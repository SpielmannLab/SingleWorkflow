#
# Read and merge several cell_data_set (Monocle3) objects following the Seurat approach
#

" Read and merge several cell_data_set (Monocle3) objects following the Seurat approach

Usage: merge_sc_seurat3_step1.R --jobname=<value> --npcs=<value> --ncores=<value> --testpcs=<logical> --infolder=<folder> --outfolder=<folder> [<file> <file>...]

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --jobname=<file>     Descriptive name for your experiment.
  --npcs=<value>       Number of principal components to use.
  --ncores=<value>     Number of threads to use.
  --testpcs=<logical>  Wheather PCs contribution to the global gene expression variance should be tested.
  --infolder=<file>    Path to the single_cell_data .rds files.  
  --outfolder=<file>   Path to results folder.

Author: 
  Cesar Prada
e-mail: 
  prada@molgen.mpg.de

"-> doc

library(docopt)
arguments <- docopt(doc, quoted_args=TRUE)
print(arguments)

# --- Dependencies

pkgs <- c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
	  'limma', 'S4Vectors', 'SingleCellExperiment',
	  'SummarizedExperiment', 'batchelor', 'devtools', 'ggplot2',
          'cowplot', "dplyr", "Seurat", "future", "sctransform",
	  'Rtsne', "phateR", "Rmagic")

if (!require("BiocManager", character.only = TRUE)) {
	install.packages("BiocManager")
	BiocManager::install()
} else {
	ipkgs <- sapply(pkgs, function(...) require(..., character.only = TRUE))
	if (any(!ipkgs)) {
		BiocManager::install(pkgs[!ipkkgs])
	} else {
		message("\n\nCool! your machine has everything is needed.\n\n")
	}
}

# Install last version of both (Beta version keep changing)

#devtools::install_github('cole-trapnell-lab/leidenbase')
#devtools::install_github('cole-trapnell-lab/monocle3', force=TRUE)


# -------------------------------
suppressMessages(library(monocle3))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(Seurat))
#suppressMessages(library(batchelor))
suppressMessages(library(sctransform))
suppressMessages(library(future))
#suppressMessages(library(dplyr))
#suppressMessages(library(ggwordcloud))
# ------------------------------

# ------- Functions
cds2mca <- function(cds) {
	mca <- CreateSeuratObject(counts = exprs(cds),
			  meta.data = data.frame(pData(cds)),
			  project = jobname)
	mca@assays$RNA@meta.features <- data.frame(fData(cds))
	mca
}

# -----------------


# --------------
# Main function
# --------------

cdsfiles <- arguments$file
inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
jobname <- arguments$jobname
ncores <- as.integer(arguments$ncores)
testpcs <- as.logical(arguments$testpcs)


# ---------------------------------------------------
# Setting multicore env for some Seurat applications
# ---------------------------------------------------
plan("multiprocess", workers = ncores)
plan()

options(future.globals.maxSize=41943040000000)

mca2tsne <- function(mca, npcs, colgroup, ncores) {
	pca <- Embeddings(object = mca@reductions$pca)
	tsne_res <- Rtsne::Rtsne(pca[, seq(npcs)], 
			 dims = 2, pca = FALSE,
			 check_duplicates = FALSE, 
			 num_threads = ncores)
	tsne_data <- data.frame(tsne_res$Y)
	row.names(tsne_data) <- rownames(pca)
	colnames(tsne_data) <- c("tSNE1", "tSNE2")
	tsne_data$id <- rownames(tsne_data)
	pdata <- mca@meta.data
	pdata$id <- rownames(pdata)
	tsne_data <- merge(tsne_data, pdata, by = "id")
	g_tsne <- ggplot(tsne_data, aes(x = tSNE1, 
					y = tSNE2, 
					color = !!rlang::sym(colgroup))) +
		    geom_point(size = 0.01) +
		    theme_classic() +
		    theme(legend.position="bottom")
	return(g_tsne)
}


prepro <- function(mca, jobname, ngenes, npcs, testpca, integration, corrected, group, ncores) {
	
	if(!integration) {
		mca <- SCTransform(mca, 
			   vars.to.regress = NULL, 
			   variable.features.n = ngenes,
			   verbose = TRUE,
			   return.only.var.genes = TRUE,
			   seed.use = 1448145
			   )
		FM <- mca@assays$SCT@scale.data # ngenes x cells
	}

	FM <- mca@assays$SCT@scale.data # ngenes x cells


	mca <- RunPCA(mca, assay = NULL, npcs = npcs,
		      rev.pca = FALSE, 
		      weight.by.var = TRUE, 
		      verbose = TRUE, ndims.print = 1:5, 
		      nfeatures.print = 30, 
		      reduction.key = "PC_",
		      seed.use = 42, 
		      approx = TRUE)

	
	if(testpca) {
		g_elbow <- ElbowPlot(mca, ndims = npcs)

		mca <- JackStraw(mca, 
				 reduction = "pca", 
				 assay = NULL, 
				 dims = npcs,
				 num.replicate = 100, 
				 prop.freq = 0.01, 
				 verbose = TRUE,
				 maxit = 1000)

		mca <- ScoreJackStraw(mca, 
				      dims = seq(npcs), 
				      score.thresh = 1e-05)
		
		g_js <- JackStrawPlot(mca, dims = seq(npcs))
	}

	mca <- RunUMAP(mca, 
		       dims = seq(npcs), 
		       min.dist = 0.75)

	mca <- FindNeighbors(mca, 
			     dims = seq(npcs), 
			     reduction = "pca",
			     verbose = TRUE)

	mca <- FindClusters(mca,
			    verbose = TRUE)

	print("hey just cluster!")

	# Color by ....
	if(group == "sample") {
		colgroup <- "sample"
	} else {
		colgroup <- "seurat_clusters"
	}

	# PHATE space & visualization >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	phate_obj <- phate(t(FM),
			  ndim = 2, 
			  knn = 5, 
			  decay = 40, 
			  n.landmark = 2000,
			  gamma = 1, 
			  t = "auto", 
			  knn.dist.method = "euclidean",
			  t.max = 100, 
			  npca = npcs,
			  n.jobs = 1 # no parallel here (...=1)
			  )

	phate_df <- data.frame(phate_obj$embedding)
	phate_df[['barcode']] <- rownames(phate_df) 

	cluster_df <- data.frame(mca@meta.data)

	phate_df <- merge(phate_df, cluster_df, by = "barcode")

	g_phate <- ggplot(data.frame(phate_df)) +
		geom_point(aes(PHATE1, PHATE2, 
			       color=!!rlang::sym(colgroup)), shape = ".") +
		theme_classic() +
		theme(legend.position="bottom")
	# <<<<<<<<<<<<<<<<<<<<<<<<<<<<< PHATE space & visualization

	g_tsne <- mca2tsne(mca, npcs, colgroup, ncores)

	g_umap <- DimPlot(mca, 
			  reduction = "umap", 
			  pt.size = 0.01, 
			  group.by=colgroup) +
		ggplot2::theme(legend.position = "bottom")

	g_pca <- DimPlot(mca, 
			 reduction = "pca", 
			 pt.size = 0.01, 
			 group.by=colgroup) +
		ggplot2::theme(legend.position = "bottom")

	g_pcafeat <- DimHeatmap(mca, 
				dims = ceiling(seq(npcs)/4), 
				cells = 500, 
				balanced = TRUE,
				fast = FALSE)

	if(corrected) {
		jobname <- paste0(jobname, '_anchoring_corrected')
	}

	title <- ggdraw() + 
		draw_label(
			jobname,
			fontface = 'bold',
			hjust = 0.5
		) +
	theme(
		plot.margin = margin(0, 0, 0, 0)
  	)

	# --- Write cell_data_set object with merged samples
	mca_file <- paste0(outputfolder, '/', jobname, '_seurat3.rds')
	saveRDS(mca, file = mca_file)
	message(paste("The cell_data_set object (Seurat3)", 
		      mca_file, "has been created."))
	# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	if(testpca) {
		return(list('mca'=mca, 
			    'umap'=g_umap, 
			    'tsne'=g_tsne, 
			    'pca'=g_pca,
			    'phate'=g_phate,
			    'elbow'=g_elbow,
			    'js'=g_js,
			    'pcafeat'=g_pcafeat,
			    'title'=title))
	} else {
		return(list('mca'=mca, 
			    'umap'=g_umap, 
			    'tsne'=g_tsne, 
			    'pca'=g_pca,
			    'phate'=g_phate,
			    'title'=title))
	}
}

# ------------
# Reading cds 
# ------------

cdss <- lapply(cdsfiles, function(...) readRDS(paste0(inputfolder, '/', ...)))

if(any(!sapply(cdss, class)=="cell_data_set")) {
	stop("Please make sure that all .rds files are cell_data_set objects.")
} 

cdsnames <- gsub("\\.rds", "", cdsfiles)
names(cdss) <- cdsnames

# 
print("This files will be SCT normalized and log10 transformed:")
print(names(cdss))
#
print("###\n###\nSeurat3 step1\n###\n###")
# From monocle3 to Seurat3 object
mcas <- lapply(cdss, cds2mca)
  
pp_mcas <- parallel::mclapply(names(mcas),
			     function(nmca) {
				     prepro(mcas[[nmca]], 
					    nmca,
					    ngenes=3000, 
					    npcs=npcs, 
					    testpca=testpcs, 
					    integration=FALSE, 
					    corrected=FALSE, 
					    group="cluster",
					    ncores=ncores #cores for tSNE
					    )
			     }, mc.cores = 2)

pdf(paste0(outputfolder, "/", jobname,
	   "_Seurat3_SCT_preprocessing_independent.pdf"), 
    width = 14, height = 4)
for(i in seq(pp_mcas)) {
	plot(
	    plot_grid(pp_mcas[[i]][['title']],
		      plot_grid(pp_mcas[[i]][['pca']], pp_mcas[[i]][['tsne']],
				pp_mcas[[i]][['umap']], pp_mcas[[i]][['phate']]+theme(legend.position="bottom"), ncol = 4),
		      ncol = 1, rel_heights = c(0.1, 1))
	)
}
dev.off()

pdf(paste0(outputfolder, "/", jobname,
	   "_Seurat3_SCT_PCmetrics_preprocessing_independent.pdf"), 
    width = 15, height = 8)
for(i in seq(pp_mcas)) {
	plot(
	    plot_grid(pp_mcas[[i]][['title']],
		      plot_grid(pp_mcas[[i]][['elbow']], pp_mcas[[i]][['js']],
				ncol = 2),
		      ncol = 1, rel_heights = c(0.1, 1))
	)
}
dev.off()

pdf(paste0(outputfolder, "/", jobname,
	   "_Seurat3_SCT_PCgenes_preprocessing_independent.pdf"), 
    width = 11, height = 52)
for(i in seq(pp_mcas)) {
	plot(
	    plot_grid(pp_mcas[[i]][['title']],
		      plot_grid(pp_mcas[[i]][['pcafeat']], ncol = 1),
		      ncol = 1, rel_heights = c(0.01, 1))
	)
}
dev.off()
