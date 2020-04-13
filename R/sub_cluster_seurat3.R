#
# Sub-cluster a merged seurat3 object based on a given colData variable
#

"  Sub-cluster a merged seurat3 object based on a given colData variable

Usage: sub_cluster_seurat3.R --jobname=<value> [--mcafile=<file>] --npcs=<value> --testpcs=<logical> --cellvar=<value> [--samplevar=<value>] --nvgenes=<value> --mincells=<value> --infolder=<folder> --outfolder=<folder> --ncores=<value> 

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --jobname=<file>     Descriptive name for your experiment.
  --mcafile=<file>     Optional. If !is.null, do not use jobname to guess mcafilename.
  --npcs=<value>       Number of PCs to consider for clustering.
  --testpcs=<logical>  Wheather the mca input object contains JackStraw-test results.
  --cellvar=<value>    colData variable to use as reference for sub-clustering.
  --samplevar=<value>  Optional. If more than one sample that need to be integrated. 
  --nvgenes=<value>    Number of variant genes to consider per sub_cluster pre-processing.
  --mincells=<value>   Min. number of cells per sample to consider in pre-processing.
  --infolder=<file>    Path to the single_cell_data .rds files.  
  --outfolder=<file>   Path to results folder.
  --ncores=<value>     Number of processors to use

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
#suppressMessages(library(sctransform))
#suppressMessages(library(future))
#suppressMessages(library(dplyr))
#suppressMessages(library(ggwordcloud))
# ------------------------------

# ---------- Functions ---------

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

# Load parameters
inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
jobname <- arguments$jobname
cellvar <- arguments$cellvar
npcs <- as.numeric(arguments$npcs)
nvgenes <- as.numeric(arguments$nvgenes)
mincells <- as.numeric(arguments$mincells)
testpcs <- as.logical(arguments$testpcs)

if(!is.null(arguments$samplevar)) {
	samplevar <- arguments$samplevar
} 

# ---------------------------------------------------
# Setting multicore env for some Seurat applications
# ---------------------------------------------------
ncores <- as.numeric(arguments$ncores)
plan("multiprocess", workers = ncores)
plan()

# --- Read cell_data_set object with merged samples
if(!is.null(arguments$mcafile)) {
	mca_file <- paste0(inputfolder, '/', arguments$mcafile)
	jobname <- gsub("\\.rds", "", arguments$mcafile)
} else {
	mca_file <- paste0(inputfolder, '/', jobname, '_seurat3_merged_clustered_resrange.rds')
}

mca <- readRDS(mca_file)

message(paste("The cell_data_set object (Seurat3)", 
	      mca_file, "has been read."))

# --- Check colData variable existance

mca@meta.data$barcode <- rownames(mca@meta.data)
cell_meta <- mca@meta.data
cellvar_levels <- unique(cell_meta[[cellvar]])

if(length(cellvar_levels) > 1) {

	print("Subcluster will be based on these cellvar levels:")	
	print(cellvar_levels)

} else {

	stop("Please provide a valid colvar containing more than one level.")

}

# --- Split mca
DefaultAssay(mca) <- "RNA" # <<<<<<<<<<<<<<<<<< critical!

mcas <- lapply(cellvar_levels, function(l) {
	l_cells <- dplyr::filter(cell_meta, !!rlang::sym(cellvar) == l)[['barcode']]
	mca[, l_cells]
})

mcas <- setNames(mcas, cellvar_levels)

options(future.globals.maxSize=41943040000000)

# Setting sub-folder to write sub mca objects

fls <- list.dirs(path = outputfolder)

if(any(grepl(jobname, fls))) {

	print(paste("The folder", jobname, "in the indicated outputfolder will be overwrite."))

} else {

	system(paste0("mkdir ", outputfolder, "/", jobname))
	print(paste("The folder", paste0(outputfolder, "/", jobname),
		    "was created."))

}

rds_outputfolder <- paste0(outputfolder, "/", jobname)

lapply(seq(mcas), function(i) {
	       #
	       print(names(mcas)[i])
	       #
	       cl_mca <- mcas[[i]]
	       #
	       cl_mca@meta.data$barcode <- rownames(cl_mca@meta.data)
	       cell_meta <- cl_mca@meta.data
	       cellvar_levels <- unique(cell_meta[[samplevar]])

	       if(length(cellvar_levels) > 1) {
		       
		       print("Subcluster will be based on these cellvar levels:")	
		       print(cellvar_levels)

		} else {

			stop("Please provide a valid samplevar containing more than one level.")
		
	       }
	       #
	       # --- Split mca in the compprising samples
	       #
	       s_cl_mcas <- lapply(cellvar_levels, function(l) {
				l_cells <- dplyr::filter(cell_meta, !!rlang::sym(samplevar) == l)[['barcode']]
				cl_mca[, l_cells]
			})
	       #
	       print(sapply(s_cl_mcas, dim))
	       #
	       # filter samples contributing with less than mincell cells 
	       s_cl_mcas <- s_cl_mcas[which(sapply(s_cl_mcas, dim)[2, ] > mincells)]
	       print(sapply(s_cl_mcas, dim))
	       # Normalization of the sample-wise and cluster specific datasets
	       s_cl_mcas <- lapply(s_cl_mcas, function(mcai) {

				 SCTransform(object = mcai,
				     verbose = TRUE,
				     return.only.var.genes = FALSE)
		})
	       #
	       print(sapply(s_cl_mcas, print))
	       #
	       # Integration removing samples effect only
	       #
	       mca.features <- SelectIntegrationFeatures(object.list = s_cl_mcas)#, 
							 #nfeatures = nvgenes)
	       #
	       print(head(mca.features))
	       #
	       s_cl_mcas <- PrepSCTIntegration(object.list = s_cl_mcas, 
					       anchor.features = mca.features,
					       verbose = TRUE)
	       #
	       # min cells - 1 might be a maximum nearest neighbors 
	       #
	       minK <- min(sapply(s_cl_mcas, dim)[2, ])-1
	       print(minK)
	       if(minK > 200) {
		       minK <- 200
	       } else if (minK > 100) {
		       minK <- 100
	       }
	       print(minK)
	       #
	       mca.anchors <- FindIntegrationAnchors(object.list = s_cl_mcas,
						     normalization.method = "SCT", 
						     anchor.features = mca.features, 
						     #k.anchor = minK,
						     k.filter = minK,
						     #dims = seq(maxDim),
						     #k.score = minK,
						     verbose = TRUE)
	       #
	       # integrate data and keep full geneset
	       #
	       i_s_cl_mca <- IntegrateData(anchorset = mca.anchors,
					   normalization.method = "SCT",
					   features.to.integrate = mca.features, 
					   verbose = TRUE)
	       #
	       # Dimensionality reduction
	       #
	       FM <- i_s_cl_mca@assays$SCT@scale.data # ngenes x cells
	       #
	       i_s_cl_mca <- RunPCA(i_s_cl_mca, 
			      assay = NULL, 
			      npcs = npcs,
			      rev.pca = FALSE, 
			      weight.by.var = TRUE, 
			      verbose = TRUE, 
			      ndims.print = 1:5, 
			      nfeatures.print = 30, 
			      reduction.key = "PC_",
			      seed.use = 42, 
			      approx = TRUE)
	       #
	       if(testpcs) {
		       #
		       g_elbow <- ElbowPlot(i_s_cl_mca,
					    ndims = npcs)
		       
		       i_s_cl_mca <- JackStraw(i_s_cl_mca, 
					       reduction = "pca", 
					       assay = NULL, 
					       dims = npcs,
					       num.replicate = 100, 
					       prop.freq = 0.01, 
					       verbose = TRUE,
					       maxit = 1000)

		       i_s_cl_mca <- ScoreJackStraw(i_s_cl_mca, 
						    dims = seq(npcs), 
						    score.thresh = 1e-05)

		       g_js <- JackStrawPlot(i_s_cl_mca, dims = seq(npcs))
		}
	       #
	       i_s_cl_mca <- RunUMAP(i_s_cl_mca, 
				     dims = seq(npcs), 
				     min.dist = 0.2)
	       # clustering
	       i_s_cl_mca <- FindNeighbors(i_s_cl_mca, 
					   dims = seq(npcs), 
					   verbose = TRUE)
	       #
	       i_s_cl_mca <- FindClusters(i_s_cl_mca,
					  #resolution = 0.07,
					  verbose = TRUE)
	       #
	       print("hey just cluster!")
	       #
	       # writing down the mca object
	       mca_file <- paste0(rds_outputfolder, '/', jobname, 
				  gsub("\\.", "_", cellvar), 
				  "_cluster", "_", names(mcas)[i], 
				  "_Seurat3_preprocessed_sample_integrated_subcluster.rds")
	       print(i_s_cl_mca)
	       saveRDS(i_s_cl_mca, file = mca_file)
	       message(paste("The cell_data_set object (Seurat3)", 
			 mca_file, "has been created."))
	       #
	       # Color by ....
	       colgroup <- samplevar
	       #

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
	       cluster_df[['barcode']] <- rownames(cluster_df)
	       
	       phate_df <- merge(phate_df, cluster_df, by = "barcode")


		# <<<<<<<<<<<<<<<<<<<<<<<<<<<<< PHATE space & visualization
		#
		# tSNE
		g_tsne <- mca2tsne(i_s_cl_mca, 
				   npcs, 
				   colgroup, 
				   ncores)
		# umap
		g_umap <- DimPlot(i_s_cl_mca, 
			  reduction = "umap", 
			  pt.size = 0.01, 
			  group.by=colgroup) +
       			ggplot2::theme(legend.position = "bottom")
		# pca
		g_pca <- DimPlot(i_s_cl_mca, 
				 reduction = "pca", 
				 pt.size = 0.01, 
				 group.by=colgroup) +
       			ggplot2::theme(legend.position = "bottom")
		# phate
		g_phate <- ggplot(data.frame(phate_df)) +
				geom_point(aes(PHATE1, PHATE2, 
					       color=!!rlang::sym(colgroup)), shape = ".") +
                    		theme_classic() +
				theme(legend.position="bottom")

		if(testpcs) {
			
			l <- list('mca'=i_s_cl_mca, 
				  'umap'=g_umap, 
				  'tsne'=g_tsne,
				  'pca'=g_pca,
				  'phate'=g_phate,
				  'elbow'=g_elbow,
				  'js'=g_js,
				  'title'=title)
		} else {

			l <- list('mca'=i_s_cl_mca, 
				  'umap'=g_umap, 
				  'tsne'=g_tsne, 
				  'pca'=g_pca,
				  'phate'=g_phate,
				  'title'=title)
		}
		#
		p_sub_mca <- l
		#
		pdf(paste0(rds_outputfolder, "/", jobname, 
			   '_',
			   gsub("\\.", "_", cellvar), 
			   '_', 
			   "cluster", "_", names(mcas)[i], 
			   "_",
			   "Seurat3_subclusters.pdf"), 
		    width = 13, height = 5)

		plot(
		     plot_grid(p_sub_mca[['pca']], p_sub_mca[['tsne']],
			       p_sub_mca[['umap']], p_sub_mca[['phate']], 
				ncol = 4)
		     )
		#
		dev.off()
		jpeg(paste0(rds_outputfolder, "/", jobname, 
			   '_',
			   gsub("\\.", "_", cellvar), 
			   '_', 
			   "cluster", "_", names(mcas)[i], 
			   "_",
			   "Seurat3_subclusters.jpeg"), 
		    height = 450, width = 1600, pointsize = 74, quality = 100)

		plot(
		     plot_grid(p_sub_mca[['pca']], p_sub_mca[['tsne']],
			       p_sub_mca[['umap']], p_sub_mca[['phate']], 
				ncol = 4)
		     )
		#
		dev.off()
})
