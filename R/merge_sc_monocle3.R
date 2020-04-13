#
# Read and merge several cell_data_set (Monocle3) objects
#

" Read and merge several cell_data_set (Monocle3) objects

Usage: merge_sc_monocle3.R --jobname=<value> --npcs=<value> --ncores=<value> --infolder=<folder> --outfolder=<folder> [<file> <file>...]

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --jobname=<file>     Descriptive name for your experiment.
  --npcs=<value>       Number of principal components to use.
  --ncores=<value>     Number of threads to use.
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
          'cowplot', "dplyr", "phateR", "Rmagic")

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
# Phate installation
#reticulate::py_install("phate", pip=TRUE)
#reticulate::py_install("magic-impute", pip=TRUE)
#devtools::install_github("KrishnaswamyLab/phateR")
#devtools::install_github("KrishnaswamyLab/MAGIC/Rmagic")


# -------------------------------
suppressMessages(library(monocle3))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(phateR))
suppressMessages(library(Rmagic))
pymagic <- reticulate::import("magic")
#suppressMessages(library(batchelor))
#suppressMessages(library(dplyr))
#suppressMessages(library(ggwordcloud))
# ------------------------------

# --------------
# Main function
# --------------

cdsfiles <- arguments$file
inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
jobname <- arguments$jobname
npcs <- as.numeric(arguments$npcs)
ncores <- as.integer(arguments$ncores)
# ------------
# Reading cds 
# ------------

# --------------------------
# Monocle3 pre-processing
# --------------------------

cdss <- lapply(cdsfiles, function(...) readRDS(paste0(inputfolder, '/', ...)))

if(any(!sapply(cdss, class)=="cell_data_set")) {
	stop("Please make sure that all .rds files are cell_data_set objects.")
} 

cdsnames <- gsub("\\.rds", "", cdsfiles)
names(cdss) <- cdsnames

# 
print("This files will be size.factor normalized and log10 transformed:")
print(names(cdss))
#

# Factor.size + Log2 normalization + dimensionality reduction (monocle apprach)

sflog <- function(cds, jobname, integration, corrected, npcs) {
	if(!corrected) {
		cds <- preprocess_cds(cds, method = "PCA", num_dim = npcs,
		      norm_method = "log", use_genes = NULL,
		      scaling = TRUE,  verbose = TRUE, 
		      cores = ncores)
		# Monocle3 normalizarion for norm_method = "log"
		FM <- SingleCellExperiment::counts(cds)
		pseudo_count <- 1 # because norm_method = "log"
		FM <- Matrix::t(Matrix::t(FM)/size_factors(cds))
		FM@x = log2(FM@x + 1)
	}

	if(!corrected) {
		premethod <- "PCA"
	} else {
		premethod <- "Aligned"
		# Monocle3 normalizarion for norm_method = "log"
		FM <- SingleCellExperiment::counts(cds)
		pseudo_count <- 1 # because norm_method = "log"
		FM <- Matrix::t(Matrix::t(FM)/size_factors(cds))
		FM@x = log2(FM@x + 1)
		cores <- MulticoreParam(ncores)
		MNNres <- fastMNN(FM, batch=colData(cds)$sample,
				  d = npcs, correct.all=TRUE, BPPARAM=cores)
		FM <- MNNres@assays@data$reconstructed
		rownames(FM) <- MNNres@rowRanges@partitioning@NAMES 
	}

	cds <- reduce_dimension(cds, max_components = 2, reduction_method = "UMAP",
				preprocess_method = premethod, umap.metric = "cosine", 
				umap.min_dist = 0.1, umap.n_neighbors = 15L, 
				umap.fast_sgd = TRUE, umap.nn_method = "annoy", 
				cores = ncores, verbose = TRUE)
	
	cds <- reduce_dimension(cds, max_components = 2, reduction_method = "tSNE",
				preprocess_method = premethod, cores = ncores, 
				verbose = TRUE)
	
	cds <- cluster_cells(cds, reduction_method="UMAP")

	cds@clusters$tSNE <- cds@clusters$UMAP
	cds@clusters$PCA <- cds@clusters$UMAP

	if(integration) {
		cl <- "sample"
	} else {
		cl <- "cluster"
	}

	g2 <- plot_cells(cds, color_cells_by=cl, 
			 reduction_method = "UMAP",
			 label_cell_groups=FALSE)


	g3 <- plot_cells(cds, color_cells_by=cl, 
			 reduction_method = "tSNE",
			 label_cell_groups=FALSE)

	g4 <- plot_cells(cds, color_cells_by=cl, 
			 reduction_method = "PCA",
			 label_cell_groups=FALSE)

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

	cluster_df <- data.frame(cds@clusters$UMAP$clusters)
	cluster_df[['barcode']] <- rownames(cluster_df)

	if(integration) {
		meta <- colData(cds)
		meta[['barcode']] <- rownames(meta)
		meta <- merge(meta, cluster_df, by = "barcode")
		phate_df <- merge(phate_df, meta, by = "barcode")
	} else {
		meta <- colData(cds)
		meta[['barcode']] <- rownames(meta)
		meta <- merge(meta, cluster_df, by = "barcode")
		phate_df <- merge(phate_df, meta, by = "barcode")
	}

	if(cl == "cluster") {
		colnames(phate_df)[which(colnames(phate_df) == "cds.clusters.UMAP.clusters")] <- "cluster"
	}

	g_phate <- ggplot(data.frame(phate_df)) +
		geom_point(aes(PHATE1, PHATE2, 
			       color=!!rlang::sym(cl)), shape = ".") +
		theme_classic()
	# <<<<<<<<<<<<<<<<<<<<<<<<<<<<< PHATE space & visualization

	if(corrected) {
		jobname <- paste0(jobname, '_MNN_corrected')
	}

	title <- ggdraw() + 
		draw_label(
			jobname,
			fontface = 'bold',
			#x = 0,
			hjust = 0.5
		) +
	theme(
		# add margin on the left of the drawing canvas,
		# so title is aligned with left edge of first plot
		plot.margin = margin(0, 0, 0, 0)
  	)

	return(list('cds'=cds, 'umap'=g2, 'tsne'=g3, 'pca'=g4, 
		    'phate'=g_phate, 'title'=title))
}

# Size.factor and log2 transformed
pp_cds <- parallel::mclapply(names(cdss), 
			     function(ncds) sflog(cdss[[ncds]], ncds, FALSE, FALSE, npcs), 
			     mc.cores = ncores)

pdf(paste0(outputfolder, "/", jobname, "_monocle3_preprocessing_independent.pdf"), 
    width = 14, height = 3)
for(i in seq(pp_cds)) {
	plot(
	    plot_grid(pp_cds[[i]][['title']],
		      plot_grid(pp_cds[[i]][['pca']], pp_cds[[i]][['tsne']],
				pp_cds[[i]][['umap']], pp_cds[[i]][['phate']], ncol = 4),
		      ncol = 1, rel_heights = c(0.1, 1))
	)
}
dev.off()

ncdss <- lapply(pp_cds, function(pcds) pcds[['cds']])
names(ncdss) <- names(cdss)

# Holding QC metrics for all samples
cds <- lapply(names(cdss), function(cdsname) {
	scds <- cdss[[cdsname]]
	oldcols <- seq(ncol(fData(scds)))[-c(1:2)]
	colnames(fData(scds))[oldcols] <- paste0(colnames(fData(scds))[oldcols], 
						'_', cdsname)
	scds
})

names(cds) <- cdsnames
cds <- combine_cds(cdss, keep_all_genes = FALSE, cell_names_unique = FALSE)
pData(cds)$sample <- as.factor(pData(cds)$sample)

# Combined samples with no MNN correction
pp_m_cds <- sflog(cds, jobname, TRUE, FALSE, npcs) # color by sample
pp_m_cds_cl <- sflog(cds, jobname, FALSE, FALSE, npcs) # color by cluster

cds <- pp_m_cds[['cds']]

# kMNN-Marioni's approach to remove the sample variance effect
cds <- align_cds(cds, preprocess_method="PCA",
		 alignment_group="sample")


# aligned / batch-corrected cds
pp_m_al_cds <- sflog(cds, jobname, TRUE, TRUE, npcs) # color by sample 
pp_m_al_cds_cl <- sflog(cds, jobname, FALSE, TRUE, npcs) # color by cluster

pdf(paste0(outputfolder, "/", jobname, "_monocle3_preprocessing_allsamples.pdf"), width = 25, height = 5)
    plot_grid(pp_m_cds[['title']],
	      plot_grid(pp_m_cds[['pca']], pp_m_cds[['tsne']],
			pp_m_cds[['umap']], pp_m_cds[['phate']], ncol = 4),
	      ncol = 1, rel_heights = c(0.1, 1)
    )
    plot_grid(pp_m_cds_cl[['title']],
	      plot_grid(pp_m_cds_cl[['pca']], pp_m_cds_cl[['tsne']],
			pp_m_cds_cl[['umap']], pp_m_cds_cl[['phate']], ncol = 4),
	      ncol = 1, rel_heights = c(0.1, 1)
    )
    plot_grid(pp_m_al_cds[['title']],
	      plot_grid(pp_m_al_cds[['pca']], pp_m_al_cds[['tsne']],
			pp_m_al_cds[['umap']], pp_m_al_cds[['phate']], ncol = 4),
	      ncol = 1, rel_heights = c(0.1, 1)
    )
    plot_grid(pp_m_al_cds_cl[['title']],
	      plot_grid(pp_m_al_cds_cl[['pca']], pp_m_al_cds_cl[['tsne']],
			pp_m_al_cds_cl[['umap']], pp_m_al_cds_cl[['phate']], ncol = 4),
	      ncol = 1, rel_heights = c(0.1, 1)
    )
dev.off()

png(paste0(outputfolder, "/", jobname, "_monocle3_preprocessing_allsamples_sample_non_aligned.png"), 
    width = 2500, height = 500)
    plot_grid(pp_m_cds[['title']],
	      plot_grid(pp_m_cds[['pca']], pp_m_cds[['tsne']],
			pp_m_cds[['umap']], pp_m_cds[['phate']], ncol = 4),
	      ncol = 1, rel_heights = c(0.1, 1)
    )
dev.off()

png(paste0(outputfolder, "/", jobname, "_monocle3_preprocessing_allsamples_cluster_non_aligned.png"), 
    width = 2500, height = 500)
   plot_grid(pp_m_cds_cl[['title']],
	      plot_grid(pp_m_cds_cl[['pca']], pp_m_cds_cl[['tsne']],
			pp_m_cds_cl[['umap']], pp_m_cds_cl[['phate']], ncol = 4),
	      ncol = 1, rel_heights = c(0.1, 1)
    )
dev.off()

png(paste0(outputfolder, "/", jobname, "_monocle3_preprocessing_allsamples_sample_aligned.png"), 
    width = 2500, height = 500)
   plot_grid(pp_m_al_cds[['title']],
	      plot_grid(pp_m_al_cds[['pca']], pp_m_al_cds[['tsne']],
			pp_m_al_cds[['umap']], pp_m_al_cds[['phate']], ncol = 4),
	      ncol = 1, rel_heights = c(0.1, 1)
    )
dev.off()

png(paste0(outputfolder, "/", jobname, "_monocle3_preprocessing_allsamples_cluster_aligned.png"), 
    width = 2500, height = 500)
    plot_grid(pp_m_al_cds_cl[['title']],
	      plot_grid(pp_m_al_cds_cl[['pca']], pp_m_al_cds_cl[['tsne']],
			pp_m_al_cds_cl[['umap']], pp_m_al_cds_cl[['phate']], ncol = 4),
	      ncol = 1, rel_heights = c(0.1, 1)
    )
dev.off()


# --- Write cell_data_set object with merged samples

cds_file <- paste0(outputfolder, '/', jobname, '_moncocle3_merged.rds')

cds <- pp_m_al_cds_cl[['cds']]
cds
saveRDS(cds, file = cds_file)

message(paste("The cell_data_set object (Monocle3)", 
	      cds_file, "has been created."))
