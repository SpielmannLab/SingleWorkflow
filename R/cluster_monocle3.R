#
# Cluster a merged cell_data_set (Monocle3) object using several resolutions
#

" Cluster a merged cell_data_set (Monocle3) object using several resolutions

Usage: cluster_monocle3.R --jobname=<value> [--cdsfile=<file>] --npcs=<value> --infolder=<folder> --outfolder=<folder> --ncores=<value> 

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --jobname=<file>     Descriptive name for your experiment.
  --cdsfile=<file>     Optional. If !is.null, do not use jobname to guess mcafilename.
  --npcs=<value>       Number of principal components to use.
  --infolder=<file>    Path to the single_cell_data .rds files.  
  --outfolder=<file>   Path to results folder.
  --ncores=<value>     Number of threads to use.

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
          'cowplot', "dplyr", "phateR", "Rmagic", "monocle3")

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
suppressMessages(library(monocle3))
pymagic <- reticulate::import("magic")
#suppressMessages(library(batchelor))
#suppressMessages(library(dplyr))
#suppressMessages(library(ggwordcloud))
# ------------------------------

# --------------
# Main function
# --------------

# ---------- Functions ---------

getcluster <- function(cds, res) {
	cds <- cluster_cells(cds, 
			     reduction_method = "Aligned",#c("UMAP", "tSNE", "PCA", "LSI", "Aligned"), 
			     k = 20, 
			     cluster_method = "louvain",
			     num_iter = 2, 
			     #partition_qval = 0.05, 
			     weight = FALSE,
			     resolution = res, 
			     random_seed = NULL,
			     verbose = TRUE)

	cds@colData$tmpcluster <- cds@clusters$Aligned$clusters
	
	g_cluster <- plot_cells(cds,
				color_cells_by="tmpcluster")
	
	g_partition <- plot_cells(cds, 
				  color_cells_by="partition", 
				  group_cells_by="partition")

	clstr <- data.frame('barcode' = names(cds@clusters$Aligned$clusters),
			    'cluster' = cds@clusters$Aligned$clusters)

	colnames(clstr)[2] <- paste0(colnames(clstr)[2], '_', 
				     gsub("\\.", "_", res))

	return(list('g_cl' = g_cluster,
		    'g_pr' = g_partition,
		    'cluster' = clstr))
}

# --- Read cell_data_set object with merged samples

if(!is.null(arguments$mcafile)) {
	cds_file <- paste0(inputfolder, '/', arguments$cdsfile)
	jobname <- gsub("\\.rds", "", arguments$cdsfile)
} else {
	cds_file <- paste0(inputfolder, '/', jobname, '_moncocle3_merged.rds')
}

cds <- readRDS(cds_file)

message(paste("The cell_data_set object (Monocle3)", 
		cds_file, "has been read."))

# clustering resolutions
res <- c(0.000001, 0.00001, 0.0001, 0.001, 0.05) 
	
cls <- mclapply(res, function(r) getcluster(cds, r), mc.cores = ncores)
 
cluster_jobname <- paste0(jobname, '_', ncol(cds), '_cells_and_', npcs, 'PCs_from_', nrow(cds), "_genes")  

# --- Write QC plots
title <- ggdraw() + 
  draw_label(
    cluster_jobname,
    fontface = 'bold',
    hjust = 0.5
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

rfile <- paste0(outputfolder, '/', cluster_jobname, '_', 'monocle3_clustering.pdf')

print(rfile)

cl_gplots <- lapply(cls, function(cl) cl[c('g_cl', 'g_pr')])
cl_gplots <- unlist(cl_gplots, recursive = FALSE)

pdf(rfile, height = 11.69, width = 8.27) 

    plot(
	 plot_grid(title, plot_grid(plotlist = cl_gplots, ncol = 2), 
	      ncol = 1, rel_heights = c(0.02, 1))
	 )

dev.off()

jpeg(paste0(outputfolder, '/', cluster_jobname, "_", "monocle3_clustering.jpeg"), 
		   height = 1980, width = 2980, pointsize = 74, quality = 100)

    plot(
    	plot_grid(title, plot_grid(plotlist = cl_gplots, ncol = 4), 
	      ncol = 1, rel_heights = c(0.02, 1))
    )

dev.off()

meta_cls <- Reduce(cbind, lapply(cls, function(cl) select(cl[['cluster']], matches("cluster"))))
meta_cls[['barcode']] <- rownames(meta_cls)
meta <- colData(cds)
meta[['barcode']] <- rownames(meta)
meta <- merge(meta, meta_cls, by = "barcode")
rownames(meta) <- meta[['barcode']]
pData(cds) <- meta

cds_file <- paste0(outputfolder, '/', jobname, '_monocle3_merged_clustered_resrange.rds')

print(cds)

saveRDS(cds, file = cds_file)

message(paste("The cell_data_set object (Monocle3)", 
	      cds_file, "has been created."))
