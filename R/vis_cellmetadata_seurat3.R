# Visualize cell metadata in UMAP space (Seurat3)

" Visualize cell metadata in UMAP space (Seurat3)

Usage: vis_cellmetadata_seurat3.R --jobname=<value> [--mcafile=<file>] --infolder=<folder> --outfolder=<folder> --ncores=<value> 

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --jobname=<file>     Descriptive name for your experiment.
  --mcafile=<file>     Optional. If !is.null, do not use jobname to guess mcafilename.
  --infolder=<file>    Path to mca .rds object file 
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
	  'Rtsne')

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


inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
jobname <- arguments$jobname

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
	mca_file <- paste0(inputfolder, '/', jobname, '_seurat3_merged_clustered.rds')
}

mca <- readRDS(mca_file)

message(paste("The cell_data_set object (Seurat3)", 
	      mca_file, "has been read."))

metavar <- colnames(mca@meta.data) 

mt <- lapply(metavar, function(colgroup) {
		if(class(mca@meta.data[[colgroup]]) %in% c("integer", "numeric")) {

			mca@meta.data[[colgroup]] <- as.numeric(mca@meta.data[[colgroup]])
		
			print("Num:")
			print(colgroup)
			print(class(mca@meta.data[[colgroup]]))

			FeaturePlot(mca, 
			  reduction = "umap", 
			  pt.size = 0.01, 
			  features=colgroup) +
		 		#ggplot2::theme(legend.position = "none") +
				ggplot2::theme(legend.position = "none") + 
				ggtitle(colgroup) +
				scale_colour_gradient(low = "grey", high = "blue")
	
		} else if(class(mca@meta.data[[colgroup]]) %in% c("character", "factor")) {
			
			mca@meta.data[[colgroup]] <- as.factor(mca@meta.data[[colgroup]])
			
			print("Fac:")
			print(colgroup)
			print(class(mca@meta.data[[colgroup]]))
	
			DimPlot(mca, 
			  reduction = "umap", 
			  pt.size = 0.01, 
			  group.by=colgroup) +
		 		#ggplot2::theme(legend.position = "none") +
				ggplot2::theme(legend.position = "none") + 
				ggtitle(colgroup)
		}
	})

rfile <- paste0(outputfolder, '/', jobname, '_', 'metadata_umap.pdf')

# Plot metadata in UMAP space << PDF >>
title <- ggdraw() + 
	  draw_label(
	    jobname,
	    fontface = 'bold',
	    hjust = 0.5
	  ) +
	  theme(
	    plot.margin = margin(0, 0, 0, 0)
	  )
pdf(rfile, height = 22.69, width = 8.27) 
    plot_grid(title, plot_grid(plotlist = mt, ncol = 2), 
	      ncol = 1, rel_heights = c(0.02, 1))
dev.off()

# Plot metadata in UMAP space << PDF >>

rfile <- paste0(outputfolder, '/', jobname, '_', 'metadata_umap.jpeg')

jpeg(rfile, width = 1980, height = 4980, pointsize = 74, quality = 100)
    plot_grid(title, plot_grid(plotlist = mt, ncol = 3), 
	      ncol = 1, rel_heights = c(0.02, 1))
dev.off()

