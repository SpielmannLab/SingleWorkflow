#
# Estimate cell type composition difference given a colData variable (Seurat3 objects only)
#

"  Estimate cell type composition difference given a colData variable (Seurat3 objects only)

Usage: diff_cellcomp.R --jobname=<value> --groupvar=<value> --resolution=<value> --infolder=<folder> --outfolder=<folder>

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --jobname=<file>     Descriptive name for your experiment.
  --groupvar=<value>   Variable to gather cells.
  --resolution=<value> Clustering resolution.
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
	  'Rtsne', 'MASS')

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
suppressMessages(library(MASS))
#suppressMessages(library(dplyr))
#suppressMessages(library(ggwordcloud))
# ------------------------------


inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
jobname <- arguments$jobname
groupvar <- arguments$groupvar
resolution <- arguments$resolution

# --- Read cell_data_set object with merged samples
mca_file <- paste0(outputfolder, '/', jobname, '_seurat3_merged_clustered_resrange.rds')

mca <- readRDS(mca_file)

message(paste("The cell_data_set object (Seurat3)", 
	      mca_file, "has been read."))

# Setting default cell identity
cluster_var <- paste0("integrated_snn_res.", resolution) 
Idents(mca) <- cluster_var

cond <- DimPlot(mca, 
	  reduction = "umap", 
	  pt.size = 0.01, 
	  split.by=groupvar) + 
		ggplot2::theme(legend.position = "bottom") + 
		ggtitle(groupvar)

cond2 <- DimPlot(mca, 
		  reduction = "umap", 
		  pt.size = 0.01, 
		  group.by=groupvar) +
			ggplot2::theme(legend.position = "bottom") + 
			ggtitle(groupvar)
# --- title
title <- ggdraw() + 
  draw_label(
    paste0(jobname, ' ', groupvar),
    fontface = 'bold',
    hjust = 0.5
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

jpeg(paste0(outputfolder, '/', jobname, "_", "pickedup_clusters_umap.jpeg"), 
		    width = 1380, height = 900, pointsize = 74, quality = 100)
    plot_grid(title, plot_grid(cond, cond2, ncol = 2), 
	      ncol = 1, rel_heights = c(0.02, 1))
dev.off()


umap <- mca@reductions$umap@cell.embeddings

mca@meta.data$barcode <- rownames(mca@meta.data)

meta_data <- mca@meta.data

groupvar_levels <- unique(meta_data[[groupvar]])

dens_plots <- lapply(groupvar_levels, function(l) {
	l_cells <- dplyr::filter(meta_data, !!rlang::sym(groupvar) == l)[['barcode']]
	x <- umap[l_cells, 'UMAP_1']
	y <- umap[l_cells, 'UMAP_2']
	kde2d(y, x, n = 1000)
})

dens_plots <- setNames(dens_plots, groupvar_levels)

# PDF 2D density plots
pdf(paste0(outputfolder, '/', groupvar, "_differential_density_umap.pdf"))
	for(l in groupvar_levels) {
		filled.contour(dens_plots[[l]], 
			       zlim = seq(0, 0.01, 0.001), 
				plot.title = title(main = l))
	}
dev.off()

# JPEG 2D density plots
for(l in groupvar_levels) {
	jpeg(paste0(outputfolder, '/', groupvar, "_", l, 
		    "_differential_density_umap.jpeg"))
		     filled.contour(dens_plots[[l]], zlim = seq(0, 0.01, 0.001), 
				plot.title = title(main = l))
	dev.off()
}

	
print(head(meta_data))

# Plots for differential cell proportion

prop_plots <- lapply(groupvar_levels, function(l) {
	ggplot(dplyr::filter(meta_data, !!rlang::sym(groupvar) == l), 
	       aes(x = !!rlang::sym(groupvar), 
		   fill = !!rlang::sym(cluster_var))) +
		     geom_bar(position = "stack", 
			      aes(y = ..count../sum(..count..))) +
		     scale_y_continuous(labels = scales::percent_format()) +
		     theme_classic()
})

jpeg(paste0(outputfolder, '/', jobname, "_", groupvar, "_cell_prop_condition.jpeg"), 
		    width = 1380, height = 900, pointsize = 74, quality = 100)
    plot_grid(title, plot_grid(plotlist = prop_plots, 
			       ncol = length(groupvar_levels)), 
	      ncol = 1, rel_heights = c(0.02, 1))
dev.off()

pdf(paste0(outputfolder, '/', jobname, "_", groupvar, "_cell_prop_condition.pdf"), 
		    width = 13, height = 9)
    plot_grid(title, plot_grid(plotlist = prop_plots, 
			       ncol = length(groupvar_levels)), 
	      ncol = 1, rel_heights = c(0.02, 1))
dev.off()


prop_plots <- lapply(groupvar_levels, function(l) {
	ggplot(dplyr::filter(meta_data, !!rlang::sym(groupvar) == l), 
	       aes(x = !!rlang::sym(cluster_var), 
		   fill = !!rlang::sym(cluster_var))) +
		     geom_bar(position = position_dodge(), 
			      aes(y = ..count../sum(..count..))) +
		     scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
		     theme_classic()
})

jpeg(paste0(outputfolder, '/', jobname, "_", groupvar, "_barplot_cell_prop_condition.jpeg"), 
		    width = 1380, height = 900, pointsize = 74, quality = 100)
    plot_grid(title, plot_grid(plotlist = prop_plots, 
			       ncol = length(groupvar_levels)), 
	      ncol = 1, rel_heights = c(0.02, 1))
dev.off()

pdf(paste0(outputfolder, '/', jobname, "_", groupvar, "_barplot_cell_prop_condition.pdf"), 
		    width = 13, height = 9)
    plot_grid(title, plot_grid(plotlist = prop_plots, 
			       ncol = length(groupvar_levels)), 
	      ncol = 1, rel_heights = c(0.02, 1))
dev.off()
