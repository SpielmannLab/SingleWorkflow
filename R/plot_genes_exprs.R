#
# Plot a given list of genes on given sc R objects 
#

" Plot a given list of genes on given sc R objects 

Usage: plot_genes_exprs.R --jobname=<value> --infolder=<folder> --outfolder=<folder> --genes=<file> --draw=<value> [<file> <file>...]

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --jobname=<file>     Descriptive name for your experiment.
  --infolder=<file>    Path to the single_cell_data .rds files.  
  --outfolder=<file>   Path to results folder.
  --genes=<file>       One column file with interesting genes. 
  --draw=<value>       Either PDF or JPEG for plot format.

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

inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
jobname <- arguments$jobname
genefile <- arguments$genes
draw <- arguments$draw
scfiles <- arguments$file

# Reading inputs
genes <- read.delim(genefile)
genes <- unique(as.character(genes[[1]]))

print("Marker genes to be plot:")
print(genes)

scfiles <- scfiles[5]

sc_objs <- mclapply(scfiles, function(...) readRDS(paste0(inputfolder, '/', ...)), mc.cores = 4)
names(sc_objs) <- gsub("\\.rds", "", scfiles)
sc_cls <- sapply(sc_objs, class)

# Monocle object visualization

cds <- sc_objs[which(sc_cls == "cell_data_set")]

if(length(genes) > 25) {
	genes <- sample(genes, 25)
}

if(draw == "PDF") {
	pdf(paste0(outputfolder, '/', jobname, ".pdf"), width = 8.27, height = 11.69)
	for(cds_name in names(cds)) {
		plot(
			plot_cells(cds[[cds_name]], genes=genes) +
				ggtitle(cds_name)
		)
	}
	dev.off()
} else if(draw == "JPEG") {
	
	for(cds_name in names(cds)) {
		jpeg(paste0(outputfolder, '/', jobname, "_", cds_name, ".jpeg"), 
		    width = 1980, height = 1980, pointsize = 74, quality = 100)
			plot(
				plot_cells(cds[[cds_name]], genes=genes) +
					ggtitle(cds_name)
			)
		dev.off()
		jpeg(paste0(outputfolder, '/', jobname, "_", cds_name, "_sample.jpeg"), 
			     width = 980, height = 980, pointsize = 74, quality = 100)
			plot(
			     plot_cells(cds[[cds_name]], color_cells_by="sample") +
				ggtitle(cds_name)
			)
		dev.off()
		pdf(paste0(outputfolder, '/', jobname, "_", cds_name, "_sample.pdf"), 
			     width = 9, height = 9)
			plot(
			     plot_cells(cds[[cds_name]], color_cells_by="sample") +
				ggtitle(cds_name)
			)
		dev.off()

	}
}

cds_name <- names(cds)
jpeg(paste0(outputfolder, '/', jobname, "_", cds_name, "_limb.jpeg"), 
		     width = 980, height = 980, pointsize = 74, quality = 100)
		plot(
		     plot_cells(cds[[cds_name]], color_cells_by="limb") +
			ggtitle(cds_name)
		)
dev.off()

# Seurat marker visualizaion
mca <- sc_objs[which(sc_cls == "Seurat")]

if(length(genes) > 25) {
	genes <- sample(genes, 25)
}
genes

ref <- read.delim("~/limb_project/data/biomart_mm_20191206.tsv")
head(ref)

genes <- as.character(filter(ref, Gene.name%in%genes)$Gene.stable.ID)

jpeg(paste0(outputfolder, '/', jobname, "_",  "seurat3.jpeg"), 
	    width = 1980, height = 1980, pointsize = 74, quality = 100)
		plot(
		     FeaturePlot(mca[[1]], 
			  reduction = "umap",
			  pt.size = 0.01, 
			  features=genes) +
				ggplot2::theme(legend.position = "none") + 
				ggtitle(jobname) +
				scale_colour_gradient(low = "grey", high = "blue")
		)
dev.off()
