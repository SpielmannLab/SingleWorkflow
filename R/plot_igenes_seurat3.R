#
# Plot interesting genes an/or pathways for seurat3 object for given cell-phenodata variable and clustering resolution ### in dev ###
#

" Plot interesting genes an/or pathways for seurat3 object for given cell-phenodata variable and clustering resolution ### in dev ###

Usage: plot_igenes_seurat3.R --jobname=<value> [--mcafile=<file>] [--genefile=<file>] [--gmtfile=<file>] --specie=<value> --cellvar=<value> --resolution=<value> --infolder=<folder> --outfolder=<folder> --ncores=<value> 

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --jobname=<file>     Descriptive name for your experiment.
  --mcafile=<file>     Optional. If !is.null, do not use jobname to guess mcafilename.
  --genefile=<file>    One column text table with genes to plot. (ENSG.../ENMUSG...).
  --gmtfile=<file>     .gmt files with gene-sets to test together.
  --specie=<value>     Either hg or mm.
  --cellvar=<value>    colData variable to use as reference for sub-clustering.
  --resolution=<value> Clustering resolution.
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
	  'Rtsne', 'data.table')

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

cellvar <- arguments$cellvar
inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
jobname <- arguments$jobname
genefile <- arguments$genefile
res <- as.numeric(arguments$resolution)
# ------------------------------------------------------------------------------
# --- Loading mca onject
# ------------------------------------------------------------------------------

ncores <- as.numeric(arguments$ncores)
plan("multiprocess", workers = ncores)
plan()

if(is.null(arguments$gmtfile) & is.null(arguments$igenefile)) {

	stop("Please provide either an igenes or gmt file.")

} else if(!is.null(arguments$gmtfile)) {

	gmtfile <- arguments$gmtfile
	gmt <- CEMiTool::read_gmt(gmtfile)

} else if(!is.null(arguments$igenfile)) {
		
	genefile <- arguments$igenfile
	igenes <- fread(genefile)
	colnames(igenes)[1] <- "gene"
	igenes[[1]] <- toupper(igenes[[1]])
}

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

# ----- set ident -------
ref_cluster <- paste0("integrated_snn_res.", res)
Idents(object = mca) <- ref_cluster 
print(table(Idents(mca)))


# Here we go...
