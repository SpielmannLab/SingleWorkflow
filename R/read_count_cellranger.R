#
# Read 10X - CellRanger output and write and SingleCellExperiment object
#

"CellRanger output to cell_data_set (Monocle3)

Usage: read_count_cellranger.R --infolder=<folder> --outfolder=<folder> --jobname=<value> --genome=<value>

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --infolder=<file>    Path to ../outs parent folder.
  --outfolder=<file>   Path to results folder.
  --jobname=<value>    Descriptive name.
  --genome=<value>     Genome string eg. mm10.

Author: 
  Cesar Prada
e-mail: 
  prada@molgen.mpg.de

"-> doc

library(docopt)
arguments<-docopt(doc, quoted_args=TRUE)
print(arguments)

# --- Dependencies

pkgs <- c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
	  'limma', 'S4Vectors', 'SingleCellExperiment',
	  'SummarizedExperiment', 'batchelor', 'devtools')

if (!requireNamespace("BiocManager", quietly = TRUE)) {
	install.packages("BiocManager")
	BiocManager::install()
} else {
	if (!requireNamespace(pkgs, quietly = TRUE)) {
		BiocManager::install(pkgs)
	} else {
		message("Cool! your machine has everything is needed")
	}
}

# Install last version of both (Beta version keep changing)

#devtools::install_github('cole-trapnell-lab/leidenbase')
#devtools::install_github('cole-trapnell-lab/monocle3')

#
supressMessages(library(monocle3))
#

# --- Parameters

inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
jobname <- arguments$jobname
`genome` <- arguments$genome
umi_cutoff <- 0

# --- Run

cds <- load_cellranger_data(inputfolder)

cds_file <- paste0(outputfolder, '/', jobname, '.rds')

saveRDS(cds, file = cds_file)

message(paste("The cell_data_set object (Monocle3)", 
	      cds_file, "has been created."))
