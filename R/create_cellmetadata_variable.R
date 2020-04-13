#
# Add cell phenodata variables (project specific) 
#

"  Add cell phenodata variables (project specific) 

Usage: create_cellmetadata.R --rdsfile=<value> --infolder=<folder> --outfolder=<folder> 

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --rdsfile=<file>
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

# -------------------------------
suppressMessages(library(monocle3))
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
# ------------------------------

inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
rdsfile <- arguments$rdsfile

# ----------------------------------------------------
#  Read cell_data_set object with cluster metadata
# ----------------------------------------------------

mca_file <- paste0(inputfolder, '/', rdsfile)

mca <- readRDS(mca_file)

message(paste("The cell_data_set object (Seurat3)", 
	      mca_file, "has been read."))

# Including meta.data relevant variables
# <<<<<<<<  Specific for the PD project
#s <- mca@meta.data$sample
#mca@meta.data$sex <- unlist(regmatches(s, regexec("_M_|_F_", s)))
#mca@meta.data$age <- unlist(regmatches(s, regexec("_\\d+_|_\\d+_", s)))
#mca@meta.data$condition <- unlist(regmatches(s, regexec("^CO|^PD", s)))
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


## <<<<<<<< Specific for the L_B6_project
#s <- mca@meta.data$sample
#mca@meta.data$limb <- unlist(regmatches(s, regexec("FLs|HLs", s)))
#mca@meta.data$stage <- unlist(regmatches(s, regexec("E\\d+_\\d", s)))
#mca@meta.data$sample <- unlist(regmatches(s, regexec("E\\d+_\\dBlack6HLs|E\\d+_\\dBlack6FLs", s)))
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

print(mca)

mca_file <- paste0(outputfolder, '/', rdsfile)
saveRDS(mca, file = mca_file)
message(paste("The cell_data_set object (Seurat3)", 
	      mca_file, "has been created."))
