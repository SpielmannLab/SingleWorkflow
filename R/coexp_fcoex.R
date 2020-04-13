#
# Coexpression analysis based on entropy and FCBF. Input magic imputed data matrix object. ## In dev ##
#

"  Coexpression analysis based on entropy and FCBF. Input imputed data matrix magic object. ## In dev ##

Usage: coexp_fcoex.R --jobname=<value> --mcafile --infolder=<folder> --outfolder=<folder> --ncores=<value> 

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --jobname=<file>     Descriptive name for your experiment.
  --mcafile=<file>     mca.rds file
  --colgroup=<value>   additional color labor group to include for reference 
  --infolder=<file>    Path to the single_cell_data .rds files.  
  --ncores=<value>     # of processors to use
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
#devtools::install_github('lubianat/fcoex')


# -------------------------------
suppressMessages(library(monocle3))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(Seurat))
#suppressMessages(library(batchelor))
suppressMessages(library(sctransform))
suppressMessages(library(future))
suppressMessages(library(fcoex))
#suppressMessages(library(dplyr))
#suppressMessages(library(ggwordcloud))
# ------------------------------


inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
jobname <- arguments$jobname
npcs <- arguments$npcs
# ---------------------------------------------------
# Setting multicore env for some Seurat applications
# ---------------------------------------------------
ncores <- as.numeric(arguments$ncores)
plan("multiprocess", workers = ncores)
plan()


## dummy example from the package
library(SingleCellExperiment) 
data("mini_pbmc3k")
targets <- colData(mini_pbmc3k)$clusters
head(targets)
exprs <- as.data.frame(assay(mini_pbmc3k, "logcounts"))
exprs[1:10, 1:10]
dim(exprs)
fc <- new_fcoex(exprs, targets)
fc <- discretize(fc)
fc <- find_cbf_modules(fc)
fc 
##


# --- Read cell_data_set object with merged samples and cell_metadata
mca_file <- paste0(outputfolder, '/', jobname, '_seurat3_merged_clustered_resrange.rds')

mca <- readRDS(mca_file)

message(paste("The cell_data_set object (Seurat3)", 
	      mca_file, "has been read."))

# --- Read .rds object with magic results. This must be produced by the using the same provided mca_file 

magic_file <- paste0(outputfolder, '/', jobname, '_magic_output.rds')
magic <- readRDS(magic_file)
message(paste("The .rds object ", 
	      magic_file, "with inpute_magic.R output has been read."))

exprs <- data.frame(t(magic$magic$result))
exprs <- exprs[colnames(magic$corr), ]
rownames(exprs) <- gsub("\\.\\d+$", "", rownames(exprs))
colnames(exprs) <- paste0("V", seq(ncol(exprs))) 
print(exprs[1:5, 1:10])
print(dim(exprs))

# --- Get cell metadata
target <- mca@meta.data
print(target[1:3, c(1:5, ncol(target))])
print(dim(target))
print(nrow(target) == ncol(exprs))
# - set cluster reference
ident_col <- colnames(target)[grep(paste0("_snn_res.", res, "$"), colnames(target))]
print(ident_col)
targets <- as.factor(target[[ident_col]])
print(head(targets))

fc <- new_fcoex(`exprs`, targets)

fc <- discretize(fc, number_of_bins = 5)

fc <- find_cbf_modules(fc,
		       n_genes = 50,
		       verbose = TRUE,
		       is_parallel = FALSE)

fc <- get_nets(fc)

save_plots(name = paste0("fcoex_", jobname), fc, 
	   force = TRUE, directory = paste0(outputfolder, "/"))

mods <- fc@module_list 


for(coex_mod in names(mods)) {
	genes <- mods[[coex_mod]]
	jpeg(paste0(outputfolder, '/', jobname, "_", "fcoex_mod_", coex_mod, "_umap.jpeg"), 
		    width = 1980, height = 1900, pointsize = 74, quality = 100)
	plot(
	     FeaturePlot(mca, 
	    	reduction = "umap", 
	    	pt.size = 0.01, 
	    	features=genes) +
			ggplot2::theme(legend.position = "none") + 
			ggtitle(coex_mod) +
			scale_colour_gradient(low = "grey", high = "blue")
	)
	dev.off()
}

saveRDS(fc, paste0(outputfolder, "/", jobname, "_fcoex_object.rds"))

modtab <- Reduce(rbind,
		 lapply(names(fc@module_list), function(mod) {
				data.frame('gene' = fc@module_list[[mod]],
					   'mod' = mod)
				})
		 )

head(modtab)
write.table(modtab, paste0(outputfolder, "/", jobname, "_fcoex_modules.tsv"),
	    sep = "\t", quote = FALSE, row.names = FALSE)
