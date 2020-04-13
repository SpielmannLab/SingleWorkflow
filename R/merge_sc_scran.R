#
# Read and merge several cell_data_set (scran) objects following the scran approach ## In dev ##
#

" Read and merge several cell_data_set (scran) objects following the Seurat approach ## In dev ##

Usage: merge_sc_scran.R --jobname=<value> --infolder=<folder> --outfolder=<folder> --ncores=<value> [<file> <file>...]

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --jobname=<file>     Descriptive name for your experiment.
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

# ------- Functions
cds2mca <- function(cds) {
	mca <- CreateSeuratObject(counts = exprs(cds),
			  meta.data = data.frame(pData(cds)),
			  project = jobname)
	mca@assays$RNA@meta.features <- data.frame(fData(cds))
	mca
}

# -----------------


# --------------
# Main function
# --------------

cdsfiles <- arguments$file
inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
jobname <- arguments$jobname



#sce <- lapply(cdss, cds2sce)


# Commom features ----
#cfeat <- unlist(lapply(sce, rownames))
#cfeat <- table(cfeat)
#cfeat <- names(which(cfeat == length(sce)))
# -------------------

#sce <- lapply(sce, function(x) x[cfeat, ])
#lapply(sce, dim)

########
#rescaled <- multiBatchNorm(sce[[1]], sce[[2]])
#rescaled <- multiBatchNorm(sapply(sce, function(x) x))
########


# ------------------
# Marionis' approach
# ------------------
#library(scater)
#library(scran)

sce <- cds2sce(cds)
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters=clusters)
print("\n\n\nSize Factor summary:")
print(summary(sizeFactors(sce)))
sce <- logNormCounts(sce)

sce <- runPCA(sce, ntop = 2000, 
	      subset_row = NULL)

sce <- runTSNE(sce, ntop = 2000,
	       dimred = "PCA")

sce <- runUMAP(sce, ntop = 2000,
	       dimred = "PCA")

sce <- runTSNE(sce, ntop = 2000,
	       dimred = "PCA")

g_pca <- plotReducedDim(sce, dimred="PCA", colour_by="sample", ncomponents=4)
g_pca <- plotReducedDim(sce, dimred="PCA", colour_by="sample")

g_tsne <- plotReducedDim(sce, dimred="TSNE", colour_by="sample")

g_umap <- plotReducedDim(sce, dimred="UMAP", colour_by="sample")


pdf(paste0(outputfolder, "/", jobname, "_scran_preprocessing.pdf"), width = 16, height = 8)
	plot_grid(plot_grid(g_pca, g_tsne, g_umap, ncol = 2))	
dev.off()


# ---------------------
# Global variant genes
# ---------------------

dec <- modelGeneVar(sce, BPPARAM = SerialParam(TRUE))
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)

# Variance by samples
dec3 <- modelGeneVar(sce, block=sce$sample)
per.block <- dec3$per.block
par(mfrow=c(3, 2))
for (i in seq_along(per.block)) {
    decX <- per.block[[i]]
    plot(decX$mean, decX$total, xlab="Mean log-expression", 
        ylab="Variance", main=names(per.block)[i])
    curve(metadata(decX)$trend(x), col="blue", add=TRUE)
}

# Get the top 10% of genes.
top.hvgs <- getTopHVGs(dec, prop=0.1)

# Get the top 2000 genes.
top.hvgs2 <- getTopHVGs(dec, n=2000)

# Get all genes with positive biological components.
top.hvgs3 <- getTopHVGs(dec, var.threshold=0)

# Get all genes with FDR below 5%.
top.hvgs4 <- getTopHVGs(dec, fdr.threshold=0.05)


