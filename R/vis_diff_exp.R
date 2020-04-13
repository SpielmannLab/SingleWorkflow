#
# Visualize differential expression results # in dev #
#

" Visualize differential expression results # in dev #

Usage: vis_diff_exp.R --jobname=<value> --groupvar=<value> --samplevar=<value> --resolution=<value> --specie=<value> --infolder=<folder> --outfolder=<folder>

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --jobname=<file>     Descriptive name for your experiment.
  --groupvar=<value>   Variable to do comparison.
  --samplevar=<value>  Sample variable. Ref for pseudo-bulk.
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
	  'Rtsne', 'MASS', 'htmlwidgets')

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
suppressMessages(library(zinbwave))
suppressMessages(library(DESeq2))
suppressMessages(library(ggfortify))
suppressMessages(library(topconfects))
suppressMessages(library(plotly))
suppressMessages(library(edgeR))

#suppressMessages(library(dplyr))
#suppressMessages(library(ggwordcloud))
# ------------------------------

jobname <- arguments$jobname
groupvar <- arguments$groupvar
samplevar <- arguments$samplevar
res <- arguments$resolution
specie <- arguments$specie
inputfolder <- arguments$infolder
outputfolder <- arguments$outfolder

# --- Read cell_data_set object with merged samples
mca_file <- paste0(outputfolder, '/', jobname, '_seurat3_merged_clustered_resrange.rds')

mca <- readRDS(mca_file)

message(paste("The cell_data_set object (Seurat3)", 
	      mca_file, "has been read."))

# Setting default cell identity based on clustering resolution
ident_col <- colnames(mca@meta.data)[grep(paste0("_snn_res.", res), colnames(mca@meta.data))]
print(table(mca@meta.data[[ident_col]]))
Idents(mca) <- ident_col

# Checking cell metadata variable
print(table(mca@meta.data[[groupvar]]))

# Checking sample metadata
print(table(mca@meta.data[[samplevar]]))

# Gene annotation
# Gene annotation
if(specie == "hg") {

	refbiomart <- data.table::fread("<human-biomart-reference-file>")

} else if(specie == "mm"){

	refbiomart <- data.table::fread("<mouse-biomart-reference-file>")

} else {

	stop("Sorry this specie isnt supported, try mm or hg please")

}


# Counts
raw_eset <- mca@assays$RNA@counts 
print(raw_eset[1:5, 1:10])
print(dim(raw_eset))

# Sample metadata
meta <- mca@meta.data
print(head(meta,2))
meta %>%
	mutate('sample' = gsub("_filtered_FILTERED_qc_seurat3", "", `sample`)) %>%
	group_by(!!rlang::sym(samplevar)) %>%
	summarise('sex' = gsub("_", "", unique(sex)),
		  'age' = as.numeric(gsub("_", "", unique(age))),
		  'condition' = unique(condition)) %>%
	data.frame -> sample_meta

sample_meta <- DataFrame(sample_meta)
sample_meta[[groupvar]] <- as.factor(sample_meta[[groupvar]])
rownames(sample_meta) <- sample_meta[[samplevar]] 
print(head(sample_meta,2))

# Cell metadata
meta <- mca@meta.data
meta[['barcode']] <- rownames(meta)
print(head(meta,2))

# Minimum number of UMI counts to consider in the DE analysis
# (UMI counts per sample) 
mingenecount <- 200

# Top X% perturbed genes
metathr <- 0.05

# DE results files pattern
method <- "DESeq2" # edgeR
de_pattern <- paste0("_deg_", method, "_pbulk_results.tsv")
de_files <- list.files(path = inputfolder, pattern = de_pattern)
de_files <- setNames(de_files, gsub(de_pattern, "", de_files))
de_files


des <- lapply(de_files, function(f) data.table::fread(paste0(inputfolder, '/', f)))

# Volcano plot
p_thr <- 0.1

lapply(names(des), function(res_name) {
	res <- des[[res_name]]
	res %>%                                                                                                                                                                                                                   
                dplyr::mutate(perturbation = ifelse(padj <= p_thr,                                                                                                                                                              
                                                    ifelse(log2FoldChange < 0, "Down", "Up"),                                                                                                                                     
                                                    "Unperturbed")) -> res
	
	gg <- ggplot(dplyr::filter(res, !is.na(padj)),

                     aes(x = log2FoldChange, y = -log10(padj),                                                                                                                                                                    
                         text = `Gene name`,                                                                                                                                                                                      
                         color = perturbation)) +

                geom_point(alpha = 0.5) +                                                                                                                                                                                         
                geom_errorbarh(aes(xmin = log2FoldChange - lfcSE,                                                                                                                                                                 
                                   xmax = log2FoldChange + lfcSE),                                                                                                                                                                
                               alpha = .2) +

                scale_color_manual(values = c("blue", "grey", "red")) +                                                                                                                                                           
                theme_classic()                                                                                                                                                                                                   
                                                                                                                                                                                                                                  
        pdf(paste0(outputfolder, "/", jobname, "_", res_name, "_DESeq2_samplewise_pbulk_volcano_perturbed.pdf"))                                                                                                                                           
                plot(gg)                                                                                                                                                                                                          
        dev.off()                                                                                                                                                                                                                 
                                                                                                                                                                                                                                  
        htmlwidgets::saveWidget(as_widget(ggplotly(gg)),                                                                                                                                                                          
            paste0(normalizePath(outputfolder),
                   "/", jobname, "_", res_name, "_DESeq2_volcanoplot_perturbed.html"))
})


de_tab <- Reduce(rbind,
		 lapply(names(des), function(nd) {
			 de_res <- des[[nd]]
			 de_res[['cluster']] <- nd
			 de_res
		  })
		 )
			 
dim(de_tab)
head(de_tab)

p_thr <- 0.1

de_tab %>%
	group_by(`cluster`) %>%
	summarise('up' = length(which(log2FoldChange > 0 & padj <= p_thr)),
		  'dw' = length(which(log2FoldChange < 0 & padj <= p_thr))) -> de_n

dim(de_n)
de_n <- reshape2::melt(de_n)
de_n[['value']] <- de_n[['value']]*ifelse(de_n[['variable']]=="dw", -1, 1)

gg_degs <- ggplot(de_n, aes(x = `cluster`, y = `value`, label = `value`)) +
		geom_point(aes(color = `variable`, size = `value`)) +
		geom_segment(aes(x=`cluster`, xend = `cluster`, y=0, yend=`value`), color="grey", alpha = 0.5) +
		scale_color_manual(values = c("#b2182b", "#2166ac")) + 
		scale_x_discrete(limits = rev(unique(de_n[['cluster']]))) +
		coord_flip() +
		geom_text() +
		theme_classic()

pdf(paste0(outputfolder, "/", jobname, "_deg_result_", method, "_n_degs.pdf"),
    width = 4.5, height = 3)
	plot(gg_degs)
dev.off()
