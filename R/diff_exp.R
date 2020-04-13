#
# Estimate expression difference given a colData variable and a clustering resolution (Seurat3 objects only) ## In dev ##
#
 
" Estimate expression difference given a colData variable and a clustering resolution (Seurat3 objects only) ## In dev ##

Usage: diff_exp.R --jobname=<value> --groupvar=<value> --samplevar=<value> --resolution=<value> --specie=<value> --infolder=<folder> --outfolder=<folder>

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

# Main function for pseudo-bulk differential expression
pbulk_diffexp <- function(raw_eset, meta, sample_meta, jobname, mingenecount, refbiomart) { # count matrix plus metadata
	# Sum counts per sample
	ss <- unique(meta[[samplevar]])

	sample_pbulk <- Reduce(cbind,
			lapply(ss, function(s) {
				cells <- rownames(meta)[which(meta[[samplevar]] == s)]
				bulk <- data.frame(rowSums(raw_eset[, cells]))
				sname <- gsub("_filtered_FILTERED_qc_seurat3", "", s)
				colnames(bulk) <- sname
				print(head(bulk, 2))
				bulk
			})
		)
	
	print("Genes x samples matrix")
	print(head(sample_pbulk))
	print(dim(sample_pbulk))

	# Check sample expression distributions
	dist_data <- reshape2::melt(sample_pbulk)

	gg_raw_dist <- ggplot(dist_data, aes(`value`)) +
		geom_density(aes(group = `variable`, fill = `variable`), alpha = .4) +
		theme_classic()

	gg_raw_log_dist <- ggplot(dist_data, aes(log10(`value` + 1))) +
		geom_density(aes(group = `variable`, fill = `variable`), alpha = .4) +
		theme_classic()

	dds <- DESeqDataSetFromMatrix(countData = sample_pbulk[, rownames(sample_meta)],
				      colData = sample_meta,
				      design = ~ condition + age)

	print(dds)

	print(quantile(rowSums(counts(dds))))
	keep <- rowSums(counts(dds)) > mingenecount
	dds <- dds[keep, ]

	# Normalization for visualization ---------
	vsd <- vst(dds, blind=FALSE)
	dim(assay(vsd))
	head(assay(vsd), 3)
	vsd_mat <- assay(vsd)
	colnames(vsd_mat) <- colnames(sample_pbulk)
	vs <- colnames(sample_meta)
	pdf(paste0(outputfolder, "/", jobname, "_DESeq2_vsd_norm_metadata_PCA.pdf"))
	for(v in vs) {
		plot(
			autoplot(prcomp(t(vsd_mat)[sample_meta[['sample']], ]), 
				 data=data.frame(sample_meta), 
				 colour = v,
				 label = TRUE, 
				 label.size = 3) +
		     	theme_classic()
		)
	}
	dev.off()

	# Check sample distributions of this normalized expression
	dist_data <- reshape2::melt(vsd_mat)

	gg_vsd_dist <- ggplot(dist_data, aes(`value`)) +
		geom_density(aes(group = `Var2`, fill = `Var2`), alpha = .4) +
		theme_classic()

	pdf(paste0(outputfolder, "/", jobname, "_DESeq2_expression_sample_pbulk_distribution.pdf"))
		plot(gg_raw_dist)
		plot(gg_raw_log_dist)
		plot(gg_vsd_dist)
	dev.off()
	# -----------------------------------------

	dds <- DESeq(object = dds, 
		     test = "Wald",
		     quiet = TRUE,
		     parallel = FALSE, 
		     BPPARAM = bpparam())

	res <- results(dds, contrast=c("condition", "PD", "CO"))

	res[['gene']] <- rownames(res)

	res <- data.frame(res)

	# Check if stable gene ID was provided with or without version
	g2test <- res[['gene']][1]

	if(grepl("\\.\\d+$", g2test)) {
		res[['gene_version']] <- res[['gene']]
		res[['gene']] <- gsub("\\.\\d+$", "", res[['gene_version']])
	} else {
		res[['gene_version']] <- res[['gene']]
	}

	res <- merge(res, refbiomart, by.x = "gene", by.y = "Gene stable ID") 
	
	res[['index']] <- seq(nrow(res))

	confects <- normal_confects(as.numeric(res$log2FoldChange),
				    se = as.numeric(res$lfcSE), 
				    fdr = 0.05, 
				    full = TRUE)

	print(head(confects$table, 3))

	res <- merge(res, 
		     dplyr::select(confects$table, c(index, `rank`)),
		     by = 'index', all = TRUE)


	irank <- quantile(res[['rank']], metathr)
	print(irank)

	res %>%
		dplyr::mutate(perturbation = ifelse(`rank` <= irank, 
						    ifelse(log2FoldChange < 0, "Down", "Up"),
						    "Unperturbed")) -> res 

	write.table(res, paste0(outputfolder, "/", jobname, "_deg_DESeq2_pbulk_results.tsv"),
		    sep = "\t", quote = FALSE, row.names = FALSE)

	# Volcano plot
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

	pdf(paste0(outputfolder, "/", jobname, "_DESeq2_samplewise_pbulk_volcano.pdf"))
		plot(gg)
	dev.off()

	htmlwidgets::saveWidget(as_widget(ggplotly(gg)), 
            paste0(normalizePath(outputfolder), 
	           "/", jobname, "_DESeq2_volcanoplot.html"))
	# -

	# MA plot
	gg <- ggplot(dplyr::filter(res, !is.na(padj)), 
		     aes(x = baseMean, y = log2FoldChange, 
			 text = `Gene name`,
			 color = perturbation)) +
		geom_point(alpha = 0.5) +
		geom_errorbar(aes(ymin = log2FoldChange - lfcSE, 
				  ymax = log2FoldChange + lfcSE), 
			       alpha = .2) +
		scale_color_manual(values = c("blue", "grey", "red")) +
		theme_classic()

	pdf(paste0(outputfolder, "/", jobname, "_DESeq2_samplewise_pbulk_MAplot.pdf"))
		plot(gg)
	dev.off()

	htmlwidgets::saveWidget(as_widget(ggplotly(gg)), 
            paste0(normalizePath(outputfolder), 
	           "/", jobname, "_DESeq2_MAplot.html"))
	# -

	vis_eset <- reshape2::melt(vsd_mat)
	vis_eset <- merge(vis_eset, sample_meta, 
			  by.x = "Var2", by.y = "sample")

	vis_eset <- merge(vis_eset, 
			  select(res, c(gene_version, `Gene name`)), 
			  by.x = "Var1", by.y = "gene_version")

	top50 <- head(arrange(filter(res, rank <= irank), 
			      rank), 50)[['gene_version']]
	print(head(top50))

	dplyr::filter(data.frame(vis_eset), 
		      Var1 %in% top50) %>%
		ggplot(aes(x = condition, y = value,
			   fill = condition, color = condition)) +
		geom_violin(alpha = .2) +
		geom_jitter(alpha = .5) +
		theme_classic() +
		facet_wrap(~`Gene.name`, nrow = 10) -> ggex

	pdf(paste0(outputfolder, "/", jobname, "_DESeq2_samplewise_pbulk_top_50gene_violin.pdf"),
		width = 20, height = 20)
		plot(ggex)
	dev.off()
}

# all cells 
pbulk_diffexp(raw_eset, meta, sample_meta, jobname, mingenecount, refbiomart)

# For each cluster independently
cls <- unique(meta[[ident_col]])

lapply(cls, function(cl) {
	       print(cl)

	       sub_meta <- dplyr::filter(meta, !!rlang::sym(ident_col) == cl)
	       print(dim(sub_meta))

	       cl_cells <- sub_meta[['barcode']]
	       rownames(sub_meta) <- cl_cells
	       print(head(sub_meta, 2))
	       print(dim(sub_meta))
	       sub_raw_eset <- raw_eset[, cl_cells]
	       print(dim(sub_raw_eset))

	       sub_jobname <- paste0(jobname, "_cluster_", cl)
	       
	       pbulk_diffexp(sub_raw_eset, sub_meta, sample_meta, sub_jobname, mingenecount, refbiomart)

})

# ----------

# differential expression using edgeR

pbulk_edger_diffexp <- function(raw_eset, meta, sample_meta, jobname, mingenecount, refbiomart) {
	# Sum counts per sample
	ss <- unique(meta[[samplevar]])

	sample_pbulk <- Reduce(cbind,
			lapply(ss, function(s) {
				cells <- rownames(meta)[which(meta[[samplevar]] == s)]
				bulk <- data.frame(rowSums(raw_eset[, cells]))
				sname <- gsub("_filtered_FILTERED_qc_seurat3", "", s)
				colnames(bulk) <- sname
				print(head(bulk, 2))
				bulk
			})
		)
	
	print("Genes x samples matrix")
	print(head(sample_pbulk,2))
	print(dim(sample_pbulk))

 	y <- DGEList(counts = sample_pbulk[, rownames(sample_meta)],
				      samples = sample_meta)

	# Checking sample library sizes
	discarded <- scater::isOutlier(y$samples$lib.size, log=TRUE, type="lower")

	if(any(discarded)) {
		print("Notice that these samples were removed from the DE analysis:")
		print(colnames(y)[discarded])
		y <- y[, !discarded]
	} else {
		print("Sample library sizes are more/less homogeneus")
		print(y$samples$lib.size)
	}
	print("Discarded?:")
	summary(discarded)

	keep <- filterByExpr(y, group=sample_meta[[groupvar]],
			     min.total.count = 200)
	table(keep)

	y <- y[keep, ]
	
	y <- calcNormFactors(y)
	y$samples

	# The actual design matrix /should it be customed for every dataset?/
	design <- model.matrix(~ as.numeric(age) + factor(condition),
			       y$samples)

	print(`design`)
	y <- estimateDisp(y, design)
	print(summary(y$trended.dispersion))

	pdf(paste0(outputfolder, "/", jobname, "_biological_component_geneVariance_edgeR.pdf"))
		plotBCV(y)
	dev.off()
	
	fit <- glmQLFit(y, design, robust=TRUE)
	print(summary(fit$var.prior))

	pdf(paste0(outputfolder, "/", jobname, "_likelihood_dispersion_edgeR.pdf"))
		plotQLDisp(fit)
	dev.off()

	res <- glmQLFTest(fit, coef=ncol(design))
	summary(decideTests(res))
	topTags(res)

	de_res <- data.frame(topTags(res, nrow(res)))
	de_res[['gene']] <- rownames(de_res)

	# Check if stable gene ID was provided with or without version
	g2test <- de_res[['gene']][1]

	if(grepl("\\.\\d+$", g2test)) {
		de_res[['gene_version']] <- de_res[['gene']]
		de_res[['gene']] <- gsub("\\.\\d+$", "", de_res[['gene_version']])
	} else {
		de_res[['gene_version']] <- de_res[['gene']]
	}

	de_res <- merge(de_res, refbiomart, 
			by.x = "gene", by.y = "Gene stable ID") 

	print(head(de_res, 3))

	de_res %>%
		mutate(pertb = abs(logFC)*(-log10(PValue))) %>%
		arrange(-pertb) -> de_res

	de_res[['rank']] <- seq(nrow(de_res))
	print(head(de_res))
	
	metathr <- 0.05
	irank <- de_res[round(nrow(de_res)*metathr), ][['rank']]

	de_res %>%
		dplyr::mutate(perturbation = ifelse(`rank` <= irank, 
						    ifelse(logFC < 0, "Down", "Up"),
						    "Unperturbed")) -> de_res 

	write.table(de_res, paste0(outputfolder, "/", jobname, "_deg_edgeR_pbulk_results.tsv"),
		    sep = "\t", quote = FALSE, row.names = FALSE)

	gg <- ggplot(de_res, 
		     aes(x = logFC, y = -log10(PValue), 
			 text = `Gene name`,
			 color = perturbation)) +
		geom_point(alpha = 0.5) +
		scale_color_manual(values = c("blue", "grey", "red")) +
		theme_classic()

	pdf(paste0(outputfolder, "/", jobname, "_edgeR_samplewise_pbulk_volcano.pdf"))
		plot(gg)
	dev.off()

	htmlwidgets::saveWidget(as_widget(ggplotly(gg)), 
            paste0(normalizePath(outputfolder), 
	           "/", jobname, "_edgeR_volcanoplot.html"))

	vis_eset <- data.frame(res$fitted.values)
	vis_eset[['gene']] <- rownames(vis_eset)
	vis_eset <- reshape2::melt(vis_eset)
	vis_eset <- merge(vis_eset, sample_meta, 
			  by.x = "variable", by.y = "sample")

	vis_eset <- merge(vis_eset, 
			  select(de_res, c(gene_version, `Gene name`)), 
			  by.x = "gene", by.y = "gene_version")

	top50 <- head(arrange(de_res, rank), 50)[['gene_version']]
	print(head(top50))

	dplyr::filter(data.frame(vis_eset), 
		      gene %in% top50) %>%
		ggplot(aes(x = condition, y = log10(value),
			   fill = condition, color = condition)) +
		geom_violin(alpha = .2) +
		geom_jitter(alpha = .5) +
		theme_classic() +
		facet_wrap(~`Gene.name`, nrow = 10) -> ggex

	pdf(paste0(outputfolder, "/", jobname, "_edgeR_samplewise_pbulk_top_50gene_violin.pdf"),
		width = 20, height = 20)
		plot(ggex)
	dev.off()
}

# all cells 
pbulk_edger_diffexp(raw_eset, meta, sample_meta, jobname, mingenecount, refbiomart)

# For each cluster independently
cls <- unique(meta[[ident_col]])

lapply(cls, function(cl) {
	       print(cl)

	       sub_meta <- dplyr::filter(meta, !!rlang::sym(ident_col) == cl)
	       print(dim(sub_meta))

	       cl_cells <- sub_meta[['barcode']]
	       rownames(sub_meta) <- cl_cells
	       print(head(sub_meta, 2))
	       print(dim(sub_meta))
	       sub_raw_eset <- raw_eset[, cl_cells]
	       print(dim(sub_raw_eset))

	       sub_jobname <- paste0(jobname, "_cluster_", cl)
	       
	       pbulk_edger_diffexp(sub_raw_eset, sub_meta, sample_meta, sub_jobname, mingenecount, refbiomart)

})


# Comparison edgeR AND DESeq2
deseq2 <- data.table::fread(paste0(outputfolder, "/", jobname, "_deg_DESeq2_pbulk_results.tsv"))
edger <- data.table::fread(paste0(outputfolder, "/", jobname, "_deg_edgeR_pbulk_results.tsv"))

r <- merge(deseq2, edger, by = 'gene_version')
print(head(r))

r %>%
	mutate(perturbation = ifelse(`perturbation.x` != "Unperturbed" & `perturbation.y` != "Unperturbed",
				     "Perturbed", "Unperturbed")) -> r 

gg <- ggplot(r, aes(x = log2FoldChange, y = logFC, 
		    text = `Gene name.x`,
		    color = `perturbation`)) +
	geom_point(alpha = .4) +
#	geom_density_2d(alpha = .5) +
	theme_classic()

pdf(paste0(outputfolder, "/", jobname, "_corr_edgeR_DESeq2_logFC_samplewise_pbulk.pdf"))
	plot(gg)
dev.off()

htmlwidgets::saveWidget(as_widget(ggplotly(gg)), 
           paste0(normalizePath(outputfolder), 
           "/", jobname, "_corr_edgeR_DESeq2_logFC_samplewise_pbulk.html"))

# for all clusters
lapply(cls, function(cl) {
		print(cl)
		sub_jobname <- paste0(jobname, "_cluster_", cl)
		deseq2 <- data.table::fread(paste0(outputfolder, "/", sub_jobname, "_deg_DESeq2_pbulk_results.tsv"))
		edger <- data.table::fread(paste0(outputfolder, "/", sub_jobname, "_deg_edgeR_pbulk_results.tsv"))

		r <- merge(deseq2, edger, by = 'gene_version')
		print(head(r))

		r %>%
			mutate(perturbation = ifelse(`perturbation.x` != "Unperturbed" & `perturbation.y` != "Unperturbed",
					     "Perturbed", "Unperturbed")) -> r 

		gg <- ggplot(r, aes(x = log2FoldChange, y = logFC, 
			    text = `Gene name.x`,
			    color = `perturbation`)) +
			geom_point(alpha = .4) +
		#	geom_density_2d(alpha = .5) +
			theme_classic()

		pdf(paste0(outputfolder, "/", sub_jobname, "_corr_edgeR_DESeq2_logFC_samplewise_pbulk.pdf"))
			plot(gg)
		dev.off()

		htmlwidgets::saveWidget(as_widget(ggplotly(gg)), 
		           paste0(normalizePath(outputfolder), 
		           "/", sub_jobname, "_corr_edgeR_DESeq2_logFC_samplewise_pbulk.html"))
})
