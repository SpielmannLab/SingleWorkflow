#
# Correlate interesting genes, raw, normalized and magic inputed ## In dev ##
#

" Correlate interesting genes, raw, normalized and magic inputed ## In dev ##

Usage: magic_corr_igenes_seurat3.R --jobname=<value> [--mcafile=<file>] [--genefile=<file>] [--gmtfile=<file>] --specie=<value> --cellvar=<value> --resolution=<value> --infolder=<folder> --outfolder=<folder> --ncores=<value> 

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
specie <- arguments$specie
gmtfile <- arguments$gmtfile
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
	colnames(igenes) <- "gene"
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


# ------------------------------------------------------------------------------
# Marker expression visualization (umap)
# ------------------------------------------------------------------------------

# Gene annotation
if(specie == "hg") {

	refbiomart <- data.table::fread("<human-biomart-reference-file>")

} else if(specie == "mm"){

	refbiomart <- data.table::fread("<mouse-biomart-reference-file>")

} else {

	stop("Sorry this specie isnt supported, try mm or hg please")

}


if(!is.null(arguments$igenfile)) {

	if(!grepl("ENSG|ENSMUSG", igenes[[1]][1])) {
		mca.markers <- merge(igenes, refbiomart, by.x = "gene", by.y = "Gene name", all.x = TRUE) %>%
			filter(!is.na(`Gene stable ID`))
		mca.markers[['Gene name']] <- mca.markers[['gene']]
		mca.markers[['gene']] <- mca.markers[['Gene stable ID']]
		g <- mca.markers[['Gene stable ID']]

	} else {
		mca.markers <- merge(igenes, refbiomart, by.x = "gene", by.y = "Gene stable ID", all.x = TRUE)
		mca.markers[['Gene name']][which(is.na(mca.markers[['Gene name']]))] <- mca.markers[['gene']][which(is.na(mca.markers[['Gene name']]))]
		g <- mca.markers[['gene']]
	}

	# SCT normalized expression 
	m <- data.matrix(mca@assays$SCT@data)
	g <- g[g %in% rownames(m)]
	print(g)
	sub_mat <- data.frame(t(m[g, ]))
	sub_mat[['barcode']] <- rownames(sub_mat)
	umap <- data.frame(mca@reductions$umap@cell.embeddings)
	umap[['barcode']] <- rownames(umap)

	sub_umap_eset <- merge(sub_mat, umap, by = "barcode")

	ggs <- lapply(g[idx], function(igene) {
		print(igene)      
		ggplot(sub_umap_eset, 
			     aes(x = UMAP_1, y = UMAP_2, color = !!rlang::sym(igene))) +
		      	geom_point(shape = ".") +
			scale_color_gradient2(low="grey80", 
					      high="blue", 
					      na.value = "grey80") +
			ggtitle(paste(dplyr::filter(mca.markers, 
						    gene==igene)[['Gene name']], 
				      collapse = '///')) +
			theme_classic()
	})


	title <- ggdraw() + 
		draw_label(
			paste0(jobname, " interesting genes"),
			fontface = 'bold',
			hjust = 0.5
		) +
	theme(
		plot.margin = margin(0, 0, 0, 0)
	)
	jpeg(paste0(outputfolder, '/', jobname, "_cool_marker_genes_names_SCT_norm_values.jpeg"), 
	    width = 980, height = 780, pointsize = 74, quality = 100)
		plot(
		    plot_grid(title,
			      plot_grid(plotlist = ggs, ncol = 3),
			      ncol = 1, rel_heights = c(0.01, 1))
		)
	dev.off()
	pdf(paste0(outputfolder, '/', jobname, "_cool_marker_genes_names_SCT_norm_values.pdf"), 
	    width = 14, height = 6)
		plot(
		    plot_grid(title,
			      plot_grid(plotlist = ggs, ncol = 3),
			      ncol = 1, rel_heights = c(0.01, 1))
		)
	dev.off()

	# violing plots
 	s <- mca@meta.data$sample
	mca@meta.data$stage <- unlist(regmatches(s, regexec("E\\d+_\\d", s)))
	meta <- data.frame(mca@meta.data)
	meta[['barcode']] <- rownames(meta)
	sub_mat <- data.frame(t(m[g, ]))





	sub_mat[['barcode']] <- rownames(sub_mat)
	

	g <- ggplot(sub_mat, aes(x = ENSMUSG00000058665, y = ENSMUSG00000099336)) +
		geom_jitter(shape = ".") +
		theme_classic()

	pdf(paste0(outputfolder, '/', jobname, "_correlation_ENSMUSG00000058665_ENSMUSG00000099336.pdf"), 
	    width = 4, height = 4)
		plot(
		     g
		)
	dev.off()

	e95 <- dplyr::filter(meta, stage == "E9_5")[['barcode']]

	g <- ggplot(sub_mat[e95,], aes(x = ENSMUSG00000058665, y = ENSMUSG00000099336)) +
		geom_jitter(shape = ".") +
		theme_classic()

	pdf(paste0(outputfolder, '/', jobname, "_correlation_ENSMUSG00000058665_ENSMUSG00000099336.pdf"), 
	    width = 4, height = 4)
		plot(
		     g
		)
	dev.off()

	e95_ect <- dplyr::filter(meta, stage == "E9_5" & integrated_snn_res.0.1 == 2)[['barcode']]

	g <- ggplot(sub_mat[e95_ect,], aes(x = ENSMUSG00000058665, y = ENSMUSG00000099336)) +
		geom_jitter(shape = ".") +
		theme_classic()

	pdf(paste0(outputfolder, '/', jobname, "_correlation_ENSMUSG00000058665_ENSMUSG00000099336.pdf"), 
	    width = 4, height = 4)
		plot(
		     g
		)
	dev.off()

	sub_mat <- merge(sub_mat, meta, by="barcode")

	g <- ggplot(sub_mat, aes(x = ENSMUSG00000058665, y = ENSMUSG00000099336)) +
		geom_jitter(shape = ".") +
		theme_classic() +
		facet_wrap(~stage*integrated_snn_res.0.1, ncol = 6)

	pdf(paste0(outputfolder, '/', jobname, "_correlation_ENSMUSG00000058665_ENSMUSG00000099336.pdf"), 
	    width = 14, height = 14)
		plot(
		     g
		)
	dev.off()



	gene_sub_mat <- melt(sub_mat, id.vars='barcode')
	gene_sub_mat <- merge(gene_sub_mat, 
			      select(meta, c(barcode, 
					     !!rlang::sym(ref_cluster),
					     !!rlang::sym(cellvar))), 
				by = 'barcode') %>%
					group_by(variable, 
						 !!rlang::sym(ref_cluster), 
						 !!rlang::sym(cellvar)) %>%
					summarize(n = length(which(`value` > 0)),
						  sum_exp = sum(`value`)) %>%
					mutate(mean_exp_cell = sum_exp/n)
			
	head(gene_sub_mat)
	dim(gene_sub_mat)


	gene_sub_mat <- merge(gene_sub_mat, refbiomart, 
			      by.x = "variable", 
			      by.y = "Gene stable ID", all.x = TRUE)

	vln_plot <- ggplot(gene_sub_mat, aes(x = stage, y = n, fill = integrated_snn_res.0.1)) +
		geom_bar(stat="identity", position = "dodge2") +
		theme_classic() +
		facet_wrap(~variable)

	vln_plot <- ggplot(gene_sub_mat, aes(x = stage, y = n, fill = integrated_snn_res.0.1)) +
		geom_bar(stat="identity", position = "dodge2") +
		theme_classic() +
		facet_wrap(~variable)

	umap <- data.frame(mca@reductions$umap@cell.embeddings)
	umap[['barcode']] <- rownames(umap)


	sub_umap_eset <- merge(sub_mat, umap, by = "barcode") %>%
		select(matches("barcode|ENSMUSG00000099336|UMAP_"))
	
	p<-"lila_lncRNA"
	sub_umap_eset[[p]] <- rowSums(select(sub_umap_eset, matches("ENSG")))
	sub_umap_eset[[paste0('mean_',p)]] <- rowMeans(select(sub_umap_eset, matches("ENSG")))

	ss_umap_eset <- select(sub_umap_eset, matches(paste0("barcode|UMAP|", p)))
	
	ss_umap_eset <- merge(ss_umap_eset, meta, by = 'barcode')
	
	aov_res <- aov(as.formula(paste0('`', p, '`', '~', ref_cluster)), data = ss_umap_eset)
			
	vln_plot <- ggplot(ss_umap_eset, 
			   aes(x = !!rlang::sym(cellvar), 
			       y = !!rlang::sym(p))) +
			geom_violin(alpha = 0.4) +
			geom_jitter(alpha = .1, shape=".") +
			theme_classic() +
			facet_wrap(as.formula(paste0('~',
						     ref_cluster)),
				nrow=1)

	jpeg(paste0(outputfolder, '/', jobname, "_marker_genes_names_SCT_norm_values_violin.jpeg"), 
	    width = 1980, height = 1980, pointsize = 74, quality = 100)
		plot(
		    plot_grid(title,
			      plot_grid(plotlist = list(vln_plot), ncol = 1),
			      ncol = 1, rel_heights = c(0.01, 1))
		)
	dev.off()


} 

keyword <- "Cytokine|Lysosome|Mitophagy|NOD-like|Cytosolic|Dopamine|Parkinson"
act_quantile <- 0.50

if(!is.null(arguments$gmtfile)) {

	# SCT normalized expression 
	m <- data.matrix(mca@assays$SCT@data)
	meta <- data.frame(mca@meta.data)
	meta[['barcode']] <- rownames(meta)
	# -------------------------

	print(head(gmt))
		
	mca.markers <- merge(gmt, refbiomart, by.x = "gene", by.y = "Gene name")

	print(head(mca.markers))

	pathways <- split(mca.markers, mca.markers[['term']])

	sub_paths <- pathways[grep(keyword, names(pathways))]
	names(sub_paths) <- gsub(' |\\.|\\(|\\)|\\-', "_", names(sub_paths))

	print(names(sub_paths))
	
	ggs <- lapply(names(sub_paths), function(p) {

		pathway <- sub_paths[[p]]
		g <- unique(pathway[['Gene stable ID']])
		g <- g[g %in% rownames(m)]

		if(length(g) < 2) {

			return(NULL)

		} else {

			sub_mat <- data.frame(t(m[g, ]))
			sub_mat[['barcode']] <- rownames(sub_mat)
			gene_sub_mat <- melt(sub_mat, id.vars='barcode')
	    		gene_sub_mat <- merge(gene_sub_mat, 
					      select(meta, c(barcode, 
							     !!rlang::sym(ref_cluster),
							     !!rlang::sym(cellvar))), 
					      by = 'barcode') %>%
				group_by(variable, !!rlang::sym(ref_cluster), !!rlang::sym(cellvar)) %>%
				summarize(n = length(which(`value` > 0)),
					  sum_exp = sum(`value`)) %>%
				mutate(mean_exp_cell = sum_exp/n)
			
			head(gene_sub_mat)

			gene_sub_mat <- merge(gene_sub_mat, refbiomart, 
					      by.x = "variable", 
					      by.y = "Gene stable ID", all.x = TRUE)

			g_wc <- ggplot(gene_sub_mat,
				       aes(label = `Gene name`,
					   size = mean_exp_cell,
					   x =  !!rlang::sym(cellvar),
					   color = !!rlang::sym(cellvar))) +
		      		ggwordcloud::geom_text_wordcloud(rm_outside = TRUE) +
				scale_size_area(max_size = 2.5) +
				scale_color_brewer(palette = "Dark2") +
				#scale_color_gradient(low = "lightgrey", high = "black") +
				theme_minimal() +
				facet_wrap(as.formula(paste0('~', ref_cluster)),
						nrow=1) +
		      		ggtitle("Gene expression in the pathway")
		
			umap <- data.frame(mca@reductions$umap@cell.embeddings)
			umap[['barcode']] <- rownames(umap)
	
			sub_umap_eset <- merge(sub_mat, umap, by = "barcode")
			sub_umap_eset[[p]] <- rowSums(select(sub_umap_eset, matches("ENSG")))
			sub_umap_eset[[paste0('mean_',p)]] <- rowMeans(select(sub_umap_eset, matches("ENSG")))

			ss_umap_eset <- select(sub_umap_eset, matches(paste0("barcode|UMAP|", p)))
		
			ss_umap_eset <- merge(ss_umap_eset, meta, by = 'barcode')
	
			aov_res <- aov(as.formula(paste0('`', p, '`', '~', ref_cluster)), data = ss_umap_eset)
			
			vln_plot <- ggplot(ss_umap_eset, 
					   aes(x = !!rlang::sym(cellvar), 
					       y = !!rlang::sym(p))) +
					geom_violin(alpha = 0.4) +
					geom_jitter(alpha = .1, shape=".") +
					theme_classic() +
					facet_wrap(as.formula(paste0('~',
								     ref_cluster)),
						nrow=1)
		      				
			exp_threshold <- quantile(ss_umap_eset[[p]], act_quantile)
			ss_umap_eset %>%
				group_by(!!rlang::sym(ref_cluster), !!rlang::sym(cellvar)) %>%
				summarise(n = n(),
					  active_cells = length(which(!!rlang::sym(p) > exp_threshold)), 
					  mean_exp = mean(!!rlang::sym(paste0('mean_',p))),
					  sum_exp = sum(!!rlang::sym(p))) %>%
				mutate(prop_cells = n / sum(n)) %>%
				mutate(prop_act_cells = active_cells / n,
				       mean_exp_in_act_cells = sum_exp / active_cells,
				       cls_cell = paste0(!!rlang::sym(ref_cluster), '_', 
							 !!rlang::sym(cellvar))) -> bubble
			

			g <- ggplot(bubble, aes(x = !!rlang::sym(cellvar), 
						y = 2, 
						group = !!rlang::sym(cellvar))) +
		      		geom_point(aes(size = prop_act_cells, color = mean_exp_in_act_cells)) +
				ylab(p) +
				theme_classic() +
				theme(#axis.title.y=element_blank(),
				      axis.text.y=element_blank(),
				      axis.ticks.y=element_blank(),
				      legend.position="bottom") +
				scale_size(range = c(min(bubble[['prop_act_cells']]*10), 
						     max(bubble[['prop_act_cells']]*10)),
							 name="% active cells") +
		      		scale_color_gradient2(low="white", high="#7f0000") +
				facet_wrap(as.formula(paste0('~',
							     ref_cluster)), 
							     nrow = 1)

			g_barr <- ggplot(bubble, aes(x = !!rlang::sym(cellvar), 
						#y = 2, 
						group = !!rlang::sym(cellvar))) +
		      		geom_bar(aes(y = prop_act_cells,
						fill = mean_exp_in_act_cells),
					       stat = "identity") +
		      		theme_classic() +
		      		theme(legend.position="none") +
				scale_fill_gradient2(low="white", high="#7f0000") +
				facet_wrap(as.formula(paste0('~',
							     ref_cluster)), 
							     nrow = 1)


			head(data.frame(bubble))
			pdf(paste0(outputfolder, '/', jobname, 
				   "_pathway_", p, "_SCT_norm_bubble.pdf"),
				width = 11, height = 11) 
				plot(
				     plot_grid(plotlist = list(g, g_barr, vln_plot, g_wc), ncol = 1)
		    		)
			dev.off()

			



			if(summary(aov_res)[[1]][[5]][1] < 0.05) {


				g <- ggplot(sub_umap_eset, 
				     aes(x = UMAP_1, y = UMAP_2, color = !!rlang::sym(p))) +
		       			geom_point(shape = ".") +
					scale_color_gradient2(low="grey80", 
					      high="blue", 
					      na.value = "grey80") +
		      			ggtitle(p) +
					theme_classic()

				return(g)				
						
			} else {
			
				NULL
	
			}
		}
	})

	title <- ggdraw() + 
		draw_label(
			paste0(jobname, " interesting pathways"),
			fontface = 'bold',
			hjust = 0.5
		) +
	theme(
		plot.margin = margin(0, 0, 0, 0)
	)


	ggs <- Filter(Negate(is.null), ggs)

	print(length(ggs))

	jpeg(paste0(outputfolder, '/', jobname, "_pathways_names_SCT_norm_values.jpeg"), 
	    width = 600*length(ggs), height = 600, pointsize = 74, quality = 100)
		plot(
		    plot_grid(title,
			      plot_grid(plotlist = ggs, ncol = length(ggs)),
			      ncol = 1, rel_heights = c(0.01, 1))
		)
	dev.off()
}
