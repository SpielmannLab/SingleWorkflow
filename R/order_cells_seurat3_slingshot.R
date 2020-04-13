#
# Order cells (pseudotime) for a seurat3 object following the slingshot approach ## IN DEV ##
#

" Order cells (pseudotime) for a seurat3 object following the slingshot approach ## IN DEV ##

Usage: order_cells_seurat3_slingshot.R --jobname=<value> [--cdsfile=<file>] --specie=<value> --infolder=<folder> --outfolder=<folder> --ncores=<value> [--res=<value>]

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --jobname=<file>     Descriptive name for your experiment.
  --cdsfile=<file>     Optional. If !is.null, do not use jobname to guess mcafilename.
  --specie=<value>     Specie e.g. either mm or hg.
  --infolder=<file>    Path to the single_cell_data .rds files.  
  --outfolder=<file>   Path to results folder.
  --ncores=<value>     Number of threads to use.
  --res=<value>        Cluster resolution.

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
          'cowplot', "dplyr", "phateR", "Rmagic", "monocle3", 
	  "slingshot", "viridis", "gam")

if (!require("BiocManager", character.only = TRUE)) {
	install.packages("BiocManager")
	BiocManager::install()
} else {
	ipkgs <- sapply(pkgs, function(...) require(..., character.only = TRUE))
	if (any(!ipkgs)) {
		BiocManager::install(pkgs[!ipkgs])
	} else {
		message("\n\nCool! your machine has everything is needed.\n\n")
	}
}

# Install last version of both (Beta version keep changing)

#devtools::install_github('cole-trapnell-lab/leidenbase')
#devtools::install_github('cole-trapnell-lab/monocle3', force=TRUE)
# Phate installation
#reticulate::py_install("phate", pip=TRUE)
#reticulate::py_install("magic-impute", pip=TRUE)
#devtools::install_github("KrishnaswamyLab/phateR")
#devtools::install_github("KrishnaswamyLab/MAGIC/Rmagic")


# -------------------------------
suppressMessages(library(monocle3))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(phateR))
suppressMessages(library(Rmagic))
suppressMessages(library(monocle3))
#pymagic <- reticulate::import("magic")
suppressMessages(library("viridis"))
suppressMessages(library(Seurat))
#suppressMessages(library(batchelor))
#suppressMessages(library(dplyr))
#suppressMessages(library(ggwordcloud))
# ------------------------------

inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
jobname <- arguments$jobname
npcs <- as.numeric(arguments$npcs)
specie <- arguments$specie
ncores <- as.numeric(arguments$ncores)
res <- as.numeric(arguments$res)

# --- Read cell_data_set object with merged samples

if(!is.null(arguments$cdsfile)) {
	cds_file <- paste0(inputfolder, '/', arguments$cdsfile)
	jobname <- gsub("_seurat3_subeset_seurat3_merged_clustered_resrange.rds", "", arguments$cdsfile)
} else {
	stop("Please provide an specific Seurat3 file --cdsfile")
}

cds <- readRDS(cds_file)

message(paste("The cell_data_set object (Monocle3)", 
		cds_file, "has been read."))

print(dim(cds))

# Getting PCA
pca <- cds@reductions$pca@cell.embeddings
print(dim(pca))
print(head(pca,2))


ident_col <- colnames(cds@meta.data)[grep(paste0("_snn_res.", res), colnames(cds@meta.data))]
print(table(cds@meta.data[[ident_col]]))

scale_mat <- cds@assays$integrated@scale.data
npcs <- 15
mclapply(unique(cds@meta.data[[ident_col]]), function(r) {

	    sub_scale <- scale_mat[, rownames(cds@meta.data)[which(cds@meta.data[[ident_col]] == r)]]
	    
	    pcas <- prcomp(t(sub_scale))

	    pca <- pcas$x[, seq(npcs)]

	    lin <- getLineages(data=pca)

	    lin
	
	jobname <- paste0(jobname, '_res_', gsub("\\.", "_", res), '_C', r)

	ptime <- slingshot(data=pca, 
#		   clusterLabels=as.character(cds@meta.data[[ident_col]]),
		   reducedDim = NULL, 
		   start.clus = NULL, 
		   end.clus = NULL,
		   dist.fun = NULL, 
		   omega = NULL, 
		   lineages = lin, 
		   shrink = TRUE,
		   extend = "y", 
		   reweight = TRUE, 
		   reassign = TRUE, 
		   thresh = 0.001,
		   maxit = 15, stretch = 2, 
		   approx_points = FALSE,
		   smoother = "smooth.spline", 
		   shrink.method = "cosine",
		   allow.breaks = TRUE)
	ptime
 	class(ptime) 
	summary(ptime)

	ps_time <- data.frame(slingshot::slingPseudotime(ptime))
	ps_time[['barcode']] <- rownames(ps_time)
	summary(ps_time) 

	df_plot <- data.frame(reducedDim(ptime))
	df_plot[['barcode']] <- rownames(df_plot)
	df_plot <- merge(df_plot, ps_time, by = "barcode")

	pseudo_meta <- merge(df_plot, cds@meta.data, by = "barcode")
	umap <- data.frame(cds@reductions$umap@cell.embeddings)
	umap[['barcode']] <- rownames(umap)

	pseudo_meta <- merge(pseudo_meta, umap, by = 'barcode')
  
	head(pseudo_meta)

	write.table(pseudo_meta, paste0(outputfolder, '/', jobname, "_meta_trajectory_pseudotime_cluster.tsv"),
					sep = "\t", quote = FALSE, row.names = FALSE)

	g <- ggplot(pseudo_meta, aes(x = UMAP_1, y = UMAP_2, color = !!rlang::sym(ident_col))) +
		geom_point(aes(shape = "."), alpha = .3) +
		theme_classic() 
	

	cs <- colnames(pseudo_meta)[grep("curve", colnames(pseudo_meta))]
	cs

	ggs <- lapply(cs, function(traj) {
		ggplot(pseudo_meta, aes(x = UMAP_1, y = UMAP_2, color = !!rlang::sym(traj))) +
			geom_point(aes(shape = "."), alpha = .3) +
			theme_classic() +
			scale_color_gradient2(low="grey", high="red")
		})
	
	pdf(paste0(outputfolder, '/', jobname, "_UMAP_trajectory_pseudotime_cluster.pdf"), height = 4, width = 4*length(ggs))
		plot(
			plot_grid(plotlist = ggs, ncol = 5)
		)
	dev.off()
	
	
	g <- ggplot(pseudo_meta, aes(x = PC1, y = PC2, color = !!rlang::sym(ident_col))) +
		geom_point(aes(shape = "."), alpha = .3) +
		theme_classic() 

	pdf(paste0(outputfolder, '/', jobname, "_PCA_cluster_trajectory.pdf"))
		plot(g)
	dev.off()

	g <- ggplot(df_plot, aes(x = PC1, y = PC2, color = curve1)) +
		geom_point(aes(shape = "."), alpha = .3) +
		theme_classic() 

	pdf(paste0(outputfolder, '/', jobname, "_PCA_pseudotime_trajectory.pdf"))
		plot(g)
	dev.off()

	colors <- colorRampPalette(magma(10))(100)
	plotcol <- colors[cut(ptime@curves$curve1$lambda, breaks=100)]

	pdf(paste0(outputfolder, '/', jobname, "_pca_pseudotime_trajectory_genes.pdf"))
		plot(reducedDim(ptime), col = plotcol, pch = 16, asp = 1, cex = 0.5)
		lines(ptime, lwd = 2, col = 'black')
	dev.off()

	print(head(df_plot))

	sct_eset <- cds@assays$integrated@scale.data

	var_genes <- names(head(sort(apply(sct_eset, 1, var), decreasing = TRUE), 1000))

	sct_eset <- data.frame(t(sct_eset[var_genes, ]))
	colnames(sct_eset) <- gsub("\\.\\d+", "", colnames(sct_eset))

	sct_eset[["barcode"]] <- rownames(sct_eset)

	sct_eset <- merge(sct_eset, df_plot, by = "barcode")

	sct_eset <- merge(sct_eset, select(pseudo_meta, 
				   matches(paste("barcode", 
						 ident_col, 
						 sep = "|"))
				   ),
			  by = "barcode")

	head(sct_eset, 2)

	gam_pval <- apply(dplyr::select(sct_eset, matches("ENSG")), 2, function(gene) {
			d <- data.frame(z = gene,
					t = sct_eset[['curve1']],
					cl = sct_eset[[ident_col]])

			suppressWarnings({
	
				  tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
	
			})
			
			g <- ggplot(d, aes(x = t, y = z)) +
				geom_point(shape = ".",aes(color = cl)) +
				geom_smooth() +
				theme_classic()
	
			p <- summary(tmp)[3][[1]][2,3]
			list('p' = p,
			     'g' = g)
		})
	
	topgenes <- names(sort(sapply(gam_pval, function(x) x$p), decreasing = FALSE))[1:75]
	
	heatdata <- cds@assays$integrated@scale.data

	rownames(heatdata) <- gsub("\\.\\d+", "", rownames(heatdata))
	
	heatdata <- heatdata[topgenes, sct_eset[['barcode']][order(sct_eset[['curve1']], na.last = NA)]]
	
	# for time, only look at the 100 most variable genes
	if(specie == "hg") {

		refbiomart <- data.table::fread("<human-biomart-reference-file>")

	} else if(specie == "mm"){

		refbiomart <- data.table::fread("<mouse-biomart-reference-file>")

	} else {

		stop("Sorry this specie isnt supported, try mm or hg please")

	}
	
	rownames(heatdata) <- sapply(rownames(heatdata), function(g) paste(filter(refbiomart, `Gene stable ID` == g)[['Gene name']], collapse = '///'))
	
	write.table(heatdata, paste0(outputfolder, '/', jobname, "_meta_trajectory_pseudotime_dynamic_mkrs_heatdata.tsv"),
				     sep = "\t", quote = FALSE, row.names = TRUE)


	# cluster annotation
	pdf(paste0(outputfolder, '/', jobname, "_dynamic_genes_heatmap.pdf"), heigh = 16, width = 16)
		heatmap(heatdata, Colv = NA, scale = "none")
	dev.off()
	# cluster annotation
	jpeg(paste0(outputfolder, '/', jobname, "_dynamic_genes_heatmap.jpeg"),
	    width = 4980, height = 4980, pointsize = 74, quality = 100)
		heatmap(heatdata, Colv = NA, scale = "none")
	dev.off()
	
	
	jpeg(paste0(outputfolder, '/', jobname, "_dynamic_marker_genes.jpeg"),
	     width = 2980, height = 3980, pointsize = 74, quality = 100)
		plot(
		     plot_grid(plotlist = lapply(names(gam_pval[topgenes]), 
						 function(x) {
							 gam_pval[[x]]$g + 
								 ggtitle(paste(dplyr::filter(refbiomart, 
											     `Gene stable ID` == x)[['Gene name']], 
									       collapse = "///"))
						 }), ncol = 10)
		     )
	dev.off()
	
	
	pdf(paste0(outputfolder, '/', jobname, "_dynamic_marker_genes.pdf"),
	     width = 29, height = 39)
		plot(
		     plot_grid(plotlist = lapply(names(gam_pval[topgenes]), 
						 function(x) {
							 gam_pval[[x]]$g + 
								 ggtitle(paste(dplyr::filter(refbiomart, 
											     `Gene stable ID` == x)[['Gene name']], 
									       collapse = "///"))
						 }), ncol = 10)
		     )
	dev.off()
}, mc.cores = ncores)
