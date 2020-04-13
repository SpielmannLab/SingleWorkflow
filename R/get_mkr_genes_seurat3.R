#
# Identify marker genes in a seurat3 object for given cell-phenodata variable
#

"  Identify marker genes in a seurat3 object for given cell-phenodata variable

Usage: get_mkr_genes_seurat3.R --jobname=<value> [--mcafile=<file>] --resolution=<value> --specie=<value> --infolder=<folder> --outfolder=<folder> --ncores=<value> 

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --jobname=<file>     Descriptive name for your experiment.
  --resolution=<value> Clustering resolution. 	
  --specie=<value>     Either hg or mm.
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


inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
jobname <- arguments$jobname
res <- as.numeric(arguments$resolution)
specie <- arguments$specie

# --------------------------------------------
# Some auxiliar functions
# --------------------------------------------

calc_simil <- function(mods, type) { # take a list of string vectors and compare them against them. It returns the -log10(Fisher.pval)
  gene.universe <- as.vector(unlist(mods))
  result <- lapply(names(mods)[-length(names(mods))], function(mod) {
    i <- which(names(mods) == mod) + 1
    similarity <- c()
    while (i <= length(mods)) {
      if(type == "fisher") {
        s <- -log10(fisher_exact_test(mods[[mod]], mods[[i]], gene.universe))
      } else if (type == "intersect") {
        s <- length(intersect(mods[[mod]], mods[[i]]))
      } else if (type == "jaccard") {
        s <- set_similarity(mods[[mod]], mods[[i]], method = "Jaccard")
      }
      similarity <- append(similarity,s)
      i <- i+1
    }
    result = cbind("DEG1" = rep(mod, length(similarity) ), 
                   "DEG2" = names(mods)[(which(names(mods)==mod)+1):length(names(mods))], 
                   "Weight" = as.numeric(similarity))
    return(result)
  })
  return(result)
}

fisher_exact_test <- function(m1, m2, universe) { # Function: Calculate the Fisher Exact Test (@Diogenes)
  isect <- length(intersect(m1, m2))
  m1uniq <- length(m1[!m1 %in% isect])
  m2uniq <- length(m2[!m2 %in% isect])
  tot <- length(universe) - length(unique(c(m1,m2)))
  mm <- c(isect, m1uniq, m2uniq, tot)
  to_fisher <- matrix(mm, nrow = 2)
  if (all(to_fisher > 0)) { 
    k <- fisher.test(to_fisher, alternative = 'greater')$p.value
    return(k)
  } else { # If any cell os the matrix is 0
    return(1)
  }
}

write_simil_matrix <- function(simil.3col) {
  simil.3col[which(simil.3col[,3]=="Inf"),3] <- max(as.numeric(simil.3col[-which(simil.3col[,3]=='Inf'),3]))
  vars <- unique(c(simil.3col[,1], simil.3col[,2]))
  nvars <- length(vars)
  simil.mat <- data.frame(matrix(0, nrow=nvars, ncol=nvars))
  colnames(simil.mat) <- vars
  row.names(simil.mat) <- vars
  for(s in 1:nrow(simil.3col)) {
    simil.mat[simil.3col[s,1], simil.3col[s,2]] <- simil.3col[s,3]
    simil.mat[simil.3col[s,2], simil.3col[s,1]] <- simil.3col[s,3]
  }
  return(simil.mat)
}


# ---------------------------------------------------
# Setting multicore env for some Seurat applications
# ---------------------------------------------------
ncores <- as.numeric(arguments$ncores)
plan("multiprocess", workers = ncores)
plan()

# ----------------------------------------------------
#  Read cell_data_set object with cluster metadata
# ----------------------------------------------------

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

ident_col <- colnames(mca@meta.data)[grep(paste0("_snn_res.", res), colnames(mca@meta.data))]
ident_col <- ident_col[1]
Idents(object = mca) <- ident_col 
print(table(Idents(mca)))

# ----------------------------------------------------
#  Visualizing the selected cluster granularity
# ----------------------------------------------------

g_cluster <- Seurat::DimPlot(mca, 
		     reduction = "umap", 
		     pt.size = 0.01) +
			     ggplot2::theme(legend.position = "bottom") + 
			     ggtitle(paste0(jobname, "_", res))

# --- title
title <- ggdraw() + 
  draw_label(
    paste0(jobname, ' res: ', res),
    fontface = 'bold',,
    hjust = 0.5
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

jpeg(paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "_", res), "_umap.jpeg"), 
     	    width = 980, height = 900, pointsize = 74, quality = 100)
    plot_grid(title, plot_grid(g_cluster, ncol = 1), 
	      ncol = 1, rel_heights = c(0.02, 1))
dev.off()

pdf(paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "_", res), "_umap.pdf"))
    plot_grid(title, plot_grid(g_cluster, ncol = 1), 
	      ncol = 1, rel_heights = c(0.02, 1))
dev.off()

# ----------------------------------------------------
# Identifying the marker genes
# ----------------------------------------------------

mca.markers.roc <- FindAllMarkers(mca,
			      only.pos = TRUE, 
			      min.pct = 0.1,
			      test.use = "roc",
			      logfc.threshold = 0.25)

print(head(mca.markers.roc, 2))

cns <- colnames(mca.markers.roc)
cns[grep("pct|cluster", cns)] <- paste0(cns[grep("pct|cluster", cns)], "_auc")
colnames(mca.markers.roc) <- cns

mca.markers.wilcox <- FindAllMarkers(mca,
			      only.pos = TRUE, 
			      min.pct = 0.1,
			      test.use = "wilcox",
			      logfc.threshold = 0.25)

print(head(mca.markers.wilcox, 2))

cns <- colnames(mca.markers.wilcox)
cns[grep("pct|cluster", cns)] <- paste0(cns[grep("pct|cluster", cns)], "_wilcox")
colnames(mca.markers.wilcox) <- cns

mca.markers <- merge(mca.markers.roc, mca.markers.wilcox, by = "gene", all = TRUE)


# Check if stable gene ID was provided with or without version
g2test <- mca.markers[['gene']][1]

if(grepl("\\.\\d+$", g2test)) {
	mca.markers[['gene_version']] <- mca.markers[['gene']]
	mca.markers[['gene']] <- gsub("\\.\\d+$", "", mca.markers[['gene_version']])
} else {
	mca.markers[['gene_version']] <- mca.markers[['gene']]
}


print(head(mca.markers, 2))


# ------------------------------------------------------------------------------
# Marker expression visualization (umap)
# ------------------------------------------------------------------------------

if(specie == "hg") {

	refbiomart <- data.table::fread("<biomart>.tsv")
	mca.markers <- merge(mca.markers, refbiomart, by.x = "gene", by.y = "Gene stable ID", all.x = TRUE) 

} else if(specie == "mm"){

	refbiomart <- data.table::fread("<biomart>.tsv")
	mca.markers <- merge(mca.markers, refbiomart, by.x = "gene", by.y = "Gene stable ID", all.x = TRUE) 

} else {

	stop("Sorry this specie isnt supported, try mm or hg please")

}

print(head(mca.markers))

write.table(mca.markers, paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "", res), "_marker_genes.tsv"),
	    sep = "\t", quote = FALSE, row.names = FALSE)


auc <- dplyr::filter(mca.markers, !is.na(cluster_auc))
auc <- split(auc, auc[['cluster_auc']]) 
auc <- lapply(auc, function(w) head(dplyr::arrange(w, -myAUC), 200))
topmarkers <- Reduce(rbind, lapply(auc, function(top) {
					   top <- dplyr::mutate(top, d = pct.1_auc - pct.2_auc)
					   head(dplyr::arrange(top, -d), 100)
	})
)


print(table(topmarkers$cluster_auc))

dim(topmarkers)

write.table(topmarkers, paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "", res),"_top100_marker_genes.tsv"),
	    sep = "\t", quote = FALSE, row.names = FALSE)

mkrs <- split(topmarkers, topmarkers$cluster_auc) 
mkrs <- Filter(Negate(is.null), mkrs)

if(any(grepl("gene_version", colnames(mkrs[[1]])))) {
	   mkrs <- lapply(mkrs, function(m) head(m$gene_version, 50))
} else {
	   mkrs <- lapply(mkrs, function(m) head(m$gene, 50))
	
}


parallel::mclapply(names(mkrs), function(cls) {

	g <- mkrs[[cls]]
	gnames <- filter(topmarkers, gene%in%g)[['Gene name']]

	ggs <- FeaturePlot(mca, 
			   features = unique(g), 
			   combine=FALSE)

	ggs <- lapply(seq(ggs), function(i) {
			      igene <- names(ggs[[i]]$data)[4]
			      ggs[[i]] + 
				      ggtitle(paste(unique(dplyr::filter(topmarkers,
							    gene_version==igene)[['Gene name']]),
				      collapse = '///')) 
			   })

	title <- ggdraw() + 
		draw_label(
			paste0(jobname, " C", cls),
			fontface = 'bold',
			hjust = 0.5
		) +
	theme(
		plot.margin = margin(0, 0, 0, 0)
  	)

	jpeg(paste0(outputfolder, '/', jobname, "_C", cls, "_res_", gsub("\\.", "", res), "_auc_marker_genes.jpeg"), 
	    width = 1980, height = 1980, pointsize = 74, quality = 100)

	plot(
	    plot_grid(title,
		      plot_grid(plotlist = ggs, ncol = 5),
		      ncol = 1, rel_heights = c(0.01, 1))
	)

	dev.off()
}, mc.cores = ncores)

# ------------------------------------------------------------------------------
# Marker expression visualization (heatmap) [AUC]
# ------------------------------------------------------------------------------

if(length(colnames(mca)) > 10000) {
	ncell <- 10000
} else {
	ncell <- length(colnames(mca))
}

cells <- sample(colnames(mca), ncell, replace=FALSE)

print(head(cells))

pdf(paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "", res), "_auc_marker_genes_heatmap.pdf"), 
     width = 18, height = 15)
plot(
DoHeatmap(
  mca,
  features = unlist(mkrs),# arrange(topmarkers, `cluster`)[['gene']],
  cells = cells,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = "SCT",
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
) + scale_fill_gradientn(colors = c("#92c5de", "#f7f7f7", "#b2182b"))
)
dev.off()

jpeg(paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "", res),"_auc_marker_genes_heatmap.jpeg"), 
 	    width = 1980, height = 980, pointsize = 74, quality = 100)
plot(
DoHeatmap(
  mca,
  features = unlist(mkrs),
  cells = cells,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = "SCT",
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
) + scale_fill_gradientn(colors = c("#92c5de", "#f7f7f7", "#b2182b"))
)
dev.off()


#
# Wilcoxon
#

wilcox <- dplyr::filter(mca.markers, !is.na(cluster_wilcox))
wilcox <- split(wilcox, wilcox[['cluster_wilcox']]) 
wilcox <- lapply(wilcox, function(w) head(dplyr::arrange(w, p_val_adj), 200))
topmarkers <- Reduce(rbind, lapply(wilcox, function(top) {
					   top <- dplyr::mutate(top, d = pct.1_wilcox - pct.2_wilcox)
					   head(dplyr::arrange(top, -d), 100)
	})
)

print(table(topmarkers$cluster_wilcox))

write.table(topmarkers, paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "", res),"_wilcox_top100_marker_genes.tsv"),
	    sep = "\t", quote = FALSE, row.names = FALSE)

mkrs <- split(topmarkers, topmarkers$cluster_auc) 
mkrs <- Filter(Negate(is.null), mkrs)

if(any(grepl("gene_version", colnames(mkrs[[1]])))) {
	   mkrs <- lapply(mkrs, function(m) head(m$gene_version, 50))
} else {
	   mkrs <- lapply(mkrs, function(m) head(m$gene, 50))
	
}


parallel::mclapply(names(mkrs), function(cls) {

	g <- mkrs[[cls]]
	gnames <- filter(topmarkers, gene%in%g)[['Gene name']]

	ggs <- FeaturePlot(mca, 
			   features = unique(g), 
			   combine=FALSE)

	ggs <- lapply(seq(ggs), function(i) {
			      igene <- names(ggs[[i]]$data)[4]
			      ggs[[i]] + 
				      ggtitle(paste(unique(dplyr::filter(topmarkers,
							    gene_version==igene)[['Gene name']]),
				      collapse = '///')) 
			   })

	title <- ggdraw() + 
		draw_label(
			paste0(jobname, " C", cls),
			fontface = 'bold',
			hjust = 0.5
		) +
	theme(
		plot.margin = margin(0, 0, 0, 0)
  	)

	jpeg(paste0(outputfolder, '/', jobname, "_C", cls, "_res_", gsub("\\.", "", res), "_wilcox_marker_genes.jpeg"), 
	    width = 1980, height = 1980, pointsize = 74, quality = 100)

	plot(
	    plot_grid(title,
		      plot_grid(plotlist = ggs, ncol = 5),
		      ncol = 1, rel_heights = c(0.01, 1))
	)

	dev.off()
}, mc.cores = ncores)


# ------------------------------------------------------------------------------
# Marker expression visualization (heatmap)
# ------------------------------------------------------------------------------

if(length(colnames(mca)) > 10000) {
	ncell <- 10000
} else {
	ncell <- length(colnames(mca))
}

cells <- sample(colnames(mca), ncell, replace=FALSE)

print(head(cells))

pdf(paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "", res), "_wilcox_marker_genes_heatmap.pdf"), 
     width = 18, height = 15)
plot(
DoHeatmap(
  mca,
  features = unlist(mkrs),# arrange(topmarkers, `cluster`)[['gene']],
  cells = cells,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = "SCT",
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
) + scale_fill_gradientn(colors = c("#92c5de", "#f7f7f7", "#b2182b"))
)
dev.off()

jpeg(paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "", res),"_wilcox_marker_genes_heatmap.jpeg"), 
 	    width = 1980, height = 980, pointsize = 74, quality = 100)
plot(
DoHeatmap(
  mca,
  features = unlist(mkrs),
  cells = cells,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = "SCT",
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
) + scale_fill_gradientn(colors = c("#92c5de", "#f7f7f7", "#b2182b"))
)
dev.off()

# -----------------------------------------------------
# Marker annotation
# -----------------------------------------------------
mca.markers %>%
	filter(!is.na(cluster_auc)) %>%
	group_by(cluster_auc) %>% 
	top_n(n = 200, wt = myAUC) -> topmarkers

topmarkers <- split(topmarkers, as.factor(topmarkers$cluster_auc))
topmarkers <- Reduce(rbind, lapply(topmarkers, function(top) head(arrange(top, -power), 50)))

mkrs <- split(topmarkers, topmarkers[['cluster_auc']])

genes <- lapply(mkrs, function(x) unique(append(toupper(x$`HGNC symbol`), toupper(x$`Gene name`))))

# Refs
refiles <- c("<path-to>PanglaoDB_markers_19_Nov_2019.tsv", 
	     "<pathto>cellmarkerDB_20191120.tsv")

names(refiles) <- c("PanglaoDB", "CellMarkerDB")

refs <- lapply(refiles, function(f) fread(f))

print("Marker gene references:")
print(lapply(refs, dim))
print(lapply(refs, function(g) head(g, 2)))

panglaodb <- refs[['PanglaoDB']] %>%
	dplyr::mutate(cell_type = paste0(`cell type`, '_', organ, `germ layer`)) %>%
	select(c(`official gene symbol`, cell_type))

panglaodb <- split(panglaodb, panglaodb[['cell_type']]) 
panglaodb <- lapply(panglaodb, function(x) x$`official gene symbol`)
print(panglaodb[1:4])

cmdb <- refs[['CellMarkerDB']] %>%
	dplyr::mutate(cell_type = paste0(cellName, '_', tissueType, '_', cellType)) %>%
	select(c(geneSymbol, cell_type))

cmdb <- split(cmdb, cmdb[['cell_type']]) 
cmdb <- lapply(cmdb, function(x) toupper(unlist(strsplit(x$geneSymbol, ", "))))
print(cmdb[1:4])

names(genes) <- paste0("C", names(genes))

igenes <- append(genes, cmdb)
smilmat <- calc_simil(igenes, "fisher")
smilmat <- Reduce(rbind, smilmat)
smilmat <- write_simil_matrix(smilmat)
mat <- smilmat[names(cmdb), names(genes)]

max_enrch <- apply(mat, 1, max)

idx <- tail(sort(max_enrch), length(genes)*4)

mat <- mat[-which(max_enrch < idx[1]), ]

dim(mat)
head(mat)

colfunc <- colorRampPalette(c("#deebf7", "#9ecae1", "#3182bd"))
pdf(paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "", res), "_hm_similarity_gene_cell_markers_fisher_CellMarkerDB.pdf"),
    width = 30, height = 15)
  heatmap(data.matrix(mat), scale = "none", Rowv=NULL, Colv=NA, col = colfunc(1000))
dev.off()

jpeg(paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "", res),
	    "_hm_similarity_gene_cell_markers_fisher_CellMarkerDB.jpeg"),
    width = 3500, height = 1900, quality = 100, pointsize=36)
       heatmap(data.matrix(mat), scale = "none", Rowv=NULL, Colv=NA, col = colfunc(1000))
dev.off()

igenes <- append(genes, panglaodb)
smilmat <- calc_simil(igenes, "fisher")
smilmat <- Reduce(rbind, smilmat)
smilmat <- write_simil_matrix(smilmat)
mat <- smilmat[names(panglaodb), names(genes)]

max_enrch <- apply(mat, 1, max)

idx <- tail(sort(max_enrch), length(genes)*4)

mat <- mat[-which(max_enrch < idx[1]), ]

dim(mat)
head(mat)

pdf(paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "", res), "_hm_similarity_gene_cell_markers_fisher_panglaoDB.pdf"),
    width = 30, height = 15)
  heatmap(data.matrix(mat), scale = "none", Rowv=NULL, Colv=NA, col = colfunc(1000))
dev.off()

jpeg(paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "", res),
	    "_hm_similarity_gene_cell_markers_fisher_panglaoDB.jpeg"),
    width = 3500, height = 1900, quality = 100, pointsize=36)
       heatmap(data.matrix(mat), scale = "none", Rowv=NULL, Colv=NA, col = colfunc(1000))
dev.off()
