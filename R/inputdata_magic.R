#
# Imputing missing data values on sparse data sets (Magic R package) ## In dev ##
#

" Imputing missing data values on sparse data sets (Magic R package) ## In dev ##

Usage: inputdata_magic.R --jobname=<value> --rdsfile=<value> --npcs=<value> --genequant=<value> --nclusters=<value> --specie=<value> --npath=<value> --infolder=<folder> --outfolder=<folder> --ncores=<value>

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --jobname=<value>    Descriptive name for your experiment.
  --rdsfile=<file>     The .rds file containing sc object.
  --npcs=<value>       Number of PCs to consider in magic.
  --genequant=<value>  Quantile of genes to discard before magic.
  --nclusters<value>   How many co-expression modules to be defined.
  --specie=<value>     Either hg or mm.
  --npath=<value>      Maximum number of significant (FDR < 0.05) pathways per geneset. 
  --infolder=<path>    Path to the single_cell_data .rds files.Either Seurat or Monocle.  
  --outfolder=<path>   Path to results folder.
  --ncores=<value>     # of processors to use

Authors: 
  Jana Henck and Cesar Prada
e-mail: 
  henck@molgen.mpg.de

"-> doc
library(docopt)
arguments <- docopt(doc, quoted_args=TRUE)
print(arguments)

# --- Dependencies

pkgs <- c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
          'limma', 'S4Vectors', 'SingleCellExperiment',
          'SummarizedExperiment', 'batchelor', 'devtools', 'ggplot2',
          'cowplot', "dplyr","Rmagic","viridis","phateR", "Rfast",
	  "ComplexHeatmap", "enrichR")

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

#install.packages('Rfast')
# --- load the library 
library(Rmagic)
library(ggplot2)
library(viridis)
library(phateR)
library(Rfast)
library(ComplexHeatmap)
library(enrichR)


# Load parameters
jobname <- arguments$jobname
rdsfile <- arguments$rdsfile
npcs <- as.numeric(arguments$npcs)
genefilterquantile <- as.numeric(arguments$genequant)
nclusters <- as.numeric(arguments$nclusters)
specie <- arguments$specie
npath <- as.numeric(arguments$npath)
inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
ncores <- as.numeric(arguments$ncores)

# --- Read the .rds file
scfile <- paste0(outputfolder, '/', rdsfile)

scobject <- readRDS(scfile)

message(paste("The single cell data set object", 
	      scfile, "has been read."))

# check scobject class

objecttype <- class(scobject)
print(objecttype)

if (objecttype == 'cell_data_set') {

	m_dataset <- scobject

	} else if (objecttype == 'Seurat') { 

		m_dataset <- scobject@assays$RNA@counts

	} else {
		stop (paste("Sorry, currently we do not support your dataset type"))
}

# Pre-procesing recomended by the authors
print(dim(m_dataset))
print(m_dataset[1:5, 1:10])

# Filter
print(quantile(rowSums(m_dataset > 0)))
keep_rows <- rowSums(m_dataset > 0) > quantile(rowSums(m_dataset > 0), genefilterquantile)

m_dataset <- m_dataset[names(which(keep_rows)), ]

m_dataset <- library.size.normalize(t(m_dataset), verbose = TRUE)

print(dim(m_dataset))
print(m_dataset[1:5, 1:10])
print(class(m_dataset))

sr_m_dataset <- log10(m_dataset + 1)
sr_m_dataset <- as(sr_m_dataset, "dgCMatrix")

print(sr_m_dataset[1:5, 1:10])
print(class(sr_m_dataset))

magic_out <-  magic(sr_m_dataset, genes = 'all_genes', knn = 5, knn.max = NULL,
       decay = 1, t = 3, npca = npcs, init = NULL, t.max = 20,
       knn.dist.method = "euclidean", verbose = TRUE, n.jobs = 10,
       seed = 123)


# -- get the top N most variable genes
# scran approach to filter genes with informative/biological variation
gvar <- scran::modelGeneVar(t(data.matrix(magic_out$result)))
vargenes <- scran::getTopHVGs(gvar, var.threshold = 0)

print(length(vargenes))

x <- data.matrix(magic_out$result)
x <- x[, vargenes]

corr <- Rfast::cora(x)

hc_out <- hclust(dist(corr), method = "ward.D")

palette = colorRampPalette(c("blue", "white", "red"))(50)

clusterCut <- cutree(hc_out, nclusters)
print(table(clusterCut))

ha = HeatmapAnnotation(cluster = as.factor(paste0('C', clusterCut[hc_out$order])))
ra = rowAnnotation(cluster = as.factor(paste0('C', clusterCut[hc_out$order])))

jpeg(paste0(outputfolder, '/', jobname, '_corrlation_log2_heatmap_hvgenes.jpeg'),
	height = 1980, width = 2280, pointsize = 74, quality = 100)
	Heatmap(corr[hc_out$order, hc_out$order],
		col=palette, 
		show_row_names=FALSE,
		show_column_names=FALSE,
		cluster_rows=FALSE, 
		cluster_columns=FALSE,
		top_annotation=ha,
		right_annotation=ra
	) 
dev.off()

jpeg(paste0(outputfolder, '/', jobname, '_corrlation_log2_clustertree_hvgenes.jpeg'),
	height = 1980, width = 5280, pointsize = 74, quality = 100)
	plot(hc_out)
dev.off()

# Functional enrichment of the co-exp modules
if(specie == "hg") {

	refbiomart <- data.table::fread("<human-biomart-reference-file>")

} else if(specie == "mm"){

	refbiomart <- data.table::fread("<mouse-biomart-reference-file>")

} else {

	stop("Sorry this specie isnt supported, try mm or hg please")

}


# Check if stable gene ID was provided with or without version
g2test <- names(clusterCut[1])

if(grepl("\\.\\d+$", g2test)) {
	names(clusterCut) <- gsub("\\.\\d+$", "", names(clusterCut))
} 

cls <- split(names(clusterCut), as.factor(clusterCut))
cls <- lapply(cls, function(g) unique(dplyr::filter(refbiomart, `Gene stable ID` %in% g)[['Gene name']]))

dbs <- listEnrichrDbs()
print(head(dbs)) 

q_dbs <- c("Human_Gene_Atlas", 
	   "Mouse_Gene_Atlas", 
	   "OMIM_Disease", 
	   "Tissue_Protein_Expression_from_Human_Proteome_Map",
	   "ARCHS4_Tissues",
	   "ARCHS4_TFs_Coexp",
	   "GO_Biological_Process_2018",
	   "GO_Cellular_Component_2018",
	   "GO_Molecular_Function_2018",
	   "KEGG_2019_Human",
	   "KEGG_2019_Mouse",
	   "ProteomicsDB_2020",
	   "Reactome_2016")

enrch <- lapply(cls, function(g) enrichr(g, q_dbs))

gg_enrch <- Reduce(rbind, 
		   lapply(names(enrch), function(cl) {
			   enr <- enrch[[cl]]
			   
			   enr <- Reduce(rbind, lapply(names(enr), function(path) {
							       enres <- enr[[path]]
							       enres[['pathway']] <- path
							       head(enres, npath)
							})
			   		)

			   enr <- dplyr::filter(enr, Adjusted.P.value < 0.05)
			   enr[['module']] <- paste0('M', cl)
			   enr
		   })
	)

ggs <- lapply(unique(gg_enrch$pathway), function(p) {
	
	print(p)

	ggplot(dplyr::filter(gg_enrch, pathway == p), 
	       aes(x = Term, y = -log10(Adjusted.P.value))) +
		geom_bar(stat = "identity", aes(fill = Combined.Score)) +
		theme_classic() +
		coord_flip() +
		ggtitle(paste(p)) +
		facet_wrap(~module, nrow = 1)

	})

pdf(paste0(outputfolder, '/', jobname, '_corexp_modules_enrichment.pdf'),
    width = 15, height = 8)
	for(i in seq(ggs)) {
		plot(ggs[[i]])
	}
dev.off()

# -- Plots and Visualization

g1 <- colnames(corr)[hc_out$order][2]
g2 <- colnames(corr)[hc_out$order][5]

m_df <- data.frame(data.matrix(sr_m_dataset))

g <- ggplot(m_df) +
  geom_point(aes(!!rlang::sym(g1), !!rlang::sym(g2))) +
  scale_color_viridis()

png(paste0(outputfolder, '/', jobname, '_coorrelation_before_magic.png'))
	plot(g)
dev.off()

g <- ggplot(magic_out) +
  geom_point(aes(!!rlang::sym(g1), !!rlang::sym(g2))) +
  scale_color_viridis()

png(paste0(outputfolder, '/', jobname, '_coorrelation_after_magic.png'))
	plot(g)
dev.off()

# ------------------------

# --- Write cell_data_set object with merged samples

mca_file <- paste0(outputfolder, '/', jobname, '_magic_output.rds')

saveRDS(list('magic' = magic_out,
	     'corr' = corr,
	     'hclust' = hc_out), 
	file = mca_file)

message(paste("The R object ", 
	      mca_file, "has been created."))

