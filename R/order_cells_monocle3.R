#
# Order cells (pseudotime) for a monocle3 object following the monocle3 approach ## In dev ##
#

" Order cells (pseudotime) for a monocle3 object following the monocle3 approach ## In dev ##

Usage: ordercells_monocle3.R --jobname=<value> [--cdsfile=<file>] --npcs=<value> --infolder=<folder> --outfolder=<folder> --ncores=<value> 

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --jobname=<file>     Descriptive name for your experiment.
  --cdsfile=<file>     Optional. If !is.null, do not use jobname to guess mcafilename.
  --npcs=<value>       Number of principal components to use.
  --infolder=<file>    Path to the single_cell_data .rds files.  
  --outfolder=<file>   Path to results folder.
  --ncores=<value>     Number of threads to use.

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
          'cowplot', "dplyr", "phateR", "Rmagic", "monocle3", "slingshot")

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
pymagic <- reticulate::import("magic")
suppressMessages(library(Seurat))
#suppressMessages(library(batchelor))
#suppressMessages(library(dplyr))
#suppressMessages(library(ggwordcloud))
# ------------------------------

# --------------
# Main function
# --------------

# ---------- Functions ---------

getcluster <- function(cds, res) {
	cds <- cluster_cells(cds, 
			     reduction_method = "Aligned",#c("UMAP", "tSNE", "PCA", "LSI", "Aligned"), 
			     k = 20, 
			     cluster_method = "louvain",
			     num_iter = 2, 
			     #partition_qval = 0.05, 
			     weight = FALSE,
			     resolution = res, 
			     random_seed = NULL,
			     verbose = TRUE)

	cds@colData$tmpcluster <- cds@clusters$Aligned$clusters
	#cds@colData$tmpcluster <- cds@clusters$PCA$clusters

	g_cluster <- plot_cells(cds,
				color_cells_by="tmpcluster")
	
	g_partition <- plot_cells(cds, 
				  color_cells_by="partition", 
				  group_cells_by="partition")

	clstr <- data.frame('barcode' = names(cds@clusters$Aligned$clusters),
			    'cluster' = cds@clusters$Aligned$clusters)

	colnames(clstr)[2] <- paste0(colnames(clstr)[2], '_', 
				     gsub("\\.", "_", res))

	return(list('g_cl' = g_cluster,
		    'g_pr' = g_partition,
		    'cluster' = clstr))
}


# ------------------------------

inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
jobname <- arguments$jobname
npcs <- as.numeric(arguments$npcs)
testpcs <- as.logical(arguments$testpcs)


# --- Read cell_data_set object with merged samples
if(!is.null(arguments$cdsfile)) {
	cds_file <- paste0(inputfolder, '/', arguments$cdsfile)
	jobname <- gsub("\\.rds", "", arguments$cdsfile)
} else {
	cds_file <- paste0(inputfolder, '/', jobname, '_monocle3_merged_clustered_resrange.rds')
}

cds <- readRDS(cds_file)


message(paste("The cell_data_set object (Monocle3)", 
		cds_file, "has been read."))

print(dim(cds))


# For cds objects after filtering with no pre-processing

cds <- preprocess_cds(cds, 
		      method = "PCA", 
		      num_dim = npcs,
		      norm_method = "log",
		      use_genes = NULL,
		      scaling = TRUE,  
		      verbose = TRUE, 
		      cores = ncores)

cds <- reduce_dimension(cds, 
		 max_components = 2, 
		 reduction_method = "UMAP",
		 preprocess_method = "PCA", 
		 umap.metric = "cosine", 
		 umap.min_dist = 0.1, 
		 umap.n_neighbors = 15L, 
		 umap.fast_sgd = FALSE, 
		 umap.nn_method = "annoy", 
		 cores = ncores, 
		 verbose = TRUE)

cds <- cluster_cells(cds, 
		     reduction_method = "PCA",#c("UMAP", "tSNE", "PCA", "LSI", "Aligned"), 
		     k = 20, 
		     cluster_method = "louvain",
		     num_iter = 10, 
		     partition_qval = 0.05, 
		     weight = FALSE,
		     resolution = 0.0005, 
		     random_seed = NULL,
		     verbose = TRUE)

cds <- cluster_cells(cds,  
		     reduction_method = "UMAP",#c("UMAP", "tSNE", "PCA", "LSI", "Aligned"), 
		     k = 20, 
		     cluster_method = "louvain",
		     num_iter = 10, 
		     partition_qval = 0.05, 
		     weight = FALSE,
		     #resolution = 0.00000001, 
		     random_seed = NULL,
		     verbose = TRUE)

cds@colData$tmpcluster <- cds@clusters$PCA$clusters

g_cluster <- plot_cells(cds,
			color_cells_by="tmpcluster")

pdf(paste0(outputfolder, "/", jobname, "cluster.pdf")
plot(g_cluster)
dev.off()

g_partition <- plot_cells(cds, 
			  color_cells_by="partition",
			  group_cells_by="partition")

pdf(paste0(outputfolder, "/", jobname, "_partition_cluster.pdf")
plot(g_partition)
dev.off()

cds <- learn_graph(cds,
		   use_partition = TRUE,
		   close_loop = TRUE,
		   learn_graph_control = NULL,
		   verbose = TRUE)

g_partition <- plot_cells(cds, 
			  color_cells_by="partition", 
			  group_cells_by="partition",
			  label_leaves=FALSE,
			  label_branch_points=FALSE)

pdf(paste0(outputfolder, "/", jobname, "_trajectory.pdf")
plot(g_partition)
dev.off()

get_earliest_principal_node <- function(cds, cell_ids){

#  cell_ids <- rownames(colData(cds))
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex

  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])

  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))]
  
  root_pr_nodes
}

# B6
cell_ids <- names(partitions(cds)[partitions(cds) == 3])
# D129
cell_ids <- names(partitions(cds)[partitions(cds) == 5])

cds <- order_cells(cds,
		   root_pr_nodes = get_earliest_principal_node(cds, cell_ids))

g_ptime <- plot_cells(cds, 
		      color_cells_by="pseudotime", 
		      label_leaves=FALSE,
		      label_branch_points=FALSE)


pdf(paste0(outputfolder, "/", jobname, "_pseudotime_cluster.pdf")
	plot(g_ptime)
dev.off()
