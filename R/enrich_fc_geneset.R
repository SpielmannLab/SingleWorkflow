#
# Gene set enrichment analysis based on a differential expression result and gene sets ## In dev ##
#

" Gene set enrichment analysis based on a differential expression result and gene sets ## In dev ##

Usage: enrich_fc_geneset.R --jobname=<value> --defile=<value> --gsfile=<value> --specie=<value> --infolder=<folder> --outfolder=<folder> --ncores=<value>

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --jobname=<value>    Descriptive name for your experiment.
  --defile=<file>     The .rds file containing sc object.
  --gsfile=<value>   How many co-expression modules to be defined.
  --specie=<value>     Either hg or mm.
  --infolder=<path>    Path to the single_cell_data .rds files.Either Seurat or Monocle.  
  --outfolder=<path>   Path to results folder.
  --ncores=<value>     # of processors to use

Authors: 
  Cesar Prada
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
	  "ComplexHeatmap", "enrichR", "fgsea", "circlize")

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


#method <- "DESeq2"
#defile <- paste0("_deg_", method, "_pbulk_results.tsv")

# Load parameters
gmtfolder <- arguments$gmtfolder
jobname <- arguments$jobname
defile <- arguments$rdsfile
gsfile <- arguments$gsfile
specie <- arguments$specie
inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
ncores <- as.numeric(arguments$ncores)
method <- arguments$method

# Differential expression files
defiles <- list.files(path = inputfolder, pattern = defile)
des <- lapply(defiles, function(f) data.table::fread(paste0(inputfolder, '/', f)))
des <- setNames(des, gsub(defile, "", defiles))
print(sapply(des, dim))

# gmt file with genesets
if(!is.null(gmtfolder)) {
	gmt_files <- list.files(path = gmtfolder) 
	gmts <- lapply(gmt_files, function(...) CEMiTool::read_gmt(paste0(gmtfolder, '/', ...)))
	gmts <- lapply(gmts, function(g) {
			       g[['gene']] <- gsub(",-?.+", "", g[['gene']])
			       g
	  })
	gmts <- setNames(gmts, gsub(".txt", "", gmt_files))
} else {
	gmts <- list(CEMiTool::read_gmt(gsfile))
	gmt_name <- unlist(strsplit(gsfile, '/'))
	gmt_name <- gmt_name[length(gmt_name)]
	names(gmts) <- gmt_name 
}


# -------------------------------------
# fGSEA per gmt and for all DE results
# -------------------------------------
fc_var <- "log2FoldChange"
p_thr <- 0.01
library(circlize)

lapply(names(gmts)[c(6:9)], function(gmt_name) {
	       fgsea_gmt(des, gmts[[gmt_name]], gmt_name, fc_var, p_thr)
	  })

fgsea_gmt <- function(des, gmt, gmt_name, fc_var, p_thr) {
	#
	gmt <- split(gmt, gmt[['term']])
	gmt <- lapply(gmt, function(g) as.character(g[['gene']]))
	print(str(gmt))
	#
	fgseares <- lapply(names(des), function(dn) {
	       print(dn)
	       d <- des[[dn]]
	       d %>%
		       dplyr::filter(`Gene name` != "") %>%
		       dplyr::filter(!is.na(`Gene name`)) %>%
		       dplyr::arrange(pvalue) %>%
		       dplyr::filter(!duplicated(`Gene name`)) %>%
		       dplyr::arrange(!!rlang::sym(fc_var)) -> fc_rank

		fc_rank <- setNames(fc_rank[[fc_var]], fc_rank[['Gene name']])
	        # 
		print(head(fc_rank,10))
	        print(tail(fc_rank,10))
		#

		fgseares <- fgsea(pathways = gmt,
				  stats = fc_rank,
				  minSize=5,
				  maxSize=500,
				  nperm=10000)

		#
		fgseares[['de_result']] <- dn
		
		# Package visualization	for top pathways
		topPathwaysUp <- fgseares[ES > 0][head(order(pval), n=10), pathway]
		topPathwaysDown <- fgseares[ES < 0][head(order(pval), n=10), pathway]
		topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

		pdf(paste0("~/fgsea_test_", gmt_name, ".pdf"))
			plotGseaTable(gmt[topPathways], 
				      fc_rank, 
				      fgseares, 
				      gseaParam = 0.5)
		dev.off()
	        # ------
	        return(fgseares)
	})
	#
	tab_res <- Reduce(rbind, fgseares)
	#
	# --Visualizing significancy
	
	tab_vis <- lapply(fgseares, function(f) {
			  p <- unique(f[['de_result']])
			  sf <- select(f, c(pathway, padj))
			  colnames(sf)[2] <- p
			  sf
	  })

	tab_vis <- Reduce(function(x, y) merge(x, y, 
					       by = "pathway", 
					       all = TRUE), 
			  tab_vis)
	#
	print(head(tab_vis, 2))
	print(dim(tab_vis))
	#
	rownames(tab_vis) <- tab_vis[['pathway']]
	tab_vis <- select(tab_vis, -pathway)
	# Heatmap
	x <- data.matrix(-log10(tab_vis))
	rownames(x) <- rownames(tab_vis)
	print(dim(x))
	print(head(x,2))
	
	col_fun = colorRamp2(c(0, 1, 2), 
			     c("#ffffd9", "#7fcdbb", "#225ea8"))
	pdf(paste0("~/sign_NES_", gmt_name, "_heatmap.pdf"), width = 6, height = 15)
	ht <- Heatmap(x, 
	      cluster_rows = FALSE, 
	      cluster_columns = FALSE,
	      col = col_fun,
	      heatmap_legend_param = list(title_position = "topcenter",
				    title = "-log10(FDR)", 
				    color_bar = "continuous", 
				    legend_direction = "horizontal",
				    legend_width = unit(5, "cm")))
		draw(ht, heatmap_legend_side = "top")
	dev.off()

	# -- Visualizing DE expression (NES)
	tab_vis <- lapply(fgseares, function(f) {
			  p <- unique(f[['de_result']])
			  sf <- filter(f, padj < p_thr) %>% 
				  select(c(pathway, NES))
			  colnames(sf)[2] <- p
			  sf
	})

	tab_vis <- Reduce(function(x, y) merge(x, y, by = "pathway", all = TRUE), tab_vis)
	head(tab_vis)

	rownames(tab_vis) <- tab_vis[['pathway']]
	tab_vis <- select(tab_vis, -pathway)

	x <- data.matrix(tab_vis)
	rownames(x) <- rownames(tab_vis)
	dim(x)
	head(x)

	col_fun = colorRamp2(c(-2, -1, 1, 2), c("#3288bd", "white", "white", "#d53e4f"))

	pdf(paste0("~/NES_", gmt_name, "_heatmap.pdf"), width = 6, height = 15)
	ht <- Heatmap(x, 
	      cluster_rows = FALSE, 
	      cluster_columns = FALSE,
	      col = col_fun,
	      heatmap_legend_param = list(title_position = "topcenter",
				    title = "Normalized Enrichment Score", 
				    color_bar = "continuous", 
				    legend_direction = "horizontal",
				    legend_width = unit(5, "cm")))
		draw(ht, heatmap_legend_side = "top")
	dev.off()
}





# -----------------------------------------------
# Fisher exact test
# Using Enrichr from Mayaan lab
# Functional enrichment
# -----------------------------------------------

# Gene annotation
if(specie == "hg") {

	refbiomart <- data.table::fread("<human-biomart-reference-file>")

} else if(specie == "mm"){

	refbiomart <- data.table::fread("<mouse-biomart-reference-file>")

} else {

	stop("Sorry this specie isnt supported, try mm or hg please")

}

p_thr <- 0.1

cls <- lapply(des, function(d) {
		      d <- dplyr::mutate(d, 'reg' = ifelse(log2FoldChange > 0, 'up', 'down')) 
		      d <- split(d, d[['reg']])
		      lapply(d, function(dd) {
				     dplyr::filter(dd, log2FoldChange != 0 & padj < p_thr)[['Gene name']] %>%
					     unique
			})
	      })

cls <- unlist(cls, recursive = FALSE)
sapply(cls, length)

cls <- cls[which(sapply(cls, length) > 0)]

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

npath <- 10
p_thr <- 0.1

gg_enrch <- Reduce(rbind, 
		   lapply(names(enrch), function(cl) {
			   enr <- enrch[[cl]]
			   
			   enr <- Reduce(rbind, lapply(names(enr), function(path) {
							       enres <- enr[[path]]
							       if(nrow(enres) > 0) {
								       enres[['pathway']] <- path
								       head(enres, npath)
							       } else {
								       NULL
							       }
							})
			   		)

			   enr <- dplyr::filter(enr, Adjusted.P.value < p_thr)
			   if(nrow(enr) > 0) {
				   enr[['module']] <- paste0(cl)
				   enr
			   } else {
				   NULL
			   }
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


pdf(paste0(outputfolder, '/', jobname, '_deg_fisher_enrichment.pdf'),
    width = 15, height = 8)
	for(i in seq(ggs)) {
		plot(ggs[[i]])
	}
dev.off()

# ----------------------------------------------------------------------------------
# Fisher exact test using R base
# Testing intersection significancy
# ----------------------------------------------------------------------------------


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

fisher_exact_test <- function(m1, m2, universe) { # Function: Calculate the Fisher Exact Test (Diogenes)
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

names(gmts)
sapply(gmts, function(g) length(unique(g$term)))

names(cls)

minus_log10_pthr <- 1.3

mclapply(names(gmts), function(gmt_name) {
		 fisher_enrch(gmts[[gmt_name]], gmt_name, cls, minus_log10_pthr)
  }, mc.cores = 9)

fisher_enrch <- function(gmt, gmt_name, cls, minus_log10_pthr) {
	#
	refs <- split(gmt, gmt[['term']])
	refs <- lapply(refs, function(f) f[['gene']])
	igenes <- append(cls, refs)
	#
	smilmat <- calc_simil(igenes, "fisher")
	smilmat <- Reduce(rbind, smilmat)
	smilmat <- write_simil_matrix(smilmat)
	mat <- cbind(data.matrix(smilmat[names(refs), names(cls)[grep("down", names(cls))]])*-1, 
		     smilmat[names(refs), names(cls)[grep("up", names(cls))]])
	mat <- mat[, names(cls)]
	#
	write.table(mat, paste0(outputfolder, "/fisher_intersect_", gmt_name, "_", jobname, "_heatmap.tsv"),
		    sep = "\t", quote = FALSE, row.names = TRUE)
	#
	max_enrch <- apply(mat, 2, max)
	print(max_enrch)
	#
	max_ids <- unique(unlist(apply(mat, 2, function(m) which(abs(as.numeric(m)) >= minus_log10_pthr))))
	#
	mat <- mat[max_ids, ]

	int_tab <- Reduce(rbind,
			 lapply(rownames(mat), function(p) {
				Reduce(rbind,
				      lapply(colnames(mat), function(gs) {
					      data.frame('pathway' = p,
							 'gene_set' = gs,
							 'intersected_genes' = paste(intersect(refs[[p]],
											       cls[[gs]]),
								   		     collapse = "///")
							 )
						})
				      )
			})
		)

	int_tab[['path_geneset']] <- paste(int_tab[['pathway']], int_tab[['gene_set']], sep = "_")

	long_mat <- mat
	long_mat[['pathway']] <- rownames(long_mat)
	long_mat <- reshape2::melt(long_mat, id.vars = "pathway")
	long_mat[['path_geneset']] <- paste(long_mat[['pathway']], long_mat[['variable']], sep = "_")

	int_tab <- merge(int_tab, long_mat, by = "path_geneset")

	gg_heat <- ggplot(int_tab, aes(x = variable,  y = pathway.x, fill = as.numeric(value), text = intersected_genes)) + 
			geom_tile() +
			theme_classic() +
			scale_fill_gradient2(high = "#b2182b",  mid = "white", low = "#2166ac") +
			theme(axis.text.x = element_text(angle = 90)) +
			scale_y_discrete(limits = rownames(mat))


	htmlwidgets::saveWidget(plotly::as_widget(plotly::ggplotly(gg_heat)), 
            paste0(normalizePath(outputfolder), 
	           "/", jobname, "_", gmt_name, 
		   "_fisher_enrch.html"))
	
	col_fun = colorRamp2(c(-2, 0, 2), 
		     c("#2166ac", "white", "#b2182b" ))

	pdf(paste0(outputfolder, "/fisher_intersect_", gmt_name, "_", jobname, "_heatmap.pdf"), width = 6, height = 15)
		ht <- Heatmap(data.matrix(mat), 
		      cluster_rows = FALSE, 
		      cluster_columns = FALSE,
		      col = col_fun,
		      heatmap_legend_param = list(title_position = "topcenter",
					    title = "-log10(FDR)", 
					    color_bar = "continuous", 
					    legend_direction = "horizontal",
					    legend_width = unit(5, "cm")))
			draw(ht, heatmap_legend_side = "top")
	dev.off()
}

names(gmts) 
minus_log10_pthr <- 6 

mclapply(names(gmts)[c(3)], function(gmt_name) {
		 vis_fisher_enrch(gmts[[gmt_name]], gmt_name, cls, minus_log10_pthr)
  }, mc.cores = 9)

vis_fisher_enrch <- function(gmt, gmt_name, cls, minus_log10_pthr) {
	#
	refs <- split(gmt, gmt[['term']])
	refs <- lapply(refs, function(f) f[['gene']])
	igenes <- append(cls, refs)
	#
	#
	mat <- data.table::fread(paste0(outputfolder, "/fisher_intersect_", gmt_name, "_", jobname, "_heatmap.tsv"))
	mat <- data.frame(mat)
	rownames(mat) <- mat[['V1']]
	mat <- mat[, -1]
	print(head(mat,1))
	#
	max_enrch <- apply(mat, 2, max)
	#print(max_enrch)
	#
	max_ids <- unique(unlist(apply(mat, 2, function(m) which(abs(as.numeric(m)) >= minus_log10_pthr))))
	#
	mat <- mat[max_ids, ]

	int_tab <- Reduce(rbind,
			 lapply(rownames(mat), function(p) {
				Reduce(rbind,
				      lapply(colnames(mat), function(gs) {
					      data.frame('pathway' = p,
							 'gene_set' = gs,
							 'intersected_genes' = paste(intersect(refs[[p]],
											       cls[[gs]]),
								   		     collapse = "///")
							 )
						})
				      )
			})
		)

	int_tab[['path_geneset']] <- paste(int_tab[['pathway']], int_tab[['gene_set']], sep = "_")

	long_mat <- mat
	long_mat[['pathway']] <- rownames(long_mat)
	long_mat <- reshape2::melt(long_mat, id.vars = "pathway")
	long_mat[['path_geneset']] <- paste(long_mat[['pathway']], long_mat[['variable']], sep = "_")

	int_tab <- merge(int_tab, long_mat, by = "path_geneset")

	gg_heat <- ggplot(int_tab, aes(x = variable,  y = pathway.x, fill = as.numeric(value), text = intersected_genes)) + 
			geom_tile() +
			theme_classic() +
			scale_fill_gradient2(high = "#b2182b",  mid = "white", low = "#2166ac") +
			theme(axis.text.x = element_text(angle = 90)) +
			scale_y_discrete(limits = rownames(mat))


	ijobname <- paste(jobname, "minus_log10_p", gsub("\\.", "_", minus_log10_pthr), sep =  "_")
	htmlwidgets::saveWidget(plotly::as_widget(plotly::ggplotly(gg_heat)), 
            paste0(normalizePath(outputfolder), 
	           "/", ijobname, "_", gmt_name, 
		   "_fisher_enrch.html"))
	
	col_fun = colorRamp2(c(-2, 0, 2), 
		     c("#2166ac", "white", "#b2182b" ))

	pdf(paste0(outputfolder, "/fisher_intersect_", gmt_name, "_", ijobname, "_heatmap.pdf"), width = 6, height = 15)
		ht <- Heatmap(data.matrix(mat), 
		      cluster_rows = FALSE, 
		      cluster_columns = FALSE,
		      col = col_fun,
		      heatmap_legend_param = list(title_position = "topcenter",
					    title = "-log10(FDR)", 
					    color_bar = "continuous", 
					    legend_direction = "horizontal",
					    legend_width = unit(5, "cm")))
			draw(ht, heatmap_legend_side = "top")
	dev.off()
}
