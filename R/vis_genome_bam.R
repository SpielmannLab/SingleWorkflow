# Visualize CellRanger results and the reference transcriptome ## IN DEV : only input parameters, functioning-core ##

" Visualize CellRanger results and the reference transcriptome ## IN DEV : only input parameters, functioning-core ##

Usage: plot_igenes_seurat3.R --jobname=<value> [--mcafile=<file>] [--genefile=<file>] [--gmtfile=<file>] --specie=<value> --cellvar=<value> --resolution=<value> --infolder=<folder> --outfolder=<folder> --ncores=<value> 

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
  --infolder=<file>    Path cellranger /outs folder with bam files
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
	  'SummarizedExperiment', 'devtools', 'ggplot2',
          'cowplot', "dplyr", "Seurat", 
	  'data.table', 'Gviz')

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


jobname <- ""
`genome` <- "mm10"
bamfile <- "possorted_genome_bam.bam"

inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
jobname <- arguments$jobname
genefile <- arguments$genefile
res <- as.numeric(arguments$resolution)

# --- Functions

extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- strsplit(gtf_attributes, "; ")
  att <- gsub("\"","",unlist(att))
  if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return('NaN')}
}



# ------------------------------------------------------------------------------
# --- Loading mca onject
# ------------------------------------------------------------------------------

jobname <- ""
bamFile <- paste0(inputfolder, '/', bamfile)
gfrom <- 120537901
gfrom <- 120737901
gto <- 120959099
gfrom <- 120846901
gto <- 120849526


repeatmasker <- UcscTrack(genome="mm10", 
			  chromosome="chr1", 
			  track="rmsk", 
			  from=gfrom, 
			  to=gto,
			  trackType="GeneRegionTrack", 
			  rstarts="genoStart", 
			  rends="genoEnd", 
			  gene="repClass",
			  symbol="repName", 
			  transcript="repName", 
			  strand="strand", 
			  fill="#8282d2", 
			  name="UCSC repeat-masker",
			  transcriptAnnotation="symbol",
			  fill = "salmon")

repeatmasker2 <- UcscTrack(genome="mm10", 
			  chromosome="chr1", 
			  track="rmsk", 
			  from=gfrom, 
			  to=gto,
			  trackType="GeneRegionTrack", 
			  rstarts="genoStart", 
			  rends="genoEnd", 
			  gene="repClass",
			  symbol="repName", 
			  transcript="repName", 
			  strand="strand", 
			  fill="#8282d2", 
			  name="UCSC repeat-masker",
			  transcriptAnnotation="gene",
			  col = NULL,
			  fill = "white")

gencode23 <- UcscTrack(genome="mm10", 
		       	  chromosome="chr1", 
			  track="knownGene", 
			  from=gfrom, 
			  to=gto,
			  trackType="GeneRegionTrack", 
			  rstarts="exonStarts", 
			  rends="exonEnds", 
			  gene="name",
			  symbol="name", 
			  transcript="name", 
			  strand="strand", 
			  fill="#8282d2", 
			  name="gencode_v23",
			  transcriptAnnotation="gene",
			  col = NULL,
			  fill = "black")



options(ucscChromosomeNames=FALSE)

gtrack <- GenomeAxisTrack()

itrack <- IdeogramTrack(genome="mm10", 
			chromosome="chr1")

alTrack <- AlignmentsTrack(bamFile, 
			   isPaired=FALSE, 
			   genome="mm10",
			   chromosome="chr1",
			   )
alTrack


dTrack4 <- DataTrack(range=bamFile, 
		     genome="mm10", 
		     type="l", 
		     name="Coverage", 
		     window=-1, 
		     chromosome="chr1")

ht <- HighlightTrack(trackList=list(repeatmasker2, repeatmasker, reftr, gtrack, dTrack4, alTrack), 
		     start=120847926, width=284, chromosome="chr1", col="white")


gtfile <- "genes.gtf"
gtf <- fread(gtfile, skip="##", header = FALSE, quote = "") %>% filter(V3 == 'exon')

head(gtf)
colnames(gtf) <- c('chromosome', 'db', 'type', 'start', 'end', 'V6', 'strand', 'V8', 'V9')

sgtf <- dplyr::filter(gtf, chromosome == "chr1" & start > gfrom & end < gto)
dim(sgtf)

sgtf[['gene']] <- unlist(lapply(sgtf[['V9']], function(g) extract_attributes(g, "gene_id")))
sgtf[['symbol']] <- unlist(lapply(sgtf[['V9']], function(g) extract_attributes(g, "gene_name")))
sgtf[['exon']] <- unlist(lapply(sgtf[['V9']], function(g) extract_attributes(g, "exon_id")))
sgtf[['transcript']] <- unlist(lapply(sgtf[['V9']], function(g) extract_attributes(g, "transcript_id")))
sgtf[['feature']] <- unlist(lapply(sgtf[['V9']], function(g) extract_attributes(g, "gene_type")))


sgtf <- sgtf %>% select(chromosome, start, end, strand, feature, gene, exon, transcript, symbol)
head(sgtf)

reftr <- GeneRegionTrack(sgtf,
			 genome="mm10",
			 transcriptAnnotation="symbol", 
			 chromosome="chr1")
head(reftr)

pdf(paste0(outputfolder, '/', jobname, "track_bam_ref.pdf"), width = 14, height = 22)
	plotTracks(list(itrack, repeatmasker2, repeatmasker, gtrack, reftr, dTrack4, alTrack), chromosome="chr1", from=gfrom, to=gto)
dev.off()

pdf(paste0(outputfolder, '/', jobname, "track_bam_ref_highlighted.pdf"), width = 14, height = 22)
	plotTracks(list(itrack, gencode23, ht), chromosome="chr1", from=gfrom, to=gto)
dev.off()
pdf(paste0(outputfolder, '/', jobname, "track_bam_ref.pdf"), width = 14, height = 22)
	plotTracks(list(itrack,  gtrack, gencode23,reftr, alTrack), chromosome="chr1", from=gfrom, to=gto)
dev.off()

