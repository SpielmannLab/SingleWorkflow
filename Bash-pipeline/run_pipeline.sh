#!/bin/bash


# ---------------------
# Reading raw data (CellRanger outs folder)
# ---------------------

jobname="<cool-descriptive-string>"
outfolder="<path-to-writing output>"
sp="hg" # "mm"

IFS=$'\n'

cr2cds=( # per cellranger count run/sample
	"Rscript read_count_cellranger.R --infolder=<path1> --outfolder=${outfolder} --jobname="" --genome=${sp}"
	"Rscript read_count_cellranger.R --infolder=<path2> --outfolder=${outfolder} --jobname="" --genome=${sp}"
	"Rscript read_count_cellranger.R --infolder=<...> --outfolder=${outfolder} --jobname="" --genome=${sp}"
	)

parallel {} ::: ${cr2cds[*]}

# ---------------------
# QC
# ---------------------

infolder="path-to-read_count_cellranger.R-output"
dataset="RAW"

qcreport=(
	"Rscript qc_sc.R --scobject=sc1.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript qc_sc.R --scobject=sc2.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript qc_sc.R --scobject=s....rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	)

parallel {} ::: ${qcreport[*]}

# ---------------------
# FILTERING
# ---------------------

minfeat=1
maxfeat="Inf"
mincellfeat=3
maxcellfeat="Inf"
mincountcell=500
maxcountcell=20000
mingenecell=200
maxgenecell=5000
pctmt=0.1
pctrb=0.1
mt="FALSE"
rb="FALSE"

filter=(
	"Rscript filter_qc.R --scobject=sc1_RAW_qc.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --minfeat=${minfeat} --maxfeat=${maxfeat} --mincellfeat=${mincellfeat} --maxcellfeat=${maxcellfeat} --mincountcell=${mincountcell} --maxcountcell=${maxcountcell} --mingenecell=${mingenecell} --maxgenecell=${maxgenecell} --pctmt=${pctmt} --pctrb=${pctrb} --mt=${mt} --rb=${rb}"
	"Rscript filter_qc.R --scobject=sc2_RAW_qc.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --minfeat=${minfeat} --maxfeat=${maxfeat} --mincellfeat=${mincellfeat} --maxcellfeat=${maxcellfeat} --mincountcell=${mincountcell} --maxcountcell=${maxcountcell} --mingenecell=${mingenecell} --maxgenecell=${maxgenecell} --pctmt=${pctmt} --pctrb=${pctrb} --mt=${mt} --rb=${rb}"
	"Rscript filter_qc.R --scobject=..._RAW_qc.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --minfeat=${minfeat} --maxfeat=${maxfeat} --mincellfeat=${mincellfeat} --maxcellfeat=${maxcellfeat} --mincountcell=${mincountcell} --maxcountcell=${maxcountcell} --mingenecell=${mingenecell} --maxgenecell=${maxgenecell} --pctmt=${pctmt} --pctrb=${pctrb} --mt=${mt} --rb=${rb}"
	)

parallel {} ::: ${filter[*]}

# -----------------------
# Duplets score Scrublet
# -----------------------

npcs=15

duplets=(
	"python3 run_scrublet.py sc1 ${infolder} ${outfolder} ${npcs}"
	"python3 run_scrublet.py sc2 ${infolder} ${outfolder} ${npcs}"
	"python3 run_scrublet.py ... ${infolder} ${outfolder} ${npcs}"
	)

source /<path-to-virtenv>/scrublet/bin/activate
parallel {} ::: ${duplets[*]}
deactivate

# ---------------------
# Filter duplets given a duplet-score Scrublet threshold
# ---------------------

dataset="FILTERED"
dupthr=0.15


filterdup=(
	"Rscript filter_duplets.R --scobject=sc1_filtered.rds --dupthr=${dupthr} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript filter_duplets.R --scobject=sc2_filtered.rds --dupthr=${dupthr} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript filter_duplets.R --scobject=sc..._filtered.rds --dupthr=${dupthr} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	)

parallel {} ::: ${filterdup[*]}

# ---------------------
# QC after filter
# ---------------------

dataset="FILTERED"

qcreport=(
	# PD
	"Rscript qc_sc.R --scobject=sc1_filtered.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript qc_sc.R --scobject=sc2_filtered.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript qc_sc.R --scobject=sc..._filtered.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	)

parallel {} ::: ${qcreport[*]}

# --------------------------------------------------------------
# ---------------------
#                            Individual samples
# ---------------------
# --------------------------------------------------------------

# Clustering
ncor=5
npcs=50 # testing the amount of PCs to use
Rscript cluster_seurat3.R --jobname=${jobname} --mcafile=sc1_filtered_FILTERED_qc_seurat3.rds --npcs=${npcs} --testpcs=TRUE --infolder=${infolder} --outfolder=${outfolder} --ncores=${ncor}
npcs=15 # chosen number of PCs from previous step
Rscript cluster_seurat3.R --jobname=${jobname} --mcafile=CR_E9_5Black6FLs_filtered_FILTERED_qc_seurat3.rds --npcs=${npcs} --testpcs=FALSE --infolder=${infolder} --outfolder=${outfolder} --ncores=${ncor}
# Annotation
mca_file_cluster="sc1_filtered_FILTERED_qc_seurat3_seurat3_merged_clustered_resrange.rds"
res=0.1
Rscript /project/zvi/cprada/src_monocle3/get_mkr_genes_seurat3.R --jobname=${jobname} --mcafile=${mca_file_cluster} --resolution=${res} --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --ncores=${ncor}

# --------------------------------------------------------------
# ---------------------
# Merge samples
# ---------------------
# --------------------------------------------------------------
ncor=3
npcs=15

# monacle3 
Rscript merge_sc_monocle3.R --jobname=${jobname} --npcs=${npcs} --ncores=${ncor} --infolder=${infolder} --outfolder=${outfolder} sc1_filtered_FILTERED_qc.rds sc2_filtered_FILTERED_qc.rds sc..._filtered_FILTERED_qc.rds 
# seurat3 step1
Rscript merge_sc_seurat3_step1.R --jobname=${jobname} --npcs=${npcs} --ncores=${ncor} --testpcs=TRUE --infolder=${infolder} --outfolder=${outfolder} sc1_filtered_FILTERED_qc.rds sc2_filtered_FILTERED_qc.rds sc..._filtered_FILTERED_qc.rds 
#
Rscript merge_sc_seurat3_step1_5.R --jobname=${jobname} --npcs=${npcs} --ncores=${ncor} --testpcs=TRUE --infolder=${infolder} --outfolder=${outfolder} sc1_filtered_FILTERED_qc_seurat3.rds sc2_filtered_FILTERED_qc_seurat3.rds sc..._filtered_FILTERED_qc_seurat.rds 
# seurat3 step2
Rscript merge_sc_seurat3_step2.R --jobname=${jobname} --npcs=${npcs} --ncores=${ncor} --testpcs=TRUE --infolder=${infolder} --outfolder=${outfolder}

# --------------------------------------------------
# Clustering
# --------------------------------------------------

npcs=50 # testing the amount of PCs to use
Rscript cluster_seurat3.R --jobname=${jobname} --npcs=${npcs} --testpcs=TRUE --infolder=${infolder} --outfolder=${outfolder} --ncores=${ncor}
# 
npcs=15 # chosen number of PCs from previous step
Rscript /project/zvi/cprada/src_monocle3/cluster_seurat3.R --jobname=${jobname} --npcs=${npcs} --testpcs=FALSE --infolder=${infolder} --outfolder=${outfolder} --ncores=${ncor}
#
Rscript /project/zvi/cprada/src_monocle3/cluster_monocle3.R --jobname=${jobname} --npcs=${npcs} --infolder=${infolder} --outfolder=${outfolder} --ncores=${ncor}

# --------------------------------------------------
# Marker gene identification and functional-enrichment
# --------------------------------------------------

res=0.01
Rscript get_mkr_genes_seurat3.R --jobname=${jobname} --resolution=${res} --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --ncores=${ncor}

# --------------------------------------------------
# Add colData variables (project specific script)
# --------------------------------------------------

rdsfile=${jobname}"_seurat3_merged_clustered_resrange.rds"
Rscript create_cellmetadata_variable.R --rdsfile=${rdsfile} --infolder=${infolder} --outfolder=${outfolder}

# --------------------------------------------------
# Differential cellular composition (two-level variables)
# --------------------------------------------------

Rscript diff_cellcomp.R --jobname=${jobname} --groupvar=<condition> --resolution=<number> --infolder=${infolder} --outfolder=${outfolder}

# --------------------------------------------------
# Differential gene expression (two-level variables)
# --------------------------------------------------

Rscript /project/zvi/cprada/src_monocle3/diff_exp.R --jobname=${jobname} --groupvar=<condition> --samplevar=<sample> --resolution=${res} --specie=${sp} --infolder=${infolder} --outfolder=${outfolder}

##################################
# Data imputation and spearman correlation analysis
##################################

Rscript /project/zvi/cprada/src_monocle3/inputdata_magic.R --jobname=${jobname} --rdsfile=${rdsfile} --npcs=${npcs} --genequant=0.25 --nclusters=7 --specie=${sp} --npath=5 --infolder=${infolder} --outfolder=${outfolder} --ncores=${ncor}

##################################
# Sub-clustering
##################################

mkdir ${outfolder}"/"${jobname}

npcs=15
nvgenes=500 # min. number most variable genes
mincell=100 # min number of cells per sample in each sub-cluster to be considered

Rscript sub_cluster_seurat3.R --jobname=${jobname} --npcs=${npcs} --testpcs=TRUE --cellvar=<ref-cell-metadata> --samplevar=sample --nvgenes=500 --mincells=100 --infolder=${infolder} --outfolder=${outfolder} --ncores=${ncor}

soutfolder=${outfolder}"/"${jobname}
sinfolder=${infolder}"/"${jobname}
ncor=2
res=0.04 # clusering resolution

mcafiles=`find ${outfolder}"/"${jobname} -name "*_Seurat3_preprocessed_sample_integrated_subcluster.rds" -printf '%P\n'`

declare -a cluster
for i in ${mcafiles[@]}
do

   f="Rscript cluster_seurat3.R --jobname=${jobname} --mcafile="${i}" --npcs=${npcs} --testpcs=FALSE --infolder=${sinfolder} --outfolder=${soutfolder} --ncores=${ncor}"
   cluster+=(${f})
   echo $f 
   
done

parallel {} ::: ${cluster[*]}

mcafiles=`find ${outfolder}"/"${jobname} -name "*_Seurat3_preprocessed_sample_integrated_subcluster_clustered_resrange.rds" -printf '%P\n'`

declare -a ann_mkrs
for i in ${mcafiles[@]}
do

   f="Rscript get_mkr_genes_seurat3.R --jobname=${jobname} --mcafile="${i}" --resolution=${res} --specie=${sp} --infolder=${sinfolder} --outfolder=${soutfolder} --ncores=${ncor}"
   ann_mkrs+=(${f})
   echo $f 
   
done

parallel {} ::: ${ann_mkrs[*]}

##################################
# pseudo-ordering sub-clusters:
##################################

ncor=5
res=0.04

mcafiles=`find ${outfolder}"/"${jobname} -name "*_Seurat3_preprocessed_sample_integrated_subcluster_clustered_resrange.rds" -printf '%P\n'`

declare -a pseudotmp
for i in ${mcafiles[@]}
do

   f="Rscript order_cells_seurat3_slingshot.R --jobname=${jobname} --cdsfile=${i} --specie=${sp} --infolder=${sinfolder} --outfolder=${soutfolder} --ncores=${ncor} --res=${res}"
   pseudotmp+=(${f})
   echo $f 
   
done

parallel {} ::: ${pseudotmp[*]}
# here we go...
