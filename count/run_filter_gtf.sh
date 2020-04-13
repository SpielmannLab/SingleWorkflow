#!/bin/bash

## init profile just once
if [ -z "$TX_HOME" ] ; then
  source /<path-to-env>/10xgenomics/profile.cellranger
fi

## basic
DATE=$(date +%Y-%m-%d)
THREADS=40
MEM=256

## sample specific

## actual call
cellranger mkgtf /<path-to>/GRCh38/gencode/gencode.v33.annotation.gtf /<path-to>/gencode.v33.annotation_filtered.gtf \
                   --attribute=gene_type:protein_coding \
                   --attribute=gene_type:lncRNA \
                   --attribute=gene_type:processed_pseudogene \
                   --attribute=gene_type:transcribed_unprocessed_pseudogene \
                   --attribute=gene_type:unprocessed_pseudogene \
                   --attribute=gene_type:polymorphic_pseudogene \
                   --attribute=gene_type:misc_RNA \
                   --attribute=gene_type:snRNA \
                   --attribute=gene_type:miRNA \
                   --attribute=gene_type:transcribed_unitary_pseudogene \
                   --attribute=gene_type:transcribed_processed_pseudogene \
                   --attribute=gene_type:snoRNA \
                   --attribute=gene_type:scaRNA \
                   --attribute=gene_type:IG_LV_gene \
                   --attribute=gene_type:IG_V_gene \
                   --attribute=gene_type:IG_V_pseudogene \
                   --attribute=gene_type:IG_D_gene \
                   --attribute=gene_type:IG_J_gene \
                   --attribute=gene_type:IG_J_pseudogene \
                   --attribute=gene_type:IG_C_gene \
                   --attribute=gene_type:IG_C_pseudogene \
                   --attribute=gene_type:TR_V_gene \
                   --attribute=gene_type:TR_V_pseudogene \
                   --attribute=gene_type:TR_D_gene \
                   --attribute=gene_type:TR_J_gene \
                   --attribute=gene_type:TR_J_pseudogene \
                   --attribute=gene_type:TR_C_gene
