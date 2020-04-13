#!/bin/bash

# turn transcripts into exons for single-nuclei mapping

awk 'BEGIN{FS="\t"; OFS="\t"} $3 == "transcript"{ $3="exon"; print}' \
       genes.gtf > genes_premrna.gtf
