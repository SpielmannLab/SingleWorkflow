#!/bin/bash

## init profile just once
if [ -z "$TX_HOME" ] ; then
  source /<path-to-env>/10xgenomics/profile.cellranger
fi

## basic
DATE=$(date +%Y-%m-%d)
THREADS=50
MEM=100

## sample specific

## call
cellranger mkref \
	--genome=GRCh38_nuclei \
	--fasta=/<path-to>/GRCh38/gencode/GRCh38.p13.genome.fa \
	--genes=/<path-to>/GRCh38/gencode/gencode.v33.annotation_filtered.gtf \
	--nthreads=$THREADS \
	--memgb=$MEM
