#!/bin/bash

# sending cellranger count jobs to the cluster (mxqsub) 

OUTDIR="/<path-to-write-results>"
REFDIR="/<path-to-ref-transcriptome>"

mxqsub --tmpdir 500G --stdout=. --stderr=. --threads=20 --memory=300000 --runtime=24h ./cellranger_count_jobscluster.sh sc1 $REFDIR <path-to-fasq-folder-sc1> $OUTDIR
mxqsub --tmpdir 500G --stdout=. --stderr=. --threads=20 --memory=300000 --runtime=24h ./cellranger_count_jobscluster.sh sc2 $REFDIR <path-to-fasq-folder-sc2> $OUTDIR
mxqsub --tmpdir 500G --stdout=. --stderr=. --threads=20 --memory=300000 --runtime=24h ./cellranger_count_jobscluster.sh sc... $REFDIR <path-to-fasq-folder-sc...> $OUTDIR
