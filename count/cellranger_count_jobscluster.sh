#! /usr/bin/bash

if [ $# -eq 4 ]; then
	SAMPLE="$1"
	TRANSCRIPTOME="$2"
	FASQ_DIR="$3"
	OUTFOLDER="$4"
else
	echo "Usage: $0 <SAMPLE_NAME> <REF_TRANSCRIP_FOLDER> <FASQ_FOLDER> <WHERE_TO_WRITE_RESULTS>"
	exit
fi

set -ve

ORIGDIR="$(/bin/pwd)"

cp -ar $FASQ_DIR/ $MXQ_JOB_TMPDIR/

cd $MXQ_JOB_TMPDIR

FASTQ_SET=$(ls)

echo $FASTQ_SET

## init profile just once
if [ -z "$TX_HOME" ] ; then
	source /<path-to-env>/10xgenomics/profile.cellranger
fi

DATE=$(date +%Y-%m-%d)
THREADS=20
MEM=200

cellranger count \
	--id=cr_count_${SAMPLE}_${DATE} \
	--fastqs=${MXQ_JOB_TMPDIR}/${FASTQ_SET} \
	--transcriptome=${TRANSCRIPTOME} \
	--localcores=${THREADS} \
	--localmem=${MEM} \
	--disable-ui

RESULT=$(ls | grep "cr_count")

cp -r $MXQ_JOB_TMPDIR/$RESULT $OUTFOLDER/
