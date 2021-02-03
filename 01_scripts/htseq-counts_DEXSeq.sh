#!/bin/bash

RBASE="/home/kylewellband/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts"
NCPUS=1
GFF="02_reference/genes.gff"
INPUT="05_aligned_bam"
OUTPUT="09_htseq-counts"

ls -1 ${INPUT}/*.bam |
perl -pe 's/(.*\/)(.*)(\.bam)/\2/' |
parallel -j $NCPUS python ${RBASE}/dexseq_count.py \
-p no \
-s no \
-f bam \
-r pos $GFF ${INPUT}/{}.bam ${OUTPUT}/{}.txt

