#!/bin/bash

RBASE="/Users/kyle/Library/R/3.5/library/DEXSeq/python_scripts"
NCPUS=8
GFF="02_reference/genes_DEXSeq.gff"
INPUT="05_aligned_bam"
OUTPUT="09_htseq-counts"

ls -1 ${INPUT}/*.sorted.out.bam |
perl -pe 's/(.*\/)(.*)(\.bam)/\2/' |
parallel -j $NCPUS ~/miniconda3/bin/python ${RBASE}/dexseq_count.py -p no -s no -f bam -r pos $GFF ${INPUT}/{}.bam ${OUTPUT}/{}.txt

