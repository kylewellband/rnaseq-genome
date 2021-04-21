#!/bin/bash

NCPUS=8
GFF="02_reference/genes.gff"
INPUT="05_aligned_bam"
OUTPUT="09_htseq-counts"

ls -1 ${INPUT}/*.sorted.out.bam |
perl -pe 's/(.*\/)(.*)(\.bam)/\2/' |
parallel -j $NCPUS "htseq-count -s no -f bam -r pos -t CDS -i Parent ${INPUT}/{}.bam $GFF > ${OUTPUT}/{}_genes.txt"

