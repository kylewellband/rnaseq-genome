#!/bin/bash
# 00_index_STAR_reference.sh
# by: Kyle Wellband
# Sept. 8, 2019
#

# Copy script as it was run
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10_log_files"
cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# Global variables
GENOME="02_reference/genome.fasta"  # Genomic reference .fasta
GFF="02_reference/genes.gff"
NCPUS=4

# Build reference
echo "Building ${GENOME} index with SJs..."

BIT_SIZE=$(python3 01_scripts/util/calcSTARGenomeBitSize.py ${GENOME})
Nsjs=$(cat 05_aligned_bam/*SJ.out.tab | cut -f1-3 | sort | uniq | wc -l) 

STAR --runMode genomeGenerate \
    --runThreadN ${NCPUS} \
    --genomeDir "02_reference/STARindexSJ/" \
    --genomeFastaFiles ${GENOME} \
    --sjdbGTFfile ${GFF} \
    --sjdbGTFtagExonParentTranscript Parent \
    --sjdbGTFfeatureExon exon \
    --limitSjdbInsertNsj ${Nsjs} \
    --sjdbFileChrStartEnd 05_aligned_bam/*SJ.out.tab \
    --genomeChrBinNbits ${BIT_SIZE} \
    --genomeSAindexNbases 12


