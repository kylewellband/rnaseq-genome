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

# Modules
module load samtools
module load star 
module load python

# Build reference
echo "Building ${GENOME} index..."

samtools faidx ${GENOME}

gatk --java-options '-Xmx4G' \
    CreateSequenceDictionary \
        -R ${GENOME} \
        -O ${GENOME%.*}.dict

BIT_SIZE=$(python3 01_scripts/util/calcSTARGenomeBitSize.py ${GENOME})

STAR --runMode genomeGenerate \
    --runThreadN ${NCPUS} \
    --genomeDir "02_reference/STARindex/" \
    --genomeFastaFiles ${GENOME} \
    --sjdbGTFfile ${GFF} \
    --sjdbGTFtagExonParentTranscript Parent \
    --sjdbGTFfeatureExon exon \
    --genomeChrBinNbits ${BIT_SIZE} \
    --genomeSAindexNbases 14


