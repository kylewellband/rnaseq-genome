#!/bin/bash
# build_hisat2_reference.sh
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
NCPUS=4

# Modules
module load hisat2 

# Build reference
echo "Building ${GENOME} index..."
hisat2-build -p $NCPUS "$GENOME" "$GENOME"


