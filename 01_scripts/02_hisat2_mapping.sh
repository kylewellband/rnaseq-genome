#!/bin/bash
# 02_hisat2_mapping.sh
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
GENOME="02_reference/genome.fasta"  # Genomic reference.fasta
INPUT="04_trimmed_reads"
OUTPUT="05_aligned_bam"
NCPUS=4

# Modules
module load hisat2 samtools

# Align reads
for file in $(ls "$INPUT"/*_R1.fastq.gz | perl -pe 's/_R[12]\.fastq\.gz//g')
do
    name=$(basename $file)
    
    echo "Aligning $file"

    hisat2 -p $NCPUS \
        -x "$GENOME" \
        -k 20 \
        --dta \
        --rg-id "$name" \
        --rg "SM:$name" \
        # for single end data: uncomment the next line and comment the following three
        # -U "$name"_R1.fastq.gz |
        --no-mixed --no-discordant \
        -1 "$name"_R1.fastq.gz \
        -2 "$name"_R2.fastq.gz |
    samtools view -b -q5 |
    samtools sort -@ $NCPUS - -o "$OUTPUT"/"$name".bam
    
    samtools index "$OUTPUT"/"$name".bam

done &> | tee $LOG_FOLDER/"$TIMESTAMP"_mapping.log
