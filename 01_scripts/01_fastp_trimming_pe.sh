#!/bin/bash

# 4 CPU
# 10 Go

# Copy script as it was run
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10_log_files"
cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# Global variables
LENGTH=80 # This is good for PE100, consider increasing for longer reads (e.g. PE125, PE150)
QUAL=25
INPUT="03_raw_data"
OUTPUT="04_trimmed_reads"
NCPUS=4

# Trim reads with fastp
for file in $(ls "$INPUT"/*_R1.fastq.gz | perl -pe 's/_R[12]\.fastq\.gz//g')
do
    name=$(basename $file)

    # Fastp
    fastp -w $NCPUS \
        -i "$file"_R1.fastq.gz \
        -I "$file"_R2.fastq.gz \
        -o $OUTPUT/"$name"_R1.fastq.gz \
        -O $OUTPUT/"$name"_R2.fastq.gz  \
        --length_required="$LENGTH" \
        --qualified_quality_phred="$QUAL" \
        --trim_poly_x \
	--correction \
        --trim_tail1=1 \
        --trim_tail2=1 \
        --json $OUTPUT/"$name".json \
        --html $OUTPUT/"$name".html  \
        --report_title="$name"report.html

done 2>&1 | tee 10_log_files/"$TIMESTAMP"_fastp.log
