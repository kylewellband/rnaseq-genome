#!/bin/bash
# 02_salmon_SA_pe.sh
# by: Kyle Wellband
# Jan. 20, 2020
#

# Copy script as it was run
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10_log_files"
cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# Global variables
INDEX="02_reference/salmon_index"
INPUT="04_trimmed_reads"
OUTPUT="07_salmon_quant"
NCPUS=4

# Modules
module load salmon

# Psuedo Align 
for file in $(ls "$INPUT"/*_R1.fastq.gz | perl -pe 's/_R1\.fastq\.gz//g')
do
    name=$(basename $file)
    
    echo -e "\nPseudo aligning $file"

    salmon quant \
        -i $INDEX \
        -l A \
        -1 $INPUT/"$name"_R1.fastq.gz \
        -2 $INPUT/"$name"_R2.fastq.gz \
        -p $NCPUS \
        --validateMappings \
	--gcBias \
	--seqBias \
	-o $OUTPUT/"$name"_quant

done | tee $LOG_FOLDER/"$TIMESTAMP"_salmon_quant.log



