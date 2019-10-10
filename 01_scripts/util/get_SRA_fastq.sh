#!/bin/bash
# get_SRA_fastq.sh
# by: Kyle Wellband
#

# Copy script as it was run
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10_log_files"
cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# Global variables
NCPUS=4
SRA_FILE="${1}"
OUTPUT="03_raw_data"

# Modules
module load SRA-toolkit

for sample in $(cat "${SRA_FILE}");
do

    echo "Fetching sample ${sample}..."
    fasterq-dump \
        --threads ${NCPUS} \
        --progress \
        -O ${OUTPUT} \
        ${sample}
    
    gzip ${OUTPUT}/${sample}.fastq

done
