#!/bin/bash
# 03_stringtie_assembly.sh
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
INPUT="05_aligned_bam"
OUTPUT="06_stringtie_gtf"
NCPUS=4

# Modules
module load stringtie

# Run StringTie
for file in $(ls "$INPUT"/*.sorted.out.bam | perl -pe 's/Aligned\.sorted\.out\.bam//g')
do
    name=$(basename $file)
    
    echo "Running stringtie for sample: ${name}"
    stringtie -p ${NCPUS} "${file}"Aligned.sorted.out.bam > $OUTPUT/"${name}".gtf
    
done

# Merge StringTie results
echo "Merging stringtie results to: 02_reference/stringtie_genes.gtf"
stringtie -p ${NCPUS} --merge ${OUTPUT}/*.gtf > 02_reference/stringtie_genes.gtf


