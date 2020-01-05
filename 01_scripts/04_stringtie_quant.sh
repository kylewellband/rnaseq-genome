#!/bin/bash
# 04_stringtie_quant.sh
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
INPUT="06_stringtie_gtf"
OUTPUT="07_ballgown"
NCPUS=4
GFF="02_reference/stringtie_genes.gtf"

# Modules
module load stringtie
module load python

# Split off to a new script for quantification...
for file in $(ls "$INPUT"/*.sorted.out.bam | perl -pe 's/Aligned\.sorted\.out\.bam//g')
do
    name=$(basename $file)
    
    echo "Running stringtie -eB for sample: ${name}"
    stringtie -p ${NCPUS} -G $GFF -eB "${file}"Aligned.sorted.out.bam > ${OUTPUT}/"${name}".ctab
    
done

ls ${OUTPUT}/*.ctab > ${OUTPUT}/ballgown_files.txt

01_scripts/util/prepDE.py -i ${OUTPUT}/ballgown_files.txt -g ${OUTPUT}/gene_count_matrix.csv -t ${OUTPUT}/transcript_count_matrix.csv

