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
GENOME="02_reference/genome.fasta"  # Genomic reference.fasta
INPUT="06_stringtie_gtf"
OUTPUT="07_ballgown"
NCPUS=4
GFF="02_reference/genes.gtf"

# Modules
module load stringtie python


# Split off to a new script for quantification...
for file in $(ls "$INPUT"/*.bam | perl -pe 's/\.bam//g')
do
    name=$(basename $file)
    
    echo "Running stringtie -eB for sample: ${name}"
    stringtie -p ${NCPUS} -G $GFF -eB "${file}" > ${OUTPUT}/"${name}".ctab
    
    echo "${name}\t${OUTPUT}/${name}.ctab" >> ${OUTPUT}/ballgown_files.txt
done

#this is included in stringtie... does this work on Manitou?
prepDE.py -i ${OUTPUT}/ballgown_files.txt -g ${OUTPUT}/gene_count_matrix.csv -t ${OUTPUT}/transcript_count_matrix.csv

