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
INPUT="05_aligned_bam"
OUTPUT="06_stringtie_gtf"
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

ls ${OUTPUT}/*.ctab | awk -F "/" '{sub(/.ctab/, "", $2); print $2 "\t" $1 "/" $2 ".ctab"}' > ${OUTPUT}/stringtie_ctab_files.txt


