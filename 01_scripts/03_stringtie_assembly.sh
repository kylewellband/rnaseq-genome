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
INPUT="05_aligned_bam"
OUTPUT="06_stringtie_gtf"
NCPUS=4
#why do I need this?
GFF="~/Projects/sfo_eqtl/01data/sfon_gene_salp_pos.gff"

# Modules
module load stringtie

# Run StringTie
for file in $(ls "$INPUT"/*.bam | perl -pe 's/\.bam//g')
do
    name=$(basename $file)
    
    echo "Running stringtie for sample: ${name}"
    stringtie -p ${NCPUS} "${file}" > $OUTPUT/"${name}".gtf
    
done

# Merge StringTie results
echo "Merging stringtie results and copying to: 02_reference/genes.gtf"
stringtie -p ${NCPUS} --merge ${OUTPUT}/*.gtf > ${OUTPUT}/merged_transcripts.gtf

cp ${OUTPUT}/merged_transcripts.gtf 02_reference/genes.gtf

