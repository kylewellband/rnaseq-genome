#!/bin/bash
# 02_STAR_mapping.sh
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
module load samtools
module load star

# Align reads
for file in $(ls "$INPUT"/*_R1.fastq.gz | perl -pe 's/_R[12]\.fastq\.gz//g')
do
    name=$(basename $file)
    
    echo "Aligning $file"

    STAR --runThreadN ${NCPUS} \
        --genomeDir ./genome \
        --readFilesIn ${INPUT}/${name}_R1.fastq.gz ${INPUT}/${name}_R2.fastq.gz \
        --readFilesCommand gunzip -c \
        --twopassMode Basic \
        --outSAMmapqUnique 60 \
        --outFileNamePrefix ${OUTPUT}/${name} \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterMultimapNmax 20 \
        --waspOutputMode SAMtag \
        --quantMode TranscriptomeSAM GeneCounts \
        --outSAMattributes NH HI AS nM vW \
        --outSAMattrRGline ID:${name} SN:${name} PL:ILLUMINA
    
    samtools index "$OUTPUT"/"$name".bam

done &> | tee $LOG_FOLDER/"$TIMESTAMP"_mapping.log
