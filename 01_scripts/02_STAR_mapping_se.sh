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

## Need to adjust settings in /etc/sysctl.conf to take advantage of shared memory
## See STAR manual for details
# STAR --runThreadN $NCPUS \
#     --genomeDir 02_reference/STARindex/ \
#     --genomeLoad LoadAndExit

# Align reads first pass
for file in $(ls "$INPUT"/*.fastq.gz | perl -pe 's/\.fastq\.gz//g')
do
    name=$(basename $file)
    
    echo -e "\nFirst Pass: Aligning $file"

    STAR --runThreadN ${NCPUS} \
        --genomeDir 02_reference/STARindex/ \
        --readFilesIn ${INPUT}/${name}.fastq.gz \
        --readFilesCommand gunzip -c \
        --twopassMode None \
        --outSAMmapqUnique 60 \
        --outFileNamePrefix ${OUTPUT}/${name} \
        --outSAMtype BAM Unsorted \
        --outFilterMultimapNmax 20 \
        --outSAMattrRGline ID:${name} SN:${name} PL:ILLUMINA

    rm "$OUTPUT"/"$name"Aligned.out.bam

done | tee $LOG_FOLDER/"$TIMESTAMP"_firstpass_mapping.log

# Align reads second pass inserting SJ on the fly
for file in $(ls "$INPUT"/*.fastq.gz | perl -pe 's/\.fastq\.gz//g')
do
    name=$(basename $file)
    
    echo -e "\nSecond Pass: Aligning $file"

    STAR --runThreadN ${NCPUS} \
        --genomeDir 02_reference/STARindex/ \
        --readFilesIn ${INPUT}/${name}.fastq.gz \
        --readFilesCommand gunzip -c \
        --twopassMode None \
        --sjdbFileChrStartEnd ${OUTPUT}/*SJ.out.tab \
        --limitSjdbInsertNsj=2000000 \
        --outSAMmapqUnique 60 \
        --outFileNamePrefix ${OUTPUT}/${name} \
        --outSAMtype BAM Unsorted \
        --outFilterMultimapNmax 20 \
        --quantMode TranscriptomeSAM GeneCounts \
        --outSAMattributes NH HI AS nM \
        --outSAMattrRGline ID:${name} SM:${name} PL:ILLUMINA

    samtools sort -@ $NCPUS "$OUTPUT"/"$name"Aligned.out.bam > "$OUTPUT"/"$name"Aligned.sorted.out.bam
    rm "$OUTPUT"/"$name"Aligned.out.bam
    samtools index "$OUTPUT"/"$name"Aligned.sorted.out.bam

done | tee $LOG_FOLDER/"$TIMESTAMP"_secondpass_mapping.log

# STAR --runThreadN $NCPUS \
#     --genomeDir 02_reference/STARindex/ \
#     --genomeLoad Remove

