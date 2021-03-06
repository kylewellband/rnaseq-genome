#!/bin/bash

GENOME="02_reference/genome.fasta"
TRANSCRIPTOME="02_reference/mrna.fasta"
NCPUS=4

grep "^>" $GENOME | cut -d " " -f 1 > 02_reference/decoys.txt
sed -i.bak -e 's/>//g' 02_reference/decoys.txt

if [[ ! -f $GENOME ]] || [[ ! -f $TRANSCRIPTOME ]]; then
    echo "Error: genome or transcriptome files not found!"
    exit
fi

cat $TRANSCRIPTOME $GENOME > 02_reference/gentrome.fasta

salmon index -t 02_reference/gentrome.fasta -d 02_reference/decoys.txt -p $NCPUS -i 02_reference/salmon_index
