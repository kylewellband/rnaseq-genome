#!/usr/bin/env Rscript

library(DEXSeq)
library(BiocParallel)

# setwd("/mnt/HDD1/coregone/EU/rnaseq-genome/")

# inDir <- "09_htseq-counts"
BPPARAM = MulticoreParam(8)

# countFiles <- list.files(inDir, pattern=".txt$", full.names=TRUE)
# basename(countFiles)

flattenedFile <- "02_reference/genes_DEXSeq.gff"

args <- commandArgs(TRUE) # args <- "sample_info.txt"

samples <- read.table(args[1], header = TRUE, stringsAsFactors = FALSE)
# samples <- samples[order(samples$fastq_name),]

countFiles <- paste0("09_htseq-counts/", samples$fastq_name, "_R1Aligned.sorted.out.txt")

sampleTable <- data.frame(row.names = paste0(samples$fastq_name),
                          condition = samples$form)

dxd <- DEXSeqDataSetFromHTSeq(
    countFiles,
    sampleData = sampleTable,
    design = ~ sample + exon + condition:exon,
    flattenedfile = flattenedFile)

dxd <- estimateSizeFactors(dxd)

dxd <- estimateDispersions(dxd, BPPARAM = BPPARAM)

plotDispEsts(dxd)

dxd <- testForDEU(dxd, BPPARAM = BPPARAM)

dxd <- estimateExonFoldChanges(dxd, fitExpToVar="condition", BPPARAM = BPPARAM)

dxr <- DEXSeqResults(dxd)

save(dxd, dxr, file = paste0(basename(args[1]), ".Rdata"))
