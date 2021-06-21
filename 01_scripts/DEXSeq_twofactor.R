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
                          lake = samples$lake,
                          condition = samples$form)
fullModel <- ~ sample + exon + condition:exon + lake:exon
reducedModel_form <- ~ sample + exon + lake:exon
reducedModel_lake <- ~ sample + exon + form:exon

dxd <- DEXSeqDataSetFromHTSeq(
    countFiles,
    sampleData = sampleTable,
    design = fullModel,
    flattenedfile = flattenedFile)

dxd <- estimateSizeFactors(dxd)

dxd <- estimateDispersions(dxd, formula = fullModel, BPPARAM = BPPARAM)

#plotDispEsts(dxd)

dxd_form <- testForDEU(dxd, fullModel = fullModel, reducedModel = reducedModel_form, BPPARAM = BPPARAM)
dxd_lake <- testForDEU(dxd, fullModel - fullModel, reducedModel = reducedModel_lake, BPPARAM = BPPARAM)
dxd_form <- estimateExonFoldChanges(dxd_form, fitExpToVar="condition", BPPARAM = BPPARAM)
dxd_lake <- estimateExonFoldChanges(dxd_lake, fitExpToVar="lake", BPPARAM = BPPARAM)
dxr_form <- DEXSeqResults(dxd_form)
dxr_lake <- DEXSeqResults(dxd_lake)

save(dxd, dxd_form, dxd_lake, dxr_form, dxr_lake, file = paste0(basename(args[1]), ".Rdata"))
