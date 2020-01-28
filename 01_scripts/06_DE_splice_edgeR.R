#!/usr/bin/env Rscript
## Working script for DML / DMR quantification using DSS

## Install and load necessary packages
for (p in c("data.table", "BiocManager", "Rsubread", "edgeR", "tidyverse", "GenomicFeatures")) {
    if (!suppressMessages(require(p, character.only = T))) {
        message(paste("Installing:", p))
        if(p %in% c("tximport", "edgeR", "GenomicFeatures")) {
            BiocManager::install(p)
        } else {
            install.packages(p, repos = "https://mirror.its.dal.ca/cran", dependencies = T)
        }
        suppressMessages(require(p, character.only = T))}
    rm(p)
}

args <- commandArgs(T)
# args <- c(1, "PE", "samples_splice.txt") ; setwd("~/Projects/tetrahymena/rnaseq-genome/")

gff <- "02_reference/genes.gff"

if (length(args) != 3)
    stop("Usage: 06_DE_splice_edgeR.R <n_cores> <PE/SR> <samples.txt>")

n_cores <- as.integer(args[1])

if (!args[2] %in% c("SR", "PE"))
    stop(paste0("option \"", args[2], "\" not recognized. Please specify \"SR\" or \"PE\""))

if (args[2] == "PE") {
    PE <- TRUE
} else {
    PE <- FALSE
}

samples <- read.table(args[3], header = T, stringsAsFactors = FALSE)

bam_files <- samples$file
names(bam_files) <- samples$sample

if (!(file.exists("07_DE_results/exon_feature_counts.rds"))) {
    fc <- featureCounts(bam_files, annot.ext = gff,
                        isGTFAnnotationFile = TRUE, GTF.featureType = "exon", GTF.attrType = "ID",
                        useMetaFeatures = FALSE, allowMultiOverlap = TRUE,
                        isPairedEnd = PE, nthreads = n_cores, strandSpecific = TRUE)
    saveRDS(fc, "07_DE_results/exon_feature_counts.rds", compress = "gzip")
} else {
    fc <- readRDS("07_DE_results/exon_feature_counts.rds")
}

# Specify unique groups and experimental design
group <- samples$group
design <- model.matrix(~group)

# Make DGE
y <- DGEList(counts = fc$counts, 
                 genes = fc$annotation, 
                 samples = samples,
                 group = group)
dim(y)
head(y$genes)

keep <- filterByExpr(y)
table(keep)

y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples

y <- estimateDisp(y, design = design)
plotBCV(y)

fit <- glmQLFit(y, design)
plotQLDisp(fit)

qlf <- glmQLFTest(fit)
topTags(qlf)

is_de <- decideTests(qlf, p.value = 0.05, lfc = log2(1.5))
summary(is_de)

fit$genes$ExonID <- fit$genes$GeneID
fit$genes$GeneID <- sub(".exon[0-9]+", "", fit$genes$GeneID)
sp <- diffSpliceDGE(fit, coef=2, geneid="GeneID", exonid="ExonID")
dim(topSpliceDGE(sp, test="gene", n=100000))
sum(topSpliceDGE(sp, test="gene", n=100000)$FDR < 0.05)

