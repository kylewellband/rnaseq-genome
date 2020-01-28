#!/usr/bin/env Rscript
## Working script for DML / DMR quantification using DSS

## Install and load necessary packages
for (p in c("data.table", "BiocManager", "tximport", "edgeR", "tidyverse", "GenomicFeatures")) {
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
# args <- c("samples.txt", "known") ; setwd("~/Projects/tetrahymena/rnaseq-genome/")
 

## Sanity checking
if (length(args) != 2)
    stop("Usage: tximport_edgeR.R <samples.txt> <known/denovo>")

if (!args[2] %in% c("known", "denovo"))
    stop(paste0("option \"", args[2], "\" not recognized. Please specify \"denovo\" or \"known\""))

if (args[2] == "denovo") {
    denovo <- TRUE
} else {
    denovo <- FALSE
}

samples <- read.table(args[1], header = T, stringsAsFactors = FALSE)

#### Set up group factor and design matrix
group <- samples$group
design <- model.matrix(~ samples$group)

if(!all(c("sample", "file", formula_parts) %in% colnames(samples)))
    stop("Samples file must contain a header row with names: \'sample\', \'file\', and factor(s).")

## Create and save, or load transcript to gene DB
if (denovo) {
    if (!file.exists("02_reference/genes_denovo.sqlite")) {
        txdb <- makeTxDbFromGFF(file = "02_reference/stringtie_genes.gtf", format = "gtf")
        saveDb(txdb, "02_reference/genes_denovo.sqlite")
    } else {
        txdb <- loadDb("02_reference/genes_denovo.sqlite")
    }
} else {
    if (!file.exists("02_reference/genes_known.sqlite")) {
        txdb <- makeTxDbFromGFF(file = "02_reference/genes.gff", format = "gff")
        saveDb(txdb, "02_reference/genes_known.sqlite")
    } else {
        txdb <- loadDb("02_reference/genes_known.sqlite")
    }
}

# files <- paste0(list.dirs("06_stringtie_gtf", full.names = TRUE, recursive = FALSE), "/t_data.ctab")
# names(files) <- sub(".*\\.", "", dirname(files))

files <- samples[,"file"]
names(files) <- samples[,"sample"]

if (denovo) {
    k <- keys(txdb, keytype = "TXNAME")
    tx2gene <- select(txdb, k, "GENEID", "TXNAME")
} else {
    tmp <- read_tsv(files[1])
    tx2gene <- tmp[, c("t_name", "gene_id")]
    rm(tmp)
}

txi <- tximport(files, type = "stringtie", tx2gene = tx2gene)

## Create DGE object
cts <- txi$counts
normMat <- txi$length

# Obtaining per-observation scaling factors for length, adjusted to avoid
# changing the magnitude of the counts.
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat

# Computing effective library sizes from scaled counts, to account for
# composition biases between samples.
eff_lib <- calcNormFactors(normCts) * colSums(normCts)

# Combining effective library sizes with the length factors, and calculating
# offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff_lib, "*")
normMat <- log(normMat)

# Creating a DGEList object for use in edgeR.
y <- DGEList(cts, group = group)
y <- scaleOffset(y, normMat)

# filtering
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]

# estimate dispersion
y <- estimateDisp(y, design, robust = TRUE)

pdf("07_DE_results/BCV_plot.pdf", width = 6, height = 6)
plotBCV(y)
dev.off()

pdf("07_DE_results/MDS_plot.pdf", width = 6, height = 6)
plotMDS(y)
dev.off()

# fit model
fit <- glmQLFit(y, design, robust = TRUE)

colnames(design)
qlf <- glmQLFTest(fit)

# Smmarize results
topTags(qlf)

summary(decideTests(qlf))

plotMD(qlf)
