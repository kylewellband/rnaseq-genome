setwd("~/Projects/sfo_eqtl/")

library(Rsubread)
library(edgeR)

gff_annot <- "~/Projects/sfo_eqtl/01data/sfon_gene_salp_pos.gff"
bam_dir <- "/mnt/HDD1/sfo_map/RNAseq"

SE_bam_files <- system(paste0("ls ", bam_dir, "/*.bam"), intern = TRUE)[-1:-4]

if (!(file.exists("~/Projects/sfo_eqtl/01data/exon_feature_counts.RDS"))) {
    fc_SE <- featureCounts(SE_bam_files, annot.ext = gff_annot,
                           isGTFAnnotationFile = TRUE, GTF.featureType = "exon", GTF.attrType = "ID",
                           useMetaFeatures = FALSE, allowMultiOverlap = TRUE)
    saveRDS(fc_SE, "~/Projects/sfo_eqtl/01data/exon_feature_counts.RDS")
} else {
    fc_SE <- readRDS("~/Projects/sfo_eqtl/01data/exon_feature_counts.RDS")
}

y.all <- DGEList(counts=fc_SE$counts, genes=fc_SE$annotation)
y <- sumTechReps(y.all, ID = sub(".bam$", "", row.names(y.all$samples)))

dim(y)
head(y$genes)

samples <- read.table("~/Dropbox/sfo_eqtl/01data/additional_file_s7_sfeq_interpretation_v1.csv", header = TRUE, sep = ",")
samples <- samples[grep("SFO", samples$fish.id),]
row.names(samples) <- sub("_", ".", samples$fish.id)
samples <- samples[sub(".bam", "", row.names(y$samples)),]

Sex <- factor(samples$sex)
keep <- filterByExpr(y, group = Sex)
table(keep)

y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples

design <- model.matrix(~Sex)

y <- estimateDisp(y, design = design)
plotBCV(y)

fit <- glmQLFit(y, design)
plotQLDisp(fit)

qlf <- glmQLFTest(fit, coef=2)
topTags(qlf)

is_de <- decideTests(qlf, p.value = 0.05, lfc = log2(1.5))
summary(is_de)

fit$genes$ExonID <- fit$genes$GeneID
fit$genes$GeneID <- sub(".exon[0-9]+", "", fit$genes$GeneID)
sp <- diffSpliceDGE(fit, coef=2, geneid="GeneID", exonid="ExonID")
dim(topSpliceDGE(sp, test="gene", n=100000))
sum(topSpliceDGE(sp, test="gene", n=100000)$FDR < 0.05)

plotSpliceDGE(sp, geneid="QSF_A1CF.2.2.mrna1", genecol="GeneID")
plotSpliceDGE(sp, geneid="QSF_ASAH2.1.1.mrna1", genecol="GeneID")
plotSpliceDGE(sp, geneid="QSF_CFB.3.11.mrna1", genecol="GeneID")
plotSpliceDGE(sp, geneid="QSF_MYLK.3.3.mrna1", genecol="GeneID")
plotSpliceDGE(sp, geneid="QSF_WDR37.3.3.mrna1", genecol="GeneID")
plotSpliceDGE(sp, geneid="QSF_DKC1.3.3.mrna1", genecol="GeneID")
plotSpliceDGE(sp, geneid="QSF_DKC1.2.3.mrna1", genecol="GeneID")
plotSpliceDGE(sp, geneid="QSF_DMXL1.1.1.mrna1", genecol="GeneID")
