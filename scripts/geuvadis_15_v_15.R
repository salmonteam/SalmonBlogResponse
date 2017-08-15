library(tximport)
dat <- read.delim("geuvadis_samples.txt", header=TRUE)
condition <- dat$Performer
run <- as.character(dat$"Comment.ENA_RUN.")

files <- file.path("geuvadis_quants/salmon",run,"quant.sf")
txi <- tximport(files, type="salmon", txOut=TRUE)

kfiles <- file.path("geuvadis_quants/kallisto",run,"abundance.tsv")
ktxi <- tximport(kfiles, type="kallisto", txOut=TRUE)

library(genefilter)
logtpm <- log2(txi$abundance + 1)
logtpm <- logtpm[rowMeans(txi$abundance) > .1,]
tt <- rowttests(logtpm, condition)
sum(p.adjust(tt$p.value, method="BH") < .01, na.rm=TRUE)

klogtpm <- log2(ktxi$abundance + 1)
klogtpm <- klogtpm[rowMeans(ktxi$abundance) > .1,]
ktt <- rowttests(klogtpm, condition)
sum(p.adjust(ktt$p.value, method="BH") < .01, na.rm=TRUE)

library(DESeq2)
df <- data.frame(condition)
dds <- DESeqDataSetFromTximport(txi, df, ~ condition)
dds <- estimateSizeFactors(dds)
keep <- rowSums(counts(dds, normalized=TRUE) >= 10) >= 5
table(keep)
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, alpha=.01)
(sal.tot <- sum(res$padj < .01, na.rm=TRUE))

kdds <- DESeqDataSetFromTximport(ktxi, df, ~ condition)
kdds <- estimateSizeFactors(kdds)
keep <- rowSums(counts(kdds, normalized=TRUE) >= 10) >= 5
table(keep)
kdds <- kdds[keep,]
kdds <- DESeq(kdds)
kres <- results(kdds, alpha=.01)
(kal.tot <- sum(kres$padj < .01, na.rm=TRUE))

txi.cfa <- tximport(files, type="salmon", txOut=TRUE, countsFromAbundance="lengthScaledTPM")
ktxi.cfa <- tximport(kfiles, type="kallisto", txOut=TRUE, countsFromAbundance="lengthScaledTPM")

limma <- function(x, design) {
  dgel <- DGEList(x$counts)
  dgel <- calcNormFactors(dgel)
  L <- min(colSums(dgel$counts)/1e6)
  keep <- rowSums(cpm(dgel) > 10/L) >= 5
  dgel <- dgel[keep,]
  v <- voom(dgel,design,plot=FALSE)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  topTable(fit, number=nrow(dgel), sort.by="none")
}

library(limma)
library(edgeR)

design <- model.matrix(~condition)
limres <- limma(txi.cfa, design)
sum(limres$adj.P.Val < .01)

klimres <- limma(ktxi.cfa, design)
sum(klimres$adj.P.Val < .01)


