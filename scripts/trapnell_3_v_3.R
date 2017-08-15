library(tximport)
srrs <- paste0("SRR",493300 + 66:71)
files <- file.path("quant/salmon",srrs,"quant.sf")
txi <- tximport(files, type="salmon", txOut=TRUE)
kfiles <- file.path("quant/kallisto",srrs,"abundance.tsv") # not bias corrected
#kfiles <- file.path("quant/kallisto_bc",srrs,"abundance.tsv") # bias corrected
ktxi <- tximport(kfiles, type="kallisto", txOut=TRUE)

library(DESeq2)
df <- data.frame(condition = factor(rep(1:2,each=3)))
dds <- DESeqDataSetFromTximport(txi, df, ~ condition)
dds <- estimateSizeFactors(dds)
keep <- rowSums(counts(dds, normalized=TRUE) >= 1) >= 3
table(keep)
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, alpha=.05)
(sal.tot <- sum(res$padj < .05, na.rm=TRUE))

kdds <- DESeqDataSetFromTximport(ktxi, df, ~ condition)
kdds <- estimateSizeFactors(kdds)
keep <- rowSums(counts(kdds, normalized=TRUE) >= 1) >= 3
table(keep)
kdds <- kdds[keep,]
kdds <- DESeq(kdds)
kres <- results(kdds, alpha=.05)
(kal.tot <- sum(kres$padj < .05, na.rm=TRUE))

sig <- rownames(res)[which(res$padj < .05)]
ksig <- rownames(kres)[which(kres$padj < .05)]
both <- length(intersect(ksig, sig))

both
sal.tot - both
kal.tot - both
