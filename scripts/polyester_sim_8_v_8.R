load("simulate.rda")
library(tximport)
library(readr)
library(Biostrings)

fasta <- "reference.fa"
txseq <- readDNAStringSet(fasta)
# only those transcripts that were simulated
load(paste0("out/out_1/sim_counts_matrix.rda"))
sub.cts <- matrix(NA,nrow=nrow(counts_matrix),ncol=16)
for (i in 1:8) {
  load(paste0("out/out_",i,"/sim_counts_matrix.rda"))
  sub.cts[,c(i,i+8)] <- counts_matrix
}
rownames(sub.cts) <- rownames(counts_matrix)
# all transcripts in the FASTA
cts <- matrix(0,nrow=length(txseq),ncol=16)
rownames(cts) <- names(txseq)
cts[match(rownames(sub.cts),rownames(cts)),] <- sub.cts
# now calculate gold standard TPMs
gold <- cts
frag.length <- 200 # needs to match with simulate_reads.R
for (i in 1:16) {
  gold[,i] <- cts[,i] / pmax(width(txseq) - frag.length + 1, 1)
}
gold <- sweep(gold, 2, 1e6/colSums(gold), `*`)

n <- 8 # samples per group
dirs <- paste0("out/out_",rep(1:n,2),"/sample_0",rep(1:2,each=n))
tpm <- list()
txi <- list()
txi.cfa <- list()
meths <- c("Salmon", "kallisto")
dir.type <- c(Salmon="salmon", kallisto="kallisto")
file.type <- c(Salmon="quant.sf", kallisto="abundance.tsv")
type <- c(Salmon="salmon", kallisto="kallisto")
for (m in meths) {
  txi[[m]] <- tximport(file.path(dirs,dir.type[m],file.type[m]),
                       type=type[m],
                       txOut=TRUE, dropInfReps=TRUE)
  txi.cfa[[m]] <- tximport(file.path(dirs,dir.type[m],file.type[m]),
                           type=type[m],
                           txOut=TRUE, dropInfReps=TRUE,
                           countsFromAbundance="lengthScaledTPM")
  tpm[[m]] <- txi[[m]]$abundance
}
tpm0 <- tpm

tpm <- tpm0
idx <- rowSums(gold) > 10
gold.mid1 <- median(apply(gold[idx,1:n], 2, median))
gold.mid2 <- median(apply(gold[idx,(n+1):(2*n)], 2, median))
for (m in meths) {
  mid1 <- median(apply(tpm[[m]][idx,1:n], 2, median))
  mid2 <- median(apply(tpm[[m]][idx,(n+1):(2*n)], 2, median))
  tpm[[m]][,1:n] <- tpm[[m]][,1:n] / mid1 * gold.mid1
  tpm[[m]][,(n+1):(2*n)] <- tpm[[m]][,(n+1):(2*n)] / mid2 * gold.mid2
}

errfun_mae <- function(x,y,pc=1) {
  idx <- x > pc
  x2 <- log2(x[idx] + pc)
  y2 <- log2(y[idx] + pc)
  median( abs(x2 - y2), na.rm=TRUE )
}
myScatter <- function(x,y,algo,pc=1,xmax=3,xlab="log10(TPM)",ylab="") {
  err <- as.data.frame(log2(y+pc)-log2(x+pc))
  colnames(err) <- 'val'
  err$col[abs(err$val) > 0.5] <- rgb(0.8, 0.0, 0.0, 0.1)
  err$col[abs(err$val) <= 0.5] <- rgb(0, 0, 0, 0.1)

  out <-  sum(abs(err$val[x>0]) > 0.5)
  fout <- out / sum(x > 0)
  
  plot(log10(x+pc),log2(y+pc)-log2(x+pc),
       xlim=c(log10(pc),5),ylim=c(-2,2),
       main=algo,xlab=xlab,ylab=ylab,cex=.6,pch=16,col=err$col)
  abline(h=0,col="red",lwd=1)
  usr <- par( "usr" )
  efrac <- sprintf("%.3f", fout)
  legend("topright", legend=efrac, text.col="red", bty="n")
}

par(mfrow=c(2,2), mar=c(5,2,3,0))
for (i in c(1,9)) {
  myScatter(gold[,i], tpm[["Salmon"]][,i], "Salmon")
  myScatter(gold[,i], tpm[["kallisto"]][,i], "kallisto")
}

err <- lapply(0:1, function(j) {
         lapply(meths, function(m) {
           sapply(1:n, function(i) {
             errfun_mae(gold[,i + j*n], tpm[[m]][,i + j*n])
                })
              })
            })
err <- do.call(c, err)
head <- paste(rep(meths, 2), rep(c("low", "high"), each=4))
tdat <- as.data.frame(unlist(lapply(err, mean)), row.names=head)
colnames(tdat) <- "Mean MAE"
kable(t(tdat), format="html", digits=4, table.attr='class="flat-table"')

library(RColorBrewer)
palette(brewer.pal(4,"Set1"))
par(mfrow=c(1,1),mar=c(7,5,1,1))
cols <- c(seq_along(meths),seq_along(meths))
par(mar=c(9,8,1,1))   # extra large bottom margin
boxplot(err, ylab="median( | log2(estimate/truth) | )",
        las=3, range=0, border=cols, names=rep(meths,2), ylim=c(0, 0.25), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
for (i in 1:(2*length(meths))) points(runif(n,i-.1,i+.1), err[[i]], col=cols[i])
abline(v = length(meths) + .5, col="gray")

gold.de <- numeric(length(txseq))
names(gold.de) <- names(txseq)
gold.de[rownames(fold_changes)] <- log2(fold_changes[,2])
condition <- factor(rep(1:2, each=n))
design <- model.matrix(~ condition)
library(limma)
# this function uses limma's lmFit for speed
# it performs simple t-tests, not moderated t-statistics / ebayes methods
ttests <- function(x, design, pc=1, coef=ncol(design)) {
  # correct for global shift across condition due to different library size
  sf <- exp(median(log(rowSums(x[,1:n]) / rowSums(x[,(n+1):(2*n)])), na.rm=TRUE))
  x[,(n+1):(2*n)] <- sf * x[,(n+1):(2*n)]
  fit <- lmFit(log2(x + pc), design)
  ordinary.t <- fit$coef / fit$stdev.unscaled / fit$sigma
  pvals <- 2*pt(abs(ordinary.t[,coef]), df=ncol(x)-ncol(design), lower.tail=FALSE)
  data.frame(dm=fit$coef[,coef], pvalues=pvals)
}


library(DESeq2)
deseq2tests <- function(x, condition) {
  dds <- DESeqDataSetFromTximport(x, data.frame(condition), ~condition)
  dds <- estimateSizeFactors(dds)
  keep <- rowSums(counts(dds, normalized=TRUE) >= 10) >= 3
  dds <- dds[keep,]
  dds <- DESeq(dds, fitType="local", minRep=Inf)
  res <- results(dds)
  df <- data.frame(dm=numeric(nrow(x$counts)), pvalues=rep(1,nrow(x$counts)))
  df$dm[keep] <- res$log2FoldChange
  df$pvalues[keep] <- res$pvalue
  df
}

library(edgeR)
library(limma)
limmatests <- function(x, design) {
  dgel <- DGEList(x$counts)
  dgel <- calcNormFactors(dgel)
  L <- min(colSums(dgel$counts)/1e6)
  keep <- rowSums(cpm(dgel) > 10/L) >= 3
  dgel <- dgel[keep,]
  v <- voom(dgel,design,plot=FALSE)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  tt <- topTable(fit, number=nrow(dgel), sort.by="none")
  df <- data.frame(dm=numeric(nrow(x$counts)), pvalues=rep(1,nrow(x$counts)))
  df$dm[keep] <- tt$logFC
  df$pvalues[keep] <- tt$P.Value
  df
}

tlist <- list()
llist <- list()
dlist <- list()
tlist[["gold"]] <- ttests(gold, design)
for (m in meths) {
  tlist[[m]] <- ttests(tpm[[m]], design)
  llist[[m]] <- limmatests(txi.cfa[[m]], design)
  dlist[[m]] <- deseq2tests(txi[[m]], condition)
}

par(mfrow=c(1,3), mar=c(5,5,1,1))
for (m in c("gold",meths)) {
  boxplot(tlist[[m]]$dm ~ gold.de, main=m)
  abline(h=-1:1, col="red")
}

sensPrecCurve <- function(tlist, gold.de, xlim=c(0,1), ylim=c(0,1), add=FALSE, lty=1, ...) {
  for (i in seq_along(tlist)) {
    ps <- 10^seq(from=-100, to=0, length=500)
    tpr <- sapply(ps, function(p) mean(tlist[[i]]$pvalues[gold.de != 0] < p, na.rm=TRUE))
    fdr <- sapply(ps, function(p) mean(gold.de[tlist[[i]]$pvalues < p] == 0, na.rm=TRUE))
    if (i == 1 & !add) {
      plot(fdr, tpr, type="l", xlim=xlim, ylim=ylim,
           cex=1.5, pch=20, col=i, lwd=1, lty=lty, ...)
    } else {
      points(fdr, tpr, type="l", cex=1.5, pch=20, col=i, lwd=1, lty=lty)
    }
  }
}

library(rafalib)
bigpar()
palette(brewer.pal(4,"Set1"))
sensPrecCurve(tlist[-1], gold.de,
         xlab="False Discovery Rate",
         ylab="Sensitivity", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
sensPrecCurve(dlist, gold.de, add=TRUE, lty=2)
sensPrecCurve(llist, gold.de, add=TRUE, lty=3)
legend("bottomright", c("Salmon t-test TPM","Salmon DESeq2",
                        "Salmon CPM limma-voom",
                        "kallisto t-test TPM","kallisto DESeq2",
                        "kallisto CPM limma-voom"),
       col=rep(1:2,each=3), lty=c(1:3,1:3), lwd=1, inset=.05, bg="white")

maxSensForPrec <- function(tlist, gold.de, thresholds=c(.01,.05,.1), low, high) {
  m <- matrix(NA, ncol=length(tlist), nrow=length(thresholds))
  for (i in seq_along(thresholds)) {
    for (j in seq_along(tlist)) {
      ps <- rev(10^seq(from=low, to=high, length=500))
      tpr <- sapply(ps, function(p) mean(tlist[[j]]$pvalues[gold.de != 0] < p,
                                         na.rm=TRUE))
      fdr <- sapply(ps, function(p) mean(gold.de[tlist[[j]]$pvalues < p] == 0,
                                         na.rm=TRUE))
      ps.idx <- which(fdr < thresholds[i])[1]
      print(log10(ps[ps.idx]))
      print((thresholds[i] - fdr[ps.idx])/thresholds[i])
      m[i,j] <- tpr[ps.idx]
    }
  }
  rownames(m) <- as.character(thresholds)
  colnames(m) <- names(tlist)
  m
}

maxSensForPrec(tlist[-1], gold.de, low=-15, high=-5)
maxSensForPrec(dlist, gold.de, low=-60, high=-10)
maxSensForPrec(llist, gold.de, low=-12, high=-6)

kable(maxSens, format="html", digits=4, table.attr='class="flat-table"')

