system("tar -xvf GEUV_quants.tgz")
library(tximport)
files <- c("GEUV_quants/ERR188140/salmon_gc/quant.sf",
           "GEUV_quants/ERR188140/kallisto/bias/abundance.tsv")
file.exists(files)
txi <- tximport(files[1], type="salmon", txOut=TRUE,
                countsFromAbundance="lengthScaledTPM")
ktxi <- tximport(files[2], type="kallisto", txOut=TRUE,
                 countsFromAbundance="lengthScaledTPM", dropInfReps=TRUE)

tpm <- cbind(txi$abundance, ktxi$abundance)
counts <- cbind(txi$counts, ktxi$counts)

maplot <- function(x,y,pts=TRUE,add=FALSE,col,ylim=c(-10,10),xlab="log10 TPM",...) {
  mu <- 0.5 * (log10(x+1) + log10(y+1))
  idx <- x > 1 & y > 1
  lfc <- log2(y) - log2(x)
  mu <- mu[idx]
  lfc <- lfc[idx]
  nbrks <- 20
  brks <- c(0,quantile(mu, probs=seq(0,1,length=nbrks-2)),max(mu)+1)
  mu.cut <- cut(mu, brks)
  qs <- sapply(levels(mu.cut), function(i) {
    quantile(lfc[mu.cut == i], c(.05, .95))
  })
  xlim <- c(.25,quantile(mu,.9999))
  if (!add) plot(mu,lfc,xlim=xlim,ylim=ylim,type="n",ylab="log2 fold change",xlab=xlab,...)
  if (pts) points(mu, lfc)
  if (!add) abline(h=0, col="grey")
  lines(.5*(brks[-1]+brks[-nbrks]),qs[2,],col=col,lwd=3)
  lines(.5*(brks[-1]+brks[-nbrks]),qs[1,],col=col,lwd=3)
}

library(rafalib)
library(RColorBrewer)

bigpar(1,2)
palette(brewer.pal(4,"Set1"))
maplot(tpm[,1], tpm[,2], col=3, ylim=c(-5,5))
maplot(counts[,1], counts[,2], col=3, xlab="log10 counts", ylim=c(-5,5))


