files <- c("just_quants/A_BGI_1/quant_bc/quant.sf",
           "just_quants/A_CNL_1/quant_bc/quant.sf")
kfiles <- c("just_quants/A_BGI_1/kquant_bc/abundance.tsv",
            "just_quants/A_CNL_1/kquant_bc/abundance.tsv")
library(tximport)
txi <- tximport(files, type="salmon", txOut=TRUE, countsFromAbundance="lengthScaledTPM")
ktxi <- tximport(kfiles, type="kallisto", txOut=TRUE, countsFromAbundance="lengthScaledTPM")

x <- txi$abundance[,1]
y <- txi$abundance[,2]

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
doleg <- function() legend("topright",c("Salmon","kallisto"),lwd=3,col=1:2,cex=2,inset=.05)
library(RColorBrewer)

bigpar(2,3)
palette(brewer.pal(4,"Set1"))
maplot(txi$abundance[,1],txi$abundance[,2],col=1,main="Salmon A-CNL-1 vs A-BGI-1")
maplot(ktxi$abundance[,1],ktxi$abundance[,2],col=2,main="kallisto A-CNL-1 vs A-BGI-1")
maplot(txi$abundance[,1],txi$abundance[,2],col=1,pts=FALSE,ylim=c(-3,3),main="5%, 95% quant. A-CNL-1 vs A-BGI-1")
maplot(ktxi$abundance[,1],ktxi$abundance[,2],col=2,add=TRUE,pts=FALSE)
doleg()
maplot(txi$counts[,1],txi$counts[,2],col=1,main="Salmon A-CNL-1 vs A-BGI-1",xlab="log10 counts")
maplot(ktxi$counts[,1],ktxi$counts[,2],col=2,main="kallisto A-CNL-1 vs A-BGI-1",xlab="log10 counts")
maplot(txi$counts[,1],txi$counts[,2],col=1,pts=FALSE,ylim=c(-3,3),
       main="5%, 95% quant. A-CNL-1 vs A-BGI-1",xlab="log10 counts")
maplot(ktxi$counts[,1],ktxi$counts[,2],col=2,add=TRUE,pts=FALSE,xlab="log10 counts")
doleg()


bigpar(2,2)
palette(brewer.pal(4,"Set1"))
maplot(txi$abundance[,1],ktxi$abundance[,1],col=3,main="kallisto vs. Salmon A-BGI-1")
maplot(txi$abundance[,2],ktxi$abundance[,2],col=3,main="kallisto vs. Salmon A-CNL-1")
maplot(txi$counts[,1],ktxi$counts[,1],col=3,main="kallisto vs. Salmon A-BGI-1",
       xlab="log10 counts")
maplot(txi$counts[,2],ktxi$counts[,2],col=3,main="kallisto vs. Salmon A-CNL-1",
       xlab="log10 counts")

