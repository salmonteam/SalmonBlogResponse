salmon <- tximport::tximport("GEUV_quants/ERR188140/salmon_gc/quant.sf",
                             type="salmon", txOut=TRUE, countsFromAbundance="lengthScaledTPM")
kallisto <- tximport::tximport("GEUV_quants/ERR188140/kallisto/bias/abundance.h5", 
                               type="kallisto", txOut=TRUE, countsFromAbundance="lengthScaledTPM")


maplot <- function(x, y, ...) 
  plot( log2(x)/2+ log2(y)/2, log2(y)-log2(x), ...)

## the correlation is 
cor(log(kallisto$counts+1), log(salmon$counts+1))

png(file="~/Desktop/ma-plot.png",
    width=800,
    height=425)

par(mfrow=c(1,2),
    cex.lab = 1.35, 
    cex.main = 1.45, 
    cex.axis = 1,
    mar = c(3, 3,  2.6, 0.5),
    mgp = c(1.5, 0.5, 0))

plot(kallisto$counts, salmon$counts,
     ylab = "Salmon_gc counts",
     xlab = "kallisto counts",
     main = "ERR188140: Scatterplot like one in blog post")
abline(0, 1, lty=2, col="grey")

maplot(kallisto$counts+1, salmon$counts+1,
       ylab = "Difference between Salmon and kallisto",
       xlab = "Average of Salmon and kallisto",
       main = "MA-plot version of plot on the left")
abline(h = c(-1,1), col="red", lty=2)
abline(h = 0, lty=2, col="grey")
dev.off()


