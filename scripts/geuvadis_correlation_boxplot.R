## The code for computing correlations 
## and code for the plot are fine and can be used directly

### read in data quickly with scan
library(matrixStats)
kallisto <- scan("GEUV_all/GEUVADIS_bias_corrected_kallisto.tsv", what="c")
salmon <- scan("GEUV_all/GEUVADIS_bias_corrected_salmon.tsv", what="c")
pd <- read.table("GEUV_all/meta.tsv", header=TRUE, stringsAsFactors = FALSE)

## convert into numeric matrices
## with log2( values + 1)

##first kallisto
n <- length(kallisto)
ncol <- min(which(!grepl("ERR", kallisto[2:1000]))) 
nrow <- n/ncol
rownames <- kallisto[seq(ncol+1, n, ncol)]
colnames <- kallisto[2:ncol]

##sanity check:
stop(!identical(colnames,pd[,1]))

kallisto  <- kallisto[-(1:ncol)]
kallisto <- kallisto[-seq(1, n-ncol, ncol)]
kallisto <- matrix(log2(as.numeric(kallisto)+1), ncol=ncol-1, nrow=nrow-1,
            byrow = TRUE)

## same for Salmon
n <- length(salmon)
ncol <- min(which(!grepl("ERR", salmon[2:1000]))) 
nrow <- n/ncol
rownames_s <- salmon[seq(ncol+1, n, ncol)]

##sanity check:
stop(!identical(rownames,rownames_s))

##sanity check:
colnames <- salmon[2:ncol]
stop(!identical(colnames,pd[,1]))

salmon  <- salmon[-(1:ncol)]
salmon <- salmon[-seq(1, n-ncol, ncol)]
salmon <- matrix(log2(as.numeric(salmon)+1), ncol=ncol-1, nrow=nrow-1,
                   byrow = TRUE)

### Compute correlations 
library(tidyverse)
dat <- tibble(
  correlation = sapply(1:ncol(salmon), function(i) cor(kallisto[,i], salmon[,i])),
  lab = as.character(pd$labs))

### which one is ERR188140
i <- which(colnames=="ERR188140")

### make the plot
dat %>% mutate(lab = reorder(lab, correlation, median)) %>%
  ggplot(aes(lab, correlation, col=lab)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(width = .1, height = 0, pch=21, cex = 0.5) +
  ylab("Correlation between Salmon and kallisto") +
  ggtitle("GEUVADIS Dataset\nPlot based on log(TPM+1)") + 
  annotate("text", x = as.numeric(pd$lab[i]), y = dat$correlation[i], label = "ERR188140") +
  theme_bw()
  
