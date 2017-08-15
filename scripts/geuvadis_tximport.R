library(tximport)
library(readr)
library(tidyverse)

d <- read.table('/mnt/scratch6/hirak/meta.tsv',header=TRUE)
bias_files_salmon <- file.path("/mnt/scratch6/hirak/sra_archive/GEUV-1/", d$sample, "bias","quant.sf")
names(bias_files_salmon) <- d$sample
bias_files_kallisto <- file.path("/mnt/scratch6/hirak/sra_archive/GEUV-1_kallisto/", d$sample, "bias","abundance.h5")
names(bias_files_kallisto) <- d$sample
txi_kallisto <- tximport(bias_files_kallisto, type="kallisto", txOut=TRUE, countsFromAbundance=c("lengthScaledTPM"), dropInfReps=TRUE)
txi_salmon <- tximport(bias_files_salmon, type="salmon", txOut=TRUE, countsFromAbundance=c("lengthScaledTPM"), dropInfReps=TRUE)
# get rid of gencode names
kallisto_rn <- rownames(txi_kallisto$counts)
m1 <- do.call(rbind, strsplit(kallisto_rn, "|", fixed=TRUE))[,1]
rownames(txi_kallisto$counts) <- m1
# Check if salmon and kallisto row names are equal
all.equal(rownames(salmon),rownames(kallisto))
# Check if salmon and kallisto column names are equal
all.equal(colnames(salmon),colnames(kallisto))

salmon <- txi_salmon$counts
kallisto <- txi_kallisto$counts

# Check if salmon and kallisto row names are equal
all.equal(rownames(salmon),rownames(kallisto))
# Check if salmon and kallisto column names are equal
all.equal(colnames(salmon),colnames(kallisto))

i <- which(colnames(salmon)=="ERR188140")

dat <- tibble(correlation = sapply(1:ncol(salmon), function(i)
cor(log2(kallisto[,i]+1), log2(salmon[,i]+1))),
lab = as.character(d$labs))

my.plot <- dat %>% mutate(lab = reorder(lab, correlation, median)) %>%
  ggplot(aes(lab, correlation, col=lab)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(width = .1, height = 0, pch=21, cex = 0.5) +
  ylab("Correlation between Salmon and kallisto") +
  ggtitle("GEUVADIS Dataset\nPlot based on log(COUNT+1)") +
  annotate("text", x = as.numeric(dat$lab[i]), y = dat$correlation[i],
label = "ERR188140") +
  theme_bw()

ggsave("geu.pdf",my.plot,width=8,height=5)