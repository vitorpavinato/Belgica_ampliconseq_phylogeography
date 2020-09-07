####################################
#    Mapped Read Depth and IQR     #
####################################

## Vitor Pavinato
## vitor.pavinato@supagro.fr
## CFAES - OSU

rm(list=ls())
ls()

sessionInfo()$running;
sessionInfo()$platform;
R.version.string;
.Platform$GUI;

## Configuration 
genomesize = 89583723 # in bp
input_file = "alignment_coverage.txt"
input_folder = "INDIR"

## Load libraries
#library(ggplot2)

## read samtools output with read depth per nucleotide position
rd <- read.table(file = paste0(input_folder,"/", input_file), header = F)
#rd <- read.table(file = "sample.txt", header = F)
colnames(rd) <- c("CHROM", "POS", paste0("L",1:dim(rd[,-c(1:2)])[2]))

## Mean mapped read depth for every sample
mmrd <- colSums(rd[,-c(1:2)])/genomesize

write.table(mmrd, file = "mean_mapped_read_depth_sample.txt",
            quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")

## Mean mapped read depth across samples
mmrd_mean <- sum(rowSums(rd[,-c(1:2)]))/genomesize/dim(rd[,-c(1:2)])[2]
write.table(mmrd_mean, file = "mean_mapped_read_depth.txt",
            quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")

## Inter-Quantile Range (IQR)

# Mean across samples
mean_quantiles <- quantile(rowMeans(rd[,-c(1:2)]), c(0.25, 0.75))
mean_iqr <- mean_quantiles[2] - mean_quantiles[1]

write.table(mean_iqr, file = "mean_IQR.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# For every sample
quantiles <- apply(rd[,-c(1:2)], 2, quantile, probs=c(0.25, 0.75))

iqr <- quantiles[2, ] - quantiles[1, ]

write.table(iqr, file = "IQR.txt",
            quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")

## Mapped read depth histogram

# Mean across samples
pdf(file = "mean_mapped_read_depth.pdf")
par(mar=c(5,5,4,1)+.1)
hist(rowMeans(rd[,-c(1:2)]), breaks = 100,
     main = "",
     ylab = "Number of Reference Bases",
     xlab = "Mapped Read Depth",
     col = "#035aa6")
dev.off()

# For every sample
pdf(file = "mapped_read_depth.pdf")
par(mar=c(5,5,4,1)+.1, mfrow=c(4, 6))
hist(rd[, paste0("L", 1)], 
     main = paste0("L", 1),
     ylab = "Number of Reference Bases",
     xlab = "Mapped Read Depth",
     col = "#035aa6")

for (i in 2:(dim(rd[,-c(1:2)])[2]-1)){
  hist(rd[, paste0("L", i)], 
       main = paste0("L", i),
       ylab = "",
       xlab = "",
       col = "#035aa6")
}

hist(rd[, paste0("L", dim(rd[,-c(1:2)])[2])], 
     main = paste0("L", dim(rd[,-c(1:2)])[2]),
     ylab = "Number of Reference Bases",
     xlab = "Mapped Read Depth",
     col = "#035aa6")
dev.off()


