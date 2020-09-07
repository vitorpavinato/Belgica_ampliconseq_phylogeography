###############################################
#      Detect selection:  Bayescan plots      #
###############################################

# Script to plot Bayescan results
# First you need to run genome scan with Bayescan

## Vitor Pavinato
## vitor.pavinato@supagro.fr
## INRA

rm(list=ls())

sessionInfo()$running;
sessionInfo()$platform;
R.version.string;
.Platform$GUI;

#install.packages(c("boa", "coda"))
library(boa)
library(coda)

recordSessionInfo <- sessionInfo();
save(recordSessionInfo, file="results/detect_selection/bayescan/sessionInfo.RData");
save.image(file = "results/detect_selection/bayescan/workspaceFile.RData");

load(file = "results/detect_selection/bayescan/workspaceFile.RData")

## Posterior distribution of FST and FIS
par.sel<-read.table("results/detect_selection/bayescan/bayescan_1245_snps_short.sel", header = T, colClasses="numeric")
rownames(par.sel) <- par.sel$thin
par.sel <- par.sel[, -c(1)]
head(par.sel)

par.sel.mcmc <- as.mcmc(par.sel)

pdf(file = "results/detect_selection/bayescan/bayescan_f_posterior_convergencey_1245_snps.pdf")
plot(par.sel.mcmc)
dev.off()

# Calculate 95% HPDI -  similar to CI
fisHPD <- boa.hpd(par.sel$Fst1, 0.05)
fstHPD <- boa.hpd(par.sel$Fst2, 0.05)

# Mean and Median for FSI and FST
fisMean <- mean(par.sel$Fst1) 
fisMedian <- median(par.sel$Fst1)

fstMean <- mean(par.sel$Fst2) 
fstMedian <- median(par.sel$Fst2)

pdf(file = "results/detect_selection/bayescan/bayescan_f_posterior_distribution_1245_snps.pdf")
par(mfrow=c(2,1))
# FIS posterior distribution
plot(density(par.sel$Fst1),xlab=expression("F"[IS]),main="Parameter posterior distribution")
abline(v = fisMedian, col = "red")
abline(v = c(fisHPD[1], fisHPD[2]), col = "grey", lty = 2) 

# FST posterior distribution
plot(density(par.sel$Fst2),xlab=expression("F"[ST]), main="")
abline(v = fstMedian, col = "red")
abline(v = c(fstHPD[1], fstHPD[2]), col = "grey", lty = 2)
dev.off()

# Table with posterior parameters estimates
bayesF <- data.frame(fis_mean=fisMean, fis_median=fisMedian, fis_HPDl=fisHPD[1], fis_HPDu=fisHPD[2],
                     fst_mean=fstMean, fst_median=fstMedian, fst_HPDl=fstHPD[1], fst_HPDu=fstHPD[2])

rownames(bayesF) <- ""

write.table(bayesF, file="results/detect_selection/bayescan/bayescan_f_posterior_estimates_1245_snps.txt", sep="\t", quote = F)


# Load this R source to plot bayescan result - it comes with bayescan version 2.1
source("~/Softwares/PopGen/BayeScan2.1/R functions/plot_R.r")

## Outlier detection
# Plot bayescan output with their function
plot_bayescan("results/detect_selection/bayescan/bayescan_fst_1245_snps.txt", FDR=0.10)

# Plot FST outliers with manual R plot
fst.outlier<-read.table("results/detect_selection/bayescan/bayescan_fst_1245_snps.txt", header=TRUE)
head(fst.outlier)
dim(fst.outlier)

snp.colors <- ifelse(fst.outlier$fst > fstHPD[2], "#d53e4f", "black")

pdf(file = "results/detect_selection/bayescan/bayescan_fstQvalue_outlier_1245_snps.pdf", width = 9, height = 5)
par(mar = c(5, 4, 4, 4) + 0.3)
plot(log10(fst.outlier$qval), fst.outlier$fst, 
     xlab=expression(log[10]("q-value")), ylab=expression(paste("Posterior F"[ST])), 
     col=adjustcolor(snp.colors, alpha.f = 0.8), pch=19, cex=1.5, 
     xlim = c(-0.05,-0.4))
#identify(fst.outlier$log10.PO., fst.outlier$fst) 
## identified loci:
## 1260 snps = 611 and 492
## 1245 snps = 596 and 477 (same locus since 596 + 15 that were removed = 611)
abline(h = c(fstHPD[2], fstMedian), col = "#969696", lty = c(2,1))
legend("topleft", legend = c(expression(paste("Upper bound ", F[ST] ," 95% HDP interval")), 
                             expression(paste("Posterior mean ", F[ST]))
                             ), 
       col = c("#969696", "#969696"),  lty = c(2,1), pch = "", bty = "n", cex = 1)
dev.off()

pdf(file = "results/detect_selection/bayescan/bayescan_fstSNPs_outlier_1245_snps.pdf", width = 9, height = 5)
par(mar = c(5, 4, 4, 4) + 0.3)
plot(1:length(rownames(fst.outlier)), fst.outlier$fst, 
     xlab = "SNPs", ylab = expression(paste("Posterior F"[ST])), 
     col=adjustcolor(snp.colors, alpha.f = 0.8), pch=19, cex=1.5)
#identify(1:length(rownames(fst.outlier)), fst.outlier$fst) 
## identified loci:
# 475  477  596  597  673  678 1090 1198
abline(h = c(fstHPD[2], fstMedian), col = "#969696", lty = c(2,1))
par(new = TRUE)
plot(1:length(rownames(fst.outlier)), fst.outlier$prob, 
     pch="", axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side=4, at = pretty(range(fst.outlier$prob)))
mtext("Posterior Probability", side=4, line=3)
legend("topright", legend = c(expression(paste("Upper bound ", F[ST] ," 95% HDP interval")), 
                             expression(paste("Posterior mean ", F[ST]))
                             ), 
col = c("#969696", "#969696"),  lty = c(2,1), pch = "", bty = "n", cex = 1)
dev.off()

# outliers identified:
# 475(490) - gi|676158749|gb|JPYR01001350.1| 7448 7449
# 477(492) - gi|676158749|gb|JPYR01001350.1|	7517 7518
# 596(611) - gi|676158379|gb|JPYR01001535.1|	415	416
# 597(612) - gi|676158379|gb|JPYR01001535.1| 1490 1491
# 673(688) - gi|676157863|gb|JPYR01001792.1| 2457 2458
# 678(693) - gi|676157863|gb|JPYR01001792.1| 2712 2713 
# 1090(1105) - gi|676155571|gb|JPYR01002934.1| 359 360
# 1198(1213) -  gi|676153183|gb|JPYR01004124.1| 114 115

# load OutFlank detected outliers
outliers.outflank <- read.table(file="results/detect_selection/outflank/outlierFinalTable.txt",header = TRUE )

fst.outlier[outliers.outflank$rownames,]

# create the vector of color for OutFlank outliers loci
snp.colors.outflank <- rep("black", dim(fst.outlier)[1])

snp.colors.outflank[outliers.outflank$rownames] <- "#31a354"

# superimpose the Bayescan outliers
snp.colors.outflank[fst.outlier$fst > fstHPD[2]] <- "#d53e4f"

pdf(file = "results/detect_selection/bayescan/bayescanOutFlank_fstSNPs_outlier_1245_snps.pdf", width = 9, height = 5)
par(mar = c(4, 4, 2, 4) + 0.3)
plot(1:length(rownames(fst.outlier)), fst.outlier$fst, 
     xlab = "SNPs", ylab = expression(paste("Posterior F"[ST])), 
     col=adjustcolor(snp.colors.outflank, alpha.f = 0.8), pch=19, cex=1.5)
abline(h = c(fstHPD[2], fstMedian), col = "#969696", lty = c(2,1))
par(new = TRUE)
plot(1:length(rownames(fst.outlier)), fst.outlier$prob, 
     pch="", axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side=4, at = pretty(range(fst.outlier$prob)))
mtext("Posterior Probability", side=4, line=3)
legend("topleft", 
       legend = c(expression(paste("Upper bound ", F[ST] ," 95% HDP interval")), 
                  expression(paste("Posterior mean ", F[ST])),
                  expression(paste("Highest posterior ", F[ST], "'s")),
                  "OutFlank outliers"
                  ), 
       col = c(rep("#969696", 2),
               adjustcolor("#d53e4f", alpha.f = 0.8),
               adjustcolor("#31a354", alpha.f = 0.8)
               ),  
       lty = c(2,1,0,0), pch = c(26,26,19,19), bty = "n", cex = 1)
dev.off()


save.image(file = "results/detect_selection/bayescan/workspaceFile.RData");

