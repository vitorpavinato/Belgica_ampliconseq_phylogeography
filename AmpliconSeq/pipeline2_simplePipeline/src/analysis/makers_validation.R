##########################################
##       SITE MEAN DEPTH AND QUAL,      ##
##   PARALOGOUS IDENTIFICATION ISSUES,  ##
##      QUAL and DEPTH RELATIONTHIP,    ##
##     ALLELE FREQUENCY DISTRIBUTION,   ##
##   GENOTYPE FREQUENCY DISTRIBUTION,   ##
##        AND MARKER VALIDATION         ##
##########################################

## Vitor Pavinato
## vitor.pavinato@supagro.fr
## CFAES - OSU

rm(list=ls())
ls()

## LOAD LIBRARIES
library(scales)

sessionInfo()$running;
sessionInfo()$platform;
R.version.string;
.Platform$GUI;

#######################
### SITE MEAN DEPTH ###
#######################

## For each site, DEPTH averaged across samples.
## Results obtained with "vcftools --site-mean-depth".
## This output file has the suffix ".ldepth.mean".

## WITHIN AMPLICON - TARGETED variants
##------------------------

# mean-site depth
targeted.snps.dp <- read.table(file = "results/coverage/calledVariants/withinAmplicons/targeted_snps.ldepth.mean",
                               na.strings = ".", header = TRUE)
dim(targeted.snps.dp)

targeted.snps.dp <- targeted.snps.dp[complete.cases(targeted.snps.dp), ]
dim(targeted.snps.dp)
# 47  4

# Mean across samples
targeted.mean_quantiles <- quantile(targeted.snps.dp$MEAN_DEPTH, c(0.25, 0.75))
targeted.mean_iqr <- targeted.mean_quantiles[2] - targeted.mean_quantiles[1]
# 524.195

# Mean across samples
#pdf(file = "raw_mean_site_depth.pdf")
par(mar=c(5,5,4,1)+.1)
hist(t(targeted.snps.dp$MEAN_DEPTH), breaks = 10, freq = T,
     main = "",
     ylab = "Number of Variant Sites",
     xlab = "Site Depth",
     col = "#035aa6")
#dev.off()

min(targeted.snps.dp$MEAN_DEPTH) # 784.143
max(targeted.snps.dp$MEAN_DEPTH) # 6066.95
mean(targeted.snps.dp$MEAN_DEPTH) # 3152.887
sqrt(var(targeted.snps.dp$MEAN_DEPTH)) # 1017.195
 
mean(targeted.snps.dp$MEAN_DEPTH) - sqrt(var(targeted.snps.dp$MEAN_DEPTH)) # 2135.692
mean(targeted.snps.dp$MEAN_DEPTH) + sqrt(var(targeted.snps.dp$MEAN_DEPTH)) # 4170.082

dens1 <- density(targeted.snps.dp$MEAN_DEPTH, na.rm = T)
plot(dens1)

## WITHIN AMPLICON - NON-TARGETED variants
##------------------------
nontargeted.snps.dp <- read.table(file = "results/coverage/calledVariants/withinAmplicons/nontargeted_snps.ldepth.mean",
                                  na.strings = ".", header = TRUE)
dim(nontargeted.snps.dp)

nontargeted.snps.dp <- nontargeted.snps.dp[complete.cases(nontargeted.snps.dp), ]
dim(nontargeted.snps.dp)
#  126   4

# Mean across samples
nontargeted.mean_quantiles <- quantile(nontargeted.snps.dp$MEAN_DEPTH, c(0.25, 0.75))
nontargeted.mean_iqr <- nontargeted.mean_quantiles[2] - nontargeted.mean_quantiles[1]
# 1986.3

# Mean across samples
#pdf(file = "raw_mean_site_depth.pdf")
par(mar=c(5,5,4,1)+.1)
hist(t(nontargeted.snps.dp$MEAN_DEPTH), breaks = 10, freq = T,
     main = "",
     ylab = "Number of Variant Sites",
     xlab = "Site Depth",
     col = "#035aa6")
#dev.off()

min(nontargeted.snps.dp$MEAN_DEPTH) # 709.25
max(nontargeted.snps.dp$MEAN_DEPTH) # 3973.95
mean(nontargeted.snps.dp$MEAN_DEPTH) # 2621.971
sqrt(var(nontargeted.snps.dp$MEAN_DEPTH)) # 1008.11

mean(nontargeted.snps.dp$MEAN_DEPTH) - sqrt(var(nontargeted.snps.dp$MEAN_DEPTH)) # 1613.861
mean(nontargeted.snps.dp$MEAN_DEPTH) + sqrt(var(nontargeted.snps.dp$MEAN_DEPTH)) # 3630.081

dens2 <- density(nontargeted.snps.dp$MEAN_DEPTH)
plot(dens2)

## OUTSIDE AMPLICON - OUTSIDE PRIMER REGION
##------------------------
outsidePrimer.snps.dp <- read.table(file = "results/coverage/calledVariants/outsideAmplicon/variants_outside_outsidePrimers.ldepth.mean",
                                    na.strings = ".", header = TRUE)

dim(outsidePrimer.snps.dp)
# 166   4

# List of identified SNPs that were probably produced by a mutation (insertion)
# that moved a bit the primer matching
slipped.primersSNPs <- c(2,3,4,15,16,29,30,31,32,33,35,36,37,41,42,64,65,76,77,117)
# 20

toBeAddedPrimersList <- outsidePrimer.snps.dp[slipped.primersSNPs, ]

# remove these from the table of nontargeted SNPs outside primer regions
# add these to the list of nontargeted SNPs within primer regions
outsidePrimer.snps.dp <- outsidePrimer.snps.dp[-slipped.primersSNPs,]
dim(outsidePrimer.snps.dp)
# 146   4

outsidePrimer.snps.dp <- outsidePrimer.snps.dp[complete.cases(outsidePrimer.snps.dp), ]
dim(outsidePrimer.snps.dp)
# 118   4

# Mean across samples
#pdf(file = "raw_mean_site_depth.pdf")
par(mar=c(5,5,4,1)+.1)
hist(t(outsidePrimer.snps.dp$MEAN_DEPTH), breaks = 10, freq = T,
     main = "",
     ylab = "Number of Variant Sites",
     xlab = "Site Depth",
     col = "#035aa6")
#dev.off()

min(outsidePrimer.snps.dp$MEAN_DEPTH) # 1.5
max(outsidePrimer.snps.dp$MEAN_DEPTH) # 380.524
mean(outsidePrimer.snps.dp$MEAN_DEPTH) # 78.74609
sqrt(var(outsidePrimer.snps.dp$MEAN_DEPTH)) # 97.62463

dens3 <- density(outsidePrimer.snps.dp$MEAN_DEPTH)
plot(dens3)

## OUTSIDE AMPLICON - WITHIN PRIMER REGION
##------------------------
withinPrimer.snps.dp <- read.table(file = "results/coverage/calledVariants/outsideAmplicon/variants_outside_withinPrimers.ldepth.mean",
                                   na.strings = ".", header = TRUE)

dim(withinPrimer.snps.dp)
# 54  4

withinPrimer.snps.dp <- rbind(withinPrimer.snps.dp, toBeAddedPrimersList)
dim(withinPrimer.snps.dp)
# 74  4

withinPrimer.snps.dp <- withinPrimer.snps.dp[complete.cases(withinPrimer.snps.dp), ]
dim(withinPrimer.snps.dp)
# 63  4

# Mean across samples
#pdf(file = "raw_mean_site_depth.pdf")
par(mar=c(5,5,4,1)+.1)
hist(t(withinPrimer.snps.dp$MEAN_DEPTH), breaks = 10, freq = T,
     main = "",
     ylab = "Number of Variant Sites",
     xlab = "Site Depth",
     col = "#035aa6")
#dev.off()

min(withinPrimer.snps.dp$MEAN_DEPTH) # 1.11111
max(withinPrimer.snps.dp$MEAN_DEPTH) # 7501
mean(withinPrimer.snps.dp$MEAN_DEPTH) # 933.3676
sqrt(var(withinPrimer.snps.dp$MEAN_DEPTH)) # 1741.158

dens4 <- density(withinPrimer.snps.dp$MEAN_DEPTH)
plot(dens4)

## Combinade coverage histograms
##------------------------------
c1 <- rgb(49,54,149, max = 255, alpha = 100, names = "dk.blue")
c2 <- rgb(171,217,233, max = 255, alpha = 90, names = "lt.blue")
c3 <- rgb(253,174,97, max = 255, alpha = 90, names = "lt.orange")
c4 <- rgb(165,0,38, max = 255, alpha = 90, names = "dk.red")

b <- min(c(targeted.snps.dp$MEAN_DEPTH, nontargeted.snps.dp$MEAN_DEPTH,
           withinPrimer.snps.dp$MEAN_DEPTH, outsidePrimer.snps.dp$MEAN_DEPTH)) # Set the minimum for the breakpoints
e <- max(c(targeted.snps.dp$MEAN_DEPTH, nontargeted.snps.dp$MEAN_DEPTH,
            withinPrimer.snps.dp$MEAN_DEPTH, outsidePrimer.snps.dp$MEAN_DEPTH)) # Set the maximum for the breakpoints
ax <- pretty(b:e, n = 25) # Make a neat vector for the breakpoints
ax

hg.targeted <- hist(targeted.snps.dp$MEAN_DEPTH, breaks = ax, plot = FALSE) 
hg.nontargeted <- hist(nontargeted.snps.dp$MEAN_DEPTH, breaks = ax, plot = FALSE) 
hg.wPrimers <- hist(withinPrimer.snps.dp$MEAN_DEPTH, breaks = ax, plot = FALSE) 
hg.oPrimers <- hist(outsidePrimer.snps.dp$MEAN_DEPTH, breaks = ax, plot = FALSE) 

pdf(file="results/coverage/calledVariants/hist_mean_site_depth_variants.pdf")
par(mar=c(5,5,4,1)+.1)
plot(hg.targeted, col = c1, ylim=c(0,55), main = "", xlab = "Site Depth") # Plot 1st histogram using a transparent color
plot(hg.nontargeted, col = c2, add = TRUE) # Add 2nd histogram using different color
plot(hg.wPrimers, col = c3, add = TRUE)
plot(hg.oPrimers, col = c4, add = TRUE)
legend("topright", legend=c("targeted SNPs within amplicon ",
                            "non-targeted SNPs within amplicon ",
                            "SNPs within 'primer' region",
                            "SNPs from multiple-aligned reads"), 
       fill=c(c1,c2,c3,c4), bty="n")
dev.off()

## Combined coverage density plots
##------------------------------
pdf(file="results/coverage/calledVariants/dens_mean_site_depth_variants.pdf")
par(mar=c(5,5,4,1)+.1)
plot(NA, xlim = range(ax), ylim = c(0,0.002), 
     xlab = "Site Depth", ylab = "Density", cex.axis=1.2, cex.lab  = 1.2)
lines(dens1 , col="#313695", lwd=2)
lines(dens2 , col="#74add1", lwd=2)
lines(dens4 , col="#f46d43", lwd=2)
lines(dens3 , col="#a50026", lwd=2)

legend("topright", legend=c("targeted SNPs within amplicon ",
                            "non-targeted SNPs within amplicon ",
                            "SNPs within 'primer' region",
                            "SNPs from multiple-aligned reads"), 
       col=c("#313695","#74add1","#f46d43","#a50026"), lwd = c(2,2,2,2), bty="n")
dev.off()


#######################
### SITE DP vs QUAL ###
###  RELATIONSHIP   ###
#######################

## Tables were obtained with bcftools query
## bcftools query -f '%CHROM\t%POS\t%REF,%ALT\t%QUAL\t%NS\t%DP\t%AF\n' snps.vcf > snps_infoTable.txt

## WITHIN AMPLICON - TARGETED variants
##------------------------

# Overall DEPTH and QUAL
targeted.snps.table <- read.table(file = "results/variantCalling/within_amplicon/targeted_snps_infoTable.txt",
                                  na.strings = ".", header = TRUE)
dim(targeted.snps.table)
# 49  7

# Keep the same SNPs as for site-mean depth analyze
targeted.snps.table <- targeted.snps.table[rownames(targeted.snps.dp), ]
dim(targeted.snps.table)
# 47  7

targeted.snps.table <- targeted.snps.table[complete.cases(targeted.snps.table), ]
dim(targeted.snps.table)
# 47  7

targeted.snps.table <- targeted.snps.table[!targeted.snps.table$CHROM == "gi|676159314|gb|JPYR01001074.1|", ] # remove mtDNA
dim(targeted.snps.table)
# 46  7

## WITHIN AMPLICON - NON-TARGETED variants
##------------------------
nontargeted.snps.table <- read.table(file = "results/variantCalling/within_amplicon/nontargeted_snps_infoTable.txt",
                                     na.strings = ".", header = TRUE)
# Overall DEPTH and QUAL
dim(nontargeted.snps.table)
# 171   7

# Keep the same SNPs as for site-mean depth analyze
nontargeted.snps.table <- nontargeted.snps.table[rownames(nontargeted.snps.dp), ]
dim(nontargeted.snps.table)
# 126   7

nontargeted.snps.table <- nontargeted.snps.table[complete.cases(nontargeted.snps.table), ]
dim(nontargeted.snps.table)
# 126   7

nontargeted.snps.table <- nontargeted.snps.table[!nontargeted.snps.table$CHROM == "gi|676159314|gb|JPYR01001074.1|", ] # remove mtDNA
dim(nontargeted.snps.table)
# 125   7

## OUTSIDE AMPLICON - OUTSIDE PRIMER REGION
##------------------------
outsidePrimer.snps.table <- read.table(file = "results/variantCalling/outside_amplicon/variants_outside_outsidePrimers_infoTable.txt",
                                       na.strings = ".", header = TRUE)

dim(outsidePrimer.snps.table)
# 166   7

# List of identified SNPs that were probably produced by a mutation (insertion)
# that moved a bit the primer matching
slipped.primersSNPs

toBeAddedPrimersList2 <- outsidePrimer.snps.table[slipped.primersSNPs, ]

# remove these from the table of nontargeted SNPs outside primer regions
# add these to the list of nontargeted SNPs within primer regions
outsidePrimer.snps.table <- outsidePrimer.snps.table[-slipped.primersSNPs,]
dim(outsidePrimer.snps.table)
# 146   7

# Keep only the same at the end of outsidePrimer.snps.dp dataset
outsidePrimer.snps.table <- outsidePrimer.snps.table[rownames(outsidePrimer.snps.dp), ]
dim(outsidePrimer.snps.table)
# 118   7

## OUTSIDE AMPLICON - WITHIN PRIMER REGION
##------------------------
withinPrimer.snps.table <- read.table(file = "results/variantCalling/outside_amplicon/variants_outside_withinPrimers_infoTable.txt",
                                      na.strings = ".", header = TRUE)

dim(withinPrimer.snps.table)
# 54  7

withinPrimer.snps.table <- rbind(withinPrimer.snps.table, toBeAddedPrimersList2)
dim(withinPrimer.snps.table)
# 74  7

withinPrimer.snps.table <- withinPrimer.snps.table[rownames(withinPrimer.snps.dp), ]
dim(withinPrimer.snps.table)
# 63  7

## Combine all the datasets and calculate 
## the DEPTH and QUAL thresholds
##----------------------------------------
combined.tables <- rbind(targeted.snps.table, nontargeted.snps.table, withinPrimer.snps.table, outsidePrimer.snps.table)

dpCutOff <- mean(combined.tables$DP) + 3*sqrt(mean(combined.tables$DP))
qualCutOff <- 2*mean(combined.tables$DP)

## CORRECT relationship between DEPTH and QUAL
c1.qualDP = ifelse(targeted.snps.table$DP > dpCutOff & targeted.snps.table$QUAL < qualCutOff, c1, "#313695")
c2.qualDP = ifelse(nontargeted.snps.table$DP > dpCutOff & nontargeted.snps.table$QUAL < qualCutOff, c2, "#74add1")
c3.qualDP = ifelse(withinPrimer.snps.table$DP > dpCutOff & withinPrimer.snps.table$QUAL < qualCutOff, c3, "#f46d43")
c4.qualDP = ifelse(outsidePrimer.snps.table$DP > dpCutOff & outsidePrimer.snps.table$QUAL < qualCutOff, c4, "#a50026")

pdf(file="results/variantCalling/qualDepth_relationship.pdf")
par(mar=c(5,5,4,1)+.1)
plot(combined.tables$DP, combined.tables$QUAL, 
     col= c(c1.qualDP,c2.qualDP,c3.qualDP,c4.qualDP), pch=19, cex=1.5,
     xlab="DEPTH", ylab="QUAL", cex.axis = 1.2,
     cex.lab = 1.2)
abline(h=qualCutOff, v=dpCutOff, col="gray", lty="dashed")
#legend("topleft", legend=c("targeted SNPs within amplicon ",
#                            "non-targeted SNPs within amplicon ",
#                            "SNPs within 'primer' region",
#                            "SNPs from multiple-aligned reads"), 
#       col=c("#313695","#74add1","#f46d43","#a50026"), pch = c(19,19,19,19), bty="n")

dev.off()

###########################
### FILTER TARGETED AND ###
###  NONTARGETED SNPs   ###
###########################

## WITHIN AMPLICON - TARGETED variants
##------------------------ 
# all SNPs regardless the relationship
dim(targeted.snps.table)
# 46 + 1 mtDNA
(47/59)*100 # 79.66102% of the targeted SNPs

length(unique(targeted.snps.table$CHROM))
# 29 + 1 mtDNA
(30/48)*100 # 62.5% of the amplicons

# keeping only SNPs within the limits of QUAL vs DEPTH relationship
final.targeted.snps <- targeted.snps.table[!(targeted.snps.table$DP > dpCutOff & targeted.snps.table$QUAL < qualCutOff), ]
dim(final.targeted.snps )
# 38 + 1 mtDNA
(39/59)*100 # 66.10169% of the targeted SNPs

length(unique(final.targeted.snps$CHROM))
# 25 + 1 mtDNA
(26/48)*100 # 54.16667% of the amplicons

## WITHIN AMPLICON - NONTARGETED variants
##------------------------
# all SNPs regardless the relationship
dim(nontargeted.snps.table)
# 125 + 1 mtDNA

# keeping only SNPs within the limits of QUAL vs DEPTH relationship
final.nontargeted.snps <- nontargeted.snps.table[!(nontargeted.snps.table$DP > dpCutOff & nontargeted.snps.table$QUAL < qualCutOff), ]

# remove mono
final.nontargeted.snps <- final.nontargeted.snps[!(final.nontargeted.snps$AF==0 | final.nontargeted.snps$AF==1), ]
dim(final.nontargeted.snps)
# 41  7

## Number of AMPLICONS with Data
##-------------------------------
# Total number of amplicon with at least one SNP - all SNPs
length(unique(c(targeted.snps.table$CHROM, nontargeted.snps.table$CHROM)))
# 35 + 1 mtDNA
(36/48)*100 # 75% of the amplicons

# Total number of amplicon with at least one SNP - QUAL vs DEPTH filtered
length(unique(c(final.targeted.snps$CHROM, final.nontargeted.snps$CHROM)))
# 30 + 1 mtDA
(31/48)*100 # 64.58333% of the amplicons + 1 mtDNA

########################
### ALLELE FREQUENCY ###
###   DISTRIBUTION   ###
########################

## FOR TARGETED and NON TARGETED SNPS WITHIN AMPLICONS
##-----------------------------------------------------

pdf(file="results/variantCalling/alternative_alleleFreq.pdf")
par(mar=c(5,5,4,1)+.1)
plot(NA, xlim = range(0,1), ylim = c(0,2.5), 
     xlab = "Alternative Allele Frequency", ylab = "Density", cex.axis=1.2, cex.lab  = 1.2)
lines(density(final.targeted.snps$AF) , col="#313695", lwd=2)
lines(density(as.numeric(final.nontargeted.snps$AF)) , col="#74add1", lwd=2)
legend("topright", legend=c("targeted SNPs",
                            "non-targeted SNPs"), 
       col=c("#313695","#74add1"), lwd = c(2,2,2,2), bty="n")
dev.off()

################################
### MARKER VALIDATION PART I ###
###     GENOTYPE AND AAF     ###
################################

## FUNCTIONS
calcGenFreq <- function(x)
{
        m = x[, -c(1:2)]
        NS <- apply(!is.na(m),  1, sum)
        countRR <- apply(m==0, 1,sum, na.rm=TRUE)
        countRA <- apply(m==1, 1,sum, na.rm=TRUE)
        countAA <- apply(m==2, 1,sum, na.rm=TRUE)
        
        
        res <- data.frame(freqRR=countRR/NS,
                          freqRA=countRA/NS,
                          freqAA=countAA/NS)
        return(cbind(x[,1:2], NS=NS, res))
}

## WITHIN AMPLICON - TARGETED variants
##------------------------------------

# Import data
targeted.genotypes <- t(read.table(file="results/variantCalling/within_amplicon/genotypes012/targeted_snps.012",
                                 header = F, na.strings = "-1", sep = "\t"))

targeted.genotypes <- targeted.genotypes[-1, ]

targeted.samples <- as.vector(read.csv(file="results/variantCalling/within_amplicon/genotypes012/targeted_snps.012.indv",
                                 header = F))

colnames(targeted.genotypes) <- targeted.samples$V1

targeted.snps <- read.table(file="results/variantCalling/within_amplicon/genotypes012/targeted_snps.012.pos",
                            header = F, sep = "\t", col.names = c("CHROM", "POS"))

targeted.genotypes.table <- cbind(targeted.snps, targeted.genotypes)
dim(targeted.genotypes.table)

# Keep only the genotypes of the SNPs QUAl vs DEPTH filtered
targeted.genotypes.table.fil <- targeted.genotypes.table[paste0(targeted.genotypes.table$CHROM, "_",targeted.genotypes.table$POS) %in% 
                                                         paste0(final.targeted.snps$CHROM, "_",final.targeted.snps$POS), ]
dim(targeted.genotypes.table.fil)
# 38 23

# Calculate Genotype Freq for each marker
targeted.genotypes.freqs <- calcGenFreq(x=targeted.genotypes.table.fil)

## Combinade coverage histograms
##------------------------------
crr1 <- rgb(8,29,88, max = 255, alpha = 90, names = "dkk.blue")
cra2 <- rgb(34,94,168, max = 255, alpha = 90, names = "blueish")
caa3 <- rgb(65,182,196, max = 255, alpha = 90, names = "green.blue")

b <- min(c(targeted.genotypes.freqs$freqRR, targeted.genotypes.freqs$freqRA,
           targeted.genotypes.freqs$freqAA)) # Set the minimum for the breakpoints
e <- max(c(targeted.genotypes.freqs$freqRR, targeted.genotypes.freqs$freqRA,
           targeted.genotypes.freqs$freqAA)) # Set the maximum for the breakpoints
ax <- pretty(b:e, n = 10) # Make a neat vector for the breakpoints
ax

hg.targeted.rr <- hist(targeted.genotypes.freqs$freqRR, breaks = ax, plot = FALSE) 
hg.targeted.ra <- hist(targeted.genotypes.freqs$freqRA, breaks = ax, plot = FALSE) 
hg.targeted.aa <- hist(targeted.genotypes.freqs$freqAA, breaks = ax, plot = FALSE) 

pdf(file="results/variantCalling/targeted_genotypeFreq.pdf")
par(mar=c(5,5,4,1)+.1)
plot(hg.targeted.rr, col = crr1, ylim=c(0,45), main = "", xlab = "Genotype Frequency", ylab = "Counts") # Plot 1st histogram using a transparent color
plot(hg.targeted.ra, col = cra2, add = TRUE) # Add 2nd histogram using different color
plot(hg.targeted.aa, col = caa3, add = TRUE)
legend("topright", legend=c("Homozygous Ref/Ref",
                            "Heterozygous Ref/Alt",
                            "Homozygous Alt/Alt"), 
       fill=c(crr1,cra2,caa3), bty="n")
dev.off()

## WITHIN AMPLICON - NON TARGETED variants
##----------------------------------------

# Import data
nontargeted.genotypes <- t(read.table(file="results/variantCalling/within_amplicon/genotypes012/nontargeted_snps.012",
                                      header = F, na.strings = "-1", sep = "\t"))

nontargeted.genotypes <- nontargeted.genotypes[-1, ]

nontargeted.samples <- as.vector(read.csv(file="results/variantCalling/within_amplicon/genotypes012/nontargeted_snps.012.indv",
                                       header = F))

colnames(nontargeted.genotypes) <- nontargeted.samples$V1

nontargeted.snps <- read.table(file="results/variantCalling/within_amplicon/genotypes012/nontargeted_snps.012.pos",
                            header = F, sep = "\t", col.names = c("CHROM", "POS"))

nontargeted.genotypes.table <- cbind(nontargeted.snps, nontargeted.genotypes)

# Keep only the genotypes of the SNPs QUAl vs DEPTH filtered
nontargeted.genotypes.table.fil <- nontargeted.genotypes.table[paste0(nontargeted.genotypes.table$CHROM, "_",nontargeted.genotypes.table$POS) %in% 
                                                               paste0(final.nontargeted.snps$CHROM, "_",final.nontargeted.snps$POS), ]
dim(nontargeted.genotypes.table.fil)
# 41 23

# Calculate Genotype Freq for each marker
nontargeted.genotypes.freqs <- calcGenFreq(x=nontargeted.genotypes.table.fil)

b <- min(c(nontargeted.genotypes.freqs$freqRR, nontargeted.genotypes.freqs$freqRA,
           nontargeted.genotypes.freqs$freqAA)) # Set the minimum for the breakpoints
e <- max(c(nontargeted.genotypes.freqs$freqRR, nontargeted.genotypes.freqs$freqRA,
           nontargeted.genotypes.freqs$freqAA)) # Set the maximum for the breakpoints
ax <- pretty(b:e, n = 10) # Make a neat vector for the breakpoints
ax

hg.nontargeted.rr <- hist(nontargeted.genotypes.freqs$freqRR, breaks = ax, plot = FALSE) 
hg.nontargeted.ra <- hist(nontargeted.genotypes.freqs$freqRA, breaks = ax, plot = FALSE) 
hg.nontargeted.aa <- hist(nontargeted.genotypes.freqs$freqAA, breaks = ax, plot = FALSE) 

pdf(file="results/variantCalling/nontargeted_genotypeFreq.pdf")
par(mar=c(5,5,4,1)+.1)
plot(hg.nontargeted.rr, col = crr1, ylim=c(0,45), main = "", xlab = "Genotype Frequency", ylab = "Counts") # Plot 1st histogram using a transparent color
plot(hg.nontargeted.ra, col = cra2, add = TRUE) # Add 2nd histogram using different color
plot(hg.nontargeted.aa, col = caa3, add = TRUE)
#legend("topright", legend=c("Homozygous Ref/Ref",
#                            "Heterozygous Ref/Alt",
#                            "Homozygous Alt/Alt"), 
#       fill=c(crr1,cra2,caa3), bty="n")
dev.off()

#################################
### MARKER VALIDATION PART II ###
###   COMPARISON WITH THE WGS ###
#################################

## WGS TARGETED SNPS
##-------------------------------

# UPLOAD GENOTYPE FILES FOR WGS TARGETED SNPS
wgs.targeted.genotypes <- t(read.table(file="/Users/vitorpavinato/Dropbox/PosDoc_OSU_2019/Belgica-genomics-NSF-NERC/WGS/markerDiscovery/results/variantCalling/vcf/genotypes012/targeted.012",
                                       header = F, na.strings = "-1", sep = "\t"))

wgs.targeted.genotypes <- wgs.targeted.genotypes[-1, ]
dim(wgs.targeted.genotypes)
# 58 23

wgs.targeted.samples <- as.vector(read.csv(file="/Users/vitorpavinato/Dropbox/PosDoc_OSU_2019/Belgica-genomics-NSF-NERC/WGS/markerDiscovery/results/variantCalling/vcf/genotypes012/targeted.012.indv",
                                       header = F))

colnames(wgs.targeted.genotypes) <- wgs.targeted.samples$V1

wgs.targeted.snps <- read.table(file="/Users/vitorpavinato/Dropbox/PosDoc_OSU_2019/Belgica-genomics-NSF-NERC/WGS/markerDiscovery/results/variantCalling/vcf/genotypes012/targeted.012.pos",
                                header = F, sep = "\t", col.names = c("CHROM", "POS"))

wgs.targeted.genotypes.table <- cbind(wgs.targeted.snps, wgs.targeted.genotypes)
dim(wgs.targeted.genotypes.table)

# KEEP ONLY THE GENOTYPES RECOVERED WITH THE AMPLICON SEQUENCING
wgs.targeted.genotypes.table.rec <- wgs.targeted.genotypes.table[paste0(wgs.targeted.genotypes.table$CHROM, "_",wgs.targeted.genotypes.table$POS) %in% 
                                                                 paste0(targeted.snps.table$CHROM, "_",targeted.snps.table$POS), ]
dim(wgs.targeted.genotypes.table.rec)
# 46 25

# KEEP ONLY THE SAMPLES ALSO GENOTYPED WITH THE AMPLICON SEQUENCING
wgs.targeted.genotypes.table.rec <- wgs.targeted.genotypes.table.rec[, colnames(wgs.targeted.genotypes.table.rec[,-c(1:2)]) %in% 
                                                                       colnames(targeted.genotypes.table[-c(1:2)])]

dim(wgs.targeted.genotypes.table.rec)
#46 23

# UPLOAD THE FILE WITH INDIVIDUAL SITE DEPTH FOR WGS TARGETED SNPS
wgs.targeted.genotypes.depth <- read.table(file="/Users/vitorpavinato/Dropbox/PosDoc_OSU_2019/Belgica-genomics-NSF-NERC/WGS/markerDiscovery/results/coverage/calledVariants/targeted_DP.txt",
                                  header = F, na.strings = "NA", sep = "\t")
dim(wgs.targeted.genotypes.depth)
# 58 25

colnames(wgs.targeted.genotypes.depth) <- c("CHROM", "POS", wgs.targeted.samples$V1)
dim(wgs.targeted.genotypes.depth)

# KEEP ONLY THE DEPTH FOR GENOTYPES RECOVERED WITH THE AMPLICON SEQUENCING
wgs.targeted.genotypes.depth.rec <- wgs.targeted.genotypes.depth[paste0(wgs.targeted.genotypes.depth$CHROM, "_",wgs.targeted.genotypes.depth$POS) %in% 
                                                                 paste0(targeted.snps.table$CHROM, "_",targeted.snps.table$POS), ]
dim(wgs.targeted.genotypes.depth.rec)
# 46 25

# KEEP ONLY THE SAMPLES ALSO GENOTYPED WITH THE AMPLICON SEQUENCING
wgs.targeted.genotypes.depth.rec <- wgs.targeted.genotypes.depth.rec[, colnames(wgs.targeted.genotypes.depth.rec[,-c(1:2)]) %in% 
                                                                       colnames(targeted.genotypes.table[-c(1:2)])]

dim(wgs.targeted.genotypes.depth.rec)
# 46 23

## AMPLICON SEQ TARGETED SNPS
##-------------------------------

# KEEP ONLY THE GENOTYPES RECOVERED WITH THE AMPLICON SEQUENCING
amplicon.targeted.genotypes.table.rec <- targeted.genotypes.table[paste0(targeted.genotypes.table$CHROM, "_",targeted.genotypes.table$POS) %in% 
                                                                  paste0(targeted.snps.table$CHROM, "_",targeted.snps.table$POS), ]
dim(amplicon.targeted.genotypes.table.rec)
# 46 23

# UPLOAD THE FILE WITH INDIVIDUAL SITE DEPTH FOR AMPLICON SEQ TARGETED SNPS
amplicon.targeted.genotypes.depth <- read.table(file="results/coverage/calledVariants/withinAmplicons/targeted_snps_DP.txt",
                                                header = F, na.strings = ".", sep = "\t")
dim(amplicon.targeted.genotypes.depth)
# 49 23

colnames(amplicon.targeted.genotypes.depth) <- c("CHROM", "POS", targeted.samples$V1)
dim(amplicon.targeted.genotypes.depth)

# KEEP ONLY THE DEPTH FOR GENOTYPES RECOVERED WITH THE AMPLICON SEQUENCING
amplicon.targeted.genotypes.depth.rec <- amplicon.targeted.genotypes.depth[paste0(amplicon.targeted.genotypes.depth$CHROM, "_",amplicon.targeted.genotypes.depth$POS) %in% 
                                                                           paste0(targeted.snps.table$CHROM, "_",targeted.snps.table$POS), ]
dim(amplicon.targeted.genotypes.depth.rec)
# 46 23

## GENOTYPING ERROR RATE AS A FUNCTION OF DEPTH
##----------------------------------------

# WGS TABLES 
dim(wgs.targeted.genotypes.table.rec) # 46 23
dim(wgs.targeted.genotypes.depth.rec) # 46 23

# AMPLICON SEQ TABLES 
dim(amplicon.targeted.genotypes.table.rec) # 46 23
dim(amplicon.targeted.genotypes.depth.rec) # 46 23

n_markers = 46
n_samples = 21

match.wgs.amplicon.genotypes <- matrix(data = NA, nrow = n_markers, ncol = n_samples)
for (i in 1:n_markers)
{
        res = wgs.targeted.genotypes.table.rec[i, -c(1:2)] == amplicon.targeted.genotypes.table.rec[i, -c(1:2)]
        match.wgs.amplicon.genotypes[i, ] = res
}

## BY MARKER
##----------
#errorByMarker <- list()
#for (i in seq_len(n_markers))
#{
#        t <- data.frame(cov=as.numeric(amplicon.targeted.genotypes.depth.rec[i,-c(1:2)]), match=match.wgs.amplicon.genotypes[i,])
#        t <- t[complete.cases(t), ]
#        
#        #obs_pred_class <- data.frame(obs=selection_1, pred=class_selection_1$model.rf$predictions, logthetaPS=logthetaPS) 
#        error <- !t$match
#        tol<-0.2
#        local_error <- array(NA,nrow(t))
#        for (j in seq_len(nrow(t))){
#                
#                distance <- abs(t$cov-t$cov[j])
#                
#                # calculate weigths from epachnikov kernel
#                nacc <- ceiling(length(distance) * tol)
#                ds   <- sort(distance)[nacc]
#                weights <- 1 - (distance/ds)^2
#                weights[which(weights<0)]<-0
#                # calculated weighthed proportion of error
#                local_error[j]<-sum(error*weights)/sum(weights)
#        } # end of j loop
#        
#        t["error"] <- local_error
#        t <- t[order(t$cov),]
#        
#        errorByMarker[[i]] <- t
#
#} # end of i loop
#
# simple plot
#plot(errorByMarker[[12]]$cov,
#     errorByMarker[[12]]$error,
#     xlim = c(1,13335), ylim = c(0,1),
#     xlab="DEPTH",
#     ylab="error rate",
#     type="l",lwd=1,col="gray") 
#
#for (i in 2:n_markers)
#{
#        lines(errorByMarker[[i]]$cov,
#              errorByMarker[[i]]$error,
#              type="l",lwd=1,col="gray")
#            
#}

## OVERALL
##----------
match.wgs.amplicon.table <- data.frame(dp_wgs=matrix(t(wgs.targeted.genotypes.depth.rec[,-c(1:2)])),
                                       dp_amp=matrix(t(amplicon.targeted.genotypes.depth.rec[,-c(1:2)])),
                                       match =matrix(t(match.wgs.amplicon.genotypes)))

match.wgs.amplicon.table <- match.wgs.amplicon.table[complete.cases(match.wgs.amplicon.table), ]

error <- !match.wgs.amplicon.table$match
tol<-0.2
local_error <- array(NA,nrow(match.wgs.amplicon.table))
for (j in seq_len(nrow(match.wgs.amplicon.table))){
        
        distance <- abs(match.wgs.amplicon.table$dp_amp-match.wgs.amplicon.table$dp_amp[j])
        
        # calculate weigths from epachnikov kernel
        nacc <- ceiling(length(distance) * tol)
        ds   <- sort(distance)[nacc]
        weights <- 1 - (distance/ds)^2
        weights[which(weights<0)]<-0
        # calculated weighthed proportion of error
        local_error[j]<-sum(error*weights)/sum(weights)
} # end of j loop

match.wgs.amplicon.table["error"] <- local_error
match.wgs.amplicon.table <- match.wgs.amplicon.table[order(match.wgs.amplicon.table$dp_amp),]

# Error rate as a function of the WGS depth
#pdf(file = "results/gentoypingErrorRate/gentypingError_byWGSdepth.pdf")
pdf(file = "results/gentoypingErrorRate/gentypingError_byAMPdepth.pdf")
par(mar=c(5,5,4,1)+.1)
#plot(match.wgs.amplicon.table$dp_wgs,
plot(match.wgs.amplicon.table$dp_amp,
      match.wgs.amplicon.table$error,
      ylim = c(0,1),
      xlab="Site Depth",
      ylab="Error Rate",
      type="l",lwd=2,col="black"
      ,
      cex.axis=1.2, cex.lab  = 1.2
      ) 
dev.off()

## COUNT THE NUMBER OF SIMILAR GENOTYPES 
##----------------------------------------

colnames(match.wgs.amplicon.genotypes) <- names(amplicon.targeted.genotypes.table.rec[,-c(1:2)])
rownames(match.wgs.amplicon.genotypes) <- seq(1,n_markers,1)

# BY MARKERS
similar.byMarkers <- apply(match.wgs.amplicon.genotypes, 1, function(x) sum(x, na.rm = TRUE))

pdf(file = "results/gentoypingErrorRate/similarities_ByMarkers.pdf")
par(mar=c(5,5,4,1)+.1)
barplot(round((similar.byMarkers/ncol(match.wgs.amplicon.genotypes))*100, 2), 
        ylim = c(0,100), ylab = "% of similar genotypes", las=2, 
        col = c(rep("#999999",length(1:7)),
                "#bd0026","#999999","#bd0026","#bd0026",
                rep("#999999",length(12:23)),
                "#bd0026","#999999","#bd0026",
                "#999999","#999999","#bd0026",
                rep("#999999",length(30:42)),
                "#bd0026","#bd0026",
                "#999999","#999999"
                ),
        cex.lab = 1.2, cex.axis = 1.2, cex.names = 1.2
        )
legend("topright", legend = c("Passed QUAL x DEPTH filter","NOT Passed QUAL x DEPTH filter"), 
       col = c("#999999","#bd0026"), pch = c(15,15), lty = c(0,0), bty = "n", cex = 1)
dev.off()

# BY SAMPLE
similar.bySample <- apply(match.wgs.amplicon.genotypes, 2, function(x) sum(x, na.rm = TRUE))

pdf(file = "results/gentoypingErrorRate/similarities_BySample.pdf")
par(mar=c(5,5,4,1)+.1)
barplot(round((similar.bySample/nrow(match.wgs.amplicon.genotypes))*100, 2), 
        ylim = c(0,100), ylab = "% of similar genotypes", las=2, col = "#999999",
        cex.lab = 1.2, cex.axis = 1.2, cex.names = 1.2)
dev.off()


#######################################
### COUNT MISSING AND HETEROZYGOTES ###
#######################################

## WGS TARGETED SNPS
##-------------------------------

# Count the number of missing data
wgs.missing.markers <- apply(wgs.targeted.genotypes.table.rec[,-c(1:2)], 1, function(x) sum(is.na(x)))/ncol(wgs.targeted.genotypes.table.rec[,-c(1:2)])
wgs.missing.sample <- apply(wgs.targeted.genotypes.table.rec[,-c(1:2)], 2, function(x) sum(is.na(x)))/nrow(wgs.targeted.genotypes.table.rec[,-c(1:2)])

# Count the number of heterozygous genotypes
wgs.heterozygous.markers <- apply(wgs.targeted.genotypes.table.rec[,-c(1:2)], 1, function(x) sum(x == 1, na.rm = TRUE))/ncol(wgs.targeted.genotypes.table.rec[,-c(1:2)])
wgs.heterozygous.sample <- apply(wgs.targeted.genotypes.table.rec[,-c(1:2)], 2, function(x) sum(x == 1, na.rm = TRUE))/nrow(wgs.targeted.genotypes.table.rec[,-c(1:2)])

## AMPLICON SEQ TARGETED SNPS
##-------------------------------

# Count the number of missing data
amplicon.missing.markers <- apply(amplicon.targeted.genotypes.table.rec[,-c(1:2)], 1, function(x) sum(is.na(x)))/ncol(amplicon.targeted.genotypes.table.rec[,-c(1:2)])
amplicon.missing.sample <- apply(amplicon.targeted.genotypes.table.rec[,-c(1:2)], 2, function(x) sum(is.na(x)))/nrow(amplicon.targeted.genotypes.table.rec[,-c(1:2)])

# Count the number of heterozygous genotypes
amplicon.heterozygous.markers <- apply(amplicon.targeted.genotypes.table.rec[,-c(1:2)], 1, function(x) sum(x == 1, na.rm = TRUE))/ncol(amplicon.targeted.genotypes.table.rec[,-c(1:2)])
amplicon.heterozygous.sample <- apply(amplicon.targeted.genotypes.table.rec[,-c(1:2)], 2, function(x) sum(x == 1, na.rm = TRUE))/nrow(amplicon.targeted.genotypes.table.rec[,-c(1:2)])

## MISSING COMPARISON

# BY MARKER
missing.comparison.markers <- cbind(wgs.missing.markers*100, amplicon.missing.markers*100)
rownames(missing.comparison.markers ) <- seq(1,n_markers)

pdf(file = "results/gentoypingErrorRate/missingData_byMarkers.pdf")
par(mar=c(5,5,4,1)+.1)
barplot(t(missing.comparison.markers), 
        beside=T, ylim=c(0,100), ylab="% of missing data", las=2, col=c("#4575b4","#d73027"),
        cex.lab = 1.2, cex.axis = 1.2, cex.names = 1.2)
legend("topright", legend = c("Targeted genotypes in WGS","Targeted genotypes in Ampliconseq"), 
       col = c("#4575b4","#d73027"), pch = c(15,15), lty = c(0,0), bty = "n", cex = 1.2)
dev.off()

# BY SAMPLE
missing.comparison.sample <- cbind(wgs.missing.sample*100, amplicon.missing.sample*100)

pdf(file = "results/gentoypingErrorRate/missingData_bySample.pdf")
par(mar=c(5,5,4,1)+.1)
barplot(t(missing.comparison.sample), 
        beside=T, ylim=c(0,100), ylab="% of missing data", las=2, col=c("#4575b4","#d73027"),
        cex.lab = 1.2, cex.axis = 1.2, cex.names = 1.2)
#legend("topright", legend = c("Targeted genotypes in WGS","Targeted genotypes in Ampliconseq"), 
#       col = c("#4575b4","#d73027"), pch = c(15,15), lty = c(0,0), bty = "n", cex = 1.2)
dev.off()

## N HETEROZYGOUS GENOTYPES COMPARISON

# BY MARKER
heterozygous.comparison.markers <- cbind(wgs.heterozygous.markers*100, amplicon.heterozygous.markers*100)
rownames(heterozygous.comparison.markers ) <- seq(1,n_markers)

pdf(file = "results/gentoypingErrorRate/heterozygous_byMarkers.pdf")
par(mar=c(5,5,4,1)+.1)
barplot(t(heterozygous.comparison.markers), 
        beside=T, ylim=c(0,100), ylab="% of Heterozygous Genotypes", las=2, col=c("#4575b4","#d73027"),
        cex.lab = 1.2, cex.axis = 1.2, cex.names = 1.2)
#legend("topright", legend = c("Targeted genotypes in WGS","Targeted genotypes in Ampliconseq"), 
#       col = c("#4575b4","#d73027"), pch = c(15,15), lty = c(0,0), bty = "n", cex = 1.2)
dev.off()

# BY SAMPLE
heterozygous.comparison.sample <- cbind(wgs.heterozygous.sample*100, amplicon.heterozygous.sample*100)

pdf(file = "results/gentoypingErrorRate/heterozygous_bySample.pdf")
par(mar=c(5,5,4,1)+.1)
barplot(t(heterozygous.comparison.sample), 
        beside=T, ylim=c(0,100), ylab="% of Heterozygous Genotypes", las=2, col=c("#4575b4","#d73027"),
        cex.lab = 1.2, cex.axis = 1.2, cex.names = 1.2)
#legend("topright", legend = c("Targeted genotypes in WGS","Targeted genotypes in Ampliconseq"), 
#       col = c("#4575b4","#d73027"), pch = c(15,15), lty = c(0,0), bty = "n", cex = 1.2)
dev.off()

## Compare the concordance of all genotype types

wgs.targeted.refGeno <- ifelse(is.na(wgs.targeted.genotypes.table.rec[,-c(1:2)]), "NA", 
                           ifelse(wgs.targeted.genotypes.table.rec[,-c(1:2)] == 0, "RR", 
                              ifelse(wgs.targeted.genotypes.table.rec[,-c(1:2)] == 2, "AA", "RA")))

amplicon.targeted.refGeno <- ifelse(is.na(amplicon.targeted.genotypes.table.rec[,-c(1:2)]), "NA", 
                                 ifelse(amplicon.targeted.genotypes.table.rec[,-c(1:2)] == 0, "RR",
                                    ifelse(amplicon.targeted.genotypes.table.rec[,-c(1:2)] == 2, "AA", "RA")))

table(amplicon.targeted.refGeno, wgs.targeted.refGeno)
#                           wgs.targeted.refGeno
#amplicon.targeted.refGeno  AA  NA  RA  RR
#                           AA  19   0  11   1
#                           NA   0   0   2   7
#                           RA  22  72 339 181
#                           RR   4  28  86 194

sum(table(amplicon.targeted.refGeno, wgs.targeted.refGeno))
# 966

((19+339+194)/sum(table(amplicon.targeted.refGeno, wgs.targeted.refGeno)))*100
# true positives = 57.1428%

