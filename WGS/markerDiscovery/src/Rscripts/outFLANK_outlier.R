##############################################
##        Detect selection: outFLANK       ##
##############################################

## Vitor Pavinato
## vitor.pavi@gmail.com

rm(list=ls())
ls()

sessionInfo()$running
sessionInfo()$platform
R.version.string
.Platform$GUI

## outFLANK installation
#install.packages("devtools", dep=TRUE)
#library(devtools)
#source("http://bioconductor.org/biocLite.R")
#biocLite("qvalue")
#install_github("whitlock/OutFLANK")

#install.packages("bigstatsr")
#library(bigstatsr)

#devtools::install_github("privefl/bigsnpr")
#library(bigsnpr)

## Packages
library(OutFLANK)
#library(bigstatsr)
#library(bigsnpr)

sessionInfo()$otherPkgs$OutFLANK$Version
#sessionInfo()$otherPkgs$bigstatsr$Version
#sessionInfo()$otherPkgs$bigsnpr$Version

recordSessionInfo <- sessionInfo()
save(recordSessionInfo, file="results/detect_selection/outflank/sessionInfo.RData")
save.image(file = "results/detect_selection/outflank/workspaceFile.RData")

load(file = "results/detect_selection/outflank/workspaceFile.RData")

## load additional files 
vcf012snps <- read.table(file = "data/datasets/filtered.012.pos", col.names = c("scaffold", "position"))

scaffoldname <- vcf012snps$scaffold
locusnames <- seq(from=1,to=length(vcf012snps$position))

levels(vcf012snps$scaffold) <- seq(from=length(unique(levels(vcf012snps$scaffold))), to=1)

# at this point there are the info for all the snps inlcuding those in the mtDNA
infosnps <- data.frame(scaffold=scaffoldname, 
                       id=as.numeric(vcf012snps$scaffold,decreasing = T), 
                       position=vcf012snps$position, 
                       locusnames=locusnames)


popnames <- c(rep("D1",11), rep("HP",12))

## load genotype files
## The easiest way to convert genotypes to outFLANK input file is to use the 012 output produced by vcftools
## for biallelic SNPs only

snptable <- read.table(file="data/datasets/filtered.012.txt", 
                       col.names = c("ind", paste(infosnps$locusnames)),
                       header = FALSE, na.strings = -1)

snptable <- snptable[,-c(1)]

# remove mtDNA SNPs info
infosnps <- infosnps[-c(417:431),]

# remove mtDNA SNPs genotypes
snptable <- snptable[, -c(417:431)]

# remove snps with NAs to run the identification of redundant snps - high r2
#noNAsnptable <- snptable[, complete.cases(t(snptable))]
#countref <- apply(noNAsnptable, 2, function(x)sum(x == 0))
#countalt <- apply(noNAsnptable, 2, function(x)sum(x == 2))
#keep_snps  <- countref < 20 & countalt < 20
#
#noNAsnptable <- noNAsnptable[, keep_snps]
#noNAinfosnps <- infosnps[infosnps$locusnames %in% row.names(t(noNAsnptable)), ]

# prepare the full snp matrix to run outflank
snptable <- as.matrix(snptable)
snptable[is.na(snptable)] <- 9

# check if NAs were replaced by 9
table(snptable)

## obtaining outFLANK input data - data without mtDNA SNPs
fstdataframe <- MakeDiploidFSTMat(SNPmat = snptable, locusNames = infosnps$locusnames, popNames = popnames)

# data check 1 - full data
pdf(file = "results/detect_selection/outflank/fstHeplot.pdf")
plot(fstdataframe$He, fstdataframe$FST, 
     col=rgb(0,0,0,0.4), pch = 19,
     xlab = expression("H"[E]), ylab = expression("F"[ST]))
dev.off()

# data check 2 - full data
pdf(file = "results/detect_selection/outflank/fstFstCorplot.pdf")
plot(fstdataframe$FST, fstdataframe$FSTNoCorr, 
     col = rgb(0,0,0,0.4), pch = 19, 
     xlab = expression("Corrected F"[ST]), ylab = expression("Uncorrected F"[ST]))
abline(0,1, col = "#3182bd", lty = 2)
dev.off()

## find quasi-independent SNPs to calibrate the Null distribution
#fakesnptable <- cbind(noNAsnptable, noNAsnptable, noNAsnptable)
#fakeinfosnps <- rbind(noNAinfosnps, noNAinfosnps, noNAinfosnps)
#
#fakegenotypes<-add_code256(big_copy(fakesnptable,type="raw"),code=bigsnpr:::CODE_012)
#snppc<-snp_autoSVD(G=fakegenotypes, 
#                   infos.chr = fakeinfosnps$id, 
#                   infos.pos = fakeinfosnps$position)
#
#which_pruned <- attr(newpc, which="subset") # Indexes of remaining SNPS after pruning
#length(which_pruned)

## find quasi-independent SNPs to calibrate the Null distribution
## load inter scaffold::snp r2 calculation - vcftools
## mtDNA SNPs already removed
interr2 <- read.table(file = "results/intra_population/linkage/r2_inter_scaffolds/filtered_interchrom_geno_r2.txt", header = T)

linterr2 <- interr2[which(interr2$R.2 <= 0.01), ]
linterr2c12 <- linterr2[!duplicated(linterr2[, c(1,2)]), c(1,2)]
colnames(linterr2c12) <- c("scaffold", "position")

## load intra scaffold::snp r2 calculation - vcftools
intrar2 <- read.table(file = "results/intra_population/linkage/r2_same_scaffold/filtered_geno_r2.txt", header = TRUE)

## remove mtDNA SNPs
intrar2 <- intrar2[-c(1288:1293), ]

lintrar2 <- intrar2[which(intrar2$R.2 <= 0.01), ]
lintrar2c12 <- lintrar2[!duplicated(lintrar2[, c(1,2)]), c(1,2)]
colnames(lintrar2c12) <- c("scaffold", "position")

lr2infosnps <- rbind(linterr2c12, lintrar2c12)
lr2infosnps <- lr2infosnps[!duplicated(lr2infosnps[,c(1,2)]), ]

infosnpsIndp <- infosnps[paste0(infosnps$scaffold,":",infosnps$position) %in% 
                           paste0(lr2infosnps$scaffold,":",lr2infosnps$position),]
head(infosnpsIndp)
head(fstdataframe)

dim(fstdataframe[infosnpsIndp$locusnames, ])


## run outFLANK
outflankout <- OutFLANK(FstDataFrame = fstdataframe[infosnpsIndp$locusnames, ], 
                        LeftTrimFraction = 0.10, RightTrimFraction = 0.10, Hmin = 0.2,
                        NumberOfSamples = 2, qthreshold = 0.05)

dim(outflankout$results)

OutFLANKResultsPlotter(outflankout, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom = FALSE,
                       RightZoomFraction = 0.05, titletext = NULL)

hist(outflankout$results$pvaluesRightTail)

table(outflankout$results$OutlierFlag)

outlierP <- pOutlierFinderChiSqNoCorr(fstdataframe, Fstbar = outflankout$FSTNoCorrbar, 
                                      dfInferred = outflankout$dfInferred, qthreshold = 0.10, Hmin=0.1)

head(outlierP)
hist(outlierP$pvaluesRightTail)
table(outlierP$OutlierFlag)

finalout <- outlierP$OutlierFlag==TRUE

## SNPs identified as outliers with Bayescan
infosnps[596,]
fstdataframe[596,]
outlierP[which(outlierP$LocusName==611),]
fstdataframe[which(fstdataframe$LocusName==611),]


## final plots
pdf(file = "results/detect_selection/outflank/fstHeplotOutliers.pdf")
plot(outlierP$He, outlierP$FST, 
     col=rgb(0,0,0,0.4), pch = 19,
     xlab = expression("H"[E]), ylab = expression("F"[ST]))
points(outlierP$He[finalout], outlierP$FST[finalout], col="#31a354", pch=19)
#points(outlierP$He[which(outlierP$LocusName==492)], outlierP$FST[which(outlierP$LocusName==492)], col="#d53e4f", pch=19) 
#points(outlierP$He[which(outlierP$LocusName==611)], outlierP$FST[which(outlierP$LocusName==611)], col="#d53e4f", pch=19) 
dev.off()

pdf(file = "results/detect_selection/outflank/fstsnpsplotOutliers.pdf")
plot(outlierP$LocusName, outlierP$FST,
     col=rgb(0,0,0,0.4), pch=19,
     xlab="Position", ylab=expression("F"[ST]))
points(outlierP$LocusName[finalout], outlierP$FST[finalout], col="#31a354", pch=19)
#points(outlierP$LocusName[which(outlierP$LocusName==492)], outlierP$FST[which(outlierP$LocusName==492)], col="#d53e4f", pch=19) 
#points(outlierP$LocusName[which(outlierP$LocusName==611)], outlierP$FST[which(outlierP$LocusName==611)], col="#d53e4f", pch=19) 
dev.off()

pdf(file = "results/detect_selection/outflank/snpsPvalueplotOutliers.pdf")
plot(outlierP$LocusName, -log10(outlierP$pvalues),
     col=rgb(0,0,0,0.4), pch=19,
     xlab="Position", ylab="-log10(p-value)")
points(outlierP$LocusName[finalout], -log10(outlierP$pvalues[finalout]), col="#31a354", pch=19)
#points(outlierP$LocusName[which(outlierP$LocusName==492)], -log10(outlierP$pvalues[which(outlierP$LocusName==492)]), col="#d53e4f", pch=19) 
#points(outlierP$LocusName[which(outlierP$LocusName==611)], -log10(outlierP$pvalues[which(outlierP$LocusName==611)]), col="#d53e4f", pch=19)
abline(h=3, col = "darkgrey", lty = 2)
dev.off()

pdf(file = "results/detect_selection/outflank/snpsQvalueplotOutliers.pdf")
plot(outlierP$LocusName, -log10(outlierP$qvalues),
     col=rgb(0,0,0,0.4), pch=19,
     xlab="Position", ylab="-log10(q-value)")
points(outlierP$LocusName[finalout], -log10(outlierP$qvalues[finalout]), col="#31a354", pch=19)
#points(outlierP$LocusName[which(outlierP$LocusName==492)], -log10(outlierP$qvalues[which(outlierP$LocusName==492)]), col="#d53e4f", pch=19) 
#points(outlierP$LocusName[which(outlierP$LocusName==611)], -log10(outlierP$qvalues[which(outlierP$LocusName==611)]), col="#d53e4f", pch=19)
dev.off()


## New outliers ##
## FDR == 10%   ##

outlierP[which(outlierP$qvalues <= 0.10 ), ]
outlocusnames <- outlierP[which(outlierP$qvalues <= 0.10 ), 1]

infosnps[infosnps$locusnames %in% outlocusnames,]

outlierTable <- data.frame(rownames=rownames(outlierP[which(outlierP$qvalues <= 0.10 ), ]),
                           infosnps[infosnps$locusnames %in% outlocusnames,], 
                           outlierP[which(outlierP$qvalues <= 0.10 ), c(2,3,9,10,11,12)])
write.table(outlierTable, file="results/detect_selection/outflank/outlierFinalTable.txt", sep="\t", quote = F)

save.image(file = "results/detect_selection/outflank/workspaceFile.RData")
