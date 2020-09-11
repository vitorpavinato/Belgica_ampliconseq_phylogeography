##########################################
##   READS MAPPED ON AND OFF TARGETS   ###
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

#############################
### READ MAPPED ON TARGET ###
#############################

targeted.samples <- as.vector(read.csv(file="results/variantCalling/within_amplicon/genotypes012/targeted_snps.012.indv",
                                       header = F))


on.targets.files <- read.table(file = paste0("results/coverage/alignedReads","/",targeted.samples[1,1],"_bedCov_onTargets.txt"),
                               header = F, na.strings = ".")

on.targets.files <- on.targets.files[,-c(5:7)]
colnames(on.targets.files) <- c("CHROM", "start", "end", targeted.samples[1,1])

n_samples = 21

for (i in 2:n_samples)
{
  on.targets.temp <- read.table(file = paste0("results/coverage/alignedReads","/",targeted.samples[i,1],"_bedCov_onTargets.txt"),
                                header = F, na.strings = ".")
  
  on.targets.files[, paste(targeted.samples[i,1])] <- on.targets.temp[,4]
  
}

dim(on.targets.files)

# RE-ARRANGE ROWS TO MATCH AMPLICON ORDER
on.targets.files <- on.targets.files[c(2,4,3,1,5,6:14,16,17,15,18,19,
                                       22:25,20,21,27:48,26),
                                     ]
# ADD THE NAMES OF THE AMPLICONS AS ROW NAMES
amplicon_names <- paste0("Amplicon_",seq(1,48,1))

rownames(on.targets.files) <- amplicon_names 

# PLOT BY AMPLICON
#pdf(file="results/coverage/alignedReads/readsByAmplicon.pdf")
#par(mar=c(6,5,3,1)+.1)
boxplot(t(on.targets.files[, -c(1:3)]), ylab = "Number of reads", xaxt="n", xlab="")
axis(1, at=seq(1:48), labels=FALSE)
text(x=seq(1:48), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=colnames(t(on.targets.files[, -c(1:3)])), srt=90, adj=1, xpd=TRUE, cex = 1)
#dev.off()

# PLOT BY SAMPLE
#pdf(file="results/coverage/alignedReads/readsBySample.pdf")
#par(mar=c(6,5,3,1)+.1)
boxplot(on.targets.files[, -c(1:3)], ylab = "Number of reads", xaxt="n", xlab="")
axis(1, at=seq(1:21), labels=FALSE)
text(x=seq(1:21), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=colnames(on.targets.files[, -c(1:3)]), srt=90, adj=1, xpd=TRUE, cex = 1)
#dev.off()


##############################
### READ MAPPED OFF TARGET ###
##############################

off.targets.files <- read.table(file = paste0("results/coverage/alignedReads","/",targeted.samples[1,1],"_bedCov_offTargets.txt"),
                                header = F, na.strings = ".")

off.targets.files <- off.targets.files[,-c(5:7)]
colnames(off.targets.files) <- c("CHROM", "start", "end", targeted.samples[1,1])

n_samples = 21

for (i in 2:n_samples)
{
  off.targets.temp <- read.table(file = paste0("results/coverage/alignedReads","/",targeted.samples[i,1],"_bedCov_offTargets.txt"),
                                 header = F, na.strings = ".")
  
  off.targets.files[, paste(targeted.samples[i,1])] <- off.targets.temp[,4]
  
}

dim(off.targets.files)

colSums(off.targets.files[,-c(1:3)])

# KEEP ONLY ROWS WITH ROW COLSUM > 100
off.targets.files.reduced <- off.targets.files[rowSums(off.targets.files[,-(1:3)]) > 100,]
dim(off.targets.files.reduced)

# ADD THE NAMES OF THE AMPLICONS AS ROW NAMES
offtargets_names <- c("JPYR01000230","JPYR01000252","JPYR01000264","JPYR01000268","JPYR01000640","JPYR01000826","JPYR01000927","JPYR01001001","JPYR01001067",
                      "JPYR01001074","JPYR01001157","JPYR01001234","JPYR01001314","JPYR01001460","JPYR01001654","JPYR01001792","JPYR01001909","JPYR01001941",
                      "JPYR01002107","JPYR01002279","JPYR01002341","JPYR01002728","JPYR01002830","JPYR01003019","JPYR01003058","JPYR01003493","JPYR01003888",
                      "JPYR01003905","JPYR01004396","JPYR01004693")

rownames(off.targets.files.reduced) <- offtargets_names

# PLOT BY AMPLICON
#pdf(file="results/coverage/alignedReads/readsByOffTargets.pdf")
#par(mar=c(6,5,3,1)+.1)
boxplot(t(off.targets.files.reduced[, -c(1:3)]), ylab = "Number of reads", xaxt="n", xlab="")
axis(1, at=seq(1:30), labels=FALSE)
text(x=seq(1:30), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=colnames(t(off.targets.files.reduced[, -c(1:3)])), srt=90, adj=1, xpd=TRUE, cex = 1)
#dev.off()

# PLOT BY SAMPLE
#pdf(file="results/coverage/alignedReads/readsBySampleOffTargets.pdf")
#par(mar=c(6,5,3,1)+.1)
boxplot(off.targets.files.reduced[, -c(1:3)], ylab = "Number of reads", xaxt="n", xlab="")
axis(1, at=seq(1:21), labels=FALSE)
text(x=seq(1:21), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=colnames(off.targets.files.reduced[, -c(1:3)]), srt=90, adj=1, xpd=TRUE, cex = 1)
#dev.off()

## COMBINED PLOT
##-----------------------

#pdf(file="results/coverage/alignedReads/readBy_combined.pdf")
par(mar=c(5,5,2,1)+.1, mfrow=c(2,2))

# PLOT BY AMPLICON - TARGETS
boxplot(t(on.targets.files[, -c(1:3)]), ylab = "Number of reads", xaxt="n", xlab="")
axis(1, at=seq(1:48), labels=FALSE)
text(x=seq(1:48), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=colnames(t(on.targets.files[, -c(1:3)])), srt=90, adj=1, xpd=TRUE, cex = 0.7)

# PLOT BY AMPLICON - OFF TARGETS
boxplot(t(off.targets.files.reduced[, -c(1:3)]), ylab = "Number of reads", xaxt="n", xlab="")
axis(1, at=seq(1:30), labels=FALSE)
text(x=seq(1:30), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=colnames(t(off.targets.files.reduced[, -c(1:3)])), srt=90, adj=1, xpd=TRUE, cex = 0.7)

# PLOT BY SAMPLE - TARGETS
boxplot(on.targets.files[, -c(1:3)], ylab = "Number of reads", xaxt="n", xlab="")
axis(1, at=seq(1:21), labels=FALSE)
text(x=seq(1:21), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=colnames(on.targets.files[, -c(1:3)]), srt=90, adj=1, xpd=TRUE, cex = 0.7)

# PLOT BY SAMPLE - OFF TARGETS
boxplot(off.targets.files.reduced[, -c(1:3)], ylab = "Number of reads", xaxt="n", xlab="")
axis(1, at=seq(1:21), labels=FALSE)
text(x=seq(1:21), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=colnames(off.targets.files.reduced[, -c(1:3)]), srt=90, adj=1, xpd=TRUE, cex = 0.7)

#dev.off()

hist(rowMeans(on.targets.files[,-c(1:3)]), breaks = 25)
hist(rowMeans(off.targets.files.reduced[,-c(1:3)]), breaks = 25)


