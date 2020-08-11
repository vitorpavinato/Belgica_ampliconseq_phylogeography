#####################################################
## 	    Population Genetics WGS        	   ##
#####################################################

## Vitor Pavinato
## vitor.pavinato@supagro.fr
## INRA

rm(list=ls())
ls()

sessionInfo()$running;
sessionInfo()$platform;
R.version.string;
.Platform$GUI;

#install.packages(c("ape","vegan","adegenet", "HardyWeinberg"), dep=TRUE)
library(ape)
library(vegan)
library(adegenet)
library(pegas)
library(hierfstat)
library(HardyWeinberg)
library(ggplot2)

recordSessionInfo <- sessionInfo();

## Read raw SNP table
# This file contains 1260 SNPs coded asvcftool 012 format:
#  0 = homozigous genotype Ref/Ref;
#  1 = heterozigous genotype Ref/Alt;
#  2 = homozigous genotype Alt/Alt;
# -1 = missing data.
#
# and additional colums: col      1: individuals names;
#                        col      2: factor::Population that each individal belongs;
#                        col 3-1262: SNPs;
#                        col   1263: factor::color

## Load data set
raw.table<-read.table("data/datasets/filtered.012.table.txt", header=TRUE, na.string="-1")
dim(raw.table)
raw.table[1:10,1:10]

# Include sample names
ind.names<-raw.table$Ind.names

# Include population names
pop.ind<-raw.table$pop
pop.names<-c("D_pop", "H_pop")

# Vector of colors
ind.colors <- c(rep("#ff7f00",11),rep("#377eb8",12))

# Take only the columns that contain genotypes
snps<-raw.table[,3:1262]
dim(snps)

# Count missing data:
# SNPs
missing.snps <- apply(snps, 2, function(x) sum(is.na(x))/length(x))
round(mean(missing.snps)*100, 2) #8.53
round(sd(missing.snps)*100, 2) #6.16
round(max(missing.snps)*100, 2) #17.39
round(min(missing.snps)*100, 2) #0

# WGS marker to loci refence list
marker2snp <-matrix(nrow=1260,ncol=2)
marker2snp[,1]<-colnames(snps)
marker2snp[,2]<-paste("SNP",1:ncol(snps),sep="")
write.table(marker2snp, file="results/markers2lociNames",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

# GLOBAL - Check for monomorphic SNPs and for MAF < 0.05
countRR <- apply(snps==0, 2,sum, na.rm=TRUE)
countAA <- apply(snps==2, 2,sum, na.rm=TRUE)
countRA <- apply(snps==1, 2,sum, na.rm=TRUE)

freq.q <- (2*(countAA) + countRA)/((2*(countRR) + countRA)+(2*(countAA) + countRA))

countHomos <- cbind(index=1:ncol(snps),countRR, countAA, q=freq.q)
dim(countHomos[which(countHomos[,4]==0), ]) 
# No SNPs fixed for the Ref allele

dim(countHomos[which(1-countHomos[,4]==0), ])
# No SNPs fixed for the Alt allele

# MAF
maf.threshold = 1/(2 * 23) #0.02173913

dim(countHomos[which(countHomos[,4] < 0.0217 | 1-countHomos[,4] < 0.0217), ])
# 0 SNP with MAF < 1/(2*N[cs])

#sum(countHomos[,4] >= 0.05 & countHomos[,4] <= 0.95) ## VER isso depois 997
#sum(countHomos[,4] >= 0.05) ## VER isso depois  1027
#sum(countHomos[,4] >= 0.0217) ## 1260

# D1 - Check for monomorphic SNPs
countRR.d1 <- apply(snps[1:11,]==0, 2,sum, na.rm=TRUE)
countAA.d1 <- apply(snps[1:11,]==2, 2,sum, na.rm=TRUE)
countRA.d1 <- apply(snps[1:11,]==1, 2,sum, na.rm=TRUE)

freq.q.d1 <- (2*(countAA.d1) + countRA.d1)/((2*(countRR.d1) + countRA.d1)+(2*(countAA.d1) + countRA.d1))

countHomos.d1 <- cbind(index=1:ncol(snps), countRR.d1, countAA.d1, q=freq.q.d1)
dim(countHomos.d1[which(countHomos.d1[,4] == 0), ]) 
# 141 SNPs fixed for the Ref allele

dim(countHomos.d1[which(1-countHomos.d1[,4] == 0), ]) 
# 27 SNPs fixed for the Ref allele

# MAF
maf.threshold.d1 = 1/(2 * 11) #0.04545455

dim(countHomos.d1[which(countHomos.d1[,4] > 0 & countHomos.d1[,4] < 0.05 | 
                                1-countHomos.d1[,4] > 0 & 1-countHomos.d1[,4] < 0.05), ])
# 34 SNPs with maf < 0.05

fixed.d1 <- rbind(countHomos.d1[which(countHomos.d1[,4] == 0), ],
                  countHomos.d1[which(1-countHomos.d1[,4] == 0), ])

write.table(fixed.d1, 
            file="results/intra_population/fixed_snps_d1.txt",
            quote=FALSE, 
            sep="\t",
            row.names=TRUE,
            col.names=TRUE)

# create list with fixed RR and AA, and maf < 0.5 for D1 population
remove.snps.d1 <- c(as.vector(countHomos.d1[which(countHomos.d1[,4] == 0), ][,1]),
                    as.vector(countHomos.d1[which(1-countHomos.d1[,4] == 0), ][,1])
                    #,as.vector(countHomos.d1[which(countHomos.d1[,4] > 0 & countHomos.d1[,4] < 0.05 | 
                    #                                      1-countHomos.d1[,4] > 0 & 1-countHomos.d1[,4] < 0.05), ][,1])
                    )

remove.snps.d1 <- sort(remove.snps.d1)

# HP - Check for monomorphic SNPs
countRR.hp <- apply(snps[12:23,]==0, 2,sum, na.rm=TRUE)
countAA.hp <- apply(snps[12:23,]==2, 2,sum, na.rm=TRUE)
countRA.hp <- apply(snps[12:23,]==1, 2,sum, na.rm=TRUE)

freq.q.hp <- (2*(countAA.hp) + countRA.hp)/((2*(countRR.hp) + countRA.hp)+(2*(countAA.hp) + countRA.hp))

countHomos.hp <- cbind(index=1:ncol(snps), countRR.hp, countAA.hp, q=freq.q.hp)
dim(countHomos.hp[which(countHomos.hp[,4] == 0), ]) 
# 208 SNPs fixed for the Ref allele

dim(countHomos.hp[which(1-countHomos.hp[,4] == 0), ]) 
# 11 SNPs fixed for the Ref allele

countHomos.hp[which(countHomos.hp[,4] == 0), ][,1]
countHomos.hp[which(1-countHomos.hp[,4] == 0), ][,1]

# MAF
maf.threshold.hp = 1/(2 * 12) #0.04166667

dim(countHomos.hp[which(countHomos.hp[,4] > 0 & countHomos.hp[,4] < 0.05 | 
                                1-countHomos.hp[,4] > 0 & 1-countHomos.hp[,4] < 0.05), ])
# 78 SNPs with maf < 0.05

fixed.hp <- rbind(countHomos.hp[which(countHomos.hp[,4] == 0), ],
                  countHomos.hp[which(1-countHomos.hp[,4] == 0), ])

write.table(fixed.hp, 
            file="results/intra_population/fixed_snps_hp.txt",
            quote=FALSE, 
            sep="\t",
            row.names=TRUE,
            col.names=TRUE)

# create list with fixed RR and AA, and maf < 0.5 for D1 population
remove.snps.hp <- c(as.vector(countHomos.hp[which(countHomos.hp[,4] == 0), ][,1]),
                    as.vector(countHomos.hp[which(1-countHomos.hp[,4] == 0), ][,1])
                    #,as.vector(countHomos.hp[which(countHomos.hp[,4] > 0 & countHomos.hp[,4] < 0.05 | 
                    #                                      1-countHomos.hp[,4] > 0 & 1-countHomos.hp[,4] < 0.05), ][,1])
)

remove.snps.hp <- sort(remove.snps.hp)


## Exploratory analysis I - genotypes and allele frequencies
##-----------------------------------------------------------

## Adegenet package
###################

# table to genlight
gen<-new("genlight", snps)

# Acessors
ploidy(gen) #ok all diploidy
nInd(gen)
nLoc(gen)
pop(gen)
other(gen)

indNames(gen)<-ind.names
locNames(gen)<-paste("SNP",1:nLoc(gen),sep="")
pop(gen)<-pop.ind
other(gen)<-ind.colors

#################
# Genotype plot #
#################
# all SNPs - removed mtDNA
pdf(file = "results/intra_population/genotypes_afs/genotype_plot_filtered_snps.pdf")
glPlot(gen[, -c(417:431)], posi="topleft")
dev.off()

###############################
# Allele frequency histograms #
###############################
# all SNPs - removed mtDNA
myFreq<-glMean(gen[, -c(417:431)])
myFreq<-c(myFreq, 1-myFreq)

pdf(file = "results/intra_population/genotypes_afs/allele_freq_distr_filtered_snps_.pdf")
hist(myFreq,
     proba=TRUE,
     col="grey",
     xlab="Allele frequencies",
     main="Distribution of allele frequencies",
     nclass=20)
dens<-density(myFreq, bw=.05)
lines(dens$x, dens$y*2,lwd=1)
dev.off()


# within population
###################

pps<-seppop(gen)

## D1 population

# remove fixed and mtDNA SNPs in D1
remove.snps.d1 <- unique(sort(c(remove.snps.d1, 417:431)))

Dpps<-pps$D_pop
myFreq.1<-glMean(Dpps[ , -remove.snps.d1])
myFreq.1<-c(myFreq.1, 1-myFreq.1)

pdf(file = "results/intra_population/genotypes_afs/allele_freq_distr_filtered_snps_d1.pdf")
hist(myFreq.1,
     proba=TRUE,
     col="grey",
     xlab="Allele frequencies",
     main="Distribution of allele frequencies",
     nclass=20)
dens.1<-density(myFreq.1, bw=.05)
lines(dens.1$x, dens.1$y*2,lwd=1)
dev.off()

## HP population

# remove fixed and mtDNA SNPs in HP
remove.snps.hp <- unique(sort(c(remove.snps.hp, 417:431)))

Hpps<-pps$H_pop
myFreq.2<-glMean(Hpps[, -remove.snps.hp])
myFreq.2<-c(myFreq.2, 1-myFreq.2)

pdf(file = "results/intra_population/genotypes_afs/allele_freq_distr_filtered_snps_hp.pdf")
hist(myFreq.2,
     proba=TRUE,
     col="grey",
     xlab="Allele frequencies",
     main="Distribution of allele frequencies",
     nclass=20)
dens.2<-density(myFreq.2, bw=.05)
lines(dens.2$x, dens.2$y*2,lwd=1)
dev.off()

########################################
#  Allele frequency spectrum and Theta #
########################################

## egglib calculated AFS
########################

## load the egglib output - without mtDNA SNPs
#egglib_summstats <- read.table(file = "results/summary_statistics/filtered_egglib_summstats.txt", header = TRUE)
#afsAll <- egglib_summstats[1, grepl("^SFS", names(egglib_summstats))]
#colnames(afsAll) <- c(1:length(afsAll))
#
## AFS plot
#pdf(file="results/intra_population/genotypes_afs/afs_global_egglib.pdf")
#barplot(as.matrix(afsAll), ylim = c(0,600), ylab = "Number of SNPs", xlab = "", col = "#999999")
#dev.off()

## manually calculated AFS
##########################

# Calculate allele frequencies
calculateAFS <- function(x,n,genes=2){
        
        g0 <- apply(x==0, 2, sum, na.rm=T)
        g1 <- apply(x==1, 2, sum, na.rm=T)
        g2 <- apply(x==2, 2, sum, na.rm=T)
        
        # calculate allele frequencies
        fp <- (2*(g0) + g1)/((2*(g0) + g1)+(2*(g2) + g1))
        fq <- (2*(g2) + g1)/((2*(g0) + g1)+(2*(g2) + g1))
        
        # define the MAF
        mafs <- pmin(fp,fq)
        
        # transform the frequencies in counts
        count_mafs = mafs*n*genes
        
        # create the bins of the folded AFS
        obs_afs <- table((cut(count_mafs,breaks = seq(from=0,to=ceiling(max(count_mafs)),by=1))))
        
        return(obs_afs)
}

# D1 AFS plot
pdf(file="results/intra_population/genotypes_afs/afs_d1.pdf")
barplot(calculateAFS(x=snps[1:11,-remove.snps.d1],n=11,genes=2), ylim = c(0,600),
        ylab = "Number of SNPs", xlab = "", names.arg=seq(1,11,1), col = "#ff7f00")
dev.off()

# HP AFS plot 
pdf(file="results/intra_population/genotypes_afs/afs_hp.pdf")
barplot(calculateAFS(x=snps[12:23,-remove.snps.hp],n=12,genes=2), ylim = c(0,600),
        ylab = "Number of SNPs", xlab = "", names.arg=seq(1,12,1), col = "#377eb8")
dev.off()

# Global AFS plot
pdf(file="results/intra_population/genotypes_afs/afs_global.pdf")
barplot(calculateAFS(x=snps[,-c(417:431)],n=23,genes=2), ylim = c(0,600),
        ylab = "Number of SNPs", xlab = "", names.arg=seq(1,23,1), col = "#999999")
dev.off()

# some thing about Theta
########################

## Number of heterozygous genotypes
## 185,000 out of 88,780,579 positions Kelley et al. 2014
round(100*(195860/88780579),1)
# 0.2% of heterozygotes genotypes

## Ne calculated with one genome, assuming mu = 8.4e-9
(195860/88780579)/(4*8.4e-9)
# 65658.13

## in 23 whole genome sequencing samples
## on average 708 sites out of 1,245 polymorphic site 
## where heterozygous calls in low coverage WGS
mean(apply(snps[,-c(417:431)] == 1, 1, sum, na.rm=TRUE))
# 708.5652

## on average 57% of polymrphic sites were heterozygotes
round(mean(apply(snps[,-c(417:431)] == 1, 1, sum, na.rm=TRUE)/1245)*100,1)
# 56.9

## But, the average number of expected heterozygotes genotypes per sample,
## based on Kelley et al. 2014 would be
0.02* 1245
# 24.9

## Ne calculated with SNPs position that are heterozygotes in 23 samples, assuming mu = 8.4e-9,
## assuming the total site as being the one found with polymirphis in the wgs by Kelley et al. 2014
((1245*0.569)/185000)/(4*8.4e-9)
# 113964.8

calculateTheta <- function(x, size){
        g1 = apply(x == 1, 1, sum, na.rm=TRUE)
        theta = g1/size
        return(theta)
}

global.theta = calculateTheta(x=snps[,-c(417:431)], 185000)
mean(global.theta)
var(global.theta)

# Ne calculated with each genome
mean(global.theta/(4*8.4e-9))
var(global.theta/(4*8.4e-9))

## Now for each population

## D1
## in 11 whole genome sequencing samples
## on average 711 sites out of 1,083 polymorphic site 
## where heterozygous calls in low coverage WGS
mean(apply(snps[1:11,-c(417:431)] == 1, 1, sum, na.rm=TRUE))
# 711.6364

## on average 65% of polymrphic sites were heterozygotes
round(mean(apply(snps[1:11,-remove.snps.d1] == 1, 1, sum, na.rm=TRUE)/1083)*100,1)
# 65.7

## Ne calculated with SNPs position that are heterozygotes in 11 samples, assuming mu = 8.4e-9,
## assuming the total site as being the one found with polymirphis in the wgs by Kelley et al. 2014
((1083*0.657)/185000)/(4*8.4e-9)
# 114467.7

d1.theta = calculateTheta(x=snps[1:11,-c(417:431)], 185000)
mean(d1.theta) #0.003846683
var(d1.theta) #2.914403e-08

# Ne calculated with each genome
mean(d1.theta/(4*8.4e-9)) #114484.6
var(d1.theta/(4*8.4e-9))

## HP
## in 12 whole genome sequencing samples
## on average 705 sites out of 1,034polymorphic site 
## where heterozygous calls in low coverage WGS
mean(apply(snps[12:23,-c(417:431)] == 1, 1, sum, na.rm=TRUE))
# 705.75

## on average 65% of polymrphic sites were heterozygotes
round(mean(apply(snps[12:23,-remove.snps.hp] == 1, 1, sum, na.rm=TRUE)/1034)*100,1)
# 68.3

## Ne calculated with SNPs position that are heterozygotes in 11 samples, assuming mu = 8.4e-9,
## assuming the total site as being the one found with polymirphis in the wgs by Kelley et al. 2014
((1034*0.683)/185000)/(4*8.4e-9)
# 113613.6

hp.theta = calculateTheta(x=snps[12:23,-c(417:431)], 185000)
mean(hp.theta) #0.003814865
var(hp.theta) #8.633375e-09

# Ne calculated with each genome
mean(hp.theta/(4*8.4e-9)) #113537.6
var(hp.theta/(4*8.4e-9))

## Exploratory analysis II - genotype counts and Hardy-Weiberg Equilibrium
##------------------------------------------------------------------------

## genotypes counts
###################

# Genotype counts function
GenotypeCount <- function(snps_table, bySNP=TRUE){
        
        if (bySNP){
                n0 <- apply(snps_table==0, 2, sum, na.rm=T)
                n1 <- apply(snps_table==1, 2, sum, na.rm=T)
                n2 <- apply(snps_table==2, 2, sum, na.rm=T)
        } else {
                n0 <- apply(snps_table==0, 1, sum, na.rm=T)
                n1 <- apply(snps_table==1, 1, sum, na.rm=T)
                n2 <- apply(snps_table==2, 1, sum, na.rm=T)
        }
        res <- as.matrix(data.frame(RR=n0, RA=n1, AA=n2))
        
        return(res)
}

# Genotype counts - by SNP
geno.counts <- GenotypeCount(snps_table = snps, bySNP = TRUE)

# plot total number of each genotype
pdf(file="results/intra_population/genotypes_afs/genotype_counts_bySNP_total.pdf")
barplot(c(sum(geno.counts[-c(417:431),1]), sum(geno.counts[-c(417:431),2]), sum(geno.counts[-c(417:431),3])),
        ylim = c(0,20000),ylab = "Number of genotypes", 
        xlab = "", names.arg=c("Ref/Ref", "Ref/Alt", "Alt/Alt"), col = "#999999")
dev.off()

# remove mtDNA from the geno.counts - for the snp-by-snp genotype plot
geno.counts.reduced <- geno.counts[-c(417:431),]

pdf(file="results/intra_population/genotypes_afs/genotype_counts_bySNP.pdf")
plot(geno.counts.reduced[1,], type="b", ylim = c(0,25), 
     col=ifelse(geno.counts.reduced[1,2] > geno.counts.reduced[1,1], "#999999", "red"),
     xaxt="n", xlab="", ylab="Number of genotypes")
for (i in 2:1245){
        lines(geno.counts.reduced[i,], type="b", 
              col=ifelse(geno.counts.reduced[i,2] > geno.counts.reduced[i,1], "#999999","red")
              )
}
axis(1, at=seq(1:3), labels=FALSE)
text(x=seq(1:3), y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]),
     labels=c("Ref/Ref", "Ref/Alt", "Alt/Alt"), srt=1, adj=0.5, xpd=TRUE)
legend("topright", legend = c("#'s of Ref/Ref > Ref/Alt", "#'s of Ref/Alt > Ref/Ref"), 
       col = c("red", "#999999"), pch = 1, lty = 1, bty = "n", cex = 1)
dev.off()

# Genotype counts - by sample
geno.counts.sample <- GenotypeCount(snps_table = snps[,-c(417:431)], bySNP = FALSE)

pdf(file="results/intra_population/genotypes_afs/genotype_counts_bySample.pdf")
plot(geno.counts.sample[1,], type="b", ylim=c(0,800), 
     col=ifelse(geno.counts.sample[1,2] > geno.counts.sample[1,1], "#999999", "red"),
     xaxt="n", xlab="", ylab="Number of genotypes")
for (i in 2:23){
        lines(geno.counts.sample[i,], type="b", 
              col=ifelse(geno.counts.sample[i,2] > geno.counts.sample[i,1], "#999999","red")
        )
}
axis(1, at=seq(1:3), labels=FALSE)
text(x=seq(1:3), y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]),
     labels=c("Ref/Ref", "Ref/Alt", "Alt/Alt"), srt=1, adj=0.5, xpd=TRUE)
legend("topright", legend = c("#'s of Ref/Ref > Ref/Alt", "#'s of Ref/Alt > Ref/Ref"), 
       col = c("red", "#999999"), pch = 1, lty = 1, bty = "n", cex = 1)
dev.off()

# HardyWeinberg package
#######################

## GLOBAL HWP test
##-----------------

## HWE test for one marker
# Chisq test
HW.test <- HWChisq(geno.counts[1,], verbose = TRUE, cc = 0)

# LR test
HW.lrtest <- HWLratio(geno.counts[1,], verbose = TRUE)

# Fisher's Exact test
HW.exacttest <- HWExact(geno.counts[1,], alternative = "two.sided", pvaluetype = "selome")

## HWE test for many markers
hwe.tests <- HWExactStats(geno.counts, alternative = "two.sided", pvaluetype = "selome")

hwe.table <- data.frame(Markers=rownames(geno.counts), Ex_pvalue=hwe.tests)

# Remove markers in the scaffold JPYR01001074.1 - mtDNA
# Markers #417 to 431
hwe.table <- hwe.table[-c(417:431), ]

# Proportion of markers with departure from HWP
alpha=0.05
# without bonferroni correction
numb.dp.pexact <- sum(hwe.table$Ex_pvalue < alpha) #710
prop.dp.pexact <- (sum(hwe.table$Ex_pvalue < alpha)/length(hwe.table$Markers))*100 # 57.02811

hwe.sign <- hwe.table[which(hwe.table$Ex_pvalue < 0.05), ]

write.table(hwe.sign,
            file      = "results/intra_population/hwe/global_pexact_sign_HardyWeinberg.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            append    = FALSE)


alpha.bc=0.05/dim(hwe.table)[1]
# with bonferroni correction:
numb.dp.pexact.bc <- sum(hwe.table$Ex_pvalue < alpha.bc) #461
prop.dp.pexact.bc <- (sum(hwe.table$Ex_pvalue < alpha.bc)/length(hwe.table$Markers))*100 # 37.02811

hwe.sign.bc <- hwe.table[which(hwe.table$Ex_pvalue < alpha.bc), ]

write.table(hwe.sign.bc,
            file      = "results/intra_population/hwe/global_pexact_sign_bonf_HardyWeinberg.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            append    = FALSE)


# Ternary plot
HWTernaryPlot(geno.counts[-c(417:431), ], vbounds = FALSE, mafbounds = TRUE, mafvalue = 0.05, 
              alpha = alpha, region = 7,curtyp = "dotted", pvaluetype = "selome")

HWTernaryPlot(geno.counts[-c(417:431), ], vbounds = FALSE, mafbounds = TRUE, mafvalue = 0.05, 
              alpha = alpha.bc, region = 7,curtyp = "dotted", pvaluetype = "selome")

# Q-Q plot
HWQqplot(geno.counts[-c(417:431), ], nsim=1000, pvaluetype = "selome", logplot = TRUE)

# Plot p-values for each SNP
plot(hwe.table$Ex_pvalue, pch=19, col=ifelse(hwe.table$Ex_pvalue < alpha.bc, "#1a1a1a", "#bababa"))
abline(h=c(alpha, alpha.bc), col=c("red", "blue"), lty=c(1,1)) # No SNP deviates from HWP 

# Population = D1
#----------------

snps.d1 <- snps[which(pop.ind=="D_pop"), ]
dim(snps.d1)

geno.counts.d1 <- GenotypeCount(snps_table = snps.d1)

hwe.tests.d1 <- HWExactStats(geno.counts.d1, alternative = "two.sided", pvaluetype = "selome")

hwe.table.d1 <- data.frame(Markers=rownames(geno.counts.d1), Ex_pvalue=hwe.tests.d1)

# Remove markers in the scaffold JPYR01001074.1 - mitocondria genome
# Markers #417 to 431
hwe.table.d1 <- hwe.table.d1[-c(417:431), ]

# Proportion of markers with departure from HWP
# without bonferroni correction
numb.dp.pexact.d1 <- sum(hwe.table.d1$Ex_pvalue < alpha) #598
prop.dp.pexact.d1 <- (sum(hwe.table.d1$Ex_pvalue < alpha)/length(hwe.table.d1$Markers))*100 # 48.03213

hwe.sign.d1 <- hwe.table.d1[which(hwe.table.d1$Ex_pvalue < 0.05), ]

write.table(hwe.sign.d1,
            file      = "results/intra_population/hwe/d1_pexact_sign_HardyWeinberg.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            append    = FALSE)

# with bonferroni correction:
numb.dp.pexact.d1.bc <- sum(hwe.table.d1$Ex_pvalue < alpha.bc) #0
prop.dp.pexact.d1.bc <- (sum(hwe.table.d1$Ex_pvalue < alpha.bc)/length(hwe.table.d1$Markers))*100 # 0

# Ternary plot
HWTernaryPlot(geno.counts.d1[-c(417:431), ], vbounds = FALSE, mafbounds = TRUE, mafvalue = 0.05, 
              alpha = alpha, region = 7,curtyp = "dotted", pvaluetype = "selome")

# Q-Q plot
HWQqplot(geno.counts.d1[-c(417:431), ], nsim=1000, pvaluetype = "selome", logplot = TRUE)

# Plot p-values for each SNP
plot(hwe.table.d1$Ex_pvalue, pch=19, col=ifelse(hwe.table.d1$Ex_pvalue < alpha.bc, "#1a1a1a", "#bababa"))
abline(h=c(alpha, alpha.bc), col=c("red", "blue"), lty=c(1,1)) # No SNP deviates from HWP 

# Population = HP
#----------------

snps.hp <- snps[which(pop.ind=="H_pop"), ]
dim(snps.hp)

geno.counts.hp <- GenotypeCount(snps_table = snps.hp)

hwe.tests.hp <- HWExactStats(geno.counts.hp, alternative = "two.sided", pvaluetype = "selome")

hwe.table.hp <- data.frame(Markers=rownames(geno.counts.hp), Ex_pvalue=hwe.tests.hp)

# Remove markers in the scaffold JPYR01001074.1 - mitocondria genome
# Markers #417 to 431
hwe.table.hp <- hwe.table.hp[-c(417:431), ]

# Proportion of markers with departure from HWP
# without bonferroni correction
numb.dp.pexact.hp <- sum(hwe.table.hp$Ex_pvalue < alpha) #629
prop.dp.pexact.hp <- (sum(hwe.table.hp$Ex_pvalue < alpha)/length(hwe.table.hp$Markers))*100 # 50.52209

hwe.sign.hp <- hwe.table.hp[which(hwe.table.hp$Ex_pvalue < 0.05), ]

write.table(hwe.sign.hp,
            file      = "results/intra_population/hwe/hp_pexact_sign_HardyWeinberg.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            append    = FALSE)

# with bonferroni correction:
numb.dp.pexact.hp.bc <- sum(hwe.table.hp$Ex_pvalue < alpha.bc) #0
prop.dp.pexact.hp.bc <- (sum(hwe.table.hp$Ex_pvalue < alpha.bc)/length(hwe.table.hp$Markers))*100 # 0

# Ternary plot
HWTernaryPlot(geno.counts.hp[-c(417:431), ], vbounds = FALSE, mafbounds = TRUE, mafvalue = 0.05, 
              alpha = alpha, region = 7,curtyp = "dotted", pvaluetype = "selome")

# Q-Q plot
HWQqplot(geno.counts.hp[-c(417:431), ], nsim=1000, pvaluetype = "selome", logplot = TRUE)

# Plot p-values for each SNP
plot(hwe.table.hp$Ex_pvalue, pch=19, col=ifelse(hwe.table.hp$Ex_pvalue < alpha.bc, "#1a1a1a", "#bababa"))
abline(h=c(alpha, alpha.bc), col=c("red", "blue"), lty=c(1,1)) # No SNP deviates from HWP 

# Proportion of PHW deviates - without bonferroni correction
prop <- c(prop.dp.pexact, prop.dp.pexact.d1, prop.dp.pexact.hp)
pops <- c("Global", "D1", "HP")

hwe.prop.table <- data.frame(pops=pops, prop=prop)
barplot(hwe.prop.table$prop, names.arg=pops, ylim = c(0,100),
        ylab = "% of departure from HWP")

# Proportion of PHW deviates - with bonferroni correction
prop.bc <- c(prop.dp.pexact.bc, prop.dp.pexact.d1.bc, prop.dp.pexact.hp.bc)
pops.bc <- c("Global", "D1", "HP")

hwe.prop.table.bc <- data.frame(pops=pops.bc, prop=prop.bc)
barplot(hwe.prop.table.bc$prop, names.arg=pops.bc, ylim = c(0,100),
        ylab = "% of departure from HWP")

# Combined Plot 1: ternary plots
pdf(file = "results/intra_population/hwe/HardyWeinberg_ternary_plots.pdf")
par(mfrow=c(2,2))

HWTernaryPlot(geno.counts[-c(417:431), ], vbounds = FALSE, mafbounds = TRUE, mafvalue = 0.05, 
              alpha = alpha, region = 7,curtyp = "dotted", pvaluetype = "selome",
              main = "'Global w/out correction'")

HWTernaryPlot(geno.counts[-c(417:431), ], vbounds = FALSE, mafbounds = TRUE, mafvalue = 0.05, 
              alpha = alpha.bc, region = 7,curtyp = "dotted", pvaluetype = "selome",
              main = "'Global with correction'")

HWTernaryPlot(geno.counts.d1[-c(417:431), ], vbounds = FALSE, mafbounds = TRUE, mafvalue = 0.05, 
              alpha = alpha, region = 7,curtyp = "dotted", pvaluetype = "selome",
              main = "'D1 w/out correction'")

HWTernaryPlot(geno.counts.hp[-c(417:431), ], vbounds = FALSE, mafbounds = TRUE, mafvalue = 0.05, 
              alpha = alpha, region = 7,curtyp = "dotted", pvaluetype = "selome",
              main = "'HP w/out correction'")

dev.off()

# Combined Plot 2: Q-Q plots
pdf(file = "results/intra_population/hwe/HardyWeinberg_qq_plots.pdf")
par(mfrow=c(2,2))

HWQqplot(geno.counts[-c(417:431), ], nsim=1000, pvaluetype = "selome", logplot = TRUE, main = "Q-Q plot - 'Global'")

HWQqplot(geno.counts.d1[-c(417:431), ], nsim=1000, pvaluetype = "selome", logplot = TRUE, main = "Q-Q plot - 'D1'")

HWQqplot(geno.counts.hp[-c(417:431), ], nsim=1000, pvaluetype = "selome", logplot = TRUE, main = "Q-Q plot - 'HP'")

barplot(hwe.prop.table$prop, names.arg=pops.bc, ylim = c(0,100),
        ylab = "% of departure from HWP", main = "% departure w/out correction")

dev.off()

## adegenet and pegas packages
##############################

# Input format 1: GENIND - converting from GENLIGHT
gen.mat <- as.matrix(gen[, -c(417:431)]) # x is a genlight object
gen.mat[gen.mat == 0] <- "1/1" # homozygote reference
gen.mat[gen.mat == 1] <- "1/2" # heterozygote
gen.mat[gen.mat == 2] <- "2/2" # homozygote alternate
gid.gen <- df2genind(gen.mat, sep = "/", ploidy = 2)

# Acessors
ploidy(gid.gen) #ok all diploidy
nInd(gid.gen)
nLoc(gid.gen)
indNames(gid.gen)
locNames(gid.gen)
pop(gid.gen) 
other(gid.gen)

pop(gid.gen) <- pop.ind
other(gid.gen) <- ind.colors
is.genind(gid.gen)

loci.gid <- as.loci(gid.gen)

# HWP test
# Global
hw.gid <- hw.test(loci.gid, B=1000)
dim(hw.gid)
hw.gid[1:10,]

# remove mtDNA SNPs
hw.gid.table <- data.frame(Markers=colnames(snps)[-c(417:431)], hw.gid)

# Proportion of markers with departure from HWP
alpha=0.05
# without bonferroni correction
numb.dp.gid <- sum(hw.gid.table[,5] < alpha) #709
prop.dp.gid <- (sum(hw.gid.table[,5] < alpha)/length(rownames(hw.gid.table)))*100 # 56.94779

hw.gid.sign <- hw.gid.table[which(hw.gid.table[,5] < alpha), ]

write.table(hw.gid.sign,
            file      = "results/intra_population/hwe/global_pexact_sign_pegas.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            append    = FALSE)


alpha.bc=0.05/dim(hw.gid)[1]
# with bonferroni correction:
numb.dp.gid.bc <- sum(hw.gid.table[,5] < alpha.bc) #587
prop.dp.gid.bc <- (sum(hw.gid.table[,5] < alpha.bc)/length(rownames(hw.gid.table)))*100 # 47.14859

hw.gid.sign.bc <- hw.gid.table[which(hw.gid.table[,5] < alpha.bc), ]

write.table(hw.gid.sign.bc,
            file      = "results/intra_population/hwe/global_pexact_sign_bonf_pegas.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            append    = FALSE)

# Q-Q plot
plot(-log10((1:length(hw.gid[,4]))/length(hw.gid[,4])), -log10(hw.gid[,4][order(hw.gid[,4])]),
     xlim = c(0,3), ylim = c(0,3),
     xlab="-log10(expected P)", ylab="-log10(observed P)", main="QQ plot")
abline(a=0, b=1, col="red")

# Plot p-values for each SNP
plot(hw.gid[,4], pch=19, col=ifelse(hw.gid[,4] < alpha.bc, "#1a1a1a", "#bababa"))
abline(h=c(alpha, alpha.bc), col=c("red", "blue"), lty=c(1,1)) 

# Population-specific HWP test - Spliting in subpopulations
pops.gid <- seppop(gid.gen)

# Population = D1
#-----------------
gid.d1<-pops.gid$D_pop
loci.gid.d1 <- as.loci(gid.d1)

# HWP test
hw.gid.d1<-hw.test(loci.gid.d1, B=1000)
dim(hw.gid.d1)
hw.gid.d1[1:10,]

hw.gid.d1.table <- data.frame(Markers=colnames(snps)[-c(417:431)], hw.gid.d1)

# Proportion of markers with departure from HWP
# without bonferroni correction
numb.dp.gid.d1 <- sum(hw.gid.d1.table[,5] < alpha) #594
prop.dp.gid.d1 <- (sum(hw.gid.d1.table[,5] < alpha)/length(rownames(hw.gid.d1.table)))*100 # 47.71084

hw.gid.sign.d1 <- hw.gid.d1.table[which(hw.gid.d1.table[,5] < alpha), ]

write.table(hw.gid.sign.d1,
            file      = "results/intra_population/hwe/d1_pexact_sign_pegas.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            append    = FALSE)


# with bonferroni correction:
numb.dp.gid.d1.bc <- sum(hw.gid.d1.table[,5] < alpha.bc) #13
prop.dp.gid.d1.bc <- (sum(hw.gid.d1.table[,5] < alpha.bc)/length(rownames(hw.gid.d1.table)))*100 # 1.285141

hw.gid.sign.d1.bc <- hw.gid.d1.table[which(hw.gid.d1.table[,5] < alpha.bc), ]

write.table(hw.gid.sign.d1.bc,
            file      = "results/intra_population/hwe/d1_pexact_sign_bonf_pegas.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            append    = FALSE)

# Q-Q plot
plot(-log10((1:length(hw.gid.d1[,4]))/length(hw.gid.d1[,4])), -log10(hw.gid.d1[,4][order(hw.gid.d1[,4])]),
     xlim = c(0,3), ylim = c(0,3),
     xlab="-log10(expected P)", ylab="-log10(observed P)", main="QQ plot")
abline(a=0, b=1, col="red")

# Plot p-values for each SNP
plot(hw.gid.d1[,4], pch=19, col=ifelse(hw.gid.d1[,4] < alpha.bc, "#1a1a1a", "#bababa"))
abline(h=c(alpha, alpha.bc), col=c("red", "blue"), lty=c(1,1)) 

# Population = HP
#------------------
gid.hp<-pops.gid$H_pop
loci.gid.hp <- as.loci(gid.hp)

# HWP test
hw.gid.hp<-hw.test(loci.gid.hp, B=1000)
dim(hw.gid.hp)
hw.gid.hp[1:10,]

hw.gid.hp.table <- data.frame(Markers=colnames(snps)[-c(417:431)], hw.gid.hp)

# Proportion of markers with departure from HWP
# without bonferroni correction
numb.dp.gid.hp <- sum(hw.gid.hp.table[,5] < alpha) #625
prop.dp.gid.hp <- (sum(hw.gid.hp.table[,5] < alpha)/length(rownames(hw.gid.hp.table)))*100 # 50.2008

hw.gid.sign.hp <- hw.gid.hp.table[which(hw.gid.hp.table[,5] < alpha), ]

write.table(hw.gid.sign.hp,
            file      = "results/intra_population/hwe/hp_pexact_sign_pegas.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            append    = FALSE)

# with bonferroni correction:
numb.dp.gid.hp.bc <- sum(hw.gid.hp.table[,5] < alpha.bc) #49
prop.dp.gid.hp.bc <- (sum(hw.gid.hp.table[,5] < alpha.bc)/length(rownames(hw.gid.hp.table)))*100 # 3.935743

hw.gid.sign.hp.bc <- hw.gid.hp.table[which(hw.gid.hp.table[,5] < alpha.bc), ]

write.table(hw.gid.sign.hp.bc,
            file      = "results/intra_population/hwe/hp_pexact_sign_bonf_pegas.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            append    = FALSE)

# Q-Q plot
plot(-log10((1:length(hw.gid.hp[,4]))/length(hw.gid.hp[,4])), -log10(hw.gid.hp[,4][order(hw.gid.hp[,4])]),
     xlim = c(0,3), ylim = c(0,3), 
     xlab="-log10(expected P)", ylab="-log10(observed P)", main="QQ plot")
abline(a=0, b=1, col="red")

# Plot p-values for each SNP
plot(hw.gid.hp[,4], pch=19, col=ifelse(hw.gid.hp[,4] < alpha.bc, "#1a1a1a", "#bababa"))
abline(h=c(alpha, alpha.bc), col=c("red", "blue"), lty=c(1,1)) 

## Exploratory analysis III - heterozygosity and genetic diversity
##-------------------------------------------------------------------------

## ADEGENET - Loci Heterozygosity
#################################

# DIFFERENCES BETWEEN HE AND HO

# GLOBAL
#----------------
global_stats<-summary(gid.gen)

global_stats_table <- data.frame(Markers=colnames(snps)[-c(417:431)], Loci=locNames(gid.gen), 
                                 HE=global_stats$Hexp, HO=global_stats$Hobs, deltaHEHO =  (global_stats$Hexp - global_stats$Hobs))

write.table(global_stats_table,
            file      = "results/intra_population/heterozygosity/global_HEHO_loci_table_adegenet.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            append    = FALSE)

pdf(file="results/intra_population/heterozygosity/global_HOHE_plot_adegenet.pdf")
plot(global_stats_table$HE ~ global_stats_table$HO, col=rgb(0,0,0,0.4), pch=19, 
     ylab="Expected heterozygosity", xlab="Observed heterozygosity", 
     xlim=c(0,1),ylim=c(0,1))
abline(a=0, b=1, lty=3, col='red', lwd=1)
dev.off()

# take the SNPs which HO < HE
write.table(global_stats_table[which(global_stats_table$HO < global_stats_table$HE), ],
            file      = "results/intra_population/heterozygosity/global_HOsmallHE_table_adegenet.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            append    = FALSE)

bartlett.test(list(global_stats_table$HE,global_stats_table$HO))
#Bartlett test of homogeneity of variances
#
#data:  list(global_stats_table$HE, global_stats_table$HO)
#Bartlett's K-squared = 655.59, df = 1, p-value < 2.2e-16

# Mean and variance
mean(global_stats_table$HO) # 0.615991
var(global_stats_table$HO) # 0.1652393

mean(global_stats_table$HE) # 0.3444635
var(global_stats_table$HE) # 0.03625761

# additional visualization plots
plot(global_stats_table$HO, col=adjustcolor("#de2d26", alpha.f = 0.8), pch=1, xlab="SNP", ylab="Heterozygosity")
points(global_stats_table$HE, col=adjustcolor("#3182bd", alpha.f = 0.6), pch=1)

pdf(file="results/intra_population/heterozygosity/global_deltaHet_plot_adegenet.pdf")
plot(global_stats_table$deltaHEHO, col="black", type="l",
     ylim=c(-0.5, 0.5), xlab="SNP", ylab=expression(paste(Delta," Heterozygosity")))
dev.off()


# Population = D1
#----------------

# remove remaining monomorphic markers - mtDNA had been removed before
remove.snps.d1.nomtdna <- remove.snps.d1[-c(54:68)]

all.d1.loci.names <- locNames(gid.d1)

keep.loci.d1 <- setdiff(all.d1.loci.names, paste0("SNP",remove.snps.d1.nomtdna))

gid.d1.poly <- gid.d1[loc = keep.loci.d1]

d1_stats<-summary(gid.d1.poly)

d1_stats_table <- data.frame(Markers=colnames(snps)[-remove.snps.d1], 
                             Loci=locNames(gid.d1.poly), 
                             HE=d1_stats$Hexp, HO=d1_stats$Hobs, 
                             deltaHEHO = (d1_stats$Hexp - d1_stats$Hobs))

write.table(d1_stats_table,
            file      = "results/intra_population/heterozygosity/d1_HEHO_loci_table_adegenet.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            append    = FALSE)

pdf(file="results/intra_population/heterozygosity/d1_HOHE_plot_adegenet.pdf")
plot(d1_stats_table$HE ~ d1_stats_table$HO, col=rgb(0,0,0,0.4), pch=19, 
     ylab="Expected heterozygosity", xlab="Observed heterozygosity", 
     xlim=c(0,1),ylim=c(0,1))
abline(a=0, b=1, lty=3, col='red', lwd=1)
dev.off()

# take the SNPs which HO < HE
write.table(d1_stats_table[which(d1_stats_table$HO < d1_stats_table$HE), ],
            file      = "results/intra_population/heterozygosity/d1_HOsmallHE_table_adegenet.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            append    = FALSE)

bartlett.test(list(d1_stats_table$HE, d1_stats_table$HO))
#Bartlett test of homogeneity of variances
#
#data:  list(d1_stats_table$HE, d1_stats_table$HO)
#Bartlett's K-squared = 689.94, df = 1, p-value < 2.2e-16

# Mean and variance
mean(d1_stats_table$HO) # 0.7154042
var(d1_stats_table$HO) # 0.1282169

mean(d1_stats_table$HE) # 0.3959183
var(d1_stats_table$HE) # 0.02381147

# additional visualization plots
plot(d1_stats_table$HO, col=adjustcolor("#de2d26", alpha.f = 0.8), pch=1, xlab="SNP", ylab="Heterozygosity")
points(d1_stats_table$HE, col=adjustcolor("#3182bd", alpha.f = 0.6), pch=1)

pdf(file="results/intra_population/heterozygosity/d1_deltaHet_plot_adegenet.pdf")
plot(d1_stats_table$deltaHEHO, col="black", type="l",
     ylim=c(-0.5, 0.5), xlab="SNP", ylab=expression(paste(Delta," Heterozygosity")))
dev.off()

# Population = HP
#----------------

# remove remaining monomorphic markers - mtDNA had been removed before
remove.snps.hp.nomtdna <- remove.snps.hp[-c(37:51)]

all.hp.loci.names <- locNames(gid.hp)

keep.loci.hp <- setdiff(all.hp.loci.names, paste0("SNP",remove.snps.hp.nomtdna))

gid.hp.poly <- gid.hp[loc = keep.loci.hp]

hp_stats<-summary(gid.hp.poly)

hp_stats_table <- data.frame(Markers=colnames(snps)[-remove.snps.hp], 
                             Loci=locNames(gid.hp.poly), 
                             HE=hp_stats$Hexp, HO=hp_stats$Hobs, 
                             deltaHEHO = (hp_stats$Hexp - hp_stats$Hobs))

write.table(hp_stats_table,
            file      = "results/intra_population/heterozygosity/hp_HEHO_loci_table_adegenet.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            append    = FALSE)

pdf(file="results/intra_population/heterozygosity/hp_HOHE_plot_adegenet.pdf")
plot(hp_stats_table$HE ~ hp_stats_table$HO, col=rgb(0,0,0,0.4), pch=19, 
     ylab="Expected Heterozygosity", xlab="Observed Heterozygosity", 
     xlim=c(0,1),ylim=c(0,1))
abline(a=0, b=1, lty=3, col='red', lwd=1)
dev.off()

write.table(hp_stats_table[which(hp_stats_table$HO < hp_stats_table$HE), ],
            file      = "results/intra_population/heterozygosity/hp_HOsmallHE_table_adegenet.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            append    = FALSE)

bartlett.test(list(hp_stats_table$HE, hp_stats_table$HO))
#Bartlett test of homogeneity of variances
#
#data:  list(hp_stats_table$HE, hp_stats_table$HO)
#Bartlett's K-squared = 648.62, df = 1, p-value < 2.2e-16

# Mean and variance
mean(hp_stats_table$HO) # 0.7386625
var(hp_stats_table$HO) # 0.1266963

mean(hp_stats_table$HE) # 0.4036168
var(hp_stats_table$HE) # 0.02386642

# additional visualization plots
plot(hp_stats_table$HO, col=adjustcolor("#de2d26", alpha.f = 0.8), pch=1, xlab="SNP", ylab="Heterozygosity")
points(hp_stats_table$HE, col=adjustcolor("#3182bd", alpha.f = 0.6), pch=1)

pdf(file="results/intra_population/heterozygosity/hp_deltaHet_plot_adegenet.pdf")
plot(hp_stats_table$deltaHEHO, col="black", type="l",
     ylim=c(-0.5, 0.5), xlab="SNP", ylab=expression(paste(Delta," Heterozygosity")))
dev.off()

## ADEGENET/GENPOP - Inbreeding and mean Heterozygosity
#######################################################

# GLOBAL
#----------------

# INDIVIDUAL INBREEDING COEFFICIENT
ind_inbreeding <- inbreeding(gid.gen, pop = NULL, res.type = "sample", N=1000, M=2000)

ind_inbreeding_vector <- do.call(rbind, lapply(ind_inbreeding, function(x) mean(x)))

ind_inbreeding_table <- data.frame(pop=pop(gid.gen), ind=indNames(gid.gen), FIS = ind_inbreeding_vector)

FIS_pop <- aggregate(FIS ~ pop, mean, data = ind_inbreeding_table)

# POPULATION'S EXPECTED HETEROZYGOSITY
HE_pop <- Hs(genind2genpop(gid.gen))

# TEST DIFFERENCES BETWEEN  HE
Hs.test(gid.d1, gid.hp, 1000)
#Monte-Carlo test
#Call: Hs.test(x = gid.d1, y = gid.hp, n.sim = 1000)
#
#Observation: 0.009188523 
#
#Based on 1000 replicates
#Simulated p-value: 0.07992008 
#Alternative hypothesis: two-sided 
#
#Std.Obs.y   Expectation      Variance 
#1.725499e+00 -3.819141e-05  2.859334e-05 

# Population = D1
#----------------

ind.inbreeding.d1 <- inbreeding(gid.d1.poly, pop = NULL, res.type = "sample", N=1000, M=2000)

inbreeding.d1.vector <- do.call(rbind, lapply(ind.inbreeding.d1, function(x) mean(x)))

FIS.d1.poly <- mean(inbreeding.d1.vector)

HE.d1.poly <- Hs(genind2genpop(gid.d1.poly))

# Population = HP
#----------------

ind.inbreeding.hp <- inbreeding(gid.hp.poly, pop = NULL, res.type = "sample", N=1000, M=2000)

inbreeding.hp.vector <- do.call(rbind, lapply(ind.inbreeding.hp, function(x) mean(x)))

FIS.hp.poly <- mean(inbreeding.hp.vector)

HE.hp.poly <- Hs(genind2genpop(gid.hp.poly))

Pop_fis_HE_table <- data.frame(Pop=popNames(gid.gen), HE = HE_pop, HE_poly = c(HE.d1.poly, HE.hp.poly), 
                                                      FIS = FIS_pop$FIS, FIS_poly = c(FIS.d1.poly, FIS.hp.poly))
write.table(Pop_fis_HE_table,
            file      = "results/intra_population/heterozygosity/populations_he_fis_table_adegenet.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            append    = FALSE)


## ADEGENET/FSTAT - Heterozygosity: loci and mean
########################################

hfstat.gid <- genind2hierfstat(gid.gen, pop=gid.gen@pop)

hfstat_stat <- basic.stats(hfstat.gid, diploid = TRUE, digits = 8) 

write.table(hfstat_stat$perloc,
            file      = "results/intra_population/heterozygosity/global_Diver_loci_table_fstat.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = TRUE,
            append    = FALSE)

write.table(hfstat_stat$overall,
            file      = "results/intra_population/heterozygosity/global_Diver_overall_table_fstat.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = TRUE,
            append    = FALSE)

hfstat_ar <- allelic.richness(hfstat.gid, min.n=NULL, diploid=TRUE)

# AR visualization
plot(hfstat_ar$Ar[,1], col="#ff7f00", xlab="SNPs", ylab="Allele richness")
points(hfstat_ar$Ar[,2], col="#377eb8")

# D1 summary table
#-----------------
summary_loci_d1_table <- data.frame(N=hfstat_stat$n.ind.samp[,1], AR=hfstat_ar$Ar[,1], 
                                    HE=hfstat_stat$Hs[,1], HO=hfstat_stat$Ho[,1], FIS=hfstat_stat$Fis[,1])

summstats.mean.d1.hfstat <- apply(summary_loci_d1_table, 2, mean, na.rm=TRUE)

summary_loci_d1_table <- summary_loci_d1_table[-remove.snps.d1.nomtdna,]

write.table(summary_loci_d1_table,
            file      = "results/intra_population/heterozygosity/d1_summary_loci_table_fstat.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = TRUE,
            append    = FALSE)

summstats.mean.d1.poly.hfstat <- apply(summary_loci_d1_table, 2, mean, na.rm=TRUE)

# HP summary table
#-----------------
summary_loci_hp_table <- data.frame(N=hfstat_stat$n.ind.samp[,2], AR=hfstat_ar$Ar[,2], 
                                    HE=hfstat_stat$Hs[,2], HO=hfstat_stat$Ho[,2], FIS=hfstat_stat$Fis[,2])

summstats.mean.hp.hfstat <- apply(summary_loci_hp_table, 2, mean, na.rm=TRUE)

summary_loci_hp_table <- summary_loci_hp_table[-remove.snps.hp.nomtdna,]

write.table(summary_loci_hp_table,
            file      = "results/intra_population/heterozygosity/hp_summary_loci_table_fstat.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = TRUE,
            append    = FALSE)

summstats.mean.hp.poly.hfstat <- apply(summary_loci_hp_table, 2, mean, na.rm=TRUE)

write.table(rbind(summstats.mean.d1.hfstat,
                  summstats.mean.hp.hfstat,
                  summstats.mean.d1.poly.hfstat,
                  summstats.mean.hp.poly.hfstat),
            file      = "results/intra_population/heterozygosity/summary_d1hp_mean_table_fstat.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = TRUE,
            append    = FALSE)

## ADEGENET/PEGAS - Other statistics
####################################

# Global
#--------
pegas_stats <- summary(loci.gid)

# Heterozygosity
pegas_H_global <- sapply(pegas_stats, function(x) heterozygosity(x$allele))

# Theta H
pegas_thetaH_global <- sapply(pegas_stats, function(x) theta.h(x$allele))

# Theta K
pegas_thetaK_global <- sapply(pegas_stats, function(x) theta.k(x$allele))

pegas_global_table <- data.frame(Het=pegas_H_global, thetaH=pegas_thetaH_global, thetaK=pegas_thetaK_global)

write.table(pegas_global_table,
            file      = "results/intra_population/heterozygosity/global_Het_thetas_loci_table_pegas.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = TRUE,
            append    = FALSE)

# D1 population
#--------------
pegas_stats_d1 <- summary(loci.gid.d1)

# Heterozygosity
pegas_H_d1 <- sapply(pegas_stats_d1, function(x) heterozygosity(x$allele))

# Theta H
pegas_thetaH_d1 <- sapply(pegas_stats_d1, function(x) theta.h(x$allele))

pegas_d1_table <- data.frame(Het=pegas_H_d1, thetaH=pegas_thetaH_d1)

# only polymorphic SNPs
loci.gid.d1.poly <- as.loci(gid.d1.poly)

pegas_stats_d1_poly <- summary(loci.gid.d1.poly)

# Heterozygosity
pegas_H_d1_poly <- sapply(pegas_stats_d1_poly, function(x) heterozygosity(x$allele))

# Theta H
pegas_thetaH_d1_poly <- sapply(pegas_stats_d1_poly, function(x) theta.h(x$allele))

pegas_d1_table_poly <- data.frame(Het=pegas_H_d1_poly, thetaH=pegas_thetaH_d1_poly)

write.table(pegas_d1_table_poly,
            file      = "results/intra_population/heterozygosity/d1_Het_thetas_loci_table_pegas.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = TRUE,
            append    = FALSE)

# HP population
#--------------
pegas_stats_hp <- summary(loci.gid.hp)

# Heterozygosity
pegas_H_hp <- sapply(pegas_stats_hp, function(x) heterozygosity(x$allele))

# Theta H
pegas_thetaH_hp <- sapply(pegas_stats_hp, function(x) theta.h(x$allele))

pegas_hp_table <- data.frame(Het=pegas_H_hp, thetaH=pegas_thetaH_hp)

# only polymorphic SNPs
loci.gid.hp.poly <- as.loci(gid.hp.poly)

pegas_stats_hp_poly <- summary(loci.gid.hp.poly)

# Heterozygosity
pegas_H_hp_poly <- sapply(pegas_stats_hp_poly, function(x) heterozygosity(x$allele))

# Theta H
pegas_thetaH_hp_poly <- sapply(pegas_stats_hp_poly, function(x) theta.h(x$allele))

pegas_hp_table_poly <- data.frame(Het=pegas_H_hp_poly, thetaH=pegas_thetaH_hp_poly)

write.table(pegas_hp_table_poly,
            file      = "results/intra_population/heterozygosity/hp_Het_thetas_loci_table_pegas.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = TRUE,
            append    = FALSE)

summstats_pegas <- rbind(apply(pegas_global_table, 2, mean, na.rm=TRUE)[c(1,3)],
                         apply(pegas_d1_table, 2, mean, na.rm=TRUE),
                         apply(pegas_d1_table_poly, 2, mean, na.rm=TRUE),
                         apply(pegas_hp_table, 2, mean, na.rm=TRUE),
                         apply(pegas_hp_table_poly, 2, mean, na.rm=TRUE))

rownames(summstats_pegas) <- c("global",
                               "d1",
                               "d1_poly",
                               "hp",
                               "hp_poly")

write.table(summstats_pegas,
            file      = "results/intra_population/heterozygosity/summary_d1hp_mean_table_pegas.txt",
            sep       = "\t",
            quote     = FALSE,
            col.names = TRUE,
            row.names = TRUE,
            append    = FALSE)

## ADEGENET - AFS
#################

# Global AFS
global.maf <- minorAllele(gid.gen)

pdf(file="results/intra_population/genotypes_afs/afs_global_adegenet.pdf")
hist(global.maf, col ="#999999", xlab = "AFS", main = "")
dev.off()

global_maf_table <- data.frame(Pop="Global", maf=global.maf)

# D1 AFS
d1.maf <- minorAllele(gid.d1.poly)

pdf(file="results/intra_population/genotypes_afs/afs_d1_adegenet.pdf")
hist(d1.maf, col ="#ff7f00", xlab = "AFS", main = "")
dev.off()

d1_maf_table <- data.frame(Pop="D1", maf=d1.maf)

# HP AFS
hp.maf <- minorAllele(gid.hp.poly)

pdf(file="results/intra_population/genotypes_afs/afs_hp_adegenet.pdf")
hist(hp.maf, col ="#377eb8", xlab = "AFS", main = "")
dev.off()

hp_maf_table <- data.frame(Pop="HP", maf=hp.maf)

maf_table <- rbind(global_maf_table,d1_maf_table,hp_maf_table)

## MAF plot
pdf(file = "results/intra_population/genotypes_afs/afs_populations_adegenet.pdf", width = 11, height = 9)
p1 <- ggplot(maf_table, aes(maf, fill = Pop)) 
p1 <- p1 + geom_histogram(alpha = 0.5, aes(y = ..density..),position = 'identity')
p1 <- p1 + xlab("MAF") + ylab("Density")
p1 <- p1 + scale_fill_manual(values=c("#999999", "#ff7f00", "#377eb8"), 
                             name="Populations: ",
                             breaks=c("Global", "D1", "HP"),
                             labels=c("Global", "D1", "HP"))

p1 <- p1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                              legend.position = "right", axis.text = element_text(size = 12),
                              axis.title=element_text(size=14))
print(p1)
dev.off()

## Exploratory analysis IV - data structure - distances, heatmaps and PCAs
##-------------------------------------------------------------------------

snps.paste<-snps
dim(snps.paste)

#rownames(snps.paste)<-ind.names
rownames(snps.paste) <- c(sprintf("%s%d", "D1-0", 1:8), sprintf("%s%d", "D1-", 10:12),
                          sprintf("%s%d", "HP-0", 1:9), sprintf("%s%d", "HP-", 10:12))
snps.paste[1:10,1:10]

colnames(snps.paste)<-paste("SNP_",1:ncol(snps.paste),sep="") # This changed from SNP.01 to SNP1

# remove mtDNA SNPs
snps.paste <- snps.paste[, -c(417:431)]
dim(snps.paste)

######################
# Euclidean Distance #
######################

dist<-vegdist(snps.paste, method = "euclidean", na.rm=TRUE)

pdf(file = "results/inter_population/distance/nj_euclidian_ind_dist_filtered_snps.pdf")
nj.d<-nj(dist)
plot(nj.d)
dev.off()

pdf(file = "results/inter_population/distance/ward_euclidian_ind_dist_filtered_snps.pdf")
ward<-hclust(dist, method = "ward.D")
plot(ward)
dev.off()

###########
# Heatmap #
###########

pdf(file = "results/inter_population/distance/heatmap_filtered_snps.pdf")
heat.obj<-heatmap(as.matrix(snps.paste))
dev.off()
#str(heat.obj)

#################################
# Principal componente analysis #
#################################

#####
## Price et al 2010 - vcov of matrix individuals x markers
#####

W<-scale(snps.paste, scale=FALSE) #centering
W[1:10,1:10]
W[is.na(W)]<-0
cov.W<-cov(t(W))
eig.result<-eigen(cov.W)
eig.vec<-eig.result$vectors
lambda<-eig.result$values

# PVE
pdf(file = "results/inter_population/pca/fraction_variance_pca_price-et-al-2010.pdf")
plot(lambda/sum(lambda),ylab="Fraction of total variance", ylim=c(0,0.5), type='b')
dev.off()
100*lambda[1]/sum(lambda) #16.39946
100*lambda[2]/sum(lambda) #6.774808

# PCA plot
pdf(file = "results/inter_population/pca/pca_price-et-al-2010.pdf")
plot(eig.vec[,1],eig.vec[,2],col=ind.colors, 
     xlim = c(-0.4, 0.4), ylim = c(-0.4, 0.4),
     xlab="PC1 (16.40%)", ylab="PC2 (6.77%)", cex=1.5, pch=19)
legend("bottomright", legend = c("D1", "HP"), 
       col = c("#ff7f00", "#377eb8"), pch = 19, bty = "n", cex = 1)
abline(v=0,h=0,col="grey",lty=3)
dev.off()

#####
## precomp package - SVD decomposition of COV matrix (Price)
#####

pc.precomp<-prcomp(cov.W, center=F, scale=F, retx=TRUE)
names(pc.precomp)

lambda.2<-pc.precomp$sdev

# PVE
pdf(file = "results/inter_population/pca/fraction_variance_pca_precomp.pdf")
plot(lambda.2/sum(lambda.2),ylab="Fraction of total variance", ylim=c(0,0.5), type='b')
dev.off()
round(100*lambda.2[1]/sum(lambda.2),2) #16.39946
round(100*lambda.2[2]/sum(lambda.2),2)#6.774808

pdf(file = "results/inter_population/pca/pca_precomp.pdf")
plot(pc.precomp$x[,1],pc.precomp$x[,2],col=ind.colors,
     xlim = c(-0.6, 0.6), ylim = c(-0.15, 0.15),
     xlab="PC1 (16.40%)",ylab="PC2 (6.77%)",cex=1.5, pch=19)
legend("bottomright", legend = c("D1", "HP"), 
       col = c("#ff7f00", "#377eb8"), pch = 19, bty = "n", cex = 1)
abline(v=0,h=0,col="grey",lty=3)
dev.off()

#####
## PCA optmized for genlight object
#####
pca<-glPca(gen)

# PVE
pdf(file = "results/inter_population/pca/fraction_variance_pca_genlight.pdf")
plot(pca$eig/sum(pca$eig),ylab="Fraction of total variance", ylim=c(0,0.5), type='b')
dev.off()
round(100*(pca$eig[1]/(sum(pca$eig))), digits = 2) #16.40
round(100*(pca$eig[2]/(sum(pca$eig))), digits = 2) #6.82

# PCA plot
pdf(file = "results/inter_population/pca/pca_genlight.pdf")
par(xpd=FALSE,oma=c(0,0,0,0), pty="m")
plot(pca$scores[,1], pca$scores[,2], col=adjustcolor(ind.colors, alpha.f = 0.8),
     xlim = c(-4, 4), ylim = c(-3, 3),
     xlab="PC 1 (16.40%)", ylab="PC 2 (6.82%)", cex=2, pch=19, asp=1)
legend("bottomright", legend = c("Population D1", "Population HP"), 
       col = adjustcolor(c("#ff7f00", "#377eb8"), alpha.f = 0.8), pch = 19, bty = "n", cex = 0.8)
abline(v=0,h=0,col="grey",lty=3)
dev.off()

## Workflow to screen and find informative CANDIDATE SNPs
##-------------------------------------------------------------------------

## GOAL: identify informative candidate SNPs to include in a target panel
## by using a population genetics approach

## FIRST PART: Population genetics-based filtering:
# 1-With a PCA analysis, identify informative SNPs by ranking the PCA loadings;
# 2-Outliers detection of non-neutral loci;
# 3-Overall and population-specific HWP test;
# 4-Remove SNPs with 100% of heterozygotes genotypes;
# 5-Remove SNPs within or close to a clusters of SNPs;
# 6-Get the final list of candidate SNPs


## FIRST PART: Population genetics-based filtering:
# use all SNPs including the ones found in mtDNA 
# to reproduce original results for the pca loadings snp filtering

## STEP 1.0 - Principal component analysis
gen_original <- new("genlight", snps)
locNames(gen_original)<-paste("SNP", 1:nLoc(gen_original),sep="")
pca_original <- glPca(gen_original)

## STEP 1.1 - PCA loadings
loadings<-loadingplot(pca_original) #threshold of 75%

pca_loadings <- data.frame(Markers=colnames(snps)[loadings$var.idx], 
                           Loci=gen_original@loc.names[loadings$var.idx],
                           index=loadings$var.idx, 
                           pca_loading=loadings$var.values)

pca_loadings <- pca_loadings[order(-pca_loadings$pca_loading), ]

write.table(pca_loadings, 
            file="results/candidate_snps/pca_loadings.txt", 
            row.names = F,
            quote = F)

# Filtered loci lists based on the ranked PCA loadings
# Handling data to access loci/scaffold name
# FIRST LIST: Highest 60th PCA loadigns SNPs
candidate_snps_1 <- c(417,492,611,896,490,897,693,688,690,641,661,898,495,612,692,795,616,649,628,
                      650,1213,687,613,1010,643,1203,9,1051,701,1105,617,124,984,647,799,421,1044,
                      736,696,648,360,744,72,935,1067,893,1210,486,794,545,453,654,644,568,773,934,
                      1030,543,891,365)

candidate_snps_list_1 <- data.frame(Markers=colnames(snps)[candidate_snps_1], 
                                    Loci=gen_original@loc.names[candidate_snps_1])

write.table(candidate_snps_list_1, 
            file="results/candidate_snps/candidate_snps_list_1.txt",
            row.names = F,
            quote = F)

# SECOND LIST: Following 30th loci - 61 to 90 PCA loadings
candidate_snps_2 <- c(424,573,747,749,572,624,656,336,1121,566,707,708,651,933,1216,992,750,625,
                      487,646,396,739,280,535,548,727,575,729,782,932,1038)

candidate_snps_list_2 <- data.frame(Markers=colnames(snps)[candidate_snps_2], 
                                    Loci=gen_original@loc.names[candidate_snps_2])

write.table(candidate_snps_list_2, 
            file="results/candidate_snps/candidate_snps_list_2.txt",
            row.names = F,
            quote = F)

## STEP 2 - Outliers detection - BAYESCAN
# see src/wgs_filtered_bayescan_charts.R

## STEP 3 - Overall and population-specific HWP test
# src/wgs_filtered_hwetests_HardyWeinberg_vcftools.R

## STEP 4 - Removal of SNPs with 100% of heterozygotes genotypes - IGV and R
# Visual inspection made manually using IGV
# Statistical inspection with HO and HE calculation

## STEP 5 - Removal of clusters of SNPs - IGV and R
# Visual inspection made manually using IGV


## Comparative SNP set analysis: in silico VALIDATION of targeted SNPs
##-------------------------------------------------------------------------

# Creating the vectors with the candidate snps indexes
targeted_idx <- c(9,10,11,12,43,74,75,91,124,135,145,155,195,209,365,421,487,611,612,613,624,625,641,642,643,644,646,647,648,
                  649,650,654,655,656,673,687,688,693,696,695,773,795,799,896,897,932,933,934,992,1010,1017,1038,
                  1051,1067,1098,1110,1121,1203,1210)

# Subseting the original snp set with the candidate snps
targeted_snps <- snps[,targeted_idx]

# remove mtDNA SNP and a targeted that did not enter in the panel - JPYR01002157.1_336
targeted_snps <- targeted_snps[,-c(16,47)]

# Working with final fluidigm snp set - filtered loci, random and genes
targeted_samples.names   <- scan(file="data/datasets/targeted.samples.txt", what="character",quiet=TRUE)
rownames(targeted_snps) <- targeted_samples.names

# markers names
markers_names <- c(sprintf("%s0%d", "Bant_tgt", 1:9), sprintf("%s%d", "Bant_tgt", 10:57))
colnames(targeted_snps) <- markers_names

## Exploratory analysis - genotypes and allele frequencies
##-----------------------------------------------------------

###############################
#  Allele frequency spectrum  #
###############################

## manually calculated AFS
##########################

# D1 AFS plot - w/out monomorphic SNPs
pdf(file="results/intra_population/genotypes_afs/afs_d1_targeted.pdf")
barplot(calculateAFS(x=targeted_snps[1:11,-c(15,22,34,36,37,43,52,54,55)],n=11, genes=2),
        ylim = c(0,20), ylab = "Number of SNPs", xlab = "", names.arg=seq(1,11,1), col = "#ff7f00")
dev.off()

# HP AFS plot - w/out monomorphic SNPs
pdf(file="results/intra_population/genotypes_afs/afs_hp_targete.pdf")
barplot(calculateAFS(x=targeted_snps[12:23,-c(5,10,13,16,17,18,19,29,30,32,33,39,48,53,56)],n=12,genes=2),
        ylim = c(0,20),ylab = "Number of SNPs", xlab = "", names.arg=seq(1,12,1), col = "#377eb8")
dev.off()

# Global AFS plot
pdf(file="results/intra_population/genotypes_afs/afs_global_targeted.pdf")
barplot(calculateAFS(x=targeted_snps, n=23,genes=2), ylim = c(0,20),
        ylab = "Number of SNPs", xlab = "", names.arg=seq(1,23,1), col = "#999999")
dev.off()

## Data structure - distances, heatmaps and PCAs
##-----------------------------------------------

######################
# Euclidean Distance #
######################

dist_targeted<-vegdist(targeted_snps, method = "euclidean", na.rm=TRUE)

pdf(file = "results/inter_population/distance/nj_euclidian_ind_dist_targeted_snps.pdf")
nj.d_targeted<-nj(dist_targeted)
plot(nj.d_targeted)
dev.off()

pdf(file = "results/inter_population/distance/ward_euclidian_ind_dist_targeted_snps.pdf")
ward_targeted<-hclust(dist_targeted, method = "ward.D")
plot(ward_targeted)
dev.off()

###########
# Heatmap #
###########

pdf(file = "results/inter_population/distance/heatmap_targeted_snps.pdf")
heat.obj_targeted<-heatmap(as.matrix(targeted_snps))
dev.off()
#str(heat.obj)

#################################
# Principal componente analysis #
#################################

#####
## Price et al 2010 - vcov of matrix individuals x markers
#####

W.targeted<-scale(targeted_snps, scale=FALSE) #centering
W.targeted[1:10,1:10]
W.targeted[is.na(W.targeted)]<-0
cov.W.targeted<-cov(t(W.targeted))
eig.result.targeted<-eigen(cov.W.targeted)
eig.vec.targeted<-eig.result.targeted$vectors
lambda.targeted<-eig.result.targeted$values

# PVE
pdf(file = "results/inter_population/pca/fraction_variance_pca_price-et-al-2010_targeted.pdf")
plot(lambda.targeted/sum(lambda.targeted),ylab="Fraction of total variance", ylim=c(0,0.5), type='b')
dev.off()
round(100*lambda.targeted[1]/sum(lambda.targeted),2) #46.10
round(100*lambda.targeted[2]/sum(lambda.targeted),2) #9.06

# PCA plot
pdf(file = "results/inter_population/pca/pca_price-et-al-2010_targeted.pdf")
plot(eig.vec.targeted[,1],eig.vec.targeted[,2],col=ind.colors, 
     xlim = c(-0.4, 0.4), ylim = c(-0.4, 0.4),
     xlab="PC1 (46.10%)", ylab="PC2 (9.06%)", cex=1.5, pch=19)
legend("bottomright", legend = c("D1", "HP"), 
       col = c("#ff7f00", "#377eb8"), pch = 19, bty = "n", cex = 1)
abline(v=0,h=0,col="grey",lty=3)
dev.off()

#####
## precomp package - SVD decomposition of COV matrix (Price)
#####

pc.precomp_targeted<-prcomp(cov.W.targeted, center=F, scale=F, retx=TRUE)
names(pc.precomp_targeted)

lambda.2.targeted<-pc.precomp_targeted$sdev

# PVE
pdf(file = "results/inter_population/pca/fraction_variance_pca_precomp_targeted.pdf")
plot(lambda.2.targeted/sum(lambda.2.targeted),ylab="Fraction of total variance", ylim=c(0,0.5), type='b')
dev.off()
round(100*lambda.2.targeted[1]/sum(lambda.2.targeted),2) #46.10
round(100*lambda.2.targeted[2]/sum(lambda.2.targeted),2) #9.06

pdf(file = "results/inter_population/pca/pca_precomp_targeted.pdf")
plot(pc.precomp_targeted$x[,1],pc.precomp_targeted$x[,2],col=ind.colors,
     xlim = c(-0.6, 0.6), ylim = c(-0.15, 0.15),
     xlab="PC1 (46.10%)",ylab="PC2 (9.06%)",cex=1.5, pch=19)
legend("bottomright", legend = c("D1", "HP"), 
       col = c("#ff7f00", "#377eb8"), pch = 19, bty = "n", cex = 1)
abline(v=0,h=0,col="grey",lty=3)
dev.off()

#####
## PCA optmized for genlight object
#####

## Adegenet package
###################

# table to genlight
gen_targeted<-new("genlight", targeted_snps)

# Acessors
ploidy(gen_targeted) <- 2
nInd(gen_targeted)
nLoc(gen_targeted)
pop(gen_targeted)
other(gen_targeted)

indNames(gen_targeted)
#pop(gen_targeted)<-pop.ind
other(gen_targeted)<-ind.colors

pca_targeted<-glPca(gen_targeted)

# PVE
pdf(file = "results/inter_population/pca/fraction_variance_pca_genlight_targeted.pdf")
plot(pca_targeted$eig/sum(pca_targeted$eig),ylab="Fraction of total variance", ylim=c(0,0.5), type='b')
dev.off()
round(100*(pca_targeted$eig[1]/(sum(pca_targeted$eig))), digits = 2) #45.56
round(100*(pca_targeted$eig[2]/(sum(pca_targeted$eig))), digits = 2) #9.05

# PCA plot
pdf(file = "results/inter_population/pca/pca_genlight_targeted.pdf")
par(xpd=FALSE,oma=c(0,0,0,0), pty="m")
plot(pca_targeted$scores[,1], pca_targeted$scores[,2], col=adjustcolor(ind.colors, alpha.f = 0.8),
     xlim = c(-4, 4), ylim = c(-3, 3),
     xlab="PC 1 (45.56%)", ylab="PC 2 (9.05%)", cex=2, pch=19, asp=1)
legend("bottomright", legend = c("Population D1", "Population HP"), 
       col = adjustcolor(c("#ff7f00", "#377eb8"), alpha.f = 0.8), pch = 19, bty = "n", cex = 0.8)
abline(v=0,h=0,col="grey",lty=3)
dev.off()

## Prepare the egglib outout file for the targeted dataset
##-------------------------------------------------------
#targeted_samples   <- scan(file="data/datasets/targeted.samples.txt", what="character",quiet=TRUE)
#
#c_header           <- c("chrom", "position","alleles")
#
#targeted.data <- read.table(file="data/datasets/targeted.txt", 
#                            header=FALSE, 
#                            col.names = c(c_header,targeted_samples), 
#                            check.names = FALSE)
#
## remove mtDNA SNP and a targeted that did not enter in the panel - JPYR01002157.1_336
#targeted.data <- targeted.data[-c(16,47),]
#
## convert dataset to transposed 0,1,2 data - BUT KEEP SAMPLES IN COLUMNS
#targeted.data.mtx <- as.matrix(targeted.data[,-c(1:3)])
#
#targeted.data.mtx[targeted.data.mtx == "./."] <- NA
#targeted.data.mtx[targeted.data.mtx == "0/0"] <- 0
#targeted.data.mtx[targeted.data.mtx == "1/1"] <- 2
#targeted.data.mtx[targeted.data.mtx == "0/1" | targeted.data.mtx == "1/0"] <- 1
#
## first create the egglib summary stats input file
#targeted.data.mtx2egglib <- targeted.data.mtx
#
#targeted.data.mtx2egglib[targeted.data.mtx2egglib == "2"] <- "22"
#targeted.data.mtx2egglib[targeted.data.mtx2egglib == "0"] <- "11"
#targeted.data.mtx2egglib[targeted.data.mtx2egglib == "1"] <- "12"
#targeted.data.mtx2egglib[is.na(targeted.data.mtx2egglib)] <- "00"
#
#colnames(targeted.data.mtx2egglib) <- c(paste0("indiv", seq(from=1, to=11, by=1), "@pop1", ""), 
#                                        paste0("indiv", seq(from=1, to=12, by=1), "@pop2", ""))
#
#targeted.data.data2egglib <- as.data.frame(targeted.data.mtx2egglib)
#
#targeted2egglib <- data.frame(targeted.data[,c(1:2)], status="S",selection="Y", alleles=targeted.data[,c(3)])
#
#targeted2egglib <- cbind(targeted2egglib, targeted.data.data2egglib)
#
#write.table(targeted2egglib, file = "data/datasets/targeted_egglib_input.txt", quote=FALSE, sep="\t", row.names = FALSE)
#
## run egglib
#egglib_run <- paste("/Users/vitorpavinato/myPython2venvs/python-egglib3.0.0b22/bin/python",
#                    paste0(getwd(), "/", "src/pythonscripts/summstats.py"),
#                    paste0("input-file=", "data/datasets/targeted_egglib_input.txt"),
#                    paste0("output-file=", "results/summary_statistics/targeted_egglib_summstats.txt"),
#                    paste0("LSS=", paste0(c("He", "Dj", "WCst"), collapse = ",")),
#                    paste0("LSSp=", paste0(c("He", "Dj"), collapse = ",")),
#                    paste0("WSS=", paste0(c("He", "Dj", "WCst", "S", "thetaW", "Pi", "D", "Da", "ZZ", "ZnS"), collapse = ",")),
#                    paste0("WSSp=", paste0(c("He", "Dj", "S", "thetaW", "Pi", "D", "ZZ", "ZnS"), collapse = ",")),
#                    paste0("GSS=", paste0(c("He", "Dj", "WCst", "S", "thetaW" , "Pi", "D", "Da", "SFS"), collapse = ",")),
#                    paste0("GSSp=",paste0(c("He", "S", "thetaW", "Pi", "D", "Da"), collapse = ",")),
#                    paste0("wspan=", 150),
#                    paste0("SFS-bins=", 20),
#                    paste0("select=", "all"))
#
#system(egglib_run)
