##########################################
##   SITE MEAN DEPTH, MAPQ & GQ QUAL    ##
##          VARIANT CALL FILES          ##
##      QUAL vs DEPTH RELATIONSHIP      ##
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

## RAW variants
##------------------------
raw.variants.dp <- read.table(file = "results/coverage/calledVariants/raw_mpileup_call.ldepth.mean",
                              header = TRUE)
dim(raw.variants.dp)
# 715721      4

# Mean across samples
mean_quantiles <- quantile(raw.variants.dp$MEAN_DEPTH, c(0.25, 0.75))
mean_iqr <- mean_quantiles[2] - mean_quantiles[1]
# 4.869561

# Mean across samples
pdf(file = "results/coverage/calledVariants/raw_mean_site_depth.pdf")
par(mar=c(5,5,4,1)+.1)
hist(t(raw.variants.dp$MEAN_DEPTH), breaks = 100,
     main = "",
     ylab = "Number of Variant Sites",
     xlab = "Site Depth",
     col = "#035aa6"
     #,
     #cex.axis = 1.2,
     #cex.lab = 1.2
     )
dev.off()

min(raw.variants.dp$MEAN_DEPTH) # 0.0434783
max(raw.variants.dp$MEAN_DEPTH) # 125
mean(raw.variants.dp$MEAN_DEPTH) # 4.170089

dim(raw.variants.dp[which(raw.variants.dp$MEAN_DEPTH > 20), ])
# 8127    4

## FILTERED variants
##------------------------
filtered.variants.dp <- read.table(file = "results/coverage/calledVariants/filtered.ldepth.mean",
                                   header = TRUE)
dim(filtered.variants.dp)
# 1260    4 

# Mean across samples
mean_quantiles.filtered <- quantile(filtered.variants.dp$MEAN_DEPTH, c(0.25, 0.75))
mean_filtered_iqr <- mean_quantiles.filtered[2] - mean_quantiles.filtered[1]
# 13.93475

# Mean across samples
pdf(file = "results/coverage/calledVariants/filtered_mean_site_depth.pdf")
par(mar=c(5,5,4,1)+.1)
hist(t(filtered.variants.dp$MEAN_DEPTH), breaks = 100, 
     main = "",
     ylab = "Number of Variant Sites",
     xlab = "Site Depth",
     col = "#035aa6"
     #,
     #cex.axis = 1.2,
     #cex.lab = 1.2
     )
dev.off()

min(filtered.variants.dp$MEAN_DEPTH) # 30
max(filtered.variants.dp$MEAN_DEPTH) # 79.2174
mean(filtered.variants.dp$MEAN_DEPTH) # 41.41094

## FILTERED TARGETED variants - including mtDNA sites
##---------------------------------------------------
filtered.targeted.variants.dp <- read.table(file = "results/coverage/calledVariants/targeted.ldepth.mean",
                                            header = TRUE) 
dim(filtered.targeted.variants.dp)
# 59  4

# Mean across samples
pdf(file = "results/coverage/calledVariants/targeted_mean_site_depth__.pdf")
hist(t(filtered.targeted.variants.dp$MEAN_DEPTH), breaks = 100, 
     main = "",
     ylab = "Number of Variant Sites",
     xlab = "Site Depth",
     col = "#035aa6",
     cex.axis = 1.5,
     cex.lab = 1.5)
dev.off()

##################
### SITE DEPTH ###
##################

## For each site, sum the depth across samples. 
## Results obtained with "vcftools --site-depth" or with bcftools query 'chr%CHROM\t%POS\t%REF,%ALT\t%DP\n'.
## This output file has the suffix ".DEPTH".

## RAW variants
##------------------------
raw.variants.totaldp  <- read.table(file = "results/coverage/calledVariants/raw_mpileup_call.vcf.DEPTH", header = F)
min(raw.variants.totaldp/23) # 0.08695652
max(raw.variants.totaldp/23) # 125.6522
mean(raw.variants.totaldp$V1/23) # 4.580311
mean(raw.variants.totaldp$V1) # 105.3472

## FILTERED variants
##------------------------
filtered.variants.totaldp  <- read.table(file = "results/coverage/calledVariants/filtered.vcf.DEPTH", header = F)
min(filtered.variants.totaldp/23) # 30.04348
max(filtered.variants.totaldp/23) # 88.95652
mean(filtered.variants.totaldp$V1/23) # 43.11798
mean(filtered.variants.totaldp$V1) # 991.7135

## FILTERED TARGETED variants - including mtDNA sites
##---------------------------------------------------
filtered.targeted.variants.totaldp  <- read.table(file = "results/coverage/calledVariants/targeted.vcf.DEPTH", header = F)
min(filtered.targeted.variants.totaldp/23) # 30.82609
max(filtered.targeted.variants.totaldp/23) # 80.30435
mean(filtered.targeted.variants.totaldp$V1/23) # 40.19381
mean(filtered.targeted.variants.totaldp$V1) # 924.4576

####################
### SITE MEAN GQ ###
####################

##  RAW variants
##------------------------
raw.variants.gq <- read.table(file = "results/variantCalling/gq/raw_variants_GQ.txt",
                              header = F)

pdf(file = "results/variantCalling/gq/raw_site_mean_gq.pdf")
par(mar=c(5,5,4,1)+.1)
hist(rowMeans(raw.variants.gq[,-c(1:3)]), breaks = 100,
     main = "",
     ylab = "Number of Variant Sites",
     xlab = "Mean GQ",
     col = "#ffa931",
     cex.axis = 1.5,
     cex.lab = 1.5)
dev.off()

## FILTERED variants
##------------------------
filtered.variants.gq <- read.table(file = "results/variantCalling/gq/filtered_GQ.txt",
                                   header = F)

pdf(file = "results/variantCalling/gq/filtered_site_mean_gq.pdf")
par(mar=c(5,5,4,1)+.1)
hist(rowMeans(filtered.variants.gq[,-c(1:3)]), breaks = 50,
     main = "",
     ylab = "Number of Variant Sites",
     xlab = "Mean GQ",
     col = "#ffa931",
     cex.axis = 1.5,
     cex.lab = 1.5)
dev.off()

## FILTERED TARGETD variants
##------------------------
targeted.variants.gq <- read.table(file = "results/variantCalling/gq/targeted_GQ.txt",
                                   header = F)

pdf(file = "results/variantCalling/gq/targeted_site_mean_gq.pdf")
par(mar=c(5,5,4,1)+.1)
hist(rowMeans(targeted.variants.gq[,-c(1:3)]), breaks = 20,
     main = "",
     ylab = "Number of Variant Sites",
     xlab = "Mean GQ",
     col = "#ffa931",
     cex.axis = 1.5,
     cex.lab = 1.5)
dev.off()

#################
### SITE QUAL ###
#################

##  RAW variants
##------------------------
raw.variants.QUAL <- apply(raw.variants.gq[, -c(1:3)], MARGIN = 1, sum)

min(raw.variants.QUAL) # 59
max(raw.variants.QUAL) # 2277
mean(raw.variants.QUAL) # 366.8173

## FILTERED variants
##------------------------
filtered.variants.QUAL <- apply(filtered.variants.gq[, -c(1:3)], MARGIN = 1, sum)

min(filtered.variants.QUAL) # 764
max(filtered.variants.QUAL) # 2270
mean(filtered.variants.QUAL) # 1516.018

## FILTERED TARGETD variants
##------------------------
targeted.variants.QUAL <- apply(targeted.variants.gq[, -c(1:3)], MARGIN = 1, sum)

min(targeted.variants.QUAL) # 826
max(targeted.variants.QUAL) # 2152
mean(targeted.variants.QUAL) # 1400.881

############
###  MQ  ###
############

## RAW variants
##------------------------
raw.variants.mq <- read.table(file = "results/variantCalling/mq/mq_raw_mpileup_call.txt",
                               header = F)

pdf(file = "results/variantCalling/mq/raw_mean_site_mq_.pdf")
hist(t(raw.variants.mq$V4), breaks = 20, 
     main = "",
     ylab = "Number of Variant Sites",
     xlab = "Site MQ",
     col = "#d92027",
     cex.axis = 1.5,
     cex.lab = 1.5)
dev.off()

raw.variants.mq.den <- density(raw.variants.mq$V4)

## FILTERED variants
##------------------------
filtered.variants.mq <- read.table(file = "results/variantCalling/mq/mq_filtered.txt",
                                     header = F)

pdf(file = "results/variantCalling/mq/filtered_mean_site_mq.pdf")
hist(t(filtered.variants.mq$V4), breaks = 20, 
     main = "",
     ylab = "Number of Variant Sites",
     xlab = "Site MQ",
     col = "#d92027",
     cex.axis = 1.5,
     cex.lab = 1.5)
dev.off()

filtered.variants.mq.den <- density(filtered.variants.mq$V4)

## FILTERED TARGETED variants
##------------------------
filtered.targeted.variants.mq <- read.table(file = "results/variantCalling/mq/mq_targeted.txt",
                                     header = F)

pdf(file = "results/variantCalling/mq/targeted_mean_site_mq.pdf")
hist(t(filtered.targeted.variants.mq$V4), breaks = 10, 
     main = "",
     ylab = "Number of Variant Sites",
     xlab = "Site MQ",
     col = "#d92027",
     cex.axis = 1.5,
     cex.lab = 1.5)
dev.off()

filtered.targeted.variants.mq.den <- density(filtered.targeted.variants.mq$V4)

########################
## PLOT DEPTH vs QUAL ##
##       FILTER 1     ##
########################

dpCutOff <- mean(raw.variants.totaldp$V1) + 3*sqrt(mean(raw.variants.totaldp$V1))
qualCutOff <- 2*mean(raw.variants.totaldp$V1)

# CORRECT relationship between DEPTH and QUAL
top.qualDP.color = ifelse(raw.variants.totaldp$V1 > dpCutOff & raw.variants.QUAL < 2*raw.variants.totaldp$V1,
                          "gray", "black")

# USED relationship between DEPTH and QUAL - High Depth & High Qual
top.qualDP.color2 = ifelse(raw.variants.totaldp$V1 > dpCutOff & raw.variants.QUAL > qualCutOff,
                           "black", "grey")

 
dpqual.table<- data.frame(DP=raw.variants.totaldp$V1, QUAL=raw.variants.QUAL, 
                          boleanDP=(raw.variants.totaldp$V1 > dpCutOff), boleanQUAL=(raw.variants.QUAL < 2*raw.variants.totaldp$V1),
                          boleanLowQ= (raw.variants.totaldp$V1 > dpCutOff) & (raw.variants.QUAL < 2*raw.variants.totaldp$V1), 
                          boleanLowQ2= !((raw.variants.totaldp$V1 > dpCutOff) & (raw.variants.QUAL > qualCutOff)),
                          boleanLowQcol=top.qualDP.color, boleanLowQcol2=top.qualDP.color2)



# DEPTH AND QUAL RELATIONSHIP
#----------------------------

# PLOT CORRECT RELATIONSHIP
par(mar=c(5,5,4,1)+.1)
plot(raw.variants.totaldp$V1, raw.variants.QUAL, 
     col= alpha(top.qualDP.color, 0.4), pch=19,
     xlab="DEPTH", ylab="QUAL", cex.axis = 1.2,
     cex.lab = 1.2)

points(filtered.variants.totaldp$V1, filtered.variants.QUAL, col="red")
points(filtered.targeted.variants.totaldp$V1, targeted.variants.QUAL, col=alpha("blue", 0.75), pch=19)
legend("bottomright", legend=c("Likely true variant", "Likely false variant", "Initially kept SNPs", "Targeted SNPs"), 
       col=c("black", "grey", "red", "blue"), pch=c(19,19,1,1), box.lty=0)

# PLOT USED RELATIONSHIP - High Depth & High Qual
par(mar=c(5,5,4,1)+.1)
plot(raw.variants.totaldp$V1, raw.variants.QUAL, 
     col= alpha(top.qualDP.color2, 0.4), pch=19,
     xlab="DEPTH", ylab="QUAL", cex.axis = 1.2,
     cex.lab = 1.2)
abline(v=dpCutOff, h=qualCutOff, col="darkgray", lty="dashed")

points(filtered.variants.totaldp$V1, filtered.variants.QUAL, col="red")
points(filtered.targeted.variants.totaldp$V1, targeted.variants.QUAL, col=alpha("blue", 0.75), pch=19)
legend("bottomright", legend=c("High DP & QUAL", "Low DP & QUAL", "Initially kept SNPs", "Targeted SNPs"), 
       col=c("black", "grey", "red", "blue"), pch=c(19,19,1,1), box.lty=0)

# MEAN SITE DEPTH PLOT CORRECTED AFTER REMOVAL OF LOW QUALITY SNPS
# MEAN SITE DEPTH PLOT WITH AND WITHOUT LOW QUALITY SNPS 
#-----------------------------------------------------------------

# PLOT CORRECT RELATIONSHIP
par(mar=c(5,5,4,1)+.1,  mfrow=c(1, 2))
hist(raw.variants.dp$MEAN_DEPTH, breaks = 100,
     main = "",
     ylab = "Number of Variant Sites",
     xlab = "Site Depth",
     col = "#035aa6")

hist(raw.variants.dp$MEAN_DEPTH[!dpqual.table$boleanLowQ], breaks = 100,
     main = "",
     ylab = "Number of Variant Sites",
     xlab = "Site Depth",
     col = "#035aa6")

write.table(x=raw.variants.dp[dpqual.table$boleanLowQ, c(1:2)], 
            file = "results/variantCalling/DepthQUAL_relationship/excludedHighDPlowQUALvar.txt",
            row.names = F, col.names = F, quote = F)

raw.variants.dp[!dpqual.table$boleanLowQ, c(1:2)]

# PLOT USED RELATIONSHIP - High Depth & High Qual
par(mar=c(5,5,4,1)+.1,  mfrow=c(1, 2))
hist(raw.variants.dp$MEAN_DEPTH, breaks = 100,
     main = "",
     ylab = "Number of Variant Sites",
     xlab = "Site Depth",
     col = "#035aa6")

hist(raw.variants.dp$MEAN_DEPTH[!dpqual.table$boleanLowQ2], breaks = 100,
     main = "",
     ylab = "Number of Variant Sites",
     xlab = "Site Depth",
     col = "#035aa6")


# MEAN SITE DEPTH HISTOGRAM WITH MAX MEAN DEPTH CUTOFF
# Filter applied after the removel of DP & QUAL SNPS
# CALLED FILTER 2
#---------------------------------------------------

raw.variants.dp.mod <- data.frame(ID=paste0(raw.variants.dp$CHROM,"_",raw.variants.dp$POS),
                                  MEAN_DEPTH=raw.variants.dp$MEAN_DEPTH,
                                  VAR_DEPT=raw.variants.dp$VAR_DEPTH)


filtered.variants.dp.mod <- data.frame(ID=paste0(filtered.variants.dp$CHROM,"_",filtered.variants.dp$POS),
                                       MEAN_DEPTH=filtered.variants.dp$MEAN_DEPTH,
                                       VAR_DEPT=filtered.variants.dp$VAR_DEPTH)


# CORRECT RELATIONSHIP
raw.variants.dp.mod.filtered <- raw.variants.dp.mod[!dpqual.table$boleanLowQ, ]

min(raw.variants.dp.mod.filtered[raw.variants.dp.mod.filtered$ID %in% filtered.variants.dp.mod$ID, "MEAN_DEPTH"]) # 30
max(raw.variants.dp.mod.filtered[raw.variants.dp.mod.filtered$ID %in% filtered.variants.dp.mod$ID, "MEAN_DEPTH"]) # 45.6087
mean(raw.variants.dp.mod.filtered[raw.variants.dp.mod.filtered$ID %in% filtered.variants.dp.mod$ID, "MEAN_DEPTH"]) # 33.61009

p1 <- hist(raw.variants.dp.mod.filtered$MEAN_DEPTH, breaks = 100)
p2 <- hist(raw.variants.dp.mod.filtered[raw.variants.dp.mod.filtered$ID %in% filtered.variants.dp.mod$ID, "MEAN_DEPTH"])

plot(p1, col="darkgray", xlim=c(0,80), ylim = c(0,15000), main = "", xlab = "Mean Site Depth")  # first histogram
plot(p2, col="red", xlim=c(0,80), add=T)  # second
abline(v=(dpCutOff/10), col="darkgreen", lty="dashed")
abline(v= 30, col="red", lty="dashed")
legend("topright", legend=c("Cutoff - DP cutoff/10", "min site depth of kept SNPs"), 
       col=c("darkgreen", "red"), lty=c("dashed", "dashed"), pch=NULL, box.lty=0)


# USED RELATIONSHIP
raw.variants.dp.mod.filtered2 <- raw.variants.dp.mod[!dpqual.table$boleanLowQ2, ]

min(filtered.variants.dp.mod$MEAN_DEPTH) # 30
max(filtered.variants.dp.mod$MEAN_DEPTH) # 79.2174
mean(filtered.variants.dp.mod$MEAN_DEPTH) # 41.41094

p1 <- hist(raw.variants.dp.mod.filtered2$MEAN_DEPTH, breaks = 100)
p2 <- hist(raw.variants.dp.mod.filtered2[raw.variants.dp.mod.filtered2$ID %in% filtered.variants.dp.mod$ID, "MEAN_DEPTH"])

plot(p1, col="darkgray", xlim=c(0,80), ylim = c(0,15000), main = "", xlab = "Mean Site Depth")  # first histogram
plot(p2, col="red", xlim=c(0,80), add=T)  # second
abline(v=(dpCutOff/10), col="darkgreen", lty="dashed")
abline(v= 30, col="red", lty="dashed")
legend("topright", legend=c("Cutoff - DP cutoff/10", "min site depth of kept SNPs"), 
       col=c("darkgreen", "red"), lty=c("dashed", "dashed"), pch=NULL, box.lty=0)
