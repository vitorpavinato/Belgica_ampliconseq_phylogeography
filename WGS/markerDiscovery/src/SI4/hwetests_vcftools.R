###########################################
#       vcftools HWP test and charts      #
###########################################

## Vitor Pavinato
## vitor.pavinato@supagro.fr
## INRA

rm(list=ls())
ls()

sessionInfo()$running
sessionInfo()$platform
R.version.string
.Platform$GUI

recordSessionInfo <- sessionInfo();
save(recordSessionInfo, file="results/intra_population/hwe/sessionInfo.RData");
save.image(file = "results/intra_population/hwe/workspaceFile.RData");

load(file = "results/intra_population/hwe/workspaceFile.RData")

# HWE exact test obtained with vcftools
hwe.global <- read.table(file = "results/intra_population/hwe/filtered_global_hwe.txt", header = T)
hwe.d1 <- read.table(file = "results/intra_population/hwe/filtered_d1_hwe.txt", header = T)
hwe.hp <- read.table(file = "results/intra_population/hwe/filtered_hp_hwe.txt", header = T)

hwe.global <- hwe.global[-c(417:431), ]
hwe.d1 <- hwe.d1[-c(417:431), ]
hwe.hp <- hwe.hp[-c(417:431), ]

## GLOBAL HWP test
hwe.data <- data.frame(Markers=paste0(hwe.global$CHR, ":", hwe.global$POS), pvalue=hwe.global$P_HWE)
summary(hwe.data$pval)

alpha=0.05
alpha.bc=0.05/1245

# without bonferroni correction
sum(hwe.data$pval < alpha) #710
prop.dp.pexact <- (sum(hwe.data$pvalue < alpha)/length(hwe.data$Markers))*100 # 57.02811

# with bonferroni correction
sum(hwe.data$pvalue < alpha.bc) #461
prop.dp.pexact.bc <- (sum(hwe.data$pvalue < alpha.bc)/length(hwe.data$Markers))*100 # 37.02811

# Population = D1
#-----------------

hwe.data.d1 <- data.frame(Markers=paste0(hwe.d1$CHR, ":", hwe.d1$POS), pvalue=hwe.d1$P_HWE)
summary(hwe.data.d1$pval)

# without bonferroni correction
sum(hwe.data.d1$pval < alpha) #598
prop.dp.pexact.d1 <- (sum(hwe.data.d1$pvalue < alpha)/length(hwe.data.d1$Markers))*100 # 48.03213

# with bonferroni correction
sum(hwe.data.d1$pvalue < alpha.bc) #0
prop.dp.pexact.d1.bc <- (sum(hwe.data.d1$pvalue < alpha.bc)/length(hwe.data.d1$Markers))*100 # 0

# Population = HP
#-----------------

hwe.data.hp <- data.frame(Markers=paste0(hwe.hp$CHR, ":", hwe.hp$POS), pvalue=hwe.hp$P_HWE)
summary(hwe.data.hp$pval)

# without bonferroni correction
sum(hwe.data.hp$pval < alpha) #629
prop.dp.pexact.hp <- (sum(hwe.data.hp$pvalue < alpha)/length(hwe.data.hp$Markers))*100 # 50.52209

# with bonferroni correction
sum(hwe.data.hp$pvalue < alpha.bc) #0
prop.dp.pexact.hp.bc <- (sum(hwe.data.hp$pvalue < alpha.bc)/length(hwe.data.hp$Markers))*100 # 0

# Proportion of PHW deviates - without bonferroni correction
prop <- c(prop.dp.pexact, prop.dp.pexact.d1, prop.dp.pexact.hp)
pops <- c("Global", "D1", "HP")

hwe.prop.table <- data.frame(pops=pops, prop=prop)

# Proportion of PHW deviates - with bonferroni correction
prop.bc <- c(prop.dp.pexact.bc, prop.dp.pexact.d1.bc, prop.dp.pexact.hp.bc)
pops.bc <- c("Global", "D1", "HP")

hwe.prop.table.bc <- data.frame(pops=pops.bc, prop=prop.bc)

# Combined Plot 2: Q-Q plots
pdf(file = "results/intra_population/hwe/vcftools_pexact_qq_plots.pdf", width = 6.5, height = 6)
par(mfrow=c(2,2))
plot(-log10((1:length(hwe.data$pval))/length(hwe.data$pval)), -log10(hwe.data$pval[order(hwe.data$pval)]),
     xlim = c(0,6), ylim = c(0,6),
     xlab="-log10(expected P)", ylab="-log10(observed P)", main="Q-Q plot - 'Global'")
abline(a=0, b=1, col="red")
mtext("(A)", side=3, line=1.5, cex=1.5, adj=-0.35)

plot(-log10((1:length(hwe.data.d1$pval))/length(hwe.data.d1$pval)), -log10(hwe.data.d1$pval[order(hwe.data.d1$pval)]),
     xlim = c(0,3), ylim = c(0,3),
     xlab="-log10(expected P)", ylab="-log10(observed P)", main="Q-Q plot - 'D1'")
abline(a=0, b=1, col="red")
mtext("(B)", side=3, line=1.5, cex=1.5, adj=-0.35)

plot(-log10((1:length(hwe.data.hp$pval))/length(hwe.data.hp$pval)), -log10(hwe.data.hp$pval[order(hwe.data.hp$pval)]),
     xlim = c(0,3), ylim = c(0,3),
     xlab="-log10(expected P)", ylab="-log10(observed P)", main="Q-Q plot - 'HP'")
abline(a=0, b=1, col="red")
mtext("(C)", side=3, line=1.5, cex=1.5, adj=-0.35)

barplot(hwe.prop.table$prop, names.arg=pops.bc, ylim = c(0,100),
        ylab = "% of departure from HWP", main = "% departure w/out correction")
mtext("(D)", side=3, line=1.5, cex=1.5, adj=-0.35)

dev.off()

save.image(file = "results/intra_population/hwe/workspaceFile.RData")



