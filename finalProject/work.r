#--------------------------------------------------------------------------------------
# The following is R code written for exploratory analysis of the FAMUS dataset. 
# 
# Written by Nick Strayer and Peter Euclide for final project for Stats 295 Fall 2014 with
# Professor Richard Single at the University of Vermont.
#
# Code outputs csv files with results of each association test for use with an online
# manhattan plot generator. Hosted at www.uvm.edu/~nstrayer/manhattanPlot
#--------------------------------------------------------------------------------------
#READ IN THE FULL FAMUSS DATASET
setwd("/Users/Nick/fall14/statGen/finalProject")
library(vioplot)
library(haplo.stats)
library(ggplot2)
library(plyr)
source("http://www.uvm.edu/~rsingle/Rdata/scripts_stat295F14.R")
fms <- otherdata("FMS_data.tsv", sep="\t")
attach(fms)
#--------------------------------------------------------------------------------------

#We assign out cutoff at the 90th percentile
cutoff = quantile(NDRM.CH, c(.95), na.rm = TRUE) 

#testing with a random snp first.
chisq.test(table(NDRM.CH > cutoff, actn3_1671064))

#Take out the b2b col
fms2 = fms[,colnames(fms) != "b2b"]

#get the names of the snps
snps = names(fms2)[5:226]

percents_NA = NULL
for (snp in snps){
  percents_NA = c(percents_NA, ( summary(fms[[snp]])["NA's"] / length(fms[[snp]])) )
}

#hist(percent_NA, breaks = 20)
vioplot(percents_NA)

#So we will throw away snps with NA percentange of higher than 70. 


#-----------------------------------------------------------------------------------------------------
# Non-Dominant Arm association:
#----------------------------------------------------------------------------------------------------
#Initialize some lists
NDRM_pvals    = NULL
NDRM_usedSnps = NULL
for (snp in snps){
  percent_NA = summary(fms2[[snp]])["NA's"] / length(fms2[[snp]])
  tabledVals = table(NDRM.CH > 100, fms2[[snp]])
  if ((percent_NA < .7) && !(0 %in% tabledVals)){
    pval = chisq.test(tabledVals)$p.value
    NDRM_pvals = c(NDRM_pvals, pval)
    NDRM_usedSnps = c(NDRM_usedSnps, snp)
  }
}
bonferPVal = .05/length(NDRM_pvals)
bonferPVal
barplot(-log10(NDRM_pvals))

#make a csv for exporting
NDRM_Out = data.frame(NDRM_usedSnps, -log10(NDRM_pvals))
write.csv(NDRM_Out, file = "NDRM_Data.csv")

#To the manhattans!
pdf("manhattan.pdf")
barplot(-log10(pvals))
dev.off()

#Here is a new comment. 

#Find the pvalue for significance with conservative bonferoni correction. 
bonferPVal = .05/length(snps)


#-----------------------------------------------------------------------------------------------------
# Dominant arm association:
#----------------------------------------------------------------------------------------------------
DRM_cutoff = quantile(DRM.CH, c(.90), na.rm = TRUE) 
DRM_pvals    = NULL
DRM_usedSnps = NULL
for (snp in snps){
  percent_NA = summary(fms2[[snp]])["NA's"] / length(fms2[[snp]])
  tabledVals = table(DRM.CH > DRM_cutoff, fms2[[snp]])
  if ((percent_NA < .7) && !(0 %in% tabledVals)){
    pval = chisq.test(tabledVals)$p.value
    DRM_pvals = c(DRM_pvals, pval)
    DRM_usedSnps = c(DRM_usedSnps, snp)
  }
}
barplot(-log10(DRM_pvals))

#make a csv for exporting
DRM_Out = data.frame(DRM_usedSnps, -log10(DRM_pvals))
write.csv(DRM_Out, file = "d3Viz/DRM_Data.csv")

#-----------------------------------------------------------------------------------------------------
# Blood Pressure association:
#----------------------------------------------------------------------------------------------------
SBP_cutoff = quantile(SBP, c(.90), na.rm = TRUE) 
SBP_pvals    = NULL
SBP_usedSnps = NULL
for (snp in snps){
  percent_NA = summary(fms2[[snp]])["NA's"] / length(fms2[[snp]])
  tabledVals = table(SBP > SBP_cutoff, fms2[[snp]])
  if ((percent_NA < .7) && !(0 %in% tabledVals)){
    pval = chisq.test(tabledVals)$p.value
    SBP_pvals = c(SBP_pvals, pval)
    SBP_usedSnps = c(SBP_usedSnps, snp)
  }
}
barplot(-log10(SBP_pvals))

#make a csv for exporting
SBP_Out = data.frame(SBP_usedSnps, -log10(SBP_pvals))
write.csv(SBP_Out, file = "d3Viz/SBP_Data.csv")

#-----------------------------------------------------------------------------------------------------
# Hardy Weinberg Manhattan Plot:
#----------------------------------------------------------------------------------------------------
HWE_pvals = NULL
HWE_usedSnps = NULL
Color = NULL
for (snp in snps){
  currentSnp = fms2[[snp]]
  if (length(unique(currentSnp)) > 2){
    genotyped    = genotype(currentSnp, sep = "")
    HWE_pval     = HWE.chisq(genotyped)$p.value
    HWE_pvals    = c(HWE_pvals, HWE_pval)
    HWE_usedSnps = c(HWE_usedSnps, snp)
#     if ( (-log10(HWE_pval) < 4.5) && (-log10(HWE_pval) > 3.8)){ #some output to check on the interesting snps
#       print(paste(snp,HWE_pval))
#       print(summary(currentSnp))
#     }
  }
}

results = data.frame(HWE_usedSnps, -log10(HWE_pvals)) #Make a dataframe of the results for plotting with ggplot

#Let's plot with ggplot because we can, and it's pretty.
ggplot(data=results, aes(x=HWE_usedSnps, y=X.log10.HWE_pvals.)) + 
  geom_bar(aes(fill=X.log10.HWE_pvals.>3),stat="identity") +
  theme(axis.text.x = element_blank()) + 
  ggtitle("Hardy Weinberg Equilibrium chi-squared for SNPs") + 
  ylab("-Log10(P-Value)") + xlab("SNP")

#make a csv for exporting
HWE_Out = data.frame(HWE_usedSnps, -log10(HWE_pvals))
write.csv(HWE_Out, file = "HWE_Data.csv")


#-----------------------------------------------------------------------------------------------------
# OR for NDM most significant SNP:
#----------------------------------------------------------------------------------------------------
#run HWE chi-square 
table(nos3_rs1799983)
snp1 <- genotype(nos3_rs1799983, sep="")
summary(snp1)
HWE.chisq(snp1)
tabledVals = table(fms2$nos3_rs1799983,NDRM.CH > 100 )
chisq.test(tabledVals)

x <- tabledVals
or.GG.TT <- (x[1,1]*x[3,2])/(x[1,2]*x[3,1])
or.GT.TT <- (x[2,1]*x[3,2])/(x[2,2]*x[3,1])
or.GG.TT  # 3.95 more likely to have muscle bellow 90th percentile if GG
or.GT.TT  # 1.99 times more likely to have muscle below 90th percentile if GT
