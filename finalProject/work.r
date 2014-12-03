#--------------------------------------------------------------------------------------
#VARIABLES
#
#Gender: this is sex, NOT GENDER
#Age
#Race
#NDRM.CH: % change in muscle mass before & after exercise training (Non-Dominant Arm)
#DRM.CH:  % change in muscle mass before & after exercise training (Dominant Arm)
#Pre.weight
#Pre.height
#pre.BMI
#SBP:     Systolic Blood Pressure
#Dom.Arm: Dominant Arm
#Post.weight
#Post.Height
#Calc.post.BMI
#
#SNPs in lots of different genes
#--------------------------------------------------------------------------------------
# Your assignment is to do an exploratory analysis of genetic and/or other factors  
# that are predictive of % change in muscle mass before & after exercise training 
# in the non-dominant arm (NDRM.CH). 

#--------------------------------------------------------------------------------------
#READ IN THE FULL FAMUSS DATASET

library(vioplot)

source("http://www.uvm.edu/~rsingle/Rdata/scripts_stat295F14.R")
fms <- otherdata("FMS_data.tsv", sep="\t")
attach(fms)
#--------------------------------------------------------------------------------------

#We assign out cutoff at the 90th percentile
cutoff = quantile(NDRM.CH, c(.9), na.rm = TRUE) 

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

#Initialize some lists
pvals = NULL
pval = NULL
numOfSnps = 0
for (snp in snps){
  
  percent_NA = summary(fms2[[snp]])["NA's"] / length(fms2[[snp]])
  tabledVals = table(NDRM.CH > 100, fms2[[snp]])
  if ((percent_NA < .7) && !(0 %in% tabledVals)){
    pval = chisq.test(tabledVals)$p.value
    pvals = c(pvals, pval)
    numOfSnps = numOfSnps + 1
    if (-log10(pval) > 2){
      print(snp)
    }
  } else {
  }
  }
bonferPVal = .05/numOfSnps
bonferPVal
barplot(-log10(pvals))


#To the manhattans!
pdf("labeled_manhattan.pdf")
barplot(-log10(pvals))
text(113,2.7,  "nos3_rs1799983", pos = 4, cex = 0.8)
text(167,2.2,  "resistin_c980g", pos = 4, cex = 0.8)
text(195,6.36, "visfatin_6947766", pos = 4, cex = 0.8)
dev.off()

#Here is a new comment. 

#Find the pvalue for significance with conservative bonferoni correction. 
bonferPVal = .05/length(snps)

sampleSnp = nos3_rs1799983

