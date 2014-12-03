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

source("http://www.uvm.edu/~rsingle/Rdata/scripts_stat295F14.R")
fms <- otherdata("FMS_data.tsv", sep="\t")
attach(fms)

#TEST CODE (to make sure the data is read in correctly):
summary(resistin_c30t)
# CC   CT NA's 
#712   21  664 
table(resistin_c30t) #NOTE: NA's not listed
# CC  CT 
#712  21 
#--------------------------------------------------------------------------------------

#We assign out cutoff at the 90th percentile
cutoff = quantile(NDRM.CH, c(.9), na.rm = TRUE) 

#testing with a random snp first.
chisq.test(table(NDRM.CH > cutoff, actn3_1671064))

#Take out the b2b col
fms2 = fms[,colnames(fms) != "b2b"]

#get the names of the snps
snps = names(fms2)[5:226]

#Initialize some lists
pvals = NULL
anomalies = NULL
for (snp in snps){
  pval = chisq.test(table(NDRM.CH > 100, fms2[[snp]]))$p.value
  
  #if there is something wierd about the pval let's exclude it. 
  if (pval < 1.205409e-10 || pval == "NaN"){ 
    anomalies = c(anomalies, snp)
    #print(snp)
  } else { 
    if ( (-log10(pval) > 2) & !(snp %in% anomalies) ){
      print(paste(snp, toString(pval)))
    }
    #if we didn't exclude it, put it into a nice list for plotting
    pvals = c(pvals, pval)
  }
}
#To the manhattans!
barplot(-log10(pvals))
text(150,3, "Here")
text(200,2.5, "Here2")
text(230,6.3, "Here3")
dev.off()

#Here is a new comment. 


