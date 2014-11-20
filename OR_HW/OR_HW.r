#READ IN THE DATA
source("http://www.uvm.edu/~rsingle/Rdata/scripts_stat295F14.R")
fms <- otherdata("FMS_small.tsv", sep="\t")
attach(fms)
library(LDheatmap)
library(haplo.stats)

############################################################################################################
#Problem 2
############################################################################################################
#part a
actn3Snp = summary(genotype(actn3_rs540874, sep=""))
actn3Snp 

#part b
#Genotypes
tableActn   = table(actn3_rs540874)
proportions = tableActn/sum(tableActn)
genoProportions = cbind(tableActn, proportions)
genoProportions

##Allele 
genotypes = names(tableActn)
A = 0
G = 0
for (genotype in genotypes){
  for (allele in strsplit(genotype, split = "")[[1]]) {
    if (allele == "A"){
      A =  as.integer(tableActn[genotype]) + A
    } else {
      G =  as.integer(tableActn[genotype]) + G
    }
  }
}

alleles = c("A", "G")
counts  = c( A ,  G)
freqs   = counts/sum(counts)
alleleFreqs = cbind(alleles, counts, freqs)
alleleFreqs

############################################################################################################
#Problem 3
############################################################################################################
actn3.Summary     = summary(genotype(actn3_1671064, sep = ""))
actn3.MajorAllele = actn3.Summary$allele.freq
actn3.MajorAllele

races = c()
majorA = c()
for (race in unique(Race) ){
  if ( is.na(race) || race == "Am Indian" ){  #The american indian population had no data. 
  } else {
    tmpDf = fms[Race == race,]
    tmpActn3.Summary    = summary(genotype(tmpDf$actn3_1671064, sep = "") )
    tmpActn3.MajorAllele = tmpActn3.Summary$allele.freq[1,2]
    races  = c(races, race )
    majorA = c(majorA, tmpActn3.MajorAllele)
  }
}
MajorAlleleFreqs = cbind(races, majorA)
MajorAlleleFreqs 


############################################################################################################
#Problem 4
############################################################################################################
# 4) Test for a significant association (using a chi-squared test) between SNPs in the esr1 gene
# and pre-test BMI (use pre.BMI>25 to define the 2 BMI groups for the test).

#Get list of snps from the esr1 gene:

p_values = NULL
esr1_snps = c()
p_values  = c()
for (snp in names(fms)){ if (grepl("esr1", snp)){ esr1_snps = c(esr1_snps,snp) } }
esr1_snps

for (snp in esr1_snps){
  tmp = na.omit(  data.frame(fms[snp], (fms["pre.BMI"] > 25))   )
  X2Test = chisq.test(table(tmp))
  p_values = c(p_values, X2Test$p.value)
  #print(chisq.test(table(tmp))) 
}
results = data.frame(esr1_snps, p_values)
results


############################################################################################################
#Problem 5
############################################################################################################
# 5) For any in (4) that are significant, summarize the results. That is, describe which genotype
# categories are most responsible for the significant outcome.

#I am over coding this to safe guard future applications
significantSnps = results[results$p_values < 0.05,]$esr1_snps
for (snp in significantSnps){
  mod.anova <- aov(    na.omit(    (fms["pre.BMI"] > 25) ~ fms[[snp]]   )     )
  mod.anova
  print(TukeyHSD(mod.anova))
}
#So GG-AA are significantly different. 

############################################################################################################
#Problem 6
############################################################################################################
# 6) For any in (4) that are significant, compute the two genotype odds ratios relating pre-test 
# BMI (dichotomized as above) and the SNP using either GG or TT as the reference group. State an 
# interpretation of your computed odds ratios. 

#attach(fms)

x <- table(esr1_rs1042717,pre.BMI>25)
x
or.AA.GG <- (x[1,1]*x[3,2])/(x[1,2]*x[3,1])
or.GA.GG <- (x[2,1]*x[3,2])/(x[2,2]*x[3,1])
or.AA.GG  
or.GA.GG

#GA is higher associated with a high BMI, but both are smaller than GG (i.e. < 1)

############################################################################################################
#Problem 7
############################################################################################################
# 7) For any in (4) that are significant, compute a T-test comparing the pre.BMI scores (not
# dichotomized, but actual scores) in the two homozygote categories. Summarize the test results.
# Include histograms and boxplots [using the boxplot() command] of the two groups to assess
# whether or not a T-test is appropriate on these data. Perform non-parametric tests if
# necessary.

#So I am just going to go with the fact I am working with esr1_rs1042717 at this point. 

AAs = pre.BMI[esr1_rs1042717 == "AA"]
GGs = pre.BMI[esr1_rs1042717 == "GG"]

par(mfrow=c(2,1))
hist(AAs)
hist(GGs)
boxplot(AAs, GGs, names = c("AA" ,"GG") )
t.test(AAs,GGs)

wilcox.test(AAs,GGs) 

summary(GGs)
############################################################################################################
#Problem 8
############################################################################################################
# Summarize the LD patterns based on r-squared and D-prime for SNPs within the esr1 gene.
# Include a matrix of LD values and a plot of each [using the LDheatmap() function]. Are there
# any notable differences between results for the two measures? If so, describe them. Which
# pairs of SNPs in the esr1 gene have significant LD? Do all of these correspond to high values
# of D-prime and/or r-squared?

rs1801132 = genotype( esr1_rs1801132 , sep = "")
rs1042717 = genotype( esr1_rs1042717 , sep = "")
rs2228480 = genotype( esr1_rs2228480 , sep = "")
rs2077647 = genotype( esr1_rs2077647 , sep = "")
rs9340799 = genotype( esr1_rs9340799 , sep = "")
rs2234693 = genotype( esr1_rs2234693 , sep = "")
snp.data  = data.frame(rs1801132,rs1042717,rs2228480,rs2077647,rs9340799,rs2234693)

#Let's get some tables of LD values!
LD(snp.data)$"R^2"
LD(snp.data)$"D'"
LD(snp.data)$"P-value"
#rs2228480 to all other SNPs are much higher for d' than for R2
#In addition the little botton triangle of rs2077647-(rs9340799, rs2234693) and rs9340799 - rs2234693 are different as well. 
#These all correspond directly with the snp pairs with significant LD

#Now let's draw some LD heatmaps
png(filename = "r2.png") #First for R-Squared
LDheatmap(snp.data, LDmeasure="r",  add.map=F, SNP.name = substr(esr1_snps,6,14))
dev.off()

png(filename = "Dprime.png") #...and for D-prime
LDheatmap(snp.data, LDmeasure="D'",  add.map=F, SNP.name = substr(esr1_snps,6,14))
dev.off()

############################################################################################################
#Problem 9
############################################################################################################
# 9) Determine if there is significant deviation from HWP for the akt1_t10726c_t12868c SNP based
# on the full dataset. Stratifying your analysis by Race and interpret your results.

akt1.geno = genotype(akt1_t10726c_t12868c,sep = "")
summary(akt1.geno)
HWE.chisq(akt1.geno)

HWE.PVals = c()
for (race in races){
  tempDf = fms[Race == race,]
  temp.geno = genotype(tempDf$akt1_t10726c_t12868c,sep = "")
  HWE.PVals = c(HWE.PVals, HWE.chisq(temp.geno)$p.value)
}
HWE.results = cbind(races, HWE.PVals)



