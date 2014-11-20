source("http://www.uvm.edu/~rsingle/Rdata/scripts_stat295F14.R")
library(haplo.stats)
hgdp <- otherdata("HGDP_AKT1.txt") 
pops.use <- c("Bantu","Biaka_Pygmies","Yoruba","Palestinian","Druze","Bedouin","Sardinian","French_Basque","French","Colombian","Karitiana")
hgdp <- hgdp[hgdp$Population %in% pops.use,] #use just a subset of populations
hgdp$Geographic.area <- factor(hgdp$Geographic.area)
hgdp$Population      <- factor(hgdp$Population)
head(hgdp)
attach(hgdp)

biaka = hgdp[Population == "Biaka_Pygmies",] #Isolate just the Baika_Pygmies for problems 1-5
snp2 = biaka$AKT1.C6024T #Put the snps into easier to type names
snp2 = gsub("C","A", gsub("T","a", snp2))
snp3 = biaka$AKT1.G2347T
snp3 = gsub("G","B", gsub("T","b", snp3))

#Problem 1
library(gmodels)
CrossTable(snp2, snp3)
#4 are predetermined and none need to be sorted out by em as we are looking at only 2 points with 2 alleles each.

#Problem 2
library(genetics)
snp2.geno = genotype(snp2,sep="") #Run the genotype function on the snps
snp3.geno = genotype(snp3,sep="")

p.A = summary(snp2.geno)$allele.freq[1,2] #Let's grab the allele frequencies and name them. 
p.a = summary(snp2.geno)$allele.freq[2,2]
p.B = summary(snp3.geno)$allele.freq[1,2]
p.b = summary(snp3.geno)$allele.freq[2,2]
p.A #let's see 'em.
p.a
p.B
p.b

#Problem 3
geno           = cbind(substr(snp2,1,1), substr(snp2,2,2),substr(snp3,1,1), substr(snp3,2,2))
haplo          = haplo.em(geno, locus.label=c("snp2", "snp3"))
haplo.RowNum   = dim(haplo$haplotype)[1] #not sure why I care about making these this form. 
colNames       = c("haplotype", "frequency")
haplo.HapNames = c()
for (i in 1:haplo.RowNum){
  haplo.HapNames = c(haplo.HapNames, paste(haplo$haplotype[i,], collapse = " "))
}
haplo.HapNameFreq           = data.frame(haplo.HapNames, haplo$hap.prob) #assemble the dataframe
colnames(haplo.HapNameFreq) = colNames                                   #Name the columns easier to understand things
haplo.HapNameFreq

#Problem 4
D      = haplo.HapNameFreq$frequency[1] - (p.a * p.b)
Dprime = D/min(p.A*p.b, p.a*p.B)
r2     = D^2/ (p.A*p.a*p.B*p.b)
prob4 = c(D,Dprime,r2)

#Problem 5 
snps = data.frame(snp2.geno,snp3.geno)
LD   = LD(snps)
prob5 = c(LD$D[3],LD$"D'"[3],LD$"R^2"[3] )

#let's compare!
values = c("D", "D'", "r^2")
compareDf = data.frame(values, prob4, prob5) 
compareDf
# looks pretty good to me! 

#Problem 6
regions   = unique(Geographic.area) #let's initialize a bunch of stuff
SNP2      = NULL 
SNP3      = NULL
tempDf    = NULL
region    = NULL
population= NULL
snp2.MAF  = NULL
snp3.MAF  = NULL
HWE.pval2 = NULL
HWE.pval3 = NULL
Dprime    = NULL
rSquared  = NULL
tempSnps  = NULL

#Release the for loops! 
for (reg in regions){                                                #start by looping over the regions
  populations  = unique(hgdp$Population[Geographic.area == reg])     #get the populations in that region
  tempDf       = hgdp[Geographic.area == reg,]                       #set up a dataframe of just the region
  for (pop in populations){                                          #now loop over the present populations
    region     = c(region, reg)
    population = c(population, pop)
    SNP2.geno  = genotype(tempDf$AKT1.C6024T[tempDf$Population == pop], sep = "") #made capital to avoid messing with the snp variables from earlier problems
    SNP3.geno  = genotype(tempDf$AKT1.G2347T[tempDf$Population == pop], sep = "")
    snp2.MAF   = c(snp2.MAF, summary(SNP2.geno)$allele.freq[2,2])
    snp3.MAF   = c(snp3.MAF, summary(SNP3.geno)$allele.freq[2,2])
    HWE.pval2  = c(HWE.pval2, HWE.chisq(SNP2.geno)$p.value)
    HWE.pval3  = c(HWE.pval3, HWE.chisq(SNP3.geno)$p.value)
    tempSnps   = data.frame(SNP2.geno,SNP3.geno)
    Dprime     = c(Dprime,  LD(tempSnps)$"D'"[3])
    rSquared   = c(rSquared,LD(tempSnps)$"R^2"[3] )
  }
}

bigDf = data.frame(region, population, snp2.MAF, snp3.MAF, HWE.pval2, HWE.pval3, Dprime, rSquared)
bigDf