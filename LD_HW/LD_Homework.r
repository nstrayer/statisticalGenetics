source("http://www.uvm.edu/~rsingle/Rdata/scripts_stat295F14.R")
library(LDheatmap)
library(haplo.stats)
hgdp <- otherdata("HGDP_AKT1.txt") 
pops.use <- c("Bantu","Biaka_Pygmies","Yoruba","Palestinian","Druze","Bedouin",
              "Sardinian","French_Basque","French","Colombian","Karitiana")
hgdp <- hgdp[hgdp$Population %in% pops.use,] #use just a subset of populations
hgdp$Geographic.area <- factor(hgdp$Geographic.area)
hgdp$Population      <- factor(hgdp$Population)
table(hgdp$Population,hgdp$Geographic.area)
attach(hgdp)

#---------------------------------------------------------------------------------
# 1) Compute 3-locus haplotype frequency estimates 
# (snp2-snp3-snp4 = AKT1.C6024T - AKT1.G2347T - AKT1.G2375A) for each of the 3 
# Central_Africa populations in the reduced dataset. Describe to what extent the 
# most common haplotypes in each of the three populations are shared. 

Populations <- unique(Population[Geographic.area == "Central_Africa"]) #let's isolate the populations
Populations

AKT1snps <- names(hgdp)[substr(names(hgdp),1,4)=="AKT1"] #taken from EM-intro.R notes

geno <- cbind(substr(AKT1.C6024T,1,1), substr(AKT1.C6024T,2,2),
              substr(AKT1.G2347T,1,1), substr(AKT1.G2347T,2,2),
              substr(AKT1.G2375A,1,1), substr(AKT1.G2375A,2,2))
geno

AKT1snps.2.4 <- AKT1snps[2:4] #We want to select the 2nd through 4th snps

geno.R1 <- geno[Population=="Biaka_Pygmies",] #Select the three Central African populations per the Populations call earlier
geno.R2 <- geno[Population=="Bantu",]
geno.R3 <- geno[Population=="Yoruba",]

haplo.R1 <- haplo.em(geno.R1, locus.label=AKT1snps.2.4) #Run the haplotype function on the populations. 
haplo.R2 <- haplo.em(geno.R2, locus.label=AKT1snps.2.4)
haplo.R3 <- haplo.em(geno.R3, locus.label=AKT1snps.2.4)

haplo.R1$hap.prob #Let's see what we got. 
haplo.R2$hap.prob
haplo.R3$hap.prob

#What we will do now is, see how many identified haplotypes there are, then construct a zipped list of the haplotype name (allele1,allele2,allele3) and 
#also its corresponding frequency. From there we will plot in a comparison barchart to allow rapid understanding of the frequency differences. 

#First we get the number of haplotypes found for each group
haplo.R1.RowNum <- dim(haplo.R1$haplotype)[1]
haplo.R2.RowNum <- dim(haplo.R2$haplotype)[1]
haplo.R3.RowNum <- dim(haplo.R3$haplotype)[1]

#Set column names for later dataframe creation:
colNames <- c("haplotype", "frequency")

#Now we fill a list with the names of the haplotypes
haplo.R1.HapNames <- c()
for (i in 1:haplo.R1.RowNum){
  haplo.R1.HapNames <- c(haplo.R1.HapNames, paste(haplo.R1$haplotype[i,], collapse = " "))
}

#bind that list of names with their corresponding frequencies
haplo.R1.HapNameFreq <- data.frame(haplo.R1.HapNames, haplo.R1$hap.prob) #assemble the dataframe
colnames(haplo.R1.HapNameFreq) <- colNames                               #Name the columns easier to undestand things
haplo.R1.HapNameFreq<- haplo.R1.HapNameFreq[order(haplo.R1.HapNameFreq$haplotype),] #Sort them alphabetically for later comparison

#rinse and repeate for the other two groups. 
haplo.R2.HapNames <- c()
for (i in 1:haplo.R2.RowNum){
  haplo.R2.HapNames <- c(haplo.R2.HapNames, paste(haplo.R2$haplotype[i,], collapse = " "))
}
haplo.R2.HapNameFreq <- data.frame(haplo.R2.HapNames, haplo.R2$hap.prob) 
colnames(haplo.R2.HapNameFreq) <- colNames
haplo.R2.HapNameFreq<- haplo.R2.HapNameFreq[order(haplo.R2.HapNameFreq$haplotype),] 

haplo.R3.HapNames <- c()
for (i in 1:haplo.R3.RowNum){
  haplo.R3.HapNames <- c(haplo.R3.HapNames, paste(haplo.R3$haplotype[i,], collapse = " "))
}
haplo.R3.HapNameFreq <- data.frame(haplo.R3.HapNames, haplo.R3$hap.prob)
colnames(haplo.R3.HapNameFreq) <- colNames
haplo.R3.HapNameFreq<- haplo.R3.HapNameFreq[order(haplo.R3.HapNameFreq$haplotype),]

#Let's plot!
png(filename = "Haplotype_Frequency_Chart.png") # print the file to the home directory
attach(mtcars)
barplot(haplo.R1.HapNameFreq$frequency,names.arg = haplo.R1.HapNameFreq$haplotype, main = "Biaka Pygmies")
barplot(haplo.R2.HapNameFreq$frequency,names.arg = haplo.R2.HapNameFreq$haplotype, main = "Bantu")
barplot(haplo.R3.HapNameFreq$frequency,names.arg = haplo.R3.HapNameFreq$haplotype, main = "Yoruba")
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# 2) Create an LDheatmap() plot based on r^2 for each of the 4 geographic regions 
# in the reduced dataset. Copy and paste these into a text editor so that they all 
# fit on one page. 

GeoRegions = unique(Geographic.area)
snp.label <- c("snp1","snp2","snp3","snp4" )
#Let's plot!
for (region in GeoRegions){
  tempDf   <- hgdp[Geographic.area == region,]
  snp1     <- genotype(tempDf$AKT1.C0756A,sep="")
  snp2     <- genotype(tempDf$AKT1.C6024T,sep="")
  snp3     <- genotype(tempDf$AKT1.G2347T,sep="")
  snp4     <- genotype(tempDf$AKT1.G2375A,sep="")
  snp.data <- data.frame(snp1,snp2,snp3,snp4)
  png(filename = paste(region, ".png", sep = ""))
  LDheatmap(snp.data, LDmeasure="r",  add.map=F, SNP.name=snp.label)
  dev.off()
} 

#------------------------------------------------------------------------------------------------------------------------
#3) Compute estimates of pairwise LD, separately for each geographic region, using 
#the D' measure for all SNPs in the AKT1 gene. Cut and paste the LD matrices for 
#each region into a text editor. To make things easier, use all data from 
#individuals in a given geographic region together for this question. 

GeoRegions = unique(Geographic.area)
snp.label <- c("snp1","snp2","snp3","snp4" )
LDs <- NULL #For problem 6

for (region in GeoRegions){
  tempDf   <- hgdp[Geographic.area == region,]
  snp1     <- genotype(tempDf$AKT1.C0756A,sep="")
  snp2     <- genotype(tempDf$AKT1.C6024T,sep="")
  snp3     <- genotype(tempDf$AKT1.G2347T,sep="")
  snp4     <- genotype(tempDf$AKT1.G2375A,sep="")
  snp.data <- data.frame(snp1,snp2,snp3,snp4)
  print(region)
  #print(LD(snp.data)$"D'") #problem 3
  print(LD(snp.data)$"r")   #problem 4
  
#   avgLD <- mean( LD(snp.data) $"r" , na.rm=T) # Problem 6 get the avg ld for the matrix
#   print(avgLD ) 
#   LDs <- c(LDs, avgLD)
} 

print("The average linkage disequilibrium over all groups")
print(mean(LDs))
print(t.test(LDs)$conf.int)

#4) Repeat question #3 using the r correlation measure of LD. 
#See three

#5) Describe the extent to which the LD measures are similar or different across 
#geographic regions. Which measure highlights any  differences that may exist to 
#a greater extent? Can you describe any general trend. Does this trend hold for 
#all NSP pairs, or only some?

#6) Do the estimates of LD within a geographic region vary between populations in 
#that region? Do they vary to roughly the same extent in all regions, or only in 
#specific regions? Use the r correlation measure to answer this question. You will 
#need to compute LD measures separately for each population in order to answer this 
#question. 

GeoRegions = unique(Geographic.area)
snp.label <- c("snp1","snp2","snp3","snp4" )
#region <- "Israel"
#pop <- "Druze"
LDs <- NULL
for (region in GeoRegions){
  tempDf   <- hgdp[Geographic.area == region,]
  Populations = unique(tempDf$Population)
  
  for (pop in Populations){
    subTempDf <- tempDf[tempDf$Population == pop,]
    snp1      <- genotype(subTempDf$AKT1.C0756A,sep="")
    snp2      <- genotype(subTempDf$AKT1.C6024T,sep="")
    snp3      <- genotype(subTempDf$AKT1.G2347T,sep="")
    snp4      <- genotype(subTempDf$AKT1.G2375A,sep="")
    snp.data  <- data.frame(snp1,snp2,snp3,snp4)
    print(paste(pop, region, sep=" in ")) #allow for easy reading
    avgLD     <- mean( LD(snp.data) $"r" , na.rm=T) # get the avg ld for the matrix
    print(avgLD)
    LDs       <- c(LDs, avgLD) #append to an array of average LDs for end meaning.
  }
} 

print("The average linkage disequilibrium over all groups")
print(mean(LDs))
print(t.test(LDs)$conf.int)


#7) use the results of the 3-locus haplotype frequency estimates from question #1 
##to help explain any patterns seen in the answer to question #6.