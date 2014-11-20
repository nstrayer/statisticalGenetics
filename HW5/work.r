source("http://www.uvm.edu/~rsingle/Rdata/scripts_stat295F14.R")
amd <- otherdata("amd_chr1_fixed.dat") 
snps.all <- names(amd)[-c(1:3)]
head(amd[,1:13]) #the first 13 columns of the amd dataset
snps.all[1:10]
attach(amd)
table(amd$s2,amd$cc)

################################################################################################
# Problem 1
################################################################################################
snp2 = table(data.frame(cc,s2))[,-1]
print(chisq.test(snp2))

# > X-squared = 1.1978, df = 2, p-value = 0.5494
      
################################################################################################
# Problem 2
################################################################################################
pvals = c()
for (snp in snps.all){
  temp  = table(data.frame(cc,amd[snp]))                #table the snp with the case control status
  if(length(colnames(temp)) == 4){ temp = temp[,-1] }   #I am pretty sure that removing the 00s is the most frustrating coding thing I have ever done. 
  pvals = c(pvals, chisq.test(temp)$p.value )           #do a chi-squared test and get the p-value
}
summary(pvals)

################################################################################################
# Problem 3
################################################################################################
#pdf("manhattanPlot.pdf")
plot(-log10(pvals), xlab = "SNPs", main = "Manhattan plot of chi-squared p-values", col = "#1b9e77")
#dev.off()

################################################################################################
# Problem 4
################################################################################################
pdf("searchingHistogram.pdf") #lets make a histogram to try and find this artifact.
hist(pvals)
#hist(pvals, breaks = 400)
arrows(.6,500,.738,115)
text(.52,520, labels = "Ah ha! Got you!")
dev.off()


#We bring back a modified problem 2 loop to fish out these snps

suspiciousSNPs = c() #make a list of the snps
everyLength = c()
suspiciousLengths = c()        #Get their number of genotypes-- I added this after checking out these snps in the second loop
for (snp in snps.all){
  temp  = table(data.frame(cc,amd[snp]))   
  numOfCols = length(colnames(temp))
  everyLength = c(everyLength,numOfCols)
  if(numOfCols == 4){ temp = temp[,-1] }    
  pval = chisq.test(temp)$p.value
  if(pval > .735 && pval < .74){           #if they fall into my eyeballed range
    suspiciousSNPs = c(suspiciousSNPs, snp)#Add them to a list of snps I am suspicious of
    suspiciousLengths = c(suspiciousLengths, numOfCols)        #Lets grab their number of genotypes to be safe. 
    } 
}
summary(suspiciousLengths)
summary(everyLength)

for (snp in suspiciousSNPs){               #look over the list of snps I just made
  print(table(data.frame(cc,amd[snp]))  )             
}

################################################################################################
# Problem 5
################################################################################################
lowPvals = c(1,2)
lowSnps = c("blah", "blaah")

for (snp in snps.all){
  temp  = table(data.frame(cc,amd[snp]))   
  numOfCols = length(colnames(temp))
  if(numOfCols == 4){ temp = temp[,-1] }    
  pval = chisq.test(temp)$p.value
  if (pval < lowPvals[2]){ #if the pvalue is smaller than the larger of the two smallest seen pvals
    lowPvals[2] = pval #swap
    lowSnps[2]  = snp
    if(lowPvals[2] < lowPvals[1]){ #if out of order (small, bigger)
      lowPvals[2] = lowPvals[1] #swap the bigger one to 2nd space
      lowSnps[2]  = lowSnps[1]
      
      lowPvals[1] = pval #and put the new smaller one at the front
      lowSnps[1]  = snp
    }
  }
}
lowSnps #we're going to be looking for snps#s between 6000 and 8000
-log10(lowPvals) #and manhattaned pvals of >5

#this ran on the very first try, I am pretty proud of that. 

################################################################################################
# Problem 6
################################################################################################

#Allelic chi-squared
################################################################################################

#SNP = s7062
SNP = s7064

alleles = c()
caseCond = c()
for (i in 1:length(SNP)){
  genotype = toString(SNP[i])
  for (allele in strsplit(genotype, split = "")[[1]]) {
    alleles = c(alleles, allele)
    caseCond = c(caseCond, cc[i])
  }
}

# results_s7062 = table(alleles, caseCond)
# chisq.test(results_s7062)

results_s7064 = table(alleles, caseCond)[-1,]
chisq.test(results_s7064)

#ORs
################################################################################################
OR_Func    = function(a,b,c,d){ return ((a*d)/(b*c))}
OR_Conf_Int = function(a,b,c,d, OR){ 
  std_er = sqrt( (1/a) + (1/b) + (1/c) + (1/d) )
  bottomLn = log(OR) - 1.96 * std_er
  topLn    = log(OR) + 1.96 * std_er
  return (c(exp(bottomLn),exp(topLn)))}

#set up the tables...
s7062_table = table(s7062,cc)     
s7064_table = table(s7064,cc)[-1,]

#find the dominant OR 
#set up a-d for 7062 dominant
s7062_dom_a = (s7062_table[1,2]+ s7062_table[2,2])
s7062_dom_b = (s7062_table[1,1]+ s7062_table[2,1])
s7062_dom_c = s7062_table[3,2]
s7062_dom_d = s7062_table[3,1]

#and for 7064 dominant
s7064_dom_a = (s7064_table[1,2]+ s7064_table[2,2])
s7064_dom_b = (s7064_table[1,1]+ s7064_table[2,1])
s7064_dom_c = s7064_table[3,2]
s7064_dom_d = s7064_table[3,1]

#Find the OR and confidence int. 
s7062_dom_OR  = OR_Func(s7062_dom_a,s7062_dom_b, s7062_dom_c,s7062_dom_d)
s7062_dom_Conf_Int = OR_Conf_Int(s7062_dom_a,s7062_dom_b, s7062_dom_c,s7062_dom_d, s7062_dom_OR)

s7064_dom_OR  = OR_Func(s7064_dom_a,s7064_dom_b, s7064_dom_c,s7064_dom_d)
s7064_dom_Conf_Int = OR_Conf_Int(s7064_dom_a,s7064_dom_b, s7064_dom_c,s7064_dom_d, s7064_dom_OR)


#now find the reccessive OR

#for 7062
s7062_rec_a = s7062_table[1,2]
s7062_rec_b = s7062_table[1,1]
s7062_rec_c = (s7062_table[3,2] + s7062_table[2,2])
s7062_rec_d = (s7062_table[3,1] + s7062_table[2,1])

#for 7064
s7064_rec_a = s7064_table[1,2]
s7064_rec_b = s7064_table[1,1]
s7064_rec_c = (s7064_table[3,2] + s7064_table[2,2])
s7064_rec_d = (s7064_table[3,1] + s7064_table[2,1])

#Now to find the ORs and Conf Intervals:
s7062_rec_OR  = OR_Func(s7062_rec_a,s7062_rec_b, s7062_rec_c,s7062_rec_d)
s7062_rec_Conf_Int = OR_Conf_Int(s7062_rec_a,s7062_rec_b, s7062_rec_c,s7062_rec_d, s7062_rec_OR)

s7064_rec_OR  = OR_Func(s7064_rec_a,s7064_rec_b, s7064_rec_c,s7064_rec_d)
s7064_rec_Conf_Int = OR_Conf_Int(s7064_rec_a,s7064_rec_b, s7064_rec_c,s7064_rec_d, s7064_rec_OR)


#Now we do it logistic style
#-------------------------
# snp = s7062
# CC_Log = cc==2

snp = factor(s7064[s7064 != "00"])
CC_Log = cc[s7064 != "00"]==2

Dom <- c(1, 1, 0) #With A as the 'variant' allele
Rec <- c(1, 0, 0) #With A as the 'variant' allele

#Dominant model 
contrasts(snp,1) <- cbind(Dom)
contrasts(snp)
Dom_model <- glm(CC_Log ~ snp, family=binomial)
Dom_OR = exp(Dom_model$coeff)

#Recessive model 
contrasts(snp,1) <- cbind(Rec)
contrasts(snp)
Rec_model <- glm(CC_Log ~ snp, family=binomial)
Rec_OR = exp(Rec_model$coeff)


#95% CI for dominant OR
Dom_m.sum <- summary(Dom_model) 
ci.l <- exp( Dom_m.sum$coeff[,1] - 1.96*(Dom_m.sum$coeff[,2]) )
ci.u <- exp( Dom_m.sum$coeff[,1] + 1.96*(Dom_m.sum$coeff[,2]) )
Dom_CI = rbind(ci.l, ci.u)

#95% CI for recessive OR
Rec_m.sum <- summary(Rec_model) 
ci.l <- exp( Rec_m.sum$coeff[,1] - 1.96*(Rec_m.sum$coeff[,2]) )
ci.u <- exp( Rec_m.sum$coeff[,1] + 1.96*(Rec_m.sum$coeff[,2]) )
Rec_CI = rbind(ci.l, ci.u)

DomResults = paste("Dominant: OR:", toString(Dom_OR[2]), " Confidence Interval: (", toString(Dom_CI[1,2]), toString(Dom_CI[2,2]), ")")
RecResults = paste("Recesive: OR:", toString(Rec_OR[2]), " Confidence Interval: (", toString(Rec_CI[1,2]), toString(Rec_CI[2,2]), ")")
print(DomResults)
print(RecResults)

# Allelic OR
################################################################################################
#We take the allele tables from part a: 

# results_s7062
# s7062_Al_OR = OR_Func(results_s7062[2,1],results_s7062[2,2],results_s7062[1,1],results_s7062[1,2])

#4.032717

#now we uncomment the part a code for 7604 and repeat:

results_s7064
s7064_Al_OR = OR_Func(results_s7064[2,1],results_s7064[2,2],results_s7064[1,1],results_s7064[1,2])
#3.909774






