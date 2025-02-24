
# If firewall bypass needed
#Sys.setenv(http_proxy='http://wwwcache.gla.ac.uk:8080')

#install.packages("genetics") 

#install.packages("HardyWeinberg")

library(genetics) 

geno_data = read.table('./data/module3-HWELDexercisegenotypedata-v3.txt', header = TRUE, stringsAsFactors = T)


GenoCount = summary(geno_data['MSNP1'])
GenoCount


for (col in colnames(geno_data)) {
  GenoCount = summary(geno_data[col])
  print(GenoCount)
}

for (col in colnames(geno_data)) {
  g_frame = geno_data[col]
  colnames(g_frame) = c('loc')
  Geno = genotype(g_frame$loc, sep='')
  print(summary(Geno))
}



get_nonhwe_snps = function(geno_data){
  snps = colnames(geno_data)
  nonhwe_snps = vector('list',length(snps))
  i = 1
  while (i <= length(snps)) {
    print(i)
    col = snps[i]
    print(col)
    
    g_frame = geno_data[col]
    colnames(g_frame) = c('loc')
    Geno = genotype(g_frame$loc, sep='')
    hwe_results = HWE.exact(Geno)
 i + 1
  }
  return(nonhwe_snps)
}

nonhwe_snps = get_nonhwe_snps(geno_data)



#test_frame = geno_data$MSNP2
test_table = table(geno_data['MSNP2'])#as.factor(as.vector(geno_data['MSNP3']))

Geno = genotype(geno_data$MSNP3, sep='')
summary(Geno)

MSNP2geno = genotype(geno_data$MSNP2, sep="") 

HWE.chisq(MSNP2geno) 

hwe_results = HWE.exact(MSNP3geno)
hwe_results$p.value

obs_counts = table(geno_data['MSNP3'])
obs_counts
#obs_counts = table(c('ct','tt','ct','tt','ct','ct'))
length(obs_counts)
names(obs_counts)[1]

obs_counts[[1]]
obs_counts[[2]]

n_obs = sum(obs_counts)

freq_ma = (2*obs_counts[[1]] + obs_counts[[2]])/(2*n_obs)

#(2*ObsCount[[1]] + ObsCount[[2]])/ (2*Nobs)

ExpCount = c(n_obs*(freq_ma)^2, 2*n_obs*freq_ma*(1-freq_ma), n_obs*(1-freq_ma)^2) 
ExpCount

ChiSqStat <- sum((obs_counts - ExpCount)^2/ExpCount) 
ChiSqStat 

SigThresholdStat <- qchisq(1-0.05,df=1) 
SigThresholdStat 