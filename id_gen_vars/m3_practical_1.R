
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
  g_frame = geno_data['MSNP3']
  colnames(g_frame) = c(loc)
  Geno = genotype(g_frame$loc, sep='')
  print(summary(Geno))
}

test_frame = geno_data$MSNP3
test_frame1 = geno_data['MSNP3']#as.factor(as.vector(geno_data['MSNP3']))

Geno = genotype(geno_data$MSNP3, sep='')
summary(Geno)