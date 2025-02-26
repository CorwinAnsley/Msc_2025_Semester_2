#install.packages('qqman')

library(qqman)

# Load data
chrom1_linear = read.table("./cw_data/cw_output_qc__chr_1.assoc.linear" , header = T)
chrom2_linear = read.table("./cw_data/cw_output_qc__chr_2.assoc.linear" , header = T) 

chrom1_logistic = read.table("./cw_data/cw_output_dichot_qc__chr_1.assoc.logistic" , header = T)
chrom2_logistic = read.table("./cw_data/cw_output_dichot_qc__chr_2.assoc.logistic" , header = T) 
# Make manhattan plot
#plot(chrom$BP, -log10(chrom$P))

# Get only genome wide significant hita

#sig_results = chrom[which(chrom$P <= 5e-8),]

# Plot with qqman

# Logistic 1
png('./plots/manhattan_logistic_1.png')
manhattan(chrom1_logistic)
# make plot
dev.off()

png('./plots/qq_logistic_1.png')
qq(chrom1_logistic$P)
# make plot
dev.off()

# Logistic 2
png('./plots/manhattan_logistic_2.png')
manhattan(chrom2_logistic)
# make plot
dev.off()

png('./plots/qq_logistic_2.png')
qq(chrom2_logistic$P)
# make plot
dev.off()