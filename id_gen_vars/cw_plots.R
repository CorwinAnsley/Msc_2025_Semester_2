install.packages('qqman')

library(qqman)

# Load data
chrom = read.table("./cw_data/cw_output_qc__chr_1.assoc.linear" , header = T) 

# Make manhattan plot
plot(chrom$BP, -log10(chrom$P))

# Get only genome wide significant hita

sig_results = chrom[which(chrom$P <= 5e-8),]

# Plot with qqman

manhattan(chrom)
qq(chrom$P)
