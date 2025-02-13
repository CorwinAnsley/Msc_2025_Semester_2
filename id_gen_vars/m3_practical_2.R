# Load data
chrom = read.table("./data/test_output_qc__1234567_1.assoc.linear" , header = T) 

# Make manhattan plot
plot(chrom$BP, -log10(chrom$P))

# Get only genome wide significant hita

sig_results = chrom[which(chrom$P <= 5e-8),]
