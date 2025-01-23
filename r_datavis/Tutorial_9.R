library(reshape2)

em_scaled = read.table("./data/em_scaled.csv", header=TRUE, row.names=1, sep= "\t")
master= read.table("./data/master.csv", header=TRUE, row.names=1, sep= "\t")
master_sig = read.table("./data/master_sig.csv", header=TRUE, row.names=1, sep= "\t")
em_symbols_sig = read.table("./data/em_symbols_sig.csv", header=TRUE, row.names=1, sep= "\t")
em_scaled = na.omit(em_scaled)
sig_genes = row.names(master_sig)
em_scaled_sig = em_scaled[sig_genes,]
em_scaled_sig = na.omit(em_scaled_sig)

hm_matrix = as.matrix(em_scaled_sig[1:100,])
#hm_matrix = melt(hm_matrix)

source("./omics_functions.R")

ggp = plot_heatmap(em_scaled_sig)

ggp
