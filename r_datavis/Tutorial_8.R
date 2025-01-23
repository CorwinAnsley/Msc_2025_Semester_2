ss = read.table("./data/sample_sheet.csv", header=TRUE, sep="\t")
master= read.table("./data/master.csv", header=TRUE, row.names=1, sep= "\t")

em_symbols = master[ , as.vector(ss$SAMPLE)]
em_scaled = data.frame(t(scale(data.frame(t(em_symbols)))))
em_scaled = na.omit(em_scaled)
  
  source("./omics_functions.R")
  
  #tnf_gene_data = get_gene_data('Tnf',em_symbols,ss$SAMPLE_GROUP,group_order=c('gut','duct','node'))
  
  candidate_genes = row.names(em_scaled)[1:10]
  ggp  = boxplot_facets(em_scaled, candidate_genes, ss$SAMPLE_GROUP, 2, 5, group_order=c('gut','duct','node'))
  ggp 

top_gene_data_m = melt(top_gene_data, id.vars='sample_group')

ggp = ggplot(tnf_gene_data,aes(x=sample_group,y=expression, fill=sample_group)) + geom_boxplot()
ggp
