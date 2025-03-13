if (!require("BiocManager", quietly = TRUE)) 
  
  install.packages("BiocManager") 



BiocManager::install("DESeq2") 

library(DESeq2)
library(ggplot2)

countTable = read.csv("./data/gene_count_matrix_nd.csv",row.names=1)
colnames(countTable) = c('tb1','tb2','tb3','wf1','wf2','wf3')

colTable = read.csv("./data/Design-Exp.csv",row.names=1)

dds = DESeqDataSetFromMatrix(countData=countTable,colData=colTable,design= ~ group)
dds

# Task C - plots

#ggp = ggplot(data=countTable, aes(x=log2fold, y=-log10(p.adj), col=diffexpr, label=symbol)) + 
#  geom_point()

