if (!require("BiocManager", quietly = TRUE)) 
  
  install.packages("BiocManager") 



BiocManager::install("DESeq2") 

library(DESeq2)
library(ggplot2)

BiocManager::install("limma")
library("limma")

countTable = read.csv("./data/gene_count_matrix_nd.csv",row.names=1)
colnames(countTable) = c('tb1','tb2','tb3','wf1','wf2','wf3')

colTable = read.csv("./data/Design-Exp.csv",row.names=1)

dds = DESeqDataSetFromMatrix(countData=countTable,colData=colTable,design= ~ group)
dds

# Filter out zeroes
notAllZero = (rowSums(counts(dds))>0)
dds = dds[notAllZero,]


dim(dds)[1]

dds = DESeq(dds)

# Task C - plots

plotDispEsts(dds)

# Perform rlog for qc
rld = rlog(dds, blind=TRUE)

# matrix of log2(raw-counts)
lgc.raw = log2(counts(dds,normalized=FALSE)+1)

# matrix of log2(normalized-counts)
lgc.norm = log2(counts(dds,normalized=TRUE)+1)

par(mfrow = c(1,3)) #split plot into 3 columns
boxplot(lgc.raw, main="raw counts")
boxplot(lgc.norm, main="norm counts")
boxplot(assay(rld),main="rlog")
