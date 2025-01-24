#install.packages("BiocManager")
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Mm.eg.db")

library(clusterProfiler)
library(org.Mm.eg.db)

em_scaled = read.table("./data/em_scaled.csv", header=TRUE, row.names=1, sep= "\t")
master= read.table("./data/master.csv", header=TRUE, row.names=1, sep= "\t")
master_sig = read.table("./data/master_sig.csv", header=TRUE, row.names=1, sep= "\t")
em_symbols_sig = read.table("./data/em_symbols_sig.csv", header=TRUE, row.names=1, sep= "\t")
em_scaled = na.omit(em_scaled)
sig_genes = row.names(master_sig)
em_scaled_sig = em_scaled[sig_genes,]
em_scaled_sig = na.omit(em_scaled_sig)
ss = read.table("./data/sample_sheet.csv", header=TRUE, sep="\t")

sig_genes_entrez = bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Task 3
ora_results = enrichGO(gene = sig_genes_entrez$ENTREZID, OrgDb = org.Mm.eg.db, readable = T, ont =
                         "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
ggp = barplot(ora_results, showCategory=10)
ggp

# Task 5
up_sig_genes = row.names(subset(master_sig, log2fold > 0))
down_sig_genes = row.names(subset(master_sig, log2fold < 0))



up_sig_genes_entrez = bitr(up_sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
ora_results = enrichGO(gene = up_sig_genes_entrez$ENTREZID, OrgDb = org.Mm.eg.db, readable = T, ont =
                         "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
ggp = barplot(ora_results, showCategory=10)
ggp

# Task 6
# we want the log2 fold change
gsea_input = master$log2fold
# add gene names the vector
names(gsea_input) = row.names(master)
# omit any NA values
gsea_input = na.omit(gsea_input)
# sort the list in decreasing order (required for clusterProfiler)
gsea_input = sort(gsea_input, decreasing = TRUE)

gse_results = gseGO(geneList=gsea_input,
                    ont ="BP",
                    keyType = "SYMBOL",
                    nPerm = 10000,
                    minGSSize = 3,
                    maxGSSize = 800,
                    pvalueCutoff = 0.05,
                    verbose = TRUE,
                    OrgDb = org.Mm.eg.db,
                    pAdjustMethod = "none")
ggp = ridgeplot(gse_results)

### Tutorial 11
ora_results = enrichGO(gene = sig_genes_entrez$ENTREZID, OrgDb = org.Mm.eg.db, readable = T, ont =
                         "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.10)

source("./omics_functions.R")

e_genes = get_enriched_genes_ora(ora_results)
#\print(length(e_genes))

ggp = boxplot_facets(em_scaled, e_genes, ss$SAMPLE_GROUP, 5, 26, group_order=c('gut','duct','node'))
ggp

e_genes_scaled = em_scaled[e_genes,]
e_genes_scaled = na.omit(e_genes_scaled)

ggp = plot_heatmap(e_genes_scaled)
ggp


