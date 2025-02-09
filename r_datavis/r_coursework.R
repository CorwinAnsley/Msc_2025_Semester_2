source("./omics_functions.R")

em = read.table("./data/EM.csv", header=TRUE, row.names=1, sep= "\t")
de_senes_vs_prolif = read.table("./data/DE_Senes_vs_Prolif.csv", header=TRUE, row.names=1, sep= "\t")
de_mtd_vs_prolif = read.table("./data/DE_Senes_MtD_vs_Prolif.csv", header=TRUE, row.names=1, sep= "\t")
de_mtd_vs_senes = read.table("./data/DE_Senes_MtD_vs_Senes.csv", header=TRUE, row.names=1, sep= "\t")
annotations = read.table("./data/Human_Background_GRCh38.p13.csv", header=TRUE, row.names=1, sep= "\t")
ss = read.table("./data/sample_sheet.csv", header=TRUE, sep="\t")

#de_mtd_vs_prolif = na.omit(de_mtd_vs_prolif)
#de_mtd_vs_prolif = na.omit(de_mtd_vs_prolif)
#de_mtd_vs_prolif = na.omit(de_mtd_vs_prolif)

de_tables = list(de_senes_vs_prolif, de_mtd_vs_senes, de_mtd_vs_prolif)
#de_tables = list(de_senes_vs_prolif)

master = load_tables(de_tables, em, annotations)

em_symbols = master[ , as.vector(ss$SAMPLE)]
#em_symbols = em_symbols[,-1]
em_scaled = data.frame(t(scale(data.frame(t(em_symbols)))))
em_scaled = na.omit(em_scaled)


sig_genes_senes_vs_prolif = get_sig_genes(master, 'p.adj_1', 'log2fold_1')
sig_genes_mtd_vs_senes = get_sig_genes(master, 'p.adj_2', 'log2fold_2')
sig_genes_mtd_vs_prolif = get_sig_genes(master, 'p.adj_3', 'log2fold_3')

em_scaled_sig_senes_vs_prolif = em_scaled[sig_genes_senes_vs_prolif,]
em_scaled_sig_mtd_vs_senes = em_scaled[sig_genes_mtd_vs_senes,]
em_scaled_sig_mtd_vs_prolif = em_scaled[sig_genes_mtd_vs_prolif,]

# Create expression density plot
ggp = expr_density_facets(em, 3, 3)
ggsave("./plots/expr_density.pdf")

# Create pca plot
ggp = pca_graph(em_scaled, ss)
ggsave("./plots/pca.pdf")

source("./omics_functions.R")

# Create senes vs prolif volcano plot
ggp = volcano_plot_df_table(master, p_column = 'p.adj_1', log2Fold_column = 'log2fold_1')
ggp
ggsave("./plots/volcano_senes_vs_prolif.pdf", width = 20, height = 20)

# Create mtd vs senes volcano plot
ggp = volcano_plot_df_table(master, p_column = 'p.adj_2', log2Fold_column = 'log2fold_2')
ggsave("./plots/volcano_mtd_vs_senes.pdf")

# Create mtd vs prolif volcano plot
ggp = volcano_plot_df_table(master, p_column = 'p.adj_3', log2Fold_column = 'log2fold_3')
ggsave("./plots/volcano_mtd_vs_prolif.pdf")

# Create heatmaps
ggp = plot_heatmap(em_scaled_sig_senes_vs_prolif)
ggsave("./plots/heatmap_senes_vs_prolif.pdf", width = 9, height = 9)

ggp = plot_heatmap(em_scaled_sig_mtd_vs_senes)
ggsave("./plots/heatmap_mtd_vs_senes.pdf", width = 9, height = 9)

ggp = plot_heatmap(em_scaled_sig_mtd_vs_prolif)
ggsave("./plots/heatmap_mtd_vs_prolif.pdf", width = 9, height = 9)

# Create heatmap rug
ggp = heatmap_rug(ss$SAMPLE_GROUP)
ggsave("./plots/heatmap_rug.pdf", width = 9, height = 1)

# Get ORA results
ora_results = get_ora_results(sig_genes_senes_vs_prolif)
ora_results_table = convert_ora_results_to_table(ora_results)