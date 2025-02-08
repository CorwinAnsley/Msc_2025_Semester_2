source("./omics_functions.R")

em = read.table("./data/EM.csv", header=TRUE, row.names=1, sep= "\t")
de_mtd_vs_prolif = read.table("./data/DE_Senes_MtD_vs_Prolif.csv", header=TRUE, row.names=1, sep= "\t")
de_mtd_vs_senes = read.table("./data/DE_Senes_MtD_vs_Senes.csv", header=TRUE, row.names=1, sep= "\t")
de_senes_vs_prolif = read.table("./data/DE_Senes_vs_Prolif.csv", header=TRUE, row.names=1, sep= "\t")
annotations = read.table("./data/Human_Background_GRCh38.p13.csv", header=TRUE, row.names=1, sep= "\t")
ss = read.table("./data/sample_sheet.csv", header=TRUE, sep="\t")

#de_mtd_vs_prolif = na.omit(de_mtd_vs_prolif)
#de_mtd_vs_prolif = na.omit(de_mtd_vs_prolif)
#de_mtd_vs_prolif = na.omit(de_mtd_vs_prolif)

de_tables = list(de_senes_vs_prolif, de_mtd_vs_senes, de_mtd_vs_prolif)

master = load_tables(de_tables, em, annotations)

em_symbols = master[ , as.vector(ss$SAMPLE)]
#em_symbols = em_symbols[,-1]
em_scaled = data.frame(t(scale(data.frame(t(em_symbols)))))
em_scaled = na.omit(em_scaled)

ggp = expr_density_facets(em, 3, 3)
#ggsave("./plots/expr_density.pdf")

source("./omics_functions.R")
ggp = pca_graph(em_scaled, ss)
#ggsave("./plots/pca.pdf")
ggp

test_df = volcano_plot_df_table(master, p_column = 'p.adj_1', log2Fold_column = 'log2fold_1')
ggp