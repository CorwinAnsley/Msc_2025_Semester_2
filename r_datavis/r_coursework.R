#### Import functions ####
source("./omics_functions.R")


#### Load data files ####
em = read.table("./data/EM.csv", header=TRUE, row.names=1, sep= "\t")
de_senes_vs_prolif = read.table("./data/DE_Senes_vs_Prolif.csv", header=TRUE, row.names=1, sep= "\t")
de_mtd_vs_prolif = read.table("./data/DE_Senes_MtD_vs_Prolif.csv", header=TRUE, row.names=1, sep= "\t")
de_mtd_vs_senes = read.table("./data/DE_Senes_MtD_vs_Senes.csv", header=TRUE, row.names=1, sep= "\t")
annotations = read.table("./data/Human_Background_GRCh38.p13.csv", header=TRUE, row.names=1, sep= "\t")
ss = read.table("./data/sample_sheet.csv", header=TRUE, sep="\t")

#de_mtd_vs_prolif = na.omit(de_mtd_vs_prolif)
#de_mtd_vs_prolif = na.omit(de_mtd_vs_prolif)
#de_mtd_vs_prolif = na.omit(de_mtd_vs_prolif)

#### Get tables and sig genes ####

de_tables = list(de_senes_vs_prolif, de_mtd_vs_senes, de_mtd_vs_prolif)
de_names = c('senes_vs_prolif', 'mtd_vs_senes', 'mtd_vs_prolif')
#de_tables = list(de_senes_vs_prolif)

master = load_master_table(de_tables, de_names, em, annotations)

#names(df)[names(df) == 'old.var.name'] <- 'new.var.name'

em_symbols = master[ , as.vector(ss$SAMPLE)]
#em_symbols = em_symbols[,-1]
em_scaled = data.frame(t(scale(data.frame(t(em_symbols)))))
em_scaled = na.omit(em_scaled)


sig_genes_senes_vs_prolif = get_sig_genes(master, 'p.adj_1', 'log2fold_1')
sig_genes_mtd_vs_senes = get_sig_genes(master, 'p.adj_2', 'log2fold_2')
sig_genes_mtd_vs_prolif = get_sig_genes(master, 'p.adj_3', 'log2fold_3')
all_sig_genes = c(sig_genes_senes_vs_prolif, sig_genes_mtd_vs_senes, sig_genes_mtd_vs_prolif)
all_sig_genes = unique(all_sig_genes)


em_scaled_sig_senes_vs_prolif = em_scaled[sig_genes_senes_vs_prolif,]
em_scaled_sig_mtd_vs_senes = em_scaled[sig_genes_mtd_vs_senes,]
em_scaled_sig_mtd_vs_prolif = em_scaled[sig_genes_mtd_vs_prolif,]
em_scaled_all_sig = em_scaled[all_sig_genes,]

#### Figure 1 ####

# Create expression density plot
ggp = expr_density_facets(em, 3, 3)
ggsave("./plots/expr_density.pdf", width = 9, height = 9)

#### Figure 2 ####

# Create pca plot
ggp = pca_graph(em_scaled, ss)
ggsave("./plots/pca.pdf", width = 9, height = 9)


#### Figure 3 ####

# Create senes vs prolif volcano plot
ggp = volcano_plot_df_table(master, p_column = 'p.adj_senes_vs_prolif', log2Fold_column = 'log2fold_senes_vs_prolif',plot_title = 'Senes vs Prolif')
ggsave("./plots/volcano_senes_vs_prolif.pdf", width = 9, height = 9)

# Create mtd vs senes volcano plot
ggp = volcano_plot_df_table(master, p_column = 'p.adj_mtd_vs_senes', log2Fold_column = 'log2fold_mtd_vs_senes',plot_title = 'SenesMtD vs Senes')
ggsave("./plots/volcano_mtd_vs_senes.pdf", width = 9, height = 9))

# Create mtd vs prolif volcano plot
ggp = volcano_plot_df_table(master, p_column = 'p.adj_mtd_vs_prolif', log2Fold_column = 'log2fold_mtd_vs_prolif',plot_title = 'SenesMtD vs Prolif')
ggsave("./plots/volcano_mtd_vs_prolif.pdf", width = 9, height = 9))

#### Figure 4 ####

# Create heatmaps
ggp = plot_heatmap(em_scaled_sig_senes_vs_prolif)
ggsave("./plots/heatmap_senes_vs_prolif.pdf", width = 9, height = 9)

ggp = plot_heatmap(em_scaled_sig_mtd_vs_senes)
ggsave("./plots/heatmap_mtd_vs_senes.pdf", width = 9, height = 9)

ggp = plot_heatmap(em_scaled_sig_mtd_vs_prolif)
ggsave("./plots/heatmap_mtd_vs_prolif.pdf", width = 9, height = 9)

ggp = plot_heatmap(em_scaled_all_sig)
ggp
ggsave("./plots/heatmap_all.pdf", width = 9, height = 9)

# Create heatmap rug
ggp = heatmap_rug(ss$SAMPLE_GROUP)
ggsave("./plots/heatmap_rug.pdf", width = 9, height = 1)

# Get ORA results
ora_results = get_ora_results(sig_genes_senes_vs_prolif)
ora_results_table = convert_ora_results_to_table(ora_results)

senes_vs_prolif_e_genes = get_enriched_genes_from_table(ora_results_table, 1)
#senes_vs_prolif_e_genes = get_enriched_genes_ora(ora_results)

senes_vs_prolif_e_genes_scaled = em_scaled[senes_vs_prolif_e_genes,]
senes_vs_prolif_e_genes_scaled = na.omit(senes_vs_prolif_e_genes_scaled)

ggp = boxplot_facets(senes_vs_prolif_e_genes_scaled, senes_vs_prolif_e_genes, ss$SAMPLE_GROUP, 5, 26, group_order=c('Prolif','Senes','Senes_MtD'))
ggp

#### Figure 5 ####

signature_1 = unique(c(sig_genes_senes_vs_prolif, sig_genes_mtd_vs_prolif))
signature_1 = master[signature_1,]
signature_1 = subset(signature_1, log2fold_senes_vs_prolif > 0 & log2fold_mtd_vs_prolif > 0)
signature_1 = row.names(signature_1)
em_scaled_signature_1 = em_scaled[signature_1,]

ggp = plot_heatmap(em_scaled_signature_1)
ggp

source("./omics_functions.R")

plot_signature(signature_1, em_scaled, ss, "signature_1")

ora_results_s1 = get_ora_results(signature_1)
ora_results_table_s1 = convert_ora_results_to_table(ora_results_s1)

ggp = barplot(ora_results_s1, showCategory=10) + my_theme + theme(legend.title = element_blank())
ggsave("./plots/barplot_signature_1.pdf", width = 9, height = 9)
ggp

  #row.names(subset(master, (sig_genes_senes_vs_prolif == TRUE &log2fold_1 > 0) & (sig_genes_mtd_vs_senes == TRUE & log2fold_2 > 0)))
