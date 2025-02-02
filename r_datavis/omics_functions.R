#install.packages("ggplot2")
#install.packages("ggrepel")

#install.packages("xtable")

#install.packages("reshape2")

#install.packages("amap")

#install.packages("BiocManager")
#BiocManager::install("clusterProfiler")
#install.packages("stringi")
#BiocManager::install("org.Mm.eg.db")

#BiocManager::install("STRINGdb")

library(ggplot2)
library(ggrepel)


library("xtable")
options(xtable.floating=FALSE)
options(xtable.timestamp="")


library(reshape2)

library(amap)

library(STRINGdb)

library(clusterProfiler)
library(org.Mm.eg.db)

my_theme = theme(
  plot.title = element_text(size=30),
  axis.text.x = element_text(size=8),
  axis.text.y = element_text(size=8),
  axis.title.x = element_text(size=12),
  axis.title.y = element_text(size=12)
)

load_tables = function(de_tables, em, annotations, symbol_column = 'SYMBOL'){
  master = merge(em, annotations,by.x=0,by.y=0)
  i = 1
  for (de in de_tables) {
    master = merge(master,de,by.x=1,by.y=0)
    colnames(master)[which(names(master) == "p")] = paste("p_",as.character(i),sep='')
    colnames(master)[which(names(master) == "p.adj")] = paste("p.adj_",as.character(i),sep='')
    colnames(master)[which(names(master) == "log2fold")] = paste("log2fold_",as.character(i),sep='')
    i = i + 1
  }
  master= na.omit(master)
  row.names(master) = master[,symbol_column]
  names(master)[1] = 'ensemble_id'
  return(master)
}

# Helper function to get data for specific gene
get_gene_data = function(gene, gene_frame, sample_groups, group_order=c() ) {
  gene_data = gene_frame[gene,]
  gene_data = t(gene_data)
  gene_data = data.frame(gene_data)
  gene_data$sample_group = sample_groups
  #names(gene_data) = c("expression","sample_group")
  
  #if length(group_order) > 0 {
  gene_data$sample_group = factor(gene_data$sample_group, levels=group_order)
  #}
  return (gene_data)
}

boxplot_facets = function(gene_frame, candidate_genes, sample_groups, nrow, ncol, group_order=c()){
  genes_to_plot = get_gene_data(candidate_genes,gene_frame,sample_groups,group_order=group_order)
  top_gene_data_m = melt(genes_to_plot, id.vars='sample_group')
  
  ggp = ggplot(top_gene_data_m,aes(x=sample_group,y=value, fill=sample_group)) +
    geom_boxplot() +
    facet_wrap(~variable, nrow=nrow, ncol=ncol) +
    my_theme
  return(ggp)
}

volcano_plot_df_table = function(df, p_max = 0.05, log2Fold_threshold = 1, name_column = "symbol", p_column = 'p.adj', log2Fold_column = 'log2Fold', symbol_labels = TRUE) {
  # adding label to de tables for up and down regulated genes
  df$diffexpr = "NO" 
  df$symbol = row.names(df)
  #df$log2Fold = df[log2Fold_column]
  df$diffexpr[df[log2Fold_column] > log2Fold_threshold & df[p_column] < p_max] = "UP"
  df$diffexpr[df[log2Fold_column] < -log2Fold_threshold  & df[p_column] < p_max] = "DOWN"
  
  df_top_ten = df[1:10,]
  
  # adding name labels for the significant genes
  df$delabel = NA
  if (symbol_labels) {
    df$delabel[df$diffexpr != "NO"] = df$symbol[df$diffexpr != "NO"]
  }
  
  ggp = ggplot(data=df, aes(x=log2fold, y=-log10(p.adj), col=diffexpr, label=symbol)) + 
    geom_point() +
    theme_minimal() +
    geom_text_repel(data = df_top_ten, max.overlaps=100, show.legend = FALSE) +
    scale_color_manual(values=c("darkcyan", "black", "darkred"), labels=c("Down-regulated","Non-significant", "Up-regulated"),name="") +
    geom_vline(xintercept=c(-log2Fold_threshold , log2Fold_threshold ), col="red") +
    geom_hline(yintercept=-log10(p_max), col="red") +
    my_theme
    #theme(legend.position="none")
  ggp
}

ma_plot_df_table = function(df, p_max = 0.05, log2Fold_threshold = 1, name_column = "symbol", p_column = 'p.adj', log2Fold_column = 'log2Fold', symbol_labels = TRUE) {
  # adding label to de tables for up and down regulated genes
  df$diffexpr = "NO" 
  #df$symbol = df[name_column]
  #df$log2Fold = df[log2Fold_column]
  df$diffexpr[df[log2Fold_column] > log2Fold_threshold & df[p_column] < p_max] = "UP"
  df$diffexpr[df[log2Fold_column] < -log2Fold_threshold  & df[p_column] < p_max] = "DOWN"
  
  # adding name labels for the significant genes
  df$delabel = NA
  if (symbol_labels) {
    df$delabel[df$diffexpr != "NO"] = df$symbol[df$diffexpr != "NO"]
  }
  
  ggp = ggplot(data=df, aes(x=log10(mean_expr), y=log2fold, col=diffexpr, label=delabel)) + 
    geom_point() +
    theme_minimal() +
    geom_text_repel(max.overlaps=100) +
    scale_color_manual(values=c("darkcyan", "black", "darkred")) +
    theme(legend.position="none")
  ggp
}

pca_graph = function(em_scaled,ss){
  pca = prcomp(t(as.matrix(sapply(em_scaled,as.numeric))))
  pca_coordinates = data.frame(pca$x)
  
  vars = apply(pca$x, 2, var)
  prop_x = round(vars["PC1"] / sum(vars),4) * 100
  prop_y = round(vars["PC2"] / sum(vars),4) * 100
  x_axis_label = paste("PC1 ", " (",prop_x, "%)",sep="")
  y_axis_label = paste("PC2 ", " (",prop_y, "%)",sep="")
  
  ggp = ggplot(pca_coordinates, aes(x=PC1, y=PC2,col=ss$SAMPLE_GROUP)) +
    geom_point() +
    geom_text_repel(aes(label=ss$SAMPLE),show.legend = FALSE) +
    scale_color_manual(values=as.vector(c("darkcyan", "black", "darkred"))) +
    xlab(x_axis_label) +
    ylab(y_axis_label)
  
  return(ggp)
}

expr_density_graph = function(table,sample){
  ggp = ggplot(table, aes(x = log10(sample+0.01))) + 
    geom_density(colour ='red', fill='coral', linewidth = 0.5, alpha = 0.5)
}

expr_density_facets = function(table, nrow, ncol){
  table_melt = melt(table)
  ggp = ggplot(table_melt, aes(x = log10(value+0.01))) + 
    geom_density(colour ='red', fill='coral', linewidth = 0.5, alpha = 0.5) +
    facet_wrap(~variable, nrow=nrow, ncol=ncol) +
    my_theme +
    labs(x = "Expression (log10)", y = "Density") +
    theme( panel.grid = element_blank())
  
  return(ggp)
}

plot_heatmap = function(em_table, dist_method = "spearman", cluster_method = "average", reorder_func = "average", by_x = FALSE){
  hm_matrix = as.matrix(em_table)
  if (by_x){
  hm_matrix = t(hm_matrix)
  }
  # Get the distances and cluster
  dist = Dist(hm_matrix,nbproc = 2, method="spearman")
  cluster = hclust(dist, method="average")
  
  # get cluster dendrogram and untangle
  tree = as.dendrogram(cluster)
  tree_reorder = reorder(tree,0,FUN=reorder_func)
  
  # Get untangled order and reorder the original matrix
  order = order.dendrogram(tree_reorder)
  hm_matrix_clustered = hm_matrix[order,]
  #hm_matrix_clustered = t(hm_matrix_clustered)
  
  # melt the matrix
  hm_matrix_clustered_melt = melt(hm_matrix_clustered)
  
  
  # Create heatmap colour palette
  pal_colours = c("purple","black","yellow")
  heatmap_palette = colorRampPalette(pal_colours)(100)
  
  ggp = ggplot(hm_matrix_clustered_melt, aes(x=Var2, y=Var1, fill=value)) + 
    geom_tile() + 
    scale_fill_gradientn(colours = heatmap_palette) +
    ylab("") +
    xlab("") +
    my_theme +
    theme(axis.text.y = element_blank(), axis.ticks=element_blank(), legend.title = element_blank(),
          legend.spacing.x = unit(0.25, 'cm'))
  return(ggp)
}

heatmap_rug = function(sample_groups){
  # rug for continuous variable
  groups_data = as.matrix(as.numeric(as.factor(sample_groups)))
  groups_data = melt(groups_data)
  # heatmap
  colours = c("red","cyan","purple")
  ggp = ggplot(groups_data, aes(x = Var1, y = Var2, fill = value)) + 
    geom_tile() + 
    scale_fill_gradientn(colours = colours) +
    geom_tile(linetype="blank") +
    labs(x = "", y = "") +
    theme(legend.position="none", legend.title = element_blank(), axis.text.x = element_blank(), axis.text.y =
            element_blank(), axis.ticks=element_blank())
  return(ggp)
}

# Returns the results of an ORA on given list of genes
get_ora_results = function(genes){
  # convert gene symbols to entrezid
  genes_entrez = bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  
  # run an ORA and store the results
  ora_results = enrichGO(gene = genes_entrez$ENTREZID, OrgDb = org.Mm.eg.db, readable = T, ont =
                           "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
  return(ora_results)
}

# Takes ora results and converts to data frame
convert_ora_results_to_table = function(ora_results) {
  gene_sets = ora_results$geneID
  description = ora_results $Description
  p_adj = ora_results$p.adjust
  
  ora_results_table = data.frame(cbind(description, gene_sets,  p_adj),row.names = 1)
  return(ora_results_table)
}

# Takes gsea results and convets to data frame
convert_gsea_results_to_table = function(gsea_results){
  description = gse_results$Description
  p_adj = gse_results$p.adjust
  NES = gse_results$NES
  gene_sets = gse_results$core_enrichment
  
  gsea_results_table = data.frame(cbind(description, gene_sets, p_adj, NES),row.names = 1)
  return(gsea_results_table)
}

# Returns list of genes from ora/gsea results table
get_enriched_genes_from_table = function(results_table, row_num){
  enriched_gene_set = as.character(results_table [row_num,1])
  candidate_genes = unlist(strsplit(enriched_gene_set, "/"))
  return(candidate_genes)
}

