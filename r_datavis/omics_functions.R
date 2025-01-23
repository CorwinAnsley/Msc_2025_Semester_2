library(ggplot2)
library(ggrepel)

#install.packages("xtable")
library("xtable")
options(xtable.floating=FALSE)
options(xtable.timestamp="")
#library(crayon)

#install.packages("reshape2")
library(reshape2)

my_theme = theme(
  plot.title = element_text(size=30),
  axis.text.x = element_text(size=8),
  axis.text.y = element_text(size=8),
  axis.title.x = element_text(size=12),
  axis.title.y = element_text(size=12)
)

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

load_tables = function(de_tables, em, annotations){
  master = merge(em, annotations,by.x=0,by.y=0)
  for (de in de_tables) {
    master = merge(master_temp, de,by.x=1,by.y=0)
  }
  return(master)
  
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