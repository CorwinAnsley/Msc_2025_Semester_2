#### Install and load dependencies ####
#install.packages("ggplot2")
#install.packages("ggrepel")

#install.packages("xtable")

#install.packages("reshape2")

#install.packages("amap")

#install.packages("BiocManager")
#install.packages("stringi")
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Mm.eg.db")

#BiocManager::install("STRINGdb")

#install.packages('stringr')

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

library('stringr')

#### Set Theme ####

# Custom theme for all plots
my_theme = theme(
  plot.title = element_text(size=30),
  axis.text.x = element_text(size=12),
  axis.text.y = element_text(size=12),
  axis.title.x = element_text(size=16),
  axis.title.y = element_text(size=16),
  legend.title = element_blank()
) + theme_minimal()

#### Data functions ####

# Returns a master table colating the data from provided em, annotations and list of de tables
load_master_table = function(de_tables, de_names, em, annotations, symbol_column = 'SYMBOL'){
  master = merge(em, annotations,by.x=0,by.y=0)
  i = 1
  # Add p, p.adj and log2fold for all de tables
  for (de in de_tables) {
    master = merge(master,de,by.x=1,by.y=0)
    colnames(master)[which(names(master) == "p")] = paste("p_",de_names[i],sep='')
    colnames(master)[which(names(master) == "p.adj")] = paste("p.adj_",de_names[i],sep='')
    colnames(master)[which(names(master) == "log2fold")] = paste("log2fold_",de_names[i],sep='')
    i = i + 1
  }
  
  # Clean up the master table and return it
  master = na.omit(master)
  row.names(master) = master[,symbol_column]
  names(master)[1] = 'ensemble_id'
  master = master[,-which(names(master) == symbol_column)]
  return(master)
}

# Returns a list of sig genes for given p (should be adjusted) and log2fold column
get_sig_genes = function(master, p_column, log2fold_column, p_threshold = 0.001, log2fold_threshold=2) {
  master['p.adj'] = master[p_column]
  master['log2fold'] = master[log2fold_column]
  master_sig = subset(master, p.adj < p_threshold)
  master_sig = subset(master, abs(log2fold) >log2fold_threshold)
  sig_genes = row.names(master_sig)
  return(sig_genes)
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

#### Boxplot functions ####
# Create facets of box plots
boxplot_facets = function(gene_frame, candidate_genes, sample_groups, nrow, ncol, group_order=c()){
  genes_to_plot = get_gene_data(candidate_genes,gene_frame,sample_groups,group_order=group_order)
  top_gene_data_m = melt(genes_to_plot, id.vars='sample_group')
  
  ggp = ggplot(top_gene_data_m,aes(x=sample_group,y=value, fill=sample_group)) +
    geom_boxplot() +
    facet_wrap(~variable, nrow=nrow, ncol=ncol) +
    my_theme
  return(ggp)
}

# Create a metagene by averaging out columns and return a boxplot of its expression
metagene_boxplot = function(gene_frame, sample_groups)
{
  gene_frame = na.omit(gene_frame)
  metagene = data.frame(colMeans(x=gene_frame))
  metagene['sample_group'] = sample_groups
  colnames(metagene) = c("expression", "sample_group")
  #metagene = t(metagene)
  #return(metagene)
  
  ggp = ggplot(metagene,aes(x=sample_group,y=expression, fill=sample_group)) +
    geom_boxplot() +
    scale_color_manual(values=as.vector(c("darkcyan", "black", "darkred")))+
    xlab("") +
    my_theme +
    theme(legend.title = element_blank())
  return(ggp)
}

#### General Plots ####

# Create a volcano plot
volcano_plot_df_table = function(df,
                                 p_column = 'p.adj', 
                                 log2Fold_column = 'log2Fold', 
                                 plot_title = '',
                                 p_threshold = 0.001, 
                                 log2fold_threshold = 2) {
  # adding label to de tables for up and down regulated genes
  df$diffexpr = "NO" 
  df$symbol = row.names(df)
  df['log2fold'] = df[log2Fold_column]
  df['p.adj'] = df[p_column]
  df$diffexpr[df[log2Fold_column] > log2fold_threshold & df[p_column] < p_threshold] = "UP"
  df$diffexpr[df[log2Fold_column] < -log2fold_threshold  & df[p_column] < p_threshold] = "DOWN"
  
  sorted_order = order(df[,'p.adj'], decreasing=FALSE)
  df = df[sorted_order,]
  
  df_sig_up = subset(df, p.adj < p_threshold & log2fold > log2fold_threshold)
  df_sig_down = subset(df, p.adj < p_threshold & log2fold < -log2fold_threshold)#
  
  df_up_top5 = df_sig_up[1:5,]
  df_down_top5 = df_sig_down[1:5,]

  ggp = ggplot(data=df, aes(x=log2fold, y=-log10(p.adj), col=diffexpr, label=symbol)) + 
    geom_point() +
    labs(title=plot_title) +
    #theme_minimal() +
    geom_text_repel(data = df_up_top5, max.overlaps=100, show.legend = FALSE) +
    geom_text_repel(data = df_down_top5, max.overlaps=100, show.legend = FALSE) +
    scale_color_manual(values=c("darkcyan", "black", "darkred"), labels=c("Down-regulated","Non-significant", "Up-regulated"),name="") +
    geom_vline(xintercept=c(-log2fold_threshold , log2fold_threshold ), col="red") +
    geom_hline(yintercept=-log10(p_threshold), col="red") +
    my_theme 
  

  return(ggp)
  #return(df)
}

# Create ma plot
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

# Create pca plot
pca_graph = function(em_scaled,ss){
  pca = prcomp(t(as.matrix(sapply(em_scaled,as.numeric))))
  pca_coordinates = data.frame(pca$x)
  
  vars = apply(pca$x, 2, var)
  prop_x = round(vars["PC1"] / sum(vars),4) * 100
  prop_y = round(vars["PC2"] / sum(vars),4) * 100
  x_axis_label = paste("PC1 ", " (",prop_x, "%)",sep="")
  y_axis_label = paste("PC2 ", " (",prop_y, "%)",sep="")
  
  ggp = ggplot(pca_coordinates, aes(x=PC1, y=PC2, colour=ss$SAMPLE_GROUP)) +
    geom_point() +
    geom_text_repel(aes(label=ss$SAMPLE),show.legend = FALSE) +
    scale_color_manual(values=as.vector(c("darkcyan", "black", "darkred"))) +
    labs(color = "Sample Group\n") +
    xlab(x_axis_label) +
    ylab(y_axis_label) +
    my_theme
  
  return(ggp)
}

# Create expression density graph
expr_density_graph = function(table,sample){
  ggp = ggplot(table, aes(x = log10(sample+0.01))) + 
    geom_density(colour ='red', fill='coral', linewidth = 0.5, alpha = 0.5)
}

# Create facets of expression density graphs
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

#### Heatmap functions ####

# plot a heatmap
plot_heatmap = function(em_table, dist_method = "spearman", cluster_method = "average", reorder_func = "average", by_x = FALSE){
  em_taBLE = na.omit(em_table)
  
  hm_matrix = as.matrix(em_table)
  if (by_x){
  hm_matrix = t(hm_matrix)
  }
  hm_matrix = na.omit(hm_matrix)
  # Get the distances and cluster
  dist = Dist(hm_matrix,nbproc = 2, method="spearman")
  #return(dist)
  dist = na.omit(dist)
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
    theme(axis.text.y = element_blank(), axis.ticks=element_blank(), legend.title = element_blank(),# axis.text.x = element_text(size=8),
          legend.spacing.x = unit(0.25, 'cm'))
  return(ggp)
}

#Create a rug for a heatmap
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
    my_theme +
    theme(legend.position="none", legend.title = element_blank(), axis.text.x = element_blank(), axis.text.y =
            element_blank(), axis.ticks=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())
  return(ggp)
}

#### Pathway and signature analysis functions ####

# Returns the results of an ORA on given list of genes
get_ora_results = function(genes){
  genes = stringr::str_to_title(genes) 
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

# Creates a heatmap, a metagene boxplot and ora pathway analysis from a list of genes or signature
plot_signature = function(gene_list, em_scaled, ss, signature_name)
{
  # get a subset of expression data from the list of genes
  em_scaled_signature = em_scaled[gene_list,]
  
  # Plot a heatmap
  ggp = plot_heatmap(em_scaled_signature)
  heatmap_filename = paste("./plots/heatmap_",signature_name,sep='')
  heatmap_filename = paste(heatmap_filename,".pdf",sep='')
  ggsave(heatmap_filename, width = 9, height = 9)
  
  #ggp = make_metagene_boxplot(gene_list, em_scaled, groups)
  #save_plot(ggp, …)
  
  # Plot a metagene boxplot
  ggp = metagene_boxplot(em_scaled_signature, ss$SAMPLE_GROUP)
  boxplot_filename = paste("./plots/boxplot_",signature_name,sep='')
  boxplot_filename = paste(boxplot_filename,".pdf",sep='')
  ggsave(boxplot_filename, width = 5, height = 5)
  
  # Run pathway analysis and create barplot
  ora_results = get_ora_results(gene_list)
  ora_results_table = convert_ora_results_to_table(ora_results)
  ggp = barplot(ora_results, showCategory=10) + labs(x = "") + my_theme + theme(legend.title = element_blank())
  barplot_filename = paste("./plots/barplot_",signature_name,sep='')
  barplot_filename = paste(barplot_filename,".pdf",sep='')
  ggsave(barplot_filename, width = 5, height = 5)
  
}
                        
