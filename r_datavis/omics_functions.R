library(ggplot2)
library(ggrepel)

#install.packages("xtable")
library("xtable")
options(xtable.floating=FALSE)
options(xtable.timestamp="")
library(crayon)

my_theme = theme(
  plot.title = element_text(size=30),
  axis.text.x = element_text(size=8),
  axis.text.y = element_text(size=8),
  axis.title.x = element_text(size=12),
  axis.title.y = element_text(size=12)
)

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