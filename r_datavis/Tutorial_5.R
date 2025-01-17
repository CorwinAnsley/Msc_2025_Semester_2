#install.packages("ggplot2")

library(ggplot2)
library(ggrepel)

em = read.table("./data/em.csv", header=TRUE, row.names=1, sep= "\t")
de = read.table("./data/de_duct_vs_gut.csv", header=TRUE, row.names=1, sep= "\t")
annotations = read.table("./data/annotations.csv", header=TRUE, row.names=1, sep= "\t")
em_symbols = read.table("./data/em_symbols.csv", header=TRUE, row.names=1, sep= "\t")
em_scaled = read.table("./data/em_scaled.csv", header=TRUE, row.names=1, sep= "\t")
master= read.table("./data/master.csv", header=TRUE, row.names=1, sep= "\t")
master_sig = read.table("./data/master_sig.csv", header=TRUE, row.names=1, sep= "\t")
em_symbols_sig = read.table("./data/em_symbols_sig.csv", header=TRUE, row.names=1, sep= "\t")
#em_scaled_sig = read.table("./data/em_scaled_sig.csv", header=TRUE, row.names=1, sep= "\t")

# Task 2

  master_sig_up = subset(master_sig, log2fold>0)
  master_sig_down = subset(master_sig, log2fold<0)
  
  master_sig_up_top5 = master_sig_up[1:5,]
  master_sig_down_top5 = master_sig_down[1:5,]
  
  ggp = ggplot(master, aes(x=log2fold,y=-log10(p.adj))) + 
    geom_point(size=0.5) + 
    geom_point(data=master_sig_up, color='red') +
    geom_point(data=master_sig_down, color='blue') +
    labs(title="Volcano", x="log2fold", y="-log10(p)") + 
    geom_vline(xintercept= -1, linetype="dashed", color = "grey", linewidth=0.5) +
    geom_vline(xintercept= 1, linetype="dashed", color = "grey", linewidth=0.5) + 
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "grey", linewidth=0.5) +
    theme_classic() + 
    #xlim(c(-20, 20)) +
    #ylim(c(0, 50)) +
    geom_text_repel(data=master_sig_up_top5, aes(label=row.names(master_sig_up_top5))) +
    geom_text_repel(data=master_sig_down_top5, aes(label=row.names(master_sig_down_top5)))
  
  ggp
  
  #master$symbol = row.names(master)
  
  source("./omics_functions.R")
  
  volcano_plot_df_table(master,log2Fold_column = 'log2fold',symbol_labels = FALSE, name_column = 0)
  ma_plot_df_table(master,log2Fold_column = 'log2fold',symbol_labels = FALSE, name_column = 0)
  #png("./data/volcano.png", height = 400, width = 400)
  #print(ggp)
  #dev.off()
  
