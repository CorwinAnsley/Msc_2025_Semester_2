install.packages("ggplot2")

#library(ggplot2)

em = read.table("./data/em.csv", header=TRUE, row.names=1, sep= "\t")
de = read.table("./data/de_duct_vs_gut.csv", header=TRUE, row.names=1, sep= "\t")
annotations = read.table("./data/annotations.csv", header=TRUE, row.names=1, sep= "\t")
em_symbols = read.table("./data/em_symbols.csv", header=TRUE, row.names=1, sep= "\t")
em_scaled = read.table("./data/em_scaled.csv", header=TRUE, row.names=1, sep= "\t")
master= read.table("./data/master.csv", header=TRUE, row.names=1, sep= "\t")
master_sig = read.table("./data/master_sig.csv", header=TRUE, row.names=1, sep= "\t")
em_symbols_sig = read.table("./data/em_symbols_sig.csv", header=TRUE, row.names=1, sep= "\t")
em_scaled_sig = read.table("./data/em_scaled_sig.csv", header=TRUE, row.names=1, sep= "\t")