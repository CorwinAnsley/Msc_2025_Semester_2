library(ggplot2)
library(ggrepel)

em = read.table("./data/em.csv", header=TRUE, row.names=1, sep= "\t")
de_duct_vs_gut = read.table("./data/de_duct_vs_gut.csv", header=TRUE, row.names=1, sep= "\t")
de_node_vs_duct
annotations = read.table("./data/annotations.csv", header=TRUE, row.names=1, sep= "\t")

