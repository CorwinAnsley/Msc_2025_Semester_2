library(ggplot2)
library(ggrepel)

ss = read.table("./data/sample_sheet.csv", header=TRUE, sep="\t")
master= read.table("./data/master.csv", header=TRUE, row.names=1, sep= "\t")

em_symbols = master[ , as.vector(ss$SAMPLE)]
em_scaled = data.frame(t(scale(data.frame(t(em_symbols)))))
em_scaled = na.omit(em_scaled)

# Task 2
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

#ggp

# Task 3
ggp = ggplot(master, aes(x = log10(gut_r1+0.01))) + 
geom_density(colour ='coral', linewidth = 0.5, alpha = 5)

ggp
