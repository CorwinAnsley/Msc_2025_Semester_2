em = read.table("./data/em.csv", header=TRUE, row.names=1, sep= "\t")
de = read.table("./data/de_duct_vs_gut.csv", header=TRUE, row.names=1, sep= "\t")
annotations = read.table("./data/annotations.csv", header=TRUE, row.names=1, sep= "\t")

# Tutorial 2

# Task 6
de['ENSMUSG00000024084',]
de[,'log2fold']
de['ENSMUSG00000045010','p']

# Task 7
new_names = c("chromosome","start","stop","name","type")
names(annotations) = new_names

# Task 8
master_temp = merge(em, annotations,by.x=0,by.y=0)
master = merge(master_temp, de,by.x=1,by.y=0)
row.names(master) = master[,'name']
names(master)[1] = 'ensemble_id'
master = master[,-14]

em_symbols = master[,-c(11:17)]
em_symbols = em_symbols[,-1]

# Tutorial 3+4

master = na.omit(master)

# Task 2
sorted_order = order(master[,'p'], decreasing=FALSE)
master = master[sorted_order,]

# Task 3
#rowMeans(em)
#rowMeans(master)
rowMeans(master[,c(2:10)])

master$mean_expr = rowMeans(master[,c(2:10)])

# Task 4
master$mlog10p = -log10(master$p)

# Task 5
master$sig = as.factor(master$p.adj < 0.05 & abs(master$log2fold) > 1.0)
em_scaled = data.frame(t(scale(data.frame(t(em_symbols)))))
em_scaled = na.omit(em_scaled)

# Task 7  

master_sig = subset(master, p.adj < 0.05)
master_sig = subset(master, abs(log2fold) >1))




