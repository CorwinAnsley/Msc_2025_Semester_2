source("./omics_functions.R")

em = read.table("./data/EM.csv", header=TRUE, row.names=1, sep= "\t")
de_mtd_vs_prolif = read.table("./data/DE_Senes_MtD_vs_Prolif.csv", header=TRUE, row.names=1, sep= "\t")
de_mtd_vs_senes = read.table("./data/DE_Senes_MtD_vs_Senes.csv", header=TRUE, row.names=1, sep= "\t")
de_senes_vs_prolif = read.table("./data/DE_Senes_vs_Prolif.csv", header=TRUE, row.names=1, sep= "\t")
annotations = read.table("./data/Human_Background_GRCh38.p13.csv", header=TRUE, row.names=1, sep= "\t")
ss = read.table("./data/sample_sheet.csv", header=TRUE, row.names=1, sep="\t")

de_tables = list(de_senes_vs_prolf, de_mtd_vs_senes, de_mtd_vs_prolif)


master = load_tables(de_tables, em, annotations)

em_symbols = master[,-c(11:17)]
em_symbols = em_symbols[,-1]

ggp = expr_density_facets(em_symbols, 3, 3)
ggp
#names(master)['p'] = 'p_1'