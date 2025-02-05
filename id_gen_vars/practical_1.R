anova_data = read.table("./data/ANOVAdata.txt", header = TRUE, stringsAsFactors = T)

class(anova_data$housing)

summary(anova_data)

levelsbree