anova_data = read.table("./data/ANOVAdata.txt", header = TRUE, stringsAsFactors = T)

class(anova_data$housing)

summary(anova_data)

levels(anova_data$breed)

quantile(anova_data$weight, prob=c(0.025,0.975))

t.test(anova_data$weight ~ anova_data$housing)

wilcox.test(anova_data$weight ~ anova_data$housing)

chisq.test(table(anova_data$housing,anova_data$breed))

## 1a
hist(anova_data$weight)

shapiro.test(anova_data$weight)

## 1b
plot(anova_data$tick, anova_data$weight, pch=15, col='red', cex.lab=1.5)

lm_1 = lm(weight~tick, data=anova_data)
summary(lm_1)
par(mfrow=c(2,2))
abline(lm(weight~tick, data=anova_data),col="blue",lwd=3) 

plot(lm_1)

par(mfrow=c(1,1))
plot(anova_data$tick, anova_data$weight, pch=15, col='red', cex.lab=1.5,
     abline(lm(weight~tick, data=anova_data),col="blue",lwd=3) )

## 2

boxplot(anova_data$weight~anova_data$breed, xlab="Breed", ylab="Weight",col="red")

aov1 = aov(weight~breed, data = anova_data)
summary(aov1)

pairwise.t.test(anova_data$weight, anova_data$breed, p.adj = "bonf")

pairwise.t.test(anova_data$weight, anova_data$breed, p.adj = "fdr")
