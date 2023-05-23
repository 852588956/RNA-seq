library(GSVA)
library(dplyr)
library(openxlsx)
library(ggpubr)
#加载数据

types <- read.table("../data/types.txt", header = T, sep = "\t")


tpm <- read.table("../data/tpm.txt", header = T, row.names = 1)
tpm <- log(tpm+0.01)
tpm <- as.matrix(tpm)

#sen的
geneset <- read.xlsx("../data/41467_2022_32552_MOESM4_ESM.xlsx")
geneset <- list(as.character(as.matrix(geneset[, 1])))

#geneset <- read.table("../data/cellAge/cellage3.tsv", header = T, sep = "\t")
#geneset <- list(as.character(as.matrix(geneset[, 2])))






result=gsva(tpm, geneset, method='gsva', kcdf="Gaussian", verbose=T)
result=data.matrix(result)
rownames(result)="SenMayo"

#画图
senmayo <- t(result)
types$senmayo <- senmayo[rownames(types), 1]
G1 <- types[which(types$treatment == "Placebo"), ]

png("GSE141595_Placebo.png", width = 2500, height = 2500, res = 284)
ggboxplot(types, x="age", y="senmayo",
          color = "treatment", add = "jitter")+
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(types$senmayo), li1netype=2)+# Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 1)+ # Add global annova p-value
  stat_compare_means(label = "p.format", method = "t.test", ref.group = "25-37 years")
dev.off()




#选出最年轻的5个和最老的5个
G2  <- G1[order(G1$senmayo), ]
G3 <- rbind(G2[which(G2$age == "25-37 years"), ][1:5, ], G2[which(G2$age == "60-79 years"), ][11:15, ])
write.table(G3,"Placebo.txt", sep = "\t", quote = F)








