suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
library(GOplot)
library("enrichplot")
#载入基因
res <- read.table("old5_vs_young5_DESeq2.txt", sep = "\t")

all <- rownames(res)
up <- rownames(res)[which(res$type=="up")]
down <- rownames(res)[which(res$type=="down")]

#上调基因
tp <- c("BP", "MF", "CC")
for(i in 1:length(tp)){
  ego.up <- enrichGO(gene = up,
                     universe = all,
                     OrgDb = org.Hs.eg.db,
                     keyType = 'SYMBOL',
                     ont = tp[i],
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05,
                     readable = TRUE)
  png(paste0(tp[i], "_GO_up.png"), width = 2500, height = 2500, res = 284)
  p <- goplot(ego.up)
  print(p)
  dev.off()
  
  png(paste0(tp[i], "_GO_up_dot.png"), width = 2500, height = 2500, res = 284)
  p <- dotplot(ego.up)
  print(p)
  dev.off()
}
#下调
for (i in 1:length(tp)) {
  ego.down <- enrichGO(gene = down,
                     universe = all,
                     OrgDb = org.Hs.eg.db,
                     keyType = 'SYMBOL',
                     ont = tp[i],
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05,
                     readable = TRUE)
  png(paste0(tp[i], "_GO_down.png"), width = 2500, height = 2500, res = 284)
  p <- goplot(ego.down)
  print(p)
  dev.off()
  
  png(paste0(tp[i], "_GO_down_dot.png"), width = 2500, height = 2500, res = 284)
  p <- dotplot(ego.down)
  print(p)
  dev.off()
}


#KEGG

all <- bitr(rownames(res), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
all <- all[,2]
up <- rownames(res)[which(res$type=="up")]
up <- bitr(up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
up <- up[,2]
down <- rownames(res)[which(res$type=="down")]
down <- bitr(down, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
down <- down[,2]



# 上调
kk.up <- enrichKEGG(gene = up,
                    organism = 'hsa',
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)
png(paste0("KEGG_up_barplot.png"), width = 2500, height = 2500, res = 284)
barplot(kk.up, showCategory = 10)
dev.off()

png(paste0("KEGG_up_dotplot.png"), width = 2500, height = 2500, res = 284)
dotplot(kk.up, showCategory = 20)
dev.off()

png(paste0("KEGG_up_plot.png"), width = 2500, height = 2500, res = 284)
cnetplot(kk.up) #直接使用enrichKEGG函数分析出的结果，ek2中ID已转换为symbol
dev.off()

#下调
kk.down <- enrichKEGG(gene = down,
                      organism = 'hsa',
                      pvalueCutoff = 0.05)

png(paste0("KEGG_down_barplot.png"), width = 2500, height = 2500, res = 284)
barplot(kk.down, showCategory = 5)
dev.off()

png(paste0("KEGG_down_dotplot.png"), width = 2500, height = 2500, res = 284)
dotplot(kk.down, showCategory = 10)
dev.off()

png(paste0("KEGG_down_plot.png"), width = 2500, height = 2500, res = 284)
cnetplot(kk.down) #直接使用enrichKEGG函数分析出的结果，ek2中ID已转换为symbol
dev.off()






