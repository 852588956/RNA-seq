library(DESeq2)

library(stringr)
#加载数据"
group <- read.table("sample_group.txt", header = T, row.names = 1, sep = "\t")
old <- rownames(group)[which(group$group == "old")]
young <- rownames(group)[which(group$group == "young")]

counts <- read.table("gene_count_matrix.csv", header = T, row.names = 1, sep = ",")

genename <- rownames(counts)
for(i in 1:length(genename)){
  genename[i] <- str_split(genename[i], "[|]")[[1]][2]
}
counts$genename <- genename
genes <- unique(genename)
#重复的基因取平均值
gene_counts <- counts[1, 1:10]
for(i in 1:length(genes)){
  gn <- genes[i]
  gs <- counts[which(counts$genename == gn), 1:10]
  gene_counts[gn, ] <- as.integer(colMeans(gs))
}
gene_counts <- gene_counts[2:nrow(gene_counts), ]


gene_counts <- gene_counts[, c(young, old)]
coldata <- data.frame(condition = factor(rep(c('young', 'old'), each = 5),
                                          levels = c('young', 'old')))


gene_counts <- na.omit(gene_counts)
#计算差异倍数列表
dds <- DESeqDataSetFromMatrix(countData = gene_counts, colData = coldata, design = ~condition)
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)

#old 相较于 young
res <- results(dds1, contrast = c('condition', 'old', 'young'))

res

res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)

write.table(res1, "old5_vs_young5_DESeq2.txt", col.names = NA, sep = "\t", quote = FALSE)


#PCA
suppressMessages(library("ggplot2"))
vsd <- vst(dds1, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition", "sizeFactor"), returnData=TRUE)
percentVar <- round(100*attr(pcaData, "percentVar"))
png("PCA.png", width = 2500, height = 2500, res = 284)
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1:", percentVar[1],"% variance")) +
  ylab(paste0("PC2:", percentVar[2],"% variance")) +
  coord_fixed()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()
#write.table(as.matrix(c(old, young)), "samples_id.txt",quote = F)



#火山图
res1 <- na.omit(res1)
res1$logp <- -log10(res1$padj)
res1 <- data.frame(res1)
res1$type <- ifelse(res1$padj < 0.05,ifelse(abs(res1$log2FoldChange) > 1, ifelse(res1$log2FoldChange < -1,'down','up'),'noSig'),'noSig')
write.table(res1, "old5_vs_young5_DESeq2.txt", sep = "\t", quote = FALSE)
table(res1$type)
png("volcano.png", width = 2500, height = 2500, res = 284)
ggplot(res1, aes(x = log2FoldChange, y = -log10(padj)))+
  geom_point(aes(color = type), alpha = 0.5, size = 2)+
  theme_bw(base_size = 16) + 
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), axis.text = element_text(color = 'black'),
        axis.title = element_text(color = 'black'))+
  scale_color_manual(name='', values = c('up'='#DA1212','noSig'='grey','down'='#3E7C17'), label = c('up'='up (1710)','noSig'='noSig (18682)','down'='down (126)')) + 
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  scale_x_continuous(breaks = c(-6,-4,-2,-1,0,1,2,4,6)) 
dev.off()
#明显差异很大的三个基因
#上调
rownames(res1[which(abs(res1$log2FoldChange)>10 & res1$log2FoldChange>0), ])
#下调
rownames(res1[which(abs(res1$log2FoldChange)>10 & res1$log2FoldChange<0), ])
























