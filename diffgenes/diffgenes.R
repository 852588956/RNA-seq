library(DESeq2)


#加载数据"
group <- read.table("../GSVA/Placebo.txt", header = T, row.names = 1, sep = "\t")
young <- rownames(group)[which(group$age == "25-37 years")]
old <- rownames(group)[which(group$age == "60-79 years")]

counts <- read.table("../data/GSE141595_counts.txt", header = T, row.names = 1, sep = "\t")
#counts <- read.table("../data/counts.txt", header = T, row.names = 1)
counts <- counts[, c(young, old)]

coldata <- data.frame(condition = factor(rep(c('young', 'old'), each = 5),
                                          levels = c('young', 'old')))


counts <- na.omit(counts)
#计算差异倍数列表
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~condition)
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)

#old 相较于 young
res <- results(dds1, contrast = c('condition', 'old', 'young'))

res

res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
write.table(res1, "old_vs_young_DESeq2.txt", col.names = NA, sep = "\t", quote = FALSE)


#PCA
suppressMessages(library("ggplot2"))
vsd <- vst(dds1, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition", "sizeFactor"), returnData=TRUE)
percentVar <- round(100*attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1:", percentVar[1],"% variance")) +
  ylab(paste0("PC2:", percentVar[2],"% variance")) +
  coord_fixed()

#write.table(as.matrix(c(old, young)), "samples_id.txt",quote = F)
