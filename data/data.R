library(GEOquery)

#将所有数据处理成两个文件，1是所有样本的表达矩阵，去掉所有有缺失的基因
#                          2是所有样本属性的矩阵，包括其来源和种类（年轻or老年or处理过的老年）


#GSE141595数据处理 
my_id <- "GSE141595"
gse <- getGEO(my_id,destdir = './', GSEMatrix = TRUE)
length(gse)
gse<-gse[[1]]
##样本临床信息
p <- pData(gse) ## print the sample information
f <- fData(gse) ## print the gene annotation
e <- exprs(gse) ## print the expression data
##第一个count文件
data1 <- read.delim("GSE141595_Osteocyte_enriched_bone_biopsies_gene_count.csv",
                    header = T, row.names = 2, sep = ",")
samples <- colnames(data1)[6:ncol(data1)]
gid <- c()
title <- c()
for(i in 1:length(samples)){
  title <- c(title, strsplit(samples[i], "_")[[1]][2])
}
for(i in 1:length(title)){
  gid <- c(gid, rownames(p)[which(p$title == title[i])])
}
colnames(data1)[6:length(data1)] <- gid


##第二个count文件
data2 <- read.delim("GSE141595_Raw_Data_GSM6246351-GSM6246365_June_13_2022.csv",
                    header = T, row.names = 2, sep = ",")
colnames(data2)[6:ncol(data2)] <- rownames(p)[31: nrow(p)]


##基因交集
gene <- intersect(rownames(data1), rownames(data2))

GSE141595 <- cbind(data1[gene, ], data2[gene, 6:ncol(data2)])

types <- p[, c("age:ch1", "treatment:ch1")]
types$ori <- "GSE141595"

GSE141595 <- GSE141595[!duplicated(GSE141595$GeneName), ]
rownames(GSE141595) <- GSE141595$GeneName



counts <- na.omit(GSE141595)

#算tpm值
genes_length <- abs(counts$Stop - counts$Start) / 1000
tpm <- counts[, 6:ncol(counts)]
for (i in 1:ncol(tpm)) {
  tpm[, i] <- tpm[, i] / genes_length
  all_reads <- sum(tpm[, i]) 
  tpm[, i] <- tpm[, i] / all_reads
}




write.table(tpm, "tpm.txt", quote = FALSE)



#整理分类
colnames(types) <- c("age", "treatment", "source")


write.table(types, "types.txt", quote = FALSE, sep = "\t")







