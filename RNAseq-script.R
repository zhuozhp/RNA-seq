# We get CSV file "gene_count.csv" from StringTie, and we performed the following code.
# install the necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

library(BiocInstaller) 
BiocManager::install(c('DESeq2','amap','BiocParallel'),
                     ask = F,update = F)

install.packages("RColorBrewer")
install.packages("gplots")
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("dplyr")
install.packages("hexbin")
install.packages("FactoMineR")
install.packages("corrplot")
install.packages("factoextra")

library("DESeq2")
library("RColorBrewer")
library("gplots")
library("amap")
library("ggplot2")
library("BiocParallel")
library("pheatmap")
library("dplyr")
library("hexbin")
library("FactoMineR")
library("corrplot")
library("factoextra")

# data import
data <- read.csv("gene_count.csv", header=T, row.names="gene_id", com='', quote='')
data <- data[rowSums(data)>2,]
data <- data[apply(data, 1, var)!=0,]
head(data, 2)
dim(data)

# conditions infomation setting
conditions <- factor(c(rep("BLM-0D_WT",3), rep("BLM-REC-1D_WT",3), rep("BLM-REC-3D_WT",3),
                       rep("BLM-0D_MT",3), rep("BLM-REC-1D_MT",3), rep("BLM-REC-3D_MT",3)), 
                     levels = c("BLM-0D_WT","BLM-REC-1D_WT","BLM-REC-3D_WT",
                                "BLM-0D_MT","BLM-REC-1D_MT","BLM-REC-3D_MT"))
colData  <- data.frame(row.names=colnames(data), conditions)
sample  <- data.frame(row.names=colnames(data), conditions)
sample
sample <- sample[match(colnames(data), rownames(sample)),, drop=F]
sample_rowname <- rownames(sample)

# DESeq2 comparison
ddsMat <- DESeqDataSetFromMatrix(countData = data, 
                                 colData = sample,  
                                 design= ~ conditions)
dds <- DESeq(ddsMat)
rld <- rlog(dds, blind = FALSE)
rlogMat <- assay(rld)

# sample rlog boxplot
n.sample <- ncol(data)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
pdf(file = "rlog_boxplot_of_each_sample.pdf")
par(cex=0.8, cex.main=2, cex.lab=1)
par(pin=c(4,4),las=2)
boxplot(rlogMat, col = cols,main="expression value")
dev.off()

normalized_counts <- DESeq2::counts(dds, normalized=TRUE)
normalized_counts_mad <- apply(normalized_counts, 1, mad)
write.csv(normalized_counts, file = "normalized_counts.csv")

rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]
write.csv(normalized_counts, file = "normalized_counts_rlog.csv")

# pearson correlation
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(50)
pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
hc <- hcluster(t(rlogMat), method="pearson")
pdf("sample_rlog_pearson_correlation.pdf", pointsize=8)
heatmap.2(pearson_cor, Rowv=as.dendrogram(hc), symm=T, 
          trace="none", col=hmcol, margins=c(12,12),
          main="The pearson correlation of each sample")
dev.off()

# pca plot
pca_data <- plotPCA(rld, intgroup = c("conditions"), 
                    returnData = T)
percentVar <- round(100 * attr(pca_data, "percentVar"))
ggplot(pca_data, aes(PC1, PC2, color = conditions)) +
  geom_point(size = 2.5) +
  xlab(paste("PC1:", percentVar[1], "% variance")) +
  ylab(paste("PC2:", percentVar[2], "% variance")) +
  labs(title="PCA plot") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", 
        legend.title = element_blank(), panel.background = element_blank(),
        panel.grid = element_blank(), legend.key = element_blank(), 
        panel.border = element_rect(colour = "gray", fill = NA))
ggsave("ggplot2_PCA_plot.pdf", width = 6, height = 5, units = "in")


# two samples comparison, take "BLM-REC-1D_WT" vs "BLM-0D_WT" as an example.
sampleA = "BLM-REC-1D_WT"
sampleB = "BLM-0D_WT"
res <- results(dds, contrast = c("conditions", 
                                 sampleA, sampleB), test = "Wald")

{
  baseA <- DESeq2::counts(dds, normalized=TRUE)[, colData$conditions == sampleA]
  
  if (is.vector(baseA)){
    baseMeanA <- as.data.frame(baseA)
  } else {
    baseMeanA <- as.data.frame(rowMeans(baseA))
  }
  colnames(baseMeanA) <- sampleA
  head(baseMeanA)
  
  baseB <- DESeq2::counts(dds, normalized=TRUE)[, colData$conditions == sampleB]
  if (is.vector(baseB)){
    baseMeanB <- as.data.frame(baseB)
  } else {
    baseMeanB <- as.data.frame(ID=rownames(baseB), rowMeans(baseB))
  }
  colnames(baseMeanB) <- sampleB
  head(baseMeanB)
  result <- cbind(baseA, baseB, baseMeanA, baseMeanB, as.data.frame(res))
  head(result, 2)
  result$baseMean <- rowMeans(cbind(baseA, baseB))
  result$padj[is.na(result$padj)] <- 1
  result <- na.omit(result)
  result$ID <- rownames(result)
  result <- result[order(result$pvalue), ]
}
comp1 <- paste("Count_matrix_",sampleA, "_vs_", sampleB, ".csv", sep="")
write.csv(result, file = comp1)


# volcano plot
library(ggrepel)
result_1 <- na.omit(as.data.frame(res))
result_1$ID <- rownames(result_1)
rownames(result_1) <- result_1$ID
result_1 <- result_1[,-1]
head(result_1)
volcano <- ggplot(data,aes(logFC,-1*log10(pvalue)))
volcano$threshold <- as.factor(ifelse(result_1$pvalue < 0.05 & 
                                        abs(result_1$log2FoldChange) >= 0.6, 
                                      ifelse(result_1$log2FoldChange >= 0.6, 
                                             'Up','Down'),'No Significant'))
gene_list <- read.csv("select.csv", header = T)
head(gene_list)
gene_list <- gene_list[,1:2]
colnames(gene_list) <- c("ID", "Gene.Name")
res_int <- as.data.frame(res[gene_list$ID,])
dim(res_int)
res_int$ID <- rownames(res_int)
res_int <- merge(res_int, anno, by="ID")
rownames(result_1) <- result_1$ID
head(res_int)
  
res_int_up <- subset(res_int, res_int$log2FoldChange>=0.6 & res_int$pvalue<=0.05)
res_int_dw <- subset(res_int, res_int$log2FoldChange<=-0.6 & res_int$pvalue<=0.05)
res_int_ns <- subset(res_int, res_int$pvalue > 0.05 | res_int$log2FoldChange<0.6 & res_int$log2FoldChange>-0.6)

  
threshold <- volcano$threshold
p <- ggplot(data = result_1, aes(x = log2FoldChange, y = -log10(pvalue), 
                                 colour=threshold)) +
  geom_point(data = result_1, aes(x = log2FoldChange, y = -log10(pvalue), 
                                  colour=threshold), alpha=0.6, size=0.95) +
  scale_color_manual(values=c("#2E2EFE","gray80","#FE2E2E")) +
  xlim(c(-10, 10)) + ylim(c(0, 15)) + 
  geom_vline(xintercept=c(-0.6,0.6),lty=4,col="darkgray",lwd=0.3) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="darkgray",lwd=0.3) +
  labs(x="log2FoldChange", y="-log10(pvalue)", cex=2.5, title=NULL) +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="none", legend.title = element_blank(),
        legend.key = element_blank(),
        panel.grid = element_blank(), 
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8))+
  geom_text_repel(data=res_int_up, aes(label=res_int_up$ID),
                  col="darkred",alpha = 1, fontface = "bold", cex=3, 
                  arrow = arrow(length=unit(0.005, "npc")),segment.alpha = 0.95,
                  force = 0.2, max.iter = 3e3, 
                  box.padding=unit(0.6, "lines"), point.padding=unit(0.1, "lines"), 
                  segment.colour = "darkred", segment.size = 0.4, nudge_x = 0.59, 
                  min.segment.length = 0.1, nudge_y = -0.3) + 
  geom_text_repel(data=res_int_dw, aes(label=res_int_dw$ID),
                  col="darkblue",alpha = 1, fontface = "bold", cex=3,
                  arrow = arrow(length=unit(0.005, "npc")),force = 0.2, max.iter = 3e3, 
                  box.padding=unit(0.6, "lines"), point.padding=unit(0.1, "lines"), 
                  segment.colour = "darkblue", segment.size = 0.4,
                  segment.alpha = 1,
                  min.segment.length = 0.2, nudge_y = 0.1, nudge_x = -0.69) +
  geom_text_repel(data=res_int_ns, aes(label=res_int_ns$ID),
                  col="black",alpha = 0.4, fontface = "bold", cex=3,
                  arrow = arrow(length=unit(0.005, "npc")),force = 0.5, max.iter = 3e3, 
                  box.padding=unit(0.3, "lines"), point.padding=unit(0.1, "lines"), 
                  segment.colour = "black", segment.size = 0.4,
                  segment.alpha = 0.4,
                  min.segment.length = 0.2, nudge_x = -0.1, nudge_y = -0.1)
  
p
  
file_base <- paste("vocano_plot_",sampleA, " vs ", sampleB, ".pdf", sep="")
ggsave(file = file_base, width = 5.4 , height = 6.1)

# install necessary packages
BiocManager::install(c('KEGG.db','clusterProfiler','org.Dm.eg.db'),
                     ask = F,update = F)
install.packages("car")

library("KEGG.db")
library("clusterProfiler")
library("car")
library("ggplot2")
library("org.Dm.eg.db")

pthres_1 <- 0.05
logFC_1 <- 0.7
res <- results(dds, contrast = c("conditions", 
                                 sampleA, sampleB), test = "Wald")
res1 <- as.data.frame(res1)
res1 <- na.omit(res1)
up <- as.vector(rownames(res1[res1$pvalue<=pthres_1 & res1$log2FoldChange>=logFC_1,]))
length(up)
dw <- as.vector(rownames(res1[res1$pvalue<=pthres_1 & res1$log2FoldChange<= -logFC_1,]))
length(dw)
gene <- as.data.frame(up)
head(gene)
gene$status <- "up"
colnames(gene) <- c("ID", "status")
gene_dw <- as.data.frame(dw)
head(gene_dw)
gene_dw$status <- "dw"
colnames(gene_dw) <- c("ID", "status")
gene <- rbind(gene, gene_dw)
rm(gene_dw)

# pathway bubble plot
ga <- list("REC-1D_WT_up" = up, "REC-1D_WT_dw"= dw)
X <- enrichDAVID(ga,idType = "FLYBASE_GENE_ID", minGSSize = 1,
                 maxGSSize = 500, annotation = "KEGG_PATHWAY", 
                 pvalueCutoff = 0.8, pAdjustMethod = "BH", qvalueCutoff = 1,
                 species = NA, david.user = "zhuozhp@mail2.sysu.edu.cn")
file_name <- paste("kegg_enrichment_", sampleA, "_vs_", sampleB, "_logFC",
                   logFC_1, "_p", pthres_1, sep = "")
pdf(file = paste(file_name, ".pdf", sep = ""))
p <- dotplot(X, showCategory=15, font.size = 10) + 
  scale_color_continuous(low='red', high='green')
p
dev.off()
