# load librarys
library(stringr)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(ComplexHeatmap)
library(ggpubr)


# read gene id and gene name 
genemap <- read.csv("/data2/ref/genome_data/mm10/mm10_vM25_geneid_geneName.txt",
                    sep="\t",
                    header = T, check.names = F, stringsAsFactors = FALSE)
genenameFreq <- data.frame(table(genemap$gene_name))
colnames(genenameFreq) <- c("geneName","Freq")
genenameFreq$geneName <- as.character(genenameFreq$geneName)

# gene name map to >=2 gene id
freq2genes <- genenameFreq$geneName[genenameFreq$Freq>=2]
for (i in 1:length(freq2genes)){
  idx <- which(genemap$gene_name==freq2genes[i])
  genemap$gene_name[idx] <- paste0(genemap$gene_id[idx],"_",genemap$gene_name[idx])
}



# set work dir
setwd("~/projects/chemi-CATI/")

files_list <- list.files(pattern = ".rsem.genes.results")
samples_list <- str_split(files_list,pattern = "\\.",simplify = T)[,1]

samples_mtx_list <- lapply(files_list, function(x){
  mtx <- read.csv(x,header = T,sep="\t",check.names = F, stringsAsFactors = F)
  return(mtx)
})
names(samples_mtx_list) <- samples_list


# TPM mat
tpm_mtx <- data.frame(lapply(samples_mtx_list,function(x)  x$TPM))
row.names(tpm_mtx) <- samples_mtx_list[[1]]$gene_id
row.names(tpm_mtx) <- genemap$gene_name[match(row.names(tpm_mtx) , genemap$gene_id)]
colnames(tpm_mtx) <-  c("Polq1","Polq2","Polq4","C4","C1","C2","C3","Polq3")

# FPKM mat
fpkm_mtx <- data.frame(lapply(samples_mtx_list,function(x)  x$FPKM))
row.names(fpkm_mtx) <- samples_mtx_list[[1]]$gene_id
row.names(fpkm_mtx) <- genemap$gene_name[match(row.names(fpkm_mtx) , genemap$gene_id)]
colnames(fpkm_mtx) <-  c("Polq1","Polq2","Polq4","C4","C1","C2","C3","Polq3")

# expected_counts
expectedCounts_mtx <- data.frame(lapply(samples_mtx_list,function(x)  x$expected_count))
row.names(expectedCounts_mtx) <- samples_mtx_list[[1]]$gene_id
row.names(expectedCounts_mtx) <- genemap$gene_name[match(row.names(expectedCounts_mtx) , genemap$gene_id)]
colnames(expectedCounts_mtx) <-  c("Polq1","Polq2","Polq4","C4","C1","C2","C3","Polq3")


write.table(tpm_mtx, file="merged_TPM_matrix.txt", sep="\t",quote=FALSE)
write.table(fpkm_mtx, file="merged_FPKM_matrix.txt", sep="\t",quote=FALSE)
write.table(expectedCounts_mtx, file="merged_expectedCounts_matrix.txt", sep="\t",quote=FALSE)


rm(list=ls())

expectedCounts_mtx <- read.csv("merged_expectedCounts_matrix.txt",sep="\t",
                               check.names = F, header = T)
tpm_mtx <- read.csv("merged_TPM_matrix.txt",sep="\t",check.names = F, header = T)
fpkm_mtx <- read.csv("merged_FPKM_matrix.txt",sep="\t",check.names = F, header = T)

expectedCounts_mtx <- expectedCounts_mtx[rowSums(expectedCounts_mtx)>0,]

# filter genes
keep_rows1 <- apply(expectedCounts_mtx[,c(1,2,3,8)],1,function(x){
  res=FALSE
  if(sum(x>0)>=2){
    res=TRUE
  }
  return(res)
})
table(keep_rows1)

keep_rows2 <- apply(expectedCounts_mtx[,c(4:7)],1,function(x){
  res=FALSE
  if(sum(x>0)>=2){
    res=TRUE
  }
  return(res)
})
table(keep_rows2)

keep_rows <- (keep_rows2 | keep_rows1)
counts_mtx <- expectedCounts_mtx[keep_rows,]

coldata <- data.frame(row.names = colnames(counts_mtx),
                      sample = colnames(counts_mtx),
                      condition = c("Polq","Polq","Polq","Ctrl",
                                    "Ctrl","Ctrl","Ctrl","Polq"))
dds <- DESeqDataSetFromMatrix(countData = round(counts_mtx),
                              colData = coldata,
                              design = ~condition)

vstd <- vst(dds, blind=TRUE)
pcaData <- plotPCA(vstd, intgroup = c("condition"), ntop = 500,returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p <- ggplot(pcaData, aes(PC1, PC2,
                    color = coldata$condition,
                    label = coldata$sample)) +
  geom_point(size=3) + #coord_fixed() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  geom_text_repel(
    nudge_x = .15,
    box.padding = 0.5,
    nudge_y = 1,
    segment.curvature = -0.1,
    segment.ncp = 3,
    segment.angle = 20,
    size = 3)+
  theme_bw()
ggsave(p, 
       file="./output/deseq2_vst_top500_pca.pdf",
       height = 3.98, width = 6.18)

dds_rlog <- rlog(dds, blind=TRUE)
pcaData <- plotPCA(dds_rlog, intgroup = c("condition"), ntop = 500,returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p <- ggplot(pcaData, aes(PC1, PC2,
                         color = coldata$condition,
                         label = coldata$sample)) +
  geom_point(size=3) + #coord_fixed() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  geom_text_repel(
    nudge_x = .15,
    box.padding = 0.5,
    nudge_y = 1,
    segment.curvature = -0.1,
    segment.ncp = 3,
    segment.angle = 20,
    size = 3)+
  theme_bw()
ggsave(p, 
       file="./output/deseq2_rlog_top500_pca.pdf",
       height = 3.98, width = 6.18)


pheatmap(cor(assay(vstd),method = "pearson"))
pheatmap(cor(assay(vstd),method = "spearman"))

pheatmap(cor(assay(dds_rlog),method = "pearson"))
pheatmap(cor(assay(dds_rlog),method = "spearman"))

dds <- DESeq(dds) 
results <- results(dds, alpha=0.05, lfcThreshold=1)  
resultsOrdered <- results[order(results$pvalue),] 
summary(results) 
result_data <- merge(as.data.frame(results),
                     as.data.frame(counts(dds, normalized=TRUE)),
                     by="row.names",
                     sort=FALSE)

result_data[which(result_data$log2FoldChange >= 1 & result_data$padj < 0.05),'sig'] <- 'up'
result_data[which(result_data$log2FoldChange <= -1 & result_data$padj < 0.05),'sig'] <- 'down'
result_data[which(abs(result_data$log2FoldChange) <= 1 | result_data$padj >= 0.05),'sig'] <- 'none'

write.table(result_data,
            file="./output/deseq2_dds_degs_ALL.xls",
            sep="\t",quote=FALSE, row.names = FALSE)


p <- ggplot(data = result_data, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(size = 1) +  
  scale_color_manual(values = c('red', 'gray', '#0000FF'), limits = c('up', 'none', 'down')) +  
  labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = 'Polq vs Ctrl', color = '') +  
  theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') + 
  geom_hline(yintercept = (-log10(0.05)), lty = 3, color = 'black') 
p
ggsave(p, file="./output/deseq2_Volcano.pdf",
       height = 4.68, width = 6.19)

res_select <- subset(result_data, sig %in% c('up', 'down'))
write.table(res_select,
            file="./output/deseq2_dds_degs_Sig.xls",
            sep="\t",quote=FALSE, row.names = FALSE)

res_select_tpm_mtx <- tpm_mtx[res_select$Row.names,]
pheatmap(log10(res_select_tpm_mtx+1))


Polq_fpkm <- cbind(coldata,fpkm=as.numeric(fpkm_mtx['Polq',]))
Polq_fpkm$fpkm <- log10(Polq_fpkm$fpkm+1)
p <- ggplot(Polq_fpkm, aes(x=condition, y=fpkm,fill=condition))+geom_boxplot(width=0.5) +
  stat_compare_means() +theme_classic() +
  theme(legend.position = "na",
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 12)) +
  xlab(label = "") + ylab(label="log10 (FPMK+1)")
ggsave(p, 
       file="./output/Polq_fpkm_boxplot.pdf",
       height = 3.98, width = 4.18)
