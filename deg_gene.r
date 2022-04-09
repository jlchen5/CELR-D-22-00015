rm(list = ls())
library(DESeq2)
library(ggplot2)

#load counts
filelist<-dir("./",pattern = "*txt")
filelist <- paste("./",filelist,sep="")
file1.count<-read.table(filelist[1],header = T,sep = "\t")


count<-matrix(ncol = length(filelist)+1,nrow = nrow(file1.count))
count<-as.data.frame(count)
count[,1:2]<-file1.count[,c(1,2)]

for(i in 2:length(filelist)){
  file.count<-read.table(filelist[i],header = T,sep = "\t")
  print(filelist[i])
  count[,i+1]<-file.count[,2]
}
head(count)
colnames(count) <- c("ensid","WT_suz12-1","Suz12KO-1","WT_suz12-2",
                     "Suz12KO-2","WT_suz12-3","Suz12KO-3",)

#id conversion
library(org.Mm.eg.db)
# 1. Convert from ensembl.gene to gene.symbol
ensembl.genes <- as.vector(count$ensid)
annots <- select(org.Mm.eg.db, ensembl.genes, 
                 columns="SYMBOL", keytype="ENSEMBL")



rownames(count) <- count[,1]
count <- count[,-1]
head(count)


run_DESeq2<-function(exprSet,group_list,contrast_list){
  (colData <- data.frame(row.names=colnames(exprSet),
                         group_list=group_list))
  dds <- DESeqDataSetFromMatrix(countData = exprSet,
                                colData = colData,
                                design = ~ group_list)
  dds <- DESeq(dds)
  res <- results(dds, 
                 contrast=c("group_list",contrast_list))
  resOrdered <- res[order(res$padj),]
  print(mcols(resOrdered))
  DEG =as.data.frame(resOrdered)
  DESeq2_DEG = na.omit(DEG)
  DESeq2_DEG$FoldChange<-2^DESeq2_DEG$log2FoldChange
  nrDEG=DESeq2_DEG[,c(2,7,5:6)]
  colnames(nrDEG)=c('log2FoldChange','FoldChange','pvalue','padj')
  nrDEG
}

exprSet<-count[,c("WT-1","Suz12KO-1","WT-2","Suz12KO-2","WT-3","Suz12KO-3")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Suz12KO","WT")
Suz12KO.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)

exprSet<-count[,c("WT-1","Suz12KO-1","WT-2","Suz12KO-2","WT-3","Suz12KO-3")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Suz12KO","WT")
Suz12KO.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)