rm(list = ls())  
options(stringsAsFactors = F)
options(scipen = 5)

################################ 1  data.frame: expr cpm repClass
#load counts
filelist<-dir("coding_genes/",pattern = "*simple")
filelist <- paste("coding_genes/",filelist,sep="")
file1.count<-read.table(filelist[1],header = T,sep = "\t")

count<-matrix(ncol = length(filelist)+1,nrow = nrow(file1.count))
count<-as.data.frame(count)
count[,1:2]<-file1.count[,c(1,3)]

for(i in 2:length(filelist)){
  file.count<-read.table(filelist[i],header = T,sep = "\t")
  print(filelist[i])
  count[,i+1]<-file.count[,3]
}
head(count)   

colnames(count) <- c("Geneid",filelist)
head(count)


colnames(count) <- gsub(pattern="coding_genes/",replacement = "",colnames(count))
colnames(count) <- gsub(pattern=".coding_genes.count.simple",replacement = "",colnames(count))
head(count)

rownames(count) <- count$Geneid
expr <- count[,2:ncol(count)]
colnames(expr) <- gsub(pattern="_rep",replacement = "-",colnames(expr))
head(expr)
################################ out
suppressMessages(library(DESeq2))
run_DESeq2<-function(exprSet,group_list,contrast_list){
  (colData <- data.frame(row.names=colnames(exprSet),
                         group_list=group_list))
  dds <- DESeqDataSetFromMatrix(countData = exprSet,
                                colData = colData,
                                design = ~ group_list)
  # 这里需要修改
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

#######################有重复样本，用Deseq2做差异分析
exprSet<-expr[,c("mmu_Eed_RNASeq-1","mmu_Eed_RNASeq-2","mESC_wt_RNASeq-1","mESC_wt_RNASeq-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("mmu_Eed_RNASeq","mESC_wt_RNASeq")
Eed_ko.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)

exprSet<-expr[,c("mmu_Ezh2_RNASeq-1","mmu_Ezh2_RNASeq-2","mESC_wt_RNASeq-1","mESC_wt_RNASeq-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("mmu_Ezh2_RNASeq","mESC_wt_RNASeq")
Ezh2_ko.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)

exprSet<-expr[,c("suz12_ko-1","suz12_ko-2","suz12_ko-3","wt_for_suz12-1","wt_for_suz12-2","wt_for_suz12-3")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("suz12_ko","wt_for_suz12")
suz12_ko.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)


DEG_merge<-c(rownames(Eed_ko.vs.WT),rownames(Ezh2_ko.vs.WT),rownames(suz12_ko.vs.WT))

DEG_all<-names(table(DEG_merge))[table(DEG_merge)==3]

DEG_all_FC<-data.frame(Eed_ko.vs.WT[DEG_all,"log2FoldChange"],
                       Ezh2_ko.vs.WT[DEG_all,"log2FoldChange"],
                       suz12_ko.vs.WT[DEG_all,"log2FoldChange"])

rownames(DEG_all_FC)<-DEG_all
colnames(DEG_all_FC)<-c("Eed_ko","Ezh2_ko","suz12_ko")
DEG_all_FC<-na.omit(DEG_all_FC)

write.csv(DEG_all_FC,file = "gene_fc.csv")




