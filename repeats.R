rm(list = ls())  
options(stringsAsFactors = F)
options(scipen = 5)

################################ 1  data.frame: expr cpm repClass
filelist<-dir("repeats/",pattern = "*simple")
filelist <- paste("repeats/",filelist,sep="")
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

expr <- count[,-1]
rownames(expr) <- count[,1]
colnames(expr)<-c("Control_for_Group2-1","Control_for_Group2-2",
                  "Control_for_Group3","Control_for_Group4",
                  "Control_for_Group1-1","Control_for_Group1-2",
                  "D0_sgDnmt1-1","D0_sgDnmt1-2","D0_sgMyc-1","D0_sgMyc-2",
                  "D1_sgDnmt1-1","D1_sgDnmt1-2","D1_sgMyc-1","D1_sgMyc-2",
                  "Group1_shSae1-1","Group1_shSae1-2",
                  "Group1_shSenp6-1","Group1_shSenp6-2",
                  "Group1_shSumo2-1","Group1_shSumo2-2",
                  "Group1_shTrim28","Group1_shUba2-1",
                  "Group1_shUba2-2","Group1_shUbe2i-1",
                  "Group1_shUbe2i-2","Group2_shChaf1a-1",
                  "Group2_shChaf1a-2","Group2_shChaf1b-1",
                  "Group2_shChaf1b-2","Group3_shAtf7ip",
                  "Group4_shEset","KAP1_KO-0",
                  "KAP1_WT-0","KAP1_WT-1","KAP1_WT-2","KAP1_WT-3",
                  "KAP1_KO-1","KAP1_KO-2","KAP1_KO-3","LIN28_WT_ES-1",
                  "LIN28_WT_ES-2","LIN28_WT_ES-3","Lin28a_ko-2","Lin28a_ko-3",
                  "Lin28b_ko-1","Lin28b_ko-2","Lin28b_ko-3","LINE1_siControl-1",
                  "LINE1_siControl-2","LINE1_siControl-3",
                  "LINE1_siDux-1","LINE1_siDux-2","LINE1_siDux-3",
                  "Max_KD-1","Max_KD-2","mESC_wt_RNASeq-1","mESC_wt_RNASeq-2",
                  "Mga_KD-1","Mga_KD-2","Mga_WT-1","Mga_WT-2",
                  "mmu_Eed_RNASeq-1","mmu_Eed_RNASeq-2",
                  "mmu_Ezh2_RNASeq-1","mmu_Ezh2_RNASeq-2",
                  "Pcgf1_KO-1","Pcgf1_KO-2","Pcgf124_KO-1","Pcgf124_KO-2",
                  "Pcgf24_KO-1","Pcgf24_KO-2","Pcgf35_KO-1","Pcgf35_KO-2",
                  "Pcgf356_KO-1","Pcgf356_KO-2","Pcgf6_KD-1","Pcgf6_KD-2",
                  "Pcgf6_KO-1","Pcgf6_KO-2",
                  "Rif1_KO-1","Rif1_KO-2",
                  "Rif1_WT-1","Rif1_WT-2",
                  "RNF2_KD-1","RNF2_KD-2",
                  "Smchd1_WT_1_WT1-1","Smchd1_WT_1_WT1-2",
                  "Smchd1_WT_2_WT2-1","Smchd1_WT_2_WT2-2",
                  "Smchd1_WT_3_WT3-1","Smchd1_WT_3_WT3-2",
                  "Smchd1_KO1-1","Smchd1_KO1-2",
                  "Smchd1_KO2-1","Smchd1_KO2-2",
                  "Smchd1_KO3-1","Smchd1_KO3-2",
                  "suz12_ko-1","suz12_ko-2","suz12_ko-3",
                  "WT-1","WT-2",
                  "wt_for_suz12-1","wt_for_suz12-2","wt_for_suz12-3",
                  "Ythdc1_cKO_4OHT-1","Ythdc1_cKO_4OHT-2",
                  "Ythdc1_cKO_DMSO-1","Ythdc1_cKO_DMSO-2")

#lib <- apply(expr,2,sum) # colculate library size
#cpm <- t(t(expr)/lib)*10^6 # colculate CPM

#check<-cpm[,c("WT-1","MaxKO-2","WT-2","Rif1KO-1","Rif1KO-2","RNF2KO-1","Pcgf6KO-1","Pcgf6KO-2","RNF2KO-2","MaxKO-1")]
#check[c("Max","Rif1","Rnf2","Pcgf6"),]


################################ out
outDir<-"repeats.out/"
if(! file_test("-d", outDir)){
  dir.create(outDir)
}


suppressMessages(library(DESeq2))
run_DESeq2<-function(exprSet,group_list,contrast_list){
  (colData <- data.frame(row.names=colnames(exprSet),
                         group_list=group_list))
  dds <- DESeqDataSetFromMatrix(countData = exprSet,
                                colData = colData,
                                design = ~ group_list)
  # 这里需要修改
  dds <- DESeq(dds)
  res <- results(dds,contrast=c("group_list",contrast_list))
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
exprSet<-expr[,c("D0_sgDnmt1-1","D0_sgDnmt1-2","D1_sgDnmt1-1","D1_sgDnmt1-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("D1_sgDnmt1","D0_sgDnmt1")
sgDnmt1.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
sgDnmt1.vs.WT_pvalue0.05<-subset(sgDnmt1.vs.WT,pvalue<0.05)

exprSet<-expr[,c("D0_sgMyc-1","D0_sgMyc-2","D1_sgMyc-1","D1_sgMyc-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("D1_sgMyc","D0_sgMyc")
sgMyc.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
sgMyc.vs.WT_pvalue0.05<-subset(sgMyc.vs.WT,pvalue<0.05)

exprSet<-expr[,c("Group1_shSae1-1","Group1_shSae1-2","Control_for_Group1-1","Control_for_Group1-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Group1_shSae1","Control_for_Group1")
shSae1.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
shSae1.vs.WT_pvalue0.05<-subset(shSae1.vs.WT,pvalue<0.05)

exprSet<-expr[,c("Group1_shSenp6-1","Group1_shSenp6-2","Control_for_Group1-1","Control_for_Group1-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Group1_shSenp6","Control_for_Group1")
shSenp6.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
shSenp6.vs.WT_pvalue0.05<-subset(shSenp6.vs.WT,pvalue<0.05)

exprSet<-expr[,c("Group1_shSumo2-1","Group1_shSumo2-2","Control_for_Group1-1","Control_for_Group1-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Group1_shSumo2","Control_for_Group1")
shSumo2.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
shSumo2.vs.WT_pvalue0.05<-subset(shSumo2.vs.WT,pvalue<0.05)

exprSet<-expr[,c("Group1_shUba2-1","Group1_shUba2-2","Control_for_Group1-1","Control_for_Group1-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Group1_shUba2","Control_for_Group1")
shUba2.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
shUba2.vs.WT_pvalue0.05<-subset(shUba2.vs.WT,pvalue<0.05)

exprSet<-expr[,c("Group1_shUbe2i-1","Group1_shUbe2i-2","Control_for_Group1-1","Control_for_Group1-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Group1_shUbe2i","Control_for_Group1")
shUbe2i.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
shUbe2i.vs.WT_pvalue0.05<-subset(shUbe2i.vs.WT,pvalue<0.05)

exprSet<-expr[,c("Group2_shChaf1a-1","Group2_shChaf1a-2","Control_for_Group2-1","Control_for_Group2-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Group2_shChaf1a","Control_for_Group2")
shChaf1a.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
shChaf1a.vs.WT_pvalue0.05<-subset(shChaf1a.vs.WT,pvalue<0.05)

exprSet<-expr[,c("Group2_shChaf1b-1","Group2_shChaf1b-2","Control_for_Group2-1","Control_for_Group2-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Group2_shChaf1b","Control_for_Group2")
shChaf1b.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
shChaf1b.vs.WT_pvalue0.05<-subset(shChaf1b.vs.WT,pvalue<0.05)

exprSet<-expr[,c("KAP1_KO-0","KAP1_KO-1","KAP1_KO-2","KAP1_KO-3","KAP1_WT-0","KAP1_WT-1","KAP1_WT-2","KAP1_WT-3")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("KAP1_KO","KAP1_WT")
KAP1_KO.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
KAP1_KO.vs.WT_pvalue0.05<-subset(KAP1_KO.vs.WT,pvalue<0.05)

exprSet<-expr[,c("Lin28a_ko-2","Lin28a_ko-3","LIN28_WT_ES-2","LIN28_WT_ES-3")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Lin28a_ko","LIN28_WT_ES")
Lin28a_ko.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
Lin28a_ko.vs.WT_pvalue0.05<-subset(Lin28a_ko.vs.WT,pvalue<0.05)

exprSet<-expr[,c("Lin28b_ko-1","Lin28b_ko-2","Lin28b_ko-3","LIN28_WT_ES-1", "LIN28_WT_ES-2","LIN28_WT_ES-3")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Lin28b_ko","LIN28_WT_ES")
Lin28b_ko.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
Lin28b_ko.vs.WT_pvalue0.05<-subset(Lin28b_ko.vs.WT,pvalue<0.05)

exprSet<-expr[,c("LINE1_siDux-1","LINE1_siDux-2","LINE1_siDux-3","LINE1_siControl-1","LINE1_siControl-2","LINE1_siControl-3")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("LINE1_siDux","LINE1_siControl")
LINE1_siDux.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
LINE1_siDux.vs.WT_pvalue0.05<-subset(LINE1_siDux.vs.WT,pvalue<0.05)

exprSet<-expr[,c("Max_KD-1","Max_KD-2","Rif1_WT-1","Rif1_WT-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Max_KD","Rif1_WT")
Max_KD.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
Max_KD.vs.WT_pvalue0.05<-subset(Max_KD.vs.WT,pvalue<0.05)

exprSet<-expr[,c("Mga_KD-1","Mga_KD-2","Mga_WT-1","Mga_WT-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Mga_KD","Mga_WT")
Mga_KD.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
Mga_KD.vs.WT_pvalue0.05<-subset(Mga_KD.vs.WT,pvalue<0.05)

exprSet<-expr[,c("Pcgf124_KO-1","Pcgf124_KO-2","WT-1","WT-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Pcgf124_KO","WT")
Pcgf124_KO.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
Pcgf124_KO.vs.WT_pvalue0.05<-subset(Pcgf124_KO.vs.WT,pvalue<0.05)

exprSet<-expr[,c("Pcgf1_KO-1","Pcgf1_KO-2","WT-1","WT-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Pcgf1_KO","WT")
Pcgf1_KO.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
Pcgf1_KO.vs.WT_pvalue0.05<-subset(Pcgf1_KO.vs.WT,pvalue<0.05)

exprSet<-expr[,c("Pcgf24_KO-1","Pcgf24_KO-2","WT-1","WT-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Pcgf24_KO","WT")
Pcgf24_KO.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
Pcgf24_KO.vs.WT_pvalue0.05<-subset(Pcgf24_KO.vs.WT,pvalue<0.05)

exprSet<-expr[,c("Pcgf356_KO-1","Pcgf356_KO-2","WT-1","WT-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Pcgf356_KO","WT")
Pcgf356_KO.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
Pcgf356_KO.vs.WT_pvalue0.05<-subset(Pcgf356_KO.vs.WT,pvalue<0.05)

exprSet<-expr[,c("Pcgf35_KO-1","Pcgf35_KO-2","WT-1","WT-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Pcgf35_KO","WT")
Pcgf35_KO.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
Pcgf35_KO.vs.WT_pvalue0.05<-subset(Pcgf35_KO.vs.WT,pvalue<0.05)

exprSet<-expr[,c("Pcgf6_KO-1","Pcgf6_KO-2","WT-1","WT-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Pcgf6_KO","WT")
Pcgf6_KO.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
Pcgf6_KO.vs.WT_pvalue0.05<-subset(Pcgf6_KO.vs.WT,pvalue<0.05)

exprSet<-expr[,c("Pcgf6_KD-1","Pcgf6_KD-2","Rif1_WT-1","Rif1_WT-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Pcgf6_KD","Rif1_WT")
Pcgf6_KD.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
Pcgf6_KD.vs.WT_pvalue0.05<-subset(Pcgf6_KD.vs.WT,pvalue<0.05)

exprSet<-expr[,c("Rif1_KO-1","Rif1_KO-2","Rif1_WT-1","Rif1_WT-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Rif1_KO","Rif1_WT")
Rif1_KO.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
Rif1_KO.vs.WT_pvalue0.05<-subset(Rif1_KO.vs.WT,pvalue<0.05)

exprSet<-expr[,c("RNF2_KD-1","RNF2_KD-2","Rif1_WT-1","Rif1_WT-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("RNF2_KD","Rif1_WT")
RNF2_KD.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
RNF2_KD.vs.WT_pvalue0.05<-subset(RNF2_KD.vs.WT,pvalue<0.05)

exprSet<-expr[,c("Smchd1_KO1-1","Smchd1_KO1-2","Smchd1_WT_1_WT1-1","Smchd1_WT_1_WT1-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Smchd1_KO1","Smchd1_WT_1_WT1")
Smchd1_KO1.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
Smchd1_KO1.vs.WT_pvalue0.05<-subset(Smchd1_KO1.vs.WT,pvalue<0.05)

exprSet<-expr[,c("Smchd1_KO2-1","Smchd1_KO2-2","Smchd1_WT_2_WT2-1","Smchd1_WT_2_WT2-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Smchd1_KO2","Smchd1_WT_2_WT2")
Smchd1_KO2.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
Smchd1_KO2.vs.WT_pvalue0.05<-subset(Smchd1_KO2.vs.WT,pvalue<0.05)

exprSet<-expr[,c("Smchd1_KO3-1","Smchd1_KO3-2","Smchd1_WT_3_WT3-1","Smchd1_WT_3_WT3-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Smchd1_KO3","Smchd1_WT_3_WT3")
Smchd1_KO3.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
Smchd1_KO3.vs.WT_pvalue0.05<-subset(Smchd1_KO3.vs.WT,pvalue<0.05)

exprSet<-expr[,c("Ythdc1_cKO_4OHT-1","Ythdc1_cKO_4OHT-2","Ythdc1_cKO_DMSO-1","Ythdc1_cKO_DMSO-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("Ythdc1_cKO_4OHT","Ythdc1_cKO_DMSO")
Ythdc1_cKO_4OHT.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
Ythdc1_cKO_4OHT.vs.WT_pvalue0.05<-subset(Ythdc1_cKO_4OHT.vs.WT,pvalue<0.05)

exprSet<-expr[,c("mmu_Eed_RNASeq-1","mmu_Eed_RNASeq-2","mESC_wt_RNASeq-1","mESC_wt_RNASeq-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("mmu_Eed_RNASeq","mESC_wt_RNASeq")
Eed_ko.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
Eed_ko.vs.WT_pvalue0.05<-subset(Eed_ko.vs.WT,pvalue<0.05)

exprSet<-expr[,c("mmu_Ezh2_RNASeq-1","mmu_Ezh2_RNASeq-2","mESC_wt_RNASeq-1","mESC_wt_RNASeq-2")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("mmu_Ezh2_RNASeq","mESC_wt_RNASeq")
Ezh2_ko.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
Ezh2_ko.vs.WT_pvalue0.05<-subset(Ezh2_ko.vs.WT,pvalue<0.05)

exprSet<-expr[,c("suz12_ko-1","suz12_ko-2","suz12_ko-3","wt_for_suz12-1","wt_for_suz12-2","wt_for_suz12-3")]
group_list<-as.factor(sub("-\\d","",colnames(exprSet)))
contrast_list<-c("suz12_ko","wt_for_suz12")
suz12_ko.vs.WT<-run_DESeq2(exprSet,group_list,contrast_list)
suz12_ko.vs.WT_pvalue0.05<-subset(suz12_ko.vs.WT,pvalue<0.05)


DEG_merge<-c(rownames(sgDnmt1.vs.WT),rownames(sgMyc.vs.WT),
             rownames(shSae1.vs.WT),rownames(shSenp6.vs.WT),
             rownames(shSumo2.vs.WT),rownames(shUba2.vs.WT),
             rownames(shUbe2i.vs.WT),rownames(shChaf1a.vs.WT),
             rownames(shChaf1b.vs.WT),rownames(KAP1_KO.vs.WT),
             rownames(Lin28a_ko.vs.WT),rownames(Lin28b_ko.vs.WT),
             rownames(LINE1_siDux.vs.WT),rownames(Max_KD.vs.WT),
             rownames(Mga_KD.vs.WT),rownames(Pcgf124_KO.vs.WT),
             rownames(Pcgf1_KO.vs.WT),rownames(Pcgf24_KO.vs.WT),
             rownames(Pcgf356_KO.vs.WT),rownames(Pcgf35_KO.vs.WT),
             rownames(Pcgf6_KO.vs.WT),rownames(Pcgf6_KD.vs.WT),
             rownames(Rif1_KO.vs.WT),rownames(RNF2_KD.vs.WT),
             rownames(Smchd1_KO1.vs.WT),rownames(Smchd1_KO2.vs.WT),
             rownames(Smchd1_KO3.vs.WT),rownames(Ythdc1_cKO_4OHT.vs.WT),
             rownames(Eed_ko.vs.WT),rownames(Ezh2_ko.vs.WT),
             rownames(suz12_ko.vs.WT))

DEG_all<-names(table(DEG_merge))[table(DEG_merge)==31]

DEG_all_FC<-data.frame(sgDnmt1.vs.WT[DEG_all,"log2FoldChange"],
                       sgMyc.vs.WT[DEG_all,"log2FoldChange"],
                       shSae1.vs.WT[DEG_all,"log2FoldChange"],
                       shSenp6.vs.WT[DEG_all,"log2FoldChange"],
                       shSumo2.vs.WT[DEG_all,"log2FoldChange"],
                       shUba2.vs.WT[DEG_all,"log2FoldChange"],
                       shUbe2i.vs.WT[DEG_all,"log2FoldChange"],
                       shChaf1a.vs.WT[DEG_all,"log2FoldChange"],
                       shChaf1b.vs.WT[DEG_all,"log2FoldChange"],
                       KAP1_KO.vs.WT[DEG_all,"log2FoldChange"],
                       Lin28a_ko.vs.WT[DEG_all,"log2FoldChange"],
                       Lin28b_ko.vs.WT[DEG_all,"log2FoldChange"],
                       LINE1_siDux.vs.WT[DEG_all,"log2FoldChange"],
                       Max_KD.vs.WT[DEG_all,"log2FoldChange"],
                       Mga_KD.vs.WT[DEG_all,"log2FoldChange"],
                       Pcgf124_KO.vs.WT[DEG_all,"log2FoldChange"],
                       Pcgf1_KO.vs.WT[DEG_all,"log2FoldChange"],
                       Pcgf24_KO.vs.WT[DEG_all,"log2FoldChange"],
                       Pcgf356_KO.vs.WT[DEG_all,"log2FoldChange"],
                       Pcgf35_KO.vs.WT[DEG_all,"log2FoldChange"],
                       Pcgf6_KO.vs.WT[DEG_all,"log2FoldChange"],
                       Pcgf6_KD.vs.WT[DEG_all,"log2FoldChange"],
                       Rif1_KO.vs.WT[DEG_all,"log2FoldChange"],
                       RNF2_KD.vs.WT[DEG_all,"log2FoldChange"],
                       Smchd1_KO1.vs.WT[DEG_all,"log2FoldChange"],
                       Smchd1_KO2.vs.WT[DEG_all,"log2FoldChange"],
                       Smchd1_KO3.vs.WT[DEG_all,"log2FoldChange"],
                       Ythdc1_cKO_4OHT.vs.WT[DEG_all,"log2FoldChange"],
                       Eed_ko.vs.WT[DEG_all,"log2FoldChange"],
                       Ezh2_ko.vs.WT[DEG_all,"log2FoldChange"],
                       suz12_ko.vs.WT[DEG_all,"log2FoldChange"])

rownames(DEG_all_FC)<-DEG_all
colnames(DEG_all_FC)<-c("sgDnmt1","sgMyc","shSae1","shSenp6","shSumo2",
                        "shUba2","shUbe2i","shChaf1a","shChaf1b",
                        "KAP1_KO","Lin28a_ko","Lin28b_ko","LINE1_siDux",
                        "Max_KD","Mga_KD","Pcgf124_KO","Pcgf1_KO",
                        "Pcgf24_KO","Pcgf356_KO","Pcgf35_KO","Pcgf6_KO",
                        "Pcgf6_KD","Rif1_KO","RNF2_KD","Smchd1_KO1",
                        "Smchd1_KO2","Smchd1_KO3","Ythdc1_cKO_4OHT",
                        "Eed_ko","Ezh2_ko","suz12_ko")
DEG_all_FC<-na.omit(DEG_all_FC)

library(pheatmap)
pheatmap(1-cor(DEG_all_FC),
         color = colorRampPalette(c("red","yellow"))(100),
         show_rownames=T,
         show_colnames=T,
         display_numbers = T,
         filename = paste0(outDir,"/1rep.same.deg.diff.pdf")
)

pheatmap(cor(DEG_all_FC),
         color = colorRampPalette(c("yellow","red"))(100),
         show_rownames=T,
         show_colnames=T,
         display_numbers = F,
         border=FALSE,
         filename = paste0(outDir,"/1rep.same.deg.cor.pdf")
)

################################无重复样本，用edgeR做差异分析
library(edgeR)
counts<-expr[,c("Group1_shTrim28","Control_for_Group1-1")]
group<-c("shTrim28","Group1")
y <- DGEList(counts=counts,group=group)
y <- calcNormFactors(y)
y_bcv <- y
bcv <- 0.01
et <- exactTest(y_bcv,dispersion = bcv^2)
gene1 <- decideTestsDGE(et,p.value = 0.05,lfc = 0)
summary(gene1)
shTrim28_vs_WT <- et$table
shTrim28_vs_WT <- subset(shTrim28_vs_WT,PValue<0.05)
DEG_all_FC$shTrim28<-shTrim28_vs_WT[rownames(DEG_all_FC),"logFC"]


counts<-expr[,c("Group3_shAtf7ip","Control_for_Group3")]
group<-c("shAtf7ip","Group3")
y <- DGEList(counts=counts,group=group)
y <- calcNormFactors(y)
y_bcv <- y
bcv <- 0.01
et <- exactTest(y_bcv,dispersion = bcv^2)
gene1 <- decideTestsDGE(et,p.value = 0.05,lfc = 0)
summary(gene1)
shAtf7ip_vs_WT <- et$table
shAtf7ip_vs_WT <- subset(shAtf7ip_vs_WT,PValue<0.05)
DEG_all_FC$shAtf7ip<-shAtf7ip_vs_WT[rownames(DEG_all_FC),"logFC"]

counts<-expr[,c("Group4_shEset","Control_for_Group4")]
group<-c("shEset","Group4")
y <- DGEList(counts=counts,group=group)
y <- calcNormFactors(y)
y_bcv <- y
bcv <- 0.01
et <- exactTest(y_bcv,dispersion = bcv^2)
gene1 <- decideTestsDGE(et,p.value = 0.05,lfc = 0)
summary(gene1)
shEset_vs_WT <- et$table
shEset_vs_WT <- subset(shEset_vs_WT,PValue<0.05)
DEG_all_FC$shEset<-shEset_vs_WT[rownames(DEG_all_FC),"logFC"]

DEG_all_FC<-na.omit(DEG_all_FC)
write.csv(DEG_all_FC,file = "deg_fc_rep.csv")


#cor.matrix<-cor.dist(as.matrix(t(DEG_all_FC)),diag=TRUE,upper=TRUE)
#aa<-as.matrix(cor.matrix)

pheatmap(cor(DEG_all_FC),
         color = colorRampPalette(c("yellow","red"))(100),
         show_rownames=T,
         show_colnames=T,
         display_numbers = F,
         border=FALSE,
         filename = paste0(outDir,"/2all.same.deg.cor.pdf")
)

pheatmap(cor(DEG_all_FC[,1:(ncol(DEG_all_FC)-1)]),
         color = colorRampPalette(c("yellow","red"))(100),
         show_rownames=T,
         show_colnames=T,
         display_numbers = F,
         border=FALSE,
         filename = paste0(outDir,"/4neg.out.same.deg.cor.pdf")
)

#################################
# data in paper
DEG_paper<-c(rownames(shSae1.vs.WT),rownames(shSenp6.vs.WT),rownames(shSumo2.vs.WT),
             rownames(shUba2.vs.WT),rownames(shUbe2i.vs.WT),rownames(shTrim28.vs.WT),
             rownames(shChaf1a.vs.WT),rownames(shChaf1b.vs.WT)
)
DEG_paper_all<-names(table(DEG_paper))[table(DEG_paper)==8]
DEG_paper_all_FC<-data.frame(
  shSae1.vs.WT[DEG_paper_all,"log2FoldChange"],
  shSenp6.vs.WT[DEG_paper_all,"log2FoldChange"],
  shSumo2.vs.WT[DEG_paper_all,"log2FoldChange"],
  shUba2.vs.WT[DEG_paper_all,"log2FoldChange"],
  shUbe2i.vs.WT[DEG_paper_all,"log2FoldChange"],
  shTrim28.vs.WT[DEG_paper_all,"log2FoldChange"],
  shChaf1a.vs.WT[DEG_paper_all,"log2FoldChange"],
  shChaf1b.vs.WT[DEG_paper_all,"log2FoldChange"])

rownames(DEG_paper_all_FC)<-DEG_paper_all
colnames(DEG_paper_all_FC)<-c("shSae1","shSenp6","shSumo2","shUba2","shUbe2i","shTrim28","shChaf1a","shChaf1b")
DEG_paper_all_FC<-na.omit(DEG_paper_all_FC)

DEG_paper_all_FC$shAtf7ip <- shAtf7ip_vs_WT[rownames(DEG_paper_all_FC),"logFC"]
DEG_paper_all_FC$shSetdb1 <- shSetdb1_vs_WT[rownames(DEG_paper_all_FC),"logFC"]

pheatmap(cor(DEG_paper_all_FC),
         color = colorRampPalette(c("yellow","red"))(100),
         show_rownames=T,
         show_colnames=T,
         display_numbers = F,
         border=FALSE,
         filename = paste0(outDir,"/3paper.same.deg.cor.pdf")
)

#########################################################################################
log2FoldChangeCutoff <- 1
Rif1_KO_DEG <- rownames(subset(Rif1_KO.vs.WT,abs(log2FoldChange)>log2FoldChangeCutoff & padj<0.05))
RNF2_KO_DEG <- rownames(subset(RNF2_KO.vs.WT,abs(log2FoldChange)>log2FoldChangeCutoff & padj<0.05))
Pcgf6_KO_DEG <- rownames(subset(Pcgf6_KO.vs.WT,abs(log2FoldChange)>log2FoldChangeCutoff & padj<0.05))
Max_KO_DEG <- rownames(subset(Max_KO.vs.WT,abs(log2FoldChange)>log2FoldChangeCutoff & padj<0.05))
shMga_DEG <- rownames(subset(shMga.vs.WT,abs(log2FoldChange)>log2FoldChangeCutoff & padj<0.05))

library(VennDiagram)
venn.diagram(list(Rif1=Rif1_KO_DEG,
                  RNF2=RNF2_KO_DEG,
                  Pcgf6=Pcgf6_KO_DEG,
                  Max=Max_KO_DEG,
                  Mga=shMga_DEG),
             margin = 0.04,
             col = rainbow(5),
             filename=paste0(outDir,"/5factor.jpg"))
venn.diagram(list(Rif1=Rif1_KO_DEG,
                  RNF2=RNF2_KO_DEG,
                  Pcgf6=Pcgf6_KO_DEG,
                  Max=Max_KO_DEG),
             col = rainbow(4),
             filename=paste0(outDir,"/4factor.jpg"))
venn.diagram(list(Rif1=Rif1_KO_DEG,
                  RNF2=RNF2_KO_DEG,
                  Pcgf6=Pcgf6_KO_DEG),
             col = rainbow(3),
             filename=paste0(outDir,"/3factor.jpg"))

library(eulerr)
vd <- euler(c(Rif1KO=length(Rif1_KO_DEG)-length(intersect(Rif1_KO_DEG,RNF2_KO_DEG)), 
              RNF2KO=length(RNF2_KO_DEG)-length(intersect(Rif1_KO_DEG,RNF2_KO_DEG)), 
              "Rif1KO&RNF2KO"=length(intersect(Rif1_KO_DEG,RNF2_KO_DEG)))
)
pdf(width=10,paste0(outDir,"/Venn_Rif1_RNF2.pdf"))
plot(vd,
     fills = list(fill = c(colors()[616], colors()[468],"red"), alpha = 0.6),
     labels = list(col = "black", font = 4), 
     edges = FALSE,
     quantities = TRUE
     )
dev.off()

vd <- euler(c(Rif1KO=length(Rif1_KO_DEG)-length(intersect(Rif1_KO_DEG,Pcgf6_KO_DEG)), 
              Pcgf6KO=length(Pcgf6_KO_DEG)-length(intersect(Rif1_KO_DEG,Pcgf6_KO_DEG)), 
              "Rif1KO&Pcgf6KO"=length(intersect(Rif1_KO_DEG,Pcgf6_KO_DEG)))
)
pdf(width=10,paste0(outDir,"/Venn_Rif1_Pcgf6.pdf"))
plot(vd,
     fills = list(fill = c(colors()[616], colors()[468],"red"), alpha = 0.6),
     labels = list(col = "black", font = 4), 
     edges = FALSE,
     quantities = TRUE
)
dev.off()

vd <- euler(c(Rif1KO=length(Rif1_KO_DEG)-length(intersect(Rif1_KO_DEG,Max_KO_DEG)), 
              MaxKO=length(Max_KO_DEG)-length(intersect(Rif1_KO_DEG,Max_KO_DEG)), 
              "Rif1KO&MaxKO"=length(intersect(Rif1_KO_DEG,Max_KO_DEG)))
)
pdf(width=10,paste0(outDir,"/Venn_Rif1_Max.pdf"))
plot(vd,
     fills = list(fill = c(colors()[616], colors()[468],"red"), alpha = 0.6),
     labels = list(col = "black", font = 4), 
     edges = FALSE,
     quantities = TRUE
)
dev.off()

vd <- euler(c(Rif1KO=length(Rif1_KO_DEG)-length(intersect(Rif1_KO_DEG,shMga_DEG)), 
              shMga=length(shMga_DEG)-length(intersect(Rif1_KO_DEG,shMga_DEG)), 
              "Rif1KO&shMga"=length(intersect(Rif1_KO_DEG,shMga_DEG)))
)
pdf(width=10,paste0(outDir,"/Venn_Rif1_Mga.pdf"))
plot(vd,
     fills = list(fill = c(colors()[616], colors()[468],"red"), alpha = 0.6),
     labels = list(col = "black", font = 4), 
     edges = FALSE,
     quantities = TRUE
)
dev.off()

#######################
intersect(Rif1_KO_DEG,RNF2_KO_DEG)
Cls<-read.table("ucsc.mm9.repeat.txt",header = T,sep = "\t")
co<-subset(Cls,repName %in% intersect(Rif1_KO_DEG,RNF2_KO_DEG))
dim(co)
pie(table(co$renamed_repFamily))

subset(co,renamed_repFamily=="ERV1")

######## anno ######
# import annotation and rename the classification
anno <- read.csv("mm9.repeatMasker.class.csv",header=T)
DEG_all_FC$RepName <- rownames(DEG_all_FC)
DEG_all_FC.anno <- as.data.frame(merge(DEG_all_FC,anno,by="RepName"))
write.csv(DEG_all_FC.anno,file = "deg_fc_rep_anno.csv")






