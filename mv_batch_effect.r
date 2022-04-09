rm(list = ls())  
options(stringsAsFactors = F)
options(scipen = 5)
suppressMessages(library(sva))
suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))

################################ 1  data.frame: expr cpm repClass
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
                  "Lin28b_ko-1","Lin28b_ko-2","Lin28b_ko-3",
                  "LINE1_siControl-1","LINE1_siControl-2","LINE1_siControl-3",
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
                  "suz12_ko-1","suz12_ko-2","suz12_ko-3","WT-1","WT-2",
                  "wt_for_suz12-1","wt_for_suz12-2","wt_for_suz12-3",
                  "Ythdc1_cKO_4OHT-1","Ythdc1_cKO_4OHT-2",
                  "Ythdc1_cKO_DMSO-1","Ythdc1_cKO_DMSO-2")
##### 导入phenomena
pheno <- read.csv("pheno.csv",header = T)

##### 未去除批次效应，PCA #####
colData <- data.frame(row.names =colnames(expr),batch=pheno$batch)
dds <- DESeqDataSetFromMatrix(countData = expr, colData = colData, design = ~ batch) 
cpm <- fpm(dds,robust = FALSE)
cpm <- as.data.frame(cpm)
vst <- vst(dds)
plotPCA(vst,intgroup="batch") 
data <- plotPCA(vst,intgroup="type2",returnData = TRUE)
data$cell.name <- rownames(data)
data.pheno <- merge(data,pheno,by="type2")
data.pheno

ggplot(data.pheno, aes(x = PC1, y=PC2)) +
  geom_point(size=4.5,stroke = 0.8,aes(color=group)) + 
  theme(plot.title=element_text(hjust=0.5),title =element_text(size=12)) +
  xlab(paste0("PC1:","27% variance"))+
  ylab(paste0("PC2:","21% variance"))

##### 去除批次效应，PCA #####
expr <- ComBat(dat = as.matrix(expr),batch=pheno$batch,mod=NULL,par.prior=T)
expr[expr < 0] <- 0
expr <- round(expr,digits = 0)
#expr <- na.omit(expr)
colData <- data.frame(row.names =colnames(expr),batch=pheno$batch)
dds <- DESeqDataSetFromMatrix(countData = expr, colData = colData, design = ~ batch) 
cpm <- fpm(dds,robust = FALSE)
cpm <- as.data.frame(cpm)
vst <- vst(dds)
plotPCA(vst,intgroup="batch") 
data <- plotPCA(vst,intgroup="type2",returnData = TRUE)
data$cell.name <- rownames(data)
data.pheno <- merge(data,pheno,by="cell.name")
data.pheno

ggplot(data.pheno, aes(x = PC1, y=PC2)) +
  geom_point(size=4.5,stroke = 0.8,aes(shape=origin,color=group )) + 
  labs(title=" ") +
  theme(plot.title=element_text(hjust=0.5),title =element_text(size=12)) +
  xlab(paste0("PC1:","36% variance"))+
  ylab(paste0("PC2:","15% variance"))





library(limma)
removeBatchEffect(x, batch=NULL, batch2=NULL, covariates=NULL,
                  design=matrix(1,ncol(x),1), ...)










