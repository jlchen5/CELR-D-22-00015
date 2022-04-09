suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))
rm(list = ls())  
options(stringsAsFactors = F)

#load counts
filelist<-dir("coding_genes.count/",pattern = "*simple")
filelist <- paste("coding_genes.count/",filelist,sep="")
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


colnames(count) <- gsub(pattern="coding_genes.count/",replacement = "",colnames(count))
colnames(count) <- gsub(pattern=".coding_genes.count.simple",replacement = "",colnames(count))
head(count)

rownames(count) <- count$Geneid
count <- count[,-1]

head(count)


# get pheno.
library("xlsx")
pheno.raw <- read.xlsx("pheno.xlsx",sheetIndex = 1)

pheno <- pheno.raw[pheno.raw$ID %in% colnames(count),] 
rownames(pheno) <- pheno$ID
count <- count[,pheno$ID]
pheno  <- pheno[colnames(count),]
pheno

table(colnames(count)==pheno$ID)
colnames(count) <- pheno$ID

#PCA
colData <- data.frame(row.names =colnames(count),type=pheno$type)

dds <- DESeqDataSetFromMatrix(countData = count, colData = colData, design = ~ type) 
cpm <- fpm(dds,robust = FALSE)
cpm <- as.data.frame(cpm)
cpm$gene_name <- rownames(cpm)

vst <- vst(dds)
plotPCA(vst,intgroup="type") 
data <- plotPCA(vst,intgroup="type",returnData = TRUE)
data$ID <- rownames(data)
data.pheno <- merge(data,pheno,by="ID")
data.pheno

ggplot(data.pheno, aes(x = PC1, y=PC2)) +
  geom_point(size=4.5,stroke = 0.8,aes(color=group )) + 
  labs(title=" ") + theme_bw()+
  theme(plot.title=element_text(hjust=0.5),title =element_text(size=12) )  +
  xlab(paste0( " PC1:  " ," 73% variance " ))+
  ylab(paste0( " PC2:  " ," 11% variance " ))



# import annotation and rename the classification
anno <- read.csv("rawData/Transposable_Elements_in_Human_Genome_20190819.csv",header=T)
anno <- anno[,c("RepName","Family","Superfamily","Class","repeatMasker.repName")]
dim(anno) ## In total, there are 15819 repeats 
# only repeats with the repeatMasker.repName can be dealt in the PCA and DEG analysis, so those repeats without repeatMasker.repName will be removed.
anno <- anno[complete.cases(anno$repeatMasker.repName),]
dim(anno) ## After filtering, 15476 repeats kept.
anno$Class <- gsub(pattern = "DNA transposon",replacement = "DNA Transposon",anno$Class)
anno[anno$RepName %in% c("(CATTC)n","(GAATG)n"),]
anno <- anno[c(-4748,-7259),] ## After removeing duplicates, 15474 repeats kept.
rownames(anno) <- anno$repeatMasker.repName
dim(anno)


dds <- DESeq(dds)
res <- results(dds,contrast = c("type","shRNF2","WT"))
res <- as.data.frame(res)
res$gene_name <- rownames(res)

cpm.res <- as.data.frame(merge(cpm,res,by="gene_name"))
cpm.res$gene_name <- rownames(res)
#cpm.res.anno <- merge(anno,cpm.res,by="repeatMasker.repName")
write.csv(cpm.res,file="genes_shRNF2.csv")









