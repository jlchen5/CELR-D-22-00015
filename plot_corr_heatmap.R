rm(list = ls())
for.all.gene <- read.csv(file = "deg_fc_gene.csv")
prc2.gene <- read.csv(file = "gene_fc.csv")
for.all.rep <- read.csv(file = "deg_fc_rep.csv")
prc2.rep <- read.csv(file = "rep_fc.csv")

prc2.gene$geneid <- prc2.gene$X

all.gene <- merge(for.all.gene,prc2.gene,by = "geneid")
all.rep <- merge(for.all.rep,prc2.rep,by = "X")

#write.csv(all.gene,file = "all.gene.fc.csv")
#write.csv(all.rep,file = "all.rep.fc.csv")


## cal corr
rm(list = ls())
gene.corr <- read.csv(file = "all.gene.fc.csv",header = T)
rep.corr <- read.csv(file = "all.rep.fc.csv",header = T)
rownames(gene.corr) <- gene.corr[,1]
rownames(rep.corr) <- rep.corr[,1]
gene.corr <- gene.corr[,-1]
rep.corr <- rep.corr[,-1]

gc <- as.data.frame(cor(gene.corr))
rc <- as.data.frame(cor(rep.corr))

pheatmap::pheatmap(gc,cluster_cols = F,border_color = 'black',cluster_rows = F)
pheatmap::pheatmap(rc,cluster_cols = F,border_color = 'black',cluster_rows = F)

