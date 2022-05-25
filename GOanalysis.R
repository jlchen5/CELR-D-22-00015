rm(list = ls())
library(data.table)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(dplyr)


target_gene_id <- unique(read.delim("rif1_pcgf6.txt",header = F))
target_gene_id <- bitr(target_gene_id$V1, fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)
target_gene_id <- target_gene_id$ENTREZID


ego <- enrichGO(OrgDb="org.Mm.eg.db",
                   gene = target_gene_id,
                   pvalueCutoff = 0.05,
                   ont = "ALL",
                   readable=TRUE)


dotplot(ego, split = "ONTOLOGY", font.size = 10, showCategory = 5) + facet_grid(ONTOLOGY ~ ., scale = "free") 

