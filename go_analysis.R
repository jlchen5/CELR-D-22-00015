rm(list = ls())
library(data.table)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(dplyr)


target_gene_id <- unique(read.delim("rif1_pcgf6.txt",header = F))


# transform id   *****
target_gene_id <- bitr(target_gene_id$V1, fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)

target_gene_id <- target_gene_id$ENTREZID



# 全局分析
ego <- enrichGO(OrgDb="org.Mm.eg.db",
                   gene = target_gene_id,
                   pvalueCutoff = 0.05,
                   ont = "ALL",
                   readable=TRUE)

#write.csv(ego_result_MF,file = "MCs_ego_result_MF.csv")
#write.csv(ego_result_CC,file = "MCs_ego_result_CC.csv")
#write.csv(ego_result_BP,file = "MCs_ego_result_BP.csv")


dotplot(ego, split = "ONTOLOGY", font.size = 10, #按照go中的MF/BP/CC进行换行输出
        showCategory = 5) + facet_grid(ONTOLOGY ~ ., scale = "free") 

