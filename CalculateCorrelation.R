rm(list = ls())
library(pheatmap)
#library(amap)
#library(ggplot2)
#library(gplots)
#library(BiocParallel)
data <- read.csv("deg_fc_gene2.csv",header = T)
rownames(data) <- data[,1]
data <- data[,-1]

# calculate correlation
data_corr <- cor(data)
data_corr <- as.data.frame(data_corr)


pheatmap(data_corr,cluster_cols  = T,fontsize_row=12,
         color = colorRampPalette(c("#4575b4","#abd9e9","#abd9e9","#ffffbf","#fee090","#fee090","#d73027"))(100),
         fontsize_col = 12,angle_col = 315)

############################################

rm(list = ls())
library(pheatmap)
data <- read.csv("res2.csv",header = T)
rownames(data) <- data[,1]
data <- data[,-1]

pheatmap(data,cluster_cols  = T,fontsize_row=12,display_numbers = T,
         color = colorRampPalette(c("#fee090","#fee090","#f46d43","#f46d43","#d73027"))(100),
         fontsize_col = 12,angle_col = 315)
