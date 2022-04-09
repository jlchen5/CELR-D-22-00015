################# 0.init environment #################
rm(list = ls())
library(pheatmap)
################# 1.read table #################
rawdata <- read.csv('density.csv',header = T)
rawdata <- data.frame(rawdata)
rownames(rawdata) <- rawdata[,2]

anno <- data.frame(rawdata[,1])
rownames(anno) <- rownames(rawdata)
colnames(anno) <- 'ERVs'

plotdata <- data.frame(rawdata[,3:6])


################# 2.plot heatmap #################
pheatmap(plotdata,scale = 'column',cluster_cols = F,border_color = 'black',
         #color = colorRampPalette(c("blue1", "white", "red1"))(100),
         cluster_rows = F,annotation_row = anno)
