#火山图绘制
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
rm(list=ls())

dataset <- read.csv("rif1_ko_rep_deg.csv",header = T)
dataset$RepName <- dataset$X

#anno.raw <- read.csv("mm9.repeatMasker.class.csv",header = T)
#dataset_anno <- as.data.frame(merge(dataset,anno.raw,by="RepName"))
#write.csv(dataset_anno,file = "Pcgf6_kd_deg_anno.csv")

# 设置pvalue和logFC的阈值
cut_off_pvalue = 0.05
cut_off_logFC = 1
# 根据阈值分别为上调基因设置‘up’，下调基因设置‘Down’，无差异设置‘Stable’，保存到change列
# 这里的change列用来设置火山图点的颜色
dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_logFC, 
                        ifelse(dataset$log2FoldChange> cut_off_logFC ,'Up','Down'),
                        'Stable')
# 绘制火山图====================================
p <- ggplot(
  #设置数据
  dataset, 
  aes(x = log2FoldChange, 
      y = -log10(pvalue), 
      colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  
  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  
  # 坐标轴
  labs(x="log2(Foldchange)",
       y="-log10(P-value)")+
  theme_bw()+
  
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )

# 这里设置logFC值大于5的差异基因来标记
dataset$label = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= 1, as.character(dataset$X),"")
p + geom_text_repel(data = dataset, aes(x = dataset$log2FoldChange, 
                                        y = -log10(dataset$pvalue), 
                                        label = label),
                    size = 3,box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.8, "lines"), 
                    segment.color = "black", 
                    show.legend = FALSE)
