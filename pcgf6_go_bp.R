rm(list = ls())
library(ggplot2)
library(dplyr)
library(stringi)
library(clusterProfiler)
library(stringr)
library(enrichplot)


GO_BP=read.csv("pcgf6_go_bp.csv",header = T)


##
ggplot(data=GO_BP,aes(x=Term,y=Count))+
  geom_point(aes(size=Count,color=-log10(PValue)))+
  scale_colour_gradient(low="blue",high="red")+
  theme_bw()+
  labs()+
  theme(axis.text.x = element_text(face = "bold",color = "gray50",angle = 90,vjust = 1,hjust = 1))


