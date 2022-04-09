rm(list = ls())

data <- read.csv("deg_fc_gene2.csv",header = T)
rownames(data) <- data[,1]
data <- data[,-1]

data <- t(data)
data.pca <- prcomp(data)
summary(data.pca)
screeplot(data.pca, npcs = 10, type = "lines")

pca_prcp_contrib <- data.pca$sdev %>% .^2 %>% {./sum(.) * 100} %>% .[1:2] %>% signif(digits = 4)
pca_prcp_contrib

plot(data.pca$x,cex = 3,main = "PCA analysis", col = c(rep("pink1",1),rep("pink2",1)),pch = c(rep(16,3),rep(16,3)))
# 添加分隔线
abline(h=0,v=0,lty=2,col="gray")
# 添加标签
text(data.pca$x,labels = rownames(data.pca$x),pos = 4,offset = 0.6,cex = 1)
# 添加图例
legend("bottomright",title = "Sample",inset = 0.01,legend = rownames(data.pca$x),col = c(rep("red",3),rep("blue",3)),pch = c(rep(16,3),rep(17,3)))

#############################
library(tidyverse)
library(magrittr)
library(glue)
library(data.table)
library(irlba)
library(ggplot2)
library(ggrepel)
library(ggsci)
library(scales)

PCAdata <- prcomp_irlba(data, n = 3, scale. = T)
PCprop <- (PCAdata$sdev)^2 %>% {round(. / sum(.), 3) * 100}

plotData <- as.data.table(PCAdata$rotation)
plotData[, id := colnames(data)][,][]

ggplot(plotData, aes(x = PC1, y = PC2)) +
  geom_point(show.legend = F,size=8) +
  geom_text_repel(aes(label = id)) +
  scale_color_d3() +
  labs(
    x = str_c("PC1 : ", PCprop[1], " % Variance"),
    y = str_c("PC2 : ", PCprop[2], " % Variance")) +
  theme_bw()



