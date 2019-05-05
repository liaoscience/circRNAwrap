setwd("E:/work_tmp/OneDrive/circRNA_tools_compare/result/time_memory/")
df1 <- read.table("time_memory_align.txt", head=T)
df1 <- data.frame(df1)
library(ggplot2)
library(RColorBrewer)
align <- ggplot(df1, aes(x = software)) +
  geom_col(aes( y = runtime, fill="darkturquoise")) +
  geom_text(aes(y = runtime, label = runtime ), fontface = "bold", vjust = 1.4, color = "black", size = 3) +
  geom_line(aes(y = memory * 10, group = 1, color = 'blackline'), size=1, linejoin="round") +
  geom_text(aes(y = memory * 10, label = round(memory, 1)), vjust = 1.4, color = "black", size = 3) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~ . / 1)) +
  scale_fill_manual('', labels = 'runtime', values = "deeppink4") +
  scale_color_manual('', labels = 'mem', values = 'black') +
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())

df2 <- read.table("time_memory_detection.txt", head=T)
df2 <- data.frame(df2)
library(ggplot2)
detection <- ggplot(df2, aes(x = software)) +
  geom_col(aes( y = runtime, fill="deepskyblue1")) +
  geom_text(aes(y = runtime, label = runtime ), fontface = "bold", vjust = 1.4, color = "black", size = 3) +
  geom_line(aes(y = memory * 10, group = 1, color = 'blackline'), size=1, linejoin="round") +
  geom_text(aes(y = memory * 10, label = round(memory, 1)), vjust = 1.4, color = "black", size = 3) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~ . / 1)) +
  scale_fill_manual('', labels = 'runtime', values = "aquamarine4") +
  scale_color_manual('', labels = 'mem', values = 'black') +
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())

df3 <- read.table("time_memory_abundance.txt", head=T)
df3 <- data.frame(df3)
library(ggplot2)
abundance <- ggplot(df3, aes(x = software)) +
  geom_col(aes( y = runtime, fill="redfill")) +
  geom_text(aes(y = runtime, label = runtime ), fontface = "bold", vjust = 1.4, color = "black", size = 3) +
  geom_line(aes(y = memory * 10, group = 1, color = 'blackline'), size=1, linejoin="round") +
  geom_text(aes(y = memory * 10, label = round(memory, 1)), vjust = 1.4, color = "black", size = 3) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~ . / 1)) +
  scale_fill_manual('', labels = 'runtime', values = "deepskyblue3") +
  scale_color_manual('', labels = 'mem', values = 'black') +
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())

library(RColorBrewer)
library("gridExtra")
p <- grid.arrange(align, detection, abundance, ncol=3, nrow=1, widths=c(4,4,4))
require(grid)

pdf( "runtime.pdf",height=3,width=12)
plot(p)
dev.off()

