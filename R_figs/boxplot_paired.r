setwd("E:/work_tmp/OneDrive/circRNA_tools_compare/result/figs/paired_boxplot/")

library(ggpubr)

T <- read.table("list.txt", header = TRUE)

ggpaired(T, x = "tools", y = "count", id = "sample",
         color = "library", line.color = "gray", line.size = 0.8, palette = "jco", title = "  different alignment tools",  xlab = "alignment tools", ylab = "count", repel=TRUE, font.label = list(size = 14, face = "bold", color ="red")) + stat_compare_means(paired = TRUE)


# paired 
ggpaired(ToothGrowth, x = "supp", y = "len",
         color = "supp", line.color = "gray", line.size = 0.4,
         palette = "jco")+
  stat_compare_means(paired = TRUE)