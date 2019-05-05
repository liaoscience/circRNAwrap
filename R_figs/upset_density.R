#
#install.packages("UpSetR")
library(UpSetR)
library(UpSetR)
library(reshape)
library(ggplot2)
library(plyr)
library(reshape)
library(plyr)
setwd("E:/work_tmp/OneDrive/circRNA_tools_compare/result/raw_result/")
dir()
sample <- "SRR1049828"
sum <- as.numeric("281")

#data <- read.table("./circ2_test/circ2.upset1.txt",header=TRUE)
data <- read.table("./SRR1049828/SRR1049828.upset1.txt",header=TRUE)
# 4 10
# 3.88 10.18
# 26 202 33 296 32 298 31 280 30 263 29 238 28 230 27 244 circ2 226 mix 281 circ 282
#pdf (file=paste("circ2_upset", ".pdf", sep=""),  width=5.65, height=2.75)
#upset(data, main.bar.color = rainbow(sum), 
      sets = c("acfs", "CIRCexplorer", "CIRCexplorer2", "circRNA_finder", "CIRI", "DCC", "find_circ", "KNIFE", "mapsplice", "standard"),#
      mb.ratio = c(0.55, 0.45),#
      order.by = "freq", #
      keep.order = TRUE, #keep.order
      show.numbers = FALSE,
      number.angles = 30, #
      point.size = 1.5, line.size = 0.3, #
      mainbar.y.label = "Genre Intersections", sets.x.label = "methods", #
      
      text.scale = c(1, 1, 1, 1, 1, 1)) #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
#dev.off()
#library(grid)
#grid.text(paste(sample,"Choosing the Top Largest Sets and Plot Formatting", sep = " ", collapse = NULL),x = 0.65, y=0.95, gp=gpar(fontsize=12))

another.plot <- function(data, x, y) {
    myplot <- (ggplot(data, aes_string(x = x)) + geom_density(aes(fill = factor(data$fre)), 
        alpha = 0.4) + theme_bw() + theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black")) + xlim(0,60) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), legend.key.size = unit(0.4, "cm")))
}
upset(data,
      sets = c("acfs", "CIRCexplorer", "CIRCexplorer2", "circRNA_finder", "CIRI", "DCC", "find_circ", "KNIFE", "mapsplice", "standard"),#
      mb.ratio = c(0.55, 0.45),#
      order.by = "freq", #
      keep.order = TRUE, #keep.order
      number.angles = 30, #
      point.size = 2, line.size = 1, #
      mainbar.y.label = "Genre Intersections", sets.x.label = "methods", #
      text.scale = c(1.3, 1.3, 1, 1, 1.5, 1), 
      attribute.plots = list(gridrows = 50, plots = list(list(plot = histogram, 
        x = "fre", color = "red", queries = F, xlim=c(1,10)), list(plot = scatter_plot, x = "fre", 
        y = "average", queries = T, xlim=c(1,10)), list(plot = another.plot, x = "average", 
        y = "fre", queries = F)), ncols = 3))
        
library(grid)
grid.text(paste(sample,"Choosing the Top Largest Sets and Plot Formatting", sep = " ", collapse = NULL),x = 0.65, y=0.95, gp=gpar(fontsize=12))