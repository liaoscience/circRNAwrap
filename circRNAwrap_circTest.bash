#!/bin/sh
#$ -N further_analysis
#$ -S /bin/bash
#$ -q all.q
##$ -e further_analysis.err
##$ -o further_analysis.log
#$ -cwd
#$ -pe make 1


#################  content, circRNAs further analysis 

#################  Circtest, DE circRNAs and compare with linear RNAs



. $1/circRNAwrap.configs
sample=$2
threads=$3
dir=$4
echo "-------------options----------"
echo "the config files direction is:  " $1
echo "------------------------------"
echo "the sample is:                  " $sample
echo "------------------------------"
echo "the threads numbers is:         " $threads 
echo "------------------------------"
echo "the work direction is:          " $dir
echo "------------------------------"




# circPrimer for circRNAs

# Circtest prepare, linear abundance, circRNA abundance
# RAISE result and DCC result for CircTest

Circ backsplice and linearsplice


perl $put_together_ws ./linearRNA | sort -r -k2 | sed 's/\t/,/g' > Linear.csv
perl $put_together_ws ./circRNA | sort -r -k2 | sed 's/\t/,/g' > Circ.csv

# R part

> install.packages("devtools")
> require(devtools)
> install_github('dieterich-lab/CircTest')
> library(CircTest)

> Circ <- read.delim('Circ.csv', header = T, as.is = T)
> Linear <- read.delim('Linear.csv', header = T, as.is = T)

> Circ_filtered <- Circ.filter(circ = Circ, linear = Linear, Nreplicates = 3, filter.sample = 3, filter.count = 5, percentage = 0.1, circle_description = 1)

# filter linear table by remaining circles
> Linear_filtered <- Linear[rownames(Circ_filtered),]

> test <- Circ.test(Circ_filtered, Linear_filtered, group=c(rep(1,3),rep(2,3)), circle_description = 1)

for (i in rownames(test$summary_table))  { 
    Circ.ratioplot(Circ_filtered, Linear_filtered, plotrow=i, groupindicator1=c(rep('Control',3),rep('Treatment',3)), 
		   lab_legend='Condition', circle_description = 1 )
 }
  
for (i in rownames(test$summary_table))  {
    Circ.lineplot(Circ_filtered, Linear_filtered, plotrow=i, groupindicator1=c(rep('Control',3),rep('Treatment',3)),
		  circle_description = 1 )
 }
 


 