#!/bin/sh
#$ -N further_analysis
#$ -S /bin/bash
#$ -q all.q
##$ -e further_analysis.err
##$ -o further_analysis.log
#$ -cwd
#$ -pe make 1


#################  content, circRNAs further analysis 
#################  circRNA view, miRNA target prediction, DE circRNAs, primer design
#################  Ularcirc, view, miRNA target, ORF
#################  Circtest, DE circRNAs and compare with linear RNAs
#################  circPrimer, primer design
#################  


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
awk '{if(/^chr/) print $1}' ./${sample}/sailfish/quant_circular/quant.sf > ${sample}.circPrimer.list

awk '{if(/^chr/ || NR==1) print $1"\t"$5}' ./${sample}/sailfish/quant_circular/quant.sf > ${sample}.circRNA.profile
awk '{if(!/^chr/) print $1"\t"$5}' ./${sample}/sailfish/quant_circular/quant.sf > ${sample}.linearRNA.profile

# Circtest prepare, linear abundance, circRNA abundance

# Ularcirc data prepare STAR output tab junctions, prepare in align and detection step.





