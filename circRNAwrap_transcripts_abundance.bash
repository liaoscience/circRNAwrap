#!/bin/sh
#$ -N abundance
#$ -S /bin/bash
#$ -q all.q
##$ -e abundance.err
##$ -o abundance.log
#$ -cwd
#$ -pe make 1


################# content, detected circRNAs, 
################# use several tools for circRNA transcript and abundance
################# 1 CIRI-AS
################# 2 RAISE
################# 3 CIRCexplorer2
################# 4 sailfish-cir


# dataset and reference

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

#### sample_1.fastq
#### sample_2.fastq


## 1. CIRI-AS version: 1.2, URL: https://sourceforge.net/projects/ciri/files/CIRI-AS/
cd $dir/${sample}

# fastq are same length, not trimmed reads
echo "CIRI-AS begin" && echo ${sample} && date
mkdir ${sample}_CIRI-AS
#time fastq-dump --split-3 $sra/${sample}.sra
time $bwa mem -t 4 -T 19 $genome ./${sample}_1.fastq ./${sample}_2.fastq > ./${sample}.bwa.sam
time perl $ciri -T 5 -I ./${sample}.bwa.sam -F $genome -A $GTF -G ./${sample}_CIRI-AS/${sample}.log -O ./${sample}_CIRI-AS/${sample}.CIRI.txt
time perl $ciri_as -S ${sample}.bwa.sam -C ./${sample}_CIRI-AS/${sample}.CIRI.txt -O ./${sample}_CIRI-AS/${sample}.CIRI.sequence -F $genome -A $GTF
rm ${sample}.bwa.sam
#rm ${sample}_1.fastq
#rm ${sample}_2.fastq

echo "CIRI-AS done" && echo ${sample} && date

## 2. RAISE
echo "RAISE begin" && echo ${sample} && date
cd $dir/${sample}
time bash $circRNAwrap/circRNAwrap_RAISE.bash ${sample} ${threads} ${dir}

echo "RAISE done" && echo ${sample} && date

## 3. CIRCexplorer2 version: 2 https://github.com/YangLab/CIRCexplorer2.git
cd $dir/${sample}
echo "circExplorer2 begin" && echo ${sample} && date
cd ${sample}_tophat
echo "circExplorer2 begin" && echo ${sample} && date
#time $CIRCexplorer2 parse -t TopHat-Fusion tophat_fusion/accepted_hits.bam > CIRCexplorer2_parse.log
#time $CIRCexplorer2 annotate -r $gene_annotation_hg19 -g $genome -b back_spliced_junction.bed -o ${sample}.CIRCexplorer2.txt > CIRCexplorer2_annotate.log
time $CIRCexplorer2 assemble -p $threads -r $gene_annotation_hg19 -m ${sample}_tophat -o assemble > CIRCexplorer2_assemble.log
time $CIRCexplorer2 denovo -r $gene_annotation_hg19 -g $genome -b back_spliced_junction.bed --abs abs --as as -m ../${sample}_tophat -d ../assemble -n $polyA_tophat -o denovo > CIRCexplorer2_denovo.log

cd $dir/${sample}
cp -r ./${sample}_tophat/as ./${sample}_CIRCexplorer2
cp -r ./${sample}_tophat/abs ./${sample}_CIRCexplorer2
cp -r ./${sample}_tophat/denovo ./${sample}_CIRCexplorer2


echo "circExplorer2 done" && echo ${sample} && date

## 4. sailfish-cir version: 1.2, URL: https://github.com/zerodel/sailfish-cir
cd $dir/${sample}
echo "sailfish_cir begin" && echo ${sample} && date
if [ ! -d ${sample}_sailfish ];then
        mkdir ${sample}_sailfish
fi
awk '{print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"_"$6"\t"$5"\t"$6}' ./circRNA_validate/${sample}.circ.txt > ./sailfish/$sample.circ.bed
time python $sailfish_cir -g $genome -a $GTF -1 ./${sample}_1.fastq -2 ./${sample}_2.fastq -o ./sailfish --bed ./sailfish/$sample.circ.bed
echo "sailfish_cir done" && echo ${sample} && date
