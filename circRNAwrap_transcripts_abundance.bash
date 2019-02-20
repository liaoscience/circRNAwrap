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
################# 2 FUCHS
################# 3 RAISE
################# 4 CIRCexplorer2
################# 5 sailfish-cir
################# 6 CIRCpseudo


# dataset and reference

sample=$1
dir=/home/lilin/workdir/data/circRNA/data
threads=$2
# fastq=/home/lilin/workdir/data/circRNA/data
#### sample_1.fastq
#### sample_2.fastq


#tools
ciri_as="/home/lilin/workdir/git/CIRI-full_v2.0/bin/CIRI_AS_v1.2/CIRI_AS_v1.2.pl"
ciri=/home/lilin/workdir/git/CIRI-full_v2.0/bin/CIRI_v2.0.6/CIRI2.pl
sailfish_cir="/home/lilin/workdir/git/sailfish-cir/sailfish_cir.py"
pseudo="/home/lilin/workdir/git/CIRCpseudo/CIRCpseudo.pl"


# reference
genome=/home/lilin/workdir/reference/gatk4/hg19_ref/ucsc.hg19.fasta
GTF=/home/lilin/workdir/reference/gatk4/hg19_ref/hg19_genes.gtf
genome=/home/lilin/workdir/reference/gatk4/hg19_ref/ucsc.hg19.fasta  # include bwa faidx file
REF="/home/lilin/workdir/reference/gatk4/hg19_ref/annotation/hg19_gene.txt"

## 1. CIRI-AS version: 1.2, URL: https://sourceforge.net/projects/ciri/files/CIRI-AS/

cd $dir/${sample}


# fastq are same length, not trimmed reads
echo "CIRI-AS begin" && echo ${sample} && date
mkdir ${sample}_CIRI-AS
time fastq-dump --split-3 $sra/${sample}.sra
time bwa mem -t 4 -T 19 $genome ./${sample}_1.fastq ./${sample}_2.fastq > ./${sample}.bwa.sam
time perl $ciri -T 5 -I ./${sample}.bwa.sam -F $genome -A $GTF -G ./${sample}_CIRI-AS/${sample}.log -O ./${sample}_CIRI-AS/${sample}.CIRI.txt
time perl $ciri_as -S ${sample}.bwa.sam -C ./${sample}_CIRI-AS/${sample}.CIRI.txt -O ./${sample}_CIRI-AS/${sample}.CIRI.sequence -F $genome -A $GTF
rm ${sample}.bwa.sam
rm ${sample}_1.fastq
rm ${sample}_2.fastq

echo "CIRI-AS done" && echo ${sample} && date

## 5. sailfish-cir version: 1.2, URL: https://github.com/zerodel/sailfish-cir
cd $dir/${sample}
echo "sailfish_cir begin" && echo ${sample} && date
mkdir sailfish
awk '{print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"_"$6"\t"$5"\t"$6}' ./RAISE/kept_7_tools.bed > ./sailfish/$sample.circ.bed
time python $sailfish_cir -g $genome -a $GTF -1 ./${sample}_1.fastq -2 ./${sample}_2.fastq -o ./sailfish --bed ./sailfish/$sample.circ.bed
echo "sailfish_cir done" && echo ${sample} && date


## 6. CIRCpseudo version: 1.2, URL: https://github.com/YangLab/CIRCpseudo
cd $dir/${sample}
echo "pseudo begin" && echo ${sample} && date
time perl $pseudo -circ ./RAISE/${sample}.circ.bed -ref $REF -genome $genome -bwaidx $genome -output ${sample}.pseudo.txt 1> pseudo.log
echo "pseudo done" && echo ${sample} && date
