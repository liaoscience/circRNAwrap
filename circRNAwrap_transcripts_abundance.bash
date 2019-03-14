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

sample=$1
dir=/home/lilin/workdir/data/circRNA/data
threads=$2
# fastq=/home/lilin/workdir/data/circRNA/data
#### sample_1.fastq
#### sample_2.fastq


# tools
ciri_as="/home/lilin/workdir/git/CIRI-full_v2.0/bin/CIRI_AS_v1.2/CIRI_AS_v1.2.pl"
ciri=/home/lilin/workdir/git/CIRI-full_v2.0/bin/CIRI_v2.0.6/CIRI2.pl
sailfish_cir="/home/lilin/workdir/git/sailfish-cir/sailfish_cir.py"


# reference
genome=/home/lilin/workdir/reference/gatk4/hg19_ref/ucsc.hg19.fasta  # include bwa faidx file
GTF=/home/lilin/workdir/reference/gatk4/hg19_ref/hg19_genes.gtf  
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

## 2. RAISE
echo "RAISE begin" && echo ${sample} && date

time bash circRNAwrap_RAISE.bash ${sample}

echo "RAISE done" && echo ${sample} && date

## 3. CIRCexplorer2 version: 2 https://github.com/YangLab/CIRCexplorer2.git

echo "circExplorer2 begin" && echo ${sample} && date
#time $CIRCexplorer2 parse -t TopHat-Fusion tophat_fusion/accepted_hits.bam > CIRCexplorer2_parse.log
#time $CIRCexplorer2 annotate -r $gene_annotation_hg19 -g $genome -b back_spliced_junction.bed -o ${sample}.CIRCexplorer2.txt > CIRCexplorer2_annotate.log
time $CIRCexplorer2 assemble -p $threads -r $gene_annotation_hg19 -m ${sample}_tophat -o assemble > CIRCexplorer2_assemble.log
time $CIRCexplorer2 denovo -r $gene_annotation_hg19 -g $genome -b back_spliced_junction.bed --abs abs --as as -m ../${sample}_tophat -d ../assemble -n $polyA_tophat -o denovo > CIRCexplorer2_denovo.log
echo "circExplorer2 done" && echo ${sample} && date

## 4. sailfish-cir version: 1.2, URL: https://github.com/zerodel/sailfish-cir
cd $dir/${sample}
echo "sailfish_cir begin" && echo ${sample} && date
mkdir sailfish
awk '{print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"_"$6"\t"$5"\t"$6}' ./RAISE/${sample}.circ.txt > ./sailfish/$sample.circ.bed
time python $sailfish_cir -g $genome -a $GTF -1 ./${sample}_1.fastq -2 ./${sample}_2.fastq -o ./sailfish --bed ./sailfish/$sample.circ.bed
echo "sailfish_cir done" && echo ${sample} && date