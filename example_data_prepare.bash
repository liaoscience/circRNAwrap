#!/bin/sh
#$ -N data
#$ -S /bin/bash
#$ -q all.q
##$ -e data.err
##$ -o data.log
#$ -cwd
#$ -pe make 1

cd /home/lilin/workdir/data/circRNA/

#sra data
for sample in SRR1049826 SRR1049827 SRR1049828 SRR1049829 SRR1049830 SRR1049831 SRR1049832 SRR1049833;do
/home/lilin/.aspera/connect/bin/ascp -i /home/lilin/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -T -l 200m anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR104/${sample}/${sample}.sra ./

#wget -c ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR778/$sample/$sample.sra
#fastq-dump --split-3 $sample.sra
#fastqc $sample""_1.fastq
#fastqc $sample""_2.fastq
done




#simulate reads
cd /home/lilin/workdir/data/circRNA/simulate/

#linear RNA 
~/art_illumina -ss HS25 -d simulate_linear -na -i hg19_genes.fa -o simulate_linear -ef -l 101 -f 50 -p -m 350 -s 10 -sp -rs data_prepare -qs -13 -qs2 -13

#filter linear unmapped reads
for sample in simulate_linear;do
hisat2 -p 20 --mm --dta -x /home/lilin/workdir/reference/gatk4/hg19_ref/hisat2/hg19/genome -1 ./${sample}1.fq -2 ./${sample}2.fq -S ${sample}.sam 2> ${sample}.map.rate

awk '{ if($2==83) print $1}' ${sample}.sam | awk '!a[$0]++' > ${sample}.map.id
awk '{if(NR%4!=0)ORS=" ";else ORS="\n"}1' ${sample}1.fq > ${sample}_1.txt && awk '{if(NR%4!=0)ORS=" ";else ORS="\n"}1' ${sample}2.fq > ${sample}_2.txt &
awk -F '[@ /]' 'NR==FNR{a[$1]=$0;next}NR>FNR{if($2 in a)print $0}' ${sample}.map.id ${sample}_1.txt | awk '{print $1"\n"$2"\n"$3"\n"$4}' > ${sample}_1.fastq && awk -F '[@ /]' 'NR==FNR{a[$1]=$0;next}NR>FNR{if($2 in a)print $0}' ${sample}.unmap.id ${sample}_2.txt | awk '{print $1"\n"$2"\n"$3"\n"$4}' > ${sample}_2.fastq
done


#circRNA
awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0" ":$0 }' "/home/lilin/workdir/git/KNIFE-1.4/circularRNApipeline_Standalone/index/hg19_junctions_scrambled.fa" | awk '{if($2!~"N" && $1~"rev") print }' > hg19_genes_knife.fa &
awk '{if (/^>chr1\|/) print $1"\n"$2}' hg19_genes_knife.fa > hg19_genes_knife1.fa
#cd-hit -i hg19_genes_knife_chr1.fa -o hg19_genes_knife_uniq.fa -c 0.99 -T 20 -M 0 
#circRNA
"/home/lilin/workdir/git/art_bin_MountRainier/art_illumina" -ss HS25 -ef -d simulate_circ -na -i hg19_genes_knife1.fa -f 30 -p -l 101 -f 20 -m 200 -s 10 -o simulate_circ

#filter linear mapped reads
for sample in simulate_circ;do
#hisat2 -p 20 --mm --dta -x /home/lilin/workdir/reference/gatk4/hg19_ref/hisat2/hg19/genome -1 ./${sample}1.fq -2 ./${sample}2.fq -S ${sample}.sam 2> ${sample}.map.rate

awk '{if($3~"*" || $2==77 || $2==69 || $2==133) print $1}' ${sample}.sam | awk '!a[$0]++' > ${sample}.unmap.id
awk '{if(NR%4!=0)ORS=" ";else ORS="\n"}1' ${sample}1.fq > ${sample}_1.txt && awk '{if(NR%4!=0)ORS=" ";else ORS="\n"}1' ${sample}2.fq > ${sample}_2.txt &
awk -F '[@ /]' 'NR==FNR{a[$1]=$0;next}NR>FNR{if($2 in a)print $0}' ${sample}.unmap.id ${sample}_1.txt | awk '{print $1"\n"$2"\n"$3"\n"$4}' > ${sample}_1.fastq && awk -F '[@ /]' 'NR==FNR{a[$1]=$0;next}NR>FNR{if($2 in a)print $0}' ${sample}.unmap.id ${sample}_2.txt | awk '{print $1"\n"$2"\n"$3"\n"$4}' > ${sample}_2.fastq
done

# the generated circRNA list and abundance reads
awk -F '[@|:]' '{if(NR%4==1) {$8=substr($8,1,1); print $2"\t"$6"\t"$4"\t"$3"\t"$5"\t"$8}}' circ_1.fastq | awk '{a[$0]++}END{for (i in a)print i"\t"a[i]*2 | "sort -k 1"}' > circRNA_standard.txt

awk -F '[@ /]' 'NR==FNR{a[$1]=$0;next}NR>FNR{if($2 in a)print $0}' ${sample}.unmap_both.id ${sample}_1.txt | awk '{print $1"\n"$2"\n"$3"\n"$4}' > ${sample}_both_1.fastq && awk -F '[@ /]' 'NR==FNR{a[$1]=$0;next}NR>FNR{if($2 in a)print $0}' ${sample}.unmap_both.id ${sample}_2.txt | awk '{print $1"\n"$2"\n"$3"\n"$4}' > ${sample}_both_2.fastq
