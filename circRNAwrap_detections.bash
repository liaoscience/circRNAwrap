#!/bin/sh
#$ -N detection
#$ -S /bin/bash
#$ -q all.q
##$ -e detection.err
##$ -o detection.log
#$ -cwd
#$ -pe make 1


################# content, detected circRNAs, 
################# use several tools for circRNA detection
################# 1 KNIFE
################# 2 find_circ
################# 3 acfs
################# 4 CIRI-full_v2
################# 5 CIRCexplorer2
################# 6 circRNA_finder
################# 7 mapsplice
################# 8 DCC
################# 8 linear RNA abundance analysis


# dataset and reference

sample=$1
dir=/home/lilin/workdir/data/circRNA/data
threads=$2
# fastq=/home/lilin/workdir/data/circRNA/data
#### sample_1.fastq
#### sample_2.fastq


# tools and pipeline
# KNIFE bowtie and bowtie2
# find_circ bowtie2
# acfs CIRI CIRCexplorer2 circRNA_finder mapsplice
KNIFE=/home/lilin/workdir/git/KNIFE-1.4/circularRNApipeline_Standalone/
find_circ=/home/lilin/workdir/git/find_circ
acfs=/home/lilin/workdir/git/acfs/ACF_MAKE.pl
ciri=/home/lilin/workdir/git/CIRI-full_v2.0/bin/CIRI_v2.0.6/CIRI2.pl
CIRCexplorer=/home/lilin/workdir/git/CIRCexplorer/
circRNA_finder=/home/lilin/workdir/git/circRNA_finder
MAPSPLICE_DIR=/home/lilin/miniconda2/envs/metawrap-env/bin/
hisat2=/home/lilin/miniconda2/bin/hisat2
samtools=/home/lilin/workdir/git/samtools-1.4/samtools
tophat2=/home/lilin/miniconda2/bin/tophat


# reference
bowtie1=/home/lilin/workdir/reference/gatk4/hg19_ref/bowtie1/ucsc.hg19
bowtie2_ref=/home/lilin/workdir/reference/gatk4/hg19_ref/bowtie2/ucsc.hg19
fasta=/home/lilin/workdir/reference/gatk4/hg19_ref/bowtie2/ucsc.hg19.fa
acfs_config=/home/lilin/workdir/git/acfs/acfs_hg19.txt
genome=/home/lilin/workdir/reference/gatk4/hg19_ref/ucsc.hg19.fasta
GTF=/home/lilin/workdir/reference/gatk4/hg19_ref/hg19_genes.gtf
STAR=/home/lilin/workdir/reference/gatk4/hg19_ref/STAR/
mapsplice_genome=/home/lilin/workdir/reference/gatk4/hg19_ref/mapsplice/
mapsplice_index=/home/lilin/workdir/reference/gatk4/hg19_ref/mapsplice/total
hisat2_reference=/home/lilin/workdir/reference/gatk4/hg19_ref/hisat2/hg19/genome

hisat2_head_hg19="/home/lilin/workdir/data/circRNA/data/hisat2.sam"


## 1. KNIFE version: 1.4, URL: https://github.com/lindaszabo/KNIFE
#bowtie-1.1.2
#bowtie2-2.2.6
#perl
#python2.6
#R 3.1
#samtools 0.1.19
cd $dir/${sample}
if [ ! -d ${sample}_KNIFE ];then
        mkdir ${sample}_KNIFE
fi
cd $KNIFE
echo "KNIFE begin" && echo ${sample} && date
time bash completeRun.sh $dir/${sample}/ complete $dir/${sample}/ ${sample}_KNIFE 13 sam 2>&1 | tee $dir/${sample}/${sample}_KNIFE.log 
echo "KNIFE done" && echo ${sample} && date

## 2. find_circ https://github.com/marvin-jens/find_circ.git
cd $dir/${sample}
echo "find_circ begin" && echo ${sample} && date
awk '{i++; if(NR%4==1) print "@"i; else{print }}' ./${sample}_1.fastq > ${sample}_1.fq
awk '{i++; if(NR%4==1) print "@"i; else{print }}' ./${sample}_2.fastq > ${sample}_2.fq
bowtie2 -p $threads --very-sensitive --phred33 --mm -M20 --score-min=C,-15,0 -x $bowtie2_ref -q -1 ./${sample}_1.fq -2 ./${sample}_2.fq > ${sample}.bowtie2.sam 2> bowtie2.log
$samtools sort -@ $threads -o ${sample}.bowtie2.bam ${sample}.bowtie2.sam
samtools view -hf 4 ${sample}.bowtie2.bam | samtools view -Sb - > ${sample}_unmapped_bowtie2.bam
$find_circ/unmapped2anchors.py ${sample}_unmapped_bowtie2.bam > ${sample}_anchor.fastq
if [ ! -d ${sample}_fc ];then
        mkdir ${sample}_fc
fi
time bowtie2 -p $threads --reorder --mm -M 20 --score-min=C,-15,0 -q -x $bowtie2_ref -U ${sample}_anchor.fastq 2> $sample.bt2_second.log | $find_circ/find_circ.py -G $genome -p $sample -s ./$sample""_fc/$sample.sites.log > ./$sample""_fc/$sample.sites.bed 2> ./$sample""_fc/$sample.sites.reads
rm ${sample}_1.fq
rm ${sample}_2.fq
rm ${sample}_anchor.fastq
echo "find_circ done" && echo ${sample} && date

## 3. acfs https://github.com/arthuryxt/acfs.git
cd $dir/${sample}
echo "acfs begin" && echo ${sample} && date
bamToFastq -i ${sample}_unmapped_bowtie2.bam -fq ${sample}_unmapped_bowtie2.fq
if [ ! -d ${sample}_acfs ];then
        mkdir ${sample}_acfs
fi
awk '{if(NR%4==2) print;}' ${sample}_unmapped_bowtie2.fq > ./${sample}_acfs/unmapped.reads
cd ./${sample}_acfs
awk '{a[$1]++}END{for(i in a){print i,a[i]}}' unmapped.reads > unmapped.group.reads #the methods is too simple to stupid
awk '{i++; print ">"i"reads\n"$1}' unmapped.group.reads > unmapped.group.fa
awk '{i++; print ">"i"reads\t"$2}' unmapped.group.reads |sed 's/>//g' > unmapped.group.expr
sed -i '1i\newid\tsample1' unmapped.group.expr

rm unmapped.reads
cp $acfs_config SPEC_acfs_${sample}.txt
sed -i s/sample1/${sample}/g SPEC_acfs_${sample}.txt
perl $acfs SPEC_acfs_${sample}.txt SPEC_acfs_${sample}.sh
time sh SPEC_acfs_${sample}.sh
rm unmap*
rm Step*
rm circle_candidates*.p1.*
rm circle_candidates_CBR.sam
rm circle_candidates_MEA.sam
rm circle_candidates_CBR.CL.sa
rm circle_candidates_CBR.CL.amb
rm circle_candidates_CBR.CL.ann
rm circle_candidates_CBR.CL.bwt
rm circle_candidates_CBR.CL.pac
rm circle_candidates_CBR.CL
rm circle_candidates_CBR.pseudo.exon.fa
rm circle_candidates_CBR.pseudo.gene.fa
rm circle_candidates_MEA.CL.amb
rm circle_candidates_MEA.CL.ann
rm circle_candidates_MEA.CL.bwt
rm circle_candidates_MEA.CL.pac
rm circle_candidates_MEA.CL.sa
rm circle_candidates_CBR.agtf
rm circle_candidates_CBR.agtf_gene
rm circle_candidates_CBR.refFlat
rm circle_candidates_MEA.CL
rm circle_candidates_MEA.pseudo.exon.fa
rm circle_candidates_MEA.pseudo.gene.fa
rm circle_candidates_MEA.agtf
rm circle_candidates_MEA.agtf_gene
rm circle_candidates_MEA.err
rm circle_candidates_MEA.ext50
rm circle_candidates_MEA.gtf
rm circle_candidates_MEA.gtf.ext50
rm circle_candidates_MEA.refFlat
rm circle_candidates_MEA.gtf_2G
rm circle_candidates_CBR
rm circle_candidates_discard
rm circle_candidates_MEA
rm circle_candidates
rm SPEC_acfs_${sample}.sh
rm SPEC_acfs_${sample}.txt
echo "acfs end" && echo ${sample} && date

# 4, CIRI-full_v2 https://sourceforge.net/projects/ciri/files/CIRI2/
cd $dir/${sample}/
echo "ciri begin" && echo ${sample} && date
time bwa mem -t $threads -T 19 $genome ./${sample}_1.fastq ./${sample}_2.fastq > ./${sample}.bwa.sam
if [ ! -d ${sample}_CIRI ];then
        mkdir ${sample}_CIRI
fi
time perl $ciri -T $threads -I ./${sample}.bwa.sam -F $genome -A $GTF -G ./${sample}_CIRI/${sample}.log -O ./${sample}_CIRI/${sample}.CIRI.txt
echo "ciri done" && echo ${sample} && date
#rm ./${sample}.bwa.sam

## 5. circexplorer version: 1.1.5, URL: https://github.com/YangLab/CIRCexplorer
cd $dir/${sample}
echo "circexplorer begin" && echo ${sample} && date
time STAR --genomeDir $STAR --readFilesIn ./${sample}_1.fastq ./${sample}_2.fastq --runThreadN $threads --chimSegmentMin 20 --chimScoreMin 1 --alignIntronMax 100000 --outFilterMismatchNmax 4 --alignTranscriptsPerReadNmax 100000 --outFilterMultimapNmax 2 --outFileNamePrefix ./$sample --outSAMtype BAM Unsorted
cd ./
if [ ! -d ${sample}_CIRCexplorer ];then
        mkdir ${sample}_CIRCexplorer
fi
time python $CIRCexplorer/circ/star_parse.py ${sample}Chimeric.out.junction ${sample}_CIRCexplorer/${sample}_junction.txt && \
time python $CIRCexplorer/circ/CIRCexplorer.py -j ${sample}_CIRCexplorer/${sample}_junction.txt -g $genome -r $CIRCexplorer/test/data/ref.txt -o ${sample}_CIRCexplorer/${sample}
echo "circexplorer done" && echo ${sample} && date


## 6. circRNA_finder version: 1.1.5, URL: https://github.com/orzechoj/circRNA_finder.git
cd $dir/${sample}
echo "circRNA_finder begin" && echo ${sample} && date
#STAR --genomeDir $STAR --readFilesIn ../${sample}_1.fastq ./${sample}_2.fastq --runThreadN 10 --chimSegmentMin 20 --chimScoreMin 1 --alignIntronMax 100000 --outFilterMismatchNmax 4 --alignTranscriptsPerReadNmax 100000 --outFilterMultimapNmax 2 --outFileNamePrefix ./${sample} --outSAMtype BAM Unsorted
#time perl $circRNA_finder/runStar1.pl --inFile1 ./${sample}_1.fastq --inFile2 ./${sample}_2.fastq --genomeDir $STAR --outPrefix ./$sample.
if [ ! -d ${sample}_circRNA_finder ];then
        mkdir ${sample}_circRNA_finder
fi
perl $circRNA_finder/postProcessStarAlignment.pl ./ ./
rm $sample.Aligned.out.sam
mkdir ${sample}_circRNA_finder
cp ${sample}s_filteredJunctions_fw.bed ./${sample}_circRNA_finder/$sample.circRNA_finder.bed
echo "circRNA_finder done" && echo ${sample} && date

## 7. mapsplice 
source activate metawrap-env
cd $dir/${sample}
echo "mapsplice begin" && echo ${sample} && date
if [ ! -d ${sample}_mapsplice ];then
        mkdir ${sample}_mapsplice
fi
OUTPUT_DIR=./${sample}_mapsplice
READ_FILE_END1=./${sample}_unmapped_bowtie2.fq
#READ_FILE_END2=../${sample}_2.fastq
time python $MAPSPLICE_DIR/mapsplice.py \
       -1 $READ_FILE_END1 \
       -c $mapsplice_genome \
       -x $mapsplice_index \
	   --qual-scale phred33 \
       -p $threads --bam --fusion-non-canonical --min-fusion-distance 200 \
	   --gene-gtf $GTF \
       -o $OUTPUT_DIR 2>log.txt
echo "mapsplice done" && echo ${sample} && date
source deactivate metawrap-env

## DCC


## 8. linear RNA abundance
cd $dir/${sample}
echo "hisat2 begin" && echo ${sample} && date
time $hisat2 -p $threads --mm --dta -x $hisat2_reference -1 ../${sample}_1.fastq -2 ../${sample}_2.fastq -S ${sample}.hisat2.sam 2> ${sample}.map.rate
awk '{ if($2==77 || $2==141 || $2==89 || $2==133 || $2==137 || $2==69) print }' ${sample}.hisat2.sam > ${sample}.unmap.sam
awk '{m=2; if($2==77 || $2==69) m=1; if($6~/\*/) print "@"$1"_"$2"_"m"\n"$10"\n\+\n"$11}' ${sample}.unmap.sam > ${sample}.unmap.fq
$samtools sort -@ $threads -o ${sample}.hisat2.bam ${sample}.hisat2.sam
$samtools rmdup ${sample}.hisat2.bam ${sample}.hisat2.rmdup.bam
rm ${sample}.hisat2.sam &
awk '{if($6~/\*/) print }' ${sample}.unmap.sam > ${sample}.unmapped.sam
cat $hisat2_head_hg19 ${sample}.unmapped.sam > ${sample}.unmap.head.sam
$samtools view -Sb ${sample}.unmap.head.sam > ${sample}.unmap.bam
echo "hisat2 end" && echo ${sample} && date

## 9. tophat alignment
cd $dir/${sample}
echo "tophat2 begin" && echo ${sample} && date
time $tophat2 -a 6 --microexon-search -m 2 -p $threads -G $GTF -o ${sample}_tophat $bowtie2_ref ../${sample}_1.fastq ../${sample}_2.fastq
time samtools view ${sample}_tophat/unmapped.bam > ${sample}_tophat/unmapped.sam

cd ${sample}_tophat
time bamToFastq -i ./unmapped.bam -fq ./unmapped.fastq
time "/home/lilin/workdir/git/tophat-2.1.0.Linux_x86_64/tophat2" -o tophat_fusion -p $threads --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search $bowtie1 ./unmapped.fastq

echo "tophat2 end" && echo ${sample} && date



## 5. circexplorer version: 1.1.5, URL: https://github.com/YangLab/CIRCexplorer
cd $dir/${sample}

gene_annotation_ucsc_hg19="/home/lilin/workdir/reference/gatk4/hg19_ref/hg19_ens.txt"
hisat2_head_hg19="/home/lilin/workdir/data/circRNA/data/hisat2.sam"
CIRCexplorer2=CIRCexplorer2
## circExplorer2 , URL: https://github.com/YangLab/CIRCexplorer2.git

cd $dir/${sample}
echo "circExplorer2 begin" && echo ${sample} && date

CIRCexplorer2 align -p $THREADS -o circExplorer2_output -G $GTF -i $bowtie1 -j $find_circ_ref ../${sample}_1.fastq ../${sample}_2.fastq
$CIRCexplorer2 parse -t STAR ${sample}.Chimeric.out.junction > CIRCexplorer2_parse.log
$CIRCexplorer2 annotate -r $gene_annotation_ucsc_hg19 -g $genome -b back_spliced_junction.bed -o ${sample}.$CIRCexplorer2.txt > CIRCexplorer2_annotate.log
$CIRCexplorer2 assemble -p 10 --remove-rRNA --max-bundle-frags=50000 -r $gene_annotation_ucsc_hg19 -m ${sample}_tophat -o assemble > CIRCexplorer2_assemble.log
CIRCexplorer2 denovo -r $gene_annotation_ucsc_hg19 -g $genome -b back_spliced_junction.bed --abs abs --as as -m ${sample}_tophat -n pAplus_tophat -o denovo > CIRCexplorer2_denovo.log

echo "circExplorer2 done" && echo ${sample} && date

done


## 8. DCC, URL: https://github.com/dieterich-lab/DCC.git

cd $dir/${sample}
rm -r ${sample}_FUCHS
if [ ! -d ${sample}_FHCHS ];then
 mkdir ${sample}_FUCHS
fi
echo "DCC begin" && echo ${sample} && date

cd ${sample}_FUCHS

time STAR --runThreadN $threads --genomeDir $STAR --outSAMtype BAM SortedByCoordinate --readFilesIn $dir/${sample}_1.fastq $dir/${sample}_2.fastq --outFileNamePrefix ${sample}_FUCHS --quantMode GeneCounts --genomeLoad NoSharedMemory --outReadsUnmapped Fastx --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 10 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.05 --outFilterMatchNminOverLread 0.7 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15 --twopassMode Basic --alignSoftClipAtReferenceEnds No --outSAMattributes NH HI AS nM NM MD jM jI XS --sjdbGTFfile $GTF

gzip ${sample}_FUCHSUnmapped.out.mate1
mv ${sample}_FUCHSUnmapped.out.mate1.gz Unmapped_out_mate1.fastq.gz
time STAR --readFilesCommand zcat --runThreadN $threads --genomeDir $STAR --outSAMtype BAM SortedByCoordinate --readFilesIn Unmapped_out_mate1.fastq.gz --outFileNamePrefix ${sample}.mate1. --quantMode GeneCounts --genomeLoad NoSharedMemory --outReadsUnmapped Fastx --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 10 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.05 --outFilterMatchNminOverLread 0.7 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15 --twopassMode Basic --alignSoftClipAtReferenceEnds No --outSAMattributes NH HI AS nM NM MD jM jI XS --sjdbGTFfile $GTF

gzip ${sample}_FUCHSUnmapped.out.mate2
mv ${sample}_FUCHSUnmapped.out.mate2.gz Unmapped_out_mate2.fastq.gz
time STAR --readFilesCommand zcat --runThreadN $threads --genomeDir $STAR --outSAMtype BAM SortedByCoordinate --readFilesIn Unmapped_out_mate2.fastq.gz --outFileNamePrefix ${sample}.mate2. --quantMode GeneCounts --genomeLoad NoSharedMemory --outReadsUnmapped Fastx --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 10 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.05 --outFilterMatchNminOverLread 0.7 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15 --twopassMode Basic --alignSoftClipAtReferenceEnds No --outSAMattributes NH HI AS nM NM MD jM jI XS --sjdbGTFfile $GTF

echo "DCC done" && echo ${sample} && date
