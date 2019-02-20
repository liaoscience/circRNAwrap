## 1. CIRI-AS version: 1.2, URL: https://sourceforge.net/projects/ciri/files/CIRI-AS/

# fastq are same length, not trimmed reads
echo "CIRI-AS begin" && echo ${sample} && date
mkdir ${sample}_CIRI-AS
#time fastq-dump --split-3 $sra/${sample}.sra
time bwa mem -t $threads -T 19 $genome ./${sample}_1.fastq ./${sample}_2.fastq > ./${sample}.bwa.sam
time perl $ciri -T $threads -I ./${sample}.bwa.sam -F $genome -A $GTF -G ./${sample}_CIRI-AS/${sample}.log -O ./${sample}_CIRI-AS/${sample}.CIRI.txt
time perl $ciri_as -S ${sample}.bwa.sam -C ./${sample}_CIRI-AS/${sample}.CIRI.txt -O ./${sample}_CIRI-AS/${sample}.CIRI.sequence -F $genome -A $GTF
rm ${sample}.bwa.sam
rm ${sample}_1.fastq
rm ${sample}_2.fastq
echo "CIRI-AS done" && echo ${sample} && date

## 2. RAISE 
echo "RAISE begin" && echo ${sample} && date

time circRNAwrap_RAISE_total.bash ${sample}

echo "RAISE done" && echo ${sample} && date

## 3. CIRCexplorer2 version: 2 https://github.com/YangLab/CIRCexplorer2.git

echo "circExplorer2 begin" && echo ${sample} && date
time $CIRCexplorer2 parse -t TopHat-Fusion tophat_fusion/accepted_hits.bam > CIRCexplorer2_parse.log
time $CIRCexplorer2 annotate -r $gene_annotation_hg19 -g $genome -b back_spliced_junction.bed -o ${sample}.CIRCexplorer2.txt > CIRCexplorer2_annotate.log
time $CIRCexplorer2 assemble -p $threads -r $gene_annotation_hg19 -m ${sample}_tophat -o assemble > CIRCexplorer2_assemble.log

time $CIRCexplorer2 denovo -r $gene_annotation_hg19 -g $genome -b back_spliced_junction.bed --abs abs --as as -m ../${sample}_tophat -d ../assemble -n /home/lilin/workdir/data/circRNA/data/polyA/SRR8083958_tophat/ -o denovo > CIRCexplorer2_denovo.log
echo "circExplorer2 done" && echo ${sample} && date



