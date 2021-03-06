# circRNAwrap pipeline, for quick pipeline, we applied the CIRI and CIRI-AS for the circRNA analysis, owning to that CIRI are sensitive to detect circRNA candidates, and CIRI-AS are useful for circRNA transcript prediction. 

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
# alignment, identification, transcript prediction, abundance estimation

cd $dir/${sample}
# fastq are same length, not trimmed reads
echo "CIRI-AS begin" && echo ${sample} && date
mkdir ${sample}_CIRI-AS
#time fastq-dump --split-3 $sra/${sample}.sra
time $bwa mem -t $threads -T 19 $genome ./${sample}_1.fastq ./${sample}_2.fastq > ./${sample}.bwa.sam
time perl $ciri -T $threads -I ./${sample}.bwa.sam -F $genome -A $GTF -G ./${sample}_CIRI-AS/${sample}.log -O ./${sample}_CIRI-AS/${sample}.CIRI.txt
time perl $ciri_as -S ${sample}.bwa.sam -C ./${sample}_CIRI-AS/${sample}.CIRI.txt -O ./${sample}_CIRI-AS/${sample}.CIRI.sequence -F $genome -A $GTF
rm ${sample}.bwa.sam
#rm ${sample}_1.fastq
#rm ${sample}_2.fastq
echo "CIRI-AS done" && echo ${sample} && date
