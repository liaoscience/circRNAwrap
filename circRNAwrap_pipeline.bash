#!/bin/sh
#$ -N detection
#$ -S /bin/bash
#$ -q all.q
##$ -e detection.err
##$ -o detection.log
#$ -cwd
#$ -pe make 1


# pipeline direction, configs file direction
# software direction
# before runing, give the configs files

# circRNAwrap pipeline

# alignment, identification, transcript prediction, abundance estimation
# usually, we combine the alignment and identification -> detection

# full pipeline
#### ./sample/
#### sample_1.fastq
#### sample_2.fastq



dir=/home/lilin/workdir/data/sra_circRNA/
threads=10
sample=SRR1049826
circRNAwrap=/home/lilin/workdir/git/circRNAwrap_v3/

nohup bash $circRNAwrap/circRNAwrap_align_detections.bash $circRNAwrap $sample $threads $dir 1> test &

# transcript prediction and abundance estimation -> estimation

bash $circRNAwrap/circRNAwrap_transcripts_abundance.bash $circRNAwrap $sample $threads $dir


# quick pipeline

bash $circRNAwrap/circRNAwrap_quick_mode.bash $circRNAwrap $sample $threads $dir


nohup parallel -j 2 "bash CIRCexplorer2_assembl2.bash {1}" ::: SRR1049832 SRR1049833 &


expand 

nohup bash /home/lilin/workdir/git/circRNAwrap_v3/circRNAwrap_transcripts_abundance.bash /home/lilin/workdir/git/circRNAwrap_v3 SRR1049826 8 /home/lilin/workdir/data/sra_circRNA/ &



change CIRCexplorer2_assemble from cufflinks to stringtie

gtfToGenePred SRR1049830.gtf SRR1049830.genePred

awk '{if($3~/transcript/) print $10,$14,$16,$18}' SRR1049830.gtf | sort -rk2 | awk '!a[$1]++' | sed 's/"//g' | sed 's/;//g' | sort -k1 > SRR1049830.gname

awk -F '[\.\t ]' 'NR==FNR{a[$1"\t"$2]=$3"\t"$4;next}NR>FNR{if($1"\t"$2 in a)print a[$1"\t"$2]"\t"$0}' SRR1049830.gname SRR1049830.genePred | awk '{if(/^[0-9]/) $2=$3; else $3=$1; $1=""; print }' | sed 's/  / /g' | sed 's/^ //g' | sed 's/ /\t/g' > SRR1049830_ref.txt

awk '{if(/^[0-9]/) $2=$3; else $3=$1; $1=""; print }'


