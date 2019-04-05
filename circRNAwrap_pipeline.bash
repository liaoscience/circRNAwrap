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



dir=/home/lilin/workdir/data/circRNA/data/
threads=2
sample=head1000
circRNAwrap=/home/lilin/workdir/git/circRNAwrap_v3/

bash $circRNAwrap/circRNAwrap_align_detections.bash $circRNAwrap $sample $threads $dir

# transcript prediction and abundance estimation -> estimation

bash $circRNAwrap/circRNAwrap_transcripts_abundance.bash $circRNAwrap $sample $threads $dir


# quick pipeline

bash $circRNAwrap/circRNAwrap_quick_mode.bash $circRNAwrap $sample $threads $dir