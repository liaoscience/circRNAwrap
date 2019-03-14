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
. $circRNAwrap/circRNAwrap.configs

# circRNAwrap pipeline

# alignment, identification, transcript prediction, abundance estimation
# usually, we combine the alignment and identification -> detection

# full pipeline

bash $circRNAwrap/circRNAwrap_align_detections.bash sample threads dir

# transcript prediction and abundance estimation -> estimation

bash $circRNAwrap/circRNAwrap_transcripts_abundance.bash sample threads dir


# quick pipeline

. $circRNAwrap/circRNAwrap.configs
bash $circRNAwrap/circRNAwrap_quick_mode.bash sample threads dir