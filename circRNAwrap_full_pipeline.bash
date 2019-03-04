# circRNAwrap pipeline

# alignment, identification, transcript prediction, abundance estimation

# usually, we combine the alignment and identification -> detection

bash circRNAwrap_align_detections.bash $1 $2

# transcript prediction and abundance estimation -> estimation

bash circRNAwrap_transcripts_abundance.bash $1 $2




