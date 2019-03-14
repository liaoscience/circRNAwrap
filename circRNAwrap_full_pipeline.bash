# circRNAwrap pipeline

# alignment, identification, transcript prediction, abundance estimation

# usually, we combine the alignment and identification -> detection

bash circRNAwrap_align_detections.bash $1 $2

# transcript prediction and abundance estimation -> estimation

bash circRNAwrap_transcripts_abundance.bash $1 $2


# parallel run
parallel -j 2 --xapply --link 'bash circRNAwrap_transcripts_abundance.bash {1} 8' ::: sample1 sample2 sample3

parallel -j 2 --xapply --link 'bash circRNAwrap_align_detections.bash {1} 8' ::: sample1 sample2 sample3

# list of samples
parallel -j 2 --xapply --link 'bash circRNAwrap_transcripts_abundance.bash {1} 8' :::: sample.txt

sample.txt: 
sample1
sample2
sample3
sample4