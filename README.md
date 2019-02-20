######################circRNAwrap is a flexible pipeline for circRNA expand analysis

#####################


circRNAwrap
==============

Contains the following files:
- prepare_files.txt
- circRNAwrap_detections.bash
- circRNAwrap_abundance_transcripts.bash

test on 8 circRNA detection tools
 - KNIFE
 - find_circ
 - CIRI2
 - CIRCexplorer2
 - mapsplice
 - acfs
 - circRNA_finder
 - DCC

circRNA sequence prediction
 - RAISE
 - CIRI-as
 - CIRCexplorer2
 
circRNA abundance estimation
 - sailfish-cir

expanded apply
 - parallel
 


These scripts have been tested on various Linux distributions. Before they can be run, make sure that the following prerequisites are installed:
 - perl
 - awk
 - STAR (versions 2.3.1)
 - samtools (1.3)
 - bwa (0.7.3)
 - hisat2 (2.0.4)
 - bowtie2 (2.2.5)
 - bedtools (2.24.0)
 - stringtie (1.3.0)
 - htseq-count (0.6.0)
 - github


To run the scripts to identify circular RNAs, first run hisat2 and STAR, once for each data set:

bash circRNAwrap_detections.bash

Next, run the post processing scripts:

bash circRNAwrap_abundance_transcripts.bash

For each library the following output files or directions are produced:

folds:
a) <lib_name>_fc: find_circ output
b) <lib_name>_acfs: acfs output
c) <lib_name>_mapsplice: mapsplice output
d) circRNA_validation: circRNA realignment details
files:
a) <lib name>_filteredJunctions.bed: circRNA_finder output
b) <lib name>_total.txt: circRNA abundance
b) <lib name>_total.txt: circRNA abundance detail

email: llsnnu@163.com
