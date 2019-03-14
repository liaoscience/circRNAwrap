######################circRNAwrap is a flexible pipeline for circRNA expand analysis

##################### from the RNA-seq data paired fastq files to circRNA list, transcripts predicted info, and circRNA abundance. 


circRNAwrap
==============

Contains the following files:
- circRNAwrap_configs.txt
- circRNAwrap_align_detections.bash
- circRNAwrap_transcript_abundance.bash

examples of this pipeline
- circRNAwrap_quick_pipeline.bash
- circRNAwrap_full_pipeline.bash


include 8 circRNA detection tools, users could add more tools if necessary
 - KNIFE
 - find_circ
 - CIRI2
 - CIRCexplorer2 and CIRCexplorer
 - mapsplice
 - acfs
 - circRNA_finder
 - DCC

include 3 circRNA sequence prediction tools
 - RAISE
 - CIRI-as
 - CIRCexplorer2
 
include 1 circRNA abundance estimation tool
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


before run the scripts, firstly, prepare the index and annotation files. owning to that we applied several different tools for circRNA identification, so we have to inistall the tools, and prepare the index for each softwares.

the detail in index and annotation, which include tools install, index and annotation prepare  
- index_and_annotation.txt
 
 

To run the scripts to identify circular RNAs, just sub the shell script:

bash circRNAwrap_align_detections.bash

Next, run the post processing scripts:

bash circRNAwrap_transcript_abundance.bash


For each library the following output files or directions are produced in full pipeline:

folds:
a) <lib_name>_KNIFE:                KNIFE output
b) <lib_name>_find_circ:            find_circ output
c) <lib_name>_acfs:                 acfs output
d) <lib_name>_CIRI:                 CIRI output, CIRI-AS output
e) <lib_name>_CIRCexplorer2:        CIRCexplorer2 output, include backsplice and alternative splice
f) <lib_name>_circRNA_finder:       circRNA_finder output
g) <lib_name>_mapsplice:            mapsplice output
h) <lib_name>_DCC:                  DCC output

i) <lib_name>_circRNA_validation:   circRNA realignment details
j) <lib_name>_RAISE:                RAISE output
k) <lib_name>_sailfish-cir:         sailfish-cir output

l) <lib_name>_sum:                  sum result
