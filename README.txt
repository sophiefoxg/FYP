#### FYP predicting enhancer promoter interactions ####

# Directory structure #

## R script EPI.R ##
uses .txt ESC promoter-other interactions
creates a genomic interactions object
creates metadata file that shows distance between interaction



## slurm scripts for RNA-seq processing ##
# RNA-esc
# using provided RNA seq data QC, trim and multi QC revealed : GSM723776 dataset to be better
#justification: ERR031629_trim is problematic:36% of reads failed QC
# star alignment to mm9 genome
# featurecounts created summary file of genes expressed
# TPM.R calculates transcripts per million
# (still need to do) execute slurm script. this runs all RNA-seq processing for the dataset

# RNA-flc
# same as esc

# chip -esc
# QCP
# align (samtools and bowtie2)
# peaks (macs2)
# chip qc 

# Hi-c
# EPI.R
# chromatin interactions using genomic interactions objects shows the distance between promoter-promoter interactions and promoter enhancer interactions for flc (and esc)




