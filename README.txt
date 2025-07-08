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
