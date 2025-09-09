# Predicting Promoter-Enhancer Interactions using Machine Learning #

## Promoter–Enhancer Interaction Analysis

This repository contains R scripts for processing Promoter Capture Hi-C (PCHi-C) data from Schoenfelder et al. (2015) to identify promoter–enhancer and promoter–promoter interactions in mouse embryonic stem cells (ESC) and fetal liver cells (FLC).
The pipeline integrates PCHi-C with ChIP-seq (H3K27ac broad peaks) and RNA-seq (TPM expression) to build annotated interaction datasets and train random forest models for predicting enhancer–promoter interactions.

## Project overview:
<img width="1094" height="853" alt="image" src="https://github.com/user-attachments/assets/333b4a74-82a8-4f50-9ef8-b8ee3125463a" />

# Pre-processing
All pre-processing steps were performed in a **MobaXterm** terminal within a **SLURM** executable environment.


## Pipeline for RNA-sequencing
 #### found in ./RNA/esc
RNA-seq data available:
ESC (GSM723776) and FLC (GSM661638) from Shen et al (Shen et al. 2012) 
MM9 reference genome UCSC
### Scripts
1fastQCP.slurm 
2star.slurm 
3featurecounts.slurm  
TPM.R
#### Requirements:
fastp v0.23.4, FastQC v0.12.1, STAR v2.7, subread/featureCounts v2.0.6



## Pipeline for analysing chip-sequencing data as histone marks in esc and flc:
#### found in ./chip/esc
Chip-seq data available at:
ESC: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29184
FLC:GSE153563
### Scripts
1QCP.slurm  
2align.slurm 
3peaks.slurm
#### Requirements:
FastQC v0.12.1, fastp v0.23.4, Bowtie2 v2.5.2, samtools v1.19.2, MACS2 v2.2.8, bedtools v2.31.1

## CpG
CpG island annotation step:done in command line
#### Requirements
bedtools v2.31.1


# Model 1

### found in ./Hi-C
### EPI.R: Identifying true EPIs and sampling negatives
Required to clean PCHi-C data from Schoenfelder et al. (2015) in ./Hi-c folder:

- `ESC_promoter_other_significant_interactions.txt`  
- `ESC_promoter_promoter_significant_interactions.txt`  
- `FLC_promoter_other_significant_interactions.txt`  
- `FLC_promoter_promoter_significant_interactions.txt`  
### Sampling
Two types of sampling:  
1. Not matching for distance  
2. Matched for distance  

##  Results for Model 1 and example of LOCO
#### Requirements
- BiocManager v1.30.25  , GenomicInteractions v1.40.0 (Harmston et al. 2015) , GenomicRanges v1.58.0 (Lawrence et al. 2013) , IRanges v2.40.0 (Lawrence et al. 2013),  rtracklayer v1.66.0 (Lawrence et al. 2009), randomForest v4.7-1.1 (Liaw & Wiener 2002)  
ggplot2 v3.5.1 (Wickham 2016), cowplot v1.1.3 (Wilke 2020), pROC v1.18.5 (Robin et al. 2011), PRROC v1.3.1 (Grau et al. 2015) 



# Model 2

#### found in: ./Hi-c
### chromatin_marks_esc.R
Overlaps genomic coordinates of Positives and negatives with histone marks found in pre-processing steps, resulting in signal values for each promoter, enhancer and "negative"

##### Required packages (plus EPI.R packages):
IRanges v2.4


# Model 3

#### found in: ./Hi-c
### motif_esc.R
Identifies motifs and calculates similarity scores between positive and negative set

##### Required packages (plus EPI.R packages):
TFBSTools v1.44, lsa v0.73, and the mm9 genome v1.4.


# Final Model
#### found in: ./Hi-c
### Final_model.R
Integrates all features and plots all models used in the report
##### Required packages = EPI.R packages


# Cross Cell Type
#### found in: ./Hi-c
### cell_type_comparison.R
Compares the final models performance between esc and flc cells
##### Required packages = EPI.R packages

# Additional Scripts:
#### found in: ./Hi-c
### Plots.R
Plots the majority of figures seen in the project report

#### found in: ./Hi-c
### M2vsM3.R
the beginning of error analysis for model 2 and 3

