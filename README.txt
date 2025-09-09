####FYP predicting enhancer promoter interactions ####

#./Hi-C
##EPI.R: Identifying true EPIs and sampling negatives
Required to clean PCHi-C data from Schoenfelder et al. (2015) in ./Hi-c folder:

- `ESC_promoter_other_significant_interactions.txt`  
- `ESC_promoter_promoter_significant_interactions.txt`  
- `FLC_promoter_other_significant_interactions.txt`  
- `FLC_promoter_promoter_significant_interactions.txt`  
##Sampling

Two types of sampling:  
1. Not matching for distance  
2. Matched for distance  

## Requirements

- BiocManager v1.30.25  , GenomicInteractions v1.40.0 (Harmston et al. 2015) , GenomicRanges v1.58.0 (Lawrence et al. 2013) , IRanges v2.40.0 (Lawrence et al. 2013),  rtracklayer v1.66.0 (Lawrence et al. 2009), randomForest v4.7-1.1 (Liaw & Wiener 2002)  
ggplot2 v3.5.1 (Wickham 2016), cowplot v1.1.3 (Wilke 2020), pROC v1.18.5 (Robin et al. 2011), PRROC v1.3.1 (Grau et al. 2015) 
# Model 1 #

## R script EPI.R ##
# Promoter–Enhancer Interaction Modeling (ESC vs FLC)
This R script builds **positive** promoter–enhancer (PE) interaction sets from Hi-C data, constructs **negative control** sets (window-based and distance-matched), annotates features (**distance**, **H3K27ac**, **CpG**, **RNA TPM**), and trains a simple **random-forest** classifier with cross-chromosome validation.  
It runs the full pipeline for **ESC** and **FLC** samples.

---

## TL;DR

### 1. Load significant Hi-C interactions
- **Promoter–other**
- **Promoter–promoter**

### 2. Call enhancers
- Identify “other” anchors overlapping **H3K27ac** peaks
- Keep only **promoter–enhancer** pairs

### 3. Annotate PE pairs
- Genomic distance
- **H3K27ac** signal (mean per overlapping peak)
- **CpG island** overlap
- **RNA TPM** for the promoter

### 4. Generate negatives
- **Window-based**: random HindIII fragments within ±1.5 Mb of the promoter, excluding known enhancers
- **Distance-matched**: random fragments sampled to match the positive distances

### 5. Train random forest (RF)
- **Features**: `log10(distance)` + `TPM`
- **Evaluation**: AUC and PR curves via **leave-one-chromosome-out** cross-validation

## Required Packages

### Core R Packages
- **randomForest** – train and evaluate random forest models  
- **pROC** – generate and plot ROC curves  
- **PRROC** – generate and plot precision–recall curves  
- **ggplot2** – general plotting  
- **cowplot** – arrange multiple plots into a grid  

### Bioconductor Packages
- **GenomicInteractions** – store, manipulate, and analyze genomic interaction data  
- **GenomicRanges** – represent genomic intervals and their annotations  
- **rtracklayer** – import/export genomic data formats (e.g., BED, GTF)


