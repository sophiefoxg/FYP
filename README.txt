#### FYP predicting enhancer promoter interactions ####

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


