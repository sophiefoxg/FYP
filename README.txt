# FYP predicting enhancer promoter interactions 

# ./Hi-C
## EPI.R: Identifying true EPIs and sampling negatives
Required to clean PCHi-C data from Schoenfelder et al. (2015) in ./Hi-c folder:

- `ESC_promoter_other_significant_interactions.txt`  
- `ESC_promoter_promoter_significant_interactions.txt`  
- `FLC_promoter_other_significant_interactions.txt`  
- `FLC_promoter_promoter_significant_interactions.txt`  
## Sampling
Two types of sampling:  
1. Not matching for distance  
2. Matched for distance  

## Model 1 and example of LOCO
## Requirements
- BiocManager v1.30.25  , GenomicInteractions v1.40.0 (Harmston et al. 2015) , GenomicRanges v1.58.0 (Lawrence et al. 2013) , IRanges v2.40.0 (Lawrence et al. 2013),  rtracklayer v1.66.0 (Lawrence et al. 2009), randomForest v4.7-1.1 (Liaw & Wiener 2002)  
ggplot2 v3.5.1 (Wickham 2016), cowplot v1.1.3 (Wilke 2020), pROC v1.18.5 (Robin et al. 2011), PRROC v1.3.1 (Grau et al. 2015) 





