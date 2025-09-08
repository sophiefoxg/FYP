#M2 vs M3
library(randomForest)
library(PRROC)
library(pROC)
#ESC

setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")
#clean and add together
#M2_esc<- read.delim("ESC_combined_M2.tsv", sep = "\t", header = TRUE)
#M3_esc<-read.delim("ESC_combined_M3.tsv", sep = "\t", header = TRUE)

head(M2_esc)
nrow(M2_esc)
nrow(M3_esc)
head(M3_esc)
#add CpG and motif to M2
#M2_esc$motif_similarity<-M3_esc$motif_similarity
#M2_esc$CpG_island<-M3_esc$CpG_island
#M2_esc$CpG_levels<-M3_esc$CpG_levels

head(M2_esc)
#combined_esc_final.df<-M2_esc

#seperate and add promoter id back in.
pe_esc.gi<-readRDS("ESC_samples_positive.rds")
negative_gi_esc<-readRDS("negative_interactions_esc.rds")

head(pe_esc.gi)
positives_esc <- subset(combined_esc_final.df, class == 1)
negatives_esc <- subset(combined_esc_final.df, class == 0)

positives_esc$promoter.id<-pe_esc.gi$anchor1.promoter.id
negatives_esc$promoter.id<-negative_gi_esc$anchor1.promoter.id

head(positives_esc)
head(negatives_esc)
#recombine
combined_esc_final.df <- rbind(positives_esc, negatives_esc)
setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")
#write.table(combined_esc_final.df, file = "ESC_combined_final.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

combined_esc_final.df <- read.delim("ESC_combined_final.tsv", sep = "\t", header = TRUE)
nrow(combined_esc_final.df)
combined_esc_final.df<-na.omit(combined_esc_final.df)
nrow(combined_esc_final.df)



## 
set.seed(42)
chromosomes <- unique(combined_esc_final.df$seqnames1)

# cross-chromosome validation
run_cv <- function(formula_obj) {
  all_labels <- c()
  all_predictions <- c()
  auc_values <- c()
  
  for (chr in chromosomes) {
    testdf  <- combined_esc_final.df[combined_esc_final.df$seqnames1 == chr, ]
    traindf <- combined_esc_final.df[combined_esc_final.df$seqnames1 != chr, ]
    
    model <- randomForest(formula_obj, data = traindf, ntree = 500)
    pred_probs <- predict(model, newdata = testdf, type = "prob")[, "1"]
    
    roc_obj <- roc(testdf$class, pred_probs)
    auc_values <- c(auc_values, auc(roc_obj))
    
    all_labels      <- c(all_labels, testdf$class)
    all_predictions <- c(all_predictions, pred_probs)
  }
  
  list(labels = all_labels, preds = all_predictions, auc = mean(auc_values))
}

#define model 2
f_Model2_esc <- as.formula(as.factor(class) ~ log10_distance + TPM +
                             H3K27acA1 + H3K27acA2 +
                             H3K4me1A1 + H3K4me1A2 +
                             H3K9acA1  + H3K9acA2 +
                             H3K4me3A1 + H3K4me3A2 +
                             H3K27me3A1 + H3K27me3A2 +
                             H3K9me3A1  + H3K9me3A2)
#define model 3
f_Model3_esc <- as.formula(as.factor(class) ~ log10_distance + TPM + motif_similarity + CpG_island)

f_Final_esc<-as.formula(as.factor(class) ~ log10_distance + TPM +
                          H3K27acA1 + H3K27acA2 +
                          H3K4me1A1 + H3K4me1A2 +
                          H3K9acA1  + H3K9acA2 +
                          H3K4me3A1 + H3K4me3A2 +
                          H3K27me3A1 + H3K27me3A2 +
                          H3K9me3A1  + H3K9me3A2 +motif_similarity+CpG_levels+CpG_island)

##run models
res_Model2_esc <- run_cv(f_Model2_esc)
res_Model3_esc <- run_cv(f_Model3_esc)
res_Final_esc<-run_cv((f_Final_esc))
##Pr curves
pr_Model2_esc <- pr.curve(
  scores.class0 = res_Model2_esc$preds[res_Model2_esc$labels == 1],
  scores.class1 = res_Model2_esc$preds[res_Model2_esc$labels == 0],
  curve = TRUE
)

pr_Model3_esc <- pr.curve(
  scores.class0 = res_Model3_esc$preds[res_Model3_esc$labels == 1],
  scores.class1 = res_Model3_esc$preds[res_Model3_esc$labels == 0],
  curve = TRUE
)
pr_Final_esc<- pr.curve(
  scores.class0 = res_Final_esc$preds[res_Final_esc$labels == 1],
  scores.class1 = res_Final_esc$preds[res_Final_esc$labels == 0],
  curve = TRUE
)
plot(pr_Model2_esc)
plot(pr_Model3_esc)
plot(pr_Final_esc)

## confusion matrices (ESC)
cutoff <- 0.5
pred_class_Model2_esc <- ifelse(res_Model2_esc$preds >= cutoff, 1, 0)
pred_class_Model3_esc <- ifelse(res_Model3_esc$preds >= cutoff, 1, 0)


print(table(True = res_Model2_esc$labels, Predicted = pred_class_Model2_esc))
print(table(True = res_Model3_esc$labels, Predicted = pred_class_Model3_esc))



## build a combined results frame
res_all <- data.frame(
  true   = res_Model2_esc$labels,  
  prob_M2 = res_Model2_esc$preds,
  prob_M3 = res_Model3_esc$preds,
  pred_M2 = pred_class_Model2_esc,
  pred_M3 = pred_class_Model3_esc
)


# add the features back in (assumes the order of rows is unchanged)
res_all <- cbind(res_all, combined_esc_final.df)
#make cpg binary
res_all$CpG_island <- ifelse(res_all$CpG_island == "CpG", 1, 0)
# False Positives
FP_M2 <- res_all[res_all$true == 0 & res_all$pred_M2 == 1, ]
FP_M3 <- res_all[res_all$true == 0 & res_all$pred_M3 == 1, ]

# False Negatives
FN_M2 <- res_all[res_all$true == 1 & res_all$pred_M2 == 0, ]
FN_M3 <- res_all[res_all$true == 1 & res_all$pred_M3 == 0, ]

# overlaps & uniques by row position
FP_both    <- res_all[res_all$true == 0 & res_all$pred_M2 == 1 & res_all$pred_M3 == 1, ]
FP_only_M2 <- res_all[res_all$true == 0 & res_all$pred_M2 == 1 & res_all$pred_M3 == 0, ]
FP_only_M3 <- res_all[res_all$true == 0 & res_all$pred_M2 == 0 & res_all$pred_M3 == 1, ]

FN_both    <- res_all[res_all$true == 1 & res_all$pred_M2 == 0 & res_all$pred_M3 == 0, ]
FN_only_M2 <- res_all[res_all$true == 1 & res_all$pred_M2 == 0 & res_all$pred_M3 == 1, ]
FN_only_M3 <- res_all[res_all$true == 1 & res_all$pred_M2 == 1 & res_all$pred_M3 == 0, ]

# quick counts
cat("False Positives: both =", nrow(FP_both),
    "only M2 =", nrow(FP_only_M2),
    "only M3 =", nrow(FP_only_M3), "\n")
cat("False Negatives: both =", nrow(FN_both),
    "only M2 =", nrow(FN_only_M2),
    "only M3 =", nrow(FN_only_M3), "\n")



# Add columns to flag correctness
res_all$M2_correct <- res_all$pred_M2 == res_all$true
res_all$M3_correct <- res_all$pred_M3 == res_all$true
# M2 correct, M3 incorrect
M2_correct_M3_wrong <- res_all[res_all$M2_correct == TRUE & res_all$M3_correct == FALSE, ]
# M3 correct, M2 incorrect
M3_correct_M2_wrong <- res_all[res_all$M3_correct == TRUE & res_all$M2_correct == FALSE, ]
# Both correct
Both_correct <- res_all[res_all$M2_correct == TRUE & res_all$M3_correct == TRUE, ]
# Both wrong
Both_wrong <- res_all[res_all$M2_correct == FALSE & res_all$M3_correct == FALSE, ]

# quick counts
cat("M2 correct, M3 wrong:", nrow(M2_correct_M3_wrong), "\n")
cat("M3 correct, M2 wrong:", nrow(M3_correct_M2_wrong), "\n")
cat("Both correct:", nrow(Both_correct), "\n")
cat("Both wrong:", nrow(Both_wrong), "\n")

# median values for key features
features <- c(
  "log10_distance","TPM","motif_similarity","CpG_island",
  "H3K27acA1","H3K27acA2",
  "H3K4me1A1","H3K4me1A2",
  "H3K9acA1","H3K9acA2",
  "H3K4me3A1","H3K4me3A2",
  "H3K27me3A1","H3K27me3A2",
  "H3K9me3A1","H3K9me3A2"
)


apply(M2_correct_M3_wrong[, features], 2, function(x) mean(as.numeric(x), na.rm=TRUE))
apply(M3_correct_M2_wrong[, features], 2, function(x) mean(as.numeric(x), na.rm=TRUE))

BW_FP  <- Both_wrong[Both_wrong$true == 0, ] 
BW_FN  <- Both_wrong[Both_wrong$true == 1, ]


BW_FP_med <- apply(BW_FP[, features, drop=FALSE], 2, function(x) median(as.numeric(x), na.rm=TRUE))
BW_FN_med <- apply(BW_FN[, features, drop=FALSE], 2, function(x) median(as.numeric(x), na.rm=TRUE))
BW_FP_mean <- apply(BW_FP[, features, drop=FALSE], 2, function(x) mean(as.numeric(x), na.rm=TRUE))
BW_FN_mean <- apply(BW_FN[, features, drop=FALSE], 2, function(x) mean(as.numeric(x), na.rm=TRUE))

cat("\nBoth wrong FP medians:\n"); print(BW_FP_med)
cat("\nBoth wrong FN medians:\n"); print(BW_FN_med)
cat("\nBoth wrong FP mean:\n"); print(BW_FP_mean)
cat("\nBoth wrong FN mean:\n"); print(BW_FN_mean)



## 2) summary stats for ALL vs BOTH_WRONG
ALL_med <- apply(res_all[, features, drop=FALSE], 2, function(x) median(as.numeric(x), na.rm=TRUE))
ALL_mean <- apply(res_all[, features, drop=FALSE], 2, function(x) mean(as.numeric(x), na.rm=TRUE))

BW_med  <- apply(Both_wrong[, features, drop=FALSE], 2, function(x) median(as.numeric(x), na.rm=TRUE))
BW_mean <- apply(Both_wrong[, features, drop=FALSE], 2, function(x) mean(as.numeric(x), na.rm=TRUE))

cat("\nALL — medians:\n"); print(ALL_med)
cat("\nBoth wrong — medians:\n"); print(BW_med)
cat("\nALL — means:\n");   print(ALL_mean)


cat("\nBoth wrong — means:\n");   print(BW_mean)

library(clusterProfiler)
library(org.Mm.eg.db)

##get the promoter IDs for each set
# flag correctness
res_all$M2_correct <- res_all$pred_M2 == res_all$true
res_all$M3_correct <- res_all$pred_M3 == res_all$true
# subsets
M2_correct_M3_wrong <- res_all[res_all$M2_correct & !res_all$M3_correct, ]
M3_correct_M2_wrong <- res_all[res_all$M3_correct & !res_all$M2_correct, ]

# extract promoter ids
genes_M2corr_M3wrong <- unique(M2_correct_M3_wrong$promoter.id)
genes_M3corr_M2wrong <- unique(M3_correct_M2_wrong$promoter.id)

##map to ENTREZ IDs
map_M2corr_M3wrong <- bitr(genes_M2corr_M3wrong,
                           fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"),
                           OrgDb=org.Mm.eg.db)
map_M3corr_M2wrong <- bitr(genes_M3corr_M2wrong,
                           fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"),
                           OrgDb=org.Mm.eg.db)

# clean duplicates and missing
map_M2corr_M3wrong <- map_M2corr_M3wrong[!duplicated(map_M2corr_M3wrong$ENSEMBL) & !is.na(map_M2corr_M3wrong$ENTREZID), ]
map_M3corr_M2wrong <- map_M3corr_M2wrong[!duplicated(map_M3corr_M2wrong$ENSEMBL) & !is.na(map_M3corr_M2wrong$ENTREZID), ]

##run GO enrichment (BP)
ego_M2corr_M3wrong <- enrichGO(
  gene          = map_M2corr_M3wrong$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

ego_M3corr_M2wrong <- enrichGO(
  gene          = map_M3corr_M2wrong$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

##make terms readable
ego_M2corr_M3wrong <- setReadable(ego_M2corr_M3wrong, OrgDb=org.Mm.eg.db, keyType="ENTREZID")
ego_M3corr_M2wrong <- setReadable(ego_M3corr_M2wrong, OrgDb=org.Mm.eg.db, keyType="ENTREZID")

## ---- check top terms ----
head(as.data.frame(ego_M2corr_M3wrong)[, c("Description")], 10)
head(as.data.frame(ego_M3corr_M2wrong)[, c("Description")], 10)

nrow(Both_wrong) # =697
length(unique(Both_wrong$promoter.id)) #494
##get promoter IDs for Both models wrong
genes_BothWrong <- unique(Both_wrong$promoter.id)


##map to ENTREZ IDs
map_BothWrong <- bitr(genes_BothWrong,
                      fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"),
                      OrgDb=org.Mm.eg.db)

# clean duplicates and missing
map_BothWrong <- map_BothWrong[!duplicated(map_BothWrong$ENSEMBL) & !is.na(map_BothWrong$ENTREZID), ]

##run GO enrichment (BP)
ego_BothWrong <- enrichGO(
  gene          = map_BothWrong$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

## make terms readable
ego_BothWrong <- setReadable(ego_BothWrong, OrgDb=org.Mm.eg.db, keyType="ENTREZID")

##check top terms
head(as.data.frame(ego_BothWrong)[, c("Description")], 10)

#### final model wrong, tells me what i could add in if i have more time


###### M2 vs M3 — FLC ######
library(randomForest)
library(PRROC)
library(pROC)

setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")
#### DONT RUN again
#load
M2_flc <- read.delim("FLC_combined_M2.tsv", sep = "\t", header = TRUE)
M3_flc <- read.delim("FLC_combined_M3.tsv", sep = "\t", header = TRUE)

# add CpG/motif from M3 into M2
M2_flc$motif_similarity <- M3_flc$motif_similarity
M2_flc$CpG_island       <- M3_flc$CpG_island
M2_flc$CpG_levels       <- M3_flc$CpG_levels

combined_flc_final.df <- M2_flc

#add promoter.id
pe_flc.gi       <- readRDS("FLC_samples_positive.rds")
negative_gi_flc <- readRDS("negative_interactions_flc.rds")

positives_flc <- subset(combined_flc_final.df, class == 1)
negatives_flc <- subset(combined_flc_final.df, class == 0)

positives_flc$promoter.id <- pe_flc.gi$anchor1.promoter.id
negatives_flc$promoter.id <- negative_gi_flc$anchor1.promoter.id

# recombine
combined_flc_final.df <- rbind(positives_flc, negatives_flc)

# write.table(combined_flc_final.df, file = "FLC_combined_final.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

combined_flc_final.df <- read.delim("FLC_combined_final.tsv", sep = "\t", header = TRUE)
combined_flc_final.df <- na.omit(combined_flc_final.df)

## rf
set.seed(42)
chromosomes <- unique(combined_flc_final.df$seqnames1)

# cross-chromosome validation 
run_cv <- function(formula_obj) {
  all_labels <- c()
  all_predictions <- c()
  auc_values <- c()
  for (chr in chromosomes) {
    testdf  <- combined_flc_final.df[combined_flc_final.df$seqnames1 == chr, ]
    traindf <- combined_flc_final.df[combined_flc_final.df$seqnames1 != chr, ]
    model <- randomForest(formula_obj, data = traindf, ntree = 500)
    pred_probs <- predict(model, newdata = testdf, type = "prob")[, "1"]
    roc_obj <- roc(testdf$class, pred_probs)
    auc_values <- c(auc_values, auc(roc_obj))
    all_labels      <- c(all_labels, testdf$class)
    all_predictions <- c(all_predictions, pred_probs)
  }
  list(labels = all_labels, preds = all_predictions, auc = mean(auc_values))
}

# define models
f_Model2_flc <- as.formula(as.factor(class) ~ log10_distance + TPM +
                             H3K27acA1 + H3K27acA2 +
                             H3K4me1A1 + H3K4me1A2 +
                             H3K9acA1  + H3K9acA2 +
                             H3K4me3A1 + H3K4me3A2 +
                             H3K27me3A1 + H3K27me3A2 +
                             H3K9me3A1  + H3K9me3A2)
#m3
f_Model3_flc <- as.formula(as.factor(class) ~ log10_distance + TPM + motif_similarity + CpG_island)
#final
f_Final_flc<-as.formula(as.factor(class) ~ log10_distance + TPM +
                          H3K27acA1 + H3K27acA2 +
                          H3K4me1A1 + H3K4me1A2 +
                          H3K9acA1  + H3K9acA2 +
                          H3K4me3A1 + H3K4me3A2 +
                          H3K27me3A1 + H3K27me3A2 +
                          H3K9me3A1  + H3K9me3A2 +motif_similarity+CpG_levels+CpG_island)
## run models 
res_Model2_flc <- run_cv(f_Model2_flc)
res_Model3_flc <- run_cv(f_Model3_flc)

##  PR curves 
pr_Model2_flc <- pr.curve(
  scores.class0 = res_Model2_flc$preds[res_Model2_flc$labels == 1],
  scores.class1 = res_Model2_flc$preds[res_Model2_flc$labels == 0],
  curve = TRUE
)
pr_Model3_flc <- pr.curve(
  scores.class0 = res_Model3_flc$preds[res_Model3_flc$labels == 1],
  scores.class1 = res_Model3_flc$preds[res_Model3_flc$labels == 0],
  curve = TRUE
)
pr_Final_flc <- pr.curve(
  scores.class0 = res_Final_flc$preds[res_Final_flc$labels == 1],
  scores.class1 = res_Final_flc$preds[res_Final_flc$labels == 0],
  curve = TRUE
)
plot(pr_Model2_flc)
plot(pr_Model3_flc)
plot(pr_Final_flc)

##confusion matrices (FLC)
cutoff <- 0.5
pred_class_Model2_flc <- ifelse(res_Model2_flc$preds >= cutoff, 1, 0)
pred_class_Model3_flc <- ifelse(res_Model3_flc$preds >= cutoff, 1, 0)

cat("Model 2 (FLC):\n")
print(table(True = res_Model2_flc$labels, Predicted = pred_class_Model2_flc))
cat("\nModel 3 (FLC):\n")
print(table(True = res_Model3_flc$labels, Predicted = pred_class_Model3_flc))

##build combined results
res_all_flc <- data.frame(
  true    = res_Model2_flc$labels,  
  prob_M2 = res_Model2_flc$preds,
  prob_M3 = res_Model3_flc$preds,
  pred_M2 = pred_class_Model2_flc,
  pred_M3 = pred_class_Model3_flc
)

# add features back in 
res_all_flc <- cbind(res_all_flc, combined_flc_final.df)
#turn Cpg binary for summary reports
res_all_flc$CpG_island <- ifelse(res_all_flc$CpG_island == "CpG", 1, 0)

#Errors
FP_M2_flc <- res_all_flc[res_all_flc$true == 0 & res_all_flc$pred_M2 == 1, ]
FP_M3_flc <- res_all_flc[res_all_flc$true == 0 & res_all_flc$pred_M3 == 1, ]
FN_M2_flc <- res_all_flc[res_all_flc$true == 1 & res_all_flc$pred_M2 == 0, ]
FN_M3_flc <- res_all_flc[res_all_flc$true == 1 & res_all_flc$pred_M3 == 0, ]

FP_both_flc    <- res_all_flc[res_all_flc$true == 0 & res_all_flc$pred_M2 == 1 & res_all_flc$pred_M3 == 1, ]
FP_only_M2_flc <- res_all_flc[res_all_flc$true == 0 & res_all_flc$pred_M2 == 1 & res_all_flc$pred_M3 == 0, ]
FP_only_M3_flc <- res_all_flc[res_all_flc$true == 0 & res_all_flc$pred_M2 == 0 & res_all_flc$pred_M3 == 1, ]

FN_both_flc    <- res_all_flc[res_all_flc$true == 1 & res_all_flc$pred_M2 == 0 & res_all_flc$pred_M3 == 0, ]
FN_only_M2_flc <- res_all_flc[res_all_flc$true == 1 & res_all_flc$pred_M2 == 0 & res_all_flc$pred_M3 == 1, ]
FN_only_M3_flc <- res_all_flc[res_all_flc$true == 1 & res_all_flc$pred_M2 == 1 & res_all_flc$pred_M3 == 0, ]

cat("False Positives (FLC): both =", nrow(FP_both_flc),
    " only M2 =", nrow(FP_only_M2_flc),
    " only M3 =", nrow(FP_only_M3_flc), "\n")
cat("False Negatives (FLC): both =", nrow(FN_both_flc),
    " only M2 =", nrow(FN_only_M2_flc),
    " only M3 =", nrow(FN_only_M3_flc), "\n")

# correct
res_all_flc$M2_correct <- res_all_flc$pred_M2 == res_all_flc$true
res_all_flc$M3_correct <- res_all_flc$pred_M3 == res_all_flc$true

M2_correct_M3_wrong_flc <- res_all_flc[ res_all_flc$M2_correct & !res_all_flc$M3_correct, ]
M3_correct_M2_wrong_flc <- res_all_flc[ res_all_flc$M3_correct & !res_all_flc$M2_correct, ]
Both_correct_flc        <- res_all_flc[ res_all_flc$M2_correct &  res_all_flc$M3_correct, ]
Both_wrong_flc          <- res_all_flc[!res_all_flc$M2_correct & !res_all_flc$M3_correct, ]

cat("M2 correct, M3 wrong (FLC):", nrow(M2_correct_M3_wrong_flc), "\n")
cat("M3 correct, M2 wrong (FLC):", nrow(M3_correct_M2_wrong_flc), "\n")
cat("Both correct (FLC):", nrow(Both_correct_flc), "\n")
cat("Both wrong (FLC):", nrow(Both_wrong_flc), "\n")

#feature summaries
features <- c(
  "log10_distance","TPM","motif_similarity","CpG_island", "H3K27acA1","H3K27acA2","H3K4me1A1","H3K4me1A2","H3K9acA1","H3K9acA2","H3K4me3A1","H3K4me3A2","H3K27me3A1","H3K27me3A2","H3K9me3A1","H3K9me3A2"
)

apply(M2_correct_M3_wrong_flc[, features, drop=FALSE], 2, function(x) mean(as.numeric(x), na.rm=TRUE))
apply(M3_correct_M2_wrong_flc[, features, drop=FALSE], 2, function(x) mean(as.numeric(x), na.rm=TRUE))

BW_FP_flc  <- Both_wrong_flc[Both_wrong_flc$true == 0, ]   # shared FP
BW_FN_flc  <- Both_wrong_flc[Both_wrong_flc$true == 1, ]   # shared FN

BW_overall_med_flc <- apply(Both_wrong_flc[, features, drop=FALSE], 2, function(x) median(as.numeric(x), na.rm=TRUE))
BW_FP_med_flc <- apply(BW_FP_flc[, features, drop=FALSE], 2, function(x) median(as.numeric(x), na.rm=TRUE))
BW_FN_med_flc <- apply(BW_FN_flc[, features, drop=FALSE], 2, function(x) median(as.numeric(x), na.rm=TRUE))
BW_FP_mean_flc <- apply(BW_FP_flc[, features, drop=FALSE], 2, function(x) mean(as.numeric(x), na.rm=TRUE))
BW_FN_mean_flc <- apply(BW_FN_flc[, features, drop=FALSE], 2, function(x) mean(as.numeric(x), na.rm=TRUE))

cat("\nBoth wrong (FLC) — overall medians:\n"); print(BW_overall_med_flc)
cat("\nShared FP (FLC) medians:\n"); print(BW_FP_med_flc)
cat("\nShared FN (FLC) medians:\n"); print(BW_FN_med_flc)
cat("\nShared FP (FLC) means:\n"); print(BW_FP_mean_flc)
cat("\nShared FN (FLC) means:\n"); print(BW_FN_mean_flc)

## overall dataswet vs BOTH_WRONG  summaries
ALL_med_flc  <- apply(res_all_flc[, features, drop=FALSE], 2, function(x) median(as.numeric(x), na.rm=TRUE))
ALL_mean_flc <- apply(res_all_flc[, features, drop=FALSE], 2, function(x) mean(as.numeric(x), na.rm=TRUE))
BW_med_flc   <- apply(Both_wrong_flc[, features, drop=FALSE], 2, function(x) median(as.numeric(x), na.rm=TRUE))
BW_mean_flc  <- apply(Both_wrong_flc[, features, drop=FALSE], 2, function(x) mean(as.numeric(x), na.rm=TRUE))

cat("\nALL (FLC) — medians:\n"); print(ALL_med_flc)
cat("\nBoth wrong (FLC) — medians:\n"); print(BW_med_flc)
cat("\nALL (FLC) — means:\n");   print(ALL_mean_flc)
cat("\nBoth wrong (FLC) — means:\n");   print(BW_mean_flc)

# GO enrichment
library(clusterProfiler)
library(org.Mm.eg.db)

# M2-correct, M3-wrong 
genes_M2corr_M3wrong_flc <- unique(M2_correct_M3_wrong_flc$promoter.id)
genes_M2corr_M3wrong_flc <- genes_M2corr_M3wrong_flc[!is.na(genes_M2corr_M3wrong_flc) & genes_M2corr_M3wrong_flc != "not_promoter"]

map_M2corr_M3wrong_flc <- bitr(genes_M2corr_M3wrong_flc,
                               fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"),
                               OrgDb=org.Mm.eg.db)
map_M2corr_M3wrong_flc <- map_M2corr_M3wrong_flc[!duplicated(map_M2corr_M3wrong_flc$ENSEMBL) & !is.na(map_M2corr_M3wrong_flc$ENTREZID), ]

ego_M2corr_M3wrong_flc <- enrichGO(
  gene          = map_M2corr_M3wrong_flc$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)
ego_M2corr_M3wrong_flc <- setReadable(ego_M2corr_M3wrong_flc, OrgDb=org.Mm.eg.db, keyType="ENTREZID")
head(as.data.frame(ego_M2corr_M3wrong_flc)[, "Description", drop=FALSE], 10)

#M3-correct, M2-wrong 
genes_M3corr_M2wrong_flc <- unique(M3_correct_M2_wrong_flc$promoter.id)
genes_M3corr_M2wrong_flc <- genes_M3corr_M2wrong_flc[!is.na(genes_M3corr_M2wrong_flc) & genes_M3corr_M2wrong_flc != "not_promoter"]

map_M3corr_M2wrong_flc <- bitr(genes_M3corr_M2wrong_flc,
                               fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"),
                               OrgDb=org.Mm.eg.db)
map_M3corr_M2wrong_flc <- map_M3corr_M2wrong_flc[!duplicated(map_M3corr_M2wrong_flc$ENSEMBL) & !is.na(map_M3corr_M2wrong_flc$ENTREZID), ]

ego_M3corr_M2wrong_flc <- enrichGO(
  gene          = map_M3corr_M2wrong_flc$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)
ego_M3corr_M2wrong_flc <- setReadable(ego_M3corr_M2wrong_flc, OrgDb=org.Mm.eg.db, keyType="ENTREZID")
head(as.data.frame(ego_M3corr_M2wrong_flc)[, "Description", drop=FALSE], 10)

#Both_wrong  3080
cat("\nBoth_wrong (FLC) size =", nrow(Both_wrong_flc), 
    " | unique promoters =", length(unique(Both_wrong_flc$promoter.id)), "\n")
#1653 unique promoters
genes_BothWrong_flc <- unique(Both_wrong_flc$promoter.id)
genes_BothWrong_flc <- genes_BothWrong_flc[!is.na(genes_BothWrong_flc) & genes_BothWrong_flc != "not_promoter"]

map_BothWrong_flc <- bitr(genes_BothWrong_flc,
                          fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"),
                          OrgDb=org.Mm.eg.db)
map_BothWrong_flc <- map_BothWrong_flc[!duplicated(map_BothWrong_flc$ENSEMBL) & !is.na(map_BothWrong_flc$ENTREZID), ]

ego_BothWrong_flc <- enrichGO(
  gene          = map_BothWrong_flc$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)
ego_BothWrong_flc <- setReadable(ego_BothWrong_flc, OrgDb=org.Mm.eg.db, keyType="ENTREZID")
head(as.data.frame(ego_BothWrong_flc)[, "Description", drop=FALSE], 10)
