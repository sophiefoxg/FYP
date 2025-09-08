# Cross cell type comparison
library(randomForest)
library(pROC)
library(PRROC)

set.seed(42)

# Data
setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")
ESC_df <- read.delim("ESC_combined_final.tsv")
FLC_df <- read.delim("FLC_combined_final.tsv")

ESC_df<-na.omit(ESC_df)
FLC_df<-na.omit(FLC_df)

## ========= ESC train → FLC test =========
shared_chr <- intersect(unique(ESC_df$seqnames1), unique(FLC_df$seqnames1))

all_labels <- c()
all_predictions <- c()
auc_values <- c()

for (chr in shared_chr) {
  # train on ESC excluding chr; test on that chr in FLC
  traindf <- ESC_df[ESC_df$seqnames1 != chr, ]
  testdf  <- FLC_df[FLC_df$seqnames1 == chr, ]
  
  model <- randomForest(
    as.factor(class) ~ log10_distance + TPM + motif_similarity + CpG_island + CpG_levels +
      H3K27acA2 + H3K27acA1 + H3K4me1A1 + H3K4me1A2 + H3K9acA1 + H3K9acA2 +
      H3K4me3A1 + H3K4me3A2 + H3K27me3A1 + H3K27me3A2 + H3K9me3A1 + H3K9me3A2,
    data = traindf, ntree = 500
  )
  
  pred_probs <- predict(model, newdata = testdf, type = "prob")[, "1"]
  roc_obj <- roc(testdf$class, pred_probs)
  
  auc_values      <- c(auc_values, auc(roc_obj))
  all_labels      <- c(all_labels, testdf$class)
  all_predictions <- c(all_predictions, pred_probs)
}

mean_auc_ESC_to_FLC <- mean(auc_values)
print(mean_auc_ESC_to_FLC)

# PR curve (renamed)
pr_ESC_to_FLC <- pr.curve(
  scores.class0 = all_predictions[all_labels == 1],
  scores.class1 = all_predictions[all_labels == 0],
  curve = TRUE
)
plot(pr_ESC_to_FLC, main = "ESC-trained → FLC-tested (PR)")



## ========= FLC train → ESC test =========
shared_chr <- intersect(unique(FLC_df$seqnames1), unique(ESC_df$seqnames1))

all_labels <- c()
all_predictions <- c()
auc_values <- c()

for (chr in shared_chr) {
  # train on FLC excluding chr; test on that chr in ESC
  traindf <- FLC_df[FLC_df$seqnames1 != chr, ]
  testdf  <- ESC_df[ESC_df$seqnames1 == chr, ]
  
  model <- randomForest(
    as.factor(class) ~ log10_distance + TPM + motif_similarity + CpG_island + CpG_levels +
      H3K27acA2 + H3K27acA1 + H3K4me1A1 + H3K4me1A2 + H3K9acA1 + H3K9acA2 +
      H3K4me3A1 + H3K4me3A2 + H3K27me3A1 + H3K27me3A2 + H3K9me3A1 + H3K9me3A2,
    data = traindf, ntree = 500
  )
  
  pred_probs <- predict(model, newdata = testdf, type = "prob")[, "1"]
  roc_obj <- roc(testdf$class, pred_probs)
  
  auc_values      <- c(auc_values, auc(roc_obj))
  all_labels      <- c(all_labels, testdf$class)
  all_predictions <- c(all_predictions, pred_probs)
}

mean_auc_FLC_to_ESC <- mean(auc_values)
print(mean_auc_FLC_to_ESC)

# ROC & PR plots
plot(roc(all_labels, all_predictions), main = "FLC-trained → ESC-tested (ROC)")

pr_FLC_to_ESC <- pr.curve(
  scores.class0 = all_predictions[all_labels == 1],
  scores.class1 = all_predictions[all_labels == 0],
  curve = TRUE
)
plot(pr_FLC_to_ESC, main = "FLC-trained → ESC-tested (PR)")




