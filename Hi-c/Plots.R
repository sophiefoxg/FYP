setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")
###motif enrichment
#install packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomicInteractions")
library(GenomicInteractions)


#setwd("C:/Users/foxgs/OneDrive/Documents/pos_neg_sets")

pe_esc.gi<-readRDS("ESC_samples_positive.rds")
negative_gi_esc<-readRDS("negative_interactions_esc.rds")
negative_nodis_esc<-readRDS("negative_interactions_no_distance_esc.rds")

pe_flc.gi <- readRDS("FLC_samples_positive.rds")
negative_gi_flc <- readRDS("negative_interactions_flc.rds")
negative_nodis_flc <- readRDS("negative_interactions_no_distance_flc.rds")


negative_gi_flc$distance<-calculateDistances(negative_gi_flc)
head(negative_gi_flc)

median(pe_esc.gi$distance)
median(pe_flc.gi$distance)

median(negative_nodis_esc$distance)
median(negative_nodis_flc$distance)

median(negative_gi_esc$distance)
median(negative_gi_flc$distance)

hist(log10(pe_esc.gi$distance), breaks=50,
     col="skyblue", main="ESC: Distance Distribution",
     xlab="log10 Distance (bp)")

hist(log10(pe_flc.gi$distance), breaks=50,
     col="salmon", main="FLC: Distance Distribution",
     xlab="log10 Distance (bp)")

head(negative_nodis_esc)
#ngeative no dis
hist(log10(negative_nodis_esc$distance), breaks=50,
     col="skyblue", main="ESC: No Distance Matching Distribution",
     xlab="log10 Distance (bp)")

hist(log10(negative_nodis_flc$distance), breaks=50,
     col="salmon", main="FLC: No Distance Matching Distribution",
     xlab="log10 Distance (bp)")



#ngeative distance matching
hist(log10(negative_gi_esc$distance), breaks=50,
     col="skyblue", main="ESC: Distance Matching Distribution",
     xlab="log10 Distance (bp)")

hist(log10(negative_gi_flc$distance), breaks=50,
     col="salmon", main="FLC: Distance Matching Distribution",
     xlab="log10 Distance (bp)")

#correlation between distance and promoter expression:
cor(pe_esc.gi$TPM, log10(pe_esc.gi$distance), use = "complete.obs", method = "pearson")

cor(pe_flc.gi$TPM, log10(pe_flc.gi$distance), use = "complete.obs", method = "pearson")


cor(pe_esc.gi$H3K27ac, pe_esc.gi$TPM, use = "complete.obs")
cor(pe_flc.gi$H3K27ac, pe_flc.gi$TPM, use = "complete.obs")

#auc and pr curves
install.packages("randomForest")
install.packages("pROC")
install.packages("PRROC")
library(randomForest)
library(pROC)
library(PRROC)

setwd("C:/Users/foxgs/OneDrive/Documents/combineddfs")
#combined df of positives and negatives
ESC_df <- read.delim("ESC_samples_distance_match.tsv")
FLC_df <- read.delim("FLC_samples_distance_match.tsv")


ESC_nodis_df <- read.delim("ESC_samples_no_distance_match.tsv")

FLC_nodis_df <- read.delim("FLC_samples_no_distance_match.tsv")

combined.df<-FLC_nodis_df

# Get unique chromosomes from data frame
chromosomes <- unique(combined.df$seqnames1)

set.seed(42)
all_labels <- c()
all_predictions <- c()
auc_values <- c()
#select test chr and train on all others
for (chr in chromosomes) {
  # Split train/test based on current chromosome
  testdf  <- combined.df[combined.df$seqnames1 == chr, ]
  traindf <- combined.df[combined.df$seqnames1 != chr, ]
  #train model on all other chromosomes  
  model <- randomForest(as.factor(class) ~ log10_distance + TPM ,
                        data = traindf, ntree = 500)
  #predict probabilities on the test chromosome
  pred_probs <- predict(model, newdata = testdf, type = "prob")[,"1"]
  roc_obj <- roc(testdf$class, pred_probs)
  auc_values <- c(auc_values, auc(roc_obj))
  
  # save predictions and labels as list for ROC
  all_labels <- c(all_labels, testdf$class)
  all_predictions <- c(all_predictions, pred_probs)
}

#average auc performance across chromosomes
mean_auc <- mean(auc_values)
print(auc_values)
#plot avg ROC using all predictions 
plot(roc(all_labels, all_predictions))
roc(all_labels, all_predictions)

# plot precision-recall curve averaged


pr <- pr.curve(scores.class0 = all_predictions[all_labels == 1],
               scores.class1 = all_predictions[all_labels == 0],
               curve = TRUE)
plot(pr)

