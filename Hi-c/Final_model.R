setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")

combined_esc_final.df <- read.delim("ESC_combined_final.tsv", sep = "\t", header = TRUE)
combined_esc_final.df<-na.omit(combined_esc_final.df)

combined_esc_final.df <- combined_esc_final.df[ , !(names(combined_esc_final.df) %in% "CpG_levels")]


head(combined_esc_final.df)
#RF plots
library(randomForest)
library(pROC)
library(PRROC)

set.seed(42)
chromosomes <- unique(combined_esc_final.df$seqnames1)

all_labels <- c()
all_predictions <- c()
auc_values <- c()
#select test chr and train on all others
for (chr in chromosomes) {
  # Split train/test based on current chromosome
  testdf  <- combined_esc_final.df[combined_esc_final.df$seqnames1 == chr, ]
  traindf <- combined_esc_final.df[combined_esc_final.df$seqnames1 != chr, ]
  #train model on all other chromosomes  
  #
  model <- randomForest(as.factor(class) ~ log10_distance + TPM  + motif_similarity+ CpG_island+H3K27acA2+H3K27acA1 +H3K4me1A1 +H3K4me1A2 +H3K9acA1+H3K9acA2+H3K4me3A1+H3K4me3A2+H3K27me3A1+H3K27me3A2+H3K9me3A1+H3K9me3A2,
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

#dev.off()
#plot avg ROC using all predictions 
plot(roc(all_labels, all_predictions))
roc(all_labels, all_predictions)

# plot precision-recall curve averaged
pr <- pr.curve(scores.class0 = all_predictions[all_labels == 1],
               scores.class1 = all_predictions[all_labels == 0],
               curve = TRUE)
plot(pr)

varImpPlot(model)
# importance
imp <- importance(model, type = 2)  # MeanDecreaseGini
imp_df <- data.frame(Feature = rownames(imp), Importance = imp[, 1])

#rename features 
imp_df$Feature <- gsub("^TPM$", "Expression", imp_df$Feature)
imp_df$Feature <- gsub("^log10_distance$", "Distance", imp_df$Feature)
imp_df$Feature <- gsub("^motif_similarity$", "Motif Similarity", imp_df$Feature)
imp_df$Feature <- gsub("^CpG_island$", "CpG binary", imp_df$Feature)
imp_df$Feature <- gsub("A2$", " (E)", imp_df$Feature)
imp_df$Feature <- gsub("A1$", " (P)", imp_df$Feature)

#sort by importance
imp_df <- imp_df[order(imp_df$Importance, decreasing = TRUE), ]

library(ggplot2)
#plot 
ggplot(imp_df, aes(x = Importance, y = reorder(Feature, Importance))) +
  geom_col(fill = "skyblue") +
  labs(title = "Feature Importance Final Model ",
       x = "Mean Decrease Gini",
       y = "Feature") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 12)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)))


#for plotting define model formulas
#baseline model
f_Base_esc <- as.formula(as.factor(class) ~ log10_distance + TPM)

#model2
f_Model2_esc <- as.formula(as.factor(class) ~ log10_distance + TPM +
                             H3K27acA1 + H3K27acA2 +
                             H3K4me1A1 + H3K4me1A2 +
                             H3K9acA1  + H3K9acA2 +
                             H3K4me3A1 + H3K4me3A2 +
                             H3K27me3A1 + H3K27me3A2 +
                             H3K9me3A1  + H3K9me3A2)

#model3
f_Model3_esc <- as.formula(as.factor(class) ~ log10_distance + TPM + motif_similarity + CpG_island)

#model final
f_Final_esc <- as.formula(as.factor(class) ~ log10_distance + TPM +
                            H3K27acA1 + H3K27acA2 +
                            H3K4me1A1 + H3K4me1A2 +
                            H3K9acA1  + H3K9acA2 +
                            H3K4me3A1 + H3K4me3A2 +
                            H3K27me3A1 + H3K27me3A2 +
                            H3K9me3A1  + H3K9me3A2 +
                            motif_similarity + CpG_island)

#final no enhancer mark
f_Final_noH3K27acA2_esc <- as.formula(as.factor(class) ~ log10_distance + TPM +
                                        H3K27acA1 +
                                        H3K4me1A1 + H3K4me1A2 +
                                        H3K9acA1  + H3K9acA2 +
                                        H3K4me3A1 + H3K4me3A2 +
                                        H3K27me3A1 + H3K27me3A2 +
                                        H3K9me3A1  + H3K9me3A2 +
                                        motif_similarity + CpG_island)

#cross chromosome validation loop
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
    
    all_labels <- c(all_labels, testdf$class)
    all_predictions <- c(all_predictions, pred_probs)
  }
  
  list(labels = all_labels, preds = all_predictions, auc = mean(auc_values))
}

res_Base_esc              <- run_cv(f_Base_esc)
res_Model2_esc            <- run_cv(f_Model2_esc)
res_Model3_esc            <- run_cv(f_Model3_esc)
res_Final_esc             <- run_cv(f_Final_esc)
res_Final_noH3K27acA2_esc <- run_cv(f_Final_noH3K27acA2_esc)

#create pr obj
pr_Base_esc <- pr.curve(scores.class0 = res_Base_esc$preds[res_Base_esc$labels == 1],
                        scores.class1 = res_Base_esc$preds[res_Base_esc$labels == 0],
                        curve = TRUE)

pr_Model2_esc <- pr.curve(scores.class0 = res_Model2_esc$preds[res_Model2_esc$labels == 1],
                          scores.class1 = res_Model2_esc$preds[res_Model2_esc$labels == 0],
                          curve = TRUE)

pr_Model3_esc <- pr.curve(scores.class0 = res_Model3_esc$preds[res_Model3_esc$labels == 1],
                          scores.class1 = res_Model3_esc$preds[res_Model3_esc$labels == 0],
                          curve = TRUE)

pr_Final_esc <- pr.curve(scores.class0 = res_Final_esc$preds[res_Final_esc$labels == 1],
                         scores.class1 = res_Final_esc$preds[res_Final_esc$labels == 0],
                         curve = TRUE)

pr_Final_noH3K27acA2_esc <- pr.curve(scores.class0 = res_Final_noH3K27acA2_esc$preds[res_Final_noH3K27acA2_esc$labels == 1],
                                     scores.class1 = res_Final_noH3K27acA2_esc$preds[res_Final_noH3K27acA2_esc$labels == 0],
                                     curve = TRUE)

#model 2 plot vs baseline
plot(pr_Model2_esc, main = "ESC: Model 2", lwd = 2, col = "yellow")
plot(pr_Base_esc,   add = TRUE, lwd = 2, col = "red")
legend("bottomleft",
       legend = c(
         sprintf("Model 2 (AUC = %.3f)",  pr_Model2_esc$auc.integral),
         sprintf("Baseline (AUC = %.3f)", pr_Base_esc$auc.integral)
       ),
       col = c("yellow","red"), lwd = 2, bty = "n")

#model 3 plot vs baseline
plot(pr_Model3_esc, main = "ESC: Model 3", lwd = 2, col = "green")
plot(pr_Base_esc,   add = TRUE, lwd = 2, col = "red")
legend("bottomleft",
       legend = c(
         sprintf("Model 3 (AUC = %.3f)",  pr_Model3_esc$auc.integral),
         sprintf("Baseline (AUC = %.3f)", pr_Base_esc$auc.integral)
       ),
       col = c("green","red"), lwd = 2, bty = "n")

#final model vs baseline vs noH3K27acA2
dev.off()
plot(pr_Final_esc,             main = "ESC: Final Model", lwd = 2, col = "green")
plot(pr_Final_noH3K27acA2_esc, add = TRUE, lwd = 2, col = "blue")
plot(pr_Base_esc,              add = TRUE, lwd = 2, col = "red")
legend("bottomleft",
       legend = c(
         sprintf("Final (AUC = %.3f)",            pr_Final_esc$auc.integral),
         sprintf("Without H3K27acA2 (AUC = %.3f)", pr_Final_noH3K27acA2_esc$auc.integral),
         sprintf("Baseline (AUC = %.3f)",         pr_Base_esc$auc.integral)
       ),
       col = c("green","blue","red"), lwd = 2, bty = "n")

#final model vs model2 vs model3 vs baseline
plot(pr_Final_esc,  main = "ESC:Final Model", lwd = 2, col = "green")
plot(pr_Final_noH3K27acA2_esc,add = TRUE, lwd = 2, col = "blue")
plot(pr_Model2_esc, add = TRUE, lwd = 2, col = "yellow")
plot(pr_Model3_esc, add = TRUE, lwd = 2, col = "orange")
plot(pr_Base_esc,   add = TRUE, lwd = 2, col = "red")
plot.new()
legend("bottomleft",
       legend = c(
         sprintf("Final (AUC = %.3f)",   pr_Final_esc$auc.integral),
         sprintf("No H3K27ac (AUC = %.3f)", pr_Final_noH3K27acA2_esc$auc.integral),
         sprintf("Model 2 (AUC = %.3f)", pr_Model2_esc$auc.integral),
         sprintf("Model 3 (AUC = %.3f)", pr_Model3_esc$auc.integral),
         sprintf("Baseline (AUC = %.3f)",pr_Base_esc$auc.integral)
       ),
       col = c("green","blue","yellow","orange","red"), lwd = 2, bty = "n")




#####FLC

setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")

combined_flc_final.df <- read.delim("FLC_combined_final.tsv", sep = "\t", header = TRUE)
combined_flc_final.df<-na.omit(combined_flc_final.df)

combined_flc_final.df <- combined_flc_final.df[ , !(names(combined_flc_final.df) %in% "CpG_levels")]
head(combined_flc_final.df)
# RF plots
library(randomForest)
library(pROC)
library(PRROC)

set.seed(42)
chromosomes <- unique(combined_flc_final.df$seqnames1)

all_labels <- c()
all_predictions <- c()
auc_values <- c()
# select test chr and train on all others
for (chr in chromosomes) {
  # Split train/test based on current chromosome
  testdf  <- combined_flc_final.df[combined_flc_final.df$seqnames1 == chr, ]
  traindf <- combined_flc_final.df[combined_flc_final.df$seqnames1 != chr, ]
  # train model on all other chromosomes  
  #
  model <- randomForest(as.factor(class) ~ log10_distance + TPM  + motif_similarity + CpG_island + H3K27acA2 + H3K27acA1 + H3K4me1A1 + H3K4me1A2 + H3K9acA1 + H3K9acA2 + H3K4me3A1 + H3K4me3A2 + H3K27me3A1 + H3K27me3A2 + H3K9me3A1 + H3K9me3A2,
                        data = traindf, ntree = 500)
  # predict probabilities on the test chromosome
  pred_probs <- predict(model, newdata = testdf, type = "prob")[, "1"]
  roc_obj <- roc(testdf$class, pred_probs)
  auc_values <- c(auc_values, auc(roc_obj))
  
  # save predictions and labels as list for ROC
  all_labels <- c(all_labels, testdf$class)
  all_predictions <- c(all_predictions, pred_probs)
}


# average auc performance across chromosomes
mean_auc <- mean(auc_values)
print(auc_values)

#dev.off()
# plot avg ROC using all predictions 
plot(roc(all_labels, all_predictions))
roc(all_labels, all_predictions)

# plot precision-recall curve averaged
pr <- pr.curve(scores.class0 = all_predictions[all_labels == 1],
               scores.class1 = all_predictions[all_labels == 0],
               curve = TRUE)
plot(pr)

varImpPlot(model)
# importance
imp <- importance(model, type = 2)  # MeanDecreaseGini
imp_df <- data.frame(Feature = rownames(imp), Importance = imp[, 1])

# rename features 
imp_df$Feature <- gsub("^TPM$", "Expression", imp_df$Feature)
imp_df$Feature <- gsub("^log10_distance$", "Distance", imp_df$Feature)
imp_df$Feature <- gsub("^motif_similarity$", "Motif Similarity", imp_df$Feature)
imp_df$Feature <- gsub("^CpG_island$", "CpG binary", imp_df$Feature)
imp_df$Feature <- gsub("A2$", " (E)", imp_df$Feature)
imp_df$Feature <- gsub("A1$", " (P)", imp_df$Feature)

# sort by importance
imp_df <- imp_df[order(imp_df$Importance, decreasing = TRUE), ]

library(ggplot2)
# plot 
ggplot(imp_df, aes(x = Importance, y = reorder(Feature, Importance))) +
  geom_col(fill = "skyblue") +
  labs(title = "Feature Importance Final Model ",
       x = "Mean Decrease Gini",
       y = "Feature") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 12)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)))
#RF plots
library(randomForest)
library(pROC)
library(PRROC)

set.seed(42)
chromosomes <- unique(combined_flc_final.df$seqnames1)

#for plotting define model formulas
#baseline model
f_Base_flc <- as.formula(as.factor(class) ~ log10_distance + TPM)

#model2
f_Model2_flc <- as.formula(as.factor(class) ~ log10_distance + TPM  +
                         H3K27acA1 + H3K27acA2 +
                         H3K4me1A1 + H3K4me1A2 + +H3K9acA1 + H3K9acA2 +
                         H3K4me3A1 + H3K4me3A2 +
                         H3K27me3A1 + H3K27me3A2 +
                         H3K9me3A1 + H3K9me3A2)
#model3
f_Model3_flc <- as.formula(as.factor(class) ~ log10_distance + TPM + motif_similarity + CpG_island)

#model final
f_Final_flc   <- as.formula(as.factor(class) ~ log10_distance + TPM  + H3K27acA1 +H3K27acA2+H3K4me1A1 +H3K4me1A2+ +H3K9acA1+H3K9acA2+H3K4me3A1+H3K4me3A2+H3K27me3A1+H3K27me3A2+H3K9me3A1+H3K9me3A2+motif_similarity+CpG_island)
#final no enhancer mark
f_Final_noH3K27acA2_flc <- as.formula(as.factor(class) ~ log10_distance + TPM +H3K27acA1 +H3K4me1A1 +H3K4me1A2+ +H3K9acA1+H3K9acA2+H3K4me3A1+H3K4me3A2+H3K27me3A1+H3K27me3A2+H3K9me3A1+H3K9me3A2+motif_similarity+CpG_island)

#cross chromosome validation loop
run_cv <- function(formula_obj) {
  all_labels <- c()
  all_predictions <- c()
  auc_values <- c()
  
  for (chr in chromosomes) {
    testdf  <- combined_flc_final.df[combined_flc_final.df$seqnames1 == chr, ]
    traindf <- combined_flc_final.df[combined_flc_final.df$seqnames1 != chr, ]
    
    model <- randomForest(formula_obj, data = traindf, ntree = 500)
    pred_probs <- predict(model, newdata = testdf, type = "prob")[,"1"]
    
    roc_obj <- roc(testdf$class, pred_probs)
    auc_values <- c(auc_values, auc(roc_obj))
    
    all_labels <- c(all_labels, testdf$class)
    all_predictions <- c(all_predictions, pred_probs)
  }
  
  list(labels = all_labels, preds = all_predictions, auc = mean(auc_values))
}

res_Base_flc               <- run_cv(f_Base_flc)
res_Model2_flc             <- run_cv(f_Model2_flc)
res_Model3_flc             <- run_cv(f_Model3_flc)
res_Final_flc              <- run_cv(f_Final_flc)
res_Final_noH3K27acA2_flc  <- run_cv(f_Final_noH3K27acA2_flc)


#create pr obj
pr_Base_flc <- pr.curve(scores.class0 = res_Base_flc$preds[res_Base_flc$labels == 1],
                        scores.class1 = res_Base_flc$preds[res_Base_flc$labels == 0],
                        curve = TRUE)

pr_Model2_flc <- pr.curve(scores.class0 = res_Model2_flc$preds[res_Model2_flc$labels == 1],
                          scores.class1 = res_Model2_flc$preds[res_Model2_flc$labels == 0],
                          curve = TRUE)

pr_Model3_flc <- pr.curve(scores.class0 = res_Model3_flc$preds[res_Model3_flc$labels == 1],
                          scores.class1 = res_Model3_flc$preds[res_Model3_flc$labels == 0],
                          curve = TRUE)

pr_Final_flc <- pr.curve(scores.class0 = res_Final_flc$preds[res_Final_flc$labels == 1],
                         scores.class1 = res_Final_flc$preds[res_Final_flc$labels == 0],
                         curve = TRUE)


pr_Final_noH3K27acA2_flc <- pr.curve(scores.class0 = res_Final_noH3K27acA2_flc$preds[res_Final_noH3K27acA2_flc$labels == 1],
                                     scores.class1 = res_Final_noH3K27acA2_flc$preds[res_Final_noH3K27acA2_flc$labels == 0],
                                     curve = TRUE)

#model 2 plot vs baseline
plot(pr_Model2_flc,            main = "FLC: Model 2", lwd = 2, col = "yellow")
plot(pr_Base_flc,             add = TRUE, lwd = 2, col = "red")
plot.new()
legend("bottomleft",
       legend = c(
         sprintf("Model 2(AUC = %.3f)",           pr_Model2_flc$auc.integral),
         sprintf("Baseline Model (AUC = %.3f)",        pr_Base_flc$auc.integral)
       ),
       col = c("yellow","red"), lwd = 2, bty = "n")


#model 3 plot vs baseline
plot(pr_Model3_flc,            main = "FLC: Model 3", lwd = 2, col = "green")
plot(pr_Base_flc,             add = TRUE, lwd = 2, col = "red")
plot.new()
legend("bottomleft",
       legend = c(
         sprintf("Model 3 (AUC = %.3f)",           pr_Model3_flc$auc.integral),
         sprintf("Baseline (AUC = %.3f)",        pr_Base_flc$auc.integral)
       ),
       col = c("green","red"), lwd = 2, bty = "n")


#all models model vs baseline vs no k27
plot(pr_Final_flc,            main = "FLC: Final Model", lwd = 2, col = "green")
plot(pr_Final_noH3K27acA2_flc,add = TRUE, lwd = 2, col = "blue")
plot(pr_Model2_flc,      add = TRUE, lwd = 2, col = "yellow")
plot(pr_Model3_flc,add = TRUE, lwd = 2, col = "orange")
plot(pr_Base_flc,             add = TRUE, lwd = 2, col = "red")
plot.new()
legend("bottomleft",
       legend = c(
         sprintf("Final Model(AUC = %.3f)",           pr_Final_flc$auc.integral),
         sprintf("Without H3K27acA2 (AUC = %.3f)", pr_Final_noH3K27acA2_flc$auc.integral),
         sprintf("Model 2 (AUC = %.3f)",  pr_Model2_flc$auc.integral),
         sprintf("Model 3 (AUC = %.3f)", pr_Model3_flc$auc.integral),
         sprintf("Baseline (AUC = %.3f)",        pr_Base_flc$auc.integral)
       ),
       col = c("green","blue","yellow","orange","red"), lwd = 2, bty = "n")
