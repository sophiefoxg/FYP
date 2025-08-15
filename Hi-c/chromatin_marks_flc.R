#M2
#chromatin marks
library(rtracklayer)
library(GenomicRanges)
library(GenomicInteractions)

setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")


pe_flc.gi<-readRDS("FLC_samples_positive.rds")
negative_gi<-readRDS("negative_interactions_flc.rds")

#enhancer marks
setwd("/mnt/clusters/admiral/data/c2007523/FYP/chip/flc/K9ac/broadpeaks")


# adding amount of acetylation: 
h3K9ac <- import("H3K9ac_signal.bed", format="BED")

#anchor one
ov <- findOverlaps(anchorOne(pe_flc.gi), h3K9ac)
signal_vals <- tapply(as.numeric(mcols(h3K9ac)$name[subjectHits(ov)]),
                      queryHits(ov), mean)
pe_flc.gi$H3K9acA1 <- 0
pe_flc.gi$H3K9acA1[as.numeric(names(signal_vals))] <- signal_vals
head(pe_flc.gi)

mean(pe_flc.gi$H3K9acA1)
sum(pe_flc.gi$H3K9acA1 > 0, na.rm = TRUE)



# check if any negatives overlap h3K9ac

#anchor 1
ov_neg <- findOverlaps(anchorOne(negative_gi), h3K9ac)
neg_signal_vals <- tapply(as.numeric(mcols(h3K9ac)$name[subjectHits(ov_neg)]),
                          queryHits(ov_neg), mean)
negative_gi$H3K9acA1 <- 0
negative_gi$H3K9acA1[as.numeric(names(neg_signal_vals))] <- neg_signal_vals
head(negative_gi)
mean(negative_gi$H3K9acA1)
sum(negative_gi$H3K9acA1 > 0, na.rm = TRUE)
##anchor 2
ov <- findOverlaps(anchorTwo(pe_flc.gi), h3K9ac)
signal_vals <- tapply(as.numeric(mcols(h3K9ac)$name[subjectHits(ov)]),
                      queryHits(ov), mean)
pe_flc.gi$H3K9acA2 <- 0
pe_flc.gi$H3K9acA2[as.numeric(names(signal_vals))] <- signal_vals
head(pe_flc.gi)

mean(pe_flc.gi$H3K9acA2)
sum(pe_flc.gi$H3K9acA2 > 0, na.rm = TRUE)
#anchor 2
ov_neg <- findOverlaps(anchorTwo(negative_gi), h3K9ac)
neg_signal_vals <- tapply(as.numeric(mcols(h3K9ac)$name[subjectHits(ov_neg)]),
                          queryHits(ov_neg), mean)
negative_gi$H3K9acA2 <- 0
negative_gi$H3K9acA2[as.numeric(names(neg_signal_vals))] <- neg_signal_vals
head(negative_gi)
mean(negative_gi$H3K9acA2)
sum(negative_gi$H3K9acA2 > 0, na.rm = TRUE)

setwd("/mnt/clusters/admiral/data/c2007523/FYP/chip/flc/K4me3/broadpeaks")

# adding amount of methylation: 
h3k4me3 <- import("H3K4me3_signal.bed", format="BED")

ov <- findOverlaps(anchorOne(pe_flc.gi), h3k4me3)
signal_vals <- tapply(as.numeric(mcols(h3k4me3)$name[subjectHits(ov)]),
                      queryHits(ov), mean)
pe_flc.gi$H3K4me3A1 <- 0
pe_flc.gi$H3K4me3A1[as.numeric(names(signal_vals))] <- signal_vals
head(pe_flc.gi)

mean(pe_flc.gi$H3K4me3A1)
sum(pe_flc.gi$H3K4me3A1 > 0, na.rm = TRUE)

# check if any negatives overlap h3k4me3
ov_neg <- findOverlaps(anchorOne(negative_gi), h3k4me3)
neg_signal_vals <- tapply(as.numeric(mcols(h3k4me3)$name[subjectHits(ov_neg)]),
                          queryHits(ov_neg), mean)
negative_gi$H3K4me3A1 <- 0
negative_gi$H3K4me3A1[as.numeric(names(neg_signal_vals))] <- neg_signal_vals
mean(negative_gi$H3K4me3A1)
sum(negative_gi$H3K4me3A1 > 0, na.rm = TRUE)

#anchor 2
ov <- findOverlaps(anchorTwo(pe_flc.gi), h3k4me3)
signal_vals <- tapply(as.numeric(mcols(h3k4me3)$name[subjectHits(ov)]),
                      queryHits(ov), mean)
pe_flc.gi$H3K4me3A2 <- 0
pe_flc.gi$H3K4me3A2[as.numeric(names(signal_vals))] <- signal_vals
head(pe_flc.gi)

mean(pe_flc.gi$H3K4me3A2)
sum(pe_flc.gi$H3K4me3A2 > 0, na.rm = TRUE)

# check if any negatives overlap h3k4me3
ov_neg <- findOverlaps(anchorTwo(negative_gi), h3k4me3)
neg_signal_vals <- tapply(as.numeric(mcols(h3k4me3)$name[subjectHits(ov_neg)]),
                          queryHits(ov_neg), mean)
negative_gi$H3K4me3A2 <- 0
negative_gi$H3K4me3A2[as.numeric(names(neg_signal_vals))] <- neg_signal_vals
mean(negative_gi$H3K4me3A2)
sum(negative_gi$H3K4me3A2 > 0, na.rm = TRUE)

setwd("/mnt/clusters/admiral/data/c2007523/FYP/chip/flc/K9me3/broadpeaks")

# adding amount of methylation: 
h3k9me3 <- import("H3K9me3_signal.bed", format="BED")
#anchor one
ov <- findOverlaps(anchorOne(pe_flc.gi), h3k9me3)
signal_vals <- tapply(as.numeric(mcols(h3k9me3)$name[subjectHits(ov)]),
                      queryHits(ov), mean)
pe_flc.gi$H3K9me3A1 <- 0
pe_flc.gi$H3K9me3A1[as.numeric(names(signal_vals))] <- signal_vals
head(pe_flc.gi)

mean(pe_flc.gi$H3K9me3A1)
sum(pe_flc.gi$H3K9me3A1 > 0, na.rm = TRUE)

# check if any negatives overlap h3k9me3
ov_neg <- findOverlaps(anchorOne(negative_gi), h3k9me3)
neg_signal_vals <- tapply(as.numeric(mcols(h3k9me3)$name[subjectHits(ov_neg)]),
                          queryHits(ov_neg), mean)
negative_gi$H3K9me3A1 <- 0
negative_gi$H3K9me3A1[as.numeric(names(neg_signal_vals))] <- neg_signal_vals
mean(negative_gi$H3K9me3A1)
sum(negative_gi$H3K9me3A1 > 0, na.rm = TRUE)
#anchor two
ov <- findOverlaps(anchorTwo(pe_flc.gi), h3k9me3)
signal_vals <- tapply(as.numeric(mcols(h3k9me3)$name[subjectHits(ov)]),
                      queryHits(ov), mean)
pe_flc.gi$H3K9me3A2 <- 0
pe_flc.gi$H3K9me3A2[as.numeric(names(signal_vals))] <- signal_vals
head(pe_flc.gi)

mean(pe_flc.gi$H3K9me3A2)
sum(pe_flc.gi$H3K9me3A2 > 0, na.rm = TRUE)

# check if any negatives overlap h3k9me3
ov_neg <- findOverlaps(anchorTwo(negative_gi), h3k9me3)
neg_signal_vals <- tapply(as.numeric(mcols(h3k9me3)$name[subjectHits(ov_neg)]),
                          queryHits(ov_neg), mean)
negative_gi$H3K9me3A2 <- 0
negative_gi$H3K9me3A2[as.numeric(names(neg_signal_vals))] <- neg_signal_vals
mean(negative_gi$H3K9me3A2)
sum(negative_gi$H3K9me3A2 > 0, na.rm = TRUE)


setwd("/mnt/clusters/admiral/data/c2007523/FYP/chip/flc/K27me3/broadpeaks")

#H3K27me3 anchor 1 expectations low because already H3K27ac
h3K27me3 <- import("H3K27me3_signal.bed", format="BED")
ov <- findOverlaps(anchorOne(pe_flc.gi), h3K27me3)
head(ov)
signal_vals <- tapply(as.numeric(mcols(h3K27me3)$name[subjectHits(ov)]),
                      queryHits(ov), mean)
pe_flc.gi$H3K27me3A1 <- 0
pe_flc.gi$H3K27me3A1[as.numeric(names(signal_vals))] <- signal_vals
head(pe_flc.gi)
mean(pe_flc.gi$H3K27me3A1)
sum(pe_flc.gi$H3K27me3A1 > 0, na.rm = TRUE)

head(pe_flc.gi)
head(negative_gi)

# check if any negatives 
ov_neg <- findOverlaps(anchorOne(negative_gi), h3K27me3)
head(ov_neg)
neg_signal_vals <- tapply(as.numeric(mcols(h3K27me3)$name[subjectHits(ov_neg)]),
                          queryHits(ov_neg), mean)
negative_gi$H3K27me3A1 <- 0
negative_gi$H3K27me3A1[as.numeric(names(neg_signal_vals))] <- neg_signal_vals
mean(negative_gi$H3K27me3A1)
sum(negative_gi$H3K27me3A1 > 0, na.rm = TRUE)

#H3K27me3 anchor 2
ov <- findOverlaps(anchorTwo(pe_flc.gi), h3K27me3)
head(ov)
signal_vals <- tapply(as.numeric(mcols(h3K27me3)$name[subjectHits(ov)]),
                      queryHits(ov), mean)
pe_flc.gi$H3K27me3A2 <- 0
pe_flc.gi$H3K27me3A2[as.numeric(names(signal_vals))] <- signal_vals
head(pe_flc.gi)

mean(pe_flc.gi$H3K27me3A2)
sum(pe_flc.gi$H3K27me3A2 > 0, na.rm = TRUE)


# check if any negatives 
ov_neg <- findOverlaps(anchorTwo(negative_gi), h3K27me3)
head(ov_neg)
neg_signal_vals <- tapply(as.numeric(mcols(h3K27me3)$name[subjectHits(ov_neg)]),
                          queryHits(ov_neg), mean)
negative_gi$H3K27me3A2 <- 0
negative_gi$H3K27me3A2[as.numeric(names(neg_signal_vals))] <- neg_signal_vals
mean(negative_gi$H3K27me3A2)
sum(negative_gi$H3K27me3A2 > 0, na.rm = TRUE)


setwd("/mnt/clusters/admiral/data/c2007523/FYP/chip/flc/K4me1/broadpeaks")

# adding amount of methylation: 
h3K4me1 <- import("H3K4me1_signal.bed", format="BED")

#anchor one
ov <- findOverlaps(anchorOne(pe_flc.gi), h3K4me1)
signal_vals <- tapply(as.numeric(mcols(h3K4me1)$name[subjectHits(ov)]),
                      queryHits(ov), mean)
pe_flc.gi$H3K4me1A1 <- 0
pe_flc.gi$H3K4me1A1[as.numeric(names(signal_vals))] <- signal_vals
head(pe_flc.gi)

mean(pe_flc.gi$H3K4me1A1)
sum(pe_flc.gi$H3K4me1A1 > 0, na.rm = TRUE)

# check if any negatives overlap h3K4me1
ov_neg <- findOverlaps(anchorOne(negative_gi), h3K4me1)
neg_signal_vals <- tapply(as.numeric(mcols(h3K4me1)$name[subjectHits(ov_neg)]),
                          queryHits(ov_neg), mean)
negative_gi$H3K4me1A1 <- 0
negative_gi$H3K4me1A1[as.numeric(names(neg_signal_vals))] <- neg_signal_vals
mean(negative_gi$H3K4me1A1)
sum(negative_gi$H3K4me1A1 > 0, na.rm = TRUE)

#anchor two
ov <- findOverlaps(anchorTwo(pe_flc.gi), h3K4me1)
signal_vals <- tapply(as.numeric(mcols(h3K4me1)$name[subjectHits(ov)]),
                      queryHits(ov), mean)
pe_flc.gi$H3K4me1A2 <- 0
pe_flc.gi$H3K4me1A2[as.numeric(names(signal_vals))] <- signal_vals
head(pe_flc.gi)

mean(pe_flc.gi$H3K4me1A2)
sum(pe_flc.gi$H3K4me1A2 > 0, na.rm = TRUE)

# check if any negatives overlap h3K4me1
ov_neg <- findOverlaps(anchorTwo(negative_gi), h3K4me1)
neg_signal_vals <- tapply(as.numeric(mcols(h3K4me1)$name[subjectHits(ov_neg)]),
                          queryHits(ov_neg), mean)
negative_gi$H3K4me1A2 <- 0
negative_gi$H3K4me1A2[as.numeric(names(neg_signal_vals))] <- neg_signal_vals
mean(negative_gi$H3K4me1A2)
sum(negative_gi$H3K4me1A2 > 0, na.rm = TRUE)
sum(negative_gi$H3K4me1A2 == 0, na.rm = TRUE)
length(negative_gi)



# 
setwd("/mnt/clusters/admiral/data/c2007523/FYP/chip/flc/K27ac/broadpeaks")
h3k27ac<-import("H3K27ac_signal.bed")

#anchor one
ov <- findOverlaps(anchorOne(pe_flc.gi), h3k27ac)
signal_vals <- tapply(as.numeric(mcols(h3k27ac)$name[subjectHits(ov)]),
                      queryHits(ov), mean)
pe_flc.gi$H3K27acA1 <- 0
pe_flc.gi$H3K27acA1[as.numeric(names(signal_vals))] <- signal_vals
head(pe_flc.gi)

mean(pe_flc.gi$H3K27acA1)
sum(pe_flc.gi$H3K27acA1 > 0, na.rm = TRUE)

# check if any negatives overlap h3k27ac
ov_neg <- findOverlaps(anchorOne(negative_gi), h3k27ac)
neg_signal_vals <- tapply(as.numeric(mcols(h3k27ac)$name[subjectHits(ov_neg)]),
                          queryHits(ov_neg), mean)
negative_gi$H3K27acA1 <- 0
negative_gi$H3K27acA1[as.numeric(names(neg_signal_vals))] <- neg_signal_vals
mean(negative_gi$H3K27acA1)
sum(negative_gi$H3K27acA1 > 0, na.rm = TRUE)

#anchor two
ov <- findOverlaps(anchorTwo(pe_flc.gi), h3k27ac)
signal_vals <- tapply(as.numeric(mcols(h3k27ac)$name[subjectHits(ov)]),
                      queryHits(ov), mean)
pe_flc.gi$H3K27acA2 <- 0
pe_flc.gi$H3K27acA2[as.numeric(names(signal_vals))] <- signal_vals
head(pe_flc.gi)

mean(pe_flc.gi$H3K27acA2)
sum(pe_flc.gi$H3K27acA2 > 0, na.rm = TRUE)

# check if any negatives overlap h3k27ac
ov_neg <- findOverlaps(anchorTwo(negative_gi), h3k27ac)
neg_signal_vals <- tapply(as.numeric(mcols(h3k27ac)$name[subjectHits(ov_neg)]),
                          queryHits(ov_neg), mean)
negative_gi$H3K27acA2 <- 0
negative_gi$H3K27acA2[as.numeric(names(neg_signal_vals))] <- neg_signal_vals
mean(negative_gi$H3K27acA2)
sum(negative_gi$H3K27acA2 > 0, na.rm = TRUE)
sum(negative_gi$H3K27acA2 == 0, na.rm = TRUE)
length(negative_gi)




#cleanup
#add a column that defines if its a negative or positive example:
pe_flc.gi$class<-1
negative_gi$class<-0
#convert to data frame (it didn't like GI objects)
positives.df<-as.data.frame(pe_flc.gi)
negatives.df<-as.data.frame(negative_gi)

head(positives.df)
#removing columns
positives.df <- positives.df[, !names(positives.df) %in% c(
  "strand1", "strand2", "counts", "loe", "anchor2.promoter.id",
  "start1", "end1", "width1",
  "start2", "end2", "width2",
  "anchor1.node.class", "anchor1.promoter.id","seqnames2", "anchor2.node.class", "CpG_island", "H3K27ac"
)]

head(positives.df)
head(negatives.df)
negatives.df <- negatives.df[, !names(negatives.df) %in% c(
  "counts", "strand1", "strand2",
  "start1", "end1", "width1",
  "start2", "end2", "width2", "seqnames2",
  "anchor1.name", "anchor1.score", "anchor2.name", "anchor2.score", "anchor1.promoter.id", "CpG_island", "H3K27ac"
)]

head(negatives.df)
head(positives.df)

#scaling distance
positives.df$log10_distance<-log10(positives.df$distance)
negatives.df$log10_distance<-log10(negatives.df$distance)
#removing old distnce column
negatives.df<-negatives.df[,!(names(negatives.df) %in% c("distance"))]
positives.df<-positives.df[,!(names(positives.df) %in% c("distance"))]
#only using TPM, distance, chr
combined.df <- rbind(positives.df, negatives.df)
head(combined.df)
#one had 0 expression so just removed it
colSums(is.na(combined.df))
combined.df<-na.omit(combined.df)


#save
setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")
write.table(combined.df, file = "FLC_combined_M2.tsv", sep = "\t", quote = FALSE,row.names = FALSE)
#read in if necessary
combined.df <- read.delim("FLC_combined_M2.tsv", sep = "\t", header = TRUE)

##### random forests ####
library(randomForest)
library(PRROC)
library(pROC)

chromosomes <- unique(combined.df$seqnames1)

set.seed(42)
all_labels <- c()
all_predictions <- c()
auc_values <- c()

for (chr in chromosomes) {
  testdf  <- combined.df[combined.df$seqnames1 == chr, ]
  traindf <- combined.df[combined.df$seqnames1 != chr, ]
  
  model <- randomForest(as.factor(class) ~ log10_distance + TPM +
                          H3K27acA2 + H3K27acA1 +
                          H3K4me1A1 + H3K4me1A2 +
                          H3K9acA1 + H3K9acA2 +
                          H3K4me3A1 + H3K4me3A2 +
                          H3K27me3A1 + H3K27me3A2 +
                          H3K9me3A1 + H3K9me3A2,
                        data = traindf, ntree = 500)
  
  pred_probs <- predict(model, newdata = testdf, type = "prob")[,"1"]
  roc_obj <- roc(testdf$class, pred_probs)
  auc_values <- c(auc_values, auc(roc_obj))
  
  all_labels <- c(all_labels, testdf$class)
  all_predictions <- c(all_predictions, pred_probs)
}

mean_auc <- mean(auc_values)
print(auc_values)

plot(roc(all_labels, all_predictions))

pr <- pr.curve(scores.class0 = all_predictions[all_labels == 1],
               scores.class1 = all_predictions[all_labels == 0],
               curve = TRUE)
plot(pr)

varImpPlot(model)
