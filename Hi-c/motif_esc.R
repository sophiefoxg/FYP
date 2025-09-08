###motif enrichment esc
#install packages
if (!require("BiocManager", quietly = TRUE))
 install.packages("BiocManager")
BiocManager::install("GenomicInteractions")
BiocManager::install("TFBSTools")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm9")
BiocManager::install("lsa")
BiocManager::install("JASPAR2022") 
#load
library(JASPAR2022)
library(lsa)
library(BSgenome.Mmusculus.UCSC.mm9)
library(TFBSTools)
library(GenomicInteractions)
library(rtracklayer)

setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")
pe_esc.gi<-readRDS("ESC_samples_positive.rds")
negative_gi_esc<-readRDS("negative_interactions_esc.rds")

#load motifs
opts   <- list(species = 10090, matrixtype = "PWM")
motifs <- getMatrixSet(JASPAR2022, opts)

#function for counting motif hits
count_hits_per_motif <- function(motifs, sequence, threshold = "85%") {
  sapply(motifs, function(motif) {
    length(searchSeq(motif, sequence, strand = "*", min.score = threshold))
  })
}

#positives esc
# ESC PROMOTERS
promoters_esc <- anchorOne(pe_esc.gi)
strand(promoters_esc) <- "+"

promoter_vecs_esc <- vector("list", length(promoters_esc))

for (i in seq_along(promoters_esc)) {
  cat("Processing ESC promoter", i, "of", length(promoters_esc), "\n")
  
  # 250 kb upstream and 500 bp downstream from TSS
  promoter <- promoters_esc[i]
  tss <- resize(promoter, width = 1, fix = "start")
  start(tss) <- start(tss) - 2500
  end(tss)   <- end(tss)   + 500
  
  promoter_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm9, tss)
  promoter_vecs_esc[[i]] <- count_hits_per_motif(motifs, promoter_seq[[1]])
}
#save
#saveRDS(promoter_vecs_esc, "esc_promoter_motif_vectors.rds")
#reload
promoter_vecs_esc<-readRDS("esc_promoter_motif_vectors.rds")
head(promoter_vecs_esc)

# ESC ENHANCERS
enhancers_esc <- anchorTwo(pe_esc.gi)

enhancer_vecs_esc <- vector("list", length(enhancers_esc))

for (i in seq_along(enhancers_esc)) {
  cat("Processing ESC enhancer", i, "of", length(enhancers_esc), "\n")
  
  enhancer <- enhancers_esc[i]
  enhancer_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm9, enhancer)
  enhancer_vecs_esc[[i]] <- count_hits_per_motif(motifs, enhancer_seq[[1]])
}
#saveRDS(enhancer_vecs_esc, "esc_enhancer_motif_vectors.rds")
enhancer_vecs_esc<-readRDS("esc_enhancer_motif_vectors.rds")


#cosin sim
cosine_pe_esc <- mapply(cosine, promoter_vecs_esc, enhancer_vecs_esc)
head(cosine_pe_esc)
mean(cosine_pe_esc)
mean(cosine_pe_esc, na.rm = TRUE)
pe_esc.gi <- readRDS("ESC_samples_positive.rds")
mcols(pe_esc.gi)$motif_similarity <- cosine_pe_esc
head(pe_esc.gi)
#saveRDS(pe_esc.gi, "positive_interactions_esc_motif.rds")




# ESC NEGATIVES



negatives_esc <- anchorTwo(negative_gi_esc)
setwd("/mnt/clusters/admiral/data/c2007523/FYP/chip/esc/K27ac/broadpeaks")
#get peaks
h3k27ac <- import("H3K27ac_signal.bed")

#set default length for negative "enhancers" as avg peak length (~270bp)
default_length <- as.integer(round(mean(width(h3k27ac), na.rm = TRUE)))
if (is.na(default_length) || default_length <= 0L) default_length <- 500L  # safety

#set all anchors to default-length (centred)
regions_to_sample <- resize(negatives_esc, width = default_length, fix = "center")

#does negative anchor have a peak?
ov <- findOverlaps(negatives_esc, h3k27ac, ignore.strand = TRUE)

#if it does, replace with the peak coordinates
if (length(ov) > 0) {
  #any overlaps
  qh <- queryHits(ov)
  #the peak information for each overlap
  sh <- subjectHits(ov)
  # peak defined
  hit <- !duplicated(qh)
  idx <- qh[hit]   
  # 
  #define the region to sample
  peaks_for_replace <- h3k27ac[sh[hit]]
  #remove metadata
  mcols(peaks_for_replace) <- NULL
#add coordinates in to regions to sample
  regions_to_sample[idx] <- peaks_for_replace
}
setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")
#10,000

n_to_process <- min(10000L, length(negatives_esc))
negatives_esc_subset <- regions_to_sample[seq_len(n_to_process)]
negative_vecs_esc <- vector("list", n_to_process)

#using regions to sample
for (i in seq_len(n_to_process)) {
  cat("Processing ESC negative", i, "of", n_to_process, "\n")
  
  negative <- negatives_esc_subset[i]                  # peak if overlapped; fallback window otherwise
  negative_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm9, negative)
  negative_vecs_esc[[i]] <- count_hits_per_motif(motifs, negative_seq[[1]])
}
#save
setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")
#saveRDS(negative_vecs_esc, "esc_negative_motif_vectors_first10000.rds")
#reload
negative_vecs_esc_1<-readRDS("esc_negative_motif_vectors_first10000.rds")
head(negative_vecs_esc)


## Chunk 10,001–20,000
start <- 10001L
end <- min(start + 10000L - 1L, length(regions_to_sample))
negatives_subset <- regions_to_sample[start:end]
negative_vecs_esc <- vector("list", length(negatives_subset))

for (i in seq_along(negatives_subset)) {
  cat("Processing ESC negative", start + i - 1L, "of", length(regions_to_sample), "\n")
  negative_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm9, negatives_subset[i])
  negative_vecs_esc[[i]] <- count_hits_per_motif(motifs, negative_seq[[1]])
}
#save
setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")
#saveRDS(negative_vecs_esc, "esc_negative_motif_vectors_rows_10001to20000.rds")
head(negative_vecs_esc)
#reload
negative_vecs_esc_2<-readRDS("esc_negative_motif_vectors_rows_10001to20000.rds")

## Chunk 20,001–30,000
start <- 20001L
end <- min(start + 10000L - 1L, length(regions_to_sample))
negatives_subset <- regions_to_sample[start:end]
negative_vecs_esc <- vector("list", length(negatives_subset))

for (i in seq_along(negatives_subset)) {
  cat("Processing ESC negative", start + i - 1L, "of", length(regions_to_sample), "\n")
  negative_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm9, negatives_subset[i])
  negative_vecs_esc[[i]] <- count_hits_per_motif(motifs, negative_seq[[1]])
}
setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")
#save
#saveRDS(negative_vecs_esc, "esc_negative_motif_vectors_rows20001to30000.rds")
#reload
negative_vecs_esc_3<-readRDS("esc_negative_motif_vectors_rows20001to30000.rds")

## Chunk 30,001–40,000
start <- 30001L
end <- min(start + 10000L - 1L, length(regions_to_sample))
negatives_subset <- regions_to_sample[start:end]
negative_vecs_esc <- vector("list", length(negatives_subset))

for (i in seq_along(negatives_subset)) {
  cat("Processing ESC negative", start + i - 1L, "of", length(regions_to_sample), "\n")
  negative_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm9, negatives_subset[i])
  negative_vecs_esc[[i]] <- count_hits_per_motif(motifs, negative_seq[[1]])
}
setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")
#save
#saveRDS(negative_vecs_esc, "esc_negative_motif_vectors_rows30001to40000.rds")
#reload
negative_vecs_esc_4<-readRDS("esc_negative_motif_vectors_rows30001to40000.rds")


## Chunk 40,001–50,000
start <- 40001L
end <- min(start + 10000L - 1L, length(regions_to_sample))
negatives_subset <- regions_to_sample[start:end]
negative_vecs_esc <- vector("list", length(negatives_subset))

for (i in seq_along(negatives_subset)) {
  cat("Processing ESC negative", start + i - 1L, "of", length(regions_to_sample), "\n")
  negative_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm9, negatives_subset[i])
  negative_vecs_esc[[i]] <- count_hits_per_motif(motifs, negative_seq[[1]])
}
#save
#saveRDS(negative_vecs_esc, "esc_negative_motif_vectors_rows40001to50000.rds")
#reload
negative_vecs_esc_5<-readRDS("esc_negative_motif_vectors_rows_40001to50000.rds")

## Chunk 50,001–60,000
start <- 50001L
end <- min(start + 10000L - 1L, length(regions_to_sample))
negatives_subset <- regions_to_sample[start:end]
negative_vecs_esc <- vector("list", length(negatives_subset))

for (i in seq_along(negatives_subset)) {
  cat("Processing ESC negative", start + i - 1L, "of", length(regions_to_sample), "\n")
  negative_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm9, negatives_subset[i])
  negative_vecs_esc[[i]] <- count_hits_per_motif(motifs, negative_seq[[1]])
}
#save
#saveRDS(negative_vecs_esc, "esc_negative_motif_vectors_rows50001to60000.rds")
#reload
negative_vecs_esc_6<-readRDS("esc_negative_motif_vectors_rows_50001to60000.rds")

setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")

# combine together
# all Anchor2 of negative set "false enhancers"
negative_vec_1<-readRDS("esc_negative_motif_vectors_first10000.rds")
negative_vec_2<-readRDS("esc_negative_motif_vectors_rows_10001to20000.rds")
negative_vec_3<-readRDS("esc_negative_motif_vectors_rows20001to30000.rds")
negative_vec_4<-readRDS("esc_negative_motif_vectors_rows30001to40000.rds")
negative_vec_5<-readRDS("esc_negative_motif_vectors_rows40001to50000.rds")
negative_vec_6<-readRDS("esc_negative_motif_vectors_rows50001to60000.rds")


head(negative_vec_1)
negative_vecs_esc <- c(
  negative_vec_1,
  negative_vec_2,
  negative_vec_3,
  negative_vec_4,
  negative_vec_5,
  negative_vec_6
)

length(negative_vecs_esc)
#save the combined object
#saveRDS(negative_vecs_esc, "esc_negative_motif_vectors.rds")
setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")

#promoter-negatives matching
negative_vecs_esc <- readRDS("esc_negative_motif_vectors.rds")
promoter_vecs_esc <- readRDS("esc_promoter_motif_vectors.rds")

#full df
pe_esc.gi <- readRDS("ESC_samples_positive.rds")
negative_gi_esc<-readRDS("negative_interactions_esc.rds")

#promoter-negatives matching
negative_vecs_esc <- readRDS("esc_negative_motif_vectors.rds")
promoter_vecs_esc <- readRDS("esc_promoter_motif_vectors.rds")

#define promoter ids using positive set ranges
promoter_lookup <- setNames(
  promoter_vecs_esc,
  paste0(seqnames(anchorOne(pe_esc.gi)), ":",
         start(anchorOne(pe_esc.gi)), "-",
         end(anchorOne(pe_esc.gi)))
)


head(promoter_lookup)
head(pe_esc.gi)
#cosne calc
# adding promoter motifs using the promoter IDs
promoter_ids <- paste0(
  seqnames(anchorOne(negative_gi_esc)), ":",
  start(anchorOne(negative_gi_esc)), "-",
  end(anchorOne(negative_gi_esc))
)

# get promoter vectors
promoter_vecs_neg_esc <- promoter_lookup[promoter_ids]
length(promoter_vecs_neg_esc)
library(lsa)
#calculate cosine similarity
cosine_neg_esc <- mapply(cosine, promoter_vecs_neg_esc, negative_vecs_esc)
head(negative_vecs_esc)

head(cosine_neg_esc)
mean(cosine_neg_esc)

mcols(negative_gi_esc)$motif_similarity <- cosine_neg_esc
#saveRDS(negative_gi_esc, "negative_interactions_esc_motif.rds")
hist(negative_gi_esc$motif_similarity)
length(negative_gi_esc)
length(pe_esc.gi)

#head(npe_esc.gi#head(negative_gi_esc)

#cosine similarity:
#promoter-enhancer
promoter_vecs_esc <- readRDS("esc_promoter_motif_vectors.rds")
enhancer_vecs_esc <- readRDS("esc_enhancer_motif_vectors.rds")
cosine_pe_esc <- mapply(cosine, promoter_vecs_esc, enhancer_vecs_esc)
length(cosine_pe_esc)

mcols(pe_esc.gi)$motif_similarity <- cosine_pe_esc
#saveRDS(pe_esc.gi, "positive_interactions_esc_motif.rds")


length(pe_esc.gi)
#summary
pe_esc.gi<-readRDS("positive_interactions_esc_motif.rds")
hist(pe_esc.gi$motif_similarity)


#CPG annotation

#add in cpg annotation
library(rtracklayer)
setwd("/mnt/clusters/admiral/data/c2007523/FYP/")
cpg<-import("CpG_islands.bed", format="BED")
promoters_ESC <- anchorOne(pe_esc.gi)
cpg_hits_ESC <- findOverlaps(promoters_ESC, cpg)

pe_esc.gi$CpG_island <- "noCpG"
pe_esc.gi$CpG_island[unique(queryHits(cpg_hits_ESC))] <- "CpG"
sum(pe_esc.gi$CpG_island=="CpG")
sum(pe_esc.gi$CpG_island=="noCpG")
# annotate negative promoters with CpG information
promoters_ESC_neg <- anchorOne(negative_gi_esc)
cpg_hits_ESC_neg <- findOverlaps(promoters_ESC_neg, cpg)
negative_gi_esc$CpG_island <- "noCpG"
negative_gi_esc$CpG_island[unique(queryHits(cpg_hits_ESC_neg))] <- "CpG"

#Cpg levels
setwd("/mnt/clusters/admiral/data/c2007523/FYP")
CpG_levels<-import("cpg_counts.bed")
ov <- findOverlaps(anchorOne(pe_esc.gi), CpG_levels)
signal_vals <- tapply(as.numeric(mcols(CpG_levels)$name[subjectHits(ov)]),
                      queryHits(ov), mean)
pe_esc.gi$CpG_levels <- 0
pe_esc.gi$CpG_levels[as.numeric(names(signal_vals))] <- signal_vals
mean(pe_esc.gi$CpG_levels)
sum(pe_esc.gi$CpG_levels > 0, na.rm = TRUE)
sum(pe_esc.gi$CpG_levels == 0, na.rm = TRUE)


ov <- findOverlaps(anchorOne(negative_gi_esc), CpG_levels)
signal_vals <- tapply(as.numeric(mcols(CpG_levels)$name[subjectHits(ov)]),
                      queryHits(ov), mean)
negative_gi_esc$CpG_levels <- 0
negative_gi_esc$CpG_levels[as.numeric(names(signal_vals))] <- signal_vals
summary(negative_gi_esc$CpG_levels)
mean(negative_gi_esc$CpG_levels)
sum(negative_gi_esc$CpG_levels > 0, na.rm = TRUE)
sum(negative_gi_esc$CpG_levels== 0, na.rm = TRUE)
####random forests

#cleanup
#add a column that defines if its a negative or positive example:
pe_esc.gi$class<-1
negative_gi_esc$class<-0
#convert to data frame (it didn't like GI objects)
positives.df<-as.data.frame(pe_esc.gi)
negatives.df<-as.data.frame(negative_gi_esc)

head(positives.df)
#removing columns
positives.df <- positives.df[, !names(positives.df) %in% c(
  "strand1", "strand2", "counts", "loe", "anchor2.promoter.id",
  "start1", "end1", "width1",
  "start2", "end2", "width2",
  "anchor1.node.class", "anchor1.promoter.id","seqnames2", "anchor2.node.class", "H3K27ac"
)]

head(positives.df)
head(negatives.df)
negatives.df <- negatives.df[, !names(negatives.df) %in% c(
  "counts", "strand1", "strand2",
  "start1", "end1", "width1",
  "start2", "end2", "width2", "seqnames2",
  "anchor1.name", "anchor1.score", "anchor2.name", "anchor2.score", "anchor1.promoter.id","H3K27ac"
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

head(combined.df)
#save
setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")
#write.table(combined.df, file = "ESC_combined_M3.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
#reload if needed
combined.df <- read.delim("ESC_combined_M3.tsv", sep = "\t", header = TRUE)

combined.df<-na.omit(combined.df)

#####random forests ####
library(randomForest)
library(PRROC)
library(pROC)
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
  #Model 3 all marks, distance, expression, motif similarity+CPG annotations
  model <- randomForest(as.factor(class) ~ log10_distance + TPM  + motif_similarity+CpG_island,
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

dev.off()
#plot avg ROC using all predictions 
plot(roc(all_labels, all_predictions))
roc(all_labels, all_predictions)

# plot precision-recall curve averaged
pr <- pr.curve(scores.class0 = all_predictions[all_labels == 1],
               scores.class1 = all_predictions[all_labels == 0],
               curve = TRUE)
plot(pr)

varImpPlot(model)

library(randomForest)
library(ggplot2)

# importance
imp <- importance(model, type = 2)  # MeanDecreaseGini
imp_df <- data.frame(Feature = rownames(imp), Importance = imp[, 1])

#rename features 
imp_df$Feature <- gsub("^TPM$", "Expression", imp_df$Feature)
imp_df$Feature <- gsub("^log10_distance$", "Distance", imp_df$Feature)
imp_df$Feature <- gsub("A2$", " (E)", imp_df$Feature)
imp_df$Feature <- gsub("A1$", " (P)", imp_df$Feature)

#sort by importance
imp_df <- imp_df[order(imp_df$Importance, decreasing = TRUE), ]

#plot 
ggplot(imp_df, aes(x = Importance, y = reorder(Feature, Importance))) +
  geom_col(fill = "skyblue") +
  labs(title = "Feature Importance Model 3 ",
       x = "Mean Decrease Gini",
       y = "Feature") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 12)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)))


# negatives
negative_motifs_esc <- combined.df[combined.df$class == 0, ]

# positives
positive_motifs_esc <- combined.df[combined.df$class == 1, ]
hist(negative_motifs_esc$motif_similarity)
hist(positive_motifs_esc$motif_similarity)

dev.off()


hist(positive_motifs_esc$motif_similarity,
     main = "ESC Positive set: Motif similarity",
     xlab = "Cosine Similarity",
     col = "lightgreen",
     border = "white",
     probability = TRUE)  # scale to density

hist(negative_motifs_esc$motif_similarity,
     main = "ESC Negative set: Motif Similarity",
     xlab = "Cosine Similarity",
     col = "lightcoral",
     border = "white",
     probability = TRUE)  # scale to density

par(mfrow = c(1,1))  # reset layout

library(ggplot2)

plot_df <- rbind(
  data.frame(motif_similarity = positive_motifs_esc$motif_similarity, class = "Positive"),
  data.frame(motif_similarity = negative_motifs_esc$motif_similarity, class = "Negative")
)

ggplot(plot_df, aes(x = motif_similarity, fill = class)) +
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.7) +
  facet_wrap(~class, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("Positive" = "lightgreen", "Negative" = "lightcoral")) +
  theme_minimal(base_size = 14) +
  labs(title = "Distribution of Motif Cosine Similarity Scores",
       x = "Cosine Similarity Score", y = "Density")
