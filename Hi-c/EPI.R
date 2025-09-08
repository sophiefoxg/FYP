### install packages and data read in ###

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicInteractions")
library(GenomicInteractions)

x1 = read.delim("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c/ESC_promoter_other_significant_interactions.txt")


head(x1)
#start bait and end bait is the first promoter with a gene associated to it, sometimes several genes
#ignore expression.quartile

#start and end in the second column will be what its interacting with, (the other) 
# raw count is the reads associated with it

##log.observed.expected based off of GTHIc binomial distribution model, exactly as descrived,  large number= more reads than expected
#### ESC ####
### data organisation and cleaning ###
#create a genomic ranges object for promoters
esc1.gr = GRanges(x1$chr.bait,IRanges(x1$start.bait, x1$end.bait))
#assign promoters 
esc1.gr$node.class= "promoter"
#assign ensembl id to promoter
esc1.gr$promoter.id= x1$Ensembl.Gene.ID

#create a genomic ranges object for others
esc2.gr = GRanges(x1$chr, IRanges(x1$start, x1$end))
esc2.gr$node.class= "other"
esc2.gr$promoter.id= NA

#removing all promoters that span several genes


head(esc1.gr)
head(esc2.gr)



#promoter and other dataset
po_esc.gi=GenomicInteractions(esc1.gr,esc2.gr, counts = x1$raw.count, loe=x1$log.observed.expected.)
head(po_esc.gi)

# distance between promoter other interactions, using anchor1 and anchor 2
summary(calculateDistances(po_esc.gi))

#^ look into max distance this is large
hist(log10(calculateDistances(po_esc.gi)), xlab = "Distance log10(Mb)", ylab= "Frequency", main = "Promoter-Other Distance (ESC)" )



#adding in hits from chip-seq data
setwd("/mnt/clusters/admiral/data/c2007523/FYP/chip/esc/broadpeaks/")
if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  BiocManager::install("rtracklayer")
}

library(rtracklayer)
library(GenomicRanges)


#broad peaks
setwd("/mnt/clusters/admiral/data/c2007523/FYP/chip/esc/broadpeaks")
esc_broad<-import("broadpeaksesc.bed")
head(esc_broad)

hits_escbroad<-findOverlaps(anchorTwo(po_esc.gi),esc_broad)

po_esc.gi$anchor2.node.class[queryHits(hits_escbroad)] = "enhancer"
table(po_esc.gi$anchor2.node.class)
# using broad peaks, 8075 enhancers found


# promoter enhancer interactions
pe_esc.gi<-po_esc.gi[po_esc.gi$anchor2.node.class== "enhancer"]
head(pe_esc.gi)
pe_esc.gi = pe_esc.gi[pe_esc.gi$anchor1.promoter.id!="not_promoter"]
table(pe_esc.gi$anchor2.node.class)

pe_esc.gi=pe_esc.gi[nchar(pe_esc.gi$anchor1.promoter.id) <= 18, ]

head(pe_esc.gi)
mcols(pe_esc.gi)$distance <- calculateDistances(pe_esc.gi)
hist(log10(calculateDistances(pe_esc.gi)),xlab = "Distance log10(Mb)", ylab= "Frequency", main = "Promoter-Enhancer Distance (ESC)" )

mean(pe_esc.gi$distance)

# adding amount of acetylation: 
#this bed file now also has the signal value from broad peak file
h3k27ac <- import("H3K27ac_signal.bed", format="BED")
head(h3k27ac)
ov <- findOverlaps(anchorTwo(pe_esc.gi), h3k27ac)

# adding signal value and overlapping as numeric object
signal_vals <- tapply(as.numeric(mcols(h3k27ac)$name[subjectHits(ov)]),
                      queryHits(ov), mean)
#adding "amount of acetylation as a column"
pe_esc.gi$H3K27ac[as.numeric(names(signal_vals))] <- signal_vals
mean(pe_esc.gi$H3K27ac)
#avg 3.2
head(pe_esc.gi)
#add in cpg annotation
setwd("/mnt/clusters/admiral/data/c2007523/FYP/")
cpg<-import("CpG_islands.bed", format="BED")
promoters_ESC <- anchorOne(pe_esc.gi)
cpg_hits_ESC <- findOverlaps(promoters_ESC, cpg)

pe_esc.gi$CpG_island <- "noCpG"
pe_esc.gi$CpG_island[unique(queryHits(cpg_hits_ESC))] <- "CpG"
sum(pe_esc.gi$CpG_island=="CpG")
sum(pe_esc.gi$CpG_island=="noCpG")
head(pe_esc.gi)
#read in transcription data
setwd("/mnt/clusters/admiral/data/c2007523/FYP/RNA/esc/")
TPMesc<-read.delim("TPMesc_table.txt", header=TRUE)
head(TPMesc)
#if TPM is higher than 0 then active promoter
matchRNA<-match(pe_esc.gi$anchor1.promoter.id,TPMesc$Geneid)
pe_esc.gi$TPM<-TPMesc$TPM[matchRNA]

setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")
#saveRDS(pe_esc.gi, file = "ESC_samples_positive.rds")
pe_esc.gi<-readRDS("ESC_samples_positive.rds")


table(pe_esc.gi$anchor2.node.class)
 length(unique(pe_esc.gi$anchor1.promoter.id))
 df <- as.data.frame(anchorTwo(pe_esc.gi))
 length(unique(paste(df$seqnames, df$start, df$end, sep = "_")))
 min(pe_esc.gi$distance)
 max(pe_esc.gi$distance)
 
x2 = read.delim("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c/ESC_promoter_promoter_significant_interactions.txt")

#for ESC promoter-promoter data
head (x2)
#promoter 1
esc3.gr = GRanges(x2$chr, IRanges(x2$start, x2$end))
esc3.gr$node.class= "promoter"
esc3.gr$promoter.id=x2$Ensembl.Gene.ID

#promoter 2
esc4.gr = GRanges(x2$chr.1, IRanges(x2$start.1, x2$end.1))
esc4.gr$node.class= "promoter"
esc4.gr$promoter.id=x2$Ensembl.Gene.ID.1

pp_esc.gi=GenomicInteractions(esc3.gr, esc4.gr, counts=x2$raw.count,loe=x2$log.observed.expected.)
head(pp_esc.gi)
promoters.gr = unique(c(esc1.gr[esc1.gr$node.class=="promoter"],
                        esc3.gr[esc3.gr$node.class=="promoter"],
                        esc4.gr[esc4.gr$node.class=="promoter"]))

head(promoters.gr)
promoters.gr = promoters.gr[promoters.gr$promoter.id!="not_promoter"]
promoters.gr=promoters.gr[nchar(promoters.gr$promoter.id) <= 18, ]



### 0 found at promoters


#distance between promoter promoter interactions df

hist(log10(calculateDistances(pp_esc.gi)),xlab = "Distance log10(Mb)", ylab= "Frequency", main = "Promoter-Promoter Distance (ESC)" )
mcols(pp_esc.gi)$distance <- calculateDistances(pp_esc.gi)



#removing all promoters that are "not_promoters"
pp_esc.gi = pp_esc.gi[pp_esc.gi$anchor1.promoter.id!="not_promoter",]
pp_esc.gi = pp_esc.gi[pp_esc.gi$anchor2.promoter.id!="not_promoter",]

head(pp_esc.gi)
#removing all promoters that span several genes


pp_esc.gi <- pp_esc.gi[nchar(pp_esc.gi$anchor1.promoter.id) <= 18, ]
pp_esc.gi<- pp_esc.gi[nchar(pp_esc.gi$anchor2.promoter.id) <=18,]

head(pp_esc.gi)

#read in tPM data
setwd("/mnt/clusters/admiral/data/c2007523/FYP/RNA/esc/")
TPMesc<-read.delim("TPMesc_table.txt", header=TRUE)
head(TPMesc)
matchRNA1<-match(pp_esc.gi$anchor1.promoter.id,TPMesc$Geneid)
matchRNA2<-match(pp_esc.gi$anchor2.promoter.id,TPMesc$Geneid)
pp_esc.gi$TPM1<-TPMesc$TPM[matchRNA1]
pp_esc.gi$TPM2<-TPMesc$TPM[matchRNA2]


pp_esc.gi$distances<-(calculateDistances(pp_esc.gi))



###making fragments
setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")

fragmentscsv<-read.csv("Fragments.csv")
head(fragmentscsv)
fragments<-import("mm9_HindIII.bed")
head(fragments) 

baits.gr <- GRanges(
  seqnames = fragmentscsv$Chr,
  ranges = IRanges(start = fragmentscsv$Start, end = fragmentscsv$End),
  strand = fragmentscsv$Strand
)


head(baits.gr)
#find overlap between hindIII and baits, taking away the promoters from the "gaps"
hits<-findOverlaps(fragments,baits.gr)
gaps.gr<-fragments[-queryHits(hits)]

head(gaps.gr)
length(gaps.gr)

#posiyive set of promoters
promoter_pos <- anchorOne(pe_esc.gi)


#removing positive set from the fragments
promoter_hits <- findOverlaps(gaps.gr, promoter_pos)
length(promoter_hits)
fragments_no_promoters <- gaps.gr[-queryHits(promoter_hits)]

##no baits or known promoter anchors
gaps_esc.gr<-fragments_no_promoters
length(gaps_esc.gr)


##### sampling negative set
#make it a loop
##### sampling negative set
# Prepare output list
negative_samples_list<-list()

# Loop through all promoters in pe_esc.gi
for (i in seq_along(pe_esc.gi)) {
  # define a promoter
  promoter<-anchorOne(pe_esc.gi[i])
  head(pe_esc.gi)
  # define 1.5mb window with a lower and upper bound 
  # from the promoters start and end site
  window<-1.5e6
  lower<-start(promoter)-window
  upper<-end(promoter)+window

  # create ranges object for Region Of Interest 
  ROI.gr<-GRanges(seqnames=seqnames(promoter), ranges=IRanges(start=lower, end=upper))
  
  # finding this range in gaps.gr (gaps.gr is the hind3bed file without the promoters)
  ROI_in_gaps<-findOverlaps(ROI.gr, gaps_esc.gr)
  ROI.gr<-gaps_esc.gr[subjectHits(ROI_in_gaps)]
  head(ROI.gr)
  
  # checking what enhancers this promoter interacts with if there are multiple
  promoter_hits<-which(overlapsAny(anchorOne(pe_esc.gi), promoter))
  interacting_enhancers<-anchorTwo(pe_esc.gi)[promoter_hits]
  
  # removing any interacting enhancers that overlap the region of interest
  remove_enhancers<-findOverlaps(ROI.gr, interacting_enhancers)
  ROI_no_enhancers<-ROI.gr[-queryHits(remove_enhancers)]
  
  # sample from this region of interest with no enhancers
  if (length(ROI_no_enhancers)>=10) {
    negative_sample<-sample(ROI_no_enhancers, 10, replace=FALSE)
  } else {
    negative_sample<-ROI_no_enhancers
  }
  
  # store full promoter info as metadata for each sampled fragment (had to change because some rows didnt have samples)
  if (length(negative_sample)>0) {
    mcols(negative_sample)$promoter_chr<-as.character(seqnames(promoter))
    mcols(negative_sample)$promoter_start<-start(promoter)
    mcols(negative_sample)$promoter_end<-end(promoter)
    mcols(negative_sample)$promoter_index<-i #adding in what row it came from
    mcols(negative_sample)$promoter_id<-mcols(pe_esc.gi[i])$anchor1.promoter.id
    mcols(negative_sample)$TPM<-mcols(pe_esc.gi[i])$TPM
    
    negative_samples_list[[i]]<-negative_sample
  }
}

negative_interactions<-do.call(c, negative_samples_list)

head(negative_interactions)
length(unique(negative_interactions$promoter_id))
#View(as.data.frame(negative_interactions))
#adding promoters into granges object
promoter_negatives.gr <- GRanges(seqnames=negative_interactions$promoter_chr,
                                 ranges=IRanges(start=negative_interactions$promoter_start,end=negative_interactions$promoter_end)
)

#adding negatives into a granges object with TPM and distance
negatives_a2.gr<-GRanges(seqnames=seqnames(negative_interactions),ranges=IRanges(start=start(negative_interactions), end=end(negative_interactions))
)

negative_gi<-GenomicInteractions(anchor1=promoter_negatives.gr, anchor2=negatives_a2.gr)

mcols(negative_gi)$TPM<-mcols(negative_interactions)$TPM
negative_gi$distance<-calculateDistances(negative_gi)
mcols(negative_gi)$anchor1.promoter.id <- negative_interactions$promoter_id
length(unique(negative_gi$anchor1.promoter.id))
# annotate negative promoters with CpG information
promoters_ESC_neg <- anchorOne(negative_gi)
cpg_hits_ESC_neg <- findOverlaps(promoters_ESC_neg, cpg)
negative_gi$CpG_island <- "noCpG"
negative_gi$CpG_island[unique(queryHits(cpg_hits_ESC_neg))] <- "CpG"

# H3K27ac info
ov_neg <- findOverlaps(anchorTwo(negative_gi), h3k27ac)

neg_signal_vals <- tapply(as.numeric(mcols(h3k27ac)$name[subjectHits(ov_neg)]),
                        queryHits(ov_neg), mean)

negative_gi$H3K27ac <- 0
negative_gi$H3K27ac[as.numeric(names(neg_signal_vals))] <- neg_signal_vals
mean(negative_gi$H3K27ac)
sum(negative_gi$H3K27ac > 0, na.rm = TRUE)
#644
head(negative_gi$H3K27ac[negative_gi$H3K27ac > 0])

mean(negative_gi$distance)
summary(negative_gi$TPM)
sum(negative_gi$CpG_island=="CpG")
sum(negative_gi$CpG_island=="noCpG")

head(negative_gi)
#saveRDS(negative_gi, file= "negative_interactions_no_distance_esc.rds")
negative_gi <- readRDS("negative_interactions_no_distance_esc.rds")

#visualise/plot
hist(log10(calculateDistances(pe_esc.gi)),xlab = "Distance log10(Mb)", ylab= "Frequency", main = "Promoter-Enhancer Distance (ESC)" )
hist(log10(calculateDistances(negative_gi)),xlab = "Distance log10(Mb)", ylab= "Frequency", main = "Promoter-Negative Distance not distance matched(ESC)" )


####random forests

BiocManager::install("randomForest")
library(randomForest)
library(ggplot2)
library(cowplot)

#add a column that defines if its a negative or positive example:
pe_esc.gi$class<-1
negative_gi$class<-0
#convert to data frame (it didn't like GI objects)
positives.df<-as.data.frame(pe_esc.gi)
negatives.df<-as.data.frame(negative_gi)

head(positives.df)
#removing columns
positives.df <- positives.df[, !names(positives.df) %in% c(
  "strand1", "strand2", "counts", "loe", "anchor2.promoter.id",
  "start1", "end1", "width1",
  "start2", "end2", "width2",
  "anchor1.node.class", "anchor1.promoter.id","seqnames2", "anchor2.node.class"
)]

head(positives.df)
head(negatives.df)
negatives.df <- negatives.df[, !names(negatives.df) %in% c(
  "counts", "strand1", "strand2",
  "start1", "end1", "width1",
  "start2", "end2", "width2", "seqnames2",
  "anchor1.name", "anchor1.score"
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
#write.table(combined.df, file = "ESC_samples_no_distance_match.tsv",  sep = "\t", quote = FALSE,  row.names = FALSE)
combined.df <- read.delim("ESC_samples_no_distance_match.tsv", sep = "\t", header = TRUE)


head(combined.df)

###edit for chromosome validation

library(randomForest)
library(pROC)

# Get unique chromosomes from data frame
chromosomes <- unique(combined.df$seqnames1)

set.seed(42)
all_labels <- c()
all_predictions <- c()
auc_values <- c()
#select test chr and train on all others
for (chr in chromosomes) {
#split into train and test based off current chromosome
  testdf  <- combined.df[combined.df$seqnames1 == chr, ]
  traindf <- combined.df[combined.df$seqnames1 != chr, ]
#Model 1 distance and expression 
    model <- randomForest(as.factor(class) ~ log10_distance + TPM, #+ CpG_island + H3K27ac,
                        data = traindf, ntree = 500)
#predict probabilities on test chr  
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
#0.9993

# plot precision-recall curve averaged
library(PRROC)

pr <- pr.curve(scores.class0 = all_predictions[all_labels == 1],
               scores.class1 = all_predictions[all_labels == 0],
               curve = TRUE)
plot(pr)
#results

#distance matched sampling
set.seed(123)
#create an object for all distances
all_distances <- mcols(pe_esc.gi)$distance
#creating empty lists
#nega1 is promoter
#nega2 is fake interaction enhancer/gap
neg_a1_list <- list(); neg_a2_list <- list()
#neg_tpm is RNA and neg id is the promoter id
neg_tpm <- c(); neg_id <- c()
#for this row
for (i in seq_along(pe_esc.gi)) {
#define the promoter
  promoter <- anchorOne(pe_esc.gi[i])
#define the enhancer
  enhancer <- anchorTwo(pe_esc.gi[i])
  #define promoter id
  pid <- mcols(pe_esc.gi[i])$anchor1.promoter.id
  # define tpm
  tpm <- mcols(pe_esc.gi[i])$TPM
  
  tmp_a1 <- list(); tmp_a2 <- list()
#sample 10
  for (j in 1:10) {
    repeat {
#random distance
      rand_dist <- sample(all_distances, 1)
# up or downstream
      sign_dir <- sample(c(-1, 1), 1)
  # the region of interest can be upstream or downstream of the promoter
      ROI.gr <- GRanges(seqnames(promoter),
                        IRanges(if(sign_dir==1) end(promoter)+rand_dist else start(promoter)-rand_dist,
                                width=width(promoter)))
# finds a region that doesnt contain promoter baits       
      ROI.gr <- gaps_esc.gr[subjectHits(findOverlaps(ROI.gr, gaps_esc.gr))]
   #if region exists
      if (length(ROI.gr) > 0) {
        #add sample
        neg <- sample(ROI.gr, 1)
  #does this promoter, interact with this enhancer from positive set
        overlap <- length(findOverlaps(promoter, anchorOne(pe_esc.gi[i]))) > 0 &&
          length(findOverlaps(neg, enhancer)) > 0
#if it's not in positive set, store it as a negative sample
                if (!overlap) { tmp_a1[[j]] <- promoter; tmp_a2[[j]] <- neg; break }
      }
    }
  #if 10 samples dont exist, stop looking
    if (length(tmp_a1) < j) break
  }
#store information for these interactions  
  if (length(tmp_a1) > 0) {
    neg_a1_list[[i]] <- do.call(c, tmp_a1)
    neg_a2_list[[i]] <- do.call(c, tmp_a2)
    neg_tpm <- c(neg_tpm, rep(tpm, length(tmp_a1)))
    neg_id <- c(neg_id, rep(pid, length(tmp_a1)))
  }
#progress report  
  cat(round(i / length(pe_esc.gi) * 100, 2), "% complete\n")
}
# put it into an interactions object, a1 is promoter, a2 is other region
negative_gi <- GenomicInteractions(anchor1 = do.call(c, neg_a1_list),anchor2 = do.call(c, neg_a2_list))
  
 #add TPM and promoter id back in                                
mcols(negative_gi)$TPM <- neg_tpm
mcols(negative_gi)$anchor1.promoter.id <- neg_id

View(as.data.frame(negative_gi))
View(as.data.frame(pe_esc.gi))



head(negative_gi)

head(negative_gi)
head(pe_esc.gi)
#double check that these interactions dont appear in my positive set
findOverlaps(anchorTwo(negative_gi), anchorTwo(pe_esc.gi))
pos_hits <- findOverlaps(
  GenomicInteractions(anchor1 = anchorOne(negative_gi), anchor2 = anchorTwo(negative_gi)),
  pe_esc.gi,
  type = "equal"
)

length(pos_hits)
#calculae distance
mcols(negative_gi)$distance <- calculateDistances(negative_gi)

# annotate negative promoters with CpG information
promoters_ESC_neg <- anchorOne(negative_gi)
cpg_hits_ESC_neg <- findOverlaps(promoters_ESC_neg, cpg)
negative_gi$CpG_island <- "noCpG"
negative_gi$CpG_island[unique(queryHits(cpg_hits_ESC_neg))] <- "CpG"

# check if any negatives overlap H3K27ac
ov_neg <- findOverlaps(anchorTwo(negative_gi), h3k27ac)
neg_signal_vals <- tapply(as.numeric(mcols(h3k27ac)$name[subjectHits(ov_neg)]),
                          queryHits(ov_neg), mean)
negative_gi$H3K27ac <- 0
negative_gi$H3K27ac[as.numeric(names(neg_signal_vals))] <- neg_signal_vals

# quick QC
length(pe_esc.gi)
length(negative_gi)
mean(negative_gi$H3K27ac)
sum(negative_gi$H3K27ac>0)
sum(pe_esc.gi$H3K27ac>0)
mean(negative_gi$distance)
summary(negative_gi$TPM)
sum(negative_gi$CpG_island=="CpG")
sum(negative_gi$CpG_island=="noCpG")
head(negative_gi)
length(unique(negative_gi$anchor1.promoter.id))

# visualize
hist(log10(calculateDistances(pe_esc.gi)), xlab = "Distance log10(Mb)",
     ylab = "Frequency", main = "Promoter-Enhancer Distance (ESC)" )
hist(log10(calculateDistances(negative_gi)), xlab = "Distance log10(Mb)",
     ylab = "Frequency", main = "Promoter-Negative Distance matched (ESC)" )

#setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")
#saveRDS(negative_gi, file="negative_interactions_esc.rds")
negative_gi <- readRDS("negative_interactions_esc.rds")

#cleanup

#add a column that defines if its a negative or positive example:
pe_esc.gi$class<-1
negative_gi$class<-0
#convert to data frame (it didn't like GI objects)
positives.df<-as.data.frame(pe_esc.gi)
negatives.df<-as.data.frame(negative_gi)
#removing columns
head(positives.df)
positives.df <- positives.df[, !names(positives.df) %in% c(
  "strand1", "strand2", "counts", "loe", "anchor2.promoter.id",
  "start1", "end1", "width1",
  "start2", "end2", "width2", "seqnames2",
  "anchor1.node.class", "anchor1.promoter.id", "anchor2.node.class"
)]

head(negatives.df)
negatives.df <- negatives.df[, !names(negatives.df) %in% c(
  "counts", "strand1", "strand2",
  "start1", "end1", "width1",
  "start2", "end2", "width2", "seqnames2","anchor2.name", "anchor2.score",
  "anchor1.promoter.id")]

head(negatives.df)
head(positives.df)

#scaling distance
positives.df$log10_distance<-log10(positives.df$distance)
negatives.df$log10_distance<-log10(negatives.df$distance)
#removing old distnce column
negatives.df<-negatives.df[,!(names(negatives.df) %in% c("distance"))]
positives.df<-positives.df[,!(names(positives.df) %in% c("distance"))]
#only using TPM, distance, 
combined.df <- rbind(positives.df, negatives.df)
head(combined.df)
#one had 0 expression so just removed it
colSums(is.na(combined.df))
combined.df<-na.omit(combined.df)

#save
setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")
#write.table(combined.df,file = "ESC_combined_M1.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

combined.df <- read.delim("ESC_combined_M1.tsv", sep = "\t", header = TRUE)
length(pe_esc.gi)

###edit for chromosome validation

library(randomForest)
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
  #Model 1 only TPM and distance 
  model <- randomForest(as.factor(class) ~ log10_distance + TPM, #+ CpG_island +H3K27ac ,
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
library(PRROC)

pr <- pr.curve(scores.class0 = all_predictions[all_labels == 1],
               scores.class1 = all_predictions[all_labels == 0],
               curve = TRUE)
plot(pr)



#what feature is most important
importance(model)
#amount of acetylation is most important
#then distance
#then TPM
#then cpg

#evaluate model
# labels
true_labels <- as.numeric(as.character(all_labels))
predicted_classes <- ifelse(all_predictions > 0.5, 1, 0)  

#add predictions and labels tocombined.df
results_df <- combined.df
stopifnot(nrow(results_df) == length(true_labels)) 
results_df$true_label <- true_labels
results_df$predicted_class <- predicted_classes
results_df$predicted_prob <- all_predictions

#classify by true label and pred class
results_df$confusion_type <- with(results_df, ifelse(true_label == 1 & predicted_class == 1, "TP",
                                                     ifelse(true_label == 0 & predicted_class == 1, "FP",
                                                            ifelse(true_label == 1 & predicted_class == 0, "FN", "TN"))))

#subset into TP, FP, FN, TN
true_positives   <- results_df[results_df$confusion_type == "TP", ]
false_positives  <- results_df[results_df$confusion_type == "FP", ]
false_negatives  <- results_df[results_df$confusion_type == "FN", ]
true_negatives   <- results_df[results_df$confusion_type == "TN", ]

#false negatives
nrow(false_negatives)
mean(false_negatives$TPM)
mean(false_negatives$log10_distance)
mean(false_negatives$H3K27ac)
sum(false_negatives$CpG_island=="CpG")
sum(false_negatives$CpG_island=="noCpG")
# false positives
nrow(false_positives)
mean(false_positives$TPM)
mean(false_positives$log10_distance)
mean(false_positives$H3K27ac)
sum(false_positives$CpG_island == "CpG")
sum(false_positives$CpG_island == "noCpG")
# true positives
nrow(true_positives)
mean(true_positives$TPM)
mean(true_positives$log10_distance)
mean(true_positives$H3K27ac)
sum(true_positives$CpG_island == "CpG")
sum(true_positives$CpG_island == "noCpG")
#true negatives
nrow(true_negatives)
mean(true_negatives$TPM)
mean(true_negatives$log10_distance)
mean(true_negatives$H3K27ac)
sum(true_negatives$CpG_island == "CpG")
sum(true_negatives$CpG_island == "noCpG")



















#####FLC

### install packages and data read in ###

x1 = read.delim("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c/FLC_promoter_other_significant_interactions.txt")

head(x1)
#start bait and end bait is the first promoter with a gene associated to it, sometimes several genes
#ignore expression.quartile

#start and end in the second column will be what its interacting with, (the other) 
# raw count is the reads associated with it

##log.observed.expected based off of GTHIc binomial distribution model, exactly as descrived,  large number= more reads than expected
#### FLC ####
### data organisation and cleaning ###
#create a genomic ranges object for promoters
flc1.gr = GRanges(x1$chr.bait,IRanges(x1$start.bait, x1$end.bait))
#assign promoters 
flc1.gr$node.class= "promoter"
#assign ensembl id to promoter
flc1.gr$promoter.id= x1$Ensembl.Gene.ID

#create a genomic ranges object for others
flc2.gr = GRanges(x1$chr, IRanges(x1$start, x1$end))
flc2.gr$node.class= "other"
flc2.gr$promoter.id= NA

head(flc1.gr)
head(flc2.gr)

#promoter and other dataset
po_flc.gi=GenomicInteractions(flc1.gr,flc2.gr, counts = x1$raw.count, loe=x1$log.observed.expected.)
head(po_flc.gi)

# distance between promoter other interactions, using anchor1 and anchor 2
summary(calculateDistances(po_flc.gi))
hist(log10(calculateDistances(po_flc.gi)), xlab = "Distance log10(Mb)", ylab= "Frequency", main = "Promoter-Other Distance (FLC)" )

#define enhancers
library(rtracklayer)
#broad peaks
setwd("/mnt/clusters/admiral/data/c2007523/FYP/chip/flc/broadpeaks/")
flc_broad<-import("flcbroad.bed")
hits_flcbroad<-findOverlaps(anchorTwo(po_flc.gi),flc_broad)
po_flc.gi$anchor2.node.class[queryHits(hits_flcbroad)] = "enhancer"
table(po_flc.gi$anchor2.node.class)

# promoter enhancer interactions
pe_flc.gi<-po_flc.gi[po_flc.gi$anchor2.node.class== "enhancer"]
pe_flc.gi = pe_flc.gi[pe_flc.gi$anchor1.promoter.id!="not_promoter"]
pe_flc.gi <- pe_flc.gi[nchar(pe_flc.gi$anchor1.promoter.id) <= 18, ]

length(pe_flc.gi)
mcols(pe_flc.gi)$distance <- calculateDistances(pe_flc.gi)
hist(log10(calculateDistances(pe_flc.gi)),xlab = "Distance log10(Mb)", ylab= "Frequency", main = "Promoter-Enhancer Distance (FLC)" )

# adding amount of acetylation: 
h3k27ac <- import("H3K27ac_signal.bed", format="BED")
ov <- findOverlaps(anchorTwo(pe_flc.gi), h3k27ac)
signal_vals <- tapply(as.numeric(mcols(h3k27ac)$name[subjectHits(ov)]),
                     queryHits(ov), mean)
pe_flc.gi$H3K27ac[as.numeric(names(signal_vals))] <- signal_vals
head(pe_flc.gi)
# add CpG annotation ###remove?
setwd("/mnt/clusters/admiral/data/c2007523/FYP/")
cpg<-import("CpG_islands.bed", format="BED")
promoters_FLC <- anchorOne(pe_flc.gi)
cpg_hits_FLC <- findOverlaps(promoters_FLC, cpg)
pe_flc.gi$CpG_island <- "noCpG"
pe_flc.gi$CpG_island[unique(queryHits(cpg_hits_FLC))] <- "CpG"
sum(pe_flc.gi$CpG_island=="noCpG")

#read in transcription data
setwd("/mnt/clusters/admiral/data/c2007523/FYP/RNA/flc/")
TPMflc<-read.delim("TPMflc_table.txt", header=TRUE)
matchRNA<-match(pe_flc.gi$anchor1.promoter.id,TPMflc$Geneid)
pe_flc.gi$TPM<-TPMflc$TPM[matchRNA]
length(pe_flc.gi)

#save

setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")
#saveRDS(pe_flc.gi, file="FLC_samples_positive.rds")
pe_flc.gi <- readRDS("FLC_samples_positive.rds")

min(pe_flc.gi$distance)
max(pe_flc.gi$distance)

# Read in FLC promoter-promoter interactions
x2_flc <- read.delim("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c/FLC_promoter_promoter_significant_interactions.txt")

# promoter 1
flc3.gr <- GRanges(x2_flc$chr, IRanges(x2_flc$start, x2_flc$end))
flc3.gr$node.class <- "promoter"
flc3.gr$promoter.id <- x2_flc$Ensembl.Gene.ID

# promoter 2
flc4.gr <- GRanges(x2_flc$chr.1, IRanges(x2_flc$start.1, x2_flc$end.1))
flc4.gr$node.class <- "promoter"
flc4.gr$promoter.id <- x2_flc$Ensembl.Gene.ID.1

# Create promoter-promoter GenomicInteractions object
pp_flc.gi <- GenomicInteractions(flc3.gr, flc4.gr,
                                 counts = x2_flc$raw.count,
                                 loe = x2_flc$log.observed.expected.)

# Keep unique promoter ranges
promoters_flc.gr <- unique(c(flc3.gr[flc3.gr$node.class == "promoter"],
                             flc4.gr[flc4.gr$node.class == "promoter"]))
promoters_flc.gr <- promoters_flc.gr[promoters_flc.gr$promoter.id != "not_promoter"]
promoters_flc.gr <- promoters_flc.gr[nchar(promoters_flc.gr$promoter.id) <= 18, ]

# Distance between promoter-promoter interactions
hist(log10(calculateDistances(pp_flc.gi)),
     xlab = "Distance log10(Mb)", ylab = "Frequency",
     main = "Promoter-Promoter Distance (FLC)")
mcols(pp_flc.gi)$distance <- calculateDistances(pp_flc.gi)

# Remove invalid promoters
pp_flc.gi <- pp_flc.gi[pp_flc.gi$anchor1.promoter.id != "not_promoter", ]
pp_flc.gi <- pp_flc.gi[pp_flc.gi$anchor2.promoter.id != "not_promoter", ]
pp_flc.gi <- pp_flc.gi[nchar(pp_flc.gi$anchor1.promoter.id) <= 18, ]
pp_flc.gi <- pp_flc.gi[nchar(pp_flc.gi$anchor2.promoter.id) <= 18, ]

# Read in FLC TPM data
setwd("/mnt/clusters/admiral/data/c2007523/FYP/RNA/flc/")
TPMflc <- read.delim("TPMflc_table.txt", header = TRUE)

# Match TPM to anchors
matchRNA1 <- match(pp_flc.gi$anchor1.promoter.id, TPMflc$Geneid)
matchRNA2 <- match(pp_flc.gi$anchor2.promoter.id, TPMflc$Geneid)
pp_flc.gi$TPM1 <- TPMflc$TPM[matchRNA1]
pp_flc.gi$TPM2 <- TPMflc$TPM[matchRNA2]

# Store distances as numeric
pp_flc.gi$distances <- calculateDistances(pp_flc.gi)



### making fragments ###
setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")
fragments<-import("mm9_HindIII.bed")
fragmentscsv<-read.csv("Fragments.csv")
baits.gr <- GRanges(seqnames = fragmentscsv$Chr,
                    ranges = IRanges(start = fragmentscsv$Start, end = fragmentscsv$End),
                    strand = fragmentscsv$Strand)

hits<-findOverlaps(fragments,baits.gr)
gaps.gr<-fragments[-queryHits(hits)]

promoter_pos <- anchorOne(pe_flc.gi)
promoter_hits <- findOverlaps(gaps.gr, promoter_pos)
fragments_no_promoters <- gaps.gr[-queryHits(promoter_hits)]
gaps_flc.gr<-fragments_no_promoters

##### sampling negative set #####
negative_samples_list<-list()
window<-1.5e6

for (i in seq_along(pe_flc.gi)) {
  promoter<-anchorOne(pe_flc.gi[i])
  lower<-start(promoter)-window
  upper<-end(promoter)+window
  ROI.gr<-GRanges(seqnames=seqnames(promoter), ranges=IRanges(start=lower, end=upper))
  ROI_in_gaps<-findOverlaps(ROI.gr, gaps_flc.gr)
  ROI.gr<-gaps_flc.gr[subjectHits(ROI_in_gaps)]
  
  promoter_hits<-which(overlapsAny(anchorOne(pe_flc.gi), promoter))
  interacting_enhancers<-anchorTwo(pe_flc.gi)[promoter_hits]
  remove_enhancers<-findOverlaps(ROI.gr, interacting_enhancers)
  ROI_no_enhancers<-ROI.gr[-queryHits(remove_enhancers)]
  
  if (length(ROI_no_enhancers)>=10) {
    negative_sample<-sample(ROI_no_enhancers, 10, replace=FALSE)
  } else {
    negative_sample<-ROI_no_enhancers
  }
  
  if (length(negative_sample)>0) {
    mcols(negative_sample)$promoter_chr<-as.character(seqnames(promoter))
    mcols(negative_sample)$promoter_start<-start(promoter)
    mcols(negative_sample)$promoter_end<-end(promoter)
    mcols(negative_sample)$promoter_index<-i
    mcols(negative_sample)$promoter_id<-mcols(pe_flc.gi[i])$anchor1.promoter.id
    mcols(negative_sample)$TPM<-mcols(pe_flc.gi[i])$TPM
    negative_samples_list[[i]]<-negative_sample
  }
}

negative_interactions<-do.call(c, negative_samples_list)

promoter_negatives.gr <- GRanges(seqnames=negative_interactions$promoter_chr,
                                 ranges=IRanges(start=negative_interactions$promoter_start,end=negative_interactions$promoter_end))
negatives_a2.gr<-GRanges(seqnames=seqnames(negative_interactions),
                         ranges=IRanges(start=start(negative_interactions), end=end(negative_interactions)))

negative_gi<-GenomicInteractions(anchor1=promoter_negatives.gr, anchor2=negatives_a2.gr)
mcols(negative_gi)$TPM<-mcols(negative_interactions)$TPM
negative_gi$distance<-calculateDistances(negative_gi)
mean(negative_gi$distance)
summary(negative_gi$TPM)

# annotate negatives with CpG ####remove?
promoters_FLC_neg <- anchorOne(negative_gi)
cpg_hits_FLC_neg <- findOverlaps(promoters_FLC_neg, cpg)
negative_gi$CpG_island <- "noCpG"
negative_gi$CpG_island[unique(queryHits(cpg_hits_FLC_neg))] <- "CpG"
sum(negative_gi$CpG_island=="noCpG")
sum(negative_gi$CpG_island=="CpG")

# add H3K27ac to negatives ###remove?
ov_neg <- findOverlaps(anchorTwo(negative_gi), h3k27ac)
neg_signal_vals <- tapply(as.numeric(mcols(h3k27ac)$name[subjectHits(ov_neg)]),
                        queryHits(ov_neg), mean)
negative_gi$H3K27ac <- 0
negative_gi$H3K27ac[as.numeric(names(neg_signal_vals))] <- neg_signal_vals

mean(negative_gi$H3K27ac)
mean(pe_flc.gi$H3K27ac)
sum(negative_gi$H3K27ac>0)

mcols(negative_gi)$anchor1.promoter.id <- negative_interactions$promoter_id

length(unique(negative_gi$anchor1.promoter.id))
head(negative_gi)


#save
setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")
#saveRDS(negative_gi, file = "negative_interactions_no_distance_flc.rds")
negative_gi <- readRDS("negative_interactions_no_distance_flc.rds")


#### cleanup

pe_flc.gi$class<-1
negative_gi$class<-0
positives.df<-as.data.frame(pe_flc.gi)
negatives.df<-as.data.frame(negative_gi)

positives.df <- positives.df[, !names(positives.df) %in% c(
  "strand1", "strand2", "counts", "loe", "anchor2.promoter.id",
  "start1", "end1", "width1",
  "start2", "end2", "width2", "seqnames2",
  "anchor1.node.class", "anchor1.promoter.id", "anchor2.node.class"
)]
negatives.df <- negatives.df[, !names(negatives.df) %in% c(
  "counts", "strand1", "strand2",
  "start1", "end1", "width1",
  "start2", "end2", "width2", "seqnames2","anchor1.promoter.id",
  "promoter_id")]

positives.df$log10_distance<-log10(positives.df$distance)
negatives.df$log10_distance<-log10(negatives.df$distance)
negatives.df<-negatives.df[,!(names(negatives.df) %in% c("distance"))]
positives.df<-positives.df[,!(names(positives.df) %in% c("distance"))]
head(positives.df)
head(negatives.df)
combined.df <- rbind(positives.df, negatives.df)
combined.df<-na.omit(combined.df)

#save
#write.table(combined.df,  file = "FLC_samples_no_distance_match.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
combined.df <- read.delim("FLC_samples_no_distance_match.tsv", sep = "\t", header = TRUE)

### cross chromosome validation with average ROC/PR ###
library(randomForest)
library(pROC)
library(PRROC)

chromosomes <- unique(combined.df$seqnames1)
set.seed(42)
all_labels <- c()
all_predictions <- c()
auc_values <- c()

for (chr in chromosomes) {
  testdf  <- combined.df[combined.df$seqnames1 == chr, ]
  traindf <- combined.df[combined.df$seqnames1 != chr, ]
  #model 1 only distance and expression
  model <- randomForest(as.factor(class) ~ log10_distance + TPM ,# +H3K27ac ,
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
roc(all_labels, all_predictions)

pr <- pr.curve(scores.class0 = all_predictions[all_labels == 1],
               scores.class1 = all_predictions[all_labels == 0],
               curve = TRUE)
plot(pr)



#### distance matching loop (FLC)

#distance match
set.seed(123)

all_distances <- mcols(pe_flc.gi)$distance
neg_a1_list <- list(); neg_a2_list <- list()
neg_tpm <- c(); neg_id <- c()

for (i in seq_along(pe_flc.gi)) {
  promoter <- anchorOne(pe_flc.gi[i])
  enhancer <- anchorTwo(pe_flc.gi[i])
  pid <- mcols(pe_flc.gi[i])$anchor1.promoter.id
  tpm <- mcols(pe_flc.gi[i])$TPM
  
  tmp_a1 <- list(); tmp_a2 <- list()
  
  for (j in 1:10) {
    repeat {
      rand_dist <- sample(all_distances, 1)
      sign_dir <- sample(c(-1, 1), 1)
      ROI.gr <- GRanges(seqnames(promoter),
                        IRanges(if(sign_dir==1) end(promoter)+rand_dist else start(promoter)-rand_dist,
                                width=width(promoter)))
      
      ROI.gr <- gaps_flc.gr[subjectHits(findOverlaps(ROI.gr, gaps_flc.gr))]
      if (length(ROI.gr) > 0) {
        neg <- sample(ROI.gr, 1)
        overlap <- length(findOverlaps(promoter, anchorOne(pe_flc.gi[i]))) > 0 &&
          length(findOverlaps(neg, enhancer)) > 0
        if (!overlap) { tmp_a1[[j]] <- promoter; tmp_a2[[j]] <- neg; break }
      }
    }
    if (length(tmp_a1) < j) break
  }
  
  if (length(tmp_a1) > 0) {
    neg_a1_list[[i]] <- do.call(c, tmp_a1)
    neg_a2_list[[i]] <- do.call(c, tmp_a2)
    neg_tpm <- c(neg_tpm, rep(tpm, length(tmp_a1)))
    neg_id <- c(neg_id, rep(pid, length(tmp_a1)))
  }
  
  cat(round(i / length(pe_flc.gi) * 100, 2), "% complete\n")
}

negative_gi <- GenomicInteractions(anchor1 = do.call(c, neg_a1_list),
                                   anchor2 = do.call(c, neg_a2_list))
mcols(negative_gi)$TPM <- neg_tpm
mcols(negative_gi)$anchor1.promoter.id <- neg_id

View(as.data.frame(negative_gi))

head(pe_flc.gi)
head(negative_gi)


# add TPM and distance back in

findOverlaps(anchorTwo(negative_gi), anchorTwo(pe_flc.gi))
pos_hits <- findOverlaps(
  GenomicInteractions(anchor1 = anchorOne(negative_gi), anchor2 = anchorTwo(negative_gi)),
  pe_flc.gi,
  type = "equal"
)

length(pos_hits)

mcols(negative_gi)$distance <- calculateDistances(negative_gi)

# annotate negative promoters with CpG information 
promoters_FLC_neg <- anchorOne(negative_gi)
cpg_hits_FLC_neg <- findOverlaps(promoters_FLC_neg, cpg)
negative_gi$CpG_island <- "noCpG"
negative_gi$CpG_island[unique(queryHits(cpg_hits_FLC_neg))] <- "CpG"

# check if any negatives overlap H3K27ac
ov_neg <- findOverlaps(anchorTwo(negative_gi), h3k27ac)
neg_signal_vals <- tapply(as.numeric(mcols(h3k27ac)$name[subjectHits(ov_neg)]),
                          queryHits(ov_neg), mean)
negative_gi$H3K27ac <- 0
negative_gi$H3K27ac[as.numeric(names(neg_signal_vals))] <- neg_signal_vals

# quick QC
head(negative_gi)
nrow(pe_flc.gi)
length(negative_gi)
sum(negative_gi$H3K27ac>0)
mean(negative_gi$distance)
summary(negative_gi$TPM)

sum(negative_gi$CpG_island=="CpG")
sum(negative_gi$CpG_island=="noCpG")


length(unique(negative_gi$anchor1.promoter.id))
nrow(pe_flc.gi)
mean(pe_flc.gi$H3K27ac)
sum(pe_flc.gi$H3K27ac>0)
mean(pe_flc.gi$distance)
summary(pe_flc.gi$TPM)
sum(pe_flc.gi$CpG_island=="CpG")
sum(pe_flc.gi$CpG_island=="noCpG")

# visualize
hist(log10(calculateDistances(pe_flc.gi)), xlab = "Distance log10(Mb)",
     ylab = "Frequency", main = "Promoter-Enhancer Distance (FLC)" )
hist(log10(calculateDistances(negative_gi)), xlab = "Distance log10(Mb)",
     ylab = "Frequency", main = "Promoter-Negative Distance matched (FLC)" )

setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")
#saveRDS(negative_gi, file= "negative_interactions_flc.rds")
negative_gi <- readRDS("negative_interactions_flc.rds")


####clean up

#add a column that defines if its a negative or positive example:
pe_flc.gi$class<-1
negative_gi$class<-0
#convert to data frame (it didn't like GI objects)
positives.df<-as.data.frame(pe_flc.gi)
negatives.df<-as.data.frame(negative_gi)
#removing columns
head(positives.df)
positives.df <- positives.df[, !names(positives.df) %in% c(
  "strand1", "strand2", "counts", "loe", "anchor2.promoter.id",
  "start1", "end1", "width1",
  "start2", "end2", "width2", "seqnames2",
  "anchor1.node.class", "anchor1.promoter.id", "anchor2.node.class"
)]

head(negatives.df)
negatives.df <- negatives.df[, !names(negatives.df) %in% c(
  "counts", "strand1", "strand2",
  "start1", "end1", "width1",
  "start2", "end2", "width2", "seqnames2","anchor2.name", "anchor2.score",
  "anchor1.promoter.id")]

head(negatives.df)
head(positives.df)
#scaling distance
positives.df$log10_distance<-log10(positives.df$distance)
negatives.df$log10_distance<-log10(negatives.df$distance)
#removing old distnce column
negatives.df<-negatives.df[,!(names(negatives.df) %in% c("distance"))]
positives.df<-positives.df[,!(names(positives.df) %in% c("distance"))]
#only using TPM, distance, 
combined.df <- rbind(positives.df, negatives.df)
head(combined.df)
#one had 0 expression so just removed it
colSums(is.na(combined.df))
combined.df<-na.omit(combined.df)

head(combined.df)
#save
setwd("/mnt/clusters/admiral/data/c2007523/FYP/Hi-c")
#write.table(combined.df,file = "FLC_combined_M1.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
combined.df <- read.delim("FLC_combined_M1.tsv", sep = "\t", header = TRUE)

###edit for chromosome validation

### cross chromosome validation with average ROC/PR ###

chromosomes <- unique(combined.df$seqnames1)
set.seed(42)
all_labels <- c()
all_predictions <- c()
auc_values <- c()

for (chr in chromosomes) {
  testdf  <- combined.df[combined.df$seqnames1 == chr, ]
  traindf <- combined.df[combined.df$seqnames1 != chr, ]
  # model 1 only distance and expression
  model <- randomForest(as.factor(class) ~ log10_distance + TPM  ,
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
roc(all_labels, all_predictions)

pr <- pr.curve(scores.class0 = all_predictions[all_labels == 1],
               scores.class1 = all_predictions[all_labels == 0],
               curve = TRUE)
plot(pr)

#
#what feature is most important
importance(model)
#amount of acetylation is most important
#then distance
#then TPM
#then cpg

#evaluate model
# labels
true_labels <- as.numeric(as.character(all_labels))
predicted_classes <- ifelse(all_predictions > 0.5, 1, 0)  

#add predictions and labels tocombined.df
results_df <- combined.df
stopifnot(nrow(results_df) == length(true_labels)) 
results_df$true_label <- true_labels
results_df$predicted_class <- predicted_classes
results_df$predicted_prob <- all_predictions

#classify by true label and pred class
results_df$confusion_type <- with(results_df, ifelse(true_label == 1 & predicted_class == 1, "TP",
                                                     ifelse(true_label == 0 & predicted_class == 1, "FP",
                                                            ifelse(true_label == 1 & predicted_class == 0, "FN", "TN"))))

#subset into TP, FP, FN, TN
true_positives   <- results_df[results_df$confusion_type == "TP", ]
false_positives  <- results_df[results_df$confusion_type == "FP", ]
false_negatives  <- results_df[results_df$confusion_type == "FN", ]
true_negatives   <- results_df[results_df$confusion_type == "TN", ]

#false negatives
nrow(false_negatives)
mean(false_negatives$TPM)
mean(false_negatives$log10_distance)
mean(false_negatives$H3K27ac)
sum(false_negatives$CpG_island=="CpG")
sum(false_negatives$CpG_island=="noCpG")
# false positives
nrow(false_positives)
mean(false_positives$TPM)
mean(false_positives$log10_distance)
mean(false_positives$H3K27ac)
sum(false_positives$CpG_island == "CpG")
sum(false_positives$CpG_island == "noCpG")
# true positives
nrow(true_positives)
mean(true_positives$TPM)
mean(true_positives$log10_distance)
mean(true_positives$H3K27ac)
sum(true_positives$CpG_island == "CpG")
sum(true_positives$CpG_island == "noCpG")
#true negatives
nrow(true_negatives)
mean(true_negatives$TPM)
mean(true_negatives$log10_distance)
mean(true_negatives$H3K27ac)
sum(true_negatives$CpG_island == "CpG")
sum(true_negatives$CpG_island == "noCpG")

negative_gi_esc_distance<-read.delim("ESC_samples_distance_match.tsv", sep = "\t", header = TRUE)
negative_gi_flc_distance<- read.delim("FLC_samples_distance_match.tsv", sep = "\t", header = TRUE)


negative_gi_flc_nodist<- read.delim("FLC_samples_no_distance_match.tsv", sep = "\t", header = TRUE)
negative_gi_esc_nodist<-read.delim("ESC_samples_no_distance_match.tsv", sep = "\t", header = TRUE)

#comparind datasets
length(unique(pe_esc.gi$anchor1.promoter.id))
length(unique(pe_flc.gi$anchor1.promoter.id))

length(unique(anchorTwo(pe_esc.gi)))
length(unique(anchorTwo(pe_flc.gi)))

#promoters per enhancer
mean(table(anchorTwo(pe_esc.gi)))  # ESC
mean(table(anchorTwo(pe_flc.gi)))  # FLC

#average enhancers per promoter
mean(table(anchorOne(pe_esc.gi)))  # ESC
mean(table(anchorOne(pe_flc.gi)))  # FLC


mean(mcols(pe_esc.gi)$distance)
mean(mcols(pe_flc.gi)$distance)


mean(pe_esc.gi$TPM, na.rm=TRUE)
mean(pe_flc.gi$TPM, na.rm=TRUE)

#correlation between TPM and disrance
cor(pe_esc.gi$TPM, log10(pe_esc.gi$distance), use="complete.obs")
cor(pe_flc.gi$TPM, log10(pe_flc.gi$distance), use="complete.obs")

#understanding expression vs promoter enhancer interactions
pe_esc.gi<-readRDS("ESC_samples_positive.rds")
negative_gi_esc<-readRDS("negative_interactions_esc.rds")
negative_no_dis_esc<-readRDS("negative_interactions_no_distance_esc.rds")

pe_flc.gi <- readRDS("FLC_samples_positive.rds")
negative_gi_flc <- readRDS("negative_interactions_flc.rds")
negative_no_dis_flc<-readRDS("negative_interactions_no_distance_flc.rds")

tpm

#tpm for interacting vs all
all_ESC <- TPMesc$TPM
esc_epi_promoters <- unique(pe_esc.gi$anchor1.promoter.id)
esc_TPM <- TPMesc$TPM[match(esc_epi_promoters, TPMesc$Geneid)]
esc_TPM[is.na(esc_TPM)] <- 0

median(esc_TPM)
median(all_ESC)
# calculate y_max for annotation
y_max <- max(log2(c(all_ESC, esc_TPM) + 1), na.rm = TRUE)

# boxplot with extra space at the top
boxplot(log2(all_ESC + 1), log2(esc_TPM + 1),
        names = c("All Promoters Expression", "EPI Promoter Expression"),
        ylab = "log2(TPM + 1)",
        main = "ESC Promoter Expression",
        col = c("grey80", "#56B4E9"),
        border = "black", outline = FALSE,
        ylim = c(0, y_max * 1.0))  # <-- add extra space

# significance line and stars
segments(1, y_max * 1.02, 2, y_max * 1.02)
text(1.5, y_max * 1.06, labels = "***", cex = 1.5)  # p < 0.001
#significance tests
wilcox.test(log2(all_ESC+1), log2(esc_TPM+1))
t.test(log2(all_ESC+1), log2(esc_TPM+1))



### FLC ###

all_FLC <- TPMflc$TPM
flc_epi_promoters <- unique(pe_flc.gi$anchor1.promoter.id)
flc_TPM <- TPMflc$TPM[match(flc_epi_promoters, TPMflc$Geneid)]
flc_TPM[is.na(flc_TPM)] <- 0

median(flc_TPM)
median(all_FLC)
# calculate y_max for annotation
y_max_flc <- max(log2(c(all_FLC, flc_TPM) + 1), na.rm = TRUE)

# boxplot with extra space at the top
boxplot(log2(all_FLC + 1), log2(flc_TPM + 1),
        names = c("All Promoters Expression", "EPI Promoter Expression"),
        ylab = "log2(TPM + 1)",
        main = "FLC Promoter Expression",
        col = c("grey80", "salmon"),
        border = "black", outline = FALSE,
        ylim = c(0, y_max_flc * 1.1))  # extra space

# significance line and stars
segments(1, y_max_flc * 1.02, 2, y_max_flc * 1.02)
text(1.5, y_max_flc * 1.06, labels = "***", cex = 1.5)  # p < 0.001

# significance tests
wilcox.test(log2(all_FLC + 1), log2(flc_TPM + 1))
t.test(log2(all_FLC + 1), log2(flc_TPM + 1))


#distance distributions
par(mfrow=c(1,2))

hist(log10(pe_esc.gi$distance), breaks=50,
     col="skyblue", main="ESC: Distance Distribution",
     xlab="log10 Distance (bp)")

hist(log10(pe_flc.gi$distance), breaks=50,
     col="salmon", main="FLC: Distance Distribution",
     xlab="log10 Distance (bp)")


