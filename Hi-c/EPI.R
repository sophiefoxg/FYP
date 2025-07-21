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
po.gi=GenomicInteractions(esc1.gr,esc2.gr, counts = x1$raw.count, loe=x1$log.observed.expected.)
head(po.gi)

# distance between promoter other interactions, using anchor1 and anchor 2
summary(calculateDistances(po.gi))

#^ look into max distance this is large
hist(log10(calculateDistances(po.gi)), xlab = "Distance log10(Mb)", ylab= "Frequency", main = "Promoter-Other Distance (ESC)" )




#adding in hits from chip-seq data
setwd("/mnt/clusters/admiral/data/c2007523/FYP/chip/esc/peaks")
if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  BiocManager::install("rtracklayer")
}

library(rtracklayer)
library(GenomicRanges)

# Import BED as GRanges
esc_chip <- import("ESC.bed", format = "BED")
head(po.gi)
hits_esc<-findOverlaps(anchorTwo(po.gi),esc_chip)
po.gi$anchor2.node.class[queryHits(hits_esc)  ]  = "enhancer"
table(po.gi$anchor2.node.class)

head(po.gi)



#broad peaks
#setwd("/mnt/clusters/admiral/data/c2007523/FYP/chip/esc/broadpeaks")
#esc_broad<-import("broadpeaksesc.bed")
#hits_escbroad<-findOverlaps(anchorTwo(po.gi),esc_broad)


#po.gi$anchor2.node.class[queryHits(hits_escbroad)] = "enhancer"
#table(po.gi$anchor2.node.class)

#7023 promoter enhancer interactions
pe.gi<-po.gi[po.gi$anchor2.node.class== "enhancer"]
head(pe.gi)
pe.gi = pe.gi[pe.gi$anchor1.promoter.id!="not_promoter"]
table(pe.gi$anchor2.node.class)
#enhancer 6374 


length(unique(pe.gi$anchor1.promoter.id))

pe.gi=pe.gi[nchar(pe.gi$anchor1.promoter.id) <= 18, ]

head(pe.gi)
mcols(pe.gi)$distance <- calculateDistances(pe.gi)
hist(calculateDistances(pe.gi))


#read in transcription data
setwd("/mnt/clusters/admiral/data/c2007523/FYP/RNA/esc/")
TPMesc<-read.delim("TPMesc_table.txt", header=TRUE)
head(TPMesc)
matchRNA<-match(pe.gi$anchor1.promoter.id,TPMesc$Geneid)
pe.gi$TPM<-TPMesc$TPM[matchRNA]

head(pe.gi)


enhancer.gr <- anchorTwo(pe.gi)  
head(enhancer.gr)
#if TPM is higher than 0 then inctive  add as feature
# add cpg for promoter type


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

pp.gi=GenomicInteractions(esc3.gr, esc4.gr, counts=x2$raw.count,loe=x2$log.observed.expected.)
head(pp.gi)
promoters.gr = unique(c(esc1.gr[esc1.gr$node.class=="promoter"],
                        esc3.gr[esc3.gr$node.class=="promoter"],
                        esc4.gr[esc4.gr$node.class=="promoter"]))

head(promoters.gr)
promoters.gr = promoters.gr[promoters.gr$promoter.id!="not_promoter"]
promoters.gr=promoters.gr[nchar(promoters.gr$promoter.id) <= 18, ]
#distance between promoter promoter interactions df

hist(log10(calculateDistances(pp.gi)),xlab = "Distance log10(Mb)", ylab= "Frequency", main = "Promoter-Promoter Distance (ESC)" )
mcols(pp.gi)$distance <- calculateDistances(pp.gi)



#removing all promoters that are "not_promoters"
pp.gi = pp.gi[pp.gi$anchor1.promoter.id!="not_promoter",]
pp.gi = pp.gi[pp.gi$anchor2.promoter.id!="not_promoter",]

head(pp.gi)
#removing all promoters that span several genes


pp.gi <- pp.gi[nchar(pp.gi$anchor1.promoter.id) <= 18, ]
pp.gi<- pp.gi[nchar(pp.gi$anchor2.promoter.id) <=18,]

head(pp.gi)

#read in tPM data
setwd("/mnt/clusters/admiral/data/c2007523/FYP/RNA/esc/")
TPMesc<-read.delim("TPMesc_table.txt", header=TRUE)
head(TPMesc)
matchRNA1<-match(pp.gi$anchor1.promoter.id,TPMesc$Geneid)
matchRNA2<-match(pp.gi$anchor2.promoter.id,TPMesc$Geneid)
pp.gi$TPM1<-TPMesc$TPM[matchRNA1]
pp.gi$TPM2<-TPMesc$TPM[matchRNA2]


View(as.data.frame(pp.gi))
pp.gi$distances<-(calculateDistances(pp.gi))

###making fragments

head(po.gi)
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

head(gaps.gr)
head(baits.gr)
#find overlap between hindIII and baits, taking away the promoters from the "gaps"
hits<-findOverlaps(fragments,baits.gr)
gaps.gr<-fragments[-queryHits(hits)]

head(gaps.gr)
length(gaps.gr)

#posiyive set of promoters
promoter_pos <- anchorOne(pe.gi)


#removing positive set from the fragments
promoter_hits <- findOverlaps(gaps.gr, promoter_pos)
length(promoter_hits)
fragments_no_promoters <- gaps.gr[-queryHits(promoter_hits)]

##no baits or known promoter anchors
gaps.gr<-fragments_no_promoters



##### sampling negative set
# define a promoter
promoter<-anchorOne(pe.gi[475])
#define 1.5mb window with a lower and upper bound 
#from the promoters start and end site
window<-1.5e6
lower<-start(promoter)-window
upper<-end(promoter)+window

# create ranges object for Region Of Interest 
ROI.gr<-GRanges(seqnames=seqnames(promoter),ranges=IRanges(start=lower, end=upper))

#finding this range in gaps.gr (gaps.gr is the hind3bed file without the promoters)
ROI_in_gaps<-findOverlaps(ROI.gr, gaps.gr)
ROI.gr<-gaps.gr[subjectHits(ROI_in_gaps)]
head(ROI.gr)

#checking what enhancers this promoter interacts with if there are multiple
promoter_hits<-which(overlapsAny(anchorOne(pe.gi), promoter))
interacting_enhancers<-anchorTwo(pe.gi)[promoter_hits]

#removing any interacting enhancers that overlap the region of interest
remove_enhancers<-findOverlaps(ROI.gr, interacting_enhancers)
ROI_no_enhancers<-ROI.gr[-queryHits(remove_enhancers)]

#sample from this region of interest with no enhancers
negative_p475<-sample(ROI_no_enhancers,10,replace=FALSE)
View(as.data.frame(negative_p475))
View(as.data.frame(ROI_no_enhancers))


#make it a loop
##### sampling negative set
# Prepare output list
negative_samples_list <- list()

# Loop through all promoters in pe.gi
for (i in seq_along(pe.gi)) {
  
  # define a promoter
  promoter <- anchorOne(pe.gi[i])
  
  # define 1.5mb window with a lower and upper bound 
  # from the promoters start and end site
  window <- 1.5e6
  lower <- start(promoter) - window
  upper <- end(promoter) + window
  
  # create ranges object for Region Of Interest 
  ROI.gr <- GRanges(seqnames = seqnames(promoter), ranges = IRanges(start = lower, end = upper))
  
  # finding this range in gaps.gr (gaps.gr is the hind3bed file without the promoters)
  ROI_in_gaps <- findOverlaps(ROI.gr, gaps.gr)
  ROI.gr <- gaps.gr[subjectHits(ROI_in_gaps)]
  head(ROI.gr)
  
  # checking what enhancers this promoter interacts with if there are multiple
  promoter_hits <- which(overlapsAny(anchorOne(pe.gi), promoter))
  interacting_enhancers <- anchorTwo(pe.gi)[promoter_hits]
  
  # removing any interacting enhancers that overlap the region of interest
  remove_enhancers <- findOverlaps(ROI.gr, interacting_enhancers)
  ROI_no_enhancers <- ROI.gr[-queryHits(remove_enhancers)]
  
  # sample from this region of interest with no enhancers
  if (length(ROI_no_enhancers) >= 10) {
    negative_sample <- sample(ROI_no_enhancers, 10, replace = FALSE)
  } else {
    negative_sample <- ROI_no_enhancers
  }
  
  # save result
  negative_samples_list[[i]] <- negative_sample
}

# Example view
View(as.data.frame(negative_samples_list[[1]]))









#create an empty list
negatives<-list()
#1.5mb around the promoter
window<-1.5e6
#for every line in pe.gi
for (i in seq_len(length(pe.gi))) {
#define the promoter
promoter<-anchorOne(pe.gi[i])
#define the distance between the promoter and enhancer
dist<-mcols(pe.gi)$distance[i]
#+ or - 1.5 mb either ide
lower<-dist-window
upper<-dist+window
#on the same chromosome, i can change this to do my test(for seperate chromosomes)
chr<-gaps.gr[as.character(seqnames(gaps.gr))==as.character(seqnames(promoter))]
#picking a fragment thats same distance to enhancer
gap_dist<-abs(start(chr)-start(promoter))
#within the distance
gap<-chr[gap_dist>=lower&gap_dist<=upper]

#loop using sample()
if (length(gap) > 0) {
sampled<-sample(gap, 1)
negative<-GenomicInteractions(promoter, sampled)
mcols(negative)$distance<-abs(start(promoter)-start(sampled))
negatives[[length(negatives)+1]]<-negative
  }
}

# make the GenomicInteractions object
negatives_list<-do.call(c, negatives)
head(negatives_list)
head(pe.gi)
length(pe.gi)
length(negatives_list)

#store distance? and expression levels, annotate with promoter enhancer
