### install packages and data read in ###

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicInteractions")
library(GenomicInteractions)

x1 = read.delim("/mnt/clusters/admiral/data/c2007523/Diss/ESC_promoter_other_significant_interactions.txt")
x2 = read.delim("/mnt/clusters/admiral/data/c2007523/Diss/ESC_promoter_promoter_significant_interactions.txt")


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


#all promoters

promoters.gr= unique(c(esc1.gr[esc1.gr$node.class == "promoter"],
                       esc3.gr[esc3.gr$node.class=="promoter"],
                       esc4.gr[esc4.gr$node.class=="promoter"]
))
#removing all promoters that are "not_promoters"
promoters.gr = promoters.gr[promoters.gr$promoter.id!="not_promoter"]


#removing all promoters that span several genes


head(esc1.gr)
head(esc2.gr)


#promoter and other dataset
po.gi=GenomicInteractions(esc1.gr,esc2.gr, counts = x1$raw.count, loe=x1$log.observed.expected.)
head(po.gi)



#promoter and promoter dataset

pp.gi=GenomicInteractions(esc3.gr, esc4.gr, counts=x2$raw.count,loe=x2$log.observed.expected.)
head(pp.gi)
# all interactions
all.gi= c(po.gi,pp.gi)


# distance between promoter other interactions, using anchor1 and anchor 2
summary(calculateDistances(po.gi))
#^ look into max distance this is large
hist(log10(calculateDistances(po.gi)), xlab = "Distance log10(Mb)", ylab= "Frequency", main = "Promoter-Other Distance (ESC)" )








#distance between promoter other interactions df
podistances=(as.data.frame(calculateDistances(po.gi)))
View(podistances)
colnames(podistances) <- "distance"

po.df= as.data.frame(po.gi)
podis.df= cbind(po.df, podistances)

View(as.data.frame(podis.df))

#removing all promoters that are "not_promoters"
podis.df = podis.df[podis.df$anchor1.promoter.id!="not_promoter",]
head(podis.df)

#removing all promoters that span several genes

nchar(podis.df$anchor1.promoter.id, type="chars")
podis.df <- podis.df[nchar(podis.df$anchor1.promoter.id) <= 18, ]
head(podis.df)

####subsetting for largest distances####

#ordering by distance (ordered by distance promoter other)
obdpo<- podis.df[order(podis.df$distance, decreasing = TRUE),]

#top100

obdpotop<-as.data.frame(obdpo[1:99,])

head(obdpotop)
View(as.data.frame(obdpotop))
top.po<-obdpotop$anchor1.promoter.id


#what do those genes do?
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
library("clusterProfiler")
library("org.Mm.eg.db")


head(top.po)

POGO <- enrichGO(gene = top.po, 
                 keyType = "ENSEMBL",
                 OrgDb = "org.Mm.eg.db", 
                 ont = "BP", 
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.05, 
)
View(as.data.frame(POGO))

POenrich<-cbind(POGO$Description,obdpotop)

head(POenrich)


#top100 pp

#ordering by distance (ordered by distance promoter promoter)
obdpp<- ppdis.df[order(ppdis.df$distance, decreasing = TRUE),]
obdpptop<-as.data.frame(obdpp[1:100,])
View(obdpp)

View(obdpptop)


#what do those genes do?
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
library("clusterProfiler")
library("org.Mm.eg.db")

topp1<-obdpptop$anchor1.promoter.id[obdpptop$anchor1.promoter.id!="not_promoter"]

head(topp1)
View(as.data.frame(topp1))
nchar(topp1, type="chars")
topp1<- topp1[nchar(topp1) <= 18 ]

PP1 <- enrichGO(gene = topp1, 
                 keyType = "ENSEMBL",
                 OrgDb = "org.Mm.eg.db", 
                 ont = "BP", 
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.05, 
)
View(as.data.frame(PP1))


topp2<-obdpptop$anchor2.promoter.id[obdpptop$anchor2.promoter.id!="not_promoter"]

head(topp2)
nchar(topp2, type="chars")
topp2<- topp2[nchar(topp2) <= 18 ]

PP2 <- enrichGO(gene = topp2, 
                keyType = "ENSEMBL",
                OrgDb = "org.Mm.eg.db", 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
)
View(as.data.frame(PP2))

