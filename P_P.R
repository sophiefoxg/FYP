if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicInteractions")
library(GenomicInteractions)
x2 = read.delim("/mnt/clusters/admiral/data/c2007523/Diss/ESC_promoter_promoter_significant_interactions.txt")

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

#distance between promoter promoter interactions df
ppdistances=(as.data.frame(calculateDistances(pp.gi)))
View(ppdistances)
colnames(ppdistances) <- "distance"


pp.df= as.data.frame(pp.gi)
ppdis.df= cbind(pp.df, ppdistances)

View(as.data.frame(ppdis.df))
# look into max distance this is large
hist(log10(calculateDistances(pp.gi)),xlab = "Distance log10(Mb)", ylab= "Frequency", main = "Promoter-Promoter Distance (ESC)" )


#distance between promoter other interactions df
ppdistances=(as.data.frame(calculateDistances(pp.gi)))
View(ppdistances)
colnames(ppdistances) <- "distance"

pp.df= as.data.frame(pp.gi)
ppdis.df= cbind(pp.df, ppdistances)

View(as.data.frame(ppdis.df))

#removing all promoters that are "not_promoters"
ppdis.df = ppdis.df[ppdis.df$anchor1.promoter.id!="not_promoter",]
ppdis.df = ppdis.df[ppdis.df$anchor2.promoter.id!="not_promoter",]
head(ppdis.df)

#removing all promoters that span several genes

nchar(ppdis.df$anchor1.promoter.id, type="chars")
ppdis.df <- ppdis.df[nchar(ppdis.df$anchor1.promoter.id) <= 18, ]

nchar(ppdis.df$anchor2.promoter.id, type="chars")
ppdis.df<- ppdis.df[nchar(ppdis.df$anchor2.promoter.id) <=18,]

head(ppdis.df)

####subsetting for largest distances####

#ordering by distance (ordered by distance promoter other)
obdpp<- ppdis.df[order(ppdis.df$distance, decreasing = TRUE),]

#top100

obdpptop<-as.data.frame(obdpp[1:100,])

head(obdpptop)
View(as.data.frame(obdpptop))

p1.v<-obdpptop$anchor1.promoter.id

p2<-obdpptop$anchor2.promoter.id

#what do those genes do?
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
library("clusterProfiler")
library("org.Mm.eg.db")


head(p1.v)

P1GO <- enrichGO(gene = p1.v, 
                 keyType = "ENSEMBL",
                 OrgDb = "org.Mm.eg.db", 
                 ont = "BP", 
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.05,
)
View(as.data.frame(P1GO))



