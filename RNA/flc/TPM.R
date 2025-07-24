#Read in feature counts (GSM was better qualiity)
setwd ("/mnt/clusters/admiral/data/c2007523/FYP")
# Load table
flcfc <- read.table("/mnt/clusters/admiral/data/c2007523/FYP/RNA/flc/featureCounts/GSM661638.featurecount.txt", header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
head(flcfc)
# Extract relevant columns
Geneid<- flcfc$Geneid
Length <- flcfc$Length
Counts <- flcfc$X.home.c2007523.mydata.c2007523.FYP.RNA.flc.star.GSM661638.unsort.Aligned.out.bam

# Reads Per Kilobase
rpk <- Counts / (Length / 1000)

# calculate scaling factor (sum of rpks / 1 million)
scalingfactor <- sum(rpk) / 1e6

# TPM
tpm <- rpk / scalingfactor

# Add to dataframe
flcfc$TPM <- tpm
TPMflc <- cbind(flcfc$Geneid, flcfc$TPM)
head(TPMesc)
colnames(TPMflc) <- c("Geneid", "TPM")
head (TPMflc)
setwd("/mnt/clusters/admiral/data/c2007523/FYP/RNA/flc")
write.table(TPMflc, file = "TPMflc_table.txt", sep = "\t", row.names = FALSE, quote = FALSE)
