#Read in feature counts (GSM was better qualiity)
setwd ("/mnt/clusters/admiral/data/c2007523/FYP")
# Load table
escfc <- read.table("/mnt/clusters/admiral/data/c2007523/FYP/RNA/esc/featureCounts/GSM723776.featurecount.txt", header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
head(escfc)
# Extract relevant columns
Geneid<- escfc$Geneid
Length <- escfc$Length
Counts <- escfc$X.home.c2007523.mydata.c2007523.FYP.RNA.esc.star.GSM723776.unsort.Aligned.out.bam

# Reads Per Kilobase
rpk <- Counts / (Length / 1000)

# calculate scaling factor (sum of rpks / 1 million)
scalingfactor <- sum(rpk) / 1e6

# TPM
tpm <- rpk / scalingfactor

# Add to dataframe
escfc$TPM <- tpm
TPMesc <- cbind(escfc$Geneid, escfc$TPM)
head(TPMesc)
colnames(TPMesc) <- c("Geneid", "TPM")
head (TPMesc)
setwd("/mnt/clusters/admiral/data/c2007523/FYP/RNA/esc")
write.table(TPMesc, file = "TPMesc_table.txt", sep = "\t", row.names = FALSE, quote = FALSE)


