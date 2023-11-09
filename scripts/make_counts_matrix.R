meta = read.csv("../ref/mirna_metadata.csv")
path = "../counts/"
samples = meta$sampleID

# Read in first sample
star_output = read.table(paste0(path,samples[1],".tsv"),sep="",header=FALSE)
# remove STAR output header lines 1-4 and take first two columns: 
# column 1 = gene ID, column 2 = unstranded counts
counts = star_output[-c(1:4),c(1,2)]

# Add in rest of the samples
for(i in 2:length(samples)){
    star_output = read.table(paste0(path,samples[i],".tsv"),sep="",header=FALSE)
    star_output = star_output[na.omit(match(counts$V1,star_output$V1)),]
    counts=cbind(counts,star_output[,"V2"])
}

# Rename columns
colnames(counts) = c("gene_id",samples)

# Write to file
write.table(counts, file="../counts/counts_matrix.tsv", sep="\t", quote=F, row.names = F)
