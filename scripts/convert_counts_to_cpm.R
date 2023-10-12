# Convert to CPM (counts per million)
cpm_func <- function(counts) {
    cpm <- apply(counts,2, function(x) (x/sum(x))*1000000)
    cpm = as.matrix(cpm)
    rownames(cpm) = rownames(counts)
    colnames(cpm) = colnames(counts)
    return(cpm)
}
counts = read.table("../counts/counts_matrix.tsv",sep="\t",header=TRUE, row.names=1)
cpm = cpm_func(counts[rowSums(counts)>0,])

# Write to file
write.table(cpm, file="../counts/cpm_matrix.tsv", sep="\t", quote=FALSE, row.names = TRUE)
                 
# Function to filter CPM > 2
cpm_greaterthan2_func <- function(cpm_0) {
    check_max_cpm = apply(cpm_0, 1, max)
    cpm_2 = cpm_0[which(check_max_cpm > 2), ]
    return(cpm_2)
}

# Apply CPM > 2 transformation
cpm_2 = cpm_greaterthan2_func(cpm)
print(colSums(cpm_2>0)) # count # miRNAs detected > 2 CPM

# Write CPM > 2 matrix to file
write.table(cpm_2, file="../counts/cpm_over2_matrix.tsv", sep="\t", quote=FALSE, row.names=TRUE)
