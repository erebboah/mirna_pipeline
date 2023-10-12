import pandas as pd

def cpm_func(counts):
    cpm = counts.apply(lambda x: (x / x.sum()) * 1e6, axis=0)
    return cpm

def cpm_greaterthan2_func(cpm_0):
    check_max_cpm = cpm_0.apply(max, axis=1)
    cpm_2 = cpm_0[check_max_cpm > 2]
    return cpm_2

# Read counts matrix
counts = pd.read_csv("../counts/counts_matrix.tsv", sep="\t", header=0, index_col=0)

# Apply CPM transformation
cpm = cpm_func(counts.loc[counts.sum(axis=1) > 0])

# Write CPM matrix to file
cpm.to_csv("../counts/cpm_matrix.tsv", sep="\t", index=True)

# Apply CPM > 2 transformation
cpm_2 = cpm_greaterthan2_func(cpm)

print((cpm_2 > 0).sum(axis=0))  # count # miRNAs detected > 2 CPM

# Write CPM > 2 matrix to file
cpm_2.to_csv("../counts/cpm_over2_matrix.tsv", sep="\t", index=True)
