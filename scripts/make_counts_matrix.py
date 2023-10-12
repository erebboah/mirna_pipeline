import pandas as pd

meta = pd.read_csv("../ref/mirna_practice_metadata.csv")
path = "../counts/"
samples = meta["sampleID"]

# Read in first sample
star_output = pd.read_table(f"{path}{samples.iloc[0]}.tsv", sep="\t", header=None)
# remove STAR output header lines 1-4 and take first two columns: 
# column 1 = gene ID, column 2 = unstranded counts
counts = star_output.iloc[4:, [0, 1]]

# Add in the rest of the samples
for i in range(1, len(samples)):
    star_output = pd.read_table(f"{path}{samples.iloc[i]}.tsv", sep="\t", header=None)
    star_output = star_output.loc[star_output[0].isin(counts[0].dropna()), :]
    counts = pd.concat([counts, star_output[1]], axis=1)

# Rename columns
counts.columns = ["gene_id"] + samples.tolist()

# Write to file
counts.to_csv("../counts/counts_matrix.tsv", sep="\t", index=False)