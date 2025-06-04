import argparse
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Read metadata and counts
meta = pd.read_csv("../ref/mirna_metadata.csv")
meta['genotype'] = meta['genotype'].astype(str)
meta['genotype'] = pd.Categorical(meta['genotype'], categories=['58', '513'], ordered=True)
meta.set_index('sampleID', inplace=True)

counts = pd.read_csv("../counts/counts_matrix_gene_name.tsv", sep="\t", index_col=0)

counts = counts[counts.sum(axis=1) > 0] # remove non-expressed microRNAs


# Create DESeq2 dataset
dds = DeseqDataSet(counts=counts.T, 
                   metadata=meta,
                   design_factors='genotype', 
                   ref_level=['genotype', '58'],
                   refit_cooks=True)

# Perform DESeq2 analysis
dds.deseq2()

# Extract statistics
stat_res = DeseqStats(dds, contrast=['genotype', '513', '58'])
stat_res.summary()

# Extract results
result_df = stat_res.results_df.reset_index()


result_df.to_csv('test.csv')