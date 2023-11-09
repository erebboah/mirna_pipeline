import argparse
import pandas as pd
import pyranges as pr
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Function to run DESeq2 analysis
def run_deseq2(counts, meta, group, filtering):
    dds = DeseqDataSet(
        counts=counts.T,
        clinical=meta,
        design_factors=group,
        refit_cooks=True
    )
    dds.deseq2()
    stat_res = DeseqStats(dds, contrast=[group] + filtering[group])
    stat_res.summary()
    
    # Get DEG results
    results = stat_res.results_df.copy(deep=True)
    results = results.reset_index()  # put gene_id in column
    
    return results

# Function to merge results with gene names
def merge_results_with_gene_names(results, mirna_gtf):
    results = pd.merge(results, mirna_gtf, on="gene_id", how="left")
    return results

# Function to save results to CSV
def save_results_to_csv(results, output_file):
    results.to_csv(f"../degs/{output_file}.csv", index=False)

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Run DESeq2 analysis and save results to CSV.')
    parser.add_argument('--sex', nargs='+', required=True, help='List of sexes to include')
    parser.add_argument('--technician', nargs='+', required=True, help='List of technicians to include')
    parser.add_argument('--timepoint', nargs='+', required=True, help='List of timepoints to include')
    parser.add_argument('--group', required=True, choices=['timepoint', 'sex'], help='Group to test')
    parser.add_argument('--output', required=True, help='Output file name')
    
    args = parser.parse_args()

    # Set up filtering dictionary
    filtering = {'sex': args.sex,
                 'timepoint': args.timepoint,
                 'technician': args.technician}
    
    print(filtering)

    # Read metadata and counts
    meta = pd.read_csv("../ref/mirna_metadata.csv")
    counts = pd.read_csv("../counts/counts_matrix.tsv", sep="\t", index_col=0)
    
    counts = counts[counts.sum(axis=1) > 0] # remove non-expressed microRNAs
    
    meta_selected = meta.copy(deep=True)
    for col in filtering.keys():
        meta_selected = meta_selected[meta_selected[col].isin(filtering[col])]

    meta_selected.index = meta_selected['sampleID'] # deseq2 needs metadata to be indexed

    counts_selected = counts.loc[:, meta_selected['sampleID']]
    counts_selected[counts_selected.columns] = counts_selected[counts_selected.columns].astype(int)
    
    # Run DESeq2 analysis
    results = run_deseq2(counts_selected, meta_selected, args.group, filtering)

    # Read GTF file into PyRanges
    mirna_gtf = pr.read_gtf('../ref/hg38_mirna.gtf', as_df=True)

    # Extract gene_id and gene_name columns
    mirna_gtf = mirna_gtf[['gene_id', 'gene_name']].drop_duplicates()
    
    # Merge results with gene names
    results = merge_results_with_gene_names(results, mirna_gtf)

    # Save results to CSV
    save_results_to_csv(results, args.output)

if __name__ == "__main__":
    main()
