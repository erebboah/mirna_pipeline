# MicroRNA-seq
## Overview
MicroRNA-seq libraries are built using the [ENCODE protocol](https://www.encodeproject.org/documents/49f43842-5ab4-4aa1-a6f4-2b1234955d93/@@download/attachment/Multiplexed%20miRNA%20Sequencing%20Library%20Generation%20Protocol-3.4.pdf) "Multiplexed miRNA Sequencing Library Generation Protocol".

Briefly, libraries are built from 300 ng to 1 ug of total RNA. Adapters are ligated to the 5’ and 3’ ends of small RNAs using their 5’ phosphate and 3’ hydroxyl groups, then the ligation product is reverse transcribed. The 5’ adapter adds a 6-nucleotide barcode from a one of 7 sets of 4 distinct barcodes used in downstream demultiplexing. The cDNA is amplified with 58 bp reverse and 55 bp forward primers containing additional 6-nucleotide barcodes added to the 3’ end. The 140 bp product containing mature microRNA (21-25 nucleotides) is size-selected using 10% TBE-Urea gel. 

Libraries are isolated from the TBE-Urea gel by agitated incubation at 70°C, 1000 RPM for 2 hours in a buffer containing 0.5 M ammonium acetate, 0.1% SDS, and 0.1 uM EDTA, then precipitated overnight in 50% isopropanol. 

Libraries are quantified using Qubit dsDNA HS Assay Kit and sequenced on an Illumina NextSeq 2000 with P1 or P2 100 cycle kits as 75 bp single-end reads (SE 75/0/6/0) to ~10 M raw read depth per library. We have good results sequencing pooled library at 850 pM (%Loading Concentration > 95). Submission to the ENCODE portal requires >5 M aligned reads, >300 microRNAs detected at >2 CPM, and a Spearman replicate correlation >0.85.

<img src="https://github.com/erebboah/mirna_pipeline/blob/master/mirna_overview.png" width="602" height="502">

## Environment requirements
#### Demultiplexing
- seqtk (`conda install -c bioconda seqtk`)
- python, with `os` and `sys` libraries installed

#### Quantification
Nothing! Can module load STAR and cutadapt on HPC.

#### Analysis
- R with `tidyverse`, `rtracklayer` (to grab gene names from miRNA GTF)
- python with `pandas`, `seaborn`, `matplotlib`, `sklearn`

## Demultiplexing

Library structure deviates from standard Illumina, so samples must be demultiplexed with custom code. Illumina NextSeq2000 SampleSheet (e.g. for our [practice run](https://github.com/erebboah/mirna_pipeline/blob/master/ref/SampleSheet_1080.csv)) needs to specify read setup but the SampleID / Index doesn't matter since we're just going to use the "undetermined" reads.

1. After sequencing is complete, transfer reads to HPC.

```
scp /home/sharing/runs/231006_VH00582_82_AACYV55M5/Analysis/1/Data/fastq/*.fastq.gz erebboah@hpc3.rcic.uci.edu:/pub/erebboah/mirna_pipeline/fastq/run_1080
```

2. Concatenate all reads to undetermined.fastq.gz and unzip, making sure you have enough space. Need 24G for a P1 100 cycle kit, x2 since you'll basically be copying it to demultiplex = 48G. If you clean up all files at the end (undetermined.fastq.gz and all the individual sample output directories, just keeping the counts directory), you'll only need ~5G, mostly for the demultiplexed fastqs.

```
cd /pub/erebboah/mirna_pipeline/fastq/
cat *.fastq.gz > undetermined.fastq.gz
gunzip undetermined.fastq.gz
```

3. Edit the "real" sample sheet. The first column is the 5' Adaptor ID (1-7) and forward primer index ID (1-12, 14, 15, 18, 21-25) separated by an underscore. The second column is your desired fastq name (sample ID). No header! Preferably matches some sample ID in your metadata. Should be in the same directory as the scripts. [Example of samplesheet.csv](https://github.com/erebboah/mirna_pipeline/blob/master/scripts/samplesheet.csv)

4. Edit [demux_mirna.sh](https://github.com/erebboah/mirna_pipeline/blob/master/scripts/demux_mirna.sh) so that you are using your HPC account (change `#SBATCH -A SEYEDAM_LAB`) and your conda environment (change `source ~/miniconda3/bin/activate seqtkpython3`). As inputs, you need fastq named `undetermined.fastq` in the fastq directory of this repo (e.g. `/pub/erebboah/mirna_pipeline/fastq/`) and your [samplesheet.csv](https://github.com/erebboah/mirna_pipeline/blob/master/scripts/samplesheet.csv) in the scripts directory (e.g. `/pub/erebboah/mirna_pipeline/scripts/`). 

5. Run demultiplexing bash script: `sbatch demux_mirna.sh`. Output is demultiplexed, gzipped fastqs with sample IDs in samplesheet.csv, in the fastq directory along with undetermined.fastq.gz.

## Quantification
### Make STAR reference - only have to do this once
1. Download microRNA GENCODE GTFs from ENCODE portal: [vM21 mouse](https://www.encodeproject.org/files/ENCFF094ICJ/) or [v29 human](https://www.encodeproject.org/files/ENCFF470CZH/).

```
wget https://www.encodeproject.org/files/ENCFF094ICJ/@@download/ENCFF094ICJ.gtf.gz
wget https://www.encodeproject.org/files/ENCFF470CZH/@@download/ENCFF470CZH.gtf.gz
mv ENCFF094ICJ.gtf.gz mm10_mirna.gtf.gz
mv ENCFF470CZH.gtf.gz hg38_mirna.gtf.gz
gunzip *.gtf.gz
```

2. Download reference sequences from the portal: [mm10 mouse](https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/) or [hg38 human](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/).

```
wget https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz

gunzip *.fasta.gz
```

3. Once you have the reference files, run STAR in genomeGenerate mode: `sbatch make_ref.sh`

### Trim and map reads
Run cutadapt to trim adapters and STAR to map. Specify human or mouse in sbatch statement. 

Inputs:
   - Demultiplexed and gzipped fastqs in `fastq`
   - STAR reference directory, e.g. `ref/mm10` or `ref/hg38`
   - samplesheet.csv in `scripts` to indicate the samples you want to map

Human: `sbatch trim_map.sh hg38`

Mouse: `sbatch trim_map.sh mm10`
  
Each sample will get its own output directory, named the same as the sample ID, e.g. `ENC4_453_SB`. Within the sample directory, there's `cutadapt` and `star` directories containing intermediate files. The actual tab-separated microRNA quantifications per sample are in the main counts directory `counts` (e.g. `/pub/erebboah/mirna_pipeline/counts/ENC4_453_SB.tsv`). Remove the sample directories to save space if you only need the final counts. 

Other than counts, you may be interested in the STAR report (e.g. `ENC4_453_SB/star/Log.final.out`) and signal files to display on the UCSC genome browser (e.g. `ENC4_453_SB/star/Signal.UniqueMultiple.str1.out.wig`). The percent uniquely mapped reads is low because the microRNAs are so short. On average I get ~43% multi-mapped reads, and pass ENCODE standards. For microRNA-seq, the number of multi-mapped reads plus unique reads are the total number of aligned reads (should be > 5M). Feel free to poke around my old [ENCODE miRNA-seq spreadsheet](https://docs.google.com/spreadsheets/d/1qcve4QnxcMVTgyIxT3ouhblCO1y5pmKZ2B9L-TCnbBg/edit#gid=898072680) for other QC stuff.

## Analysis
1. Concatenate counts per sample into a counts matrix. Example code in [R](https://github.com/erebboah/mirna_pipeline/blob/master/scripts/make_counts_matrix.R) and [python](https://github.com/erebboah/mirna_pipeline/blob/master/scripts/make_counts_matrix.py) using data from our practice run.
2. Convert counts to CPM (counts per million) to normalize for library depth: [R](https://github.com/erebboah/mirna_pipeline/blob/master/scripts/convert_counts_to_cpm.R) and [python](https://github.com/erebboah/mirna_pipeline/blob/master/scripts/convert_counts_to_cpm.py)
3. Principal component analysis: [R](https://github.com/erebboah/mirna_pipeline/blob/master/scripts/make_pca.R) and [python](https://github.com/erebboah/mirna_pipeline/blob/master/scripts/make_pca.py)
4. Differential microRNA expression analysis between conditions using pyDeseq2: python
5. Heatmap and volcano plot of differentially expressed microRNAs: R and python

## Summary
1. `sbatch demux_mirna.sh`
2. `sbatch make_ref.sh` (only need to run once to generate references!)
3. `sbatch trim_map.sh GRCh38` (e.g. for human samples)
4. Analysis in R and/or python.
