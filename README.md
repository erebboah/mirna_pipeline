# MicroRNA-seq Analysis
## Overview
MicroRNA-seq libraries are built using the [ENCODE protocol](https://www.encodeproject.org/documents/49f43842-5ab4-4aa1-a6f4-2b1234955d93/@@download/attachment/Multiplexed%20miRNA%20Sequencing%20Library%20Generation%20Protocol-3.4.pdf) "Multiplexed miRNA Sequencing Library Generation Protocol".

Briefly, libraries are built from 300 ng to 1 ug of total RNA. Adapters are ligated to the 5’ and 3’ ends of small RNAs using their 5’ phosphate and 3’ hydroxyl groups, then the ligation product is reverse transcribed. The 5’ adapter adds a 6-nucleotide barcode from a one of 7 sets of 4 distinct barcodes used in downstream demultiplexing. The cDNA is amplified with 58 bp reverse and 55 bp forward primers containing additional 6-nucleotide barcodes added to the 3’ end. The 140 bp product containing mature microRNA (21-25 nucleotides) is size-selected using 10% TBE-Urea gel. 

Libraries are isolated from the TBE-Urea gel by agitated incubation at 70°C, 1000 RPM for 2 hours in a buffer containing 0.5 M ammonium acetate, 0.1% SDS, and 0.1 uM EDTA, then precipitated overnight in 50% isopropanol. 

Resulting microRNA-seq libraries are quantified using Qubit dsDNA HS Assay Kit and sequenced on an Illumina NextSeq 2000 with P2 100 cycle kits as 75 bp single-end reads (SE 75/0/6/0) to ~10 M raw read depth per library. Submission to the ENCODE portal requires >5 M aligned reads, >300 microRNAs detected at >2 CPM, and a Spearman replicate correlation >0.85.

<img src="https://github.com/erebboah/mirna_pipeline/blob/master/mirna_overview.png" width="602" height="502">

## Environment requirements
Run the demultiplexing code in an environment with
- seqtk
- python, with `os` and `sys` libraries installed 

## Demultiplexing
Library structure deviates from standard Illumina, so samples must be demultiplexed with custom code. Illumina NextSeq2000 SampleSheet (e.g. for our [practice run](https://github.com/erebboah/mirna_pipeline/blob/master/ref/SampleSheet_1080.csv)) needs to specify read setup but the SampleID / Index doesn't matter since we're just going to use the "undetermined" reads.

1. After sequencing is complete, transfer reads to HPC.

```
scp /home/sharing/runs/231006_VH00582_82_AACYV55M5/Analysis/1/Data/fastq/*.fastq.gz erebboah@hpc3.rcic.uci.edu:/pub/erebboah/mirna_pipeline/fastq/run_1080
```

2. Concatenate all reads to undetermined.fastq.gz and unzip, making sure you have enough space. Need 24G for a P1 100 cycle kit. 

```
cd /pub/erebboah/mirna_pipeline/fastq/
cat *.fastq.gz > undetermined.fastq.gz
gunzip undetermined.fastq.gz
```

3. Edit the "real" sample sheet. The first column is the 5' Adaptor ID (1-7) and forward primer index ID (1-12, 14, 15, 18, 21-25) separated by an underscore. The second column is your desired fastq name (sample ID). No header! Preferably matches some sample ID in your metadata. Should be in the same directory as the scripts. [Example of samplesheet.csv](https://github.com/erebboah/mirna_pipeline/blob/master/scripts/samplesheet.csv)

4. Edit [demux_mirna.sh](https://github.com/erebboah/mirna_pipeline/blob/master/scripts/demux_mirna.sh) so that you are using your HPC account (change `#SBATCH -A SEYEDAM_LAB`) and your conda environment (change `source ~/miniconda3/bin/activate seqtkpython3`). As inputs, you need fastq named `undetermined.fastq` in the fastq directory of this repo (e.g. `/pub/erebboah/mirna_pipeline/fastq/`) and your [samplesheet.csv](https://github.com/erebboah/mirna_pipeline/blob/master/scripts/samplesheet.csv) in the scripts directory (e.g. `/pub/erebboah/mirna_pipeline/scripts/`). 

5. Run demultiplexing bash script: `sbatch demux_mirna.sh`. Output is demultiplexed fastqs named as in samplesheet.csv in the fastq directory.

## Quantification
### Make STAR reference - only have to do once
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
Run cutadapt to trim adapters and STAR to map. Specify human or mouse in sbatch call. 

Human: `sbatch trim_map.sh GRCh38`

Mouse: `sbatch trim_map.sh mm10`

Inputs should have been generated in previous steps. You need:
   - Demultiplexed and gzipped fastqs in `fastq`
   - STAR reference directory, e.g. `ref/mm10` or `ref/hg38`
   - samplesheet.csv in `scripts`
  
Each sample will get its own output directory, named the same as the sample ID, e.g. `ENC4_453_NM`. Within the sample directory, there's `cutadapt` and `star` directories containing intermediate files. The actual tab-separated microRNA quantifications per sample are in the main counts directory `counts` (e.g. `/pub/erebboah/mirna_pipeline/counts/ENC4_453_NM.tsv`). Feel free to remove the sample directories to save space if you only need the final counts.

## Analysis
1. Concatenate counts per sample into a counts matrix. Example code in R and python from our practice run.
2. Convert counts to CPM (counts per million) to normalize for library depth: R and python
3. Principal component analysis: R and python
4. Differential microRNA expression analysis between conditions using pyDeseq2: python
5. Heatmap and volcano plot of differentially expressed microRNAs: R and python

## Summary
1. `sbatch demux_mirna.sh`
2. `sbatch make_ref.sh` (only need to run once to generate references!)
3. `sbatch trim_map.sh GRCh38` (e.g. for human samples)
4. Analysis in R and/or python.
