# MicroRNA-seq Analysis
## Overview
MicroRNA-seq libraries are built using the [ENCODE protocol](https://www.encodeproject.org/documents/49f43842-5ab4-4aa1-a6f4-2b1234955d93/@@download/attachment/Multiplexed%20miRNA%20Sequencing%20Library%20Generation%20Protocol-3.4.pdf) "Multiplexed miRNA Sequencing Library Generation Protocol".

Briefly, libraries are built from 400 ng of total RNA. Adapters are ligated to the 5’ and 3’ ends of small RNA using its 5’ phosphate and 3’ hydroxyl groups, then the ligation product is reverse transcribed. The 5’ adapter adds a 6-nucleotide barcode from a one of 7 sets of 4 distinct barcodes used in downstream demultiplexing. The cDNA is amplified with 58 bp reverse and 55 bp forward primers containing additional 6-nucleotide barcodes added to the 3’ end. The 140 bp product containing mature microRNA (21-25 nucleotides) is size-selected using 10% TBE-Urea gel. 

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

3. Edit the "real" sample sheet. The first column is the 5' Adaptor ID (1-7) and forward primer index ID (1-12, 14, 15, 18, 21-25) separated by an underscore. The second column is your desired fastq name (sample ID). No header! Preferably matches some sample ID in your metadata. Should be in the same folder as the scripts.

[Example of samplesheet.csv](https://github.com/erebboah/mirna_pipeline/blob/master/scripts/samplesheet.csv)

Edit `[demux_mirna.sh](https://github.com/erebboah/mirna_pipeline/blob/master/scripts/demux_mirna.sh)` so that you are using your HPC account (change `#SBATCH -A SEYEDAM_LAB`) and your conda environment (change `source ~/miniconda3/bin/activate seqtkpython3`). As inputs, you need fastq named `undetermined.fastq` in the fastq directory of this repo (e.g. `/pub/erebboah/mirna_pipeline/fastq/`) and your [samplesheet.csv](https://github.com/erebboah/mirna_pipeline/blob/master/scripts/samplesheet.csv) in the scripts directory (e.g. `/pub/erebboah/mirna_pipeline/scripts/`).

## Quantification
1. Make STAR reference. Download microRNA GENCODE GTFs from ENCODE portal: [vM21 mouse](https://www.encodeproject.org/files/ENCFF094ICJ/) or [v29 human](https://www.encodeproject.org/files/ENCFF470CZH/).

```
cd ../ref
wget https://www.encodeproject.org/files/ENCFF094ICJ/@@download/ENCFF094ICJ.gtf.gz
wget https://www.encodeproject.org/files/ENCFF470CZH/@@download/ENCFF470CZH.gtf.gz
mv ENCFF094ICJ.gtf.gz mm10_vM21_mirna.gtf
mv ENCFF470CZH.gtf.gz hg38_v29_mirna.gtf
```

Download reference sequences from the portal: [mm10 mouse](https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/) or [hg38 human](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/).

```
wget https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
```

Once you have the reference files, run STAR in genomeGenerate mode

## Downstream analysis
