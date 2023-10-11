# MicroRNA-seq Analysis
## Overview
MicroRNA-seq libraries are built using the [ENCODE protocol](https://www.encodeproject.org/documents/49f43842-5ab4-4aa1-a6f4-2b1234955d93/@@download/attachment/Multiplexed%20miRNA%20Sequencing%20Library%20Generation%20Protocol-3.4.pdf) "Multiplexed miRNA Sequencing Library Generation Protocol".

Briefly, libraries are built from 400 ng of total RNA. Adapters are ligated to the 5’ and 3’ ends of small RNA using its 5’ phosphate and 3’ hydroxyl groups, then the ligation product is reverse transcribed using SuperScript II Reverse Transcriptase (Invitrogen cat. #18064-071). The 5’ adapter adds a 6-nucleotide barcode from a one of 7 sets of 4 distinct barcodes used in downstream demultiplexing. The cDNA is amplified using Phusion PCR master mix (NEB cat. #M0531S) with 58 bp reverse and 55 bp forward primers containing additional 6-nucleotide barcodes added to the 3’ end. The 140 bp product containing mature microRNA (21-25 nucleotides) is size-selected using 10% TBE-Urea gel (BioRad cat. #456-6033). 

Libraries are isolated from the TBE-Urea gel by agitated incubation at 70°C, 1000 RPM for 2 hours in a buffer containing 0.5 M ammonium acetate (Ambion cat. #AM9070G), 0.1% SDS (Sigma cat. #L6026-50G), and 0.1 uM EDTA (Ambion cat. #AM9261), then precipitated overnight in 50% isopropanol. 

Resulting microRNA-seq libraries are quantified using Qubit dsDNA HS Assay Kit (Thermo cat. #Q32854) and sequenced on an Illumina NextSeq 2000 with P2 100 cycle kits (Illumina cat. #20046811) as 75 bp single-end reads to around 10 M raw read depth per library. Submission to the ENCODE portal requires >5 M aligned reads, >300 microRNAs detected at >2 CPM, and a Spearman replicate correlation >0.85.

<img src="https://github.com/erebboah/mirna_pipeline/blob/master/mirna_overview.png" width="602" height="502">


## Demultiplexing
Library structure deviates from standard Illumina, so samples must be demultiplexed with custom code. Illumina NextSeq2000 SampleSheet (e.g. for our [practice run](https://github.com/erebboah/mirna_pipeline/blob/master/ref/SampleSheet_1080.csv)) needs to specify read setup but the SampleID / Index doesn't matter since we're just going to use the "undetermined" reads.

1. After sequencing is complete, transfer reads to HPC.

```
scp /home/sharing/runs/231006_VH00582_82_AACYV55M5/Analysis/1/Data/fastq/*.fastq.gz erebboah@hpc3.rcic.uci.edu:/pub/erebboah/mirna_pipeline/fastq/run_1080
```

2. Cat all reads to undetermined.fastq.gz and unzip, making sure you have enough space. Need 24G of space for a P1 100 cycle kit.

```
cd /pub/erebboah/mirna_pipeline/fastq/run_1080
cat *.fastq.gz > undetermined.fastq.gz
gunzip undetermined.fastq.gz
```

## Quantification


## Downstream analysis
