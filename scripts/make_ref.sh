#!/bin/bash
#SBATCH --job-name=ref    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=16         ## number of cores the job needs
#SBATCH --output=ref-%J.out ## output log file
#SBATCH --error=ref-%J.err ## error log file

module load star

STAR \
	--runThreadN 16 \
	--runMode genomeGenerate \
	--genomeDir ../ref/mm10/ \
	--sjdbGTFfile ../ref/mm10_mirna.gtf \
	--sjdbOverhang 1 \
	--genomeFastaFiles ../ref/mm10_no_alt_analysis_set_ENCODE.fasta

gzip ../ref/mm10_no_alt_analysis_set_ENCODE.fasta

STAR \
        --runThreadN 16 \
        --runMode genomeGenerate \
        --genomeDir ../ref/hg38/ \
        --sjdbGTFfile ../ref/hg38_mirna.gtf \
        --sjdbOverhang 1 \
        --genomeFastaFiles ../ref/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

gzip ../ref/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
