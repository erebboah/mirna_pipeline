#!/bin/bash
#SBATCH --job-name=ref    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p free               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=16         ## number of cores the job needs
#SBATCH --output=ref-%J.out ## output log file
#SBATCH --error=ref-%J.err ## error log file

module load star

STAR \
	--runThreadN 16 \
	--runMode genomeGenerate \
	--genomeDir ../ref/hg38/ \
	--sjdbGTFfile ../ref/hg38_mirna.gtf \
	--sjdbOverhang 1 \
	--genomeFastaFiles ../ref/hg38.fa

gzip ../ref/hg38.fa

