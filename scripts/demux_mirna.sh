#!/bin/bash
#SBATCH --job-name=mirna_demux    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p free               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=8         ## number of cores the job needs
#SBATCH --output=mirna_demux-%J.out ## output log file
#SBATCH --error=mirna_demux-%J.err ## error log file

python3 dual_demux.py ../fastq/undetermined.fastq.gz ../fastq
