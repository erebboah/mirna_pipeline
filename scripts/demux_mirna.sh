#!/bin/bash
#SBATCH --job-name=mirna_demux    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=8         ## number of cores the job needs
#SBATCH --output=mirna_demux-%J.out ## output log file
#SBATCH --error=mirna_demux-%J.err ## error log file

source ~/miniconda3/bin/activate seqtkpython3 # replace with your conda env name containing seqtk and python3

mkdir ../fastq/temp
module load cutadapt
module load ucsc-tools

##### Figure out which reads correspond to our samples

python3 five_prime_demux.py ../fastq/undetermined.fastq ../fastq/temp

samples=(1 2 3 4 5 6 7)

for sample in "${samples[@]}"; do
    cutadapt -a ATCTCGTATGCCGTCTTCTGCTT -e 0.25 -m 15 --match-read-wildcards --untrimmed-output=../fastq/temp/${sample}_NO3AD.fastq ../fastq/temp/${sample}.fastq > ../fastq/temp/${sample}_trimmed.fastq
    rm ../fastq/temp/${sample}_NO3AD.fastq 
    rm ../fastq/temp/${sample}.fastq

    mkdir -p ../fastq/temp/${sample}
    python3 three_prime_demux.py ../fastq/temp/${sample}_trimmed.fastq ../fastq/temp/${sample}

    # Rename files in the sample directory
    for f in ../fastq/temp/${sample}/*; do
        new_name="${sample}_$(basename "$f")"
        mv -- "$f" "../fastq/temp/${new_name}"
    done
done

##### Actually demultiplex original fastq

for file in `awk -F "," '{print $1}' samplesheet.csv`; do
	awk 'NR%4==1' ../fastq/temp/${file}.fastq | awk -F " " '{print $1}' | sed 's/@//' > ../fastq/temp/${file}.names
	seqtk subseq undetermined.fastq ../fastq/temp/${file}.names > ../fastq/${file}.fastq
	
	# print number of reads 
	awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' ../fastq/${file}.fastq
	
	# validate format using ucsc-tools
	validateFiles -chromInfo=../ref/mm10.chrom.sizes -type=fastq ../fastq/${file}.fastq
done

## rename fastqs according to samplesheet.csv
sed 's|^|mv -vi "../fastq/|; s|, |.fastq" "../fastq/|; s/$/.fastq";/' < samplesheet.csv | bash -

## gzip final files
gzip ../fastq/ENC* ##### TEMPORARY ###########

## remove temp files
rm -r ../fastq/temp
