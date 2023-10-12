#!/bin/bash
#SBATCH --job-name=map    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --array=1-12              ## number of tasks to launch (number of samples)
#SBATCH --cpus-per-task=16         ## number of cores the job needs
#SBATCH --output=map-%J.out ## output log file
#SBATCH --error=map-%J.err ## error log file

inpath="/pub/erebboah/mirna_pipeline/"
genome=$1 # either hg38 or mm10

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <genome>"
    exit 1
fi

file=$inpath"scripts/samplesheet.csv"
prefix=`head -n $SLURM_ARRAY_TASK_ID  ${file} | tail -n 1 |  awk -F ", " '{print $2}'`
adapter=$(head -n $SLURM_ARRAY_TASK_ID ${file} | tail -n 1 | awk -F "," '{gsub(/_/,"",$1); print substr($1, 1, 1)}')
five_prime_adapters=$inpath"/ref/five_prime_adapter_set"$adapter".fasta"
three_prime_adapters=$inpath"/ref/three_prime_adapter.fasta"

module load cutadapt
module load star

mkdir ${inpath}${prefix}
mkdir ${inpath}${prefix}/cutadapt

cutadapt \
	-a file:${three_prime_adapters} \
	-e 0.25 \
	--match-read-wildcards \
	--untrimmed-output=${inpath}${prefix}/cutadapt/${prefix}_NO3AD.fastq \
	${inpath}temp/${prefix}.fastq.gz \
	| cutadapt \
	-e 0.34 \
	--match-read-wildcards \
	--no-indels \
	-m 15 \
	-O 6 \
	-n 1 \
	-g file:${five_prime_adapters} \
	--untrimmed-output=${inpath}${prefix}/cutadapt/${prefix}_NO5AD.fastq \
	--too-short-output=${inpath}${prefix}/cutadapt/${prefix}_SHORT_FAIL.fastq \
	- > ${inpath}${prefix}/cutadapt/${prefix}_trim.fastq

mkdir ${inpath}${prefix}/star
mkdir ${inpath}/counts

STAR \
	--genomeDir ${inpath}ref/${genome} \
	--readFilesIn ${inpath}${prefix}/cutadapt/${prefix}_trim.fastq \
	--sjdbGTFfile ${inpath}ref/${genome}_mirna.gtf \
	--runThreadN 8 \
	--alignEndsType EndToEnd \
	--outFilterMismatchNmax 1 \
	--outFilterMultimapScoreRange 0 \
	--quantMode TranscriptomeSAM GeneCounts \
	--outReadsUnmapped Fastx \
	--outSAMtype BAM SortedByCoordinate \
	--outFilterMultimapNmax 10 \
	--outSAMunmapped Within \
	--outFilterScoreMinOverLread 0 \
	--outFilterMatchNminOverLread 0 \
	--outFilterMatchNmin 16 \
	--alignSJDBoverhangMin 1000 \
	--alignIntronMax 1 \
	--outWigType wiggle \
	--outWigStrand Stranded \
	--outWigNorm RPM \
	--outFileNamePrefix ${inpath}${prefix}/star/

mv ${inpath}${prefix}/star/ReadsPerGene.out.tab ${inpath}/counts/${prefix}.tsv


