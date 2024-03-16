#!/bin/bash 

#SBATCH --job-name=cutadapt
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --mem=1gb
#SBATCH --output=logs/cutadapt_%J.out
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kristen.wells-wrasman@cuanschutz.edu

mkdir -p logs

# Load in modules
module load cutadapt

# Update with the path to your fastq files
fastq_dir="/pl/active/Anschutz_BDC/resources/tutorials/alpine_guide/data"
fastq1="${fastq_dir}/sample1_S26_L001_R1_001.fastq.gz"
fastq2="${fastq_dir}/sample1_S26_L001_R2_001.fastq.gz"

cutadapt_dir="cutadapt"
trim_fastq1="${cutadapt_dir}/sample1_S26_L001_R1_001.fastq.gz"
trim_fastq2="${cutadapt_dir}/sample1_S26_L001_R2_001.fastq.gz"

mkdir -p $cutadapt_dir

cutadapt -m 20 \
    -a 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA' \
    -A 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'  \
    -o $trim_fastq1 -p $trim_fastq2 \
    $fastq1 $fastq2