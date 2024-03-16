#!/bin/bash 

#SBATCH --job-name=fastqc
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --mem=1gb
#SBATCH --output=logs/fastqc_%J.out
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kristen.wells-wrasman@cuanschutz.edu

mkdir -p logs

# Load in modules
module load fastqc

# Update with the path to your fastq files
fastq_dir="/pl/active/Anschutz_BDC/resources/tutorials/alpine_guide/data"
fastq1="${fastq_dir}/sample1_S26_L001_R1_001.fastq.gz"
fastq2="${fastq_dir}/sample1_S26_L001_R2_001.fastq.gz"
outdir="fastqc"
mkdir -p $outdir

fastqc $fastq1 $fastq2 --outdir $outdir