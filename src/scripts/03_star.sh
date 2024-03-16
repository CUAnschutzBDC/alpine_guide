#!/bin/bash 

#SBATCH --job-name=star
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --mem=50gb
#SBATCH --output=logs/star_%J.out
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kristen.wells-wrasman@cuanschutz.edu

mkdir -p logs

# Load in modules
module load star/2.7.10b
module load samtools

# Update with the path to your fastq files
cutadapt_dir="cutadapt"
trim_fastq1="${cutadapt_dir}/sample1_S26_L001_R1_001.fastq.gz"
trim_fastq2="${cutadapt_dir}/sample1_S26_L001_R2_001.fastq.gz"
genome="/pl/active/Anschutz_BDC/resources/ref/indices/star/mouse/GRCm38"
gtf_file="/pl/active/Anschutz_BDC/resources/ref/annotation/mouse/GRCm38/gtf/Mus_musculus.GRCm38.96.gtf"
star_dir="star"
mkdir -p $star_dir

STAR \
    --runThreadN 1 \
    --genomeDir $genome \
    --sjdbGTFfile $gtf_file \
    --readFilesIn $trim_fastq1 $trim_fastq2 \
    --readFilesCommand zcat \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix $star_dir/

samtools sort ${star_dir}/Aligned.out.bam -T /scratch/alpine/$USER > ${star_dir}/Aligned.sorted.out.bam
rm ${star_dir}/Aligned.out.bam
samtools index ${star_dir}/Aligned.sorted.out.bam