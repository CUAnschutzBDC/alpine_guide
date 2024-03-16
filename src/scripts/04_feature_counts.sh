#!/bin/bash 

#SBATCH --job-name=featureCounts
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --mem=1gb
#SBATCH --output=logs/featureCounts%J.out
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kristen.wells-wrasman@cuanschutz.edu

mkdir -p logs

# Activate conda environment
conda activate infer_experiment

# Update with the path to your fastq files
gtf_file="/pl/active/Anschutz_BDC/resources/ref/annotation/mouse/GRCm38/gtf/Mus_musculus.GRCm38.96.gtf"
star_dir="star"
featureCount_dir="featureCounts"
mkdir -p $featureCount_dir

featureCounts \
    --extraAttributes 'gene_name,gene_biotype' \
    -s 2 -p -B \
    -a ${gtf_file} \
    -o ${featureCount_dir}/s1_countsOutput \
    ${star_dir}/Aligned.sorted.out.bam
