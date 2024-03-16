#!/bin/bash 

#SBATCH --job-name=featureCounts
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --mem=1gb
#SBATCH --output=logs/featureCounts%J.out
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kristen.wells-wrasman@cuanschutz.edu

python countTable.py \
    -f featureCounts/s1_countsOutput