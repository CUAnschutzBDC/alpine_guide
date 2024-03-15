#!/bin/bash 

#SBATCH --job-name=test_job
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=1gb
#SBATCH --output=logs/test_%J.out
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kristen.wells-wrasman@cuanschutz.edu

mkdir -p logs

echo "This is the first step"
sleep 60
echo "Finished"