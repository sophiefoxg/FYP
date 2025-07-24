#!/bin/bash
#SBATCH --job-name=pipeline
#SBATCH --output=pipeline.out
#SBATCH --error=pipeline.err
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --mail-user=Fox-GmuerSE@cardiff.ac.uk
#SBATCH --mail-type=END,FAIL

# Define variables
workdir=$(pwd)

# Submit first job
jid1=$(sbatch 1fastQCP.slurm | awk '{print $4}')

# Submit second job dependent on first
jid2=$(sbatch --dependency=afterok:$jid1 2star.slurm | awk '{print $4}')

# Submit third job dependent on second
sbatch --dependency=afterok:$jid2 3featurecounts.slurm

