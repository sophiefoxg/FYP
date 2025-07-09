#!/bin/bash
#SBATCH --job-name=rnaseq
#SBATCH --output=rnaseq.out
#SBATCH --error=rnaseq.err
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --mail-user=Fox-GmuerSE@cardiff.ac.uk
#SBATCH --mail-type=END,FAIL

# Define variables
workdir=$(pwd)


sbatch 1fastQCP.slurm

