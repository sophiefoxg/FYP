#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --output=fastqcp.out
#SBATCH --error=fastqcp.err
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
rawdir=${workdir}

# Create directories if they don't exist


rawqcdir=${workdir}/rawqc
trimdir=${workdir}/trimdata
trimqcdir=${workdir}/trimqc

mkdir -p ${rawqcdir} ${trimdir} ${trimqcdir}
# Load required programs
module load fastqc
module load fastp

# Quality control for raw data using FastQC
for f in ${rawdir}/*.fastq.gz
do
 base=$(basename $f | sed 's/\.fastq\.gz//')
  fastqc -t 2 ${f} -o ${rawqcdir}
done

# Data trimming using fastp
for f in ${rawdir}/*.fastq.gz
do
 base=$(basename $f | sed 's/\.fastq\.gz//')
  fastp -q 20 -u 10 --cut_right \
   -i ${rawdir}/${base}.fastq.gz \
    -o ${trimdir}/${base}_trim.fq.gz \
        -j ${trimdir}/${base}_trim.json \
         -h ${trimdir}/${base}_trim.html
done

# Quality control for trimmed data using FastQC
for f in ${trimdir}/*_trim.fq.gz
do
    base=$(basename $f | sed 's/_trim\.fq\.gz//')
    fastqc -t 2 ${trimdir}/${base}_trim.fq.gz  -o ${trimqcdir}
done

# Generate a MultiQC report
module load py-multiqc
multiqc ${trimqcdir} -o ${workdir}
module unload py-multiqc

# Unload programs
module unload fastqc
module unload fastp

echo "QC checks and trimming complete."
