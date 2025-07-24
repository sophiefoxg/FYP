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
for f in ${rawdir}/*_1.fastq.gz
do
 base=$(basename $f | sed 's/_1\.fastq\.gz//')
  fastqc -t 2 ${rawdir}/${base}_1.fastq.gz -o ${rawqcdir}
   fastqc -t 2 ${rawdir}/${base}_2.fastq.gz -o ${rawqcdir}
done

# Data trimming using fastp
for f in ${rawdir}/*_1.fastq.gz
do
 base=$(basename $f | sed 's/_1\.fastq\.gz//')
  fastp -q 20 -u 10 --cut_right \
   -i ${rawdir}/${base}_1.fastq.gz \
    -I ${rawdir}/${base}_2.fastq.gz \
     -o ${trimdir}/${base}_trim_R1.fq.gz \
     -O ${trimdir}/${base}_trim_R2.fq.gz \
        -j ${trimdir}/${base}_trim.json \
         -h ${trimdir}/${base}_trim.html
done

# Quality control for trimmed data using FastQC
for f in ${trimdir}/*_trim_R1.fq.gz
do
    base=$(basename $f | sed 's/_trim_R1\.fq\.gz//')
    fastqc -t 2 ${trimdir}/${base}_trim_R1.fq.gz ${trimdir}/${base}_trim_R2.fq.gz -o ${trimqcdir}
done

# Generate a MultiQC report
module load py-multiqc
multiqc ${trimqcdir} -o ${workdir}
module unload py-multiqc

# Unload programs
module unload fastqc
module unload fastp

echo "QC checks and trimming complete."
