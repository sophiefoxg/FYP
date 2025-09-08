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

#Define variables
workdir=$(pwd)
rawdir=${workdir}
rawqcdir=${workdir}/rawqc
trimdir=${workdir}/trimdata

mkdir -p ${rawqcdir} ${trimdir}

# Load tools
module load fastqc
module load fastp
module load py-multiqc
# Run FastQC on raw data
for f in ${rawdir}/*.fastq.gz
do
    fastqc -t 4 $f -o ${qcdir}
done

# Optional trimming with fastp (for single-end)
for f in ${rawdir}/*.fastq.gz
do
    base=$(basename $f .fastq.gz)
    fastp -i $f \
          -o ${trimdir}/${base}_trim.fastq.gz \
          -j ${trimdir}/${base}_trim.json \
          -h ${trimdir}/${base}_trim.html
done
# Run FastQC on trimmed data
for f in ${trimdir}/*_trim.fastq.gz
do
    fastqc -t 4 $f -o ${qcdir}
done

# Generate MultiQC report
multiqc ${qcdir} -o ${workdir}

# Unload modules
module unload fastqc
module unload fastp
module unload py-multiqc
