#!/bin/bash
#SBATCH --job-name=bowtie
#SBATCH --output=bowtie.out
#SBATCH --error=bowtie.err
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --mail-user=Fox-GmuerSE@cardiff.ac.uk
#SBATCH --mail-type=END,FAIL

module load bowtie2/2.5.2-hg77ifp

workdir=$(pwd)
trimdir=${workdir}/trimdata
mkdir ${workdir}/aligned
aligndir=${workdir}/aligned
genomedir=~/mydata/c2007523/FYP/genome

# align with bowtie 2

for f in ${trimdir}/*_trim.fastq.gz
do
    base=$(basename $f _trim.fastq.gz)

 bowtie2 -p 8 -q --local \
-x ${genomedir}/mm9_index \
-U ${trimdir}/${base}_trim.fastq.gz \
-S ${aligndir}/${base}_trim.sam

done

module unload bowtie2/2.5.2-hg77ifp

#samtools sam to bam

module load samtools/1.19.2-i77kweo

for f in ${aligndir}/*_trim.sam
do
    base=$(basename $f _trim.sam)

    # Convert SAM to BAM and sort
    samtools view -@ 2 -bS ${aligndir}/${base}_trim.sam | \
    samtools sort -@ 2 -o ${aligndir}/${base}.sorted.bam

    # Index the sorted BAM
    samtools index ${aligndir}/${base}.sorted.bam
done

module unload samtools/1.19.2-i77kweo

