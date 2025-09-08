#!/bin/bash
#SBATCH --job-name=star
#SBATCH --output=star.out
#SBATCH --error=star.err
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --mail-user=Fox-GmuerSE@cardiff.ac.uk
#SBATCH --mail-type=END,FAIL

module load star/2.7.11b-clv6cdk


# varibles
workdir=$(pwd)

rawdir=${workdir}
trimdir=${workdir}/trimdata
mkdir ${workdir}/star
stardir=${workdir}/star

genomedir=~/mydata/c2007523/FYP/genome/


for f in ${trimdir}/*_trim.fq.gz
do
    base=$(basename $f | sed 's/_trim\.fq\.gz//')


# map reads to genome
STAR   --outMultimapperOrder Random \
       --outSAMmultNmax 1 \
       --runThreadN 8 \
       --runMode alignReads \
       --outSAMtype BAM Unsorted \
       --quantMode GeneCounts \
       --readFilesCommand zcat \
       --outFileNamePrefix ${stardir}/${base}-unsort. \
       --genomeDir ${genomedir} \
       --readFilesIn ${trimdir}/${base}_trim.fq.gz 

done


module unload star/2.7.11b-clv6cdk

module load fastqc
fastqc -t 2 ${stardir}/*-unsort.Aligned.out.bam -o ${stardir}
module unload fastqc

module load py-multiqc
multiqc ${stardir}/*.zip -o ${stardir}
module unload py-multiqc

