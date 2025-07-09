#!/bin/bash
#SBATCH --job-name=fc
#SBATCH --output=fc.out
#SBATCH --error=fc.err
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --mail-user=Fox-GmuerSE@cardiff.ac.uk
#SBATCH --mail-type=END,FAIL

# Write jobscript to output file (good for reproducibility)
cat $0

# Load some modules
module load subread/2.0.6-qefcod7

##Directories
workdir=$(pwd)

rawdir=${workdir}
trimdir=${workdir}/trimdata
stardir=${workdir}/star
genomedir=~/mydata/c2007523/FYP/genome
markdir=${workdir}/markdup
mkdir ${workdir}/featureCounts
fcdir=${workdir}/featureCounts

for f in ${workdir}/*.fastq.gz
do
R1=$(basename $f | cut -f1 -d.)
base=$(basename $f | sed 's/.fastq\.gz//')

featureCounts \
        -T 4 -F GTF -t exon -g gene_id \
	-a ${genomedir}/mm9.ensGene.gtf  \
	-o ${fcdir}/${base}.featurecount.txt \
	${stardir}/${base}*.bam


done
 module unload subread/2.0.6-qefcod7
