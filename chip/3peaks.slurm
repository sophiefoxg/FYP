#!/bin/bash
#SBATCH --job-name=peak
#SBATCH --output=peaks.out
#SBATCH --error=peaks.err
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

aligndir=${workdir}/aligned
peakdir=${workdir}/peaks
mkdir -p ${peakdir}

module load py-macs2/2.2.8-sca46fg
 

macs2 callpeak -t ${aligndir}/GSM851278_H3K27ac_mESC.sorted.bam \
	-c ${aligndir}/GSM723020_Input_mESC.sorted.bam \
	-f BAM -g 1.3e+8 \
	-n ESC_Acetylation_Marks \
	--outdir ${peakdir}

module unload py-macs2/2.2.8-sca46fg
