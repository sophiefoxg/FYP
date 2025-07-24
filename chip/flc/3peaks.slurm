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
broadpeakdir=${workdir}/broadpeaks
mkdir -p ${broadpeakdir}
peaksdir=${workdir}/peaks
mkdir -p ${peaksdir}

module load py-macs2/2.2.8-sca46fg
 

macs2 callpeak --broad -t ${aligndir}/GSM851302_H3K27ac_embryonic-liver.sorted.bam \
	-c ${aligndir}/GSM851303_Input_embryonic-liver.sorted.bam \
	-f BAM -g 1.3e+8 \
	-n FLC_Acetylation_Marks_broad \
	--outdir ${broadpeakdir}

macs2 callpeak -t ${aligndir}/GSM851302_H3K27ac_embryonic-liver.sorted.bam \
        -c ${aligndir}/GSM851303_Input_embryonic-liver.sorted.bam \
        -f BAM -g 1.3e+8 \
        -n FLC_Acetylation_Marks \
        --outdir ${peaksdir}


module unload py-macs2/2.2.8-sca46fg
