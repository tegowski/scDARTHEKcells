#!/bin/bash
#SBATCH -a 1-1254

BAMS=$WORKDIR/singlebams
file=$(ls WT*dupMarked.bam | sed -n ${SLURM_ARRAY_TASK_ID}p)
STEM=$(basename "$file" .dupmarked.bam)
mkdir $WORKDIR/matrix

perl $WORKDIR/software/parseBAM.pl -i $file -o $WORKDIR/matrix/$STEM.matrix --minCoverage 3 --verbose --removeDuplicates