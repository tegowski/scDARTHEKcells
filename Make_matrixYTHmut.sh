#!/bin/bash
#SBATCH -a 1-1254

BAMS=$WORKDIR/singlebams
file=$APOYTHmutmerge.bam
STEM=$(basename "$file" .bam)
mkdir $WORKDIR/matrix

perl $WORKDIR/software/parseBAM.pl -i $file -o $WORKDIR/matrix/$STEM.matrix --minCoverage 3 --verbose --removeDuplicates