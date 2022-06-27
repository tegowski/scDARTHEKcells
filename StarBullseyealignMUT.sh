#!/bin/bash
#SBATCH --mem=60G
#SBATCH -c 12
#SBATCH -p scavenger

#Merge all MUT bamfiles
samtools merge MUTCellsmerged.bam MUT*duprm.bam

#Then use Bullseye (parseBAM.pl) to make nucleotide matrix
WORKDIR="/your/path/here/"
BAMS=$WORKDIR/singlebams
STEM=$(basename "$file" Aligned.sortedByCoord.out.bam.duprm.bam)
mkdir $WORKDIR/matrix

perl $WORKDIR/software/parseBAM.pl \
-i MUTCellsmerged.bam \
-o $WORKDIR/matrix/MUTCellsmerged.matrix \
--minCoverage 5 \
--verbose \
--removeDuplicates
