#!/bin/bash
#SBATCH -p scavenger
#SBATCH -c 8
#SBATCH --mem=60G
#SBATCH -a 1-301


file=$(ls WT*Aligned.sortedByCoord.out.bam | sed -n ${SLURM_ARRAY_TASK_ID}p)
STEM=$(basename "$file" .Aligned.sortedByCoord.out.bam)

echo "$[SLURM_ARRAY_TASK_ID]"
echo "Starting File: $file"

#First, use samtools to remove optical duplicates
samtools sort -n -o $STEM.name.bam $file
samtools fixmate -r -m $STEM.name.bam $STEM.fix.bam
samtools sort -o $STEM.possort.bam $STEM.fix.bam
samtools markdup -r $STEM.possort.bam $STEM.duprm.bam

rm $STEM.name.bam $STEM.fix.bam $STEM.possort.bam

#Then use Bullseye (parseBAM.pl) to make nucleotide matrix
WORKDIR="/your/path/here"
BAMS=$WORKDIR/singlebams
file=$(ls WT*duprm.bam | sed -n ${SLURM_ARRAY_TASK_ID}p)
STEM=$(basename "$file" .duprm.bam)
mkdir $WORKDIR/matrix

perl $WORKDIR/software/parseBAM.pl \
-i $file \
-o $WORKDIR/matrix/$STEM.matrix \
--minCoverage 5 \
--verbose \
--removeDuplicates
