#!/bin/bash
#SBATCH -a 1-601
#SBATCH -c 8
#SBATCH --mem=60G
#SBATCH -p scavenger

#Define path to your working directory
WORKDIR="/work/mrt41/Ex96_SMARTseq/test"

file=$(ls *_1.fastq | sed -n ${SLURM_ARRAY_TASK_ID}p) ## this is done for paired end libraries ending with _1 and _2 for each pair.
STEM=$(basename "$file" _1.fastq)

BAM=$WORKDIR/singlebams
mkdir -p $BAM

echo "aligning $STEM"

STAR --runMode alignReads \
--genomeDir $WORKDIR/STARindex \
--runThreadN 8 \
--readFilesIn ${STEM}_1.fastq ${STEM}_2.fastq \
--outSAMtype BAM SortedByCoordinate \
--outFilterMismatchNoverReadLmax 0.06 \
--outFileNamePrefix $BAM/${STEM} \

#This next steps marks PCR duplicates in the BAM file so they can be ignored during the next step
STAR --runMode inputAlignmentsFromBAM \
--runThreadN 4 \
--inputBAMfile $BAM/${STEM}.Aligned.sortedByCoord.out.bam \
--bamRemoveDuplicatesType UniqueIdentical \
--outFileNamePrefix $BAM/${STEM}.dupMarked.bam
