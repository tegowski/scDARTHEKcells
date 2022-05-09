#!/bin/bash
#SBATCH -a 1-1667

file=$(ls *_1.fastq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p) ## this is done for paired end libraries ending with _1 and _2 for each pair.
STEM=$(basename "$file" _1.fastq.gz)
BAM=$WORKDIR/singlebams
mkdir -p $BAM

echo "aligning $STEM"

STAR --runMode alignReads \
--genomeDir $WORKDIR/Starindex \
--runThreadN 8 \
--readFilesIn $dir/${STEM}_1.fastq.gz $dir/${STEM}_2.fastq.gz \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outFilterMismatchNoverReadLmax 0.06 \
--outFileNamePrefix $BAM/${STEM}. \

#This next steps marks PCR duplicates in the BAM file so they can be ignored during the next step
STAR --runMode inputAlignmentsFromBAM \
--runThreadN 4 \
--inputBAMfile $BAM/${STEM}.Aligned.sortedByCoord.out.bam \
--bamRemoveDuplicatesType UniqueIdentical \
--outFileNamePrefix $BAM/${STEM}.dupMarked.bam