#!/bin/bash
#SBATCH -a 1-5

#This script should be run from the $WORKDIR/YTH directory
#Replace $WORKDIR with the path to your directory
$WORKDIR="/path/to/directory"

#Defines list of fastq files for each batch to be aligned
filelist=$(cat WTfilelist${SLURM_ARRAY_TASK_ID}.txt)

file=$(ls $FASTQ/*_1.fastq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
STEM=$(basename "$file" _1.fastq.gz)

#Ignore these two lines for STARsolo
#file=$(ls $FASTQ/*_1.fastq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
#STEM=$(basename "$file" _1.fastq.gz)

#Run STARsolo
STAR --genomeDir $WORKDIR/STARindex/ \
--runThreadN 8 \
--readFilesIn $filelist \
--readFilesCommand zcat \
--outFileNamePrefix  $WORKDIR/aligned/WT${SLURM_ARRAY_TASK_ID} \
--outSAMtype BAM Unsorted \
--outFilterMismatchNoverReadLmax 0.06 \
--outSAMattrRGline ID: \
--soloType SmartSeq \
--readFilesManifest $WORKDIR/software/WTmanifest${SLURM_ARRAY_TASK_ID}.txt \
--soloUMIdedup Exact \
--soloStrand Unstranded

