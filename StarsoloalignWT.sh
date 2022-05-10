#!/bin/bash
#SBATCH -a 1-2
#SBATCH -c 8
#SBATCH --mem=60G

#This script should be run from the $WORKDIR/YTH directory
#Replace $WORKDIR with the path to your directory
WORKDIR="/your/path/here"
mkdir -p $WORKDIR/aligned

#Defines list of fastq files for each batch to be aligned
filelist=$(cat $WORKDIR/software/WTfilelist${SLURM_ARRAY_TASK_ID}.txt)

#file=$(ls *_1.fastq | sed -n ${SLURM_ARRAY_TASK_ID}p)
#STEM=$(basename "$file" _1.fastq)

#Ignore these two lines for STARsolo
#file=$(ls $FASTQ/*_1.fastq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
#STEM=$(basename "$file" _1.fastq.gz)

#Run STARsolo
STAR --genomeDir $WORKDIR/STARindex/ \
--runThreadN 8 \
--readFilesIn $filelist \
--outFileNamePrefix  $WORKDIR/aligned/WT${SLURM_ARRAY_TASK_ID} \
--outSAMtype BAM Unsorted \
--outFilterMismatchNoverReadLmax 0.06 \
--outSAMattrRGline ID: \
--soloType SmartSeq \
--readFilesManifest $WORKDIR/software/WTmanifest${SLURM_ARRAY_TASK_ID}.txt \
--soloUMIdedup Exact \
--soloStrand Unstranded

