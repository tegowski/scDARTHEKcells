#!/bin/bash
#SBATCH -a 1-2
#SBATCH -p scavenger
#SBATCH -c 8
#SBATCH --mem=60G

#This script should be run from the $WORKDIR/YTH directory
#Replace $WORKDIR with the path to your directory
WORKDIR="/work/mrt41/Ex96_SMARTseq/test"
mkdir -p $WORKDIR/aligned

#Defines list of fastq files for each batch to be aligned
filelist=$(cat $WORKDIR/software/mutfilelist${SLURM_ARRAY_TASK_ID}.txt)

#file=$(ls *_1.fastq | sed -n ${SLURM_ARRAY_TASK_ID}p)
#STEM=$(basename "$file" _1.fastq)

#Ignore these two lines for STARsolo
#file=$(ls $FASTQ/*_1.fastq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
#STEM=$(basename "$file" _1.fastq.gz)

#Run STARsolo
STAR --genomeDir $WORKDIR/STARindex/ \
--runThreadN 8 \
--readFilesIn $filelist \
--outFileNamePrefix  $WORKDIR/aligned/MUT${SLURM_ARRAY_TASK_ID} \
--outSAMtype BAM Unsorted \
--outFilterMismatchNoverReadLmax 0.06 \
--outSAMattrRGline ID: \
--soloType SmartSeq \
--readFilesManifest $WORKDIR/software/mutmanifest${SLURM_ARRAY_TASK_ID}.txt \
--soloUMIdedup Exact \
--soloStrand Unstranded

