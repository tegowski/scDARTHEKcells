#!/bin/bash
#SBATCH -a 1-2

#Replace $WORKDIR with the path to your directory
$WORKDIR="/path/to/directory"

#Defines list of fastq files for each batch to be aligned
filelist=$(cat mutfilelist${SLURM_ARRAY_TASK_ID}.txt)

#Ignore these two lines for STARsolo
#file=$(ls $FASTQ/*_1.fastq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
#STEM=$(basename "$file" _1.fastq.gz)

#Run STARsolo
STAR --genomeDir $WORKDIR/STARindex/ \
--runThreadN 8 \
--readFilesIn $filelist \
--readFilesCommand zcat \
--outFileNamePrefix  $WORKDIR/aligned/mut${SLURM_ARRAY_TASK_ID} \
--outSAMtype BAM Unsorted \
--outFilterMismatchNoverReadLmax 0.06 \
--outSAMattrRGline ID: \
--soloType SmartSeq \
--readFilesManifest $WORKDIR/software/mutmanifest${SLURM_ARRAY_TASK_ID}.txt \
--soloUMIdedup Exact \
--soloStrand Unstranded