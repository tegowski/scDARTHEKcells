#!/bin/bash
#SBATCH -a 1-1667

WORKDIR="/your/path/here"
FASTQ=$WORKDIR/rawfiles #input directory
TRIM=$WORKDIR/trim #output directory
mkdir -p $TRIM

file=$(ls $FASTQ/*_1.fastq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)  ## the ls a basename argument will change depending on how your files are named
STEM=$(basename "$file" _1.fastq.gz)

R1="$STEM"_1.fastq.gz # reassign full name to Read1
R2="$STEM"_2.fastq.gz # assign name to Read 2

flexbar -r $FASTQ/$R1 -p FASTQ/$R2 -aa Nextera --zip-output GZ --threads 16 -qf i1.8 -t $TRIM/$STEM
