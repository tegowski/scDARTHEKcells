#!/bin/bash
#SBATCH -a 1-832

file=$(ls *1.fastq | sed -n ${SLURM_ARRAY_TASK_ID}p)
STEM=$(basename "$file" _1.fastq)

mv ${STEM}_1.fastq WTCell${SLURM_ARRAY_TASK_ID}_1.fastq 
mv ${STEM}_2.fastq WTCell${SLURM_ARRAY_TASK_ID}_2.fastq
