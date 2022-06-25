  GNU nano 2.9.8                                                                             duprm.sh                                                                              Modified  
#!/bin/bash
#SBATCH -p scavenger
#SBATCH -c 8
#SBATCH --mem=60G
#SBATCH -a 1-300


file=$(ls MUT*Aligned.sortedByCoord.out.bam | sed -n ${SLURM_ARRAY_TASK_ID}p)
STEM=$(basename "$file" .Aligned.sortedByCoord.out.bam)

echo "$[SLURM_ARRAY_TASK_ID]"
echo "Starting File: $file"

#First, use samtools to remove optical duplicates
samtools sort -n -o $STEM.name.bam $file
samtools fixmate -r -m $STEM.name.bam $STEM.fix.bam
samtools sort -o $STEM.possort.bam $STEM.fix.bam
samtools markdup -r $STEM.possort.bam $STEM.duprm.bam

rm $STEM.name.bam $STEM.fix.bam $STEM.possort.bam
