#!/bin/bash
#SBATCH --mem=24G
#SBATCH -p scavenger
#SBATCH -c 8
#SBATCH -a 1-301

annotation_file=/datacommons/meyerlab/userdata/matt_tego/gtf/HG38/hg38refFlat_nochr.gtf
WORKDIR="/your/path/here/"
bed=$WORKDIR/bed
file=$(ls WT*matrix.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
STEM=$(basename "$file" .matrix.gz)
mkdir -p $bed

perl $WORKDIR/software/Find_edit_site.pl \
        --annotationFile $annotation_file \
        --EditedMatrix $file \
        --controlMatrix MUTCellsmerged.matrix.gz \
        --editType C2U \
        --minEdit 10 \
        --maxEdit 95 \
        --editFoldThreshold 1.5 \
        --EditedMinCoverage 20 \
        --ControlMinCoverage 10 \
        --MinEditSites 2 \
        --cpu 8 \
        --outfile $bed/$STEM.bed \
        --bed6 \
        --verbose
