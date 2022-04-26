#!/bin/bash
#SBATCH --mem=60G
#SBATCH -c 4

#Please define $WORKDIR to the path of the parent directory for this project
$WORKDIR="/your/path/here"

STAR --runThreadN 4 \
--runMode genomeGenerate \
--genomeDir $WORKDIR/STARindex \
--genomeFastaFiles $WORKDIR/genomefasta/* \
--sjdbGTFfile $WORKDIR/gtf/Homo_sapiens.GRCh38.98.gtf \
--sjdbOverhang 49 \
--genomeSAindexNbases 3