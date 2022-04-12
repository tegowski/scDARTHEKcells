# scDART in HEK293T cells

## **Overview**
This contains information necessary for reproducing the identification of m<sup>6</sup>A 

**Prerequisites**
1) Download raw fastq files
2) Generate manifest file
3) Download flexbar (v3.0.3) and STAR (v2.7.5a), preferably into conda environment
4) Download Seurat (v3.2 or later)
5) Download genome and annotation files
6) Create STAR index for genome alignment


***1) Download raw fastq files***
Raw sequencing files are stored at SRA under the accession number GSE180954.
Files with SRA run numbers between SRR15268425-SRR15268461 and SRR15268889-SRR15270105 are APOBEC1-YTH-expressing, while those between SRR15268461-SRR15268876 are APOBEC1-YTHmut-expressing.
```bash
for (( i = 68425; i <= 68461; i++ ))
  do
fastq-dump --outdir $WORKDIR/YTH --split-files SRR152$i
  done

for (( i = 68889; i <= 70105; i++ ))
  do
fastq-dump --outdir $WORKDIR/YTH --split-files SRR152$i
  done
  
for (( i = 68461; i <= 68876; i++ ))
  do
fastq-dump --outdir $WORKDIR/YTHmut --split-files SRR152$i
  done
```
Rename all the APOBEC1-YTH-expressing files to include a cell number and transgene expression status
```bash

```
Rename all the APOBEC1-YTHmut-expressing files to include a cell number and transgene expression status
```bash

```
Then put them all in the same rawfiles directory
```bash

```
***2) Make manifest file***
The manifest file is required for using STARsolo. This file is tab separated and contains 3 columns. 1) read1-file 2) read2-file 3) read-group.

```bash
ls *WT

```


***3) Download flexbar and STAR***
If using conda environment.
```bash
conda install -c bioconda flexbar=3.0.3
conda install -c bioconda STAR=2.7.5a
```

***4) Download Seurat***
In R 3.6 or later.
```R
install.packages("Seurat")
```
For more installation instructions for Seurat, see https://satijalab.org/seurat/articles/install.html
For more detailed tutorials for Seurat, see https://satijalab.org/seurat/articles/get_started.html

***5) Download genome and annotation files***
To download genome files:
```bash
#Navigate to $WORKDIR/genomefasta
for (( i = 1; i <= 2; i++ ))
  do
wget http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.$i.fa.gz
  done
wget http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.X.fa.gz
wget http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.Y.fa.gz
wget http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.MT.fa.gz
wget http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.nonchromosomal.fa.gz
```
To download gtf annotation file:
```bash
#Navigate to $WORKDIR/gtf
wget Homo_sapiens.GRCh38.98.gtf.gz
```

***6) Create STAR index***
```bash
STAR --runThreadN 4 \
--runMode genomeGenerate \
--genomeDir $WORKDIR/STARindex \
--genomeFastaFiles $WORKDIR/genomefasta/* \
--sjdbGTFfile $WORKDIR/gtf/Homo_sapiens.GRCh38.98.gtf \
--sjdbOverhang 49 \
--genomeSAindexNbases 3
```

## **Data Processing**
###**Adapter trimming and genome alignment on Linux/Unix system**
The first two steps after obtaining fastq files occur on a Linux/Unix system
1) Trimming of adapter sequences
2) Genome alignment and feature counting

***1) Trimming of adapter sequences***
Flexbar (3.0.3) is used for trimming adapter sequences on reads. Since Nextera-compatible indexing primers were used, we can use the -aa Nextera option. If analyzing other adapter sequences, you may need to specify your adapter sequences in a fasta file. It is suggested to submit as array jobs as there are so many files. If using a system with slurm manager, the #SBATCH -a option can be used. This should be changed if using other systems.

```bash
#!/bin/bash
#SBATCH -a 1-1667

FASTQ=$WORKDIR/rawfiles #input directory
TRIM=$WORKDIR/trim #output directory
mkdir -p $TRIM

file=$(ls $FASTQ/*_1.fastq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)  ## the ls a basename argument will change depending on how your files are named
STEM=$(basename "$file" _1.fastq.gz)

R1="$STEM"_1.fastq.gz # reassign full name to Read1
R2="$STEM"_2.fastq.gz # assign name to Read 2

flexbar -r $FASTQ/$R1 -p FASTQ/$R2 -aa Nextera --zip-output GZ --threads 16 -qf i1.8 -t $TRIM/$STEM
```

***2) Genome alignment and feature counting***
This step can take alot of time if all cells are processed in one job. There therefore it is recommended to break up the files into several "batches" that can be run simultaneously. Therefore 5 different groups of "WT" or APOBEC1-YTH-expressing cells and 3 different groups of "mut" or APOBEC1-YTHmut-expressing cells can be run at a time as separate jobs. These can be re-integrated later in Seurat.

```bash
#!/bin/bash
#SBATCH -a 1-1667

file=$(ls $FASTQ/*_1.fastq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)  ## the ls a basename argument will change depending on how your files are named
STEM=$(basename "$file" _1.fastq.gz)

STAR --genomeDir $WORKDIR/STARindex/ \
--runThreadN 8 \
--readFilesIn [List of fastq.gz files in batch. These must match the files specified in the manifest being used] \
--readFilesCommand zcat \
--outFileNamePrefix  $WORKDIR/aligned/$STEM \
--outSAMtype BAM Unsorted \
--outFilterMismatchNoverReadLmax 0.06 \
--outSAMattrRGline ID: \
--soloType SmartSeq \
--readFilesManifest $WORKDIR/rawfiles/[manifest file name] \
--soloUMIdedup Exact \
--soloStrand Unstranded
```
###**Using Seurat to perform QC and eliminate low quality cells**
Move the output files from STARsolo (barcodes.tsv, features.tsv,matrix.mtx) to a new directory $WORKDIR/Data so that there are subdirectories for each batch that was run for the genome alignment ($WORKDIR/Data/WT1/, $WORKDIR/Data/WT2/, $WORKDIR/Data/WT3/, etc.). Then open an R session and run the "1_SeuratQCandcellfiltering.R" script. This will take all the files, create a Seurat Object, perform basic QC and generate graphs. THen it filter cells based on several criteria: at least 1,000,000 reads, at least 9,000 genes, a log10genes/reads ratio of 0.58, less than 10% mitochondrial RNAs. Then it will generate the QC graphs on the filtered dataset.


