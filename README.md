# scDARTHEKcells

**Overview**
This contains information necessary for reproducing the identification of m<sup>6</sup>A 

**Prerequisites**
1) Download raw fastq files
2) Download flexbar (v3.0.3) and STAR (v2.7.5a), preferably into conda environment
3) Download Seurat (v3.2 or later)
4) Download genome and annotation files
5) Create STAR index for genome alignment


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
Then rename all the files to include a cell number and transgene expression status
```bash

```
Then put them all in the same rawfiles directory
```bash

```

***2) Download flexbar and STAR***
If using conda environment.
```bash
conda install -c bioconda flexbar=3.0.3
conda install -c bioconda STAR=2.7.5a
```

***3) Download Seurat***
In R 3.6 or later.
```R
install.packages("Seurat")
```
For more installation instructions for Seurat, see https://satijalab.org/seurat/articles/install.html
For more detailed tutorials for Seurat, see https://satijalab.org/seurat/articles/get_started.html

***4) Download genome and annotation files***
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

***5) Create STAR index***


