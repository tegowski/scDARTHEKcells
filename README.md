# scDART in HEK293T cells

## **Overview**
This contains information necessary for reproducing the identification of m<sup>6</sup>A. In this example, 300 APOBEC1-YTH and 300 APOBEC1-YTHmut expressing cells will be analyzed. Please replace all instances of $WORKDIR with path to your working directory. All code files should be downloaded to a directory $WORKDIR/software. $WORKDIR is the full path to the working directory being used for this analysis. Any instance of $WORKDIR present in code should either be defined as or replaced with the full path to the working directory.

**Prerequisites**
1) Download raw fastq files
2) Download and install flexbar (v3.0.3) and STAR (v2.7.5a), preferably into conda environment
3) Download Seurat (v3.2 or later)
4) Download genome and annotation files
5) Create STAR index for genome alignment


***1) Download raw fastq files***
Raw sequencing files are stored at SRA under the accession number GSE180954.
Files with SRA run numbers between SRR15268425-SRR15268461 and SRR15268889-SRR15270105 are APOBEC1-YTH-expressing, while those between SRR15268461-SRR15268876 are APOBEC1-YTHmut-expressing.
Downloading could take several days, and the files are nearly 2TB of memory in total.
```bash
#In the $WORKDIR/software directory, run the script to download SMART-seq single-cell fastq files
sbatch fastqdump.sh
```
Rename all the APOBEC1-YTH-expressing files to include a cell number and transgene expression status.
Move the renameYTH.sh script to $WORKDIR/YTH and submit the job
```bash
sbatch renameYTH.sh
```
Rename all the APOBEC1-YTHmut-expressing files to include a cell number and transgene expression status.
Move the renameYTHmut.sh script to $WORKDIR/YTHmut and submit the job
```bash
sbatch renameYTHmut.sh
```
Once you are **sure** the files have been renamed and moved to the $WORKDIR/rawfiles directory, you can delete the originals. The scripts above copy the files, preventing the need to redownload in the case of an error.
```bash
rm -r $WORKDIR/YTH
rm -r $WORKDIR/YTHmut
```

***2) Install Bullseye, flexbar, and STAR***
It is easiest to install all requisite materials in a conda environment if you don't have root access. Please see https://github.com/mflamand/Bullseye for detailed instructions for installing Bullseye and requisite Perl modules.


First, download Bullseye scripts from https://github.com/mflamand/Bullseye and move them to $WORKDIR/software

```bash
#Navigate to $WORKDIR/software
#Create conda environment and install required Perl modules for Bullseye

conda env create -f bullseye.yml
conda activate Bullseye
wget http://www.cpan.org/authors/id/T/TO/TODDR/XML-Parser-2.46.tar.gz 
tar -xf XML-Parser-2.46.tar.gz
cd XML-Parser-2.46 
perl Makefile.PL EXPATLIBPATH=$CONDA_PREFIX/lib EXPATINCPATH=$CONDA_PREFIX/include
make
make install

cpanm Bio::DB::Fasta
cpanm Text::NSP
cpanm Array::IntSpan
cpanm MCE

#Install flexbar and STAR
conda install -c bioconda flexbar=3.0.3
conda install -c bioconda STAR=2.7.5c
conda install -c r r=4.0.0
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
for (( i = 1; i <= 22; i++ ))
  do
wget http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.$i.fa.gz
  done
wget http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.X.fa.gz
wget http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.Y.fa.gz
wget http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.MT.fa.gz
wget http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.nonchromosomal.fa.gz
gunzip *
```
To download gtf annotation file:
```bash
#Navigate to $WORKDIR/gtf
mkdir $WORKDIR/gtf
cd $WORKDIR/gtf
wget 'http://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz'
gunzip *
```

***5) Create STAR index***
```bash
cd $WORKDIR/software
sbatch STARindex.sh
```

## **Data Processing**
###**Adapter trimming and genome alignment on Linux/Unix system**
The first two steps after obtaining fastq files occur on a Linux/Unix system
1) Trimming of adapter sequences
2) Genome alignment and feature counting

***1) Trimming of adapter sequences***
Flexbar (3.0.3) is used for trimming adapter sequences on reads. Since Nextera-compatible indexing primers were used, we can use the -aa Nextera option. If analyzing other adapter sequences, you may need to specify your adapter sequences in a fasta file. It is suggested to submit as array jobs as there are so many files. If using a system with slurm manager, the #SBATCH -a option can be used. This should be changed if using other systems.

```bash
cd $WORKDIR/software
sbatch flexbar.sh
```

***2) Genome alignment and feature counting***
This step can take alot of time if all cells are processed in one job. There therefore it is recommended to break up the files into several "batches" that can be run simultaneously. Therefore 5 different groups of "WT" or APOBEC1-YTH-expressing cells and 3 different groups of "mut" or APOBEC1-YTHmut-expressing cells can be run at a time as separate jobs. These can be re-integrated later in Seurat.

Copy the scripts to the rawfiles directory and then submit them.

```bash
cp $WORKDIR/software/StarsoloalignWT.sh $WORKDIR/rawfiles/
cp $WORKDIR/software/StarsoloalignMUT.sh $WORKDIR/rawfiles/
cd $WORKDIR/rawfiles
sbatch StarsoloalignWT.sh
sbatch StarsoloalignMUT.sh
```
###**Using Seurat to perform QC and eliminate low quality cells**
The output files will be in a directory $WORKDIR/aligned/$STEMSolo.out/Gene/filtered. Move the output files from STARsolo (barcodes.tsv, features.tsv,matrix.mtx) to a new directory $WORKDIR/Data so that there are subdirectories for each batch that was run for the genome alignment ($WORKDIR/Data/WT1/, $WORKDIR/Data/WT2/, $WORKDIR/Data/WT3/, etc.). Then open an R session and run the "1_SeuratQCandcellfiltering.R" script. This will take all the files, create a Seurat Object, perform basic QC and generate graphs. THen it filter cells based on several criteria: at least 1,000,000 reads, at least 9,000 genes, a log10genes/reads ratio of 0.58, less than 10% mitochondrial RNAs. Then it will generate the QC graphs on the filtered dataset.

###Using Bullseye to identify m<sup>6</sup> methylation
*Please see https://github.com/mflamand/Bullseye for instructions on installation of Bullseye*

There are 4 Major Steps
1) Genome alignment
2) Bamfile parsing
3) Identification of C-to-U mutations
4) Filtering of results for high confidence sites

***1) Genome alignment***
For the analysis of methylation, it is easiest to align the sequences to the genome again using STAR, but not the STARsolo option. This will generate a separate bamfile for each cell. It is easiest to do this as an array, so that many files be aligned simultaneously. The #SBATCH -a option allows for this with the SLURM managment system. If you are not using SLURM, you may need to change the command. If these are run sequencially the job will be a very long time.
```bash
cp $WORKDIR/software/STARBullseyealignWT.sh $WORKDIR/rawfiles/
cp $WORKDIR/software/STARBullseyealignMUT.sh $WORKDIR/rawfiles/
cd $WORKDIR/rawfiles
sbatch STARBullseyealignWT.sh
sbatch STARBullseyealignMUT.sh
```


***2) Parse bamfiles***
First, parse the bamfiles of each APOBEC1-YTH-expressing cell. This step uses samtools to eliminate optical duplicates from the bamfile before using Bullseye to parse the reads and determine the nucleotide composition at every location with a minimum coverage.
```bash
cp $WORKDIR/software/MakematrixYTH.sh $WORKDIR/singlebams
cd $WORKDIR/singlebams
sbatch MakematrixYTH.sh
```

Optical duplicates are removed from each APOBEC1-YTHmut-expressing cell.
```bash
cp $WORKDIR/software/RemoveduplicatesMUT.sh $WORKDIR/singlebams
cd $WORKDIR/singlebams
sbatch RemoveduplicatesMUT.sh
```

Then the APOBEC1-YTHmut bamfiles are merged and that merged bamfile is parsed. This creates an average nucleotide representation matrix across all APOBEC1-YTHmut cells sequenced.
```bash
cp $WORKDIR/software/MakematrixYTHmut.sh $WORKDIR/singlebams
cd $WORKDIR/singlebams
sbatch Make_matrixYTHmut.sh
```

***3) Find C-to-U mutations***
The next step is to compare the C-to-U (C-to-T in DNA sequencing) mutations found in each APOBEC1-YTH cell to the mutation rate in the average APOBEC1-YTHmut cells. By using thresholds for minimum coverage, minimum/maximum C2U editing percentage, and a fold-change in editing over APOBEC1-YTHmut cells, we can filter out noise due to sequencing errors/SNPs and many APOBEC1 (non-YTH)-mediated editing events.

This step requires a refFlat annotation file. Please use the annotation file provided here if repeating this dataset. If using a different (non-human) dataset, please see https://github.com/mflamand/Bullseye for instructions on how to generate a refFlat gtf file for your organism of interest.

```bash
cp $WORKDIR/software/find_RNA_edit_sites.sh $WORKDIR/matrix
cd $WORKDIR/matrix
sbatch find_RNA_edit_sites.sh
```

***4) Add cell ID to bed files and make one file with all sites in all cells***
This step will take the Cell ID (e.g. WTCell1) contained within the bed file names and add it as a column to the each of the bed files.

```bash
cp $WORKDIR/software/addCellIDs.sh $WORKDIR/bed
cd $WORKDIR/bed
sbatch addCellIDs.sh
```

Then concatenate all the bedfiles together

```bash
cat 2_WTCell*.bed > Sitesinallcells.bed
```

***5) Filter for sites that occur in R-A-C motif (G/A-A-C).*** 
This enriches for sites with C2U editing adjacent to m6A site.

```bash
cp $WORKDIR/software/RACfilter.sh $WORKDIR/bed
cd $WORKDIR/bed
sbatch RACfilter.sh
```

***6) Eliminate sites that were found to have editing in HEK293T cells when APOBEC1 was overexpressed by itself.***
This should only eliminate a very small proportion of sites and only offers a modest improvement in site quality. If using a system other than HEK293T cells, this step can be omitted.

```bash
cd $WORKDIR/bed
R
```
Note: replace $WORKDIR in the below script with your working directory path

```R
> library(tidyverse)
> apoonly <- read.table("DARTseq_APOonly_hg382.bed", col.names=c("chr","start","end"))
> head(apoonly)
> dartsites <- read.table("Sitesinallcells.3nt.RAC.bed")
> head(dartsites)
> #Add position column to compare with
> apoonly <- apoonly %>% mutate(pos = paste0(chr,"-",start))
> dartsites <- dartsites %>% mutate(pos = paste0(V1,"-",V2))
> #Filter dartsites for only those that don't appear in apoonly
> dartsitesfilt <- dartsites %>% filter(!pos %in% apoonly$pos)
> write.table(dartsitesfilt,file="Sitesinallcells_APOfilt.bed",col.names=F,row.names=F,quote=F,sep="\t")
```



