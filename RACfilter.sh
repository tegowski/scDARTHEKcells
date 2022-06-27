#!/bin/bash
#SBATCH -p scavenger
#SBATCH --mem=50G

echo -e 'Input bed file'
head $1.bed
#Remove chr from first column
#awk '{print substr($1,4), $2, $3, $4, $5, $6}' $1.bed6 > $1.nochr.bed
#echo -e '\nchr removed:'
#head -5 $1.nochr.bed
#add 2nts upstream of C to T site if + strand, if - add downstream
awk '{ if ($6 == "+")
        print $1"\t"$2-2"\t"$3"\t"$4"\t"$5"\t"$6;
else
        print $1"\t"$2"\t"$3+2"\t"$4"\t"$5"\t"$6
}' $1.bed > $1.3nt.bed
echo -e '\n2nts added 5-prime of C:'
head -5 $1.3nt.bed
#get sequence at those sites
bedtools getfasta -bedOut -fo $1.3nt.sequence.bed -s -fi /datacommons/meyerlab/userdata/matt_tego/genomes/hg38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fasta -bed $1.3nt.bed

#Add sequence to tab
paste $1.bed $1.3nt.sequence.bed > $1.3ntwithseq.bed
echo -e '\n3mer context of edited C:'
head $1.3ntwithseq.bed
#retrieve lines that have RAC motif from pos file
awk 'BEGIN{IGNORECASE=1}{if($13~"GAC|AAC") {print $0}}' $1.3ntwithseq.bed > $1.3nt.RAC.bed
awk 'BEGIN{IGNORECASE=1}{if($13!~"GAC|AAC") {print $0}}' $1.3ntwithseq.bed > $1.3nt.nonRAC.bed
rm $1.3nt.sequence.bed $1.nochr.bed $1.3nt.bed $1.3ntwithseq.bed
