#!/bin/bash
#SBATCH --mem=4G

for i in ls WTCell*.bed;
do
echo "processing $i"
awk '{print $1"\t"$2"\t"$3"\t"FILENAME}' $i > $i.modbed
cut -f 4 $i.modbed > $i.Ids.txt
awk -F '[-_]' 'OFS="\t" {print $2$3$4}' $i.Ids.txt > $i.Ids2.txt
paste $i $i.Ids2.txt > 2_$i.bed
rm $i.modbed $i.Ids.txt $i.Ids2.txt;
done
