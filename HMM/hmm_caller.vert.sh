#!/bin/bash




bam=$1
fai=$2
filter=$3


samtools view -@ 2 $bam | samToBed /dev/stdin/ --useH --flag   > samBed.$bam.bed

if [filter==1];then
	cat samBed.$bam.bed | python filter.subread.py |sort -k1,1 -k2,2n >samBed.$bam.bed.f
else
	mv samBed.$bam.bed samBed.$bam.bed.f
fi



intersectBed -c -a <(bedtools makewindows -b <( awk 'BEGIN{OFS="\t"}{print $1,0,$2}' $fai) -w 100) -b samBed.$bam.bed.f > coverage.$bam.bed

awk 'BEGIN{OFS="\t";c=0;sum=0;} sum=sum+$4;c=c+1;END{print sum/c;}' coverage.$bam.bed |tail -1> mean.$bam.txt

mean=$(cat mean.$bam.txt)
echo $mean

#viterbi2.cpp
for r in ` cat $fai|cut -f 1`;do echo $r; viterbi <(grep -w $r coverage.$bam.bed |cut -f 4 ) $mean $bam.$r 1; done

for r in ` cat $fai|cut -f 1`;do echo $r; cat $bam.$r.viterout.txt >> $bam.viterout.txt;done

paste <(cat coverage.$bam.bed) <(cat $bam.viterout.txt) > $bam.viterout.bed

#plot every 50000 points ~ 5MB per plot
Rscript plot.HMM.noclip.R $bam.viterout.bed $bam 50000 mean.$bam.txt
