#!/bin/bash




file=$1
bam='echo $file |tr "/" "\n"|tail -1'
fai=$2
filter=$3 #1 or 0 for clr subread


samtools view -@ 2 $file | samToBed /dev/stdin/ --useH --flag   > samBed.$bam.bed

if [filter==1];then
	cat samBed.$bam.bed | python filter.subread.py |sort -k1,1 -k2,2n >samBed.$bam.bed.f
else
	mv samBed.$bam.bed samBed.$bam.bed.f
fi



intersectBed -c -a <(bedtools makewindows -b <( awk 'BEGIN{OFS="\t"}{print $1,0,$2}' $fai) -w 100 |intersectBed -v -a stdin -b lc.bed) -b samBed.$bam.bed.f > coverage.$bam.bed

awk 'BEGIN{OFS="\t";c=0;sum=0;} sum=sum+$4;c=c+1;END{print sum/c;}' coverage.$bam.bed |tail -1> mean.$bam.txt

mean=$(cat mean.$bam.txt)
echo $mean

#first run
for r in ` cat $fai|cut -f 1`;do echo $r; viterbi <(grep -w $r coverage.$bam.bed |cut -f 4 ) 1 pre.$bam.$r 1 $mean; done
for r in ` cat $fai|cut -f 1`;do echo $r; cat pre.$bam.$r.viterout.txt >> pre.$bam.viterout.txt;done

paste <(cat coverage.$bam.bed) <(cat pre.$bam.viterout.txt) > pre.$bam.viterout.bed

bash viter.to_call.sh pre.$bam.viterout.bed 

cat DUPcalls.masked_CN.pre.$bam.viterout.bed | awk '$4==3' -| awk 'BEGIN{OFS="\t";c=0;sum=0;} sum=sum+($3-$2);c=c+1;END{print sum/(c*100);}' - |tail -1 >scaler.$bam.txt
scaler=$(cat scaler.$bam.txt)
echo $scaler
 
#scaled run
for r in ` cat $fai|cut -f 1`;do echo $r; viterbi <(grep -w $r coverage.$bam.bed |cut -f 4 ) $scaler $bam.$r 1 $mean ; done
for r in ` cat $fai|cut -f 1`;do echo $r; cat $bam.$r.viterout.txt >> $bam.viterout.txt;done

paste <(cat coverage.$bam.bed) <(cat $bam.viterout.txt) > $bam.viterout.bed

bash viter.to_call.sh $bam.viterout.bed

rm *pre.*

#plot every 50000 points ~ 5MB per plot
Rscript plot.HMM.noclip.R $bam.viterout.bed $bam 50000 mean.$bam.txt
