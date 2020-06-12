#!/usr/bin/bash


file=$1
bed='echo $file|tr "/" "\n"|tail -1'
cat $file | awk '
BEGIN{OFS="\t";s=0;e=0;cns=0;o=3;}
{
    if(NR==1)
    {
        s=$1"\t"$2; e=$3; cns=$5;
        chr = $1;
    }
    else if ($5==cns && chr==$1)
    {
        e=$3;
        o=0;
    }
    else
    {
        print s,e,cns;
        s=$1"\t"$2;
        e=$3;
        cns=$5;
        chr=$1;
        o=1;
    }
}
END{
if(o==0){print s,e,cns;}
} ' |intersectBed -v -a stdin -b hg38.telomere.extended.bed > ini_call.$bed

awk '$4<2' ini_call.$bed > DELcalls.$bed

intersectBed -wa -wb -a <( awk '$4>2' ini_call.$bed |cut -f 1-3) -b repeatMask.merged.bed |sort -k1,1 -k2,2n | python repeatMask.py | groupBy -g 1,2,3,10 -c 9| awk 'BEGIN{OFS="\t"} $6=$5/$4' >DUPcalls.rep_int.$bed

intersectBed -v -a <( awk '$4>2' ini_call.$bed |cut -f 1-3) -b repeatMask.merged.bed|awk 'BEGIN{OFS="\t"}{print$1,$2,$3,$3-$2,0,0}'>>DUPcalls.rep_int.$bed

intersectBed -wa -wb -a <(awk '$6<0.8' DUPcalls.rep_int.$bed) -b <( awk '$4>2' ini_call.$bed ) |awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$10,$6;}' > DUPcalls.masked_CN.$bed
