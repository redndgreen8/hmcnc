import os
import tempfile
import subprocess
import os.path
import json
import argparse


# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)
RD = "/home/cmb-16/mjc/rdagnew/summerproj"

#config("hmm_caller.json")
#configfile: "sd_analysis.json"



#ap = argparse.ArgumentParser(description=".")
#ap.add_argument("--subread", help="require subread filtering step,0/1", required=True)

filt=config["subread"]
#args=ap.parse_args()

#filt=args.subread


outdir="hmm"


#hg38.fai for alignments  # needs to mask using lc.bed 

fai= open(config["asm"]+".fai")
contigs = [l.split()[0].strip().replace("|","_") for l in fai]
bamt = config["bam"]
bam = bamt.split("/")[-1]
prefix_bam = os.path.splitext(bam)[0]

rule all:
    input:
        bed="hmm/cov.bed",
        no_subreadbed="hmm/cov.no_subread.bed",
        bins="hmm/coverage.bins.bed",
        avg="hmm/mean_cov.txt",
        vitterout=expand("hmm/{ctg}.viterout.txt", ctg=contigs),
        cn=expand("hmm/copy_number.{ctg}.bed", ctg=contigs),
        allCN="hmm/copy_number.tsv",
        plot=expand("hmm/{bm}.noclip.pdf",bm=prefix_bam),
	scale=("hmm/scaler.txt"),

rule MakeCovBed:
    input:
        bam=config["bam"],
    output:
        bed="hmm/cov.bed",
    params:
        rd=RD,
    shell:"""
mkdir -p hmm
samtools view -@ 2 {input.bam} | {params.rd}/samToBed /dev/stdin/ --useH --flag   > {output.bed}
"""

rule FilterSubreads:
    input:
        bed="hmm/cov.bed",
    output:
        covbed="hmm/cov.no_subread.bed",
    params:
        sd=SD,
        f=filt,
    shell:"""

if [{params.f} == 1];then
    cat {input.bed} | python filter.subread.py |sort -k1,1 -k2,2n > {output.covbed}
else
    mv {input.bed} {output.covbed}
fi
"""

rule MakeIntersect:
    input:
        bed="hmm/cov.no_subread.bed",
    output:
        bins="hmm/coverage.bins.bed",
    params:
        rd=RD,
        asm=config["asm"],
    shell:"""

intersectBed -c -a <(bedtools makewindows -b <( awk 'BEGIN{{OFS="\t"}}{{print $1,0,$2}}' {params.asm}.fai) -w 100) -b {input.bed} > {output.bins}
"""

rule GetMeanCoverage:
    input:
        bins="hmm/coverage.bins.bed",
    output:
        avg="hmm/mean_cov.txt",
    params:
        rd=RD,
    shell:"""
awk 'BEGIN{{OFS="\t";c=0;sum=0;}} sum=sum+$4;c=c+1;END{{print sum/c;}}' {input.bins} |tail -1> {output.avg}
"""

rule PreRunVitter:
    input:
        avg="hmm/mean_cov.txt",
        bins="hmm/coverage.bins.bed",
    output:
        cov="hmm/pre.{contig}.viterout.txt",
    params:
        rd=RD,
        contig_prefix="{contig}",
    shell:"""
mean=$(cat {input.avg})
touch {output.cov}
./viterbi <(grep -w {wildcards.contig} {input.bins} |cut -f 4 ) 1 hmm/pre.{params.contig_prefix} 1 $mean

"""

rule PreorderVitter:
    input:
        CopyNumber="hmm/pre.{contig}.viterout.txt",
        bins="hmm/coverage.bins.bed",
    output:
        cn="hmm/pre.copy_number.{contig}.bed",
    shell:"""

paste <(cat {input.bins}|grep -w {wildcards.contig} ) <(cat {input.CopyNumber}  ) > {output.cn}

"""

rule PrecombineVitter:
    input:
        allCopyNumberBED=expand("hmm/pre.copy_number.{contig}.bed", contig=contigs),
    output:
        allCN="hmm/pre.copy_number.tsv",
	#dupCall="hmm/pre.dup_call.tsv",
	scale="hmm/scaler.txt",
    shell:"""

cat {input.allCopyNumberBED} > {output.allCN}
pref='echo {output.allCN} |tr "/" "\n"|tail-1'
bash viter.to_call.sh {output.allCN}
cat DUPcalls.masked_CN.$pref.viterout.bed | awk '$4==3' -| awk 'BEGIN{{OFS="\t";c=0;sum=0;}} sum=sum+($3-$2);c=c+1;END{{print sum/c;}}' - |tail -1 > {output.scale}



"""



rule RunVitter:
    input:
        avg="hmm/mean_cov.txt",
        bins="hmm/coverage.bins.bed",
	scale="hmm/scaler.txt",
    output:
        cov="hmm/{contig}.viterout.txt",
    params:
        rd=RD,
        contig_prefix="{contig}",
    shell:"""
mean=$(cat {input.avg})
scaler=$(cat {input.scale})
touch {output.cov}
./viterbi <(grep -w {wildcards.contig} {input.bins} |cut -f 4 ) $scaler hmm/pre.{params.contig_prefix} 1 $mean

"""

rule orderVitter:
    input:
        CopyNumber="hmm/{contig}.viterout.txt",
        bins="hmm/coverage.bins.bed",
    output:
        cn="hmm/copy_number.{contig}.bed",
    shell:"""

paste <(cat {input.bins}|grep -w {wildcards.contig} ) <(cat {input.CopyNumber}  ) > {output.cn}

"""

rule combineVitter:
    input:
        allCopyNumberBED=expand("hmm/copy_number.{contig}.bed", contig=contigs),
    output:
        allCN="hmm/copy_number.tsv",
    shell:"""

cat {input.allCopyNumberBED} > {output.allCN}
rm hmm/pre.*

"""


rule PlotBins:
    input:
        allCN="hmm/copy_number.tsv",
        #aln=config["aln"],
        avg="hmm/mean_cov.txt",
    output:
        plot="hmm/{prefix_bam}.noclip.pdf",
    params:
        rd=RD,
        genome_prefix="{prefix_bam}",
    shell:"""

#plot every 50000 points ~ 5 MegaBases
Rscript plot.HMM.noclip.R {input.allCN} {params.genome_prefix} 50000 {input.avg}
touch {output.plot}
"""
