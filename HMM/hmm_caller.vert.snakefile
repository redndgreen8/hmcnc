import os
import tempfile
import subprocess
import os.path
import json
import argparse


# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)
RD = "/home/cmb-16/mjc/rdagnew/summerproj/hmm/HMM"

#config("hmm_caller.json")
#configfile: "sd_analysis.json"



#ap = argparse.ArgumentParser(description=".")
#ap.add_argument("--subread", help="require subread filtering step,0/1", required=True)

filt=int(config["subread"])
#args=ap.parse_args()

#filt=args.subread
threads=config["t"]
mapQ=config["MQ"]
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
        previtterout=expand("hmm/{ctg}.viterout.txt.pre", ctg=contigs),
        vitterout=expand("hmm/{ctg}.viterout.txt", ctg=contigs),
        pcn=expand("hmm/copy_number.{ctg}.bed.pre", ctg=contigs),
        cn=expand("hmm/copy_number.{ctg}.bed", ctg=contigs),
        pallCN="hmm/pre.copy_number.tsv",
        allCN="hmm/copy_number.tsv",    
        plot=expand("hmm/{bm}.noclip.pdf",bm=prefix_bam),
	scale="hmm/scaler.txt",

rule MakeCovBed:
    input:
        bam=config["bam"],
    output:
        bed="hmm/cov.bed",
    params:
        rd=RD,
    shell:"""
mkdir -p hmm
/home/cmb-16/mjc/rdagnew/anaconda3/envs/sam_env/bin/samtools view -@ 2 {input.bam} | {params.rd}/samToBed /dev/stdin/ --useH --flag   > {output.bed}
"""
if filt==0:
    rule FilterSubreads0:
        input:
            bed="hmm/cov.bed",
        output:
            covbed="hmm/cov.no_subread.bed",
        params:
            rd=RD,
        shell:"""
    	cp {input.bed} {output.covbed}
	"""
else:
    rule FilterSubreads1:
        input:
            bed="hmm/cov.bed",
        output:
            covbed="hmm/cov.no_subread.bed",
        params:
            rd=RD,
        shell:"""
       cat {input.bed} | python {params.rd}/filter.subread.py --mapq {mapQ} |sort -k1,1 -k2,2n > {output.covbed}
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

intersectBed -c -a <(bedtools makewindows -b <( awk 'BEGIN{{OFS="\t"}}{{print $1,0,$2}}' {params.asm}.fai) -w 100|intersectBed -v -a stdin -b {params.rd}/lc.bed) -b {input.bed} > {output.bins}
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
        pcov="hmm/{contig}.viterout.txt.pre",
    params:
        rd=RD,
        contig_prefix="{contig}",
    shell:"""
mean=$(cat {input.avg})
{params.rd}/viterbi <(grep -w {wildcards.contig} {input.bins} |cut -f 4 ) 1 hmm/pre.{params.contig_prefix} 1 $mean
mv hmm/pre.{params.contig_prefix}.viterout.txt {output.pcov}
"""

rule PreorderVitter:
    input:
        pCopyNumber="hmm/{contig}.viterout.txt.pre",
        bins="hmm/coverage.bins.bed",
    output:
        pcn="hmm/copy_number.{contig}.bed.pre",
    shell:"""

paste <(cat {input.bins}|grep -w {wildcards.contig} ) <(cat {input.pCopyNumber}  ) > {output.pcn}

"""

rule PrecombineVitter:
    input:
        pallCopyNumberBED=expand("hmm/copy_number.{contig}.bed.pre", contig=contigs),
    output:
        pallCN="hmm/pre.copy_number.tsv",
    params:
        rd=RD,
    shell:"""
cat {input.pallCopyNumberBED} > {output.pallCN}
"""

rule scaler:
    input:
        pallCN="hmm/pre.copy_number.tsv",
    output:
        scale="hmm/scaler.txt",
    params:
        rd=RD,
    shell:"""
pref=$(echo {input.pallCN} |tr "/" "\\n"|tail -1)
echo $pref
bash {params.rd}/viter.to_call.sh {input.pallCN}
cat DUPcalls.masked_CN.$pref | awk '$4==3' -| awk 'BEGIN{{OFS="\t";c=0;sum=0;}} sum=sum+($3-$2);c=c+1;END{{print sum/(c*100);}}' - |tail -1 > {output.scale}

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
{params.rd}/viterbi <(grep -w {wildcards.contig} {input.bins} |cut -f 4 ) $scaler hmm/{params.contig_prefix} 1 $mean
touch {output.cov}

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
bash {params.rd}/viter.to_call.sh {input.allCN}
#plot every 50000 points ~ 5 MegaBases
#Rscript {params.rd}/plot.HMM.noclip.R {input.allCN} {params.genome_prefix} 50000 {input.avg}
touch {output.plot}

"""
