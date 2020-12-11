import os
import tempfile
import subprocess
import os.path
import json

# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)


RD = config["scr"]

#config("hmm_caller.json")


configfile: config["configfile"]

outdir=config["outdir"]+"/hmm"

fai= open(config["bam"]+".fai")
contigs = [l.split()[0].strip().replace("|","_") for l in fai]
bamt = config["bam"]
bamm = bamt.split("/")[-1]
prefix_bam = os.path.splitext(bamm)[0]




cov=config["coverage"]
sub=config["subread"]
mq=config["mq"]


if config["repeatMask"] not None:
    rep=config["repeatMask"]



rule all:
    input:
        bed="{outdir}/cov.bed",
        no_subreadbed="{outdir}/cov.no_subread.bed",
        bins="{outdir}/coverage.bins.bed.gz",
        avg="{outdir}/mean_cov.txt",
        vitterout=expand("{outdir}/{ctg}.viterout.txt", ctg=contigs),
        cn=expand("{outdir}/copy_number.{ctg}.bed", ctg=contigs),
        allCN="{outdir}/copy_number.tsv",
        plot=expand("{outdir}/{bm}.noclip.pdf",bm=prefix_bam),

rule MakeCovBed:
    input:
        bam=config["bam"],
    output:
        bed=protected("{outdir}/cov.bed"),
    params:
        rd=RD,
        out={outdir},
    shell:"""
mkdir -p {params.out}
samtools view -q 10 -F 2304 -@ 3 {input.bam} | {params.rd}/samToBed /dev/stdin/ --useH --flag   > {output.bed}
"""

if sub==0:
    rule FilterSubreads0:
        input:
            bed="{outdir}/cov.bed",
        output:
            covbed="{outdir}/cov.no_subread.bed",
        params:
            rd=RD,
        shell:"""
    	cp {input.bed} {output.covbed}
	"""
else:
    rule FilterSubreads1:
        input:
            bed="{outdir}/cov.bed",
        output:
            covbed="{outdir}/cov.no_subread.bed",
        params:
            sd=SD,
        shell:"""
    {params.sd}/RemoveRedundantSubreads.py --input {input.bed} |sort -k1,1 -k2,2n > {output.covbed}

           cat {input.bed} | python {params.rd}/filter.subread.py --mapq {mapQ} |sort -k1,1 -k2,2n > {output.covbed}
    """

rule MakeIntersect:
    input:
        bed="{outdir}/cov.no_subread.bed",
    output:
        windows="{outdir}/ref_windows.bed"
        bins=protected("{outdir}/coverage.bins.bed.gz"),
    params:
        rd=RD,
        asm={bamt},
    shell:"""
bedtools makewindows -b <( awk 'BEGIN {{OFS="\\t"}} {{print $1,0,$2}} ' {params.asm}.fai) -w 100 | intersectBed -v -a stdin -b {params.rd}/lc.bed {params.rd}/hg38.telomere.extended.bed {params.rd}/annotation/gap.bed {params.rd}/annotation/centromere.bed |bedtools sort > {output.windows}

intersectBed -sorted -c -a {output.windows} -b {input.bed}| bgzip -c > {output.bins}
tabix -C {output.bins}
"""
### annotation?



if cov is not None:
    rule GetMeanCoverage:
        input:
            bins="{outdir}/coverage.bins.bed.gz",
        output:
            avg="{outdir}/mean_cov.txt",
        params:
            rd=RD,
            cove=cov,
        shell:"""
    echo {cove}> {output.avg}
    """
else:
    rule GetMeanCoverage:
        input:
            bins="{outdir}/coverage.bins.bed.gz",
        output:
            avg="{outdir}/mean_cov.txt",
        params:
            rd=RD,
        shell:"""
    zcat {input.bins} | awk 'BEGIN{{OFS="\\t";c=0;sum=0;}} sum=sum+$4;c=c+1;END{{print sum/c;}}' | tail -1> {output.avg}
    """

rule RunVitter:
    input:
        avg="{outdir}/mean_cov.txt",
        bins="{outdir}/coverage.bins.bed.gz",
    output:
        cov=protected("{outdir}/{contig}.viterout.txt"),
    params:
        rd=RD,
        contig_prefix="{contig}",
    shell:"""
mean=$(cat {input.avg})
touch {output.cov}
tabix {input.bins} {wildcards.contig} | cut -f 4 | \
  {params.rd}/{outdir}/viterbi  /dev/stdin $mean {outdir}/{params.contig_prefix}

"""
## viterbi params ?



rule orderVitter:
    input:
        CopyNumber="{outdir}/{contig}.viterout.txt",
        bins="{outdir}/coverage.bins.bed.gz",
    output:
        cn="{outdir}/copy_number.{contig}.bed",
    shell:"""

paste <(tabix  {input.bins} {wildcards.contig} ) <(cat {input.CopyNumber}  ) > {output.cn}

"""

rule combineVitter:
    input:
        allCopyNumberBED=expand("{outdir}/copy_number.{contig}.bed", contig=contigs),
    output:
        allCN="{outdir}/copy_number.tsv",
    run:
        import sys
        sys.stderr.write("Writing " + str(output.allCN) + "\n")
        cn=open(output.allCN,'w')
        for fn in input.allCopyNumberBED:
            sys.stderr.write("writing " + str(fn) + "\n")
            f=open(fn)
            l=f.readlines()
            cn.write("".join(l))
        cn.close()


rule call:
    input:
        allCN="hmm/copy_number.tsv",
        #aln=config["aln"],
        avg="hmm/mean_cov.txt",
    output:
        call="DUPcalls.masked_CN.copy_number.tsv",
    params:
        rd=RD,
        genome_prefix=prefix_bam,
    shell:"""
bash {params.rd}/viter.to_call.sh {input.allCN} {params.rd}
touch {output.call}
}
"""
rule repeatMask:

rule compositeCall:
    input:
        call="DUPcalls.masked_CN.copy_number.tsv",
    output:
        compCall="DUPcalls.FINAL.composite.bed",
    params:
        rd=RD,
        asm=config["asm"],
    shell:"""
mergeBed -i <(cat {input.call})|sort -k1,1 -k2,2n |intersectBed -wb -b stdin -a {input.call} |groupBy -g 6,7,8 -c 2,4,5 -o collapse,collapse,collapse|awk '$3-$2 >15000'|intersectBed -wa -v -a stdin  -b <(slopBed -b 50000 -g {params.asm}.fai -i {params.rd}/annotation/gap.bed)  {params.rd}/annotation/centromere.bed >{output.compCall}
"""

rule repeatMask:







rule PlotBins:
    input:
        allCN="{outdir}/copy_number.tsv",
        #aln=config["aln"],
        avg="{outdir}/mean_cov.txt",
    output:
        plot="{outdir}/{prefix_bam}.noclip.pdf",
    params:
        rd=RD,
        genome_prefix="{prefix_bam}",
    shell:"""
touch {output.plot}
#plot every 50000 points ~ 5MB
#Rscript {params.rd}/plot.HMM.noclip.R {input.allCN} {params.genome_prefix} 50000 {input.avg}
touch {output.plot}
"""
