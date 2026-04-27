"""
HMCNC Benchmarking + IGV Track Pipeline

Runs hmcnclip on each sample, filters calls via filter_finalize_call.sh, and
evaluates against GIAB truth VCFs using truvari. Run once per development phase.

Usage:
    conda activate hmcnc

    # Just run hmcnclip (no truth VCFs needed):
    snakemake -s benchmarks/benchmark.smk --configfile benchmarks/config.yaml \
        --config phase=phase1 -j3 calls

    # Run + filter with repeat/exclude masks:
    snakemake -s benchmarks/benchmark.smk --configfile benchmarks/config.yaml \
        --config phase=phase1 -j3 filter

    # Full benchmark with truvari (needs benchmarks/truth/ populated):
    snakemake -s benchmarks/benchmark.smk --configfile benchmarks/config.yaml \
        --config phase=phase1 -j3

Config keys:
    phase       - label for output directory (e.g. "phase0", "phase1")
    ref         - hg38 chr22 FASTA path (must have .fai)
    truth_dir   - directory with <SAMPLE>.chr22.truth.vcf.gz (+.tbi)
    bam_dir     - directory with <SAMPLE>.chr22.bam (+.bai)
    binary      - path to hmcnclip binary
    excl_bed    - exclusion regions BED (default: HMM/annotation/hg38.region_to_EXCLUDE.bed)
    repeat_bed  - repeat mask BED      (default: HMM/annotation/hg38.repeatMask.merged.bed)
    extra_args  - extra hmcnclip CLI args (default: "")
"""

import json
import os
import textwrap

# ── Config ────────────────────────────────────────────────────────────────────

PHASE      = config.get("phase",      "phase_test")
REF        = config.get("ref",        "/Users/red/repos/MethSmoothEval/data/annotations/hg38.chr22.fa")
TRUTH_DIR  = config.get("truth_dir",  "benchmarks/truth")
BAM_DIR    = config.get("bam_dir",    "data/chr22_test")
BINARY     = config.get("binary",     "src/hmcnclip")
EXCL_BED   = config.get("excl_bed",   "HMM/annotation/hg38.region_to_EXCLUDE.bed")
REPEAT_BED = config.get("repeat_bed", "HMM/annotation/hg38.repeatMask.merged.bed")
FILTER_SH  = config.get("filter_sh",  "HMM/filter_finalize_call.sh")

SAMPLES    = config.get("samples",    ["HG002", "HG003", "HG004"])
OUTDIR     = f"results/{PHASE}"
EXTRA_ARGS = config.get("extra_args", "")

# ── Targets ──────────────────────────────────────────────────────────────────

# Default: full pipeline including truvari (requires truth VCFs)
rule all:
    input:
        f"{OUTDIR}/summary.tsv",
        expand(f"{OUTDIR}/truvari/{{sample}}/summary.json", sample=SAMPLES),

# Just hmcnclip, no post-processing
rule calls:
    input:
        expand(f"{OUTDIR}/{{sample}}.vcf", sample=SAMPLES),

# hmcnclip + filter_finalize_call.sh (no truth VCFs needed)
rule filter:
    input:
        expand(f"{OUTDIR}/filtered_final.{{sample}}.bed", sample=SAMPLES),

# filter + IGV tracks (colored BED, bedGraphs, session XML)
rule tracks:
    input:
        expand(f"{OUTDIR}/igv/{{sample}}/{{sample}}.session.xml", sample=SAMPLES),


# ── Step 1: Run hmcnclip ─────────────────────────────────────────────────────

rule run_hmcnclip:
    input:
        bam = f"{BAM_DIR}/{{sample}}.chr22.bam",
        bai = f"{BAM_DIR}/{{sample}}.chr22.bam.bai",
        ref = REF,
    output:
        vcf   = f"{OUTDIR}/{{sample}}.vcf",
        bed   = f"{OUTDIR}/{{sample}}.bed",
        param = f"{OUTDIR}/{{sample}}.param.1",
        cov   = f"{OUTDIR}/{{sample}}.cov.bed",
        clip  = f"{OUTDIR}/{{sample}}.clip.bed",
        log   = f"{OUTDIR}/{{sample}}.log",
    params:
        prefix = f"{OUTDIR}/{{sample}}",
        extra  = EXTRA_ARGS,
    shell:
        """
        {BINARY} -a {input.bam} --output-all {params.prefix} {params.extra} {input.ref} \
            2>&1 | tee {output.log}
        """


# ── Step 2: Filter + finalize calls ──────────────────────────────────────────
# Runs HMM/filter_finalize_call.sh which:
#   - Keeps PASS calls with CN != 2
#   - Removes excluded regions (centromeres, etc.)
#   - Annotates with repeat overlap fraction
#   - Filters duplications where merged repeat fraction >= 0.8
# Output: filtered_final.<sample>.bed

rule filter_calls:
    input:
        bed    = f"{OUTDIR}/{{sample}}.bed",
        excl   = EXCL_BED,
        repeat = REPEAT_BED,
    output:
        final      = f"{OUTDIR}/filtered_final.{{sample}}.bed",
        filter_bed = f"{OUTDIR}/{{sample}}.bed.filter.bed",
        merge_bed  = f"{OUTDIR}/{{sample}}.bed.merge.filter.bed",
    params:
        outdir  = OUTDIR,
        base    = "{sample}.bed",
        script  = os.path.abspath(FILTER_SH),
        excl    = os.path.abspath(EXCL_BED),
        repeat  = os.path.abspath(REPEAT_BED),
    shell:
        """
        cd {params.outdir}
        bash {params.script} {params.base} {params.excl} {params.repeat}
        # script writes filtered_final.<base> in cwd; rename to match output
        mv filtered_final.{params.base} filtered_final.$(basename {wildcards.sample}).bed
        """


# ── Step 3: IGV Tracks ───────────────────────────────────────────────────────
# Produces per-sample igv/ directory with:
#   <sample>.cnv.bed.gz       colored BED9 (CNV calls, itemRgb by CN state)
#   <sample>.cov.bedgraph.gz  coverage per 100 bp bin
#   <sample>.clip.bedgraph.gz raw clip count per 100 bp bin
#   <sample>.session.xml      IGV session loading all three tracks

# CN → RGB color mapping (deletions = reds, duplications = blues)
CN_COLORS = {
    0: "204,0,0",     # CN0  homozygous del  dark red
    1: "255,140,0",   # CN1  het del         orange
    3: "100,149,237", # CN3  het dup         cornflower blue
    4: "0,0,205",     # CN4                  medium blue
    5: "0,0,139",     # CN5                  dark blue
    6: "75,0,130",    # CN6                  indigo
}
DEFAULT_COLOR = "150,150,150"

rule make_cnv_bed:
    """Convert filtered calls to colored BED9, sort, bgzip, tabix."""
    input:
        bed = f"{OUTDIR}/filtered_final.{{sample}}.bed",
    output:
        gz  = f"{OUTDIR}/igv/{{sample}}/{{sample}}.cnv.bed.gz",
        tbi = f"{OUTDIR}/igv/{{sample}}/{{sample}}.cnv.bed.gz.tbi",
    run:
        import subprocess, gzip, os
        cn_colors = {int(k): v for k, v in CN_COLORS.items()}
        os.makedirs(os.path.dirname(output.gz), exist_ok=True)
        tmp = output.gz.replace(".gz", ".tmp.bed")
        with open(input.bed) as fin, open(tmp, "w") as fout:
            for line in fin:
                f = line.rstrip().split("\t")
                chrom, start, end, cn = f[0], f[1], f[2], int(f[3])
                rgb   = cn_colors.get(cn, DEFAULT_COLOR)
                name  = f"CN{cn}"
                score = "0"
                fout.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t.\t{start}\t{end}\t{rgb}\n")
        subprocess.run(f"sort -k1,1 -k2,2n {tmp} | bgzip -c > {output.gz}", shell=True, check=True)
        subprocess.run(f"tabix -p bed {output.gz}", shell=True, check=True)
        os.remove(tmp)


rule make_cov_bedgraph:
    """Extract coverage column from cov.bed → bedGraph, sort, bgzip, tabix."""
    input:
        cov = f"{OUTDIR}/{{sample}}.cov.bed",
    output:
        gz  = f"{OUTDIR}/igv/{{sample}}/{{sample}}.cov.bedgraph.gz",
        tbi = f"{OUTDIR}/igv/{{sample}}/{{sample}}.cov.bedgraph.gz.tbi",
    shell:
        """
        mkdir -p $(dirname {output.gz})
        sort -k1,1 -k2,2n {input.cov} | bgzip -c > {output.gz}
        tabix -p bed {output.gz}
        """


rule make_clip_bedgraph:
    """Extract raw clip count (col4) from clip.bed → bedGraph, sort, bgzip, tabix."""
    input:
        clip = f"{OUTDIR}/{{sample}}.clip.bed",
    output:
        gz  = f"{OUTDIR}/igv/{{sample}}/{{sample}}.clip.bedgraph.gz",
        tbi = f"{OUTDIR}/igv/{{sample}}/{{sample}}.clip.bedgraph.gz.tbi",
    shell:
        """
        mkdir -p $(dirname {output.gz})
        awk 'OFS="\\t"{{print $1,$2,$3,$4}}' {input.clip} \
        | sort -k1,1 -k2,2n | bgzip -c > {output.gz}
        tabix -p bed {output.gz}
        """


rule make_igv_session:
    """Write an IGV session XML that loads the three tracks for this sample."""
    input:
        cnv_gz  = f"{OUTDIR}/igv/{{sample}}/{{sample}}.cnv.bed.gz",
        cov_gz  = f"{OUTDIR}/igv/{{sample}}/{{sample}}.cov.bedgraph.gz",
        clip_gz = f"{OUTDIR}/igv/{{sample}}/{{sample}}.clip.bedgraph.gz",
    output:
        xml = f"{OUTDIR}/igv/{{sample}}/{{sample}}.session.xml",
    run:
        # Use absolute paths so the session works regardless of working directory
        cnv_path  = os.path.abspath(input.cnv_gz)
        cov_path  = os.path.abspath(input.cov_gz)
        clip_path = os.path.abspath(input.clip_gz)
        sample    = wildcards.sample

        xml = textwrap.dedent(f"""\
            <?xml version="1.0" encoding="UTF-8" standalone="no"?>
            <Session genome="hg38" locus="chr22" version="8">
                <Resources>
                    <Resource path="{cnv_path}"/>
                    <Resource path="{cov_path}"/>
                    <Resource path="{clip_path}"/>
                </Resources>
                <Panel height="200" name="CNV" width="1200">
                    <Track attributeKey="{sample}.cnv.bed.gz" clazz="org.broad.igv.track.FeatureTrack"
                           displayMode="EXPANDED" fontSize="10" id="{cnv_path}"
                           name="CNV calls ({sample})" visible="true"/>
                </Panel>
                <Panel height="150" name="Coverage" width="1200">
                    <Track attributeKey="{sample}.cov.bedgraph.gz" clazz="org.broad.igv.track.DataSourceTrack"
                           autoScale="true" color="0,0,180" fontSize="10" id="{cov_path}"
                           name="Coverage ({sample})" renderer="BAR_CHART" visible="true"/>
                </Panel>
                <Panel height="100" name="Clipping" width="1200">
                    <Track attributeKey="{sample}.clip.bedgraph.gz" clazz="org.broad.igv.track.DataSourceTrack"
                           autoScale="true" color="180,0,0" fontSize="10" id="{clip_path}"
                           name="Clip signal ({sample})" renderer="BAR_CHART" visible="true"/>
                </Panel>
            </Session>
        """)
        with open(output.xml, "w") as f:
            f.write(xml)


# ── Step 4: Sort + bgzip + index VCF ─────────────────────────────────────────

rule prepare_vcf:
    input:
        vcf = f"{OUTDIR}/{{sample}}.vcf",
    output:
        gz  = f"{OUTDIR}/{{sample}}.vcf.gz",
        tbi = f"{OUTDIR}/{{sample}}.vcf.gz.tbi",
    shell:
        """
        bcftools sort -O z -o {output.gz} {input.vcf}
        tabix -p vcf {output.gz}
        """


# ── Step 4: Truvari bench ─────────────────────────────────────────────────────

rule truvari_bench:
    input:
        vcf   = f"{OUTDIR}/{{sample}}.vcf.gz",
        tbi   = f"{OUTDIR}/{{sample}}.vcf.gz.tbi",
        truth = f"{TRUTH_DIR}/{{sample}}.chr22.truth.vcf.gz",
        trtbi = f"{TRUTH_DIR}/{{sample}}.chr22.truth.vcf.gz.tbi",
        ref   = REF,
    output:
        summary = f"{OUTDIR}/truvari/{{sample}}/summary.json",
        tp_base = f"{OUTDIR}/truvari/{{sample}}/tp-base.vcf.gz",
        tp_comp = f"{OUTDIR}/truvari/{{sample}}/tp-comp.vcf.gz",
        fp      = f"{OUTDIR}/truvari/{{sample}}/fp.vcf.gz",
        fn      = f"{OUTDIR}/truvari/{{sample}}/fn.vcf.gz",
    params:
        outdir = f"{OUTDIR}/truvari/{{sample}}",
    shell:
        """
        rm -rf {params.outdir}
        truvari bench \
            -b {input.truth} \
            -c {input.vcf} \
            -o {params.outdir} \
            --pick multi \
            -r 500 \
            --svtype DEL DUP \
            --passonly \
            -p 0.0
        """


# ── Step 5: Summary table ─────────────────────────────────────────────────────

rule summarize:
    input:
        jsons = expand(f"{OUTDIR}/truvari/{{sample}}/summary.json", sample=SAMPLES),
    output:
        tsv = f"{OUTDIR}/summary.tsv",
    run:
        rows = []
        for sample, jpath in zip(SAMPLES, input.jsons):
            with open(jpath) as fh:
                d = json.load(fh)
            rows.append({
                "phase":     PHASE,
                "sample":    sample,
                "TP":        d.get("TP-base", d.get("TP", "?")),
                "FP":        d.get("FP", "?"),
                "FN":        d.get("FN", "?"),
                "precision": round(d.get("precision", float("nan")), 4),
                "recall":    round(d.get("recall",    float("nan")), 4),
                "f1":        round(d.get("f1",        float("nan")), 4),
            })

        header = ["phase","sample","TP","FP","FN","precision","recall","f1"]
        with open(output.tsv, "w") as out:
            out.write("\t".join(header) + "\n")
            for r in rows:
                out.write("\t".join(str(r[h]) for h in header) + "\n")

        import sys
        print("\t".join(header), file=sys.stderr)
        for r in rows:
            print("\t".join(str(r[h]) for h in header), file=sys.stderr)
