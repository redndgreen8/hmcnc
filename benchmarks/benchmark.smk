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

    # Prepare CNV truth VCF (DEL + INS→DUP, once per truth set):
    snakemake -s benchmarks/benchmark.smk --configfile benchmarks/config.yaml -j1 prepare_truth

    # Full benchmark with truvari (needs benchmarks/truth/ populated):
    snakemake -s benchmarks/benchmark.smk --configfile benchmarks/config.yaml \
        --config phase=phase1 -j3 bench

Config keys:
    phase         - label for output directory (e.g. "phase0", "phase1")
    ref           - hg38 chr22 FASTA path (must have .fai)
    truth_dir     - directory with truth VCFs
    bam_dir       - directory with <SAMPLE>.chr22.bam (+.bai)
    binary        - path to hmcnclip binary
    excl_bed      - exclusion regions BED (default: HMM/annotation/hg38.region_to_EXCLUDE.bed)
    repeat_bed    - repeat mask BED      (default: HMM/annotation/hg38.repeatMask.merged.bed)
    extra_args    - extra hmcnclip CLI args (default: "")
    min_sv_size   - minimum SV size for benchmarking in bp (default: 1000)
    bench_samples - samples with truth VCFs available (default: ["HG002"])

Truth VCF preparation notes:
    DEL truth: PASS records with SVLEN >= min_sv_size from HG002_SVs_Tier1_v0.6.vcf.gz
    DUP truth: INS records with REPTYPE=DUP and SVLEN >= min_sv_size, with END expanded
               to POS+SVLEN and SVTYPE renamed to DUP. These represent tandem duplications
               that appear as insertions in SV callers but as elevated coverage in CNV callers.
               (chr22: ~48 events >= 1kb, ~11 >= 5kb — sparse but directionally informative)
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

SAMPLES       = config.get("samples",       ["HG002", "HG003", "HG004"])
BENCH_SAMPLES = config.get("bench_samples", ["HG002"])   # samples with truth VCFs
MIN_SV_SIZE   = config.get("min_sv_size",   1000)        # bp; filters both truth and calls
OUTDIR        = f"results/{PHASE}"
EXTRA_ARGS    = config.get("extra_args", "")

# ── Targets ──────────────────────────────────────────────────────────────────

# Default: full pipeline including truvari (requires truth VCFs)
rule all:
    input:
        f"{OUTDIR}/summary.tsv",
        expand(f"{OUTDIR}/truvari/{{sample}}/summary.json", sample=BENCH_SAMPLES),

# truvari + summary (alias for all, explicit name for CLI use)
rule bench:
    input:
        f"{OUTDIR}/summary.tsv",
        expand(f"{OUTDIR}/truvari/{{sample}}/summary.json", sample=BENCH_SAMPLES),

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


# ── Step 4a: Prepare CNV truth VCF ───────────────────────────────────────────
# Builds HG002.chr22.cnv.truth.vcf.gz from Tier1 SV VCF:
#   DEL: PASS records with SVLEN >= MIN_SV_SIZE
#   DUP: INS records with REPTYPE=DUP and SVLEN >= MIN_SV_SIZE, END expanded to
#        POS+SVLEN so the record spans the duplicated region for overlap matching.
# Run once: snakemake -s benchmarks/benchmark.smk prepare_truth -j1

rule prepare_truth:
    input:
        vcf = f"{TRUTH_DIR}/HG002.chr22.truth.vcf.gz",
        tbi = f"{TRUTH_DIR}/HG002.chr22.truth.vcf.gz.tbi",
    output:
        vcf = f"{TRUTH_DIR}/HG002.chr22.cnv.truth.vcf.gz",
        tbi = f"{TRUTH_DIR}/HG002.chr22.cnv.truth.vcf.gz.tbi",
    params:
        min_size = MIN_SV_SIZE,
    run:
        import subprocess, re, os

        min_size = int(params.min_size)
        tmp = output.vcf.replace(".vcf.gz", ".tmp.vcf")

        header = subprocess.run(
            ["bcftools", "view", "-h", input.vcf],
            capture_output=True, text=True, check=True,
        ).stdout

        records = subprocess.run(
            ["bcftools", "view", "-H", input.vcf],
            capture_output=True, text=True, check=True,
        ).stdout

        kept_del = kept_dup = 0
        with open(tmp, "w") as out:
            out.write(header)
            for line in records.splitlines():
                f = line.split("\t")
                if len(f) < 8:
                    continue
                info = f[7]
                svtype_m = re.search(r"SVTYPE=(\w+)", info)
                svlen_m  = re.search(r"SVLEN=(-?\d+)", info)
                if not svtype_m or not svlen_m:
                    continue
                svtype = svtype_m.group(1)
                svlen  = abs(int(svlen_m.group(1)))
                if svlen < min_size:
                    continue

                if svtype == "DEL" and f[6] == "PASS":
                    out.write(line + "\n")
                    kept_del += 1
                elif svtype == "INS":
                    reptype_m = re.search(r"REPTYPE=(\w+)", info)
                    if not reptype_m or reptype_m.group(1) != "DUP":
                        continue
                    new_end = int(f[1]) + svlen
                    f[4]  = "<DUP>"
                    f[6]  = "PASS"
                    info  = re.sub(r"SVTYPE=INS", "SVTYPE=DUP", info)
                    info  = re.sub(r"END=\d+",    f"END={new_end}", info)
                    f[7]  = info
                    out.write("\t".join(f) + "\n")
                    kept_dup += 1

        subprocess.run(["bcftools", "sort", "-O", "z", "-o", output.vcf, tmp], check=True)
        subprocess.run(["tabix", "-p", "vcf", output.vcf], check=True)
        os.remove(tmp)
        print(f"prepare_truth: {kept_del} DEL + {kept_dup} DUP (from INS REPTYPE=DUP) >= {min_size}bp")


# ── Step 4b: Sort + bgzip + index VCF ────────────────────────────────────────

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


# ── Step 5: Truvari bench ─────────────────────────────────────────────────────
# Uses prepared CNV truth (DEL + INS→DUP). Only runs for BENCH_SAMPLES.
# DEL precision/recall: direct comparison.
# DUP precision/recall: approximate — truth DUPs are INS records with REPTYPE=DUP
#   whose coordinates have been expanded; not all DUP calls will have a truth match.

rule truvari_bench:
    input:
        vcf   = f"{OUTDIR}/{{sample}}.vcf.gz",
        tbi   = f"{OUTDIR}/{{sample}}.vcf.gz.tbi",
        truth = f"{TRUTH_DIR}/HG002.chr22.cnv.truth.vcf.gz",
        trtbi = f"{TRUTH_DIR}/HG002.chr22.cnv.truth.vcf.gz.tbi",
        ref   = REF,
    output:
        summary = f"{OUTDIR}/truvari/{{sample}}/summary.json",
        tp_base = f"{OUTDIR}/truvari/{{sample}}/tp-base.vcf.gz",
        tp_comp = f"{OUTDIR}/truvari/{{sample}}/tp-comp.vcf.gz",
        fp      = f"{OUTDIR}/truvari/{{sample}}/fp.vcf.gz",
        fn      = f"{OUTDIR}/truvari/{{sample}}/fn.vcf.gz",
    params:
        outdir   = f"{OUTDIR}/truvari/{{sample}}",
        min_size = MIN_SV_SIZE,
    shell:
        """
        rm -rf {params.outdir}
        truvari bench \
            -b {input.truth} \
            -c {input.vcf} \
            -o {params.outdir} \
            --pick multi \
            -r 500 \
            --passonly \
            --size-min {params.min_size} \
            -p 0.0
        """


# ── Step 6: Summary table ─────────────────────────────────────────────────────

rule summarize:
    input:
        jsons = expand(f"{OUTDIR}/truvari/{{sample}}/summary.json", sample=BENCH_SAMPLES),
    output:
        tsv = f"{OUTDIR}/summary.tsv",
    run:
        rows = []
        for sample, jpath in zip(BENCH_SAMPLES, input.jsons):
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
