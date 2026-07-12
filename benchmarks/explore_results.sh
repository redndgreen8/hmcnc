#!/usr/bin/env bash
# =============================================================================
# explore_results.sh — Cross-phase benchmarking and concordance exploration
#
# NOTE: Ground truth quality not yet validated. These commands are for
# exploratory analysis only. Results should be interpreted cautiously until
# the truth set is confirmed appropriate for CNV benchmarking on chr22.
#
# Prerequisites:
#   conda activate hmcnc   (bedtools, bcftools, tabix, bgzip, truvari)
#   Run from repo root: bash benchmarks/explore_results.sh
# =============================================================================

set -euo pipefail
cd "$(dirname "$0")/.."   # always run from repo root

TRUTH_VCF="benchmarks/truth/HG002.chr22.cnv.truth.vcf.gz"
PHASES="phase1.1 phase2 phase3 phase4 phase5 phase5.6"


# -----------------------------------------------------------------------------
# 1. IGV TRACKS — build per-phase visualization sessions
#    Opens in IGV: File → Open Session → results/<phase>/igv/HG002/HG002.session.xml
# -----------------------------------------------------------------------------
build_igv_tracks() {
    local phase="${1:-phase5}"
    echo "=== Building IGV tracks for $phase ==="
    snakemake -s benchmarks/benchmark.smk --configfile benchmarks/config.yaml \
        --config phase="$phase" -j3 tracks
    echo "Session: results/$phase/igv/HG002/HG002.session.xml"
}


# -----------------------------------------------------------------------------
# 2. TRUVARI — precision/recall against HG002 chr22 CNV truth
#    Truth: benchmarks/truth/HG002.chr22.cnv.truth.vcf.gz
#    Matching: position ±500 bp, no sequence similarity required (-p 0.0)
#    NOTE: truth quality unverified — treat results as preliminary
# -----------------------------------------------------------------------------
run_truvari_all_phases() {
    echo "=== Running truvari for all phases ==="
    for phase in $PHASES; do
        echo "--- $phase ---"
        snakemake -s benchmarks/benchmark.smk --configfile benchmarks/config.yaml \
            --config phase="$phase" -j3 bench
    done
}

# Collect all summary.tsv files into one table
collect_summaries() {
    echo "=== Collecting truvari summaries ==="
    echo -e "phase\tsample\tTP\tFP\tFN\tprecision\trecall\tf1"
    for phase in $PHASES; do
        f="results/$phase/summary.tsv"
        [ -f "$f" ] && tail -n +2 "$f"
    done
}


# -----------------------------------------------------------------------------
# 3. RECIPROCAL OVERLAP — bedtools-based, threshold sweep
#    More transparent than truvari; useful for exploring sensitivity to RO cutoff
#    NOTE: truth VCF used directly as -b; works if bcftools-indexed
# -----------------------------------------------------------------------------
reciprocal_overlap_sweep() {
    local sample="${1:-HG002}"
    local phase="${2:-phase5}"
    local calls_del="results/$phase/del.filtered_final.${sample}.bed"
    local calls_dup="results/$phase/dup.filtered_final.${sample}.bed"

    echo "=== Reciprocal overlap sweep: $sample $phase ==="
    echo "Truth: $TRUTH_VCF"
    echo ""

    for thresh in 0.1 0.25 0.5 0.75 0.9; do
        del_hits=$(bedtools intersect -a "$calls_del" -b "$TRUTH_VCF" \
                       -f "$thresh" -r -u 2>/dev/null | wc -l | tr -d ' ')
        dup_hits=$(bedtools intersect -a "$calls_dup" -b "$TRUTH_VCF" \
                       -f "$thresh" -r -u 2>/dev/null | wc -l | tr -d ' ')
        echo "RO >= $thresh:  DEL hits = $del_hits  DUP hits = $dup_hits"
    done
}


# -----------------------------------------------------------------------------
# 4. CROSS-PHASE CONCORDANCE — which calls are shared between phases
#    Useful for understanding what each phase gains or loses relative to a baseline
# -----------------------------------------------------------------------------
cross_phase_concordance() {
    local sample="${1:-HG002}"
    local base_phase="${2:-phase1.1}"
    local compare_phase="${3:-phase5}"
    local thresh="${4:-0.5}"

    local base="results/$base_phase/filtered_final.${sample}.bed"
    local comp="results/$compare_phase/filtered_final.${sample}.bed"

    echo "=== Cross-phase concordance: $sample ==="
    echo "Base:    $base_phase  ($(wc -l < "$base") calls)"
    echo "Compare: $compare_phase  ($(wc -l < "$comp") calls)"
    echo "RO threshold: $thresh"
    echo ""

    shared=$(bedtools intersect -a "$comp" -b "$base" -f "$thresh" -r -u | wc -l | tr -d ' ')
    only_comp=$(bedtools intersect -a "$comp" -b "$base" -f "$thresh" -r -v | wc -l | tr -d ' ')
    only_base=$(bedtools intersect -a "$base" -b "$comp" -f "$thresh" -r -v | wc -l | tr -d ' ')

    echo "Shared (in both):           $shared"
    echo "Only in $compare_phase:     $only_comp"
    echo "Only in $base_phase:        $only_base"
}

# Compare all phases to phase1.1 baseline
cross_phase_concordance_all() {
    local sample="${1:-HG002}"
    echo "=== Cross-phase concordance vs phase1.1 baseline ($sample) ==="
    echo -e "phase\tshared\tonly_in_phase\tonly_in_baseline"
    for phase in $PHASES; do
        base="results/phase1.1/filtered_final.${sample}.bed"
        comp="results/$phase/filtered_final.${sample}.bed"
        [ ! -f "$comp" ] && continue
        shared=$(bedtools intersect -a "$comp" -b "$base" -f 0.5 -r -u | wc -l | tr -d ' ')
        only_comp=$(bedtools intersect -a "$comp" -b "$base" -f 0.5 -r -v | wc -l | tr -d ' ')
        only_base=$(bedtools intersect -a "$base" -b "$comp" -f 0.5 -r -v | wc -l | tr -d ' ')
        echo -e "$phase\t$shared\t$only_comp\t$only_base"
    done
}


# -----------------------------------------------------------------------------
# 5. CALL COUNT SUMMARY — raw counts across all phases (no truth needed)
# -----------------------------------------------------------------------------
call_count_summary() {
    echo "=== Filtered call counts (chr22 only) ==="
    echo -e "phase\tsample\tDEL\tDUP\ttotal"
    for d in results/phase*/; do
        phase=$(basename "$d")
        for s in HG002 HG003 HG004; do
            delf="$d/del.filtered_final.${s}.bed"
            dupf="$d/dup.filtered_final.${s}.bed"
            [ ! -f "$delf" ] || [ ! -f "$dupf" ] && continue
            del=$(wc -l < "$delf" | tr -d ' ')
            dup=$(wc -l < "$dupf" | tr -d ' ')
            tot=$((del + dup))
            echo -e "$phase\t$s\t$del\t$dup\t$tot"
        done
    done
}


# =============================================================================
# MAIN — uncomment the functions you want to run
# =============================================================================

# call_count_summary

# build_igv_tracks phase5

# run_truvari_all_phases
# collect_summaries

# reciprocal_overlap_sweep HG002 phase5
# reciprocal_overlap_sweep HG002 phase1.1

# cross_phase_concordance HG002 phase1.1 phase5
# cross_phase_concordance_all HG002

echo "explore_results.sh: no functions called — uncomment the ones you want in MAIN."
