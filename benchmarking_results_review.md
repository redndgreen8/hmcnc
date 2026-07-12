# Benchmarking Results Review

> [!IMPORTANT]
> **All results are chr22 only.** The BAMs in `data/chr22_test/` are chr22 extracts from HG002/HG003/HG004.
> The reference is `hg38.chr22.fa`. No genome-wide runs have been done. Call counts are therefore
> not comparable to whole-genome CNV callers and should be interpreted as a development/tuning
> benchmark, not a production result.

## What Exists

**Truth VCFs** are present in `benchmarks/truth/`:
- `HG002.chr22.cnv.truth.vcf.gz` — CNV-specific subset (small, ~14KB)
- `HG002.chr22.truth.vcf.gz` — Full chr22 SV truth (~574KB)
- `HG002_GRCh38_CMRG_SV_v1.00.vcf.gz` — CMRG medical gene SV set
- `HG002_SVs_Tier1_v0.6.vcf.gz` — GIAB SV Tier 1 (GRCh37)
- `HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` — SNV/indel benchmark (wrong type)

> [!WARNING]
> **No truvari/precision-recall results exist yet.** All results below are raw and filtered call counts only — no sensitivity or PPV numbers have been computed against truth.

---

## Filtered Call Counts by Phase

Counts from `del.filtered_final.<SAMPLE>.bed` and `dup.filtered_final.<SAMPLE>.bed`
(post exclusion-region + repeat-mask filter). **Chr22 only** — chr22 has ~51 Mb of sequence
but only ~37 Mb is non-centromeric/accessible; real CNV density here is lower than genome-wide.

### Core Phases (Phase 0 → Phase 5)

| Phase | HG002 DEL | HG002 DUP | HG003 DEL | HG003 DUP | HG004 DEL | HG004 DUP |
|-------|-----------|-----------|-----------|-----------|-----------|-----------|
| phase0 (Poisson clip) | 349 | 573 | 369 | 461 | 387 | 531 |
| phase1 (ZINB MOM) | 349 | 573 | 386 | 596 | 387 | 531 |
| **phase1.1** (`--epsi-weight 1.0`) | **92** | **186** | **90** | **164** | **117** | **176** |
| phase2 (directional clips) | 92 | 186 | 90 | 164 | 117 | 176 |
| phase3 (NR dispersion) | 92 | 186 | 90 | 164 | 117 | 176 |
| phase4 (directional transitions) | 92 | 186 | 90 | 164 | 117 | 176 |
| phase5 (directional ZINB M-step) | 92 | 186 | 90 | 164 | 117 | 176 |

**Key observation**: The `--epsi-weight 1.0` penalty introduced in Phase 1.1 is the dominant effect,
cutting calls by ~70% (922 → 278 for HG002). Phases 2–5 (ZINB/directional clipping improvements)
leave call counts essentially unchanged, because coverage depth still dominates the HMM.

---

## Phase 5.x Sub-Experiments

A series of undocumented sub-experiments ran after Phase 5, exploring parameter space:

| Sub-phase | HG002 DEL | HG002 DUP | Notes |
|-----------|-----------|-----------|-------|
| phase5.01 | 92 | 186 | Baseline (same as phase5) |
| phase5.1 | 29 | 15 | Very different — likely higher epsi-weight or clip params |
| phase5.11 | 29 | 15 | Same as 5.1 |
| phase5.2 | 92 | 186 | Back to baseline |
| phase5.21 | 29 | 15 | Same low-call config as 5.1 |
| phase5.3 | 92 | 186 | Baseline |
| phase5.4 | 92 | 186 | Baseline (HG003/HG004 +1 DUP vs phase5) |
| **phase5.5** | **204** | **38** | Inverted ratio — likely clip-dominated run |
| phase5.6 | 92 | 186 | Near-baseline |

The **phase5.5** result is striking: DEL=204, DUP=38 — the DEL/DUP ratio is completely inverted from
all other phases. Clip params show `clipMean=13.8, clipPhi=0.1` vs typical `clipMean≈3, clipPhi≈1`,
suggesting the clip model dominated transitions and inverted the calling pattern.

The **phase5.1/5.11/5.21** family (DEL=29, DUP=15) represents a very aggressive filter or penalty
setting — these would need truth comparison to know if this is better or worse precision/recall.

---

## What the Numbers Mean (Without Truth)

Without truvari results, the only interpretable signal is:

1. **Phase 0 vs Phase 1.1 drop** (922 → 278 filtered): the epsi-weight penalty is very effective at
   suppressing short noisy calls. Whether it's too aggressive (lost true positives) or right-sized
   requires truth comparison.

2. **DEL/DUP ratio consistency**: Phases 2–5/5.01/5.2/5.3/5.4/5.6 all give DEL≈90–117,
   DUP≈164–186 — a consistent ~2:1 DUP:DEL ratio. Phase 5.5 is an outlier worth investigating.

3. **Trio consistency**: HG002, HG003, HG004 (Ashkenazim trio) show similar call counts across all
   phases, which is a good sanity check — the samples are from the same family and expected to have
   overlapping CNV profiles.

---

## Next Step: Run Truvari

Truth VCF is present for HG002 chr22. The benchmark pipeline supports it:

```bash
# Full benchmark with truvari (truth VCF already in benchmarks/truth/)
snakemake -s benchmarks/benchmark.smk --configfile benchmarks/config.yaml \
  --config phase=phase5 -j3 all
```

This will produce precision/recall/F1 in a `summary.json` inside the results directory.

> [!NOTE]
> Only HG002 truth is present. HG003/HG004 truth VCFs would need to be downloaded separately
> to get trio-level benchmarking.
