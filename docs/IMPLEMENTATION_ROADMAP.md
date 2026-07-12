# HMCNC Implementation Roadmap

## Overview

HMCNC (Hidden Markov Copy Number Caller) detects copy number variants using an HMM with multiple evidence layers:
- **Read depth** (coverage per 100bp bin)
- **Clipping signatures** (soft-clipped reads indicating breakpoints)
- **SNV allele frequencies** (B-allele frequency support)

This document outlines the phased implementation plan to incorporate ZINB clipping models and directional clipping evidence as specified in the mathematical definitions.

---

## Pre-Phase: Testing Infrastructure [COMPLETED]

### Achievements
1. **Build system fixed** for macOS (Makefile `-Wl,-rpath`)
2. **Conda environment** `hmcnc` with g++ 15.2.0, htslib, samtools
3. **CLI extensions** for single-chromosome testing:
   - `--wg-mean`, `--wg-var` - Override coverage statistics
   - `--wg-clip-mean`, `--wg-clip-var` - Override clipping statistics
   - `--stats-only` - Extract stats and exit
4. **Parameter file format** updated with `clipMean`, `clipVar`
5. **Test data** extracted to `data/chr22_test/` (HG002, HG003, HG004)

### Workflow
```bash
# Step 1: Extract WG stats
hmcnclip ref.fa -a full_genome.bam --stats-only -P wg_stats.param

# Step 2: Run on single chromosome with WG stats
hmcnclip chr22.fa -a chr22.bam \
  --wg-mean 35 --wg-var 80 \
  --wg-clip-mean 1 --wg-clip-var 3 \
  -o chr22.vcf
```

---

## Phase 0: Baseline Validation [COMPLETED]

### Results
- Poisson clip model, ~1200 raw calls/sample, results in `results/phase0/`
- Filtered call counts (after exclusion + repeat filter):

| Sample | DEL | DUP | Total |
|--------|-----|-----|-------|
| HG002  | 349 | 573 | 922   |
| HG003  | 386 | 596 | 982   |
| HG004  | 387 | 531 | 918   |

---

## Phase 1: ZINB for Clipping [COMPLETED]

### Mathematical Model
Zero-Inflated Negative Binomial replaces Poisson clip model:

```
P(c | state) = π · I(c=0) + (1-π) · NB(c; μ, φ)
```

### Implementation
- `LgZINB()` added to `src/hmmcnc.cpp`
- π, φ estimated via Method of Moments from chr22 clip data
- `clipPi`, `clipPhi` added to parameter file (backward-compatible)
- Results in `results/phase1/`

**ZINB Parameters (chr22, ~30× coverage):**
- clipPi ≈ 0.997 (99.7% of bins have zero clips)
- clipPhi ≈ 0.29 (highly overdispersed NB component)

**Raw call counts:**

| Sample | DEL | DUP | Total |
|--------|-----|-----|-------|
| HG002  | 435 | 790 | 1225  |
| HG003  | 487 | 832 | 1319  |
| HG004  | 485 | 768 | 1253  |

**Filtered call counts:**

| Sample | DEL | DUP | Total |
|--------|-----|-----|-------|
| HG002  | 346 | 556 | 902   |
| HG003  | 382 | 590 | 972   |
| HG004  | 382 | 511 | 893   |

---

## Phase 1.1: Emission Calibration [COMPLETED]

### Problem
The NB emission model for CN=3 has higher variance (`var ∝ CN`) than CN=1, giving CN=3
a lower peak log-probability. This creates an asymmetric calling penalty:
- `lepsi23_nb ≈ −0.87` nats/bin at the CN=2/3 boundary
- `lepsi21_nb ≈ −1.07` nats/bin at the CN=1/2 boundary

### Changes

**Asymmetric transition initialization (`InitParams`):**
- From CN=2: CN=1 exits use `beta1 = 100 × lepsi21_nb` (more negative)
- From CN=2: CN=3 exits use `beta = 100 × lepsi23_nb` (less negative)
- Calibrates both directions to require ~100 boundary-level bins
- Effect on final calls is negligible (Baum-Welch overwrites the prior)

**Epsi emission penalty (`--epsi-weight`):**
- Adds per-bin NB LLR penalty to non-diploid state emissions
- CN=i offset = `(beta_dir / 100) × epsiWeight × |i − 2|`
- `epsiWeight=1.0` (default): full penalty, ~100-bin threshold
- `epsiWeight=0.0`: disabled (reverts to phase1 behavior)
- Pending precision/recall evaluation before choosing production default

**Call count sweep (HG002, raw):**

| `--epsi-weight` | DEL | DUP | Total |
|----------------|-----|-----|-------|
| 0.0            | 435 | 790 | 1225  |
| 0.25           | 330 | 651 | 981   |
| 0.5            | 252 | 553 | 805   |
| 1.0            | 175 | 353 | 528   |

**Results in `results/phase1.1/` (generated with `--epsi-weight 1.0`)**

**Files modified:**
- `src/hmmcnc.cpp` — asymmetric `InitParams`, epsi emission penalty block
- `include/hmcnc.h` — `epsiWeight` field in `Parameters`

**Pending:** GIAB truth VCFs needed in `benchmarks/truth/` to evaluate precision/recall
and select the production default for `--epsi-weight`.

---

## Phase 2: Directional Clipping (P_L, P_R) [COMPLETED]

### Concept
Separate left-clips and right-clips per bin:
- **Left clips (P_L)**: leading H/S, recorded at alignment `startpos` → right-hand CNV boundary (CN-increasing: j > i)
  - deletion→normal (1→2) and normal→duplication (2→3)
- **Right clips (P_R)**: trailing H/S, recorded at alignment `endpos` → left-hand CNV boundary (CN-decreasing: j < i)
  - normal→deletion (2→1) and duplication→normal (3→2)

### Expected Patterns
| State Transition               | Direction | Expected Clips |
|--------------------------------|-----------|----------------|
| Normal → Deletion  (2→1, CN↓) | j < i     | High P_R       |
| Deletion → Normal  (1→2, CN↑) | j > i     | High P_L       |
| Normal → Duplication (2→3, CN↑)| j > i     | High P_L       |
| Duplication → Normal (3→2, CN↓)| j < i     | High P_R       |

### Implementation

**BAM processing**: `ParseChrom` collects supplementary alignments (flag `BAM_FSUPPLEMENTARY`)
into `suppReads` and processes them separately from primary reads.
- **Supplementary reads** (minimap2/pbmm2 split-read output): leading H/S → `leftClipBins` at
  `pos`; trailing H/S → `rightClipBins` at `endpos`. This is the primary SV signal for long-read
  aligners — at real breakpoints pbmm2/minimap2 emit supplementaries rather than primary soft-clips.
- **Primary reads**: leading S → `leftClipBins`; trailing S → `rightClipBins`. Retained for
  completeness but dominated by pericentromeric noise in minimap2 data.

Both sources feed the same `clipBins`/`leftClipBins`/`rightClipBins` arrays. Secondary alignments
(flag `BAM_FSECONDARY`) are skipped entirely.

**ZINB parameters**: separate MOM estimation for L and R clips via `estimateZINB` lambda.
Logged to stderr as `Left clip — mean/pi/phi` and `Right clip — mean/pi/phi`.

**Directional posteriors**: `PnL/PclL/PnR/PclR` computed per-bin per-contig using L/R ZINB
params. Stored for Phase 4 transition modification — not yet wired into HMM.

**Clip bed format** (7 columns): `chr  start  end  leftClips  rightClips  Pn  Pcl`

**Call counts** (identical to Phase 1.1 — HMM logic unchanged):

| Sample | DEL | DUP | Total |
|--------|-----|-----|-------|
| HG002  | 175 | 353 | 528   |
| HG003  | 160 | 305 | 465   |
| HG004  | 182 | 322 | 504   |

**Files modified:**
- `src/hmmcnc.cpp` — `ThreadInfo` L/R pointers, `ParseChrom` routing, allocation, thread wiring, `estimateZINB` lambda, `PnL/PclL/PnR/PclR` computation
- `src/hmcnc_io.cpp` — `WriteClipBed` updated to write L/R columns
- `include/hmcnc.h` — `WriteClipBed` declaration updated

---

## Phase 3: Newton-Raphson for Dispersion [COMPLETED]

### Algorithm
Estimate dispersion parameter φ using Newton-Raphson (3 iterations) nested in the BW M-step:

```
φ_{k+1} = φ_k - L'(φ_k) / L''(φ_k)
```

Where:
- L'(φ) = Σ_{st,c>0} γ_{st,c} [ψ(φ+c) − ψ(φ) + log(φ/(φ+ν)) + (ν−c)/(φ+ν)]
- L''(φ) = Σ_{st,c>0} γ_{st,c} [ψ'(φ+c) − ψ'(φ) + 1/φ − 1/(φ+ν) − (ν−c)/(φ+ν)²]
- ψ = digamma, ψ' = trigamma (via Boost)

### Implementation

**`UpdateZINBPhi`** (new function): 3-iteration NR, clamped to φ≥0.01.

**ZINB histogram** (`zinbClipHist[state][count]`): γ-weighted per-bin per-count accumulation across all contigs in E-step (inside `ThreadedBWE` mutex).

**M-step block** (in BW outer loop after `BaumWelchM`): updates `clipPi`, `clipMean` from diploid-state histogram, calls `UpdateZINBPhi`, recomputes `Pn/Pcl`.

**NR convergence**: φ converges from MOM prior (~0.29) to ~0.85–1.3 (tighter than MOM because posterior γ concentrates on true breakpoint bins).

**Call counts** (HMM calls unchanged — main coverage model dominates):

| Sample | DEL | DUP | Total |
|--------|-----|-----|-------|
| HG002  | 175 | 352 | 527   |
| HG003  | 160 | 305 | 465   |
| HG004  | 182 | 323 | 505   |

**Files modified:**
- `src/hmmcnc.cpp` — `UpdateZINBPhi`, `zinbClipHist` accumulation in `ThreadedBWE`, ZINB M-step in BW outer loop
- `include/hmcnc.h` — `UpdateZINBPhi` declaration

---

## Phase 4: Directional Transition Modification [COMPLETED]

### Design
At each bin k, the effective log-transition probability for state i→j is:

```
log a_ij(k) = log-sum-exp(covCovTransP[i][j] + pN_ij(k),
                           clipCovCovTransP[i][j] + pCl_ij(k))
```

Where the clip posterior pair is selected by CN direction:
- j > i (CN increasing: deletion→normal, normal→duplication): use PnL[k], PclL[k]
- j < i (CN decreasing: normal→deletion, duplication→normal): use PnR[k], PclR[k]
- j == i (self-transition): use Pn[k], Pcl[k]

Same direction logic applied in the BW expected-transition accumulation.

### Implementation

**New `ForwardBackwards` overload**: takes `Pn, Pcl, PnL, PclL, PnR, PclR`; uses a `clipPair` lambda to select per-(i,j) at each bin.

**`BaumWelchEOnChrom`**: updated signature to accept all 6 posteriors; uses `dirPn`/`dirPcl` lambdas for both F-B call and expected-transition accumulation.

**`ThreadInfo`**: added `nL, clL, nR, clR` pointers.

**`ThreadedBWE`**: passes `(*nL)[curSeq]`, etc. to `BaumWelchEOnChrom`.

**Call counts** (directional posteriors are subtle — coverage still dominates):

| Sample | DEL | DUP | Total |
|--------|-----|-----|-------|
| HG002  | 175 | 353 | 528   |
| HG003  | 160 | 306 | 466   |
| HG004  | 183 | 322 | 505   |

**Files modified:**
- `src/hmmcnc.cpp` — new directional `ForwardBackwards`, updated `BaumWelchEOnChrom`, `ThreadInfo`, `ThreadedBWE`
- `include/hmcnc.h` — new `ForwardBackwards` and `BaumWelchEOnChrom` declarations

---

## Phase 5: Directional ZINB M-step [COMPLETED]

### Design
Extend the Phase 3 combined ZINB M-step to update directional ZINB parameters
(πL, μL, φL, πR, μR, φR) separately using γ-weighted directional clip histograms.
These parameters were MOM-initialized in Phase 2 and held fixed through Phase 4.

The spec (`bw_modified.pdf`) adds clip evidence to the Q-function and calls for
EM updates of πL, πR, ϵL, ϵR (zero-inflation and clip rates) in the M-step.

### New E-step accumulators
`zinbLeftHist[nStates][count]` and `zinbRightHist[nStates][count]` accumulate
γ_t(i)-weighted left/right clip counts per state alongside the existing `zinbClipHist`.

### New M-step block
After the combined ZINB update, an `updateDirZINB` lambda runs NR (via `UpdateZINBPhi`)
for each direction, updates π and mean from the diploid-state histogram, and recomputes
PnL/PclL/PnR/PclR with the refreshed parameters for the next BW iteration.

### HG002 chr22 convergence (2 BW iterations before early-stop)
| Iter | clipPi | clipMean | clipPhi | πL | μL | φL | πR | μR | φR |
|------|--------|----------|---------|-----|----|----|-----|----|----|
| 1 | 0.9986 | 3.37 | 0.845 | 0.9993 | 2.92 | 1.043 | 0.9993 | 3.44 | 0.829 |
| 2 | 0.9987 | 2.98 | 0.985 | 0.9993 | 2.61 | 1.188 | 0.9993 | 3.04 | 0.998 |

### Call counts (effectively unchanged — coverage still dominates)

| Sample | DEL | DUP | Total |
|--------|-----|-----|-------|
| HG002  | 175 | 353 | 528   |
| HG003  | 160 | 306 | 466   |
| HG004  | 183 | 322 | 505   |

**Files modified:**
- `src/hmmcnc.cpp` — `ThreadInfo` L/R hist pointers, E-step directional histogram accumulation, directional ZINB M-step + PnL/PclL/PnR/PclR recompute

---

## Phase 6: Integration & Validation

### Testing
- [ ] End-to-end testing on HG002/HG003/HG004 chr22
- [ ] Full genome testing
- [ ] Comparison against truth sets (GIAB CNV calls)
- [ ] Performance benchmarking

### Validation Metrics
- Sensitivity (recall) for deletions/duplications
- Precision (PPV)
- Breakpoint accuracy (within N bp)
- Runtime and memory usage

### Documentation
- [ ] Update CLAUDE.md with new features
- [ ] Document parameter tuning guidelines
- [ ] Create example workflows

---

## Reference: PDF Specifications

Mathematical specifications are in `definitions/`:
- `cnv.pdf` - Core HMM with ZINB and directional clipping
- `bw_modified.pdf` - Baum-Welch with clip evidence
- `cnv (dispersion).pdf` - Dispersion parameter estimation
- `cnv(ZINB_disp).pdf` - Newton-Raphson for ZINB dispersion
- `cnv (BW and NR).pdf` - Nested optimization guide
