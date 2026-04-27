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

## Phase 2: Directional Clipping (P_L, P_R)

### Concept
Separate left-clips and right-clips per bin:
- **Left clips (P_L)**: Reads clipped on left side → deletion START, duplication END
- **Right clips (P_R)**: Reads clipped on right side → deletion END, duplication START

### Expected Patterns
| State Transition | Expected Clips |
|-----------------|----------------|
| Normal → Deletion | High P_L |
| Deletion → Normal | High P_R |
| Normal → Duplication | High P_R |
| Duplication → Normal | High P_L |

### Implementation Tasks
- [ ] Modify clip counting in BAM processing to separate L/R
- [ ] Add `leftClipBins`, `rightClipBins` vectors
- [ ] Separate ZINB parameters for each direction
- [ ] Update clip bed I/O format

### Files to Modify
- `src/hmmcnc.cpp` - Separate clip counting
- `src/hmcnc_io.cpp` - Directional clip I/O
- `include/hmcnc.h` - Data structure updates

---

## Phase 3: Newton-Raphson for Dispersion

### Algorithm
Estimate dispersion parameter φ using Newton-Raphson iteration:

```
φ_{k+1} = φ_k - f(φ_k) / f'(φ_k)
```

Where f(φ) is derived from the ZINB log-likelihood gradient.

### Implementation Tasks
- [ ] Implement digamma/trigamma functions (or use Boost)
- [ ] Add NR iteration for φ estimation
- [ ] Set convergence criteria (tolerance, max iterations)
- [ ] Nest within Baum-Welch M-step

### Files to Modify
- `src/hmmcnc.cpp` - Add NR optimization
- `include/hmcnc.h` - Add helper function declarations

---

## Phase 4: Transition Modification with Clipping

### Formula
Modify transition probabilities based on clipping evidence:

```
a_{t,ij} = a_{b,ij} × (1 + α × m_ij(P_L, P_R))
```

Where:
- `a_{b,ij}` = base transition probability
- `α` = clipping influence weight
- `m_ij(P_L, P_R)` = direction-aware modifier function

### Modifier Function
```
m_ij(P_L, P_R) = {
  P_L - P_R  if transition expects left clips
  P_R - P_L  if transition expects right clips
  0          otherwise
}
```

### Implementation Tasks
- [ ] Define state transition expectations matrix
- [ ] Implement `m_ij()` function
- [ ] Add per-position transition modification
- [ ] Tune α parameter

### Files to Modify
- `src/hmmcnc.cpp` - Transition modification logic

---

## Phase 5: Baum-Welch with Clip Evidence

### Dual Forward-Backward
Run two forward-backward passes:
1. **Neutral**: Using base transitions `a_{b,ij}`
2. **Clipped**: Using modified transitions `a_{t,ij}`

### Posterior Combination
Combine posteriors with position-dependent weights:

```
γ_final(t,i) = w_n(t) × γ_neutral(t,i) + w_c(t) × γ_clipped(t,i)
```

Where weights depend on local clipping evidence.

### Implementation Tasks
- [ ] Implement dual forward-backward passes
- [ ] Add posterior combination logic
- [ ] Update expected sufficient statistics calculation
- [ ] Modify `BaumWelchEOnChrom()` and `BaumWelchM()`

### Files to Modify
- `src/hmmcnc.cpp` - Core algorithm changes
- `include/hmcnc.h` - Updated function signatures

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
