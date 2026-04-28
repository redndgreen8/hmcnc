# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

HMCNC (Hidden Markov Copy Number Caller) is a bioinformatics tool for calling copy number variations (CNVs) from sequencing data. It uses Hidden Markov Models with read depth analysis and optional SNV/clipping signatures.

## Documentation

- `docs/IMPLEMENTATION_ROADMAP.md` - Phased implementation plan for ZINB/directional clipping
- `docs/TECHNICAL_REFERENCE.md` - Codebase architecture and API reference
- `docs/MATHEMATICAL_MODELS.md` - HMM, ZINB, and algorithm specifications
- `definitions/*.pdf` - Mathematical specification documents

## Build Commands

### Makefile (Primary)
```bash
# Activate conda environment
conda activate hmcnc

# Build all targets
cd src && make

# Outputs: hmcnclip, viterbi, samToBed
```

### Meson (Alternative)
```bash
cd hmcnc && mkdir build && cd build && meson .. && ninja
```

**Dependencies:** boost, htslib (>=1.4), gxx, pthread

### Conda Environment Setup
```bash
conda create -n hmcnc
conda activate hmcnc
conda install -c conda-forge gxx_osx-arm64 htslib boost samtools bedtools tabix
```

## Running the Tool

### Basic Usage
```bash
# From BAM file
hmcnclip reference.fa -a alignments.bam -o output.vcf

# From pre-computed coverage bed
hmcnclip reference.fa -b coverage.bed -s snvs.txt -o output.vcf
```

### Single-Chromosome Testing with WG Stats
```bash
# Step 1: Extract stats from full genome
hmcnclip ref.fa -a full_genome.bam --stats-only -P wg_stats.param

# Step 2: Run on single chromosome with those stats
hmcnclip chr22.fa -a chr22.bam \
  --wg-mean 35 --wg-var 80 \
  --wg-clip-mean 1 --wg-clip-var 3 \
  -o chr22.vcf
```

### Key CLI Options
```
Input:
  -a FILE    BAM file (calculate depth on the fly)
  -b FILE    Pre-computed coverage bed
  -s FILE    SNV file
  -p FILE    Parameter file (skip training)
  -l FILE    Clipping signature file

Statistics Override:
  --wg-mean FLOAT       Haploid mean coverage
  --wg-var FLOAT        Coverage variance
  --wg-clip-mean FLOAT  Mean clip count per bin
  --wg-clip-var FLOAT   Clip variance
  --stats-only          Output stats and exit

Depth Calculation:
  --epsi-weight FLOAT   Emission penalty weight [1.0]. Per-bin NB LLR penalty
                        applied to non-diploid states. 1.0 = 100-bin threshold,
                        0.5 = 200-bin threshold, 0.0 = disabled.

Output:
  -o FILE    VCF output
  -P FILE    Trained parameter file
  -B FILE    Coverage bed output
  --bed FILE BED format output
```

## Running Tests

### Test Data
```
tests/data/           - chr22 reference and coverage bed
data/chr22_test/      - Extracted chr22 BAMs from HG002/HG003/HG004
```

### Test Script
```bash
./tests/test_single_chrom.sh
```

### Meson Tests (if using meson build)
```bash
meson test -v
```

## Architecture

### Core Algorithm Flow

1. **Coverage calculation** - Multi-threaded BAM processing into 100bp bins
2. **SNV extraction** (optional) - Variant calling from alignments
3. **Clipping analysis** (optional) - Structural variant signatures from read clipping
4. **HMM training** - Baum-Welch (EM) algorithm for parameter estimation
5. **HMM inference** - Viterbi algorithm for copy number state sequence
6. **Output** - VCF generation with optional parameter/coverage files

### Key Source Files

- `src/hmmcnc.cpp` - Main algorithm: Forward-Backward, Baum-Welch, emission calculations, BAM processing (~2700 lines)
- `src/hmcnc_io.cpp` - File I/O for BED, FAI, parameter files, SNVs
- `src/main.cpp` - CLI entry point
- `src/viterbi.cpp` - Standalone Viterbi decoder
- `src/SamToBed.cpp` - BAM to coverage bed converter
- `include/hmcnc.h` - Core data structures (SNV, Interval, Parameters)

### Coverage Models

- **Negative Binomial** (default, `-m nb`) - Better handles overdispersion
- **Poisson** (`-m pois`) - Simpler model

### Constants (in hmmcnc.cpp)

- `BIN_LENGTH = 100` - Coverage bin size in bases
- `MIN_CLIP_LENGTH = 500` - Minimum clipping to flag structural variants
- `MIN_MAPQ = 10` - Mapping quality threshold
- `MAX_CN = 6` - Maximum copy number states (0-6)

### Python Implementations

`src/hmcnc_3.py` - ZINB and directional clipping prototype (reference for C++ implementation)

### Workflow Orchestration

`HMM/caller.smk` - Snakemake pipeline for full CNV calling workflow with annotation

## Current Development Status

**Complete:**
- Pre-Phase: testing infrastructure, CLI extensions, test data (HG002/HG003/HG004 chr22)
- Phase 0: Baseline validation — Poisson clip model, results in `results/phase0/`
- Phase 1: ZINB clipping — `LgZINB()` replaces Poisson; pi/phi via MOM; results in `results/phase1/`
- Phase 1.1: Emission calibration — asymmetric transition init + `--epsi-weight` penalty
- Phase 2: Directional clipping — separate L/R clip bins + ZINB params; results in `results/phase2/`
- Phase 3: Newton-Raphson dispersion — `UpdateZINBPhi` in BW M-step; results in `results/phase3/`
- Phase 4: Directional transition modification — direction-aware clip posteriors in F-B and BW E-step; results in `results/phase4/`

**Phase ZINB Parameters (BW-updated, chr22):**
- MOM init: clipPi ≈ 0.997, clipPhi ≈ 0.29
- After NR in BW loop: clipPhi converges to ~0.85–1.3 (tighter, γ-weighted)
- Left clip: pi ≈ 0.999, phi ≈ 0.37; Right clip: pi ≈ 0.999, phi ≈ 0.26 (MOM estimates)

**Directional Convention:**
- Left clips (frontClip, S at CIGAR start): CN-decreasing transitions (j < i)
- Right clips (backClip, S at CIGAR end): CN-increasing transitions (j > i)
- Self-transitions (j == i): combined Pn/Pcl

**Benchmarking Pipeline:**
```bash
# Run samples for a given phase (no truth VCFs needed):
snakemake -s benchmarks/benchmark.smk --configfile benchmarks/config.yaml \
  --config phase=phase4 -j3 calls

# With custom epsi-weight:
snakemake -s benchmarks/benchmark.smk --configfile benchmarks/config.yaml \
  --config phase=phase4 extra_args="--epsi-weight 0.5" -j3 filter

# Full benchmark with truvari (needs benchmarks/truth/<SAMPLE>.chr22.truth.vcf.gz):
snakemake -s benchmarks/benchmark.smk --configfile benchmarks/config.yaml \
  --config phase=phase4 -j3 all
```

**Note:** hg38 chr22 reference is at `/Users/red/repos/MethSmoothEval/data/annotations/hg38.chr22.fa`.
The symlinks in `data/chr22_test/` are broken. Reference goes **last** in the command line.

**Next Phase:** Phase 5 — dual forward-backward passes with clip-evidence posteriors (`bw_modified.pdf`)

See `docs/IMPLEMENTATION_ROADMAP.md` for full roadmap.
