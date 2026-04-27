# HMCNC Technical Reference

## Architecture Overview

### Core Components

```
include/
└── hmcnc.h          # Data structures, function declarations

src/
├── main.cpp         # CLI entry point
├── hmmcnc.cpp       # Core HMM algorithms (2700+ lines)
├── hmcnc_io.cpp     # Parameter/data I/O
├── viterbi.cpp      # Standalone Viterbi decoder
└── SamToBed.cpp     # BAM to coverage bed converter
```

### Data Flow

```
BAM file
    ↓
Coverage calculation (100bp bins)
    ↓
Clip signature extraction
    ↓
SNV extraction (optional)
    ↓
HMM parameter estimation (Baum-Welch)
    ↓
Viterbi decoding → Copy number states
    ↓
Interval merging → CNV calls
    ↓
VCF output
```

---

## Key Data Structures

### Parameters (hmcnc.h:46-82)
```cpp
struct Parameters {
  // Input files
  std::string bamFileName;
  std::string covBedInFileName;
  std::string clipInFileName;
  std::string paramInFile;

  // Output files
  std::string outFileName;
  std::string paramOutFile;

  // WG stats overrides (for single-chrom testing)
  double wgMean = -1;       // -1 means estimate from data
  double wgVar = -1;
  double wgClipMean = -1;
  double wgClipVar = -1;
  bool statsOnly = false;

  // Model settings
  MODEL_TYPE model = NEG_BINOM;  // POIS or NEG_BINOM
  int nproc = 4;
};
```

### Interval (hmcnc.h:15-30)
```cpp
struct Interval {
  int start;
  int end;
  int copyNumber;
  float averageCoverage;
  double pVal;
  std::string filter;
  int distanceToFrontClip;
  int distanceToEndClip;
  int nFrontClip;
  int nEndClip;
};
```

### SNV (hmcnc.h:32-44)
```cpp
struct SNV {
  int pos;
  char refNuc;
  char altNuc;
  int ref;  // reference allele count
  int alt;  // alternate allele count
};
```

---

## HMM Model

### States
7 copy number states (s0-s6) representing CN 0-6:
- s0: Homozygous deletion (CN=0)
- s1: Heterozygous deletion (CN=1)
- s2: Normal diploid (CN=2)
- s3: Single copy gain (CN=3)
- s4-s6: Higher amplifications

### Emissions

**Coverage (Negative Binomial)**
```cpp
double LgNegBinom(int cn, int cov, float Hmean, float Hvar);
```
- Mean scales with copy number: `μ_cn = Hmean × cn`
- Variance from whole-genome estimate

**Clipping (current simple model)**
- Binary: clip present or absent
- Plans: ZINB with directional separation

### Transitions

Two transition matrices:
1. `covCovTransP` - Base transitions (neutral)
2. `clipCovCovTransP` - Modified transitions (with clip evidence)

Current implementation uses high self-transition probability with small transition rates between adjacent states.

---

## Core Algorithms

### Forward-Backward (hmcnc.h:125-139)
```cpp
double ForwardBackwards(
    const std::vector<double> &startP,
    const std::vector<std::vector<double>> &covCovTransP,
    const std::vector<std::vector<double>> &clipCovCovTransP,
    const std::vector<std::vector<double>> &emisP,
    const std::vector<int> &obs,
    std::vector<std::vector<double>> &f,  // forward probs
    std::vector<std::vector<double>> &b,  // backward probs
    std::vector<double> &Pn,              // neutral posterior
    std::vector<double> &Pcl              // clipped posterior
);
```

### Baum-Welch E-Step (hmcnc.h:88-97)
```cpp
double BaumWelchEOnChrom(
    const std::vector<double> &startP,
    std::vector<std::vector<double>> &covCovTransP,
    std::vector<std::vector<double>> &clipCovCovTransP,
    std::vector<std::vector<double>> &emisP,
    std::vector<int> &obs,
    std::vector<std::vector<double>> &f,
    std::vector<std::vector<double>> &b,
    std::vector<std::vector<double>> &expCovCovTransP,
    std::vector<std::vector<double>> &expEmisP,
    std::vector<double> &Pn,
    std::vector<double> &Pcl
);
```

### Baum-Welch M-Step (hmcnc.h:99-110)
```cpp
void BaumWelchM(
    const std::vector<double> &startP,
    const std::vector<std::vector<double>> &transP,
    const std::vector<std::vector<double>> &emisP,
    const std::vector<std::vector<std::vector<double>>> &binoP,
    int model,
    const std::vector<long> &stateTotCov,
    const std::vector<long> &stateNCov,
    const std::vector<std::vector<double>> &expTransP,
    std::vector<std::vector<double>> &expEmisP,
    std::vector<std::vector<double>> &covCovPrior,
    std::vector<std::vector<double>> &updateTransP,
    std::vector<std::vector<double>> &updateEmisP
);
```

### Viterbi Decoding
Implemented in `StorePosteriorMaxIntervals()` - finds most likely state sequence and extracts intervals.

---

## Parameter File Format

```
nStates    7
covMean    35
covVar     80
clipMean   1
clipVar    3
maxState   6
maxCov     500
startP
-1.94591
-1.94591
...
transP     7  7
[7x7 log-probability matrix]
clipTransP 7  7
[7x7 log-probability matrix]
emisP      7  500
[7x500 log-probability matrix]
```

All probabilities stored in log space.

---

## I/O Functions

### Reading
```cpp
void ReadFai(filename, contigNames, contigLengths);
void ReadCoverage(filename, contigNames, covBins);
void ReadSNVs(filename, contigNames, snvs);
void ReadParameterFile(filename, nStates, covMean, covVar,
                       clipMean, clipVar, maxState, maxCov,
                       startP, transP, clipTransP, emisP);
```

### Writing
```cpp
void WriteCovBed(filename, contigNames, covBins);
void WriteClipBed(filename, contigNames, covBins, Pn, Pcl);
void WriteParameterFile(filename, nStates, covMean, covVar,
                        clipMean, clipVar, maxState, maxCov,
                        startP, transP, clipTransP, emisP);
void WriteVCF(out, refName, sampleName, contigNames,
              contigLengths, intervals, writeFail);
void WriteBed(intervals, filename, contigNames);
```

---

## Build System

### Makefile
```makefile
all: samToBed viterbi hmcnclip

hmcnclip: main.cpp hmmcnc.o hmcnc_io.o
    g++ -g2 -I $(CONDA_PREFIX)/include $^ -o $@ \
        -L $(CONDA_PREFIX)/lib -lhts -lpthread -ffast-math \
        -Wl,-rpath,$(CONDA_PREFIX)/lib
```

### Dependencies
- htslib (BAM/CRAM reading)
- Boost (string algorithms)
- pthread

### Conda Environment
```bash
conda create -n hmcnc
conda activate hmcnc
conda install -c conda-forge \
    gxx_osx-arm64 htslib boost samtools bedtools
```

---

## CLI Reference

```
hmcnclip [OPTIONS] reference

Positionals:
  reference FILE REQUIRED     FASTA reference file

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

Output:
  -o FILE    VCF output
  -B FILE    Coverage bed output
  -P FILE    Trained parameter file
  -L FILE    Clip signature output
  --bed FILE BED format output
  --sample   Sample name in VCF

Options:
  -t INT     Threads (default: 4)
  -m TEXT    Model: pois or nb (default: nb)
  -C TEXT    Run on single chromosome
  -M         Merge consecutive bins
```

---

## Log Probability Utilities

All HMM computations use log probabilities for numerical stability:

```cpp
double PairSumOfLogP(double a, double b);  // log(exp(a) + exp(b))
double SumOfLogP(const std::vector<double> &vals);  // log-sum-exp
```
