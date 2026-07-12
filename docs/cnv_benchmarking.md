# CNV Caller Benchmarking Guide

## Overview

Benchmark a CNV/SV caller against two orthogonal truth sets:
- **Platinum Pedigree** (NA12878, CEPH-1463) — pedigree-validated, long-read truth set
- **GIAB** (HG002/HG003/HG004, Ashkenazim trio) — multi-technology integrated truth set
- **Cancer GIAB** (HG008, Pancreatic Cancer) — somatic structural variant and CNV benchmarks

Reference genome: **GRCh38**

---

## Dependencies

```bash
# Install required tools
conda create -n cnv-bench python=3.10
conda activate cnv-bench

conda install -c bioconda -c conda-forge \
  truvari \
  bcftools \
  samtools \
  hap.py \
  bedtools \
  awscli \
  curl \
  wget
```

---

## 1. Download Truth Sets

### 1a. Platinum Pedigree (NA12878)

```bash
mkdir -p truth/platinum_pedigree

# SV truth set (use for CNV benchmarking)
aws s3 cp --no-sign-request \
  s3://platinum-pedigree-data/variants/merged_sv_truthset/GRCh38/merged_hg38.svs.sort.oa.vcf.gz \
  truth/platinum_pedigree/

aws s3 cp --no-sign-request \
  s3://platinum-pedigree-data/variants/merged_sv_truthset/GRCh38/merged_hg38.svs.sort.oa.vcf.gz.tbi \
  truth/platinum_pedigree/

# SV truth set with tandem repeats excluded (cleaner for CNV benchmarking)
aws s3 cp --no-sign-request \
  s3://platinum-pedigree-data/variants/merged_sv_truthset/GRCh38/merged_hg38.svs.TRexclusion.sort.oa.vcf.gz \
  truth/platinum_pedigree/

aws s3 cp --no-sign-request \
  s3://platinum-pedigree-data/variants/merged_sv_truthset/GRCh38/merged_hg38.svs.TRexclusion.sort.oa.vcf.gz.tbi \
  truth/platinum_pedigree/

# High-confidence regions BED
aws s3 cp --no-sign-request \
  s3://platinum-pedigree-data/variants/small_variant_truthset/GRCh38/hq_regions_final.bed.gz \
  truth/platinum_pedigree/

# Latest versioned truthset (v1.2 recommended)
aws s3 cp --no-sign-request --recursive \
  s3://platinum-pedigree-data/truthset_v1.2/ \
  truth/platinum_pedigree/truthset_v1.2/
```

### 1b. GIAB — HG002 SV Truth Set (v0.6, GRCh37)

> Note: GIAB SV v0.6 is GRCh37. Liftover to GRCh38 if needed (see Section 4).

```bash
mkdir -p truth/giab/hg002

FTPDIR=ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6

curl -O --output-dir truth/giab/hg002 ${FTPDIR}/HG002_SVs_Tier1_v0.6.vcf.gz
curl -O --output-dir truth/giab/hg002 ${FTPDIR}/HG002_SVs_Tier1_v0.6.vcf.gz.tbi
curl -O --output-dir truth/giab/hg002 ${FTPDIR}/HG002_SVs_Tier1_v0.6.bed
```

### 1c. GIAB — HG002/3/4 Small Variant Benchmarks (v4.2.1, GRCh38)

```bash
BASE=https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio

# HG002
mkdir -p truth/giab/hg002
curl -O --output-dir truth/giab/hg002 \
  ${BASE}/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
curl -O --output-dir truth/giab/hg002 \
  ${BASE}/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
curl -O --output-dir truth/giab/hg002 \
  ${BASE}/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed

# HG003
mkdir -p truth/giab/hg003
curl -O --output-dir truth/giab/hg003 \
  ${BASE}/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
curl -O --output-dir truth/giab/hg003 \
  ${BASE}/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
curl -O --output-dir truth/giab/hg003 \
  ${BASE}/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed

# HG004
mkdir -p truth/giab/hg004
curl -O --output-dir truth/giab/hg004 \
  ${BASE}/HG004_NA24143_mother/NISTv4.2.1/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
curl -O --output-dir truth/giab/hg004 \
  ${BASE}/HG004_NA24143_mother/NISTv4.2.1/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
curl -O --output-dir truth/giab/hg004 \
  ${BASE}/HG004_NA24143_mother/NISTv4.2.1/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
```

### 1d. Cancer GIAB — HG008 Somatic SV/CNV Benchmarks (v0.5)

```bash
mkdir -p truth/giab/hg008

FTPDIR=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/NIST_HG008-T_somatic-stvar-CNV_DraftBenchmark_V0.5-20260318

# Download somatic SV/CNV draft benchmark V0.5
curl -O --output-dir truth/giab/hg008 ${FTPDIR}/HG008-T_somatic-stvar-CNV_DraftBenchmark_V0.5.vcf.gz
curl -O --output-dir truth/giab/hg008 ${FTPDIR}/HG008-T_somatic-stvar-CNV_DraftBenchmark_V0.5.vcf.gz.tbi
```

---

## 2. Download BAMs

### 2a. Platinum Pedigree BAMs (HiFi, GRCh38)

```bash
mkdir -p bams/platinum_pedigree

# NA12878 cell line — held out from truth set, recommended for benchmarking
aws s3 cp --no-sign-request \
  s3://platinum-pedigree-data/data/hifi/mapped/GRCh38/NA12878-cell-line-revio.GRCh38.haplotagged.bam \
  bams/platinum_pedigree/

aws s3 cp --no-sign-request \
  s3://platinum-pedigree-data/data/hifi/mapped/GRCh38/NA12878-cell-line-revio.GRCh38.haplotagged.bam.bai \
  bams/platinum_pedigree/

# List all available HiFi BAMs across the family
aws s3 ls --no-sign-request s3://platinum-pedigree-data/data/hifi/mapped/GRCh38/

# Illumina and ONT BAMs also available
aws s3 ls --no-sign-request s3://platinum-pedigree-data/data/illumina/
aws s3 ls --no-sign-request s3://platinum-pedigree-data/data/ont/
```

### 2b. GIAB BAMs (HiFi, GRCh38)

```bash
mkdir -p bams/giab

BASE=ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio

# HG002 — HiFi CCS (15kb+20kb chemistry2, GRCh38)
curl -O --output-dir bams/giab \
  ${BASE}/HG002_NA24385_son/PacBio_CCS_15kb_20kb_chemistry2/GRCh38/HG002.GRCh38.pbmm2.bam
curl -O --output-dir bams/giab \
  ${BASE}/HG002_NA24385_son/PacBio_CCS_15kb_20kb_chemistry2/GRCh38/HG002.GRCh38.pbmm2.bam.bai

# HG003 / HG004 HiFi — check FTP for latest aligned BAMs
# ${BASE}/HG003_NA24149_father/
# ${BASE}/HG004_NA24143_mother/

# Latest Revio data (2024Q4, all three samples)
# https://downloads.pacbcloud.com/public/revio/2024Q4/WGS/GIAB_trio/
```

---

## 3. Prepare Truth VCFs — Filter to CNV Types

CNVs are a subset of SVs. Filter truth sets to `DEL` and `DUP` only:

```bash
# Platinum Pedigree — CNV-only truth VCF
bcftools view -i 'SVTYPE="DEL" || SVTYPE="DUP"' \
  truth/platinum_pedigree/merged_hg38.svs.sort.oa.vcf.gz \
  -O z -o truth/platinum_pedigree/NA12878_CNV_truth.vcf.gz
bcftools index -t truth/platinum_pedigree/NA12878_CNV_truth.vcf.gz

# GIAB HG002 — CNV-only truth VCF
bcftools view -i 'SVTYPE="DEL" || SVTYPE="DUP"' \
  truth/giab/hg002/HG002_SVs_Tier1_v0.6.vcf.gz \
  -O z -o truth/giab/hg002/HG002_CNV_truth.vcf.gz
bcftools index -t truth/giab/hg002/HG002_CNV_truth.vcf.gz
```

---

## 4. (Optional) Liftover GIAB v0.6 to GRCh38

```bash
# Download chain file and GRCh37 reference if needed
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# Use bcftools or CrossMap
pip install CrossMap
CrossMap.py vcf hg19ToHg38.over.chain.gz \
  truth/giab/hg002/HG002_CNV_truth.vcf.gz \
  ref/GRCh38.fa \
  truth/giab/hg002/HG002_CNV_truth_GRCh38.vcf
bgzip truth/giab/hg002/HG002_CNV_truth_GRCh38.vcf
bcftools index -t truth/giab/hg002/HG002_CNV_truth_GRCh38.vcf.gz
```

---

## 5. Run Your CNV Caller

Replace the block below with your caller's actual command:

```bash
# Example — replace with your caller
mkdir -p calls/

# Platinum Pedigree / NA12878
your_cnv_caller \
  --bam bams/platinum_pedigree/NA12878-cell-line-revio.GRCh38.haplotagged.bam \
  --ref ref/GRCh38.fa \
  --out calls/NA12878_calls.vcf.gz

# GIAB / HG002
your_cnv_caller \
  --bam bams/giab/HG002.GRCh38.pbmm2.bam \
  --ref ref/GRCh38.fa \
  --out calls/HG002_calls.vcf.gz
```

---

## 6. Benchmark with Truvari

Truvari is the standard tool for SV/CNV benchmarking. Use `--pick multi` for CNV callers that may report overlapping events.

```bash
mkdir -p results/

# --- Platinum Pedigree / NA12878 ---
truvari bench \
  -b truth/platinum_pedigree/NA12878_CNV_truth.vcf.gz \
  --includebed truth/platinum_pedigree/hq_regions_final.bed.gz \
  -c calls/NA12878_calls.vcf.gz \
  -o results/NA12878_truvari/ \
  --passonly \
  --pick multi \
  -r 500 \
  -p 0.0 \
  --svtype DEL DUP

# --- GIAB / HG002 ---
truvari bench \
  -b truth/giab/hg002/HG002_CNV_truth.vcf.gz \
  --includebed truth/giab/hg002/HG002_SVs_Tier1_v0.6.bed \
  -c calls/HG002_calls.vcf.gz \
  -o results/HG002_truvari/ \
  --passonly \
  --pick multi \
  -r 500 \
  -p 0.0 \
  --svtype DEL DUP
```

### Key Truvari Parameters

| Parameter | Value | Notes |
|-----------|-------|-------|
| `-r` | `500` | Max breakend distance (bp). Increase to 1000 for less precise callers |
| `-p` | `0.0` | Size similarity threshold. Set to `0.7` for stricter matching |
| `--passonly` | flag | Only evaluate PASS calls. Remove if your caller doesn't set FILTER |
| `--pick` | `multi` | Allow multiple base matches per call. Use `ac` for genotype-aware |
| `--svtype` | `DEL DUP` | Restrict to CNV types only |

---

## 7. Parse and Summarise Results

```bash
# Print summary for each benchmark
echo "=== Platinum Pedigree / NA12878 ===" 
cat results/NA12878_truvari/summary.json | python3 -m json.tool

echo "=== GIAB / HG002 ==="
cat results/HG002_truvari/summary.json | python3 -m json.tool

# Quick precision/recall/F1 table
python3 - <<'EOF'
import json, glob

results = {
    "NA12878 (Platinum Pedigree)": "results/NA12878_truvari/summary.json",
    "HG002 (GIAB)"               : "results/HG002_truvari/summary.json",
}

print(f"{'Sample':<35} {'Precision':>10} {'Recall':>10} {'F1':>10} {'TP':>8} {'FP':>8} {'FN':>8}")
print("-" * 95)
for label, path in results.items():
    with open(path) as f:
        s = json.load(f)
    print(f"{label:<35} {s['precision']:>10.4f} {s['recall']:>10.4f} {s['f1']:>10.4f} "
          f"{s['TP-call']:>8} {s['FP']:>8} {s['FN']:>8}")
EOF
```

---

## 8. Stratify by SV Size

```bash
# Split truth and calls by size bin, then re-benchmark
for SIZE_MIN in 50 500 5000; do
  SIZE_MAX=$((SIZE_MIN * 10))
  LABEL="${SIZE_MIN}-${SIZE_MAX}bp"

  bcftools view -i "ABS(SVLEN)>=${SIZE_MIN} && ABS(SVLEN)<${SIZE_MAX}" \
    truth/platinum_pedigree/NA12878_CNV_truth.vcf.gz \
    -O z -o truth/platinum_pedigree/NA12878_CNV_truth_${LABEL}.vcf.gz
  bcftools index -t truth/platinum_pedigree/NA12878_CNV_truth_${LABEL}.vcf.gz

  truvari bench \
    -b truth/platinum_pedigree/NA12878_CNV_truth_${LABEL}.vcf.gz \
    --includebed truth/platinum_pedigree/hq_regions_final.bed.gz \
    -c calls/NA12878_calls.vcf.gz \
    -o results/NA12878_${LABEL}/ \
    --passonly --pick multi -r 500 -p 0.0
done
```

---

## 9. Stratify by Genomic Region (Optional)

```bash
# Evaluate performance inside segmental duplications (harder regions)
wget -O segdups_GRCh38.bed.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz

# Intersect truth with segdups
bedtools intersect -a truth/platinum_pedigree/NA12878_CNV_truth.vcf.gz \
  -b segdups_GRCh38.bed.gz -header \
  | bgzip > truth/platinum_pedigree/NA12878_CNV_truth_segdups.vcf.gz
bcftools index -t truth/platinum_pedigree/NA12878_CNV_truth_segdups.vcf.gz

truvari bench \
  -b truth/platinum_pedigree/NA12878_CNV_truth_segdups.vcf.gz \
  -c calls/NA12878_calls.vcf.gz \
  -o results/NA12878_segdups/ \
  --passonly --pick multi -r 500 -p 0.0
```

---

## Directory Structure

```
cnv-benchmark/
├── truth/
│   ├── platinum_pedigree/
│   │   ├── merged_hg38.svs.sort.oa.vcf.gz
│   │   ├── merged_hg38.svs.TRexclusion.sort.oa.vcf.gz
│   │   ├── NA12878_CNV_truth.vcf.gz          # filtered DEL/DUP
│   │   ├── hq_regions_final.bed.gz
│   │   └── truthset_v1.2/
│   └── giab/
│       ├── hg002/
│       │   ├── HG002_SVs_Tier1_v0.6.vcf.gz
│       │   ├── HG002_SVs_Tier1_v0.6.bed
│       │   └── HG002_CNV_truth.vcf.gz        # filtered DEL/DUP
│       ├── hg003/
│       ├── hg004/
│       └── hg008/
├── bams/
│   ├── platinum_pedigree/
│   │   └── NA12878-cell-line-revio.GRCh38.haplotagged.bam
│   └── giab/
│       └── HG002.GRCh38.pbmm2.bam
├── calls/
│   ├── NA12878_calls.vcf.gz
│   └── HG002_calls.vcf.gz
├── results/
│   ├── NA12878_truvari/
│   └── HG002_truvari/
└── ref/
    └── GRCh38.fa
```

---

## References

- Kronenberg et al. (2025). *The Platinum Pedigree: a long-read benchmark for genetic variants.* Nature Methods 22, 1669–1676. https://doi.org/10.1038/s41592-025-02750-y
- Zook et al. (2020). *A robust benchmark for germline structural variant detection.* Nature Biotechnology. https://doi.org/10.1038/s41587-020-0538-8
- Zook et al. (2021). *Benchmarking challenging small variants with linked and long reads.* Cell Genomics.
- Genome in a Bottle Consortium (2026). *A complete human pancreatic cancer genome.* bioRxiv. https://doi.org/10.64898/2026.05.01.722316
- Data: `s3://platinum-pedigree-data/` (open, no credentials needed)
- Data: `ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/`
- Data: NHGRI AnVIL Project (for HGSVC, T2T/CHM13, and cancer long reads) - `https://anvilproject.org/`
- Data: The Cancer Genome Atlas (TCGA) - `https://gdc.cancer.gov/`
- Data: Human Pangenome Reference Consortium (HPRC) - `https://humanpangenome.org/`
- GitHub: https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Datasets

### Additional Datasets and Open Registries
- **1KG-ONT-VIENNA panel:** Medium coverage ONT sequencing data for 1,019 samples from the 1000 Genomes Project collection, including SVs and haplotype context.
- **Genome in a Bottle on AWS (`giab`):** Full NIST reference genomes mirrored on AWS.
- **1000 Genomes Phase 3 Reanalysis with DRAGEN (`ilmn-dragen-1kgp`):** Includes high-confidence small variant, CNV, and SV calls on 1KGP samples.
- **GATK Structural Variation (SV) Data (`gatk-sv-data`):** Data needed to run Broad's SV discovery pipeline for Illumina short-read WGS data.
- **Epigenomes of the HPRC Release 2 (`hprc-epigenome`):** High-quality phased genome assemblies from over 200 individuals with comprehensive functional epigenomic data.
- **Allen Ivy Glioblastoma Atlas:** Gene expression and tissue imaging for glioblastoma human brain tumors.
- **Beat Acute Myeloid Leukemia (AML) 1.0:** Collaborative research sequencing data for AML.
- **Clinical Trial Sequencing Project - Diffuse Large B-Cell Lymphoma:** Recurrent genetic alterations (including deletions and amplifications/CNVs) in DLBCL.
- **Garvan Institute Long Read Sequencing Benchmark Data (`gtgseq`):** A benchmark resource of GIAB samples (HG001, HG002) sequenced heavily on Oxford Nanopore PromethION.
- **Somatic Mosaicism across Human Tissues (`smaht`):** NIH consortium characterizing somatic variation (mosaicism) in normal tissues, heavily utilizing high-accuracy PacBio HiFi data.
- **Nanopore Reference Human Genome (`nanopore`):** The original reference standard human genome (GM12878) using MinION nanopore sequencing.
