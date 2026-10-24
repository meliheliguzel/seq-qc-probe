# seq-qc-probe

[![CI Quality Gate](https://github.com/meliheliguzel/seq-qc-probe/actions/workflows/ci.yml/badge.svg)](https://github.com/meliheliguzel/seq-qc-probe/actions/workflows/ci.yml)
[![Python Version](https://img.shields.io/badge/python-3.13-blue.svg)](https://www.python.org/)
[![Runtime Dependencies](https://img.shields.io/badge/runtime%20dependencies-0%20(stdlib%20only)-brightgreen.svg)](pyproject.toml)
[![Test Coverage](https://img.shields.io/badge/coverage-97.2%25-brightgreen.svg)](pyproject.toml)
[![Docker: Multi-stage](https://img.shields.io/badge/docker-multi--stage%20debian--slim-blue)](Dockerfile)
[![Package Manager](https://img.shields.io/badge/managed%20by-uv-DE5FE9.svg)](https://github.com/astral-sh/uv)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**`seq-qc-probe`** is a high-throughput, single-pass streaming pre-flight validation engine and quality gatekeeper for Next-Generation Sequencing (NGS) data (**FASTQ**, **VCF 4.2**, and **SAM/BAM**). Implemented in pure Python 3.13 with **zero external runtime dependencies**, it operates with strictly **$O(1)$ constant memory** with respect to input file size.

Engineered for automated bioinformatics pipelines (Nextflow, Snakemake, Slurm, AWS Batch, GCP Life Sciences), `seq-qc-probe` verifies format invariants, mathematical consistency, and biological quality thresholds in seconds. It halts execution via deterministic UNIX exit codes (`0`, `1`, `2`) before pipelines commit expensive compute resources (alignment, assembly, variant calling) to corrupted or substandard files.

---

## 🎯 The Problem: Pre-Flight Gatekeeping vs. Late Pipeline Failure

In automated NGS production workflows, input and intermediate files routinely suffer from data corruption introduced during network transfer, cloud storage synchronization, or interrupted preprocessing:

- **Truncated Files (Premature EOF):** Incomplete transfers via rsync, S3, or FTP that cut files mid-record without triggering OS-level I/O errors.
- **Record & Coordinate Desynchronization:** Extraneous whitespace or missing columns that shift record boundaries, or non-positive coordinates ($POS \le 0$) in VCF files.
- **CIGAR & Sequence Length Mismatches:** Alignment records where CIGAR operation sums deviate from the query sequence length.
- **Severe Base Quality Degradation:** Batches with severe cycle decay or average Phred scores falling below project quality gates ($Q < 20$ or $Q < 30$).

### Pipeline Economics: Late Failure vs. Early Detection

| Scenario | Without `seq-qc-probe` | With `seq-qc-probe` |
| :--- | :--- | :--- |
| **Corrupted 50 GB FASTQ** | BWA-MEM processes data across 32 cores for 4 hours, then crashes near EOF. **Wasted: ~128 core-hours.** | Fails pre-flight gate in seconds. **Workflow halts immediately (Exit 1).** |
| **Malformed VCF Header** | Downstream annotation (VEP/bcftools) fails 2 hours into a cluster run due to a missing `#CHROM` column. | Verifies format compliance and mandatory columns prior to job dispatch. |
| **Invalid CIGAR in SAM** | Variant caller encounters illegal CIGAR operations mid-chromosome, producing partial calls or crashing silently. | Evaluates CIGAR-to-read sequence length equivalence during stream or upon write. |

---

## ⚖️ Functional Scope & Boundaries

To ensure architectural clarity, `seq-qc-probe` operates within explicit, well-defined boundaries:

### What `seq-qc-probe` DOES:
- ✅ **FASTQ Structural & Quality Gating (`check`):** Enforces 4-line record invariants, sequence/quality length congruence ($|S| \equiv |Q|$), canonical IUPAC nucleotide alphabet (`A, C, G, T, N`), Phred+33 ASCII boundaries, base-weighted global mean Phred calculation, threshold distributions ($Q \ge 20, Q \ge 30$, GC%, positional cycle quality), adapter screening, and `--min-q` gating.
- ✅ **VCF 4.2 Integrity Validation (`vcf-check`):** Validates the mandatory `##fileformat=VCFv4.2` header, verifies the 8 mandatory column headers (`#CHROM` to `INFO`), enforces positive 1-based coordinates ($POS \ge 1$), validates REF/ALT base sequences, and checks format-sample token cardinality.
- ✅ **SAM / Streaming BAM Validation (`sam-check`):** Validates all 11 mandatory tab-delimited columns, enforces mathematical equivalence between read-consuming CIGAR operations and read sequence length ($\sum (M + I + S + = + X) == \text{len}(SEQ)$), respects non-query operations ($D, N, H, P$), validates MAPQ range ($0 \le \text{MAPQ} \le 255$), and computes $MAPQ=0$ multimapping percentage and mean MAPQ.
- ✅ **Universal Transparent Streaming:** Streams uncompressed text, gzip archives (`.gz`), and standard input (`-i -`) natively without temporary disk decompression.
- ✅ **Automation & Audit-Ready Output:** Deterministic UNIX exit codes (`0` = Passed, `1` = Validation Failure, `2` = I/O Error) and machine-readable JSON reports (`-o`, `--json`).

### What `seq-qc-probe` DOES NOT Do:
- ❌ **It is NOT an aligner:** Does not map reads to references (use BWA-MEM, Bowtie2, Minimap2, or STAR).
- ❌ **It is NOT a variant caller:** Does not identify SNVs, indels, or call genotypes (use GATK, DeepVariant, or FreeBayes).
- ❌ **It is NOT a read modifier or trimmer:** Operates as an immutable, read-only inspector; never writes modified FASTQ sequences to disk (use fastp or Trimmomatic).
- ❌ **It is NOT an exploratory visual dashboard:** Generates machine-readable JSON and concise terminal summaries, not interactive HTML plots (use FastQC or MultiQC for visual reports).
- ❌ **It does NOT decode binary BAM/CRAM directly:** Operates on standard SAM text streams. Binary BAM files are validated via standard samtools streaming: `samtools view -h sample.bam | seq-qc-probe sam-check -i -`.

---

## 📐 Mathematical Invariants & Formulation

### 1. Phred+33 Base-Weighted Global Mean Quality (FASTQ)
Under Sanger and Illumina 1.8+ standards, Phred score $Q = -10 \log_{10}(P) \iff Q_i = \text{ord}(c_i) - 33$. To eliminate the mathematical skew of averaging read-level averages across variable-length reads, the global mean is calculated across all validated bases:

$$\bar{Q}_{\text{global}} = \frac{\sum_{k=1}^{R} \sum_{i=1}^{L_k} Q_{k, i}}{\sum_{k=1}^{R} L_k}$$

### 2. CIGAR Sequence Length Invariant (SAM v1.6)
According to SAM v1.6 specification, operations $M$ (alignment match/mismatch), $I$ (insertion to reference), $S$ (soft clip), $=$ (sequence match), and $X$ (sequence mismatch) consume read sequence bases. Deletions ($D$), introns ($N$), hard clips ($H$), and padding ($P$) consume reference or are clipped, and do **not** consume read bases:

$$\sum_{op \in \{M, I, S, =, X\}} \text{length}(op) == \text{len}(SEQ)$$

### 3. MAPQ Multimapping Metric (SAM)
Reads with $\text{MAPQ} = 0$ map with equal likelihood to multiple genomic loci:

$$\text{Multimapping Rate} = \frac{N_{\text{MAPQ}=0}}{N_{\text{total evaluated reads}}}$$

---

## ⚡ Memory Complexity: Strictly $O(1)$ w.r.t File Size

The engine parses records line-by-line using Python streaming generators. Single-pass accumulators record metrics in-flight and discard records immediately:
- **Flat Memory Footprint:** Traced peak memory remains $< 35\text{ KB}$ whether validating $500$ reads or $50,000,000$ reads.
- **Empirical Verification:** Enforced in the automated test suite using Python's standard `tracemalloc` library ([tests/test_metrics.py](tests/test_metrics.py)).

---

## 💻 Installation

Requires **Python 3.13 or newer**.

```bash
# Using uv (Recommended)
uv sync

# Using standard pip
pip install .
```

---

## 🚀 CLI Reference & Commands

```bash
seq-qc-probe <SUBCOMMAND> -i <FILE_OR_STREAM> [OPTIONS]
```

### Common Flags Across Subcommands

| Flag | Type | Default | Description |
| :--- | :---: | :---: | :--- |
| `-i`, `--input <PATH>` | `str` | *Required* | Path to input file (`.fastq`, `.vcf`, `.sam`, `.gz`) or `-` for stdin. |
| `--max-errors <INT>` | `int` | `10` | Maximum validation errors reported before early termination. |
| `-q`, `--quiet` | `flag` | `False` | Suppresses console output; returns exit code only. |
| `--json` | `flag` | `False` | Emits structured JSON validation report to stdout. |
| `-o`, `--output-json <PATH>` | `str` | `None` | Writes structured JSON validation report directly to file. |
| `-h`, `--help` | `flag` | - | Displays usage instructions and available options. |

---

### 1. FASTQ Inspection & Quality Gate (`check`)

```bash
seq-qc-probe check -i sample.fastq.gz --min-q 28.0 -o qc_report.json
```

*FASTQ-Specific Option:*
- `--min-q`, `--min-quality <FLOAT>` (Default: `0.0`): Enforces a minimum global base-weighted Phred quality threshold.

---

### 2. VCF 4.2 Structural Validation (`vcf-check`)

```bash
seq-qc-probe vcf-check -i variants.vcf.gz -o vcf_report.json
```

Validates `##fileformat=VCFv4.2`, verifies 8 mandatory column headers, checks coordinate validity ($POS \ge 1$), and inspects REF/ALT allele bases.

---

### 3. SAM / BAM Alignment Validation (`sam-check`)

```bash
# Validate SAM file directly
seq-qc-probe sam-check -i alignment.sam.gz -o sam_report.json

# Stream directly from BAM via samtools
samtools view -h sample.bam | seq-qc-probe sam-check -i - --max-errors 5
```

Validates 11 mandatory columns, verifies CIGAR length equivalence against `SEQ`, and computes MAPQ distribution and $MAPQ=0$ multimapping rates.

---

## 📊 Standardized JSON Output Contract

Every subcommand produces structured JSON adhering to a unified schema:

```json
{
  "schema_version": "1.0.0",
  "meta": {
    "tool": "seq-qc-probe",
    "version": "0.2.0",
    "command": "sam-check",
    "timestamp": "2026-10-24T13:43:19Z",
    "input_file": "alignment.sam.gz",
    "thresholds": {
      "max_errors": 10
    }
  },
  "validation_gate": {
    "status": "PASSED",
    "passed": true,
    "exit_code": 0,
    "total_errors": 0,
    "error_distribution": {},
    "errors": []
  },
  "mapq_metrics": {
    "total_reads": 100000,
    "mapq_zero_count": 3120,
    "mapq_zero_rate": 0.0312,
    "multimapping_rate": 0.0312,
    "mean_mapq": 48.75,
    "mapq_distribution": {
      "60": 85000,
      "0": 3120
    }
  },
  "summary": {
    "total_records_evaluated": 100000,
    "total_errors": 0,
    "is_valid": true
  }
}
```

---

## 🚦 UNIX Exit Codes Contract

Deterministic exit codes allow direct conditional routing in shell scripts and orchestrators:

| Exit Code | Classification | Description |
| :---: | :--- | :--- |
| **`0`** | **SUCCESS** | File structure is valid and quality thresholds are fully satisfied. |
| **`1`** | **VALIDATION FAILURE** | Structural corruption detected (malformed header, length mismatch, CIGAR error, invalid base) OR quality score $< \text{--min-q}$. |
| **`2`** | **I/O ERROR** | Input file not found, permission denied, or invalid CLI parameters. |

### Shell Integration Example

```bash
#!/usr/bin/env bash
set -euo pipefail

# Pre-flight gate before alignment
if seq-qc-probe check -i raw_reads.fastq.gz --min-q 28.0 -q; then
    echo "[PASSED] Pre-flight gate satisfied. Launching alignment..."
    bwa mem -t 16 ref.fa raw_reads.fastq.gz | samtools sort -o aligned.bam
else
    echo "[FAILED] Corrupted or substandard data detected. Aborting." >&2
    exit 1
fi
```

### Nextflow DSL2 Integration Example

```groovy
process PREFLIGHT_QC {
    tag "$meta.id"
    container 'seq-qc-probe:latest'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path(reads), emit: validated_reads
    path "${meta.id}.qc.json"   , emit: json_report

    script:
    """
    seq-qc-probe check \\
        --input "${reads}" \\
        --min-q 28.0 \\
        --output-json "${meta.id}.qc.json"
    """
}
```

---

## 🧪 Testing & Verification

Every build is validated against an exact analytical ground truth dataset ([tests/fixtures/benchmark_truth.fastq](tests/fixtures/benchmark_truth.fastq)):

| Metric | Analytical Truth | Engine Result | Status |
| :--- | :---: | :---: | :---: |
| Evaluated Reads | $4$ | $4$ | Verified |
| Total Bases | $80$ | $80$ | Verified |
| Global Mean Phred Score | $25.0000$ | $25.0000$ | Exact |
| Q20 Ratio (% bases $\ge Q20$) | $75.00\%$ | $75.00\%$ | Exact |
| Q30 Ratio (% bases $\ge Q30$) | $50.00\%$ | $50.00\%$ | Exact |
| Global GC Content | $46.25\%$ | $46.25\%$ | Exact |
| Adapter Contamination | $25.00\%$ | $25.00\%$ | Exact |

```bash
# Execute complete unit and integration test suite (124 tests)
uv run pytest

# Execute test suite enforcing >= 90% coverage
uv run pytest --cov=seq_qc_probe --cov-report=term-missing --cov-fail-under=90
```

---

## 📄 License

This project is licensed under the [MIT License](LICENSE).
