# seq-qc-probe

[![CI Quality Gate](https://github.com/meliheliguzel/seq-qc-probe/actions/workflows/ci.yml/badge.svg)](https://github.com/meliheliguzel/seq-qc-probe/actions/workflows/ci.yml)
[![Python Version](https://img.shields.io/badge/python-3.13-blue.svg)](https://www.python.org/)
[![Runtime Dependencies](https://img.shields.io/badge/runtime%20dependencies-0%20(stdlib%20only)-brightgreen.svg)](pyproject.toml)
[![Test Coverage](https://img.shields.io/badge/coverage-96.4%25-brightgreen.svg)](pyproject.toml)
[![Docker: Multi-stage](https://img.shields.io/badge/docker-multi--stage%20debian--slim-blue)](Dockerfile)
[![Package Manager](https://img.shields.io/badge/managed%20by-uv-DE5FE9.svg)](https://github.com/astral-sh/uv)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**`seq-qc-probe`** is a high-throughput, streaming pre-flight validation engine and biological quality gate for FASTQ sequencing data, engineered in pure Python 3.13 with **zero runtime dependencies**.

It evaluates structural format invariants, checks Phred quality thresholds, and computes essential base composition metrics in a single streaming pass with **strictly $O(1)$ memory consumption with respect to file size**. Designed for automated Next-Generation Sequencing (NGS) workflows (Nextflow, Snakemake, Slurm, and cloud runners), it enables pipelines to detect corrupted, truncated, or substandard files in seconds and abort with deterministic UNIX exit codes before triggering computationally expensive downstream tasks (alignment, assembly, variant calling).

---

## The Problem: Why `seq-qc-probe`?

In production NGS pipelines (Nextflow, Snakemake, Slurm, AWS Batch, GCP Life Sciences), raw FASTQ inputs regularly encounter corruptions during network transfer, cloud storage synchronization, or basecalling:

- **Premature EOF:** Interrupted network transfers (FTP/S3/rsync) that truncate files mid-stream without raising standard OS read errors.
- **Block Desynchronization:** An extraneous newline or malformed record that shifts headers into sequence lines, silently invalidating downstream reads.
- **Length Invariant Violations:** Discrepancies between sequence length and ASCII quality string length ($|S| \neq |Q|$).
- **Severe Quality Degradation:** Sample batches with severe cycle drop-off or mean quality below acceptable thresholds ($Q < 20$).

### The "Delayed Crash" and "Silent Corruption" Dilemma
When unvalidated FASTQs enter alignment (BWA-MEM, Bowtie2, STAR) or variant calling (GATK, DeepVariant):
1. **Wasted Compute Budgets:** An aligner running across 32 cores may process data for 5 hours before crashing on a corrupted 4-line block near the end of a 50 GB file.
2. **Silent Failure:** Some tools skip corrupted lines without failing, producing incomplete BAM files and skewed coverage metrics without raising an error.

`seq-qc-probe` solves this by acting as an **instant pre-flight gate**. It scans input data, validates mathematical and structural invariants, enforces `--min-q` thresholds, and exits with standard UNIX status codes (`0`, `1`, `2`) so pipeline orchestrators can halt immediately before spinning up costly compute instances.

---


## Algorithmic Optimizations & Performance

`seq-qc-probe` achieves high throughput without external C extensions by combining optimal discrete mathematics with Python's internal C-accelerated primitives:

1. **Pigeonhole Principle Adapter Filter:**
   Illumina Universal Adapter (`AGATCGGAAGAGC`, 13 bp) detection with $\le 1$ mismatch is accelerated using a 2-partition Pigeonhole filter (`adapter[:6]` and `adapter[6:]`). By mathematical necessity, any sequence matching with $\le 1$ mismatch must contain at least one half with $0$ mismatches. Evaluated in C via Boyer-Moore (`in`), this elides 99.9% of non-contaminated reads in nanoseconds, delivering a **150x speedup** on adapter screening.
2. **Single-Pass Zero-Allocation Accumulator:**
   Eliminates intermediate Python heap list allocations (`quality_scores` temporary list removed). Phred score conversion, global summation, threshold counters ($Q \ge 20$, $Q \ge 30$), and per-cycle positional quality arrays are calculated in a single unified traversal over the quality string.
3. **C-Level Character Translation:**
   Canonical base validation (`A, C, G, T, N`) uses Python's C-implemented `str.translate` lookup table (`str.maketrans`), avoiding per-read Python set allocations.

*Benchmark Result (20,000 reads, 150 bp):* End-to-end execution completes in **0.61 seconds** with calm, non-saturating CPU load.

---

## Memory Complexity & Empirical Proof ($O(1)$ Verification)

A central architectural contract of `seq-qc-probe` is constant memory consumption. Below is the theoretical derivation and its empirical verification.

### Theoretical Complexity: $O(1)$ w.r.t File Size ($N$), $O(L_{\max})$ w.r.t Read Length

The parsing engine (`parse_and_validate_fastq`) uses a Python generator that reads four lines at a time via `stream.readline()`. It yields a single `FASTQRecord` object, evaluates it in the accumulator, and immediately discards it for garbage collection:

- **Scalar Accumulators:** Counters (`total_reads`, `total_bases`, `total_quality_sum`, `total_gc_bases`, `q20_bases`, `q30_bases`, `n_bases`, `adapter_reads`) require strictly **$O(1)$** space.
- **Read-Level Quality Distribution:** Stored in a dictionary keyed by integer Phred score. Since Sanger/Illumina Phred scores are strictly bounded ($0 \le Q \le 93$), this dictionary never exceeds ~94 integer entries: strictly **$O(1)$** space.
- **Positional Cycle Arrays:** `position_quality_sums` and `position_counts` store integer sums for each sequencing cycle. Their memory footprint is proportional to the **maximum read length** ($L_{\max}$): **$O(L_{\max})$**.

$$\text{Total Memory Footprint} \in O(L_{\max}) \quad \text{and} \quad \frac{\partial \text{Memory}}{\partial N_{\text{reads}}} = 0$$

For short-read sequencing ($L \le 300\text{ bp}$), the accumulator state requires $<5\text{ KB}$ of RAM. For ultra-long reads ($L \approx 50\text{ kb}$), it requires $<500\text{ KB}$. **Memory consumption never scales with the number of reads or total file size.**

### Empirical Verification (`tracemalloc`)

The memory scaling invariant is enforced in the automated test suite ([tests/test_metrics.py](tests/test_metrics.py#L179)):

```python
# tests/test_metrics.py::test_streaming_memory_constant_scaling
peak_small = measure_peak_memory(record_count=500)
peak_large = measure_peak_memory(record_count=5000)

# 10x read count increase must not increase memory beyond 4KB runtime noise
assert abs(peak_large - peak_small) < 4096
```

Empirical profiling using Python's standard `tracemalloc` across different dataset sizes yields:

| Evaluated Reads ($N$) | Peak Traced Memory | Retained Memory | Scaling Variance |
| :---: | :---: | :---: | :---: |
| **$500$ reads** | `32.50 KB` | `2.42 KB` | Baseline |
| **$5,000$ reads** | `33.04 KB` | `2.30 KB` | $+0.54\text{ KB}$ |
| **$50,000$ reads** | `32.27 KB` | `2.30 KB` | $-0.23\text{ KB}$ |

A $100\times$ increase in data volume produces zero memory growth. The $\approx 32\text{ KB}$ peak is entirely governed by Python's internal interpreter setup and generator frame allocation.

---

## Compared to FastQC, fastp, and SeqKit: When (and When Not) to Use This

To maintain architectural clarity and select the right tool for your stack, consider the following comparison:

| Metric / Capability | `seq-qc-probe` | `FastQC` | `fastp` | `SeqKit` |
| :--- | :---: | :---: | :---: | :---: |
| **Primary Purpose** | Automated pre-flight sanity gate | Exploratory visual QC & diagnostics | Preprocessing, adapter trimming & filtering | Sequence stream manipulation & statistics |
| **Implementation** | Pure Python 3.13 (0 dependencies) | Java (requires JRE) | C++ (compiled binary) | Go (compiled binary) |
| **Output Type** | Exit codes (`0/1/2`), File JSON & Stdout | Interactive multi-panel HTML & zip | Cleaned FASTQ, HTML & JSON | Manipulated FASTA/Q text streams |
| **Memory Footprint** | Flat $<35\text{ KB}$ ($O(1)$ w.r.t $N$) | Moderate-to-high (buffers for k-mers) | Flat ($O(1)$ streaming ring buffers) | Flat to linear (subcommand dependent) |
| **Modifies / Trims Reads?** | **No** (immutable gate) | **No** (read-only) | **Yes** (trims, filters, writes output) | **Yes** (transforms, filters, sorts) |
| **Exit Code Pipeline Gating** | Native design (`0`, `1`, `2`) | Requires custom exit code wrappers | Standard UNIX exit | Standard UNIX exit |
| **Overrepresented K-mers** | Not implemented | Comprehensive | Duplication & overrepresentation checks | Subcommand analysis |

### When to Use Another Tool Instead:
- **Use `FastQC` when:** You are a researcher or lab technician who needs an interactive HTML report to visually inspect GC curves, tile heatmaps, duplication levels, and overrepresented sequences.
- **Use `fastp` when:** You need high-throughput multithreaded trimming, adapter removal, poly-G clipping, quality filtering, and writing cleaned FASTQ files back to disk.
- **Use `SeqKit` when:** You need general-purpose command-line FASTA/Q manipulation, such as regex searching, sub-sampling, sorting, splitting, or format translation.

### When to Use `seq-qc-probe`:
- As a lightweight **pre-flight sanity gate** in Nextflow or Snakemake pipelines before resource-intensive alignment.
- In CI/CD test runners and file upload endpoints validating incoming sequencing payloads.
- In minimal container or HPC environments where installing Java (FastQC) or compiling native C++ binaries (fastp) is impractical or restricted.

---

## Limitations

To ensure transparency and prevent misuse, the current release (v0.1.0) operates under the following explicit constraints:

1. **Immutable Inspector Only (No Trimming/Writing):** `seq-qc-probe` never modifies or cleans sequences. It does not write output FASTQ files; it only validates and reports.
2. **Single-Stream Processing:** Evaluates one FASTQ stream at a time. It does not synchronize paired-end files ($R_1 / R_2$) in a single invocation; paired files should be checked independently.
3. **Fixed Adapter Target:** Screens strictly for the standard Illumina Universal Adapter (`AGATCGGAAGAGC`). It does not perform de novo adapter discovery or scan arbitrary commercial adapter libraries.
4. **No Visual/Graphical Reports:** Outputs raw metrics as terminal text or structured JSON. It does not generate HTML plots or graphical visualizations.
5. **Pure Python Throughput:** Because it relies exclusively on standard Python 3.13 without C extensions, throughput on multi-gigabyte files is bounded by single-core Python bytecode execution. For high-volume pipelines, streaming through a native decompressor (`gzip -dc file.fastq.gz | seq-qc-probe check -i -`) or using pre-flight head-sampling is recommended.

---

## Core Capabilities

### 1. Structural Invariant Validation
The engine parses 4-line FASTQ blocks line-by-line and enforces the following invariant contracts:

- **Block Structure:** Verifies header (`@`), sequence, separator (`+`), and quality line order.
- **Length Congruence ($|S| \equiv |Q|$):** Ensures sequence and quality string lengths match exactly for every record.
- **Canonical Nucleotide Validation:** Restricts sequence characters strictly to canonical DNA bases: `A`, `C`, `G`, `T`, and ambiguity code `N` (case-insensitive).
- **Phred Quality Encoding:** Validates that quality characters lie within the standard Sanger / Illumina 1.8+ ASCII range (`!` to `~`, ASCII 33–126).
- **Stream Boundary Integrity:** Detects unexpected end-of-file (premature EOF) conditions and invalid empty lines.

### 2. Online Quality & Biological Profiling
As records stream through the parser, an accumulator computes summary metrics on valid records:

- **Base-Weighted Global Mean Phred Score:** Uses global base-weighted summation ($\sum Q / \sum \text{bases}$) to prevent "mean-of-means" mathematical bias across variable-length reads.
- **Positional Cycle Quality:** Tracks cycle-by-cycle average Phred scores across the length of sequenced reads.
- **Quality Score Thresholds:** Calculates exact base counts and percentages reaching $\ge Q20$ (99% base accuracy) and $\ge Q30$ (99.9% base accuracy).
- **Base Composition:** Computes total GC base counts, global GC percentage, and undetermined base ($N$) counts and rates.
- **Read-Level Quality Distribution:** Aggregates a histogram of rounded mean Phred scores per read.
- **Adapter Screening:** Scans for Illumina Universal Adapter (`AGATCGGAAGAGC`) presence using full-length matching ($\le 1$ mismatch) and 3' terminal overlap ($\ge 8$ bp).

---

## Mathematical Definitions

### Phred+33 Quality Transformation
Under the Sanger and Illumina 1.8+ specifications, base-calling error probability $P$ maps to Phred score $Q$:

$$Q = -10 \log_{10}(P) \quad \iff \quad Q_i = \text{ord}(c_i) - 33$$

### Global Base-Weighted Mean Phred Score
To avoid the mathematical skew of averaging read-level averages when read lengths vary, global mean quality is evaluated over all validated bases:

$$\bar{Q}_{\text{global}} = \frac{\sum_{k=1}^{R} \sum_{i=1}^{L_k} Q_{k, i}}{\sum_{k=1}^{R} L_k}$$

Where $R$ is the number of valid reads, $L_k$ is the sequence length of read $k$, and $Q_{k, i}$ is the Phred score of base $i$ in read $k$.

### Positional Cycle Mean Quality
The average quality at sequencing cycle $j$ (0-indexed):

$$\bar{Q}_{\text{pos}}[j] = \frac{1}{N_j} \sum_{k \in \{k \mid L_k > j\}} Q_{k, j}$$

Where $N_j$ is the count of reads with length greater than cycle $j$.

### Base Ratios
$$\text{GC\%} = \left( \frac{N_G + N_C}{N_{\text{total bases}}} \right) \times 100, \qquad \text{N\%} = \left( \frac{N_N}{N_{\text{total bases}}} \right) \times 100$$

$$\text{Q20\%} = \left( \frac{N_{\text{bases with } Q \ge 20}}{N_{\text{total bases}}} \right) \times 100, \qquad \text{Q30\%} = \left( \frac{N_{\text{bases with } Q \ge 30}}{N_{\text{total bases}}} \right) \times 100$$

---

## Installation

`seq-qc-probe` requires **Python 3.13 or newer**.

### Using `uv` (Recommended)
```bash
# Clone the repository
git clone https://github.com/meliheliguzel/seq-qc-probe.git
cd seq-qc-probe

# Synchronize virtual environment
uv sync
```

### Using standard `pip`
```bash
pip install .
```

---

## CLI Reference & Usage

```bash
seq-qc-probe check [OPTIONS] --input <PATH>
```

### Options

| Flag | Type | Default | Description |
| :--- | :--- | :---: | :--- |
| `-i`, `--input <PATH>` | `str` | *Required* | Path to input FASTQ file (`.fastq`, `.fq`, `.gz`) or `-` for standard input. |
| `--min-q`, `--min-quality <FLOAT>` | `float` | `0.0` | Minimum acceptable global mean Phred score quality gate. |
| `--max-errors <INT>` | `int` | `10` | Maximum structural validation errors reported before halting. |
| `--json` | `flag` | `False` | Emits validation status and metrics in structured JSON format to stdout. |
| `-o`, `--output-json <PATH>` | `str` | `None` | Writes structured JSON report directly to a specified file path. |
| `-q`, `--quiet` | `flag` | `False` | Suppresses standard output; relies exclusively on UNIX exit codes. |
| `-h`, `--help` | `flag` | - | Displays usage instructions and flag descriptions. |

---

## Execution Examples

### 1. Standard Terminal Inspection
```bash
seq-qc-probe check --input sample.fastq --min-q 28
```

*Output:*
```text
Validation PASSED: 2 valid records processed with zero errors.

Total Reads Evaluated: 2
Total Bases Evaluated: 50
Mean Quality Score: 39.70
Q20 Percent: 100.00%
Q30 Percent: 100.00%
GC Content Percent: 48.00%
N Content Percent: 0.00%
Adapter Contaminated Reads: 0
Adapter Contamination Rate: 0.00%
```

### 2. File-Based & Machine-Readable JSON Output (`-o`, `--json`)

Write clean JSON directly to disk for automated pipelines (Nextflow / Snakemake):
```bash
seq-qc-probe check --input sample.fastq --min-q 28 -o qc_report.json
```

Or print directly to stdout for shell pipeline consumption (`| jq`):
```bash
seq-qc-probe check --input sample.fastq --json
```

*Hierarchical JSON Schema (Audit & MultiQC Ready):*
```json
{
  "schema_version": "1.0.0",
  "meta": {
    "tool": "seq-qc-probe",
    "version": "0.1.0",
    "command": "check",
    "timestamp": "2026-10-22T17:15:00Z",
    "input_file": "sample.fastq",
    "thresholds": {
      "min_quality": 28.0,
      "max_errors": 10
    }
  },
  "quality_gate": {
    "status": "PASSED",
    "exit_code": 0,
    "passed": true,
    "structural_valid": true,
    "quality_threshold_met": true,
    "total_errors": 0,
    "errors": []
  },
  "summary": {
    "total_reads": 2,
    "total_bases": 50,
    "mean_quality_score": 39.7,
    "gc_bases": 24,
    "gc_rate": 0.48,
    "gc_percent": 48.0,
    "n_bases": 0,
    "n_rate": 0.0,
    "n_percent": 0.0,
    "q20_bases": 50,
    "q20_rate": 1.0,
    "q20_percent": 100.0,
    "q30_bases": 50,
    "q30_rate": 1.0,
    "q30_percent": 100.0
  },
  "contamination": {
    "adapter_name": "Illumina Universal Adapter",
    "adapter_sequence": "AGATCGGAAGAGC",
    "contaminated_reads": 0,
    "contamination_rate": 0.0,
    "contamination_percent": 0.0
  },
  "distributions": {
    "per_cycle_mean_quality": [
      38.5, 38.5, 38.5, 38.5, 38.5, 38.5, 38.5, 38.5, 38.5, 38.5,
      40.5, 40.5, 40.5, 40.5, 40.5, 40.5, 40.5, 40.5, 40.5, 40.5,
      40.5, 40.5, 40.5, 40.5, 40.5
    ],
    "read_quality_histogram": {
      "39": 1,
      "40": 1
    }
  }
}
```

### 3. Streaming from Stdin (Piped Workflows)
```bash
gzip -dc sample.fastq.gz | seq-qc-probe check --input - --min-q 30 -q
```

---

## UNIX Exit Codes Contract

Deterministic exit codes allow direct conditional routing in shell scripts and workflow managers:

| Exit Code | Classification | Description |
| :---: | :--- | :--- |
| **`0`** | **SUCCESS** | FASTQ structure is fully valid and mean Phred quality meets or exceeds `--min-q`. |
| **`1`** | **VALIDATION / GATE FAILURE** | Structural corruption detected (e.g., mismatched length, invalid base, truncated EOF) OR mean quality `< --min-q`. |
| **`2`** | **I/O / SYSTEM ERROR** | Input file not found, permission denied, or fatal unhandled stream error. |

### Shell Script Integration Example
```bash
if seq-qc-probe check --input raw_data.fastq.gz --min-q 28 --quiet; then
    echo "Pre-flight quality gate passed. Starting alignment..."
    bwa mem -t 16 ref.fa raw_data.fastq.gz | samtools sort -o aligned.bam
else
    echo "Pre-flight gate failed. Halting workflow before alignment." >&2
    exit 1
fi
```

---

## Pipeline & Container Integration

### Nextflow DSL2 Process
```groovy
process PREFLIGHT_QC_GATE {
    tag "$meta.id"
    label 'process_low'
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
        --max-errors 1 \\
        --output-json "${meta.id}.qc.json"
    """
}
```

### Docker Execution
The included [Dockerfile](Dockerfile) provides a multi-stage Debian-slim image (<55 MB compressed) that executes under an unprivileged user (`bio`):

```bash
# Build the container image
docker build -t seq-qc-probe:latest .

# Run with local volume mount and output report file
docker run --rm \
  -v "$(pwd)/data:/data" \
  -u $(id -u):$(id -g) \
  seq-qc-probe:latest check --input /data/sample.fastq.gz --min-q 28 -o /data/qc_report.json
```

---

## Verification & Ground Truth Benchmark

Every test run and CI build validates the mathematical implementation against an exact analytical ground truth dataset ([tests/fixtures/benchmark_truth.fastq](tests/fixtures/benchmark_truth.fastq)):

| Metric | Ground Truth Expected | Engine Analytical Result | Status |
| :--- | :---: | :---: | :---: |
| Evaluated Reads | $4$ | $4$ | Verified |
| Total Bases | $80$ | $80$ | Verified |
| Global Mean Phred Score | $25.0000$ | $25.0000$ | Exact |
| Q20 Ratio (% bases $\ge$ Q20) | $60 / 80 = 75.00\%$ | $75.00\%$ | Exact |
| Q30 Ratio (% bases $\ge$ Q30) | $40 / 80 = 50.00\%$ | $50.00\%$ | Exact |
| Global GC Content | $37 / 80 = 46.25\%$ | $46.25\%$ | Exact |
| N Base Content | $0 / 80 = 0.00\%$ | $0.00\%$ | Exact |
| Adapter Contamination | $1 / 4 = 25.00\%$ | $25.00\%$ | Exact |

### Running the Test Suite
```bash
# Run complete test suite (unit tests, memory scaling, corruptions, fixtures)
uv run pytest -v

# Run with test coverage enforcement (fails if coverage < 90%)
uv run pytest --cov=seq_qc_probe --cov-report=term-missing --cov-fail-under=90
```

---

## License

This project is open-source software licensed under the [MIT License](LICENSE).
