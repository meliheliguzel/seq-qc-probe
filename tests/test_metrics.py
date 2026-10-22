from pathlib import Path

import pytest

from seq_qc_probe.metrics import (
    ILLUMINA_UNIVERSAL_ADAPTER,
    SequenceQCAccumulator,
    analyze_fastq_file,
    is_adapter_contaminated,
)
from seq_qc_probe.models import FASTQRecord

FIXTURES_DIR = Path(__file__).parent / "fixtures"


def test_analyze_valid_fastq_metrics():
    """Test that analyzing a valid FASTQ file returns correct metrics."""
    valid_path = FIXTURES_DIR / "valid.fastq"
    summary, errors = analyze_fastq_file(str(valid_path))

    assert len(errors) == 0
    assert summary.total_reads_evaluated == 2
    assert summary.total_bases_evaluated == 50
    assert summary.mean_quality_score == pytest.approx(39.70, rel=1e-2)
    assert len(summary.positions_average_quality) == 25
    assert summary.positions_average_quality[0] == 38.5
    assert summary.positions_average_quality[-1] == 40.5
    assert summary.q20_percent == 100.0
    assert summary.q30_percent == 100.0
    assert summary.n_content_percent == 0.0


def test_variable_read_lengths_base_weighted_average():
    """Test that global mean quality is base-weighted rather than mean-of-means."""
    acc = SequenceQCAccumulator()

    # R1: 2 bp, Phred 40 each (sum = 80)
    acc.accumulate(FASTQRecord("@R1", "AC", "+", "II", 1))
    # R2: 4 bp, Phred 0, 0, 40, 40 (sum = 80)
    acc.accumulate(FASTQRecord("@R2", "ACGT", "+", "!!II", 5))

    summary = acc.summarize()

    assert summary.total_reads_evaluated == 2
    assert summary.total_bases_evaluated == 6
    assert summary.mean_quality_score == pytest.approx(26.6667, rel=1e-3)
    assert summary.positions_average_quality == [20.0, 20.0, 40.0, 40.0]


def test_empty_accumulator_returns_zero_metrics():
    """Test that an empty SequenceQCAccumulator returns zero metrics."""
    acc = SequenceQCAccumulator()
    summary = acc.summarize()

    assert summary.total_reads_evaluated == 0
    assert summary.total_bases_evaluated == 0
    assert summary.mean_quality_score == 0.0
    assert summary.positions_average_quality == []
    assert summary.read_quality_distribution == {}
    assert summary.q20_percent == 0.0
    assert summary.q30_percent == 0.0
    assert summary.n_content_percent == 0.0


def test_read_quality_distribution():
    """Test that read quality distribution is computed correctly."""
    acc = SequenceQCAccumulator()

    acc.accumulate(FASTQRecord("@R1", "A", "+", "!", 1))
    acc.accumulate(FASTQRecord("@R2", "A", "+", "I", 5))

    summary = acc.summarize()

    assert summary.mean_quality_score == 20.0
    assert summary.read_quality_distribution == {0: 1, 40: 1}


def test_corrupted_records_excluded_from_metrics():
    """Test that corrupted records are excluded from the metrics."""
    corrupt_path = FIXTURES_DIR / "err_length_mismatch.fastq"
    summary, errors = analyze_fastq_file(str(corrupt_path))

    assert len(errors) == 1
    assert summary.total_reads_evaluated == 1
    assert summary.mean_quality_score == 40.0
    assert summary.read_quality_distribution == {40: 1}


def test_gc_content_and_adapter_contamination():
    """Test that GC content and adapter contamination metrics are computed correctly."""
    acc = SequenceQCAccumulator()

    acc.accumulate(FASTQRecord("@R1", "ACGT", "+", "IIII", 1))
    acc.accumulate(FASTQRecord("@R2", "GGGGCGCCC" + ILLUMINA_UNIVERSAL_ADAPTER, "+", "I" * 22, 5))

    summary = acc.summarize()

    assert summary.total_reads_evaluated == 2
    assert summary.mean_quality_score == 40.0
    assert summary.read_quality_distribution == {40: 2}
    assert summary.adapter_contaminated_reads == 1
    assert summary.positions_average_quality[0] == 40.0
    assert summary.gc_content_percent == pytest.approx(69.23, rel=1e-2)
    assert summary.adapter_contamination_rate == pytest.approx(50.0, rel=1e-2)


def test_advanced_adapter_detection():
    """Test that 3' partial matches and 1-mismatch adapter occurrences are detected."""
    # 1. Exact 13bp match
    assert is_adapter_contaminated("NNNN" + ILLUMINA_UNIVERSAL_ADAPTER + "NNNN")
    # 2. 1-mismatch match: change 1 base in adapter
    one_mismatch_adapter = "AGATCGGAAGACC"  # G->C at position 12
    assert is_adapter_contaminated("NNNN" + one_mismatch_adapter + "NNNN")
    # 3. 3' partial overlap: read ends with 9bp prefix of adapter
    assert is_adapter_contaminated("ACGTACGT" + ILLUMINA_UNIVERSAL_ADAPTER[:9])
    # 4. Clean read without adapter
    assert not is_adapter_contaminated("ACGTACGTACGTACGT")


def test_q20_q30_and_n_content_metrics():
    """Test that Q20, Q30, and N-content percentages are correctly calculated."""
    acc = SequenceQCAccumulator()

    # R1: sequence with 2 Ns, quality has Q=10 ('+'), Q=20 ('5'), Q=30 ('?'), Q=40 ('I')
    acc.accumulate(FASTQRecord("@R1", "ACNN", "+", "+5?I", 1))

    summary = acc.summarize()

    assert summary.total_bases_evaluated == 4
    assert summary.n_content_percent == 50.0  # 2 Ns out of 4 bases
    # Q scores: 10, 20, 30, 40
    # Q >= 20: 3 (20, 30, 40) -> 75%
    # Q >= 30: 2 (30, 40) -> 50%
    assert summary.q20_percent == 75.0
    assert summary.q30_percent == 50.0


def test_real_variable_crlf_metrics():
    """Test metrics calculation on real-world CRLF variable length FASTQ fixture."""
    file_path = FIXTURES_DIR / "real_variable_crlf.fastq"
    summary, errors = analyze_fastq_file(str(file_path))

    assert len(errors) == 0
    assert summary.total_reads_evaluated == 3
    assert summary.total_bases_evaluated == 100
    # Expected weighted mean: (30*40 + 45*37 + 25*41) / 100 = 38.90
    assert summary.mean_quality_score == pytest.approx(38.90, rel=1e-3)
    assert len(summary.positions_average_quality) == 45
    assert summary.positions_average_quality[0] == pytest.approx((40 + 37 + 41) / 3, rel=1e-3)


def test_real_illumina_pe_metrics():
    """Test metrics on physical gzipped fixture with uncalled N bases and Q2 '#'."""
    file_path = FIXTURES_DIR / "real_illumina_pe.fastq.gz"
    summary, errors = analyze_fastq_file(str(file_path))

    assert len(errors) == 0
    assert summary.total_reads_evaluated == 2
    assert summary.total_bases_evaluated == 120
    # 4 Ns out of 120 bases = 3.333%
    assert summary.n_content_percent == pytest.approx(4 / 120 * 100, rel=1e-3)
    # 116 bases with Q=40 or Q=41 (>= 20 and >= 30) -> 116/120 = 96.667%
    assert summary.q20_percent == pytest.approx(116 / 120 * 100, rel=1e-3)
    assert summary.q30_percent == pytest.approx(116 / 120 * 100, rel=1e-3)


def test_adapter_3prime_partial_fixture():
    """Test detection of realistic 3' partial read-through adapter overlap."""
    file_path = FIXTURES_DIR / "adapter_3prime_partial.fastq"
    summary, errors = analyze_fastq_file(str(file_path))

    assert len(errors) == 0
    assert summary.total_reads_evaluated == 2
    assert summary.adapter_contaminated_reads == 1
    assert summary.adapter_contamination_rate == pytest.approx(50.0, rel=1e-3)


def test_streaming_memory_constant_scaling(tmp_path):
    """Verifies empirical O(1) memory complexity with respect to read count (N).

    Scaling read volume by 10x must not scale peak memory consumption.
    """
    import tracemalloc

    from seq_qc_probe.parser import parse_and_validate_fastq

    def measure_peak_memory(record_count: int) -> int:
        test_file = tmp_path / f"test_{record_count}.fastq"
        with open(test_file, "w", encoding="utf-8") as f:
            for i in range(record_count):
                f.write(f"@READ_{i}\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n")

        tracemalloc.start()
        accumulator = SequenceQCAccumulator()
        for item in parse_and_validate_fastq(test_file):
            if isinstance(item, FASTQRecord):
                accumulator.accumulate(item)

        _, peak_bytes = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        test_file.unlink()
        return peak_bytes

    peak_small = measure_peak_memory(500)
    peak_large = measure_peak_memory(5000)

    # 10x read count increase must result in less than 4KB memory variance (O(1) w.r.t N)
    memory_growth = abs(peak_large - peak_small)
    assert memory_growth < 4096, f"Memory scaled with read count: growth was {memory_growth} bytes"
