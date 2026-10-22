from pathlib import Path

import pytest

from seq_qc_probe.metrics import analyze_fastq_file

FIXTURES_DIR = Path(__file__).parent / "fixtures"


def test_ground_truth_exact_mathematical_validation():
    """Validates the analytical precision of all metrics on benchmark_truth.fastq.

    Analytical Ground Truth:
    - Total Reads: 4 (each exactly 20 bp)
    - Total Bases: 80 bp
    - Read 1 Mean Phred: 'I' -> ord('I') - 33 = 40.0
    - Read 2 Mean Phred: '?' -> ord('?') - 33 = 30.0
    - Read 3 Mean Phred: '5' -> ord('5') - 33 = 20.0
    - Read 4 Mean Phred: '+' -> ord('+') - 33 = 10.0
    - Global Mean Quality (Mean Q): (40 + 30 + 20 + 10) / 4 = 25.00
    - Cycle / Positional Means: (40 + 30 + 20 + 10) / 4 = 25.00 across all 20 positions
    - GC Content:
        * Read 1: 0 GC
        * Read 2: 7 GC (5 G + 2 C in AGATCGGAAGAGC)
        * Read 3: 20 GC (10 C + 10 G)
        * Read 4: 10 GC (5 C + 5 G)
        * Total: 37 GC / 80 Bases = 46.25%
    - Adapter Contamination:
        * Read 2 contains the exact sequence AGATCGGAAGAGC.
        * 1 / 4 = 25.00%
    """
    benchmark_file = FIXTURES_DIR / "benchmark_truth.fastq"
    summary, errors = analyze_fastq_file(str(benchmark_file))

    # 1. Structural validation errors must be zero
    assert len(errors) == 0

    # 2. Total reads and bases evaluation
    assert summary.total_reads_evaluated == 4
    assert summary.total_bases_evaluated == 80

    # 3. Phred mean validation
    assert summary.mean_quality_score == pytest.approx(25.0, abs=1e-5)

    # 4. Positional cycle averages (all 20 positions must be exactly 25.0)
    assert len(summary.positions_average_quality) == 20
    assert all(
        pos_q == pytest.approx(25.0, abs=1e-5) for pos_q in summary.positions_average_quality
    )

    # 5. Biological metrics (GC%, adapter contamination rate, Q20/Q30/N)
    assert summary.gc_content_percent == pytest.approx(46.25, abs=1e-4)
    assert summary.adapter_contaminated_reads == 1
    assert summary.adapter_contamination_rate == pytest.approx(25.0, abs=1e-4)
    # Reads 1, 2, 3 have Q >= 20 -> 60/80 = 75.0%
    assert summary.q20_percent == pytest.approx(75.0, abs=1e-4)
    # Reads 1, 2 have Q >= 30 -> 40/80 = 50.0%
    assert summary.q30_percent == pytest.approx(50.0, abs=1e-4)
    assert summary.n_content_percent == pytest.approx(0.0, abs=1e-4)

    assert summary.read_quality_distribution == {10: 1, 20: 1, 30: 1, 40: 1}
