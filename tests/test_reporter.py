from pathlib import Path

from seq_qc_probe.checker import FASTQChecker
from seq_qc_probe.reporter import QCReporter

FIXTURES_DIR = Path(__file__).parent / "fixtures"


def test_fastq_checker_direct_call():
    """Verify FASTQChecker executes programmatically and returns CheckResult."""
    valid_file = FIXTURES_DIR / "valid.fastq"
    checker = FASTQChecker(input_path=valid_file, min_quality=30.0)
    result = checker.run()

    assert result.status == "PASSED"
    assert result.exit_code == 0
    assert result.total_records == 2
    assert result.total_errors == 0
    assert result.summary.total_reads_evaluated == 2


def test_qc_reporter_hierarchical_structure(tmp_path: Path):
    """Verify QCReporter generates correct hierarchical schema and writes file."""
    valid_file = FIXTURES_DIR / "valid.fastq"
    checker = FASTQChecker(input_path=valid_file, min_quality=25.0)
    result = checker.run()

    reporter = QCReporter(result)
    report_dict = reporter.to_dict()

    assert "schema_version" in report_dict
    assert "meta" in report_dict
    assert "quality_gate" in report_dict
    assert "summary" in report_dict
    assert "contamination" in report_dict
    assert "distributions" in report_dict

    # Check metrics
    assert report_dict["summary"]["gc_bases"] == 24
    assert report_dict["summary"]["q20_bases"] == 50
    assert report_dict["summary"]["q30_bases"] == 50

    # Test file writing
    out_file = tmp_path / "custom" / "report.json"
    reporter.write_json(out_file)
    assert out_file.exists()
