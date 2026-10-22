import io
import json
import sys
from pathlib import Path

import pytest

from seq_qc_probe.cli import main

FIXTURES_DIR = Path(__file__).parent / "fixtures"


def test_cli_check_valid_file(capsys: pytest.CaptureFixture):
    """Test CLI check subcommand on a valid FASTQ file returns exit code 0 and clean output."""
    valid_file = FIXTURES_DIR / "valid.fastq"
    with pytest.raises(SystemExit) as exc:
        main(["check", "-i", str(valid_file), "--min-q", "30"])

    assert exc.value.code == 0
    captured = capsys.readouterr()
    assert "Validation PASSED" in captured.out
    assert "Total Reads Evaluated: 2" in captured.out
    assert "Mean Quality Score: 39.70" in captured.out
    assert "GC Content Percent: 48.00%" in captured.out


def test_cli_missing_required_input(capsys: pytest.CaptureFixture):
    """Test CLI fails with exit code 2 when required -i/--input flag is missing."""
    with pytest.raises(SystemExit) as exc:
        main(["check"])

    assert exc.value.code == 2
    captured = capsys.readouterr()
    assert "required" in captured.err or "error" in captured.err


def test_cli_file_not_found(capsys: pytest.CaptureFixture):
    """Test CLI exits with code 2 when specified input file does not exist."""
    with pytest.raises(SystemExit) as exc:
        main(["check", "-i", "/path/to/nonexistent_file.fastq"])

    assert exc.value.code == 2
    captured = capsys.readouterr()
    assert "Error: File not found" in captured.err


def test_cli_corrupted_file_exits_with_code_1(capsys: pytest.CaptureFixture):
    """Test CLI exits with code 1 when structural validation errors are detected."""
    corrupted_file = FIXTURES_DIR / "err_length_mismatch.fastq"
    with pytest.raises(SystemExit) as exc:
        main(["check", "-i", str(corrupted_file)])

    assert exc.value.code == 1
    captured = capsys.readouterr()
    assert "MISMATCHED_LENGTHS" in captured.err
    assert "Validation FAILED" in captured.err


def test_cli_min_quality_threshold_failure(capsys: pytest.CaptureFixture):
    """Test CLI exits with code 1 when file passes structural checks but fails --min-q."""
    benchmark_file = FIXTURES_DIR / "benchmark_truth.fastq"
    with pytest.raises(SystemExit) as exc:
        main(["check", "-i", str(benchmark_file), "--min-q", "30"])

    assert exc.value.code == 1
    captured = capsys.readouterr()
    assert "is below the minimum threshold" in captured.err


def test_cli_max_errors_aborts(capsys: pytest.CaptureFixture):
    """Test CLI aborts processing when maximum error threshold is reached."""
    corrupted_file = FIXTURES_DIR / "err_bad_header.fastq"
    with pytest.raises(SystemExit) as exc:
        main(["check", "-i", str(corrupted_file), "--max-errors", "1"])

    assert exc.value.code == 1
    captured = capsys.readouterr()
    assert "Maximum error threshold (1) reached" in captured.err


def test_cli_quiet_mode_suppresses_output(capsys: pytest.CaptureFixture):
    """Test CLI --quiet flag suppresses standard output and error."""
    valid_file = FIXTURES_DIR / "valid.fastq"
    with pytest.raises(SystemExit) as exc:
        main(["check", "-i", str(valid_file), "-q"])

    assert exc.value.code == 0
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_cli_json_output(capsys: pytest.CaptureFixture):
    """Test CLI --json flag outputs valid JSON summary."""
    valid_file = FIXTURES_DIR / "valid.fastq"
    with pytest.raises(SystemExit) as exc:
        main(["check", "-i", str(valid_file), "--json"])

    assert exc.value.code == 0
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert payload["status"] == "PASSED"
    assert payload["total_records"] == 2
    assert payload["total_errors"] == 0
    assert payload["total_reads_evaluated"] == 2
    assert payload["total_bases_evaluated"] == 50
    assert payload["mean_quality_score"] == 39.7
    assert payload["gc_content_percent"] == 48.0


def test_cli_stdin_pipe(capsys: pytest.CaptureFixture):
    """Test CLI reading from standard input ('-i -')."""
    raw_content = "@READ_STDIN\nACGT\n+\nIIII\n"
    original_stdin = sys.stdin
    try:
        sys.stdin = io.StringIO(raw_content)
        with pytest.raises(SystemExit) as exc:
            main(["check", "-i", "-", "--json"])
    finally:
        sys.stdin = original_stdin

    assert exc.value.code == 0
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert payload["status"] == "PASSED"
    assert payload["total_reads_evaluated"] == 1
    assert payload["total_bases_evaluated"] == 4


def test_cli_output_json_file(tmp_path: Path):
    """Test CLI -o/--output-json flag writes structured JSON report to a file."""
    valid_file = FIXTURES_DIR / "valid.fastq"
    output_file = tmp_path / "reports" / "summary.json"

    with pytest.raises(SystemExit) as exc:
        main(["check", "-i", str(valid_file), "-o", str(output_file)])

    assert exc.value.code == 0
    assert output_file.is_file()

    with open(output_file, encoding="utf-8") as f:
        data = json.load(f)

    # Check hierarchical schema
    assert data["schema_version"] == "1.0.0"
    assert data["meta"]["tool"] == "seq-qc-probe"
    assert data["meta"]["version"] == "0.1.0"
    assert data["quality_gate"]["status"] == "PASSED"
    assert data["quality_gate"]["passed"] is True
    assert data["summary"]["total_reads"] == 2
    assert data["summary"]["total_bases"] == 50
    assert data["summary"]["q20_bases"] == 50
    assert data["summary"]["q20_rate"] == 1.0
    assert data["summary"]["gc_bases"] == 24
    assert data["contamination"]["adapter_name"] == "Illumina Universal Adapter"
    assert len(data["distributions"]["per_cycle_mean_quality"]) == 25


def test_cli_output_json_with_quiet(tmp_path: Path, capsys: pytest.CaptureFixture):
    """Test CLI writes output JSON file even when --quiet flag suppresses terminal output."""
    valid_file = FIXTURES_DIR / "valid.fastq"
    output_file = tmp_path / "quiet_report.json"

    with pytest.raises(SystemExit) as exc:
        main(["check", "-i", str(valid_file), "-o", str(output_file), "-q"])

    assert exc.value.code == 0
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""
    assert output_file.is_file()
