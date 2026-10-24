import io
import json
from pathlib import Path

import pytest

from seq_qc_probe.cli import main

FIXTURES_DIR = Path(__file__).parent / "fixtures"


def test_cli_sam_valid_file(capsys: pytest.CaptureFixture):
    """Test CLI sam-check on a compliant SAM file returns exit code 0."""
    valid_file = FIXTURES_DIR / "valid_alignment.sam"
    with pytest.raises(SystemExit) as exc:
        main(["sam-check", "-i", str(valid_file)])

    assert exc.value.code == 0
    captured = capsys.readouterr()
    assert "Validation PASSED" in captured.out
    assert "Total Records Evaluated: 16" in captured.out
    assert "MAPQ & Multimapping Metrics:" in captured.out
    assert "Mean MAPQ: 41.06" in captured.out


def test_cli_sam_missing_required_input(capsys: pytest.CaptureFixture):
    """Test CLI sam-check fails with exit code 2 when required -i/--input flag is missing."""
    with pytest.raises(SystemExit) as exc:
        main(["sam-check"])

    assert exc.value.code == 2
    captured = capsys.readouterr()
    assert "required" in captured.err or "error" in captured.err


def test_cli_sam_file_not_found(capsys: pytest.CaptureFixture):
    """Test CLI sam-check exits with code 2 when input file does not exist."""
    with pytest.raises(SystemExit) as exc:
        main(["sam-check", "-i", "/path/to/nonexistent_file.sam"])

    assert exc.value.code == 2
    captured = capsys.readouterr()
    assert "Error: File not found" in captured.err


def test_cli_sam_file_not_found_quiet(capsys: pytest.CaptureFixture):
    """Test CLI sam-check exits with code 2 and suppresses error output when quiet."""
    with pytest.raises(SystemExit) as exc:
        main(["sam-check", "-i", "/path/to/nonexistent_file.sam", "-q"])

    assert exc.value.code == 2
    captured = capsys.readouterr()
    assert captured.err == ""
    assert captured.out == ""


def test_cli_sam_cigar_mismatch_exits_with_code_1(capsys: pytest.CaptureFixture):
    """Test CLI sam-check exits with code 1 when CIGAR mismatch errors are detected."""
    corrupted_file = FIXTURES_DIR / "err_cigar_mismatch.sam"
    with pytest.raises(SystemExit) as exc:
        main(["sam-check", "-i", str(corrupted_file)])

    assert exc.value.code == 1
    captured = capsys.readouterr()
    assert "Validation FAILED" in captured.err
    assert "CIGAR_LENGTH_MISMATCH" in captured.err


def test_cli_sam_quiet_mode_suppresses_output(capsys: pytest.CaptureFixture):
    """Test CLI sam-check quiet mode suppresses all stdout and stderr."""
    corrupted_file = FIXTURES_DIR / "err_cigar_mismatch.sam"
    with pytest.raises(SystemExit) as exc:
        main(["sam-check", "-i", str(corrupted_file), "-q"])

    assert exc.value.code == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_cli_sam_json_output(capsys: pytest.CaptureFixture):
    """Test CLI sam-check outputs JSON format to stdout."""
    valid_file = FIXTURES_DIR / "valid_alignment.sam"
    with pytest.raises(SystemExit) as exc:
        main(["sam-check", "-i", str(valid_file), "--json"])

    assert exc.value.code == 0
    captured = capsys.readouterr()
    data = json.loads(captured.out)
    assert data["meta"]["command"] == "sam-check"
    assert data["validation_gate"]["status"] == "PASSED"
    assert data["validation_gate"]["passed"] is True
    assert data["summary"]["total_records_evaluated"] == 16
    assert data["mapq_metrics"]["total_reads"] == 16
    assert data["mapq_metrics"]["mean_mapq"] == 41.06


def test_cli_sam_output_json_file(tmp_path: Path):
    """Test CLI sam-check writes JSON report to a file path."""
    valid_file = FIXTURES_DIR / "valid_alignment.sam"
    out_json = tmp_path / "sam_out.json"

    with pytest.raises(SystemExit) as exc:
        main(["sam-check", "-i", str(valid_file), "-o", str(out_json)])

    assert exc.value.code == 0
    assert out_json.exists()
    data = json.loads(out_json.read_text())
    assert data["validation_gate"]["passed"] is True
    assert data["mapq_metrics"]["total_reads"] == 16


def test_cli_sam_output_json_with_quiet(tmp_path: Path, capsys: pytest.CaptureFixture):
    """Test CLI sam-check writes JSON file even when quiet is enabled."""
    valid_file = FIXTURES_DIR / "valid_alignment.sam"
    out_json = tmp_path / "sam_quiet.json"

    with pytest.raises(SystemExit) as exc:
        main(["sam-check", "-i", str(valid_file), "-o", str(out_json), "-q"])

    assert exc.value.code == 0
    assert out_json.exists()
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_cli_sam_stdin_pipe(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture):
    """Test CLI sam-check accepts input piped from stdin ('-')."""
    valid_content = (
        "@HD\tVN:1.6\tSO:coordinate\n"
        "@SQ\tSN:chr1\tLN:248956422\n"
        "read_01\t0\tchr1\t1000\t60\t10M\t*\t0\t0\tACGTACGTAC\tFFFFFFFFFF\n"
    )
    monkeypatch.setattr("sys.stdin", io.StringIO(valid_content))

    with pytest.raises(SystemExit) as exc:
        main(["sam-check", "-i", "-"])

    assert exc.value.code == 0
    captured = capsys.readouterr()
    assert "Validation PASSED" in captured.out
    assert "Total Records Evaluated: 1" in captured.out


def test_cli_sam_gzip_file(capsys: pytest.CaptureFixture):
    """Test CLI sam-check on a gzip-compressed valid SAM file."""
    valid_gz = FIXTURES_DIR / "valid_alignment.sam.gz"
    with pytest.raises(SystemExit) as exc:
        main(["sam-check", "-i", str(valid_gz)])

    assert exc.value.code == 0
    captured = capsys.readouterr()
    assert "Validation PASSED" in captured.out
    assert "Total Records Evaluated: 16" in captured.out


def test_cli_sam_max_errors_flag(capsys: pytest.CaptureFixture):
    """Test CLI sam-check respects --max-errors threshold."""
    corrupted_file = FIXTURES_DIR / "err_cigar_mismatch.sam"
    with pytest.raises(SystemExit) as exc:
        main(["sam-check", "-i", str(corrupted_file), "--max-errors", "1", "--json"])

    assert exc.value.code == 1
    captured = capsys.readouterr()
    data = json.loads(captured.out)
    assert data["validation_gate"]["total_errors"] == 1
    assert data["meta"]["thresholds"]["max_errors"] == 1


def test_handle_sam_check_direct_call():
    """Test handle_sam_check function directly with parsed arguments."""
    from seq_qc_probe.cli import build_parser, handle_sam_check

    parser = build_parser()
    valid_file = FIXTURES_DIR / "valid_alignment.sam"
    args = parser.parse_args(["sam-check", "-i", str(valid_file), "-q"])
    exit_code = handle_sam_check(args)
    assert exit_code == 0
