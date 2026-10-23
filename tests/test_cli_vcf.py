import io
import json
from pathlib import Path

import pytest

from seq_qc_probe.cli import main

FIXTURES_DIR = Path(__file__).parent / "fixtures"


def test_cli_vcf_valid_file(capsys: pytest.CaptureFixture):
    """Test CLI vcf-check on a compliant VCF file returns exit code 0."""
    valid_file = FIXTURES_DIR / "valid_vcf42.vcf"
    with pytest.raises(SystemExit) as exc:
        main(["vcf-check", "-i", str(valid_file)])

    assert exc.value.code == 0
    captured = capsys.readouterr()
    assert "Validation PASSED" in captured.out
    assert "Total Records Evaluated: 6" in captured.out


def test_cli_vcf_missing_required_input(capsys: pytest.CaptureFixture):
    """Test CLI vcf-check fails with exit code 2 when required -i/--input flag is missing."""
    with pytest.raises(SystemExit) as exc:
        main(["vcf-check"])

    assert exc.value.code == 2
    captured = capsys.readouterr()
    assert "required" in captured.err or "error" in captured.err


def test_cli_vcf_file_not_found(capsys: pytest.CaptureFixture):
    """Test CLI vcf-check exits with code 2 when input file does not exist."""
    with pytest.raises(SystemExit) as exc:
        main(["vcf-check", "-i", "/path/to/nonexistent_file.vcf"])

    assert exc.value.code == 2
    captured = capsys.readouterr()
    assert "Error: File not found" in captured.err


def test_cli_vcf_file_not_found_quiet(capsys: pytest.CaptureFixture):
    """Test CLI vcf-check exits with code 2 and suppresses error output when quiet."""
    with pytest.raises(SystemExit) as exc:
        main(["vcf-check", "-i", "/path/to/nonexistent_file.vcf", "-q"])

    assert exc.value.code == 2
    captured = capsys.readouterr()
    assert captured.err == ""
    assert captured.out == ""


def test_cli_vcf_corrupted_file_exits_with_code_1(capsys: pytest.CaptureFixture):
    """Test CLI vcf-check exits with code 1 when errors are detected."""
    corrupted_file = FIXTURES_DIR / "err_bad_coordinate.vcf"
    with pytest.raises(SystemExit) as exc:
        main(["vcf-check", "-i", str(corrupted_file)])

    assert exc.value.code == 1
    captured = capsys.readouterr()
    assert "Validation FAILED" in captured.err
    assert "INVALID_COORDINATE" in captured.err


def test_cli_vcf_quiet_mode_suppresses_output(capsys: pytest.CaptureFixture):
    """Test CLI vcf-check quiet mode suppresses all stdout and stderr."""
    corrupted_file = FIXTURES_DIR / "err_bad_coordinate.vcf"
    with pytest.raises(SystemExit) as exc:
        main(["vcf-check", "-i", str(corrupted_file), "-q"])

    assert exc.value.code == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_cli_vcf_json_output(capsys: pytest.CaptureFixture):
    """Test CLI vcf-check outputs JSON format to stdout."""
    valid_file = FIXTURES_DIR / "valid_vcf42.vcf"
    with pytest.raises(SystemExit) as exc:
        main(["vcf-check", "-i", str(valid_file), "--json"])

    assert exc.value.code == 0
    captured = capsys.readouterr()
    data = json.loads(captured.out)
    assert data["meta"]["command"] == "vcf-check"
    assert data["validation_gate"]["status"] == "PASSED"
    assert data["validation_gate"]["passed"] is True
    assert data["summary"]["total_records_evaluated"] == 6


def test_cli_vcf_output_json_file(tmp_path: Path):
    """Test CLI vcf-check writes JSON report to a file path."""
    valid_file = FIXTURES_DIR / "valid_vcf42.vcf"
    out_json = tmp_path / "vcf_out.json"

    with pytest.raises(SystemExit) as exc:
        main(["vcf-check", "-i", str(valid_file), "-o", str(out_json)])

    assert exc.value.code == 0
    assert out_json.exists()
    data = json.loads(out_json.read_text())
    assert data["validation_gate"]["passed"] is True


def test_cli_vcf_output_json_with_quiet(tmp_path: Path, capsys: pytest.CaptureFixture):
    """Test CLI vcf-check writes JSON file even when quiet is enabled."""
    valid_file = FIXTURES_DIR / "valid_vcf42.vcf"
    out_json = tmp_path / "vcf_quiet.json"

    with pytest.raises(SystemExit) as exc:
        main(["vcf-check", "-i", str(valid_file), "-o", str(out_json), "-q"])

    assert exc.value.code == 0
    assert out_json.exists()
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_cli_vcf_stdin_pipe(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture):
    """Test CLI vcf-check accepts input piped from stdin ('-')."""
    valid_content = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\trs1\tA\tT\t30\tPASS\t.\n"
    )
    monkeypatch.setattr("sys.stdin", io.StringIO(valid_content))

    with pytest.raises(SystemExit) as exc:
        main(["vcf-check", "-i", "-"])

    assert exc.value.code == 0
    captured = capsys.readouterr()
    assert "Validation PASSED" in captured.out


def test_cli_vcf_gzip_file(capsys: pytest.CaptureFixture):
    """Test CLI vcf-check on a gzip-compressed valid VCF file."""
    valid_gz = FIXTURES_DIR / "valid_vcf42.vcf.gz"
    with pytest.raises(SystemExit) as exc:
        main(["vcf-check", "-i", str(valid_gz)])

    assert exc.value.code == 0
    captured = capsys.readouterr()
    assert "Validation PASSED" in captured.out
    assert "Total Records Evaluated: 6" in captured.out


def test_cli_vcf_max_errors_flag(capsys: pytest.CaptureFixture):
    """Test CLI vcf-check respects --max-errors threshold."""
    corrupted_file = FIXTURES_DIR / "err_bad_coordinate.vcf"
    with pytest.raises(SystemExit) as exc:
        main(["vcf-check", "-i", str(corrupted_file), "--max-errors", "1", "--json"])

    assert exc.value.code == 1
    captured = capsys.readouterr()
    data = json.loads(captured.out)
    assert data["validation_gate"]["total_errors"] == 1
    assert data["meta"]["thresholds"]["max_errors"] == 1


def test_handle_vcf_check_direct_call():
    """Test handle_vcf_check function directly with parsed arguments."""
    from seq_qc_probe.cli import build_parser, handle_vcf_check

    parser = build_parser()
    valid_file = FIXTURES_DIR / "valid_vcf42.vcf"
    args = parser.parse_args(["vcf-check", "-i", str(valid_file), "-q"])
    exit_code = handle_vcf_check(args)
    assert exit_code == 0
