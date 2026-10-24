import io
import json
from pathlib import Path

import pytest

from seq_qc_probe.checker import FASTQChecker
from seq_qc_probe.models import (
    AlignmentErrorKind,
    AlignmentValidationError,
    AlignmentValidationSummary,
    MAPQSummary,
    VCFErrorKind,
    VCFValidationError,
    VCFValidationSummary,
)
from seq_qc_probe.reporter import AlignmentReporter, BaseReporter, QCReporter, VCFReporter
from seq_qc_probe.vcf_validator import VCFValidator

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


def test_base_reporter_abstract_instantiation():
    """Verify BaseReporter cannot be instantiated directly without implementations."""
    with pytest.raises(TypeError):
        BaseReporter()  # type: ignore[abstract]


def test_base_reporter_timestamp():
    """Verify BaseReporter.get_iso_timestamp produces an ISO 8601 UTC string."""
    ts = BaseReporter.get_iso_timestamp()
    assert ts.endswith("Z")
    assert "T" in ts


def test_vcf_reporter_valid_file(tmp_path: Path):
    """Verify VCFReporter formats successful validation results into dictionary and text."""
    vcf_file = FIXTURES_DIR / "valid_vcf42.vcf"
    validator = VCFValidator(vcf_file)
    summary = validator.validate()

    reporter = VCFReporter(summary, input_file=vcf_file, max_errors=10)
    data = reporter.to_dict()

    assert data["schema_version"] == "1.0.0"
    assert data["meta"]["command"] == "vcf-check"
    assert data["meta"]["input_file"] == str(vcf_file)
    assert data["validation_gate"]["status"] == "PASSED"
    assert data["validation_gate"]["passed"] is True
    assert data["validation_gate"]["exit_code"] == 0
    assert data["validation_gate"]["total_errors"] == 0
    assert data["summary"]["total_records_evaluated"] == 6

    # Test text summary
    text_summary = reporter.format_text_summary()
    assert "Validation PASSED" in text_summary
    assert "Total Records Evaluated: 6" in text_summary

    # Test print_terminal_summary to stdout
    stdout_buf = io.StringIO()
    stderr_buf = io.StringIO()
    reporter.print_terminal_summary(stdout=stdout_buf, stderr=stderr_buf)
    assert "Validation PASSED" in stdout_buf.getvalue()
    assert stderr_buf.getvalue() == ""

    # Test JSON output writing
    json_path = tmp_path / "vcf_report.json"
    reporter.write_json(json_path)
    assert json_path.exists()
    loaded = json.loads(json_path.read_text())
    assert loaded["validation_gate"]["passed"] is True


def test_vcf_reporter_failed_file():
    """Verify VCFReporter properly formats error distribution and messages."""
    summary = VCFValidationSummary(
        total_records_evaluated=10,
        total_errors=2,
        error_distribution={
            VCFErrorKind.INVALID_COORDINATE: 1,
            VCFErrorKind.INVALID_BASE: 1,
        },
        errors=[
            VCFValidationError(
                line_number=5,
                error_kind=VCFErrorKind.INVALID_COORDINATE,
                message="Invalid POS value: 0",
                column="POS",
            ),
            VCFValidationError(
                line_number=8,
                error_kind=VCFErrorKind.INVALID_BASE,
                message="Invalid REF base: X",
                column="REF",
            ),
        ],
    )

    reporter = VCFReporter(summary, input_file="bad.vcf")
    data = reporter.to_dict()

    assert data["validation_gate"]["status"] == "FAILED"
    assert data["validation_gate"]["passed"] is False
    assert data["validation_gate"]["exit_code"] == 1
    assert data["validation_gate"]["total_errors"] == 2
    assert data["validation_gate"]["error_distribution"] == {
        "INVALID_COORDINATE": 1,
        "INVALID_BASE": 1,
    }
    assert len(data["validation_gate"]["errors"]) == 2

    text = reporter.format_text_summary()
    assert "Validation FAILED" in text
    assert "INVALID_COORDINATE: 1" in text
    assert "INVALID_BASE: 1" in text
    assert "[Line 5] [INVALID_COORDINATE] Invalid POS value: 0 (Column: POS)" in text
    assert "[Line 8] [INVALID_BASE] Invalid REF base: X (Column: REF)" in text

    # Test print_terminal_summary sends to stderr on failure
    stdout_buf = io.StringIO()
    stderr_buf = io.StringIO()
    reporter.print_terminal_summary(stdout=stdout_buf, stderr=stderr_buf)
    assert stdout_buf.getvalue() == ""
    assert "Validation FAILED" in stderr_buf.getvalue()


def test_vcf_reporter_print_json_stream():
    """Verify VCFReporter.print_json serializes correctly to a given stream."""
    summary = VCFValidationSummary(
        total_records_evaluated=5,
        total_errors=0,
        error_distribution={},
        errors=[],
    )
    reporter = VCFReporter(summary, input_file="test.vcf")
    buf = io.StringIO()
    reporter.print_json(stream=buf)

    payload = json.loads(buf.getvalue())
    assert payload["summary"]["total_records_evaluated"] == 5
    assert payload["validation_gate"]["status"] == "PASSED"


def test_vcf_reporter_failed_without_errors_list():
    """Verify VCFReporter handles failed summary when error list is empty."""
    summary = VCFValidationSummary(
        total_records_evaluated=0,
        total_errors=1,
        error_distribution={},
        errors=[],
    )
    reporter = VCFReporter(summary, input_file="empty.vcf")
    text = reporter.format_text_summary()
    assert "Validation FAILED: Found 1 error(s)" in text
    assert "Errors:" not in text


def test_qc_reporter_failed_output():
    """Verify QCReporter terminal output for failed FASTQ validation."""
    corrupted_file = FIXTURES_DIR / "err_bad_header.fastq"
    checker = FASTQChecker(input_path=corrupted_file, min_quality=30.0)
    result = checker.run()

    reporter = QCReporter(result)
    assert reporter.format_text_summary() == ""

    stdout_buf = io.StringIO()
    stderr_buf = io.StringIO()
    reporter.print_terminal_summary(stdout=stdout_buf, stderr=stderr_buf)
    assert stdout_buf.getvalue() == ""
    assert "Validation FAILED" in stderr_buf.getvalue()


def test_alignment_reporter_valid_file(tmp_path: Path):
    """Verify AlignmentReporter formats a passing SAM validation outcome."""
    mapq = MAPQSummary(
        total_reads=16,
        mapq_zero_count=3,
        mapq_zero_rate=0.1875,
        mean_mapq=41.06,
        mapq_distribution={60: 7, 0: 3, 40: 2, 50: 1, 42: 1, 35: 1, 30: 1},
    )
    summary = AlignmentValidationSummary(
        total_records_evaluated=16,
        total_errors=0,
        error_distribution={},
        errors=[],
        mapq_summary=mapq,
    )

    reporter = AlignmentReporter(summary, input_file="valid.sam", max_errors=10)
    data = reporter.to_dict()

    # Meta
    assert data["schema_version"] == "1.0.0"
    assert data["meta"]["command"] == "sam-check"
    assert data["meta"]["input_file"] == "valid.sam"
    assert data["meta"]["thresholds"]["max_errors"] == 10

    # Validation Gate
    assert data["validation_gate"]["status"] == "PASSED"
    assert data["validation_gate"]["passed"] is True
    assert data["validation_gate"]["exit_code"] == 0
    assert data["validation_gate"]["total_errors"] == 0
    assert data["validation_gate"]["errors"] == []

    # MAPQ Metrics
    assert data["mapq_metrics"]["total_reads"] == 16
    assert data["mapq_metrics"]["mapq_zero_count"] == 3
    assert data["mapq_metrics"]["mapq_zero_rate"] == 0.1875
    assert data["mapq_metrics"]["multimapping_rate"] == 0.1875
    assert data["mapq_metrics"]["mean_mapq"] == 41.06
    assert data["mapq_metrics"]["mapq_distribution"][60] == 7

    # Text summary
    text = reporter.format_text_summary()
    assert "Validation PASSED: SAM alignment file is valid." in text
    assert "Total Records Evaluated: 16" in text
    assert "MAPQ & Multimapping Metrics:" in text
    assert "Mean MAPQ: 41.06" in text
    assert "MAPQ = 0 (Multimapped): 3 (18.75%)" in text

    # Terminal summary to stdout
    stdout_buf = io.StringIO()
    stderr_buf = io.StringIO()
    reporter.print_terminal_summary(stdout=stdout_buf, stderr=stderr_buf)
    assert "Validation PASSED" in stdout_buf.getvalue()
    assert stderr_buf.getvalue() == ""

    # JSON export
    json_path = tmp_path / "sam_report.json"
    reporter.write_json(json_path)
    assert json_path.exists()
    payload = json.loads(json_path.read_text())
    assert payload["validation_gate"]["passed"] is True
    assert payload["mapq_metrics"]["mapq_distribution"]["60"] == 7


def test_alignment_reporter_failed_file():
    """Verify AlignmentReporter properly formats error breakdown and individual errors."""
    mapq = MAPQSummary(
        total_reads=8,
        mapq_zero_count=2,
        mapq_zero_rate=0.25,
        mean_mapq=35.62,
        mapq_distribution={60: 4, 0: 2, 50: 1, 30: 1},
    )
    summary = AlignmentValidationSummary(
        total_records_evaluated=8,
        total_errors=2,
        error_distribution={
            AlignmentErrorKind.CIGAR_LENGTH_MISMATCH: 1,
            AlignmentErrorKind.INVALID_CIGAR: 1,
        },
        errors=[
            AlignmentValidationError(
                line_number=7,
                error_kind=AlignmentErrorKind.CIGAR_LENGTH_MISMATCH,
                message="Expected 100 bp from CIGAR '100M', but got 92 bp",
                qname="read02_short",
                cigar="100M",
                expected_length=100,
                actual_length=92,
                column="CIGAR",
            ),
            AlignmentValidationError(
                line_number=12,
                error_kind=AlignmentErrorKind.INVALID_CIGAR,
                message="Invalid CIGAR syntax: '100Z'",
                qname="read05_bad",
                cigar="100Z",
                column="CIGAR",
            ),
        ],
        mapq_summary=mapq,
    )

    reporter = AlignmentReporter(summary, input_file="corrupted.sam")
    data = reporter.to_dict()

    assert data["validation_gate"]["status"] == "FAILED"
    assert data["validation_gate"]["passed"] is False
    assert data["validation_gate"]["exit_code"] == 1
    assert data["validation_gate"]["total_errors"] == 2
    assert data["validation_gate"]["error_distribution"] == {
        "CIGAR_LENGTH_MISMATCH": 1,
        "INVALID_CIGAR": 1,
    }
    assert len(data["validation_gate"]["errors"]) == 2

    text = reporter.format_text_summary()
    assert "Validation FAILED: Found 2 error(s) across 8 record(s)." in text
    assert "MAPQ = 0 (Multimapped): 2 (25.00%)" in text
    assert "CIGAR_LENGTH_MISMATCH: 1" in text
    assert "INVALID_CIGAR: 1" in text
    assert (
        "[Line 7] [CIGAR_LENGTH_MISMATCH] (Read: read02_short, CIGAR: 100M, Column: CIGAR)" in text
    )

    # Verify terminal summary goes to stderr on failure
    stdout_buf = io.StringIO()
    stderr_buf = io.StringIO()
    reporter.print_terminal_summary(stdout=stdout_buf, stderr=stderr_buf)
    assert stdout_buf.getvalue() == ""
    assert "Validation FAILED" in stderr_buf.getvalue()


def test_alignment_reporter_print_json_stream():
    """Verify AlignmentReporter.print_json serializes correctly to a given stream."""
    summary = AlignmentValidationSummary(
        total_records_evaluated=10,
        total_errors=0,
        error_distribution={},
        errors=[],
    )
    reporter = AlignmentReporter(summary, input_file="stream.sam")
    buf = io.StringIO()
    reporter.print_json(stream=buf)

    payload = json.loads(buf.getvalue())
    assert payload["summary"]["total_records_evaluated"] == 10
    assert payload["validation_gate"]["status"] == "PASSED"


def test_alignment_reporter_failed_without_errors_list():
    """Verify AlignmentReporter handles failed summary when error list is empty."""
    summary = AlignmentValidationSummary(
        total_records_evaluated=0,
        total_errors=1,
        error_distribution={},
        errors=[],
    )
    reporter = AlignmentReporter(summary, input_file="empty.sam")
    text = reporter.format_text_summary()
    assert "Validation FAILED: Found 1 error(s)" in text
    assert "Errors:" not in text
