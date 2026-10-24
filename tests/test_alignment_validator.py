"""Unit tests for AlignmentValidator, SAMRecordParser, and CIGARValidator."""

import io
from pathlib import Path

from seq_qc_probe.alignment_validator import AlignmentValidator, CIGARValidator, validate_cigar
from seq_qc_probe.models import AlignmentErrorKind

FIXTURES_DIR = Path(__file__).parent / "fixtures"


def test_valid_sam_parses_completely():
    """Verify that a compliant SAM file passes validation with zero errors and accurate MAPQ."""
    sam_file = FIXTURES_DIR / "valid_alignment.sam"
    validator = AlignmentValidator(sam_file)
    summary = validator.validate()

    assert summary.is_valid is True
    assert summary.total_errors == 0
    assert summary.total_records_evaluated == 16
    assert len(summary.errors) == 0

    # MAPQ assertions
    mapq = summary.mapq_summary
    assert mapq.total_reads == 16
    assert mapq.mapq_zero_count == 3
    assert mapq.mapq_zero_rate == 0.1875
    assert mapq.multimapping_rate == 0.1875
    assert mapq.mean_mapq == 41.06
    assert mapq.mapq_distribution[60] == 7
    assert mapq.mapq_distribution[0] == 3


def test_valid_sam_gzip_parsing():
    """Verify that gzip-compressed SAM files (.sam.gz) are parsed transparently."""
    sam_gz = FIXTURES_DIR / "valid_alignment.sam.gz"
    validator = AlignmentValidator(sam_gz)
    summary = validator.validate()

    assert summary.is_valid is True
    assert summary.total_errors == 0
    assert summary.total_records_evaluated == 16


def test_err_cigar_mismatch_sam():
    """Verify that CIGAR length mismatches are detected accurately with correct metrics."""
    sam_file = FIXTURES_DIR / "err_cigar_mismatch.sam"
    validator = AlignmentValidator(sam_file)
    summary = validator.validate()

    assert summary.is_valid is False
    assert summary.total_errors == 6
    assert summary.total_records_evaluated == 8
    assert summary.error_distribution[AlignmentErrorKind.CIGAR_LENGTH_MISMATCH] == 6

    # Verify MAPQ metrics are collected from all records regardless of CIGAR status
    mapq = summary.mapq_summary
    assert mapq.total_reads == 8
    assert mapq.mapq_zero_count == 2
    assert mapq.mapq_zero_rate == 0.25


def test_err_invalid_cigar_syntax():
    """Verify that syntactically invalid CIGAR strings trigger INVALID_CIGAR errors."""
    sam_file = FIXTURES_DIR / "err_invalid_cigar.sam"
    validator = AlignmentValidator(sam_file)
    summary = validator.validate()

    assert summary.is_valid is False
    assert summary.total_errors == 4
    assert summary.total_records_evaluated == 5
    assert summary.error_distribution[AlignmentErrorKind.INVALID_CIGAR] == 4


def test_err_malformed_sam():
    """Verify that malformed records, bad POS, and out-of-range MAPQ are detected."""
    sam_file = FIXTURES_DIR / "err_malformed_sam.sam"
    validator = AlignmentValidator(sam_file)
    summary = validator.validate()

    assert summary.is_valid is False
    assert summary.total_errors == 4
    assert summary.total_records_evaluated == 5
    assert summary.error_distribution[AlignmentErrorKind.MALFORMED_RECORD] == 1
    assert summary.error_distribution[AlignmentErrorKind.INVALID_MAPQ] == 2
    assert summary.error_distribution[AlignmentErrorKind.INVALID_POSITION] == 1


def test_empty_sam_file_detection(tmp_path: Path):
    """Verify that an empty SAM file triggers PREMATURE_EOF."""
    empty_file = tmp_path / "empty.sam"
    empty_file.write_text("")

    validator = AlignmentValidator(empty_file)
    summary = validator.validate()

    assert summary.is_valid is False
    assert summary.error_distribution[AlignmentErrorKind.PREMATURE_EOF] == 1
    assert summary.errors[0].line_number == 1
    assert "empty" in summary.errors[0].message.lower()


def test_header_after_alignment_record_fails(tmp_path: Path):
    """Verify that header lines '@' appearing after alignment records trigger INVALID_HEADER."""
    sam_file = tmp_path / "misplaced_header.sam"
    content = (
        "@HD\tVN:1.6\n"
        "read01\t0\tchr1\t100\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tFFFFFFFFFF\n"
        "@CO\tlate header comment\n"
    )
    sam_file.write_text(content)

    validator = AlignmentValidator(sam_file)
    summary = validator.validate()

    assert summary.is_valid is False
    assert summary.error_distribution[AlignmentErrorKind.INVALID_HEADER] == 1
    assert summary.errors[0].line_number == 3


def test_mismatched_seq_qual_lengths(tmp_path: Path):
    """Verify that SEQ length != QUAL length triggers MISMATCHED_LENGTHS."""
    sam_file = tmp_path / "seq_qual_mismatch.sam"
    content = "@HD\tVN:1.6\nread01\t0\tchr1\t100\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tFFF\n"
    sam_file.write_text(content)

    validator = AlignmentValidator(sam_file)
    summary = validator.validate()

    assert summary.is_valid is False
    assert summary.error_distribution[AlignmentErrorKind.MISMATCHED_LENGTHS] == 1


def test_invalid_flag_detection(tmp_path: Path):
    """Verify that negative or non-integer FLAG triggers INVALID_FLAG."""
    sam_file = tmp_path / "bad_flag.sam"
    content = (
        "@HD\tVN:1.6\n"
        "read01\t-16\tchr1\t100\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tFFFFFFFFFF\n"
        "read02\tNOT_INT\tchr1\t200\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tFFFFFFFFFF\n"
    )
    sam_file.write_text(content)

    validator = AlignmentValidator(sam_file)
    summary = validator.validate()

    assert summary.is_valid is False
    assert summary.error_distribution[AlignmentErrorKind.INVALID_FLAG] == 2


def test_max_errors_threshold_aborts_early():
    """Verify that validation stops when max_errors threshold is reached."""
    sam_file = FIXTURES_DIR / "err_cigar_mismatch.sam"
    validator = AlignmentValidator(sam_file, max_errors=2)
    summary = validator.validate()

    assert summary.is_valid is False
    assert summary.total_errors == 2
    assert len(summary.errors) == 2


def test_stdin_pipe_support(monkeypatch):
    """Verify that passing '-' reads from sys.stdin transparently."""
    content = "@HD\tVN:1.6\nread01\t0\tchr1\t100\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tFFFFFFFFFF\n"
    fake_stdin = io.StringIO(content)
    monkeypatch.setattr("sys.stdin", fake_stdin)

    validator = AlignmentValidator("-")
    summary = validator.validate()

    assert summary.is_valid is True
    assert summary.total_records_evaluated == 1


def test_validate_cigar_direct_function():
    """Verify standalone validate_cigar function logic on diverse operations."""
    # Test class method directly
    res_cls = CIGARValidator.validate("10M", 10)
    assert res_cls.is_valid is True

    # Matches and insertions consume read bases: 10M + 2I + 8M = 20 bp
    res = validate_cigar("10M2I8M", 20)
    assert res.is_valid is True
    assert res.expected_length == 20

    # Deletions do not consume read bases: 10M + 5D + 10M = 20 bp
    res = validate_cigar("10M5D10M", 20)
    assert res.is_valid is True
    assert res.expected_length == 20

    # Hard clipping does not consume read bases: 15H + 10M = 10 bp
    res = validate_cigar("15H10M", 10)
    assert res.is_valid is True
    assert res.expected_length == 10

    # Soft clipping consumes read bases: 5S + 10M + 5S = 20 bp
    res = validate_cigar("5S10M5S", 20)
    assert res.is_valid is True
    assert res.expected_length == 20

    # Unavailable CIGAR '*'
    res = validate_cigar("*", 50)
    assert res.is_valid is True

    # Bad syntax
    res = validate_cigar("20Z", 20)
    assert res.is_valid is False
    assert res.error_kind == AlignmentErrorKind.INVALID_CIGAR


def test_validate_record_direct_method():
    """Verify validate_record direct invocation updates validator summary."""
    validator = AlignmentValidator("dummy.sam")
    record = "read01\t0\tchr1\t100\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tFFFFFFFFFF"
    validator.validate_record(record, line_number=5)

    assert validator.summary.total_records_evaluated == 1
    assert validator.summary.is_valid is True
    assert validator.summary.mapq_summary.total_reads == 1


def test_invalid_pnext_and_tlen_detection(tmp_path: Path):
    """Verify that negative PNEXT or non-integer TLEN trigger MALFORMED_RECORD."""
    sam_file = tmp_path / "bad_pnext_tlen.sam"
    content = (
        "@HD\tVN:1.6\n"
        "read01\t0\tchr1\t100\t60\t10M\t*\t-5\t0\tAAAAAAAAAA\tFFFFFFFFFF\n"
        "read02\t0\tchr1\t200\t60\t10M\t*\t0\tNOT_INT\tAAAAAAAAAA\tFFFFFFFFFF\n"
    )
    sam_file.write_text(content)

    validator = AlignmentValidator(sam_file)
    summary = validator.validate()

    assert summary.is_valid is False
    assert summary.error_distribution[AlignmentErrorKind.MALFORMED_RECORD] == 2


def test_empty_lines_in_sam_are_ignored(tmp_path: Path):
    """Verify that blank lines are safely skipped during processing."""
    sam_file = tmp_path / "with_blank_lines.sam"
    content = "@HD\tVN:1.6\n\nread01\t0\tchr1\t100\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tFFFFFFFFFF\n\n"
    sam_file.write_text(content)

    validator = AlignmentValidator(sam_file)
    summary = validator.validate()

    assert summary.is_valid is True
    assert summary.total_records_evaluated == 1


def test_empty_qname_and_rname_detection(tmp_path: Path):
    """Verify that empty QNAME or RNAME trigger MISSING_FIELD error."""
    sam_file = tmp_path / "empty_qname_rname.sam"
    content = (
        "@HD\tVN:1.6\n"
        "\t0\tchr1\t100\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tFFFFFFFFFF\n"
        "read02\t0\t\t100\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tFFFFFFFFFF\n"
    )
    sam_file.write_text(content)

    validator = AlignmentValidator(sam_file)
    summary = validator.validate()

    assert summary.is_valid is False
    assert summary.error_distribution[AlignmentErrorKind.MISSING_FIELD] == 2
    assert summary.errors[0].column == "QNAME"
    assert summary.errors[1].column == "RNAME"


def test_whitespace_only_sam_file_triggers_eof(tmp_path: Path):
    """Verify that a file containing only newlines and spaces triggers PREMATURE_EOF."""
    sam_file = tmp_path / "whitespace_only.sam"
    sam_file.write_text("\n\n   \n\n")

    validator = AlignmentValidator(sam_file)
    summary = validator.validate()

    assert summary.is_valid is False
    assert summary.error_distribution[AlignmentErrorKind.PREMATURE_EOF] == 1
