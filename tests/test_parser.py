import gzip
from pathlib import Path

import pytest

from seq_qc_probe.models import FASTQErrorKind, FASTQRecord, FASTQValidationError
from seq_qc_probe.parser import parse_and_validate_fastq

FIXTURES_DIR = Path(__file__).parent / "fixtures"


def test_valid_fastq_parses_completely():
    """Test that a valid FASTQ file is parsed correctly and yields FASTQRecord objects."""
    file_path = FIXTURES_DIR / "valid.fastq"
    results = list(parse_and_validate_fastq(str(file_path)))

    assert len(results) == 2
    assert all(isinstance(item, FASTQRecord) for item in results)

    first_record = results[0]
    assert isinstance(first_record, FASTQRecord)
    assert first_record.header == "@READ_01_CLEAN"
    assert first_record.sequence == "AGCTTAGCTAGCTACCTAGCTAGCT"
    assert first_record.separator == "+"
    assert first_record.quality == "IIIIIIIIIIIIIIIIIIIIIIIII"
    assert first_record.line_number == 1


@pytest.mark.parametrize(
    "fixture_name, expected_error_kind, expected_line",
    [
        ("err_bad_header.fastq", FASTQErrorKind.INVALID_HEADER, 1),
        ("err_bad_separator.fastq", FASTQErrorKind.INVALID_SEPARATOR, 3),
        ("err_length_mismatch.fastq", FASTQErrorKind.MISMATCHED_LENGTHS, 8),
        ("err_truncated_eof.fastq", FASTQErrorKind.PREMATURE_EOF, 8),
        ("err_truncated_at_seq.fastq", FASTQErrorKind.PREMATURE_EOF, 7),
        ("err_empty_line.fastq", FASTQErrorKind.EMPTY_LINE, 5),
    ],
)
def test_synthetic_error_modes(
    fixture_name: str,
    expected_error_kind: FASTQErrorKind,
    expected_line: int,
):
    """Test that synthetic error modes produce the expected error types and line numbers."""
    file_path = FIXTURES_DIR / fixture_name
    results = list(parse_and_validate_fastq(str(file_path)))

    errors = [item for item in results if isinstance(item, FASTQValidationError)]
    assert len(errors) >= 1

    target_error = errors[0]
    assert target_error.error_kind == expected_error_kind
    assert target_error.line_number == expected_line


def test_invalid_nucleotide_sequence(tmp_path: Path):
    """Test that a FASTQ file with an invalid base yields an INVALID_SEQUENCE error."""
    corrupted_content = (
        "@READ_INVALID_BASE\nAGCTXAGCTAGCTACCTAGCTAGCT\n+\nIIIIIIIIIIIIIIIIIIIIIIIII\n"
    )
    test_file = tmp_path / "corrupt_base.fastq"
    test_file.write_text(corrupted_content, encoding="utf-8")

    results = list(parse_and_validate_fastq(str(test_file)))
    errors = [item for item in results if isinstance(item, FASTQValidationError)]

    assert len(errors) == 1
    assert errors[0].error_kind == FASTQErrorKind.INVALID_SEQUENCE
    assert errors[0].line_number == 2
    assert errors[0].record is not None
    assert errors[0].record.header == "@READ_INVALID_BASE"


def test_empty_sequence_record_yields_error(tmp_path: Path):
    """Test that an empty sequence yields an INVALID_SEQUENCE error."""
    content = "@EMPTY_SEQ\n\n+\n\n"
    test_file = tmp_path / "empty_seq.fastq"
    test_file.write_text(content, encoding="utf-8")

    results = list(parse_and_validate_fastq(str(test_file)))
    errors = [item for item in results if isinstance(item, FASTQValidationError)]

    assert len(errors) >= 1
    assert errors[0].error_kind == FASTQErrorKind.INVALID_SEQUENCE


def test_invalid_phred_quality_character(tmp_path: Path):
    """Test that quality score characters below ASCII 33 or above 126 yield INVALID_QUALITY."""
    # ASCII 32 is space (' ')
    content = "@BAD_QUAL\nACGT\n+\nII I\n"
    test_file = tmp_path / "bad_qual.fastq"
    test_file.write_text(content, encoding="utf-8")

    results = list(parse_and_validate_fastq(str(test_file)))
    errors = [item for item in results if isinstance(item, FASTQValidationError)]

    assert len(errors) == 1
    assert errors[0].error_kind == FASTQErrorKind.INVALID_QUALITY
    assert errors[0].line_number == 4
    assert errors[0].record is not None


def test_empty_file_yields_nothing(tmp_path: Path):
    """Test that an empty FASTQ file yields no records or errors without mutating fixtures."""
    empty_file = tmp_path / "empty.fastq"
    empty_file.touch()

    results = list(parse_and_validate_fastq(str(empty_file)))
    assert results == []


def test_gzip_fastq_parsing(tmp_path: Path):
    """Test that gzipped FASTQ files (.fastq.gz) are decompressed and parsed transparently."""
    content = b"@GZ_READ_1\nACGT\n+\nIIII\n@GZ_READ_2\nTGCA\n+\nJJJJ\n"
    gz_file = tmp_path / "sample.fastq.gz"
    with gzip.open(gz_file, "wb") as f:
        f.write(content)

    results = list(parse_and_validate_fastq(str(gz_file)))
    assert len(results) == 2
    assert all(isinstance(r, FASTQRecord) for r in results)
    assert results[0].header == "@GZ_READ_1"
    assert results[1].header == "@GZ_READ_2"


def test_mixed_stream_recovers_after_corruption(tmp_path: Path):
    """Test that stream continues parsing subsequent records after an error."""
    mixed_content = (
        "@READ_VALID_1\nACGT\n+\nIIII\n"
        "@READ_CORRUPT_LEN\nACGT\n+\nII\n"
        "@READ_VALID_2\nTGCA\n+\nJJJJ\n"
    )
    mixed_file = tmp_path / "mixed.fastq"
    mixed_file.write_text(mixed_content, encoding="utf-8")

    results = list(parse_and_validate_fastq(str(mixed_file)))

    assert len(results) == 3
    assert isinstance(results[0], FASTQRecord)
    assert results[0].header == "@READ_VALID_1"

    assert isinstance(results[1], FASTQValidationError)
    assert results[1].error_kind == FASTQErrorKind.MISMATCHED_LENGTHS

    assert isinstance(results[2], FASTQRecord)
    assert results[2].header == "@READ_VALID_2"


def test_real_variable_crlf_parsing():
    """Test parsing a real-world style file with CRLF line endings and variable lengths."""
    file_path = FIXTURES_DIR / "real_variable_crlf.fastq"
    results = list(parse_and_validate_fastq(str(file_path)))

    assert len(results) == 3
    assert all(isinstance(r, FASTQRecord) for r in results)

    # Read 1: 30 bp
    assert len(results[0].sequence) == 30
    assert results[0].header.startswith("@M00123:1:")

    # Read 2: 45 bp with lowercase soft-masked bases
    assert len(results[1].sequence) == 45
    assert results[1].sequence.startswith("acgt")

    # Read 3: 25 bp, ended at EOF without a trailing newline
    assert len(results[2].sequence) == 25
    assert results[2].quality == "JJJJJJJJJJJJJJJJJJJJJJJJJ"


def test_real_illumina_pe_gzip_fixture():
    """Test parsing a physical gzipped FASTQ fixture with Illumina CASAVA 1.8+ headers."""
    file_path = FIXTURES_DIR / "real_illumina_pe.fastq.gz"
    results = list(parse_and_validate_fastq(str(file_path)))

    assert len(results) == 2
    assert all(isinstance(r, FASTQRecord) for r in results)

    # Read 1 has uncalled N bases and Q2 '#'
    assert "NNNN" in results[0].sequence
    assert "####" in results[0].quality
    assert results[0].header == "@M00967:43:000000000-A3J0G:1:1101:14495:1817 1:N:0:ATCACG"
