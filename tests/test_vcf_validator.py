import io
from pathlib import Path

import pytest

from seq_qc_probe.models import VCFErrorKind
from seq_qc_probe.vcf_validator import VCFValidator

FIXTURES_DIR = Path(__file__).parent / "fixtures"


def test_valid_vcf42_parses_completely():
    """Verify that a compliant VCF 4.2 file passes validation with zero errors."""
    vcf_file = FIXTURES_DIR / "valid_vcf42.vcf"
    validator = VCFValidator(vcf_file)
    summary = validator.validate()

    assert summary.is_valid is True
    assert summary.total_errors == 0
    assert summary.total_records_evaluated == 6
    assert len(summary.errors) == 0


def test_valid_vcf42_gzip_parsing():
    """Verify that gzip-compressed VCF 4.2 files are parsed transparently."""
    vcf_gz = FIXTURES_DIR / "valid_vcf42.vcf.gz"
    validator = VCFValidator(vcf_gz)
    summary = validator.validate()

    assert summary.is_valid is True
    assert summary.total_errors == 0
    assert summary.total_records_evaluated == 6


def test_missing_fileformat_detection():
    """Verify that a VCF missing ##fileformat=VCFv4.2 as line 1 is flagged."""
    vcf_file = FIXTURES_DIR / "err_missing_fileformat.vcf"
    validator = VCFValidator(vcf_file)
    summary = validator.validate()

    assert summary.is_valid is False
    assert summary.error_distribution[VCFErrorKind.MISSING_FILEFORMAT] == 1
    assert summary.errors[0].line_number == 1
    assert "fileformat" in summary.errors[0].message.lower()


def test_missing_column_header_detection():
    """Verify that a header missing mandatory columns triggers INVALID_HEADER_COLUMNS."""
    vcf_file = FIXTURES_DIR / "err_missing_column.vcf"
    validator = VCFValidator(vcf_file)
    summary = validator.validate()

    assert summary.is_valid is False
    assert summary.error_distribution[VCFErrorKind.INVALID_HEADER_COLUMNS] >= 1


def test_bad_coordinates_detection():
    """Verify that POS <= 0 or non-integer POS values are detected as INVALID_COORDINATE."""
    vcf_file = FIXTURES_DIR / "err_bad_coordinate.vcf"
    validator = VCFValidator(vcf_file)
    summary = validator.validate()

    assert summary.is_valid is False
    assert summary.error_distribution[VCFErrorKind.INVALID_COORDINATE] == 3
    assert summary.total_records_evaluated == 3

    error_lines = [e.line_number for e in summary.errors]
    assert error_lines == [3, 4, 5]


def test_bad_ref_alt_bases_detection():
    """Verify that invalid IUPAC bases or dot in REF are flagged as INVALID_BASE."""
    vcf_file = FIXTURES_DIR / "err_bad_ref_alt.vcf"
    validator = VCFValidator(vcf_file)
    summary = validator.validate()

    assert summary.is_valid is False
    assert summary.error_distribution[VCFErrorKind.INVALID_BASE] == 4
    assert summary.total_records_evaluated == 4


def test_malformed_meta_lines_detection():
    """Verify that syntax errors in meta lines (missing '=', unclosed '<') are flagged."""
    vcf_file = FIXTURES_DIR / "err_malformed_meta.vcf"
    validator = VCFValidator(vcf_file)
    summary = validator.validate()

    assert summary.is_valid is False
    assert summary.error_distribution[VCFErrorKind.MALFORMED_LINE] == 2


def test_empty_file_detection(tmp_path: Path):
    """Verify that an empty file produces a PREMATURE_EOF error."""
    empty_vcf = tmp_path / "empty.vcf"
    empty_vcf.write_text("")

    validator = VCFValidator(empty_vcf)
    summary = validator.validate()

    assert summary.is_valid is False
    assert summary.error_distribution[VCFErrorKind.PREMATURE_EOF] == 1


def test_max_errors_threshold_aborts_early():
    """Verify that validation aborts early when max_errors threshold is reached."""
    vcf_file = FIXTURES_DIR / "err_bad_ref_alt.vcf"
    validator = VCFValidator(vcf_file, max_errors=2)
    summary = validator.validate()

    assert summary.total_errors == 2
    assert len(summary.errors) == 2


def test_format_without_samples_fails(tmp_path: Path):
    """Verify that a 9-column header with FORMAT but no samples is invalid."""
    bad_vcf = tmp_path / "format_no_samples.vcf"
    bad_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"
        "chr1\t100\trs1\tA\tT\t30\tPASS\t.\tGT\n"
    )
    validator = VCFValidator(bad_vcf)
    summary = validator.validate()

    assert summary.is_valid is False
    assert summary.error_distribution[VCFErrorKind.INVALID_HEADER_COLUMNS] >= 1


def test_column_count_mismatch(tmp_path: Path):
    """Verify that records with differing column counts produce INVALID_RECORD."""
    bad_vcf = tmp_path / "mismatch_cols.vcf"
    bad_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\trs1\tA\tT\n"  # Only 5 columns
    )
    validator = VCFValidator(bad_vcf)
    summary = validator.validate()

    assert summary.is_valid is False
    assert summary.error_distribution[VCFErrorKind.INVALID_RECORD] == 1


def test_meta_line_after_header_fails(tmp_path: Path):
    """Verify that meta-lines appearing after #CHROM header trigger MALFORMED_LINE."""
    bad_vcf = tmp_path / "meta_after_header.vcf"
    bad_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        "chr1\t100\trs1\tA\tT\t30\tPASS\t.\n"
    )
    validator = VCFValidator(bad_vcf)
    summary = validator.validate()

    assert summary.is_valid is False
    assert summary.error_distribution[VCFErrorKind.MALFORMED_LINE] >= 1


def test_stdin_pipe_support(monkeypatch: pytest.MonkeyPatch):
    """Verify that VCF validation works when reading from standard input ('-')."""
    valid_content = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\trs1\tA\tT\t30\tPASS\t.\n"
    )
    monkeypatch.setattr("sys.stdin", io.StringIO(valid_content))

    validator = VCFValidator("-")
    summary = validator.validate()

    assert summary.is_valid is True
    assert summary.total_records_evaluated == 1


def test_meta_only_file_flags_missing_header(tmp_path: Path):
    """Verify that a file with only meta lines triggers missing header error."""
    vcf = tmp_path / "meta_only.vcf"
    vcf.write_text('##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description="Passed">\n')
    s = VCFValidator(vcf).validate()
    assert s.is_valid is False
    assert s.error_distribution[VCFErrorKind.INVALID_HEADER_COLUMNS] >= 1


def test_first_line_header_without_fileformat(tmp_path: Path):
    """Verify that a file starting directly with #CHROM handles header and missing format."""
    vcf = tmp_path / "no_format_header.vcf"
    vcf.write_text(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\trs1\tA\tT\t30\tPASS\t.\n"
    )
    s = VCFValidator(vcf).validate()
    assert s.is_valid is False
    assert s.error_distribution[VCFErrorKind.MISSING_FILEFORMAT] == 1
    assert s.total_records_evaluated == 1


def test_first_line_record_without_meta(tmp_path: Path):
    """Verify that a file starting directly with a record handles data line."""
    vcf = tmp_path / "record_first.vcf"
    vcf.write_text("chr1\t100\trs1\tA\tT\t30\tPASS\t.\n")
    s = VCFValidator(vcf).validate()
    assert s.is_valid is False
    assert s.total_records_evaluated == 1


def test_empty_key_in_meta(tmp_path: Path):
    """Verify that a meta line like ##=val produces a MALFORMED_LINE error."""
    vcf = tmp_path / "empty_meta_key.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##=bad_val\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    s = VCFValidator(vcf).validate()
    assert s.is_valid is False
    assert s.error_distribution[VCFErrorKind.MALFORMED_LINE] >= 1


def test_empty_sample_identifier_fails(tmp_path: Path):
    """Verify that whitespace/empty sample identifier is flagged."""
    vcf = tmp_path / "empty_sample.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t   \n"
    )
    s = VCFValidator(vcf).validate()
    assert s.is_valid is False
    assert s.error_distribution[VCFErrorKind.INVALID_HEADER_COLUMNS] >= 1


def test_sample_without_format_tag(tmp_path: Path):
    """Verify that 10th column present without FORMAT at 9th column fails."""
    vcf = tmp_path / "sample_no_format.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tNOTFORMAT\tSAMPLE1\n"
    )
    s = VCFValidator(vcf).validate()
    assert s.is_valid is False
    assert s.error_distribution[VCFErrorKind.INVALID_HEADER_COLUMNS] >= 1


def test_empty_lines_and_comments_in_data_section(tmp_path: Path):
    """Verify that empty lines or stray # in data section trigger MALFORMED_LINE."""
    vcf = tmp_path / "stray_lines.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\trs1\tA\tT\t30\tPASS\t.\n"
        "\n"
        "#comment_in_data\n"
        "chr1\t200\trs2\tC\tG\t30\tPASS\t.\n"
    )
    s = VCFValidator(vcf).validate()
    assert s.is_valid is False
    assert s.error_distribution[VCFErrorKind.MALFORMED_LINE] == 2


def test_empty_alt_allele_triggers_invalid_base(tmp_path: Path):
    """Verify that empty allele in ALT (e.g. 'A,,C') triggers INVALID_BASE."""
    vcf = tmp_path / "empty_alt_allele.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\trs1\tA\tA,,C\t30\tPASS\t.\n"
    )
    s = VCFValidator(vcf).validate()
    assert s.is_valid is False
    assert s.error_distribution[VCFErrorKind.INVALID_BASE] >= 1
