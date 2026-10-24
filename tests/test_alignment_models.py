"""Unit tests for alignment data models in seq_qc_probe.models."""

from seq_qc_probe.models import (
    AlignmentErrorKind,
    AlignmentRecord,
    AlignmentValidationError,
    AlignmentValidationSummary,
    CIGARValidationResult,
    MAPQSummary,
)


def test_alignment_error_kind_enum():
    """Verify that AlignmentErrorKind contains the required error types."""
    assert AlignmentErrorKind.CIGAR_LENGTH_MISMATCH is not None
    assert AlignmentErrorKind.INVALID_CIGAR is not None
    assert AlignmentErrorKind.MALFORMED_RECORD is not None
    assert AlignmentErrorKind.MISSING_FIELD is not None
    assert AlignmentErrorKind.INVALID_MAPQ is not None
    assert AlignmentErrorKind.INVALID_FLAG is not None
    assert AlignmentErrorKind.INVALID_POSITION is not None
    assert AlignmentErrorKind.INVALID_REFERENCE_NAME is not None
    assert AlignmentErrorKind.INVALID_HEADER is not None
    assert AlignmentErrorKind.MISMATCHED_LENGTHS is not None
    assert AlignmentErrorKind.PREMATURE_EOF is not None
    assert AlignmentErrorKind.INVALID_ALIGNMENT is not None


def test_alignment_validation_error_attributes():
    """Verify that AlignmentValidationError stores structural error attributes accurately."""
    err = AlignmentValidationError(
        line_number=10,
        error_kind=AlignmentErrorKind.CIGAR_LENGTH_MISMATCH,
        message="Length mismatch",
        record="read01\t0\tchr1\t100\t60\t100M\t*\t0\t0\tACGT\tFFFF",
        qname="read01",
        cigar="100M",
        expected_length=100,
        actual_length=95,
        column="CIGAR",
    )

    assert err.line_number == 10
    assert err.error_kind == AlignmentErrorKind.CIGAR_LENGTH_MISMATCH
    assert err.message == "Length mismatch"
    assert err.qname == "read01"
    assert err.cigar == "100M"
    assert err.expected_length == 100
    assert err.actual_length == 95
    assert err.column == "CIGAR"
    assert err.record is not None


def test_cigar_validation_result():
    """Verify CIGARValidationResult represents CIGAR calculation outcomes properly."""
    res_valid = CIGARValidationResult(
        is_valid=True,
        cigar="100M",
        expected_length=100,
        actual_length=100,
        operations=[(100, "M")],
    )
    assert res_valid.is_valid is True
    assert res_valid.expected_length == 100
    assert res_valid.actual_length == 100
    assert res_valid.operations == [(100, "M")]
    assert res_valid.error_kind is None

    res_err = CIGARValidationResult(
        is_valid=False,
        cigar="100M",
        expected_length=100,
        actual_length=80,
        operations=[(100, "M")],
        error_kind=AlignmentErrorKind.CIGAR_LENGTH_MISMATCH,
        message="Expected 100 bp from CIGAR '100M', but sequence length is 80 bp",
    )
    assert res_err.is_valid is False
    assert res_err.expected_length == 100
    assert res_err.actual_length == 80
    assert res_err.error_kind == AlignmentErrorKind.CIGAR_LENGTH_MISMATCH


def test_mapq_summary_accumulation_and_metrics():
    """Verify MAPQSummary accumulates reads and calculates multimapping & mean scores."""
    mapq_summary = MAPQSummary()
    assert mapq_summary.total_reads == 0
    assert mapq_summary.mapq_zero_count == 0
    assert mapq_summary.mapq_zero_rate == 0.0
    assert mapq_summary.multimapping_rate == 0.0
    assert mapq_summary.mean_mapq == 0.0

    # Add reads: 2x MAPQ 60, 2x MAPQ 0
    mapq_summary.add_read(60)
    mapq_summary.add_read(60)
    mapq_summary.add_read(0)
    mapq_summary.add_read(0)

    assert mapq_summary.total_reads == 4
    assert mapq_summary.mapq_zero_count == 2
    assert mapq_summary.mapq_zero_rate == 0.5
    assert mapq_summary.multimapping_rate == 0.5
    assert mapq_summary.mean_mapq == 30.0
    assert mapq_summary.mapq_distribution == {60: 2, 0: 2}


def test_mapq_summary_empty_recalculate():
    """Verify recalculate on empty MAPQSummary resets rates to 0.0 safely."""
    mapq = MAPQSummary()
    mapq.recalculate()
    assert mapq.mapq_zero_rate == 0.0
    assert mapq.mean_mapq == 0.0


def test_mapq_summary_manual_population_recalculate():
    """Verify recalculate with manually populated distribution computes accurate rates."""
    mapq = MAPQSummary(
        total_reads=10,
        mapq_zero_count=4,
        mapq_distribution={60: 6, 0: 4},
    )
    mapq.recalculate()
    assert mapq.mapq_zero_rate == 0.4
    assert mapq.mean_mapq == 36.0


def test_alignment_validation_summary_defaults_and_validation():
    """Verify AlignmentValidationSummary can be constructed with defaults and validates state."""
    summary = AlignmentValidationSummary()
    assert summary.total_records_evaluated == 0
    assert summary.total_errors == 0
    assert summary.is_valid is True
    assert len(summary.errors) == 0
    assert isinstance(summary.mapq_summary, MAPQSummary)

    # Add an error
    err = AlignmentValidationError(
        line_number=1,
        error_kind=AlignmentErrorKind.CIGAR_LENGTH_MISMATCH,
        message="Mismatch",
    )
    summary.total_errors += 1
    summary.errors.append(err)
    assert summary.is_valid is False


def test_alignment_record_dataclass():
    """Verify AlignmentRecord creates a structured immutable SAM record."""
    rec = AlignmentRecord(
        qname="read01",
        flag=0,
        rname="chr1",
        pos=1000,
        mapq=60,
        cigar="100M",
        rnext="*",
        pnext=0,
        tlen=0,
        seq="A" * 100,
        qual="F" * 100,
        tags={"NM": 0},
        line_number=6,
    )
    assert rec.qname == "read01"
    assert rec.flag == 0
    assert rec.rname == "chr1"
    assert rec.pos == 1000
    assert rec.mapq == 60
    assert rec.cigar == "100M"
    assert rec.seq == "A" * 100
    assert rec.qual == "F" * 100
    assert rec.tags == {"NM": 0}
    assert rec.line_number == 6
