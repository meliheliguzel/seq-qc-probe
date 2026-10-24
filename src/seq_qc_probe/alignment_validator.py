"""SAM alignment record, CIGAR string, and MAPQ metrics validator."""

import re
from pathlib import Path
from typing import Any

from .io import open_stream
from .models import (
    AlignmentErrorKind,
    AlignmentRecord,
    AlignmentValidationError,
    AlignmentValidationSummary,
    CIGARValidationResult,
    MAPQSummary,
)

# CIGAR specification according to SAM v1.6
CIGAR_SYNTAX_PATTERN = re.compile(r"^(\d+[MIDNSHP=X])+$")
CIGAR_OPERATION_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")

QUERY_CONSUMING_OPS: frozenset[str] = frozenset("MIS=X")

REFERENCE_CONSUMING_OPS: frozenset[str] = frozenset("MDN=X")

MANDATORY_COLUMN_COUNT = 11


class CIGARValidator:
    """Validates and evaluates CIGAR strings against sequence lengths."""

    @staticmethod
    def validate(cigar: str, seq_length: int | None = None) -> CIGARValidationResult:
        """Validates CIGAR syntax and mathematical equivalence with sequence length."""
        if cigar == "*":
            return CIGARValidationResult(
                is_valid=True,
                cigar="*",
                expected_length=seq_length or 0,
                actual_length=seq_length or 0,
                operations=[],
            )

        if not CIGAR_SYNTAX_PATTERN.match(cigar):
            return CIGARValidationResult(
                is_valid=False,
                cigar=cigar,
                expected_length=0,
                actual_length=seq_length or 0,
                operations=[],
                error_kind=AlignmentErrorKind.INVALID_CIGAR,
                message=(
                    f"Invalid CIGAR syntax: '{cigar}' does not match pattern '^(\\d+[MIDNSHP=X])+$'"
                ),
            )

        raw_ops = CIGAR_OPERATION_PATTERN.findall(cigar)
        operations = [(int(length), op) for length, op in raw_ops]

        expected_length = sum(length for length, op in operations if op in QUERY_CONSUMING_OPS)

        if seq_length is not None and expected_length != seq_length:
            msg = (
                f"CIGAR length mismatch: operations '{cigar}' consume "
                f"{expected_length} bp of read sequence, but SEQ length is {seq_length} bp"
            )
            return CIGARValidationResult(
                is_valid=False,
                cigar=cigar,
                expected_length=expected_length,
                actual_length=seq_length,
                operations=operations,
                error_kind=AlignmentErrorKind.CIGAR_LENGTH_MISMATCH,
                message=msg,
            )

        return CIGARValidationResult(
            is_valid=True,
            cigar=cigar,
            expected_length=expected_length,
            actual_length=seq_length if seq_length is not None else expected_length,
            operations=operations,
        )


def validate_cigar(cigar: str, seq_length: int | None = None) -> CIGARValidationResult:
    """Convenience functional interface for CIGAR validation."""
    return CIGARValidator.validate(cigar, seq_length)


class SAMRecordParser:
    """Parses and validates individual SAM alignment records adhering to Single Responsibility."""

    def __init__(self, cigar_validator: CIGARValidator | None = None) -> None:
        self.cigar_validator = cigar_validator or CIGARValidator()

    def parse_and_validate(
        self, line: str, line_number: int
    ) -> tuple[AlignmentRecord | None, list[AlignmentValidationError], int | None]:
        """Parses a SAM line and validates all 11 mandatory fields.

        Returns:
            A tuple of (parsed_record_or_None, list_of_validation_errors, valid_mapq_or_None).
        """
        errors: list[AlignmentValidationError] = []
        columns = line.split("\t")

        if len(columns) < MANDATORY_COLUMN_COUNT:
            errors.append(
                AlignmentValidationError(
                    line_number=line_number,
                    error_kind=AlignmentErrorKind.MALFORMED_RECORD,
                    message=(
                        f"Record has {len(columns)} column(s); at least "
                        f"{MANDATORY_COLUMN_COUNT} mandatory tab-delimited columns required"
                    ),
                    record=line,
                )
            )
            return None, errors, None

        qname = columns[0]
        flag_str = columns[1]
        rname = columns[2]
        pos_str = columns[3]
        mapq_str = columns[4]
        cigar = columns[5]
        rnext = columns[6]
        pnext_str = columns[7]
        tlen_str = columns[8]
        seq = columns[9]
        qual = columns[10]

        if not qname:
            errors.append(
                AlignmentValidationError(
                    line_number=line_number,
                    error_kind=AlignmentErrorKind.MISSING_FIELD,
                    message="QNAME cannot be empty",
                    record=line,
                    column="QNAME",
                )
            )

        if not rname:
            errors.append(
                AlignmentValidationError(
                    line_number=line_number,
                    error_kind=AlignmentErrorKind.MISSING_FIELD,
                    message="RNAME cannot be empty",
                    record=line,
                    column="RNAME",
                )
            )

        flag = self._validate_flag(flag_str, line_number, qname, line, errors)

        pos = self._validate_pos(pos_str, line_number, qname, line, errors)

        mapq = self._validate_mapq(mapq_str, line_number, qname, line, errors)

        pnext = self._validate_integer_field(
            pnext_str, "PNEXT", line_number, qname, line, errors, min_value=0
        )

        tlen = self._validate_integer_field(tlen_str, "TLEN", line_number, qname, line, errors)

        seq_length = None if seq == "*" else len(seq)
        self._validate_seq_qual(seq, qual, line_number, qname, cigar, line, errors)

        self._validate_cigar_field(cigar, seq_length, line_number, qname, line, errors)

        # Parse optional tags if any
        tags: dict[str, Any] = {}
        if len(columns) > MANDATORY_COLUMN_COUNT:
            for tag_entry in columns[MANDATORY_COLUMN_COUNT:]:
                parts = tag_entry.split(":", 2)
                if len(parts) == 3:
                    tags[parts[0]] = parts[2]

        record: AlignmentRecord | None = None
        if (
            flag is not None
            and pos is not None
            and mapq is not None
            and pnext is not None
            and tlen is not None
        ):
            record = AlignmentRecord(
                qname=qname,
                flag=flag,
                rname=rname,
                pos=pos,
                mapq=mapq,
                cigar=cigar,
                rnext=rnext,
                pnext=pnext,
                tlen=tlen,
                seq=seq,
                qual=qual,
                tags=tags,
                line_number=line_number,
            )

        return record, errors, mapq

    def _validate_flag(
        self,
        flag_str: str,
        line_number: int,
        qname: str,
        line: str,
        errors: list[AlignmentValidationError],
    ) -> int | None:
        try:
            flag = int(flag_str)
            if flag < 0:
                raise ValueError("FLAG cannot be negative")
            return flag
        except ValueError:
            errors.append(
                AlignmentValidationError(
                    line_number=line_number,
                    error_kind=AlignmentErrorKind.INVALID_FLAG,
                    message=f"FLAG must be a non-negative integer: '{flag_str}'",
                    record=line,
                    qname=qname,
                    column="FLAG",
                )
            )
            return None

    def _validate_pos(
        self,
        pos_str: str,
        line_number: int,
        qname: str,
        line: str,
        errors: list[AlignmentValidationError],
    ) -> int | None:
        try:
            pos = int(pos_str)
            if pos < 0:
                raise ValueError("POS cannot be negative")
            return pos
        except ValueError:
            errors.append(
                AlignmentValidationError(
                    line_number=line_number,
                    error_kind=AlignmentErrorKind.INVALID_POSITION,
                    message=f"POS coordinate must be a non-negative integer: '{pos_str}'",
                    record=line,
                    qname=qname,
                    column="POS",
                )
            )
            return None

    def _validate_mapq(
        self,
        mapq_str: str,
        line_number: int,
        qname: str,
        line: str,
        errors: list[AlignmentValidationError],
    ) -> int | None:
        try:
            mapq = int(mapq_str)
            if mapq < 0 or mapq > 255:
                errors.append(
                    AlignmentValidationError(
                        line_number=line_number,
                        error_kind=AlignmentErrorKind.INVALID_MAPQ,
                        message=f"MAPQ score out of range [0, 255]: '{mapq_str}'",
                        record=line,
                        qname=qname,
                        column="MAPQ",
                    )
                )
                return None
            return mapq
        except ValueError:
            errors.append(
                AlignmentValidationError(
                    line_number=line_number,
                    error_kind=AlignmentErrorKind.INVALID_MAPQ,
                    message=f"MAPQ score must be an integer between 0 and 255: '{mapq_str}'",
                    record=line,
                    qname=qname,
                    column="MAPQ",
                )
            )
            return None

    def _validate_integer_field(
        self,
        val_str: str,
        field_name: str,
        line_number: int,
        qname: str,
        line: str,
        errors: list[AlignmentValidationError],
        min_value: int | None = None,
    ) -> int | None:
        try:
            val = int(val_str)
            if min_value is not None and val < min_value:
                raise ValueError(f"{field_name} must be >= {min_value}")
            return val
        except ValueError:
            errors.append(
                AlignmentValidationError(
                    line_number=line_number,
                    error_kind=AlignmentErrorKind.MALFORMED_RECORD,
                    message=f"{field_name} must be an integer: '{val_str}'",
                    record=line,
                    qname=qname,
                    column=field_name,
                )
            )
            return None

    def _validate_seq_qual(
        self,
        seq: str,
        qual: str,
        line_number: int,
        qname: str,
        cigar: str,
        line: str,
        errors: list[AlignmentValidationError],
    ) -> None:
        if seq != "*" and qual != "*" and len(seq) != len(qual):
            errors.append(
                AlignmentValidationError(
                    line_number=line_number,
                    error_kind=AlignmentErrorKind.MISMATCHED_LENGTHS,
                    message=(
                        f"Length of QUAL ({len(qual)}) does not match length of SEQ ({len(seq)})"
                    ),
                    record=line,
                    qname=qname,
                    cigar=cigar,
                    expected_length=len(seq),
                    actual_length=len(qual),
                    column="QUAL",
                )
            )

    def _validate_cigar_field(
        self,
        cigar: str,
        seq_length: int | None,
        line_number: int,
        qname: str,
        line: str,
        errors: list[AlignmentValidationError],
    ) -> None:
        cigar_res = self.cigar_validator.validate(cigar, seq_length)
        if not cigar_res.is_valid and cigar_res.error_kind:
            errors.append(
                AlignmentValidationError(
                    line_number=line_number,
                    error_kind=cigar_res.error_kind,
                    message=cigar_res.message or "Invalid CIGAR",
                    record=line,
                    qname=qname,
                    cigar=cigar,
                    expected_length=cigar_res.expected_length,
                    actual_length=cigar_res.actual_length,
                    column="CIGAR",
                )
            )


class AlignmentValidator:
    """Validates SAM format files and streams for structural, CIGAR, and MAPQ integrity."""

    def __init__(
        self,
        file_path: str | Path,
        max_errors: int | None = None,
        parser: SAMRecordParser | None = None,
    ) -> None:
        self.file_path = file_path
        self.max_errors = max_errors
        self.parser = parser or SAMRecordParser()
        self.seen_records: bool = False
        self.summary = AlignmentValidationSummary(
            total_records_evaluated=0,
            total_errors=0,
            error_distribution={kind: 0 for kind in AlignmentErrorKind},
            errors=[],
            mapq_summary=MAPQSummary(),
        )

    def validate(self) -> AlignmentValidationSummary:
        """Validates the SAM file or stream with a memory-efficient single-pass pipeline."""
        stream, should_close = open_stream(self.file_path)
        has_content = False

        try:
            for line_number, raw_line in enumerate(stream, start=1):
                line = raw_line.rstrip("\r\n")

                if not line.strip():
                    continue

                has_content = True

                if self.max_errors and self.summary.total_errors >= self.max_errors:
                    break

                self._process_line(line, line_number)

            if not has_content:
                self._record_error(
                    line_number=1,
                    error_kind=AlignmentErrorKind.PREMATURE_EOF,
                    message="SAM file is empty",
                )
        finally:
            if should_close:
                stream.close()

        return self.summary

    def validate_record(self, record: str, line_number: int = 1) -> None:
        """Validates a single alignment record directly and updates the summary."""
        self._process_record_line(record, line_number)

    def _process_line(self, line: str, line_number: int) -> None:
        """Processes a single line, delegating to header or record handler."""
        if line.startswith("@"):
            self._process_header_line(line, line_number)
        else:
            self._process_record_line(line, line_number)

    def _process_header_line(self, line: str, line_number: int) -> None:
        """Validates SAM header lines, ensuring headers precede all alignment records."""
        if self.seen_records:
            self._record_error(
                line_number=line_number,
                error_kind=AlignmentErrorKind.INVALID_HEADER,
                message="SAM header line '@' found after alignment records",
                record=line,
                column="HEADER",
            )

    def _process_record_line(self, line: str, line_number: int) -> None:
        """Validates a single alignment record and records errors and MAPQ metrics."""
        self.seen_records = True
        self.summary.total_records_evaluated += 1

        _, errors, mapq = self.parser.parse_and_validate(line, line_number)

        for err in errors:
            self._record_error(
                line_number=err.line_number,
                error_kind=err.error_kind,
                message=err.message,
                record=err.record,
                qname=err.qname,
                cigar=err.cigar,
                expected_length=err.expected_length,
                actual_length=err.actual_length,
                column=err.column,
            )

        if mapq is not None and 0 <= mapq <= 255:
            self.summary.mapq_summary.add_read(mapq)

    def _record_error(
        self,
        line_number: int,
        error_kind: AlignmentErrorKind,
        message: str,
        record: str | None = None,
        qname: str | None = None,
        cigar: str | None = None,
        expected_length: int | None = None,
        actual_length: int | None = None,
        column: str | None = None,
    ) -> None:
        """Records a validation error in the summary."""
        self.summary.total_errors += 1
        self.summary.error_distribution[error_kind] += 1
        self.summary.errors.append(
            AlignmentValidationError(
                line_number=line_number,
                error_kind=error_kind,
                message=message,
                record=record,
                qname=qname,
                cigar=cigar,
                expected_length=expected_length,
                actual_length=actual_length,
                column=column,
            )
        )
