"""VCF 4.2 format, header, and column validator."""

from pathlib import Path

from .io import open_stream
from .models import VCFErrorKind, VCFValidationError, VCFValidationSummary

MANDATORY_COLUMNS: list[str] = [
    "#CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO",
]
OPTIONAL_FORMAT_COLUMN: str = "FORMAT"
VALID_BASES: frozenset[str] = frozenset("ACGTNacgtn")


class VCFValidator:
    """Validates VCF 4.2 files for structural correctness and biological consistency."""

    def __init__(
        self,
        file_path: str | Path,
        max_errors: int | None = None,
    ) -> None:
        self.file_path = file_path
        self.max_errors = max_errors
        self.header_columns: list[str] | None = None
        self.header_seen: bool = False
        self.summary = VCFValidationSummary(
            total_records_evaluated=0,
            total_errors=0,
            error_distribution={kind: 0 for kind in VCFErrorKind},
            errors=[],
        )

    def validate(self) -> VCFValidationSummary:
        """Validates the VCF file with a streamlined, single-pass pipeline."""
        stream, should_close = open_stream(self.file_path)
        total_lines = 0

        try:
            for line_number, raw_line in enumerate(stream, start=1):
                total_lines = line_number
                line = raw_line.rstrip("\r\n")

                if self.max_errors and self.summary.total_errors >= self.max_errors:
                    break

                self._process_line(line, line_number)

            if total_lines == 0:
                self._record_error(
                    line_number=1,
                    error_kind=VCFErrorKind.PREMATURE_EOF,
                    message="VCF file is empty",
                )
            elif not self.header_seen and (
                not self.max_errors or self.summary.total_errors < self.max_errors
            ):
                self._record_error(
                    line_number=total_lines,
                    error_kind=VCFErrorKind.INVALID_HEADER_COLUMNS,
                    message="Missing mandatory #CHROM header line in VCF file",
                )
        finally:
            if should_close:
                stream.close()

        return self.summary

    def _process_line(self, line: str, line_number: int) -> None:
        """Routes each line to its corresponding domain validator based on file structure."""
        if line_number == 1:
            self.fileformat_line_control(line, line_number)
            if not line.startswith("##"):
                self._route_non_meta_line(line, line_number)
            return

        if line.startswith("##"):
            self.other_meta_line_control(line, line_number)
        elif line.startswith("#"):
            self.header_line_control(line, line_number)
        else:
            self.record_line_control(line, line_number)

    def _route_non_meta_line(self, line: str, line_number: int) -> None:
        """Handles non-meta lines appearing unexpectedly at the very first line."""
        if line.startswith("#"):
            self.header_line_control(line, line_number)
        elif line.strip():
            self.record_line_control(line, line_number)

    def fileformat_line_control(self, line: str, line_number: int) -> bool:
        """Validates the mandatory ##fileformat=VCFv4.2 specification on the initial line."""
        if line.strip() != "##fileformat=VCFv4.2":
            msg = (
                f"Missing or incorrect fileformat line: "
                f"expected '##fileformat=VCFv4.2', found '{line[:40]}'"
            )
            self._record_error(
                line_number=line_number,
                error_kind=VCFErrorKind.MISSING_FILEFORMAT,
                message=msg,
            )
            return False
        return True

    def other_meta_line_control(self, line: str, line_number: int) -> bool:
        """Validates general meta-information lines (e.g. ##FILTER, ##INFO, ##contig)."""
        if self.header_seen:
            self._record_error(
                line_number=line_number,
                error_kind=VCFErrorKind.MALFORMED_LINE,
                message="Meta-information line '##' found after #CHROM header",
            )
            return False

        content = line[2:]
        if "=" not in content:
            self._record_error(
                line_number=line_number,
                error_kind=VCFErrorKind.MALFORMED_LINE,
                message=f"Malformed meta line missing '=' key-value delimiter: '{line[:40]}'",
            )
            return False

        key, value = content.split("=", 1)
        if not key.strip():
            self._record_error(
                line_number=line_number,
                error_kind=VCFErrorKind.MALFORMED_LINE,
                message=f"Meta line has empty key: '{line[:40]}'",
            )
            return False

        if value.startswith("<") and not value.endswith(">"):
            self._record_error(
                line_number=line_number,
                error_kind=VCFErrorKind.MALFORMED_LINE,
                message=f"Unclosed structured meta line attribute value: '{line[:40]}'",
            )
            return False

        return True

    def header_line_control(self, line: str, line_number: int) -> bool:
        """Validates the #CHROM header line and sample column structure."""
        if self.header_seen:
            self._record_error(
                line_number=line_number,
                error_kind=VCFErrorKind.MALFORMED_LINE,
                message="Duplicate header line encountered",
            )
            return False

        line_parts = line.split("\t")
        self.header_columns = line_parts
        self.header_seen = True

        if len(line_parts) < 8:
            self._record_error(
                line_number=line_number,
                error_kind=VCFErrorKind.INVALID_HEADER_COLUMNS,
                message=f"Header line has insufficient columns: {len(line_parts)}",
            )
            return False

        if line_parts[:8] != MANDATORY_COLUMNS:
            mismatched = [
                f"{exp}!={got}"
                for exp, got in zip(MANDATORY_COLUMNS, line_parts[:8], strict=False)
                if exp != got
            ]
            self._record_error(
                line_number=line_number,
                error_kind=VCFErrorKind.INVALID_HEADER_COLUMNS,
                message=f"Header columns do not match VCF 4.2 standard: {', '.join(mismatched)}",
            )
            return False

        if len(line_parts) == 9:
            self._record_error(
                line_number=line_number,
                error_kind=VCFErrorKind.INVALID_HEADER_COLUMNS,
                message="FORMAT column present without any sample genotype columns",
            )
            return False

        if len(line_parts) >= 10:
            if line_parts[8] != OPTIONAL_FORMAT_COLUMN:
                msg = (
                    f"Column 9 must be 'FORMAT' when sample columns exist, found '{line_parts[8]}'"
                )
                self._record_error(
                    line_number=line_number,
                    error_kind=VCFErrorKind.INVALID_HEADER_COLUMNS,
                    message=msg,
                    column=line_parts[8],
                )
                return False

            for sample in line_parts[9:]:
                if not sample.strip():
                    self._record_error(
                        line_number=line_number,
                        error_kind=VCFErrorKind.INVALID_HEADER_COLUMNS,
                        message="Empty sample identifier found in header",
                    )
                    return False

        return True

    def record_line_control(self, line: str, line_number: int) -> bool:
        """Validates a record line against column count, coordinates, and nucleotide alleles."""
        if not self.header_seen:
            self._record_error(
                line_number=line_number,
                error_kind=VCFErrorKind.INVALID_HEADER_COLUMNS,
                message="Data record found before #CHROM header line",
            )
            self.header_seen = True

        if not line.strip():
            self._record_error(
                line_number=line_number,
                error_kind=VCFErrorKind.MALFORMED_LINE,
                message="Empty line encountered in data section",
            )
            return False

        if line.startswith("#"):
            self._record_error(
                line_number=line_number,
                error_kind=VCFErrorKind.MALFORMED_LINE,
                message="Unexpected header or comment line in record section",
            )
            return False

        self.summary.total_records_evaluated += 1

        line_parts = line.split("\t")
        if len(line_parts) < 8:
            self._record_error(
                line_number=line_number,
                error_kind=VCFErrorKind.INVALID_RECORD,
                message=f"Record line has insufficient columns: {len(line_parts)}",
            )
            return False

        if self.header_columns is not None:
            expected = len(self.header_columns)
            if len(line_parts) != expected:
                self._record_error(
                    line_number=line_number,
                    error_kind=VCFErrorKind.INVALID_RECORD,
                    message=f"Record line has {len(line_parts)} columns, expected {expected}",
                )
                return False

        _chrom, pos, _id, ref, alt, _qual, _filter, _info = line_parts[:8]

        try:
            pos_val = int(pos)
            if pos_val <= 0:
                self._record_error(
                    line_number=line_number,
                    error_kind=VCFErrorKind.INVALID_COORDINATE,
                    message=f"Invalid POS value: {pos} (must be a 1-based positive integer)",
                    column="POS",
                )
        except ValueError:
            self._record_error(
                line_number=line_number,
                error_kind=VCFErrorKind.INVALID_COORDINATE,
                message=f"Invalid non-integer POS value: {pos}",
                column="POS",
            )

        if not ref or ref == "." or not all(base in VALID_BASES for base in ref):
            self._record_error(
                line_number=line_number,
                error_kind=VCFErrorKind.INVALID_BASE,
                message=f"Invalid REF base: {ref}",
                column="REF",
            )

        alt_alleles = alt.split(",")
        for allele in alt_alleles:
            if not allele:
                self._record_error(
                    line_number=line_number,
                    error_kind=VCFErrorKind.INVALID_BASE,
                    message="Empty ALT allele found",
                    column="ALT",
                )
            elif allele in (".", "*") or (
                allele.startswith("<") and allele.endswith(">") and len(allele) > 2
            ):
                continue
            elif not all(base in VALID_BASES for base in allele):
                self._record_error(
                    line_number=line_number,
                    error_kind=VCFErrorKind.INVALID_BASE,
                    message=f"Invalid ALT base: {allele}",
                    column="ALT",
                )

        return True

    def _record_error(
        self,
        line_number: int,
        error_kind: VCFErrorKind,
        message: str,
        column: str | None = None,
    ) -> None:
        """Records a validation error in the summary."""
        self.summary.total_errors += 1
        self.summary.error_distribution[error_kind] += 1
        self.summary.errors.append(
            VCFValidationError(
                line_number=line_number,
                error_kind=error_kind,
                message=message,
                column=column,
            )
        )
