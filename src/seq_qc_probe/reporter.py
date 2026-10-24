import json
import sys
import textwrap
from abc import ABC, abstractmethod
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, TextIO

from .models import AlignmentValidationSummary, CheckResult, VCFValidationSummary
from .version import __version__


class BaseReporter(ABC):
    """Abstract base reporter providing standard JSON serialization and I/O utilities."""

    @abstractmethod
    def to_dict(self) -> dict[str, Any]:
        """Builds a structured, JSON-serializable dictionary report."""
        ...

    @abstractmethod
    def format_text_summary(self) -> str:
        """Formats human-readable terminal output."""
        ...

    @abstractmethod
    def print_terminal_summary(
        self, stdout: TextIO | None = None, stderr: TextIO | None = None
    ) -> None:
        """Prints terminal summary and failure notifications."""
        ...

    @staticmethod
    def get_iso_timestamp() -> str:
        """Returns standardized ISO 8601 UTC timestamp."""
        return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")

    def print_json(self, stream: TextIO | None = None) -> None:
        """Prints formatted JSON report to the given stream."""
        out = sys.stdout if stream is None else stream
        out.write(json.dumps(self.to_dict(), indent=2) + "\n")

    def write_json(self, file_path: str | Path) -> None:
        """Writes formatted JSON report directly to a file."""
        target = Path(file_path)
        target.parent.mkdir(parents=True, exist_ok=True)
        with open(target, "w", encoding="utf-8") as f:
            f.write(json.dumps(self.to_dict(), indent=2) + "\n")


class QCReporter(BaseReporter):
    """Handles serialization and formatting of FASTQ quality control results."""

    def __init__(self, result: CheckResult) -> None:
        self.result = result

    def to_dict(self) -> dict[str, Any]:
        """Builds a comprehensive, hierarchical JSON-serializable report dictionary."""
        res = self.result
        summary = res.summary

        # Calculate normalized rates (0.0 to 1.0)
        q20_rate = (
            round(summary.q20_percent / 100.0, 4) if summary.total_bases_evaluated > 0 else 0.0
        )
        q30_rate = (
            round(summary.q30_percent / 100.0, 4) if summary.total_bases_evaluated > 0 else 0.0
        )
        gc_rate = (
            round(summary.gc_content_percent / 100.0, 4)
            if summary.total_bases_evaluated > 0
            else 0.0
        )
        n_rate = (
            round(summary.n_content_percent / 100.0, 6)
            if summary.total_bases_evaluated > 0
            else 0.0
        )
        contamination_rate = (
            round(summary.adapter_contamination_rate / 100.0, 4)
            if summary.total_reads_evaluated > 0
            else 0.0
        )

        timestamp = self.get_iso_timestamp()

        return {
            "schema_version": "1.0.0",
            "meta": {
                "tool": "seq-qc-probe",
                "version": __version__,
                "command": "check",
                "timestamp": timestamp,
                "input_file": res.input_file,
                "thresholds": {
                    "min_quality": res.min_quality,
                    "max_errors": res.max_errors,
                },
            },
            "quality_gate": {
                "status": res.status,
                "exit_code": res.exit_code,
                "passed": res.exit_code == 0,
                "structural_valid": res.total_errors == 0,
                "quality_threshold_met": summary.mean_quality_score >= res.min_quality,
                "total_errors": res.total_errors,
                "errors": [
                    {
                        "line_number": err.line_number,
                        "error_kind": err.error_kind.name,
                        "message": err.message,
                    }
                    for err in res.errors
                ],
            },
            "summary": {
                "total_reads": summary.total_reads_evaluated,
                "total_bases": summary.total_bases_evaluated,
                "mean_quality_score": round(summary.mean_quality_score, 2),
                "gc_bases": summary.gc_bases,
                "gc_rate": gc_rate,
                "gc_percent": round(summary.gc_content_percent, 2),
                "n_bases": summary.n_bases,
                "n_rate": n_rate,
                "n_percent": round(summary.n_content_percent, 4),
                "q20_bases": summary.q20_bases,
                "q20_rate": q20_rate,
                "q20_percent": round(summary.q20_percent, 2),
                "q30_bases": summary.q30_bases,
                "q30_rate": q30_rate,
                "q30_percent": round(summary.q30_percent, 2),
            },
            "contamination": {
                "adapter_name": "Illumina Universal Adapter",
                "adapter_sequence": "AGATCGGAAGAGC",
                "contaminated_reads": summary.adapter_contaminated_reads,
                "contamination_rate": contamination_rate,
                "contamination_percent": round(summary.adapter_contamination_rate, 2),
            },
            "distributions": {
                "per_cycle_mean_quality": [round(q, 2) for q in summary.positions_average_quality],
                "read_quality_histogram": summary.read_quality_distribution,
            },
            # Top-level backwards compatibility fields for v0.1.0 consumers
            "status": res.status,
            "total_records": res.total_records,
            "total_errors": res.total_errors,
            "total_reads_evaluated": summary.total_reads_evaluated,
            "total_bases_evaluated": summary.total_bases_evaluated,
            "mean_quality_score": round(summary.mean_quality_score, 2),
            "q20_percent": round(summary.q20_percent, 2),
            "q30_percent": round(summary.q30_percent, 2),
            "gc_content_percent": round(summary.gc_content_percent, 2),
            "n_content_percent": round(summary.n_content_percent, 2),
            "adapter_contaminated_reads": summary.adapter_contaminated_reads,
            "adapter_contamination_rate": round(summary.adapter_contamination_rate, 2),
            "positions_average_quality": [round(q, 2) for q in summary.positions_average_quality],
            "read_quality_distribution": summary.read_quality_distribution,
        }

    def format_text_summary(self) -> str:
        """Formats standard human-readable terminal output for FASTQ QC."""
        res = self.result
        summary = res.summary

        if res.status == "PASSED":
            return textwrap.dedent(f"""\
                Validation PASSED: {res.total_records} valid records processed with zero errors.

                Total Reads Evaluated: {summary.total_reads_evaluated}
                Total Bases Evaluated: {summary.total_bases_evaluated}
                Mean Quality Score: {summary.mean_quality_score:.2f}
                Q20 Percent: {summary.q20_percent:.2f}%
                Q30 Percent: {summary.q30_percent:.2f}%
                GC Content Percent: {summary.gc_content_percent:.2f}%
                N Content Percent: {summary.n_content_percent:.2f}%
                Adapter Contaminated Reads: {summary.adapter_contaminated_reads}
                Adapter Contamination Rate: {summary.adapter_contamination_rate:.2f}%
            """)
        return ""

    def print_terminal_summary(
        self, stdout: TextIO | None = None, stderr: TextIO | None = None
    ) -> None:
        """Prints terminal summary and failure notifications for FASTQ QC."""
        out = sys.stdout if stdout is None else stdout
        err = sys.stderr if stderr is None else stderr

        res = self.result
        summary = res.summary

        if res.status == "PASSED":
            out.write(self.format_text_summary())
        else:
            if res.total_errors > 0:
                err.write(
                    f"\nValidation FAILED: Found {res.total_errors} error(s) "
                    f"across {res.total_records} valid records.\n"
                )
            if summary.mean_quality_score < res.min_quality:
                err.write(
                    f"\nValidation FAILED: Mean quality score {summary.mean_quality_score:.2f} "
                    f"is below the minimum threshold of {res.min_quality}.\n"
                )


class VCFReporter(BaseReporter):
    """Handles serialization and formatting of VCF validation results."""

    def __init__(
        self,
        summary: VCFValidationSummary,
        input_file: str | Path = "",
        max_errors: int | None = None,
    ) -> None:
        self.summary = summary
        self.input_file = str(input_file)
        self.max_errors = max_errors

    def to_dict(self) -> dict[str, Any]:
        """Builds a structured JSON-serializable dictionary report for VCF validation."""
        s = self.summary
        error_dist = {kind.name: count for kind, count in s.error_distribution.items() if count > 0}
        return {
            "schema_version": "1.0.0",
            "meta": {
                "tool": "seq-qc-probe",
                "version": __version__,
                "command": "vcf-check",
                "timestamp": self.get_iso_timestamp(),
                "input_file": self.input_file,
                "thresholds": {
                    "max_errors": self.max_errors,
                },
            },
            "validation_gate": {
                "status": "PASSED" if s.is_valid else "FAILED",
                "passed": s.is_valid,
                "exit_code": 0 if s.is_valid else 1,
                "total_errors": s.total_errors,
                "error_distribution": error_dist,
                "errors": [
                    {
                        "line_number": err.line_number,
                        "error_kind": err.error_kind.name,
                        "message": err.message,
                        "column": err.column,
                    }
                    for err in s.errors
                ],
            },
            "summary": {
                "total_records_evaluated": s.total_records_evaluated,
                "total_errors": s.total_errors,
                "is_valid": s.is_valid,
            },
        }

    def format_text_summary(self) -> str:
        """Formats standard human-readable terminal output for VCF validation."""
        s = self.summary
        if s.is_valid:
            return textwrap.dedent(f"""\
                Validation PASSED: VCF 4.2 file is valid.
                Total Records Evaluated: {s.total_records_evaluated}
                Total Errors: 0
            """)

        lines = [
            f"Validation FAILED: Found {s.total_errors} error(s) "
            f"across {s.total_records_evaluated} record(s).",
            "",
            "Error Breakdown:",
        ]
        for kind, count in s.error_distribution.items():
            if count > 0:
                lines.append(f"  - {kind.name}: {count}")

        if s.errors:
            lines.append("")
            lines.append("Errors:")
            for err in s.errors:
                col_info = f" (Column: {err.column})" if err.column else ""
                lines.append(
                    f"  [Line {err.line_number}] [{err.error_kind.name}] {err.message}{col_info}"
                )
        return "\n".join(lines) + "\n"

    def print_terminal_summary(
        self, stdout: TextIO | None = None, stderr: TextIO | None = None
    ) -> None:
        """Prints terminal summary to stdout or stderr depending on pass/fail."""
        out = sys.stdout if stdout is None else stdout
        err = sys.stderr if stderr is None else stderr

        if self.summary.is_valid:
            out.write(self.format_text_summary())
        else:
            err.write(self.format_text_summary())


class AlignmentReporter(BaseReporter):
    """Handles serialization and formatting of SAM/BAM alignment validation and MAPQ results."""

    def __init__(
        self,
        summary: AlignmentValidationSummary,
        input_file: str | Path = "",
        max_errors: int | None = None,
    ) -> None:
        self.summary = summary
        self.input_file = str(input_file)
        self.max_errors = max_errors

    def to_dict(self) -> dict[str, Any]:
        """Builds a structured JSON-serializable dictionary report for SAM validation."""
        s = self.summary
        error_dist = {kind.name: count for kind, count in s.error_distribution.items() if count > 0}
        mapq = s.mapq_summary

        return {
            "schema_version": "1.0.0",
            "meta": {
                "tool": "seq-qc-probe",
                "version": __version__,
                "command": "sam-check",
                "timestamp": self.get_iso_timestamp(),
                "input_file": self.input_file,
                "thresholds": {
                    "max_errors": self.max_errors,
                },
            },
            "validation_gate": {
                "status": "PASSED" if s.is_valid else "FAILED",
                "passed": s.is_valid,
                "exit_code": 0 if s.is_valid else 1,
                "total_errors": s.total_errors,
                "error_distribution": error_dist,
                "errors": [
                    {
                        "line_number": err.line_number,
                        "error_kind": err.error_kind.name,
                        "message": err.message,
                        "qname": err.qname,
                        "cigar": err.cigar,
                        "expected_length": err.expected_length,
                        "actual_length": err.actual_length,
                        "column": err.column,
                    }
                    for err in s.errors
                ],
            },
            "mapq_metrics": {
                "total_reads": mapq.total_reads,
                "mapq_zero_count": mapq.mapq_zero_count,
                "mapq_zero_rate": mapq.mapq_zero_rate,
                "multimapping_rate": mapq.multimapping_rate,
                "mean_mapq": mapq.mean_mapq,
                "mapq_distribution": mapq.mapq_distribution,
            },
            "summary": {
                "total_records_evaluated": s.total_records_evaluated,
                "total_errors": s.total_errors,
                "is_valid": s.is_valid,
            },
        }

    def format_text_summary(self) -> str:
        """Formats standard human-readable terminal output for SAM validation."""
        s = self.summary
        mapq = s.mapq_summary
        pct_zero = mapq.mapq_zero_rate * 100.0

        mapq_block = textwrap.dedent(f"""\
            MAPQ & Multimapping Metrics:
              Total Reads: {mapq.total_reads}
              Mean MAPQ: {mapq.mean_mapq:.2f}
              MAPQ = 0 (Multimapped): {mapq.mapq_zero_count} ({pct_zero:.2f}%)
        """).strip()

        if s.is_valid:
            header_block = textwrap.dedent(f"""\
                Validation PASSED: SAM alignment file is valid.
                Total Records Evaluated: {s.total_records_evaluated}
                Total Errors: 0
            """).strip()
            return f"{header_block}\n\n{mapq_block}\n"

        lines = [
            f"Validation FAILED: Found {s.total_errors} error(s) "
            f"across {s.total_records_evaluated} record(s).",
            "",
            mapq_block,
            "",
            "Error Breakdown:",
        ]
        for kind, count in s.error_distribution.items():
            if count > 0:
                lines.append(f"  - {kind.name}: {count}")

        if s.errors:
            lines.append("")
            lines.append("Errors:")
            for err in s.errors:
                context_parts: list[str] = []
                if err.qname:
                    context_parts.append(f"Read: {err.qname}")
                if err.cigar:
                    context_parts.append(f"CIGAR: {err.cigar}")
                if err.column:
                    context_parts.append(f"Column: {err.column}")
                ctx_info = f" ({', '.join(context_parts)})" if context_parts else ""
                lines.append(
                    f"  [Line {err.line_number}] [{err.error_kind.name}]{ctx_info} {err.message}"
                )

        return "\n".join(lines) + "\n"

    def print_terminal_summary(
        self, stdout: TextIO | None = None, stderr: TextIO | None = None
    ) -> None:
        """Prints terminal summary to stdout or stderr depending on pass/fail."""
        out = sys.stdout if stdout is None else stdout
        err = sys.stderr if stderr is None else stderr

        if self.summary.is_valid:
            out.write(self.format_text_summary())
        else:
            err.write(self.format_text_summary())
