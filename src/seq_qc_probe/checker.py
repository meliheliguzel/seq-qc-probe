import sys
from pathlib import Path

from .metrics import SequenceQCAccumulator
from .models import CheckResult, FASTQRecord, FASTQValidationError
from .parser import parse_and_validate_fastq


class FASTQChecker:
    """Orchestrates FASTQ validation and quality gate enforcement."""

    def __init__(
        self,
        input_path: str | Path,
        min_quality: float = 0.0,
        max_errors: int = 10,
        quiet: bool = False,
    ) -> None:
        self.input_path = str(input_path)
        self.min_quality = min_quality
        self.max_errors = max_errors
        self.quiet = quiet

    def run(self) -> CheckResult:
        """Executes validation stream, accumulates metrics, and evaluates quality gate."""
        if self.input_path != "-":
            input_file = Path(self.input_path)
            if not input_file.is_file():
                return CheckResult(
                    input_file=self.input_path,
                    min_quality=self.min_quality,
                    max_errors=self.max_errors,
                    total_records=0,
                    total_errors=0,
                    errors=[],
                    summary=SequenceQCAccumulator().summarize(),
                    status="FAILED",
                    exit_code=2,
                    error_message=f"Error: File not found: {self.input_path}",
                )

        total_records = 0
        total_errors = 0
        captured_errors: list[FASTQValidationError] = []
        accumulator = SequenceQCAccumulator()

        try:
            for item in parse_and_validate_fastq(self.input_path):
                if isinstance(item, FASTQValidationError):
                    total_errors += 1
                    if len(captured_errors) < self.max_errors:
                        captured_errors.append(item)

                    if not self.quiet and total_errors <= self.max_errors:
                        sys.stderr.write(
                            f"[{item.error_kind.name}] Line {item.line_number}: {item.message}\n"
                        )
                    if total_errors >= self.max_errors:
                        if not self.quiet:
                            sys.stderr.write(
                                f"Aborted: Maximum error threshold ({self.max_errors}) reached.\n"
                            )
                        break
                elif isinstance(item, FASTQRecord):
                    total_records += 1
                    accumulator.accumulate(item)

            summary = accumulator.summarize()

        except Exception as exc:
            return CheckResult(
                input_file=self.input_path,
                min_quality=self.min_quality,
                max_errors=self.max_errors,
                total_records=total_records,
                total_errors=total_errors,
                errors=captured_errors,
                summary=accumulator.summarize(),
                status="FAILED",
                exit_code=2,
                error_message=f"Fatal error while processing stream: {exc}",
            )

        is_failed = total_errors > 0 or summary.mean_quality_score < self.min_quality
        exit_code = 1 if is_failed else 0
        status = "FAILED" if is_failed else "PASSED"

        return CheckResult(
            input_file=self.input_path,
            min_quality=self.min_quality,
            max_errors=self.max_errors,
            total_records=total_records,
            total_errors=total_errors,
            errors=captured_errors,
            summary=summary,
            status=status,
            exit_code=exit_code,
        )
