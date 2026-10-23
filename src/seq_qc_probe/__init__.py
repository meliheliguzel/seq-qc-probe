"""seq-qc-probe: High-throughput FASTQ quality control and structural integrity validator."""

from .checker import FASTQChecker
from .metrics import SequenceQCAccumulator, analyze_fastq_file
from .models import (
    CheckResult,
    FASTQErrorKind,
    FASTQRecord,
    FASTQValidationError,
    SequenceQCMetricsSummary,
    VCFErrorKind,
    VCFValidationError,
    VCFValidationSummary,
)
from .parser import parse_and_validate_fastq
from .reporter import QCReporter
from .vcf_validator import VCFValidator
from .version import __version__

__all__ = [
    "CheckResult",
    "FASTQChecker",
    "FASTQErrorKind",
    "FASTQRecord",
    "FASTQValidationError",
    "QCReporter",
    "SequenceQCAccumulator",
    "SequenceQCMetricsSummary",
    "VCFErrorKind",
    "VCFValidationError",
    "VCFValidationSummary",
    "VCFValidator",
    "__version__",
    "analyze_fastq_file",
    "parse_and_validate_fastq",
]
