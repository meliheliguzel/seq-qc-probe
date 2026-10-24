"""seq-qc-probe: High-throughput FASTQ quality control and structural integrity validator."""

from .alignment_validator import AlignmentValidator, CIGARValidator, validate_cigar
from .checker import FASTQChecker
from .metrics import SequenceQCAccumulator, analyze_fastq_file
from .models import (
    AlignmentErrorKind,
    AlignmentRecord,
    AlignmentValidationError,
    AlignmentValidationSummary,
    CheckResult,
    CIGARValidationResult,
    FASTQErrorKind,
    FASTQRecord,
    FASTQValidationError,
    MAPQSummary,
    SequenceQCMetricsSummary,
    VCFErrorKind,
    VCFValidationError,
    VCFValidationSummary,
)
from .parser import parse_and_validate_fastq
from .reporter import AlignmentReporter, BaseReporter, QCReporter, VCFReporter
from .vcf_validator import VCFValidator
from .version import __version__

__all__ = [
    "AlignmentErrorKind",
    "AlignmentRecord",
    "AlignmentReporter",
    "AlignmentValidationError",
    "AlignmentValidationSummary",
    "AlignmentValidator",
    "BaseReporter",
    "CIGARValidationResult",
    "CIGARValidator",
    "CheckResult",
    "FASTQChecker",
    "FASTQErrorKind",
    "FASTQRecord",
    "FASTQValidationError",
    "MAPQSummary",
    "QCReporter",
    "SequenceQCAccumulator",
    "SequenceQCMetricsSummary",
    "VCFErrorKind",
    "VCFReporter",
    "VCFValidationError",
    "VCFValidationSummary",
    "VCFValidator",
    "__version__",
    "analyze_fastq_file",
    "parse_and_validate_fastq",
    "validate_cigar",
]
