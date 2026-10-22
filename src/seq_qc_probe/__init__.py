"""seq-qc-probe: High-throughput FASTQ quality control and structural integrity validator."""

from .metrics import SequenceQCAccumulator, analyze_fastq_file
from .models import FASTQErrorKind, FASTQRecord, FASTQValidationError, SequenceQCMetricsSummary
from .parser import parse_and_validate_fastq

__version__ = "0.1.0"

__all__ = [
    "FASTQErrorKind",
    "FASTQRecord",
    "FASTQValidationError",
    "SequenceQCAccumulator",
    "SequenceQCMetricsSummary",
    "__version__",
    "analyze_fastq_file",
    "parse_and_validate_fastq",
]
