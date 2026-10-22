from dataclasses import dataclass
from enum import Enum, auto


class FASTQErrorKind(Enum):
    """Error types for FASTQ parsing and validation."""

    INVALID_HEADER = auto()
    INVALID_SEPARATOR = auto()
    INVALID_SEQUENCE = auto()
    INVALID_QUALITY = auto()
    MISMATCHED_LENGTHS = auto()
    PREMATURE_EOF = auto()
    EMPTY_LINE = auto()


@dataclass(frozen=True)
class FASTQRecord:
    """Represents a single valid FASTQ record."""

    header: str
    sequence: str
    separator: str
    quality: str
    line_number: int


@dataclass(frozen=True)
class FASTQValidationError:
    """Represents a structural validation error in a FASTQ record or stream."""

    line_number: int
    error_kind: FASTQErrorKind
    message: str
    record: FASTQRecord | None = None


@dataclass
class SequenceQCMetricsSummary:
    """A summary of the Sequence QC metrics for a FASTQ file."""

    total_reads_evaluated: int
    mean_quality_score: float
    positions_average_quality: list[float]
    read_quality_distribution: dict[int, int]
    gc_content_percent: float
    adapter_contaminated_reads: int
    adapter_contamination_rate: float
    total_bases_evaluated: int = 0
    q20_percent: float = 0.0
    q30_percent: float = 0.0
    n_content_percent: float = 0.0
