from dataclasses import dataclass, field
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
    gc_bases: int = 0
    q20_bases: int = 0
    q30_bases: int = 0
    n_bases: int = 0


@dataclass(frozen=True)
class CheckResult:
    """Encapsulates the complete outcome of a sequence validation run."""

    input_file: str
    min_quality: float
    max_errors: int
    total_records: int
    total_errors: int
    errors: list[FASTQValidationError]
    summary: SequenceQCMetricsSummary
    status: str
    exit_code: int
    error_message: str | None = None


class VCFErrorKind(Enum):
    """Error types for VCF parsing and validation."""

    INVALID_RECORD = auto()
    MISSING_FIELD = auto()
    INVALID_GENOTYPE = auto()
    PREMATURE_EOF = auto()
    MISSING_FILEFORMAT = auto()
    INVALID_HEADER_COLUMNS = auto()
    INVALID_COORDINATE = auto()
    INVALID_BASE = auto()
    MALFORMED_LINE = auto()


@dataclass(frozen=True)
class VCFValidationError:
    """Represents a structural validation error in a VCF record or stream."""

    line_number: int
    error_kind: VCFErrorKind
    message: str
    column: str | None = None


@dataclass
class VCFValidationSummary:
    """A summary of the VCF validation metrics for a VCF file."""

    total_records_evaluated: int
    total_errors: int
    error_distribution: dict[VCFErrorKind, int]
    errors: list[VCFValidationError] = field(default_factory=list)

    @property
    def is_valid(self) -> bool:
        """Returns True if there are zero errors encountered."""
        return self.total_errors == 0
