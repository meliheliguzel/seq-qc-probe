from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Any


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


class AlignmentErrorKind(Enum):
    """Error types for alignment parsing and validation."""

    CIGAR_LENGTH_MISMATCH = auto()
    INVALID_CIGAR = auto()
    MALFORMED_RECORD = auto()
    MISSING_FIELD = auto()
    INVALID_MAPQ = auto()
    INVALID_FLAG = auto()
    INVALID_POSITION = auto()
    INVALID_REFERENCE_NAME = auto()
    INVALID_HEADER = auto()
    MISMATCHED_LENGTHS = auto()
    PREMATURE_EOF = auto()
    INVALID_ALIGNMENT = auto()


@dataclass(frozen=True)
class AlignmentRecord:
    """Represents a single parsed SAM alignment record."""

    qname: str
    flag: int
    rname: str
    pos: int
    mapq: int
    cigar: str
    rnext: str
    pnext: int
    tlen: int
    seq: str
    qual: str
    tags: dict[str, Any] = field(default_factory=dict)
    line_number: int | None = None


@dataclass(frozen=True)
class CIGARValidationResult:
    """Represents the outcome of validating a CIGAR string against read sequence length."""

    is_valid: bool
    cigar: str
    expected_length: int
    actual_length: int
    operations: list[tuple[int, str]] = field(default_factory=list)
    error_kind: AlignmentErrorKind | None = None
    message: str | None = None


@dataclass(frozen=True)
class AlignmentValidationError:
    """Represents a structural validation error in an alignment record or stream."""

    line_number: int
    error_kind: AlignmentErrorKind
    message: str
    record: str | None = None
    qname: str | None = None
    cigar: str | None = None
    expected_length: int | None = None
    actual_length: int | None = None
    column: str | None = None


@dataclass
class MAPQSummary:
    """A summary of the MAPQ metrics for an alignment file."""

    total_reads: int = 0
    mapq_zero_count: int = 0
    mapq_zero_rate: float = 0.0
    mean_mapq: float = 0.0
    mapq_distribution: dict[int, int] = field(default_factory=dict)
    _total_mapq_sum: int = field(default=0, repr=False)

    @property
    def multimapping_rate(self) -> float:
        """Ratio of reads with MAPQ == 0 (multimapped reads)."""
        return self.mapq_zero_rate

    def add_read(self, mapq: int) -> None:
        """Records a single alignment read's MAPQ value and updates metrics."""
        self.total_reads += 1
        if mapq == 0:
            self.mapq_zero_count += 1
        self.mapq_distribution[mapq] = self.mapq_distribution.get(mapq, 0) + 1
        self._total_mapq_sum += mapq
        self.mapq_zero_rate = round(self.mapq_zero_count / self.total_reads, 4)
        self.mean_mapq = round(self._total_mapq_sum / self.total_reads, 2)

    def recalculate(self) -> None:
        """Recalculates mean_mapq and mapq_zero_rate based on total_reads and distribution."""
        if self.total_reads > 0:
            self.mapq_zero_rate = round(self.mapq_zero_count / self.total_reads, 4)
            self._total_mapq_sum = sum(q * count for q, count in self.mapq_distribution.items())
            self.mean_mapq = round(self._total_mapq_sum / self.total_reads, 2)
        else:
            self._total_mapq_sum = 0
            self.mapq_zero_rate = 0.0
            self.mean_mapq = 0.0


@dataclass
class AlignmentValidationSummary:
    """A summary of the alignment validation metrics for an alignment file."""

    total_records_evaluated: int = 0
    total_errors: int = 0
    error_distribution: dict[AlignmentErrorKind, int] = field(default_factory=dict)
    errors: list[AlignmentValidationError] = field(default_factory=list)
    mapq_summary: MAPQSummary = field(default_factory=MAPQSummary)

    @property
    def is_valid(self) -> bool:
        """Returns True if there are zero errors encountered."""
        return self.total_errors == 0
