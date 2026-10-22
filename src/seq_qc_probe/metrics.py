from .models import FASTQRecord, FASTQValidationError, SequenceQCMetricsSummary
from .parser import parse_and_validate_fastq

ILLUMINA_UNIVERSAL_ADAPTER: str = "AGATCGGAAGAGC"


def is_adapter_contaminated(
    seq: str,
    adapter: str = ILLUMINA_UNIVERSAL_ADAPTER,
    min_overlap: int = 8,
    max_mismatches: int = 1,
) -> bool:
    """Detects Illumina adapter contamination.

    Supports partial 3' overlap or full match with mismatch tolerance.
    """
    seq_upper = seq.upper()
    adapter_len = len(adapter)
    seq_len = len(seq_upper)

    if seq_len >= adapter_len:
        for i in range(seq_len - adapter_len + 1):
            sub = seq_upper[i : i + adapter_len]
            mismatches = sum(1 for a, b in zip(sub, adapter, strict=True) if a != b)
            if mismatches <= max_mismatches:
                return True

    max_prefix = min(adapter_len - 1, seq_len)
    return any(seq_upper.endswith(adapter[:k]) for k in range(max_prefix, min_overlap - 1, -1))


class SequenceQCAccumulator:
    """Accumulates Phred quality and biological metrics from FASTQ records in streaming mode."""

    total_reads: int
    total_bases: int
    total_quality_sum: int
    read_quality_distribution: dict[int, int]
    position_quality_sums: list[int]
    position_counts: list[int]
    total_gc_bases: int
    q20_bases: int
    q30_bases: int
    n_bases: int
    adapter_reads: int

    def __init__(self) -> None:
        self.total_reads = 0
        self.total_bases = 0
        self.total_quality_sum = 0
        self.read_quality_distribution = {}
        self.position_quality_sums = []
        self.position_counts = []
        self.total_gc_bases = 0
        self.q20_bases = 0
        self.q30_bases = 0
        self.n_bases = 0
        self.adapter_reads = 0

    def accumulate(self, record: FASTQRecord) -> float:
        """Accumulate quality metrics from a single FASTQ record."""
        if not record.quality:
            return 0.0

        self.total_reads += 1
        quality_scores = [ord(char) - 33 for char in record.quality]

        seq_upper = record.sequence.upper()
        self.total_bases += len(seq_upper)
        self.total_gc_bases += seq_upper.count("G") + seq_upper.count("C")
        self.n_bases += seq_upper.count("N")

        for score in quality_scores:
            self.total_quality_sum += score
            if score >= 20:
                self.q20_bases += 1
            if score >= 30:
                self.q30_bases += 1

        if is_adapter_contaminated(seq_upper):
            self.adapter_reads += 1

        read_mean: float = sum(quality_scores) / len(quality_scores)
        rounded_mean = round(read_mean)
        self.read_quality_distribution[rounded_mean] = (
            self.read_quality_distribution.get(rounded_mean, 0) + 1
        )

        for i, score in enumerate(quality_scores):
            if i >= len(self.position_quality_sums):
                self.position_quality_sums.append(score)
                self.position_counts.append(1)
            else:
                self.position_quality_sums[i] += score
                self.position_counts[i] += 1

        return read_mean

    def summarize(self) -> SequenceQCMetricsSummary:
        """Summarize the accumulated metrics into a SequenceQCMetricsSummary."""
        if self.total_reads == 0:
            return SequenceQCMetricsSummary(
                total_reads_evaluated=0,
                total_bases_evaluated=0,
                mean_quality_score=0.0,
                positions_average_quality=[],
                read_quality_distribution={},
                gc_content_percent=0.0,
                q20_percent=0.0,
                q30_percent=0.0,
                n_content_percent=0.0,
                adapter_contaminated_reads=0,
                adapter_contamination_rate=0.0,
            )

        main_average_quality = (
            (self.total_quality_sum / self.total_bases) if self.total_bases > 0 else 0.0
        )

        positions_average_quality = [
            pos_sum / count
            for pos_sum, count in zip(self.position_quality_sums, self.position_counts, strict=True)
        ]

        gc_content_percent = (
            (self.total_gc_bases / self.total_bases * 100) if self.total_bases > 0 else 0.0
        )
        q20_percent = (self.q20_bases / self.total_bases * 100) if self.total_bases > 0 else 0.0
        q30_percent = (self.q30_bases / self.total_bases * 100) if self.total_bases > 0 else 0.0
        n_content_percent = (self.n_bases / self.total_bases * 100) if self.total_bases > 0 else 0.0
        adapter_contamination_rate = (
            (self.adapter_reads / self.total_reads * 100) if self.total_reads > 0 else 0.0
        )

        return SequenceQCMetricsSummary(
            total_reads_evaluated=self.total_reads,
            total_bases_evaluated=self.total_bases,
            mean_quality_score=main_average_quality,
            positions_average_quality=positions_average_quality,
            read_quality_distribution=dict(sorted(self.read_quality_distribution.items())),
            gc_content_percent=gc_content_percent,
            q20_percent=q20_percent,
            q30_percent=q30_percent,
            n_content_percent=n_content_percent,
            adapter_contaminated_reads=self.adapter_reads,
            adapter_contamination_rate=adapter_contamination_rate,
        )


def analyze_fastq_file(
    file_path: str,
    max_errors: int = 1000,
) -> tuple[SequenceQCMetricsSummary, list[FASTQValidationError]]:
    """Analyze a FASTQ file and return a summary of Sequence QC metrics.

    Caps memory usage by storing at most max_errors validation errors.
    """
    accumulator = SequenceQCAccumulator()
    validation_errors: list[FASTQValidationError] = []

    for record in parse_and_validate_fastq(file_path):
        if isinstance(record, FASTQValidationError):
            if len(validation_errors) < max_errors:
                validation_errors.append(record)
        else:
            accumulator.accumulate(record)

    return (accumulator.summarize(), validation_errors)
