from .models import FASTQRecord, FASTQValidationError, SequenceQCMetricsSummary
from .parser import parse_and_validate_fastq

ILLUMINA_UNIVERSAL_ADAPTER: str = "AGATCGGAAGAGC"


def is_adapter_contaminated(
    seq: str,
    adapter: str = ILLUMINA_UNIVERSAL_ADAPTER,
    min_overlap: int = 8,
    max_mismatches: int = 1,
) -> bool:
    """Detects Illumina adapter contamination using fast C-level filters and Pigeonhole principle.

    Supports partial 3' overlap or full match with mismatch tolerance.
    """
    seq_upper = seq.upper()
    adapter_len = len(adapter)
    seq_len = len(seq_upper)

    if seq_len >= adapter_len:
        if adapter in seq_upper:
            return True

        mid = adapter_len // 2
        part1 = adapter[:mid]
        part2 = adapter[mid:]
        if max_mismatches > 1 or (part1 in seq_upper or part2 in seq_upper):
            for i in range(seq_len - adapter_len + 1):
                sub = seq_upper[i : i + adapter_len]
                mismatches = 0
                for a, b in zip(sub, adapter, strict=True):
                    if a != b:
                        mismatches += 1
                        if mismatches > max_mismatches:
                            break
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
        """Accumulate quality metrics from a single FASTQ record in a single pass."""
        if not record.quality:
            return 0.0

        self.total_reads += 1

        seq_upper = record.sequence.upper()
        self.total_bases += len(seq_upper)
        self.total_gc_bases += seq_upper.count("G") + seq_upper.count("C")
        self.n_bases += seq_upper.count("N")

        if is_adapter_contaminated(seq_upper):
            self.adapter_reads += 1

        read_qual_sum = 0
        q20_add = 0
        q30_add = 0
        pos_sums = self.position_quality_sums
        pos_counts = self.position_counts
        pos_len = len(pos_sums)
        q_len = len(record.quality)

        for i, char in enumerate(record.quality):
            score = ord(char) - 33
            read_qual_sum += score
            if score >= 20:
                q20_add += 1
                if score >= 30:
                    q30_add += 1

            if i < pos_len:
                pos_sums[i] += score
                pos_counts[i] += 1
            else:
                pos_sums.append(score)
                pos_counts.append(1)
                pos_len += 1

        self.total_quality_sum += read_qual_sum
        self.q20_bases += q20_add
        self.q30_bases += q30_add

        read_mean: float = read_qual_sum / q_len
        rounded_mean = round(read_mean)
        self.read_quality_distribution[rounded_mean] = (
            self.read_quality_distribution.get(rounded_mean, 0) + 1
        )

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
                gc_bases=0,
                q20_bases=0,
                q30_bases=0,
                n_bases=0,
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
            gc_bases=self.total_gc_bases,
            q20_bases=self.q20_bases,
            q30_bases=self.q30_bases,
            n_bases=self.n_bases,
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
