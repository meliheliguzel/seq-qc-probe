import gzip
import sys
from collections.abc import Iterator
from pathlib import Path
from typing import TextIO

from .models import FASTQErrorKind, FASTQRecord, FASTQValidationError

VALID_BASES = frozenset("ACGTNacgtn")
_INVALID_BASE_TRANSLATOR = str.maketrans("", "", "ACGTNacgtn")


def _open_stream(file_path: str | Path) -> tuple[TextIO, bool]:
    """Opens a file path or stdin stream, handling gzip transparently.

    Returns the stream and a boolean indicating whether it should be closed upon exit.
    """
    path_str = str(file_path)
    if path_str == "-":
        return sys.stdin, False

    path = Path(path_str)
    with open(path, "rb") as f:
        magic = f.read(2)

    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace"), True
    return open(path, encoding="utf-8", errors="replace"), True


def parse_and_validate_fastq(
    file_path: str | Path,
) -> Iterator[FASTQRecord | FASTQValidationError]:
    """Parses a FASTQ file and yields valid FASTQRecord or FASTQValidationError objects.

    Enforces strict 4-line block invariant to prevent parser desynchronization.
    Memory footprint is O(1) streaming. Supports plain text, gzip (.gz), and stdin ('-').
    """
    stream, should_close = _open_stream(file_path)
    line_number = 0

    try:
        while True:
            raw_header = stream.readline()
            if not raw_header:
                break

            line_number += 1

            if raw_header in ("\n", "\r\n"):
                yield FASTQValidationError(
                    line_number=line_number,
                    error_kind=FASTQErrorKind.EMPTY_LINE,
                    message=f"Empty line encountered at line {line_number}",
                )
                continue

            start_line = line_number
            header = raw_header.rstrip("\r\n")

            raw_seq = stream.readline()
            if not raw_seq:
                yield FASTQValidationError(
                    line_number=line_number + 1,
                    error_kind=FASTQErrorKind.PREMATURE_EOF,
                    message=f"Unexpected end of file after header at line {line_number + 1}",
                )
                break
            line_number += 1
            sequence = raw_seq.rstrip("\r\n")

            raw_sep = stream.readline()
            if not raw_sep:
                yield FASTQValidationError(
                    line_number=line_number + 1,
                    error_kind=FASTQErrorKind.PREMATURE_EOF,
                    message=f"Unexpected end of file after sequence at line {line_number + 1}",
                )
                break
            line_number += 1
            separator = raw_sep.rstrip("\r\n")

            raw_qual = stream.readline()
            if not raw_qual:
                yield FASTQValidationError(
                    line_number=line_number + 1,
                    error_kind=FASTQErrorKind.PREMATURE_EOF,
                    message=f"Unexpected end of file after separator at line {line_number + 1}",
                )
                break
            line_number += 1
            quality = raw_qual.rstrip("\r\n")

            # --- (INVARIANT VALIDATION) ---
            rec = FASTQRecord(
                header=header,
                sequence=sequence,
                separator=separator,
                quality=quality,
                line_number=start_line,
            )
            has_error = False

            if not header.startswith("@"):
                yield FASTQValidationError(
                    line_number=start_line,
                    error_kind=FASTQErrorKind.INVALID_HEADER,
                    message=f"Invalid FASTQ header at line {start_line}: {header}",
                    record=rec,
                )
                has_error = True

            if len(sequence) == 0:
                yield FASTQValidationError(
                    line_number=start_line + 1,
                    error_kind=FASTQErrorKind.INVALID_SEQUENCE,
                    message=f"Empty sequence at line {start_line + 1}",
                    record=rec,
                )
                has_error = True
            elif bool(sequence.translate(_INVALID_BASE_TRANSLATOR)):
                yield FASTQValidationError(
                    line_number=start_line + 1,
                    error_kind=FASTQErrorKind.INVALID_SEQUENCE,
                    message=f"Invalid sequence at line {start_line + 1}: {sequence}",
                    record=rec,
                )
                has_error = True

            if not separator.startswith("+"):
                yield FASTQValidationError(
                    line_number=start_line + 2,
                    error_kind=FASTQErrorKind.INVALID_SEPARATOR,
                    message=f"Invalid FASTQ separator at line {start_line + 2}: {separator}",
                    record=rec,
                )
                has_error = True

            if len(sequence) != len(quality):
                yield FASTQValidationError(
                    line_number=start_line + 3,
                    error_kind=FASTQErrorKind.MISMATCHED_LENGTHS,
                    message=(
                        f"Sequence length ({len(sequence)}) and quality length "
                        f"({len(quality)}) mismatch at record starting at line {start_line}"
                    ),
                    record=rec,
                )
                has_error = True

            if quality and (min(quality) < "!" or max(quality) > "~"):
                yield FASTQValidationError(
                    line_number=start_line + 3,
                    error_kind=FASTQErrorKind.INVALID_QUALITY,
                    message=f"Invalid Phred quality character detected at line {start_line + 3}",
                    record=rec,
                )
                has_error = True

            if not has_error:
                yield rec

    finally:
        if should_close:
            stream.close()
