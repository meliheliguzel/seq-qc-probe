"""Input/Output streaming utilities for file and standard input handling."""

import gzip
import sys
from pathlib import Path
from typing import TextIO


def open_stream(file_path: str | Path) -> tuple[TextIO, bool]:
    """Opens a file path or stdin stream, handling gzip transparently.

    Args:
        file_path: Path to the target file, or '-' for standard input.

    Returns:
        A tuple of (stream, should_close). The boolean flag indicates
        whether the caller should close the stream upon exit (False for stdin).
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
