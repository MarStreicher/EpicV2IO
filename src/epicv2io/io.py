"""Shared helpers for reading beta-value text files."""

from __future__ import annotations

import csv
from pathlib import Path


def infer_separator(path: Path, sample_bytes: int = 64_384) -> str:
    """Infer a comma, tab, or semicolon delimiter from a text file."""
    with path.open("r", encoding="utf-8", errors="replace") as f:
        sample = f.read(sample_bytes)
    return csv.Sniffer().sniff(sample, delimiters=",\t;").delimiter
