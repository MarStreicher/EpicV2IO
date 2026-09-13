"""Shared helpers for reading beta-value text files."""

from __future__ import annotations

from pathlib import Path

_SEPARATOR_BY_SUFFIX = {
    ".csv": ",",
    ".txt": "\t",
}


def infer_separator(path: Path) -> str:
    """Return the delimiter implied by the betas file extension.

    - ``.csv`` → comma (``,``)
    - ``.txt`` → tab (``\\t``)
    """
    suffix = path.suffix.lower()
    try:
        return _SEPARATOR_BY_SUFFIX[suffix]
    except KeyError as exc:
        supported = ", ".join(sorted(_SEPARATOR_BY_SUFFIX))
        raise ValueError(
            f"Unsupported betas file extension {suffix!r} for {path}; "
            f"expected one of: {supported}"
        ) from exc
