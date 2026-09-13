from pathlib import Path

import pytest

from epicv2io.io import infer_separator


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("betas.csv", ","),
        ("betas.CSV", ","),
        ("betas.txt", "\t"),
        ("betas.TXT", "\t"),
    ],
)
def test_infer_separator_from_extension(name: str, expected: str) -> None:
    assert infer_separator(Path(name)) == expected


def test_infer_separator_rejects_unknown_extension() -> None:
    with pytest.raises(ValueError, match="Unsupported betas file extension"):
        infer_separator(Path("betas.tsv"))
