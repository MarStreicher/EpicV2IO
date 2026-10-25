from pathlib import Path

import pandas as pd
import pytest

from epicv2io.io import (
    filter_cg_rows,
    infer_separator,
    keep_sample_columns,
    probes_to_samples,
)


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


def test_filter_cg_rows_drops_non_cg_prefixes() -> None:
    frame = pd.DataFrame(
        {"Sample_A": [0.1, 0.2, 0.3]},
        index=["cg1", "ch1", "rs1"],
    )

    assert filter_cg_rows(frame).index.tolist() == ["cg1"]


def test_keep_sample_columns_drops_annotation_and_controls() -> None:
    frame = pd.DataFrame(
        {
            "Sample_A": [0.1],
            "CHR": ["chr1"],
            "Control_note": ["x"],
        }
    )

    assert keep_sample_columns(frame).columns.tolist() == ["Sample_A"]


def test_to_samples_by_features_sets_axis_names() -> None:
    frame = pd.DataFrame({"Sample_A": [0.1]}, index=pd.Index(["cg1"], name="IlmnID"))

    wide = probes_to_samples(frame)

    assert wide.index.tolist() == ["Sample_A"]
    assert wide.columns.tolist() == ["cg1"]
    assert wide.index.name == "sample_id"
    assert wide.columns.name == "IlmnID"
