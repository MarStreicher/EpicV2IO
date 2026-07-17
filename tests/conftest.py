from pathlib import Path

import pandas as pd
import pytest


@pytest.fixture
def example_betas_csv(tmp_path: Path) -> Path:
    """Create a small semicolon-delimited beta file for one test run."""
    betas = pd.DataFrame(
        {
            "Sample_A": [0.10, 0.20, 0.30, 0.40],
            "Sample_B": [0.15, 0.25, 0.35, 0.45],
            "CHR": ["chr1", "chrX", "chr2", "chr3"],
            "Control_note": ["keep out"] * 4,
        },
        index=pd.Index(
            [
                "cg00000029_TC11",
                "cg00000103_BC11",
                "ch00000001",
                "rs00000001",
            ]
        ),
    )
    path = tmp_path / "example_betas.csv"
    betas.to_csv(path, sep=";", index_label="IlmnID")
    return path


@pytest.fixture
def invalid_betas_csv(tmp_path: Path) -> Path:
    """Create a small invalid semicolon-delimited beta file for one test run."""
    betas = pd.DataFrame(
        {
            "Sample_A": [0.10, 0.20, 0.30, 0.40],
            "Sample_B": [0.15, 0.25, 0.35, 0.45],
            "CHR": ["chr1", "chrX", "chr2", "chr3"],
            "Control_note": ["keep out"] * 4,
        },
        index=pd.Index(
            [
                "ch00000029_TC11",
                "rs00000103_BC11",
                "ch00000001",
                "rs00000001",
            ]
        ),
    )
    path = tmp_path / "example_betas.csv"
    betas.to_csv(path, sep=";", index_label="IlmnID")
    return path


@pytest.fixture
def example_manifest() -> pd.DataFrame:
    """Return a minimal manifest with rows exercising different filters."""
    return pd.DataFrame(
        {
            "IlmnID": [
                "cg00000029_TC11",
                "cg00000103_BC11",
                "cg_not_present_TC11",
            ],
            "CHR": ["chr1", "chrX", "chrM"],
            "MismatchPos": ["N", "Y", "N"],
            "MissingPos": ["N", "N", "Y"],
            "CH_BLAT": ["N", "N", "Y"],
            "CH_WGBS_evidence": ["N", "N", "N"],
        }
    )
