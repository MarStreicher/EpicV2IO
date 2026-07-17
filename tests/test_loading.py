from pathlib import Path
from typing import Optional

import pandas as pd
import pytest

from epicv2io import BetasLoader
from epicv2io import betas_loader


@pytest.mark.parametrize(
    ("delimiter", "suffix"),
    [
        (";", ".csv"),
        ("\t", ".tsv"),
        (",", ".csv"),
    ],
    ids=["semicolon", "tab", "comma"],
)
def test_init_infers_delimiter(
    delimiter: str,
    suffix: str,
    example_betas_csv: Path,
    example_manifest: pd.DataFrame,
    tmp_path: Path,
) -> None:
    """The loader should infer the delimiter and preserve probe IDs as its index."""
    input_path = tmp_path / f"example_betas{suffix}"

    contents = example_betas_csv.read_text().replace(";", delimiter)
    input_path.write_text(contents)

    loader = BetasLoader(input_path, manifest=example_manifest)

    assert loader.raw_data.index.tolist() == [
        "cg00000029_TC11",
        "cg00000103_BC11",
        "ch00000001",
        "rs00000001",
    ]
    assert loader.raw_data.loc["cg00000029_TC11", "Sample_A"] == pytest.approx(0.10)
    assert loader.manifest is not example_manifest


def test_invalid_probe_index_error(
    invalid_betas_csv: Path, example_manifest: pd.DataFrame
) -> None:
    with pytest.raises(ValueError, match="Expected probe IDs"):
        BetasLoader(invalid_betas_csv, manifest=example_manifest)


def test_keep_only_cg_probes(
    example_betas_csv: Path, example_manifest: pd.DataFrame
) -> None:
    loader = BetasLoader(example_betas_csv, manifest=example_manifest)
    data = loader.load_data()

    assert data.index.tolist() == ["cg00000029_TC11", "cg00000103_BC11"]


def test_load_data_keeps_only_sample_columns(
    example_betas_csv: Path, example_manifest: pd.DataFrame
) -> None:
    loader = BetasLoader(example_betas_csv, manifest=example_manifest)

    data = loader.load_data()

    assert data.columns.tolist() == ["Sample_A", "Sample_B"]
    assert loader.valid_data.columns.tolist() == [
        "Sample_A",
        "Sample_B",
        "CHR",
        "Control_note",
    ]


def test_load_data_wide_sets_axes_and_transposes_values(
    example_betas_csv: Path, example_manifest: pd.DataFrame
) -> None:
    loader = BetasLoader(example_betas_csv, manifest=example_manifest)

    data = loader.load_data_wide()

    assert data.index.tolist() == ["Sample_A", "Sample_B"]
    assert data.columns.tolist() == ["cg00000029_TC11", "cg00000103_BC11"]
    assert data.index.name == "sample_id"
    assert data.columns.name == "IlmnID"
    assert data.loc["Sample_B", "cg00000103_BC11"] == pytest.approx(0.25)


def test_get_manifest_cpgs_missing_in_betas_returns_sorted_ids(
    example_betas_csv: Path, example_manifest: pd.DataFrame
) -> None:
    manifest = pd.concat(
        [
            example_manifest,
            pd.DataFrame(
                {
                    "IlmnID": ["cg_aaa_TC11", "rs_not_a_cpg"],
                    "CHR": ["chr2", "chr2"],
                    "MismatchPos": ["N", "N"],
                    "MissingPos": ["N", "N"],
                    "CH_BLAT": ["N", "N"],
                    "CH_WGBS_evidence": ["N", "N"],
                }
            ),
        ],
        ignore_index=True,
    )
    loader = BetasLoader(example_betas_csv, manifest=manifest)

    assert loader.get_manifest_cpgs_missing_in_betas() == [
        "cg_aaa_TC11",
        "cg_not_present_TC11",
    ]


@pytest.mark.parametrize(
    ("flag", "manifest_column", "flag_value", "chromosome"),
    [
        ("exclude_mismatch_pos", "MismatchPos", " y ", "chr1"),
        ("exclude_missing_pos", "MissingPos", "TRUE", "chr1"),
        ("exclude_ch_blat", "CH_BLAT", "1", "chr1"),
        ("exclude_ch_wgbs_evidence", "CH_WGBS_evidence", "Y", "chr1"),
        ("exclude_sex_chromosomes", None, None, " chrX "),
        ("exclude_mitochondrial", None, None, "chrM"),
    ],
)
def test_each_manifest_exclusion_rule_removes_matching_probe(
    flag: str,
    manifest_column: Optional[str],
    flag_value: object,
    chromosome: str,
    example_betas_csv: Path,
    example_manifest: pd.DataFrame,
) -> None:
    manifest = example_manifest.copy()
    manifest.loc[0, "CHR"] = chromosome
    manifest.loc[1, "CHR"] = "chr2"
    if manifest_column is not None:
        manifest.loc[0, manifest_column] = flag_value
        manifest.loc[1, manifest_column] = "N"
    loader = BetasLoader(example_betas_csv, manifest=manifest)

    data = loader.load_data(**{flag: True})

    assert data.index.tolist() == ["cg00000103_BC11"]
    assert loader.load_candidate_cg_count == 2
    assert loader.load_exclusion_counts_by_rule == {flag: 1}
    assert loader.load_exclusion_union_count == 1


def test_overlapping_exclusion_rules_count_unique_removed_probes(
    example_betas_csv: Path, example_manifest: pd.DataFrame
) -> None:
    manifest = example_manifest.copy()
    manifest["MismatchPos"] = "N"
    manifest.loc[0, ["MismatchPos", "MissingPos"]] = "Y"
    manifest.loc[1, "MissingPos"] = "Y"
    loader = BetasLoader(example_betas_csv, manifest=manifest)

    data = loader.load_data(exclude_mismatch_pos=True, exclude_missing_pos=True)

    assert data.empty
    assert loader.load_exclusion_counts_by_rule == {
        "exclude_mismatch_pos": 1,
        "exclude_missing_pos": 2,
    }
    assert loader.load_exclusion_union_count == 2


def test_load_data_without_filters_resets_exclusion_statistics(
    example_betas_csv: Path, example_manifest: pd.DataFrame
) -> None:
    loader = BetasLoader(example_betas_csv, manifest=example_manifest)
    loader.load_data(exclude_missing_pos=True)

    loader.load_data()

    assert loader.load_exclusion_counts_by_rule == {}
    assert loader.load_exclusion_union_count == 0


def test_manifest_filter_requires_overlap_with_betas(
    example_betas_csv: Path, example_manifest: pd.DataFrame
) -> None:
    manifest = example_manifest.copy()
    manifest["IlmnID"] = ["cg_absent_1", "cg_absent_2", "cg_absent_3"]
    loader = BetasLoader(example_betas_csv, manifest=manifest)

    with pytest.raises(ValueError, match="No overlap"):
        loader.load_data(exclude_missing_pos=True)


def test_loader_uses_a_copy_of_manifest_and_normalises_ids(
    example_betas_csv: Path, example_manifest: pd.DataFrame
) -> None:
    manifest = example_manifest.copy()
    manifest["IlmnID"] = pd.Series([1, 2, 3], dtype="int64")

    loader = BetasLoader(example_betas_csv, manifest=manifest)

    assert loader.manifest["IlmnID"].tolist() == ["1", "2", "3"]
    assert manifest["IlmnID"].tolist() == [1, 2, 3]


def test_loader_uses_packaged_manifest_by_default(
    example_betas_csv: Path,
    example_manifest: pd.DataFrame,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls = []

    def fake_load_peters_manifest() -> pd.DataFrame:
        calls.append(True)
        return example_manifest.copy()

    monkeypatch.setattr(
        betas_loader, "load_peters_manifest", fake_load_peters_manifest
    )

    loader = BetasLoader(example_betas_csv)

    assert calls == [True]
    assert loader.manifest.equals(example_manifest)

