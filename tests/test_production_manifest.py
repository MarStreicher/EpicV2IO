from pathlib import Path

import pandas as pd
import pytest

from epicv2io import BetasLoader, CgProbeId
from epicv2io.manifest import MANIFEST_COLUMNS, load_peters_manifest


@pytest.fixture(scope="module")
def production_manifest() -> pd.DataFrame:
    return load_peters_manifest()


@pytest.mark.integration
def test_packaged_manifest_has_expected_schema_and_counts(
    production_manifest: pd.DataFrame,
) -> None:
    assert production_manifest.columns.tolist() == list(MANIFEST_COLUMNS)
    assert len(production_manifest) == 936_166

    cg_rows = production_manifest["IlmnID"].astype(str).str.startswith("cg")
    mismatch_rows = (
        production_manifest["MismatchPos"]
        .astype(str)
        .str.strip()
        .str.upper()
        .isin({"Y", "TRUE", "1"})
    )

    assert cg_rows.sum() == 933_252
    assert (cg_rows & mismatch_rows).sum() == 7_507


@pytest.mark.integration
def test_all_production_cpg_ids_are_parseable(
    production_manifest: pd.DataFrame,
) -> None:
    cg_ids = production_manifest.loc[
        production_manifest["IlmnID"].astype(str).str.startswith("cg"), "IlmnID"
    ]

    invalid_ids = [probe_id for probe_id in cg_ids if CgProbeId.parse(probe_id) is None]

    assert invalid_ids == []


@pytest.mark.integration
def test_loader_uses_packaged_manifest_for_exclusions(
    production_manifest: pd.DataFrame, tmp_path: Path
) -> None:
    manifest = production_manifest.copy()
    cg_rows = manifest["IlmnID"].astype(str).str.startswith("cg")
    mismatch_rows = (
        manifest["MismatchPos"]
        .astype(str)
        .str.strip()
        .str.upper()
        .isin({"Y", "TRUE", "1"})
    )
    mismatch_id = manifest.loc[cg_rows & mismatch_rows, "IlmnID"].iloc[0]
    retained_id = manifest.loc[cg_rows & ~mismatch_rows, "IlmnID"].iloc[0]
    path = tmp_path / "production_manifest_smoke.tsv"
    pd.DataFrame(
        {"Sample_A": [0.1, 0.2]},
        index=pd.Index([mismatch_id, retained_id], name="IlmnID"),
    ).to_csv(path, sep="\t")

    loader = BetasLoader(path)
    data = loader.load_data(exclude_mismatch_pos=True)

    assert data.index.tolist() == [retained_id]
    assert loader.load_exclusion_counts_by_rule == {"exclude_mismatch_pos": 1}
