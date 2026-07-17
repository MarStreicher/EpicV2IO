from pathlib import Path
from typing import Optional, Sequence

import pandas as pd
import pytest

from epicv2io import manifest


@pytest.mark.parametrize(
    ("columns", "expected_columns"),
    [
        (None, list(manifest.MANIFEST_COLUMNS)),
        (("IlmnID", "CHR"), ["IlmnID", "CHR"]),
        ([], []),
    ],
)
def test_load_peters_manifest_reads_packaged_parquet_with_requested_columns(
    columns: Optional[Sequence[str]],
    expected_columns: list[str],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    package_root = tmp_path / "epicv2io"
    resource_dir = package_root / "resources"
    resource_dir.mkdir(parents=True)
    parquet_path = resource_dir / "peters_epicv2_manifest.parquet"
    parquet_path.write_bytes(b"parquet-placeholder")
    observed = {}
    expected = pd.DataFrame({"IlmnID": ["cg1"]})

    monkeypatch.setattr(manifest.resources, "files", lambda package: package_root)

    def fake_read_parquet(file_object, *, columns):
        observed["contents"] = file_object.read()
        observed["columns"] = columns
        return expected

    monkeypatch.setattr(manifest.pd, "read_parquet", fake_read_parquet)

    result = manifest.load_peters_manifest(columns)

    assert result is expected
    assert observed == {
        "contents": b"parquet-placeholder",
        "columns": expected_columns,
    }
