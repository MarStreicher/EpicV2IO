from importlib import resources
from typing import Sequence, Optional

import pandas as pd

_MANIFEST_FILENAME = "peters_epicv2_manifest.parquet"

MANIFEST_COLUMNS: tuple[str, ...] = (
    "IlmnID",
    "CHR",
    "MismatchPos",
    "MissingPos",
    "CH_BLAT",
    "CH_WGBS_evidence",
)


def load_peters_manifest(
    columns: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    usecols = list(columns) if columns is not None else list(MANIFEST_COLUMNS)
    with (
        resources.files("epicv2io")
        .joinpath("resources", _MANIFEST_FILENAME)
        .open("rb") as f
    ):
        return pd.read_parquet(f, columns=usecols)
