"""Load the packaged Peters et al. EPICv2 probe annotation subset.

The Parquet file shipped under ``epicv2io/resources/`` is a *derived*
column subset of the extended EPICv2 manifest from Peters et al. (2024).

Source / provenance
-------------------
- Paper (CC BY 4.0):
  Peters TJ, Meyer B, Ryan L, et al. Characterisation and reproducibility
  of the HumanMethylationEPIC v2.0 BeadChip for DNA methylation profiling.
  *BMC Genomics* (2024) 25:251.
  https://doi.org/10.1186/s12864-024-10027-5
- Supplementary CSV: Additional file 4 of that article.
- Bioconductor / AnnotationHub: package ``EPICv2manifest``, record
  ``AH116484`` (Artistic-2.0).

Only these columns are retained for EpicV2IO:
``IlmnID``, ``CHR``, ``MismatchPos``, ``MissingPos``, ``CH_BLAT``,
``CH_WGBS_evidence``.

When using this manifest in publications, cite the Peters et al. (2024)
BMC Genomics paper (and prefer also acknowledging ``EPICv2manifest`` /
AH116484). See the package README / NOTICE for full attribution.
"""

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
    """Load the packaged Peters et al. EPICv2 annotation subset.

    Returns a DataFrame with the requested columns (default: all
    ``MANIFEST_COLUMNS``). This is third-party annotation data distributed
    as a derived subset under the paper's CC BY 4.0 terms; cite Peters et
    al., *BMC Genomics* (2024) 25:251 when you use it.
    """
    usecols = list(columns) if columns is not None else list(MANIFEST_COLUMNS)
    with (
        resources.files("epicv2io")
        .joinpath("resources", _MANIFEST_FILENAME)
        .open("rb") as f
    ):
        return pd.read_parquet(f, columns=usecols)
