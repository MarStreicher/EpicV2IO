from __future__ import annotations

from pathlib import Path
from typing import cast

import pandas as pd

_SEPARATOR_BY_SUFFIX = {
    ".csv": ",",
    ".txt": "\t",
}

NON_SAMPLE_COLUMNS = [
    "TargetID",
    "CHR",
    "MAPINFO",
    "STRAND",
    "UCSC_REFGENE_NAME",
    "UCSC_REFGENE_ACCESSION",
    "UCSC_REFGENE_GROUP",
    "UCSC_CPG_ISLANDS_NAME",
    "RELATION_TO_UCSC_CPG_ISLAND",
    "PHANTOM4_ENHANCERS",
    "PHANTOM5_ENHANCERS",
    "DMR",
    "450K_ENHANCER",
    "REGULATORY_FEATURE_NAME",
    "REGULATORY_FEATURE_GROUP",
    "GENCODEBASICV12_NAME",
    "GENCODEBASICV12_ACCESSION",
    "GENCODEBASICV12_GROUP",
    "GENCODECOMPV12_NAME",
    "GENCODECOMPV12_ACCESSION",
    "GENCODECOMPV12_GROUP",
    "DNASE_HYPERSENSITIVITY_NAME",
    "DNASE_HYPERSENSITIVITY_EVIDENCE_COUNT",
    "TFBS_NAME",
    "TFBS_EVIDENCE_COUNT",
    "METHYL450_LOCI",
    "SNP_ID",
    "SNP_DISTANCE",
    "SNP_MINORALLELEFREQUENCY",
    "CHR_HG38",
    "START_HG38",
    "END_HG38",
    "STRAND_HG38",
]


def infer_separator(path: Path) -> str:
    suffix = path.suffix.lower()
    try:
        return _SEPARATOR_BY_SUFFIX[suffix]
    except KeyError as exc:
        supported = ", ".join(sorted(_SEPARATOR_BY_SUFFIX))
        raise ValueError(
            f"Unsupported betas file extension {suffix!r} for {path}; "
            f"expected one of: {supported}"
        ) from exc


def read_betas_matrix(path: Path) -> pd.DataFrame:
    raw_data = pd.read_csv(
        path, sep=infer_separator(path), low_memory=False, index_col=0
    )
    if not raw_data.index.astype(str).str.startswith("cg").any():
        raise ValueError(
            f"Expected probe IDs in index, got: {raw_data.index[:5].tolist()}"
        )
    raw_data.index = raw_data.index.astype(str)
    raw_data.index.name = "IlmnID"
    return raw_data


def filter_cg_rows(df: pd.DataFrame) -> pd.DataFrame:
    return cast(pd.DataFrame, df.loc[df.index.str.startswith("cg")])


def keep_sample_columns(df: pd.DataFrame) -> pd.DataFrame:
    sample_columns = [
        c
        for c in df.columns
        if c not in NON_SAMPLE_COLUMNS and not str(c).startswith("Control_")
    ]
    return cast(pd.DataFrame, df.loc[:, sample_columns])


def probes_to_samples(df: pd.DataFrame) -> pd.DataFrame:
    wide = df.T
    wide.index.name = "sample_id"
    wide.columns.name = "IlmnID"
    return wide
