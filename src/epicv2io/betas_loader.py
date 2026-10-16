from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from .exclusions import (
    ExclusionFlags,
    exclusion_sets_by_rule,
    print_exclusion_summary,
)
from .io import (
    NON_SAMPLE_COLUMNS,
    filter_cg_rows,
    keep_sample_columns,
    read_betas_matrix,
    probes_to_samples,
)
from .manifest import load_peters_manifest

__all__ = ["BetasLoader", "NON_SAMPLE_COLUMNS"]


class BetasLoader:
    """
    Load EPICv2 betas and shipped Peters manifest, intersect on IlmnID, keep cg probes only.
    """

    def __init__(self, path: Path, manifest: Optional[pd.DataFrame] = None) -> None:
        self.path_betas = path
        self._manifest_override = manifest
        self._raw_data: Optional[pd.DataFrame] = None
        self._manifest: Optional[pd.DataFrame] = None

        self.valid_data = pd.DataFrame()

        self.load_candidate_cg_count: int = 0
        self.load_exclusion_counts_by_rule: Dict[str, int] = {}
        self.load_exclusion_union_count: int = 0

    def _load_manifest(self) -> pd.DataFrame:
        if self._manifest_override is None:
            manifest = load_peters_manifest()
        else:
            manifest = self._manifest_override.copy()
        manifest["IlmnID"] = manifest["IlmnID"].astype(str)
        return manifest

    @property
    def raw_data(self) -> pd.DataFrame:
        if self._raw_data is None:
            self._raw_data = read_betas_matrix(self.path_betas)
        return self._raw_data

    @property
    def manifest(self) -> pd.DataFrame:
        if self._manifest is None:
            self._manifest = self._load_manifest()
        return self._manifest

    def get_manifest_cpgs_missing_in_betas(self) -> List[str]:
        """
        Return IlmnIDs of manifest cg probes that are not present in the betas index.
        """
        ilmn = self.manifest["IlmnID"]
        manifest_cg = set(ilmn[ilmn.str.startswith("cg")])
        betas_index = set(self.raw_data.index.astype(str))
        return sorted(manifest_cg - betas_index)

    def load_data(self, **flags: bool) -> pd.DataFrame:
        """
        Load cg probes present in both betas and the shipped Peters manifest.
        """
        exclusions = ExclusionFlags(**flags)
        df_cg = filter_cg_rows(self.raw_data)
        self.load_candidate_cg_count = int(len(df_cg))

        if exclusions.any_enabled():
            per_rule, excluded = exclusion_sets_by_rule(
                self.manifest, df_cg.index, exclusions
            )
            self.load_exclusion_counts_by_rule = {
                key: len(ids) for key, ids in per_rule.items()
            }
            self.load_exclusion_union_count = len(excluded)
            print_exclusion_summary(self.load_candidate_cg_count, per_rule, excluded)
            self.valid_data = df_cg.loc[~df_cg.index.isin(excluded)]
        else:
            self.load_exclusion_counts_by_rule = {}
            self.load_exclusion_union_count = 0
            print_exclusion_summary(self.load_candidate_cg_count, {}, set())
            self.valid_data = df_cg

        print(self.valid_data.shape)
        return keep_sample_columns(self.valid_data)

    def load_data_wide(self, **flags: bool) -> pd.DataFrame:
        """
        Convenience wrapper around ``load_data()`` that returns a samples × features matrix.
        """
        return self.to_samples_by_features(self.load_data(**flags))

    @staticmethod
    def to_samples_by_features(df: pd.DataFrame) -> pd.DataFrame:
        """
        Convert a probes × samples matrix into samples × features.
        """
        return probes_to_samples(df)
