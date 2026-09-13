import pandas as pd

from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, cast

from .io import infer_separator
from .manifest import load_peters_manifest

_TRUTHY_MANIFEST_FLAGS = frozenset({"Y", "TRUE", "1"})

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


class BetasLoader:
    """
    Load EPICv2 betas and shipped Peters manifest, intersect on IlmnID, keep cg probes only.
    """

    def __init__(self, path: Path, manifest: Optional[pd.DataFrame] = None) -> None:
        self.path_betas = path

        sep = infer_separator(self.path_betas)
        self.raw_data = pd.read_csv(
            self.path_betas, sep=sep, low_memory=False, index_col=0
        )
        if not self.raw_data.index.astype(str).str.startswith("cg").any():
            raise ValueError(
                f"Expected probe IDs in index, got: {self.raw_data.index[:5].tolist()}"
            )

        self.raw_data.index = self.raw_data.index.astype(str)
        self.raw_data.index.name = "IlmnID"

        self.manifest = load_peters_manifest() if manifest is None else manifest.copy()
        self.manifest["IlmnID"] = self.manifest["IlmnID"].astype(str)

        self.valid_data = pd.DataFrame()

        self.load_candidate_cg_count: int = 0
        self.load_exclusion_counts_by_rule: Dict[str, int] = {}
        self.load_exclusion_union_count: int = 0

    @staticmethod
    def _truthy_manifest_mask(values: pd.Series) -> pd.Series:
        normalised = values.astype(str).str.strip().str.upper()
        return normalised.isin(_TRUTHY_MANIFEST_FLAGS)

    def _get_mismatch_pos(self, df: pd.DataFrame) -> List[str]:
        """cg probes with MismatchPos flagged (Illumina vs SeSAMe location discrepancy)."""
        ilmn = cast(pd.Series, df["IlmnID"]).astype(str)
        mismatch = cast(pd.Series, df["MismatchPos"])
        m = ilmn.str.startswith("cg") & self._truthy_manifest_mask(mismatch)
        return ilmn[m].tolist()

    def _manifest_exclusion_sets_by_rule(
        self,
        candidate_ids: pd.Index,
        *,
        exclude_mismatch_pos: bool,
        exclude_missing_pos: bool,
        exclude_ch_blat: bool,
        exclude_ch_wgbs_evidence: bool,
        exclude_sex_chromosomes: bool,
        exclude_mitochondrial: bool,
    ) -> Tuple[Dict[str, Set[str]], Set[str]]:
        """
        For each enabled exclusion flag, IlmnIDs among candidates matching that rule.
        Returns (per_rule_sets, union). Per-rule counts can overlap; union is unique removed.
        """
        cand_set = set(candidate_ids.astype(str))
        intersection_cg = self.manifest[
            self.manifest["IlmnID"].isin(cast(pd.Series, cand_set))
        ]
        ilmn = intersection_cg["IlmnID"]
        per_rule: Dict[str, Set[str]] = {}

        if intersection_cg.empty:
            raise ValueError(
                "No overlap between the loaded betas probes and the Peters et al. manifest (IlmnID)."
            )

        def _flag_rule(column: str) -> Set[str]:
            mask = self._truthy_manifest_mask(cast(pd.Series, intersection_cg[column]))
            return set(ilmn[mask])

        if exclude_mismatch_pos:
            per_rule["exclude_mismatch_pos"] = _flag_rule("MismatchPos")
        if exclude_missing_pos:
            per_rule["exclude_missing_pos"] = _flag_rule("MissingPos")
        if exclude_ch_blat:
            per_rule["exclude_ch_blat"] = _flag_rule("CH_BLAT")
        if exclude_ch_wgbs_evidence:
            per_rule["exclude_ch_wgbs_evidence"] = _flag_rule("CH_WGBS_evidence")

        # Only normalise CHR when a chromosome-based rule is actually requested.
        if exclude_sex_chromosomes or exclude_mitochondrial:
            chr_u = (
                cast(pd.Series, intersection_cg["CHR"])
                .astype(str)
                .str.strip()
                .str.upper()
            )
            if exclude_sex_chromosomes:
                per_rule["exclude_sex_chromosomes"] = set(
                    ilmn[chr_u.isin(("CHRX", "CHRY"))]
                )
            if exclude_mitochondrial:
                per_rule["exclude_mitochondrial"] = set(ilmn[chr_u == "CHRM"])

        union_excluded: Set[str] = set()
        for s in per_rule.values():
            union_excluded |= s
        return per_rule, union_excluded

    def _print_exclusion_summary(
        self,
        candidate_count: int,
        per_rule: Dict[str, Set[str]],
        union_excluded: Set[str],
    ) -> None:
        print(f"Candidate cg CpGs (betas ∩ manifest, cg rows): {candidate_count:,}")
        if not per_rule:
            print("  No manifest exclusion flags enabled; 0 CpGs removed by rules.")
            return
        rule_order = (
            "exclude_mismatch_pos",
            "exclude_missing_pos",
            "exclude_ch_blat",
            "exclude_ch_wgbs_evidence",
            "exclude_sex_chromosomes",
            "exclude_mitochondrial",
        )
        for name in rule_order:
            if name not in per_rule:
                continue
            n = len(per_rule[name])
            print(f"  {name}: {n:,} CpGs match this rule (among candidates)")
        print(
            f"  Union (unique CpGs removed by any enabled rule): {len(union_excluded):,}"
        )
        print(f"  Remaining after union: {candidate_count - len(union_excluded):,}")

    def get_manifest_cpgs_missing_in_betas(self) -> List[str]:
        """
        Return IlmnIDs of manifest cg probes that are not present in the betas index.
        """
        ilmn = self.manifest["IlmnID"]
        manifest_cg = set(ilmn[ilmn.str.startswith("cg")])
        betas_index = set(self.raw_data.index.astype(str))
        return sorted(manifest_cg - betas_index)

    def _filter_cg_rows(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Removal of all non cg rows in the data frame.

        cg - Methylation at CpG sites
        ch - Methylation at CpH sites
        nv - Nucleotide variants (SNPs without dbSNP rsIDs)
        rs - SNPs with dbSNP rsIDs
        """
        return cast(pd.DataFrame, df.loc[df.index.str.startswith("cg")])

    def _keep_sample_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        sample_columns = [
            c
            for c in df.columns
            if c not in NON_SAMPLE_COLUMNS and not c.startswith("Control_")
        ]
        return cast(pd.DataFrame, df.loc[:, sample_columns])

    @staticmethod
    def to_samples_by_features(df: pd.DataFrame) -> pd.DataFrame:
        """
        Convert a probe/features × samples matrix into samples × features.
        """
        wide = df.T
        wide.index.name = "sample_id"
        wide.columns.name = "IlmnID"
        return wide

    def _prepare_cg_frame(
        self,
        *,
        exclude_mismatch_pos: bool = False,
        exclude_missing_pos: bool = False,
        exclude_ch_blat: bool = False,
        exclude_ch_wgbs_evidence: bool = False,
        exclude_sex_chromosomes: bool = False,
        exclude_mitochondrial: bool = False,
    ) -> pd.DataFrame:
        df_cg = self._filter_cg_rows(self.raw_data)

        self.load_candidate_cg_count = int(len(df_cg))
        any_manifest_filter = any(
            (
                exclude_mismatch_pos,
                exclude_missing_pos,
                exclude_ch_blat,
                exclude_ch_wgbs_evidence,
                exclude_sex_chromosomes,
                exclude_mitochondrial,
            )
        )

        if any_manifest_filter:
            per_rule, excluded = self._manifest_exclusion_sets_by_rule(
                df_cg.index,
                exclude_mismatch_pos=exclude_mismatch_pos,
                exclude_missing_pos=exclude_missing_pos,
                exclude_ch_blat=exclude_ch_blat,
                exclude_ch_wgbs_evidence=exclude_ch_wgbs_evidence,
                exclude_sex_chromosomes=exclude_sex_chromosomes,
                exclude_mitochondrial=exclude_mitochondrial,
            )
            self.load_exclusion_counts_by_rule = {
                k: len(v) for k, v in per_rule.items()
            }
            self.load_exclusion_union_count = len(excluded)
            self._print_exclusion_summary(
                self.load_candidate_cg_count, per_rule, excluded
            )
            self.valid_data = df_cg.loc[~df_cg.index.isin(excluded)]
        else:
            self.load_exclusion_counts_by_rule = {}
            self.load_exclusion_union_count = 0
            self._print_exclusion_summary(self.load_candidate_cg_count, {}, set())
            self.valid_data = df_cg

        print(self.valid_data.shape)

        out = self._keep_sample_columns(self.valid_data)
        return out

    def load_data(
        self,
        *,
        exclude_mismatch_pos: bool = False,
        exclude_missing_pos: bool = False,
        exclude_ch_blat: bool = False,
        exclude_ch_wgbs_evidence: bool = False,
        exclude_sex_chromosomes: bool = False,
        exclude_mitochondrial: bool = False,
    ) -> pd.DataFrame:
        """
        Load cg probes present in both betas and the shipped Peters manifest.
        """
        return self._prepare_cg_frame(
            exclude_mismatch_pos=exclude_mismatch_pos,
            exclude_missing_pos=exclude_missing_pos,
            exclude_ch_blat=exclude_ch_blat,
            exclude_ch_wgbs_evidence=exclude_ch_wgbs_evidence,
            exclude_sex_chromosomes=exclude_sex_chromosomes,
            exclude_mitochondrial=exclude_mitochondrial,
        )

    def load_data_wide(
        self,
        *,
        exclude_mismatch_pos: bool = False,
        exclude_missing_pos: bool = False,
        exclude_ch_blat: bool = False,
        exclude_ch_wgbs_evidence: bool = False,
        exclude_sex_chromosomes: bool = False,
        exclude_mitochondrial: bool = False,
    ) -> pd.DataFrame:
        """
        Convenience wrapper around ``load_data()`` that returns a samples × features matrix.
        """
        return self.to_samples_by_features(
            self.load_data(
                exclude_mismatch_pos=exclude_mismatch_pos,
                exclude_missing_pos=exclude_missing_pos,
                exclude_ch_blat=exclude_ch_blat,
                exclude_ch_wgbs_evidence=exclude_ch_wgbs_evidence,
                exclude_sex_chromosomes=exclude_sex_chromosomes,
                exclude_mitochondrial=exclude_mitochondrial,
            )
        )
