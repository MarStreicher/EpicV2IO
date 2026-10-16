from __future__ import annotations

from dataclasses import astuple, dataclass
from typing import Dict, Set, Tuple, cast

import pandas as pd

_TRUTHY_MANIFEST_FLAGS = frozenset({"Y", "TRUE", "1"})

_FLAG_COLUMNS = (
    ("exclude_mismatch_pos", "MismatchPos"),
    ("exclude_missing_pos", "MissingPos"),
    ("exclude_ch_blat", "CH_BLAT"),
    ("exclude_ch_wgbs_evidence", "CH_WGBS_evidence"),
)

_RULE_ORDER = (
    "exclude_mismatch_pos",
    "exclude_missing_pos",
    "exclude_ch_blat",
    "exclude_ch_wgbs_evidence",
    "exclude_sex_chromosomes",
    "exclude_mitochondrial",
)


@dataclass(frozen=True)
class ExclusionFlags:
    exclude_mismatch_pos: bool = False
    exclude_missing_pos: bool = False
    exclude_ch_blat: bool = False
    exclude_ch_wgbs_evidence: bool = False
    exclude_sex_chromosomes: bool = False
    exclude_mitochondrial: bool = False

    def any_enabled(self) -> bool:
        return any(astuple(self))


def _truthy_manifest_mask(values: pd.Series) -> pd.Series:
    normalised = values.astype(str).str.strip().str.upper()
    return normalised.isin(_TRUTHY_MANIFEST_FLAGS)


def exclusion_sets_by_rule(
    manifest: pd.DataFrame,
    candidate_ids: pd.Index,
    flags: ExclusionFlags,
) -> Tuple[Dict[str, Set[str]], Set[str]]:
    cand_set = set(candidate_ids.astype(str))
    intersection_cg = manifest[manifest["IlmnID"].isin(cand_set)]
    ilmn = intersection_cg["IlmnID"]
    per_rule: Dict[str, Set[str]] = {}

    if intersection_cg.empty:
        raise ValueError(
            "No overlap between the loaded betas probes and the Peters et al. manifest (IlmnID)."
        )

    def _flag_rule(column: str) -> Set[str]:
        mask = _truthy_manifest_mask(cast(pd.Series, intersection_cg[column]))
        return set(ilmn[mask])

    for attr, column in _FLAG_COLUMNS:
        if getattr(flags, attr):
            per_rule[attr] = _flag_rule(column)

    if flags.exclude_sex_chromosomes or flags.exclude_mitochondrial:
        chr_u = (
            cast(pd.Series, intersection_cg["CHR"]).astype(str).str.strip().str.upper()
        )
        if flags.exclude_sex_chromosomes:
            per_rule["exclude_sex_chromosomes"] = set(
                ilmn[chr_u.isin(("CHRX", "CHRY"))]
            )
        if flags.exclude_mitochondrial:
            per_rule["exclude_mitochondrial"] = set(ilmn[chr_u == "CHRM"])

    union_excluded: Set[str] = set()
    for ids in per_rule.values():
        union_excluded |= ids
    return per_rule, union_excluded


def print_exclusion_summary(
    candidate_count: int,
    per_rule: Dict[str, Set[str]],
    union_excluded: Set[str],
) -> None:
    print(f"Candidate cg CpGs (betas ∩ manifest, cg rows): {candidate_count:,}")
    if not per_rule:
        print("  No manifest exclusion flags enabled; 0 CpGs removed by rules.")
        return
    for name in _RULE_ORDER:
        if name not in per_rule:
            continue
        n = len(per_rule[name])
        print(f"  {name}: {n:,} CpGs match this rule (among candidates)")
    print(f"  Union (unique CpGs removed by any enabled rule): {len(union_excluded):,}")
    print(f"  Remaining after union: {candidate_count - len(union_excluded):,}")
