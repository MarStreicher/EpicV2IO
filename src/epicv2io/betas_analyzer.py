from pathlib import Path

import numpy as np
import pandas as pd
import sys

from .io import infer_separator
from .manifest import load_peters_manifest

_PROBE_TYPE_PREFIXES = ("cg", "ch", "rs", "nv", "control")
_SEX_CHROMOSOMES = ("chrX", "chrY")


def _print_counts(
    title: str,
    categories: tuple[str, ...],
    counts: pd.Series,
    n_rows: int,
    *,
    extra_labels: tuple[str, ...] = (),
) -> None:
    print(title)
    print("-" * 40)
    for label in categories:
        c = int(counts.get(label, 0))
        pct = 100.0 * c / n_rows if n_rows else 0.0
        print(f"  {label:10} {c:>10,}  ({pct:5.2f}%)")
    for label in extra_labels:
        c = int(counts.get(label, 0))
        if c:
            pct = 100.0 * c / n_rows if n_rows else 0.0
            print(f"  {label:10} {c:>10,}  ({pct:5.2f}%)")
    print("-" * 40)
    print(f"  {'total':10} {n_rows:>10,}  (100.00%)")
    print()


class BetasAnalyzer:
    def __init__(self, betas_path: Path) -> None:
        self._betas_path = betas_path
        if not betas_path.is_file():
            print(f"Error: betas file not found: {betas_path}", file=sys.stderr)
            sys.exit(1)

    def get_header_and_sample_count(self) -> tuple[list[str], int]:
        sep = infer_separator(self._betas_path)
        empty_frame = pd.read_csv(self._betas_path, sep=sep, nrows=0)
        header = empty_frame.columns.tolist()
        return (header, len(header))

    def _load_manifest(self) -> pd.DataFrame:
        manifest = load_peters_manifest(columns=("IlmnID", "CHR"))
        return manifest.astype(str)

    def summarise(self) -> None:
        header_columns, n_columns = self.get_header_and_sample_count()
        print(f"Betas file: {self._betas_path}")
        print(f"Header (first 10 columns): {header_columns[:10]}...")
        print(f"Number of columns (samples + index): {n_columns}")
        print(f"Number of sample columns: {n_columns - 1}")
        print()

        sep = infer_separator(self._betas_path)
        probes = pd.read_csv(
            self._betas_path,
            sep=sep,
            usecols=[0],
            index_col=0,
            dtype=str,
        )
        ids = probes.index.str.strip().str.lower()
        n_rows = len(ids)

        conditions = [ids.str.startswith(prefix) for prefix in _PROBE_TYPE_PREFIXES]
        probe_type = pd.Series(
            np.select(conditions, _PROBE_TYPE_PREFIXES, default="other"),
            index=ids,
            name="probe_type",
        )

        manifest = self._load_manifest()
        unique_chr = sorted(manifest["CHR"].dropna().astype(str).unique().tolist())

        manifest["IlmnID"] = manifest["IlmnID"].str.strip().str.lower()
        probe_to_chr = manifest.set_index("IlmnID")["CHR"]
        chr_series = ids.to_series().map(probe_to_chr).fillna("undefined")

        type_counts = probe_type.value_counts()
        chr_counts = chr_series.value_counts()

        print(f"Total probe rows (index size): {n_rows:,}")
        print()
        _print_counts(
            "Probe type counts (from index prefix):",
            _PROBE_TYPE_PREFIXES,
            type_counts,
            n_rows,
            extra_labels=("other",),
        )
        _print_counts(
            "Probe chr counts (from manifest):",
            tuple(unique_chr),
            chr_counts,
            n_rows,
            extra_labels=("undefined",),
        )

        print("Distribution of autosomal vs sex chromosomes (from manifest):")
        print("-" * 40)
        total_sex = int(chr_counts.reindex(_SEX_CHROMOSOMES, fill_value=0).sum())
        pct = 100.0 * total_sex / n_rows if n_rows else 0.0
        print(f"  {'sex':10} {total_sex:>10,}  ({pct:5.2f}%)")
        total_autosomal = n_rows - total_sex
        pct = 100.0 * total_autosomal / n_rows if n_rows else 0.0
        print(f"  {'autosomal':10} {total_autosomal:>10,}  ({pct:5.2f}%)")
        print("-" * 40)
        print(f"  {'total':10} {n_rows:>10,}  (100.00%)")
        print()
