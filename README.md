<p align="right">
  <img src="https://raw.githubusercontent.com/MarStreicher/EpicV2IO/main/assets/epicv2io_logo.png" alt="epicv2io logo" width="250">
</p>

# EpicV2IO

I/O library for Illumina MethylationEPIC v2 (EPICv2) beta arrays.

EpicV2IO loads beta files, retains cg probes, optionally applies QC exclusions from the shipped Peters et al. probe annotation manifest, and provides helpers for probe-ID parsing and file summaries. It is designed as a lightweight dependency for analysis packages (for example [betaSieve](https://github.com/MarStreicher/betaSieve)) and for ML or custom workflows that need clean EPICv2 matrices.

---

## Features

- Load EPICv2 beta matrices (CSV / TSV / semicolon-delimited; delimiter inferred)
- Ship a compact Peters et al. EPICv2 manifest as Parquet
- Optional QC exclusions via manifest flags (mismatch / missing position, cross-hybridisation, sex / mitochondrial chromosomes)
- Parse EPICv2 `cg…_…` probe identifiers into site, design, and replicate components
- Summarise probe-type and chromosome composition of a betas file

---

## Installation

### From GitHub

```bash
pip install git+https://github.com/MarStreicher/EpicV2IO.git
```

### Conda environment

```bash
conda env create -f environment.yml
conda activate env_epicv2io
```

### Requirements

- Python ≥ 3.9
- pandas ≥ 1.5
- pyarrow ≥ 10
- openpyxl

Dependencies are installed automatically with `pip`.

---

## Data and Annotation Resources

EpicV2IO uses a compact Parquet representation of the supplementary EPICv2 manifest published by Peters et al. Only the columns required for duplicate-probe filtering and annotation are retained.

The package distributes a preprocessed Parquet file to reduce installation size and improve loading performance.

| Column               | Role                                                                 |
| -------------------- | -------------------------------------------------------------------- |
| `IlmnID`           | Probe identifier (join key with betas index)                         |
| `CHR`              | Chromosome (`chr1` … `chr22`, `chrX`, `chrY`, `chrM`, …) |
| `MismatchPos`      | Illumina vs SeSAMe location discrepancy flag                         |
| `MissingPos`       | Missing / unresolved genomic position flag                           |
| `CH_BLAT`          | Cross-hybridisation evidence (BLAT)                                  |
| `CH_WGBS_evidence` | Cross-hybridisation evidence (WGBS)                                  |

Source in AnnotationHub: **AH116484** (Peters EPICv2 manifest). The rebuild script is `scripts/build_peters_manifest_parquet.R`.

### Purpose

Provide a small, versioned annotation table that travels with the package so loaders and downstream tools do not need a separate Bioconductor / AnnotationHub download at runtime.

### Example

```python
from epicv2io.manifest import load_peters_manifest, MANIFEST_COLUMNS

manifest = load_peters_manifest()
print(manifest.shape)
print(list(manifest.columns))
print(manifest.head(3))
```

Example result:

```text
(937690, 6)
['IlmnID', 'CHR', 'MismatchPos', 'MissingPos', 'CH_BLAT', 'CH_WGBS_evidence']

            IlmnID   CHR MismatchPos MissingPos CH_BLAT CH_WGBS_evidence
0  cg00000029_TC11  chr16           N          N       N                N
1  cg00000109_TC21   chr3           N          N       N                N
2  cg00000155_BC11  chr12           N          N       N                N
```

Subset columns if you only need chromosome mapping:

```python
chr_map = load_peters_manifest(columns=["IlmnID", "CHR"])
```

---

## Quick Start — `BetasLoader`

### Purpose

Read a beta-value matrix, keep **cg** probes only, drop annotation / control columns, and optionally remove probes that match enabled Peters QC flags. Returns either a probes × samples frame or a samples × probes (wide) frame for ML-style workflows.

**Expected input**

- Rows: probe IDs (`IlmnID`), including at least some IDs starting with `cg`
- Columns: samples (plus optional annotation columns such as `CHR`, which are dropped)
- Values: beta values in \[0, 1\]
- Delimiter: comma, tab, or semicolon (auto-detected)

### Example

```python
from pathlib import Path
from epicv2io import BetasLoader

loader = BetasLoader(Path("betas.csv"))

# probes × samples (cg only; no QC flags)
betas = loader.load_data()
print(betas.shape)
print(betas.iloc[:3, :3])

# samples × probes
betas_wide = loader.load_data_wide(
    exclude_mismatch_pos=True,
    exclude_missing_pos=True,
    exclude_ch_blat=True,
    exclude_ch_wgbs_evidence=True,
    exclude_sex_chromosomes=True,
    exclude_mitochondrial=True,
)
print(betas_wide.shape)
print(betas_wide.iloc[:2, :3])
```

Console summary printed during a filtered load (counts are dataset-dependent):

```text
Candidate cg CpGs (betas ∩ manifest, cg rows): 850,000
  exclude_mismatch_pos: 1,234 CpGs match this rule (among candidates)
  exclude_missing_pos: 567 CpGs match this rule (among candidates)
  exclude_ch_blat: 8,901 CpGs match this rule (among candidates)
  exclude_ch_wgbs_evidence: 2,345 CpGs match this rule (among candidates)
  exclude_sex_chromosomes: 19,000 CpGs match this rule (among candidates)
  exclude_mitochondrial: 12 CpGs match this rule (among candidates)
  Union (unique CpGs removed by any enabled rule): 28,500
  Remaining after union: 821,500
(821500, 48)
```

Example matrix result (`load_data`):

```text
                 Sample_A  Sample_B  Sample_C
IlmnID
cg00000029_TC11      0.12      0.15      0.11
cg00000109_TC21      0.88      0.91      0.87
cg00000155_BC11      0.45      0.44      0.48
```

Example wide result (`load_data_wide`):

```text
IlmnID     cg00000029_TC11  cg00000109_TC21  cg00000155_BC11
sample_id
Sample_A              0.12             0.88             0.45
Sample_B              0.15             0.91             0.44
```

After loading, exclusion statistics are available on the loader:

```python
print(loader.load_candidate_cg_count)
print(loader.load_exclusion_counts_by_rule)
print(loader.load_exclusion_union_count)
```

Manifest cg probes absent from the betas file:

```python
missing = loader.get_manifest_cpgs_missing_in_betas()
print(len(missing), missing[:5])
```

### Exclusion flags

All default to `False`. Truthy manifest values are `Y`, `TRUE`, and `1` (case-insensitive).

| Argument                     | Manifest column / rule         |
| ---------------------------- | ------------------------------ |
| `exclude_mismatch_pos`     | `MismatchPos`                |
| `exclude_missing_pos`      | `MissingPos`                 |
| `exclude_ch_blat`          | `CH_BLAT`                    |
| `exclude_ch_wgbs_evidence` | `CH_WGBS_evidence`           |
| `exclude_sex_chromosomes`  | `CHR` in `chrX` / `chrY` |
| `exclude_mitochondrial`    | `CHR` == `chrM`            |

When any exclusion flag is enabled, at least one probe ID must overlap the shipped manifest (`IlmnID`); otherwise a `ValueError` is raised.

---

## `CgProbeId`

### Purpose

Parse EPICv2 cg probe names of the form `cgDDDD_XXXX` into components used for duplicate grouping (same site / same design / replicate suffix).

Pattern: `^cg(\d+)_([TB][CO][12])(\d+)$`

| Component        | Meaning                                   | Example (`cg00000029_TC11`) |
| ---------------- | ----------------------------------------- | ----------------------------- |
| `site_id`      | CpG site                                  | `cg00000029`                |
| `design_type`  | Design block`[TB][CO][12]`              | `TC1`                       |
| `design_id`    | Site + design (exact-duplicate group key) | `cg00000029_TC1`            |
| `replicate_id` | Numeric replicate suffix                  | `1`                         |

### Example

```python
from epicv2io import CgProbeId

probe = "cg00000029_TC11"
print(CgProbeId.parse(probe))
print(CgProbeId.site_id(probe))
print(CgProbeId.design_id(probe))
print(CgProbeId.design_type(probe))
print(CgProbeId.replicate_id(probe))
```

Example result:

```text
('cg00000029', 'cg00000029_TC1', 'TC1', '1')
cg00000029
cg00000029_TC1
TC1
1
```

Non-matching IDs return `None`:

```python
print(CgProbeId.parse("rs12345678"))  # None
```

---

## `BetasAnalyzer`

### Purpose

Quick inventory of a betas file: sample/column counts, probe-type mix from ID prefixes (`cg`, `ch`, `rs`, `nv`, `control`), and chromosome distribution via the Peters manifest. Useful before a full load.

**Note:** the analyzer currently expects a **tab**-delimited betas file.

### Example

```python
from pathlib import Path
from epicv2io import BetasAnalyzer

BetasAnalyzer(Path("betas.txt")).summarise()
```

Example result:

```text
Betas file: betas.txt
Header (first 10 columns): ['IlmnID', 'Sample_A', 'Sample_B', ...]...
Number of columns (samples + index): 49
Number of sample columns: 48

Total probe rows (index size): 937,055

Probe type counts (from index prefix):
----------------------------------------
  cg            850,000  (90.71%)
  ch             50,000  ( 5.34%)
  rs             20,000  ( 2.13%)
  nv             10,000  ( 1.07%)
  control         5,000  ( 0.53%)
  other           2,055  ( 0.22%)
----------------------------------------
  total         937,055  (100.00%)

Probe chr counts (from manifest):
----------------------------------------
  chr1           75,000  ( 8.00%)
  ...
  undefined      12,000  ( 1.28%)
----------------------------------------
  total         937,055  (100.00%)

Distribution of autosomal vs sex chromosomes (from manifest):
----------------------------------------
  sex            25,000  ( 2.67%)
  autosomal     912,055  (97.33%)
----------------------------------------
  total         937,055  (100.00%)
```

---

## Citation

If you use EpicV2IO in published work, please cite:

> Streicher M. EpicV2IO: I/O library for reading Illumina MethylationEPIC v2 BeadChip array files. GitHub repository. 2026.

Please also cite the EPICv2 annotation resource:

> Peters TJ, Buckberry S, Hogg K, et al. Characterisation and recommendations for analysis of the Illumina MethylationEPIC v2 BeadChip. *Clinical Epigenetics*. 2024.

---

## License

MIT License.
