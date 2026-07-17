from pathlib import Path

import pandas as pd
import pytest

from epicv2io import BetasAnalyzer
from epicv2io import betas_analyzer


def test_init_exits_with_error_when_betas_file_is_missing(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    missing_path = tmp_path / "missing.tsv"

    with pytest.raises(SystemExit) as exc_info:
        BetasAnalyzer(missing_path)

    assert exc_info.value.code == 1
    assert capsys.readouterr().err == f"Error: betas file not found: {missing_path}\n"


def test_get_header_and_sample_count_reads_tab_delimited_header(
    tmp_path: Path,
) -> None:
    path = tmp_path / "betas.tsv"
    path.write_text("IlmnID\tSample_A\tSample_B\ncg1\t0.1\t0.2\n")

    assert BetasAnalyzer(path).get_header_and_sample_count() == (
        ["IlmnID", "Sample_A", "Sample_B"],
        3,
    )


def test_load_manifest_requests_only_needed_columns(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "betas.tsv"
    path.write_text("IlmnID\n")
    source = pd.DataFrame({"IlmnID": [1], "CHR": [2]})
    observed = {}

    def fake_load_peters_manifest(*, columns):
        observed["columns"] = columns
        return source

    monkeypatch.setattr(
        betas_analyzer, "load_peters_manifest", fake_load_peters_manifest
    )

    result = BetasAnalyzer(path)._load_manifest()

    assert observed["columns"] == ("IlmnID", "CHR")
    assert result.iloc[0].tolist() == ["1", "2"]


def test_summarise_reports_probe_and_chromosome_counts(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    path = tmp_path / "betas.tsv"
    path.write_text(
        "IlmnID\tSample_A\tSample_B\n"
        " cg1 \t0.1\t0.2\n"
        "CH2\t0.1\t0.2\n"
        "rs3\t0.1\t0.2\n"
        "nv4\t0.1\t0.2\n"
        "control_5\t0.1\t0.2\n"
        "other6\t0.1\t0.2\n"
    )
    analyzer = BetasAnalyzer(path)
    manifest = pd.DataFrame(
        {
            "IlmnID": ["CG1", "ch2", "rs3"],
            "CHR": ["chr1", "chrX", "chrY"],
        }
    )
    monkeypatch.setattr(analyzer, "_load_manifest", lambda: manifest)

    analyzer.summarise()

    output = capsys.readouterr().out
    assert f"Betas file: {path}" in output
    assert "Number of columns (samples + index): 3" in output
    assert "Number of sample columns: 2" in output
    assert "Total probe rows (index size): 6" in output
    for probe_type in ("cg", "ch", "rs", "nv", "control", "other"):
        assert f"  {probe_type:10}          1  (16.67%)" in output
    assert "  undefined           3  (50.00%)" in output
    assert "  sex                 2  (33.33%)" in output
    assert "  autosomal           4  (66.67%)" in output


def test_summarise_handles_file_with_no_probe_rows(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    path = tmp_path / "empty.tsv"
    path.write_text("IlmnID\tSample_A\n")
    analyzer = BetasAnalyzer(path)
    monkeypatch.setattr(
        analyzer,
        "_load_manifest",
        lambda: pd.DataFrame({"IlmnID": [], "CHR": []}, dtype=str),
    )

    analyzer.summarise()

    output = capsys.readouterr().out
    assert "Total probe rows (index size): 0" in output
    assert "  sex                 0  ( 0.00%)" in output
    assert "  autosomal           0  ( 0.00%)" in output
