from pathlib import Path

from epicv2io.io import infer_separator


def test_infer_separator_honours_sample_byte_limit(tmp_path: Path) -> None:
    path = tmp_path / "betas.txt"
    path.write_text("IlmnID;A\ncg1;0.1\nignored,comma\n")

    assert infer_separator(path, sample_bytes=17) == ";"
