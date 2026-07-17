from typing import Optional

import pytest

from epicv2io import CgProbeId


@pytest.mark.parametrize(
    ("probe_id", "expected"),
    [
        (
            "cg00000029_TC11",
            ("cg00000029", "cg00000029_TC1", "TC1", "1"),
        ),
        (
            " cg123_BO220 ",
            ("cg123", "cg123_BO2", "BO2", "20"),
        ),
        (
            "cg9_TO110",
            ("cg9", "cg9_TO1", "TO1", "10"),
        ),
    ],
)
def test_parse_valid_probe_ids(
    probe_id: str, expected: tuple[str, str, str, str]
) -> None:
    assert CgProbeId.parse(probe_id) == expected


@pytest.mark.parametrize(
    "probe_id",
    [
        "CG00000029_TC11",
        "cg00000029_XC11",
        "cg00000029_TC31",
        "cg00000029_TC1",
        "cg_TC11",
        "cg00000029_TC11_extra",
        "ch00000029_TC11",
        "",
    ],
)
def test_parse_rejects_invalid_probe_ids(probe_id: str) -> None:
    assert CgProbeId.parse(probe_id) is None


def test_component_accessors_return_parsed_values() -> None:
    probe_id = "cg00000029_TC110"

    assert CgProbeId.site_id(probe_id) == "cg00000029"
    assert CgProbeId.design_id(probe_id) == "cg00000029_TC1"
    assert CgProbeId.design_type(probe_id) == "TC1"
    assert CgProbeId.replicate_id(probe_id) == "10"


@pytest.mark.parametrize(
    "accessor",
    [
        CgProbeId.site_id,
        CgProbeId.design_id,
        CgProbeId.design_type,
        CgProbeId.replicate_id,
    ],
)
def test_component_accessors_return_none_for_invalid_probe(accessor) -> None:
    assert accessor("not-a-probe") is None


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("cg000123_TC11", "000123"),
        (" CG42-anything", "42"),
        ("cg7", "7"),
        ("ch123", None),
        ("prefix-cg123", None),
        (123, None),
    ],
)
def test_site_digits_is_case_insensitive_and_accepts_stringable_values(
    value: str, expected: Optional[str]
) -> None:
    assert CgProbeId.site_digits(value) == expected


def test_parse_results_are_cached() -> None:
    CgProbeId.parse.cache_clear()

    first = CgProbeId.parse("cg00000029_TC11")
    second = CgProbeId.parse("cg00000029_TC11")

    assert first is second
    assert CgProbeId.parse.cache_info().hits == 1
