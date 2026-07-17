from __future__ import annotations

import re
from functools import lru_cache
from typing import Optional, Tuple

# (site_id, design_id, design_type, replicate_id) from ``parse``.
CgProbeParseResult = Tuple[str, str, str, str]

_CG_PROBE_PATTERN = re.compile(r"^cg(\d+)_([TB][CO][12])(\d+)$")
_CG_SITE_PATTERN = re.compile(r"^cg(\d+)", re.IGNORECASE)


class CgProbeId:
    """
    Parse cg probe identifiers of the form cgDDDD_XXXX.
    """

    @staticmethod
    @lru_cache(maxsize=None)
    def parse(probe_id: str) -> Optional[CgProbeParseResult]:
        """
        Parse ``^cg(\\d+)_([TB][CO][12])(\\d+)$`` into a 4-tuple.

        - site_id: digits after ``cg`` (same target CpG).
        - design_id: ``site_id`` + ``design_type``; stable key for exact-duplicate grouping.
        - design_type: the three-character block ``[TB][CO][12]`` after the underscore.
        - replicate_id: numeric replicate suffix after ``design_type`` (one or more digits;
          e.g. ``1``, ``10``, ``20``).

        Returns ``None`` if the string does not match the cg probe pattern.
        """
        m = _CG_PROBE_PATTERN.match(probe_id.strip())
        if not m:
            return None

        digits = m.group(1)
        site_id = f"cg{digits}"
        design_type = m.group(2)
        replicate_id = m.group(3)
        design_id = f"{site_id}_{design_type}"

        return (site_id, design_id, design_type, replicate_id)

    @staticmethod
    def site_id(probe_id: str) -> Optional[str]:
        """Return the site part (DDDD) of cgDDDD_XXXX, or ``None`` if not a cg probe."""
        parsed = CgProbeId.parse(probe_id)
        return parsed[0] if parsed else None

    @staticmethod
    def design_id(probe_id: str) -> Optional[str]:
        """Return ``site_id + design_type``; same value ⇒ same exact-duplicate group."""
        parsed = CgProbeId.parse(probe_id)
        return parsed[1] if parsed else None

    @staticmethod
    def design_type(probe_id: str) -> Optional[str]:
        """Return the ``[TB][CO][12]`` trio after the underscore, or ``None`` if not a cg probe."""
        parsed = CgProbeId.parse(probe_id)
        return parsed[2] if parsed else None

    @staticmethod
    def replicate_id(probe_id: str) -> Optional[str]:
        """Return the numeric replicate suffix (digits after ``design_type``)."""
        parsed = CgProbeId.parse(probe_id)
        return parsed[3] if parsed else None

    @staticmethod
    def site_digits(s: str) -> Optional[str]:
        """Return the digits after ``cg`` (case-insensitive), or ``None`` if not a cg probe."""
        m = _CG_SITE_PATTERN.match(str(s).strip())
        return m.group(1) if m else None
