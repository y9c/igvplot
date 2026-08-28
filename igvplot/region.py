"""Genomic interval parsing and representation.

Coordinates are stored in the same convention as pysam / BAM files:
**0-based, half-open** (``start`` inclusive, ``end`` exclusive). Human-readable
region strings (``chr1:1,000-2,000``) use 1-based, inclusive coordinates and
are converted on parse.
"""
from __future__ import annotations

import re
from typing import Union

__all__ = ["Region"]


class Region:
    """A genomic interval on a single chromosome/contig.

    Parameters
    ----------
    chrom: str
        Chromosome / contig name.
    start: int
        0-based inclusive start position.
    end: int
        0-based exclusive end position.
    """

    __slots__ = ("chrom", "start", "end")

    def __init__(self, chrom: str, start: int, end: int):
        start = int(start)
        end = int(end)
        if not chrom:
            raise ValueError("chrom must be a non-empty string")
        if start < 0 or end < 0:
            raise ValueError("coordinates must be non-negative")
        if end < start:
            raise ValueError(f"end ({end}) is smaller than start ({start})")
        self.chrom = str(chrom)
        self.start = start
        self.end = end

    # ------------------------------------------------------------------ #
    # properties
    # ------------------------------------------------------------------ #
    @property
    def length(self) -> int:
        """Number of reference bases in the interval."""
        return self.end - self.start

    def __repr__(self) -> str:
        return f"Region({self.chrom!r}, {self.start}, {self.end})"

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, Region)
            and (self.chrom, self.start, self.end) == (other.chrom, other.start, other.end)
        )

    def __hash__(self) -> int:
        return hash((self.chrom, self.start, self.end))

    # ------------------------------------------------------------------ #
    # constructors
    # ------------------------------------------------------------------ #
    @classmethod
    def from_string(cls, text: str) -> "Region":
        """Parse a human-readable region.

        Accepts strings such as ``chr1:1,000-2,000`` or ``11:1000-2000``
        (1-based, inclusive). Commas and whitespace around coordinates are
        ignored.
        """
        cleaned = re.sub(r"[\s_,]", "", text).strip()
        match = re.match(r"^(?P<chrom>[A-Za-z0-9_.]+):(?P<start>\d+)-(?P<end>\d+)$", cleaned)
        if not match:
            raise ValueError(
                f"Cannot parse region string {text!r}; expected 'chr1:1000-2000'"
            )
        start1, end1 = int(match.group("start")), int(match.group("end"))
        if end1 < start1:
            raise ValueError(f"Region {text!r} has end < start")
        return cls(match.group("chrom"), start1 - 1, end1)

    @classmethod
    def from_any(cls, region: Union["Region", str, tuple, list]) -> "Region":
        """Coerce a string, a ``(chrom, start, end)`` tuple/list, or a
        :class:`Region` into a :class:`Region`. Tuples/lists are 0-based,
        half-open (pysam convention)."""
        if isinstance(region, cls):
            return region
        if isinstance(region, str):
            return cls.from_string(region)
        if isinstance(region, (tuple, list)):
            if len(region) != 3:
                raise ValueError(
                    "a tuple/list region must have exactly 3 items: (chrom, start, end)"
                )
            chrom, start, end = region
            return cls(chrom, int(start), int(end))
        raise TypeError(f"Cannot build Region from {type(region).__name__}")

    @classmethod
    def centered(cls, chrom: str, pos: int, window: int = 100) -> "Region":
        """Build a region centred on a 1-based ``pos`` with ``window`` bases on
        each side (the smvplot 'chr:center' convention)."""
        pos, window = int(pos), int(window)
        center0 = pos - 1  # 0-based
        start = max(0, center0 - window)
        return cls(chrom, start, center0 + window + 1)

    # ------------------------------------------------------------------ #
    # geometry
    # ------------------------------------------------------------------ #
    def overlaps(self, start: int, end: int) -> bool:
        """True if the half-open interval [start, end) overlaps this region."""
        return not (end <= self.start or start >= self.end)
