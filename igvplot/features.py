"""Load gene / feature annotation tracks, rendered by ``dna_features_viewer``.

``dna_features_viewer`` provides clean, publication-style gene/feature arrows
with automatic label layout. We keep the annotation in **global** reference
coordinates (0-based), crop it to the requested region, and plot it on a
matplotlib axis shared with the read/coverage axes so everything lines up.
"""
from __future__ import annotations

import os
from os import fspath
from typing import List, Optional, Union

from dna_features_viewer import BiopythonTranslator, GraphicFeature, GraphicRecord

from .region import Region

__all__ = ["load_features", "crop_to_region", "FeaturesError"]


class FeaturesError(Exception):
    """Raised when an annotation source cannot be loaded."""


def load_features(
    source: Union[str, "object", GraphicRecord, List[GraphicFeature]],
    region: Optional[Region] = None,
) -> GraphicRecord:
    """Load an annotation source and return a dna_features_viewer
    :class:`GraphicRecord` in **global** coordinates.

    Parameters
    ----------
    source:
        - a ``.gb`` / ``.genbank`` file path,
        - a Biopython ``SeqRecord``,
        - an existing :class:`GraphicRecord` (used as-is),
        - a list of :class:`GraphicFeature` (assumed to be in global
          coordinates), or
        - a ``.gff`` / ``.gff3`` / ``.gtf`` path (requires the optional
          ``bcbio-gff`` dependency).
    region:
        If given, the record is cropped to this (0-based) interval first.
    """
    if isinstance(source, os.PathLike):
        source = fspath(source)
    if isinstance(source, GraphicRecord):
        record = source
    elif isinstance(source, (str, bytes)):
        record = _load_from_path(fspath(source))
    elif isinstance(source, list):
        # assume global coordinates; sequence_length irrelevant for a region
        record = GraphicRecord(
            first_index=0,
            sequence_length=int(1e9),
            features=source,
        )
    elif _looks_like_biopython_record(source):
        record = BiopythonTranslator().translate_record(source)
    else:
        raise FeaturesError(
            f"Unsupported annotation source of type {type(source).__name__}; "
            "expected a .gb/.gff path, a SeqRecord, a GraphicRecord or a list "
            "of GraphicFeature."
        )

    if region is not None:
        record = crop_to_region(record, region)
    return record


def _load_from_path(path: str) -> GraphicRecord:
    lower = path.lower()
    if lower.endswith((".gb", ".gbk", ".genbank")):
        return BiopythonTranslator().translate_record(path)
    if lower.endswith((".gff", ".gff3", ".gtf")):
        return _load_from_gff(path)
    # Unknown extension -> try as genbank, then as GFF.
    try:
        return BiopythonTranslator().translate_record(path)
    except Exception:
        return _load_from_gff(path)


def _load_from_gff(path: str) -> GraphicRecord:
    try:
        from BCBio import GFF
    except ImportError as exc:  # pragma: no cover
        raise FeaturesError(
            "Reading GFF/GTF requires the 'bcbio-gff' package: "
            "pip install 'igvplot[gff]'"
        ) from exc

    with open(path) as handle:
        try:
            record = next(GFF.parse(handle))
        except StopIteration:
            raise FeaturesError(f"No records found in GFF file {path!r}")
    return BiopythonTranslator().translate_record(record)


def _looks_like_biopython_record(source) -> bool:
    return hasattr(source, "features") and hasattr(source, "seq")


def crop_to_region(record: GraphicRecord, region: Region) -> GraphicRecord:
    """Crop a GraphicRecord to a 0-based interval, tolerating out-of-range.

    DnaFeaturesViewer's ``crop`` keeps features in global coordinates and sets
    ``first_index`` to the window start, so the plotted arrows align exactly
    with our read/coverage axes.
    """
    first_index = record.first_index or 0
    last_index = first_index + record.sequence_length
    if region.start >= last_index or region.end <= first_index:
        return GraphicRecord(
            first_index=region.start,
            sequence_length=region.length,
            features=[],
        )
    start = max(first_index, min(region.start, last_index))
    end = max(first_index, min(region.end, last_index))
    try:
        return record.crop((start, end))
    except Exception:
        return GraphicRecord(first_index=start, sequence_length=end - start, features=[])
