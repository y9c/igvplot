"""BigWig coverage source (optional).

Primary path uses ``pyBigWig`` (fast, in-process, indexed random access). When
``pyBigWig`` is unavailable it falls back to the UCSC ``bigWigToBedGraph``
binary via subprocess and slices the region. If neither is present, a helpful
error is raised.

Install the preferred backend with::

    pip install pybigwig          # microtask
    # or, via conda/bioconda:
    conda install -c bioconda pybigwig

The ucsc ``bigWigToBedGraph`` binary is also available in bioconda:
``ucsc-bigwigtobedgraph``.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from os import fspath

import numpy as np

from .region import Region

__all__ = ["BigWigUnavailableError", "read_bigwig_coverage", "coverage_from_bedgraph"]


class BigWigUnavailableError(RuntimeError):
    """Raised when no BigWig backend (pyBigWig or bigWigToBedGraph) is found."""


def coverage_from_bedgraph(
    bedgraph_text: str, region: Region, n: int, default: float = 0.0
) -> np.ndarray:
    """Parse bedGraph text (or an iterable of intervals) into a per-base array
    aligned to ``region`` (0-based, half-open).

    Each line is ``chrom<TAB>start<TAB>end<TAB>value`` (0-based, half-open
    intervals in bedGraph). Non-overlapping / out-of-region intervals are
    ignored; positions not covered get ``default``.
    """
    depths = np.full(n, default, dtype=np.float64)
    lines = bedgraph_text.splitlines() if isinstance(bedgraph_text, str) else bedgraph_text
    for line in lines:
        if isinstance(line, (list, tuple)):
            # iterable of intervals: (chrom, start, end, value)
            if len(line) < 4:
                continue
            chrom, s, e, val = line[0], int(line[1]), int(line[2]), float(line[3])
        else:
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                chrom, s, e, val = parts[0], int(parts[1]), int(parts[2]), float(parts[3])
            except (ValueError, IndexError):
                continue
        if chrom != region.chrom:
            continue
        s = max(s, region.start)
        e = min(e, region.end)
        if e > s:
            depths[s - region.start : e - region.start] = val
    return depths


def _read_bigwig_with_pybigwig(path: str, region: Region) -> np.ndarray:
    import pyBigWig  # noqa: F401  (imported lazily)

    n = region.length
    bw = pyBigWig.open(fspath(path))
    try:
        chrom = region.chrom
        if chrom not in bw.chroms():
            raise KeyError(chrom)
        vals = bw.values(chrom, region.start, region.end, numpy=True)
    finally:
        bw.close()
    if vals is None:
        return np.zeros(n, dtype=np.float64)
    vals = np.asarray(vals, dtype=np.float64)
    if vals.shape[0] != n:
        out = np.zeros(n, dtype=np.float64)
        out[: min(n, vals.shape[0])] = vals[: min(n, vals.shape[0])]
        return out
    return np.nan_to_num(vals)


def _read_bigwig_with_ucsc_tools(path: str, region: Region) -> np.ndarray:
    exe = shutil.which("bigWigToBedGraph")
    if not exe:
        raise BigWigUnavailableError(
            "BigWig coverage requires either 'pyBigWig' (pip install pybigwig) "
            "or the UCSC 'bigWigToBedGraph' binary (conda install -c bioconda "
            "ucsc-bigwigtobedgraph)."
        )
    with tempfile.NamedTemporaryFile(suffix=".bedGraph", delete=False) as tmp:
        out_path = tmp.name
    try:
        proc = subprocess.run(
            [exe, path, out_path], capture_output=True, text=True, check=False
        )
        if proc.returncode != 0:
            raise BigWigUnavailableError(
                f"bigWigToBedGraph failed: {proc.stderr.strip()}"
            )
        with open(out_path) as fh:
            data = fh.read()
    finally:
        try:
            os.unlink(out_path)
        except OSError:
            pass
    return coverage_from_bedgraph(data, region, region.length)


def read_bigwig_coverage(path: str, region) -> np.ndarray:
    """Read per-base coverage from a BigWig over ``region``.

    Returns a numpy array of length ``region.length`` aligned to the 0-based,
    half-open region. Falls back from pyBigWig to the UCSC tools.
    """
    region = Region.from_any(region)
    path = fspath(path)
    try:
        return _read_bigwig_with_pybigwig(path, region)
    except ImportError:
        pass
    return _read_bigwig_with_ucsc_tools(path, region)
