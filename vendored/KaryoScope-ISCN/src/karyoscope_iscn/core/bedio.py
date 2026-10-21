"""Minimal BED reader for cytoband / region annotation files.

Sample annotation BEDs have the layout ``col1  start  end  label`` where
``col1`` is the *scaffold/contig* name of the sample assembly and ``label`` is
the matched reference feature (a cytoband such as ``13q21.1`` or a region
feature). The reference ideogram uses the same layout with ``col1`` being the
reference chromosome.

This reader is deliberately dependency-free (it does not import karyoplot) so
the core logic and its tests run without the plotting stack. It transparently
reads ``.gz`` files.
"""

from __future__ import annotations

import gzip
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from karyoscope_iscn.exceptions import InputError


@dataclass(frozen=True)
class BedInterval:
    """One BED record: ``chrom`` is col1 (scaffold or reference chrom)."""

    chrom: str
    start: int
    end: int
    name: str


@contextmanager
def open_text(path: Path) -> Iterator[TextIO]:
    """Open a text file, transparently decompressing ``.gz``."""
    # This generator IS the context manager; it closes `handle` in `finally`.
    if path.suffix == ".gz":
        handle: TextIO = gzip.open(path, "rt")  # noqa: SIM115
    else:
        handle = path.open("rt")
    try:
        yield handle
    finally:
        handle.close()


def read_bed(path: Path, *, min_columns: int = 4) -> Iterator[BedInterval]:
    """Yield :class:`BedInterval` records from a (optionally gzipped) BED file.

    Raises :class:`InputError` on a missing file or a row with too few columns
    or non-integer coordinates.
    """
    if not path.is_file():
        raise InputError(f"BED file not found: {path}")
    with open_text(path) as handle:
        for i, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split("\t")
            if len(parts) < min_columns:
                raise InputError(
                    f"{path}:{i}: expected >={min_columns} columns, got {len(parts)}: {line!r}"
                )
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError as e:
                raise InputError(f"{path}:{i}: non-integer coordinates: {line!r}") from e
            yield BedInterval(chrom=parts[0], start=start, end=end, name=parts[3])


def load_bed(path: Path, *, min_columns: int = 4) -> list[BedInterval]:
    """Eagerly read a BED file into a list."""
    return list(read_bed(path, min_columns=min_columns))


def iter_bed_filtered(
    path: Path,
    *,
    names: set[str] | None = None,
    chroms: set[str] | None = None,
) -> Iterator[BedInterval]:
    """Yield only records whose name and/or scaffold match -- fast over huge BEDs.

    A cheap substring pre-check (``token in line``) rejects the vast majority of
    lines before the expensive split / int-parse / object construction, so
    extracting (say) the SST1 intervals on one scaffold from a ~20 M-line
    scaffolded BED costs a substring scan per line instead of 20 M ``BedInterval``
    allocations. The substring test can false-positive (a value matching in
    another column), so an exact column check confirms each survivor. Malformed
    rows are skipped (this is a tolerant fast path, not the strict ``read_bed``).
    """
    if not path.is_file():
        raise InputError(f"BED file not found: {path}")
    with open_text(path) as handle:
        for raw in handle:
            # Require a substring hit in EACH provided group (name AND chrom), so a
            # line survives the cheap pass only if it can match the exact filters.
            if names is not None and not any(n in raw for n in names):
                continue
            if chroms is not None and not any(c in raw for c in chroms):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, name = parts[0], parts[3]
            if chroms is not None and chrom not in chroms:
                continue
            if names is not None and name not in names:
                continue
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            yield BedInterval(chrom=chrom, start=start, end=end, name=name)


def group_by_scaffold(intervals: list[BedInterval]) -> dict[str, list[BedInterval]]:
    """Group intervals by ``chrom`` (scaffold), each sorted by start.

    Preserves first-appearance order of scaffolds.
    """
    groups: dict[str, list[BedInterval]] = {}
    for iv in intervals:
        groups.setdefault(iv.chrom, []).append(iv)
    for ivs in groups.values():
        ivs.sort(key=lambda iv: (iv.start, iv.end))
    return groups
