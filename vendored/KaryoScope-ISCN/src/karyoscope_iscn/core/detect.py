"""Structural-rearrangement detection from scaffold segments.

A sample scaffold that maps to a single reference chromosome is normal; one that
maps to two or more is a structural rearrangement (translocation / derivative
chromosome). This module extracts breakpoints (with the cytoband on each side)
and the per-scaffold structural events, then consolidates reciprocal pairs.

Centromere-dependent classification (dicentric, Robertsonian) is layered on top
in a later module once centromere calls are available; what lives here is the
deterministic, cytoband-only structural skeleton.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import pairwise

from karyoscope_iscn.core.segments import Segment

#: Event kinds emitted here. Centromere-derived kinds are added downstream.
TRANSLOCATION = "translocation"


@dataclass(frozen=True)
class Breakpoint:
    """A junction on a scaffold between two reference chromosomes.

    ``pos`` is the scaffold coordinate of the junction; ``up_band`` /
    ``down_band`` are the cytobands flanking it (the basis for ISCN notation).
    """

    scaffold: str
    pos: int
    up_chrom: str
    down_chrom: str
    up_band: str | None
    down_band: str | None


@dataclass
class StructuralEvent:
    """A structural rearrangement detected on one scaffold."""

    kind: str
    scaffold: str
    chromosomes: list[str]
    breakpoints: list[Breakpoint] = field(default_factory=list)

    def chromosome_set(self) -> frozenset[str]:
        return frozenset(self.chromosomes)


def find_breakpoints(scaffold: str, segments: list[Segment]) -> list[Breakpoint]:
    """Return breakpoints at each junction between differently-mapped segments.

    ``segments`` must be in scaffold order. The flanking bands are the band
    nearest the junction on each side (the last band of the upstream segment and
    the first band of the downstream segment), or ``None`` when that side
    carried only bare-chromosome/filler labels.
    """
    breakpoints: list[Breakpoint] = []
    for up, down in pairwise(segments):
        if up.chrom == down.chrom:
            continue
        breakpoints.append(
            Breakpoint(
                scaffold=scaffold,
                pos=down.start,
                up_chrom=up.chrom,
                down_chrom=down.chrom,
                up_band=up.bands[-1] if up.bands else None,
                down_band=down.bands[0] if down.bands else None,
            )
        )
    return breakpoints


def distinct_chromosomes(segments: list[Segment]) -> list[str]:
    """Return the distinct reference chromosomes a scaffold maps to, in order."""
    seen: list[str] = []
    for s in segments:
        if s.chrom not in seen:
            seen.append(s.chrom)
    return seen


def detect_structural_event(
    scaffold: str,
    segments: list[Segment],
) -> StructuralEvent | None:
    """Classify a scaffold's structural state from its segments.

    Returns ``None`` for a normal scaffold (zero or one chromosome). A scaffold
    spanning two or more chromosomes yields a :data:`TRANSLOCATION` event with
    its breakpoints. (Dicentric / Robertsonian refinement happens downstream,
    once centromere calls are available.)
    """
    chroms = distinct_chromosomes(segments)
    if len(chroms) < 2:
        return None
    return StructuralEvent(
        kind=TRANSLOCATION,
        scaffold=scaffold,
        chromosomes=chroms,
        breakpoints=find_breakpoints(scaffold, segments),
    )


@dataclass
class ReciprocalGroup:
    """Two scaffolds whose structural events involve the same chromosome pair."""

    chromosomes: frozenset[str]
    events: list[StructuralEvent]


def consolidate_reciprocal(events: list[StructuralEvent]) -> list[ReciprocalGroup]:
    """Group events that share an identical chromosome set into reciprocal groups.

    Two derivative scaffolds carrying the same chromosome pair are the two
    products of one reciprocal translocation. Only groups with two or more
    events are returned; singletons are not reciprocal.
    """
    by_set: dict[frozenset[str], list[StructuralEvent]] = {}
    for ev in events:
        by_set.setdefault(ev.chromosome_set(), []).append(ev)
    return [
        ReciprocalGroup(chromosomes=chrom_set, events=evs)
        for chrom_set, evs in by_set.items()
        if len(evs) >= 2
    ]
