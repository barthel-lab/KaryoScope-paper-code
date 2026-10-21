"""Segment a sample scaffold into contiguous reference-chromosome runs.

Given a scaffold's cytoband-labeled intervals, group them into maximal runs that
map to a single reference chromosome. A scaffold that switches from one
chromosome to another (e.g. chr13 -> chr14) yields multiple segments -- the raw
material for translocation / Robertsonian detection. Sentinel fillers
(``categorized`` / ``novel``) and unassigned gaps extend the current segment
without breaking it; bare-chromosome labels extend their chromosome's segment
but contribute no arm/band.

All coordinates are scaffold coordinates. Band ordering and reference
coordinates are *not* used here -- this module only partitions the scaffold by
chromosome identity.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from karyoscope_iscn.core.bedio import BedInterval
from karyoscope_iscn.core.cytoband import parse_band
from karyoscope_iscn.reference import Reference


@dataclass
class Segment:
    """A maximal run of one scaffold that maps to a single reference chromosome."""

    scaffold: str
    chrom: str
    start: int
    end: int
    covered_bp: int = 0
    arms: set[str] = field(default_factory=set)
    bands: list[str] = field(default_factory=list)

    @property
    def span(self) -> int:
        return self.end - self.start


def segment_scaffold(
    scaffold: str,
    intervals: list[BedInterval],
    reference: Reference,
) -> list[Segment]:
    """Partition one scaffold's intervals into per-chromosome segments.

    ``intervals`` should be sorted by start (see
    :func:`karyoscope_iscn.core.bedio.group_by_scaffold`). Intervals whose label
    does not resolve to a reference chromosome (``categorized``/``novel`` or an
    unknown token) extend the current segment as filler.
    """
    from karyoscope_iscn.core.cytoband import chrom_for_token

    segments: list[Segment] = []
    current: Segment | None = None

    for iv in intervals:
        band = parse_band(iv.name)
        chrom = chrom_for_token(band.chrom_token, reference)

        if chrom is None:
            # Filler (gap / centromere / telomere / unresolved): extend current.
            if current is not None:
                current.end = max(current.end, iv.end)
            continue

        if current is not None and current.chrom == chrom:
            current.end = max(current.end, iv.end)
        else:
            current = Segment(scaffold=scaffold, chrom=chrom, start=iv.start, end=iv.end)
            segments.append(current)

        current.covered_bp += iv.end - iv.start
        if band.arm is not None:
            current.arms.add(band.arm)
            current.bands.append(band.raw)

    return segments


def merge_small_segments(
    segments: list[Segment],
    *,
    min_bp: int,
) -> list[Segment]:
    """Drop segments whose covered bp is below ``min_bp``.

    Tiny segments are typically noise (a handful of cross-mapping k-mers). They
    are removed rather than merged, since a real chromosomal segment carries
    substantial sequence. Returns a new list; order is preserved.
    """
    return [s for s in segments if s.covered_bp >= min_bp]


def is_acrocentric_short_arm(segment: Segment, reference: Reference) -> bool:
    """True if a segment is *only* an acrocentric chromosome's short (p) arm.

    The acrocentric short arms (the rDNA/satellite p-arms of chr13/14/15/21/22) are
    nearly homologous and recombine across the human population. A k-mer index built
    from the single CHM13 individual therefore maps short-arm sequence from one
    acrocentric onto another, so a scaffold segment that is *entirely* an acrocentric
    p-arm is shared acrocentric sequence -- expected, not a rearrangement partner.

    Deliberately strict, because this filter's output is also a chromosome's
    copy-number evidence: the chromosome must be acrocentric *and* the segment must
    carry p-arm bands only. Anything reaching the q-arm is kept. For deciding
    rearrangement *partners*, use the broader
    :func:`is_acrocentric_shared_region`, which also covers the proximal q11 band.
    """
    info = reference.chromosomes.get(segment.chrom)
    return bool(info and info.acrocentric and segment.arms == {"p"})


def is_acrocentric_shared_region(segment: Segment, reference: Reference) -> bool:
    """True if a segment is acrocentric sequence shared across the acrocentrics.

    Broader than :func:`is_acrocentric_short_arm`: the shared region does not stop
    at the centromere, so the proximal ``q11`` band -- pericentromeric satellite and
    segmental duplication -- counts too. HG02841's ``der(14)t(14;22)(?;q11.21)``
    rests entirely on chr22's q11.21 band, which the p-only test let through.

    **Use this for rearrangement-partner decisions only, never to filter the
    segment list that feeds copy number.** Those segments are a chromosome's
    evidence that it is present at all: dropping centromere-only and short-arm
    chr22 contigs cost GM03786 its chr22 copy number and produced a spurious
    ``-22``. A partner reaching beyond ``q11`` (any real Robertsonian, which
    carries a whole long arm) is not shared sequence and is kept.
    """
    info = reference.chromosomes.get(segment.chrom)
    if not (info and info.acrocentric):
        return False
    if segment.arms == {"p"}:
        return True
    if not segment.bands:
        return False
    return all(_is_short_arm_or_proximal(b) for b in segment.bands)


def _is_short_arm_or_proximal(label: str) -> bool:
    """True for an acrocentric p-arm band or the pericentromeric ``q11`` band."""
    band = parse_band(label)
    if band.arm == "p":
        return True
    return band.arm == "q" and (band.sub or "").startswith("11")


def drop_acrocentric_short_arms(
    segments: list[Segment],
    reference: Reference,
) -> list[Segment]:
    """Drop segments that are only an acrocentric short arm (see
    :func:`is_acrocentric_short_arm`).

    Removes shared-acrocentric cross-mapping so it is not mistaken for a
    translocation partner. Returns a new list; order is preserved. Callers should
    re-run :func:`merge_adjacent_segments` afterwards in case removing a spurious
    segment leaves two same-chromosome pieces adjacent.
    """
    return [s for s in segments if not is_acrocentric_short_arm(s, reference)]


def merge_adjacent_segments(segments: list[Segment]) -> list[Segment]:
    """Coalesce consecutive segments that map to the same chromosome.

    After :func:`merge_small_segments` removes a noise segment sitting between
    two pieces of the same chromosome, those pieces should read as one segment.
    Merging unions their span, covered bp, arms, and bands (in order). Returns a
    new list.
    """
    merged: list[Segment] = []
    for seg in segments:
        if merged and merged[-1].chrom == seg.chrom:
            prev = merged[-1]
            prev.end = max(prev.end, seg.end)
            prev.start = min(prev.start, seg.start)
            prev.covered_bp += seg.covered_bp
            prev.arms |= seg.arms
            prev.bands.extend(seg.bands)
        else:
            merged.append(
                Segment(
                    scaffold=seg.scaffold,
                    chrom=seg.chrom,
                    start=seg.start,
                    end=seg.end,
                    covered_bp=seg.covered_bp,
                    arms=set(seg.arms),
                    bands=list(seg.bands),
                )
            )
    return merged
