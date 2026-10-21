"""Symmetric cumulative k-mer breakpoint refinement.

A second, finer breakpoint estimate that *composes on top of* the
content-dominance transition interval from :mod:`karyoscope_iscn.core.refine`.
Given that interval extended by a margin on each side (the analysis window), it
computes, at fine resolution, two cumulative-fraction curves from the
chromosome-identity signal:

* the **left-cumulative fraction** for one fusion partner -- the proportion of
  that partner's chromosome-specific content in the window lying at or before
  each position (rises 0 -> 1 left to right), and
* the **right-cumulative fraction** for the other partner -- the proportion
  lying at or after each position (falls 1 -> 0 left to right).

Their sum is the **combined score** (theoretical max 2 at a perfectly clean
junction). The position of the maximum is the **inferred breakpoint**, and the
**plateau** is the interval over which the score stays within ``1 - threshold``
of its maximum (the *99% plateau width* at the default ``threshold = 0.99``).

This is *not* a confidence interval -- it is the width of the score plateau, a
descriptive measure of how sharply the junction is localised.

Chromosome identity comes from a pluggable resolver shared with
:mod:`~karyoscope_iscn.core.refine` (see its ``PIVOTS``), so this runs on the
scaffolded (pre-bin) BED of whichever feature set the pivot names: cytoband
labels by default, or the ``chromosome`` feature set directly. The two modules
must always be given the SAME pivot -- refine supplies the transition interval
this analysis then sharpens, so a mismatch would sharpen one label space's
interval using another's evidence. Pure Python, no array dependency: the curves
are short (window / resolution bins).
"""

from __future__ import annotations

from dataclasses import dataclass, field

from karyoscope_iscn.core.bedio import BedInterval
from karyoscope_iscn.core.detect import Breakpoint
from karyoscope_iscn.core.refine import DEFAULT_PIVOT, ChromResolver, resolver_for
from karyoscope_iscn.reference import Reference

#: Fine resolution (bp) for the cumulative curves.
DEFAULT_RESOLUTION = 100
#: Margin (bp) added to each side of the transition interval to form the window.
DEFAULT_WINDOW_PAD = 1_000_000
#: Plateau threshold: positions scoring within ``1 - threshold`` of the maximum.
DEFAULT_PLATEAU_THRESHOLD = 0.99


@dataclass(frozen=True)
class CumulativeRefinement:
    """The cumulative-k-mer breakpoint estimate for one junction.

    ``positions`` and the curve lists are parallel (one entry per resolution
    bin) and are kept for the density-panel visualisation; the scalar fields are
    what the report records.
    """

    inferred_breakpoint: int
    plateau_start: int
    plateau_end: int
    max_score: float
    threshold: float
    window_start: int
    window_end: int
    resolution: int
    positions: list[int] = field(default_factory=list)
    left_fraction: list[float] = field(default_factory=list)
    right_fraction: list[float] = field(default_factory=list)
    combined_score: list[float] = field(default_factory=list)

    @property
    def plateau_width(self) -> int:
        return self.plateau_end - self.plateau_start


def _coverage(
    intervals: list[BedInterval],
    chrom: str,
    chrom_of: ChromResolver,
    win_start: int,
    n_bins: int,
    resolution: int,
) -> list[int]:
    """Per-bin bp coverage of intervals resolving to ``chrom`` over the window."""
    cov = [0] * n_bins
    win_end = win_start + n_bins * resolution
    for iv in intervals:
        if iv.end <= win_start or iv.start >= win_end:
            continue
        if chrom_of(iv.name) != chrom:
            continue
        s = max(0, (iv.start - win_start) // resolution)
        e = min(n_bins, -(-(iv.end - win_start) // resolution))  # ceil division
        for i in range(s, e):
            bin_start = win_start + i * resolution
            overlap = min(iv.end, bin_start + resolution) - max(iv.start, bin_start)
            if overlap > 0:
                cov[i] += overlap
    return cov


def _coverage_label(
    intervals: list[BedInterval], label: str, win_start: int, n_bins: int, resolution: int
) -> list[int]:
    """Per-bin bp coverage of intervals whose name is exactly ``label``."""
    cov = [0] * n_bins
    win_end = win_start + n_bins * resolution
    for iv in intervals:
        if iv.end <= win_start or iv.start >= win_end or iv.name != label:
            continue
        s = max(0, (iv.start - win_start) // resolution)
        e = min(n_bins, -(-(iv.end - win_start) // resolution))  # ceil division
        for i in range(s, e):
            bin_start = win_start + i * resolution
            overlap = min(iv.end, bin_start + resolution) - max(iv.start, bin_start)
            if overlap > 0:
                cov[i] += overlap
    return cov


def _crossover(
    cov_a: list[int],
    cov_b: list[int],
    win_start: int,
    n_bins: int,
    resolution: int,
    plateau_threshold: float,
) -> CumulativeRefinement | None:
    """Symmetric left/right cumulative-fraction crossover from two coverage curves.

    The shared core of :func:`cumulative_breakpoint` (chromosome sides) and
    :func:`cumulative_breakpoint_between` (band sides): orient the curves, build
    the rising/falling cumulative fractions, and return the max-score position
    (inferred breakpoint) and its plateau.
    """
    # Orient: the side with more content in the left half is the "left" one.
    mid = n_bins // 2
    a_is_left = sum(cov_a[:mid]) >= sum(cov_b[:mid])
    left_cov, right_cov = (cov_a, cov_b) if a_is_left else (cov_b, cov_a)

    total_left = sum(left_cov)
    total_right = sum(right_cov)
    if total_left == 0 or total_right == 0:
        return None

    # left-cumulative: proportion at or before each bin (rises 0 -> 1).
    left_fraction: list[float] = []
    running = 0
    for v in left_cov:
        running += v
        left_fraction.append(running / total_left)

    # right-cumulative: proportion at or after each bin (falls 1 -> 0).
    right_fraction = [0.0] * n_bins
    running = 0
    for i in range(n_bins - 1, -1, -1):
        running += right_cov[i]
        right_fraction[i] = running / total_right

    positions = [win_start + i * resolution for i in range(n_bins)]
    combined = [left_fraction[i] + right_fraction[i] for i in range(n_bins)]

    max_score = max(combined)
    cutoff = max_score * plateau_threshold
    plateau_idx = [i for i, c in enumerate(combined) if c >= cutoff]
    plateau_start = positions[plateau_idx[0]]
    plateau_end = positions[plateau_idx[-1]] + resolution
    # Inferred breakpoint = midpoint of the maximal-score plateau (within 1e-9).
    top = [i for i, c in enumerate(combined) if c >= max_score - 1e-9]
    inferred = (positions[top[0]] + positions[top[-1]]) // 2

    return CumulativeRefinement(
        inferred_breakpoint=inferred,
        plateau_start=plateau_start,
        plateau_end=plateau_end,
        max_score=max_score,
        threshold=plateau_threshold,
        window_start=win_start,
        window_end=win_start + n_bins * resolution,
        resolution=resolution,
        positions=positions,
        left_fraction=left_fraction,
        right_fraction=right_fraction,
        combined_score=combined,
    )


def cumulative_breakpoint(
    bp: Breakpoint,
    intervals: list[BedInterval],
    reference: Reference,
    window: tuple[int, int],
    *,
    resolution: int = DEFAULT_RESOLUTION,
    plateau_threshold: float = DEFAULT_PLATEAU_THRESHOLD,
    pivot: str = DEFAULT_PIVOT,
    chrom_of: ChromResolver | None = None,
) -> CumulativeRefinement | None:
    """Infer a *translocation* breakpoint by cumulative k-mer analysis over ``window``.

    The two sides are the fusion partners' **chromosomes**. ``intervals`` are the
    scaffolded (pre-bin) records for the breakpoint's scaffold, in whichever label
    space ``pivot`` names -- pass the same pivot (and the same BED) used for
    :func:`~karyoscope_iscn.core.refine.refine_breakpoint`. ``chrom_of`` overrides
    ``pivot`` with a ready-made resolver. Returns ``None`` when either partner has
    no chromosome-specific content in the window.
    """
    chrom_a, chrom_b = bp.up_chrom, bp.down_chrom
    if not chrom_a or not chrom_b or chrom_a == chrom_b:
        return None

    win_start = max(0, window[0])
    n_bins = (window[1] - win_start) // resolution
    if n_bins < 1:
        return None

    resolve = chrom_of or resolver_for(pivot, reference)
    cov_a = _coverage(intervals, chrom_a, resolve, win_start, n_bins, resolution)
    cov_b = _coverage(intervals, chrom_b, resolve, win_start, n_bins, resolution)
    return _crossover(cov_a, cov_b, win_start, n_bins, resolution, plateau_threshold)


def cumulative_breakpoint_between(
    intervals: list[BedInterval],
    label_a: str,
    label_b: str,
    window: tuple[int, int],
    *,
    resolution: int = DEFAULT_RESOLUTION,
    plateau_threshold: float = DEFAULT_PLATEAU_THRESHOLD,
) -> CumulativeRefinement | None:
    """Infer an *intra-chromosomal* junction by cumulative k-mer analysis.

    Same machinery as :func:`cumulative_breakpoint`, but the two sides are two
    **cytoband labels** (``label_a`` / ``label_b``) rather than chromosomes -- so
    it localises a deletion/inversion/duplication junction where one band gives
    way to the other on the same chromosome. ``window`` is a *local* span around
    the junction (the flanking band runs' scaffold extent), so distant repeats of
    the same band don't bias the crossover. Returns ``None`` when a side has no
    content in the window.
    """
    if not label_a or not label_b or label_a == label_b:
        return None
    win_start = max(0, window[0])
    n_bins = (window[1] - win_start) // resolution
    if n_bins < 1:
        return None
    cov_a = _coverage_label(intervals, label_a, win_start, n_bins, resolution)
    cov_b = _coverage_label(intervals, label_b, win_start, n_bins, resolution)
    return _crossover(cov_a, cov_b, win_start, n_bins, resolution, plateau_threshold)
