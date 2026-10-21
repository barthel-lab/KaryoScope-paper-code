"""Breakpoint refinement from the finest-resolution scaffolded annotation.

A translocation breakpoint detected from the *binned* cytoband BED only knows
the junction to within a bin (the segment boundary, :attr:`Breakpoint.pos`).
The flanking interval emitted into ``_regions.tsv`` is, by default, just a fixed
symmetric pad around that estimate. This module sharpens both: given the
**scaffolded (pre-bin) cytoband BED** -- the finest resolution we retain -- it
locates where one flanking chromosome's content gives way to the other and
returns the midpoint of the intervening ambiguous zone as the refined centre,
plus that zone as a tight interval.

The algorithm is the content-dominance method ported from the legacy
``KaryoScope_assembly_contig_zoom_plot.py`` zoom plotter, de-hardcoded for the
cytoband pivot: chromosome identity comes from parsing each interval's cytoband
label (:func:`~karyoscope_iscn.core.cytoband.parse_band` ->
:func:`~karyoscope_iscn.core.cytoband.chrom_for_token`) rather than from
``chrN_specific`` chromosome-featureset names.

Refinement is pure analysis (no rendering). The report runs it by default so
``_regions.tsv`` carries exact coordinates and the zoom plot can consume them
directly without re-deriving anything.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Callable
from dataclasses import dataclass

from karyoscope_iscn.core.bedio import BedInterval
from karyoscope_iscn.core.cytoband import chrom_for_token, parse_band
from karyoscope_iscn.core.detect import Breakpoint
from karyoscope_iscn.reference import Reference

#: Aggregation window for the coarse dominance pass.
DEFAULT_DOMINANCE_WINDOW = 100_000
#: Bin size for the sharp-transition drill-down (finer than any band boundary,
#: but coarse enough that scattered cross-mapped specks don't move the call).
DEFAULT_FINE_BIN = 10_000
#: A window counts as dominated only if the winning chromosome covers at least
#: this fraction of the *window* (not just of chromosome-specific content).
DEFAULT_MIN_WINDOW_FRACTION = 0.25
#: ...and at least this fraction of the two flanking chromosomes' content.
DEFAULT_SPECIFICITY = 0.50
#: Symmetric pad used as the interval fallback when no transition is resolved.
DEFAULT_FALLBACK_PAD = 500_000
#: How far (in dominance windows) the last upstream-dominated window may sit PAST
#: the first downstream-dominated one before the scan is judged to have found no
#: transition. A small overlap is the normal ambiguous zone; a large one means the
#: two chromosomes' content is interleaved, or the coarse ``bp.pos`` that
#: :func:`_orient` keys off is wrong and upstream/downstream came out backwards.
DEFAULT_MAX_DOMINANCE_OVERLAP_WINDOWS = 2


@dataclass(frozen=True)
class RefinedBreakpoint:
    """A breakpoint's sharpened centre and ambiguous-zone interval.

    ``refined`` is ``False`` when the content-dominance scan could not resolve a
    transition and the values fall back to ``pos`` and a symmetric pad.
    """

    center: int
    start: int
    end: int
    refined: bool


#: Maps a BED interval's label to a reference chromosome, or ``None`` when the
#: label carries no chromosome identity (an internal hierarchy node, ``novel``,
#: or a band from an unknown chromosome). Only labels resolving to one of the
#: two flanking chromosomes contribute to the dominance scan, so ``None`` is a
#: normal, silent outcome.
ChromResolver = Callable[[str], "str | None"]


def cytoband_resolver(reference: Reference) -> ChromResolver:
    """Resolve chromosome identity by parsing a cytoband label (``13q12`` -> ``chr13``)."""
    cache: dict[str, str | None] = {}

    def resolve(name: str) -> str | None:
        if name not in cache:
            cache[name] = chrom_for_token(parse_band(name).chrom_token, reference)
        return cache[name]

    return resolve


def chromosome_resolver(reference: Reference) -> ChromResolver:
    """Resolve chromosome identity from a ``chromosome`` feature-set label.

    That feature set labels each k-mer with the chromosome it is specific to, so
    the label *is* the chromosome name (``chr13``). Internal hierarchy nodes
    (``autosome``, ``acrocentric``, ``categorized``) and ``novel`` carry no
    chromosome identity and resolve to ``None``.

    This is the pivot the published Robertsonian breakpoint analysis used; the
    cytoband pivot is the same dominance method over a different label space.
    """
    cache: dict[str, str | None] = {}
    known = reference.chromosomes

    def resolve(name: str) -> str | None:
        if name not in cache:
            cache[name] = name if name in known else None
        return cache[name]

    return resolve


#: Pivot name -> resolver factory. ``cytoband`` preserves prior behaviour.
PIVOTS: dict[str, Callable[[Reference], ChromResolver]] = {
    "cytoband": cytoband_resolver,
    "chromosome": chromosome_resolver,
}
DEFAULT_PIVOT = "cytoband"


def resolver_for(pivot: str, reference: Reference) -> ChromResolver:
    """Return the resolver for ``pivot``; raises on an unknown name."""
    try:
        factory = PIVOTS[pivot]
    except KeyError:
        raise ValueError(f"unknown refine pivot {pivot!r}; known: {sorted(PIVOTS)}") from None
    return factory(reference)


def _window_content(
    intervals: list[BedInterval],
    chrom1: str,
    chrom2: str,
    chrom_of: ChromResolver,
    window: int,
) -> dict[int, dict[str, int]]:
    """Aggregate per-chromosome covered bp into ``window``-sized bins.

    Only intervals resolving to ``chrom1`` or ``chrom2`` contribute; an interval
    spanning several windows is split proportionally by overlap.
    """
    content: dict[int, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for iv in intervals:
        chrom = chrom_of(iv.name)
        if chrom != chrom1 and chrom != chrom2:
            continue
        w_start = (iv.start // window) * window
        while w_start < iv.end:
            w_end = w_start + window
            overlap = min(iv.end, w_end) - max(iv.start, w_start)
            if overlap > 0:
                content[w_start][chrom] += overlap
            w_start += window
    return content


def _orient(
    content: dict[int, dict[str, int]],
    chrom1: str,
    chrom2: str,
    pivot: int,
) -> tuple[str, str]:
    """Return ``(upstream_chrom, downstream_chrom)`` by content before ``pivot``."""
    c1_before = sum(c.get(chrom1, 0) for w, c in content.items() if w < pivot)
    c2_before = sum(c.get(chrom2, 0) for w, c in content.items() if w < pivot)
    if c1_before >= c2_before:
        return chrom1, chrom2
    return chrom2, chrom1


def _dominant(
    content: dict[str, int],
    chrom: str,
    other: str,
    window: int,
    min_window_fraction: float,
    specificity: float,
) -> bool:
    """True when ``chrom`` dominates a window over ``other``."""
    c_bp = content.get(chrom, 0)
    total = c_bp + content.get(other, 0)
    return c_bp >= window * min_window_fraction and total > 0 and c_bp / total >= specificity


def _drilldown_center(
    intervals: list[BedInterval],
    upstream: str,
    downstream: str,
    chrom_of: ChromResolver,
    lo: int,
    hi: int,
    fine_bin: int,
) -> int | None:
    """Find the sharp transition centre at ``fine_bin`` resolution in ``[lo, hi)``.

    Each fine bin is won by whichever flanking chromosome covers the most bp in
    it. The downstream chromosome's start is only honoured once it wins two
    consecutive bins (a real block), guarding against isolated cross-mapped
    specks. Returns the midpoint of the last-upstream / first-downstream wins,
    or ``None`` if no clean transition is found.
    """
    bins: dict[int, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for iv in intervals:
        chrom = chrom_of(iv.name)
        if chrom != upstream and chrom != downstream:
            continue
        if iv.end <= lo or iv.start >= hi:
            continue
        sb = max(iv.start, lo) // fine_bin
        eb = (min(iv.end, hi) - 1) // fine_bin
        for b_idx in range(sb, eb + 1):
            bs = b_idx * fine_bin
            be = bs + fine_bin
            overlap = min(iv.end, be) - max(iv.start, bs)
            if overlap > 0:
                bins[b_idx][chrom] += overlap

    last_up_pos: int | None = None
    first_down_pos: int | None = None
    prev_was_down = False
    for b_idx in sorted(bins):
        winner = max(bins[b_idx].items(), key=lambda kv: kv[1])[0]
        pos = b_idx * fine_bin
        if winner == upstream:
            last_up_pos = pos + fine_bin
        elif winner == downstream and first_down_pos is None and prev_was_down:
            first_down_pos = pos - fine_bin
        prev_was_down = winner == downstream

    # No consecutive downstream block: accept a directly-adjacent first win.
    if first_down_pos is None and last_up_pos is not None:
        for b_idx in sorted(bins):
            winner = max(bins[b_idx].items(), key=lambda kv: kv[1])[0]
            pos = b_idx * fine_bin
            if winner == downstream and pos <= last_up_pos + fine_bin:
                first_down_pos = pos
                break

    if last_up_pos is not None and first_down_pos is not None:
        return (last_up_pos + first_down_pos) // 2
    return None


def refine_breakpoint(
    bp: Breakpoint,
    intervals: list[BedInterval],
    reference: Reference,
    *,
    window: int = DEFAULT_DOMINANCE_WINDOW,
    fine_bin: int = DEFAULT_FINE_BIN,
    min_window_fraction: float = DEFAULT_MIN_WINDOW_FRACTION,
    specificity: float = DEFAULT_SPECIFICITY,
    fallback_pad: int = DEFAULT_FALLBACK_PAD,
    max_overlap_windows: int = DEFAULT_MAX_DOMINANCE_OVERLAP_WINDOWS,
    pivot: str = DEFAULT_PIVOT,
    chrom_of: ChromResolver | None = None,
) -> RefinedBreakpoint:
    """Sharpen one breakpoint's centre and interval from scaffolded intervals.

    ``intervals`` are the scaffolded (pre-bin) cytoband BED records **for this
    breakpoint's scaffold** (already filtered to ``bp.scaffold``). The two
    flanking chromosomes are ``bp.up_chrom`` / ``bp.down_chrom``.

    Falls back to ``bp.pos`` with a symmetric ``fallback_pad`` interval (and
    ``refined=False``) whenever the transition cannot be resolved.
    """
    fallback = RefinedBreakpoint(
        center=bp.pos,
        start=max(0, bp.pos - fallback_pad),
        end=bp.pos + fallback_pad,
        refined=False,
    )
    chrom1, chrom2 = bp.up_chrom, bp.down_chrom
    if not chrom1 or not chrom2 or chrom1 == chrom2:
        return fallback

    if chrom_of is None:
        chrom_of = resolver_for(pivot, reference)
    content = _window_content(intervals, chrom1, chrom2, chrom_of, window)
    if not content:
        return fallback

    upstream, downstream = _orient(content, chrom1, chrom2, bp.pos)

    last_upstream_end: int | None = None
    first_downstream_start: int | None = None
    for w in sorted(content):
        c = content[w]
        if _dominant(c, upstream, downstream, window, min_window_fraction, specificity):
            last_upstream_end = w + window
        if first_downstream_start is None and _dominant(
            c, downstream, upstream, window, min_window_fraction, specificity
        ):
            first_downstream_start = w

    if last_upstream_end is None or first_downstream_start is None:
        return fallback

    # No clean transition: upstream dominance persisting well past the start of
    # downstream dominance is interleaved content, not a junction. Taking
    # min()/max() anyway yields an "interval" spanning the whole contig with a
    # meaningless midpoint centre -- still flagged refined -- which is how
    # HG01175's dic(6;20) reported [0, 172,200,000] and a cumulative breakpoint
    # OUTSIDE its own interval. Report unrefined instead.
    if last_upstream_end - first_downstream_start > max_overlap_windows * window:
        return fallback

    start = min(last_upstream_end, first_downstream_start)
    end = max(last_upstream_end, first_downstream_start)
    gap = first_downstream_start - last_upstream_end

    if gap > window:
        # Wide ambiguous zone -- the midpoint is the best centre estimate.
        center = (last_upstream_end + first_downstream_start) // 2
    else:
        # Sharp/narrow transition -- drill down to place the centre precisely.
        lo = max(0, last_upstream_end - window)
        hi = first_downstream_start + 20 * window
        center = _drilldown_center(intervals, upstream, downstream, chrom_of, lo, hi, fine_bin)
        if center is None:
            center = (last_upstream_end + first_downstream_start) // 2

    return RefinedBreakpoint(center=center, start=start, end=end, refined=True)
