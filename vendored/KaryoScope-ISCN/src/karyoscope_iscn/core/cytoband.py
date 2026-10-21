"""The cytoband label model.

A cytoband annotation BED labels each interval of a *sample* scaffold with the
CHM13 cytoband it matches, e.g. ``13q21.1``. The label encodes everything the
ISCN report needs about reference identity:

* the **chromosome** token (``13``),
* the **arm** (``q``), and
* the **sub-band** (``21.1``).

Some labels are not bands: ``categorized`` (resolved to the genome but not a
band), ``novel`` (k-mer absent from the index), and the bare chromosome token
(``13``, resolved to a chromosome but not an arm/band). These are represented
with ``arm``/``sub`` (and, for the sentinels, ``chrom_token``) set to ``None``.

Parsing is pure string work over the ISCN naming convention -- it does not
depend on any reference. Use :func:`chrom_for_token` to map a token to a
reference chromosome name (``13`` -> ``chr13``).
"""

from __future__ import annotations

from dataclasses import dataclass

from karyoscope_iscn.reference import Reference

#: Non-band sentinel labels emitted by the annotator.
CATEGORIZED = "categorized"
NOVEL = "novel"
SENTINELS = frozenset({CATEGORIZED, NOVEL})


@dataclass(frozen=True)
class Band:
    """A parsed cytoband label."""

    raw: str
    chrom_token: str | None
    arm: str | None
    sub: str | None

    @property
    def is_sentinel(self) -> bool:
        """True for ``categorized`` / ``novel`` (no chromosome)."""
        return self.raw in SENTINELS

    @property
    def has_arm(self) -> bool:
        """True when the label resolves to a specific arm (and thus a band)."""
        return self.arm is not None

    @property
    def is_chromosomal(self) -> bool:
        """True when the label maps to a chromosome (band or bare chromosome)."""
        return self.chrom_token is not None


def _split_arm(label: str) -> tuple[str, str | None, str | None]:
    """Split a band label into (chrom_token, arm, sub) on the first p/q."""
    for i, ch in enumerate(label):
        if ch in "pq":
            token = label[:i]
            sub = label[i + 1 :] or None
            return token, ch, sub
    # No arm letter: the whole label is the chromosome token (bare chromosome).
    return label, None, None


def parse_band(label: str) -> Band:
    """Parse a cytoband label into a :class:`Band`.

    Sentinels (``categorized``/``novel``) yield a band with all fields ``None``.
    A bare chromosome token (``"13"``) yields ``chrom_token`` set with no arm.
    """
    label = label.strip()
    if label in SENTINELS or not label:
        return Band(raw=label, chrom_token=None, arm=None, sub=None)
    token, arm, sub = _split_arm(label)
    return Band(raw=label, chrom_token=token or None, arm=arm, sub=sub)


def chrom_for_token(token: str | None, reference: Reference) -> str | None:
    """Map a band chromosome token (``"13"``) to a reference chromosome (``"chr13"``).

    Matches by normalizing reference names (dropping a leading ``chr``) against
    the token, so it works whether or not the reference uses a ``chr`` prefix.
    Returns ``None`` for a ``None`` token or one absent from the reference.
    """
    if token is None:
        return None
    for name in reference.chrom_names():
        normalized = name[3:] if name.lower().startswith("chr") else name
        if normalized == token or name == token:
            return name
    return None
