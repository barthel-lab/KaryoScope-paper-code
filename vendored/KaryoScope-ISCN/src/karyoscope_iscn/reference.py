"""Precomputed reference data (chromosome sizes, centromeres, acrocentric flags,
canonical cytoband ideogram).

These artifacts replace the formerly-hardcoded ``CHM13_CENTROMERES`` /
``CHM13_CHROM_SIZES`` / ``ACROCENTRIC_CHROMS`` tables. They are *derived once*
from stable reference inputs (the assembly ``.fai`` for sizes, the smoothed
CHM13 cytoband reference BED for centromere/arm boundaries and band
coordinates, and the database's chromosome hierarchy for acrocentric grouping)
by the ``precompute`` command, then committed under
``karyoscope_iscn/data/references/<name>/`` so the pipeline never re-runs
KaryoScope on the reference.

A reference directory contains:

* ``chromosomes.tsv`` -- one row per chromosome: size, p-arm end, q-arm start,
  centromere midpoint, acrocentric flag.
* ``cytobands.bed`` -- the canonical band ideogram (``chrom  start  end  band``).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from karyoscope_iscn.exceptions import ReferenceError

#: Packaged reference data shipped with the wheel.
DEFAULT_REFERENCES_DIR = Path(__file__).parent / "data" / "references"

CHROMOSOMES_FILE = "chromosomes.tsv"
#: Optional per-chromosome mappable-fraction table (see the file for provenance).
MAPPABILITY_FILE = "mappability.tsv"
IDEOGRAM_FILE = "cytobands.bed"

_CHROM_HEADER = ("chrom", "size", "p_end", "q_start", "centromere", "acrocentric")


@dataclass(frozen=True)
class ChromInfo:
    """Reference facts for one chromosome (all coordinates 0-based bp)."""

    chrom: str
    size: int
    p_end: int | None
    q_start: int | None
    centromere: int | None
    acrocentric: bool

    def arm_of(self, pos: int) -> str | None:
        """Return ``'p'`` or ``'q'`` for ``pos`` relative to the centromere.

        Returns ``None`` when the centromere is unknown (cannot decide).
        """
        if self.centromere is None:
            return None
        return "p" if pos < self.centromere else "q"


@dataclass
class Reference:
    """A loaded reference: per-chromosome facts plus the optional ideogram path."""

    name: str
    chromosomes: dict[str, ChromInfo]
    ideogram_path: Path | None = None
    #: Per-chromosome mappable fraction (see ``mappability.tsv``). Copy number
    #: divides expected_bp by this. Missing chromosome -> 1.0, so a reference with
    #: no table behaves exactly as before.
    mappability: dict[str, float] = field(default_factory=dict)

    def mappable_fraction(self, chrom: str) -> float:
        """Fraction of ``chrom`` the k-mer index can reach; 1.0 when uncalibrated."""
        return self.mappability.get(chrom, 1.0)

    def chrom_names(self) -> list[str]:
        """Chromosome names in file order."""
        return list(self.chromosomes)

    def info(self, chrom: str) -> ChromInfo:
        try:
            return self.chromosomes[chrom]
        except KeyError:
            raise ReferenceError(
                f"chromosome {chrom!r} is not in reference {self.name!r}; "
                f"known: {self.chrom_names()!r}"
            ) from None

    def size(self, chrom: str) -> int:
        return self.info(chrom).size

    def centromere(self, chrom: str) -> int | None:
        return self.info(chrom).centromere

    def acrocentrics(self) -> set[str]:
        return {c for c, i in self.chromosomes.items() if i.acrocentric}

    def arm_of(self, chrom: str, pos: int) -> str | None:
        return self.info(chrom).arm_of(pos)


def parse_chromosomes_tsv(path: Path) -> dict[str, ChromInfo]:
    """Parse a ``chromosomes.tsv`` into ``{chrom: ChromInfo}``."""
    if not path.is_file():
        raise ReferenceError(f"reference chromosomes file not found: {path}")
    lines = path.read_text().splitlines()
    if not lines:
        raise ReferenceError(f"reference chromosomes file is empty: {path}")
    header = tuple(h.strip() for h in lines[0].split("\t"))
    if header != _CHROM_HEADER:
        raise ReferenceError(
            f"{path} has unexpected header {list(header)!r}, expected {list(_CHROM_HEADER)!r}"
        )

    def _opt_int(value: str) -> int | None:
        value = value.strip()
        return None if value in ("", "NA", "None") else int(value)

    out: dict[str, ChromInfo] = {}
    for i, raw in enumerate(lines[1:], start=2):
        if not raw.strip():
            continue
        parts = raw.split("\t")
        if len(parts) != len(_CHROM_HEADER):
            raise ReferenceError(
                f"{path}:{i}: expected {len(_CHROM_HEADER)} columns, got {len(parts)}"
            )
        chrom = parts[0].strip()
        out[chrom] = ChromInfo(
            chrom=chrom,
            size=int(parts[1]),
            p_end=_opt_int(parts[2]),
            q_start=_opt_int(parts[3]),
            centromere=_opt_int(parts[4]),
            acrocentric=parts[5].strip() in ("1", "true", "True"),
        )
    return out


def load_reference(name: str, references_dir: Path | None = None) -> Reference:
    """Load a committed reference by ``name`` (e.g. ``"chm13v2.0"``).

    Raises :class:`ReferenceError` if the reference directory or its
    ``chromosomes.tsv`` is missing.
    """
    base = references_dir or DEFAULT_REFERENCES_DIR
    ref_dir = base / name
    if not ref_dir.is_dir():
        available = sorted(p.name for p in base.iterdir() if p.is_dir()) if base.is_dir() else []
        raise ReferenceError(
            f"reference {name!r} not found under {base}. Available: {available!r}. "
            "Generate it with `karyoscope-iscn precompute`."
        )
    chromosomes = parse_chromosomes_tsv(ref_dir / CHROMOSOMES_FILE)
    ideogram = ref_dir / IDEOGRAM_FILE
    return Reference(
        name=name,
        chromosomes=chromosomes,
        ideogram_path=ideogram if ideogram.is_file() else None,
        mappability=parse_mappability_tsv(ref_dir / MAPPABILITY_FILE),
    )


def parse_mappability_tsv(path: Path) -> dict[str, float]:
    """Parse an optional ``mappability.tsv`` into ``{chrom: factor}``.

    Absent file -> empty mapping -> every chromosome scores 1.0, i.e. the
    uncalibrated behaviour. Blank lines and ``#`` comments are skipped.
    """
    if not path.is_file():
        return {}
    out: dict[str, float] = {}
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or line.startswith("chrom\t"):
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        try:
            factor = float(parts[1])
        except ValueError:
            continue
        if factor <= 0:
            raise ReferenceError(f"{path}: non-positive mappability {parts[1]!r} for {parts[0]!r}")
        out[parts[0]] = factor
    return out


def available_references(references_dir: Path | None = None) -> list[str]:
    """List committed reference names."""
    base = references_dir or DEFAULT_REFERENCES_DIR
    if not base.is_dir():
        return []
    return sorted(p.name for p in base.iterdir() if (p / CHROMOSOMES_FILE).is_file())
