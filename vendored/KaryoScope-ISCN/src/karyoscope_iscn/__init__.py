"""KaryoScope-ISCN: ISCN cytogenetic reports and karyotype figures.

A consumer of KaryoScope cytoband/region annotations that produces ISCN-style
cytogenetic reports and publication karyotype figures. Feature definitions
(chromosomes, cytobands, regions, ...) are loaded dynamically from the k-mer
database ``hierarchy.tsv``; nothing about the reference is hardcoded.
"""

from karyoscope_iscn._version import __version__

__all__ = ["__version__"]
