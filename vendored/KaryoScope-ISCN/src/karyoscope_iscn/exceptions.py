"""Exception hierarchy for KaryoScope-ISCN.

All errors raised intentionally by this package derive from
:class:`KaryoscopeIscnError`, so callers (and the CLI) can catch the whole
family with a single ``except``.
"""

from __future__ import annotations


class KaryoscopeIscnError(Exception):
    """Base class for all KaryoScope-ISCN errors."""


class DatabaseError(KaryoscopeIscnError):
    """A k-mer database is missing, malformed, or inconsistent."""


class ReferenceError(KaryoscopeIscnError):
    """Precomputed reference data is missing or could not be derived."""


class KaryoscopeRunError(KaryoscopeIscnError):
    """A shelled-out ``karyoscope`` invocation failed."""


class InputError(KaryoscopeIscnError):
    """A required input file is missing or malformed."""
