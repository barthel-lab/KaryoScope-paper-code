"""Core business logic for KaryoScope-ISCN (no CLI awareness).

Database/reference loading, BED merging, ISCN report generation, and figure
rendering live here so they can be unit-tested and reused independently of the
``click`` command layer.
"""
