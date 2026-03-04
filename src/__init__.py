"""
EncoMPASS: tools to build and analyze the EncoMPASS membrane protein database.

Subpackages
-----------
- encompass.pipeline           : high-level pipeline orchestration
- encompass.sources            : fetching and combining data from OPM, PDB, UniProt, etc.
- encompass.struct_comparisons : structural comparisons and DB-wide analyses
- encompass.symmetry           : symmetry detection and related tools
- encompass.site_db            : export into the PostgreSQL database used by the site
- encompass.utils              : small, generic helper utilities
- encompass.data               : packaged reference tables and templates
"""

# Optional: version; keep in sync with your packaging metadata
__version__ = "0.1.0"

# Optional: expose a convenient top-level entrypoint if you have one
try:
    from .pipeline.run_encompass import run_full_pipeline
except ImportError:
    # During early refactors / partial installs, it's OK if this isn't available
    run_full_pipeline = None

__all__ = []
