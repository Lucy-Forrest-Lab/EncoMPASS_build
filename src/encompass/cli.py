# src/encompass/cli.py

from __future__ import annotations

import sys
from typing import List


def main(argv: List[str] | None = None) -> None:
    """
    Top-level CLI for EncoMPASS.

    Subcommands:
      - structures       : run the full structure pipeline
                         (wraps encompass.pipeline.run_encompass.main)
      - symmetry       : run a single symmetry step
                         (wraps encompass.symmetry.run_symmetry_step.main)
      - site-db        : run site database utilities / demo
                         (wraps encompass.site_db.main.main)
      - init-templates : copy template files to a user directory
                         (wraps encompass.pipeline.templates.main)
    """
    if argv is None:
        argv = sys.argv[1:]

    if not argv:
        print(
            "Usage:\n"
            "  encompass structures [structures-args...]\n"
            "  encompass symmetry [symmetry-args...]\n"
            "  encompass site-db [site-db-args...]\n"
            "  encompass init-templates [--dest PATH]\n\n"
            "Run 'encompass <subcommand> -h' for details on that subcommand.\n"
        )
        sys.exit(1)

    cmd, *rest = argv

    if cmd == "structures":
        # Full structure pipeline
        from .pipeline import run_encompass as run_enc
        run_enc.main(rest)

    elif cmd == "symmetry":
        # Single symmetry step
        from .symmetry import run_symmetry_step as rss
        rss.main(rest)

    elif cmd in ("site-db", "sitedb"):
        # Site database tools / demo
        from .site_db import main as site_main
        site_main.main(rest)

    elif cmd == "init-templates":
        # Copy packaged templates into a user-editable directory
        from .pipeline import templates as tmpl
        tmpl.main(rest)

    elif cmd in ("-h", "--help", "help"):
        print(
            "Usage:\n"
            "  encompass structures [structures-args...]\n"
            "  encompass symmetry [symmetry-args...]\n"
            "  encompass site-db [site-db-args...]\n"
            "  encompass init-templates [--dest PATH]\n\n"
            "Run 'encompass <subcommand> -h' for details on that subcommand.\n"
        )
        sys.exit(0)

    else:
        print(f"Unknown subcommand: {cmd!r}\n", file=sys.stderr)
        print(
            "Available subcommands:\n"
            "  structures\n"
            "  symmetry\n"
            "  site-db\n"
            "  init-templates\n",
            file=sys.stderr,
        )
        sys.exit(1)
