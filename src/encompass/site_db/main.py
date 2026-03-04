# src/encompass/site_db/main.py

from __future__ import annotations

from .encompassService import menu, create_pdb


def main(argv: list[str] | None = None) -> None:
    """
    Entry point for the site-db subcommand.

    For now, this just runs create_pdb() as a simple test/demo.
    Later you could extend this to parse arguments (e.g. DB URL,
    table selection, etc.) and/or expose the interactive menu().
    """
    # Simple behavior for now: create a sample PDB entry
    create_pdb()

    # If you later want an interactive CLI, you could do:
    # while True:
    #     choice = menu()
    #     if choice == "1":
    #         create_pdb()
    #     elif choice == "5":
    #         break
    #     else:
    #         print("Invalid choice. Try again.")


if __name__ == "__main__":
    main()

