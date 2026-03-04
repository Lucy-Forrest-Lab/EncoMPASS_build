# src/encompass/pipeline/templates.py

from __future__ import annotations

import shutil
from pathlib import Path
from importlib.resources import files, as_file
import argparse


# Package layout constants
_DATA_PKG = "encompass.data"
_TEMPLATES_SUBDIR = "templates"
_REFERENCE_SUBDIR = "reference"

# Placeholders used inside data/templates/instructions.txt
_INSTR_DEL_PLACEHOLDER = "<path to file with pdb codes that should be ignored, e.g. delete_list.txt>"
_INSTR_FORCE_PLACEHOLDER = "<path to file that insists on particular source for the structure of a pdb code, e.g. replace_list.txt>"
_INSTR_DST_PLACEHOLDER = "<path to data structure template, e.g. str_data_entry_current.json>"


def _copy_other_templates(templates_dir, dest: Path, overwrite: bool = False) -> None:
    """
    Copy all template files except instructions.txt from templates_dir into dest.
    """
    for entry in templates_dir.iterdir():
        if not entry.is_file():
            continue
        if entry.name == "instructions.txt":
            continue  # handled separately

        target = dest / entry.name
        if target.exists() and not overwrite:
            # Don't clobber user changes
            continue

        with as_file(entry) as src_path:
            shutil.copy(src_path, target)


def _write_instructions_with_paths(templates_dir, dest: Path, overwrite: bool = False) -> None:
    """
    Read the packaged instructions.txt template, substitute placeholders
    with the actual paths to delete_list.txt, replace_list.txt and
    str_data_entry_current.json inside the installed package, and write
    the result into dest.
    """
    instr_res = templates_dir.joinpath("instructions.txt")
    if not instr_res.is_file():
        # Nothing to do if there is no instructions template
        return

    target = dest / "instructions.txt"
    if target.exists() and not overwrite:
        # Don't overwrite an existing user-edited instructions file
        return

    data_root = files(_DATA_PKG)
    ref_dir = data_root.joinpath(_REFERENCE_SUBDIR)

    del_res = ref_dir.joinpath("delete_list.txt")
    forcedb_res = ref_dir.joinpath("replace_list.txt")
    dst_res = ref_dir.joinpath("str_data_entry_current.json")

    # Use as_file so this works even if the package is installed as a wheel/zip
    with as_file(instr_res) as instr_path, \
         as_file(del_res) as del_path, \
         as_file(forcedb_res) as forcedb_path, \
         as_file(dst_res) as dst_path:

        text = Path(instr_path).read_text()

        text = text.replace(_INSTR_DEL_PLACEHOLDER, str(del_path))
        text = text.replace(_INSTR_FORCE_PLACEHOLDER, str(forcedb_path))
        text = text.replace(_INSTR_DST_PLACEHOLDER, str(dst_path))

    target.write_text(text)


def copy_templates(dest_dir: str | Path = ".", overwrite: bool = False) -> None:
    """
    Copy all template files from encompass.data.templates into dest_dir.

    - All templates except instructions.txt are copied verbatim.
    - instructions.txt is treated as a template: placeholder strings
      for delete_list.txt, replace_list.txt and str_data_entry_current.json
      are replaced with the actual paths to those reference files inside
      the installed package.

    Existing files are not overwritten unless overwrite=True.
    """
    dest = Path(dest_dir)
    dest.mkdir(parents=True, exist_ok=True)

    data_root = files(_DATA_PKG)
    templates_dir = data_root.joinpath(_TEMPLATES_SUBDIR)

    _copy_other_templates(templates_dir, dest, overwrite=overwrite)
    _write_instructions_with_paths(templates_dir, dest, overwrite=overwrite)


def main(argv: list[str] | None = None) -> None:
    """
    CLI entrypoint used by `encompass init-templates`.

    Examples
    --------
    Copy templates into the current directory:

        encompass init-templates

    Copy into a custom directory (created if needed):

        encompass init-templates --dest config/ --overwrite
    """
    parser = argparse.ArgumentParser(
        description=(
            "Copy EncoMPASS template files into a user-editable directory. "
            "instructions.txt is updated to point at packaged reference files."
        )
    )
    parser.add_argument(
        "--dest",
        default=".",
        help="Destination directory (default: current directory).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing files in the destination.",
    )
    args = parser.parse_args(argv)

    copy_templates(dest_dir=args.dest, overwrite=args.overwrite)
