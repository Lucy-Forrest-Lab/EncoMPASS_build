# EncoMPASS

EncoMPASS is a pipeline for building, analyzing, and serving the **EncoMPASS membrane protein structure database**. It:

- Retrieves and normalizes structural data from PDB / OPM and related sources  
- Builds a curated EncoMPASS repository on disk  
- Runs large-scale **structure comparisons** and **symmetry analysis** (CE-Symm, SymD, QuatSymm, AnaNaS, MSSD)  
- Produces data structures that can be exported as XML (for the legacy website) or loaded into a **PostgreSQL** database (via `site_db`) for a modern web front-end.

---

## Project layout

This repository is now a Python package with a `src/` layout:

- `src/encompass/`
  - `pipeline/` – core repository & database build pipeline  
    - `run_encompass.py` – orchestrates the multi-stage build  
    - `initialize_repository.py`, `config.py`, `supporting_functions.py`, …
  - `sources/` – code for pulling data from source databases (PDB, UniProt, OPM, etc.)
  - `struct_comparisons/` – structure comparison / analysis pipeline (neighbors, plots, etc.)
  - `symmetry/` – symmetry pipeline (CE-Symm, SymD, QuatSymm, AnaNaS, MSSD, transfer)  
  - `site_db/` – tools for building the data structures used by the EncoMPASS website  
    - `create_data_struct.py` – builds a complete `web_data` structure from an EncoMPASS repository  
    - `create_xml_from_struct.py` – optionally renders XML from `web_data`  
    - `models.py`, `database.py`, `encompassService.py`, `dao.py` – SQLAlchemy models & DB utilities
  - `data/`
    - `reference/` – reference txt/json files used by the pipeline  
    - `templates/` – templates that can be copied and edited by users (e.g. instructions file)
  - `utils/` – assorted utilities and validation scripts
- `scripts/` – thin shell wrappers for batch / cluster runs (symmetry submissions, plotting helpers, etc.)
- `tests/` – Python tests (`encompass_tests.py`)

At runtime, EncoMPASS also expects an `encompass.env` file (environment variables) and a pipeline **instructions file** describing paths, database locations, and other options.

---

## Installation

EncoMPASS is currently intended to be installed **from source**.

```bash
# Clone this repository
git clone https://github.com/Lucy-Forrest-Lab/EncoMPASS.git
cd EncoMPASS

# Create and activate a virtual environment
python -m venv .venv
source .venv/bin/activate

# Install the package (plus optional extras)
python -m pip install --upgrade pip
pip install -e ".[all]"
````

> **Note:** The `.[all]` extra installs Python-side dependencies for:
>
> * the main pipeline
> * symmetry & analysis code
> * the site_db / SQLAlchemy integration.

External tools (PPM, MUSCLE, FrTM-Align, CE-Symm, SymD, QuatSymm, AnaNaS, etc.) are **not** installed by this package and must be available in your environment, typically via modules or the EncoMPASS container.

---

## Command-line interface (CLI)

Installing the package exposes a top-level `encompass` command that provides a small CLI wrapper around the main pipeline entry points.

Typical usage pattern:

```bash
# Show CLI help
encompass --help

# Run the main EncoMPASS pipeline
encompass run-pipeline \
    --main-path /path/to/encompass_repository \
    --instr-file EncoMPASS_options_relative_to_main.txt

# Run a symmetry update step (formerly single_str_update.py)
encompass run-symmetry-step \
    --db /path/to/EncoMPASS.db \
    --instr-file /path/to/instructions.txt \
    --step cesymm \
    --label cesymm_beta_example \
    --locusts-dir /path/to/locusts_dir \
    --str-type all

# Copy editable templates (instructions, etc.) into the current directory
encompass init-templates

# Build the aggregated web_data structure (for XML or PostgreSQL)
encompass build-site-data \
    --main-path /path/to/encompass_repository \
    --instr-file EncoMPASS_options_relative_to_main.txt \
    --output web_data.pkl
```

> Exact subcommand names/flags are defined in `src/encompass/cli.py`.
> The examples above match the intended usage: `run-pipeline`, `run-symmetry-step`, `init-templates`, and `build-site-data`.

You can still run modules directly via `python -m` if you prefer:

```bash
python -m encompass.pipeline.run_encompass -h
python -m encompass.symmetry.run_symmetry_step -h
python -m encompass.site_db.create_data_struct -h
python -m encompass.site_db.create_xml_from_struct -h
```

---

## Configuration & templates

Two key configuration files:

1. **`encompass.env`**
   Defines environment variables used by the pipeline and symmetry jobs (e.g. `ENC_DB`, `ENC_DB_INSTRUCT`, `LOCUSTS_TMP`, `PYTHON_ENV`, etc.).
   Your shell wrappers in `scripts/` typically `source` this file before running Python modules.

2. **Instruction file (template in `data/templates/instructions.txt`)**
   Controls:

   * database path
   * location of input/output folders
   * reference files (e.g. `delete_list.txt`, `replace_list.txt`, `str_data_entry_current.json`)
   * paths to external tools (PPM, MUSCLE, etc.)

The `encompass init-templates` command can be used to copy the templates into the current working directory (or a specified directory), where you can then edit them for your specific deployment.

Reference files (e.g. `define_locations.txt`, `delete_list.txt`, `deletion_codes.txt`, `str_data_entry_current.json`) now live under:

```text
src/encompass/data/reference/
```

and are resolved via the pipeline configuration, rather than by assuming they sit next to the Python scripts.

---

## Site database (`site_db`)

The `encompass.site_db` package provides tooling to export the processed EncoMPASS data into formats suitable for the public web interface:

* **Step 1: build a canonical `web_data` structure**

  ```bash
  encompass build-site-data \
      --main-path /path/to/encompass_repository \
      --instr-file EncoMPASS_options_relative_to_main.txt \
      --output web_data.pkl
  ```

  This reads the repository (str_data, analysis outputs, symmetry results, inferred symmetry, etc.), and builds a single `web_data` dict with **all information needed** by both:

  * the legacy XML site, and
  * a PostgreSQL-backed site using SQLAlchemy models.

* **Step 2 (optional): generate XML**

  ```bash
  python -m encompass.site_db.create_xml_from_struct \
      -w web_data.pkl \
      -o EncoMPASS
  ```

  This recreates the XML files that the original legacy site uses, but using only `web_data` as input.

* **Step 3: populate PostgreSQL**

  The modules in `encompass.site_db` (`models.py`, `database.py`, `encompassService.py`, `dao.py`) define SQLAlchemy models and helper functions for loading the same information into a PostgreSQL database. This is intended for powering a modern web front-end that mirrors everything previously exposed via XML.

---

## Dependencies & containers

For reproducible environments, we recommend using the existing container definition:

* Public container (from 2023-12-19) is available at:
  [https://github.com/Lucy-Forrest-Lab/EncoMPASS-containers-deps](https://github.com/Lucy-Forrest-Lab/EncoMPASS-containers-deps)

It includes:

* **PPM v2.0** to insert structures into the membrane when OPM does not have the desired biological assembly

  * [PPM website](https://opm.phar.umich.edu/ppm_server2_cgopm)
* **MUSCLE v3.8.31** for sequence alignments
* **FrTM-Align** for structure alignments

  * [FrTM-Align website](https://sites.gatech.edu/cssb/fr-tm-align/)
* Symmetry tools:

  * CE-Symm v2.2.3
  * QuatSymm v2.2.3
  * SymD v1.61 and v1.3w
  * AnaNaS v1.1

Python-side dependencies (e.g. `numpy`, `pandas`, `biopython`, `requests`, `sqlalchemy`, etc.) are specified in `pyproject.toml` and installed via `pip`.

---

## Release notes

### v1.1.0

* Refactored code into a **standard Python package** under `src/encompass/`
* Added a top-level **`encompass` CLI** with subcommands for:

  * running the main pipeline
  * running symmetry steps
  * initializing templates
  * building site data (`web_data`)
* Updated to be compatible with **Python 3.12**
* Introduced `site_db`:

  * aggregation into a single `web_data` object
  * optional XML generation from `web_data`
  * SQLAlchemy models for PostgreSQL export
* Fixed bugs in:

  * handling MUSCLE (updated to MUSCLE v5)
  * output folder specification in `complete_information`
* Added wrapper code `run_encompass.py` to allow dataset compilation to be run in stages
* Updated API calls to OPM and PDB to match current web services (as of 2025)
* Updated to newer versions of PPM (configuration-dependent)

### v1.0.0

#### Structure retrieval

* Added information about processing and decision-making steps to the header of each structure

#### Structure alignment

* TMs of all sequence-related chains are considered when deciding which comparisons to make
* 1 & 2 TM chains have a different set of rules from larger chains, including a condition on the size of the domains on either side of the membrane

#### Symmetry algorithms used

* CE-Symm v2.2.3
* QuatSymm v2.2.3
* SymD v1.6
* AnaNaS v1.1

#### Multi-step Symmetry Selection (MSSD)

* Integrated QuatSymm into the MSSD procedure. QuatSymm results are post-processed to guess the specific repeat range; the output is only used if the resulting symmetry has comparable RMSD and TM-score to the one reported by QuatSymm.
* Quaternary symmetries with only 1 TM chain in a repeat are now considered acceptable and are reported.

---

## Contributors

* Antoniya A. Aleksandrova
* Edoardo Sarti
* Lucy R. Forrest

If you use EncoMPASS in your work, please cite:

- Aleksandrova AA, Sarti E, Forrest LR. EncoMPASS: An encyclopedia of membrane proteins analyzed by structure and symmetry. *Structure* 32(4):492–504.e4 (2024). [https://doi.org/10.1016/j.str.2024.01.011](https://doi.org/10.1016/j.str.2024.01.011)

- Sarti E, Aleksandrova AA, Ganta SK, Yavatkar AS, Forrest LR. EncoMPASS: an online database for analyzing structure and symmetry in membrane proteins. *Nucleic Acids Research* 47(D1):D315–D321 (2019). [https://doi.org/10.1093/nar/gky952](https://doi.org/10.1093/nar/gky952)


