import os
import sys
import time
import argparse
import subprocess
import pickle as pkl

from encompass.pipeline.initialize_repository import initialize_repository
from encompass.sources.combine_sources import read_checkpoint
from encompass.symmetry.symmetry_exec_functions import cesymm_alignment, symd_alignment


def eprint(*args, **kwargs) -> None:
    """Print to stderr."""
    print(*args, file=sys.stderr, **kwargs)


def format_time(timestamp: float) -> str:
    """Return a human-readable timestamp string."""
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(timestamp))


def _split_semicolon_list(value: str | None) -> list[str]:
    """
    Turn a semicolon-separated string like 'a;b;c;' into ['a', 'b', 'c'].
    Empty entries and 'na' are dropped.
    """
    if not value:
        return []
    items = [v.strip() for v in value.split(";")]
    return [v for v in items if v and v.lower() != "na"]


def _build_prefixed_paths(
    names: list[str],
    web_prefix: str = "",
    fs_prefix: str | None = None,
    check_fs: bool = False,
) -> list[dict]:
    """
    Given file *names* and optional web / filesystem prefixes, build a list of:
        {"name": ..., "web_path": ..., "disk_path": ...}

    If check_fs is True, disk_path is only kept if the file actually exists.
    """
    paths: list[dict] = []
    for name in names:
        if not name:
            continue
        web_path = (web_prefix or "") + name if web_prefix is not None else None
        disk_path = None
        if fs_prefix:
            candidate = fs_prefix + name
            if (not check_fs) or os.path.isfile(candidate):
                disk_path = candidate
        paths.append(
            {
                "name": name,
                "web_path": web_path,
                "disk_path": disk_path,
            }
        )
    return paths


# ---------------------------------------------------------------------------
# Deletion & neighbor info
# ---------------------------------------------------------------------------

def prepare_deletion_info(locations: dict) -> tuple[list[tuple[str, str]], dict[str, str]]:
    """
    Read deletion codes/messages from system files and return
    (deletion_codes, deletion_messages).
    """
    deletion_codes: list[tuple[str, str]] = []
    print("Storing deletion codes from", locations["SYSFILES"]["delcodes"])
    with open(locations["SYSFILES"]["delcodes"]) as dcf:
        for line in dcf:
            if not line.strip():
                continue
            dc, c = line.split("\t")
            deletion_codes.append((dc, c.strip()))

    deletion_messages: dict[str, str] = {}
    print("Storing deletion messages from", locations["SYSFILES"]["delmsg"])
    with open(locations["SYSFILES"]["delmsg"]) as dmf:
        for line in dmf:
            if not line.strip():
                continue
            c, m = line.split("\t")
            deletion_messages[c.strip()] = m

    return deletion_codes, deletion_messages


def prepare_neighbor_info(locations: dict):
    """
    Compute neighbor counts (seq/str/total) from files on disk using wc -l.

    Returns:
        neighbor_counts: (seqneighs, strneighs, totneighs)
        filename_examples: (seqneighs_fn_example, strneighs_fn_example, totneighs_fn_example)
    """
    start_time = time.time()

    def _wc_counts(dir_path: str) -> tuple[dict[str, int], str]:
        cmd = (
            "for x in `ls {dir}`; do wc -l {dir}/${{x}}; done"
            .format(dir=dir_path)
        )
        out = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True,
        ).stdout

        counts: dict[str, int] = {}
        example_name: str | None = None

        if out is not None:
            for raw in out:
                parts = raw.decode("utf8").split()
                if not parts or "total" in parts[1]:
                    continue
                n = int(parts[0])
                fname = os.path.basename(parts[1])
                # filenames are like seqneigh_1abc_A.txt; last 10..4 is "1abc_A"
                counts[fname[-10:-4]] = n
                if example_name is None:
                    example_name = fname

        if example_name is None:
            example_name = ""

        return counts, example_name

    seqneighs, seqneighs_fn_templ = _wc_counts(locations["FSYSPATH"]["seqneighs"])
    print("obtained seqneigh neighbor counts; example", seqneighs_fn_templ)

    strneighs, strneighs_fn_templ = _wc_counts(locations["FSYSPATH"]["strneighs"])
    print("obtained strneigh neighbor counts; example", strneighs_fn_templ)

    totneighs, totneighs_fn_templ = _wc_counts(locations["FSYSPATH"]["totneighs"])
    print("obtained totneigh neighbor counts; example", totneighs_fn_templ)

    elapsed_time = time.time() - start_time
    print("Neighbor counting took {:.2f} seconds".format(elapsed_time))

    neighbor_counts = (seqneighs, strneighs, totneighs)
    filename_examples = (seqneighs_fn_templ, strneighs_fn_templ, totneighs_fn_templ)
    return neighbor_counts, filename_examples


# ---------------------------------------------------------------------------
# Symmetry helpers (consensus + alignments + MSSD summaries)
# ---------------------------------------------------------------------------

def consensus_symmetry(symm_dic: dict, prot_class: str, transfer_list=None) -> tuple[str, str]:
    """
    Compute consensus symmetry order as in the original mssd logic.

    Returns:
        (reported_order, consensus_score)
    """
    cs, mssd, transfer, reported = "na", "na", "na", "na"

    if "cesymm" in symm_dic:
        cs = symm_dic["cesymm"].get("symmetry_order", "na")
        reported = cs

    if transfer_list:
        # transfer_list is expected to be a list of dicts
        transfer = transfer_list[0].get("template_order", "na")
        reported = transfer

    if "selected" in symm_dic:
        sel = symm_dic["selected"]["source"].split(";")
        key = sel[0]
        mssd = symm_dic.get(key, {}).get("symmetry_order", "na")
        reported = mssd

    if prot_class == "beta":
        consensus = "1.0"
    elif transfer == "na":
        # 2-way consensus between cs and mssd
        consensus = [cs, mssd].count(cs)
        consensus = str(float("{:.1f}".format(consensus / 2.0)))
    else:
        # 3-way consensus between cs, mssd, transfer
        consensus = [cs, mssd, transfer].count(cs)
        if [cs, mssd, transfer].count(mssd) > consensus:
            consensus = [cs, mssd, transfer].count(mssd)
        consensus = str(float("{:.1f}".format(consensus / 3.0)))

    return reported, consensus


def get_cesymm_alignment(
    struct: str,
    xml_type: str,
    symm_dic_pdbi: dict,
    locations: dict,
    verbose: int = 1,
) -> str:
    """
    Return the CE-Symm alignment string for a structure or chain.

    struct    : "1abc" for whole, or "1abc_A" for chain
    xml_type  : "whole" or "chains"
    symm_dic_pdbi : entry from the symmetry results dict for this struct
    locations : EncoMPASS locations dict
    """
    if "cesymm" not in symm_dic_pdbi:
        return ""

    ces = symm_dic_pdbi["cesymm"]

    # 1) Preferred: alignment directly in the dict
    if ces.get("alignment"):
        if verbose:
            print("\t\tusing CE-Symm alignment from dictionary")
        return ces["alignment"]

    # 2) Backwards-compatibility: typo key
    if ces.get("alignmnet"):
        if verbose:
            print("\t\tusing CE-Symm alignment from 'alignmnet' key")
        return ces["alignmnet"]

    # 3) Fallback: try to reconstruct from FASTA on disk
    if xml_type == "whole":
        fasta_name = struct[0:4]
    else:  # "chains"
        fasta_name = struct

    fasta_root = (
        locations["FSYSPATH"]["main"]
        + locations["FSYS"]["cesymm"]
        + struct[0:4]
        + "/"
        + fasta_name
    )

    if verbose:
        print("\t\tNo CE-Symm alignment in dict, trying fasta:", fasta_root + ".fasta")

    if os.path.isfile(fasta_root + ".fasta"):
        aln = cesymm_alignment(
            locations["FSYSPATH"]["main"],
            locations["FSYS"]["cesymm"] + struct[0:4] + "/" + fasta_name,
        )
        if not aln and verbose:
            print("\t\tAlignment file exists but is empty of alignment:", fasta_root + ".fasta")
        return aln or ""
    else:
        if verbose:
            print("\tWARNING, expected CE-Symm alignment file doesn't exist:", fasta_root + ".fasta")
        return ""


def get_symd_alignment(
    struct: str,
    xml_type: str,
    symm_dic_pdbi: dict,
    locations: dict,
    verbose: int = 1,
) -> str:
    """
    Return the SymD alignment string for a structure or chain.
    """
    if "symd" not in symm_dic_pdbi:
        return ""

    sd = symm_dic_pdbi["symd"]

    # 1) Alignment in dictionary
    if sd.get("alignment"):
        if verbose:
            print("\t\tusing SymD alignment from dictionary")
        return sd["alignment"]

    if sd.get("alignmnet"):
        if verbose:
            print("\t\tusing SymD alignment from 'alignmnet' key")
        return sd["alignmnet"]

    # 2) Fallback: fasta-derived
    if xml_type == "whole":
        fasta_name = struct[0:4]
    else:
        fasta_name = struct

    fasta_root = (
        locations["FSYSPATH"]["main"]
        + locations["FSYS"]["symd"]
        + struct[0:4]
        + "/"
        + fasta_name
        + "-best"
    )

    if verbose:
        print("\t\tNo SymD alignment in dict, trying fasta:", fasta_root + ".fasta")

    if os.path.isfile(fasta_root + ".fasta"):
        aln = symd_alignment(
            locations["FSYSPATH"]["main"],
            locations["FSYS"]["symd"] + struct[0:4] + "/" + fasta_name,
        )
        if not aln and verbose:
            print("\t\tAlignment file exists but is empty of alignment:", fasta_root + ".fasta")
        return aln or ""
    else:
        if verbose:
            print("\tWARNING, expected SymD alignment file doesn't exist:", fasta_root + ".fasta")
        return ""


def get_mssd_alignments(symm_dic_pdbi: dict, locations: dict, verbose: int = 1) -> list[dict]:
    """
    Return a list of {method, alignment} dicts for the selected MSSD components.
    """
    if "selected" not in symm_dic_pdbi:
        return []

    sel = symm_dic_pdbi["selected"]["source"].strip(";").split(";")
    alignments: list[dict] = []

    for s in sel:
        if not s:
            continue

        entry = symm_dic_pdbi.get(s, {})
        if not entry:
            continue

        if "alignment" in entry:
            aln = entry["alignment"]
        elif "alignmnet" in entry:
            aln = entry["alignmnet"]
        else:
            # Fallback to cesymm_alignment using files_dir/files_key,
            # same as mssd_xml currently does.
            aln = cesymm_alignment(
                locations["FSYSPATH"]["main"] + entry["files_dir"],
                entry["files_key"],
            )
        alignments.append({"method": s, "alignment": aln})

    return alignments


# ---------------------------------------------------------------------------
# Symmetry file-path bundles (stored fully in web_data)
# ---------------------------------------------------------------------------

def build_cesymm_file_info(
    pdbi: str,
    monomer: bool,
    symm_dic_pdbi: dict,
    locations: dict,
) -> dict:
    """
    Build a CE-Symm file/paths bundle for a given structure.

    Returns a dict with:
        - symmetry numbers (order, repeats, etc.)
        - image / pml / jmol file names AND resolved paths
        - download file paths for stdout/xml/axes/fasta
    """
    if "cesymm" not in symm_dic_pdbi:
        return {}

    ce = symm_dic_pdbi["cesymm"]
    fsys = locations["FSYS"]
    fspath = locations["FSYSPATH"]

    # For whole-structure CE-Symm paths, XML uses 'symm_whole' vs 'symm_chains'
    webdir = "symm_whole" if monomer else "symm_chains"

    image_names = _split_semicolon_list(ce.get("image_files"))
    pml_names = _split_semicolon_list(ce.get("pml_files"))
    jmol_names = _split_semicolon_list(ce.get("jmol_files"))

    image_paths = _build_prefixed_paths(
        image_names,
        web_prefix=fsys.get(webdir + "pngs", ""),
        fs_prefix=fspath.get(webdir + "pngs"),
        check_fs=True,
    )
    pml_paths = _build_prefixed_paths(
        pml_names,
        web_prefix=fsys.get(webdir + "pngs", ""),
        fs_prefix=fspath.get(webdir + "pngs"),
        check_fs=True,
    )
    jmol_paths = _build_prefixed_paths(
        jmol_names,
        web_prefix=fsys.get(webdir + "jsons", ""),
        fs_prefix=fspath.get(webdir + "jsons"),
        check_fs=True,
    )

    # Download CE-Symm files (stdout/xml/axes/fasta), as in cesymm_xml()
    download_suffixes = ["_stdout.out", "_output.xml", ".axes", ".fasta"]
    downloads: list[dict] = []
    for suffix in download_suffixes:
        if monomer:
            disk = fspath["cesymm"] + pdbi[0:4] + "/" + pdbi[0:4] + suffix
        else:
            disk = fspath["cesymm"] + pdbi[0:4] + "/" + pdbi + suffix
        if os.path.isfile(disk):
            web = disk.replace(fspath["cesymm"], fsys["cesymm"])
            downloads.append(
                {
                    "suffix": suffix,
                    "disk_path": disk,
                    "web_path": web,
                }
            )

    return {
        "symmetry_order": ce.get("symmetry_order"),
        "repeats_number": ce.get("repeats_number"),
        "symmetry_levels": ce.get("symmetry_levels"),
        "repeat_length": ce.get("repeat_length"),
        "aligned_length": ce.get("aligned_length"),
        "refined_rmsd": ce.get("refined_rmsd"),
        "refined_tmscore": ce.get("refined_tmscore"),
        "unrefined_rmsd": ce.get("unrefined_rmsd"),
        "unrefined_tmscore": ce.get("unrefined_tmscore"),
        "seed": ce.get("seed"),
        "repeats": ce.get("repeats"),
        "image_files": image_names,
        "image_paths": image_paths,
        "pml_files": pml_names,
        "pml_paths": pml_paths,
        "jmol_files": jmol_names,
        "jmol_paths": jmol_paths,
        "downloads": downloads,
    }


def build_symd_file_info(
    pdbi: str,
    monomer: bool,
    symm_dic_pdbi: dict,
    locations: dict,
) -> dict:
    """
    Build SymD file/paths bundle. Pattern follows symd_xml().

    - Download: -info.txt, -trfm.pdb, -best.fasta (if they exist)
    - Image / pml / jmol names as stored; these are already "full" in the XML.
    """
    if "symd" not in symm_dic_pdbi:
        return {}

    sd = symm_dic_pdbi["symd"]
    fsys = locations["FSYS"]
    fspath = locations["FSYSPATH"]

    download_suffixes = ["-info.txt", "-trfm.pdb", "-best.fasta"]
    downloads: list[dict] = []
    for suffix in download_suffixes:
        if monomer:
            disk = fspath["symd"] + pdbi[0:4] + "/" + pdbi[0:4] + suffix
        else:
            disk = fspath["symd"] + pdbi[0:4] + "/" + pdbi + suffix
        if os.path.isfile(disk):
            web = disk.replace(fspath["symd"], fsys["symd"])
            downloads.append(
                {
                    "suffix": suffix,
                    "disk_path": disk,
                    "web_path": web,
                }
            )

    image_names = _split_semicolon_list(sd.get("image_files"))
    pml_names = _split_semicolon_list(sd.get("pml_files"))
    jmol_names = _split_semicolon_list(sd.get("jmol_files") or sd.get("jmol_symd"))

    return {
        "symmetry_order": sd.get("symmetry_order"),
        "image_files": image_names,
        "pml_files": pml_names,
        "jmol_files": jmol_names,
        "downloads": downloads,
    }


def build_analysis_file_info(
    pdbi: str,
    monomer: bool,
    symm_dic_pdbi: dict,
    locations: dict,
) -> list[dict]:
    """
    Build a list of "analysis_files" entries, one per selected method
    (cesymm, symd, quatsymm, etc.), matching the spirit of mssd_xml().

    Each entry includes:
        - method name
        - relative path from DB root (for downloads)
        - download paths (stdout/xml/axes/fasta)
        - per-method image / pml / jmol / super-pml names and resolved paths
    """
    if "selected" not in symm_dic_pdbi:
        return []

    sel = symm_dic_pdbi["selected"]["source"].strip(";").split(";")
    sel = [s for s in sel if s]

    fsys = locations["FSYS"]
    fspath = locations["FSYSPATH"]
    db_root_fs = fspath["main"].rstrip("/")
    db_name = os.path.basename(db_root_fs)

    results: list[dict] = []

    for method in sel:
        entry = symm_dic_pdbi.get(method, {})
        files_key = entry.get("files_key", "")

        files_dir_path = entry.get("files_dir_path") or entry.get("files_dir") or ""
        rel_path = ""
        disk_prefix = None
        web_prefix = None

        if files_dir_path:
            # In mssd_xml(), rel_path is computed by splitting on "<DBNAME>/"
            if f"{db_name}/" in files_dir_path:
                rel_path = files_dir_path.split(f"{db_name}/", 1)[1]
                if not rel_path.endswith("/"):
                    rel_path += "/"
                disk_prefix = db_root_fs + "/" + rel_path
                web_prefix = rel_path  # relative to DB root
            else:
                disk_prefix = files_dir_path
                web_prefix = ""

        downloads: list[dict] = []
        for suffix in ["_stdout.out", "_output.xml", ".axes", ".fasta"]:
            if disk_prefix and files_key:
                disk = os.path.join(disk_prefix, files_key + suffix)
                if os.path.isfile(disk):
                    web_path = (web_prefix or "") + files_key + suffix
                    downloads.append(
                        {
                            "suffix": suffix,
                            "disk_path": disk,
                            "web_path": web_path,
                        }
                    )

        # Per-method images (as in mssd_xml); webdir depends on method & monomer
        xml_type = "whole" if monomer else "chains"
        if method == "cesymm":
            webdir = "symm_" + xml_type
        else:
            webdir = "analysis_" + xml_type

        image_names = _split_semicolon_list(entry.get("image_files"))
        pml_names = _split_semicolon_list(entry.get("pml_files"))
        jmol_names = _split_semicolon_list(entry.get("jmol_files"))
        super_pml_names = _split_semicolon_list(entry.get("super_pml_files"))

        image_paths = _build_prefixed_paths(
            image_names,
            web_prefix=fsys.get(webdir + "pngs", ""),
            fs_prefix=fspath.get(webdir + "pngs"),
            check_fs=True,
        )
        pml_paths = _build_prefixed_paths(
            pml_names,
            web_prefix=fsys.get(webdir + "pngs", ""),
            fs_prefix=fspath.get(webdir + "pngs"),
            check_fs=True,
        )
        t3dmol_paths = _build_prefixed_paths(
            jmol_names,
            web_prefix=fsys.get(webdir + "jsons", ""),
            fs_prefix=fspath.get(webdir + "jsons"),
            check_fs=True,
        )
        super_pml_paths = _build_prefixed_paths(
            super_pml_names,
            web_prefix=fsys.get(f"analysis_{xml_type}super", ""),
            fs_prefix=fspath.get(f"analysis_{xml_type}super"),
            check_fs=True,
        )

        results.append(
            {
                "method": method,
                "files_key": files_key,
                "rel_path": rel_path,
                "downloads": downloads,
                "image_files": image_names,
                "image_paths": image_paths,
                "pml_files": pml_names,
                "pml_paths": pml_paths,
                "jmol_files": jmol_names,
                "jmol_paths": t3dmol_paths,
                "super_pml_files": super_pml_names,
                "super_pml_paths": super_pml_paths,
            }
        )

    return results


def build_transfer_file_info(
    inferred_dic_list: list[dict],
    locations: dict,
    xml_type: str,
) -> dict:
    """
    Build a transfer-symmetry bundle from the inferred_dic_list.

    Captures:
        - Download paths for axes / alnres (relative + full disk)
        - Image / pml / 3dmol / super-pml names AND resolved web/disk paths
    """
    if not inferred_dic_list:
        return {}

    fsys = locations["FSYS"]
    fspath = locations["FSYSPATH"]

    webdir = "transfer_" + xml_type

    # DownloadTransfer: axes_file and alnres_file (paths relative to DB root)
    downloads: list[dict] = []
    for entry in inferred_dic_list:
        for key in ("axes_file", "alnres_file"):
            rel_path = entry.get(key)
            if not rel_path:
                continue
            disk_path = os.path.join(fspath["main"], rel_path)
            if os.path.isfile(disk_path):
                downloads.append(
                    {
                        "kind": key,
                        "rel_path": rel_path,
                        "disk_path": disk_path,
                    }
                )

    # TransferImage: aggregate names from all entries, then build full paths
    image_names = _split_semicolon_list(
        ";".join(e.get("image_files", "") for e in inferred_dic_list)
    )
    pml_names = _split_semicolon_list(
        ";".join(e.get("pml_files", "") for e in inferred_dic_list)
    )
    jmol_names = _split_semicolon_list(
        ";".join(e.get("jmol_files", "") for e in inferred_dic_list)
    )
    super_pml_names = _split_semicolon_list(
        ";".join(e.get("super_pml_files", "") for e in inferred_dic_list)
    )

    image_paths = _build_prefixed_paths(
        image_names,
        web_prefix=fsys.get(webdir + "pngs", ""),
        fs_prefix=fspath.get(webdir + "pngs"),
        check_fs=True,
    )
    pml_paths = _build_prefixed_paths(
        pml_names,
        web_prefix=fsys.get(webdir + "pngs", ""),
        fs_prefix=fspath.get(webdir + "pngs"),
        check_fs=True,
    )
    t3dmol_paths = _build_prefixed_paths(
        jmol_names,
        web_prefix=fsys.get(webdir + "jsons", ""),
        fs_prefix=fspath.get(webdir + "jsons"),
        check_fs=True,
    )
    super_pml_paths = _build_prefixed_paths(
        super_pml_names,
        web_prefix=fsys.get(webdir + "super", ""),
        fs_prefix=fspath.get(webdir + "super"),
        check_fs=True,
    )

    return {
        "downloads": downloads,
        "image_files": image_names,
        "image_paths": image_paths,
        "pml_files": pml_names,
        "pml_paths": pml_paths,
        "jmol_files": jmol_names,
        "jmol_paths": t3dmol_paths,
        "super_pml_files": super_pml_names,
        "super_pml_paths": super_pml_paths,
    }


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_str_data(locations: dict) -> dict:
    """
    Load str_data dict from cache/template using read_checkpoint.
    """
    str_data_filename = locations["FSYSPATH"]["cache"] + "str_data_completegen.pkl"
    str_template_filename = locations["SYSFILES"]["data_structure_template"]
    print("Reading", str_data_filename)
    print("Reading", str_template_filename)
    str_data = read_checkpoint(str_data_filename, str_template_filename)
    print("Finished reading str_data at", format_time(time.time()))
    return str_data


def load_mssd_data(locations: dict) -> dict:
    """
    Load mssd_data (symmetry_results.pkl).
    """
    mssd_data_filename = locations["FSYSPATH"]["mssd"] + "symmetry_results.pkl"
    print("Reading", mssd_data_filename)
    mssd_data = pkl.load(open(mssd_data_filename, "rb"))
    print("Finished reading mssd_data at", format_time(time.time()))
    return mssd_data


def load_neighbor_counts(locations: dict):
    """
    Load or compute neighbor counts.

    Returns:
        neighbor_counts: (seqneighs, strneighs, totneighs)
        neighbor_example_fns: (seqneighs_fn_example, strneighs_fn_example, totneighs_fn_example)
    """
    neighbor_counts_filename = locations["FSYSPATH"]["cache"] + "neighborcounts.pkl"
    print("Preparing neighbor counts")
    if os.path.exists(neighbor_counts_filename):
        print("WARNING: Checkpoint neighbor count file found, reading", neighbor_counts_filename)
        with open(neighbor_counts_filename, "rb") as f:
            (all_counts, example_filenames) = pkl.load(f)
        print("\tloaded", neighbor_counts_filename)
        return all_counts, example_filenames
    else:
        (seqneighs, strneighs, totneighs), (seqneighs_fn_example,
                                            strneighs_fn_example,
                                            totneighs_fn_example) = prepare_neighbor_info(locations)
        with open(neighbor_counts_filename, "wb") as f:
            pkl.dump(
                ((seqneighs, strneighs, totneighs),
                 (seqneighs_fn_example, strneighs_fn_example, totneighs_fn_example)),
                f,
            )
            print("Checkpoint neighbor count data saved:", neighbor_counts_filename)
        return (seqneighs, strneighs, totneighs), (seqneighs_fn_example,
                                                   strneighs_fn_example,
                                                   totneighs_fn_example)


# ---------------------------------------------------------------------------
# Path builders (structure / chain) – these store *final paths* in web_data
# ---------------------------------------------------------------------------

def build_structure_paths(pdbi: str, locations: dict) -> dict:
    """Return file/url paths for the whole structure (non-symmetry)."""
    fsys = locations["FSYS"]
    return {
        "files": {
            "pdb_file":   f"{fsys['whole']}{pdbi}_enc.pdb",
            "image_whole": f"{fsys['images_figs_whole']}{pdbi}.png",
            "t3dmol_json": f"{fsys['3dmol_whole']}{pdbi}.json",
        },
        "external_urls": {
            "pdb":   f"http://www.rcsb.org/pdb/explore/explore.do?structureId={pdbi.upper()}",
            "opm":   f"http://opm.phar.umich.edu/protein.php?pdbid={pdbi}",
            "pdbtm": f"http://pdbtm.enzim.hu/data/database/{pdbi[1:3]}/{pdbi}.trpdb.gz",
        },
    }


def build_chain_paths(pdbi: str, tmch: str, str_data_pdbi: dict, locations: dict) -> dict:
    """Return file/url paths for a specific chain (non-symmetry)."""
    fsys = locations["FSYS"]
    pdbi_ch = f"{pdbi}_{tmch}"

    struct = str_data_pdbi["ENCOMPASS"]["structure"]
    redund = struct.get("redundant_chains", {})

    if tmch in redund:
        repr_ch = redund[tmch]
        if repr_ch in redund:
            repr_ch = redund[repr_ch]
    else:
        repr_ch = tmch

    return {
        "files": {
            "pdb_file":            f"{fsys['whole']}{pdbi}_enc.pdb",
            "image_chain":         f"{fsys['images_figs_chain']}{pdbi_ch}.png",
            "t3dmol_json":         f"{fsys['3dmol_whole']}{pdbi}.json",
            "seq_neighbors_file":  f"{fsys['seqneighs']}seqneigh_{pdbi_ch}.txt",
            "str_neighbors_file":  f"{fsys['strneighs']}strneigh_{pdbi_ch}.txt",
            "tot_neighbors_file":  f"{fsys['totneighs']}totneigh_{pdbi_ch}.txt",
            "resdistance_image":   f"{fsys['distributions_figs']}distr_{pdbi_ch}.png",
            "resdistance_data":    f"{fsys['distributions_data']}distr_{pdbi}_{repr_ch}.txt",
            "structure_alignments": f"{fsys['stralns']}stralns_{pdbi}_{repr_ch}.txt.gz",
            "estimators":          f"{fsys['estimators']}est_{pdbi}_{repr_ch}.txt",
        },
        "plots": {
            "radial_distribution_image": f"{fsys['polar_figs']}p_{pdbi_ch}.png",
            "radial_distribution_map":   f"{fsys['polar_maps']}{pdbi_ch}_cmap.txt",
            "density_plot_image":        f"{fsys['densityscatter_figs']}ds_{pdbi_ch}.png",
            "density_plot_map":          f"{fsys['densityscatter_maps']}{pdbi_ch}_cmap.txt",
            "topology_image":            f"{fsys['topologies_figs']}top_{pdbi_ch}.jpeg",
        },
        "external_urls": {
            "pdb":   f"http://www.rcsb.org/pdb/explore/explore.do?structureId={pdbi.upper()}",
            "opm":   f"http://opm.phar.umich.edu/protein.php?pdbid={pdbi}",
            "pdbtm": f"http://pdbtm.enzim.hu/data/database/{pdbi[1:3]}/{pdbi}.trpdb.gz",
        },
    }


# ---------------------------------------------------------------------------
# Main web_data builder
# ---------------------------------------------------------------------------

def build_web_data(
    locations: dict,
    str_data: dict,
    symm_data: dict,
    neighbor_counts: tuple[dict, dict, dict],
    deletion_info: tuple[list, dict],
    verbose: int = 1,
) -> dict:
    """
    Aggregate all information needed for the website (and for XML / SQL export)
    into a single, self-contained ``web_data`` dictionary.

    The goal is that *all* paths stored here are **relative** paths (web paths)
    that can be written directly into XML or SQL, and that we do **not** depend
    on any external location dictionaries when rendering later.

    High-level shape (per PDB code ``pdbi``)::

        web_data[pdbi] = {
            "pdb_code": "1a0s",
            "title": "...",              # ENCOMPASS name
            "class": "alpha" | "beta" | ...,
            "monomeric": True/False,

            "structure": {
                "method": "...",
                "resolution": 2.4,
                "size": 413,
                "n_chains": 1,
                "n_tmchains": 1,
                "chain_uniprot_accs": ["P22340", ...],
                "general_structure": {
                    "pdb_file": "database/selection/whole_structs/1a0s_enc.pdb",
                    "image_file": "analysis/analysis/images/images_figs/chains/1a0s_P.png",
                    "t3dmol_file": "webarchive/whole_structs_json/1a0s.json",
                },
                "urls": {
                    "pdb": "...",
                    "opm": "...",
                    "pdbtm": "...",
                },
            },

            "neighbors": {  # structure-level summary (optional)
                "sequence_neighbors": int | None,
                "structure_neighbors": int | None,
                "total_neighbors": int | None,
            },

            "deletion": {
                "code": str | None,
                "message": str | None,
            },

            "symmetry": {
                "cesymm_alignment_whole": str,   # <Aln>…</Aln> block, or ""
                "symd_alignment_whole": str,
                "analysis_alignments": [         # MSSD-style alignments, if any
                    {"method": "cesymm", "alignment": "..."},
                    {"method": "symd", "alignment": "..."},
                    ...
                ],
                "consensus_by_chain": {
                    "1a0s_P": {
                        "reported_order": "C1",
                        "consensus_score": 1.0,
                        "symmetry_order_cesymm": "C1" | "na",
                        "symmetry_order_mssd": "C12" | "na",
                        "symmetry_order_transfer": "C3" | "na",
                    },
                    ...
                },
                # Optional per-PDB transfer-symmetry information (when available)
                "transfer": {
                    "image_files": [ "webarchive/transfer_chains_pngs/1a0s_P_transfer_1.png", ... ],
                    "pml_files":   [ "webarchive/transfer_chains_pngs/pymol_script_1a0s_P_transfer_1.pml", ... ],
                    "jmol_files":  [ "webarchive/transfer_chains_json/3dmol_1a0s_P_transfer_1.json", ... ],
                    "super_pml_files": [ "webarchive/transfer_superposition_pml/superposition_1a0s_P_transfer_1.pml", ... ],
                },
            },

            "chains": {
                "P": {
                    "pdb_chain_id": "1a0s_P",
                    "size": 413,
                    "sequence": "SGFEFH...",
                    "tm_domains": [
                        (75, 87),
                        (105, 124),
                        ...
                    ],
                    "uniprot_accs": ["P22340", ...],

                    "neighbors": {
                        "sequence_neighbors": 8,
                        "structure_neighbors": 393,
                        "total_neighbors": 393,
                    },
                    "neighbor_files": {
                        "sequence_neighbors_file": "analysis/analysis/neighbor_lists/seq_neighbors/seqneigh_1a0s_P.txt",
                        "structure_neighbors_file": "analysis/analysis/neighbor_lists/str_neighbors/strneigh_1a0s_P.txt",
                        "total_neighbors_file": "analysis/analysis/neighbor_lists/tot_neighbors/totneigh_1a0s_P.txt",
                    },
                    "images": {
                        "radial_distribution": "analysis/analysis/polarplots/polarplots_figs/p_1a0s_P.png",
                        "radial_distribution_map": "analysis/analysis/polarplots/polarplots_maps/1a0s_P_cmap.txt",
                        "density_plot": "analysis/analysis/densityscatter/densityscatter_figs/ds_1a0s_P.png",
                        "density_plot_map": "analysis/analysis/densityscatter_maps/1a0s_P_cmap.txt",
                        "res_distance_distribution": "analysis/analysis/distributions/distributions_figs/distr_1a0s_P.png",
                        "topology": "analysis/analysis/topologies/topologies_figs/top_1a0s_P.jpeg",
                    },
                    "additional_files": {
                        "res_distance_distribution": "analysis/analysis/distributions/distributions_data/distr_1a0s_P.txt",
                        "structure_alignments": "database/alignments/stralns/stralns_1a0s_P.txt.gz",
                        "estimators": "analysis/analysis/estimators/est_1a0s_P.txt",
                    },

                    "symmetry_consensus": { ... },
                    "cesymm_alignment": str,
                    "symd_alignment": str,
                },
                ...
            },
        }
    """

    (seqneighs, strneighs, totneighs) = neighbor_counts
    deletion_codes, deletion_messages = deletion_info  # deletion_codes currently unused

    web_data: dict = {}

    for pdbi, str_data_pdbi in str_data.items():
        if verbose:
            print(pdbi)

        # --- choose symmetry dictionary for this PDB (beta barrels may have their own pickle) ---
        symm_results = symm_data
        if str_data_pdbi["ENCOMPASS"]["class"] == "beta":
            beta_dictionary_filename = (
                locations["FSYSPATH"]["images_cesymm_symd_log"]
                + pdbi
                + "/"
                + pdbi
                + "_dic.pkl"
            )
            if os.path.isfile(beta_dictionary_filename):
                if verbose:
                    print("\tReading beta barrel pickle for", pdbi, beta_dictionary_filename)
                with open(beta_dictionary_filename, "rb") as f:
                    symm_results = pkl.load(f)

        # --- monomer vs multimer ---
        num_tm_ch = len(
            [c for c in str_data_pdbi["ENCOMPASS"]["structure"]["ktmchains"] if c != "-"]
        )
        monomer = num_tm_ch == 1
        if verbose:
            print("\n", pdbi, num_tm_ch, "chains", monomer)

        # --- structure-level metadata ---

        enc_struct = str_data_pdbi["ENCOMPASS"]["structure"]
        enc_meta = str_data_pdbi["ENCOMPASS"]

        title = enc_meta.get("name", "")
        prot_class = enc_meta.get("class", "")
        resolution = enc_meta.get("resolution", None)
        if resolution in ("None", "", "na", None):
            resolution_val = None
        else:
            try:
                resolution_val = float(resolution)
            except Exception:
                resolution_val = None

        size = enc_meta.get("n_of_aas", None)
        n_chains = enc_meta.get("n_of_chains", None)
        n_tmchains = enc_meta.get("n_of_tmchains", None)

        # structure-level UniProt stoichiometry
        stoich = enc_struct.get("UniProt_stoichiometry", [])

        # Normalize to a list
        if isinstance(stoich, str):
            codes = [stoich]
        elif isinstance(stoich, (list, tuple)):
            codes = list(stoich)
        else:
            codes = []

        chain_uniprot_accs = [c for c in codes if isinstance(c, str) and c.strip()]


        # basic URLs (as in original XML)
        pdb_url = f"http://www.rcsb.org/pdb/explore/explore.do?structureId={pdbi.upper()}"
        opm_url = f"http://opm.phar.umich.edu/protein.php?pdbid={pdbi}"
        pdbtm_url = f"http://pdbtm.enzim.hu/data/database/{pdbi[1:3]}/{pdbi}.trpdb.gz"

        # structure-level files (whole structure)
        struct_paths = {
            "pdb_file": f"{locations['FSYS']['whole']}{pdbi}_enc.pdb",
            "image_file": f"{locations['FSYS']['images_figs_whole']}{pdbi}.png",
            "t3dmol_file": f"{locations['FSYS']['3dmol_whole']}{pdbi}.json",
        }

        # --- base container for this structure ---

        web_struct = {
            "pdb_code": pdbi,
            "title": title,
            "class": prot_class,
            "monomeric": monomer,
            "structure": {
                "method": str_data_pdbi.get("FROM_PDB", {}).get("experimental_method", ""),
                "resolution": resolution_val,
                "size": size,
                "n_chains": n_chains,
                "n_tmchains": n_tmchains,
                "chain_uniprot_accs": chain_uniprot_accs,
                "general_structure": struct_paths,
                "urls": {
                    "pdb": pdb_url,
                    "opm": opm_url,
                    "pdbtm": pdbtm_url,
                },
            },
            "chains": {},
            "symmetry": {
                "consensus_by_chain": {},
                "cesymm_alignment_whole": "",
                "symd_alignment_whole": "",
                "analysis_alignments": [],
                "transfer": {
                    "image_files": [],
                    "pml_files": [],
                    "jmol_files": [],
                    "super_pml_files": [],
                },
            },
            "neighbors": {
                # you *can* aggregate to structure-level if desired; for now leave None
                "sequence_neighbors": None,
                "structure_neighbors": None,
                "total_neighbors": None,
            },
            "deletion": {},
        }

        # --- deletion information (as in original create_xml) ---
        del_key = str_data_pdbi.get("delete_keyword")
        if del_key:
            web_struct["deletion"]["code"] = del_key
            web_struct["deletion"]["message"] = deletion_messages.get(del_key, "").replace('"', "").strip()
        else:
            web_struct["deletion"]["code"] = None
            web_struct["deletion"]["message"] = None

        # --- whole-structure symmetry (for monomers) ---

        if pdbi in symm_results:
            if verbose:
                print("\tSymmetry results for", pdbi, "found in dictionary")
            symm_dic_pdbi = symm_results[pdbi]
            web_struct["symmetry"]["cesymm_alignment_whole"] = get_cesymm_alignment(
                struct=pdbi,
                xml_type="whole",
                symm_dic_pdbi=symm_dic_pdbi,
                locations=locations,
                verbose=verbose,
            )
            web_struct["symmetry"]["symd_alignment_whole"] = get_symd_alignment(
                struct=pdbi,
                xml_type="whole",
                symm_dic_pdbi=symm_dic_pdbi,
                locations=locations,
                verbose=verbose,
            )

            # MSSD-style analysis (if present)
            analysis_alignments = get_mssd_alignments(
                symm_dic_pdbi=symm_dic_pdbi,
                locations=locations,
                verbose=verbose,
            )
            web_struct["symmetry"]["analysis_alignments"] = analysis_alignments

            # Transfer-symmetry paths (if stored in symm_dic_pdbi)
            transfer_list = symm_dic_pdbi.get("transfer_list", [])
            if transfer_list:
                img_paths = []
                pml_paths = []
                jmol_paths = []
                super_pml_paths = []
                for entry in transfer_list:
                    # these should already be *relative* paths
                    img_paths.extend(entry.get("image_files", []))
                    pml_paths.extend(entry.get("pml_files", []))
                    jmol_paths.extend(entry.get("jmol_files", []))
                    super_pml_paths.extend(entry.get("super_pml_files", []))

                web_struct["symmetry"]["transfer"]["image_files"] = img_paths
                web_struct["symmetry"]["transfer"]["pml_files"] = pml_paths
                web_struct["symmetry"]["transfer"]["jmol_files"] = jmol_paths
                web_struct["symmetry"]["transfer"]["super_pml_files"] = super_pml_paths

        # --- per-chain aggregation ---

        consensus_by_chain: dict = {}
        redundant_chains = enc_struct.get("redundant_chains", {})

        for ich, tmch in enumerate(enc_struct["ktmchains"]):
            if tmch == "-":
                continue

            pdbi_ch = f"{pdbi}_{tmch}"
            chain_info = enc_struct["chains"][ich]

            # transfer dictionary (per chain), if present
            transf_dict_path = (
                locations["FSYSPATH"]["transfer"] + pdbi_ch + "/" + pdbi_ch + "_transfer.pkl"
            )
            if os.path.isfile(transf_dict_path):
                if verbose:
                    print("\tReading transfer pickle for", pdbi_ch, transf_dict_path)
                transfer_dic = pkl.load(open(transf_dict_path, "rb"))[pdbi_ch]
            else:
                if verbose:
                    print("\tNo transfer data file available for", pdbi_ch)
                transfer_dic = None

            # choose symmetry entry for this chain:
            #   prefer explicit chain entry; otherwise for monomer fall back to whole-structure entry
            if pdbi_ch in symm_results:
                chain_symm_entry = symm_results[pdbi_ch]
            elif monomer and pdbi in symm_results:
                chain_symm_entry = symm_results[pdbi]
            else:
                chain_symm_entry = {}

            # chain-level consensus (as in original code)
            chain_consensus = {}
            if chain_symm_entry:
                reported_order, consensus_score = consensus_symmetry(
                    chain_symm_entry,
                    str_data_pdbi["ENCOMPASS"]["class"],
                    transfer_list=transfer_dic,
                )
                chain_consensus = {
                    "reported_order": reported_order,
                    "consensus_score": consensus_score,
                }

            consensus_by_chain[pdbi_ch] = chain_consensus

            # chain-level alignments
            cesymm_alignment_str = ""
            symd_alignment_str = ""
            if chain_symm_entry:
                if "cesymm" in chain_symm_entry:
                    cesymm_alignment_str = get_cesymm_alignment(
                        struct=pdbi_ch,
                        xml_type="chains",
                        symm_dic_pdbi=chain_symm_entry,
                        locations=locations,
                        verbose=verbose,
                    )
                if "symd" in chain_symm_entry:
                    symd_alignment_str = get_symd_alignment(
                        struct=pdbi_ch,
                        xml_type="chains",
                        symm_dic_pdbi=chain_symm_entry,
                        locations=locations,
                        verbose=verbose,
                    )

            # --- chain-level neighbors, images, additional files ---

            seq_n = seqneighs.get(pdbi_ch, None)
            str_n = strneighs.get(pdbi_ch, None)
            tot_n = totneighs.get(pdbi_ch, None)

            neighbor_files = {
                "sequence_neighbors_file": f"{locations['FSYS']['seqneighs']}seqneigh_{pdbi_ch}.txt",
                "structure_neighbors_file": f"{locations['FSYS']['strneighs']}strneigh_{pdbi_ch}.txt",
                "total_neighbors_file": f"{locations['FSYS']['totneighs']}totneigh_{pdbi_ch}.txt",
            }

            images = {
                "radial_distribution": f"{locations['FSYS']['polar_figs']}p_{pdbi_ch}.png",
                "radial_distribution_map": f"{locations['FSYS']['polar_maps']}{pdbi_ch}_cmap.txt",
                "density_plot": f"{locations['FSYS']['densityscatter_figs']}ds_{pdbi_ch}.png",
                "density_plot_map": f"{locations['FSYS']['densityscatter_maps']}{pdbi_ch}_cmap.txt",
                "res_distance_distribution": f"{locations['FSYS']['distributions_figs']}distr_{pdbi_ch}.png",
                "topology": f"{locations['FSYS']['topologies_figs']}top_{pdbi_ch}.jpeg",
            }

            # representative chain for files stored once per redundancy group
            ch = tmch
            if ch in redundant_chains:
                repr_ch = redundant_chains[ch]
                if repr_ch in redundant_chains:
                    repr_ch = redundant_chains[repr_ch]
            else:
                repr_ch = ch

            additional_files = {
                "res_distance_distribution": f"{locations['FSYS']['distributions_data']}distr_{pdbi}_{repr_ch}.txt",
                "structure_alignments": f"{locations['FSYS']['stralns']}stralns_{pdbi}_{repr_ch}.txt.gz",
                "estimators": f"{locations['FSYS']['estimators']}est_{pdbi}_{repr_ch}.txt",
            }

            # sequence and TM domains from str_data
            sequence = chain_info.get("sequence", "")
            tm_domains = []
            tm_extrema = chain_info.get("TM_regions", {}).get("TM_regions_extrema", [])
            for tm in tm_extrema:
                try:
                    start = tm[0][0][0]
                    end = tm[1][0][0]
                    tm_domains.append((start, end))
                except Exception:
                    if isinstance(tm, (list, tuple)) and len(tm) >= 2:
                        tm_domains.append((tm[0], tm[1]))

            # chain-level UniProt accession(s)
            chain_uniprot = chain_info.get("UniProt_acc")
            if isinstance(chain_uniprot, str):
                chain_uniprot_accs = [chain_uniprot] if chain_uniprot else []
            elif isinstance(chain_uniprot, (list, tuple)):
                chain_uniprot_accs = [u for u in chain_uniprot if u]
            else:
                chain_uniprot_accs = []

            chain_entry = {
                "pdb_chain_id": pdbi_ch,
                "size": chain_info.get("n_of_residues"),
                "sequence": sequence,
                "tm_domains": tm_domains,
                "uniprot_accs": chain_uniprot_accs,
                "neighbors": {
                    "sequence_neighbors": seq_n,
                    "structure_neighbors": str_n,
                    "total_neighbors": tot_n,
                },
                "neighbor_files": neighbor_files,
                "images": images,
                "additional_files": additional_files,
                "symmetry_consensus": chain_consensus,
                "cesymm_alignment": cesymm_alignment_str,
                "symd_alignment": symd_alignment_str,
            }

            web_struct["chains"][tmch] = chain_entry

        # store consensus map at structure level
        web_struct["symmetry"]["consensus_by_chain"] = consensus_by_chain

        # save into web_data
        web_data[pdbi] = web_struct

    return web_data


# ---------------------------------------------------------------------------
# High-level pipeline + CLI
# ---------------------------------------------------------------------------

def build_web_data_from_repository(
    main_path: str = "",
    instr_filename_path: str = "",
    verbose: int = 1,
):
    """
    Full pipeline:
    1. initialize_repository(...) → options, locations
    2. load_str_data(locations)
    3. load_mssd_data(locations)
    4. load_neighbor_counts(locations)
    5. prepare_deletion_info(locations)
    6. build_web_data(...)
    """
    options, locations = initialize_repository(main_path, instr_filename_path)
    str_data = load_str_data(locations)
    mssd_data = load_mssd_data(locations)
    neighbor_counts, neighbor_examples = load_neighbor_counts(locations)
    deletion_info = prepare_deletion_info(locations)

    web_data = build_web_data(
        locations, str_data, mssd_data, neighbor_counts, deletion_info, verbose=verbose
    )
    return web_data, (options, locations)


def main(argv: list[str] | None = None) -> int:
    """
    Simple CLI entrypoint to build the web_data structure from an EncoMPASS repository.

    Example:
        python -m encompass.site_db.create_data_struct \\
            --main-path /path/to/EncoMPASS/db \\
            --instr-file EncoMPASS_options_relative_to_main.txt \\
            --output web_data.pkl
    """
    parser = argparse.ArgumentParser(
        description="Build aggregated web/DB data structure from an EncoMPASS repository."
    )
    parser.add_argument(
        "--main-path",
        "-m",
        dest="main_path",
        default="",
        help="Main EncoMPASS repository path (if empty, initialize_repository will infer it).",
    )
    parser.add_argument(
        "--instr-file",
        "-i",
        dest="instr_filename_path",
        default="",
        help=(
            "Instruction file name or path. If a bare filename is given and main-path is set, "
            "it is interpreted relative to main-path. If empty, initialize_repository() "
            "will use its default behavior."
        ),
    )
    parser.add_argument(
        "--output",
        "-o",
        dest="output",
        default="web_data.pkl",
        help="Output pickle file to write web_data to (default: web_data.pkl).",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        dest="verbose",
        type=int,
        default=1,
        help="Verbosity level (0 = quiet, 1 = normal, >1 = extra logging).",
    )

    args = parser.parse_args(argv)

    print("Building web_data from repository...")
    print("  main_path          :", args.main_path or "(empty / use initialize_repository default)")
    print("  instr_filename_path:", args.instr_filename_path or "(empty / use initialize_repository default)")
    print("  output             :", args.output)
    print("  verbose            :", args.verbose)

    web_data, (options, locations) = build_web_data_from_repository(
        main_path=args.main_path,
        instr_filename_path=args.instr_filename_path,
        verbose=args.verbose,
    )

    # Save result
    with open(args.output, "wb") as f:
        pkl.dump(web_data, f)
    print(f"web_data written to {args.output}")

    # Small summary for sanity-checking
    n_structs = len(web_data)
    print(f"Number of structures in web_data: {n_structs}")
    if n_structs:
        example_keys = list(web_data.keys())[:5]
        print("Example PDB codes:", ", ".join(example_keys))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
