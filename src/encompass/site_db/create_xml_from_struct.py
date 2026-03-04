from __future__ import annotations

import argparse
import os
import pickle as pkl
import xml.etree.ElementTree as ET
from typing import Any, Dict, List


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def indent(elem: ET.Element, level: int = 0) -> None:
    """
    Pretty-print helper for xml.etree.ElementTree.
    Mutates the tree in-place to add indentation whitespace.
    """
    i = "\n" + level * "\t"
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "\t"
        for child in elem:
            indent(child, level + 1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


def insert_alignment_xml(parent: ET.Element, tag_name: str, alignment_text: str) -> None:
    """
    Insert an alignment snippet (which already contains <Aln> tags etc.)
    as real XML children under a parent element <tag_name>.

    If parsing fails, we fall back to putting the raw text inside the tag.
    """
    alignment_text = (alignment_text or "").strip()
    elem = ET.SubElement(parent, tag_name)

    if not alignment_text:
        # empty alignment → nothing more to do
        return

    # Wrap the snippet so it's valid XML
    wrapped = f"<{tag_name}>{alignment_text}</{tag_name}>"
    try:
        temp_root = ET.fromstring(wrapped)
    except ET.ParseError:
        # Fallback: just text
        elem.text = alignment_text
        return

    # Move children from temp_root into elem
    for child in temp_root:
        elem.append(child)


def set_chain_uniprot_attrs(elem: ET.Element, accs: List[str]) -> None:
    """
    Populate something like:
      <ChainUniProtAccs code1="P22340" code2="Q99999" />
    from a list of accessions.
    """
    non_empty = [a for a in accs if a]
    for idx, acc in enumerate(non_empty, start=1):
        elem.set(f"code{idx}", acc)


# ---------------------------------------------------------------------------
# Whole-structure XML
# ---------------------------------------------------------------------------

def build_whole_structure_xml(root: ET.Element, pdb_code: str, entry: Dict[str, Any]) -> None:
    """
    Create the <Structure> element for a whole structure (1 per PDB code)
    and attach it to the 'whole' XML root, using the new web_data schema.
    """
    struct_elem = ET.SubElement(root, "Structure", ID=pdb_code)

    title_text = entry.get("title", "") or ""
    prot_class = entry.get("class", "") or ""
    structure = entry.get("structure", {}) or {}
    chains = entry.get("chains", {}) or {}

    # ----- Header -----
    header = ET.SubElement(struct_elem, "Header")
    title = ET.SubElement(header, "Title")
    title.text = title_text

    # ----- General -----
    general = ET.SubElement(struct_elem, "General")
    ET.SubElement(general, "PDBCode").text = pdb_code
    ET.SubElement(general, "Class").text = prot_class

    method = structure.get("method", "") or ""
    ET.SubElement(general, "Method").text = method

    size_elem = ET.SubElement(general, "Size")
    size_val = structure.get("size")
    if size_val is not None:
        size_elem.set("value", str(size_val))
    else:
        size_elem.set("value", "0")

    # Structure-level UniProt accs
    chain_uniprot_accs = structure.get("chain_uniprot_accs") or []
    cu_elem = ET.SubElement(general, "ChainUniProtAccs")
    set_chain_uniprot_attrs(cu_elem, chain_uniprot_accs)

    # Total TM domains (sum over chains)
    tm_total = 0
    for ch_entry in chains.values():
        tm_total += len(ch_entry.get("tm_domains") or [])
    ET.SubElement(general, "TMDomains").text = str(tm_total)

    # GeneralStructure: PDB file, image, t3dmol
    general_structure = structure.get("general_structure") or {}
    gs = ET.SubElement(general, "GeneralStructure")
    ET.SubElement(gs, "DownloadPDBFile").text = general_structure.get("pdb_file", "") or ""
    ET.SubElement(gs, "ImageFile").text = general_structure.get("image_file", "") or ""
    ET.SubElement(gs, "T3dmolFile").text = general_structure.get("t3dmol_file", "") or ""

    urls = structure.get("urls") or {}
    ET.SubElement(general, "PDB_URL").text = urls.get("pdb", "") or ""
    ET.SubElement(general, "OPM_URL").text = urls.get("opm", "") or ""
    ET.SubElement(general, "PDBTM_URL").text = urls.get("pdbtm", "") or ""

    res_elem = ET.SubElement(general, "Resolution")
    res_val = structure.get("resolution")
    res_elem.text = "" if res_val is None else str(res_val)

    # ----- StructureInformation -----
    sinfo = ET.SubElement(struct_elem, "StructureInformation")
    ET.SubElement(sinfo, "Member").text = pdb_code
    ET.SubElement(sinfo, "Class").text = prot_class
    ET.SubElement(sinfo, "TMDomains").text = str(tm_total)

    struct_neighbors = entry.get("neighbors") or {}
    for key, tag in [
        ("sequence_neighbors", "SequenceNeighbors"),
        ("structure_neighbors", "StructureNeighbors"),
        ("total_neighbors", "TotalNeighbors"),
    ]:
        val = struct_neighbors.get(key)
        el = ET.SubElement(sinfo, tag)
        el.text = "" if val is None else str(val)

    # ----- Deletion info (if any) -----
    deletion = entry.get("deletion") or {}
    if deletion.get("code") or deletion.get("message"):
        del_elem = ET.SubElement(struct_elem, "Deletion")
        if deletion.get("code"):
            ET.SubElement(del_elem, "Code").text = str(deletion["code"])
        if deletion.get("message"):
            ET.SubElement(del_elem, "Message").text = deletion["message"]

    # ----- Symmetry (whole-structure alignments only) -----
    symm = entry.get("symmetry") or {}
    symm_elem = ET.SubElement(struct_elem, "Symmetry")

    cesymm_whole = symm.get("cesymm_alignment_whole") or ""
    symd_whole = symm.get("symd_alignment_whole") or ""

    if cesymm_whole:
        insert_alignment_xml(symm_elem, "CE-SymmAlignment", cesymm_whole)
    if symd_whole:
        insert_alignment_xml(symm_elem, "SymdAlignment", symd_whole)

    # You could also add a SymmetryAnalysis section based on symm["analysis_alignments"]
    # if you decide to expose MSSD-style alignments at whole-structure level.


# ---------------------------------------------------------------------------
# Chain-level XML
# ---------------------------------------------------------------------------

def build_chain_structure_xml(
    root: ET.Element,
    pdb_code: str,
    struct_entry: Dict[str, Any],
    chain_id: str,
    chain_entry: Dict[str, Any],
) -> None:
    """
    Create the <Structure> element for a specific chain (PDB_chain)
    and attach it to the 'chains' XML root, using the new web_data schema.
    """
    struct_chain_id = f"{pdb_code}_{chain_id}"
    struct_elem = ET.SubElement(root, "Structure", ID=struct_chain_id)

    title_text = struct_entry.get("title", "") or ""
    prot_class = struct_entry.get("class", "") or ""
    structure = struct_entry.get("structure", {}) or {}

    sequence = chain_entry.get("sequence") or ""
    tm_domains = chain_entry.get("tm_domains") or []
    chain_uniprot_accs = chain_entry.get("uniprot_accs") or []

    # ----- Header -----
    header = ET.SubElement(struct_elem, "Header")
    title = ET.SubElement(header, "Title")
    title.text = title_text

    # ----- General -----
    gen = ET.SubElement(struct_elem, "General")
    ET.SubElement(gen, "PDBCode").text = pdb_code
    ET.SubElement(gen, "ChainID").text = chain_id
    ET.SubElement(gen, "Class").text = prot_class

    method = structure.get("method", "") or ""
    ET.SubElement(gen, "Method").text = method

    # Chain size
    size_elem = ET.SubElement(gen, "Size")
    size_val = chain_entry.get("size")
    if size_val is None and sequence:
        size_val = len(sequence)
    size_elem.set("value", str(size_val or 0))

    # ChainUniProtAccs (per chain)
    cu_elem = ET.SubElement(gen, "ChainUniProtAccs")
    set_chain_uniprot_attrs(cu_elem, chain_uniprot_accs)

    # TM domains (count + details)
    ET.SubElement(gen, "TMDomains").text = str(len(tm_domains))

    seq_elem = ET.SubElement(gen, "Sequence", name=f">{pdb_code}_{chain_id}:")
    seq_text_elem = ET.SubElement(seq_elem, "Seq")
    seq_text_elem.text = sequence

    for idx, (start, end) in enumerate(tm_domains):
        tm_elem = ET.SubElement(gen, "TMDomain", ID=str(idx))
        tm_range = ET.SubElement(tm_elem, "TMDRange")
        tm_range.text = f"{start} - {end}"

    # GeneralStructure (whole-structure files reused)
    general_structure = structure.get("general_structure") or {}
    gs = ET.SubElement(gen, "GeneralStructure")
    ET.SubElement(gs, "DownloadPDBFile").text = general_structure.get("pdb_file", "") or ""
    ET.SubElement(gs, "ImageFile").text = general_structure.get("image_file", "") or ""
    ET.SubElement(gs, "T3dmolFile").text = general_structure.get("t3dmol_file", "") or ""

    urls = structure.get("urls") or {}
    ET.SubElement(gen, "PDB_URL").text = urls.get("pdb", "") or ""
    ET.SubElement(gen, "OPM_URL").text = urls.get("opm", "") or ""
    ET.SubElement(gen, "PDBTM_URL").text = urls.get("pdbtm", "") or ""

    res_elem = ET.SubElement(gen, "Resolution")
    res_val = structure.get("resolution")
    res_elem.text = "" if res_val is None else str(res_val)

    # ----- StructureInformation -----
    sinfo = ET.SubElement(struct_elem, "StructureInformation")
    ET.SubElement(sinfo, "Member").text = struct_chain_id
    ET.SubElement(sinfo, "Class").text = prot_class
    ET.SubElement(sinfo, "TMDomains").text = str(len(tm_domains))

    neigh = chain_entry.get("neighbors") or {}
    ET.SubElement(
        sinfo, "SequenceNeighbors"
    ).text = "" if neigh.get("sequence_neighbors") is None else str(neigh.get("sequence_neighbors"))
    ET.SubElement(
        sinfo, "StructureNeighbors"
    ).text = "" if neigh.get("structure_neighbors") is None else str(neigh.get("structure_neighbors"))
    ET.SubElement(
        sinfo, "TotalNeighbors"
    ).text = "" if neigh.get("total_neighbors") is None else str(neigh.get("total_neighbors"))

    # Neighbor files
    neighbor_files = chain_entry.get("neighbor_files") or {}
    ET.SubElement(sinfo, "SequenceNeighborsFile").text = neighbor_files.get(
        "sequence_neighbors_file", ""
    )
    ET.SubElement(sinfo, "StructureNeighborsFile").text = neighbor_files.get(
        "structure_neighbors_file", ""
    )
    ET.SubElement(sinfo, "TotalNeighborsFile").text = neighbor_files.get(
        "total_neighbors_file", ""
    )

    # Images block
    images = chain_entry.get("images") or {}
    images_elem = ET.SubElement(sinfo, "Images")
    ET.SubElement(images_elem, "RadialDistribution").text = images.get(
        "radial_distribution", ""
    )
    ET.SubElement(images_elem, "RadialDistributionMap").text = images.get(
        "radial_distribution_map", ""
    )
    ET.SubElement(images_elem, "AlignmentsInDensityPlot").text = images.get(
        "density_plot", ""
    )
    ET.SubElement(images_elem, "AlignmentsInDensityPlotMap").text = images.get(
        "density_plot_map", ""
    )
    ET.SubElement(images_elem, "ResDistanceDistribution").text = images.get(
        "res_distance_distribution", ""
    )
    ET.SubElement(images_elem, "Topology").text = images.get("topology", "")

    # AdditionalFiles block
    additional_files = chain_entry.get("additional_files") or {}
    add_elem = ET.SubElement(sinfo, "AdditionalFiles")
    ET.SubElement(add_elem, "ResDistanceDistribution").text = additional_files.get(
        "res_distance_distribution", ""
    )
    ET.SubElement(add_elem, "StructureAlignments").text = additional_files.get(
        "structure_alignments", ""
    )
    ET.SubElement(add_elem, "Estimators").text = additional_files.get(
        "estimators", ""
    )

    # ----- Symmetry (per chain) -----
    symm_elem = ET.SubElement(struct_elem, "Symmetry")

    chain_consensus = chain_entry.get("symmetry_consensus") or {}
    ET.SubElement(symm_elem, "ReportedOrder").text = str(
        chain_consensus.get("reported_order", "na")
    )
    ET.SubElement(symm_elem, "ConsensusScore").text = str(
        chain_consensus.get("consensus_score", "na")
    )

    # Alignments (chain-level CE-Symm / SymD) with proper <Aln> children
    ce_aln = chain_entry.get("cesymm_alignment") or ""
    sd_aln = chain_entry.get("symd_alignment") or ""

    if ce_aln:
        insert_alignment_xml(symm_elem, "CE-SymmAlignment", ce_aln)
    if sd_aln:
        insert_alignment_xml(symm_elem, "SymdAlignment", sd_aln)

    # Transfer / SymmetryAnalysis sections from the original XML are not
    # reproduced in full here, because we no longer keep all those metrics
    # and file paths in web_data. If needed later, we can extend web_data
    # and then mirror those tags too.


# ---------------------------------------------------------------------------
# Top-level XML creator
# ---------------------------------------------------------------------------

def create_xml_from_web_data(
    web_data: Dict[str, Any],
    whole_xml_path: str,
    chains_xml_path: str,
) -> None:
    """
    Build 'whole' and 'chains' XML files from the aggregated web_data
    dictionary produced by create_data_struct.py.
    """
    # Whole structures XML
    whole_root = ET.Element("Database")
    for pdb_code, entry in sorted(web_data.items()):
        build_whole_structure_xml(whole_root, pdb_code, entry)

    indent(whole_root)
    whole_tree = ET.ElementTree(whole_root)
    whole_dir = os.path.dirname(whole_xml_path)
    if whole_dir:
        os.makedirs(whole_dir, exist_ok=True)
    whole_tree.write(whole_xml_path, encoding="utf-8", xml_declaration=True)

    # Chain XML
    chains_root = ET.Element("Database")
    for pdb_code, entry in sorted(web_data.items()):
        chains = entry.get("chains", {}) or {}
        for chain_id, chain_entry in sorted(chains.items()):
            build_chain_structure_xml(chains_root, pdb_code, entry, chain_id, chain_entry)

    indent(chains_root)
    chains_tree = ET.ElementTree(chains_root)
    chains_dir = os.path.dirname(chains_xml_path)
    if chains_dir:
        os.makedirs(chains_dir, exist_ok=True)
    chains_tree.write(chains_xml_path, encoding="utf-8", xml_declaration=True)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Create EncoMPASS XML files (whole/chains) from the aggregated "
            "web_data structure produced by create_data_struct.py."
        )
    )
    parser.add_argument(
        "--web-data",
        "-w",
        dest="web_data_path",
        required=True,
        help="Path to web_data pickle (output of create_data_struct.py).",
    )
    parser.add_argument(
        "--output-prefix",
        "-o",
        dest="output_prefix",
        required=True,
        help=(
            "Prefix for the output XML files. "
            "Two files will be created: "
            "<prefix>_whole.xml and <prefix>_chains.xml."
        ),
    )

    args = parser.parse_args(argv)

    print("Loading web_data from:", args.web_data_path)
    with open(args.web_data_path, "rb") as f:
        web_data = pkl.load(f)

    whole_xml_path = args.output_prefix + "_whole.xml"
    chains_xml_path = args.output_prefix + "_chains.xml"

    print("Writing whole XML to:", whole_xml_path)
    print("Writing chains XML to:", chains_xml_path)
    create_xml_from_web_data(web_data, whole_xml_path, chains_xml_path)
    print("Done.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
