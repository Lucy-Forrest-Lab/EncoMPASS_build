import analysis_parallel
import pickle as pkl
import time
from initialize_repository import *
from combine_sources import *
childdir = os.path.dirname(os.path.abspath(__file__)) + "/symmetry"
sys.path.append(childdir)
from symmetry_exec_functions import cesymm_alignment
from symmetry_exec_functions import symd_alignment


def eprint(*args, **kwargs): # enables printing to standard err
    print(*args, file=sys.stderr, **kwargs)


def format_time(timestamp):
    return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(timestamp))


def read_exclusions(delete_list):
    exclude_pdbs = set({})
    if delete_list:
        with open(delete_list, 'r') as f:
            for line in f:
                if not line.startswith("#"):
                    exclude_pdbs.add(line[0:4])
                    #print("Excluding:", line[0:4])
    return exclude_pdbs

def squash_webarchive(options, locations, generate=True):

    # Create directory where to place the squashfs archive
    out_dir = options['PATHS'][('waout', 'wa_out_dir')]
    if not os.path.exists(options['PATHS'][('waout', 'wa_out_dir')]):
        os.mkdir(options['PATHS'][('waout', 'wa_out_dir')])
    squash_fn = out_dir + '/' + options['RUN'][('code', 'main_code')] + '.squashfs' 
    if os.path.exists(squash_fn):
        os.remove(squash_fn)
    p = os.getcwd()
    if squash_fn[0] != "/":
        squash_fn = p + "/" + squash_fn

    # Check the existence of the include list
    if not os.path.exists(locations['SYSFILES']['squashin']):
        print("ERROR: please provide include list!", locations['SYSFILES']['squashin'])
        exit(1)

    # Read include list and add files and directories to be copied
    incldir_set = set()
    inclfile_set = set()
    with open(locations['SYSFILES']['squashin']) as sf:
        for line in sf:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split('\t')
            k = fields[0].strip()
            if k not in locations['FSYS']:
                print("ERROR: {0} not a valid keyword for the locations dictionary!".format(k))
                exit(1)
            if len(fields) == 2:
                filts = fields[1].strip()
                os.chdir(locations['FSYSPATH']['main'])
                outs = subprocess.run(["ls", locations['FSYS'][k[0]]] + filts.split(), stdout=subprocess.PIPE).stdout.decode('utf8').split()
                for o in outs:
                    print(o)
                    inclfile_set.add(o)
                os.chdir(p)
            elif len(fields) == 1:
                incldir_set.add(locations['FSYS'][k])
    if not incldir_set:
        print("ERROR: Empty SQUASHFS list!")
        exit(1)

    # Invokes mksquashfs from appropriate singularity container. Squashfs version must include the --no-strip flag
    # NOTE: for this to work, process must be placed in the main directory of the database to be squashed!
    print("Changing directory to", locations['FSYSPATH']['main'])
    os.chdir(locations['FSYSPATH']['main'])
    print("Running mksquashfs")
    print(" ".join(["singularity", "exec", options['PATHS'][('container', 'singularity_container')], "mksquashfs"] + sorted(list(incldir_set)) + sorted(list(inclfile_set)) + [squash_fn, "-no-strip"]))
    subprocess.run(
        ["singularity", "exec", options['PATHS'][('container', 'singularity_container')], "mksquashfs"] + sorted(
            list(incldir_set)) + sorted(list(inclfile_set)) + [squash_fn, "-no-strip"])

    print("Changing directory back to", p)
    os.chdir(p)

    # Example command for mounting squashfs
    # NOTE: your system could not allow mounting images like this!
    #print("EXAMPLE mount dir.sqsh /mnt/dir -t squashfs -o loop")

# xml attribute generator
def value_enumerator(xml_entry_name, dict_value, xml_key, tab_num=4, dict_value_prefix="", check_dir=False):
    tabs = '\t'*tab_num
    xml_txt = '{0}<{1} '.format(tabs, xml_entry_name)
    if '&amp;' in dict_value:
        xml_txt += '{0}{1}="{2}{3}" '.format(xml_key, "1", dict_value_prefix, str(dict_value))
    else:
        a: object
        for i,a in enumerate(dict_value.strip(';').split(';')):
            if check_dir:
                if not os.path.isfile(check_dir + str(a)):
                    if str(a):
                        print(dict_value.strip(';').split(';'))
                        full_path = check_dir + str(a)
                        print(f"WARNING: No image found for {full_path} in {xml_entry_name}")
                        continue
            if str(a) == "na":
                continue
            # e.g. <ImageFile name="webarchive/symm_chains_pngs/1a0s_P_symd.png" />
            xml_txt += '{0}{1}="{2}{3}" '.format(xml_key, str(i+1), dict_value_prefix, str(a))
    xml_txt += '/>\n'
    return xml_txt


def format_alignments_from_string(aln):
    if aln == "":
        return ""
    xml_txt = ""
    for a in aln.strip('\n').split('\n'):
        xml_txt += '\t' + a + "\n"
    return xml_txt


def cesymm_empty(locations, struct):
    txt = '\t\t\t<CE-Symm>\n'
    txt += '\t\t\t\t<CE-SymmOrder />\n'
    txt += '\t\t\t\t<CE-SymmAngle />\n'
    txt += '\t\t\t\t<CE-SymmTranslation />\n'
    txt += '\t\t\t\t<CE-SymmRMSD />\n'
    txt += '\t\t\t\t<CE-SymmTMScore />\n'
    txt += '\t\t\t\t<CE-SymmCoverage />\n'
    txt += '\t\t\t\t<CE-SymmAlignment>\n'
    txt += '\t\t\t\t</CE-SymmAlignment>\n'
    txt += '\t\t\t\t<AdditionalCE-Symm>\n'
    txt += '\t\t\t\t\t<CE-SymmNumRepeats />\n'
    txt += '\t\t\t\t\t<CE-SymmLevels />\n'
    txt += '\t\t\t\t\t<CE-SymmRepLength />\n'
    txt += '\t\t\t\t\t<CE-SymmAlignedResidues />\n'
    txt += '\t\t\t\t\t<CE-SymmUnrefinedRMSD />\n'
    txt += '\t\t\t\t\t<CE-SymmUnrefinedTMScore />\n'
    txt += '\t\t\t\t\t<CE-SymmSeed />\n'
    txt += '\t\t\t\t\t<CE-SymmRepeats />\n'
    txt += '\t\t\t\t</AdditionalCE-Symm>\n'
    txt += '\t\t\t\t<DownloadCE-Symm>\n'
    for suffix in ["_stdout.out", "_output.xml", ".axes", ".fasta"]:
        if os.path.isfile(locations['FSYSPATH']['cesymm'] + struct + suffix):
            txt += '\t\t\t\t\t<Path> {0}{1}{2} </Path>\n'.format(locations['FSYS']['cesymm'], struct, suffix)
    txt += '\t\t\t\t</DownloadCE-Symm>\n'
    txt += '\t\t\t\t<CE-SymmImage>\n'
    txt += '\t\t\t\t\t<ImageFile />\n'
    txt += '\t\t\t\t\t<PmlFile />\n'
    txt += '\t\t\t\t\t<T3dmolFile />\n'
    txt += '\t\t\t\t</CE-SymmImage>\n'
    txt += '\t\t\t</CE-Symm>\n'
    return txt

def symd_empty(locations, struct):
    txt = '\t\t\t<SymD>\n'
    txt += '\t\t\t\t<SymdOrder />\n'
    txt += '\t\t\t\t<SymdAngle />\n'
    txt += '\t\t\t\t<SymdTranslation />\n'
    txt += '\t\t\t\t<SymdRMSD />\n'
    txt += '\t\t\t\t<SymdTMScore />\n'
    txt += '\t\t\t\t<SymdCoverage />\n'
    txt += '\t\t\t\t<SymdAlignment>\n'
    txt += '\t\t\t\t</SymdAlignment>\n'
    txt += '\t\t\t\t<AdditionalSymD>\n'
    txt += '\t\t\t\t\t<SymdIS />\n'
    txt += '\t\t\t\t\t<SymdISAngle />\n'
    txt += '\t\t\t\t\t<SymdISTranslation />\n'
    txt += '\t\t\t\t\t<SymdZTMScore />\n'
    txt += '\t\t\t\t\t<SymdAlignedResidues />\n'
    txt += '\t\t\t\t</AdditionalSymD>\n'
    txt += '\t\t\t\t<DownloadSymD>\n'
    for suffix in ["-info.txt", "-trfm.pdb", "-best.fasta"]:
        if os.path.isfile(locations['FSYSPATH']['symd'] + struct + suffix):
            txt += '\t\t\t\t\t<Path> {0}{1}{2} </Path>\n'.format(locations['FSYS']['symd'], struct, suffix)
    txt += '\t\t\t\t</DownloadSymD>\n'
    txt += '\t\t\t\t<SymDImage>\n'
    txt += '\t\t\t\t\t<ImageFile />\n'
    txt += '\t\t\t\t\t<PmlFile />\n'
    txt += '\t\t\t\t\t<T3dmolFile />\n'
    txt += '\t\t\t\t</SymDImage>\n'
    txt += '\t\t\t</SymD>\n'
    return txt

def mssd_empty():
    # Symmetry analysis
    txt = '\t\t\t<AnalysisChainOrder />\n'
    txt += '\t\t\t<AnalysisOrder name="C1" />\n'
    txt += '\t\t\t<AnalysisAngle angle1="null" />\n'
    txt += '\t\t\t<AnalysisTranslation transl1="null" />\n'
    txt += '\t\t\t<AnalysisRMSD value="null" />\n'
    txt += '\t\t\t<AnalysisTMScore value="null" />\n'
    txt += '\t\t\t<AnalysisCoverage value="0.0" />\n'
    txt += '\t\t\t<AnalysisNumRepeats value="1" />\n'
    txt += '\t\t\t<AnalysisLevels value="0" />\n'
    txt += '\t\t\t<AnalysisRepLength value="0" />\n'
    txt += '\t\t\t<AnalysisAlignedResidues value="0" />\n'
    txt += '\t\t\t<AnalysisQuaternaryInternal type1="Quaternary" />\n'
    txt += '\t\t\t<AnalysisTopology />\n'
    txt += '\t\t\t<AnalysisAxisAngle />\n'
    txt += '\t\t\t<AnalysisRepeats />\n'
    txt += '\t\t\t<AnalysisAlignment>\n'
    txt += '\t\t\t</AnalysisAlignment>\n'
    txt += '\t\t\t<DownloadAnalysis>\n'
    txt += '\t\t\t</DownloadAnalysis>\n'
    txt += '\t\t\t<AnalysisImage>\n'
    txt += '\t\t\t\t<ImageFile />\n'
    txt += '\t\t\t\t<PmlFile />\n'
    txt += '\t\t\t\t<T3dmolFile />\n'
    txt += '\t\t\t\t<SuperPmlFile />\n'
    txt += '\t\t\t</AnalysisImage>\n'
    return txt


def transfer_empty():
    txt = '\t\t\t<TransferTemplate />\n'
    txt += '\t\t\t<TransferOrder name="C1" />\n'
    txt += '\t\t\t<TransferAngle angle1="null" />\n'
    txt += '\t\t\t<TransferTranslation transl1="null" />\n'
    txt += '\t\t\t<TransferRMSD value="null" />\n'
    txt += '\t\t\t<TransferTMScore value="null" />\n'
    txt += '\t\t\t<TransferCoverage value="0.0" />\n'
    txt += '\t\t\t<TransferNumRepeats value="1" />\n'
    txt += '\t\t\t<TransferLevels value="0" />\n'
    txt += '\t\t\t<TransferRepLength value="0" />\n'
    txt += '\t\t\t<TransferAlignedResidues value="0" />\n'
    txt += '\t\t\t<TransferTopology />\n'
    txt += '\t\t\t<TransferAxisAngle />\n'
    txt += '\t\t\t<TransferRepeats />\n'
    txt += '\t\t\t<DownloadTransfer>\n'
    txt += '\t\t\t</DownloadTransfer>\n'
    txt += '\t\t\t<TransferImage>\n'
    txt += '\t\t\t\t<ImageFile />\n'
    txt += '\t\t\t\t<PmlFile />\n'
    txt += '\t\t\t\t<T3dmolFile />\n'
    txt += '\t\t\t\t<SuperPmlFile />\n'
    txt += '\t\t\t</TransferImage>\n'
    return txt


def cesymm_xml(locations, struct, symm_dic_pdbi, xml_type, monomeric=False):
    # Symmetry
    txt = '\t\t\t<CE-Symm Algorithm="jCE-symm" Version="2.2" >\n'
    txt += '\t\t\t\t<CE-SymmOrder name="{0}" />\n'.format(symm_dic_pdbi['cesymm']['symmetry_order'])
    txt += value_enumerator("CE-SymmAngle", symm_dic_pdbi['cesymm']['unit_angle'], "angle", tab_num=4)
    txt += value_enumerator("CE-SymmTranslation", symm_dic_pdbi['cesymm']['unit_translation'], "transl", tab_num=4)
    txt += '\t\t\t\t<CE-SymmRMSD value="{0}" />\n'.format(symm_dic_pdbi['cesymm']['refined_rmsd'])
    txt += '\t\t\t\t<CE-SymmTMScore value="{0}" />\n'.format(symm_dic_pdbi['cesymm']['refined_tmscore'])
    txt += '\t\t\t\t<CE-SymmCoverage value="{0}" />\n'.format(symm_dic_pdbi['cesymm']['coverage'])
    txt += '\t\t\t\t<CE-SymmAlignment>\n'
    print("\tPrinting alignment for", struct, xml_type)
    if 'alignment' in symm_dic_pdbi['cesymm'] and symm_dic_pdbi['cesymm']['alignment']:
        print("\t\tusing keyword alignment")
        aln = format_alignments_from_string(symm_dic_pdbi['cesymm']['alignment'])
    elif 'alignmnet' in symm_dic_pdbi['cesymm'] and symm_dic_pdbi['cesymm']['alignmnet']:
        print("\t\tusing keyword alignmnet")
        aln = format_alignments_from_string(symm_dic_pdbi['cesymm']['alignmnet']) 
    else:
        if xml_type == 'whole':
            fasta_name = struct[0:4]
        elif xml_type == 'chains':
            fasta_name = struct
        fasta_root = locations['FSYSPATH']['main'] + locations['FSYS']['cesymm'] + struct[0:4] + "/" + fasta_name
        print("\t\tNo entry for CE-Symm alignment in dictionary, trying fasta in", fasta_root, ".fasta")
        if os.path.isfile(fasta_root + '.fasta'):
            aln = cesymm_alignment(locations['FSYSPATH']['main'], locations['FSYS']['cesymm'] + struct[0:4] + '/' + fasta_name)
            if not aln:
                print("\t\tAlignment file exists but is empty of alignment: ", fasta_root, ".fasta")
        else:
            print("\tWARNING, expected alignment file doesn't exist: ", fasta_root, ".fasta")
            aln = ""

    txt += aln
    txt += '\t\t\t\t</CE-SymmAlignment>\n'
    txt += '\t\t\t\t<AdditionalCE-Symm>\n'
    txt += '\t\t\t\t\t<CE-SymmNumRepeats value="{0}" />\n'.format(symm_dic_pdbi['cesymm']['repeats_number'])
    txt += '\t\t\t\t\t<CE-SymmLevels value="{0}" />\n'.format(symm_dic_pdbi['cesymm']['symmetry_levels'])
    txt += '\t\t\t\t\t<CE-SymmRepLength value="{0}" />\n'.format(symm_dic_pdbi['cesymm']['repeat_length'])
    txt += '\t\t\t\t\t<CE-SymmAlignedResidues value="{0}" />\n'.format(symm_dic_pdbi['cesymm']['aligned_length'])
    txt += '\t\t\t\t\t<CE-SymmUnrefinedRMSD value="{0}" />\n'.format(symm_dic_pdbi['cesymm']['unrefined_rmsd'])
    txt += '\t\t\t\t\t<CE-SymmUnrefinedTMScore value="{0}" />\n'.format(symm_dic_pdbi['cesymm']['unrefined_tmscore'])
    txt += '\t\t\t\t\t<CE-SymmSeed value="{0}" />\n'.format(symm_dic_pdbi['cesymm']['seed'])
    txt += value_enumerator("CE-SymmRepeats", symm_dic_pdbi['cesymm']['repeats'], "name", tab_num=5)
    txt += '\t\t\t\t</AdditionalCE-Symm>\n'
    txt += '\t\t\t\t<DownloadCE-Symm>\n'

    for suffix in ["_stdout.out", "_output.xml", ".axes", ".fasta"]:
        path = locations['FSYSPATH']['cesymm'] + struct[0:4] + "/" + struct + suffix
        if os.path.isfile(path):
            txt += '\t\t\t\t\t<Path> {0} </Path>\n'.format(locations['FSYS']['cesymm'] + struct[0:4] + "/" + struct + suffix)
    txt += '\t\t\t\t</DownloadCE-Symm>\n'
    txt += '\t\t\t\t<CE-SymmImage>\n'
    if monomeric:
        xml_type = "whole"
    webdir = "symm_" + xml_type
    if 'image_files' not in symm_dic_pdbi['cesymm'] or not symm_dic_pdbi['cesymm']['image_files']:
        print("\tWARNING: ", struct, "does not have image_files, pml_files and jmol_files")
        txt += '\t\t\t\t\t<ImageFile />\n'
        txt += '\t\t\t\t\t<PmlFile />\n'
        txt += '\t\t\t\t\t<T3dmolFile />\n'
    else:
        print("\t\tWriting paths for CE-Symm image_files, pml_files, jmol_files", webdir, locations['FSYS'][webdir + 'pngs'], symm_dic_pdbi['cesymm']['image_files'])
        txt += value_enumerator("ImageFile", symm_dic_pdbi['cesymm']['image_files'], "name", tab_num=5, dict_value_prefix=locations['FSYS'][webdir + 'pngs'], check_dir=locations['FSYSPATH'][webdir + 'pngs'])
        txt += value_enumerator("PmlFile", symm_dic_pdbi['cesymm']['pml_files'], "name", tab_num=5, dict_value_prefix=locations['FSYS'][webdir + 'pngs'], check_dir=locations['FSYSPATH'][webdir + 'pngs'])
        txt += value_enumerator("T3dmolFile", symm_dic_pdbi['cesymm']['jmol_files'], "name", tab_num=5, dict_value_prefix=locations['FSYS'][webdir + 'jsons'], check_dir=locations['FSYSPATH'][webdir + 'jsons'])
    txt += '\t\t\t\t</CE-SymmImage>\n'
    txt += '\t\t\t</CE-Symm>\n'
    return txt


def symd_xml(locations, struct, symm_dic_pdbi, xml_type, monomeric=False):
    txt = '\t\t\t<SymD Algorithm="SymD" Version="1.6" >\n'
    if symm_dic_pdbi['symd']['symmetry_order'] == 'na':
        txt += '\t\t\t\t<SymdOrder value="1" />\n'
        print("WARNING: Correcting na in SymdOrder for", struct, xml_type)
    else:
        txt += '\t\t\t\t<SymdOrder value="{0}" />\n'.format(symm_dic_pdbi['symd']['symmetry_order'])
    txt += '\t\t\t\t<SymdAngle value="{0}" />\n'.format(symm_dic_pdbi['symd']['unit_angle'])
    txt += '\t\t\t\t<SymdTranslation value="{0}" />\n'.format(symm_dic_pdbi['symd']['unit_translation'])
    txt += '\t\t\t\t<SymdRMSD value="{0}" />\n'.format(symm_dic_pdbi['symd']['rmsd'])
    txt += '\t\t\t\t<SymdTMScore value="{0}" />\n'.format(symm_dic_pdbi['symd']['tmscore'])
    txt += '\t\t\t\t<SymdCoverage value="{0}" />\n'.format(symm_dic_pdbi['symd']['coverage'])
    txt += '\t\t\t\t<SymdAlignment>\n'
    print("\tPrinting alignment for", struct, xml_type)
    if 'alignment' in symm_dic_pdbi['cesymm'] and symm_dic_pdbi['cesymm']['alignment']:
        print("\t\tusing keyword alignment")
        aln = format_alignments_from_string(symm_dic_pdbi['symd']['alignment'])
    elif 'alignmnet' in symm_dic_pdbi['cesymm'] and symm_dic_pdbi['cesymm']['alignmnet']:
        print("\t\tusing keyword alignmnet")
        aln = format_alignments_from_string(symm_dic_pdbi['symd']['alignmnet'])  # fix when dic key corrected!
    else:
        if xml_type == 'whole':
            fasta_name = struct[0:4]
        elif xml_type == 'chains':
            fasta_name = struct
        fasta_root = locations['FSYSPATH']['main'] + locations['FSYS']['symd'] + struct[0:4] + "/" + fasta_name + '-best'
        print("\t\tBut no entry for symd alignment in dictionary, so trying fasta: ", fasta_root, ".fasta")
        if os.path.isfile(fasta_root + '.fasta'):
            aln = symd_alignment(locations['FSYSPATH']['main'], locations['FSYS']['symd'] + struct[0:4] + '/' + fasta_name)
            if not aln:
                print("\t\tAlignment file exists but is empty of alignment: ", fasta_root, ".fasta")
        else:
            print("\tWARNING: expecting alignment file, but doesn't exist: ", fasta_root, ".fasta")
            aln = ""

    txt += aln
    txt += '\t\t\t\t</SymdAlignment>\n'
    txt += '\t\t\t\t<AdditionalSymD>\n'
    txt += '\t\t\t\t\t<SymdIS value="{0}" />\n'.format(symm_dic_pdbi['symd']['initial_shift'])
    txt += '\t\t\t\t\t<SymdISAngle value="{0}" />\n'.format(symm_dic_pdbi['symd']['is_angle'])
    txt += '\t\t\t\t\t<SymdISTranslation value="{0}" />\n'.format(symm_dic_pdbi['symd']['is_translation'])
    if not isinstance(symm_dic_pdbi['symd']['z_tmscore'], float):
        symm_dic_pdbi['symd']['z_tmscore'] = 0
    txt += '\t\t\t\t\t<SymdZTMScore value="{0}" />\n'.format(symm_dic_pdbi['symd']['z_tmscore'])
    txt += '\t\t\t\t\t<SymdAlignedResidues value="{0}" />\n'.format(symm_dic_pdbi['symd']['aligned_length'])
    txt += '\t\t\t\t</AdditionalSymD>\n'
    txt += '\t\t\t\t<DownloadSymD>\n'
    for suffix in ["-info.txt", "-trfm.pdb", "-best.fasta"]:
        path = locations['FSYSPATH']['symd'] + struct[0:4] + "/" + struct + suffix
        if os.path.isfile(path):
            txt += '\t\t\t\t\t<Path> {0} </Path>\n'.format(locations['FSYS']['symd'] + struct[0:4] + "/" + struct + suffix)
    txt += '\t\t\t\t</DownloadSymD>\n'
    txt += '\t\t\t\t<SymDImage>\n'
    if 'image_files' not in symm_dic_pdbi['symd']:
        symm_dic_pdbi['symd']['image_files'], symm_dic_pdbi['symd']['pml_files'], symm_dic_pdbi['symd']['jmol_files'] = "", "", ""
    txt += '\t\t\t\t\t<ImageFile name="{0}" />\n'.format(symm_dic_pdbi['symd']['image_files'])
    txt += '\t\t\t\t\t<PmlFile name="{0}" />\n'.format(symm_dic_pdbi['symd']['pml_files'])
    txt += '\t\t\t\t\t<T3dmolFile name="{0}" />\n'.format(symm_dic_pdbi['symd']['jmol_symd']) # change this dictionary keyword for consistency to "jmol_files"
    txt += '\t\t\t\t</SymDImage>\n'
    txt += '\t\t\t</SymD>\n'
    return txt


def mssd_xml(locations, symm_dic_pdbi, xml_type):
    # Symmetry analysis
    sel = symm_dic_pdbi['selected']['source'].strip(';').split(';')
    txt = '\t\t\t<AnalysisChainOrder value="{0}" />\n'.format(symm_dic_pdbi[sel[0]]['chains'])
    txt += '\t\t\t<AnalysisOrder name="{0}" />\n'.format(' &amp; '.join([symm_dic_pdbi[s]['symmetry_order'] for s in sel]))
    txt += value_enumerator("AnalysisAngle", ' &amp; '.join([str(symm_dic_pdbi[s]['unit_angle']).strip(';') for s in sel]), "angle", tab_num=3)
    txt += value_enumerator("AnalysisTranslation", ' &amp; '.join([str(symm_dic_pdbi[s]['unit_translation']).strip(';') for s in sel]), "transl", tab_num=3)
    txt += '\t\t\t<AnalysisRMSD value="{0}" />\n'.format(' &amp; '.join([str(symm_dic_pdbi[s]['refined_rmsd']).strip(';') for s in sel]))
    txt += '\t\t\t<AnalysisTMScore value="{0}" />\n'.format(' &amp; '.join([str(symm_dic_pdbi[s]['refined_tmscore']) for s in sel]))
    txt += '\t\t\t<AnalysisCoverage value="{0}" />\n'.format(' &amp; '.join([str(symm_dic_pdbi[s]['coverage']) for s in sel]))
    txt += '\t\t\t<AnalysisNumRepeats value="{0}" />\n'.format(' &amp; '.join([str(symm_dic_pdbi[s]['repeats_number']).strip(';') for s in sel]))
    txt += '\t\t\t<AnalysisLevels value="{0}" />\n'.format(' &amp; '.join([str(symm_dic_pdbi[s]['symmetry_levels']).strip(';') for s in sel]))
    txt += '\t\t\t<AnalysisRepLength value="{0}" />\n'.format(' &amp; '.join([str(symm_dic_pdbi[s]['repeat_length']).strip(';') for s in sel]))
    txt += '\t\t\t<AnalysisAlignedResidues value="{0}" />\n'.format(' &amp; '.join([str(symm_dic_pdbi[s]['aligned_length']).strip(';') for s in sel]))
    txt += value_enumerator("AnalysisQuaternaryInternal", ' &amp; '.join([symm_dic_pdbi[s]['symmetry_type'].strip(';') for s in sel]), "type", tab_num=3)
    txt += value_enumerator("AnalysisTopology", ' &amp; '.join([symm_dic_pdbi[s]['topology'].strip(';') for s in sel]), "type", tab_num=3)
    txt += value_enumerator("AnalysisAxisAngle", ' &amp; '.join([symm_dic_pdbi[s]['axis_angle_with_membrane_normal'].strip(';') for s in sel]), "angle", tab_num=3)
    txt += value_enumerator("AnalysisRepeats", ' &amp; '.join([symm_dic_pdbi[s]['repeats'].strip(';') for s in sel]), "name", tab_num=3)
    txt += '\t\t\t<AnalysisAlignment>\n'
    for s in sel:
        if 'alignment' in symm_dic_pdbi[s]:
            aln = symm_dic_pdbi[s]['alignment']
        elif 'alignmnet' in symm_dic_pdbi[s]:
            aln = symm_dic_pdbi[s]['alignmnet']
        else:
            aln = cesymm_alignment(locations['FSYSPATH']['main'] + symm_dic_pdbi[s]['files_dir'], symm_dic_pdbi[s]['files_key'])
        txt += aln
    txt += '\t\t\t</AnalysisAlignment>\n'
    txt += '\t\t\t<DownloadAnalysis>\n'
    for s in sel:
        for suffix in ["_stdout.out", "_output.xml", ".axes", ".fasta"]:
            # we are making this more complicated than it needs to be because quatsymm dictionary has a wrong files_dir value
            db_name = locations['FSYSPATH']['main'].strip('/').split('/')[-1]
            print(f"db_name: {db_name}, symm_dic_pdbi, {symm_dic_pdbi[s]['files_dir_path'].split(db_name + '/')} ")
            rel_path = symm_dic_pdbi[s]['files_dir_path'].split(db_name + "/")[1]
            path = locations['FSYSPATH']['main'] + rel_path + symm_dic_pdbi[s]['files_key'] + suffix
            print(path)
            if os.path.isfile(path):
                txt += '\t\t\t\t<Path> {0} </Path>\n'.format(rel_path + symm_dic_pdbi[s]['files_key'] + suffix)
    txt += '\t\t\t</DownloadAnalysis>\n'
    txt += '\t\t\t<AnalysisImage>\n'
    image_files, pml_files, t3dmol_files, superpml_files = "", "", "", ""
    #print(f"symmdic: {symm_dic_pdbi}")
    print(f"sel: {sel}")
    for s in sel:
        if s == 'cesymm':
            webdir = 'symm_' + xml_type
        else:
            webdir = 'analysis_' + xml_type
        print(f"s: {s}, webdir: {webdir}, xml: {xml_type}")
        image_files += locations['FSYS'][webdir + 'pngs'] + symm_dic_pdbi[s]['image_files'].strip(';') + ";"
        print(f"image_files: {image_files}")
        pml_files += locations['FSYS'][webdir + 'pngs'] + symm_dic_pdbi[s]['pml_files'].strip(';') + ";"
        print(f"pml_files: {pml_files}")
        t3dmol_files += locations['FSYS'][webdir + 'jsons'] + symm_dic_pdbi[s]['jmol_files'].strip(';') + ";"
        print(f"t3dmol files: {t3dmol_files}")
        print(locations['FSYS']['analysis_' + xml_type + 'super'])
        superpml_files += locations['FSYS']['analysis_' + xml_type + 'super'] + symm_dic_pdbi[s]['super_pml_files'].strip(';') + ";"
        print(f"superpml files: {superpml_files}")

    txt += value_enumerator("ImageFile", image_files.strip(';'), "name", tab_num=4, dict_value_prefix='')
    txt += value_enumerator("PmlFile", pml_files.strip(';'), "name", tab_num=4)
    txt += value_enumerator("T3dmolFile", t3dmol_files.strip(';'), "name", tab_num=4)
    txt += value_enumerator("SuperPmlFile", superpml_files.strip(';'), "name", tab_num=4)
    txt += '\t\t\t</AnalysisImage>\n'

    return txt


def transfer_xml(inferred_dic_list, locations, xml_type):
    txt = '\t\t\t<TransferTemplate ID="{0}" />\n'.format(inferred_dic_list[0]['template_structure'])
    txt += '\t\t\t<TransferOrder name="{0}" />\n'.format(inferred_dic_list[0]['template_order'])
    txt += value_enumerator("TransferAngle", ' &amp; '.join([inferred_dic_list[s]['unit_angle'].strip(';') for s in range(len(inferred_dic_list))]), "angle", tab_num=3)
    txt += value_enumerator("TransferTranslation", ' &amp; '.join([inferred_dic_list[s]['unit_translation'].strip(';') for s in range(len(inferred_dic_list))]), "transl", tab_num=3)
    txt += '\t\t\t<TransferRMSD value="{0}" />\n'.format(' &amp; '.join([inferred_dic_list[s]['refined_rmsd'].strip(';') for s in range(len(inferred_dic_list))]))
    txt += '\t\t\t<TransferTMScore value="{0}" />\n'.format(' &amp; '.join([str(inferred_dic_list[s]['refined_tmscore']) for s in range(len(inferred_dic_list))]))
    txt += '\t\t\t<TransferCoverage value="{0}" />\n'.format(' &amp; '.join([str(inferred_dic_list[s]['coverage']) for s in range(len(inferred_dic_list))]))
    txt += '\t\t\t<TransferNumRepeats value="{0}" />\n'.format(' &amp; '.join([str(inferred_dic_list[s]['repeats_number']).strip(';') for s in range(len(inferred_dic_list))]))
    txt += '\t\t\t<TransferLevels value="{0}" />\n'.format(' &amp; '.join([inferred_dic_list[s]['symmetry_levels'].strip(';') for s in range(len(inferred_dic_list))]))
    txt += '\t\t\t<TransferRepLength value="{0}" />\n'.format(' &amp; '.join([str(inferred_dic_list[s]['repeat_length']).strip(';') for s in range(len(inferred_dic_list))]))
    txt += '\t\t\t<TransferAlignedResidues value="{0}" />\n'.format(' &amp; '.join([str(inferred_dic_list[s]['aligned_length']).strip(';') for s in range(len(inferred_dic_list))]))
    txt += value_enumerator("TransferTopology", ' &amp; '.join([inferred_dic_list[s]['topology'].strip(';') for s in range(len(inferred_dic_list))]), "type", tab_num=3)
    txt += value_enumerator("TransferAxisAngle", ' &amp; '.join([inferred_dic_list[s]['axis_angle_with_membrane_normal'].strip(';') for s in range(len(inferred_dic_list))]), "angle", tab_num=3)
    txt += value_enumerator("TransferRepeats", ' &amp; '.join([inferred_dic_list[s]['repeats'].strip(';') for s in range(len(inferred_dic_list))]), "name", tab_num=3)
    txt += '\t\t\t<DownloadTransfer>\n'
    for s in range(len(inferred_dic_list)):
        for path in [inferred_dic_list[s]['axes_file'], inferred_dic_list[s]['alnres_file']]:
            if os.path.isfile(locations['FSYSPATH']['main'] + path):
                txt += '\t\t\t\t<Path> {0} </Path>\n'.format(path)
            else:
                print(f"{path} doesn't exist.")
    txt += '\t\t\t</DownloadTransfer>\n'
    txt += '\t\t\t<TransferImage>\n'
    webdir = "transfer_" + xml_type
    txt += value_enumerator("ImageFile", ';'.join([inferred_dic_list[s]['image_files'].strip(';') for s in range(len(inferred_dic_list))]), "name", tab_num=4, dict_value_prefix=locations['FSYS'][webdir + 'pngs'], check_dir=locations['FSYSPATH'][webdir + 'pngs'])
    txt += value_enumerator("PmlFile", ';'.join([inferred_dic_list[s]['pml_files'].strip(';') for s in range(len(inferred_dic_list))]), "name", tab_num=4, dict_value_prefix=locations['FSYS'][webdir + 'pngs'], check_dir=locations['FSYSPATH'][webdir + 'pngs'])
    txt += value_enumerator("T3dmolFile", ';'.join([inferred_dic_list[s]['jmol_files'].strip(';') for s in range(len(inferred_dic_list))]), "name", tab_num=4, dict_value_prefix=locations['FSYS'][webdir + 'jsons'], check_dir=locations['FSYSPATH'][webdir + 'jsons'])
    txt += value_enumerator("SuperPmlFile", ';'.join([inferred_dic_list[s]['super_pml_files'].strip(';') for s in range(len(inferred_dic_list))]), "name", tab_num=4, dict_value_prefix=locations['FSYS'][webdir + 'super'], check_dir=locations['FSYSPATH'][webdir + 'super'])
    txt += '\t\t\t</TransferImage>\n'
    return txt


def symmetry_xml(locations, struct, prot_class, xml_type ='whole', symm_dic_pdbi={}, transfer_dic=None, monomeric=False):
    txt = '\t\t<Symmetry>\n'
    if 'cesymm' in symm_dic_pdbi and symm_dic_pdbi['cesymm']['coverage'] != "na":
        print("\tFilling ce-symm")
        txt += cesymm_xml(locations, struct, symm_dic_pdbi, xml_type, monomeric)
    else:
        print("\tNot filling ce-symm")
        txt += cesymm_empty(locations, struct)
    if 'symd' in symm_dic_pdbi and symm_dic_pdbi['symd']['coverage'] != "na":
        print("\tFilling symd")
        txt += symd_xml(locations, struct, symm_dic_pdbi, xml_type, monomeric)
    else:
        print("\tNot filling symd")
        txt += symd_empty(locations, struct)
    txt += '\t\t</Symmetry>\n'
    txt += '\t\t<SymmetryAnalysis>\n'
    print("\t\tAt MSSD analysis section")

    if prot_class == "beta":
        txt += '\t\t\t<Info text="No symmetry analysis available for beta-barrel transmembrane proteins." />\n'
    elif 'selected' in symm_dic_pdbi:
        print("\tFound selected tag in symm_dic_pdbi", xml_type)
        print("\tSelected features:", symm_dic_pdbi['selected'])
        if monomeric:
            if xml_type == "whole" and 'quatinfo' in symm_dic_pdbi['selected']:
                print("\tPrint quaternary info message for monomer in xml for xml_type", xml_type)
                txt += '\t\t\t<Info text="{0}" />\n'.format(symm_dic_pdbi['selected']['quatinfo'])
                txt += mssd_empty()
            elif xml_type == "chains" and 'info' in symm_dic_pdbi['selected']:
                print("\tPrint chain info message for monomer for xml_type", xml_type)
                txt += '\t\t\t<Info text="{0}" />\n'.format(symm_dic_pdbi['selected']['info'])
                txt += mssd_empty()
            else:
                print("\tMonomer mssd, doing mssd as whole in xml for", xml_type)
                txt += mssd_xml(locations, symm_dic_pdbi, xml_type="whole")
        elif not monomeric and xml_type == "chains":
            if 'info' in symm_dic_pdbi['selected']:
                print("\tPrint chain info message for chain in complex for xml_type", xml_type)
                txt += '\t\t\t<Info text="{0}" />\n'.format(symm_dic_pdbi['selected']['info'])
                txt += mssd_empty()
            else:
                print("\tChain mssd result", xml_type)
                txt += mssd_xml(locations, symm_dic_pdbi, xml_type=xml_type)
        elif xml_type == "whole":
            if 'quatinfo' in symm_dic_pdbi['selected']:
                print("\tPrint quaternary info message for whole struct in xml for xml_type", xml_type)
                txt += '\t\t\t<Info text="{0}" />\n'.format(symm_dic_pdbi['selected']['quatinfo'])
                txt += mssd_empty()
            else:
                print("\tWhole struct mssd result", xml_type)
                txt += mssd_xml(locations, symm_dic_pdbi, xml_type=xml_type)

    elif 'selected' not in symm_dic_pdbi:
        print("\tNo selected tag in symm_dic_pdbi")
        #print(symm_dic_pdbi)
        if xml_type == "chain" or monomeric:
            print("\t\tuse unprogrammed info text because xml_type is chain:", xml_type, "or monomer:", monomeric)
            txt += '\t\t\t<Info text="{0}" />\n'.format("No symmetry found in the membrane-bound region during analysis.")
        else:
            print("\t\tfor TM complex")
            txt += '\t\t\t<Info text="{0}" />\n'.format("No quaternary symmetry found in the membrane-bound region during analysis. For internal symmetry, please check the specific chain page.")
        txt += mssd_empty()

    txt += '\t\t</SymmetryAnalysis>\n'
    if xml_type == 'chains':
        txt += '\t\t<Transfer>\n'
        if transfer_dic:
            txt += transfer_xml(transfer_dic, locations, xml_type)
        else:
            txt += '\t\t\t<Info text="No data for inferred symmetry available." />\n'
            txt += transfer_empty()
        txt += '\t\t</Transfer>\n'

    return txt


def consensus_symmetry(symm_dic, prot_class, transfer_list=None):
    cs, mssd, transfer, reported = "na", "na", "na", "na"
    if 'cesymm' in symm_dic:
        cs = symm_dic['cesymm']['symmetry_order']
        reported = cs
    if transfer_list:
        transfer = transfer_list[0]['template_order']
        reported = transfer
    if 'selected' in symm_dic:
        sel = symm_dic['selected']['source'].split(';')
        mssd = symm_dic[sel[0]]['symmetry_order']
        reported = mssd

    consensus = [cs, mssd, transfer].count(cs)
    if [cs, mssd, transfer].count(mssd) > consensus:
        consensus = [cs, mssd, transfer].count(mssd)
    if prot_class == 'beta':
        consensus = '1.0'
    elif transfer == 'na':
        consensus = [cs, mssd].count(cs)
        consensus = str(float('%.1f' % (consensus/2)))
    else:
        consensus = [cs, mssd, transfer].count(cs)
        if [cs, mssd, transfer].count(mssd) > consensus:
            consensus = [cs, mssd, transfer].count(mssd)
        consensus = str(float('%.1f' % (consensus/3)))
    return reported, consensus


def struct_xml_entry(pdbi, str_data_pdbi, symm_dic_pdbi, locations, neighbors, deletion_info, consensus, monomeric=False):

    # Neighbors
    seqneighs, strneighs, totneighs = neighbors

    # Template
    # Opening
    txt =  '\t<Structure ID="{0}">\n'.format(pdbi)
   
    deletion_codes, deletion_messages = deletion_info

    deltxt = ''
    del_key = str_data_pdbi['delete_keyword']
    if del_key:
        print("\tThis delete keyword is present:", del_key, "so only printing Info statement for", pdbi)
        deltxt = "The requested structure is not in the EncoMPASS database " + deletion_messages[del_key].replace('"', '').strip()
    for dci, dc in enumerate(deletion_codes):
        if deletion_codes[dci][0] in [x['location'] for x in str_data_pdbi['PASSPORT']]:
            deltxt = "The requested structure is not in the EncoMPASS database " + deletion_messages[del_key].replace('"', '').strip()
    if deltxt:
        txt += '\t\t<Info> text="{0}" </Info>\n'.format(deltxt)
        txt += '\t</Structure>\n'
        return txt, False

    # Header
    txt += '\t\t<Header>\n\t\t\t<Title> {0} </Title>\n\t\t</Header>\n'.format(str_data_pdbi['ENCOMPASS']['name'])

    # General
    txt += '\t\t<General>\n'
    txt += '\t\t\t<PDBCode> {0} </PDBCode>\n'.format(pdbi)
    txt += '\t\t\t<Method> {0} </Method>\n'.format(str_data_pdbi['FROM_PDB']['experimental_method'])
    if str_data_pdbi['ENCOMPASS']['resolution'] == "None":
        txt += '\t\t\t<Resolution> null </Resolution>\n'
    else:
        txt += '\t\t\t<Resolution> {0} </Resolution>\n'.format(str_data_pdbi['ENCOMPASS']['resolution'])

    txt += '\t\t\t<Size value="{0}" />\n'.format(str_data_pdbi['ENCOMPASS']['n_of_aas'])
    txt += '\t\t\t<Chains> {0} </Chains>\n'.format(str_data_pdbi['ENCOMPASS']['n_of_chains'])
    txt += '\t\t\t<TMChains> {0} </TMChains>\n'.format(str_data_pdbi['ENCOMPASS']['n_of_tmchains'])
    # account for empty first chains - LF
    codes = str_data_pdbi['ENCOMPASS']['structure']['UniProt_stoichiometry']
    non_empty_codes = [code for code in codes if code.strip()]  # Exclude empty codes
    if non_empty_codes:
        txt += '\t\t\t<UniProtAccs {0} />\n'.format(
            ' '.join(['code' + str(i + 1) + '=' + '"' + x + '"' for i, x in enumerate(non_empty_codes)]))
    else:
        txt += '\t\t\t<UniProtAccs />\n'

    ## Summary chain info in whole struture xml
    for ich, tmch in enumerate(str_data_pdbi['ENCOMPASS']['structure']['ktmchains']):
        if tmch == '-':
            continue
        txt += '\t\t\t<Chain ID="{0}">\n'.format(tmch)
        txt += '\t\t\t\t<NTMDomains> {0} </NTMDomains>\n'.format(len(str_data_pdbi['ENCOMPASS']['structure']['chains'][ich]['TM_regions']['TM_regions_extrema']))
        txt += '\t\t\t\t<Sequence name=">{0}_{1}:">\n'.format(pdbi, tmch)
        txt += '\t\t\t\t\t<Seq> {0} </Seq>\n'.format(str_data_pdbi['ENCOMPASS']['structure']['chains'][ich]['sequence'])
        txt += '\t\t\t\t</Sequence>\n'
        
        ### TM Domains
        for iextr, extr in enumerate(str_data_pdbi['ENCOMPASS']['structure']['chains'][ich]['TM_regions']['TM_regions_extrema']):
            txt += '\t\t\t\t<TMDomain ID="{0}">\n'.format(iextr)
            txt += '\t\t\t\t\t<TMDRange> {0} - {1} </TMDRange>\n'.format(extr[0][0][0], extr[1][0][0])
            txt += '\t\t\t\t</TMDomain>\n'

        txt += '\t\t\t</Chain>\n'

    txt += '\t\t\t<GeneralStructure>\n'
    txt += '\t\t\t\t<DownloadPDBFile> {0} </DownloadPDBFile>\n'.format(locations['FSYS']['whole'] + pdbi + '_enc.pdb')
    txt += '\t\t\t\t<ImageFile> {0} </ImageFile>\n'.format(locations['FSYS']['images_figs_whole'] + pdbi + '.png')
    txt += '\t\t\t\t<T3dmolFile> {0} </T3dmolFile>\n'.format(locations['FSYS']['3dmol_whole'] + pdbi + '.json')
    txt += '\t\t\t</GeneralStructure>\n'
    txt += '\t\t\t<PDB_URL> http://www.rcsb.org/pdb/explore/explore.do?structureId={0} </PDB_URL>\n'.format(pdbi.upper())
    txt += '\t\t\t<OPM_URL> http://opm.phar.umich.edu/protein.php?pdbid={0} </OPM_URL>\n'.format(pdbi)
    txt += '\t\t\t<PDBTM_URL> http://pdbtm.enzim.hu/data/database/{0}/{1}.trpdb.gz </PDBTM_URL>\n'.format(pdbi[1:3], pdbi)
    txt += '\t\t</General>\n'

    # StructureInformation
    txt += '\t\t<StructureInformation>\n'
    for ich, tmch in enumerate(str_data_pdbi['ENCOMPASS']['structure']['ktmchains']):
        if tmch == '-':
            continue
        pdbi_ch = pdbi + '_' + tmch
        #print(pdbi_ch)
        txt += '\t\t\t<ChainInformation ID="{0}_{1}">\n'.format(pdbi, tmch)
        txt += '\t\t\t\t<Member> {0}_{1} </Member>\n'.format(pdbi, tmch)
        txt += '\t\t\t\t<Class> {0} </Class>\n'.format(str_data_pdbi['ENCOMPASS']['class'])
        txt += '\t\t\t\t<TMdomains> {0} </TMdomains>\n'.format(len(str_data_pdbi['ENCOMPASS']['structure']['chains'][ich]['TM_regions']['TM_regions_extrema']))
        txt += '\t\t\t\t<SequenceNeighbors> {0} </SequenceNeighbors>\n'.format(seqneighs[pdbi_ch])
        txt += '\t\t\t\t<StructureNeighbors> {0} </StructureNeighbors>\n'.format(strneighs[pdbi_ch])
        txt += '\t\t\t\t<TotalNeighbors> {0} </TotalNeighbors>\n'.format(totneighs[pdbi_ch])
        ### TONI
        txt += '\t\t\t\t<Symmetry> {0} </Symmetry>\n'.format(consensus[pdbi_ch]['reported_order'])
        txt += '\t\t\t\t<SymmConsensus> {0} </SymmConsensus>\n'.format(consensus[pdbi_ch]['consensus_score'])
        txt += '\t\t\t</ChainInformation>\n'
    txt += '\t\t</StructureInformation>\n'

    # Symmetry
    txt += symmetry_xml(locations, pdbi, str_data_pdbi['ENCOMPASS']['class'], xml_type ='whole', symm_dic_pdbi=symm_dic_pdbi, monomeric=monomeric)

    txt += '\t</Structure>\n'

    return txt, True


def chain_xml_entry(pdbi, ch, str_data_pdbi, symm_dic_pdbi, locations, neighbor_counts, transfer_list=None, monomeric=False):
    seqneighs, strneighs, totneighs = neighbor_counts
    ich = str_data_pdbi['ENCOMPASS']['structure']['kchains'].index(ch)
    pdbi_ch = pdbi + "_" + ch

    txt =  '\t<Structure ID="{0}_{1}">\n'.format(pdbi, ch)
    txt += '\t\t<Header>\n\t\t\t<Title> {0} </Title>\n\t\t</Header>\n'.format(str_data_pdbi['ENCOMPASS']['name'])
    txt += '\t\t<General>\n'
    txt += '\t\t\t<PDBCode> {0} </PDBCode>\n'.format(pdbi)
    txt += '\t\t\t<Method> {0} </Method>\n'.format(str_data_pdbi['FROM_PDB']['experimental_method'])
    txt += '\t\t\t<Size value="{0}" />\n'.format(str_data_pdbi['ENCOMPASS']['structure']['chains'][ich]['n_of_residues'])
    ### WARNING this line is not general enough - TODO, figure out what else it needs
    txt += '\t\t\t<ChainUniProtAccs {0} />\n'.format(' '.join(['code'+str(i+1)+'='+'"'+ x.split("'")[1] +'"' for i,x in enumerate([str_data_pdbi['ENCOMPASS']['structure']['chains'][ich]['UniProt_acc']]) if len(x.split("'"))>1]))
    ### /WARNING
    txt += '\t\t\t<TMDomains> {0} </TMDomains>\n'.format(len(str_data_pdbi['ENCOMPASS']['structure']['chains'][ich]['TM_regions']['TM_regions_extrema']))
    txt += '\t\t\t<Sequence name=">{0}_{1}:">\n'.format(pdbi, ch)
    txt += '\t\t\t\t<Seq> {0} </Seq>\n'.format(str_data_pdbi['ENCOMPASS']['structure']['chains'][ich]['sequence'])
    txt += '\t\t\t</Sequence>\n'

    ### TM Domains
    for iextr, extr in enumerate(
            str_data_pdbi['ENCOMPASS']['structure']['chains'][ich]['TM_regions']['TM_regions_extrema']):
        txt += '\t\t\t\t<TMDomain ID="{0}">\n'.format(iextr)
        txt += '\t\t\t\t\t<TMDRange> {0} - {1} </TMDRange>\n'.format(extr[0][0][0], extr[1][0][0])
        txt += '\t\t\t\t</TMDomain>\n'

    txt += '\t\t\t<GeneralStructure>\n'
    txt += '\t\t\t\t<DownloadPDBFile> {0} </DownloadPDBFile>\n'.format(locations['FSYS']['whole'] + pdbi + '_enc.pdb')
    txt += '\t\t\t\t<ImageFile> {0} </ImageFile>\n'.format(locations['FSYS']['images_figs_chain'] + pdbi_ch + '.png')
    txt += '\t\t\t\t<T3dmolFile> {0} </T3dmolFile>\n'.format(locations['FSYS']['3dmol_whole'] + pdbi + '.json')
    txt += '\t\t\t</GeneralStructure>\n'
    txt += '\t\t\t<PDB_URL> http://www.rcsb.org/pdb/explore/explore.do?structureId={0} </PDB_URL>\n'.format(pdbi.upper())
    txt += '\t\t\t<OPM_URL> http://opm.phar.umich.edu/protein.php?pdbid={0} </OPM_URL>\n'.format(pdbi)
    txt += '\t\t\t<PDBTM_URL> http://pdbtm.enzim.hu/data/database/{0}/{1}.trpdb.gz </PDBTM_URL>\n'.format(pdbi[1:3], pdbi)
    if str_data_pdbi['ENCOMPASS']['resolution'] == "None":
        txt += '\t\t\t<Resolution> null </Resolution>\n'
    else:
        txt += '\t\t\t<Resolution> {0} </Resolution>\n'.format(str_data_pdbi['ENCOMPASS']['resolution'])
    txt += '\t\t</General>\n'


    # StructureInformation
    txt += '\t\t<StructureInformation>\n'
    txt += '\t\t\t<Member> {0}_{1} </Member>\n'.format(pdbi, ch)
    txt += '\t\t\t<Class> {0} </Class>\n'.format(str_data_pdbi['ENCOMPASS']['class'])
    txt += '\t\t\t<TMdomains> {0} </TMdomains>\n'.format(len(str_data_pdbi['ENCOMPASS']['structure']['chains'][ich]['TM_regions']['TM_regions_extrema']))
    txt += '\t\t\t<SequenceNeighbors> {0} </SequenceNeighbors>\n'.format(seqneighs[pdbi_ch])
    txt += '\t\t\t<StructureNeighbors> {0} </StructureNeighbors>\n'.format(strneighs[pdbi_ch])
    txt += '\t\t\t<TotalNeighbors> {0} </TotalNeighbors>\n'.format(totneighs[pdbi_ch])
    txt += '\t\t\t<SequenceNeighborsFile> {0} </SequenceNeighborsFile>\n'.format(locations['FSYS']['seqneighs'] + 'seqneigh_' + pdbi_ch + '.txt')
    txt += '\t\t\t<StructureNeighborsFile> {0} </StructureNeighborsFile>\n'.format(locations['FSYS']['strneighs'] + 'strneigh_' + pdbi_ch + '.txt')
    txt += '\t\t\t<TotalNeighborsFile> {0} </TotalNeighborsFile>\n'.format(locations['FSYS']['totneighs'] + 'totneigh_' + pdbi_ch + '.txt')
    txt += '\t\t\t<Images>\n'
    txt += '\t\t\t\t<RadialDistribution> {0} </RadialDistribution>\n'.format(locations['FSYS']['polar_figs'] + 'p_' + pdbi_ch + '.png')
    txt += '\t\t\t\t<RadialDistributionMap> {0} </RadialDistributionMap>\n'.format(locations['FSYS']['polar_maps'] + pdbi_ch + '_cmap.txt')
    txt += '\t\t\t\t<AlignmentsInDensityPlot> {0} </AlignmentsInDensityPlot>\n'.format(locations['FSYS']['densityscatter_figs'] + 'ds_' + pdbi_ch + '.png')
    txt += '\t\t\t\t<AlignmentsInDensityPlotMap> {0} </AlignmentsInDensityPlotMap>\n'.format(locations['FSYS']['densityscatter_maps'] + pdbi_ch + '_cmap.txt')
    txt += '\t\t\t\t<ResDistanceDistribution> {0} </ResDistanceDistribution>\n'.format(locations['FSYS']['distributions_figs'] + 'distr_' + pdbi_ch + '.png')
    txt += '\t\t\t\t<Topology> {0} </Topology>\n'.format(locations['FSYS']['topologies_figs'] + 'top_' + pdbi_ch + '.jpeg')
    txt += '\t\t\t</Images>\n'
    txt += '\t\t\t<AdditionalFiles>\n'
    # Structure alignments are present only if the chain is a representative one
    if ch in str_data[pdbi]['ENCOMPASS']['structure']['redundant_chains']:
        repr_ch = str_data[pdbi]['ENCOMPASS']['structure']['redundant_chains'][ch]
        if repr_ch in str_data[pdbi]['ENCOMPASS']['structure']['redundant_chains']:
            repr_ch = str_data[pdbi]['ENCOMPASS']['structure']['redundant_chains'][repr_ch]
    else:
        repr_ch = ch
    txt += '\t\t\t\t<ResDistanceDistribution> {0} </ResDistanceDistribution>\n'.format(locations['FSYS']['distributions_data'] + 'distr_' + pdbi + '_' + repr_ch + '.txt')
    txt += '\t\t\t\t<StructureAlignments> {0} </StructureAlignments>\n'.format(locations['FSYS']['stralns'] + 'stralns_' + pdbi + '_' + repr_ch + '.txt.gz') ###!
    txt += '\t\t\t\t<Estimators> {0} </Estimators>\n'.format(locations['FSYS']['estimators'] + 'est_' + pdbi + '_' + repr_ch + '.txt')
    txt += '\t\t\t</AdditionalFiles>\n'
    txt += '\t\t</StructureInformation>\n'

    print("\tMonomer =", monomeric)
    txt += symmetry_xml(locations, pdbi_ch, str_data_pdbi['ENCOMPASS']['class'], xml_type ='chains', symm_dic_pdbi=symm_dic_pdbi, transfer_dic=transfer_list, monomeric=monomeric)

    txt +=  '\t</Structure>\n'

    return txt


def create_template(locations, str_data, neighbor_counts, neighborlist_example_filenames, symm_data, v):

    seqneigh_counts, strneigh_counts, totneigh_counts = neighbor_counts
    seqneighs_fn_eg, strneighs_fn_eg, totneighs_fn_eg = neighborlist_example_filenames

    deletion_info = prepare_deletion_info(locations)

    wtxt = ctxt = '<Database>\n'
    for pdbi in str_data:

        print(pdbi)
        symm_results = symm_data # contains the symmetry information for the pdb and its chains
        if str_data[pdbi]['ENCOMPASS']['class'] == 'beta':
            beta_dictionary_filename = locations['FSYSPATH']['images_cesymm_symd_log'] + pdbi + "/" + pdbi + "_dic.pkl"
            if os.path.isfile(beta_dictionary_filename):
                print("\tReading beta barrel pickle for", pdbi, beta_dictionary_filename) if v else None
                with open(beta_dictionary_filename, "rb") as pkl_file:
                    symm_results = pkl.load(pkl_file)

        num_tm_ch = len([c for c in str_data[pdbi]['ENCOMPASS']['structure']['ktmchains'] if c != '-'])
        monomer = False  # initialize
        if num_tm_ch == 1:
            monomer = True
        print("\n", pdbi, num_tm_ch, "chains", monomer) if v else None

        chain_txt = ""
        consensus = {}
        for ich, tmch in enumerate(str_data[pdbi]['ENCOMPASS']['structure']['ktmchains']):
            if tmch == '-':
                continue
            pdbi_ch = pdbi + '_' + tmch

            transf_dict_path = locations['FSYSPATH']['transfer'] + pdbi_ch + "/" + pdbi_ch + "_transfer.pkl"
            if os.path.isfile(transf_dict_path):
                print("\tReading transfer pickle for", pdbi_ch, transf_dict_path) if v else None
                transfer_dic = pkl.load(open(transf_dict_path, "rb"))[pdbi_ch]
            else:
                print("\tNo transfer data file available for", pdbi_ch) if v else None
                transfer_dic = None

            if pdbi_ch in symm_results:
                print("\tComputing consensus for chain", pdbi_ch) if v else None
                reported_order, consensus_score = consensus_symmetry(symm_results[pdbi_ch], str_data[pdbi]['ENCOMPASS']['class'], transfer_list=transfer_dic)
            if monomer and pdbi in symm_results:
                print("\tComputing consensus for monomer", pdbi_ch) if v else None
                reported_order, consensus_score = consensus_symmetry(symm_results[pdbi], str_data[pdbi]['ENCOMPASS']['class'], transfer_list=transfer_dic)
            else:
                print("\tUpdating consensus for no symm_results", pdbi_ch) if v else None
                reported_order, consensus_score = consensus_symmetry({}, str_data[pdbi]['ENCOMPASS']['class'], transfer_list=transfer_dic)

            consensus[pdbi_ch] = {}
            consensus[pdbi_ch]['reported_order'] = reported_order
            consensus[pdbi_ch]['consensus_score'] = consensus_score
            print("\t\tFound reported_order, consensus_score", reported_order, consensus_score) if v else None

            if pdbi_ch not in seqneighs:
                seqneighs[pdbi_ch] = 0
                create_empty_neighborlist_file(locations, pdbi_ch, seqneigh_counts, seqneighs_fn_eg, "seq", v)
            if pdbi_ch not in strneighs:
                strneighs[pdbi_ch] = 0
                create_empty_neighborlist_file(locations, pdbi_ch, strneigh_counts, strneighs_fn_eg, "str", v)
            if pdbi_ch not in totneighs:
                totneighs[pdbi_ch] = 0
                create_empty_neighborlist_file(locations, pdbi_ch, totneigh_counts, totneighs_fn_eg, "tot", v)

            if pdbi_ch in symm_results:
                print("\tCreating chain xml data with symmetry results", pdbi_ch) if v else None
                ctxt += chain_xml_entry(pdbi, tmch, str_data[pdbi], symm_results[pdbi_ch], locations, neighbor_counts, transfer_list=transfer_dic, monomeric=monomer)
            elif monomer and pdbi in symm_results:            # for monomers use whole structure data, if it exists
                print("\tCreating chain xml data for monomer or single-TM chain using symmetry results for whole", pdbi, pdbi_ch) if v else None
                ctxt += chain_xml_entry(pdbi, tmch, str_data[pdbi], symm_results[pdbi], locations, neighbor_counts, transfer_list=transfer_dic, monomeric=monomer)
            else:
                print("\tCreating chain xml data without symmetry results", pdbi_ch) if v else None
                ctxt += chain_xml_entry(pdbi, tmch, str_data[pdbi], {}, locations, neighbor_counts, transfer_list=transfer_dic, monomeric=monomer)

        print("Whole chain:", pdbi)
        #  are creating the whole xml after the chain xml because we need the symmetry consensus score for each chain
        if pdbi in symm_results:
            print("\tCreating whole structure xml text with symmetry", pdbi) if v else None
            wtxt0, ok = struct_xml_entry(pdbi, str_data[pdbi], symm_results[pdbi], locations, neighbor_counts, deletion_info, consensus, monomeric=monomer)
        else:
            print("\tCreating whole structure xml text without symmetry", pdbi) if v else None
            wtxt0, ok = struct_xml_entry(pdbi, str_data[pdbi], {}, locations, neighbor_counts, deletion_info, consensus, monomeric=monomer)
        wtxt += wtxt0
        if not ok: # Don't continue if protein was in deletion set
            continue
        ctxt += chain_txt

    wtxt += '</Database>\n'
    ctxt += '</Database>\n'
    return wtxt, ctxt


def prepare_deletion_info(locations):

    deletion_codes = []
    print("Storing deletion codes from", locations['SYSFILES']['delcodes'])
    with open(locations['SYSFILES']['delcodes']) as dcf:
        for line in dcf:
            if not line.strip():
                continue
            dc, c = line.split('\t')
            #print(dc, c)
            deletion_codes.append((dc, c.strip()))

    deletion_messages = {}
    print("Storing deletion messages from", locations['SYSFILES']['delmsg'])
    with open(locations['SYSFILES']['delmsg']) as dmf:
        for line in dmf:
            if not line.strip():
                continue
            c, m = line.split('\t')
            #print(c.strip(), m)
            deletion_messages[c.strip()] = m
    deletion_info = (deletion_codes, deletion_messages)

    return deletion_info


def prepare_neighbor_info(locations):
    start_time = time.time()
    seqneighs = {os.path.basename(x[1][-10:-4]): x[0] for x in [y.decode('utf8').split() for y in subprocess.Popen(
        "for x in `ls {0}`; do wc -l {0}/${{x}}; done".format(locations['FSYSPATH']['seqneighs']),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).stdout] if "total" not in x[1]}
    print("obtained seqneigh neighbor counts")

    seqneighs_fn_templ = [os.path.basename(x[1]) for x in [y.decode('utf8').split() for y in subprocess.Popen(
        "for x in `ls {0}`; do wc -l {0}/${{x}}; done".format(locations['FSYSPATH']['seqneighs']),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).stdout] if "total" not in x[1]][0]
    print("filename example is", seqneighs_fn_templ)
    strneighs = {os.path.basename(x[1][-10:-4]): x[0] for x in [y.decode('utf8').split() for y in subprocess.Popen(
        "for x in `ls {0}`; do wc -l {0}/${{x}}; done".format(locations['FSYSPATH']['strneighs']),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).stdout] if "total" not in x[1]}
    print("obtained strneighs neighbor counts")

    strneighs_fn_templ = [os.path.basename(x[1]) for x in [y.decode('utf8').split() for y in subprocess.Popen(
        "for x in `ls {0}`; do wc -l {0}/${{x}}; done".format(locations['FSYSPATH']['strneighs']),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).stdout] if "total" not in x[1]][0]
    print("filename example is", strneighs_fn_templ)
    totneighs = {os.path.basename(x[1][-10:-4]): x[0] for x in [y.decode('utf8').split() for y in subprocess.Popen(
        "for x in `ls {0}`; do wc -l {0}/${{x}}; done".format(locations['FSYSPATH']['totneighs']),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).stdout] if "total" not in x[1]}
    print("obtained totneighs neighbor counts")

    totneighs_fn_templ = [os.path.basename(x[1]) for x in [y.decode('utf8').split() for y in subprocess.Popen(
        "for x in `ls {0}`; do wc -l {0}/${{x}}; done".format(locations['FSYSPATH']['totneighs']),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).stdout] if "total" not in x[1]][0]
    print("filename example is", totneighs_fn_templ)

    end_time = time.time()
    elapsed_time = end_time - start_time

    print("Neighbor counting took {:.2f} seconds".format(elapsed_time))
    neighbor_counts = seqneighs, strneighs, totneighs
    filename_examples = seqneighs_fn_templ, strneighs_fn_templ, totneighs_fn_templ
    return neighbor_counts, filename_examples


def create_empty_neighborlist_file(locations, pdbi_ch, neighbor_type, neighs_fn_templ, label, v):

    filename = locations['FSYSPATH'][neighbor_type] + neighs_fn_templ[:-10] + pdbi_ch + neighs_fn_templ[-4:]
    os.system("touch {0}".format(filename))
    print("EMPTY NEIGHBOR FILE CREATED:", label, filename) if v else None


def test_xml(xml_filename):
    import xml.etree.ElementTree as ET

    tree = ET.parse(xml_filename)
    root = tree.getroot()
    for str_cell in root.findall("Database/Structure"):
        for ch_cell in str_cell.findall("General/Chain"):
            print(ch_cell.find("NTMDomains").text)
    pass    


def webarchive(options, locations, str_data, all_neighs, all_neighs_fn_templ, table, symm_data, v):

    out_dir = options['PATHS'][('waout', 'wa_out_dir')]
    wa_whole_fn = out_dir + '/' + options['RUN'][('code', 'main_code')] + '_whole.xml' 
    wa_chains_fn = out_dir + '/' + options['RUN'][('code', 'main_code')] + '_chains.xml'
    print("Opening FILES to write: ", wa_whole_fn, wa_chains_fn)
    wtxt, ctxt = create_template(locations, str_data, all_neighs, all_neighs_fn_templ, symm_data, v)
    print("Creating FILES: ", wa_whole_fn, wa_chains_fn)
    with open(wa_whole_fn, 'w') as wf:
        wf.write(wtxt)
    with open(wa_chains_fn, 'w') as cf:
        cf.write(ctxt)


def store_neighbor_counts():
    global f, seqneighs, strneighs, totneighs, seqneighs_fn_example, strneighs_fn_example, totneighs_fn_example
    print("Preparing neighbor counts")
    neighbor_counts_filename = locations['FSYSPATH']['cache'] + 'neighborcounts.pkl'
    if os.path.exists(neighbor_counts_filename):
        print("WARNING: Checkpoint neighbor count file found, reading", neighbor_counts_filename)
        with open(neighbor_counts_filename, 'rb') as f:
            (all_counts, example_filenames) = pkl.load(f)
        print("\tloaded", neighbor_counts_filename)
        (seqneighs, strneighs, totneighs) = all_counts
        (seqneighs_fn_example, strneighs_fn_example, totneighs_fn_example) = example_filenames
        print(example_filenames)
    else:
        (seqneighs, strneighs, totneighs), (seqneighs_fn_example, strneighs_fn_example, totneighs_fn_example) = prepare_neighbor_info(locations)
        with open(neighbor_counts_filename, 'wb') as f:
            pickle.dump(
                ((seqneighs, strneighs, totneighs), (seqneighs_fn_example, strneighs_fn_example, totneighs_fn_example)),f)
            print("Checkpoint neighbor count data saved for future iterations:", neighbor_counts_filename)


def read_data_structure():
    global str_data
    str_data_filename = locations['FSYSPATH']['cache'] + 'str_data_straln.pkl'
    str_template_filename = locations['SYSFILES']['data_structure_template']
    print("Reading", str_data_filename)
    print("Reading", str_template_filename)
    str_data = read_checkpoint(str_data_filename, str_template_filename)

    #str_data = pkl.load(open("/data/TMB-CSB/Forrest/0_working/encompass_project/checks/create_xml_tmp/str_data_failed_stucts.pkl", "rb"))
    load_time = time.time()
    print("Finished reading at", format_time(load_time))


def read_mssd_data():
    global mssd_data
    mssd_data_filename = locations['FSYSPATH']['mssd'] + "symmetry_results.pkl"
    print("Reading", mssd_data_filename)
    mssd_data = pkl.load(open(mssd_data_filename, "rb"))
    load_time = time.time()
    print("Finished at", format_time(load_time))

if __name__ == "__main__":

    start_time = time.time()
    print("Running create_xml.py at ", format_time(start_time))
    options, locations = initialize_repository()
    verbose_level = 1

    read_data_structure()
    # only uncomment next statement if don't want eliminated entries to have Info statements
    #str_data = {x : str_data[x] for x in str_data if not ('eliminated' in str_data[x]['status'])}

    if options['PARAMETERS'][('squash', 'make_squash_filesystem')]:
        print("Creating squash filesystem")
        squash_webarchive(options, locations, generate=True) # TODO - exclude files from exclusions list?
    else:
        print("Not creating squash filesystem, because squash is set to False")

    read_mssd_data()
    store_neighbor_counts()

    webarchive(options, locations, str_data, (seqneighs, strneighs, totneighs),
               (seqneighs_fn_example, strneighs_fn_example, totneighs_fn_example), {}, mssd_data, verbose_level)

    end_time = time.time()
    print("Finished at", format_time(end_time))

