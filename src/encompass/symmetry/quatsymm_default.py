# Name: quatsymm.py
# Author: Antoniya A. Aleksandrova
# Date: 9 May 2022
# Updated:
# Language: Python 3.8
"""
Executes QuatSymm on a single structure
"""

from encompass.pipeline.initialize_repository import initialize_repository
from encompass.symmetry.symmetry_exec_functions import *

start_time = time.time()


def quatsymm_default(locations: dict, options: dict, inputf: str, pdb_suff: str, run_id: str = "0000") -> str:
    """
    Submit a structure to QuatSymm

    :param locations: dictionary with the relevant paths (log, pdb, wkdir, output)
    :param options:
    :param inputf: name of input structure
    :param pdb_suff: file extension used in file name of structure
    :param run_id: code identifying this specific submission
    :return: name of input structure (to denote complete status)
    """
    # Set options

    log_file = open(locations['FSYSPATH']['quatsymm_log'] + inputf + '_qsymm.log', 'w')
    save_out = sys.stdout
    save_err = sys.stderr
    sys.stdout = sys.stderr = log_file

    print("Run id: %s" % run_id)

    if not os.path.isfile(locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb'):
        raise SystemExit("{0} does not exist.".format(locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb'))
    wkdir = locations['FSYSPATH']['symmtemp'] + "quatsymm/" + inputf + "/"
    if not os.path.isdir(wkdir):
        os.mkdir(wkdir)
    out_dir = locations['FSYSPATH']['quatsymm'] + inputf + "/"
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    quatsymm_exe = options['ALL']['sigquatsymm']
    # quatsymm_exe = "/data/Encompass/software/quatsymm-2.2.3-SNAPSHOT/runQuatSymm.sh"
    #quatsymm_exe = "/data/TMB-CSB/Aleksandrova/software/quatsymm-2.2.3-SNAPSHOT/runQuatSymm.sh"

    suffix = "_qsymm"

    tm_archive = tm_archive = locations['SYSFILES']['tm_archive']
    a = pkl.load(open(tm_archive, 'rb'))
    chain = ';'.join(sorted(a[inputf]['tmchains']))
    print(chain, locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb')
    num = strip_tm_chains(wkdir, inputf, locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb', chain)
    os.rename(wkdir + inputf + "_tmp.pdb", wkdir + inputf + ".pdb")
    oriented = wkdir + inputf + ".pdb"

    quatsymm(quatsymm_exe, out_dir, inputf, oriented, suffix)

    # Clean-up
    print("{0} completed.".format(inputf))
    # shutil.rmtree(wkdir)
    print("Time[s]:%s" % (time.time() - start_time))
    print("Date: %s" % time.asctime())

    sys.stdout = save_out
    sys.stderr = save_err
    log_file.close()

    return inputf


if __name__ == "__main__":
    if len(sys.argv) < 5:
        raise SystemExit("syntax: %s <locations_path> <pdbcode> <pdb_suffix> <run_id>" % sys.argv[0])
    else:
        locations_path = sys.argv[1].strip()
        struct = sys.argv[2][0:4]
        pdb_extension = sys.argv[3].strip()
        run_name = sys.argv[4].strip()

    #    currentdir = os.path.dirname(os.path.realpath(__file__))
    #    parentdir = os.path.dirname(currentdir)
    #    sys.path.append(parentdir)
    # locations_path = os.path.expanduser('~') + "/" + "EncoMPASS_1.1_2018_03_22/"
    opt, loc = initialize_repository(main_path=locations_path,
                                     instr_filename_path="EncoMPASS_options_relative_to__instruction_file_2023_12_13_shulien.txt")
    # pdb_suff = '_enc'
    if not os.path.isdir(loc['FSYSPATH']['symmtemp'] + "quatsymm/"):
        os.mkdir(loc['FSYSPATH']['symmtemp'] + "quatsymm/")
    quatsymm_default(loc, opt, struct, pdb_extension, run_id=run_name)
