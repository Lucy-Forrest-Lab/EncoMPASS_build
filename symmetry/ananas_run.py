# Name: ananas_run.py
# Author: Antoniya A. Aleksandrova
# Date: 30 Jan 2019
# Language: Python 3.5
# Description: Executes an AnAnaS run with default options and a protein complex

import os
import sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from initialize_repository import initialize_repository
from symmetry_exec_functions import *
import time

start_time = time.time()


def ananas_run(locations, inputf, pdb_suff, run_id="0000"):

    # Set options
    log_file = open(locations['FSYSPATH']['ananas_log'] + inputf + '_ananas.log', 'w')
    save_out = sys.stdout
    save_err = sys.stderr
    sys.stdout = sys.stderr = log_file

    if not os.path.isfile(locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb'):
        raise SystemExit("{0} does not exist.".format(locations['FSYSPATH']['whole'] + pdb_suff + '.pdb'))
    wkdir = locations['FSYSPATH']['symmtemp'] + "ananas/"

    ananas_path = options['ALL']['sigcesymm']

    print("In: " + inputf + " submitted")
    out_dir = locations['FSYSPATH']['ananas']
    tm_archive = locations['FSYSPATH']['symmetry'] + '.tm_archive.pkl'
    a = pkl.load(open(tm_archive, 'rb'))

    # Find TM chains
    chain = ""
    chain = ';'.join(sorted(a[inputf]['tmchains'])) + ';'  # tm chains as a string A;B;C;...
    if chain == "":
        raise SystemExit("{0} has no TM chains.".format(locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb'))
    else:
        print(chain)
    if len(chain) > 2:
        # Run on whole protein with TM chains only
        num = strip_tm_chains(wkdir, inputf, locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb', chain)
        os.rename(wkdir + inputf + "_tmp.pdb", wkdir + inputf + ".pdb")
        oriented = wkdir + inputf + ".pdb"
        if num == 0:
            raise SystemExit("{0} has 0 TM atoms.".format(locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb'))

        ananas(wkdir, ananas_path, out_dir, inputf, oriented)
        os.remove(wkdir + inputf + ".pdb")
    else:
        print("Single chain: structure will be skipped.")

    # Clean-up
    print("{0} completed.".format(inputf))
    print("Time[s]:%s" % (time.time() - start_time))
    print("Date: %s" % time.asctime())
    print("Run id: %s" % run_id)
    sys.stdout = save_out
    sys.stderr = save_err
    log_file.close()
    return inputf


if __name__ == "__main__":
    if len(sys.argv) < 5:
        raise SystemExit("syntax: %s <locations_path> <pdb_id:e.g 2wcd> <pdb_suff> <run_id>" % sys.argv[0])
    else:
        locations_path = sys.argv[1].strip()
        inputf = sys.argv[2].strip()[0:4]
        pdb_suff = sys.argv[3].strip()
        run_id = sys.argv[4].strip()

    options, locations = initialize_repository(main_path=locations_path, instr_filename_path="EncoMPASS_options_relative_to__instruction_file.txt")

    ananas_run(locations, inputf, pdb_suff)
