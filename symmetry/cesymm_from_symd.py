# Name: cesymm_from_symd.py
# Author: Antoniya A. Aleksandrova
# Language: Python 3.5+
# Description: Submit fragments of a protein for analysis with CE-Symm based on the results from SymD

import os
import sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from initialize_repository import initialize_repository
from symmetry_exec_functions import *

start_time = time.time()

def cesymm_from_symd(locations, options, inputf, pdb_suff, mem_req_limit = 20000, run_id="0000"):
    pdb = inputf[0:4]
    chain = inputf[5:6]
    seq = inputf[7:]

    log_file = open(locations['FSYSPATH']['cesymm_from_symd_log'] + inputf + '_fragment.log', 'w')
    save_out = sys.stdout
    save_err = sys.stderr
    sys.stdout = sys.stderr = log_file

    # Set options
    if os.path.isfile(locations['FSYSPATH']['whole'] + pdb + pdb_suff +'.pdb')==False:
        raise SystemExit("{0} does not exist.".format(locations['FSYSPATH']['whole'] + pdb + pdb_suff +'.pdb'))
    wkdir = locations['FSYSPATH']['symmtemp'] + "cesymm_from_symd/"  + inputf + "/"
    if not os.path.isdir(wkdir):
        os.mkdir(wkdir)

    print("In: " + inputf+" submitted")

    out_dir=locations['FSYSPATH']['cesymm_from_symd'] + inputf + "/"
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    oriented=locations['FSYSPATH']['whole']+ pdb + pdb_suff + ".pdb"
    winsize=8
    scorethreshold=0.4
    minreplen=15
    maxsymlev=0
    minlength=20 #20 residues
    maxrmsd=99

    flds=seq.split('-')
    first_resid=flds[0]
    last_resid=flds[1]

    num=strip_tm_chains(wkdir, pdb + "_" + chain, locations['FSYSPATH']['whole'] + pdb + pdb_suff + '.pdb', chain)  #extract chain
    os.rename(wkdir + pdb + "_" + chain + "_tmp.pdb", wkdir + pdb + "_" + chain + ".pdb")
    oriented=wkdir + pdb + "_" + chain + ".pdb"
    if num == 0:
        raise SystemExit("{0} has 0 TM atoms.".format(locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb'))
    if num > mem_req_limit:
        cesymm_dir = options['ALL']['sigcesymm']
    else:
        cesymm_dir = options['ALL']['sigcesymmsmall']
    symd_dir=options['ALL']['sigsymd']
    residues_extract(wkdir, pdb + "." + chain, oriented, first_resid, last_resid) # extract sequence
    os.remove(wkdir + pdb + "_" + chain + ".pdb")
    oriented = wkdir + pdb + "." + chain + "_" + first_resid + "-" + last_resid + ".pdb"

    print("In: CE-Symm")
    ce_symm(wkdir,cesymm_dir,out_dir,pdb+"."+chain+"_"+first_resid+"-"+last_resid,oriented,winsize,scorethreshold,minreplen,maxsymlev,maxrmsd)
    num=num_atoms_pdb(oriented)
    print("In: Symd")
    symd(wkdir,symd_dir,out_dir,pdb+"."+chain+"_"+first_resid+"-"+last_resid,oriented,num)

    #Clean-up
    print("{0} completed.".format(inputf))
    os.remove(oriented)
    shutil.rmtree(wkdir)
    print("Time[s]:%s" % (time.time() - start_time))
    print("Date: %s" % time.asctime())
    print("Run id: %s" % run_id)
    sys.stdout = save_out
    sys.stderr = save_err
    log_file.close()
    return inputf

if __name__ == "__main__":
    if len(sys.argv) < 5:
        raise SystemExit("syntax: %s <locations_path> <pdbcode_chain_seq> <pdb_suffix> <run_id>" % sys.argv[0])
    else:
        locations_path = sys.argv[1].strip()
        inputf = sys.argv[2]
        pdb_suff = sys.argv[3].strip()
        run_id = sys.argv[4].strip()

    options, locations = initialize_repository(main_path=locations_path, instr_filename_path="EncoMPASS_options_relative_to__instruction_file_benchmark.txt")

    cesymm_from_symd(locations, options, inputf, pdb_suff, mem_req_limit=20000, run_id = run_id)