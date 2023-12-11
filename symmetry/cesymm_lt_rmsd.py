# Name: cesymm_lt_rmsd.py
# Author: Antoniya A. Aleksandrova
# Language: Python 3.5+
"""
Executes a CE-Symm run with unrefinedscore=0.2 and maxrmsd (=2.5 or 3) on each chain in the database
"""

import os
import time
import sys
import shutil
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from initialize_repository import initialize_repository
from symmetry_exec_functions import *
start_time = time.time()

def cesymm_lt_rmsd(locations, options, inputf, pdb_suff, maxrmsd, run_id="0000"):
    maxrmsd = float(maxrmsd)
    if maxrmsd.is_integer():
        rmsd_descriptor = '_rmsd_' + str(int(maxrmsd)) # _rmsd_3
    else:
        rmsd_descriptor = '_rmsd_' + str(int(maxrmsd)) + '_' + str(int((maxrmsd - int(maxrmsd))*10)) # _rmsd_2_5

    log_file = open(locations['FSYSPATH']['cesymm' + rmsd_descriptor + '_log'] + inputf + rmsd_descriptor + '.log', 'w')
    save_out = sys.stdout
    save_err = sys.stderr
    sys.stdout = sys.stderr = log_file

    # Set options
    if os.path.isfile(locations['FSYSPATH']['whole']+inputf+pdb_suff+'.pdb')==False:
        raise SystemExit("{0} does not exist.".format(locations['FSYSPATH']['whole']+inputf+pdb_suff+'.pdb'))
    wkdir = locations['FSYSPATH']['symmtemp'] + "cesymm" + rmsd_descriptor + "/" + inputf + "/"
    if not os.path.isdir(wkdir):
        os.mkdir(wkdir)
    cesymm_exe = options['ALL']['sigcesymmsmall']

    out_dir = locations['FSYSPATH']['cesymm' + rmsd_descriptor] + inputf + "/"
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    tm_archive = locations['FSYSPATH']['symmetry'] + '.tm_archive.pkl'
    a = pkl.load(open(tm_archive, 'rb'))
    winsize=16
    scorethreshold=0.2
    minreplen=15
    maxsymlev=0
    minlength=20 # minimum number of residues required for CE-Symm submission

    print("In: "+inputf+" submitted")

    # Find TM chains
    chain = ""
    if 'tmchains' in a[inputf] and len(a[inputf]['tmchains'])>0:
        chain = ';'.join(sorted(a[inputf]['tmchains'])) + ';'  # tm chains as a string A;B;C;..
        print(chain)
    else:
        raise SystemExit("{0} has no TM chains.".format(locations['FSYSPATH']['whole']+inputf+pdb_suff+'.pdb'))

    if len(chain)>2:
        monomer='N'
    else:
        monomer='Y'
    all_chain=int(len(chain)/2)  # the number of tm chains

    # Run on whole protein with TM chains only

    num=strip_tm_chains(wkdir,inputf,locations['FSYSPATH']['whole']+inputf+pdb_suff+'.pdb',chain)
    os.rename(wkdir+inputf+"_tmp.pdb",wkdir+inputf+".pdb")
    oriented=wkdir+inputf+".pdb"
    if num==0:
        raise SystemExit("{0} has 0 TM atoms.".format(locations['FSYSPATH']['whole']+inputf+'_enc.pdb'))

    if monomer=='Y':
        ce_symm(wkdir,cesymm_exe,out_dir,inputf,oriented,winsize,scorethreshold,minreplen,maxsymlev,maxrmsd)

    # Run on each TM chain of oligomer
    if monomer=='N':
        # Check the order of chains
        oriented_old=oriented

        # Get each chain
        for i in range(0,all_chain):
            ch_id=chain[i*2]
            print("In: Chain "+ch_id)
            oriented=oriented_old
            get_chain(oriented,ch_id,inputf,wkdir)
            oriented=wkdir+inputf+"_"+ch_id+".pdb"
            res_num=resnum(oriented)
            if res_num>minlength:
                ce_symm(wkdir,cesymm_exe,out_dir,inputf+"_"+ch_id,oriented,winsize,scorethreshold,minreplen,maxsymlev,maxrmsd)
            else:
                print("Out: Chain {0} has too few residues ({1}).".format(inputf+"_"+ch_id,str(res_num)))

    #Clean-up
    print("{0} completed.".format(inputf))
    shutil.rmtree(wkdir)
    print("Time[s]:%s" % (time.time() - start_time))
    print("Date: %s" % time.asctime())
    print("Run id: %s" % run_id)
    sys.stdout = save_out
    sys.stderr = save_err
    log_file.close()
    return inputf

if __name__ == "__main__":
    if len(sys.argv) < 6:
        raise SystemExit("syntax: %s <locations_path> <pdb_options:e.g 2wcd_ordered_df> <pdb_suff> <run_id> <max rmsd>" % sys.argv[0])
    else:
        locations_path = sys.argv[1].strip()
        inputf = sys.argv[2].strip()
        pdb_suff = sys.argv[3].strip()
        run_id = sys.argv[4].strip()
        maxrmsd = sys.argv[5].strip()

    options, locations = initialize_repository(main_path=locations_path, instr_filename_path="EncoMPASS_options_relative_to__instruction_file_benchmark.txt")

    cesymm_lt_rmsd(locations, options, inputf, pdb_suff, maxrmsd, run_id)
