# Name: cesymmm_quat_minlen.py
# Author: Antoniya A. Aleksandrova
# Language: Python 3.5+
# Description: Run a structure with CE-Symm and options, eg. from quat_symm_minlen_small_pdbs.log

import os
import sys
import shutil
import argparse
import time
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from initialize_repository import initialize_repository
from symmetry_exec_functions import *
start_time = time.time()

def cesymm_quat_minlen(locations, options, inputf_options, pdb_suff, mem_req_limit = 20000, run_id="0000"):
    log_file = open(locations['FSYSPATH']['cesymm_quat_minlen_log'] + inputf_options + '_quat.log', 'w')
    save_out = sys.stdout
    save_err = sys.stderr
    sys.stdout = sys.stderr = log_file
    print("Run id: %s" % run_id)

    # Set options
    if not os.path.isfile(locations['FSYSPATH']['cesymm_quat_minlen_log'] + 'quat_symm_minlen_pdbs.pkl'):
        raise SystemExit("Dictionary with running values " + locations['FSYSPATH']['cesymm_quat_minlen_log'] + 'quat_symm_minlen_pdbs.pkl not found.\n')
    inp_dic_full = pkl.load(open(locations['FSYSPATH']['cesymm_quat_minlen_log'] + 'quat_symm_minlen_pdbs.pkl','rb'))
    inp_dic = inp_dic_full[inputf_options]
    wkdir = locations['FSYSPATH']['symmtemp'] + "cesymm_quat_minlen/" + inp_dic['name'] + "/"
    if os.path.isdir(wkdir)==False:
        os.mkdir(wkdir)
    chain_ordered_dir=locations['FSYSPATH']['symmtemp'] + "chain_ordered_str/"
    out_dir=locations['FSYSPATH']['cesymm_quat_minlen'] + inp_dic['name'] + "/"
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    cesymm_exe = options['ALL']['sigcesymm']
    cesymm_exe_small = options['ALL']['sigcesymmsmall']

    tm_archive = locations['FSYSPATH']['symmetry'] + '.tm_archive.pkl'
    a = pkl.load(open(tm_archive, 'rb'))
    suffix="_quat"

    print("In: " + inp_dic['name'] + " submitted")
    inputf = inp_dic['name']
    pdb=inputf[0:4]

    if '_lt_mlev' in inputf or '_df_mlev' in inputf:
        inputf=inputf[:-8]
        inputf_options = inputf + '_' + inp_dic['mode']
    elif '_mlev' in inputf:
        raise SystemExit("mlev should be accompanied by df or lt")
    elif '_lt' in inputf or '_df' in inputf:
        inputf = inputf[:-3]

    if inp_dic['mode']=='df':
        winsize=8
        scorethreshold=0.4
    elif inp_dic['mode']=='lt':
        winsize=16
        scorethreshold=0.2
    minreplen=15 # this is the default value; keep this default to avoid issues, eg. 2vpw
    order=0 # no order enforced
    seed=inp_dic['seed']
    maxsymlev=int(inp_dic['required_symm_levels'])

    if a[pdb]['natoms'] > mem_req_limit:
        cesymm_path = cesymm_exe
    else:
        cesymm_path = cesymm_exe_small

    try:
        if len(inputf)==4: # standard pdb
            chain = ';'.join(sorted(a[inputf]['tmchains'])) + ';'  # tm chains as a string A;B;C;..
            num=strip_tm_chains(wkdir,pdb,locations['FSYSPATH']['whole']+pdb+pdb_suff+'.pdb',chain)
            os.rename(wkdir+pdb+"_tmp.pdb", wkdir+inputf+".pdb")
            oriented=wkdir+inputf+".pdb"
            ce_symm_minlen(cesymm_path,out_dir,inputf_options,oriented,winsize,scorethreshold,minreplen,maxsymlev,order,seed,suffix)
            os.remove(oriented)
        if len(inputf)>6: # re-ordered chains pdb
            oriented=chain_ordered_dir+inputf+".pdb"
            ce_symm_minlen(cesymm_path,out_dir,inputf_options,oriented,winsize,scorethreshold,minreplen,maxsymlev,order,seed,suffix)
        print(inp_dic['name'] +" completed.")
    except:
        print(sys.exc_info()[0])
        print("Error with " + inp_dic['name'])

    shutil.rmtree(wkdir)
    print("Time[s]:%s" % (time.time() - start_time))
    print("Date: %s" % time.asctime())
    sys.stdout = save_out
    sys.stderr = save_err
    log_file.close()
    return inp_dic['name']

if __name__ == "__main__":
    if len(sys.argv) < 5:
        raise SystemExit("syntax: %s <locations_path> <pdb_options:e.g 2wcd_ordered_df> <pdb_suff> <run_id>" % sys.argv[0])
    else:
        locations_path = sys.argv[1].strip()
        inputf = sys.argv[2].strip()
        pdb_suff = sys.argv[3].strip()
        run_id = sys.argv[4].strip()

    options, locations = initialize_repository(main_path=locations_path, instr_filename_path="EncoMPASS_options_relative_to__instruction_file_benchmark.txt")

    cesymm_quat_minlen(locations, options, inputf, pdb_suff, mem_req_limit = 20000, run_id = run_id)