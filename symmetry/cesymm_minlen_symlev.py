# Name: cesymmm_minlen_symlev.py
# Author: Antoniya A. Aleksandrova
# Language: Python 3.5+
# Description:


import os
import sys
import shutil
import time
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from initialize_repository import initialize_repository
from symmetry_exec_functions import *
start_time = time.time()

def cesymm_minlen_symlev(locations, options, inputf_options, pdb_suff, mem_req_limit=20000, run_id="0000"):

    log_file = open(locations['FSYSPATH']['cesymm_minlen_log'] + inputf_options + '_minlen.log', 'w')
    save_out = sys.stdout
    save_err = sys.stderr
    sys.stdout = sys.stderr = log_file

    # Set options
    if not os.path.isfile(locations['FSYSPATH']['cesymm_minlen_log'] + 'minlen_violations_pdbs.pkl'):
        raise SystemExit("Dictionary with running values " + locations['FSYSPATH']['cesymm_minlen_log'] + 'minlen_violations_pdbs.pkl not found.\n')
    inp_dic_full = pkl.load(open(locations['FSYSPATH']['cesymm_minlen_log'] + 'minlen_violations_pdbs.pkl','rb'))
    inp_dic = inp_dic_full[inputf_options]
    wkdir = locations['FSYSPATH']['symmtemp'] + "cesymm_minlen/" + inp_dic['name'] + "/"
    if not os.path.isdir(wkdir):
        os.mkdir(wkdir)
    chain_ordered_dir=locations['FSYSPATH']['symmtemp']+"chain_ordered_str/"
    out_dir=locations['FSYSPATH']['cesymm_minlen'] + inp_dic['name'] + "/"
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    cesymm_exe = options['ALL']['sigcesymm']
    cesymm_exe_small = options['ALL']['sigcesymmsmall']

    tm_archive = locations['FSYSPATH']['symmetry'] + '.tm_archive.pkl'
    a = pkl.load(open(tm_archive, 'rb'))
    suffix="_mlev"
    order=0

    print("In: " + inp_dic['name'] + " submitted")
    pdb=inp_dic['name'][0:4]
    inputf = inp_dic['name'][:-3] # eg. 2wcd_ordered
    inputf_options = inp_dic['name'] # eg. 2wcd_ordered_df
    if '_df' in inputf_options:
        winsize=8
        scorethreshold=0.4
    elif '_lt' in inputf_options:
        winsize=16
        scorethreshold=0.2
    else:
        raise SystemExit('In input "%s" there is no specified CE-Symm mode (df/lt)' %(inputf_options))

    minreplen=15 # this is the default value; keep this default to avoid issues, eg. 2vpw
    seed=inp_dic['seed']
    maxsymlev=int(inp_dic['cesymm_symm_levels'])-1
    try:
        if len(inputf)==6: # chains
            chain=inputf[5:6]
            num=strip_tm_chains(wkdir,inputf,locations['FSYSPATH']['whole']+pdb+pdb_suff+'.pdb',chain)
            if num > mem_req_limit:
                cesymm_dir = cesymm_exe
            else:
                cesymm_dir = cesymm_exe_small
            os.rename(wkdir+inputf+"_tmp.pdb",wkdir+inputf+".pdb")
            oriented=wkdir+inputf+".pdb"
            ce_symm_minlen(cesymm_dir,out_dir,inputf_options,oriented,winsize,scorethreshold,minreplen,maxsymlev,order,seed,suffix)
            os.remove(oriented)
        else:
            if len(inputf)==4: # standard pdb
                # Find TM chains
                chain = ';'.join(sorted(a[inputf]['tmchains'])) + ';'  # tm chains as a string A;B;C;..
                num=strip_tm_chains(wkdir,pdb,locations['FSYSPATH']['whole']+pdb+pdb_suff+'.pdb',chain)
                if num > mem_req_limit:
                    cesymm_dir = cesymm_exe
                else:
                    cesymm_dir = cesymm_exe_small
                os.rename(wkdir+pdb+"_tmp.pdb",wkdir+inputf+".pdb")
                oriented=wkdir+inputf+".pdb"
                ce_symm_minlen(cesymm_dir,out_dir,inputf_options,oriented,winsize,scorethreshold,minreplen,maxsymlev,order,seed,suffix)
                os.remove(oriented)
            if len(inputf)>6: # re-ordered chains pdb
                oriented=chain_ordered_dir+inputf+".pdb"
                if a[pdb]['natoms'] > mem_req_limit:
                    cesymm_dir = cesymm_exe
                else:
                    cesymm_dir = cesymm_exe_small
                ce_symm_minlen(cesymm_dir,out_dir,inputf_options,oriented,winsize,scorethreshold,minreplen,maxsymlev,order,seed,suffix)
        print(inputf_options+" completed.")
    except:
        print(sys.exc_info()[0])
        print("Error with "+inputf_options + ': ' + sys.exc_info()[0])

    shutil.rmtree(wkdir)
    print("Time[s]:%s" % (time.time() - start_time))
    print("Date: %s" % time.asctime())
    print("Run id: %s" % run_id)
    sys.stdout = save_out
    sys.stderr = save_err
    log_file.close()
    return inputf_options

if __name__ == "__main__":
    if len(sys.argv) < 5:
        raise SystemExit("syntax: %s <locations_path> <pdb_options:e.g 2wcd_ordered_df> <pdb_suff> <run_id>" % sys.argv[0])
    else:
        locations_path = sys.argv[1].strip()
        inputf = sys.argv[2].strip()
        pdb_suff = sys.argv[3].strip()
        run_id = sys.argv[4].strip()

    options, locations = initialize_repository(main_path=locations_path, instr_filename_path="EncoMPASS_options_relative_to__instruction_file_benchmark.txt")

    cesymm_minlen_symlev(locations, options, inputf, pdb_suff, mem_req_limit = 20000, run_id = run_id)
