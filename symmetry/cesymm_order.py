# Name: cesymmm_order.py
# Author: Antoniya A. Aleksandrova
# Language: Python 3.5+
# Description:

import os
import sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from initialize_repository import initialize_repository
from symmetry_exec_functions import *
import shutil
import time
start_time = time.time()


def cesymm_order(locations, options, out_struct_name, pdb_suff, mem_req_limit=20000, run_id="0000"): # inp_dic['name'] is for example "3lbw_ordered_df";
    """

    :param locations: dictionary with directories of the database
    :param pdb_suff: eg. _enc
    :param out_struct_name: eg. 2wcd_ordered_df
    :param mem_req_limit: the threshold number of atoms which should trigger a switch to a version of runCESymm.sh with larger Java heap allowance
    :param run_id: unique identifier of the run
    :return: the name of the completed structure (out_struct_name)
    """
    log_file = open(locations['FSYSPATH']['cesymm_order_log'] + out_struct_name + '_ord.log', 'w')
    save_out = sys.stdout
    save_err = sys.stderr
    sys.stdout = sys.stderr = log_file

    # Set options
    if not os.path.isfile(locations['FSYSPATH']['cesymm_order_log'] + 'cs_symd_order.pkl'):
        raise SystemExit("Dictionary with running values " + locations['FSYSPATH']['cesymm_order_log'] + 'cs_symd_order.pkl not found.\n')
    inp_dic_full = pkl.load(open(locations['FSYSPATH']['cesymm_order_log'] + 'cs_symd_order.pkl','rb'))
    inp_dic = inp_dic_full[out_struct_name]
    wkdir=locations['FSYSPATH']['symmtemp'] + "cesymm_order/"
    if not os.path.isdir(wkdir):
        os.mkdir(wkdir)
    #if not os.path.isdir(wkdir + "temp_structures"):
    #    os.mkdir(wkdir + "temp_structures")

    print("In: " + inp_dic['name'] + " submitted")
    chain_ordered_dir=locations['FSYSPATH']['symmtemp']+"chain_ordered_str/"
    out_dir=locations['FSYSPATH']['cesymm_order'] + inp_dic['name'] + "/"
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    tm_archive = locations['FSYSPATH']['symmetry'] + '.tm_archive.pkl'
    a = pkl.load(open(tm_archive, 'rb'))
    suffix="_ord"

    pdb = inp_dic['name'][0:4] # eg. 2wcd
    inputf = inp_dic['name'][:-3] # eg. 2wcd_ordered
    if '_df' in out_struct_name:
        winsize=8
        scorethreshold=0.4
    elif '_lt' in out_struct_name:
        winsize=16
        scorethreshold=0.2
    else:
        raise SystemExit('In input %s there is no specified CE-Symm mode (df/lt)' % (out_struct_name))

    minreplen = 15 # this is the default value; keep this default to avoid issues, eg. 2vpw
    seed = inp_dic['seed']
    maxsymlev = 0 # default
    order = int(inp_dic['symd_order'])
    cesymm_exe = options['ALL']['sigcesymm']
    cesymm_exe_small = options['ALL']['sigcesymmsmall']
    #cesymm_exe ="/data/TMB-CSB/Aleksandrova/software/cesymm-2.2.2/runCESymm.sh"
    #cesymm_exe_small ="/data/TMB-CSB/Aleksandrova/software/cesymm-2.2.2/runCESymm.sh"

    try:
        if len(inputf)==6: # chains
            if os.path.isfile(out_dir + out_struct_name + suffix + "_stdout.out"):
                pass
            chain=inputf[5:6]
            num=strip_tm_chains(wkdir,inputf,locations['FSYSPATH']['whole'] + pdb + pdb_suff + '.pdb',chain)
            if num > mem_req_limit:
                cesymm_dir = cesymm_exe
            else:
                cesymm_dir = cesymm_exe_small
            os.rename(wkdir + inputf + "_tmp.pdb", wkdir + inputf + ".pdb")
            oriented = wkdir + inputf + ".pdb"
            ce_symm_minlen(cesymm_dir, out_dir, out_struct_name, oriented, winsize, scorethreshold, minreplen, maxsymlev, order, seed, suffix)
            os.remove(oriented)
        else:
            if len(inputf)==4: # standard pdb
                if os.path.isfile(out_dir + out_struct_name + suffix + "_stdout.out"):
                    pass
                # Find TM chains
                chain=""
                chain = ';'.join(sorted(a[inputf]['tmchains'])) + ';'  # tm chains as a string A;B;C;..
                num=strip_tm_chains(wkdir,pdb,locations['FSYSPATH']['whole'] + pdb + pdb_suff + '.pdb',chain)
                if num > mem_req_limit:
                    cesymm_dir = cesymm_exe
                else:
                    cesymm_dir = cesymm_exe_small
                os.rename(wkdir+pdb+"_tmp.pdb",wkdir+inputf+".pdb")
                oriented=wkdir+inputf+".pdb"
                ce_symm_minlen(cesymm_dir, out_dir, out_struct_name, oriented, winsize, scorethreshold, minreplen, maxsymlev, order, seed, suffix)
                os.remove(oriented)
            if len(inputf)>6: # re-ordered chains pdb
                if os.path.isfile(out_dir + out_struct_name + suffix + "_stdout.out"):
                    pass
                if a[pdb]['natoms'] > mem_req_limit:
                    cesymm_dir = cesymm_exe
                else:
                    cesymm_dir = cesymm_exe_small
                oriented=chain_ordered_dir+inputf+".pdb"
                ce_symm_minlen(cesymm_dir, out_dir, out_struct_name, oriented, winsize, scorethreshold, minreplen, maxsymlev, order, seed, suffix)
        print(out_struct_name + " completed.")
    except:
        print(sys.exc_info()[0])
        print("Error with " + out_struct_name + ': ' + sys.exc_info()[0])

    #shutil.rmtree(wkdir) # this causes issues when program is run in parallel
    print("Time[s]:%s" % (time.time() - start_time))
    print("Date: %s" % time.asctime())
    print("Run id: %s" % run_id)
    sys.stdout = save_out
    sys.stderr = save_err
    log_file.close()
    return out_struct_name

if __name__ == "__main__":
    if len(sys.argv) < 5:
        raise SystemExit("syntax: %s <locations_path> <pdb_options:e.g 2wcd_ordered_df> <pdb_suff> <run_id>" % sys.argv[0])
    else:
        locations_path = sys.argv[1].strip()
        inputf = sys.argv[2].strip()
        pdb_suff = sys.argv[3].strip()
        run_id = sys.argv[4].strip()

    options, locations = initialize_repository(main_path=locations_path, instr_filename_path="EncoMPASS_options_relative_to__instruction_file_benchmark.txt")

    cesymm_order(locations, options, inputf, pdb_suff, mem_req_limit = 20000, run_id = run_id)



