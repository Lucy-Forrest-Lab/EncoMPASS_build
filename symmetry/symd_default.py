# Name: symd_default.py
# Author: Antoniya A. Aleksandrova
# Language: Python 3.5+
# Description: Executes a SymD 1.6 run with a single protein, taking into consideration guessed chain order

import os
import sys
import numpy as np
from symmetry_exec_functions import *
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from initialize_repository import initialize_repository
import time
start_time = time.time()


def symd_default(locations, options, inputf, pdb_suff, run_id="0000"):

    # Set options
    log_file = open(locations['FSYSPATH']['symd_log'] + inputf + '_symd.log', 'w')
    save_out = sys.stdout
    save_err = sys.stderr
    sys.stdout = sys.stderr = log_file
    
    print("Run id: %s" % run_id)
    print("In: "+inputf+" submitted")
    if not os.path.isfile(locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb'):
        raise SystemExit("{0} does not exist.".format(locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb'))
    wkdir=locations['FSYSPATH']['symmtemp']+"symd/"
    if not os.path.isdir(wkdir):
        os.mkdir(wkdir)
    if not os.path.isdir(wkdir + "temp_structures"):
        os.mkdir(wkdir + "temp_structures")
    chain_ordered_dir=locations['FSYSPATH']['symmtemp']+"chain_ordered_str/"
    if not os.path.isdir(chain_ordered_dir):
        os.mkdir(chain_ordered_dir)
    symd_exe=options['ALL']['sigsymd']
    print(symd_exe)
    out_dir=locations['FSYSPATH']['symd'] + inputf +"/"
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    oriented=locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb'
    tm_archive = locations['FSYSPATH']['symmetry'] + '.tm_archive.pkl'
    minlength=20 #20 residues
    # Get the membrane z-coordinates
    limit_bottom,limit_top=membrane_limits(oriented)

    # Find TM chains
    chain=""
    a = pkl.load(open(tm_archive, 'rb'))
    chain = ';'.join(sorted(a[inputf]['tmchains'])) + ';'  # tm chains as a string A;B;C;...

    if chain=="":
        raise SystemExit("{0} has no TM chains.".format(locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb'))
    else:
        print(chain)

    if len(chain)>2:
        monomer='N'
    else:
        monomer='Y'
    all_chain=int(len(chain)/2)  # the number of tm chains

    # Run on whole protein with TM chains only
    num=strip_tm_chains(wkdir,inputf,locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb',chain)
    os.rename(wkdir+inputf+"_tmp.pdb",wkdir+inputf+".pdb")
    oriented=wkdir+inputf+".pdb"
    if num==0:
        raise SystemExit("{0} has 0 TM atoms.".format(locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb'))

    res_num = resnum(oriented) #SymD 1.6w limits runs to 12000 residues
    if res_num > 12000:
        print("Skipping because there are more than 12000 residues in the structure ({0})".format(str(res_num)))
    else:
        symd(wkdir,symd_exe,out_dir,inputf,oriented,num)

    # Run on each TM chain of oligomer
    if monomer=='N':
        # Check the order of chains
        oriented_old=oriented
        if all_chain>3 and res_num < 12001:
            # Center-of-mass algorithm
            order=find_chain_order(inputf,oriented,limit_top,limit_bottom,chain)
            order_com_symd=""
            if order!="" and order!=chain and (order not in chain*2) and order!=chain[::-1]:
                print("In: CoM chain order: "+order)
                oriented=wkdir+inputf+"_ordered.pdb"
                symd(wkdir,symd_exe,out_dir,inputf+'_ordered',oriented,num)
                # Find the order from the SymD results on the center-of-mass (com) ordered chains
                order_com_symd=symd_chain_order(wkdir, out_dir, inputf+'_ordered', oriented, chain)
                if os.path.isfile(wkdir+inputf+'_ordered_sd_order.pdb'):
                  os.rename(wkdir+inputf+'_ordered'+'_sd_order.pdb',wkdir+inputf+'_com_sd.pdb')
                ##
                os.rename(wkdir+inputf+"_ordered.pdb", chain_ordered_dir+inputf+"_ordered.pdb")
                oriented=oriented_old
            if os.path.isfile(wkdir+inputf+"_ordered.pdb"):
                os.rename(wkdir+inputf+"_ordered.pdb",wkdir+"temp_structures/"+inputf+"_ordered.pdb")
            # Symd algorithm
            order_symd=symd_chain_order(wkdir, out_dir, inputf, oriented, chain)
            if order_symd!="" and order_symd!=order and order_symd!=chain and (order_symd not in order*2) and (order_symd not in chain*2) and (order_symd!=order[:-1][::-1]+';') and (order_symd!=chain[:-1][::-1]+';'):
                oriented=wkdir+inputf+"_sd_order.pdb"
                num_tmp=num_atoms_pdb(oriented)
                print("In: Symd chain order: "+order_symd)
                symd(wkdir,symd_exe,out_dir,inputf+"_sd_order",oriented,num_tmp)
                oriented=oriented_old
                os.rename(wkdir+inputf+"_sd_order.pdb",chain_ordered_dir+inputf+"_sd_order.pdb")
            if os.path.isfile(wkdir+inputf+"_sd_order.pdb"):
                os.rename(wkdir+inputf+"_sd_order.pdb",wkdir+"temp_structures/"+inputf+"_sd_order.pdb")
            # SymD algorithm on CoM results
            if order_com_symd!="" and order_com_symd!=order_symd and order_com_symd!=order and order_com_symd!=chain and (order_com_symd not in order*2) and (order_com_symd not in chain*2) and (order_com_symd not in order_symd*2) and (order_com_symd!=order_symd[:-1][::-1]+';') and (order_com_symd!=order[:-1][::-1]+';') and (order_com_symd!=chain[:-1][::-1]+';'):
                oriented=wkdir+inputf+"_com_sd.pdb"
                num_tmp=num_atoms_pdb(oriented)
                print("In: CoM and Symd chain order: "+order_com_symd)
                symd(wkdir,symd_exe,out_dir,inputf+"_com_sd",oriented,num_tmp)
                oriented=oriented_old
                os.rename(wkdir+inputf+"_com_sd.pdb",chain_ordered_dir+inputf+"_com_sd.pdb")
            if os.path.isfile(wkdir+inputf+"_com_sd.pdb"):
                os.rename(wkdir+inputf+"_com_sd.pdb",wkdir+"temp_structures/"+inputf+"_com_sd.pdb")
        # Get each chain
        for i in range(0,all_chain):
            ch_id=chain[i*2]
            print("In: Chain "+ch_id)
            oriented=oriented_old
            get_chain(oriented,ch_id,inputf,wkdir)
            oriented=wkdir+inputf+"_"+ch_id+".pdb"
            res_num=resnum(oriented)
            if res_num>minlength:
                num=num_atoms_pdb(oriented)
                symd(wkdir,symd_exe,out_dir,inputf+"_"+ch_id,oriented,num)
            else:
                print("Out: Chain {0} has too few residues ({1}).".format(inputf+"_"+ch_id,str(res_num)))
            os.rename(wkdir+inputf+"_"+ch_id+".pdb",wkdir+"temp_structures/"+inputf+"_"+ch_id+".pdb")

    # Clean-up
    print("{0} completed.".format(inputf))
    os.remove(wkdir + inputf + ".pdb")
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
        inputf = sys.argv[2][0:4]
        pdb_suff = sys.argv[3].strip()
        run_id = sys.argv[4].strip()

    currentdir = os.path.dirname(os.path.realpath(__file__))
    parentdir = os.path.dirname(currentdir)
    sys.path.append(parentdir)
    options, locations = initialize_repository(main_path=locations_path)
    if not os.path.isdir(locations['FSYSPATH']['symmtemp']+"symd/"):
        os.mkdir(locations['FSYSPATH']['symmtemp']+"symd/")
    if not os.path.isdir(locations['FSYSPATH']['symmtemp']+"symd/temp_structures"):
        os.mkdir(locations['FSYSPATH']['symmtemp']+"symd/temp_structures")

    symd_default(locations, options, inputf, pdb_suff, run_id = run_id)