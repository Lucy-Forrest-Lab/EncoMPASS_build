# Name: select_pdbs_cesymm_from_symd.py
# Author: Antoniya A. Aleksandrova
# Language: Python 3
# Description: Find which fragments of proteins need to be further analysed with CE-Symm based on the results from SymD; both CE-Symm default and SymD runs must be complete

import os
import sys
sys.path.append(os.path.expandvars('$MAIN'))
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from initialize_repository import initialize_repository
from symmetry_exec_functions import *

def select_cesymm_from_symd(locations, pdb_suff):
    wkdir = locations['FSYSPATH']['symmtemp'] + "cesymm_from_symd/"
    if os.path.isdir(wkdir) == False:
        os.mkdir(wkdir)
    tm_archive = locations['FSYSPATH']['symmetry'] + '.tm_archive.pkl'
    a = pkl.load(open(tm_archive, 'rb'))
    output=open(locations['FSYSPATH']['symmtemp'] + "cesymm_from_symd/fragments_to_run.log","w")
    cesymm_from_symd_set = set()

    symd_completed, failed = check_completion(locations, sorted(list(a)), 'symd')
    cesymm_completed, failed = check_completion(locations, sorted(list(a)), 'cesymm')

    for pdb in sorted(list(a)):
        print(pdb)
        symd_out_dir = locations['FSYSPATH']['symd'] + pdb + "/"
        cesymm_out_dir = locations['FSYSPATH']['cesymm'] + pdb + "/"
        chains=';'.join(sorted(a[pdb]['tmchains'])) + ';'
        if os.path.isfile(locations['FSYSPATH']['whole'] + pdb + pdb_suff + '.pdb')==False:
            # change when STATUS flag becomes reliable
            print("{0} does not exist.".format(locations['FSYSPATH']['whole'] + pdb + pdb_suff + '.pdb'))
            continue
            #raise SystemExit("{0} does not exist.".format(locations['FSYSPATH']['whole'] + pdb + pdb_suff + '.pdb'))

        # Check whether the SymD/CE-Symm analysis on the pdb was completed as determined by check_completion.py
        if pdb not in symd_completed:
            print("{0}: SymD analysis incomplete".format(pdb))
            continue

        if pdb not in cesymm_completed:
            print("{0}: CE-Symm analysis incomplete".format(pdb))
            continue

        for ind in range(0,len(chains)//2):
            ch=chains[ind*2]
            if len(chains)==2:
                inputf=pdb
            else:
                inputf=pdb+"_"+ch
            print(inputf)
            # Get only TM chains (strip PDB)
            num=strip_tm_chains(wkdir,inputf,locations['FSYSPATH']['whole'] + pdb + pdb_suff + '.pdb',ch)
            os.rename(wkdir+inputf+"_tmp.pdb",wkdir+inputf+".pdb")
            oriented=wkdir+inputf+".pdb"

            # Load SymD data
            symmlist = {}
            symd_dic=symd_data(symd_out_dir,inputf)
            symmlist[inputf]={}
            symmlist[inputf]['symd']={}
            symmlist[inputf]['symd']=symd_dic

            #Load CE-Symm data
            cesymm_dic = cesymm_data(cesymm_out_dir,inputf)
            symmlist[inputf]['cesymm']={}
            symmlist[inputf]['cesymm']=cesymm_dic

            try:
                if symmlist[inputf]['cesymm']['repeats_number']!='na':
                    e=int(symmlist[inputf]['cesymm']['repeats_number'])
                e=int(symmlist[inputf]['symd']['symmetry_order'])
                e=float(symmlist[inputf]['symd']['coverage'])
            except ValueError:
                continue
            if symmlist[inputf]['cesymm']['repeats_number']!='na' and int(symmlist[inputf]['cesymm']['repeats_number'])<2 and symmlist[inputf]['symd']['topology']!='na' and int(symmlist[inputf]['symd']['symmetry_order'])<9 and float(symmlist[inputf]['symd']['coverage'])<0.5:
                # Find out whether the repeats are overlapping
                first_resid, last_resid, match_fr, match_lr, pdb_first, pdb_last=symd_sequence(symd_out_dir, inputf, oriented)
                x=sorted([int(first_resid),int(last_resid)])
                xs=set(range(x[0],x[1]+1))
                if int(match_fr)<=int(match_lr):
                    yr=range(int(match_fr),int(match_lr)+1)
                else:
                    yr=[i for j in (range(int(pdb_first),int(match_fr)+1), range(int(match_fr),int(pdb_last)+1)) for i in j]
                overlap=len(xs.intersection(yr))
                inputf=pdb+"_"+ch
                output.write(inputf+"_"+first_resid+"-"+last_resid+"\n")
                cesymm_from_symd_set.add(inputf+"_"+first_resid+"-"+last_resid)
                if overlap!=0 and (int(pdb_last)-int(last_resid)+1)>80:
                    fr=int(last_resid)+1
                    first_resid=str(fr)
                    last_resid=pdb_last
                    output.write(inputf+"_"+first_resid+"-"+last_resid+"\n")
                    cesymm_from_symd_set.add(inputf + "_" + first_resid + "-" + last_resid)

            os.remove(oriented)
    pkl.dump(cesymm_from_symd_set, open(locations['FSYSPATH']['cesymm_from_symd_log'] + 'fragments_to_run.pkl', 'wb'))
    return cesymm_from_symd_set

if __name__ == "__main__":
    locations_path = os.path.expandvars("ENC_DB")
    options, locations = initialize_repository(main_path=locations_path)
    pdb_suff = '_enc'
    cesymm_from_symd_set = select_cesymm_from_symd(locations, pdb_suff)