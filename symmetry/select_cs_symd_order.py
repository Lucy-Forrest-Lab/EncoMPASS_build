# Name: select_cs_symd_order.py
# Author: Antoniya A. Aleksandrova
# Language: Python 3.5+
# Description: Compares the symmetry order of SymD and CE-Symm to figure out if CE-Symm should be run with a different order (eg: 4j7c_I)
# Prerequisites: Chain (internal) SymD and CE-Symm default results


import sys
import os
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from initialize_repository import initialize_repository
from symmetry_exec_functions import *

# Locate key directories and files
def select_cs_symd_order(locations):
    tm_archive = locations['FSYSPATH']['symmetry'] + '.tm_archive.pkl'
    out=open(locations['FSYSPATH']['cesymm_order_log'] + 'cs_symd_order.log','w')
    out.write("#Input\tCE-Symm Seed\tSymd Order\tCE-Symm Symmetry Levels\tCE-Symm Repeat Length\n")
    a = pkl.load(open(tm_archive, 'rb'))
    cs_symd_order_dic = {}

    for pdb in sorted(list(a)):
        symd_out_dir = locations['FSYSPATH']['symd'] + pdb + "/"
        cesymm_out_dir = locations['FSYSPATH']['cesymm'] + pdb + "/"
        tm_chains_list = sorted(a[pdb]['tmchains'])
        prot_type = a[pdb]['class']
        mode = '_df'
        if prot_type=='alpha':
            inputf = pdb
            all_names = ['', '_ordered', '_com_sd', '_sd_order']   # Handle chain re-ordering
            for name in all_names:
                if os.path.isfile(cesymm_out_dir + inputf + name + "_stdout.out"):
                    cesymm = cesymm_data(cesymm_out_dir, inputf + name)
                    symd = symd_data(symd_out_dir, inputf + name)
                    if symd['symmetry_order']!='na' and cesymm['repeats_number']!='na' and int(cesymm['repeats_number'])>1 and 10>int(symd['symmetry_order'])>3 and float(symd['coverage'])>0.7 and int(cesymm['symmetry_levels'])==1 and int(symd['symmetry_order']) != int(cesymm['repeats_number'])  and int(symd['symmetry_order']) % int(cesymm['repeats_number'])==0 and int(cesymm['repeat_length'])/(int(symd['symmetry_order']) / int(cesymm['repeats_number']))>40:
                       print(inputf+name, 'increase to ',symd['symmetry_order'])
                       out.write(inputf + name + mode + "\t" + cesymm['seed'] + "\t" + str(symd['symmetry_order'])+ "\t" + str(cesymm['symmetry_levels']) + "\t" + str(cesymm['repeat_length']) + "\n")
                       cs_symd_order_dic[inputf + name + mode] = {'name': inputf + name + mode,'seed': cesymm['seed'], 'symd_order': str(symd['symmetry_order']), 'cesymm_symm_levels': str(cesymm['symmetry_levels']), 'cesymm_rep_length': str(cesymm['repeat_length']), 'nres': a[pdb]['nres'], 'natoms': a[pdb]['natoms']}
                    if symd['symmetry_order']!='na' and cesymm['repeats_number']!='na' and int(cesymm['repeats_number'])>1 and 20>int(symd['symmetry_order'])>3 and float(symd['coverage'])>0.8 and int(cesymm['symmetry_levels'])>1 and (int(symd['symmetry_order']) == int(cesymm['repeats_number'])  or int(symd['symmetry_order']) % int(cesymm['repeats_number'])==0 or int(cesymm['repeats_number']) % int(symd['symmetry_order']) == 0) :
                       print(inputf+name, 'forced to ',symd['symmetry_order'])
                       out.write(inputf + name + mode + "\t" + cesymm['seed'] + "\t" + str(symd['symmetry_order']) + "\t" + str(cesymm['symmetry_levels']) + "\t" + str(cesymm['repeat_length']) + "\n")
                       cs_symd_order_dic[inputf + name + mode] = {'name': inputf + name + mode, 'seed': cesymm['seed'], 'symd_order': str(symd['symmetry_order']), 'cesymm_symm_levels': str(cesymm['symmetry_levels']), 'cesymm_rep_length': str(cesymm['repeat_length']), 'nres': a[pdb]['nres'], 'natoms': a[pdb]['natoms']}


            for ch in tm_chains_list:
                inputf=pdb+'_'+ch
                if os.path.isfile(cesymm_out_dir + inputf + ".axes"):
                   cesymm = cesymm_data(cesymm_out_dir, inputf)
                   symd = symd_data(symd_out_dir,inputf)
                   if symd['symmetry_order']!='na' and cesymm['repeats_number']!='na' and int(cesymm['repeats_number'])>1 and 10>int(symd['symmetry_order'])>3 and float(symd['coverage'])>0.4 and int(cesymm['symmetry_levels'])==1 and int(symd['symmetry_order']) != int(cesymm['repeats_number'])  and int(symd['symmetry_order']) % int(cesymm['repeats_number'])==0:
                      print(inputf, symd['symmetry_order'])
                      out.write(inputf + name + mode + "\t" + cesymm['seed'] + "\t" + str(symd['symmetry_order'])+ "\t" + str(cesymm['symmetry_levels']) + "\t" + str(cesymm['repeat_length']) + "\n")
                      cs_symd_order_dic[inputf + name + mode] = {'name': inputf + name + mode, 'seed': cesymm['seed'],
                                                      'symd_order': str(symd['symmetry_order']),
                                                      'cesymm_symm_levels': str(cesymm['symmetry_levels']),
                                                      'cesymm_rep_length': str(cesymm['repeat_length']),
                                                      'nres': a[pdb][ch]['nres'],
                                                      'natoms': a[pdb][ch]['natoms']}
    pkl.dump(cs_symd_order_dic, open(locations['FSYSPATH']['cesymm_order_log'] + 'cs_symd_order.pkl','wb'))
    return cs_symd_order_dic


if __name__ == "__main__":
    locations_path = os.path.expandvars("ENC_DB")
    options, locations = initialize_repository(main_path=locations_path)
    cs_symd_order_dic = select_cs_symd_order(locations)