# Name: sort_symmetries.py
# Author: Antoniya A. Aleksandrova
# Language: Python 3.5+
# Date: 30 Jan 2017
# Updated: April 2023
# Description: Sort the selected symmetries in the Chain Symmetry Analysis section from best to worst 
#              Assume that, for example, D10>D2>C2>H>R and, 
#              if the order is the same, then the number of aligned residues should be considered, 
#              and if that's the same too, the total number of residues in the structure should be considered

import sys; print('Python %s on %s' % (sys.version, sys.platform))
import pickle as pkl
import os
import re
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from initialize_repository import initialize_repository


def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split('(\d+)', text) ]

def symm_order_for_mixed_sym(symm_dic):
    """
    Assign a symmetry order for cases with multiple independent symmetries

    :param symm_dic:
    :return:
    """
    sel = symm_dic['selected']['source']
    if ";" in sel:
        orders = [symm_dic[x]['symmetry_order'] for x in sel.strip(";").split(";")]
        orders.sort(key=natural_keys, reverse=True)
        if all(x == orders[0] for x in orders)==True:
            trump_order = orders[0]
        elif 'D' in orders:
            trump_order = [k for i,k in enumerate(orders) if 'D' in k][0]
        elif 'C' in orders:
            trump_order = [k for i,k in enumerate(orders) if 'C' in k][0]
        elif 'H' in orders:
            trump_order = [k for i,k in enumerate(orders) if 'H' in k][0]
        elif 'R' in orders:
            trump_order = [k for i,k in enumerate(orders) if 'R' in k][0]
        aligned_len = sum([int(symm_dic[x]['aligned_length']) for x in sel.strip(";").split(";")])
        size = [int(symm_dic[x]['size']) for x in sel.strip(";").split(";")][0]
        chain = [symm_dic[x]['chains'].strip(";") for x in sel.strip(";").split(";")][0]
    else:
        trump_order = symm_dic[sel]['symmetry_order']
        aligned_len = symm_dic[sel]['aligned_length']
        size = symm_dic[sel]['size']
        chain = symm_dic[sel]['chains'].strip(";")
    return trump_order, aligned_len, size, chain


def rank_by_symmetry(locations):
    results = pkl.load(open(locations['FSYSPATH']['mssd'] + "symmetry_results.pkl","rb"))
    #tm_archive = locations['FSYSPATH']['symmetry'] + '.tm_archive.pkl'
    symm_chains_list = []
    symm_complex_list = []
    for pdb in results.keys():
        if (len(pdb) == 4 and len([k for k in results if pdb in k]) == 1) or len(pdb)==6: #single-chain cases
            if 'selected' in results[pdb] and 'info' not in results[pdb]['selected']:
                order, alignment_len, prot_size, chain = symm_order_for_mixed_sym(results[pdb])
                # split order in a character and a digit for sorting
                if len(order) > 1:
                    digits = int(order[1:])
                else:
                    digits = ""
                symm_chains_list.append([pdb[0:4] + "_" + chain, order[0], digits, alignment_len, prot_size])
            elif 'selected' not in results[pdb]:
                print("%s was not run through MSSD, possibly because it is a beta-barrel or was too small for symmetry analysis."%(pdb))
        else:  # handle complexes
            if 'selected' in results[pdb] and 'quatinfo' not in results[pdb]['selected']:
                order, alignment_len, prot_size, chain = symm_order_for_mixed_sym(results[pdb])
                # split order in a character and a digit for sorting
                if len(order) > 1:
                    digits = int(order[1:])
                else:
                    digits = ""
                symm_complex_list.append([pdb, order[0], digits, alignment_len, prot_size])
            elif 'selected' not in results[pdb]:
                print("%s was not run through MSSD, possibly because it is a beta-barrel or was too small for symmetry analysis."%(pdb))

    def sorting(symm_list):
        pseudosorted_symm_list=[x for x in sorted(symm_list, key = lambda x : (x[1], x[2],x[3], x[4]), reverse=True)]
        custom_order = ['D', 'C', 'H', 'R']
        sorted_symm_list = []
        for c in custom_order:
            for l in pseudosorted_symm_list:
                if l[1] == c:
                    sorted_symm_list.append(l)
        return sorted_symm_list

    sorted_chains_symm_list = sorting(symm_chains_list)
    sorted_complex_symm_list = sorting(symm_complex_list)
#    print([x for x in sorted_symm_list if x[0] == '2a65_A'])
    pkl.dump(sorted_chains_symm_list, open(locations['FSYSPATH']['mssd'] + "ranked_chain_symmetries.pkl", "wb"))
    pkl.dump(sorted_complex_symm_list, open(locations['FSYSPATH']['mssd'] + "ranked_complex_symmetries.pkl", "wb"))

    return sorted_chains_symm_list, sorted_complex_symm_list

if __name__ == "__main__":
    locations_path = os.path.expandvars('$ENC_DB')
    options, locs = initialize_repository(main_path=locations_path, instr_filename_path="EncoMPASS_options_relative_to__instruction_file.txt")
    lst_chains, lst_complex = rank_by_symmetry(locs)

