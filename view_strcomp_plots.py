# Lucy Forrest 2023-06-29
# to select out chains that have a given number of TMs
# so can visualize their plots
# reads two arguments; nTMs pkl file and count of TMs to get list for

import os
import sys
import pickle as pkl

pkl_file: str = sys.argv[1]
tm_to_consider = int(sys.argv[2])

f = os.path.exists(pkl_file)
if not f: 
    print("File does not exist:", pkl_file)
    exit(1)
p = pkl.load(open(pkl_file, 'rb'))

print("Looking for chains with", tm_to_consider)

for PDB_identifier, chain_nTMs in p.items():
    for chain in chain_nTMs:

        if chain_nTMs[chain] == tm_to_consider:
            print("ds_", PDB_identifier, "_", chain, ".png")
