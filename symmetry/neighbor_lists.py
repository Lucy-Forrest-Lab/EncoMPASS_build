"""

"""
import sys; print('Python %s on %s' % (sys.version, sys.platform))
import pickle as pkl
import os
import glob
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from initialize_repository import initialize_repository

def combine_pairwise_neighbor_lists(fdir, suffix = "pairwise_neighbors.pkl"):
	all_pairs = []
	for fname in glob.glob(fdir + "*" + suffix):
		struct_pairs_list = pkl.load(open(fname,"rb"))
		for p in struct_pairs_list:
			all_pairs.append(p)
	return all_pairs


def merge(lsts: list) -> list:
	"""
	Connects pairs of related structures into clusters of related structures
	e.g. [['1a0s_P', '1a0s_Q'], ['1a0s_P', '1a0s_R'], ['1a0s_P', '1a0t_P'], ...] into
	[frozenset({'1a0s_Q', '1a0s_Q', '1a0s_P', '1a0s_P', '1a0t_P',...}, frozenset({'1mpn_B'...})]

	"""

	sets = [set(lst) for lst in lsts if lst]
	merged = 1
	while merged:
		merged = 0
		results = []
		while sets:
			common, rest = sets[0], sets[1:]
			sets = []
			for x in rest:
				if x.isdisjoint(common):
					sets.append(x)
				else:
					merged = 1
					common |= x # common becomes the union of x and common
			results.append(frozenset([i for i in common]))
		sets = results
	print("Merged chain neighbors into sets.\n")
	return sets

def check_unity(structlist, neighbor_sets):
	for struct in structlist:
		clust_num = [1 for x in neighbor_sets if struct in x]
		if len(clust_num) == 0:
			print("WARNING: %s does not appear in any neighbor cluster!"%(struct))
		elif sum(clust_num) > 1 :
			print("WARNING: %s appears in multiple neighbor clusters. This will be a problem for transfer." %(struct))

def write_chain_neighbors(slist, chain_neigh, outname):
	# Collects neighbors for each structure
	# Writes it down on a file
	neighd = {}
	f = open(outname, 'w')
	f.write("#chain-name	neighbors\n")
	for struct in slist:
		sfound = []
		for s in chain_neigh:
			if struct in s:
				sfound = s
				break
		field1 = ''.join(sorted([x+';' for x in sfound if x != struct]))
		if field1:
			field1 = field1[:-1]
		else:
			field1 = 'None'
		f.write('{0}\t{1}\n'.format(struct, field1))
		neighd[struct] = [x for x in sfound if x != struct]
	f.close()
	# neighd = {'1a0s_P': ['1a0s_R', '1oh2_Q', '1oh2_R', '1a0t_R', '1oh2_P', '1a0s_Q', '1a0t_P', '1a0t_Q'], '1a0s_Q': ['1a0s_R', '1oh2_Q', '1oh2_R', '1a0t_R', '1oh2_P', '1a0t_P', '1a0s_P', '1a0t_Q'], ...]
	print("Neighboring chains file completed.\n")

def choose_chain_representative_based_on_symm(struct_list, symm_rank_list, neighbors_list, tm_arch = False):
	"""

	:param struct_list:
	:param symm_rank_list:
	:param neighbors_list:
	:param tm_arch:
	:return:
	"""

	sorted_neighbors = {} # include structures with no neighbors and no symmetry in the list
	sorted_neighbors_for_transfer = {} # excludes structures with no neighbors and no symmetry in the representative
	no_symm_rank = max(100000,len(symm_rank_list))
	for struct in struct_list:
		clust_num = [1 for x in neighbors_list if struct in x]
		if len(clust_num) == 0:
			neighbors_list.append(frozenset({struct}))

	for cluster in neighbors_list:
		clust_symm = []
		for struct in cluster:
			unk = 0
			if tm_arch:
				if tm_arch[struct[0:4]][struct[5:6]]['unk'] == True:
					unk = 1
			pair = [struct, no_symm_rank, unk]
			for rank, symm in enumerate(symm_rank_list):
				if symm[0] == struct:
					pair = [struct, rank, unk]
					break
			clust_symm.append(pair)

		clust_symm.sort(key=lambda x: (x[2], x[1], x[0]))  # sort by no-unique-resids, symmetry rank in ascending order and then alphabetically
		if clust_symm[0][1] != no_symm_rank and len(clust_symm) > 1:
			sorted_neighbors_for_transfer[clust_symm[0][0]] = [x[0] for x in clust_symm[1:]]
		sorted_neighbors[clust_symm[0][0]] = [x[0] for x in clust_symm[1:]]

	return sorted_neighbors, sorted_neighbors_for_transfer

def write_transfer_list(transfer_neighbors, outname):
	with open(outname, "w") as out:
		out.write("#parent\tchildren\n")
		for x in transfer_neighbors:
			out.write(f"{x}\t{';'.join(transfer_neighbors[x])}\n")



if __name__ == "__main__":
	locations_path=os.path.expandvars("$ENC_DB")  # path to the EncoMPASS database
	date_label = "4-10-2023"

	options, locs = initialize_repository(main_path=locations_path, instr_filename_path="EncoMPASS_options_relative_to__instruction_file.txt")
	pairwise_data_dir = locs['FSYSPATH']['neighbors']
	mode = "pairwise_neighbors"
	generate = True
	if generate == True:
		chain_pairs = combine_pairwise_neighbor_lists(pairwise_data_dir, suffix = mode + ".pkl")
		chain_clusters = merge(chain_pairs)
		pkl.dump(chain_clusters, open(pairwise_data_dir + "all_chain_clusters_" + mode + ".pkl", "wb"))
	else:
		chain_clusters = pkl.load(open(pairwise_data_dir + "all_chain_clusters_" + mode +".pkl", "rb"))
	tm_archive = pkl.load(open(locs['FSYSPATH']['symmetry'] + '.tm_archive.pkl', 'rb'))
	ch_list = [pdb + "_" + ch for pdb in tm_archive.keys() for ch in tm_archive[pdb]['tmchains']]
	write_chain_neighbors(ch_list, chain_clusters, locs['FSYSPATH']['neighbors'] + mode + "_4_each_chain_seqid0.85_tmscore_0.6.txt")
	ranked_symm = pkl.load(open(locs['FSYSPATH']['mssd'] + "ranked_chain_symmetries_" + date_label + ".pkl", "rb"))
	unique_chains, transfer_dic = choose_chain_representative_based_on_symm(ch_list, ranked_symm, chain_clusters, tm_archive)
	# unique_chains contains structures with no neighbors and no symmetry
	pkl.dump(transfer_dic, open(pairwise_data_dir + "parent_children_symmetry-ranked_neighbors.pkl", 'wb'))
	write_transfer_list(transfer_dic, pairwise_data_dir + "parent_children_symmetry-ranked_neighbors.txt")