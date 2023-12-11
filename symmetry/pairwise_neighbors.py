import sys
import os
import pickle as pkl
import time
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from initialize_repository import initialize_repository

start_time = time.time()

def extract_neighbors_for_struct(structname, fname, seqid_thresh = 0.85, tmscore_thresh = 0.6):
	"""
	Find all the neighbors of a given structure (they have to satisfy the requirements in both directions).
	:param structname: list of chains to consider, e.g. ['1okc_A','2a65_A','2a65_B',...]
	:param fname: name of the alignment results table
	:param seqid_thresh: MUSCLE sequence id threshold for neighbors
	:param tmscore_thresh: FrTMAlign TM-score threshold for neighbors
	:return: list of pairs of neighboring structures
	"""
	log_file = open(locations['FSYSPATH']['neighbors_log'] + structname + '_pairwise_neighbors.log', 'w')
	save_out = sys.stdout
	save_err = sys.stderr
	sys.stdout = sys.stderr = log_file

	print("Extracting neighbors for %s.\n"%(structname))
	neighbors = []
	with open(fname, 'r') as summary_table:
		for line in summary_table:
			if line.startswith("#"):
				continue
			pdbi1, ch1, pdbi2, ch2, alnlen, seq_seqid, str_seqid, rmsd, tmscore = line.split()
			if (structname == pdbi1 + "_" + ch1 or structname == pdbi2 + "_" + ch2) \
					and pdbi1 + "_" + ch1 != pdbi2 + "_" + ch2:
				if float(seq_seqid) > seqid_thresh or float(tmscore) > tmscore_thresh: # change this line if MUSCLE seq id is in different column
					neighbors.append([pdbi1 + "_" + ch1, pdbi2 + "_" + ch2])

	pairwise_neighbors = [x for x in neighbors if [x[1], x[0]] in neighbors]

	print("%s pairwise neighbor extraction completed.\n" % structname)
	print("Time[s]:%s" % (time.time() - start_time))
	print("Date: %s" % time.asctime())
	sys.stdout = save_out
	sys.stderr = save_err
	log_file.close()

	return pairwise_neighbors, neighbors


if __name__ == "__main__":
	if len(sys.argv) < 5:
		raise SystemExit("syntax: %s <locations_path> <pdb_options:e.g 2wcd_ordered_df> <pdb_suff> <run_id>" % sys.argv[0])
	else:
		locations_path = sys.argv[1].strip()
		inputf = sys.argv[2].strip()
		pdb_suff = sys.argv[3].strip()
		run_id = sys.argv[4].strip()

	options, locations = initialize_repository(main_path=locations_path, instr_filename_path="EncoMPASS_options_relative_to__instruction_file.txt")
	alignments_table = locations['SYSFILES']['summarytable']
	outdir = locations['FSYSPATH']['neighbors']
	outfile_pairwise = outdir + inputf + "_pairwise_neighbors.pkl"
	outfile_one_direction = outdir + inputf + "_one_direction_neighbors.pkl"
	pairwise_lst, lst = extract_neighbors_for_struct(inputf, alignments_table)
	pkl.dump(pairwise_lst, open(outfile_pairwise, "wb"))
	pkl.dump(lst, open(outfile_one_direction, "wb"))
