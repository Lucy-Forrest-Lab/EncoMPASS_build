# Name: quatsymm_support.py
# Author: Antoniya A. Aleksandrova
# Language: Python 3.8+
"""
Post-processing QuatSymm results
"""

import time
import shutil
from symmetry_exec_functions import *
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from initialize_repository import initialize_repository
from typing import List, Tuple, Any


def read_quatsymm_out(quatsymm_out_dir, struct):

	# Name    Size    Subunits        Stoichiometry   Pseudostoichiometry     Symmetry        Local   Method  SymmRMSD        SymmTMscore
	# /data/local/.../quatsymm/2h8p/2h8p.pdb       8       [A, G, E, C, B, H, F, D]        A4B4    false   C4      false   ROTATION        0.00    1.00
	quatsymm_dic = {'size': 0, 'subunits': 'na', 'stoichiometry': 'na', 'pseudosymmetry': 'na', 'symmetry': 'na', 'local': 'na', 'method': 'na',
					'rmsd': 'na', 'tmscore': 'na'}
	if os.path.isfile(quatsymm_out_dir + struct + ".out") and os.path.getsize(quatsymm_out_dir + struct + ".fasta") > 5:
		result = open(quatsymm_out_dir + struct + ".out", "r")
		for ln in result:
			if struct[0:4] in ln:
				fields = ln.split('\t')
				if len(fields) == 10:
					quatsymm_dic = {'size': int(fields[1]), 'subunits': fields[2], 'stoichiometry': fields[3], 'pseudosymmetry': fields[4],
									'symmetry': fields[5], 'local': fields[6], 'method': fields[7], 'rmsd': float(fields[8]), 'tmscore': float(fields[9])}
				else:
					print("Warning: Irregularity in the file" + quatsymm_out_dir + struct + ".out. Number of columns is not 10.")
	else:
		print("No proper QuatSymm out file exists at ", quatsymm_out_dir + struct + ".out.")
	return quatsymm_dic

def prep_struct(locations, wkdir, struct, pdb_suff):
	"""

	:param locations:
	:param wkdir:
	:param struct:
	:param pdb_suff:
	:return:
	"""
	tm_archive = locations['FSYSPATH']['symmetry'] + '.tm_archive.pkl'
	a = pkl.load(open(tm_archive, 'rb'))
	chain = ';'.join(sorted(a[struct]['tmchains']))
	print(chain, locations['FSYSPATH']['whole']+struct+pdb_suff+'.pdb')
	num=strip_tm_chains(wkdir,struct,locations['FSYSPATH']['whole']+struct+pdb_suff+'.pdb',chain)
	os.rename(wkdir+struct+"_tmp.pdb", wkdir+struct+".pdb")
	oriented=wkdir+struct+".pdb"
	return oriented

def read_quatsymm_fasta(quatsymm_out_dir, struct, oriented):

	cesm=list(SeqIO.parse(quatsymm_out_dir + struct +".fasta", "fasta"))
	reference=parse_structure(oriented)[0]
	resseq_A = get_pdb_sequence_with_chains(reference) # creates a list of the form [(2, 'D', 'A'), (3, 'Q', 'A'), (4, 'A', 'A'), (5, 'L', 'A'), (6, 'S', 'A')...]

def match_chain_names(structure: Any) -> Tuple[dict, dict]:
	# 1. Get chain names in actual structure
	struct_chains = structure.get_chains()
	# 2. Guess chain names based on order and alphabet
	alphabet = letters_range()
	ch_dic_struct = {}
	ch_dict_qsym = {}
	for i,x in enumerate(struct_chains):
		ch_dic_struct[x.id] = alphabet[i]  # the correspondence between structure chains (keys) and QuatSymm chain names (values)
		ch_dict_qsym[alphabet[i]] = x.id  # the correspondence between QuatSymm chain names (keys) and structure chains (values)
	return ch_dic_struct, ch_dict_qsym


def quatsymm_chain_ranges(qsymm_fasta: str, structure: Any) -> Tuple[List, List]:
	"""
	Load QuatSymm fasta file and modify the record to include the pdb-relevant chain names and
	residue ranges corresponding to each alignment

	:param qsymm_fasta: the filename of the QuatSymm-generated repeat alignment
	:param structure: a BioPython loaded pdb (model level)
	:return: qsym: a modified list of the aligned QuatSymm repeats with id corresponding to the QuatSymm chain name, name corresponding to the
			pdb file chain name, and description containing the residue range of the sequence,
			e.g. [SeqRecord(seq=Seq('SALHWRAAGAATVLLVIVLLAGSYLAVLAERGAPGAQLITYPRALWWACETATTVGY'), id='_A', name='C', description='22-78', dbxrefs=[]),...]
			qsym_resid: the corresponding residue list of the repeats, e.g. [[(22,'S','C','',''), (21,'A','C','',''),...], ...]
	"""
	qsym = list(SeqIO.parse(qsymm_fasta, "fasta"))
	struct_chains, qsym_chains = match_chain_names(structure)

	# Ensure that the CE-Symm alignment sequence is part of the guessed chain & get residue ranges
	chains_seq = get_chains_seq(structure) # What happens when there is no chain name?
	chains_seq_resids = get_pdb_sequence_with_chains(structure)

	qsym_resid = [[] for i in qsym]
	for i,entry in enumerate(qsym):
		qsym[i].name = qsym_chains[entry.id.strip('_')]
		for ch in chains_seq:
			if ch.id == qsym[i].name:
				qsym_seq_clean = ''.join(c.upper() for c in entry.seq if c not in '/-')
				for n in range(3,len(qsym_seq_clean)):
					flag = ch.seq.find(qsym_seq_clean[0:n])
					if flag == -1:
						print("Difference! The QuatSymm sequence has only matched the pdb sequence up to the last residue here:", qsym_seq_clean[0:n])
						break

				ind_first = ch.seq.find(qsym_seq_clean)
				if ind_first == -1:
					raise SystemExit("Error: QuatSymm repeat sequence %s is not part of chain %s for %s."%(entry.id, ch.id, qsymm_fasta))
				ind_last = ind_first + len(qsym_seq_clean) - 1
				resid_count = 0
				for r in chains_seq_resids:
					if r[2] == ch.id:
						if resid_count == ind_first:
							resid_first = r[0]
							aln_count = 0
						elif resid_count == ind_last:
							resid_last = r[0]
						if resid_count >= ind_first and resid_count <= ind_last:
							while entry.seq[aln_count] == '-':
								qsym_resid[i].append(('-','-','-','-','-'))
								aln_count +=1

							assert entry.seq[aln_count].upper() == r[1], print("Mismatch between chain sequence and QuatSymm alignment sequence")
							qsym_resid[i].append((r[0], entry.seq[aln_count], r[2], r[3], r[4]))
							aln_count +=1
						resid_count +=1
				qsym[i].description = str(resid_first) + "-" + str(resid_last)


	return qsym, qsym_resid

def chains_to_repeats(qsym_aln, qsym_resid, symm_order):
	"""

	:param qsym_aln: a list of the aligned QuatSymm repeats generated by the program in a BioPython SeqRecord form,
			e.g. [SeqRecord(seq=Seq('SALHWRAAGAATVLLVIVLLAGSYLAVLAERGAPGAQLITYPRALWWACETATTVGY'), id='_A', name='C', description='22-78', dbxrefs=[]),...]

	:param qsym_resid:

	:param symm_order: the symmetry order; used for calculating the number of expected repeats

	:return: 1) a re-ordered list of the alignment where each element of the list corresponds to one repeat, e.g.
	[[SeqRecord(seq=Seq('SALHWRAAGAA'), id='_A', name='C', description='22-78', dbxrefs=[]), SeqRecord(seq=Seq('DLYPVTLWGRL'), id='_B', name='D', description='80-122', dbxrefs=[])],
	[SeqRecord(seq=Seq('SALHWRAAGAA'), id='_G', name='O', description='22-78', dbxrefs=[]), SeqRecord(seq=Seq('DLYPVTL'), id='_H', name='P', description='80-122', dbxrefs=[])],
	[SeqRecord(seq=Seq('SALHWRAAGAA'), id='_E', name='K', description='22-78', dbxrefs=[]), SeqRecord(seq=Seq('DLYPVTL'), id='_F', name='L', description='80-122', dbxrefs=[])]],
	and 2) a similar list but with the residue ids instead:
	[[[(22, 'S', 'C', ' '), (23, 'A', 'C', ' '), ...],[(80, 'D', 'D', ' '), (81, 'L', 'D', ' '),...]],
	[[(22, 'S', 'O', ' '), (23, 'A', 'O', ' '), ...],[(80, 'D', 'P', ' '), (81, 'L', 'P', ' '),...]],[[...]]]
	3) a similar list but in a form compatible with PyMOL selections: [[(22,'C'),(78,'C'),(80,'D'),(122,'D')],[(22,'O'),(78,'O'),(80,'P'),(122,'P)],[..]]
	4) "repeats": an axes-file compatible format, e.g. (C_22-78;O_27-78;K_27-78)(D_80-122;P_80-122;L_80-122)
	"""
	try:
		dim = int(''.join([s for s in symm_order if s.isdigit()]))
	except:
		"Number of repeats will be estimated based on // in sequence alignment rather than order."
		dim = len(qsym_aln)//len([1 for i in qsym_aln if '//' in i.seq])
	reps = [[] for _ in range(dim)]
	reps_sel = [[] for _ in range(len(qsym_aln))]
	reps_sel_layers = [[] for _ in range(dim)]
	reps_resid = [[] for _ in range(dim)]
	reps_ax = [[] for _ in range(len(qsym_aln)//dim)]
	repeats_transf=[[-1 for j in range(0, len(qsym_aln))] for i in range(0,1)]  # same number of elements as symmetry levels, eg. [[-1, -1, -1, -1], [-1, -1, -1, -1]]
	unique_repeats = []

	k=0
	m=0
	for i,r in enumerate(qsym_aln):
		reps[k].append(r)
		first_digit=r.description[0]  # in case the residues has a negative id (i.e. -1)
		beg=r.description[1:].split('-')[0]
		beg=first_digit+beg
		end=r.description[1:].split('-')[-1]
		ind = m + k*(len(qsym_aln)//dim)
		reps_sel[ind].append((int(beg), r.name))
		reps_sel[ind].append((int(end), r.name))
		reps_sel_layers[k].append((int(beg), r.name))
		reps_sel_layers[k].append((int(end), r.name))
		unique_repeats.append(r.name + '_' + r.description)
		reps_resid[k].append(qsym_resid[i])
		reps_ax[m].append(r.name + "_" + r.description)
		repeats_transf[0][i] = k
		k+=1
		if k == dim:
			k = 0
			m+=1


	repeats = ""
	for i, r in enumerate(reps_ax):
		repeats = repeats + "(" + ";".join(r) + ")"
	repeats = repeats + ""

	return reps, reps_resid, reps_sel, repeats, reps_ax, repeats_transf, reps_sel_layers, unique_repeats

def modified_fasta(reps, pdb, fname):
	"""
	Write a EncoMPASS-compliant FASTA file of the repeats from QuatSymm
	:param reps: the aligned sequences with each element of the list corresponds to one repeat, e.g.
	[[SeqRecord(seq=Seq('SALHWRAAGAA'), id='_A', name='C', description='22-78'), SeqRecord(seq=Seq('DLYPVTLWGRL'), id='_B', name='D', description='80-122')],
	[SeqRecord(seq=Seq('SALHWRAAGAA'), id='_G', name='O', description='22-78'), SeqRecord(seq=Seq('DLYPVTL'), id='_H', name='P', description='80-122')],
	[SeqRecord(seq=Seq('SALHWRAAGAA'), id='_E', name='K', description='22-78'), SeqRecord(seq=Seq('DLYPVTL'), id='_F', name='L', description='80-122')]]

	:param pdb: the PDB id to be included in the sequence tag in the FASTA file
	:param fname: the name of the FASTA file that should be produced
	:return:
	"""
	f = open(fname, "w")
	aligned = 0
	for r in reps:
		name = [c.name + "_" + c.description for c in r]
		f.write(">" + pdb.upper() + "." + ";".join(name) + "\n")
		seq = [str(c.seq.strip("/")) for c in r]
		repeat = [s for s in ''.join(seq) if s != '-' and s!= '/']
		aligned += len(repeat)
		rep_len = len("".join(seq))
		f.write("/".join(seq) + "\n")
	f.close()
	return rep_len, aligned

def get_chains_seq(structure):
	"""
	Retrieves the AA sequence from a PDB structure for each chain in a BioPython SeqRecord format.
	It's a list that looks like [SeqRecord(seq=Seq('SALHWRAAGAATVLLVIVLLAGSYLAVLAERGAPGAQLITYPRALWWACETATTVGY'), id='A', name='A', description='A'), SeqRecord(seq=Seq('SALHWRAAGAATVLLVIVLLAGSYLAVLAERGAPGAQLITYPRALWWACETATTVGY'), id='B', name='B', description='B')]

	:param structure: a model from a BioPython parsed pdb file
	:return: seq_list, a list of chains and their sequences in a SeqRecord format
	"""

	from Bio.SeqRecord import SeqRecord
	seq_list = []
	for chain in structure:
		seq = ''.join([aa3to1.get(residue.resname, 'X') for residue in chain if is_aa(residue) and residue.has_id('CA')])
		seq_list.append(SeqRecord(seq=Seq(seq), id = chain.id, name = chain.id, description = chain.id))
	return seq_list

def letters_range():
	"""
	Creates a list of all possible 1- and 2-letter combinations ordered alphabetically, i.e. ['A', 'B',..., 'AA', 'BA',..]
	:return: ['A', 'B',..., 'AA', 'BA',..]
	"""

	import string
	alphabet = list(string.ascii_uppercase)
	r = alphabet
	for l in alphabet:
		r = r + [c+l for c in alphabet]
	return r

#if dihedral, this method won't work; check that symmetry is not dihedral; e.g. of dihedral in PDB: 2wcd
def create_aln_map(reps_resid): #[[[A],[B]],[[C],[D]],[[F],[G]]] -> [[A,C,F],[B,D,G]]
	"""
	:param reps_resid: [[[A],[B]],[[C],[D]],[[F],[G]]]
	:return: aln_map:  [[A,C,F],[B,D,G]], e.g. [[(3, 'N', ' ', ' '), (3, 'E', ' ', ' ')], [(4, 'N', ' ', ' '), (4, 'E', ' ', ' ')]...]; note that residue names are not included
	"""
	pairs = [item for sublist in reps_resid[0] for item in sublist] # flatten reps_resid[0]
	aln_map = [[] for i in pairs] # [[],[]]
	for num, r in enumerate(reps_resid): # num=0, r=[[A],[B]]
		repeat = [item for sublist in r for item in sublist] # [A, B]
		map_count = 0
		#for map_count, repeat in enumerate(r):
		for c in repeat:
			if c[1] == '-' or c[1].islower():
				aln_map[map_count].append(('-','-','-','-'))
			else:
				aln_map[map_count].append((c[0], c[2], c[3], c[4]))
			map_count +=1
	return aln_map

# Test with 7phq for gaps

def symm_type(order):
	"""
	Guess whether symmetry is closed and order based on the letter in the symmetry order
	:param order: e.g. C22 or H2
	:return: "CLOSED" or "OPEN"
	"""
	order_letter = ''.join([s for s in order if not s.isdigit()])
	if order_letter == "C" or order_letter == "D":
		type = "CLOSED"
	else:
		type = "OPEN"

	return type



if __name__ == "__main__":
	if len(sys.argv) < 3:
		raise SystemExit("syntax: %s <QuatSymm results directory> <structure name>" % sys.argv[0])
	else:
		outdir = sys.argv[1].strip()
		inputf = sys.argv[2].strip()
	qsymm_results =  read_quatsymm_out(outdir, inputf)
	print("Qsymm results: ", qsymm_results)
	if qsymm_results['symmetry'] == 'C1' or 'D' in qsymm_results['symmetry'] or qsymm_results['symmetry'] == 'na':
		print("No single-level symmetry found.")
	else:
		summary = {}
		summary['source'] = 'quatsymm'
		summary['inputf'] = inputf[0:4]
		summary['raw_data'] = qsymm_results
		summary['chains'] = qsymm_results['subunits'].strip('[]').replace(', ', ';') + ";"
		summary['symmetry_order'] = qsymm_results['symmetry']
		summary['symmetry_levels'] = '1'  # we don't support multiple levels with QuatSymm
		summary['repeat_level'] = '1;'  # we don't support multiple levels with QuatSymm
		summary['internal_symmetry'] = 'No'
		summary['quaternary_symmetry'] = 'Yes'
		summary['symmetry_type'] = "Quaternary;"

		qsymm=list(SeqIO.parse(outdir + inputf + ".fasta", "fasta"))

		locations_path = os.path.expandvars("$ENC_DB")  # path to the EncoMPASS
		options, locations = initialize_repository(main_path=locations_path, instr_filename_path="EncoMPASS_options_relative_to__instruction_file.txt")
		wkdir = locations['FSYSPATH']['symmtemp'] + "quatsymm/"
		if os.path.isdir(wkdir) == False:
			os.mkdir(wkdir)
		oriented = prep_struct(locations, wkdir, inputf[0:4], "_enc")
		summary['oriented'] = oriented

		reference=parse_structure(oriented)[0]
		resseq = get_pdb_sequence_with_chains(reference) # creates a list of the form [(2, 'D', 'A'), (3, 'Q', 'A'), (4, 'A', 'A'), (5, 'L', 'A'), (6, 'S', 'A')...]
		from Bio.SeqRecord import SeqRecord
		seq_list = []
		for chain in reference:
			seq = ''.join([aa3to1.get(residue.resname, 'X') for residue in chain if is_aa(residue) and residue.has_id('CA')])
			seq_list.append(SeqRecord(seq=Seq(seq), id = chain.id, name = chain.id, description = chain.id))
		chains_seq = match_chain_names(reference)


		modified_qsym, mod_qsym_resid = quatsymm_chain_ranges(outdir + inputf + ".fasta", reference)
		print("Chains matched:\n")
		for i in modified_qsym:
			print(i.id,i.name,i.description)
		reps, reps_resid, repeat_selection_flat, repeats, repeats_ax, repeats_transf, repeat_selection, unique_reps = chains_to_repeats(modified_qsym, mod_qsym_resid, symm_order= qsymm_results['symmetry'])
		print("repeats:", repeats)
		summary['repeat_selection'] = repeat_selection
		summary['repeats'] = repeats.replace(";",",") + ";"
		summary['repeats_num'] = str(len(reps))
		repeat_len, aligned_res = modified_fasta(reps, inputf[0:4], outdir + inputf + "_enc.fasta")
		summary['alignment'] = ''
		summary['repeat_length'] = repeat_len
		summary['aligned_length'] = aligned_res
		quat_aln_map = create_aln_map(reps_resid)
		type = symm_type(qsymm_results['symmetry'])
		summary['closed_opened'] = type + ";"
		tmscore,rmsd,angle,translation,pt1,pt2,prot_len,axis_angle,axis_type, rot_matrix, transl_vector = cesymm_rmsd_tmscore(oriented, quat_aln_map, type)
		print("TM-score: %f, RMSD: %f, Angle: %.2f, Protein size: %d, Axis angle: %.2f, Axis type: %s\n"%(tmscore,rmsd,angle, prot_len,axis_angle,axis_type))
		summary['refined_tmscore'] = tmscore
		summary['refined_rmsd'] = str(rmsd)
		summary['unit_angle'] = str(angle) + ";"
		summary['unit_translation'] = str(translation) + ";"
		summary['size'] = prot_len
		summary['coverage'] = str(round(aligned_res/prot_len,2))
		summary['axis_angle_with_membrane_normal'] = str(axis_angle) + ";"
		summary['topology'] = axis_type



		axes_path=create_axes_file(outdir, inputf, angle, translation, pt1, pt2, repeats, len(quat_aln_map[0]), "_enc", type)

		print("Repeat selection: ", repeat_selection)
		pymol = "/data/TMB-CSB/apps/CentOS7-LabLinux/pymol/2.3.4/bin/pymol -c"
		summary['image_key'] = inputf #e.g. inputf + "_analysis"
		summary['files_key'] = inputf + "_enc"
		summary['files_dir'] = locations['FSYS']['quatsymm'] + inputf[0:4] + "/"
		summary['files_dir_path'] = outdir
		image_files, pml_files, jmol_files = cesymm_images(wkdir, outdir, summary['files_key'], repeat_selection_flat,
														   oriented, reference, 1,
														   pymol, outdir, outdir, summary['image_key'])

		summary['image_files'] = image_files
		summary['pml_files'] = pml_files
		summary['jmol_files'] = jmol_files

		super_pml_files = pymol_superposition_script(inputf + "_enc", rot_matrix, transl_vector, repeats_ax, repeat_selection, unique_reps, outdir, repeats_transf, inputf, level=0)
		summary['super_pml_files'] = super_pml_files
		summary['alignment'] = cesymm_alignment(summary['files_dir_path'], summary['files_key'])


		# create pkl; make sure files are saved in appropriate places
		print("Final dictionary: ", summary)
		pkl.dump(summary, open(summary['files_dir_path'] + summary['files_key'] + ".pkl", 'wb'))