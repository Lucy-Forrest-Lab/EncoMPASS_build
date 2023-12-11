# quatsymm_complete.py
# Author: AAA
# Date: 6 Dec 2022
"""
Runs QuatSymm for a given structure and completes all the information necessary for EncoMPASS MSSD analysis.
"""

from symmetry_exec_functions import *
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from initialize_repository import initialize_repository
from quatsymm_default import *
from quatsymm_support import *
from Bio.SeqRecord import SeqRecord

start_time = time.time()

def quatsymm_complete(locations, options, inputf, pdb_suff, run_id="0000"):

	# Run QuatSymm
	quatsymm_default(locations, options, inputf, pdb_suff, run_id)

	log_file = open(locations['FSYSPATH']['quatsymm_log'] + inputf + '_qsymm.log', 'a')
	save_out = sys.stdout
	save_err = sys.stderr
	sys.stdout = sys.stderr = log_file

	print("Initiating EncoMPASS modifications (Run id: %s)." % run_id)
	pymol = "pymol -c"
	outdir = locations['FSYSPATH']['quatsymm'] + inputf + "/"
	wkdir = locations['FSYSPATH']['symmtemp'] + "quatsymm/" + inputf[0:4] + "/"

	qsymm_results =  read_quatsymm_out(outdir, inputf + "_qsymm")
	print("Quatsymm results: ", qsymm_results)
	if qsymm_results['symmetry'] == 'C1' or 'D' in qsymm_results['symmetry'] or qsymm_results['symmetry'] == 'na':
		print("No single-level symmetry found.")
		if os.path.isdir(wkdir):
			shutil.rmtree(wkdir)
	else:
		summary = {}
		summary['source'] = 'quatsymm'
		summary['inputf'] = inputf[0:4]
		summary['raw_data'] = qsymm_results
		summary['symmetry_order'] = qsymm_results['symmetry']
		summary['symmetry_levels'] = '1'  # we don't support multiple levels with QuatSymm
		summary['repeat_level'] = '1;'  # we don't support multiple levels with QuatSymm
		summary['internal_symmetry'] = 'No'
		summary['quaternary_symmetry'] = 'Yes'
		summary['symmetry_type'] = "Quaternary;"

		inputf = inputf + "_qsymm"
		qsymm=list(SeqIO.parse(outdir + inputf + ".fasta", "fasta"))
		if os.path.isdir(wkdir) == False:
			os.mkdir(wkdir)

		oriented = prep_struct(locations, wkdir, inputf[0:4], pdb_suff)
		summary['oriented'] = oriented

		reference=parse_structure(oriented)[0]
		resseq = get_pdb_sequence_with_chains(reference)  # creates a list of the form [(2, 'D', 'A'), (3, 'Q', 'A'), (4, 'A', 'A'), (5, 'L', 'A'), (6, 'S', 'A')...]

		seq_list = []
		for chain in reference:
			seq = ''.join([aa3to1.get(residue.resname, 'X') for residue in chain if is_aa(residue) and residue.has_id('CA')])
			seq_list.append(SeqRecord(seq=Seq(seq), id = chain.id, name = chain.id, description = chain.id))
		chains_str, chains_qsym = match_chain_names(reference)
		chain_order = [chains_qsym[x] for x in qsymm_results['subunits'].strip('[]').split(', ')]
		summary['chains'] = ";".join(chain_order) + ";"
		summary['raw_data']['chain_match_qsym_to_pdb'] = chains_qsym


		modified_qsym, mod_qsym_resid = quatsymm_chain_ranges(outdir + inputf + ".fasta", reference)
		print("\nChains matched:")
		for i in modified_qsym:
			print(i.id,i.name,i.description)
		reps, reps_resid, repeat_selection_flat, repeats, repeats_ax, repeats_transf, repeat_selection, unique_reps = chains_to_repeats(modified_qsym, mod_qsym_resid, symm_order= qsymm_results['symmetry'])
		print("repeats:", repeats)
		summary['repeat_selection'] = repeat_selection
		summary['repeat_selection_images'] = repeat_selection_flat
		summary['repeats'] = repeats.replace(";",",") + ";"
		summary['repeats_number'] = str(len(reps))
		repeat_len, aligned_res = modified_fasta(reps, inputf[0:4], outdir + inputf + pdb_suff + ".fasta")
		summary['alignment'] = ''
		summary['repeat_length'] = repeat_len
		summary['aligned_length'] = aligned_res
		quat_aln_map = create_aln_map(reps_resid)
		type = symm_type(qsymm_results['symmetry'])
		summary['closed_opened'] = type + ";"
		tmscore,rmsd,angle,translation,pt1,pt2,prot_len,axis_angle,axis_type, rot_matrix, transl_vector = cesymm_rmsd_tmscore(oriented, quat_aln_map, type)
		print("TM-score: %.2f, RMSD: %.2f, Angle: %.2f, Protein size: %d, Axis angle: %.2f, Axis type: %s\n"%(tmscore,rmsd,angle, prot_len,axis_angle,axis_type))
		summary['refined_tmscore'] = tmscore
		summary['refined_rmsd'] = str(rmsd)
		summary['unit_angle'] = str(angle) + ";"
		summary['unit_translation'] = str(translation) + ";"
		summary['size'] = prot_len
		summary['coverage'] = str(round(aligned_res/prot_len,2))
		summary['axis_angle_with_membrane_normal'] = str(axis_angle) + ";"
		summary['topology'] = axis_type



		axes_path=create_axes_file(outdir, inputf, angle, translation, pt1, pt2, repeats, len(quat_aln_map[0]), pdb_suff, type)
		# r.sort(key=lambda item: (len(item), item)) # sort r by length and alphabetically, i.e. A, B, C..., AA, BA, CA, .., BA,..

		summary['image_key'] = inputf + pdb_suff  # e.g. inputf + "_analysis"
		summary['files_key'] = inputf + pdb_suff
		summary['files_dir'] = outdir[len(locations['FSYSPATH']['main']):]
		summary['files_dir_path'] = outdir
		image_files, pml_files, jmol_files = cesymm_images(wkdir, outdir, summary['files_key'], repeat_selection_flat,
														   oriented, reference, 1,
														   pymol, locations['FSYSPATH']['quatsymm_pngs'], locations['FSYSPATH']['quatsymm_jsons'], summary['image_key'])

		summary['image_files'] = image_files
		summary['pml_files'] = pml_files
		summary['jmol_files'] = jmol_files

		super_pml_files = pymol_superposition_script(inputf + pdb_suff, rot_matrix, transl_vector, repeats_ax, repeat_selection_flat, unique_reps, locations['FSYSPATH']['quatsymm_super'], repeats_transf, inputf, level=0)
		summary['super_pml_files'] = super_pml_files
		summary['alignment'] = cesymm_alignment(summary['files_dir_path'], summary['files_key'])

		# create pkl; make sure files are saved in appropriate places
		print("\nFinal dictionary:\n ", summary)
		pkl.dump(summary, open(summary['files_dir_path'] + summary['files_key'] + ".pkl", 'wb'))
		shutil.rmtree(wkdir)

	# Clean-up
	print("{0} analysis completed.".format(inputf))
	print("Time[s]:%s" % (time.time() - start_time))
	print("Date: %s" % time.asctime())

	sys.stdout = save_out
	sys.stderr = save_err
	log_file.close()


if __name__ == "__main__":
	if len(sys.argv) < 5:
		raise SystemExit("syntax: %s <locations_path> <pdbcode> <pdb_suffix> <run_id>" % sys.argv[0])
	else:
		locations_path = sys.argv[1].strip()
		inputf = sys.argv[2][0:4]
		pdb_suff = sys.argv[3].strip()
		run_id = sys.argv[4].strip()

	options, locations = initialize_repository(main_path=locations_path, instr_filename_path="EncoMPASS_options_relative_to__instruction_file_benchmark.txt")
	if os.path.isdir(locations['FSYSPATH']['symmtemp'] + "quatsymm/") == False:
		os.mkdir(locations['FSYSPATH']['symmtemp'] + "quatsymm/")
	quatsymm_complete(locations, options, inputf, pdb_suff, run_id = run_id)