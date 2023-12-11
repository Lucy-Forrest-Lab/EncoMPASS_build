# Author: Antoniya A. Aleksandrova


from symmetry_exec_functions import *
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from initialize_repository import initialize_repository

def transfer_symmetry_image(locations, template, target, pdb_suff, pymol=None, run_id="0000"):
	log_file = open(locations['FSYSPATH']['images_transfer_log'] + target + '_image.log', 'w')
	save_out = sys.stdout
	save_err = sys.stderr
	sys.stdout = sys.stderr = log_file
	print("Run id: %s" % run_id)
	print(f"Image generation for symmetry inferred from {template} to {target}.")
	ray = 1 # whether ray-tracing should be used in PyMOL image generation (0=no, 1=yes)
	if not pymol:
		pymol = "pymol -c"
	pic_dir_chain=locations['FSYSPATH']['transfer_chainspngs']
	if not os.path.isdir(pic_dir_chain + template + "/"):
		os.mkdir(pic_dir_chain + template + "/")
	jmol_dir_chain=locations['FSYSPATH']['transfer_chainsjsons']

	wkdir = locations['FSYSPATH']['symmtemp'] + "transfer/" + target + "/"
	if not os.path.isdir(wkdir):
		os.mkdir(wkdir)

	transfer_out_dir = locations['FSYSPATH']['transfer'] + target + "/"
	if not os.path.isfile(transfer_out_dir + target + "_transfer.pkl"):
		print(f"The inferred symmetry results for this target were not completed since {transfer_out_dir + target + '_transfer.pkl'} doesn't exist.")
		return {}
	transfer_dic = pkl.load(open(transfer_out_dir + target + "_transfer.pkl", "rb"))
	oriented_opm = locations['FSYSPATH']['whole'] + target[0:4] + pdb_suff + '.pdb'
	strip_tm_chains(wkdir, target, oriented_opm, target[5:6])
	os.rename(wkdir + target + "_tmp.pdb", wkdir + target + ".pdb")
	oriented_target = wkdir + target + ".pdb"
	reference = parse_structure(oriented_target)[0]
	for i,trans_dic_tmp in enumerate(transfer_dic[target]):
		image_key = target + "_transfer"
		if len(transfer_dic[target]) > 1:
			image_key += "_" + str(i + 1)
		image_files, pml_files, jmol_files = cesymm_images(wkdir,transfer_out_dir,trans_dic_tmp['files_key'],trans_dic_tmp['repeats_selection'],oriented_target,reference,ray,pymol,pic_dir_chain,jmol_dir_chain,image_key)
		print(image_files, pml_files)
		transfer_dic[target][i]['image_files'] = image_files
		transfer_dic[target][i]['pml_files'] = pml_files
		transfer_dic[target][i]['jmol_files'] = jmol_files
		for f in image_files.strip(";").split(";"):
			force_symlink(os.path.relpath(locations['FSYSPATH']['transfer_chainspngs'] + f, locations['FSYSPATH']['transfer_chainspngs'] + template + "/"), locations['FSYSPATH']['transfer_chainspngs'] + template + "/" + f)
			print('Creating link {} -> {}'.format(locations['FSYSPATH']['transfer_chainspngs'] + template + "/" + f, os.path.relpath(locations['FSYSPATH']['transfer_chainspngs'] + f, locations['FSYSPATH']['transfer_chainspngs'] + template + "/")))

	pkl.dump(transfer_dic, open(transfer_out_dir + target + "_transfer.pkl", "wb"))
	# Clean-up
	print("{0} completed.".format(target))
	shutil.rmtree(wkdir)
	print("Time[s]:%s" % (time.time() - start_time))
	print("Date: %s" % time.asctime())
	sys.stdout = save_out
	sys.stderr = save_err
	log_file.close()
	return transfer_dic

if __name__ == "__main__":
	if len(sys.argv) < 6:
		raise SystemExit("syntax: %s <locations_path> <template> <target> <pdb_suff> <run_id>" % sys.argv[0])
	else:
		locations_path = sys.argv[1].strip()
		template = sys.argv[2].strip()
		target = sys.argv[3].strip()
		pdb_suff = sys.argv[4].strip()
		run_id = sys.argv[5].strip()
		if len(sys.argv) > 6:
			pymol_path = sys.argv[6].strip()
		else:
			pymol_path = None

	options, locs = initialize_repository(main_path=locations_path, instr_filename_path="EncoMPASS_options_relative_to__instruction_file.txt")
	inferred_symm_dic = transfer_symmetry_image(locs, template, target, pdb_suff, pymol=pymol_path, run_id=run_id)
