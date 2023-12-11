# Name: inferred_symmetry.py
# Author: Antoniya A. Aleksandrova
# Language: Python 3.8+
# Date: 30 Jan 2017 (previously version transfer_execution.py)
# Updated: May 2023
# Description:

import sys
import os
import pickle as pkl
import time
from symmetry_exec_functions import *
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from initialize_repository import initialize_repository

start_time = time.time()

def write_transfer_out_files(inferred_symm_dic, target_out_dir, target):
	with open(target_out_dir + target + "_transfer.out","w") as output:
		for k in ['template_order', 'coverage', 'size', 'unit_angle','unit_translation', 'refined_tmscore',
				  'refined_rmsd', 'repeats_number', 'topology', 'axis_angle_with_membrane_normal', 'symmetry_levels',
				  'aligned_length', 'repeat_length', 'template_keys']:
			value = '&'.join([inferred_symm_dic[i][k].strip(';') for i in range(len(inferred_symm_dic))])
			output.write(f"{k}\t{value.strip(';')}\n")
		output.write(f"template_structure\t{inferred_symm_dic[0]['template_structure']}\n")

	print(f"Wrote {target_out_dir + target + '_transfer.out'}.")
	parents = target_out_dir + target + "_transfer.parent"
	output=open(parents,"w")
	output.write(f"{'&'.join([inferred_symm_dic[i]['template_keys'] for i in range(len(inferred_symm_dic))])}")
	output.close()
	print(f"Wrote {parents}.")

def write_children_file(template_out_dir, template, target):
	children = template_out_dir + template[0:6] + "_transfer.children"
	if os.path.isfile(children):
		with open(children, "r") as f:
			prev_structures=f.readline().strip()
		output=open(children,"w")
		output.write(prev_structures+","+target)
		output.close()
	else:
		output=open(children,"w")
		output.write(target)
		output.close()
	print(f"Wrote {children}.")


def transfer_exec(template, target, locations, suffix, frtm_exec_path="/EncoMPASS/frtmalign/frtmalign", run_id="0000"):
	# Set options
	log_file = open(locations['FSYSPATH']['transfer_log'] + target + '_transfer.log', 'w')
	save_out = sys.stdout
	save_err = sys.stderr
	sys.stdout = sys.stderr = log_file
	print("Run id: %s" % run_id)

	print(f"Inferring symmetry from template {template} to target {target}.\n")
	symm_dic = pkl.load(open(locs['FSYSPATH']['mssd'] + "symmetry_results.pkl", "rb"))
	if len([k for k in symm_dic if template in k]) == 0:
		template_pdb = template[0:4]
	else:
		template_pdb = template

	if 'selected' not in symm_dic[template_pdb] or 'info' in symm_dic[template_pdb]['selected']:
		print("%s has no known symmetry.\n" % template)
		return
	sel = symm_dic[template_pdb]['selected']['source']
	transfer_out_dir = locations['FSYSPATH']['transfer'] + target + "/"
	if not os.path.isdir(transfer_out_dir):
		os.mkdir(transfer_out_dir)
	transfer_out_dir_template = locations['FSYSPATH']['transfer'] + template + "/"
	if not os.path.isdir(transfer_out_dir_template):
		os.mkdir(transfer_out_dir_template)
	wkdir = locations['FSYSPATH']['symmtemp'] + "transfer/" + target + "/"
	if not os.path.isdir(wkdir):
		os.mkdir(wkdir)
	super_pml_dir = locations['FSYSPATH']['transfer_chainssuper']

	if not os.path.isfile(locations['FSYSPATH']['whole'] + template[0:4] + suffix + '.pdb'):
		print("{0} does not exist.".format(locations['FSYSPATH']['whole'] + template[0:4] + suffix + '.pdb'))
		return
	oriented_opm = locations['FSYSPATH']['whole'] + template[0:4] + suffix + '.pdb'
	num = strip_tm_chains(wkdir, template, oriented_opm, template[5:6])
	os.rename(wkdir + template + "_tmp.pdb", wkdir + template + ".pdb")
	oriented_template = wkdir + template + ".pdb"

	if os.path.isfile(transfer_out_dir + target + "_transfer.out") or target == template: #skips targets that have already had their transfer
		print('Transfer already exists.')
		return
	if not os.path.isfile(locations['FSYSPATH']['whole'] + target[0:4] + suffix + '.pdb'):
		print("{0} does not exist.".format(locations['FSYSPATH']['whole'] + target[0:4] + suffix + '.pdb'))
		return
	oriented_opm = locations['FSYSPATH']['whole'] + target[0:4] + suffix + '.pdb'
	num = strip_tm_chains(wkdir, target, oriented_opm, target[5:6])
	os.rename(wkdir + target + "_tmp.pdb", wkdir + target + ".pdb")
	oriented_target = wkdir + target + ".pdb"

	t = time.localtime()
	frtmaln_file = locations['FSYSPATH']['strseqalns'] + 'strseqalns_' + template + ".txt"
	if os.path.isfile(frtmaln_file):
		with open(frtmaln_file, "r") as f:
			if f"INIT {template} {target}" not in f.read():
				print(f"Non-symmetrically aligned chains: {template} - {target}")
				frtmaln_file = locations['FSYSPATH']['strseqalns'] + 'strseqalns_' + target + ".txt"

	trans_dic = {target: []}
	for i, s in enumerate(sel.strip(";").split(";")):
		template_symm_data_dir = locations['FSYSPATH']['main'] + symm_dic[template_pdb][s]['files_dir']
		print("template files location: ", template_symm_data_dir)
		template_files_key = symm_dic[template_pdb][s]['files_key']
		template_symm_order = symm_dic[template_pdb][s]['symmetry_order']
		print("template symmetry order: ", template_symm_order)
		image_key = target + "_transfer"
		if len(sel.strip(";").split(";")) > 1:
			image_key += "_" + str(i + 1)
		trans_dic_tmp = transfer_symmetry(wkdir, frtm_exec_path, template_symm_data_dir, transfer_out_dir, oriented_template, oriented_target, template_files_key, target, template_symm_order, frtmaln_file, super_pml_dir, image_key)

		if trans_dic_tmp != {}:
			trans_dic_tmp['files_dir'] = locations['FSYS']['transfer'] + target + "/"
			trans_dic_tmp['template_keys'] = template_files_key
			trans_dic_tmp['template_structure'] = template
			trans_dic_tmp['date'] = f'{t.tm_year}-{t.tm_mon}-{t.tm_mday}'
			trans_dic_tmp['database'] = locations['FSYSPATH']['main'].split('/')[-2]
			trans_dic_tmp['axes_file'] = trans_dic_tmp['axes_file'].split(trans_dic_tmp['database'] +"/")[1]
			trans_dic_tmp['alnres_file'] = trans_dic_tmp['alnres_file'].split(trans_dic_tmp['database'] +"/")[1]
		if trans_dic_tmp=={}:
			print(f"Alignment doesn't include the symmetric regions found for {symm_dic[template_pdb][s]['files_key']}.")
			continue
		trans_dic[target].append(trans_dic_tmp)
		write_transfer_out_files(trans_dic[target], transfer_out_dir, target)
	write_children_file(transfer_out_dir_template, template, target)
	pkl.dump(trans_dic,open(transfer_out_dir + target + "_transfer.pkl", "wb"))
	print(trans_dic)

	# Clean-up
	print("{0} completed.".format(target))
	shutil.rmtree(wkdir)
	print("Time[s]:%s" % (time.time() - start_time))
	print("Date: %s" % time.asctime())
	sys.stdout = save_out
	sys.stderr = save_err
	log_file.close()
	return trans_dic


if __name__ == "__main__":
	if len(sys.argv) < 7:
		raise SystemExit("syntax: %s <locations_path> <template> <target> <pdb_suffix> <run_id> <FrTMAlign exec. path>" % sys.argv[0])
	else:
		locations_path = sys.argv[1].strip()
		template_chain = sys.argv[2]
		target_chain = sys.argv[3]
		pdb_suff = sys.argv[4].strip()
		run_id = sys.argv[5].strip()
		frtmpath = sys.argv[6]
	print('frtmpath ', frtmpath)
	options, locs = initialize_repository(main_path=locations_path, instr_filename_path="EncoMPASS_options_relative_to__instruction_file.txt")
	wkdir = locs['FSYSPATH']['symmtemp']+"transfer/"
	if not os.path.isdir(wkdir):
		os.mkdir(wkdir)
	transfer_results = transfer_exec(template_chain, target_chain, locs, pdb_suff, frtmpath)

