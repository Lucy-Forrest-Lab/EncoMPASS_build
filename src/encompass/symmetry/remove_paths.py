"""
Remove any reference to specific paths from publicly released files
"""
import subprocess
from encompass.pipeline.initialize_repository import initialize_repository
from encompass.symmetry.symmetry_exec_functions import *


def find_paths_to_change(fname, db_name):
	with open(fname, "r") as f:
		substitution_list = set({})
		for line in f:
			if db_name in line:
				print("From the previous database: ", fname)
				fields = line.split()
				revealing_paths = [x for x in fields if db_name in x]
				assert len(
					revealing_paths) < 2, f"Too many information-revealing paths in one line: file {fname}, line {line}"
				revealing_path = revealing_paths[0]
				if "pdb" in revealing_path:
					struct = revealing_path.split("/")[-1]
					substitution_list.add('/'.join(revealing_path.split("/")[:-1]))
				else:
					print("WARNING: uknown info-revealing path: ", revealing_path)
					substitution_list.add(revealing_path)
	return substitution_list


def replace_path(fname, substitution_list):
	print("substitution list" , substitution_list)
	for path in substitution_list:
		command = ["sed", "-i", "-e", 's|' + path + '|PATH|g', fname]
		print(command)
		subprocess.run(command, check=True)


def modify_file(fname, db_name):
	subs_paths = find_paths_to_change(fname, db_name)
	if len(subs_paths) > 0:
		replace_path(fname, subs_paths)
		print("Modified ", fname)
	return


if __name__ == "__main__":
	if len(sys.argv) < 2:
		raise SystemExit("syntax: %s <locations_path>" % sys.argv[0])
	else:
		locations_path = sys.argv[1].strip()
	options, locations = initialize_repository(main_path=locations_path,
											   instr_filename_path="EncoMPASS_options_relative_to__instruction_file_2023_12_13_shulien.txt")


	# database_name = locations['FSYSPATH']['main'].strip('/').split('/')[-1]
	database_name = "encompass_OPM20231213_PDB20231213"
	print(f"database name: {database_name}")
	# note that as of now, we're not changing the container-based path /EncoMPASS/symd1.61/src/symd1.61-linux
	# or it's leslie equivalent
	folders_to_change = ["cesymm", "cesymm_from_symd", "cesymm_low_thr", "cesymm_minlen", "cesymm_order",
	 					 "cesymm_quat_minlen", "cesymm_rmsd_2_5", "cesymm_rmsd_3", "quatsymm", "mssd", "symd"]
	for f in folders_to_change:
		print(f"folder: {f}")
		dir_list = [dir for dir in glob.glob(locations['FSYSPATH'][f] + "*") if os.path.isdir( dir) and ('jsons' not in dir and 'pngs' not in dir and 'super' not in dir)]
		for dir in dir_list:
			file_list = [file for file in glob.glob(dir + "/*") if '.pkl' not in file]
			for f in file_list:
				modify_file(f, database_name)

	#modify_file(file, database_name)
