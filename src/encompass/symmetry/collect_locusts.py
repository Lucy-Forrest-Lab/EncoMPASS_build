from encompass.symmetry.run_locusts import *

def collect_results(locations_path, locations, job_code, stepname, local_dir):
	fn_opt = default_fn_options(options)
	pdb_list = ['1okc', '2a65']
	command_template = 'source ~/encompass/bin/activate; python cesymm_default.py <arg0> <arg1> <arg2> <arg3> <arg4>'
	run_locusts(fn_opt, locations_path, locations, job_code, stepname, local_dir, pdb_list, command_template, opt_arg = '', gather = True)
	return

if __name__ == "__main__":
	from encompass.pipeline.initialize_repository import initialize_repository
	locations_path = "/data/local/encompass/dev/toni_full_dec2020_beforeFrTM"
	options, locations = initialize_repository(main_path=locations_path, instr_filename_path="EncoMPASS_options_relative_to__instruction_file.txt")
	#command = 'source ~/encompass/bin/activate; python cesymm_default.py <arg0> <arg1> <arg2> <arg3> <arg4>'
	date = "05-01-2022"
	str_type = "all"
	stepname = "cesymm_rmsd_3"
	locusts_tmp = "/data/local/encompass/symmetry/locusts_" + stepname + "_" + str_type + "/"
	job_code = stepname + "_" + str_type + "_" + date

	collect_results(locations_path, locations, job_code, stepname, locusts_tmp)
	#run_locusts(fn_opt, locations_path, locations, 'encompass-06-01-2021', 'cesymm', locusts_tmp, pdb_list, command)
