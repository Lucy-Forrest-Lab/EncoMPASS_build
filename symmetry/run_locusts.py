import collections
import locusts.swarm
import sys
import os
sys.path.append(os.path.expandvars('$ENCSOURCESYM'))

def write_locusts_parfile(fn_options, locusts_filename, job_code, gather):
	text = (
		"### Generic\n"
		"{0:30}\t{1}\n\n"
	).format('run_on_hpc', fn_options['run_on_hpc'])

	if fn_options['run_on_hpc']:
		text += "### HPC\n"

		for k in ['host_name', 'partition', 'requested_nodes', 'cpus_per_node', 'constraint',
				  'min_stack_per_core', 'exclusive', 'walltime', 'waiting_time', 'memory',
				  'memory_per_cpu', 'data_transfer_protocol', 'nodewise_scratch_folder',
				  'nodewise_scratch_memory', 'email_address', 'extra_outer_statements']:
			if fn_options[k]:
				text += "{0:30}\t{1}\n".format(k, fn_options[k])
		text += "{0:30}\t{1}\n".format('hpc_exec_dir', fn_options[
			'hpc_exec_path'] + '/' + job_code)  # EncoMPASS hpc_exec_path is the parent of each locusts hpc_exec_dir

		if fn_options['hpc_singularity'] != '':
			text += "\n### Singularity\n"
			for k in ['singularity', 'singularity_container', 'singularity_modload']:
				text += "{0:30}\t{1}\n".format(k, fn_options['hpc_' + k])
	else:
		text += (
			"### Local multithreading\n"
			"{0:30}\t{1}\n\n"
		).format('number_of_processors', fn_options['number_of_processors'])
		if fn_options['singularity'] != '':
			text += "\n### Singularity\n"
			for k in ['singularity', 'singularity_container', 'singularity_modload', 'extra_outer_statements']:
				text += "{0:30}\t{1}\n".format(k, fn_options[k])

	if gather == True:
		text += "{0:30}\t{1}\n".format("only_gather", "True")

	with open(locusts_filename, 'w') as outf:
		outf.write(text)

	return

def default_fn_options(options):
	fn_opt = collections.OrderedDict()
	fn_opt['run_on_hpc'] = 'True'
	fn_opt['host_name'] = options['ALL']['host_name']
	fn_opt['partition'] = 'norm'
	fn_opt['requested_nodes'] = '8'
	fn_opt['cpus_per_node'] = '28'
	fn_opt['exclusive'] = 'True'
	fn_opt['memory'] = '247g'
	fn_opt['memory_per_cpu'] = ''
	fn_opt['constraint'] = 'x2695|x2680'
	fn_opt['min_stack_per_core'] = '1' # the minimum number of tasks that a core needs to complete (queue on it)
	fn_opt['walltime'] = '3-00:00:00'
	fn_opt['waiting_time'] = '600' # time [s] between tasks check on HPC (requires ssh login)
	fn_opt['hpc_exec_path'] = options['ALL']['hpc_exec_path']
	fn_opt['local_shared_dir'] = ''
	fn_opt['nodewise_scratch_folder'] = ''
	fn_opt['nodewise_scratch_memory'] = ''
	fn_opt['email_address'] = 'antoniya.aleksandrova@nih.gov'
	fn_opt['extra_outer_statements'] = 'export SINGULARITY_BINDPATH="$SINGULARITY_BINDPATH,/data/Encompass/"'
	fn_opt['number_of_processors'] = options['ALL']['number_of_processors'] #for local runs
	fn_opt['singularity'] = options['ALL']['singularity']
	fn_opt['singularity_modload'] = options['ALL']['singularity_modload']
	fn_opt['singularity_container'] = options['ALL']['singularity_container']
	fn_opt['hpc_singularity'] = options['ALL']['hpc_singularity']
	fn_opt['hpc_singularity_modload'] = options['ALL']['hpc_singularity_modload']
	fn_opt['hpc_singularity_container'] = options['ALL']['hpc_singularity_container']
	fn_opt['data_transfer_protocol'] = options['ALL']['data_transfer_protocol']

	return fn_opt

def write_locusts_fst(locusts_fstname, dbname, stepname):
	f = open(locusts_fstname, 'w')
	text = '''\
#WORKDIR : {wkdir}

{db} : *
	database :
		selection :
			whole_structs : *
	symmetries :
		{step} : **
		temp : **
	logs :
		symmetry : *
			.tm_archive.pkl
			{step} : **
	! summary_table.txt
encsource : **
	! __pycache__\
	'''.format(wkdir = 'encsource/symmetry/', db = dbname, step = stepname)
	f.write(text)

def write_locusts_cesymm_symd_fst(locusts_fstname, dbname, stepname):
	f = open(locusts_fstname, 'w')
	text = '''\
#WORKDIR : {wkdir}

{db} : *
	database :
		selection :
			whole_structs : *
	symmetries :
		cesymm : **
		symd: **
		temp : **
	logs :
		symmetry : *
			cesymm : **
			symd : **
			{step} : **
			.tm_archive.pkl
	webarchive: **
	! summary_table.txt
encsource : **
	! __pycache__\
	'''.format(wkdir = 'encsource/symmetry/', db = dbname, step = stepname)
	f.write(text)

def write_locusts_mssd_fst(locusts_fstname, dbname):
	f = open(locusts_fstname, 'w')
	text = '''\
#WORKDIR : {wkdir}

{db} : *
	database :
		selection :
			whole_structs : *
	symmetries : **
	logs :
		symmetry : **
			.tm_archive.pkl
	webarchive : **
	! summary_table.txt
encsource : **
	! __pycache__\
	'''.format(wkdir = 'encsource/symmetry/', db = dbname)
	f.write(text)

def write_locusts_neighbors_fst(locusts_fstname, dbname, stepname):
	f = open(locusts_fstname, 'w')
	text = '''\
#WORKDIR : {wkdir}

{db} : 
	EncoMPASS_options_relative_to__instruction_file.txt
	summary_table.txt
	symmetries :
		{step} : **
	logs :
		symmetry :
			.tm_archive.pkl
			{step} : **
encsource : **
	! __pycache__\
	'''.format(wkdir = 'encsource/symmetry/', db = dbname, step = stepname)
	f.write(text)

def write_locusts_transfer_fst(locusts_fstname, dbname, stepname):
	f = open(locusts_fstname, 'w')
	text = f'''\
#WORKDIR : encsource/symmetry/

{dbname} : 
	EncoMPASS_options_relative_to__instruction_file.txt
	database :
		selection :
			whole_structs : *
		alignments :
			str_seqalns : **
	symmetries : 
		cesymm : **
		cesymm_from_symd : **
		cesymm_low_thr : **
		cesymm_minlen : **
		cesymm_order : **
		cesymm_rmsd_2_5 : **
		cesymm_rmsd_3 : **
		mssd : **
		symd : **
		temp :
			{stepname} : **
		{stepname} : **			
	logs :
		symmetry : 
			{stepname} : **
			.tm_archive.pkl
	webarchive : 
		{stepname}_chains_json : **
		{stepname}_chains_super : **
encsource : **
	! __pycache__\

	'''
	f.write(text)
	f.close()

def write_locusts_images_transfer_fst(locusts_fstname, dbname, stepname):
	f = open(locusts_fstname, 'w')
	text = f'''\
#WORKDIR : encsource/symmetry/

{dbname} : 
	EncoMPASS_options_relative_to__instruction_file.txt
	database :
		selection :
			whole_structs : *
	symmetries :
		{stepname} : **
		transfer : **
		temp : 
			transfer : **
	logs :
		symmetry : *
			.tm_archive.pkl
			{stepname} : **
	webarchive : 
		transfer_chains_json : **
		transfer_chains_pngs : **
		transfer_chains_super : **
encsource : **
	! __pycache__\
	'''
	f.write(text)

def run_locusts(options, locations_path, locations, job_code, stepname, local_dir, pdb_list, command_template, opt_arg = [''], gather = False): #, str_data, str_set, only_analysis=False
	"""
	Launches Locusts in a manner compatible with the symmetry runs

	:param options: note that these options are the fn_options, not the database ones
	:param locations_path: path to the EncoMPASS database version being used
	:param locations: dictionary of locations within the EncoMPASS database
	:param job_code: the name of the HPC folder and the label attached to each log file from the run
	:param stepname: the name of the procedure being submitted (e.g. 'cesymm_low_thr')
	:param local_dir: the folder on the local machine that can be used to store intermediary locusts files
	:param pdb_list: list of structures that will be analyzed
	:param command_template: as required by locusts; e.g. "python cesymm_default.py <arg0> <arg1> <arg2> <arg3> <arg4>"
	:param opt_arg: can be left empty (default) or given a values such as 3 for the maxrmsd in cesymm_lt_rmsd runs
	:param gather: use locusts only to collect completed results from biowulf back to the local machine
	:return:
	"""

	db_loc = locations_path.split('/')
	parname = locations['FSYSPATH']['locusts_log'] + stepname + '_locusts.par'
	print(options, parname)
	write_locusts_parfile(options, parname, job_code, gather)
	fstname = locations['FSYSPATH']['locusts_log'] + stepname + '_locusts.fst'
	if stepname == "mssd":
		write_locusts_mssd_fst(fstname, db_loc[-1])
	elif stepname == "images_cesymm_symd":
		write_locusts_cesymm_symd_fst(fstname, db_loc[-1], stepname)
	elif stepname == "neighbors":
		write_locusts_neighbors_fst(fstname, db_loc[-1], stepname)
	elif stepname == "transfer":
		write_locusts_transfer_fst(fstname, db_loc[-1], stepname)
	elif stepname == "images_transfer":
		write_locusts_images_transfer_fst(fstname, db_loc[-1], stepname)
	else:
		write_locusts_fst(fstname, db_loc[-1], stepname)

	batch_job_code = stepname + "Run"  # A unique identifier of your choice
	env_root = "/".join(db_loc[:-1]) + "/"  # the root of the tree is .fst
	if options['run_on_hpc']:
		run_locations_path = options['hpc_exec_path'] + job_code + "/" + batch_job_code + "/" + db_loc[-1]
	else:
		run_locations_path = locations_path
	print(run_locations_path)
	pdb_suff = '_sb'
	if "transfer" in stepname:
		argument_list = [[run_locations_path, template, target, pdb_suff, job_code] + opt_arg for (template,target) in pdb_list]
	else:
		print("opt_arg: ", opt_arg, type(opt_arg))
		argument_list = [[run_locations_path, pdb, pdb_suff, job_code] + opt_arg for pdb in pdb_list]
	#print('arg_list:', argument_list)
	locusts.swarm.launch(
		code=batch_job_code,
		cmd=command_template,
		args=argument_list,
		envroot=env_root,
		envfs=fstname,
		parf=parname,
		locdir=local_dir
	 )

	print(run_locusts.__name__ + ' completed.')


if __name__ == "__main__":
	from initialize_repository import initialize_repository
	locations_path = os.path.expandvars("ENC_DB")
	options, locations = initialize_repository(main_path=locations_path)
	
	fn_opt = default_fn_options(options)
	locusts_tmp = '~/test_locusts/' # tmp dir for locusts on local machine
	pdb_list = ['6kig', '6dr2']
	fn_opt['requested_nodes'] = '2'
	fn_opt['cpus_per_node'] = '1'
	fn_opt['walltime'] = '3-00:00:00'
	command = 'source ~/encompass/bin/activate; python cesymm_default.py <arg0> <arg1> <arg2> <arg3> <arg4>'
	run_locusts(fn_opt, locations_path, locations, 'encompass-06-01-2021', 'cesymm', locusts_tmp, pdb_list, command)
	
	