import collections
import locusts.swarm
import sys
import os
import shutil
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
				  'nodewise_scratch_space', 'email_address', 'extra_outer_statements']:
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
			for k in ['singularity', 'singularity_container', 'singularity_modload']:
				text += "{0:30}\t{1}\n".format(k, fn_options[k])

	if gather == True:
		text += "{0:30}\t{1}\n".format("only_gather", "True")

	with open(locusts_filename, 'w') as outf:
		outf.write(text)

	return

def default_fn_options(options):
   fn_opt = collections.OrderedDict()
   fn_opt['run_on_hpc'] = options['ALL']['run_on_hpc']
   fn_opt['host_name'] = options['ALL']['host_name']
   fn_opt['partition'] = options['ALL']['partition']
   fn_opt['requested_nodes'] = options['ALL']['nodes']
   fn_opt['cpus_per_node'] = options['ALL']['cpus_per_node']
   fn_opt['exclusive'] = options['ALL']['exclusive_nodes']
   fn_opt['memory'] = str(options['ALL']['nodemem']) + 'g' if str(options['ALL']['nodemem']).strip() else ""
   fn_opt['memory_per_cpu'] = ''
   fn_opt['constraint'] = ''
   fn_opt['min_stack_per_core'] = '1' # PER CORE! not per node.
   fn_opt['walltime'] = options['ALL']['walltime']
   fn_opt['waiting_time'] = '600' # time [s] between tasks check on HPC (requires ssh login)
   fn_opt['hpc_exec_path'] = options['ALL']['hpc_exec_path']
   fn_opt['local_shared_dir'] = ''
   fn_opt['nodewise_scratch_folder'] = options['ALL']['nodescratch']
   fn_opt['nodewise_scratch_space'] = options['ALL']['nodescratchspace']
   fn_opt['email_address'] = options['ALL']['email_address']
   fn_opt['extra_outer_statements'] = options['ALL']['extra_outer_statements'] 
   fn_opt['singularity'] = options['ALL']['singularity']
   fn_opt['singularity_modload'] = options['ALL']['singularity_modload']
   fn_opt['singularity_container'] = options['ALL']['singularity_container']
   fn_opt['hpc_singularity'] = options['ALL']['hpc_singularity']
   fn_opt['hpc_singularity_modload'] = options['ALL']['hpc_singularity_modload']
   fn_opt['hpc_singularity_container'] = options['ALL']['hpc_singularity_container']
   fn_opt['data_transfer_protocol'] = options['ALL']['data_transfer_protocol']

   return fn_opt


def local_fn_options(options):
   fn_opt = collections.OrderedDict()
   fn_opt['run_on_hpc'] = False
   fn_opt['number_of_processors'] = options['ALL']['number_of_processors'] #for local runs
   fn_opt['singularity'] = ''
   fn_opt['singularity_modload'] = options['ALL']['singularity_modload']
   fn_opt['singularity_container'] = options['ALL']['singularity_container']

   return fn_opt



def write_locusts_fst(locusts_fstname, db_name, an_path, workdir, not_profiles=False):
    # CONVENTION: There is an EncoMPASS symlink just outside the main folder of the database
    # NB: workdir must be in the database
    # NB: an_path must be a FSYS path

    an_path_list = [x for x in an_path.split('/') if x]

    f = open(locusts_fstname, 'w')
    text = f'''\
#WORKDIR : {workdir}

EncoMPASS : *
	! lucy_enc_env
	! __pycache__
{db_name} :'''
    for ix, x in enumerate(an_path_list):
        text += '\n' + '\t'*(ix+1) + x + ' : '
    text += '**\n'

    if not_profiles:
        text += '\t'*(ix+1) + '! distributions/distributions_data\n'

    f.write(text)


def write_locusts_fst_leslie(locusts_fstname, env_root, workdir, an_path):

    an_path = rebase_path(an_path, env_root)
    an_path_list = [x for x in an_path.split('/') if x]

    workdir = rebase_path(workdir, env_root)

    f = open(locusts_fstname, 'w')
    text = f'''\
#WORKDIR : {workdir}

EncoMPASS : *.py
'''
    for ix, x in enumerate(an_path_list):
        text += '\n' + '\t'*(ix) + x + ' : '
    text += '**\n'

    f.write(text)


def write_locusts_fst_shulien(locusts_fstname, env_root, workdir, an_path):

    an_path = rebase_path(an_path, env_root)
    an_path_list = [x for x in an_path.split('/') if x]

    workdir = rebase_path(workdir, env_root)

    f = open(locusts_fstname, 'w')
    text = f'''#WORKDIR : {workdir}

code :
	Edo :
		EncoMPASS : *.py
'''
    for ix, x in enumerate(an_path_list):
        text += '\n' + '\t\t\t'*(ix) + x + ' : '
    text += '**\n'

    f.write(text)


def rebase_path(path, new_base, permissive=False):
    trailslash = ''
    if path.strip()[-1] == '/':
        trailslash = '/'
    if new_base.count('/') == 1 and new_base[-1] == '/':
        new_base = new_base[:-1]
    if '/' not in new_base:
        path_list = [x for x in path.split('/') if x]
        trailslash = ''
        if path.strip()[-1] == '/':
            trailslash = '/'
        if new_base not in path_list:
            print("ERROR: new_base not in path_list")
            print(new_base)
            print(path_list)
            if permissive: 
                return False
            else:
                exit(1)
        if len([x for x in path.split('/') if x == new_base]) > 1:
            print("new_base appears multiple times in path")
            print(new_base)
            print(path)
            if permissive:
                return False
            else:
                exit(1)
        else:
            path_list = path_list[path_list.index(new_base)+1:]
        rebased_path = '/'.join(path_list) + trailslash
    else:
        if new_base not in path:
            print("new_base not in path")
            print(new_base)
            print(path)
            if permissive:
                return False
            else:
                exit(1)
        if path.index(new_base) != 0 :
            print("if new_base is a path, path must start with it")
            print(new_base)
            print(path)
            if permissive:
                return False
            else:
                exit(1)
        rebased_path = path.replace(new_base, '')
        rebased_path = '/'.join([x for x in rebased_path.split('/') if x]) + trailslash
    if not rebased_path.replace('/', ''):
        rebased_path = ''
    return rebased_path


def run_locusts(options, locations, workdir, job_code, stepname, local_dir, pdb_lists, command_template, fixed_args = [], gather = False): #, str_data, str_set, only_analysis=False
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
    :param gather: use locusts only to cocllect completed results from biowulf back to the local machine
    :return:
    """
    # NB: workdir must be an absolute path!
    
    # The env_root is one level up from the Encompass database main directory
    # NB: in env_root there MUST be a symlink to the script folder, and it MUST be named "EncoMPASS"
    env_root = '/'.join([x for i, x in enumerate(locations['FSYSPATH']['main'].split('/')) if x != "" or i == 0][:-1])

    an_path = locations['FSYS']['analysis'] #parent
    ancache_path = locations['FSYSPATH']['ancache'] #child
    if os.path.exists(ancache_path):
        shutil.rmtree(ancache_path)
    os.mkdir(ancache_path)

    parname = ancache_path + stepname + '_locusts.par'
    write_locusts_parfile(options, parname, job_code, gather)
    
    locusts_fstname = ancache_path + stepname + '_locusts.fst'
    db_name = [x for x in locations['FSYSPATH']['main'].split('/') if x][-1]
    print("DB NAME", db_name)
    print("WORKDIR", workdir)
    print("ANPATH", an_path)
    rebased_workdir = rebase_path(workdir, env_root)
    print("REBASED_WORKDIR", rebased_workdir)
    write_locusts_fst(locusts_fstname, db_name, an_path, rebased_workdir)
    workdir_prepend = ''.join(['../' for x in rebased_workdir if x == '/'])    

    print("rebased_workdir", rebased_workdir)
    print("workdir_prepend", workdir_prepend)

    batch_job_code = stepname + "Run"  # A unique identifier of your choice

    print("LISTS", len(pdb_lists))
    print("ALL", len(sum(pdb_lists, [])))

    fn_list = []
    for il, pdb_list in enumerate(pdb_lists):
        fn_local = ancache_path + 'an_list_{0}'.format(str(il).zfill(3))
        fn_runtime = workdir_prepend + rebase_path(ancache_path, env_root) + 'an_list_{0}'.format(str(il).zfill(3))
        with open(fn_local, 'w') as f:
            for pdb_ch in pdb_list:
                f.write(pdb_ch+'\n')
        fn_list.append(fn_runtime)

    print("FNS", len(fn_list))

    print("ROOT FST", env_root)
    print("PARNAME", parname)
    print("FSNAME", locusts_fstname)

    rebased_fixed_args = [workdir_prepend + rebase_path(a, env_root) for a in fixed_args]
    argument_list = [[*rebased_fixed_args, fn] for fn in fn_list]

    command_template = " ".join([workdir_prepend + x if x.startswith('EncoMPASS') else x for x in command_template.split()])
    print("LOCDIR", local_dir)
    print("# ARGLIST", argument_list)
    print("# will use command:")
    print("# COMMAND", command_template)
    print("# launching locusts now")

    locusts.swarm.launch(
        code=batch_job_code,
        cmd=command_template,
        args=argument_list,
        envroot=env_root,
        envfs=locusts_fstname,
        parf=parname,
        locdir=local_dir
    )
	
    print(run_locusts.__name__ + ' completed.')


if __name__ == "__main__":
    tests = [
       ("/absolute/path/rebase_here/a/b/c/", '/absolute/path/rebase_here/', 'a/b/c/'),
       ("/absolute/path/rebase_here/a/b/c/", 'rebase_here', 'a/b/c/'),
       ("/absolute/path/rebase_here/a/b/c/", 'base_here', False),
       ("/absolute/path/rebase_here/a/b/c/", 'rebase_here/', 'a/b/c/'),
       ("/absolute/path/rebase_here/a/b/c/", '/rebase_here/', False),
       ("relative/path/rebase_here/a/b/c/", 'relative/path/rebase_here/', 'a/b/c/'),
       ("relative/path/rebase_here/a/b/c/", 'rebase_here', 'a/b/c/'),
       ("relative/path/rebase_here/a/b/c/", 'base_here', False),
       ("relative/path/rebase_here/a/b/c/", 'rebase_here/', 'a/b/c/'),
       ("relative/path/rebase_here/a/b/c/", '/rebase_here/', False),
       ("/absolute/path/ambiguous/path/a/b/c/", 'path', False),
       ("relative/path/ambiguous/path/a/b/c/", 'path', False),
       ("/absolute/path/ambiguous/path/a/b/c/", '/absolute/path/ambiguous/path', 'a/b/c/'),
       ("relative/path/ambiguous/path/a/b/c/", 'relative/path/ambiguous/path/', 'a/b/c/'),
    ]

    for path_to_rebase, new_base, exp_result in tests:
        print("PATH_TO_REBASE", path_to_rebase)
        print("NEW_BASE", new_base)
        res = rebase_path(path_to_rebase, new_base, permissive=True)
        print("EXPECTED_RESULT", exp_result, "REBASED_PATH", res)
        print("OUTCOME", exp_result == res, "\n\n\n")
