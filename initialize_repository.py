# Name: initialize_repository.py
# Version: EncoMPASS 2.0
# Language: python3
# Description:
# Author: Edoardo Sarti
# Date: 2018/10/28

import argparse
import collections
import datetime
from supporting_functions import *

class FooAction(argparse.Action):
	def __call__(self, parser, namespace, values, option_string=None):
		setattr(namespace, self.dest, values)
		setattr(namespace, self.dest+'_nondefault', True)
class FooAction2(argparse.Action):
	def __call__(self, parser, namespace, values, option_string=None):
		setattr(namespace, self.dest, True)


def load_options(correct_types=False, thr_log_status='ERROR'):
    # Option sections
    options = collections.OrderedDict()
    options['GENERAL'] = collections.OrderedDict() 
    options['PATHS'] = collections.OrderedDict()
    options['RUN'] = collections.OrderedDict()
    options['PARAMETERS'] = collections.OrderedDict()
    options['EXECUTION'] = collections.OrderedDict()
    options['ALL'] = collections.OrderedDict()

    # GENERAL options *must* be set (no defaults)
    options['GENERAL'][('s', 'instruction_filename')] = ''
    options['GENERAL'][('', 'out_instruction_filename')] = ''
    options['GENERAL'][('main', 'main_directory')] = ''	# If it does not exist, it will create one

    # PATH options *must* be set (no defaults)
    options['PATHS'][('dst', 'data_structure_template')] = os.path.dirname(os.path.abspath(__file__)) + '/str_data_entry_current.json'
    options['PATHS'][('ppm', 'ppm_path')] = ''
    options['PATHS'][('muscle', 'muscle_path')] = ''
    options['PATHS'][('cesymm', 'cesymm_path')] = ''
    options['PATHS'][('symd', 'symd_path')] = ''
    options['PATHS'][('frtmalign', 'frtmalign_path')] = ''
    options['PATHS'][('pymol', 'pymol_path')] = ''
    options['PATHS'][('dssp', 'dssp_path')] = ''
    options['PATHS'][('mpiq', 'mpiq_path')] = ''
    options['PATHS'][('sing', 'singularity')] = ''
    options['PATHS'][('singmod', 'singularity_modload')] = ''
    options['PATHS'][('container', 'singularity_container')] = ''
    options['PATHS'][('hpcsing', 'hpc_singularity')] = ''
    options['PATHS'][('hpcsingmod', 'hpc_singularity_modload')] = ''
    options['PATHS'][('hpccontainer', 'hpc_singularity_container')] = ''
    options['PATHS'][('hpcppm', 'hpc_ppm_path')] = ''
    options['PATHS'][('hpcmuscle', 'hpc_muscle_path')] = ''
    options['PATHS'][('hpccesymm', 'hpc_cesymm_path')] = '/data/Encompass/software/cesymm-2.2.3-SNAPSHOT/runCESymm.sh'
    options['PATHS'][('hpccesymmsmall', 'hpc_cesymm_small_path')] = '/data/Encompass/software/cesymm-2.2.3-SNAPSHOT/runCESymm_small.sh'
    options['PATHS'][('hpcquatsymm', 'hpc_quatsymm_path')] = ''
    options['PATHS'][('hpcsymd', 'hpc_symd_path')] = ''
    options['PATHS'][('hpcfrtmalign', 'hpc_frtmalign_path')] = ''
    options['PATHS'][('hpcdssp', 'hpc_dssp_path')] = ''
    options['PATHS'][('sigppm', 'sig_ppm_path')] = '/EncoMPASS/ppm/immers'
    options['PATHS'][('sigmuscle', 'sig_muscle_path')] = '/EncoMPASS/muscle3.8.31/muscle3.8.31_i86linux64'
    options['PATHS'][('sigcesymm', 'sig_cesymm_path')] = '/EncoMPASS/cesymm-2.1.0/runCESymm.sh'
    options['PATHS'][('sigcesymmsmall', 'sig_cesymm_small_path')] = '/EncoMPASS/cesymm-2.1.0/runCESymm_small.sh'
    options['PATHS'][('sigsymd', 'sig_symd_path')] = '/EncoMPASS/symd1.61/src/symd1.61-linux'
    options['PATHS'][('sigquatsymm', 'sig_quatsymm_path')] = '/EncoMPASS/quatsymm-2.2.3-SNAPSHOT/runQuatSymm.sh'
    options['PATHS'][('sigfrtmalign', 'sig_frtmalign_path')] = '/EncoMPASS/frtmalign/frtmalign'
    options['PATHS'][('sigdssp', 'sig_dssp_path')] = '/EncoMPASS/dssp-3.1.4/bin/mkdssp'
    options['PATHS'][('waout', 'wa_out_dir')] = ''
    options['PATHS'][('wa', 'webarchive')] = ''
    options['PATHS'][('me', 'missing_entries')] = ''
    options['PATHS'][('', 'ppm_reslib_path')] = ''

    
    # RUN options *must* be set, and have a default
    options['RUN'][('hpc', 'run_on_hpc')] = 'False'
    options['RUN'][('np', 'number_of_processors')] = '1'
    options['RUN'][('code', 'main_code')] = 'ENC' + datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d') + 'XXXX'
    options['RUN'][('host', 'host_name')] = ''
    options['RUN'][('partition', 'hpc_partition')] = ''
    options['RUN'][('constraint', 'hpc_constraint_processors')] = ''
    options['RUN'][('walltime', 'hpc_walltime')] = '24:00:00'
    options['RUN'][('hpctpc', 'min_stack_per_core')] = '10'
    options['RUN'][('hpcdir', 'hpc_exec_path')] = ''
    options['RUN'][('shdir', 'local_shared_dir')] = ''
    options['RUN'][('nodescratch', 'nodewise_scratch_folder')] = ''
    options['RUN'][('nodescratchspace', 'nodewise_scratch_space')] = ''
    options['RUN'][('nodescratchpath', 'nodewise_scratch_path')] = ''
    options['RUN'][('email', 'email_address')] = ''
    options['RUN'][('extracmds', 'extra_outer_statements')] = ''
    options['RUN'][('nodes', 'requested_nodes')] = ''
    options['RUN'][('nodemem', 'memory_per_node')] = ''
    options['RUN'][('excl', 'exclusive_nodes')] = 'False'
    options['RUN'][('cpn', 'cpus_per_node')] = ''
    options['RUN'][('dtp', 'data_transfer_protocol')] = ''
    options['RUN'][('ol', 'offline')] = 'False'
    options['RUN'][('fud', 'force_update')] = 'False'
    options['RUN'][('', 'download_date_threshold')] = '30'
    options['RUN'][('repo', 'from_repo')] = ''
    options['RUN'][('pdbmirror', 'pdb_local_mirror')] = ''
    options['RUN'][('mmcifmirror', 'mmcif_local_mirror')] = ''   
    options['RUN'][('hpcactipy', 'hpc_python_activate_command')] = ''

 
    # PARAMETERS options *must* be set, and have a default
    options['PARAMETERS'][('res_thr', 'structure_resolution_threshold')] = ''	# If 0, it is None
    options['PARAMETERS'][('hole_thr', 'structure_hole_threshold')] = '100'	# number of residues
    options['PARAMETERS'][('holefrac_thr', 'structure_hole_fraction_threshold')] = '0.3'
    options['PARAMETERS'][('seqid_thr', 'sequence_identity_threshold')] = '0.85'
    options['PARAMETERS'][('tmscore_thr', 'tmscore_threshold')] = '0.6'
    options['PARAMETERS'][('no_nmr', 'no_nmr')] = 'False'
    options['PARAMETERS'][('no_cryoEM', 'no_cryoEM')] = 'False'
    options['PARAMETERS'][('no_theoretical', 'no_theoretical')] = 'False'
    options['PARAMETERS'][('', 'hole_thr')] = '100'
    options['PARAMETERS'][('', 'hole_frac_thr')] = '0.3'
    options['PARAMETERS'][('', 'disordered_thr')] = '100' # Deactivated
    options['PARAMETERS'][('', 'intra_persub_clashes_thr')] = '100'
    options['PARAMETERS'][('', 'inter_persub_clashes_thr')] = '30'
    options['PARAMETERS'][('', 'clashes_thr')] = '100'
    options['PARAMETERS'][('', 'clash_frac_thr')] = '0.1'
    options['PARAMETERS'][('', 'clash_distance')] = '2.5'
    options['PARAMETERS'][('', 'with_UNK')] = 'False'
    options['PARAMETERS'][('', 'consensus_signature_PDB')] = 'True'
    options['PARAMETERS'][('', 'use_PDB_format')] = 'True'
    options['PARAMETERS'][('squash', 'make_squash_filesystem')] = 'False'
    
    # EXECUTION options are mainly for debugging
    options['EXECUTION'][('debug', 'debug')] = 'False'
    options['EXECUTION'][('paraldebug', 'parallel_debug')] = 'False'
    options['EXECUTION'][('del_orig_instr', 'delete_original_instructions')] = 'False'
    options['EXECUTION'][('no_related', 'neglect_opm_related_structures')] = 'False'
    options['EXECUTION'][('locp', 'location_permissive')] = 'True'
    options['EXECUTION'][('del', 'deletion_entries')] = ''
    options['EXECUTION'][('repl', 'replace_entries')] = ''
    options['EXECUTION'][('add', 'additional_entries')] = ''
    options['EXECUTION'][('forcedb', 'force_db_choice')] = ''
    options['EXECUTION'][('redo', 'recompute_list')] = '[]'
    options['EXECUTION'][('update', 'update_all')] = 'False'
    options['EXECUTION'][('sample', 'sample_structures')] = '[]'
    options['EXECUTION'][('', 'execute_list')] = '[]'
    options['EXECUTION'][('ext', 'external_list')] = '[]'
    options['EXECUTION'][('adofn', 'analysis_do_only_fn')] = ''
    options['EXECUTION'][('', 'no_neighbor_lists_analysis')] = 'False'
    options['EXECUTION'][('', 'filter_entries_analysis')] = 'False'    
    options['EXECUTION'][('', 'batch_size_analysis')] = '10'
    options['EXECUTION'][('', 'no_structurewisetable_analysis')] = 'True'
    options['EXECUTION'][('', 'quick_test_analysis')] = 'False'
    options['EXECUTION'][('', 'sort_by_chain_size_analysis')] = 'False'
    options['EXECUTION'][('', 'sort_by_num_neighbors_analysis')] = 'True'


    for x in options:
        if x == 'ALL':
            continue
        for y in options[x]:
            y1, y2 = y
            options['ALL'][y1] = (x, y)
            options['ALL'][y2] = (x, y)
    
    if correct_types:
        new_options = {'ALL' : {}}
        for x in options:
            if x != 'ALL':
                new_options[x] = {}
                for y in options[x]:
                    new_options[x][y] = string_decode_element(options[x][y], permissive=True)
                    y1, y2 = y
                    new_options['ALL'][y1] = (x, y)
                    new_options['ALL'][y2] = (x, y)
        options = new_options
	
    return options	


def command_line_parser(thr_log_status='ERROR'):
	this_name = command_line_parser.__name__

	options = load_options(thr_log_status=thr_log_status)

	# 1. Parsing from command line with argparse
	parser = argparse.ArgumentParser()

	falsedef = set()
	for s_name, l_name in set([options['ALL'][x][1] for x in options['ALL']]):
		if not s_name:
			continue
		default = options[options['ALL'][l_name][0]][options['ALL'][l_name][1]]
		if default == 'False':
			parser.add_argument('-'+s_name, '--'+l_name, action='store_const', default=default, const='True')
			falsedef.add(l_name)
		else:
			parser.add_argument('-'+s_name, '--'+l_name, nargs='?', default=default, action=FooAction)

	parsed, unknown = parser.parse_known_args()


	# 2. Transfer into options dictionary, restoring correct variable types
	options = {}
	for x in parsed.__dict__:
		if not hasattr(parsed, x+'_nondefault') and x not in falsedef:
			continue
		if x in parsed.__dict__ and parsed.__dict__=='':
			options[x] = ''
			continu
		if x in falsedef:
			optx = string_decode_element(parsed.__dict__[x], permissive=True)
			if optx:
				options[x] = True
		else:
			options[x] = string_decode_element(parsed.__dict__[x], permissive=True)


	deb_msg = ('DEBUG', this_name, "These arguments were parsed from the command line: {0}".format(options))
	print_log(deb_msg)

	return options


def instruction_file_parser(instruction_filename, thr_log_status='ERROR'):
	this_name = instruction_file_parser.__name__

	if not instruction_filename:
		crt_msg = ('CRITICAL', this_name, "An instruction file must be specified to run EncoMPASS 2.0")
		print_log(crt_msg)

	def_options = load_options(thr_log_status=thr_log_status)
	options = {}
	field0, field1 = "", ""
	with open(instruction_filename, 'r') as instruction_file:
		for l in instruction_file:
			if not l.strip() or l.startswith("#"):
				continue
			if "#" in l:
				line = l[:l.index("#")]
			else:
				line = l
			# -- from here on, the loop variable is "line" -- #
			if not line.strip():
				continue
			# If line starts with spaces, it is a continuation line (field0 is empty)
			if line.lstrip() != line:
				fields = ["", " ".join([x for x in line.split()]).strip()]
			else:
				fields = [line.split()[0], " ".join([x for x in line.split()[1:]]).strip()]
			if fields[0]:
				if field1.strip():
					options[field0] = string_decode_element(field1.strip())
					field1 = ""
				# Do not parse (presumably deprecated) instruction filename from instruction file...
				elif fields[0] == "s":
					continue
				# If field0 not found in option dictionary, interrupt
				elif fields[0] not in def_options['ALL']:
					wfl_msg = ('CRITICAL', this_name, 'option {0} from instruction file is not valid'.format(fields[0]))
					print_log(wfl_msg)
				field0 = fields[0]
			# If field1 is not present, the default value will be applied
			if not fields[1].strip():
				wfl_msg = ('NOTICE', this_name, 'option {0} from instruction file does not have an attribute: the command line or default value will be applied'.format(field0))
				print_log(wfl_msg)
			else:
				field1 += fields[1].lstrip()
		if field1:
			options[field0] = string_decode_element(field1.strip())
			field1 = ""

	return options


def instruction_file_writer(options, out_instruction_filename, thr_log_status='ERROR'):
	this_name = instruction_file_writer.__name__

	out_instruction_file = open(out_instruction_filename, 'w')
	ts = time.time()
	out_instruction_file.write("# AUTOMATICALLY GENERATED INSTRUCTION FILE {0}\n# Reports all the options with which EncoMPASS is being run\n\n".format(datetime.datetime.fromtimestamp(ts).strftime('%Y/%m/%d %H:%M:%S')))
	for sector in options:
		if sector == 'ALL':
			continue
		out_instruction_file.write("\n# {0}\n".format(sector))
		for s_key, l_key in options[sector]:
			sk = ''+s_key
			if not sk:
				sk = l_key
			is_string = type(options[sector][(s_key, l_key)]) == str
			phrase = str(options[sector][(s_key, l_key)]).split()
			if is_string and options[sector][(s_key, l_key)].strip():
				phrase[0] = "'" + phrase[0]
				phrase[-1] += "'"
			s = ''
			line = " ".join(phrase)
			for nc, c in enumerate(line):
				s += c
				if nc != 0 and len(s) % 80 == 0:
					out_instruction_file.write("{0:30}{1:80}\n".format(sk, s))
					s = ''
					sk = ''
			if s:
				out_instruction_file.write("{0:30}{1:80}# {2}\n".format(sk, s, "".join([x+" " for x in l_key.split('_')])))
	out_instruction_file.close()


def check_options(command_line_options, instruction_file_options, thr_log_status='ERROR'):
	this_name = check_options.__name__

	def_options = load_options(correct_types=True, thr_log_status=thr_log_status)
	options = collections.OrderedDict()
	for x in def_options['ALL']:
		f0 = def_options['ALL'][x][0]
		f1 = def_options['ALL'][x][1]
		if f0 not in options:
			options[f0] = collections.OrderedDict()
		options[f0][f1] = def_options[f0][f1]
	for x in instruction_file_options:
		if x not in def_options['ALL']:
			wfl_msg = ('CRITICAL', this_name, 'option {0} from instruction file is unknown'.format(x))
			print_log(wfl_msg)
		if_tp = type(instruction_file_options[x])
		df_tp = type(def_options[def_options['ALL'][x][0]][def_options['ALL'][x][1]])
		if if_tp != df_tp and def_options[def_options['ALL'][x][0]][def_options['ALL'][x][1]]:
			wtp_msg = ('CRITICAL', this_name, 'in instruction file, option {0} is of type {1}, but it should be of type {2}'.format(x, if_tp, df_tp))
			print_log(wtp_msg)
		else:
			f0 = def_options['ALL'][x][0]
			f1 = def_options['ALL'][x][1]
			if f0 not in options:
				options[f0] = collections.OrderedDict()
			options[f0][f1] = instruction_file_options[x]
	for x in command_line_options:
		if x not in def_options['ALL']:
			wfl_msg = ('CRITICAL', this_name, 'option {0} from command line is unknown'.format(x))
			print_log(wfl_msg)
		cl_tp = type(command_line_options[x])
		df_tp = type(def_options[def_options['ALL'][x][0]][def_options['ALL'][x][1]])
		if cl_tp != df_tp:
			wtp_msg = ('CRITICAL', this_name, 'in command line, option {0} is of type {1}, but it should be of type {2}'.format(x, cl_tp, df_tp))
			print_log(wtp_msg)
		else:
			if x in instruction_file_options and command_line_options[x] != instruction_file_options[x]:
				chg_msg = ('NOTICE', this_name, 'option {0} defined in instruction file as {1} was overrun by command line value {2}'.format(x, instruction_file_options[x], command_line_options[x]))
				print_log(chg_msg)
			f0 = def_options['ALL'][x][0]
			f1 = def_options['ALL'][x][1]
			if f0 not in options:
				options[f0] = collections.OrderedDict()
			options[f0][f1] = command_line_options[x]

	# Change paths if singularity container is used
	# WARNING: you still need a PPM path to retrieve the res.lib file!
	options['PATHS'][('', 'ppm_reslib_path')] = options['PATHS'][('ppm', 'ppm_path')]
	if options['PATHS'][('sing', 'singularity')] and options['PATHS'][('container', 'singularity_container')]:
		prefix = options['PATHS'][('sing', 'singularity')] + ' exec ' + options['PATHS'][('container', 'singularity_container')] + ' '
		for x, y in [('ppm', 'ppm_path'), ('muscle', 'muscle_path'), ('cesymm', 'cesymm_path'), ('symd', 'symd_path'), ('frtmalign', 'frtmalign_path'), ('dssp', 'dssp_path')]:
			sigx, sigy = "sig"+x, "sig_" + y
			options['PATHS'][(x, y)] = prefix +  options['PATHS'][(sigx, sigy)]

	if options['RUN'][('hpc', 'run_on_hpc')]:
		if not options['RUN'][('hpcdir', 'hpc_exec_path')]:
			wtp_msg = ('CRITICAL', this_name, 'run_on_hpc option is active, but no hpc exec path was inserted. Please use the option -hpcdir or --hpc_exec_path')
			print_log(wtp_msg)
		elif ((options['RUN'][('hpcdir', 'hpc_exec_path')] and options['RUN'][('hpcdir', 'hpc_exec_path')][0] != "/") and ((not options['RUN'][('hpclsf', 'local_shared_folder')]) or options['RUN'][('hpclsf', 'local_shared_folder')][0] != "/")):
			wtp_msg = ('CRITICAL', this_name, 'the option hpc_exec_path only takes absolute paths. Please insert a valid absolute path')
			print_log(wtp_msg)

	options['ALL'] = collections.OrderedDict()
	for x in options:
		if x == 'ALL':
			continue
		for y in options[x]:
			options['ALL'][y[0]] = options[x][y]
			options['ALL'][y[1]] = options[x][y]
	
	return options


def read_locations_file(locations_filename, options=collections.OrderedDict(), thr_log_status='ERROR'):
	this_name = read_locations_file.__name__

	locations = collections.OrderedDict()
	with open(locations_filename, 'r') as locations_file:
		for line in locations_file:
			if not line.strip() or line.strip().startswith('#'):
				continue
			fields = line.split(' : ')
			if len(fields) == 1:
				sector = fields[0].strip()
				locations[sector] = collections.OrderedDict()
				ancestors = ['']
			else:
				key = fields[0]
				ntabs = 0
				for c in key:
					if c == '\t':
						ntabs += 1
					else:
						break
				key = key.strip()
				if len(ancestors) == ntabs + 1:
					ancestors[-1] = key
				elif len(ancestors) == ntabs:
					ancestors.append(key)
				elif len(ancestors) > ntabs + 1:
					ancestors = ancestors[:ntabs+1]
					ancestors[-1] = key
					
				subfields = fields[1].split(' + ')
				value = ''
				if len(ancestors) > 1:
					value += locations[sector][ancestors[-2]]
				if len(subfields) > 1:
					for nsf, sf in enumerate(subfields):
						if nsf != len(subfields) - 1:
							value += locations['FSYSPATH'][sf.strip()]
						else:
							value += sf.strip()
				else:
					c = '/'
					if sector != 'EXTERNAL':
						value += fields[1].strip() + '/'
					else:
						value += fields[1].strip()
				locations[sector][key] = value

	if options:
		if not 'EXTERNAL' in locations:
			locations['EXTERNAL'] = collections.OrderedDict()
		for x in options['PATHS']:
			locations['EXTERNAL'][x[1]] = options['PATHS'][x]

	return locations


def bind_locations(locations, options, main_path="", thr_log_status='ERROR'):
    this_name = bind_locations.__name__

    if not (options['ALL']['main'] or main_path):
        mai_msg = ('CRITICAL', this_name, 'Please enter the main directory of the project')
        print_log(mai_msg)
    else:
        # If main path is relative, make it absolute
        if main_path:
            mp = main_path
        else:
            mp = options['ALL']['main']
        if mp[0] != '/':
            cwd = os.getcwd()
            if cwd[-1] != '/':
                cwd += '/'
            options['ALL']['main'] = cwd + mp
            options['GENERAL'][('main', 'main_directory')] = cwd + mp 
            print_log(('NOTICE', this_name, 'The main location is a relative path. It will be transcribed as the following absolute path: {0}'.format(options['GENERAL'][('main', 'main_directory')])))

    if 'PATHS' in options:
        if not 'EXTERNAL' in locations:
            locations['EXTERNAL'] = collections.OrderedDict()
        for x in options['PATHS']:
            locations['EXTERNAL'][x[1]] = options['PATHS'][x]

    # Write absolute paths
    for x in locations:
        for y in locations[x]:
            if '<main>' in locations[x][y]:
                if main_path:
                    locations[x][y] = locations[x][y].replace('<main>', main_path)
                else:
                    locations[x][y] = locations[x][y].replace('<main>', options['ALL']['main'])
            elif '<src>' in locations[x][y]:
                locations[x][y] = locations[x][y].replace('<src>', os.path.dirname(os.path.abspath(__file__)))

    # Check the existence of the directories
    if not main_path:
        if os.path.exists(locations['FSYSPATH']['main']):
            for y in locations['FSYSPATH']:
                if not os.path.exists(locations['FSYSPATH'][y]):
                    if options['EXECUTION'][('locp', 'location_permissive')]:
                        loc_msg = ('NOTICE', this_name, "Directory {0} is not present in this repository and will be created".format(locations['FSYSPATH'][y]))
                        #print("MKDIR", locations['FSYSPATH'][y])
                        try:
                            os.mkdir(locations['FSYSPATH'][y])
                        except:
                            print("JUMP MKDIR", locations['FSYSPATH'][y])
                    else:
                        loc_msg = ('ERROR', this_name, "Directory {0} is not present in this repository. For adding new locations on-the-run, please activate the -locp / --location_permissive flag".format(locations['FSYSPATH'][y]))
                    print_log(loc_msg)
        else:
            for y in locations['FSYSPATH']:
                if not os.path.exists(locations['FSYSPATH'][y]):
                    os.mkdir(locations['FSYSPATH'][y])

    
    # Change back to relative paths
    locations['FSYS'] = {}
    for x in locations['FSYSPATH']:
        locations['FSYS'][x] = locations['FSYSPATH'][x].replace(locations['FSYSPATH']['main'], "")
    locations['LOCSYSFILES'] = {}
    for x in locations['SYSFILES']:
        locations['LOCSYSFILES'][x] = locations['SYSFILES'][x].replace(locations['FSYSPATH']['main'], "")

    return locations

def initialize_repository(main_path="", instr_filename_path=""):
    this_name = initialize_repository.__name__

    # Read options from command line
    if not main_path:
        cl_opt = command_line_parser()
        instr_filename_path = cl_opt['instruction_filename']
    elif instr_filename_path == "":
        cl_opt = {}
        instr_filename_paths = glob.glob(main_path + "/EncoMPASS_options_relative_to*")
        if len(instr_filename_paths) > 1:
            crt_msg = ('CRITICAL', this_name, "More than one instruction file in the main path {0}.\nPlease only retain the wanted one.".format(main_path))
            print_log(crt_msg) 
        instr_filename_path = instr_filename_paths[0]
    else:
        cl_opt = {}
        instr_filename_path = main_path + "/" + instr_filename_path

    print("INSTR", instr_filename_path, "\n")

    # Read options from instruction file
    if_opt = instruction_file_parser(instr_filename_path)

    # Merge options and default values and create the "options" structure
    options = check_options(cl_opt, if_opt)
    config.log_filename = 'ENCOMPASS_DEBUG_' + options['RUN'][('code', 'main_code')] + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S") + '.txt'

    # Creates file system and 
    locations = read_locations_file(os.path.dirname(os.path.abspath(__file__)) + '/define_locations.txt', options=options)
    locations = bind_locations(locations, options, main_path=main_path)
    if not main_path:
        if not os.path.exists(options['PATHS'][('dst', 'data_structure_template')]):
            nocp_msg = ('CRITICAL', this_name, 'inexistent data structure template file: {0}'.format(options['PATHS'][('dst', 'data_structure_template')]))
            print_log(nocp_msg)
        shutil.copyfile(options['PATHS'][('dst', 'data_structure_template')], locations['SYSFILES']['data_structure_template'])
	
    # Write reference summary instruction file
    opt_s = options['GENERAL'][('s', 'instruction_filename')]
    options['GENERAL'][('', 'out_instruction_filename')] = options['ALL']['main'] + '/EncoMPASS_options_relative_to__' + os.path.basename(opt_s)

    if not main_path:
        instruction_file_writer(options, options['GENERAL'][('', 'out_instruction_filename')])

    if options['EXECUTION'][('debug', 'debug')]:
        for optcl in options:
            for k in options[optcl]:
                print(optcl, k, options[optcl][k])
        print("-----------------------------------------------------------------------\n\n")

    # Defines global variables in the config module
    config.str_data_template_fn = options['PATHS'][('dst', 'data_structure_template')]

    return options, locations

if __name__ == "__main__":
    options, locations = initialize_repository()
    print("AFTER INIT", config.log_filename)
