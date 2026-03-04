import multiprocessing
import datetime
import argparse
import shutil
import os
import sys
from pathlib import Path
from typing import Dict

from encompass.symmetry.symmetry_exec_functions import *
from encompass.symmetry.symd_default import *
from encompass.symmetry.cesymm_default import cesymm_default
from encompass.symmetry.cesymm_lower_threshold import cesymm_low_thr
from encompass.symmetry.select_pdbs_cesymm_from_symd import *
from encompass.symmetry.cesymm_from_symd import cesymm_from_symd
from encompass.symmetry.select_cs_symd_order import *
from encompass.symmetry.cesymm_order import cesymm_order
from encompass.symmetry.select_minlen_violations_pdbs import select_minlen_violations
from encompass.symmetry.cesymm_minlen_symlev import cesymm_minlen_symlev
from encompass.symmetry.select_quat_symm_minlen import select_quat_symm_minlen
from encompass.symmetry.cesymm_quat_minlen import cesymm_quat_minlen
from encompass.symmetry.cesymm_lt_rmsd import cesymm_lt_rmsd
from encompass.symmetry.ananas_run import ananas_run
from encompass.symmetry.cesymm_symd_pngs import dic_for_cesymm_symd
from encompass.symmetry.results_selection import results_selection
from encompass.symmetry.run_locusts import *


# Multiprocessing runners
def cesymm_order_mps_runner(data):
    locations, pdb_suff, flds = data
    print(("Analyzing %s with CE-Symm order" % flds[0]))
    cesymm_order_success = cesymm_order(locations, pdb_suff, flds)
    return cesymm_order_success

def mps_runner(data):
    return data[0](*data[1]) # this is, for example, equivalent to saying "return cesymm_order(locations, pdb_suff, flds)"

def mps_parallelize(np, data):

    pool = multiprocessing.Pool(processes=np)
    pool_outputs = pool.map(mps_runner, data)
    pool.close()
    pool.join()
    return pool_outputs

def cesymm_big_small_run(pdbs_to_run, fn_opt, locations_path, locations, runs_label, stepname, funcname, locusts_tmp, opt_arg = [''], str_type = 'none', run_local=False):
    """
    Cesymm_big_small_run takes a list of structures and a procedure code, divides the structures into small and big
    and submits them to locusts separately
    :param pdbs_to_run: dictionary with the structures to be submitted and each of their properties such as number of atoms or residues
    :param fn_opt: dictionary with default options for the HPC or local runs such as wall time, partition, etc.
    :param locations_path: path to the EncoMPASS database version being used
    :param locations: dictionary of locations within the EncoMPASS database
    :param runs_label: the name of the HPC folder and the label attached to each log file from the run
    :param stepname: the name of the procedure being submitted (e.g. 'cesymm_low_thr')
    :param funcname: the name of the python script that will be executed by locusts (e.g. 'cesymm_lower_threshold.py')
    :param locusts_tmp: the temporary folder on the local computer where locusts can write files during transfer
    :param opt_arg: optional argument, can be empty or carry a value such as max rmsd; default: []
    :param str_type: sets which types of structures should be submitted; options: 'all','giant','big','small','none'; default: 'none'
    :return:
    """

    # Handle local runs
    if run_local:
        command = 'source ~/venv/encompass/bin/activate; python ' + funcname + ' <arg0> <arg1> <arg2> <arg3> <arg4>'
        print("\nThe following structures will be submitted LOCALLY to CE-Symm (" + str(len(pdbs_to_run)) + "):")
        for pdb in pdbs_to_run:
            print(pdb + ' ', end='')
        print("\nRunning " + stepname + " (local) ...\n")
        run_locusts(fn_opt, locations_path, locations, runs_label, stepname, locusts_tmp, pdbs_to_run, command, opt_arg=opt_arg)
        return

    command = 'source ~/encompass/bin/activate; python ' + funcname + ' <arg0> <arg1> <arg2> <arg3> <arg4>'
    # Initiate giant structures (> 8000 residues)
    pdbs_to_run_giant = {x: pdbs_to_run[x] for x in pdbs_to_run if pdbs_to_run[x]['nres'] > 8000}
    pdbs_to_run_giant = sorted(pdbs_to_run_giant.keys(), key=lambda x: (pdbs_to_run_giant[x]['nres']), reverse=True)
    print("\nThe following GIANT structures will be submitted to CE-Symm (" + str(len(pdbs_to_run_giant)) + "):")
    for pdb in pdbs_to_run_giant:
        print(pdb + ' ', end='')
    fn_opt['waiting_time'] = '900'
    if len(pdbs_to_run_giant) > 0 and (str_type == 'all' or str_type == 'giant'):
        fn_opt['waiting_time'] = '900'
        fn_opt['partition'] = "unlimited"
        fn_opt['memory'] = '121g'
        fn_opt['constraint'] = ''
        fn_opt['cpus_per_node'] = str(min(5, len(pdbs_to_run_giant)))  # we're setting a max of 2 jobs on a node at a time due to mem. req.
        fn_opt['requested_nodes'] = max(min(10, len(pdbs_to_run_giant) // int(fn_opt['cpus_per_node'])), 1)
        print('\nrequested nodes: ', fn_opt['requested_nodes'])
        fn_opt['walltime'] = '20-00:00:00' # in theory, we can leave this empty since default is unlimited
        print("\nRunning " + stepname + " (giant proteins) ...\n")
        run_locusts(fn_opt, locations_path, locations, runs_label, stepname, locusts_tmp, pdbs_to_run_giant, command, opt_arg=opt_arg)


    # Initiate big structures
    pdbs_to_run_big = {x: pdbs_to_run[x] for x in pdbs_to_run if pdbs_to_run[x]['natoms'] > 20000 and pdbs_to_run[x]['nres'] < 8001}
    #pdbs_to_run_big = {x: pdbs_to_run[x] for x in pdbs_to_run if pdbs_to_run[x]['natoms'] > 20000}
    pdbs_to_run_big = sorted(pdbs_to_run_big.keys(), key=lambda x: (pdbs_to_run_big[x]['natoms']), reverse=True)
    print("\nThe following BIG structures will be submitted to CE-Symm (" + str(len(pdbs_to_run_big)) + "):")
    for pdb in pdbs_to_run_big:
        print(pdb + ' ', end='')
    if len(pdbs_to_run_big) > 0 and (str_type == 'all' or str_type == 'big'):
        fn_opt['partition'] = 'norm'
        fn_opt['cpus_per_node'] = str(min(7, len(pdbs_to_run_big)))  # we're setting a max of 7 jobs on a node at a time due to mem. req. (limited to max 30GB per job)
        fn_opt['requested_nodes'] = max(min(50, len(pdbs_to_run_big) // int(fn_opt['cpus_per_node'])), 1)
        print('\nrequested nodes: ', fn_opt['requested_nodes'])
        fn_opt['walltime'] = '10-00:00:00'
        print("\nRunning " + stepname + " (big proteins) ...\n")
        run_locusts(fn_opt, locations_path, locations, runs_label, stepname, locusts_tmp, pdbs_to_run_big, command, opt_arg=opt_arg)

    # Initiate small structures
    pdbs_to_run_small = {x: pdbs_to_run[x] for x in pdbs_to_run if pdbs_to_run[x]['natoms'] <= 20000}
    pdbs_to_run_small = sorted(pdbs_to_run_small.keys(), key=lambda x: (pdbs_to_run_small[x]['natoms']), reverse=True)
    print("\nThe following SMALL structures will be submitted to CE-Symm (" + str(len(pdbs_to_run_small)) + "):")
    for pdb in pdbs_to_run_small:
        print(pdb + ' ', end='')
    if len(pdbs_to_run_small) > 0 and (str_type == 'all' or str_type == 'small'):
        fn_opt['partition'] = 'norm'
        fn_opt['cpus_per_node'] = str(min(28, len(pdbs_to_run_small)))
        fn_opt['requested_nodes'] = max(min(50, len(pdbs_to_run_small) // int(fn_opt['cpus_per_node'])), 1)
        print('\nrequested nodes: ', fn_opt['requested_nodes'])
        fn_opt['walltime'] = '10-00:00:00'
        print("\nRunning " + stepname + " (small proteins) ...\n")
        run_locusts(fn_opt, locations_path, locations, runs_label, stepname, locusts_tmp, pdbs_to_run_small, command, opt_arg=opt_arg)

    return

def clean_sym_results(locations, wipe_list, procedure=None):
    """
    Remove all log files and symmetry analysis files for a given list of pdbs. If no specific symmetry procedure is specified,
    all will be removed.

    :param locations: the loaded dictionary of database paths
    :param wipe_list: list of structures to remove
    :param procedure: symmetry analysis step, e.g. "cesymm"
    :return:
    """
    for pdb in wipe_list:
        if procedure:
            proc_names = [procedure]
        else:
            proc_names = ['symd','cesymm','cesymm_low_thr','cesymm_from_symd', 'cesymm_order', 'cesymm_minlen',
                          'cesymm_quat_minlen','cesymm_rmsd_2_5', 'cesymm_rmsd_3', 'ananas', 'quatsymm',
                          'images_cesymm_symd', 'mssd', 'transfer']
        for proc_name in proc_names:
            logs = glob.glob(locations['FSYSPATH'][proc_name + '_log'] + pdb + "*")
            if proc_name in locations['FSYSPATH']:
                syms = glob.glob(locations['FSYSPATH'][proc_name] + pdb + "*")
            files = logs + syms
            for f in files:
                if os.path.isdir(f):
                    print("Removing %s" % f)
                    shutil.rmtree(f)
                elif os.path.isfile(f):
                    print("Removing %s" % f)
                    os.remove(f)
                else:
                    print("%s was not removed because it is neither file nor folder." % f)
            if len(files) == 0:
                print("Nothing removed for %s in %s." % (pdb, proc_name))
        print("\n")
    return



def symm_update(locations_path, wkdir, input_file, pdb_suff, str_type, runs_label, locusts_tmp, step_name='test', instr_path='', generate=True, clean=None, run_local=False):
    from encompass.pipeline.initialize_repository import initialize_repository
    options, locations = initialize_repository(main_path=locations_path, instr_filename_path=instr_path)
    #print(locations)
    encompass_archive_name = 'str_data_completegen.pkl'
    # encompass_archive_name = 'str_data_straln.pkl'

    #Read from a checkpoint example
    #from supporting_functions import read_checkpoint
    #str_data = read_checkpoint('/data/local/encompass/dev/toni_full_dec2020_beforeFrTM/logs/cache/str_data_completegen.pkl',locations['SYSFILES']['data_structure_template'])
    #print(str_data.keys())
    #print(str_data['1okc'])
    #return str_data

    #options, filters, locations = filesystem_info(locations_path)
    # Make sure tm_archive is created
    #create_tm_archive(locations, wkdir)

    if clean:
        if os.path.isfile(clean):
            clean_list = []
            with open(clean, 'r') as f:
                while (line := f.readline().rstrip()):
                    clean_list.append(line)
            print("All results will be cleaned for : ", clean_list)
            clean_sym_results(locations, clean_list, procedure='')
        else:
            raise SystemExit("% does not exist.\n" % clean)
        return

    if input_file:
        pdb_list = set()
        with open(input_file, 'r') as lst:
            for line in lst:
                if not line.startswith('#'):
                    pdb_list.add(line[0:4])

    if step_name == 'test':
        print("Test will run in ", locations['FSYSPATH']['symmtemp'])
        create_tm_archive_new(locations, locations['FSYSPATH']['symmtemp'], encompass_archive_name)
        return "Test completed."
    print('Will run the following procedures: ', step_name)
    if not os.path.isfile(locations['SYSFILES']['tm_archive']):
        create_tm_archive_new(locations, locations['FSYSPATH']['symmtemp'], encompass_archive_name)
    print('Using TM archive ' + locations['SYSFILES']['tm_archive'])
    tm_archive = pkl.load(open(locations['SYSFILES']['tm_archive'], 'rb'))
    # Remove structures that should be ignored because of 1) no TMs; 2) no residues; 3) UNK residues
    ignore = [pdb_key for pdb_key in tm_archive.keys() if (len(tm_archive[pdb_key]['tmchains']) == 0 or 'nres' not in tm_archive[pdb_key])] #or tm_archive[pdb_key]['unk']==True)]
    for pdb_key in ignore:
        print(f"Removing {pdb_key} from tm_archive.")
        tm_archive.pop(pdb_key)
        if input_file and pdb_key in pdb_list:
            pdb_list.remove(pdb_key)

    if not input_file:
        pdb_list = set(tm_archive.keys())
#    pdb_list = {'1a0s', '1a0t', '1a91', '1af6', '1afo', '1aig', '1aij', '1ap9', '1ar1', '1at9', '1aty', '1ay2', '1bcc',
#                '1be3', '1bgy', '1bh3', '1bl8', '1bm1', '1brd', '1brr', '1brx', '1bt9', '1bxw', '1by3', '1by5', '1c0v',
#                '1c17', '1c3w', '1c8r', '1c8s', '1c99', '1cwq', '1ds8', '1dv3', '1dv6', '1dx7', '1dxr', '1dze', '1dzo',
#                '1e0p', '1e12', '1e14', '1e54', '1e6d', '1e7p', '1ehk', '1ek9', '1eys', '1ezv', '1f4z', '1f50', '1f6g',
#                '1f6n', '1f88', '1fbb', '1fbk', '1fcp', '1fdm', '1fe1', '1fep'}

    print("Symmetry analysis will proceed with %d proteins." % len(pdb_list))

    #np = int(options['number_of_procs'])
    np = multiprocessing.cpu_count()
    mem_req_limit = 20000
    np_high_mem = 3
    fail_flag = 0
    print("Number of processors in local machine: ", np)

    # Set default locusts options
    fn_opt = default_fn_options(options)
    fn_opt['email_address'] = ""
    if run_local:
        fn_opt['run_on_hpc'] = False
        fn_opt['hpc_singularity'] = ""
        #fn_opt['singularity'] = ""
        fn_opt['singularity'] = "singularity"
        #fn_opt['singularity_modload'] = ""
        #fn_opt['singularity_modload'] = "module load singularity"
        print(f"Local singularity contaner to be used: {fn_opt['singularity_container']}")
        fn_opt['extra_outer_statements'] = 'export SINGULARITY_BINDPATH="$SINGULARITY_BINDPATH,/data/TMB-CSB/,/home/aleksandrovaaa/,/data/local/,/usr/local/"'
        # for local runs, use: 'export SINGULARITY_BINDPATH="$SINGULARITY_BINDPATH,/data/TMB-CSB/,/home/aleksandrovaaa/,/data/local/,/usr/local/"'
        fn_opt['python'] = "/data/local/encompass/venv_encompass_3.12/bin/python"
    else:
        print(f"HPC singularity container will be used: {fn_opt['hpc_singularity_container']}")
        fn_opt['python'] = options['hpcactipy'].replace("activate", "python")
        print(f"Using HPC python environment {fn_opt['python']}")
    if os.path.isdir(locusts_tmp)==False and str_type!="none":
        print("Creating folder " + locusts_tmp)
        os.mkdir(locusts_tmp)

    # Step 1
    if step_name == 'symd' or step_name == 'all': #or step_name == 'cesymm_from_symd' or step_name == 'all':

        symd_completed, failed = check_completion(locations, pdb_list, 'symd')
        pdbs_to_run = {x: tm_archive[x] for x in pdb_list.difference(symd_completed)}
        pdbs_to_run = sorted(pdbs_to_run.keys(), key=lambda x: (pdbs_to_run[x]['natoms']), reverse=True)
        print("The following structures will be submitted to SymD (" + str(len(pdbs_to_run)) + "):")
        for pdb in pdbs_to_run:
            print(pdb + ' ', end = '')

        if str_type != "none":
            if not os.path.isdir(locations['FSYSPATH']['symmtemp']+"symd/"):
                os.mkdir(locations['FSYSPATH']['symmtemp']+"symd/")
            if not os.path.isdir(locations['FSYSPATH']['symmtemp']+"symd/temp_structures"):
                os.mkdir(locations['FSYSPATH']['symmtemp']+"symd/temp_structures")
            if run_local:
                command = f"{fn_opt['python']} symd_default.py <arg0> <arg1> <arg2> <arg3>"
            else:
                fn_opt['requested_nodes'] = '40'
                fn_opt['walltime'] = '10-00:00:00'
                fn_opt['waiting_time'] = '900'
                command = "{fn_opt['python']} symd_default.py <arg0> <arg1> <arg2> <arg3>"
            if len(pdbs_to_run) > 0:
                print("\nRunning SymD ...")
                run_locusts(fn_opt, locations_path, locations, runs_label, 'symd', locusts_tmp, pdbs_to_run, command)
            else:
                print("\nNothing to submit to SymD.")

        result, negresult = check_completion(locations, pdb_list, 'symd')
        if len(pdb_list.difference(symd_completed).difference(result)) > 0:
            fail_flag = 1
        print("\nSymD completed for " + str(len(result)) + " structures. SymD failed for the following structures (" + str(len(negresult)) + "):")
        for pdb in negresult:
            print(pdb + ' ', end = '')

        symd_successful = symd_completed.union(result)

    # Step 2
    if step_name == 'cesymm' or step_name == 'all': #or step_name == 'cesymm_from_symd' or step_name == 'all':
        if not os.path.isdir(locations['FSYSPATH']['symmtemp']+"cesymm_default/"):
            os.mkdir(locations['FSYSPATH']['symmtemp']+"cesymm_default/")

        cesymm_completed, failed = check_completion(locations, pdb_list, 'cesymm')

        if step_name == 'cesymm':
            pdbs_to_run = {x: tm_archive[x] for x in pdb_list.difference(cesymm_completed)}
        else:
            pdbs_to_run = {x: tm_archive[x] for x in symd_successful.difference(cesymm_completed)}

        # Use local hpc exec to run cesymm
        #fn_opt['hpc_singularity'] = ""
        #fn_opt['singularity'] = ""
        if str_type != "none" and len(pdbs_to_run)>0:
            cesymm_big_small_run(pdbs_to_run, fn_opt, locations_path, locations, runs_label, 'cesymm', 'cesymm_default.py',
                                locusts_tmp, str_type=str_type, run_local=run_local)


        result, negresult = check_completion(locations, pdb_list, 'cesymm')
        if len(pdb_list.difference(cesymm_completed).difference(result)) > 0:
            fail_flag = 1
        print("\nCE-Symm completed for " + str(len(result)) + " structures. CE-Symm failed for the following structures (" + str(len(negresult)) + "):")
        for pdb in negresult:
            print(pdb + ' ', end = '')



    # Step 3
    if step_name == 'cesymm_low_thr' or step_name == 'all':
        if not os.path.isdir(locations['FSYSPATH']['symmtemp']+"cesymm_low_thr/"):
            os.mkdir(locations['FSYSPATH']['symmtemp']+"cesymm_low_thr/")
        if not os.path.isdir(locations['FSYSPATH']['symmtemp']+"cesymm_low_thr/temp_structures"):
            os.mkdir(locations['FSYSPATH']['symmtemp']+"cesymm_low_thr/temp_structures")

        """
        struct = "6cp6"
        cesymm_low_thr(locations, options, struct, pdb_suff, mem_req_limit = 20000, run_id = runs_label)
        with open(locations['FSYSPATH']['cesymm_low_thr_log'] + struct + 'cesymm_lt.log', 'r') as f:
            print(f.read())
        """

        completed, failed = check_completion(locations, pdb_list, 'cesymm_low_thr')

        if step_name == 'cesymm_low_thr':
            pdbs_to_run = {x: tm_archive[x] for x in pdb_list.difference(completed)}
        else:
            pdbs_to_run = {x: tm_archive[x] for x in symd_successful.difference(completed)}


        if str_type != "none" and len(pdbs_to_run) > 0:
            cesymm_big_small_run(pdbs_to_run, fn_opt, locations_path, locations, runs_label, 'cesymm_low_thr', 'cesymm_lower_threshold.py',
                                 locusts_tmp, str_type = str_type, run_local=run_local)

        result, negresult = check_completion(locations, pdb_list, 'cesymm_low_thr')
        if len(pdb_list.difference(completed).difference(result)) > 0:
            fail_flag = 1
        print("CE-Symm Lower Threshold failed for the following structures (" + str(len(negresult)) + "):")
        for pdb in negresult:
            print(pdb + ' ', end = '')



    # Step 4
    ## Runs for ~ 1:40 hr
    if step_name == 'cesymm_from_symd' or step_name == 'all':
        if not os.path.isdir(locations['FSYSPATH']['symmtemp'] + "cesymm_from_symd/"):
            os.mkdir(locations['FSYSPATH']['symmtemp'] + "cesymm_from_symd/")

        """
        #struct = "5d92_C_-43-204"
        struct = "5twv_F_1163-1577"
        cesymm_from_symd(locations, options, struct, pdb_suff, mem_req_limit = 20000, run_id = runs_label)
        with open(locations['FSYSPATH']['cesymm_from_symd_log'] + struct + '_fragment.log', 'r') as f:
            print(f.read())
        """

        if generate == True:
            cesymm_from_symd_set = select_cesymm_from_symd(locations, pdb_suff)
        else:
            cesymm_from_symd_set = pkl.load(open(locations['FSYSPATH']['cesymm_from_symd_log'] + 'fragments_to_run.pkl', 'rb'))
        cesymm_from_symd_set_run = {x for x in cesymm_from_symd_set if x[0:4] in pdb_list}
        completed, failed = check_completion(locations, cesymm_from_symd_set_run, 'cesymm_from_symd')
        pdbs_to_run = {x: tm_archive[x[0:4]][x[5:6]] for x in cesymm_from_symd_set_run.difference(completed)}

        if len(pdbs_to_run) > 0 and str_type != 'none':
            cesymm_big_small_run(pdbs_to_run, fn_opt, locations_path, locations, runs_label, 'cesymm_from_symd', 'cesymm_from_symd.py',
                                locusts_tmp, str_type = str_type, run_local=run_local)

        print("\n\nDONE\n\n")
        result, negresult = check_completion(locations, cesymm_from_symd_set_run, 'cesymm_from_symd')
        if len(cesymm_from_symd_set_run.difference(completed).difference(result)) > 0:
            fail_flag = 1
        print("\nCE-Symm from SymD completed for " + str(len(result)) + " structures. CE-Symm from SymD failed for the following structures (" + str(len(negresult)) + "):")
        for pdb in negresult:
            print(pdb + ' ', end = '')


    # Step 5
    if step_name == 'cesymm_order' or step_name == 'all':
        #Big structures run took ~21 hrs, small structures took
        if not os.path.isdir(locations['FSYSPATH']['symmtemp'] + "cesymm_order/"):
            os.mkdir(locations['FSYSPATH']['symmtemp'] + "cesymm_order/")
        if generate == True:
            cs_symd_order_dic = select_cs_symd_order(locations)
        else:
            cs_symd_order_dic = pkl.load(open(locations['FSYSPATH']['cesymm_order_log'] + 'cs_symd_order.pkl','rb'))

        cs_symd_order_run = {x for x in cs_symd_order_dic if x[0:4] in pdb_list}
        print(cs_symd_order_dic)
        completed, failed = check_completion(locations, cs_symd_order_run, 'cesymm_order')
        print("Completed", completed, "Failed", failed)
        """
        cesymm_order(locations, options, '6dra_df', pdb_suff, mem_req_limit = 20000, run_id = runs_label)
        with open(locations['FSYSPATH']['cesymm_order_log'] + '6dra_df' + '_ord.log', 'r') as f:
            print(f.read())
        """
        pdbs_to_run = {x: cs_symd_order_dic[x] for x in cs_symd_order_run.difference(completed)}
        print("Will run:", pdbs_to_run)
        if str_type != 'none' and len(pdbs_to_run) > 0:
            cesymm_big_small_run(pdbs_to_run, fn_opt, locations_path, locations, runs_label, 'cesymm_order', 'cesymm_order.py',
                                locusts_tmp, str_type = str_type, run_local=run_local)

        result, negresult = check_completion(locations, cs_symd_order_run, 'cesymm_from_symd')
        if len(cs_symd_order_run.difference(completed).difference(result)) > 0:
            fail_flag = 1
        print("CE-Symm with enforced order failed for the following structures (" + str(len(negresult)) + "):")
        for pdb in negresult:
            print(pdb + ' ', end = '')



    
    # Step 6
    if step_name == 'cesymm_minlen' or step_name == 'all':
        if not os.path.isdir(locations['FSYSPATH']['symmtemp'] + "cesymm_minlen/"):
            os.mkdir(locations['FSYSPATH']['symmtemp'] + "cesymm_minlen/")
        if generate == True:
            minlen_dic = select_minlen_violations(locations, pdb_suff, tm_archive)
        else:
            minlen_dic = pkl.load(open(locations['FSYSPATH']['cesymm_minlen_log'] + 'minlen_violations_pdbs.pkl','rb'))

        minlen_run = {x for x in minlen_dic if x[0:4] in pdb_list}
        completed, failed = check_completion(locations, minlen_run, 'cesymm_minlen')
        print(completed, failed, set(minlen_dic.keys()))
        #'3v3c_lt', '6nfv_ordered_lt', '3m75_A_df'

        #cesymm_minlen_symlev(locations, options, '3v3c_lt', pdb_suff, mem_req_limit = 20000, run_id = runs_label)
        #with open(locations['FSYSPATH']['cesymm_minlen_log'] + '3v3c_lt' + '_minlen.log', 'r') as f:
        #    print(f.read())

        pdbs_to_run = {x: minlen_dic[x] for x in minlen_run.difference(completed)}

        if str_type != 'none' and len(pdbs_to_run) > 0:
            cesymm_big_small_run(pdbs_to_run, fn_opt, locations_path, locations, runs_label, 'cesymm_minlen', 'cesymm_minlen_symlev.py',
                                 locusts_tmp, str_type = str_type, run_local=run_local)

        result, negresult = check_completion(locations, minlen_run, 'cesymm_minlen')
        if len(minlen_run.difference(completed).difference(result)) > 0:
            fail_flag = 1
        print("CE-Symm minlen failed for the following structures (" + str(len(negresult)) + "):")
        for pdb in negresult:
            print(pdb + ' ', end = '')


    # Step 7
    if step_name == 'cesymm_quat_minlen' or step_name == 'all':
        # Generation of selected structurs took ~44 hrs.
        # Big structures took ~27 hrs and small structures took ~5 hrs (total ~33 hrs)
        if not os.path.isdir(locations['FSYSPATH']['symmtemp'] + "cesymm_quat_minlen/"):
            os.mkdir(locations['FSYSPATH']['symmtemp'] + "cesymm_quat_minlen/")
        if generate == True:
            quat_minlen_dic = select_quat_symm_minlen(locations, pdb_suff, tm_archive)
        else:
            quat_minlen_dic = pkl.load(open(locations['FSYSPATH']['cesymm_quat_minlen_log'] + 'quat_symm_minlen_pdbs.pkl','rb'))
        quat_minlen_run = {x for x in quat_minlen_dic if x[0:4] in pdb_list}
        completed, failed = check_completion(locations, quat_minlen_run, 'cesymm_quat_minlen')
        print(completed, failed, set(quat_minlen_dic.keys()))

        """
        cesymm_quat_minlen_symlev(locations, options, '6drc_df', pdb_suff, mem_req_limit = 20000, run_id = runs_label)
        with open(locations['FSYSPATH']['cesymm_quat_minlen_log'] + '6drc_df' + '_quat.log', 'r') as f:
	        print(f.read())
        """
        pdbs_to_run = {x: quat_minlen_dic[x] for x in quat_minlen_run.difference(completed)}
        if str_type != 'none' and len(pdbs_to_run) > 0:
            cesymm_big_small_run(pdbs_to_run, fn_opt, locations_path, locations, runs_label, 'cesymm_quat_minlen', 'cesymm_quat_minlen.py',
                             locusts_tmp, str_type = str_type, run_local=run_local)

        result, negresult = check_completion(locations, quat_minlen_run, 'cesymm_quat_minlen')
        if len(quat_minlen_run.difference(completed).difference(result)) > 0:
            fail_flag = 1
        print("CE-Symm Quat minlen failed for the following structures (" + str(len(negresult)) + "):")
        for pdb in negresult:
            print(pdb + ' ', end = '')


    # Step 8
    if step_name == 'cesymm_rmsd_2_5' or step_name == 'all':
        maxrmsd = 2.5
        if not os.path.isdir(locations['FSYSPATH']['symmtemp'] + "cesymm_rmsd_2_5/"):
            os.mkdir(locations['FSYSPATH']['symmtemp'] + "cesymm_rmsd_2_5/")

        """
        struct = "6cp6"
        cesymm_lt_rmsd(locations, options, struct, pdb_suff, maxrmsd, run_id = runs_label)
        with open(locations['FSYSPATH']['cesymm_rmsd_2_5_log'] + struct + '_rmsd_2_5.log', 'r') as f:
            print(f.read())
        """

        completed, failed = check_completion(locations, pdb_list, 'cesymm_rmsd_2_5')
        pdbs_to_run = {x: tm_archive[x] for x in pdb_list.difference(completed)}
        #print(pdbs_to_run['6l7o']['nres'])
        # Assign the number of residues of the whole structure to be
        # the number of residues of the largest chain since this procedure
        # only runs on chains
        for x in pdbs_to_run:
            max_ch, max_natoms = '', 0
            for ch in pdbs_to_run[x]['tmchains']:
                if pdbs_to_run[x][ch]['natoms'] > max_natoms:
                    max_ch = ch
                    max_natoms = pdbs_to_run[x][ch]['natoms']
            pdbs_to_run[x]['nres'] = pdbs_to_run[x][max_ch]['nres']
            pdbs_to_run[x]['natoms'] = pdbs_to_run[x][max_ch]['natoms']

        if str_type != 'none' and len(pdbs_to_run) > 0:
            cesymm_big_small_run(pdbs_to_run, fn_opt, locations_path, locations, runs_label, 'cesymm_rmsd_2_5', 'cesymm_lt_rmsd.py',
                                 locusts_tmp, opt_arg = ['2.5'], str_type = str_type, run_local=run_local)

        result, negresult = check_completion(locations, pdb_list, 'cesymm_rmsd_2_5')
        if len(pdb_list.difference(completed).difference(result)) > 0:
            fail_flag = 1
        print("CE-Symm Lower Threshold RMSD 2.5 failed for the following structures (" + str(len(negresult)) + "):")
        for pdb in negresult:
            print(pdb + ' ', end = '')


    # Step 9
    if step_name == 'cesymm_rmsd_3' or step_name == 'all':
        # Took 22.7 hrs
        maxrmsd = 3
        if not os.path.isdir(locations['FSYSPATH']['symmtemp'] + "cesymm_rmsd_3/"):
            os.mkdir(locations['FSYSPATH']['symmtemp'] + "cesymm_rmsd_3/")

        """
        struct = "6cp6"
        cesymm_lt_rmsd(locations, options, struct, pdb_suff, maxrmsd, run_id = runs_label)
        with open(locations['FSYSPATH']['cesymm_rmsd_3_log'] + struct + '_rmsd_3.log', 'r') as f:
            print(f.read())
        """

        completed, failed = check_completion(locations, pdb_list, 'cesymm_rmsd_3')
        pdbs_to_run = {x: tm_archive[x] for x in pdb_list.difference(completed)}
        #print(pdbs_to_run['6l7o']['nres'])
        # Assign the number of residues of the whole structure to be
        # the number of residues of the largest chain since this procedure
        # only runs on chains
        for x in pdbs_to_run:
            max_ch, max_natoms = '', 0
            for ch in pdbs_to_run[x]['tmchains']:
                if pdbs_to_run[x][ch]['natoms'] > max_natoms:
                    max_ch = ch
                    max_natoms = pdbs_to_run[x][ch]['natoms']
            pdbs_to_run[x]['nres'] = pdbs_to_run[x][max_ch]['nres']
            pdbs_to_run[x]['natoms'] = pdbs_to_run[x][max_ch]['natoms']

        if str_type != 'none' and len(pdbs_to_run) > 0:
            cesymm_big_small_run(pdbs_to_run, fn_opt, locations_path, locations, runs_label, 'cesymm_rmsd_3', 'cesymm_lt_rmsd.py',
                                 locusts_tmp, opt_arg = ['3'], str_type = str_type, run_local=run_local)

        result, negresult = check_completion(locations, pdb_list, 'cesymm_rmsd_3')
        if len(pdb_list.difference(completed).difference(result)) > 0:
            fail_flag = 1
        print("CE-Symm Lower Threshold RMSD 3 failed for the following structures (" + str(len(negresult)) + "):")
        for pdb in negresult:
            print(pdb + ' ', end = '')


    # Step 10
    if step_name == 'quatsymm' or step_name == 'all':
        if os.path.isdir(locations['FSYSPATH']['symmtemp'] + "quatsymm/") == False:
            os.mkdir(locations['FSYSPATH']['symmtemp'] + "quatsymm/")
        completed, failed = check_completion(locations, pdb_list, 'quatsymm')
        print("Completed structures before run: ", completed, "\nFailed structures before run: ", failed)
        pdbs_to_run = pdb_list.difference(completed)

        if run_local and len(pdbs_to_run) > 0:
            fn_opt['singularity'] = ""
            fn_opt['singularity_container'] = ""
            #fn_opt['extra_outer_statements'] = ""
            fn_opt['extra_outer_statements'] = "module load pymol"
            command = "source ~/venv/encompass/bin/activate; python quatsymm_complete.py <arg0> <arg1> <arg2> <arg3>"
            run_locusts(fn_opt, locations_path, locations, runs_label, 'quatsymm', locusts_tmp, pdb_list, command, opt_arg = [''], gather = False)
        elif len(pdbs_to_run) > 0:
            fn_opt['partition'] = 'quick,norm'
            fn_opt['memory'] = '58g'
            fn_opt['constraint'] = ''
            fn_opt['cpus_per_node'] = 1
            fn_opt['requested_nodes'] = '40'
            fn_opt['walltime'] = '0-03:30:00'
            fn_opt['waiting_time'] = '600'
            fn_opt['extra_outer_statements'] = fn_opt['extra_outer_statements'] + "; module load pymol/2.5.0"
            command = "source ~/encompass/bin/activate; python quatsymm_complete.py <arg0> <arg1> <arg2> <arg3>"
            run_locusts(fn_opt, locations_path, locations, runs_label, 'quatsymm', locusts_tmp, pdb_list, command, opt_arg = [''], gather = False)


        result, negresult = check_completion(locations, pdbs_to_run, 'quatsymm')
        if len(pdbs_to_run.difference(completed).difference(result)) > 0:
            fail_flag = 1
        print("QuatSymm failed for the following structures (" + str(len(negresult)) + "):")
        for pdb in negresult:
            print(pdb + ' ', end = '')


    # Step X
    if step_name == 'ananas' or step_name == "all":
        if os.path.isdir(locations['FSYSPATH']['symmtemp'] + "ananas/") == False:
            os.mkdir(locations['FSYSPATH']['symmtemp'] + "ananas/")
        ananas_completed, failed = check_completion(locations, pdb_list, 'ananas')
        data = []
        print("The following structures will be submitted to AnAnaS:")
        for pdb in pdb_list.difference(ananas_completed):
            data.append((ananas_run, (locations, pdb, pdb_suff)))
            print(pdb + ' ', end='')

        if len(data) > 0:
            print("\nRunning AnAnaS ...")
            result = mps_parallelize(np, data)
            print('AnAnaS completed: ', result)
    #    if len(pdb_list.difference(ananas_completed).difference(result)) > 0:
    #        fail_flag = 1

    # Step 11: Data & images for CE-Symm and SymD
    # Takes ~ 27 hrs
    if step_name == "images_cesymm_symd" or step_name == "all":
        if os.path.isdir(locations['FSYSPATH']['symmtemp'] + "webarchive/") == False:
            os.mkdir(locations['FSYSPATH']['symmtemp'] + "webarchive/")
        completed, failed = check_completion(locations, pdb_list, 'images_cesymm_symd')
        print(completed, failed)
        pdbs_to_run = pdb_list.difference(completed)


        fn_opt['singularity'] = ""
        fn_opt['singularity_container'] = ""
        fn_opt['extra_outer_statements'] = ""
        if run_local and len(pdbs_to_run) > 0:
            command = "source ~/venv/encompass/bin/activate; python cesymm_symd_pngs.py <arg0> <arg1> <arg2> <arg3>"
            run_locusts(fn_opt, locations_path, locations, runs_label, 'images_cesymm_symd', locusts_tmp, pdb_list, command, opt_arg = [''], gather = False)
        elif len(pdbs_to_run) > 0:
            fn_opt['partition'] = 'quick,norm'
            fn_opt['memory'] = '58g'
            fn_opt['constraint'] = ''
            fn_opt['cpus_per_node'] = 1
            fn_opt['requested_nodes'] = '40'
            fn_opt['walltime'] = '0-03:30:00'
            fn_opt['waiting_time'] = '600'
            fn_opt['hpc_singularity'] = ""
            fn_opt['singularity'] = ""
            fn_opt['extra_outer_statements'] = "module load pymol/2.5.0"
            command = "source ~/encompass/bin/activate; python cesymm_symd_pngs.py <arg0> <arg1> <arg2> <arg3>"
            run_locusts(fn_opt, locations_path, locations, runs_label, 'images_cesymm_symd', locusts_tmp, pdb_list, command, opt_arg = [''], gather = False)
        print('The CE-Symm and SymD images were completed.')

    # Step 12: MSSD results selection
    if step_name == 'mssd' or step_name == 'all':
        if os.path.isdir(locations['FSYSPATH']['symmtemp'] + "webarchive/") == False:
            os.mkdir(locations['FSYSPATH']['symmtemp'] + "webarchive/")
        if not os.path.isdir(locations['FSYSPATH']['symmtemp']  + "webarchive/mssd/"):
            os.mkdir(locations['FSYSPATH']['symmtemp']  + "webarchive/mssd/")

        wkdir = locations['FSYSPATH']['symmtemp']  + "selected/"
        if os.path.isdir(wkdir) == False:
            os.mkdir(wkdir)

        ### Clean up before running
        if fail_flag == 0:
            outs_list = [locations['SYSFILES']['mssd_symmetries_whole'],
                         locations['SYSFILES']['mssd_symmetries_chains']]

            for out in outs_list:
                if os.path.isfile(out):
                    olds_wkdir = wkdir + 'old/'
                    if os.path.isdir(olds_wkdir) == False:
                        os.mkdir(olds_wkdir)
                    t = os.path.getmtime(out)
                    d = str(datetime.date.fromtimestamp(t))
                    out_name = out.split('/')[-1]
                    os.rename(out, olds_wkdir + out_name[:-4] + '_' + d + out_name[-4:])

            selected_out = locations['FSYSPATH']['mssd']
            if os.path.isdir(selected_out) == False:
                os.mkdir(selected_out)
            completed, failed = check_completion(locations, pdb_list, 'mssd')
            print("Completed: ", completed, "\n\nFailed: ", failed)
            pdbs_to_run = pdb_list.difference(completed)

            fn_opt['singularity'] = ""
            fn_opt['singularity_container'] = ""
            fn_opt['extra_outer_statements'] = ""
            if run_local:
                command = "source ~/venv/encompass/bin/activate; python results_selection.py <arg0> <arg1> <arg2> <arg3> <arg4>"
                run_locusts(fn_opt, locations_path, locations, runs_label, 'mssd', locusts_tmp, pdb_list, command, opt_arg = [''], gather = False)
            else:
                fn_opt['partition'] = 'norm'
                fn_opt['memory'] = '58g'
                fn_opt['constraint'] = ''
                fn_opt['cpus_per_node'] = 1
                fn_opt['requested_nodes'] = '60'
                fn_opt['walltime'] = '0-08:45:00'
                fn_opt['waiting_time'] = '600'
                fn_opt['hpc_singularity'] = ""
                fn_opt['singularity'] = ""
                fn_opt['extra_outer_statements'] = "module load pymol/2.5.0; module load muscle/5.1"
                command = "source ~/encompass/bin/activate; python results_selection.py <arg0> <arg1> <arg2> <arg3> <arg4>"
                run_locusts(fn_opt, locations_path, locations, runs_label, 'mssd', locusts_tmp, pdb_list, command, opt_arg = [''], gather = False)
            print('Result selection completed.')


            """
            symmdic = {}
            for count, inputf in enumerate(sorted(list(pdbs_to_run))):
                print("#%d of %d proteins (%d perc., pdb %s)" % (count, len(pdbs_to_run), 100*count//len(pdbs_to_run), inputf))
                symmdic.update(results_selection(locations, options, pdb_suff, inputf))
            print('Result selection completed.')
            pkl.dump(symmdic, open(wkdir + 'symmetry_results_' + runs_label + '.pkl','wb'))
            print('Results file written and dictionary saved at ' + wkdir + 'symmetry_results_' + runs_label + '.pkl.')
            """

        else:
            print('Failed runs - not ready for result selection.')

    # Step 13: Combine MSSD results
    if step_name == "mssd_combine_dic" or step_name == "all":
        completed, failed = check_completion(locations, pdb_list, 'mssd')
        print("Completed (", len(completed), ") :", completed, "\n\nFailed (", len(failed), "): ", failed)

        results = {}
        for inputf in completed:
            print(inputf)
            print(locations['FSYSPATH']['mssd_log']  + inputf + "/" + inputf + "_dic.pkl")
            if os.path.isfile(locations['FSYSPATH']['mssd_log']  + inputf + "/" + inputf + "_dic.pkl"):
                d = pkl.load(open(locations['FSYSPATH']['mssd_log']  + inputf + "/" + inputf + "_dic.pkl","rb"))
                #for struct in d.keys():
                #    d[struct]['name'] = tm_archive[struct[0:4]]['name']
                results.update(d)

        #t = time.localtime()
        #date = "%s-%s-%s" % (t.tm_mon, t.tm_mday, t.tm_year)

        pkl.dump(results, open(locations['FSYSPATH']['mssd'] + "symmetry_results.pkl","wb"))

    # Step 14: Finding the sequence and structural neighbors of each structure
    if step_name == 'neighbors' or step_name == 'all':
        if not os.path.isdir(locations['FSYSPATH']['neighbors']):
            os.mkdir(locations['FSYSPATH']['neighbors'])
        ch_list = [pdb + "_" + ch for pdb in tm_archive.keys() for ch in tm_archive[pdb]['tmchains']]

        completed, failed = check_completion(locations, ch_list, 'neighbors')
        print("Completed before start: ", completed, "\n\nFailed before start: ", failed)
        #'3v3c_lt', '6nfv_ordered_lt', '3m75_A_df'

        #cesymm_minlen_symlev(locations, options, '3v3c_lt', pdb_suff, mem_req_limit = 20000, run_id = runs_label)
        #with open(locations['FSYSPATH']['cesymm_minlen_log'] + '3v3c_lt' + '_minlen.log', 'r') as f:
        #    print(f.read())

        chains_to_run = {x for x in set(ch_list).difference(completed)}
        # Use local hpc exec to run cesymm
        fn_opt['partition'] = 'quick,norm'
        fn_opt['cpus_per_node'] = str(min(28, len(chains_to_run)))
        fn_opt['requested_nodes'] = max(min(50, len(chains_to_run) // int(fn_opt['cpus_per_node'])), 1)
        fn_opt['walltime'] = '0-03:55:00'
        fn_opt['waiting_time'] = '300'
        fn_opt['hpc_singularity'] = ""
        fn_opt['singularity'] = ""
        fn_opt['extra_outer_statements'] = "export OMP_NUM_THREADS=1; ulimit -u 4096"
        command = "source ~/encompass/bin/activate; python pairwise_neighbors.py <arg0> <arg1> <arg2> <arg3>"
        run_locusts(fn_opt, locations_path, locations, runs_label, 'neighbors', locusts_tmp, list(chains_to_run), command, opt_arg = [], gather = False)

        result, negresult = check_completion(locations, ch_list, 'neighbors')
        if len(set(ch_list).difference(completed).difference(result)) > 0:
            fail_flag = 1
        print("The pairwise neighbors procedure failed for the following structures (" + str(len(negresult)) + "):")
        for pdb in negresult:
            print(pdb + ' ', end = '')

    # Step 15: Finding the sequence and structural neighbors of each structure
    if step_name == 'sort_neighbors' or step_name == 'all':
        # Rank MSSD symmetries
        from .sort_symmetries import rank_by_symmetry
        ranked_symm_chains, ranked_symm_complex = rank_by_symmetry(locations)
        print("Finished sorting symmetry.")

        # Combine results from pairwise neighbors
        from .neighbor_lists import (
            combine_pairwise_neighbor_lists,
            merge,
            write_chain_neighbors,
            choose_chain_representative_based_on_symm,
            write_transfer_list,
        )
        pairwise_data_dir = locations['FSYSPATH']['neighbors']
        mode = "pairwise_neighbors"
        generate = True
        if generate == True:
            chain_pairs = combine_pairwise_neighbor_lists(pairwise_data_dir, suffix = mode + ".pkl")
            chain_clusters = merge(chain_pairs)
            pkl.dump(chain_clusters, open(pairwise_data_dir + "all_chain_clusters_" + mode + ".pkl", "wb"))
        else:
            chain_clusters = pkl.load(open(pairwise_data_dir + "all_chain_clusters_" + mode +".pkl", "rb"))
        print("Clustered neighbors.")
        ch_list = [pdb + "_" + ch for pdb in tm_archive.keys() for ch in tm_archive[pdb]['tmchains']]
        #check_unity(ch_list, chain_clusters)
        write_chain_neighbors(ch_list, chain_clusters, locations['FSYSPATH']['neighbors'] + mode + "_4_each_chain_seqid0.85_tmscore_0.6.txt")
        # Choose representative chain based on symmetry
        unique_chains, transfer_dic = choose_chain_representative_based_on_symm(ch_list, ranked_symm_chains, chain_clusters, tm_archive)
        #unique_chains contains structures with no neighbors and no symmetry
        pkl.dump(transfer_dic, open(pairwise_data_dir + "parent_children_symmetry-ranked_neighbors.pkl", 'wb'))
        print("Selected representative chains.")
        write_transfer_list(transfer_dic, pairwise_data_dir + "parent_children_symmetry-ranked_neighbors.txt")

    # Step 16: Inferred symmetry
    if step_name == 'transfer' or step_name == 'all':
        parent_children_dic = pkl.load(open(locations['FSYSPATH']['neighbors'] + "parent_children_symmetry-ranked_neighbors.pkl", 'rb'))
        ch_list = [pdb + "_" + ch for pdb in tm_archive.keys() for ch in tm_archive[pdb]['tmchains']]
        completed, failed = check_completion(locations, ch_list, 'transfer')
        print("Completed (", len(completed), ") :", completed, "\n\nFailed (", len(failed), "): ", failed)
        frtm_path = "/EncoMPASS/frtmalign/frtmalign"
        template_target_list = []
        for template in parent_children_dic:
            template_target_list = template_target_list + [[template, target] for target in parent_children_dic[template] if not target in completed]
        print(f"{len(template_target_list)} template-target pairs will be submitted. The total number of templates is {len(parent_children_dic)}.")

        fn_opt['partition'] = 'quick,norm'
        fn_opt['memory'] = '40g'
        fn_opt['constraint'] = ''
        fn_opt['cpus_per_node'] = 1
        fn_opt['requested_nodes'] = '60'
        fn_opt['walltime'] = '0-03:45:00'
        fn_opt['waiting_time'] = '600'

        command = "source ~/encompass/bin/activate; python inferred_symmetry.py <arg0> <arg1> <arg2> <arg3> <arg4> <arg5>" # double check this!!!!
        run_locusts(fn_opt, locations_path, locations, runs_label, 'transfer', locusts_tmp, template_target_list, command, opt_arg = [frtm_path], gather = True)
        result, negresult = check_completion(locations, ch_list, 'transfer')
        print(f"The transfer procedure completed {len(result)} structures and failed for the following structures {len(negresult)}:")
        # for pdb in negresult:
        #     print(pdb + ' ', end = '')
        # print('\nInferred symmetries completed.')

    # Step 17: Inferred symmetry images
    if step_name == "images_transfer" or step_name == "all":
        if os.path.isdir(locations['FSYSPATH']['symmtemp'] + "transfer/") == False:
            os.mkdir(locations['FSYSPATH']['symmtemp'] + "transfer/")
        parent_children_dic = pkl.load(open(locations['FSYSPATH']['neighbors'] + "parent_children_symmetry-ranked_neighbors.pkl", 'rb'))
        print(len(parent_children_dic))
        ch_list = [pdb + "_" + ch for pdb in tm_archive.keys() for ch in tm_archive[pdb]['tmchains']]
        completed, failed = check_completion(locations, ch_list, 'images_transfer')
        template_target_list = []
        for template in parent_children_dic:
            template_target_list = template_target_list + [[template, target] for target in parent_children_dic[template] if not target in completed]
        print(f"{len(template_target_list)} template-target pairs will be submitted.")
        pymol_path = '"pymol -c"'
        fn_opt['partition'] = 'quick,norm'
        fn_opt['memory'] = '58g'
        fn_opt['constraint'] = ''
        fn_opt['cpus_per_node'] = 1
        fn_opt['requested_nodes'] = '40'
        fn_opt['walltime'] = '0-03:30:00'
        fn_opt['waiting_time'] = '600'
        fn_opt['hpc_singularity'] = ""
        fn_opt['singularity'] = ""
        fn_opt['extra_outer_statements'] = "module load pymol/2.5.0"
        command = "source ~/encompass/bin/activate; python inferred_symmetry_images.py <arg0> <arg1> <arg2> <arg3> <arg4> <arg5>"
        run_locusts(fn_opt, locations_path, locations, runs_label, 'images_transfer', locusts_tmp, template_target_list, command, opt_arg = [pymol_path], gather = False)
        result, negresult = check_completion(locations, ch_list, 'images_transfer')

        print(f"The transfer images generation procedure completed {len(result)} structures and failed for the following structures {len(negresult)}:")
        for pdb in negresult:
            print(pdb + ' ', end = '')
        print('\nTransfer images step completed.')

    # Step 13: Combine transfer results
    if step_name == "transfer_combine_dic" or step_name == "all":
        ch_list = [pdb + "_" + ch for pdb in tm_archive.keys() for ch in tm_archive[pdb]['tmchains']]
        completed, failed = check_completion(locations, ch_list, 'transfer')
        print("Completed (", len(completed), ") :", completed, "\n\nFailed (", len(failed), "): ", failed)

        results = {}
        for inputf in completed:
            print(inputf)
            dic_path = locations['FSYSPATH']['transfer']  + inputf + "/" + inputf + "_transfer.pkl"
            if os.path.isfile(dic_path):
                d = pkl.load(open(dic_path,"rb"))
                results.update(d)

        pkl.dump(results, open(locations['FSYSPATH']['transfer'] + "transfer_results.pkl","wb"))

    return


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description = 'Execute a single step or all of the symmetry analysis. '
                                                   'Possible steps are: symd, cesymm, cesymm_low_thr, cesymm_from_symd, '
                                                   'cesymm_order, cesymm_minlen, cesymm_quat_minlen, cesymm_rmsd_2_5, '
                                                   'cesymm_rmsd_3, ananas, quatsymm, mssd, mssd_combine_dic, neighbors, '
                                                   'sort_neighbors, transfer, images_transfer, transfer_combine_dic, '
                                                   'all, test, clean (creates tm_dictionary only).'
                                                   'Runs can be further broken down based on types of proteins: '
                                                   'small (<=20 000 atoms), big (>20 000 atoms and <=8000 residues or giant (>8000 residues)')
    parser.add_argument('-step', type=str, required = True, choices=['symd','cesymm','cesymm_low_thr','cesymm_from_symd', 'cesymm_order', 'cesymm_minlen', 'cesymm_quat_minlen','cesymm_rmsd_2_5', 'cesymm_rmsd_3', 'ananas', 'quatsymm', 'images_cesymm_symd', 'mssd', 'mssd_combine_dic', 'neighbors', 'sort_neighbors', 'transfer', 'images_transfer', 'transfer_combine_dic', 'all', 'test', 'clean'], help = 'which step of the symmetry analysis should be performed')
    parser.add_argument('-label', nargs='?', required = True, help = 'name of the run and, if relevant, of the HPC folder')
    parser.add_argument('-locusts_dir', nargs='?', required = True, help = 'the directory that locusts should use for temporary storage')


    # Optional
    parser.add_argument('-db', nargs='?', default='/data/local/encompass/dev/toni_full_dec2020_beforeFrTM', help = 'path to the database version to use')
    parser.add_argument('-instr_file', nargs='?', help = 'instruction file for the already built database, e.g. xxx_relative_to_xxx.txt')
    parser.add_argument('-pdb_suff', type=str, nargs='?', default = "_enc", help = 'suffix in the name of the pdb files, eg. "_enc"')
    parser.add_argument('-str_type', type=str, default='none', choices=['small', 'big', 'giant', 'all', 'none'], help = 'which structures should be analyzed (small, big, giant, all, none)')
    parser.add_argument('-wkdir', nargs='?', help = 'where output should be saved if not in the database', default = "")
    parser.add_argument('-in_file', nargs='?', help = 'a file with a list of pdbs to be analyzed (one per line)', default = False)
    parser.add_argument('-clean_list', nargs='?', help = 'a filename with a list of pdbs (one on each line), which should be purged from the symmetry results')
    parser.add_argument('--generate', help = 'if present, all selection algorithms will be run to generate list of pdbs to run for the steps', action = "store_true", default = False)
    parser.add_argument('--locusts-debug', help = 'if present, will pass the debug signal to locusts',action = "store_true")
    parser.add_argument('--local', help = 'whether we want to run locally or on HPC', action = "store_true", default = False)

    return parser

def main(argv: list[str] | None = None) -> None:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    args.locusts_dir = os.path.join(args.locusts_dir, '') # add trailing slash
    args.db = args.db[:-1] if args.db[-1] == "/" else args.db # remove trailing slash
    print("instr file", args.instr_file)
    print("Initiating symmetry analysis with the following options:\n")
    print("Database path: {}\nResults working directory: {}\nInput file with pdbs: {}\nPdb suffix: {}\nStep name: {}\nRun label: {}\nLocusts temporary directory: {}\nLocal run: {}\n".format(args.db, args.wkdir, args.in_file, args.pdb_suff, args.step, args.label, args.locusts_dir, args.local))
    symm_update(locations_path=args.db,
                wkdir=args.wkdir,
                input_file=args.in_file,
                pdb_suff=args.pdb_suff,
                step_name=args.step,
                str_type=args.str_type,
                runs_label=args.label,
                locusts_tmp=args.locusts_dir,
                instr_path=args.instr_file,
                generate = args.generate,
                clean = args.clean_list,
                run_local=args.local)


if __name__ == "__main__":
    main()