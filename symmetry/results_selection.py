# Name: results_selection.py
# Author: Antoniya A. Aleksandrova
# Language: Python 3.5+
# Description: Select the symmetry option that best fits MSSD criteria


import sys
import os
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from symmetry_exec_functions import *
from initialize_repository import initialize_repository
import glob

def all_options_write(file, inputf, dic):
    file.write("%s\t%-23s\t\t\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
        inputf, dic['source'], dic['chains'],
        str(dic['repeats_number']), str(dic['symmetry_levels']),
        dic['symmetry_order'], dic['topology'],
        dic['symmetry_type'], str(dic['unit_angle']),
        str(dic['unit_translation']), str(dic['refined_tmscore']),
        str(dic['refined_rmsd']), str(dic['repeat_length']),
        str(dic['aligned_length']), str(dic['size']),
        str(dic['coverage'])))
    return

def selected_write(file, inputf, dic, status):
    file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
        inputf, dic['chains'], str(dic['repeats_number']),
        str(dic['symmetry_levels']), dic['symmetry_order'],
        dic['topology'], dic['symmetry_type'], str(dic['unit_angle']),
        str(dic['unit_translation']), str(dic['refined_tmscore']),
        str(dic['refined_rmsd']), str(dic['repeat_length']),
        str(dic['aligned_length']), str(dic['size']),
        str(dic['coverage']), dic['source'], status))
    return


def fragments_check(struct, ch, ch_pdb, sel, tm_arch, cs_data, locations, symmdic, oriented, doubletmhelices, ray, pymol, want_image, pic_dir, jmol_dir, super_dir, out2, image_files, pml_files, jmol_files, super_pml_files, tm_score_super):
    ### Is there information in the fragments?
    print('frag inputs:', struct, ch, ch_pdb)
    structid = struct + "_" + ch_pdb if ch_pdb != '' else struct
    fragments = {}
    cs_lt_out_dir = locations['FSYSPATH']['cesymm_low_thr'] + struct + "/"
    cs_from_symd_dir = locations['FSYSPATH']['cesymm_from_symd']
    # We start with CE-Symm Lower Threshold fragments because:
    #    a) fragments are run only if CE-Symm default has produced no result
    #    b) fragments are run with the lower-threshold only if the run with
    #       lower threshold on the whole protein produced symmetry that left
    #       unaligned more than 80 successive amino acids which could be
    #       considered a fragment
    # Therefore, it is not possible for these fragments to have any overlap with the repeats in the selected symmetry
    for f in glob.glob(cs_lt_out_dir + struct + "." + ch + "*.axes"):
        f_name = f[len(cs_lt_out_dir):-5]
        fragment_dic = cs_data(cs_lt_out_dir, f_name)
        if fragment_dic['repeats_number'] != 'na' and int(fragment_dic['repeats_number']) > 1 and int(
                fragment_dic['repeat_length']) >= doubletmhelices and \
                int(fragment_dic['aligned_length']) >= 2 * doubletmhelices:
            reference = parse_structure(oriented)[0]
            repeat_selection = get_repeat_resid(cs_lt_out_dir, f_name, reference)
            fragment_dic['repeat_selection'] = repeat_selection
            fragment_dic['files_key'] = f_name
            fragment_dic['files_dir'] = locations['FSYS']['cesymm_low_thr'] + struct + "/"
            fragment_dic['files_dir_path'] = locations['FSYSPATH']['cesymm_low_thr'] + struct + "/"
            fragment_dic['oriented'] = oriented
            fragment_dic['chains'] = ch + ";"
            fragment_dic['source'] = 'fragment_cesymm_lt'
            tmscore, rmsd, angle, prot_len = partial_repeats_scores(fragment_dic['files_dir_path'], oriented, f_name)
            fragment_dic['refined_tmscore'] = tmscore
            fragment_dic['coverage'] = float('%.2f' % (fragment_dic['aligned_length'] / prot_len))
            fragment_dic['size'] = prot_len
            fragments['clt'] = fragment_dic
            fragments['source'] = 'clt;'
            all_options_write(out2, struct + "_" + ch, fragment_dic)

        # Next we investigate the SymD-suggested fragments and check if they overlap with the selected symmetry or
    for i, f in enumerate(glob.glob(cs_from_symd_dir + structid + "*")):
        fragments_out_dir = os.path.join(f, '') # add a /
        f_name = f[len(cs_from_symd_dir):]
        fragments_out_dir_rel = locations['FSYS']['cesymm_from_symd'] + f_name + "/"
        f_name = f_name[0:4] + "." + f_name[5:]
        fragment_dic = cs_data(fragments_out_dir, f_name)
        if fragment_dic['repeats_number'] != 'na' and int(fragment_dic['repeats_number']) > 1 and int(
                fragment_dic['repeat_length']) >= doubletmhelices and \
                int(fragment_dic['aligned_length']) >= 2 * doubletmhelices:
            reference = parse_structure(oriented)[0]
            repeat_selection = get_repeat_resid(fragments_out_dir, f_name, reference)
            overlap_sel = repeats_overlap(repeat_selection, symmdic[structid][sel][
                'repeat_selection'])  # overlap with the selected symmetry
            if 'clt' in fragments:
                overlap = repeats_overlap(repeat_selection, fragments['clt']['repeat_selection'])  # overlap
                if overlap == 'independent' and overlap_sel == 'independent':
                    fragment_dic['repeat_selection'] = repeat_selection
                    fragment_dic['files_key'] = f_name
                    fragment_dic['files_dir'] = fragments_out_dir_rel
                    fragment_dic['files_dir_path'] = fragments_out_dir
                    fragment_dic['oriented'] = oriented
                    fragment_dic['chains'] = ch + ";"
                    fragment_dic['source'] = 'fragment_symd'
                    tmscore, rmsd, angle, prot_len = partial_repeats_scores(fragment_dic['files_dir_path'], oriented,
                                                                            f_name)
                    fragment_dic['refined_tmscore'] = tmscore
                    fragment_dic['coverage'] = float('%.2f' % (fragment_dic['aligned_length'] / prot_len))
                    fragment_dic['size'] = prot_len
                    fragments['symd_' + str(i)] = fragment_dic
                    fragments['source'] = fragments['source'] + 'symd_' + str(i) + ';'
                    all_options_write(out2, structid, fragment_dic)


            elif overlap_sel == 'independent':
                fragment_dic['repeat_selection'] = repeat_selection
                fragment_dic['files_key'] = f_name
                fragment_dic['files_dir'] = fragments_out_dir_rel
                fragment_dic['files_dir_path'] = fragments_out_dir
                fragment_dic['oriented'] = oriented
                fragment_dic['chains'] = ch + ";"
                fragment_dic['source'] = 'fragment_symd'
                tmscore, rmsd, angle, prot_len = partial_repeats_scores(fragment_dic['files_dir_path'], oriented,
                                                                        f_name)
                fragment_dic['refined_tmscore'] = tmscore
                fragment_dic['coverage'] = float('%.2f' % (fragment_dic['aligned_length'] / prot_len))
                fragment_dic['size'] = prot_len
                fragments['symd_' + str(i)] = fragment_dic
                fragments['source'] = fragments.get('source', '') + 'symd_' + str(i) + ';'
                all_options_write(out2, structid, fragment_dic)



    if 'source' in fragments:

        print('sel', sel, symmdic[structid][sel])
        doms, crossings = symmetry_tm_domains(tm_arch, symmdic[structid][sel]['repeat_selection'], struct, ch_pdb)
        source = fragments['source'][:-1].split(';')
        symm_location_frag = 'out'
        for i in source:
            doms_tmp, crossings_tmp = symmetry_tm_domains(tm_arch, fragments[i]['repeat_selection'], struct, ch_pdb)
            if doms_tmp == 'out':
                re.sub(i + ';', '', fragments['source'])  # removes fragments that have their repeats outside the membrane
            else:
                symm_location_frag = 'in'
                fragments[i]['rep_has_tm_domains'] = doms_tmp
                fragments[i]['rep_tm_domain_num'] = crossings_tmp
        if (symmdic[structid][sel]['repeats_number'] == 'na' or int(symmdic[structid][sel]['repeats_number']) == 1
            or int(symmdic[structid][sel]['aligned_length']) < 2 * doubletmhelices or doms == 'out') \
                and symm_location_frag == 'in':
            source = fragments['source'][:-1].split(';')
            if len(source) == 1:
                s = source[0]
                symmdic[structid]['selected']['source'] = s
            else:
                symmdic[structid]['selected']['source'] = fragments['source']
            for s in source:
                symmdic[structid][s] = fragments[s]
                reference = parse_structure(oriented)[0]
                image_key = structid + "_" + s + "_analysis"
                if want_image == 1:
                    symmdic[structid][s]['image_files'], \
                    symmdic[structid][s]['pml_files'], \
                    symmdic[structid][s]['jmol_files'] = cesymm_images(wkdir, fragments[s]['files_dir_path'], fragments[s]['files_key'],
                                                                     fragments[s]['repeat_selection'], oriented,
                                                                     reference, ray, pymol, pic_dir, jmol_dir,
                                                                     image_key)
                    symmdic[structid][s]['super_pml_files'], tm_score_super = repeats_superposition(fragments[s]['files_dir_path'], oriented,
                                                                                                  fragments[s]['files_key'], super_dir,
                                                                                                  image_key, fragments[s]['repeat_selection'])
                    symmdic[structid][s]['image_key'] = image_key


        elif symm_location_frag == 'in':
            symmdic[structid][sel]['rep_has_tm_domains'] = doms
            symmdic[structid][sel]['rep_tm_domain_num'] = crossings
            symmdic[structid]['selected']['source'] = symmdic[structid]['selected']['source'] + ';' + fragments['source']
            s = fragments['source'].strip(';')
            symmdic[structid][s] = fragments[s]
            reference = parse_structure(oriented)[0]
            image_key = structid + "_" + s +  "_analysis"
            if want_image == 1:
                symmdic[structid][s]['image_files'], \
                symmdic[structid][s]['pml_files'], \
                symmdic[structid][s]['jmol_files'] = cesymm_images(wkdir, fragments[s]['files_dir_path'], fragments[s]['files_key'],
                                                                 fragments[s]['repeat_selection'], oriented,
                                                                 reference, ray, pymol, pic_dir, jmol_dir,
                                                                 image_key)
                symmdic[structid][s]['super_pml_files'], tm_score_super = repeats_superposition(fragments[s]['files_dir_path'], oriented,
                                                                                              fragments[s]['files_key'], super_dir,
                                                                                              image_key, fragments[s]['repeat_selection'])
                symmdic[structid][s]['image_key'] = image_key

    return symmdic, fragments, tm_score_super

def results_selection(locations, options, pdb_suff, inputf, run_id="0000", want_image = 1, ray = 0):
    """
    This function executes the multi-step symmetry selection (MSSD) analysis on a given protein structure.

    :param locations: dictionary with all the relevant paths to folders, executables and lists
    :param options:
    :param pdb_suff: the suffix attached to the pdb code in the name of the structure (e.g. _enc)
    :param inputf: the pdb code of the structure (e.g. 1okc)
    :param run_id: an identifier for this run
    :param want_image: should images be produced? (default: 1, i.e. yes)
    :param ray: should PyMOL images be ray-traced? (default: 0, i.e. no)
    :return: a dictionary with all the symmetry information for the structure and its TM chains
    """
    if not os.path.isdir(locations['FSYSPATH']['mssd_log'] + inputf):
        os.mkdir(locations['FSYSPATH']['mssd_log'] + inputf)


    log_file = open(locations['FSYSPATH']['mssd_log'] + inputf + "/" + inputf + '_run.log', 'w')
    save_out = sys.stdout
    save_err = sys.stderr
    sys.stdout = sys.stderr = log_file

    print("Run id: %s" % run_id)

    if os.path.isfile(locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb') == False:
        raise SystemExit("{0} does not exist.".format(locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb'))
    tm_archive = locations['FSYSPATH']['symmetry'] + '.tm_archive.pkl'
    oriented_opm = locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb'
    symd_out_dir = locations['FSYSPATH']['symd'] + inputf + "/"
    cesymm_out_dir = locations['FSYSPATH']['cesymm'] + inputf + "/"
    cesymm_lt_out_dir = locations['FSYSPATH']['cesymm_low_thr'] + inputf + "/"
    cesymm_minlen_out_dir = locations['FSYSPATH']['cesymm_minlen'] # followed by eg 4f35_df
    cesymm_quat_minlen_out_dir = locations['FSYSPATH']['cesymm_quat_minlen'] # followed by eg 4f35_df
    quatsymm = locations['FSYSPATH']['quatsymm'] + inputf + "/"
    cesymm_from_symd_dir = locations['FSYSPATH']['cesymm_from_symd'] # followed by eg 4f35_df

    chain_ordered_dir = locations['FSYSPATH']['symmtemp'] + "chain_ordered_str/"
    pymol = "pymol -c"
    muscle_exe = "muscle"
    selected_out = locations['FSYSPATH']['mssd']  + inputf + "/"
    if not os.path.isdir(selected_out):
        os.mkdir(selected_out)
    wkdir = locations['FSYSPATH']['symmtemp']  + "mssd/"
    pic_dir_whole = locations['FSYSPATH']['analysis_wholepngs']
    pic_dir_chain = locations['FSYSPATH']['analysis_chainspngs']
    jmol_dir_whole = locations['FSYSPATH']['analysis_wholejsons']
    jmol_dir_chain = locations['FSYSPATH']['analysis_chainsjsons']
    super_dir_whole = locations['FSYSPATH']['analysis_wholesuper']
    super_dir_chain = locations['FSYSPATH']['analysis_chainssuper']


    doubletmhelices = 40
    pdb = inputf[0:4] # inputf should already be a 4-letter code
    print(inputf)
    # Find TM chains
    a = pkl.load(open(tm_archive, 'rb'))
    if len(a[inputf]['tmchains']) > 0:
        chains_str = ';'.join(sorted(a[inputf]['tmchains'])) + ';'  # tm chains as a string A;B;C;...
        chains = sorted(a[inputf]['tmchains'])
        print(chains_str)
        type = a[inputf]['class']
    else:
        raise SystemExit("{0} has no TM chains.".format(locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb'))

    if type == 'beta':
        print("{0} is a beta-barrel so MSSD will not analyze it.".format(inputf))
        return {}
    if a[inputf]['unk'] == True:
        print("{0} contains UNK residues, which may cause errors.".format(inputf))
    #    return {}


    #######
    out1_name = locations['FSYSPATH']['mssd_log'] + inputf + "/" + inputf + "_selected.txt"
    out2_name = locations['FSYSPATH']['mssd_log'] + inputf + "/" + inputf + "_all_options_list.txt"
    out3_name = locations['FSYSPATH']['mssd_log'] + inputf + "/" + inputf + "_changes_from_default.txt"

    out1 = open(out1_name, "w")
    out2 = open(out2_name, "w")
    out3 = open(out3_name, "w")  # Record statistics on how much the analysis changes CE-Symm results
    out3.write(inputf + "\t")

    cesymm_quat_symm = 0
    discard_small = 0
    discard_nonmemb = 0

    # Prepare structure
    num = strip_tm_chains(wkdir, inputf, oriented_opm, chains)  # should not strip them in order because they were not run that way
    oriented = wkdir + inputf + "_tmp.pdb"
    pdb_chain_order = chain_id(oriented)

    ########## Quaternary Symmetry
    # Create a dictionary of all the CE-Symm and SymD symmetry features (does not require a loaded structure)
    if pdb_chain_order != chains_str:
        chains_str = pdb_chain_order
        chains = pdb_chain_order[:-1].split(';')
    symmlist = {}
    if os.path.isfile(locations['FSYSPATH']['images_cesymm_symd_log'] + inputf + "/" + inputf + "_dic.pkl"):
        symmlist = pkl.load(open(locations['FSYSPATH']['images_cesymm_symd_log'] + inputf + "/" + inputf + "_dic.pkl", "rb"))
        symmlist[inputf]['type'] = type
    else:
        print("No dictionary from the cesymm/symd data and images step. Will attempt to read in the raw data now.\n")
        symmlist[inputf] = {'type': type}
        cesymm_dic = cesymm_data(cesymm_out_dir, inputf)
        symmlist[inputf]['cesymm'] = {}
        symmlist[inputf]['cesymm'] = cesymm_dic
        symd_dic = symd_data(symd_out_dir, inputf)
        symmlist[inputf]['symd'] = {}
        symmlist[inputf]['symd'] = symd_dic

    symmlist[inputf]['cesymm']['chains'] = chains_str
    symmlist[inputf]['cesymm']['source'] = 'cesymm'
    symmlist[inputf]['cesymm']['oriented'] = oriented
    all_options_write(out2, inputf, symmlist[inputf]['cesymm'])

    cesymm_dic = {}
    cesymm_dic = cesymm_data(cesymm_lt_out_dir, inputf)
    symmlist[inputf]['cesymm_lt'] = {}
    symmlist[inputf]['cesymm_lt'] = cesymm_dic
    symmlist[inputf]['cesymm_lt']['chains'] = chains_str
    symmlist[inputf]['cesymm_lt']['source'] = 'cesymm_lt'
    symmlist[inputf]['cesymm_lt']['oriented'] = oriented

    symmlist[inputf]['selected'] = {}
    symmlist[inputf]['selected']['source'] = 'cesymm'

    limit_bottom, limit_top = membrane_limits(oriented_opm)

    ### Is there information in the default run?
    if symmlist[inputf]['cesymm']['repeats_number'] == 'na':
        images = ""
        pml = ""
        pass
    elif int(symmlist[inputf]['cesymm']['repeats_number']) > 1:
        symmlist[inputf]['selected']['source'] = 'cesymm'
        reference = parse_structure(oriented)[0]
        repeat_selection = get_repeat_resid(cesymm_out_dir, inputf, reference)
        symm_location = is_symmetry_in_membrane(limit_bottom, limit_top, repeat_selection, oriented)
        doms, crossings = symmetry_tm_domains(tm_archive, repeat_selection, inputf, '')
        symmlist[inputf]['cesymm']['rep_has_tm_domains'] = doms
        symmlist[inputf]['cesymm']['rep_tm_domain_num'] = crossings
        symmlist[inputf]['cesymm']['symm_within_membrane'] = symm_location
        if doms == 'in':
            cesymm_quat_symm = 1
        elif symm_location == 'out':
            discard_nonmemb = 1
        else:
            discard_small = 1

    # Account for multi-chain cases where CE-Symm declares no symmetry even though it's unrefined TM-score is 1
    elif len(chains) > 1 \
            and (int(symmlist[inputf]['cesymm']['repeats_number']) == 1
                 or int(symmlist[inputf]['cesymm']['repeat_length']) < doubletmhelices) \
            and (float(symmlist[inputf]['cesymm']['unrefined_tmscore']) > 0.90
                 or (float(symmlist[inputf]['symd']['tmscore']) > 0.90
                     and float(symmlist[inputf]['symd']['coverage']) > 0.90)):
        # exmaple: 1yew: A;B;C;E;F;G;I;J;K; angle=120
        # find the angle of symmetry
        k = 0
        if symmlist[inputf]['cesymm']['unit_angle'] != 'na' \
                and len(symmlist[inputf]['cesymm']['unit_angle'][:-1].split(';')) == 1 \
                and float(symmlist[inputf]['cesymm']['unit_angle'].strip(';')) != 0 \
                and float(symmlist[inputf]['cesymm']['unrefined_tmscore']) > 0.90:
            k = round(360 / float(symmlist[inputf]['cesymm']['unit_angle'].strip(';')))
        if (k != 0 and len(chains) % k != 0) or k == 0 or round(k) == 1:
            k = round(360 / float(symmlist[inputf]['symd']['unit_angle']))  # k=3
            if len(chains) % k != 0 or round(k) == 1:
                k = round(360 / float(symmlist[inputf]['symd']['is_angle']))
        # divide chains in repeats
        if len(chains) % k == 0 and round(k) != 1:
            n = int(len(chains) / k)  # n=3
            repeats_list = [chains[i:i + n] for i in range(0, len(chains), n)]  # repeats_list=[[A,B,C],[E,F,G],[I,J,K]]
            cs_bug_fix_dic = bug_fixer(muscle_exe, wkdir, inputf, oriented, repeats_list, "_analysis")
            symmlist[inputf]['cesymm_bug_fix'] = cs_bug_fix_dic
            if len(cs_bug_fix_dic) > 0:
                cs_bug_fix_dic['chains'] = chains_str
                cs_bug_fix_dic['source'] = 'cesymm_bug_fix'
                cs_bug_fix_dic['oriented'] = oriented
                symmlist[inputf]['cesymm_bug_fix'] = cs_bug_fix_dic
                symmlist[inputf]['selected']['source'] = 'cesymm_bug_fix'
                all_options_write(out2, inputf, symmlist[inputf]['cesymm_bug_fix'])


    all_names = ['_ordered', '_com_sd', '_sd_order']
    for name in all_names:
        # Handle chain re-ordering
        if os.path.isfile(cesymm_out_dir + inputf + name + "_stdout.out"):
            name_dic = 'cesymm' + name
            symmlist[inputf][name_dic] = {}
            symmlist[inputf][name_dic] = cesymm_data(cesymm_out_dir, inputf + name)
            symmlist[inputf][name_dic]['source'] = name_dic
            oriented_tmp = chain_ordered_dir + inputf + name + ".pdb"
            symmlist[inputf][name_dic]['oriented'] = oriented_tmp
            chains_str_tmp = find_chains(oriented_tmp)
            chains_tmp = chains_str_tmp[:-1].split(';')
            symmlist[inputf][name_dic]['chains'] = chains_str_tmp
            print(chains_tmp)

            # Fix the CE-Symm bug with refinement
            symd_dic = 'symd' + name
            symmlist[inputf][symd_dic] = {}
            symmlist[inputf][symd_dic] = symd_data(symd_out_dir, inputf + name)
            if symmlist[inputf][name_dic]['repeats_number'] != 'na' \
                    and len(chains_tmp) > 1 \
                    and (int(symmlist[inputf][name_dic]['repeats_number']) == 1
                         or int(symmlist[inputf][name_dic]['repeat_length']) < doubletmhelices) \
                    and (float(symmlist[inputf][name_dic]['unrefined_tmscore']) > 0.90
                         or (float(symmlist[inputf][symd_dic]['tmscore']) > 0.90
                             and float(symmlist[inputf][symd_dic]['coverage']) > 0.90)):
                # exmaple: 1yew: A;B;C;E;F;G;I;J;K; angle=120
                # find the angle of symmetry
                k = 0
                if symmlist[inputf][name_dic]['unit_angle'] != 'na' \
                        and len(symmlist[inputf][name_dic]['unit_angle'][:-1].split(';')) == 1 \
                        and float(symmlist[inputf][name_dic]['unit_angle'].strip(';')) != 0 \
                        and float(symmlist[inputf][name_dic]['unrefined_tmscore']) > 0.90:
                    k = round(360 / float(symmlist[inputf][name_dic]['unit_angle'].strip(';')))
                if (k != 0 and len(chains_tmp) % k != 0) or k == 0 or round(k) == 1:
                    k = round(360 / float(symmlist[inputf][symd_dic]['unit_angle']))  # k=3
                    if len(chains_tmp) % k != 0 or round(k) == 1:
                        k = round(360 / float(symmlist[inputf][symd_dic]['is_angle']))
                # divide chains_tmp in repeats
                if len(chains_tmp) % k == 0 and round(k) != 1:
                    n = int(len(chains_tmp) / k)  # n=3
                    repeats_list = [chains_tmp[i:i + n] for i in
                                    range(0, len(chains_tmp), n)]  # repeats_list=[[A,B,C],[E,F,G],[I,J,K]]
                    cs_bug_fix_dic = bug_fixer(muscle_exe, wkdir, inputf + name, oriented_tmp, repeats_list, "_analysis")
                    if len(cs_bug_fix_dic) > 0:
                        cs_bug_fix_dic['chains'] = chains_str_tmp
                        symmlist[inputf][name_dic]['chains'] = chains_str_tmp
                        all_options_write(out2, inputf, symmlist[inputf][name_dic])
                        symmlist[inputf][name_dic] = cs_bug_fix_dic
                        symmlist[inputf][name_dic]['oriented'] = oriented_tmp
                        symmlist[inputf][name_dic]['source'] = name_dic + '_bug_fix'

            # Were the analyzed chains just a subset of the TM chains?
            additions = ""
            for ch in chains:  # K;3;2;J;F;4;1;G;B;R;A;L;H;I; but only B;A; appear in the symmetry order (2wsc)
                if ch not in chains_str_tmp:
                    additions = additions + ch + ";"

            # Re-calculate scores if only a subset of chains was analyzed by CE-Symm
            if len(additions) != 0 and symmlist[inputf][name_dic]['repeats_number'] != 'na' and int(
                    symmlist[inputf][name_dic]['repeats_number']) > 1:
                chains_str_tmp = chains_str_tmp + additions
                chains_tmp = chains_str_tmp[:-1].split(';')
                symmlist[inputf][name_dic]['chains'] = chains_str_tmp
                num = strip_tm_chains_in_order(wkdir, inputf + name, oriented_opm, chains_tmp)
                oriented_tmp = wkdir + inputf + name + "_tmp.pdb"
                if symmlist[inputf][name_dic]['source'] == name_dic + '_bug_fix':
                    tm_score, rmsd, angle, prot_len = partial_repeats_scores(wkdir, oriented_tmp,
                                                                             symmlist[inputf][name_dic]['files_key'])
                else:
                    tm_score, rmsd, angle, prot_len = partial_repeats_scores(cesymm_out_dir, oriented_tmp, inputf + name)
                if abs(rmsd - float(symmlist[inputf][name_dic]['refined_rmsd'])) < 0.2:
                    symmlist[inputf][name_dic]['refined_tmscore'] = tm_score
                    symmlist[inputf][name_dic]['coverage'] = float(
                        '%.2f' % (symmlist[inputf][name_dic]['aligned_length'] / prot_len))
                    symmlist[inputf][name_dic]['size'] = prot_len
                    symmlist[inputf][name_dic]['oriented'] = oriented_tmp

            all_options_write(out2, inputf, symmlist[inputf][name_dic])

            # Compare with the current best
            sel = symmlist[inputf]['selected']['source']
            #   if no symmetry was found before but now we have one (with only quaternary symmetries), select it
            if (symmlist[inputf][sel]['repeats_number'] == 'na'
                or int(symmlist[inputf][sel]['repeats_number']) == 1
                or (int(symmlist[inputf][sel]['repeat_length']) < doubletmhelices and int(
                        symmlist[inputf][name_dic]['repeat_length']) > doubletmhelices)
                or (symmlist[inputf][sel]['internal_symmetry'] == 'Yes'
                    and symmlist[inputf][name_dic]['internal_symmetry'] == 'No')) \
                    and symmlist[inputf][name_dic]['repeats_number'] != 'na' \
                    and int(symmlist[inputf][name_dic]['repeats_number']) > 1 \
                    and int(symmlist[inputf][name_dic]['repeat_length']) > doubletmhelices:
                symmlist[inputf]['selected']['source'] = name_dic

            #  otherwise, select the current symmetry if either it both has more symmetry levels and a repeat length >doubletmhelices
            # or both has the same number of levels and greater coverage, provided they both have or do not have internal symmetry
            #  or has fewer levels, same symmetry order and greater or same coverage (eg. 3k0g)
            elif symmlist[inputf][sel]['repeats_number'] != 'na' \
                    and int(symmlist[inputf][sel]['repeats_number']) > 1 \
                    and symmlist[inputf][name_dic]['repeats_number'] != 'na' \
                    and int(symmlist[inputf][name_dic]['repeats_number']) > 1:
                if int(symmlist[inputf][sel]['symmetry_levels']) < int(symmlist[inputf][name_dic]['symmetry_levels']) \
                        and symmlist[inputf][sel]['internal_symmetry'] == symmlist[inputf][name_dic]['internal_symmetry'] \
                        and float(symmlist[inputf][sel]['coverage']) <= float(symmlist[inputf][name_dic]['coverage']) \
                        and int(symmlist[inputf][sel]['repeats_number']) == int(
                    symmlist[inputf][name_dic]['repeats_number']):
                    pass
                # potential error next (elif vs if)
                elif (int(symmlist[inputf][sel]['symmetry_levels']) < int(symmlist[inputf][name_dic]['symmetry_levels'])
                    and (int(symmlist[inputf][name_dic]['repeat_length']) >= doubletmhelices
                         and symmlist[inputf][name_dic]['internal_symmetry'] == symmlist[inputf][sel]['internal_symmetry']
                         and float(symmlist[inputf][sel]['coverage']) <= float(symmlist[inputf][name_dic]['coverage']))) \
                        or (
                        int(symmlist[inputf][sel]['symmetry_levels']) == int(symmlist[inputf][name_dic]['symmetry_levels'])
                        and float(symmlist[inputf][sel]['coverage']) < float(symmlist[inputf][name_dic]['coverage'])
                        and symmlist[inputf][name_dic]['internal_symmetry'] == symmlist[inputf][sel]['internal_symmetry']) \
                        or (
                        int(symmlist[inputf][sel]['symmetry_levels']) > int(symmlist[inputf][name_dic]['symmetry_levels'])
                        and symmlist[inputf][sel]['symmetry_order'] == symmlist[inputf][name_dic]['symmetry_order']
                        and float(symmlist[inputf][sel]['coverage']) <= float(symmlist[inputf][name_dic]['coverage']))\
                        or (int(symmlist[inputf][sel]['symmetry_levels']) == int(symmlist[inputf][name_dic]['symmetry_levels'])
                        and float(symmlist[inputf][sel]['coverage']) == float(symmlist[inputf][name_dic]['coverage'])
                        and int(symmlist[inputf][sel]['repeats_number']) < int(symmlist[inputf][name_dic]['repeats_number'])
                        and symmlist[inputf][name_dic]['internal_symmetry'] == "No"):
                    symmlist[inputf]['selected']['source'] = name_dic
                    # Clean up
        #    if len(additions)!=0:
        #        os.remove(wkdir+inputf+"_ordered_tmp.pdb")

    ### Is there information in the lower threshold run?
    if symmlist[inputf]['cesymm_lt']['repeats_number'] == 'na' or int(symmlist[inputf]['cesymm_lt']['repeats_number']) == 1:
        all_options_write(out2, inputf, symmlist[inputf]['cesymm_lt'])
        pass
    elif int(symmlist[inputf]['cesymm_lt']['repeats_number']) > 1:
        # Compare with the current best
        sel = symmlist[inputf]['selected']['source']
        #   if no symmetry was found before but now we have one with repeat length at least doubletmhelices, select it
        if (symmlist[inputf][sel]['repeats_number'] == 'na'
            or int(symmlist[inputf][sel]['repeats_number']) == 1
            or (int(symmlist[inputf][sel]['repeat_length']) < doubletmhelices and int(
                    symmlist[inputf]['cesymm_lt']['repeat_length']) > doubletmhelices)
            or (symmlist[inputf][sel]['internal_symmetry'] == 'Yes'
                and symmlist[inputf]['cesymm_lt']['internal_symmetry'] == 'No')) \
                and symmlist[inputf]['cesymm_lt']['repeats_number'] != 'na' \
                and int(symmlist[inputf]['cesymm_lt']['repeats_number']) > 1 \
                and int(symmlist[inputf]['cesymm_lt']['repeat_length']) > doubletmhelices:
            symmlist[inputf]['selected']['source'] = 'cesymm_lt'

        #  otherwise, select the current symmetry if either it both has more symmetry levels and a repeat length >doubletmhelices
        # or both has the same number of levels and greater coverage
        #  or has fewer levels, same symmetry order and greater or same coverage (eg. 3k0g)
        elif symmlist[inputf][sel]['repeats_number'] != 'na' \
                and int(symmlist[inputf][sel]['repeats_number']) > 1 \
                and symmlist[inputf]['cesymm_lt']['repeats_number'] != 'na' \
                and int(symmlist[inputf]['cesymm_lt']['repeats_number']) > 1:
            if int(symmlist[inputf][sel]['symmetry_levels']) < int(symmlist[inputf]['cesymm_lt']['symmetry_levels']) \
                    and symmlist[inputf][sel]['internal_symmetry'] == symmlist[inputf]['cesymm_lt']['internal_symmetry'] \
                    and float(symmlist[inputf][sel]['coverage']) <= float(symmlist[inputf]['cesymm_lt']['coverage']) \
                    and int(symmlist[inputf][sel]['repeats_number']) == int(
                symmlist[inputf]['cesymm_lt']['repeats_number']):
                pass
            elif (int(symmlist[inputf][sel]['symmetry_levels']) < int(symmlist[inputf]['cesymm_lt']['symmetry_levels'])
                  and (int(symmlist[inputf]['cesymm_lt']['repeat_length']) >= doubletmhelices
                       and symmlist[inputf]['cesymm_lt']['internal_symmetry'] == symmlist[inputf][sel]['internal_symmetry']
                       and float(symmlist[inputf][sel]['coverage']) <= float(symmlist[inputf]['cesymm_lt']['coverage']))) \
                    or (
                    int(symmlist[inputf][sel]['symmetry_levels']) == int(symmlist[inputf]['cesymm_lt']['symmetry_levels'])
                    and float(symmlist[inputf][sel]['coverage']) < float(symmlist[inputf]['cesymm_lt']['coverage'])
                    and symmlist[inputf]['cesymm_lt']['internal_symmetry'] == symmlist[inputf][sel]['internal_symmetry']) \
                    or (
                    int(symmlist[inputf][sel]['symmetry_levels']) > int(symmlist[inputf]['cesymm_lt']['symmetry_levels'])
                    and symmlist[inputf][sel]['symmetry_order'] == symmlist[inputf]['cesymm_lt']['symmetry_order']
                    and float(symmlist[inputf][sel]['coverage']) <= float(symmlist[inputf]['cesymm_lt']['coverage'])) \
                    or (int(symmlist[inputf][sel]['symmetry_levels']) == int(symmlist[inputf]['cesymm_lt']['symmetry_levels'])
                        and float(symmlist[inputf][sel]['coverage']) == float(symmlist[inputf]['cesymm_lt']['coverage'])
                        and int(symmlist[inputf][sel]['repeats_number']) < int(symmlist[inputf]['cesymm_lt']['repeats_number'])
                        and symmlist[inputf]['cesymm_lt']['internal_symmetry'] == "No"):
                symmlist[inputf]['selected']['source'] = 'cesymm_lt'
        all_options_write(out2, inputf, symmlist[inputf]['cesymm_lt'])

    # Account for multi-chain cases where CE-Symm declares no symmetry even though it's unrefined TM-score is 1
    elif len(chains) > 1 \
            and (int(symmlist[inputf]['cesymm_lt']['repeats_number']) == 1
                 or int(symmlist[inputf]['cesymm_lt']['repeat_length']) < doubletmhelices) \
            and (float(symmlist[inputf]['cesymm_lt']['unrefined_tmscore']) > 0.90
                 or (float(symmlist[inputf]['symd']['tmscore']) > 0.90
                     and float(symmlist[inputf]['symd']['coverage']) > 0.90)):
        # exmaple: 1yew: A;B;C;E;F;G;I;J;K; angle=120
        # find the angle of symmetry
        k = 0
        if symmlist[inputf]['cesymm_lt']['unit_angle'] != 'na' \
                and len(symmlist[inputf]['cesymm_lt']['unit_angle'][:-1].split(';')) == 1 \
                and float(symmlist[inputf]['cesymm_lt']['unit_angle'].strip(';')) != 0 \
                and float(symmlist[inputf]['cesymm_lt']['unrefined_tmscore']) > 0.90:
            k = round(360 / float(symmlist[inputf]['cesymm_lt']['unit_angle'].strip(';')))
        # divide chains in repeats
        if k != 0 and len(chains) % k == 0 and round(k) != 1:
            n = int(len(chains) / k)  # n=3
            repeats_list = [chains[i:i + n] for i in range(0, len(chains), n)]  # repeats_list=[[A,B,C],[E,F,G],[I,J,K]]
            cs_bug_fix_dic = bug_fixer(muscle_exe, wkdir, inputf, oriented, repeats_list, "_analysis")
            symmlist[inputf]['cesymm_lt_bug_fix'] = cs_bug_fix_dic
            if len(cs_bug_fix_dic) > 0:
                cs_bug_fix_dic['chains'] = chains_str
                cs_bug_fix_dic['source'] = 'cesymm_lt_bug_fix'
                cs_bug_fix_dic['oriented'] = oriented
                symmlist[inputf]['cesymm_lt_bug_fix'] = cs_bug_fix_dic
                symmlist[inputf]['selected']['source'] = 'cesymm_lt_bug_fix'
                all_options_write(out2, inputf, symmlist[inputf]['cesymm_lt_bug_fix'])


    all_names = ['_ordered', '_com_sd', '_sd_order']
    for name in all_names:
        # Handle chain re-ordering
        if os.path.isfile(cesymm_lt_out_dir + inputf + name + "_stdout.out"):
            name_dic = 'cesymm_lt' + name
            symmlist[inputf][name_dic] = {}
            symmlist[inputf][name_dic] = cesymm_data(cesymm_lt_out_dir, inputf + name)
            symmlist[inputf][name_dic]['source'] = name_dic
            oriented_tmp = chain_ordered_dir + inputf + name + ".pdb"
            symmlist[inputf][name_dic]['oriented'] = oriented_tmp
            chains_str_tmp = find_chains(oriented_tmp)
            chains_tmp = chains_str_tmp[:-1].split(';')
            symmlist[inputf][name_dic]['chains'] = chains_str_tmp
            print(chains_tmp)

            # Fix the CE-Symm bug with refinement
            symd_dic = 'symd' + name
            symmlist[inputf][symd_dic] = {}
            symmlist[inputf][symd_dic] = symd_data(symd_out_dir, inputf + name)
            if symmlist[inputf][name_dic]['repeats_number'] != 'na' \
                    and len(chains_tmp) > 1 \
                    and (int(symmlist[inputf][name_dic]['repeats_number']) == 1
                         or int(symmlist[inputf][name_dic]['repeat_length']) < doubletmhelices) \
                    and (float(symmlist[inputf][name_dic]['unrefined_tmscore']) > 0.90
                         or (float(symmlist[inputf][symd_dic]['tmscore']) > 0.90
                             and float(symmlist[inputf][symd_dic]['coverage']) > 0.90)):
                # exmaple: 1yew: A;B;C;E;F;G;I;J;K; angle=120
                # find the angle of symmetry
                k = 0
                if symmlist[inputf][name_dic]['unit_angle'] != 'na' \
                        and len(symmlist[inputf][name_dic]['unit_angle'][:-1].split(';')) == 1 \
                        and float(symmlist[inputf][name_dic]['unit_angle'].strip(';')) != 0 \
                        and float(symmlist[inputf][name_dic]['unrefined_tmscore']) > 0.90:
                    k = round(360 / float(symmlist[inputf][name_dic]['unit_angle'].strip(';')))
                if (k != 0 and len(chains_tmp) % k != 0) or k == 0 or round(k) == 1:
                    k = round(360 / float(symmlist[inputf][symd_dic]['unit_angle']))  # k=3
                    if len(chains_tmp) % k != 0 or round(k) == 1:
                        k = round(360 / float(symmlist[inputf][symd_dic]['is_angle']))
                # divide chains_tmp in repeats
                if len(chains_tmp) % k == 0 and round(k) != 1:
                    n = int(len(chains_tmp) / k)  # n=3
                    repeats_list = [chains_tmp[i:i + n] for i in
                                    range(0, len(chains_tmp), n)]  # repeats_list=[[A,B,C],[E,F,G],[I,J,K]]
                    cs_bug_fix_dic = bug_fixer(muscle_exe, wkdir, inputf + name, oriented_tmp, repeats_list, "_analysis")
                    if len(cs_bug_fix_dic) > 0:
                        cs_bug_fix_dic['chains'] = chains_str_tmp
                        symmlist[inputf][name_dic]['chains'] = chains_str_tmp
                        all_options_write(out2, inputf, symmlist[inputf][name_dic])
                        symmlist[inputf][name_dic] = cs_bug_fix_dic
                        symmlist[inputf][name_dic]['oriented'] = oriented_tmp
                        symmlist[inputf][name_dic]['source'] = name_dic + '_bug_fix'

            # Were the analyzed chains just a subset of the TM chains?
            additions = ""
            for ch in chains:  # K;3;2;J;F;4;1;G;B;R;A;L;H;I; but only B;A; appear in the symmetry order (2wsc)
                if ch not in chains_str_tmp:
                    additions = additions + ch + ";"

            # Re-calculate scores if only a subset of chains was analyzed by CE-Symm
            if len(additions) != 0 and symmlist[inputf][name_dic]['repeats_number'] != 'na' and int(
                    symmlist[inputf][name_dic]['repeats_number']) > 1:
                chains_str_tmp = chains_str_tmp + additions
                chains_tmp = chains_str_tmp[:-1].split(';')
                symmlist[inputf][name_dic]['chains'] = chains_str_tmp
                num = strip_tm_chains_in_order(wkdir, inputf + name, oriented_opm, chains_tmp)
                oriented_tmp = wkdir + inputf + name + "_tmp.pdb"
                print("name_dic to bug fix", name_dic)
                if symmlist[inputf][name_dic]['source'] == name_dic + '_bug_fix':
                    tm_score, rmsd, angle, prot_len = partial_repeats_scores(wkdir, oriented_tmp,
                                                                             symmlist[inputf][name_dic]['files_key'])
                else:
                    tm_score, rmsd, angle, prot_len = partial_repeats_scores(cesymm_lt_out_dir, oriented_tmp, inputf + name)
                #            tm_score, rmsd, angle, prot_len=partial_repeats_scores(cesymm_lt_out_dir, oriented_tmp, inputf+name)
                if abs(rmsd - float(symmlist[inputf][name_dic]['refined_rmsd'])) < 0.2:
                    symmlist[inputf][name_dic]['refined_tmscore'] = tm_score
                    symmlist[inputf][name_dic]['coverage'] = float(
                        '%.2f' % (symmlist[inputf][name_dic]['aligned_length'] / prot_len))
                    symmlist[inputf][name_dic]['size'] = prot_len
                    symmlist[inputf][name_dic]['oriented'] = oriented_tmp

            all_options_write(out2, inputf, symmlist[inputf][name_dic])

            # Compare with the current best
            sel = symmlist[inputf]['selected']['source']
            print("here:", symmlist[inputf]['selected'])
            print(symmlist[inputf][name_dic])
            #   if no symmetry was found before but now we have one, select it
            if (symmlist[inputf][sel]['repeats_number'] == 'na'
                or int(symmlist[inputf][sel]['repeats_number']) == 1
                or (int(symmlist[inputf][sel]['repeat_length']) < doubletmhelices and int(
                        symmlist[inputf][name_dic]['repeat_length']) > doubletmhelices)
                or (symmlist[inputf][sel]['internal_symmetry'] == 'Yes'
                    and symmlist[inputf][name_dic]['internal_symmetry'] == 'No')) \
                    and symmlist[inputf][name_dic]['repeats_number'] != 'na' \
                    and int(symmlist[inputf][name_dic]['repeats_number']) > 1 \
                    and int(symmlist[inputf][name_dic]['repeat_length']) > doubletmhelices \
                    and symmlist[inputf][name_dic]['internal_symmetry'] == 'No':
                symmlist[inputf]['selected']['source'] = name_dic

            #  otherwise, select the current symmetry if either it both has more symmetry levels and a repeat length >doubletmhelices and greater coverage
            # or both has the same number of levels and greater coverage
            # or has fewer levels, same symmetry order (eg. C8) and greater or same coverage (eg. 3k0g)
            # or has the same levels, same coverage, but more repeats and is quaternary (e.g. 3ux4, 1k4c)
            elif symmlist[inputf][sel]['repeats_number'] != 'na' \
                    and int(symmlist[inputf][sel]['repeats_number']) > 1 \
                    and symmlist[inputf][name_dic]['repeats_number'] != 'na' \
                    and int(symmlist[inputf][name_dic]['repeats_number']) > 1:
                if int(symmlist[inputf][sel]['symmetry_levels']) < int(symmlist[inputf][name_dic]['symmetry_levels']) \
                        and symmlist[inputf][sel]['internal_symmetry'] == symmlist[inputf][name_dic]['internal_symmetry'] \
                        and float(symmlist[inputf][sel]['coverage']) <= float(symmlist[inputf][name_dic]['coverage']) \
                        and int(symmlist[inputf][sel]['repeats_number']) == int(
                    symmlist[inputf][name_dic]['repeats_number']):
                    pass
                elif (int(symmlist[inputf][sel]['symmetry_levels']) < int(symmlist[inputf][name_dic]['symmetry_levels'])
                      and (int(symmlist[inputf][name_dic]['repeat_length']) >= doubletmhelices
                           and symmlist[inputf][name_dic]['internal_symmetry'] == symmlist[inputf][sel]['internal_symmetry']
                           and float(symmlist[inputf][sel]['coverage']) <= float(symmlist[inputf][name_dic]['coverage']))) \
                        or (
                        int(symmlist[inputf][sel]['symmetry_levels']) == int(symmlist[inputf][name_dic]['symmetry_levels'])
                        and float(symmlist[inputf][sel]['coverage']) < float(symmlist[inputf][name_dic]['coverage'])
                        and symmlist[inputf][name_dic]['internal_symmetry'] == symmlist[inputf][sel]['internal_symmetry']) \
                        or (
                        int(symmlist[inputf][sel]['symmetry_levels']) > int(symmlist[inputf][name_dic]['symmetry_levels'])
                        and symmlist[inputf][sel]['symmetry_order'] == symmlist[inputf][name_dic]['symmetry_order']
                        and float(symmlist[inputf][sel]['coverage']) <= float(symmlist[inputf][name_dic]['coverage'])) \
                        or (int(symmlist[inputf][sel]['symmetry_levels']) == int(symmlist[inputf][name_dic]['symmetry_levels'])
                            and float(symmlist[inputf][sel]['coverage']) == float(symmlist[inputf][name_dic]['coverage'])
                            and int(symmlist[inputf][sel]['repeats_number']) < int(symmlist[inputf][name_dic]['repeats_number'])
                            and symmlist[inputf][name_dic]['internal_symmetry'] == "No"):
                    symmlist[inputf]['selected']['source'] = name_dic
                    # Clean up
        #    if len(additions)!=0:
        #        os.remove(wkdir+inputf+"_ordered_tmp.pdb")

    # Handle corrections to the number of levels / minimum repeat length AND to the number of levels for quaternary symmetry AND enforced order
    all_names = ['_df', '_ordered_df', '_sd_order_df', '_com_sd_df', '_lt', '_ordered_lt', '_sd_order_lt', '_com_sd_lt']
    modifier_dic = {'_mlev': 'cesymm_minlen', '_quat': 'cesymm_quat_minlen', '_ord': 'cesymm_order'}
    for modifier in modifier_dic:
        for name in all_names:
            cs_out_dir = locations['FSYSPATH'][modifier_dic[modifier]] + inputf + name + "/"
            cs_out_dir_rel = locations['FSYS'][modifier_dic[modifier]] + inputf + name + "/"
            if os.path.isfile(cs_out_dir + inputf + name + modifier + "_stdout.out"):
                name_dic = modifier[1:] + name
                symmlist[inputf][name_dic] = {}
                symmlist[inputf][name_dic] = cesymm_data(cs_out_dir, inputf + name + modifier)
                symmlist[inputf][name_dic]['source'] = name_dic
                if name[:-3] == '':
                    oriented_tmp = oriented
                else:
                    oriented_tmp = chain_ordered_dir + inputf + name[:-3] + ".pdb"
                if 'lt' in name:
                    tag = 'cesymm_lt' + name[:-3]
                else:
                    tag = 'cesymm' + name[:-3]

                symmlist[inputf][name_dic]['oriented'] = oriented_tmp
                chains_str_tmp = find_chains(oriented_tmp)
                chains_tmp = chains_str_tmp[:-1].split(';')
                symmlist[inputf][name_dic]['chains'] = chains_str_tmp
                print(chains_tmp)
                symmlist[inputf][name_dic]['files_key'] = inputf + name + modifier
                symmlist[inputf][name_dic]['files_dir'] = cs_out_dir_rel
                symmlist[inputf][name_dic]['files_dir_path'] = cs_out_dir

                # Were the analyzed chains just a subset of the TM chains?
                additions = ""
                for ch in chains:  # K;3;2;J;F;4;1;G;B;R;A;L;H;I; but only B;A; appear in the symmetry order (2wsc)
                    if ch not in chains_str_tmp:
                        additions = additions + ch + ";"
                # Re-calculate scores if only a subset of chains was analyzed by CE-Symm
                if len(additions) != 0 and symmlist[inputf][name_dic]['repeats_number'] != 'na' and int(
                        symmlist[inputf][name_dic]['repeats_number']) > 1:
                    chains_str_tmp = chains_str_tmp + additions
                    chains_tmp = chains_str_tmp[:-1].split(';')
                    symmlist[inputf][name_dic]['chains'] = chains_str_tmp
                    num = strip_tm_chains_in_order(wkdir, inputf + name, oriented_opm, chains_tmp)
                    oriented_tmp = wkdir + inputf + name + "_tmp.pdb"
                    tm_score, rmsd, angle, prot_len = partial_repeats_scores(cs_out_dir, oriented_tmp,
                                                                             inputf + name + modifier)
                    if abs(rmsd - float(symmlist[inputf][name_dic]['refined_rmsd'])) < 0.2:
                        symmlist[inputf][name_dic]['refined_tmscore'] = tm_score
                        symmlist[inputf][name_dic]['coverage'] = float(
                            '%.2f' % (symmlist[inputf][name_dic]['aligned_length'] / prot_len))
                        symmlist[inputf][name_dic]['size'] = prot_len
                        symmlist[inputf][name_dic]['oriented'] = oriented_tmp

                all_options_write(out2, inputf, symmlist[inputf][name_dic])

                # Compare with the current best
                sel = symmlist[inputf]['selected']['source']
                #   if no symmetry was found before but now we have one that has a repeat with at least doubletmhelices a.a., select it
                if (symmlist[inputf][sel]['repeats_number'] == 'na'
                    or int(symmlist[inputf][sel]['repeats_number']) == 1
                    or (int(symmlist[inputf][sel]['repeat_length']) < doubletmhelices and int(
                            symmlist[inputf][name_dic]['repeat_length']) > doubletmhelices)
                    or (symmlist[inputf][sel]['internal_symmetry'] == 'Yes'
                        and symmlist[inputf][name_dic]['internal_symmetry'] == 'No')) \
                        and symmlist[inputf][name_dic]['repeats_number'] != 'na' \
                        and int(symmlist[inputf][name_dic]['repeats_number']) > 1 \
                        and int(symmlist[inputf][name_dic]['repeat_length']) > doubletmhelices:
                    symmlist[inputf]['selected']['source'] = name_dic

                # otherwise, select the new symmetry if it was a correction for the currently selected one
                #            elif sel==tag and int(symmlist[inputf][name_dic]['repeats_number'])>1:
                #                symmlist[inputf]['selected']['source']=name_dic

                #  otherwise, select the current symmetry if either it both has more symmetry levels and a repeat length >doubletmhelices
                # or both has the same number of levels and greater coverage
                # or has fewer levels, same symmetry order and greater or same coverage (eg. 3k0g)
                # or has the same levels, same coverage, but more repeats and is quaternary (e.g. 3ux4, 1k4c)
                elif symmlist[inputf][sel]['repeats_number'] != 'na' \
                        and int(symmlist[inputf][sel]['repeats_number']) > 1 \
                        and symmlist[inputf][name_dic]['repeats_number'] != 'na' \
                        and int(symmlist[inputf][name_dic]['repeats_number']) > 1:
                    if int(symmlist[inputf][sel]['symmetry_levels']) < int(symmlist[inputf][name_dic]['symmetry_levels']) \
                            and symmlist[inputf][sel]['internal_symmetry'] == symmlist[inputf][name_dic][
                        'internal_symmetry'] \
                            and float(symmlist[inputf][sel]['coverage']) <= float(symmlist[inputf][name_dic]['coverage']) \
                            and int(symmlist[inputf][sel]['repeats_number']) == int(
                        symmlist[inputf][name_dic]['repeats_number']):
                        pass
                    elif (int(symmlist[inputf][sel]['symmetry_levels']) < int(symmlist[inputf][name_dic]['symmetry_levels'])
                          and (int(symmlist[inputf][name_dic]['repeat_length']) >= doubletmhelices
                               and symmlist[inputf][sel]['internal_symmetry'] == symmlist[inputf][name_dic][
                                   'internal_symmetry']
                               and float(symmlist[inputf][sel]['coverage']) <= float(
                                        symmlist[inputf][name_dic]['coverage']))) \
                            or (int(symmlist[inputf][sel]['symmetry_levels']) == int(
                        symmlist[inputf][name_dic]['symmetry_levels'])
                                and float(symmlist[inputf][sel]['coverage']) < float(symmlist[inputf][name_dic]['coverage'])
                                and symmlist[inputf][sel]['internal_symmetry'] == symmlist[inputf][name_dic][
                                    'internal_symmetry']) \
                            or (int(symmlist[inputf][sel]['symmetry_levels']) > int(
                        symmlist[inputf][name_dic]['symmetry_levels'])
                                and symmlist[inputf][sel]['symmetry_order'] == symmlist[inputf][name_dic]['symmetry_order']
                                and float(symmlist[inputf][sel]['coverage']) <= float(
                                symmlist[inputf][name_dic]['coverage'])) \
                            or (int(symmlist[inputf][sel]['symmetry_levels']) == int(symmlist[inputf][name_dic]['symmetry_levels'])
                                and float(symmlist[inputf][sel]['coverage']) == float(symmlist[inputf][name_dic]['coverage'])
                                and int(symmlist[inputf][sel]['repeats_number']) < int(symmlist[inputf][name_dic]['repeats_number'])
                                and symmlist[inputf][name_dic]['internal_symmetry'] == "No"):
                        symmlist[inputf]['selected']['source'] = name_dic

    # Handle the RMSD option
    maxrmsd_options = ['cesymm_rmsd_2_5', 'cesymm_rmsd_3']
    for rmsd_option in maxrmsd_options:
        cesymm_rmsd_out_dir = locations['FSYSPATH'][rmsd_option] + inputf + "/"
        if os.path.isfile(cesymm_rmsd_out_dir + inputf + "_stdout.out"):
            name_dic = rmsd_option[7:]
            symmlist[inputf][name_dic] = {}
            symmlist[inputf][name_dic] = cesymm_data(cesymm_rmsd_out_dir, inputf)
            symmlist[inputf][name_dic]['source'] = name_dic
            symmlist[inputf][name_dic]['oriented'] = oriented
            chains_str_tmp = find_chains(oriented)
            chains_tmp = chains_str_tmp[:-1].split(';')
            symmlist[inputf][name_dic]['chains'] = chains_str_tmp
            symmlist[inputf][name_dic]['files_key'] = inputf
            symmlist[inputf][name_dic]['files_dir'] = locations['FSYS'][rmsd_option] + inputf + "/"
            symmlist[inputf][name_dic]['files_dir_path'] = locations['FSYSPATH'][rmsd_option]  + inputf + "/"
            all_options_write(out2, inputf, symmlist[inputf][name_dic])


            # Compare with the current best
            sel = symmlist[inputf]['selected']['source']
            #   if no symmetry was found before but now we have one that has a repeat with at least doubletmhelices a.a., select it
            if (symmlist[inputf][sel]['repeats_number'] == 'na'
                or int(symmlist[inputf][sel]['repeats_number']) == 1
                or (int(symmlist[inputf][sel]['repeat_length']) < doubletmhelices
                    and int(symmlist[inputf][name_dic]['repeat_length']) > doubletmhelices)
                or (symmlist[inputf][sel]['internal_symmetry'] == 'Yes'
                    and symmlist[inputf][name_dic]['internal_symmetry'] == 'No')) \
                    and symmlist[inputf][name_dic]['repeats_number'] != 'na' \
                    and int(symmlist[inputf][name_dic]['repeats_number']) > 1 \
                    and int(symmlist[inputf][name_dic]['repeat_length']) > doubletmhelices:
                symmlist[inputf]['selected']['source'] = name_dic

            #  otherwise, select the current symmetry if either it both has more symmetry levels and a repeat length > doubletmhelices
            # or both has the same number of levels and greater coverage
            #  or has fewer levels, same symmetry order and greater or same coverage (eg. 3k0g)
            elif symmlist[inputf][sel]['repeats_number'] != 'na' \
                    and int(symmlist[inputf][sel]['repeats_number']) > 1 \
                    and symmlist[inputf][name_dic]['repeats_number'] != 'na' \
                    and int(symmlist[inputf][name_dic]['repeats_number']) > 1:
                if int(symmlist[inputf][sel]['symmetry_levels']) < int(symmlist[inputf][name_dic]['symmetry_levels']) \
                        and symmlist[inputf][sel]['internal_symmetry'] == symmlist[inputf][name_dic]['internal_symmetry'] \
                        and float(symmlist[inputf][sel]['coverage']) <= float(symmlist[inputf][name_dic]['coverage']) \
                        and int(symmlist[inputf][sel]['repeats_number']) == int(
                    symmlist[inputf][name_dic]['repeats_number']):
                    pass
                elif (int(symmlist[inputf][sel]['symmetry_levels']) < int(symmlist[inputf][name_dic]['symmetry_levels'])
                      and (int(symmlist[inputf][name_dic]['repeat_length']) >= doubletmhelices
                           and symmlist[inputf][sel]['internal_symmetry'] == symmlist[inputf][name_dic]['internal_symmetry']
                           and float(symmlist[inputf][sel]['coverage']) <= float(symmlist[inputf][name_dic]['coverage']))) \
                        or (int(symmlist[inputf][sel]['symmetry_levels']) == int(symmlist[inputf][name_dic]['symmetry_levels'])
                            and float(symmlist[inputf][sel]['coverage']) < float(symmlist[inputf][name_dic]['coverage'])
                            and symmlist[inputf][sel]['internal_symmetry'] == symmlist[inputf][name_dic]['internal_symmetry']) \
                        or (int(symmlist[inputf][sel]['symmetry_levels']) > int(symmlist[inputf][name_dic]['symmetry_levels'])
                            and symmlist[inputf][sel]['symmetry_order'] == symmlist[inputf][name_dic]['symmetry_order']
                            and float(symmlist[inputf][sel]['coverage']) <= float(symmlist[inputf][name_dic]['coverage'])) \
                        or (int(symmlist[inputf][sel]['symmetry_levels']) == int(symmlist[inputf][name_dic]['symmetry_levels'])
                            and float(symmlist[inputf][sel]['coverage']) == float(symmlist[inputf][name_dic]['coverage'])
                            and int(symmlist[inputf][sel]['repeats_number']) < int(symmlist[inputf][name_dic]['repeats_number'])
                            and symmlist[inputf][name_dic]['internal_symmetry'] == "No"):
                    symmlist[inputf]['selected']['source'] = name_dic

    # Handle QuatSymm
    ## Load dictionary of QuatSymm results if available
    if os.path.isfile(locations['FSYSPATH']['quatsymm'] + inputf + "/" + inputf + "_qsymm_enc.pkl"):
        symmlist[inputf]['quatsymm'] = pkl.load(open(locations['FSYSPATH']['quatsymm'] + inputf + "/" + inputf + "_qsymm_enc.pkl",'rb'))
        # symmlist[inputf]['quatsymm']['repeats_number'] = symmlist[inputf]['quatsymm']['repeats_num']
        all_options_write(out2, inputf, symmlist[inputf]['quatsymm'])
        # Compare with current best
        sel = symmlist[inputf]['selected']['source']
        ## We assume that since sequence identity is considered in QuatSymm, single TM repeats represent
        ## functionally-relevant symmetries; we also assume that since only TM chains are considered, we
        ## don't need to check for the number of transmembrane segments.
        dtm = float(symmlist[inputf]['quatsymm']['refined_tmscore']) - float(symmlist[inputf]['quatsymm']['raw_data']['tmscore'])
        drmsd = float(symmlist[inputf]['quatsymm']['refined_rmsd']) - float(symmlist[inputf]['quatsymm']['raw_data']['rmsd'])
        print("check1:", (np.abs(dtm) <= 0.1 or np.abs(drmsd) <= 1))
        print("check2:", int(symmlist[inputf][sel]['symmetry_levels']) >= int(symmlist[inputf]['quatsymm']['symmetry_levels']))
        print("check3:", float(symmlist[inputf][sel]['coverage']) <= float(symmlist[inputf]['quatsymm']['coverage']))
        if (np.abs(dtm) <= 0.1 or np.abs(drmsd) <= 1) \
                and ((symmlist[inputf][sel]['repeats_number'] == 'na' or
                      int(symmlist[inputf][sel]['repeats_number']) == 1 or
                      int(symmlist[inputf][sel]['repeat_length']) < doubletmhelices) or #and int(symmlist[inputf]['quatsymm']['repeat_length']) > doubletmhelices)) or
                     (int(symmlist[inputf][sel]['symmetry_levels']) >= int(symmlist[inputf]['quatsymm']['symmetry_levels'])
                      and float(symmlist[inputf][sel]['coverage']) <= float(symmlist[inputf]['quatsymm']['coverage']))
                     or (int(symmlist[inputf][sel]['symmetry_levels']) == int(symmlist[inputf]['quatsymm']['symmetry_levels'])
                         and float(symmlist[inputf][sel]['coverage']) == float(symmlist[inputf]['quatsymm']['coverage'])
                         and int(symmlist[inputf][sel]['repeats_number']) < int(symmlist[inputf]['quatsymm']['repeats_number']))
                     or symmlist[inputf][sel]['internal_symmetry'] == 'Yes'):
            symmlist[inputf]['selected']['source'] = 'quatsymm'
            print("Quatsymm selected.")

    sel = symmlist[inputf]['selected']['source']

    # Image creation

    image_files = ""
    pml_files = ""
    jmol_files = ""
    super_pml_files = ""
    tm_score_super = "na"

    if 'cesymm_lt' in symmlist[inputf][sel]['source']:
        if 'bug_fix' in symmlist[inputf][sel]['source']:
            reference = parse_structure(symmlist[inputf][sel]['oriented'])[0]
            repeat_selection = symmlist[inputf][sel]['repeat_selection']
            files_key = symmlist[inputf][sel]['files_key']
            symmlist[inputf][sel]['files_dir'] = locations['FSYS']['mssd'] + inputf + "/"
            symmlist[inputf][sel]['files_dir_path'] = locations['FSYSPATH']['mssd'] + inputf + "/"
            image_key = inputf + "_analysis"
            if int(symmlist[inputf][sel]['aligned_length']) >= 2 * doubletmhelices:
                if want_image == 1:
                    symmlist[inputf][sel]['image_files'], \
                    symmlist[inputf][sel]['pml_files'], \
                    symmlist[inputf][sel]['jmol_files'] = cesymm_images(wkdir, wkdir, files_key, repeat_selection,
                                                                       symmlist[inputf][sel]['oriented'], reference, ray,
                                                                       pymol, pic_dir_whole, jmol_dir_whole, image_key)
                    symmlist[inputf][sel]['super_pml_files'], tm_score_super = repeats_superposition(wkdir, symmlist[inputf][sel]['oriented'],
                                                                        symmlist[inputf][sel]['files_key'], super_dir_whole,
                                                                        image_key,
                                                                        symmlist[inputf][sel]['repeat_selection'])
                    symmlist[inputf][sel]['image_key'] = image_key
                os.rename(wkdir + files_key + ".axes", selected_out + image_key + ".axes")
                os.rename(wkdir + files_key + ".fasta", selected_out + image_key + ".fasta")
                symmlist[inputf][sel]['files_key'] = image_key
        elif 'ordered' in symmlist[inputf][sel]['source']:
            reference = parse_structure(symmlist[inputf][sel]['oriented'])[0]
            repeat_selection = get_repeat_resid(cesymm_lt_out_dir, inputf + "_ordered", reference)
            symmlist[inputf][sel]['repeat_selection'] = repeat_selection
            symmlist[inputf][sel]['files_key'] = inputf + "_ordered"
            symmlist[inputf][sel]['files_dir'] = locations['FSYS']['cesymm_low_thr'] + inputf + "/"
            symmlist[inputf][sel]['files_dir_path'] = locations['FSYSPATH']['cesymm_low_thr'] + inputf + "/"
            image_key = inputf + "_analysis"
            if int(symmlist[inputf][sel]['aligned_length']) >= 2 * doubletmhelices and want_image == 1:
                symmlist[inputf][sel]['image_files'], \
                symmlist[inputf][sel]['pml_files'], \
                symmlist[inputf][sel]['jmol_files'] = cesymm_images(wkdir, cesymm_lt_out_dir, inputf + "_ordered",
                                                                   repeat_selection, symmlist[inputf][sel]['oriented'],
                                                                   reference, ray, pymol, pic_dir_whole, jmol_dir_whole,
                                                                   image_key)
                symmlist[inputf][sel]['super_pml_files'], tm_score_super = repeats_superposition(symmlist[inputf][sel]['files_dir_path'],
                                                                        symmlist[inputf][sel]['oriented'],
                                                                        symmlist[inputf][sel]['files_key'], super_dir_whole,
                                                                        image_key,
                                                                        symmlist[inputf][sel]['repeat_selection'])
                symmlist[inputf][sel]['image_key'] = image_key

            # remember to remove any path indication from files and to record downloadable files
            # include pointer to normal cesymm files
        elif 'sd_order' in symmlist[inputf][sel]['source']:
            reference = parse_structure(symmlist[inputf][sel]['oriented'])[0]
            repeat_selection = get_repeat_resid(cesymm_lt_out_dir, inputf + "_sd_order", reference)
            symmlist[inputf][sel]['repeat_selection'] = repeat_selection
            symmlist[inputf][sel]['files_key'] = inputf + "_sd_order"
            symmlist[inputf][sel]['files_dir'] = locations['FSYS']['cesymm_low_thr'] + inputf + "/"
            symmlist[inputf][sel]['files_dir_path'] = locations['FSYSPATH']['cesymm_low_thr'] + inputf + "/"
            image_key = inputf + "_analysis"
            if int(symmlist[inputf][sel]['aligned_length']) >= 2 * doubletmhelices and want_image == 1:
                symmlist[inputf][sel]['image_files'], \
                symmlist[inputf][sel]['pml_files'], \
                symmlist[inputf][sel]['jmol_files'] = cesymm_images(wkdir, cesymm_lt_out_dir, inputf + "_sd_order",
                                                                   repeat_selection, symmlist[inputf][sel]['oriented'],
                                                                   reference, ray, pymol, pic_dir_whole, jmol_dir_whole,
                                                                   image_key)
                symmlist[inputf][sel]['super_pml_files'], tm_score_super = repeats_superposition(symmlist[inputf][sel]['files_dir_path'],
                                                                        symmlist[inputf][sel]['oriented'],
                                                                        symmlist[inputf][sel]['files_key'], super_dir_whole,
                                                                        image_key,
                                                                        symmlist[inputf][sel]['repeat_selection'])
                symmlist[inputf][sel]['image_key'] = image_key

        elif 'com_sd' in symmlist[inputf][sel]['source']:
            reference = parse_structure(symmlist[inputf][sel]['oriented'])[0]
            repeat_selection = get_repeat_resid(cesymm_lt_out_dir, inputf + "_com_sd", reference)
            symmlist[inputf][sel]['repeat_selection'] = repeat_selection
            symmlist[inputf][sel]['files_key'] = inputf + "_com_sd"
            symmlist[inputf][sel]['files_dir'] = locations['FSYS']['cesymm_low_thr'] + inputf + "/"
            symmlist[inputf][sel]['files_dir_path'] = locations['FSYSPATH']['cesymm_low_thr'] + inputf + "/"
            image_key = inputf + "_analysis"
            if int(symmlist[inputf][sel]['aligned_length']) >= 2 * doubletmhelices and want_image == 1:
                symmlist[inputf][sel]['image_files'], \
                symmlist[inputf][sel]['pml_files'], \
                symmlist[inputf][sel]['jmol_files'] = cesymm_images(wkdir, cesymm_lt_out_dir, inputf + "_com_sd",
                                                                   repeat_selection, symmlist[inputf][sel]['oriented'],
                                                                   reference, ray, pymol, pic_dir_whole, jmol_dir_whole,
                                                                   image_key)
                symmlist[inputf][sel]['super_pml_files'], tm_score_super = repeats_superposition(symmlist[inputf][sel]['files_dir_path'],
                                                                        symmlist[inputf][sel]['oriented'],
                                                                        symmlist[inputf][sel]['files_key'], super_dir_whole,
                                                                        image_key,
                                                                        symmlist[inputf][sel]['repeat_selection'])
                symmlist[inputf][sel]['image_key'] = image_key
        else:
            reference = parse_structure(symmlist[inputf][sel]['oriented'])[0]
            repeat_selection = get_repeat_resid(cesymm_lt_out_dir, inputf, reference)
            symmlist[inputf][sel]['repeat_selection'] = repeat_selection
            symmlist[inputf][sel]['files_key'] = inputf
            symmlist[inputf][sel]['files_dir'] = locations['FSYS']['cesymm_low_thr'] + inputf + "/"
            symmlist[inputf][sel]['files_dir_path'] = locations['FSYSPATH']['cesymm_low_thr'] + inputf + "/"
            image_key = inputf + "_analysis"
            if int(symmlist[inputf][sel]['aligned_length']) >= 2 * doubletmhelices and want_image == 1:
                symmlist[inputf][sel]['image_files'], \
                symmlist[inputf][sel]['pml_files'], \
                symmlist[inputf][sel]['jmol_files'] = cesymm_images(wkdir, cesymm_lt_out_dir, inputf, repeat_selection,
                                                                   symmlist[inputf][sel]['oriented'], reference, ray, pymol,
                                                                   pic_dir_whole, jmol_dir_whole, image_key)
                symmlist[inputf][sel]['super_pml_files'], tm_score_super = repeats_superposition(symmlist[inputf][sel]['files_dir_path'],
                                                                        symmlist[inputf][sel]['oriented'],
                                                                        symmlist[inputf][sel]['files_key'], super_dir_whole,
                                                                        image_key,
                                                                        symmlist[inputf][sel]['repeat_selection'])
                symmlist[inputf][sel]['image_key'] = image_key


    elif 'cesymm' in symmlist[inputf][sel]['source']:
        if 'bug_fix' in symmlist[inputf][sel]['source']:
            reference = parse_structure(symmlist[inputf][sel]['oriented'])[0]
            repeat_selection = symmlist[inputf][sel]['repeat_selection']
            files_key = symmlist[inputf][sel]['files_key']
            symmlist[inputf][sel]['files_dir'] = locations['FSYS']['mssd']  + inputf + "/"
            symmlist[inputf][sel]['files_dir_path'] = locations['FSYSPATH']['mssd'] + inputf + "/"
            image_key = inputf + "_analysis"
            if int(symmlist[inputf][sel]['aligned_length']) >= 2 * doubletmhelices:
                if want_image == 1:
                    symmlist[inputf][sel]['image_files'], \
                    symmlist[inputf][sel]['pml_files'], \
                    symmlist[inputf][sel]['jmol_files'] = cesymm_images(wkdir, wkdir, files_key, repeat_selection,
                                                                       symmlist[inputf][sel]['oriented'], reference, ray,
                                                                       pymol, pic_dir_whole, jmol_dir_whole, image_key)
                    symmlist[inputf][sel]['super_pml_files'], tm_score_super = repeats_superposition(wkdir, symmlist[inputf][sel]['oriented'],
                                                                        symmlist[inputf][sel]['files_key'], super_dir_whole,
                                                                        image_key,
                                                                        symmlist[inputf][sel]['repeat_selection'])
                    symmlist[inputf][sel]['image_key'] = image_key
                os.rename(wkdir + files_key + ".axes", selected_out + image_key + ".axes")
                os.rename(wkdir + files_key + ".fasta", selected_out + image_key + ".fasta")
                symmlist[inputf][sel]['files_key'] = image_key
        elif 'ordered' in symmlist[inputf][sel]['source']:
            reference = parse_structure(symmlist[inputf][sel]['oriented'])[0]
            repeat_selection = get_repeat_resid(cesymm_out_dir, inputf + "_ordered", reference)
            symmlist[inputf][sel]['repeat_selection'] = repeat_selection
            symmlist[inputf][sel]['files_key'] = inputf + "_ordered"
            symmlist[inputf][sel]['files_dir'] = locations['FSYS']['cesymm'] + inputf + "/"
            symmlist[inputf][sel]['files_dir_path'] = locations['FSYSPATH']['cesymm'] + inputf + "/"
            image_key = inputf + "_analysis"
            if int(symmlist[inputf][sel]['aligned_length']) >= 2 * doubletmhelices and want_image == 1:
                symmlist[inputf][sel]['image_files'], \
                symmlist[inputf][sel]['pml_files'], \
                symmlist[inputf][sel]['jmol_files'] = cesymm_images(wkdir, cesymm_out_dir, inputf + "_ordered",
                                                                   repeat_selection, symmlist[inputf][sel]['oriented'],
                                                                   reference, ray, pymol, pic_dir_whole, jmol_dir_whole,
                                                                   image_key)
                symmlist[inputf][sel]['super_pml_files'], tm_score_super = repeats_superposition(symmlist[inputf][sel]['files_dir_path'],
                                                                        symmlist[inputf][sel]['oriented'],
                                                                        symmlist[inputf][sel]['files_key'], super_dir_whole,
                                                                        image_key,
                                                                        symmlist[inputf][sel]['repeat_selection'])
                symmlist[inputf][sel]['image_key'] = image_key

        elif 'sd_order' in symmlist[inputf][sel]['source']:
            reference = parse_structure(symmlist[inputf][sel]['oriented'])[0]
            repeat_selection = get_repeat_resid(cesymm_out_dir, inputf + "_sd_order", reference)
            symmlist[inputf][sel]['repeat_selection'] = repeat_selection
            symmlist[inputf][sel]['files_key'] = inputf + "_sd_order"
            symmlist[inputf][sel]['files_dir'] = locations['FSYS']['cesymm'] + inputf + "/"
            symmlist[inputf][sel]['files_dir_path'] = locations['FSYSPATH']['cesymm'] + inputf + "/"
            image_key = inputf + "_analysis"
            if int(symmlist[inputf][sel]['aligned_length']) >= 2 * doubletmhelices and want_image == 1:
                symmlist[inputf][sel]['image_files'], \
                symmlist[inputf][sel]['pml_files'], \
                symmlist[inputf][sel]['jmol_files'] = cesymm_images(wkdir, cesymm_out_dir, inputf + "_sd_order",
                                                                   repeat_selection, symmlist[inputf][sel]['oriented'],
                                                                   reference, ray, pymol, pic_dir_whole, jmol_dir_whole,
                                                                   image_key)
                symmlist[inputf][sel]['super_pml_files'], tm_score_super = repeats_superposition(symmlist[inputf][sel]['files_dir_path'],
                                                                        symmlist[inputf][sel]['oriented'],
                                                                        symmlist[inputf][sel]['files_key'], super_dir_whole,
                                                                        image_key,
                                                                        symmlist[inputf][sel]['repeat_selection'])
                symmlist[inputf][sel]['image_key'] = image_key

        elif 'com_sd' in symmlist[inputf][sel]['source']:
            symmlist[inputf][sel]['files_key'] = inputf + "_com_sd"
            symmlist[inputf][sel]['files_dir'] = locations['FSYS']['cesymm'] + inputf + "/"
            symmlist[inputf][sel]['files_dir_path'] = locations['FSYSPATH']['cesymm'] + inputf + "/"
            reference = parse_structure(symmlist[inputf][sel]['oriented'])[0]
            repeat_selection = get_repeat_resid(symmlist[inputf][sel]['files_dir_path'], inputf + "_com_sd", reference)
            symmlist[inputf][sel]['repeat_selection'] = repeat_selection
            image_key = inputf + "_analysis"
            if int(symmlist[inputf][sel]['aligned_length']) >= 2 * doubletmhelices and want_image == 1:
                symmlist[inputf][sel]['image_files'], \
                symmlist[inputf][sel]['pml_files'], \
                symmlist[inputf][sel]['jmol_files'] = cesymm_images(wkdir, symmlist[inputf][sel]['files_dir'], inputf + "_com_sd",
                                                                   repeat_selection, symmlist[inputf][sel]['oriented'],
                                                                   reference, ray, pymol, pic_dir_whole, jmol_dir_whole,
                                                                   image_key)
                symmlist[inputf][sel]['super_pml_files'], tm_score_super = repeats_superposition(symmlist[inputf][sel]['files_dir_path'],
                                                                        symmlist[inputf][sel]['oriented'],
                                                                        symmlist[inputf][sel]['files_key'], super_dir_whole,
                                                                        image_key,
                                                                        symmlist[inputf][sel]['repeat_selection'])
                symmlist[inputf][sel]['image_key'] = image_key
        else:
            symmlist[inputf][sel]['files_key'] = inputf
            symmlist[inputf][sel]['files_dir'] = locations['FSYS']['cesymm'] + inputf + "/"
            symmlist[inputf][sel]['files_dir_path'] = locations['FSYSPATH']['cesymm'] + inputf + "/"
            reference = parse_structure(symmlist[inputf][sel]['oriented'])[0]
            repeat_selection = get_repeat_resid(symmlist[inputf][sel]['files_dir_path'], inputf, reference)
            symmlist[inputf][sel]['repeat_selection'] = repeat_selection
            if symmlist[inputf][sel]['aligned_length'] != "na" and int(symmlist[inputf][sel]['aligned_length']) >= 2 * doubletmhelices:
                symmlist[inputf][sel]['image_key'] = inputf + "_cesymm"
                lst = glob.glob(locations['FSYSPATH']['symm_wholepngs'] + inputf + "_cesymm_*.png")
                # lst=sorted(lst, key=lambda name: int(name[len(locations['FSYSPATH']['symm_wholepngs']+inputf+'_cesymm_'):-4]))
                for f in lst:
                    image_files = image_files + f[len(locations['FSYSPATH']['symm_wholepngs']):] + ";"
                    force_symlink(os.path.relpath(f, locations['FSYSPATH']['analysis_wholepngs']), locations['FSYSPATH']['analysis_wholepngs'] + f[len(locations['FSYSPATH']['symm_wholepngs']):])
                    print('Creating link {} -> {}'.format(locations['FSYSPATH']['analysis_wholepngs'] + f[len(locations['FSYSPATH']['symm_wholepngs']):], os.path.relpath(f, locations['FSYSPATH']['analysis_wholepngs'])))
                symmlist[inputf][sel]['image_files'] = image_files
                lst = glob.glob(locations['FSYSPATH']['symm_wholepngs'] + 'pymol_script_' + inputf + "_cesymm_*.pml")
                # lst=sorted(lst, key=lambda name: int(name[len(locations['FSYSPATH']['symm_wholepngs']+'pymol_script_'+inputf+'_cesymm_'):-4]))
                for f in lst:
                    pml_files = pml_files + f[len(locations['FSYSPATH']['symm_wholepngs']):] + ";"
                symmlist[inputf][sel]['pml_files'] = pml_files
                lst = glob.glob(locations['FSYSPATH']['symm_wholejsons'] + '3dmol_' + inputf + "_cesymm_*.json")
                # lst=sorted(lst, key=lambda name: int(name[len(locations['FSYSPATH']['symm_wholepngs']+'pymol_script_'+inputf+'_cesymm_'):-4]))
                for f in lst:
                    jmol_files = jmol_files + f[len(locations['FSYSPATH']['symm_wholejsons']):] + ";"
                symmlist[inputf][sel]['jmol_files'] = jmol_files
                if len(symmlist[inputf][sel]['repeat_selection']) > 1:
                    symmlist[inputf][sel]['super_pml_files'], tm_score_super = repeats_superposition(symmlist[inputf][sel]['files_dir_path'],
                                                                            symmlist[inputf][sel]['oriented'],
                                                                            symmlist[inputf][sel]['files_key'],
                                                                            super_dir_whole, inputf,
                                                                            symmlist[inputf][sel]['repeat_selection'])

                else:
                    symmlist[inputf][sel]['super_pml_files'], tm_score_super = '', 'na'

    elif 'mlev' in symmlist[inputf][sel]['source']:
        reference = parse_structure(symmlist[inputf][sel]['oriented'])[0]
        repeat_selection = get_repeat_resid(symmlist[inputf][sel]['files_dir_path'], symmlist[inputf][sel]['files_key'], reference)
        symmlist[inputf][sel]['repeat_selection'] = repeat_selection
        image_key = inputf + "_analysis"
        if int(symmlist[inputf][sel]['aligned_length']) >= 2 * doubletmhelices and want_image == 1:
            symmlist[inputf][sel]['image_files'], \
            symmlist[inputf][sel]['pml_files'], \
            symmlist[inputf][sel]['jmol_files'] = cesymm_images(wkdir, symmlist[inputf][sel]['files_dir_path'],
                                                               symmlist[inputf][sel]['files_key'], repeat_selection,
                                                               symmlist[inputf][sel]['oriented'], reference, ray, pymol,
                                                               pic_dir_whole, jmol_dir_whole, image_key)
            symmlist[inputf][sel]['super_pml_files'], tm_score_super = repeats_superposition(symmlist[inputf][sel]['files_dir_path'],
                                                                    symmlist[inputf][sel]['oriented'],
                                                                    symmlist[inputf][sel]['files_key'], super_dir_whole,
                                                                    image_key, symmlist[inputf][sel]['repeat_selection'])
            symmlist[inputf][sel]['image_key'] = image_key

    elif 'quatsymm' in symmlist[inputf][sel]['source']:
        print("Moving quatsymm images to webarchive.")
        for f in symmlist[inputf]['quatsymm']['image_files'].strip(";").split(";"):
            shutil.copyfile(locations['FSYSPATH']['quatsymm'] + "pngs/" + f, locations['FSYSPATH']['analysis_wholepngs'] + f)
        for f in symmlist[inputf]['quatsymm']['pml_files'].strip(";").split(";"):
            shutil.copyfile(locations['FSYSPATH']['quatsymm'] + "pngs/" + f, locations['FSYSPATH']['analysis_wholepngs'] + f)
        for f in symmlist[inputf]['quatsymm']['jmol_files'].strip(";").split(";"):
            shutil.copyfile(locations['FSYSPATH']['quatsymm'] + "jsons/" + f, locations['FSYSPATH']['analysis_wholejsons'] + f)
        for f in symmlist[inputf]['quatsymm']['super_pml_files'].strip(";").split(";"):
            shutil.copyfile(locations['FSYSPATH']['quatsymm'] + "super/" + f, locations['FSYSPATH']['analysis_wholesuper'] + f)

    elif 'quat' in symmlist[inputf][sel]['source']:
        reference = parse_structure(symmlist[inputf][sel]['oriented'])[0]
        repeat_selection = get_repeat_resid(symmlist[inputf][sel]['files_dir_path'], symmlist[inputf][sel]['files_key'], reference)
        symmlist[inputf][sel]['repeat_selection'] = repeat_selection
        image_key = inputf + "_analysis"
        if int(symmlist[inputf][sel]['aligned_length']) >= 2 * doubletmhelices and want_image == 1:
            symmlist[inputf][sel]['image_files'], \
            symmlist[inputf][sel]['pml_files'], \
            symmlist[inputf][sel]['jmol_files'] = cesymm_images(wkdir, symmlist[inputf][sel]['files_dir_path'],
                                                               symmlist[inputf][sel]['files_key'], repeat_selection,
                                                               symmlist[inputf][sel]['oriented'], reference, ray, pymol,
                                                               pic_dir_whole, jmol_dir_whole, image_key)
            symmlist[inputf][sel]['super_pml_files'], tm_score_super = repeats_superposition(symmlist[inputf][sel]['files_dir_path'],
                                                                    symmlist[inputf][sel]['oriented'],
                                                                    symmlist[inputf][sel]['files_key'], super_dir_whole,
                                                                    image_key, symmlist[inputf][sel]['repeat_selection'])
            symmlist[inputf][sel]['image_key'] = image_key

    elif 'ord' in symmlist[inputf][sel]['source']:
        reference = parse_structure(symmlist[inputf][sel]['oriented'])[0]
        repeat_selection = get_repeat_resid(symmlist[inputf][sel]['files_dir_path'], symmlist[inputf][sel]['files_key'],
                                            reference)
        symmlist[inputf][sel]['repeat_selection'] = repeat_selection
        image_key = inputf + "_analysis"
        if int(symmlist[inputf][sel]['aligned_length']) >= 2 * doubletmhelices and want_image == 1:
            symmlist[inputf][sel]['image_files'], \
            symmlist[inputf][sel]['pml_files'], \
            symmlist[inputf][sel]['jmol_files'] = cesymm_images(wkdir, symmlist[inputf][sel]['files_dir_path'],
                                                               symmlist[inputf][sel]['files_key'], repeat_selection,
                                                               symmlist[inputf][sel]['oriented'], reference, ray, pymol,
                                                               pic_dir_whole, jmol_dir_whole, image_key)
            symmlist[inputf][sel]['super_pml_files'], tm_score_super = repeats_superposition(symmlist[inputf][sel]['files_dir_path'],
                                                                    symmlist[inputf][sel]['oriented'],
                                                                    symmlist[inputf][sel]['files_key'], super_dir_whole,
                                                                    image_key, symmlist[inputf][sel]['repeat_selection'])
            symmlist[inputf][sel]['image_key'] = image_key

    elif 'rmsd' in symmlist[inputf][sel]['source']:
        reference = parse_structure(symmlist[inputf][sel]['oriented'])[0]
        repeat_selection = get_repeat_resid(symmlist[inputf][sel]['files_dir_path'], symmlist[inputf][sel]['files_key'],
                                            reference)
        symmlist[inputf][sel]['repeat_selection'] = repeat_selection
        image_key = inputf + "_analysis"
        if int(symmlist[inputf][sel]['aligned_length']) >= 2 * doubletmhelices and want_image == 1:
            symmlist[inputf][sel]['image_files'], \
            symmlist[inputf][sel]['pml_files'], \
            symmlist[inputf][sel]['jmol_files'] = cesymm_images(wkdir, symmlist[inputf][sel]['files_dir_path'],
                                                               symmlist[inputf][sel]['files_key'], repeat_selection,
                                                               symmlist[inputf][sel]['oriented'], reference, ray, pymol,
                                                               pic_dir_whole, jmol_dir_whole, image_key)
            symmlist[inputf][sel]['super_pml_files'], tm_score_super = repeats_superposition(symmlist[inputf][sel]['files_dir_path'],
                                                                    symmlist[inputf][sel]['oriented'],
                                                                    symmlist[inputf][sel]['files_key'], super_dir_whole,
                                                                    image_key, symmlist[inputf][sel]['repeat_selection'])
            symmlist[inputf][sel]['image_key'] = image_key

    ### Is there information in the fragments?
    if len(chains) == 1:
        symmlist, fragments, tm_score_super = fragments_check(inputf, chains[0], '', sel, tm_archive, cesymm_data, locations, symmlist, oriented, doubletmhelices, ray, pymol, want_image, pic_dir_whole, jmol_dir_whole, super_dir_whole, out2, image_files, pml_files, jmol_files, super_pml_files, tm_score_super)

    sel = symmlist[inputf]['selected']['source']

    overall_sel = sel
    print("Selection: ", sel)
    if sel != 'cesymm' and cesymm_quat_symm == 0:
        out3.write("enrichment;")
    elif sel == 'cesymm' and discard_small == 1:
        out3.write("discard_small;")
    elif sel == 'cesymm' and discard_nonmemb == 1:
        out3.write("discard_nonmemb;")
    elif sel != 'cesymm' and cesymm_quat_symm == 1:
        out3.write("precision;")

    decision = "Accepted"
    internal_sym = "No"


    # Check if there is at least one repeat in the membrane
    if (';' not in sel) and (
            symmlist[inputf][sel]['repeats_number'] == 'na' or int(symmlist[inputf][sel]['repeats_number']) == 1):
        doms, crossings = 'out', []
        symmlist[inputf][sel]['rep_has_tm_domains'] = doms
        symmlist[inputf][sel]['rep_tm_domain_num'] = crossings
        internal_sym = symmlist[inputf][sel]['internal_symmetry']

    elif (';' not in sel):
        doms, crossings = symmetry_tm_domains(tm_archive, symmlist[inputf][sel]['repeat_selection'], inputf, '')
        symmlist[inputf][sel]['rep_has_tm_domains'] = doms
        symmlist[inputf][sel]['rep_tm_domain_num'] = crossings
        internal_sym = symmlist[inputf][sel]['internal_symmetry']
    else:
        source = sel[:-1].split(';')
        if symmlist[inputf][source[0]]['internal_symmetry'] == 'Yes' or symmlist[inputf][source[1]]['internal_symmetry'] == 'Yes':
            internal_sym = "Yes"
            decision = "Rejected"
        selected_write(out1, inputf, symmlist[inputf][source[0]], decision)
        selected_write(out1, inputf, symmlist[inputf][source[1]], decision)
        if len(chains) == 1 and internal_sym == "Yes":
            selected_write(out1, inputf + "_" + str(chains[0]), symmlist[inputf][source[0]], "Accepted")
            selected_write(out1, inputf + "_" + str(chains[0]), symmlist[inputf][source[1]], "Accepted")


    # Record Final Info Decisions
    if internal_sym == 'Yes':
        symmlist[inputf]['selected']['quatinfo'] = "No quaternary symmetry found in the membrane-bound region during analysis. For internal symmetry, please check the specific chain page."
        decision = "Rejected"

    if (';' not in sel) and (symmlist[inputf][sel]['repeats_number'] == 'na' or int(symmlist[inputf][sel]['repeats_number']) == 1
                             or (int(symmlist[inputf][sel]['aligned_length']) < 2 * doubletmhelices and sel !="quatsymm") or (doms == 'out' and sum(crossings) < 3 and sel != "quatsymm")):
        symmlist[inputf]['selected']['info'] = "No symmetry found in the membrane-bound region during analysis."
        symmlist[inputf]['selected']['quatinfo'] = "No quaternary symmetry found in the membrane-bound region during analysis. For internal symmetry, please check the specific chain page."
        for redundant in glob.glob(locations['FSYSPATH']['analysis_wholepngs'] + inputf + '_analysis*.png'):
            if os.path.islink(redundant):
                os.unlink(redundant)
            else:
                os.remove(redundant)
        for redundant in glob.glob(
            locations['FSYSPATH']['analysis_wholepngs'] + 'pymol_script_' + inputf + '_analysis*.pml'):
            os.remove(redundant)
        decision = "Rejected"
        selected_write(out1, inputf, symmlist[inputf][sel], decision)
    elif (';' not in sel):
        selected_write(out1, inputf, symmlist[inputf][sel], decision)
        if len(chains) == 1 and internal_sym == "Yes":
            selected_write(out1, inputf + "_" + str(chains[0]), symmlist[inputf][sel], "Accepted")


        ################################################# Many Chains ######################
    if len(chains) > 1:
        cesymm_int_symm = 0
        discard_small = 0
        discard_nonmemb = 0
        for ch in chains:
            inputf = pdb + "_" + ch
            print(inputf)
            if os.path.isfile(cesymm_out_dir + inputf + ".axes"):
                num = strip_tm_chains_in_order(wkdir, inputf, oriented_opm, ch)
                oriented = wkdir + inputf + "_tmp.pdb"

                if os.path.isfile(cesymm_minlen_out_dir + inputf + "_df_mlev.axes"):
                    symmlist[inputf] = {}
                    symmlist[inputf]['cesymm'] = {}
                    cesymm_dic = cesymm_data(cesymm_minlen_out_dir, inputf + "_df_mlev")
                    cesymm_dic['source'] = 'mlev'
                    cesymm_dic['files_key'] = inputf + "_df_mlev"
                    cesymm_dic['files_dir'] = locations['FSYS']['cesymm_minlen'] + pdb + "/"
                    cesymm_dic['files_dir_path'] = locations['FSYSPATH']['cesymm_minlen'] + pdb + "/"
                    symd_dic = symd_data(symd_out_dir, inputf)
                    symmlist[inputf]['symd'] = {}
                    symmlist[inputf]['symd'] = symd_dic
                elif os.path.isfile(locations['FSYSPATH']['images_cesymm_symd_log'] + pdb + "/" + pdb + "_dic.pkl"):
                    cesymm_dic = symmlist[inputf]['cesymm']
                    cesymm_dic['source'] = 'cesymm'
                    cesymm_dic['files_key'] = inputf
                    cesymm_dic['files_dir'] = locations['FSYS']['cesymm'] + pdb + "/"
                    cesymm_dic['files_dir_path'] = locations['FSYSPATH']['cesymm'] + pdb + "/"

                else:
                    symmlist[inputf] = {}
                    symmlist[inputf]['cesymm'] = {}
                    cesymm_dic = cesymm_data(cesymm_out_dir, inputf)
                    cesymm_dic['source'] = 'cesymm'
                    cesymm_dic['files_key'] = inputf
                    cesymm_dic['files_dir'] = locations['FSYS']['cesymm'] + pdb + "/"
                    cesymm_dic['files_dir_path'] = locations['FSYSPATH']['cesymm'] + pdb + "/"
                    symd_dic = symd_data(symd_out_dir, inputf)
                    symmlist[inputf]['symd'] = {}
                    symmlist[inputf]['symd'] = symd_dic

                symmlist[inputf]['cesymm'] = cesymm_dic
                symmlist[inputf]['cesymm']['chains'] = ch + ';'
                symmlist[inputf]['cesymm']['oriented'] = oriented
                all_options_write(out2, inputf, symmlist[inputf]['cesymm'])

                cesymm_dic = {}
                symmlist[inputf]['cesymm_lt'] = {}
                if os.path.isfile(cesymm_minlen_out_dir + inputf + "_lt_mlev.axes"):
                    cesymm_dic = cesymm_data(cesymm_minlen_out_dir, inputf + "_lt_mlev")
                    cesymm_dic['source'] = 'mlev_lt'
                    cesymm_dic['files_key'] = inputf + "_lt_mlev"
                    cesymm_dic['files_dir'] = locations['FSYS']['cesymm_minlen'] + pdb + "/"
                    cesymm_dic['files_dir_path'] = locations['FSYSPATH']['cesymm_minlen'] + pdb + "/"
                else:
                    cesymm_dic = cesymm_data(cesymm_lt_out_dir, inputf)
                    cesymm_dic['source'] = 'cesymm_lt'
                    cesymm_dic['files_key'] = inputf
                    cesymm_dic['files_dir'] = locations['FSYS']['cesymm_low_thr'] + pdb + "/"
                    cesymm_dic['files_dir_path'] = locations['FSYSPATH']['cesymm_low_thr'] + pdb + "/"
                symmlist[inputf]['cesymm_lt'] = cesymm_dic
                symmlist[inputf]['cesymm_lt']['chains'] = ch + ';'
                symmlist[inputf]['cesymm_lt']['oriented'] = oriented



                symmlist[inputf]['selected'] = {}
                symmlist[inputf]['selected']['source'] = 'cesymm'

                ### Is there information in the default run?
                if symmlist[inputf]['cesymm']['repeats_number'] == 'na':
                    pass
                elif int(symmlist[inputf]['cesymm']['repeats_number']) > 1:
                    symmlist[inputf]['selected']['source'] = 'cesymm'
                    reference = parse_structure(oriented)[0]
                    repeat_selection = get_repeat_resid(symmlist[inputf]['cesymm']['files_dir_path'],
                                                        symmlist[inputf]['cesymm']['files_key'], reference)
                    doms, crossings = symmetry_tm_domains(tm_archive, repeat_selection, pdb, ch)
                    symm_location = is_symmetry_in_membrane(limit_bottom, limit_top, repeat_selection, oriented)
                    symmlist[inputf]['cesymm']['rep_has_tm_domains'] = doms
                    symmlist[inputf]['cesymm']['rep_tm_domain_num'] = crossings
                    symmlist[inputf]['cesymm']['symm_within_membrane'] = symm_location
                    if doms == 'in':
                        cesymm_int_symm = 1
                    elif symm_location == 'out':
                        discard_nonmemb = 1
                    else:
                        discard_small = 1

                ### Is there information in the lower threshold run?
                if symmlist[inputf]['cesymm_lt']['repeats_number'] == 'na' or int(
                        symmlist[inputf]['cesymm_lt']['repeats_number']) == 1:
                    all_options_write(out2, inputf, symmlist[inputf]['cesymm_lt'])
                elif int(symmlist[inputf]['cesymm_lt']['repeats_number']) > 1:
                    # Compare with the current best
                    sel = symmlist[inputf]['selected']['source']

                    #   if no symmetry was found before but now we have one, select it
                    if (symmlist[inputf][sel]['repeats_number'] == 'na'
                        or int(symmlist[inputf][sel]['repeats_number']) == 1) \
                            and symmlist[inputf]['cesymm_lt']['repeats_number'] != 'na' \
                            and int(symmlist[inputf]['cesymm_lt']['repeats_number']) > 1 \
                            and int(symmlist[inputf]['cesymm_lt']['repeat_length']) > doubletmhelices:
                        symmlist[inputf]['selected']['source'] = 'cesymm_lt'

                    #  otherwise, select the current symmetry if either it both has more symmetry levels and a repeat length >doubletmhelices or both has the same number of levels and greater coverage
                    elif symmlist[inputf][sel]['repeats_number'] != 'na' \
                            and int(symmlist[inputf][sel]['repeats_number']) > 1 \
                            and symmlist[inputf]['cesymm_lt']['repeats_number'] != 'na' \
                            and int(symmlist[inputf]['cesymm_lt']['repeats_number']) > 1:
                        if (int(symmlist[inputf][sel]['symmetry_levels']) < int(
                                symmlist[inputf]['cesymm_lt']['symmetry_levels'])
                            and (int(symmlist[inputf]['cesymm_lt']['repeat_length']) >= doubletmhelices
                                 and float(symmlist[inputf][sel]['coverage']) < float(
                                            symmlist[inputf]['cesymm_lt']['coverage']))) \
                                or (int(symmlist[inputf][sel]['symmetry_levels']) == int(
                            symmlist[inputf]['cesymm_lt']['symmetry_levels'])
                                    and float(symmlist[inputf][sel]['coverage']) < float(
                                    symmlist[inputf]['cesymm_lt']['coverage'])):
                            symmlist[inputf]['selected']['source'] = 'cesymm_lt'
                    all_options_write(out2, inputf, symmlist[inputf]['cesymm_lt'])

                ### Is there information in the maxrmsd runs?
                maxrmsd_options = ['cesymm_rmsd_2_5', 'cesymm_rmsd_3']
                for rmsd_option in maxrmsd_options:
                    cesymm_dic = {}
                    cesymm_rmsd_out_dir = locations['FSYSPATH'][rmsd_option] + pdb + "/"
                    if os.path.isfile(cesymm_rmsd_out_dir + inputf + ".axes"):
                        cesymm_dic = cesymm_data(cesymm_rmsd_out_dir, inputf)
                        name_dic = rmsd_option[7:]
                        cesymm_dic['source'] = name_dic
                        cesymm_dic['files_key'] = inputf
                        cesymm_dic['files_dir'] = locations['FSYS'][rmsd_option] + pdb + "/"
                        cesymm_dic['files_dir_path'] = locations['FSYSPATH'][rmsd_option] + pdb + "/"
                        symmlist[inputf][name_dic] = cesymm_dic
                        symmlist[inputf][name_dic]['chains'] = ch + ';'
                        symmlist[inputf][name_dic]['oriented'] = oriented

                        if symmlist[inputf][name_dic]['repeats_number'] == 'na' or int(symmlist[inputf][name_dic]['repeats_number']) == 1:
                            all_options_write(out2, inputf, symmlist[inputf][name_dic])
                        elif int(symmlist[inputf][name_dic]['repeats_number']) > 1:
                            # Compare with the current best
                            sel = symmlist[inputf]['selected']['source']

                            #   if no symmetry was found before but now we have one, select it
                            if (symmlist[inputf][sel]['repeats_number'] == 'na'
                                or int(symmlist[inputf][sel]['repeats_number']) == 1) \
                                    and symmlist[inputf][name_dic]['repeats_number'] != 'na' \
                                    and int(symmlist[inputf][name_dic]['repeats_number']) > 1 \
                                    and int(symmlist[inputf][name_dic]['repeat_length']) > doubletmhelices:
                                symmlist[inputf]['selected']['source'] = name_dic

                            # otherwise, select the current symmetry if either it both has more symmetry levels and a repeat length >doubletmhelices
                            # or both has the same number of levels and greater coverage
                            elif symmlist[inputf][sel]['repeats_number'] != 'na' \
                                    and int(symmlist[inputf][sel]['repeats_number']) > 1 \
                                    and symmlist[inputf][name_dic]['repeats_number'] != 'na' \
                                    and int(symmlist[inputf][name_dic]['repeats_number']) > 1:

                                if (int(symmlist[inputf][sel]['symmetry_levels']) < int(symmlist[inputf][name_dic]['symmetry_levels'])
                                    and (int(symmlist[inputf][name_dic]['repeat_length']) >= doubletmhelices
                                         and float(symmlist[inputf][sel]['coverage']) < float(symmlist[inputf][name_dic]['coverage']))) \
                                        or (int(symmlist[inputf][sel]['symmetry_levels']) == int(symmlist[inputf][name_dic]['symmetry_levels'])
                                            and float(symmlist[inputf][sel]['coverage']) < float(symmlist[inputf][name_dic]['coverage'])):
                                    symmlist[inputf]['selected']['source'] = name_dic
                            all_options_write(out2, inputf, symmlist[inputf][name_dic])


                ### Are there results in the minimum length/ maximum symmetry levels run?

                sel = symmlist[inputf]['selected']['source']

                # Image creation
                if 'cesymm_lt' in symmlist[inputf][sel]['source']:
                    reference = parse_structure(symmlist[inputf][sel]['oriented'])[0]
                    repeat_selection = get_repeat_resid(symmlist[inputf]['cesymm_lt']['files_dir_path'],
                                                        symmlist[inputf]['cesymm_lt']['files_key'], reference)
                    symmlist[inputf][sel]['repeat_selection'] = repeat_selection
                    image_key = inputf + "_analysis"
                    if int(symmlist[inputf][sel]['aligned_length']) >= 2 * doubletmhelices and want_image == 1:
                        symmlist[inputf][sel]['image_files'], \
                        symmlist[inputf][sel]['pml_files'], \
                        symmlist[inputf][sel]['jmol_files'] = cesymm_images(wkdir,
                                                                           symmlist[inputf]['cesymm_lt']['files_dir_path'],
                                                                           symmlist[inputf]['cesymm_lt']['files_key'],
                                                                           repeat_selection,
                                                                           symmlist[inputf][sel]['oriented'], reference,
                                                                           ray, pymol, pic_dir_chain, jmol_dir_chain,
                                                                           image_key)
                        symmlist[inputf][sel]['super_pml_files'], tm_score_super = repeats_superposition(
                            symmlist[inputf]['cesymm_lt']['files_dir_path'], symmlist[inputf][sel]['oriented'],
                            symmlist[inputf]['cesymm_lt']['files_key'], super_dir_chain, image_key,
                            symmlist[inputf]['cesymm_lt']['repeat_selection'])
                        symmlist[inputf][sel]['image_key'] = image_key

                elif 'rmsd' in symmlist[inputf][sel]['source']:
                    reference = parse_structure(symmlist[inputf][sel]['oriented'])[0]
                    repeat_selection = get_repeat_resid(symmlist[inputf][sel]['files_dir_path'],
                                                        symmlist[inputf][sel]['files_key'], reference)
                    symmlist[inputf][sel]['repeat_selection'] = repeat_selection
                    image_key = inputf + "_analysis"
                    if int(symmlist[inputf][sel]['aligned_length']) >= 2 * doubletmhelices and want_image == 1:
                        symmlist[inputf][sel]['image_files'], \
                        symmlist[inputf][sel]['pml_files'], \
                        symmlist[inputf][sel]['jmol_files'] = cesymm_images(wkdir, symmlist[inputf][sel]['files_dir_path'],
                                                                           symmlist[inputf][sel]['files_key'],
                                                                           repeat_selection,
                                                                           symmlist[inputf][sel]['oriented'], reference,
                                                                           ray, pymol, pic_dir_chain, jmol_dir_chain,
                                                                           image_key)
                        symmlist[inputf][sel]['super_pml_files'], tm_score_super = repeats_superposition(symmlist[inputf][sel]['files_dir_path'],
                                                                                symmlist[inputf][sel]['oriented'],
                                                                                symmlist[inputf][sel]['files_key'],
                                                                                super_dir_chain, image_key,
                                                                                symmlist[inputf][sel]['repeat_selection'])
                        symmlist[inputf][sel]['image_key'] = image_key

                else:
                    reference = parse_structure(symmlist[inputf][sel]['oriented'])[0]
                    repeat_selection = get_repeat_resid(symmlist[inputf]['cesymm']['files_dir_path'],
                                                        symmlist[inputf]['cesymm']['files_key'], reference)
                    symmlist[inputf][sel]['repeat_selection'] = repeat_selection
                    if symmlist[inputf][sel]['aligned_length'] != "na" and int(symmlist[inputf][sel]['aligned_length']) >= 2 * doubletmhelices and symmlist[inputf][sel][
                        'source'] == 'cesymm':
                        if 'image_files' not in symmlist[inputf][sel]:
                            symmlist[inputf][sel]['image_key'] = inputf + "_cesymm"
                            lst = glob.glob(locations['FSYSPATH']['symm_chainspngs'] + inputf + "_cesymm_*.png")
                            for f in lst:
                                image_files = image_files + f[len(locations['FSYSPATH']['symm_chainspngs']):] + ";"
                                #os.symlink(f, locations['FSYSPATH']['analysis_chainspngs'] + f[len(locations['FSYSPATH']['symm_chainspngs']):])
                                force_symlink(os.path.relpath(f, locations['FSYSPATH']['analysis_chainspngs']), locations['FSYSPATH']['analysis_chainspngs'] + f[len(locations['FSYSPATH']['symm_chainspngs']):])
                                print('Creating link {} -> {}'.format(locations['FSYSPATH']['analysis_chainspngs'] + f[len(locations['FSYSPATH']['symm_chainspngs']):], os.path.relpath(f, locations['FSYSPATH']['analysis_chainspngs'])))
                            symmlist[inputf][sel]['image_files'] = image_files
                            lst = glob.glob(
                                locations['FSYSPATH']['symm_chainspngs'] + 'pymol_script_' + inputf + "_cesymm_*.pml")
                            for f in lst:
                                pml_files = pml_files + f[len(locations['FSYSPATH']['symm_chainspngs']):] + ";"
                            symmlist[inputf][sel]['pml_files'] = pml_files
                            lst = glob.glob(locations['FSYSPATH']['symm_chainsjsons'] + '3dmol_' + inputf + "_cesymm_*.json")
                            # lst=sorted(lst, key=lambda name: int(name[len(locations['FSYSPATH']['symm_wholepngs']+'pymol_script_'+inputf+'_cesymm_'):-4]))
                            for f in lst:
                                jmol_files = jmol_files + f[len(locations['FSYSPATH']['symm_chainsjsons']):] + ";"
                            symmlist[inputf][sel]['jmol_files'] = jmol_files
                        else:
                            lst = glob.glob(locations['FSYSPATH']['symm_chainspngs'] + inputf + "_cesymm_*.png")
                            for f in lst:
                                force_symlink(os.path.relpath(f, locations['FSYSPATH']['analysis_chainspngs']), locations['FSYSPATH']['analysis_chainspngs'] + f[len(locations['FSYSPATH']['symm_chainspngs']):])
                                print('Creating link {} -> {}'.format(locations['FSYSPATH']['analysis_chainspngs'] + f[len(locations['FSYSPATH']['symm_chainspngs']):], os.path.relpath(f, locations['FSYSPATH']['analysis_chainspngs'])))
                        if len(repeat_selection) > 1:
                            symmlist[inputf][sel]['super_pml_files'], tm_score_super = repeats_superposition(
                                symmlist[inputf]['cesymm']['files_dir_path'], symmlist[inputf][sel]['oriented'],
                                symmlist[inputf]['cesymm']['files_key'], super_dir_chain, inputf + "_analysis",
                                repeat_selection)
                        else:
                            symmlist[inputf][sel]['super_pml_files'], tm_score_super = '', 'na'
                    elif symmlist[inputf][sel]['aligned_length'] != "na" and int(symmlist[inputf][sel]['aligned_length']) >= 2 * doubletmhelices:  # account for mlev cases
                        image_key = inputf + "_analysis"
                        symmlist[inputf][sel]['image_key'] = image_key
                        symmlist[inputf][sel]['image_files'], \
                        symmlist[inputf][sel]['pml_files'], \
                        symmlist[inputf][sel]['jmol_files'] = cesymm_images(wkdir,
                                                                           symmlist[inputf]['cesymm']['files_dir_path'],
                                                                           symmlist[inputf]['cesymm']['files_key'],
                                                                           repeat_selection,
                                                                           symmlist[inputf][sel]['oriented'], reference,
                                                                           ray, pymol, pic_dir_chain, jmol_dir_chain,
                                                                           image_key)
                        if len(repeat_selection) > 1:
                            symmlist[inputf][sel]['super_pml_files'], tm_score_super = repeats_superposition(
                                symmlist[inputf]['cesymm']['files_dir_path'], symmlist[inputf][sel]['oriented'],
                                symmlist[inputf]['cesymm']['files_key'], super_dir_chain, image_key, repeat_selection)
                        else:
                            symmlist[inputf][sel]['super_pml_files'], tm_score_super = '', 'na'
                        # print(pic_dir_chain,image_files)

                ### Is there information in the fragments?
                symmlist, fragments, tm_score_super = fragments_check(pdb, ch, ch, sel, tm_archive, cesymm_data, locations, symmlist, oriented, doubletmhelices, ray, pymol, want_image, pic_dir_chain, jmol_dir_chain, super_dir_chain, out2, image_files, pml_files, jmol_files, super_pml_files, tm_score_super)
                #so far we only have an image and a result written out for fragments if there is only one fragment

                sel = symmlist[inputf]['selected']['source']
                print("Selection: ", sel)
                decision = "Accepted"

                    # Check if there is at least one repeat in the membrane
                if (';' not in sel) and (symmlist[inputf][sel]['repeats_number'] == 'na' or int(
                        symmlist[inputf][sel]['repeats_number']) == 1):
                    doms, crossings = 'out', []
                    symmlist[inputf][sel]['rep_has_tm_domains'] = doms
                    symmlist[inputf][sel]['rep_tm_domain_num'] = crossings
                elif (';' not in sel):
                    doms, crossings = symmetry_tm_domains(tm_archive, symmlist[inputf][sel]['repeat_selection'], pdb, ch)
                    symmlist[inputf][sel]['rep_has_tm_domains'] = doms
                    symmlist[inputf][sel]['rep_tm_domain_num'] = crossings
                else:
                    source = sel[:-1].split(';')
                    for s in source:
                        selected_write(out1, inputf, symmlist[inputf][s], decision)

                ###
                if (';' not in sel) and (symmlist[inputf][sel]['repeats_number'] == 'na' or int(
                        symmlist[inputf][sel]['repeats_number']) == 1 or int(
                    symmlist[inputf][sel]['aligned_length']) < 2 * doubletmhelices or doms == 'out'):
                    decision = "Rejected"
                    symmlist[inputf]['selected']['info'] = "No symmetry found in the membrane-bound region during analysis."
                    for redundant in glob.glob(locations['FSYSPATH']['analysis_chainspngs'] + inputf + '_analysis*.png'):
                        if os.path.islink(redundant):
                            os.unlink(redundant)
                        else:
                            os.remove(redundant)
                    for redundant in glob.glob(
                            locations['FSYSPATH']['analysis_chainspngs'] + 'pymol_script_' + inputf + '_analysis*.pml'):
                        os.remove(redundant)
                    selected_write(out1, inputf, symmlist[inputf][sel], decision)
                elif (';' not in sel):
                    selected_write(out1, inputf, symmlist[inputf][sel], decision)
                ###
            if sel != 'cesymm' and cesymm_int_symm == 0:
                out3.write("enrichment;")
            elif sel == 'cesymm' and discard_small == 1:
                out3.write("discard_small;")
            elif sel == 'cesymm' and discard_nonmemb == 1:
                out3.write("discard_nonmemb;")
            elif sel != 'cesymm' and cesymm_int_symm == 1:
                out3.write("precision;")

    for f in glob.glob(wkdir + pdb + "*.pdb"):
        os.remove(f)
    for f in glob.glob(wkdir + pdb + "*.axes"):
        os.remove(f)
    for f in glob.glob(wkdir + pdb + "*.fasta"):
        os.remove(f)
    out1.close()
    out2.close()
    out3.write("\n")
    out3.close()
    print(symmlist)
    pkl.dump(symmlist, open(locations['FSYSPATH']['mssd_log'] + inputf[0:4] + "/" + inputf[0:4] + "_dic.pkl", "wb"))
    print(inputf[0:4] + " MSSD analysis completed.")
    print("Time[s]:%s" % (time.time() - start_time))
    print("Date: %s" % time.asctime())

    sys.stdout = save_out
    sys.stderr = save_err
    log_file.close()

    return symmlist

if __name__ == "__main__":
    if len(sys.argv) < 5:
        raise SystemExit("syntax: %s <locations_path> <pdb_options:e.g 2wcd_ordered_df> <pdb_suff> <run_id>" % sys.argv[0])
    else:
        locations_path = sys.argv[1].strip()
        inputf = sys.argv[2].strip()
        pdb_suff = sys.argv[3].strip()
        run_id = sys.argv[4].strip()

    options, locations = initialize_repository(main_path=locations_path, instr_filename_path="EncoMPASS_options_relative_to__instruction_file.txt")


    if not os.path.isdir(locations['FSYSPATH']['symmtemp']  + "webarchive/"):
        os.mkdir(locations['FSYSPATH']['symmtemp']  + "webarchive/")
    if not os.path.isdir(locations['FSYSPATH']['symmtemp']  + "webarchive/mssd/"):
        os.mkdir(locations['FSYSPATH']['symmtemp']  + "webarchive/mssd/")

    wkdir = locations['FSYSPATH']['symmtemp']  + "mssd/"
    if os.path.isdir(wkdir) == False:
        os.mkdir(wkdir)
    selected_out = locations['FSYSPATH']['mssd']
    if os.path.isdir(selected_out) == False:
        os.mkdir(selected_out)

    symmlist = results_selection(locations, options, pdb_suff, inputf, run_id = run_id)
