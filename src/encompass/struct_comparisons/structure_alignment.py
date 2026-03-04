from encompass.pipeline.supporting_functions import *

def get_redundancy_references(str_data_pdbi, pdbi, only_tmchains=False):
    redref = {}  # e.g. {"A" : {"A", "B", "C"}, "D" : {"D", "E"}, "H" : {"H"}}
    #print("REDUNDANCY", pdbi, str_data_pdbi['ENCOMPASS']['structure']['redundant_chains'])
    for ch in str_data_pdbi['ENCOMPASS']['structure']['redundant_chains']:
        if only_tmchains and ch not in str_data_pdbi['ENCOMPASS']['structure']['ktmchains']:
            continue
        ref_ch = str_data_pdbi['ENCOMPASS']['structure']['redundant_chains'][ch]
        if ref_ch in str_data_pdbi['ENCOMPASS']['structure']['redundant_chains']:
            print("WARNING:", pdbi, ref_ch, "is both a reference and a duplicate. This is a fatal error. EXIT")
            exit(1)
        if ref_ch not in redref:
            redref[ref_ch] = {ref_ch}
        redref[ref_ch].add(ch)
    key = 'kchains' if not only_tmchains else 'ktmchains'
    for ch in [x for x in str_data_pdbi['ENCOMPASS']['structure'][key] if x != "-"]:
        if ch not in str_data_pdbi['ENCOMPASS']['structure']['redundant_chains']:
            if ch not in redref:
                redref[ch] = {ch}
    return redref


def write_all_scheduled_alns(exelist, pdbich_nte_list, scheduledalns_fn):
    # BASED on the fact that pdbich_nte_list has each reference followed by its duplicates
    c = 0
    with open(scheduledalns_fn, "w") as f:
        f.write("PDB1\tCHAIN1\tPDB2\tCHAIN2\tISBETA1\tISBETA2\tNTM1\tNTM2\tDUPL1\tDUPL2\tTMSET1\tTMSET2\n")
        for ientry1, entry1 in enumerate(pdbich_nte_list):
            pdbi_ch1, cl1, ntm_list1, allLS1, dupl1 = entry1 
            pdbi1, ch1 = pdbi_ch1[:4], pdbi_ch1[5]
            b1 = True if cl1 == "beta" else False
            tmset1 = ",".join(sorted([str(x) for x in set(ntm_list1)]))
            ntm1_0 = ntm_list1[0]
            if not dupl1:
                iref1 = ientry1

            # Add comparisons with the redundant chains relative to pdbi_ch1 (but not pdbi_ch1 itself)
            iplus = 0
            while len(pdbich_nte_list) > iref1+iplus and (iplus==0 or pdbich_nte_list[iref1+iplus][4]):
                if iref1+iplus == ientry1:
                    iplus += 1
                    continue
                pdbi_ch2, cl2, ntm_list2, allLS2, dupl2 = pdbich_nte_list[iref1+iplus]
                pdbi2, ch2 = pdbi_ch2[:4], pdbi_ch2[5]
                b2 = True if cl2 == "beta" else False
                tmset2 = ",".join(sorted([str(x) for x in set(ntm_list2)]))
                ntm2_0 = ntm_list2[0]
                f.write(f"{pdbi1}\t{ch1}\t{pdbi2}\t{ch2}\t{b1}\t{b2}\t{ntm1_0}\t{ntm2_0}\t{dupl1}\t{dupl2}\t{tmset1}\t{tmset2}\n")
                iplus += 1
                c += 1

            # Add all other comparisons
            for ientry2 in exelist[iref1]:
                iplus = 0
                while iplus == 0 or (len(pdbich_nte_list) > ientry2+iplus and pdbich_nte_list[ientry2+iplus][4]):
                    pdbi_ch2, cl2, ntm_list2, allLS2, dupl2 = pdbich_nte_list[ientry2+iplus]
                    pdbi2, ch2 = pdbi_ch2[:4], pdbi_ch2[5]
                    b2 = True if cl2 == "beta" else False
                    tmset2 = ",".join(sorted([str(x) for x in set(ntm_list2)]))
                    ntm2_0 = ntm_list2[0]
                    f.write(f"{pdbi1}\t{ch1}\t{pdbi2}\t{ch2}\t{b1}\t{b2}\t{ntm1_0}\t{ntm2_0}\t{dupl1}\t{dupl2}\t{tmset1}\t{tmset2}\n")
                    iplus += 1
                    c += 1
        print("RESULTING ALIGNMENTS (INCLUDED REDUNDANT CHAINS):", c)

def comparison_lists(options, locations, str_data):

    # Effective structures that will be compared
    eff_str_data = {x for x in str_data if 'eliminated' not in str_data[x]['status']}
    print("EFFECTIVE WHOLE ENTRIES CONSIDERED:", len(eff_str_data))
    sys.stdout.flush()

    # Check on TM chains: eliminate TM chains with 0 TM regions from ktmchains in str_data
    # TODO: Move this piece of code to previous library
    pdbich_nte_list = []
    for pdbi in eff_str_data:
        for ich, ch in [(i, x) for i, x in enumerate(str_data[pdbi]['ENCOMPASS']['structure']['ktmchains']) if x != '-']:
            ntm = len(str_data[pdbi]['ENCOMPASS']['structure']['chains'][ich]['TM_regions']['TM_regions_extrema'])
            if not ntm:
                print("WARNING:", pdbi, ch, "is a TM chain with 0 TM regions, and it will be removed. But this should have been caught before this point!! EXIT")
                exit(1)
                str_data[pdbi]['ENCOMPASS']['structure']['ktmchains'][ich] = '-'

    # Create main information lookup list
    # STRUCTURE: contains reference and then its duplicates
    # EXAMPLE: if duplicates in 1abc are ("A", "C", "E") and ("B", "D"), the order will be:  1abc_A, 1abc_C, 1abc_E, 1abc_B, 1abc_D
    print("SELECTION LOOP 1", time.ctime())
    sys.stdout.flush()
    nondup_set = set()
    for pdbi in sorted(list(eff_str_data)):
        cl = str_data[pdbi]['ENCOMPASS']['class']
        local_lists = {}
        
        # Get the reference chains in the pdb entry (example: {"A" : {"A", "B", "C"}, "D" : {"D", "E"}})
        redundancy_references = get_redundancy_references(str_data[pdbi], pdbi, only_tmchains=True) 
        

        """
        # For eliminating references of references (which is an error and should be corrected at the source!)
        redundancy_list_d2 = {}
        for ref_ch in redundancy_list_within_entry:
            for ch in redundancy_list_within_entry[ref_ch]:
                if ch not in redundancy_list_d2:
                    redundancy_list_d2[ch] = redundancy_list_within_entry[ref_ch]
        redundancy_list_within_entry = redundancy_list_d2
        """

        
        for ch in sorted(list(redundancy_references)): 
            pdbi_ch = pdbi + '_' + ch
            # Get list of indices corresponding to the duplicate chains of ch (ch included!)
            ichlist = sorted([str_data[pdbi]['ENCOMPASS']['structure']['kchains'].index(x) for x in redundancy_references[ch]])
            print("ICHLIST", ichlist)
            N_TM_regions = []
            LENGTH_stats = []
            for x in ichlist:
                N_TM_regions.append(len(str_data[pdbi]['ENCOMPASS']['structure']['chains'][x]['TM_regions']['TM_regions_extrema']))
                print(pdbi, x, str_data[pdbi]['ENCOMPASS']['structure']['chains'][x]['TM_regions']['Nterm_length'], str_data[pdbi]['ENCOMPASS']['structure']['chains'][x]['TM_regions']['middle_linkers_length'], str_data[pdbi]['ENCOMPASS']['structure']['chains'][x]['TM_regions']['Cterm_length'])
                NTL = int(str_data[pdbi]['ENCOMPASS']['structure']['chains'][x]['TM_regions']['Nterm_length'])
                MLL = int(str_data[pdbi]['ENCOMPASS']['structure']['chains'][x]['TM_regions']['middle_linkers_length'])
                CTL = int(str_data[pdbi]['ENCOMPASS']['structure']['chains'][x]['TM_regions']['Cterm_length'])
                LENGTH_stats.append((NTL, MLL, CTL))
                xch = str_data[pdbi]['ENCOMPASS']['structure']['kchains'][x]
                if xch != ch:
                    duplicate_of = ch
                    if duplicate_of not in local_lists:
                        local_lists[duplicate_of] = []
                    pdbi_x = pdbi + "_" + xch
                    local_lists[duplicate_of].append((pdbi_x, cl, N_TM_regions, LENGTH_stats, duplicate_of))
                elif xch == ch:
                    duplicate_of = None
                    nondup_set.add(pdbi_ch)
                    if ch not in local_lists:
                        local_lists[ch] = []
                    local_lists[ch].append((pdbi_ch, cl, N_TM_regions, LENGTH_stats, duplicate_of))
                print("DUPLICATE OF", duplicate_of)
        for c in sorted([x for x in local_lists]):
            print("LOCAL_LISTS", pdbi_ch, local_lists)
            pdbich_nte_list += sorted(local_lists[c], key= lambda x : x[0])

    # Filter main information list:
    #   check if there are all references of duplicates
    print("CHECK REFERENCE FOR EACH DUPLICATE", time.ctime())
    sys.stdout.flush()
    pdbich_nte_list_copy = []
    for pdbi_ch, cl, N_TM_regions, tripl, duplicate_of in pdbich_nte_list:
        if duplicate_of and pdbi_ch[:4] + "_" + duplicate_of not in nondup_set:
             print("NO REFERENCE FOR", pdbi_ch, ":", pdbi_ch[:4] + "_" + duplicate_of)
        else:
             pdbich_nte_list_copy.append((pdbi_ch, cl, N_TM_regions, tripl, duplicate_of))
    pdbich_nte_list = pdbich_nte_list_copy

    # Save main information list
    #   pickle version
    pickle.dump(pdbich_nte_list, open(locations['FSYSPATH']['cache'] + "pdbich_nte_list.pkl", "wb"))
    #   human readable version
    with open(locations['FSYSPATH']['cache'] + "entry_infos.tsv", "w") as f:
        for pdbi_ch1, cl1, ntm_list1, allLS1, dupl1 in pdbich_nte_list:
            pdbi1, ch1 = pdbi_ch1[:4], pdbi_ch1[5:]
            tmset1 = ",".join(sorted([str(x) for x in set(ntm_list1)]))
            bcl1 = True if cl1 == "beta" else False
            ntm1_0 = ntm_list1[0]
            f.write(f"{pdbi1}\t{ch1}\t{bcl1}\t{ntm1_0}\t{dupl1}\t{tmset1}\n")

    print("EFFECTIVE CHAINS CONSIDERED:", len(pdbich_nte_list), "ALL_VS_ALL WOULD INCLUDE THIS MANY COMPARISONS:", len(pdbich_nte_list)**2-len(pdbich_nte_list))
    pdbich_nte_list_nodupl = [x for x in pdbich_nte_list if not x[4]]
    print("EFFECTIVE CHAINS CONSIDERED (NO DUPL):", len(pdbich_nte_list_nodupl), "ALL_VS_ALL WOULD INCLUDE THIS MANY ACTUAL COMPARISONS:", len(pdbich_nte_list_nodupl)**2-len(pdbich_nte_list_nodupl)) 

    # Pair selection
    print("SELECTION LOOP 2", time.ctime())
    sys.stdout.flush()

    local_stats = {
        'bb' : 0,
        'minTM/maxTM >= 3/4' : 0,
        '3TM' : 0,
        '2TM-ext' : 0,
        '1TM-ext' : 0
    }

    exelist = []
    for ientry1, entry1 in enumerate(pdbich_nte_list):
        exelist.append(set())

    for ientry1, entry1 in enumerate(pdbich_nte_list):
        pdbi_ch1, cl1, ntm_list1, allLS1, dupl1 = entry1
        if dupl1:
            continue
        pdbi1, ch1 = pdbi_ch1[:4], pdbi_ch1[5:]
        tmset1 = ",".join(sorted([str(x) for x in set(ntm_list1)]))
        for ientry2, entry2 in enumerate(pdbich_nte_list):
            pdbi_ch2, cl2, ntm_list2, allLS2, dupl2 = entry2
            if dupl2:
                continue
            pdbi2, ch2 = pdbi_ch2[:4], pdbi_ch2[5:]
            tmset2 = ",".join(sorted([str(x) for x in set(ntm_list2)]))

            # Filter: same pdbi and chain
            if pdbi_ch1 == pdbi_ch2:
                #print("DISCARD", pdbi_ch1, pdbi_ch2, "SAME")
                continue
            # Filter: A-B comparison
            if cl1 != cl2:
                #print("DISCARD", pdbi_ch1, pdbi_ch2, "DIFFERENT CLASS", cl1, cl2)
                continue

            # Selection criteria
            # B-B comparisons: all
            if cl1 == 'beta':
                #print("OK", pdbi_ch1, pdbi_ch2, "BETA")
                local_stats['bb'] += 1
                exelist[ientry1].add(ientry2)
                exelist[ientry2].add(ientry1)
            # Ai-Aj (max(i,j)>3 and min(i,j)>3) comparisons: min/max >= 3/4 [i.e. 3/4, 4/5, 4/6,... but not 3/5, 3/6, 4/6, ...]
            elif min(max(ntm_list1), max(ntm_list2)) > 2 and max(max(ntm_list1), max(ntm_list2)) > 3:
                if max(min(ntm_list1), min(ntm_list2))*0.75 <= min(max(ntm_list1), max(ntm_list2)):
                    #print("OK", pdbi_ch1, pdbi_ch2, "Ai-Aj", ntm_list1[0], ntm_list2[0], max(min(ntm_list1), min(ntm_list2))*0.75 <= min(max(ntm_list1), max(ntm_list2)))
                    local_stats['minTM/maxTM >= 3/4'] += 1
                    exelist[ientry1].add(ientry2)
                    exelist[ientry2].add(ientry1)
                #else:
                    #print("DISCARD", pdbi_ch1, pdbi_ch2, "Ai-Aj NO OPTION", ntm_list1[0], ntm_list2[0], max(min(ntm_list1), min(ntm_list2))*0.75 <= min(max(ntm_list1), max(ntm_list2)))
            # A3-A3 comparisons: all
            elif 3 in ntm_list1 and 3 in ntm_list2:
                local_stats['3TM'] += 1
                exelist[ientry1].add(ientry2)
                exelist[ientry2].add(ientry1)
            # A1-A1 comparisons: no big extrema or |max(1a,1b)-max(2a,2b)| < min(max(1a,1b), max(2a,2b))
            # (the difference between the biggest in each chain must not exceed the smallest of the two)
            elif 1 in ntm_list1 and 1 in ntm_list2:
                found = False
                for ls1 in allLS1:
                    for ls2 in allLS2:
                        ntl1, _, ctl1 = ls1
                        ntl2, _, ctl2 = ls2
                        if (max(ntl1, ctl1) < 100 and max(ntl2, ctl2) < 100) or (abs(max(ntl1, ctl1) - max(ntl2, ctl2)) < min(max(ntl1, ctl1), max(ntl2,  ctl2))):
                            #print("OK", pdbi_ch1, pdbi_ch2, "A1-A1", ls1, ls2, "({0} AND {1}) OR {2}".format(max(ntl1, ctl1) < 100, max(ntl2, ctl2) < 100, abs(max(ntl1, ctl1) - max(ntl2, ctl2)) < min(max(ntl1, ctl1), max(ntl2,  ctl2))))
                            local_stats['1TM-ext'] += 1
                            exelist[ientry1].add(ientry2)
                            exelist[ientry2].add(ientry1)
                            found = True
                            break
                    if found:
                        break
                #if not found:
                    #print("DISCARD", pdbi_ch1, pdbi_ch2, "A1-A1 NO OPTION", allLS1, allLS2, "max(ntl1, ctl1) < 100 and max(ntl2, ctl2) < 100) or (abs(max(ntl1, ctl1) - max(ntl2, ctl2)) < min(max(ntl1, ctl1), max(ntl2,  ctl2))")
            # A2-A2 comparisons : no big extrema or (A1-A1 condition and (no big intra or similar intras))
            elif 2 in ntm_list1 and 2 in ntm_list2:
                found = False
                for ls1 in allLS1:
                    for ls2 in allLS2:
                        ntl1, mll1, ctl1 = ls1
                        ntl2, mll2, ctl2 = ls2
                        if (
                                ((max(ntl1, ctl1) < 100 and max(ntl2, ctl2) < 100) 
                                or (abs(max(ntl1, ctl1) - max(ntl2, ctl2)) < min(max(ntl1, ctl1), max(ntl2, ctl2))))
                                and ((mll1 < 100 and mll2 < 100) 
                                or (abs(mll1 - mll2) < min(mll1, mll2)))
                                ):
                            #print("OK", pdbi_ch1, pdbi_ch2, "A2-A2", ls1, ls2, "({0} OR {1}) AND ({2} OR {3})".format(max(ntl1, ctl1) < 100 and max(ntl2, ctl2) < 100, abs(max(ntl1, ctl1) - max(ntl2, ctl2)) < min(max(ntl1, ctl1), max(ntl2, ctl2)), mll1 < 100 and mll2 < 100, abs(mll1 - mll2) < min(mll1, mll2)))
                            local_stats['2TM-ext'] += 1
                            exelist[ientry1].add(ientry2)
                            exelist[ientry2].add(ientry1)
                            found = True
                            break
                    if found:
                        break
                #if not found:
                    #print("DISCARD", pdbi_ch1, pdbi_ch2, "A2-A2 NO OPTION", allLS1, allLS2, "((max(ntl1, ctl1) < 100 and max(ntl2, ctl2) < 100) or (abs(max(ntl1, ctl1) - max(ntl2, ctl2)) < min(max(ntl1, ctl1), max(ntl2, ctl2)))) and ((mll1 < 100 and mll2 < 100) or (abs(mll1 - mll2) < min(mll1, mll2)))")
            #else:
                #print("DISCARD", pdbi_ch1, pdbi_ch2, "NO OPTION", ntm_list1, ntm_list2)

    for i in range(len(exelist)):
        exelist[i] = tuple(exelist[i])

    print("EXELIST", exelist)
    print("SCHEDULED ALIGNMENT TO LAUNCH:", functools.reduce(lambda count, l: count + len(l), exelist, 0)) 
    print(time.ctime())
    pickle.dump(exelist, open(locations['FSYSPATH']['cache'] + "exelist.pkl", "wb")) 

    """
    # Remove cases where the reference chain is not in the data structure anymore because, for example, it had structural errors
    print("ORPHANS?")
    sys.stdout.flush()
    dellines = []
    for i, t in enumerate(pre_exelist):
        pdbi1, ch1, pdbi2, ch2, isb1, isb2, ntm_list1, ntm_list2, dupl1, dupl2, _, _ = t
        if (dupl1 and pdbi1 + "_" + dupl1 not in nondup_set[pdbi2 + "_" + ch2]) or (dupl2 and pdbi2 + "_" + dupl2 not in nondup_set[pdbi1 + "_" + ch1]):
            dellines.append(i)
            print("DELETE ORPHAN", t)
            sys.stdout.flush()
     """

    write_all_scheduled_alns(exelist, pdbich_nte_list, locations['SYSFILES']['H_scheduledalns'])

    # Show statistics
    for x in sorted(local_stats):
        print("STATS", "criterion", x, local_stats[x], local_stats[x])
    print("STATS", "criterion", "Finished", "time", time.ctime(), "\n")
    sys.stdout.flush()
    print(time.ctime())

    gc.collect()

    return exelist, pdbich_nte_list


def structure_alignment(options, locations, str_data, entry_exelist, pdbich_nte_list, only_table=False, only_gather=False, split_exelist=False): # more_list=[], less_list=[], entry_exelist=None, only_table=False, only_gather=False, split_exelist=False):
    this_name = 'structure_alignment'
#    external_filename = options['output_tab']

    #less_list = list(set(less_list) - set(more_list))

    # BOH
    already_processed = []
    ex_list = {}

    if split_exelist:
        exelists = []
#        tmranges = [[1], [2], [3], [4,5], [6,7,8,9], [x for x in range(10, 100)]]
        tmranges = [[1,2,3],[4,5,6],[x for x in range(7, 100)]]
        for tmrange in tmranges:
            el = []
            for ientry, entry in enumerate(pdbich_nte_list):
                pdbi_ch, cl, ntm_list, allLS, dupl = entry
                if max(ntm_list) in tmrange:
                    el.append(entry_exelist[ientry])
            exelists.append(el)
    else:
        exelists = [entry_exelist]

    ip1 = -1
    is_first = True
    for ic, exelist in enumerate(exelists):
        batch_job_code = options['RUN'][('code', 'main_code')] + "FrTMAlign" + str(ic).zfill(3)
        output_dir = locations['FSYSPATH']['cache'] + 'output_{0}/'.format(batch_job_code)
        outputs = ['stats_<id>.txt', 'seq_seqalns_<id>.txt', 'str_seqalns_<id>.txt', 'stralns_<id>.txt.gz'] 
 
        if not only_table:
            # Instruct locusts
            if options['PATHS'][('sing', 'singularity')] and options['PATHS'][('container', 'singularity_container')]:
                frtm_path = options['PATHS'][('sigfrtmalign', 'sig_frtmalign_path')]
                muscle_path = options['PATHS'][('sigmuscle', 'sig_muscle_path')]
            else:
                if options['RUN'][('hpc', 'run_on_hpc')]:
                    frtm_path = options['PATHS'][('hpcfrtmalign', 'hpc_frtmalign_path')]
                    muscle_path = options['PATHS'][('hpcmuscle', 'hpc_muscle_path')]
                else:
                    frtm_path = options['PATHS'][('frtmalign', 'frtmalign_path')]
                    muscle_path = options['PATHS'][('muscle', 'muscle_path')]

            # 0. Create locusts parameter file
            parameter_file = locations['FSYSPATH']['logs'] + 'FrTMAlign_locusts.par'
            write_locusts_parfile(options, parameter_file, options['RUN'][('code', 'main_code')] + '_FrTMAlign', only_gather=only_gather)

            # .1 Create input and output dir
            # Input (do not rewrite if only_gather)
            input_dir = locations['FSYSPATH']['cache'] + 'input_{0}/'.format(batch_job_code)
            if not only_gather:
                if os.path.exists(input_dir):
                    shutil.rmtree(input_dir)
                os.mkdir(input_dir)
    
            # Output (erase in any case)
            if not os.path.exists(output_dir):
                os.mkdir(output_dir)
    
            # Script (do not rewrite if only_gather)
            seqid_py = input_dir + "seqid.py"
            if not only_gather:
               shutil.copyfile(locations['SYSFILES']['seqid_py'], seqid_py)
        
            # .2 Create the local script and list of files        
            pdbiset = set()
            entrylist_fn = output_dir + "/entry_list.txt"
            entrylist_f = open(entrylist_fn, 'w')
            for ilist1, list1 in enumerate(exelist):
                if not list1:
                    continue
                ip1 += 1
                pdbi_ch1, cl1, ntm_list1, allLS1, dupl1 = pdbich_nte_list[ilist1]
                pdbiset.add(pdbi_ch1)
                entrylist_f.write("{0:10d}\t{1:6s}\n".format(ip1, pdbi_ch1))
                # If only_gather, only fill pdbiset
                if only_gather:
                    continue
                locscr_filename = input_dir + "straln_exe_{0}.sh".format(ip1)
                with open(locscr_filename, 'w') as lsf:
                    lsf.write('PRESHARED="${1}"\n')
                    lsf.write('SHARED=`dirname ${PRESHARED}`\n')
                    lsf.write('pdb1={0}\n'.format(pdbi_ch1))
                    if options['ALL']['nodescratchpath']:
                        lsf.write('THIS=$(pwd); mkdir {0}; cp strlist_{1}.txt {0}; cd {0}\n'.format(options['ALL']['nodescratch']+'/'+str(ip1), ip1))
                        lsf.write('cp ${SHARED}/*/seqid.py .\n')
                        lsf.write('for pdb2 in `cat strlist_{0}.txt`\n'.format(ip1))
                        lsf.write('do\n')
                        lsf.write('    cp ${SHARED}/*/${pdb1}_enc.pdb ${SHARED}/*/${pdb2}_enc.pdb .\n\n')
                        lsf.write('done\n')
                    lsf.write('for pdb2 in `cat strlist_{0}.txt`\n'.format(ip1))
                    lsf.write('do\n')
                    if not options['ALL']['nodescratchpath']:
                        lsf.write('    cp ${SHARED}/*/${pdb1}_enc.pdb ${SHARED}/*/${pdb2}_enc.pdb .\n')
                    lsf.write('    {0} ${{pdb1}}_enc.pdb ${{pdb2}}_enc.pdb -o prestraln_{1}.txt > out_tmp_{1}.txt\n'.format(frtm_path, ip1))
                    lsf.write("    sed -i 's/\\x0/X/g' out_tmp_{0}.txt\n".format(ip1))
                    lsf.write('    echo "INIT ${{pdb1}} ${{pdb2}}" >> stralns_{0}.txt\n'.format(ip1))
                    lsf.write('    grep "ATOM\|TER" prestraln_{0}.txt >> stralns_{0}.txt\n'.format(ip1))
                    lsf.write('    echo "INIT ${{pdb1}} ${{pdb2}}" >> str_seqalns_{0}.txt\n'.format(ip1))
                    lsf.write('    echo "INIT ${{pdb1}} ${{pdb2}}" >> seq_seqalns_{0}.txt\n'.format(ip1))
                    lsf.write('    grep -A3 "(\\\":" out_tmp_{0}.txt | awk -v pdb1=${{pdb1}} -v pdb2=${{pdb2}} \'NR==2{{print ">"pdb1; print}} NR==4{{print ">"pdb2; print}}\' > str_seq_prealn_{0}.fa\n'.format(ip1))
                    lsf.write('    {0} -in str_seq_prealn_{1}.fa -out seq_seq_prealn_{1}.fa\n'.format(muscle_path, ip1))
                    lsf.write('    SEQSEQID=`python3 seqid.py seq_seq_prealn_{0}.fa | awk \'{{print $3}}\'`\n'.format(ip1))
                    lsf.write('    cat seq_seq_prealn_{0}.fa >> seq_seqalns_{0}.txt\n'.format(ip1))
                    lsf.write('    STRSEQID=`python3 seqid.py str_seq_prealn_{0}.fa | awk \'{{print $3}}\'`\n'.format(ip1))
                    lsf.write('    cat str_seq_prealn_{0}.fa >> str_seqalns_{0}.txt\n'.format(ip1))
                    lsf.write('    echo -n "${{pdb1}} ${{pdb2}} " >> stats_{0}.txt\n'.format(ip1))
                    lsf.write('    grep "Aligned length" out_tmp_{0}.txt | awk \'BEGIN{{FS=","}}{{print $1, "=", $2, "=", $3}}\' | awk -v sesid=${{SEQSEQID}} -v stsid=${{STRSEQID}} \'BEGIN{{FS="="}}{{printf "%6d\t%8.4f\t%8.4f\t%8.4f\t%8.4f", $2, sesid, stsid, $4, $6}}\' >> stats_{0}.txt\n'.format(ip1))
                    lsf.write('    awk \'BEGIN{{print ""}}\' >> stats_{0}.txt\n'.format(ip1))
                    if not options['ALL']['nodescratchpath']:
                        lsf.write('    rm ${pdb1}_enc.pdb ${pdb2}_enc.pdb\n')
                    lsf.write('done\n')
                    lsf.write('rm *.pdb\n')
                    #lsf.write('grep "INIT\|ATOM\|TER" stralns_{0}.txt > stralns_{0}.txt\n'.format(ip1))
                    lsf.write('gzip stralns_{0}.txt\n'.format(ip1))
                    lsf.write('rm stralns_{0}.txt\n'.format(ip1))
                    if options['ALL']['nodescratchpath']:
                        lsf.write('cp stralns_{0}.txt.gz seq_seqalns_{0}.txt str_seqalns_{0}.txt stats_{0}.txt $THIS\n'.format(ip1))
                strlist_filename = input_dir + "strlist_{0}.txt".format(ip1)
                with open(strlist_filename, 'w') as strlsf:
                    for ientry2 in exelist[ilist1]:
                        pdbi_ch2, cl2, ntm_list2, allLS2, dupl2 = pdbich_nte_list[ientry2]
                        strlsf.write("{0}\n".format(pdbi_ch2))
                        pdbiset.add(pdbi_ch2)
   
            entrylist_f.close()
 
            # Copy pdb structures in input
            if not only_gather:
                for pdbi_ch in pdbiset:
                    shutil.copyfile(locations['FSYSPATH']['chains'] + '{0}_enc.pdb'.format(pdbi_ch), input_dir + '{0}_enc.pdb'.format(pdbi_ch))
        
            # .3 Instruct locusts
            specific_inputs = ['strlist_<id>.txt', 'straln_exe_<id>.sh']
            command_template = 'cp <shared>pyscr . ; S=`dirname <shared>{0} `; bash straln_exe_<id>.sh ${{S}}'.format([x for x in pdbiset][0])
            shared_inputs = ['{0}:{0}_enc.pdb'.format(x) for x in pdbiset] + ['pyscr:seqid.py']
        
            # .4 Launch
            locusts.swarm.launch(
                indir=input_dir,
                outdir=output_dir,
                code=batch_job_code,
                spcins=specific_inputs,
                shdins=shared_inputs,
                outs=outputs,
                cmd=command_template,
                parf=parameter_file
            )

        gc.collect()
        is_first = False

    lengths = {}
    for pdbi in str_data:
        for i, x in enumerate(str_data[pdbi]['ENCOMPASS']['structure']['ktmchains']):
            lengths[pdbi+"_"+x] = len(str_data[pdbi]['ENCOMPASS']['structure']['chains'][i]['residues'])
    is_first = True
    outputs = ['stats_<id>.txt', 'seqseqalns_<id>.txt', 'str_seqalns_<id>.txt', 'stralns_<id>.txt.gz']
    for ic, exelist in enumerate(exelists):
        batch_job_code = options['RUN'][('code', 'main_code')] + "FrTMAlign" + str(ic).zfill(3)
        output_dir = locations['FSYSPATH']['cache'] + "output_" + batch_job_code + "/"
        move_and_filter(locations, output_dir, str_data)
        table = make_new_table(locations, exelist, pdbich_nte_list, str_data, lengths, output_dir, outputs, append=(not is_first))
        is_first = False
    return table


def move_and_filter(locations, output_dir, str_data):

    # id -> pdbi
    entrylist_fn = output_dir + "/entry_list.txt"
    entryd = {}
    with open(entrylist_fn) as f:
        for line in f:
            i, pdbi = [x.strip() for x in line.split("\t")]
            entryd[i] = pdbi

    # subpaths list from output.log
    subpaths = []
    with open(output_dir + 'output.log') as outf:
        for line in outf:
            if (not line.startswith("stats_")) and line.split()[1] == "present":
                subpaths.append(line.split()[2])

    # copy
    for subpath in subpaths:
        bn = os.path.basename(subpath)
        the_id = bn.split("_")[-1].split(".")[0]
        pdbi_ch = entryd[the_id]
        pdbi = pdbi_ch[:4]
        ch = pdbi_ch[5]

        # Recover all redundant chains
        pdbilist = [pdbi_ch]
        redundancy_references = get_redundancy_references(str_data[pdbi], pdbi, only_tmchains=True)
        if ch in redundancy_references:
            for rch in redundancy_references[ch]:
                pdbilist.append(pdbi + "_" + rch)

        # Copy once for each redundant chain
        for pdbi_ch in pdbilist:
            new_bn = bn.replace("_"+the_id+".", "_"+pdbi_ch+".")
            if "seq_seqaln" in new_bn:
                shutil.copyfile(subpath, locations['FSYSPATH']['seqseqalns'] + "/" + new_bn)
            if "str_seqaln" in new_bn:
                shutil.copyfile(subpath, locations['FSYSPATH']['strseqalns'] + "/" + new_bn)
            if "straln" in new_bn:
                shutil.copyfile(subpath, locations['FSYSPATH']['stralns'] + "/" + new_bn)
    

def make_new_table(locations, exelist, pdbich_nte_list, str_data, lengths, output_dir, outputs, equivalence={}, allow_incomplete=False, append=False):
    this_name = make_new_table.__name__

    print("OUTPUT", output_dir)
    summary_table_fname = locations['SYSFILES']['summarytable']

    subpaths = []
    with open(output_dir + 'output.log') as outf:
        for line in outf:
            if line.startswith("stats_") and line.split()[1] == "present":
                subpaths.append(line.split()[2])

    all_redundancy_references = {}
    for pdbi in str_data:
        all_redundancy_references[pdbi] = get_redundancy_references(str_data[pdbi], pdbi, only_tmchains=True)
        print("REDUNDANCY", pdbi, all_redundancy_references[pdbi])

    idx_pdbich_nte_list = {}
    for ientry, entry in enumerate(pdbich_nte_list):
        pdbi_ch, cl, ntm_list, allLS, dupl = entry
        idx_pdbich_nte_list[pdbi_ch] = ientry

    print("Create table")
    print("fill table with straln data")
    not_completed = set()
    not_completed_unknown = 0
    pairs_with_alignments = set()
    mode = "a" if append else "w"
    fout = open(summary_table_fname, mode)
    #fout.write(f"PDB1\tCHAIN1\tPDB2\tCHAIN2\tALN_LEN\tSEQ_SEQID\tSTR_SEQID\tRMSD\tTMSCORE\n")
    for subpath in subpaths[:]:
        print("SUBPATH", subpath)
        first_time_as_primary_pdb = True
        flag = True
        with open(subpath) as f:
            for line in f:
                fields = line.split()
                if len(fields) == 7:
                    pdbc1, pdbc2, alnlen, seqseqid, strseqid, rmsd, tmscore = fields
                    pdbi1, ch1 = pdbc1.split("_")
                    pdbi2, ch2 = pdbc2.split("_")
                    if flag:
                        print(f"CHECK: {pdbc1}-{pdbc2}")
                    pairs_with_alignments.add((pdbc1, pdbc2))
                    refch1 = ch1 if ch1 not in str_data[pdbi1]['ENCOMPASS']['structure']['redundant_chains'] else str_data[pdbi1]['ENCOMPASS']['structure']['redundant_chains'][ch1]
                    refch2 = ch2 if ch2 not in str_data[pdbi2]['ENCOMPASS']['structure']['redundant_chains'] else str_data[pdbi2]['ENCOMPASS']['structure']['redundant_chains'][ch2]
                    refpdbc1 = pdbi1 + "_" + refch1
                    refpdbc2 = pdbi2 + "_" + refch2
                    if idx_pdbich_nte_list[refpdbc1] not in exelist[idx_pdbich_nte_list[refpdbc2]]:
                        print((pdbc1, pdbc2), "will not be included in summary_table.txt")
                        continue
                else:
                    if len(fields) == 2:
                        pdbc1, pdbc2 = fields
                        pdbi1, ch1 = pdbc1.split("_")
                        pdbi2, ch2 = pdbc2.split("_")
                        pairs_with_alignments.add((pdbc1, pdbc2))
                        not_completed.add((pdbc1, pdbc2))
                    else:
                        not_completed_unknown += 1
                    continue

                if first_time_as_primary_pdb:
                    print("SELF")
                    if ch1 in all_redundancy_references[pdbi1]:
                        for rch1 in all_redundancy_references[pdbi1][ch1]:
                            if pdbc1 not in lengths:
                                print(f"WARNING: NO LENGTH FOR {pdbc1}")
                                lengths[pdbc1] = 0
                            for rch2 in all_redundancy_references[pdbi1][ch1]:
                                if rch1 != rch2:
                                    print(pdbi1, rch1, pdbi1, rch2)
                                    fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(pdbi1, rch1, pdbi1, rch2, lengths[pdbc1], "1.0000", "1.0000", "0.0000", "1.0000"))

                            """
                            print(pdbi1, ch1, pdbi1, rch1)
                            print(pdbi1, rch1, pdbi1, ch1)
                            fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(pdbi1, ch1, pdbi1, rch1, lengths[pdbc1], "1.0000", "1.0000", "0.0000", "1.0000"))
                            if flag:
                                print(f"CHECK: first_time_as_primary with {pdbc2}, wrote {pdbi1} {ch1} {pdbi1} {rch1}")
                            pairs_with_alignments.add((pdbc1, pdbi1 + "_" + rch1))
                            fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(pdbi1, rch1, pdbi1, ch1, lengths[pdbc1], "1.0000", "1.0000", "0.0000", "1.0000"))
                            if flag:
                                print(f"CHECK: first_time_as_primary with {pdbc2}, wrote {pdbi1} {rch1} {pdbi1} {ch1}")
                            pairs_with_alignments.add((pdbi1 + "_" + rch1, pdbc1))
                            for rch2 in all_redundancy_references[pdbi1][ch1]:
                                if rch1 != rch2:
                                    print(pdbi1, rch1, pdbi1, rch2)
                                    fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(pdbi1, rch1, pdbi1, rch2, lengths[pdbc1], "1.0000", "1.0000", "0.0000", "1.0000"))
                            """
                    first_time_as_primary_pdb = False

                """
                print("NON-SELF (COMPUTED)")
                print(pdbi1, ch1, pdbi2, ch2)
                fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(pdbi1, ch1, pdbi2, ch2, alnlen, seqseqid, strseqid, rmsd, tmscore))
                if flag:
                    print(f"CHECK: wrote {pdbi1} {ch1} {pdbi2} {ch2}")
                pairs_with_alignments.add((pdbc1, pdbc2))
                """

                if ch1 in all_redundancy_references[pdbi1]:
                    for rch1 in all_redundancy_references[pdbi1][ch1]:
                        if ch2 in all_redundancy_references[pdbi2]:
                            for rch2 in all_redundancy_references[pdbi2][ch2]:
                                print(pdbi1, ch1, pdbi2, rch2)
                                fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(pdbi1, rch1, pdbi2, rch2, alnlen, seqseqid, strseqid, rmsd, tmscore))
                                pairs_with_alignments.add((pdbi1 + "_" + rch1, pdbi2 + "_" + rch2))
                """
                print("NON-SELF, COPIES OF 2")
                if ch2 in all_redundancy_references[pdbi2]:
                    for rch2 in all_redundancy_references[pdbi2][ch2]:
                        print(pdbi1, ch1, pdbi2, rch2)
                        fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(pdbi1, ch1, pdbi2, rch2, alnlen, seqseqid, strseqid, rmsd, tmscore))
                        pairs_with_alignments.add((pdbc1, pdbi2 + "_" + rch2))

                print("NON-SELF, COPIES OF 1")
                if ch1 in all_redundancy_references[pdbi1]:
                    for rch1 in all_redundancy_references[pdbi1][ch1]:
                        print(pdbi1, rch1, pdbi2, ch2)
                        fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(pdbi1, rch1, pdbi2, ch2, alnlen, seqseqid, strseqid, rmsd, tmscore))
                        pairs_with_alignments.add((pdbi1 + "_" + rch1, pdbc2))
                        print("NON-SELF, COPIES OF 2 OF COPIES OF 1")
                        if ch2 in all_redundancy_references[pdbi2]:
                            for rch2 in all_redundancy_references[pdbi2][ch2]:
                                print(pdbi1, rch1, pdbi2, rch2)
                                fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(pdbi1, rch1, pdbi2, rch2, alnlen, seqseqid, strseqid, rmsd, tmscore))
                                pairs_with_alignments.add((pdbi1 + "_" + rch1, pdbi2 + "_" + rch2))
                """

    tmp_fn = locations['FSYSPATH']['main'] + "tmp.txt"

    for p in not_completed:
        print("Not completed:", p)
    print("Total KNOWN not completed", len(not_completed))
    print("Total UNKNOWN not completed", not_completed_unknown)
    #print("Scheduled alignments that were omitted from the summary table (perhaps due to a difference between the scheduling table used for the run and for this analysis): ", selection.difference(pairs_with_alignments), " total: ", len(selection.difference(pairs_with_alignments)))

    return {}

def complete_straln(exelist_name, summary_name):
    exelist = pd.read_csv(exelist_name, sep="\t", header=0)#, dtype={'PDB1' : 'category', 'CHAIN1' : 'category', 'PDB2' : 'category', 'CHAIN2' : 'category', 'ISBETA1' : 'bool', 'ISBETA2' : 'bool', 'NTM1' : 'uint8', 'NTM2' : 'uint8', 'DUPLICATE_OF1' : 'category', 'DUPLICATE_OF2' : 'category', 'TM_SET_1' : 'category', 'TM_SET_2' : 'category'})
    print(exelist)
    exeset = set((exelist['PDB1']+exelist['CHAIN1']+exelist['PDB2']+exelist['CHAIN2']).unique())
    summarylist = pd.read_csv(summary_name, sep="\t", header=None, names=["PDB1", "CHAIN1", "PDB2", "CHAIN2", "ALN_LEN", "SEQ_SEQID", "STR_SEQID", "RMSD", "TMSCORE"])
    summaryset = set((summarylist['PDB1']+summarylist['CHAIN1']+summarylist['PDB2']+summarylist['CHAIN2']).unique())
    diffset = exeset - summaryset
    print("WARNING: Comparisons that are in the scheduled alignment but not the summary table (perhaps the scheduled alignment tabel changed between run and analysis): ",len(diffset), [x[0:4]+"_"+x[4:5]+"--"+x[5:9]+"_"+x[9:10] for x in diffset])

"""
def move_and_filter(in_fn, out_fn_template):
    with open(in_fn) as f:
        for line in f:
            _, pdbi_ch, _ = line.split()
            break
    out_fn = out_fn_template.replace("XXX", pdbi_ch)
    with open(in_fn) as f, open(out_fn, "w") as fo:
        copy_flag = False
        for line in f:
            if line.startswith("INIT"):
                if len(line.split()) != 3:
                    print("FORMAT ERROR", line.strip())
                    continue
                _, itself, pdbi_ch_add = line.split()
                if (itself, pdbi_ch_add) in selection:
                    #print("COPIED", (itself, pdbi_ch_add))
                    copy_flag = True
                else:
                    #print("DISCARDED", (itself, pdbi_ch_add))
                    copy_flag = False
            if copy_flag:
                fo.write(line)
    if pdbi_ch in red:
        for pdbi_ch_add in red[pdbi_ch]:
            shutil.copyfile(out_fn, out_fn_template.replace("XXX", pdbi_ch_add))
"""

if __name__ == "__main__":
    from encompass.pipeline.supporting_functions import *
    from encompass.pipeline.initialize_repository import initialize_repository
    from encompass.sources.combine_sources import *
    from encompass.sources.complete_information import *

    options, locations = initialize_repository()

    ### REPAIR STEP : UNCOMMENT THIS SECTION
    #str_data = read_checkpoint(locations['FSYSPATH']['cache'] + 'str_data_completegen.pkl', locations['SYSFILES']['data_structure_template'])
    #str_data = find_redundant_chains(locations, str_data) #REPAIR STEP
    #write_checkpoint(str_data, locations['FSYSPATH']['cache'] + 'str_data_straln.pkl')
    ###

    #test_list = ['1a0s', '1a0t']
    str_data = read_checkpoint(locations['FSYSPATH']['cache'] + 'str_data_completegen.pkl', locations['SYSFILES']['data_structure_template'])
    #print("THIS IS A TEST!!")
    #str_data = {x : str_data[x] for x in str_data if x in test_list}

    print(">>>exelist")
    sys.stdout.flush()
    exelist, pdbich_nte_list = comparison_lists(options, locations, str_data)

    print(">>>straln")
    sys.stdout.flush()
    print("LOAD EXELIST")
    exelist = pickle.load(open(locations['FSYSPATH']['cache'] + "exelist.pkl", "rb"))
    print("LOAD PDBICH")
    pdbich_nte_list = pickle.load(open(locations['FSYSPATH']['cache'] + "pdbich_nte_list.pkl", "rb"))
    print("EXECUTE STRALN")
    structure_alignment(options, locations, str_data, exelist, pdbich_nte_list, only_table=True, only_gather=False)
    #structure_alignment(options, locations, str_data, exelist, pdbich_nte_list, only_table=False, only_gather=False)# more_list=[], less_list=[], entry_exelist=exelist, only_table=False, only_gather=False, split_exelist=False)
    complete_straln(exelist_name=locations['SYSFILES']['H_scheduledalns'], summary_name=locations['SYSFILES']['summarytable'])
