from supporting_functions import *

def simpler_comparison_lists(options, locations, str_data):

    # Effective structures that will be compared
    eff_str_data = {x for x in str_data if 'eliminated' not in str_data[x]['status']}
    print("EFFECTIVE WHOLE ENTRIES CONSIDERED:", len(eff_str_data))

    pdbich_nte_list = []
    # Look at each TM chain: if it # of TM resghions is 0, remove chain from TM chains 
    # WARNING: This modifies str_data!
    for pdbi in eff_str_data:
        for ich, ch in [(i, x) for i, x in enumerate(str_data[pdbi]['ENCOMPASS']['structure']['ktmchains']) if x != '-']:
            N_TM_regions = len(str_data[pdbi]['ENCOMPASS']['structure']['chains'][ich]['TM_regions']['TM_regions_extrema'])
            if not N_TM_regions:
                str_data[pdbi]['ENCOMPASS']['structure']['ktmchains'][ich] = '-'

    # Loop on all TM chains (see above) 
    nondup_set = set()
    for pdbi in sorted(list(eff_str_data)):#[:200] + ["3cn5", "2rdd"]:  #HERE!!!!
        cl = str_data[pdbi]['ENCOMPASS']['class']
        local_lists = {}

        redundancy_list_within_entry = {}  # e.g. {'A': ['A'], 'C': ['C'], 'E': ['E']}
        for ch in [x for x in str_data[pdbi]['ENCOMPASS']['structure']['ktmchains'] if x != '-']:
            if ch in str_data[pdbi]['ENCOMPASS']['structure']['redundant_chains']:
                ref_ch = str_data[pdbi]['ENCOMPASS']['structure']['redundant_chains'][ch]
            else:
                ref_ch = ch
            if ref_ch not in redundancy_list_within_entry:
                redundancy_list_within_entry[ref_ch] = []
            redundancy_list_within_entry[ref_ch].append(ch)
        # For eliminating references of references (which is an error and should be corrected at the source!)
        redundancy_list_d2 = {}
        for ref_ch in redundancy_list_within_entry:
            for ch in redundancy_list_within_entry[ref_ch]:
                if ch not in redundancy_list_d2:
                    redundancy_list_d2[ch] = redundancy_list_within_entry[ref_ch]
        redundancy_list_within_entry = redundancy_list_d2

        for ch in redundancy_list_within_entry: #ich, ch in [(i, x) for i, x in enumerate(str_data[pdbi]['ENCOMPASS']['structure']['ktmchains']) if x != '-']:
            pdbi_ch = pdbi + '_' + ch
            ichlist = [str_data[pdbi]['ENCOMPASS']['structure']['kchains'].index(x) for x in redundancy_list_within_entry[ch]]
            N_TM_regions = []
            LENGTH_stats = []
            for x in ichlist:
                N_TM_regions.append(len(str_data[pdbi]['ENCOMPASS']['structure']['chains'][x]['TM_regions']['TM_regions_extrema']))
                NTL = int(str_data[pdbi]['ENCOMPASS']['structure']['chains'][x]['TM_regions']['Nterm_length'])
                MLL = int(str_data[pdbi]['ENCOMPASS']['structure']['chains'][x]['TM_regions']['middle_linkers_length'])
                CTL = int(str_data[pdbi]['ENCOMPASS']['structure']['chains'][x]['TM_regions']['Cterm_length'])
                LENGTH_stats.append((NTL, MLL, CTL))
            if ch in str_data[pdbi]['ENCOMPASS']['structure']['redundant_chains'] and str_data[pdbi]['ENCOMPASS']['structure']['redundant_chains'][ch] in str_data[pdbi]['ENCOMPASS']['structure']['ktmchains']:
                duplicate_of = str_data[pdbi]['ENCOMPASS']['structure']['redundant_chains'][ch]
                if duplicate_of not in local_lists:
                    local_lists[duplicate_of] = []
                local_lists[duplicate_of].append((pdbi_ch, cl, N_TM_regions, LENGTH_stats, duplicate_of))
            else:
                duplicate_of = None
                nondup_set.add(pdbi_ch)
                if ch not in local_lists:
                    local_lists[ch] = []
                local_lists[ch].append((pdbi_ch, cl, N_TM_regions, LENGTH_stats, duplicate_of))
        for c in local_lists:
            pdbich_nte_list += local_lists[c]

    pdbich_nte_list_copy = []
    for pdbi_ch, cl, N_TM_regions, tripl, duplicate_of in pdbich_nte_list:
        if duplicate_of and pdbi_ch[:4] + "_" + duplicate_of not in nondup_set:
             print("NO ORIGINAL!", pdbi_ch[:4] + "_" + duplicate_of, pdbi_ch)
        else:
             pdbich_nte_list_copy.append((pdbi_ch, cl, N_TM_regions, tripl, duplicate_of))

    pdbich_nte_list = pdbich_nte_list_copy
    pickle.dump(pdbich_nte_list, open(locations['FSYSPATH']['cache'] + "pdbich_nte_list.pkl", "wb"))

    pre_exelist = []
    ind_exelist = []
    local_stats = {
        'bb' : 0,
        'minTM/maxTM >= 3/4' : 0,
        '3TM' : 0,
        '2TM-ext' : 0,
        '1TM-ext' : 0
    }


    for pdbi_ch1, cl1, ntm1, allLS1, dupl1 in pdbich_nte_list:
        pdbi1, ch1 = pdbi_ch1[:4], pdbi_ch1[5:]
        tmset1 = ",".join(sorted([str(x) for x in set(ntm1)]))
        for pdbi_ch2, cl2, ntm2, allLS2, dupl2 in pdbich_nte_list:
            pdbi2, ch2 = pdbi_ch2[:4], pdbi_ch2[5:]
            tmset2 = ",".join(sorted([str(x) for x in set(ntm2)]))
            print(pdbi_ch1, pdbi_ch2)
            if pdbi_ch1 == pdbi_ch2:
                print("DISCARD", pdbi_ch1, pdbi_ch2, "SAME")
                continue
            # A-B comparisons: none
            if cl1 != cl2:
                print("DISCARD", pdbi_ch1, pdbi_ch2, "DIFFERENT CLASS", cl1, cl2)
                continue

            p1p2 = pdbi_ch1 + '__' + pdbi_ch2
            # B-B comparisons: all
            if cl1 == 'beta':
                print("OK", pdbi_ch1, pdbi_ch2, "BETA")
                local_stats['bb'] += 1
                pre_exelist.append((pdbi1, ch1, pdbi2, ch2, True, True, ntm1[0], ntm2[0], dupl1, dupl2, tmset1, tmset2))
                ind_exelist.append(p1p2)
            # Ai-Aj (max(i,j)>3 and min(i,j)>3) comparisons: min/max >= 3/4 [i.e. 3/4, 4/5, 4/6,... but not 3/5, 3/6, 4/6, ...]
            elif min(max(ntm1), max(ntm2)) > 2 and max(max(ntm1), max(ntm2)) > 3:
                if max(min(ntm1), min(ntm2))*0.75 <= min(max(ntm1), max(ntm2)):
                    print("OK", pdbi_ch1, pdbi_ch2, "Ai-Aj", ntm1[0], ntm2[0], max(min(ntm1), min(ntm2))*0.75 <= min(max(ntm1), max(ntm2)))
                    local_stats['minTM/maxTM >= 3/4'] += 1
                    pre_exelist.append((pdbi1, ch1, pdbi2, ch2, False, False, ntm1[0], ntm2[0], dupl1, dupl2, tmset1, tmset2))
                    ind_exelist.append(p1p2)
                else:
                    print("DISCARD", pdbi_ch1, pdbi_ch2, "Ai-Aj NO OPTION", ntm1[0], ntm2[0], max(min(ntm1), min(ntm2))*0.75 <= min(max(ntm1), max(ntm2)))
            # A3-A3 comparisons: all
            elif 3 in ntm1 and 3 in ntm2:
                local_stats['3TM'] += 1
                pre_exelist.append((pdbi1, ch1, pdbi2, ch2, False, False, ntm1[0], ntm2[0], dupl1, dupl2, tmset1, tmset2))
                ind_exelist.append(p1p2)
            # A1-A1 comparisons: no big extrema or |max(1a,1b)-max(2a,2b)| < min(max(1a,1b), max(2a,2b))
            # (the difference between the biggest in each chain must not exceed the smallest of the two)
            elif 1 in ntm1 and 1 in ntm2:
                found = False
                for ls1 in allLS1:
                    for ls2 in allLS2:
                        ntl1, _, ctl1 = ls1
                        ntl2, _, ctl2 = ls2
                        if (max(ntl1, ctl1) < 100 and max(ntl2, ctl2) < 100) or (abs(max(ntl1, ctl1) - max(ntl2, ctl2)) < min(max(ntl1, ctl1), max(ntl2,  ctl2))):
                            print("OK", pdbi_ch1, pdbi_ch2, "A1-A1", ls1, ls2, "({0} AND {1}) OR {2}".format(max(ntl1, ctl1) < 100, max(ntl2, ctl2) < 100, abs(max(ntl1, ctl1) - max(ntl2, ctl2)) < min(max(ntl1, ctl1), max(ntl2,  ctl2))))
                            local_stats['1TM-ext'] += 1
                            pre_exelist.append((pdbi1, ch1, pdbi2, ch2, False, False, ntm1[0], ntm2[0], dupl1, dupl2, tmset1, tmset2))
                            ind_exelist.append(p1p2)
                            found = True
                            break
                    if found:
                        break
                if not found:
                    print("DISCARD", pdbi_ch1, pdbi_ch2, "A1-A1 NO OPTION", allLS1, allLS2, "max(ntl1, ctl1) < 100 and max(ntl2, ctl2) < 100) or (abs(max(ntl1, ctl1) - max(ntl2, ctl2)) < min(max(ntl1, ctl1), max(ntl2,  ctl2))")
            # A2-A2 comparisons : no big extrema or (A1-A1 condition and (no big intra or similar intras))
            elif 2 in ntm1 and 2 in ntm2:
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
                            print("OK", pdbi_ch1, pdbi_ch2, "A2-A2", ls1, ls2, "({0} OR {1}) AND ({2} OR {3})".format(max(ntl1, ctl1) < 100 and max(ntl2, ctl2) < 100, abs(max(ntl1, ctl1) - max(ntl2, ctl2)) < min(max(ntl1, ctl1), max(ntl2, ctl2)), mll1 < 100 and mll2 < 100, abs(mll1 - mll2) < min(mll1, mll2)))
                            local_stats['2TM-ext'] += 1
                            pre_exelist.append((pdbi1, ch1, pdbi2, ch2, False, False, ntm1[0], ntm2[0], dupl1, dupl2, tmset1, tmset2))
                            ind_exelist.append(p1p2)
                            found = True
                            break
                    if found:
                        break
                if not found:
                    print("DISCARD", pdbi_ch1, pdbi_ch2, "A2-A2 NO OPTION", allLS1, allLS2, "((max(ntl1, ctl1) < 100 and max(ntl2, ctl2) < 100) or (abs(max(ntl1, ctl1) - max(ntl2, ctl2)) < min(max(ntl1, ctl1), max(ntl2, ctl2)))) and ((mll1 < 100 and mll2 < 100) or (abs(mll1 - mll2) < min(mll1, mll2)))")
            else:
                print("DISCARD", pdbi_ch1, pdbi_ch2, "NO OPTION", ntm1, ntm2)

    nondup_set = {}
    for pdbi1, ch1, pdbi2, ch2, isb1, isb2, ntm1, ntm2, dupl1, dupl2, _, _ in pre_exelist:
        if pdbi2 + "_" + ch2 not in nondup_set:
            nondup_set[pdbi2 + "_" + ch2] = set()
        if pdbi1 + "_" + ch1 not in nondup_set:
            nondup_set[pdbi1 + "_" + ch1] = set()
        if not dupl1:
            nondup_set[pdbi2 + "_" + ch2].add(pdbi1 + "_" + ch1)
            nondup_set[pdbi1 + "_" + ch1].add(pdbi1 + "_" + ch1)
        if not dupl2:
            nondup_set[pdbi1 + "_" + ch1].add(pdbi2 + "_" + ch2)
            nondup_set[pdbi2 + "_" + ch2].add(pdbi2 + "_" + ch2)

    # Remove cases where the reference chain is not in the data structure anymore because, for example, it had structural errors
    dellines = []
    for i, t in enumerate(pre_exelist):
        pdbi1, ch1, pdbi2, ch2, isb1, isb2, ntm1, ntm2, dupl1, dupl2, _, _ = t
        if (dupl1 and pdbi1 + "_" + dupl1 not in nondup_set[pdbi2 + "_" + ch2]) or (dupl2 and pdbi2 + "_" + dupl2 not in nondup_set[pdbi1 + "_" + ch1]):
            dellines.append(i)
            print("DELETE ORPHAN", t)

    exelist = pd.DataFrame(pre_exelist, columns=['PDB1', 'CHAIN1', 'PDB2', 'CHAIN2', 'ISBETA1', 'ISBETA2', 'NTM1', 'NTM2', 'DUPLICATE_OF1', 'DUPLICATE_OF2', 'TM_SET_1', 'TM_SET_2'])#, dtype={'PDB1' : 'category', 'CHAIN1' : 'category', 'PDB2' : 'category', 'CHAIN2' : 'category', 'ISBETA1' : 'bool', 'ISBETA2' : 'bool', 'NTM1' : 'uint8', 'NTM2' : 'uint8', 'DUPLICATE_OF1' : 'category', 'DUPLICATE_OF' : 'category'})  # Main data structure
    exelist = exelist.astype({'PDB1' : 'category', 'CHAIN1' : 'category', 'PDB2' : 'category', 'CHAIN2' : 'category', 'ISBETA1' : 'bool', 'ISBETA2' : 'bool', 'NTM1' : 'uint8', 'NTM2' : 'uint8', 'DUPLICATE_OF1' : 'category', 'DUPLICATE_OF2' : 'category', 'TM_SET_1' : 'category', 'TM_SET_2' : 'category'})

    exelist.drop(exelist.index[[dellines]], inplace=True)

    exelist.to_csv(locations['SYSFILES']['H_scheduledalns'], sep="\t")

    # Show statistics
    for x in sorted(local_stats):
        print("STATS", "criterion", x, local_stats[x], local_stats[x])
    print("STATS", "criterion", "Finished", "time", time.ctime(), "\n")

    gc.collect()

    return exelist


def simple_comparison_lists(options, locations, str_data):
    """Returns the pairs to be compared, in form of a dict of dicts.
    For each structure there is a dict of structures, each key corresponding
    to a 4uple class-#TEs-pdbich1-pdbich2
    """

    eff_str_data = {x : str_data[x] for x in str_data if 'eliminated' not in str_data[x]['status']}

    # Get the number of TEs in each chain and cvreate list chain-#TEs
    pdbich_nte_list = []
    for pdbi in eff_str_data:
        for iTEs, ch in enumerate(eff_str_data[pdbi]['ENCOMPASS']['kTEs']):
            print(pdbi, ch, str_data[pdbi]['ENCOMPASS']['structure']['equivalent_chains'])
            if ch in str_data[pdbi]['ENCOMPASS']['structure']['equivalent_chains']:
                duplicate_of = str_data[pdbi]['ENCOMPASS']['structure']['equivalent_chains'][ch]
            else:
                duplicate_of = None
            TEs = eff_str_data[pdbi]['ENCOMPASS']['TEs'][iTEs]
            nte = 0
            tmlist = []
            for iTE, TE in enumerate(TEs['segments']):
                if TE["is_TM"] == 'True':
                    tmlist.append(iTE)
                nte += 1
            if tmlist:
                pdbich_nte_list.append((pdbi + '_' + ch, tmlist, duplicate_of))

    # Compile the list of pairs to be compared
    pre_exelist = []
    ind_exelist = []
    for pdbi_ch1, tml1, dupl1 in pdbich_nte_list:
        ntm1 = len(tml1)
        pdbi1 = pdbi_ch1[:4]
        ch1 = pdbi_ch1[5]
        ich1 = str_data[pdbi1]['ENCOMPASS']['kTEs'].index(ch1)
        TE_ch1 = str_data[pdbi1]['ENCOMPASS']['TEs'][ich1]
        for pdbi_ch2, tml2, dupl2 in pdbich_nte_list:
            if pdbi_ch2 == pdbi_ch1:
                print("DISCARD", pdbi_ch1, pdbi_ch2, "SAME")
                continue
            ntm2 = len(tml2)
            pdbi2 = pdbi_ch2[:4]
            ch2 = pdbi_ch2[5]
            ich2 = str_data[pdbi2]['ENCOMPASS']['kTEs'].index(ch2)
            TE_ch2 = str_data[pdbi2]['ENCOMPASS']['TEs'][ich2]
            p1p2 = pdbi_ch1 + '__' + pdbi_ch2
            # A-B comparisons: none
            if str_data[pdbi1]['ENCOMPASS']['class'] != str_data[pdbi2]['ENCOMPASS']['class']:
                print("DISCARD", pdbi_ch1, pdbi_ch2, "NOT THE SAME CLASS", str_data[pdbi1]['ENCOMPASS']['class'], str_data[pdbi2]['ENCOMPASS']['class'])
                continue
            # B-B comparisons: all
            if str_data[pdbi1]['ENCOMPASS']['class'] == 'beta':
                pre_exelist.append((pdbi1, ch1, pdbi2, ch2, True, True, ntm1, ntm2, dupl1, dupl2))
                ind_exelist.append(p1p2)
            # Ai-Aj (max(i,j)>3 and min(i,j)>3) comparisons: min/max >= 3/4 [i.e. 3/4, 4/5, 4/6,... but not 3/5, 3/6, 4/6, ...]
            elif min(ntm1, ntm2) > 2 and max(ntm1, ntm2) > 3:
                if max(ntm1, ntm2)*0.75 <= min(ntm1, ntm2):
                    pre_exelist.append((pdbi1, ch1, pdbi2, ch2, False, False, ntm1, ntm2, dupl1, dupl2))
                    ind_exelist.append(p1p2)
            # A3-A3 comparisons: all
            elif ntm1 == 3 and ntm2 == 3:
                pre_exelist.append((pdbi1, ch1, pdbi2, ch2, False, False, ntm1, ntm2, dupl1, dupl2))
                ind_exelist.append(p1p2)
            # A1-A1 comparisons: no big extrema or |max(1a,1b)-max(2a,2b)| < min(max(1a,1b), max(2a,2b))
            # (the difference between the biggest in each chain must not exceed the smallest of the two)
            elif ntm1 == 1 and ntm2 == 1:
                ich1 = str_data[pdbi1]['ENCOMPASS']['structure']['kchains'].index(ch1)
                ich2 = str_data[pdbi2]['ENCOMPASS']['structure']['kchains'].index(ch2)
                res_pdbc1 = [x for x in str_data[pdbi1]['ENCOMPASS']['structure']['chains'][ich1]['kresidues']]
                res_pdbc2 = [x for x in str_data[pdbi2]['ENCOMPASS']['structure']['chains'][ich2]['kresidues']]
                seg_pdbc1_s, seg_pdbc1_e = TE_ch1['segments'][tml1[0]]['residues_num_type3'][0][0], TE_ch1['segments'][tml1[0]]['residues_num_type3'][-1][0]
                seg_pdbc2_s, seg_pdbc2_e = TE_ch2['segments'][tml2[0]]['residues_num_type3'][0][0], TE_ch2['segments'][tml2[0]]['residues_num_type3'][-1][0]
                dom1a = len([x for x in res_pdbc1 if x[0] < seg_pdbc1_s])
                dom1b = len([x for x in res_pdbc1 if x[0] < seg_pdbc1_e])
                dom2a = len([x for x in res_pdbc2 if x[0] < seg_pdbc2_s])
                dom2b = len([x for x in res_pdbc2 if x[0] < seg_pdbc2_e])
                if (max(dom1a, dom1b) < 100 and max(dom2a, dom2b) < 100) or (abs(max(dom1a, dom1b) - max(dom2a, dom2b)) < min(max(dom1a, dom1b), max(dom2a, dom2b))):
                    pre_exelist.append((pdbi1, ch1, pdbi2, ch2, False, False, ntm1, ntm2, dupl1, dupl2))
                    ind_exelist.append(p1p2)
            # A2-A2 comparisons : no big extrema or (A1-A1 condition and (no big intra or similar intras))
            elif ntm1 == 2 and ntm2 == 2:
                res_pdbc1 = [x for x in str_data[pdbi1]['ENCOMPASS']['structure']['chains'][ich1]['kresidues']]
                res_pdbc2 = [x for x in str_data[pdbi2]['ENCOMPASS']['structure']['chains'][ich2]['kresidues']]
                seg_pdbc1_1s, seg_pdbc1_2e = TE_ch1['segments'][tml1[0]]['residues_num_type3'][0][0], TE_ch1['segments'][tml1[1]]['residues_num_type3'][-1][0]
                seg_pdbc2_1s, seg_pdbc2_2e = TE_ch2['segments'][tml2[0]]['residues_num_type3'][0][0], TE_ch2['segments'][tml2[1]]['residues_num_type3'][-1][0]
                seg_pdbc1_1e, seg_pdbc1_2s = TE_ch1['segments'][tml1[0]]['residues_num_type3'][-1][0], TE_ch1['segments'][tml1[1]]['residues_num_type3'][0][0]
                seg_pdbc2_1e, seg_pdbc2_2s = TE_ch2['segments'][tml2[0]]['residues_num_type3'][-1][0], TE_ch2['segments'][tml2[1]]['residues_num_type3'][0][0]
                dom1a = len([x for x in res_pdbc1 if x[0] < seg_pdbc1_1s])
                dom1b = len([x for x in res_pdbc1 if x[0] < seg_pdbc1_2e])
                dom2a = len([x for x in res_pdbc2 if x[0] < seg_pdbc2_1s])
                dom2b = len([x for x in res_pdbc2 if x[0] < seg_pdbc2_2e])
                dom1intra = len([x for x in res_pdbc1 if x[0] > seg_pdbc1_1e and x[0] < seg_pdbc1_2s])
                dom2intra = len([x for x in res_pdbc2 if x[0] > seg_pdbc2_1e and x[0] < seg_pdbc2_2s])
                if (
                        ((max(dom1a, dom1b) < 100 and max(dom2a, dom2b) < 100) 
                        or (abs(max(dom1a, dom1b) - max(dom2a, dom2b)) < min(max(dom1a, dom1b), max(dom2a, dom2b))))
                        and ((dom1intra < 100 and dom2intra < 100) 
                        or (abs(dom1intra - dom2intra) < min(dom1intra, dom2intra)))
                        ):
                    pre_exelist.append((pdbi1, ch1, pdbi2, ch2, False, False, ntm1, ntm2, dupl1, dupl2))
                    ind_exelist.append(p1p2)
            else:
                print("DISCARD", pdbi_ch1, pdbi_ch2, "NO OPTION")

    exelist = pd.DataFrame(pre_exelist, columns=['PDB1', 'CHAIN1', 'PDB2', 'CHAIN2', 'ISBETA1', 'ISBETA2', 'NTM1', 'NTM2', 'DUPLICATE_OF1', 'DUPLICATE_OF2', 'TM_SET_1', 'TM_SET_2'])
    exelist = exelist.astype({'PDB1' : 'category', 'CHAIN1' : 'category', 'PDB2' : 'category', 'CHAIN2' : 'category', 'ISBETA1' : 'bool', 'ISBETA2' : 'bool', 'NTM1' : 'uint8', 'NTM2' : 'uint8', 'DUPLICATE_OF1' : 'category', 'DUPLICATE_OF2' : 'category', 'TM_SET_1' : 'category', 'TM_SET_2' : 'category'})

    exelist.to_csv(locations['SYSFILES']['H_scheduledalns'], sep="\t")
    return exelist

def structure_alignment(options, locations, str_data, more_list=[], less_list=[], entry_exelist=None, only_table=False, only_gather=False, split_exelist=False):
    this_name = 'structure_alignment'

    less_list = list(set(less_list) - set(more_list))

    already_processed = []
    ex_list = {}

    # If a scheduling file was created, load it
    if type(entry_exelist) == type(None) and os.path.exists(locations['SYSFILES']['H_scheduledalns']):
        entry_exelist = pd.read_csv(locations['SYSFILES']['H_scheduledalns'], sep="\t", index_col=0, dtype={'PDB1' : 'category', 'CHAIN1' : 'category', 'PDB2' : 'category', 'CHAIN2' : 'category', 'ISBETA1' : 'bool', 'ISBETA2' : 'bool', 'NTM1' : 'uint8', 'NTM2' : 'uint8', 'DUPLICATE_OF1' : 'category', 'DUPLICATE_OF2' : 'category', 'TM_SET_1' : 'category', 'TM_SET_2' : 'category'})

    print("EXELIST")
    print(entry_exelist, len(entry_exelist))

    if split_exelist:
        exelists = []
        tmranges = [[1,2,3],[4,5,6],[x for x in range(7, 100)]]
        for tmrange in tmranges:
            exelists.append(entry_exelist.loc[entry_exelist['NTM1'].isin(tmrange)])
    else:
        exelists = [entry_exelist]

    is_first = True
    for ic, exelist in enumerate(exelists):
        batch_job_code = options['RUN'][('code', 'main_code')] + "FrTMAlign" + str(ic).zfill(3)
        output_dir = locations['FSYSPATH']['cache'] + 'output_{0}/'.format(batch_job_code)
        outputs = ['stats_<id>.txt', 'seq_seqalns_<id>.txt', 'str_seqalns_<id>.txt', 'stralns_<id>.ATOM.txt.gz'] 
 
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
            entrylist = pd.unique(exelist[(exelist["DUPLICATE_OF1"].isnull()) & (exelist["DUPLICATE_OF2"].isnull())][['PDB1', 'CHAIN1']].agg('_'.join, axis=1))
            entrylist_fn = output_dir + "/entry_list.txt"
            entrylist_f = open(entrylist_fn, 'w')
            for ip1, pdbi1 in enumerate(entrylist):
                pdbiset.add(pdbi1)
                entrylist_f.write("{0:10d}\t{1:6s}\n".format(ip1, pdbi1))
                # If only_gather, only fill pdbiset
                if only_gather:
                    continue
                locscr_filename = input_dir + "straln_exe_{0}.sh".format(ip1)
                with open(locscr_filename, 'w') as lsf:
                    lsf.write('PRESHARED="${1}"\n')
                    lsf.write('SHARED=`dirname ${PRESHARED}`\n')
                    lsf.write('pdb1={0}\n'.format(pdbi1))
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
                        lsf.write('    cp ${SHARED}/${pdb1}_enc.pdb ${SHARED}/${pdb2}_enc.pdb .\n')
                    lsf.write('    {0} ${{pdb1}}_enc.pdb ${{pdb2}}_enc.pdb -o prestraln_{1}.txt > out_tmp_{1}.txt\n'.format(frtm_path, ip1))
                    lsf.write("    sed -i 's/\\x0/X/g' out_tmp_{0}.txt\n".format(ip1))
                    lsf.write('    echo "INIT ${{pdb1}} ${{pdb2}}" >> stralns_{0}.txt\n'.format(ip1))
                    lsf.write('    cat prestraln_{0}.txt >> stralns_{0}.txt\n'.format(ip1))
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
                    lsf.write('grep "INIT\|ATOM\|TER" stralns_{0}.txt > stralns_{0}.ATOM.txt\n'.format(ip1))
                    lsf.write('gzip stralns_{0}.ATOM.txt\n'.format(ip1))
                    if options['ALL']['nodescratchpath']:
                        lsf.write('cp stralns_{0}.ATOM.txt.gz seq_seqalns_{0}.txt str_seqalns_{0}.txt stats_{0}.txt $THIS\n'.format(ip1))
                strlist_filename = input_dir + "strlist_{0}.txt".format(ip1)
                with open(strlist_filename, 'w') as strlsf:
                    sel = exelist[(exelist['PDB1'] == pdbi1[:4]) & (exelist['CHAIN1'] == pdbi1[5])]
                    for pdbi2 in pd.unique(sel[(sel["DUPLICATE_OF1"].isnull()) & (sel["DUPLICATE_OF2"].isnull())][['PDB2', 'CHAIN2']].agg('_'.join, axis=1)):
                        strlsf.write("{0}\n".format(pdbi2))
                        pdbiset.add(pdbi2)
   
            entrylist_f.close()
 
            # Copy pdb structures in input
            if not only_gather:
                for pdbi in pdbiset:
                    shutil.copyfile(locations['FSYSPATH']['chains'] + '{0}_enc.pdb'.format(pdbi), input_dir + '{0}_enc.pdb'.format(pdbi))
        
            # .3 Instruct locusts
            specific_inputs = ['strlist_<id>.txt', 'straln_exe_<id>.sh']
            command_template = 'cp <shared>pyscr . ; S=`dirname <shared>{0} `; bash straln_exe_<id>.sh ${{S}}'.format(entrylist[0])
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
    outputs = ['stats_<id>.txt', 'seqseqalns_<id>.txt', 'str_seqalns_<id>.txt', 'stralns_<id>.ATOM.txt.gz']
    for ic, exelist in enumerate(exelists):
        batch_job_code = options['RUN'][('code', 'main_code')] + "FrTMAlign" + str(ic).zfill(3)
        output_dir = locations['FSYSPATH']['cache'] + "output_" + batch_job_code + "/"
        #if not only_table:
        move_and_filter(locations, output_dir, exelist)
        table = make_new_table(locations, exelist, lengths, output_dir, outputs, append=(not is_first))
        is_first = False
    return table


def move_and_filter(locations, output_dir, exelist):

    # pdbi_ch -> [redundantchains]
    redundancy = list_redundant_chains(exelist)

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

        # Recover all redundant chains
        pdbilist = [pdbi_ch]
        if pdbi_ch in redundancy:
            for rch in redundancy[pdbi_ch]:
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
    

def list_redundant_chains(el):  # el is the exelist pandas framework
    redundancy = {}
    for line in el.itertuples():
        ind, pdb1, ch1, pdb2, ch2, isbeta1, isbeta2, ntm1, ntm2, d1, d2, _, _ = line
        if d1 == d1:  # excludes NaN values since NaN != NaN is True
            pdbc1 = pdb1 + "_" + d1
            if pdbc1 not in redundancy:
                redundancy[pdbc1] = []
            if ch1 not in redundancy[pdbc1]:
                redundancy[pdbc1].append(ch1)
        if d2 == d2:
            pdbc2 = pdb2 + "_" + d2
            if pdbc2 not in redundancy:
                redundancy[pdbc2] = []
            if ch2 not in redundancy[pdbc2]:
                redundancy[pdbc2].append(ch2)
    return redundancy


def make_new_table(locations, exelist, lengths, output_dir, outputs, equivalence={}, allow_incomplete=False, append=False):
    this_name = make_new_table.__name__

    print("OUTPUT", output_dir)
    summary_table_fname = locations['SYSFILES']['summarytable']

    subpaths = []
    with open(output_dir + 'output.log') as outf:
        for line in outf:
            if line.startswith("stats_") and line.split()[1] == "present":
                subpaths.append(line.split()[2])

    print('exelist selection that will be applied for creating summary_table.txt')
    selection = set(zip(exelist['PDB1'].astype(str) + "_" + exelist['CHAIN1'].astype(str), exelist['PDB2'].astype(str) + "_" + exelist['CHAIN2'].astype(str)))

    redexelist = exelist[(~exelist["DUPLICATE_OF1"].isnull()) | (~exelist["DUPLICATE_OF2"].isnull())]
    print("N OF REDUNDANT CHAINS", len(redexelist))

    redundancy = list_redundant_chains(redexelist)

    print("Create table")
    table = pd.DataFrame(columns=exelist.columns.tolist() + ["ALN_LEN", "SEQ_SEQID", "STR_SEQID", "RMSD", "TMSCORE"])
    table = table.astype({'PDB1' : 'category', 'CHAIN1' : 'category', 'PDB2' : 'category', 'CHAIN2' : 'category', 'ISBETA1' : 'bool', 'ISBETA2' : 'bool', 'NTM1' : 'uint8', 'NTM2' : 'uint8', 'DUPLICATE_OF1' : 'category', 'DUPLICATE_OF2' : 'category', 'TM_SET_1' : 'category', 'TM_SET_2' : 'category', "ALN_LEN" : "uint16", "SEQ_SEQID" : "float32", "STR_SEQID" : "float32", "RMSD" : "float32", "TMSCORE" : "float32"})

    print("fill table with straln data")
    not_completed = set()
    not_completed_unknown = 0
    pairs_with_alignments = set()
    mode = "a" if append else "w"
    fout = open(summary_table_fname, mode)
    for subpath in subpaths:
        print("SUBPATH", subpath)
        first_time_as_primary_pdb = True
        flag = False
        with open(subpath) as f:
            for line in f:
                fields = line.split()
                if len(fields) == 7:
                    pdbc1, pdbc2, alnlen, seqseqid, strseqid, rmsd, tmscore = fields
                    if flag:
                        print(f"CHECK: {pdbc1}-{pdbc2}")
                    pairs_with_alignments.add((pdbc1, pdbc2))
                    if (pdbc1, pdbc2) not in selection:
                        print((pdbc1, pdbc2), "will not be included in summary_table.txt")
                        continue
                else:
                    if len(fields) == 2:
                        pdbc1, pdbc2 = fields
                        pairs_with_alignments.add((pdbc1, pdbc2))
                        not_completed.add((pdbc1, pdbc2))
                    else:
                        not_completed_unknown += 1
                    continue
                pdbi1, ch1 = pdbc1.split("_")
                pdbi2, ch2 = pdbc2.split("_")
                if first_time_as_primary_pdb:
                    if pdbc1 in redundancy:
                        for rch1 in redundancy[pdbc1]:
                            fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(pdbi1, ch1, pdbi1, rch1, lengths[pdbc1], "1.0000", "1.0000", "0.0000", "1.0000"))
                            if flag:
                                print(f"CHECK: first_time_as_primary with {pdbc2}, wrote {pdbi1} {ch1} {pdbi1} {rch1}")
                            pairs_with_alignments.add((pdbc1, pdbi1 + "_" + rch1))
                            fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(pdbi1, rch1, pdbi1, ch1, lengths[pdbc1], "1.0000", "1.0000", "0.0000", "1.0000"))
                            if flag:
                                print(f"CHECK: first_time_as_primary with {pdbc2}, wrote {pdbi1} {rch1} {pdbi1} {ch1}")
                            pairs_with_alignments.add((pdbi1 + "_" + rch1, pdbc1))
                            for rch2 in redundancy[pdbc1]:
                                if rch1 != rch2:
                                    fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(pdbi1, rch1, pdbi1, rch2, lengths[pdbc1], "1.0000", "1.0000", "0.0000", "1.0000"))
                    first_time_as_primary_pdb = False
                fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(pdbi1, ch1, pdbi2, ch2, alnlen, seqseqid, strseqid, rmsd, tmscore))
                if flag:
                    print(f"CHECK: wrote {pdbi1} {ch1} {pdbi2} {ch2}")
                pairs_with_alignments.add((pdbc1, pdbc2))
                if pdbc2 in redundancy:
                    for rch2 in redundancy[pdbc2]:
                        fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(pdbi1, ch1, pdbi2, rch2, alnlen, seqseqid, strseqid, rmsd, tmscore))
                        pairs_with_alignments.add((pdbc1, pdbi2 + "_" + rch2))
                if pdbc1 in redundancy:
                    for rch1 in redundancy[pdbc1]:
                        fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(pdbi1, rch1, pdbi2, ch2, alnlen, seqseqid, strseqid, rmsd, tmscore))
                        pairs_with_alignments.add((pdbi1 + "_" + rch1, pdbc2))
                        if pdbc2 in redundancy:
                            for rch2 in redundancy[pdbc2]:
                                fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(pdbi1, rch1, pdbi2, rch2, alnlen, seqseqid, strseqid, rmsd, tmscore))
                                pairs_with_alignments.add((pdbi1 + "_" + rch1, pdbi2 + "_" + rch2))
    tmp_fn = locations['FSYSPATH']['main'] + "tmp.txt"

    for p in not_completed:
        print("Not completed:", p)
    print("Total KNOWN not completed", len(not_completed))
    print("Total UNKNOWN not completed", not_completed_unknown)
    print("Scheduled alignments that were omitted from the summary table (perhaps due to a difference between the scheduling table used for the run and for this analysis): ", selection.difference(pairs_with_alignments), " total: ", len(selection.difference(pairs_with_alignments)))

    return {}

def complete_straln(exelist_name, summary_name):
    exelist = pd.read_csv(exelist_name, sep="\t", index_col=0, dtype={'PDB1' : 'category', 'CHAIN1' : 'category', 'PDB2' : 'category', 'CHAIN2' : 'category', 'ISBETA1' : 'bool', 'ISBETA2' : 'bool', 'NTM1' : 'uint8', 'NTM2' : 'uint8', 'DUPLICATE_OF1' : 'category', 'DUPLICATE_OF2' : 'category', 'TM_SET_1' : 'category', 'TM_SET_2' : 'category'})
    exeset = set((exelist['PDB1']+exelist['CHAIN1']+exelist['PDB2']+exelist['CHAIN2']).unique())
    summarylist = pd.read_csv(summary_name, sep="\t", header=None, names=["PDB1", "CHAIN1", "PDB2", "CHAIN2", "ALN_LEN", "SEQ_SEQID", "STR_SEQID", "RMSD", "TMSCORE"])
    summaryset = set((summarylist['PDB1']+summarylist['CHAIN1']+summarylist['PDB2']+summarylist['CHAIN2']).unique())
    diffset = exeset - summaryset
    print("WARNING: Comparisons that are in the scheduled alignment but not the summary table (perhaps the scheduled alignment tabel changed between run and analysis): ",len(diffset), [x[0:4]+"_"+x[4:5]+"--"+x[5:9]+"_"+x[9:10] for x in diffset])


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
                    copy_flag = True
                else:
                    copy_flag = False
            if copy_flag:
                fo.write(line)
    if pdbi_ch in red:
        for pdbi_ch_add in red[pdbi_ch]:
            shutil.copyfile(out_fn, out_fn_template.replace("XXX", pdbi_ch_add))


if __name__ == "__main__":
    from supporting_functions import *
    from initialize_repository import *
    from combine_sources import *
    from complete_information import *

    options, locations = initialize_repository()

    str_data = read_checkpoint(locations['FSYSPATH']['cache'] + 'str_data_completegen.pkl', locations['SYSFILES']['data_structure_template'])

    print(">>>exelist")
    exelist = simpler_comparison_lists(options, locations, str_data)

    print(">>>straln")
    structure_alignment(options, locations, str_data, more_list=[], less_list=[], entry_exelist=exelist, only_table=True, only_gather=False, split_exelist=False)
    #complete_straln(exelist_name=locations['SYSFILES']['H_scheduledalns'], summary_name=locations['SYSFILES']['summarytable'])
