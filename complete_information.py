# Name: 
# Language: python3
# Libraries:
# Description: 
# Author: Edoardo Sarti
# Date: 

from supporting_functions import *

def calc_seqid_up(aln, bn, en, seq, unp_seq, fmt='biopython'):
    # aln: [intervals1, intervals2]
    # intervals: (interval1, interval2, ..., intervalN)
    # interval: (begin, end)
    # intervals1 and intervals2 contain the same number of intervals, and corresponding ones
    # have the same length (end - begin), and are in fact the matched-mismatched regions

    opts = [[x[1]-1 - x[0], x[1]-1 - bn, en - x[0], en - bn] for x in aln.aligned[1]]
    norm = sum([min([max(0,x+1) for x in y]) for y in opts])
    seqid = 0
    strb, stre = -1, -1
    if norm != 0:
        for i in range(len(aln.aligned[0])):
            ix, iy = aln.aligned[0][i], aln.aligned[1][i]
            if strb == -1:
                strb = ix[0]
            stre = ix[1]
            for j in range(ix[1] - ix[0]):
                if seq[ix[0]+j] == unp_seq[iy[0]+j] and iy[0]+j >= bn and iy[0]+j <= en:
                    seqid += 1
        seqid /= norm
    return seqid, strb, stre


def map_retrieve(ids2map, source_fmt='ACC+ID', target_fmt='ACC', output_fmt='tab'):
    '''
    Retrieves ACC from IDs
    Supported output formats:
    - 'tab' : returns a dictionary of ID -> ACC codes
    - 'fasta' : returns a dictionary of ID -> sequence
    '''

    BASE = 'http://www.uniprot.org'
    KB_ENDPOINT = '/uniprot/'
    TOOL_ENDPOINT = '/uploadlists/'
    lenids = len(ids2map)

    if output_fmt in ['tab', 'fasta']:
        nope = {}
    else:
        nope = []

    if type(ids2map) != list:
        ids2map = ids2map.split()

    REPEAT = 5
    rep = 0
    results = {'results' : {}}
    while (not results['results']) and rep<REPEAT: 
        job_id = uniprot_requests.submit_id_mapping(
            from_db="UniProtKB_AC-ID", to_db="UniProtKB", ids=ids2map
        )

        if uniprot_requests.check_id_mapping_results_ready(job_id):
            link = uniprot_requests.get_id_mapping_results_link(job_id)
            results = uniprot_requests.get_id_mapping_results_search(link)
        rep += 1

    d = {}
    if output_fmt == 'tab':
        for r in results['results']:
            d[r['from']] = r['to']['primaryAccession']
    elif output_fmt == 'fasta':
        for r in results['results']:
            if 'sequence' in r['to']:
                d[r['from']] = r['to']['sequence']['value']
    return d


def uniprot_sifts(locations, str_data):
    # SIFTS file must be done like this:
    # PDB CHAIN ACC BEGINRES_NUM ENDRES_NUM BEGINRES_PDBIDX ENDRES_PDBIDX BEGINRES_UNIPROT ENDRES_UNIPROT
    up_data = {}
    up_set = set()

    if not os.path.exists(locations['SYSFILES']['sifts_uniprot']):
        url = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz"
        local_filename = locations['SYSFILES']['sifts_uniprot']
        reiterated_gzip_download(url, local_filename)

    with open(locations['SYSFILES']['sifts_uniprot']) as f:
        for line in f:
            if line.strip().startswith("#"):
                continue
            pdbi, ch, acc, begnum, endnum, begpdb, endpdb, begp, endp = line.split()
            if pdbi in str_data:
                if pdbi not in up_data:
                    up_data[pdbi] = []
                if acc not in up_data[pdbi]:
                    up_data[pdbi].append((acc, int(begp)-1, int(endp)-1))
                up_set.add(acc)
    return up_data, up_set


def assign_uniprot_acc(locations, str_data, eq):
    """Needs:
    - biopython
    - str_data[pdbi]['ENCOMPASS']['structure']['chains'][ich]['sequence']
    - str_data[pdbi]['ENCOMPASS']['structure']['chains'][ich]['TM_regions']['TM_regions_extrema'] 
    """
    EXTR = 1000000

    # Statistics
    local_stats = {
        'no_conversion' : [],
        'up_codes' : [],
        'short_sequence' : [],
        'antibody' : [],
        'best_just_below_thr' : [],
        'no_up_codes' : [],
        'best_much_below_thr' : [],
        '#codes' : {},
        'no_up_TM' : [],
        'up_sifts_only' : [],
        'up_sifts&json' : []
    }

    # Searches in PDB json files and compiles a UniProt list dictionary
    #  that can also accommodate for starting and ending points of thee struct chain
    uniprot_codes = set()
    uniprot_dict = {} # UniProt list dictionary for 3-tuples (acc, init, end)   
    for pdbi in str_data:
        pdbjson_fn = locations['FSYSPATH']['PDBjsons'] + pdbi + '.json'
        uniprot_code_list = []
        if os.path.exists(pdbjson_fn):
            with codecs.open(pdbjson_fn, 'r+', encoding="utf-8") as json_file:
                js = json.load(json_file, encoding='utf-8')
            for x in js['polymer_entities']:
                ul = x['rcsb_polymer_entity_container_identifiers']['uniprot_ids']
                if ul:
                    uniprot_codes |= set(ul)
                    uniprot_code_list += ul

        # init and end are defaulted at some extreme value
        uniprot_dict[pdbi] = [(acc, -EXTR, EXTR) for acc in sorted(list(set(uniprot_code_list)))]

    # Look into the SIFTS database for completing the data
    up_data_sifts, up_set_sifts = uniprot_sifts(locations, str_data)

    # Merge the UniProt codes found
    uniprot_codes |= up_set_sifts

    # If the code has been found in SIFTS, add init and end
    for pdbi in up_data_sifts:
        if pdbi not in uniprot_dict:
            local_stats['up_sifts_only'].append(pdbi)
            uniprot_dict[pdbi] = up_data_sifts[pdbi]
        else:
            # If the code has been found in SIFTS, delete the current entry
            #  and replace it with the SIFTS one
            for acc, bn, en in up_data_sifts[pdbi]:
                if (acc, -EXTR, EXTR) in uniprot_dict[pdbi]:
                    local_stats['up_sifts&json'].append(pdbi)
                    delid = uniprot_dict[pdbi].index((acc, -EXTR, EXTR))
                    del uniprot_dict[pdbi][delid]
                if (acc, bn, en) not in uniprot_dict[pdbi]:
                    uniprot_dict[pdbi].append((acc, bn, en))  # Replace in any case

    # Retrieve the sequences corresponding to the uniprot codes
    binsize = 100
    unp_seqs = {} 
    for i in range(0, len(uniprot_codes)//binsize + 1):
        partial_l = sorted(list(uniprot_codes))[binsize*i:binsize*(i+1)]
        d = map_retrieve(partial_l, output_fmt='fasta')
        for e in d:
            unp_seqs[e] = d[e]

    aligner = init_pairwise_alignment()

    SEQID_THR = 0.95
    for pdbi in [x for x in str_data]:
        for ich, ch in enumerate(str_data[pdbi]['ENCOMPASS']['structure']['kchains']):
            resids = str_data[pdbi]['ENCOMPASS']['structure']['chains'][ich]['residues']
            seq = "".join([from3to1(x[1]) for x in resids])
            ticket = []
            if len(seq) <= 15:
                local_stats['short_sequence'].append((pdbi, ch))
            if 'antibody' in str_data[pdbi]['ENCOMPASS']['name'].lower()\
                    or 'antibody' in str_data[pdbi]['FROM_PDB']['name'].lower():
                local_stats['antibody'].append((pdbi, ch))

            # Couple the chain with pertinent UniProt accs
            #  If the alignment between the chain and the UniProt sequence has seqid>0.95,
            #  then add the acc to the list. If the list remains empty and the best seqid is >0.9,
            #  then add the corresponding acc.
            right_acc_name = [] # List of UniProt codes associated to a chain
            bestid, bestaln, bestname, bestextr, bestupextr = 0, [], '', (-1, -1), (-1, -1)
            status = {}
            for name, bn, en in uniprot_dict[pdbi]:
                try:
                    unp_seq = unp_seqs[name]  # name can not be in unp_seqs if seq has not been retrieved
                    result = aligner.align(seq, unp_seq)  # biopython can fail for whatever reason
                    if not result:  # when result is called, biopython can still fail for too many optimal solutions
                        status[pdbi+"_"+ch] = 'BIOPYTHON ERROR - NULL'
                        continue
                except:
                    if name not in unp_seqs:
                        status[pdbi+"_"+ch] = 'UNIPROT ERROR - NO SEQUENCE'
                    else:
                        status[pdbi+"_"+ch] = 'BIOPYTHON ERROR - FAILED'
                    continue
                if result:
                    aln = result[0]
                    seqid, strb, stre = calc_seqid_up(aln, bn, en, seq, unp_seq, fmt='biopython')
                if seqid > SEQID_THR:
                    right_acc_name.append((name, (bn, en), (resids[strb][0][0], resids[stre-1][0][0])))
                if seqid > bestid:
                    bestid = seqid
                    bestaln = aln
                    bestname = name
                    if strb == -1 or stre == -1:
                        bestextr = (-1, -1)
                    else:
                        bestextr = (resids[strb][0][0], resids[stre-1][0][0])
                    bestupextr = (bn, en)
            if not right_acc_name:
                if pdbi+"_"+ch not in status:
                    if bestid > 0.9:
                        right_acc_name.append((bestname, bestupextr, bestextr))
                        local_stats['best_just_below_thr'].append((pdbi, ch, bestname))
                    elif not uniprot_dict[pdbi]:
                        local_stats['no_up_codes'].append((pdbi, ch))
                    else:
                        local_stats['best_much_below_thr'].append((pdbi, ch, bestname))
            ran = ",".join([x[0] for x in right_acc_name])
            str_data[pdbi]['ENCOMPASS']['structure']['chains'][ich]['UniProt_acc'] = right_acc_name # This is a list
            if ran not in str_data[pdbi]['ENCOMPASS']['structure']['UniProt_stoichiometry']:
                str_data[pdbi]['ENCOMPASS']['structure']['UniProt_stoichiometry'][ran] = []
            str_data[pdbi]['ENCOMPASS']['structure']['UniProt_stoichiometry'][ran].append(ch)
            n = len(right_acc_name)
            if n not in local_stats['#codes']:
                local_stats['#codes'][n] = []
            local_stats['#codes'][n].append(pdbi)

    # If some equivalent chain has the UniProt code, transfer it
    # 1) Equivalence dictionary already present in variable eq
    # 2) Transfer UniProt codes between equivalent chains 
    # NEEDS TM analysis!
    for pdbi in sorted([x for x in str_data]):
        for ich, ch in enumerate(sorted(str_data[pdbi]['ENCOMPASS']['structure']['kchains'])):
            ENCchain = str_data[pdbi]['ENCOMPASS']['structure']['chains'][ich]
            if len(ENCchain['TM_regions']['TM_regions_extrema'])>0\
                    and (not ENCchain['UniProt_acc']):
                pdbi_ch1 = pdbi + '_' + ch
                if pdbi_ch1 not in eq:
                    continue
                for pdbi_ch2 in eq[pdbi_ch1]:
                    pdbi2, ch2 = pdbi_ch2.split('_')
                    if ch2 in str_data[pdbi2]['ENCOMPASS']['structure']['kchains']:
                        ich2 = str_data[pdbi2]['ENCOMPASS']['structure']['kchains'].index(ch2)
                    else:
                        continue
                    ENCchain2 = str_data[pdbi2]['ENCOMPASS']['structure']['chains'][ich2]
                    if ENCchain2['UniProt_acc']:
                        ENCchain['UniProt_acc'] = ENCchain2['UniProt_acc']
                        break
            if len(ENCchain['TM_regions']['TM_regions_extrema'])>0 and (not ENCchain['UniProt_acc']):
                local_stats['no_up_TM'].append((pdbi, ch))

    # Show statistics
    for x in sorted(local_stats):
        if type(local_stats[x]) is list:
            print("STATS", "assign_uniprot_acc", x, len(local_stats[x]), local_stats[x])
        elif type(local_stats[x]) is dict:
            for y in sorted(local_stats[x]):
                print("STATS", "assign_uniprot_acc", x, y, len(local_stats[x][y]), local_stats[x][y])
    print("STATS", "assign_uniprot_acc", "Finished", "time", time.ctime(), "\n")

    return str_data 


def aln_and_seqid(seq1, seq2, sn1, sn2, muscle_path='', cache_path=''):
    """Calculcates pairwise alignment and sequence identity
    Low-level function: independent of EncDS
    """

    if type(seq1) == list:
        seq1 = ''.join([x for x in seq1])
    if type(seq2) == list:
        seq2 = ''.join([x for x in seq2])

    # If MUSCLE will be used, cache_path must be specified
    if (muscle_path and not cache_path) or (cache_path and not muscle_path):
        print("ERROR1010")
        exit(1)

    # MUSCLE is the preferred method
    muscle_failed = False
    if muscle_path:
        # Cache name with random label
        signature = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
        muscle_in_filename = cache_path + 'muscle_in_' + signature + '.fa'
        muscle_out_filename = cache_path + 'muscle_out_' + signature + '.fa'
        muscle_stdout = cache_path + 'muscle_' + signature + '.out'
        muscle_stderr = cache_path + 'muscle_' + signature + '.err'
        # Write input fasta file
        with open(muscle_in_filename, 'w') as f:
            f.write('>{0}\n{1}\n>{2}\n{3}\n'.format(sn1, seq1, sn2, seq2))
        
        # Execute same MUSCLE job 3 times before giving up
        muscle_plist = muscle_path.split()
        i = 0
        while (not os.path.exists(muscle_out_filename)) and i<3:
            if i>0:
                print("Muscle output {0} not found. Try number {1}".format(sn1+":"+sn2, i+1))
            try:
                subprocess.check_call(muscle_plist + ['-in', muscle_in_filename, '-out', muscle_out_filename], stdout=open(muscle_stdout, 'w'), stderr=open(muscle_stderr, 'w'))
            except:
                print("Muscle failed on {0}".format(sn1+":"+sn2))
            i += 1 

        # Check if, after at most 3 tries, the output is there
        # If so, parse out put and create alignment
        if os.path.exists(muscle_out_filename):
            with open(muscle_out_filename) as f:
                aseqs = []
                aseq = ''
                for line in f:
                    if not line:
                        continue
                    if line.startswith('>'):
                        if aseq:
                            aseqs.append(aseq)
                        aseq = ''
                        continue
                    aseq += line.strip()
                if aseq:
                    aseqs.append(aseq)
                aln = list(zip(aseqs[0], aseqs[1]))
            os.remove(muscle_in_filename)
            os.remove(muscle_out_filename)
            os.remove(muscle_stdout)
            os.remove(muscle_stderr)
        else:
            muscle_failed = True

    # Calculate seqid
    norm = len([x for x in aln if '-' not in x])
    if norm == 0:
        seqid = 0
    else:
        seqid = len([x for x in aln if x[0] == x[1]])/norm

    return aln, seqid


def find_copies_of_chains(coord_dict_pdbi, analysis_pdbi, pdbi, pdb_data_pdbi={}, pdb_bioassembly_group=-1): #, thr_log_status="ERROR"):
	"""
	Checks if chains with coordinates are exact copies
	pdb_bioassembly_group : assembly number to consider (-1 equals "all")
	"""

	# annotate the copies in this way: [list_of_equals1, list_of_equals2, ...] 
	# Note: in case of an OPM sturcure, this data cannot be trusted. The only way is to do the seq_from_structure approach
	this_name = find_copies_of_chains.__name__

	if pdb_bioassembly_group != -1:
		if not pdb_data_pdbi:
			print("ERROR")
			exit(1)
		elif pdb_bioassembly_group < 0 or len(pdb_data_pdbi['biological_assemblies'].show_list(quiet=True)) <= pdb_bioassembly_group:
			print("ERROR2")
			print(pdb_bioassembly_group)
			print(pdb_data_pdbi['biological_assemblies'].show_list())
			exit(1)

	chains_to_consider = []
	if pdb_bioassembly_group != -1:
		for transf in pdb_data_pdbi['biological_assemblies'][pdb_bioassembly_group]['transformations']:
			for c in transf['chains']:
				if c not in chains_to_consider:
					chains_to_consider.append(c)
	else:
		chains_to_consider = [x for x in coord_dict_pdbi['COORDS']]

	if len(chains_to_consider) < 2:
		if len(chains_to_consider) == 0:
			print("ERROR2020")
			exit(1)
		return False, [chains_to_consider] # list of one list

	seqs = {}
	# Make a list of sequences that have to be compared. In case of a bioassembly, only take the ones mentioned in the bioassembly
	# if [1 for x in str_data_pdbi['PASSPORT'] if 'chain_noba' in x['location'] or 'chain_banoc' in x['location']]: (use this code for the "not in agreement" part)
	for c in chains_to_consider:
		seqs[c] = [from3to1(x[1]) for x in coord_dict_pdbi['COORDS'][c]]

	inds_to_consider = [analysis_pdbi['seqid_legend'].index(c) for c in chains_to_consider]
	seqids = np.array(analysis_pdbi['seqid_mtx'])

	# Make a list of lists of chains having a Seq. Id. >= 95%
	pairs = set()
	sn = pdbi + "_"
	for ic1, (c1, i1) in enumerate(list(zip(chains_to_consider, inds_to_consider))):
		for c2, i2 in list(zip(chains_to_consider, inds_to_consider))[ic1:]:
			if i1 > -1 and i2 > -1 and seqids[i1, i2] >= 0.95:
				pairs.add((c1, c2))

	copies_list = sorted(merge(pairs, sorted_lists=True))
	for c in sorted(list(seqs.keys())):
		c_found = False
		for cl in copies_list:
			if c in cl:
				c_found = True
				break
		if not c_found:
			for cli, cl in enumerate(copies_list):
				if sorted([c, cl[0]])[0] == c:
					copies_list = copies_list[:cli] + [c] + copies_list[cli:]
					break
	has_copies_of_chains = False
	for c in copies_list:
		if len(c) > 1:
			has_copies_of_chains = True
			break

	return has_copies_of_chains, copies_list


def generate_full_structure(str_data_pdbi, pdbi, coords={}, bioassembly=0): #, thr_log_status="ERROR"):
    this_name = generate_full_structure.__name__

    IDmtx = np.array([[1,0,0],[0,1,0],[0,0,1]])
    IDvec = np.array([0,0,0])

    # 1. Build the np matrices for rotation and translation
    rotmatrices = []
    tarrays = []
    translrots = str_data_pdbi['FROM_PDB']['translrots']
    ktranslrots = str_data_pdbi['FROM_PDB']['ktranslrots']
    for nmat in range(len(translrots)):
        rotmatrices.append(np.array(translrots[nmat]['matrix']))
        tarrays.append(np.array(translrots[nmat]['vector']))

    # 2. Generate structure 
    if not coords:
        coords, _ = retrieve_coords_from_CIF(str_data_pdbi['FROM_PDB']['cache_coordinate_file'])

    chain_sets = {}	# Sets of repeated chains (their union could be smaller than the set of all chains contained in the pdb)
    used_chains = set()	# Names of already used chains, and initialization
    transfs = str_data_pdbi['FROM_PDB']['biological_assemblies'][bioassembly]['transformations']
    for x in range(len(transfs)):
        for c in transfs[x]['chains']:
            used_chains.add(c)
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890'	# Alphabetical order used to assign names to new chains
    alphabetical = list(alphabet)
    for x in alphabet:
        alphabetical += ["{0}{1}".format(x,y) for y in alphabet]

    new_coords = {}
    for nop in range(len(transfs)):
        ops = transfs[nop]['operators']  # Operators are individual rototranslation functions. They can combine multiple rototranslators. Example: [(M1,v1),(M2,v2)] is the function v -> (M2*((M1*v)+v1))+v2
        for chain in [x for x in coords if x in transfs[nop]['chains']]:  # DANGER
            for iop in range(len(ops)): # each chain might be applied to multiple operators
                op = ops[iop]
                is_identity = False
                # If it is the identity operator, do not change chain names
                if len(op) == 1:
                    tr = op[0]
                    itr = ktranslrots.index(tr)
                    if (rotmatrices[itr] == IDmtx.tolist()).all() and (tarrays[itr] == IDvec.tolist()).all():  # DANGER
                        used_chains.add(chain)
                        if nop not in chain_sets:
                            chain_sets[(nop,iop)] = []
                        chain_sets[(nop,iop)].append(chain)
                        new_chain = chain
                        new_coords[new_chain] = coords[chain]
                        is_identity = True
                # If it isn't, change chain name
                if not is_identity:
                    for chid in alphabetical:
                        if chid not in used_chains and chid not in coords:
                            # Choose a new chain name
                            used_chains.add(chid)
                            if nop not in chain_sets:
                                chain_sets[(nop,iop)] = []
                            chain_sets[(nop,iop)].append(chid)
                            new_chain = chid

                            # Initialize the coords of the new chain
                            new_coords[new_chain] = {}
                            for resid_name in coords[chain]:
                                new_coords[new_chain][resid_name] = {}
                            break
                    # Apply operator, transf by transf
                    for tr in op: # each non-ID operator might be a combination of transformations
                        itr = ktranslrots.index(tr)
                        for resid_name in coords[chain]:
                            resid, resname = resid_name
                            if resname == "DUM":
                                continue
                            for atid in coords[chain][resid_name]:
                                new_coords[new_chain][resid_name][atid] = np.dot(rotmatrices[itr], coords[chain][resid_name][atid]) + tarrays[itr]
                        coords[new_chain] = new_coords[new_chain]
    return new_coords, chain_sets


def calc_consensus_PDB(str_data, pdbi, pdb_bioassembly_group=0):
    bioassembly = 0
    PDB_coords = retrieve_coords_from_CIF(str_data[pdbi]['FROM_PDB']['cache_coordinate_file'])
    has_copies_of_chains_PDB, copies_ch_PDB = find_copies_of_chains(PDB_coords, str_data[pdbi]['FROM_PDB']['analysis'], pdbi, pdb_data_pdbi=str_data[pdbi]['FROM_PDB'], pdb_bioassembly_group=pdb_bioassembly_group)
    transfs = str_data[pdbi]['FROM_PDB']['biological_assemblies'][bioassembly]['transformations']
    chain_signature_PDB = find_gcd([len(x) for x in copies_ch_PDB]) * len(transfs)
    return chain_signature_PDB


def calc_rotation_axis(rotmtx):
    axis_of_rotation = np.array([rotmtx[2][1] - rotmtx[1][2], rotmtx[0][2] - rotmtx[2][0], rotmtx[1][0] - rotmtx[0][1]])
    if np.linalg.norm(axis_of_rotation) == 0 and sum([rotmtx[x][x] for x in range(3)]) == -1:
        aor = []
        for x in range(3):
            if rotmtx[x][x] == 1:
                aor.append(1)
            else:
                aor.append(0)
        axis_of_rotation = np.array(aor)
    if np.linalg.norm(axis_of_rotation) > config.epsilon:
        axis_of_rotation = np.array(axis_of_rotation)/np.linalg.norm(axis_of_rotation)
    return axis_of_rotation


def select_biological_assembly(data):
    this_name = select_biological_assembly.__name__

    local_stats = {
        'add_identity_transf' : [],
        'sameseq' : [],
        '1chain0transf' : [],
        'PDB!=OPM' : [],
        'PDB?=OPM' : [],
        'weirdtransf' : [],
        'biomatrix_error' : [],
        'undecided' : [],
        'unsafe' : [],
        'safe' : [],
        'null_rotax' : []
    }

    # It ONLY takes pdbis for which the safety is UNDECIDED
    options, locations, str_data = data
    bioassembly = 0

    # Evaluates safety of OPM groups instead of single structures
    if options['PARAMETERS'][('', 'consensus_signature_PDB')]:
        consensus_signature_PDB = {}
        visited_pdbi = set()
        for pdbi in str_data:
            if pdbi in visited_pdbi:
                continue
            visited_pdbi.add(pdbi)
            pdb_signatures = [] # List for studying the consensus
            # Calculates PDB signature as: (greatest common divisor of copies of present chains) * (biomatrix transformations)
            # Example: a complex A,B,C,D,E,F where A,B are identical and C,D,E,F are identical has gcd = 2 (there are 2 full copies of the A,2*C unit)
            # If to that complex we add the non-identical chain G, theh gcd=1 (there is only one copy of the unit 2*A,4*C,G)
            if not ({'eliminated', 'pdb_eliminated', 'opm_eliminated'} & set(str_data[pdbi]['status'])):
                pdb_signatures.append(calc_consensus_PDB(str_data, pdbi, pdb_bioassembly_group=0))
            else:
                str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_jump', pdbi, "This entry has been deemed very problematic, and thus will not undergo comparative biological assembly analyses."))
            for sec_pdbi in [x['pdbid'] for x in str_data[pdbi]['FROM_OPM']['secondary_representations'] if x['pdbid'] in str_data]:
                visited_pdbi.add(sec_pdbi)
                if not ({'eliminated', 'pdb_eliminated'} & set(str_data[sec_pdbi]['status'])):
                    pdb_signatures.append(calc_consensus_PDB(str_data, sec_pdbi, pdb_bioassembly_group=0))
                else:
                    str_data[sec_pdbi]['PASSPORT'].append(passport_entry(this_name+'_jump2', pdbi, "This entry has been deemed very problematic, and thus will not undergo comparative biological assembly analyses."))
            if not pdb_signatures: # It means the entries were eliminated
                continue
    
            # Sort signatures by number of recurrences and form consensus
            sigd = {i:pdb_signatures.count(i) for i in pdb_signatures}
            invsigd = {}
            for i in sigd:
                if sigd[i] not in invsigd:
                    invsigd[sigd[i]] = []
                invsigd[sigd[i]].append(i)
            maxi = sorted([i for i in invsigd])[-1]
    
            # If two signatures are equally voted, then put a signature value that will be catched hereafter (-1)
            multiplicity = len(invsigd[maxi])
            if multiplicity > 1:
                val = -1
            else:
                val = invsigd[maxi][0]
            consensus_signature_PDB[pdbi] = val
    
            # Update consensus for the whole group
            for sec_pdbi in [x['pdbid'] for x in str_data[pdbi]['FROM_OPM']['secondary_representations'] if x['pdbid'] in str_data]:
                consensus_signature_PDB[sec_pdbi] = val # can be -1 !! Catch it!

    # Decision				
    decsn = {'SAFE' : {}, 'UNSAFE' : {}}
    status = {}	# Shortcut for str_data[pdbi]['status']
    passport = {} # Shortcut for str_data[pdbi]['PASSPORT']
    to_be_checked = False
    IDmtx = np.array([[1,0,0],[0,1,0],[0,0,1]])
    IDvec = np.array([0,0,0])
    for pdbi in str_data:
        if config.debug:
            print("CBA", pdbi)

        status[pdbi] = str_data[pdbi]['status']
        passport[pdbi] = str_data[pdbi]['PASSPORT']

        # If there is the PDB structure, the transformations are taken from there. Otherwise, it is the identity
        PDB_coords = retrieve_coords_from_CIF(str_data[pdbi]['FROM_PDB']['cache_coordinate_file'])
        PDB_seqs = get_sequence_from_CIF(str_data[pdbi]['FROM_PDB']['cache_coordinate_file'])
        OPMpdb_coords = retrieve_coords_from_CIF(str_data[pdbi]['FROM_OPM']['cache_coordinate_file'])
        
        # Preparation: fill in missing transformation sets
        # Case 1: broken PDB but working OPM. Transformations = identity
        if ('pdb_eliminated' in status[pdbi] and 'opm_eliminated' not in status[pdbi]) or\
                ('pdb_eliminated' not in status[pdbi] and len(str_data[pdbi]['FROM_PDB']['biological_assemblies'][bioassembly]['transformations']) == 0):
            from_pdb = str_data[pdbi]['FROM_PDB']

            trd = FixedDict(from_pdb['translrots'].get_fdict_template())
            trd['matrix'] = IDmtx.tolist()
            trd['vector'] = IDvec.tolist()
            from_pdb['translrots'].append(trd)
            from_pdb['ktranslrots'].append('1')

        # Get transformation set
        # Case 1: broken PDB but working OPM. Transformations = identity
        if ('pdb_eliminated' in status[pdbi] and 'opm_eliminated' not in status[pdbi]) or\
                ('pdb_eliminated' not in status[pdbi] and len(str_data[pdbi]['FROM_PDB']['biological_assemblies'][bioassembly]['transformations']) == 0):
            #has_copies_of_chains_OPM, copies_ch_OPM = find_copies_of_chains(OPMpdb_coords, str_data[pdbi]['FROM_OPM']['analysis'], pdbi, pdb_data_pdbi=str_data[pdbi]['FROM_PDB'], pdb_bioassembly_group=-1)

            from_pdb = str_data[pdbi]['FROM_PDB']

            trd = FixedDict(from_pdb['translrots'].get_fdict_template())
            trd['matrix'] = IDmtx.tolist()
            trd['vector'] = IDvec.tolist()
            from_pdb['translrots'].append(trd)
            from_pdb['ktranslrots'].append('1')

            bafd = FixedDict(from_pdb['biological_assemblies'].get_fdict_template())
            tfd = FixedDict(bafd['transformations'].get_fdict_template())
            tfd['chains'] = sorted([c for c in OPMpdb_coords])
            tfd['operator'] = ['1']
            bafd['transformations'].append(tfd)
            from_pdb['biological_assemblies'].append(bafd)
            passport[pdbi].append(passport_entry(this_name+'_topmc', pdbi, "This structure does not have reliable biological assembly records coming from the PDB database. The coordinate file coming from the OPM database will be used, along with the implicit identity transformation matrix."))

            local_stats['add_identity_transf'].append(pdbi)
        # Case 2: working PDB and broken OPM. Take transformations from PDB
        if 'opm_eliminated' in status[pdbi] and 'pdb_eliminated' not in status[pdbi]:
            passport[pdbi].append(passport_entry(this_name+'_tpdbc', pdbi, "This structure has reliable biological assembly records coming from the PDB database: the related transformation matrices will be used to generate the whole biological assembly"))

        # Safety decision
        PDB_conditions = {
            'sameseq' : None,
            '1chain0transf' : None,
            'PDB!=OPM' : None,
            'weirdtransf' : None
        }

        # Condition 1: in PDB, two or more chains with coordinates have the exact same sequence
        seqset = {PDB_seqs[ch] for ch in PDB_seqs}
        if len(seqset) < len(PDB_seqs):
            PDB_conditions['sameseq'] = True
            local_stats['sameseq'].append(pdbi)
        elif len(seqset) == len(PDB_seqs):
            PDB_conditions['sameseq'] = False
        else:
            print_log((
                "ERROR",
                this_name,
                ("sequence set: {0}\n"
                    "chains with coordinates list: {1}")\
                    .format("\n".join(list(seqset)), PDB_seqs)
            ))

        # Condition 2: in PDB, only 1 chain and the identity transf
        translrots = str_data[pdbi]['FROM_PDB']['translrots']
        #print(pdbi, [ch for ch in PDB_seqs])
        #print(pdbi, len(PDB_seqs), len(translrots), translrots[0]['matrix'], translrots[0]['vector'], translrots[0]['matrix'] == IDmtx.tolist(), translrots[0]['vector'] == IDvec.tolist())
        if len(PDB_seqs) == 1 and len(translrots) == 1 and translrots[0]['matrix'] == IDmtx.tolist() and translrots[0]['vector'] == IDvec.tolist():
            PDB_conditions['1chain0transf'] = True
            local_stats['1chain0transf'].append(pdbi)
        else:
            PDB_conditions['1chain0transf'] = False

        # Condition 3: OPM and PDB signatures differ
        # Case 1: both PDB and OPM working.
        # If the structure is from OPM, nothing in the structure will be changed or used for changing other structures
        # BUT it gets checked for consistency of copies of chains. If it has the same number of copies of all chains as the one suggested by PDB, ok, otherwise
        # we trust PDB better than OPM and we update the OPM status in monitored.
        if not ({'eliminated', 'pdb_eliminated', 'opm_eliminated'} & set(status[pdbi])):
            # The signature indicates the greatest common divisor among the repeated chains.
            # ex: the signature of [[A, B, C], [D, E, F], [G, H], [I, L]] is 1, of [[A, B, C, K], [D, E, F, J], [G, H], [I, L]] is 2, of [[A, B, C], [D, E, F], [G, H], [I]] is 1)
            # Ideally, the PDB signature should always be 1 * n_of_transformations
            # and the OPM signature should coincide with that, demonstrating an agreement on the BA.
            transfs = str_data[pdbi]['FROM_PDB']['biological_assemblies'][0]['transformations']

            # OPM signature
            OPMpdb_coords = retrieve_coords_from_CIF(str_data[pdbi]['FROM_OPM']['cache_coordinate_file'])
            has_copies_of_chains_OPM, copies_ch_OPM = find_copies_of_chains(OPMpdb_coords, str_data[pdbi]['FROM_OPM']['analysis'], pdbi, pdb_data_pdbi=str_data[pdbi]['FROM_PDB'], pdb_bioassembly_group=-1)
            chain_signature_OPM = find_gcd([len(x) for x in copies_ch_OPM])	

            # PDB signature
            # Catch problematic group signatures and revert, for this case, to individually-calculated signature
            if options['PARAMETERS'][('', 'consensus_signature_PDB')] and consensus_signature_PDB[pdbi] != -1:
                chain_signature_PDB = consensus_signature_PDB[pdbi]
                print("IF")
            else:
                PDB_coords = retrieve_coords_from_CIF(str_data[pdbi]['FROM_PDB']['cache_coordinate_file'])
                has_copies_of_chains_PDB, copies_ch_PDB = find_copies_of_chains(PDB_coords, str_data[pdbi]['FROM_PDB']['analysis'], pdbi, pdb_data_pdbi=str_data[pdbi]['FROM_PDB'], pdb_bioassembly_group=0)
#               transfs = str_data[pdbi]['FROM_PDB']['biological_assemblies'][0]['transformations']
                chain_signature_PDB = find_gcd([len(x) for x in copies_ch_PDB]) * len(transfs)
                print("ELSE")

            # OPM <-> PDB signatures
            sig_msg = ('NOTICE', this_name, "The structure {0} is characterized by the following chain signatures:\nPDB: {1}\nOPM: {2}".format(pdbi, chain_signature_PDB, chain_signature_OPM))
            print_log(sig_msg)
            if chain_signature_OPM != chain_signature_PDB:
                PDB_conditions['PDB!=OPM'] = True
                #status[pdbi].append('opm_monitored')
                #status[pdbi].append('pdb_unsafe')
                #passport[pdbi].append(passport_entry(this_name+'_chsigneq', pdbi, "This entry displays inconsistencies regarding its quaternary structure. The coordinate file coming from the OPM database is thus considered problematic, and the entry is classified as unreliable (it will not be used to amend other structures)".format(chain_signature_PDB, chain_signature_OPM)))
                passport[pdbi].append(passport_entry(this_name+'_chsigneq_followup', pdbi, "The structure displays two different subunit signatures: PDB {0}; OPM {1}.".format(chain_signature_PDB, chain_signature_OPM)))
            else:
                PDB_conditions['PDB!=OPM'] = False
        else:
            local_stats['PDB?=OPM'].append(pdbi)

        # Condition 4: Translation parallel to axis of rotation
        transfs = str_data[pdbi]['FROM_PDB']['biological_assemblies'][bioassembly]['transformations']
        t_opers = []
        for t in transfs:
            t_opers.append(t['operators'])
        trs = set()
        #print(t_opers)
        for ops in t_opers:
            for op in ops:
                trs |= set(op)

        delta = 0.1
        ktranslrots = str_data[pdbi]['FROM_PDB']['ktranslrots']
        for tr in trs:
            itr = ktranslrots.index(tr)
            if translrots[itr]['matrix'] == IDmtx.tolist():
                continue
            rotax = calc_rotation_axis(translrots[itr]['matrix'])
            if np.linalg.norm(rotax) < config.epsilon:
                local_stats['null_rotax'].append((pdbi, itr, translrots[itr]['matrix']))
                #print("NULL ROTATION AXIS", pdbi, itr, translrots[itr]['matrix'], rotax)
            trv = np.array(translrots[itr]['vector'])
            if not (trv==np.zeros(3)).all():
                trv /= np.linalg.norm(trv)
            if abs(np.dot(rotax, trv)) > 1 - delta:
                PDB_conditions['weirdtransf'] = True
                local_stats['weirdtransf'].append(pdbi)
                break
        if PDB_conditions['weirdtransf'] is None:
            PDB_conditions['weirdtransf'] = False

        # Decision on conditions 1-4
        is_safe = True
        if PDB_conditions['sameseq']: # None is caught before
            is_safe = False
            passport[pdbi].append(passport_entry(
                this_name+'_sameseq',
                pdbi,
                ("This entry's PDB coordinate file contains two or more "
                    "identical sequences.")
            ))
        else:
            passport[pdbi].append(passport_entry(
                this_name+'_nosameseq',
                pdbi,
                ("This entry's PDB coordinate file contains only unique "
                    "sequences")
            ))

        if PDB_conditions['1chain0transf']: # None is not possible
            is_safe = False
            passport[pdbi].append(passport_entry(
                this_name+'_1chain0transf',
                pdbi,
                ("This entry's PDB coordinate file contains only one "
                    "subunit and no bioassembly transformation "
                    "(excecpt the identity).")
            ))
        else:
            passport[pdbi].append(passport_entry(
                this_name+'_no1chain0transf',     
                pdbi,       
                ("This entry's PDB coordinate file contains either "
                    "multiple subunits of at least one nontrivial "
                    "biological assembly transformation.")
            ))

        if PDB_conditions['PDB!=OPM'] is None:
            passport[pdbi].append(passport_entry(
                this_name+'_skipPDBvsOPM',
                pdbi,
                ("The comparison between PDB and OPM biological assemblies "
                    "could not be done.")
            ))
            is_safe = None
        elif PDB_conditions['PDB!=OPM']:
            is_safe = False
            passport[pdbi].append(passport_entry(
                this_name+'_PDB!=OPM',
                pdbi,
                ("The biological assemblies in PDB and OPM records "
                "are not identical.")
            ))
        else:
            passport[pdbi].append(passport_entry(
                this_name+'_PDB==OPM',
                pdbi,
                ("The biological assemblies in PDB and OPM records "
                "are identical.")
            ))
        if PDB_conditions['weirdtransf']: # None is not possible
            is_safe = False
            passport[pdbi].append(passport_entry(
                this_name+'_weirdtransf',
                pdbi,
                ("The PDB biological assembly transformations contain at "
                    "least one translation that is parallel to the axis of "
                    "the associated rotation")
            ))
        else:
            passport[pdbi].append(passport_entry(
                this_name+'_noweirdtransf',
                pdbi,
                ("The PDB biological assembly transformations do not contain "
                    "translations that are parallel to the axis of "
                    "the associated rotation")
            ))
        if is_safe is None or is_safe:
            # Conndition 5: Detached subunits
            # Generate full structure
            full_coords_PDB, equivalent_chain_sets = generate_full_structure(str_data[pdbi], pdbi, coords=PDB_coords['COORDS'])
            passport[pdbi].append(passport_entry(
                this_name+'_geneqch', 
                pdbi, 
                ("The full biological assembly as described in the coordinate "
                    "file coming from the PDB database was generated from the "
                    "corresponding transformation matrices, resulting in the "
                    "following sets of equivalent chains: {0}")\
                    .format(equivalent_chain_sets)
            ))

        # Last check: if SAFE, then try to generate. Then mark as SAFE or UNSAFE
        # If structure wasn't deemed unsafe, it is because the OPM structure agrees with the PDB one, and thus we chose the OPM structure over the PDB one.
        # We now have to compile a list of safe structures to remedy the unsafe ones: we start from the OPM-chosen and we see if their full PDB structure is actually reliable.
        # If it isn't, that's ok, if it is, they are chosen as SAFE.
        if 'pdb_unsafe' not in status[pdbi] and 'pdb_safety_not_available' not in status[pdbi]:
            # Check on the generation of the full PDB structure, to see if the transformations are really ok
            if not PDB_coords:
                PDB_coords = retrieve_coords_from_CIF(str_data[pdbi]['FROM_PDB']['cache_coordinate_file'])
            # Generate full structure
            full_coords_PDB, equivalent_chain_sets = generate_full_structure(str_data[pdbi], pdbi, coords=PDB_coords['COORDS'])
            passport[pdbi].append(passport_entry(
                this_name+'_geneqch', 
                pdbi, 
                ("The full biological assembly as described in the coordinate "
                    "file coming from the PDB database was generated from the "
                    "corresponding transformation matrices, resulting in the "
                    "following sets of equivalent chains: {0}")\
                    .format(equivalent_chain_sets)
            ))


            # Check for errors
            pass_entry, biomatrix_error = check_PDB_biomatrix_errors(options, equivalent_chain_sets, full_coords_PDB, pdbi)
            
            # Conndition 5: Detached subunits
            if biomatrix_error:
                is_safe = False
                local_stats['biomatrix_error'].append(pdbi)
            passport[pdbi].append(pass_entry)

        if is_safe is None:
            status[pdbi].append('pdb_undecided')
            passport[pdbi].append(passport_entry(
                this_name+'_undecided',
                pdbi,
                ("This entry does not present biological assembly faults "
                "but due to incomplete records it will not be taken as "
                "a reference for correcting biological assembly defects "
                "in other entries")
            ))
            local_stats['undecided'].append(pdbi)
        elif not is_safe:
            has_copies_of_chains_PDB, copies_ch_PDB = find_copies_of_chains(PDB_coords, str_data[pdbi]['FROM_PDB']['analysis'], pdbi, pdb_data_pdbi=str_data[pdbi]['FROM_PDB'], pdb_bioassembly_group=0)
            decsn['UNSAFE'][pdbi] = [sorted(x)[0] for x in copies_ch_PDB]
            status[pdbi].append('pdb_unsafe')
            passport[pdbi].append(passport_entry(
                this_name+'_unsafe',
                pdbi,
                ("This entry presents biological assembly defects and will "
                    "undergo a correction procedure")
            ))
            local_stats['unsafe'].append(pdbi)
        else:
            decsn['SAFE'][pdbi] = str_data[pdbi]['FROM_PDB']
            status[pdbi].append('pdb_safe')
            passport[pdbi].append(passport_entry(
                this_name+'_safe',
                pdbi,
                ("This entry does not present biological assembly defects "
                "and its records are complete: it will be used as a reference "
                "for correcting biological assembly defects in other entries")
            ))
            local_stats['safe'].append(pdbi)

    return decsn, status, passport, [x for x in str_data], local_stats


def PDBTM_download(locations, pdbi):
	pdb_url = 'http://pdbtm.enzim.hu/data/database/{0}/{1}.trpdb.gz'.format(pdbi[1:3], pdbi)
	pdb_filename = locations['FSYSPATH']['PDBTMpdbs'] + pdbi + '_PDBTM.pdb'
	pdb_downloaded = reiterated_gzip_download(pdb_url, pdb_filename)

	xml_url = 'http://pdbtm.enzim.hu/data/database/{0}/{1}.xml'.format(pdbi[1:3], pdbi)
	xml_filename = locations['FSYSPATH']['PDBTMxmls'] + pdbi + '_PDBTM.xml'
	xml_downloaded = reiterated_simple_download(xml_url, xml_filename)

	return pdb_downloaded and xml_downloaded


def get_sequence_from_CIF(cache_path):
	this_name = get_sequence_from_CIF.__name__
	with open(cache_path) as cache_file:
		loop_descr = False
		descriptors = []
		column_res = -1
		column_ch = -1
		old_resid = -99999
		sequence = {}
		for line in cache_file:
			if not line.strip() or line.startswith("#"):
				continue
			if line.startswith("loop"):
				loop_descr = True
				continue
			if loop_descr:
				if line.startswith("_"):
					descriptors.append(line.strip())
				else:
					loop_descr = False
					if "_atom_site.label_comp_id" in descriptors:
						column_res = descriptors.index("_atom_site.label_comp_id")
					if "_atom_site.label_asym_id" in descriptors:
						column_ch = descriptors.index("_atom_site.label_asym_id")
					if "_atom_site.label_seq_id" in descriptors:
						column_resid = descriptors.index("_atom_site.label_seq_id")
				continue
			else:
				if column_res < 0 or column_ch < 0:
					continue
				resid = int(line.split()[column_resid])
				ch = line.split()[column_ch]
				if ch not in sequence:
					sequence[ch] = ''
					old_resid = -999999
				if resid > old_resid:
					sequence[ch] += from3to1(line.split()[column_res])
					old_resid = resid
	if not sequence:
		print_log((
			"ERROR",
			this_name,
			"No sequences found in {0}".format(cache_path)
		))
	for ch in sequence:
		if not sequence[ch]:
			print_log((
				"ERROR",
				this_name,
				"Empty sequence for chain {0} in {1}".format(ch, cache_path)
			))
	return sequence

def apply_rotation(c_dict, R):
	new_c_dict = {'COORDS' : {}, 'DUMMIES' : []}
	for ch in c_dict['COORDS']:
		new_c_dict['COORDS'][ch] = {}
		for resid_name in c_dict['COORDS'][ch]:
			new_c_dict['COORDS'][ch][resid_name] = {}
			for at in c_dict['COORDS'][ch][resid_name]:
				new_c_dict['COORDS'][ch][resid_name][at] = np.dot(R, c_dict['COORDS'][ch][resid_name][at])
	for at, c in c_dict['DUMMIES']:
		new_c_dict['DUMMIES'].append((at, tuple(np.dot(R, c))))

	return new_c_dict


def add_DUM_to_coord_dict(c_dict, z):
	gcx, gcy, n = 0, 0, 0
	for ch in c_dict['COORDS']:
#		print("LEN SUBD", len(c_dict['COORDS'][ch]))
		for resid_name in c_dict['COORDS'][ch]:
			if 'CA' not in c_dict['COORDS'][ch][resid_name]:
				print("CA not there", c_dict['COORDS'][ch][resid_name])
				continue
			cx, cy, cz = c_dict['COORDS'][ch][resid_name]['CA']
			gcx, gcy = gcx + cx, gcy + cy
			n += 1
	gcx, gcy = gcx/n, gcy/n

	dmax = 0
	for ch in c_dict['COORDS']:
		for resid_name in c_dict['COORDS'][ch]:
			if 'CA' not in c_dict['COORDS'][ch][resid_name]:
				continue
			cx, cy, cz = c_dict['COORDS'][ch][resid_name]['CA']
			d = math.sqrt((gcx - cx)**2 + (gcy - cy)**2)
			if d > dmax:
				dmax = d
	
	for i in range(int(gcx-dmax-10), int(gcx+dmax+11), 2):
		for j in range(int(gcy-dmax-10), int(gcy+dmax+11), 2):
			d = math.sqrt((gcx-i)**2 + (gcy-j)**2)
			if d >= dmax+10:
				continue
			c_dict['DUMMIES'].append(("O", (i, j, z)))
			c_dict['DUMMIES'].append(("N", (i, j, -z)))
	return c_dict


def PDBTM_to_ENC(PDBTM_dict, str_data_pdbi, pdbi, cache_filename=''):
	this_name = PDBTM_to_ENC.__name__
	x, y, z = -9999, -9999, -9999
	with open(str_data_pdbi['FROM_PDBTM']['xml_file']) as xml_f:
		for line in xml_f:
			fields = line.split()
			if fields[0] == '<NORMAL':
				x = float(fields[1].split('=')[1][1:-1])
				y = float(fields[2].split('=')[1][1:-1])
				z = float(fields[3].split('=')[1][1:-3])
				break
	if x == -9999:
		print("NORMAL NOT FOUND")
		str_data_pdbi['PASSPORT'].append(passport_entry(this_name+'_rel', pdbi, "This structure as recorded in the coordinate file from the PDBTM databaase is not provided with information regarding its insertion in the model lipid bilayer."))
		if str_data_pdbi['FROM_PDBTM']['cache_coordinate_file']:
			write_mmCIF(str_data_pdbi, PDBTM_dict, str_data_pdbi['FROM_PDBTM']['cache_coordinate_file'], LOC_pad)
			str_data_pdbi['status'].append('pdbtm_approved')
			str_data_pdbi['status'].append('pdbtm_monitored')
		else:
			str_data_pdbi['status'].append('pdbtm_eliminated')
		return PDBTM_dict, str_data_pdbi

	if abs(x) > 0.001 or abs(y) > 0.001:
		print("NORMAL NOT ALIGNED ON Z", x, y, z)

	normal = np.array([x, y, z])
	z = np.linalg.norm(normal)
	angle, axis, R = find_rotation(normal, np.array([0,0,z]))
	PDBTM_dict = apply_rotation(PDBTM_dict, R)
	PDBTM_dict = add_DUM_to_coord_dict(PDBTM_dict, z)

	if cache_filename:
		write_mmCIF(str_data_pdbi, PDBTM_dict, cache_filename, LOC_pad)
		str_data_pdbi['present_cache'] = cache_filename
	return PDBTM_dict, str_data_pdbi
 

def build_by_comparison(data):
	options, locations, str_data_pdbi, pdbi = data
	scan_data = (options, locations, str_data_pdbi)

	str_data_pdbi['FROM_PDBTM']['coordinate_file'] = locations['FSYSPATH']['PDBTMpdbs'] + pdbi + '_PDBTM.pdb'
	str_data_pdbi['FROM_PDBTM']['cache_coordinate_file'] = locations['FSYSPATH']['cache'] + pdbi + '_PDBTM.cif'
	str_data_pdbi['FROM_PDBTM']['xml_file'] = locations['FSYSPATH']['PDBTMxmls'] + pdbi + '_PDBTM.xml'

	downloaded = PDBTM_download(locations, pdbi)
	if downloaded:
		PDBTM_dict, PDBTM_report = parse_PDB_coords(str_data_pdbi['FROM_PDBTM']['coordinate_file'], scan=True, options=options) 
		if PDBTM_dict:
			str_data_pdbi = check_version(options, [PDBTM_dict, PDBTM_report, 'pdbtm'], str_data_pdbi, pdbi)
			PDBTM_dict, str_data_pdbi = PDBTM_to_ENC(PDBTM_dict, str_data_pdbi, pdbi, cache_filename=str_data_pdbi['FROM_PDBTM']['cache_coordinate_file'])
		else:
			str_data_pdbi['PASSPORT'].append(passport_entry(this_name+'_rel', pdbi, "The structural and/or format inconsistencies of this entry could not be resolved with the help of the coordinate file from the PDBTM database. This entry is scheduled for deletion"))
			str_data_pdbi['status'].append('eliminated')
			str_data[pdbi]['delete_keyword'] = 'Unclear'

	return str_data_pdbi, [pdbi]

def decision(locations, str_data):
	this_name = decision.__name__

	LOC_pad = locations['SYSFILES']['pdbatomdict']

	run_in_PPM = set()
	n_opm, n_pdb, n_pdbtm, n_pdbtm_ppm, n_del = 0, 0, 0, 0, 0
	for pdbi in str_data:
		print("DECISION", pdbi)
		print(str_data[pdbi]['status'])
		# If eliminated, there should already be the ENC file with the comments. Is it there??
		if 'eliminated' in str_data[pdbi]['status']:
			n_del += 1
			print_log(('NOTICE', this_name, "del {0}".format(pdbi)))
			str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_strdel', pdbi, "This entry will not be reflected in the EncoMPASS database"))
			write_ENC(str_data[pdbi], {'COORDS' : {}, 'DUMMIES' : {}}, locations['FSYSPATH']['ENCpdbs'] + pdbi + '_enc.pdb', locations['SYSFILES']['ignoredpasscodes'], LOC_pad)
			continue
		if 'pdbtm_approved' in str_data[pdbi]['status'] and 'pdbtm_eliminated' not in str_data[pdbi]['status']:
			str_data[pdbi]['present_cache'] = str_data[pdbi]['FROM_PDBTM']['cache_coordinate_file']
			run_in_PPM.add(pdbi)
			n_pdbtm_ppm += 1
			print_log(('NOTICE', this_name, "pdbtm -> PPM {0}".format(pdbi)))
			str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_pdbtminppm', pdbi, "This structure reflects its description in the PDBTM database, and will be analyzed with the PPM software in order to insert it in a model lipid bilayer"))
			str_data[pdbi]['status'].append('run_in_PPM')
		# If OPM ok, then go on with it
		elif 'opm_eliminated' not in str_data[pdbi]['status'] and 'opm_monitored' not in str_data[pdbi]['status']:
			str_data[pdbi]['present_cache'] = str_data[pdbi]['FROM_OPM']['cache_coordinate_file']
			c_dict_pdbi = retrieve_coords_from_CIF(str_data[pdbi]['present_cache'])
			str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_strok', pdbi, "This structure reflects its description in the OPM database, and has passed all structural checks"))
			write_ENC(str_data[pdbi], c_dict_pdbi, str_data[pdbi]['ENCOMPASS']['coordinate_file'], locations['SYSFILES']['ignoredpasscodes'], LOC_pad)
			n_opm += 1
			print_log(('NOTICE', this_name, "opm {0}".format(pdbi)))
		# If not, if PDB ok, go on and add to PPM
		elif 'pdb_eliminated' not in str_data[pdbi]['status'] and 'pdb_monitored' not in str_data[pdbi]['status']: #and 'pdb_unsafe' not in str_data[pdbi]['status']:
			str_data[pdbi]['present_cache'] = str_data[pdbi]['FROM_PDB']['cache_coordinate_file']
			run_in_PPM.add(pdbi)
			n_pdb += 1
			print_log(('NOTICE', this_name, "pdb -> PPM {0}".format(pdbi)))
			str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_runinppm', pdbi, "This structure reflects its description in the PDB database, and will be analyzed with the PPM software in order to insert it in a model lipid bilayer"))
			str_data[pdbi]['status'].append('run_in_PPM')
		elif 'virus' in str_data[pdbi]['FROM_OPM']['superfamily']['name'].lower() or 'viral' in str_data[pdbi]['FROM_OPM']['superfamily']['name'].lower()\
                    or 'virus' in str_data[pdbi]['FROM_OPM']['name'].lower() or 'viral' in str_data[pdbi]['FROM_OPM']['name'].lower()\
                    or 'virus' in str_data[pdbi]['FROM_PDB']['name'].lower() or 'virus' in str_data[pdbi]['FROM_PDB']['title'].lower()\
                    or 'viral' in str_data[pdbi]['FROM_PDB']['name'].lower() or 'viral' in str_data[pdbi]['FROM_PDB']['title'].lower():
			str_data[pdbi]['status'].append('opm_rescued')
			str_data[pdbi]['present_cache'] = str_data[pdbi]['FROM_OPM']['cache_coordinate_file']
			c_dict_pdbi = retrieve_coords_from_CIF(str_data[pdbi]['present_cache'])
			str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_isviral', pdbi, "This is a viral membrane protein, and constitutes an exception to the PDB/OPM consistency criteria. The OPM description is chosen"))
			write_ENC(str_data[pdbi], c_dict_pdbi, str_data[pdbi]['ENCOMPASS']['coordinate_file'], locations['SYSFILES']['ignoredpasscodes'], LOC_pad)
			n_opm += 1
			print_log(('NOTICE', this_name, "opm {0}".format(pdbi)))
		elif 'pdb_monitored' in str_data[pdbi]['status']:
			str_data[pdbi]['status'].append('opm_rescued')
			str_data[pdbi]['present_cache'] = str_data[pdbi]['FROM_OPM']['cache_coordinate_file']
			c_dict_pdbi = retrieve_coords_from_CIF(str_data[pdbi]['present_cache'])
			str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_pdbsus', pdbi, "The structure is hardly classifiable. Because of inconsistencies in its PDB version, the OPM version will be considered viable"))
			write_ENC(str_data[pdbi], c_dict_pdbi, str_data[pdbi]['ENCOMPASS']['coordinate_file'], locations['SYSFILES']['ignoredpasscodes'], LOC_pad)
			n_opm += 1
			print_log(('NOTICE', this_name, "opm {0}".format(pdbi)))
		else:
			print_log(('NOTICE', this_name, "unc {0}".format(pdbi)))
			str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_strunc', pdbi, "This structure has defied any previous classification. It will be manually assessed."))
	n_unc = len(str_data) - n_opm - n_pdb - n_pdbtm - n_pdbtm_ppm - n_del

	print_log(('NOTICE', this_name, "Statistics: from OPM: {0}, from PDB (to be inserted in membrane): {1}, from PDBTM: {2}, from PDBTM (to be reinserted in membrane): {3}, deleted: {4}, unclassified: {5}, TOTAL: {6}".format(n_opm, n_pdb, n_pdbtm, n_pdbtm_ppm, n_del, n_unc, len(str_data))))

	return str_data, run_in_PPM


def cleanup(locations, str_data):
	this_name = cleanup.__name__

	for pdbi in str_data:
		if 'eliminated' in str_data[pdbi]['status']:
			str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_cleanupdel', pdbi, "This entry will not be reflected in the EncoMPASS database"))
			write_ENC(str_data[pdbi], {}, "", locations['SYSFILES']['ignoredpasscodes'])	
		if not str_data[pdbi]['delete_keyword']:
			str_data[pdbi]['delete_keyword'] = 'Unknown'
	return


def parallel_generate_full_structure(options, locations, str_data,
        generate_only=None, thr_log_status="ERROR", cycle_0=True, cycle_1=True):
    """Generates the complete EncoMPASS-grade structure starting from the 
    PDB structure. (The OPM structure is already ok).
    cycle_1 : Service option for skipping MUSCLE and biological assembly
    """

    this_name = parallel_generate_full_structure.__name__

    # The actual set of structures is the one that has not been already eliminated
    str_set = {pdbi for pdbi in str_data if 'eliminated' not in str_data[pdbi]['status']}

    # Run MUSCLE sequence alignments
    if cycle_0: # Service option
        str_data = run_MUSCLE(options, locations, str_data, str_set, only_analysis=False) # False in production
        write_checkpoint(str_data, locations['FSYSPATH']['cache'] + 'str_data_pre2a.pkl')
    else:	
        str_data = read_checkpoint(locations['FSYSPATH']['cache'] + 'str_data_pre2a.pkl', 
            locations['SYSFILES']['data_structure_template'])
	
    # Step 1
    biological_assembly = {
        'SAFE' : {}, # contains str_data entries
        'UNSAFE' : {} # contains only pdb names and the selected chain names
    }

    # Step 1a: prepare independent inputs
    data = []
    pre_data = []
    for pdbi in str_data:
        status = str_data[pdbi]['status']
        if 'eliminated' in status or 'pdb_eliminated' in status:
            continue
        # Simple version
        if pdbi in str_set:
            if options['PARAMETERS'][('', 'consensus_signature_PDB')]: # Option for grouping judgements over OPM representation groups 
                pre_data.append(pdbi)
            else:
                data.append((options, locations, str_data[pdbi]))

    # Data is submitted in batches corresponding to OPM representative/related groups
    if options['PARAMETERS'][('', 'consensus_signature_PDB')]:
        ids_list, secrep_d = read_OPM_representation_chart(locations['SYSFILES']['OPMreprchart'])
        for pdbi in secrep_d:
            group_data = {}
            if pdbi not in pre_data:
                continue
            group_data[pdbi] = str_data[pdbi]
            for sec_pdbi in secrep_d[pdbi]:
                if sec_pdbi not in pre_data:
                    continue
                group_data[sec_pdbi] = str_data[sec_pdbi]
            data.append((options, locations, group_data))

    # Cycle 1: Evaluate biological assembly of those entries that are neither SAFE nor UNSAFE
    # At their first appearence in EncoMPASS, all structures lack these keywords

    numproc = int(options['RUN'][('np', 'number_of_processors')])	# Option used ONLY for SINGLE machines (no clusters!)	
    pool_outputs = []
    if cycle_1:
        for d in data:
            o = select_biological_assembly(d)
            pool_outputs.append(o)
	
        # Step 3a: collect outputs
        # The routine can also delete entries in OPM_data, but we have to retain the data on OPM-parsed structures
        local_stats = {}
        for s in pool_outputs:
            bioassembly_chunk, status_chunk, passport_chunk, structs_chunk, ls = s
            for k in ls:
                if k not in local_stats:
                    local_stats[k] = []
                local_stats[k] += ls[k]
            if not bioassembly_chunk:
                continue
            for k in biological_assembly:	# over SAFE, UNSAFE
                for kk in bioassembly_chunk[k]: # over pdbi 
                    biological_assembly[k][kk] = bioassembly_chunk[k][kk]
            for pdbi in structs_chunk:
                if pdbi in status_chunk: 
                    if (not generate_only or pdbi in generate_only):
                        str_data[pdbi]['status'] = status_chunk[pdbi]
                    else:
                        if pdbi in str_data:
                            str_data[pdbi]['status'].append('eliminated')
                            str_data[pdbi]['delete_keyword'] = 'Unclear'
            for pdbi in structs_chunk:
                if pdbi in passport_chunk:
                    str_data[pdbi]['PASSPORT'] = passport_chunk[pdbi]
        # Show statistics
        for x in sorted(local_stats):
            if type(local_stats[x]) is list:
                print("STATS", "select_biological_assembly", x, len(local_stats[x]), local_stats[x])
            elif type(local_stats[x]) is dict:
                for y in sorted(local_stats[x]):
                    print("STATS", "select_biological_assembly", x, y, len(local_stats[x][y]), local_stats[x][y])
        print("STATS", "select_biological_assembly", "Finished", "time", time.ctime(), "\n")

	
        write_checkpoint(str_data, locations['FSYSPATH']['cache'] + 'str_data_2a.pkl')

    str_data = read_checkpoint(locations['FSYSPATH']['cache'] + 'str_data_2a.pkl', locations['SYSFILES']['data_structure_template']) 

    LOC_pad = locations['SYSFILES']['pdbatomdict']

    local_stats = {
        'unsafe_and_opmelim' : [],
        'unsafe_and_opmgood' : [],
        'unsafe' : [],
        'pdbtm_download_ok' : [],
        'pdbtm_parsed_ok' : [],
        'pdbtm_no_parsed' : [],
        'pdbtm_download_no' : []
    }

    # Cycle Abis: remedy unsafe set with OPM or PDBTM
    for pdbi in biological_assembly['UNSAFE']:
        local_stats['unsafe'].append(pdbi)
        if 'opm_eliminated' not in str_data[pdbi]['status']:
            str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_pass', pdbi, "Given the previous warnings over possibile inconsistencies of the PDB structure, the OPM structure is chosen for representing this entry"))
            local_stats['unsafe_and_opmgood'].append(pdbi)
        else:
            local_stats['unsafe_and_opmelim'].append(pdbi)
            str_data[pdbi]['FROM_PDBTM']['coordinate_file'] = locations['FSYSPATH']['PDBTMpdbs'] + pdbi + '_PDBTM.pdb'
            str_data[pdbi]['FROM_PDBTM']['cache_coordinate_file'] = locations['FSYSPATH']['cache'] + pdbi + '_PDBTM.cif'
            str_data[pdbi]['FROM_PDBTM']['xml_file'] = locations['FSYSPATH']['PDBTMxmls'] + pdbi + '_PDBTM.xml'
    
            downloaded = PDBTM_download(locations, pdbi)
            if downloaded:
                local_stats['pdbtm_download_ok'].append(pdbi)
                PDBTM_dict, PDBTM_report = parse_PDB_coords(str_data[pdbi]['FROM_PDBTM']['coordinate_file'], scan=True, options=options) 
                if PDBTM_dict:
                    local_stats['pdbtm_parsed_ok'].append(pdbi)
                    str_data[pdbi] = check_version(options, [PDBTM_dict, PDBTM_report, 'pdbtm'], str_data[pdbi], pdbi)
                    write_mmCIF(str_data[pdbi], PDBTM_dict, str_data[pdbi]['FROM_PDBTM']['cache_coordinate_file'], LOC_pad)
                    str_data[pdbi]['status'].append('pdbtm_approved')
                else:
                    local_stats['pdbtm_no_parsed'].append(pdbi)
                    str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_rel', pdbi, "The structural and/or format inconsistencies of this entry could not be resolved with the help of the coordinate file from the PDBTM database. This entry is scheduled for deletion"))
                    str_data[pdbi]['status'].append('eliminated')
                    str_data[pdbi]['delete_keyword'] = 'Unclear'
            else:
                local_stats['pdbtm_download_no'].append(pdbi)

    # Show statistics
    for x in sorted(local_stats):
        if type(local_stats[x]) is list:
            print("STATS", "select_biological_assembly_BIS", x, len(local_stats[x]), local_stats[x])
        elif type(local_stats[x]) is dict:
            for y in sorted(local_stats[x]):
                print("STATS", "select_biological_assembly_BIS", x, y, len(local_stats[x][y]), local_stats[x][y])
    print("STATS", "select_biological_assembly_BIS", "Finished", "time", time.ctime(), "\n")
            
    str_d, run_in_PPM = decision(locations, str_data)
    return str_data, run_in_PPM
		

def find_rotation(v_in, v_out):
	I = np.zeros((3,3))
	I[0,0] = 1
	I[1,1] = 1
	I[2,2] = 1

	v_in = v_in/np.linalg.norm(v_in)
	v_out = v_out/np.linalg.norm(v_out)
	axis_of_rotation = np.cross(v_in, v_out)
	if np.linalg.norm(axis_of_rotation) < 0.01:
		return 0, [0, 0, 0], I
	angle_of_rotation = np.arcsin(np.linalg.norm(axis_of_rotation))
	axis_of_rotation = axis_of_rotation/np.linalg.norm(axis_of_rotation)

	Z = np.zeros((3,3))
	Z[0,1] = - axis_of_rotation[2]
	Z[0,2] = axis_of_rotation[1]
	Z[1,2] = - axis_of_rotation[0]
	Z[1,0] = - Z[0,1]
	Z[2,0] = - Z[0,2]
	Z[2,1] = - Z[1,2]

	I = np.zeros((3,3))
	I[0,0] = 1
	I[1,1] = 1
	I[2,2] = 1

	R = I + math.sin(angle_of_rotation)*Z + (1 - math.cos(angle_of_rotation))*np.linalg.matrix_power(Z,2)

	return angle_of_rotation, axis_of_rotation, R


def check_PDB_biomatrix_errors(options, equivalent_chain_sets, coords, pdbi, thr_log_status="ERROR"):
	this_name = check_PDB_biomatrix_errors.__name__

	print(pdbi, equivalent_chain_sets)
	report = scan_struct(options, coords=coords)

	# Check also if there are clashes. If clashes, error
	chain_sets = equivalent_chain_sets
	detached_units = []
	error = False
	if len(chain_sets) > 1:
		for cs1 in chain_sets:
			contact_found = False
			for cs2 in chain_sets:
				if cs1 == cs2:
					continue
				for cid1 in chain_sets[cs1]:
					for cid2 in chain_sets[cs2]:
						if report['chains_in_contact'][(cid1, cid2)]:
							contact_found = True
							break
					if contact_found:
						break
				if contact_found:
					break
			if not contact_found:
				print("DETACHED SET OF CHAINS", chain_sets[cs1])
				detached_units.append(chain_sets[cs1])
				error = True

	if detached_units:
		pass_entry = passport_entry(this_name+'_detached', pdbi, "There are detached replicas")
		for du in detached_units:
			det_err = ('NOTICE', this_name, "Chains {0} form a detached replica".format(", ".join([x for x in du])))
			print_log(det_err)#, thr_log_status=thr_log_status, log_filename=log_filename)
			error = True
	if not error:
		pass_entry = passport_entry(this_name+'_nostrerr', pdbi, "No structural errors detected in the generation of the full structure")
	return pass_entry, error


def distribute_items_in_fixed_length_list(n_items, list_length):
	return [((list_length-1)//n_items)+1]*(((list_length-1)%n_items)+1) + [list_length//n_items]*(n_items - list_length%n_items)


def run_MUSCLE(options, locations, str_data, str_set, only_analysis=False):
    """Runs MUSCLE over all pairs of chains in complexes contained in str_set set
    """

    def format_mtx(mtx_d, pdbi, antp_str):
        antp = antp_str.get_fdict_template()
        alpha = sorted(list({p[0] for p in mtx_d}))
        d = len(alpha)
        seqid_mtx = np.ones((d, d))
        for i1, pc1 in enumerate(alpha):
            for i2, pc2 in list(enumerate(alpha))[i1+1:]:
                if (pc1, pc2) not in mtx_d:
                    pdbi_c1 = pdbi + "_" + pc1
                    pdbi_c2 = pdbi + "_" + pc2
                    print_log((
                        "ERROR",
                        this_name,
                        ("Output number {0} of job {1} misses the comparison "
                            "between {2} and {3}")\
                            .format(index, batch_job_code, pdbi_c1, pdbi_c2)
                        ))
                seqid_mtx[i1, i2] = seqid_mtx[i2, i1] = mtx_d[(pc1, pc2)]
                antp['seqid_mtx'] = seqid_mtx.tolist()
                antp['seqid_legend'] = alpha
        return antp

    this_name = run_MUSCLE.__name__

    if not str_set:
        return str_data

    if options['PATHS'][('sing', 'singularity')] and options['PATHS'][('container', 'singularity_container')]:
        muscle_path = options['PATHS'][('sigmuscle', 'sig_muscle_path')]
    else:
        if options['RUN'][('hpc', 'run_on_hpc')]:
            muscle_path = options['PATHS'][('hpcmuscle', 'hpc_muscle_path')]
        else:
            muscle_path = options['PATHS'][('muscle', 'muscle_path')]

    for db in ['PDB', 'OPM']:
        # Instruct locusts
        batch_job_code = options['RUN'][('code', 'main_code')] + db + "_MUSCLE"
        output_dir = locations['FSYSPATH']['cache'] + 'output_{0}/'.format(batch_job_code)

        if not only_analysis:
            # .0 Create locusts parameter file
            parameter_file = locations['FSYSPATH']['logs'] + 'MUSCLE_locusts.par'
            write_locusts_parfile(options, parameter_file, batch_job_code)

            # .1 Create input and output dir
            input_dir = locations['FSYSPATH']['cache'] + 'input_{0}/'.format(batch_job_code)
            if os.path.exists(input_dir):
                shutil.rmtree(input_dir)
            os.mkdir(input_dir)
            if os.path.exists(output_dir):
                shutil.rmtree(output_dir)
            os.mkdir(output_dir)
            seqid_py = input_dir + "seqid.py"
            shutil.copyfile(locations['SYSFILES']['seqid_py'], seqid_py)
    
    
            # .2 Create input files:
            # .2.1-2 MUSCLE seqlist and command files
            command_file_content = (
                '# Argument1 : <id>\n'
                'SEQS_FILENAME="${{1}}_seqpairs.txt"\n'
                'TMPLIST_FILENAME="${{1}}.tmp.txt"\n'
                'LOG_FILENAME="${{1}}.log.txt"\n'
                'OUT_FILENAME="${{1}}_seqids.txt"\n'
                'rm ${{TMPLIST_FILENAME}}\n'
                "awk -v lfn=seqlist.tmp.txt '\n"
                '$1=="INIT"{{if (filename) {{close(filename)}}; filename = $2 "_" $3 ".fa"; printf "" > filename; next}}\n'
                "{{print >> filename; print filename >> lfn}}' ${{SEQS_FILENAME}}\n"
                "sort -nk1,2 seqlist.tmp.txt | uniq > ${{TMPLIST_FILENAME}}\n"
                "for x in `cat ${{TMPLIST_FILENAME}}`\n"
                "do\n"
                '\t{0} -in ${{x}} -out ${{x%".fa"}}.out.fa\n'
                '\tpython3 seqid.py ${{x%".fa"}}.out.fa >> ${{OUT_FILENAME}}\n'
                "done &> ${{LOG_FILENAME}}\n"
                "rm ${{TMPLIST_FILENAME}}\n"
            ).format(muscle_path)

            tmplist = []
            spf = ""
            for pdbi in sorted(str_set):
                if len(tmplist) == 10:
                    index = tmplist[0]+"-"+tmplist[-1]
                    tmplist = []
                    seq_pairs_filename = input_dir + "{0}_seqpairs.txt".format(index)
                    with open(seq_pairs_filename, 'w') as spff:
                        spff.write(spf)
                    command_filename = input_dir + "{0}_cmd.sh".format(index)
                    with open(command_filename, 'w') as cf:
                        cf.write(command_file_content)
                    subprocess.call(["chmod", "777", command_filename])
                    spf = ""
    
                db_cache = str_data[pdbi]['FROM_'+db]['cache_coordinate_file']
                if os.path.exists(db_cache):
                    db_coords = retrieve_coords_from_CIF(db_cache)
                    seqs = {}
                    for c in db_coords['COORDS']:
                        seqs[c] = [from3to1(x[1]) for x in db_coords['COORDS'][c]]
                    sn = pdbi + "_"
                    chlist = sorted(list(seqs.keys()))
                    nc = len(chlist)
                    seq_pairs_filename = input_dir + "{0}_seqpairs.txt".format(pdbi)
                    for ic1, c1 in enumerate(chlist):
                        sn1 = sn + c1
                        for ii, c2 in enumerate(chlist[ic1+1:]):
                            ic2 = ic1 + 1 + ii
                            sn2 = sn + c2
                            spf += "INIT\t{0}\t{1}\n".format(sn1, sn2)
                            spf += ">{0}\n{1}\n>{2}\n{3}\n".format(sn1, "".join(seqs[c1]), sn2, "".join(seqs[c2]))
                else:
                    print_log((
                        'ERROR',
                        this_name,
                        'Database cache file {0} of entry {1} does not exist'\
                            .format(db_cache, pdbi)
                    ))
                    print("ERROR", db_cache, pdbi)
                tmplist.append(pdbi)
            if tmplist:
                index = tmplist[0]+"-"+tmplist[-1]
                seq_pairs_filename = input_dir + "{0}_seqpairs.txt".format(index)
                with open(seq_pairs_filename, 'w') as spff:
                    spff.write(spf)
                command_filename = input_dir + "{0}_cmd.sh".format(index)
                with open(command_filename, 'w') as cf:
                    cf.write(command_file_content)
                subprocess.call(["chmod", "777", command_filename])

            # .3 Instruct locusts
            specific_inputs = ['<id>_seqpairs.txt', '<id>_cmd.sh']
            shared_inputs = ['pyscr:seqid.py']
            exec_filename = locations['FSYSPATH']['cache'] + db + '_MUSCLE_exec_instructions.txt'
            command_template = "cp <shared>pyscr . ; bash <id>_cmd.sh <id>"
            outputs = ['<id>_seqids.txt']
   
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

        # .5 Collect
        # Read output.log and get the directory of the output file for each index in str_set (3rd field)
        # For each of these output files read the content and make one square matrix and one index list for each pdbi
        # Then write these two data in the appropriate str_data space
        outlog = output_dir + 'output.log'
        fromdb = "FROM_{0}".format(db)
        strcount = set() 
        seqid_legend = [] # Chain names in order
        seqid_linmtx = [] # Linear seqid matrix

        # Main output "directory" file
        with open(outlog) as olf:
            for line in olf:
                if not line.strip():
                    continue
                fields = line.split()
                index = fields[0].split(".")[0]
                presence = fields[1]

                # Check whether all output parts are present
                if presence != "present":
                    print_log((
                        'CRITICAL',
                        this_name,
                        ("Output number {0} of job {1} is not present, "
                            "EncoMPASS halts").format(index, batch_job_code)
                        ))

                # Parses the results files
                address = fields[2]
                mtx_d = {}
                prev_pdbi = ""
                with open(address) as adf:
                    for adline in adf:
                        if not adline.strip():
                            continue
                        pdbi_c1, pdbi_c2, seqid = adline.split()
                        (pdbi, c1), (_, c2) =\
                            pdbi_c1.split("_"), pdbi_c2.split("_")
                        #print(pdbi_c1, pdbi_c2, seqid)

                        # Records NC*NC seqid matrix
                        if pdbi not in strcount:
                            if mtx_d: 
                                str_data[prev_pdbi][fromdb]['analysis'] =\
                                    FixedDict(format_mtx(mtx_d, prev_pdbi, str_data[pdbi][fromdb]['analysis']))
                                #str_data[prev_pdbi][fromdb]['analysis'].show_dictionary()
                            prev_pdbi = pdbi
                            mtx_d = {}
                            strcount.add(pdbi)

                        # Records present line
                        if (c1, c2) not in mtx_d:
                            mtx_d[(c1, c2)] = mtx_d[(c2, c1)] = seqid
                        else:
                            print_log((
                                'ERROR',
                                this_name,
                                ("In output file {0}, chain pair {1} "
                                    "is repeated").format(address,
                                    (c1, c2))
                                ))

                    # Records last NC*NC seqid matrix
                    if mtx_d:
                        str_data[pdbi][fromdb]['analysis'] =\
                            FixedDict(format_mtx(mtx_d, prev_pdbi, str_data[pdbi][fromdb]['analysis']))

    for pdbi in strcount:
        for db in ['PDB', 'OPM']:
            fromdb = 'FROM_'+db
            if not str_data[pdbi][fromdb]['analysis']['seqid_legend']:
                db_cache = str_data[pdbi][fromdb]['cache_coordinate_file']
                if os.path.exists(db_cache):
                    db_coords = retrieve_coords_from_CIF(db_cache)
                    chainlist = [c for c in db_coords['COORDS']]
                    if len(chainlist) > 1:
                        print("ERROR: missing chains in MUSCLE output", pdbi, chainlist, str_data[pdbi][fromdb]['analysis']['seqid_legend'])
                        exit(1)
                    else:
                        antp = str_data[pdbi][fromdb]['analysis'].get_fdict_template()
                        antp['seqid_mtx'] = [1.0]
                        antp['seqid_legend'] = chainlist[0]
                else:
                    if db.lower() + '_eliminated' in str_data[pdbi]['status']:
                        print("OK: missing coordinate file taken already into account", pdbi, db_cache)
                    else:
                        print("ERROR: missing coordinate file", pdbi, db_cache)
                        exit(1)

    return str_data


def run_PPM(options, locations, str_set, str_data, only_analysis=False):
    this_name = run_PPM.__name__

    LOC_pad = locations['SYSFILES']['pdbatomdict']

    if options['PATHS'][('sing', 'singularity')] and options['PATHS'][('container', 'singularity_container')]:
        ppm_path = options['PATHS'][('sigppm', 'sig_ppm_path')]
    else:
        if options['RUN'][('hpc', 'run_on_hpc')]:
            ppm_path = options['PATHS'][('hpcppm', 'hpc_ppm_path')]
        else:
            ppm_path = options['PATHS'][('ppm', 'ppm_path')]

    # Instruct locusts
    batch_job_code = options['RUN'][('code', 'main_code')] + "PPM"
    output_dir = locations['FSYSPATH']['cache'] + 'output_{0}/'.format(batch_job_code)

    if str_set and not only_analysis:
        # .0 Create locusts parameter file
        parameter_file = locations['FSYSPATH']['logs'] + 'PPM_locusts.par'
        write_locusts_parfile(options, parameter_file, options['RUN'][('code', 'main_code')] + '_PPM')

        # .1 Create input and output dir
        input_dir = locations['FSYSPATH']['cache'] + 'input_{0}/'.format(batch_job_code)
        if os.path.exists(input_dir):
            shutil.rmtree(input_dir)
        os.mkdir(input_dir)
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.mkdir(output_dir)
    
        # .2 Create input files:
        # .2.1 PPM input spec
        for pdbi in str_set:
            inpfname = input_dir + pdbi + '.inp'
            with open(inpfname, 'w') as inpf:
                inpf.write(" 0 in  {0}_temp.pdb\n".format(pdbi))
    
        # .2.2 temporary PDB file
        for pdbi in str_set:
            tmp_filename = input_dir + pdbi + '_temp.pdb'
            cif_filename = str_data[pdbi]['present_cache']
            if not cif_filename:
                msg = 'Structure {0} has no designated present cache file'.format(pdbi)
                print_log('CRITICAL', this_name, msg)
            c_dict_pdbi = retrieve_coords_from_CIF(cif_filename)
            write_ENC(str_data[pdbi], c_dict_pdbi, tmp_filename, locations['SYSFILES']['ignoredpasscodes'], LOC_pad)
    
        # .2.3 res.lib file must be present in each clean environment
        if options['PATHS'][('sing', 'singularity')] and options['PATHS'][('container', 'singularity_container')]:
            cmd = options['PATHS'][('sing', 'singularity')] + ' exec ' + options['PATHS'][('container', 'singularity_container')] + ' cp  /EncoMPASS/ppm/res.lib ' + input_dir
        else:
            cmd = 'cp ' + os.path.dirname(options['PATHS'][('', 'ppm_reslib_path')]) + '/res.lib ' + input_dir
        p = subprocess.Popen(cmd.split(), stdout=open('/dev/null', 'w'), stderr=open('/dev/null', 'w'))
        p.wait() 
    
        # .3 Instruct locusts
        specific_inputs = ['<id>_temp.pdb', '<id>.inp']
        shared_inputs = ['rl:res.lib']
        outputs = ['<id>_temp_datapar1', '<id>_temp_datasub1', '<id>_tempout.pdb']
        command_template = 'cp <shared>rl . ; {0} < {1} > {2}'.format(ppm_path, '<id>.inp', '<id>_log.out')
    
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
    outputs = ['<id>_temp_datapar1', '<id>_temp_datasub1', '<id>_tempout.pdb']

    # .5 Collect
    output_paths = {}
    if str_set:
        with open(output_dir + "output.log") as of:
            for line in of:
                fields = line.split()
                output_paths[fields[0]] = (fields[1], fields[2])

    # Analysis
    for pdbi in str_set:
        failed = False
        for outfname in outputs:
            if output_paths[outfname.replace('<id>', pdbi)][0] != "present":
                failed = True
                break
        if failed:
            str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_fppm', pdbi, "During the analysis with the PPM software, the process has failed because of unexpected structural or format inconsistencies. This entry is scheduled for deletion"))
            str_data[pdbi]['status'].append('eliminated')
            str_data[pdbi]['delete_keyword'] = 'Unclear'
            continue

        # Output PDB file
        ppmout_path = output_paths['<id>_tempout.pdb'.replace('<id>', pdbi)][1]
        PPM_dict, _ = parse_PDB_coords(ppmout_path, with_UNK=options['PARAMETERS'][('', 'with_UNK')])
        if PPM_dict and 'COORDS' in PPM_dict and PPM_dict['COORDS']:
            str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_ppmok', pdbi, "The software PPM has successfully completed the analysis of this structure."))
        else:
            str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_ppmerr', pdbi, "The software PPM could not complete the analysis for this structure. This entry is scheduled for deletion"))
            str_data[pdbi]['status'].append('eliminated')
            str_data[pdbi]['delete_keyword'] = 'Unclear'
            continue
        str_data[pdbi]['present_cache'] = locations['FSYSPATH']['TMPcoords'] + pdbi + '_fromPPM.cif'
        write_mmCIF(str_data[pdbi], PPM_dict, str_data[pdbi]['present_cache'], LOC_pad)
        write_ENC(str_data[pdbi], PPM_dict, str_data[pdbi]['ENCOMPASS']['coordinate_file'], locations['SYSFILES']['ignoredpasscodes'], LOC_pad)
        if not PPM_dict['DUMMIES']:
            print("ERROR: NO DUM ATOMS")
            exit(1)
        str_data[pdbi]['ENCOMPASS']['lim_sup'] = abs(PPM_dict['DUMMIES'][0][1][2])
        str_data[pdbi]['ENCOMPASS']['lim_inf'] = -str_data[pdbi]['ENCOMPASS']['lim_sup']

        # Output segments
        with open(output_paths['<id>_temp_datasub1'.replace('<id>', pdbi)][1]) as datasub_file:
            for line in datasub_file:
                if not line.strip():
                    continue
                subud = FixedDict(str_data[pdbi]['FROM_OPM']['subunits'].get_fdict_template())
                macrofields = line.split(';')
                protein_letter = macrofields[1]
                subud['tilt'] = macrofields[2]
                subud['segment'] = []
                for seg in macrofields[3].split('(')[1:]:
                    n = seg.split('-')
                    if len(n) == 2:
                        n1, n2 = seg.split('-')
                        n1 = int(''.join(n1.split()))
                        n2 = int(''.join(n2.split(')')[0].split()))
                    elif len(n) == 1:
                        n1, n2 = n[0][:n[0].index(')')].split()
                    else:
                        print_log('CRITICAL', this_name, 'Structure {0} has shown inconsistencies in the definition of TM domains'.format(pdbi))
                    subud['segment'].append(tuple([n1, n2]))
                str_data[pdbi]['FROM_OPM']['ksubunits'].append(protein_letter)
                str_data[pdbi]['FROM_OPM']['subunits'].append(subud)
            if not str_data[pdbi]['FROM_OPM']['ksubunits']:
                str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_notm', pdbi, "The software PPM has not found any transmembrane domain inside this structure. This entry is scheduled for deletion"))
                str_data[pdbi]['status'].append('eliminated')
                str_data[pdbi]['delete_keyword'] = 'Monotopic'
                continue

    return str_data


def make_structurewise_table(str_data, filename):
    df = [[
        "PDB", 
        "CHAIN", 
        "UNIPROT_ACC", 
        "CHAINS_IN_COMPLEX", 
        "TMCHAINS_IN_COMPLEX", 
        "CHAIN_LEN", 
        "TOPOL_ELMTS", 
        "TM_REGIONS", 
        "TM_COVERAGE", 
        "EXT_N_TERM", 
        "EXT_MIDDLE", 
        "EXT_C_TERM", 
        "SS_COV", 
        "SS_COV_IN_TM_REG"
    ]]
    dftypes = {
        "PDB" : "category",
        "CHAIN" : "category",
        "UNIPROT_ACC" : "category",
        "CHAINS_IN_COMPLEX" : "uint8",
        "TMCHAINS_IN_COMPLEX" : "uint8",
        "TOPOL_ELMTS" : "uint8",
        "TM_REGIONS" : "uint8",
        "TM_COVERAGE" : "uint16",
        "EXT_N_TERM" : "uint16",
        "EXT_MIDDLE" : "uint16",
        "EXT_C_TERM" : "uint16",
        "SS_COV" : "uint16",
        "SS_COV_IN_TM_REG" : "uint16"
    }

    for pdbi in str_data: 
        chs = len(str_data[pdbi]['ENCOMPASS']['structure']['kchains'])
        tmchs = len([x for x in str_data[pdbi]['ENCOMPASS']['structure']['ktmchains'] if x != '-'])
        print(df[0])
        for ich, ch in enumerate(str_data[pdbi]['ENCOMPASS']['structure']['kchains']):
            chfd = str_data[pdbi]['ENCOMPASS']['structure']['chains'][ich]
            if type(chfd['UniProt_acc']) == str:
                try:
                    unp_acc = ",".join([x[0] for x in ast.literal_eval(chfd['UniProt_acc'])])
                except:
                    unp_acc = chfd['UniProt_acc']
            elif type(chfd['UniProt_acc']) == list:
                unp_acc = ",".join([x[0] for x in chfd['UniProt_acc']])
            else:
                unp_acc = ""
            ch_len = len(chfd['sequence'])
            nTEs = len(chfd["TEs"])
            nTMs = len(chfd['TM_regions']['TM_regions_extrema'])
            ch_tm_cov = chfd['TM_regions']['TM_coverage'] if chfd['TM_regions']['TM_coverage'] else 0
            len_nterm = chfd['TM_regions']['Nterm_length'] if chfd['TM_regions']['Nterm_length'] and nTMs else 0
            len_middle = chfd['TM_regions']['middle_linkers_length'] if chfd['TM_regions']['middle_linkers_length'] and nTMs else 0
            len_cterm = chfd['TM_regions']['Cterm_length'] if chfd['TM_regions']['Cterm_length'] and nTMs else 0
            ss_res = chfd['ss_coverage'] if chfd['ss_coverage'] else 0
            ch_ss_cov_in_tm = chfd['TM_regions']['ss_coverage_of_TMareas'] if chfd['TM_regions']['ss_coverage_of_TMareas'] else 0
            print(pdbi, ch, unp_acc, chs, tmchs, ch_len, nTEs, nTMs, ch_tm_cov, len_nterm, len_middle, len_cterm, ss_res, ch_ss_cov_in_tm)
            df.append([pdbi, ch, unp_acc, chs, tmchs, ch_len, nTEs, nTMs, ch_tm_cov, len_nterm, len_middle, len_cterm, ss_res, ch_ss_cov_in_tm])
                
    structwise_table = pd.DataFrame(df[1:], columns=df[0]).astype(dftypes)
    structwise_table.to_csv(filename, sep="\t", index=False)  # DOES NOT WRITE THE INDEX


def post_insertion_analysis(options, locations, str_data):
    this_name = post_insertion_analysis.__name__

    for pdbi in str_data:
        if 'eliminated' in str_data[pdbi]['status']:
            continue

        tmchs = 0
        chs = 0
        # Topological elements
        topological_elements, stats = define_TM_regions(options, locations, str_data[pdbi]['ENCOMPASS']['coordinate_file'], str_data[pdbi], pdbi, c_dict={}, cache_coord_filename='', span=None)

        #print(pdbi, "TOPOL ELEM", topological_elements)
        #print(stats)

        ch_tm_cov = 0
        ch_ss_cov_in_tm = 0
        thereareTEs = False
        for ich, ch in enumerate(str_data[pdbi]['ENCOMPASS']['structure']['kchains']):
            ch_tm_cov = 0
            ch_ss_cov_in_tm = 0
            chs += 1
            tei = 0
            tmi = 0
            aftertm = 0
            chosen_st_tmrs = []
            ilinker = 0
            if ch not in topological_elements:
                str_data[pdbi]['ENCOMPASS']['structure']['ktmchains'].append('-')
                str_data[pdbi]['ENCOMPASS']['structure']['chains'][ich]['is_TM'] = False
                continue
            ch_len, ss_res, len_linkers = stats["STATS_CHAINS"][ch]
            len_tmlinkers = [len_linkers[0]]
            for te, st_tmr in zip(topological_elements[ch], stats["STATS_TM_REGIONS"][ch][1:]):
                te_weights, te_seg = te
                te_weight, tm_weight = te_weights
                first_tmres, last_tmres, len_tmreg, len_ssintm = st_tmr
                te_segtype, te_segid = te_seg[0], te_seg[1:]
                if te_weight < 0.5:
                    len_tmlinkers[-1] += len_tmreg + len_linkers[ilinker+1]
                    temsg = ('NOTICE', this_name, "Structure {0} chain {1} TE {2} is discarded: TE_score {3}, TE_segtype {4}, TE_residues {5}".format(pdbi, ch, tei, te_weight, te_segtype, [x[0][0] for x in te_segid]))
                    print_log(temsg)
                    continue

                tei += 1

                if tm_weight > 0:
                    tmi += 1
                    is_tm = True
                    anf = [x[0][0] for x in te_segid]
                    # First and last TM residues can be outside TMTE!
                    if first_tmres[0][0] in anf:
                        beforetm = len(anf[:anf.index(first_tmres[0][0])])
                    else:
                        beforetm = 0
                    if last_tmres[0][0] in anf:
                        aftertm = len(anf[anf.index(last_tmres[0][0]):])
                    else:
                        aftertm = 0
                    len_tmlinkers[-1] += beforetm
                    len_tmlinkers.append(len_linkers[ilinker+1] + aftertm)
                    ch_tm_cov += len_tmreg
                    ch_ss_cov_in_tm += len_ssintm
                    chosen_st_tmrs.append(st_tmr)
                else:
                    is_tm = False
                    len_tmlinkers[-1] += len_tmreg + len_linkers[ilinker+1]
                ilinker += 1

                temsg = ('NOTICE', this_name, "Structure {0} chain {1} TE {2}: is_TM {3}, TE_score {4}, TM_score {5}, TE_segtype {6}, TE_residues {7}".format(pdbi, ch, tei, is_tm, te_weight, tm_weight, te_segtype, [x[0][0] for x in te_segid]))
                print_log(temsg)

                tet = FixedDict(str_data[pdbi]['ENCOMPASS']['structure']['chains'][ich]['TEs'].get_fdict_template())
                tet['ss_type'] = te_segtype
                tet['residues_term'] = [te_segid[0], te_segid[-1]]
                tet['TE_score'] = te_weight
                tet['TM_score'] = tm_weight
                tet['is_TM'] = is_tm
                str_data[pdbi]['ENCOMPASS']['structure']['chains'][ich]['TEs'].append(tet)
                thereareTEs = True
            if tmi:
                tmchs += 1
                str_data[pdbi]['ENCOMPASS']['structure']['ktmchains'].append(ch)
                str_data[pdbi]['ENCOMPASS']['structure']['chains'][ich]['is_TM'] = True
            else:
                str_data[pdbi]['ENCOMPASS']['structure']['ktmchains'].append('-')
                str_data[pdbi]['ENCOMPASS']['structure']['chains'][ich]['is_TM'] = False

            len_nterm = len_tmlinkers[0]
            len_cterm = len_tmlinkers[-1] + aftertm
            len_middle = sum(len_tmlinkers[1:-1])

            str_data[pdbi]['ENCOMPASS']['structure']['chains'][ich]['ss_coverage'] = ss_res
            tmfd = str_data[pdbi]['ENCOMPASS']['structure']['chains'][ich]['TM_regions']
            tmfd["TM_coverage"] = ch_tm_cov
            tmfd["Nterm_length"] = len_nterm
            tmfd["middle_linkers_length"] = len_middle
            tmfd["Cterm_length"] = len_cterm
            tmfd["ss_coverage_of_TMareas"] = ch_ss_cov_in_tm
            tm_extrema = []
            for st_tmr in chosen_st_tmrs:
                first_tmres, last_tmres, len_tmreg, len_ssintm = st_tmr
                tm_extrema.append((first_tmres, last_tmres))
            tmfd["TM_regions_extrema"] = tm_extrema

            str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_tei', pdbi, "Number of Topological Elements found in chain {0}: {1}".format(ch, tei)))
            str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_tmi', pdbi, "Number of TM Topological Elements found in chain {0}: {1}".format(ch, tmi)))
        if not thereareTEs:
            str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_note', pdbi, "No Topological Elements have been found in the transmembrane subunits of this structure. Its insertion in the lipid bilayer has not been confirmed : this entry is scheduled for deletion"))
            str_data[pdbi]['status'].append('eliminated')
            str_data[pdbi]['delete_keyword'] = 'Monotopic'

    for pdbi in str_data:
        pdbi_no_TM = True
        for ich, ch in [(i, x) for i, x in enumerate(str_data[pdbi]['ENCOMPASS']['structure']['ktmchains']) if x != '-']:
            N_TM_regions = len(str_data[pdbi]['ENCOMPASS']['structure']['chains'][ich]['TM_regions']['TM_regions_extrema'])
            if not N_TM_regions:
                str_data[pdbi]['ENCOMPASS']['structure']['ktmchains'][ich] = '-'    
                str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_0tmr', pdbi, "No TM regions in chain {0}".format(ch)))
            else:
                pdbi_no_TM = False
        if pdbi_no_TM:
            str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_0tmr', pdbi, "This entry has no TM regions"))
            if 'eliminated' not in str_data[pdbi]['status']:
                str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_0tmr', pdbi, "This entry has no TM regions"))
                str_data[pdbi]['status'].append('eliminated')
                str_data[pdbi]['delete_keyword'] = 'Monotopic'

    for pdbi in str_data:
        if str_data[pdbi]['ENCOMPASS']['structure']['ktmchains']:
            str_data[pdbi]['ENCOMPASS']['n_of_chains'] = len(str_data[pdbi]['ENCOMPASS']['structure']['ktmchains'])
            str_data[pdbi]['ENCOMPASS']['n_of_tmchains'] = len([x for x in str_data[pdbi]['ENCOMPASS']['structure']['ktmchains'] if x != '-'])
        

    make_structurewise_table(str_data, locations['SYSFILES']['structurewise_table'])

    write_checkpoint(str_data, locations['FSYSPATH']['cache'] + 'str_data_5.pkl')

    return str_data


def orthogonal_distance_regression(data):
	# The ODR line will pass from the centroid (can be proved)
	datamean = data.mean(axis=0)

	# Do an SVD on the mean-centered data. vv[0] will contain the direction of the best fit line
	uu, dd, vv = np.linalg.svd(data - datamean)

	return datamean, vv[0]


def calculate_trace(dseg, ch_coords):
	"""
	Careful: this function contains implicit thresholds!
	No trace will be output if the segment is not at least as long as the window for that type of secondary structure
	"""
	alpha_types = ['G', 'H', 'I']
	alpha_intervals = [3, 4, 5]	# These are thresholds
	beta_types = ['B', 'E']
	beta_intervals = [3, 3]		# These are thresholds
	ss_type = dseg[0]
	trace = []
	if len(dseg[1:]) < 3:
		return trace
	if ss_type in alpha_types:
		alpha_thr = alpha_intervals[alpha_types.index(ss_type)]
		for nr in range(len(dseg[1:-alpha_thr])):
			points = []
			vectors = []
			for nnr in range(nr, nr+alpha_thr):
				resid_name = dseg[1+nnr]
				v1, v2 = None, None
				if resid_name in ch_coords:
					for atname in ch_coords[resid_name]:
						if atname == 'O':
							v2 = np.array(ch_coords[resid_name][atname])
						else:
							points.append(ch_coords[resid_name][atname])
							if atname == 'C':
								v1 = np.array(ch_coords[resid_name][atname])
				
					if type(v1) == type(None) or type(v2) == type(None):
						print("warning: no C or O in residue", resid_name)
					else:
						vectors.append(list(v2 - v1))
			points = np.array(points)
			centroid = points.mean(axis=0)
			vectors = np.array(vectors)
			meanCOdir = vectors.mean(axis=0)
			meanCOdir /= np.linalg.norm(meanCOdir)
			trace.append((centroid, meanCOdir))
	elif ss_type in beta_types:
		beta_thr = beta_intervals[beta_types.index(ss_type)]
		for nr in range(len(dseg[1:-beta_thr])):
			points = []
			for nnr in range(nr, nr+beta_thr):
				resid_name = dseg[1+nnr]
				for atname in ch_coords[resid_name]:
					points.append(ch_coords[resid_name][atname])
			points = np.array(points)
			centroid, line_vec = orthogonal_distance_regression(points)
			trace.append((centroid, line_vec))
	return trace


def define_TM_regions(options, locations, coord_filename, str_data_pdbi, pdbi, c_dict={}, cache_coord_filename='', span=None):
    this_name = define_TM_regions.__name__

    # Obtain DSSP segments
    dssp_out_filename = locations['FSYSPATH']['cache'] + pdbi + '_dssp.txt'
    tmp_coord_filename = locations['FSYSPATH']['cache'] + pdbi + '_dssp.pdb'
    subprocess.Popen(['grep', '-v', 'DUM', coord_filename], stdout=open(tmp_coord_filename, 'w'))
    dssp_segments, Tres = make_DSSP(options['PATHS'][('dssp', 'dssp_path')], tmp_coord_filename, dssp_out_filename, pdbi, leave_original_dssp=True)
    dssp_to_ab = {'H':'A', 'G':'A', 'I':'A', 'B':'B', 'E':'B'}

    if not dssp_segments:
        return {}, {}

    # Coordinates are needed. If not provided, try to retrieve them
    if not c_dict:
        if cache_coord_filename:
            c_dict = retrieve_coords_from_CIF(cache_coord_filename)
        else:
            c_dict, _ = parse_PDB_coords(coord_filename, with_UNK=options['PARAMETERS'][('', 'with_UNK')])

    # Dummy atoms are needed.  
    if not c_dict['DUMMIES']:
        msg = passport_entry(this_name+'_tmnodum', pdbi, "The reference coordinate file has no DUM atoms. The TM regions cannot be defined without this information")
        return {}, {}

    lim_sup = abs(c_dict['DUMMIES'][0][1][2])
    lim_inf = -lim_sup

    stats_tmrs = {}
    stats_tmchs = {}

    # Calculate attributes of relevant DSSP segments
    topological_elements = {}
    for ch in dssp_segments:
        stats_tmrs[ch] = [['FIRST_TMRES', 'LAST_TMRES', 'LEN_TMREG', 'LEN_SSINTM']]
        dssp_segments_refined = []
        traces = []
        ss_res = 0
        for dseg in dssp_segments[ch]:
            # Calculate trace
            # CAUTION: contains implicit ss-type-dependent thresholds!!
            # G = 3; H = 4; I = 5; B = 3; E = 3.
            ss_res += len(dseg) - 1
            traces.append(calculate_trace(dseg, c_dict['COORDS'][ch]))

        jump_next = False
        for ntrace, trace in enumerate(traces[:-1]):
            dseg1 = dssp_segments[ch][ntrace][:]
            dseg2 = dssp_segments[ch][ntrace+1][:]
            tr1 = traces[ntrace]
            tr2 = traces[ntrace+1]
            if jump_next:
                dseg1 = dssp_segments_refined[-1][:]
                dseg2 = dssp_segments[ch][ntrace+1][:]
                tr1 = calculate_trace(dseg1, c_dict['COORDS'][ch])
                tr2 = traces[ntrace+1]
            # Merging segments
            # First, extend segments with Turns T
            dsegp1, dsegp2 = dseg1[:], dseg2[:]
            if dssp_to_ab[dsegp1[0]] == dssp_to_ab[dsegp2[0]] == "A":
                for ir in range(dseg1[-1][0][0]+1, dseg2[1][0][0]):
                    missed = True
                    if ch in Tres:
                        for x in Tres[ch]:
                            if x[0][0] == ir:
                                dsegp1.append(x)
                                missed = False
                                break
                    if missed:
                        break
                for ir in range(dseg2[1][0][0]-1, dseg1[-1][0][0], -1):
                    missed = True
                    if ch in Tres:
                        for x in Tres[ch]:
                           if x[0][0] == ir:
                               dsegp2 = [dsegp2[0]] + [x] + dsegp2[1:]
                               missed = False
                               break
                    if missed:
                        break

            linker = []
            for r in c_dict['COORDS'][ch]:
                if r[0][0] > dsegp1[-1][0][0] and r[0][0] < dsegp2[1][0][0]:
                    linker.append(r)
            linker = sorted(linker, key= lambda x: x[0][0])

            # Segments must be of the same DSSP type and closer than 4 residues, or of the same alpha/beta type, and closer than 3 residues
            if (dsegp1[0] == dsegp2[0] and dsegp2[1][0][0] - dsegp1[-1][0][0] < 4) or (dssp_to_ab[dsegp1[0]] == dssp_to_ab[dsegp2[0]] and dsegp2[1][0][0] - dsegp1[-1][0][0] <= 2):
                short_segments = False
                # Check if at least one of the two is short (less than 3 (G, B, E) or 4 (H) or 5 (I))
                # If not, calculate angle between traces
                if (len(tr1) == 0) or (len(tr2) == 0):
                    short_segments = True
                else:
                    ang = math.acos(np.dot(tr1[-1][1], tr2[0][1]))/math.pi*180
                # If segments are short and close or they are not too angled, merge them
                to_be_merged = False
                if (short_segments and dsegp2[1][0][0] - dsegp1[-1][0][0] <= 2) or (not short_segments and abs(ang) < 90):
                    to_be_merged = True
                elif not short_segments:
                    dist_from_0_sign_1s = (lim_inf+(lim_sup-lim_inf)/2 - c_dict['COORDS'][ch][dsegp1[1]]['CA'][2])/((lim_sup-lim_inf)/2)
                    dist_from_0_sign_2e = (lim_inf+(lim_sup-lim_inf)/2 - c_dict['COORDS'][ch][dsegp2[-1]]['CA'][2])/((lim_sup-lim_inf)/2)
                    if abs(dist_from_0_sign_1s) >= 1 and abs(dist_from_0_sign_2e) >= 1:
                        maxdepth = 0
                        for r in linker:
                            dist_from_0_sign_r = (lim_inf+(lim_sup-lim_inf)/2 - c_dict['COORDS'][ch][r]['CA'][2])/((lim_sup-lim_inf)/2)
                            initmem = math.copysign(1, dist_from_0_sign_1s)
                            depth = initmem*(initmem - dist_from_0_sign_r)
                            if maxdepth < depth:
                                maxdepth = depth
                        if maxdepth < 5/3:
                            to_be_merged = True
                if to_be_merged:
                    if jump_next:
                        added_seg = []
                        for ix in range(dssp_segments_refined[-1][-1][0][0]+1, dsegp2[1][0][0]):
                            for resid_name in c_dict['COORDS'][ch]:
                                resid, resname = resid_name
                                if ix == resid[0] and resid not in [dssp_segments_refined[-1][-1][0], dsegp2[1][0]]:
                                    added_seg.append((resid, resname))
                        dssp_segments_refined[-1] += added_seg + dsegp2[1:]
                    else:
                        added_seg = []
                        for ix in range(dsegp1[-1][0][0]+1, dsegp2[1][0][0]):
                            for resid_name in c_dict['COORDS'][ch]:
                                resid, resname = resid_name
                                if ix == resid[0] and resid not in [dsegp1[-1][0], dsegp2[1][0]]:
                                    added_seg.append((resid, resname))
                        dssp_segments_refined.append(dsegp1 + added_seg + dsegp2[1:])
                    jump_next = True
                    continue
            if not jump_next:
                dssp_segments_refined.append(dseg1)
            jump_next = False
        if not jump_next:
            if dssp_segments[ch]:
                dssp_segments_refined.append(dssp_segments[ch][-1])
            else:
                continue
				
        linkers = []
        new_dssp_segments_refined = []
        if dssp_segments_refined:
            nterminal = []
            for r in c_dict['COORDS'][ch]:
                if r[0][0] < dssp_segments_refined[0][1][0][0]:
                    nterminal.append(r)
            nterminal = sorted(nterminal, key= lambda x: x[0][0])
            linkers.append(nterminal)
        else:
            seq = []
            for r in c_dict['COORDS'][ch]:
                seq.append(r)
            linkers.append(seq)
        for ndseg, dseg in enumerate(dssp_segments_refined[:-1]):
            # Calculate trace
            # CAUTION: contains implicit ss-type-dependent thresholds!!
            # G = 3; H = 4; I = 5; B = 3; E = 3.
            trace = calculate_trace(dseg, c_dict['COORDS'][ch])
            if not trace:
                continue
            linker = []
            for r in c_dict['COORDS'][ch]:
                if r[0][0] > dseg[-1][0][0] and r[0][0] < dssp_segments_refined[ndseg+1][1][0][0]:
                    linker.append(r)
            linker = sorted(linker, key= lambda x: x[0][0])
            linkers.append(linker)
            new_dssp_segments_refined.append(dseg)
        if dssp_segments_refined:
            new_dssp_segments_refined.append(dssp_segments_refined[-1])
            cterminal = []
            for r in c_dict['COORDS'][ch]:
                if r[0][0] > dssp_segments_refined[-1][-1][0][0]:
                    cterminal.append(r)
            cterminal = sorted(cterminal, key= lambda x: x[0][0])
            linkers.append(cterminal)
        else:
            continue

        linkers = [[]]
        ids_old = 0
        for r in c_dict['COORDS'][ch]:
            inds = False
            for ids, ds in enumerate(new_dssp_segments_refined): 
                if r in ds:
                    inds = True
                    if ids_old != ids:
                        linkers.append([])
                        ids_old = ids
                    break
            if inds:
                if not linkers[-1]:
                    continue
                linkers.append([])
            else:
                linkers[-1].append(r)
                

        dssp_segments[ch] = new_dssp_segments_refined[:]

        # Once the TEs are defined, we have to characterize them
        topological_elements[ch] = []
        len_linkers = [0]
        last_tmres = ((-999, ''), '')
        for ndseg, dseg in enumerate(dssp_segments[ch]):

            attributes = []
            # Calculate TE weight:
            # It is 1 if element enters the membrane
            # It is 0 < TEweight < 1 if element is close to the membrane
            # It is 0 otherwise
            TE_weight = 0
            TM_weight = 0
            TE_norm = 0
            TM_norm = 0
            initmem = math.copysign(1,c_dict['COORDS'][ch][dseg[1]]['CA'][2] - (lim_inf+(lim_sup-lim_inf)/2))
            intm_res = []
            beflink, aftlink = False, False
            for rn in range(len(dseg[1:])):
                if dseg[1+rn] not in c_dict['COORDS'][ch]:   # Why this?
                    continue
                # Min distance from center of membrane
                dist_from_0_sign = (c_dict['COORDS'][ch][dseg[1+rn]]['CA'][2] - (lim_sup+lim_inf)/2)/((lim_sup-lim_inf)/2)
                dist_from_0 = abs(dist_from_0_sign)
                if dist_from_0 < 1:
                    TE_weight += 1
                    # TM criterion: if depth is less than 25% of the membrane, it is a penalty
                    # otherwise a positive factor
                    TM_norm += 1
                    intm_res.append(dseg[1+rn])
                    if rn == 0:
                        beflink = True
                    elif rn == len(dseg[1:])-1:
                        aftlink = True
                    depth = initmem*(initmem - dist_from_0_sign)
                    TM_weight += (depth - 0.5)
                elif dist_from_0 < 1.4:
                    TE_weight += 1 - (dist_from_0 - 1)/0.4

                if dist_from_0 < 2:
                    TE_norm += 1
            len_ssintm = len(intm_res)

            # Check linker before TE
            subtrtolink = 0
            if beflink:
                for r in reversed(linkers[ndseg]):
                    dist_from_0_sign = (c_dict['COORDS'][ch][r]['CA'][2] - (lim_sup+lim_inf)/2)/((lim_sup-lim_inf)/2)
                    if abs(dist_from_0_sign) > 1 or r == last_tmres:
                        break
                    else:
                        depth = initmem*(initmem - dist_from_0_sign)
                        TM_weight += (depth - 0.5)
                        TM_norm += 1
                        intm_res.append(r)
                        subtrtolink += 1
            len_linkers[-1] += len(linkers[ndseg]) - subtrtolink

            # Check linker after TE
            subtrtolink = 0
            if aftlink:
                for r in linkers[ndseg+1]: # There will always be one more linker: the (maybe empty) cterminal. If it is not there, the sequence does not have segs
                    dist_from_0_sign = (c_dict['COORDS'][ch][r]['CA'][2] - (lim_sup+lim_inf)/2)/((lim_sup-lim_inf)/2)
                    if abs(dist_from_0_sign) > 1: # The first residue that goes out of the membrane makes the count stop
                        break
                    else:
                        depth = initmem*(initmem - dist_from_0_sign)
                        TM_weight += (depth - 0.5)
                        TM_norm += 1
                        intm_res.append(r)
                        subtrtolink += 1
            len_linkers.append(len(linkers[ndseg+1]) - subtrtolink)

            if TE_norm != 0:
                if TM_norm != 0:
                    TM_weight /= TM_norm
                TE_weight /= TE_norm

            intm_res = sorted(intm_res)
            if intm_res:
                first_tmres, last_tmres, len_tmreg = intm_res[0], intm_res[-1], len(intm_res)
            else:
                first_tmres, last_tmres, len_tmreg = 0, 0, 0
            stats_tmrs[ch].append([first_tmres, last_tmres, len_tmreg, len_ssintm])
            topological_elements[ch].append(((TE_weight, TM_weight), [dssp_to_ab[dseg[0]]] + dseg[1:]))
        stats_tmchs[ch] = [len(c_dict['COORDS'][ch]), ss_res, len_linkers]

    return topological_elements, {"STATS_CHAINS" : stats_tmchs, "STATS_TM_REGIONS" : stats_tmrs}


def division_into_chains(options, locations, pdbi, str_data_pdbi, str_filename, out_dir, is_ciff=False):
	this_name = division_into_chains.__name__

	LOC_pad = locations['SYSFILES']['pdbatomdict']

	if is_ciff:
		c_dict = retrieve_coords_from_CIF(str_filename)
	else:
		c_dict, _ = parse_PDB_coords(str_filename, with_UNK=options['PARAMETERS'][('', 'with_UNK')])


	if not (c_dict['COORDS'] and c_dict['DUMMIES']):
		det_err = ('NOTICE', this_name, "Structure {0} has a cache CIF file which is not up to date (does not contain COORD or DUMMIES) and for this reason it will not appear in EncoMPASS (TO BE CORRECTED).".format(pdbi))
		return []

	z = abs(c_dict['DUMMIES'][0][1][2])

	ch_list = []
	for ch in c_dict['COORDS']:
		ch_list.append(ch)
		chain_filename = out_dir + pdbi + '_' + ch + '_enc.pdb'
		chain_dict = add_DUM_to_coord_dict({'COORDS' : {ch : c_dict['COORDS'][ch]}, 'DUMMIES' : []}, z)
		write_ENC(str_data_pdbi, chain_dict, chain_filename, locations['SYSFILES']['ignoredpasscodes'], LOC_pad, no_pass=True)
	return ch_list


def structure_sorter(options, locations, str_data):
	def fill_str_data(coords, str_data_pdbi):
		n_of_aas = 0
		for ch in coords:
			chfd = FixedDict(str_data_pdbi['ENCOMPASS']['structure']['chains'].get_fdict_template())
			str_data_pdbi['ENCOMPASS']['structure']['kchains'].append(ch)
			chfd['residues'] = list(coords[ch])
			chfd['sequence'] = "".join([from3to1(x[1]) for x in coords[ch]])
			chfd['n_of_residues'] = len(coords[ch])
			n_of_aas += len(coords[ch])
			str_data_pdbi['ENCOMPASS']['structure']['chains'].append(chfd)
		str_data_pdbi['ENCOMPASS']['n_of_aas'] = n_of_aas
		str_data_pdbi['ENCOMPASS']['n_of_chains'] = len(coords)
		return str_data_pdbi

	LOC_pad = locations['SYSFILES']['pdbatomdict']

	this_name = structure_sorter.__name__
	# Division into chains
	for pdbi in str_data:
		if 'eliminated' in str_data[pdbi]['status']:
			write_ENC(str_data[pdbi], {}, str_data[pdbi]['ENCOMPASS']['coordinate_file'], locations['SYSFILES']['ignoredpasscodes'], LOC_pad, only_pass=True)
			continue
		whole_filename = locations['FSYSPATH']['whole'] + pdbi + '_enc.pdb'
		if str_data[pdbi]['present_cache'] and os.path.exists(str_data[pdbi]['present_cache']):
			c_dict = retrieve_coords_from_CIF(str_data[pdbi]['present_cache'])
			if c_dict['COORDS'] and c_dict['DUMMIES'] and ('eliminated' not in str_data[pdbi]['status']):
				ch_list = division_into_chains(options, locations, pdbi, str_data[pdbi], str_data[pdbi]['present_cache'], locations['FSYSPATH']['chains'], is_ciff=True)
				write_ENC(str_data[pdbi], c_dict, str_data[pdbi]['ENCOMPASS']['coordinate_file'], locations['SYSFILES']['ignoredpasscodes'], LOC_pad)
				write_ENC(str_data[pdbi], c_dict, whole_filename, locations['SYSFILES']['ignoredpasscodes'], LOC_pad)
				str_data[pdbi]['ENCOMPASS']['coordinate_file'] = whole_filename
				str_data[pdbi] = fill_str_data(c_dict['COORDS'], str_data[pdbi])
			else:
				str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_unkerr1', pdbi, "This entry encountered some critical error during the structural analyses. It will not be reflected in the EncoMPASS database."))
				str_data[pdbi]['status'].append('eliminated')
				str_data[pdbi]['delete_keyword'] = 'Unclear'
				write_ENC(str_data[pdbi], {}, whole_filename, locations['SYSFILES']['ignoredpasscodes'], LOC_pad, only_pass=True)
		elif str_data[pdbi]['ENCOMPASS']['coordinate_file'] and os.path.exists(str_data[pdbi]['ENCOMPASS']['coordinate_file']):
			ch_list = division_into_chains(options, locations, pdbi, str_data[pdbi], str_data[pdbi]['ENCOMPASS']['coordinate_file'], locations['FSYSPATH']['chains'])
			if not ch_list:
				str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_unkerr2', pdbi, "This entry encountered some critical error during the structural analyses. It will not be reflected in the EncoMPASS database."))
				str_data[pdbi]['status'].append('eliminated')
				str_data[pdbi]['delete_keyword'] = 'Unclear'
				write_ENC(str_data[pdbi], {}, str_data[pdbi]['ENCOMPASS']['coordinate_file'], locations['SYSFILES']['ignoredpasscodes'], LOC_pad, only_pass=True)
				continue
			ENC_dict, _ = parse_PDB_coords(output_dir + '<index>_tempout.pdb'.replace('<index>', pdbi), with_UNK=options['PARAMETERS'][('', 'with_UNK')])
			write_ENC(str_data[pdbi], ENC_dict, whole_filename, locations['SYSFILES']['ignoredpasscodes'], LOC_pad)
			shutil.copyfile(str_data[pdbi]['ENCOMPASS']['coordinate_file'], whole_filename)
			str_data[pdbi]['ENCOMPASS']['coordinate_file'] = whole_filename
			str_data[pdbi] = fill_str_data(c_dict['COORDS'], str_data[pdbi])
		else:
			str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_unkerr3', pdbi, "This entry encountered some critical error during the structural analyses. It will not be reflected in the EncoMPASS database."))
			str_data[pdbi]['status'].append('eliminated')
			str_data[pdbi]['delete_keyword'] = 'Unclear'
			write_ENC(str_data[pdbi], {}, str_data[pdbi]['ENCOMPASS']['coordinate_file'], locations['SYSFILES']['ignoredpasscodes'], LOC_pad, only_pass=True)

	write_checkpoint(str_data, locations['FSYSPATH']['cache'] + 'str_data_4.pkl')
	return str_data
			

def static_figures(options, locations, str_data):
	for pdbi in str_data: 
		pdb_filename = locations['FSYSPATH']['ENCpdbs'] + pdbi + '_enc.pdb'
		if not os.path.exists(pdb_filename):
			print("OHIBO", pdb_filename)
			continue
		pml_filename = locations['FSYSPATH']['stfig'] + pdbi + '.pml'
		fig_filename = locations['FSYSPATH']['stfig'] + pdbi + '.png'
		with open(pml_filename, 'w') as pml_file:
			pml_file.write(create_whole_png_text(pdb_filename, fig_filename))
		os.system("cd {0}; {1} -ucq {2}; cd -".format(locations['FSYSPATH']['stfig'], options['PATHS'][('pymol', 'pymol_path')], pml_filename))

def find_redundant_chains(locations, str_data):
    copies = {}
    epsilon = 0.5
    for pdbi in str_data:
        if 'eliminated' in str_data[pdbi]['status']:
            continue
        samech = []
        partiallists = []
        insame = set()
        allch = str_data[pdbi]['ENCOMPASS']['structure']['kchains']
        str_data[pdbi]['ENCOMPASS']['structure']['redundant_chains'] = {}
        for ich1, ch1 in enumerate(allch):
            print(pdbi, ch1)
            seq1 = str_data[pdbi]['ENCOMPASS']['structure']['chains'][ich1]['sequence']
            chain_pdbfile1 = locations['FSYSPATH']['chains'] + pdbi + '_' + ch1 + '_enc.pdb'
            for ich2, ch2 in [(ix, x) for ix, x in enumerate(allch)][ich1+1:]:
                seq2 = str_data[pdbi]['ENCOMPASS']['structure']['chains'][ich2]['sequence']
                if seq1 == seq2:
                    print(ch2, "seq1 == seq2")
                    print(ch1, ich1, seq1)
                    print(ch2, ich2, seq2)
                    chain_pdbfile2 = locations['FSYSPATH']['chains'] + pdbi + '_' + ch2 + '_enc.pdb'
                    ref_str = PDBParser(PERMISSIVE=True, QUIET=True, get_header=False).get_structure(pdbi, chain_pdbfile1)
                    alt_str = PDBParser(PERMISSIVE=True, QUIET=True, get_header=False).get_structure(pdbi, chain_pdbfile2)
                    rmsd = fast_superimpose(ref_str[0][ch1], alt_str[0][ch2])
                    print(ch2, "rmsd", rmsd, "eps", epsilon)
                    print(ch1, ch2, rmsd)
                    if rmsd < epsilon:
                        samech.append((ch1, ch2))
                        insame.add(ch1)
                        insame.add(ch2)
                        if ch2 not in str_data[pdbi]['ENCOMPASS']['structure']['redundant_chains']:
                            if ch1 in str_data[pdbi]['ENCOMPASS']['structure']['redundant_chains']:
                                true_ref = str_data[pdbi]['ENCOMPASS']['structure']['redundant_chains'][ch1]
                                print(pdbi, ch2, "pointing to a duplicate", ch1, "redirected to true ref", true_ref)
                            else:
                                true_ref = ch1
                            str_data[pdbi]['ENCOMPASS']['structure']['redundant_chains'][ch2] = true_ref
                            print(str_data[pdbi]['ENCOMPASS']['structure']['redundant_chains'])
            partiallists = sorted(merge(samech, sorted_lists=True))
        copies[pdbi] = sorted(merge(samech, sorted_lists=True)) + [[x] for x in allch if x not in insame]
        print("COPIES", pdbi)
        print(copies[pdbi])
    return str_data, copies


def fast_superimpose(ref_chain, alt_chain):
    # Checks if two chains are in the exact same conformation
    ref_atoms = []
    alt_atoms = []
    for ref_res, alt_res in zip(ref_chain, alt_chain):
        try:
            ref_atoms.append(ref_res['CA'])                
            alt_atoms.append(alt_res['CA'])
        except:
            continue
    super_imposer = Superimposer()
    super_imposer.set_atoms(ref_atoms, alt_atoms)
    super_imposer.apply(alt_chain.get_atoms())
    return super_imposer.rms


def generate_chain_pdb_files(options, locations, str_data, do_first_parse=True, thr_log_status="ERROR"):
	this_name = generate_chain_pdb_files.__name__

	str_data, run_in_PPM = parallel_generate_full_structure(options, locations, str_data)
	write_checkpoint(str_data, locations['FSYSPATH']['cache'] + 'str_data_2b.pkl')

	# PPM runner
	str_data = read_checkpoint(locations['FSYSPATH']['cache'] + 'str_data_2b.pkl', locations['SYSFILES']['data_structure_template'])
	run_in_PPM = set()
	for pdbi in str_data:
		if 'run_in_PPM' in str_data[pdbi]['status']:
			run_in_PPM.add(pdbi)
	#print("RUN IN PPM", run_in_PPM)
	str_data = run_PPM(options, locations, run_in_PPM, str_data, only_analysis=False)

	return str_data


def run_PPM_wrapper(options, locations, str_data, only_analysis=False):
    # Compile list of structures to be run
    run_in_PPM = set()
    for pdbi in str_data:
        if 'run_in_PPM' in str_data[pdbi]['status']:
            run_in_PPM.add(pdbi)
    print_log((
        "NOTICE", 
        run_PPM_wrapper.__name__,
        f"Structures to be run in PPM: {run_in_PPM}"
    ))

    # Run PPM
    str_data = run_PPM(options, locations, run_in_PPM, str_data, only_analysis=only_analysis)
    write_checkpoint(str_data, locations['FSYSPATH']['cache'] + 'str_data_3.pkl')
    return str_data


if __name__ == "__main__":
    from supporting_functions import *
    from initialize_repository import *
    from combine_sources import *

    options, locations = initialize_repository()
    ids, secd = mainlist_from_OPM(options, locations)
    str_data = scan_OPM(options, locations, ids)
    str_data = scan_PDB(options, locations, str_data=str_data)
    str_data = combine_sources(options, locations, str_data)
    #str_data = read_checkpoint(locations['FSYSPATH']['cache'] + 'str_data_1.pkl', locations['SYSFILES']['data_structure_template'])
    str_data, run_in_PPM = parallel_generate_full_structure(options, locations, str_data, cycle_0=True)
    #write_checkpoint(str_data, locations['FSYSPATH']['cache'] + 'str_data_2b.pkl')
    #str_data = read_checkpoint(locations['FSYSPATH']['cache'] + 'str_data_2b.pkl', locations['SYSFILES']['data_structure_template'])
    str_data = run_PPM_wrapper(options, locations, str_data, only_analysis=False)
    #str_data = read_checkpoint(locations['FSYSPATH']['cache'] + 'str_data_3.pkl', locations['SYSFILES']['data_structure_template'])
    str_data = structure_sorter(options, locations, str_data)
    #str_data = read_checkpoint(locations['FSYSPATH']['cache'] + 'str_data_4.pkl', locations['SYSFILES']['data_structure_template'])
    str_data = post_insertion_analysis(options, locations, str_data)
    str_data, unique_chains = find_redundant_chains(locations, str_data)
    str_data = assign_uniprot_acc(locations, str_data, unique_chains)
    
    LOC_pad = locations['SYSFILES']['pdbatomdict']
    for pdbi in str_data:
        if 'eliminated' in str_data[pdbi]['status']:
            only_pass = True
            c_dict_pdbi = {'COORDS' : {}, 'DUMMIES' : {}}
        else:
            only_pass = False
            c_dict_pdbi = retrieve_coords_from_CIF(str_data[pdbi]['present_cache'])
        write_ENC(str_data[pdbi], c_dict_pdbi, locations['FSYSPATH']['ENCpdbs'] + pdbi + '_enc.pdb', locations['SYSFILES']['ignoredpasscodes'], LOC_pad, only_pass=only_pass)
        write_ENC(str_data[pdbi], c_dict_pdbi, str_data[pdbi]['ENCOMPASS']['coordinate_file'], locations['SYSFILES']['ignoredpasscodes'], LOC_pad, only_pass=only_pass)
    write_checkpoint(str_data, locations['FSYSPATH']['cache'] + 'str_data_completegen.pkl')

    make_structurewise_table(str_data, locations['SYSFILES']['structurewise_table'])
    str_data = read_checkpoint(locations['FSYSPATH']['cache'] + 'str_data_completegen.pkl', locations['SYSFILES']['data_structure_template'])
