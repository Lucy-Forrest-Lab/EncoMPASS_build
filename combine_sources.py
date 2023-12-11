# Name: combine_sources.py
# Language: python3
# Author: Edoardo Sarti
# Date: Apr 03 2022

from supporting_functions import *

def pdb_query(path="pdb_query.json"):
    """Returns the json query for the PDB API (2021-)
    """
    with open(path) as pf:
        return pf.read()

def check_update_repo(options, locations, record=False, 
        thr_log_status='ERROR'):
    """Determines whether the repository must be updated (i.e. new structures 
    must be retrieved)
    """

    this_name = check_update_repo.__name__

    # Options and locations aliases
    OPT_download_date_threshold =\
        int(options['RUN'][('', 'download_date_threshold')]) # Days from last update
    OPT_force_update = options['RUN'][('fud', 'force_update')]
    LOC_OPMlastddate = locations['SYSFILES']['OPMlastddate']
    LOC_main = locations['FSYSPATH']['main']

    # Checks current time and date
    update_repo = False
    now_datetime = datetime.datetime.now()
    now_timestamp = time.mktime(now_datetime.timetuple()) # Epoch time in seconds

    # If "force update" option is present, return True
    if OPT_force_update:
        print_log((
            'NOTICE', 
            this_name, 
            f"Force update contents of the folder {LOC_main}"))
        update_repo = True
    # If time from last update > download_date_threshold option, return True
    elif os.path.exists(LOC_OPMlastddate):
        # Time span conversion into days
        with open(LOC_OPMlastddate) as ldd_file:
            for line in ldd_file:
                last_download_timestamp = line.split()[0] # Epoch time in seconds
                len_ldt = len(last_download_timestamp)
                last_download_datetime = line[len_ldt-1:].strip()
        time_span_days = int((int(last_download_timestamp) - int(now_timestamp))
            /(24*3600))
        # If time from last update > download_date_threshold opt, return True
        if time_span_days > OPT_download_date_threshold:
            print_log((
                'NOTICE', 
                this_name, 
                ("Update contents of the folder {0}\nLast update: {1}"
                    "\nDeadline: {2} days"
                ).format(
                    LOC_main, 
                    last_download_datetime,
                    OPT_download_date_threshold)
            ))
            update_repo = True
    # Else, it must be the first update. 
    else:
        print_log((
            'NOTICE', 
            this_name, 
            ("First update contents of the folder {0} (or file {1} "
                "not found)").format(
                LOC_main, 
                LOC_OPMlastddate)
        ))
        update_repo = True

    # If the date is not recorded, the update will be issued again next time
    if record:
        print_log((
            'NOTICE', 
            this_name, 
            "Update complete, new date recorded"))
        with open(LOC_OPMlastddate, 'w') as repo_date_file:
            repo_date_file.write("{0} {1}\n".format(
                int(now_timestamp), 
                now_datetime.strftime('%Y-%m-%d %H:%M:%S')))
    return update_repo


def mainlist_from_OPM(options, locations, str_data={}, 
        use_pkl=False, # Deprecated?
        thr_log_status='ERROR'):
    """Returns the list of primary and secondary OPM entries to be considered
    during this instance of EncoMPASS
    """

    this_name = mainlist_from_OPM.__name__

    # Options and locations aliases
    OPT_from_repo = options['RUN'][('repo', 'from_repo')]
    OPT_offline = options['RUN'][('ol', 'offline')]
    OPT_sample_structures = options['EXECUTION'][('sample', 'sample_structures')]
    OPT_debug = options['EXECUTION'][('debug', 'debug')] 
    LOC_OPMjsons = locations['FSYSPATH']['OPMjsons']
    LOC_OPMreprchart = locations['SYSFILES']['OPMreprchart']

    # Statistics
    local_stats = {
        'primary' : [],
        'secondary' : [],
        'associated_secondary' : [],
        'selected_primary' : [],
        'selected_secondary' : [],
        'temp_primary' : [],
        'sample_not_found' : []
    }

    # Check valid option setting
    if OPT_offline and (not OPT_from_repo):
        print_log((
            'CRITICAL',
            this_name,
            ("Offline option is on: please use from_repo option to specify a "
                "local repository where to take the OPM entry list from")
        ))

    # Check whether OPM must be queried
    update_repo = check_update_repo(options, locations)

    # Download/retrieve main json file with PDB codes of primary structures
    mainjson_filename = LOC_OPMjsons + 'main.json'
    repofilename = OPT_from_repo + LOC_OPMjsons + 'main.json'
    if OPT_from_repo and os.path.exists(repofilename):
        shutil.copyfile(repofilename, mainjson_filename)
    else:
        res = requests.get(
            # The .../types/1/... specifies this is the list for type-1 OPM structures (Transmembrane)
            "https://lomize-group-opm.herokuapp.com/types/1/primary_structures?pageSize=100000",
            params={},
            headers={"Accept" : "application/json"}
        )
        with open(mainjson_filename, 'w') as f:
            data = res.json()
            f.write(json.dumps(data))
            local_stats['primary'] = [x['pdbid'] for x in data['objects']]

    if not os.path.exists(mainjson_filename):
        print_log((
            'CRITICAL',
            this_name,
            "Could not retrieve the list of primary OPM structures"
        ))
                
    #  Extract primary entry info from json file
    if OPT_debug:
        print("Reading JSON file", mainjson_filename)
    with codecs.open(mainjson_filename, 'r+', encoding="utf-8") as mainjson_f:
        data = json.load(mainjson_f, encoding='utf-8')
    ids = [(x['id'], x['pdbid'], False) for x in data['objects']] # Primary (<int>ID, <str>PDBID, <bool>is_temp) list 
    ids_d = {x['id'] : x['pdbid'] for x in data['objects']} # Primary <int>ID-><str>PDBID dict (internal use)
    revids_d = {x['pdbid'] : x['id'] for x in data['objects']} # Primary <str>PDBID-><int>ID dict (internal use)

    # Download/retrieve secondary structure json file
    secjson_filename = LOC_OPMjsons + 'secondary.json'
    repofilename = OPT_from_repo + LOC_OPMjsons + 'secondary.json'
    if OPT_from_repo and os.path.exists(repofilename):
        shutil.copyfile(repofilename, secjson_filename)
    else:
        if not OPT_offline:
            res = requests.get(
                "https://lomize-group-opm.herokuapp.com/secondary_representations?pageSize=1000000",
                params={},
                headers={"Accept" : "application/json"}
            )
            with open(secjson_filename, 'w') as f:
                data = res.json()
                f.write(json.dumps(data))
                local_stats['secondary'] = [x['pdbid'] for x in data['table']['objects']]
    if not os.path.exists(secjson_filename):
        print_log((
            'CRITICAL',
            this_name,
            "Could not retrieve the list of secondary OPM structures"
        ))

    #  Extract secondary entries info from json file
    with codecs.open(secjson_filename, 'r+', encoding="utf-8") as secjson_f:
        data = json.load(secjson_f, encoding='utf-8')
    secrep = [] # Secondary <str>PDBID list (internal use)
    secrep_d = {} # Secondary <str>PDBID->[<str>PDBID, ...] dict
    revsecrep_d = {} # Secondary <str>PDBID-><int>ID dict (internal use)
    for x in data['table']['objects']:
        if x['primary_structure_id'] not in ids_d:
            #print("CONT", x['primary_structure_id'])
            continue
        if ids_d[x['primary_structure_id']] not in secrep_d:
            secrep_d[ids_d[x['primary_structure_id']]] = []
        if x['pdbid'] not in local_stats['associated_secondary']:
            local_stats['associated_secondary'].append(x['pdbid'])
        secrep.append(x['pdbid'])
        secrep_d[ids_d[x['primary_structure_id']]].append(x['pdbid'])
        revsecrep_d[x['pdbid']] = x['primary_structure_id']

    # Only choose the structures in the sample list. Primary structures 
    # corresponding to cited secondary ones will be added and then removed 
    # after the OPM data retrieval
    if OPT_sample_structures:
        new_ids = [] # Primary (<int>ID, <str>PDBID, <bool>is_temp) list
        new_secrep_d = {} # Secondary <int>ID->[<str>PDBID, ...] dict
        ids1 = [x[1] for x in ids]
        for entry in OPT_sample_structures:
            if entry[:4] in ids1:
                # If a temporary primary entry was added because of an orphan secondary entry,
                # but later in the list the same entry is read as a rightful entry, switch off
                # the is_temp bit
                if (revids_d[entry[:4]], entry[:4], True) in new_ids:
                    idx = new_ids.index((revids_d[entry[:4]], entry[:4], True))
                    new_ids[idx] = (revids_d[entry[:4]], entry[:4], False)
                else:
                    new_ids.append((revids_d[entry[:4]], entry[:4], False))
                # If it is a "+" entry, add all the secondary representations
                if len(entry) == 5 and entry[4] == "+":
                    if entry[:4] in secrep_d: # entry is absent if it does not have secondary reprs
                        new_secrep_d[entry[:4]] = secrep_d[entry[:4]]
                    else:
                        print_log((
                            'WARNING',
                            this_name,
                            ("Entry {0} was included in sample_structures, yet "
                                "{1} does not have secondary representations")\
                                .format(entry, entry[:4])
                        ))
            elif entry[:4] in secrep:
                repr_id = revsecrep_d[entry[:4]]
                # Add the primary entry corresponding to the orphan secondary entry
                if (repr_id, ids_d[repr_id], False) not in new_ids:
                    # This entry is only temporary and should not be considered as an EncoMPASS entry
                    new_ids.append((repr_id, ids_d[repr_id], True)) 
                if repr_id not in new_secrep_d:
                    new_secrep_d[ids_d[repr_id]] = []
                new_secrep_d[ids_d[repr_id]].append(entry[:4])
            else:
                local_stats['sample_not_found'].append(entry)
        local_stats['temp_primary'] = [x[1] for x in new_ids if x[2] is True]
        ids = new_ids
        secrep_d = new_secrep_d

    local_stats['selected_primary'] = [x[1] for x in ids]
    for x in secrep_d:
        local_stats['selected_secondary'] += secrep_d[x]

    # Write representation chart of all the entries considered in this EncoMPASS repository
    write_OPM_representation_chart(ids, secrep_d, LOC_OPMreprchart)
    
    # Show statistics
    for x in sorted(local_stats):
        print("STATS", "mainlist_from_OPM", x, len(local_stats[x]), local_stats[x])
    print("STATS", "mainlist_from_OPM", "Finished", "time", time.ctime(), "\n")

    # List of primary entries [(<int>ID, <str>PDBID, <bool>is_temp), ...] and dictionary <str>PDBID->[<str>PDBID, ...] linking each primary ID to list of secondary PDBIDs
    return ids, secrep_d 


def scan_OPM(options, locations, entries, str_data={}, use_pkl=False, 
        thr_log_status='ERROR'):

    this_name = scan_OPM.__name__

    # Options and locations aliases
    OPT_dst = options['PATHS'][('dst', 'data_structure_template')]
    OPT_from_repo = options['RUN'][('repo', 'from_repo')]
    OPT_offline = options['RUN'][('ol', 'offline')]
    OPT_debug = options['EXECUTION'][('debug', 'debug')]
    OPT_sample_structures =\
        options['EXECUTION'][('sample', 'sample_structures')]
    LOC_OPM = locations['FSYSPATH']['OPM']
    LOC_OPMjsons = locations['FSYSPATH']['OPMjsons']
    LOC_OPMpdbs = locations['FSYSPATH']['OPMpdbs']
    LOC_UniProt = locations['FSYS']['UniProt']
    LOC_uniprot_all = locations['SYSFILES']['uniprot_all']
    LOC_OPMreprchart = locations['SYSFILES']['OPMreprchart']

    # Statistics
    local_stats = {
        'json_not_found' : [],
        'no_segments' : [],
        'no_alphabeta' : [],
        'OPMpdb_not_found' : []
    }

    # Update repository
    update_repo = check_update_repo(options, locations)


    # Download/retrieval of OPM infos about primary entries
    print("scan_OPM", "Download json files", "time", time.ctime())
    primary_structures_data = {}
    for i, pdbi, is_temp in sorted(entries, key=lambda x:x[1]):
        if OPT_debug:
            print(pdbi)
        json_filename = LOC_OPMjsons + pdbi + '.json'
        repofilename = OPT_from_repo + LOC_OPMjsons + pdbi + '.json'
        if OPT_from_repo and os.path.exists(repofilename):
            shutil.copyfile(repofilename, json_filename)
        else:
            if ((not OPT_offline) and
                    (update_repo or (not os.path.exists(json_filename)))):
                res = requests.get(
                    f"https://lomize-group-opm.herokuapp.com/primary_structures/{i}",
                    params={},
                    headers={"Accept" : "application/json"}
                )
                time.sleep(0.1)
                if not os.path.exists(json_filename):
                    with open(json_filename, 'w') as f:
                        f.write(json.dumps(res.json()))

        if not os.path.exists(json_filename):
            print_log((
                'ERROR', 
                this_name,
                f"Could not find file {json_filename}"
            ))
            local_stats['json_not_found'].append(pdbi)
        else:
            with codecs.open(json_filename, 'r+', encoding='utf-8')\
                    as json_file:
                primary_structures_data[(i, pdbi, is_temp)] = json.load(json_file, 
                    encoding='utf-8')

    # Download OPM coordinate files
    print("scan_OPM", "Download coordinate files", "time", time.ctime())
    for i, pdbi, is_temp in sorted(primary_structures_data, key=lambda x:x[1]):
        if options['EXECUTION'][('debug', 'debug')]:
            print(pdbi)
        if pdbi not in str_data:
            str_data[pdbi] = define_str_data_entry(OPT_dst)

        # Update status
        str_data[pdbi]['status'].append('initial_opm')

        # EncoMPASS name is set to OPM name
        str_data[pdbi]['ENCOMPASS']['name'] =\
            primary_structures_data[(i, pdbi, is_temp)]['name']
        str_data[pdbi]['PASSPORT'].append(passport_entry(
            this_name+'_encname_opm', 
            pdbi, 
            "The name of this structure follows the OPM nomenclature"
        ))

        # Download/retrieve the OPM coordinate files
        pdb_path = LOC_OPMpdbs + pdbi + '.pdb'
        repofilename = OPT_from_repo + LOC_OPMpdbs + pdbi + '.pdb'
        if OPT_from_repo and os.path.exists(repofilename):
            shutil.copyfile(repofilename, pdb_path)
        else:
            if not os.path.exists(pdb_path) and not OPT_offline: 
                download_result = reiterated_simple_download(
                    f'https://opm-assets.storage.googleapis.com/pdb/{pdbi}.pdb', 
                    pdb_path
                )

        # Remove 0-segment structures
        if int(primary_structures_data[(i, pdbi, is_temp)]['subunit_segments']) == 0:
            str_data[pdbi]['PASSPORT'].append(passport_entry(
                this_name+'_zero_tmsegments', 
                pdbi, 
                ("According to OPM, this structure has 0 TM segments, i.e. "
                    "it is not a transmembrane protein structure. The entry "
                    "will not be considered in EncoMPASS.")
            ))
            str_data[pdbi]['status'].append('eliminated')
            str_data[pdbi]['delete_keyword'] = 'Monotopic'
            local_stats['no_segments'].append(pdbi)

        # Sort attributes
        for x in str_data[pdbi]['FROM_OPM']:
            # Copy all superfamily sub-attributes
            if 'superfamily' in x:
                for y in str_data[pdbi]['FROM_OPM']['superfamily']:
                    str_data[pdbi]['FROM_OPM']['superfamily'][y] =\
                        primary_structures_data[(i, pdbi, is_temp)]['family']['superfamily'][y]
            # Class renaming: alpha-beta-other
            # TO DO add a statistic here
            elif 'class' in x:
                classtype = primary_structures_data[(i, pdbi, is_temp)]['family']\
                    ['superfamily']['classtype']['name']
                str_data[pdbi]['FROM_OPM']['classtype'] =\
                    primary_structures_data[(i, pdbi, is_temp)]['family']\
                        ['superfamily']['classtype']['name']
                if classtype == 'Alpha-helical polytopic':
                    enc_class = 'alpha'
                elif classtype == 'Beta-barrel transmembrane':
                    enc_class = 'beta'
                else:
                    local_stats['no_alphabeta'].append(pdbi)
                    enc_class = 'other'
                str_data[pdbi]['ENCOMPASS']['class'] = enc_class
                str_data[pdbi]['PASSPORT'].append(passport_entry(
                    this_name+'_encclass_opm',
                    pdbi,
                    ("This structure is classified according to the OPM "
                        "topology classification, in the category: {0}"
                    ).format(enc_class)
                ))
            # Check whether the OPM coordinate file is there
            elif x == 'coordinate_file':
                if os.path.exists(pdb_path):
                    str_data[pdbi]['FROM_OPM']['coordinate_file'] = pdb_path
                else:
                    local_stats['OPMpdb_not_found'].append(pdbi)
                    print_log((
                        'ERROR', 
                        this_name, 
                        "OPM pdb file {0} was not found in {1}".format(
                            pdbi, locations['FSYSPATH']['OPMpdbs'])
                    ))
                    if ((not str_data[pdbi]['status']) or 
                            str_data[pdbi]['status'][-1] != 'opm_eliminated'):
                        str_data[pdbi]['status'].append('opm_eliminated')
                    if ((not str_data[pdbi]['status']) or 
                            str_data[pdbi]['status'][-1] != 'run_ppm'):
                        str_data[pdbi]['status'].append('run_ppm')
            # Copy other dictionaries
            elif (isinstance(str_data[pdbi]['FROM_OPM'][x], FixedDict) and 
                    x != 'analysis'):
                for y in str_data[pdbi]['FROM_OPM'][x]:
                    str_data[pdbi]['FROM_OPM'][x][y] =\
                        primary_structures_data[(i, pdbi, is_temp)][x][y]
            # Copy other lists
            elif isinstance(str_data[pdbi]['FROM_OPM'][x], FixedList):
                # If there are no secondary representation, cycle
                if (x == 'secondary_representations' and 
                        (not primary_structures_data[(i, pdbi, is_temp)][x])):
                    continue
                # Chains as described in OPM
                # The 'subunit' key in OPM json contains the infos regarding TM segments
                if x == 'subunits':
                    for subunit in primary_structures_data[(i, pdbi, is_temp)][x]:
                        subunits_key = subunit['protein_letter']
                        str_data[pdbi]['FROM_OPM']['ksubunits'].append(
                            subunits_key)
                        tempd = FixedDict(
                            str_data[pdbi]['FROM_OPM']['subunits']\
                            .get_fdict_template()
                        )
                        for y in tempd:
                            # Fill OPM segments
                            if y == "segment":
                                tempd[y] = []
                                for seg in subunit[y].split('(')[1:]:
                                    n = seg.split('-')
                                    if len(n) == 2:
                                        n1, n2 = seg.split('-')
                                        n1 = int(''.join(n1.split()))
                                        n2 = int(''.join(
                                            n2.split(')')[0].split())
                                        )
                                    elif len(n) == 1:
                                        n1, n2 = n[0][:n[0].index(')')].split()
                                    tempd[y].append(tuple([n1, n2]))
                                #print("OPM SEGMENTS", pdbi, subunits_key, tempd[y])
                            else:
                                tempd[y] = subunit[y]
                        str_data[pdbi]['FROM_OPM']['subunits']\
                            .append(tempd)
                        #print("WRITE", tempd.show_dictionary())
                # Otherwise, OPM json and the str_data have the same keys and data can be copied
                else:
                    for il, l in enumerate(primary_structures_data[(i, pdbi, is_temp)][x]):
                        tempd = FixedDict(
                            str_data[pdbi]['FROM_OPM'][x].get_fdict_template())
                        for y in tempd:
                            tempd[y] = primary_structures_data[(i, pdbi, is_temp)][x][il][y]
                        str_data[pdbi]['FROM_OPM'][x].append(tempd)
            # Copy other non-dictionary, non-list data
            else:
                if x in primary_structures_data[(i, pdbi, is_temp)]:
                    str_data[pdbi]['FROM_OPM'][x] =\
                        primary_structures_data[(i, pdbi, is_temp)][x]
        # Check the primary representation case
        str_data[pdbi]['FROM_OPM']['primary_representation'] = True

        # Secondary representations
        if len(str_data[pdbi]['FROM_OPM']['secondary_representations']) > 0:
            # Loop on the codes found in the OPM json page
            sec_pdbis = [x['pdbid'] for x in\
                primary_structures_data[(i, pdbi, is_temp)]['secondary_representations']]
            for sec_pdbi in sec_pdbis:
                # If they are not among the codes to sample, cycle
                if (OPT_sample_structures and 
                        sec_pdbi not in OPT_sample_structures):
                    continue
                # Retrieve the OPM coordinate file
                sec_pdb_path = LOC_OPMpdbs + sec_pdbi + '.pdb'
                repofilename = OPT_from_repo + LOC_OPMpdbs + sec_pdbi + '.pdb'
                if OPT_from_repo and os.path.exists(repofilename):
                    shutil.copyfile(repofilename, sec_pdb_path)
                else:
                    if not OPT_offline and not os.path.exists(sec_pdb_path):
                        download_result = reiterated_simple_download(
                            f'https://opm-assets.storage.googleapis.com/pdb/{sec_pdbi}.pdb', 
                            sec_pdb_path
                        )
                # Copy the str_data entry of the pdbi entry into sec_pdbi -> DANGEROUS: CHECK ALL PATHS
                str_data[sec_pdbi] = FixedDict(str_data[pdbi])


                if str_data[pdbi]['ENCOMPASS']['class'] == 'other':
                    local_stats['no_alphabeta'].append(sec_pdbi)
                
                # OPM coordinate file analysis
                if os.path.exists(sec_pdb_path):
                    str_data[sec_pdbi]['FROM_OPM']\
                        ['coordinate_file'] = sec_pdb_path
                    str_data[sec_pdbi]['FROM_OPM']\
                        ['primary_representation'] = False

                    # Add as secondary representations all other secondary 
                    # representations of its primary representation except
                    # itself, plus the primary representation
                    str_data[sec_pdbi]['FROM_OPM']['secondary_representations'] =\
                        str_data[sec_pdbi]['FROM_OPM']\
                        .show_dictionary(quiet=True)['secondary_representations']
                    for iixx, xx in enumerate(str_data[sec_pdbi]['FROM_OPM']\
                            ['secondary_representations']):
                        if xx['pdbid'] == sec_pdbi:
                            iitodel = iixx
                            break
                    del(str_data[sec_pdbi]['FROM_OPM']\
                        ['secondary_representations'][iitodel])

                    elm = str_data[sec_pdbi]['FROM_OPM']\
                        ['secondary_representations'].get_fdict_template()
                    for ii, dd, tt in entries:
                        if dd == pdbi:
                            elm['id'] = ii
                            break
                    elm['pdbid'] = pdbi
                    str_data[sec_pdbi]['FROM_OPM']\
                        ['secondary_representations'].append(FixedDict(elm))
                else:
                    local_stats['OPMpdb_not_found'].append(sec_pdbi)
                    print_log((
                        'ERROR',
                        this_name,
                        "OPM coordinate file {0} was not found in {1}"\
                            .format(sec_pdbi, LOC_OPMpdbs)
                    ))
                    if ((not str_data[sec_pdbi]['status']) or
                            str_data[sec_pdbi]['status'][-1] != 'opm_eliminated'):
                        str_data[sec_pdbi]['status'].append('opm_eliminated')
                    if ((not str_data[sec_pdbi]['status']) or
                            str_data[sec_pdbi]['status'][-1] != 'run_ppm'):
                        str_data[sec_pdbi]['status'].append('run_ppm')

        # If it does not have subunits, put in run_ppm list
        if len(str_data[pdbi]['FROM_OPM']['subunits']) == 0:
                str_data[pdbi]['status'].append('run_ppm')
        # Remove primary entries that were added only in order to describe some orphan secondary entry
        if is_temp:
            del(str_data[pdbi])
            print_log((
                'NOTICE',
                this_name,
                "Entry {0} was only a proxy, it will be deleted".format(pdbi)
            ))

    # If from_repo, copy uniprot data from repository
    if OPT_from_repo:
        repofilename = OPT_from_repo + LOC_UniProt + os.path.basename(LOC_uniprot_all)
        print(repofilename, LOC_uniprot_all)
        shutil.copyfile(repofilename, LOC_uniprot_all)

    # Make pickle
    str_data_pkl = LOC_OPM + '.str_data_after_OPM_read.pkl'
    pickle.dump(str_data, open(str_data_pkl, 'wb'))
    # Record current date 
    check_update_repo(options, locations, record=True)

    # Show statistics
    for x in sorted(local_stats):
        print("STATS", "scan_OPM", x, len(local_stats[x]), local_stats[x])
    print("STATS", "scan_OPM", "Finished", "time", time.ctime(), "\n")
    return str_data


def scan_PDB(options, locations, str_data={}, pdb_list=[], use_pkl=False, thr_log_status='ERROR'):
    """scan_PDB
    Collects information from the RCSB GraphQL API and records them in the data
    structure (FROM_PDB key)
    """

    this_name = scan_PDB.__name__

    # Options and locations aliases
    local_stats = {
        'no_record' : [],
        'no_str' : [],
        'no_title' : [],
        'no_expmethod' : [],
        'no_resolution' : [],
        'no_depo_date' : [],
        'no_entities' : []
    }

    if not pdb_list:
        pdb_list = sorted([x for x in str_data])

    # Pickle
    str_data_pkl = locations['FSYSPATH']['PDB'] + '.str_data_after_PDB_read.pkl'
    if use_pkl and os.path.exists(str_data_pkl):
        str_data = pickle.load(open(str_data_pkl, 'rb'))

    # Obsolete PDB entries
    obsolete = set()
    obsolete_filename = locations['FSYSPATH']['PDBjsons'] + 'obsolete_entries.json'

    #  Get them from repo or server
    repofilename = options['RUN'][('repo', 'from_repo')] + locations['FSYS']['PDBjsons'] + 'obsolete_entries.json'
    if options['RUN'][('repo', 'from_repo')] and os.path.exists(repofilename):
        shutil.copyfile(repofilename, obsolete_filename)
    else:
        pass

    # Json info records
    temparch = {}  # Keys 'coordinate_file' and 'header' do not come from json but from structure

    # pdbis -> batches
    batch_size = 100
    batches = []
    with open(locations['FSYSPATH']['cache'] + 'PDBall.txt', 'w') as f:
        for ipdbi, pdbi in enumerate(pdb_list):
            if ipdbi%batch_size == 0:
                batches.append([])
            batches[-1].append(pdbi.upper())
            f.write(pdbi.upper()+'\n')

    # query by batches
    url = 'https://data.rcsb.org/graphql?'
    headers = []
    print("scan_PDB", "Querying PDB by batches", "time", time.ctime())
    for ibatch, batch in enumerate(batches):
        if options['EXECUTION'][('debug', 'debug')]:
            print("Batch", ibatch)
        # Batch query
        json_query = pdb_query(path=locations['SYSFILES']['pdbquery'])
        json_query = json_query.replace("<PDBI>", ", ".join(['"'+x+'"' for x in batch]))
        if options['EXECUTION'][('debug', 'debug')]:
            print("JSON query file", json_query)
        r = requests.post(url, json={'query': json_query})
        js = json.loads(r.text)
        if options['EXECUTION'][('debug', 'debug')]:
            print("JSON data", js)
        ebatch = []
        for e in js['data']['entries']:
            ebatch.append(e['rcsb_entry_container_identifiers']['entry_id'].lower())

        # Download and parse individual structures
        for pdbi in batch:
            if options['EXECUTION'][('debug', 'debug')]:
                print(pdbi)

            if pdbi.lower() not in ebatch:
                local_stats['no_record'].append(pdbi.lower())
                continue
            ipdbi = ebatch.index(pdbi.lower())
            with open(locations['FSYSPATH']['PDBjsons'] + pdbi.lower() + '.json', 'w') as jf:
                json.dump(js['data']['entries'][ipdbi], jf)
            temparch[pdbi.lower()] = js['data']['entries'][ipdbi]
            pdb_filename = locations['FSYSPATH']['PDBpdbs'] + pdbi.lower() + '.pdb'
            mmcif_filename = locations['FSYSPATH']['PDBmmcifs'] + pdbi.lower() + '.cif'

            # Download coordinate file
            if options['RUN'][('mmcifmirror', 'mmcif_local_mirror')]:
                mirrormmcif_fn = options['RUN']['mmcifmirror', 'mmcif_local_mirror'] + pdbi.upper() + '.cif'
                mirrormmcif_fn2 = options['RUN']['mmcifmirror', 'mmcif_local_mirror'] + pdbi.lower() + '.cif'
                if os.path.exists(mirrormmcif_fn):
                    shutil.copyfile(mirrormmcif_fn, mmcif_filename)
                elif os.path.exists(mirrormmcif_fn2):
                    shutil.copyfile(mirrormmcif_fn2, mmcif_filename)
                else:
                    str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_no_pdbmirror', pdbi, "This entry was not found in the pdb mirror"))
            elif options['RUN'][('pdbmirror', 'pdb_local_mirror')]:
                mirrorpdb_fn = options['RUN']['pdbmirror', 'pdb_local_mirror'] + '/' + pdbi.lower()[1:3] + '/pdb' + pdbi.lower() + '.ent.gz'
                if os.path.exists(mirrorpdb_fn):
                    shutil.copyfile(mirrorpdb_fn, pdb_filename+'.gz')
                os.system('gunzip ' + pdb_filename + '.gz')
            elif options['RUN'][('repo', 'from_repo')]:
                repofilename = options['RUN'][('repo', 'from_repo')] + locations['FSYS']['PDBpdbs'] + pdbi.lower() + '.pdb'
                if os.path.exists(repofilename):
                    shutil.copyfile(repofilename, pdb_filename)
            else:
                download_result = reiterated_simple_download('https://files.rcsb.org/download/{0}.cif'.format(pdbi.upper()), mmcif_filename)
                download_result_2 = reiterated_simple_download('https://files.rcsb.org/download/{0}.pdb'.format(pdbi.upper()), pdb_filename)

            # Parse header
            if os.path.exists(mmcif_filename) or os.path.exists(mmcif_filename.replace(".gz", "")+"gz"):
                temparch[pdbi.lower()]['coordinate_file'] = mmcif_filename
                header, hreport = parse_mmCIF_header_wrap(mmcif_filename)
                temparch[pdbi.lower()]['header'] = header
                headers.append(temparch[pdbi.lower()]['header'])  # USED NOW! For deciding alpha-beta-other class
            elif os.path.exists(pdb_filename) or os.path.exists(pdb_filename.replace(".gz", "")+"gz"):
                temparch[pdbi.lower()]['coordinate_file'] = pdb_filename
                temparch[pdbi.lower()]['header'], hreport = parse_PDB_header(pdb_filename) # not used so far
            else:
                local_stats['no_str'].append(pdbi.lower())

    mmCIF_header_stats(locations['FSYSPATH']['PDBstats'], headers)
            
    # Fill FROM_PDB with info collected in temparch and with those in the PDB header
    print("scan_PDB", "Fill data structure", "time", time.ctime())
    for pdbi in pdb_list:
        if options['EXECUTION'][('debug', 'debug')]:
            print(pdbi)
        if pdbi in local_stats['no_record']:
            str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_no_graphql', pdbi, "The GraphQL entry for {0} was not found. The PDB entry will be disregarded".format(pdbi)))
            str_data[pdbi]['status'].append('pdb_eliminated')
            continue
        if pdbi in local_stats['no_str']:
            str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_no_coordfile', pdbi, "The PDB coordinate file for {0} was not found. The PDB entry will be disregarded".format(pdbi)))
            str_data[pdbi]['status'].append('pdb_eliminated')
            continue
        from_pdb = FixedDict(str_data[pdbi]['FROM_PDB'].get_fdict_template())  # Template FROM_PDB structure

        from_pdb['coordinate_file'] = temparch[pdbi.lower()]['coordinate_file']
        if not from_pdb['coordinate_file']:
            str_data[pdbi]['status'].append('pdb_eliminated')

        from_pdb['ss_class'] = temparch[pdbi]['header']['ss_class']

        if not str_data[pdbi]['status'] or str_data[pdbi]['status'][-1] != 'initial_pdb':
            str_data[pdbi]['status'].append('initial_pdb')

        # General info
        from_pdb['pdb_id'] = pdbi
        if temparch[pdbi]['struct'] is not None:
            from_pdb['name'] = temparch[pdbi]['struct']['title']
        else:
            wnores_msg = ('WARNING', this_name, f"PDB json entry {pdbi} does not have the 'struct' key. Name will be set to None")
            print_log(wnores_msg)
            from_pdb['name'] = None
            local_stats['no_title'].append(pdbi)

        if temparch[pdbi]['exptl'] is not None:
            from_pdb['experimental_method'] = ",".join(x['method'] for x in temparch[pdbi]['exptl'])
            wnores_msg = ('WARNING', this_name, f"PDB json entry {pdbi} does not have the 'exptl' key. Experimental method will be set to None")
            print_log(wnores_msg)
        else:
            from_pdb['experimental_method'] = None
            local_stats['no_expmethod'].append(pdbi)

        if temparch[pdbi]['refine']: # It is a list and can be empty
            # If no 'refine', there are 3 keywords missing, not just 'PDB_resolution'
            from_pdb['resolution'] = temparch[pdbi]['refine'][0]['ls_d_res_high']
        else:
            if temparch[pdbi]['header']['resolution']:
                from_pdb['resolution'] = temparch[pdbi]['header']['resolution']
            else:
                local_stats['no_resolution'].append(pdbi)
            wnores_msg = ('WARNING', this_name, f"PDB json entry {pdbi} does not have the 'refine' key. Resolution will be set to None")
            print_log(wnores_msg)
            from_pdb['resolution'] = None
        if temparch[pdbi]['pdbx_database_status'] is not None:
            print(temparch[pdbi]['pdbx_database_status'])
            from_pdb['deposition_date'] = temparch[pdbi]['pdbx_database_status']['recvd_initial_deposition_date']
        else:
            local_stats['no_depo_date'].append(pdbi)
            wnores_msg = ('WARNING', this_name, f"PDB json entry {pdbi} does not have the 'database_PDB_rev' key. Deposition date will be set to None")
            print_log(wnores_msg)
            from_pdb['deposition_date'] = None

        # Polymers
        new2oldch = {}
        polychains = set()
        alreadythere = set()
        if temparch[pdbi]['polymer_entities'] is not None:
            for ip, p in enumerate(temparch[pdbi]['polymer_entities']):
                chains = temparch[pdbi]['polymer_entities'][ip]['entity_poly']['pdbx_strand_id']
                for inst in p['polymer_entity_instances']:
                    if inst['rcsb_polymer_entity_instance_container_identifiers']['auth_asym_id'] in polychains: # HETATMs of the same chain are now considered a separate chain, we don't want to have them
                        continue
                    if inst['rcsb_polymer_entity_instance_container_identifiers']['auth_asym_id'] in chains:
                        new2oldch[inst['rcsb_polymer_entity_instance_container_identifiers']['asym_id']] = inst['rcsb_polymer_entity_instance_container_identifiers']['auth_asym_id']
                        polychains.add(inst['rcsb_polymer_entity_instance_container_identifiers']['auth_asym_id'])
                    else:
                        exit(1)
                for ch in sorted(list(polychains)):
                    poly = FixedDict(from_pdb['polymers'].get_fdict_template())
                    poly['chain'] = ch
                    poly['sequence'] = temparch[pdbi]['polymer_entities'][ip]['entity_poly']['pdbx_seq_one_letter_code_can']
                    if type(temparch[pdbi]['polymer_entities'][ip]['rcsb_polymer_entity_container_identifiers']['uniprot_ids']) == list:
                        poly['uniprot_acc'] =  ",".join(temparch[pdbi]['polymer_entities'][ip]['rcsb_polymer_entity_container_identifiers']['uniprot_ids'])
                    if type(temparch[pdbi]['polymer_entities'][ip]['rcsb_entity_source_organism']) == list:
                        poly['taxonomy'] = ",".join([x['ncbi_scientific_name'] if x['ncbi_scientific_name'] else "" for x in temparch[pdbi]['polymer_entities'][ip]['rcsb_entity_source_organism']])
                    from_pdb['polymers'].append(poly)
        else:
            local_stats['no_entities'].append(pdbi)
            wnores_msg = ('WARNING', this_name, f"PDB json entry {pdbi} does not have the 'polymer_entities' key. Polymer section will be empty")
            print_log(wnores_msg)

        # Biological assembly
        if temparch[pdbi]['assemblies'] is not None:
            # Assemblies
            # NB: an assembly is a way to form a complex. Different assemblies are found when different hypotheses on the
            #  biological assembly have been formulated. Usually, we consider the first assembly (but we record all)
            oper_tr = FixedDict(from_pdb['translrots'].get_fdict_template())
            for ib, b in enumerate(temparch[pdbi]['assemblies']):  
                ba = FixedDict(from_pdb['biological_assemblies'].get_fdict_template())
                ba['method'] = temparch[pdbi]['assemblies'][ib]['pdbx_struct_assembly']['method_details']
                ba['description'] = temparch[pdbi]['assemblies'][ib]['pdbx_struct_assembly']['details']
                # Transformations
                # NB: a transformation is a SET of matrices that applies on the SAME SET of chains. An assembly can have more than one transformation
                # (such as 3dh4), because different matrices can apply to different chains and still contribute giving the same complex
                for itr in range(len(temparch[pdbi]['assemblies'][ib]['pdbx_struct_assembly_gen'])):
                    tr = FixedDict(ba['transformations'].get_fdict_template())
                    tr['chains'] = [new2oldch[x] for x in temparch[pdbi]['assemblies'][ib]['pdbx_struct_assembly_gen'][itr]['asym_id_list'] if x in new2oldch] # List of chains(new name) for each assembly matrix)
                    oper_list, oper_keys = read_transf_matrices_json(temparch[pdbi]['assemblies'][ib]['pdbx_struct_oper_list'])
                    opexps = read_oper_expression(temparch[pdbi]['assemblies'][ib]['pdbx_struct_assembly_gen'][itr]['oper_expression'], oper_keys)
                    #print("READ_OPER_EXPRESSION", pdbi, ib, itr, temparch[pdbi]['assemblies'][ib]['pdbx_struct_assembly_gen'][itr]['oper_expression'], oper_keys, opexps)
                    # Rototranslations are individual pairs of rotation matrix and translation vectors
                    #  They are the building blocks of operators
                    # NB: first the matrix will be applied, then the vector.
                    # NB2: in json files, translrots are listed as depending on each operator, but in mmCIF files they are not
                    #  (i.e. their id is UNIQUE throughout the file, it does not reset at each operator change
                    #  Here, we keep the more synthetic (and correct?) mmCIF notation
                    for ioper, oper in enumerate(oper_list):
                        op = FixedDict(from_pdb['translrots'].get_fdict_template())
                        op['matrix'] = oper['matrix']
                        op['vector'] = oper['vector']
                        # Append ONLY if ID is not already there
                        if oper_keys[ioper] not in from_pdb['ktranslrots']:
                            from_pdb['ktranslrots'].append(oper_keys[ioper])
                            from_pdb['translrots'].append(op)
                    # Operators are individual rototranslation functions. They can combine multiple rototranslators
                    # Example: [(M1,v1),(M2,v2)] is the function v -> (M2*((M1*v)+v1))+v2
                    # NB: there exists always M*,v* such that [(M1,v1),...,(Mn,vn)] = [(M*,v*)], but I don't know how to calculate it
                    tr['operators'] = opexps
                    ba['transformations'].append(tr)
                from_pdb['biological_assemblies'].append(ba)
        else:
            ba = FixedDict(from_pdb['biological_assemblies'].get_fdict_template())
            tr = FixedDict(ba['transformations'].get_fdict_template())
            tr['chains'] = sorted(list(polychains))
            tr['operators'] = [["1"]]
            ba['transformations'].append(tr)
            from_pdb['biological_assemblies'].append(ba)

            op = FixedDict(from_pdb['translrots'].get_fdict_template())
            op['matrix'] = [[1,0,0],[0,1,0],[0,0,1]]
            op['vector'] = [0, 0, 0]
            from_pdb['translrots'].append(op)
            from_pdb['ktranslrots'].append("1")

        str_data[pdbi]['FROM_PDB'] = from_pdb

    for x in sorted(local_stats):
        print("STATS", "scan_PDB", x, len(local_stats[x]), local_stats[x])
    pickle.dump(str_data, open(str_data_pkl, 'wb'))
    print("STATS", "scan_OPM", "Finished", "time", time.ctime(), "\n")
    return str_data


def check_version(options, tvec, str_data_pdbi, pdbi, tm_chains=[], thr_log_status="ERROR"):
    """Checks that 
    """
    this_name = check_version.__name__

    pdb_dict, report, prov = tvec
    if options['EXECUTION'][('debug', 'debug')]:
        print("check_version", pdbi, prov, "time", time.ctime())

    hole_thr = options['PARAMETERS'][('', 'hole_thr')]
    hole_frac_thr = options['PARAMETERS'][('', 'hole_frac_thr')]
    disordered_thr = options['PARAMETERS'][('', 'disordered_thr')]
    inter_persub_clashes_thr = options['PARAMETERS'][('', 'inter_persub_clashes_thr')]
    intra_persub_clashes_thr = options['PARAMETERS'][('', 'intra_persub_clashes_thr')]

    elim_tag, no_biopy, approval = prov + '_eliminated', prov + '_no_biopy', prov + '_approved'
    if not pdb_dict or ('COORDS' in pdb_dict and not pdb_dict['COORDS']):
        str_data_pdbi['PASSPORT'].append(passport_entry(this_name+'_no_coords', pdbi, "The coordinate file from the {0} database was not found".format(elim_tag[:3].upper())))
        if not str_data_pdbi['status'] or str_data_pdbi['status'][-1] != elim_tag:
            str_data_pdbi['status'].append(elim_tag)
        return str_data_pdbi

    # Check consistency with PDB format (if consensus_signature_PDB option is there)
    if options['PARAMETERS'][('', 'consensus_signature_PDB')] and not report["pdb_format"]:
        str_data_pdbi['PASSPORT'].append(passport_entry(this_name+'no_pdbfmt', pdbi, "Coordinate file from the {0} database cannot be translated in PDB format (common causes: >100000 atoms or multiple-character chain names). The coordinate file will not be considered in EncoMPASS"))
        if not str_data_pdbi['status'] or str_data_pdbi['status'][-1] != elim_tag:
            str_data_pdbi['status'].append(elim_tag)


    # Chain-wise checks
    for chain in sorted([x for x in pdb_dict['COORDS']]):
        # Check for holes
        if len(report['holes'][chain]) > hole_thr:
            str_data_pdbi['PASSPORT'].append(passport_entry(this_name+'_chain_hole', pdbi, "Subunit {0} of the structure as recorded in the coordinate file from the {1} database contains at least one undetermined segment longer than {2} amino acids. The coordinate file is labeled as problematic.".format(chain, prov.upper(), hole_thr)))
            if (not tm_chains) or (chain in tm_chains):
                if not str_data_pdbi['status'] or str_data_pdbi['status'][-1] != elim_tag:
                    str_data_pdbi['status'].append(elim_tag)
        # Check for holes (percent)
        if report['holes_frac'][chain] > hole_frac_thr:
            #print(pdbi, chain, report['holes'][chain])
            str_data_pdbi['PASSPORT'].append(passport_entry(this_name+'_chain_holefrc', pdbi, "Subunit {0} of the structure as recorded in the coordinate file from the {1} database contains undetermined regions amounting to more than {2}% of its sequence length. The coordinate file is labeled as problematic.".format(chain, prov.upper(), hole_frac_thr*100)))
            if (not tm_chains) or (chain in tm_chains):
                if not str_data_pdbi['status'] or str_data_pdbi['status'][-1] != elim_tag:
                    str_data_pdbi['status'].append(elim_tag)
        # Check for disordered regions
        if len(report['disordered'][chain]) > disordered_thr:
            str_data_pdbi['PASSPORT'].append(passport_entry(this_name+'_chain_disres', pdbi, 'Subunit {0} of the structure as recorded in the coordinate file from the {1} database contains disordered or repeated residue indexes. The coordinate file is labeled as problematic.'.format(chain, prov.upper())))
            if (not tm_chains) or (chain in tm_chains):
                if not str_data_pdbi['status'] or str_data_pdbi['status'][-1] != elim_tag:
                    str_data_pdbi['status'].append(elim_tag)
        # Check for unknown residues (UNK or other strange types, which were included as UNK)
        if report['UNK'][chain]:
            if options['PARAMETERS'][('', 'with_UNK')]:
                neg = ''
            else:
                neg = 'not'
            str_data_pdbi['PASSPORT'].append(passport_entry(this_name+'_UNK_resnames', pdbi, 'Subunit {0} of the structure as recorded in the coordinate file from the {1} database contains unknown residues (labeled as UNK in the original coordinate file). In compliance with the EncoMPASS guidelines, these residues will {2} be kept into account in the model'.format(chain, prov.upper(), neg)))

        # Check for chain-wise clashes
        interctot = 0
        for chain2 in sorted([x for x in pdb_dict['COORDS']]):
            if chain == chain2:
                continue
            interctot += len(report['clashes'][(chain, chain2)])
        if interctot > inter_persub_clashes_thr or len(report['clashes'][(chain, chain2)]) > intra_persub_clashes_thr:
            if interctot > inter_persub_clashes_thr:
                str_data_pdbi['PASSPORT'].append(passport_entry(this_name+'_interclash', pdbi, 'Subunit {0} of the structure as recorded in the coordinate file from the {1} database has too many steric clashes with the other subunits ({2} > {3})'.format(chain, prov.upper(), interctot, inter_persub_clashes_thr)))
            else:
                str_data_pdbi['PASSPORT'].append(passport_entry(this_name+'_interclash', pdbi, 'Subunit {0} of the structure as recorded in the coordinate file from the {1} database contains too many steric clashes ({2} > {3})'.format(chain, prov.upper(), len(report['clashes'][(chain, chain2)]), intra_persub_clashes_thr)))
            if not str_data_pdbi['status'] or str_data_pdbi['status'][-1] != elim_tag:
                str_data_pdbi['status'].append(elim_tag)

    if elim_tag not in str_data_pdbi['status'] and no_biopy not in str_data_pdbi['status']:
        str_data_pdbi['status'].append(approval)

    # Check inter-chain errors
    if report['merged_chains']:
        str_data_pdbi['PASSPORT'].append(passport_entry(this_name+'_conch', pdbi, 'The structure as recorded in the coordinate file from the {0} database contains concatenating chains. Subunit names have been modified: {1}'.format(prov.upper(), " ".join([",".join([x for x in m[1:]])+" -> "+m[0] for m in report['merged_chains']]))))

    # Check MSE->MET replacement
    if report['mse_res']:
        rep_str = ", ".join("chain {0} resID {1}{2}".format(ch, resid[0], resid[1]) for resid, ch in report['mse_res'])
        str_data_pdbi['PASSPORT'].append(passport_entry(this_name+'_mse2met', pdbi, 'The following MSE residues were replaced by MET residues: ' + rep_str))
    return str_data_pdbi


def check_whole_structures(data):
    """
    Decides whether to use the OPM (preferred) or PDB structure, or to condemn the entry
    """

    this_name = check_whole_structures.__name__

    # Decompress input
    options, locations, str_data, in_parallel, thr_log_status = data    # in_parallel is a boolean that tells if this calculation is run in parallel

    # Options and locations aliases
    LOC_pad = locations['SYSFILES']['pdbatomdict']

    tic = time.time()
    new_OPM_data = {}

    for pdbi in sorted([k for k in str_data]):
        # Parses PDB as it is downloaded from the PDB
        # with scan == True: clash control, hole control

        # Filter
        if 'eliminated' in str_data[pdbi]['status']:
            continue

        str_data[pdbi]['FROM_PDB']['cache_coordinate_file'] = locations['FSYSPATH']['TMPcoords'] + pdbi+'_PDB.cif'
        str_data[pdbi]['FROM_OPM']['cache_coordinate_file'] = locations['FSYSPATH']['TMPcoords'] + pdbi+'_OPM.cif'
        PDB_dict, PDB_report = parse_mmCIF_coords_wrap(str_data[pdbi]['FROM_PDB']['coordinate_file'], scan=True, options=options)
        if PDB_dict and 'COORDS' in PDB_dict and PDB_dict['COORDS']:
            write_mmCIF(str_data[pdbi], PDB_dict, str_data[pdbi]['FROM_PDB']['cache_coordinate_file'], LOC_pad)

        OPMpdb_dict, OPM_report = parse_PDB_coords(str_data[pdbi]['FROM_OPM']['coordinate_file'], scan=True, options=options)	# Can be null!
        if OPMpdb_dict and 'COORDS' in OPMpdb_dict and OPMpdb_dict['COORDS']:
            write_mmCIF(str_data[pdbi], OPMpdb_dict, str_data[pdbi]['FROM_OPM']['cache_coordinate_file'], LOC_pad)

        # Checks on OPM and PDB structure (trust OPM more than PDB)
        for tvec in [(OPMpdb_dict, OPM_report, 'opm'), (PDB_dict, PDB_report, 'pdb')]:
            tmch = str_data[pdbi]['FROM_OPM']['ksubunits'] if tvec[2] == 'opm' else []
            str_data[pdbi] = check_version(options, tvec, str_data[pdbi], pdbi, tm_chains=tmch)

        # Check OPM-PDB elimination from previous checks
        if 'pdb_eliminated' in str_data[pdbi]['status'] and 'opm_eliminated' in str_data[pdbi]['status']:
            if not str_data[pdbi]['status'] or str_data[pdbi]['status'][-1] != 'eliminated':
                str_data[pdbi]['status'].append('eliminated')
                str_data[pdbi]['delete_keyword'] = 'Unclear'
            str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_nocoord', pdbi, 'There are no reliable coordinate files for this structure. The entry is scheduled for deletion.'))
        pdbswitch = False	# Flag for changing the reference to PDB. True when OPM was deemed incorrect or dangerous in the following post-processing
        # In case OPM is not eliminated
        if 'eliminated' not in str_data[pdbi]['status'] and 'opm_eliminated' not in str_data[pdbi]['status']:
            # Decides whether to switch to PDB or to keep the OPM structure
            pdbswitch = False
            for sub_ch in str_data[pdbi]['FROM_OPM']['ksubunits']:
                if not OPMpdb_dict or sub_ch not in OPMpdb_dict['COORDS']:
                    if not OPMpdb_dict:
                        str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_noopmcoord', pdbi, 'There is no OPM coordinate files for this structure.'))
                    else:
                        str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_noopmch', pdbi, 'In the OPM coordinate file, chain {0} described on the corresponding OPM web page is missing. The OPM coordinate file will not be considered anymore.'.format(sub_ch)))
                    pdbswitch = True
                    str_data[pdbi]['status'].append('opm_eliminated')
                    break
            # Switch to PDB
            if pdbswitch:
                # Check whether the switch is possible, or whether the structure must be eliminated
                if 'pdb_eliminated' in str_data[pdbi]['status']:
                    if not str_data[pdbi]['status'] or str_data[pdbi]['status'][-1] != 'eliminated':
                        str_data[pdbi]['status'].append('eliminated')
                        str_data[pdbi]['delete_keyword'] = 'Unclear'
                    str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_nocoord2', pdbi, 'There are no reliable coordinate files for this structure. The entry is scheduled for deletion.'))
                    pdbswitch = False
            else:
                # Maintain OPM
                str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_coord_opm', pdbi, 'Coordinate file taken from OPM'))
#                str_data[pdbi]['ENCOMPASS']['coordinate_file'] = str_data[pdbi]['FROM_OPM']['coordinate_file']
                str_data[pdbi]['present_cache'] = str_data[pdbi]['FROM_OPM']['cache_coordinate_file']

        # Check PDB-only condemnation from previous checks
        if 'eliminated' not in str_data[pdbi]['status'] and 'pdb_eliminated' not in str_data[pdbi]['status']:
            bio_asb = str_data[pdbi]['FROM_PDB']['biological_assemblies']
            # If biological assembly is present, check consistency between biological assemblies and coordinates
            if len(bio_asb) > 0:
                # If a chain with coordinates is not in biological assembly, the PDB becomes monitored
                for ch in PDB_dict['COORDS']:
                    notin = True
                    for transf in bio_asb[0]['transformations']:
                        if ch in transf['chains']:
                            notin = False
                            break
                    if notin:
                        str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_chain_noba', pdbi, 'Subunit {0} of the structure as recorded in the coordinate file from the PDB database is not mentioned in the chosen set of transformation matrices for generating the biological assembly. This structure is classified as very problematic, and will be discarded should a solid record from another database (OPM, PDBTM) not be found.'.format(ch)))
                        if not str_data[pdbi]['status'] or str_data[pdbi]['status'][-1] != 'pdb_eliminated':
                            str_data[pdbi]['status'].append('pdb_eliminated')
                # If a chain in biological assembly is not found to have coordinates, the PDB becomes monitored
                for transf in bio_asb[0]['transformations']:
                    for ch in transf['chains']:
                        if ch not in PDB_dict['COORDS']:
                            str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_chain_banoc', pdbi, 'Subunit {0} of the structure as recorded in the coordinate file from the PDB database is mentioned in the chosen set of transformation matrices for generating the biological assembly, yet its coordinates are not present. This structure is classified as very problematic, and will be discarded should a solid record from another database (OPM, PDBTM) not be found.'.format(ch)))
                            if not str_data[pdbi]['status'] or str_data[pdbi]['status'][-1] != 'pdb_eliminated':
                                str_data[pdbi]['status'].append('pdb_eliminated')
            # If biological assembly is not present, check if there is one or more chains. 
            else:
                if not str_data[pdbi]['status'] or str_data[pdbi]['status'][-1] != 'pdb_monitored':
                    str_data[pdbi]['status'].append('pdb_monitored')
                if len(PDB_dict['COORDS']) == 1:
                    str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_1ch_noba', pdbi, 'The structure as recorded in the coordinate file from the PDB database contains a single subunit and does not have a biological assembly record'))
                else:
                    str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_2pch_noba', pdbi, 'The structure as recorded in the coordinate file from the PDB database contains more than one subunit and does not have a biological assembly record'))

            # If pdbswitch is present, change to PDB, else maintain the OPM structure
            if pdbswitch:
                str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_coord_pdb', pdbi, 'Basing on the evidences of the previous checks, the coordinate file taken to represent this EncoMPASS entry is taken from the PDB database'))
#                str_data[pdbi]['ENCOMPASS']['coordinate_file'] = str_data[pdbi]['FROM_PDB']['coordinate_file'] # WARNING HERE. WHY WAS THIS COMMENTED?
                str_data[pdbi]['present_cache'] = str_data[pdbi]['FROM_PDB']['cache_coordinate_file']

    toc = time.time()
    return str_data, [pdbi for pdbi in str_data]


def parallel_homogeneous_parsing(options, locations, str_data, thr_log_status="ERROR"):
    '''
    For each entry of the main dictionary, it parses the corresponding coordinate file(s) and tries to correct format and structure errors
    It is parallelized (active if the option 'number_of_procs' is found)
    '''

    # Call to check_whole_structures

    if options['RUN'][('np', 'number_of_processors')] > 1:
        # Step 1: prepare independent inputs
        data = []
        for pdbi in sorted(str_data):
            if ('eliminated' in str_data[pdbi]['status'] or
                'forced_pdb' in str_data[pdbi]['status'] or
                'forced_opm' in str_data[pdbi]['status']):
                continue 
            data.append((options, locations, {pdbi : str_data[pdbi]}, True, thr_log_status))


        # Step 2: run multiprocess
        np = int(options['RUN'][('np', 'number_of_processors')])
        pool_outputs = []
        if options['EXECUTION'][('debug', 'debug')] and not options['EXECUTION'][('paraldebug', 'parallel_debug')]:
            print("check_whole_structures (serial version)", "time", time.ctime())
            for d in data:
                print(",".join([x for x in d[2]]))
                o = check_whole_structures(d)
                pool_outputs.append(o)
        else:
            if options['EXECUTION'][('debug', 'debug')]:
                print("check_whole_structures (parallel version)", "time", time.ctime())
            pool = multiprocessing.Pool(processes=np)
            pool_outputs = pool.map(check_whole_structures, data)
            pool.close()
            pool.join

        # Step 3: collect outputs
        allrPPM = set()
        for s in pool_outputs:
           if not s:
               continue
           str_dict, pdbis = s
           pdbi = pdbis[0]
           if pdbi in str_dict:
               str_data[pdbi] = copy.deepcopy(str_dict[pdbi])
    else:
        str_data, tmp = check_whole_structures((options, locations, str_data, False, thr_log_status))

    print("check_whole_structures", "Finished", "time", time.ctime(), "\n")
    return str_data


def combine_sources(options, locations, str_data):
    """Combines OPM and PDB
    """
    this_name = combine_sources.__name__

    # Options and locations aliases
    LOC_pad = locations['SYSFILES']['pdbatomdict']

    accepted_techniques = ['SOLUTION NMR', 'SOLID-STATE NMR', 'ELECTRON MICROSCOPY', 'ELECTRON CRYSTALLOGRAPHY', 'THEORETICAL MODEL', 'FIBER DIFFRACTION', 'X-RAY DIFFRACTION']

    local_stats = {
        'rescued_alphabeta' : [],
        'no_alphabeta_2' : [], # This number can be less than no_alphabeta - rescued_alohabeta because some entries could have been eliminated
        'resolution_filtered' : [],
        'unclass_experiment_filtered' : [],
        'nmr_experiment_filtered' : [],
        'cryoem_experiment_filtered' : [],
        'theo_experiment_filtered' : [],
        'hand_filtered' : [],
        'hand_db' : [],
        'eliminated' : [],
        'opm_eliminated' : [],
        'pdb_eliminated' : []
    }

    deletion_fn = options['EXECUTION'][('del', 'deletion_entries')]
    # ASSUMES format: tsv 2 columns. 1: PDBID, 2: comment
    del_d = {}
    if os.path.exists(deletion_fn):
        with open(deletion_fn) as delf:
            for line in delf:
                tlist = [x.strip() for x in line.split("\t")]
                if len(tlist) != 3:
                    print(tlist)
                    print_log((
                        'CRITICAL',
                        this_name,
                        "File {0} badly formatted. Line:\n{1}"\
                            .format(deletion_fn, line)
                    ))
                pdbi, keyword, comment = [x.strip() for x in tlist]
                del_d[pdbi] = (keyword, comment)

    replace_fn = options['EXECUTION'][('repl', 'replace_entries')]
    # ASSUMES format: tsv 2 columns. 1: PDBID, 2: path
    
    forcedb_fn = options['EXECUTION'][('forcedb', 'force_db_choice')]
    # ASSUMES format: tsv 2 columns. 1: PDBID, 2: OPM or PDB
    forcedb_d = {}
    if os.path.exists(forcedb_fn):
        with open(forcedb_fn) as delf:
            for line in delf:
                tlist = [x.strip() for x in line.split("\t")]
                if len(tlist) != 3:
                    print(tlist)
                    print_log((
                        'CRITICAL',
                        this_name,
                        "File {0} badly formatted. Line:\n{1}"\
                            .format(deletion_fn, line)
                    ))
                pdbi, db, motivation = [x.strip() for x in tlist]
                if db in ['OPM', 'PDB']:
                    forcedb_d[pdbi] = (db, motivation)
                else:
                    print_log((
                        'ERROR',
                        this_name,
                        ("PDB code {0} was forced to take info from "
                             "non-existing database {1}")\
                             .format(pdbi, db)
                    ))


    print("combine_sources", "combine inputs", "time", time.ctime())
    for pdbi in sorted(str_data):
        if options['EXECUTION'][('debug', 'debug')]:
            print(pdbi)

        # Name
        if str_data[pdbi]['FROM_OPM']['name']:
            str_data[pdbi]['ENCOMPASS']['name'] = str_data[pdbi]['FROM_OPM']['name']
        elif str_data[pdbi]['FROM_PDB']['name']:
            str_data[pdbi]['ENCOMPASS']['name'] = str_data[pdbi]['FROM_PDB']['name']
        elif str_data[pdbi]['FROM_PDB']['title']:
            str_data[pdbi]['ENCOMPASS']['name'] = str_data[pdbi]['FROM_PDB']['title']

        str_data[pdbi]['ENCOMPASS']['coordinate_file'] = locations['FSYSPATH']['ENCpdbs'] + pdbi + '_enc.pdb'

        # Class (already assigned in from_OPM, reviewed here)
        if str_data[pdbi]['ENCOMPASS']['class'] == 'other':
            if str_data[pdbi]['FROM_PDB']['ss_class'] == 'all-alpha':
                str_data[pdbi]['ENCOMPASS']['class'] = 'alpha'
                local_stats['rescued_alphabeta'].append(pdbi)
            if str_data[pdbi]['FROM_PDB']['ss_class'] == 'all-beta':
                str_data[pdbi]['ENCOMPASS']['class'] = 'beta'
                local_stats['rescued_alphabeta'].append(pdbi)
        if str_data[pdbi]['ENCOMPASS']['class'] == 'other':
            local_stats['no_alphabeta_2'].append(pdbi)

        # Resolution
        if str_data[pdbi]['FROM_OPM']['resolution']:
            str_data[pdbi]['ENCOMPASS']['resolution'] = str_data[pdbi]['FROM_PDB']['resolution']
            str_data[pdbi]['ENCOMPASS']['structure']['quality']['metric'] = 'resolution'
            str_data[pdbi]['ENCOMPASS']['structure']['quality']['value'] = str_data[pdbi]['FROM_PDB']['resolution']
        else:
            str_data[pdbi]['ENCOMPASS']['resolution'] = 'N/A'	# Not Applicable

        # Resolution FILTERS
        if options['PARAMETERS'][('res_thr', 'structure_resolution_threshold')]:
            if str_data[pdbi]['ENCOMPASS']['resolution'] == 'N/A':
                if not str_data[pdbi]['status'] or str_data[pdbi]['status'][-1] != 'eliminated':
                    str_data[pdbi]['status'].append('eliminated')
                    str_data[pdbi]['delete_keyword'] = 'Filtered'
                str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_nores', pdbi, "The coordinate file from the PDB database does not contain any information regarding the resolution of the structure, and therefore does not comply with the EncoMPASS guidelines (filter 'resolution <= {0}'). This entry is scheduled for deletion.".format(options['PARAMETERS'][('res_thr', 'structure_resolution_threshold')])))
                local_stats['resolution_filtered'].append(pdbi)
            elif float(str_data[pdbi]['ENCOMPASS']['resolution']) > float(options['PARAMETERS'][('res_thr', 'structure_resolution_threshold')]):
                if not str_data[pdbi]['status'] or str_data[pdbi]['status'][-1] != 'eliminated':
                    str_data[pdbi]['status'].append('eliminated')
                    str_data[pdbi]['delete_keyword'] = 'Filtered'
                str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_nores', pdbi, "The structure as recorded in the coordinate file from the PDB database has a resolution of {0}, and therefore does not comply with the EncoMPASS guidelines (filter 'resolution <= {1})'. This entry is scheduled for deletion.".format(str_data[pdbi]['ENCOMPASS']['resolution'], options['PARAMETERS'][('res_thr', 'structure_resolution_threshold')])))
                local_stats['resolution_filtered'].append(pdbi)

        # Experiment
        if str_data[pdbi]['FROM_PDB']['experimental_method']:
            str_data[pdbi]['ENCOMPASS']['structure']['experiment'] = str_data[pdbi]['FROM_PDB']['experimental_method']

        # Experiment FILTERS
        #  Unrecognized experimental techniques
        if ";" in str_data[pdbi]['ENCOMPASS']['structure']['experiment']:
            str_data[pdbi]['ENCOMPASS']['structure']['experiment'] = str_data[pdbi]['ENCOMPASS']['structure']['experiment'][:str_data[pdbi]['ENCOMPASS']['structure']['experiment'].index(";")]
        if str_data[pdbi]['ENCOMPASS']['structure']['experiment'] and str_data[pdbi]['ENCOMPASS']['structure']['experiment'] not in accepted_techniques and not [1 for x in accepted_techniques if x in str_data[pdbi]['ENCOMPASS']['structure']['experiment']]:
            if not str_data[pdbi]['status'] or str_data[pdbi]['status'][-1] != 'eliminated':
                str_data[pdbi]['status'].append('eliminated')
                str_data[pdbi]['delete_keyword'] = 'Filtered'
            str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_noexpdta', pdbi, "This structure as recorded in the coordinate file from the PDB database has been determined with the unclassified experimental technique: {0}. This entry is scheduled for deletion".format(str_data[pdbi]['ENCOMPASS']['structure']['experiment'])))
            local_stats['unclass_experiment_filtered'].append(pdbi)
            continue

        #  Note: there would be also FIBER DIFFRACTION, but a filter for that is not implemented yet (resolution is defined)
        if options['PARAMETERS'][('no_nmr', 'no_nmr')] and 'NMR' in str_data[pdbi]['ENCOMPASS']['structure']['experiment']:	# SOLUTION NMR or SOLID-STATE NMR
            if not str_data[pdbi]['status'] or str_data[pdbi]['status'][-1] != 'eliminated':
                str_data[pdbi]['status'].append('eliminated')
                str_data[pdbi]['delete_keyword'] = 'Filtered'
            str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_nmrout', pdbi, "This structure as recorded in the coordinate file from the PDB database has been determined with the experimental technique: NMR. Such technique is not compliant with the EncoMPASS guidelines (filter 'no NMR structures'). This entry is scheduled for deletion."))
            local_stats['nmr_experiment_filtered'].append(pdbi)
            continue
        elif options['PARAMETERS'][('no_cryoEM', 'no_cryoEM')] and 'ELECTRON' in str_data[pdbi]['ENCOMPASS']['structure']['experiment']:	# ELECTRON MICROSCOPY or ELECTRON CRYSTALLOGRAPHY
            if not str_data[pdbi]['status'] or str_data[pdbi]['status'][-1] != 'eliminated':
                str_data[pdbi]['status'].append('eliminated')
                str_data[pdbi]['delete_keyword'] = 'Filtered'
            str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_cryoout', pdbi, "This structure as recorded in the coordinate file from the PDB database has been determined with the experimental technique: Cryo-EM. Such technique is not compliant with the EncoMPASS guidelines (filter 'no cryo-EM structures'). This entry is scheduled for deletion."))
            local_stats['cryoem_experiment_filtered'].append(pdbi)
            continue
        elif options['PARAMETERS'][('no_theoretical', 'no_theoretical')] and 'THEORETICAL' in str_data[pdbi]['ENCOMPASS']['structure']['experiment']:	# THEORETICAL MODEL
            if not str_data[pdbi]['status'] or str_data[pdbi]['status'][-1] != 'eliminated':
                str_data[pdbi]['status'].append('eliminated')
                str_data[pdbi]['delete_keyword'] = 'Filtered'
            str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_theoout', pdbi, "This structure as recorded in the coordinate file from the PDB database is a theoretical structure and therefore does not comply with the EncoMPASS guidelines (filter 'no theoretical structures'). This entry is scheduled for deletion."))
            local_stats['theo_experiment_filtered'].append(pdbi)
            continue

        # Other quality metrics (in case not resolution)
        # ...

        if pdbi in forcedb_d:
            str_data[pdbi]['status'].append('forcedb_'+forcedb_d[pdbi][0].lower())
            str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_forced', pdbi, "This entry has been manually set to take structural information from the {0} databse. The reason is: {1}".format(forcedb_d[pdbi][0], forcedb_d[pdbi][1])))
            if forcedb_d[pdbi][0] == "OPM":
                str_data[pdbi]['status'].append('pdb_eliminated')
            else:
                str_data[pdbi]['status'].append('opm_eliminated')
            str_data[pdbi]['present_cache'] = str_data[pdbi]['FROM_'+forcedb_d[pdbi][0].upper()]['cache_coordinate_file']
            local_stats['hand_db'].append((pdbi, forcedb_d[pdbi][0]))
        if pdbi in del_d:
            str_data[pdbi]['status'].append('eliminated')
            str_data[pdbi]['PASSPORT'].append(passport_entry(this_name+'_delled', pdbi, "This entry has been manually removed from the database. The reason is: {0}".format(del_d[pdbi][1])))
            str_data[pdbi]['delete_keyword'] = del_d[pdbi][0]
            local_stats['hand_filtered'].append(pdbi)

    str_data = parallel_homogeneous_parsing(options, locations, str_data, thr_log_status="ERROR")

    print("combine_sources", "write ENC structures", "time", time.ctime())
    for pdbi in sorted(str_data):
        if options['EXECUTION'][('debug', 'debug')]:
            print(pdbi)
        tmpOPMcif = locations['FSYSPATH']['TMPcoords'] + pdbi + '_OPM.cif'
        tmpPDBcif = locations['FSYSPATH']['TMPcoords'] + pdbi + '_PDB.cif'
        if os.path.exists(tmpOPMcif):
            c_dict_pdbi = retrieve_coords_from_CIF(tmpOPMcif)
        elif os.path.exists(tmpPDBcif):
            c_dict_pdbi = retrieve_coords_from_CIF(tmpPDBcif)
        else:
            c_dict_pdbi = {'COORDS' : {}, 'DUMMIES' : []}

        str_data_pdbi = str_data[pdbi]
        if 'eliminated' in str_data[pdbi]['status']:
            local_stats['eliminated'].append(pdbi)
            write_ENC(str_data_pdbi, c_dict_pdbi, locations['FSYSPATH']['ENCpdbs'] + pdbi + '_enc.pdb', locations['SYSFILES']['ignoredpasscodes'], LOC_pad)
        else:
            write_ENC(str_data_pdbi, c_dict_pdbi, locations['FSYSPATH']['ENCpdbs'] + pdbi + '_enc.tmp.pdb', locations['SYSFILES']['ignoredpasscodes'], LOC_pad)

        if 'pdb_eliminated' in str_data[pdbi]['status']:
            local_stats['pdb_eliminated'].append(pdbi)
        elif 'opm_eliminated' in str_data[pdbi]['status']:
            local_stats['opm_eliminated'].append(pdbi)

    for x in sorted(local_stats):
        print("STATS", "combine_sources", x, len(local_stats[x]), local_stats[x])

    write_checkpoint(str_data, locations['FSYSPATH']['cache'] + 'str_data_1.pkl')
    print("combine_sources", "Finished", "time", time.ctime(), "\n")
    return str_data


def write_ENC_structures(locations, str_data):

    # Locations aliases
    LOC_pad = locations['SYSFILES']['pdbatomdict']
    LOC_TMPcoords = locations['FSYSPATH']['TMPcoords']
    LOC_ENCpdbs = locations['FSYSPATH']['ENCpdbs']

    for pdbi in str_data:
        str_data_pdbi = str_data[pdbi]
        tmpOPMcif = LOC_TMPcoords + pdbi + '_OPM.cif'
        tmpPDBcif = LOC_TMPcoords + pdbi + '_PDB.cif'
        if os.path.exists(tmpOPMcif) and 'opm_eliminated' not in str_data_pdbi['status']:
            c_dict_pdbi = retrieve_coords_from_CIF(tmpOPMcif)
        elif os.path.exists(tmpPDBcif) and 'pdb_eliminated' not in str_data_pdbi['status']:
            c_dict_pdbi = retrieve_coords_from_CIF(tmpPDBcif)
        else:
            c_dict_pdbi = {'COORDS' : {}, 'DUMMIES' : []}
        write_ENC(str_data_pdbi, c_dict_pdbi, LOC_ENCpdbs + pdbi + '_enc.pdb', locations['SYSFILES']['ignoredpasscodes'], LOC_pad)


if __name__ == '__main__':
    import initialize_repository
    options, locations = initialize_repository.initialize_repository()
    ids, secd = mainlist_from_OPM(options, locations)
    str_data = scan_OPM(options, locations, ids)
    str_data = scan_PDB(options, locations, str_data=str_data)
    str_data = combine_sources(options, locations, str_data)
