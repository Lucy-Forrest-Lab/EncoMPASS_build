# Name: parse_source_databases.py
# Language: python3
# Libraries:
# Description:
# Author: Edoardo Sarti, Lucy Forrest
# Date: 2025-02-23

from encompass.pipeline.supporting_functions import *


def uniprot_sifts(locations, str_data):
	# SIFTS file must be done like this:
	# PDB CHAIN ACC BEGINRES_NUM ENDRES_NUM BEGINRES_PDBIDX ENDRES_PDBIDX BEGINRES_UNIPROT ENDRES_UNIPROT
	uniprot_data = {}
	uniprot_set = set()

	if not os.path.exists(locations['SYSFILES']['sifts_uniprot']):
		#url = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz"
		url = locations['EXTERNAL']['uniprot_url']
		local_filename = locations['SYSFILES']['sifts_uniprot']
		iterate_download(url, local_filename, True) # zipped

	with open(locations['SYSFILES']['sifts_uniprot']) as f:
		for line in f:
			if line.strip().startswith("#"):
				continue
			pdbi, ch, acc, begnum, endnum, begpdb, endpdb, begp, endp = line.split()
			if pdbi in str_data:
				if pdbi not in uniprot_data:
					uniprot_data[pdbi] = []
				if acc not in uniprot_data[pdbi]:
					uniprot_data[pdbi].append((acc, int(begp)-1, int(endp)-1))
				uniprot_set.add(acc)
	return uniprot_data, uniprot_set


def mainlist_from_OPM(options, locations, str_data={},
					  use_pkl=False,  # Deprecated?
					  thr_log_status='ERROR'):
	"""Returns the list of primary and secondary OPM entries to be considered
	during this instance of EncoMPASS
	"""

	this_name = mainlist_from_OPM.__name__
	print("Get main list from OPM", this_name)

	# Options and locations aliases
	OPT_from_repo = options['RUN'][('repo', 'from_repo')]
	OPT_offline = options['RUN'][('ol', 'offline')]
	OPT_sample_structures = options['EXECUTION'][('sample', 'sample_structures')]
	OPT_debug = options['EXECUTION'][('debug', 'debug')]
	LOC_OPMjsons = locations['FSYSPATH']['OPMjsons']
	LOC_OPMreprchart = locations['SYSFILES']['OPMreprchart']

	# Statistics
	local_stats = {
		'primary': [],
		'secondary': [],
		'associated_secondary': [],
		'selected_primary': [],
		'selected_secondary': [],
		'temp_primary': [],
		'sample_not_found': []
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
            locations['EXTERNAL'][opm_primary_url],
            #"https://lomize-group-opm.herokuapp.com/types/1/primary_structures?pageSize=100000",
			params={},
			headers={"Accept": "application/json"}
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
	ids = [(x['id'], x['pdbid'], False) for x in data['objects']]  # Primary (<int>ID, <str>PDBID, <bool>is_temp) list
	ids_d = {x['id']: x['pdbid'] for x in data['objects']}  # Primary <int>ID-><str>PDBID dict (internal use)
	revids_d = {x['pdbid']: x['id'] for x in data['objects']}  # Primary <str>PDBID-><int>ID dict (internal use)

	# Download/retrieve secondary structure json file
	secjson_filename = LOC_OPMjsons + 'secondary.json'
	repofilename = OPT_from_repo + LOC_OPMjsons + 'secondary.json'
	if OPT_from_repo and os.path.exists(repofilename):
		shutil.copyfile(repofilename, secjson_filename)
	else:
		if not OPT_offline:
			res = requests.get(
                locations['EXTERNAL'][opm_secondary_url],
				#"https://lomize-group-opm.herokuapp.com/secondary_representations?pageSize=1000000",
				params={},
				headers={"Accept": "application/json"}
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
	secrep = []  # Secondary <str>PDBID list (internal use)
	secrep_d = {}  # Secondary <str>PDBID->[<str>PDBID, ...] dict
	revsecrep_d = {}  # Secondary <str>PDBID-><int>ID dict (internal use)
	for x in data['table']['objects']:
		if x['primary_structure_id'] not in ids_d:
			# print("CONT", x['primary_structure_id'])
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
		new_ids = []  # Primary (<int>ID, <str>PDBID, <bool>is_temp) list
		new_secrep_d = {}  # Secondary <int>ID->[<str>PDBID, ...] dict
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
					if entry[:4] in secrep_d:  # entry is absent if it does not have secondary reprs
						new_secrep_d[entry[:4]] = secrep_d[entry[:4]]
					else:
						print_log((
							'WARNING',
							this_name,
							("Entry {0} was included in sample_structures, yet "
							 "{1} does not have secondary representations") \
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
	print("Start scanning OPM", this_name)

	# Options and locations aliases
	OPT_dst = options['PATHS'][('dst', 'data_structure_template')]
	OPT_from_repo = options['RUN'][('repo', 'from_repo')]
	OPT_offline = options['RUN'][('ol', 'offline')]
	OPT_debug = options['EXECUTION'][('debug', 'debug')]
	OPT_sample_structures = \
		options['EXECUTION'][('sample', 'sample_structures')]
	LOC_OPM = locations['FSYSPATH']['OPM']
	LOC_OPMjsons = locations['FSYSPATH']['OPMjsons']
	LOC_OPMpdbs = locations['FSYSPATH']['OPMpdbs']
	LOC_UniProt = locations['FSYS']['UniProt']
	LOC_uniprot_all = locations['SYSFILES']['uniprot_all']
	LOC_OPMreprchart = locations['SYSFILES']['OPMreprchart']

	# Statistics
	local_stats = {
		'json_not_found': [],
		'no_segments': [],
		'no_alphabeta': [],
		'OPMpdb_not_found': []
	}

	# Update repository
	update_repo = check_update_repo(options, locations)

	# Download/retrieval of OPM infos about primary entries
	print("scan_OPM", "Download json files", "time", time.ctime())
	primary_structures_data = {}
	for i, pdbi, is_temp in sorted(entries, key=lambda x: x[1]):
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
					f"{locations['EXTERNAL']['opm_primary_root']}/{i}",
					#f"https://lomize-group-opm.herokuapp.com/primary_structures/{i}",
					params={},
					headers={"Accept": "application/json"}
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
			with codecs.open(json_filename, 'r+', encoding='utf-8') \
					as json_file:
				primary_structures_data[(i, pdbi, is_temp)] = json.load(json_file,
																		encoding='utf-8')

	# Download OPM coordinate files
	print("scan_OPM", "Download coordinate files", "time", time.ctime())
	for i, pdbi, is_temp in sorted(primary_structures_data, key=lambda x: x[1]):
		if options['EXECUTION'][('debug', 'debug')]:
			print(pdbi)
		if pdbi not in str_data:
			str_data[pdbi] = define_str_data_entry(OPT_dst)

		# Update status
		str_data[pdbi]['status'].append('initial_opm')

		# EncoMPASS name is set to OPM name
		str_data[pdbi]['ENCOMPASS']['name'] = \
			primary_structures_data[(i, pdbi, is_temp)]['name']
		str_data[pdbi]['PASSPORT'].append(passport_entry(
			this_name + '_encname_opm',
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
				download_result = iterate_download(
					f"{locations['EXTERNAL']['opm_coord_root']}/{pdbi}.pdb",
					#f'https://opm-assets.storage.googleapis.com/pdb/{pdbi}.pdb',
					pdb_path
				)

		# Remove 0-segment structures
		if int(primary_structures_data[(i, pdbi, is_temp)]['subunit_segments']) == 0:
			str_data[pdbi]['PASSPORT'].append(passport_entry(
				this_name + '_zero_tmsegments',
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
					str_data[pdbi]['FROM_OPM']['superfamily'][y] = \
						primary_structures_data[(i, pdbi, is_temp)]['family']['superfamily'][y]
			# Class renaming: alpha-beta-other
			# TO DO add a statistic here
			elif 'class' in x:
				classtype = primary_structures_data[(i, pdbi, is_temp)]['family'] \
					['superfamily']['classtype']['name']
				str_data[pdbi]['FROM_OPM']['classtype'] = \
					primary_structures_data[(i, pdbi, is_temp)]['family'] \
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
					this_name + '_encclass_opm',
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
					str_data[pdbi]['FROM_OPM'][x][y] = \
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
							str_data[pdbi]['FROM_OPM']['subunits'] \
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
								# print("OPM SEGMENTS", pdbi, subunits_key, tempd[y])
							else:
								tempd[y] = subunit[y]
						str_data[pdbi]['FROM_OPM']['subunits'] \
							.append(tempd)
						# print("WRITE", tempd.show_dictionary())
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
					str_data[pdbi]['FROM_OPM'][x] = \
						primary_structures_data[(i, pdbi, is_temp)][x]
		# Check the primary representation case
		str_data[pdbi]['FROM_OPM']['primary_representation'] = True

		# Secondary representations
		if len(str_data[pdbi]['FROM_OPM']['secondary_representations']) > 0:
			# Loop on the codes found in the OPM json page
			all_sec_pdbi = [x['pdbid'] for x in \
							primary_structures_data[(i, pdbi, is_temp)]['secondary_representations']]
			for sec_pdbi in all_sec_pdbi:
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
						download_result = iterate_download(
							f"{locations['EXTERNAL']['opm_coord_root']}/{sec_pdbi}.pdb",
							#f'https://opm-assets.storage.googleapis.com/pdb/{sec_pdbi}.pdb',
							sec_pdb_path
						)
				# Copy the str_data entry of the pdbi entry into sec_pdbi -> DANGEROUS: CHECK ALL PATHS
				str_data[sec_pdbi] = FixedDict(str_data[pdbi])

				if str_data[pdbi]['ENCOMPASS']['class'] == 'other':
					local_stats['no_alphabeta'].append(sec_pdbi)

				# OPM coordinate file analysis
				if os.path.exists(sec_pdb_path):
					#                    str_data[sec_pdbi]['FROM_OPM'] =\
					#                        FixedDict(str_data[pdbi]['FROM_OPM'])   # Copy everything, yes! But then adapt POV
					str_data[sec_pdbi]['FROM_OPM'] \
						['coordinate_file'] = sec_pdb_path
					str_data[sec_pdbi]['FROM_OPM'] \
						['primary_representation'] = False

					# Add as secondary representations all other secondary
					# representations of its primary representation except
					# itself, plus the primary representation
					str_data[sec_pdbi]['FROM_OPM']['secondary_representations'] = \
						str_data[sec_pdbi]['FROM_OPM'] \
							.show_dictionary(quiet=True)['secondary_representations']
					for secondary_entry, xx in enumerate(str_data[sec_pdbi]['FROM_OPM'] \
																 ['secondary_representations']):
						if xx['pdbid'] == sec_pdbi:
							sec_to_del = secondary_entry
							break
					del (str_data[sec_pdbi]['FROM_OPM']['secondary_representations'][sec_to_del])
					# print("deleting secondary rep", sec_to_del)

					elm = str_data[sec_pdbi]['FROM_OPM'] \
						['secondary_representations'].get_fdict_template()
					for ii, dd, tt in entries:
						if dd == pdbi:
							elm['id'] = ii
							break
					elm['pdbid'] = pdbi
					str_data[sec_pdbi]['FROM_OPM'] \
						['secondary_representations'].append(FixedDict(elm))
				else:
					local_stats['OPMpdb_not_found'].append(sec_pdbi)
					print_log((
						'ERROR',
						this_name,
						"OPM coordinate file {0} was not found in {1}" \
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
			del (str_data[pdbi])
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
	#str_data_pkl = LOC_OPM + '.str_data_after_OPM_read.pkl'
	#pickle.dump(str_data, open(str_data_pkl, 'wb'))
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
	print("Start scanning PDB", this_name)

	# Options and locations aliases
	local_stats = {
		'no_record': [],
		'no_mmcif': [],
		'empty_mmcif': [],
		'corrupted_mmcif': [],
		'no_pdb': [],
		'empty_pdb': [],
		'corrupted_pdb': [],
		'no_str': [],
		'no_title': [],
		'no_expmethod': [],
		'no_resolution': [],
		'no_depo_date': [],
		'no_entities': []
	}

	if not pdb_list:
		pdb_list = sorted([x for x in str_data])

	# move management of pickles to wrapper script
	# Pickle
	#str_data_pkl = locations['FSYSPATH']['PDB'] + '.str_data_after_PDB_read.pkl'
	#if use_pkl and os.path.exists(str_data_pkl):
	#	str_data = pickle.load(open(str_data_pkl, 'rb'))

	# Obsolete PDB entries
	obsolete = set()
	obsolete_filename = locations['FSYSPATH']['PDBjsons'] + 'obsolete_entries.json'

	#  Get them from repo or server
	repofilename = options['RUN'][('repo', 'from_repo')] + locations['FSYS']['PDBjsons'] + 'obsolete_entries.json'
	if options['RUN'][('repo', 'from_repo')] and os.path.exists(repofilename):
		shutil.copyfile(repofilename, obsolete_filename)
	else:
		# TO DO: find obsolete
		pass

	#  TO DO: make a loop for filtering obsolete entries ('eliminated')

	# Json info records
	temparch = {}  # Keys 'coordinate_file' and 'header' do not come from json but from structure

	# pdbis -> batches
	batch_size = 100
	batches = []
	with open(locations['FSYSPATH']['cache'] + 'PDBall.txt', 'w') as f:
		for ipdbi, pdbi in enumerate(pdb_list):
			if ipdbi % batch_size == 0:
				batches.append([])
			batches[-1].append(pdbi.upper())
			f.write(pdbi.upper() + '\n')

	# query by batches
	url = 'https://data.rcsb.org/graphql?'
	headers = []
	print("scan_PDB", "Querying PDB by batches", "time", time.ctime())
	for ibatch, batch in enumerate(batches):
		if options['EXECUTION'][('debug', 'debug')]:
			print("Batch", ibatch)
		# Batch query
		json_query = pdb_query(path=locations['SYSFILES']['pdbquery'])
		json_query = json_query.replace("<PDBI>", ", ".join(['"' + x + '"' for x in batch]))
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
					str_data[pdbi]['PASSPORT'].append(
						passport_entry(this_name + '_no_pdbmirror', pdbi, "This entry was not found in the pdb mirror"))
			elif options['RUN'][('pdbmirror', 'pdb_local_mirror')]:
				mirrorpdb_fn = options['RUN']['pdbmirror', 'pdb_local_mirror'] + '/' + pdbi.lower()[
																					   1:3] + '/pdb' + pdbi.lower() + '.ent.gz'
				if os.path.exists(mirrorpdb_fn):
					shutil.copyfile(mirrorpdb_fn, pdb_filename + '.gz')
				os.system('gunzip ' + pdb_filename + '.gz')
			elif options['RUN'][('repo', 'from_repo')]:
				repofilename = options['RUN'][('repo', 'from_repo')] + locations['FSYS'][
					'PDBpdbs'] + pdbi.lower() + '.pdb'
				if os.path.exists(repofilename):
					shutil.copyfile(repofilename, pdb_filename)
			else:
				download_result = iterate_download(
					"{locations['EXTERNAL']['pdbdb_mmcif']}/{0}-assembly1.cif.gz".format(pdbi.lower()),
					#'https://files.rcsb.org/pub/pdb/data/structures/all/mmCIF/{0}.cif.gz'.format(pdbi.lower()),
					mmcif_filename, check="bio.mmcif", gzipped=True)
				download_result_2 = iterate_download(
					"{locations['EXTERNAL']['pdbdb_pdb']}/{0}.pdb1.gz".format(pdbi.lower()),
					#'https://files.rcsb.org/pub/pdb/data/structures/all/pdb/pdb{0}.ent.gz'.format(pdbi.lower()),
					pdb_filename, check="bio.pdb", gzipped=True)

			# Check if file is empty or corrupted
			pdb_ok, mmcif_ok = True, True
			parser = MMCIFParser()  # TODO check functional
			if not os.path.exists(mmcif_filename):
				mmcif_ok = False
				local_stats['no_mmcif'].append(pdbi.lower())
			elif os.path.getsize(mmcif_filename) == 0:
				mmcif_ok = False
				local_stats['empty_mmcif'].append(pdbi.lower())
			else:
				try:
					structure = parser.get_structure("try", mmcif_filename)
				except:
					mmcif_ok = False
					local_stats['corrupted_mmcif'].append(pdbi.lower())

			parser = PDBParser()  # TODO - check functional
			if not os.path.exists(pdb_filename):
				pdb_ok = False
				local_stats['no_pdb'].append(pdbi.lower())
			elif os.path.getsize(pdb_filename) == 0:
				pdb_ok = False
				local_stats['empty_pdb'].append(pdbi.lower())
			else:
				try:
					structure = parser.get_structure("running pdbparser on:", pdb_filename)
				except:
					pdb_ok = False
					local_stats['corrupted_pdb'].append(pdbi.lower())

			# Parse header
			if mmcif_ok:
				temparch[pdbi.lower()]['coordinate_file'] = mmcif_filename
				header, hreport = parse_mmCIF_header_wrap(mmcif_filename)
				temparch[pdbi.lower()]['header'] = header
				headers.append(temparch[pdbi.lower()]['header'])  # USED NOW! For deciding alpha-beta-other class
			elif pdb_ok and False:  # SKIP THIS OPTION FOR NOW - PARSER IS OBSOLETE
				temparch[pdbi.lower()]['coordinate_file'] = pdb_filename
				temparch[pdbi.lower()]['header'], hreport = parse_PDB_header(pdb_filename)  # not used so far
			else:
				local_stats['no_str'].append(pdbi.lower())

	mmCIF_header_stats(locations['FSYSPATH']['PDBstats'], headers)

	# Fill FROM_PDB with info collected in temparch and with those in the PDB header
	print("scan_PDB", "Fill data structure", "time", time.ctime())
	for pdbi in pdb_list:
		if options['EXECUTION'][('debug', 'debug')]:
			print(pdbi)
		if pdbi in local_stats['no_record']:
			str_data[pdbi]['PASSPORT'].append(passport_entry(this_name + '_no_graphql', pdbi,
															 "The GraphQL entry for {0} was not found. The PDB entry will be disregarded".format(
																 pdbi)))
			str_data[pdbi]['status'].append('pdb_eliminated')
			continue
		if pdbi in local_stats['no_str']:
			str_data[pdbi]['PASSPORT'].append(passport_entry(this_name + '_no_coordfile', pdbi,
															 "The PDB coordinate file for {0} was not found. The PDB entry will be disregarded".format(
																 pdbi)))
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
			wnores_msg = (
			'WARNING', this_name, f"PDB json entry {pdbi} does not have the 'struct' key. Name will be set to None")
			print_log(wnores_msg)
			from_pdb['name'] = None
			local_stats['no_title'].append(pdbi)

		if temparch[pdbi]['exptl'] is not None:
			from_pdb['experimental_method'] = ",".join(x['method'] for x in temparch[pdbi]['exptl'])
			wnores_msg = ('WARNING', this_name,
						  f"PDB json entry {pdbi} does not have the 'exptl' key. Experimental method will be set to None")
			print_log(wnores_msg)
		else:
			from_pdb['experimental_method'] = None
			local_stats['no_expmethod'].append(pdbi)

		if temparch[pdbi]['refine']:  # It is a list and can be empty
			# If no 'refine', there are 3 keywords missing, not just 'PDB_resolution'
			from_pdb['resolution'] = temparch[pdbi]['refine'][0]['ls_d_res_high']
		else:
			if temparch[pdbi]['header']['resolution']:
				from_pdb['resolution'] = temparch[pdbi]['header']['resolution']
			else:
				local_stats['no_resolution'].append(pdbi)
			wnores_msg = ('WARNING', this_name,
						  f"PDB json entry {pdbi} does not have the 'refine' key. Resolution will be set to None")
			print_log(wnores_msg)
			from_pdb['resolution'] = None
		if temparch[pdbi]['pdbx_database_status'] is not None:
			print(temparch[pdbi]['pdbx_database_status'])
			from_pdb['deposition_date'] = temparch[pdbi]['pdbx_database_status']['recvd_initial_deposition_date']
		else:
			local_stats['no_depo_date'].append(pdbi)
			wnores_msg = ('WARNING', this_name,
						  f"PDB json entry {pdbi} does not have the 'database_PDB_rev' key. Deposition date will be set to None")
			print_log(wnores_msg)
			from_pdb['deposition_date'] = None

		# TODO: database references: from where does the sequence come from?
		# Polymers
		new2oldch = {}
		polychains = set()
		alreadythere = set()
		# print(pdbi)
		if temparch[pdbi]['polymer_entities'] is not None:
			for ip, p in enumerate(temparch[pdbi]['polymer_entities']):
				chains = temparch[pdbi]['polymer_entities'][ip]['entity_poly']['pdbx_strand_id']
				for inst in p['polymer_entity_instances']:
					if inst['rcsb_polymer_entity_instance_container_identifiers'][
						'auth_asym_id'] in polychains:  # HETATMs of the same chain are now considered a separate chain, we don't want to have them
						continue
					if inst['rcsb_polymer_entity_instance_container_identifiers']['auth_asym_id'] in chains:
						new2oldch[inst['rcsb_polymer_entity_instance_container_identifiers']['asym_id']] = \
						inst['rcsb_polymer_entity_instance_container_identifiers']['auth_asym_id']
						polychains.add(inst['rcsb_polymer_entity_instance_container_identifiers']['auth_asym_id'])
					else:
						exit(1)
				for ch in sorted(list(polychains)):
					poly = FixedDict(from_pdb['polymers'].get_fdict_template())
					poly['chain'] = ch
					poly['sequence'] = temparch[pdbi]['polymer_entities'][ip]['entity_poly'][
						'pdbx_seq_one_letter_code_can']
					if type(temparch[pdbi]['polymer_entities'][ip]['rcsb_polymer_entity_container_identifiers'][
								'uniprot_ids']) == list:
						poly['uniprot_acc'] = ",".join(
							temparch[pdbi]['polymer_entities'][ip]['rcsb_polymer_entity_container_identifiers'][
								'uniprot_ids'])
					if type(temparch[pdbi]['polymer_entities'][ip]['rcsb_entity_source_organism']) == list:
						poly['taxonomy'] = ",".join(
							[x['ncbi_scientific_name'] if x['ncbi_scientific_name'] else "" for x in
							 temparch[pdbi]['polymer_entities'][ip]['rcsb_entity_source_organism']])
					from_pdb['polymers'].append(poly)
			# print("PDB CHAIN NAMES : AUTHOR CHAIN NAMES", pdbi, new2oldch)
		else:
			local_stats['no_entities'].append(pdbi)
			wnores_msg = ('WARNING', this_name,
						  f"PDB json entry {pdbi} does not have the 'polymer_entities' key. Polymer section will be empty")
			print_log(wnores_msg)

		# Biological assembly
		#        print(temparch[pdbi])
		#        print(temparch[pdbi]['assemblies'])
		if temparch[pdbi]['assemblies'] is not None:
			# Assemblies
			# NB: an assembly is a way to form a complex. Different assemblies are found when different hypotheses on the
			#  biological assembly have been formulated. Usually, we consider the first assembly (but we record all)
			oper_tr = FixedDict(from_pdb['translrots'].get_fdict_template())
			for ib, b in enumerate(temparch[pdbi]['assemblies']):
				ba = FixedDict(from_pdb['biological_assemblies'].get_fdict_template())
				#                print(temparch[pdbi]['assemblies'][ib])
				#                print(temparch[pdbi]['assemblies'][ib]['pdbx_struct_assembly'])
				#                print(temparch[pdbi]['assemblies'][ib]['pdbx_struct_assembly']['method_details'])
				#                print(type(ba['method']), type(temparch[pdbi]['assemblies'][ib]['pdbx_struct_assembly']['method_details']))
				ba['method'] = temparch[pdbi]['assemblies'][ib]['pdbx_struct_assembly']['method_details']
				ba['description'] = temparch[pdbi]['assemblies'][ib]['pdbx_struct_assembly']['details']
				# Transformations
				# NB: a transformation is a SET of matrices that applies on the SAME SET of chains. An assembly can have more than one transformation
				# (such as 3dh4), because different matrices can apply to different chains and still contribute giving the same complex
				for itr in range(len(temparch[pdbi]['assemblies'][ib]['pdbx_struct_assembly_gen'])):
					tr = FixedDict(ba['transformations'].get_fdict_template())
					tr['chains'] = [new2oldch[x] for x in
									temparch[pdbi]['assemblies'][ib]['pdbx_struct_assembly_gen'][itr]['asym_id_list'] if
									x in new2oldch]  # List of chains(new name) for each assembly matrix
					oper_list, oper_keys = read_transf_matrices_json(
						temparch[pdbi]['assemblies'][ib]['pdbx_struct_oper_list'])
					opexps = read_oper_expression(
						temparch[pdbi]['assemblies'][ib]['pdbx_struct_assembly_gen'][itr]['oper_expression'], oper_keys)
					# print("READ_OPER_EXPRESSION", pdbi, ib, itr, temparch[pdbi]['assemblies'][ib]['pdbx_struct_assembly_gen'][itr]['oper_expression'], oper_keys, opexps)
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
			op['matrix'] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
			op['vector'] = [0, 0, 0]
			from_pdb['translrots'].append(op)
			from_pdb['ktranslrots'].append("1")

		# TO DO: compare with biological assembly in PDB file and see if there are discrepancies between the two records
		# WARNING: that record is not correctly parsed! Make sure it is!

		str_data[pdbi]['FROM_PDB'] = from_pdb

	for x in sorted(local_stats):
		print("STATS", "scan_PDB", x, len(local_stats[x]), local_stats[x])
	#pickle.dump(str_data, open(str_data_pkl, 'wb'))
	print("STATS", "scan_OPM", "Finished", "time", time.ctime(), "\n")
	return str_data


def PDBTM_download(locations, pdbi):

	pdb_url = locations['EXTERNAL']['pdbtm_root'] + {1} + ".trpdb".format(pdbi[1:3], pdbi)
	#pdb_url = 'http://pdbtm.unitmp.org/api/v1/entry/{1}.trpdb'.format(pdbi[1:3], pdbi)
	pdb_filename = locations['FSYSPATH']['PDBTMpdbs'] + pdbi + '_PDBTM.pdb'
	pdb_downloaded = iterate_download(pdb_url, pdb_filename, gzipped=False)

	xml_url = locations['EXTERNAL']['pdbtm_root'] + {1} + ".xml".format(pdbi[1:3], pdbi)
	#xml_url = 'http://pdbtm.unitmp.org/api/v1/entry/{1}.xml'.format(pdbi[1:3], pdbi)
	xml_filename = locations['FSYSPATH']['PDBTMxmls'] + pdbi + '_PDBTM.xml'
	xml_downloaded = iterate_download(xml_url, xml_filename, gzipped=False)

	return pdb_downloaded and xml_downloaded


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


def pdb_query(path="pdb_query.json"):
	"""Returns the json query for the PDB API (2021-)
	"""
	with open(path) as pf:
		return pf.read()


def test_boh():
    options = Options()
    options.headless = True
#    service_log_path = "{0}/chromedriver.log".format(os.getcwd())
#    service_args = ['--verbose']
#    driver = webdriver.Chrome(options=options, service_log_path=service_log_path, service_args=service_args)
    driver = webdriver.Firefox(options=options)

    driver.get("https://opm.phar.umich.edu/")

    print("DONE")
    exit(1)


def test_for_redundancy():
    from encompass.pipeline.initialize_repository import initialize_repository
    options, locations = initialize_repository()
    ids, secd = mainlist_from_OPM(options, locations)
    print(ids, secd)
    str_data = scan_OPM(options, locations, ids)


def main():
	return


# ensures that main is only executed when the script is run directly, not when imported as a module
if __name__ == "__main__":
	main()
