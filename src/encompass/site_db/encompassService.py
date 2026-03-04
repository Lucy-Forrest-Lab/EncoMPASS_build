import pdb
from turtle import title
from encompass.site_db.database.database import SessionLocal, Base, engine
from encomapss.sited_db.dao.dao import ( save_pdbwhole, get_wholepdb, get_all_wholepdbs, update_wholepdb, delete_pdb, get_wholepdbByPDBCode, save_pdbChain)
from encomapss.site_db.models.models import *


def create_pdb():
    
    db = SessionLocal()

    print("hi")
    pdb_obj = PdbStructureWhole(title="Sample PDB", chains=2, pdb_code="1a0s")

    print(pdb_obj)
    child1 = PdbChainInformation(class_id="1a0s_P", ps_class="beta", ps_member="1a0s_P", tm_domains=4)
    child2 = PdbChainInformation(class_id="1a0s_Q", ps_class="beta", ps_member="1a0s_P", tm_domains=4)

    pdb_obj.pdb_chain_information.append(child1)
    pdb_obj.pdb_chain_information.append(child2)

    ceSymmClassification = PdbCeSymmClassification(algo_version = "", algorithm = "", symm_angle = "", symm_coverage = 0.00, symm_image_file = "", symm_levels = 0, symm_num_repeats = 0, symm_order = "", symm_pml_file = "", symm_repeats = "", symm_rmsd = 0.00, symm_tm_score = 0.00, symm_translation = "", symm_repeat_length = 0, json_file_path_3dmol = "", symm_aligned_residues = 0, symm_seed = 0, symm_unrefined_rmsd = 0.00, symm_unrefined_tm_score = 0.00)


    symmAlignedResidueList = PcscSymmAlignedResidueList(ccsc_symm_aligned_residue= 0.00)

    symmAngleList = PcscSymmAngleList(symm_angle=0.00)

    symmSeedList =PcscSymmSeedList(pcsc_symm_seed= 0)

    symmTranslationList = PcscSymmTranslationList(symm_translation = 0.0)

    symmUnrefinedRmsdList = PcscSymmUnrefinedRmsdList(pcsc_symm_unrefined_rmsd=0.00)

    symmUnrefinedTmScoreList = PcscSymmUnrefinedTmScoreList(pcsc_symm_unrefined_tm_score=0.00)

    symmAlignments = PdbSymmAlignments(sa_name = "", sa_seq="", sa_seq2="")

    symmalignDownloads = PdbSymmalignDownloads(file_path="")

    ceSymmClassification.pcsc_symm_aligned_residue_list.append(symmAlignedResidueList)
    ceSymmClassification.pcsc_symm_angle_list.append(symmAngleList)
    ceSymmClassification.pcsc_symm_seed_list.append(symmSeedList)
    ceSymmClassification.pcsc_symm_translation_list.append(symmTranslationList)
    ceSymmClassification.pcsc_symm_unrefined_rmsd_list.append(symmUnrefinedRmsdList)
    ceSymmClassification.pcsc_symm_unrefined_tm_score_list.append(symmUnrefinedTmScoreList)
    ceSymmClassification.pdb_symm_alignments.append(symmAlignments)
    ceSymmClassification.pdb_symmalign_downloads.append(symmalignDownloads)

    pdb_obj.pdb_ce_symm_classification.append(ceSymmClassification)

    ceSymdClassification = PdbCeSymdClassification( algo_version = "", algorithm = "", symd_angle = 0.00, symd_coverage = 0.00, symd_image_file = "", symd_is = 0.00, symd_is_angle = 0.00, symd_is_translation = 0.00, symd_order = 0, symd_pml_file = "", symd_rmsd = 0.00, symd_tm_score = 0.00, symd_translation = 0.00, symd_ztm_score = 0.00, json_file_path_3dmol = "", symd_aligned_residues = 0)

    symdAlignments = PdbSymdAlignments(sa_name = "", sa_seq = "", sa_seq2 ="", sa_seq3="", sa_seq4="")

    symdalignDownloads = PdbSymdalignDownloads(file_path="")

    ceSymdClassification.pdb_symd_alignments.append(symdAlignments)
    ceSymdClassification.pdb_symdalign_downloads.append(symdalignDownloads)

    pdb_obj.pdb_ce_symd_classification.append(ceSymdClassification)

    ntmchaindetails = PdbNtmChainsDetails(chain_id="P", ntm_domains=18, seq_name="1aos_P", seq="SGFEFHGYARSGVIMNDSGASTKSGAYITPAGETGGAIGRLGNQADTYVEMNLEHKQTLDNGATTRFKVMVADGQTSYNDWTASTSDLNVRQAFVELGNLPTFAGPFKGSTLWAGKRFDRDNFDIHWIDSDVVFLAGTGGGIYDVKWNDGLRSNFSLYGRNFGDIDDSSNSVQNYILTMNHFAGPLQMMVSGLRAKDNDERKDSNGNLAKGDAANTGVHALLGLHNDSFYGLRDGSSKTALLYGHGLGAEVKGIGSDGALRPGADTWRIASYGTTPLSENWSVAPAMLAQRSKDRYADGDSYQWATFNLRLIQAINQNFALAYEGSYQYMDLKPEGYNDRQAVNGSFYKLTFAPTFKVGSIGDFFSRPEIRFYTSWMDWSKKLNNYASDDALGSDGFNSGGEWSFGVQMETWF")

    ntmChainDomains = PdbNtmChainDomains(tm_domain_id= 1,  tmd_range= "")

    ntmchaindetails.pdb_ntm_chain_domains.append(ntmChainDomains)

    pdb_obj.pdb_ntm_chains_details.append(ntmchaindetails)

    uniport=PdbUniprotAcList(uniprot_access_code="P22340")

    pdb_obj.pdb_uniprot_ac_list.append(uniport)

    symmetryAnalysis = PdbSymmetryAnalysis(psa_chain_order= "C;M;D;N;E;O;G;Q;J;T;", psa_order="C2", psa_angle= "180.0", psa_translation= 0.0, psa_rmsd=0.04, psa_tmscore= 0.93, psa_coverage= 0.93, psa_num_repeats= 2, psa_levels= 1, psa_rep_length= 860, psa_quaternary_internal= "Quaternary", psa_topology= "Parallel", psa_repeats= "(C_1-36,O_1-36)", psa_axis_angle=-1.45, psa_aligned_residues= 1720)
    
    saAlignments = PdbAnalysisAlignments(aa_name="", aa_seq="", aa_seq2="")

    saDownloads = PdbAnalysisDownloads(file_path="")

    saImages = PdbAnalysisImages(image_file_and_path="", pml_file_and_path="", json_file_path_3dmol="", super_pml_file_and_path="")

    saAlignedResiduesList = PsaAlignedResiduesList(psa_aligned_residue=0)
    saAngleList = PsaAngleList(psa_angle=0.00)

    saAxisAngleList = PsaAxisAngleList(psa_axis_angle=0.00)

    saRmsdList = PsaRmsdList(psa_rmsd=0.00)

    saTmscoreList = PsaTmscoreList(psa_tmscore=0.00)

    saTranslationList = PsaTranslationList(psa_translation=0.00)

    symmetryAnalysis.pdb_analysis_downloads.append(saDownloads)

    symmetryAnalysis.pdb_analysis_alignments.append(saAlignments)

    symmetryAnalysis.pdb_analysis_images.append(saImages)

    symmetryAnalysis.psa_aligned_residues_list.append(saAlignedResiduesList)

    symmetryAnalysis.psa_angle_list.append(saAngleList)

    symmetryAnalysis.psa_axis_angle_list.append(saAxisAngleList)

    symmetryAnalysis.psa_rmsd_list.append(saRmsdList)

    symmetryAnalysis.psa_tmscore_list.append(saTmscoreList)

    symmetryAnalysis.psa_translation_list.append(saTranslationList)


    pdb_obj.pdb_symmetry_analysis.append(symmetryAnalysis)
    # ----------------------------------------------
    # Pass PDB object to DAO function
    # ----------------------------------------------
    saved_obj = save_pdbwhole(db, pdb_obj)

    print("Whole PDB saved:", saved_obj.psw_id)
    #print("Children:", saved_obj.chainInformation)

    db.close()
    print("Conection closed")


def create_pdb_chain():
    
    db = SessionLocal()

    pdb_obj = get_wholepdbByPDBCode(db, "1a0s") 

    print(pdb_obj)

    chain_obj = ChainStructure(tm_domains=3, opm_url="", pdb_code = "", pdb_url="", seq_name="", sequence = "", struct_pdb_download = "", struct_gif_path = "", struct_id = "1a0s_P", title = "", chain_size = 3, json_file_path_3dmol = '', image_file_path = "", resolution = 0.0, psw_id = pdb_obj.psw_id, pdbtm_url = "", method = "", chain_uniprot_access_codes = "")

    chainStructClassification = ChainStructClassification(ps_class = "", ps_member = "", tm_domains = 3, seq_neighbors = 3, struct_neighbors = 2, total_neighbors = 4, seq_neighbors_file_path = "", struct_neighbors_file_path = "", total_neighbors_file_path = "")
    
    chain_obj.chain_struct_classification.append(chainStructClassification)

    addntlFiles = ChainStructAddntlFiles(file_type = "", file_path = "")

    chain_obj.chain_struct_addntl_files.append(addntlFiles)

    structImgdatDownloads = ChainStructImgdatDownloads(file_type = "", file_path = "", map_path = "")

    chain_obj.chain_struct_imgdat_downloads.append(structImgdatDownloads)

    chainUniprot = ChainUniprotAcList(uniprot_access_code = "")

    chain_obj.chain_uniprot_ac_list.append(chainUniprot)

    chainCeSymdClassification = ChainCeSymdClassification(algo_version = "", algorithm = "", symd_angle = 0.0, symd_coverage = 0.0, symd_image_file = "", symd_is = 0.0, symd_is_angle = 0.0, symd_is_translation = 0.0, symd_order = 0.0, symd_pml_file = "", symd_rmsd = 0.0, symd_tm_score = 0.0, symd_translation = 0.0, symd_ztm_score = 0.0, json_file_path_3dmol = "", symd_aligned_residues = 0)

    symdAlignments = ChainSymdAlignments( sa_name = "", sa_seq = "", sa_seq2 = "")

    symdAlignDownloads = ChainSymdalignDownloads(file_path = "")

    chainCeSymdClassification.chain_symd_alignments.append(symdAlignments)

    chainCeSymdClassification.chain_symdalign_downloads.append(symdAlignDownloads)

    chain_obj.chain_ce_symd_classification.append(chainCeSymdClassification)


    chainCeSymmClassification = ChainCeSymmClassification(algo_version = "", algorithm = "", symm_angle = "", symm_coverage = 0.0, symm_image_file = "", symm_levels = 0, symm_num_repeats = 0, symm_order = "", symm_pml_file = "", symm_repeats = "", symm_rmsd = 0.0, symm_tm_score = 0.0, symm_translation = "", symm_repeat_length = 0, json_file_path_3dmol = "", symm_unrefined_rmsd = 0.0, symm_seed = 0, symm_unrefined_tm_score = 0.0, symm_aligned_residues = 0)

    symmAngle = CcscSymmAngleList(symm_angle = 0.0)

    symmTranslationList = CcscSymmTranslationList(symm_translation = 0.0)
    symmAlignments = ChainSymmAlignments(sa_name = "", sa_seq ="", sa_seq2="")
    symmalignDownloads = ChainSymmalignDownloads(file_path = "")

    chainCeSymmClassification.ccsc_symm_angle_list.append(symmAngle)
    chainCeSymmClassification.ccsc_symm_translation_list.append(symmTranslationList)
    chainCeSymmClassification.chain_symm_alignments.append(symmAlignments)
    chainCeSymmClassification.chain_symmalign_downloads.append(symmalignDownloads)

    chain_obj.chain_ce_symm_classification.append(chainCeSymmClassification)

    chainTransfer = ChainTransfer(ct_template = "", ct_order = "", ct_angle = "", ct_translation = "", ct_rmsd = "", ct_tmscore = "", ct_coverage = "", ct_num_repeats = "", ct_rep_length = "", ct_topology = "", ct_axis_angle = "", ct_repeats = "", ct_levels = "", ct_aligned_residues = "", ct_info = "")

    transferDownloads = ChainTransferDownloads(file_path = "")

    transferImages = ChainTransferImages(image_file_and_path = "", pml_file_and_path = "", json_file_path_3dmol = "", super_pml_file_and_path = "")
    alignedResidues = CtAlignedResiduesList(ct_aligned_residues = 0.0)
    transferAngle = CtAngleList(ct_angle = 0.0)
    transferAxisAngle = CtAxisAngleList(ct_axis_angle = 0.0)
    transferCoverage = CtCoverageList(ct_coverage = 0.0)
    transferLevel = CtLevelsList(ct_levels = 0)
    transferNumRepeat = CtNumRepeatsList(ct_num_repeats = 0)
    transferRepLength = CtRepLengthList(ct_rep_length = 0)
    transferRmsd = CtRmsdList(ct_rmsd = 0.0)
    transferTmscore = CtTmscoreList(ct_tmscore = 0.0)
    transferTranslation = CtTranslationList(ct_translation = 0.0)

    chainTransfer.chain_transfer_downloads.append(transferDownloads)
    chainTransfer.chain_transfer_images.append(transferImages)
    chainTransfer.ct_aligned_residues_list.append(alignedResidues)
    chainTransfer.ct_angle_list.append(transferAngle)
    chainTransfer.ct_axis_angle_list.append(transferAxisAngle)
    chainTransfer.ct_coverage_list.append(transferCoverage)
    chainTransfer.ct_levels_list.append(transferLevel)
    chainTransfer.ct_num_repeats_list.append(transferNumRepeat)
    chainTransfer.ct_rep_length_list.append(transferRepLength)
    chainTransfer.ct_rmsd_list.append(transferRmsd)
    chainTransfer.ct_tmscore_list.append(transferTmscore)
    chainTransfer.ct_translation_list.append(transferTranslation)

    chain_obj.chain_transfer.append(chainTransfer)


    ## Code to be added for ChainCeSymmClassification,  ChainCeSymdClassification, ChainSymmetryAnalysis tables  data 

    saved_obj = save_pdbChain(db, chain_obj)

    print("Chain PDB Details saved:", saved_obj.ps_id)

    db.close()
    print("Conection closed")

