from typing import Optional
import datetime

from sqlalchemy import Column, Date, DateTime, Double, ForeignKeyConstraint, Integer, PrimaryKeyConstraint, Sequence, String, Table, Text, text
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship

class Base(DeclarativeBase):
    pass


class ChainStructure(Base):
    __tablename__ = 'chain_structure'
    __table_args__ = (
        PrimaryKeyConstraint('ps_id', name='chain_structure_pk'),
        {'schema': 'protdb3'}
    )

    ps_id: Mapped[int] = mapped_column(Integer, Sequence('psid_seq', schema='protdb3'), primary_key=True)
    tm_domains: Mapped[Optional[int]] = mapped_column(Integer)
    opm_url: Mapped[Optional[str]] = mapped_column(String(255))
    pdb_code: Mapped[Optional[str]] = mapped_column(String(255))
    pdb_url: Mapped[Optional[str]] = mapped_column(String(255))
    seq_name: Mapped[Optional[str]] = mapped_column(String(255))
    sequence: Mapped[Optional[str]] = mapped_column(String(4000))
    struct_pdb_download: Mapped[Optional[str]] = mapped_column(String(255))
    struct_gif_path: Mapped[Optional[str]] = mapped_column(String(255))
    struct_id: Mapped[Optional[str]] = mapped_column(String(255))
    title: Mapped[Optional[str]] = mapped_column(String(255))
    chain_size: Mapped[Optional[int]] = mapped_column(Integer)
    json_file_path_3dmol: Mapped[Optional[str]] = mapped_column(String(255))
    image_file_path: Mapped[Optional[str]] = mapped_column(String(255))
    resolution: Mapped[Optional[float]] = mapped_column(Double(53))
    psw_id: Mapped[Optional[int]] = mapped_column(Integer)
    pdbtm_url: Mapped[Optional[str]] = mapped_column(String(255))
    method: Mapped[Optional[str]] = mapped_column(String(55))
    chain_uniprot_access_codes: Mapped[Optional[str]] = mapped_column(String(255))

    chain_ce_symd_classification: Mapped[list['ChainCeSymdClassification']] = relationship('ChainCeSymdClassification', back_populates='ps')
    chain_ce_symm_classification: Mapped[list['ChainCeSymmClassification']] = relationship('ChainCeSymmClassification', back_populates='ps')
    chain_struct_addntl_files: Mapped[list['ChainStructAddntlFiles']] = relationship('ChainStructAddntlFiles', back_populates='ps')
    chain_struct_classification: Mapped[list['ChainStructClassification']] = relationship('ChainStructClassification', back_populates='ps')
    chain_struct_imgdat_downloads: Mapped[list['ChainStructImgdatDownloads']] = relationship('ChainStructImgdatDownloads', back_populates='ps')
    chain_symmetry_analysis: Mapped[list['ChainSymmetryAnalysis']] = relationship('ChainSymmetryAnalysis', back_populates='ps')
    chain_tm_domains: Mapped[list['ChainTmDomains']] = relationship('ChainTmDomains', back_populates='ps')
    chain_transfer: Mapped[list['ChainTransfer']] = relationship('ChainTransfer', back_populates='ps')
    chain_uniprot_ac_list: Mapped[list['ChainUniprotAcList']] = relationship('ChainUniprotAcList', back_populates='ps')


class DbVersionRec(Base):
    __tablename__ = 'db_version_rec'
    __table_args__ = (
        PrimaryKeyConstraint('dv_id', name='db_version_rec_pk'),
        {'schema': 'protdb3'}
    )

    dv_id: Mapped[int] = mapped_column(Integer, Sequence('dvid_seq', schema='protdb3'), primary_key=True)
    db_version: Mapped[Optional[str]] = mapped_column(String(255))
    version_dt: Mapped[Optional[datetime.date]] = mapped_column(Date)
    db_file_name: Mapped[Optional[str]] = mapped_column(String(255))
    json_file_name: Mapped[Optional[str]] = mapped_column(String(255))
    comments: Mapped[Optional[str]] = mapped_column(String(255))
    db_sql_file_name: Mapped[Optional[str]] = mapped_column(String(55))
    xml_data_file_name: Mapped[Optional[str]] = mapped_column(String(55))


t_employee = Table(
    'employee', Base.metadata,
    Column('emp_name', Text, nullable=False),
    Column('salary', Integer, nullable=False),
    schema='protdb3'
)


class PdbStructureWhole(Base):
    __tablename__ = 'pdb_structure_whole'
    __table_args__ = (
        PrimaryKeyConstraint('psw_id', name='pdb_structure_whole_pk'),
        {'schema': 'protdb3'}
    )

    psw_id: Mapped[int] = mapped_column(Integer, Sequence('pswid_seq', schema='protdb3'), primary_key=True)
    chains: Mapped[Optional[int]] = mapped_column(Integer)
    opm_url: Mapped[Optional[str]] = mapped_column(String(255))
    pdb_code: Mapped[Optional[str]] = mapped_column(String(255))
    pdb_url: Mapped[Optional[str]] = mapped_column(String(255))
    struct_pdb_download: Mapped[Optional[str]] = mapped_column(String(255))
    struct_gif_path: Mapped[Optional[str]] = mapped_column(String(255))
    struct_id: Mapped[Optional[str]] = mapped_column(String)
    title: Mapped[Optional[str]] = mapped_column(String(255))
    struct_size: Mapped[Optional[int]] = mapped_column(Integer)
    image_file_path: Mapped[Optional[str]] = mapped_column(String(255))
    json_file_path_3dmol: Mapped[Optional[str]] = mapped_column(String(255))
    resolution: Mapped[Optional[float]] = mapped_column(Double(53))
    pdbtm_url: Mapped[Optional[str]] = mapped_column(String(255))
    tm_chains: Mapped[Optional[int]] = mapped_column(Integer)
    method: Mapped[Optional[str]] = mapped_column(String(55))
    uniprot_access_codes: Mapped[Optional[str]] = mapped_column(String(600))
    psw_info: Mapped[Optional[str]] = mapped_column(String(1000))

    pdb_ce_symd_classification: Mapped[list['PdbCeSymdClassification']] = relationship('PdbCeSymdClassification', back_populates='psw')
    pdb_ce_symm_classification: Mapped[list['PdbCeSymmClassification']] = relationship('PdbCeSymmClassification', back_populates='psw')
    pdb_chain_information: Mapped[list['PdbChainInformation']] = relationship('PdbChainInformation', back_populates='psw')
    pdb_ntm_chains_details: Mapped[list['PdbNtmChainsDetails']] = relationship('PdbNtmChainsDetails', back_populates='psw')
    pdb_symmetry_analysis: Mapped[list['PdbSymmetryAnalysis']] = relationship('PdbSymmetryAnalysis', back_populates='psw')
    pdb_uniprot_ac_list: Mapped[list['PdbUniprotAcList']] = relationship('PdbUniprotAcList', back_populates='psw')


class TestPdbStructureWhole(Base):
    __tablename__ = 'test_pdb_structure_whole'
    __table_args__ = (
        PrimaryKeyConstraint('psw_id', name='test_pdb_structure_whole_pk'),
        {'schema': 'protdb3'}
    )

    psw_id: Mapped[int] = mapped_column(Integer, Sequence('pswid_seq', schema='protdb3'), primary_key=True)
    chains: Mapped[Optional[int]] = mapped_column(Integer)
    opm_url: Mapped[Optional[str]] = mapped_column(String(255))
    pdb_code: Mapped[Optional[str]] = mapped_column(String(255))
    pdb_url: Mapped[Optional[str]] = mapped_column(String(255))
    struct_pdb_download: Mapped[Optional[str]] = mapped_column(String(255))
    struct_gif_path: Mapped[Optional[str]] = mapped_column(String(255))
    struct_id: Mapped[Optional[str]] = mapped_column(String)
    title: Mapped[Optional[str]] = mapped_column(String(255))
    struct_size: Mapped[Optional[int]] = mapped_column(Integer)
    image_file_path: Mapped[Optional[str]] = mapped_column(String(255))
    json_file_path_3dmol: Mapped[Optional[str]] = mapped_column(String(255))
    resolution: Mapped[Optional[float]] = mapped_column(Double(53))
    pdbtm_url: Mapped[Optional[str]] = mapped_column(String(255))
    tm_chains: Mapped[Optional[int]] = mapped_column(Integer)
    method: Mapped[Optional[str]] = mapped_column(String(55))
    uniprot_access_codes: Mapped[Optional[str]] = mapped_column(String(600))
    psw_info: Mapped[Optional[str]] = mapped_column(String(1000))

    test_pdb_chain_information: Mapped[list['TestPdbChainInformation']] = relationship('TestPdbChainInformation', back_populates='psw')


t_user_feedback_rec = Table(
    'user_feedback_rec', Base.metadata,
    Column('ufr_id', Integer, Sequence('ufrid_seq', schema='protdb3'), nullable=False),
    Column('enc_section', String(85)),
    Column('name', String(85)),
    Column('email', String(85)),
    Column('organization', String(255)),
    Column('comments', String(4000)),
    Column('structure_identifier', String(55)),
    Column('pdb_code', String(55)),
    Column('dt_submitted', DateTime),
    schema='protdb3'
)


t_v_chain_structure_chain_info = Table(
    'v_chain_structure_chain_info', Base.metadata,
    Column('psw_id', Integer),
    Column('pdb_code', String(255)),
    Column('title', String(255)),
    Column('chain_size', Integer),
    Column('seq_name', String(255)),
    Column('sequence', String(4000)),
    Column('chain_uniprot_access_codes', String(255)),
    Column('tm_domains', Integer),
    Column('seq_neighbors', Integer),
    Column('struct_neighbors', Integer),
    Column('total_neighbors', Integer),
    Column('symmetry', String(255)),
    Column('symm_consensus', String(255)),
    schema='protdb3'
)


t_v_chain_structure_classification = Table(
    'v_chain_structure_classification', Base.metadata,
    Column('pdb_code', String(255)),
    Column('title', String(255)),
    Column('resolution', Double(53)),
    Column('t1_tm_domains', Integer),
    Column('chain_size', Integer),
    Column('seq_name', String(255)),
    Column('sequence', String(4000)),
    Column('ps_class', String(255)),
    Column('t2_tm_domains', Integer),
    Column('seq_neighbors', Integer),
    Column('struct_neighbors', Integer),
    Column('total_neighbors', Integer),
    schema='protdb3'
)


class ChainCeSymdClassification(Base):
    __tablename__ = 'chain_ce_symd_classification'
    __table_args__ = (
        ForeignKeyConstraint(['ps_id'], ['protdb3.chain_structure.ps_id'], name='chain_ce_symd_classification_fk'),
        PrimaryKeyConstraint('csd_id', name='chain_ce_symd_classification_pk'),
        {'schema': 'protdb3'}
    )

    csd_id: Mapped[int] = mapped_column(Integer, Sequence('ccsdid_seq', schema='protdb3'), primary_key=True)
    algo_version: Mapped[Optional[str]] = mapped_column(String(255))
    algorithm: Mapped[Optional[str]] = mapped_column(String(255))
    symd_angle: Mapped[Optional[float]] = mapped_column(Double(53))
    symd_coverage: Mapped[Optional[float]] = mapped_column(Double(53))
    symd_image_file: Mapped[Optional[str]] = mapped_column(String(1000))
    symd_is: Mapped[Optional[float]] = mapped_column(Double(53))
    symd_is_angle: Mapped[Optional[float]] = mapped_column(Double(53))
    symd_is_translation: Mapped[Optional[float]] = mapped_column(Double(53))
    symd_order: Mapped[Optional[float]] = mapped_column(Double(53))
    symd_pml_file: Mapped[Optional[str]] = mapped_column(String(1000))
    symd_rmsd: Mapped[Optional[float]] = mapped_column(Double(53))
    symd_tm_score: Mapped[Optional[float]] = mapped_column(Double(53))
    symd_translation: Mapped[Optional[float]] = mapped_column(Double(53))
    symd_ztm_score: Mapped[Optional[float]] = mapped_column(Double(53))
    ps_id: Mapped[Optional[int]] = mapped_column(Integer)
    json_file_path_3dmol: Mapped[Optional[str]] = mapped_column(String(1000))
    symd_aligned_residues: Mapped[Optional[int]] = mapped_column(Integer)

    ps: Mapped[Optional['ChainStructure']] = relationship('ChainStructure', back_populates='chain_ce_symd_classification')
    chain_symd_alignments: Mapped[list['ChainSymdAlignments']] = relationship('ChainSymdAlignments', back_populates='csd')
    chain_symdalign_downloads: Mapped[list['ChainSymdalignDownloads']] = relationship('ChainSymdalignDownloads', back_populates='csd')


class ChainCeSymmClassification(Base):
    __tablename__ = 'chain_ce_symm_classification'
    __table_args__ = (
        ForeignKeyConstraint(['ps_id'], ['protdb3.chain_structure.ps_id'], name='chain_ce_symm_classification_fk'),
        PrimaryKeyConstraint('csc_id', name='chain_ce_symm_classification_pk'),
        {'schema': 'protdb3'}
    )

    csc_id: Mapped[int] = mapped_column(Integer, Sequence('cscmid_seq', schema='protdb3'), primary_key=True)
    ps_id: Mapped[int] = mapped_column(Integer, nullable=False)
    algo_version: Mapped[Optional[str]] = mapped_column(String(255))
    algorithm: Mapped[Optional[str]] = mapped_column(String(255))
    symm_angle: Mapped[Optional[str]] = mapped_column(String(255))
    symm_coverage: Mapped[Optional[float]] = mapped_column(Double(53))
    symm_image_file: Mapped[Optional[str]] = mapped_column(String(1000))
    symm_levels: Mapped[Optional[int]] = mapped_column(Integer)
    symm_num_repeats: Mapped[Optional[int]] = mapped_column(Integer)
    symm_order: Mapped[Optional[str]] = mapped_column(String(255))
    symm_pml_file: Mapped[Optional[str]] = mapped_column(String(1000))
    symm_repeats: Mapped[Optional[str]] = mapped_column(String(1000))
    symm_rmsd: Mapped[Optional[float]] = mapped_column(Double(53))
    symm_tm_score: Mapped[Optional[float]] = mapped_column(Double(53))
    symm_translation: Mapped[Optional[str]] = mapped_column(String(255))
    symm_repeat_length: Mapped[Optional[int]] = mapped_column(Integer)
    json_file_path_3dmol: Mapped[Optional[str]] = mapped_column(String(1000))
    symm_unrefined_rmsd: Mapped[Optional[float]] = mapped_column(Double(53))
    symm_seed: Mapped[Optional[int]] = mapped_column(Integer)
    symm_unrefined_tm_score: Mapped[Optional[float]] = mapped_column(Double(53))
    symm_aligned_residues: Mapped[Optional[int]] = mapped_column(Integer)

    ps: Mapped['ChainStructure'] = relationship('ChainStructure', back_populates='chain_ce_symm_classification')
    ccsc_symm_aligned_residue_list: Mapped[list['CcscSymmAlignedResidueList']] = relationship('CcscSymmAlignedResidueList', back_populates='csc')
    ccsc_symm_angle_list: Mapped[list['CcscSymmAngleList']] = relationship('CcscSymmAngleList', back_populates='csc')
    ccsc_symm_seed_list: Mapped[list['CcscSymmSeedList']] = relationship('CcscSymmSeedList', back_populates='csc')
    ccsc_symm_translation_list: Mapped[list['CcscSymmTranslationList']] = relationship('CcscSymmTranslationList', back_populates='csc')
    ccsc_symm_unrefined_rmsd_list: Mapped[list['CcscSymmUnrefinedRmsdList']] = relationship('CcscSymmUnrefinedRmsdList', back_populates='csc')
    ccsc_symm_unrefined_tm_score_list: Mapped[list['CcscSymmUnrefinedTmScoreList']] = relationship('CcscSymmUnrefinedTmScoreList', back_populates='csc')
    chain_symm_alignments: Mapped[list['ChainSymmAlignments']] = relationship('ChainSymmAlignments', back_populates='csc')
    chain_symmalign_downloads: Mapped[list['ChainSymmalignDownloads']] = relationship('ChainSymmalignDownloads', back_populates='csc')


class ChainStructAddntlFiles(Base):
    __tablename__ = 'chain_struct_addntl_files'
    __table_args__ = (
        ForeignKeyConstraint(['ps_id'], ['protdb3.chain_structure.ps_id'], name='chain_struct_addntl_files_fk'),
        PrimaryKeyConstraint('csaf_id', name='chain_struct_addntl_files_pk'),
        {'schema': 'protdb3'}
    )

    csaf_id: Mapped[int] = mapped_column(Integer, Sequence('csafid_seq', schema='protdb3'), primary_key=True)
    ps_id: Mapped[Optional[int]] = mapped_column(Integer)
    file_type: Mapped[Optional[str]] = mapped_column(String(65))
    file_path: Mapped[Optional[str]] = mapped_column(String(125))

    ps: Mapped[Optional['ChainStructure']] = relationship('ChainStructure', back_populates='chain_struct_addntl_files')


class ChainStructClassification(Base):
    __tablename__ = 'chain_struct_classification'
    __table_args__ = (
        ForeignKeyConstraint(['ps_id'], ['protdb3.chain_structure.ps_id'], name='chain_struct_classification_fk'),
        PrimaryKeyConstraint('sc_id', name='chain_struct_classification_pk'),
        {'schema': 'protdb3'}
    )

    sc_id: Mapped[int] = mapped_column(Integer, Sequence('scid_seq', schema='protdb3'), primary_key=True)
    ps_class: Mapped[Optional[str]] = mapped_column(String(255))
    ps_member: Mapped[Optional[str]] = mapped_column(String(255))
    tm_domains: Mapped[Optional[int]] = mapped_column(Integer)
    ps_id: Mapped[Optional[int]] = mapped_column(Integer)
    seq_neighbors: Mapped[Optional[int]] = mapped_column(Integer)
    struct_neighbors: Mapped[Optional[int]] = mapped_column(Integer)
    total_neighbors: Mapped[Optional[int]] = mapped_column(Integer)
    seq_neighbors_file_path: Mapped[Optional[str]] = mapped_column(String(100))
    struct_neighbors_file_path: Mapped[Optional[str]] = mapped_column(String(100))
    total_neighbors_file_path: Mapped[Optional[str]] = mapped_column(String(100))

    ps: Mapped[Optional['ChainStructure']] = relationship('ChainStructure', back_populates='chain_struct_classification')


class ChainStructImgdatDownloads(Base):
    __tablename__ = 'chain_struct_imgdat_downloads'
    __table_args__ = (
        ForeignKeyConstraint(['ps_id'], ['protdb3.chain_structure.ps_id'], name='chain_struct_imgdat_downloads_fk'),
        PrimaryKeyConstraint('cs_id', name='chain_struct_imgdat_downloads_pk'),
        {'schema': 'protdb3'}
    )

    cs_id: Mapped[int] = mapped_column(Integer, Sequence('csid_seq', schema='protdb3'), primary_key=True)
    ps_id: Mapped[Optional[int]] = mapped_column(Integer)
    file_type: Mapped[Optional[str]] = mapped_column(String(255))
    file_path: Mapped[Optional[str]] = mapped_column(String(255))
    map_path: Mapped[Optional[str]] = mapped_column(String(255))

    ps: Mapped[Optional['ChainStructure']] = relationship('ChainStructure', back_populates='chain_struct_imgdat_downloads')


class ChainSymmetryAnalysis(Base):
    __tablename__ = 'chain_symmetry_analysis'
    __table_args__ = (
        ForeignKeyConstraint(['ps_id'], ['protdb3.chain_structure.ps_id'], name='chain_symmetry_analysis_fk'),
        PrimaryKeyConstraint('csa_id', name='chain_symmetry_analysis_pk'),
        {'schema': 'protdb3'}
    )

    csa_id: Mapped[int] = mapped_column(Integer, Sequence('csaid_seq', schema='protdb3'), primary_key=True)
    ps_id: Mapped[Optional[int]] = mapped_column(Integer)
    csa_order: Mapped[Optional[str]] = mapped_column(String(55))
    csa_angle: Mapped[Optional[str]] = mapped_column(String(55))
    csa_translation: Mapped[Optional[str]] = mapped_column(String(55))
    csa_rmsd: Mapped[Optional[str]] = mapped_column(String(255))
    csa_tmscore: Mapped[Optional[str]] = mapped_column(String(255))
    csa_coverage: Mapped[Optional[str]] = mapped_column(String(255))
    csa_num_repeats: Mapped[Optional[str]] = mapped_column(String(255))
    csa_levels: Mapped[Optional[str]] = mapped_column(String(255))
    csa_rep_length: Mapped[Optional[str]] = mapped_column(String(255))
    csa_quaternary_internal: Mapped[Optional[str]] = mapped_column(String(55))
    csa_topology: Mapped[Optional[str]] = mapped_column(String(400))
    csa_repeats: Mapped[Optional[str]] = mapped_column(String(400))
    csa_info: Mapped[Optional[str]] = mapped_column(String(255))
    csa_axis_angle: Mapped[Optional[str]] = mapped_column(String(400))
    csa_aligned_residues: Mapped[Optional[str]] = mapped_column(String(400))
    csa_chain_order: Mapped[Optional[str]] = mapped_column(String(55))

    ps: Mapped[Optional['ChainStructure']] = relationship('ChainStructure', back_populates='chain_symmetry_analysis')
    chain_analysis_alignments: Mapped[list['ChainAnalysisAlignments']] = relationship('ChainAnalysisAlignments', back_populates='csa')
    chain_analysis_downloads: Mapped[list['ChainAnalysisDownloads']] = relationship('ChainAnalysisDownloads', back_populates='csa')
    chain_analysis_images: Mapped[list['ChainAnalysisImages']] = relationship('ChainAnalysisImages', back_populates='csa')
    csa_aligned_residues_list: Mapped[list['CsaAlignedResiduesList']] = relationship('CsaAlignedResiduesList', back_populates='csa')
    csa_angle_list: Mapped[list['CsaAngleList']] = relationship('CsaAngleList', back_populates='csa')
    csa_axis_angle_list: Mapped[list['CsaAxisAngleList']] = relationship('CsaAxisAngleList', back_populates='csa')
    csa_coverage_list: Mapped[list['CsaCoverageList']] = relationship('CsaCoverageList', back_populates='csa')
    csa_levels_list: Mapped[list['CsaLevelsList']] = relationship('CsaLevelsList', back_populates='csa')
    csa_num_repeats_list: Mapped[list['CsaNumRepeatsList']] = relationship('CsaNumRepeatsList', back_populates='csa')
    csa_rep_length_list: Mapped[list['CsaRepLengthList']] = relationship('CsaRepLengthList', back_populates='csa')
    csa_rmsd_list: Mapped[list['CsaRmsdList']] = relationship('CsaRmsdList', back_populates='csa')
    csa_tmscore_list: Mapped[list['CsaTmscoreList']] = relationship('CsaTmscoreList', back_populates='csa')
    csa_translation_list: Mapped[list['CsaTranslationList']] = relationship('CsaTranslationList', back_populates='csa')


class ChainTmDomains(Base):
    __tablename__ = 'chain_tm_domains'
    __table_args__ = (
        ForeignKeyConstraint(['ps_id'], ['protdb3.chain_structure.ps_id'], name='chain_tm_domains_fk'),
        PrimaryKeyConstraint('nd_id', name='chain_tm_domains_pk'),
        {'schema': 'protdb3'}
    )

    nd_id: Mapped[int] = mapped_column(Integer, Sequence('ndid_seq', schema='protdb3'), primary_key=True)
    tm_domain_id: Mapped[Optional[int]] = mapped_column(Integer)
    tmd_range: Mapped[Optional[str]] = mapped_column(String(255))
    ps_id: Mapped[Optional[int]] = mapped_column(Integer)

    ps: Mapped[Optional['ChainStructure']] = relationship('ChainStructure', back_populates='chain_tm_domains')


class ChainTransfer(Base):
    __tablename__ = 'chain_transfer'
    __table_args__ = (
        ForeignKeyConstraint(['ps_id'], ['protdb3.chain_structure.ps_id'], name='chain_transfer_fk'),
        PrimaryKeyConstraint('ct_id', name='chain_transfer_pk'),
        {'schema': 'protdb3'}
    )

    ct_id: Mapped[int] = mapped_column(Integer, Sequence('ctid_seq', schema='protdb3'), primary_key=True)
    ps_id: Mapped[Optional[int]] = mapped_column(Integer)
    ct_template: Mapped[Optional[str]] = mapped_column(String(55))
    ct_order: Mapped[Optional[str]] = mapped_column(String(55))
    ct_angle: Mapped[Optional[str]] = mapped_column(String(55))
    ct_translation: Mapped[Optional[str]] = mapped_column(String(55))
    ct_rmsd: Mapped[Optional[str]] = mapped_column(String(55))
    ct_tmscore: Mapped[Optional[str]] = mapped_column(String(55))
    ct_coverage: Mapped[Optional[str]] = mapped_column(String(55))
    ct_num_repeats: Mapped[Optional[str]] = mapped_column(String(55))
    ct_rep_length: Mapped[Optional[str]] = mapped_column(String(55))
    ct_topology: Mapped[Optional[str]] = mapped_column(String(55))
    ct_axis_angle: Mapped[Optional[str]] = mapped_column(String(55))
    ct_repeats: Mapped[Optional[str]] = mapped_column(String(255))
    ct_levels: Mapped[Optional[str]] = mapped_column(String(55))
    ct_aligned_residues: Mapped[Optional[str]] = mapped_column(String(255))
    ct_info: Mapped[Optional[str]] = mapped_column(String(255))

    ps: Mapped[Optional['ChainStructure']] = relationship('ChainStructure', back_populates='chain_transfer')
    chain_transfer_downloads: Mapped[list['ChainTransferDownloads']] = relationship('ChainTransferDownloads', back_populates='ct')
    chain_transfer_images: Mapped[list['ChainTransferImages']] = relationship('ChainTransferImages', back_populates='ct')
    ct_aligned_residues_list: Mapped[list['CtAlignedResiduesList']] = relationship('CtAlignedResiduesList', back_populates='ct')
    ct_angle_list: Mapped[list['CtAngleList']] = relationship('CtAngleList', back_populates='ct')
    ct_axis_angle_list: Mapped[list['CtAxisAngleList']] = relationship('CtAxisAngleList', back_populates='ct')
    ct_coverage_list: Mapped[list['CtCoverageList']] = relationship('CtCoverageList', back_populates='ct')
    ct_levels_list: Mapped[list['CtLevelsList']] = relationship('CtLevelsList', back_populates='ct')
    ct_num_repeats_list: Mapped[list['CtNumRepeatsList']] = relationship('CtNumRepeatsList', back_populates='ct')
    ct_rep_length_list: Mapped[list['CtRepLengthList']] = relationship('CtRepLengthList', back_populates='ct')
    ct_rmsd_list: Mapped[list['CtRmsdList']] = relationship('CtRmsdList', back_populates='ct')
    ct_tmscore_list: Mapped[list['CtTmscoreList']] = relationship('CtTmscoreList', back_populates='ct')
    ct_translation_list: Mapped[list['CtTranslationList']] = relationship('CtTranslationList', back_populates='ct')


class ChainUniprotAcList(Base):
    __tablename__ = 'chain_uniprot_ac_list'
    __table_args__ = (
        ForeignKeyConstraint(['ps_id'], ['protdb3.chain_structure.ps_id'], name='chain_uniprot_ac_list_fk'),
        PrimaryKeyConstraint('cual_id', name='chain_uniprot_ac_list_pk'),
        {'schema': 'protdb3'}
    )

    cual_id: Mapped[int] = mapped_column(Integer, Sequence('cualid_seq', schema='protdb3'), primary_key=True)
    ps_id: Mapped[int] = mapped_column(Integer, nullable=False)
    uniprot_access_code: Mapped[str] = mapped_column(String(20), nullable=False)

    ps: Mapped['ChainStructure'] = relationship('ChainStructure', back_populates='chain_uniprot_ac_list')


class PdbCeSymdClassification(Base):
    __tablename__ = 'pdb_ce_symd_classification'
    __table_args__ = (
        ForeignKeyConstraint(['psw_id'], ['protdb3.pdb_structure_whole.psw_id'], name='pdb_ce_symd_classification_fk'),
        PrimaryKeyConstraint('csd_id', name='pdb_ce_symd_classification_pk'),
        {'schema': 'protdb3'}
    )

    csd_id: Mapped[int] = mapped_column(Integer, Sequence('csdid_seq', schema='protdb3'), primary_key=True)
    algo_version: Mapped[Optional[str]] = mapped_column(String(255))
    algorithm: Mapped[Optional[str]] = mapped_column(String(255))
    symd_angle: Mapped[Optional[float]] = mapped_column(Double(53))
    symd_coverage: Mapped[Optional[float]] = mapped_column(Double(53))
    symd_image_file: Mapped[Optional[str]] = mapped_column(String(1800))
    symd_is: Mapped[Optional[float]] = mapped_column(Double(53))
    symd_is_angle: Mapped[Optional[float]] = mapped_column(Double(53))
    symd_is_translation: Mapped[Optional[float]] = mapped_column(Double(53))
    symd_order: Mapped[Optional[int]] = mapped_column(Integer)
    symd_pml_file: Mapped[Optional[str]] = mapped_column(String(1800))
    symd_rmsd: Mapped[Optional[float]] = mapped_column(Double(53))
    symd_tm_score: Mapped[Optional[float]] = mapped_column(Double(53))
    symd_translation: Mapped[Optional[float]] = mapped_column(Double(53))
    symd_ztm_score: Mapped[Optional[float]] = mapped_column(Double(53))
    psw_id: Mapped[Optional[int]] = mapped_column(Integer)
    json_file_path_3dmol: Mapped[Optional[str]] = mapped_column(String(1800))
    symd_aligned_residues: Mapped[Optional[int]] = mapped_column(Integer)

    psw: Mapped[Optional['PdbStructureWhole']] = relationship('PdbStructureWhole', back_populates='pdb_ce_symd_classification')
    pdb_symd_alignments: Mapped[list['PdbSymdAlignments']] = relationship('PdbSymdAlignments', back_populates='csd')
    pdb_symdalign_downloads: Mapped[list['PdbSymdalignDownloads']] = relationship('PdbSymdalignDownloads', back_populates='csd')


class PdbCeSymmClassification(Base):
    __tablename__ = 'pdb_ce_symm_classification'
    __table_args__ = (
        ForeignKeyConstraint(['psw_id'], ['protdb3.pdb_structure_whole.psw_id'], name='pdb_ce_symm_classification_fk'),
        PrimaryKeyConstraint('csc_id', name='pdb_ce_symm_classification_pk'),
        {'schema': 'protdb3'}
    )

    csc_id: Mapped[int] = mapped_column(Integer, Sequence('cscid_seq', schema='protdb3'), primary_key=True)
    algo_version: Mapped[Optional[str]] = mapped_column(String(255))
    algorithm: Mapped[Optional[str]] = mapped_column(String(255))
    symm_angle: Mapped[Optional[str]] = mapped_column(String(3000))
    symm_coverage: Mapped[Optional[float]] = mapped_column(Double(53))
    symm_image_file: Mapped[Optional[str]] = mapped_column(String(3000))
    symm_levels: Mapped[Optional[int]] = mapped_column(Integer)
    symm_num_repeats: Mapped[Optional[int]] = mapped_column(Integer)
    symm_order: Mapped[Optional[str]] = mapped_column(String(255))
    symm_pml_file: Mapped[Optional[str]] = mapped_column(String(3000))
    symm_repeats: Mapped[Optional[str]] = mapped_column(String(3000))
    symm_rmsd: Mapped[Optional[float]] = mapped_column(Double(53))
    symm_tm_score: Mapped[Optional[float]] = mapped_column(Double(53))
    symm_translation: Mapped[Optional[str]] = mapped_column(String(3000))
    psw_id: Mapped[Optional[int]] = mapped_column(Integer)
    symm_repeat_length: Mapped[Optional[int]] = mapped_column(Integer)
    json_file_path_3dmol: Mapped[Optional[str]] = mapped_column(String(3000))
    symm_aligned_residues: Mapped[Optional[int]] = mapped_column(Integer)
    symm_seed: Mapped[Optional[int]] = mapped_column(Integer)
    symm_unrefined_rmsd: Mapped[Optional[float]] = mapped_column(Double(53))
    symm_unrefined_tm_score: Mapped[Optional[float]] = mapped_column(Double(53))

    psw: Mapped[Optional['PdbStructureWhole']] = relationship('PdbStructureWhole', back_populates='pdb_ce_symm_classification')
    pcsc_symm_aligned_residue_list: Mapped[list['PcscSymmAlignedResidueList']] = relationship('PcscSymmAlignedResidueList', back_populates='csc')
    pcsc_symm_angle_list: Mapped[list['PcscSymmAngleList']] = relationship('PcscSymmAngleList', back_populates='csc')
    pcsc_symm_seed_list: Mapped[list['PcscSymmSeedList']] = relationship('PcscSymmSeedList', back_populates='csc')
    pcsc_symm_translation_list: Mapped[list['PcscSymmTranslationList']] = relationship('PcscSymmTranslationList', back_populates='csc')
    pcsc_symm_unrefined_rmsd_list: Mapped[list['PcscSymmUnrefinedRmsdList']] = relationship('PcscSymmUnrefinedRmsdList', back_populates='csc')
    pcsc_symm_unrefined_tm_score_list: Mapped[list['PcscSymmUnrefinedTmScoreList']] = relationship('PcscSymmUnrefinedTmScoreList', back_populates='csc')
    pdb_symm_alignments: Mapped[list['PdbSymmAlignments']] = relationship('PdbSymmAlignments', back_populates='csc')
    pdb_symmalign_downloads: Mapped[list['PdbSymmalignDownloads']] = relationship('PdbSymmalignDownloads', back_populates='csc')


class PdbChainInformation(Base):
    __tablename__ = 'pdb_chain_information'
    __table_args__ = (
        ForeignKeyConstraint(['psw_id'], ['protdb3.pdb_structure_whole.psw_id'], name='pdb_chain_information_fk'),
        PrimaryKeyConstraint('pcc_id', name='pdb_chain_information_pk'),
        {'schema': 'protdb3'}
    )

    pcc_id: Mapped[int] = mapped_column(Integer, Sequence('pccid_seq', schema='protdb3'), primary_key=True)
    class_id: Mapped[Optional[str]] = mapped_column(String(255))
    ps_class: Mapped[Optional[str]] = mapped_column(String(255))
    ps_member: Mapped[Optional[str]] = mapped_column(String(255))
    tm_domains: Mapped[Optional[int]] = mapped_column(Integer)
    psw_id: Mapped[Optional[int]] = mapped_column(Integer)
    seq_neighbors: Mapped[Optional[int]] = mapped_column(Integer)
    struct_neighbors: Mapped[Optional[int]] = mapped_column(Integer)
    total_neighbors: Mapped[Optional[int]] = mapped_column(Integer)
    symmetry: Mapped[Optional[str]] = mapped_column(String(255))
    symm_consensus: Mapped[Optional[str]] = mapped_column(String(255))
    seq_neighbors_file_path: Mapped[Optional[str]] = mapped_column(String(150))
    struct_neighbors_file_path: Mapped[Optional[str]] = mapped_column(String(150))
    total_neighbors_file_path: Mapped[Optional[str]] = mapped_column(String(150))

    psw: Mapped[Optional['PdbStructureWhole']] = relationship('PdbStructureWhole', back_populates='pdb_chain_information')


class PdbNtmChainsDetails(Base):
    __tablename__ = 'pdb_ntm_chains_details'
    __table_args__ = (
        ForeignKeyConstraint(['psw_id'], ['protdb3.pdb_structure_whole.psw_id'], name='pdb_ntm_chains_details_fk'),
        PrimaryKeyConstraint('ncd_id', name='pdb_ntm_chains_details_pk'),
        {'schema': 'protdb3'}
    )

    ncd_id: Mapped[int] = mapped_column(Integer, Sequence('ncdid_seq', schema='protdb3'), primary_key=True)
    chain_id: Mapped[Optional[str]] = mapped_column(String(255))
    ntm_domains: Mapped[Optional[int]] = mapped_column(Integer)
    seq: Mapped[Optional[str]] = mapped_column(String(4000))
    seq_name: Mapped[Optional[str]] = mapped_column(String(255))
    psw_id: Mapped[Optional[int]] = mapped_column(Integer)

    psw: Mapped[Optional['PdbStructureWhole']] = relationship('PdbStructureWhole', back_populates='pdb_ntm_chains_details')
    pdb_ntm_chain_domains: Mapped[list['PdbNtmChainDomains']] = relationship('PdbNtmChainDomains', back_populates='ncd')


class PdbSymmetryAnalysis(Base):
    __tablename__ = 'pdb_symmetry_analysis'
    __table_args__ = (
        ForeignKeyConstraint(['psw_id'], ['protdb3.pdb_structure_whole.psw_id'], name='pdb_symmetry_analysis_fk'),
        PrimaryKeyConstraint('psa_id', name='pdb_symmetry_analysis_pk'),
        {'schema': 'protdb3'}
    )

    psa_id: Mapped[int] = mapped_column(Integer, Sequence('psaid_seq', schema='protdb3'), primary_key=True)
    psw_id: Mapped[int] = mapped_column(Integer, nullable=False)
    psa_chain_order: Mapped[Optional[str]] = mapped_column(String(255))
    psa_order: Mapped[Optional[str]] = mapped_column(String(255))
    psa_angle: Mapped[Optional[str]] = mapped_column(String(255))
    psa_translation: Mapped[Optional[str]] = mapped_column(String(255))
    psa_rmsd: Mapped[Optional[float]] = mapped_column(Double(53))
    psa_tmscore: Mapped[Optional[float]] = mapped_column(Double(53))
    psa_coverage: Mapped[Optional[float]] = mapped_column(Double(53))
    psa_num_repeats: Mapped[Optional[int]] = mapped_column(Integer)
    psa_levels: Mapped[Optional[int]] = mapped_column(Integer)
    psa_rep_length: Mapped[Optional[int]] = mapped_column(Integer)
    psa_quaternary_internal: Mapped[Optional[str]] = mapped_column(String(255))
    psa_topology: Mapped[Optional[str]] = mapped_column(String(255))
    psa_repeats: Mapped[Optional[str]] = mapped_column(String(1000))
    psa_info: Mapped[Optional[str]] = mapped_column(String(1000))
    psa_axis_angle: Mapped[Optional[str]] = mapped_column(String(255))
    psa_aligned_residues: Mapped[Optional[int]] = mapped_column(Integer)

    psw: Mapped['PdbStructureWhole'] = relationship('PdbStructureWhole', back_populates='pdb_symmetry_analysis')
    pdb_analysis_alignments: Mapped[list['PdbAnalysisAlignments']] = relationship('PdbAnalysisAlignments', back_populates='psa')
    pdb_analysis_downloads: Mapped[list['PdbAnalysisDownloads']] = relationship('PdbAnalysisDownloads', back_populates='psa')
    pdb_analysis_images: Mapped[list['PdbAnalysisImages']] = relationship('PdbAnalysisImages', back_populates='psa')
    psa_aligned_residues_list: Mapped[list['PsaAlignedResiduesList']] = relationship('PsaAlignedResiduesList', back_populates='psa')
    psa_angle_list: Mapped[list['PsaAngleList']] = relationship('PsaAngleList', back_populates='psa')
    psa_axis_angle_list: Mapped[list['PsaAxisAngleList']] = relationship('PsaAxisAngleList', back_populates='psa')
    psa_rmsd_list: Mapped[list['PsaRmsdList']] = relationship('PsaRmsdList', back_populates='psa')
    psa_tmscore_list: Mapped[list['PsaTmscoreList']] = relationship('PsaTmscoreList', back_populates='psa')
    psa_translation_list: Mapped[list['PsaTranslationList']] = relationship('PsaTranslationList', back_populates='psa')


class PdbUniprotAcList(Base):
    __tablename__ = 'pdb_uniprot_ac_list'
    __table_args__ = (
        ForeignKeyConstraint(['psw_id'], ['protdb3.pdb_structure_whole.psw_id'], name='pdb_uniprot_ac_list_fk'),
        PrimaryKeyConstraint('pual_id', name='pdb_uniprot_ac_list_pk'),
        {'schema': 'protdb3'}
    )

    pual_id: Mapped[int] = mapped_column(Integer, Sequence('pualid_seq', schema='protdb3'), primary_key=True)
    psw_id: Mapped[int] = mapped_column(Integer, nullable=False)
    uniprot_access_code: Mapped[str] = mapped_column(String, nullable=False)

    psw: Mapped['PdbStructureWhole'] = relationship('PdbStructureWhole', back_populates='pdb_uniprot_ac_list')


class TestPdbChainInformation(Base):
    __tablename__ = 'test_pdb_chain_information'
    __table_args__ = (
        ForeignKeyConstraint(['psw_id'], ['protdb3.test_pdb_structure_whole.psw_id'], name='test_pdb_chain_information_fk'),
        PrimaryKeyConstraint('pcc_id', name='test_pdb_chain_information_pk'),
        {'schema': 'protdb3'}
    )

    pcc_id: Mapped[int] = mapped_column(Integer, Sequence('pccid_seq', schema='protdb3'), primary_key=True)
    class_id: Mapped[Optional[str]] = mapped_column(String(255))
    ps_class: Mapped[Optional[str]] = mapped_column(String(255))
    ps_member: Mapped[Optional[str]] = mapped_column(String(255))
    tm_domains: Mapped[Optional[int]] = mapped_column(Integer)
    psw_id: Mapped[Optional[int]] = mapped_column(Integer)
    seq_neighbors: Mapped[Optional[int]] = mapped_column(Integer)
    struct_neighbors: Mapped[Optional[int]] = mapped_column(Integer)
    total_neighbors: Mapped[Optional[int]] = mapped_column(Integer)
    symmetry: Mapped[Optional[str]] = mapped_column(String(255))
    symm_consensus: Mapped[Optional[str]] = mapped_column(String(255))
    seq_neighbors_file_path: Mapped[Optional[str]] = mapped_column(String(150))
    struct_neighbors_file_path: Mapped[Optional[str]] = mapped_column(String(150))
    total_neighbors_file_path: Mapped[Optional[str]] = mapped_column(String(150))

    psw: Mapped[Optional['TestPdbStructureWhole']] = relationship('TestPdbStructureWhole', back_populates='test_pdb_chain_information')


class CcscSymmAlignedResidueList(Base):
    __tablename__ = 'ccsc_symm_aligned_residue_list'
    __table_args__ = (
        ForeignKeyConstraint(['csc_id'], ['protdb3.chain_ce_symm_classification.csc_id'], name='ccsc_symm_aligned_residue_list_fk'),
        PrimaryKeyConstraint('csarl_id', name='ccsc_symm_aligned_residue_list_pk'),
        {'schema': 'protdb3'}
    )

    csarl_id: Mapped[int] = mapped_column(Integer, Sequence('csarlid_seq', schema='protdb3'), primary_key=True)
    csc_id: Mapped[Optional[int]] = mapped_column(Integer)
    ccsc_symm_aligned_residue: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    csc: Mapped[Optional['ChainCeSymmClassification']] = relationship('ChainCeSymmClassification', back_populates='ccsc_symm_aligned_residue_list')


class CcscSymmAngleList(Base):
    __tablename__ = 'ccsc_symm_angle_list'
    __table_args__ = (
        ForeignKeyConstraint(['csc_id'], ['protdb3.chain_ce_symm_classification.csc_id'], name='ccsc_symm_angle_list_fk'),
        PrimaryKeyConstraint('csal_id', name='ccsc_symm_angle_list_pk'),
        {'schema': 'protdb3'}
    )

    csal_id: Mapped[int] = mapped_column(Integer, Sequence('csalid1_seq', schema='protdb3'), primary_key=True)
    csc_id: Mapped[Optional[int]] = mapped_column(Integer)
    symm_angle: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    csc: Mapped[Optional['ChainCeSymmClassification']] = relationship('ChainCeSymmClassification', back_populates='ccsc_symm_angle_list')


class CcscSymmSeedList(Base):
    __tablename__ = 'ccsc_symm_seed_list'
    __table_args__ = (
        ForeignKeyConstraint(['csc_id'], ['protdb3.chain_ce_symm_classification.csc_id'], name='ccsc_symm_seed_list_fk'),
        PrimaryKeyConstraint('cssl_id', name='ccsc_symm_seed_list_pk'),
        {'schema': 'protdb3'}
    )

    cssl_id: Mapped[int] = mapped_column(Integer, Sequence('csslid_seq', schema='protdb3'), primary_key=True)
    csc_id: Mapped[Optional[int]] = mapped_column(Integer)
    ccsc_symm_seed: Mapped[Optional[int]] = mapped_column(Integer)

    csc: Mapped[Optional['ChainCeSymmClassification']] = relationship('ChainCeSymmClassification', back_populates='ccsc_symm_seed_list')


class CcscSymmTranslationList(Base):
    __tablename__ = 'ccsc_symm_translation_list'
    __table_args__ = (
        ForeignKeyConstraint(['csc_id'], ['protdb3.chain_ce_symm_classification.csc_id'], name='ccsc_symm_translation_list_fk'),
        PrimaryKeyConstraint('cstl_id', name='ccsc_symm_translation_list_pk'),
        {'schema': 'protdb3'}
    )

    cstl_id: Mapped[int] = mapped_column(Integer, Sequence('ccstlid1_seq', schema='protdb3'), primary_key=True)
    csc_id: Mapped[Optional[int]] = mapped_column(Integer)
    symm_translation: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    csc: Mapped[Optional['ChainCeSymmClassification']] = relationship('ChainCeSymmClassification', back_populates='ccsc_symm_translation_list')


class CcscSymmUnrefinedRmsdList(Base):
    __tablename__ = 'ccsc_symm_unrefined_rmsd_list'
    __table_args__ = (
        ForeignKeyConstraint(['csc_id'], ['protdb3.chain_ce_symm_classification.csc_id'], name='ccsc_symm_unrefined_rmsd_list_fk'),
        PrimaryKeyConstraint('csurl_id', name='ccsc_symm_unrefined_rmsd_list_pk'),
        {'schema': 'protdb3'}
    )

    csurl_id: Mapped[int] = mapped_column(Integer, Sequence('csurlid_seq', schema='protdb3'), primary_key=True)
    csc_id: Mapped[Optional[int]] = mapped_column(Integer)
    ccsc_symm_unrefined_rmsd: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    csc: Mapped[Optional['ChainCeSymmClassification']] = relationship('ChainCeSymmClassification', back_populates='ccsc_symm_unrefined_rmsd_list')


class CcscSymmUnrefinedTmScoreList(Base):
    __tablename__ = 'ccsc_symm_unrefined_tm_score_list'
    __table_args__ = (
        ForeignKeyConstraint(['csc_id'], ['protdb3.chain_ce_symm_classification.csc_id'], name='ccsc_symm_unrefined_tm_score_list_fk'),
        PrimaryKeyConstraint('csutsl_id', name='ccsc_symm_unrefined_tm_score_list_pk'),
        {'schema': 'protdb3'}
    )

    csutsl_id: Mapped[int] = mapped_column(Integer, Sequence('csutslid_seq', schema='protdb3'), primary_key=True)
    csc_id: Mapped[Optional[int]] = mapped_column(Integer)
    ccsc_symm_unrefined_tm_score: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    csc: Mapped[Optional['ChainCeSymmClassification']] = relationship('ChainCeSymmClassification', back_populates='ccsc_symm_unrefined_tm_score_list')


class ChainAnalysisAlignments(Base):
    __tablename__ = 'chain_analysis_alignments'
    __table_args__ = (
        ForeignKeyConstraint(['csa_id'], ['protdb3.chain_symmetry_analysis.csa_id'], name='chain_analysis_alignments_fk'),
        PrimaryKeyConstraint('caa_id', name='chain_analysis_alignments_pk'),
        {'schema': 'protdb3'}
    )

    caa_id: Mapped[int] = mapped_column(Integer, Sequence('caaid1_seq', schema='protdb3'), primary_key=True)
    aa_name: Mapped[Optional[str]] = mapped_column(String(1500))
    aa_seq: Mapped[Optional[str]] = mapped_column(String(4000))
    csa_id: Mapped[Optional[int]] = mapped_column(Integer)
    aa_seq2: Mapped[Optional[str]] = mapped_column(String(4000))

    csa: Mapped[Optional['ChainSymmetryAnalysis']] = relationship('ChainSymmetryAnalysis', back_populates='chain_analysis_alignments')


class ChainAnalysisDownloads(Base):
    __tablename__ = 'chain_analysis_downloads'
    __table_args__ = (
        ForeignKeyConstraint(['csa_id'], ['protdb3.chain_symmetry_analysis.csa_id'], name='chain_analysis_downloads_fk'),
        PrimaryKeyConstraint('cad_id', name='chain_analysis_downloads_pk'),
        {'schema': 'protdb3'}
    )

    cad_id: Mapped[int] = mapped_column(Integer, Sequence('cadid_seq', schema='protdb3'), primary_key=True)
    csa_id: Mapped[Optional[int]] = mapped_column(Integer)
    file_path: Mapped[Optional[str]] = mapped_column(String(85))

    csa: Mapped[Optional['ChainSymmetryAnalysis']] = relationship('ChainSymmetryAnalysis', back_populates='chain_analysis_downloads')


class ChainAnalysisImages(Base):
    __tablename__ = 'chain_analysis_images'
    __table_args__ = (
        ForeignKeyConstraint(['csa_id'], ['protdb3.chain_symmetry_analysis.csa_id'], name='chain_analysis_images_fk'),
        PrimaryKeyConstraint('cai_id', name='chain_analysis_images_pk'),
        {'schema': 'protdb3'}
    )

    cai_id: Mapped[int] = mapped_column(Integer, Sequence('caiid_seq', schema='protdb3'), primary_key=True)
    csa_id: Mapped[Optional[int]] = mapped_column(Integer)
    image_file_and_path: Mapped[Optional[str]] = mapped_column(String(1000))
    pml_file_and_path: Mapped[Optional[str]] = mapped_column(String(1000))
    json_file_path_3dmol: Mapped[Optional[str]] = mapped_column(String(1000))
    super_pml_file_and_path: Mapped[Optional[str]] = mapped_column(String(1000))

    csa: Mapped[Optional['ChainSymmetryAnalysis']] = relationship('ChainSymmetryAnalysis', back_populates='chain_analysis_images')


class ChainSymdAlignments(Base):
    __tablename__ = 'chain_symd_alignments'
    __table_args__ = (
        ForeignKeyConstraint(['csd_id'], ['protdb3.chain_ce_symd_classification.csd_id'], name='chain_symd_alignments_fk'),
        PrimaryKeyConstraint('sa_id', name='chain_symd_alignments_pk'),
        {'schema': 'protdb3'}
    )

    sa_id: Mapped[int] = mapped_column(Integer, Sequence('csyaid1_seq', schema='protdb3'), primary_key=True)
    sa_name: Mapped[Optional[str]] = mapped_column(String(255))
    sa_seq: Mapped[Optional[str]] = mapped_column(String(4000))
    csd_id: Mapped[Optional[int]] = mapped_column(Integer)
    sa_seq2: Mapped[Optional[str]] = mapped_column(String(4000))

    csd: Mapped[Optional['ChainCeSymdClassification']] = relationship('ChainCeSymdClassification', back_populates='chain_symd_alignments')


class ChainSymdalignDownloads(Base):
    __tablename__ = 'chain_symdalign_downloads'
    __table_args__ = (
        ForeignKeyConstraint(['csd_id'], ['protdb3.chain_ce_symd_classification.csd_id'], name='chain_symdalign_downloads_fk'),
        PrimaryKeyConstraint('csdd_id', name='chain_symdalign_downloads_pk'),
        {'schema': 'protdb3'}
    )

    csdd_id: Mapped[int] = mapped_column(Integer, Sequence('csddid_seq', schema='protdb3'), primary_key=True)
    csd_id: Mapped[Optional[int]] = mapped_column(Integer)
    file_path: Mapped[Optional[str]] = mapped_column(String(255))

    csd: Mapped[Optional['ChainCeSymdClassification']] = relationship('ChainCeSymdClassification', back_populates='chain_symdalign_downloads')


class ChainSymmAlignments(Base):
    __tablename__ = 'chain_symm_alignments'
    __table_args__ = (
        ForeignKeyConstraint(['csc_id'], ['protdb3.chain_ce_symm_classification.csc_id'], name='chain_symm_alignments_fk'),
        PrimaryKeyConstraint('sma_id', name='chain_symm_alignments_pk'),
        {'schema': 'protdb3'}
    )

    sma_id: Mapped[int] = mapped_column(Integer, Sequence('csmaid_seq', schema='protdb3'), primary_key=True)
    sa_name: Mapped[Optional[str]] = mapped_column(String(255))
    sa_seq: Mapped[Optional[str]] = mapped_column(String(4000))
    csc_id: Mapped[Optional[int]] = mapped_column(Integer)
    sa_seq2: Mapped[Optional[str]] = mapped_column(String(4000))

    csc: Mapped[Optional['ChainCeSymmClassification']] = relationship('ChainCeSymmClassification', back_populates='chain_symm_alignments')


class ChainSymmalignDownloads(Base):
    __tablename__ = 'chain_symmalign_downloads'
    __table_args__ = (
        ForeignKeyConstraint(['csc_id'], ['protdb3.chain_ce_symm_classification.csc_id'], name='chain_symmalign_downloads_fk'),
        PrimaryKeyConstraint('csad_id', name='chain_symmalign_downloads_pk'),
        {'schema': 'protdb3'}
    )

    csad_id: Mapped[int] = mapped_column(Integer, Sequence('csadid_seq', schema='protdb3'), primary_key=True)
    csc_id: Mapped[Optional[int]] = mapped_column(Integer)
    file_path: Mapped[Optional[str]] = mapped_column(String(255))

    csc: Mapped[Optional['ChainCeSymmClassification']] = relationship('ChainCeSymmClassification', back_populates='chain_symmalign_downloads')


class ChainTransferDownloads(Base):
    __tablename__ = 'chain_transfer_downloads'
    __table_args__ = (
        ForeignKeyConstraint(['ct_id'], ['protdb3.chain_transfer.ct_id'], name='chain_transfer_downloads_fk'),
        PrimaryKeyConstraint('ctd_id', name='chain_transfer_downloads_pk'),
        {'schema': 'protdb3'}
    )

    ctd_id: Mapped[int] = mapped_column(Integer, Sequence('ctdid_seq', schema='protdb3'), primary_key=True)
    ct_id: Mapped[Optional[int]] = mapped_column(Integer)
    file_path: Mapped[Optional[str]] = mapped_column(String(255))

    ct: Mapped[Optional['ChainTransfer']] = relationship('ChainTransfer', back_populates='chain_transfer_downloads')


class ChainTransferImages(Base):
    __tablename__ = 'chain_transfer_images'
    __table_args__ = (
        ForeignKeyConstraint(['ct_id'], ['protdb3.chain_transfer.ct_id'], name='chain_transfer_images_fk'),
        PrimaryKeyConstraint('ctm_id', name='chain_transfer_images_pk'),
        {'schema': 'protdb3'}
    )

    ctm_id: Mapped[int] = mapped_column(Integer, Sequence('ctmid_seq', schema='protdb3'), primary_key=True)
    ct_id: Mapped[Optional[int]] = mapped_column(Integer)
    image_file_and_path: Mapped[Optional[str]] = mapped_column(String(1000))
    pml_file_and_path: Mapped[Optional[str]] = mapped_column(String(1000))
    json_file_path_3dmol: Mapped[Optional[str]] = mapped_column(String(1000))
    super_pml_file_and_path: Mapped[Optional[str]] = mapped_column(String(1000))

    ct: Mapped[Optional['ChainTransfer']] = relationship('ChainTransfer', back_populates='chain_transfer_images')


class CsaAlignedResiduesList(Base):
    __tablename__ = 'csa_aligned_residues_list'
    __table_args__ = (
        ForeignKeyConstraint(['csa_id'], ['protdb3.chain_symmetry_analysis.csa_id'], name='csa_aligned_residues_list_fk'),
        PrimaryKeyConstraint('carl_id', name='csa_aligned_residues_list_pk'),
        {'schema': 'protdb3'}
    )

    carl_id: Mapped[int] = mapped_column(Integer, Sequence('carlid_seq', schema='protdb3'), primary_key=True)
    csa_id: Mapped[Optional[int]] = mapped_column(Integer)
    csa_aligned_residue: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    csa: Mapped[Optional['ChainSymmetryAnalysis']] = relationship('ChainSymmetryAnalysis', back_populates='csa_aligned_residues_list')


class CsaAngleList(Base):
    __tablename__ = 'csa_angle_list'
    __table_args__ = (
        ForeignKeyConstraint(['csa_id'], ['protdb3.chain_symmetry_analysis.csa_id'], name='csa_angle_list_fk'),
        PrimaryKeyConstraint('cal_id', name='csa_angle_list_pk'),
        {'schema': 'protdb3'}
    )

    cal_id: Mapped[int] = mapped_column(Integer, Sequence('calid_seq', schema='protdb3'), primary_key=True)
    csa_id: Mapped[Optional[int]] = mapped_column(Integer)
    csa_angle: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    csa: Mapped[Optional['ChainSymmetryAnalysis']] = relationship('ChainSymmetryAnalysis', back_populates='csa_angle_list')


class CsaAxisAngleList(Base):
    __tablename__ = 'csa_axis_angle_list'
    __table_args__ = (
        ForeignKeyConstraint(['csa_id'], ['protdb3.chain_symmetry_analysis.csa_id'], name='csa_axis_angle_list_fk'),
        PrimaryKeyConstraint('caal_id', name='csa_axis_angle_list_pk'),
        {'schema': 'protdb3'}
    )

    caal_id: Mapped[int] = mapped_column(Integer, Sequence('caalid_seq', schema='protdb3'), primary_key=True)
    csa_id: Mapped[Optional[int]] = mapped_column(Integer)
    csa_axis_angle: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    csa: Mapped[Optional['ChainSymmetryAnalysis']] = relationship('ChainSymmetryAnalysis', back_populates='csa_axis_angle_list')


class CsaCoverageList(Base):
    __tablename__ = 'csa_coverage_list'
    __table_args__ = (
        ForeignKeyConstraint(['csa_id'], ['protdb3.chain_symmetry_analysis.csa_id'], name='csa_coverage_list_fk'),
        PrimaryKeyConstraint('csacl_id', name='csa_coverage_list_pk'),
        {'schema': 'protdb3'}
    )

    csacl_id: Mapped[int] = mapped_column(Integer, Sequence('csaclid_seq', schema='protdb3'), primary_key=True)
    csa_id: Mapped[Optional[int]] = mapped_column(Integer)
    csa_coverage: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    csa: Mapped[Optional['ChainSymmetryAnalysis']] = relationship('ChainSymmetryAnalysis', back_populates='csa_coverage_list')


class CsaLevelsList(Base):
    __tablename__ = 'csa_levels_list'
    __table_args__ = (
        ForeignKeyConstraint(['csa_id'], ['protdb3.chain_symmetry_analysis.csa_id'], name='csa_levels_list_fk'),
        PrimaryKeyConstraint('csls_id', name='csa_levels_list_pk'),
        {'schema': 'protdb3'}
    )

    csls_id: Mapped[int] = mapped_column(Integer, Sequence('cslsid_seq', schema='protdb3'), primary_key=True)
    csa_id: Mapped[Optional[int]] = mapped_column(Integer)
    csa_levels: Mapped[Optional[int]] = mapped_column(Integer, server_default=text('0'))

    csa: Mapped[Optional['ChainSymmetryAnalysis']] = relationship('ChainSymmetryAnalysis', back_populates='csa_levels_list')


class CsaNumRepeatsList(Base):
    __tablename__ = 'csa_num_repeats_list'
    __table_args__ = (
        ForeignKeyConstraint(['csa_id'], ['protdb3.chain_symmetry_analysis.csa_id'], name='csa_num_repeats_list_fk'),
        PrimaryKeyConstraint('csnr_id', name='csa_num_repeats_list_pk'),
        {'schema': 'protdb3'}
    )

    csnr_id: Mapped[int] = mapped_column(Integer, Sequence('csnrid_seq', schema='protdb3'), primary_key=True)
    csa_id: Mapped[Optional[int]] = mapped_column(Integer)
    csa_num_repeats: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    csa: Mapped[Optional['ChainSymmetryAnalysis']] = relationship('ChainSymmetryAnalysis', back_populates='csa_num_repeats_list')


class CsaRepLengthList(Base):
    __tablename__ = 'csa_rep_length_list'
    __table_args__ = (
        ForeignKeyConstraint(['csa_id'], ['protdb3.chain_symmetry_analysis.csa_id'], name='csa_rep_length_list_fk'),
        PrimaryKeyConstraint('csrll_id', name='csa_rep_length_list_pk'),
        {'schema': 'protdb3'}
    )

    csrll_id: Mapped[int] = mapped_column(Integer, Sequence('csrllid_seq', schema='protdb3'), primary_key=True)
    csa_id: Mapped[Optional[int]] = mapped_column(Integer)
    csa_rep_length: Mapped[Optional[int]] = mapped_column(Integer, server_default=text('0'))

    csa: Mapped[Optional['ChainSymmetryAnalysis']] = relationship('ChainSymmetryAnalysis', back_populates='csa_rep_length_list')


class CsaRmsdList(Base):
    __tablename__ = 'csa_rmsd_list'
    __table_args__ = (
        ForeignKeyConstraint(['csa_id'], ['protdb3.chain_symmetry_analysis.csa_id'], name='csa_rmsd_list_fk'),
        PrimaryKeyConstraint('csrl_id', name='csa_rmsd_list_pk'),
        {'schema': 'protdb3'}
    )

    csrl_id: Mapped[int] = mapped_column(Integer, Sequence('csrlid_seq', schema='protdb3'), primary_key=True)
    csa_id: Mapped[Optional[int]] = mapped_column(Integer)
    csa_rmsd: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    csa: Mapped[Optional['ChainSymmetryAnalysis']] = relationship('ChainSymmetryAnalysis', back_populates='csa_rmsd_list')


class CsaTmscoreList(Base):
    __tablename__ = 'csa_tmscore_list'
    __table_args__ = (
        ForeignKeyConstraint(['csa_id'], ['protdb3.chain_symmetry_analysis.csa_id'], name='csa_tmscore_list_fk'),
        PrimaryKeyConstraint('cstl_id', name='csa_tmscore_list_pk'),
        {'schema': 'protdb3'}
    )

    cstl_id: Mapped[int] = mapped_column(Integer, Sequence('cstlid_seq', schema='protdb3'), primary_key=True)
    csa_id: Mapped[Optional[int]] = mapped_column(Integer)
    csa_tmscore: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    csa: Mapped[Optional['ChainSymmetryAnalysis']] = relationship('ChainSymmetryAnalysis', back_populates='csa_tmscore_list')


class CsaTranslationList(Base):
    __tablename__ = 'csa_translation_list'
    __table_args__ = (
        ForeignKeyConstraint(['csa_id'], ['protdb3.chain_symmetry_analysis.csa_id'], name='csa_translation_list_fk'),
        PrimaryKeyConstraint('ctl_id', name='csa_translation_list_pk'),
        {'schema': 'protdb3'}
    )

    ctl_id: Mapped[int] = mapped_column(Integer, Sequence('ctlid_seq', schema='protdb3'), primary_key=True)
    csa_id: Mapped[Optional[int]] = mapped_column(Integer)
    csa_translation: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    csa: Mapped[Optional['ChainSymmetryAnalysis']] = relationship('ChainSymmetryAnalysis', back_populates='csa_translation_list')


class CtAlignedResiduesList(Base):
    __tablename__ = 'ct_aligned_residues_list'
    __table_args__ = (
        ForeignKeyConstraint(['ct_id'], ['protdb3.chain_transfer.ct_id'], name='ct_aligned_residues_list_fk'),
        PrimaryKeyConstraint('car_id', name='ct_aligned_residues_list_pk'),
        {'schema': 'protdb3'}
    )

    car_id: Mapped[int] = mapped_column(Integer, Sequence('carid_seq', schema='protdb3'), primary_key=True)
    ct_id: Mapped[Optional[int]] = mapped_column(Integer)
    ct_aligned_residues: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    ct: Mapped[Optional['ChainTransfer']] = relationship('ChainTransfer', back_populates='ct_aligned_residues_list')


class CtAngleList(Base):
    __tablename__ = 'ct_angle_list'
    __table_args__ = (
        ForeignKeyConstraint(['ct_id'], ['protdb3.chain_transfer.ct_id'], name='ct_angle_list_fk'),
        PrimaryKeyConstraint('cta_id', name='ct_angle_list_pk'),
        {'schema': 'protdb3'}
    )

    cta_id: Mapped[int] = mapped_column(Integer, Sequence('ctaid_seq', schema='protdb3'), primary_key=True)
    ct_id: Mapped[Optional[int]] = mapped_column(Integer)
    ct_angle: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    ct: Mapped[Optional['ChainTransfer']] = relationship('ChainTransfer', back_populates='ct_angle_list')


class CtAxisAngleList(Base):
    __tablename__ = 'ct_axis_angle_list'
    __table_args__ = (
        ForeignKeyConstraint(['ct_id'], ['protdb3.chain_transfer.ct_id'], name='ct_axis_angle_list_fk'),
        PrimaryKeyConstraint('ctal_id', name='ct_axis_angle_list_pk'),
        {'schema': 'protdb3'}
    )

    ctal_id: Mapped[int] = mapped_column(Integer, Sequence('ctalid_seq', schema='protdb3'), primary_key=True)
    ct_id: Mapped[Optional[int]] = mapped_column(Integer)
    ct_axis_angle: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    ct: Mapped[Optional['ChainTransfer']] = relationship('ChainTransfer', back_populates='ct_axis_angle_list')


class CtCoverageList(Base):
    __tablename__ = 'ct_coverage_list'
    __table_args__ = (
        ForeignKeyConstraint(['ct_id'], ['protdb3.chain_transfer.ct_id'], name='ct_coverage_list_fk'),
        PrimaryKeyConstraint('ccl_id', name='ct_coverage_list_pk'),
        {'schema': 'protdb3'}
    )

    ccl_id: Mapped[int] = mapped_column(Integer, Sequence('cclid_seq', schema='protdb3'), primary_key=True)
    ct_id: Mapped[Optional[int]] = mapped_column(Integer)
    ct_coverage: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    ct: Mapped[Optional['ChainTransfer']] = relationship('ChainTransfer', back_populates='ct_coverage_list')


class CtLevelsList(Base):
    __tablename__ = 'ct_levels_list'
    __table_args__ = (
        ForeignKeyConstraint(['ct_id'], ['protdb3.chain_transfer.ct_id'], name='ct_levels_list_fk'),
        PrimaryKeyConstraint('cls_id', name='ct_levels_list_pk'),
        {'schema': 'protdb3'}
    )

    cls_id: Mapped[int] = mapped_column(Integer, Sequence('clsid_seq', schema='protdb3'), primary_key=True)
    ct_id: Mapped[Optional[int]] = mapped_column(Integer)
    ct_levels: Mapped[Optional[int]] = mapped_column(Integer, server_default=text('0.0'))

    ct: Mapped[Optional['ChainTransfer']] = relationship('ChainTransfer', back_populates='ct_levels_list')


class CtNumRepeatsList(Base):
    __tablename__ = 'ct_num_repeats_list'
    __table_args__ = (
        ForeignKeyConstraint(['ct_id'], ['protdb3.chain_transfer.ct_id'], name='ct_num_repeats_list_fk'),
        PrimaryKeyConstraint('cnr_id', name='ct_num_repeats_list_pk'),
        {'schema': 'protdb3'}
    )

    cnr_id: Mapped[int] = mapped_column(Integer, Sequence('cnrid_seq', schema='protdb3'), primary_key=True)
    ct_id: Mapped[Optional[int]] = mapped_column(Integer)
    ct_num_repeats: Mapped[Optional[int]] = mapped_column(Integer, server_default=text('0.0'))

    ct: Mapped[Optional['ChainTransfer']] = relationship('ChainTransfer', back_populates='ct_num_repeats_list')


class CtRepLengthList(Base):
    __tablename__ = 'ct_rep_length_list'
    __table_args__ = (
        ForeignKeyConstraint(['ct_id'], ['protdb3.chain_transfer.ct_id'], name='ct_rep_length_list_fk'),
        PrimaryKeyConstraint('crll_id', name='ct_rep_length_list_pk'),
        {'schema': 'protdb3'}
    )

    crll_id: Mapped[int] = mapped_column(Integer, Sequence('crllid_seq', schema='protdb3'), primary_key=True)
    ct_id: Mapped[Optional[int]] = mapped_column(Integer)
    ct_rep_length: Mapped[Optional[int]] = mapped_column(Integer, server_default=text('0.0'))

    ct: Mapped[Optional['ChainTransfer']] = relationship('ChainTransfer', back_populates='ct_rep_length_list')


class CtRmsdList(Base):
    __tablename__ = 'ct_rmsd_list'
    __table_args__ = (
        ForeignKeyConstraint(['ct_id'], ['protdb3.chain_transfer.ct_id'], name='ct_rmsd_list_fk'),
        PrimaryKeyConstraint('crl_id', name='ct_rmsd_list_pk'),
        {'schema': 'protdb3'}
    )

    crl_id: Mapped[int] = mapped_column(Integer, Sequence('crlid_seq', schema='protdb3'), primary_key=True)
    ct_id: Mapped[Optional[int]] = mapped_column(Integer)
    ct_rmsd: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    ct: Mapped[Optional['ChainTransfer']] = relationship('ChainTransfer', back_populates='ct_rmsd_list')


class CtTmscoreList(Base):
    __tablename__ = 'ct_tmscore_list'
    __table_args__ = (
        ForeignKeyConstraint(['ct_id'], ['protdb3.chain_transfer.ct_id'], name='ct_tmscore_list_fk'),
        PrimaryKeyConstraint('csl_id', name='ct_tmscore_list_pk'),
        {'schema': 'protdb3'}
    )

    csl_id: Mapped[int] = mapped_column(Integer, Sequence('cslid_seq', schema='protdb3'), primary_key=True)
    ct_id: Mapped[Optional[int]] = mapped_column(Integer)
    ct_tmscore: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    ct: Mapped[Optional['ChainTransfer']] = relationship('ChainTransfer', back_populates='ct_tmscore_list')


class CtTranslationList(Base):
    __tablename__ = 'ct_translation_list'
    __table_args__ = (
        ForeignKeyConstraint(['ct_id'], ['protdb3.chain_transfer.ct_id'], name='ct_translation_list_fk'),
        PrimaryKeyConstraint('ctt_id', name='ct_translation_list_pk'),
        {'schema': 'protdb3'}
    )

    ctt_id: Mapped[int] = mapped_column(Integer, Sequence('cttid_seq', schema='protdb3'), primary_key=True)
    ct_id: Mapped[Optional[int]] = mapped_column(Integer)
    ct_translation: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    ct: Mapped[Optional['ChainTransfer']] = relationship('ChainTransfer', back_populates='ct_translation_list')


class PcscSymmAlignedResidueList(Base):
    __tablename__ = 'pcsc_symm_aligned_residue_list'
    __table_args__ = (
        ForeignKeyConstraint(['csc_id'], ['protdb3.pdb_ce_symm_classification.csc_id'], name='pcsc_symm_aligned_residue_list_fk'),
        PrimaryKeyConstraint('psarl_id', name='pcsc_symm_aligned_residue_list_pk'),
        {'schema': 'protdb3'}
    )

    psarl_id: Mapped[int] = mapped_column(Integer, Sequence('psarlid_seq', schema='protdb3'), primary_key=True)
    csc_id: Mapped[Optional[int]] = mapped_column(Integer)
    ccsc_symm_aligned_residue: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    csc: Mapped[Optional['PdbCeSymmClassification']] = relationship('PdbCeSymmClassification', back_populates='pcsc_symm_aligned_residue_list')


class PcscSymmAngleList(Base):
    __tablename__ = 'pcsc_symm_angle_list'
    __table_args__ = (
        ForeignKeyConstraint(['csc_id'], ['protdb3.pdb_ce_symm_classification.csc_id'], name='pcsc_symm_angle_list_fk'),
        PrimaryKeyConstraint('sal_id', name='pcsc_symm_angle_list_pk'),
        {'schema': 'protdb3'}
    )

    sal_id: Mapped[int] = mapped_column(Integer, Sequence('salid_seq', schema='protdb3'), primary_key=True)
    csc_id: Mapped[Optional[int]] = mapped_column(Integer)
    symm_angle: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    csc: Mapped[Optional['PdbCeSymmClassification']] = relationship('PdbCeSymmClassification', back_populates='pcsc_symm_angle_list')


class PcscSymmSeedList(Base):
    __tablename__ = 'pcsc_symm_seed_list'
    __table_args__ = (
        ForeignKeyConstraint(['csc_id'], ['protdb3.pdb_ce_symm_classification.csc_id'], name='pcsc_symm_seed_list_fk'),
        PrimaryKeyConstraint('pssl_id', name='pcsc_symm_seed_list_pk'),
        {'schema': 'protdb3'}
    )

    pssl_id: Mapped[int] = mapped_column(Integer, Sequence('psslid_seq', schema='protdb3'), primary_key=True)
    csc_id: Mapped[Optional[int]] = mapped_column(Integer)
    pcsc_symm_seed: Mapped[Optional[int]] = mapped_column(Integer)

    csc: Mapped[Optional['PdbCeSymmClassification']] = relationship('PdbCeSymmClassification', back_populates='pcsc_symm_seed_list')


class PcscSymmTranslationList(Base):
    __tablename__ = 'pcsc_symm_translation_list'
    __table_args__ = (
        ForeignKeyConstraint(['csc_id'], ['protdb3.pdb_ce_symm_classification.csc_id'], name='pcsc_symm_translation_list_fk'),
        PrimaryKeyConstraint('sil_id', name='pcsc_symm_translation_list_pk'),
        {'schema': 'protdb3'}
    )

    sil_id: Mapped[int] = mapped_column(Integer, Sequence('silid_seq', schema='protdb3'), primary_key=True)
    csc_id: Mapped[Optional[int]] = mapped_column(Integer)
    symm_translation: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    csc: Mapped[Optional['PdbCeSymmClassification']] = relationship('PdbCeSymmClassification', back_populates='pcsc_symm_translation_list')


class PcscSymmUnrefinedRmsdList(Base):
    __tablename__ = 'pcsc_symm_unrefined_rmsd_list'
    __table_args__ = (
        ForeignKeyConstraint(['csc_id'], ['protdb3.pdb_ce_symm_classification.csc_id'], name='pcsc_symm_unrefined_rmsd_list_fk'),
        PrimaryKeyConstraint('psurl_id', name='pcsc_symm_unrefined_rmsd_list_pk'),
        {'schema': 'protdb3'}
    )

    psurl_id: Mapped[int] = mapped_column(Integer, Sequence('psurlid_seq', schema='protdb3'), primary_key=True)
    csc_id: Mapped[Optional[int]] = mapped_column(Integer)
    pcsc_symm_unrefined_rmsd: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    csc: Mapped[Optional['PdbCeSymmClassification']] = relationship('PdbCeSymmClassification', back_populates='pcsc_symm_unrefined_rmsd_list')


class PcscSymmUnrefinedTmScoreList(Base):
    __tablename__ = 'pcsc_symm_unrefined_tm_score_list'
    __table_args__ = (
        ForeignKeyConstraint(['csc_id'], ['protdb3.pdb_ce_symm_classification.csc_id'], name='pcsc_symm_unrefined_tm_score_list_fk'),
        PrimaryKeyConstraint('psutsl_id', name='pcsc_symm_unrefined_tm_score_list_pk'),
        {'schema': 'protdb3'}
    )

    psutsl_id: Mapped[int] = mapped_column(Integer, Sequence('psutslid_seq', schema='protdb3'), primary_key=True)
    csc_id: Mapped[Optional[int]] = mapped_column(Integer)
    pcsc_symm_unrefined_tm_score: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    csc: Mapped[Optional['PdbCeSymmClassification']] = relationship('PdbCeSymmClassification', back_populates='pcsc_symm_unrefined_tm_score_list')


class PdbAnalysisAlignments(Base):
    __tablename__ = 'pdb_analysis_alignments'
    __table_args__ = (
        ForeignKeyConstraint(['psa_id'], ['protdb3.pdb_symmetry_analysis.psa_id'], name='pdb_analysis_alignments_fk'),
        PrimaryKeyConstraint('paa_id', name='pdb_analysis_alignments_pk'),
        {'schema': 'protdb3'}
    )

    paa_id: Mapped[int] = mapped_column(Integer, Sequence('paaid_seq', schema='protdb3'), primary_key=True)
    aa_name: Mapped[Optional[str]] = mapped_column(String(1500))
    aa_seq: Mapped[Optional[str]] = mapped_column(String(4000))
    psa_id: Mapped[Optional[int]] = mapped_column(Integer)
    aa_seq2: Mapped[Optional[str]] = mapped_column(String(4000))

    psa: Mapped[Optional['PdbSymmetryAnalysis']] = relationship('PdbSymmetryAnalysis', back_populates='pdb_analysis_alignments')


class PdbAnalysisDownloads(Base):
    __tablename__ = 'pdb_analysis_downloads'
    __table_args__ = (
        ForeignKeyConstraint(['psa_id'], ['protdb3.pdb_symmetry_analysis.psa_id'], name='pdb_analysis_downloads_fk'),
        PrimaryKeyConstraint('pad_id', name='pdb_analysis_downloads_pk'),
        {'schema': 'protdb3'}
    )

    pad_id: Mapped[int] = mapped_column(Integer, Sequence('padid_seq', schema='protdb3'), primary_key=True)
    psa_id: Mapped[Optional[int]] = mapped_column(Integer)
    file_path: Mapped[Optional[str]] = mapped_column(String(255))

    psa: Mapped[Optional['PdbSymmetryAnalysis']] = relationship('PdbSymmetryAnalysis', back_populates='pdb_analysis_downloads')


class PdbAnalysisImages(Base):
    __tablename__ = 'pdb_analysis_images'
    __table_args__ = (
        ForeignKeyConstraint(['psa_id'], ['protdb3.pdb_symmetry_analysis.psa_id'], name='pdb_analysis_images_fk'),
        PrimaryKeyConstraint('pai_id', name='pdb_analysis_images_pk'),
        {'schema': 'protdb3'}
    )

    pai_id: Mapped[int] = mapped_column(Integer, Sequence('paiid_seq', schema='protdb3'), primary_key=True)
    psa_id: Mapped[Optional[int]] = mapped_column(Integer)
    image_file_and_path: Mapped[Optional[str]] = mapped_column(String(1000))
    pml_file_and_path: Mapped[Optional[str]] = mapped_column(String(1000))
    json_file_path_3dmol: Mapped[Optional[str]] = mapped_column(String(1000))
    super_pml_file_and_path: Mapped[Optional[str]] = mapped_column(String(1000))

    psa: Mapped[Optional['PdbSymmetryAnalysis']] = relationship('PdbSymmetryAnalysis', back_populates='pdb_analysis_images')


class PdbNtmChainDomains(Base):
    __tablename__ = 'pdb_ntm_chain_domains'
    __table_args__ = (
        ForeignKeyConstraint(['ncd_id'], ['protdb3.pdb_ntm_chains_details.ncd_id'], name='pdb_ntm_chain_domains_fk'),
        PrimaryKeyConstraint('ntd_id', name='pdb_ntm_chain_domains_pk'),
        {'schema': 'protdb3'}
    )

    ntd_id: Mapped[int] = mapped_column(Integer, Sequence('ntdid_seq', schema='protdb3'), primary_key=True)
    tm_domain_id: Mapped[Optional[int]] = mapped_column(Integer)
    tmd_range: Mapped[Optional[str]] = mapped_column(String)
    ncd_id: Mapped[Optional[int]] = mapped_column(Integer)

    ncd: Mapped[Optional['PdbNtmChainsDetails']] = relationship('PdbNtmChainsDetails', back_populates='pdb_ntm_chain_domains')


class PdbSymdAlignments(Base):
    __tablename__ = 'pdb_symd_alignments'
    __table_args__ = (
        ForeignKeyConstraint(['csd_id'], ['protdb3.pdb_ce_symd_classification.csd_id'], name='pdb_symd_alignments_fk'),
        PrimaryKeyConstraint('sda_id', name='pdb_symd_alignments_pk'),
        {'schema': 'protdb3'}
    )

    sda_id: Mapped[int] = mapped_column(Integer, Sequence('sdaid_seq', schema='protdb3'), primary_key=True)
    sa_name: Mapped[Optional[str]] = mapped_column(String(255))
    sa_seq: Mapped[Optional[str]] = mapped_column(String(4000))
    csd_id: Mapped[Optional[int]] = mapped_column(Integer)
    sa_seq2: Mapped[Optional[str]] = mapped_column(String(4000))
    sa_seq3: Mapped[Optional[str]] = mapped_column(String(4000))
    sa_seq4: Mapped[Optional[str]] = mapped_column(String(4000))

    csd: Mapped[Optional['PdbCeSymdClassification']] = relationship('PdbCeSymdClassification', back_populates='pdb_symd_alignments')


class PdbSymdalignDownloads(Base):
    __tablename__ = 'pdb_symdalign_downloads'
    __table_args__ = (
        ForeignKeyConstraint(['csd_id'], ['protdb3.pdb_ce_symd_classification.csd_id'], name='pdb_symdalign_downloads_fk'),
        PrimaryKeyConstraint('ssdd_id', name='pdb_symdalign_downloads_pk'),
        {'schema': 'protdb3'}
    )

    ssdd_id: Mapped[int] = mapped_column(Integer, Sequence('ssddid_seq', schema='protdb3'), primary_key=True)
    csd_id: Mapped[Optional[int]] = mapped_column(Integer)
    file_path: Mapped[Optional[str]] = mapped_column(String(255))

    csd: Mapped[Optional['PdbCeSymdClassification']] = relationship('PdbCeSymdClassification', back_populates='pdb_symdalign_downloads')


class PdbSymmAlignments(Base):
    __tablename__ = 'pdb_symm_alignments'
    __table_args__ = (
        ForeignKeyConstraint(['csc_id'], ['protdb3.pdb_ce_symm_classification.csc_id'], name='pdb_symm_alignments_fk'),
        PrimaryKeyConstraint('sma_id', name='pdb_symm_alignments_pk'),
        {'schema': 'protdb3'}
    )

    sma_id: Mapped[int] = mapped_column(Integer, Sequence('smaid_seq', schema='protdb3'), primary_key=True)
    sa_name: Mapped[Optional[str]] = mapped_column(String(255))
    sa_seq: Mapped[Optional[str]] = mapped_column(String(4000))
    csc_id: Mapped[Optional[int]] = mapped_column(Integer)
    sa_seq2: Mapped[Optional[str]] = mapped_column(String(4000))

    csc: Mapped[Optional['PdbCeSymmClassification']] = relationship('PdbCeSymmClassification', back_populates='pdb_symm_alignments')


class PdbSymmalignDownloads(Base):
    __tablename__ = 'pdb_symmalign_downloads'
    __table_args__ = (
        ForeignKeyConstraint(['csc_id'], ['protdb3.pdb_ce_symm_classification.csc_id'], name='pdb_symmalign_downloads_fk'),
        PrimaryKeyConstraint('ssad_id', name='pdb_symmalign_downloads_pk'),
        {'schema': 'protdb3'}
    )

    ssad_id: Mapped[int] = mapped_column(Integer, Sequence('ssadid_seq', schema='protdb3'), primary_key=True)
    csc_id: Mapped[Optional[int]] = mapped_column(Integer)
    file_path: Mapped[Optional[str]] = mapped_column(String(255))

    csc: Mapped[Optional['PdbCeSymmClassification']] = relationship('PdbCeSymmClassification', back_populates='pdb_symmalign_downloads')


class PsaAlignedResiduesList(Base):
    __tablename__ = 'psa_aligned_residues_list'
    __table_args__ = (
        ForeignKeyConstraint(['psa_id'], ['protdb3.pdb_symmetry_analysis.psa_id'], name='psa_aligned_residues_list_fk'),
        PrimaryKeyConstraint('parl_id', name='psa_aligned_residues_list_pk'),
        {'schema': 'protdb3'}
    )

    parl_id: Mapped[int] = mapped_column(Integer, Sequence('parlid_seq', schema='protdb3'), primary_key=True)
    psa_id: Mapped[Optional[int]] = mapped_column(Integer)
    psa_aligned_residue: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    psa: Mapped[Optional['PdbSymmetryAnalysis']] = relationship('PdbSymmetryAnalysis', back_populates='psa_aligned_residues_list')


class PsaAngleList(Base):
    __tablename__ = 'psa_angle_list'
    __table_args__ = (
        ForeignKeyConstraint(['psa_id'], ['protdb3.pdb_symmetry_analysis.psa_id'], name='psa_angle_list_fk'),
        PrimaryKeyConstraint('pal_id', name='psa_angle_list_pk'),
        {'schema': 'protdb3'}
    )

    pal_id: Mapped[int] = mapped_column(Integer, Sequence('palid_seq', schema='protdb3'), primary_key=True)
    psa_id: Mapped[int] = mapped_column(Integer, nullable=False)
    psa_angle: Mapped[Optional[float]] = mapped_column(Double(53))

    psa: Mapped['PdbSymmetryAnalysis'] = relationship('PdbSymmetryAnalysis', back_populates='psa_angle_list')


class PsaAxisAngleList(Base):
    __tablename__ = 'psa_axis_angle_list'
    __table_args__ = (
        ForeignKeyConstraint(['psa_id'], ['protdb3.pdb_symmetry_analysis.psa_id'], name='psa_axis_angle_list_fk'),
        PrimaryKeyConstraint('paal_id', name='psa_axis_angle_list_pk'),
        {'schema': 'protdb3'}
    )

    paal_id: Mapped[int] = mapped_column(Integer, Sequence('paalid_seq', schema='protdb3'), primary_key=True)
    psa_id: Mapped[int] = mapped_column(Integer, nullable=False)
    psa_axis_angle: Mapped[float] = mapped_column(Double(53), nullable=False)

    psa: Mapped['PdbSymmetryAnalysis'] = relationship('PdbSymmetryAnalysis', back_populates='psa_axis_angle_list')


class PsaRmsdList(Base):
    __tablename__ = 'psa_rmsd_list'
    __table_args__ = (
        ForeignKeyConstraint(['psa_id'], ['protdb3.pdb_symmetry_analysis.psa_id'], name='psa_rmsd_list_fk'),
        PrimaryKeyConstraint('psrl_id', name='psa_rmsd_list_pk'),
        {'schema': 'protdb3'}
    )

    psrl_id: Mapped[int] = mapped_column(Integer, Sequence('psrlid_seq', schema='protdb3'), primary_key=True)
    psa_id: Mapped[Optional[int]] = mapped_column(Integer)
    psa_rmsd: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    psa: Mapped[Optional['PdbSymmetryAnalysis']] = relationship('PdbSymmetryAnalysis', back_populates='psa_rmsd_list')


class PsaTmscoreList(Base):
    __tablename__ = 'psa_tmscore_list'
    __table_args__ = (
        ForeignKeyConstraint(['psa_id'], ['protdb3.pdb_symmetry_analysis.psa_id'], name='psa_tmscore_list_fk'),
        PrimaryKeyConstraint('pstl_id', name='psa_tmscore_list_pk'),
        {'schema': 'protdb3'}
    )

    pstl_id: Mapped[int] = mapped_column(Integer, Sequence('pstlid_seq', schema='protdb3'), primary_key=True)
    psa_id: Mapped[Optional[int]] = mapped_column(Integer)
    psa_tmscore: Mapped[Optional[float]] = mapped_column(Double(53), server_default=text('0.0'))

    psa: Mapped[Optional['PdbSymmetryAnalysis']] = relationship('PdbSymmetryAnalysis', back_populates='psa_tmscore_list')


class PsaTranslationList(Base):
    __tablename__ = 'psa_translation_list'
    __table_args__ = (
        ForeignKeyConstraint(['psa_id'], ['protdb3.pdb_symmetry_analysis.psa_id'], name='psa_translation_list_fk'),
        PrimaryKeyConstraint('ptl_id', name='psa_translation_list_pk'),
        {'schema': 'protdb3'}
    )

    ptl_id: Mapped[int] = mapped_column(Integer, Sequence('ptlid_seq', schema='protdb3'), primary_key=True)
    psa_id: Mapped[int] = mapped_column(Integer, nullable=False)
    psa_translation: Mapped[Optional[float]] = mapped_column(Double(53))

    psa: Mapped['PdbSymmetryAnalysis'] = relationship('PdbSymmetryAnalysis', back_populates='psa_translation_list')
