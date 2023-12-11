# Name: cesymm_symd_pngs.py
# Author: Antoniya A. Aleksandrova
# Language: Python 3.5+
# Date: 16 Jun 2022

"""
Create the images and dictionary with the data from CE-Symm and SymD that should be
displayed in the first symmetry section of the website EncoMPASS

Usage: python cesymm_symd_pngs.py <locations> <pdb code> <pdb_suff> <run_id>
"""

import sys
import os
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from symmetry_exec_functions import *
from initialize_repository import initialize_repository

def dic_for_cesymm_symd(locations, pdb_suff, inputf, ray=1, run_id = "0000"):
    """
    Produce images and data dictionary for symd and ce-symm default results

    :param locations: dictionary with all essential paths
    :param pdb_suff: suffix attached to the pdb code (e.g. "_enc" in "1okc_enc.pdb")
    :param inputf: pdb code
    :param ray: turn on/off (1/0) PyMOL ray tracing
    :return:
    """
    if not os.path.isdir(locations['FSYSPATH']['images_cesymm_symd_log'] + inputf):
        os.mkdir(locations['FSYSPATH']['images_cesymm_symd_log'] + inputf)
    log_file = open(locations['FSYSPATH']['images_cesymm_symd_log'] + inputf + "/" + inputf + '_run.log', 'w')
    save_out = sys.stdout
    save_err = sys.stderr
    sys.stdout = sys.stderr = log_file

    print("Run id: %s" % run_id)
    if os.path.isfile(locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb')==False:
        raise SystemExit("{0} does not exist.".format(locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb'))
    oriented_opm=locations['FSYSPATH']['whole']+ inputf + pdb_suff + '.pdb'
    tm_archive = locations['FSYSPATH']['symmetry'] + '.tm_archive.pkl'
    dir = locations['FSYSPATH']['symmtemp'] + "webarchive/"
    symd_out_dir=locations['FSYSPATH']['symd'] + inputf + "/"
    cesymm_out_dir=locations['FSYSPATH']['cesymm'] + inputf + "/"
    pymol = "pymol -c" # path to PyMOL with command line option
    pic_dir_whole=locations['FSYSPATH']['symm_wholepngs']
    pic_dir_chain=locations['FSYSPATH']['symm_chainspngs']
    print(locations['FSYSPATH']['symm_wholejsons'])
    jmol_dir_whole=locations['FSYSPATH']['symm_wholejsons']
    jmol_dir_chain=locations['FSYSPATH']['symm_chainsjsons']
    cs_files_path=locations['FSYS']['cesymm']
    symd_files_path=locations['FSYS']['symd']

    # Find TM chains
    a = pkl.load(open(tm_archive, 'rb'))
    if len(a[inputf]['tmchains']) > 0:
        chains = sorted(a[inputf]['tmchains'])
        print(chains)

    else:
        raise SystemExit("{0} has no TM chains.".format(locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb'))


    # Create a dictionary of all the CE-Symm and SymD symmetry features (does not require a loaded structure)
    symmlist = {}
    symmlist[inputf]={}
    cesymm_dic = cesymm_data(cesymm_out_dir,inputf)
    symmlist[inputf]['cesymm']={}
    symmlist[inputf]['cesymm']=cesymm_dic
    symd_dic=symd_data(symd_out_dir,inputf)
    symmlist[inputf]['symd']={}
    symmlist[inputf]['symd']=symd_dic

    # Create CE-Symm symmetry image

    num=strip_tm_chains(dir,inputf,oriented_opm,chains)
    tmp_ref=parse_structure(dir+inputf+"_tmp.pdb")[0]
    image_key=inputf+"_cesymm"
    try:
        if int(symmlist[inputf]['cesymm']['symmetry_levels'])>0 and symmlist[inputf]['cesymm']['repeats_number']!='na' and int(symmlist[inputf]['cesymm']['repeats_number'])>1:
          repeat_selection=get_repeat_resid(cesymm_out_dir,inputf,tmp_ref)
          image_files, pml_files, jmol_files=cesymm_images(dir,cesymm_out_dir,inputf,repeat_selection,oriented_opm,tmp_ref,ray,pymol,pic_dir_whole,jmol_dir_whole,image_key)
          symmlist[inputf]['cesymm']['image_files'] = image_files
          symmlist[inputf]['cesymm']['pml_files'] = pml_files
          symmlist[inputf]['cesymm']['jmol_files'] = jmol_files

    except:
      print(sys.exc_info()[0])
      symmlist[inputf]['cesymm']['image_files'] = ""
      symmlist[inputf]['cesymm']['pml_files'] = ""
      symmlist[inputf]['cesymm']['jmol_files'] = ""


    # Create SymD symmetry image
    if os.path.isfile(symd_out_dir+inputf+"-trfm.pdb"):
        aligned_ids, aligned_chains=get_symd_aligned_resid(symd_out_dir,inputf,tmp_ref)
        symd_images(dir,symd_out_dir,inputf,oriented_opm,ray,pymol,pic_dir_whole,jmol_dir_whole,aligned_ids,aligned_chains)
        symmlist[inputf]['symd']['image_files'] = locations['FSYS']['symm_wholepngs']+inputf+'_symd.png'
        symmlist[inputf]['symd']['pml_files'] = locations['FSYS']['symm_wholepngs']+'pymol_script_'+inputf+'_symd.pml'
        symmlist[inputf]['symd']['jmol_symd'] = locations['FSYS']['symm_wholejsons']+'3dmol_'+inputf+'_symd.json'
    else:
        symmlist[inputf]['symd']['image_files'] = ""
        symmlist[inputf]['symd']['pml_files'] = ""
        symmlist[inputf]['symd']['jmol_symd'] = ""

    os.remove(dir+inputf+"_tmp.pdb")
    # Start writing to file for whole PDB
    try:
      symmlist[inputf]['cesymm']['alignment'] = cesymm_alignment(cesymm_out_dir,inputf)
    except:
      symmlist[inputf]['cesymm']['alignment'] = ""

    try:
      symmlist[inputf]['symd']['alignment'] = symd_alignment(symd_out_dir,inputf)
    except:
      symmlist[inputf]['symd']['alignment']= ""


    ###### Deal with chains ######
    pdb = inputf
    if len(chains)>1:
        for chain in chains:
            inputf = pdb + "_" + chain
            num=strip_tm_chains(dir,inputf,oriented_opm,chain)
            oriented=dir + inputf + "_tmp.pdb"
            tmp_ref=parse_structure(dir + inputf + "_tmp.pdb")[0]
            symmlist[inputf]={}
            cesymm_dic = cesymm_data(cesymm_out_dir,inputf)
            symmlist[inputf]['cesymm'] = {}
            symmlist[inputf]['cesymm'] = cesymm_dic

            symd_dic=symd_data(symd_out_dir,inputf)
            symmlist[inputf]['symd'] = {}
            symmlist[inputf]['symd'] = symd_dic

            # Create CE-Symm symmetry image
            image_key=inputf + "_cesymm"
            try:
              if int(symmlist[inputf]['cesymm']['symmetry_levels'])>0 and symmlist[inputf]['cesymm']['repeats_number']!='na' and int(symmlist[inputf]['cesymm']['repeats_number'])>1:
                repeat_selection=get_repeat_resid(cesymm_out_dir,inputf,tmp_ref)
                image_files, pml_files, jmol_files=cesymm_images(dir,cesymm_out_dir,inputf,repeat_selection,oriented,tmp_ref,ray,pymol,pic_dir_chain,jmol_dir_chain,image_key)
                symmlist[inputf]['cesymm']['image_files'] = image_files
                symmlist[inputf]['cesymm']['pml_files'] = pml_files
                symmlist[inputf]['cesymm']['jmol_files'] = jmol_files
            except:
                symmlist[inputf]['cesymm']['image_files'] = ""
                symmlist[inputf]['cesymm']['pml_files'] = ""
                symmlist[inputf]['cesymm']['jmol_files'] = ""
            # Create SymD symmetry image
            if os.path.isfile(symd_out_dir+inputf+"-trfm.pdb"):
                aligned_ids, aligned_chains=get_symd_aligned_resid(symd_out_dir,inputf,tmp_ref)
                symd_images(dir,symd_out_dir,inputf,oriented,ray,pymol,pic_dir_chain,jmol_dir_chain,aligned_ids,aligned_chains)
                symmlist[inputf]['symd']['image_files'] = locations['FSYS']['symm_chainspngs'] + inputf + "_symd.png"
                symmlist[inputf]['symd']['pml_files'] = locations['FSYS']['symm_chainspngs']+'pymol_script_'+inputf+'_symd.pml'
                symmlist[inputf]['symd']['jmol_symd'] = locations['FSYS']['symm_chainsjsons']+'3dmol_'+inputf+'_symd.json'
            else:
                symmlist[inputf]['symd']['image_files'] = ""
                symmlist[inputf]['symd']['pml_files'] = ""
                symmlist[inputf]['symd']['jmol_symd'] = ""
            os.remove(dir+inputf+"_tmp.pdb")

            try:
                symmlist[inputf]['cesymm']['alignment']=cesymm_alignment(cesymm_out_dir,inputf+"_"+chain)
            except:
                symmlist[inputf]['cesymm']['alignment'] =''

            try:
                symmlist[inputf]['symd']['alignment']=symd_alignment(symd_out_dir,inputf+'_'+chain)
            except:
                symmlist[inputf]['symd']['alignment'] = ""

    pkl.dump(symmlist, open(locations['FSYSPATH']['images_cesymm_symd_log'] + inputf[0:4] + "/" + inputf[0:4] + "_dic.pkl", "wb"))
    print(inputf[0:4] + "  data and image generation completed.")
    print("Time[s]:%s" % (time.time() - start_time))
    print("Date: %s" % time.asctime())

    sys.stdout = save_out
    sys.stderr = save_err
    log_file.close()
    return symmlist

if __name__ == "__main__":
    if len(sys.argv) < 5:
        raise SystemExit("syntax: %s <locations_path> <pdbcode> <pdb_suff> <run_id>" % sys.argv[0])
    else:
        locations_path = sys.argv[1].strip()
        inputf = sys.argv[2][0:4]
        pdb_suff = sys.argv[3].strip()
        run_id = sys.argv[4].strip()

    options, locations = initialize_repository(main_path=locations_path, instr_filename_path="EncoMPASS_options_relative_to__instruction_file_benchmark.txt")
    dir=locations['FSYSPATH']['symmtemp']+"webarchive/"
    if os.path.isdir(dir)==False:
        os.mkdir(dir)
    print(dic_for_cesymm_symd(locations, pdb_suff, inputf, run_id = run_id))