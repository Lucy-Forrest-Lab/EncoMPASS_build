# Name: select_quat_symm_minlen.py
# Author: Antoniya A. Aleksandrova
# Language: Python 3.5+
# Description: Finds results for all protein complexes with more than 1 chain, which have multiple levels that are not all with quaternary symmetry, and notes them down for a later run with CE-Symm and restricted levels option


import sys
import os
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from initialize_repository import initialize_repository
from symmetry_exec_functions import *

def select_quat_symm_minlen(locations, pdb_suff, tm_dic=False):
    tm_archive = locations['FSYSPATH']['symmetry'] + '.tm_archive.pkl'
    if tm_dic==False:
        tm_dic = pkl.load(open(tm_archive, 'rb'))
    chain_ordered_dir=locations['FSYSPATH']['symmtemp'] +"chain_ordered_str/"
    wkdir=locations['FSYSPATH']['symmtemp'] + "cesymm_quat_minlen/"
    if os.path.isdir(wkdir)==False:
        os.mkdir(wkdir)
    out1=open(locations['FSYSPATH']['cesymm_quat_minlen_log']+"quat_symm_minlen_small_pdbs.log","w")
    out2=open(locations['FSYSPATH']['cesymm_quat_minlen_log']+"quat_symm_minlen_big_pdbs.log","w")
    out1.write("#Input\tCE-Symm Seed\tCE-Symm Quat. Symmetry Levels\tCE-Symm Repeat Length\n")
    out2.write("#Input\tCE-Symm Seed\tCE-Symm Quat. Symmetry Levels\tCE-Symm Repeat Length\n")
    select_quat_dic = {}

    tm_dic = pkl.load(open(tm_archive, 'rb'))
    for inputf in sorted(list(tm_dic)):
        print(inputf)
        cesymm_out_dir = locations['FSYSPATH']['cesymm'] + inputf + "/"
        cesymm_lt_out_dir = locations['FSYSPATH']['cesymm_low_thr']  + inputf + "/"
        if os.path.isfile(locations['FSYSPATH']['whole']+inputf+pdb_suff+'.pdb')==False:
            continue
        chains=sorted(tm_dic[inputf]['tmchains'])
        type=tm_dic[inputf]['class']
        if type!='alpha' or len(chains)<2:
            continue
        symmlist = {}
        symmlist[inputf]={}
        # Whole pdb
        oriented_opm=locations['FSYSPATH']['whole']+inputf+pdb_suff+'.pdb'
        limit_bottom,limit_top=membrane_limits(oriented_opm)
        num=strip_tm_chains(wkdir,inputf,oriented_opm,chains) #should not strip them in order because they were not run that way
        if num<20000:
            out=out1
        else:
            out=out2

        all_suffixes=['_com_sd','_ordered','_sd_order', '', '_com_sd_lt','_ordered_lt','_sd_order_lt', '_lt']
        for suff in all_suffixes:
            if '_lt' in suff:
                cs_dir=cesymm_lt_out_dir
                name = suff[:-3]
                tag = 'lt'
            else:
                cs_dir=cesymm_out_dir
                name = suff
                tag = 'df'
            if os.path.isfile(cs_dir+inputf+name+"_stdout.out"):
                cesymm_dic = cesymm_data(cs_dir,inputf+name)
                symmlist[inputf]['cesymm']={}
                symmlist[inputf]['cesymm']=cesymm_dic
                if name=='':
                    oriented_tmp=wkdir+inputf+"_tmp.pdb"
                else:
                    oriented_tmp=chain_ordered_dir+inputf+name+".pdb"
                if symmlist[inputf]['cesymm']['repeats_number']=='na':
                    pass
                elif int(symmlist[inputf]['cesymm']['repeats_number']) > 1:
                    reference = parse_structure(oriented_tmp)[0]
                    repeat_selection = get_repeat_resid(cs_dir,inputf+name,reference)
                    symm_location = is_symmetry_in_membrane(limit_bottom,limit_top,repeat_selection,oriented_tmp)
                    doms, crossings = symmetry_tm_domains(tm_archive, repeat_selection, inputf, '')
                    if int(symmlist[inputf]['cesymm']['symmetry_levels']) > 1 and (int(symmlist[inputf]['cesymm']['repeat_length'])>=40 and (doms=='in' and symm_location=='in')):
                        levels, level_type, axes_per_level, order_of_level, level_symm_type=cesymm_levels(inputf+name,cs_dir)
                        if 'Quaternary' in level_symm_type:
                            num_levels=level_symm_type.count('Quaternary')
                            out.write(inputf+name+"\t"+tag+"\t"+symmlist[inputf]['cesymm']['seed']+"\t"+str(num_levels)+"\t"+str(symmlist[inputf]['cesymm']['repeat_length'])+"\n")
                            select_quat_dic[inputf + name + '_' + tag] = {'name': inputf + name + '_' + tag,
                                                                          'mode' : tag,
                                                                          'seed': symmlist[inputf]['cesymm']['seed'],
                                                                          'required_symm_levels': num_levels,
                                                                          'cesymm_rep_length': symmlist[inputf]['cesymm']['repeat_length'],
                                                                          'natoms': tm_dic[inputf]['natoms'],
                                                                          'nres': tm_dic[inputf]['nres']}
                        if 'Quaternary' in level_symm_type and 'Internal' in level_symm_type:
                            print(inputf+' - internal symmetry can be extracted here\n')



        all_names=['_com_sd_df','_ordered_df','_sd_order_df', '_df', '_com_sd_lt','_ordered_lt','_sd_order_lt', '_lt']
        for name in all_names:
            cesymm_minlen_out_dir=locations['FSYSPATH']['cesymm_minlen'] + inputf + name + "/"
            if os.path.isfile(cesymm_minlen_out_dir+inputf+name+"_mlev_stdout.out"):
                symmlist[inputf]['cesymm']={}
                symmlist[inputf]['cesymm']=cesymm_data(cesymm_minlen_out_dir,inputf+name+"_mlev")
                if 'lt' in name:
                    tag = 'lt'
                elif 'df' in name:
                    tag = 'df'
                print(tag)
                if name[:-3] == '':
                    oriented_tmp = wkdir + inputf + "_tmp.pdb"
                else:
                    oriented_tmp = chain_ordered_dir + inputf + name[:-3] + ".pdb"
                reference = parse_structure(oriented_tmp)[0]
                repeat_selection = get_repeat_resid(cesymm_minlen_out_dir,inputf+name+"_mlev",reference)
                symm_location = is_symmetry_in_membrane(limit_bottom,limit_top,repeat_selection,oriented_tmp)
                doms, crossings = symmetry_tm_domains(tm_archive, repeat_selection, inputf, '')
                if int(symmlist[inputf]['cesymm']['symmetry_levels'])>1 and (int(symmlist[inputf]['cesymm']['repeat_length'])>=40 and (doms=='in' and symm_location=='in')):
                    levels, level_type, axes_per_level, order_of_level, level_symm_type=cesymm_levels(inputf+name+"_mlev",cesymm_minlen_out_dir)
                    if 'Quaternary' in level_symm_type:
                        num_levels=level_symm_type.count('Quaternary')
                        out.write(inputf+name+"_mlev"+"\t"+tag+"\t"+symmlist[inputf]['cesymm']['seed']+"\t"+str(num_levels)+"\t"+str(symmlist[inputf]['cesymm']['repeat_length'])+"\n")
                        select_quat_dic[inputf + name + '_mlev'] = {'name': inputf + name + '_mlev',
                                                                    'mode' : tag,
                                                                    'seed': symmlist[inputf]['cesymm']['seed'],
                                                                    'required_symm_levels': num_levels,
                                                                    'cesymm_rep_length': symmlist[inputf]['cesymm']['repeat_length'],
                                                                    'natoms': tm_dic[inputf]['natoms'],
                                                                    'nres': tm_dic[inputf]['nres']}


        os.remove(wkdir+inputf+"_tmp.pdb")
    out.close()
    pkl.dump(select_quat_dic, open(locations['FSYSPATH']['cesymm_quat_minlen_log'] + "quat_symm_minlen_pdbs.pkl", "wb"))
    return select_quat_dic

if __name__ == "__main__":
    locations_path = os.path.expandvars("$ENC_DB")
    options, locations = initialize_repository(main_path=locations_path)
    pdb_suff = '_enc'
    quat_symm_minlen_dic = select_quat_symm_minlen(locations, pdb_suff)