# Name: select_minlen_violations_pdbs.py
# Author: Antoniya A. Aleksandrova
# Language: Python 3.5+
# Date: 4 Mar 2017
# Updated: 4 Mar 2017, 7 Apr 2022
# Description: Finds results for all proteins and chains, which have multiple levels and their repeat lengths are below 2 TM helices, and notes them down for a later run with CE-Symm and restricted levels option


import sys
import os
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from initialize_repository import initialize_repository
from symmetry_exec_functions import *
import pickle as pkl

def select_minlen_violations(locations, pdb_suff, tm_dic=False):
    tm_archive = locations['FSYSPATH']['symmetry'] + '.tm_archive.pkl'
    if tm_dic==False:
        tm_dic = pkl.load(open(tm_archive, 'rb'))

    chain_ordered_dir=locations['FSYSPATH']['symmtemp'] + "chain_ordered_str/"
    wkdir=locations['FSYSPATH']['symmtemp'] +"cesymm_minlen/"
    if os.path.isdir(wkdir)==False:
        os.mkdir(wkdir)
    out=open(locations['FSYSPATH']['cesymm_minlen_log']+"minlen_violations_pdbs.log","w")
    out.write("#Input\tCE-Symm Seed\tCE-Symm Symmetry Levels\tCE-Symm Repeat Length\n")
    select_minlen_dic = {}

    for inputf in sorted(list(tm_dic)):
        print(inputf)
        cesymm_out_dir = locations['FSYSPATH']['cesymm'] + inputf + "/"
        cesymm_lt_out_dir = locations['FSYSPATH']['cesymm_low_thr']  + inputf + "/"
        if os.path.isfile(locations['FSYSPATH']['whole']+inputf+ pdb_suff +'.pdb')==False:
            continue
        chains=sorted(tm_dic[inputf]['tmchains'])
        type=tm_dic[inputf]['class']
        if type!='alpha':
            continue
        symmlist = {}
        symmlist[inputf]={}
        # Whole pdb
        oriented_opm=locations['FSYSPATH']['whole']+inputf+ pdb_suff +'.pdb'
        limit_bottom,limit_top=membrane_limits(oriented_opm)
        ### Default case
        strip_tm_chains(wkdir,inputf,oriented_opm,chains) #should not strip them in order because they were not run that way
        oriented=wkdir+inputf+"_tmp.pdb"

        cases = [[cesymm_out_dir, '_df'], [cesymm_lt_out_dir, '_lt']]
        for case in cases:
            cesymm_dic = cesymm_data(case[0],inputf)
            symmlist[inputf]['cesymm']={}
            symmlist[inputf]['cesymm']=cesymm_dic

            if symmlist[inputf]['cesymm']['repeats_number']=='na':
                pass
            elif int(symmlist[inputf]['cesymm']['repeats_number'])>1:
                reference=parse_structure(oriented)[0]
                repeat_selection=get_repeat_resid(case[0],inputf,reference)
                symm_location=is_symmetry_in_membrane(limit_bottom,limit_top,repeat_selection,oriented)
                doms=symmetry_tm_domains(tm_archive, repeat_selection, inputf, '')
                if int(symmlist[inputf]['cesymm']['symmetry_levels'])>1 and (int(symmlist[inputf]['cesymm']['repeat_length'])<40 or (doms=='out' and symm_location=='in')):
                    print(doms, symm_location)
                    out.write(inputf + case[1] + "\t" + symmlist[inputf]['cesymm']['seed']+"\t"+str(symmlist[inputf]['cesymm']['symmetry_levels'])+"\t"+str(symmlist[inputf]['cesymm']['repeat_length'])+"\n")
                    select_minlen_dic[inputf + case[1]]= {'name': inputf + case[1],
                                                          'seed': symmlist[inputf]['cesymm']['seed'],
                                                          'cesymm_symm_levels': symmlist[inputf]['cesymm']['symmetry_levels'],
                                                          'cesymm_rep_length': symmlist[inputf]['cesymm']['repeat_length'],
                                                          'nres': tm_dic[inputf]['nres'],
                                                          'natoms': tm_dic[inputf]['natoms']}

            all_names = ['_ordered', '_com_sd', '_sd_order']
            for name in all_names:
                if os.path.isfile(case[0]+inputf+name+"_stdout.out"):
                    cesymm_dic = cesymm_data(case[0],inputf+name)
                    symmlist[inputf]['cesymm']={}
                    symmlist[inputf]['cesymm']=cesymm_dic
                    oriented_tmp=chain_ordered_dir+inputf+name+".pdb"

                    if symmlist[inputf]['cesymm']['repeats_number']=='na':
                        pass
                    elif int(symmlist[inputf]['cesymm']['repeats_number'])>1:
                        reference=parse_structure(oriented_tmp)[0]
                        repeat_selection=get_repeat_resid(case[0],inputf+name,reference)
                        symm_location=is_symmetry_in_membrane(limit_bottom,limit_top,repeat_selection,oriented_tmp)
                        doms=symmetry_tm_domains(tm_archive, repeat_selection, inputf, '')
                        if int(symmlist[inputf]['cesymm']['symmetry_levels'])>1 and (int(symmlist[inputf]['cesymm']['repeat_length'])<40 or (doms=='out' and symm_location=='in')):
                            print(doms, symm_location)
                            out.write(inputf+name+case[1]+"\t"+symmlist[inputf]['cesymm']['seed']+"\t"+str(symmlist[inputf]['cesymm']['symmetry_levels'])+"\t"+str(symmlist[inputf]['cesymm']['repeat_length'])+"\n")
                            select_minlen_dic[inputf + name + case[1]] = {'name': inputf + name + case[1],
                                                                          'seed': symmlist[inputf]['cesymm']['seed'],
                                                                          'cesymm_symm_levels': symmlist[inputf]['cesymm']['symmetry_levels'],
                                                                          'cesymm_rep_length': symmlist[inputf]['cesymm']['repeat_length'],
                                                                          'nres': tm_dic[inputf]['nres'],
                                                                          'natoms': tm_dic[inputf]['natoms']}

        os.remove(oriented)

        # Chains
        if len(chains)>1:
            pdb=inputf
            for ch in chains:
                inputf=pdb+"_"+ch
                strip_tm_chains_in_order(wkdir, inputf, oriented_opm, ch)
                oriented = wkdir + inputf + "_tmp.pdb"
                cases = [[cesymm_out_dir, '_df'], [cesymm_lt_out_dir, '_lt']]
                for case in cases:
                    if os.path.isfile(case[0]+inputf+".axes"):
                        cesymm_dic = cesymm_data(case[0],inputf)
                        symmlist[inputf]={}
                        symmlist[inputf]['cesymm']={}
                        symmlist[inputf]['cesymm']=cesymm_dic

                        if symmlist[inputf]['cesymm']['repeats_number']=='na':
                            pass
                        elif int(symmlist[inputf]['cesymm']['repeats_number'])>1:
                            reference=parse_structure(oriented)[0]
                            repeat_selection=get_repeat_resid(case[0],inputf,reference)
                            symm_location=is_symmetry_in_membrane(limit_bottom,limit_top,repeat_selection,oriented)
                            doms=symmetry_tm_domains(tm_archive, repeat_selection, pdb, ch)
                            if int(symmlist[inputf]['cesymm']['symmetry_levels'])>1 and (int(symmlist[inputf]['cesymm']['repeat_length'])<40 or (doms=='out' and symm_location=='in')):
                                print(doms, symm_location)
                                out.write(inputf+case[1]+"\t"+symmlist[inputf]['cesymm']['seed']+"\t"+str(symmlist[inputf]['cesymm']['symmetry_levels'])+"\t"+str(symmlist[inputf]['cesymm']['repeat_length'])+"\n")
                                select_minlen_dic[inputf + case[1]] = {'name': inputf + case[1],
                                                                       'seed': symmlist[inputf]['cesymm']['seed'],
                                                                       'cesymm_symm_levels': symmlist[inputf]['cesymm']['symmetry_levels'],
                                                                       'cesymm_rep_length': symmlist[inputf]['cesymm']['repeat_length'],
                                                                       'nres': tm_dic[pdb][ch]['nres'],
                                                                       'natoms': tm_dic[pdb][ch]['natoms']}

                os.remove(oriented)

    out.close()
    pkl.dump(select_minlen_dic, open(locations['FSYSPATH']['cesymm_minlen_log'] + 'minlen_violations_pdbs.pkl','wb'))
    return select_minlen_dic

if __name__ == "__main__":
    locations_path = os.path.expandvars("ENC_DB")
    options, locations = initialize_repository(main_path=locations_path)
    pdb_suff = '_enc'
    minlen_violations_dic = select_minlen_violations(locations, pdb_suff)