# Name: cesymm_default.py
# Author: Antoniya A. Aleksandrova
# Language: Python 3.5+
# Description: Executes a CE-Symm run with default options and a single protein, taking into consideration guessed chain order


import time
from symmetry_exec_functions import *
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from initialize_repository import initialize_repository


start_time = time.time()


def cesymm_default(locations, options, inputf, pdb_suff, mem_req_limit = 20000, run_id="0000"):

    # Set options
    log_file = open(locations['FSYSPATH']['cesymm_log'] + inputf + '_cesymm_df.log', 'w')
    save_out = sys.stdout
    save_err = sys.stderr
    sys.stdout = sys.stderr = log_file
    print("Run id: %s" % run_id)

    if not os.path.isfile(locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb'):
        raise SystemExit("{0} does not exist.".format(locations['FSYSPATH']['whole'] + pdb_suff + '.pdb'))
    wkdir = locations['FSYSPATH']['symmtemp'] + "cesymm_default/"
    if not os.path.isdir(wkdir):
        os.mkdir(wkdir)

    print("In: " + inputf + " submitted")

    chain_ordered_dir = locations['FSYSPATH']['symmtemp'] + "chain_ordered_str/"
    out_dir = locations['FSYSPATH']['cesymm'] + inputf +"/"
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    tm_archive = locations['FSYSPATH']['symmetry'] + '.tm_archive.pkl'
    winsize = 8
    scorethreshold = 0.4
    minreplen = 15  # default options for CE-Symm
    maxsymlev = 0  # default options for CE-Symm
    minlength = 20  # 20 residues
    maxrmsd = 99

    a = pkl.load(open(tm_archive, 'rb'))

    # Find TM chains
    chain = ""
    chain = ';'.join(sorted(a[inputf]['tmchains'])) + ';'  # tm chains as a string A;B;C;...

    if chain == "":
        raise SystemExit("{0} has no TM chains.".format(locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb'))
    else:
        print(chain)

    if len(chain) > 2:
        monomer = 'N'
    else:
        monomer = 'Y'
    all_chain = int(len(chain) / 2)  # the number of tm chains

    # Run on whole protein with TM chains only
    num = strip_tm_chains(wkdir, inputf, locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb', chain)

    if num > mem_req_limit:
        cesymm_path = options['ALL']['sigcesymm']
    else:
        cesymm_path = options['ALL']['sigcesymmsmall']

    # clear_disordered_atoms(wkdir+inputf+"_tmp.pdb")
    os.rename(wkdir + inputf + "_tmp.pdb", wkdir + inputf + ".pdb")
    oriented = wkdir + inputf + ".pdb"
    if num == 0:
        raise SystemExit("{0} has 0 TM atoms.".format(locations['FSYSPATH']['whole'] + inputf + pdb_suff + '.pdb'))

    ce_symm(wkdir, cesymm_path, out_dir, inputf, oriented, winsize, scorethreshold, minreplen, maxsymlev, maxrmsd)

    # Run on each TM chain of oligomer
    if monomer == 'N':
        # Check the order of chains
        oriented_old = oriented
        if all_chain > 3:
            #        order=find_chain_order(inputf,oriented,limit_top,limit_bottom,chain)
            #        if order!=chain:
            if os.path.isfile(chain_ordered_dir + inputf + "_ordered.pdb"):
                print("In: CoM chain order")
                oriented = chain_ordered_dir + inputf + "_ordered.pdb"
                ce_symm(wkdir, cesymm_path, out_dir, inputf + '_ordered', oriented, winsize, scorethreshold, minreplen,
                        maxsymlev, maxrmsd)
                #            os.rename(wkdir+inputf+"_ordered.pdb", wkdir+"temp_structures/"+inputf+"_ordered.pdb")
                oriented = oriented_old
            if os.path.isfile(chain_ordered_dir + inputf + "_sd_order.pdb"):
                print("In: Symd chain order")
                oriented = chain_ordered_dir + inputf + "_sd_order.pdb"
                ce_symm(wkdir, cesymm_path, out_dir, inputf + "_sd_order", oriented, winsize, scorethreshold, minreplen,
                        maxsymlev, maxrmsd)
                oriented = oriented_old
            if os.path.isfile(chain_ordered_dir + inputf + "_com_sd.pdb"):
                print("In: CoM and Symd chain order")
                oriented = chain_ordered_dir + inputf + "_com_sd.pdb"
                ce_symm(wkdir, cesymm_path, out_dir, inputf + "_com_sd", oriented, winsize, scorethreshold, minreplen,
                        maxsymlev, maxrmsd)
                oriented = oriented_old

                # Get each chain
        for i in range(0, all_chain):
            ch_id = chain[i * 2]
            print("In: Chain " + ch_id)
            oriented = oriented_old
            get_chain(oriented, ch_id, inputf, wkdir)
            oriented = wkdir + inputf + "_" + ch_id + ".pdb"
            res_num = resnum(oriented)
            if res_num > minlength:
                ce_symm(wkdir, cesymm_path, out_dir, inputf + "_" + ch_id, oriented, winsize, scorethreshold, minreplen,
                        maxsymlev, maxrmsd)
            else:
                print("Out: Chain {0} has too few residues ({1}).".format(inputf + "_" + ch_id, str(res_num)))
            #        os.rename(wkdir+inputf+"_"+ch_id+".pdb",wkdir+"temp_structures/"+inputf+"_"+ch_id+".pdb")
            os.remove(wkdir + inputf + "_" + ch_id + ".pdb")

    # Clean-up
    print("{0} completed.".format(inputf))
    os.remove(wkdir + inputf + ".pdb")
    print("Time[s]:%s" % (time.time() - start_time))
    print("Date: %s" % time.asctime())
    sys.stdout = save_out
    sys.stderr = save_err
    log_file.close()
    return inputf


if __name__ == "__main__":
    if len(sys.argv) < 5:
        raise SystemExit("syntax: %s <locations_path> <pdbcode> <pdb_suffix> <run_id>" % sys.argv[0])
    else:
        locations_path = sys.argv[1].strip()
        inputf = sys.argv[2][0:4]
        pdb_suff = sys.argv[3].strip()
        run_id = sys.argv[4].strip()

    options, locations = initialize_repository(main_path=locations_path, instr_filename_path="EncoMPASS_options_relative_to__instruction_file_benchmark.txt")
    if not os.path.isdir(locations['FSYSPATH']['symmtemp']+"cesymm_default/"):
        os.mkdir(locations['FSYSPATH']['symmtemp']+"cesymm_default/")

    cesymm_default(locations, options, inputf, pdb_suff, mem_req_limit = 20000, run_id = run_id)
