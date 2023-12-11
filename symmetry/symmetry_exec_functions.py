# Name: symmetry_exec_functions.py
# Author: Antoniya Aleksandrova
# Language: Python 3.5+
# All the functions used to generate the symmetry information in EnocMPASS
from __future__ import with_statement

import errno
import glob
import sys
import os
import numpy as np
import pickle as pkl
import json
import subprocess
import operator
import shutil
import re
import time

start_time = time.time()

import argparse
from Bio.PDB import PDBParser, FastMMCIFParser, Superimposer, PDBIO, Selection, Select
from Bio.PDB.Polypeptide import is_aa

from Bio import pairwise2
from Bio.Align import substitution_matrices
from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.PDB.vectors import *
from itertools import combinations


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    if all([x == 0 for x in vector]) == True:
        print("Warning: the zero vector is being processed by the unit_vector function.\n")
        return vector
    else:
        return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.dot(v1_u, v2_u))
    if np.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return np.pi
    return angle


# Atomic masses were extracted from VMD and supplemented with the wikipedia info for the ones not present in the structure

def center_of_mass(pdbfile, chain, limit_top, limit_bottom):
    """
    Calculates center of mass of a protein and/or ligand structure.
    Returns:
        center (list): List of float coordinates [x,y,z] that represent the
        center of mass (precision 3).
    Input: file.pdb chain "upper z-coord. limit of the membrane" "lower z-coord. limit of the membrane"
    """
    ATOMIC_WEIGHTS = {'BE': 9.012182, 'BA': 137.327, 'BH': 270, 'BI': 208.9804, 'BK': 247, 'BR': 79.904, 'RU': 101.07,
                      'RE': 186.207, 'RF': 267, 'LU': 174.9668, 'RA': 226, 'RB': 85.4678, 'RN': 222, 'RH': 102.9055,
                      'H': 1.0080000162124634, 'P': 30.973762, 'GE': 72.63, 'GD': 157.25, 'GA': 69.723, 'UUT': 286,
                      'OS': 190.23, 'HS': 269, 'C': 12.01099967956543, 'HO': 164.93032, 'HF': 178.49,
                      'HG': 1.0080000162124634, 'HE': 1.0080000162124634, 'PR': 140.90765, 'PT': 195.084, 'PU': 244,
                      'PB': 30.9737606048584, 'PA': 30.9737606048584, 'PD': 106.42, 'PO': 209, 'PM': 145, 'ZN': 65.38,
                      'K': 39.0983, 'O': 15.99899959564209, 'S': 32.060001373291016, 'W': 183.84, 'EU': 151.964,
                      'ZR': 91.224, 'ER': 167.259, 'MD': 258, 'MG': 24.30500030517578, 'MO': 95.96, 'MN': 54.938045,
                      'MT': 278, 'U': 238.02891, 'FR': 223, 'FE': 55.845, 'FL': 289, 'FM': 257, 'NI': 58.6934,
                      'NO': 259, 'NA': 22.98976928, 'NB': 92.90638, 'ND': 144.242, 'NE': 14.006999969482422, 'RG': 281,
                      'ES': 252, 'NP': 237, 'B': 10.81, 'CO': 58.933195, 'CN': 285, 'CM': 247, 'CL': 35.45,
                      'CA': 12.01099967956543, 'CF': 251, 'CE': 12.01099967956543, 'N': 14.006999969482422,
                      'V': 50.9415, 'CS': 132.9054519, 'CR': 51.9961, 'CU': 63.546, 'SR': 87.62, 'UUP': 288, 'UUS': 294,
                      'SI': 28.085, 'SN': 118.71, 'SM': 150.36, 'SC': 44.955912, 'SB': 121.76, 'SG': 32.060001373291016,
                      'SE': 78.96, 'YB': 173.054, 'DB': 268, 'DY': 162.5, 'DS': 281, 'LA': 138.90547, 'F': 18.9984032,
                      'LI': 6.94, 'LV': 293, 'TL': 204.38, 'TM': 168.93421, 'LR': 262, 'TH': 232.03806, 'TI': 47.867,
                      'TE': 127.6, 'TB': 158.92535, 'TC': 98, 'TA': 180.94788, 'AC': 227, 'AG': 107.8682,
                      'I': 126.90447, 'IR': 192.217, 'AM': 243, 'AL': 26.9815386, 'AS': 74.9216, 'AR': 39.948,
                      'AU': 196.966569, 'AT': 210, 'IN': 114.818, 'Y': 88.90585, 'CD': 12.01099967956543, 'XE': 131.293,
                      'OCT': 15.99899959564209}

    center = [None, None, None]

    with open(pdbfile, 'r') as pdb:

        # extract coordinates [ [x1,y1,z1], [x2,y2,z2], ... ]
        coordinates = []
        masses = []
        for line in pdb:
            if line.startswith("ATOM") and line[21:22].strip() == chain and limit_bottom < float(
                    line[46:54]) <= limit_top:
                coordinates.append([float(line[30:38]),  # x_coord
                                    float(line[38:46]),  # y_coord
                                    float(line[46:54])  # z_coord
                                    ])
                element_name = line[13:17].strip()
                if element_name not in ATOMIC_WEIGHTS:
                    element_name = line.split()[2].strip()[0]
                masses.append(ATOMIC_WEIGHTS[element_name])

        assert len(coordinates) == len(masses)

        # calculate relative weight of every atomic mass
        total_mass = sum(masses)
        weights = [float(atom_mass / total_mass) for atom_mass in masses]

        # calculate center of mass
        center = [sum([coordinates[i][j] * weights[i]
                       for i in range(len(weights))]) for j in range(3)]
        center_rounded = [round(center[i], 3) for i in range(3)]
        return center


def list_duplicates_of(seq, item):
    start_at = -1
    locs = []
    while True:
        try:
            loc = seq.index(item, start_at + 1)
        except ValueError:
            break
        else:
            locs.append(loc)
            start_at = loc
    return locs

def force_symlink(target, link_name):
    try:
        os.symlink(target, link_name)
    except OSError as e:
        if e.errno == errno.EEXIST:
            os.remove(link_name)
            os.symlink(target, link_name)
        else:
            raise e

def ce_symm(wkdir, cesymm_exe, out_dir, inputf, oriented, winsize, unrefinedtmscorethreshold, minlen, symmlevels,
            maxrmsd):
    ##### Submit file to CE-Sym version 2.0 #####
    os.chdir(wkdir)
    run_data = []
    seeds = [3, 5, 10]
    variability = '0'
    tolerance = 5
    win = "-winsize=" + str(winsize)
    unrefinedtmscore = "-unrefinedscorethreshold=" + str(unrefinedtmscorethreshold)
    minreplen = "-minlen=" + str(minlen)  # minimum length of symmetric subunit
    maxsymlev = "-symmlevels=" + str(symmlevels)  # maximum number of allowed symmetry levels
    maxrmsd = "-maxrmsd=" + str(maxrmsd)  # maximum RMSD at which to stop alignmnet optimization
    variability_file = open(out_dir + inputf + ".variability", "w")
    for run in range(1, 4):
        xmlfile = "-xml=" + wkdir + inputf + "_output_r" + str(run) + ".xml"
        rnd = "-rndseed=" + str(seeds[run - 1])
        fasta = "-fasta=" + wkdir + inputf + "_r" + str(run) + ".fasta"
        axes = "-axes=" + wkdir + inputf + "_r" + str(run) + ".axes"

        cmd = cesymm_exe.split() + ["-J"] + ["-stats=-"] + [rnd] + [oriented] + [xmlfile] + [fasta] + [axes] + [win] + [
            unrefinedtmscore] + [minreplen] + [maxsymlev] + [maxrmsd]
        cesymm = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        screen, err = cesymm.communicate()
        p_tmp = screen.decode().split('\n')
        p = p_tmp[len(p_tmp) - 2]
        variability_file.write("Seed {0}\n".format(str(rnd)))
        variability_file.write(p + "\n")
        if os.path.isfile(xmlfile[5:]) and os.path.getsize(xmlfile[5:]) > 5:
            out = open(wkdir + "/" + inputf + "_stdout_r" + str(run) + ".out", "w")
            out.write(p)
            out.close()
            fields = p.split()
            ind = 0
            try:
                int(fields[4])
                ce_order = fields[
                    2]  # 18
            except ValueError:
                ind = -1
                ce_order = "na"
            num_repeats = fields[1]  # 17
            tmScore = float(fields[ind + 10])  # 24; refined
            ce_rmsd = float(fields[ind + 11])  # 25; refined
            repeat_length = int(fields[ind + 12])
            sym_level = int(fields[ind + 4])  # 20
            core_length = int(fields[ind + 13])  # 29
            length = int(fields[ind + 14])  # 30
            coverage = fields[ind + 15]  # 31
            cce_type = "na"
            angle = "na"
            internal = "N"
            quaternary = "N"
            axis_angle = "na"
            detected_symmetry = 'na'
            if sym_level > 0:
                cce_type = ""
                angle = ""
                f = open(axes[6:], 'r')
                for line in f:
                    if inputf in line:
                        flds = line.split()
                        angle = angle + flds[4] + ";"
                        pt1 = [float(x) for x in flds[6].split(",")]
                        pt1 = np.array(pt1)
                        pt2 = [float(x) for x in flds[7].split(",")]
                        pt2 = np.array(pt2)
                        u_tmp = pt1 - pt2
                        u = u_tmp / np.linalg.norm(u_tmp)
                        axis_angle = round(angle_between(u, (0, 0, 1)) * 180 / np.pi, 2)
                        if np.abs(np.dot(u, (0, 0, 1))) > 0.5:
                            cce_type = cce_type + "orthogonal"
                        else:
                            cce_type = cce_type + "parallel"
                        rep = flds[8].split(")")
                        chain = [rep[0][pos + 1] for pos, char in enumerate(rep[0]) if
                                 char == '.']  # finds all chain names in first symmetry bracket (they are preceded by .)
                        if chain.count(chain[0]) != len(
                                chain):  # if all repeats are from the same chain, no quaternary structure symmetry was detected
                            detected_symmetry = "Quat"
                            quaternary = 'Y'
                        else:
                            detected_symmetry = "Int"
                            internal = 'Y'
                        cce_type = cce_type + "-" + detected_symmetry + ";"

                f.close()
            run_data.append(
                [ce_order, sym_level, ce_rmsd, tmScore, cce_type, num_repeats, internal, quaternary, angle, axis_angle,
                 detected_symmetry, core_length, length, coverage, repeat_length])


        else:
            ce_order = 'na'
            sym_level = 'na'
            ce_rmsd = 'na'
            tmScore = 'na'
            cce_type = 'na'
            num_repeats = 'na'
            internal = "na"
            quaternary = "na"
            angle = "na"
            axis_angle = "na"
            detected_symmetry = 'na'
            coverage = 'na'
            length = 'na'
            core_length = 0
            variability = 'na'
            repeat_length = 0
            if not os.path.isfile(xmlfile[5:]):
                print("Error running CE-Symm - xmlfile does not exist")
                break
            else:
                out = open(wkdir + "/" + inputf + "_stdout_r" + str(run) + ".out", "w")
                out.write(p)
                out.close()
                run_data.append(
                    [ce_order, sym_level, ce_rmsd, tmScore, cce_type, num_repeats, internal, quaternary, angle,
                     axis_angle, detected_symmetry, core_length, length, coverage, repeat_length])

    if len(run_data) > 2:
        ind = 0
        print("Core length in each run: ", run_data[0][11], run_data[1][11], run_data[2][11])
        # Compare number of repeats and core length
        # The goal here is not to discard (filter) results but to make sure that if there is an option with at least two TM crossings in a repeat, it is selected
        if run_data[0][5] == run_data[1][5] and run_data[0][5] == run_data[2][5] and np.abs(
                int(run_data[0][11]) - int(run_data[1][11])) < tolerance and np.abs(
                int(run_data[0][11]) - int(run_data[2][11])) < tolerance and np.abs(
                int(run_data[2][11]) - int(run_data[1][11])) < tolerance:
            variability = '0'
            cores = [run_data[0][11], run_data[1][11], run_data[2][11]]
            ind = np.argmax(cores)
        else:
            variability = '1'
            rep_len = [run_data[0][14], run_data[1][14], run_data[2][14]]  # length of repeats
            repeats = [run_data[0][5], run_data[1][5], run_data[2][5]]  # number of repeats
            cores = [run_data[0][11], run_data[1][11], run_data[2][11]]  # core length

            if (all(i >= 40 for i in rep_len) or all(
                    i < 40 for i in rep_len)) == False:  # Check for repeats that have at least two TM crossings
                for i, r in enumerate(rep_len):
                    if r < 40:
                        repeats[i] = 0
                        cores[i] = 0

            ind = np.argmax(repeats)
            if repeats.count(repeats[0]) == len(repeats):
                ind = np.argmax(cores)
        """
        ce_order=run_data[ind][0]
        sym_level=run_data[ind][1]
        ce_rmsd=run_data[ind][2]
        tmScore=run_data[ind][3]
        cce_type=run_data[ind][4]
        num_repeats=run_data[ind][5]
        internal=run_data[ind][6]
        quaternary=run_data[ind][7]
        angle=run_data[ind][8]
        axis_angle=run_data[ind][9]
        detected_symmetry=run_data[ind][10]
        core_length=run_data[ind][11]
        length=run_data[ind][12]
        coverage=run_data[ind][13]
        """

        src = wkdir + inputf + "_output_r" + str(ind + 1) + ".xml"
        dst = out_dir + inputf + "_output.xml"
        shutil.copyfile(src, dst)
        src = wkdir + inputf + "_stdout_r" + str(ind + 1) + ".out"
        dst = out_dir + inputf + "_stdout.out"
        shutil.copyfile(src, dst)
        src = wkdir + inputf + "_r" + str(ind + 1) + ".fasta"
        dst = out_dir + inputf + ".fasta"
        shutil.copyfile(src, dst)
        src = wkdir + inputf + "_r" + str(ind + 1) + ".axes"
        dst = out_dir + inputf + ".axes"
        shutil.copyfile(src, dst)
        [os.remove(wkdir + inputf + "_output_r" + str(i) + ".xml") for i in range(1, 4)]
        [os.remove(wkdir + inputf + "_stdout_r" + str(i) + ".out") for i in range(1, 4)]
        [os.remove(wkdir + inputf + "_r" + str(i) + ".fasta") for i in range(1, 4)]
        [os.remove(wkdir + inputf + "_r" + str(i) + ".axes") for i in range(1, 4)]
        seed_file = open(out_dir + inputf + ".seed", "w")
        seed_file.write(str(seeds[ind]))
        seed_file.close()
        variability_file.write("Variability: " + variability + "\n")
        variability_file.close()
        print("Out: CE-Symm procedure completed successfully.")
    else:
        print("Incomplete run: only {0} runs were completed successfully.".format(len(run_data)))

    os.chdir(wkdir)
    return


def ce_symm_minlen(cesymm_exe, out_dir, inputf, oriented, winsize, unrefinedtmscorethreshold, minlen, symmlevels, order,
                   seed, suffix):
    ##### Submit file to CE-Sym version 2.0 #####
    win = "-winsize=" + str(winsize)
    unrefinedtmscore = "-unrefinedscorethreshold=" + str(unrefinedtmscorethreshold)
    minreplen = "-minlen=" + str(minlen)  # minimum length of symmetric subunit
    maxsymlev = "-symmlevels=" + str(symmlevels)  # maximum number of allowed symmetry levels
    order = "-order=" + str(order)  # enforced order
    xmlfile = "-xml=" + out_dir + inputf + suffix + "_output.xml"
    rnd = "-rndseed=" + str(seed)
    fasta = "-fasta=" + out_dir + inputf + suffix + ".fasta"
    axes = "-axes=" + out_dir + inputf + suffix + ".axes"
    seed_file = open(out_dir + inputf + suffix + ".seed", "w")
    seed_file.write(str(seed))
    seed_file.close()

    cmd = cesymm_exe.split() + ["-J"] + ["-stats=-"] + [rnd] + [oriented] + [xmlfile] + [fasta] + [axes] + [win] + [
        unrefinedtmscore] + [minreplen] + [maxsymlev] + [order]
    cesymm = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    screen, err = cesymm.communicate()
    p_tmp = screen.decode().split('\n')
    p = p_tmp[len(p_tmp) - 2]
    if os.path.isfile(xmlfile[5:]):
        out = open(out_dir + "/" + inputf + suffix + "_stdout.out", "w")
        out.write(p)
        out.close()

    else:
        print("Error running CE-Symm - xmlfile does not exist")

    return


def symd(wkdir, symd_exe, out_dir, inputf, oriented, atoms):
    ##### Submit file to SymD #####
    os.chdir(wkdir)
    cmd = symd_exe.split() + ["-mxna"] + [str(atoms)] + [oriented]
    symd = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p, err = symd.communicate()
    if err.decode() == '' and os.path.isfile(wkdir + inputf + "-info.txt"):
        os.rename(wkdir + inputf + "-best.fasta", out_dir + inputf + "-best.fasta")
        os.rename(wkdir + inputf + "-info.txt", out_dir + inputf + "-info.txt")
        os.rename(wkdir + inputf + "-trfm.pdb", out_dir + inputf + "-trfm.pdb")
        print("Out: SymD procedure completed successfully.")
    else:
        print("Error running SymD")
        print(err.decode())

    return


def ananas(wkdir, ananas_exe, out_dir, inputf, oriented):
    """
    Submit file to AnAnaS and store the output
    """
    output = open(wkdir + inputf + "_ananas.out", "w")
    json_out = out_dir + inputf + "_ananas.json"
    fnull = open(os.devnull, 'w')
    cmd = ananas_exe.split() + [oriented] + ["-a"] + ["-p"] + ["--json"] + [json_out]
    p = subprocess.Popen(cmd, stdout=output, stderr=fnull)
    p.wait()
    fnull.close()
    output.close()
    if os.path.isfile(json_out):
        os.rename(wkdir + inputf + "_ananas.out", out_dir + inputf + "_ananas.out")
        print("Out: AnAnaS procedure completed successfully.")
    else:
        print("Error running AnAnaS")
    return inputf


def quatsymm(quatsymm_exe, out_dir, inputf, oriented, suffix, threads='2', minseqlen='', minseqid='', minseqcov='',
             minstructcov='', maxclustrmd='', clustmethod='', maxsymmrmsd=''):
    ##### Submit file to QuatSymm version 2.2.3 #####
    # Note that QuatSymm doesn't work with HETATM DUM entries
    # ./runQuatSymm.sh -J 2nwx --stats=2nwx.out --fasta=2nwx.fasta --threads=2
    opts = {'--threads': threads, '--minSeqLen': minseqlen, '--minSeqId': minseqid, '--minSequenceCoverage': minseqcov,
            '--minStructureCoverage': minstructcov, '--maxClustRmsd': maxclustrmd, '--clustMethod': clustmethod,
            '--maxSymmRmsd': maxsymmrmsd}
    cmd = quatsymm_exe.split() + ["-J"] + ["--stats=" + out_dir + inputf + suffix + ".out"] \
          + [oriented] + ["--fasta=" + out_dir + inputf + suffix + ".fasta"]
    for entry in opts:
        if opts[entry] != '':
            cmd = cmd + [entry + "=" + str(opts[entry])]

    quatsymm = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    screen, err = quatsymm.communicate()
    tmp_err = err.decode().strip().split(
        '\n')  # handling the failed attempts to download amino acid cif files when there is no internet

    #    if err.decode()=='' and os.path.isfile(out_dir + inputf + suffix + ".out"):
    if os.path.isfile(out_dir + inputf + suffix + ".out") and (
            err.decode() == '' or all(['Could not download' in ln for ln in tmp_err])):
        print("Out: QuatSymm procedure completed successfully.")

    else:
        print("Error running QuatSymm - out does not exist")
        print(err.decode())
    return


def get_chain(oriented, ch_id, inputf, out_dir):
    """
    We are assuming that the pdb we are passing to this functions has been
    cleaned so that all important HETATM have been renamed to ATOM
    (see strip_pdb function for details)
    """
    f = open(oriented, 'r')
    o = open(out_dir + inputf + "_" + ch_id + ".pdb", "w")
    f.seek(0)
    for line in f:
        if line.startswith("ATOM") and line[21:22] == ch_id:
            o.write(line)
    o.write("END\n")
    f.close()
    o.close()


def chain_id(oriented):
    # obtain the ids for all chains in case the pdb is not listed in the OPM
    f = open(oriented, 'r')
    chain = ''
    ch_old = ''
    for line in f:
        if line.startswith("ATOM"):
            ch_id = line[21:22]
            if ch_old != ch_id:
                chain = chain + ch_id + ";"
                ch_old = ch_id
    return chain


def halve_chain(oriented, inputf):
    f = open(oriented, 'r')
    o = open(inputf + "_h.pdb", 'w')
    num = 0
    for line in f:
        if line.startswith("ATOM"):
            num = num + 1
    f.seek(0)
    count = 0
    for line in f:
        if line.startswith("ATOM") and count <= num / 2:
            count = count + 1
            o.write(line)
    o.write("END\n")
    f.close()
    o.close()


def num_atoms_pdb(oriented):
    """
    Finds the number of atoms in a pdb file. 
    Important as a tool to check if number of atoms in pdb exceeds 10,000 
    and therefore will crash a SymD run 
    """
    f = open(oriented, 'r')
    num = 0
    for line in f:
        if line.startswith("ATOM"):
            num = num + 1
    f.close()
    return num


def order_pdbs_by_size(pdb_list):
    return


def opm_sql(opm_dir, inputf):
    ###### Working with OPM SQL file
    sql = open(opm_dir + "OPM06-05-16.sql", "r")
    flag = 0
    chain = ''
    name = ''
    family_name = ''
    opm_id = inputf
    family = 'not applicable'
    #  Check whether the pdb id is related to another structure or is representative
    for line in sql:
        if "-- Dumping data for table `protein`" in line:
            flag = 2
        if ("-- Dumping data for table `relatedproteins" in line) and flag == 2:
            flag = 6
        if flag == 2 and inputf in line:
            flds = line.split(',')
            family = flds[1].strip()
            flds = line.split("'")
            name = '"' + flds[3] + '"'
            flag = 0
        if flag == 6 and inputf in line:
            flds = line.split("'")
            if inputf in flds[3]:
                opm_id = flds[1]
                break
    sql.seek(0)
    flag = 0
    #  Get subunits (chains) for the protein
    for line in sql:
        if "-- Dumping data for table `protein`" in line:
            flag = 2
        if flag == 2 and opm_id in line:
            flds = line.split(',')
            family = flds[1].strip()
            flds = line.split("'")
            name = '"' + flds[3] + '"'
            flag = 0
        if "-- Dumping data for table `subunits`" in line:
            flag = 1
        if flag == 1 and opm_id in line:
            flds = line.split()[3]
            chain = chain + flds[1:][:-2] + ";"
        if flag == 1 and '-- -' in line:
            break
    sql.seek(0)
    #  Get name and family of the protein
    flag = 0
    superfml = 'not_applicable'
    for line in sql:
        if "-- Dumping data for table `family`" in line:
            flag = 1
        if flag == 1 and line.startswith("(" + family + ','):
            index = line.split()[1]
            superfml = index[:-5] + "'"
            flag = 0
        if "-- Dumping data for table `superfamily`" in line:
            flag = 2
        if flag == 2 and superfml in line:
            fld = line.split("'")[3]
            family_name = '"' + fld + '"'
            flag = 0
    sql.close()
    return chain, name, family_name


def symd_sequence(dir, inputf, oriented):
    RESIDUES = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
                'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y', 'MSE': 'M', 'UNK': '*'}
    alignment = open(dir + inputf + "-best.fasta", "r")
    flag = 0
    sequence = ""
    sequence2 = ""
    for line in alignment:
        if line.startswith(">") and flag == 1:
            flag = 2
        if line.startswith(">") and flag == 0:
            flag = 1
        if not line.startswith(">") and flag == 1:
            sequence = sequence + line.strip()
        if not line.startswith(">") and flag == 2:
            sequence2 = sequence2 + line.strip()

    flag = 0
    count = 0
    first_cap = 0
    last_cap = 0
    for ind, c in enumerate(sequence):
        if c == '-':
            count = count + 1
        if (c.isupper() or c == '*') and flag == 0:
            first_cap = ind - count
            first_cap_symd = ind
            flag = 1
        if (c.isupper() or c == '*') and flag == 1:
            last_cap = ind - count
            last_cap_symd = ind
            last_let = c
    protein = open(oriented, "r")
    ind = 0
    print("oriented: ", oriented, first_cap, last_cap)
    for line in protein:
        if line.startswith("ATOM") and line[13:16].strip() == "CA":
            if ind == first_cap:
                resid_first_cap = line[22:27].strip()
            if ind == last_cap:
                resid_last_cap = line[22:27].strip()
            ind = ind + 1

    protein.seek(0)
    pdbseq = ""
    flag = 0
    for line in protein:
        if line.startswith("ATOM") and line[13:16].strip() == "CA":
            resid = line[22:27]
            if flag == 1 and resid != old_resid:
                pdbseq = pdbseq + RESIDUES[line[17:20]]
                old_resid = resid
            elif flag == 0:
                flag = 1
                old_resid = resid
                pdbseq = pdbseq + RESIDUES[line[17:20]]
    tag = ""
    i = 0
    for c in sequence2[first_cap_symd:]:
        if c != "-":
            tag = tag + c
            i = i + 1
        if i == 10:
            break
    tag = tag.upper()
    first_cap = pdbseq.find(tag)
    if first_cap == -1:  # this takes into account a case when the first_cap will be at the very end of the PDB sequence
        circular = pdbseq + pdbseq
        first_cap = circular.find(tag)
    tag = ""
    i = 0
    temp = sequence2[0:last_cap_symd + 1]
    for c in temp[::-1]:
        if c != "-":
            tag = c + tag
            i = i + 1
        if i == 10:
            break
    tag = tag.upper()
    add = len(
        tag) - 1  # this takes into account a case when the matched sequence is very short and there's nothing after it
    last_cap = pdbseq.find(tag) + add

    protein.seek(0)
    ind = 0
    flag = 0
    for line in protein:
        if line.startswith("ATOM") and line[13:16].strip() == "CA":
            if flag == 0:
                pdb_first = line[22:27].strip()
                flag = 1
            if ind == first_cap:
                match_first_cap = line[22:27].strip()
            if ind == last_cap:
                match_last_cap = line[22:27].strip()
            pdb_last = line[22:27].strip()
            ind = ind + 1

    return resid_first_cap, resid_last_cap, match_first_cap, match_last_cap, pdb_first, pdb_last


def residues_extract(wkdir, inputf, oriented, first_resid, last_resid):
    f = open(oriented, 'r')
    o = open(wkdir + inputf + "_" + first_resid + "-" + last_resid + ".pdb", 'w')
    resids = range(int(first_resid), int(last_resid) + 1)
    for line in f:
        # if line.startswith("ATOM") and int(first_resid)<=int(re.search(r'\d+',line[22:27]).group())<=int(last_resid):
        if line.startswith("ATOM") and int(re.search(r'\d+', line[22:27]).group()) in resids:
            o.write(line)
    o.write("END\n")
    f.close()
    o.close()


def strip_pdb(dir, inputf, oriented_enc):
    f = open(oriented_enc, 'r')
    o = open(dir + "/" + inputf + "_tmp.pdb", "w")
    LINELEM = "{:76s}{:>2s}\n"
    for line in f:
        if (line.startswith("ATOM") or line.startswith("TER ")) and from3to1(line[17:20].strip()) != 0:
            if line[76:78].strip() != '':
                o.write(line)
            else:
                if line.startswith("TER "):
                    o.write(line)
                else:
                    atom = line[12:16].strip()
                    elem = atom[0]
                    o.write(LINELEM.format(line[0:76], elem))
        if line.startswith("TER\n"):
            o.write(line)
        if line.startswith("HETATM") and line[17:20] == 'MSE':
            if line[76:78].strip() != '':
                o.write("ATOM  " + line[6:])
            else:
                atom = line[12:16].strip()
                elem = atom[0]
                o.write("ATOM  " + line[6:76] + " " + elem + "\n")

    o.write("END\n")
    f.close()
    o.close()


def tm_chains(dir, inputf, oriented_tmp, tm_chains):
    f = open(oriented_tmp, 'r')
    o = open(dir + "/" + inputf + ".pdb", "w")
    chains = tm_chains.split(';')
    for line in f:
        if (line.startswith("ATOM") or line.startswith("TER ")) and line[21:22] in chains:
            o.write(line)
        if line.startswith("END"):
            o.write(line)
    f.close()
    o.close()
    os.remove(oriented_tmp)


def membrane_limits(oriented):
    f = open(oriented, 'r')
    limit = []
    flag = 0
    for line in f:
        if (line.startswith("HETATM") or line.startswith("ATOM")) and line[17:20] == 'DUM' and line[
                                                                                               12:16].strip() == 'N':
            limit.append(float(line[46:54]))
            flag = flag + 1
        if (line.startswith("HETATM") or line.startswith("ATOM")) and line[17:20] == 'DUM' and line[
                                                                                               12:16].strip() == 'O':
            limit.append(float(line[46:54]))
            flag = flag + 1
        if flag == 2:
            break
    f.close()
    limit = sorted(limit)
    return limit[0], limit[1]


def order_pdb_chains(inputf, oriented, ordered_chains):
    f = open(oriented, 'r')
    outputf = open(inputf + "_ordered.pdb", "w")
    for i in range(0, len(ordered_chains)):
        for line in f:
            if (line.startswith("ATOM") or line.startswith("TER ")) and line[21:22].strip() == ordered_chains[i]:
                outputf.write(line)
        f.seek(0)
    f.close()
    outputf.write("END\n")
    outputf.close()


def find_chain_order(inputf, oriented, limit_top, limit_bottom, chain):
    com_chains = []
    all_chain = len(chain) // 2
    for i in range(0, all_chain):
        ch_id = chain[i * 2]
        com = center_of_mass(oriented, ch_id, limit_top, limit_bottom)
        com_chains.append(com)

    dist_between_com = []
    for i in range(0, len(com_chains)):
        temp_dist = []
        for j in range(0, len(com_chains)):
            temp_dist.append(np.linalg.norm(np.array(com_chains[i]) - np.array(com_chains[j])))
        dist_between_com.append(temp_dist)

    max_dist = np.max(dist_between_com)
    first_index, value = max(enumerate(dist_between_com), key=operator.itemgetter(1))
    first_chain = chain[first_index * 2]
    sorted_dist = []
    for i in range(0, len(com_chains)):
        sorted_dist.append(sorted(dist_between_com[i]))

    neighbors = []
    neighbor_ind = []
    for i in range(0, len(com_chains)):
        tmp = []
        tmp_indeces = []
        for j in range(1, len(com_chains)):
            lst_of_ind = list_duplicates_of(dist_between_com[i], sorted_dist[i][j])
            if len(lst_of_ind) == 1:
                ind = lst_of_ind[0]
                tmp.append(chain[ind * 2])
                tmp_indeces.append(ind)
            else:
                for ind in lst_of_ind:
                    if chain[ind * 2] not in tmp:
                        tmp.append(chain[ind * 2])
                        tmp_indeces.append(ind)
                        break

        neighbors.append(tmp)
        neighbor_ind.append(tmp_indeces)

    order = [first_chain]
    ind = first_index
    for i in range(1, len(com_chains)):
        for j in range(0, len(com_chains) - 1):
            if neighbors[ind][j] not in order:
                next = neighbors[ind][j]
                break
        order.append(next)
        ind = chain.find(next) // 2

    order_pdb_chains(inputf, oriented, order)
    new_chain = ''
    for i in range(0, len(order)):
        new_chain = new_chain + order[i] + ";"
    if len(new_chain.split(';')) == len(set(new_chain.split(';'))):
        return new_chain
    else:
        print(
            "Warning: repeated chain in the predicted order (%s). Please check that the membrane limits are selected correctly in case of multiple membranes (%f, %f)." % (
            new_chain, limit_bottom, limit_top))
        return ""


def resnum(oriented):
    f = open(oriented, "r")
    res_num = 0
    old_resid = -1
    for line in f:
        if line.startswith("ATOM") and line[13:16].strip() == "CA" and int(
                re.search(r'\d+', line[22:27]).group()) != old_resid:
            res_num += 1
            old_resid = int(re.search(r'\d+', line[22:27]).group())
    return res_num


def from3to1_general(resname):
    f3t1 = {'ALA': 'A',
            'ARG': 'R',
            'ASN': 'N',
            'ASP': 'D',
            'CYS': 'C',
            'GLN': 'Q',
            'GLU': 'E',
            'GLY': 'G',
            'HIS': 'H',
            'ILE': 'I',
            'LEU': 'L',
            'LYS': 'K',
            'MET': 'M',
            'PHE': 'F',
            'PRO': 'P',
            'SER': 'S',
            'THR': 'T',
            'TRP': 'W',
            'TYR': 'Y',
            'VAL': 'V',
            'MSE': 'M'}

    if resname in list(f3t1.keys()):
        return f3t1[resname]
    else:
        return '0'


def find_chains(oriented):
    f = open(oriented, "r")
    chains = ""
    for line in f:
        if line.startswith("ATOM") and line[21:22].strip() not in chains:
            chains = chains + line[21:22].strip() + ";"
    f.close()
    return chains


######### Alignment functions 
def get_pdb_sequence(structure):
    """
    Retrieves the AA sequence from a PDB structure.
    """

    _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, 'X'))
    seq = [_aainfo(r) for r in structure.get_residues() if is_aa(r)]
    return seq


def get_pdb_sequence_with_chains(structure):
    """
    Retrieves the AA sequence from a PDB structure. It's a list that looks like [(5, 'R', 'A'), (6, 'E', 'A'), (7, 'H', 'A'), (8, 'W', 'A'),...]
    """

    _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, 'X'), r.get_parent().get_id(), r.id[0], r.id[2])
    seq = [_aainfo(r) for r in structure.get_residues() if (is_aa(r) and r.has_id('CA'))]
    return seq


def calculate_identity(sequenceA, sequenceB):
    """
    Returns the percentage of identical characters between two sequences.
    Assumes the sequences are aligned.
    """

    sa, sb, sl = sequenceA, sequenceB, len(sequenceA)
    matches = [sa[i] == sb[i] for i in range(sl)]
    seq_id = (100 * sum(matches)) / sl

    gapless_sl = sum([1 for i in range(sl) if (sa[i] != '-' and sb[i] != '-')])
    gap_id = (100 * sum(matches)) / gapless_sl
    return (seq_id, gap_id)


def align_sequences(structA, structB, **kwargs):
    """
    Performs a global pairwise alignment between two sequences
    using the BLOSUM62 matrix and the Needleman-Wunsch algorithm
    as implemented in Biopython. Returns the alignment, the sequence
    identity and the residue mapping between both original sequences.
    """

    matrix = kwargs.get('matrix', substitution_matrices.load("BLOSUM62"))
    gap_open = kwargs.get('gap_open', -10.0)
    gap_extend = kwargs.get('gap_extend', -0.5)
    trim_ends = kwargs.get('trim_ends', True)

    resseq_A = get_pdb_sequence(structA)
    resseq_B = get_pdb_sequence(structB)

    sequence_A = ''.join([i[1] for i in resseq_A])
    sequence_B = ''.join([i[1] for i in resseq_B])
    alns = pairwise2.align.globalds(sequence_A, sequence_B,
                                    matrix, gap_open, gap_extend,
                                    penalize_end_gaps=(False, False))

    best_aln = alns[0]
    aligned_A, aligned_B, score, begin, end = best_aln

    # Equivalent residue numbering
    # Relative to reference
    mapping = {}
    aa_i_A, aa_i_B = 0, 0
    for aln_i, (aa_aln_A, aa_aln_B) in enumerate(zip(aligned_A, aligned_B)):
        if aa_aln_A == '-':
            if aa_aln_B != '-':
                aa_i_B += 1
        elif aa_aln_B == '-':
            if aa_aln_A != '-':
                aa_i_A += 1
        else:
            assert resseq_A[aa_i_A][1] == aa_aln_A
            assert resseq_B[aa_i_B][1] == aa_aln_B
            mapping[resseq_A[aa_i_A][0]] = resseq_B[aa_i_B][0]
            aa_i_A += 1
            aa_i_B += 1

    # Gapless alignment
    # def _trimmer(sequence):
    #     """Returns indices of first and last ungapped position"""

    #     leading = [i for (i, aa) in enumerate(sequence) if aa != '-'][0]
    #     trailing = [i for (i, aa) in enumerate(sequence[::-1]) if aa != '-'][0]

    #     trailing = len(sequence) - trailing
    #     return (leading, trailing)

    # lead_A, trail_A = _trimmer(aligned_A)
    # lead_B, trail_B = _trimmer(aligned_B)

    # lead = max(lead_A, lead_B)
    # trail = min(trail_A, trail_B)
    # trim_aln_A = aligned_A[lead:trail]
    # trim_aln_B = aligned_B[lead:trail]
    # mismatch = ''.join(['+' if a!=b else ' ' for (a,b) in zip(trim_aln_A, trim_aln_B)])

    # Calculate (gapless) sequence identity
    seq_id, g_seq_id = calculate_identity(aligned_A, aligned_B)
    return ((aligned_A, aligned_B), seq_id, g_seq_id, mapping)
    # return ((trim_aln_A, trim_aln_B, mismatch), seq_id, g_seq_id, mapping)


def align_sequences_with_chains(wkdir, structA, structB, **kwargs):
    """
    Performs a global pairwise alignment between two sequences
    using the BLOSUM62 matrix and the Needleman-Wunsch algorithm
    as implemented in Biopython. Returns the alignment, the sequence
    identity and the residue mapping between both original sequences.
    Gaps are not included in the residue map.
    """

    matrix = kwargs.get('matrix', matlist.blosum62)
    gap_open = kwargs.get('gap_open', -10.0)
    gap_extend = kwargs.get('gap_extend', -0.5)
    trim_ends = kwargs.get('trim_ends', True)
    muscle_path = kwargs.get('muscle_path', "muscle")

    resseq_A = get_pdb_sequence_with_chains(structA)
    resseq_B = get_pdb_sequence_with_chains(structB)
    #    print(resseq_B)

    sequence_A = ''.join([i[1] for i in resseq_A])
    sequence_B = ''.join([i[1] for i in resseq_B])

    try:
        if os.path.isfile(wkdir + "out.fasta"):
            os.remove(wkdir + "out.fasta")
        if os.path.isfile(wkdir + "tmp.fasta"):
            os.remove(wkdir + "tmp.fasta")
        f = open(wkdir + "tmp.fasta", "w")
        f.write(">seqA\n")
        f.write(sequence_A + "\n")
        f.write(">seqB\n")
        f.write(sequence_B + "\n")
        f.close()
        fnull = open(os.devnull, 'w')
        p = subprocess.Popen(muscle_path.split() + ['-in', wkdir + "tmp.fasta", '-out', wkdir + "out.fasta"],
                             stdout=fnull, stderr=fnull)
        p.wait()
        fnull.close()
        tmp = SeqIO.to_dict(SeqIO.parse(wkdir + "out.fasta", "fasta"))
        aligned_A = tmp['seqA'].seq
        aligned_B = tmp['seqB'].seq
        os.remove(wkdir + "out.fasta")
        os.remove(wkdir + "tmp.fasta")
    except:
        print(sys.exc_info()[0])
        print("Muscle procedure failed. Reverting to BioPython...")
        alns = pairwise2.align.globalds(sequence_A, sequence_B,
                                        matrix, gap_open, gap_extend,
                                        penalize_end_gaps=(False, False))

        best_aln = alns[0]
        aligned_A, aligned_B, score, begin, end = best_aln

    #    print(aligned_A)
    #    print(aligned_B)

    # Equivalent residue numbering
    # Relative to reference
    mapping = {}
    aa_i_A, aa_i_B = 0, 0
    for aln_i, (aa_aln_A, aa_aln_B) in enumerate(zip(aligned_A, aligned_B)):
        if aa_aln_A == '-':
            if aa_aln_B != '-':
                aa_i_B += 1
        elif aa_aln_B == '-':
            if aa_aln_A != '-':
                aa_i_A += 1
        else:
            assert resseq_A[aa_i_A][1] == aa_aln_A
            assert resseq_B[aa_i_B][1] == aa_aln_B
            mapping[(resseq_A[aa_i_A][0], resseq_A[aa_i_A][2], resseq_A[aa_i_A][3], resseq_A[aa_i_A][4])] = (
            resseq_B[aa_i_B][0], resseq_B[aa_i_B][2], resseq_B[aa_i_B][3], resseq_B[aa_i_B][4])
            aa_i_A += 1
            aa_i_B += 1

    # Calculate (gapless) sequence identity
    seq_id, g_seq_id = calculate_identity(aligned_A, aligned_B)
    return ((aligned_A, aligned_B), seq_id, g_seq_id, mapping)


def frtmaln(wkdir, template, target, frtmalign_path):
    """
    Performs a Fr-TMAlign alignmnet
    """
    # dir=os.getcwd()

    # os.chdir(wkdir)
    if not os.path.isfile(template):
        raise SystemExit(f"Template {template} doesn't exist.")
    if not os.path.isfile(target):
        raise SystemExit(f"Target {target} doesn't exist.")
    template_pdb = template.split('/')[-1:][0]
    target_pdb = target.split('/')[-1:][0]
    if template != wkdir + template_pdb:
        shutil.copy2(template, wkdir + template_pdb)
    if target != wkdir + target_pdb:
        shutil.copy2(target, wkdir + target_pdb)
    out = wkdir + "frtm.out"
    output = open(out, "w")
    if ".sif" in frtmalign_path or "singularity" in frtmalign_path:
        if shutil.which("singularity") is None:
            raise SystemExit("Load/install singularity and export SINGULARITY_BINDPATH before running.")

    fnull = open(os.devnull, 'w')
    p = subprocess.Popen(frtmalign_path.split() + [template_pdb, target_pdb], stdout=output, stderr=fnull, cwd=wkdir)
    p.wait()
    fnull.close()
    output.close()
    f = open(out, "r")
    o = open(wkdir + "frtm.fasta", "w")
    flag = 0
    for line in f:
        if flag == 3:
            o.write(line)
        if flag == 1 and ':' not in line:
            o.write(line)
        elif ':' in line and flag == 1:
            flag = 2
        if flag == 2 and ':' not in line:
            o.write(">" + target_pdb[:-4] + "\n")
            flag = 3
            o.write(line)
        if "denotes the residue pairs of distance" in line:
            flag = 1
            o.write(">" + template_pdb[:-4] + "\n")
    o.close()
    f.close()
    os.remove(out)
    if template != wkdir + template_pdb:
        os.remove(wkdir + template_pdb)
    if target != wkdir + target_pdb:
        os.remove(wkdir + target_pdb)
    return

def extract_frtmalign_fasta(fname, outname, pdb1, pdb2):
    """
    Compatible with the 2018 version of EncoMPASS
    :param fname: file that contains a list of alignments and statistics
    :param outname: name of the extracted fasta file we want to create
    :param pdb1: name of one of the structures in the aligned pair
    :param pdb2: name of the other structure in the aligned pair
    :return: 0 if the structures were never found, 1 if the alignment was incomplete, 2 if all the information was found
    """
    f = open(fname, 'r')
    flag = 0
    out = open(outname, "w")
    for line in f:
        if flag == 1:
            out.write(line)
        if line.startswith('BEGIN') and pdb1 in line and pdb2 in line:
            flag = 1
        if flag == 1 and line.startswith('RMSD'):
            flag = 2
        if flag == 2:
            break
    f.close()
    out.close()
    return flag

def extract_frtm_alignment(fname, outname, pdb1, pdb2):
    """
    Compatible with 2022+ version of EncoMPASS
    :param fname:
    :param outname:
    :param pdb1:
    :param pdb2:
    :return:
    """
    out = open(outname, "w")
    flag = 0
    with open(fname, "r") as f:
        for line in f:
            if flag == 1:
                out.write(line)
            if line.startswith("INIT") and flag == 1:
                flag = 2
            if line.startswith("INIT") and pdb1 in line and pdb2 in line:
                flag = 1
    out.close()
    return flag

def frtmalign_sequences_with_chains(wkdir, structA, structB, file, oriented_template, oriented_target, template, target,
                                    **kwargs):
    """
    Uses previously performed Fr-TMalign fasta to get the sequence alignment and residue mapping between files
    """
    frtmalign_path = kwargs.get('frtmalign_path', "frtmalign")
    #flag = extract_frtmalign_fasta(file, wkdir + "frtm.fasta", target, template)
    flag = extract_frtm_alignment(file, wkdir + "frtm.fasta", target, template)
    if flag != 2:
        print("Fr-TMalign alignment doesn't exist. Attempting alignment ...")
        print(wkdir)
        frtmaln(wkdir, oriented_template, oriented_target, frtmalign_path)
        print("Alignmnet completed")

    tmp = SeqIO.to_dict(SeqIO.parse(wkdir + "frtm.fasta", "fasta"))
    # print(tmp)
    aligned_A = tmp[template].seq
    aligned_B = tmp[target].seq
    os.remove(wkdir + "frtm.fasta")

    resseq_A = get_pdb_sequence_with_chains(structA)
    resseq_B = get_pdb_sequence_with_chains(structB)

    # print(aligned_A)
    # print(resseq_A)
    # print(aligned_B)
    # print(resseq_B)
    # Equivalent residue numbering
    # Relative to reference
    mapping = {}
    aa_i_A, aa_i_B = 0, 0
    for aln_i, (aa_aln_A, aa_aln_B) in enumerate(zip(aligned_A, aligned_B)):
        if aa_aln_A == '-':
            if aa_aln_B != '-':
                aa_i_B += 1
        elif aa_aln_B == '-':
            if aa_aln_A != '-':
                aa_i_A += 1
        else:
            # print(resseq_A[aa_i_A], aa_aln_A, " ",end="")
            # print(resseq_B[aa_i_B], aa_aln_B)
            curr_aa_i_A = aa_i_A
            curr_aa_i_B = aa_i_B
            if resseq_A[aa_i_A][1] != aa_aln_A and (
                    resseq_A[aa_i_A][1] == aligned_A[(aln_i + 1) % (len(aligned_A) - 1)]):
                aa_i_B += 1
                print('pass')
                continue
            while (aa_i_A < len(resseq_A) and ((curr_aa_i_A + 2) > aa_i_A or resseq_A[aa_i_A][1] == 'M')) and \
                    resseq_A[aa_i_A][1] != aa_aln_A:
                print('skipping residue ', resseq_A[aa_i_A], 'from template')
                aa_i_A += 1
                print(aa_i_A, aa_i_B, " counters")
            if (aa_i_A == len(resseq_A) or (curr_aa_i_A + 2) == aa_i_A) and resseq_A[aa_i_A][1] != aa_aln_A:
                aa_i_A = curr_aa_i_A + 1
                aa_i_B += 1
                print(aa_i_A, aa_i_B, " counters")
                continue
            if resseq_B[aa_i_B][1] != aa_aln_B and (
                    resseq_B[aa_i_B][1] == aligned_B[(aln_i + 1) % (len(aligned_B) - 1)]):
                aa_i_A += 1
                print('pass')
                continue
            while (aa_i_B < len(resseq_B) and ((curr_aa_i_B + 2) > aa_i_B or resseq_B[aa_i_B][1] == 'M')) and \
                    resseq_B[aa_i_B][1] != aa_aln_B:
                print('skipping residue ', resseq_B[aa_i_B], 'from target')
                aa_i_B += 1
                print(aa_i_A, aa_i_B, " counters")
            if (aa_i_B == len(resseq_B) or (curr_aa_i_B + 2) == aa_i_B) and resseq_B[aa_i_B][1] != aa_aln_B:
                aa_i_B = curr_aa_i_B
                aa_i_A += 1
                print(aa_i_A, aa_i_B, " counters")
                continue

            assert resseq_A[aa_i_A][1] == aa_aln_A
            assert resseq_B[aa_i_B][1] == aa_aln_B
            mapping[(resseq_A[aa_i_A][0], resseq_A[aa_i_A][2], resseq_A[aa_i_A][3], resseq_A[aa_i_A][4])] = (
            resseq_B[aa_i_B][0], resseq_B[aa_i_B][2], resseq_B[aa_i_B][3], resseq_B[aa_i_B][4])
            aa_i_A += 1
            aa_i_B += 1

    # Calculate (gapless) sequence identity
    seq_id, g_seq_id = calculate_identity(aligned_A, aligned_B)
    return ((aligned_A, aligned_B), seq_id, g_seq_id, mapping)


def align_sequences_with_chains_and_gaps(wkdir, structA, structB, **kwargs):
    """
    Performs a global pairwise alignment between two sequences
    using the BLOSUM62 matrix and the Needleman-Wunsch algorithm
    as implemented in Biopython. Returns the alignment and the 
    residue mapping between both original sequences. Gaps are also
    included in the resude map.
    """

    matrix = kwargs.get('matrix', matlist.blosum62)
    gap_open = kwargs.get('gap_open', -10.0)
    gap_extend = kwargs.get('gap_extend', -0.5)
    trim_ends = kwargs.get('trim_ends', True)
    muscle_path = kwargs.get('muscle_path', "muscle")

    resseq_A = get_pdb_sequence_with_chains(structA)
    resseq_B = get_pdb_sequence_with_chains(structB)

    sequence_A = ''.join([i[1] for i in resseq_A])
    sequence_B = ''.join([i[1] for i in resseq_B])

    try:
        if os.path.isfile(wkdir + "out.fasta"):
            os.remove(wkdir + "out.fasta")
        if os.path.isfile(wkdir + "tmp.fasta"):
            os.remove(wkdir + "tmp.fasta")
        f = open(wkdir + "tmp.fasta", "w")
        f.write(">seqA\n")
        f.write(sequence_A + "\n")
        f.write(">seqB\n")
        f.write(sequence_B + "\n")
        f.close()
        fnull = open(os.devnull, 'w')
        p = subprocess.Popen(muscle_path.split() + ['-in', wkdir + "tmp.fasta", '-out', wkdir + "out.fasta"],
                             stdout=fnull, stderr=fnull)
        p.wait()
        fnull.close()
        tmp = SeqIO.to_dict(SeqIO.parse(wkdir + "out.fasta", "fasta"))
        aligned_A = tmp['seqA'].seq
        aligned_B = tmp['seqB'].seq
        os.remove(wkdir + "out.fasta")
        os.remove(wkdir + "tmp.fasta")
    except:
        print("Muscle procedure failed. Reverting to BioPython...")
        alns = pairwise2.align.globalds(sequence_A, sequence_B,
                                        matrix, gap_open, gap_extend,
                                        penalize_end_gaps=(False, False))

        best_aln = alns[0]
        aligned_A, aligned_B, score, begin, end = best_aln

    mapping = [[], []]
    aa_i_A, aa_i_B = 0, 0
    for aln_i, (aa_aln_A, aa_aln_B) in enumerate(zip(aligned_A, aligned_B)):
        if aa_aln_A == '-':
            if aa_aln_B != '-':
                mapping[0].append(('-', '-'))
                mapping[1].append((resseq_B[aa_i_B][0], resseq_B[aa_i_B][2]))
                aa_i_B += 1
        elif aa_aln_B == '-':
            if aa_aln_A != '-':
                mapping[0].append((resseq_A[aa_i_A][0], resseq_A[aa_i_A][2]))
                mapping[1].append(('-', '-'))
                aa_i_A += 1
        else:
            assert resseq_A[aa_i_A][1] == aa_aln_A
            assert resseq_B[aa_i_B][1] == aa_aln_B
            mapping[0].append((resseq_A[aa_i_A][0], resseq_A[aa_i_A][2]))
            mapping[1].append((resseq_B[aa_i_B][0], resseq_B[aa_i_B][2]))
            aa_i_A += 1
            aa_i_B += 1

    return ((aligned_A, aligned_B), mapping)


def parse_structure(spath):
    """Parses a PDB/cif structure"""

    if not os.path.isfile(spath):
        return IOError('File not found: {0}'.format(spath))

    if spath.endswith(('pdb', 'ent')):
        parser = PDBParser(QUIET=True)
    elif spath.endswith('cif'):
        parser = FastMMCIFParser()
    else:
        raise Exception('Format not supported ({0}). Must be .pdb/.ent or .cif'.format(spath))

    sname = os.path.basename(spath.split('.')[0])
    return parser.get_structure(sname, spath)


def cesymm_related_repeats(inputf, cesymm_dir):
    axes = open(cesymm_dir + inputf + ".axes", "r")
    flag = 0
    repeats_all = []
    repeats_type = []
    repeats_level = []
    axes_per_level = []
    order_of_level = []
    for line in axes:
        if inputf[0:4] in line:
            flag = 1
            flds = line.split()
            order = int(flds[1])
            #            repeats_all.append([])
            repeats = flds[8].strip(')/(/\n')
            if ')(' in repeats:
                repeats = repeats.split(')(')
                repeats = [k.split(';') for k in repeats]
            else:
                repeats = [repeats.split(';')]
            if len(repeats_all) < order:
                repeats_all.append(repeats)
                repeats_type.append(flds[2])
                repeats_level.append(order)
                axes_per_level.append(1)
                order_of_level.append(flds[3])
            else:
                for k in repeats:
                    repeats_all[order - 1].append(k)
                axes_per_level[order - 1] += 1
                if flds[2] != repeats_type[order - 1]:
                    print(("Warning: Combining %s with %s symmetries!") % (repeats_type[order - 1], flds[2]))
    if flag == 1:
        return repeats_all, repeats_type, repeats_level, axes_per_level, order_of_level
    else:
        raise SystemExit("Error: %s has no detected symmetry\n" % inputf)


def cesymm_levels(inputf, cesymm_dir):
    axes = open(cesymm_dir + inputf + ".axes", "r")
    flag = 0
    level_type = []  # Open/Closed
    levels = []  # 1, 2, ...
    axes_per_level = []
    order_of_level = []  # C2, H2, ...
    level_symm_type = []  # Quaternal/Internal

    for line in axes:
        if inputf[0:4] in line:
            flag = 1
            flds = line.split()
            level = int(flds[1])
            #            repeats_all.append([])
            if level not in levels:
                rep = flds[8].split(")")
                chain = [rep[0][pos + 1] for pos, char in enumerate(rep[0]) if
                         char == '.']  # finds all chain names in first symmetry bracket (they are preceded by .)
                if chain.count(chain[0]) != len(
                        chain):  # if all repeats are from the same chain, no quaternary structure symmetry was detected
                    level_symm_type.append("Quaternary")
                else:
                    level_symm_type.append("Internal")
                level_type.append(flds[2])
                levels.append(level)
                axes_per_level.append(1)
                order_of_level.append(flds[3])

            else:
                axes_per_level[level - 1] += 1
                if flds[2] != level_type[level - 1]:
                    print("Warning: Combining %s with %s symmetries!" % (repeats_type[order - 1], flds[2]))

    if flag == 1:
        return levels, level_type, axes_per_level, order_of_level, level_symm_type
    else:
        raise SystemExit("Error: %s has no detected symmetry\n" % inputf)


def find_tm_chains(inputf, opm_dir, flag):
    # Turn flag to 0 to use Edo's dictionary str_info; Turn to 1 when the structure is not part of Edo's dictionary
    if flag == 0:
        str_info = read_data(opm_dir + "opm_archive.txt")
        tm_ch = str_info[inputf]['tmchains']
    else:
        chain = ""
        #    if downloaded==1:
        #      chain=chain_id(oriented)
        ppm_results = open(opm_dir + "all_opm_ppm_results_corrected.dat", "r")
        for line in ppm_results:
            if line.startswith(inputf):
                flds = line.split()
                chain = flds[1]
        ppm_results.close()
        if chain == "":
            tm_ch = []
            raise SystemExit(inputf + " has no known TM chains.")
        else:
            tm_ch = chain.split(';')[:-1]
    return (tm_ch)


def strip_tm_chains(wkdir, inputf, oriented_enc, chains):
    """
    Use for transfer of symmetry
    Extract only the transmembrane chains atoms that are properly specified and order them as they appear in the PDB file.
    Note that chains can be either an array/list or a string.
    """
    print(oriented_enc, chains)
    f = open(oriented_enc, 'r')
    altloc = ' '
    flag = 0
    # Choose the altloc with highest occupancy if multiple present
    for line in f:
        if line.startswith("ATOM") and line[12:16].strip() == 'CA' and line[16:17] != ' ' and (
                float(line[54:60]) > 0.5 or flag == 0):
            altloc = line[16:17]
            flag = 1
    f.seek(0)
    # print(altloc)
    print(wkdir + inputf + "_tmp.pdb")
    o = open(wkdir + inputf + "_tmp.pdb", "w")
    num = 0
    LINELEM = "{:76s}{:>2s}\n"
    atomnames = []
    old_resid = '-88'
    for line in f:
        if (line.startswith("ATOM") or line.startswith("TER ")) and (
                line[21:22] in chains or chains == "") and from3to1_general(line[17:20].strip()) != 0:
            if old_resid != '-88' and line[22:27] == old_resid and line[12:16] in atomnames or (
                    line[16:17] != altloc and line[16:17] != ' '):  # sort out disordered atoms
                continue
            elif (old_resid == '-88' or line[22:27] != old_resid):
                old_resid = line[22:27]
                atomnames = []
                # print(line)
            atomnames.append(line[12:16])

            if line[76:78].strip() != '':  # ensure that the element symbol is included
                o.write(line)
                num += 1
            else:
                if line.startswith("TER "):
                    o.write(line)
                else:
                    atom = line[12:16].strip()
                    elem = atom[0]
                    o.write(LINELEM.format(line[0:76], elem))
                    num += 1
        if line.startswith("HETATM") and line[17:20] == 'MSE' and (line[21:22] in chains or chains == ""):
            if old_resid != '-88' and line[22:27] == old_resid and line[12:16] in atomnames or (
                    line[16:17] != altloc and line[16:17] != ' '):
                continue
            elif (old_resid == '-88' or line[22:27] != old_resid):
                old_resid = line[22:27]
                atomnames = []
            atomnames.append(line[12:16])

            if line[76:78].strip() != '':
                o.write("ATOM  " + line[6:])
                num += 1
            else:
                atom = line[12:16].strip()
                elem = atom[0]
                o.write("ATOM  " + line[6:76] + " " + elem + "\n")
                num += 1

    o.write("END\n")
    f.close()
    o.close()
    return num


def strip_tm_chains_in_order(wkdir, inputf, oriented_enc, chains):
    """
    Extract only the transmembrane chains atoms that are properly specified and order the chains as requested in the variable chains.
    Note that chains is an array here, eg. ['A','B','C'].
    """
    f = open(oriented_enc, 'r')
    altloc = ' '
    flag = 0
    for line in f:
        if line.startswith("ATOM") and line[12:16].strip() == 'CA' and line[16:17] != ' ' and (
                float(line[54:60]) > 0.5 or flag == 0):
            altloc = line[16:17]
    f.seek(0)
    o = open(wkdir + inputf + "_tmp.pdb", "w")
    num = 0
    LINELEM = "{:76s}{:>2s}\n"
    for chain in chains:
        f.seek(0)
        atomnames = []
        old_resid = '-88'
        for line in f:
            if (line.startswith("ATOM") or line.startswith("TER ")) and line[21:22] == chain and from3to1_general(
                    line[17:20].strip()) != 0:
                if old_resid != '-88' and line[22:27] == old_resid and line[12:16] in atomnames or (
                        line[16:17] != altloc and line[16:17] != ' '):  # sort out disordered atoms
                    continue
                elif (old_resid == '-88' or line[22:27] != old_resid):
                    old_resid = line[22:27]
                    atomnames = []
                atomnames.append(line[12:16])

                if line[76:78].strip() != '':  # ensure that the element symbol is included
                    o.write(line)
                    num += 1
                else:
                    if line.startswith("TER "):
                        o.write(line)
                    else:
                        atom = line[12:16].strip()
                        elem = atom[0]
                        o.write(LINELEM.format(line[0:76], elem))
                        num += 1
            if line.startswith("HETATM") and line[17:20] == 'MSE' and line[21:22] == chain:
                if old_resid != '-88' and line[22:27] == old_resid and line[12:16] in atomnames or (
                        line[16:17] != altloc and line[16:17] != ' '):
                    continue
                elif (old_resid == '-88' or line[22:27] != old_resid):
                    old_resid = line[22:27]
                    atomnames = []
                atomnames.append(line[12:16])

                if line[76:78].strip() != '':
                    o.write("ATOM  " + line[6:])
                    num += 1
                else:
                    atom = line[12:16].strip()
                    elem = atom[0]
                    o.write("ATOM  " + line[6:76] + " " + elem + "\n")
                    num += 1

    o.write("END\n")
    f.close()
    o.close()
    return num


def cesymm_last_resid_repeat(inputf, oriented, std_file):
    """
    Finds the last residue that has been aligned in a single chain. Returns the difference 
    between that residue and the last residue in the chain. This is useful if we want to 
    know whether a significant portion of the protein is not aligned and might contain an 
    independent symmetry.
    """
    stdout = open(std_file, "r")
    num_rep = 0
    order = 'C1'
    sym_level = 0
    for line in stdout:
        if inputf in line:
            flds = line.split()
            num_rep = int(flds[1])
            order = flds[2].strip()
            try:
                sym_level = int(flds[4])
            except ValueError:
                sym_level = int(flds[3])
                order = ""
            try:
                repeats = flds[16]
            except:
                repeats = ""
            break
    stdout.close()
    if num_rep > 1 and order != 'C1' and sym_level > 0 and repeats != "":
        flds = repeats.split(';')
        last_resids = []
        for entry in flds:
            res = entry.split('-')
            last_resids.append(int(res[1]))
        last_res_id = max(last_resids)
        f = open(oriented, 'r')
        res_num = 0
        old_resid = -1
        last_id = 0
        for line in f:
            if line.startswith("ATOM") and line[13:16].strip() == "CA" and int(
                    re.search(r'\d+', line[22:27]).group()) > last_res_id and int(
                    re.search(r'\d+', line[22:27]).group()) != old_resid:
                res_num += 1
                old_resid = int(re.search(r'\d+', line[22:27]).group())
    else:
        res_num = 0
        old_resid = 0
        last_res_id = 0
    return res_num, old_resid, last_res_id


def symd_chain_order(wkdir, symd_out_dir, inputf, oriented, chains):
    """
    We use the results from a SymD run to order the chains in a multi-chain complex
    * chains is a string such as 'A;B;C;'
    """
    s_template = parse_structure(oriented)

    resseq_A = get_pdb_sequence_with_chains(
        s_template)  # creates a list of the form [(2, 'D', 'A'), (3, 'Q', 'A'), (4, 'A', 'A'), (5, 'L', 'A'), (6, 'S', 'A')...]
    symd_fasta = SeqIO.to_dict(SeqIO.parse(symd_out_dir + inputf + "-best.fasta", "fasta"))

    # Find the initial shift
    flds = symd_fasta[inputf + '-original'].description.split('=')[1]
    init_shift = int(flds.split('(')[0])
    init_shift = (len(resseq_A) + init_shift) % len(resseq_A)

    # Re-create the SymD alignments with the chain names instead of the amino acid names
    i = 0
    orig = ''
    for c in symd_fasta[inputf + '-original'].seq:
        if c != '-':
            orig += resseq_A[i][2]
            i += 1
        if c == '-':
            orig += c

    i = 0
    perm = ''
    for c in symd_fasta[inputf + '-permuted'].seq:
        if c != '-':
            ind = (i + init_shift) % len(resseq_A)
            perm += resseq_A[ind][2]
            i += 1
        if c == '-':
            perm += c

    # Determine how many chains we should have in each repeat unit
    symd_info = open(symd_out_dir + inputf + "-info.txt", "r")
    for line in symd_info:
        if line.startswith("! Derived unit angle:"):
            angle = float(line.split(':')[1])
            break
    symd_info.close()

    if angle != 0:
        order = round(360 / angle)
        chain_unit = len(chains) // (2 * order)
    else:
        chain_unit = len(chains)

    # For each chain find out where it is supposed to go based on the frequency of alignments with another chain
    ch = []  # list of chains
    aligned_ch = []  # list of aligned chains to the ch list
    counts = []  # how many times that particular alignment is encountered

    for i, c in enumerate(orig):
        if c != '-' and symd_fasta[inputf + '-original'].seq[i].isupper() and c not in ch:
            ch.append(c)
            ind = ch.index(c)
            if perm[i] != '-':
                aligned_ch.append([perm[i]])
                counts.append([1])
            else:
                aligned_ch.append([])
                counts.append([0])
        elif c != '-' and perm[i] != '-' and symd_fasta[inputf + '-permuted'].seq[i].isupper() and perm[i] not in \
                aligned_ch[ind]:
            aligned_ch[ind].append(perm[i])
            counts[ind].append(1)
        elif c != '-' and perm[i] != '-' and symd_fasta[inputf + '-permuted'].seq[i].isupper():
            subind = aligned_ch[ind].index(perm[i])
            counts[ind][subind] += 1

    # Choose the corresponding chain to each chain based on high frequency
    matches = []
    for i, j in enumerate(counts):
        index, value = max(enumerate(j), key=operator.itemgetter(1))
        if value != 0 and value > 5:
            matches.append(aligned_ch[i][index])
        else:
            matches.append("")

    # Form a list of pairs of connected chains
    lst = []
    for i, c in enumerate(ch):
        if matches[i] != "":
            lst.append([c, matches[i]])
        else:
            lst.append([c])
    #    limit=len(lst)
    limit = 10

    # Link all connected chains
    count = 0
    while len(lst) != chain_unit and count < limit:
        count += 1
        for i, elem in enumerate(lst):
            for k in elem:
                if len(lst) < i + 1:
                    break
                for ind, j in enumerate(lst):
                    if ind != i and k in j:
                        pos = j.index(k)
                        elem = elem + j[(pos + 1):len(j)]
                        elem = j[0:pos] + elem
                        if i < ind:
                            lst[i] = elem
                            del lst[ind]
                        else:
                            lst[ind] = elem
                        break

    # Remove repeating chains from a sublist that often occur in a cyclical symmetry
    for i, elem in enumerate(lst):
        seen = set()
        lst[i] = [x for x in elem if x not in seen and not seen.add(x)]

    # Remove remaining chains that have no partner (single chains)
    refined = []
    for i, elem in enumerate(lst):
        if len(elem) != 1:
            refined.append(elem)
    lst = refined

    if len(lst) != 0:
        # Format the final order in which the chains should appear in the pdb
        n = len(lst[0])
        order = []
        if np.all([len(x) == n for x in lst]) == True:
            for i in range(0, n):
                for item in lst:
                    if item[i] not in order:
                        order = order + [item[i]]
            order_pdb_chains(wkdir + inputf + "_symd", oriented, order)
            os.rename(wkdir + inputf + "_symd_ordered.pdb", wkdir + inputf + "_sd_order.pdb")
            return ';'.join(order) + ';'

        else:
            print("The order is inconsistent: ", lst)
            return ""


    else:
        print("SymD detected slip symmetry.")
        return ""

    #### Reading CE-Symm and SymD information ####


def cesymm_data(cesymm_out_dir, inputf):
    """
    Crawl through the CE-Symm result files and collect all the necessary information for a structure
    """
    ce_order = "na"
    sym_level = "na"
    ce_rmsd = "na"
    tmScore = "na"
    cce_type = "na"
    num_repeats = "na"
    internal = "na"
    quaternary = "na"
    angle = "na"
    axis_angle = "na"
    detected_symmetry = "na"
    coverage = "na"
    length = "na"
    core_length = "na"
    translation = "na"
    unrefined_tmscore = "na"
    unrefined_rmsd = "na"
    seed = "na"
    symm_type = "na"
    repeats = "na"
    repeat_length = "na"
    closed_opened = "na"
    repeat_level = "na"
    if os.path.isfile(cesymm_out_dir + inputf + ".axes") and os.path.isfile(
            cesymm_out_dir + inputf + "_stdout.out") and os.path.getsize(cesymm_out_dir + inputf + "_stdout.out") > 5:
        stdout = open(cesymm_out_dir + inputf + "_stdout.out", "r")
        for ln in stdout:
            if inputf[0:4] in ln:
                fields = ln.split()
                num_repeats = fields[1]  # 17
                ce_order = fields[2]  # 18
                ind = 0
                try:
                    sym_level = int(fields[4])
                except ValueError:
                    ind = -1
                    sym_level = int(fields[ind + 4])  # 20
                    ce_order = 'na'
                closed_opened = fields[ind + 5]
                tmScore = float(fields[ind + 10])  # 24
                ce_rmsd = float(fields[ind + 11])  # 25
                unrefined_tmscore = float(fields[ind + 8])
                unrefined_rmsd = float(fields[ind + 9])
                repeat_length = int(fields[ind + 12])
                core_length = int(fields[ind + 13])  # 29
                length = int(fields[ind + 14])  # 30
                coverage = fields[ind + 15]  # 31
                if length > 0 and float(
                        coverage) - core_length / length > 0.01:  # correct errors in CE-Symm coverage calculation
                    coverage = str(round(core_length / length, 2))
                cce_type = "na"
                angle = fields[ind + 6]
                internal = "No"
                quaternary = "No"
                axis_angle = "na"
                detected_symmetry = "na"
                translation = fields[ind + 7]
                with open(cesymm_out_dir + inputf + ".seed", "r") as f:
                    seed = f.readline().strip()
                symm_type = "na"
                repeats = "na"
                repeat_level = "na"
                if sym_level > 0:
                    cce_type = ""
                    angle = ""
                    axis_angle = ""
                    translation = ""
                    symm_type = ""
                    repeats = ""
                    closed_opened = ""
                    repeat_level = ""
                    f = open(cesymm_out_dir + inputf + ".axes", "r")
                    for line in f:
                        if inputf[0:4] in line:
                            flds = line.split()
                            closed_opened = closed_opened + flds[2] + ';'
                            angle = angle + flds[4] + ";"
                            translation = translation + flds[5] + ";"
                            repeats = repeats + flds[8].replace(";", ",") + ";"
                            repeat_level = repeat_level + flds[1] + ";"
                            pt1 = [float(x) for x in flds[6].split(",")]
                            pt1 = np.array(pt1)
                            pt2 = [float(x) for x in flds[7].split(",")]
                            pt2 = np.array(pt2)
                            u_tmp = pt1 - pt2
                            u = unit_vector(u_tmp)
                            axis_angle_tmp = angle_between(u, (0, 0, 1)) * 180 / np.pi
                            if axis_angle_tmp > 90:
                                axis_angle_tmp = axis_angle_tmp - 180
                            axis_angle = axis_angle + str(float('%.2f' % round(axis_angle_tmp, 2))) + ';'
                            if np.abs(np.dot(u, (0, 0, 1))) > 0.5:
                                cce_type = cce_type + "Parallel;"
                            else:
                                cce_type = cce_type + "Antiparallel;"
                            rep = flds[8].split(")")
                            chain = [rep[0][pos + 1] for pos, char in enumerate(rep[0]) if
                                     char == '.']  # finds all chain names in first symmetry bracket (they are preceded by .)
                            if chain.count(chain[0]) != len(
                                    chain):  # if all repeats are from the same chain, no quaternary structure symmetry was detected
                                detected_symmetry = "Quaternary"
                                quaternary = 'Yes'
                            else:
                                detected_symmetry = "Internal"
                                internal = 'Yes'
                            symm_type = symm_type + detected_symmetry + ";"
                    f.close()
                    repeats = repeats.replace(inputf[0:4].upper() + ".", "")


    else:
        print("CE-Symm files missing for %s" % inputf)

    cesymm_dic = {}
    cesymm_dic['size'] = length
    cesymm_dic['coverage'] = coverage
    cesymm_dic['unit_angle'] = angle
    cesymm_dic['unit_translation'] = translation
    cesymm_dic['refined_tmscore'] = str(tmScore)
    cesymm_dic['refined_rmsd'] = str(ce_rmsd)
    cesymm_dic['unrefined_tmscore'] = str(unrefined_tmscore)
    cesymm_dic['unrefined_rmsd'] = str(unrefined_rmsd)
    cesymm_dic['repeats_number'] = num_repeats
    cesymm_dic['topology'] = cce_type
    cesymm_dic['symmetry_order'] = ce_order
    cesymm_dic['symmetry_levels'] = str(sym_level)
    cesymm_dic['aligned_length'] = core_length
    cesymm_dic['seed'] = seed
    cesymm_dic['internal_symmetry'] = internal
    cesymm_dic['quaternary_symmetry'] = quaternary
    cesymm_dic['axis_angle_with_membrane_normal'] = axis_angle
    cesymm_dic['symmetry_type'] = symm_type
    cesymm_dic['repeats'] = repeats
    cesymm_dic['repeat_length'] = repeat_length
    cesymm_dic['closed_opened'] = closed_opened
    cesymm_dic['repeat_level'] = repeat_level
    return cesymm_dic


def symd_data(symd_out_dir, inputf):
    """
    Crawl through the SymD result files and collect all the necessary information for a structure
    """
    if os.path.isfile(symd_out_dir + inputf + "-info.txt") and os.path.isfile(symd_out_dir + inputf + "-trfm.pdb"):
        f = open(symd_out_dir + inputf + '-info.txt', 'r')
        for line in f:
            if line.startswith("! Protein Size:"):
                flds = line.split()
                size = int(flds[3])
            if line.startswith("! Best IS:"):
                flds = line.split(':')
                best_is = flds[1].strip()
            if line.startswith("! N-aligned at best IS:"):
                flds = line.split(':')
                aligned = int(flds[1])
                align_per = str(round(1.0 * aligned / size, 2))
            if line.startswith("! TM/Nr at best IS:"):
                flds = line.split(':')
                symd_tm = float(flds[1])
                symd_tm = str(round(symd_tm, 2))
            if line.startswith("! Z(TsC), Z(Ts), Z(TM) at best IS:"):
                flds = line.split(',')
                z_tm = flds[-1].strip()
            if line.startswith("! RMSD at best IS:"):
                flds = line.split(':')
                symd_rmsd = float(flds[1])
                symd_rmsd = str(round(symd_rmsd, 3))
            if line.startswith("! Highest-scoring angle:"):
                flds = line.split()
                symd_angle = float(flds[3])
            if line.startswith("! Highest-scoring p-transl:"):
                flds = line.split()
                symd_translation = flds[3].strip()
            if line.startswith("! Derived unit angle:"):
                flds = line.split()
                symd_unit_angle = float(flds[4])
                if symd_unit_angle != 0:
                    symd_order = int(round(360.0 / symd_unit_angle))
                else:
                    symd_order = "na"
            if line.startswith("! Derived unit p-transl:"):
                flds = line.split()
                symd_unit_translation = flds[4].strip()
                break
        f.close()
        f = open(symd_out_dir + inputf + '-trfm.pdb', 'r')
        counter = 0
        flag = 0
        pts = []
        for line in f:
            if line.startswith("MODEL        3"):
                flag = 1
            if flag == 1 and 0 <= counter < 2 and line.startswith("ATOM") and line[13:16].strip() == 'CA':
                pts.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                counter += 1
        f.close()
        if len(pts) == 2:
            u = np.array(pts[0]) - np.array(pts[1])
            u = unit_vector(u)
            axis_angle = angle_between(u, (0, 0, 1)) * 180 / np.pi
            if axis_angle > 90:
                axis_angle = axis_angle - 180
            axis_angle = str(float('%.2f' % axis_angle))
            if np.abs(np.dot(u, (0, 0, 1))) > 0.5:
                symd_type = "orthogonal"
            else:
                symd_type = "parallel"
        else:
            symd_order = symd_rmsd = symd_angle = symd_type = align_per = aligned = symd_tm = z_tm = symd_unit_translation = \
                symd_translation = symd_unit_angle = size = best_is = axis_angle = "na"

    else:
        print("Missing SymD files for %s" % inputf)
        symd_order = "na"
        symd_rmsd = "na"
        symd_angle = "na"
        symd_type = "na"
        align_per = "na"
        aligned = "na"
        symd_tm = "na"
        z_tm = "na"
        symd_unit_translation = "na"
        symd_translation = "na"
        symd_unit_angle = "na"
        size = "na"
        best_is = "na"
        axis_angle = "na"

    symd_dic = {}
    symd_dic['size'] = str(size)
    symd_dic['coverage'] = str(align_per)
    symd_dic['unit_angle'] = str(symd_unit_angle)
    symd_dic['unit_translation'] = symd_unit_translation
    symd_dic['tmscore'] = symd_tm
    symd_dic['rmsd'] = symd_rmsd
    symd_dic['initial_shift'] = best_is
    symd_dic['z_tmscore'] = z_tm
    symd_dic['is_angle'] = str(symd_angle)
    symd_dic['is_translation'] = symd_translation
    symd_dic['aligned_length'] = str(aligned)
    symd_dic['symmetry_order'] = str(symd_order)
    symd_dic['topology'] = symd_type
    symd_dic['axis_angle_with_membrane_normal'] = axis_angle
    return symd_dic


def ananas_data(out_dir, inputf, oriented_struct):
    fname = out_dir + inputf + "_ananas"
    ananas_dic = {}
    if os.path.isfile(fname + ".json") and len(open(fname + ".json").readlines()) > 2 and os.path.isfile(
            fname + ".out"):
        # Figure out the names and residue ranges of the chains and connect these to the AnAnaS numbering
        resseq_A = get_pdb_sequence_with_chains(
            oriented_struct)  # creates a list of the form [(2, 'D', 'A'), (3, 'Q', 'A'), (4, 'A', 'A'), (5, 'L', 'A'), (6, 'S', 'A')...]
        chain_resids = {}  # eg. {0: ('A', 2, 334), 1:('B', 5, 175)}
        chains = []
        [chains.append(i[2]) for i in resseq_A if i[2] not in chains]
        for i, chain in enumerate(chains):
            residues = list(filter(lambda x: x[2] == chain, resseq_A))
            chain_resids[i] = (chain, residues[0], residues[-1])

        # Read info from AnAnaS
        info = open(fname + ".out", 'r')
        for line in info:
            if line.startswith("Number of chains read"):
                chains_num = int(line.strip().split()[-1])
        info.close()
        with open(fname + ".json") as jf:
            data = json.load(jf)
        for i, symm in enumerate(data):
            if 'c' in symm['group'] and (int(symm['group'][1:]) > chains_num or (
                    int(symm['group'][1:]) < chains_num and chains_num % int(
                    symm['group'][1:]) != 0)):  # e.g. c9 to indicate symmetry group C9
                print('WARNING: Inconsistent chain number or cyclic order')
                continue
            elif 'd' in symm['group'] and (int(symm['group'][1:]) * 2 > chains_num or (
                    int(symm['group'][1:]) * 2 < chains_num and chains_num % int(
                    symm['group'][1:]) * 2 != 0)):  # e.g. d6 to indicate symmetry group D6
                print('WARNING: Inconsistent chain number or dihedral order')
                continue
            ananas_dic['symmetry_order'] = symm['group'].upper()  # the global order (e.g. D6)
            ananas_dic['avg_rmsd'] = symm['Average_RMSD']
            ananas_dic['symmetries'] = []

            # Find one (the first) highest order permutation
            for j, t in enumerate(symm['transforms']):
                level = {}
                if t['ORDER'] == int(ananas_dic['symmetry_order'][1:]):
                    level['order'] = t['ORDER']

                    # estimate topology and angle with membrane
                    u = np.array(t['P0']) - np.array(t['P1'])
                    level['topology'], axis_angle = determine_repeat_topology(u)

                    # Trace which chain corresponds to which and form the symmetry-related tuples
                    repeats_tuples = []  # eg. [[0,2],[1,3]]
                    for m in range(len(t['permutation'])):
                        tup = []
                        if len([True for k in repeats_tuples if (m in k or t['permutation'][m] in k)]) == 0:
                            tup.append(m)
                            ind = t['permutation'][m]
                            while ind not in tup:
                                tup.append(ind)
                                ind = t['permutation'][ind]
                            repeats_tuples.append(tup)
                    level['repeats_chain_indexes'] = np.transpose(repeats_tuples)
                    level['repeats'] = [
                        [(chain_resids[m][0], chain_resids[m][1][0], chain_resids[m][2][0]) for m in ind] for ind in
                        level[
                            'repeats_chain_indexes']]  # eg. [[('A', 1, 453), ('B', 44, 295)], [('C', 1, 453), ('D', 44, 295)]], i.e. [Rep1, Rep2]
                    repeats_text = [[chain_resids[m][0] + '_' + str(chain_resids[m][1][0])
                                     + '-' + str(chain_resids[m][2][0]) for m in ind] for ind in repeats_tuples]
                    level['repeats_text'] = ''.join(['(' + ','.join(m) + ')' for m in
                                                     repeats_text])  # CE-Symm compatible format, i.e. (A_1-453,C_1-453)(B_44-295,D_44-295)

                    level['symm_units'] = 'Quaternary'
                    level['group'] = 'Cyclic'
                    level['permutation'] = t['permutation']
                    ananas_dic['symmetries'].append(level)
                    break
            # for dihedral symmetries, extract the two-fold symmetries
            if ananas_dic['symmetry_order'][0] == 'D':
                chains_per_repeat = chains_num // (int(ananas_dic['symmetry_order'][1:]) * 2)
                two_fold_symms = np.array(ananas_dic['symmetries'][0]['repeats_chain_indexes'])
                two_fold_symms = two_fold_symms.reshape(two_fold_symms.shape[0],
                                                        two_fold_symms.shape[1] // chains_per_repeat, chains_per_repeat)
                for correspondence in two_fold_symms:
                    for t in symm['transforms']:
                        if [t['permutation'][m] for m in correspondence[0]] == list(correspondence[1]):
                            level = {}
                            level['order'] = t['ORDER']
                            assert level['order'] == len(correspondence)

                            # estimate topology and angle with membrane
                            u = np.array(t['P0']) - np.array(t['P1'])
                            level['topology'], axis_angle = determine_repeat_topology(u)
                            level['repeats'] = [
                                [(chain_resids[m][0], chain_resids[m][1][0], chain_resids[m][2][0]) for m in r] for r in
                                correspondence]
                            repeats_text = [[chain_resids[m][0] + '_' + str(chain_resids[m][1][0])
                                             + '-' + str(chain_resids[m][2][0]) for m in r] for r in correspondence]
                            repeats_text = np.transpose(repeats_text)
                            level['repeats_text'] = ''.join(
                                ['(' + ','.join(m) + ')' for m in repeats_text])  # CE-Symm compatible format
                            level['symm_units'] = 'Quaternary'
                            level['group'] = 'Cyclic'
                            level['permutation'] = t['permutation']
                            ananas_dic['symmetries'].append(level)
                            break
            break
        print(ananas_dic)

    else:
        if os.path.isfile(fname + ".json") and len(open(fname + ".json").readlines()) <= 2:
            print("AnAnaS found no symmetry %s" % inputf)
        else:
            print("Missing AnAnaS files for %s" % inputf)
    if len(ananas_dic) == 0:
        ananas_dic = {'symmetry_order': None, 'avg_rmsd': None, 'symmetries': [
            {'order': None, 'group': None, 'topology': None, 'repeats': None, 'symm_units': None}]}

    return ananas_dic


def determine_repeat_topology(u):
    u = unit_vector(u)
    axis_angle = angle_between(u, (0, 0, 1)) * 180 / np.pi
    if axis_angle > 90:
        axis_angle = axis_angle - 180
    if np.abs(np.dot(u, (0, 0, 1))) > 0.5:
        return "Parallel", str(float('%.2f' % axis_angle))
    else:
        return "Antiparallel", str(float('%.2f' % axis_angle))


def symd_quat_str(symd_out_dir, inputf, oriented, chain):
    # Determines whether SymD has found a quaternary or internal symmetry
    f = open(oriented, 'r')
    chain_len = []
    #	print inputf
    #	print chain
    all_chain = len(chain) // 2
    for i in range(0, all_chain):
        ch_id = chain[i * 2]
        #		print ch_id
        f.seek(0)
        num = 0
        old_resid = -1
        for line in f:
            if line.startswith("ATOM") and line[21:22].strip() == ch_id and line[12:16].strip() == 'CA' and int(
                    re.search(r'\d+', line[
                                      22:27]).group()) != old_resid:  # What happens in cases where there are two residues with 50% occupancy?
                num += 1
                old_resid = int(re.search(r'\d+', line[22:27]).group())
        chain_len.append(num)
    f.close()

    #	print chain_len
    alignment = open(symd_out_dir + inputf + "-best.fasta", "r")
    flag = 0
    sequence = ""
    for line in alignment:  # Extract alignement sequence from SymD
        if line.startswith(">") and flag == 1:
            flag = 2
        if line.startswith(">") and flag == 0:
            flag = 1
        if not line.startswith(">") and flag == 1:
            sequence = sequence + line.strip()
    sequence_clean = ""
    for c in sequence:  # Get rid of all the gaps
        if c != '-':
            sequence_clean = sequence_clean + c
    alignment.close()
    #	print len(sequence_clean)
    ### Check how many aligned residues there are in each chain
    big_let = 0
    total_big = 0  # records the total number of big letters in the sequence
    aligned_resid = []  # records the number of big letters in each chain
    count = -1  # counts how many residues before a given chain
    i = 0  # counts the chain index
    for ind, c in enumerate(sequence_clean):
        if c.isupper():
            big_let += 1
        if ind == count + chain_len[i]:
            count = ind
            #			print count
            i += 1
            total_big = total_big + big_let
            aligned_resid.append(big_let)
            big_let = 0
    if len(aligned_resid) != len(chain_len):
        aligned_resid.append(big_let)
        total_big = total_big + big_let
    #	print ind
    #	print aligned_resid
    #	print total_big
    #	print "###"
    if np.any([x > total_big - 6 for x in aligned_resid]) == True:
        quat = 'no'
    else:
        quat = 'yes'
    #	os.remove(dir+"/"+inputf+"_tmp.pdb")
    #	print quat
    return quat


def cesymm_alignment(cesymm_out_dir, inputf):
    """
    Read the alignments from CE-Symm and remove the // characters
    
    """
    aln = open(cesymm_out_dir + inputf + ".fasta", "r")
    cs_aln = ""
    for line in aln:
        if line.startswith('>'):
            cs_aln = cs_aln + '\t\t\t\t<Aln name="' + line.strip() + '" >\n'
        elif not line.startswith('//'):
            cs_aln = cs_aln + '\t\t\t\t\t<Seq> ' + line.strip() + ' </Seq>\n'
            cs_aln = cs_aln + '\t\t\t\t</Aln>\n'
    aln.close()
    return cs_aln


def symd_alignment(symd_out_dir, inputf):
    """
    Remove all the new line characteres from the -best.fasta files and read 
    each alignment as a string
    
    """
    aln = open(symd_out_dir + inputf + "-best.fasta", "r")
    seq = ""
    symd_aln = ""
    for line in aln:
        if line.startswith('>'):
            if seq != "":
                symd_aln = symd_aln + '\t\t\t\t\t<Seq> ' + seq + ' </Seq>\n'
                symd_aln = symd_aln + '\t\t\t\t</Aln>\n'
                seq = ""
            symd_aln = symd_aln + '\t\t\t\t<Aln name="' + line.strip('\n') + '" >\n'

        else:
            seq = seq + line.strip()
    symd_aln = symd_aln + '\t\t\t\t\t<Seq> ' + seq + ' </Seq>\n'
    symd_aln = symd_aln + '\t\t\t\t</Aln>\n'
    aln.close()
    return symd_aln


def rgb2hex(r, g, b):
    hex = "#{:02x}{:02x}{:02x}".format(r, g, b)
    return hex


def scaled_rgb2hex(r, g, b):
    r = int(255 * r)
    g = int(255 * g)
    b = int(255 * b)
    hex = "#{:02x}{:02x}{:02x}".format(r, g, b)
    return hex


def colorblind_palette():
    """
    Adapted from:
    https://pymolwiki.org/index.php/Colorblindfriendly
    """
    cb_colors = (("vermillion", (0.835, 0.369, 0.000), rgb2hex(213, 94, 0)),  # (213,  94,   0)
                 ("blue", (0.000, 0.447, 0.698), rgb2hex(0, 114, 178)),  # (  0, 114, 178)
                 ("yellow", (0.941, 0.894, 0.259), rgb2hex(240, 228, 66)),  # (240, 228,  66),
                 ("bluish_green", (0.000, 0.620, 0.451), rgb2hex(0, 158, 115)),  # (  0, 158, 115)
                 ("orange", (0.902, 0.624, 0.000), rgb2hex(230, 159, 0)),  # (230, 159,   0)
                 ("sky_blue", (0.337, 0.706, 0.914), rgb2hex(86, 180, 233)),  # ( 86, 180, 233)
                 ("reddish_purple", (0.800, 0.475, 0.655), rgb2hex(204, 121, 167)))  # (204, 121, 167)

    cb_names = []
    for i in cb_colors:
        cb_names.append("cb_" + i[0])
    return cb_colors, cb_names


def additional_colors_palette():
    pml_colors = (('density', (0.1, 0.1, 0.6), scaled_rgb2hex(0.1, 0.1, 0.6)),
                  ('firebrick', (0.698, 0.13, 0.13), scaled_rgb2hex(0.698, 0.13, 0.13)),
                  ('tv_yellow', (1.0, 1.0, 0.2), scaled_rgb2hex(1.0, 1.0, 0.2)),
                  ('purple', (0.75, 0.0, 0.75), scaled_rgb2hex(0.75, 0.0, 0.75)),
                  ('tv_orange', (1.0, 0.55, 0.15), scaled_rgb2hex(1.0, 0.55, 0.15)),
                  ('forest', (0.2, 0.6, 0.2), scaled_rgb2hex(0.2, 0.6, 0.2)),
                  ('deepteal', (0.1, 0.6, 0.6), scaled_rgb2hex(0.1, 0.6, 0.6)),
                  ('marine', (0.0, 0.5, 1.0), scaled_rgb2hex(0.0, 0.5, 1.0)),
                  ('chartreuse', (0.5, 1.0, 0.0), scaled_rgb2hex(0.5, 1.0, 0.0)),
                  ('raspberry', (0.7, 0.3, 0.4), scaled_rgb2hex(0.7, 0.3, 0.4)),
                  ('wheat', (0.99, 0.82, 0.65), scaled_rgb2hex(0.99, 0.82, 0.65)),
                  ('hotpink', (1.0, 0.0, 0.5), scaled_rgb2hex(1.0, 0.0, 0.5)),
                  ('aquamarine', (0.5, 1.0, 1.0), scaled_rgb2hex(0.5, 1.0, 1.0)),
                  ('brightorange', (1.0, 0.7, 0.2), scaled_rgb2hex(1.0, 0.7, 0.2)),
                  ('brown', (0.65, 0.32, 0.17), scaled_rgb2hex(0.65, 0.32, 0.17)),
                  ('green', (0.0, 1.0, 0.0), scaled_rgb2hex(0.0, 1.0, 0.0)),
                  ('red', (1.0, 0.0, 0.0), scaled_rgb2hex(1.0, 0.0, 0.0)),
                  ('lightpink', (1.0, 0.75, 0.87), scaled_rgb2hex(1.0, 0.75, 0.87)))

    pml_names = []
    for i in pml_colors:
        pml_names.append(i[0])
    return pml_colors, pml_names


def cesymm_images(wkdir, cesymm_out_dir, inputf, repeat_selection, oriented, reference, ray, pymol, pic_dir, jmol_dir,
                  image_key):
    pml_files = ""
    image_files = ""
    jmol_files = ""
    if os.path.isfile(cesymm_out_dir + inputf + ".axes"):
        f = open(cesymm_out_dir + inputf + ".axes", "r")
        sym_num = 0
        for line in f:
            if inputf[0:4] in line:
                sym_num += 1
                flds = line.split()
                pt1 = [float(x) for x in flds[6].split(",")]
                pt1 = np.array(pt1)
                pt2 = [float(x) for x in flds[7].split(",")]
                pt2 = np.array(pt2)
                u_tmp = pt1 - pt2
                u = u_tmp / np.linalg.norm(u_tmp)
                # Scale the axis
                atoms = Selection.unfold_entities(reference, 'A')
                flag = 0
                scale = 0
                while flag == 0:  # move a plane normal to the axis upwards and check if there are any atoms above the plane
                    ref_pt = pt1 + u * scale
                    scale = scale + 2
                    i = 0
                    for atom in atoms:
                        i += 1
                        vec = atom.get_coord() - ref_pt
                        if np.dot(u, vec) > 0:
                            break
                        if i == len(list(atoms)):
                            upper_limit = ref_pt
                            flag = 1
                flag = 0
                scale = 0
                while flag == 0:  # move a plane normal to the axis downwards and check if there are any atoms below the plane
                    ref_pt = pt1 + u * scale
                    scale = scale - 2
                    i = 0
                    for atom in atoms:
                        i += 1
                        vec = atom.get_coord() - ref_pt
                        if np.dot(u, vec) < 0:
                            break
                        if i == len(list(atoms)):
                            lower_limit = ref_pt
                            flag = 1

                repeats = flds[8].strip(')/(/\n')
                repeats = repeats.split(')(')
                for j, ent in enumerate(repeats):
                    repeats[j] = ent.split(';')

                ### Create PyMol session
                cb_colors, cb_names = colorblind_palette()
                pml_colors, pml_names = additional_colors_palette()
                colors = cb_names + pml_names
                out = open(wkdir + "pymol_script_" + image_key + ".pml", "w")
                out.write("load " + oriented + ", " + inputf + "\n")
                out.write(
                    "set antialias,1\nset line_smooth, 1\nset depth_cue, 0.5\nset specular, 0\nset surface_quality, 1\n")
                out.write("set cartoon_sampling, 14\nset ribbon_sampling, 10\nset ray_trace_fog, 0.5\n")
                out.write("viewport 350,350\n")
                out.write("set ignore_case, off\n")
                for c in cb_colors:
                    out.write(
                        "set_color cb_" + c[0] + "=[" + str(c[1][0]) + "," + str(c[1][1]) + "," + str(c[1][2]) + "]\n")
                out.write("hide everything,*\n")
                out.write("show cartoon\n")
                out.write("color grey70, " + inputf + "\n")
                out.write("set cartoon_transparency, 0.5, " + inputf + "\n")
                out.write("bg_color white\nset opaque_background, on\n")
                if np.any([s != s for s in upper_limit]) == False or np.any([s != s for s in lower_limit]) == False:
                    out.write("pseudoatom pt1, pos=[" + str(upper_limit[0]) + "," + str(upper_limit[1]) + "," + str(
                        upper_limit[2]) + "]\n")
                    out.write("pseudoatom pt2, pos=[" + str(lower_limit[0]) + "," + str(lower_limit[1]) + "," + str(
                        lower_limit[2]) + "]\n")
                    out.write("distance axis, /pt1, /pt2\n")
                    out.write("set dash_gap, 0\n")
                    out.write("hide label, axis\n")
                    out.write("set dash_radius, 0.400\n")
                    out.write("color black, axis\n")
                out.write("rotate [1,0,0], angle=-90\n")
                # Find the correct residues to select for the repeat that is connected with each particular symmetry axis;
                # Note that we are not selecting each CE-Symm defined repeat but rather the repeats that correspond to each axis
                # So in (2A65.A_9-235;2A65.B_9-235)(2A65.A_242-441;2A65.B_242-441)
                # Repeat 1 is 2A65.A_9-235 and 2A65.A_242-441
                # Repeat 2 is 2A65.B_9-235 and 2A65.B_242-441
                k = 0
                rep_num = 0
                print(repeat_selection)
                for j in range(0, len(repeats[0])):
                    # sel_text=""
                    for s in range(0, len(repeats)):
                        sel_text = ""
                        entry = repeats[s][j]
                        try:
                            ch = entry.split('.')[1][0]
                        except:
                            ch = entry.split('_')[0]
                        beg = entry.split('_')[1]
                        first_digit = beg[0]  # in case the residues has a negative id (i.e. -1)
                        beg = beg[1:].split('-')[0]
                        beg = first_digit + beg
                        locator = [i for i, sublist in enumerate(repeat_selection) if ((int(beg), ch) == sublist[0])]
                        # print(locator)
                        for i in range(0, int(len(repeat_selection[locator[0]]) / 2)):
                            sel_text = sel_text + "(chain " + repeat_selection[locator[0]][i * 2][
                                1] + " and resid " + str(repeat_selection[locator[0]][i * 2][0]) + "-" + str(
                                repeat_selection[locator[0]][i * 2 + 1][0]) + ")"
                            if i < len(repeat_selection[locator[0]]) / 2 - 1:
                                sel_text = sel_text + " or "
                        # if s<(len(repeats)-1):
                        #    sel_text=sel_text+" or "
                        out.write("create repeat" + str(rep_num + 1) + ", " + inputf + " and " + sel_text + "\n")
                        out.write("color " + colors[(k % len(colors))] + ",repeat" + str(rep_num + 1) + "\n")
                        rep_num += 1
                    k += 1
                out.write("deselect\n")
                out.write("zoom " + inputf + "\n")
                out.write("cd " + pic_dir + "\n")
                out.write("png " + image_key + "_" + str(sym_num) + ".png, ray=" + str(ray) + "\n")
                out.write("quit")
                out.close()
                os.system(pymol + " " + wkdir + "pymol_script_" + image_key + ".pml")
                os.remove(wkdir + "pymol_script_" + image_key + ".pml")
                # Create gif file
                os.chdir(pic_dir)
                image_files = image_files + image_key + '_' + str(sym_num) + '.png;'

                # Create user-friendly pymol script
                out = open(pic_dir + "pymol_script_" + image_key + "_" + str(sym_num) + ".pml", "w")
                pml_files = pml_files + "pymol_script_" + image_key + "_" + str(sym_num) + ".pml;"
                if '_' not in inputf[5:7]:
                    input_name = inputf[0:4]
                else:
                    input_name = inputf[0:6]
                out.write("load " + input_name + ".pdb, " + inputf[0:4] + "\n")
                out.write("set antialias,1\nset line_smooth, 1\nset depth_cue, 0.5\nset specular, 0\n")
                out.write("set cartoon_sampling, 14\nset ribbon_sampling, 10\nset ray_trace_fog, 0.5\n")
                out.write("set ignore_case, off\n")
                out.write("hide everything,*\n")
                for c in cb_colors:
                    out.write(
                        "set_color cb_" + c[0] + "=[" + str(c[1][0]) + "," + str(c[1][1]) + "," + str(c[1][2]) + "]\n")
                out.write("show cartoon\n")
                out.write("color grey70, " + inputf[0:4] + "\n")
                out.write("set cartoon_transparency, 0.5, " + inputf[0:4] + "\n")
                out.write("bg_color white\nset opaque_background, off\n")
                if np.any([s != s for s in upper_limit]) == False or np.any([s != s for s in lower_limit]) == False:
                    out.write("pseudoatom pt1, pos=[" + str(upper_limit[0]) + "," + str(upper_limit[1]) + "," + str(
                        upper_limit[2]) + "]\n")
                    out.write("pseudoatom pt2, pos=[" + str(lower_limit[0]) + "," + str(lower_limit[1]) + "," + str(
                        lower_limit[2]) + "]\n")
                    out.write("distance axis, /pt1, /pt2\n")
                    out.write("set dash_gap, 0\n")
                    out.write("hide label, axis\n")
                    out.write("set dash_radius, 0.4\n")
                    out.write("color black, axis\n")
                else:
                    print("NaN values in axis\n")
                out.write("rotate [1,0,0], angle=-90\n")
                k = 0
                rep_num = 0
                for j in range(0, len(repeats[0])):
                    # sel_text=""
                    for s in range(0, len(repeats)):
                        sel_text = ""
                        entry = repeats[s][j]
                        try:
                            ch = entry.split('.')[1][0]
                        except:
                            ch = entry.split('_')[0]
                        beg = entry.split('_')[1]
                        first_digit = beg[0]
                        beg = beg[1:].split('-')[0]
                        beg = first_digit + beg
                        locator = [i for i, sublist in enumerate(repeat_selection) if ((int(beg), ch) == sublist[0])]
                        for i in range(0, int(len(repeat_selection[locator[0]]) / 2)):
                            sel_text = sel_text + "(chain " + repeat_selection[locator[0]][i * 2][
                                1] + " and resid " + str(repeat_selection[locator[0]][i * 2][0]) + "-" + str(
                                repeat_selection[locator[0]][i * 2 + 1][0]) + ")"
                            if i < len(repeat_selection[locator[0]]) / 2 - 1:
                                sel_text = sel_text + " or "
                        out.write("create repeat" + str(rep_num + 1) + ", " + inputf[0:4] + " and " + sel_text + "\n")
                        out.write("color " + colors[(k % len(colors))] + ",repeat" + str(rep_num + 1) + "\n")
                        rep_num += 1
                    k += 1
                out.write("deselect\n")
                out.write("zoom " + inputf[0:4] + "\n")
                out.close()
                os.chdir(wkdir)

                # Create 3dmol script
                out = open(jmol_dir + "3dmol_" + image_key + "_" + str(sym_num) + ".json", "w")
                jmol_files = jmol_files + "3dmol_" + image_key + "_" + str(sym_num) + ".json;"
                colors_hex = [c[2] for c in cb_colors] + [c[2] for c in pml_colors]
                out.write(
                    '{\n\t"rotate":[{"axis":"x", "angle":-90},{"axis":"y", "angle":0},{"axis":"z", "angle":0}],\n\t"chains": [\n')

                k = 0
                rep_num = 0
                rep_colors = []  # list where each cell corresponds to a different color  [[('A','3-20'),('B','1-15')],[('A','35-51')]]
                for j in range(0, len(repeats[0])):
                    # sel_text=""
                    rep_colors.append([])
                    for s in range(0, len(repeats)):
                        sel_text = ""
                        entry = repeats[s][j]
                        try:
                            ch = entry.split('.')[1][0]
                        except:
                            ch = entry.split('_')[0]
                        beg = entry.split('_')[1]
                        first_digit = beg[0]
                        beg = beg[1:].split('-')[0]
                        beg = first_digit + beg
                        locator = [i for i, sublist in enumerate(repeat_selection) if ((int(beg), ch) == sublist[0])]
                        for i in range(0, int(len(repeat_selection[locator[0]]) / 2)):
                            rep_colors[-1].append((repeat_selection[locator[0]][i * 2][1],
                                                   str(repeat_selection[locator[0]][i * 2][0]) + "-" + str(
                                                       repeat_selection[locator[0]][i * 2 + 1][0])))
                        rep_num += 1
                    k += 1
                chain_ordered = []  # create a list of the form [['A',('3-20','#ff3301'),('35-51','#d00322')],['B',('1-15','#ff3301')]]
                for j, r in enumerate(rep_colors):
                    for sel in r:
                        if len([chain_ordered[k] for k in range(0, len(chain_ordered)) if
                                sel[0] in chain_ordered[k]]) == 0:  # if the chain is not yet in chain_ordered
                            chain_ordered.append([sel[0], (sel[1], colors_hex[(j % len(colors_hex))])])
                        else:
                            ind = [k for k in range(0, len(chain_ordered)) if sel[0] in chain_ordered[k]]
                            assert len(ind) == 1
                            chain_ordered[ind[0]].append((sel[1], colors_hex[(j % len(colors_hex))]))
                for j, r in enumerate(chain_ordered):
                    out.write('\t\t{"chain":"' + r[0] + '",\n')
                    out.write('\t\t\t"color":"",\n')
                    out.write('\t\t\t"resid":[')
                    for k in range(1, len(r)):
                        out.write('{"pos":"' + r[k][0] + '", "color":"' + r[k][1] + '"}')
                        if k < len(r) - 1:
                            out.write(',')
                    out.write(']\n\t\t\t}')
                    if j < len(chain_ordered) - 1:
                        out.write(',\n')
                    else:
                        out.write('\n')

                out.write('\t\t],\n')
                if np.any([s != s for s in upper_limit]) == False or np.any([s != s for s in lower_limit]) == False:
                    out.write(
                        '\t"axis": {"pt1": {"x":' + str(upper_limit[0]) + ',"y":' + str(upper_limit[1]) + ',"z":' + str(
                            upper_limit[2]) + '}, ')
                    out.write('"pt2": {"x":' + str(lower_limit[0]) + ',"y":' + str(lower_limit[1]) + ',"z":' + str(
                        lower_limit[2]) + '},')
                    out.write('"axisradius": 0.4}\n')

                else:
                    out.write('\t"axis":""\n')
                out.write('}')
                out.close()

    return image_files, pml_files, jmol_files


def get_repeat_resid(cesymm_out_dir, inputf, reference):
    """
    CE-Symm lists repeats in a way that makes it unclear where the repeat ends in multi-chain structures. 
    For example, 2A65.A_9-235 indicates that the repeat starts with residue 9 in chain A but the end might 
    be residue 235 of chain C, encompassing chain B on the way. Therefore, we need to double-check the 
    repeat listed in the fasta alignments with the sequence from the pdb file and determine what ranges of 
    residues constitute each repeat. Here, algn_map is an array of the type 
    [[(9,'A'),(235,'A'),(1,'B'),(235,'B'),(1,'C'),(235,'C')],[...]], so that within each [] we have the ranges 
    for one particular repeat. 
     
    """
    resseq_A = get_pdb_sequence_with_chains(
        reference)  # creates a list of the form [(2, 'D', 'A'), (3, 'Q', 'A'), (4, 'A', 'A'), (5, 'L', 'A'), (6, 'S', 'A')...]
    cesm = list(SeqIO.parse(cesymm_out_dir + inputf + ".fasta", "fasta"))
    repeat_start = []
    algn_map = [[] for i in range(len(cesm))]
    for num, r in enumerate(cesm):
        # print(r.seq)
        r_id = r.id[len(inputf[0:4]):]
        r_id = r_id.split('_')
        ch = r_id[0][-1]
        # print(ch)
        repeat = r_id[1].strip()
        first_digit = repeat[0]  # handling negative residues
        repeat = repeat[1:].split('-')
        repeat[0] = first_digit + repeat[0]

        #    	res_start.append(repeat[0]) #this is the resid of the residue that starts the repeat
        repeat_start = [i for i, sublist in enumerate(resseq_A) if (int(repeat[0]) == sublist[0] and ch == sublist[
            2])]  # finds the index of the residue that starts the repeat
        #        print(repeat_start)
        count = 0
        small_lett = 0
        for i, c in enumerate(r.seq):
            if c != '-' and c != '/':
                if c.isupper() and resseq_A[repeat_start[0] + count][1] == c and count == 0:
                    algn_map[num].append((resseq_A[repeat_start[0] + count][0], resseq_A[repeat_start[0] + count][2]))
                #              elif resseq_A[repeat_start[0]+count][2]!=resseq_A[repeat_start[0]+count-1][2] and c.isupper() and resseq_A[repeat_start[0]+count][1] == c:
                #                algn_map[num].append((resseq_A[repeat_start[0]+count-1][0],resseq_A[repeat_start[0]+count-1][2]))
                #                algn_map[num].append((resseq_A[repeat_start[0]+count][0],resseq_A[repeat_start[0]+count][2]))
                elif resseq_A[repeat_start[0] + count][2] != resseq_A[repeat_start[0] + count - small_lett - 1][
                    2] and c.isupper() and resseq_A[repeat_start[0] + count][1] == c:
                    algn_map[num].append((resseq_A[repeat_start[0] + count - small_lett - 1][0],
                                          resseq_A[repeat_start[0] + count - small_lett - 1][2]))
                    algn_map[num].append((resseq_A[repeat_start[0] + count][0], resseq_A[repeat_start[0] + count][2]))
                elif len(r.seq) == i + 1 or r.seq[
                    i + 1] == '/':  # note that small letters at the endo of repeats are not counted by CE-Symm in the repeat id (it's a bug!) - they should be
                    algn_map[num].append((resseq_A[repeat_start[0] + count][0], resseq_A[repeat_start[0] + count][2]))
                    break
                count += 1
                if c.isupper():
                    small_lett = 0
                else:
                    small_lett += 1
            if c == '-' and (len(r.seq) == i + 1 or r.seq[i + 1] == '/') and len(
                    algn_map[num]) % 2 == 1:  # handles cases like ...KAF--
                algn_map[num].append(
                    (resseq_A[repeat_start[0] + count - 1][0], resseq_A[repeat_start[0] + count - 1][2]))
        if len(algn_map[num]) % 2 == 1:
            algn_map[num].append(algn_map[num][-1])
    return algn_map


def get_symd_aligned_resid(symd_out_dir, inputf, reference):
    resseq_A = get_pdb_sequence_with_chains(
        reference)  # creates a list of the form [(2, 'D', 'A'), (3, 'Q', 'A'), (4, 'A', 'A'), (5, 'L', 'A'), (6, 'S', 'A')...]
    symd_fasta = SeqIO.to_dict(SeqIO.parse(symd_out_dir + inputf + "-best.fasta", "fasta"))
    i = 0
    aligned = []
    aligned_chains = []
    end = '0'
    flag = 0
    for c in symd_fasta[inputf + '-original'].seq:
        # print(c, end="")
        if c != '-':
            if c.isupper():
                if resseq_A[i][2] not in aligned_chains:
                    # if end!='0':
                    #    print(ind, str(resseq_A[i-1][0]), c)
                    if end != '0' and last_char.isupper():
                        aligned[ind].append(str(resseq_A[i - 1][0]))
                    # if end!='0' and len(aligned[ind])%2==1: # provisions for a single capital letter in the beginning of a chain sequence
                    #    print('here')
                    #    aligned[ind].append(aligned[ind][-1])
                    aligned_chains.append(resseq_A[i][2])
                    if i == (
                            len(resseq_A) - 1):  # this provisions for a single capital letter at the very end of the sequence
                        aligned.append([str(resseq_A[i][0]), str(resseq_A[i][0])])
                    else:
                        aligned.append([str(resseq_A[i][0])])
                        ind = aligned_chains.index(resseq_A[i][2])
                        end = str(resseq_A[i][0])
                    flag = 1
                else:
                    ind = aligned_chains.index(resseq_A[i][2])
                    end = str(resseq_A[i][0])
                    flag = 1
                    if i == (len(resseq_A) - 1):
                        aligned[ind].append(end)
                    if last_char.islower():  # accounts for unaligned residues (breaks) in a single chain
                        aligned[ind].append(end)
                        flag = 1
            i += 1
            last_char = c
        if c.islower() and end != '0' and flag == 1:
            aligned[ind].append(end)
            flag = 0
    for i, x in enumerate(aligned):  # handles a single capital letter in the beginning of a chain sequence, eg. 3tdo
        if len(x) % 2 == 1:
            aligned[i].append(x[-1])
    print(aligned, aligned_chains)
    return aligned, aligned_chains  # a list of the residue ids and a list of the chains of all the capital letter residues in the alignment
    #                                eg: aligned=[[2,10,20,33],[1,1,5,18]] aligned_chains=['A','L']


def symd_images(wkdir, symd_out_dir, inputf, oriented, ray, pymol, pic_dir, jmol_dir, aligned, aligned_chains):
    f = open(symd_out_dir + inputf + '-trfm.pdb', 'r')
    print('aligned: ', aligned_chains, aligned)
    counter = 0
    flag = 0
    pts = []
    for line in f:
        if line.startswith("MODEL        3"):
            flag = 1
        if flag == 1 and 0 <= counter < 3 and line.startswith("ATOM") and line[13:16].strip() == 'CA':
            pts.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            counter += 1
    f.close()
    out = open(wkdir + "pymol_script_" + inputf + "_symd.pml", "w")
    out.write("load " + oriented + ", " + inputf + "\n")  # the original pdb is displayed, not the model one
    out.write('dss\n')
    out.write("set antialias, 1\nset line_smooth, 1.00\nset depth_cue, 0.5\nset specular, 0\nset surface_quality, 1\n")
    out.write("set cartoon_sampling, 14\nset ribbon_sampling, 10\nset ray_trace_fog, 0.5\n")
    out.write("viewport 350,350\n")
    out.write("set ignore_case, off\n")
    out.write("hide everything,*\n")
    out.write("show cartoon, " + inputf + "\n")
    out.write("color grey70, " + inputf + "\n")
    out.write("set cartoon_transparency, 0.5, " + inputf + "\n")
    out.write("set all_states, on\n")
    if len(aligned_chains) != 0:
        out.write("pseudoatom pt1, pos=[" + str(pts[0][0]) + "," + str(pts[0][1]) + "," + str(pts[0][2]) + "]\n")
        out.write("pseudoatom pt2, pos=[" + str(pts[2][0]) + "," + str(pts[2][1]) + "," + str(pts[2][2]) + "]\n")
        out.write("distance axis, /pt1, /pt2\n")
        out.write("set dash_gap, 0\n")
        out.write("hide label, axis\n")
        out.write("set dash_radius, 0.400\n")
        out.write("color black, axis\n")
    out.write("bg_color white\nset opaque_background, on\n")
    out.write("rotate [1,0,0], angle=-90\n")
    if len(aligned_chains) != 0:
        sel_text = "("
        for i, ch in enumerate(aligned_chains):
            sel_text = sel_text + "(chain " + ch + " and (resid "
            for k in range(0, len(aligned[i]) // 2):
                sel_text = sel_text + str(aligned[i][k * 2]) + "-" + str(aligned[i][k * 2 + 1])
                if k != len(aligned[i]) // 2 - 1:
                    sel_text = sel_text + " or resid "
            sel_text = sel_text + "))"
            if i < (len(aligned_chains) - 1):
                sel_text = sel_text + " or "
        sel_text = sel_text + ")"
        out.write("create aligned, " + inputf + " and " + sel_text + "\n")
        out.write("show cartoon, aligned\n")
        out.write("color marine, aligned\n")
    out.write("deselect\n")
    out.write("zoom " + inputf + "\n")
    out.write("cd " + pic_dir + "\n")
    out.write("png " + inputf + "_symd.png, ray=" + str(ray) + "\n")
    out.write("quit")
    out.close()
    os.system(pymol + " " + wkdir + "pymol_script_" + inputf + "_symd.pml")
    os.remove(wkdir + "pymol_script_" + inputf + "_symd.pml")
    os.chdir(pic_dir)

    out = open(pic_dir + "pymol_script_" + inputf + "_symd.pml", "w")
    out.write("load " + inputf + ".pdb, " + inputf + "\n")
    out.write('dss\n')
    out.write("set antialias, 1\nset line_smooth, 1.00\nset depth_cue, 0.5\nset specular, 0\nset surface_quality, 1\n")
    out.write("set cartoon_sampling, 14\nset ribbon_sampling, 10\nset ray_trace_fog, 0.5\n")
    out.write("viewport 350,350\n")
    out.write("set ignore_case, off\n")
    out.write("hide everything,*\n")
    out.write("show cartoon, " + inputf + "\n")
    out.write("color grey70, " + inputf + "\n")
    out.write("set cartoon_transparency, 0.5, " + inputf + "\n")
    out.write("set all_states, on\n")
    if len(aligned_chains) != 0:
        out.write("pseudoatom pt1, pos=[" + str(pts[0][0]) + "," + str(pts[0][1]) + "," + str(pts[0][2]) + "]\n")
        out.write("pseudoatom pt2, pos=[" + str(pts[2][0]) + "," + str(pts[2][1]) + "," + str(pts[2][2]) + "]\n")
        out.write("distance axis, /pt1, /pt2\n")
        out.write("set dash_gap, 0\n")
        out.write("hide label, axis\n")
        out.write("set dash_radius, 0.400\n")
        out.write("color black, axis\n")
    out.write("bg_color white\nset opaque_background, on\n")
    out.write("rotate [1,0,0], angle=-90\n")
    if len(aligned_chains) != 0:
        sel_text = "("
        for i, ch in enumerate(aligned_chains):
            sel_text = sel_text + "(chain " + ch + " and (resid "
            for k in range(0, len(aligned[i]) // 2):
                sel_text = sel_text + str(aligned[i][k * 2]) + "-" + str(aligned[i][k * 2 + 1])
                if k != len(aligned[i]) // 2 - 1:
                    sel_text = sel_text + " or resid "
            sel_text = sel_text + "))"
            if i < (len(aligned_chains) - 1):
                sel_text = sel_text + " or "
        sel_text = sel_text + ")"
        out.write("create aligned, " + inputf + " and " + sel_text + "\n")
        out.write("show cartoon, aligned\n")
        out.write("color marine, aligned\n")
    out.write("deselect\n")
    out.write("zoom " + inputf + "\n")
    out.close()
    os.chdir(wkdir)

    # Create 3dmol script
    out = open(jmol_dir + "3dmol_" + inputf + "_symd" + ".json", "w")
    out.write(
        '{\n\t"rotate":[{"axis":"x", "angle":-90},{"axis":"y", "angle":0},{"axis":"z", "angle":0}],\n\t"chains": [\n')
    if len(aligned_chains) != 0:
        for i, ch in enumerate(aligned_chains):
            out.write('\t\t{"chain":"' + ch + '",\n')
            out.write('\t\t\t"color":"",\n')
            out.write('\t\t\t"resid":[')
            for k in range(0, len(aligned[i]) // 2):
                out.write(
                    '{"pos":"' + str(aligned[i][k * 2]) + '-' + str(aligned[i][k * 2 + 1]) + '", "color":"#0080ff"}')
                if k < (len(aligned[i]) // 2 - 1):
                    out.write(',')
            out.write(']\n\t\t\t}')
            if i < len(aligned_chains) - 1:
                out.write(',\n')
            else:
                out.write('\n')
        out.write('\t\t],\n')

        out.write(
            '\t"axis": {"pt1": {"x":' + str(pts[0][0]) + ',"y":' + str(pts[0][1]) + ',"z":' + str(pts[0][2]) + '}, ')
        out.write('"pt2": {"x":' + str(pts[2][0]) + ',"y":' + str(pts[2][1]) + ',"z":' + str(pts[2][2]) + '},')
        out.write('"axisradius": 0.4}\n')
    else:
        out.write('\t"axis":""\n')
    out.write('}')
    out.close()
    return


#### Symmetry Transfer
def tm_score_calc(mobi_ca, transformed_ca, lmin):
    """Takes the lists of vector coordinates of the initial and transformed structures and calculates a TM score"""
    d0 = 1.24 * (lmin - 15) ** (1 / 3.0) - 1.8
    tm_score_elem = 0

    for i in range(0, len(mobi_ca)):
        tm_score_elem += 1.0 / (1 + pow(np.linalg.norm(mobi_ca[i] - transformed_ca[i]) / d0, 2))
    #        print(np.linalg.norm(mobi_ca[i]-transformed_ca[i]))
    return tm_score_elem / lmin


def rmsd_score(original_ca, transformed_ca):
    """Calculates the contribution to the RMSD from a given set of coordinates"""
    rmsd_elem = 0
    for i in range(0, len(transformed_ca)):
        rmsd_elem = rmsd_elem + pow(np.linalg.norm(transformed_ca[i] - original_ca[i]), 2)
    rmsd_elem = rmsd_elem / len(transformed_ca)
    return np.sqrt(rmsd_elem)


def get_repeats_alignment(wkdir, oriented, repeats,
                          muscle_path="muscle"):
    """

    :param wkdir: directory where temporary files can be saved
    :param oriented: path to pdb file
    :param repeats: e.g. [[A,B,C],[E,F,G],[I,J,K]]
    :param muscle_path: MUSCLE executable
    :return:
    """
    from random import randint
    id = str(randint(50, 2250))
    inf = open(wkdir + "tmp-" + id + ".fasta", "w")
    resseq = []

    for i, rep in enumerate(repeats):
        num = strip_tm_chains_in_order(wkdir, 'rep' + str(i) + "-" + id, oriented, rep)
        struct_tmp = parse_structure(wkdir + 'rep' + str(i) + "-" + id + '_tmp.pdb')[0]
        resseq_tmp = get_pdb_sequence_with_chains(struct_tmp)
        resseq.append(resseq_tmp)
        sequence_tmp = ''.join([i[1] for i in resseq_tmp])
        inf.write(">Seq" + str(i) + "\n")
        inf.write(sequence_tmp + "\n")
    inf.close()
    fnull = open(os.devnull, 'w')
    p = subprocess.Popen(
        muscle_path.split() + ['-in', wkdir + "tmp-" + id + ".fasta", '-out', wkdir + "out-" + id + ".fasta"],
        stdout=fnull, stderr=fnull)
    # p=subprocess.Popen(muscle_path.split() + ['-align', wkdir+"tmp.fasta", '-output', wkdir+"out.fasta"], stdout=fnull, stderr=fnull)
    p.wait()
    fnull.close()
    alns = sorted(SeqIO.parse(wkdir + "out-" + id + ".fasta", "fasta"), key=lambda x: int(x.id[3:]))
    os.remove(wkdir + "out-" + id + ".fasta")
    os.remove(wkdir + "tmp-" + id + ".fasta")

    aln_map = [[] for i in range(len(alns[0]))]
    for num, r in enumerate(alns):
        count = 0  # sequence count
        map_count = 0  # array count
        for c in r.seq:
            ind = int(r.id[len("Seq"):])
            if c == '-':
                aln_map[map_count].append(('-', '-', '-', '-'))
            else:
                assert resseq[ind][count][1] == c
                aln_map[map_count].append(
                    (resseq[ind][count][0], resseq[ind][count][2], resseq[ind][count][3], resseq[ind][count][4]))
                count += 1
            map_count += 1

    for i in range(0, len(repeats)):
        os.remove(wkdir + 'rep' + str(i) + '-' + id + '_tmp.pdb')

    return aln_map


def partial_repeats_scores(cesymm_out_dir, oriented, inputf):
    s_template = parse_structure(oriented)
    s_mobile = parse_structure(oriented)
    template = s_template[0]
    mobile = s_mobile[0]

    resseq_A = get_pdb_sequence_with_chains(
        template)  # creates a list of the form [(2, 'D', 'A'), (3, 'Q', 'A'), (4, 'A', 'A'), (5, 'L', 'A'), (6, 'S', 'A')...]
    cesm = SeqIO.parse(cesymm_out_dir + inputf + ".fasta", "fasta")
    unique_repeats = []
    for i in list(cesm):
        unique_repeats.append(i.id)
    cesm = SeqIO.to_dict(SeqIO.parse(cesymm_out_dir + inputf + ".fasta", "fasta"))
    repeats_all, repeats_type, repeats_level, axes_per_level, order_of_level = cesymm_related_repeats(inputf,
                                                                                                      cesymm_out_dir)
    print("partial_repeats:", cesymm_out_dir, repeats_all)
    # repeats_all=[[['1PV6.A_6-94', '1PV6.A_219-309'],['1PV6.A_101-187', '1PV6.A_310-399']],[['1PV6.A_6-94', '1PV6.A_101-187'],['1PV6.A_219-309', '1PV6.A_310-399']]]
    # Example: repeats_all=[[['1PV6.A_6-94', '1PV6.A_219-309'],['1PV6.A_101-187', '1PV6.A_310-399']],[['1PV6.A_6-94', '1PV6.A_101-187'],['1PV6.A_219-309', '1PV6.A_310-399']]]
    for ind, symm in enumerate(repeats_all):
        positions = sum([sum(1 for c in cesm[symm[k][1]].seq if c != '/') for k in range(0, len(symm))])
        algn_map = [[] for i in range(positions)]
        map_ind = 0
        for sub_symm in symm:  # sub_symm=['1PV6.A_6-94', '1PV6.A_219-309']
            for k in sub_symm:  # k='1PV6.A_6-94'
                r = cesm[k]
                r_id = k[len(inputf[0:4]):].split('_')
                ch = r_id[0][-1]
                repeat = r_id[1].strip()
                first_digit = repeat[0]  # handling negative residues
                repeat = repeat[1:].split('-')
                repeat[0] = first_digit + repeat[0]
                repeat_start = [i for i, sublist in enumerate(resseq_A) if (int(repeat[0]) in sublist and ch == sublist[
                    2])]  # finds the index of the residue that starts the repeat
                count = 0
                for i, c in enumerate(r.seq):  # Last edited
                    if c != '/':
                        if c != '-' and c.isupper() and resseq_A[repeat_start[0] + count][1] == c:
                            algn_map[i + map_ind].append((resseq_A[repeat_start[0] + count][0],
                                                          resseq_A[repeat_start[0] + count][2],
                                                          resseq_A[repeat_start[0] + count][3],
                                                          resseq_A[repeat_start[0] + count][4]))
                            count += 1
                        if c == '-':
                            algn_map[i + map_ind].append(('-', '-', '-', '-'))
                        if c.islower():
                            algn_map[i + map_ind].append(('-', '-', '-', '-'))
                            count += 1
            map_ind += sum(1 for c in cesm[k].seq if c != '/')

        rotations = 2
        tm_score = 0
        rmsd = 0
        for resi in algn_map:
            if all(elem == ('-', '-', '-', '-') for elem in resi) == False:
                rotations = len(resi)
                break

        # Create a list of atoms to be paired during superposition for obtaining the axis of rotation
        fixed_ca_list, mobile_ca_list = [], []
        fixed_rl = []
        mobile_rl = []
        fixed_coord = []

        if repeats_type[ind] == 'CLOSED':
            for k in range(0, rotations):
                for i, resi in enumerate(algn_map):
                    index = (rotations - 1 + k) % rotations
                    if resi[k] != ('-', '-', '-', '-') and resi[index] != ('-', '-', '-','-'):  # and s_mobile[0][resi[k][1]][(resi[k][2],resi[k][0],resi[k][3])]['CA'].get_altloc()==s_mobile[0][resi[index][1]][(resi[index][2],resi[index][0],resi[index][3])]['CA'].get_altloc():
                        fixed_ca_list.append(s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'])
                        fixed_coord.append(s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])][
                                               'CA'].get_vector().get_array())
                        mobile_ca_list.append(
                            s_mobile[0][resi[index][1]][(resi[index][2], resi[index][0], resi[index][3])]['CA'])
                        fixed_rl.append(s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])].get_id()[1])
                        mobile_rl.append(
                            s_mobile[0][resi[index][1]][(resi[index][2], resi[index][0], resi[index][3])].get_id()[1])
        elif repeats_type[ind] == 'OPEN':
            for k in range(1, rotations):
                for i, resi in enumerate(algn_map):
                    index = k - 1
                    if resi[k] != ('-', '-', '-', '-') and resi[index] != ('-', '-', '-', '-') and \
                            s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'].get_altloc() == \
                            s_mobile[0][resi[index][1]][(resi[index][2], resi[index][0], resi[index][3])][
                                'CA'].get_altloc():
                        fixed_ca_list.append(s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'])
                        fixed_coord.append(s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])][
                                               'CA'].get_vector().get_array())
                        mobile_ca_list.append(
                            s_mobile[0][resi[index][1]][(resi[index][2], resi[index][0], resi[index][3])]['CA'])
                        fixed_rl.append(s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])].get_id()[1])
                        mobile_rl.append(
                            s_mobile[0][resi[index][1]][(resi[index][2], resi[index][0], resi[index][3])].get_id()[1])

        # Superimpose matching residues
        si = Superimposer()
        si.set_atoms(fixed_ca_list, mobile_ca_list)
        T = np.array(np.transpose(si.rotran[0]))
        transl = si.rotran[1]
        ang, rotax, screw_transl, axis_pos = get_rotation_axis(si.rotran[0], transl)
        if len(axis_pos) == 2 and axis_pos == 'na':
            center = sum(fixed_coord) / len(fixed_coord)
            proj = np.dot(center, rotax)
            center_on_axis = np.array([proj * i for i in rotax])
            axis_pos = center - center_on_axis
        #    ax2=list(np.array(list(rotax))*10+transl)
        ax2 = axis_pos + 10 * np.array(rotax)
        axis_angle = round(angle_between(list(rotax), (0, 0, 1)) * 180 / np.pi, 2)
        if np.abs(np.dot(list(rotax), (0, 0, 1))) > 0.5:
            axis_type = "Parallel;"
        else:
            axis_type = "Antiparallel;"

        # Apply transformation to coordinates
        aligned_coord_res = []
        for i, resi in enumerate(algn_map):
            if np.any([resi[k] != ('-', '-', '-', '-') for k in range(0, len(resi))]):
                aligned_coord_res.append([])
                for k in range(0, len(resi)):
                    if resi[k] != ('-', '-', '-', '-'):
                        aligned_coord_res[-1].append(resi[k])
                        pt = list(s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'].get_vector())
                        for j in range(0, k):
                            pt = np.dot((pt - transl), T)
                        s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'].set_coord(pt)
                    else:
                        aligned_coord_res[-1].append(('-', '-', '-', '-'))

    #### Calculate the RMSD & TM-score ####
    unique_repeats = list(unique_repeats)
    positions = sum(1 for c in cesm[unique_repeats[0]].seq if c != '/')  # ind=1, positions=97
    algn_map = [[] for i in range(positions)]

    # Load the alignments of all repeats of the template against each other
    for k in unique_repeats:
        r = cesm[k]
        r_id = k[len(inputf[0:4]):].split('_')
        ch = r_id[0][-1]
        repeat = r_id[1].strip()
        first_digit = repeat[0]  # handling negative residues
        repeat = repeat[1:].split('-')
        repeat[0] = first_digit + repeat[0]
        repeat_start = [i for i, sublist in enumerate(resseq_A) if (int(repeat[0]) in sublist and ch == sublist[
            2])]  # finds the index of the residue that starts the repeat
        count = 0
        for i, c in enumerate(r.seq):  # Last edited
            if c != '/':
                if c != '-' and c.isupper() and resseq_A[repeat_start[0] + count][1] == c:
                    algn_map[i].append((resseq_A[repeat_start[0] + count][0], resseq_A[repeat_start[0] + count][2],
                                        resseq_A[repeat_start[0] + count][3], resseq_A[repeat_start[0] + count][4]))
                    count += 1

                if c == '-':
                    algn_map[i].append(('-', '-', '-', '-'))
                if c.islower():
                    algn_map[i].append(('-', '-', '-', '-'))
                    count += 1

    num_repeats = len(unique_repeats)
    # Find the length of the target protein
    prot_len = 0
    prot_len = len(resseq_A)

    #    rmsd=0
    tm_score = 0
    #    n=0
    dummy = list(range(0, num_repeats))  # a list of the N numbers representing the N unique repeats
    combo_list = list(combinations(dummy, 2))  # a combination of 2 elements from N
    for combo in combo_list:
        k = combo[0]
        ind = combo[1]
        repeat1 = []
        repeat2 = []
        for i, resi in enumerate(algn_map):
            # Check if the template alignment pair is valid and that both residues have a correspondence in the target
            if resi[k] != ('-', '-', '-', '-') and resi[ind] != ('-', '-', '-', '-') and \
                    s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'].get_altloc() == \
                    s_mobile[0][resi[ind][1]][(resi[ind][2], resi[ind][0], resi[ind][3])][
                        'CA'].get_altloc():  # checks that they are coordinates rather than ['-']
                #                rmsd=rmsd+pow(np.linalg.norm(s_mobile[0][resi[k][1]][(resi[k][2],resi[k][0],resi[k][3])]['CA'].get_vector()-s_mobile[0][resi[ind][1]][(resi[ind][2],resi[ind][0],resi[ind][3])]['CA'].get_vector()),2)
                repeat1.append(
                    s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'].get_vector().get_array())
                repeat2.append(s_mobile[0][resi[ind][1]][(resi[ind][2], resi[ind][0], resi[ind][3])][
                                   'CA'].get_vector().get_array())
            #                n+=1
        tm_score = tm_score + tm_score_calc(repeat1, repeat2, prot_len)

    tm_score = tm_score * num_repeats / len(combo_list)
    #    rmsd=np.sqrt(rmsd/n)
    angle = ang * 180.0 / np.pi

    # Store coordinates for RMSD procedure
    coord_repeats = []
    altlocs = []
    for i, resi in enumerate(algn_map):
        # Check if the template alignment pair is valid and that both residues have a correspondence in the target
        coord_repeats.append([])
        altlocs.append([])
        for k in range(0, len(resi)):
            if resi[k] != ('-', '-', '-', '-'):
                coord_repeats[i].append(
                    s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'].get_vector().get_array())
                altlocs[i].append(s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'].get_altloc())
            else:
                coord_repeats[i].append(None)
                altlocs[i].append(None)

    tmp = list(map(list, zip(*coord_repeats)))
    coord_repeats = tmp
    tmp = list(map(list, zip(*altlocs)))
    altlocs = tmp
    rmsd = get_rmsd_cesymm(coord_repeats, altlocs)  # Java RMSD

    return float('%.2f' % tm_score), float('%.2f' % rmsd), float('%.2f' % angle), prot_len


def cesymm_rmsd_tmscore(oriented, algn_map, repeats_type):
    """
    Apply the CE-Symm procedure for obtaining TM-score and RMSD given a single-level symmetry.
    Useful in case we have a list of repeats and want to see what the score is.
    """
    s_target = parse_structure(oriented)
    s_mobile = parse_structure(oriented)
    target = s_target[0]
    mobile = s_mobile[0]

    rotations = 2
    tm_score = 0
    rmsd = 0
    for resi in algn_map:
        if all(elem == ('-', '-', '-', '-') for elem in resi) == False:
            rotations = len(resi)
            break

    # Create a list of atoms to be paired during superposition for obtaining the axis of rotation
    fixed_ca_list, mobile_ca_list = [], []
    fixed_rl = []
    mobile_rl = []
    rot_matrix = []
    transl_vector = []
    fixed_coord = []

    if repeats_type == 'CLOSED':
        for k in range(0, rotations):
            for i, resi in enumerate(algn_map):
                index = (rotations - 1 + k) % rotations
                if resi[k] != ('-', '-', '-', '-') and resi[index] != ('-', '-', '-', '-') and \
                        s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'].get_altloc() == \
                        s_mobile[0][resi[index][1]][(resi[index][2], resi[index][0], resi[index][3])][
                            'CA'].get_altloc():
                    fixed_ca_list.append(s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'])
                    fixed_coord.append(
                        s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'].get_vector().get_array())
                    mobile_ca_list.append(
                        s_mobile[0][resi[index][1]][(resi[index][2], resi[index][0], resi[index][3])]['CA'])
                    fixed_rl.append(s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])].get_id()[1])
                    mobile_rl.append(
                        s_mobile[0][resi[index][1]][(resi[index][2], resi[index][0], resi[index][3])].get_id()[1])
    elif repeats_type == 'OPEN':
        for k in range(1, rotations):
            for i, resi in enumerate(algn_map):
                index = k - 1
                if resi[k] != ('-', '-', '-', '-') and resi[index] != ('-', '-', '-', '-') and \
                        s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'].get_altloc() == \
                        s_mobile[0][resi[index][1]][(resi[index][2], resi[index][0], resi[index][3])][
                            'CA'].get_altloc():
                    fixed_ca_list.append(s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'])
                    fixed_coord.append(
                        s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'].get_vector().get_array())
                    mobile_ca_list.append(
                        s_mobile[0][resi[index][1]][(resi[index][2], resi[index][0], resi[index][3])]['CA'])
                    fixed_rl.append(s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])].get_id()[1])
                    mobile_rl.append(
                        s_mobile[0][resi[index][1]][(resi[index][2], resi[index][0], resi[index][3])].get_id()[1])

    # print(fixed_rl)
    # print(mobile_rl)
    # Superimpose matching residues
    si = Superimposer()
    si.set_atoms(fixed_ca_list, mobile_ca_list)
    T = np.array(np.transpose(si.rotran[0]))
    rot_matrix.append(T)
    transl = si.rotran[1]
    transl_vector.append(transl)
    ang, rotax, screw_transl, axis_pos = get_rotation_axis(si.rotran[0], transl)
    if len(axis_pos) == 2 and axis_pos == 'na':
        center = sum(fixed_coord) / len(fixed_coord)
        proj = np.dot(center, rotax)
        center_on_axis = np.array([proj * i for i in rotax])
        axis_pos = center - center_on_axis
    #    ax2=list(np.array(list(rotax))*10+transl)
    ax2 = axis_pos + 10 * np.array(rotax)
    axis_angle = round(angle_between(list(rotax), (0, 0, 1)) * 180 / np.pi, 2)
    if axis_angle > 90:
        axis_angle = 90 - axis_angle
    if np.abs(np.dot(list(rotax), (0, 0, 1))) > 0.5:
        axis_type = "Parallel;"
    else:
        axis_type = "Antiparallel;"

    # Apply transformation to coordinates
    aligned_coord_res = []
    for i, resi in enumerate(algn_map):
        if np.any([resi[k] != ('-', '-', '-', '-') for k in range(0, len(resi))]):
            aligned_coord_res.append([])
            for k in range(0, len(resi)):
                if resi[k] != ('-', '-', '-', '-'):
                    aligned_coord_res[-1].append(resi[k])
                    pt = list(s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'].get_vector())
                    for j in range(0, k):
                        pt = np.dot((pt - transl), T)
                    s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'].set_coord(pt)
                else:
                    aligned_coord_res[-1].append(('-', '-', '-', '-'))

    #### Calculate the RMSD & TM-score ####
    num_repeats = len(algn_map[0])
    # Find the length of the target protein
    prot_len = 0
    resseq_B = get_pdb_sequence_with_chains(target)
    prot_len = len(resseq_B)

    # rmsd=0
    tm_score = 0
    # n=0
    dummy = list(range(0, num_repeats))  # a list of the N numbers representing the N unique repeats
    combo_list = list(combinations(dummy, 2))  # a combination of 2 elements from N
    for combo in combo_list:
        k = combo[0]
        ind = combo[1]
        repeat1 = []
        repeat2 = []
        for i, resi in enumerate(algn_map):
            # Check if the template alignment pair is valid and that both residues have a correspondence in the target
            if resi[k] != ('-', '-', '-', '-') and resi[ind] != ('-', '-', '-', '-') and \
                    s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'].get_altloc() == \
                    s_mobile[0][resi[ind][1]][(resi[ind][2], resi[ind][0], resi[ind][3])][
                        'CA'].get_altloc():  # checks that they are coordinates rather than ['-']
                # rmsd=rmsd+pow(np.linalg.norm(s_mobile[0][resi[k][1]][(resi[k][2],resi[k][0],resi[k][3])]['CA'].get_vector()-s_mobile[0][resi[ind][1]][(resi[ind][2],resi[ind][0],resi[ind][3])]['CA'].get_vector()),2)
                repeat1.append(
                    s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'].get_vector().get_array())
                repeat2.append(s_mobile[0][resi[ind][1]][(resi[ind][2], resi[ind][0], resi[ind][3])][
                                   'CA'].get_vector().get_array())
                # n+=1
        # print(repeat1, repeat2)
        tm_score = tm_score + tm_score_calc(repeat1, repeat2, prot_len)

    tm_score = tm_score * num_repeats / len(combo_list)
    # rmsd=np.sqrt(rmsd/n)
    angle = ang * 180.0 / np.pi
    # Store coordinates for RMSD procedure
    coord_repeats = []
    altlocs = []
    for i, resi in enumerate(algn_map):
        # Check if the template alignment pair is valid and that both residues have a correspondence in the target
        coord_repeats.append([])
        altlocs.append([])
        for k in range(0, len(resi)):
            if resi[k] != ('-', '-', '-', '-'):
                coord_repeats[i].append(
                    s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'].get_vector().get_array())
                altlocs[i].append(s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'].get_altloc())
            else:
                coord_repeats[i].append(None)
                altlocs[i].append(None)

    tmp = list(map(list, zip(*coord_repeats)))
    coord_repeats = tmp
    tmp = list(map(list, zip(*altlocs)))
    altlocs = tmp
    rmsd = get_rmsd_cesymm(coord_repeats, altlocs)  # Java RMSD

    #    transl_magnitude=np.linalg.norm(transl)

    return float('%.2f' % tm_score), float('%.2f' % rmsd), float('%.2f' % angle), float(
        '%.2f' % screw_transl), ax2, axis_pos, prot_len, axis_angle, axis_type, rot_matrix, transl_vector


def create_fasta_file(out_dir, inputf, oriented, aln_map, structure, image_key):
    output = open(out_dir + inputf + image_key + ".fasta", "w")
    structure = parse_structure(oriented)[0]
    repeats = "("
    residue_count = 0
    for i in range(0, len(aln_map[0])):
        seq = ""
        tag = ">" + inputf[0:4].upper() + "."
        repeats = repeats + inputf[0:4].upper() + "." # e.g. 3wmm
        end = -1
        for j in range(0, len(aln_map)):
            if end == -1 and aln_map[j][i][0] != '-':
                tag = tag + aln_map[j][i][1] + "_" + str(aln_map[j][i][0])
                repeats = repeats + aln_map[j][i][1] + "_" + str(aln_map[j][i][0])
            if aln_map[j][i][0] != '-':
                seq = seq + aa3to1.get(
                    structure[aln_map[j][i][1]][(aln_map[j][i][2], aln_map[j][i][0], aln_map[j][i][3])].resname, 'X')
                residue_count += 1
                end = j
            if aln_map[j][i][0] == '-':
                seq = seq + '-'
        tag = tag + "-" + str(aln_map[end][i][0])
        #        tag=tag+"-"+aln_map[end][i][1]+"_"+str(aln_map[end][i][0])         # if we want to include the name of the chain for the last residue, i.e. 2a65.A_1-B_5 rather than 2a65.A_1-5
        #        repeats=repeats+"-"+aln_map[end][i][1]+"_"+str(aln_map[end][i][0])+";"
        repeats = repeats + "-" + str(aln_map[end][i][0]) + ";"
        output.write(tag + "\n")
        output.write(seq + "\n")
    rep_length = len(seq)
    repeats = repeats[:-1] + ")"
    output.close()
    return out_dir + inputf + image_key + ".fasta", repeats, len(aln_map[0]), residue_count, rep_length


def create_axes_file(out_dir, inputf, angle, transl, pt1, pt2, repeats, num_reps, image_key, sym_type):
    output = open(out_dir + inputf + image_key + ".axes", "w")
    output.write(
        "Name    SymmLevel       SymmType        SymmOrder       RotationAngle   ScrewTranslation        Point1  Point2  AlignedRepeats\n")
    output.write("%s\t%d\t%s\t%d\t%.2f\t%.2f\t%.2f,%.2f,%.2f\t%.2f,%.2f,%.2f\t%s\n" % (
    inputf[0:4], 1, sym_type, num_reps, angle, transl, pt1[0], pt1[1], pt1[2], pt2[0], pt2[1], pt2[2], repeats))
    return out_dir + inputf + image_key + ".axes"


def bug_fixer(muscle_path, wkdir, inputf, oriented, repeats_list, files_key):
    bug_fix_dic = {}
    aln_map = get_repeats_alignment(wkdir, oriented, repeats_list, muscle_path)
    tmscore, rmsd, angle, translation, pt1, pt2, prot_len, axis_angle, axis_type, rot_matrix, transl_vector = cesymm_rmsd_tmscore(
        oriented, aln_map, 'CLOSED')
    print(tmscore, rmsd, angle)
    if tmscore > 0.90:
        try:
            fasta_path, repeats, num_repeats, aligned, rep_length = create_fasta_file(wkdir, inputf, oriented, aln_map,
                                                                                      oriented, files_key)
            axes_path = create_axes_file(wkdir, inputf, angle, translation, pt1, pt2, repeats, num_repeats, files_key,
                                         'CLOSED')
            reference = parse_structure(oriented)[0]
            repeat_selection = get_repeat_resid(wkdir, inputf + files_key,
                                                reference)  # the chain order matters for this one
            # images,pml=cesymm_images(wkdir,wkdir,inputf+files_key,repeat_selection,oriented,reference,ray,pymol,pic_dir,files_key)
            bug_fix_dic['coverage'] = aligned / prot_len
            bug_fix_dic['size'] = prot_len
            bug_fix_dic['unit_angle'] = angle
            bug_fix_dic['unit_translation'] = translation
            bug_fix_dic['refined_tmscore'] = tmscore
            bug_fix_dic['refined_rmsd'] = rmsd
            bug_fix_dic['repeats_number'] = num_repeats
            bug_fix_dic['topology'] = axis_type
            bug_fix_dic['symmetry_order'] = 'C' + str(num_repeats)
            bug_fix_dic['symmetry_levels'] = 1
            bug_fix_dic['aligned_length'] = aligned
            bug_fix_dic['repeat_length'] = rep_length
            bug_fix_dic['internal_symmetry'] = 'No'
            bug_fix_dic['quaternary_symmetry'] = 'Yes'
            bug_fix_dic['axis_angle_with_membrane_normal'] = str(float('%.2f' % axis_angle)) + ';'
            bug_fix_dic['symmetry_type'] = 'Quaternary;'
            bug_fix_dic['repeats'] = repeats.replace(";", ",") + ";"
            # bug_fix_dic['repeats'] = (repeats.replace(";", ",") + ";").replace(inputf[0:4].upper()+".", "")
            bug_fix_dic['files_key'] = inputf + files_key
            bug_fix_dic['repeat_selection'] = repeat_selection
            bug_fix_dic['closed_opened'] = 'CLOSED;'
            bug_fix_dic['repeat_level'] = '1;'
        except:
            print(sys.exc_info())
            pass
    else:
        tmscore, rmsd, angle, translation, pt1, pt2, prot_len, axis_angle, axis_type, rot_matrix, transl_vector = cesymm_rmsd_tmscore(
            oriented, aln_map, 'OPEN')
        print(tmscore, rmsd, angle)
        if tmscore > 0.98:
            try:
                fasta_path, repeats, num_repeats, aligned, rep_length = create_fasta_file(wkdir, inputf, oriented,
                                                                                          aln_map, oriented, files_key)
                axes_path = create_axes_file(wkdir, inputf, angle, translation, pt1, pt2, repeats, num_repeats,
                                             files_key, 'OPEN')
                reference = parse_structure(oriented)[0]
                repeat_selection = get_repeat_resid(wkdir, inputf + files_key,
                                                    reference)  # the chain order matters for this one
                # files_key="_analysis"
                # images,pml=cesymm_images(wkdir,wkdir,inputf+files_key,repeat_selection,oriented,reference,ray,pymol,pic_dir,files_key)
                bug_fix_dic['coverage'] = aligned / prot_len
                bug_fix_dic['size'] = prot_len
                bug_fix_dic['unit_angle'] = angle
                bug_fix_dic['unit_translation'] = translation
                bug_fix_dic['refined_tmscore'] = tmscore
                bug_fix_dic['refined_rmsd'] = rmsd
                bug_fix_dic['repeats_number'] = num_repeats
                bug_fix_dic['topology'] = axis_type
                bug_fix_dic['symmetry_order'] = 'R'
                bug_fix_dic['symmetry_levels'] = 1
                bug_fix_dic['aligned_length'] = aligned
                bug_fix_dic['repeat_length'] = rep_length
                bug_fix_dic['internal_symmetry'] = 'No'
                bug_fix_dic['quaternary_symmetry'] = 'Yes'
                bug_fix_dic['axis_angle_with_membrane_normal'] = str(float('%.2f' % axis_angle)) + ';'
                bug_fix_dic['symmetry_type'] = 'Quaternary;'
                bug_fix_dic['repeats'] = repeats.replace(";", ",") + ";"
                #bug_fix_dic['repeats'] = (repeats.replace(";", ",") + ";").replace(inputf[0:4].upper()+".", "")
                bug_fix_dic['files_key'] = files_key
                bug_fix_dic['repeat_selection'] = repeat_selection
                bug_fix_dic['closed_opened'] = 'OPEN;'
                bug_fix_dic['repeat_level'] = '1;'
            except:
                print(sys.exc_info())
                pass
    return bug_fix_dic


def rotation_axis_big_angle(rot, transl, c, theta):
    """
    This is only for cases where the angle is at least 5 degrees different from 0
    rot is the 3x3 rotation matrix
    transl is the translation vector
    c is the cosine of the rotation angle (cos(theta))
    
    We are following the method descirbed in Kim et al. 2010 (SI) and implemented in 
    CE-Symm's function calculateRotationalAxis
     
    """
    sum = 0
    rot_ax = []

    for i in range(0, 3):
        rot_ax.append(np.sqrt(rot[i, i] - c))  # v_i=np.sqrt(R[i,i]-c)
        sum += rot_ax[i] * rot_ax[i]
    for i in range(0, 3):
        rot_ax[i] = rot_ax[i] / np.sqrt(sum)

    # Determine the sign
    d0 = rot[2, 1] - rot[1, 2]  # =2u[0]*sin(theta]
    d1 = rot[0, 2] - rot[2, 0]  # =2u[1]*sin(theta]
    d2 = rot[1, 0] - rot[0, 1]  # =2u[2]*sin(theta]

    s12 = rot[2, 1] + rot[1, 2]  # =2*u[1]*u[2]*(1-cos(theta]]
    s02 = rot[0, 2] + rot[2, 0]  # =2*u[0]*u[2]*(1-cos(theta]]
    s01 = rot[1, 0] + rot[0, 1]  # =2*u[0]*u[1]*(1-cos(theta]]

    # Take the biggest d for the sign to ensure numerical stability
    if np.abs(d0) < np.abs(d1):  # not d0
        if np.abs(d1) < np.abs(d2):  # d2
            if d2 >= 0:  # u[2] positive
                if s02 < 0:
                    rot_ax[0] = -rot_ax[0]
                if s12 < 0:
                    rot_ax[1] = -rot_ax[1]
            else:  # u[2] negative
                rot_ax[2] = -rot_ax[2]
                if s02 >= 0:
                    rot_ax[0] = -rot_ax[0]
                if s12 >= 0:
                    rot_ax[1] = -rot_ax[1]

        else:  # d1
            if (d1 >= 0):  # u[1] positive
                if s01 < 0:
                    rot_ax[0] = -rot_ax[0]
                if s12 < 0:
                    rot_ax[2] = -rot_ax[2]
            else:  # u[1] negative
                rot_ax[1] = -rot_ax[1]
                if s01 >= 0:
                    rot_ax[0] = -rot_ax[0]
                if s12 >= 0:
                    rot_ax[2] = -rot_ax[2]
    else:  # not d1
        if (np.abs(d0) < np.abs(d2)):  # d2
            if (d2 >= 0):  # u[2] positive
                if s02 < 0:
                    rot_ax[0] = -rot_ax[0]
                if s12 < 0:
                    rot_ax[1] = -rot_ax[1]
            else:  # u[2] negative
                rot_ax[2] = -rot_ax[2]
                if s02 >= 0:
                    rot_ax[0] = -rot_ax[0]
                if s12 >= 0:
                    rot_ax[1] = -rot_ax[1]
        else:  # d0
            if (d0 >= 0):  # u[0] positive
                if s01 < 0:
                    rot_ax[1] = -rot_ax[1]
                if s02 < 0:
                    rot_ax[2] = -rot_ax[2]
            else:  # u[0] negative
                rot_ax[0] = -rot_ax[0]
                if s01 >= 0:
                    rot_ax[1] = -rot_ax[1]
                if s02 >= 0:
                    rot_ax[2] = -rot_ax[2]

    scale = np.dot(transl, np.array(rot_ax))
    screw_transl_vec = scale * np.array(rot_ax)
    perp_transl_vec = transl - screw_transl_vec
    h = np.cross(perp_transl_vec, np.array(rot_ax)) / (2.0 * np.tan(theta / 2.0))
    axis_pos = h + 0.5 * perp_transl_vec

    return rot_ax, np.abs(scale), axis_pos


def rotation_axis_small_angle(rot, transl):
    rot_ax = transl / np.linalg.norm(transl)
    scale = np.dot(transl, np.array(rot_ax))
    axis_pos = 'na'
    return rot_ax, np.abs(scale), axis_pos


def get_rotation_axis(rotation, transl):
    c = (np.trace(rotation) - 1) / 2.0  # =cos(theta)
    # c is sometimes slightly out of the [-1,1] range due to numerical instabilities
    if (-1 - 1e-8 < c and c < -1):
        c = -1
    if (1 + 1e-8 > c and c > 1):
        c = 1
    if (-1 > c or c > 1):
        raise SystemExit("Not a valid rotation matrix")

    theta = np.arccos(c)
    min_angle = 5 * np.pi / 180
    if theta < min_angle:
        rot_ax, screw_transl, axis_pos = rotation_axis_small_angle(rotation, transl)
    else:
        rot_ax, screw_transl, axis_pos = rotation_axis_big_angle(rotation, transl, c, theta)
    return theta, rot_ax, screw_transl, axis_pos


def get_rmsd_cesymm(coord_repeats, altlocs):
    ### CE-Symm RMSD
    sumSqDist = 0
    comparisons = 0

    for r1 in range(0, len(coord_repeats)):
        for c in range(0, len(coord_repeats[r1])):
            refAtom = coord_repeats[r1][c]
            refAltloc = altlocs[r1][c]
            if refAtom is None:
                continue

            nonNullSqDist = 0
            nonNullLength = 0
            for r2 in range(r1 + 1, len(coord_repeats)):
                atom = coord_repeats[r2][c]
                atom_altloc = altlocs[r2][c]
                if atom is not None and refAltloc == atom_altloc:
                    nonNullSqDist += (np.linalg.norm(refAtom - atom)) ** 2
                    nonNullLength += 1

            if nonNullLength > 0:
                comparisons += 1
                sumSqDist += 1.0 * nonNullSqDist / nonNullLength
    return np.sqrt(sumSqDist / comparisons)


def tm_score_calc_cs(mobi_ca, transformed_ca, l1, l2):
    lmin = min(l1, l2)
    """Takes the lists of vector coordinates of the initial and transformed structures and calculates a TM score"""
    d0 = 1.24 * (lmin - 15) ** (1 / 3.0) - 1.8
    tm_score_elem = 0

    for i in range(0, len(mobi_ca)):
        tm_score_elem += 1.0 / (1 + pow(np.linalg.norm(mobi_ca[i] - transformed_ca[i]) / d0, 2))
    #        print(np.linalg.norm(mobi_ca[i]-transformed_ca[i]))
    return tm_score_elem / lmin


def get_tmscore_cesymm(coord_repeats):
    sumTM = 0
    comparisons = 0
    for r1 in range(0, len(coord_repeats)):
        for r2 in range(r1 + 1, len(coord_repeats)):
            ln = len(coord_repeats[r1])
            ref = [None] * ln
            aln = [None] * ln
            nonNullLen = 0
            for c in range(0, ln):
                if (coord_repeats[r1][c] is None) == False and (coord_repeats[r2][c] is None) == False:
                    ref[nonNullLen] = coord_repeats[r1][c]
                    aln[nonNullLen] = coord_repeats[r2][c]
                    nonNullLen += 1
            if nonNullLen < ln:
                ref = ref[0:nonNullLen]
                aln = aln[0:nonNullLen]
            sumTM += tm_score_calc_cs(ref, aln, len(coord_repeats[r1]), len(coord_repeats[r2]))
            comparisons += 1

    return sumTM / comparisons


def transfer_symmetry(wkdir, frtm_path, cesymm_dir, transfer_out_dir, oriented_template, oriented_target, inputf,
                      target_pdb, ce_order, aln_file, super_pml_dir, image_key):

    """
    Take the repeats found in the template and find their relationships (i.e. RMSD/TM-score) in the target
    """

    chain = ""
    chain2 = ""
    # Load the template
    s_template = parse_structure(oriented_template)
    if chain == "":
        template = s_template[0]
    else:
        try:
            template = s_template[0][chain]
        except KeyError:
            raise Exception('Chain {0} not found in template structure'.format(chain))

    # Load the target and its movable copy mobile
    s_target = parse_structure(oriented_target)
    s_mobile = parse_structure(oriented_target)
    if chain2 == "":
        target = s_target[0]
        mobile = s_mobile[0]
    else:
        try:
            target = s_target[0][chain2]
            mobile = s_mobile[0][chain2]
        except KeyError:
            raise Exception('Chain {0} not found in target structure'.format(chain2))

    # Check if alignment exists
    if os.path.isfile(aln_file) == False:
        print("Alignment file does not exist.")
        aln_file = wkdir + "dummy.txt"
        f = open(aln_file, "w")
        f.write("no alignmnet")
        f.close()

    resseq_A = get_pdb_sequence_with_chains(
        template)  # creates a list of the form [(2, 'D', 'A'), (3, 'Q', 'A'), (4, 'A', 'A'), (5, 'L', 'A'), (6, 'S', 'A')...]
    cesm = SeqIO.parse(cesymm_dir + inputf + ".fasta", "fasta")
    unique_repeats = []
    for i in list(cesm):
        unique_repeats.append(i.id)
        parent_id = i.id[0:4].lower() + "_" + i.id[5:6]
    cesm = SeqIO.to_dict(SeqIO.parse(cesymm_dir + inputf + ".fasta", "fasta"))
    repeats_all, repeats_type, repeats_levels, axes_per_level, order_of_level = cesymm_related_repeats(inputf,
                                                                                                       cesymm_dir)

    # Align sequences to get mapping between residues
    aln, seq_id, gapless_id, res_map = frtmalign_sequences_with_chains(wkdir, template, target, aln_file,
                                                                       oriented_template, oriented_target, parent_id,
                                                                       target_pdb,
                                                                       frtmalign_path=frtm_path)
    # print(res_map)
    # Create a list of the residue ids of the target structure
    target_res_list = []
    for res in target.get_residues():
        if is_aa(res.get_resname()):
            target_res_list.append((res.get_id()[1], res.get_parent().get_id(), res.get_id()[0], res.get_id()[2]))
    # print(target_res_list)
    target_sym_descriptors = []
    transforms = []
    transl_vectors = []
    repeats_transf = [[-1 for j in range(0, len(unique_repeats))] for i in range(0,
                                                                                 len(repeats_all))]  # same number of elements as symmetry levels, eg. [[-1, -1, -1, -1], [-1, -1, -1, -1]]

    # template_ca_list, target_ca_list = [], []
    # for temp_res in res_map:
    #    template_ca_list.append(s_template[0][temp_res[1]][temp_res[0]]['CA'])
    #    target_ca_list.append(s_target[0][res_map[temp_res][1]][temp_res[0]]['CA'])

    # repeats_all=[[['1PV6.A_6-94', '1PV6.A_219-309'],['1PV6.A_101-187', '1PV6.A_310-399']],[['1PV6.A_6-94', '1PV6.A_101-187'],['1PV6.A_219-309', '1PV6.A_310-399']]]
    # Example: repeats_all=[[['1PV6.A_6-94', '1PV6.A_219-309'],['1PV6.A_101-187', '1PV6.A_310-399']],[['1PV6.A_6-94', '1PV6.A_101-187'],['1PV6.A_219-309', '1PV6.A_310-399']]]
    mode = ""
    for ind, symm in enumerate(repeats_all):
        positions = sum([sum(1 for c in cesm[symm[k][1]].seq if c != '/') for k in range(0, len(symm))])
        algn_map = [[] for i in range(positions)]
        map_ind = 0
        for sub_symm in symm:  # sub_symm=['1PV6.A_6-94', '1PV6.A_219-309']
            for l, k in enumerate(sub_symm):  # k='1PV6.A_6-94'
                r = cesm[k]
                r_id = k[len(inputf[0:4]):].split('_')
                ch = r_id[0][-1]
                repeat = r_id[1].strip()
                first_digit = repeat[0]  # handling negative residues
                repeat = repeat[1:].split('-')
                repeat[0] = first_digit + repeat[0]
                #    	res_start.append(repeat[0]) #this is the resid of the residue that starts the repeat
                repeat_start = [i for i, sublist in enumerate(resseq_A) if (int(repeat[0]) in sublist and ch == sublist[
                    2])]  # finds the index of the residue that starts the repeat
                #    print(repeat_start)
                count = 0
                for i, c in enumerate(r.seq):  # Last edited
                    if c != '/':
                        if c != '-' and c.isupper() and resseq_A[repeat_start[0] + count][1] != c:
                            # print(c, resseq_A[repeat_start[0]+count])
                            algn_map[i + map_ind].append(('-', '-', '-', '-'))
                        if c != '-' and c.isupper() and resseq_A[repeat_start[0] + count][1] == c:
                            #    		    print(c,resseq_A[repeat_start[0]+count][1],resseq_A[repeat_start[0]+count][0],"\n")
                            algn_map[i + map_ind].append((resseq_A[repeat_start[0] + count][0],
                                                          resseq_A[repeat_start[0] + count][2],
                                                          resseq_A[repeat_start[0] + count][3],
                                                          resseq_A[repeat_start[0] + count][4]))
                            count += 1
                        #           if resseq_A[repeat_start[0]+count][2] not in chains:
                        #               chains.append(resseq_A[repeat_start[0]+count][2])
                        if c == '-':
                            algn_map[i + map_ind].append(('-', '-', '-', '-'))
                        if c.islower():
                            algn_map[i + map_ind].append(('-', '-', '-', '-'))
                            count += 1
                u = unique_repeats.index(k)
                repeats_transf[ind][u] = l
            map_ind += sum(1 for c in cesm[k].seq if c != '/')

        # print(chains)

        # print(algn_map)
        rotations = 2
        tm_score = 0
        rmsd = 0
        for resi in algn_map:
            if all(elem == ('-', '-', '-', '-') for elem in resi) == False:
                rotations = len(resi)
                break

        # Create a list of atoms to be paired during superposition for obtaining the axis of rotation
        fixed_ca_list, mobile_ca_list = [], []
        fixed_rl = []
        mobile_rl = []
        fixed_coord = []

        if repeats_type[ind] == 'CLOSED':
            for k in range(0, rotations):
                for i, resi in enumerate(algn_map):
                    index = (rotations - 1 + k) % rotations
                    # print(k, index, resi[k], resi[index])
                    if resi[k] != ('-', '-', '-', '-') and resi[index] != ('-', '-', '-', '-') and resi[
                        k] in res_map and res_map[resi[k]] in target_res_list and resi[index] in res_map and res_map[
                        resi[index]] in target_res_list and s_mobile[0][res_map[resi[k]][1]][
                        (res_map[resi[k]][2], res_map[resi[k]][0], res_map[resi[k]][3])]['CA'].get_altloc() == \
                            s_mobile[0][res_map[resi[index]][1]][
                                (res_map[resi[index]][2], res_map[resi[index]][0], res_map[resi[index]][3])][
                                'CA'].get_altloc():
                        fixed_ca_list.append(s_mobile[0][res_map[resi[k]][1]][
                                                 (res_map[resi[k]][2], res_map[resi[k]][0], res_map[resi[k]][3])]['CA'])
                        fixed_coord.append(s_mobile[0][res_map[resi[k]][1]][
                                               (res_map[resi[k]][2], res_map[resi[k]][0], res_map[resi[k]][3])][
                                               'CA'].get_vector().get_array())
                        mobile_ca_list.append(s_mobile[0][res_map[resi[index]][1]][(
                        res_map[resi[index]][2], res_map[resi[index]][0], res_map[resi[index]][3])]['CA'])
                        fixed_rl.append(s_mobile[0][res_map[resi[k]][1]][
                                            (res_map[resi[k]][2], res_map[resi[k]][0], res_map[resi[k]][3])].get_id()[
                                            1])
                        mobile_rl.append(s_mobile[0][res_map[resi[index]][1]][(
                        res_map[resi[index]][2], res_map[resi[index]][0], res_map[resi[index]][3])].get_id()[1])
        elif repeats_type[ind] == 'OPEN':
            for k in range(1, rotations):
                for i, resi in enumerate(algn_map):
                    index = k - 1
                    if resi[k] != ('-', '-', '-', '-') and resi[index] != ('-', '-', '-', '-') and resi[
                        k] in res_map and res_map[resi[k]] in target_res_list and resi[index] in res_map and res_map[
                        resi[index]] in target_res_list and s_mobile[0][res_map[resi[k]][1]][
                        (res_map[resi[k]][2], res_map[resi[k]][0], res_map[resi[k]][3])]['CA'].get_altloc() == \
                            s_mobile[0][res_map[resi[index]][1]][
                                (res_map[resi[index]][2], res_map[resi[index]][0], res_map[resi[index]][3])][
                                'CA'].get_altloc():
                        fixed_ca_list.append(s_mobile[0][res_map[resi[k]][1]][
                                                 (res_map[resi[k]][2], res_map[resi[k]][0], res_map[resi[k]][3])]['CA'])
                        fixed_coord.append(s_mobile[0][res_map[resi[k]][1]][
                                               (res_map[resi[k]][2], res_map[resi[k]][0], res_map[resi[k]][3])][
                                               'CA'].get_vector().get_array())
                        mobile_ca_list.append(s_mobile[0][res_map[resi[index]][1]][(
                        res_map[resi[index]][2], res_map[resi[index]][0], res_map[resi[index]][3])]['CA'])
                        fixed_rl.append(s_mobile[0][res_map[resi[k]][1]][
                                            (res_map[resi[k]][2], res_map[resi[k]][0], res_map[resi[k]][3])].get_id()[
                                            1])
                        mobile_rl.append(s_mobile[0][res_map[resi[index]][1]][(
                        res_map[resi[index]][2], res_map[resi[index]][0], res_map[resi[index]][3])].get_id()[1])

        if len(fixed_ca_list) == 0 or len(mobile_ca_list) == 0:
            return {}

        # Superimpose matching residues
        si = Superimposer()
        si.set_atoms(fixed_ca_list, mobile_ca_list)
        T = np.array(np.transpose(si.rotran[0]))
        transl = si.rotran[1]
        ang, rotax, screw_transl, axis_pos = get_rotation_axis(si.rotran[0], transl)
        if len(axis_pos) == 2 and axis_pos == 'na':
            center = sum(fixed_coord) / len(fixed_coord)
            proj = np.dot(center, rotax)
            center_on_axis = np.array([proj * i for i in rotax])
            axis_pos = center - center_on_axis
        #       print("Rotation angle: ",ang*180.0/np.pi,"Axis: ",rotax)
        #       print("Transformation RMS is", si.rms)
        axis_angle = round(angle_between(list(rotax), (0, 0, 1)) * 180 / np.pi, 2)
        angle = ang * 180 / np.pi
        transforms.append(T)
        transl_vectors.append(transl)

        # Apply transformation to coordinates
        aligned_coord_res = []
        for i, resi in enumerate(algn_map):
            if np.any([resi[k] != ('-', '-', '-', '-') for k in range(0, len(resi))]):
                aligned_coord_res.append([])
                for k in range(0, len(resi)):
                    if resi[k] != ('-', '-', '-', '-') and resi[k] in res_map and res_map[resi[k]] in target_res_list:
                        aligned_coord_res[-1].append(res_map[resi[k]])
                        pt = list(s_mobile[0][res_map[resi[k]][1]][
                                      (res_map[resi[k]][2], res_map[resi[k]][0], res_map[resi[k]][3])][
                                      'CA'].get_vector())
                        for j in range(0, k):
                            pt = np.dot((pt - transl), T)
                        s_mobile[0][res_map[resi[k]][1]][
                            (res_map[resi[k]][2], res_map[resi[k]][0], res_map[resi[k]][3])]['CA'].set_coord(pt)
                    else:
                        aligned_coord_res[-1].append(('-', '-', '-', '-'))
        ax2 = axis_pos + 10 * np.array(rotax)
        target_sym_descriptors.append(
            [float('%.2f' % angle), float('%.2f' % screw_transl), transl, ax2, axis_pos, axis_angle])

    #### Calculate the RMSD & TM-score ####
    unique_repeats = list(unique_repeats)
    positions = sum(1 for c in cesm[unique_repeats[0]].seq if c != '/')  # ind=1, positions=97
    algn_map = [[] for i in range(positions)]
    #    coord_repeats=[[] for i in range(positions)]

    # Load the alignments of all repeats of the template against each other
    for k in unique_repeats:
        r = cesm[k]
        r_id = k[len(inputf[0:4]):].split('_')
        ch = r_id[0][-1]
        repeat = r_id[1].strip()
        first_digit = repeat[0]  # handling negative residues
        repeat = repeat[1:].split('-')
        repeat[0] = first_digit + repeat[0]
        #    	res_start.append(repeat[0]) #this is the resid of the residue that starts the repeat
        repeat_start = [i for i, sublist in enumerate(resseq_A) if (int(repeat[0]) in sublist and ch == sublist[
            2])]  # finds the index of the residue that starts the repeat
        #    print(repeat_start)
        count = 0
        for i, c in enumerate(r.seq):  # Last edited
            if c != '/':
                if c != '-' and c.isupper() and resseq_A[repeat_start[0] + count][1] != c:
                    algn_map[i].append(('-', '-', '-', '-'))
                if c != '-' and c.isupper() and resseq_A[repeat_start[0] + count][1] == c:
                    #    		    print(c,resseq_A[repeat_start[0]+count][1],resseq_A[repeat_start[0]+count][0],"\n")
                    algn_map[i].append((resseq_A[repeat_start[0] + count][0], resseq_A[repeat_start[0] + count][2],
                                        resseq_A[repeat_start[0] + count][3], resseq_A[repeat_start[0] + count][4]))
                    count += 1
                #                resi=(resseq_A[repeat_start[0]+count][0],resseq_A[repeat_start[0]+count][2],resseq_A[repeat_start[0]+count][3],resseq_A[repeat_start[0]+count][4])
                #                if resi in res_map and res_map[resi] in target_res_list:
                #                    resit=res_map[resi]
                #                    coord_repeats[i].append(s_mobile[0][resit[1]][(resit[2],resit[0],resit[3])]['CA'].get_vector().get_array())
                #                else:
                #                    coord_repeats[i].append(None)
                #           if resseq_A[repeat_start[0]+count][2] not in chains:
                #               chains.append(resseq_A[repeat_start[0]+count][2])
                if c == '-':
                    algn_map[i].append(('-', '-', '-', '-'))
                #                coord_repeats[i].append(None)
                if c.islower():
                    algn_map[i].append(('-', '-', '-', '-'))
                    #                coord_repeats[i].append(None)
                    count += 1

    num_repeats = len(unique_repeats)
    # Find the length of the target protein
    prot_len = 0
    resseq_B = get_pdb_sequence_with_chains(target)
    prot_len = len(resseq_B)

    rmsd = 0
    tm_score = 0
    # n=0
    dummy = list(range(0, len(unique_repeats)))  # a list of the N numbers representing the N unique repeats
    combo_list = list(combinations(dummy, 2))  # a combination of 2 elements from N
    for combo in combo_list:
        k = combo[0]
        ind = combo[1]
        repeat1 = []
        repeat2 = []
        for i, resi in enumerate(algn_map):
            # Check if the template alignment pair is valid and that both residues have a correspondence in the target
            if resi[k] != ('-', '-', '-', '-') and resi[ind] != ('-', '-', '-', '-') and resi[k] in res_map and res_map[
                resi[k]] in target_res_list and resi[ind] in res_map and res_map[resi[ind]] in target_res_list and \
                    s_mobile[0][res_map[resi[k]][1]][(res_map[resi[k]][2], res_map[resi[k]][0], res_map[resi[k]][3])][
                        'CA'].get_altloc() == s_mobile[0][res_map[resi[ind]][1]][
                (res_map[resi[ind]][2], res_map[resi[ind]][0], res_map[resi[ind]][3])][
                'CA'].get_altloc():  # checks that they are coordinates rather than ['-']
                # print("[",res_map[resi[k]],",",res_map[resi[ind]],"] ",str(pow(np.linalg.norm(s_mobile[0][res_map[resi[k]][1]][(res_map[resi[k]][2],res_map[resi[k]][0],res_map[resi[k]][3])]['CA'].get_vector()-s_mobile[0][res_map[resi[ind]][1]][res_map[resi[ind]][0]]['CA'].get_vector()),2)))
                # rmsd=rmsd+pow(np.linalg.norm(s_mobile[0][res_map[resi[k]][1]][(res_map[resi[k]][2],res_map[resi[k]][0],res_map[resi[k]][3])]['CA'].get_vector()-s_mobile[0][res_map[resi[ind]][1]][res_map[resi[ind]][0]]['CA'].get_vector()),2)
                repeat1.append(
                    s_mobile[0][res_map[resi[k]][1]][(res_map[resi[k]][2], res_map[resi[k]][0], res_map[resi[k]][3])][
                        'CA'].get_vector().get_array())
                repeat2.append(s_mobile[0][res_map[resi[ind]][1]][
                                   (res_map[resi[ind]][2], res_map[resi[ind]][0], res_map[resi[ind]][3])][
                                   'CA'].get_vector().get_array())
                # n+=1
        # print(repeat1, repeat2)
        tm_score = tm_score + tm_score_calc(repeat1, repeat2, prot_len)
    tm_score = tm_score * num_repeats / len(combo_list)
    # rmsd=np.sqrt(rmsd/n)

    # Store coordinates for RMSD procedure
    coord_repeats = []
    altlocs = []
    for i, resi in enumerate(algn_map):
        # Check if the template alignment pair is valid and that both residues have a correspondence in the target
        coord_repeats.append([])
        altlocs.append([])
        for k in range(0, len(resi)):
            if resi[k] != ('-', '-', '-', '-') and resi[k] in res_map and res_map[resi[k]] in target_res_list:
                resit = res_map[resi[k]]
                coord_repeats[i].append(
                    s_mobile[0][resit[1]][(resit[2], resit[0], resit[3])]['CA'].get_vector().get_array())
                altlocs[i].append(s_mobile[0][resit[1]][(resit[2], resit[0], resit[3])]['CA'].get_altloc())
            else:
                coord_repeats[i].append(None)
                altlocs[i].append(None)

    tmp = list(map(list, zip(*coord_repeats)))
    coord_repeats = tmp
    tmp = list(map(list, zip(*altlocs)))
    altlocs = tmp
    rmsd = get_rmsd_cesymm(coord_repeats, altlocs)  # Java RMSD
    #    print("\nEstimated RMSD is: ",rmsd)
    #    print("Estimated TM-score is: ",tm_score)
    # os.remove(oriented_enc)
    # if os.path.isfile(oriented_enc_2):
    #    os.remove(oriented_enc_2)

    #### Find the leaves of the tree (all derivative axes) - does not work for dihedral cases with order >2
    algn_map = [[] for i in range(positions)]
    for k in unique_repeats:
        r = cesm[k]
        r_id = k.split('_')
        ch = r_id[0][-1]
        repeat = r_id[1].strip()
        first_digit = repeat[0]  # handling negative residues
        repeat = repeat[1:].split('-')
        repeat[0] = first_digit + repeat[0]

        repeat_start = [i for i, sublist in enumerate(resseq_A) if (int(repeat[0]) in sublist and ch == sublist[
            2])]  # finds the index of the residue that starts the repeat
        count = 0
        for i, c in enumerate(r.seq):
            if c != '/':
                if c != '-' and c.isupper() and resseq_A[repeat_start[0] + count][1] != c:
                    algn_map[i].append(('-', '-', '-', '-'))
                if c != '-' and c.isupper() and resseq_A[repeat_start[0] + count][1] == c:
                    algn_map[i].append((resseq_A[repeat_start[0] + count][0], resseq_A[repeat_start[0] + count][2],
                                        resseq_A[repeat_start[0] + count][3], resseq_A[repeat_start[0] + count][4]))
                    count += 1
                if c == '-':
                    algn_map[i].append(('-', '-', '-', '-'))
                if c.islower():
                    algn_map[i].append(('-', '-', '-', '-'))
                    count += 1
    aligned = 0
    target_algn_map = []
    for i, resi in enumerate(algn_map):
        # Check if the template alignment pair is valid and that both residues have a correspondence in the target
        target_algn_map.append([])
        for k in range(0, len(resi)):
            if resi[k] != ('-', '-', '-', '-') and resi[k] in res_map and res_map[resi[k]] in target_res_list:
                target_algn_map[i].append(res_map[resi[k]])
                aligned += 1
            else:
                target_algn_map[i].append(('-', '-', '-', '-'))
    tmp = []  # Remove columns that are entirely made up of gaps
    for i, resi in enumerate(target_algn_map):
        if all(x == ('-', '-', '-', '-') for x in resi) == False:
            tmp.append(resi)
    target_algn_map = tmp
    del tmp

    target_repeats = ['' for i in range(0, len(target_algn_map[0]))]
    target_repeats_selection = []  # [[(9,'A'),(235,'A'),(1,'B'),(235,'B'),(1,'C'),(235,'C')],[...]],
    count = 0
    beg_resid_target = 10000
    end_resid_target = -10
    target_chain = ''
    for j in range(0, len(target_algn_map[0])):
        flag = 0
        for i, resi in enumerate(target_algn_map):
            if flag == 0 and resi[j] != ('-', '-', '-', '-') and resi.count(('-', '-', '-', '-')) <= len(resi) - 2:
                target_repeats[j] = target_pdb[0:4] + '.' + resi[j][1] + "_" + str(resi[j][0])
                target_repeats_selection.append([(resi[j][0], resi[j][1])])
                if resi[j][0] < beg_resid_target:
                    beg_resid_target = resi[j][0]
                    target_chain = resi[j][1]
                flag = 1
            if flag == 1 and resi[j] != ('-', '-', '-', '-') and resi.count(('-', '-', '-', '-')) <= len(resi) - 2:
                last = resi[j]
                if last[0] > end_resid_target:
                    end_resid_target = last[0]
        target_repeats[j] = target_repeats[j] + '-' + str(last[0])
        target_repeats_selection[j].append((last[0], last[1]))
    #    print(target_repeats)
    target_repeats_all = []
    target_repeats_type = repeats_type
    # Example: repeats_all=[[['1PV6.A_6-94', '1PV6.A_219-309'],['1PV6.A_101-187', '1PV6.A_310-399']],[['1PV6.A_6-94', '1PV6.A_101-187'],['1PV6.A_219-309', '1PV6.A_310-399']]]
    for ind, symm in enumerate(repeats_all):
        target_repeats_all.append([])
        for i, sub_symm in enumerate(symm):  # sub_symm=['1PV6.A_6-94', '1PV6.A_219-309']
            target_repeats_all[ind].append([])
            for k in sub_symm:
                pos = unique_repeats.index(k)
                target_repeats_all[ind][i].append(target_repeats[pos])

    print('target_repeats_all: ', target_repeats_all)
    super_pml_files = ""
    for ind, symm in enumerate(target_repeats_all):
        # Create PyMOL superposition scripts
        cb_colors, cb_names = colorblind_palette()
        pml_colors, pml_names = additional_colors_palette()
        colors = cb_names + pml_names
        out = open(super_pml_dir + "superposition_" + image_key + "_level_" + str(ind + 1) + ".pml", "w")
        super_pml_files = super_pml_files + "superposition_" + image_key + "_level_" + str(ind + 1) + ".pml;"
        out.write("load " + target_pdb[0:4] + ".pdb, " + target_pdb + "\n")
        out.write("set antialias,1\nset line_smooth, 1\nset depth_cue, 0.5\nset specular, 0\nset surface_quality, 1\n")
        out.write("set cartoon_sampling, 14\nset ribbon_sampling, 10\nset ray_trace_fog, 0.5\n")
        out.write("viewport 350,350\n")
        out.write("set ignore_case, off\n")
        for c in cb_colors:
            out.write("set_color cb_" + c[0] + "=[" + str(c[1][0]) + "," + str(c[1][1]) + "," + str(c[1][2]) + "]\n")
        out.write("hide everything,*\n")
        out.write("show cartoon\n")
        out.write("bg_color white\nset opaque_background, on\n")
        out.write("import numpy as np\n")

        for m in range(0, ind + 1):
            out.write("transl" + str(m) + "=np.array([" + str(transl_vectors[m][0]) + "," + str(
                transl_vectors[m][1]) + "," + str(transl_vectors[m][2]) + "])\n")
            out.write(
                "T" + str(m) + "=np.array([[" + str(transforms[m][0][0]) + "," + str(transforms[m][0][1]) + "," + str(
                    transforms[m][0][2]) + "],[" + str(transforms[m][1][0]) + "," + str(
                    transforms[m][1][1]) + "," + str(transforms[m][1][2]) + "],[" + str(
                    transforms[m][2][0]) + "," + str(transforms[m][2][1]) + "," + str(transforms[m][2][2]) + "]])\n")

        k = 0
        for j in range(0, len(symm[0])):
            for s in range(0, len(symm)):
                sel_text = ""
                entry = symm[s][j]
                print(entry)
                ch = entry.split('.')[1][0]
                beg = entry.split('_')[1]
                first_digit = beg[0]  # in case the residues has a negative id (i.e. -1)
                beg = beg[1:].split('-')[0]
                beg = first_digit + beg
                locator = [i for i, sublist in enumerate(target_repeats_selection) if ((int(beg), ch) == sublist[0])]
                for i in range(0, int(len(target_repeats_selection[locator[0]]) / 2)):
                    sel_text = sel_text + "(chain " + target_repeats_selection[locator[0]][i * 2][
                        1] + " and resid " + str(target_repeats_selection[locator[0]][i * 2][0]) + "-" + str(
                        target_repeats_selection[locator[0]][i * 2 + 1][0]) + ")"
                    if i < len(target_repeats_selection[locator[0]]) / 2 - 1:
                        sel_text = sel_text + " or "
                # if s<(len(repeats)-1):
                #    sel_text=sel_text+" or "
                u = unique_repeats.index(repeats_all[ind][s][j])
                out.write("create repeat" + str(u + 1) + ", " + target_pdb + " and " + sel_text + "\n")
                out.write("color " + colors[(k % len(colors))] + ",repeat" + str(u + 1) + "\n")
                for m in range(0, ind + 1):
                    for p in range(0, repeats_transf[m][u]):
                        out.write(
                            "alter_state 1, repeat" + str(u + 1) + ", (x,y,z)=np.dot((np.array([x,y,z])-transl" + str(
                                m) + "),T" + str(m) + ")\n")

            k += 1
        out.write("deselect\n")
        out.write("rotate [1,0,0], angle=-90\n")
        out.write("zoom " + target_pdb + "\n")
        out.write("delete " + target_pdb + "\n")
        #        out.write("cd "+transfer_out_dir+"\n")
        #        out.write("png "+target_pdb+"_transfer_super_"+str(ind)+".png\n")
        #        out.write("quit")
        out.close()

    #    print(target_repeats_all, repeats_type, repeats_levels, axes_per_level)
    target_repeats_cyclic = []
    for ind, symm in enumerate(target_repeats_all):
        pairs_per_axis = len(symm) / axes_per_level[ind]
        target_repeats_cyclic.append([])
        str_rep = ""
        for i, sub_symm in enumerate(symm):
            for j, k in enumerate(sub_symm):
                if k[4:5] == '_':
                    sub_symm[j] = k[0:4].upper() + '.' + k[5:]
            if i % pairs_per_axis != 0:
                str_rep = str_rep + "(" + ";".join(sub_symm) + ")"
            elif i != 0:
                target_repeats_cyclic[ind].append(str_rep)
                str_rep = "(" + ";".join(sub_symm) + ")"
            elif i == 0:
                str_rep = "(" + ';'.join(sub_symm) + ")"
        target_repeats_cyclic[ind].append(str_rep)
    #    print(transforms)
    # Build tree of transformations (only binary tree supported at this stage)
    s = [np.array([i]) for i in reversed(range(0, len(transforms)))]
    all_transforms = [[]]
    for t in s:
        all_transforms_tmp = []
        for i, ax in enumerate(all_transforms):
            all_transforms_tmp.append(list(ax) + list(t))
        all_transforms = all_transforms + all_transforms_tmp

    all_transforms = all_transforms[1:]

    levels = []
    half = len(all_transforms) + 1
    for i in range(0, len(s)):
        subs = all_transforms[0:(half - 1)]
        half = int((len(subs) + 1) / 2)
        spacing = len(subs) + 1
        tmp = [subs[half - 1]]
        jump = spacing
        while half + jump <= len(all_transforms):
            tmp.append(all_transforms[half - 1 + jump])
            jump += spacing
        levels.append(tmp)

        # Calculate axes from tree
    all_axes = []
    for ind, l in enumerate(levels):
        axes_tmp = []
        for sub_ind, j in enumerate(l):
            for i, t in enumerate(j):
                if i == 0:
                    pt = target_sym_descriptors[t][3]
                    pt2 = target_sym_descriptors[t][4]
                if i > 0:
                    transl = target_sym_descriptors[t][2]
                    T = transforms[t]
                    pt = np.dot((pt - transl), T)
                    pt2 = np.dot((pt2 - transl), T)
            descriptor = [target_sym_descriptors[ind][0], target_sym_descriptors[ind][1],
                          target_sym_descriptors[ind][2], pt, pt2, target_sym_descriptors[ind][5],
                          target_repeats_cyclic[ind][sub_ind]]
            axes_tmp.append(descriptor)
        all_axes.append(axes_tmp)

    # Address cases with independent symmetries
    if any(char.isdigit() for char in inputf[6:]):
        """
        beg_resid=inputf[7:].split('-')[0]
        end_resid=inputf[7:].split('-')[1]
        chain=inputf[5:6]
        beg_resid_target=res_map[(int(beg_resid),chain,' ',' ')][0]
        end_resid_target=res_map[(int(end_resid),chain,' ',' ')][0]
        """
        target_pdb = target_pdb + '_' + str(beg_resid_target) + '-' + str(end_resid_target)

    # Print axes
    axout = transfer_out_dir + target_pdb + "_transfer.axes"
    output = open(axout, "w")
    closed_opened = ""
    angle = ""
    translation = ""
    axis_angle_with_membrane_normal = ""
    topology = ""
    repeats_web = ""
    super_pml_files_all = ""
    output.write(
        "Name    SymmLevel       SymmType        SymmOrder       RotationAngle   ScrewTranslation        Point1  Point2  AlignedRepeats\n")
    for ind, symm in enumerate(all_axes):
        for i, sub_symm in enumerate(symm):
            pt1 = all_axes[ind][i][3]
            pt2 = all_axes[ind][i][4]
            rotax = unit_vector(np.array(pt2) - np.array(pt1))
            axis_angle = round(angle_between(list(rotax), (0, 0, 1)) * 180 / np.pi, 2)
            axis_angle_with_membrane_normal = axis_angle_with_membrane_normal + str(float('%.2f' % axis_angle)) + ";"
            if np.abs(np.dot(rotax, (0, 0, 1))) > 0.5:
                topology = topology + "Parallel;"
            else:
                topology = topology + "Antiparallel;"
            angle = angle + str(all_axes[ind][i][0]) + ";"
            translation = translation + str(all_axes[ind][i][1]) + ";"
            closed_opened = closed_opened + repeats_type[ind] + ";"
            output.write("%s\t%d\t%s\t%d\t%.3f\t%.3f\t%.3f,%.3f,%.3f\t%.3f,%.3f,%.3f\t%s\n" % (
            target_pdb, repeats_levels[ind], repeats_type[ind], int(order_of_level[ind]), all_axes[ind][i][0],
            all_axes[ind][i][1], pt1[0], pt1[1], pt1[2], pt2[0], pt2[1], pt2[2], all_axes[ind][i][6]))
            repeats_web = repeats_web + all_axes[ind][i][6].replace(";", ",") + ";"
            name = [k for k in super_pml_files.split(';')[:-1] if 'level_' + str(repeats_levels[ind]) + '.pml' in k][0]
            super_pml_files_all = super_pml_files_all + name.strip() + ';'
    repeats_web = repeats_web.replace(target_pdb.upper() + '.', '')
    output.close()

    # Print alignment
    alnres = transfer_out_dir + target_pdb + "_transfer.alnres"
    output = open(alnres, "w")
    tags = []
    seqs = []
    resids = []
    for ind, rep in enumerate(target_repeats):
        if rep[4:5] == '_':
            rep = rep[0:4].upper() + '.' + rep[5:]
        tags.append(">" + rep)
        seq = []
        for r in target_algn_map:
            if r[ind] != ('-', '-', '-', '-'):
                for i in resseq_B:
                    if r[ind][0] == i[0] and r[ind][1] == i[2]:
                        seq.append(i[1])
            else:
                seq.append('-')
        seqs.append(''.join(seq))
        resid = ["{:<3s}".format(str(i[ind][0])) for i in target_algn_map]
        resids.append(','.join(resid))
        output.write(tags[ind] + "\n")
        output.write(seqs[ind] + "\n")
    output.write("\n")
    for ind in range(0, len(target_repeats)):
        output.write(tags[ind] + ":resid\n")
        output.write(resids[ind] + "\n")
    output.close()
    rep_length = len(target_algn_map)
    num_repeats = len(target_repeats)
    symmetry_levels = len(repeats_levels)
    if os.path.isfile(wkdir + "dummy.txt"):
        os.remove(wkdir + "dummy.txt")
    trans_dic = {}
    trans_dic['template_order'] = ce_order
    trans_dic['coverage'] = str(float('%.2f' % (aligned / prot_len)))
    trans_dic['size'] = str(prot_len)
    trans_dic['unit_angle'] = angle
    trans_dic['unit_translation'] = translation
    trans_dic['refined_tmscore'] = str(float('%.2f' % tm_score))
    trans_dic['refined_rmsd'] = str(float('%.2f' % rmsd))
    trans_dic['repeats_number'] = str(num_repeats)
    trans_dic['topology'] = topology
    trans_dic['axis_angle_with_membrane_normal'] = axis_angle_with_membrane_normal
    trans_dic['symmetry_levels'] = str(symmetry_levels)
    trans_dic['aligned_length'] = str(aligned)
    trans_dic['repeat_length'] = str(rep_length)
    trans_dic['template_structure'] = parent_id
    trans_dic['files_key'] = target_pdb + "_transfer"
    trans_dic['files_dir_path'] = transfer_out_dir
    trans_dic['axes_file'] = axout
    trans_dic['alnres_file'] = alnres
    #    trans_dic['out_file']=out
    #    trans_dic['parents_file']=parents
    trans_dic['repeats'] = repeats_web
    trans_dic['repeats_selection'] = target_repeats_selection
    trans_dic['super_pml_files'] = super_pml_files_all
    #    return float('%.2f'%tm_score), float('%.2f'%rmsd), axout, alnres, out, parent
    return trans_dic


def is_symmetry_in_membrane(limit_bottom, limit_top, repeat_selection, oriented):
    """
    Checks whether there are at least 40 a.a. from any repeat that are in the membrane
    """
    f = open(oriented, "r")
    count = 0
    for i in range(0, len(repeat_selection)):
        count = 0
        for j in range(0, int(len(repeat_selection[i]) / 2)):
            f.seek(0)
            for line in f:
                if line.startswith("ATOM") and line[21:22].strip() == repeat_selection[i][j * 2][1] and \
                        repeat_selection[i][j * 2 + 1][0] >= int(re.search(r'\d+', line[22:27]).group()) >= \
                        repeat_selection[i][j * 2][0] and limit_bottom <= float(line[46:54]) <= limit_top:
                    count = count + 1
            if count >= 40:
                break
        if count >= 40:
            break
    f.close()
    if count >= 40:
        return 'in'
    else:
        return 'out'


def symmetry_tm_domains(tm_archive, repeat_selection, pdb, chain):
    """
    Checks whether any repeat has at least two transmembrane crossings 
    repeat_selction=[[(9,'A'),(235,'A'),(1,'B'),(235,'B'),(1,'C'),(235,'C')],[...]], so that within each [] we have the ranges 
    domains=[['A',['1-20','30-45',..]],['B',['5-23','30-42',...]]]
    """
    f = pkl.load(open(tm_archive, 'rb'))
    print(repeat_selection)
    domains = []
    for entry in f:
        if entry == pdb:
            if chain in f[pdb]:
                # old EncoMPASS archives and MemSTATS use tmdoms, new EncoMPASS archives use TM_regions
                #print("check here:", f[pdb][chain])
                if 'tmdoms' in f[pdb][chain]:
                    tmdoms = [str(i[0]) + '-' + str(i[1]) for i in f[pdb][chain]['tmdoms']]
                else:
                    tmdoms = [str(i[0][0][0]) + '-' + str(i[1][0][0]) for i in
                              f[pdb][chain]['TM_regions']['TM_regions_extrema']]
                domains.append([chain, tmdoms])
            if chain == "":
                for ch in f[pdb]['tmchains']:
                    if 'tmdoms' in f[pdb][ch]:
                        tmdoms = [str(i[0]) + '-' + str(i[1]) for i in f[pdb][ch]['tmdoms']]
                    else:
                        tmdoms = [str(i[0][0][0]) + '-' + str(i[1][0][0]) for i in
                                  f[pdb][ch]['TM_regions']['TM_regions_extrema']]
                    domains.append([ch, tmdoms])
    #    print(domains)
    #    print(repeat_selection)
    all_crossings = []
    crossings = 0
    tot_crossings = 0
    for i in range(0, len(repeat_selection)):
        crossings = 0
        for j in range(0, int(len(repeat_selection[i]) / 2)):
            reps = set(range(repeat_selection[i][j * 2][0], repeat_selection[i][j * 2 + 1][0] + 1))
            doms = set()
            for ind, ch_doms in enumerate(domains):
                if ch_doms[0] == repeat_selection[i][j * 2][1]:
                    for dom in ch_doms[1]:
                        rg = dom.split('-')
                        if len(rg) == 2:
                            doms = set(range(int(rg[0]), int(rg[1]) + 1))
                        overlap = len(reps.intersection(doms))
                        if overlap > 5:
                            crossings += 1
        tot_crossings = tot_crossings + crossings
        all_crossings.append(crossings)
    print("All crossings: ", all_crossings, "Total crossings: ", tot_crossings)

    if any(s > 1 for s in all_crossings) == True and tot_crossings > 3:
        return 'in', all_crossings
    else:
        return 'out', all_crossings


def repeats_overlap(new_selection, accepted_selection):
    # we assume that we are dealing with a single-chain protein/fragmet
    print('Current function name:', sys._getframe().f_code.co_name)
    print(new_selection, accepted_selection)
    overlap = 0
    for i in range(0, len(new_selection)):
        xs = set(range(new_selection[i][0][0], new_selection[i][1][0] + 1))
        for j in range(0, len(accepted_selection)):
            yr = set(range(accepted_selection[j][0][0], accepted_selection[j][1][0] + 1))
            overlap = len(xs.intersection(yr))
            if overlap != 0:
                break
        if overlap != 0:
            break

    if overlap != 0:
        print("New and accepted repeats overlap.")
        return 'overlap'
    else:
        print("New and accepted repeats are independent.")
        return 'independent'


def clear_disordered_atoms(oriented):
    s = parse_structure(oriented)
    io = PDBIO()

    class NotDisordered(Select):
        def accept_atom(self, atom):
            return (not atom.is_disordered() or atom.get_altloc() == 'A') and atom.get_parent().get_id()[2] == ' '

    io = PDBIO()
    io.set_structure(s)
    io.save(oriented, select=NotDisordered())
    return


def repeats_superposition(cesymm_out_dir, oriented, inputf, super_dir, image_key, repeats_selection):
    """ Create a PyMOL script that superimposes the repeats of each symmetry level """
    s_mobile = parse_structure(oriented)
    mobile = s_mobile[0]

    super_pml_files = ""
    # repeats_selection=get_repeat_resid(cesymm_out_dir,inputf,mobile)
    resseq_A = get_pdb_sequence_with_chains(
        mobile)  # creates a list of the form [(2, 'D', 'A'), (3, 'Q', 'A'), (4, 'A', 'A'), (5, 'L', 'A'), (6, 'S', 'A')...]
    cesm = SeqIO.parse(cesymm_out_dir + inputf + ".fasta", "fasta")
    unique_repeats = []
    for i in list(cesm):
        unique_repeats.append(i.id)
    cesm = SeqIO.to_dict(SeqIO.parse(cesymm_out_dir + inputf + ".fasta", "fasta"))
    repeats_all, repeats_type, repeats_level, axes_per_level, order_of_level = cesymm_related_repeats(inputf,
                                                                                                      cesymm_out_dir)
    rot_matrix = []
    transl_vector = []
    repeats_transf = [[-1 for j in range(0, len(unique_repeats))] for i in range(0,
                                                                                 len(repeats_all))]  # same number of elements as symmetry levels, eg. [[-1, -1, -1, -1], [-1, -1, -1, -1]]

    # Example: repeats_all=[[['1PV6.A_6-94', '1PV6.A_219-309'],['1PV6.A_101-187', '1PV6.A_310-399']],[['1PV6.A_6-94', '1PV6.A_101-187'],['1PV6.A_219-309', '1PV6.A_310-399']]]
    for ind, symm in enumerate(repeats_all):
        tag = [t for t in cesm for k in range(0, len(symm)) if
               symm[k][1] in t]  # get the fasta file identifier of each sequence in the repeat
        positions = sum([sum(1 for c in cesm[k].seq if c != '/') for k in set(tag)])
        # positions=sum([sum(1 for c in cesm[symm[k][1]].seq if c!='/') for k in range(0,len(symm))])
        algn_map = [[] for i in range(0, positions)]
        map_ind = 0
        for sub_symm in symm:  # sub_symm=['1PV6.A_6-94', '1PV6.A_219-309']
            for l, k in enumerate(sub_symm):  # k='1PV6.A_6-94'
                try:
                    r = cesm[k]
                    r_id = k[len(inputf[0:4]):].split('_')
                    ch = r_id[0][-1]
                except:
                    tag = [t for t in cesm if k in t]
                    r = cesm[tag[0]]
                    r_id = k.split('_')
                    ch = r_id[0]
                repeat = r_id[1].strip()
                first_digit = repeat[0]  # handling negative residues
                repeat = repeat[1:].split('-')
                repeat[0] = first_digit + repeat[0]
                repeat_start = [i for i, sublist in enumerate(resseq_A) if (int(repeat[0]) in sublist and ch == sublist[
                    2])]  # finds the index of the residue that starts the repeat
                count = 0
                for i, c in enumerate(r.seq):  # Last edited
                    if c != '/':
                        if c != '-' and c.isupper() and resseq_A[repeat_start[0] + count][1] == c:
                            algn_map[i + map_ind].append((resseq_A[repeat_start[0] + count][0],
                                                          resseq_A[repeat_start[0] + count][2],
                                                          resseq_A[repeat_start[0] + count][3],
                                                          resseq_A[repeat_start[0] + count][4]))
                            count += 1
                        if c == '-':
                            algn_map[i + map_ind].append(('-', '-', '-', '-'))
                        if c.islower():
                            algn_map[i + map_ind].append(('-', '-', '-', '-'))
                            count += 1
                try:
                    u = unique_repeats.index(k)
                except:
                    u = unique_repeats.index(inputf[0:4].upper() + '.' + k)
                repeats_transf[ind][u] = l
            try:
                map_ind += sum(1 for c in cesm[k].seq if c != '/')
            except:
                map_ind += sum(1 for c in cesm[inputf[0:4].upper() + '.' + k].seq if c != '/')

        rotations = 2
        tm_score = 0
        rmsd = 0
        for resi in algn_map:
            if all(elem == ('-', '-', '-', '-') for elem in resi) == False:
                rotations = len(resi)
                break

        # Create a list of atoms to be paired during superposition for obtaining the axis of rotation
        fixed_ca_list, mobile_ca_list = [], []
        fixed_rl = []
        mobile_rl = []

        if repeats_type[ind] == 'CLOSED':
            for k in range(0, rotations):
                for i, resi in enumerate(algn_map):
                    index = (rotations - 1 + k) % rotations
                    if resi[k] != ('-', '-', '-', '-') and resi[index] != ('-', '-', '-','-'):  # and s_mobile[0][resi[k][1]][(resi[k][2],resi[k][0],resi[k][3])]['CA'].get_altloc()==s_mobile[0][resi[index][1]][(resi[index][2],resi[index][0],resi[index][3])]['CA'].get_altloc():
                        fixed_ca_list.append(s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'])
                        mobile_ca_list.append(
                            s_mobile[0][resi[index][1]][(resi[index][2], resi[index][0], resi[index][3])]['CA'])
                        fixed_rl.append(s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])].get_id()[1])
                        mobile_rl.append(
                            s_mobile[0][resi[index][1]][(resi[index][2], resi[index][0], resi[index][3])].get_id()[1])
        elif repeats_type[ind] == 'OPEN':
            for k in range(1, rotations):
                for i, resi in enumerate(algn_map):
                    index = k - 1
                    if resi[k] != ('-', '-', '-', '-') and resi[index] != ('-', '-', '-', '-') and \
                            s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'].get_altloc() == \
                            s_mobile[0][resi[index][1]][(resi[index][2], resi[index][0], resi[index][3])][
                                'CA'].get_altloc():
                        fixed_ca_list.append(s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'])
                        mobile_ca_list.append(
                            s_mobile[0][resi[index][1]][(resi[index][2], resi[index][0], resi[index][3])]['CA'])
                        fixed_rl.append(s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])].get_id()[1])
                        mobile_rl.append(
                            s_mobile[0][resi[index][1]][(resi[index][2], resi[index][0], resi[index][3])].get_id()[1])

        # Superimpose matching residues
        si = Superimposer()
        si.set_atoms(fixed_ca_list, mobile_ca_list)
        T = np.array(np.transpose(si.rotran[0]))
        transl = si.rotran[1]
        rot_matrix.append(T)
        transl_vector.append(transl)

        # Apply transformation to coordinates
        aligned_coord_res = []
        for i, resi in enumerate(algn_map):
            if np.any([resi[k] != ('-', '-', '-', '-') for k in range(0, len(resi))]):
                aligned_coord_res.append([])
                for k in range(0, len(resi)):
                    if resi[k] != ('-', '-', '-', '-'):
                        aligned_coord_res[-1].append(resi[k])
                        pt = list(s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'].get_vector())
                        for j in range(0, k):
                            pt = np.dot((pt - transl), T)
                        s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'].set_coord(pt)
                    else:
                        aligned_coord_res[-1].append(('-', '-', '-', '-'))

        # Generate PyMOL superposition script
        # repeats_all=[[['1PV6.A_6-94', '1PV6.A_219-309'],['1PV6.A_101-187', '1PV6.A_310-399']],[['1PV6.A_6-94', '1PV6.A_101-187'],['1PV6.A_219-309', '1PV6.A_310-399']]]

        cb_colors, cb_names = colorblind_palette()
        pml_colors, pml_names = additional_colors_palette()
        colors = cb_names + pml_names
        out = open(super_dir + "superposition_" + image_key + "_level_" + str(ind + 1) + ".pml", "w")
        super_pml_files = super_pml_files + "superposition_" + image_key + "_level_" + str(ind + 1) + ".pml;"
        out.write("load " + inputf[0:4] + ".pdb, " + inputf + "\n")
        out.write("set antialias,1\nset line_smooth, 1\nset depth_cue, 0.5\nset specular, 0\nset surface_quality, 1\n")
        out.write("set cartoon_sampling, 14\nset ribbon_sampling, 10\nset ray_trace_fog, 0.5\n")
        out.write("viewport 350,350\n")
        out.write("set ignore_case, off\n")
        for c in cb_colors:
            out.write("set_color cb_" + c[0] + "=[" + str(c[1][0]) + "," + str(c[1][1]) + "," + str(c[1][2]) + "]\n")
        out.write("hide everything,*\n")
        out.write("show cartoon\n")
        out.write("bg_color white\nset opaque_background, on\n")
        out.write("import numpy as np\n")
        # Find the correct residues to select for the repeat that is connected with each particular symmetry axis;
        # Note that we are not selecting each CE-Symm defined repeat but rather the repeats that correspond to each axis
        # So in (2A65.A_9-235;2A65.B_9-235)(2A65.A_242-441;2A65.B_242-441)
        # Repeat 1 is 2A65.A_9-235 and 2A65.A_242-441
        # Repeat 2 is 2A65.B_9-235 and 2A65.B_242-441
        # repeats_all=[[['1PV6.A_6-94', '1PV6.A_219-309'],['1PV6.A_101-187', '1PV6.A_310-399']],[['1PV6.A_6-94', '1PV6.A_101-187'],['1PV6.A_219-309', '1PV6.A_310-399']]]

        for m in range(0, ind + 1):
            out.write("transl" + str(m) + "=np.array([" + str(transl_vector[m][0]) + "," + str(
                transl_vector[m][1]) + "," + str(transl_vector[m][2]) + "])\n")
            out.write(
                "T" + str(m) + "=np.array([[" + str(rot_matrix[m][0][0]) + "," + str(rot_matrix[m][0][1]) + "," + str(
                    rot_matrix[m][0][2]) + "],[" + str(rot_matrix[m][1][0]) + "," + str(
                    rot_matrix[m][1][1]) + "," + str(rot_matrix[m][1][2]) + "],[" + str(
                    rot_matrix[m][2][0]) + "," + str(rot_matrix[m][2][1]) + "," + str(rot_matrix[m][2][2]) + "]])\n")

        k = 0
        for j in range(0, len(symm[0])):
            for s in range(0, len(symm)):
                sel_text = ""
                entry = symm[s][j]
                try:
                    ch = entry.split('.')[1][0]
                except:
                    ch = entry.split('_')[0]
                beg = entry.split('_')[1]
                first_digit = beg[0]  # in case the residues has a negative id (i.e. -1)
                beg = beg[1:].split('-')[0]
                beg = first_digit + beg
                locator = [i for i, sublist in enumerate(repeats_selection) if ((int(beg), ch) == sublist[0])]
                for i in range(0, int(len(repeats_selection[locator[0]]) / 2)):
                    sel_text = sel_text + "(chain " + repeats_selection[locator[0]][i * 2][1] + " and resid " + str(
                        repeats_selection[locator[0]][i * 2][0]) + "-" + str(
                        repeats_selection[locator[0]][i * 2 + 1][0]) + ")"
                    if i < len(repeats_selection[locator[0]]) / 2 - 1:
                        sel_text = sel_text + " or "
                # if s<(len(repeats)-1):
                #    sel_text=sel_text+" or "
                try:
                    u = unique_repeats.index(entry)
                except:
                    u = unique_repeats.index(inputf[0:4].upper() + "." + entry)
                out.write("create repeat" + str(u + 1) + ", " + inputf + " and " + sel_text + "\n")
                out.write("color " + colors[(k % len(colors))] + ",repeat" + str(u + 1) + "\n")
                for m in range(0, ind + 1):
                    for p in range(0, repeats_transf[m][u]):
                        out.write(
                            "alter_state 1, repeat" + str(u + 1) + ", (x,y,z)=np.dot((np.array([x,y,z])-transl" + str(
                                m) + "),T" + str(m) + ")\n")

            k += 1
        out.write("deselect\n")
        out.write("rotate [1,0,0], angle=-90\n")
        out.write("zoom " + inputf + "\n")
        out.write("delete " + inputf + "\n")
        #        out.write("cd "+super_dir+"\n")
        #        out.write("png "+image_key+"_super_"+str(ind)+".png\n")
        #        out.write("quit")
        out.close()
    #        os.system(pymol+" "+super_dir+"level_"+str(ind)+"_superposition_"+image_key+".pml")

    #### Calculate the RMSD & TM-score ####
    unique_repeats = list(unique_repeats)
    positions = sum(1 for c in cesm[unique_repeats[0]].seq if c != '/')  # ind=1, positions=97
    algn_map = [[] for i in range(positions)]

    # Load the alignments of all repeats of the template against each other
    for k in unique_repeats:
        try:
            r = cesm[k]
            r_id = k[len(inputf[0:4]):].split('_')
            ch = r_id[0][-1]
        except:
            tag = [t for t in cesm if k in t]
            r = cesm[tag[0]]
            r_id = k.split('_')
            ch = r_id[0]
        repeat = r_id[1].strip()
        first_digit = repeat[0]  # handling negative residues
        repeat = repeat[1:].split('-')
        repeat[0] = first_digit + repeat[0]
        repeat_start = [i for i, sublist in enumerate(resseq_A) if (int(repeat[0]) in sublist and ch == sublist[
            2])]  # finds the index of the residue that starts the repeat
        count = 0
        for i, c in enumerate(r.seq):  # Last edited
            if c != '/':
                if c != '-' and c.isupper() and resseq_A[repeat_start[0] + count][1] == c:
                    algn_map[i].append((resseq_A[repeat_start[0] + count][0], resseq_A[repeat_start[0] + count][2],
                                        resseq_A[repeat_start[0] + count][3], resseq_A[repeat_start[0] + count][4]))
                    count += 1

                if c == '-':
                    algn_map[i].append(('-', '-', '-', '-'))
                if c.islower():
                    algn_map[i].append(('-', '-', '-', '-'))
                    count += 1

    num_repeats = len(unique_repeats)
    # Find the length of the target protein
    prot_len = 0
    prot_len = len(resseq_A)

    #    rmsd=0
    tm_score = 0
    #    n=0
    dummy = list(range(0, num_repeats))  # a list of the N numbers representing the N unique repeats
    combo_list = list(combinations(dummy, 2))  # a combination of 2 elements from N
    for combo in combo_list:
        k = combo[0]
        ind = combo[1]
        repeat1 = []
        repeat2 = []
        for i, resi in enumerate(algn_map):
            # Check if the template alignment pair is valid and that both residues have a correspondence in the target
            if resi[k] != ('-', '-', '-', '-') and resi[ind] != ('-', '-', '-', '-') and \
                    s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'].get_altloc() == \
                    s_mobile[0][resi[ind][1]][(resi[ind][2], resi[ind][0], resi[ind][3])][
                        'CA'].get_altloc():  # checks that they are coordinates rather than ['-']
                #                rmsd=rmsd+pow(np.linalg.norm(s_mobile[0][resi[k][1]][(resi[k][2],resi[k][0],resi[k][3])]['CA'].get_vector()-s_mobile[0][resi[ind][1]][(resi[ind][2],resi[ind][0],resi[ind][3])]['CA'].get_vector()),2)
                repeat1.append(
                    s_mobile[0][resi[k][1]][(resi[k][2], resi[k][0], resi[k][3])]['CA'].get_vector().get_array())
                repeat2.append(s_mobile[0][resi[ind][1]][(resi[ind][2], resi[ind][0], resi[ind][3])][
                                   'CA'].get_vector().get_array())
            #                n+=1
        tm_score = tm_score + tm_score_calc(repeat1, repeat2, prot_len)

    tm_score = tm_score * num_repeats / len(combo_list)

    return super_pml_files, tm_score


def pymol_superposition_script(inputf, rot_matrix, transl_vector, repeats, repeats_selection, unique_repeats, super_dir,
                               repeats_transf, image_key, level=0):
    # Generate PyMOL superposition script
    # repeats_all=[[['1PV6.A_6-94', '1PV6.A_219-309'],['1PV6.A_101-187', '1PV6.A_310-399']],[['1PV6.A_6-94', '1PV6.A_101-187'],['1PV6.A_219-309', '1PV6.A_310-399']]]
    super_pml_files = ""
    cb_colors, cb_names = colorblind_palette()
    pml_colors, pml_names = additional_colors_palette()
    colors = cb_names + pml_names
    out = open(super_dir + "superposition_" + image_key + "_level_" + str(level + 1) + ".pml", "w")
    super_pml_files = super_pml_files + "superposition_" + image_key + "_level_" + str(level + 1) + ".pml;"
    out.write("load " + inputf[0:4] + ".pdb, " + inputf + "\n")
    out.write("set antialias,1\nset line_smooth, 1\nset depth_cue, 0.5\nset specular, 0\nset surface_quality, 1\n")
    out.write("set cartoon_sampling, 14\nset ribbon_sampling, 10\nset ray_trace_fog, 0.5\n")
    out.write("viewport 350,350\n")
    out.write("set ignore_case, off\n")
    for c in cb_colors:
        out.write("set_color cb_" + c[0] + "=[" + str(c[1][0]) + "," + str(c[1][1]) + "," + str(c[1][2]) + "]\n")
    out.write("hide everything,*\n")
    out.write("show cartoon\n")
    out.write("bg_color white\nset opaque_background, on\n")
    out.write("import numpy as np\n")
    # Find the correct residues to select for the repeat that is connected with each particular symmetry axis;
    # Note that we are not selecting each CE-Symm defined repeat but rather the repeats that correspond to each axis
    # So in (2A65.A_9-235;2A65.B_9-235)(2A65.A_242-441;2A65.B_242-441)
    # Repeat 1 is 2A65.A_9-235 and 2A65.A_242-441
    # Repeat 2 is 2A65.B_9-235 and 2A65.B_242-441
    # repeats_all=[[['1PV6.A_6-94', '1PV6.A_219-309'],['1PV6.A_101-187', '1PV6.A_310-399']],[['1PV6.A_6-94', '1PV6.A_101-187'],['1PV6.A_219-309', '1PV6.A_310-399']]]

    for m in range(0, level + 1):
        out.write(
            "transl" + str(m) + "=np.array([" + str(transl_vector[m][0]) + "," + str(transl_vector[m][1]) + "," + str(
                transl_vector[m][2]) + "])\n")
        out.write("T" + str(m) + "=np.array([[" + str(rot_matrix[m][0][0]) + "," + str(rot_matrix[m][0][1]) + "," + str(
            rot_matrix[m][0][2]) + "],[" + str(rot_matrix[m][1][0]) + "," + str(rot_matrix[m][1][1]) + "," + str(
            rot_matrix[m][1][2]) + "],[" + str(rot_matrix[m][2][0]) + "," + str(rot_matrix[m][2][1]) + "," + str(
            rot_matrix[m][2][2]) + "]])\n")

    k = 0
    for j in range(0, len(repeats[0])):
        for s in range(0, len(repeats)):
            sel_text = ""
            entry = repeats[s][j]
            try:
                ch = entry.split('.')[1][0]
            except:
                ch = entry.split('_')[0]
            beg = entry.split('_')[1]
            first_digit = beg[0]  # in case the residues has a negative id (i.e. -1)
            beg = beg[1:].split('-')[0]
            beg = first_digit + beg
            locator = [i for i, sublist in enumerate(repeats_selection) if ((int(beg), ch) == sublist[0])]
            for i in range(0, int(len(repeats_selection[locator[0]]) / 2)):
                sel_text = sel_text + "(chain " + repeats_selection[locator[0]][i * 2][1] + " and resid " + str(
                    repeats_selection[locator[0]][i * 2][0]) + "-" + str(
                    repeats_selection[locator[0]][i * 2 + 1][0]) + ")"
                if i < len(repeats_selection[locator[0]]) / 2 - 1:
                    sel_text = sel_text + " or "
            # if s<(len(repeats)-1):
            #    sel_text=sel_text+" or "
            # u=unique_repeats.index(entry)
            try:
                u = unique_repeats.index(entry)
            except:
                u = unique_repeats.index(inputf[0:4].upper() + "." + entry)
            # u = locator[0]
            out.write("create repeat" + str(u + 1) + ", " + inputf + " and " + sel_text + "\n")
            out.write("color " + colors[(k % len(colors))] + ",repeat" + str(u + 1) + "\n")
            for m in range(0, level + 1):
                for p in range(0, repeats_transf[m][u]):
                    out.write("alter_state 1, repeat" + str(u + 1) + ", (x,y,z)=np.dot((np.array([x,y,z])-transl" + str(
                        m) + "),T" + str(m) + ")\n")

        k += 1
    out.write("deselect\n")
    out.write("rotate [1,0,0], angle=-90\n")
    out.write("zoom " + inputf + "\n")
    out.write("delete " + inputf + "\n")
    #        out.write("cd "+super_dir+"\n")
    #        out.write("png "+image_key+"_super_"+str(ind)+".png\n")
    #        out.write("quit")
    out.close()
    #        os.system(pymol+" "+super_dir+"level_"+str(ind)+"_superposition_"+image_key+".pml")
    return super_pml_files


def create_tm_archive(locations, wkdir):
    print("[INFO]: Creating a slimmer TM archive...")
    opm_archive_path = sorted(glob.glob(locations['FSYSPATH']['main'] + '.opm_archive_*.pkl'))[-1]
    print("[INFO]: Full archive to be used: ", opm_archive_path)
    opm_archive = pkl.load(open(opm_archive_path, 'rb'))
    tm_archive = {}
    for pdb in sorted(list(opm_archive)):
        print(pdb)
        tm_archive[pdb] = {'class': '', 'tmchains': {}}
        tm_archive[pdb]['class'] = opm_archive[pdb]['class']
        tm_archive[pdb]['tmchains'] = opm_archive[pdb]['tmchains']
        for chain in opm_archive[pdb]['tmchains']:
            tm_archive[pdb][chain] = opm_archive[pdb][chain]
            pdb_path = glob.glob(locations['FSYSPATH']['whole'] + pdb + "*.pdb")
            if len(pdb_path) > 0:
                num = strip_tm_chains(wkdir, pdb, pdb_path[0], chain)
                os.remove(wkdir + pdb + "_tmp.pdb")
                tm_archive[pdb][chain]['natoms'] = num
            else:
                tm_archive[pdb][chain]['natoms'] = 0

    pkl.dump(tm_archive, open(locations['FSYSPATH']['symmetry'] + '.tm_archive.pkl', 'wb'))
    print("[INFO]: TM archive generated - " + locations['FSYSPATH']['symmetry'] + '.tm_archive.pkl')
    return


def create_tm_archive_new(locations, wkdir, arch_name):
    currentdir = os.path.dirname(os.path.realpath(__file__))
    parentdir = os.path.dirname(currentdir)
    sys.path.append(parentdir)
    from supporting_functions import read_checkpoint

    print("[INFO]: Creating a slimmer TM archive...")
    encompass_archive_path = locations['FSYSPATH']['cache'] + arch_name
    print("[INFO]: Full archive to be used: ", encompass_archive_path)
    encompass_archive = read_checkpoint(encompass_archive_path, locations['SYSFILES']['data_structure_template'])
    tm_archive = {}
    for pdb in sorted(list(encompass_archive)):
        print(pdb)
        tm_archive[pdb] = {'class': '', 'tmchains': [], 'natoms': 0, 'name' : ''}
        tm_archive[pdb]['class'] = encompass_archive[pdb]['ENCOMPASS']['class']
        tm_archive[pdb]['name'] = encompass_archive[pdb]['ENCOMPASS']['name']
        pdb_path = glob.glob(locations['FSYSPATH']['whole'] + pdb + "*.pdb")
        tm_archive[pdb]['tmchains'] = [ch for ch in encompass_archive[pdb]['ENCOMPASS']['structure']['ktmchains'] if
                                       ch != '-']
        if len(pdb_path) > 0:
            num = strip_tm_chains(wkdir, pdb, pdb_path[0], tm_archive[pdb]['tmchains'])
            res_num = resnum(wkdir + pdb + "_tmp.pdb")
            os.remove(wkdir + pdb + "_tmp.pdb")
            tm_archive[pdb]['natoms'] = num
            tm_archive[pdb]['nres'] = res_num
            tm_archive[pdb]['unk'] = False
        for ichain, chain in [(ich, ch) for ich, ch in
                              enumerate(encompass_archive[pdb]['ENCOMPASS']['structure']['ktmchains']) if ch != '-']:
            tm_archive[pdb][chain] = encompass_archive[pdb]['ENCOMPASS']['structure']['chains'][ichain].show_dictionary(
                quiet=True, deeptransparency=True)
            tm_archive[pdb][chain]['natoms'] = 0
            tm_archive[pdb][chain]['nres'] = 0
            tm_archive[pdb][chain]['unk'] = False
            if len(pdb_path) > 0:
                num = strip_tm_chains(wkdir, pdb, pdb_path[0], chain)
                res_num = resnum(wkdir + pdb + "_tmp.pdb")
                os.remove(wkdir + pdb + "_tmp.pdb")
                tm_archive[pdb][chain]['natoms'] = num
                tm_archive[pdb][chain]['nres'] = res_num
                if 'X' in tm_archive[pdb][chain]['sequence']:
                    tm_archive[pdb][chain]['unk'] = True
                    tm_archive[pdb]['unk'] = True

    pkl.dump(tm_archive, open(locations['FSYSPATH']['symmetry'] + '.tm_archive.pkl', 'wb'))
    print("[INFO]: TM archive generated - " + locations['FSYSPATH']['symmetry'] + '.tm_archive.pkl')
    return


# Completion check
def check_completion(locations, pdb_list, proc_name):
    failed = set()
    completed = set()
    error_log = open('error_' + proc_name + '.txt', 'w')
    error_log.write("#PDB_code\n")
    for pdb in pdb_list:
        if os.path.isdir(locations['FSYSPATH'][proc_name + '_log'] + pdb):
            logs = glob.glob(locations['FSYSPATH'][proc_name + '_log'] + pdb + "/" + pdb + "*.log")
        else:
            logs = glob.glob(locations['FSYSPATH'][proc_name + '_log'] + pdb + "*.log")
        if len(logs) == 1:
            f = open(logs[0], 'r')
        elif len(logs) > 1:
            print('Warning: multiple possible log files to check completion for %s' % pdb)
            similarity = [abs(len(s.split('/')[-1]) - len(pdb)) for s in logs]
            ind = similarity.index(min(similarity))
            f = open(logs[ind], 'r')
            print('Using log file %s' % logs[ind])
        else:
            continue
        flag = 0
        for line in f:
            if "Error" in line or "Traceback" in line:
                flag = 1
                error_log.write(pdb + "\t" + line.strip() + "\n")
                failed.add(pdb)
                break
            if line.startswith(pdb) and (" completed." in line) and flag == 0:
                completed.add(pdb)
            if (".pdb does not exist." in line) and flag == 0:
                completed.add(pdb)
        f.close()

    return completed, failed
