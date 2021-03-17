#!/home/li/ANACONDA/bin/python

import os
import sys
# import argparse
import numpy as np


def check_outcar(args):

    print("===============Checking OUTCAR==============")

    if not os.path.isfile("OUTCAR"):
        print("Error: No OUTCAR Found !")
        exit(1)

    with open("OUTCAR", "r") as f:
        outcar = f.readlines()
    for i in range(len(outcar)):
        outcar[i] = outcar[i].strip().strip("\n")

    lookup_rlx = "reached required accuracy"
    lookup_scf = "aborting loop because EDIFF is reached"

    rlx_conv = None
    for line in outcar[::-1]:
        if "reached required accuracy" in line:
            rlx_conv = True
            break

    scf_conv = None
    for line in outcar[::-1]:
        if "aborting loop because EDIFF is reached" in line:
            scf_conv = True
            break

    if rlx_conv is not None:
        if rlx_conv:
            print("Relaxation Converged !")
        else:
            print("Relaxation NOT Converged !")
            exit(1)

    if scf_conv is not None:
        if scf_conv:
            print("Electronic SCF Converged !")
        else:
            print("Electronic SCF NOT Converged !")
            exit(1)

    if "energy" in args: get_energy(outcar)
    if "magmom" in args: get_magmom(outcar)


def get_energy(outcar):

    for i, line in enumerate(outcar[::-1]):
        if "E-fermi" in line:
            efermi = float(line.split(":")[1].split()[0])
            print("Efermi = %.6f [eV]" % efermi)
            break

    for i, line in enumerate(outcar[::-1]):
        if "TOTEN" in line:
            etotal = float(line.split("=")[1].split()[0])
            print("Etotal = %.6f [eV]" % etotal)
            break

def get_magmom(outcar):

    for i, line in enumerate(outcar[::-1]):
        if "number of ions" in line:
            natoms = int(line.strip().split("NIONS =")[1])
            print("Number of atoms in the cell: ", natoms)
            break

    mspnx = mspny = mspnz = None
    for i, line in enumerate(outcar):
        if "magnetization (x)" in line:
            mspnx = outcar[i+4:i+4+natoms]
        if "magnetization (y)" in line:
            mspny = outcar[i+4:i+4+natoms]
        if "magnetization (z)" in line:
            mspnz = outcar[i+4:i+4+natoms]

    mspnx_dict = {}
    mspnx_list = []
    if mspnx is not None:
        for i, line in enumerate(mspnx):
            line = line.strip().split()
            mx = round(float(line[-1]), 3)
            mspnx_list.append(mx)

    mspny_dict = {}
    mspny_list = []
    if mspny is not None:
        for i, line in enumerate(mspny):
            line = line.strip().split()
            my = round(float(line[-1]), 3)
            mspny_list.append(my)

    mspnz_dict = {}
    mspnz_list = []
    if mspnz is not None:
        for i, line in enumerate(mspnz):
            line = line.strip().split()
            mz = round(float(line[-1]), 3)
            mspnz_list.append(mz)

    if mspnx is not None: print("Magnetization (x) = ", mspnx_list)
    if mspny is not None: print("Magnetization (y) = ", mspny_list)
    if mspnz is not None: print("Magnetization (z) = ", mspnz_list)


def main():

    args = sys.argv[1:] # a list of args

    if args == []:
        print("Check VASP Outputs with Arguments below")
        print("Mandatory Arguments:")
        print("outcar    --->    for OUTCAR checking to see if a relax or scf is converged")
        print("Optional Arguments:")
        print("energy    --->    for OUTCAR checking with energy information printed out")
        print("magmom    --->    for OUTCAR checking with magmom information printed out")
        print("How to Use:")
        print("python check_output.py outcar [energy] [magmom]")
        exit(1)

    if "outcar" in args: check_outcar(args[1::])

if __name__ == "__main__":
   main()
