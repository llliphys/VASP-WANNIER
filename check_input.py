#!/home/li/ANACONDA/bin/python

import os
import sys
import numpy as np

def check_incar(args):

    print("===============Checking INCAR==============")

    if not os.path.isfile("INCAR"):
        print("Error: No INCAR Found !")
        exit(1)

    with open("INCAR", "r") as f:
        incar = f.readlines()
    for i in range(len(incar)):
        incar[i] = incar[i].strip().strip("\n")

    for line in incar:
        if not args:
            print(line)
        else:
            for arg in args:
                if arg.upper() in line:
                    print(line)


def check_poscar(args):

    print("===============Checking POSCAR==============")

    if not os.path.isfile("POSCAR"):
        print("Error: No POSCAR Found !")
        eixt(1)

    with open("POSCAR", "r") as f:
        poscar = f.readlines()
    for i in range(len(poscar)):
        poscar[i] = poscar[i].strip().strip("\n")

    if "lattice" in args:
        for line in poscar[0:7]: print(line)
    if "position" in args:
        for line in poscar[7::]: print(line)
    if not args:
        for line in poscar[0::]: print(line)

def check_contcar(args):

    print("===============Checking CONTCAR==============")

    if not os.path.isfile("CONTCAR"):
        print("Error: No CONTCAR Found !")
        eixt(1)

    with open("CONTCAR", "r") as f:
        poscar = f.readlines()
    for i in range(len(poscar)):
        poscar[i] = poscar[i].strip().strip("\n")

    if "lattice" in args:
        for line in poscar[0:7]: print(line)
    if "position" in args:
        for line in poscar[7::]: print(line)
    if not args:
        for line in poscar[0::]: print(line)

def check_potcar(args):

    print("===============Checking POTCAR==============")

    if not os.path.isfile("POTCAR"):
        print("Error: No POTCAR Found !")
        exit(1)

    with open("POTCAR", "r") as f:
        potcar = f.readlines()
    for i in range(len(potcar)):
        potcar[i] = potcar[i].strip().strip("\n")

    for line in potcar:
        if "titel" in args and "TITEL" in line: print(line)
        if "enmax" in args and "ENMAX" in line: print(line)
        if not args and "TITEL" in line: print(line)
        if not args and "ENMAX" in line: print(line)

def check_kpoints(args):

    print("===============Checking KPOINTS==============")

    if not os.path.isfile("KPOINTS"):
        print("Error: No KPOINTS Found !")
        exit(1)

    with open("KPOINTS", "r") as f:
        kpoints = f.readlines()
    for i in range(len(kpoints)):
        kpoints[i] = kpoints[i].strip().strip("\n")
    for line in kpoints:
        if not args: print(line)

def main():

    args = sys.argv[1:] # a list of args

    if args == []:
        print("Check VASP Inputs with Arguments below")
        print("Mandatory Arguments:")
        print("incar     --->    for INCAR checking")
        print("poscar    --->    for POSCAR checking")
        print("potcar    --->    for POTCAR checking")
        print("kpoints   --->    for KPOINTS checking")
        print("contcar   --->    for CONTCAR checking")
        print("Optional Arguments:")
        print("nsw/...   --->    for INCAR checking with given tag information printed out")
        print("lattice   --->    for POSCAR checking with lattice information printed out")
        print("position  --->    for POSCAR checking with atomic coordinates printed out")
        print("encutmax  --->    for POTCAR checking with ENMAX information printed out")
        exit(1)

    if "incar" in args: check_incar(args[1::])
    if "poscar" in args: check_poscar(args[1::])
    if "potcar" in args: check_potcar(args[1::])
    if "kpoints" in args: check_kpoints(args[1::])
    if "contcar" in args: check_kpoints(args[1::])

if __name__ == "__main__":
   main()
