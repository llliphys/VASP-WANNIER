#!/home/li/ANACONDA/bin/python

import os
import sys
import argparse
import numpy as np
# import ase.io.vasp
# from ase import io
import ase

def modify_cell(infile=None, direct=None, vacuum=None, numcell=None):

    print("===============Modifying POSCAR==============")

    if infile is None:
        infile = "POSCAR"

    print(infile)

    if not os.path.isfile(infile):
        print("Specify a CIF or POSCAR file")
        exit(1)

    # if direct is None:
    #     direct = False
    # else:
    #     direct = True
    #
    # if vacuum is None:
    #     vacuum = vacuum
    # else:
    #     vacuum = 20

    # if numcell is None:
    #     numcell = (1, 1, 1)
    # else:
    #     numcell = numcell

    atoms = ase.io.read(infile)
    atoms.center(axis=2, vacuum=vacuum/2.0)
    # atoms.write("POSCAR-UC", format="vasp", vasp5=True, direct=direct)

    product = numcell[0] * numcell[1] * numcell[2]

    if product == 1:
        ase.io.vasp.write_vasp("POSCAR-UC", atoms*numcell, label="Supercell", vasp5=True, direct=direct, sort=True)
    else:
        ase.io.vasp.write_vasp("POSCAR-SC", atoms*numcell, label="Supercell", vasp5=True, direct=direct, sort=True)

def main():

    args = sys.argv[1:]

    arg = argparse.ArgumentParser() #(add_help=True)
    arg.add_argument('-i', type=str, default=None, help='POSCAR or CIF file (type=str)')
    arg.add_argument('-x', type=int, default=1, help='Supercell size along x (type=int)')
    arg.add_argument('-y', type=int, default=1, help='Supercell size along y (type=int)')
    arg.add_argument('-z', type=int, default=1, help='Supercell size along z (type=int)')
    arg.add_argument('-d', type=str, default='true', help='Direct or Cartesian (type=str)')
    arg.add_argument('-v', type=float, default=20, help='Vacuum length (type=float)')

    args = arg.parse_args(args)

    modify_cell(infile=args.i, numcell=(args.x, args.y, args.z), direct=args.d=="true", vacuum=args.v)

if __name__ == "__main__":
   main()
