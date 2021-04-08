#!/home/li/ANACONDA/bin/python


import os
import re
import sys
import math
import numpy as np
import scipy as sp
from itertools import chain
import matplotlib.pyplot as plt

# try:
#     from pandas import read_table
#     use_pandas = True
# except ImportError:
#     use_pandas = False

class WorkFunction:
    """
    A python module to extract from a VASP calculation the fermi energy, the planar (averaged)
    potential, the vacuum potential, and the work function of a slab system with thick vacuum.
    To this end, one has to provide two VASP output files, i.e., OUTCAR and LOCPOT (LVHAR has
    to be set .TRUE. in the INCAR file to obtain the LOCPOT file) .

    Some part of this module is taken from https://github.com/WMD-group/MacroDensity but with
    substantial modification for easy-to-use purposes. Pls feel free to accomodate your owns.

    Example usage is given in the end of this module, e.g. calculating thework function and
    plotting the planar averaged potential along a specified direction (axis = "x"/"y"/"z").
    """
    def __init__(self, wrkdir=None):

        self.wrkdir = wrkdir

    def GetFermilevel(self, outcar="OUTCAR"):

        outcar = open(outcar, "r").read()
        match = re.findall(r"E-fermi\s*:\s*(-?\d+.\d+)", outcar)[-1]
        efermi = float(match)
        self.efermi = efermi


    def GetPotential(self, locpot="LOCPOT"):

        with open(locpot, "r") as f:
            fread = f.read()
            _ = f.readline()
            scale_factor = float(f.readline())

            lattice = np.zeros(shape=(3, 3))
            for row in range(3):
                lattice[row] = [float(x) for x in f.readline().split()]
            lattice = lattice * scale_factor

            num_of_species = len(f.readline().split())
            num_of_atoms_per_specie = [int(x) for x in f.readline().split()]
            num_of_atoms = sum(num_of_atoms_per_specie)

            _ = f.readline()

            coordinate = np.zeros(shape=(num_of_atoms, 3))
            for atom_i in range(num_of_atoms):
                coordinate[atom_i] = [float(x) for x in f.readline().split()]

            # Convert to Cartesian Coordinates
            for i in range(coordinate.shape[0]):
                v = coordinate[i]
                v = lattice.T.dot(v)
                coordinate[i] = v

            _ = f.readline()

            nx, ny, nz = [int(x) for x in f.readline().split()]

            # regex = r"\s*" + r"\s*".join(["(\d+\.\d+E[+-]\d+)" for i in range(5)])
            # pot = re.findall(regex, fread)
            # pot = [[float(v) for v in vs] for vs in data]
            # pot = np.array(pot).flatten()

            # Extract 1D-flattened local potential
            pot = (f.readline().split() for i in range(int(math.ceil(nx * ny * nz / 5))))
            pot = np.fromiter(chain.from_iterable(pot), float)

            # Rearrange the 1D-flattened potential into a (nx, ny, nz)-shaped array
            potential = np.zeros((nx, ny, nz))
            l = 0
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        potential[i, j, k] = pot[l]
                        l = l + 1

            self.coordinate = coordinate
            self.meshgrid = (nx, ny, nz)
            self.potential = potential


    def PlanarAverage(self, axis="z"):

        pot = self.potential
        nx, ny, nz = self.meshgrid
        coord = self.coordinate

        if axis == "x":
            x = coord[:, 0]
            position = np.linspace(x.min(), x.max(), nx)
            pot_x_plane = np.zeros(shape=(ny, nz))
            plane_average = np.zeros(shape=(nx))
            for ix in range(nx):
                pot_x_plane[:, :] = pot[ix, :, :]
                plane_average[ix] = pot_x_plane.mean()

        if axis == "y":
            y = coord[:, 1]
            position = np.linspace(y.min(), y.max(), ny)
            pot_y_plane = np.zeros(shape=(nx, nz))
            plane_average = np.zeros(shape=(ny))
            for iy in range(ny):
                pot_y_plane[:, :] = pot[:, iy, :]
                plane_average[iy] = pot_y_plane.mean()

        if axis == "z":
            z = coord[:, 2]
            position = np.linspace(z.min(), z.max(), nz)
            pot_z_plane = np.zeros(shape=(nx, ny))
            plane_average = np.zeros(shape=(nz))
            for iz in range(nz):
                pot_z_plane[:, :] = pot[:, :, iz]
                plane_average[iz] = pot_z_plane.mean()

        self.coords = position
        self.planar = plane_average
        self.enevac = plane_average.max()

    def PlotPotential(self, ax=None, grid=True, show=True):

        if ax is None: fig, ax = plt.subplots(1, 1, sharex=True)

        ax.plot(self.coords, self.planar, lw=3)
        ax.set_xlim(self.coords.min(), self.coords.max())
        ax.set_facecolor((0.95,0.95,0.95))
        if grid: ax.grid(True)
        ax.set_xlabel("$\mathrm{Position} \ (\mathrm{\AA})$", fontsize=30)
        ax.set_ylabel("$\mathrm{Potential} \ (\mathrm{eV})$", fontsize=30)
        if show: plt.show()


def main():

    c = WorkFunction()
    c.GetFermilevel(outcar="OUTCAR") # path for OUTCAR
    c.GetPotential(locpot="LOCPOT")  # path for LOCPOT
    c.PlanarAverage(axis="z") # potential along z direction
    print("Work Function = %.3f eV" % (c.enevac - c.efermi))
    # c.PlotPotential(show=True) # Plot e.g. Planar Potential


if __name__ == "__main__":

    main()
