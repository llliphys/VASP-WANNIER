import os
import re
import argparse
import numpy as np
import argparse

# check pymatgen installed or not
try:
    import pymatgen
    pymatgen_install = True
except ImportError:
    pymatgen_install = False

# check ase installed or not
try:
    import ase
    ase_install = True
except ImportError:
    ase_install = False

def GetPoscar():
    """
    ToDo: Use ASE to get atom symbols, indices and positions:  
    if ase_install:
        from ase import io
    else:
        print("Warning: ASE is not yet installed on your machine !")
        exit(1)
    """

    if not os.path.exists("POSCAR"):
        print("No POSCAR Exist")
        exit(1)

    with open("POSCAR", "r") as f:
        poscar = f.readlines()

    avec = poscar[2]
    bvec = poscar[3]
    cvec = poscar[4]

    ax, ay, az = [float(s) for s in avec.strip().split()]
    bx, by, bz = [float(s) for s in bvec.strip().split()]
    cx, cy, cz = [float(s) for s in cvec.strip().split()]

    avec = np.array([ax, ay, az])
    bvec = np.array([bx, by, bz])
    cvec = np.array([cx, cy, cz])

    atoms_lst = poscar[5].strip().split()
    natom_str = poscar[6].strip().split()
    natom_arr = np.array([int(i) for i in natom_str])
    natoms = np.sum(natom_arr)

    atoms_all = []
    idexs_all = []
    idx = -1
    for i in range(len(atoms_lst)):
        atom_i = atoms_lst[i]
        for j in range(natom_arr[i]):
            atoms_all.append(atom_i)
            idx += 1
            idexs_all.append(idx)
    atoms_all = np.array(atoms_all)
    idexs_all = np.array(idexs_all)

    cords_all = []
    for i in np.arange(8, 8+natoms):
        line = poscar[i]
        x, y, z = [float(s) for s in line.strip().split()]
        # convert direct to cartesian coordinates
        Rvec = avec * x + bvec * y + cvec * z
        x, y, z = Rvec[0], Rvec[1], Rvec[2]
        cords_all.append([x, y, z])
    cords_all = np.vstack(cords_all)

    return idexs_all, atoms_all, cords_all



def MakeClean():

    inputs = ["INCAR", "POSCAR", "POTCAR", "KPOINTS"]

    for f in os.listdir("."):
        if f not in inputs: 
            os.system("rm -rf %s" % f)


def GetEfermi():

    outcar = open("OUTCAR", "r").read()
    match = re.findall(r"E-fermi\s*:\s*(-?\d+.\d+)", outcar)[-1]
    Efermi = float(match)

    return Efermi


def GetNbands():

    outcar = open("OUTCAR", "r").read()
    match = re.findall(r"NBANDS\s*=\s*(\d+)", outcar)[-1]
    Nbands = int(match)

    return Nbands


def MakeKpath(bands_kpath=None, bands_nkpts=None):

    if pymatgen_install:
        from pymatgen.io.vasp.inputs import Kpoints
        from pymatgen.core import Structure
        from pymatgen.symmetry.bandstructure import HighSymmKpath
    else:
        print("Warning: PYMATGEN is not yet installed on your machine !")
        exit(1)

    if not os.path.isfile("POSCAR"):
        print("No POSCAR present for making a kpath for Wannier band !")
        exit(1)

    if bands_nkpts is None: 
        bands_nkpts = 100

    kpath = bands_kpath
    nkpts = bands_nkpts

    if kpath is not None:
        kpath = [[str(s) for s in kpath]]
        ibz = HighSymmKpath(Structure.from_file("POSCAR"))
        ibz.kpath['path'] = kpath
        Kpoints.automatic_linemode(nkpts, ibz).write_file("KPATH")
    else:
        ibz = HighSymmKpath(Structure.from_file("POSCAR"))
        Kpoints.automatic_linemode(nkpts, ibz).write_file("KPATH")

    with open("KPATH", "r") as f:
        klines = f.readlines()

    regex = r"\s*(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*\!\s*\\*([A-Z]*)[a-z]*"

    file = open("kpath.win", "w")
    file.write("\n")
    file.write("begin kpoint_path \n")
    i = 0
    for line in klines:
        find = re.match(regex, line)
        if find != None:
            i = i + 1
            kx, ky, kz, s = find.groups()
            kx, ky, kz = float(kx), float(ky), float(kz)
            if i % 2 != 0: file.write("%s  %.9f  %.9f  %.9f  " % (s, kx, ky, kz))
            if i % 2 == 0: file.write("%s  %.9f  %.9f  %.9f \n" % (s, kx, ky, kz))
    file.write("end kpoint_path \n")
    file.write("\n")
    file.close()
        
        
def GrepText(text, lines):

    regex = r"\s*(%s)\s*" % text

    has_no_text = True
    for line in lines:
        if re.match(regex, line): #text in line:
            has_no_text = False

    return has_no_text

    # Alternative: Find All Matched Data Info
    # kpoints = open("KPATH", "r").read()
    # match = re.findall(regex, kpoints)
    # type(match) = a list of tuples


class PythWann:

    """
    - A Python module for carrying out a series of Wannier calculations including e.g.
    "Vasp2Wannier", "WannierFitting", "BerryCurvature", "AnomalousHallConductivity",
    "SurfaceStates", which are reprensted here by task="WAN", "FIT", "BC", "AHC", "SS",
    respectively.

    - Before doing any calculations, one has always to specify two working directories,
    one of which contains the input files for Wannier calculations and the other contains
    the output files from the calculations. The two directories for inputs and outputs are
    reprensted by scrdir and objdir, respectively. Note that scrdir and objdir can be the
    same or different directories, depending on specific calculations.

    - A rootdir can also be specified to indicate the root directory that contains the two
    subdirectories inside it, i.e., root/srcdir and root/objdir.

    """
    
    def __init__(self, task="WAN", topdir=None, wrkdir=None, srcdir=None, objdir=None):

        os.chdir(topdir)
        
        if wrkdir is None:
            wrkdir = os.getcwd()
        if srcdir is None:
            print("No INPUT Directory for WANNIER Calculation !")
            exit(1)
        if objdir is None:
            print("No OUTPUT Directory for WANNIER Calculation !")
            exit(1)    
                    
        self.task = task
        self.wrkdir = os.path.abspath(wrkdir)
        self.srcdir = os.path.join(self.wrkdir, srcdir)
        self.objdir = os.path.join(self.wrkdir, objdir)

        if not os.path.exists(self.srcdir): 
            os.mkdir(self.srcdir)
        
        if not os.path.exists(self.objdir): 
            os.mkdir(self.objdir)

    def MakeInput(self, 
                  icharg=None, istart=None, nkpts=None, gamma=True, \
                  nbands=None, npar=None, kpar=None, \
                  lsoc=False, lgct=True, kmesh_tol=None, search_shells=None, \
                  proj=None, num_wann=None, num_bands=None, exclude_bands=None, \
                  dis_win_min=None, dis_win_max=None, dis_froz_min=None, dis_froz_max=None, \
                  dis_num_iter=None, num_iter=None, \
                  write_hr=False, write_xyz=False, \
                  bands_plot=False, bands_num_points=None, bands_kpath=None, \
                  kpath=False, kpath_task=None, kpath_bands_colour=None, kpath_num_points=None, \
                  fermi_surface_plot=False, fermi_surface_num_points=None, \
                  fermi_energy=None, fermi_shift=None, \
                  kslice=False, kslice_task=None, kslice_2dkmesh=None, \
                  kslice_b1=None, kslice_b2=None, kslice_corner=None, \
                  berry=False, berry_task=None, berry_kmesh=None, berry_curv_unit=None, \
                  fermi_energy_min=None, fermi_energy_max=None, fermi_energy_step=None, \
                  kubo_freq_min=None, kubo_freq_max=None, kubo_freq_step=None, \
                  kubo_smr_type=None,  kubo_smr_fixed_en_width=None,\
                  make_clean=True, \
                  seed="wannier90"):

        self.lsoc = lsoc
        self.seed = seed
                
        win= dict()

        if lsoc: win.update(spinors="true")
        if lgct: win.update(guiding_centres="true")
        # if proj is not None: win.update(projections=proj)
        if num_wann is not None: win.update(num_wann=num_wann)
        if nbands is not None: win.update(num_bands=nbands)
        if exclude_bands is not None:
            win.update(exclude_bands=exclude_bands)
            if "," in exclude_bands:
                nbands_exclude = 0
                for s1 in exclude_bands.split(","):
                    s2 = s1.split("-")
                    nbands_exclude += int(s2[1]) - int(s2[0]) + 1
            else:
                s2 = exclude_bands.split("-")
                nbands_exclude = int(s2[1]) - int(s2[0]) + 1
            if nbands is not None: win.update(num_bands=nbands-nbands_exclude)
        if kmesh_tol is not None: win.update(kmesh_tol=kmesh_tol)
        if search_shells is not None: win.update(search_shells=search_shells)
        if dis_win_min is not None: win.update(dis_win_min=dis_win_min)
        if dis_win_max is not None: win.update(dis_win_max=dis_win_max)
        if dis_froz_min is not None: win.update(dis_froz_min=dis_froz_min)
        if dis_froz_max is not None: win.update(dis_froz_max=dis_froz_max)
        if num_iter is not None: win.update(num_iter=num_iter)
        if dis_num_iter is not None: win.update(dis_num_iter=dis_num_iter)
        if bands_plot: win.update(bands_plot="true")
        if bands_num_points is not None: win.update(bands_num_points=bands_num_points)
        if write_hr: win.update(write_hr="true")
        if write_xyz: win.update(write_xyz="true")
        # if restart is not None: win.update(restart=restart)
        if fermi_surface_plot: win.update(fermi_surface_plot="true")
        if fermi_surface_num_points is not None: win.update(fermi_surface_num_points=fermi_surface_num_points)
        # if fermi_energy is not None: win.update(fermi_energy=fermi_energy)
        if kpath: win.update(kpath="true")
        if berry_curv_unit is not None: win.update(berry_curv_unit=berry_curv_unit)
        if kpath_task is not None: win.update(kpath_task=kpath_task)
        if kpath_bands_colour is not None: win.update(kpath_bands_colour=kpath_bands_colour)
        if kpath_num_points is not None: win.update(kpath_num_points=kpath_num_points)
        if kslice: win.update(kslice="true")
        if kslice_task is not None: win.update(kslice_task=kslice_task)
        if kslice_2dkmesh is not None: win.update(kslice_2dkmesh=kslice_2dkmesh)
        if kslice_corner is not None: win.update(kslice_corner=kslice_corner)
        if kslice_b1 is not None: win.update(kslice_b1=kslice_b1)
        if kslice_b2 is not None: win.update(kslice_b2=kslice_b2)
        if berry: win.update(berry="true")
        if berry_task is not None: win.update(berry_task=berry_task)
        if berry_kmesh is not None: win.update(berry_kmesh=berry_kmesh)
        if fermi_energy_min is not None: win.update(fermi_energy_min=fermi_energy_min)
        if fermi_energy_max is not None: win.update(fermi_energy_max=fermi_energy_max)
        if fermi_energy_step is not None: win.update(fermi_energy_step=fermi_energy_step)  
        if kubo_freq_min is not None: win.update(kubo_freq_min=kubo_freq_min)   
        if kubo_freq_max is not None: win.update(kubo_freq_max=kubo_freq_max)   
        if kubo_freq_step is not None: win.update(kubo_freq_step=kubo_freq_step)   
        if kubo_smr_type is not None: win.update(kubo_smr_type=kubo_smr_type) 
        if kubo_smr_fixed_en_width is not None: win.update(kubo_smr_fixed_en_width=kubo_smr_fixed_en_width)
                
        os.chdir(self.objdir)
        if make_clean: MakeClean()
        os.chdir(self.wrkdir)
            
        if self.task.upper() == "WAN":

            os.chdir(self.srcdir)
            if not os.path.isfile("CHGCAR") and not os.path.isfile("WAVECAR"):
                raise SystemExit("No CHGCAR or WAVECAR for VASP2WANN Calculation !")
            os.chdir(self.wrkdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/INCAR %s/INCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/POSCAR %s/POSCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/POTCAR %s/POTCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/KPOINTS %s/KPOINTS" % (self.srcdir, self.objdir))
                if istart is not None:
                    os.system("cp %s/WAVECAR %s/WAVECAR" % (self.srcdir, self.objdir))
                if icharg is not None: 
                    os.system("cp %s/CHGCAR %s/CHGCAR" % (self.srcdir, self.objdir))            
            
        if self.task.upper() == "CHK":

            os.chdir(self.srcdir)
            if not os.path.isfile("%s.win"%seed):
                raise SystemExit("No INPUTS (.win) for Wannier Checking Calculation !")
            if not os.path.isfile("%s.amn"%seed):
                raise SystemExit("No INPUTS (.amn) for Wannier Checking Calculation !")
            if not os.path.isfile("%s.mmn"%seed):
                raise SystemExit("No INPUTS (.mmn) for Wannier Checking Calculation !")
            if not os.path.isfile("%s.eig"%seed):
                raise SystemExit("No INPUTS (.eig) for Wannier Checking Calculation !")
            os.chdir(self.wrkdir)

            os.chdir(self.srcdir)
            if fermi_energy is None: fermi_energy = GetEfermi()
            os.system("rm -rf EFERMI; echo %.3f >> EFERMI" % fermi_energy)
            os.chdir(self.wrkdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/POSCAR %s/POSCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/EFERMI %s/EFERMI" % (self.srcdir, self.objdir))
                os.system("cp %s/%s.win %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.amn %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.mmn %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.eig %s/" % (self.srcdir, seed, self.objdir))
            
            
        if self.task.upper() == "FIT":

            os.chdir(self.srcdir)
            if not os.path.isfile("%s.win"%seed):
                raise SystemExit("No INPUTS (.win) for Wannier Fitting Calculation !")
            if not os.path.isfile("%s.amn"%seed):
                raise SystemExit("No INPUTS (.amn) for Wannier Fitting Calculation !")
            if not os.path.isfile("%s.mmn"%seed):
                raise SystemExit("No INPUTS (.mmn) for Wannier Fitting Calculation !")
            if not os.path.isfile("%s.eig"%seed):
                raise SystemExit("No INPUTS (.eig) for Wannier Fitting Calculation !")
            os.chdir(self.wrkdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/POSCAR %s/POSCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/EFERMI %s/EFERMI" % (self.srcdir, self.objdir))
                os.system("cp %s/%s.win %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.amn %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.mmn %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.eig %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.chk %s/" % (self.srcdir, seed, self.objdir))

            os.chdir(self.objdir)
            MakeKpath(bands_kpath=bands_kpath)
            os.chdir(self.wrkdir)


        if self.task.upper() == "FSF":

            os.chdir(self.srcdir)
            if not os.path.isfile("%s.win"%seed):
                raise SystemExit("No INPUTS (.win) for Fermi Surface Calculation !")
            if not os.path.isfile("%s.amn"%seed):
                raise SystemExit("No INPUTS (.amn) for Fermi Surface Calculation !")
            if not os.path.isfile("%s.mmn"%seed):
                raise SystemExit("No INPUTS (.mmn) for Fermi Surface Calculation !")
            if not os.path.isfile("%s.eig"%seed):
                raise SystemExit("No INPUTS (.eig) for Fermi Surface Calculation !")
            if not os.path.isfile("%s.chk"%seed):
                raise SystemExit("No INPUTS (.chk) for Fermi Surface Calculation !")
            os.chdir(self.wrkdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/EFERMI %s/EFERMI" % (self.srcdir, self.objdir))
                os.system("cp %s/%s.win %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.amn %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.mmn %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.eig %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.chk %s/" % (self.srcdir, seed, self.objdir))


        if self.task.upper() == "BC1":

            os.chdir(self.srcdir)
            if not os.path.isfile("%s.win"%seed):
                raise SystemExit("No INPUTS (.win) for 1D Berry Curvature Calculation !")
            if not os.path.isfile("%s.amn"%seed):
                raise SystemExit("No INPUTS (.amn) for 1D Berry Curvature Calculation !")
            if not os.path.isfile("%s.mmn"%seed):
                raise SystemExit("No INPUTS (.mmn) for 1D Berry Curvature Calculation !")
            if not os.path.isfile("%s.eig"%seed):
                raise SystemExit("No INPUTS (.eig) for 1D Berry Curvature Calculation !")
            if not os.path.isfile("%s.chk"%seed):
                raise SystemExit("No INPUTS (.chk) for 1D Berry Curvature Calculation !")
            os.chdir(self.wrkdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/EFERMI %s/EFERMI" % (self.srcdir, self.objdir))
                os.system("cp %s/%s.win %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.amn %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.mmn %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.eig %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.chk %s/" % (self.srcdir, seed, self.objdir))
      
            
        if self.task.upper() == "BC2":

            os.chdir(self.srcdir)
            if not os.path.isfile("%s.win"%seed):
                raise SystemExit("No INPUTS (.win) for 2D Berry Curvature Calculation !")
            if not os.path.isfile("%s.amn"%seed):
                raise SystemExit("No INPUTS (.amn) for 2D Berry Curvature Calculation !")
            if not os.path.isfile("%s.mmn"%seed):
                raise SystemExit("No INPUTS (.mmn) for 2D Berry Curvature Calculation !")
            if not os.path.isfile("%s.eig"%seed):
                raise SystemExit("No INPUTS (.eig) for 2D Berry Curvature Calculation !")
            if not os.path.isfile("%s.chk"%seed):
                raise SystemExit("No INPUTS (.chk) for 2D Berry Curvature Calculation !")
            os.chdir(self.wrkdir)
            
            if self.objdir != self.srcdir:
                os.system("cp %s/EFERMI %s/EFERMI" % (self.srcdir, self.objdir))
                os.system("cp %s/%s.win %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.amn %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.mmn %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.eig %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.chk %s/" % (self.srcdir, seed, self.objdir))


        if self.task.upper() == "AHC":

            os.chdir(self.srcdir)
            if not os.path.isfile("%s.win"%seed):
                raise SystemExit("No INPUTS (.win) for Anomalous Hall Calculation !")
            if not os.path.isfile("%s.amn"%seed):
                raise SystemExit("No INPUTS (.amn) for Anomalous Hall Calculation !")
            if not os.path.isfile("%s.mmn"%seed):
                raise SystemExit("No INPUTS (.mmn) for Anomalous Hall Calculation !")
            if not os.path.isfile("%s.eig"%seed):
                raise SystemExit("No INPUTS (.eig) for Anomalous Hall Calculation !")
            if not os.path.isfile("%s.chk"%seed):
                raise SystemExit("No INPUTS (.chk) for Anomalous Hall Calculation !")
            os.chdir(self.wrkdir)
            
            if self.objdir != self.srcdir:
                os.system("cp %s/EFERMI %s/EFERMI" % (self.srcdir, self.objdir))
                os.system("cp %s/%s.win %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.amn %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.mmn %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.eig %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.chk %s/" % (self.srcdir, seed, self.objdir))

        if self.task.upper() == "OPT":

            os.chdir(self.srcdir)
            if not os.path.isfile("%s.win"%seed):
                raise SystemExit("No INPUTS (.win) for Optical Hall Calculation !")
            if not os.path.isfile("%s.amn"%seed):
                raise SystemExit("No INPUTS (.amn) for Optical Hall Calculation !")
            if not os.path.isfile("%s.mmn"%seed):
                raise SystemExit("No INPUTS (.mmn) for Optical Hall Calculation !")
            if not os.path.isfile("%s.eig"%seed):
                raise SystemExit("No INPUTS (.eig) for Optical Hall Calculation !")
            if not os.path.isfile("%s.chk"%seed):
                raise SystemExit("No INPUTS (.chk) for Optical Hall Calculation !")
            os.chdir(self.wrkdir)
            
            if self.objdir != self.srcdir:
                os.system("cp %s/%s.win %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.amn %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.mmn %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.eig %s/" % (self.srcdir, seed, self.objdir))
                os.system("cp %s/%s.chk %s/" % (self.srcdir, seed, self.objdir))
                
        os.chdir(self.objdir)

        if self.task.upper() == "WAN":
            
            with open("INCAR", "r") as f:
                incar = f.readlines()
                
            incar = [line for line in incar if "ALGO" not in line]
            incar = [line for line in incar if "NELM" not in line]
            incar = [line for line in incar if "NELMIN" not in line]                    
            incar = [line for line in incar if "LCHARG" not in line]
            incar = [line for line in incar if "LWAVE" not in line]
            incar = [line for line in incar if "ICHARG" not in line]
            incar = [line for line in incar if "ISTART" not in line]
            incar = [line for line in incar if "LORBIT" not in line]
            incar = [line for line in incar if "PRECFOCK" not in line]
            incar = [line for line in incar if "LHFCALC" not in line]
            incar = [line for line in incar if "HFSCREEN" not in line]
            incar = [line for line in incar if "AEXX" not in line]
            incar = [line for line in incar if "ENCUTGW" not in line]
            incar = [line for line in incar if "OMEGAMAX" not in line]
            incar = [line for line in incar if "NOMEGA" not in line]
            incar = [line for line in incar if "LOPTICS" not in line]
            incar = [line for line in incar if "CSHIFT" not in line]
            incar = [line for line in incar if "NEDOS" not in line]
            incar = [line for line in incar if "KPAR" not in line]
            incar = [line for line in incar if "NPAR" not in line]
            
            if npar is not None: 
                incar.append("NPAR = {} \n".format(npar))
            if kpar is not None: 
                incar.append("KPAR = {} \n".format(kpar))

            if istart is not None:
                incar.append("ISTART = 1 \n")
                incar.append("ALGO = None \n")
                incar.append("NELM = 1 \n")
            if icharg is not None: 
                incar.append("ICHARG = 11 \n")
                incar.append("ALGO = Normal \n")
                incar.append("NELM = 200 \n")    
                            
            incar.append("LWAVE = .FALSE. \n")
            incar.append("LCHARG = .FALSE. \n")
            incar.append("LWANNIER90 = .TRUE.")

            with open("INCAR", "w") as f:
                f.writelines(incar)
                
            index, atoms, cords = GetPoscar()

            file = open("%s.win" % seed, "w")
            
            # if "dis_win_min" in win:
            #     win["dis_win_min"] += fermi_energy

            # if "dis_win_max" in win:
            #     win["dis_win_max"] += fermi_energy

            # if "dis_froz_min" in win:
            #     win["dis_froz_min"] += fermi_energy

            # if "dis_froz_max" in win:
            #     win["dis_froz_max"] += fermi_energy

            file.write("begin projections \n")
            if proj.lower() == "random":
                file.write("%s \n" % proj)
            else:
                for key, val in proj.items():
                    for i in index:
                        xi, yi, zi = cords[i]
                        if key.isdigit() and i + 1 == int(key):
                            file.write("c = %.7f,%.7f,%.7f: %s \n" % (xi, yi, zi, val))
                    if key.isalpha():
                        file.write("%s: %s \n" % (key, val))
            file.write("end projections \n")

            for key, val in win.items():
                file.write("{} = {} \n".format(key, val))

            file.close()


        if self.task.upper() == "CHK":

            with open("%s.win" % seed, "r") as f:
                wannier = f.readlines()

            lookup_list = ["begin projections", "end projections", \
                           "guiding_centres", "spinors", "kmesh_tol", \
                           "num_wann", "num_bands", "exclude_bands", \
                           "dis_win_min", "dis_win_max", \
                           "dis_froz_min", "dis_froz_max", \
                           "dis_num_iter", "num_iter"]

            index = []
            for i, line in enumerate(wannier):
                 for lookup in lookup_list:
                     regex = r"\s*(%s)\s*" % lookup
                     if re.match(regex, line):
                        index.append(i)
                        
            k = max(index)
            wannier_select = wannier[0:k+1]
            wannier_remain = wannier[k+1::]

            fermi_energy = float(open("EFERMI", "r").read())
            
            if "dis_win_min" in win:
                win["dis_win_min"] += fermi_energy

            if "dis_win_max" in win:
                win["dis_win_max"] += fermi_energy

            if "dis_froz_min" in win:
                win["dis_froz_min"] += fermi_energy

            if "dis_froz_max" in win:
                win["dis_froz_max"] += fermi_energy

            for k, v in win.items():
                has_no_text = GrepText(k, wannier_select)
                if has_no_text:
                    wannier_select.append("{} = {} \n".format(k, v))

            for i, line in enumerate(wannier_select):
                for k, v in win.items():
                    if k in line: wannier_select[i] = "{} = {} \n".format(k, v)

            wannier = wannier_select + wannier_remain
            with open("%s.win" % seed, "w") as f:
                f.writelines(wannier)


        if self.task.upper() == "FIT":

            with open("%s.win" % seed, "r") as f:
                wannier = f.readlines()

            wannier.append("\n")

            for k, v in win.items():
                has_no_text = GrepText(k, wannier)
                if has_no_text:
                    wannier.append("{} = {} \n".format(k, v))

            with open("%s.win" % seed, "w") as f:
                f.writelines(wannier)
                
            os.system("cat kpath.win >> %s.win \n" % self.seed)
            
            
        if self.task.upper() == "FSF":

            with open("%s.win" % seed, "r") as f:
                wannier = f.readlines()

            fermi_energy = float(open("EFERMI", "r").read())

            if fermi_shift is None: fermi_shift = 0.0

            win.update(fermi_energy=fermi_energy+fermi_shift)

            wannier.append("\n")

            for k, v in win.items():
                has_no_text = GrepText(k, wannier)
                if has_no_text:
                    wannier.append("{} = {} \n".format(k, v))

            with open("%s.win" % seed, "w") as f:
                f.writelines(wannier)


        if self.task.upper() == "BC1":

            with open("%s.win" % seed, "r") as f:
                wannier = f.readlines()

            fermi_energy = float(open("EFERMI", "r").read())

            if fermi_shift is None: fermi_shift = 0.0
            
            win.update(fermi_energy=fermi_energy+fermi_shift)

            wannier.append("\n")

            for k, v in win.items():
                has_no_text = GrepText(k, wannier)
                if has_no_text:
                    wannier.append("{} = {} \n".format(k, v))

            with open("%s.win" % seed, "w") as f:
                f.writelines(wannier)

            os.system("cat kpath.win >> %s.win \n" % self.seed)
            
            
        if self.task.upper() == "BC2":

            with open("%s.win" % seed, "r") as f:
                wannier = f.readlines()
            
            fermi_energy = float(open("EFERMI", "r").read())
            
            if fermi_shift is None: fermi_shift = 0.0
            if berry_kmesh is None: berry_kmesh = (25, 25, 25)
   
            berry_kmesh = " ".join([str(k) for k in berry_kmesh]) 
            
            win.update(berry_kmesh=berry_kmesh)
            win.update(fermi_energy=fermi_energy+fermi_shift)

            wannier.append("\n")

            for k, v in win.items():
                has_no_text = GrepText(k, wannier)
                if has_no_text:
                    wannier.append("{} = {} \n".format(k, v))

            with open("%s.win" % seed, "w") as f:
                f.writelines(wannier)


        if self.task.upper() == "AHC":

            with open("%s.win" % seed, "r") as f:
                wannier = f.readlines()

            fermi_energy = float(open("EFERMI", "r").read())

            if berry_kmesh is None: berry_kmesh = (25, 25, 25)
            if fermi_energy_min is None: fermi_energy_min = -2.5
            if fermi_energy_max is None: fermi_energy_max = 2.5
            if fermi_energy_step is not None: fermi_energy_step = 0.001
            
            berry_kmesh = " ".join([str(k) for k in berry_kmesh]) 
            
            # win.update("fermi_energy_min = {}".format(fermi_energy_min + fermi_energy))
            # win.update("fermi_energy_max = {}".format(fermi_energy_max + fermi_energy))
            # win.update("fermi_energy_step = {}".format(fermi_energy_step))

            win.update(berry_kmesh=berry_kmesh)
            win.update(fermi_energy_min=fermi_energy_min+fermi_energy)
            win.update(fermi_energy_max=fermi_energy_max+fermi_energy)
            win.update(fermi_energy_step=fermi_energy_step)

            wannier.append("\n")

            for k, v in win.items():
                has_no_text = GrepText(k, wannier)
                if has_no_text:
                    wannier.append("{} = {} \n".format(k, v))

            with open("%s.win" % seed, "w") as f:
                f.writelines(wannier)
                

        if self.task.upper() == "OPT":

            with open("%s.win" % seed, "r") as f:
                wannier = f.readlines()

            if berry_kmesh is None: berry_kmesh = (25, 25, 25)
            if fermi_energy is None: fermi_energy = 0.0 
            if fermi_shift is None: fermi_shift = 0.0
            if kubo_freq_min is None: kubo_freq_min = 0.0
            if kubo_freq_max is None: kubo_freq_max = 10.0
            if kubo_freq_step is None: kubo_freq_step = 0.001
            if kubo_smr_type is None: kubo_smr_type = "gauss"
            if kubo_smr_fixed_en_width is None: kubo_smr_fixed_en_width = 0.01
    
            berry_kmesh = " ".join([str(k) for k in berry_kmesh]) 
            
            # win.update("kubo_freq_min = {}".format(kubo_freq_min))
            # win.update("kubo_freq_max = {}".format(kubo_freq_max))
            # win.update("kubo_freq_step = {}".format(kubo_freq_step))
            # win.update("kubo_smr_type = {}".format(kubo_smr_type))
            # win.update("kubo_smr_fixed_en_width = {}".format(kubo_smr_fixed_en_width))

            win.update(berry_kmesh=berry_kmesh)
            win.update(fermi_energy=fermi_energy+fermi_shift)
            win.update(kubo_freq_min=kubo_freq_min)
            win.update(kubo_freq_max=kubo_freq_max)
            win.update(kubo_freq_step=kubo_freq_step)
            win.update(kubo_smr_type=kubo_smr_type)
            win.update(kubo_smr_fixed_en_width=kubo_smr_fixed_en_width)
            
            wannier.append("\n")

            for k, v in win.items():
                has_no_text = GrepText(k, wannier)
                if has_no_text:
                    wannier.append("{} = {} \n".format(k, v))

            with open("%s.win" % seed, "w") as f:
                f.writelines(wannier)
                
        os.chdir(self.wrkdir)


    def Submission(self, 
                   queue="LOCAL", 
                   jobname=None, 
                   ntasks=None, 
                   lprocs=None, 
                   nnodes=None, 
                   ncores=None,
                   walltime=None, 
                   prefix=None):

        self.queue = queue

        os.chdir(self.objdir)
        
        if queue is None: queue = "LOCAL"
        if jobname is None: jobname = self.task
        # if ntasks is None: ntasks = 1 # Mumber of MPI process
        # if lprocs is None: lprocs = 1 # Number of OMP threads
        # if nnodes is None: nnodes = 1 # Number of Nodes
        # if ncores is None: ncores = 1 # Number of cores per Node
        if walltime is None: walltime = "01:00:00"
        if prefix is None: prefix = "mpirun"

        if nnodes is not None and ncores is not None: 
            ntasks = int(nnodes * ncores)
            
        if queue.upper() == "LOCAL":

            runjob = open("runjob", "w")
            
            runjob.write("#!/bin/bash \n")
            runjob.write("\n")
            
            self.RunJobWrite(runjob, ntasks=ntasks, prefix=prefix)

            runjob.close()

        if queue.upper() == "SLURM":
                        
            runjob = open("runjob", "w")
            
            runjob.write("#!/bin/bash \n")
            runjob.write("\n")
            
            runjob.write('#SBATCH --job-name={} \n'.format(jobname))
            if ntasks is not None and lprocs is not None:
                runjob.write('#SBATCH --ntasks={} --cpus-per-task={} \n'.format(ntasks, lprocs))
            if nnodes is not None and ncores is not None:
                runjob.write('#SBATCH --nodes={} --ntasks-per-node={} \n'.format(nnodes, ncores))
            runjob.write('#SBATCH --time={} \n'.format(walltime))
            runjob.write('#SBATCH --mail-user=longlong.li@uantwerpen.be \n')
            runjob.write('#SBATCH -o {}-%j.out \n'.format(jobname))
            runjob.write('#SBATCH -e {}-%j.err \n'.format(jobname))

            runjob.write("\n")

            runjob.write("module load --purge \n")
            runjob.write("module load calcua/2020a \n")
            runjob.write("module load intel/2020a \n")
            runjob.write("module load VASP/5.4.4-intel-2020a-Wannier90-2.1.0 \n")
            
            runjob.write("\n")

            self.RunJobWrite(runjob, ntasks=ntasks, prefix=prefix)
                            
            runjob.close()

        if queue.upper() == "PBS":
                        
            runjob = open("runjob", "w")
            
            runjob.write("#!/bin/bash \n")
            runjob.write("\n")

            runjob.write('#PBS -N {} \n'.format(jobname))
            if ntasks is not None and lprocs is not None:
                runjob.write('#PBS -L tasks={}:lprocs={} \n'.format(ntasks, lprocs))
            if nnodes is not None and ncores is not None:
                runjob.write('#PBS -l nodes={}:ppn={} \n'.format(nnodes, ncores))              
            runjob.write('#PBS -l walltime={} \n'.format(walltime))
            runjob.write('#PBS -M longlong.li@uantwerpen.be \n')

            runjob.write("\n")
            runjob.write("cd $PBS_O_WORKDIR \n")
            runjob.write("\n")
                                    
            runjob.write("module load --purge \n")
            runjob.write("module load calcua/2020a \n")
            runjob.write("module load intel/2020a \n")
            runjob.write("export PATH=$VSC_DATA/VASP2WANN/bin:$PATH \n")
            
            runjob.write("\n")
            
            self.RunJobWrite(runjob, ntasks=ntasks, prefix=prefix)
            
            runjob.close()
            
        if self.queue.upper() == "LOCAL": 
            os.system("bash runjob")
            
        if self.queue.upper() == "SLURM": 
            os.system("sbatch runjob")
            
        if self.queue.upper() == "PBS": 
            os.system("qsub runjob")  

        os.chdir(self.wrkdir)
        

    def RunJobWrite(self, runjob, ntasks=1, prefix="mpirun"):
        
        if self.task.upper() == "WAN":
            if self.lsoc:
                if prefix == "mpirun":
                    runjob.write("mpirun -np {} vasp_ncl \n".format(ntasks))
                if prefix == "srun":
                    runjob.write("srun -n {} vasp_ncl \n".format(ntasks))
            else:
                if prefix == "mpirun":
                    runjob.write("mpirun -np {} vasp_std \n".format(ntasks))
                if prefix == "srun":
                    runjob.write("srun -n {} vasp_std \n".format(ntasks))
                    
        if self.task.upper() == "CHK":
            runjob.write("wannier90.x {} \n".format(self.seed))
            
        if self.task.upper() == "FIT":
            runjob.write("wannier90.x {} \n".format(self.seed))

        if self.task.upper() == "FSF":
            runjob.write("wannier90.x {} \n".format(self.seed))

        if self.task.upper() == "BC1":
            if prefix == "mpirun":
                runjob.write("mpirun -np {} postw90.x {} \n".format(ntasks, self.seed))
            if prefix == "srun":
                runjob.write("srun -n {} postw90.x {} \n".format(ntasks, self.seed))

        if self.task.upper() == "BC2":
            if prefix == "mpirun":
                runjob.write("mpirun -np {} postw90.x {} \n".format(ntasks, self.seed))
            if prefix == "srun":
                runjob.write("srun -n {} postw90.x {} \n".format(ntasks, self.seed))
                
        if self.task.upper() == "AHC":
            if prefix == "mpirun":
                runjob.write("mpirun -np {} postw90.x {} \n".format(ntasks, self.seed))
            if prefix == "srun":
                runjob.write("srun -n {} postw90.x {} \n".format(ntasks, self.seed))
                
        if self.task.upper() == "OPT":
            if prefix == "mpirun":
                runjob.write("mpirun -np {} postw90.x {} \n".format(ntasks, self.seed))
            if prefix == "srun":
                runjob.write("srun -n {} postw90.x {} \n".format(ntasks, self.seed))
                
                
class WanFitMP(PythWann):

    def __init__(self, wrkdir=None, srcdir=None, objdir=None):

        self.srcdir = srcdir
        self.objdir = objdir
        self.wrkdir = wrkdir

    def FitWrap(self, kwargs):

        self.RunFit(**kwargs)

    def RunFit(self, dis_win_min=None, dis_win_max=None, dis_froz_min=None, dis_froz_max=None, num_iter=None, dis_num_iter=None, fermi_energy=None, bands_kpath=None):

        WRKDIR = self.wrkdir
        WANDIR = self.srcdir

        if dis_win_min is not None and dis_win_max is None:
            WAMDIR = "%s_WMIN_%.1f" % (WANDIR, dis_win_min)
        if dis_win_min is None and dis_win_max is not None:
            WAMDIR = "%s_WMAX_%.1f" % (WANDIR, dis_win_max)
        if dis_win_min is not None and dis_win_max is not None:
            WAMDIR = "%s_WMIN_%.1f_WMAX_%.1f" % (WAN, dis_win_min, dis_win_max)
        if dis_froz_min is not None and dis_froz_max is None:
            WAMDIR = "%s_FMIN_%.1f" % (WANDIR, dis_froz_min)
        if dis_froz_min is None and dis_froz_max is not None:
            WAMDIR = "%s_FMAX_%.1f" % (WANDIR, dis_froz_max)
        if dis_froz_min is not None and dis_froz_max is not None:
            WAMDIR = "%s_FMIN_%.1f_FMAX_%.1f" % (WANDIR, dis_froz_min, dis_froz_max)

        # WAM = self.objdir + "/" + WAM # os.path.join(self.objdir, WAM)

        Job = PythWann(task="FIT", wrkdir=WRKDIR, srcdir=WANDIR, objdir=WAMDIR)
        Job.MakeInput(dis_win_min=dis_win_min, dis_win_max=dis_win_max, dis_froz_min=dis_froz_min, dis_froz_max=dis_froz_max, num_iter=num_iter, dis_num_iter=dis_num_iter, fermi_energy=fermi_energy, bands_kpath=bands_kpath)
        Job.Submission(queue="LOCAL")

    def RunPool(self, params=None):

        if len(params) > mp.cpu_count():
            pool = mp.Pool(processes=mp.cpu_count())
        else:
            pool = mp.Pool(processes=len(params))
        pool.map(self.FitWrap, params)
        pool.close()
        pool.join()


# if __name__ == "__main__":
#     main()

    # import PythWann

    # PWD = os.getcwd()
    
    # Job = PythWann(task="WAN", topdir=PWD, wrkdir="MoS2", srcdir="SCF", objdir="WAN")
    # Job.MakeInput(icharg=1, lsoc=True, guiding_centres=True, projections={"Mo": "d"})
    # Job.Submission(queue="local", ntasks=40) # queue = local/slurm/pbs

    # Job = PythWann(task="FIT", topdir=PWD, wrkdir="MoS2", srcdir="WAN", objdir="WAN")
    # Job.MakeInput(bands_plot=True, bands_num_points=100, write_hr=True)
    # Job.Submission(queue="local", ntasks=1) # queue = local/slurm/pbs

    # Job = PythWann(task="AHC", topdir=PWD, wrkdir="MoS2", srcdir="WAN", objdir="WAN")
    # job.MakeInput(berry=True, berry_task="ahc", berry_kmesh="300 300 10")
    # Job.Submission(queue="local", ntasks=40) # queue = local/slurm/pbs
