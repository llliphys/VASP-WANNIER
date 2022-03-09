import os
import re
import yaml
import numpy as np

# check pymatgen installed or not
try:
    import pymatgen
    pymatgen_install = True
except ImportError:
    pymatgen_install = False

try:
    import ase
    ase_install = True
except ImportError:
    ase_install = False
    
def DeleteLarge():

    inputs = ["WAVECAR", "CHGCAR", "WAVEDER"]

    for f in os.listdir("."):
        if f in inputs: os.system("rm -rf %s" % f)
        
def MakeClean():

    inputs = ["INCAR", "POSCAR", "POTCAR", "KPOINTS"]

    for f in os.listdir("."):
        if f not in inputs: os.system("rm -rf %s" % f)

    for f in os.listdir("."):
        if "POSCAR-" in f: os.system("rm -rf %s" % f)


def GetLattice(lattice=None):

    avec = lattice[2].strip().strip("\n")
    bvec = lattice[3].strip().strip("\n")
    cvec = lattice[4].strip().strip("\n")

    avec = np.array([float(s) for s in avec.split()])
    bvec = np.array([float(s) for s in bvec.split()])
    cvec = np.array([float(s) for s in cvec.split()])

    a = np.sqrt(np.sum(avec**2))
    b = np.sqrt(np.sum(bvec**2))
    c = np.sqrt(np.sum(cvec**2))

    cosab = avec.dot(bvec) / (np.linalg.norm(avec) * np.linalg.norm(bvec))
    cosbc = bvec.dot(cvec) / (np.linalg.norm(bvec) * np.linalg.norm(cvec))
    cosac = avec.dot(cvec) / (np.linalg.norm(avec) * np.linalg.norm(cvec))

    return avec, bvec, cvec, a, b, c


def ModifyCell(poscar=None, contcar=None):

    # motivation: to remedy poscar shape after relaxation
    # 1. read the initial poscar to get the initial lattice vectors
    # 2. calculate the angles between the initial lattice vectors
    # 3. modify the contcar lattce vectors according to the initial angles

    if poscar is None: poscar = "POSCAR"
    if contcar is None: contcar = "CONTCAR"

    with open(poscar, "r") as f: poscar = f.readlines()
    with open(contcar, "r") as f: contcar = f.readlines()

    avec_1, bvec_1, cvec_1, a_1, b_1, c_1 = GetLattice(lattice=poscar)
    avec_2, bvec_2, cvec_2, a_2, b_2, c_2 = GetLattice(lattice=contcar)

    avec = avec_1 * (a_2 / a_1)
    bvec = bvec_1 * (b_2 / b_1)
    cvec = cvec_1 * (c_2 / c_1)

    ax, ay, az = avec
    bx, by, bz = bvec
    cx, cy, cz = cvec

    contcar[2] = "  %.15f  %.15f  %.15f \n" % (ax, ay, az)
    contcar[3] = "  %.15f  %.15f  %.15f \n" % (bx, by, bz)
    contcar[4] = "  %.15f  %.15f  %.15f \n" % (cx, cy, cz)

    natom_str = poscar[6].strip().split()
    natom_arr = np.array([int(i) for i in natom_str])
    natoms = np.sum(natom_arr)

    for i in np.arange(8, 8+natoms):
        line = poscar[i].strip().strip("\n")
        x, y, z = [float(s) for s in line.split()]
        contcar[i] = "  %.15f  %.15f  %.15f \n" %(x, y, z)

    with open(infile, "w") as f:
        f.writelines(contcar)

def MakeKmesh(kmin=-0.1, kmax=0.1, nkpts=21, kz=0):  # for 2D only

    """
    This function generates a 2D Kmesh for e.g. a Fermi surface calculation.
    """

    kxs = np.linspace(kmin, kmax, nkpts)
    kys = np.linspace(kmin, kmax, nkpts)

    kmesh = open("KPOINTS", "w")
    kmesh.write("2D KMESH\n")
    kmesh.write("%d\n" % (nkpts * nkpts))
    kmesh.write("Cartesian\n")
    for kx in kxs:
        for ky in kys:
            kmesh.write("%.12f  %.12f  %.12f  %.1f\n" % (kx, ky, kz, 1.0))
    kmesh.close()


def MakeKpath(kpath=None, nkpts=None):

    """
    This function uses PYMATGEN to generate a kpath for a band structure calculation.
    """

    if pymatgen_install:
        from pymatgen.io.vasp.inputs import Kpoints
        from pymatgen.core import Structure
        from pymatgen.symmetry.bandstructure import HighSymmKpath
    else:
        print("Warning: PYMATGEN is not yet installed on your machine!")
        exit(1)

    if nkpts is None: nkpts = 100

    if kpath is not None:
        kpath = [[str(s) for s in kpath]]
        ibz = HighSymmKpath(Structure.from_file("POSCAR"))
        ibz.kpath['path'] = kpath
        Kpoints.automatic_linemode(nkpts, ibz).write_file("KPOINTS")
    else:
        ibz = HighSymmKpath(Structure.from_file("POSCAR"))
        Kpoints.automatic_linemode(nkpts, ibz).write_file("KPOINTS")

    with open("KPOINTS", "r") as f:
        klines = f.readlines()

    regex = r"\s*(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*\!\s*\\*([A-Z]*)[a-z]*"

    file = open("kpath.win", "w")
    file.write("\n")
    file.write("bands_plot = true \n")
    file.write("bands_num_points = %g \n" % nkpts)
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

def CheckConv(infile="OUTCAR"):

    """
    This function checks whether a DFT calculation is converged.
    """

    if not os.path.exists(infile):
        print("Error: No %s Exist !" % infile)
        exit(1)

    with open(infile, "r") as f:
        outcar = f.readlines()

    converge = False

    rlx_conv = False
    for line in outcar[::-1]:
        if "reached required accuracy" in line:
            rlx_conv = True
            break

    scf_conv = False
    for line in outcar[::-1]:
        if "aborting loop because EDIFF is reached" in line:
            scf_conv = True
            break

    if rlx_conv: 
        converge = rlx_conv
    
    if scf_conv: 
        converge = scf_conv
    
    return converge


def GetEnergy(infile="OSZICAR"):

    if not os.path.exists(infile):
        print("Error: No %s Exist !" % infile)
        exit(1)

    with open(infile, "r") as f:
        flines = f.readlines()

    energy = float(flines[-1].strip().split()[2])

    return energy

def GetNkpts(infile="OUTCAR"):

    if not os.path.exists(infile):
        print("Error: No %s Exist !" % infile)
        exit(1)

    with open(from_file, "r") as f:
        fread = f.read()

    regex = r"\s*k-points\s+(NKPTS\s*=\s+\d+)"
    match = re.findall(regex, fread)[-1]
    nkpts = int(match.split("=")[1])

    return nkpts

def GetCPUTime(infile="OUTCAR"):

    if not os.path.exists(infile):
        print("Error: No %s Exist !" % infile)
        exit(1)

    with open(infile, "r") as f:
        flines = f.readlines()

    for line in flines[::-1]:
        if "Elapsed time (sec)" in line:
            cputime = float(line.strip().split(":")[1])
            break

    return cputime


def GetMagmom(infile="OUTCAR"):

    if not os.path.exists(infile):
        print("Error: No %s Exist !" % infile)
        exit(1)

    with open(infile, "r") as f:
        outcar = f.readlines()

    natoms = None
    for i, line in enumerate(outcar[::-1]):
        if "number of ions" in line:
            natoms = int(line.strip().split("NIONS =")[1])
            print("Number of atoms in the cell: ", natoms)
            break

    mspn = None
    for i, line in enumerate(outcar):
        if "magnetization (x)" in line:
            mspn = outcar[i+4:i+4+natoms]

    magmom = []
    if mspn is not None:
        for i, line in enumerate(mspn):
            line = line.strip().split()
            ms = round(float(line[-1]), 6)
            magmom.append(ms)

    return magmom


def GetEfermi(infile="OUTCAR"):

    if not os.path.exists(infile):
        print("Error: No %s Exist !" % infile)
        exit(1)

    outcar = open(infile, "r").read()
    match = re.findall(r"E-fermi\s*:\s*(-?\d+.\d+)", outcar)[-1]
    efermi = float(match)

    return efermi

def GetNbands(infile="OUTCAR"):

    if not os.path.exists(infile):
        print("Error: No %s Exist !" % infile)
        exit(1)

    outcar = open(infile, "r").read()
    match = re.findall(r"NBANDS\s*=\s*(\d+)", outcar)[-1]
    nbands = int(match)

    return nbands


def GetPhBand(infile="band.yaml", outfile="PhononBand.dat", units="THz"):
    
    """
    Here units can be chosen among THz, meV, and 1/cm
    1 THz = 4.13 meV = 72 1/cm (to be checked !!!!!!!)
    """
        
    with open(infile, "r") as bandfile:
        data = yaml.safe_load(bandfile)
        
    phonon = data["phonon"]
    natoms = data["natom"]
    nkpts = len(phonon)

    lines = []
    for i in range(nkpts):
        distance = phonon[i]["distance"]
        band = phonon[i]["band"]
        nbands = len(band)
        lines.append("%.12f\t" % distance)
        for j in range(nbands):
            eigval = band[j]["frequency"]
            lines.append("%.12f\t" % eigval)
        lines.append("\n" )

    with open(outfile, "w") as f:
        f.writelines(lines)
            
                    
class PythVasp:

    """
    - A Python class for carrying out a series of DFT-VASP calculations including e.g.
    "Relaxation", "Static SCF", "Band Structure", "Density of States", "Charge Density"
    "Fermi Surface", "Hybrid Functional", "GW Calculation", "BSE calculation", which are
    reprensted here by task="RLX", "CON", "SCF", "BND", "DOS", "CHG", "FSF", "HSE", "GWA",
    and "BSE", respectively. Here task="CON" refers to as continuing a calculation on top
    of a preconverged one.

    - Before doing any calculations, one has always to specify two working directories,
    one of which contains the input files for VASP calculations and the other contains
    the output files from the calculations. The two directories for inputs and outputs are
    labeled by scrdir and objdir, respectively. Note that scrdir and objdir can be the same
    or different directories, depending on specific calculations. Optionally, a rootdir can
    also be specified to indicate the root directory which contains the two subdirectories
    i.e., rootdir/srcdir and rootdir/objdir.

    - A child class based on this parent class is implemented for testing the convergence
    against various important DFT parameters (e.g. ENCUT, SIGMA, KPOINTS). This testing can
    be carried out by setting task="RLX" or "SCF" with the latter for a single-point testing.

    - Example scripts are given in the end to show the basic lines of how to use these modules.

    """

    def __init__(self, task="RLX", topdir=None, wrkdir=None, srcdir=None, objdir=None):

        os.chdir(topdir)
        
        if wrkdir is None:
            wrkdir = os.getcwd()
        if srcdir is None:
            print("No INPUT Directory for VASP Calculation !")
            exit(1)
        if objdir is None:
            print("No OUTPUT Directory for VASP Calculation !")
            exit(1)    
                    
        self.task = task
        self.wrkdir = os.path.abspath(wrkdir)
        self.srcdir = os.path.join(self.wrkdir, srcdir)
        self.objdir = os.path.join(self.wrkdir, objdir)
        
        if not os.path.exists(self.srcdir): 
            os.mkdir(self.srcdir)
        
        if not os.path.exists(self.objdir): 
            os.mkdir(self.objdir)

    def MakeInput(self, encut=None, prec=None, symprec=None, \
                  nsw=None, isif=None, ibrion=None, potim=None, ediffg=None, \
                  algo=None, ialgo=None, nelm=None, nelmin=None, ediff=None, \
                  icharg=None, istart=None, lorbit=None, \
                  ismear=None, sigma=None, nedos=None, \
                  nelect=None, \
                  isym=None, ivdw=None, \
                  ispin=None, magmom=None, \
                  lsoc=False, saxis=None, lorbmom=False, \
                  lwave=False, lcharg=False, \
                  lhfcalc=False, lmodelhf=False, hfscreen=None, aexx=None, precfock=None, \
                  lvtot=False, lvhar=False, \
                  lreal=False, addgrid=False, \
                  ldau=False, utype=None, ldaul=None, ldauu=None, ldauj=None, \
                  amix=None, bmix=None, amix_mag=None, bmix_mag=None, \
                  loptics=False, cshift=None, omegamax=None, \
                  nbands=None, encutgw=None, nomega=None, \
                  qpoint=None, nbandso=None, nbandsv=None, nbseeig=None, scissor=None, \
                  ldipol=False, efield=None, idipol=None, dipol=None, \
                  lpard=False, nbmod=None, eint=None, iband=None, kpuse=None, \
                  lsepk=False, lsepb=False, \
                  nkpts=None, gamma=True, shift=None, kpath=None, \
                  make_kmesh=None, kmin=None, kmax=None, \
                  kptfile=None, \
                  poscar=None, potcar=None, optcell=None, \
                  npar=None, kpar=None, ncore=None, \
                  lwannier=False, wanfile=None, \
                  celldim=None, \
                  make_clean=False, \
                  delete_large=False):

        self.lsoc = lsoc
            
        incar = dict() # OrderedDict()

        if encut is not None: incar.update(ENCUT=encut)
        if prec is not None: incar.update(PREC=prec)
        if symprec is not None: incar.update(SYMPREC=symprec)
        if precfock is not None: incar.update(PRECFOCK=precfock)
        if nsw is not None: incar.update(NSW=nsw)
        if isif is not None: incar.update(ISIF=isif)
        if ibrion is not None: incar.update(IBRION=ibrion)
        if ediffg is not None: incar.update(EDIFFG=ediffg)
        if potim is not None: incar.update(POTIM=potim)
        if algo is not None: incar.update(ALGO=algo)
        if ialgo is not None: incar.update(IALGO=ialgo)
        if nelm is not None: incar.update(NELM=nelm)
        if nelmin is not None: incar.update(NELMIN=nelmin)
        if ediff is not None: incar.update(EDIFF=ediff)
        if icharg is not None: incar.update(ICHARG=icharg)
        if istart is not None: incar.update(ISTART=istart)
        if lorbit is not None: incar.update(LORBIT=lorbit)
        if ismear is not None: incar.update(ISMEAR=ismear)
        if sigma is not None: incar.update(SIGMA=sigma)
        if nedos is not None: incar.update(NEDOS=nedos)
        if nelect is not None: incar.update(NELECT=nelect)
        if ivdw is not None: incar.update(IVDW=ivdw)
        if isym is not None: incar.update(ISYM=isym)
        if lwave: incar.update(LWAVE=".TRUE.")
        if lcharg: incar.update(LCHARG=".TRUE.")
        if not lwave: incar.update(LWAVE=".FALSE.")
        if not lcharg: incar.update(LCHARG=".FALSE.")
        if lhfcalc: incar.update(LHFCALC=".TRUE.")
        if lmodelhf: incar.update(LMODELHF=".TRUE.")
        if hfscreen is not None: incar.update(HFSCREEN=hfscreen)
        if aexx is not None: incar.update(AEXX=aexx)
        if lvtot: incar.update(LVTOT=".TRUE.")
        if lvhar: incar.update(LVHAR=".TRUE.")
        if lreal: incar.update(LREAL="Auto")
        if addgrid: incar.update(ADDGRID=".TRUE.")
        if loptics: incar.update(LOPTICS=".TRUE.")
        if cshift is not None: incar.update(CSHIFT=cshift)
        if omegamax is not None: incar.update(OMEGAMAX=omegamax)
        if ldau: incar.update(LDAU=".TRUE.")
        if utype is not None: incar.update(LDAUTYPE=utype)
        if ldaul is not None and "2" in ldaul:
            incar.update(LDAUPRINT=2, LMAXMIX=4)
        if ldaul is not None and "3" in ldaul:
            incar.update(LDAUPRINT=2, LMAXMIX=6)
        if ldaul is not None: incar.update(LDAUL=ldaul)
        if ldauu is not None: incar.update(LDAUU=ldauu)
        if ldauj is not None: incar.update(LDAUJ=ldauj)
        if amix is not None: incar.update(AMIX=amix)
        if bmix is not None: incar.update(BMIX=bmix)
        if amix_mag is not None: incar.update(AMIX_MAG=amix_mag)
        if bmix_mag is not None: incar.update(BMIX_MAG=bmix_mag)
        if lsoc: incar.update(LSORBIT=".TRUE.")
        if lorbmom: incar.update(LORBMOM=".TRUE.")
        if ispin is not None: incar.update(ISPIN=ispin)
        if magmom is not None: incar.update(MAGMOM=magmom)
        if saxis is not None: incar.update(SAXIS=saxis)
        if nbands is not None: incar.update(NBANDS=nbands)
        if encutgw is not None: incar.update(ENCUTGW=encutgw)
        if nomega is not None: incar.update(NOMEGA=nomega)
        if nbandso is not None: incar.update(NBANDSO=nbandso)
        if nbandsv is not None: incar.update(NBANDSV=nbandsv)
        if nbseeig is not None: incar.update(NBSEEIG=nbseeig)
        if qpoint is not None: incar.update(KPOINTS_BSE=qpoint)
        if scissor is not None: incar.update(SCISSOR=scissor)
        if efield is not None: incar.update(EFIELD=efield)
        if idipol is not None: incar.update(IDIPOL=idipol)
        if ldipol: incar.update(LDIPOL=".TRUE.")
        if dipol is not None: incar.update(DIPOL=dipol)
        if lpard: incar.update(LPARD=".TRUE.")
        if nbmod is not None: incar.update(NBMOD=nbmod)
        if eint is not None: incar.update(EINT=eint)
        if iband is not None: incar.update(IBAND=iband)
        if kpuse is not None: incar.update(KPUSE=kpuse)
        if lsepk: incar.update(LSEPK=".TRUE.")
        if lsepb: incar.update(LSEPB=".TRUE.")
        if kpar is not None: incar.update(KPAR=kpar)
        if npar is not None: incar.update(NPAR=npar)
        if ncore is not None: incar.update(NCORE=ncore)
        if lwannier: incar.update(LWANNIER90_RUN=".TRUE.")

        # win= dict()

        # if lsoc: win.update(spinors="true")
        # if lgct: win.update(guiding_centres="true")
        # if proj is None: win.update(projections="random")
        
        if kptfile is not None: 
            self.kptfile = os.path.join(self.wrkdir, kptfile)

        if wanfile is not None: 
            self.wanfile = os.path.join(self.wrkdir, wanfile)    
    
        # if nbands is None and self.task.upper() == "SCF":
        #     os.chdir(self.srcdir)
        #     nbands = GetNbands()
        #     os.chdir(self.wrkdir)
        #     incar.update(NBANDS=nbands)
        #     if lsoc: incar.update(NBANDS=int(2*nbands))

        # if nbands is None and self.task.upper() == "BND":
        #     os.chdir(self.srcdir)
        #     nbands = GetNbands()
        #     os.chdir(self.wrkdir)
        #     incar.update(NBANDS=int(2*nbands))

        # if nbands is None and self.task.upper() == "DOS":
        #     os.chdir(self.srcdir)
        #     nbands = GetNbands()
        #     os.chdir(self.wrkdir)
        #     incar.update(NBANDS=int(2*nbands))
            
        # if nbands is None and self.task.upper() == "GWA":
        #     os.chdir(self.srcdir)
        #     nbands = GetNbands()
        #     os.chdir(self.wrkdir)
        #     incar.update(NBANDS=int(4*nbands))
                        
        # if magmom is None and self.task.upper() == "SCF":
        #    os.chdir(self.srcdir)
        #    magmom = GetMagmom()
        #    os.chdir(self.wrkdir)
        #    if magmom != []:
        #        if not lsoc: magmom = "  ".join([str(m) for m in magmom])
        #        if lsoc: magmom = "  ".join(["0  0  %.3f" % m for m in magmom])
        #        incar.update(MAGMOM=magmom)
        
        os.chdir(self.objdir)

        if make_clean: MakeClean()
            
        with open("INCAR", "w") as f:
            for key, val in incar.items():
                f.write("{} = {} \n".format(key, val))
                
        if optcell is not None:
            with open("OPTCELL", "w") as f:
                for ijk in optcell:
                    f.write("{}{}{} \n".format(ijk[0], ijk[1], ijk[2]))
                    
        with open("KPOINTS", "w") as f:
            f.write("K-Mesh \n")
            f.write("0 \n")
            f.write("Gamma \n") if gamma else f.write("MP \n")
            if type(nkpts) == list or type(nkpts) == tuple:
                f.write("{} {} {} \n".format(nkpts[0], nkpts[1], nkpts[2]))
            if type(nkpts) == int:
                f.write("{} {} {} \n".format(nkpts, nkpts, nkpts))
            if nkpts is None:
                f.write("{} {} {} \n".format(10, 10, 10))
            if shift is None:
                f.write("0 0 0 \n")
            else:
                f.write("{} {} {} \n".format(shift[0], shift[1], shift[2]))
                
        if kptfile is not None:
            os.system("cp %s KPOINTS" % self.kptfile)

        if wanfile is not None:
            os.system("cp %s wannier90.win" % self.wanfile)
                  
        os.chdir(self.wrkdir)

        if self.task.upper() == "RLX":

            if poscar is None: 
                poscar = "POSCAR"
            if potcar is None: 
                potcar = "POTCAR"

            os.chdir(self.srcdir)
            if not os.path.isfile(poscar):
                raise SystemExit("No POSCAR for DFT Relaxation !")
            if not os.path.isfile(potcar):
                raise SystemExit("No POTCAR for DFT Relaxation !")
            os.chdir(self.wrkdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/%s %s/POSCAR" % (self.srcdir, poscar, self.objdir))
                os.system("cp %s/%s %s/POTCAR" % (self.srcdir, potcar, self.objdir))
                
        if self.task.upper() == "CON":

            os.chdir(self.srcdir)
            converge = CheckConv()
            if not converge: 
                raise SystemExit("Relaxation is NOT Converged !")
            os.chdir(self.wrkdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/CONTCAR %s/POSCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/POTCAR %s/POTCAR" % (self.srcdir, self.objdir))
                if "ICHARG" in incar:
                    os.system("cp %s/CHGCAR %s/CHGCAR" % (self.srcdir, self.objdir))
                if "ISTART" in incar:
                    os.system("cp %s/WAVECAR %s/WAVECAR" % (self.srcdir, self.objdir))
            else:
                os.system("cp %s/CONTCAR %s/POSCAR" % (self.srcdir, self.objdir))

        if self.task.upper() == "SCF":

            os.chdir(self.srcdir)
            converge = CheckConv()
            if not converge: 
                raise SystemExit("Relaxation is NOT Converged !")
            os.chdir(self.wrkdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/CONTCAR %s/POSCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/POTCAR %s/POTCAR" % (self.srcdir, self.objdir))
                if "ICHARG" in incar:
                    os.system("cp %s/CHGCAR %s/CHGCAR" % (self.srcdir, self.objdir))
                if "ISTART" in incar:
                   os.system("cp %s/WAVECAR %s/WAVECAR" % (self.srcdir, self.objdir))

            
        if self.task.upper() == "CHG":

            os.chdir(self.srcdir)
            converge = CheckConv()
            if not converge: 
                raise SystemExit("Static SCF is NOT Converged !")
            os.chdir(self.wrkdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/CONTCAR %s/POSCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/POTCAR %s/POTCAR" % (self.srcdir, self.objdir))
                if "ICHARG" in incar:
                    os.system("cp %s/CHGCAR %s/CHGCAR" % (self.srcdir, self.objdir))
                if "ISTART" in incar:
                   os.system("cp %s/WAVECAR %s/WAVECAR" % (self.srcdir, self.objdir))

            
        if self.task.upper() == "FSF":

            os.chdir(self.srcdir)
            converge = CheckConv()
            if not converge: 
                raise SystemExit("Static SCF is NOT Converged !")
            os.chdir(self.wrkdir)

            os.chdir(self.objdir)
            if make_kmesh: MakeKmesh(kmin=kmin, kmax=kmax, nkpts=nkpts)
            os.chdir(self.wrkdir)
            
            if self.objdir != self.srcdir:
                os.system("cp %s/CONTCAR %s/POSCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/POTCAR %s/POTCAR" % (self.srcdir, self.objdir))
                if "ICHARG" in incar:
                    os.system("cp %s/CHGCAR %s/CHGCAR" % (self.srcdir, self.objdir))
                if "ISTART" in incar:
                   os.system("cp %s/WAVECAR %s/WAVECAR" % (self.srcdir, self.objdir))


        if self.task.upper() == "DOS":

            os.chdir(self.srcdir)
            converge = CheckConv()
            if not converge: 
                raise SystemExit("Static SCF is NOT Converged !")
            os.chdir(self.wrkdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/POSCAR %s/POSCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/POTCAR %s/POTCAR" % (self.srcdir, self.objdir))
                if "ICHARG" in incar:
                    os.system("cp %s/CHGCAR %s/CHGCAR" % (self.srcdir, self.objdir))
                if "ISTART" in incar:
                   os.system("cp %s/WAVECAR %s/WAVECAR" % (self.srcdir, self.objdir))


        if self.task.upper() == "BND":

            os.chdir(self.srcdir)
            converge = CheckConv()
            if not converge: 
                raise SystemExit("Static SCF is NOT Converged !")
            fermi_energy = GetEfermi()
            os.system("rm -rf EFERMI; echo %.3f >> EFERMI" % fermi_energy)
            os.chdir(self.wrkdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/POSCAR %s/POSCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/POTCAR %s/POTCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/EFERMI %s/EFERMI" % (self.srcdir, self.objdir))
                if lhfcalc: 
                    os.system("cp %s/IBZKPT %s/IBZKPT" % (self.srcdir, self.objdir))
                if "ICHARG" in incar:
                    os.system("cp %s/CHGCAR %s/CHGCAR" % (self.srcdir, self.objdir))
                if "ISTART" in incar:
                   os.system("cp %s/WAVECAR %s/WAVECAR" % (self.srcdir, self.objdir))

            os.chdir(self.objdir)
            
            if kpath is not None:
                MakeKpath(kpath=kpath, nkpts=nkpts)
                
            if kptfile is not None:
                with open("KPOINTS", "r") as f: flines = f.readlines()
                if nkpts is None: nkpts = 100
                flines[1] = "%g \n" % nkpts
                with open("KPOINTS", "w") as f: f.writelines(flines)
                
            if lhfcalc:
                # Reuse the irreducible BZ KPOINTS
                # from a previous standard DFT run.
                ibz_lines = open("IBZKPT").readlines()
                linemode_lines = open("KPOINTS").readlines()

                abs_path = []
                i = 4
                while i < len(linemode_lines):
                    start_kpt = linemode_lines[i].split()
                    end_kpt = linemode_lines[i+1].split()
                    increments = [
                        (float(end_kpt[0]) - float(start_kpt[0])) / nkpts,
                        (float(end_kpt[1]) - float(start_kpt[1])) / nkpts,
                        (float(end_kpt[2]) - float(start_kpt[2])) / nkpts]
                    abs_path.append(start_kpt[:3] + ["0", start_kpt[4]])
                    for n in range(1, nkpts):
                        abs_path.append(
                            [str(float(start_kpt[0]) + increments[0] * n),
                                str(float(start_kpt[1]) + increments[1] * n),
                                str(float(start_kpt[2]) + increments[2] * n), "0"])
                    abs_path.append(end_kpt[:3] + ["0", end_kpt[4]])
                    i += 3
                n_linemode_kpts = len(abs_path)

                with open("KPOINTS", "w") as kpts:
                    kpts.write("KPATH-HSE-Band \n")
                    kpts.write("{} \n".format(n_ibz_kpts + n_linemode_kpts))
                    kpts.write("Reciprocal \n")
                    for line in ibz_lines[3:]:
                        kpts.write(line)
                    for point in abs_path:
                        kpts.write("{} \n".format(" ".join(point)))
                        
            os.chdir(self.wrkdir)
        

        if self.task.upper() == "OPT":

            os.chdir(self.srcdir)
            converge = CheckConv()
            if not converge: 
                raise SystemExit("Static SCF is NOT Converged !")
            os.chdir(self.wrkdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/POSCAR %s/POSCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/POTCAR %s/POTCAR" % (self.srcdir, self.objdir))
                if "ICHARG" in incar:
                    os.system("cp %s/CHGCAR %s/CHGCAR" % (self.srcdir, self.objdir))
                if "ISTART" in incar:
                   os.system("cp %s/WAVECAR %s/WAVECAR" % (self.srcdir, self.objdir))
        
               
        if self.task.upper() == "GWA":

            os.chdir(self.srcdir)
            converge = CheckConv()
            if not converge: 
                raise SystemExit("Static SCF is NOT Converged !")
            if not os.path.isfile("WAVECAR"):
                raise SystemExit("No WAVECAR for GW Calculation !")
            if not os.path.isfile("WAVEDER"):
                raise SystemExit("No WAVEDER for GW Calculation !")
            os.chdir(self.wrkdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/POSCAR %s/POSCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/POTCAR %s/POTCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/WAVECAR %s/WAVECAR" % (self.srcdir, self.objdir))
                os.system("cp %s/WAVEDER %s/WAVEDER" % (self.srcdir, self.objdir))

        if self.task.upper() == "WAN":

            os.chdir(self.srcdir)
            if not os.path.isfile("CHGCAR") and not os.path.isfile("WAVECAR"):
                raise SystemExit("No CHGCAR or WAVECAR for WAN Calculation !")
            fermi_energy = GetEfermi()
            os.system("rm -rf EFERMI; echo %.3f >> EFERMI" % fermi_energy)
            os.chdir(self.wrkdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/POSCAR %s/POSCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/POTCAR %s/POTCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/KPOINTS %s/KPOINTS" % (self.srcdir, self.objdir))
                os.system("cp %s/EFERMI %s/EFERMI" % (self.srcdir, self.objdir))
                if icharg is not None:
                    os.system("cp %s/CHGCAR %s/CHGCAR" % (self.srcdir, self.objdir))
                if istart is not None:
                    os.system("cp %s/WAVECAR %s/WAVECAR" % (self.srcdir, self.objdir))

            os.chdir(self.objdir)
            if kpath is not None:
                MakeKpath(kpath=kpath, nkpts=nkpts)
            if wanfile is None:
                file = open("wannier90.win", "w")
                file.write("begin projections \n")
                file.write("random \n")
                file.write("end projections \n")
                file.write("guiding_centres = true \n")
                file.close()
            if os.path.isfile("kpath.win"):
                os.system("cat kpath.win >> wannier90.win")
            os.chdir(self.wrkdir)
            
        if self.task.upper() == "IPA":

            os.chdir(self.srcdir)
            if not os.path.isfile("CHGCAR") and not os.path.isfile("WAVECAR"):
                raise SystemExit("No CHGCAR or WAVECAR for IPA Calculation !")
            os.chdir(self.wrkdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/POSCAR %s/POSCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/POTCAR %s/POTCAR" % (self.srcdir, self.objdir))
                if "ICHARG" in incar:
                    os.system("cp %s/CHGCAR %s/CHGCAR" % (self.srcdir, self.objdir))
                if "ISTART" in incar:
                   os.system("cp %s/WAVECAR %s/WAVECAR" % (self.srcdir, self.objdir))

        if self.task.upper() == "BSE":

            os.chdir(self.srcdir)
            
            has_no_wfull = True
            for f in os.listdir("."):
                if f.endswith(".tmp") and f.startswith("WFULL"): 
                    has_no_wfull = False
            if has_no_wfull:
                raise SystemExit("No WFULL000*.tmp for BSE Calculation !")
            
            has_no_w000 = True
            for f in os.listdir("."):
                if f.endswith(".tmp") and f.startswith("W000"): 
                    has_no_w000 = False
            if has_no_w000:
                raise SystemExit("No W000*.tmp for BSE Calculation !")
            
            if self.objdir != self.srcdir:
                os.system("cp %s/POSCAR %s/POSCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/POTCAR %s/POTCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/WAVECAR %s/WAVECAR" % (self.srcdir, self.objdir))
                os.system("cp %s/WAVEDER %s/WAVEDER" % (self.srcdir, self.objdir))
                os.system("cp %s/*.tmp %s/" % (self.srcdir, self.objdir))

        if self.task.upper() == "PHN":

            os.chdir(self.srcdir)
            converge = CheckConv()
            if not converge: 
                raise SystemExit("Static SCF is NOT Converged !")
            magmom = GetMagmom()
            os.chdir(self.wrkdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/CONTCAR %s/POSCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/POTCAR %s/POTCAR" % (self.srcdir, self.objdir))
                if "ICHARG" in incar:
                    os.system("cp %s/CHGCAR %s/CHGCAR" % (self.srcdir, self.objdir))
                if "ISTART" in incar:
                   os.system("cp %s/WAVECAR %s/WAVECAR" % (self.srcdir, self.objdir))

            os.chdir(self.objdir)
            
            celldim = " ".join([str(s) for s in celldim])
            magmom = " ".join([str(s) for s in magmom])
            
            self.celldim = celldim
            self.magmom = magmom
            
            if magmom == "":
                os.system("phonopy -d --dim='%s'" % celldim)
            else:
                os.system("phonopy -d --dim='%s' --magmom='%s'" % (celldim, magmom))
                with open("MAGMOM", "r") as f: magmom = f.readline()
                with open("INCAR", "r") as f: incar = f.readlines()
                for i, line in enumerate(incar):
                    if "MAGMOM" in line: 
                        incar[i] = "MAGMOM = %s \n" % magmom.split("=")[1]
                with open("INCAR", "w") as f: f.writelines(incar)
            os.system("cp POSCAR POSCAR-UC")
            os.system("cp SPOSCAR POSCAR-SC")
            os.system("mv SPOSCAR POSCAR")
            
            os.chdir(self.wrkdir)

    def PhonopyPP(self, kpath="GMKG", nkpts=101, kmesh=(10,10,10), outfile="PhononBand.dat"):

        os.chdir(self.objdir)

        os.system("phonopy --fc vasprun.xml")

        with open("POSCAR-UC", "r") as f: 
            poscar = f.readlines()
            
        atomtyp = poscar[5].strip().split()
        atomtyp = " ".join([str(s) for s in atomtyp])
        
        celldim, magmom = self.celldim, self.magmom

        if self.task.upper() == "DOS":
            
            kmesh = " ".join([str(s) for s in kmesh])
            f = open("dos.conf", "w")
            f.write("ATOM_NAME = %s \n" % atomtyp)
            f.write("DIM = %s \n" % celldim)
            f.write("MP = %s \n" % kmesh)
            # f.write("DOS_RANGE = 0 10 0.001 # EMIN, EMAX, DE \n")
            # f.write("PDOS = 1 2, 3 4 5 # EMIN, EMAX, DE \n")

        if self.task.upper() == "BND":

            MakeKpath(poscar="POSCAR-UC", kpath=kpath, nkpts=nkpts, kpoints="KPATH")
            
            with open("KPATH", "r") as f: klines = f.readlines()
            
            regex = r"\s*(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*\!\s*\\*([A-Z]*)[a-z]*"

            f = open("band.conf", "w")
            f.write("ATOM_NAME = %s \n" % atomtyp)
            f.write("DIM = %s \n" % celldim)
            f.write("BAND = ")
            i = 0
            for line in klines:
                find = re.match(regex, line)
                if find != None:
                    i = i + 1
                    kx, ky, kz, s = find.groups()
                    kx, ky, kz = float(kx), float(ky), float(kz)
                    if i % 2 == 1: f.write("%.9f  %.9f  %.9f  " % (kx, ky, kz))
            if i % 2 == 0: f.write("%.9f  %.9f  %.9f  " % (kx, ky, kz))
            f.write("\n")
            # f.write("BAND_POINTS = {} \n".format(nkpts))
            f.write("NPOINTS = {} \n".format(nkpts))
            f.write("BAND_LABELS = ")
            i = 0
            for line in klines:
                find = re.match(regex, line)
                if find != None:
                    i = i + 1
                    kx, ky, kz, s = find.groups()
                    if i % 2 != 0: f.write("%s  " % s)
            if i % 2 == 0: f.write("%s  " % s)
            f.write("\n")
            f.write("EIGENVECTORS = .TRUE. \n")
            f.write("FORCE_CONSTANTS = READ \n")
            f.close()

            if magmom == "":
                os.system("phonopy --dim='%s' -c POSCAR-UC -p band.conf" % celldim)
            else:
                os.system("phonopy --dim='%s' -c POSCAR-UC -p band.conf --magmom='%s'" % (celldim, magmom))

            # os.system("phonopy-bandplot --gnuplot > %s" % outfile)
            GetPhBand(infile="band.yaml", outfile=outfile, units="THz")
        
        os.chdir(self.wrkdir)
        

    def MakeShift(self, 
                  poscar="POSCAR",
                  shift_in_abc=None, 
                  shift_in_xyz=None, 
                  shift_in_top=False, 
                  shift_in_bot=False,
                  shift_in_both=True,
                  is_direct=True, 
                  vacuum=20):
    
        """
        shift_in_abc (a tuple of three floats): shift in direct coordinates. Defaults to None.
        shift_in_xyz (a tuple of three floats): shift in cartesian coordinates. Defaults to None.
        shift_on_which (str): for which poscar lattice is shifted. Defaults to None.
        strain_on_which (str): for which poscar stran is applied. Defaults to None.
        in_direct (bool): whether to write poscar in direct coordinates. Defaults to True.
        vacuum (float): the size of vacuum to eliminate the spurious interaction. Defaults to 20.
        limitations: only support bilayer structures placed in the center of unit cell and 
        separated by an interlayer distance.
        """

        if ase_install:
            import ase.io as io
        else:
            print("Warning: ASE is not yet installed on your machine!")
            exit(1)

        os.chdir(self.objdir)
        
        os.system("cp %s %s_NO_SHIFT" % (poscar, poscar))
                
        atoms = io.read(poscar)
        cell = atoms.get_cell()[:] 
        positions = atoms.get_positions()
        cvector = cell[-1][-1]
        zcoords = positions[:, 2]
        
        # positions_top_layer = positions[zcoords > cvector/2, :]
        # positions_bot_layer = positions[zcoords < cvector/2, :]

        shift_xyz = np.zeros(3)
        if shift_in_abc:
            shift_xyz += cell.T.dot(np.asarray(shift_in_abc))
        if shift_in_xyz:
            shift_xyz += np.asarray(shift_in_xyz)
            
        if shift_in_top:
            positions[zcoords > cvector/2, 0] += shift_xyz[0]
            positions[zcoords > cvector/2, 1] += shift_xyz[1]
            positions[zcoords > cvector/2, 2] += shift_xyz[2]
        if shift_in_bot: 
            positions[zcoords < cvector/2, 0] += shift_xyz[0]
            positions[zcoords < cvector/2, 1] += shift_xyz[1]
            positions[zcoords < cvector/2, 2] += shift_xyz[2]
        if shift_in_both:
            positions[zcoords > cvector/2, 0] += shift_xyz[0]/2
            positions[zcoords > cvector/2, 1] += shift_xyz[1]/2
            positions[zcoords > cvector/2, 2] += shift_xyz[2]/2
            positions[zcoords < cvector/2, 0] -= shift_xyz[0]/2
            positions[zcoords < cvector/2, 1] -= shift_xyz[1]/2
            positions[zcoords < cvector/2, 2] -= shift_xyz[2]/2
            
        atoms.set_positions(positions)
        atoms.center(vacuum=0.5*vacuum, axis=2)
                
        if is_direct:
            atoms.write(poscar, format="vasp", vasp5=True, direct=True)
        else:
            atoms.write(poscar, format="vasp", vasp5=True, direct=False)
        
        os.chdir(self.wrkdir)
        
        
    def ApplyStrain(self, poscar="POSCAR", epsilon=0):
        
        if ase_install:
            import ase.io as io
        else:
            print("Warning: ASE is not yet installed on your machine!")
            exit(1)

        os.chdir(self.objdir)
        
        os.system("cp %s %s_NO_STRAIN" % (poscar, poscar))
        
        atoms = io.read(poscar)   
        cell = atoms.get_cell()[:]
        a, b, c = np.sqrt(np.sum(cell**2, axis=1))
        ab = (a + b) / 2 * (1 + epsilon)
        avec = (ab, 0, 0)
        bvec = (ab*(-1/2), ab*np.sqrt(3)/2, 0) 
        cvec = (0, 0, c)
        cell[0] = np.asarray(avec)
        cell[1] = np.asarray(bvec)
        cell[2] = np.asarray(cvec)
        atoms.set_cell(cell, scale_atoms=True)
        atoms.center(vacuum=10, axis=2)
        atoms.write(poscar, format="vasp", vasp5=True, direct=True)
    
        os.chdir(self.wrkdir)

    def Submission(self, 
                   queue=None, 
                   jobname=None, 
                   ntasks=None, 
                   lprocs=None, 
                   nnodes=None, 
                   ncores=None,
                   walltime=None,
                   prefix=None):

        """
        Submission script to run calculations on LOCAL/PBS/SLURM systems.
        Arguments:
        queue = "LOCAL" or "SLURM" or "PBS"
        jobname = Name of a submitted job: e.g. "relaxation"
        ntasks = Number of MPI processes requested: e.g. 16 
        walltime = Maximally allowed CPU time: "00:30:00"
        """

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
            
            if self.lsoc:
                runjob.write("%s %g vasp_ncl \n" % (prefix, ntasks))
            else:
                runjob.write("%s %g vasp_std \n" % (prefix, ntasks))
                    
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
            
            if self.lsoc:
                runjob.write("srun -n %g vasp_ncl \n" % ntasks)
            else:
                runjob.write("srun -n %g vasp_std \n" % ntasks)
            
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
            # runjob.write('#PBS -o {}.out \n'.format(jobname))
            # runjob.write('#PBS -e {}.err \n'.format(jobname))
            
            runjob.write("\n")
            runjob.write("cd $PBS_O_WORKDIR \n")
            runjob.write("\n")

            runjob.write("module load --purge \n")
            # runjob.write("module load calcua/2020a \n")
            # runjob.write("module load intel/2020a \n")
            # runjob.write("export PATH=$VSC_DATA/VASP2WANN/bin:$PATH \n")
            runjob.write("module load intel/2018a \n")
            runjob.write("module load VASP/5.4.4-intel-2018a-Wannier90-1.2 \n")
            
            runjob.write("\n")

            if self.lsoc:
                runjob.write("mpirun -np %g vasp_ncl \n" % ntasks)
            else:
                runjob.write("mpirun -np %g vasp_std \n" % ntasks)

            runjob.close()

        if queue.upper() == "LOCAL": 
            os.system("bash runjob")
            
        if queue.upper() == "SLURM": 
            os.system("sbatch runjob")
            
        if queue.upper() == "PBS": 
            os.system("qsub runjob")  

        os.chdir(self.wrkdir)


class TestConv(PythVasp):

    def __init__(self, task="SCF", topdir=None, wrkdir=None, srcdir=None, objdir=None):

        os.chdir(topdir)
        
        if wrkdir is None:
            wrkdir = os.getcwd()
        if srcdir is None:
            print("No INPUT Directory for VASP Calculation !")
            exit(1)
        if objdir is None:
            print("No OUTPUT Directory for VASP Calculation !")
            exit(1)    
                    
        self.task = task
        self.topdir = topdir
        self.wrkdir = os.path.abspath(wrkdir)
        self.srcdir = os.path.join(self.wrkdir, srcdir)
        self.objdir = os.path.join(self.wrkdir, objdir)

        if not os.path.exists(self.srcdir): 
            os.mkdir(self.srcdir)
        
        if not os.path.exists(self.objdir): 
            os.mkdir(self.objdir)
            
    def MakeInput(self, **kwargs):

        input_dict = dict()

        for k, v in kwargs.items():
            input_dict[k] = v

        input_dict.update(nsw=0)
        input_dict.update(ibrion=-1)
        input_dict.update(ediffg="-1E-2")
        input_dict.update(ediff="1E-6")

        self.input_dict = input_dict


    def Scheduler(self, queue="LOCAL", jobname=None, ntasks=None, walltime=None):

        self.queue = queue
        self.jobname = jobname
        self.ntasks = ntasks
        self.walltime = walltime

    def LoopParam(self, param_name="encut", param_list=[100,200,300]):

        TASK = self.task
        TOPDIR = self.topdir
        WRKDIR = self.wrkdir
        SRCDIR = self.srcdir
        OBJDIR = self.objdir
        input_dict = self.input_dict

        for param in param_list:
            objdir = OBJDIR + "/" + "{}".format(param)
            if "nkpts" in param_name:
                param_str = "x".join([str(s) for s in param])
                objdir = OBJDIR + "/" + "{}".format(param_str)
            if "encut" in param_name: input_dict.update(encut=param)
            if "nkpts" in param_name: input_dict.update(nkpts=param)
            if "sigma" in param_name: input_dict.update(sigma=param)
            if "npar" in param_name: input_dict.update(npar=param)
            if "kpar" in param_name: input_dict.update(kpar=param)
            super().__init__(task=TASK, topdir=TOPDIR, wrkdir=WRKDIR, srcdir=SRCDIR, objdir=objdir)
            super().MakeInput(**input_dict)
            super().Submission(queue=self.queue, jobname=self.jobname, ntasks=self.ntasks, walltime=self.walltime)

        # confirm directories 
        self.wrkdir = WRKDIR
        self.srcdir = SRCDIR
        self.objdir = OBJDIR


    def CollectData(self, param_name="encut", param_list=[100,200,300]):

        WRKDIR = self.wrkdir
        OBJDIR = self.objdir

        energy_list = []
        nkpts_list = []
        cputime_list = []
        for param in param_list:
            objdir = OBJDIR + "/" + "{}".format(param)
            if "nkpts" in param_name:
                param_str = "x".join([str(s) for s in param])
                objdir = OBJDIR + "/" + "{}".format(param_str)
            os.chdir(objdir)
            conv = CheckConv(infile="OUTCAR")
            if conv["scf"]:
                energy = GetEnergy(infile="OUTCAR")
                nkpts = GetNkpts(infile="OUTCAR")
                cputime = GetCPUTime(infile="OUTCAR")
                energy_list.append(energy)
                nkpts_list.append(nkpts)
                cputime_list.append(cputime)
            os.chdir(WRKDIR)

        if "nkpts" in param_name: 
            param_list = nkpts_list

        os.chdir(OBJDIR)
        with open("summary.dat", "w") as f:
            f.write("#%s   TotalEnergy (eV)   NumberOfKpoints   CPURunTime (h)" % param_name.upper())
            for param, energy, nkpts, cputime in zip(param_list, energy_list, nkpts_list, cputime_list):
                f.write("{}   {}   {}   {} \n".format(param, energy, nkpts, cputime/3600))
        with open("plot_data.py", "w") as f:
            f.write("import numpy as np \n")
            f.write("import matplotlib.pyplot as plt \n")
            f.write("\n")
            f.write("data = np.loadtxt('summary.dat') \n")
            f.write("plt.plot(data[:, 0], data[:, 1], 'bo-', lw=2, ms=10) \n")
            f.write("plt.plot(data[:, 0], data[:, 2], 'bo-', lw=2, ms=10) \n")
            f.write("plt.plot(data[:, 0], data[:, 3], 'bo-', lw=2, ms=10) \n")
            f.write("plt.show() \n")
        os.chdir(WRKDIR)


#if __name__ == "__main__":

    # import PythVasp
    
    # PWD = os.getcwd()

    # Job = PythVasp(task="RLX", topdir=PWD, wrkdir="MoS2", srcdir="INP", objdir="RLX")
    # Job.MakeInput(encut=500, nsw=100, ibrion=2, isif=2, ediffg="-1E-2", algo="Normal",
    # nelm=200, ediff="1E-6", ismear=0, sigma=0.05, nkpts=(12,12,12), gamma=True)
    # Job.Submission(queue="local", ncores=40) # queue = local/slurm/pbs

    # Job = PythVasp(task="PHN", topdir=PWD, wrkdir="MoS2", srcdir="RLX", objdir="PHN")
    # Job.MakeInput(encut=500, ibrion=8, ialgo=38, ediff="1E-8", ismear=0, sigma=0.05,
    # nkpts=(12,12,12), gamma=True, lreal=False, addgrid=True, celldim="221")
    # Job.Submission(queue="local", ncores=40) # queue = local/slurm/pbs
