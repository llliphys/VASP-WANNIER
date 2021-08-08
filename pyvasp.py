import os
import re
import numpy as np


# check pymatgen installed or not
try:
    import pymatgen
    pymatgen_install = True
except ImportError:
    pymatgen_install = False


def GetEnergy(infile="OSZICAR"):

    with open(infile, "r") as f:
        flines = f.readlines()

    energy = float(flines[-1].strip().split()[2])

    return energy

def GetNkpts(infile="OUTCAR"):

    with open(infile, "r") as f:
        fread = f.read()

    regex = r"\s*k-points\s+(NKPTS\s*=\s+\d+)"
    match = re.findall(regex, fread)[-1]
    nkpts = int(match.split("=")[1])

    return nkpts

def GetCPUTime(infile="OUTCAR"):

    with open(infile, "r") as f:
        flines = f.readlines()

    for line in flines[::-1]:
        if "Elapsed time (sec)" in line:
            cputime = float(line.strip().split(":")[1])
            break

    return cputime


def GetMagmom(infile="OUTCAR"):

    with open(infile, "r") as f:
        outcar = f.readlines()

    natoms = None
    for i, line in enumerate(outcar):
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

    outcar = open(infile, "r").read()
    match = re.findall(r"E-fermi\s*:\s*(-?\d+.\d+)", outcar)[-1]
    Efermi = float(match)

    return efermi

def GetNbands(infile="OUTCAR"):

    outcar = open(infile, "r").read()
    match = re.findall(r"NBANDS\s*=\s*(\d+)", outcar)[-1]
    nbands = int(match)

    return nbands


def MakeClean(chgcar=False, wavecar=False):

    inputs = ["INCAR", "POSCAR", "POTCAR", "KPOINTS", "OPTCELL"]
    if chgcar: inputs = ["INCAR", "POSCAR", "POTCAR", "KPOINTS", "OPTCELL", "CHGCAR"]
    if wavecar: inputs = ["INCAR", "POSCAR", "POTCAR", "KPOINTS", "OPTCELL", "WAVECAR"]
    if chgcar and wavecar: inputs = ["INCAR", "POSCAR", "POTCAR", "KPOINTS", "OPTCELL", "CHGCAR", "WAVECAR"]

    for f in os.listdir("."):
        if f not in inputs: os.system("rm -rf %s" % f)


def MakeKmesh2D(kx1=0, kx2=1, nkx=101, ky1=0, ky2=1, nky=101, kz=0):

    """
    This function generates a 2D Kmesh for e.g. a Fermi surface calculation.
    """

    kxs = np.linspace(kx1, kx2, nkx)
    kys = np.linspace(ky1, ky2, nky)

    kmesh = open("KPOINTS", "w")
    kmesh.write("Nkx-by-Nky-Kmesh\n")
    kmesh.write("%d\n" % (nkx * nky))
    kmesh.write("Cartesian\n")
    for kx in kxs:
        for ky in kys:
            kmesh.write("%.12f  %.12f  %.12f  %.1f\n" % (kx, ky, kz, 1.0))
    kmesh.close()


def MakeKpath(poscar=None, kpath=None, nkpts=None, kpoints=None):

    if pymatgen_install:
        from pymatgen.io.vasp.inputs import Kpoints
        from pymatgen.core import Structure
        from pymatgen.symmetry.bandstructure import HighSymmKpath
    else:
        print("Warning: Pymatgen is not yet installed on your machine!")
        print("Either install it using conda or provide KPOINTS by hand!")
        exit(1)

    if poscar is None: poscar = "POSCAR"
    if nkpts is None: nkpts = 100
    if kpoints is None: kpoints = "KPOINTS"

    if kpath is not None:
        kpath = [[str(s) for s in kpath]]
        ibz = HighSymmKpath(Structure.from_file(poscar))
        ibz.kpath['path'] = kpath
        Kpoints.automatic_linemode(nkpts, ibz).write_file(kpoints)
    else:
        ibz = HighSymmKpath(Structure.from_file(poscar))
        Kpoints.automatic_linemode(nkpts, ibz).write_file(kpoints)



class PythVasp:

    """
    - A Python class for carrying out a series of DFT-VASP calculations including e.g.
    "Relaxation", "Static SCF", "Band Structure", "Density of States", "Charge Density"
    "Fermi Surface", "Hybrid Functional", "GW Calculation", "BSE calculation", which are 
    reprensted by task="RLX", "SCF", "BND", "DOS", "CHG", "FSF", "HSE", "GWA", and "BSE", 
    respectively. This class module can be used as an automative workflow to implement
    creating inputs from sratch and restarting from existing calculations to submitting
    calculations to popular job manager systems such as PBS and SLURM.

    - Before doing any calculations, one has always to specify two working directories,
    one of which contains the input files for VASP calculations and the other contains
    the output files from the calculations. The two directories for inputs and outputs
    are labeled by scrdir and objdir, respectively. Note that scrdir and objdir can be
    the same or different directories, depending on specific calculations.

    - A rootdir can also be specified to indicate the root directory that contains the two
    subdirectories inside it, i.e., rootdir/srcdir and rootdir/objdir.

    - A child class based on this parent class is implemented for testing the convergence
    against various important DFT parameters (e.g. ENCUT, SIGMA, KPOINTS).

    """

    def __init__(self, task="RLX", rootdir=None, srcdir=None, objdir=None):

        self.task = task
        self.srcdir = srcdir
        self.objdir = objdir
        self.rootdir = rootdir

        if self.rootdir is None: self.rootdir = os.getcwd()

        if self.srcdir is None or self.objdir is None:
            raise SystemExit("Error: No INPUTS for VASP Calculation !")
        else:
            if not os.path.exists(self.srcdir): os.mkdir(self.srcdir)
            if not os.path.exists(self.objdir): os.mkdir(self.objdir)

    def MakeInput(self, encut=None, prec=None, symprec=None, \
                  nsw=None, isif=None, ibrion=None, potim=None, ediffg=None, \
                  algo=None, ialgo=None, nelm=None, nelmin=None, ediff=None, \
                  icharg=None, istart=None, lorbit=None, \
                  ismear=None, sigma=None, nedos=None, \
                  nelect=None, \
                  ivdw=None, isym=None, \
                  ispin=None, magmom=None, nupdown=None, \
                  lsoc=False, saxis=None, lorbmom=None, \
                  lwave=None, lcharg=None, \
                  lhfcalc=None, hfscreen=None, aexx=None, precfock=None, \
                  lvtot=None, lvhar=None, \
                  lreal=None, addgrid=None, \
                  ldau=None, utype=None, ldaul=None, ldauu=None, ldauj=None, \
                  amix=None, bmix=None, amix_mag=None, bmix_mag=None, \
                  nbands=None, loptics=None, cshift=None, omegamax=None, \
                  encutgw=None, nomega=None, nbandso=None, nbandsv=None, \
                  efield=None, idipol=None, ldipol=None, dipol=None, \
                  lpard=None, nbmod=None, eint=None, iband=None, kpuse=None, \
                  nkpts=None, gamma=True, shift=None, kpath=None, kptdir=None, \
                  poscar=None, potcar=None, optcell=None, \
                  npar=None, kpar=None, ncore=None, \
                  celldim=None, \
                  make_clean=None, \
                  make_2dkmesh=None):

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
        if lsoc: incar.update(ISYM=0, LSORBIT=".TRUE.", LORBMOM=".TRUE.")
        if ispin is not None: incar.update(ISPIN=ispin)
        if magmom is not None: incar.update(MAGMOM=magmom)
        if nupdown is not None: incar.update(NUPDOWN=nupdown)
        if saxis is not None: incar.update(SAXIS=saxis)
        if nbands is not None: incar.update(NBANDS=nbands)
        if encutgw is not None: incar.update(ENCUTGW=encutgw)
        if nomega is not None: incar.update(NOMEGA=nomega)
        if nbandso is not None: incar.update(NBANDSO=nbandso)
        if nbandsv is not None: incar.update(NBANDSV=nbandsv)
        if efield is not None: incar.update(EFIELD=efield)
        if idipol is not None: incar.update(IDIPOL=idipol)
        if ldipol is not None: incar.update(LDIPOL=ldipol)
        if dipol is not None: incar.update(DIPOL=dipol)
        if lpard: incar.update(LPARD=".TRUE.")
        if nbmod is not None: incar.update(NBMOD=nbmod)
        if eint is not None: incar.update(EINT=eint)
        if iband is not None: incar.update(IBAND=iband)
        if kpuse is not None: incar.update(KPUSE=kpuse)
        if kpar is not None: incar.update(KPAR=kpar)
        if npar is not None: incar.update(NPAR=npar)
        if ncore is not None: incar.update(NCORE=ncore)

        chgcar = wavecar = False
        if "ICHARG" in incar: chgcar = True
        if "ISTART" in incar: wavecar = True

        if nbands is None and self.task.upper() == "CON":
            os.chdir(self.srcdir)
            nbands = GetNbands()
            os.chdir(self.rootdir)
            if not lsoc: incar.update(NBANDS=nbands)
            if lsoc: incar.update(NBANDS=int(2*nbands))

        if nbands is None and self.task.upper() == "SCF":
            os.chdir(self.srcdir)
            nbands = GetNbands()
            os.chdir(self.rootdir)
            if not lsoc: incar.update(NBANDS=nbands)
            if lsoc: incar.update(NBANDS=int(2*nbands))

        # if magmom is None and self.task.upper() == "SCF":
        #    os.chdir(self.srcdir)
        #    magmom = GetMagmom()
        #    os.chdir(self.rootdir)
        #    if magmom != []:
        #        if not lsoc: magmom = "  ".join([str(m) for m in magmom])
        #        if lsoc: magmom = "  ".join(["0  0  %.3f" % m for m in magmom])
        #        incar.update(MAGMOM=magmom)

        os.chdir(self.objdir)

        with open("INCAR", "w") as f:
            for key, val in incar.items():
                f.write("{} = {} \n".format(key, val))

        if optcell is not None:
            # example format: optcell = [[1, 1, 1], [1, 1, 1], [1, 1, 1]]
            with open("OPTCELL", "w") as f:
                for ijk in optcell:
                    f.write("{}{}{} \n".format(ijk[0], ijk[1], ijk[2]))

        with open("KPOINTS", "w") as f:
            f.write("K-Mesh \n")
            f.write("0 \n")
            f.write("Gamma \n") if gamma else f.write("MP \n")
            if type(nkpts) == tuple or type(nkpts) == list:
                f.write("{} {} {} \n".format(nkpts[0], nkpts[1], nkpts[2]))
            if type(nkpts) == int:
                f.write("{} {} {} \n".format(nkpts, nkpts, nkpts))
            if shift is None:
                f.write("0 0 0 \n")
            else:
                f.write("{} {} {} \n".format(shift[0], shift[1], shift[2]))

        os.chdir(self.rootdir)

        if self.task.upper() == "RLX":

            if poscar is None: poscar = "POSCAR"
            if potcar is None: potcar = "POTCAR"

            os.chdir(self.srcdir)
            if not os.path.isfile(poscar) or not os.path.isfile(potcar):
                raise SystemExit("Error: No POSCAR and/or POTCAR for DFT Relaxation !")
            os.chdir(self.rootdir)

            os.chdir(self.objdir)
            if make_clean: MakeClean(chgcar=chgcar, wavecar=wavecar)
            os.chdir(self.rootdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/%s %s/POSCAR" % (self.srcdir, poscar, self.objdir))
                os.system("cp %s/%s %s/POTCAR" % (self.srcdir, potcar, self.objdir))

        if self.task.upper() == "CON":

            os.chdir(self.srcdir)
            if not os.path.isfile("CONTCAR"):
                raise SystemExit("Error: No CONTCAR for DFT Continuation !")
            os.chdir(self.rootdir)

            os.system("cp %s/CONTCAR %s/POSCAR" % (self.srcdir, self.objdir))

            os.chdir(self.objdir)
            if make_clean: MakeClean(chgcar=chgcar, wavecar=wavecar)
            os.chdir(self.rootdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/POTCAR %s/POTCAR" % (self.srcdir, self.objdir))
                if "ICHARG" in incar:
                    os.system("cp %s/CHGCAR %s/CHGCAR" % (self.srcdir, self.objdir))
                if "ISTART" in incar:
                   os.system("cp %s/WAVECAR %s/WAVECAR" % (self.srcdir, self.objdir))

        if self.task.upper() == "SCF":

            os.chdir(self.srcdir)
            if not os.path.isfile("CONTCAR"):
                raise SystemExit("Error: No CONTCAR for SCF Calculation !")
            os.chdir(self.rootdir)

            os.system("cp %s/CONTCAR %s/POSCAR" % (self.srcdir, self.objdir))

            os.chdir(self.objdir)
            if make_clean: MakeClean(chgcar=chgcar, wavecar=wavecar)
            os.chdir(self.rootdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/POTCAR %s/POTCAR" % (self.srcdir, self.objdir))
                if "ICHARG" in incar:
                    os.system("cp %s/CHGCAR %s/CHGCAR" % (self.srcdir, self.objdir))
                if "ISTART" in incar:
                   os.system("cp %s/WAVECAR %s/WAVECAR" % (self.srcdir, self.objdir))

        if self.task.upper() == "CHG":

            os.chdir(self.srcdir)
            if not os.path.isfile("CONTCAR"):
                raise SystemExit("Error: No CONTCAR for CHG Calculation !")
            os.chdir(self.rootdir)

            os.chdir(self.objdir)
            if make_clean: MakeClean(chgcar=chgcar, wavecar=wavecar)
            os.chdir(self.rootdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/POSCAR %s/POSCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/POTCAR %s/POTCAR" % (self.srcdir, self.objdir))
                if "ICHARG" in incar:
                    os.system("cp %s/CHGCAR %s/CHGCAR" % (self.srcdir, self.objdir))
                if "ISTART" in incar:
                   os.system("cp %s/WAVECAR %s/WAVECAR" % (self.srcdir, self.objdir))

        if self.task.upper() == "FSF":

            os.chdir(self.srcdir)
            if not os.path.isfile("CONTCAR"):
                raise SystemExit("Error: No CONTCAR for FSF Calculation !")
            os.chdir(self.rootdir)

            os.chdir(self.objdir)
            if make_clean: MakeClean(chgcar=chgcar, wavecar=wavecar)
            if make_2dkmesh: MakeKmesh2D(kx1=-1, kx2=1, nkx=51, ky1=-1, ky2=1, nky=51, kz=0)
            os.chdir(self.rootdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/POSCAR %s/POSCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/POTCAR %s/POTCAR" % (self.srcdir, self.objdir))
                if "ICHARG" in incar:
                    os.system("cp %s/CHGCAR %s/CHGCAR" % (self.srcdir, self.objdir))
                if "ISTART" in incar:
                   os.system("cp %s/WAVECAR %s/WAVECAR" % (self.srcdir, self.objdir))

        if self.task.upper() == "BND" or self.task.upper() == "DOS":

            os.chdir(self.srcdir)
            if not os.path.isfile("CHGCAR"):
                raise SystemExit("Error: No CHGCAR for BND or DOS Calculation !")
            os.chdir(self.rootdir)

            os.chdir(self.objdir)
            if make_clean: MakeClean(chgcar=chgcar, wavecar=wavecar)
            os.chdir(self.rootdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/POSCAR %s/POSCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/POTCAR %s/POTCAR" % (self.srcdir, self.objdir))
                if lhfcalc: os.system("cp %s/IBZKPT %s/IBZKPT" % (self.srcdir, self.objdir))
                if "ICHARG" in incar:
                    os.system("cp %s/CHGCAR %s/CHGCAR" % (self.srcdir, self.objdir))
                if "ISTART" in incar:
                   os.system("cp %s/WAVECAR %s/WAVECAR" % (self.srcdir, self.objdir))

            os.chdir(self.objdir)
            if self.task.upper() == "BND":
                if kpath is not None:
                    MakeKpath(kpath=kpath, nkpts=nkpts)
                if kptdir is not None:
                    os.system("cp %s/%s/KPOINTS KPOINTS" % (self.rootdir, kptdir))
                    with open("KPOINTS", "r") as f:
                        flines = f.readlines()
                    flines[1] = "%g \n" % nkpts
                    with open("KPOINTS", "w") as f:
                        f.writelines(flines)
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
                        kpts.write("K-Mesh For HSE Band \n")
                        kpts.write("{} \n".format(n_ibz_kpts + n_linemode_kpts))
                        kpts.write("Reciprocal \n")
                        for line in ibz_lines[3:]:
                            kpts.write(line)
                        for point in abs_path:
                            kpts.write("{} \n".format(" ".join(point)))
            os.chdir(self.rootdir)

        if self.task.upper() == "GWA":

            os.chdir(self.srcdir)
            if not os.path.isfile("WAVECAR") or not os.path.isfile("WAVEDER"):
                raise SystemExit("Error: No WAVECAR and/or WAVEDER for GW Calculation !")
            os.chdir(self.rootdir)

            os.chdir(self.objdir)
            if make_clean: MakeClean(chgcar=chgcar, wavecar=wavecar)
            os.chdir(self.rootdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/POSCAR %s/POSCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/POTCAR %s/POTCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/WAVECAR %s/WAVECAR" % (self.srcdir, self.objdir))
                os.system("cp %s/WAVEDER %s/WAVEDER" % (self.srcdir, self.objdir))

        if self.task.upper() == "IPA":

            os.chdir(self.srcdir)
            if not os.path.isfile("CHGCAR") or not os.path.isfile("WAVECAR"):
                raise SystemExit("Error: No CHGCAR and/or WAVECAR for Optical Calculation !")
            os.chdir(self.rootdir)

            os.chdir(self.objdir)
            if make_clean: MakeClean(chgcar=chgcar, wavecar=wavecar)
            os.chdir(self.rootdir)

            if self.objdir != self.srcdir:
                os.system("cp %s/POSCAR %s/POSCAR" % (self.srcdir, self.objdir))
                os.system("cp %s/POTCAR %s/POTCAR" % (self.srcdir, self.objdir))
                if "ICHARG" in incar:
                    os.system("cp %s/CHGCAR %s/CHGCAR" % (self.srcdir, self.objdir))
                if "ISTART" in incar:
                   os.system("cp %s/WAVECAR %s/WAVECAR" % (self.srcdir, self.objdir))

        if self.task.upper() == "BSE":

            os.chdir(self.srcdir)
            has_no_wfile = True
            for f in os.listdir("."):
                if f.endswith(".tmp"): has_no_wfile = False
            if has_no_wfile:
                raise SystemExit("Error: No WFULL000*.tmp and W000*.tmp for BSE Calculation !")
            os.chdir(self.rootdir)

            if make_clean: os.system("rm -rf %s" % self.objdir)

            os.system("cp -r %s %s" % (self.srcdir, self.objdir))

        if self.task.upper() == "PHN":

            os.chdir(self.srcdir)
            if not os.path.isfile("CONTCAR"):
                raise SystemExit("Error: No CONTCAR for PHONON Calculation !")
            magmom = GetMagmom(infile="OUTCAR")
            os.chdir(self.rootdir)

            os.chdir(self.objdir)
            if make_clean: MakeClean(chgcar=chgcar, wavecar=wavecar)
            os.chdir(self.rootdir)

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
            self.magmom = magmom
            if magmom == "":
                os.system("phonopy -d --dim='%s'" % celldim)
            else:
                os.system("phonopy -d --dim='%s' --magmom='%s'" % (celldim, magmom))
                with open("MAGMOM", "r") as f: magmom = f.readline()
                with open("INCAR", "r") as f: incar = f.readlines()
                for i, line in enumerate(incar):
                    if "MAGMOM" in line: incar[i] = "MAGMOM = %s \n" % magmom.split("=")[1]
                with open("INCAR", "w") as f: f.writelines(incar)
            os.system("cp POSCAR POSCAR-UC")
            os.system("cp SPOSCAR POSCAR-SC")
            os.system("mv SPOSCAR POSCAR")
            os.chdir(self.rootdir)

    def PhonopyPP(self, celldim=None, kpath=None, nkpts=None, outfile="phonon_band.dat"):

        os.chdir(self.objdir)

        os.system("phonopy --fc vasprun.xml")

        with open("POSCAR-UC", "r") as f: poscar = f.readlines()
        atomtyp = poscar[5].strip().split()

        atomtyp = " ".join([str(s) for s in atomtyp])
        celldim = " ".join([str(s) for s in celldim])
        magmom = self.magmom

        if self.task.upper() == "BND":

            MakeKpath(poscar="POSCAR-UC", kpath=kpath, kpoints="KPATH")
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
                    if i % 2 != 0: f.write("%.9f  %.9f  %.9f  " % (kx, ky, kz))
            if i % 2 == 0: f.write("%.9f  %.9f  %.9f  " % (kx, ky, kz))
            f.write("\n")
            f.write("BAND_POINTS = {} \n".format(nkpts))
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
            f.write("FORCE_CONSTANTS = READ \n")
            f.close()

            if magmom == "":
                os.system("phonopy --dim='%s' -c POSCAR-UC -p band.conf" % celldim)
            else:
                os.system("phonopy --dim='%s' -c POSCAR-UC -p band.conf --magmom='%s'" % (celldim, magmom))

            os.system("phonopy-bandplot --gnuplot > %s" % outfile)

        os.chdir(self.rootdir)


    def Submission(self, queue=None, jobname=None, \
                   nnodes=None, ncores=None, walltime=None, \
                   bindir=None, vasp5=None, vasp6=None):

        """
        Submission script to run VASP calculations on either LOCAL or SLURM systems.
        Arguments:
        queue = "Local" or "Slurm" (will include PBS schedule)
        jobname = Name of a submitted job: "slab relaxation"
        nndoes = Number of computer nodes requested: 2 (nodes)
        ncores = Number of cores per node requested: 32 (cores)
        walltime = Mximally allowed CPU time: "00:30:00"
        bindir = User specified VASP bin directory: /usr/bin
        """

        self.queue = queue

        os.chdir(self.objdir)

        if nnodes is None: nnodes = 1
        if ncores is None: ncores = mp.cpu_count()

        if queue.upper() == "LOCAL":
            runjob = open('runjob', 'w')
            runjob.write('#!/bin/bash \n')
            runjob.write("\n\n")
            #runjob.write("export PATH=$HOME/ANACONDA/bin:$PATH \n")
            # if vasp5: runjob.write("export PATH=$HOME/VASP.5.4.4/bin:$PATH \n")
            # if vasp6: runjob.write("export PATH=$HOME/VASP.6.1.1/bin:$PATH \n")
            # runjob.write("\n\n")
            # runjob.write("module unload compiler mpi \n")
            # runjob.write("module load compiler mpi/openmpi/4.0.2 \n")
            # runjob.write("\n\n")
            if bindir is None:
                runjob.write("export PATH=$HOME/VASP2WANN/bin:$PATH \n")
                runjob.write("\n\n")
                if self.lsoc:
                    runjob.write('mpirun -np %g vasp_ncl \n' % ncores)
                else:
                    runjob.write('mpirun -np %g vasp_std \n' % ncores)
            if bindir is not None:
                runjob.write('mpirun -np %g %s \n' % (ncores, bindir))
            runjob.close()

        if queue.upper() == "SLURM":
            ntasks = int(nnodes * ncores)
            runjob = open('runjob', 'w')
            runjob.write('#!/bin/bash \n')
            runjob.write('#SBATCH --job-name={} \n'.format(jobname))
            runjob.write('#SBATCH -o {}.out \n'.format(jobname))
            runjob.write('#SBATCH -e {}.err \n'.format(jobname))
            runjob.write('#SBATCH --nodes={} \n'.format(nnodes))
            runjob.write('#SBATCH --ntasks={} \n'.format(ntasks))
            runjob.write('#SBATCH --time {} \n'.format(walltime))
            runjob.write("\n\n")
            #runjob.write("export PATH=$HOME/ANACONDA/bin:$PATH \n")
            runjob.write("export PATH=$HOME/VASP2WANN/bin:$PATH \n")
            runjob.write("export PATH=$HOME/WANNIER90/V3.1:$PATH \n")
            runjob.write("\n\n")
            runjon.write("module load slurm_setup \n")
            if vasp5: runjob.write("module load vasp \n")
            if vasp6: runjob.write("module load vasp/6.1.2 \n")
            runjob.write("\n\n")
            if LSORBIT == ".TRUE.":
                if vasp5: runjob.write("vasp5 -n $SLURM_NTASKS -s full \n")
                if vasp6: runjob.write("vasp6 -n $SLURM_NTASKS -s ncl \n")
            else:
                if vasp5: runjob.write("vasp5 -n $SLURM_NTASKS \n")
                if vasp6: runjob.write("vasp6 -n $SLURM_NTASKS \n")
            # runjob.write('cd $SLURM_SUBMIT_DIR \n')
            # runjob.write('module load intel/2016.0.109 \n')
            # runjob.write('module load openmpi/1.10.1 \n')
            # runjob.write('module load vasp/5.4.1 \n')
            # runjob.write('mpirun {} \n'.format(binary))
            runjob.close()

        if self.queue.upper() == "LOCAL": os.system("bash runjob")
        if self.queue.upper() == "SLURM": os.system("sbatch runjob")

        os.chdir(self.rootdir)


class TestConv(PythVasp):

    def __init__(self, task="CON", rootdir=None, srcdir=None, objdir=None):

        self.task = task
        self.srcdir = srcdir
        self.objdir = objdir
        self.rootdir = rootdir

        if rootdir is None: self.rootdir = os.getcwd()

        if self.srcdir is None or self.objdir is None:
            raise SystemExit("Error: No INPUTS for VASP Calculation !")
        else:
            if not os.path.exists(self.srcdir): os.mkdir(self.srcdir)
            if not os.path.exists(self.objdir): os.mkdir(self.objdir)

    def MakeInput(self, **kwargs):

        input_dict = dict()

        for k, v in kwargs.items():
            input_dict[k] = v

        input_dict.update(nsw=0)
        input_dict.update(ibrion=-1)
        input_dict.update(ediffg="-1E-2")
        input_dict.update(ediff="1E-6")

        self.input_dict = input_dict

    def Scheduler(self, queue="LOCAL", jobname=None, \
                  nnodes=None, ncores=None, walltime=None, \
                  vasp5=None, vasp6=None, bindir=None):

        self.queue = queue
        self.jobname = jobname
        self.nnodes = nnodes
        self.ncores = ncores
        self.walltime = walltime
        self.vasp5 = vasp5
        self.vasp6 = vasp6
        self.bindir = bindir

    def LoopParam(self, param_name="encut", param_list=[100,200,300]):

        TASK = self.task
        ROOTDIR = self.rootdir
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
            super().__init__(task=TASK, rootdir=ROOTDIR, srcdir=SRCDIR, objdir=objdir)
            super().MakeInput(**input_dict)
            super().Submission(queue=self.queue, jobname=self.jobname, \
                               nnodes=self.nnodes, ncores=self.ncores, walltime=self.walltime, \
                               vasp5=self.vasp5, vasp6=self.vasp6, bindir=self.bindir)

        self.rootdir = ROOTDIR
        self.srcdir = SRCDIR
        self.objdir = OBJDIR


    def CollectData(self, param_name="encut", param_list=[100,200,300], from_file="OSZICAR"):

        ROOTDIR = self.rootdir
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
            energy = GetEnergy(infile="OUTCAR")
            energy_list.append(energy)
            cputime = GetCPUTime(infile="OUTCAR")
            cputime_list.append(cputime)
            nkpts = GetNkpts(infile="OUTCAR")
            nkpts_list.append(nkpts)
            os.chdir(ROOTDIR)

        if "nkpts" in param_name: param_list = nkpts_list

        os.chdir(OBJDIR)
        with open("summary.dat", "w") as f:
            for param, energy, cputime in zip(param_list, energy_list, cputime_list):
                f.write("{}  {} {} \n".format(param, energy, cputime))
        with open("plot_data.py", "w") as f:
            f.write("import numpy as np \n")
            f.write("import matplotlib.pyplot as plt \n")
            f.write("\n")
            f.write("data = np.loadtxt('summary.dat') \n")
            f.write("plt.plot(data[:, 0], data[:, 1], 'bo-', lw=2, ms=10) \n")
            f.write("plt.plot(data[:, 0], data[:, 2], 'bo-', lw=2, ms=10) \n")
            f.write("plt.show() \n")
        os.chdir(ROOTDIR)
