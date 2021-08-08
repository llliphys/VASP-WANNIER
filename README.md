### Intro
Here are some python scripts used to check the VASP input/output files: INCAR, POSCAR, POTCAR, KPOINTS, OSZICAR, OUTCAR, and so on. To check out how to use them, just simply command `python script_name.py` or `python script_name.py -h` on your terminal to have a first glance.

- `pycheck.py`: simply merge the above two scripts into a single script (more compact).

- `pycell.py`: manipulate CIF or POSCAR file, e.g., making supercell, varying vacuum size, etc.

- `workfunc.py`: calculate the work function and plot the planar averaged potental of a slab system.

- `pyvasp.py`: A Python class for carrying out a series of DFT-VASP calculations including e.g. "Relaxation", "Static SCF", "Band Structure", "Density of States", "Charge Density", "Fermi Surface", "Hybrid Functional", "GW Calculation", "BSE calculation", which are reprensted by task="RLX", "SCF", "BND", "DOS", "CHG", "FSF", "HSE", "GWA", and "BSE", respectively. This class module can be used as an automative workflow to implement creating inputs from sratch and restarting from existing calculations to submitting calculations to popular job manager systems such as PBS and SLURM.

- `to be continued`: more scripts are yet to come ...


### Install
The simpliest way to use these scripts is to put them into a directory and then export the path of this directory as `PYTHONPATH`.
```
export PYTHONPATH=/the/path/of/the/directory:${PYTHONPATH}
```

Prior to using these python scripts, the following prerequisites need to be installed:
- numpy: https://numpy.org/
- ase: https://wiki.fysik.dtu.dk/ase/

and the following can also be installed but are optional for pycell.py and pycheck.py.
- scipy: https://www.scipy.org/
- matplotlib: https://matplotlib.org/


### Usage
```
python pycheck.py incar [poscar] [potcar] [outcar]
```

```
python pycell.py -i POSCAR -x 2 -y 2 -z 1 -v 20
Here `-x 2 -y 2 -z 1` builds a supercell of 2x2x1 and `-v 20` sets the vacuum size to 20 A.
````

```
workfunc.py is a Python Class which implements the following instance methods: 
wf = WorkFunction()
wf.GetFermilevel(outcar="OUTCAR")
wf.GetPotential(locpot="LOCPOT")
wf.PlanarAverage(axis="z")
wf.PlotPotential(show=True)
```

```
pyvasp.py is a Python class which implements an antomative workflow from creating inputs 
    from sratch and restarting from existing calculations to submitting calculations to 
    popular job manager systems such as PBS and SLURM.

from pyvasp import PythVasp

INP = "INPUT" # containing POSCAR and POTCAR

# Perform Relaxation
RLX = "RELAX" # working directory for RELAXATION
Job = PythVasp(task="RLX", srcdir=INP, objdir=RLX)
Job.MakeInput(encut=600, ispin=2, nsw=500, isif=3, ibrion=2, ediffg="-1E-2", algo="Normal", nelm=100, nelmin=8, ediff="1E-6", ismear=0, sigma=0.1, isym=2, lwave=True, lcharg=True, lorbit=11, nkpts=(8,8,2), gamma=True, kpar=4, npar=10, poscar="INPUT/POSCAR", potcar="INPUT/POTCAR", make_clean=True)
Job.Submission(queue="SLURM", ncores=40)

# Perform Static SCF 
SCF = "STATIC" # working directory for STATIC SCF
Job = PythVasp(task="SCF", srcdir=RLX, objdir=SCF)
Job.MakeInput(encut=600, ispin=2, icharg=1, lorbit=11, algo="Normal", nelm=300, nelmin=8, ediff="1E-6", ismear=0, sigma=0.1, lwave=False, lcharg=True, nkpts=(8,8,2), gamma=True, kpar=4, npar=10, make_clean=True)
Job.Submission(queue="SLURM", ncores=40)
