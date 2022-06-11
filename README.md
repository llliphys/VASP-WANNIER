### Intro
Here are some python scripts to pre- and post-process DFT-VASP and VASP2WANNIER calculations. To check them out, just simply command `python script_name.py` or `python script_name.py -h` on your terminal to have a first insight.

- `pycheck.py`: Inspect VASP input and output files, e.g. checking if input setup is correct.

- `pycell.py`: Manipulate CIF or POSCAR file, e.g., making supercell, varying vacuum size, etc.

- `workfunc.py`: Calculate the work function and plot the planar averaged potental of a slab system.

- `pythvasp.py`: A Python class module implementing a series of DFT-VASP calculations including e.g. "Relaxation", "Static SCF", "Band Structure", "Density of States", "Charge Density", "Fermi Surface", "Hybrid Functional", "GW Calculation", "BSE calculation", which are launched by setting task="RLX", "SCF", "BND", "DOS", "CHG", "FSF", "HSE", "GWA", and "BSE", respectively. This class module allows to create inputs from sratch, to restart from existing calculations, and to submit calculations in PBS and SLURM scripts.

- `pythwann.py`: A Python class module implementing a series of VASP2WANNIER calculations including e.g. "Wannierization", "Band Interplolation", "Berry Curvature", "Anomalous Hall Conductivity", which are launched by setting task="WAN", "BND", "BC", and "AHC". This class module allows to create inputs from sratch, to restart from existing calculations, and to submit calculations in PBS and SLURM scripts.

- `to be continued`: More scripts are yet to come ...


### Install
The simpliest way to use these scripts is to put them into a directory and then export the path of this directory as `PYTHONPATH`.
```
export PYTHONPATH=/the/path/of/the/directory:${PYTHONPATH}
```

To run these scripts smoothly, the following prerequisites need to be installed:
- numpy: https://numpy.org/
- ase: https://wiki.fysik.dtu.dk/ase/
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
Here below presents a basic usage of pythvasp.py:

from pythvasp import PythVasp

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
