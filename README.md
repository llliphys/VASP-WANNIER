### Intro
Here are some python scripts used to check the VASP input/output files: INCAR, POSCAR, POTCAR, KPOINTS, OSZICAR, OUTCAR, and so on. To check out how to use them, just simply command `python script_name.py` or `python script_name.py -h` on your terminal to have a first glance.

- `pycheck.py`: simply merge the above two scripts into a single script (more compact).

- `pycell.py`: manipulate CIF or POSCAR file, e.g., making supercell, varying vacuum size, etc.

- `workfunc.py`: calculate the work function and plot the planar averaged potental of a slab system.

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
python pycell.py -i POSCAR -x 2 -y 2 -z 1 -v 20
```
Here `-x 2 -y 2 -z 1` builds a supercell of
the size 2x2x1 and `-v 20` sets the vacuum size to 20 A.

```
workfunc.py is a Python Class which implements the following instance methods: 
wf = WorkFunction()
wf.GetFermilevel(outcar="OUTCAR")
wf.GetPotential(locpot="LOCPOT")
wf.PlanarAverage(axis="z")
wf.PlotPotential(show=True)
```
