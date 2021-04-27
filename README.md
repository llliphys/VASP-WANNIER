### Intro
Here are some python scripts used to check the VASP input/output files: INCAR, POSCAR, POTCAR, KPOINTS, OSZICAR, OUTCAR, EIGENVAL, DOSCAR, and so on. To check out how to use them, just simply command `python script_name.py` or `python script_name.py -h` to have a first glance.

- `check_input.py`: check INCAR, POSCAR, POTCAR, and KPOINTS to see each setup is correct.

- `check_output.py`: check OUTCAR to see if a relaxation or scf calculation is converged.

- `pycheck.py`: simply merge the above two scripts into a single script (more compact).

- `pycell.py`: manipulate CIF or POSCAR file, e.g., making supercell, varying vacuum size, etc.

- `workfunc.py`: calculate the work function and plot the planar averaged potental of a slab system.

- `to be continued`: more scripts are yet to come ...

---

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

---

### Usage
```
python check_input.py incar [ediff] [isif]
python check_input.py poscar [lattice] [position]
python check_input.py potcar [titel] [encut]
python check_output.py outcar [energy] [magmom]
python pycheck.py incar [poscar] [potcar] [outcar]
python pycell.py -i POSCAR -x 2 -y 2 -z 1 -v 20
```
