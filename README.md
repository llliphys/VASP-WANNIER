### Intro
---
Here are some python scripts used to check the VASP input/output files: INCAR, POSCAR, POTCAR, KPOINTS, OSZICAR, OUTCAR, EIGENVAL, DOSCAR, and so on. To check out how to use them, just simply command `python script.py` or `python script.py -h`to have a first glance.

- `check_input.py`: check INCAR, POSCAR, POTCAR, and KPOINTS to see each setup is correct.

- `check_output.py`: check OUTCAR to see if a relaxation or scf calculation is converged.

- `pycheck.py`: simply merge the above two scripts into a single script (more compact).

- `pycell.py`: manipulate CIF or POSCAR file, e.g. making supercell, varying vacuum lenght.

- `to be continued`: more scripts are yet to come ...

### Install
---
The simpliest way is to put these scripts into a directory and then export the path of this directoy as `PYTHONPATH`.
```
export PYTHONPATH=/the/path/of/the/directory:${PYTHONPATH}
```

### Usage
---
```
python check_input.py incar [ediff] [isif]
python check_input.py poscar [lattice] [position]
python check_input.py potcar [titel] [encut]
python check_output.py outcar [energy] [magmom]
python pycheck.py incar/poscar/potcar/outcar
python pycell.py -i POSCAR -x 2 -y 2 -z 1 -v 20
```
