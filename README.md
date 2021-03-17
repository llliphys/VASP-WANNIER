### VaspScripts
__________________
Here are some python scripts used to check the VASP input/output files: INCAR, POSCAR, POTCAR, KPOINTS, OSZICAR, OUTCAR, EIGENVAL, DOSCAR, and so on. To check out how to use them, just simply command `python script.py` to have a first glance.

- `check_input.py`: check INCAR, POSCAR, POTCAR, and KPOINTS to see each setup is correct.

- `check_output.py`: check OUTCAR to see if a relaxation or scf calculation is converged.

- `to be continued`: more scripts are yet to come ...

### Installation
__________________
The simpliest way is to put these scripts into a directory and then export the path of this directoy as `PYTHONPATH`.
```
export PYTHONPATH=/the/path/of/the/directory:${PYTHONPATH}
```
