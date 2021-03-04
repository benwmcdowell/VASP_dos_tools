# VASP_dos_tools
Tools for evaluating the site projected and total density of states (DOS) from a VASP calculation

Included files:
- plot_dos.py | plots the site projected DOS from the DOSCAR of a VASP calculation
                options:
                      -a, --atomnums: plot the site projected DOS for the atoms at specified indices (ie. 0,4,355)
                      -t, --types: plot the site projected DOS for atoms of specified types (ie. Ag, C)
                      -i, --integrated: plot the integrated DOS instead, integrated between specified energies relative to Ef (ie 0,3 integrates from the Fermi level to 3 eV above the Fermi level)
