# pofrcom
Calculates p(r) for the centre of mass of a functional group from LAMMPS output

To compile, use the make command with the included makefile or:
$ gfortran -Wall -Wextra -Wconversion -pedantic -ffpe-trap=zero,overflow,underflow -O3 RDF_com.f08 -o pofrcom

Requires:
Masses.in is a two collumn, space separated file with atom type on column #1 and mass on column #2.
Params.in has all the following configurable settings (first 12 columns are option, anything after that is not read):
lammpstrj_filename
Number of columns in lammpstrj_filename
Number of steps in lammpstrj_filename
Ignore first X steps in lammpstrj_filename
Number of atoms in target molecule
Number of FGs in target molecule
Number of FG types in target molecule
Size of FG in atoms
Atomtype of first atom in FG
Boxlength (needs ot be changed to be read from input)
SolventNum (not implemented)
Size of histogram bin
Lowest bin
Number of bins
Logical switch to calculate 2D (XY) p(r) (T/F)
