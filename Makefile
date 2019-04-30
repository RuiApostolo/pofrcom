######################################################
#
# Mixed MPI/OpenMP free- and fixed-format Fortran code
#
# Instructions:
#	- Choose a compiler and compile options
#	- "make clean"
#	- "make"
#
# If adding source files, make sure they are included
# in the "F77SRC" or "F08SRC" variables below. Add
# the appropriate dependencies in the bottom section
# (with care). There is no need to edit the targets
# section.
#
# Author: K.D'Mellow, Dec 2008.
#
######################################################

######################################################
#
# Build directory where compilation output files will
# be stored
#
######################################################
BUILD_DIR=./


# Executable

 EXE=    pr_com


######################################################
#
# Compiler and Flags Section
# (use only one compiler)
#
######################################################
#FC=     pgf90
#FC=     pgf90 -fast -O4
#FC=     pgf90 -fast -O4 -Minline=levels:1 -Minfo=inline
#FC=     pgf90 -fast -Mipa=fast,inline -Minfo=inline -tp core2-64
#FC=     pgf90 -fast -Mipa=fast,inline -Minfo=inline -tp x64
#FC=     pgf90 -fast -Mipa=fast,inline -tp x64
#FC=     pgf90 -O2 -Munroll=c:1 -Mnoframe -Mlre  -Mvect=sse -Mscalarsse -Mcache_align -Mflushz
#FC=     pgf90 -O2 -Munroll=c:1 -Mnoframe -Mlre -Mscalarsse -Mcache_align -Mflushz -Mbounds
#FC=     pgf90
#FC=     mpif90 -O4
#FC=     xlf90 -O4 -qfree
#FC=     ifort -O2
#FC=     f95 -O2
#FC=     gfortran -fbounds-check
#FC=     gfortran -O4  #-g -ffpe-trap=invalid,zero,overflow
#FC=	 gfortran -Wall -Wextra -Wconversion -pedantic -g -fbacktrace -ffpe-trap=zero,overflow,underflow
FC=	 gfortran -Wall -Wextra -Wconversion -pedantic -ffpe-trap=zero,overflow,underflow -O3

# C Preprocessing Directives 
CPPFLAGS= $(DEBUG)

#OpenMP flags
#OMPFLAGS= -fopenmpcd 
#OMPFLAGS= -mp -Minfo=mp
#OMPFLAGS= -mp
#OMPFLAGS = -fopenmp # For gfortran
#OMPFLAGS -openmp # For ifort

#Debugging
#DEBUGGING= -g -fbacktrace -fcheckall

#Any kind of profiling
#PROFILING= -pg

#
# library output for packaging up external routines.
#
#LIBRARYROOT=mylib
#LIBRARY=lib$(LIBRARYROOT).a
#LIB= -L. -l$(LIBRARYROOT)

FFLAGS=     $(DEBUGGING) $(PROFILING) $(OMPFLAGS) $(CPPFLAGS)
# Or if just want to preprocess
#FFLAGS = -F -g $(PROFILING) $(OMPFLAGS) $(CPPFLAGS)

LFLAGS= $(FFLAGS) $(LIB)



######################################################
#
# Source Files Section - add files as appropriate
#
######################################################
#Makefile is integrated into make dependencies, so specify its name
MF=Makefile

#Fortran 77 source files
F77SRC=

#Fortran 90 source files
F08SRC= \
        RDF_com.f08


#######################################################
#
# Makefile Targets: No need to edit this section
#
#######################################################
.SUFFIXES:
.SUFFIXES: .f08 .f .mod .o .inc

OBJ=    $(F78SRC:.f=.o) $(F08SRC:.f08=.o)
SRC=    $(F77SRC) $(F08SRC) $(INC)

.f08.o:
	$(FC) $(FFLAGS) -c $<

.f08.mod:
	$(FC) $(FFLAGS) -c $<

.f.o:
	$(FC) $(FFLAGS) -c $<

all:    $(EXE)

$(EXE): $(OBJ)
	$(FC) $(LFLAGS) -o $@ $(OBJ)

$(OBJ): $(MF) $(INC)

tar:
	tar cvf $(EXE).tar $(MF) $(SRC)

#libraries: my_libexternals.o
#	ar -r $(LIBRARY) $<


clean:
	rm -f $(OBJ) $(EXE) *.oo *.ipa *.ipo *.mod *.out core



######################################################
#
# Dependencies section. Modify this with care as
# appropriate to satisfy use-based dependencies within
# modules.
#
######################################################
# RDF_com.o: modules.o 

