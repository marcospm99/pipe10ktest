INSTDIR		= ./install/
PROGDIR		= ./program/
UTILDIR		= ./utils/
UTIL		= prim2matlab

TRANSFORM	= fftw3
MODSOBJ		= type.o wksp.o iohdf5.o meshs.o mpi.o parameters.o \
		  timestep.o transform.o variables.o velocity.o sta.o

#COMPILER	= g95 -C  
#COMPFLAGS	= -cpp -c -O3
#COMPILER	= pgf90 #-C
#COMPFLAGS	= -Mpreprocess -c -fast #-mcmodel=medium
#COMPILER	= pathf90 #-pg -C
#COMPFLAGS	= -cpp -c -O3 -OPT:Ofast -march=opteron -fno-second-underscore
#COMPILER	= gfortran
#COMPFLAGS	= -ffree-line-length-none -x f95-cpp-input -c -O3 \
		  -I/usr/include \
                  #-C #-pg


COMPILER	= h5pfc 
##COMPILER	= ifort
##COMPFLAGS	= -g -cpp -c -O3 -heap-arrays 1024 -fpe0 -check all -warn unused -debug extended
### -mcmodel=medium seems to be neeessary because all the memory is static! 
### DBFLAGS = -g -warn unused -debug extended -check all -fPIC -fpe0 

FLAGS           = -mcmodel=medium -cpp -O3 -qmkl
##LIBRARYF  =  -I/usr/include/ -qmkl -I$(MKLROOT)/include/fftw
##HDF5_LINK=-I/dss/dsshome1/lrz/sys/spack/release/22.2.1/opt/skylake_avx512/hdf5/1.10.7-intel-v7t5g3t/include
##HDF5_F90_LIB = -L/dss/dsshome1/lrz/sys/spack/release/22.2.1/opt/skylake_avx512/hdf5/1.10.7-intel-v7t5g3t/lib/ -lhdf5_hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5

COMPFLAGS = -c $(FLAGS) $(DBFLAGS) $(LIBRARYF) $(LIBRARYF) $(HDF5_LINK)
LIBS      =  cheby.o
LFLAGS    = $(FLAGS) $(DBFLAGS) $(LIBRARYF) $(LIBRARYF)


#------------------------------------------------------------------------
all : 	$(MODSOBJ) $(PROGDIR)main.f90
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)main.f90
#	$(COMPILER) -list=cosita.lst -o ./pipe10k main.o $(MODSOBJ) $(LIBS) -g -fpe0 -check all -warn unused -debug extended
	$(COMPILER) -o ./pipe10k main.o $(MODSOBJ) $(LIBS) $(LFLAGS) $(HDF5_F90_LIB)

install : pipe10k
	if test ! -d $(INSTDIR); then mkdir -p $(INSTDIR); fi
	mv ./$@ $(INSTDIR)
	date > $(INSTDIR)/main.info
	echo $(HOSTNAME) >> $(INSTDIR)/main.info
	pwd >> $(INSTDIR)/main.info
	echo $(COMPILER) $(COMPFLAGS) >> $(INSTDIR)/main.info
	grep "define _N" parallel.h >> $(INSTDIR)/main.info
	cut -d! -f1 $(PROGDIR)parameters.f90 | grep = | \
	   cut -d: -f3  >> $(INSTDIR)/main.info

util : 	$(MODSOBJ) $(UTILDIR)/$(UTIL).f90
	$(COMPILER) $(COMPFLAGS) $(UTILDIR)/$(UTIL).f90
	$(COMPILER) -o ./$(UTIL).out $(UTIL).o $(MODSOBJ) $(LIBS)

#------------------------------------------------------------------------
clean :
	rm -f *.o *.mod *.d *.il core *.out

#------------------------------------------------------------------------
type.o : $(PROGDIR)type.f90 parameters.o parallel.h
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)type.f90

wksp.o : $(PROGDIR)wksp.f90 type.o parallel.h
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)wksp.f90

iohdf5.o : $(PROGDIR)iohdf5.F90 velocity.o parallel.h
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)iohdf5.F90

meshs.o : $(PROGDIR)meshs.f90 parameters.o mpi.o
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)cheby.f
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)meshs.f90

mpi.o : $(PROGDIR)mpi.f90 parallel.h
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)mpi.f90

parameters.o : $(PROGDIR)parameters.f90 parallel.h
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)parameters.f90

timestep.o : $(PROGDIR)timestep.f90 variables.o
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)timestep.f90

transform.o : $(PROGDIR)transform.$(TRANSFORM).f90 variables.o
	$(COMPILER) $(COMPFLAGS) -o transform.o \
	$(PROGDIR)transform.$(TRANSFORM).f90

variables.o : $(PROGDIR)variables.f90 wksp.o meshs.o
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)variables.f90

velocity.o : $(PROGDIR)velocity.f90 timestep.o transform.o
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)velocity.f90

sta.o : $(PROGDIR)sta.f90 timestep.o transform.o
			$(COMPILER) $(COMPFLAGS) $(PROGDIR)sta.f90