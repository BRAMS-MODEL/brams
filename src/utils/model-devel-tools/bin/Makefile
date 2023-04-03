DSRC1=../../../brams/cuparm
DSRC2=../src
F_COMP=pgf90
F_OPTS= -O0 -Mbounds -traceback #-r8  -Mbounds -traceback 
#F_COMP=gfortran 
#F_OPTS=-O0 -fcheck=bounds -g -fcheck=all -finit-real=nan \
# -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,denormal #-Wall
#  -finit-real=zero -Wall

LOADER=$(F_COMP)
LOADER_OPTS=-v #-Bstatic -Mbounds -traceback #-Mlarge_arrays
ARCHIVE=ar r

#LIBS =	-L/scratchin/grupos/catt-brams/shared/libs/intel/netcdf-4.1.1/lib -lnetcdf
#LIBS =	-L/Users/saulofreitas/localbin/libs/pgi/netcdf_new/lib -lnetcdf

# Compiler commands.

F_COMMAND = $(F_COMP) -c $(F_OPTS) $(INCLUDES)

# Define archive and executable names.

BASE=gf
EXE=$(BASE).x
ARC=$(BASE).a

# Define source.
F_SRC=  $(DSRC2)/module_gate.f90 \
        $(DSRC1)/Phys_const.f90 \
        $(DSRC1)/module_cu_gf_v5.1.f90 \
        $(DSRC2)/module-from-1d-oldversion.f90 \
        $(DSRC1)/module_cu_gd_fim.f90 \
	$(DSRC2)/1d-gf-test.f90 


# Define targets.

all: $(EXE)

$(EXE): $(ARC) 
	$(LOADER) -o $(EXE) $(LOADER_OPTS) $(ARC) $(LIBS)
###	ln -fs $(EXE) $(BASE)
	rm -f r*.o

$(ARC): $(F_SRC)
	$(F_COMMAND) $(?)
	$(ARCHIVE) $(ARC) *.o
#	$(F_COMMAND) $(DSRC)/main.f90
#	rm -f *.o


install:
	ln -fs $(pwd)/$(EXE) ../$(BASE)
#	ln -fs `pwd`/$(EXE) ../test/$(BASE)

clean:
	rm -f $(ARC) $(EXE) $(BASE) *.o  ___* core* *.mod
