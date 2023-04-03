F90=gfortran
F90FLAGS=
INCLUDES=-I../utils/include/
F90LOAD= $(F90FLAGS)

EXEC1=TimeStampReport.exe
OBJ1=TimeStampReport.o
EXEC2=TimeStampSummary.exe
OBJ2=TimeStampSummary.o

all:	$(EXEC1) $(EXEC2)

$(EXEC1):	$(OBJ1)
	$(F90) -o $(EXEC1) $(INCLUDES) $(F90LOAD) $(OBJ1)

$(EXEC2):	$(OBJ2)
	$(F90) -o $(EXEC2) $(INCLUDES) $(F90LOAD) $(OBJ2)

TimeStampReport.o:	TimeStampReport.f90 ../utils/include/i8.h

TimeStampSummary.o:	TimeStampSummary.f90 ../utils/include/i8.h

clean:
	rm -f *.o *.mod $(EXEC1) $(EXEC2)


.SUFFIXES	:	.o .f90

.f90.o:
	$(F90) -c $(F90FLAGS) $(INCLUDES) $*.f90



