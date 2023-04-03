#

include paths.mk
include include.mk.$(OPT)

#
INCLUDES=$(GRIBINC) $(NETCDFINC) $(SRCINCLUDES) $(BRAMSINCLUDES)
LIBS=$(GRIBLIB) $(NETCDFCLIB) $(NETCDFFLIB)

#

F_COMMAND = $(FC) -c $(FC_OPTS) $(INCLUDES)

#

include objects.mk


all: $(EXE)
	@echo "---------------"
	@echo $(EXE) done!
	@echo "---------------"

$(EXE): $(OBJ) $(MAIN_OBJ)
	$(FLOADER) $(FLOADER_OPTS) -o $(EXE) $(MAIN_OBJ) $(OBJ)	$(LIBS)

clean:
	rm -fr *.o *.mod $(EXE)

execP:
	/opt/scali/bin/mpirun -machinefile host -np 12 $(EXE)

execS:
	./$(EXE)


include depend.mk
