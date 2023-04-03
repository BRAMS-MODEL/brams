REV=1.0
EXE=../bin/lq_$(REV)

### IMPORTANT ###
#*** dont forget to use 'make clean' after changing any of these options ***

BRAMSINCLUDES=-I$(BRAMS_INCLUDE)

FC           =gfortran
FC_OPTS      = -ffree-form -ffree-line-length-none
FLOADER      =gfortran

