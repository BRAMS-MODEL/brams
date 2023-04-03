REV=1.0
EXE=../bin/geraSM_$(REV)

### IMPORTANT ###
#*** dont forget to use 'make clean' after changing any of these options ***

GRIBLIB=-L$(GRIB2_PATH)/lib -lwgrib2_api -lwgrib2
GRIBINC=-I$(GRIB2_PATH)/include -I$(GRIB2_PATH)/lib -I../src
NETCDFCLIB=-L$(NETCDFC_PATH)/lib -lnetcdf
NETCDFFLIB=-L$(NETCDFF_PATH)/lib -lnetcdff
NETCDFINC=-I$(NETCDFF_PATH)/include -I$(NETCDFC_PATH)/include
BRAMSINCLUDES=-I$(BRAMS_INCLUDE)

FC           =gfortran
FC_OPTS      = -DUSEPARALELL=$(USEPARALELL) -DUSEGRIB=$(USEGRIB) -ffree-form -ffree-line-length-none
FLOADER      =gfortran

