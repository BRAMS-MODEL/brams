# Box-Matrix Makefile
# Author: Luiz Flávio Rodrigues (LFR)
# Date: 16Set2011
# Using F90 instead f77 

PROG =	aero

SRCS =	aero_driver.f90 aero_testsubs.f90 aero_tracer_com.f90 \
	TRAMP_config.F90 TRAMP_actv.f90 TRAMP_coag.f90 TRAMP_depv.f90 \
	TRAMP_diam.f90 TRAMP_dicrete.f90 TRAMP_init.f90 \
	TRAMP_issoropia.f90 TRAMP_isofwd.f90 TRAMP_isorev.f90 TRAMP_matrix.f90 \
	TRAMP_npf.f90 TRAMP_param_GISS.f90 TRAMP_quad.f90 TRAMP_setup.f90 \
	TRAMP_subs.f90 TRAMP_thermo_isorr.f90 TRAMP_isrpia.f90

OBJS =	aero_driver.o aero_testsubs.o aero_tracer_com.o \
	TRAMP_config.o TRAMP_actv.o TRAMP_coag.o TRAMP_depv.o TRAMP_diam.o \
	TRAMP_dicrete.o TRAMP_init.o TRAMP_issoropia.o \
	TRAMP_isofwd.o TRAMP_isorev.o TRAMP_matrix.o TRAMP_npf.o \
	TRAMP_param_GISS.o TRAMP_quad.o TRAMP_setup.o TRAMP_subs.o \
	TRAMP_thermo_isorr.o TRAMP_isrpia.o

LIBS =	

#Please, select the compilers:
CC = icc
CFLAGS = -O0
FC = ifort
FFLAGS = -O0
F90 = ifort

# LFR: Select Mechanism from 1 to 8 TRACERS_AMP_M[1-8]
# LFR: Please, see TRAMP_config.F90 to know about each one
F90FLAGS = -O0 -DTRACERS_AMP_M1 
# LFR: You may include others options using -D, like:
#-DTRACERS_ON
#-DTRACERS_WATER
#-DTRACERS_DRYDEP
#-DTRDIAG_WETDEPO
#-DTRACERS_AMP

LDFLAGS = -s

FCOMP= $(F90) $(F90FLAGS) -c
all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .F90 .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

#Dependencies:

aero_driver.o: TRAMP_actv.o TRAMP_coag.o TRAMP_config.o TRAMP_diam.o \
	TRAMP_dicrete.o TRAMP_init.o TRAMP_npf.o \
	TRAMP_param_GISS.o TRAMP_setup.o TRAMP_subs.o
        
aero_testsubs.o: TRAMP_actv.o TRAMP_config.o TRAMP_depv.o TRAMP_diam.o \
	TRAMP_param_GISS.o TRAMP_setup.o
        
aero_tracer_com.o: TRAMP_config.o TRAMP_param_GISS.o 

TRAMP_actv.o: TRAMP_config.o TRAMP_param_GISS.o

TRAMP_coag.o: TRAMP_config.o TRAMP_param_GISS.o TRAMP_setup.o

TRAMP_config.o : ./TRAMP_config.F90 TRAMP_param_GISS.o 
	$(FCOMP) $(<F:.F90=.F90)

TRAMP_depv.o: TRAMP_config.o TRAMP_param_GISS.o

TRAMP_diam.o: TRAMP_config.o TRAMP_param_GISS.o TRAMP_setup.o aero_tracer_com.o

TRAMP_dicrete.o: TRAMP_coag.o TRAMP_param_GISS.o

TRAMP_init.o: TRAMP_config.o TRAMP_dicrete.o TRAMP_param_GISS.o TRAMP_setup.o

TRAMP_matrix.o: TRAMP_actv.o TRAMP_coag.o TRAMP_config.o TRAMP_depv.o \
	TRAMP_diam.o TRAMP_npf.o TRAMP_setup.o TRAMP_subs.o
        
TRAMP_npf.o: TRAMP_param_GISS.o TRAMP_setup.o

TRAMP_setup.o: TRAMP_config.o TRAMP_param_GISS.o

TRAMP_subs.o: TRAMP_config.o TRAMP_param_GISS.o TRAMP_setup.o

TRAMP_thermo_isorr.o: TRAMP_param_GISS.o

TRAMP_issoropia.o: TRAMP_isrpia.o

TRAMP_isofwd.o : TRAMP_isrpia.o

TRAMP_isorev.o: TRAMP_isrpia.o
