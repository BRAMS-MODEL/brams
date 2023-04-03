!f90
!############################# Change Log ##################################
! 4.3.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!  Mission Research Corporation / *ASTeR Division
!###########################################################################

!PARAMETER's retirados do programa principal: ramspost_A.f90
INTEGER, PARAMETER :: maxpatch = 5
INTEGER, PARAMETER :: maxgx    = 800  ! maxgx=int(1.2*float(nxpmax)),maxgy=int(1.2*float(nypmax))
INTEGER, PARAMETER :: maxgy    = 800

!------------------------------------------------------------------------------------------------
!  Set maximum values of PARAMETERs:

INTEGER, PARAMETER :: MAXGRDS  = 4    ! MAXGRDS - Maximum number of grids
                                      !
INTEGER, PARAMETER :: NXPMAX   = 800  ! NXPMAX  - Maximum number of points
                                      !           in x-direction
INTEGER, PARAMETER :: NYPMAX   = 800  ! NYPMAX  - Maximum number of points
                                      !           in y-direction
INTEGER, PARAMETER :: NZPMAX   = 120  ! NZPMAX  - Maximum number of points
                                      !           in z-direction
INTEGER, PARAMETER :: NZGMAX   = 12   ! NZGMAX  - Maximum number of soil levels
                                      !
INTEGER, PARAMETER :: MAXSCLR  = 50   ! MAXSCLR - Maximum number of additional
                                      !           scalars
INTEGER, PARAMETER :: MAXHP    = 1000 ! MAXHP   - Maximum number of u, v, OR t
                                      !           points in a single vertical
                                      !           level interpolated from 
                                      !           opposite hemispheric grid

!------------------------------------------------------------------------------------------------
!  Set MAXDIM to the largest of NXPMAX,NYPMAX,NZPMAX+10,NZGMAX

INTEGER, PARAMETER :: MAXDIM  = 800

!------------------------------------------------------------------------------------------------
!  maxmach is the max number of processors that can be used in a parallel run

INTEGER, PARAMETER :: maxmach   = 1

INTEGER, PARAMETER :: nxyzpm    = nzpmax*nxpmax*nypmax,                    &
                      maxdimp   = maxdim+1, nstyp   = 12,    nvtyp=30,     &
                      nkeep     = 90,       nke     = nkeep, nintgm=12,    &
                      maxsndg   = 200,      maxvarf = 200,   maxsstf = 200

INTEGER, PARAMETER :: maxsched  = 200
INTEGER, PARAMETER :: maxschent = 5
INTEGER, PARAMETER :: maxfiles  = 2000
INTEGER, PARAMETER :: nplmax    = 50

!------------------------------------------------------------------------------------------------
!PARAMETER's retirados das funções: calccine, calcnce, calccape, calcnpe,
!tempvirtual, razaodemistura, potencial, potencialeq, tparcela: ramspost_D.f90
REAL,    PARAMETER :: epsi     = 0.62198

!------------------------------------------------------------------------------------------------
!PARAMETER's retirados da função: integrando: ramspost_C.f90
REAL,    PARAMETER :: ra       = 287.04

!------------------------------------------------------------------------------------------------
!PARAMETER's retirados da subroutina R_POLARST, ll_xy: ramspost_C.f90, polarst.f90
REAL,    PARAMETER :: erad     = 6367000.

!------------------------------------------------------------------------------------------------
!PARAMETER's retirados da subroutina SEAPRS_0: ramspost_B.f90
REAL,    PARAMETER :: R        = 287.04
REAL,    PARAMETER :: G        = 9.8
REAL,    PARAMETER :: GAMMA    = 287.04
REAL,    PARAMETER :: TC       = 273.16+17.5 ! T CRITICAL IN PSFC/PSLV
REAL,    PARAMETER :: PCONST   = 100.

!------------------------------------------------------------------------------------------------
!PARAMETER's retirados da subrotina ge_to_xy: ramspost_A.f90
REAL,    PARAMETER :: rt       = 6367000.00
!------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------------
!PARAMETER's retirados da subrotina ll_xy,xy_ll: polarst.f90
REAL,    PARAMETER :: erad2    = 1.2734e7
REAL,    PARAMETER :: pi180    = 3.14159265/180.
!------------------------------------------------------------------------------------------------




