!f90
!############################# Change Log ##################################
! 4.3.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!  Mission Research Corporation / *ASTeR Division
!###########################################################################

! Estrutura data retirada do Programa Principal do arquivo ramspost_A.f90
CHARACTER(LEN = 3), DIMENSION(12), PARAMETER :: cmo = (/'jan','feb','mar',     & 
                                                        'apr','may','jun',     & 
                                                        'jul','aug','sep',     &
                                                        'oct','nov','dec'/)

! Estrutura data retirada da subrotina rams_anal_init do arquivo ramspost_B.f90
INTEGER,            DIMENSION(13), PARAMETER :: mondays =  (/31,28,31,30,      & 
                                                             31,30,31,31,      &
                                                             30,31,30,31,31/)

! Estrutura data retirada da subrotina rams_comp_slmstf do arquivo ramspost_B.f90

REAL,               DIMENSION(12), PARAMETER :: slmsts0  = (/0.395,0.410,0.435, &
                                                             0.485,0.451,0.420, &
                                                             0.477,0.476,0.426, &
                                                             0.492,0.482,0.863/)

! Estrutura data retirada da subrotina R_TRANSFM do arquivo ramspost_C.f90

REAL,                              PARAMETER :: TERDEV = 0.
