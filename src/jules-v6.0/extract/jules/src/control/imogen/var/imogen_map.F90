#if !defined(UM_JULES)
!******************************COPYRIGHT**************************************
! (c) Centre for Ecology and Hydrology. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237] 
!******************************COPYRIGHT**************************************

MODULE imogen_map

USE imogen_constants, ONLY: n_imogen_land

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Module for Various mapping routines for JULES-IMOGEN simulations
!
! Code Owner: Please refer to ModuleLeaders.txt
!             This file belongs in IMOGEN
!
! Code Description:
!   Language: Fortran 90.
!   
!-----------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Module variables
!------------------------------------------------------------------------------
INTEGER, DIMENSION(n_imogen_land) :: sgindinv

CONTAINS

SUBROUTINE get_imogen_map(imogenOrderFile)
!------------------------------------------------------------------------------
! Module imports
!------------------------------------------------------------------------------
USE io_constants, ONLY: imogen_unit

USE model_grid_mod, ONLY: latitude,longitude

USE ancil_info, ONLY: land_pts

USE jules_fields_mod, ONLY: ainfo

USE theta_field_sizes, ONLY: t_i_length

USE imogen_constants, ONLY: n_imogen_land

IMPLICIT NONE

CHARACTER(LEN=*) ::                                                           &
  imogenOrderFile ! Filename to read IMOGEN points order from

!------------------------------------------------------------------------------
! Local variable declarations
!------------------------------------------------------------------------------
INTEGER :: indlat,indlon,i,j,l,map1(96,56)
INTEGER :: sgjind(land_pts),sgind(land_pts)

OPEN(imogen_unit, FILE=imogenOrderFile,                                       &
                       STATUS='old', POSITION='rewind', ACTION='read')
READ(imogen_unit,*) map1
CLOSE(imogen_unit)

DO l = 1,land_pts
  j = (ainfo%land_index(l) - 1) / t_i_length + 1
  i = ainfo%land_index(l) - (j-1) * t_i_length

  indlat = INT((latitude(i,j) + 55.1) / 2.5) + 1
  indlon = INT((longitude(i,j) + 180.1) / 3.75) + 1
  sgind(l) = map1(indlon, 57 - indlat)
  sgjind(l) = (indlat-1) * 96 + indlon
END DO

sgindinv(:) = 0
DO i = 1,n_imogen_land
  DO j = 1,land_pts
    IF (i == sgind(j)) sgindinv(i) = j
  END DO
END DO

END SUBROUTINE get_imogen_map

END MODULE imogen_map
#endif
