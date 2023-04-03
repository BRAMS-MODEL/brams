#if defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

MODULE um_latlon_mod

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   This module provides equivalent latitude/longitude definitions to 
!   standalone used for running river routing with MetUM
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Module variables
!-----------------------------------------------------------------------------

REAL, ALLOCATABLE :: longitude(:,:)  ! True longitude values on local task
REAL, ALLOCATABLE :: latitude(:,:)   ! True latitude values on local task

! Index of land points on global (full domain) grid, defined on land points
INTEGER, ALLOCATABLE :: global_land_index(:)
! Number of land points on the global (full domain) grid
INTEGER :: global_land_pts

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UM_LATLON_MOD'

CONTAINS

!-----------------------------------------------------------------------------
! Description:
!   Set global (full domain) latitude, longitude and mask variables,
!   using UM module variables. These are required for regridding when
!   using RFM and TRIP river routing schemes
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
SUBROUTINE um_latlon()

USE theta_field_sizes,    ONLY: t_i_length, t_j_length
USE trignometric_mod,     ONLY: true_longitude,true_latitude
USE conversions_mod,      ONLY: recip_pi_over_180
USE UM_ParVars,           ONLY: glsize 
USE Field_Types,          ONLY: fld_type_p

USE atm_land_sea_mask, ONLY: atmos_landmask, atmos_number_of_landpts

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! Working variables
INTEGER :: i,j,l   ! Loop index (l: land point index)

! Mask of land points on global (full domain) grid
LOGICAL :: global_land_mask(glsize(1,fld_type_p),glsize(2,fld_type_p))

!Error reporting
CHARACTER(LEN=*), PARAMETER  :: RoutineName = 'UM_LATLON'

!Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Define JULES module variable from UM variable
global_land_pts = atmos_number_of_landpts

! Allocate the index of land points on global grid
ALLOCATE(global_land_index(global_land_pts))

! Allocate module variables on local task
ALLOCATE(latitude(t_i_length, t_j_length))
ALLOCATE(longitude(t_i_length, t_j_length))

! Set local task latitude and longitude based on UM module variables
DO j = 1,t_j_length
  DO i = 1,t_i_length
    latitude(i,j)  =  recip_pi_over_180 * true_latitude(i,j)
    longitude(i,j) =  recip_pi_over_180 * true_longitude(i,j)
  END DO
END DO

! Translate UM variable atmos_landmask from 1D vector to 2D field
DO l = 1, glsize(1,fld_type_p) * glsize(2,fld_type_p)
  global_land_mask( MOD( l - 1, glsize(1,fld_type_p) ) + 1,                   &
                    (l-1) / glsize(1,fld_type_p) + 1 ) = atmos_landmask(l)
END DO

! Set land points indexing on full domain grid
l = 1
DO j = 1,glsize(2,fld_type_p)
  DO i = 1,glsize(1,fld_type_p)
    IF ( global_land_mask(i,j) ) THEN
      global_land_index(l) = (j-1) * glsize(1,fld_type_p) + i
      l = l + 1
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE um_latlon

END MODULE um_latlon_mod
#endif
