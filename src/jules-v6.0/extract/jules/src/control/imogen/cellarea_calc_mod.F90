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
MODULE cellarea_calc_mod

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CELLAREA_CALC_MOD'

CONTAINS

SUBROUTINE cellarea_calc(land_pts, lat, darea)
!-----------------------------------------------------------------------------
! Description:
!   This code corresponds tothe output of HADCM3, with a longitudinal
!   spacing of 3.75 degrees (96 points around a given latitudinal "band"
!   and a latitude spacing of 2.5 degrees with a top and bottom box
!   of 1.25 degrees.
! 
! Code Owner: Please refer to ModuleLeaders.txt
!             This file belongs in IMOGEN
!
! Written by: C. Huntingford (November 1999)
! 
! Code Description:
!   Language: Fortran 90.
!   
!-----------------------------------------------------------------------------

USE conversions_mod, ONLY: pi, pi_over_180
USE planet_constants_mod, ONLY: rad_earth => planet_radius

USE parkind1, ONLY: jprb, jpim
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                        &
  land_pts
               ! IN Number of land points

REAL, INTENT(IN) ::                                                           &
  lat(land_pts)
               ! IN Latitute (degrees)

REAL, INTENT(INOUT) ::                                                        &
  darea(land_pts)
               ! INOUT Gridbox area (m2)

INTEGER ::                                                                    &
  l

REAL ::                                                                       &
  lat_max,                                                                    &
               ! Latitute of "top" of GCM box (degrees)
  lat_min,                                                                    &
               ! Latitute of "bottom" of GCM box (degrees)
  lat_max_rad,                                                                &
               ! Latitute of "top" of GCM box (rad)
  lat_min_rad,                                                                &
               ! Latitute of "bottom" of GCM box (rad)
  area_earth
               ! Surface area of the earth


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CELLAREA_CALC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------

area_earth = 4.0 * pi * (rad_earth**2)

DO l = 1,land_pts
  IF (lat(l) >= 88.75) THEN
    lat_max = 90.0
    lat_min = 88.75
  ELSE IF (lat(l) <= -88.75) THEN
    lat_max = -88.75
    lat_min = -90.0
  ELSE
    lat_max = lat(l) + 1.25
    lat_min = lat(l) - 1.25
  END IF

  lat_max_rad = lat_max * pi_over_180
  lat_min_rad = lat_min * pi_over_180

  darea(l) =                                                                  &
    (2 * pi * (rad_earth**2) * (SIN(lat_max_rad) - SIN(lat_min_rad))) / 96.0

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE cellarea_calc
END MODULE cellarea_calc_mod
#endif

