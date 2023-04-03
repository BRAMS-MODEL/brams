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

SUBROUTINE diffcarb_land_co2(land_pts, d_land_atmos_co2, lat, dctot_co2, conv)

USE cellarea_calc_mod, ONLY: cellarea_calc

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Calculates the global change in atmospheric CO2 due to
!   grid box changes in total carbon content.
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

INTEGER, INTENT(IN) ::                                                        &
  land_pts
               ! IN Number of land points

REAL, INTENT(IN) ::                                                           &
  lat(land_pts),                                                              &
               ! IN Latitute (degrees)
  dctot_co2(land_pts),                                                        &
               ! IN Change in total surface gridbox CO2
               !    content during a period "YEAR_CO2" (kg C/m2)
  conv
               ! IN Converts global emission of C (Gt/yr)
               !    into atmospheric CO2 (ppm/GtC)


REAL, INTENT(INOUT) ::                                                        &
  d_land_atmos_co2
               ! OUT Change in atmospheric CO2 concentration
               !     as a result of land-atmosphere feedbacks
               !     between calls to SCENARIO (ppm/"YEAR_CO2")


REAL ::                                                                       &
  land_gain,                                                                  &
               ! Total gain in carbon by land (kg C)
  darea(land_pts)
               ! WORK Gridbox area (m2)

INTEGER ::                                                                    &
  l            ! Looping parameter



land_gain = 0.0

CALL cellarea_calc(land_pts, lat, darea)

DO l = 1,land_pts
  land_gain = land_gain + darea(l) * dctot_co2(l)
END DO

d_land_atmos_co2 = -(land_gain / 1.0e12) * conv

RETURN

END SUBROUTINE diffcarb_land_co2
#endif
