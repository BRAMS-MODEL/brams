!******************************COPYRIGHT**************************************
! (c) Centre for Ecology and Hydrology. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237] 
!******************************COPYRIGHT**************************************

SUBROUTINE radf_non_co2(                                                      &
  year,q,nyr_non_co2,file_non_co2,file_non_co2_vals                           &
)

USE imogen_constants, ONLY: nyr_max
USE io_constants, ONLY: imogen_unit

USE logging_mod, ONLY: log_fatal

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Calculates the radiative forcing due to non-CO2 GHGs. This version
!   is consistent with the HadCM3 GHG run (AAXZE).
!
!   Written by Peter Cox (Sept 1998)
!   Adjusted by Chris Huntingford (Dec 1999) to include
!   other scenarios copied from file.
!
! Code Owner: Please refer to ModuleLeaders.txt
!             This file belongs in IMOGEN
!
! Code Description:
!   Language: Fortran 90.
!   
!-----------------------------------------------------------------------------

INTEGER ::                                                                    &
  year,                                                                       &
           ! IN Julian year.
  i,                                                                          &
           ! WORK Loop counter.
  nyr_non_co2
           ! IN Number of years for which NON_CO2
           !    forcing is prescribed.

CHARACTER(LEN=180) ::                                                         &
  file_non_co2_vals   ! IN File of non-co2 radiative forcings

LOGICAL ::                                                                    &
  file_non_co2        ! .T. if non_co2 forcings are read in from
                      ! a data file

REAL :: q             ! OUT non-CO2 radiative forcing (W/m2).

!-----------------------------------------------------------------------
! Local parameters
!-----------------------------------------------------------------------

INTEGER ::                                                                    &
  years(nyr_max)      ! Years for which radiative is
                      ! prescribed.

REAL ::                                                                       &
  q_non_co2(nyr_max),                                                         &
                      ! Specified radiative forcings (W/m2).
  growth_rate         ! Growth rate in the forcing after the
                      ! last prescribed value (%/pa).
PARAMETER (growth_rate = 0.0)

DATA years   / 1859,  1875,  1890,  1900,  1917,                              &
               1935,  1950,  1960,  1970,  1980,                              &
               1990,  2005,  2020,  2030,  2040,                              &
               2050,  2060,  2070,  2080,  2090,                              &
               2100,  679 * 2100 /

DATA q_non_co2 / 0.0344,0.0557,0.0754,0.0912,0.1176,                          &
                0.1483,0.1831,0.2387,0.3480,0.4987,                           &
                0.6627,0.8430,0.9225,0.9763,1.0575,                           &
                1.1486,1.2316,1.3025,1.3604,1.4102,                           &
                1.4602,679 * 1.4602 /




IF ( .NOT. file_non_co2 .AND. nyr_non_co2 /= 21) THEN
  CALL log_fatal("RADF_NON_CO2", 'Reset value of nyr_non_co2')
END IF

IF (nyr_non_co2 > nyr_max) THEN
  CALL log_fatal("RADF_NON_CO2", 'nyr_non_co2 too large')
END IF


!-----------------------------------------------------------------------
! File of non-co2 forcings read in if required
!-----------------------------------------------------------------------
IF (file_non_co2) THEN
  OPEN(imogen_unit, FILE=file_non_co2_vals,                                   &
                    STATUS='old', POSITION='rewind', ACTION='read')
  DO i = 1,nyr_non_co2
    READ(imogen_unit,*) years(i),q_non_co2(i)
  END DO
  CLOSE(imogen_unit)
END IF

!-----------------------------------------------------------------------
! Now calculate the non_co2 forcing
!-----------------------------------------------------------------------
IF (year < years(1)) THEN
  q = 0.0
ELSE IF (year > years(nyr_non_co2)) THEN
  q = q_non_co2(nyr_non_co2)                                                  &
    * ((1.0+0.01 * growth_rate)**(year - years(nyr_non_co2)))
ELSE
  DO i = 1,nyr_non_co2-1
    IF ((year >= years(i)) .AND. (year <= years(i+1))) THEN
      q = q_non_co2(i) + (year - years(i))                                    &
        * (q_non_co2(i+1) - q_non_co2(i)) / (years(i+1) - years(i))
    END IF
  END DO
END IF

RETURN

END SUBROUTINE radf_non_co2
