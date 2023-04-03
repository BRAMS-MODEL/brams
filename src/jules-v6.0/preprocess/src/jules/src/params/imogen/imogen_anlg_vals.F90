!******************************COPYRIGHT**************************************
! (c) Centre for Ecology and Hydrology. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237] 
!******************************COPYRIGHT**************************************

MODULE imogen_anlg_vals

USE missing_data_mod, ONLY: rmdi

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Default parameters and variables required for the imogen analogue model. 
!     Values can be set in the imogen.nml
!
! Code Owner: Please refer to ModuleLeaders.txt
!             This file belongs in IMOGEN
!
! Code Description:
!   Language: Fortran 90.
!   
!-----------------------------------------------------------------------------

REAL ::                                                                       &
  q2co2 = 3.74,                                                               &
              ! Radiative forcing due to doubling CO2 (W/m2
  f_ocean = 0.711,                                                            &
              ! Fractional coverage of the ocean
  kappa_o = 383.8,                                                            &
              ! Ocean eddy diffusivity (W/m/K)
  lambda_l = 0.52,                                                            &
              ! Inverse of climate sensitivity
              !   over land (W/m2/K)
  lambda_o = 1.75,                                                            &
              ! Inverse of climate sensitivity
              !   over ocean (W/m2/K)
  mu = 1.87,                                                                  &
              ! Ratio of land to ocean temperature
              !   anomalies
  t_ocean_init = 289.28,                                                      &
              ! Initial ocean temperature (K)
  diff_frac_const_imogen = rmdi
              ! Fraction of SW radiation that is diffuse for IMOGEN

INTEGER ::                                                                    &
  nyr_non_co2 = 21
              ! Number of years for which NON_CO2 forcing
              ! is prescribed.

CHARACTER(LEN=180) ::                                                         &
  dir_patt = '',                                                              &
              ! Directory containing the patterns
  dir_clim = '',                                                              &
              ! Directory containing initialising climatology.
  dir_anom = ''
              ! Directory containing prescribed anomalies

LOGICAL ::                                                                    &
  file_non_co2 = .FALSE.
              ! If true, then non-CO2 radiative forcings ar
              ! contained within a file.


NAMELIST  / imogen_anlg_vals_list/ q2co2,f_ocean,kappa_o,                     &
                                 lambda_l,lambda_o,mu,                        &
                                 t_ocean_init,dir_patt,                       &
                                 dir_clim,nyr_non_co2,                        &
                                 file_non_co2,dir_anom,                       &
                                 diff_frac_const_imogen

END MODULE imogen_anlg_vals
