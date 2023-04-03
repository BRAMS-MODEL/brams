!******************************COPYRIGHT**************************************
! (c) Centre for Ecology and Hydrology. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237] 
!******************************COPYRIGHT**************************************

MODULE imogen_progs

USE imogen_constants, ONLY: n_olevs,nfarray

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Module for the prognostic variables required for JULES-IMOGEN simulations
!
! Code Owner: Please refer to ModuleLeaders.txt
!             This file belongs in IMOGEN
!
! Code Description:
!   Language: Fortran 90.
!   
!-----------------------------------------------------------------------------

REAL ::                                                                       &
  co2_ppmv,                                                                   &
            ! Atmospheric CO2 concentration (ppmv)
  co2_start_ppmv,                                                             &
            ! Atmospheric CO2 concentration at start of year (ppmv)
  co2_change_ppmv,                                                            &
            ! Change in CO2 between restarts (ppmv)
  ch4_ppbv
            ! Atmospheric CH2 concentration (ppbv)

REAL ::                                                                       &
  dtemp_o(n_olevs),                                                           &
            ! Ocean temperature anomalies (K)
  fa_ocean(nfarray),                                                          &
            ! CO2 fluxes from the atmosphere to the
            ! ocean (ie positive upwards) (ppm/m2/yr)
  d_land_atmos_co2,                                                           &
            ! Change in atmospheric CO2 concentration due to
            ! land co2 feedbacks (ppm/year).
  d_land_atmos_ch4,                                                           &
            ! Change in global total CH4 flux from the land to the atmosphere 
            ! w.r.t a specified reference flux, fch4_ref (kgC/year)
  d_ocean_atmos,                                                              &
            ! Change in atmospheric CO2 concentration due to
            ! ocean feedbacks (ppm/year)
  c_emiss_out
            ! Prescribed emissions in Gt C per year

INTEGER ::                                                                    &
  seed_wg(4)
            ! Seeding number for weather generator
            ! random number

END MODULE imogen_progs
