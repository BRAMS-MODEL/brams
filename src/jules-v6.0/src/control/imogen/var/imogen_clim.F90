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

MODULE imogen_clim

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Module for the climate data for JULES-IMOGEN simulations
!
! Code Owner: Please refer to ModuleLeaders.txt
!             This file belongs in IMOGEN
!
! Code Description:
!   Language: Fortran 90.
!   
!-----------------------------------------------------------------------------


! Maximum and minimum latitude and longitude for the climatology
! variables
REAL :: latmin_clim,latmax_clim,longmin_clim,longmax_clim

!-----------------------------------------------------------------
! Driving "control" climatology
!-----------------------------------------------------------------
REAL, ALLOCATABLE ::                                                          &
  t_clim(:,:),                                                                &
              ! Control climate temperature
  rainfall_clim(:,:),                                                         &
              ! Control climate rainfall
  snowfall_clim(:,:),                                                         &
              ! Control climate snowfall
  rh15m_clim(:,:),                                                            &
              ! Control climate relative humidity
  uwind_clim(:,:),                                                            &
              ! Control climate u-wind
  vwind_clim(:,:),                                                            &
              ! Control climate v-wind
  dtemp_clim(:,:),                                                            &
              ! Control climate diurnal Temp
  pstar_ha_clim(:,:),                                                         &
              ! Control climate pressure
  sw_clim(:,:),                                                               &
              ! Control climate shortwave radiation
  lw_clim(:,:),                                                               &
              ! Control climate longwave radiation
  f_wet_clim(:,:)
              ! Control climate fraction we

REAL, ALLOCATABLE ::                                                          &
  lat(:),                                                                     &
              ! Latitudinal position of land
  long(:),                                                                    &
              ! Longitudinal position of land
  dctot_co2(:),                                                               &
              ! Change in total surface gridbox CO2
              ! content during a period "YEAR_CO2" (kg C/m2)
  dctot_ch4(:)
              ! Total surface gridbox CH4 flux to the atmosphere
              ! during a period "YEAR_CO2" (kg C/m2)

END MODULE imogen_clim
#endif
