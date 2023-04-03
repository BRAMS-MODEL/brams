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

SUBROUTINE drdat(                                                             &
  iyear,gpoints,t_anom_dat,precip_anom_dat,rh15m_anom_dat,                    &
  uwind_anom_dat,vwind_anom_dat,dtemp_anom_dat,                               &
  pstar_ha_anom_dat,sw_anom_dat,lw_anom_dat,dir_anom,                         &
  longmin_dat,latmin_dat,longmax_dat,latmax_dat,mm,drive_month                &
)

USE io_constants, ONLY: imogen_unit

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads prescribed driving data for the IMOGEN simulations
!
! Code Owner: Please refer to ModuleLeaders.txt
!             This file belongs in IMOGEN
! Written by: P. Cox (1997) and C. Huntingford (2004)
! Code Description:
!   Language: Fortran 90.
!   
!-----------------------------------------------------------------------------

INTEGER ::                                                                    &
  iyear,                                                                      &
              ! WORK Year of interest.
  gpoints,                                                                    &
              ! IN Number of land points.
  mm,                                                                         &
              ! IN Number months in year.
  im          ! WORK Month loop parameter

CHARACTER(LEN=4) ::                                                           &
  drive_year,                                                                 &
              ! IN Label for year of driving data.
  drive_month(12)
              ! Labels for month of driving data.

CHARACTER(LEN=180) ::                                                         &
  dir_anom,                                                                   &
              ! IN Directory containing anomalies.
  driver_anom
              ! IN Directory containing anomalies.

!-----------------------------------------------------------------------
! Climatological forcing variables.
!-----------------------------------------------------------------------
REAL ::                                                                       &
  latmin_dat,latmax_dat,                                                      &
             ! IN Latitudinal limits of the area
  longmin_dat,longmax_dat,                                                    &
             ! IN Longitudinal limits of the area
  lat(gpoints),                                                               &
             ! WORK Latitude (degrees).
  long(gpoints),                                                              &
             ! WORK Longitude (degrees).
  dtemp_anom_dat(gpoints,mm),                                                 &
             ! OUT Diurnal temperature range (K).
  lw_anom_dat(gpoints,mm),                                                    &
             ! OUT Downward surface longwave radiation (W/m).
  pstar_ha_anom_dat(gpoints,mm),                                              &
             ! OUT Surface pressure (Pa).
  rainfall_anom_dat(gpoints,mm),                                              &
             ! WORK Rainfall rate (mm/day).
  rh15m_anom_dat(gpoints,mm),                                                 &
             ! OUT Relative humidity at 1.5m (%).
  snowfall_anom_dat(gpoints,mm),                                              &
             ! WORK Snowfall rate (mm/day).
  sw_anom_dat(gpoints,mm),                                                    &
             ! OUT Downward surface shortwave radiation (W/m2).
  t_anom_dat(gpoints,mm),                                                     &
             ! OUT Air temperature (K).
  precip_anom_dat(gpoints,mm),                                                &
             ! OUT Precipitation (mm/day)
  uwind_anom_dat(gpoints,mm),                                                 &
             ! OUT Windspeed u-components (m/s).
  vwind_anom_dat(gpoints,mm)
             ! OUT Windspeed v-components (m/s).

!-----------------------------------------------------------------------
! Loop counters.
!-----------------------------------------------------------------------
INTEGER :: l ! Loop counters

!-----------------------------------------------------------------------
! Convert year to string.
!-----------------------------------------------------------------------
WRITE(drive_year,'(I4)') iyear

DO im = 1,mm
  driver_anom = TRIM(dir_anom) // drive_year // drive_month(mm)

  OPEN(imogen_unit, FILE=driver_anom,                                         &
                    STATUS='old', POSITION='rewind', ACTION='read')

  READ(imogen_unit,*) longmin_dat,latmin_dat,longmax_dat,latmax_dat
  DO l = 1,gpoints
    READ(imogen_unit,*) long(l),lat(l),t_anom_dat(l,im),                      &
                        rh15m_anom_dat(l,im),uwind_anom_dat(l,im),            &
                        vwind_anom_dat(l,im),lw_anom_dat(l,im),               &
                        sw_anom_dat(l,im),dtemp_anom_dat(l,im),               &
                        rainfall_anom_dat(l,im),                              &
                        snowfall_anom_dat(l,im),pstar_ha_anom_dat(l,im)

    precip_anom_dat(l,im) = rainfall_anom_dat(l,im)                           &
                          + snowfall_anom_dat(l,im)
  END DO
END DO
CLOSE(imogen_unit)

RETURN

END SUBROUTINE drdat
#endif
