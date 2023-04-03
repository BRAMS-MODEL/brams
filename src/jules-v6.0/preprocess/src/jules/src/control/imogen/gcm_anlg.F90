!******************************COPYRIGHT**************************************
! (c) Centre for Ecology and Hydrology. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237] 
!******************************COPYRIGHT**************************************

SUBROUTINE gcm_anlg(                                                          &
  q,land_pts,t_anom_am,precip_anom_am,rh15m_anom_am,                          &
  uwind_anom_am,vwind_anom_am,dtemp_anom_am,pstar_ha_anom_am,                 &
  sw_anom_am,lw_anom_am,n_olevs,dir_patt,f_ocean,kappa_o,                     &
  lambda_l,lambda_o,mu,dtemp_o,longmin_am,latmin_am,longmax_am,               &
  latmax_am,mm                                                                &
)

USE io_constants, ONLY: imogen_unit

USE imogen_map, ONLY: sgindinv

USE imogen_constants, ONLY: drive_month,n_imogen_land

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   GCM Analouge Model. Uses pattern scaling to supply meteorological anomalies
!   based on the global temperature anomaly.
!
! Code Owner: Please refer to ModuleLeaders.txt
!             This file belongs in IMOGEN
!
! Written by: C.Huntingford & P.Cox, 1998
! Code Description:
!   Language: Fortran 90.
!   
!-----------------------------------------------------------------------------

INTEGER ::                                                                    &
  land_pts,                                                                   &
                 ! IN Number of land points.
  im,                                                                         &
                 ! IN Loop over months
  n_olevs,                                                                    &
                 ! IN Number of ocean thermal layers.
  mm             ! In Number of months in a year

CHARACTER(LEN=180) ::                                                         &
  dir_patt,                                                                   &
                 ! IN Directory containing anomaly patterns.
  driver_patt    ! WORK File containing anomaly
                 !      patterns from analogue model.

REAL ::                                                                       &
  lat_am(land_pts),                                                           &
                  ! Latitude read from file.
  long_am(land_pts)
                  !Longitude read from file.

!-----------------------------------------------------------------
! Climatological forcing variables.
!-----------------------------------------------------------------
REAL ::                                                                       &
  f_ocean,                                                                    &
                 ! IN Fractional coverage of the ocean.
  kappa_o,                                                                    &
                 ! IN Ocean eddy diffusivity (W/m/K).
  lambda_l,                                                                   &
                 ! IN Inverse climate sensitivity over
                 !    land (W/m2/K).
  lambda_o,                                                                   &
                 ! IN Inverse climate sensitivity over
                 !    ocean (W/m2/K).
  mu,                                                                         &
                 ! IN Ratio of land to ocean
                 !    temperature anomalies.
  dtemp_o(n_olevs),                                                           &
                 ! INOUT Ocean mean temperature anomaly (K).
  dtemp_l,                                                                    &
                 ! OUT Land mean temperature anomaly (K).
  dtemp_anom_am(land_pts,mm),                                                 &
                 ! OUT Diurnal temperature range (K).
  lw_anom_am(land_pts,mm),                                                    &
                 ! OUT Downward surface longwave
                 !     radiation anomaly (W/m).
  rainfall_anom_am(land_pts,mm),                                              &
                 ! OUT Rainfall rate anomaly (mm/day)
  rh15m_anom_am(land_pts,mm),                                                 &
                 ! OUT Relative humidity at 1.5m anom
  snowfall_anom_am(land_pts,mm),                                              &
                 ! OUT Snowfall rate anomaly (mm/day)
  sw_anom_am(land_pts,mm),                                                    &
                 ! OUT Downward surface shortwave
                 !     radiation anomaly (W/m2).
  t_anom_am(land_pts,mm),                                                     &
                 ! OUT Air temperature anomaly (K).
  precip_anom_am(land_pts,mm),                                                &
                 ! OUT Precipitation (snowfall
                 !     plus rainfall) anomaly (mm/day)
  uwind_anom_am(land_pts,mm),                                                 &
                 ! OUT Wind speed anomaly (m/s).
  vwind_anom_am(land_pts,mm),                                                 &
                 ! OUT Wind speed anomaly (m/s).
  latmin_am,latmax_am,                                                        &
                 ! WORK Latitudinal limits of the area
                 !      (degrees).
  longmin_am,longmax_am,                                                      &
                 ! WORK Longitudinal limits of the area
                 !      (degrees).
  pstar_ha_anom_am(land_pts,mm),                                              &
                 ! OUT Surface pressure (hPa).
  q
                 ! WORK Increase in radiative forcing (W/m2)

!-----------------------------------------------------------------
! Anomaly patterns scaled to land mean temperature anomalies.
!-----------------------------------------------------------------
REAL ::                                                                       &
  ddtemp_day_pat,                                                             &
                 ! WORK Diurnal temperature range (.).
  dlw_c_pat,                                                                  &
                 ! WORK Downward surface longwave
                 !      radiation (W/m2/K).
  dpstar_c_pat,                                                               &
                 ! WORK Surface pressure (hPa/K).
  drainfall_pat,                                                              &
                 ! WORK Rainfall rate (mm/day/K).
  drh15m_pat,                                                                 &
                 ! WORK Relative humidity at 1.5m (%/K).
  dsnowfall_pat,                                                              &
                 ! WORK Snowfall rate (mm/day/K).
  dsw_c_pat,                                                                  &
                 ! WORK Surface shortwave radiation (W/m2/K).
  dt_c_pat,                                                                   &
                 ! WORK Air temperature (.).
  duwind_pat,dvwind_pat
                 ! WORK Windspeed components (m/s/K).

!-----------------------------------------------------------------
! Loop counters.
!-----------------------------------------------------------------
INTEGER ::                                                                    &
  i,l    ! WORK


!-----------------------------------------------------------------
! Loop over months
!-----------------------------------------------------------------
DO im = 1,mm

  !-----------------------------------------------------------------
  ! Calculate new area mean temperature anomalies
  !-----------------------------------------------------------------
  IF (im == 1) THEN
    CALL delta_temp(                                                          &
      n_olevs,f_ocean,kappa_o,lambda_l,lambda_o,mu,q,                         &
      dtemp_l,dtemp_o                                                         &
    )
  END IF

  !-----------------------------------------------------------------
  ! Define the anomaly patterns and read the header
  !-----------------------------------------------------------------
  driver_patt = TRIM(dir_patt) // drive_month(im)

  OPEN(imogen_unit, FILE=driver_patt,                                         &
                    STATUS='old', POSITION='rewind', ACTION='read')

  READ(imogen_unit,*) longmin_am,latmin_am,longmax_am,latmax_am

  !-----------------------------------------------------------------
  ! Read in initial climatology and then define the new climate data.
  !-----------------------------------------------------------------
  DO i = 1,n_imogen_land
    IF (sgindinv(i) > 0) THEN
      l = sgindinv(i)

      READ(imogen_unit,*)                                                     &
        long_am(l),lat_am(l),dt_c_pat,drh15m_pat,duwind_pat,                  &
        dvwind_pat,dlw_c_pat,dsw_c_pat,ddtemp_day_pat,                        &
        drainfall_pat,dsnowfall_pat,dpstar_c_pat

      t_anom_am(l,im) = dt_c_pat * dtemp_l
      rh15m_anom_am(l,im) = drh15m_pat * dtemp_l
      uwind_anom_am(l,im) = duwind_pat * dtemp_l
      vwind_anom_am(l,im) = dvwind_pat * dtemp_l
      rainfall_anom_am(l,im) = drainfall_pat * dtemp_l
      snowfall_anom_am(l,im) = dsnowfall_pat * dtemp_l
      dtemp_anom_am(l,im) = ddtemp_day_pat * dtemp_l
      lw_anom_am(l,im) = dlw_c_pat * dtemp_l
      sw_anom_am(l,im) = dsw_c_pat * dtemp_l
      pstar_ha_anom_am(l,im) = dpstar_c_pat * dtemp_l

      precip_anom_am(l,im) = rainfall_anom_am(l,im) +                         &
                             snowfall_anom_am(l,im)
    ELSE
      READ(imogen_unit,*)
    END IF
  END DO     !End of loop over land points

  CLOSE(imogen_unit)
END DO     !End of loop over months

RETURN
END SUBROUTINE gcm_anlg
