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

SUBROUTINE imogen_update_clim(progs)

USE model_time_mod, ONLY: current_time

USE ancil_info, ONLY: dim_cs1, dim_cslayer

USE io_constants, ONLY: imogen_unit

USE datetime_mod, ONLY: secs_in_day

USE aero, ONLY: co2_mmr

USE ancil_info, ONLY: land_pts, nsoilt

USE imogen_progs, ONLY:                                                       &
  co2_ppmv,co2_start_ppmv,dtemp_o

USE imogen_run, ONLY:                                                         &
  include_co2,c_emissions,file_scen_co2_ppmv,anom,anlg,                       &
  co2_init_ppmv,include_non_co2,file_non_co2_vals,wgen,                       &
  land_feed_co2, land_feed_ch4

USE imogen_time, ONLY:                                                        &
  mm,md,nsdmax,step_day

USE imogen_anlg_vals, ONLY:                                                   &
  q2co2,nyr_non_co2,file_non_co2,dir_patt,f_ocean,kappa_o,                    &
  lambda_l,lambda_o,mu,dir_anom

USE imogen_constants, ONLY:                                                   &
  n_olevs, drive_month

USE imogen_clim, ONLY:                                                        &
!   Scalars
    latmin_clim,latmax_clim,longmin_clim,longmax_clim,                        &
!   Arrays
    t_clim,rainfall_clim,snowfall_clim,rh15m_clim,uwind_clim,                 &
    vwind_clim,dtemp_clim,pstar_ha_clim,sw_clim,lw_clim,                      &
    f_wet_clim, lat, long, dctot_co2

USE imogen_drive_vars, ONLY:                                                  &
  t_out,conv_rain_out,conv_snow_out,ls_rain_out,ls_snow_out,                  &
  qhum_out,wind_out,pstar_out,sw_out,lw_out

USE logging_mod, ONLY: log_fatal

USE jules_print_mgr, ONLY:                                                    &
  jules_message,                                                              &
  jules_print

!TYPE definitions
USE prognostics, ONLY: progs_type
USE jules_fields_mod, ONLY: trifctltype

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Updates the climatology and meteorology based on the global temperature
!   anomaly. The Global temperautre anomaly is calculated within as the a
!   function of the atmospheric composition, i.e. the CO2 concentration plus the
!   non-CO2 components.
!
! Code Owner: Please refer to ModuleLeaders.txt
!             This file belongs in IMOGEN
!
! Code Description:
!   Language: Fortran 90.
!
!-----------------------------------------------------------------------------

! Arguments

!TYPES containing the data
TYPE(progs_type), INTENT(IN OUT) :: progs

INTEGER :: l, n, nn, m ! Loop counter

INTEGER ::                                                                    &
  tally_co2_file,                                                             &
                 !WORK If used, checks CO2 value available in
                 !     file "FILE_SCEN_CO2_PPMV"
  yr_co2_file,                                                                &
                 !WORK If used, reads in years available in
                 !     file "FILE_SCEN_CO2_PPMV"
  kode           ! Used to hold IOSTAT while reading file

REAL ::                                                                       &
  co2_file_ppmv  ! If used, is prescribed CO2 value for read
                 ! in years in "FILE_SCEN_CO2_PPMV"

REAL ::                                                                       &
  q_co2,                                                                      &
               !WORK Radiative forcing due to CO2 (W/m2)
  q_ch4,                                                                      &
               !WORK Radiative forcing due to changes in CH4 concentration
               !     from reference projection (W/m2)
  q_non_co2,                                                                  &
               !WORK Radiative forcing due to non-CO2 (W/m2)
  q            !WORK Total radiative forcing, both CO2 and NON CO2

REAL :: latmin,latmax,longmin,longmax
               !WORK Max and min lat and long for files

! Dummy variables for weather generator which is not available in
! this version
REAL ::                                                                       &
  precip_wg(land_pts,mm,md),                                                  &
  tmin_wg(land_pts,mm,md),                                                    &
  tmax_wg(land_pts,mm,md),                                                    &
  sw_wg(land_pts,mm,md),                                                      &
  rh15m_wg(land_pts,mm,md)

!-----------------------------------------------------------------
! Variables to hold calculated anomalies
!-----------------------------------------------------------------
REAL ::                                                                       &
  t_anom(land_pts,mm),                                                        &
                ! WORK Temperature anomalies (K)
  precip_anom(land_pts,mm),                                                   &
                ! WORK Precip anomalies (mm/day)
  rh15m_anom(land_pts,mm),                                                    &
                ! WORK Relative humidity anomalies
  uwind_anom(land_pts,mm),                                                    &
                ! WORK u-wind anomalies (m/s)
  vwind_anom(land_pts,mm),                                                    &
                ! WORK v-wind anomalies (m/s)
  dtemp_anom(land_pts,mm),                                                    &
                ! WORK Diurnal Temperature (K)
  pstar_ha_anom(land_pts,mm),                                                 &
                ! WORK Pressure anomalies (hPa)
  sw_anom(land_pts,mm),                                                       &
                ! WORK Shortwave radiation anomalie
  lw_anom(land_pts,mm)
                ! WORK Longwave radiation anomalies

!------------------------------------------------------------------------
! Initialisation
q_co2     = 0.0
q_non_co2 = 0.0
q_ch4     = 0.0

t_anom(:,:)        = 0.0
precip_anom(:,:)   = 0.0
rh15m_anom(:,:)    = 0.0
uwind_anom(:,:)    = 0.0
vwind_anom(:,:)    = 0.0
dtemp_anom(:,:)    = 0.0
pstar_ha_anom(:,:) = 0.0
sw_anom(:,:)       = 0.0
lw_anom(:,:)       = 0.0

WRITE(jules_message,*) 'Updating IMOGEN climate'
CALL jules_print('imogen_update_clim',jules_message)


! Capture CO2 concentration at beginning of year
IF (include_co2) THEN
  co2_start_ppmv = co2_ppmv
END IF


! Hydrology 20th Century simulations (note also check for this run in
! subroutine IMOGEN_CONFIRMED_RUN which includes more stringent checks)
! OR analogue model simulations with CO2 prescribed.
IF (( .NOT. c_emissions) .AND. include_co2) THEN
  ! This works by reading in a file of CO2 concentrations, and checks that
  ! year is represented
  OPEN(imogen_unit, FILE=file_scen_co2_ppmv,                                  &
                      STATUS='old', POSITION='rewind', ACTION='read')

  tally_co2_file = 0
  kode = 0

  DO WHILE ( .TRUE.)
    READ(imogen_unit,FMT=*,IOSTAT = kode) yr_co2_file,co2_file_ppmv
    ! Check for end of file (or just any error) and exit loop if found
    ! iostat > 0 means illegal data
    ! iostat < 0 means end of record or end of file
    IF (kode /= 0) EXIT

    IF (yr_co2_file == current_time%year) THEN
      co2_ppmv = co2_file_ppmv
      tally_co2_file = tally_co2_file + 1
      ! We have found the correct year so exit loop
      EXIT
    END IF
  END DO

  CLOSE(imogen_unit)

  ! Check that value has been found.
  IF (tally_co2_file /= 1)                                                    &
    CALL log_fatal("imogen_update_clim", 'CO2 value not found in file')
END IF

! Now calculate the added monthly anomalies, either from analogue model
! or prescribed directly.
IF (anom) THEN
  IF (anlg) THEN
    ! This call is to the GCM analogue model. It is prescribed CO2
    ! concentration, and calculates non-CO2, and returns total change in
    ! radiative forcing, Q. Recall that the AM has a "memory" through
    ! DTEMP_0 - ie the ocean temperatures. Note that in this version of the
    ! code, the AM is updated yearly.

    ! Calculate the CO2 forcing
    IF (include_co2) THEN
      CALL radf_co2(co2_ppmv,co2_init_ppmv,q2co2,q_co2)
    END IF

    ! Calculate the non CO2 forcing
    IF (include_non_co2) THEN
      CALL radf_non_co2(                                                      &
        current_time%year,q_non_co2,nyr_non_co2,file_non_co2,                 &
        file_non_co2_vals                                                     &
      )
    END IF

    IF (land_feed_ch4) THEN
      CALL radf_ch4(q_ch4)
    END IF
    ! Calculate the total forcing
    q = q_co2 + q_non_co2 + q_ch4

    ! Call the GCM analogue model that responds to this forcing
    IF (include_co2 .OR. include_non_co2) THEN
      CALL gcm_anlg(                                                          &
        q,land_pts,t_anom,precip_anom,rh15m_anom,uwind_anom,                  &
        vwind_anom,dtemp_anom,pstar_ha_anom,sw_anom,lw_anom,                  &
        n_olevs,dir_patt,f_ocean,kappa_o,lambda_l,lambda_o,                   &
        mu,dtemp_o,longmin,latmin,longmax,latmax,mm                           &
      )

      ! Check driving files are compatible.
      IF ((ABS(longmin_clim - longmin) >= 1.0e-6) .OR.                        &
         (ABS(latmin_clim - latmin) >= 1.0e-6) .OR.                           &
         (ABS(longmax_clim - longmax) >= 1.0e-6) .OR.                         &
         (ABS(latmax_clim - latmax) >= 1.0e-6)) THEN
        CALL log_fatal("imogen_update_clim",                                  &
                       'Driving files are incompatible')
      END IF
    END IF

    ! Option where anomalies are prescribed directly not using the analogue
  ELSE
    CALL drdat(                                                               &
      current_time%year,land_pts,t_anom,precip_anom,rh15m_anom,               &
      uwind_anom,vwind_anom,dtemp_anom,pstar_ha_anom,sw_anom,                 &
      lw_anom,dir_anom,longmin,latmin,longmax,latmax,mm,                      &
      drive_month                                                             &
    )

    ! Check driving files are compatible.
    IF ((ABS(longmin_clim - longmin) >= 1.0e-6) .OR.                          &
       (ABS(latmin_clim - latmin) >= 1.0e-6) .OR.                             &
       (ABS(longmax_clim - longmax) >= 1.0e-6) .OR.                           &
       (ABS(latmax_clim - latmax) >= 1.0e-6)) THEN
      CALL log_fatal("imogen_update_clim",                                    &
                     'Driving files are incompatible')
    END IF
  END IF
ELSE
  ! Set anomalies to zero.
  t_anom(:,:) = 0.0
  sw_anom(:,:) = 0.0
  lw_anom(:,:) = 0.0
  pstar_ha_anom(:,:) = 0.0
  rh15m_anom(:,:) = 0.0
  precip_anom(:,:) = 0.0
  uwind_anom(:,:) = 0.0
  vwind_anom(:,:) = 0.0
  dtemp_anom(:,:) = 0.0
END IF         !End of where anomalies are calculated

! Now incorporate anomalies with climate data.
! At this point, have climatology, WG if called and anomalies of either
! "_AM" or "_DAT". Now calculate the daily values of the driving data.
! This is calculated using subroutine CLIM_CALC
CALL clim_calc(                                                               &
  land_pts,wgen,mm,md,t_clim,sw_clim,lw_clim,                                 &
  pstar_ha_clim,rh15m_clim,rainfall_clim,snowfall_clim,                       &
  uwind_clim,vwind_clim,dtemp_clim,f_wet_clim,tmin_wg,                        &
  tmax_wg,sw_wg,rh15m_wg,precip_wg,t_anom,sw_anom,lw_anom,                    &
  pstar_ha_anom,rh15m_anom,precip_anom,uwind_anom,                            &
  vwind_anom,dtemp_anom,sw_out,t_out,lw_out,conv_rain_out,                    &
  conv_snow_out,ls_rain_out,ls_snow_out,pstar_out,wind_out,                   &
  qhum_out,nsdmax,step_day,progs%seed_rain,secs_in_day,lat,long)

! Compute current land carbon as dctot_co2 for use in imogen_update_carb.
IF (land_feed_co2) THEN
  DO l = 1,land_pts
    dctot_co2(l) = trifctltype%cv_gb(l)
  END DO

  DO l = 1,land_pts
    DO m = 1,nsoilt
      DO n = 1,dim_cs1
        DO nn = 1,dim_cslayer
          dctot_co2(l) = dctot_co2(l) + progs%cs_pool_soilt(l,m,nn,n)
        END DO
      END DO
    END DO
  END DO
END IF

!Unit conversion ppm to mmr (g/g): mol mass co2 / mol mass dry air * 1e-6
co2_mmr = co2_ppmv * 44.0 / 28.97 * 1.0e-6


RETURN

END SUBROUTINE imogen_update_clim
#endif
