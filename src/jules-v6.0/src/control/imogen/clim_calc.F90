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

SUBROUTINE clim_calc(                                                         &
  gpoints,wgen,mm,md,t_clim,sw_clim,lw_clim,                                  &
  pstar_ha_clim,rh15m_clim,rainfall_clim,snowfall_clim,                       &
  uwind_clim,vwind_clim,dtemp_clim,f_wet_clim,tmin_wg,tmax_wg,                &
  sw_wg,rh15m_wg,precip_wg,t_anom,sw_anom,lw_anom,                            &
  pstar_ha_anom,rh15m_anom,precip_anom,uwind_anom,vwind_anom,                 &
  dtemp_anom,sw_out,t_out,lw_out,conv_rain_out,conv_snow_out,                 &
  ls_rain_out,ls_snow_out,pstar_out,wind_out,qhum_out,nsdmax,                 &
  step_day,seed_rain,sec_day,lat,long)

USE missing_data_mod, ONLY: rmdi

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Calculates the daily meteorology based on the monthly climatology and 
!      anomaly patterns
!
! Code Owner: Please refer to ModuleLeaders.txt
!             This file belongs in IMOGEN
!
! Code Description:
!   Language: Fortran 90.
!   
!-----------------------------------------------------------------------------

INTEGER ::                                                                    &
  mm,                                                                         &
           !IN Number of months
  md,                                                                         &
           !IN Number of days in (GCM) month
  step_day,                                                                   &
           !IN Number of daily timesteps of IMPACTS_MODEL
  istep,                                                                      &
           !Looping parameter over suib-daily periods
  sec_day,                                                                    &
           !IN Number of seconds in each day
  nsdmax,                                                                     &
           !IN Maximum number of possible subdaily increments
  l,j,k,                                                                      &
           !Loop parameters
  gpoints,                                                                    &
           !Number of points, not including Antartica
  seed_rain(4)
           !WORK Seeding number for subdaily rainfall.

LOGICAL ::                                                                    &
  wgen     !IN Is the weather generator switched on.

REAL ::                                                                       &
  t_anom(gpoints,mm),                                                         &
           ! Temperature anomalies fro
  precip_anom(gpoints,mm),                                                    &
           ! Precip anomalies from AM
  rh15m_anom(gpoints,mm),                                                     &
           ! Relative humidity anomali
  uwind_anom(gpoints,mm),                                                     &
           ! u-wind anomalies from AM
  vwind_anom(gpoints,mm),                                                     &
           ! v-wind anomalies from AM
  dtemp_anom(gpoints,mm),                                                     &
           ! Diurnal Temperature (K)
  pstar_ha_anom(gpoints,mm),                                                  &
           ! Pressure anomalies from A
  sw_anom(gpoints,mm),                                                        &
           ! Shortwave radiation anoma
  lw_anom(gpoints,mm),                                                        &
           ! Longwave radiation anomal
  lat(gpoints),                                                               &
           ! Latitudinal position of l
  long(gpoints)
           ! Longitudinal position of

! Driving "control" climatology
REAL ::                                                                       &
  t_clim(gpoints,mm),                                                         &
           !IN Control climate temperature
  rainfall_clim(gpoints,mm),                                                  &
           !IN Control climate rainfall
  snowfall_clim(gpoints,mm),                                                  &
           !IN Control climate snowfall
  rh15m_clim(gpoints,mm),                                                     &
           !IN Control climate relative humidity
  uwind_clim(gpoints,mm),                                                     &
           !IN Control climate u-wind
  vwind_clim(gpoints,mm),                                                     &
           !IN Control climate v-wind
  dtemp_clim(gpoints,mm),                                                     &
           !IN Control climate diurnal Tem
  pstar_ha_clim(gpoints,mm),                                                  &
           !IN Control climate pressure
  sw_clim(gpoints,mm),                                                        &
           !IN Control climate shortwave radiation
  lw_clim(gpoints,mm),                                                        &
           !IN Control climate longwave radiation
  f_wet_clim(gpoints,mm)
           !IN Control climate fraction we

! Output from the weather generator when called
REAL ::                                                                       &
  precip_wg(gpoints,mm,md),                                                   &
           ! Daily precipitation
  tmin_wg(gpoints,mm,md),                                                     &
           ! Daily minimum temperature
  tmax_wg(gpoints,mm,md),                                                     &
           ! Daily maximum temperature
  sw_wg(gpoints,mm,md),                                                       &
           ! Daily shortwave radiation
  rh15m_wg(gpoints,mm,md)
           ! Daily relative humidity

! Create "local" (in time) values of the arrays below
REAL ::                                                                       &
  t_out_local(gpoints,nsdmax),                                                &
           !WORK temperature (K)
  conv_rain_out_local(gpoints,nsdmax),                                        &
           !WORK temperature (mm/day)
  ls_rain_out_local(gpoints,nsdmax),                                          &
           !WORK temperature (mm/day)
  ls_snow_out_local(gpoints,nsdmax),                                          &
           !WORK temperature (mm/day)
  qhum_out_local(gpoints,nsdmax),                                             &
           !WORK humidity (kg/kg)
  wind_out_local(gpoints,nsdmax),                                             &
           !WORK wind  (m/s)
  pstar_out_local(gpoints,nsdmax),                                            &
           !WORK pressure (Pa)
  sw_out_local(gpoints,nsdmax),                                               &
           !WORK shortwave radiation
  lw_out_local(gpoints,nsdmax)
           !WORK longwave radiation

! Create fine temperal resolution year of climatology (to be used by imp
! studies or DGVMs).
REAL ::                                                                       &
  t_out(gpoints,mm,md,nsdmax),                                                &
           !OUT Calculated temperature
  conv_rain_out(gpoints,mm,md,nsdmax),                                        &
           !OUT Calculated convective rainfall (mm/day)
  conv_snow_out(gpoints,mm,md,nsdmax),                                        &
           !OUT Calculated convective rainfall (mm/day)
  ls_rain_out(gpoints,mm,md,nsdmax),                                          &
           !OUT Calculated large scale rainfall (mm/day)
  ls_snow_out(gpoints,mm,md,nsdmax),                                          &
           !OUT Calculated large scale snowfall (mm/day)
  qhum_out(gpoints,mm,md,nsdmax),                                             &
           !OUT Calculated humidity
  wind_out(gpoints,mm,md,nsdmax),                                             &
           !OUT Calculated wind  (m/s)
  pstar_out(gpoints,mm,md,nsdmax),                                            &
           !OUT Calculated pressure
  sw_out(gpoints,mm,md,nsdmax),                                               &
           !OUT Calculated shortwave radiation
  lw_out(gpoints,mm,md,nsdmax)
           !OUT Calculated longwave radiation

! Variables for daily climatology
REAL ::                                                                       &
  t_daily(gpoints,mm,md),                                                     &
           !WORK Calculated temperature
  precip_daily(gpoints,mm,md),                                                &
           !WORK Calculated temperature
  uwind_daily(gpoints,mm,md),                                                 &
           !WORK Calculated "u"-wind
  vwind_daily(gpoints,mm,md),                                                 &
           !WORK Calculated "v"-wind
  wind_daily(gpoints,mm,md),                                                  &
           !WORK Calculated wind  (m/s)
  dtemp_daily(gpoints,mm,md),                                                 &
           !WORK Calculated diurnal Temperature
  pstar_daily(gpoints,mm,md),                                                 &
           !WORK Calculated pressure (Pa)
  sw_daily(gpoints,mm,md),                                                    &
           !WORK Calculated shortwave radiation
  lw_daily(gpoints,mm,md)
           !WORK Calculated longwave radiation

! And "local" (in time) values of the DAILY variables.
REAL ::                                                                       &
  t_daily_local(gpoints),                                                     &
           !WORK Calculated temperature (K)
  precip_daily_local(gpoints),                                                &
           !WORK Calculated precip
  rh15m_daily_local(gpoints),                                                 &
           !WORK Calculated humidity (%)
  wind_daily_local(gpoints),                                                  &
           !WORK Calculated wind  (m/s)
  dtemp_daily_local(gpoints),                                                 &
           !WORK Calculated diurnal Temperature
  pstar_daily_local(gpoints),                                                 &
           !WORK Calculated pressure (Pa)
  sw_daily_local(gpoints),                                                    &
           !WORK Calculated shortwave radiation
  lw_daily_local(gpoints)
           !WORK Calculated longwave radiation

REAL ::                                                                       &
  rh15m_daily(gpoints,mm,md)
           !WORK Calculated relative humidity

REAL ::                                                                       &
  tlocal,                                                                     &
           !WORK Local value of T_DAILY in
  plocal   !WORK Local value of PSTAR_DAIL

!---------------------------------------------------------------------
! Variables required to split rainfall up so that it rains roughly
! the correct no. of days/month when weather generator is switched off.
!    INTEGER ::                                                    &
!      NO_RAINDAY,                                                 &
               ! WORK No. of rainy days in the month.
!      INT_RAINDAY
               ! WORK Rain day interval.
!      IC_RAINDAY
               ! WORK No. of rain days counter.
!    REAL ::                                                       &
!      TOT_RAIN ! WORK Total rain counter.
!---------------------------------------------------------------------


! Calculate monthly means and add anomalies.
! Anomalies will be zero, if set ANOM=.F.
DO j = 1,mm                 !Loop over the months
  DO l = 1,gpoints         !Loop over land points
    DO k = 1,md           !Loop over the days

      IF (wgen) THEN
        !     Temperature (K)
        t_daily(l,j,k) = (0.5 *                                               &
                          (tmin_wg(l,j,k) + tmax_wg(l,j,k))                   &
                         ) + t_anom(l,j)

        !     Shortwave radiation (W/m2)
        sw_daily(l,j,k) = sw_wg(l,j,k) + sw_anom(l,j)

        !     Relative humidity (kg/kg)
        rh15m_daily(l,j,k) = rh15m_wg(l,j,k) + rh15m_anom(l,j)
        dtemp_daily(l,j,k) = (tmax_wg(l,j,k) - tmin_wg(l,j,k))                &
                           + dtemp_anom(l,j)

        !     Precip : At present, the Hadley Centre GCM AM has only PRECIP in t
        !     column (ie include snowfall).
        precip_daily(l,j,k) = precip_wg(l,j,k) + precip_anom(l,j)
      ELSE
        t_daily(l,j,k)     = t_clim(l,j) + t_anom(l,j)
        sw_daily(l,j,k)    = sw_clim(l,j) + sw_anom(l,j)
        rh15m_daily(l,j,k) = rh15m_clim(l,j) + rh15m_anom(l,j)
        dtemp_daily(l,j,k) = dtemp_clim(l,j) + dtemp_anom(l,j)

        precip_daily(l,j,k) = rainfall_clim(l,j)                              &
                            + snowfall_clim(l,j)                              &
                            + precip_anom(l,j)

        !CH-EDITTED OUT LINES BELOW BUT WILL BE RE_IMPLEMENTED WITH DOCUMENTATIO
        ! To correct the no. of rain days per month if WGEN is off:
        !                 IF (K == 1)THEN
        !                    NO_RAINDAY=NINT(F_WET_CLIM(L,J)*MD)
        !                    IF (NO_RAINDAY <  1)NO_RAINDAY=1
        !                    IF (NO_RAINDAY >  MD)NO_RAINDAY=MD
        !                    INT_RAINDAY=INT(FLOAT(MD)/FLOAT(NO_RAINDAY))
        !                    IF (INT_RAINDAY >  MD)INT_RAINDAY=MD
        !                    IC_RAINDAY=0
        !                    TOT_RAIN=0.0
        !                 END IF
        !
        !                 IF (MOD(K,INT_RAINDAY) == 0
        !    &                 .AND.IC_RAINDAY <  NO_RAINDAY)THEN
        !                    PRECIP_DAILY(L,J,K) = (RAINFALL_CLIM(L,J)+
        !    &                 SNOWFALL_CLIM(L,J)+PRECIP_ANOM(L,J))
        !    &                    *FLOAT(MD)/FLOAT(NO_RAINDAY)
        !                    IC_RAINDAY=IC_RAINDAY+1
        !                 ELSE
        !                    PRECIP_DAILY(L,J,K) = 0.0
        !                 END IF
        !                 TOT_RAIN=TOT_RAIN+PRECIP_DAILY(L,J,K)
        !CH-END OF EDITTING OUT
      END IF

      !     Make sure precip anomalies do not produce negative rainfall
      precip_daily(l,j,k) = MAX(precip_daily(l,j,k), 0.0)

      !     Pressure (Pa) - Include unit conversion from HPa to Pa
      pstar_daily(l,j,k) =                                                    &
                100.0 * (pstar_ha_clim(l,j) + pstar_ha_anom(l,j))

      !     Check on humidity bounds
      rh15m_daily(l,j,k) = MIN(rh15m_daily(l,j,k), 100.0)
      rh15m_daily(l,j,k) = MAX(rh15m_daily(l,j,k), 0.0)
      !     Now convert humidity units into required (g/kg)
      tlocal = t_daily(l,j,k)
      plocal = pstar_daily(l,j,k)

      !     Check to make sure anomalies do not produce negative values
      dtemp_daily(l,j,k) = MAX(dtemp_daily(l,j,k), 0.0)

      !     Longwave radiation (W/m2)
      lw_daily(l,j,k) = lw_clim(l,j) + lw_anom(l,j)

      !     Wind
      uwind_daily(l,j,k) = uwind_clim(l,j) + uwind_anom(l,j)
      vwind_daily(l,j,k) = vwind_clim(l,j) + vwind_anom(l,j)
      wind_daily(l,j,k)  = SQRT(                                              &
        (uwind_daily(l,j,k)**2) + (vwind_daily(l,j,k)**2)                     &
      )

      !     Check to make sure anomalies do not produce zero windspeed (ie
      !     below measurement level).
      wind_daily(l,j,k) = MAX(wind_daily(l,j,k), 0.01)
    END DO
  END DO
END DO

!     Disaggregate down to sub-daily (ie TSTEP values)
!     This calls subroutine DAY_CALC which for each day converts values
!     "_daily" to "_out".

!     Variables going in (with units) are SW(W/m2), Precip (mm/day), Tem
!     (K), DTemp (K), LW (W/m2), PSTAR(Pa), Wind (m/2) and QHUM (kg/kg)

!     Variables coming out are subdaily estimates of above variables, ex
!     for DTEMP (which no longer has meaning), and temperature dependent
!     the splitting of precipitation back into LS_CONV, LS_SNOW and
!     Convective rainfall (there is no convective snow, and so below, th
!     set of have a zero value).

DO j = 1,mm                !Loop over the months
  DO k = 1,md              !Loop over the days
    DO l = 1,gpoints
      sw_daily_local(l)     = sw_daily(l,j,k)
      precip_daily_local(l) = precip_daily(l,j,k)
      t_daily_local(l)      = t_daily(l,j,k)
      dtemp_daily_local(l)  = dtemp_daily(l,j,k)
      lw_daily_local(l)     = lw_daily(l,j,k)
      pstar_daily_local(l)  = pstar_daily(l,j,k)
      wind_daily_local(l)   = wind_daily(l,j,k)
      rh15m_daily_local(l)  = rh15m_daily(l,j,k)
    END DO

    CALL day_calc(                                                            &
      gpoints,step_day,md,sw_daily_local,precip_daily_local,                  &
      t_daily_local,dtemp_daily_local,lw_daily_local,                         &
      pstar_daily_local,wind_daily_local,rh15m_daily_local,                   &
      sw_out_local,t_out_local,lw_out_local,                                  &
      conv_rain_out_local,ls_rain_out_local,                                  &
      ls_snow_out_local,pstar_out_local,wind_out_local,                       &
      qhum_out_local,j,k,lat,long,sec_day,nsdmax,seed_rain                    &
    )

    ! Finalise value and set unused output values as rmdi as a precaution.
    DO l = 1,gpoints
      DO istep = 1,step_day
        sw_out(l,j,k,istep)        = sw_out_local(l,istep)
        t_out(l,j,k,istep)         = t_out_local(l,istep)
        lw_out(l,j,k,istep)        = lw_out_local(l,istep)
        conv_rain_out(l,j,k,istep) =                                          &
                                 conv_rain_out_local(l,istep)
        conv_snow_out(l,j,k,istep) = 0.0
        ls_rain_out(l,j,k,istep)   = ls_rain_out_local(l,istep)
        ls_snow_out(l,j,k,istep)   = ls_snow_out_local(l,istep)
        pstar_out(l,j,k,istep)     = pstar_out_local(l,istep)
        wind_out(l,j,k,istep)      = wind_out_local(l,istep)
        qhum_out(l,j,k,istep)      = qhum_out_local(l,istep)
      END DO

      DO istep = step_day+1,nsdmax
        sw_out(l,j,k,istep)        = rmdi
        t_out(l,j,k,istep)         = rmdi
        lw_out(l,j,k,istep)        = rmdi
        conv_rain_out(l,j,k,istep) = rmdi
        conv_snow_out(l,j,k,istep) = rmdi
        ls_rain_out(l,j,k,istep)   = rmdi
        ls_snow_out(l,j,k,istep)   = rmdi
        pstar_out(l,j,k,istep)     = rmdi
        wind_out(l,j,k,istep)      = rmdi
        qhum_out(l,j,k,istep)      = rmdi
      END DO
    END DO
  END DO                  !End of loop over days
END DO                    !End of loop over months

RETURN

END SUBROUTINE clim_calc
#endif
