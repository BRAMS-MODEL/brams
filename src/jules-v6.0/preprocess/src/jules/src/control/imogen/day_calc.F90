!******************************COPYRIGHT**************************************
! (c) Centre for Ecology and Hydrology. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237] 
!******************************COPYRIGHT**************************************
SUBROUTINE day_calc(                                                          &
  pointsm,step_day,day_mon,sw_l,precip_l,t_l,dtemp_day_l,lw_l,                &
  pstar_l,wind_l,rh15m_l,sw_sd,t_sd,lw_sd,conv_rain_sd,                       &
  ls_rain_sd,ls_snow_sd,pstar_sd,wind_sd,qhum_sd,month,iday,                  &
  lat,long,sec_day,nsdmax,seed_rain                                           &
)

USE conversions_mod, ONLY: rhour_per_day, rsec_per_hour, pi
USE update_mod, ONLY: dur_ls_rain,dur_conv_rain,dur_ls_snow
USE datetime_mod, ONLY: secs_in_hour

USE qsat_mod, ONLY: qsat

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   This routine calculates sub-daily variability
!
! Code Owner: Please refer to ModuleLeaders.txt
!             This file belongs in IMOGEN
!
! Written by: C. Huntingford (April 2001) - based on earlier version by P. Cox
! Code Description:
!   Language: Fortran 90.
!   
!-----------------------------------------------------------------------------

INTEGER ::                                                                    &
  day_mon,                                                                    &
          !IN "Number of days in a month" (day)
  pointsm,                                                                    &
          !IN Maximum number of points in grid.
  step_day,                                                                   &
          !IN WORK Calculated number of timesteps per day
  sec_day,                                                                    &
          !IN Number of seconds in day (sec)
  nsdmax,                                                                     &
          !IN Maximum possible number of sub-daily timesteps.
  seed_rain(4),                                                               &
          !IN Seeding numbers required to disaggregate
          !   the rainfall
  n_tally !WORK Counting up number of precipitation periods.

!-----------------------------------------------------------------------
! Single day arrays of incoming variables
!-----------------------------------------------------------------------
REAL ::                                                                       &
  sw_l(pointsm),                                                              &
          ! IN Daily values for downward shortwave
          !    radiation (W/m2)
  precip_l(pointsm),                                                          &
          ! In Daily values of rain+snow (mm/day)
  t_l(pointsm),                                                               &
          ! IN Daily values of temperature (K)
  dtemp_day_l(pointsm),                                                       &
          ! IN Daily values of diurnal temperature range (K)
  lw_l(pointsm),                                                              &
          ! IN Daily values of surface longwave
          !    radiation (W/m**2).
  pstar_l(pointsm),                                                           &
          ! IN Daily values of pressure at level 1 (Pa).
  wind_l(pointsm),                                                            &
          ! IN Daily values of wind speed (m/s)
  rh15m_l(pointsm)
          ! Humidity (%)

!-----------------------------------------------------------------------
! Single day (global) arrays of variables above, but for up to hourly
! periods. NOTE: DTEMP_DAY_SD does not exist as T_SD combines T_L
! and DTEMP_DAY_L
!-----------------------------------------------------------------------

REAL ::                                                                       &
  sw_sd(pointsm,nsdmax),                                                      &
          ! OUT Sub-daily values for downward
          !     shortwave radiation (W/m2)
  t_sd(pointsm,nsdmax),                                                       &
          ! OUT Sub-daily values of temperature (K)
  lw_sd(pointsm,nsdmax),                                                      &
          ! OUT Sub-daily values of surface
          !     longwave radiation (W/m**2).
  pstar_sd(pointsm,nsdmax),                                                   &
          ! OUT Sub-daily values of pressure at level 1 (Pa).
  wind_sd(pointsm,nsdmax),                                                    &
          ! OUT Sub-daily values of wind speed (m/s)
  qhum_sd(pointsm,nsdmax),                                                    &
          ! OUT Sub_daily humidity deficit calculated
          !     from RH15M (kg/kg)
  conv_rain_sd(pointsm,nsdmax),                                               &
          ! OUT Sub-daily convective rain (mm/day)
  ls_rain_sd(pointsm,nsdmax),                                                 &
          ! OUT Sub-daily large scale rain (mm/day)
  ls_snow_sd(pointsm,nsdmax)
          ! OUT Sub-daily large scale snow (mm/day)

INTEGER ::                                                                    &
  n_event(pointsm,nsdmax),                                                    &
          ! 1: if rains/snows during timestep period
          ! 0: otherwise
  n_event_local(nsdmax)
          ! As N_EVENT, but for each gridpoint

REAL ::                                                                       &
  prec_loc(nsdmax)
          ! Temporary value of rainfall for each gridbox.

REAL ::                                                                       &
  t_sd_local(pointsm),                                                        &
          ! WORK Intermediate calculation of temperature (K)
  pstar_sd_local(pointsm),                                                    &
          ! WORK Intermediate calculation of pressure (Pa)
  qs_sd_local(pointsm)
          ! WORK Saturated humidity deficit associated with
          !      T_SD_LOCAL and PSTAR_SD_LOCAL

INTEGER ::                                                                    &
  istep,                                                                      &
          ! WORK Loop over sub-daily periods
  month,                                                                      &
          ! IN Month of interest
  iday,                                                                       &
          ! IN Day number since beginning of month
  daynumber,                                                                  &
          ! WORK Daynumber since beginning of year
  l       ! WORK Loop over land points

REAL ::                                                                       &
  sun(pointsm,nsdmax),                                                        &
          ! WORK Normalised solar radiation
  time_max(pointsm),                                                          &
          ! WORK Time (UTC) at which temperature is maximum (hours)
  lat(pointsm),                                                               &
          ! IN Latitude (degrees)
  long(pointsm),                                                              &
          ! IN Longitude (degrees)
  time_day,                                                                   &
          ! WORK Time of day (hours)
  timestep,                                                                   &
          ! WORK Timestep (seconds)
  random_num_sd
          ! WORK Random number associated with rai

REAL ::                                                                       &
  temp_conv,                                                                  &
          ! WORK Temperature above which rainfall is convective (K)
  temp_snow
          ! WORK Temperature below which snow occurs (K)
PARAMETER(temp_conv = 293.15)
PARAMETER(temp_snow = 275.15)

REAL ::                                                                       &
  init_hour_conv_rain,                                                        &
          !WORK Start of convective rain event (hour)
  init_hour_ls_rain,                                                          &
          !WORK Start of large scale rain event (hour)
  init_hour_ls_snow,                                                          &
          !WORK Start of large scale snow event (hour)
  dur_conv_rain_in_hours,                                                     &
          !WORK Start of convective rain event (hour)
  dur_ls_rain_in_hours,                                                       &
          !WORK Duration of large scale rain event (hour)
  dur_ls_snow_in_hours,                                                       &
          !WORK Start of large scale snow event (hour)
  end_hour_conv_rain,                                                         &
          !WORK End of convective rain event (hour)
  end_hour_ls_rain,                                                           &
          !WORK End of large scale rain event (hour)
  end_hour_ls_snow,                                                           &
          !WORK End of large scale snow event (hour)
  hourevent,                                                                  &
          !WORK Local variable giving hours during diurnal
          !     period for checking if precip. event occurs
  max_precip_rate,                                                            &
          !WORK Maximum allowed precip. rate allowed within
          !     each sub-daily timestep (mm/day). (This only
          !     applies when STEP_DAY >= 2)
  period_len
          !WORK Length of period (hr)


!-----------------------------------------------------------------------
! Calculate the maximum precipitation rate. It is noted that 58 mm/day
! over 8 timesteps, and where all fell within a single 3 hour period
! caused numerical issues for MOSES. This corresponded to a rate of
! 464 mm/day during the 3-hour period. Hence, place a limit of 350
! mm/day.
!-----------------------------------------------------------------------
PARAMETER(max_precip_rate = 350.0)

!-----------------------------------------------------------------------
! First check whether sub-daily calculations are required.
!-----------------------------------------------------------------------

!      DUR_CONV_RAIN_IN_HOURS = 2.0
!      DUR_LS_RAIN_IN_HOURS = 5.0
!      DUR_LS_SNOW_IN_HOURS = 5.0

!      DUR_CONV_RAIN_IN_HOURS = 6.0
!      DUR_LS_RAIN_IN_HOURS = 1.0
!      DUR_LS_SNOW_IN_HOURS = 1.0

dur_conv_rain_in_hours = dur_conv_rain / secs_in_hour
dur_ls_rain_in_hours   = dur_ls_rain   / secs_in_hour
dur_ls_snow_in_hours   = dur_ls_snow   / secs_in_hour

period_len = rhour_per_day / REAL(step_day)

!-----------------------------------------------------------------------
! First check whether sub-daily calculations are required.
!-----------------------------------------------------------------------
IF (step_day >= 2) THEN

  !-----------------------------------------------------------------------
  ! Ensure that the durations are at least as long as a time period for
  ! the model to prevent solution "falling through gaps"
  !-----------------------------------------------------------------------
  IF (dur_conv_rain_in_hours <= period_len)                                   &
    dur_conv_rain_in_hours = period_len + 1.0e-6
  IF (dur_ls_rain_in_hours <= period_len)                                     &
    dur_ls_rain_in_hours = period_len + 1.0e-6
  IF (dur_ls_snow_in_hours <= period_len)                                     &
    dur_ls_snow_in_hours = period_len + 1.0e-6

  timestep = REAL(sec_day) / REAL(step_day)

  !-----------------------------------------------------------------------
  ! Calculate the diurnal cycle in the SW radiation
  !-----------------------------------------------------------------------
  daynumber = INT((month-1.0) * REAL(day_mon)) + iday
  CALL sunny(                                                                 &
    daynumber,step_day,pointsm,1990,lat,long,sun,time_max                     &
  )

  !-----------------------------------------------------------------------
  ! Loop over timesteps
  !-----------------------------------------------------------------------
  DO istep = 1,step_day

    !-----------------------------------------------------------------------
    ! Calculate timestep values of the driving data.
    !-----------------------------------------------------------------------
    time_day = (REAL(istep) - 0.5) * timestep

    DO l = 1,pointsm
      t_sd(l,istep) = t_l(l) + 0.5 * dtemp_day_l(l) * COS(2.0 * pi            &
                      * (time_day - rsec_per_hour * time_max(l))              &
                      / REAL(sec_day))
      lw_sd(l,istep) = lw_l(l)                                                &
                     * (4.0 * t_sd(l,istep) / t_l(l) - 3.0)
      sw_sd(l,istep) = sw_l(l) * sun(l,istep)
    END DO

    !-----------------------------------------------------------------------
    ! Calculate timestep values of the driving data that is not split up
    ! into diurnal behaviour.
    !-----------------------------------------------------------------------
    DO l = 1,pointsm
      pstar_sd(l,istep) = pstar_l(l)
      wind_sd(l,istep)  = wind_l(l)
    END DO

    !-----------------------------------------------------------------------
    ! Check that humidity value is not greater than QSAT (but otherwise, QHU
    ! is not split up into diurnal behaviour).
    !-----------------------------------------------------------------------
    DO l = 1,pointsm
      pstar_sd_local(l) = pstar_sd(l,istep)
      t_sd_local(l)     = t_sd(l,istep)
    END DO

    CALL qsat(qs_sd_local,t_sd_local,pstar_sd_local,pointsm)

    DO l = 1,pointsm
      qhum_sd(l,istep) = 0.01 * rh15m_l(l) * qs_sd_local(l)
    END DO

  END DO   ! End of timestep loop within the individual days.

  !-----------------------------------------------------------------------
  ! Calculate daily rainfall disaggregation
  !-----------------------------------------------------------------------
  DO l = 1,pointsm

    !-----------------------------------------------------------------------
    ! Precipitation is split into four components,
    ! these being large scale rain, convective rain, large scale snow,
    ! convective snow. Call random number generator for different
    ! durations.
    !-----------------------------------------------------------------------
    CALL rndm(random_num_sd,seed_rain)

    !-----------------------------------------------------------------------
    ! Calculate type of precipitation. The decision is based purely up
    ! mean daily temperature, T_L. The cutoffs are:
    !
    ! Convective scale rain (duration CONV_RAIN_DUR): T_L > 20.0oC
    ! Large scale rain (duration LS_RAIN_DUR) : 20.0oC > T_L > 2oC
    ! Large scale snow (duration LS_SNOW_DUR) : T_L < 2oC
    ! Convective snow - IGNORED
    !
    !-----------------------------------------------------------------------
    ! Initialise arrays
    DO istep = 1,step_day
      conv_rain_sd(l,istep) = 0.0
      ls_rain_sd(l,istep)   = 0.0
      ls_snow_sd(l,istep)   = 0.0
      n_event(l,istep)      = 0
      n_event_local(istep)  = 0
      prec_loc(istep)       = 0.0
    END DO

    ! Calculate rainfall disaggregation.

    ! Start with convective rain
    ! (temperatures based upon mean daily temperature)

    ! First check if warm enough for convective rain
    IF (t_l(l) >= temp_conv) THEN

      init_hour_conv_rain = random_num_sd                                     &
                            * (rhour_per_day - dur_conv_rain_in_hours)
      end_hour_conv_rain = init_hour_conv_rain + dur_conv_rain_in_hours

      n_tally = 0

      DO istep = 1,step_day
        hourevent = (REAL(istep) - 0.5) * period_len
        IF (hourevent >= init_hour_conv_rain .AND.                            &
           hourevent < end_hour_conv_rain) THEN
          n_event(l,istep) = 1
          n_tally = n_tally + 1
        END IF
      END DO

      DO istep = 1,step_day
        IF (n_event(l,istep) == 1) THEN      !Rains on this day
          conv_rain_sd(l,istep) =                                             &
                    (REAL(step_day) / REAL(n_tally)) * precip_l(l)
          prec_loc(istep) = conv_rain_sd(l,istep)
          n_event_local(istep) = n_event(l,istep)
        END IF
      END DO

      ! Check that no convective rain periods
      ! exceed MAX_PRECIP_RATE, or if so,
      ! then redistribute. The variable that is redistributed is local
      ! variable PREC_LOC - CONV_RAIN_SD is then set to this after the
      ! call to REDIS.

      CALL redis(                                                             &
        nsdmax,step_day,max_precip_rate,prec_loc,                             &
        n_event_local,n_tally                                                 &
      )

      DO istep = 1,step_day
        conv_rain_sd(l,istep) = prec_loc(istep)
      END DO

      ! Now look at large scale rainfall components
    ELSE IF (t_l(l) < temp_conv .AND. t_l(l) >= temp_snow) THEN

      init_hour_ls_rain = random_num_sd                                       &
                          * (rhour_per_day - dur_ls_rain_in_hours)
      end_hour_ls_rain = init_hour_ls_rain + dur_ls_rain_in_hours

      n_tally = 0

      DO istep = 1,step_day
        hourevent = (REAL(istep) - 0.5) * period_len
        IF (hourevent >= init_hour_ls_rain .AND.                              &
           hourevent < end_hour_ls_rain) THEN
          n_event(l,istep) = 1
          n_tally = n_tally + 1
        END IF
      END DO

      DO istep = 1,step_day
        IF (n_event(l,istep) == 1) THEN      !Rains on this day
          ls_rain_sd(l,istep) =                                               &
                    (REAL(step_day) / REAL(n_tally)) * precip_l(l)
          prec_loc(istep) = ls_rain_sd(l,istep)
          n_event_local(istep) = n_event(l,istep)
        END IF
      END DO

      ! Check that no large scale rain periods
      ! exceed MAX_PRECIP_RATE, or if so,
      ! then redistribute.

      CALL redis(                                                             &
        nsdmax,step_day,max_precip_rate,prec_loc,                             &
        n_event_local,n_tally                                                 &
      )

      DO istep = 1,step_day
        ls_rain_sd(l,istep) = prec_loc(istep)
      END DO

      ! Now look at large scale snow components
    ELSE

      init_hour_ls_snow = random_num_sd                                       &
                          * (rhour_per_day - dur_ls_snow_in_hours)
      end_hour_ls_snow = init_hour_ls_snow + dur_ls_snow_in_hours

      n_tally = 0

      DO istep = 1,step_day
        hourevent = (REAL(istep) - 0.5) * period_len
        IF (hourevent >= init_hour_ls_snow .AND.                              &
           hourevent < end_hour_ls_snow) THEN
          n_event(l,istep) = 1
          n_tally = n_tally + 1
        END IF
      END DO

      DO istep = 1,step_day
        IF (n_event(l,istep) == 1) THEN       !Rains on this da
          ls_snow_sd(l,istep) =                                               &
                    (REAL(step_day) / REAL(n_tally)) * precip_l(l)
          prec_loc(istep) = ls_snow_sd(l,istep)
          n_event_local(istep) = n_event(l,istep)
        END IF
      END DO

      ! Check that no large scale snow periods exceed MAX_PRECIP_RATE,
      ! or if so, then redistribute.
      CALL redis(                                                             &
        nsdmax,step_day,max_precip_rate,prec_loc,                             &
        n_event_local,n_tally                                                 &
      )

      DO istep = 1,step_day
        ls_snow_sd(l,istep) = prec_loc(istep)
      END DO

    END IF

  END DO        ! End of large loop over different land points.
  !                    in calculation of different rainfall behaviours

ELSE           ! Now case where no subdaily variation (TSTEP=1)

  DO l = 1,pointsm
    pstar_sd_local(l) = pstar_sd(l,1)
    t_sd_local(l) = t_sd(l,1)
  END DO

  CALL qsat(qs_sd_local,t_sd_local,pstar_sd_local,pointsm)

  DO l = 1,pointsm
    sw_sd(l,1) = sw_l(l)
    t_sd(l,1) = t_l(l)
    lw_sd(l,1) = lw_l(l)
    pstar_sd(l,1) = pstar_l(l)
    wind_sd(l,1) = wind_l(l)
    qhum_sd(l,1) = 0.01 * rh15m_l(l) * qs_sd_local(l)

    IF (t_l(l) >= temp_conv) THEN
      conv_rain_sd(l,1) = precip_l(l)
    ELSE IF (t_l(l) < temp_conv .AND. t_l(l) >= temp_snow) THEN
      ls_rain_sd(l,1) = precip_l(l)
    ELSE
      ls_snow_sd(l,1) = precip_l(l)
    END IF

  END DO

END IF !End of loop to chose whether sub-daily is required

RETURN

END SUBROUTINE day_calc
