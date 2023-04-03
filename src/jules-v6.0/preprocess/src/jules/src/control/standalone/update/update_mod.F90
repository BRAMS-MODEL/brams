! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE update_mod

USE missing_data_mod, ONLY: imdi, rmdi

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module variables
!-----------------------------------------------------------------------------
LOGICAL ::                                                                    &
  l_imogen = .FALSE.
    ! Switch for using IMOGEN to provide driving data
    ! This is read in JULES DRIVE

INTEGER ::                                                                    &
  io_precip_type,                                                             &
    !  Flag indicating how precipitation is input
    !      1 = total precipitation is read in
    !      2 = values for total rainfall and total snowfall are read in
    !      3 = values for large-scale rainfall, convective rainfall and
    !          total snowfall are read in
    !      4 = values for convective rainfall, large-scale rainfall,
    !          convective snowfall and large-scale snowfall are read in
  io_rad_type,                                                                &
    !  Flag indicating how radiation is input
    !      1 downward longwave and downward shortwave fluxes provided
    !      2 net all wavelength flux (downward is positive) and downward
    !        shortwave flux are provided
    !      3 net downward longwave and net downward shortwave are provided
    !      4 downward longwave and net downward shortwave are provided
  precip_disagg_method = imdi
    !  Flag indicating how precip will be disaggregated
    !      1 no disaggregation
    !      2 imogen method
    !      3 as imogen method except no redistribution when precip is
    !        over max_precip_rate
    !      4 as above, but wet timesteps are distributed randomly through
    !        day, rather than in sequence

INTEGER, ALLOCATABLE :: prescribed_sthuf_levels(:)
    ! Stores the indices of the sm_levels which have prescribed sthuf.

LOGICAL ::                                                                    &
  io_wind_speed,                                                              &
    !   T means that the windspeed is input
    !   F means 2 components of wind are input
  use_diff_rad,                                                               &
    !   T means diffuse radiation is input
    !   F means diffuse radiation is set to a constant fraction
    !     (diff_frac_const)
  l_daily_disagg = .FALSE.,                                                   &
    !   T means expect daily driving data and disaggregate it
    !   F means diaggregationonly happens if l_imogen=T
  l_disagg_const_rh = .FALSE.,                                                &
    !   Switch controlling sub-daily disaggregation of humidity.
    !   T means keep relative humidity constant.
    !   F means keep specific humidity constant.
  l_perturb_driving = .FALSE.
    !   Switch indicating whether any perturbations should be applied to 
    !   temperature and precip driving data.

REAL ::                                                                       &
  bl_height       = 1000.0,                                                   &
    !   Height above surface of top of boundary layer (m).
    !   This value can be overridden via prescribed data.
  t_for_snow      = rmdi,                                                     &
    !   air temperature (K) at or below which precipitation is assumed to
    !   be snow
  t_for_con_rain  = rmdi,                                                     &
    !   air temperature (K) at or above which rainfall is taken to
    !   be convective (rather than large-scale) in nature
  diff_frac_const = rmdi,                                                     &
    !   a constant value for fraction of radiation that is diffuse
  dur_ls_rain     = rmdi,                                                     &
    !   duration of large scale rain event (seconds)
  dur_conv_rain   = rmdi,                                                     &
    !   duration of convective rain event (seconds)
  dur_ls_snow     = rmdi,                                                     &
    !   duration of large scale snow event (seconds)
  dur_conv_snow   = rmdi,                                                     &
    !   duration of convective snow event (seconds)
  temperature_abs_perturbation = rmdi,                                        &
    !   absolute perturbation to apply to temperature driving data                      
  precip_rel_perturbation = rmdi
    !   relative perturbation to apply to precip driving data

LOGICAL ::                                                                    &
  have_prescribed_veg = .FALSE.,                                              &
    !   T - a vegetation variable has been prescribed that requires a call
    !       to sparm after update
    !   F - no such variables have been prescribed
  have_prescribed_sthuf = .FALSE.
    !   T - sthuf is prescribed
    !   F - sthuf is given at initialisation

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UPDATE_MOD'

CONTAINS

! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE update_derived_variables(crop_vars,psparms,ainfo,urban_param,progs,&
                                    jules_vars)

!Use in relevant subroutines
USE sparm_mod,            ONLY: sparm

!Use in relevant variables
USE datetime_mod,         ONLY: secs_in_day
USE model_time_mod,       ONLY: current_time, timestep_len
USE jules_vegetation_mod, ONLY: l_croprotate, frac_min
USE jules_irrig_mod,      ONLY: l_irrig_dmd
USE ancil_info,           ONLY: land_pts, nsurft, surft_pts, nsoilt, soil_pts
USE theta_field_sizes,    ONLY: t_i_length
USE forcing,              ONLY: pstar_ij, u_1_ij, v_1_ij, u_0_ij, v_0_ij,     &
                                con_rain_ij, con_snow_ij, ls_rain_ij,         &
                                ls_snow_ij, sw_down_ij, lw_down_ij,           &
                                qw_1_ij, tl_1_ij
USE imogen_drive_vars,    ONLY: pstar_out, wind_out, conv_rain_out,           &
                                ls_rain_out, ls_snow_out, sw_out, lw_out,     &
                                qhum_out, t_out, conv_snow_out

USE string_utils_mod, ONLY: to_string
USE jules_surface_types_mod, ONLY: ncpft, nnpft
USE jules_soil_mod,       ONLY: dzsoil, sm_levels

USE water_constants_mod,  ONLY: rho_water

USE freeze_soil_mod,      ONLY: freeze_soil

USE logging_mod, ONLY: log_warn, log_fatal

!TYPE definitions
USE crop_vars_mod, ONLY: crop_vars_type
USE p_s_parms, ONLY: psparms_type
USE ancil_info,    ONLY: ainfo_type
USE urban_param_mod, ONLY: urban_param_type
USE prognostics, ONLY: progs_type
USE jules_vars_mod, ONLY: jules_vars_type

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Updates variables that are derived from those given in time-varying files
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
!Arguments
!TYPES containing field data (IN OUT)
TYPE(crop_vars_type), INTENT(IN OUT) :: crop_vars
TYPE(psparms_type), INTENT(IN OUT) :: psparms
TYPE(ainfo_type), INTENT(IN OUT) :: ainfo
TYPE(urban_param_type), INTENT(IN OUT) :: urban_param
TYPE(progs_type), INTENT(IN OUT) :: progs
TYPE(jules_vars_type), INTENT(IN OUT) :: jules_vars

! Work variables
INTEGER :: insd  ! Timestep in day - used to index IMOGEN arrays

INTEGER :: i,j,l,m,n  ! Index variables
LOGICAL :: reset_done  ! Indicates if a reset of frac to frac_min was
                         ! performed
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------


IF ( l_imogen ) THEN
  !-------------------------------------------------------------------------------
  ! If IMOGEN is enabled, copy the correct timestep of the current climatology
  ! into the driving variables
  !-------------------------------------------------------------------------------
  !Get the timestep in the day that we are on
  insd = (current_time%time / timestep_len) + 1

  DO l = 1,land_pts
    j = (ainfo%land_index(l) - 1) / t_i_length + 1
    i = ainfo%land_index(l) - (j-1) * t_i_length

    pstar_ij(i,j)    = pstar_out(l,current_time%month,current_time%day,insd)
    u_1_ij(i,j)      = wind_out(l,current_time%month,current_time%day,insd)
    v_1_ij(i,j)      = 0.0
    u_0_ij(i,j)      = 0.0
    v_0_ij(i,j)      = 0.0
    con_rain_ij(i,j) =                                                        &
      conv_rain_out(l,current_time%month,current_time%day,insd)               &
      / REAL(secs_in_day)
    con_snow_ij(i,j) =                                                        &
      conv_snow_out(l,current_time%month,current_time%day,insd)               &
      / REAL(secs_in_day)
    ls_rain_ij(i,j)  =                                                        &
      ls_rain_out(l,current_time%month,current_time%day,insd)                 &
      / REAL(secs_in_day)
    ls_snow_ij(i,j)  =                                                        &
      ls_snow_out(l,current_time%month,current_time%day,insd)                 &
      / REAL(secs_in_day)
    sw_down_ij(i,j)  = sw_out(l,current_time%month,current_time%day,insd)
    lw_down_ij(i,j)  = lw_out(l,current_time%month,current_time%day,insd)
    qw_1_ij(i,j)     = qhum_out(l,current_time%month,current_time%day,insd)
    tl_1_ij(i,j)     = t_out(l,current_time%month,current_time%day,insd)
  END DO

ELSE
  !-----------------------------------------------------------------------------
  ! Otherwise, update the main driving variables based on what was given in the
  ! input file(s) - the update of those variables has already taken place
  !-----------------------------------------------------------------------------
  ! Initialise wind based on io_wind_speed
  ! This is just a case of initialising variables that are not set via file
  ! u_1 is always set from file, v_1 is only set from file if using both
  ! wind components
  IF ( io_wind_speed ) THEN
    v_1_ij(:,:) = 0.0
  END IF
  u_0_ij(:,:) = 0.0
  v_0_ij(:,:) = 0.0

  !------------------------------------------------------------------------------
  ! Apply any perturbations
  IF ( l_perturb_driving ) THEN
    tl_1_ij(:,:) = tl_1_ij(:,:) + temperature_abs_perturbation
    con_rain_ij(:,:) = precip_rel_perturbation * con_rain_ij(:,:)
    con_snow_ij(:,:) = precip_rel_perturbation * con_snow_ij(:,:)
    ls_rain_ij(:,:)  = precip_rel_perturbation * ls_rain_ij(:,:)
    ls_snow_ij(:,:)  = precip_rel_perturbation * ls_snow_ij(:,:)
  END IF

  IF (l_daily_disagg) THEN
    !-----------------------------------------------------------------------------
    ! If we are using the disaggregator, we need to disaggregate the precip at
    ! the start of every day, then impose a diurnal cycle
    !-----------------------------------------------------------------------------
    IF ( current_time%time == 0 ) THEN
      CALL update_precip_variables()
      CALL fill_disaggregated_precip_arrays(progs)
    END IF

    CALL impose_diurnal_cycle()
  ELSE
    !-----------------------------------------------------------------------------
    ! Otherwise, just update the precip variables based on the io_precip_type
    !-----------------------------------------------------------------------------
    CALL update_precip_variables()
  END IF

  !-----------------------------------------------------------------------------
  ! Radiation variables are updated in CONTROL, since the new albedos need to be
  ! known before the update takes place
  !-----------------------------------------------------------------------------
END IF

!-------------------------------------------------------------------------------
! Copy information to U, V and T grids (assume that att grids are the same)
!-------------------------------------------------------------------------------
jules_vars%u_0_p_ij(:,:) = u_0_ij(:,:)
jules_vars%v_0_p_ij(:,:) = v_0_ij(:,:)
jules_vars%u_1_p_ij(:,:) = u_1_ij(:,:)
jules_vars%v_1_p_ij(:,:) = v_1_ij(:,:)


!-----------------------------------------------------------------------------
! Update any variables dependent on variables being prescribed that are not
! driving variables
!-----------------------------------------------------------------------------
IF ( have_prescribed_veg ) THEN
  WHERE ( progs%lai_pft(:,:) <= 0.01 )
    progs%lai_pft = 0.01
  END WHERE

  CALL sparm(land_pts, nsurft, surft_pts,                                     &
             ainfo%surft_index, ainfo%frac_surft,                             &
             progs%canht_pft, progs%lai_pft, psparms%z0m_soil_gb,             &
             psparms%catch_snow_surft, psparms%catch_surft, psparms%z0_surft, &
             psparms%z0h_bare_surft, urban_param%ztm_gb)

      !infiltration_rate does not need to be called because frac has not been
      !changed. This should be reviewed if a new parametrisation is added.

END IF

IF ( have_prescribed_sthuf ) THEN

  !-----------------------------------------------------------------------------
  ! Calculate soil moisture content from wetness
  !-----------------------------------------------------------------------------
  DO i = 1,sm_levels
    IF ( ANY(prescribed_sthuf_levels == i) ) THEN
      ! Use the sthuf_soilt that was just read in from the prescribed data file
      progs%smcl_soilt(:,:,i) = rho_water * dzsoil(i) *                       &
                          jules_vars%sthuf_soilt(:,:,i) *                     &
                          psparms%smvcst_soilt(:,:,i)
    ELSE
      ! Use the sthu_soilt and sthu_soilt values calculated on the previous
      ! timestep
      progs%smcl_soilt(:,:,i) = rho_water * dzsoil(i) *                       &
                    (psparms%sthu_soilt(:,:,i) + psparms%sthf_soilt(:,:,i)) * &
                    psparms%smvcst_soilt(:,:,i)
    END IF
  END DO

  !-----------------------------------------------------------------------------
  ! Calculate frozen and unfrozen fractions of soil moisture.
  ! freeze_soil assumes bexp_soilt, sathh_soilt, smvcst_soilt are
  ! constant in soil column, so loop around soil layers to allow
  ! depth varying soil properties and to maintain compatibility with the UM
  !-----------------------------------------------------------------------------
  DO i = 1,sm_levels
    DO m = 1,nsoilt
      CALL freeze_soil (land_pts,1,psparms%bexp_soilt(:,m,i),dzsoil(i:i),     &
                        psparms%sathh_soilt(:,m,i),progs%smcl_soilt(:,m,i),   &
                        progs%t_soil_soilt(:,m,i),psparms%smvcst_soilt(:,m,i),&
                        psparms%sthu_soilt(:,m,i),psparms%sthf_soilt(:,m,i))
    END DO
  END DO
END IF
!-----------------------------------------------------------------------------
!   If using l_croprotate ensure that fractions of
!   crop PFTs are not below minimum. Only do this over land points.
!   This ensures that each crop fraction is initialised properly so it can be
!   used later. If a fraction is zero at the start it cannot be used later.
!-----------------------------------------------------------------------------

IF (l_croprotate) THEN
  ! Set up a flag to see if any points were reset to frac_min
  reset_done = .FALSE.

  DO j = 1,soil_pts
    i = ainfo%soil_index(j)

    IF ( ANY( ainfo%frac_surft(i,nnpft+1:nnpft + ncpft) < frac_min ) ) THEN
      ! Reset all small values. Renormalisation is done later,
      ! but will fail if frac_min is sufficiently large.
      ! We only reset crop PFT tiles
      DO n = 1,ncpft
        WHERE ( ainfo%frac_surft(:,n + nnpft) < frac_min )
          ainfo%frac_surft(:,n + nnpft) = frac_min
        END WHERE
      END DO
      reset_done = .TRUE.
    END IF
  END DO

  IF ( reset_done )                                                           &
    CALL log_warn("update_derived_variables",                                 &
                  "frac < frac_min at one or more points - reset to " //      &
                  "frac_min at those points for crop tiles")

  !-----------------------------------------------------------------------------
  ! Check that ainfo%frac_surft sums to 1.0 (with a bit of leeway).
  !-----------------------------------------------------------------------------
  DO i = 1,land_pts
    IF ( ABS( SUM(ainfo%frac_surft(i,:)) - 1.0 ) >= 1.0e-2 )                  &
      ! If the discrepancy is big enough, bail
      CALL log_fatal("update_derived_variables",                              &
                     "frac does not sum to 1 at point " //                    &
                     TRIM(to_string(i)))
      ! Ignore any discrepancy below the threshold completely
  END DO
END IF
!-----------------------------------------------------------------------------
! Update irrigation fractions
!-----------------------------------------------------------------------------
IF ( l_irrig_dmd ) THEN
  CALL assign_irrig_fraction(crop_vars,ainfo)
  CALL update_irrig_variables(crop_vars,psparms,ainfo)
END IF

RETURN

END SUBROUTINE update_derived_variables


SUBROUTINE update_precip_variables()

USE forcing, ONLY: tl_1_ij, ls_snow_ij, ls_rain_ij, con_rain_ij, con_snow_ij
USE jules_surface_mod, ONLY: l_point_data

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Update precipitation variables
!   What we need to do depends on io_precip_type
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

SELECT CASE ( io_precip_type )
CASE ( 1 )
  ! Total precipitation is given in file and stored in ls_rain_ij

  ! First partition into rain and snow - all snow is assumed to be large-scale
  ls_snow_ij(:,:) = 0.0
  WHERE ( tl_1_ij(:,:) <= t_for_snow )
    ls_snow_ij = ls_rain_ij
    ls_rain_ij = 0.0
  END WHERE

  ! Now ls_rain_ij contains total rainfall
  ! If using point data assume all rain is large-scale
  con_rain_ij(:,:) = 0.0
  IF ( .NOT. l_point_data ) THEN
    ! Otherwise partition into convective and large scale based on t_for_con_rain
    WHERE ( tl_1_ij(:,:) >= t_for_con_rain )
      con_rain_ij = ls_rain_ij
      ls_rain_ij  = 0.0
    END WHERE
  END IF

  ! Convective snow is assumed to be 0 unless explicitly provided
  con_snow_ij(:,:) = 0.0

CASE ( 2 )
  ! Total rainfall given in file and currently stored in ls_rain_ij

  ! If using point data assume all rain is large-scale
  con_rain_ij(:,:) = 0.0
  IF ( .NOT. l_point_data ) THEN
    ! Otherwise partition into convective and large scale based on t_for_con_rain
    WHERE ( tl_1_ij(:,:) >= t_for_con_rain )
      con_rain_ij = ls_rain_ij
      ls_rain_ij  = 0.0
    END WHERE
  END IF

  ! Total snow was given in file, and all assumed to be large scale

  ! Convective snow is assumed to be 0 unless explicitly provided
  con_snow_ij(:,:) = 0.0

CASE ( 3 )
  ! Convective and large-scale rain given in file

  ! Total snow was given in file, and all assumed to be large-scale

  ! Convective snow is assumed to be 0 unless explicitly provided
  con_snow_ij(:,:) = 0.0

CASE ( 4 )
  ! There is nothing to be done for case 4, since all components are input from
  ! file
END SELECT

END SUBROUTINE update_precip_variables


SUBROUTINE fill_disaggregated_precip_arrays(progs)

USE disaggregated_precip, ONLY: ls_rain_disagg, con_rain_disagg,              &
                                 ls_snow_disagg, con_snow_disagg

USE forcing, ONLY:  ls_rain_ij, con_rain_ij, ls_snow_ij, con_snow_ij

USE model_time_mod, ONLY: timestep_len, timesteps_in_day

USE datetime_mod, ONLY: secs_in_day, hours_in_day, secs_in_hour

USE theta_field_sizes, ONLY: t_i_length, t_j_length

!TYPE definitions
USE prognostics, ONLY: progs_type

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!     Disaggregates a day mean precip into a precip event
!     Used only when l_daily_disagg = T
!
! Method:
!     Same as used in day_calc within IMOGEN
!
! Code Author: Karina Williams
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Arguments

!TYPES containing the data
TYPE(progs_type), INTENT(IN OUT) :: progs

INTEGER, PARAMETER :: n_prec = 4
CHARACTER(LEN=7) :: list_precip_str(n_prec)
CHARACTER(LEN=7) :: precip_str

INTEGER :: i, t_i, t_j
REAL :: dur_hours
REAL,PARAMETER :: small = 1.0e-6
REAL :: timestep_len_hours
REAL :: max_precip_rate
INTEGER :: n_tally

REAL,ALLOCATABLE ::    prec_3d(:,:,:)
REAL,ALLOCATABLE ::    prec_2d(:,:)
INTEGER,ALLOCATABLE :: n_event_local(:)
REAL,ALLOCATABLE ::    prec_1d(:)


!-----------------------------------------------------------------------------


IF ( .NOT. l_daily_disagg) RETURN

timestep_len_hours = REAL(timestep_len) / REAL(secs_in_hour)

!-----------------------------------------------------------------------------
! Calculate the maximum precipitation rate. It is noted that 58 mm/day
! over 8 timesteps, and where all fell within a single 3 hour period
! caused numerical issues for MOSES. This corresponded to a rate of
! 464 mm/day during the 3-hour period. Hence, place a limit of 350
! mm/day.
!-----------------------------------------------------------------------------
max_precip_rate = 350.0 / REAL(secs_in_day) ! in mm s-1 i.e. kg m-2 s-1

list_precip_str = (/ 'lscrain','conrain','lscsnow','consnow' /)

ALLOCATE(n_event_local(timesteps_in_day))
ALLOCATE(prec_1d(timesteps_in_day))
ALLOCATE(prec_2d(t_i_length, t_j_length))
ALLOCATE(prec_3d(t_i_length, t_j_length, timesteps_in_day))

DO i = 1,n_prec

  precip_str = list_precip_str(i)

  SELECT CASE( precip_str )
  CASE ( 'lscrain' )
    dur_hours = dur_ls_rain / REAL(secs_in_hour)
    prec_2d(:,:) = ls_rain_ij(:,:)

  CASE ( 'conrain' )
    dur_hours = dur_conv_rain / REAL(secs_in_hour)
    prec_2d(:,:) = con_rain_ij(:,:)

  CASE ( 'lscsnow' )
    dur_hours = dur_ls_snow / REAL(secs_in_hour)
    prec_2d(:,:) = ls_snow_ij(:,:)

  CASE ( 'consnow' )
    dur_hours = dur_conv_snow / REAL(secs_in_hour)
    prec_2d(:,:) = con_snow_ij(:,:)
  END SELECT

  !-----------------------------------------------------------------------------
  ! Ensure that the durations are at least as long as a time period for
  ! the model to prevent solution "falling through gaps"
  !-----------------------------------------------------------------------------

  IF ( dur_hours <= timestep_len_hours ) THEN
    dur_hours = timestep_len_hours + small
  END IF

  DO t_j = 1,t_j_length
    DO t_i = 1,t_i_length

      CALL fill_n_event_local(dur_hours, timestep_len_hours,                  &
                              precip_disagg_method, n_event_local, n_tally)

      ! recall prec_2d and prec_3d are rates
      prec_3d(t_i,t_j,:) =  REAL(n_event_local(:)) *                          &
                 (REAL(timesteps_in_day) / REAL(n_tally)) * prec_2d(t_i,t_j)

      IF ( precip_disagg_method == 2 ) THEN
        ! Check that no convective rain periods
        ! exceed max_precip_rate, or if so,
        ! then redistribute.

        prec_1d(:) = prec_3d(t_i,t_j,:)

        CALL redis(                                                           &
          timesteps_in_day, timesteps_in_day, max_precip_rate, prec_1d,       &
          n_event_local, n_tally                                              &
        )

        prec_3d(t_i,t_j,:) =  prec_1d(:)

      END IF

    END DO
  END DO

  SELECT CASE( precip_str )
  CASE ( 'lscrain' )
    ls_rain_disagg(:,:,:)  = prec_3d(:,:,:)

  CASE ( 'conrain' )
    con_rain_disagg(:,:,:) = prec_3d(:,:,:)

  CASE ( 'lscsnow' )
    ls_snow_disagg(:,:,:)  = prec_3d(:,:,:)

  CASE ( 'consnow' )
    con_snow_disagg(:,:,:) = prec_3d(:,:,:)
  END SELECT
END DO

IF ( precip_disagg_method > 1 ) THEN
  ! Update the stored seed in case a dump is written this timestep
  CALL RANDOM_SEED( get = progs%seed_rain )
END IF

DEALLOCATE(prec_3d)
DEALLOCATE(prec_2d)
DEALLOCATE(prec_1d)
DEALLOCATE(n_event_local)

END SUBROUTINE fill_disaggregated_precip_arrays


!###############################################################################
!###############################################################################


SUBROUTINE fill_n_event_local(dur_hours, timestep_len_hours,                  &
                              precip_disagg_method, n_event_local, n_tally)

USE model_time_mod, ONLY: timesteps_in_day
USE datetime_mod, ONLY: hours_in_day

IMPLICIT NONE

! Arguments
REAL, INTENT(IN) :: dur_hours
REAL, INTENT(IN) :: timestep_len_hours
INTEGER, INTENT(IN) :: precip_disagg_method

INTEGER, INTENT(OUT) :: n_event_local(:)
INTEGER, INTENT(OUT) :: n_tally

! Work variables
REAL :: init_hour_prec
REAL :: end_hour_prec
REAL :: hourevent
REAL :: random_num_sd

INTEGER, ALLOCATABLE :: loc(:)
INTEGER :: istep
INTEGER :: j


!-----------------------------------------------------------------------------

! Initialise arrays
n_event_local(:)  = 0

SELECT CASE( precip_disagg_method )
CASE (1)
  n_event_local(:)  = 1
  n_tally = SUM(n_event_local(:))

CASE (2, 3)
  CALL RANDOM_NUMBER(random_num_sd)

  init_hour_prec = random_num_sd * ( REAL(hours_in_day) - dur_hours )
  end_hour_prec = init_hour_prec + dur_hours

  DO istep = 1,timesteps_in_day
    hourevent = (REAL(istep) - 0.5) * timestep_len_hours

    IF ( hourevent >= init_hour_prec .AND. hourevent < end_hour_prec) THEN
      n_event_local(istep) = 1
    END IF
  END DO

  n_tally = SUM(n_event_local(:))

CASE (4)
  n_tally = NINT( dur_hours / REAL(hours_in_day) * timesteps_in_day )

  ALLOCATE(loc(n_tally))

  DO istep = 1, n_tally
    loc(istep) = istep
  END DO

  DO istep = n_tally+1, timesteps_in_day
    CALL RANDOM_NUMBER(random_num_sd) ! 0 <= random_num_sd < 1

    j = INT (random_num_sd * REAL(istep) ) + 1 ! 1 <= integer random number <= istep

    IF ( j <= n_tally ) THEN
      loc(j) = istep
    END IF
  END DO

  DO istep = 1, n_tally
    n_event_local(loc(istep)) = 1
  END DO

  DEALLOCATE(loc)

END SELECT

END SUBROUTINE fill_n_event_local


!###############################################################################
!###############################################################################

SUBROUTINE impose_diurnal_cycle()

USE model_time_mod, ONLY: current_time, timestep_len,                         &
                           timesteps_in_day

USE forcing, ONLY:  tl_1_ij, lw_down_ij, sw_down_ij, qw_1_ij, pstar_ij,       &
                     ls_rain_ij, con_rain_ij, ls_snow_ij, con_snow_ij,        &
                     diurnal_temperature_range_ij, diff_rad_ij

USE model_grid_mod, ONLY: latitude, longitude

USE datetime_mod, ONLY: secs_in_day, secs_in_hour, l_360, l_leap

USE datetime_utils_mod, ONLY: day_of_year, days_in_year


USE theta_field_sizes, ONLY: t_i_length, t_j_length

USE conversions_mod, ONLY: pi

USE disaggregated_precip, ONLY: ls_rain_disagg, con_rain_disagg,              &
                                 ls_snow_disagg, con_snow_disagg

USE qsat_mod, ONLY: qsat

IMPLICIT NONE

REAL, ALLOCATABLE :: tl_1_nodiurnalcycle(:,:)

REAL :: qsat_1d(t_i_length * t_j_length)
REAL :: qsat_1d_day(t_i_length * t_j_length)
!      Saturated specific humidity given daily mean temperature (kg kg-1).
REAL :: qsat_2d(t_i_length,t_j_length)
REAL, ALLOCATABLE :: sun(:,:)
REAL :: rh_2d_day(t_i_length,t_j_length)
!     Daily mean relative humidity (1).
REAL :: sun_tstep(t_i_length,t_j_length)
REAL :: time_max_1d(t_i_length * t_j_length)
REAL :: time_max_2d(t_i_length,t_j_length)

INTEGER :: n_points
INTEGER :: daynumber
INTEGER :: x_tofday

!-----------------------------------------------------------------------------
! Description:
!     Imposes a diurnal cycle on top of the forcing variables tl_1_ij,
!     lw_down_ij and sw_down_ij
!     Copy ls_rain_ij, con_rain_ij, ls_snow_ij, con_snow_ij from
!     disaggregated_precip variables
!     Makes sure qw_1_ij isn't above qsat (if l_disagg_const_rh=F), or imposes
!     constant relative humidity (if l_disagg_const_rh=T)
!
!      Used only when l_daily_disagg = T
!
! Method:
!      Same as used in day_calc within IMOGEN
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IF ( .NOT. l_daily_disagg ) RETURN

ALLOCATE( tl_1_nodiurnalcycle( t_i_length, t_j_length ))

n_points = t_i_length * t_j_length

tl_1_nodiurnalcycle(:,:) = tl_1_ij(:,:)

x_tofday = current_time%time / timestep_len + 1 ! x_tofday = 1 at midnight

ALLOCATE(sun(n_points,timesteps_in_day))

!sunny calls solpos, which assumes 360 day years
daynumber  = day_of_year(                                                     &
  current_time%year, current_time%month, current_time%day, l_360, l_leap      &
)
daynumber  = NINT( REAL(daynumber) * 360.0                                    &
                / REAL(days_in_year(current_time%year, l_360, l_leap)))

CALL sunny( daynumber, timesteps_in_day, n_points, current_time%year,         &
            RESHAPE(latitude(:,:),  (/ n_points /)),                          &
            RESHAPE(longitude(:,:), (/ n_points /)),                          &
            sun, time_max_1d )

time_max_2d(:,:)   = RESHAPE(time_max_1d(:),                                  &
                     (/ t_i_length, t_j_length /))

tl_1_ij(:,:) = tl_1_nodiurnalcycle(:,:) +                                     &
               0.5 * diurnal_temperature_range_ij(:,:) *                      &
               COS( 2.0 * pi * ( current_time%time - secs_in_hour *           &
                    time_max_2d(:,:) )                                        &
               / secs_in_day )

IF ( l_disagg_const_rh ) THEN
  !     Calculate relative humidity given daily mean temperature.

  CALL qsat(qsat_1d_day,                                                      &
          RESHAPE(tl_1_nodiurnalcycle(:,:), (/ n_points /)),                  &
          RESHAPE(pstar_ij(:,:),(/ n_points /)), n_points)

  rh_2d_day(:,:) = MIN( qw_1_ij(:,:) /                                        &
                 RESHAPE( qsat_1d_day(:), (/ t_i_length, t_j_length /) ),     &
                 1.0 )
END IF

!   Calulate saturated humidity given current temperature.

CALL qsat(qsat_1d,                                                            &
          RESHAPE(tl_1_ij(:,:), (/ n_points /)),                              &
          RESHAPE(pstar_ij(:,:),(/ n_points /)), n_points)

sun_tstep(:,:) = RESHAPE(sun(:, x_tofday),                                    &
                   (/ t_i_length, t_j_length /))
qsat_2d(:,:)   = RESHAPE(qsat_1d(:),                                          &
                   (/ t_i_length, t_j_length /))

lw_down_ij(:,:) = lw_down_ij(:,:) *                                           &
                  (4.0 * tl_1_ij(:,:) / tl_1_nodiurnalcycle(:,:) - 3.0)

sw_down_ij(:,:) = sw_down_ij(:,:) * sun_tstep

IF ( use_diff_rad ) THEN
  diff_rad_ij(:,:) = diff_rad_ij(:,:) * sun_tstep
END IF

IF ( l_disagg_const_rh ) THEN
  !     Constant relative humidity.
  qw_1_ij(:,:) = rh_2d_day(:,:) * qsat_2d(:,:)
ELSE
  !     Constant specific humidity.
  qw_1_ij(:,:) = MIN( qw_1_ij(:,:), qsat_2d(:,:) )
END IF

ls_rain_ij( :,:)  =  ls_rain_disagg( :,:, x_tofday)
con_rain_ij(:,:)  =  con_rain_disagg(:,:, x_tofday)
ls_snow_ij( :,:)  =  ls_snow_disagg( :,:, x_tofday)
con_snow_ij(:,:)  =  con_snow_disagg(:,:, x_tofday)

DEALLOCATE(tl_1_nodiurnalcycle)
DEALLOCATE(sun)

END SUBROUTINE impose_diurnal_cycle
!******************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT****************************************

SUBROUTINE calc_downward_rad(sea_ice_albedo, photosynth_act_rad,ainfo,progs,  &
                             coast,jules_vars)

!Module Imports

USE theta_field_sizes, ONLY: t_i_length, t_j_length

USE ancil_info, ONLY:                                                         &
  nsurft, land_pts

USE forcing, ONLY:                                                            &
  lw_down_ij, sw_down_ij, diff_rad_ij

USE fluxes, ONLY:                                                             &
  land_albedo_ij, alb_surft, sw_surft

USE csigma, ONLY:                                                             &
  sbcon

USE jules_sea_seaice_mod, ONLY:                                               &
  nice_use

USE jules_radiation_mod, ONLY: wght_alb

USE atm_fields_bounds_mod, ONLY: tdims

!TYPE definitions
USE ancil_info, ONLY: ainfo_type
USE prognostics, ONLY: progs_type
USE coastal,       ONLY: coastal_type
USE jules_vars_mod, ONLY: jules_vars_type

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Contains the calculations for downward radiation for standalone JULES
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

!Arguments
REAL, INTENT(IN) :: sea_ice_albedo(t_i_length,t_j_length,4)  ! Sea ice albedo

REAL, INTENT(OUT)::                                                           &
photosynth_act_rad(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!                                     ! Net downward shortwave radiation
!                                     !  in band 1 (w/m2).

!TYPES containing the data
TYPE(ainfo_type), INTENT(IN OUT) :: ainfo
TYPE(progs_type), INTENT(IN OUT) :: progs
TYPE(coastal_type), INTENT(IN OUT) :: coast
TYPE(jules_vars_type), INTENT(IN OUT) :: jules_vars

!Local Variables
INTEGER :: i,j,n,k,l

!End of header

!-------------------------------------------------------------------------------
!   Change radiation to be downward components if not using io_rad_type=1
!   NOTE this assumes that the point is either 100% land or 100% sea-ice.
!   1 downward fluxes provided
!   2 net all wavelength flux (downward is positive and is stored in downward
!     longwave variable) and downward shortwave fluxes are provided
!   3 net downward fluxes are provided (in downward variables)
!   4 downward longwave and net downward shortwave are provided
!   One day we would probably like to do this in driveUpdate or the likes, but
!   at present we don't
!   have all the required variables/masks there (I think).
!-------------------------------------------------------------------------------

!Convert net total radiation to net longwave. Net total is currently stored
!in lw_down_ij. Use the average of diffuse albedos in visible and NIR on land.
IF ( io_rad_type == 2 ) THEN
  DO j = 1,t_j_length
    DO i = 1,t_i_length
      IF ( ainfo%land_mask(i,j) ) THEN
        lw_down_ij(i,j) = lw_down_ij(i,j) - sw_down_ij(i,j) *                 &
                          (1.0 -                                              &
                           (wght_alb(1) * land_albedo_ij(i,j,1) +             &
                            wght_alb(2) * land_albedo_ij(i,j,2) +             &
                            wght_alb(3) * land_albedo_ij(i,j,3) +             &
                            wght_alb(4) * land_albedo_ij(i,j,4) ))

      ELSE
        !Assume not using CICE scheme so albedo same for all rad bands
        lw_down_ij(i,j) = lw_down_ij(i,j) -                                   &
          sw_down_ij(i,j) * (1.0 - sea_ice_albedo(i,j,1))
      END IF
    END DO
  END DO
END IF   !  io_rad_type

!Convert shortwave from net to downward.
!Net flux is currently stored in sw_down_ij.
!Use the average of diffuse albedos in visible and NIR on land.
IF ( io_rad_type == 3 .OR. io_rad_type == 4 ) THEN
  DO j = 1,t_j_length
    DO i = 1,t_i_length
      IF ( ainfo%land_mask(i,j) ) THEN
        sw_down_ij(i,j) = sw_down_ij(i,j) /                                   &
                          (1.0 -                                              &
                           (wght_alb(1) * land_albedo_ij(i,j,1) +             &
                            wght_alb(2) * land_albedo_ij(i,j,2)  +            &
                            wght_alb(3) * land_albedo_ij(i,j,3)  +            &
                            wght_alb(4) * land_albedo_ij(i,j,4) ))

      ELSE
        !Assume not using CICE scheme so albedo same for all rad bands
        sw_down_ij(i,j) = sw_down_ij(i,j) / ( 1.0 - sea_ice_albedo(i,j,1) )
      END IF
    END DO
  END DO
END IF   !  io_rad_type

!Convert longwave from net to downward. Net longwave is currently stored in
!lw_down_ij.
IF ( io_rad_type == 2 .OR. io_rad_type == 3 ) THEN
  DO j = 1,t_j_length
    DO i = 1,t_i_length
      IF ( .NOT. ainfo%land_mask(i,j) ) THEN
        DO n = 1,nice_use
          lw_down_ij(i,j) = lw_down_ij(i,j) +                                 &
                            ainfo%ice_fract_ncat_sicat(i,j,n) * sbcon         &
                            * coast%tstar_sice_sicat(i,j,n)**4.0
        END DO
      END IF
    END DO
  END DO
  DO n = 1,nsurft
    DO l = 1,land_pts
      j = ( ainfo%land_index(l) - 1 ) / t_i_length + 1
      i = ainfo%land_index(l) - ( j - 1 ) * t_i_length
      lw_down_ij(i,j) = lw_down_ij(i,j) + ainfo%frac_surft(l,n) * sbcon *     &
                        progs%tstar_surft(l,n)**4.0
    END DO
  END DO
END IF   !  io_rad_type


! Now we know sw_down_ij, we can update diff_frac if required
IF ( use_diff_rad ) THEN
  k = 0
  DO j = 1,t_j_length
    DO i = 1,t_i_length
      k = k + 1
      IF ( sw_down_ij(i,j) > 1.0 ) THEN
        jules_vars%diff_frac(k) = diff_rad_ij(i,j) / sw_down_ij(i,j)
        jules_vars%diff_frac(k) = MIN( 1.0, jules_vars%diff_frac(k) )
      ELSE
        jules_vars%diff_frac(k) = 0.0
      END IF
    END DO
  END DO
ELSE
  jules_vars%diff_frac(:) = diff_frac_const
END IF

!Calculate net SW radiation on tiles.
!Use the average of diffuse albedos in visible and NIR.
DO n = 1,nsurft
  DO l = 1,land_pts
    j = ( ainfo%land_index(l) - 1 ) / t_i_length + 1
    i = ainfo%land_index(l) - (j-1) * t_i_length
    sw_surft(l,n) = ( 1.0 - ( wght_alb(1) * alb_surft(l,n,1) +                &
                             wght_alb(2) * alb_surft(l,n,2)  +                &
                             wght_alb(3) * alb_surft(l,n,3)  +                &
                             wght_alb(4) * alb_surft(l,n,4)                   &
                   ) ) * sw_down_ij(i,j)

  END DO
END DO

!Calculate photosynthetically active radiation (PAR).
photosynth_act_rad(:,:) = 0.5 * sw_down_ij(:,:)

RETURN
END SUBROUTINE calc_downward_rad
! *****************************COPYRIGHT***************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use
! and distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT***************************************

SUBROUTINE assign_irrig_fraction (crop_vars,ainfo)

USE logging_mod, ONLY: log_fatal

USE string_utils_mod, ONLY: to_string

USE ancil_info, ONLY:                                                         &
land_pts, nsurft, nsoilt

USE jules_surface_types_mod, ONLY: ntype

USE jules_irrig_mod, ONLY:                                                    &
!  imported logicals with intent(in)
   set_irrfrac_on_irrtiles,                                                   &
   frac_irrig_all_tiles, irrtiles, nirrtile


!TYPE definitions
USE crop_vars_mod, ONLY: crop_vars_type
USE ancil_info,    ONLY: ainfo_type

IMPLICIT NONE

!Arguments
!TYPES containing field data (IN OUT)
TYPE(crop_vars_type), INTENT(IN OUT) :: crop_vars
TYPE(ainfo_type), INTENT(IN OUT) :: ainfo

!Local Variables
INTEGER ::                                                                    &
 m,l,k,n
                       ! Loop indices

INTEGER ::   alltiles(nsurft)

! Local arrays
REAL ::                                                                       &
 frac_irr_rest_soilt(land_pts,nsoilt)
                      ! Remaining irrigation fraction to assign to tiles

!-------------------------------------------------------------------------------
!   Update irrigation ancillary
!-------------------------------------------------------------------------------
! create index of all tiles (veg+non-veg)
DO m = 1,nsurft
  alltiles(m) = m
END DO

! Loop over land points and assign irrigation fraction
! note frac_irr_all has shape (land_points, year) or (land_points, 1)

IF ( set_irrfrac_on_irrtiles ) THEN

  crop_vars%frac_irr_old_soilt(:,:) = crop_vars%frac_irr_soilt(:,:)
  crop_vars%frac_irr_soilt(:,:) = 0.0

  DO k = 1,nirrtile ! loop over irrigated pfts
    n = irrtiles(k)
    !Set the current soil tile (see notice above)
    IF (nsoilt == 1) THEN
      !There is only 1 soil tile
      m = 1
    ELSE ! nsoilt == nsurft
      !Soil tiles map directly on to surface tiles
      m = k
    END IF

    DO l = 1,land_pts
      crop_vars%frac_irr_soilt(l,m) = crop_vars%frac_irr_soilt(l,m) +         &
                             crop_vars%irrfrac_irrtiles(l,1) *                &
                             ainfo%frac_surft(l, n)
    END DO
  END DO
ELSE
  DO m = 1,nsoilt
    DO l = 1,land_pts
      crop_vars%frac_irr_old_soilt(l,m) = crop_vars%frac_irr_soilt(l,m)
      crop_vars%frac_irr_soilt(l,m) = crop_vars%frac_irr_all(l,1)
    END DO
  END DO
END IF

WHERE (crop_vars%frac_irr_soilt < 1.0e-3) crop_vars%frac_irr_soilt = 0.0

!-------------------------------------------------------------------------------
!   Assign irrigation fraction to each tile
!   hadrd - this was originally done in physiol.f90
!-------------------------------------------------------------------------------
IF ( .NOT. frac_irrig_all_tiles .AND. .NOT. set_irrfrac_on_irrtiles ) THEN
  ! assign irrigated fraction to pre-defined tiles ONLY
  ! in decreasing order of importance
  crop_vars%frac_irr_surft(:,:) = 0.0
  frac_irr_rest_soilt(:,:) = crop_vars%frac_irr_soilt(:,:)

  DO l = 1,land_pts
    DO k = 1,nirrtile ! loop over irrigated pfts

      !==============================================================================
      !**NOTICE REGARDING SOIL TILING**
      !
      !The following section facilitates the use of soil tiling. As implemented,
      !there are two soil tiling options:
      !
      !nsoilt == 1
      !Operate as with a single soil tile, functionally identical to JULES upto
      ! at least vn4.7 (Oct 2016)
      ! This means that a soilt variable being passed 'up' to the surface is
      ! broadcast to the surft variable (with weighting by frac if requred)
      !
      !nsoilt > 1
      !Operate with nsoilt = nsurft, with a direct mapping between them
      ! This means that a soilt variable being passed 'up' to the surface is simply
      ! copied into the surft variable
      !
      ! This will need to be refactored for other tiling approaches. This note
      ! will be replicated elsewhere in the code as required
      !
      !These comments apply until **END NOTICE REGARDING SOIL TILING**
      !==============================================================================

              !Set the current soil tile (see notice above)
      IF (nsoilt == 1) THEN
        !There is only 1 soil tile
        m = 1
      ELSE ! nsoilt == nsurft
        !Soil tiles map directly on to surface tiles
        m = k
      END IF !nsoilt

      !==============================================================================
      !**END NOTICE REGARDING SOIL TILING**
      !==============================================================================

      ! for each tile, check if the index corresponds with index in irrtiles
      DO n = 1,ntype
        IF ( alltiles(n) == irrtiles(k) ) THEN
          ! assign (remaining) irrigated fraction to this tile
          crop_vars%frac_irr_surft(l,n)      = MIN(frac_irr_rest_soilt(l,m),  &
                                         ainfo%frac_surft(l,n) )
          frac_irr_rest_soilt(l,m) = frac_irr_rest_soilt(l,m)                 &
                                     - crop_vars%frac_irr_surft(l,n)
          ! check for negative remaining frac_irr_soilt
          IF (frac_irr_rest_soilt(l,m) < 0.0) THEN
            CALL log_fatal("update_ancil_irrig",                              &
              "Error in assigning irrigated fraction at point " //            &
              to_string(l) // " - " //                                        &
              "irrigated fraction: " //                                       &
              to_string(crop_vars%frac_irr_soilt(l,m)) // ", " //             &
              "tile fraction: " // to_string(ainfo%frac_surft(l,n)) // ", " //&
              "irrigated fraction assigned:"//                                &
              to_string(crop_vars%frac_irr_surft(l,n)))
          END IF
        END IF
      END DO ! ntype
    END DO ! irrtilenames
  END DO ! land_pts

  !to prevent negative fractions (this can happen due to numerical inaccuracy)
  crop_vars%frac_irr_surft(:,:) = MAX(crop_vars%frac_irr_surft(:,:),0.0)

  ! in sf_evap, frac_irr_surft is used as a multiplier to frac
  WHERE ( ainfo%frac_surft > 0.0 )
    crop_vars%frac_irr_surft = crop_vars%frac_irr_surft / ainfo%frac_surft
  END WHERE

ELSE ! if frac_irrig_all_tiles = TRUE
  ! assign irrigation fraction to all tiles
  ! to reproduce the original results
  crop_vars%frac_irr_surft(:,:) = 0.0
  DO n = 1,nsurft

    !==============================================================================
    !**NOTICE REGARDING SOIL TILING**
    !
    !The following section facilitates the use of soil tiling. As implemented,
    !there are two soil tiling options:
    !
    !nsoilt == 1
    !Operate as with a single soil tile, functionally identical to JULES upto
    ! at least vn4.7 (Oct 2016)
    ! This means that a soilt variable being passed 'up' to the surface is
    ! broadcast to the surft variable (with weighting by frac if requred)
    !
    !nsoilt > 1
    !Operate with nsoilt = nsurft, with a direct mapping between them
    ! This means that a soilt variable being passed 'up' to the surface is simply
    ! copied into the surft variable
    !
    ! This will need to be refactored for other tiling approaches. This note
    ! will be replicated elsewhere in the code as required
    !
    !These comments apply until **END NOTICE REGARDING SOIL TILING**
    !==============================================================================

          !Set the current soil tile (see notice above)
    IF (nsoilt == 1) THEN
      !There is only 1 soil tile
      m = 1
    ELSE ! nsoilt == nsurft
      !Soil tiles map directly on to surface tiles
      m = n
    END IF !nsoilt

    !==============================================================================
    !**END NOTICE REGARDING SOIL TILING**
    !==============================================================================

    WHERE ( crop_vars%frac_irr_soilt(:,m) > EPSILON(1.0) )
      crop_vars%frac_irr_surft(:,n) = crop_vars%frac_irr_soilt(:,m)
    END WHERE
  END DO

END IF ! .not. frac_irrig_all_tiles

END SUBROUTINE assign_irrig_fraction
! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use
! and distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT****************************************

SUBROUTINE update_irrig_variables (crop_vars,psparms,ainfo)

USE logging_mod, ONLY: log_info

USE model_time_mod, ONLY: current_time

USE datetime_mod, ONLY: datetime_to_string

USE ancil_info, ONLY: land_pts, nsoilt

USE jules_surface_types_mod, ONLY: ncpft, nnpft

USE jules_irrig_mod, ONLY: irr_crop, irr_crop_dvimax

USE jules_soil_mod, ONLY: sm_levels

!TYPE definitions
USE crop_vars_mod, ONLY: crop_vars_type
USE p_s_parms, ONLY: psparms_type
USE ancil_info,    ONLY: ainfo_type

IMPLICIT NONE

!Arguments
!TYPES containing field data (IN OUT)
TYPE(crop_vars_type), INTENT(IN OUT) :: crop_vars
TYPE(psparms_type), INTENT(IN OUT) :: psparms
TYPE(ainfo_type), INTENT(IN OUT) :: ainfo

INTEGER ::                                                                    &
 l,n,m

CHARACTER(LEN=50) ::                                                          &
 frac_irr_str

! Local arrays
REAL ::                                                                       &
 sthu_nir(land_pts,sm_levels)
                      ! soil moisture content in non-irrigated fraction as
                      ! a fraction of saturation

!-------------------------------------------------------------------------------
!   When using JULES-crop model, adjust irrigation fraction
!   to those tiles that have suitable development index
!-------------------------------------------------------------------------------
!   When irr_crop=2, JULES-crop is used to determine irrigation period.
!   However, within a grid box the development index can be different for each
!   crop tile while irrigation is applied to the (single) irrigated fraction of
!   that grid box.
!   Here frac_irr_surft is set to zero in those tiles where the dvi_cpft is below
!   the threshold and irrigation is only applied to those tiles with suitable
!   dvi_cpft, i.e. the actual irrigation fraction grows or shrinks according to
!   which crop tiles have the right dvi_cpft. Note that frac_irr_soilt can also be
!   allocated to non-crop tiles, these will (currently) only be irrigated when
!   any crop has a dvi_cpft above the threshold - see subroutine irrig_dmd.

IF ( irr_crop == irr_crop_dvimax ) THEN
  DO l = 1,land_pts
    DO n = 1,ncpft

      !==============================================================================
      !**NOTICE REGARDING SOIL TILING**
      !
      !The following section facilitates the use of soil tiling. As implemented,
      !there are two soil tiling options:
      !
      !nsoilt == 1
      !Operate as with a single soil tile, functionally identical to JULES upto
      ! at least vn4.7 (Oct 2016)
      ! This means that a soilt variable being passed 'up' to the surface is
      ! broadcast to the surft variable (with weighting by frac if requred)
      !
      !nsoilt > 1
      !Operate with nsoilt = nsurft, with a direct mapping between them
      ! This means that a soilt variable being passed 'up' to the surface is simply
      ! copied into the surft variable
      !
      ! This will need to be refactored for other tiling approaches. This note
      ! will be replicated elsewhere in the code as required
      !
      !These comments apply until **END NOTICE REGARDING SOIL TILING**
      !==============================================================================
              !Set the current soil tile (see notice above)
      IF (nsoilt == 1) THEN
        !There is only 1 soil tile
        m = 1
      ELSE ! nsoilt == nsurft
        !Soil tiles map directly on to surface tiles, remembering we're
        !working on crop tiles
        m = nnpft + n
      END IF !nsoilt
      !==============================================================================
      !**END NOTICE REGARDING SOIL TILING**
      !==============================================================================

      IF ( crop_vars%dvi_cpft(l,n) <= -1 ) THEN
        !         no suitable dvi_cpft - set frac_irr_soilt for this crop tile to zero
        !             and subtract from overall irrigation fraction
        crop_vars%frac_irr_soilt(l,m) = crop_vars%frac_irr_soilt(l,m)         &
                    - ainfo%frac_surft(l,n + nnpft)                           &
                    * crop_vars%frac_irr_surft(l,n + nnpft)
        crop_vars%frac_irr_soilt(l,m) =                                       &
                    MAX( crop_vars%frac_irr_soilt(l,m), 0.0 )
        crop_vars%frac_irr_surft(l,n + nnpft) = 0.0
      END IF
    END DO
  END DO
END IF

!-------------------------------------------------------------------------------
!   Update soil moisture content in irrigated fraction (sthu_irr_soilt)
!-------------------------------------------------------------------------------
!   If frac_irr expands, the 'added fraction' does not have the moisture content
!   of the irrigated fraction, but of the non-irrigated fraction (which is not
!   a prognostic). Total gridbox moisture content should remain the same
!   Conversely, when frac_irr shrinks, the moisture content in the non-irrig
!   fraction should become higher, but in the (remaining) irrigated fraction
!   (and total gridbox) it remains the same, so no need to cater for that here

DO m = 1, nsoilt
  DO l = 1,land_pts
    IF ( crop_vars%frac_irr_soilt(l,m) > crop_vars%frac_irr_old_soilt(l,m)) THEN
      DO n = 1,sm_levels

        !Note the hard-coded index of 1 for soil tiling.
        sthu_nir(l,n) = (psparms%sthu_soilt(l,m,n)                            &
                         - crop_vars%frac_irr_old_soilt(l,m)                  &
                         * crop_vars%sthu_irr_soilt(l,m,n))                   &
                        / (1.0 - crop_vars%frac_irr_old_soilt(l,m))

        crop_vars%sthu_irr_soilt(l,m,n) =                                     &
          (crop_vars%frac_irr_old_soilt(l,m) *                                &
           crop_vars%sthu_irr_soilt(l,m,n) +                                  &
          (crop_vars%frac_irr_soilt(l,m) - crop_vars%frac_irr_old_soilt(l,m)) &
          * sthu_nir(l,n)) / crop_vars%frac_irr_soilt(l,m)
      END DO
    END IF
  END DO
END DO

!------------------------------------------------------------------------------
! Write to stdout irrigation fraction for each gridbox.
!-------------------------------------------------------------------------------
WRITE(frac_irr_str, "(' Range of frac irrigation =',8f5.2,' to ',8f5.2)")     &
                                          MINVAL(crop_vars%frac_irr_soilt),   &
                                          MAXVAL(crop_vars%frac_irr_soilt)

IF (current_time%day == 1) THEN
  CALL log_info("update_ancil_irrig",                                         &
                '### NB The ranges below include any ice points. ###')
  CALL log_info("update_ancil_irrig",                                         &
                datetime_to_string(current_time) // frac_irr_str)
END IF


END SUBROUTINE update_irrig_variables

END MODULE update_mod
