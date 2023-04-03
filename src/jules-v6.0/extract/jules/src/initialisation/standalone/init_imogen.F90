#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_imogen(nml_dir,progs_data, trifctltype)

USE io_constants, ONLY: max_file_name_len, imogen_unit, namelist_unit

USE datetime_mod, ONLY: secs_in_day, l_360

USE string_utils_mod, ONLY: to_string

USE model_time_mod, ONLY: main_run_start, is_spinup, timestep_len

USE dump_mod, ONLY: read_dump

USE update_mod, ONLY: l_imogen, diff_frac_const

USE ancil_info, ONLY: land_pts

USE imogen_constants, ONLY: drive_month

USE imogen_map, ONLY: n_imogen_land, sgindinv, get_imogen_map

USE imogen_time, ONLY: step_day, nsdmax, mm, md

USE imogen_run, ONLY: imogen_run_list, file_points_order, wgen, c_emissions,  &
                      include_co2, include_non_co2, land_feed_co2,            &
                      land_feed_ch4, ocean_feed, anlg, anom, nyr_emiss,       &
                      file_scen_emits, ch4_init_ppbv, co2_init_ppmv,          &
                      initialise_from_dump, dump_file

USE imogen_anlg_vals, ONLY: imogen_anlg_vals_list, dir_clim,                  &
                            diff_frac_const_imogen

USE imogen_clim, ONLY: t_clim, rainfall_clim, snowfall_clim, rh15m_clim,      &
                       uwind_clim, vwind_clim, dtemp_clim, pstar_ha_clim,     &
                     sw_clim, lw_clim, f_wet_clim, lat, dctot_co2, dctot_ch4, &
                       longmin_clim, latmin_clim, longmax_clim, latmax_clim,  &
                       long

USE imogen_drive_vars, ONLY: t_out, conv_rain_out, conv_snow_out, ls_rain_out,&
                             ls_snow_out, qhum_out, wind_out, pstar_out,      &
                             sw_out, lw_out

USE imogen_io_vars, ONLY: nyr_max, yr_emiss, c_emiss

USE imogen_progs, ONLY: ch4_ppbv, co2_ppmv, co2_change_ppmv, dtemp_o,         &
                        fa_ocean, seed_wg

USE aero, ONLY: co2_mmr

USE logging_mod, ONLY: log_info, log_fatal

USE errormessagelength_mod, ONLY: errormessagelength

!TYPE definitions
USE prognostics, ONLY: progs_data_type
USE trifctl, ONLY: trifctl_type

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises IMOGEN parameters and arrays
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Arguments
CHARACTER(LEN=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                         ! namelists
TYPE(trifctl_type), INTENT(IN OUT) :: trifctltype

!TYPES containing the data
TYPE(progs_data_type), INTENT(IN OUT) :: progs_data

! Work variables
CHARACTER(LEN=max_file_name_len) :: file_clim  ! File to read IMOGEN
                                               ! climatology from

INTEGER :: i,j,l,n  ! Index variables

INTEGER :: error, error_sum  ! Error indicator
CHARACTER(LEN=errormessagelength) :: iomessage


!-----------------------------------------------------------------------------


IF ( .NOT. l_imogen ) RETURN


! Open the IMOGEN namelist file
OPEN(namelist_unit, FILE=(TRIM(nml_dir) // '/' // 'imogen.nml'),              &
               STATUS='old', POSITION='rewind', ACTION='read', IOSTAT = error,&
               IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_imogen",                                               &
                 "Error opening namelist file imogen.nml " //                 &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")


! There are three namelists to read from this file
CALL log_info("init_imogen", "Reading IMOGEN_RUN_LIST namelist...")
READ(namelist_unit, NML = imogen_run_list, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_imogen",                                               &
                 "Error reading namelist IMOGEN_RUN_LIST " //                 &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

CALL log_info("init_imogen", "Reading IMOGEN_ANLG_VALS_LIST namelist...")
READ(namelist_unit, NML = imogen_anlg_vals_list, IOSTAT = error,              &
     IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_imogen",                                               &
                 "Error reading namelist IMOGEN_ANLG_VALS_LIST " //           &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")


! Close the namelist file
CLOSE(namelist_unit, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_imogen",                                               &
                 "Error closing namelist file imogen.nml " //                 &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")



!-----------------------------------------------------------------------------
! Check that the setup of JULES time variables is compatible with IMOGEN and
! set up IMOGEN time variables
!-----------------------------------------------------------------------------
! -> Force 360 day year
IF ( .NOT. l_360 )                                                            &
  CALL log_fatal("init_imogen", "360 day year must be used with IMOGEN")

! -> Check no spinup has been requested (a full IMOGEN experiment consists of
!    3 JULES runs with the first 2 acting as a really long spinup)
IF ( is_spinup )                                                              &
  CALL log_fatal("init_imogen", "IMOGEN runs do not require spinup")

! -> Run must start at 00:00 on 1st Jan for some year
!    Since we know there is no spinup, we can check this just by checking
!    the main run start time
IF ( main_run_start%time /= 0 .OR. main_run_start%day /= 1 .OR.               &
                                   main_run_start%month /= 1 )                &
  CALL log_fatal("init_imogen",                                               &
                 "IMOGEN runs must start at 00:00 on 1st Jan for some year")

! -> Set steps per day and check it is not more than the maximum number of
!    timesteps per day allowed by IMOGEN
!    Note that we already know that secs_in_day MOD timestep_len is 0
!    since this is enforced in init_time
step_day = secs_in_day / timestep_len
IF (step_day > nsdmax)                                                        &
  CALL log_fatal("init_imogen", "Too many timesteps per day")

!-----------------------------------------------------------------------------
! Find corresponding land sites from the imogen grid 'sgind'.
!-----------------------------------------------------------------------------
CALL get_imogen_map(file_points_order)

!-----------------------------------------------------------------------------
! Allocate imogen arrays
!-----------------------------------------------------------------------------
error_sum = 0
ALLOCATE(t_clim(land_pts,mm), stat = error )
error_sum = error_sum + error
ALLOCATE(rainfall_clim(land_pts,mm), stat = error )
error_sum = error_sum + error
ALLOCATE(snowfall_clim(land_pts,mm), stat = error )
error_sum = error_sum + error
ALLOCATE(rh15m_clim(land_pts,mm), stat = error )
error_sum = error_sum + error
ALLOCATE(uwind_clim(land_pts,mm), stat = error )
error_sum = error_sum + error
ALLOCATE(vwind_clim(land_pts,mm), stat = error )
error_sum = error_sum + error
ALLOCATE(dtemp_clim(land_pts,mm), stat = error )
error_sum = error_sum + error
ALLOCATE(pstar_ha_clim(land_pts,mm), stat = error )
error_sum = error_sum + error
ALLOCATE(sw_clim(land_pts,mm), stat = error )
error_sum = error_sum + error
ALLOCATE(lw_clim(land_pts,mm), stat = error )
error_sum = error_sum + error
ALLOCATE(f_wet_clim(land_pts,mm), stat = error )
error_sum = error_sum + error

ALLOCATE(t_out(land_pts,mm,md,nsdmax), stat = error )
error_sum = error_sum + error
ALLOCATE(conv_rain_out(land_pts,mm,md,nsdmax), stat = error )
error_sum = error_sum + error
ALLOCATE(conv_snow_out(land_pts,mm,md,nsdmax), stat = error )
error_sum = error_sum + error
ALLOCATE(ls_rain_out(land_pts,mm,md,nsdmax), stat = error )
error_sum = error_sum + error
ALLOCATE(ls_snow_out(land_pts,mm,md,nsdmax), stat = error )
error_sum = error_sum + error
ALLOCATE(qhum_out(land_pts,mm,md,nsdmax), stat = error )
error_sum = error_sum + error
ALLOCATE(wind_out(land_pts,mm,md,nsdmax), stat = error )
error_sum = error_sum + error
ALLOCATE(pstar_out(land_pts,mm,md,nsdmax), stat = error )
error_sum = error_sum + error
ALLOCATE(sw_out(land_pts,mm,md,nsdmax), stat = error )
error_sum = error_sum + error
ALLOCATE(lw_out(land_pts,mm,md,nsdmax), stat = error )
error_sum = error_sum + error

ALLOCATE(lat(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE(long(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE(dctot_co2(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE(dctot_ch4(land_pts), stat = error )
error_sum = error_sum + error

! Check for error.
IF ( error_sum /= 0 ) THEN
  CALL log_fatal("init_imogen", "Error allocating IMOGEN arrays")
ELSE
  ! Initialise.
  t_clim(:,:)        = 0.0
  rainfall_clim(:,:) = 0.0
  snowfall_clim(:,:) = 0.0
  rh15m_clim(:,:)    = 0.0
  uwind_clim(:,:)    = 0.0
  vwind_clim(:,:)    = 0.0
  dtemp_clim(:,:)    = 0.0
  pstar_ha_clim(:,:) = 0.0
  sw_clim(:,:)       = 0.0
  lw_clim(:,:)       = 0.0
  f_wet_clim(:,:)    = 0.0

  t_out(:,:,:,:)         = 0.0
  conv_rain_out(:,:,:,:) = 0.0
  conv_snow_out(:,:,:,:) = 0.0
  ls_rain_out(:,:,:,:)   = 0.0
  ls_snow_out(:,:,:,:)   = 0.0
  qhum_out(:,:,:,:)      = 0.0
  wind_out(:,:,:,:)      = 0.0
  pstar_out(:,:,:,:)     = 0.0
  sw_out(:,:,:,:)        = 0.0
  lw_out(:,:,:,:)        = 0.0

  lat(:)   = 0.0
  long(:)  = 0.0
  dctot_co2(:) = 0.0
  dctot_ch4(:) = 0.0
END IF

! Weather generator is not available at present
IF ( wgen )                                                                   &
  CALL log_fatal("init_imogen", "Weather generator not available at present")

!-----------------------------------------------------------------------------
! Check that the configurations proposed are valid.
! This contains the runs that are allowed on the "decision" tree given under
! directory "plots" and called "imogen.jpg"
! This routine also checks whether the selected configuration has been
!  confirmed. Unconfirmed configurations return a WARNING message, but do not
!  break the run.
!-----------------------------------------------------------------------------
CALL imogen_check(                                                            &
  c_emissions, include_co2, include_non_co2, land_feed_co2,                   &
  land_feed_ch4, ocean_feed, anlg, anom, wgen)

!-----------------------------------------------------------------------------
! Now get the emissions/CO2 concentrations.
! At present only coded for reading in a file of emission/CO2 concentration
! imogen projects of the future climate, with carbon cycle feedbacks) or
! in a file of CO2 concentrations for "Hydrology 20th Century" simulation
!-----------------------------------------------------------------------------
IF (nyr_emiss > nyr_max) THEN
  CALL log_fatal("init_imogen", "Time series of emissions is " //             &
                 "longer than allocated array"                 //             &
                 "check nyr_emiss and nyr_max")
END IF

IF ( c_emissions .AND. include_co2 .AND. anom .AND. anlg ) THEN
  OPEN(imogen_unit, FILE=file_scen_emits,                                     &
                    STATUS='old', POSITION='rewind', ACTION='read')

  DO n = 1,nyr_emiss
    READ(imogen_unit,*) yr_emiss(n),c_emiss(n)
  END DO

  CLOSE(imogen_unit)
END IF

!-----------------------------------------------------------------------
! Read in monthly control climate data.
!-----------------------------------------------------------------------
DO j = 1,mm
  file_clim = TRIM(dir_clim) // drive_month(j)

  OPEN(imogen_unit, FILE=file_clim,                                           &
                    STATUS='old', POSITION='rewind', ACTION='read')

  READ(imogen_unit,*) longmin_clim,latmin_clim,longmax_clim, latmax_clim

  DO i = 1,n_imogen_land
    IF (sgindinv(i) > 0) THEN
      l = sgindinv(i)
      READ(imogen_unit,*) long(l), lat(l), t_clim(l,j), rh15m_clim(l,j),      &
                          uwind_clim(l,j), vwind_clim(l,j), lw_clim(l,j),     &
                          sw_clim(l,j), dtemp_clim(l,j),                      &
                          rainfall_clim(l,j), snowfall_clim(l,j),             &
                          pstar_ha_clim(l,j), f_wet_clim(l,j)
    ELSE
      READ(imogen_unit,*)
    END IF
  END DO

  CLOSE(imogen_unit)
END DO

!-----------------------------------------------------------------------------
! Set up the initial conditions
!-----------------------------------------------------------------------------
IF (land_feed_ch4) THEN
  ch4_ppbv = ch4_init_ppbv
END IF

co2_ppmv = co2_init_ppmv

IF ( include_co2 ) co2_change_ppmv = 0.0

IF ( anlg ) THEN
  dtemp_o(:) = 0.0

  IF ( include_co2 .AND. ocean_feed .AND. c_emissions ) fa_ocean(:)=0.0
END IF

! In IMOGEN, the seed has 4 components
ALLOCATE( progs_data%seed_rain(4) )

! Initiate seeding values for subdaily rainfall
progs_data%seed_rain(1) = 9465
progs_data%seed_rain(2) = 1484
progs_data%seed_rain(3) = 3358
progs_data%seed_rain(4) = 8350
! Initiate seeding values for the weather generator
IF ( wgen ) THEN
  seed_wg(1) = 5810
  seed_wg(2) = 5575
  seed_wg(3) = 5817
  seed_wg(4) = 9119
END IF

trifctltype%cv_gb(:) = 0.0

! Override the IMOGEN variables from given dump file if we have been asked to
IF ( initialise_from_dump ) THEN
  CALL read_dump(dump_file, (/ 'co2_ppmv       ', 'co2_change_ppmv',          &
                               'dtemp_o        ', 'fa_ocean       ',          &
                               'seed_rain      ', 'cv             ' /))

  IF (land_feed_ch4) THEN
    CALL read_dump(dump_file, (/ 'ch4_ppbv' /) )
  END IF
END IF

! Unit conversion ppm to mmr (g/g): mol mass co2 / mol mass dry air * 1e-6
co2_mmr = co2_ppmv * 44.0 / 28.97 * 1.0e-6

!-------------------------------------------------------------------------------
! Set diff_frac_const from IMOGEN input
!-------------------------------------------------------------------------------
! Metadata does not allow "OR" triggering so in order to preserve the original
! diff_frac_const functionality it is necessary for IMOGEN to have its own value
diff_frac_const = diff_frac_const_imogen

RETURN

END SUBROUTINE init_imogen
#endif
