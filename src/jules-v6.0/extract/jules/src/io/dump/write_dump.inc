#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

SUBROUTINE write_dump()

USE jules_fields_mod, ONLY: crop_vars, psparms, ainfo, trif_vars, progs

!Science variables
! TYPE Definitions
USE jules_fields_mod, ONLY: toppdm, soilecosse, trifctltype

USE model_grid_mod, ONLY:                                                     &
  global_land_pts, grid_area_ij, latitude, longitude

USE ancil_info, ONLY:                                                         &
  dim_cs1, dim_soil_n_pool, nsurft, nsoilt, dim_cslayer

USE fluxes, ONLY:                                                             &
  t_growth_gb

USE jules_surface_types_mod, ONLY:                                            &
  npft, ntype, ncpft

USE jules_irrig_mod, ONLY: irrtiles, frac_irrig_all_tiles, nirrtile,          &
  set_irrfrac_on_irrtiles

USE imogen_progs, ONLY:                                                       &
  co2_ppmv, co2_change_ppmv, dtemp_o, fa_ocean, ch4_ppbv

USE jules_snow_mod, ONLY:                                                     &
  nsmax

USE jules_soil_mod, ONLY:                                                     &
  sm_levels, ns_deep

USE jules_soil_biogeochem_mod, ONLY:                                          &
  dim_ch4layer

USE fire_mod, ONLY:                                                           &
  fire_prog

USE metstats_mod, ONLY:                                                       &
  metstats_prog

USE aero, ONLY:                                                               &
  co2_mmr

USE jules_rivers_mod, ONLY:                                                   &
  l_rivers, rivers_lat_rp, rivers_lon_rp, rivers_sto_rp, rfm_surfstore_rp,    &
  rfm_substore_rp, rfm_flowin_rp, rfm_bflowin_rp

USE jules_water_resources_mod, ONLY:                                          &
  conveyance_loss, irrig_eff, sfc_water_frac

USE lake_mod,                 ONLY:                                           &
  lake_depth_gb, lake_fetch_gb, lake_t_mean_gb, lake_t_mxl_gb,                &
  lake_h_mxl_gb, lake_t_ice_gb, lake_h_ice_gb, lake_shape_factor_gb

!Others
USE io_constants, ONLY:                                                       &
  max_file_name_len, max_dim_var, mode_write

USE parallel_mod, ONLY:                                                       &
  is_master_task, gather_land_field

USE model_interface_mod, ONLY:                                                &
  identifier_len

USE dictionary_mod, ONLY:                                                     &
  dict, dict_create, dict_set, dict_get, dict_has_key, dict_free

USE string_utils_mod, ONLY:                                                   &
  to_string

USE file_mod, ONLY:                                                           &
  file_handle, file_open, file_def_dim, file_def_var,                         &
  file_enddef, file_write_var, file_close

USE output_mod, ONLY:                                                         &
  output_dir, run_id

USE model_time_mod, ONLY:                                                     &
  current_time, is_spinup, spinup_cycle

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Writes a dump file for the current timestep
!   Note that the writing of the dump is done by the master task with the
!   values gathered from other tasks. This means that dumps written with
!   different amounts of tasks should be interchangable.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Local parameters.
LOGICAL, PARAMETER :: l_reading_false = .FALSE.
  ! A value of .FALSE. that is passed to argument l_reading of subroutine
  ! get_dim_info to show that it is being called in connection with writing
  ! (not reading) a dump.

! Work variables
CHARACTER(LEN=max_file_name_len) :: file_name
                                    ! The filename to use for the dump file
CHARACTER(LEN=max_file_name_len) :: dt_string
                                    ! The datetime string to use in the file
                                    ! name
CHARACTER(LEN=identifier_len) :: identifiers(max_var_dump)
                                    ! The model identifiers for the variables
                                    ! to put in the dump
CHARACTER(LEN=identifier_len) :: vars_from_ancil(max_var_dump)
                             ! The variable identifiers of the ancil
                             ! variables (not used in this subroutine)

TYPE(file_handle) :: FILE  ! The dump file

INTEGER :: nvars  ! The number of variables we are processing
INTEGER :: nvars_from_ancil

! Variables used when defining dimensions
TYPE(dict) :: file_dim_ids  ! Dictionary of the dimensions that have been
                            ! defined
                            ! Maps dim_name => dim_id

INTEGER :: ndims  ! The number of levels dims for the current variable
CHARACTER(LEN=max_sdf_name_len) :: dim_names(max_dim_var)
                  ! The levels dimension names for the current variable
INTEGER :: dim_sizes(max_dim_var)
                  ! The sizes of the levels dims for the current variable
INTEGER :: dim_ids(max_dim_var)
                  ! The ids in file of the levels dims for the current
                  ! variable
INTEGER :: var_ids(max_var_dump)
                      ! The ids of the variables in the dump file

INTEGER :: i, j, m, n, s  ! Loop counters

REAL :: frac_irrig_all_tiles_real
REAL :: set_irrfrac_irrtiles_real
REAL :: irrtiles_real(npft)
REAL :: nirrtile_real

! Arrays to hold global land points version of data gathered in master task
! before writing
REAL, ALLOCATABLE :: global_data_1d(:)  ! For data with no vertical levels
REAL, ALLOCATABLE :: global_data_2d(:,:)   ! With one vertical level
REAL, ALLOCATABLE :: global_data_3d(:,:,:) ! With two "vertical" levels
                                        ! I.E. snow variables or soil C/N
REAL, ALLOCATABLE :: global_data_4d(:,:,:,:)

LOGICAL, PARAMETER :: l_output_mode = .TRUE.

!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Get the list of identifiers that we are going to output
!-----------------------------------------------------------------------------
CALL required_vars_for_configuration(nvars, identifiers,                      &
                                     nvars_from_ancil, vars_from_ancil,       &
                                     l_output_mode)

! Add latitude and longitude to the list to help offline inspection of the
! dump file, ie they are diagnostics, not prognostics or ancillaries.
! Lat & lon will not be read in from the dump to prevent confusion with the
! model grid namelists

nvars = nvars + 1
identifiers(nvars) = 'latitude'
nvars = nvars + 1
identifiers(nvars) = 'longitude'

! Similarly, add latitude and longitude of river points, if needed.
IF ( l_rivers ) THEN
  nvars = nvars + 1
  identifiers(nvars) = 'rivers_lat_rp'
  nvars = nvars + 1
  identifiers(nvars) = 'rivers_lon_rp'
END IF

!-----------------------------------------------------------------------------
! In the master task only, we open a new file and define the required
! dimensions and variables
!-----------------------------------------------------------------------------
IF ( is_master_task() ) THEN
  !---------------------------------------------------------------------------
  ! Generate the file name that we want to use and open the file
  !---------------------------------------------------------------------------
  ! File name starts with run id + indicator of a dump file
  file_name = TRIM(run_id) // ".dump."

  ! Include the current spinup cycle if there is one
  IF ( is_spinup )                                                            &
    file_name = TRIM(file_name) //                                            &
                "spin" // TRIM(to_string(spinup_cycle)) // "."

  ! Then current date and time
  WRITE(dt_string, '(I4.4,I2.2,I2.2)') current_time%year,                     &
                                       current_time%month,                    &
                                       current_time%day
  dt_string = TRIM(dt_string) // "." // TRIM(to_string(current_time%time))
  file_name = TRIM(file_name) // TRIM(dt_string)

  ! Add the extension based on dump format
  SELECT CASE ( dump_format )
  CASE ( format_ascii )
    file_name = TRIM(file_name) // ".asc"

  CASE ( format_ncdf )
    file_name = TRIM(file_name) // ".nc"

  CASE DEFAULT
    CALL log_fatal("write_dump",                                              &
                   "Unrecognised file format - " // TRIM(dump_format))
  END SELECT

  ! Prepend the output directory
  file_name = TRIM(output_dir) // "/" // TRIM(file_name)

  ! We use the lowest level file API here, as we don't want to impose a grid
  FILE=file_open(file_name, mode_write)

  !---------------------------------------------------------------------------
  ! Create the dimensions and variables
  !---------------------------------------------------------------------------
  file_dim_ids = dict_create(max_dim_dump, INT(1))

  DO i = 1,nvars

    !------------------------------------------------------------------------
    ! Get information about the dimensions used by the variable.
    ! The argument l_reading_false shows that we are writing (not reading) a
    ! dump.
    !------------------------------------------------------------------------
    CALL get_dim_info( l_reading_false, identifiers(i), ndims,  dim_sizes,    &
                       dim_names )

    !-------------------------------------------------------------------------
    ! Define the dimensions if they have not already been defined
    ! We use a dictionary to keep track of defined dimension ids
    !
    ! At the same time, gather up the dimension ids needed by the current
    ! variable.
    !-------------------------------------------------------------------------
    DO j = 1,ndims
      ! If it has not yet been defined, define the dimension, storing its id
      IF ( .NOT. dict_has_key(file_dim_ids, dim_names(j)) )                   &
        CALL dict_set(                                                        &
          file_dim_ids, dim_names(j),                                         &
          file_def_dim(FILE, dim_names(j), dim_sizes(j))                      &
        )

      ! Get the dimension id from the dict and add it to the list for this
      ! variable.
      CALL dict_get(file_dim_ids, dim_names(j), dim_ids(j))
    END DO

    !-------------------------------------------------------------------------
    ! Define the variable, saving the id in the file for later
    !-------------------------------------------------------------------------
    var_ids(i) = file_def_var(FILE, identifiers(i), dim_ids(1:ndims),         &
                              .FALSE.)

  END DO

  !---------------------------------------------------------------------------
  ! We have finished defining things
  !---------------------------------------------------------------------------
  CALL file_enddef(FILE)
  CALL dict_free(file_dim_ids)

END IF  ! MASTER TASK


!-----------------------------------------------------------------------------
! Gather data from other tasks and write it to file
!-----------------------------------------------------------------------------
! Allocate the global data arrays
IF ( is_master_task() ) THEN
  ALLOCATE(global_data_1d(global_land_pts))
  ALLOCATE(global_data_2d(global_land_pts, MAX(npft, sm_levels,               &
                                               nsurft, ntype, ns_deep,        &
                                               nsoilt, dim_ch4layer)))
  ALLOCATE(global_data_3d(global_land_pts,                                    &
                          MAX(nsurft, dim_cslayer, nsoilt),                   &
                          MAX(nsmax, dim_cs1, dim_soil_n_pool, sm_levels)))
  ALLOCATE(global_data_4d(global_land_pts, nsoilt,                            &
                          MAX(nsurft, dim_cslayer),                           &
                          MAX(nsmax, dim_cs1, dim_soil_n_pool, sm_levels)))
ELSE
  ALLOCATE(global_data_1d(1))
  ALLOCATE(global_data_2d(1,1))
  ALLOCATE(global_data_3d(1,1,1))
END IF

DO i = 1,nvars
  ! Gather the variables into a global array to write to file
  ! Note that gather_land_field can only gather one land_pts array at a time,
  ! so to gather variables with multiple levels we must loop

  CALL log_info("write_dump", identifiers(i))

  SELECT CASE ( identifiers(i) )
  CASE ( 'gs' )
    CALL gather_land_field(progs%gs_gb, global_data_1d)

    !Case if nsoilt == 1, so it is OK to hardwire the 2nd dimension to 1
  CASE ( 'sthzw' )
    CALL gather_land_field(toppdm%sthzw_soilt(:,1), global_data_1d)

  CASE ( 'toppdm%sthzw_soilt' )
    DO m = 1,nsoilt
      CALL gather_land_field(toppdm%sthzw_soilt(:,m), global_data_2d(:,m))
    END DO

    !Case if nsoilt == 1, so it is OK to hardwire the 2nd dimension to 1
  CASE ( 'zw' )
    CALL gather_land_field(toppdm%zw_soilt(:,1), global_data_1d)

  CASE ( 'toppdm%zw_soilt' )
    DO m = 1,nsoilt
      CALL gather_land_field(toppdm%zw_soilt(:,m), global_data_2d(:,m))
    END DO

  CASE ( 'cv' )
    CALL gather_land_field(trifctltype%cv_gb, global_data_1d)

  CASE ( 'frac_agr_prev' )
    CALL gather_land_field(progs%frac_agr_prev_gb, global_data_1d)

  CASE ( 'frac_past_prev' )
    CALL gather_land_field(progs%frac_past_prev_gb, global_data_1d)

  CASE ( 'wood_prod_fast' )
    CALL gather_land_field(progs%wood_prod_fast_gb, global_data_1d)

  CASE ( 'wood_prod_med' )
    CALL gather_land_field(progs%wood_prod_med_gb, global_data_1d)

  CASE ( 'wood_prod_slow' )
    CALL gather_land_field(progs%wood_prod_slow_gb, global_data_1d)

  CASE ( 'lake_fetch_gb' )
    CALL gather_land_field(lake_fetch_gb, global_data_1d)

  CASE ( 'lake_t_mean_gb' )
    CALL gather_land_field(lake_t_mean_gb, global_data_1d)

  CASE ( 'lake_t_mxl_gb' )
    CALL gather_land_field(lake_t_mxl_gb, global_data_1d)

  CASE ( 'lake_h_mxl_gb' )
    CALL gather_land_field(lake_h_mxl_gb, global_data_1d)

  CASE ( 'lake_t_ice_gb' )
    CALL gather_land_field(lake_t_ice_gb, global_data_1d)

  CASE ( 'lake_h_ice_gb' )
    CALL gather_land_field(lake_h_ice_gb, global_data_1d)

  CASE ( 'lake_shape_factor_gb' )
    CALL gather_land_field(lake_shape_factor_gb, global_data_1d)

  CASE ( 'latitude' )
    CALL gather_land_field(latitude, global_data_1d)

  CASE ( 'longitude' )
    CALL gather_land_field(longitude, global_data_1d)

    !Case if nsoilt == 1, so it is OK to hardwire the 2nd dimension to 1
  CASE ( 'n_inorg' )
    DO n = 1,dim_cslayer
      CALL gather_land_field(progs%n_inorg_soilt_lyrs(:,1,n),                 &
                             global_data_2d(:,n))
    END DO

  CASE ( 'n_inorg_soilt' )
    DO m = 1,nsoilt
      DO n = 1,dim_cslayer
        CALL gather_land_field(progs%n_inorg_soilt_lyrs(:,m,n),               &
                               global_data_3d(:,m,n))
      END DO
    END DO

  CASE ( 'substr_ch4' )
    DO n = 1,dim_ch4layer
      CALL gather_land_field(progs%substr_ch4(:,n),                           &
                             global_data_2d(:,n))
    END DO

  CASE ( 'mic_ch4' )
    DO n = 1,dim_ch4layer
      CALL gather_land_field(progs%mic_ch4(:,n),                              &
                             global_data_2d(:,n))
    END DO

  CASE ( 'mic_act_ch4' )
    DO n = 1,dim_ch4layer
      CALL gather_land_field(progs%mic_act_ch4(:,n),                          &
                             global_data_2d(:,n))
    END DO

  CASE ( 'acclim_ch4' )
    DO n = 1,dim_ch4layer
      CALL gather_land_field(progs%acclim_ch4(:,n),                           &
                             global_data_2d(:,n))
    END DO

  CASE ( 'canht' )
    DO n = 1,npft
      CALL gather_land_field(progs%canht_pft(:,n), global_data_2d(:,n))
    END DO

  CASE ( 'lai' )
    DO n = 1,npft
      CALL gather_land_field(progs%lai_pft(:,n), global_data_2d(:,n))
    END DO

    !Case if nsoilt == 1, so it is OK to hardwire the 2nd dimension to 1
  CASE ( 'cs' )
    DO n = 1,dim_cs1
      DO m = 1,dim_cslayer
        CALL gather_land_field(progs%cs_pool_soilt(:,1,m,n),                  &
                               global_data_3d(:,m,n))
      END DO
    END DO

  CASE ( 'cs_soilt' )
    DO s = 1,nsoilt
      DO n = 1,dim_cs1
        DO m = 1,dim_cslayer
          CALL gather_land_field(progs%cs_pool_soilt(:,s,m,n),                &
                                 global_data_4d(:,s,m,n))
        END DO
      END DO
    END DO

  CASE ( 'cropdvi' )
    DO n = 1,ncpft
      CALL gather_land_field(crop_vars%dvi_cpft(:,n), global_data_2d(:,n))
    END DO

  CASE ( 'croprootc' )
    DO n = 1,ncpft
      CALL gather_land_field(crop_vars%rootc_cpft(:,n), global_data_2d(:,n))
    END DO

  CASE ( 'cropharvc' )
    DO n = 1,ncpft
      CALL gather_land_field(crop_vars%harvc_cpft(:,n), global_data_2d(:,n))
    END DO

  CASE ( 'cropreservec' )
    DO n = 1,ncpft
      CALL gather_land_field(crop_vars%reservec_cpft(:,n), global_data_2d(:,n))
    END DO

  CASE ( 'croplai' )
    DO n = 1,ncpft
      CALL gather_land_field(crop_vars%croplai_cpft(:,n), global_data_2d(:,n))
    END DO

  CASE ( 'cropcanht' )
    DO n = 1,ncpft
      CALL gather_land_field(crop_vars%cropcanht_cpft(:,n), global_data_2d(:,n))
    END DO

  CASE ( 'ns' )
    DO n = 1,dim_cs1
      DO m = 1,dim_cslayer
        CALL gather_land_field(progs%ns_pool_gb(:,m,n), global_data_3d(:,m,n))
      END DO
    END DO

    !Case if nsoilt == 1, so it is OK to hardwire the 2nd dimension to 1
  CASE ( 'sthuf' )
    ! sthuf is held in sthu until it is processed
    DO n = 1,sm_levels
      CALL gather_land_field(psparms%sthu_soilt(:,1,n) +                      &
                             psparms%sthf_soilt(:,1,n),                       &
                             global_data_2d(:,n))
    END DO

  CASE ( 'sthuf_soilt' )
    ! sthuf is held in sthu until it is processed
    DO m = 1,nsoilt
      DO n = 1,sm_levels
        CALL gather_land_field(psparms%sthu_soilt(:,m,n) +                    &
                               psparms%sthf_soilt(:,m,n),                     &
                               global_data_3d(:,m,n))
      END DO
    END DO

    !Case if nsoilt == 1, so it is OK to hardwire the 2nd dimension to 1
  CASE ( 't_soil' )
    DO n = 1,sm_levels
      CALL gather_land_field(progs%t_soil_soilt(:,1,n), global_data_2d(:,n))
    END DO

  CASE ( 't_soil_soilt' )
    DO m = 1,nsoilt
      DO n = 1,sm_levels
        CALL gather_land_field(progs%t_soil_soilt(:,m,n), global_data_3d(:,m,n))
      END DO
    END DO

  CASE ( 'tsoil_deep' )
    DO n = 1,ns_deep
      CALL gather_land_field(progs%tsoil_deep_gb(:,n), global_data_2d(:,n))
    END DO

    !Case if nsoilt == 1, so it is OK to hardwire the 2nd dimension to 1
  CASE ( 'sthu_irr' )
    DO n = 1,sm_levels
      CALL gather_land_field(crop_vars%sthu_irr_soilt(:,1,n),                 &
                             global_data_2d(:,n))
    END DO

  CASE ( 'sthu_irr_soilt' )
    DO m = 1,nsoilt
      DO n = 1,sm_levels
        CALL gather_land_field(crop_vars%sthu_irr_soilt(:,m,n),               &
                               global_data_3d(:,m,n))
      END DO
    END DO

  CASE ( 'canopy' )
    DO n = 1,nsurft
      CALL gather_land_field(progs%canopy_surft(:,n), global_data_2d(:,n))
    END DO

  CASE ( 'nsnow' )
    DO n = 1,nsurft
      CALL gather_land_field(REAL(progs%nsnow_surft(:,n)), global_data_2d(:,n))
    END DO

  CASE ( 'rgrain' )
    DO n = 1,nsurft
      CALL gather_land_field(progs%rgrain_surft(:,n), global_data_2d(:,n))
    END DO

  CASE ( 'rho_snow' )
    DO n = 1,nsurft
      CALL gather_land_field(progs%rho_snow_grnd_surft(:,n),                  &
                             global_data_2d(:,n))
    END DO

  CASE ( 'snow_tile' )
    DO n = 1,nsurft
      CALL gather_land_field(progs%snow_surft(:,n), global_data_2d(:,n))
    END DO

  CASE ( 'snow_depth' )
    DO n = 1,nsurft
      CALL gather_land_field(progs%snowdepth_surft(:,n), global_data_2d(:,n))
    END DO

  CASE ( 'snow_grnd' )
    DO n = 1,nsurft
      CALL gather_land_field(progs%snow_grnd_surft(:,n), global_data_2d(:,n))
    END DO

  CASE ( 'tstar_tile' )
    DO n = 1,nsurft
      CALL gather_land_field(progs%tstar_surft(:,n), global_data_2d(:,n))
    END DO

  CASE ( 'tsurf_elev_surft' )
    DO n = 1,nsurft
      CALL gather_land_field(progs%tsurf_elev_surft(:,n), global_data_2d(:,n))
    END DO

  CASE ( 'rgrainl' )
    DO n = 1,nsmax
      DO m = 1,nsurft
        CALL gather_land_field(progs%rgrainl_surft(:,m,n),                    &
                               global_data_3d(:,m,n))
      END DO
    END DO

  CASE ( 'snow_ds' )
    DO n = 1,nsmax
      DO m = 1,nsurft
        CALL gather_land_field(progs%ds_surft(:,m,n), global_data_3d(:,m,n))
      END DO
    END DO

  CASE ( 'snow_ice' )
    DO n = 1,nsmax
      DO m = 1,nsurft
        CALL gather_land_field(progs%sice_surft(:,m,n), global_data_3d(:,m,n))
      END DO
    END DO


  CASE ( 'snow_liq' )
    DO n = 1,nsmax
      DO m = 1,nsurft
        CALL gather_land_field(progs%sliq_surft(:,m,n), global_data_3d(:,m,n))
      END DO
    END DO


  CASE ( 'tsnow' )
    DO n = 1,nsmax
      DO m = 1,nsurft
        CALL gather_land_field(progs%tsnow_surft(:,m,n), global_data_3d(:,m,n))
      END DO
    END DO

    ! Cases for metstats variables
  CASE ('temp_max_00h_r')
    CALL gather_land_field(metstats_prog(:)%temp_max_00h%run,                 &
                           global_data_1d)

  CASE ('temp_ave_00h_r')
    CALL gather_land_field(metstats_prog(:)%temp_ave_00h%run,                 &
                           global_data_1d)

  CASE ('prec_tot_00h_r')
    CALL gather_land_field(metstats_prog(:)%prec_tot_00h%run,                 &
                           global_data_1d)

  CASE ('prec_tot_12h_r')
    CALL gather_land_field(metstats_prog(:)%prec_tot_12h%run,                 &
                           global_data_1d)

  CASE ('rhum_min_00h_r')
    CALL gather_land_field(metstats_prog(:)%rhum_min_00h%run,                 &
                           global_data_1d)

  CASE ('dewp_ave_00h_r')
    CALL gather_land_field(metstats_prog(:)%dewp_ave_00h%run,                 &
                           global_data_1d)

  CASE ('wind_ave_00h_r')
    CALL gather_land_field(metstats_prog(:)%wind_ave_00h%run,                 &
                           global_data_1d)

  CASE ('temp_max_00h')
    CALL gather_land_field(metstats_prog(:)%temp_max_00h%fin,                 &
                           global_data_1d)

  CASE ('temp_ave_00h')
    CALL gather_land_field(metstats_prog(:)%temp_ave_00h%fin,                 &
                           global_data_1d)

  CASE ('temp_ave_nday')
    CALL gather_land_field(metstats_prog(:)%temp_ave_nday%fin,                &
                           global_data_1d)

  CASE ('temp_pnt_12h')
    CALL gather_land_field(metstats_prog(:)%temp_pnt_12h%fin,                 &
                           global_data_1d)

  CASE ('prec_tot_00h')
    CALL gather_land_field(metstats_prog(:)%prec_tot_00h%fin,                 &
                           global_data_1d)

  CASE ('prec_tot_12h')
    CALL gather_land_field(metstats_prog(:)%prec_tot_12h%fin,                 &
                           global_data_1d)

  CASE ('rhum_min_00h')
    CALL gather_land_field(metstats_prog(:)%rhum_min_00h%fin,                 &
                           global_data_1d)

  CASE ('rhum_pnt_12h')
    CALL gather_land_field(metstats_prog(:)%rhum_pnt_12h%fin,                 &
                           global_data_1d)

  CASE ('dewp_ave_00h')
    CALL gather_land_field(metstats_prog(:)%dewp_ave_00h%fin,                 &
                           global_data_1d)

  CASE ('wind_ave_00h')
    CALL gather_land_field(metstats_prog(:)%wind_ave_00h%fin,                 &
                           global_data_1d)

  CASE ('wind_pnt_12h')
    CALL gather_land_field(metstats_prog(:)%wind_pnt_12h%fin,                 &
                           global_data_1d)

    ! Cases for Fire variables
  CASE ( 'fire_mcarthur_r_dr' )
    CALL gather_land_field(fire_prog(:)%mcarthur%r_dr, global_data_1d)

  CASE ( 'fire_mcarthur_n_dr' )
    CALL gather_land_field(fire_prog(:)%mcarthur%n_dr, global_data_1d)

  CASE ( 'fire_canadian_ffmc' )
    CALL gather_land_field(fire_prog(:)%canadian%ffmc, global_data_1d)

  CASE ( 'fire_canadian_ffmc_mois' )
    CALL gather_land_field(fire_prog(:)%canadian%ffmc_mois,                   &
                           global_data_1d)

  CASE ( 'fire_canadian_dmc' )
    CALL gather_land_field(fire_prog(:)%canadian%dmc, global_data_1d)

  CASE ( 'fire_canadian_dc' )
    CALL gather_land_field(fire_prog(:)%canadian%dc, global_data_1d)

  CASE ( 'fire_nesterov' )
    CALL gather_land_field(fire_prog(:)%nesterov%findex, global_data_1d)

    ! ECOSSE variables
    !Case if nsoilt == 1, so it is OK to hardwire the 2nd dimension to 1
  CASE ( 'n_soil' )
    DO n = 1,dim_soil_n_pool
      DO m = 1,dim_cslayer
        CALL gather_land_field(soilecosse%n_soil_pool_soilt(:,1,m,n),         &
                               global_data_3d(:,m,n))
      END DO
    END DO

  CASE ( 'n_soil_soilt' )
    DO s = 1,nsoilt
      DO n = 1,dim_soil_n_pool
        DO m = 1,dim_cslayer
          CALL gather_land_field(progs%cs_pool_soilt(:,s,m,n),                &
                                 global_data_4d(:,s,m,n))
        END DO
      END DO
    END DO

    ! Since each task runs its own version of IMOGEN, we just use the values
    ! from the master task.
  CASE ( 'co2_ppmv', 'co2_change_ppmv', 'dtemp_o', 'fa_ocean',                &
         'seed_rain', 'ch4_ppbv' )
    ! Nothing to do

  CASE ( 'rivers_lat_rp', 'rivers_lon_rp', 'rivers_sto_rp',                   &
         'rfm_surfstore_rp', 'rfm_substore_rp',                               &
         'rfm_flowin_rp', 'rfm_bflowin_rp' )
    ! Nothing to do

    !-------------------------------------------------------------------------
    ! Ancillary variables
    !-------------------------------------------------------------------------

    ! latlon ancil namelist
  CASE ( 'grid_area' )
    CALL gather_land_field(grid_area_ij, global_data_1d)

    !Frac ancil namelist
  CASE ( 'frac' )
    DO n = 1,ntype
      CALL gather_land_field(ainfo%frac_surft(:,n), global_data_2d(:,n))
    END DO

    ! Vegetation properties ancil namelist
  CASE ( 't_growth_gb' )
    CALL gather_land_field(t_growth_gb, global_data_1d)

    !Soil properties ancil namelist
    !Cases if nsoilt == 1, so it is OK to hardwire the 2nd dimension to 1
  CASE ( 'b      ')
    DO n = 1,sm_levels
      CALL gather_land_field(psparms%bexp_soilt(:,1,n), global_data_2d(:,n))
    END DO

  CASE ( 'b_soilt')
    DO m = 1,nsoilt
      DO n = 1,sm_levels
        CALL gather_land_field(psparms%bexp_soilt(:,m,n), global_data_3d(:,m,n))
      END DO
    END DO

  CASE ( 'sathh  ')
    DO n = 1,sm_levels
      CALL gather_land_field(psparms%sathh_soilt(:,1,n), global_data_2d(:,n))
    END DO

  CASE ( 'sathh_soilt')
    DO m = 1,nsoilt
      DO n = 1,sm_levels
        CALL gather_land_field(psparms%sathh_soilt(:,m,n),                    &
                               global_data_3d(:,m,n))
      END DO
    END DO

  CASE ( 'satcon ')
    DO n = 1,sm_levels
      CALL gather_land_field(psparms%satcon_soilt(:,1,n), global_data_2d(:,n))
    END DO

  CASE ( 'satcon_soilt')
    DO m = 1,nsoilt
      DO n = 1,sm_levels
        CALL gather_land_field(psparms%satcon_soilt(:,m,n),                   &
                               global_data_3d(:,m,n))
      END DO
    END DO

  CASE ( 'sm_sat ')
    DO n = 1,sm_levels
      CALL gather_land_field(psparms%smvcst_soilt(:,1,n), global_data_2d(:,n))
    END DO

  CASE ( 'sm_sat_soilt')
    DO m = 1,nsoilt
      DO n = 1,sm_levels
        CALL gather_land_field(psparms%smvcst_soilt(:,m,n),                   &
                               global_data_3d(:,m,n))
      END DO
    END DO

  CASE ( 'sm_crit')
    DO n = 1,sm_levels
      CALL gather_land_field(psparms%smvccl_soilt(:,1,n), global_data_2d(:,n))
    END DO

  CASE ( 'sm_crit_soilt')
    DO m = 1,nsoilt
      DO n = 1,sm_levels
        CALL gather_land_field(psparms%smvccl_soilt(:,m,n),                   &
                               global_data_3d(:,m,n))
      END DO
    END DO

  CASE ( 'sm_wilt')
    DO n = 1,sm_levels
      CALL gather_land_field(psparms%smvcwt_soilt(:,1,n),                     &
                             global_data_2d(:,n))
    END DO

  CASE ( 'sm_wilt_soilt')
    DO m = 1,nsoilt
      DO n = 1,sm_levels
        CALL gather_land_field(psparms%smvcwt_soilt(:,m,n),                   &
                               global_data_3d(:,m,n))
      END DO
    END DO

  CASE ( 'hcap   ')
    DO n = 1,sm_levels
      CALL gather_land_field(psparms%hcap_soilt(:,1,n), global_data_2d(:,n))
    END DO

  CASE ( 'hcap_soilt')
    DO m = 1,nsoilt
      DO n = 1,sm_levels
        CALL gather_land_field(psparms%hcap_soilt(:,m,n), global_data_3d(:,m,n))
      END DO
    END DO

  CASE ( 'hcon   ' )
    DO n = 1,sm_levels
      CALL gather_land_field(psparms%hcon_soilt(:,1,n), global_data_2d(:,n))
    END DO

  CASE ( 'hcon_soilt')
    DO m = 1,nsoilt
      DO n = 1,sm_levels
        CALL gather_land_field(psparms%hcon_soilt(:,m,n), global_data_3d(:,m,n))
      END DO
    END DO

  CASE ( 'albsoil' )
    CALL gather_land_field(psparms%albsoil_soilt(:,1), global_data_1d)

  CASE ( 'albsoil_soilt' )
    DO m = 1,nsoilt
      CALL gather_land_field(psparms%albsoil_soilt(:,m), global_data_2d(:,m))
    END DO

  CASE ( 'clay' )
    DO n = 1,dim_cslayer
      CALL gather_land_field(psparms%clay_soilt(:,1,n),global_data_2d(:,n))
    END DO

  CASE ( 'clay_soilt' )
    DO m = 1,nsoilt
      DO n = 1,dim_cslayer
        CALL gather_land_field(psparms%clay_soilt(:,m,n),                     &
                               global_data_3d(:,m,n))
      END DO
    END DO

  CASE ( 'soil_ph' )
    DO n = 1,dim_cslayer
      CALL gather_land_field(psparms%soil_ph_soilt(:,1,n),global_data_2d(:,n))
    END DO

  CASE ( 'soil_ph_soilt' )
    DO m = 1,nsoilt
      DO n = 1,dim_cslayer
        CALL gather_land_field(psparms%soil_ph_soilt(:,m,n),                  &
                               global_data_3d(:,m,n))
      END DO
    END DO

    !Topmodel ancillaries namelist
    !Case if nsoilt == 1, so it is OK to hardwire the 2nd dimension to 1
  CASE ( 'fexp   ' )
    CALL gather_land_field(toppdm%fexp_soilt(:,1), global_data_1d)

  CASE ( 'toppdm%fexp_soilt' )
    DO m = 1,nsoilt
      CALL gather_land_field(toppdm%fexp_soilt(:,m), global_data_2d(:,m))
    END DO

    !Case if nsoilt == 1, so it is OK to hardwire the 2nd dimension to 1
  CASE ( 'ti_mean' )
    CALL gather_land_field(toppdm%ti_mean_soilt(:,1), global_data_1d)

  CASE ( 'toppdm%ti_mean_soilt' )
    DO m = 1,nsoilt
      CALL gather_land_field(toppdm%ti_mean_soilt(:,m), global_data_2d(:,m))
    END DO

    !Case if nsoilt == 1, so it is OK to hardwire the 2nd dimension to 1
  CASE ( 'ti_sig ' )
    CALL gather_land_field(toppdm%ti_sig_soilt(:,1), global_data_1d)

  CASE ( 'toppdm%ti_sig_soilt' )
    DO m = 1,nsoilt
      CALL gather_land_field(toppdm%ti_sig_soilt(:,m), global_data_2d(:,m))
    END DO

    !Agric ancillaries namelist
  CASE ( 'frac_agr' )
    CALL gather_land_field(trifctltype%frac_agr_gb, global_data_1d)

  CASE ( 'frac_past' )
    CALL gather_land_field(trif_vars%frac_past_gb, global_data_1d)

    !Crop props ancillaries namelist
  CASE ( 'cropsowdate       ' )
    DO n = 1,ncpft
      CALL gather_land_field(crop_vars%sow_date_cpft(:,n), global_data_2d(:,n))
    END DO

  CASE ( 'croplatestharvdate' )
    DO n = 1,ncpft
      CALL gather_land_field(crop_vars%latestharv_date_cpft(:,n),             &
                             global_data_2d(:,n))
    END DO

  CASE ( 'cropttveg         ' )
    DO n = 1,ncpft
      CALL gather_land_field(crop_vars%tt_veg_cpft(:,n), global_data_2d(:,n))
    END DO

  CASE ( 'cropttrep         ')
    DO n = 1,ncpft
      CALL gather_land_field(crop_vars%tt_rep_cpft(:,n), global_data_2d(:,n))
    END DO

    !Irrigation ancillaries namelist
  CASE ( 'frac_irrig' )
    CALL gather_land_field(crop_vars%frac_irr_all(:,1), global_data_1d)

  CASE ( 'irrfrac_irrtiles' )
    CALL gather_land_field(crop_vars%irrfrac_irrtiles(:,1), global_data_1d)

  CASE ( 'frac_irr_all_tiles' )
    ! Convert to a real
    IF ( frac_irrig_all_tiles ) THEN
      frac_irrig_all_tiles_real = 1.0
    ELSE
      frac_irrig_all_tiles_real = 0.0
    END IF

  CASE ( 'irrtiles' )
    irrtiles_real = REAL(irrtiles(1:npft))

  CASE ( 'nirrtile' )
    nirrtile_real = REAL(nirrtile)

  CASE ( 'set_irrfrac_on_irrtiles' )
    ! Convert to a real
    IF ( set_irrfrac_on_irrtiles ) THEN
      set_irrfrac_irrtiles_real = 1.0
    ELSE
      set_irrfrac_irrtiles_real = 0.0
    END IF

    ! Water resources properties ancil namelist
  CASE ( 'conveyance_loss' )
    CALL gather_land_field(conveyance_loss, global_data_1d)

  CASE ( 'irrig_eff' )
    CALL gather_land_field(irrig_eff, global_data_1d)

  CASE ( 'sfc_water_frac' )
    CALL gather_land_field(sfc_water_frac, global_data_1d)

    !CO2 ancil namelist
  CASE ( 'co2_mmr' )
    !Nothing to do

    !FLake ancil namelist
  CASE ( 'lake_depth_gb' )
    CALL gather_land_field(lake_depth_gb, global_data_1d)

  CASE DEFAULT
    CALL log_fatal("write_dump",                                              &
                   "No code to gather variable for dump - " //                &
                   TRIM(identifiers(i)))
  END SELECT

  !---------------------------------------------------------------------------
  ! In the master task, write the global data to file
  !---------------------------------------------------------------------------
  IF ( is_master_task() ) THEN

    SELECT CASE ( identifiers(i) )
      ! If it is a land_pts array with no levels associated,
      ! write the global_data_1d array
    CASE ( 'gs', 'sthzw', 'zw', 'cv', 'frac_agr_prev', 'frac_past_prev',      &
           'wood_prod_fast', 'wood_prod_med', 'wood_prod_slow',               &
           'temp_max_00h_r', 'temp_ave_00h_r', 'prec_tot_00h_r',              &
           'prec_tot_12h_r', 'rhum_min_00h_r', 'dewp_ave_00h_r',              &
           'wind_ave_00h_r', 'temp_max_00h',   'temp_ave_00h',                &
           'temp_ave_nday',                                                   &
           'temp_pnt_12h',   'prec_tot_00h',   'prec_tot_12h',                &
           'rhum_min_00h',   'rhum_pnt_12h',   'dewp_ave_00h',                &
           'wind_ave_00h',   'wind_pnt_12h',                                  &
           'fire_mcarthur_r_dr', 'fire_mcarthur_n_dr',                        &
           'fire_canadian_ffmc', 'fire_canadian_ffmc_mois',                   &
           'fire_canadian_dmc',  'fire_canadian_dc',                          &
           'fire_nesterov', 'lake_fetch_gb', 'lake_t_mean_gb',                &
           'lake_t_mxl_gb', 'lake_h_mxl_gb', 'lake_t_ice_gb',                 &
           'lake_h_ice_gb', 'lake_shape_factor_gb', 'latitude', 'longitude')
      CALL file_write_var(FILE, var_ids(i), global_data_1d)

    CASE ( 'toppdm%sthzw_soilt', 'toppdm%zw_soilt' )
      CALL file_write_var(FILE, var_ids(i), global_data_2d(:,1:nsoilt))

      ! If it is a variable with one levels dimension, write the appropriate
      ! number of levels to global_data_2d.
    CASE ( 'canht', 'lai' )
      CALL file_write_var(FILE, var_ids(i), global_data_2d(:,1:npft))

    CASE ( 'cropdvi', 'croprootc', 'cropharvc', 'cropreservec',               &
           'croplai', 'cropcanht' )
      CALL file_write_var(FILE, var_ids(i), global_data_2d(:,1:ncpft))

    CASE ( 'sthuf', 't_soil', 'sthu_irr' )
      CALL file_write_var(FILE, var_ids(i),                                   &
                          global_data_2d(:,1:sm_levels))

    CASE ( 'sthuf_soilt', 't_soil_soilt', 'sthu_irr_soilt' )
      CALL file_write_var(FILE, var_ids(i),                                   &
                          global_data_3d(:,1:nsoilt,1:sm_levels))

    CASE ( 'n_inorg')
      CALL file_write_var(FILE, var_ids(i),                                   &
                          global_data_2d(:,1:dim_cslayer))

    CASE ( 'n_inorg_soilt')
      CALL file_write_var(FILE, var_ids(i),                                   &
                          global_data_3d(:,1:nsoilt,1:dim_cslayer))

    CASE ( 'substr_ch4','mic_ch4','mic_act_ch4','acclim_ch4' )
      CALL file_write_var(FILE, var_ids(i),                                   &
                          global_data_2d(:,1:dim_ch4layer))

    CASE ( 'tsoil_deep' )
      CALL file_write_var(FILE, var_ids(i), global_data_2d(:,1:ns_deep))

    CASE ( 'canopy', 'nsnow', 'rgrain', 'rho_snow', 'snow_tile',              &
           'snow_depth', 'snow_grnd', 'tstar_tile', 'tsurf_elev_surft' )
      CALL file_write_var(FILE, var_ids(i), global_data_2d(:,1:nsurft))

      ! Snow and soil C/N variables with 2 levels dimensions.
    CASE ( 'rgrainl', 'snow_ds', 'snow_ice', 'snow_liq', 'tsnow' )
      CALL file_write_var(FILE, var_ids(i),                                   &
                          global_data_3d(:,1:nsurft,1:nsmax))
    CASE ( 'cs','ns' )
      CALL file_write_var(FILE, var_ids(i),                                   &
                          global_data_3d(:,1:dim_cslayer,1:dim_cs1))

    CASE ( 'cs_soilt' )
      CALL file_write_var(FILE, var_ids(i),                                   &
                          global_data_4d(:,1:nsoilt,1:dim_cslayer,            &
                                         1:dim_cs1))
    CASE ( 'n_soil' )
      CALL file_write_var(FILE, var_ids(i),                                   &
                          global_data_3d(:,1:dim_cslayer,1:dim_soil_n_pool))

    CASE ( 'n_soil_soilt' )
      CALL file_write_var(FILE, var_ids(i),                                   &
                          global_data_4d(:,1:nsoilt,1:dim_cslayer,            &
                                         1:dim_soil_n_pool))

      ! Cases for IMOGEN variables
      ! Each task runs its own version of IMOGEN - we just write the master
      ! task's versions.
    CASE ( 'co2_ppmv' )
      CALL file_write_var(FILE, var_ids(i), co2_ppmv)

    CASE ( 'co2_change_ppmv' )
      CALL file_write_var(FILE, var_ids(i), co2_change_ppmv)

    CASE ( 'dtemp_o' )
      CALL file_write_var(FILE, var_ids(i), dtemp_o)

    CASE ( 'fa_ocean' )
      CALL file_write_var(FILE, var_ids(i), fa_ocean)

    CASE ( 'seed_rain' )
      CALL file_write_var(FILE, var_ids(i), REAL(progs%seed_rain))

    CASE ( 'ch4_ppbv' )
      CALL file_write_var(FILE, var_ids(i), ch4_ppbv)

      ! Cases for river routing variables
    CASE ( 'rivers_lat_rp' )
      CALL file_write_var(FILE, var_ids(i), rivers_lat_rp)

    CASE ( 'rivers_lon_rp' )
      CALL file_write_var(FILE, var_ids(i), rivers_lon_rp)

    CASE ( 'rivers_sto_rp' )
      CALL file_write_var(FILE, var_ids(i), rivers_sto_rp)

    CASE ( 'rfm_surfstore_rp' )
      CALL file_write_var(FILE, var_ids(i), rfm_surfstore_rp)

    CASE ( 'rfm_substore_rp' )
      CALL file_write_var(FILE, var_ids(i), rfm_substore_rp)

    CASE ( 'rfm_flowin_rp' )
      CALL file_write_var(FILE, var_ids(i), rfm_flowin_rp)

    CASE ( 'rfm_bflowin_rp' )
      CALL file_write_var(FILE, var_ids(i), rfm_bflowin_rp)

      !-----------------------------------------------------------------------
      ! Ancillary variables
      !-----------------------------------------------------------------------

      ! latlon ancil namelist
    CASE ( 'grid_area' )
      CALL file_write_var(FILE, var_ids(i), global_data_1d)

      !Frac ancil namelist
    CASE ( 'frac' )
      CALL file_write_var(FILE, var_ids(i), global_data_2d(:,1:ntype))

      ! Vegetation properties ancil namelist
    CASE ( 't_growth_gb' )
      CALL file_write_var(FILE, var_ids(i), global_data_1d)

      !Soil properties ancil namelist
    CASE ( 'b      ', 'sathh  ', 'satcon ', 'sm_sat ', 'sm_crit',             &
           'sm_wilt', 'hcap   ', 'hcon   ' )
      CALL file_write_var(FILE, var_ids(i), global_data_2d(:,1:sm_levels))

    CASE ( 'b_soilt', 'sathh_soilt', 'satcon_soilt', 'sm_sat_soilt',          &
           'sm_crit_soilt', 'sm_wilt_soilt', 'hcap_soilt', 'hcon_soilt' )
      CALL file_write_var(FILE, var_ids(i),                                   &
                          global_data_3d(:,1:nsoilt,1:sm_levels))

    CASE ( 'albsoil' )
      CALL file_write_var(FILE, var_ids(i), global_data_1d)

    CASE ( 'albsoil_soilt' )
      CALL file_write_var(FILE, var_ids(i), global_data_2d(:,1:nsoilt))

    CASE ( 'clay', 'soil_ph' )
      CALL file_write_var(FILE, var_ids(i),                                   &
           global_data_2d(:,1:dim_cslayer))

    CASE ( 'clay_soilt', 'soil_ph_soilt' )
      CALL file_write_var(FILE, var_ids(i),                                   &
           global_data_3d(:,1:nsoilt,1:dim_cslayer))

      !Topmodel ancil namelist
    CASE ( 'fexp   ', 'ti_mean', 'ti_sig ' )
      CALL file_write_var(FILE, var_ids(i), global_data_1d)

    CASE ( 'toppdm%fexp_soilt', 'toppdm%ti_mean_soilt', 'toppdm%ti_sig_soilt' )
      CALL file_write_var(FILE, var_ids(i), global_data_2d(:,1:nsoilt))

      !Agric ancil namelist
    CASE ( 'frac_agr' )
      CALL file_write_var(FILE, var_ids(i), global_data_1d)

      !Crop props ancillaries namelist
    CASE ( 'cropsowdate', 'cropttveg  ', 'cropttrep  ','croplatestharvdate')
      CALL file_write_var(FILE, var_ids(i), global_data_2d(:,1:ncpft))

      !Irrigation ancillaries namelist
    CASE ( 'frac_irrig' )
      CALL file_write_var(FILE, var_ids(i), global_data_1d)

    CASE ( 'irrfrac_irrtiles' )
      CALL file_write_var(FILE, var_ids(i), global_data_1d)

    CASE ( 'frac_irr_all_tiles' )
      CALL file_write_var(FILE, var_ids(i), frac_irrig_all_tiles_real)

    CASE ( 'set_irrfrac_on_irrtiles' )
      CALL file_write_var(FILE, var_ids(i), set_irrfrac_irrtiles_real)

    CASE ( 'irrtiles' )
      CALL file_write_var(FILE, var_ids(i), irrtiles_real)

    CASE ( 'nirrtile' )
      CALL file_write_var(FILE, var_ids(i), nirrtile_real)

      ! Water resources properties ancil namelist
    CASE ( 'conveyance_loss' )
      CALL file_write_var(FILE, var_ids(i), global_data_1d)

    CASE ( 'irrig_eff' )
      CALL file_write_var(FILE, var_ids(i), global_data_1d)

    CASE ( 'sfc_water_frac' )
      CALL file_write_var(FILE, var_ids(i), global_data_1d)

      !CO2 ancil namelist
    CASE ( 'co2_mmr' )
      CALL file_write_var(FILE, var_ids(i), co2_mmr)

    CASE ( 'lake_depth_gb' )
      CALL file_write_var(FILE, var_ids(i), global_data_1d)
    CASE DEFAULT
      CALL log_fatal("write_dump",                                            &
                     "Unrecognised variable for dump - " //                   &
                     TRIM(identifiers(i)))
    END SELECT
  END IF  ! MASTER TASK
END DO

! We are done with the file and dictionaries
IF ( is_master_task() ) CALL file_close(FILE)

IF ( ALLOCATED(global_data_1d) ) DEALLOCATE(global_data_1d)
IF ( ALLOCATED(global_data_2d) ) DEALLOCATE(global_data_2d)
IF ( ALLOCATED(global_data_3d) ) DEALLOCATE(global_data_3d)
IF ( ALLOCATED(global_data_4d) ) DEALLOCATE(global_data_4d)

RETURN

END SUBROUTINE write_dump
#endif
