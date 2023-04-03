! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE dump_mod

USE io_constants, ONLY: format_len, format_ascii, format_ncdf,                &
                         max_sdf_name_len

USE logging_mod, ONLY: log_info, log_fatal

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------------

INTEGER, PARAMETER :: max_dim_dump = 14
INTEGER, PARAMETER :: max_var_dump = 84

CHARACTER(LEN=max_sdf_name_len), PARAMETER ::                                 &
  land_dim_name    = "land",                                                  &
  pft_dim_name     = "pft",                                                   &
  cpft_dim_name    = "cpft",                                                  &
  sc_pool_dim_name = "scpool",                                                &
  soil_n_pool_dim_name = "snpool",                                            &
  sc_layer_dim_name = "sclayer",                                              &
  ch4layer_dim_name  = "ch4layer",                                            &
  snow_dim_name    = "snow",                                                  &
  soil_dim_name    = "soil",                                                  &
  tile_dim_name    = "tile",                                                  &
  type_dim_name    = "type",                                                  &
  soilt_dim_name   = "soilt",                                                 &
  scalar_dim_name  = "scalar",                                                &
  nolevs_dim_name  = "olevs",                                                 &
  nfarray_dim_name = "nfarray",                                               &
  seed_dim_name    = "seed",                                                  &
  bedrock_dim_name = "bedrock",                                               &
  p_rivers_dim_name = "p_rivers"

CHARACTER(LEN=format_len), PARAMETER :: dump_format = format_ncdf

  !Create a defined type to store the flags for which ancils are read
  !from the dump file
TYPE ancil_flags
  LOGICAL :: latlon
  LOGICAL :: frac
  LOGICAL :: vegetation_props
  LOGICAL :: soil_props
  LOGICAL :: top
  LOGICAL :: agric
  LOGICAL :: crop_props
  LOGICAL :: irrig
  LOGICAL :: rivers_props
  LOGICAL :: water_resources_props
  LOGICAL :: co2
  LOGICAL :: flake
END TYPE

TYPE(ancil_flags) :: ancil_dump_read

!-----------------------------------------------------------------------------
! Visibility declarations
!-----------------------------------------------------------------------------
PRIVATE
PUBLIC  max_var_dump, required_vars_for_configuration,                        &
        read_dump, write_dump,                                                &
        ancil_dump_read

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


SUBROUTINE required_vars_for_configuration(nvars, identifiers,                &
                                           nvars_from_ancil, vars_from_ancil, &
                                           l_output_mode,                     &
                                           total_snow, include_imogen_vars,   &
                                           read_or_write_dump)

USE update_mod, ONLY: l_imogen, have_prescribed_sthuf

USE jules_hydrology_mod, ONLY: l_top

USE jules_radiation_mod, ONLY: l_snow_albedo, l_embedded_snow

USE jules_soil_biogeochem_mod, ONLY: soil_model_ecosse, soil_model_rothc,     &
                                     soil_bgc_model, l_ch4_microbe

USE jules_soil_ecosse_mod, ONLY: l_soil_N

USE jules_vegetation_mod, ONLY: can_model, can_rad_mod,                       &
                                l_triffid, l_phenol, l_veg_compete,           &
                                l_crop, l_landuse,                            &
                                l_nitrogen, l_prescsow, l_trif_crop,          &
                                photo_acclim_model, photo_adapt,              &
                                l_croprotate

USE jules_water_resources_mod, ONLY: l_water_irrigation, l_water_resources

USE jules_irrig_mod, ONLY: l_irrig_dmd

USE jules_snow_mod, ONLY: nsmax

USE jules_surface_mod, ONLY: l_elev_land_ice, l_flake_model

USE jules_soil_mod, ONLY: l_bedrock, sm_levels, l_tile_soil

USE fire_mod,     ONLY: fire_cntl, l_fire

USE metstats_mod, ONLY: metstats_flag, l_metstats

USE update_mod, ONLY: l_daily_disagg, precip_disagg_method,                   &
                      prescribed_sthuf_levels

USE jules_rivers_mod, ONLY: l_rivers, i_river_vn, rivers_trip, rivers_rfm

USE logging_mod, ONLY: log_warn

USE imogen_run, ONLY: land_feed_ch4

USE prognostics, ONLY:                                                        &
  l_broadcast_soilt


IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Returns the identifiers of the prognostic variables for the current
!   configuration of the model
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
INTEGER, INTENT(OUT) :: nvars  ! The number of variables
CHARACTER(LEN=*), INTENT(OUT) :: identifiers(:)
                               ! The model identifiers of the required
                               ! variables

INTEGER, INTENT(OUT)          :: nvars_from_ancil  !Number of ancil variables
CHARACTER(LEN=*), INTENT(OUT) :: vars_from_ancil(:)
                               ! The model identifiers of the required
                               ! ancil variables

LOGICAL, INTENT(IN)           :: l_output_mode
                               ! True if the intention of the call to this
                               ! routine is to get the variables to output
                               ! False if for input
                               ! This enables the broadcasting of non-soil
                               ! tiled initial conditions for soil tiling
                               ! runs

LOGICAL, INTENT(IN), OPTIONAL :: total_snow
                               ! T - return only the variables required if
                               !     snow layer values are to be set by the
                               !     model
                               ! F - return all snow variables required
                               !     by the current configuration
LOGICAL, INTENT(IN), OPTIONAL :: include_imogen_vars
                               ! T - include IMOGEN prognostics in the
                               !     returned list if they are needed
                               ! F - do not include IMOGEN prognostics in
                               !     the returned list, even if they are
                               !     needed
LOGICAL, INTENT(IN), OPTIONAL :: read_or_write_dump
                               ! T - include all variables that are needed
                               !     when reading or writing to dump files
                               ! F - don't include variables that are only
                               !     inputted or outputted in dump files

!-----------------------------------------------------------------------------
! Local scalar parameters.
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER ::                                                &
  RoutineName = 'REQUIRED_VARS_FOR_CONFIGURATION'
    ! The name of this routine.

! Work variables
LOGICAL :: total_snow_local  ! Local version of total_snow
                             ! Defaults to FALSE if not present (i.e.
                             ! return all required snow vars)
LOGICAL :: inc_imogen_local  ! Local version of include_imogen_vars
                             ! Defaults to TRUE if not present, i.e.
                             ! return IMOGEN vars if they are required
LOGICAL :: read_or_write_dump_local
                             ! Local version of read_or_write_dump
                             ! Defaults to TRUE if not present
LOGICAL :: l_append_soilt
                             ! Controls whether vars are appended with _soilt
INTEGER :: num_prescribed_sthuf_levels
                             ! Number of levels in the prescibed sthuf field


!-----------------------------------------------------------------------------

! Check for presence of optional arguments
total_snow_local          = .FALSE.
inc_imogen_local          = .TRUE.
read_or_write_dump_local  = .TRUE.

! Decide whether we need to append _soilt to variables
l_append_soilt = .FALSE.
IF (l_tile_soil) THEN
  IF (l_output_mode) THEN
    l_append_soilt = .TRUE.
  ELSE
    IF ( .NOT. l_broadcast_soilt) THEN
      l_append_soilt = .TRUE.
    END IF
  END IF
END IF

IF ( PRESENT(total_snow) )                                                    &
  total_snow_local         = total_snow
IF ( PRESENT(include_imogen_vars) )                                           &
  inc_imogen_local         = include_imogen_vars
IF ( PRESENT(read_or_write_dump) )                                            &
  read_or_write_dump_local = read_or_write_dump

! First, set up the array with variables that are required with every
! configuration.
! Subroutine add_to_list is called with the l_append_soilt for a variable that
! supports soil tiles.
nvars = 0

! Variables without soil tiling.
CALL add_to_list( 'canopy',     nvars, identifiers )
! Note that for total_snow = T, snow_surft is the only required snow variable
! (in addition to rgrain, depending on the snow albedo switches).
CALL add_to_list( 'snow_tile',  nvars, identifiers )
CALL add_to_list( 'tstar_tile', nvars, identifiers )

! Variables with soil tiling.
CALL add_to_list( 'cs',         nvars, identifiers, l_append_soilt )
CALL add_to_list( 't_soil',     nvars, identifiers, l_append_soilt )

!-----------------------------------------------------------------------------
! Variables that are needed only if certain options set.
!-----------------------------------------------------------------------------
IF ( can_rad_mod == 1 ) THEN
  CALL add_to_list( 'gs', nvars, identifiers )
END IF

IF (ALLOCATED(prescribed_sthuf_levels)) THEN
  num_prescribed_sthuf_levels = SIZE(prescribed_sthuf_levels)
ELSE
  num_prescribed_sthuf_levels = 0
END IF

IF ( ( .NOT. have_prescribed_sthuf ) .OR.                                     &
     ( num_prescribed_sthuf_levels < sm_levels) ) THEN
  CALL add_to_list( 'sthuf', nvars, identifiers, l_append_soilt )
END IF

IF ( l_phenol ) THEN
  ! With phenology on, LAI is prognostic for all PFTs
  CALL add_to_list( 'lai', nvars, identifiers )
ELSE IF ( l_crop ) THEN
  ! Otherwise, if the crop scheme is on, LAI is prognostic for crop PFTs only
  CALL add_to_list( 'croplai', nvars, identifiers )
END IF

IF ( l_triffid ) THEN
  ! With TRIFFID on, canopy height is prognostic for all PFTs
  CALL add_to_list( 'canht', nvars, identifiers )
ELSE IF ( l_crop ) THEN
  ! Otherwise, if the crop scheme is on, canopy height is prognostic for crop
  ! PFTs only
  CALL add_to_list( 'cropcanht', nvars, identifiers )
END IF

IF ( l_triffid .AND. l_landuse ) THEN
  ! With TRIFFID on,  agricultural frac and WP pools are required for LUC
  CALL add_to_list( 'frac_agr_prev',  nvars, identifiers )
  CALL add_to_list( 'wood_prod_fast', nvars, identifiers )
  CALL add_to_list( 'wood_prod_med',  nvars, identifiers )
  CALL add_to_list( 'wood_prod_slow', nvars, identifiers )
  IF (l_trif_crop) THEN
    CALL add_to_list( 'frac_past_prev', nvars, identifiers )
  END IF
END IF

IF ( soil_bgc_model == soil_model_rothc .AND. l_nitrogen ) THEN
  CALL add_to_list( 'n_inorg', nvars, identifiers, l_append_soilt )
  CALL add_to_list( 'ns',      nvars, identifiers, l_append_soilt )
END IF

IF ( soil_bgc_model == soil_model_ecosse .AND. l_soil_N ) THEN
  CALL add_to_list( 'n_soil', nvars, identifiers, l_append_soilt )
END IF

IF ( l_ch4_microbe ) THEN
  CALL add_to_list( 'substr_ch4', nvars, identifiers, l_append_soilt )
  CALL add_to_list( 'mic_ch4', nvars, identifiers, l_append_soilt )
  CALL add_to_list( 'mic_act_ch4', nvars, identifiers, l_append_soilt )
  CALL add_to_list( 'acclim_ch4', nvars, identifiers, l_append_soilt )
END IF

! TOPMODEL variables.
IF ( l_top ) THEN
  ! Wetness in deep layer and depth to water table..
  CALL add_to_list( 'sthzw', nvars, identifiers, l_append_soilt )
  CALL add_to_list( 'zw',    nvars, identifiers, l_append_soilt )
END IF

! Additional crop scheme prognostics, required if crop scheme is on
IF ( l_crop ) THEN
  CALL add_to_list( 'cropdvi',      nvars, identifiers )
  CALL add_to_list( 'croprootc',    nvars, identifiers )
  CALL add_to_list( 'cropharvc',    nvars, identifiers )
  CALL add_to_list( 'cropreservec', nvars, identifiers )
END IF

IF ( l_irrig_dmd ) THEN
  CALL add_to_list( 'sthu_irr', nvars, identifiers, l_append_soilt )
END IF

! Additional deep soil temperature if bedrock is on
IF ( l_bedrock ) THEN
  CALL add_to_list( 'tsoil_deep', nvars, identifiers )
END IF

! River storage if river routing on and using TRIP
IF ( l_rivers .AND. ( i_river_vn == rivers_trip ) ) THEN
  IF ( read_or_write_dump_local ) THEN
    CALL add_to_list( 'rivers_sto_rp', nvars, identifiers )
  ELSE
    CALL log_warn( RoutineName,                                               &
                  "rivers_sto_rp will be initialised to zero.")
  END IF
END IF

! Surface and subsurface stores and flows if river routing on and using RFM
IF ( l_rivers .AND. (i_river_vn == rivers_rfm) ) THEN
  IF ( read_or_write_dump_local ) THEN
    CALL add_to_list( 'rfm_surfstore_rp', nvars, identifiers )
    CALL add_to_list( 'rfm_substore_rp',  nvars, identifiers )
    CALL add_to_list( 'rfm_flowin_rp',    nvars, identifiers )
    CALL add_to_list( 'rfm_bflowin_rp',   nvars, identifiers )
  ELSE
    CALL log_warn( RoutineName,                                               &
                  "RFM river prognostics will be initialised to zero.")
  END IF
END IF

! FLake prognostics
IF (l_flake_model) THEN
  nvars = nvars + 1
  identifiers(nvars) = 'lake_fetch_gb'

  nvars = nvars + 1
  identifiers(nvars) = 'lake_t_mean_gb'

  nvars = nvars + 1
  identifiers(nvars) = 'lake_t_mxl_gb'

  nvars = nvars + 1
  identifiers(nvars) = 'lake_h_mxl_gb'

  nvars = nvars + 1
  identifiers(nvars) = 'lake_t_ice_gb'

  nvars = nvars + 1
  identifiers(nvars) = 'lake_h_ice_gb'

  nvars = nvars + 1
  identifiers(nvars) = 'lake_shape_factor_gb'
END IF

!-----------------------------------------------------------------------------
! Add metstats prognostic variables if switched on
!-----------------------------------------------------------------------------
IF ( l_metstats ) THEN
  IF (metstats_flag%temp_max_00h) THEN
    CALL add_to_list( 'temp_max_00h',   nvars, identifiers )
    CALL add_to_list( 'temp_max_00h_r', nvars, identifiers )
  END IF

  IF (metstats_flag%temp_ave_00h) THEN
    CALL add_to_list( 'temp_ave_00h',   nvars, identifiers )
    CALL add_to_list( 'temp_ave_00h_r', nvars, identifiers )
  END IF

  IF (metstats_flag%temp_ave_nday) THEN
    CALL add_to_list( 'temp_ave_nday', nvars, identifiers )
  END IF

  IF (metstats_flag%temp_pnt_12h) THEN
    CALL add_to_list( 'temp_pnt_12h', nvars, identifiers )
  END IF

  IF (metstats_flag%prec_tot_00h) THEN
    CALL add_to_list( 'prec_tot_00h',   nvars, identifiers )
    CALL add_to_list( 'prec_tot_00h_r', nvars, identifiers )
  END IF

  IF (metstats_flag%prec_tot_12h) THEN
    CALL add_to_list( 'prec_tot_12h',   nvars, identifiers )
    CALL add_to_list( 'prec_tot_12h_r', nvars, identifiers )
  END IF

  IF (metstats_flag%rhum_min_00h) THEN
    CALL add_to_list( 'rhum_min_00h',   nvars, identifiers )
    CALL add_to_list( 'rhum_min_00h_r', nvars, identifiers )
  END IF

  IF (metstats_flag%rhum_pnt_12h) THEN
    CALL add_to_list( 'rhum_pnt_12h', nvars, identifiers )
  END IF

  IF (metstats_flag%dewp_ave_00h) THEN
    CALL add_to_list( 'dewp_ave_00h',   nvars, identifiers )
    CALL add_to_list( 'dewp_ave_00h_r', nvars, identifiers )
  END IF

  IF (metstats_flag%wind_ave_00h) THEN
    CALL add_to_list( 'wind_ave_00h',   nvars, identifiers )
    CALL add_to_list( 'wind_ave_00h_r', nvars, identifiers )
  END IF

  IF (metstats_flag%wind_pnt_12h) THEN
    CALL add_to_list( 'wind_pnt_12h', nvars, identifiers )
  END IF
END IF

!-----------------------------------------------------------------------------
! Add fire prognostic variables if switched on
!-----------------------------------------------------------------------------
IF ( l_fire ) THEN
  IF (fire_cntl%mcarthur%flag) THEN
    CALL add_to_list( 'fire_mcarthur_r_dr', nvars, identifiers )
    CALL add_to_list( 'fire_mcarthur_n_dr', nvars, identifiers )
  END IF

  IF (fire_cntl%canadian%flag) THEN
    CALL add_to_list( 'fire_canadian_ffmc',      nvars, identifiers )
    CALL add_to_list( 'fire_canadian_ffmc_mois', nvars, identifiers )
    CALL add_to_list( 'fire_canadian_dmc',       nvars, identifiers )
    CALL add_to_list( 'fire_canadian_dc',        nvars, identifiers )
  END IF

  IF (fire_cntl%nesterov%flag) THEN
    CALL add_to_list( 'fire_nesterov', nvars, identifiers )
  END IF
END IF !Fire

!-----------------------------------------------------------------------------
! Add IMOGEN variables if IMOGEN is on and they are requested
!-----------------------------------------------------------------------------
IF ( l_imogen .AND. inc_imogen_local ) THEN
  CALL add_to_list( 'co2_ppmv',        nvars, identifiers )
  CALL add_to_list( 'co2_change_ppmv', nvars, identifiers )
  CALL add_to_list( 'dtemp_o',         nvars, identifiers )
  CALL add_to_list( 'fa_ocean',        nvars, identifiers )
  CALL add_to_list( 'seed_rain',       nvars, identifiers )

  ! This should possibly be under an l_triffid switch, but is here for now
  ! so as not to change behaviour for non-TRIFFID runs
  CALL add_to_list( 'cv', nvars, identifiers )

  IF (land_feed_ch4) THEN
    CALL add_to_list( 'ch4_ppbv', nvars, identifiers )
  END IF
END IF

IF ( l_daily_disagg .AND. read_or_write_dump_local .AND.                      &
     (precip_disagg_method > 1) ) THEN
  CALL add_to_list( 'seed_rain', nvars, identifiers )
END IF

!-----------------------------------------------------------------------------
! Work out what snow variables are required
! If total_snow = T, only snow_surft is required, and is always required
! so if already in the list
! We just need to add the other required snow variables based on the scheme
! enabled if total_snow = F
!-----------------------------------------------------------------------------
IF ( .NOT. total_snow_local ) THEN
  ! Snow variables not specifically for the multi-layer model
  CALL add_to_list( 'rho_snow',   nvars, identifiers )
  CALL add_to_list( 'snow_depth', nvars, identifiers )

  IF ( can_model == 4 ) THEN
    CALL add_to_list( 'snow_grnd', nvars, identifiers )
  END IF

  ! Variables for the multi-layer snow model.
  IF ( nsmax > 0 ) THEN
    CALL add_to_list( 'nsnow',    nvars, identifiers )
    CALL add_to_list( 'snow_ds',  nvars, identifiers )
    CALL add_to_list( 'snow_ice', nvars, identifiers )
    CALL add_to_list( 'snow_liq', nvars, identifiers )
    CALL add_to_list( 'tsnow',    nvars, identifiers )

    IF ( l_snow_albedo .OR. l_embedded_snow ) THEN
      CALL add_to_list( 'rgrainl', nvars, identifiers )
    END IF
  END IF   ! nsmax

END IF  ! total_snow

! Snow surface grain size.
! Note this is required in addition to the layered variable rgrainl if
! using layered snow, because the surface grain size is used for radiative
! calculations. (But as rgrain is copied from rgrainl, possibly we
! should only require rgrainl for the layered model?)
IF ( l_snow_albedo .OR. l_embedded_snow ) THEN
  CALL add_to_list( 'rgrain', nvars, identifiers )
END IF

IF ( l_elev_land_ice) THEN
  CALL add_to_list( 'tsurf_elev_surft', nvars, identifiers )
END IF

!-----------------------------------------------------------------------------
! Work out which ancillary fields need to be in the dump file
! The logic here is slightly different from normal prognostic variables
!
! There are several cases to consider:
! -The var is not needed
! -The var should be populated by ancil and written to the dump (ie write only
!  as far as this module is concerned)
! -The var should be both read from the initial condition namelist and written
!  to the dump
!
! To get the information about which ancils will be read from file, we use
! ancil_dump_read%xyz:  Defaulting to false, set in init_ancillaries_mod
!-----------------------------------------------------------------------------

!Initialise
nvars_from_ancil = 0

!-----------------------------------------------------------------------------
! Set up read/write details for latlon ancil.
! Note that latitude and longitude are not handled here.
!-----------------------------------------------------------------------------
IF ( l_water_irrigation ) THEN
  CALL add_to_list( 'grid_area', nvars, identifiers )

  ! If not reading from dump, indicate variables come from an ancil file.
  IF ( .NOT. ancil_dump_read%latlon) THEN
    CALL add_to_list( 'grid_area', nvars_from_ancil, vars_from_ancil,         &
                                     l_append_soilt )
  END IF
END IF

!-----------------------------------------------------------------------------
!Set up read/write details for frac ancil
!This is a special case: with competing vegetation on, frac is prognostic

! frac should always be written to the dump
CALL add_to_list( 'frac', nvars, identifiers )

IF ( l_veg_compete ) THEN
  ! With competing veg on, frac is a prognostic and should always be read
  ! as an initial condition
  ancil_dump_read%frac  = .TRUE.
END IF

! If not reading from dump, indicate variables come from an ancil file.
IF ( .NOT. ancil_dump_read%frac ) THEN
  CALL add_to_list( 'frac', nvars_from_ancil, vars_from_ancil )
END IF

!-----------------------------------------------------------------------------
! Set up read/write details for vegetation ancil.

IF ( photo_acclim_model == photo_adapt ) THEN
  CALL add_to_list( 't_growth_gb', nvars, identifiers )

  ! If not reading from dump, indicate variables come from an ancil file.
  IF ( .NOT. ancil_dump_read%vegetation_props) THEN
    CALL add_to_list( 't_growth_gb', nvars_from_ancil, vars_from_ancil,       &
                                     l_append_soilt )
  END IF

END IF ! photo_acclim_model

!---------------------------------------------------------------------------
! Set up read/write details for soil properties ancils
! Most soil properties are always needed.
CALL add_to_list( 'b',       nvars, identifiers, l_append_soilt )
CALL add_to_list( 'sathh',   nvars, identifiers, l_append_soilt )
CALL add_to_list( 'satcon',  nvars, identifiers, l_append_soilt )
CALL add_to_list( 'sm_sat',  nvars, identifiers, l_append_soilt )
CALL add_to_list( 'sm_crit', nvars, identifiers, l_append_soilt )
CALL add_to_list( 'sm_wilt', nvars, identifiers, l_append_soilt )
CALL add_to_list( 'hcap',    nvars, identifiers, l_append_soilt )
CALL add_to_list( 'hcon',    nvars, identifiers, l_append_soilt )
CALL add_to_list( 'albsoil', nvars, identifiers, l_append_soilt )

IF ( soil_bgc_model == soil_model_rothc .OR.                                  &
     soil_bgc_model == soil_model_ecosse) THEN
  CALL add_to_list( 'clay', nvars, identifiers, l_append_soilt )
END IF

IF ( soil_bgc_model == soil_model_ecosse ) THEN
  CALL add_to_list( 'soil_ph', nvars, identifiers, l_append_soilt )
END IF

! If not reading from dump, indicate variables come from an ancil file.
IF ( .NOT. ancil_dump_read%soil_props ) THEN
  CALL add_to_list( 'b',       nvars_from_ancil, vars_from_ancil,             &
                               l_append_soilt )
  CALL add_to_list( 'sathh',   nvars_from_ancil, vars_from_ancil,             &
                               l_append_soilt )
  CALL add_to_list( 'satcon',  nvars_from_ancil, vars_from_ancil,             &
                               l_append_soilt )
  CALL add_to_list( 'sm_sat',  nvars_from_ancil, vars_from_ancil,             &
                               l_append_soilt )
  CALL add_to_list( 'sm_crit', nvars_from_ancil, vars_from_ancil,             &
                               l_append_soilt )
  CALL add_to_list( 'sm_wilt', nvars_from_ancil, vars_from_ancil,             &
                               l_append_soilt )
  CALL add_to_list( 'hcap',    nvars_from_ancil, vars_from_ancil,             &
                               l_append_soilt )
  CALL add_to_list( 'hcon',    nvars_from_ancil, vars_from_ancil,             &
                               l_append_soilt )
  CALL add_to_list( 'albsoil', nvars_from_ancil, vars_from_ancil,             &
                               l_append_soilt )
  IF ( soil_bgc_model == soil_model_rothc .OR.                                &
       soil_bgc_model == soil_model_ecosse ) THEN
    CALL add_to_list( 'clay',  nvars_from_ancil, vars_from_ancil,             &
                               l_append_soilt )
  END IF
  IF ( soil_bgc_model == soil_model_ecosse ) THEN
    CALL add_to_list( 'soil_ph', nvars_from_ancil,                            &
                                 vars_from_ancil, l_append_soilt )
  END IF
END IF

!----------------------------------------------------------------------------
!Set up read/write details for topmodel ancils.
!Depends on the value of l_top

IF ( l_top ) THEN
  CALL add_to_list( 'fexp',    nvars, identifiers, l_append_soilt )
  CALL add_to_list( 'ti_mean', nvars, identifiers, l_append_soilt )
  CALL add_to_list( 'ti_sig',  nvars, identifiers, l_append_soilt )

  ! If not reading from dump, indicate variables come from an ancil file.
  IF ( .NOT. ancil_dump_read%top) THEN
    CALL add_to_list( 'fexp',    nvars_from_ancil,vars_from_ancil,            &
                                 l_append_soilt )
    CALL add_to_list( 'ti_mean', nvars_from_ancil,vars_from_ancil,            &
                                 l_append_soilt )
    CALL add_to_list( 'ti_sig',  nvars_from_ancil,vars_from_ancil,            &
                                 l_append_soilt )
  END IF
END IF !l_top

!----------------------------------------------------------------------------
!Set up read/write details for agric ancil namelist
!This is always needed

CALL add_to_list( 'frac_agr', nvars, identifiers )

! If not reading from dump, indicate variables come from an ancil file.
IF ( .NOT. ancil_dump_read%agric ) THEN
  CALL add_to_list( 'frac_agr', nvars_from_ancil, vars_from_ancil )
END IF

!---------------------------------------------------------------------------
!Set up read/write details for crop properites ancil namelist
!The vars required depend on a couple of switches

IF ( l_crop ) THEN
  CALL add_to_list( 'cropttveg', nvars, identifiers )
  CALL add_to_list( 'cropttrep', nvars, identifiers )
  IF ( l_prescsow ) THEN  !we additionally need the crop sowing date
    CALL add_to_list( 'cropsowdate', nvars, identifiers )
    !When l_croprotate is TRUE we need the latest possible harvest date but
    !it is also possible to use croplatestharvdate when l_croprotate is FALSE
    CALL add_to_list( 'croplatestharvdate', nvars, identifiers )
  END IF

  ! If not reading from dump, indicate variables come from an ancil file.
  IF ( .NOT. ancil_dump_read%crop_props ) THEN
    CALL add_to_list( 'cropttveg', nvars_from_ancil, vars_from_ancil )
    CALL add_to_list( 'cropttrep', nvars_from_ancil, vars_from_ancil )
    IF ( l_prescsow ) THEN  !we additionally need the crop sowing date
      CALL add_to_list( 'cropsowdate', nvars_from_ancil, vars_from_ancil )
      !When l_croprotate is TRUE we need the latest possible harvest date but
      !it is also possible to use croplatestharvdate when l_croprotate is
      !FALSE
      CALL add_to_list( 'croplatestharvdate', nvars_from_ancil,               &
                                              vars_from_ancil )
    END IF
  END IF

END IF !l_crop

  !---------------------------------------------------------------------------
  !Set up read/write details for irrigation ancil namelist
  !The vars required depend on l_irrig_dmd
IF ( l_irrig_dmd ) THEN
  CALL add_to_list( 'frac_irrig',         nvars, identifiers )
  CALL add_to_list( 'frac_irr_all_tiles', nvars, identifiers )
  !These only matter if frac_irrig_all_tiles is false, but can't gurantee that
  !that its value will be known when reading that variable from the dump
  CALL add_to_list( 'irrtiles',           nvars, identifiers )
  CALL add_to_list( 'nirrtile',           nvars, identifiers )
  CALL add_to_list( 'set_irrfrac_on_irrtiles', nvars, identifiers )
  CALL add_to_list( 'irrfrac_irrtiles',   nvars, identifiers )

  ! If not reading from dump, indicate variables come from an ancil file.
  IF ( .NOT. ancil_dump_read%irrig ) THEN
    CALL add_to_list( 'frac_irrig', nvars_from_ancil, vars_from_ancil )
    CALL add_to_list( 'frac_irr_all_tiles', nvars_from_ancil,                 &
                                            vars_from_ancil )
    !These only matter if frac_irrig_all_tiles is false, but can't gurantee that
    !that its value will be known when reading that variable from the dump
    CALL add_to_list( 'irrtiles',  nvars_from_ancil, vars_from_ancil )
    CALL add_to_list( 'nirrtile',  nvars_from_ancil, vars_from_ancil )
    CALL add_to_list( 'set_irrfrac_on_irrtiles', nvars_from_ancil,            &
                                                 vars_from_ancil )
    CALL add_to_list( 'irrfrac_irrtiles',        nvars_from_ancil,            &
                                                 vars_from_ancil )
  END IF

END IF

!Rivers need further work before being included
!IF ( ancil_dump_read%rivers_props ) THEN
!END IF

!---------------------------------------------------------------------------
!Set up read/write details for water resources properites ancil namelist.

IF ( l_water_resources ) THEN

  ! Variables that are always required.
  CALL add_to_list( 'conveyance_loss', nvars, identifiers )
  CALL add_to_list( 'sfc_water_frac', nvars, identifiers )

  ! If not reading from dump, indicate variables come from an ancil file.
  IF ( .NOT. ancil_dump_read%water_resources_props ) THEN
    CALL add_to_list( 'conveyance_loss', nvars_from_ancil, vars_from_ancil )
    CALL add_to_list( 'sfc_water_frac',  nvars_from_ancil, vars_from_ancil )
  END IF

  ! Variables that depend on a switch.
  ! Note that grid_area is handled in the latlon section above.
  IF ( l_water_irrigation ) THEN
    CALL add_to_list( 'irrig_eff', nvars, identifiers )
    ! If not reading from dump, indicate variables come from an ancil file.
    IF ( .NOT. ancil_dump_read%water_resources_props ) THEN
      CALL add_to_list( 'irrig_eff',  nvars_from_ancil, vars_from_ancil )
    END IF
  END IF

END IF  !  l_water_resources

!---------------------------------------------------------------------------
!Set up read/write details for co2 ancil namelist
!Always required, even if it just carries around the default value set in
!the aero module
CALL add_to_list( 'co2_mmr', nvars, identifiers )

! If not reading from dump, indicate variable comes from another source.
IF ( .NOT. ancil_dump_read%co2 ) THEN
  CALL add_to_list( 'co2_mmr', nvars_from_ancil, vars_from_ancil )
END IF

  !-----------------------------------------------------------------------------
  !Set up read/write details for flake ancil namelist
  !Required if l_flake_model = True, even if it just carries around the default value

IF ( l_flake_model ) THEN
  nvars = nvars + 1
  identifiers(nvars) = 'lake_depth_gb'

  IF ( .NOT. ancil_dump_read%flake ) THEN
    nvars_from_ancil = nvars_from_ancil + 1
    vars_from_ancil(nvars_from_ancil) = 'lake_depth_gb'
  END IF

END IF !l_flake_model

RETURN

END SUBROUTINE required_vars_for_configuration

!#############################################################################

SUBROUTINE add_to_list( identifier, nvars, varlist, l_append )

! Adds an identifier to a list, optionally appending a suffix.

USE logging_mod, ONLY: log_fatal

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with intent(in)
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), INTENT(IN) ::                                               &
  identifier
    ! An identifier for a variable.

!-----------------------------------------------------------------------------
! Arguments with intent(inout)
!-----------------------------------------------------------------------------

INTEGER, INTENT(INOUT) ::                                                     &
  nvars
    ! The number of values set in the list.

CHARACTER(LEN=*), INTENT(INOUT) ::                                            &
  varlist(:)
    ! A list of identifiers.

!-----------------------------------------------------------------------------
! Optional arguments with intent(in)
!-----------------------------------------------------------------------------
LOGICAL, INTENT(IN), OPTIONAL ::                                              &
  l_append
    ! Flag indicating if a suffix should be attached to the identifier.

!-----------------------------------------------------------------------------
! Local scalar parameters.
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER ::                                                &
  RoutineName = 'ADD_TO_LIST'
    ! The name of this routine.

!-----------------------------------------------------------------------------
! Local scalar variables.
!-----------------------------------------------------------------------------
LOGICAL ::                                                                    &
  l_append_soilt_local
    ! A local version of the optional argument l_append.

! end of header
!-----------------------------------------------------------------------------

! Set a local variables depending on the optional argument.
IF ( PRESENT(l_append) ) THEN
  l_append_soilt_local = l_append
ELSE
  l_append_soilt_local = .FALSE.
END IF

!-----------------------------------------------------------------------------
! Append to the list.
!-----------------------------------------------------------------------------
IF ( nvars + 1 > SIZE( varlist ) ) THEN
  CALL log_fatal( RoutineName,                                                &
                  "Too many values. Increase size of variable." )
  STOP
ELSE
  nvars          = nvars + 1
  varlist(nvars) = identifier
  IF ( l_append_soilt_local ) THEN
    varlist(nvars) = TRIM(varlist(nvars)) // '_soilt'
  END IF
END IF

RETURN
END SUBROUTINE add_to_list
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

SUBROUTINE read_dump(file_name, identifiers)

USE jules_fields_mod, ONLY: crop_vars, psparms, toppdm, ainfo, trif_vars,     &
                            soilecosse, progs, progs_data, trifctltype,       &
                            jules_vars


USE model_grid_mod, ONLY:                                                     &
  global_land_pts, grid_area_ij

USE ancil_info, ONLY:                                                         &
  land_pts, dim_cs1, dim_soil_n_pool, nsurft, nsoilt, dim_cslayer, nsoilt

USE fluxes, ONLY:                                                             &
  t_growth_gb

USE jules_surface_types_mod, ONLY:                                            &
  npft, ntype, ncpft

USE imogen_constants, ONLY:                                                   &
  n_olevs, nfarray

USE prognostics, ONLY:                                                        &
  l_broadcast_soilt


USE jules_irrig_mod, ONLY: irrtiles, frac_irrig_all_tiles, nirrtile,          &
  set_irrfrac_on_irrtiles

USE imogen_progs, ONLY:                                                       &
  co2_ppmv, co2_change_ppmv, dtemp_o, fa_ocean, ch4_ppbv

USE jules_snow_mod, ONLY:                                                     &
  nsmax

USE jules_soil_mod, ONLY:                                                     &
  sm_levels, ns_deep, l_tile_soil

USE jules_soil_biogeochem_mod, ONLY:                                          &
  dim_ch4layer

USE fire_mod, ONLY:                                                           &
  fire_prog

USE metstats_mod, ONLY:                                                       &
  metstats_prog

USE jules_rivers_mod, ONLY:                                                   &
  rivers_sto_rp, rfm_surfstore_rp,                                            &
  rfm_substore_rp, rfm_flowin_rp, rfm_bflowin_rp

USE jules_water_resources_mod, ONLY:                                          &
  conveyance_loss, irrig_eff, sfc_water_frac

USE aero, ONLY:                                                               &
  co2_mmr

!Others
USE mpi, ONLY: mpi_real, mpi_comm_world, mpi_logical, mpi_integer
USE io_constants, ONLY:                                                       &
  mode_read, max_dim_var

USE parallel_mod, ONLY:                                                       &
  master_task_id, is_master_task, scatter_land_field

USE string_utils_mod, ONLY:                                                   &
  to_string

USE file_mod, ONLY:                                                           &
  file_handle, file_open, file_introspect,                                    &
  file_inquire_dim, file_inquire_var, file_read_var,                          &
  file_close

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Check that the given file is a JULES dump compatible with the current
!   run, and read the given identifiers from it.
!   Note that the reading of the dump is done by the master task and the
!   results scattered to other tasks. This means that dumps written with
!   different amounts of tasks should be interchangable.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
CHARACTER(LEN=*) :: file_name  ! The dump file
CHARACTER(LEN=*) :: identifiers(:)  ! The model identifiers for the
                                    ! variables to define

! Local parameters.
LOGICAL, PARAMETER :: l_reading_true = .TRUE.
  ! A value of .TRUE. that is passed to argument l_reading of subroutine
  ! get_dim_info to show that it is being called in connection with reading
  ! (rather than writing) a dump.

! Work variables
TYPE(file_handle) :: FILE
  ! The opened file

INTEGER :: nvars
  ! The number of variables we are processing

INTEGER :: dim_size_file
  ! The size of the dimension currently being
  ! processed in the file

! Used when defining dimensions and variables
INTEGER :: ndims
  ! The number of dimensions the current variable has
CHARACTER(LEN=max_sdf_name_len) :: dim_names(max_dim_var)
  ! The dimension names the current variable should use
INTEGER :: dim_sizes(max_dim_var)
  ! The dimension sizes for the current variable
INTEGER :: dim_ids(max_dim_var)
  ! The dimension ids for the current variable as
  ! calculated from file_inquire_dim
LOGICAL :: is_record_dim
  ! Detects if the current dimension is a record dim
INTEGER :: ndims_file
  ! The number of dimensions the variable has in file
  ! Compared to ndims above for each variable
INTEGER :: dim_ids_file(max_dim_var)
  ! The ids of the dimensions the variable has in the file
  ! Compared to dim_ids above to verify the variable has the
  ! correct dimensions

INTEGER :: var_ids(SIZE(identifiers))
  ! The ids of the variables in the dump file

LOGICAL :: is_record_var
  ! Indicates if a variable uses the record dimension

LOGICAL :: l_read_from_dump
  ! Used to bypass checking the dimensions of ancil variables that are
  ! not to be read from the dump.

REAL :: frac_irrig_all_tiles_real
REAL :: set_irrfrac_irrtiles_real
REAL :: irrtiles_real(npft)
REAL :: nirrtile_real

! Real versions of integer valued variables
REAL :: nsnow_real(land_pts, nsurft)
REAL, ALLOCATABLE :: seed_rain_real(:)
! The number of elements of the "seed_rain" array
INTEGER :: seed_rain_len

INTEGER :: i, j, m, n, s  ! Loop counters

INTEGER :: error  ! Variable to collect MPI errors - most MPI
                  ! implementations bail on error, so this is not checked.


! Arrays to hold global land points version of data read in master task
! before scattering
REAL, ALLOCATABLE :: global_data_1d(:)     ! Data with no vertical levels
REAL, ALLOCATABLE :: global_data_2d(:,:)   ! With one vertical level
REAL, ALLOCATABLE :: global_data_3d(:,:,:) ! With two "vertical" levels
                                             ! e.g. snow variables
REAL, ALLOCATABLE :: global_data_4d(:,:,:,:) ! With 3 "vertical" levels

!-----------------------------------------------------------------------------

nvars = SIZE(identifiers)
IF (ALLOCATED(progs_data%seed_rain)) THEN
  seed_rain_len = SIZE(progs_data%seed_rain)
ELSE
  seed_rain_len = 0
END IF

!-----------------------------------------------------------------------------
! In the master task only, we open the file and check that the correct
! dimensions exist and are of a size compatible with this run
!-----------------------------------------------------------------------------
IF ( is_master_task() ) THEN
  !-----------------------------------------------------------------------------
  ! We use the lowest level file API here, as we don't want to impose the input
  ! grid
  !-----------------------------------------------------------------------------
  FILE=file_open(file_name, mode_read)

  ! We want to auto-detect the dimensions and variables in the file
  CALL file_introspect(FILE)

  DO i = 1,nvars

    !-----------------------------------------------------------------------
    ! Get information about the dimensions used by the variable.
    ! The argument l_reading_true shows that we are reading (not writing) a
    ! dump.
    !-----------------------------------------------------------------------
    CALL get_dim_info( l_reading_true, identifiers(i), ndims,  dim_sizes,     &
                       dim_names, l_read_from_dump )

    !-------------------------------------------------------------------------
    ! Check the dimensions exist and have the correct size
    !-------------------------------------------------------------------------
    IF ( l_read_from_dump ) THEN
      DO j = 1,ndims
        ! Retrive information about the dimension from the file we store the id
        ! for use outside this loop
        CALL file_inquire_dim(                                                &
          FILE, dim_names(j), dim_ids(j), dim_size_file, is_record_dim        &
        )

        ! Check that we found a dimension
        IF ( dim_ids(j) < 0 )                                                 &
          CALL log_fatal("read_dump",                                         &
                         "Could not find expected dimension '" //             &
                         TRIM(dim_names(j)) // "' in dump file")

        ! Check that the dimension is not a record dimension (there shouldn't
        ! be one in dump files).
        IF ( is_record_dim )                                                  &
          CALL log_fatal("read_dump",                                         &
                         "Dimension '" // TRIM(dim_names(j)) // "' is a " //  &
                         "record dimension - should not exist in dump file")

        ! Check that the dimension has the correct size
        IF ( dim_size_file /= dim_sizes(j) )                                  &
          CALL log_fatal("read_dump",                                         &
                         "Dimension '" // TRIM(dim_names(j)) // "' has " //   &
                         "size incompatible with current run (required: " //  &
                         TRIM(to_string(dim_sizes(j))) // ", found: " //      &
                         TRIM(to_string(dim_size_file)) // ")")
      END DO  ! dims

      !-----------------------------------------------------------------------
      ! Check that the variable exists and has the correct dimensions
      !-----------------------------------------------------------------------
      ! Retrieve information about the variable from the file
      CALL file_inquire_var(                                                  &
        FILE, identifiers(i), var_ids(i), ndims_file, dim_ids_file,           &
        is_record_var                                                         &
      )

      ! Check that we found a variable
      IF ( var_ids(i) < 1 )                                                   &
        CALL log_fatal("read_dump",                                           &
                       "Failed to find requested variable '" //               &
                       TRIM(identifiers(i)) // "' in dump file")

      ! Check that the number of dimensions match
      IF ( ndims_file /= ndims )                                              &
        CALL log_fatal("read_dump",                                           &
                       "Variable '" // TRIM(identifiers(i)) // "' has " //    &
                       "incorrect number of dimensions in dump file (" //     &
                       "expected: " // TRIM(to_string(ndims)) // ", " //      &
                       "found: " // TRIM(to_string(ndims_file)) // ")")

      ! Check that the dimension ids match
      IF ( .NOT. ALL(dim_ids(1:ndims) == dim_ids_file(1:ndims)) )             &
        CALL log_fatal("read_dump",                                           &
                       "Variable '" // TRIM(identifiers(i)) // "' has " //    &
                       "incorrect dimensions in dump file")

    END IF  !  l_read_from_dump

  END DO  ! vars

END IF  ! MASTER TASK

!-----------------------------------------------------------------------------
! Set the requested variables from the file
!
! This is done by reading the value of the variable on global land points
! in the master task, then scattering it to the other tasks
!
! We assume that if the file passed all the checks on dimensions above, then
! it will be fine to fill variables here (i.e. we don't check the dimensions
! associated with the variables)
!-----------------------------------------------------------------------------
! Allocate the global data arrays
IF ( is_master_task() ) THEN
  ALLOCATE(global_data_1d(global_land_pts))
  ALLOCATE(global_data_2d(global_land_pts, MAX(npft, dim_cs1, sm_levels,      &
                                               nsurft, ntype, ns_deep,        &
                                               nsoilt, dim_ch4layer)))
  ALLOCATE(global_data_3d(global_land_pts,                                    &
                          MAX(nsurft, dim_cslayer, nsoilt),                   &
                          MAX(nsmax, dim_cs1, sm_levels)))
  ALLOCATE(global_data_4d(global_land_pts, nsoilt,                            &
                          MAX(nsurft, dim_cslayer),                           &
                          MAX(nsmax, dim_cs1)))
ELSE
  ALLOCATE(global_data_1d(1))
  ALLOCATE(global_data_2d(1,1))
  ALLOCATE(global_data_3d(1,1,1))
  ALLOCATE(global_data_4d(1,1,1,1))
END IF

DO i = 1,nvars

  !---------------------------------------------------------------------------
  ! In the master task, read the global data
  !---------------------------------------------------------------------------
  IF ( is_master_task() ) THEN
    SELECT CASE ( identifiers(i) )

      !----------------------------------------------------------------------
      ! If it is a land_pts array with no levels associated, read into the
      ! global_data_1d array.
      !----------------------------------------------------------------------
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
           'fire_nesterov')
      CALL file_read_var(FILE, var_ids(i), global_data_1d)

      !----------------------------------------------------------------------
      ! If it is a variable with one or more levels,
      ! read the appropriate number of levels into global_data_2d/3d/4d.
      !----------------------------------------------------------------------
    CASE ( 'toppdm%sthzw_soilt', 'toppdm%zw_soilt' )
      CALL file_read_var(FILE, var_ids(i), global_data_2d(:,1:nsoilt))

    CASE ( 'canht', 'lai' )
      CALL file_read_var(FILE, var_ids(i), global_data_2d(:,1:npft))

    CASE ( 'cropdvi', 'croprootc', 'cropharvc', 'cropreservec',               &
           'croplai', 'cropcanht' )
      CALL file_read_var(FILE, var_ids(i), global_data_2d(:,1:ncpft))

    CASE ( 'sthuf', 't_soil', 'sthu_irr' )
      CALL file_read_var(FILE, var_ids(i),                                    &
                         global_data_2d(:,1:sm_levels))

    CASE ( 'sthuf_soilt', 't_soil_soilt', 'sthu_irr_soilt')
      CALL file_read_var(FILE, var_ids(i),                                    &
                         global_data_3d(:,1:nsoilt,1:sm_levels))

    CASE ( 'n_inorg' )
      CALL file_read_var(FILE, var_ids(i),                                    &
                         global_data_2d(:,1:dim_cslayer))

    CASE ( 'n_inorg_soilt' )
      CALL file_read_var(FILE, var_ids(i),                                    &
                         global_data_3d(:,1:nsoilt,1:dim_cslayer))

    CASE ( 'substr_ch4','mic_ch4','mic_act_ch4','acclim_ch4' )
      CALL file_read_var(FILE, var_ids(i),                                    &
                         global_data_2d(:,1:dim_ch4layer))

    CASE ( 'tsoil_deep' )
      CALL file_read_var(FILE, var_ids(i), global_data_2d(:,1:ns_deep))

    CASE ( 'canopy', 'nsnow', 'rgrain', 'rho_snow', 'snow_tile',              &
           'snow_depth', 'snow_grnd', 'tstar_tile', 'tsurf_elev_surft' )
      CALL file_read_var(FILE, var_ids(i), global_data_2d(:,1:nsurft))

    CASE ( 'rgrainl', 'snow_ds', 'snow_ice', 'snow_liq', 'tsnow' )
      CALL file_read_var(FILE, var_ids(i),                                    &
                         global_data_3d(:,1:nsurft,1:nsmax))

    CASE ( 'cs','ns' )
      CALL file_read_var(FILE, var_ids(i),                                    &
                         global_data_3d(:,1:dim_cslayer,1:dim_cs1))

    CASE ( 'cs_soilt', 'ns_soilt' )
      CALL file_read_var(FILE, var_ids(i),                                    &
                         global_data_4d(:,1:nsoilt,1:dim_cslayer,1:dim_cs1))

    CASE ( 'n_soil' )
      CALL file_read_var(FILE, var_ids(i),                                    &
                         global_data_3d(:,1:dim_cslayer,1:dim_soil_n_pool))

    CASE ( 'n_soil_soilt' )
      CALL file_read_var(FILE, var_ids(i),                                    &
                         global_data_4d(:,1:nsoilt,1:dim_cslayer,             &
                                        1:dim_soil_n_pool))

      ! Cases for IMOGEN variables
      ! Each task runs its own version of IMOGEN - these variables are
      ! broadcast to all tasks below.
    CASE ( 'co2_ppmv' )
      CALL file_read_var(FILE, var_ids(i), co2_ppmv)

    CASE ( 'co2_change_ppmv' )
      CALL file_read_var(FILE, var_ids(i), co2_change_ppmv)

    CASE ( 'dtemp_o' )
      CALL file_read_var(FILE, var_ids(i), dtemp_o)

    CASE ( 'fa_ocean' )
      CALL file_read_var(FILE, var_ids(i), fa_ocean)

    CASE ( 'seed_rain' )
      IF (ALLOCATED(progs_data%seed_rain)) THEN
        ALLOCATE(seed_rain_real(seed_rain_len))
        CALL file_read_var(FILE, var_ids(i), seed_rain_real)
      ELSE
        CALL log_fatal("read_dump",                                           &
                       "seed_rain needs to be allocated before use")
      END IF


    CASE ( 'ch4_ppbv' )
      CALL file_read_var(FILE, var_ids(i), ch4_ppbv)

      ! River routing variables
    CASE ( 'rivers_sto_rp' )
      CALL file_read_var(FILE, var_ids(i), rivers_sto_rp)

    CASE ( 'rfm_surfstore_rp' )
      CALL file_read_var(FILE, var_ids(i), rfm_surfstore_rp)

    CASE ( 'rfm_substore_rp' )
      CALL file_read_var(FILE, var_ids(i), rfm_substore_rp)

    CASE ( 'rfm_flowin_rp' )
      CALL file_read_var(FILE, var_ids(i), rfm_flowin_rp)

    CASE ( 'rfm_bflowin_rp' )
      CALL file_read_var(FILE, var_ids(i), rfm_bflowin_rp)

      !----------------------------------------------------------------------
      ! Ancillary variables
      !----------------------------------------------------------------------

      ! latlon ancil namelist
    CASE ( 'grid_area' )
      IF ( ancil_dump_read%latlon ) THEN
        CALL file_read_var(FILE, var_ids(i), global_data_1d)
      END IF

      !Frac ancil namelist
    CASE ( 'frac' )
      IF ( ancil_dump_read%frac ) THEN
        CALL file_read_var(FILE, var_ids(i), global_data_2d(:,1:ntype))
      END IF

      ! Vegetation properties ancil namelist
    CASE ( 't_growth_gb' )
      IF ( ancil_dump_read%vegetation_props ) THEN
        CALL file_read_var(FILE, var_ids(i), global_data_1d)
      END IF

      !Soil properties ancil namelist
    CASE ( 'b      ', 'sathh  ', 'satcon ', 'sm_sat ', 'sm_crit',             &
           'sm_wilt', 'hcap   ', 'hcon   ' )
      IF ( ancil_dump_read%soil_props ) THEN
        CALL file_read_var(FILE, var_ids(i),                                  &
                           global_data_2d(:,1:sm_levels))
      END IF

    CASE ( 'b_soilt', 'sathh_soilt', 'satcon_soilt', 'sm_sat_soilt',          &
           'sm_crit_soilt', 'sm_wilt_soilt', 'hcap_soilt', 'hcon_soilt' )
      IF ( ancil_dump_read%soil_props ) THEN
        CALL file_read_var(FILE, var_ids(i),                                  &
                           global_data_3d(:,1:nsoilt,1:sm_levels))
      END IF

    CASE ( 'albsoil' )
      IF ( ancil_dump_read%soil_props ) THEN
        CALL file_read_var(FILE, var_ids(i), global_data_1d)
      END IF

    CASE ( 'albsoil_soilt' )
      IF ( ancil_dump_read%soil_props ) THEN
        CALL file_read_var(FILE, var_ids(i), global_data_2d(:,1:nsoilt))
      END IF

    CASE ( 'clay', 'soil_ph' )
      IF ( ancil_dump_read%soil_props ) THEN
        CALL file_read_var(FILE, var_ids(i),                                  &
                           global_data_2d(:,1:dim_cslayer))
      END IF

    CASE ( 'clay_soilt', 'soil_ph_soilt' )
      IF ( ancil_dump_read%soil_props ) THEN
        CALL file_read_var(FILE, var_ids(i),                                  &
                           global_data_3d(:,1:nsoilt,1:dim_cslayer))
      END IF

      !Topmodel ancil namelist
    CASE ( 'fexp   ', 'ti_mean', 'ti_sig ' )
      IF ( ancil_dump_read%top ) THEN
        CALL file_read_var(FILE, var_ids(i), global_data_1d)
      END IF

    CASE ( 'toppdm%fexp_soilt', 'toppdm%ti_mean_soilt', 'toppdm%ti_sig_soilt' )
      IF ( ancil_dump_read%top ) THEN
        CALL file_read_var(FILE, var_ids(i), global_data_2d(:,1:nsoilt))
      END IF

      !Agric ancil namelist
    CASE ( 'frac_agr', 'frac_past' )
      IF ( ancil_dump_read%agric ) THEN
        CALL file_read_var(FILE, var_ids(i), global_data_1d)
      END IF

      !Crop props ancillaries namelist
    CASE ( 'cropsowdate', 'cropttveg  ', 'cropttrep  ','croplatestharvdate')
      IF ( ancil_dump_read%crop_props ) THEN
        CALL file_read_var(FILE, var_ids(i), global_data_2d(:,1:ncpft))
      END IF

      !Irrigation ancillaries namelist
    CASE ( 'frac_irrig', 'irrfrac_irrtiles' )
      IF ( ancil_dump_read%irrig ) THEN
        CALL file_read_var(FILE, var_ids(i), global_data_1d)
      END IF

    CASE ( 'frac_irr_all_tiles' )
      IF ( ancil_dump_read%irrig ) THEN
        CALL file_read_var(FILE, var_ids(i), frac_irrig_all_tiles_real)
      END IF

    CASE ( 'irrtiles' )
      IF ( ancil_dump_read%irrig ) THEN
        CALL file_read_var(FILE, var_ids(i), irrtiles_real)
      END IF

    CASE ( 'nirrtile' )
      IF ( ancil_dump_read%irrig ) THEN
        CALL file_read_var(FILE, var_ids(i), nirrtile_real)
      END IF

    CASE ( 'set_irrfrac_on_irrtiles' )
      IF ( ancil_dump_read%irrig ) THEN
        CALL file_read_var(FILE, var_ids(i), set_irrfrac_irrtiles_real)
      END IF

      ! Water resources properties ancil namelist
    CASE ( 'conveyance_loss', 'irrig_eff', 'sfc_water_frac' )
      IF ( ancil_dump_read%vegetation_props ) THEN
        CALL file_read_var(FILE, var_ids(i), global_data_1d)
      END IF

      !CO2 ancil namelist
    CASE ( 'co2_mmr' )
      IF ( ancil_dump_read%co2 ) THEN
        CALL file_read_var(FILE, var_ids(i), co2_mmr)
      END IF

    CASE DEFAULT
      CALL log_fatal("read_dump",                                             &
                     "Unexpected variable in dump - " //                      &
                     TRIM(identifiers(i)))
    END SELECT
  END IF  ! MASTER TASK

  !---------------------------------------------------------------------------
  ! Now scatter the variables into their final destinations.
  ! Note that scatter_land_field can only scatter one land_pts array at a time
  ! so to scatter variables with multiple levels we must loop
  !---------------------------------------------------------------------------
  SELECT CASE ( identifiers(i) )
  CASE ( 'gs' )
    CALL scatter_land_field(global_data_1d, progs%gs_gb)

  CASE ( 'sthzw' )
    IF ( l_tile_soil .AND. l_broadcast_soilt ) THEN
      DO m = 1, nsoilt
        CALL scatter_land_field(global_data_1d, toppdm%sthzw_soilt(:,m))
      END DO
    ELSE !Case if nsoilt == 1, so OK to hardwire the 2nd dimension to 1
      CALL scatter_land_field(global_data_1d, toppdm%sthzw_soilt(:,1))
    END IF

  CASE ( 'toppdm%sthzw_soilt' )
    DO m = 1, nsoilt
      CALL scatter_land_field(global_data_2d(:,m), toppdm%sthzw_soilt(:,m))
    END DO

  CASE ( 'zw' )
    IF ( l_tile_soil .AND. l_broadcast_soilt ) THEN
      DO m = 1, nsoilt
        CALL scatter_land_field(global_data_1d, toppdm%zw_soilt(:,m))
      END DO
    ELSE !Case if nsoilt == 1, so OK to hardwire the 2nd dimension to 1
      CALL scatter_land_field(global_data_1d, toppdm%zw_soilt(:,1))
    END IF

  CASE ( 'toppdm%zw_soilt' )
    DO m = 1,nsoilt
      CALL scatter_land_field(global_data_2d(:,m), toppdm%zw_soilt(:,m))
    END DO

  CASE ( 'cv' )
    CALL scatter_land_field(global_data_1d, trifctltype%cv_gb)

  CASE ( 'frac_agr_prev' )
    CALL scatter_land_field(global_data_1d, progs%frac_agr_prev_gb)

  CASE ( 'frac_past_prev' )
    CALL scatter_land_field(global_data_1d, progs%frac_past_prev_gb)

  CASE ( 'wood_prod_fast' )
    CALL scatter_land_field(global_data_1d, progs%wood_prod_fast_gb)

  CASE ( 'wood_prod_med' )
    CALL scatter_land_field(global_data_1d, progs%wood_prod_med_gb)

  CASE ( 'wood_prod_slow' )
    CALL scatter_land_field(global_data_1d, progs%wood_prod_slow_gb)

  CASE ( 'n_inorg' )
    IF ( l_tile_soil .AND. l_broadcast_soilt ) THEN
      DO m = 1, nsoilt
        DO n = 1,dim_cslayer
          CALL scatter_land_field(global_data_2d(:,n),                        &
                                  progs%n_inorg_soilt_lyrs(:,m,n))
        END DO
      END DO
    ELSE !Case if nsoilt == 1, so OK to hardwire the 2nd dimension to 1
      DO n = 1,dim_cslayer
        CALL scatter_land_field(global_data_2d(:,n),                          &
                                progs%n_inorg_soilt_lyrs(:,1,n))
      END DO
    END IF

  CASE ( 'n_inorg_soilt' )
    DO m = 1, nsoilt
      DO n = 1,dim_cslayer
        CALL scatter_land_field(global_data_3d(:,m,n),                        &
                                progs%n_inorg_soilt_lyrs(:,m,n))
      END DO
    END DO

  CASE ( 'substr_ch4' )
    DO n = 1,dim_ch4layer
      CALL scatter_land_field(global_data_2d(:,n), progs%substr_ch4(:,n))
    END DO

  CASE ( 'mic_ch4' )
    DO n = 1,dim_ch4layer
      CALL scatter_land_field(global_data_2d(:,n), progs%mic_ch4(:,n))
    END DO

  CASE ( 'mic_act_ch4' )
    DO n = 1,dim_ch4layer
      CALL scatter_land_field(global_data_2d(:,n), progs%mic_act_ch4(:,n))
    END DO

  CASE ( 'acclim_ch4' )
    DO n = 1,dim_ch4layer
      CALL scatter_land_field(global_data_2d(:,n), progs%acclim_ch4(:,n))
    END DO

  CASE ( 'canht' )
    DO n = 1,npft
      CALL scatter_land_field(global_data_2d(:,n), progs%canht_pft(:,n))
    END DO

  CASE ( 'lai' )
    DO n = 1,npft
      CALL scatter_land_field(global_data_2d(:,n), progs%lai_pft(:,n))
    END DO

  CASE ( 'cropdvi' )
    DO n = 1,ncpft
      CALL scatter_land_field(global_data_2d(:,n), crop_vars%dvi_cpft(:,n))
    END DO

  CASE ( 'croprootc' )
    DO n = 1,ncpft
      CALL scatter_land_field(global_data_2d(:,n), crop_vars%rootc_cpft(:,n))
    END DO

  CASE ( 'cropharvc' )
    DO n = 1,ncpft
      CALL scatter_land_field(global_data_2d(:,n), crop_vars%harvc_cpft(:,n))
    END DO

  CASE ( 'cropreservec' )
    DO n = 1,ncpft
      CALL scatter_land_field(global_data_2d(:,n), crop_vars%reservec_cpft(:,n))
    END DO

  CASE ( 'croplai' )
    DO n = 1,ncpft
      CALL scatter_land_field(global_data_2d(:,n), crop_vars%croplai_cpft(:,n))
    END DO

  CASE ( 'cropcanht' )
    DO n = 1,ncpft
      CALL scatter_land_field(global_data_2d(:,n), crop_vars%cropcanht_cpft(:,n))
    END DO

  CASE ( 'cs' )
    IF ( l_tile_soil .AND. l_broadcast_soilt ) THEN
      DO s = 1,nsoilt !using s rather than m as usual
        DO n = 1,dim_cs1
          DO m = 1,dim_cslayer
            CALL scatter_land_field(global_data_3d(:,m,n),                    &
                                    progs%cs_pool_soilt(:,s,m,n))
          END DO
        END DO
      END DO
    ELSE !Case if nsoilt == 1, so OK to hardwire the 2nd dimension to 1
      DO n = 1,dim_cs1
        DO m = 1,dim_cslayer
          CALL scatter_land_field(global_data_3d(:,m,n),                      &
                                  progs%cs_pool_soilt(:,1,m,n))
        END DO
      END DO
    END IF

  CASE ( 'cs_soilt' )
    DO s = 1,nsoilt !using s rather than m as usual
      DO n = 1,dim_cs1
        DO m = 1,dim_cslayer
          CALL scatter_land_field(global_data_4d(:,s,m,n),                    &
                                  progs%cs_pool_soilt(:,s,m,n))
        END DO
      END DO
    END DO

  CASE ( 'ns' )
    DO n = 1,dim_cs1
      DO m = 1,dim_cslayer
        CALL scatter_land_field(global_data_3d(:,m,n),                        &
                                progs%ns_pool_gb(:,m,n))
      END DO
    END DO

  CASE ( 'sthuf' )
    ! sthuf is held in sthu until it is processed
    IF ( l_tile_soil .AND. l_broadcast_soilt ) THEN
      DO m = 1,nsoilt
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_2d(:,n), jules_vars%sthuf_soilt(:,m,n))
        END DO
      END DO
    ELSE !Case if nsoilt == 1, so OK to hardwire the 2nd dimension to 1
      DO n = 1,sm_levels
        CALL scatter_land_field(global_data_2d(:,n), jules_vars%sthuf_soilt(:,1,n))
      END DO
    END IF

  CASE ( 'sthuf_soilt' )
    ! sthuf is held in sthu until it is processed
    DO m = 1,nsoilt
      DO n = 1,sm_levels
        CALL scatter_land_field(global_data_3d(:,m,n), jules_vars%sthuf_soilt(:,m,n))
      END DO
    END DO

  CASE ( 't_soil' )
    IF ( l_tile_soil .AND. l_broadcast_soilt ) THEN
      DO m = 1,nsoilt
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_2d(:,n),                        &
                                  progs%t_soil_soilt(:,m,n))
        END DO
      END DO
    ELSE !Case if nsoilt == 1, so OK to hardwire the 2nd dimension to 1
      DO n = 1,sm_levels
        CALL scatter_land_field(global_data_2d(:,n), progs%t_soil_soilt(:,1,n))
      END DO
    END IF

  CASE ( 't_soil_soilt' )
    DO m = 1,nsoilt
      DO n = 1,sm_levels
        CALL scatter_land_field(global_data_3d(:,m,n),                        &
                                progs%t_soil_soilt(:,m,n))
      END DO
    END DO

  CASE ( 'tsoil_deep' )
    DO n = 1,ns_deep
      CALL scatter_land_field(global_data_2d(:,n), progs%tsoil_deep_gb(:,n))
    END DO

  CASE ( 'sthu_irr' )
    IF ( l_tile_soil .AND. l_broadcast_soilt ) THEN
      DO m = 1,nsoilt
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_2d(:,n),                        &
                                  crop_vars%sthu_irr_soilt(:,m,n))
        END DO
      END DO
    ELSE !Case if nsoilt == 1, so OK to hardwire the 2nd dimension to 1
      DO n = 1,sm_levels
        CALL scatter_land_field(global_data_2d(:,n),                          &
                                crop_vars%sthu_irr_soilt(:,1,n))
      END DO
    END IF

  CASE ( 'sthu_irr_soilt' )
    DO m = 1,nsoilt
      DO n = 1,sm_levels
        CALL scatter_land_field(global_data_3d(:,m,n),                        &
                                crop_vars%sthu_irr_soilt(:,m,n))
      END DO
    END DO

  CASE ( 'canopy' )
    DO n = 1,nsurft
      CALL scatter_land_field(global_data_2d(:,n), progs%canopy_surft(:,n))
    END DO

  CASE ( 'nsnow' )
    DO n = 1,nsurft
      CALL scatter_land_field(global_data_2d(:,n), nsnow_real(:,n))
    END DO
    progs%nsnow_surft(:,:) = NINT(nsnow_real(:,:))

  CASE ( 'rgrain' )
    DO n = 1,nsurft
      CALL scatter_land_field(global_data_2d(:,n), progs%rgrain_surft(:,n))
    END DO

  CASE ( 'rho_snow' )
    DO n = 1,nsurft
      CALL scatter_land_field(global_data_2d(:,n),                            &
                              progs%rho_snow_grnd_surft(:,n))
    END DO

  CASE ( 'snow_tile' )
    DO n = 1,nsurft
      CALL scatter_land_field(global_data_2d(:,n), progs%snow_surft(:,n))
    END DO

  CASE ( 'snow_depth' )
    DO n = 1,nsurft
      CALL scatter_land_field(global_data_2d(:,n),                            &
                              progs%snowdepth_surft(:,n))
    END DO

  CASE ( 'snow_grnd' )
    DO n = 1,nsurft
      CALL scatter_land_field(global_data_2d(:,n), progs%snow_grnd_surft(:,n))
    END DO

  CASE ( 'tstar_tile' )
    DO n = 1,nsurft
      CALL scatter_land_field(global_data_2d(:,n), progs%tstar_surft(:,n))
    END DO

  CASE ( 'tsurf_elev_surft' )
    DO n = 1,nsurft
      CALL scatter_land_field(global_data_2d(:,n),                            &
                              progs%tsurf_elev_surft(:,n))
    END DO

  CASE ( 'rgrainl' )
    DO n = 1,nsmax
      DO m = 1,nsurft
        CALL scatter_land_field(global_data_3d(:,m,n),                        &
                                progs%rgrainl_surft(:,m,n))
      END DO
    END DO

  CASE ( 'snow_ds' )
    DO n = 1,nsmax
      DO m = 1,nsurft
        CALL scatter_land_field(global_data_3d(:,m,n), progs%ds_surft(:,m,n))
      END DO
    END DO

  CASE ( 'snow_ice' )
    DO n = 1,nsmax
      DO m = 1,nsurft
        CALL scatter_land_field(global_data_3d(:,m,n),                        &
                                progs%sice_surft(:,m,n))
      END DO
    END DO


  CASE ( 'snow_liq' )
    DO n = 1,nsmax
      DO m = 1,nsurft
        CALL scatter_land_field(global_data_3d(:,m,n),                        &
                                progs%sliq_surft(:,m,n))
      END DO
    END DO


  CASE ( 'tsnow' )
    DO n = 1,nsmax
      DO m = 1,nsurft
        CALL scatter_land_field(global_data_3d(:,m,n),                        &
                                progs%tsnow_surft(:,m,n))
      END DO
    END DO

    ! Fire and metstats
  CASE ( 'temp_max_00h_r' )
    CALL scatter_land_field(global_data_1d,                                   &
                            metstats_prog(:)%temp_max_00h%run)

  CASE ( 'temp_ave_00h_r' )
    CALL scatter_land_field(global_data_1d,                                   &
                            metstats_prog(:)%temp_ave_00h%run)

  CASE ( 'prec_tot_00h_r' )
    CALL scatter_land_field(global_data_1d,                                   &
                            metstats_prog(:)%prec_tot_00h%run)

  CASE ( 'prec_tot_12h_r' )
    CALL scatter_land_field(global_data_1d,                                   &
                            metstats_prog(:)%prec_tot_12h%run)

  CASE ( 'rhum_min_00h_r' )
    CALL scatter_land_field(global_data_1d,                                   &
                            metstats_prog(:)%rhum_min_00h%run)

  CASE ( 'dewp_ave_00h_r' )
    CALL scatter_land_field(global_data_1d,                                   &
                            metstats_prog(:)%dewp_ave_00h%run)

  CASE ( 'wind_ave_00h_r' )
    CALL scatter_land_field(global_data_1d,                                   &
                            metstats_prog(:)%wind_ave_00h%run)

  CASE ( 'temp_max_00h' )
    CALL scatter_land_field(global_data_1d,                                   &
                            metstats_prog(:)%temp_max_00h%fin)

  CASE ( 'temp_ave_00h' )
    CALL scatter_land_field(global_data_1d,                                   &
                            metstats_prog(:)%temp_ave_00h%fin)

  CASE ( 'temp_ave_nday' )
    CALL scatter_land_field(global_data_1d,                                   &
                            metstats_prog(:)%temp_ave_nday%fin)

  CASE ( 'temp_pnt_12h' )
    CALL scatter_land_field(global_data_1d,                                   &
                            metstats_prog(:)%temp_pnt_12h%fin)

  CASE ( 'prec_tot_00h' )
    CALL scatter_land_field(global_data_1d,                                   &
                            metstats_prog(:)%prec_tot_00h%fin)

  CASE ( 'prec_tot_12h' )
    CALL scatter_land_field(global_data_1d,                                   &
                            metstats_prog(:)%prec_tot_12h%fin)

  CASE ( 'rhum_min_00h' )
    CALL scatter_land_field(global_data_1d,                                   &
                            metstats_prog(:)%rhum_min_00h%fin)

  CASE ( 'rhum_pnt_12h' )
    CALL scatter_land_field(global_data_1d,                                   &
                            metstats_prog(:)%rhum_pnt_12h%fin)

  CASE ( 'dewp_ave_00h' )
    CALL scatter_land_field(global_data_1d,                                   &
                            metstats_prog(:)%dewp_ave_00h%fin)

  CASE ( 'wind_ave_00h' )
    CALL scatter_land_field(global_data_1d,                                   &
                            metstats_prog(:)%wind_ave_00h%fin)

  CASE ( 'wind_pnt_12h' )
    CALL scatter_land_field(global_data_1d,                                   &
                            metstats_prog(:)%wind_pnt_12h%fin)

    ! Fire module variables- land points only
  CASE ( 'fire_mcarthur_r_dr' )
    CALL scatter_land_field(global_data_1d, fire_prog(:)%mcarthur%r_dr)

  CASE ( 'fire_mcarthur_n_dr' )
    CALL scatter_land_field(global_data_1d, fire_prog(:)%mcarthur%n_dr)

  CASE ( 'fire_canadian_ffmc' )
    CALL scatter_land_field(global_data_1d, fire_prog(:)%canadian%ffmc)

  CASE ( 'fire_canadian_ffmc_mois' )
    CALL scatter_land_field(global_data_1d,                                   &
                            fire_prog(:)%canadian%ffmc_mois)

  CASE ( 'fire_canadian_dmc' )
    CALL scatter_land_field(global_data_1d, fire_prog(:)%canadian%dmc)

  CASE ( 'fire_canadian_dc' )
    CALL scatter_land_field(global_data_1d, fire_prog(:)%canadian%dc)

  CASE ( 'fire_nesterov' )
    CALL scatter_land_field(global_data_1d, fire_prog(:)%nesterov%findex)

    ! ECOSSE variables.
  CASE ( 'n_soil' )
    IF ( l_tile_soil .AND. l_broadcast_soilt ) THEN
      DO s = 1,nsoilt !using s rather than m as usual
        DO n = 1,dim_soil_n_pool
          DO m = 1,dim_cslayer
            CALL scatter_land_field(global_data_3d(:,m,n),                    &
                                    soilecosse%n_soil_pool_soilt(:,s,m,n))
          END DO
        END DO
      END DO
    ELSE !Case if nsoilt == 1, so OK to hardwire the 2nd dimension to 1
      DO n = 1,dim_soil_n_pool
        DO m = 1,dim_cslayer
          CALL scatter_land_field(global_data_3d(:,m,n),                      &
                                  soilecosse%n_soil_pool_soilt(:,1,m,n))
        END DO
      END DO
    END IF

  CASE ( 'n_soil_soilt' )
    DO s = 1,nsoilt !using s rather than m as usual
      DO n = 1,dim_soil_n_pool
        DO m = 1,dim_cslayer
          CALL scatter_land_field(global_data_4d(:,s,m,n),                    &
                                  soilecosse%n_soil_pool_soilt(:,s,m,n))
        END DO
      END DO
    END DO

    ! IMOGEN variables are just broadcast to all tasks
  CASE ( 'co2_ppmv' )
    CALL mpi_bcast(co2_ppmv, 1, mpi_real,                                     &
                   master_task_id, mpi_comm_world, error)

  CASE ( 'co2_change_ppmv' )
    CALL mpi_bcast(co2_change_ppmv, 1, mpi_real,                              &
                   master_task_id, mpi_comm_world, error)

  CASE ( 'dtemp_o' )
    CALL mpi_bcast(dtemp_o, n_olevs, mpi_real,                                &
                   master_task_id, mpi_comm_world, error)

  CASE ( 'fa_ocean' )
    CALL mpi_bcast(fa_ocean, nfarray, mpi_real,                               &
                   master_task_id, mpi_comm_world, error)

  CASE ( 'seed_rain' )
    IF (ALLOCATED(progs_data%seed_rain)) THEN
      ! The _real version should be allocated on the "master" task but not
      ! the other tasks.
      IF ( .NOT. ALLOCATED(seed_rain_real)) THEN
        ALLOCATE(seed_rain_real(seed_rain_len))
      END IF
      ! The "master" task has read the field into seed_rain_real
      CALL mpi_bcast(seed_rain_real, seed_rain_len, mpi_real,                 &
                     master_task_id, mpi_comm_world, error)
      progs_data%seed_rain(:) = NINT(seed_rain_real(:))
      DEALLOCATE(seed_rain_real)
    ELSE
      CALL log_fatal("read_dump",                                             &
          "Tried to read seed_rain but internal array was not allocated")
    END IF
  CASE ( 'ch4_ppbv' )
    CALL mpi_bcast(ch4_ppbv, 1, mpi_real,                                     &
                   master_task_id, mpi_comm_world, error)

    ! River routing variables
  CASE ( 'rivers_sto_rp', 'rfm_surfstore_rp', 'rfm_substore_rp',              &
         'rfm_flowin_rp', 'rfm_bflowin_rp' )
    ! nothing to do

    !-------------------------------------------------------------------------
    ! Ancillary variables
    !-------------------------------------------------------------------------

    ! latlon ancil namelist
  CASE ( 'grid_area' )
    IF ( ancil_dump_read%latlon ) THEN
      CALL scatter_land_field(global_data_1d, grid_area_ij)
    END IF

    !Frac ancil namelist
  CASE ( 'frac' )
    IF ( ancil_dump_read%frac ) THEN
      DO n = 1,ntype
        CALL scatter_land_field(global_data_2d(:,n), ainfo%frac_surft(:,n))
      END DO
    END IF

    ! Vegetation properties ancil namelist
  CASE ( 't_growth_gb' )
    IF ( ancil_dump_read%vegetation_props ) THEN
      CALL scatter_land_field(global_data_1d, t_growth_gb)
    END IF

    !Soil properties ancil namelist
    !Cases if nsoilt == 1, so it is OK to hardwire the 2nd dimension to 1

  CASE ( 'b      ' )
    IF ( ancil_dump_read%soil_props ) THEN
      IF ( l_tile_soil .AND. l_broadcast_soilt ) THEN
        DO m = 1,nsoilt
          DO n = 1,sm_levels
            CALL scatter_land_field(global_data_2d(:,n),                      &
                                    psparms%bexp_soilt(:,m,n))
          END DO
        END DO
      ELSE !Case if nsoilt == 1, so OK to hardwire the 2nd dimension to 1
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_2d(:,n),                        &
                                  psparms%bexp_soilt(:,1,n))
        END DO
      END IF
    END IF

  CASE ( 'b_soilt' )
    IF ( ancil_dump_read%soil_props ) THEN
      DO m = 1,nsoilt
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_3d(:,m,n),                      &
                                  psparms%bexp_soilt(:,m,n))
        END DO
      END DO
    END IF

  CASE ( 'sathh  ' )
    IF ( ancil_dump_read%soil_props ) THEN
      IF ( l_tile_soil .AND. l_broadcast_soilt ) THEN
        DO m = 1,nsoilt
          DO n = 1,sm_levels
            CALL scatter_land_field(global_data_2d(:,n),                      &
                                    psparms%sathh_soilt(:,m,n))
          END DO
        END DO
      ELSE !Case if nsoilt == 1, so OK to hardwire the 2nd dimension to 1
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_2d(:,n),                        &
                                  psparms%sathh_soilt(:,1,n))
        END DO
      END IF
    END IF

  CASE ( 'sathh_soilt' )
    IF ( ancil_dump_read%soil_props ) THEN
      DO m = 1,nsoilt
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_3d(:,m,n),                      &
                                  psparms%sathh_soilt(:,m,n))
        END DO
      END DO
    END IF

  CASE ( 'satcon ' )
    IF ( ancil_dump_read%soil_props ) THEN
      IF ( l_tile_soil .AND. l_broadcast_soilt ) THEN
        DO m = 1,nsoilt
          DO n = 1,sm_levels
            CALL scatter_land_field(global_data_2d(:,n),                      &
                                    psparms%satcon_soilt(:,m,n))
          END DO
        END DO
      ELSE !Case if nsoilt == 1, so OK to hardwire the 2nd dimension to 1
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_2d(:,n),                        &
                                  psparms%satcon_soilt(:,1,n))
        END DO
      END IF
    END IF

  CASE ( 'satcon_soilt' )
    IF ( ancil_dump_read%soil_props ) THEN
      DO m = 1,nsoilt
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_3d(:,m,n),                      &
                                  psparms%satcon_soilt(:,m,n))
        END DO
      END DO
    END IF

  CASE ( 'sm_sat ' )
    IF ( ancil_dump_read%soil_props ) THEN
      IF ( l_tile_soil .AND. l_broadcast_soilt ) THEN
        DO m = 1,nsoilt
          DO n = 1,sm_levels
            CALL scatter_land_field(global_data_2d(:,n),                      &
                                    psparms%smvcst_soilt(:,m,n))
          END DO
        END DO
      ELSE !Case if nsoilt == 1, so OK to hardwire the 2nd dimension to 1
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_2d(:,n),                        &
                                  psparms%smvcst_soilt(:,1,n))
        END DO
      END IF
    END IF

  CASE ( 'sm_sat_soilt' )
    IF ( ancil_dump_read%soil_props ) THEN
      DO m = 1,nsoilt
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_3d(:,m,n),                      &
                                  psparms%smvcst_soilt(:,m,n))
        END DO
      END DO
    END IF

  CASE ( 'sm_crit' )
    IF ( ancil_dump_read%soil_props ) THEN
      IF ( l_tile_soil .AND. l_broadcast_soilt ) THEN
        DO m = 1,nsoilt
          DO n = 1,sm_levels
            CALL scatter_land_field(global_data_2d(:,n),                      &
                                    psparms%smvccl_soilt(:,m,n))
          END DO
        END DO
      ELSE !Case if nsoilt == 1, so OK to hardwire the 2nd dimension to 1
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_2d(:,n),                        &
                                  psparms%smvccl_soilt(:,1,n))
        END DO
      END IF
    END IF

  CASE ( 'sm_crit_soilt' )
    IF ( ancil_dump_read%soil_props ) THEN
      DO m = 1,nsoilt
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_3d(:,m,n),                      &
                                  psparms%smvccl_soilt(:,m,n))
        END DO
      END DO
    END IF

  CASE ( 'sm_wilt' )
    IF ( ancil_dump_read%soil_props ) THEN
      IF ( l_tile_soil .AND. l_broadcast_soilt ) THEN
        DO m = 1,nsoilt
          DO n = 1,sm_levels
            CALL scatter_land_field(global_data_2d(:,n),                      &
                                    psparms%smvcwt_soilt(:,m,n))
          END DO
        END DO
      ELSE !Case if nsoilt == 1, so OK to hardwire the 2nd dimension to 1
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_2d(:,n),                        &
                                  psparms%smvcwt_soilt(:,1,n))
        END DO
      END IF
    END IF

  CASE ( 'sm_wilt_soilt' )
    IF ( ancil_dump_read%soil_props ) THEN
      DO m = 1,nsoilt
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_3d(:,m,n),                      &
                                  psparms%smvcwt_soilt(:,m,n))
        END DO
      END DO
    END IF

  CASE ( 'hcap   ' )
    IF ( ancil_dump_read%soil_props ) THEN
      IF ( l_tile_soil .AND. l_broadcast_soilt ) THEN
        DO m = 1,nsoilt
          DO n = 1,sm_levels
            CALL scatter_land_field(global_data_2d(:,n),                      &
                                    psparms%hcap_soilt(:,m,n))
          END DO
        END DO
      ELSE !Case if nsoilt == 1, so OK to hardwire the 2nd dimension to 1
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_2d(:,n),                        &
                                  psparms%hcap_soilt(:,1,n))
        END DO
      END IF
    END IF

  CASE ( 'hcap_soilt' )
    IF ( ancil_dump_read%soil_props ) THEN
      DO m = 1,nsoilt
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_3d(:,m,n),                      &
                                  psparms%hcap_soilt(:,m,n))
        END DO
      END DO
    END IF

  CASE ( 'hcon   ' )
    IF ( ancil_dump_read%soil_props ) THEN
      IF ( l_tile_soil .AND. l_broadcast_soilt ) THEN
        DO m = 1,nsoilt
          DO n = 1,sm_levels
            CALL scatter_land_field(global_data_2d(:,n),                      &
                                    psparms%hcon_soilt(:,m,n))
          END DO
        END DO
      ELSE !Case if nsoilt == 1, so OK to hardwire the 2nd dimension to 1
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_2d(:,n),                        &
                                  psparms%hcon_soilt(:,1,n))
        END DO
      END IF
    END IF

  CASE ( 'hcon_soilt' )
    IF ( ancil_dump_read%soil_props ) THEN
      DO m = 1,nsoilt
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_3d(:,m,n),                      &
                                  psparms%hcon_soilt(:,m,n))
        END DO
      END DO
    END IF

  CASE ( 'albsoil' )
    IF ( ancil_dump_read%soil_props ) THEN
      IF ( l_tile_soil .AND. l_broadcast_soilt ) THEN
        DO m = 1,nsoilt
          CALL scatter_land_field(global_data_1d, psparms%albsoil_soilt(:,m))
        END DO
      ELSE !Case if nsoilt == 1, so OK to hardwire the 2nd dimension to 1
        CALL scatter_land_field(global_data_1d, psparms%albsoil_soilt(:,1))
      END IF
    END IF

  CASE ( 'albsoil_soilt' )
    IF ( ancil_dump_read%soil_props ) THEN
      DO m = 1,nsoilt
        CALL scatter_land_field(global_data_2d(:,m), psparms%albsoil_soilt(:,m))
      END DO
    END IF

  CASE ( 'clay' )
    IF ( ancil_dump_read%soil_props ) THEN
      DO n = 1,dim_cslayer
        CALL scatter_land_field(global_data_2d(:,n),                          &
                                psparms%clay_soilt(:,1,n))
      END DO
    END IF

  CASE ( 'clay_soilt' )
    IF ( ancil_dump_read%soil_props ) THEN
      DO m = 1,nsoilt
        DO n = 1,dim_cslayer
          CALL scatter_land_field(global_data_3d(:,m,n),                      &
                                  psparms%clay_soilt(:,m,n))
        END DO
      END DO
    END IF

  CASE ( 'soil_ph' )
    IF ( ancil_dump_read%soil_props ) THEN
      DO n = 1,dim_cslayer
        CALL scatter_land_field(global_data_2d(:,n),                          &
                                psparms%soil_ph_soilt(:,1,n))
      END DO
    END IF

  CASE ( 'soil_ph_soilt' )
    IF ( ancil_dump_read%soil_props ) THEN
      DO m = 1,nsoilt
        DO n = 1,dim_cslayer
          CALL scatter_land_field(global_data_3d(:,m,n),                      &
                                  psparms%soil_ph_soilt(:,m,n))
        END DO
      END DO
    END IF

    !Topmodel ancil namelist
  CASE ( 'fexp   ' )
    IF ( ancil_dump_read%top ) THEN
      IF ( l_tile_soil .AND. l_broadcast_soilt ) THEN
        DO m = 1,nsoilt
          CALL scatter_land_field(global_data_1d, toppdm%fexp_soilt(:,m))
        END DO
      ELSE !Case if nsoilt == 1, so OK to hardwire the 2nd dimension to 1
        CALL scatter_land_field(global_data_1d, toppdm%fexp_soilt(:,1))
      END IF
    END IF

  CASE ( 'toppdm%fexp_soilt' )
    IF ( ancil_dump_read%top ) THEN
      DO m = 1,nsoilt
        CALL scatter_land_field(global_data_2d(:,m), toppdm%fexp_soilt(:,m))
      END DO
    END IF

    !Case if nsoilt == 1, so it is OK to hardwire the 2nd dimension to 1
  CASE ( 'ti_mean' )
    IF ( ancil_dump_read%top ) THEN
      IF ( l_tile_soil .AND. l_broadcast_soilt ) THEN
        DO m = 1,nsoilt
          CALL scatter_land_field(global_data_1d, toppdm%ti_mean_soilt(:,m))
        END DO
      ELSE !Case if nsoilt == 1, so OK to hardwire the 2nd dimension to 1
        CALL scatter_land_field(global_data_1d, toppdm%ti_mean_soilt(:,1))
      END IF
    END IF

  CASE ( 'toppdm%ti_mean_soilt' )
    IF ( ancil_dump_read%top ) THEN
      DO m = 1,nsoilt
        CALL scatter_land_field(global_data_2d(:,m), toppdm%ti_mean_soilt(:,m))
      END DO
    END IF

    !Case if nsoilt == 1, so it is OK to hardwire the 2nd dimension to 1
  CASE ( 'ti_sig ' )
    IF ( ancil_dump_read%top ) THEN
      IF ( l_tile_soil .AND. l_broadcast_soilt ) THEN
        DO m = 1,nsoilt
          CALL scatter_land_field(global_data_1d, toppdm%ti_sig_soilt(:,m))
        END DO
      ELSE !Case if nsoilt == 1, so OK to hardwire the 2nd dimension to 1
        CALL scatter_land_field(global_data_1d, toppdm%ti_sig_soilt(:,1))
      END IF
    END IF

  CASE ( 'toppdm%ti_sig_soilt' )
    IF ( ancil_dump_read%top ) THEN
      DO m = 1,nsoilt
        CALL scatter_land_field(global_data_2d(:,m), toppdm%ti_sig_soilt(:,m))
      END DO
    END IF

    !Agric ancil namelist
  CASE ( 'frac_agr' )
    IF ( ancil_dump_read%agric ) THEN
      CALL scatter_land_field(global_data_1d, trifctltype%frac_agr_gb)
    END IF

  CASE ( 'frac_past' )
    IF ( ancil_dump_read%agric ) THEN
      CALL scatter_land_field(global_data_1d, trif_vars%frac_past_gb)
    END IF

    !Crop props ancillaries namelist
  CASE ( 'cropsowdate       ' )
    IF ( ancil_dump_read%crop_props ) THEN
      DO n = 1,ncpft
        CALL scatter_land_field(global_data_2d(:,n),                          &
                                crop_vars%sow_date_cpft(:,n))
      END DO
    END IF

  CASE ( 'croplatestharvdate' )
    IF ( ancil_dump_read%crop_props ) THEN
      DO n = 1,ncpft
        CALL scatter_land_field(global_data_2d(:,n),                          &
                                crop_vars%latestharv_date_cpft(:,n))
      END DO
    END IF

  CASE ( 'cropttveg         ' )
    IF ( ancil_dump_read%crop_props ) THEN
      DO n = 1,ncpft
        CALL scatter_land_field(global_data_2d(:,n),                          &
                                crop_vars%tt_veg_cpft(:,n))
      END DO
    END IF

  CASE ( 'cropttrep         ' )
    IF ( ancil_dump_read%crop_props ) THEN
      DO n = 1,ncpft
        CALL scatter_land_field(global_data_2d(:,n),                          &
                                crop_vars%tt_rep_cpft(:,n))
      END DO
    END IF

    !Irrigation ancillaries namelist
  CASE ( 'frac_irrig' )
    IF ( ancil_dump_read%irrig ) THEN
      CALL scatter_land_field(global_data_1d,                                 &
                              crop_vars%frac_irr_all(:,1))
    END IF

  CASE ( 'frac_irr_all_tiles' )
    IF ( ancil_dump_read%irrig ) THEN
      IF ( frac_irrig_all_tiles_real > 0.5 ) THEN
        frac_irrig_all_tiles = .TRUE.
      ELSE
        frac_irrig_all_tiles = .FALSE.
      END IF
      CALL mpi_bcast(frac_irrig_all_tiles, 1, mpi_logical,                    &
           master_task_id, mpi_comm_world, error)
    END IF

  CASE ( 'irrfrac_irrtiles' )
    IF ( ancil_dump_read%irrig ) THEN
      CALL scatter_land_field(global_data_1d, crop_vars%irrfrac_irrtiles(:,1))
    END IF

  CASE ( 'set_irrfrac_on_irrtiles' )
    IF ( ancil_dump_read%irrig ) THEN
      IF ( set_irrfrac_irrtiles_real > 0.5 ) THEN
        set_irrfrac_on_irrtiles = .TRUE.
      ELSE
        set_irrfrac_on_irrtiles = .FALSE.
      END IF
      CALL mpi_bcast(set_irrfrac_on_irrtiles, 1, mpi_logical,                 &
           master_task_id, mpi_comm_world, error)
    END IF

  CASE ( 'irrtiles' )
    IF ( ancil_dump_read%irrig ) THEN
      CALL mpi_bcast(NINT(irrtiles_real), npft, mpi_integer,                  &
           master_task_id, mpi_comm_world, error)
    END IF

  CASE ( 'nirrtile' )
    IF ( ancil_dump_read%irrig ) THEN
      CALL mpi_bcast(NINT(nirrtile_real), 1, mpi_integer,                     &
           master_task_id, mpi_comm_world, error)
    END IF

    ! Water resources properties ancil namelist
  CASE ( 'conveyance_loss' )
    IF ( ancil_dump_read%water_resources_props ) THEN
      CALL scatter_land_field(global_data_1d, conveyance_loss)
    END IF

  CASE ( 'irrig_eff' )
    IF ( ancil_dump_read%water_resources_props ) THEN
      CALL scatter_land_field(global_data_1d, irrig_eff)
    END IF

  CASE ( 'sfc_water_frac' )
    IF ( ancil_dump_read%water_resources_props ) THEN
      CALL scatter_land_field(global_data_1d, sfc_water_frac)
    END IF

    !CO2 ancil namelist
  CASE ( 'co2_mmr' )
    IF ( ancil_dump_read%co2 ) THEN
      CALL mpi_bcast(co2_mmr, 1, mpi_real,                                    &
           master_task_id, mpi_comm_world, error)
    END IF

  CASE DEFAULT
    CALL log_fatal("read_dump",                                               &
                   "Unexpected variable in dump - " //                        &
                   TRIM(identifiers(i)))
  END SELECT

END DO

! We are done with the file
IF ( is_master_task() ) CALL file_close(FILE)

IF ( ALLOCATED(global_data_1d) ) DEALLOCATE(global_data_1d)
IF ( ALLOCATED(global_data_2d) ) DEALLOCATE(global_data_2d)
IF ( ALLOCATED(global_data_3d) ) DEALLOCATE(global_data_3d)
IF ( ALLOCATED(global_data_4d) ) DEALLOCATE(global_data_4d)

RETURN

END SUBROUTINE read_dump
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
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

SUBROUTINE get_dim_info( l_reading, identifier, ndims, dim_sizes, dim_names,  &
                         l_read_from_dump )

!-----------------------------------------------------------------------------
! Description:
!   Given a variable's identifier, return the number, sizes and names of its
!   dimensions. Used when reading or writing a dump file.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!-----------------------------------------------------------------------------

USE jules_fields_mod, ONLY: progs

USE ancil_info, ONLY:                                                         &
  dim_cs1, dim_soil_n_pool, nsurft, nsoilt, dim_cslayer

USE imogen_constants, ONLY:                                                   &
  n_olevs, nfarray

USE jules_rivers_mod, ONLY:                                                   &
  np_rivers

USE jules_soil_mod, ONLY:                                                     &
  sm_levels, ns_deep

USE jules_soil_biogeochem_mod, ONLY:                                          &
  dim_ch4layer

USE jules_snow_mod, ONLY:                                                     &
  nsmax

USE jules_surface_types_mod, ONLY:                                            &
  npft, ntype, ncpft

USE model_grid_mod, ONLY:                                                     &
  global_land_pts

USE string_utils_mod, ONLY:                                                   &
  to_string

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(IN).
!-----------------------------------------------------------------------------
LOGICAL, INTENT(IN) ::                                                        &
  l_reading
    ! Flag indicating where this routine has been called from.
    ! T means the call is from read_dump.
    ! F means the call is from write_dump.

CHARACTER(LEN=*), INTENT(IN) ::                                               &
  identifier  ! The identifier.

!-----------------------------------------------------------------------------
! Arguments with INTENT(OUT).
!-----------------------------------------------------------------------------
INTEGER, INTENT(OUT) ::                                                       &
  ndims,                                                                      &
    ! Number of dimensions for the variable.
  dim_sizes(:)
    ! Size of each dimension.

CHARACTER(LEN=*), INTENT(OUT) ::                                              &
  dim_names(:)
    ! Name of each dimension.

!-----------------------------------------------------------------------------
! Optional arguments.
!-----------------------------------------------------------------------------
LOGICAL, INTENT(OUT), OPTIONAL ::                                             &
  l_read_from_dump
    ! Flag passed to read_dump to indicate that a variable should be read from
    ! the dump file. This is TRUE for most variables, FALSE only for ancillary
    ! variables that were read from ancillary files (not the dump).

!-----------------------------------------------------------------------------
! If this is call from read_dump, check we have the optional argument.
! If it is present, set the default to TRUE.
!-----------------------------------------------------------------------------
IF ( l_reading ) THEN
  IF ( .NOT. PRESENT(l_read_from_dump) ) THEN
    CALL log_fatal("get_dim_info",                                            &
                   "l_read_from_dump must be present when called from " //    &
                   "read_dump.")
  ELSE
    l_read_from_dump = .TRUE.
  END IF
END IF

!-----------------------------------------------------------------------------

SELECT CASE ( identifier )

CASE ( 'gs', 'cv', 'frac_agr_prev', 'frac_past_prev',                         &
       'wood_prod_fast', 'wood_prod_med',  'wood_prod_slow',                  &
       'temp_max_00h_r', 'temp_ave_00h_r', 'prec_tot_00h_r',                  &
       'prec_tot_12h_r', 'rhum_min_00h_r', 'dewp_ave_00h_r',                  &
       'wind_ave_00h_r', 'temp_max_00h',   'temp_ave_00h',   'temp_ave_nday', &
       'temp_pnt_12h',   'prec_tot_00h',   'prec_tot_12h',                    &
       'rhum_min_00h',   'rhum_pnt_12h',   'dewp_ave_00h',                    &
       'wind_ave_00h',   'wind_pnt_12h',                                      &
       'fire_mcarthur_r_dr', 'fire_mcarthur_n_dr',                            &
       'fire_canadian_ffmc', 'fire_canadian_ffmc_mois',                       &
       'fire_canadian_dmc',  'fire_canadian_dc',                              &
       'fire_nesterov', 'lake_fetch_gb', 'lake_t_mean_gb',                    &
       'lake_t_mxl_gb', 'lake_h_mxl_gb', 'lake_t_ice_gb',                     &
       'lake_h_ice_gb', 'lake_shape_factor_gb', 'latitude', 'longitude')
  ndims = 1
  dim_names(1) = land_dim_name
  dim_sizes(1) = global_land_pts

CASE ( 'sthzw', 'zw' )
  ndims = 1
  dim_names(1) = land_dim_name
  dim_sizes(1) = global_land_pts

CASE ( 'sthzw_soilt', 'zw_soilt' )
  ndims = 2
  dim_names(1:ndims) = (/ land_dim_name, soilt_dim_name /)
  dim_sizes(1:ndims) = (/ global_land_pts, nsoilt /)

CASE ( 'canht', 'lai' )
  ndims = 2
  dim_names(1:ndims) = (/ land_dim_name, pft_dim_name /)
  dim_sizes(1:ndims) = (/ global_land_pts, npft /)

CASE ( 'cropdvi', 'croprootc', 'cropharvc', 'cropreservec',                   &
       'croplai', 'cropcanht' )
  ndims = 2
  dim_names(1:ndims) = (/ land_dim_name, cpft_dim_name /)
  dim_sizes(1:ndims) = (/ global_land_pts, ncpft /)

CASE ( 'cs', 'ns' )
  ndims = 3
  dim_names(1:ndims) = (/ land_dim_name, sc_layer_dim_name,                   &
                          sc_pool_dim_name /)
  dim_sizes(1:ndims) = (/ global_land_pts, dim_cslayer, dim_cs1 /)

CASE ( 'n_soil' )
  ndims = 3
  dim_names(1:ndims) = (/ land_dim_name, sc_layer_dim_name,                   &
                          soil_n_pool_dim_name /)
  dim_sizes(1:ndims) = (/ global_land_pts, dim_cslayer, dim_soil_n_pool /)

CASE ( 'cs_soilt', 'ns_soilt' )
  ndims = 4
  dim_names(1:ndims) = (/ land_dim_name, soilt_dim_name,                      &
                          sc_layer_dim_name, sc_pool_dim_name /)
  dim_sizes(1:ndims) = (/ global_land_pts, nsoilt, dim_cslayer, dim_cs1 /)

CASE ( 'n_soil_soilt' )
  ndims = 4
  dim_names(1:ndims) = (/ land_dim_name, soilt_dim_name,                      &
                          sc_layer_dim_name, soil_n_pool_dim_name /)
  dim_sizes(1:ndims) = (/ global_land_pts, nsoilt, dim_cslayer,               &
                          dim_soil_n_pool /)

CASE ( 'rgrainl', 'snow_ds', 'snow_ice', 'snow_liq', 'tsnow' )
  ndims = 3
  dim_names(1:ndims) = (/ land_dim_name, tile_dim_name, snow_dim_name /)
  dim_sizes(1:ndims) = (/ global_land_pts, nsurft, nsmax /)

CASE ( 'sthuf', 't_soil', 'sthu_irr' )
  ndims = 2
  dim_names(1:ndims) = (/ land_dim_name, soil_dim_name /)
  dim_sizes(1:ndims) = (/ global_land_pts, sm_levels /)

CASE ( 'sthuf_soilt', 't_soil_soilt', 'sthu_irr_soilt' )
  ndims = 3
  dim_names(1:ndims) = (/ land_dim_name, soilt_dim_name, soil_dim_name /)
  dim_sizes(1:ndims) = (/ global_land_pts, nsoilt, sm_levels /)

CASE ( 'n_inorg' )
  ndims = 2
  dim_names(1:ndims) = (/ land_dim_name, sc_layer_dim_name /)
  dim_sizes(1:ndims) = (/ global_land_pts, dim_cslayer /)

CASE ( 'n_inorg_soilt' )
  ndims = 3
  dim_names(1:ndims) = (/ land_dim_name, soilt_dim_name,                      &
                          sc_layer_dim_name /)
  dim_sizes(1:ndims) = (/ global_land_pts, nsoilt, dim_cslayer /)

CASE ( 'substr_ch4', 'mic_ch4', 'mic_act_ch4', 'acclim_ch4' )
  ndims = 2
  dim_names(1:ndims) = (/ land_dim_name, ch4layer_dim_name /)
  dim_sizes(1:ndims) = (/ global_land_pts, dim_ch4layer /)

CASE ( 'tsoil_deep' )
  ndims = 2
  dim_names(1:ndims) = (/ land_dim_name, bedrock_dim_name /)
  dim_sizes(1:ndims) = (/ global_land_pts, ns_deep /)

CASE ( 'canopy', 'nsnow', 'rgrain', 'rho_snow', 'snow_tile',                  &
       'snow_depth', 'snow_grnd', 'tstar_tile', 'tsurf_elev_surft' )
  ndims = 2
  dim_names(1:ndims) = (/ land_dim_name, tile_dim_name /)
  dim_sizes(1:ndims) = (/ global_land_pts, nsurft /)

  ! Cases for IMOGEN variables
CASE ( 'co2_ppmv', 'co2_change_ppmv', 'ch4_ppbv' )
  ! Scalar variables are represented by an array of dimension 1.
  ndims = 1
  dim_names(1) = scalar_dim_name
  dim_sizes(1) = 1

CASE ( 'dtemp_o' )
  ndims = 1
  dim_names(1) = nolevs_dim_name
  dim_sizes(1) = n_olevs

CASE ( 'fa_ocean' )
  ndims = 1
  dim_names(1) = nfarray_dim_name
  dim_sizes(1) = nfarray

CASE ( 'seed_rain' )
  ndims = 1
  dim_names(1) = seed_dim_name
  dim_sizes(1) = SIZE(progs%seed_rain)

  ! River routing variables.
CASE ( 'rivers_lat_rp', 'rivers_lon_rp',                                      &
       'rivers_sto_rp', 'rfm_surfstore_rp', 'rfm_substore_rp',                &
       'rfm_flowin_rp', 'rfm_bflowin_rp' )
  ndims = 1
  dim_names(1) = p_rivers_dim_name
  dim_sizes(1) = np_rivers

  !---------------------------------------------------------------------------
  ! Ancillary variables.
  ! If this is a call from write_dump (l_reading=F) or a call from
  ! read_dump that indicates the variable is to be read from the dump
  ! (l_reading=T and, e.g. ancil_dump_read%frac=T) we return the
  ! information about the dimensions.
  ! If this is call about a variable that was read from an ancillary
  ! file (l_reading=T and, e.g. ancil_dump_read%frac=F) we set
  ! l_read_from_dump=F to show that the field should not be read from the
  ! dump.
  ! The decision about which of these paths should be followed is made by the
  ! function need_dims.
  !---------------------------------------------------------------------------

  ! latlon ancil namelist
CASE ( 'grid_area' )
  IF ( need_dims( l_reading, ancil_dump_read%latlon ) ) THEN
    ndims = 1
    dim_names(1) = land_dim_name
    dim_sizes(1) = global_land_pts
  ELSE
    l_read_from_dump = .FALSE.
  END IF

  ! Frac ancil namelist
CASE ( 'frac' )
  IF ( need_dims( l_reading, ancil_dump_read%frac ) ) THEN
    ndims = 2
    dim_names(1:ndims) = (/ land_dim_name, type_dim_name /)
    dim_sizes(1:ndims) = (/ global_land_pts, ntype /)
  ELSE
    l_read_from_dump = .FALSE.
  END IF

  ! Vegetation properties ancil namelist
CASE ( 't_growth_gb' )
  IF ( need_dims( l_reading, ancil_dump_read%vegetation_props ) ) THEN
    ndims = 1
    dim_names(1) = land_dim_name
    dim_sizes(1) = global_land_pts
  ELSE
    l_read_from_dump = .FALSE.
  END IF

  ! Soil properties ancil namelist
CASE ( 'b      ', 'sathh  ', 'satcon ', 'sm_sat ', 'sm_crit',                 &
       'sm_wilt', 'hcap   ', 'hcon   ' )
  IF ( need_dims( l_reading, ancil_dump_read%soil_props ) ) THEN
    ndims = 2
    dim_names(1:ndims) = (/ land_dim_name, soil_dim_name /)
    dim_sizes(1:ndims) = (/ global_land_pts, sm_levels /)
  ELSE
    l_read_from_dump = .FALSE.
  END IF

CASE ( 'b_soilt', 'sathh_soilt', 'satcon_soilt', 'sm_sat_soilt',              &
       'sm_crit_soilt', 'sm_wilt_soilt', 'hcap_soilt', 'hcon_soilt' )
  IF ( need_dims( l_reading, ancil_dump_read%soil_props ) ) THEN
    ndims = 3
    dim_names(1:ndims) = (/ land_dim_name,soilt_dim_name,soil_dim_name /)
    dim_sizes(1:ndims) = (/ global_land_pts, nsoilt, sm_levels /)
  ELSE
    l_read_from_dump = .FALSE.
  END IF

CASE ( 'albsoil' )
  IF ( need_dims( l_reading, ancil_dump_read%soil_props ) ) THEN
    ndims = 1
    dim_names(1) = land_dim_name
    dim_sizes(1) = global_land_pts
  ELSE
    l_read_from_dump = .FALSE.
  END IF

CASE ( 'albsoil_soilt' )
  IF ( need_dims( l_reading, ancil_dump_read%soil_props ) ) THEN
    ndims = 2
    dim_names(1:ndims) = (/ land_dim_name, soilt_dim_name /)
    dim_sizes(1:ndims) = (/ global_land_pts, nsoilt /)
  ELSE
    l_read_from_dump = .FALSE.
  END IF

CASE ( 'clay', 'soil_ph' )
  IF ( need_dims( l_reading, ancil_dump_read%soil_props ) ) THEN
    ndims = 2
    dim_names(1:ndims) = (/ land_dim_name, sc_layer_dim_name /)
    dim_sizes(1:ndims) = (/ global_land_pts, dim_cslayer /)
  ELSE
    l_read_from_dump = .FALSE.
  END IF

CASE ( 'clay_soilt', 'soil_ph_soilt' )
  IF ( need_dims( l_reading, ancil_dump_read%soil_props ) ) THEN
    ndims = 3
    dim_names(1:ndims) = (/ land_dim_name, soilt_dim_name,                    &
                      sc_layer_dim_name /)
    dim_sizes(1:ndims) = (/ global_land_pts, nsoilt, dim_cslayer /)
  ELSE
    l_read_from_dump = .FALSE.
  END IF

  ! Topmodel ancil namelist
CASE ( 'fexp   ', 'ti_mean', 'ti_sig ' )
  IF ( need_dims( l_reading, ancil_dump_read%top ) ) THEN
    ndims = 1
    dim_names(1) = land_dim_name
    dim_sizes(1) = global_land_pts
  ELSE
    l_read_from_dump = .FALSE.
  END IF

CASE ( 'fexp_soilt', 'ti_mean_soilt', 'ti_sig_soilt' )
  IF ( need_dims( l_reading, ancil_dump_read%top ) ) THEN
    ndims = 2
    dim_names(1:ndims) = (/ land_dim_name, soilt_dim_name /)
    dim_sizes(1:ndims) = (/ global_land_pts, nsoilt /)
  ELSE
    l_read_from_dump = .FALSE.
  END IF

  ! Agric ancillaries namelist
CASE ( 'frac_agr', 'frac_past' )
  IF ( need_dims( l_reading, ancil_dump_read%agric ) ) THEN
    ndims = 1
    dim_names(1) = land_dim_name
    dim_sizes(1) = global_land_pts
  ELSE
    l_read_from_dump = .FALSE.
  END IF

  ! Crop props ancillaries namelist
CASE ( 'cropsowdate', 'cropttveg  ', 'cropttrep  ','croplatestharvdate' )
  IF ( need_dims( l_reading, ancil_dump_read%crop_props ) ) THEN
    ndims = 2
    dim_names(1:ndims) = (/ land_dim_name, cpft_dim_name /)
    dim_sizes(1:ndims) = (/ global_land_pts, ncpft /)
  ELSE
    l_read_from_dump = .FALSE.
  END IF

  ! Irrigation ancillaries namelist
CASE ( 'frac_irrig', 'irrfrac_irrtiles' )
  IF ( need_dims( l_reading, ancil_dump_read%irrig ) ) THEN
    ndims = 1
    dim_names(1) = land_dim_name
    dim_sizes(1) = global_land_pts
  ELSE
    l_read_from_dump = .FALSE.
  END IF

CASE ( 'frac_irr_all_tiles', 'nirrtile', 'set_irrfrac_on_irrtiles' )
  IF ( need_dims( l_reading, ancil_dump_read%irrig ) ) THEN
    ! Scalar variables are represented by an array of dimension 1.
    ndims = 1
    dim_names(1) = scalar_dim_name
    dim_sizes(1) = 1
  ELSE
    l_read_from_dump = .FALSE.
  END IF

CASE ( 'irrtiles' )
  IF ( need_dims( l_reading, ancil_dump_read%irrig ) ) THEN
    ndims = 1
    dim_names(1) = pft_dim_name
    dim_sizes(1) = npft
  ELSE
    l_read_from_dump = .FALSE.
  END IF

  ! Water resource properties ancil namelist
CASE ( 'conveyance_loss', 'irrig_eff', 'sfc_water_frac' )
  IF ( need_dims( l_reading, ancil_dump_read%water_resources_props ) ) THEN
    ndims = 1
    dim_names(1) = land_dim_name
    dim_sizes(1) = global_land_pts
  ELSE
    l_read_from_dump = .FALSE.
  END IF

  ! CO2 ancil namelist
CASE ( 'co2_mmr' )
  IF ( need_dims( l_reading, ancil_dump_read%co2 ) ) THEN
    ! Scalar variables are represented by an array of dimension 1.
    ndims = 1
    dim_names(1) = scalar_dim_name
    dim_sizes(1) = 1
  ELSE
    l_read_from_dump = .FALSE.
  END IF

  !FLake ancil namelist
CASE ( 'lake_depth_gb' )
  ndims = 1
  dim_names(1) = land_dim_name
  dim_sizes(1) = global_land_pts

CASE DEFAULT
  CALL log_fatal("get_dim_info",                                              &
                 "Unrecognised variable: " // TRIM( identifier ) //           &
                 " l_reading=" // to_string(l_reading) )

END SELECT

RETURN

END SUBROUTINE get_dim_info

!#############################################################################
!#############################################################################

FUNCTION need_dims( l_reading, l_read_ancil_from_dump ) RESULT( l_need_dims )

!-----------------------------------------------------------------------------
! Function to indicate whether the calling routine needs to get information
! about the dimensions of an ancillary variable during the reading or writing
! of a dump file.
!-----------------------------------------------------------------------------

IMPLICIT NONE

LOGICAL, INTENT(IN) ::                                                        &
  l_reading,                                                                  &
    ! T means we are reading a dump file.
    ! F means we are writing a dump file.
  l_read_ancil_from_dump
    ! T means this variable should be read from the dump (rather than having
    ! already been read from another file). Only used when l_reading=T.

! Function result.
LOGICAL ::                                                                    &
  l_need_dims
    ! T indicates that the calling routine should get information about the
    !   dimensions of the current ancillary variable.
    ! F indicates that no information is required (because this ancillary
    !   variable has already been read from a file).

!-----------------------------------------------------------------------------
! If this is a call while writing a dump (l_reading=F), or a call while
! reading a dump that indicates this is an ancillary variable that is to be
! read from the dump (l_reading=T and l_read_ancil_from_dump=T) we set
! l_need_dims=T to show that the program will need to get information about
! the variable's dimensions.
! If this is a call while reading a dump but for a variable that was
! previously read from an ancillary file (l_reading=T and
! l_read_ancil_from_dump=F) we set
! l_need_dims=F to show that there is no need to get information about
! the variable's dimensions.
!-----------------------------------------------------------------------------
IF ( ( l_reading .AND. l_read_ancil_from_dump ) .OR.  .NOT. l_reading ) THEN
  l_need_dims = .TRUE.
ELSE
  l_need_dims = .FALSE.
END IF

END FUNCTION need_dims


END MODULE dump_mod
