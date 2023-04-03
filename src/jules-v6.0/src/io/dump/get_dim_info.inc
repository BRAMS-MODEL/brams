#if !defined(UM_JULES)
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

#endif
