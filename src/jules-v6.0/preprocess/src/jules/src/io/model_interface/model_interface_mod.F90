! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE model_interface_mod

USE io_constants, ONLY: max_sdf_name_len, max_attr_val_len

USE logging_mod, ONLY: log_warn, log_fatal

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   This module provides the interface between the IO routines and the
!   model variables
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------------
! Length for identifiers of model variables
INTEGER, PARAMETER :: identifier_len = 29

! The name of each levels dimension in output files
CHARACTER(LEN=max_sdf_name_len), PARAMETER ::                                 &
  bl_level_dim_name_out = 'bllevel',                                          &
  pft_dim_name_out     = 'pft',                                               &
  cpft_dim_name_out    = 'cpft',                                              &
  nvg_dim_name_out     = 'nvg',                                               &
  type_dim_name_out    = 'type',                                              &
  tile_dim_name_out    = 'tile',                                              &
  soilt_dim_name_out   = 'soilt',                                             &
  snow_dim_name_out    = 'snow',                                              &
  soil_dim_name_out    = 'soil',                                              &
  scpool_dim_name_out  = 'scpool',                                            &
  soil_n_pool_dim_name_out  = 'snpool',                                       &
  sclayer_dim_name_out = 'sclayer',                                           &
  tracer_dim_name_out  = 'tracer',                                            &
  ch4layer_dim_name_out = 'ch4layer',                                         &
  bedrock_dim_name_out = 'bedrock'


!-----------------------------------------------------------------------------
! Module variables
!-----------------------------------------------------------------------------
! The name of each levels dimension in input files
CHARACTER(LEN=max_sdf_name_len) ::                                            &
  bl_level_dim_name = 'bllevel',                                              &
  pft_dim_name     = 'pft',                                                   &
  cpft_dim_name    = 'cpft',                                                  &
  nvg_dim_name     = 'nvg',                                                   &
  type_dim_name    = 'type',                                                  &
  tile_dim_name    = 'tile',                                                  &
  soilt_dim_name   = 'soilt',                                                 &
  snow_dim_name    = 'snow',                                                  &
  soil_dim_name    = 'soil',                                                  &
  scpool_dim_name  = 'scpool',                                                &
  soil_n_pool_dim_name  = 'snpool',                                           &
  sclayer_dim_name = 'sclayer',                                               &
  tracer_dim_name  = 'tracer',                                                &
  ch4layer_dim_name = 'ch4layer',                                             &
  bedrock_dim_name = 'bedrock'


! The size of each possible levels dimension.
! Initialising these to something "unphysical" means that in future we can
! recognise if any dimension that is added to the code has not been set
! correctly (i.e. it will retain this obviously-wrong initial value).
INTEGER ::                                                                    &
  bl_level_dim_size    = -1,                                                  &
  pft_dim_size         = -1,                                                  &
  cpft_dim_size        = -1,                                                  &
  nvg_dim_size         = -1,                                                  &
  type_dim_size        = -1,                                                  &
  tile_dim_size        = -1,                                                  &
  soilt_dim_size       = -1,                                                  &
  snow_dim_size        = -1,                                                  &
  soil_dim_size        = -1,                                                  &
  scpool_dim_size      = -1,                                                  &
  soil_n_pool_dim_size = -1,                                                  &
  sclayer_dim_size     = -1,                                                  &
  tracer_dim_size      = -1,                                                  &
  ch4layer_dim_size    = -1,                                                  &
  bedrock_dim_size     = -1

! Due to how CABEL sets variables differently to JULES, these numbers
! are currently necassary.
INTEGER ::                                                                    &
  cable_tile_dim_size  = 17,                                                  &
  cable_soil_dim_size  = 6,                                                   &
  cable_snow_dim_size  = 3


!-----------------------------------------------------------------------------
! Information about the actual variables available for output
!-----------------------------------------------------------------------------
! Constants for the different 'types' of variable
INTEGER, PARAMETER ::                                                         &
  var_type_bl_level       = 0,                                                &
  var_type_surface        = var_type_bl_level+1,                              &
  var_type_pft            = var_type_surface+1,                               &
  var_type_cpft           = var_type_pft+1,                                   &
  var_type_nvg            = var_type_cpft+1,                                  &
  var_type_type           = var_type_nvg+1,                                   &
  var_type_surft          = var_type_type+1,                                  &
  var_type_soilt          = var_type_surft+1,                                 &
  var_type_soilt_soil     = var_type_soilt+1,                                 &
  var_type_snow           = var_type_soilt_soil+1,                            &
  var_type_soil           = var_type_snow+1,                                  &
  var_type_scpool         = var_type_soil+1,                                  &
  var_type_soil_n_pool    = var_type_scpool+1,                                &
  var_type_bedrock        = var_type_soil_n_pool+1,                           &
  var_type_rp             = var_type_bedrock+1,                               &
  var_type_sclayer        = var_type_rp+1,                                    &
  var_type_soilt_sclayer  = var_type_sclayer+1,                               &
  var_type_soilt_sclayer_scpool = var_type_soilt_sclayer+1,                   &
  var_type_tracer         = var_type_soilt_sclayer_scpool+1,                  &
  var_type_ch4layer       = var_type_tracer+1,                                &
  var_type_pft_sclayer    = var_type_ch4layer+1,                              &
  var_type_cable_soil     = var_type_pft_sclayer+1,                           &
  var_type_cable_snow     = var_type_cable_soil+1,                            &
  var_type_cable1L_snow   = var_type_cable_snow + 1

! Derived type to contain metadata about model variables
TYPE var_metadata

  CHARACTER(LEN=identifier_len) :: identifier
                            ! The string identifier of the variable

  INTEGER :: var_type  ! The type of the variable - must be one of the above

  CHARACTER(LEN=max_attr_val_len) :: long_name
                            ! The value of the long_name attribute
  CHARACTER(LEN=max_attr_val_len) :: units
                            ! The value of the units attribute

END TYPE var_metadata

! Array holding the metadata for all model variables that we can use for input
! or output. The CABLE land surface model adds 10 prognostics for tiled
! soil/snow prognostics
! or output
INTEGER, PARAMETER :: n_vars = 604
TYPE(var_metadata) :: metadata(n_vars)

! Include the metadata DATA statement
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!-----------------------------------------------------------------------------
! DATA statements for variable metadata
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Metadata for latitude
!-----------------------------------------------------------------------------
DATA metadata(1) / var_metadata(                                              &
! String identifier
    'latitude',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox latitude",                                                       &
! Units
    "degrees"                                                                 &
  ) /
!-----------------------------------------------------------------------------
! Metadata for longitude
!-----------------------------------------------------------------------------
DATA metadata(2) / var_metadata(                                              &
! String identifier
    'longitude',                                                              &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox longitude",                                                      &
! Units
    "degrees"                                                                 &
  ) /
!-----------------------------------------------------------------------------
! Metadata for land_fraction
!-----------------------------------------------------------------------------
DATA metadata(3) / var_metadata(                                              &
! String identifier
    'land_fraction',                                                          &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for surf_hgt
!-----------------------------------------------------------------------------
DATA metadata(4) / var_metadata(                                              &
! String identifier
    'surf_hgt',                                                               &
! Variable type
    var_type_surft,                                                           &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for z_land_land
!-----------------------------------------------------------------------------
DATA metadata(5) / var_metadata(                                              &
! String identifier
    'z_land_land',                                                            &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for frac
!-----------------------------------------------------------------------------
DATA metadata(6) / var_metadata(                                              &
! String identifier
    'frac',                                                                   &
! Variable type
    var_type_type,                                                            &
! Long name
    "Fractional cover of each surface type",                                  &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for albsoil
!-----------------------------------------------------------------------------
DATA metadata(7) / var_metadata(                                              &
! String identifier
    'albsoil',                                                                &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for clay
!-----------------------------------------------------------------------------
DATA metadata(8) / var_metadata(                                              &
! String identifier
    'clay',                                                                   &
! Variable type
    var_type_sclayer,                                                         &
! Long name
    "Soil clay content",                                                      &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for albobs_sw
!-----------------------------------------------------------------------------
DATA metadata(9) / var_metadata(                                              &
! String identifier
    'albobs_sw',                                                              &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for albobs_vis
!-----------------------------------------------------------------------------
DATA metadata(10) / var_metadata(                                             &
! String identifier
    'albobs_vis',                                                             &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for albobs_nir
!-----------------------------------------------------------------------------
DATA metadata(11) / var_metadata(                                             &
! String identifier
    'albobs_nir',                                                             &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for b
!-----------------------------------------------------------------------------
DATA metadata(12) / var_metadata(                                             &
! String identifier
    'b',                                                                      &
! Variable type
    var_type_soil,                                                            &
! Long name
    "Gridbox Brooks-Corey exponent for each soil layer",                      &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for b_const_z
!-----------------------------------------------------------------------------
DATA metadata(13) / var_metadata(                                             &
! String identifier
    'b_const_z',                                                              &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sathh
!-----------------------------------------------------------------------------
DATA metadata(14) / var_metadata(                                             &
! String identifier
    'sathh',                                                                  &
! Variable type
    var_type_soil,                                                            &
! Long name
    "Gridbox saturated soil water pressure for each soil layer",              &
! Units
    "m"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sathh_const_z
!-----------------------------------------------------------------------------
DATA metadata(15) / var_metadata(                                             &
! String identifier
    'sathh_const_z',                                                          &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for satcon
!-----------------------------------------------------------------------------
DATA metadata(16) / var_metadata(                                             &
! String identifier
    'satcon',                                                                 &
! Variable type
    var_type_soil,                                                            &
! Long name
    "Gridbox saturated hydraulic conductivity for each soil layer",           &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for satcon_const_z
!-----------------------------------------------------------------------------
DATA metadata(17) / var_metadata(                                             &
! String identifier
    'satcon_const_z',                                                         &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sm_sat
!-----------------------------------------------------------------------------
DATA metadata(18) / var_metadata(                                             &
! String identifier
    'sm_sat',                                                                 &
! Variable type
    var_type_soil,                                                            &
! Long name
    "Gridbox volumetric moisture content at saturation for each soil layer",  &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sm_sat_const_z
!-----------------------------------------------------------------------------
DATA metadata(19) / var_metadata(                                             &
! String identifier
    'sm_sat_const_z',                                                         &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sm_crit
!-----------------------------------------------------------------------------
DATA metadata(20) / var_metadata(                                             &
! String identifier
    'sm_crit',                                                                &
! Variable type
    var_type_soil,                                                            &
! Long name
    "Gridbox volumetric moisture content at critical point for each soil layer", &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sm_crit_const_z
!-----------------------------------------------------------------------------
DATA metadata(21) / var_metadata(                                             &
! String identifier
    'sm_crit_const_z',                                                        &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sm_wilt
!-----------------------------------------------------------------------------
DATA metadata(22) / var_metadata(                                             &
! String identifier
    'sm_wilt',                                                                &
! Variable type
    var_type_soil,                                                            &
! Long name
    "Gridbox  volumetric moisture content at wilting point for each soil layer", &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sm_wilt_const_z
!-----------------------------------------------------------------------------
DATA metadata(23) / var_metadata(                                             &
! String identifier
    'sm_wilt_const_z',                                                        &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for hcap
!-----------------------------------------------------------------------------
DATA metadata(24) / var_metadata(                                             &
! String identifier
    'hcap',                                                                   &
! Variable type
    var_type_soil,                                                            &
! Long name
    "Gridbox dry soil heat capacity for each soil layer",                     &
! Units
    "J K-1 m-3"                                                               &
  ) /
!-----------------------------------------------------------------------------
! Metadata for hcap_const_z
!-----------------------------------------------------------------------------
DATA metadata(25) / var_metadata(                                             &
! String identifier
    'hcap_const_z',                                                           &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for hcon
!-----------------------------------------------------------------------------
DATA metadata(26) / var_metadata(                                             &
! String identifier
    'hcon',                                                                   &
! Variable type
    var_type_soil,                                                            &
! Long name
    "Gridbox dry soil thermal conductivity for each soil layer",              &
! Units
    "W m-1 K-1"                                                               &
  ) /
!-----------------------------------------------------------------------------
! Metadata for hcon_const_z
!-----------------------------------------------------------------------------
DATA metadata(27) / var_metadata(                                             &
! String identifier
    'hcon_const_z',                                                           &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fexp
!-----------------------------------------------------------------------------
DATA metadata(28) / var_metadata(                                             &
! String identifier
    'fexp',                                                                   &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for ti_mean
!-----------------------------------------------------------------------------
DATA metadata(29) / var_metadata(                                             &
! String identifier
    'ti_mean',                                                                &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for ti_sig
!-----------------------------------------------------------------------------
DATA metadata(30) / var_metadata(                                             &
! String identifier
    'ti_sig',                                                                 &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for wrr
!-----------------------------------------------------------------------------
DATA metadata(31) / var_metadata(                                             &
! String identifier
    'wrr',                                                                    &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for hwr
!-----------------------------------------------------------------------------
DATA metadata(32) / var_metadata(                                             &
! String identifier
    'hwr',                                                                    &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for hgt
!-----------------------------------------------------------------------------
DATA metadata(33) / var_metadata(                                             &
! String identifier
    'hgt',                                                                    &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for ztm
!-----------------------------------------------------------------------------
DATA metadata(34) / var_metadata(                                             &
! String identifier
    'ztm',                                                                    &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for disp
!-----------------------------------------------------------------------------
DATA metadata(35) / var_metadata(                                             &
! String identifier
    'disp',                                                                   &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for albwl
!-----------------------------------------------------------------------------
DATA metadata(36) / var_metadata(                                             &
! String identifier
    'albwl',                                                                  &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for albrd
!-----------------------------------------------------------------------------
DATA metadata(37) / var_metadata(                                             &
! String identifier
    'albrd',                                                                  &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for emisw
!-----------------------------------------------------------------------------
DATA metadata(38) / var_metadata(                                             &
! String identifier
    'emisw',                                                                  &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for emisr
!-----------------------------------------------------------------------------
DATA metadata(39) / var_metadata(                                             &
! String identifier
    'emisr',                                                                  &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for canopy
!-----------------------------------------------------------------------------
DATA metadata(40) / var_metadata(                                             &
! String identifier
    'canopy',                                                                 &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile surface/canopy water for snow-free land tiles",                     &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cs
!-----------------------------------------------------------------------------
DATA metadata(41) / var_metadata(                                             &
! String identifier
    'cs',                                                                     &
! Variable type
    var_type_scpool,                                                          &
! Long name is dynamic - set in get_var_attrs
    "",                                                                       &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for ns
!-----------------------------------------------------------------------------
DATA metadata(42) / var_metadata(                                             &
! String identifier
    'ns',                                                                     &
! Variable type
    var_type_scpool,                                                          &
! Long name
    "Gridbox soil organic Nitrogen in each pool (DPM,RPM,bio,hum)",           &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for immob_n
!-----------------------------------------------------------------------------
DATA metadata(43) / var_metadata(                                             &
! String identifier
    'immob_n',                                                                &
! Variable type
    var_type_scpool,                                                          &
! Long name
    "Gridbox soil Nitrogen immobilisation in each pool (DPM,RPM,bio,hum)",    &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for immob_n_pot
!-----------------------------------------------------------------------------
DATA metadata(44) / var_metadata(                                             &
! String identifier
    'immob_n_pot',                                                            &
! Variable type
    var_type_scpool,                                                          &
! Long name
    "Gridbox soil potential Nitrogen immobilisation in each pool (DPM,RPM,bio,hum)", &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for minl_n
!-----------------------------------------------------------------------------
DATA metadata(45) / var_metadata(                                             &
! String identifier
    'minl_n',                                                                 &
! Variable type
    var_type_scpool,                                                          &
! Long name
    "Gridbox soil Nitrogen mineralisation in each pool (DPM,RPM,bio,hum)",    &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for minl_n_pot
!-----------------------------------------------------------------------------
DATA metadata(46) / var_metadata(                                             &
! String identifier
    'minl_n_pot',                                                             &
! Variable type
    var_type_scpool,                                                          &
! Long name
    "Gridbox soil potential Nitrogen mineralisation in each pool (DPM,RPM,bio,hum)", &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for resp_s_diag
!-----------------------------------------------------------------------------
DATA metadata(47) / var_metadata(                                             &
! String identifier
    'resp_s_diag',                                                            &
! Variable type
    var_type_scpool,                                                          &
! Long name
    "Gridbox soil respiration carbon flux in each pool (DPM,RPM,bio,hum)",    &
! Units
    "kg m-2 per 360days"                                                      &
  ) /

!-----------------------------------------------------------------------------
! Metadata for resp_s_pot_diag
!-----------------------------------------------------------------------------
DATA metadata(48) / var_metadata(                                             &
! String identifier
    'resp_s_pot_diag',                                                        &
! Variable type
    var_type_scpool,                                                          &
! Long name
    "Gridbox potential soil respiration carbon flux in each pool (DPM,RPM,bio,hum)", &
! Units
    "kg m-2 per 360days"                                                      &
  ) /

!-----------------------------------------------------------------------------
! Metadata for soil_cn
!-----------------------------------------------------------------------------
DATA metadata(49) / var_metadata(                                             &
! String identifier
    'soil_CN',                                                                &
! Variable type
    var_type_scpool,                                                          &
! Long name
    "Gridbox Soil C:N in each pool (DPM,RPM,bio,hum)",                        &
! Units
    ":"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for soil_cn_gb
!-----------------------------------------------------------------------------
DATA metadata(50) / var_metadata(                                             &
! String identifier
    'soil_CN_gb',                                                             &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox Soil C:N",                                                       &
! Units
    ":"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for resp_s_pot_diag_gb
!-----------------------------------------------------------------------------
DATA metadata(51) / var_metadata(                                             &
! String identifier
    'resp_s_pot_diag_gb',                                                     &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean potential soil respiration",                                &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for resp_s_diag_gb
!-----------------------------------------------------------------------------
DATA metadata(52) / var_metadata(                                             &
! String identifier
    'resp_s_diag_gb',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean soil respiration",                                          &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for minl_n_pot_gb
!-----------------------------------------------------------------------------
DATA metadata(53) / var_metadata(                                             &
! String identifier
    'minl_n_pot_gb',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean potential soil N mineralisation",                           &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for minl_n_gb
!-----------------------------------------------------------------------------
DATA metadata(54) / var_metadata(                                             &
! String identifier
    'minl_n_gb',                                                              &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean potential soil N mineralisation",                           &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for immob_pot_gb
!-----------------------------------------------------------------------------
DATA metadata(55) / var_metadata(                                             &
! String identifier
    'immob_n_pot_gb',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean potential soil N immobilisation",                           &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for immob_n_gb
!-----------------------------------------------------------------------------
DATA metadata(56) / var_metadata(                                             &
! String identifier
    'immob_n_gb',                                                             &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean potential soil N immobilisation",                           &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for ns_gb
!-----------------------------------------------------------------------------
DATA metadata(57) / var_metadata(                                             &
! String identifier
    'ns_gb',                                                                  &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox soil organic Nitrogen (total)",                                  &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fN
!-----------------------------------------------------------------------------
DATA metadata(58) / var_metadata(                                             &
! String identifier
    'fN',                                                                     &
! Variable type
    var_type_sclayer,                                                         &
! Long name
    "Gridbox Nitrogen availability rate modifier",                            &
! Units
    ":"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for gs
!-----------------------------------------------------------------------------
DATA metadata(59) / var_metadata(                                             &
! String identifier
    'gs',                                                                     &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox surface conductance to evaporation",                             &
! Units
    "m s-1"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snow_tile
!-----------------------------------------------------------------------------
DATA metadata(60) / var_metadata(                                             &
! String identifier
    'snow_tile',                                                              &
! Variable type
    var_type_surft,                                                           &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sthuf
!-----------------------------------------------------------------------------
DATA metadata(61) / var_metadata(                                             &
! String identifier
    'sthuf',                                                                  &
! Variable type
    var_type_soil,                                                            &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for t_soil
!-----------------------------------------------------------------------------
DATA metadata(62) / var_metadata(                                             &
! String identifier
    't_soil',                                                                 &
! Variable type
    var_type_soil,                                                            &
! Long name
    "Gridbox sub-surface temperature of each layer",                          &
! Units
    "K"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for tsoil_deep
!-----------------------------------------------------------------------------
DATA metadata(63) / var_metadata(                                             &
! String identifier
    'tsoil_deep',                                                             &
! Variable type
    var_type_bedrock,                                                         &
! Long name
    "Temperature in bedrock layers",                                          &
! Units
    "K"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for tstar_tile
!-----------------------------------------------------------------------------
DATA metadata(64) / var_metadata(                                             &
! String identifier
    'tstar_tile',                                                             &
! Variable type
    var_type_surft,                                                           &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lai
!-----------------------------------------------------------------------------
DATA metadata(65) / var_metadata(                                             &
! String identifier
    'lai',                                                                    &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT leaf area index",                                                    &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for canht
!-----------------------------------------------------------------------------
DATA metadata(66) / var_metadata(                                             &
! String identifier
    'canht',                                                                  &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT canopy height",                                                      &
! Units
    "m"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sthzw
!-----------------------------------------------------------------------------
DATA metadata(67) / var_metadata(                                             &
! String identifier
    'sthzw',                                                                  &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Soil wetness in deep LSH/TOPMODEL layer",                                &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for zw
!-----------------------------------------------------------------------------
DATA metadata(68) / var_metadata(                                             &
! String identifier
    'zw',                                                                     &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean depth to water table",                                      &
! Units
    "m"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for rgrain
!-----------------------------------------------------------------------------
DATA metadata(69) / var_metadata(                                             &
! String identifier
    'rgrain',                                                                 &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile snow surface grain size",                                           &
! Units
    "microns"                                                                 &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cv
!-----------------------------------------------------------------------------
DATA metadata(70) / var_metadata(                                             &
! String identifier
    'cv',                                                                     &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean vegetation carbon at end of model timestep.",               &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for rho_snow
!-----------------------------------------------------------------------------
DATA metadata(71) / var_metadata(                                             &
! String identifier
    'rho_snow',                                                               &
! Variable type
    var_type_surft,                                                           &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snow_depth
!-----------------------------------------------------------------------------
DATA metadata(72) / var_metadata(                                             &
! String identifier
    'snow_depth',                                                             &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile snow depth (on ground)",                                            &
! Units
    "m"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snow_grnd
!-----------------------------------------------------------------------------
DATA metadata(73) / var_metadata(                                             &
! String identifier
    'snow_grnd',                                                              &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile snow on ground below canopy (snow_grnd)",                           &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for nsnow
!-----------------------------------------------------------------------------
DATA metadata(74) / var_metadata(                                             &
! String identifier
    'nsnow',                                                                  &
! Variable type
    var_type_surft,                                                           &
! Long name
   "Number of snow layers",                                                   &
! Units
   "1"                                                                        &
 ) /
!-----------------------------------------------------------------------------
! Metadata for snow_ds
!-----------------------------------------------------------------------------
DATA metadata(75) / var_metadata(                                             &
! String identifier
    'snow_ds',                                                                &
! Variable type
    var_type_snow,                                                            &
! Long name
    "Depth of snow layers for each tile",                                     &
! Units
    "m"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snow_ice
!-----------------------------------------------------------------------------
DATA metadata(76) / var_metadata(                                             &
! String identifier
    'snow_ice',                                                               &
! Variable type
    var_type_snow,                                                            &
! Long name
    "Ice mass in snow layers for each tile",                                  &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snow_liq
!-----------------------------------------------------------------------------
DATA metadata(77) / var_metadata(                                             &
! String identifier
    'snow_liq',                                                               &
! Variable type
    var_type_snow,                                                            &
! Long name
    "Liquid mass in snow layers for each tile",                               &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for tsnow
!-----------------------------------------------------------------------------
DATA metadata(78) / var_metadata(                                             &
! String identifier
    'tsnow',                                                                  &
! Variable type
    var_type_snow,                                                            &
! Long name
    "Temperature of snow layers for each tile",                               &
! Units
    "K"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for rgrainl
!-----------------------------------------------------------------------------
DATA metadata(79) / var_metadata(                                             &
! String identifier
    'rgrainl',                                                                &
! Variable type
    var_type_snow,                                                            &
! Long name
    "Grain size in snow layers for each tile",                                &
! Units
    "microns"                                                                 &
  ) /
!-----------------------------------------------------------------------------
! Metadata for pstar
!-----------------------------------------------------------------------------
DATA metadata(80) / var_metadata(                                             &
! String identifier
    'pstar',                                                                  &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox surface pressure",                                               &
! Units
    "Pa"                                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for q
!-----------------------------------------------------------------------------
DATA metadata(81) / var_metadata(                                             &
! String identifier
    'q',                                                                      &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for t
!-----------------------------------------------------------------------------
DATA metadata(82) / var_metadata(                                             &
! String identifier
    't',                                                                      &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for rad_net
!-----------------------------------------------------------------------------
DATA metadata(83) / var_metadata(                                             &
! String identifier
    'rad_net',                                                                &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Surface net radiation of land points",                                   &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lw_net
!-----------------------------------------------------------------------------
DATA metadata(84) / var_metadata(                                             &
! String identifier
    'lw_net',                                                                 &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox net downward longwave radiation at surface",                     &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sw_net
!-----------------------------------------------------------------------------
DATA metadata(85) / var_metadata(                                             &
! String identifier
    'sw_net',                                                                 &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox net downward shortwave radiation at surface",                    &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lw_down
!-----------------------------------------------------------------------------
DATA metadata(86) / var_metadata(                                             &
! String identifier
    'lw_down',                                                                &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox surface downward LW radiation",                                  &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sw_down
!-----------------------------------------------------------------------------
DATA metadata(87) / var_metadata(                                             &
! String identifier
    'sw_down',                                                                &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox surface downward SW radiation",                                  &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for diff_rad
!-----------------------------------------------------------------------------
DATA metadata(88) / var_metadata(                                             &
! String identifier
    'diff_rad',                                                               &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for precip
!-----------------------------------------------------------------------------
DATA metadata(89) / var_metadata(                                             &
! String identifier
    'precip',                                                                 &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox precipitation rate",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for tot_rain
!-----------------------------------------------------------------------------
DATA metadata(90) / var_metadata(                                             &
! String identifier
    'tot_rain',                                                               &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for tot_snow
!-----------------------------------------------------------------------------
DATA metadata(91) / var_metadata(                                             &
! String identifier
    'tot_snow',                                                               &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for con_rain
!-----------------------------------------------------------------------------
DATA metadata(92) / var_metadata(                                             &
! String identifier
    'con_rain',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox convective rainfall",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for ls_rain
!-----------------------------------------------------------------------------
DATA metadata(93) / var_metadata(                                             &
! String identifier
    'ls_rain',                                                                &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox large-scale rainfall",                                           &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for con_snow
!-----------------------------------------------------------------------------
DATA metadata(94) / var_metadata(                                             &
! String identifier
    'con_snow',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox convective snowfall",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for ls_snow
!-----------------------------------------------------------------------------
DATA metadata(95) / var_metadata(                                             &
! String identifier
    'ls_snow',                                                                &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox large-scale snowfall",                                           &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for wind
!-----------------------------------------------------------------------------
DATA metadata(96) / var_metadata(                                             &
! String identifier
    'wind',                                                                   &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox wind speed",                                                     &
! Units
    "m s-1"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for u
!-----------------------------------------------------------------------------
DATA metadata(97) / var_metadata(                                             &
! String identifier
    'u',                                                                      &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for v
!-----------------------------------------------------------------------------
DATA metadata(98) / var_metadata(                                             &
! String identifier
    'v',                                                                      &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for z1_tq_in
!-----------------------------------------------------------------------------
DATA metadata(99) / var_metadata(                                             &
! String identifier
    'z1_tq_in',                                                               &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for diurnal temperature range
!-----------------------------------------------------------------------------
DATA metadata(100) / var_metadata(                                            &
! String identifier
    'dt_range',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Diurnal temperature range",                                              &
! Units
    "K"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for albedo_land
!-----------------------------------------------------------------------------
DATA metadata(101) / var_metadata(                                            &
! String identifier
    'albedo_land',                                                            &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox albedo (as used for net shortwave calculation)",                 &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for NDVI_land
!-----------------------------------------------------------------------------
DATA metadata(102) / var_metadata(                                            &
! String identifier
    'NDVI_land',                                                              &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox NDVI (using sum of direct and diffuse for (NIR-VIS)/(NIR+VIS))", &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for canopy_gb
!-----------------------------------------------------------------------------
DATA metadata(103) / var_metadata(                                            &
! String identifier
    'canopy_gb',                                                              &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox canopy water content",                                           &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cs_gb
!-----------------------------------------------------------------------------
DATA metadata(104) / var_metadata(                                            &
! String identifier
    'cs_gb',                                                                  &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox soil carbon (total)",                                            &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for depth_frozen
!-----------------------------------------------------------------------------
DATA metadata(105) / var_metadata(                                            &
! String identifier
    'depth_frozen',                                                           &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox depth of frozen ground at surface (temperature definition)",     &
! Units
    "m"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for depth_frozen_sthf
!-----------------------------------------------------------------------------
DATA metadata(106) / var_metadata(                                            &
! String identifier
    'depth_frozen_sthf',                                                      &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox depth of frozen ground at surface (soil moisture definition)",   &
! Units
    "m"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for depth_unfrozen
!-----------------------------------------------------------------------------
DATA metadata(107) / var_metadata(                                            &
! String identifier
    'depth_unfrozen',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox depth of unfrozen ground at surface (temeprature definition)",   &
! Units
    "m"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for depth_unfrozen_sthf
!-----------------------------------------------------------------------------
DATA metadata(108) / var_metadata(                                            &
! String identifier
    'depth_unfrozen_sthf',                                                    &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox depth of unfrozen ground at surface (soil moisture definition)", &
! Units
    "m"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for drain
!-----------------------------------------------------------------------------
DATA metadata(109) / var_metadata(                                            &
! String identifier
    'drain',                                                                  &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Drainage from bottom (nshyd) soil layer",                                &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for elake
!-----------------------------------------------------------------------------
DATA metadata(110) / var_metadata(                                            &
! String identifier
    'elake',                                                                  &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean evaporation from lakes",                                    &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for emis_gb
!-----------------------------------------------------------------------------
DATA metadata(111) / var_metadata(                                            &
! String identifier
    'emis_gb',                                                                &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox emissivity",                                                     &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fch4_wetl
!-----------------------------------------------------------------------------
DATA metadata(112) / var_metadata(                                            &
! String identifier
    'fch4_wetl',                                                              &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Scaled methane flux from wetland fraction (for use in atmos chem)",      &
! Units
    "10^-9 kg m-2 s-1"                                                        &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fch4_wetl_cs
!-----------------------------------------------------------------------------
DATA metadata(113) / var_metadata(                                            &
! String identifier
    'fch4_wetl_cs',                                                           &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Scaled methane flux from wetland fraction (cs substrate)",               &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fch4_wetl_npp
!-----------------------------------------------------------------------------
DATA metadata(114) / var_metadata(                                            &
! String identifier
    'fch4_wetl_npp',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Scaled methane flux from wetland fraction (npp substrate)",              &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fch4_wetl_resps
!-----------------------------------------------------------------------------
DATA metadata(115) / var_metadata(                                            &
! String identifier
    'fch4_wetl_resps',                                                        &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Scaled methane flux from wetland fraction (soil resp substrate)",        &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fsat
!-----------------------------------------------------------------------------
DATA metadata(116) / var_metadata(                                            &
! String identifier
    'fsat',                                                                   &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Surface saturated fraction",                                             &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fsmc_gb
!-----------------------------------------------------------------------------
DATA metadata(117) / var_metadata(                                            &
! String identifier
    'fsmc_gb',                                                                &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox soil moisture availability factor (beta)",                       &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fwetl
!-----------------------------------------------------------------------------
DATA metadata(118) / var_metadata(                                            &
! String identifier
    'fwetl',                                                                  &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Wetland fraction at end of model timestep.",                             &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for gpp_gb
!-----------------------------------------------------------------------------
DATA metadata(119) / var_metadata(                                            &
! String identifier
    'gpp_gb',                                                                 &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox gross primary productivity",                                     &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for hf_snow_melt
!-----------------------------------------------------------------------------
DATA metadata(120) / var_metadata(                                            &
! String identifier
    'hf_snow_melt',                                                           &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox snowmelt heat flux",                                             &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for land_index
!-----------------------------------------------------------------------------
DATA metadata(121) / var_metadata(                                            &
! String identifier
    'land_index',                                                             &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Index (gridbox number) of land points",                                  &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lice_index
!-----------------------------------------------------------------------------
DATA metadata(122) / var_metadata(                                            &
! String identifier
    'lice_index',                                                             &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Index (gridbox number) of land ice points",                              &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lit_c_mean
!-----------------------------------------------------------------------------
DATA metadata(123) / var_metadata(                                            &
! String identifier
    'lit_c_mean',                                                             &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean carbon litter",                                             &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lw_up
!-----------------------------------------------------------------------------
DATA metadata(124) / var_metadata(                                            &
! String identifier
    'lw_up',                                                                  &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox surface upward LW radiation of land points",                     &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for npp_gb
!-----------------------------------------------------------------------------
DATA metadata(125) / var_metadata(                                            &
! String identifier
    'npp_gb',                                                                 &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox net primary productivity prior to N limitation",                 &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for qbase
!-----------------------------------------------------------------------------
DATA metadata(126) / var_metadata(                                            &
! String identifier
    'qbase',                                                                  &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Baseflow (lateral subsurface runoff)",                                   &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for qbase_zw
!-----------------------------------------------------------------------------
DATA metadata(127) / var_metadata(                                            &
! String identifier
    'qbase_zw',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Baseflow from deep LSH/TOPMODEL layer",                                  &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for resp_p_gb
!-----------------------------------------------------------------------------
DATA metadata(128) / var_metadata(                                            &
! String identifier
    'resp_p_gb',                                                              &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox plant respiration",                                              &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for resp_s_gb
!-----------------------------------------------------------------------------
DATA metadata(129) / var_metadata(                                            &
! String identifier
    'resp_s_gb',                                                              &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox soil respiration (total)",                                       &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for resp_s_dr_out
!-----------------------------------------------------------------------------
DATA metadata(130) / var_metadata(                                            &
! String identifier
    'resp_s_dr_out',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean soil respiration for driving TRIFFID",                      &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_root
!-----------------------------------------------------------------------------
DATA metadata(131) / var_metadata(                                            &
! String identifier
    'n_root',                                                                 &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT Root N content scaled by LAI_BAL in sf_stom",                        &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_stem
!-----------------------------------------------------------------------------
DATA metadata(132) / var_metadata(                                            &
! String identifier
    'n_stem',                                                                 &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT Stem N content scaled by LAI_BAL in sf_stom",                        &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_leaf
!-----------------------------------------------------------------------------
DATA metadata(133) / var_metadata(                                            &
! String identifier
    'n_leaf',                                                                 &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT Leaf N content scaled by LAI in sf_stom",                            &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for pc_s
!-----------------------------------------------------------------------------
DATA metadata(134) / var_metadata(                                            &
! String identifier
    'pc_s',                                                                   &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT net carbon available for spreading",                                 &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for runoff
!-----------------------------------------------------------------------------
DATA metadata(135) / var_metadata(                                            &
! String identifier
    'runoff',                                                                 &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox runoff rate",                                                    &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sat_excess_roff
!-----------------------------------------------------------------------------
DATA metadata(136) / var_metadata(                                            &
! String identifier
    'sat_excess_roff',                                                        &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Saturation excess surface ('Dunne') runoff",                             &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for smc_avail_top
!-----------------------------------------------------------------------------
DATA metadata(137) / var_metadata(                                            &
! String identifier
    'smc_avail_top',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name is dynamic - set in get_var_attrs
    "",                                                                       &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for smc_avail_tot
!-----------------------------------------------------------------------------
DATA metadata(138) / var_metadata(                                            &
! String identifier
    'smc_avail_tot',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox available moisture in soil column",                              &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for smc_tot
!-----------------------------------------------------------------------------
DATA metadata(139) / var_metadata(                                            &
! String identifier
    'smc_tot',                                                                &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox total soil moisture in column",                                  &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snomlt_sub_htf
!-----------------------------------------------------------------------------
DATA metadata(140) / var_metadata(                                            &
! String identifier
    'snomlt_sub_htf',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox sub-canopy snowmelt heat flux",                                  &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snow_can_gb
!-----------------------------------------------------------------------------
DATA metadata(141) / var_metadata(                                            &
! String identifier
    'snow_can_gb',                                                            &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox snow on canopy",                                                 &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snow_depth_gb
!-----------------------------------------------------------------------------
DATA metadata(142) / var_metadata(                                            &
! String identifier
    'snow_depth_gb',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Snow depth (on ground)",                                                 &
! Units
    "m"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snow_frac
!-----------------------------------------------------------------------------
DATA metadata(143) / var_metadata(                                            &
! String identifier
    'snow_frac',                                                              &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox snow-covered fraction of land points",                           &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snow_grnd_gb
!-----------------------------------------------------------------------------
DATA metadata(144) / var_metadata(                                            &
! String identifier
    'snow_grnd_gb',                                                           &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox snow below canopy (snow_grnd)",                                  &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snow_ice_gb
!-----------------------------------------------------------------------------
DATA metadata(145) / var_metadata(                                            &
! String identifier
    'snow_ice_gb',                                                            &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox ice content of snow on ground",                                  &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snow_liq_gb
!-----------------------------------------------------------------------------
DATA metadata(146) / var_metadata(                                            &
! String identifier
    'snow_liq_gb',                                                            &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox liquid content of snow on ground",                               &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snow_melt_gb
!-----------------------------------------------------------------------------
DATA metadata(147) / var_metadata(                                            &
! String identifier
    'snow_melt_gb',                                                           &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox rate of snowmelt",                                               &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for soil_index
!-----------------------------------------------------------------------------
DATA metadata(148) / var_metadata(                                            &
! String identifier
    'soil_index',                                                             &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Index (gridbox number) of soil points",                                  &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sub_surf_roff
!-----------------------------------------------------------------------------
DATA metadata(149) / var_metadata(                                            &
! String identifier
    'sub_surf_roff',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox sub-surface runoff",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for surf_roff
!-----------------------------------------------------------------------------
DATA metadata(150) / var_metadata(                                            &
! String identifier
    'surf_roff',                                                              &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox surface runoff",                                                 &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for swet_liq_tot
!-----------------------------------------------------------------------------
DATA metadata(151) / var_metadata(                                            &
! String identifier
    'swet_liq_tot',                                                           &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox unfrozen soil moisture as fraction of saturation",               &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for swet_tot
!-----------------------------------------------------------------------------
DATA metadata(152) / var_metadata(                                            &
! String identifier
    'swet_tot',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox soil moisture as fraction of saturation",                        &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for tfall
!-----------------------------------------------------------------------------
DATA metadata(153) / var_metadata(                                            &
! String identifier
    'tfall',                                                                  &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox throughfall",                                                    &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for trad
!-----------------------------------------------------------------------------
DATA metadata(154) / var_metadata(                                            &
! String identifier
    'trad',                                                                   &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox effective radiative temperature (assuming emissivity=1)",        &
! Units
    "K"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for c_veg
!-----------------------------------------------------------------------------
DATA metadata(155) / var_metadata(                                            &
! String identifier
    'c_veg',                                                                  &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT total carbon content of the vegetation at the end of model timestep.", &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fapar
!-----------------------------------------------------------------------------
DATA metadata(156) / var_metadata(                                            &
! String identifier
    'fapar',                                                                  &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT Fraction of Absorbed Photosynthetically Active Radiation",           &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fao_et0
!-----------------------------------------------------------------------------
DATA metadata(157) / var_metadata(                                            &
! String identifier
    'fao_et0',                                                                &
! Variable type
    var_type_surface,                                                         &
! Long name
    "FAO Penman-Monteith evapotranspiration for reference crop",              &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for flux_o3_stom
!-----------------------------------------------------------------------------
DATA metadata(158) / var_metadata(                                            &
! String identifier
    'flux_o3_stom',                                                           &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Flux of O3 to stomata",                                                  &
! Units
    "mol m-2 s-1"                                                             &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fsmc
!-----------------------------------------------------------------------------
DATA metadata(159) / var_metadata(                                            &
! String identifier
    'fsmc',                                                                   &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT soil moisture availability factor (beta)",                           &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for g_leaf
!-----------------------------------------------------------------------------
DATA metadata(160) / var_metadata(                                            &
! String identifier
    'g_leaf',                                                                 &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT leaf turnover rate",                                                 &
! Units
    "per 360days"                                                             &
  ) /
!-----------------------------------------------------------------------------
! Metadata for g_leaf_day
!-----------------------------------------------------------------------------
DATA metadata(161) / var_metadata(                                            &
! String identifier
    'g_leaf_day',                                                             &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT mean leaf turnover rate for input to PHENOL",                        &
! Units
    "per 360days"                                                             &
  ) /
!-----------------------------------------------------------------------------
! Metadata for g_leaf_dr_out
!-----------------------------------------------------------------------------
DATA metadata(162) / var_metadata(                                            &
! String identifier
    'g_leaf_dr_out',                                                          &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT mean leaf turnover rate for driving TRIFFID",                        &
! Units
    "per 360days"                                                             &
  ) /
!-----------------------------------------------------------------------------
! Metadata for g_leaf_phen
!-----------------------------------------------------------------------------
DATA metadata(163) / var_metadata(                                            &
! String identifier
    'g_leaf_phen',                                                            &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT mean leaf turnover rate over phenology period",                      &
! Units
    "per 360days"                                                             &
  ) /
!-----------------------------------------------------------------------------
! Metadata for gpp
!-----------------------------------------------------------------------------
DATA metadata(164) / var_metadata(                                            &
! String identifier
    'gpp',                                                                    &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT gross primary productivity",                                         &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lai_phen
!-----------------------------------------------------------------------------
DATA metadata(165) / var_metadata(                                            &
! String identifier
    'lai_phen',                                                               &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT leaf area index after phenology",                                    &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lit_c
!-----------------------------------------------------------------------------
DATA metadata(166) / var_metadata(                                            &
! String identifier
    'lit_c',                                                                  &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT carbon Litter",                                                      &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for npp_dr_out
!-----------------------------------------------------------------------------
DATA metadata(167) / var_metadata(                                            &
! String identifier
    'npp_dr_out',                                                             &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT mean NPP for driving TRIFFID/N model ",                              &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for npp
!-----------------------------------------------------------------------------
DATA metadata(168) / var_metadata(                                            &
! String identifier
    'npp',                                                                    &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT net primary productivity prior to N limitation",                     &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for o3_exp_fac
!-----------------------------------------------------------------------------
DATA metadata(169) / var_metadata(                                            &
! String identifier
    'o3_exp_fac',                                                             &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Ozone exposure factor",                                                  &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for resp_p
!-----------------------------------------------------------------------------
DATA metadata(170) / var_metadata(                                            &
! String identifier
    'resp_p',                                                                 &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT plant respiration",                                                  &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for resp_w_dr_out
!-----------------------------------------------------------------------------
DATA metadata(171) / var_metadata(                                            &
! String identifier
    'resp_w_dr_out',                                                          &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT mean wood respiration carbon flux for driving TRIFFID",              &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for resp_l
!-----------------------------------------------------------------------------
DATA metadata(172) / var_metadata(                                            &
! String identifier
    'resp_l',                                                                 &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT leaf maintenance respiration carbon flux (in sf_stom)",              &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for resp_r
!-----------------------------------------------------------------------------
DATA metadata(173) / var_metadata(                                            &
! String identifier
    'resp_r',                                                                 &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT root maintenance respiration carbon flux (in sf_stom)",              &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for resp_w
!-----------------------------------------------------------------------------
DATA metadata(174) / var_metadata(                                            &
! String identifier
    'resp_w',                                                                 &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT wood maintenance respiration carbon flux (in sf_stom)",              &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for resp_s
!-----------------------------------------------------------------------------
DATA metadata(175) / var_metadata(                                            &
! String identifier
    'resp_s',                                                                 &
! Variable type
    var_type_scpool,                                                          &
! Long name is dynamic - set in get_var_attrs
    "",                                                                       &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cosz
!-----------------------------------------------------------------------------
DATA metadata(176) / var_metadata(                                            &
! String identifier
    'cosz',                                                                   &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Cosine of the zenith angle",                                             &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for diff_frac
!-----------------------------------------------------------------------------
DATA metadata(177) / var_metadata(                                            &
! String identifier
    'diff_frac',                                                              &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox fraction of radiation that is diffuse",                          &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for ecan_gb
!-----------------------------------------------------------------------------
DATA metadata(178) / var_metadata(                                            &
! String identifier
    'ecan_gb',                                                                &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean evaporation from canopy/surface store",                     &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for ei_gb
!-----------------------------------------------------------------------------
DATA metadata(179) / var_metadata(                                            &
! String identifier
    'ei_gb',                                                                  &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox sublimation from lying snow or sea-ice",                         &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for et_stom
!-----------------------------------------------------------------------------
DATA metadata(180) / var_metadata(                                            &
! String identifier
    'et_stom',                                                                &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile stomatal transpiration",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for et_stom_gb
!-----------------------------------------------------------------------------
DATA metadata(181) / var_metadata(                                            &
! String identifier
    'et_stom_gb',                                                             &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox stomatal transpiration",                                         &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for esoil_gb
!-----------------------------------------------------------------------------
DATA metadata(182) / var_metadata(                                            &
! String identifier
    'esoil_gb',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox surface evapotranspiration from soil moisture store",            &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fqw_gb
!-----------------------------------------------------------------------------
DATA metadata(183) / var_metadata(                                            &
! String identifier
    'fqw_gb',                                                                 &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox moisture flux from surface",                                     &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for ftl_gb
!-----------------------------------------------------------------------------
DATA metadata(184) / var_metadata(                                            &
! String identifier
    'ftl_gb',                                                                 &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox surface sensible heat flux",                                     &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for land_albedo_1
!-----------------------------------------------------------------------------
DATA metadata(185) / var_metadata(                                            &
! String identifier
    'land_albedo_1',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox band 1 albedo (direct beam visible)",                            &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for land_albedo_2
!-----------------------------------------------------------------------------
DATA metadata(186) / var_metadata(                                            &
! String identifier
    'land_albedo_2',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox band 2 albedo (diffuse visible)",                                &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for land_albedo_3
!-----------------------------------------------------------------------------
DATA metadata(187) / var_metadata(                                            &
! String identifier
    'land_albedo_3',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox band 3 albedo (direct beam NIR)",                                &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for land_albedo_4
!-----------------------------------------------------------------------------
DATA metadata(188) / var_metadata(                                            &
! String identifier
    'land_albedo_4',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox band 4 albedo (diffuse NIR)",                                    &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for latent_heat
!-----------------------------------------------------------------------------
DATA metadata(189) / var_metadata(                                            &
! String identifier
    'latent_heat',                                                            &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox surface latent heat flux",                                       &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for q1p5m_gb
!-----------------------------------------------------------------------------
DATA metadata(190) / var_metadata(                                            &
! String identifier
    'q1p5m_gb',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox specific humidity at 1.5m height",                               &
! Units
    "kg kg-1"                                                                 &
  ) /
!-----------------------------------------------------------------------------
! Metadata for qw1
!-----------------------------------------------------------------------------
DATA metadata(191) / var_metadata(                                            &
! String identifier
    'qw1',                                                                    &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox specific humidity (total water content)",                        &
! Units
    "kg kg-1"                                                                 &
  ) /
!-----------------------------------------------------------------------------
! Metadata for rainfall
!-----------------------------------------------------------------------------
DATA metadata(192) / var_metadata(                                            &
! String identifier
    'rainfall',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox rainfall rate",                                                  &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snomlt_surf_htf
!-----------------------------------------------------------------------------
DATA metadata(193) / var_metadata(                                            &
! String identifier
    'snomlt_surf_htf',                                                        &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox heat flux used for surface melting of snow",                     &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snowfall
!-----------------------------------------------------------------------------
DATA metadata(194) / var_metadata(                                            &
! String identifier
    'snowfall',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox snowfall rate",                                                  &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snow_mass_gb
!-----------------------------------------------------------------------------
DATA metadata(195) / var_metadata(                                            &
! String identifier
    'snow_mass_gb',                                                           &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox snowmass (on land this is lying_snow=snow_tile+snow_grnd)",      &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for surf_ht_flux_gb
!-----------------------------------------------------------------------------
DATA metadata(196) / var_metadata(                                            &
! String identifier
    'surf_ht_flux_gb',                                                        &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox net downward heat flux at surface over land and sea-ice fraction", &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for t1p5m_gb
!-----------------------------------------------------------------------------
DATA metadata(197) / var_metadata(                                            &
! String identifier
    't1p5m_gb',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox temperature at 1.5m height",                                     &
! Units
    "K"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for taux1
!-----------------------------------------------------------------------------
DATA metadata(198) / var_metadata(                                            &
! String identifier
    'taux1',                                                                  &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox westerly component of surface wind stress",                      &
! Units
    "N m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for tauy1
!-----------------------------------------------------------------------------
DATA metadata(199) / var_metadata(                                            &
! String identifier
    'tauy1',                                                                  &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox southerly component of surface wind stress",                     &
! Units
    "N m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for tl1
!-----------------------------------------------------------------------------
DATA metadata(200) / var_metadata(                                            &
! String identifier
    'tl1',                                                                    &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox ice/liquid water temperature",                                   &
! Units
    "K"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for tstar_gb
!-----------------------------------------------------------------------------
DATA metadata(201) / var_metadata(                                            &
! String identifier
    'tstar_gb',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox surface temperature",                                            &
! Units
    "K"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for u1
!-----------------------------------------------------------------------------
DATA metadata(202) / var_metadata(                                            &
! String identifier
    'u1',                                                                     &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox westerly wind component",                                        &
! Units
    "m s-1"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for u10m
!-----------------------------------------------------------------------------
DATA metadata(203) / var_metadata(                                            &
! String identifier
    'u10m',                                                                   &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox westerly wind component at 10 m height",                         &
! Units
    "m s-1"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for v1
!-----------------------------------------------------------------------------
DATA metadata(204) / var_metadata(                                            &
! String identifier
    'v1',                                                                     &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox southerly wind component",                                       &
! Units
    "m s-1"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for v10m
!-----------------------------------------------------------------------------
DATA metadata(205) / var_metadata(                                            &
! String identifier
    'v10m',                                                                   &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox southerly wind component at 10 m height",                        &
! Units
    "m s-1"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for ext
!-----------------------------------------------------------------------------
DATA metadata(206) / var_metadata(                                            &
! String identifier
    'ext',                                                                    &
! Variable type
    var_type_soil,                                                            &
! Long name
    "Gridbox extraction of water from each soil layer",                       &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for smcl
!-----------------------------------------------------------------------------
DATA metadata(207) / var_metadata(                                            &
! String identifier
    'smcl',                                                                   &
! Variable type
    var_type_soil,                                                            &
! Long name
    "Gridbox moisture content of each soil layer",                            &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for soil_wet
!-----------------------------------------------------------------------------
DATA metadata(208) / var_metadata(                                            &
! String identifier
    'soil_wet',                                                               &
! Variable type
    var_type_soil,                                                            &
! Long name
    "Gridbox total moisture content of each soil layer, as fraction of saturation", &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sthf
!-----------------------------------------------------------------------------
DATA metadata(209) / var_metadata(                                            &
! String identifier
    'sthf',                                                                   &
! Variable type
    var_type_soil,                                                            &
! Long name
    "Gridbox frozen moisture content of each soil layer as a fraction of saturation", &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sthu
!-----------------------------------------------------------------------------
DATA metadata(210) / var_metadata(                                            &
! String identifier
    'sthu',                                                                   &
! Variable type
    var_type_soil,                                                            &
! Long name
    "Gridbox unfrozen moisture content of each soil layer as a fraction of saturation", &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for alb_tile_1
!-----------------------------------------------------------------------------
DATA metadata(211) / var_metadata(                                            &
! String identifier
    'alb_tile_1',                                                             &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile land albedo, waveband 1",                                           &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for alb_tile_2
!-----------------------------------------------------------------------------
DATA metadata(212) / var_metadata(                                            &
! String identifier
    'alb_tile_2',                                                             &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile land albedo, waveband 2",                                           &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for alb_tile_3
!-----------------------------------------------------------------------------
DATA metadata(213) / var_metadata(                                            &
! String identifier
    'alb_tile_3',                                                             &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile land albedo, waveband 3",                                           &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for alb_tile_4
!-----------------------------------------------------------------------------
DATA metadata(214) / var_metadata(                                            &
! String identifier
    'alb_tile_4',                                                             &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile land albedo, waveband 4",                                           &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for anthrop_heat
!-----------------------------------------------------------------------------
DATA metadata(215) / var_metadata(                                            &
! String identifier
    'anthrop_heat',                                                           &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Anthropogenic heat flux for each tile",                                  &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for catch
!-----------------------------------------------------------------------------
DATA metadata(216) / var_metadata(                                            &
! String identifier
    'catch',                                                                  &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile surface/canopy water capacity of snow-free land tiles",             &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for ecan
!-----------------------------------------------------------------------------
DATA metadata(217) / var_metadata(                                            &
! String identifier
    'ecan',                                                                   &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile evaporation from canopy/surface store for snow-free land tiles",    &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for ei
!-----------------------------------------------------------------------------
DATA metadata(218) / var_metadata(                                            &
! String identifier
    'ei',                                                                     &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile sublimation from lying snow for land tiles",                        &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for emis
!-----------------------------------------------------------------------------
DATA metadata(219) / var_metadata(                                            &
! String identifier
    'emis',                                                                   &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile emissivity",                                                        &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for esoil
!-----------------------------------------------------------------------------
DATA metadata(220) / var_metadata(                                            &
! String identifier
    'esoil',                                                                  &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile surface evapotranspiration from soil moisture store for snow-free land tiles", &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fqw
!-----------------------------------------------------------------------------
DATA metadata(221) / var_metadata(                                            &
! String identifier
    'fqw',                                                                    &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile surface moisture flux for land tiles",                              &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for ftl
!-----------------------------------------------------------------------------
DATA metadata(222) / var_metadata(                                            &
! String identifier
    'ftl',                                                                    &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile surface sensible heat flux for land tiles",                         &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for gc
!-----------------------------------------------------------------------------
DATA metadata(223) / var_metadata(                                            &
! String identifier
    'gc',                                                                     &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile surface conductance to evaporation for land tiles",                 &
! Units
    "m s-1"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for le
!-----------------------------------------------------------------------------
DATA metadata(224) / var_metadata(                                            &
! String identifier
    'le',                                                                     &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile surface latent heat flux for land tiles",                           &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for q1p5m
!-----------------------------------------------------------------------------
DATA metadata(225) / var_metadata(                                            &
! String identifier
    'q1p5m',                                                                  &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile specific humidity at 1.5m over land tiles",                         &
! Units
    "kg kg-1"                                                                 &
  ) /
!-----------------------------------------------------------------------------
! Metadata for rad_net_tile
!-----------------------------------------------------------------------------
DATA metadata(226) / var_metadata(                                            &
! String identifier
    'rad_net_tile',                                                           &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile surface net radiation",                                             &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sw_surft
!-----------------------------------------------------------------------------
DATA metadata(227) / var_metadata(                                            &
! String identifier
    'sw_surft',                                                               &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile surface shortwave net radiation",                                   &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lw_down_surft
!-----------------------------------------------------------------------------
DATA metadata(228) / var_metadata(                                            &
! String identifier
    'lw_down_surft',                                                          &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile surface longwave downwelling radiation",                            &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lw_up_surft
!-----------------------------------------------------------------------------
DATA metadata(229) / var_metadata(                                            &
! String identifier
    'lw_up_surft',                                                            &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile surface longwave upwelling radiation",                              &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snow_can_melt
!-----------------------------------------------------------------------------
DATA metadata(230) / var_metadata(                                            &
! String identifier
    'snow_can_melt',                                                          &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile melt of canopy snow",                                               &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snow_can
!-----------------------------------------------------------------------------
DATA metadata(231) / var_metadata(                                            &
! String identifier
    'snow_can',                                                               &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile canopy snow",                                                       &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snow_grnd_rho
!-----------------------------------------------------------------------------
DATA metadata(232) / var_metadata(                                            &
! String identifier
    'snow_grnd_rho',                                                          &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile density of snow on ground",                                         &
! Units
    "kg m-3"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snow_ground
!-----------------------------------------------------------------------------
DATA metadata(233) / var_metadata(                                            &
! String identifier
    'snow_ground',                                                            &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile snow on ground (snow_tile or snow_grnd)",                           &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snow_ice_tile
!-----------------------------------------------------------------------------
DATA metadata(234) / var_metadata(                                            &
! String identifier
    'snow_ice_tile',                                                          &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile total ice content in snow on ground",                               &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snow_liq_tile
!-----------------------------------------------------------------------------
DATA metadata(235) / var_metadata(                                            &
! String identifier
    'snow_liq_tile',                                                          &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile total liquid content in snow on ground",                            &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snow_mass
!-----------------------------------------------------------------------------
DATA metadata(236) / var_metadata(                                            &
! String identifier
    'snow_mass',                                                              &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile lying snow (snow_tile+snow_grnd)",                                  &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snow_melt
!-----------------------------------------------------------------------------
DATA metadata(237) / var_metadata(                                            &
! String identifier
    'snow_melt',                                                              &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile snow melt (melt_tile)",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for surf_ht_flux
!-----------------------------------------------------------------------------
DATA metadata(238) / var_metadata(                                            &
! String identifier
    'surf_ht_flux',                                                           &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Downward heat flux for each tile",                                       &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snow_soil_htf
!-----------------------------------------------------------------------------
DATA metadata(239) / var_metadata(                                            &
! String identifier
    'snow_soil_htf',                                                          &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Downward heat flux after snowpack to subsurface",                        &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snice_smb_surft
!-----------------------------------------------------------------------------
DATA metadata(240) / var_metadata(                                            &
! String identifier
    'snice_smb_surft',                                                        &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Rate of change of snowpack mass on tiles",                               &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snice_m_surft
!-----------------------------------------------------------------------------
DATA metadata(241) / var_metadata(                                            &
! String identifier
    'snice_m_surft',                                                          &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Total internal melt rate of snowpack on tiles",                          &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snice_freez_surft
!-----------------------------------------------------------------------------
DATA metadata(242) / var_metadata(                                            &
! String identifier
    'snice_freez_surft',                                                      &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Total internal refreezing rate in snowpack on tiles",                    &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snice_runoff_surft
!-----------------------------------------------------------------------------
DATA metadata(243) / var_metadata(                                            &
! String identifier
    'snice_runoff_surft',                                                     &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Net rate of liquid leaving snowpack on tiles",                           &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snice_sicerate_surft
!-----------------------------------------------------------------------------
DATA metadata(244) / var_metadata(                                            &
! String identifier
    'snice_sicerate_surft',                                                   &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Rate of change of solid mass in snowpack on tiles",                      &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for snice_sliqrate_surft
!-----------------------------------------------------------------------------
DATA metadata(245) / var_metadata(                                            &
! String identifier
    'snice_sliqrate_surft',                                                   &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Rate of change of liquid mass in snowpack on tiles",                     &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for surf_ht_store
!-----------------------------------------------------------------------------
DATA metadata(246) / var_metadata(                                            &
! String identifier
    'surf_ht_store',                                                          &
! Variable type
    var_type_surft,                                                           &
! Long name
    "C*(dT/dt) for each tile",                                                &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for t1p5m
!-----------------------------------------------------------------------------
DATA metadata(247) / var_metadata(                                            &
! String identifier
    't1p5m',                                                                  &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile temperature at 1.5m over land tiles",                               &
! Units
    "K"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for tstar
!-----------------------------------------------------------------------------
DATA metadata(248) / var_metadata(                                            &
! String identifier
    'tstar',                                                                  &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile surface temperature",                                               &
! Units
    "K"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for tsurf_elev_surft
!-----------------------------------------------------------------------------
DATA metadata(249) / var_metadata(                                            &
! String identifier
    'tsurf_elev_surft',                                                       &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Temperature of elevated subsurface tiles",                               &
! Units
    "K"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for z0
!-----------------------------------------------------------------------------
DATA metadata(250) / var_metadata(                                            &
! String identifier
    'z0',                                                                     &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile surface roughness",                                                 &
! Units
    "m"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for z0h
!-----------------------------------------------------------------------------
DATA metadata(251) / var_metadata(                                            &
! String identifier
    'z0h',                                                                    &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile surface thermal roughness",                                         &
! Units
    "m"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for tile_index
!-----------------------------------------------------------------------------
DATA metadata(252) / var_metadata(                                            &
! String identifier
    'tile_index',                                                             &
! Variable type
    var_type_type,                                                            &
! Long name
    "Index (gridbox number) of land points with each surface type",           &
! Units
    "m"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for co2_mmr
!-----------------------------------------------------------------------------
DATA metadata(253) / var_metadata(                                            &
! String identifier
    'co2_mmr',                                                                &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Concentration of atmospheric CO2, expressed as a mass mixing ratio.",    &
! Units
    ""                                                                        &
  ) /
!-----------------------------------------------------------------------------
! Metadata for ozone
!-----------------------------------------------------------------------------
DATA metadata(254) / var_metadata(                                            &
! String identifier
    'ozone',                                                                  &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for isoprene_gb
!-----------------------------------------------------------------------------
DATA metadata(255) / var_metadata(                                            &
! String identifier
    'isoprene_gb',                                                            &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean isoprene emission flux",                                    &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for isoprene
!-----------------------------------------------------------------------------
DATA metadata(256) / var_metadata(                                            &
! String identifier
    'isoprene',                                                               &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT isoprene emission flux",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for terpene_gb
!-----------------------------------------------------------------------------
DATA metadata(257) / var_metadata(                                            &
! String identifier
    'terpene_gb',                                                             &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean (mono-)terpene emission flux",                              &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for terpene
!-----------------------------------------------------------------------------
DATA metadata(258) / var_metadata(                                            &
! String identifier
    'terpene',                                                                &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT (mono-)terpene emission flux",                                       &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for methanol_gb
!-----------------------------------------------------------------------------
DATA metadata(259) / var_metadata(                                            &
! String identifier
    'methanol_gb',                                                            &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean methanol emission flux",                                    &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for methanol
!-----------------------------------------------------------------------------
DATA metadata(260) / var_metadata(                                            &
! String identifier
    'methanol',                                                               &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT methanol emission flux",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for acetone_gb
!-----------------------------------------------------------------------------
DATA metadata(261) / var_metadata(                                            &
! String identifier
    'acetone_gb',                                                             &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean acetone emission flux",                                     &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for acetone
!-----------------------------------------------------------------------------
DATA metadata(262) / var_metadata(                                            &
! String identifier
    'acetone',                                                                &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT acetone emission flux",                                              &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for croprootc
!-----------------------------------------------------------------------------
DATA metadata(263) / var_metadata(                                            &
! String identifier
    'croprootc',                                                              &
! Variable type
    var_type_cpft,                                                            &
! Long name
    "Root carbon of crop at the end of model timestep",                       &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cropharvc
!-----------------------------------------------------------------------------
DATA metadata(264) / var_metadata(                                            &
! String identifier
    'cropharvc',                                                              &
! Variable type
    var_type_cpft,                                                            &
! Long name
    "Carbon in harvested parts of crop at the end of model timestep",         &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cropreservec
!-----------------------------------------------------------------------------
DATA metadata(265) / var_metadata(                                            &
! String identifier
    'cropreservec',                                                           &
! Variable type
    var_type_cpft,                                                            &
! Long name
    "Carbon in crop stem reserve pool at the end of model timestep",          &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cropdvi
!-----------------------------------------------------------------------------
DATA metadata(266) / var_metadata(                                            &
! String identifier
    'cropdvi',                                                                &
! Variable type
    var_type_cpft,                                                            &
! Long name
    "Developmental index of crops at the end of model timestep",              &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cropyield
!-----------------------------------------------------------------------------
DATA metadata(267) / var_metadata(                                            &
! String identifier
    'cropyield',                                                              &
! Variable type
    var_type_cpft,                                                            &
! Long name
    "Yield carbon at the end of model timestep",                              &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for harvest_trigger
!-----------------------------------------------------------------------------
DATA metadata(268) / var_metadata(                                            &
! String identifier
    'harvest_trigger',                                                        &
! Variable type
    var_type_cpft,                                                            &
! Long name
    "Trigger condition for crop harvest",                                     &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for harvest_counter
!-----------------------------------------------------------------------------
DATA metadata(269) / var_metadata(                                            &
! String identifier
    'harvest_counter',                                                        &
! Variable type
    var_type_cpft,                                                            &
! Long name
    "Counter for crop harvest",                                               &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cropstemc
!-----------------------------------------------------------------------------
DATA metadata(270) / var_metadata(                                            &
! String identifier
    'cropstemc',                                                              &
! Variable type
    var_type_cpft,                                                            &
! Long name
    "Carbon in stem parts of crop at the end of model timestep",              &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cropleafc
!-----------------------------------------------------------------------------
DATA metadata(271) / var_metadata(                                            &
! String identifier
    'cropleafc',                                                              &
! Variable type
    var_type_cpft,                                                            &
! Long name
    "Carbon in leaf parts of crop at the end of model timestep",              &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for croplai
!-----------------------------------------------------------------------------
DATA metadata(272) / var_metadata(                                            &
! String identifier
    'croplai',                                                                &
! Variable type
    var_type_cpft,                                                            &
! Long name
    "Leaf area index of crop at the end of model timestep",                   &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cropcanht
!-----------------------------------------------------------------------------
DATA metadata(273) / var_metadata(                                            &
! String identifier
    'cropcanht',                                                              &
! Variable type
    var_type_cpft,                                                            &
! Long name
    "Canopy height of crop at the end of model timestep",                     &
! Units
    "m"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cropsowdate
!-----------------------------------------------------------------------------
DATA metadata(274) / var_metadata(                                            &
! String identifier
    'cropsowdate',                                                            &
! Variable type
    var_type_cpft,                                                            &
! Long name
    "Day of year when crop is sown",                                          &
! Units
    "day of year"                                                             &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cropttveg
!-----------------------------------------------------------------------------
DATA metadata(275) / var_metadata(                                            &
! String identifier
    'cropttveg',                                                              &
! Variable type
    var_type_cpft,                                                            &
! Long name
    "Thermal time between emergence and flowering",                           &
! Units
    "degree days"                                                             &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cropttrep
!-----------------------------------------------------------------------------
DATA metadata(276) / var_metadata(                                            &
! String identifier
    'cropttrep',                                                              &
! Variable type
    var_type_cpft,                                                            &
! Long name
    "Thermal time between flowering and maturity/harvest",                    &
! Units
    "degree days"                                                             &
  ) /
!-----------------------------------------------------------------------------
! Neutral winds
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Metadata for u10m_n
!-----------------------------------------------------------------------------
DATA metadata(277) / var_metadata(                                            &
! String identifier
    'u10m_n',                                                                 &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox equivalent neutral westerly wind component at 10 m height",      &
! Units
    "m s-1"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for v10m_n
!-----------------------------------------------------------------------------
DATA metadata(278) / var_metadata(                                            &
! String identifier
    'v10m_n',                                                                 &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox equivalent neutral southerly wind component at 10 m height",     &
! Units
    "m s-1"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for mu10m_n
!-----------------------------------------------------------------------------
DATA metadata(279) / var_metadata(                                            &
! String identifier
    'mu10m_n',                                                                &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox westerly pseudostress referenced to 10 m height",                &
! Units
    "m2 s-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for mv10m_n
!-----------------------------------------------------------------------------
DATA metadata(280) / var_metadata(                                            &
! String identifier
    'mv10m_n',                                                                &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox southerly pseudostress referenced to 10 m height",               &
! Units
    "m2 s-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sthu_irr
!-----------------------------------------------------------------------------
DATA metadata(281) / var_metadata(                                            &
! String identifier
    'sthu_irr',                                                               &
! Variable type
    var_type_soil,                                                            &
! Long name
    "Gridbox wetness of each soil layer over irrigation",                     &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for frac_irrig
!-----------------------------------------------------------------------------
DATA metadata(282) / var_metadata(                                            &
! String identifier
    'frac_irrig',                                                             &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Irrigated fraction",                                                     &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for irrDaysDiag
!-----------------------------------------------------------------------------
DATA metadata(283) / var_metadata(                                            &
! String identifier
    'irrDaysDiag',                                                            &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Number of days on which irrigation is applied",                          &
! Units
    "days"                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for irrig_water
!-----------------------------------------------------------------------------
DATA metadata(284) / var_metadata(                                            &
! String identifier
    'irrig_water',                                                            &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Water applied as irrigation",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for frac_agr
!-----------------------------------------------------------------------------
DATA metadata(285) / var_metadata(                                            &
! String identifier
    'frac_agr',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox agricultural/crop fraction",                                     &
! Units
    "-"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for frac_agr_prev
!-----------------------------------------------------------------------------
DATA metadata(286) / var_metadata(                                            &
! String identifier
    'frac_agr_prev',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox agricultural/crop fraction from previous land use update",       &
! Units
    "-"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for frac_past
!-----------------------------------------------------------------------------
DATA metadata(287) / var_metadata(                                            &
! String identifier
    'frac_past',                                                              &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox pasture fraction",                                               &
! Units
    "-"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for frac_past_prev
!-----------------------------------------------------------------------------
DATA metadata(288) / var_metadata(                                            &
! String identifier
    'frac_past_prev',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox pasture fraction from previous land use update",                 &
! Units
    "-"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for Wood_prod_fast
!-----------------------------------------------------------------------------
DATA metadata(289) / var_metadata(                                            &
! String identifier
    'wood_prod_fast',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Wood Products Fast C Pool",                                              &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for Wood_prod_med
!-----------------------------------------------------------------------------
DATA metadata(290) / var_metadata(                                            &
! String identifier
    'wood_prod_med',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Wood Products Medium C Pool",                                            &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for Wood_prod_slow
!-----------------------------------------------------------------------------
DATA metadata(291) / var_metadata(                                            &
! String identifier
    'wood_prod_slow',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Wood Products Slow C Pool",                                              &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for WP_fast_in
!-----------------------------------------------------------------------------
DATA metadata(292) / var_metadata(                                            &
! String identifier
    'WP_fast_in',                                                             &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Wood Products Fast C Pool Input",                                        &
! Units
    "kg m-2 per 360 days"                                                     &
  ) /
!-----------------------------------------------------------------------------
! Metadata for WP_med_in
!-----------------------------------------------------------------------------
DATA metadata(293) / var_metadata(                                            &
! String identifier
    'WP_med_in',                                                              &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Wood Products Medium C Pool Input",                                      &
! Units
    "kg m-2 per 360 days"                                                     &
  ) /

!-----------------------------------------------------------------------------
! Metadata for WP_slow_in
!-----------------------------------------------------------------------------
DATA metadata(294) / var_metadata(                                            &
! String identifier
    'WP_slow_in',                                                             &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Wood Products Slow C Pool Input",                                        &
! Units
    "kg m-2 per 360 days"                                                     &
  ) /

!-----------------------------------------------------------------------------
! Metadata for WP_fast_out
!-----------------------------------------------------------------------------
DATA metadata(295) / var_metadata(                                            &
! String identifier
    'WP_fast_out',                                                            &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Wood Products Fast C Pool Output",                                       &
! Units
    "kg m-2 per 360 days"                                                     &
  ) /
!-----------------------------------------------------------------------------
! Metadata for WP_med_out
!-----------------------------------------------------------------------------
DATA metadata(296) / var_metadata(                                            &
! String identifier
    'WP_med_out',                                                             &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Wood Products Medium C Pool Output",                                     &
! Units
    "kg m-2 per 360 days"                                                     &
  ) /
!-----------------------------------------------------------------------------
! Metadata for WP_slow_out
!-----------------------------------------------------------------------------
DATA metadata(297) / var_metadata(                                            &
! String identifier
    'WP_slow_out',                                                            &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Wood Products Slow C Pool Output",                                       &
! Units
    "kg m-2 per 360 days"                                                     &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cnsrv_carbon_veg2
!-----------------------------------------------------------------------------
DATA metadata(298) / var_metadata(                                            &
! String identifier
    'cnsrv_carbon_veg2',                                                      &
! Variable type
    var_type_surface,                                                         &
! Long name
    "error in land carbon conservation in veg2 routine",                      &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cnsrv_carbon_triffid
!-----------------------------------------------------------------------------
DATA metadata(299) / var_metadata(                                            &
! String identifier
    'cnsrv_carbon_triffid',                                                   &
! Variable type
    var_type_surface,                                                         &
! Long name
    "error in land carbon conservation in triffid routine",                   &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cnsrv_veg_triffid
!-----------------------------------------------------------------------------
DATA metadata(300) / var_metadata(                                            &
! String identifier
    'cnsrv_veg_triffid',                                                      &
! Variable type
    var_type_surface,                                                         &
! Long name
    "error in vegetation carbon conservation in triffid routine",             &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cnsrv_soil_triffid
!-----------------------------------------------------------------------------
DATA metadata(301) / var_metadata(                                            &
! String identifier
    'cnsrv_soil_triffid',                                                     &
! Variable type
    var_type_surface,                                                         &
! Long name
    "error in soil carbon conservation in triffid routine",                   &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cnsrv_prod_triffid
!-----------------------------------------------------------------------------
DATA metadata(302) / var_metadata(                                            &
! String identifier
    'cnsrv_prod_triffid',                                                     &
! Variable type
    var_type_surface,                                                         &
! Long name
    "error in wood product carbon conservation in triffid routine",           &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lit_c_orig
!-----------------------------------------------------------------------------
DATA metadata(303) / var_metadata(                                            &
! String identifier
    'lit_c_orig',                                                             &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT carbon Litter including LU and fire",                                &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lit_n_orig
!-----------------------------------------------------------------------------
DATA metadata(304) / var_metadata(                                            &
! String identifier
    'lit_n_orig',                                                             &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT nitrogen Litter including LU and fire",                              &
! Units
    "kg N  m-2 per 360days"                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_inorg
!-----------------------------------------------------------------------------
DATA metadata(305) / var_metadata(                                            &
! String identifier
    'n_inorg',                                                                &
! Variable type
    var_type_sclayer,                                                         &
! Long name
    "Inorganic N Pool",                                                       &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_loss
!-----------------------------------------------------------------------------
DATA metadata(306) / var_metadata(                                            &
! String identifier
    'n_loss',                                                                 &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Nitrogen loss term (fixed fraction of N_INORG)",                         &
! Units
    "kg m-2 360 days"                                                         &
  ) /
!-----------------------------------------------------------------------------
! Metadata for dpm_ratio
!-----------------------------------------------------------------------------
DATA metadata(307) / var_metadata(                                            &
! String identifier
    'dpm_ratio',                                                              &
! Variable type
    var_type_surface,                                                         &
! Long name
    "DPM:RPM ratio of overall litter input (one value per GB)",               &
! Units
    "fraction (limited to between 0 and 1)"                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lit_c_ag
!-----------------------------------------------------------------------------
DATA metadata(308) / var_metadata(                                            &
! String identifier
    'lit_c_ag',                                                               &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT carbon Litter from LU/Agriculture",                                  &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for root_abandon
!-----------------------------------------------------------------------------
DATA metadata(309) / var_metadata(                                            &
! String identifier
    'root_abandon',                                                           &
! Variable type
    var_type_pft,                                                             &
! Long name
    "carbon flux from roots abandoned during landuse change to soil",         &
! Units
    "kg C  m-2 per 360days"                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for root_abandon_n
!-----------------------------------------------------------------------------
DATA metadata(310) / var_metadata(                                            &
! String identifier
    'root_abandon_n',                                                         &
! Variable type
    var_type_pft,                                                             &
! Long name
    "nitrogen flux from roots abandoned during landuse change to soil",       &
! Units
    "kg N  m-2 per 360days"                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for harvest
!-----------------------------------------------------------------------------
DATA metadata(311) / var_metadata(                                            &
! String identifier
    'harvest',                                                                &
! Variable type
    var_type_pft,                                                             &
! Long name
    "flux of carbon to product pools due to harvest",                         &
! Units
    "kg C  m-2 per 360days"                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for harvest_gb
!-----------------------------------------------------------------------------
DATA metadata(312) / var_metadata(                                            &
! String identifier
    'harvest_gb',                                                             &
! Variable type
    var_type_surface,                                                         &
! Long name
    "gridbox mean flux of carbon to product pools due to harvest",            &
! Units
    "kg C  m-2 per 360days"                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for harvest_n
!-----------------------------------------------------------------------------
DATA metadata(313) / var_metadata(                                            &
! String identifier
    'harvest_n',                                                              &
! Variable type
    var_type_pft,                                                             &
! Long name
    "flux of nitrogen to atmosphere due to harvest",                          &
! Units
    "kg N  m-2 per 360days"                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for harvest_n_gb
!-----------------------------------------------------------------------------
DATA metadata(314) / var_metadata(                                            &
! String identifier
    'harvest_n_gb',                                                           &
! Variable type
    var_type_surface,                                                         &
! Long name
    "gridbox mean flux of nitrogen to product pools due to harvest",          &
! Units
    "kg N  m-2 per 360days"                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_fertiliser
!-----------------------------------------------------------------------------
DATA metadata(315) / var_metadata(                                            &
! String identifier
    'n_fertiliser',                                                           &
! Variable type
    var_type_pft,                                                             &
! Long name
    "N addition from fertiliser",                                             &
! Units
    "kg N  m-2 per 360days"                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_fertiliser_gb
!-----------------------------------------------------------------------------
DATA metadata(316) / var_metadata(                                            &
! String identifier
    'n_fertiliser_gb',                                                        &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox N addition from fertiliser",                                     &
! Units
    "kg N  m-2 per 360days"                                                   &
  ) /
 !-----------------------------------------------------------------------------
! Metadata for lit_n_ag
!-----------------------------------------------------------------------------
DATA metadata(317) / var_metadata(                                            &
! String identifier
    'lit_n_ag',                                                               &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT nitrogen loss due to LU/Agriculture",                                &
! Units
    "kg N  m-2 per 360days"                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for river routing latitude (2d)
!-----------------------------------------------------------------------------
DATA metadata(318) / var_metadata(                                            &
! String identifier
    'latitude_2d',                                                            &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Routing gridbox latitude",                                               &
! Units
    "degrees"                                                                 &
  ) /
!-----------------------------------------------------------------------------
! Metadata for river routing longitude (2d)
!-----------------------------------------------------------------------------
DATA metadata(319) / var_metadata(                                            &
! String identifier
    'longitude_2d',                                                           &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Routing gridbox latitude",                                               &
! Units
    "degrees"                                                                 &
  ) /
!-----------------------------------------------------------------------------
! Metadata for river direction
!-----------------------------------------------------------------------------
DATA metadata(320) / var_metadata(                                            &
! String identifier
    'direction',                                                              &
! Variable type
    var_type_surface,                                                         &
! Long name
    "River routing direction flag",                                           &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for river sequence
!-----------------------------------------------------------------------------
DATA metadata(321) / var_metadata(                                            &
! String identifier
    'sequence',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "River routing sequence flag",                                            &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for river drainage area
!-----------------------------------------------------------------------------
DATA metadata(322) / var_metadata(                                            &
! String identifier
    'area',                                                                   &
! Variable type
    var_type_surface,                                                         &
! Long name
    "River routing drainage area",                                            &
! Units
    "number of grid cells"                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for river routing latitude (1d)
!-----------------------------------------------------------------------------
DATA metadata(323) / var_metadata(                                            &
! String identifier
    'rivers_ygrid',                                                           &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Routing gridbox latitude",                                               &
! Units
    "degrees"                                                                 &
  ) /
!-----------------------------------------------------------------------------
! Metadata for river routing longitude (1d)
!-----------------------------------------------------------------------------
DATA metadata(324) / var_metadata(                                            &
! String identifier
    'rivers_xgrid',                                                           &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Routing gridbox longitude",                                              &
! Units
    "degrees"                                                                 &
  ) /
!-----------------------------------------------------------------------------
! Metadata for river storage on river points
!-----------------------------------------------------------------------------
DATA metadata(325) / var_metadata(                                            &
! String identifier
    'rivers_sto_rp',                                                          &
! Variable type
    var_type_rp,                                                              &
! Long name
    "River storage on river points",                                          &
! Units
    "kg"                                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for river flow
!-----------------------------------------------------------------------------
DATA metadata(326) / var_metadata(                                            &
! String identifier
    'rflow',                                                                  &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Routing gridbox flow",                                                   &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for river routing grid runoff
!-----------------------------------------------------------------------------
DATA metadata(327) / var_metadata(                                            &
! String identifier
    'rrun',                                                                   &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Routing gridbox runoff",                                                 &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for Metstat temp_max_00h_r
!-----------------------------------------------------------------------------
DATA metadata(328) / var_metadata(                                            &
! String identifier
    'temp_max_00h_r',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Maximum temperature to 00hrs local- running total",                      &
! Units
    "Kelvin"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for Metstat temp_ave_00h_r
!-----------------------------------------------------------------------------
DATA metadata(329) / var_metadata(                                            &
! String identifier
    'temp_ave_00h_r',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Average temperature to 00hrs local- running total",                      &
! Units
    "Kelvin"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for Metstat prec_tot_00h_r
!-----------------------------------------------------------------------------
DATA metadata(330) / var_metadata(                                            &
! String identifier
    'prec_tot_00h_r',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Total precipitation to 00hrs local- running total",                      &
! Units
    "kg/ms2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for Metstat prec_tot_12h_r
!-----------------------------------------------------------------------------
DATA metadata(331) / var_metadata(                                            &
! String identifier
    'prec_tot_12h_r',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Total precipitation to 12hrs local- running total",                      &
! Units
    "kg/m2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for Metstat rhum_min_00h_r
!-----------------------------------------------------------------------------
DATA metadata(332) / var_metadata(                                            &
! String identifier
    'rhum_min_00h_r',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Minimum relative humidity to 00hrs local- running total",                &
! Units
    "%"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for Metstat dewp_ave_00h_r
!-----------------------------------------------------------------------------
DATA metadata(333) / var_metadata(                                            &
! String identifier
    'dewp_ave_00h_r',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Average dewpoint to 00 hrs local- running total",                        &
! Units
    "Kelvin"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for Metstat wind_ave_00h_r
!-----------------------------------------------------------------------------
DATA metadata(334) / var_metadata(                                            &
! String identifier
    'wind_ave_00h_r',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Average wind speed to 00hrs local- running total",                       &
! Units
    "m/s"                                                                     &
  ) /
!-----------------------------------------------------------------------------
! Metadata for Metstat temp_max_00h
!-----------------------------------------------------------------------------
DATA metadata(335) / var_metadata(                                            &
! String identifier
    'temp_max_00h',                                                           &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Maximum temperature to 00hrs local",                                     &
! Units
    "Kelvin"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for Metstat temp_ave_00h
!-----------------------------------------------------------------------------
DATA metadata(336) / var_metadata(                                            &
! String identifier
    'temp_ave_00h',                                                           &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Average temperature to 00hrs local",                                     &
! Units
    "Kelvin"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for Metstat temp_pnt_12h
!-----------------------------------------------------------------------------
DATA metadata(337) / var_metadata(                                            &
! String identifier
    'temp_pnt_12h',                                                           &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Timestep temperature at last 12hrs local",                               &
! Units
    "Kelvin"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for Metstat prec_tot_00h
!-----------------------------------------------------------------------------
DATA metadata(338) / var_metadata(                                            &
! String identifier
    'prec_tot_00h',                                                           &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Total precipitation to 00hrs local",                                     &
! Units
    "kg/ms2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for Metstat prec_tot_12h
!-----------------------------------------------------------------------------
DATA metadata(339) / var_metadata(                                            &
! String identifier
    'prec_tot_12h',                                                           &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Total precipitation to 12hrs local",                                     &
! Units
    "kg/m2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for Metstat rhum_min_00h
!-----------------------------------------------------------------------------
DATA metadata(340) / var_metadata(                                            &
! String identifier
    'rhum_min_00h',                                                           &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Minimum relative humidity to 00hrs local",                               &
! Units
    "%"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for Metstat rhum_pnt_12h
!-----------------------------------------------------------------------------
DATA metadata(341) / var_metadata(                                            &
! String identifier
    'rhum_pnt_12h',                                                           &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Timestep relative humidity at last 12hrs local",                         &
! Units
    "%"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for Metstat dewp_ave_00h
!-----------------------------------------------------------------------------
DATA metadata(342) / var_metadata(                                            &
! String identifier
    'dewp_ave_00h',                                                           &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Average dewpoint to 00 hrs local",                                       &
! Units
    "Kelvin"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for Metstat wind_ave_00h
!-----------------------------------------------------------------------------
DATA metadata(343) / var_metadata(                                            &
! String identifier
    'wind_ave_00h',                                                           &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Average wind speed to 00hrs local",                                      &
! Units
    "m/s"                                                                     &
  ) /
!-----------------------------------------------------------------------------
! Metadata for Metstat wind_pnt_12h
!-----------------------------------------------------------------------------
DATA metadata(344) / var_metadata(                                            &
! String identifier
    'wind_pnt_12h',                                                           &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Timestep wind speed at last 12hrs local",                                &
! Units
    "m/s"                                                                     &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_mcarthur
!-----------------------------------------------------------------------------
DATA metadata(345) / var_metadata(                                            &
! String identifier
    'fire_mcarthur',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "McArthur FFDI",                                                          &
! Units
    "None"                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_mcarthur_r_dr
!-----------------------------------------------------------------------------
DATA metadata(346) / var_metadata(                                            &
! String identifier
    'fire_mcarthur_r_dr',                                                     &
! Variable type
    var_type_surface,                                                         &
! Long name
    "McArthur FFDI- Drought R term",                                          &
! Units
    "None"                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_mcarthur_n_dr
!-----------------------------------------------------------------------------
DATA metadata(347) / var_metadata(                                            &
! String identifier
    'fire_mcarthur_n_dr',                                                     &
! Variable type
    var_type_surface,                                                         &
! Long name
    "McArthur FFDI- Drought N term",                                          &
! Units
     "None"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_canadian
!-----------------------------------------------------------------------------
DATA metadata(348) / var_metadata(                                            &
! String identifier
    'fire_canadian',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Canadian FWI",                                                           &
! Units
      "None"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_canadian_ffmc
!-----------------------------------------------------------------------------
DATA metadata(349) / var_metadata(                                            &
! String identifier
    'fire_canadian_ffmc',                                                     &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Canadian FWI- FFMC",                                                     &
! Units
    "None"                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_canadian_ffmc_mois
!-----------------------------------------------------------------------------
DATA metadata(350) / var_metadata(                                            &
! String identifier
    'fire_canadian_ffmc_mois',                                                &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Canadian FWI- FFMC Moisture",                                            &
! Units
    "None"                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_canadian_dmc
!-----------------------------------------------------------------------------
DATA metadata(351) / var_metadata(                                            &
! String identifier
    'fire_canadian_dmc',                                                      &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Canadian FWI- DMC",                                                      &
! Units
    "None"                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_canadian_dc
!-----------------------------------------------------------------------------
DATA metadata(352) / var_metadata(                                            &
! String identifier
    'fire_canadian_dc',                                                       &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Canadian FWI- DC",                                                       &
! Units
    "None"                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_canadian_isi
!-----------------------------------------------------------------------------
DATA metadata(353) / var_metadata(                                            &
! String identifier
    'fire_canadian_isi',                                                      &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Canadian FWI- ISI",                                                      &
! Units
    "None"                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_canadian_bui
!-----------------------------------------------------------------------------
DATA metadata(354) / var_metadata(                                            &
! String identifier
    'fire_canadian_bui',                                                      &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Canadian FWI- BUI",                                                      &
! Units
    "None"                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_nesterov
!-----------------------------------------------------------------------------
DATA metadata(355) / var_metadata(                                            &
! String identifier
    'fire_nesterov',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Nesterov Index",                                                         &
! Units
    "None"                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for leafC
!-----------------------------------------------------------------------------
DATA metadata(356) / var_metadata(                                            &
! String identifier
    'leafC',                                                                  &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Carbon in LEAF biomass on PFTs",                                         &
! Units
    "kg/m2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for rootC
!-----------------------------------------------------------------------------
DATA metadata(357) / var_metadata(                                            &
! String identifier
    'rootC',                                                                  &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Carbon in ROOT biomass on PFTs",                                         &
! Units
    "kg/m2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for stemC
!-----------------------------------------------------------------------------
DATA metadata(358) / var_metadata(                                            &
! String identifier
    'stemC',                                                                  &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Carbon in STEM biomass on PFTs",                                         &
! Units
    "kg/m2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for woodC
!-----------------------------------------------------------------------------
DATA metadata(359) / var_metadata(                                            &
! String identifier
    'woodC',                                                                  &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Carbon in WOODY biomass on PFTs",                                        &
! Units
    "kg/m2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for dleaf
!-----------------------------------------------------------------------------
DATA metadata(360) / var_metadata(                                            &
! String identifier
    'dleaf',                                                                  &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Increment to carbon in LEAF biomass on PFTs",                            &
! Units
    "kg/m2/TRIFFID PERIOD"                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for droot
!-----------------------------------------------------------------------------
DATA metadata(361) / var_metadata(                                            &
! String identifier
    'droot',                                                                  &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Increment to carbon in ROOT biomass on PFTs",                            &
! Units
    "kg/m2/TRIFFID PERIOD"                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for dwood
!-----------------------------------------------------------------------------
DATA metadata(362) / var_metadata(                                            &
! String identifier
    'dwood',                                                                  &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Increment to carbon in WOODY biomass on PFTs",                           &
! Units
    "kg/m2/TRIFFID PERIOD"                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for leaf_litC
!-----------------------------------------------------------------------------
DATA metadata(363) / var_metadata(                                            &
! String identifier
    'leaf_litC',                                                              &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Litter CARBON due to leaf turnover",                                     &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for root_litC
!-----------------------------------------------------------------------------
DATA metadata(364) / var_metadata(                                            &
! String identifier
    'root_litC',                                                              &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Litter CARBON due to root turnover",                                     &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for wood_litC
!-----------------------------------------------------------------------------
DATA metadata(365) / var_metadata(                                            &
! String identifier
    'wood_litC',                                                              &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Litter CARBON due to wood turnover",                                     &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for leaf_litN
!-----------------------------------------------------------------------------
DATA metadata(366) / var_metadata(                                            &
! String identifier
    'leaf_litN',                                                              &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Litter NITROGEN due to leaf turnover",                                   &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for exudates
!-----------------------------------------------------------------------------
DATA metadata(367) / var_metadata(                                            &
! String identifier
    'exudates',                                                               &
! Variable type
    var_type_pft,                                                             &
! Long name
    "EXUDATES on PFTs: excess C not assimilable into plant due lack of N availability", &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for exudates_gb
!-----------------------------------------------------------------------------
DATA metadata(368) / var_metadata(                                            &
! String identifier
    'exudates_gb',                                                            &
! Variable type
    var_type_surface,                                                         &
! Long name
    "EXUDATES (GBM): excess C not assimilable into plant due lack of N availability", &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for npp after N limitation
!-----------------------------------------------------------------------------
DATA metadata(369) / var_metadata(                                            &
! String identifier
    'npp_n',                                                                  &
! Variable type
    var_type_pft,                                                             &
! Long name
    "NPP on PFTs post N-limitation",                                          &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for npp_n_gb
!-----------------------------------------------------------------------------
DATA metadata(370) / var_metadata(                                            &
! String identifier
    'npp_n_gb',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "NPP (GBM) post N-limitation",                                            &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lai_bal
!-----------------------------------------------------------------------------
DATA metadata(371) / var_metadata(                                            &
! String identifier
    'lai_bal',                                                                &
! Variable type
    var_type_pft,                                                             &
! Long name
    "LAI in balanced growth state (seasonal maximum LAI)",                    &
! Units
    "m2/m2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for litterN
!-----------------------------------------------------------------------------
DATA metadata(372) / var_metadata(                                            &
! String identifier
    'litterN',                                                                &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Local Littter N production",                                             &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_uptake_growth
!-----------------------------------------------------------------------------
DATA metadata(373) / var_metadata(                                            &
! String identifier
    'n_uptake_growth',                                                        &
! Variable type
    var_type_pft,                                                             &
! Long name
    "N TAKEN UP FOR GROWTH of existing plant biomass",                        &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_demand_growth
!-----------------------------------------------------------------------------
DATA metadata(374) / var_metadata(                                            &
! String identifier
    'n_demand_growth',                                                        &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Demand for N for GROWTH of existing plant biomass",                      &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_gas
!-----------------------------------------------------------------------------
DATA metadata(375) / var_metadata(                                            &
! String identifier
    'n_gas',                                                                  &
! Variable type
    var_type_sclayer,                                                         &
! Long name
    "Mineralised N Gas Emissions",                                            &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_leach
!-----------------------------------------------------------------------------
DATA metadata(376) / var_metadata(                                            &
! String identifier
    'n_leach',                                                                &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Leached N Losses",                                                       &
! Units
    "kg/m2/s"                                                                 &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_fix_gb
!-----------------------------------------------------------------------------
DATA metadata(377) / var_metadata(                                            &
! String identifier
    'n_fix_gb',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Fixed N",                                                                &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_fix
!-----------------------------------------------------------------------------
DATA metadata(378) / var_metadata(                                            &
! String identifier
    'n_fix',                                                                  &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Fixed N on PFTs",                                                        &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_demand_spread
!-----------------------------------------------------------------------------
DATA metadata(379) / var_metadata(                                            &
! String identifier
    'n_demand_spread',                                                        &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Demand for N for SPREADING plant across gridbox",                        &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_demand_lit
!-----------------------------------------------------------------------------
DATA metadata(380) / var_metadata(                                            &
! String identifier
    'n_demand_lit',                                                           &
! Variable type
    var_type_pft,                                                             &
! Long name
    "N demand of litter: N lost in leaf, wood, root biomass",                 &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for deposition_n
!-----------------------------------------------------------------------------
DATA metadata(381) / var_metadata(                                            &
! String identifier
    'deposition_n',                                                           &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Nitrogen deposition",                                                    &
! Units
    "kg/m2/s"                                                                 &
  ) /
!-----------------------------------------------------------------------------
! Metadata for root_litN
!-----------------------------------------------------------------------------
DATA metadata(382) / var_metadata(                                            &
! String identifier
    'root_litN',                                                              &
! Variable type
    var_type_pft,                                                             &
! Long name
    "N lost as LITTER due to ROOT turnover",                                  &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_veg
!-----------------------------------------------------------------------------
DATA metadata(383) / var_metadata(                                            &
! String identifier
    'n_veg',                                                                  &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Plant NITROGEN content (PFTs): N_LEAF+N_ROOT+N_WOOD from C equivalents", &
! Units
    "kg/m2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_veg_gb
!-----------------------------------------------------------------------------
DATA metadata(384) / var_metadata(                                            &
! String identifier
    'n_veg_gb',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Plant NITROGEN content (GBM): N_LEAF+N_ROOT+N_WOOD from C equivalents",  &
! Units
    "kg/m2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for dcveg
!-----------------------------------------------------------------------------
DATA metadata(385) / var_metadata(                                            &
! String identifier
    'dcveg',                                                                  &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Increment in vegetation C content (PFTs)",                               &
! Units
    "kg/m2/TRIFFID TIMESTEP"                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for dcveg_gb
!-----------------------------------------------------------------------------
DATA metadata(386) / var_metadata(                                            &
! String identifier
    'dcveg_gb',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Increment in vegetation C content (GBM)",                                &
! Units
    "kg/m2/TRIFFID TIMESTEP"                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for dnveg
!-----------------------------------------------------------------------------
DATA metadata(387) / var_metadata(                                            &
! String identifier
    'dnveg',                                                                  &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Increment in vegetation N content (PFTs)",                               &
! Units
    "kg/m2/TRIFFID TIMESTEP"                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for dnveg_gb
!-----------------------------------------------------------------------------
DATA metadata(388) / var_metadata(                                            &
! String identifier
    'dnveg_gb',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Increment in vegetation N content (GBM)",                                &
! Units
    "kg/m2/TRIFFID TIMESTEP"                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_demand_gb
!-----------------------------------------------------------------------------
DATA metadata(389) / var_metadata(                                            &
! String identifier
    'n_demand_gb',                                                            &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Demand for N in total (GBM)",                                            &
! Units
    "kg/m2/360 DAYS"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_demand
!-----------------------------------------------------------------------------
DATA metadata(390) / var_metadata(                                            &
! String identifier
    'n_demand',                                                               &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Demand for N in total (PFTs)",                                           &
! Units
    "kg/m2/360 DAYS"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for wood_litN
!-----------------------------------------------------------------------------
DATA metadata(391) / var_metadata(                                            &
! String identifierN
    'wood_litN',                                                              &
! Variable type
    var_type_pft,                                                             &
! Long name
    "N lost as LITTER due to WOOD turnover",                                  &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for litterC
!-----------------------------------------------------------------------------
DATA metadata(392) / var_metadata(                                            &
! String identifier
    'litterC',                                                                &
! Variable type
    var_type_pft,                                                             &
! Long name
    "LITTER CARBON as the sum of biomass lost through LEAF, ROOT and WOOD turnover", &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lit_n
!-----------------------------------------------------------------------------
DATA metadata(393) / var_metadata(                                            &
! String identifier
    'lit_n',                                                                  &
! Variable type
    var_type_pft,                                                             &
! Long name
    "NITROGEN LITTER (PFTs)",                                                 &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lit_n_t
!-----------------------------------------------------------------------------
DATA metadata(394) / var_metadata(                                            &
! String identifier
    'lit_n_t',                                                                &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Total NITROGEN LITTER",                                                  &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_uptake
!-----------------------------------------------------------------------------
DATA metadata(395) / var_metadata(                                            &
! String identifier
    'n_uptake',                                                               &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Total NITROGEN UPTAKE by plants (PFTs)",                                 &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_uptake_gb
!-----------------------------------------------------------------------------
DATA metadata(396) / var_metadata(                                            &
! String identifier
    'n_uptake_gb',                                                            &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Total NITROGEN UPTAKE by plants (GBM)",                                  &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for burnt_area_gb (INFERNO)
!-----------------------------------------------------------------------------
DATA metadata(397) / var_metadata(                                            &
! String identifier
    'burnt_area_gb',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean burnt area fraction",                                       &
! Units
    "s-1"                                                                     &
  ) /
!-----------------------------------------------------------------------------
! Metadata for burnt_area
!-----------------------------------------------------------------------------
DATA metadata(398) / var_metadata(                                            &
! String identifier
    'burnt_area',                                                             &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT burnt area fraction",                                                &
! Units
    "s-1"                                                                     &
  ) /
!-----------------------------------------------------------------------------
! Metadata for emitted_carbon_gb
!-----------------------------------------------------------------------------
DATA metadata(399) / var_metadata(                                            &
! String identifier
    'emitted_carbon_gb',                                                      &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean emitted carbon",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for emitted_carbon
!-----------------------------------------------------------------------------
DATA metadata(400) / var_metadata(                                            &
! String identifier
    'emitted_carbon',                                                         &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT emitted carbon",                                                     &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for emitted_carbon_DPM
!-----------------------------------------------------------------------------
DATA metadata(401) / var_metadata(                                            &
! String identifier
    'emitted_carbon_DPM',                                                     &
! Variable type
    var_type_surface,                                                         &
! Long name
    "DPM emitted carbon",                                                     &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for emitted_carbon_RPM
!-----------------------------------------------------------------------------
DATA metadata(402) / var_metadata(                                            &
! String identifier
    'emitted_carbon_RPM',                                                     &
! Variable type
    var_type_surface,                                                         &
! Long name
    "RPM emitted carbon",                                                     &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_CO2_gb
!-----------------------------------------------------------------------------
DATA metadata(403) / var_metadata(                                            &
! String identifier
    'fire_em_CO2_gb',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean CO2 emission from fires",                                   &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_CO2
!-----------------------------------------------------------------------------
DATA metadata(404) / var_metadata(                                            &
! String identifier
    'fire_em_CO2',                                                            &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT CO2 emission from fires",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_CO2_DPM
!-----------------------------------------------------------------------------
DATA metadata(405) / var_metadata(                                            &
! String identifier
    'fire_em_CO2_DPM',                                                        &
! Variable type
    var_type_surface,                                                         &
! Long name
    "DPM CO2 emission from fires",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_CO2_RPM
!-----------------------------------------------------------------------------
DATA metadata(406) / var_metadata(                                            &
! String identifier
    'fire_em_CO2_RPM',                                                        &
! Variable type
    var_type_surface,                                                         &
! Long name
    "RPM CO2 emission from fires",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_CO_gb
!-----------------------------------------------------------------------------
DATA metadata(407) / var_metadata(                                            &
! String identifier
    'fire_em_CO_gb',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean CO emission from fires",                                    &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_CO
!-----------------------------------------------------------------------------
DATA metadata(408) / var_metadata(                                            &
! String identifier
    'fire_em_CO',                                                             &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT CO emission from fires",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_CO_DPM
!-----------------------------------------------------------------------------
DATA metadata(409) / var_metadata(                                            &
! String identifier
    'fire_em_CO_DPM',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "DPM CO emission from fires",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_CO_RPM
!-----------------------------------------------------------------------------
DATA metadata(410) / var_metadata(                                            &
! String identifier
    'fire_em_CO_RPM',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "RPM CO emission from fires",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_CH4_gb
!-----------------------------------------------------------------------------
DATA metadata(411) / var_metadata(                                            &
! String identifier
    'fire_em_CH4_gb',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean CH4 emission from fires",                                   &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_CH4
!-----------------------------------------------------------------------------
DATA metadata(412) / var_metadata(                                            &
! String identifier
    'fire_em_CH4',                                                            &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Gridbox mean CH4 emission from fires",                                   &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_CH4_DPM
!-----------------------------------------------------------------------------
DATA metadata(413) / var_metadata(                                            &
! String identifier
    'fire_em_CH4_DPM',                                                        &
! Variable type
    var_type_surface,                                                         &
! Long name
    "DPM CH4 emission from fires",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_CH4_RPM
!-----------------------------------------------------------------------------
DATA metadata(414) / var_metadata(                                            &
! String identifier
    'fire_em_CH4_RPM',                                                        &
! Variable type
    var_type_surface,                                                         &
! Long name
    "RPM CH4 emission from fires",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_NOx_gb
!-----------------------------------------------------------------------------
DATA metadata(415) / var_metadata(                                            &
! String identifier
    'fire_em_NOx_gb',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean NOx emission from fires",                                   &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_NOx
!-----------------------------------------------------------------------------
DATA metadata(416) / var_metadata(                                            &
! String identifier
    'fire_em_NOx',                                                            &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Gridbox mean NOx emission from fires",                                   &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_NOx_DPM
!-----------------------------------------------------------------------------
DATA metadata(417) / var_metadata(                                            &
! String identifier
    'fire_em_NOx_DPM',                                                        &
! Variable type
    var_type_surface,                                                         &
! Long name
    "DPM NOx emission from fires",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_NOx_RPM
!-----------------------------------------------------------------------------
DATA metadata(418) / var_metadata(                                            &
! String identifier
    'fire_em_NOx_RPM',                                                        &
! Variable type
    var_type_surface,                                                         &
! Long name
    "RPM NOx emission from fires",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_SO2_gb
!-----------------------------------------------------------------------------
DATA metadata(419) / var_metadata(                                            &
! String identifier
    'fire_em_SO2_gb',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean SO2 emission from fires",                                   &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_SO2
!-----------------------------------------------------------------------------
DATA metadata(420) / var_metadata(                                            &
! String identifier
    'fire_em_SO2',                                                            &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Gridbox mean SO2 emission from fires",                                   &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_SO2_DPM
!-----------------------------------------------------------------------------
DATA metadata(421) / var_metadata(                                            &
! String identifier
    'fire_em_SO2_DPM',                                                        &
! Variable type
    var_type_surface,                                                         &
! Long name
    "DPM SO2 emission from fires",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_SO2_RPM
!-----------------------------------------------------------------------------
DATA metadata(422) / var_metadata(                                            &
! String identifier
    'fire_em_SO2_RPM',                                                        &
! Variable type
    var_type_surface,                                                         &
! Long name
    "RPM SO2 emission from fires",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_OC_gb
!-----------------------------------------------------------------------------
DATA metadata(423) / var_metadata(                                            &
! String identifier
    'fire_em_OC_gb',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean OC emission from fires",                                    &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_OC
!-----------------------------------------------------------------------------
DATA metadata(424) / var_metadata(                                            &
! String identifier
    'fire_em_OC',                                                             &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT OC emission from fires",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_OC_DPM
!-----------------------------------------------------------------------------
DATA metadata(425) / var_metadata(                                            &
! String identifier
    'fire_em_OC_DPM',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "DPM OC emission from fires",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_OC_RPM
!-----------------------------------------------------------------------------
DATA metadata(426) / var_metadata(                                            &
! String identifier
    'fire_em_OC_RPM',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "RPM OC emission from fires",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_BC_gb
!-----------------------------------------------------------------------------
DATA metadata(427) / var_metadata(                                            &
! String identifier
    'fire_em_BC_gb',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean BC emission from fires",                                    &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_BC
!-----------------------------------------------------------------------------
DATA metadata(428) / var_metadata(                                            &
! String identifier
    'fire_em_BC',                                                             &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT BC emission from fires",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_BC_DPM
!-----------------------------------------------------------------------------
DATA metadata(429) / var_metadata(                                            &
! String identifier
    'fire_em_BC_DPM',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "DPM BC emission from fires",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fire_em_BC_RPM
!-----------------------------------------------------------------------------
DATA metadata(430) / var_metadata(                                            &
! String identifier
    'fire_em_BC_RPM',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "RPM BC emission from fires",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for pop_den
!-----------------------------------------------------------------------------
DATA metadata(431) / var_metadata(                                            &
! String identifier
    'pop_den',                                                                &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Population Density",                                                     &
! Units
    "ppl/km2"                                                                 &
  ) /
!-----------------------------------------------------------------------------
! Metadata for flash_rate
!-----------------------------------------------------------------------------
DATA metadata(432) / var_metadata(                                            &
! String identifier
    'flash_rate',                                                             &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Cloud to Ground lightning flash rate",                                   &
! Units
    "flashes/km2"                                                             &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lit_c_fire
!-----------------------------------------------------------------------------
DATA metadata(433) / var_metadata(                                            &
! String identifier
    'lit_c_fire',                                                             &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT carbon litter due to fire",                                          &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for burnt_carbon_dpm
!-----------------------------------------------------------------------------
DATA metadata(434) / var_metadata(                                            &
! String identifier
    'burnt_carbon_dpm',                                                       &
! Variable type
    var_type_surface,                                                         &
! Long name
    "DPM burnt carbon due to fire",                                           &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lit_n_fire
!-----------------------------------------------------------------------------
DATA metadata(435) / var_metadata(                                            &
! String identifier
    'lit_n_fire',                                                             &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT nitrogen litter due to fire",                                        &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for burnt_carbon_rpm
!-----------------------------------------------------------------------------
DATA metadata(436) / var_metadata(                                            &
! String identifier
    'burnt_carbon_rpm',                                                       &
! Variable type
    var_type_surface,                                                         &
! Long name
    "RPM burnt carbon due to fire",                                           &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for veg_c_fire_emission_gb
!-----------------------------------------------------------------------------
DATA metadata(437) / var_metadata(                                            &
! String identifier
    'veg_c_fire_emission_gb',                                                 &
! Variable type
    var_type_surface,                                                         &
! Long name
    "GBM carbon emissions from fire",                                         &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for veg_c_fire_emission_pft
!-----------------------------------------------------------------------------
DATA metadata(438) / var_metadata(                                            &
! String identifier
    'veg_c_fire_emission_pft',                                                &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT carbon emissions from fire",                                         &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for g_burn_pft
!-----------------------------------------------------------------------------
DATA metadata(439) / var_metadata(                                            &
! String identifier
    'g_burn_pft',                                                             &
! Variable type
    var_type_pft,                                                             &
! Long name
    "Fire disturbance rate",                                                  &
! Units
    "m2/m2/360days"                                                           &
  ) /
!-----------------------------------------------------------------------------
! Metadata for surface storage on river points
!-----------------------------------------------------------------------------
DATA metadata(440) / var_metadata(                                            &
! String identifier
    'rfm_surfstore_rp',                                                       &
! Variable type
    var_type_rp,                                                              &
! Long name
    "Surface storage on river points",                                        &
! Units
    "m3"                                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sub-surface storage on river points
!-----------------------------------------------------------------------------
DATA metadata(441) / var_metadata(                                            &
! String identifier
    'rfm_substore_rp',                                                        &
! Variable type
    var_type_rp,                                                              &
! Long name
    "Sub-surface storage on river points",                                    &
! Units
    "m3"                                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for surface inflow on river points
!-----------------------------------------------------------------------------
DATA metadata(442) / var_metadata(                                            &
! String identifier
    'rfm_flowin_rp',                                                          &
! Variable type
    var_type_rp,                                                              &
! Long name
    "Surface inflow on river points",                                         &
! Units
    "m3"                                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sub-surface inflow on river points
!-----------------------------------------------------------------------------
DATA metadata(443) / var_metadata(                                            &
! String identifier
    'rfm_bflowin_rp',                                                         &
! Variable type
    var_type_rp,                                                              &
! Long name
    "Sub-surface inflow on river points",                                     &
! Units
    "m3"                                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cnsrv_nitrogen_triffid
!-----------------------------------------------------------------------------
DATA metadata(444) / var_metadata(                                            &
! String identifier
    'cnsrv_nitrogen_triffid',                                                 &
! Variable type
    var_type_surface,                                                         &
! Long name
    "error in land nitrogen conservation in triffid routine",                 &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cnsrv_vegN_triffid
!-----------------------------------------------------------------------------
DATA metadata(445) / var_metadata(                                            &
! String identifier
    'cnsrv_vegN_triffid',                                                     &
! Variable type
    var_type_surface,                                                         &
! Long name
    "error in vegetation nitrogen conservation in triffid routine",           &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cnsrv_soilN_triffid
!-----------------------------------------------------------------------------
DATA metadata(446) / var_metadata(                                            &
! String identifier
    'cnsrv_soilN_triffid',                                                    &
! Variable type
    var_type_surface,                                                         &
! Long name
    "error in soil nitrogen conservation in triffid routine",                 &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cnsrv_n_inorg_triffid
!-----------------------------------------------------------------------------
DATA metadata(447) / var_metadata(                                            &
! String identifier
    'cnsrv_n_inorg_triffid',                                                  &
! Variable type
    var_type_surface,                                                         &
! Long name
    "error in inorganic nitrogen conservation in triffid routine",            &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for root_abandon_gb
!-----------------------------------------------------------------------------
DATA metadata(448) / var_metadata(                                            &
! String identifier
    'root_abandon_gb',                                                        &
! Variable type
    var_type_surface,                                                         &
! Long name
    "gridbox mean C flux from roots to soil due to landuse change",           &
! Units
    "kg C  m-2 per 360days"                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for root_abandon_n_gb
!-----------------------------------------------------------------------------
DATA metadata(449) / var_metadata(                                            &
! String identifier
    'root_abandon_n_gb',                                                      &
! Variable type
    var_type_surface,                                                         &
! Long name
    "gridbox mean N flux from roots to soil due to landuse change",           &
! Units
    "kg N  m-2 per 360days"                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for resp_s_to_atmos_gb
!-----------------------------------------------------------------------------
DATA metadata(450) / var_metadata(                                            &
! String identifier
    'resp_s_to_atmos_gb',                                                     &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean soil respiration carbon flux to atmosphere",                &
! Units
    "kg m-2 per 360days"                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_gas_gb
!-----------------------------------------------------------------------------
DATA metadata(451) / var_metadata(                                            &
! String identifier
    'n_gas_gb',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "gridbox mean mineralised N Gas Emissions",                               &
! Units
    "kg/m2/360 days"                                                          &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_inorg_gb
!-----------------------------------------------------------------------------
DATA metadata(452) / var_metadata(                                            &
! String identifier
    'n_inorg_gb',                                                             &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox total N in soil inorganic species",                              &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Variables available with any multi-pool soil C model (e.g. RothC, ECOSSE).
!-----------------------------------------------------------------------------
! Soil prognostic variables and single pools from multi-pool prognostic
! variables.
!-----------------------------------------------------------------------------
! Metadata for c_bio
!-----------------------------------------------------------------------------
DATA metadata(453) / var_metadata(                                            &
! String identifier
    'c_bio',                                                                  &
! Variable type
    var_type_sclayer,                                                         &
! Long name
    "C in soil biomass",                                                      &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for c_dpm
!-----------------------------------------------------------------------------
DATA metadata(454) / var_metadata(                                            &
! String identifier
    'c_dpm',                                                                  &
! Variable type
    var_type_sclayer,                                                         &
! Long name
    "C in decomposable plant material",                                       &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for c_hum
!-----------------------------------------------------------------------------
DATA metadata(455) / var_metadata(                                            &
! String identifier
    'c_hum',                                                                  &
! Variable type
    var_type_sclayer,                                                         &
! Long name
    "C in humus",                                                             &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for c_rpm
!-----------------------------------------------------------------------------
DATA metadata(456) / var_metadata(                                            &
! String identifier
    'c_rpm',                                                                  &
! Variable type
    var_type_sclayer,                                                         &
! Long name
    "C in resistant plant material",                                          &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_bio
!-----------------------------------------------------------------------------
DATA metadata(457) / var_metadata(                                            &
! String identifier
    'n_bio',                                                                  &
! Variable type
    var_type_sclayer,                                                         &
! Long name
    "N in soil biomass",                                                      &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_dpm
!-----------------------------------------------------------------------------
DATA metadata(458) / var_metadata(                                            &
! String identifier
    'n_dpm',                                                                  &
! Variable type
    var_type_sclayer,                                                         &
! Long name
    "N in decomposable plant material",                                       &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_hum
!-----------------------------------------------------------------------------
DATA metadata(459) / var_metadata(                                            &
! String identifier
    'n_hum',                                                                  &
! Variable type
    var_type_sclayer,                                                         &
! Long name
    "N in humus",                                                             &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_rpm
!-----------------------------------------------------------------------------
DATA metadata(460) / var_metadata(                                            &
! String identifier
    'n_rpm',                                                                  &
! Variable type
    var_type_sclayer,                                                         &
! Long name
    "N in resistant plant material",                                          &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Gridbox totals of multi-pool soil prognostic variables.
!-----------------------------------------------------------------------------
! Metadata for c_bio_gb
!-----------------------------------------------------------------------------
DATA metadata(461) / var_metadata(                                            &
! String identifier
    'c_bio_gb',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "C in soil biomass, gridbox total",                                       &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for c_dpm_gb
!-----------------------------------------------------------------------------
DATA metadata(462) / var_metadata(                                            &
! String identifier
    'c_dpm_gb',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "C in decomposable plant material, gridbox total",                        &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for c_hum_gb
!-----------------------------------------------------------------------------
DATA metadata(463) / var_metadata(                                            &
! String identifier
    'c_hum_gb',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "C in soil humus, gridbox total",                                         &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for c_rpm_gb
!-----------------------------------------------------------------------------
DATA metadata(464) / var_metadata(                                            &
! String identifier
    'c_rpm_gb',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "C in resistant plant material, gridbox total",                           &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_bio_gb
!-----------------------------------------------------------------------------
DATA metadata(465) / var_metadata(                                            &
! String identifier
    'n_bio_gb',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "N in soil biomass, gridbox total",                                       &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_dpm_gb
!-----------------------------------------------------------------------------
DATA metadata(466) / var_metadata(                                            &
! String identifier
    'n_dpm_gb',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "N in decomposable plant material, gridbox total",                        &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_hum_gb
!-----------------------------------------------------------------------------
DATA metadata(467) / var_metadata(                                            &
! String identifier
    'n_hum_gb',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "N in soil humus, gridbox total",                                         &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_rpm_gb
!-----------------------------------------------------------------------------
DATA metadata(468) / var_metadata(                                            &
! String identifier
    'n_rpm_gb',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "N in resistant plant material, gridbox total",                           &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Gridbox totals of derived variables for multi-pool soil models.
!-----------------------------------------------------------------------------
! Metadata for n_soil_gb
!-----------------------------------------------------------------------------
DATA metadata(469) / var_metadata(                                            &
! String identifier
    'n_soil_gb',                                                              &
! Variable type
    var_type_surface,                                                         &
! Long name
    "N in soil, all forms, gridbox total",                                    &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! ECOSSE: components of prognostic variables.
!-----------------------------------------------------------------------------
! Metadata for n_amm.
!-----------------------------------------------------------------------------
DATA metadata(470) / var_metadata(                                            &
! String identifier
    'n_amm',                                                                  &
! Variable type
    var_type_sclayer,                                                         &
! Long name
    "N in soil ammonium",                                                     &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_nit
!-----------------------------------------------------------------------------
DATA metadata(471) / var_metadata(                                            &
! String identifier
    'n_nit',                                                                  &
! Variable type
    var_type_sclayer,                                                         &
! Long name
    "N in soil nitrate",                                                      &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Gridbox totals of components of ECOSSE prognostic variables.
!-----------------------------------------------------------------------------
! Metadata for n_amm_gb
!-----------------------------------------------------------------------------
DATA metadata(472) / var_metadata(                                            &
! String identifier
    'n_amm_gb',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "N in soil ammonium, gridbox total",                                      &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_nit_gb
!-----------------------------------------------------------------------------
DATA metadata(473) / var_metadata(                                            &
! String identifier
    'n_nit_gb',                                                               &
! Variable type
    var_type_surface,                                                         &
! Long name
    "N in soil nitrate, gridbox total",                                       &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Soil ancillaries on soil carbon layers.
!-----------------------------------------------------------------------------
! Metadata for soil_ph
!-----------------------------------------------------------------------------
DATA metadata(474) / var_metadata(                                            &
! String identifier
    'soil_ph',                                                                &
! Variable type
    var_type_sclayer,                                                         &
! Long name
    "Soil pH",                                                                &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for soil_ph_const_z
!-----------------------------------------------------------------------------
DATA metadata(475) / var_metadata(                                            &
! String identifier
    'soil_ph_const_z',                                                        &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name
! and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for slope
!-----------------------------------------------------------------------------
DATA metadata(476) / var_metadata(                                            &
! String identifier
    'slope',                                                                  &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Terrain slope",                                                          &
! Units
    "degrees"                                                                 &
  ) /
!-----------------------------------------------------------------------------
! Metadata for b_soilt
!-----------------------------------------------------------------------------
DATA metadata(477) / var_metadata(                                            &
! String identifier
    'b_soilt',                                                                &
! Variable type
    var_type_soilt_soil,                                                      &
! Long name
    "Soil tiled Brooks-Corey exponent for each soil layer",                   &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for b_const_z_soilt
!-----------------------------------------------------------------------------
DATA metadata(478) / var_metadata(                                            &
! String identifier
    'b_const_z_soilt',                                                        &
! Variable type
    var_type_soilt,                                                           &
! Variable is not available for output, so give dummy values for
! long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sathh_soilt
!-----------------------------------------------------------------------------
DATA metadata(479) / var_metadata(                                            &
! String identifier
    'sathh_soilt',                                                            &
! Variable type
    var_type_soilt_soil,                                                      &
! Long name
    "Soil tiled saturated soil water pressure for each soil layer",           &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sathh_const_z_soilt
!-----------------------------------------------------------------------------
DATA metadata(480) / var_metadata(                                            &
! String identifier
    'sathh_const_z_soilt',                                                    &
! Variable type
    var_type_soilt,                                                           &
! Variable is not available for output, so give dummy values for
! long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for satcon_soilt
!-----------------------------------------------------------------------------
DATA metadata(481) / var_metadata(                                            &
! String identifier
    'satcon_soilt',                                                           &
! Variable type
    var_type_soilt_soil,                                                      &
! Long name
    "Soil tiled saturated hydraulic conductivity for each soil layer",        &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for satcon_const_z_soilt
!-----------------------------------------------------------------------------
DATA metadata(482) / var_metadata(                                            &
! String identifier
    'satcon_const_z_soilt',                                                   &
! Variable type
    var_type_soilt,                                                           &
! Variable is not available for output, so give dummy values for
! long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sm_sat_soilt
!-----------------------------------------------------------------------------
DATA metadata(483) / var_metadata(                                            &
! String identifier
    'sm_sat_soilt',                                                           &
! Variable type
    var_type_soilt_soil,                                                      &
! Long name
    "Soil tiled volumetric moisture content at saturation" //                 &
    "for each soil layer",                                                    &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sm_sat_const_z_soilt
!-----------------------------------------------------------------------------
DATA metadata(484) / var_metadata(                                            &
! String identifier
    'sm_sat_const_z_soilt',                                                   &
! Variable type
    var_type_soilt,                                                           &
! Variable is not available for output, so give dummy values for
! long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sm_crit_soilt
!-----------------------------------------------------------------------------
DATA metadata(485) / var_metadata(                                            &
! String identifier
    'sm_crit_soilt',                                                          &
! Variable type
    var_type_soilt_soil,                                                      &
! Long name
    "Soil tiled volumetric moisture content at critical point" //             &
    "for each soil layer",                                                    &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sm_crit_const_z_soilt
!-----------------------------------------------------------------------------
DATA metadata(486) / var_metadata(                                            &
! String identifier
    'sm_crit_const_z_soilt',                                                  &
! Variable type
    var_type_soilt,                                                           &
! Variable is not available for output, so give dummy values for
! long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sm_wilt_soilt
!-----------------------------------------------------------------------------
DATA metadata(487) / var_metadata(                                            &
! String identifier
    'sm_wilt_soilt',                                                          &
! Variable type
    var_type_soilt_soil,                                                      &
! Long name
    "Soil tiled volumetric moisture content at wilting point" //              &
    "for each soil layer",                                                    &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sm_wilt_const_z_soilt
!-----------------------------------------------------------------------------
DATA metadata(488) / var_metadata(                                            &
! String identifier
    'sm_wilt_const_z_soilt',                                                  &
! Variable type
    var_type_soilt,                                                           &
! Variable is not available for output, so give dummy values for
! long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for hcap_soilt
!-----------------------------------------------------------------------------
DATA metadata(489) / var_metadata(                                            &
! String identifier
    'hcap_soilt',                                                             &
! Variable type
    var_type_soilt_soil,                                                      &
! Long name
    "Soil tiled dry soil heat capacity for each soil layer",                  &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for hcap_const_z_soilt
!-----------------------------------------------------------------------------
DATA metadata(490) / var_metadata(                                            &
! String identifier
    'hcap_const_z_soilt',                                                     &
! Variable type
    var_type_soilt,                                                           &
! Variable is not available for output, so give dummy values for
! long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for hcon_soilt
!-----------------------------------------------------------------------------
DATA metadata(491) / var_metadata(                                            &
! String identifier
    'hcon_soilt',                                                             &
! Variable type
    var_type_soilt_soil,                                                      &
! Long name
    "Soil tiled dry soil thermal conductivity for each soil layer",           &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for hcon_const_z_soilt
!-----------------------------------------------------------------------------
DATA metadata(492) / var_metadata(                                            &
! String identifier
    'hcon_const_z_soilt',                                                     &
! Variable type
    var_type_soilt,                                                           &
! Variable is not available for output, so give dummy values for
! long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fexp_soilt
!-----------------------------------------------------------------------------
DATA metadata(493) / var_metadata(                                            &
! String identifier
    'fexp_soilt',                                                             &
! Variable type
    var_type_soilt,                                                           &
! Variable is not available for output, so give dummy values for
! long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for ti_mean_soilt
!-----------------------------------------------------------------------------
DATA metadata(494) / var_metadata(                                            &
! String identifier
    'ti_mean_soilt',                                                          &
! Variable type
    var_type_soilt,                                                           &
! Variable is not available for output, so give dummy values for
! long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for ti_sig_soilt
!-----------------------------------------------------------------------------
DATA metadata(495) / var_metadata(                                            &
! String identifier
    'ti_sig_soilt',                                                           &
! Variable type
    var_type_soilt,                                                           &
! Variable is not available for output, so give dummy values for
! long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for flammability
!-----------------------------------------------------------------------------
DATA metadata(496) / var_metadata(                                            &
! String identifier
    'flammability',                                                           &
! Variable type
    var_type_pft,                                                             &
! Long name
    "PFT flammability",                                                       &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for g_burn_gb
!-----------------------------------------------------------------------------
DATA metadata(497) / var_metadata(                                            &
! String identifier
    'g_burn_gb',                                                              &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Fire disturbance rate by gridbox",                                       &
! Units
    "m2/m2/360days"                                                           &
  ) /
!-----------------------------------------------------------------------------
! Metadata for soil_ph_soilt
!-----------------------------------------------------------------------------
DATA metadata(498) / var_metadata(                                            &
! String identifier
    'soil_ph_soilt',                                                          &
! Variable type
    var_type_soilt_sclayer,                                                   &
! Long name
    "Soil pH on soil tiles",                                                  &
! Units
    "1"                                                                       &
  ) /  
!-----------------------------------------------------------------------------
! Metadata for drain_soilt
!-----------------------------------------------------------------------------
DATA metadata(499) / var_metadata(                                            &
! String identifier
    'drain_soilt',                                                            &
! Variable type
    var_type_soilt,                                                           &
! Long name
    "Drainage from bottom (nshyd) soil layer on soil tiles",                  &
! Units
    "kg m-2 s-1"                                                              &
  ) /  
!-----------------------------------------------------------------------------
! Metadata for fch4_wetl_cs_soilt
!-----------------------------------------------------------------------------
DATA metadata(500) / var_metadata(                                            &
! String identifier
    'fch4_wetl_cs_soilt',                                                     &
! Variable type
    var_type_soilt,                                                           &
! Long name
    "Scaled methane flux from wetland fraction (cs substrate) on soil tiles", &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fch4_wetl_npp_soilt
!-----------------------------------------------------------------------------
DATA metadata(501) / var_metadata(                                            &
! String identifier
    'fch4_wetl_npp_soilt',                                                    &
! Variable type
    var_type_soilt,                                                           &
! Long name
    "Scaled methane flux from wetland fraction (npp substrate) on soil tiles",&
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fch4_wetl_resps_soilt
!-----------------------------------------------------------------------------
DATA metadata(502) / var_metadata(                                            &
! String identifier
    'fch4_wetl_resps_soilt',                                                  &
! Variable type
    var_type_soilt,                                                           &
! Long name
    "Scaled methane flux from wetland fraction (soil resp substrate) " //     &
    "on soil tiles",                                                          &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fsat_soilt
!-----------------------------------------------------------------------------
DATA metadata(503) / var_metadata(                                            &
! String identifier
    'fsat_soilt',                                                             &
! Variable type
    var_type_soilt,                                                           &
! Long name
    "Surface saturated fraction on soil tiles",                               &
! Units
    "1"                                                                       &
  ) /  
!-----------------------------------------------------------------------------
! Metadata for fwetl_soilt
!-----------------------------------------------------------------------------
DATA metadata(504) / var_metadata(                                            &
! String identifier
    'fwetl_soilt',                                                            &
! Variable type
    var_type_soilt,                                                           &
! Long name
    "Wetland fraction at end of model timestep on soil tiles.",               &
! Units
    "1"                                                                       &
  ) /  
!-----------------------------------------------------------------------------
! Metadata for qbase_soilt
!-----------------------------------------------------------------------------
DATA metadata(505) / var_metadata(                                            &
! String identifier
    'qbase_soilt',                                                            &
! Variable type
    var_type_soilt,                                                           &
! Long name
    "Baseflow (lateral subsurface runoff) on soil tiles",                     &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for qbase_zw_soilt
!-----------------------------------------------------------------------------
DATA metadata(506) / var_metadata(                                            &
! String identifier
    'qbase_zw_soilt',                                                         &
! Variable type
    var_type_soilt,                                                           &
! Long name
    "Baseflow from deep LSH/TOPMODEL layer on soil tiles",                    &
! Units
    "kg m-2 s-1"                                                              &
  ) /  
!-----------------------------------------------------------------------------
! Metadata for sat_excess_roff_soilt
!-----------------------------------------------------------------------------
DATA metadata(507) / var_metadata(                                            &
! String identifier
    'sat_excess_roff_soilt',                                                  &
! Variable type
    var_type_soilt,                                                           &
! Long name
    "Saturation excess surface ('Dunne') runoff on soil tiles",               &
! Units
    "kg m-2 s-1"                                                              &
  ) /  
!-----------------------------------------------------------------------------
! Metadata for zw_soilt
!-----------------------------------------------------------------------------
DATA metadata(508) / var_metadata(                                            &
! String identifier
    'zw_soilt',                                                               &
! Variable type
    var_type_soilt,                                                           &
! Long name
    "Mean depth to water table on soil tiles",                                &
! Units
    "m"                                                                       &
  ) /  
!-----------------------------------------------------------------------------
! Metadata for ext_soilt
!-----------------------------------------------------------------------------
DATA metadata(509) / var_metadata(                                            &
! String identifier
    'ext_soilt',                                                              &
! Variable type
    var_type_soilt_soil,                                                      &
! Long name
    "Extraction of water from each soil layer on soil tiles",                 &
! Units
    "kg m-2 s-1"                                                              &
  ) /  
!-----------------------------------------------------------------------------
! Metadata for smcl_soilt
!-----------------------------------------------------------------------------
DATA metadata(510) / var_metadata(                                            &
! String identifier
    'smcl_soilt',                                                             &
! Variable type
    var_type_soilt_soil,                                                      &
! Long name
    "Moisture content of each soil layer on soil tiles",                      &
! Units
    "kg m-2"                                                                  &
  ) /  
!-----------------------------------------------------------------------------
! Metadata for sthf_soilt
!-----------------------------------------------------------------------------
DATA metadata(511) / var_metadata(                                            &
! String identifier
    'sthf_soilt',                                                             &
! Variable type
    var_type_soilt_soil,                                                      &
! Long name
    "Frozen moisture content of each soil layer as a fraction of " //         &
    "saturation on soil tiles",                                               &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sthu_soilt
!-----------------------------------------------------------------------------
DATA metadata(512) / var_metadata(                                            &
! String identifier
    'sthu_soilt',                                                             &
! Variable type
    var_type_soilt_soil,                                                      &
! Long name
    "Unfrozen moisture content of each soil layer as a fraction of "  //      &
    "saturation on soil tiles",                                               &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for t_soil_soilt
!-----------------------------------------------------------------------------
DATA metadata(513) / var_metadata(                                            &
! String identifier
    't_soil_soilt',                                                           &
! Variable type
    var_type_soilt_soil,                                                      &
! Long name
    "Sub-surface temperature of each layer on soil tiles",                    &
! Units
    "K"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for frac_irrig_soilt
!-----------------------------------------------------------------------------
DATA metadata(514) / var_metadata(                                            &
! String identifier
    'frac_irrig_soilt',                                                       &
! Variable type
    var_type_soilt,                                                           &
! Long name
    "Irrigated fraction on soil tiles",                                       &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sthu_irr_soilt
!-----------------------------------------------------------------------------
DATA metadata(515) / var_metadata(                                            &
! String identifier
    'sthu_irr_soilt',                                                         &
! Variable type
    var_type_soilt_soil,                                                      &
! Long name
    "Wetness of each soil layer over irrigation on soil tiles",               &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for croplatestharvdate
!-----------------------------------------------------------------------------
DATA metadata(516) / var_metadata(                                            &
! String identifier
    'croplatestharvdate',                                                     &
! Variable type
    var_type_cpft,                                                            &
! Long name
    "Latest day of year when crop can be harvested",                          &
! Units
    "day of year"                                                             &
  ) /
!-----------------------------------------------------------------------------
! Metadata for ftemp
!-----------------------------------------------------------------------------
DATA metadata(517) / var_metadata(                                            &
! String identifier
    'ftemp',                                                                  &
! Variable type
    var_type_soil,                                                            &
! Long name
    "Soil respiration rate modifier due to soil temperature",                 &
! Units
    "-"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fsth
!-----------------------------------------------------------------------------
DATA metadata(518) / var_metadata(                                            &
! String identifier
    'fsth',                                                                   &
! Variable type
    var_type_soil,                                                            &
! Long name
    "Soil respiration rate modifier due to soil moisture",                    &
! Units
    "-"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for fprf
!-----------------------------------------------------------------------------
DATA metadata(519) / var_metadata(                                            &
! String identifier
    'fprf',                                                                   &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Soil respiration rate modifier due to vegetation cover",                 &
! Units
    "-"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for frac_fplain_lp
!-----------------------------------------------------------------------------
DATA metadata(520) / var_metadata(                                            &
! String identifier
    'frac_fplain_lp',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Overbank inundation area as a fraction of gridcell area",                &
! Units
    "(fraction of gridcell)"                                                  &
  )                                                                           &
/
!-----------------------------------------------------------------------------
! Metadata for logn_mean
!-----------------------------------------------------------------------------
DATA metadata(521) / var_metadata(                                            &
! String identifier
    'logn_mean',                                                              &
! Variable type
    var_type_surface,                                                         &
! Long name
    "River routing logn_mean",                                                &
! Units
    "ln(metres)"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for logn_stdev
!-----------------------------------------------------------------------------
DATA metadata(522) / var_metadata(                                            &
! String identifier
    'logn_stdev',                                                             &
! Variable type
    var_type_surface,                                                         &
! Long name
    "River routing logn_stdev",                                               &
! Units
    "ln(metres)"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for clay_const_z
!-----------------------------------------------------------------------------
DATA metadata(523) / var_metadata(                                            &
! String identifier
    'clay_const_z',                                                           &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for clay_soilt
!-----------------------------------------------------------------------------
DATA metadata(524) / var_metadata(                                            &
! String identifier
    'clay_soilt',                                                             &
! Variable type
    var_type_soilt_sclayer,                                                   &
! Long name
    "Soil tiled soil clay content",                                           &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for clay_const_z_soilt
!-----------------------------------------------------------------------------
DATA metadata(525) / var_metadata(                                            &
! String identifier
    'clay_const_z_soilt',                                                     &
! Variable type
    var_type_soilt,                                                           &
! Variable is not available for output, so give dummy values for
! long name and units
    "", ""                                                                    &
  ) /
!-----------------------------------------------------------------------------
! Metadata for plant_input_c_gb.
!-----------------------------------------------------------------------------
DATA metadata(526) / var_metadata(                                            &
! String identifier
    'plant_input_c_gb',                                                       &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Plant input of carbon to soil by litterfall",                            &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for plant_input_n_gb.
!-----------------------------------------------------------------------------
DATA metadata(527) / var_metadata(                                            &
! String identifier
    'plant_input_n_gb',                                                       &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Plant input of nitrogen to soil by litterfall",                          &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for imogen d_land_atmos_co2
!-----------------------------------------------------------------------------
DATA metadata(528) / var_metadata(                                            &
! String identifier
    'd_land_atmos_co2',                                                       &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Change in atmos CO2 conc from land-atmos feedbacks in IMOGEN" //         &
    "(repeated global value)",                                                &
! Units
    "ppm/year"                                                                &
  ) /
!-----------------------------------------------------------------------------
! Metadata for imogen d_ocean_atmos
!-----------------------------------------------------------------------------
DATA metadata(529) / var_metadata(                                            &
! String identifier
    'd_ocean_atmos',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Change in atmos CO2 conc from ocean feedbacks in IMOGEN" //              &
     "(repeated global value)",                                               &
! Units
    "ppm/year"                                                                &
  ) /
!-----------------------------------------------------------------------------
! Metadata for c_emiss_out
!-----------------------------------------------------------------------------
DATA metadata(530) / var_metadata(                                            &
! String identifier
    'c_emiss_out',                                                            &
! Variable type
    var_type_surface,                                                         &
! Long name
    "prescribed carbon emissions in IMOGEN" //                                &
    "(repeated global value)",                                                &
! Units
    "Gt/year"                                                                 &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_denitrif_gb.
!-----------------------------------------------------------------------------
DATA metadata(531) / var_metadata(                                            &
! String identifier
    'n_denitrif_gb',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox rate of denitrification, expressed as N",                        &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_nitrif_gb.
!-----------------------------------------------------------------------------
DATA metadata(532) / var_metadata(                                            &
! String identifier
    'n_nitrif_gb',                                                            &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox rate of nitrification, expressed as N",                          &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for no_soil_gb.
!-----------------------------------------------------------------------------
DATA metadata(533) / var_metadata(                                            &
! String identifier
    'no_soil_gb',                                                             &
! Variable type
    var_type_surface,                                                         &
! Long name
    "N in NO flux from soil to atmosphere",                                   &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n2_denitrif_gb.
!-----------------------------------------------------------------------------
DATA metadata(534) / var_metadata(                                            &
! String identifier
    'n2_denitrif_gb',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "N in N2 lost from soil during denitrification",                          &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n2o_denitrif_gb.
!-----------------------------------------------------------------------------
DATA metadata(535) / var_metadata(                                            &
! String identifier
    'n2o_denitrif_gb',                                                        &
! Variable type
    var_type_surface,                                                         &
! Long name
    "N in N2O lost during denitrification",                                   &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n2o_nitrif_gb.
!-----------------------------------------------------------------------------
DATA metadata(536) / var_metadata(                                            &
! String identifier
    'n2o_nitrif_gb',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "N in N2O lost during nitrification, including partial nitrification",    &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n2o_part_nitrif_gb.
!-----------------------------------------------------------------------------
DATA metadata(537) / var_metadata(                                            &
! String identifier
    'n2o_part_nitrif_gb',                                                     &
! Variable type
    var_type_surface,                                                         &
! Long name
    "N in N2O lost by partial nitrification",                                 &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n2o_soil_gb.
!-----------------------------------------------------------------------------
DATA metadata(538) / var_metadata(                                            &
! String identifier
    'n2o_soil_gb',                                                            &
! Variable type
    var_type_surface,                                                         &
! Long name
    "N in N2O flux from soil to atmosphere",                                  &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_leach_amm_gb.
!-----------------------------------------------------------------------------
DATA metadata(539) / var_metadata(                                            &
! String identifier
    'n_leach_amm_gb',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox N lost through leaching of soil ammonium",                       &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_leach_nit_gb.
!-----------------------------------------------------------------------------
DATA metadata(540) / var_metadata(                                            &
! String identifier
    'n_leach_nit_gb',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox N lost through leaching of soil nitrate",                        &
! Units
    "kg m-2 s-1"                                                              &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lai_gb.
!-----------------------------------------------------------------------------
DATA metadata(541) / var_metadata(                                            &
! String identifier
    'lai_gb',                                                                 &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox leaf area index",                                                &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for apar.
!-----------------------------------------------------------------------------
DATA metadata(542) / var_metadata(                                            &
! String identifier
    'apar',                                                                   &
! Variable type
    var_type_surface,                                                         &
! Long name
    "PFT Absorbed Photosynthetically Active Radiation",                       &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for apar_gb.
!-----------------------------------------------------------------------------
DATA metadata(543) / var_metadata(                                            &
! String identifier
    'apar_gb',                                                                &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox Absorbed Photosynthetically Active Radiation",                   &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for bl_height.
!-----------------------------------------------------------------------------
DATA metadata(544) / var_metadata(                                            &
! String identifier
    'bl_height',                                                              &
! Variable type
    var_type_surface,                                                         &
! Variable is not available for output, but giving long name and units for
! clarity.
! Long name
    "Height above surface of top of boundary layer",                          &
! Units
    "m"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for level_separation.
!-----------------------------------------------------------------------------
DATA metadata(545) / var_metadata(                                            &
! String identifier
    'level_separation',                                                       &
! Variable type
    var_type_bl_level,                                                        &
! Variable is not available for output, but giving long name and units for
! clarity.
! Long name
    "Vertical separation between atmospheric boundary layer levels",          &
! Units
    "m"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for tracer.
!-----------------------------------------------------------------------------
DATA metadata(546) / var_metadata(                                            &
! String identifier
    'tracer_field',                                                           &
! Variable type
    var_type_tracer,                                                          &
! Long name
    "Surface chemical tracer concentrations as mass mixing ratios",           &
! Units
    "kg kg-1"                                                                 &
  ) /
!-----------------------------------------------------------------------------
! Metadata for tau_gb
!-----------------------------------------------------------------------------
DATA metadata(547) / var_metadata(                                            &
! String identifier
    'tau_gb',                                                                 &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Gridbox mean surface wind stress",                                       &
! Units
    "N m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for tau
!-----------------------------------------------------------------------------
DATA metadata(548) / var_metadata(                                            &
! String identifier
    'tau',                                                                    &
! Variable type
    var_type_surft,                                                           &
! Long name
    "Tile surface stress for land tiles",                                     &
! Units
    "N m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_inorg_avail_pft
!-----------------------------------------------------------------------------
DATA metadata(549) / var_metadata(                                            &
! String identifier
    'n_inorg_avail_pft',                                                      &
! Variable type
    var_type_pft_sclayer,                                                     &
! Long name
    "Inorganic N pool that is available for plant uptake",                    &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for metstat temp_ave_nday
!-----------------------------------------------------------------------------
DATA metadata(550) / var_metadata(                                            &
! String identifier
    'temp_ave_nday',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name - this is dynamic and is reset in get_var_attrs.
! This variable is not available as a diagnostic.
    "Average temperature (exp filter) over previous n_day_photo_acclim days", &
! Units
    "Kelvin"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for t_growth_gb
!-----------------------------------------------------------------------------
DATA metadata(551) / var_metadata(                                            &
! String identifier
    't_growth_gb',                                                            &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Growth (average) temperature",                                           &
! Units
    "K"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for n_soil
!-----------------------------------------------------------------------------
DATA metadata(552) / var_metadata(                                            &
! String identifier
    'n_soil',                                                                 &
! Variable type
    var_type_soil_n_pool,                                                     &
! Long name
    "Nitrogen in soil pools",                                                 &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cropfrac
!-----------------------------------------------------------------------------
DATA metadata(553) / var_metadata(                                            &
! String identifier
    'cropfrac',                                                               &
! Variable type
    var_type_cpft,                                                            &
! Long name
    "Fractional cover of each crop type",                                     &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for irrigation fraction on irrigated tiles
!-----------------------------------------------------------------------------
DATA metadata(554) / var_metadata(                                            &
! String identifier
    'irrfrac_irrtiles',                                                       &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Irrigated fraction on irrigated tiles",                                  &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for SoilTemp_CABLE 
!-----------------------------------------------------------------------------
DATA metadata(555) / var_metadata(                                            &
! String identifier
    'SoilTemp_CABLE',                                                         &
! Variable type
    var_type_cable_soil,                                                      &
! Long name
    "PFT soil temperature in layers",                                         &
! Units
    "K"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for FrozenSoilFrac_CABLE 
!-----------------------------------------------------------------------------
DATA metadata(556) / var_metadata(                                            &
! String identifier
    'FrozenSoilFrac_CABLE',                                                   &
! Variable type
    var_type_cable_soil,                                                      &
! Long name
    "PFT frozen soil in layers",                                              &
! Units
    "fraction"                                                                &
  ) /
!-----------------------------------------------------------------------------
! Metadata for SoilMoisture_CABLE 
!-----------------------------------------------------------------------------
DATA metadata(557) / var_metadata(                                            &
! String identifier
    'SoilMoisture_CABLE',                                                     &
! Variable type
    var_type_cable_soil,                                                      &
! Long name
    "PFT soil moisture in layers",                                            &
! Units
    ""                                                                        &
  ) /
!-----------------------------------------------------------------------------
! Metadata for SnowDepth_CABLE 
!-----------------------------------------------------------------------------
DATA metadata(558) /  var_metadata(                                           &
! String identifier
    'SnowDepth_CABLE',                                                        &
! Variable type
    var_type_cable_snow,                                                      &
! Long name
    "PFT snow depth in layers",                                               &
! Units
    ""                                                                        &
  ) /
!-----------------------------------------------------------------------------
! Metadata for SnowMass_CABLE 
!-----------------------------------------------------------------------------
DATA metadata(559) /  var_metadata(                                           &
! String identifier
    'SnowMass_CABLE',                                                         &
! Variable type
    var_type_cable_snow,                                                      &
! Long name
    "PFT snow mass in layers",                                                &
! Units
    ""                                                                        &
  ) /
!-----------------------------------------------------------------------------
! Metadata for SnowDensity_CABLE 
!-----------------------------------------------------------------------------
DATA metadata(560) /  var_metadata(                                           &
! String identifier
    'SnowDensity_CABLE',                                                      &
! Variable type
    var_type_cable_snow,                                                      &
! Long name
    "PFT snow density in layers",                                             &
! Units
    ""                                                                        &
  ) /
!-----------------------------------------------------------------------------
! Metadata for SnowTemp_CABLE 
!-----------------------------------------------------------------------------
DATA metadata(561) /  var_metadata(                                           &
! String identifier
    'SnowTemp_CABLE',                                                         &
! Variable type
    var_type_cable_snow,                                                      &
! Long name
    "PFT snow temperature in layers",                                         &
! Units
    "m"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for ThreeLayerSnowFlag_CABLE 
!-----------------------------------------------------------------------------
DATA metadata(562) /  var_metadata(                                           &
! String identifier
    'ThreeLayerSnowFlag_CABLE',                                               &
! Variable type
    var_type_cable1l_snow,                                                    &
! Long name
    "3 layer snow flag for tile",                                             &
! Units
    "m"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for SnowAge_CABLE 
!-----------------------------------------------------------------------------
DATA metadata(563) /  var_metadata(                                           &
! String identifier
    'SnowAge_CABLE',                                                          &
! Variable type
    var_type_cable1l_snow,                                                    &
! Long name
    "Snow Age for tile",                                                      &
! Units
    ""                                                                        &
  ) /
!-----------------------------------------------------------------------------
! Metadata for OneLyrSnowDensity_CABLE 
!-----------------------------------------------------------------------------
DATA metadata(564) /  var_metadata(                                           &
! String identifier
    'OneLyrSnowDensity_CABLE',                                                &
! Variable type
    VAR_TYPE_CABLE1l_SNOW,                                                    &
! Long name
    "PFT snow density in aggregate 1 layer",                                  &
! Units
    ""                                                                        &
  ) /
!-----------------------------------------------------------------------------
! Metadata for albsoil_soilt
!-----------------------------------------------------------------------------
DATA metadata(565) / var_metadata(                                            &
! String identifier
    'albsoil_soilt',                                                          &
! Variable type
    var_type_soilt,                                                           &
! Variable is not available for output, so give dummy values for long name
! and units.
! Long name
    "",                                                                       &
! Units
    ""                                                                        &
  ) /
!-----------------------------------------------------------------------------
! Metadata for cs_soilt
!-----------------------------------------------------------------------------
DATA metadata(566) / var_metadata(                                            &
! String identifier
    'cs_soilt',                                                               &
! Variable type
    var_type_soilt_sclayer_scpool,                                            &
! Long name
    "Soil carbon, on soil tiles",                                             &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sthuf_soilt
!-----------------------------------------------------------------------------
DATA metadata(567) / var_metadata(                                            &
! String identifier
    'sthuf_soilt',                                                            &
! Variable type
    var_type_soilt_soil,                                                      &
! Variable is not available for output, so give dummy values for long name
! and units.
! Long name
    "",                                                                       &
! Units
    ""                                                                        &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sthzw_soilt
!-----------------------------------------------------------------------------
DATA metadata(568) / var_metadata(                                            &
! String identifier
    'sthzw_soilt',                                                            &
! Variable type
    var_type_soilt,                                                           &
! Long name
    "Soil wetness in deep LSH/TOPMODEL layer, on soil tiles",                 &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for upward_sw_1
! N.B. Diagnostics upward_sw_1, ..., upward_sw_4 are meaningful only if
! appropriate values of wght_alb are set.
!-----------------------------------------------------------------------------
DATA metadata(569) / var_metadata(                                            &
! String identifier
    'upward_sw_1',                                                            &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Upward SW flux in band 1 (direct beam visible)",                         &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for upward_sw_2
!-----------------------------------------------------------------------------
DATA metadata(570) / var_metadata(                                            &
! String identifier
    'upward_sw_2',                                                            &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Upward SW flux in band 2 (diffuse visible)",                             &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for upward_sw_3
!-----------------------------------------------------------------------------
DATA metadata(571) / var_metadata(                                            &
! String identifier
    'upward_sw_3',                                                            &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Upward SW flux in band 3 (direct beam NIR)",                             &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for upward_sw_4
!-----------------------------------------------------------------------------
DATA metadata(572) / var_metadata(                                            &
! String identifier
    'upward_sw_4',                                                            &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Upward SW flux in band 4 (diffuse NIR)",                                 &
! Units
    "W m-2"                                                                   &
  ) /
!-----------------------------------------------------------------------------
! Metadata for demand_domestic
!-----------------------------------------------------------------------------
DATA metadata(573) / var_metadata(                                            &
! String identifier
    'demand_domestic',                                                        &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Demand for water for domestic use",                                      &
! Units
    "kg"                                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for demand_environment
!-----------------------------------------------------------------------------
DATA metadata(574) / var_metadata(                                            &
! String identifier
    'demand_environment',                                                     &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Demand for water for environmental use",                                 &
! Units
    "kg"                                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for demand_industry
!-----------------------------------------------------------------------------
DATA metadata(575) / var_metadata(                                            &
! String identifier
    'demand_industry',                                                        &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Demand for water for industrial use",                                    &
! Units
    "kg"                                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for demand_irrigation
!-----------------------------------------------------------------------------
DATA metadata(576) / var_metadata(                                            &
! String identifier
    'demand_irrigation',                                                      &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Demand for water for irrigation use",                                    &
! Units
    "kg"                                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for demand_livestock
!-----------------------------------------------------------------------------
DATA metadata(577) / var_metadata(                                            &
! String identifier
    'demand_livestock',                                                       &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Demand for water for livestock use",                                     &
! Units
    "kg"                                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for demand_transfers
!-----------------------------------------------------------------------------
DATA metadata(578) / var_metadata(                                            &
! String identifier
    'demand_transfers',                                                       &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Demand for water for transfers",                                         &
! Units
    "kg"                                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for demand_rate_domestic
!-----------------------------------------------------------------------------
DATA metadata(579) / var_metadata(                                            &
! String identifier
    'demand_rate_domestic',                                                   &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Demand rate for water for domestic use",                                 &
! Units
    "kg s-1"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for demand_rate_industry
!-----------------------------------------------------------------------------
DATA metadata(580) / var_metadata(                                            &
! String identifier
    'demand_rate_industry',                                                   &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Demand rate for water for industrial use",                               &
! Units
    "kg s-1"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for demand_rate_livestock
!-----------------------------------------------------------------------------
DATA metadata(581) / var_metadata(                                            &
! String identifier
    'demand_rate_livestock',                                                  &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Demand rate for water for livestock use",                                &
! Units
    "kg s-1"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for demand_rate_transfers
!-----------------------------------------------------------------------------
DATA metadata(582) / var_metadata(                                            &
! String identifier
    'demand_rate_transfers',                                                  &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Demand rate for water for transfers",                                    &
! Units
    "kg s-1"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for water_demand
!-----------------------------------------------------------------------------
DATA metadata(583) / var_metadata(                                            &
! String identifier
    'water_demand',                                                           &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Demand for water across all sectors",                                    &
! Units
    "kg"                                                                      &
  ) /
!-----------------------------------------------------------------------------
! Metadata for conveyance_loss
!-----------------------------------------------------------------------------
DATA metadata(584) / var_metadata(                                            &
! String identifier
    'conveyance_loss',                                                        &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Conveyancy loss: fraction of abstracted water lost between source " //   &
    "and user",                                                               &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for irrig_eff
!-----------------------------------------------------------------------------
DATA metadata(585) / var_metadata(                                            &
! String identifier
    'irrig_eff',                                                              &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Irrigation efficiency",                                                  &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for sfc_water_frac
!-----------------------------------------------------------------------------
DATA metadata(586) / var_metadata(                                            &
! String identifier
    'sfc_water_frac',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Fraction of demand to be met from surface water",                        &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lake_depth
!-----------------------------------------------------------------------------
DATA metadata(587) / var_metadata(                                            &
! String identifier
    'lake_depth',                                                             &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Lake depth",                                                             &
! Units
    "m"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lake_fetch_gb
!-----------------------------------------------------------------------------
DATA metadata(588) / var_metadata(                                            &
! String identifier
    'lake_fetch_gb',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Typical wind fetch",                                                     &
! Units
    "m"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lake_t_mean_gb
!-----------------------------------------------------------------------------
DATA metadata(589) / var_metadata(                                            &
! String identifier
    'lake_t_mean_gb',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Lake mean temperature",                                                  &
! Units
    "K"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lake_t_mxl_gb
!-----------------------------------------------------------------------------
DATA metadata(590) / var_metadata(                                            &
! String identifier
    'lake_t_mxl_gb',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Lake mixed-layer temperature",                                           &
! Units
    "K"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lake_t_ice_gb
!-----------------------------------------------------------------------------
DATA metadata(591) / var_metadata(                                            &
! String identifier
    'lake_t_ice_gb',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Temperature at upper boundary of lake ice",                              &
! Units
    "K"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lake_h_mxl_gb
!-----------------------------------------------------------------------------
DATA metadata(592) / var_metadata(                                            &
! String identifier
    'lake_h_mxl_gb',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "lake mixed-layer thickness",                                             &
! Units
    "m"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lake_h_ice_gb
!-----------------------------------------------------------------------------
DATA metadata(593) / var_metadata(                                            &
! String identifier
    'lake_h_ice_gb',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Lake ice thickness",                                                     &
! Units
    "m"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lake_shape_factor_gb
!-----------------------------------------------------------------------------
DATA metadata(594) / var_metadata(                                            &
! String identifier
    'lake_shape_factor_gb'   ,                                                &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Thermocline shape factor",                                               &
! Units
    ""                                                                        &
  ) /
!-----------------------------------------------------------------------------
! Metadata for g_dt_gb
!-----------------------------------------------------------------------------
DATA metadata(595) / var_metadata(                                            &
! String identifier
    'g_dt_gb',                                                                &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Ground heat flux over delta T for lakes",                                &
! Units
    "W m-2 K-1"                                                               &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lake_t_sfc_gb
!-----------------------------------------------------------------------------
DATA metadata(596) / var_metadata(                                            &
! String identifier
    'lake_t_sfc_gb',                                                          &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Surface temperature of the lake",                                        &
! Units
    "K"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lake_t_snow_gb
!-----------------------------------------------------------------------------
DATA metadata(597) / var_metadata(                                            &
! String identifier
    'lake_t_snow_gb',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Temperature of the air-snow interface",                                  &
! Units
    "K"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lake_h_snow_gb
!-----------------------------------------------------------------------------
DATA metadata(598) / var_metadata(                                            &
! String identifier
    'lake_h_snow_gb',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Snow thickness used in FLake calculations",                              &
! Units
    "m"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for lake_albedo_gb
!-----------------------------------------------------------------------------
DATA metadata(599) / var_metadata(                                            &
! String identifier
    'lake_albedo_gb',                                                         &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Lake albedo",                                                            &
! Units
    ""                                                                        &
  ) /
!-----------------------------------------------------------------------------
! Metadata for substr_ch4.
!-----------------------------------------------------------------------------
DATA metadata(600) / var_metadata(                                            &
! String identifier
    'substr_ch4',                                                             &
! Variable type
    var_type_ch4layer,                                                        &
! Long name
    "Dissolved substrate in wetland fraction",                                &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for mic_ch4.
!-----------------------------------------------------------------------------
DATA metadata(601) / var_metadata(                                            &
! String identifier
    'mic_ch4',                                                                &
! Variable type
    var_type_ch4layer,                                                        &
! Long name
    "Methanogenic biomass in wetland fraction",                               &
! Units
    "kg m-2"                                                                  &
  ) /
!-----------------------------------------------------------------------------
! Metadata for mic_act_ch4.
!-----------------------------------------------------------------------------
DATA metadata(602) / var_metadata(                                            &
! String identifier
    'mic_act_ch4',                                                            &
! Variable type
    var_type_ch4layer,                                                        &
! Long name
    "Methanogenic biomass activity level in wetland fraction",                &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for acclim_ch4.
!-----------------------------------------------------------------------------
DATA metadata(603) / var_metadata(                                            &
! String identifier
    'acclim_ch4',                                                             &
! Variable type
    var_type_ch4layer,                                                        &
! Long name
    "Acclimation modifier for methanogen calculations",                       &
! Units
    "1"                                                                       &
  ) /
!-----------------------------------------------------------------------------
! Metadata for grid_area
!-----------------------------------------------------------------------------
DATA metadata(604) / var_metadata(                                            &
! String identifier
    'grid_area',                                                              &
! Variable type
    var_type_surface,                                                         &
! Long name
    "Area of each gridbox",                                                   &
! Units
    "m2"                                                                      &
  ) /


!-----------------------------------------------------------------------------
! Visibility declarations
!-----------------------------------------------------------------------------
PRIVATE
PUBLIC                                                                        &
! Parameters
    identifier_len, bl_level_dim_name_out,                                    &
    pft_dim_name_out, cpft_dim_name_out, nvg_dim_name_out, type_dim_name_out, &
    tile_dim_name_out, soilt_dim_name_out, snow_dim_name_out,                 &
    soil_dim_name_out, scpool_dim_name_out, soil_n_pool_dim_name_out,         &
    bedrock_dim_name_out, sclayer_dim_name_out, tracer_dim_name_out,          &
    ch4layer_dim_name_out,                                                    &
! Variables
    bl_level_dim_name,                                                        &
    pft_dim_name, cpft_dim_name, nvg_dim_name, type_dim_name, tile_dim_name,  &
    soilt_dim_name, snow_dim_name, soil_dim_name, scpool_dim_name,            &
    soil_n_pool_dim_name, bedrock_dim_name, sclayer_dim_name,                 &
    tracer_dim_name, bl_level_dim_size,                                       &
    pft_dim_size, cpft_dim_size, nvg_dim_size,                                &
    type_dim_size, tile_dim_size, soilt_dim_size, snow_dim_size,              &
    soil_dim_size, scpool_dim_size, soil_n_pool_dim_size, bedrock_dim_size,   &
    sclayer_dim_size, tracer_dim_size, ch4layer_dim_size,                     &
    cable_tile_dim_size, cable_snow_dim_size, cable_soil_dim_size,            &
! Routines for changing between string and integer identifiers
    get_var_id, get_string_identifier,                                        &
! Variable inquiry routines
    get_var_levs_dims, get_var_attrs,                                         &
! Routines for getting and setting values
    extract_var, populate_var,                                                &
! Routine to check metadata
    check_variable_metadata


CONTAINS


! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.

SUBROUTINE check_variable_metadata

! Check that the metadata for input/output variables do not include
! duplicate names.

IMPLICIT NONE

! Local scalar parameters.
CHARACTER(LEN=*), PARAMETER ::                                                &
  RoutineName = 'CHECK_VARIABLE_METADATA'

! Local scalar variables.
INTEGER :: i  !  Loop counter.

!-----------------------------------------------------------------------------
! It is possible that duplicate names have been coded for variables.
! Duplicates that are coded for input or output will be detected during
! compilation because the SELECT statements in populate_var and extract_var
! will not work. But it is possible that one variable can be included in the
! input code and the other in the output code (or, unusually, included in
! neither!) in which case the duplication will not be detected, with
! unexpected results during run time.
! Here we guard against that possibility by checking the metadata for
! duplicates.
! It is possible that this code will prevent an acceptable run (if that run
! does not need to access a duplicate variable), but in general it is better
! to test and prevent duplicate names.
!-----------------------------------------------------------------------------
DO i = 1, n_vars

  ! Check if any later variable has the same name.
  IF ( ANY( metadata(i+1:n_vars)%identifier == metadata(i)%identifier ) )     &
    THEN

    CALL log_fatal(RoutineName,                                               &
                   "Duplicate name=" // TRIM(metadata(i)%identifier ) //      &
                   ". Code needs to be changed!")
  END IF

END DO

RETURN

END SUBROUTINE check_variable_metadata
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION get_var_id(identifier) RESULT(var_id)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Given a string identifier for a model variable, returns the integer id for
!   that variable
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
CHARACTER(LEN=*), INTENT(IN) :: identifier  ! Identifies the model variable

! Return type
INTEGER :: var_id  ! The integer id for the variable


! Work variables
INTEGER :: i  ! Loop index


!-----------------------------------------------------------------------------


! The returned variable id is the index in the metadata array, so we just
! have to search for the given identifier in the metadata entries
var_id = -1

DO i = 1,n_vars
  IF ( metadata(i)%identifier == identifier ) THEN
    var_id = i
    RETURN
  END IF
END DO

! If we failed to find the identifier, error out
IF ( var_id < 1 )                                                             &
  CALL log_fatal("get_var_id",                                                &
                 "Unrecognised variable identifier - '" //                    &
                 TRIM(identifier) // "'. See docs for available variables")

RETURN

END FUNCTION get_var_id
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION get_string_identifier(var_id) RESULT(identifier)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Given an integer id for a model variable, returns the string identifier
!   for that variable
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
INTEGER, INTENT(IN) :: var_id  ! The integer id for the variable

! Return type
CHARACTER(LEN=identifier_len) :: identifier  ! The string identifier for
                                             ! the model variable


!-----------------------------------------------------------------------------


! var_id is an index in the metadata array
identifier = metadata(var_id)%identifier

RETURN

END FUNCTION get_string_identifier
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

SUBROUTINE get_var_levs_dims(var_id, ndims, dim_names_in, dim_names_out, dim_sizes)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Given an identifer for a model variable, returns the levels dimensions
!   that it uses
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
INTEGER, INTENT(IN) :: var_id  ! Identifies the model variable
INTEGER, INTENT(OUT), OPTIONAL :: ndims
                               ! The number of levels dimensions that the
                               ! variable has
CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: dim_names_in(:)
                               ! The dimension names to use for the variable
                               ! in input files
CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: dim_names_out(:)
                               ! The dimension names to use for the variable
                               ! in output files
INTEGER, INTENT(OUT), OPTIONAL :: dim_sizes(:)
                               ! The dimension sizes to for the variable


!-----------------------------------------------------------------------------


! var_id is an index into the metadata array

! We decide what to return based on the 'variable type' of the requested variable
SELECT CASE ( metadata(var_id)%var_type )
CASE ( var_type_surface )
  ! Just indicate that there are no vertical dimensions
  IF ( PRESENT(ndims) ) ndims = 0

  ! Extra dimension cases for different sizes
CASE ( var_type_bl_level )
  IF ( PRESENT(ndims) ) ndims = 1
  IF ( PRESENT(dim_names_in) ) dim_names_in(1) = bl_level_dim_name
  IF ( PRESENT(dim_names_out) ) dim_names_out(1) = bl_level_dim_name_out
  IF ( PRESENT(dim_sizes) ) dim_sizes(1) = bl_level_dim_size

CASE ( var_type_pft )
  IF ( PRESENT(ndims) ) ndims = 1
  IF ( PRESENT(dim_names_in) ) dim_names_in(1) = pft_dim_name
  IF ( PRESENT(dim_names_out) ) dim_names_out(1) = pft_dim_name_out
  IF ( PRESENT(dim_sizes) ) dim_sizes(1) = pft_dim_size

CASE ( var_type_cpft )
  IF ( PRESENT(ndims) ) ndims = 1
  IF ( PRESENT(dim_names_in) ) dim_names_in(1) = cpft_dim_name
  IF ( PRESENT(dim_names_out) ) dim_names_out(1) = cpft_dim_name_out
  IF ( PRESENT(dim_sizes) ) dim_sizes(1) = cpft_dim_size

CASE ( var_type_nvg )
  IF ( PRESENT(ndims) ) ndims = 1
  IF ( PRESENT(dim_names_in) ) dim_names_in(1) = nvg_dim_name
  IF ( PRESENT(dim_names_out) ) dim_names_out(1) = nvg_dim_name_out
  IF ( PRESENT(dim_sizes) ) dim_sizes(1) = nvg_dim_size

CASE ( var_type_type )
  IF ( PRESENT(ndims) ) ndims = 1
  IF ( PRESENT(dim_names_in) ) dim_names_in(1) = type_dim_name
  IF ( PRESENT(dim_names_out) ) dim_names_out(1) = type_dim_name_out
  IF ( PRESENT(dim_sizes) ) dim_sizes(1) = type_dim_size

CASE ( var_type_surft )
  IF ( PRESENT(ndims) ) ndims = 1
  IF ( PRESENT(dim_names_in) ) dim_names_in(1) = tile_dim_name
  IF ( PRESENT(dim_names_out) ) dim_names_out(1) = tile_dim_name_out
  IF ( PRESENT(dim_sizes) ) dim_sizes(1) = tile_dim_size

CASE ( var_type_soilt )
  IF ( PRESENT(ndims) ) ndims = 1
  IF ( PRESENT(dim_names_in) ) dim_names_in(1) = soilt_dim_name
  IF ( PRESENT(dim_names_out) ) dim_names_out(1) = soilt_dim_name_out
  IF ( PRESENT(dim_sizes) ) dim_sizes(1) = soilt_dim_size

CASE ( var_type_soil )
  IF ( PRESENT(ndims) ) ndims = 1
  IF ( PRESENT(dim_names_in) ) dim_names_in(1) = soil_dim_name
  IF ( PRESENT(dim_names_out) ) dim_names_out(1) = soil_dim_name_out
  IF ( PRESENT(dim_sizes) ) dim_sizes(1) = soil_dim_size

CASE ( var_type_scpool )
  IF ( PRESENT(ndims) ) ndims = 2
  IF ( PRESENT(dim_names_in) )                                                &
    dim_names_in(1:2) = (/ sclayer_dim_name, scpool_dim_name /)
  IF ( PRESENT(dim_names_out) )                                               &
    dim_names_out(1:2) = (/ sclayer_dim_name_out, scpool_dim_name_out /)
  IF ( PRESENT(dim_sizes) )                                                   &
    dim_sizes(1:2) = (/ sclayer_dim_size, scpool_dim_size /)

CASE ( var_type_soil_n_pool )
  ! Variables with soil (C) layers and soil N pools.
  IF ( PRESENT(ndims) ) ndims = 2
  IF ( PRESENT(dim_names_in) )                                                &
    dim_names_in(1:2) = (/ sclayer_dim_name, soil_n_pool_dim_name /)
  IF ( PRESENT(dim_names_out) )                                               &
    dim_names_out(1:2) = (/ sclayer_dim_name_out, soil_n_pool_dim_name_out /)
  IF ( PRESENT(dim_sizes) )                                                   &
    dim_sizes(1:2) = (/ sclayer_dim_size, soil_n_pool_dim_size /)

CASE ( var_type_sclayer )
  IF ( PRESENT(ndims) ) ndims = 1
  IF ( PRESENT(dim_names_in) ) dim_names_in(1) = sclayer_dim_name
  IF ( PRESENT(dim_names_out) ) dim_names_out(1) = sclayer_dim_name_out
  IF ( PRESENT(dim_sizes) ) dim_sizes(1) = sclayer_dim_size

CASE ( var_type_ch4layer )
  IF ( PRESENT(ndims) ) ndims = 1
  IF ( PRESENT(dim_names_in) ) dim_names_in(1) = ch4layer_dim_name
  IF ( PRESENT(dim_names_out) ) dim_names_out(1) = ch4layer_dim_name_out
  IF ( PRESENT(dim_sizes) ) dim_sizes(1) = ch4layer_dim_size

CASE ( var_type_tracer )
  IF ( PRESENT(ndims) )         ndims            = 1
  IF ( PRESENT(dim_names_in) )  dim_names_in(1)  = tracer_dim_name
  IF ( PRESENT(dim_names_out) ) dim_names_out(1) = tracer_dim_name_out
  IF ( PRESENT(dim_sizes) )     dim_sizes(1)     = tracer_dim_size

CASE ( var_type_bedrock )
  IF ( PRESENT(ndims) ) ndims = 1
  IF ( PRESENT(dim_names_in) ) dim_names_in(1) = bedrock_dim_name
  IF ( PRESENT(dim_names_out) ) dim_names_out(1) = bedrock_dim_name_out
  IF ( PRESENT(dim_sizes) ) dim_sizes(1) = bedrock_dim_size

CASE ( var_type_snow )
  ! Snow variables have 2 levels dimensions
  IF ( PRESENT(ndims) ) ndims = 2
  IF ( PRESENT(dim_names_in) )                                                &
    dim_names_in(1:2) = (/ tile_dim_name, snow_dim_name /)
  IF ( PRESENT(dim_names_out) )                                               &
    dim_names_out(1:2) = (/ tile_dim_name_out, snow_dim_name_out /)
  IF ( PRESENT(dim_sizes) )                                                   &
    dim_sizes(1:2) = (/ tile_dim_size, snow_dim_size /)

CASE ( var_type_soilt_soil )
  ! Variables with soil tiles and levels
  IF ( PRESENT(ndims) ) ndims = 2
  IF ( PRESENT(dim_names_in) )                                                &
    dim_names_in(1:2) = (/ soilt_dim_name, soil_dim_name /)
  IF ( PRESENT(dim_names_out) )                                               &
    dim_names_out(1:2) = (/ soilt_dim_name_out, soil_dim_name_out /)
  IF ( PRESENT(dim_sizes) )                                                   &
    dim_sizes(1:2) = (/ soilt_dim_size, soil_dim_size /)

CASE ( var_type_soilt_sclayer )
  ! Variables with soil tiles and soil C levels
  IF ( PRESENT(ndims) ) ndims = 2
  IF ( PRESENT(dim_names_in) )                                                &
    dim_names_in(1:2) = (/ soilt_dim_name, sclayer_dim_name /)
  IF ( PRESENT(dim_names_out) )                                               &
    dim_names_out(1:2) = (/ soilt_dim_name_out, sclayer_dim_name_out /)
  IF ( PRESENT(dim_sizes) )                                                   &
    dim_sizes(1:2) = (/ soilt_dim_size, sclayer_dim_size /)

CASE ( var_type_soilt_sclayer_scpool )
  ! Variables with soil tiles, soil C levels, and soil C pools.
  IF ( PRESENT(ndims) ) ndims = 3
  IF ( PRESENT(dim_names_in) )                                                &
    dim_names_in(1:3) = (/ soilt_dim_name, sclayer_dim_name,                  &
                           scpool_dim_name /)
  IF ( PRESENT(dim_names_out) )                                               &
    dim_names_out(1:3) = (/ soilt_dim_name_out, sclayer_dim_name_out,         &
                            scpool_dim_name_out  /)
  IF ( PRESENT(dim_sizes) )                                                   &
    dim_sizes(1:3) = (/ soilt_dim_size, sclayer_dim_size, scpool_dim_size /)

CASE ( var_type_pft_sclayer )
  IF ( PRESENT(ndims) ) ndims = 2
  IF ( PRESENT(dim_names_in) )                                                &
    dim_names_in(1:2) = (/ pft_dim_name, sclayer_dim_name /)
  IF ( PRESENT(dim_names_out) )                                               &
    dim_names_out(1:2) = (/ pft_dim_name_out, sclayer_dim_name_out /)
  IF ( PRESENT(dim_sizes) )                                                   &
    dim_sizes(1:2) = (/ pft_dim_size, sclayer_dim_size /)

  ! These CABLE var types can piggy back on existing types if so ordained
CASE ( var_type_cable_soil )
  ! CABLE soil variables have 2 levels dimensions
  IF ( PRESENT(ndims) ) ndims = 2
  IF ( PRESENT(dim_names_in) )                                                &
    dim_names_in(1:2) = (/ tile_dim_name, soil_dim_name /)
  IF ( PRESENT(dim_names_out) )                                               &
    dim_names_out(1:2) = (/ tile_dim_name_out, soil_dim_name_out /)
  IF ( PRESENT(dim_sizes) )                                                   &
    dim_sizes(1:2) = (/ cable_tile_dim_size, cable_soil_dim_size /)
 
CASE ( var_type_cable_snow )
  ! CABLE snow variables have 2 levels dimensions
  IF ( PRESENT(ndims) ) ndims = 2
  IF ( PRESENT(dim_names_in) )                                                &
    dim_names_in(1:2) = (/ tile_dim_name, snow_dim_name /)
  IF ( PRESENT(dim_names_out) )                                               &
    dim_names_out(1:2) = (/ tile_dim_name_out, snow_dim_name_out /)
  IF ( PRESENT(dim_sizes) )                                                   &
    dim_sizes(1:2) = (/ cable_tile_dim_size, cable_snow_dim_size /)
 
CASE ( var_type_cable1l_snow )
  ! CABLE snow variables per tile only have 1 level dimensions
  IF ( PRESENT(ndims) ) ndims = 1
  IF ( PRESENT(dim_names_in) )                                                &
    dim_names_in = tile_dim_name
  IF ( PRESENT(dim_names_out) )                                               &
    dim_names_out = tile_dim_name_out
  IF ( PRESENT(dim_sizes) )                                                   &
    dim_sizes = cable_tile_dim_size

CASE DEFAULT
  ! Unrecognised variable type
  CALL log_fatal("get_var_levs_dims",                                         &
                 "Unrecognised variable type - check metadata")
END SELECT

RETURN

END SUBROUTINE get_var_levs_dims
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE get_var_attrs(var_id, int_attributes, real_attributes,             &
                         char_attributes)

USE io_constants, ONLY: max_attr_var, mdi

USE dictionary_mod, ONLY: dict, dict_create, dict_set

USE string_utils_mod, ONLY: to_string

USE jules_soil_biogeochem_mod, ONLY:                                          &
!  imported scalar parameters
     soil_model_ecosse, soil_model_rothc, soil_model_1pool,                   &
!  imported scalars (IN)
     soil_bgc_model

USE jules_soil_mod, ONLY: zsmc

USE jules_vegetation_mod, ONLY: n_day_photo_acclim

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Given an identifer for a model variable, returns the attributes that
!   should be defined on output files
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
INTEGER, INTENT(IN) :: var_id  ! Identifies the model variable
TYPE(dict), INTENT(OUT) :: int_attributes   ! Attributes with integer values
TYPE(dict), INTENT(OUT) :: real_attributes  ! Attributes with real values
TYPE(dict), INTENT(OUT) :: char_attributes  ! Attributes with character
                                            ! values


! Work variables
CHARACTER(LEN=identifier_len) :: identifier  ! The string identifier for the
                                             ! variable


!-----------------------------------------------------------------------------


! Create the dictionary objects
int_attributes = dict_create(max_attr_var, INT(1))
real_attributes = dict_create(max_attr_var, REAL(1.0))
char_attributes = dict_create(max_attr_var, "char")


! Get the string identifier
identifier = get_string_identifier(var_id)


! All variables except lat and lon will have a missing data attribute and
! a coordinate attribute.
! lat and lon do not need these since they are interpreted as coordinate
! variables as per the CF convention.
IF ( identifier /= 'latitude' .AND. identifier /= 'longitude' ) THEN
  CALL dict_set(real_attributes, "missing_value", mdi)
  CALL dict_set(real_attributes, "_FillValue", mdi)
  CALL dict_set(char_attributes, "coordinates", "latitude longitude")
END IF


! Add the standard "long_name" and "units" attributes from the metadata array
! var_id is an index into the metadata array
CALL dict_set(char_attributes, "long_name", metadata(var_id)%long_name)
CALL dict_set(char_attributes, "units", metadata(var_id)%units)


! Use a select statement to add extra attributes/override current values if
! required. The SELECT statement uses the string identifier to select a CASE
! to avoid being dependent on the implementation of integer ids.
SELECT CASE ( identifier )
CASE ( 'latitude' )
  CALL dict_set(char_attributes, "standard_name", "latitude")

CASE ( 'longitude' )
  CALL dict_set(char_attributes, "standard_name", "longitude")

CASE ( 'frac' )
  CALL dict_set(real_attributes, "valid_min", 0.0)
  CALL dict_set(real_attributes, "valid_max", 1.0)

CASE ( 'smc_avail_top' )
  ! Overwrite long_name to use a value that is calculated using zsmc.
  CALL dict_set(char_attributes, "long_name",                                 &
                "Gridbox available moisture in top " //                       &
                TRIM(to_string(zsmc)) // "m of soil")

CASE ( 'cs' )
  ! Overwrite long_name to use a value that specifies the soil C pools.
  SELECT CASE ( soil_bgc_model )
  CASE ( soil_model_ecosse, soil_model_rothc )
    CALL dict_set(char_attributes, "long_name",                               &
                  "Gridbox soil carbon in each pool (DPM,RPM,bio,hum)")
  CASE ( soil_model_1pool) 
    CALL dict_set(char_attributes, "long_name",                               &
                  "Gridbox soil carbon (single pool)")
  END SELECT

CASE ( 'resp_s' )
  ! Overwrite long_name to use a value that specifies the soil C pools.
  SELECT CASE ( soil_bgc_model )
  CASE ( soil_model_ecosse, soil_model_rothc )
    CALL dict_set(char_attributes, "long_name",                               &
      "Gridbox soil respiration from each pool (DPM,RPM,bio,hum)")
  CASE ( soil_model_1pool) 
    CALL dict_set(char_attributes, "long_name",                               &
      "Gridbox soil respiration (single pool)")
  END SELECT

CASE ( 'temp_ave_nday' )
  ! Overwrite long_name to include length of average.
  CALL dict_set(char_attributes, "long_name",                                 &
                "Average temperature (exp filter) over previous " //          &
                TRIM(to_string(n_day_photo_acclim)) // " days")

END SELECT

RETURN

END SUBROUTINE get_var_attrs
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

FUNCTION extract_var(var_id) RESULT(cube)

USE jules_fields_mod, ONLY: crop_vars, psparms, toppdm, fire_vars, ainfo,     &
                            trif_vars, soilecosse, progs, trifctltype,        &
                            jules_vars

  !Use subroutines
! TYPE Definitions
USE jules_fields_mod, ONLY: toppdm

USE gridbox_mean_mod, ONLY:                                                   &
  surftiles_to_gbm, soiltiles_to_gbm

!Use variables
!Science variables
USE water_constants_mod, ONLY:                                                &
  rho_water, tm

USE csigma, ONLY:                                                             &
  sbcon

USE jules_surface_types_mod, ONLY:                                            &
  npft, lake, ncpft

USE theta_field_sizes, ONLY:                                                  &
  t_i_length, t_j_length

USE jules_surface_mod, ONLY:                                                  &
  l_aggregate, i_aggregate_opt

USE jules_water_resources_mod, ONLY:                                          &
  conveyance_loss, demand_accum, irrig_eff, l_water_domestic,                 &
  l_water_environment, l_water_industry, l_water_irrigation,                  &
  l_water_livestock, l_water_transfers, sfc_water_frac, use_domestic,         &
  use_environment, use_industry, use_irrigation, use_livestock, use_transfers

USE jules_radiation_mod, ONLY:                                                &
  wght_alb

USE jules_vegetation_mod, ONLY:                                               &
  l_vegdrag_surft

USE c_z0h_z0m, ONLY:                                                          &
  z0h_z0m

USE model_grid_mod, ONLY:                                                     &
  grid_area_ij, latitude, longitude

USE ancil_info, ONLY:                                                         &
  dim_cs1, dim_cslayer, land_pts, nsurft, soil_pts, surft_pts, nsoilt

USE soil_ecosse_vars_mod, ONLY:                                               &
  i_amm, i_nit

USE fluxes, ONLY:                                                             &
  alb_surft, fqw_surft, hf_snow_melt_gb, sub_surf_roff_gb,                    &
  surf_roff_gb, snomlt_sub_htf_gb, snow_melt_gb,                              &
  tot_tfall_gb, ecan_ij, ei_ij, esoil_ij_soilt,                               &
  land_albedo_ij, surf_ht_flux_ij,                                            &
  melt_surft, anthrop_heat_surft, emis_surft,                                 &
  ext_soilt, fsmc_pft, ftl_surft,le_surft, radnet_surft,                      &
  surf_ht_store_surft, surf_htf_surft, ecan_surft,                            &
  ei_surft, esoil_surft, rflow_gb, rrun_gb,                                   &
  snow_soil_htf, sw_surft, t_growth_gb, z0m_surft, z0h_surft

USE gridmean_fluxes, ONLY:                                                    &
  fqw_1_ij,ftl_1_ij,taux_1_ij,tauy_1_ij

USE veg_param, ONLY:                                                          &
  secs_per_360days

USE jules_deposition_mod, ONLY:                                               &
  tracer_field

USE jules_soil_biogeochem_mod, ONLY: soil_bgc_model, soil_model_ecosse,       &
  soil_model_rothc

USE jules_soil_ecosse_mod, ONLY: dt_soilc

USE jules_soil_mod, ONLY:                                                     &
  dzsoil, zsmc, sm_levels, ns_deep

USE forcing, ONLY:                                                            &
  con_rain_ij, con_snow_ij, ls_rain_ij, ls_snow_ij,                           &
  lw_down_ij, pstar_ij, qw_1_ij, sw_down_ij, tl_1_ij,                         &
  u_1_ij, v_1_ij, diurnal_temperature_range_ij

USE sf_diags_mod, ONLY:                                                       &
  sf_diag

USE jules_snow_mod, ONLY:                                                     &
  canSnowTile

USE overbank_inundation_mod, ONLY:                                            &
  frac_fplain_lp

USE ozone_vars, ONLY:                                                         &
  flux_o3_pft, fo3_pft

USE bvoc_vars,                ONLY:                                           &
  isoprene_gb, isoprene_pft, terpene_gb , terpene_pft,                        &
  methanol_gb, methanol_pft, acetone_gb, acetone_pft

USE aero, ONLY:                                                               &
  co2_mmr

USE fire_mod,  ONLY:                                                          &
  fire_prog, fire_diag

USE imogen_progs, ONLY:                                                       &
  d_land_atmos_co2, d_ocean_atmos, c_emiss_out

USE lake_mod,   ONLY:                                                         &
  lake_depth_gb, lake_fetch_gb, lake_t_mean_gb, lake_t_mxl_gb, lake_t_ice_gb, &
  lake_h_mxl_gb, lake_h_ice_gb, lake_shape_factor_gb, g_dt_gb, lake_t_sfc_gb, &
  lake_t_snow_gb, lake_h_snow_gb, lake_albedo_gb

!Others
USE io_constants, ONLY:                                                       &
  mdi
USE logging_mod, ONLY:                                                        &
  log_fatal
USE data_cube_mod, ONLY:                                                      &
  data_cube, cube_from_array,                                                 &
  cube_free, operator (*), operator (-)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Given an identifer for a model variable, returns the data currently
!   associated with that variable as a cube
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
INTEGER, INTENT(IN) :: var_id
                     ! Identifies the variable to extract data from


! Return
TYPE(data_cube) :: cube  ! The extracted data


! Work variables
TYPE(data_cube) :: cube_land  ! Workspace cube for land data
                              ! This is required so that it can be deallocated
                              ! to avoid memory leaks

TYPE(data_cube) :: upward_cube, emis_cube, downward_cube, lw_down_cube
                              ! Work cubes for use in calculation of lw_net
                              ! Required to avoid memory leaks

REAL :: workspace_land(land_pts)  ! Used as a space for calculations
                                  ! before creating a cube

REAL :: workspace_surft(land_pts, nsurft)
                                ! Used in calculation of tile variables

REAL :: workspace_levs(land_pts, sm_levels)
                                ! Used in calculation of gridbox mean (GBM)
                                ! on soil levels

REAL :: workspace_cpft(land_pts, ncpft)
                                ! Used in calculation of crop tile variables

REAL :: workspace_cs(land_pts, dim_cslayer)
                                ! Used in calculation of soil C layer
                                ! variables.

REAL :: ones(land_pts, sm_levels)  ! An array full of ones to aid with
                                   ! calculation of first frozen/unfrozen
                                   ! layer
INTEGER :: layer(land_pts)  ! Used in calculation of first frozen/unfrozen
                            ! soil layer

REAL :: sum_frac(land_pts)  ! Used in calculation of fsmc_gb - the sum
                            ! of frac_surft over all pft tiles

REAL :: dz  ! Used in calculation of available soil moisture
            ! Size of current soil layer
REAL :: ztop  ! Used in calculation of available soil moisture
              ! Depth to top of current soil layer

INTEGER :: i,j,k,l,n,m  ! Index variables

REAL :: co2_tmp(land_pts) ! Used to out put the 1-d co2_mmr variable

!-----------------------------------------------------------------------------

! Initialise data to missing data value
workspace_land(:)     = mdi
workspace_surft(:,:)  = mdi
workspace_cpft(:,:)   = mdi
workspace_levs(:,:)   = mdi
workspace_cs(:,:)     = mdi
ones(:,:)             = 1.0

! We use the string identifier to search for CASEs in the SELECT, rather than
! being sensitive to the implementation of integer variable ids
SELECT CASE ( get_string_identifier(var_id) )
CASE ( 'latitude' )
  cube = cube_from_array(latitude)

CASE ( 'longitude' )
  cube = cube_from_array(longitude)

CASE ( 'frac' )
  cube_land = cube_from_array(ainfo%frac_surft)
  cube = map_from_land(cube_land)

CASE ( 'b' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(psparms%bexp_soilt(:,:,k),ainfo)
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for b with nsoilt>1: " //     &
                   "calculation is unphysical. Use b_soilt instead.")
  END IF

CASE ( 'b_soilt' )
  cube_land = cube_from_array(psparms%bexp_soilt)
  cube = map_from_land(cube_land)

CASE ( 'clay' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, dim_cslayer
      workspace_levs(:,k) = soiltiles_to_gbm(psparms%clay_soilt(:,:,k),ainfo)
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:dim_cslayer))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for clay with nsoilt>1: " //  &
                   "calculation is unphysical. Use psparms%clay_soilt instead.")
  END IF

CASE ( 'clay_soilt' )
  cube_land = cube_from_array(psparms%clay_soilt)
  cube = map_from_land(cube_land)

CASE ( 'sathh' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(psparms%sathh_soilt(:,:,k),ainfo)
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for sathh with nsoilt>1: " // &
                  "calculation is unphysical. Use psparms%sathh_soilt instead.")
  END IF

CASE ( 'sathh_soilt' )
  cube_land = cube_from_array(psparms%sathh_soilt)
  cube = map_from_land(cube_land)

CASE ( 'satcon' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(psparms%satcon_soilt(:,:,k),ainfo)
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for satcon with nsoilt>1: " //&
                 "calculation is unphysical. Use psparms%satcon_soilt instead.")
  END IF

CASE ( 'satcon_soilt' )
  cube_land = cube_from_array(psparms%satcon_soilt)
  cube = map_from_land(cube_land)

CASE ( 'sm_sat' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(psparms%smvcst_soilt(:,:,k),ainfo)
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for sm_sat with nsoilt>1: " //&
                   "calculation is unphysical. Use sm_sat_soilt instead.")
  END IF

CASE ( 'sm_sat_soilt' )
  cube_land = cube_from_array(psparms%smvcst_soilt)
  cube = map_from_land(cube_land)

CASE ( 'sm_crit' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(psparms%smvccl_soilt(:,:,k),ainfo)
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for sm_crit with nsoilt>1:" //&
                   " calculation is unphysical. Use sm_crit_soilt instead.")
  END IF

CASE ( 'sm_crit_soilt' )
  cube_land = cube_from_array(psparms%smvccl_soilt)
  cube = map_from_land(cube_land)

CASE ( 'sm_wilt' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(psparms%smvcwt_soilt(:,:,k),ainfo)
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for sm_wilt with nsoilt>1:" //&
                   " calculation is unphysical. Use sm_wilt_soilt instead.")
  END IF

CASE ( 'sm_wilt_soilt' )
  cube_land = cube_from_array(psparms%smvcwt_soilt)
  cube = map_from_land(cube_land)

CASE ( 'hcap' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(psparms%hcap_soilt(:,:,k),ainfo)
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for hcap with nsoilt>1: " //  &
                   "calculation is unphysical. Use psparms%hcap_soilt instead.")
  END IF

CASE ( 'hcap_soilt' )
  cube_land = cube_from_array(psparms%hcap_soilt)
  cube = map_from_land(cube_land)

CASE ( 'hcon' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(psparms%hcon_soilt(:,:,k),ainfo)
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for hcon with nsoilt>1: " //  &
                   "calculation is unphysical. Use psparms%hcon_soilt instead.")
  END IF

CASE ( 'hcon_soilt' )
  cube_land = cube_from_array(psparms%hcon_soilt)
  cube = map_from_land(cube_land)

  !-----------------------------------------------------------------------------
  ! Soil ancillaries on soil carbon layers.
  !-----------------------------------------------------------------------------
CASE ( 'soil_ph' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, dim_cslayer
      workspace_cs(:,k) = soiltiles_to_gbm(psparms%soil_ph_soilt(:,:,k),ainfo)
    END DO
    cube_land = cube_from_array(workspace_cs(:,1:dim_cslayer))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal( "extract_var",                                            &
                    "Cannot output GB mean value for soil_ph "           //   &
                    "with nsoilt>1: calculation is unphysical. "         //   &
                    "Use psparms%soil_ph_soilt instead." )
  END IF

CASE ( 'soil_ph_soilt' )
  cube_land = cube_from_array(psparms%soil_ph_soilt)
  cube = map_from_land(cube_land)

CASE ( 'albedo_land' )
  ! Calculate the albedo as used in subroutine control when calculating the net
  ! shortwave on tiles
  ! Here we take the average of diffuse albedos in VIS and NIR
  cube_land = cube_from_array(                                                &
    surftiles_to_gbm( (wght_alb(1) * alb_surft(:,:,1) +                       &
                   wght_alb(2) * alb_surft(:,:,2) +                           &
                   wght_alb(3) * alb_surft(:,:,3) +                           &
                   wght_alb(4) * alb_surft(:,:,4)), ainfo )                   &
  )
  cube = map_from_land(cube_land)

CASE ( 'canopy_gb' )
  !     Don't use the canopy_gb variable, as this is calculated part-way through
  !     a timestep.
  cube_land = cube_from_array(surftiles_to_gbm( progs%canopy_surft, ainfo ))
  cube = map_from_land(cube_land)

CASE ( 'cs_gb' )
  IF ( nsoilt == 1) THEN
    !Case for a single soil tile
    cube_land = cube_from_array(SUM(sum(progs%cs_pool_soilt(:,1,:,:), 3), 2))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output cs_gb for soil tiled runs")
  END IF

CASE ( 'cs_soilt' )
  cube_land = cube_from_array(progs%cs_pool_soilt)
  cube = map_from_land(cube_land)

CASE ( 'cv' )
  cube_land = cube_from_array(trifctltype%cv_gb)
  cube = map_from_land(cube_land)

CASE ( 'depth_frozen' )
  IF ( nsoilt == 1) THEN
    ! Get frozen depth from surface
    ! Start by assuming that there are no frozen layers anywhere,
    !and so frozen depth is 0
    workspace_land(:) = 0.0

    ! Get the first unfrozen layer for every land point using MINLOC
    ! We give an array of ones as the array to take mins from,
    ! so that we get the minimum index at which the mask is true
    layer(:) = MINLOC(ones, 2, progs%t_soil_soilt(:,1,:) >= tm)

    ! If the layer found above is 0 for any point, that means
    !no unfrozen layers were found and hence the whole column is frozen
    WHERE ( layer(:) == 0 )
      workspace_land(:) = SUM(dzsoil(:))
    END WHERE

    ! At points where the layer found above is 1, that means
    ! no frozen layers so we can ignore those points
    ! Otherwise, interpolate to estimate depth of zero degC isotherm
    ! Do all points with the same unfrozen layer at once
    DO i = 2,sm_levels
      WHERE ( layer(:) == i )
        workspace_land(:) = SUM(dzSoil(1:i-1)) + dzSoil(i) *                  &
                          (tm - progs%t_soil_soilt(:,1,i-1)) /                &
                          (progs%t_soil_soilt(:,1,i) -                        &
                           progs%t_soil_soilt(:,1,i-1))
      END WHERE
    END DO
    ! Map the result onto the full grid
    cube_land = cube_from_array(workspace_land)
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for depth_frozen: " //        &
                   "soil tiled diagnostic not yet implemented.")
  END IF

CASE ( 'depth_frozen_sthf' )
  IF ( nsoilt == 1) THEN
    ! Get frozen depth from surface in terms on frozen/unfrozen soil moisture
    ! Method similar in construction to depth_frozen
    workspace_land(:) = 0.0
    layer(:) = MINLOC(ones, 2, psparms%sthf_soilt(:,1,:) == 0.0)
    WHERE ( layer(:) == 0 )
      workspace_land(:) = SUM(dzsoil(:)) - dzsoil(sm_levels) *                &
            (1 - psparms%sthf_soilt(:,1,sm_levels) /                          &
            (psparms%sthf_soilt(:,1,sm_levels) +                              &
            psparms%sthu_soilt(:,1,sm_levels) -                               &
            psparms%sthu_min_soilt(:,1,sm_levels)))
      ! The bottom layer can be below freezing temp but only partially frozen.
    END WHERE
    DO i = 2,sm_levels
      WHERE ( layer(:) == i )
        workspace_land(:) = SUM(dzSoil(1:i-1)) - dzsoil(i-1) *                &
                (1 - psparms%sthf_soilt(:,1,i-1) /                            &
                (psparms%sthf_soilt(:,1,i-1) +                                &
                 psparms%sthu_soilt(:,1,i-1) - psparms%sthu_min_soilt(:,1,i-1)))
      END WHERE
    END DO
    ! Map the result onto the full grid
    cube_land = cube_from_array(workspace_land)
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for depth_frozen_sthf: " //   &
                   "soil tiled diagnostic not yet implemented.")
  END IF

CASE ( 'depth_unfrozen' )
  IF ( nsoilt == 1) THEN
    ! Get unfrozen depth from surface
    ! See depth_frozen for method description, but swapping frozen for unfrozen
    ! and vica-versa
    workspace_land(:) = 0.0
    layer(:) = MINLOC(ones, 2, progs%t_soil_soilt(:,1,:) < tm)
    WHERE ( layer(:) == 0 )
      workspace_land(:) = SUM(dzsoil(:))
    END WHERE
    DO i = 2,sm_levels
      WHERE ( layer(:) == i )
        workspace_land(:) = SUM(dzSoil(1:i-1)) + dzSoil(i) *                  &
                          (tm - progs%t_soil_soilt(:,1,i-1)) /                &
                          (progs%t_soil_soilt(:,1,i) -                        &
                           progs%t_soil_soilt(:,1,i-1))
      END WHERE
    END DO
    ! Map the result onto the full grid
    cube_land = cube_from_array(workspace_land)
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for depth_unfrozen: "  //     &
                   "soil tiled diagnostic not yet implemented.")
  END IF

CASE ( 'depth_unfrozen_sthf' )
  IF ( nsoilt == 1) THEN
    ! Get unfrozen depth from surface in terms of frozen/unfrozen
    ! soil moisture
    ! See depth_frozen for method, but swapping frozen for unfrozen
    ! & vica-versa
    workspace_land(:) = 0.0
    layer(:) = MINLOC(ones, 2, psparms%sthf_soilt(:,1,:) > 0.0)
    WHERE ( layer(:) == 0 )
      workspace_land(:) = SUM(dzsoil(:))
    END WHERE
    DO i = 1,sm_levels           ! HERE WE CAN NO LONGER IGNORE i=1.
      WHERE ( layer(:) == i )
        workspace_land(:) = SUM(dzSoil(1:i)) - dzSoil(i) *                    &
                            psparms%sthf_soilt(:,1,i) /                       &
                            (psparms%sthf_soilt(:,1,i) +                      &
                             psparms%sthu_soilt(:,1,i) -                      &
                             psparms%sthu_min_soilt(:,1,i))
      END WHERE
    END DO
    ! Map the result onto the full grid
    cube_land = cube_from_array(workspace_land)
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for depth_unfrozen_sthf: " // &
                   "soil tiled diagnostic not yet implemented.")
  END IF

CASE ( 'drain' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(toppdm%drain_soilt,ainfo))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for drain with nsoilt>1: " // &
                   "calculation is unphysical. Use toppdm%drain_soilt instead.")
  END IF

CASE ( 'drain_soilt' )
  cube_land = cube_from_array(toppdm%drain_soilt)
  cube = map_from_land(cube_land)

CASE ( 'elake' )
  workspace_land(:) = 0.0
  IF ( .NOT. l_aggregate .AND. lake > 0 )                                     &
    workspace_land(:) = fqw_surft(:,lake) * ainfo%frac_surft(:,lake)
  cube_land = cube_from_array(workspace_land)
  cube = map_from_land(cube_land)

CASE ( 'emis_gb' )
  cube_land = cube_from_array(surftiles_to_gbm(emis_surft, ainfo))
  cube = map_from_land(cube_land)

CASE ( 'fch4_wetl' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(toppdm%fch4_wetl_soilt,ainfo))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for fch4_wetl with "  //      &
                   "nsoilt>1. Use toppdm%fch4_wetl_soilt instead.")
  END IF

CASE ( 'fch4_wetl_soilt' )
  cube_land = cube_from_array(toppdm%fch4_wetl_soilt)
  cube = map_from_land(cube_land)

CASE ( 'fch4_wetl_cs' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(toppdm%fch4_wetl_cs_soilt,   &
                                                 ainfo))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for fch4_wetl_cs with " //    &
                   " nsoilt>1. Use toppdm%fch4_wetl_cs_soilt instead.")
  END IF

CASE ( 'fch4_wetl_cs_soilt' )
  cube_land = cube_from_array(toppdm%fch4_wetl_cs_soilt)
  cube = map_from_land(cube_land)

CASE ( 'fch4_wetl_npp' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(toppdm%fch4_wetl_npp_soilt,  &
                                                 ainfo))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for fch4_wetl_npp with " //   &
                   "nsoilt>1. Use toppdm%fch4_wetl_npp_soilt instead.")
  END IF

CASE ( 'fch4_wetl_npp_soilt' )
  cube_land = cube_from_array(toppdm%fch4_wetl_npp_soilt)
  cube = map_from_land(cube_land)

CASE ( 'fch4_wetl_resps' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(toppdm%fch4_wetl_resps_soilt,&
                                                 ainfo))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for fch4_wetl_resps with " // &
                   "nsoilt>1. Use toppdm%fch4_wetl_resps_soilt instead.")
  END IF

CASE ( 'fch4_wetl_resps_soilt' )
  cube_land = cube_from_array(toppdm%fch4_wetl_resps_soilt)
  cube = map_from_land(cube_land)

CASE ( 'fsat' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(toppdm%fsat_soilt,ainfo))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for fsat with nsoilt>1. " //  &
                   "Use toppdm%fsat_soilt instead.")
  END IF

CASE ( 'fsat_soilt' )
  cube_land = cube_from_array(toppdm%fsat_soilt)
  cube = map_from_land(cube_land)

CASE ( 'fsmc_gb' )
  ! Calculate gridbox mean over PFTs.
  ! Calculate the weighted sum over pfts
  workspace_land(:) = SUM(fsmc_pft * ainfo%frac_surft(:,1:npft), 2)
  sum_frac(:) = SUM(ainfo%frac_surft(:,1:npft), 2)
  ! Normalise to the vegetation fraction
  WHERE ( sum_frac > EPSILON(1.0) )
    workspace_land(:) = workspace_land(:) / sum_frac(:)
  ELSEWHERE
    ! Where there is no veg, set an impossible value
    workspace_land(:) = mdi
  END WHERE
  ! Copy the result onto the full grid
  cube_land = cube_from_array(workspace_land)
  cube = map_from_land(cube_land)

CASE ( 'fwetl' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(toppdm%fwetl_soilt,ainfo))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for fwetl with nsoilt>1. " // &
                   "Use toppdm%fwetl_soilt instead.")
  END IF

CASE ( 'fwetl_soilt' )
  cube_land = cube_from_array(toppdm%fwetl_soilt)
  cube = map_from_land(cube_land)

CASE ( 'gpp_gb' )
  cube_land = cube_from_array(trifctltype%gpp_gb)
  cube = map_from_land(cube_land)

CASE ( 'lai_gb' )
  ! n.b. lai_gb is calculated here to make sure it is
  ! consistent with lai_pft after lai_pft is updated
  workspace_land(:) = SUM(progs%lai_pft * ainfo%frac_surft(:,1:npft), 2)
  cube_land = cube_from_array(workspace_land)
  cube = map_from_land(cube_land)

CASE ( 'et_stom_gb' )
  cube = cube_from_array(sf_diag%et_stom_ij)

CASE ( 'et_stom' )
  cube_land = cube_from_array(sf_diag%et_stom_surft)
  cube = map_from_land(cube_land)

CASE ( 'fprf' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(sf_diag%fprf)
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output fprf with tiled soil." )
  END IF

CASE ( 'fsth' )
  IF ( nsoilt == 1) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = sf_diag%fsth(:,k)
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output fsth with tiled soil." )
  END IF

CASE ( 'ftemp' )
  IF ( nsoilt == 1) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = sf_diag%ftemp(:,k)
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output ftemp with tiled soil."  )
  END IF

CASE ( 'gs' )
  cube_land = cube_from_array(progs%gs_gb)
  cube = map_from_land(cube_land)

CASE ( 'hf_snow_melt' )
  cube_land = cube_from_array(hf_snow_melt_gb)
  cube = map_from_land(cube_land)

CASE ( 'land_index' )
  cube_land = cube_from_array(REAL(ainfo%land_index))
  cube = map_from_land(cube_land)

CASE ( 'lice_index' )
  cube_land = cube_from_array(REAL(ainfo%lice_index))
  cube = map_from_land(cube_land)

CASE ( 'lit_c_mean' )
  cube_land = cube_from_array(trifctltype%lit_c_mn_gb)
  cube = map_from_land(cube_land)

CASE ( 'lw_net' )
  ! Calculate gridbox mean upwards longwave
  cube_land = cube_from_array(                                                &
    sbcon * surftiles_to_gbm(emis_surft * progs%tstar_surft**4, ainfo)        &
  )
  upward_cube = map_from_land(cube_land)
  CALL cube_free(cube_land)

  ! Calculate gridbox mean downward longwave
  cube_land = cube_from_array(surftiles_to_gbm(emis_surft, ainfo))
  emis_cube = map_from_land(cube_land)
  lw_down_cube = cube_from_array(lw_down_ij)
  downward_cube = emis_cube * lw_down_cube

  ! Now calculate the net flux
  cube = downward_cube - upward_cube

  ! Free work cubes
  CALL cube_free(upward_cube)
  CALL cube_free(emis_cube)
  CALL cube_free(lw_down_cube)
  CALL cube_free(downward_cube)

CASE ( 'lw_up' )
  cube_land = cube_from_array(                                                &
    sbcon * surftiles_to_gbm(emis_surft * progs%tstar_surft**4, ainfo)        &
  )
  cube = map_from_land(cube_land)

CASE ( 'npp_gb' )
  cube_land = cube_from_array(trifctltype%npp_gb)
  cube = map_from_land(cube_land)

CASE ( 'qbase' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(toppdm%qbase_soilt,ainfo))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for qbase with nsoilt>1. " // &
                   "Use toppdm%qbase_soilt instead.")
  END IF

CASE ( 'qbase_soilt' )
  cube_land = cube_from_array(toppdm%qbase_soilt)
  cube = map_from_land(cube_land)

CASE ( 'qbase_zw' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(toppdm%qbase_zw_soilt,ainfo))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for qbase_zw with " //        &
                   "nsoilt>1. Use toppdm%qbase_zw_soilt instead.")
  END IF

CASE ( 'qbase_zw_soilt' )
  cube_land = cube_from_array(toppdm%qbase_zw_soilt)
  cube = map_from_land(cube_land)

CASE ( 'rad_net' )
  cube_land = cube_from_array(surftiles_to_gbm(radnet_surft, ainfo))
  cube = map_from_land(cube_land)

CASE ( 'resp_p_gb' )
  cube_land = cube_from_array(trifctltype%resp_p_gb)
  cube = map_from_land(cube_land)

CASE ( 'resp_s_gb' )
  IF ( nsoilt == 1) THEN
    !Case for a single soil tile
    cube_land = cube_from_array(SUM(                                          &
                                SUM(trifctltype%resp_s_soilt(:,1,:,:), 3), 2))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output resp_s_gb for soil tiled runs.")
  END IF

CASE ( 'resp_s_dr_out' )
  ! HACK: We only output the total respiration for now
  cube_land = cube_from_array(SUM(trifctltype%resp_s_dr_out_gb(:,:,5), 2))
  cube = map_from_land(cube_land)

CASE ( 'resp_s_diag' )
  cube_land = cube_from_array(trif_vars%resp_s_diag_gb(:,:,1:4))
  cube = map_from_land(cube_land)

CASE ( 'resp_s_pot_diag' )
  cube_land = cube_from_array(trif_vars%resp_s_pot_diag_gb(:,:,1:4))
  cube = map_from_land(cube_land)

CASE ( 'immob_n' )
  cube_land = cube_from_array(trif_vars%immob_n_gb(:,:,1:4))
  cube = map_from_land(cube_land)

CASE ( 'immob_n_pot' )
  cube_land = cube_from_array(trif_vars%immob_n_pot_gb(:,:,1:4))
  cube = map_from_land(cube_land)

CASE ( 'minl_n' )
  cube_land = cube_from_array(trif_vars%minl_n_gb(:,:,1:4))
  cube = map_from_land(cube_land)

CASE ( 'minl_n_pot' )
  cube_land = cube_from_array(trif_vars%minl_n_pot_gb(:,:,1:4))
  cube = map_from_land(cube_land)

CASE ( 'ns' )
  SELECT CASE ( soil_bgc_model )
  CASE ( soil_model_rothc )
    cube_land = cube_from_array(progs%ns_pool_gb(:,:,:))
  CASE ( soil_model_ecosse )
    IF (nsoilt == 1) THEN
      ! Take the organic pools.
      cube_land = cube_from_array(soilecosse%n_soil_pool_soilt(:,1,:,1:dim_cs1))
    ELSE
      CALL log_fatal("extract_var",                                           &
                     "Cannot output ns when nsoilt > 1")
    END IF
  END SELECT
  cube = map_from_land(cube_land)

CASE ( 'fN' )
  cube_land = cube_from_array(trif_vars%fn_gb(:,:))
  cube = map_from_land(cube_land)

CASE ( 'soil_CN' )
  IF ( nsoilt == 1) THEN
    !Case for a single soil tile
    cube_land = cube_from_array(progs%cs_pool_soilt(:,1,:,:) /                &
                                progs%ns_pool_gb(:,:,:))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output soil_CN for soil tiled runs.")
  END IF

CASE ( 'resp_s_diag_gb' )
  cube_land = cube_from_array(SUM(trif_vars%resp_s_diag_gb(:,:,5), 2))
  cube = map_from_land(cube_land)

CASE ( 'resp_s_pot_diag_gb' )
  cube_land = cube_from_array(SUM(trif_vars%resp_s_pot_diag_gb(:,:,5), 2))
  cube = map_from_land(cube_land)

CASE ( 'immob_n_gb' )
  cube_land = cube_from_array(SUM(trif_vars%immob_n_gb(:,:,5), 2))
  cube = map_from_land(cube_land)

CASE ( 'immob_n_pot_gb' )
  cube_land = cube_from_array(SUM(trif_vars%immob_n_pot_gb(:,:,5), 2))
  cube = map_from_land(cube_land)

CASE ( 'minl_n_gb' )
  cube_land = cube_from_array(SUM(trif_vars%minl_n_gb(:,:,5), 2))
  cube = map_from_land(cube_land)

CASE ( 'minl_n_pot_gb' )
  cube_land = cube_from_array(SUM(trif_vars%minl_n_pot_gb(:,:,5), 2))
  cube = map_from_land(cube_land)

CASE ( 'ns_gb' )
  SELECT CASE ( soil_bgc_model )
  CASE ( soil_model_rothc )
    cube_land = cube_from_array(SUM(sum(progs%ns_pool_gb(:,:,1:4),dim = 3), 2))
  CASE ( soil_model_ecosse )
    IF (nsoilt == 1) THEN
      ! Sum over levels and pools.
      cube_land = cube_from_array(SUM(                                        &
                                SUM(soilecosse%n_soil_pool_soilt(:,1,:,:),3),2))
    ELSE
      CALL log_fatal("extract_var",                                           &
                     "Cannot output ns_gb when nsoilt > 1")
    END IF
  END SELECT
  cube = map_from_land(cube_land)

CASE ( 'soil_CN_gb' )
  IF ( nsoilt == 1) THEN
    cube_land =                                                               &
      cube_from_array(SUM(sum(progs%cs_pool_soilt(:,1,:,1:4),dim = 3), 2) /   &
      SUM(sum(progs%ns_pool_gb(:,:,1:4),dim = 3), 2))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output soil_CN_gb for soil tiled runs.")
  END IF

CASE ( 'runoff' )
  cube_land = cube_from_array(sub_surf_roff_gb(:) + surf_roff_gb(:))
  cube = map_from_land(cube_land)

CASE ( 'sat_excess_roff' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(toppdm%dun_roff_soilt,ainfo))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for sat_excess_roff with " // &
                   "nsoilt>1. Use sat_excess_roff_soilt instead.")
  END IF

CASE ( 'sat_excess_roff_soilt' )
  cube_land = cube_from_array(toppdm%dun_roff_soilt)
  cube = map_from_land(cube_land)

CASE ( 'smc_avail_top' )
  IF ( nsoilt == 1) THEN
    workspace_land(:) = 0.0
    ! We maintain the depth of the top of the current layer as we go down
    ! through the soil column
    ztop = 0.0
    DO k = 1,sm_levels
      ! If the top of this layer is below where we want to calculate to,
      !we are done
      IF ( ztop >= zsmc ) EXIT

      ! Calculate the amount of this layer that we want to take into account
      dz = dzsoil(k)
      ! If the layer contains the calculation limit, update dz
      IF ( ztop + dz > zsmc ) dz = zsmc - ztop

      ! Add on the contribution for this layer
      DO j = 1,soil_pts
        i = ainfo%soil_index(j)
        workspace_land(i) = workspace_land(i) + rho_water * dz *              &
                            MAX(0.0,                                          &
                            psparms%sthu_soilt(i,1,k)                         &
                            * psparms%smvcst_soilt(i,1,k)                     &
                            - psparms%smvcwt_soilt(i,1,k))
      END DO
      ztop = ztop + dzsoil(k)
    END DO
    ! Map the result onto the full grid
    cube_land = cube_from_array(workspace_land)
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for"                   //     &
                   "smc_avail_top with nsoilt>1: "                     //     &
                   "soil-tiled diagnostic not yet implemented.")
  END IF

CASE ( 'smc_avail_tot' )
  IF ( nsoilt == 1) THEN
    ! This is the same as smc_avail_top, but for whole column
    workspace_land(:) = 0.0
    DO k = 1,sm_levels
      DO j = 1,soil_pts
        i = ainfo%soil_index(j)
        workspace_land(i) = workspace_land(i) + rho_water * dzsoil(k) *       &
                            MAX(0.0,                                          &
                            psparms%sthu_soilt(i,1,k) *                       &
                            psparms%smvcst_soilt(i,1,k) -                     &
                            psparms%smvcwt_soilt(i,1,k))
      END DO
    END DO
    ! Map the result onto the full grid
    cube_land = cube_from_array(workspace_land)
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for"                   //     &
                   "smc_avail_tot with nsoilt>1: "                     //     &
                   "soil-tiled diagnostic not yet implemented.")
  END IF

CASE ( 'smc_tot' )
  workspace_land(:) = 0.0
  IF ( nsoilt == 1) THEN
    DO k = 1,sm_levels
      DO j = 1,soil_pts
        i = ainfo%soil_index(j)
        workspace_land(i) = workspace_land(i) + rho_water * dzsoil(k) *       &
                            MAX(0.0,                                          &
                            (psparms%sthu_soilt(i,1,k) +                      &
                            psparms%sthf_soilt(i,1,k)) *                      &
                            psparms%smvcst_soilt(i,1,k))
      END DO
    END DO
  ELSE
    !We can't take a gbm of psparms%smvcst_soilt so add up manually
    DO k = 1,sm_levels
      DO m = 1,nsoilt
        DO j = 1,soil_pts
          i = ainfo%soil_index(j)
          workspace_land(i) = workspace_land(i) + (ainfo%frac_soilt(i,m) *    &
                              (rho_water * dzsoil(k) *                        &
                              MAX(0.0,                                        &
                              (psparms%sthu_soilt(i,m,k) +                    &
                              psparms%sthf_soilt(i,m,k)) *                    &
                              psparms%smvcst_soilt(i,m,k))))
        END DO !j
      END DO !m
    END DO !k
  END IF
  ! Map the result onto the full grid
  cube_land = cube_from_array(workspace_land)
  cube = map_from_land(cube_land)

CASE ( 'snomlt_sub_htf' )
  cube_land = cube_from_array(snomlt_sub_htf_gb)
  cube = map_from_land(cube_land)

CASE ( 'snow_can_gb' )
  ! Only include tiles where canopy snow model is used
  cube_land = cube_from_array(surftiles_to_gbm(progs%snow_surft, ainfo, canSnowTile))
  cube = map_from_land(cube_land)

CASE ( 'snow_depth_gb' )
  cube_land = cube_from_array(surftiles_to_gbm(progs%snowdepth_surft, ainfo))
  cube = map_from_land(cube_land)

CASE ( 'snow_frac' )
  ! Sum frac over tiles with snow.
  workspace_land(:) = 0.0
  IF ( l_aggregate ) THEN
    WHERE ( progs%snow_surft(:,1) + progs%snow_grnd_surft(:,1) > EPSILON(1.0) )
      workspace_land(:) = 1.0
    END WHERE
  ELSE
    workspace_land(:) = SUM(ainfo%frac_surft, 2, progs%snow_surft +           &
                            progs%snow_grnd_surft > EPSILON(1.0))
  END IF
  cube_land = cube_from_array(workspace_land)
  cube = map_from_land(cube_land)

CASE ( 'snow_grnd_gb' )
  ! Only include tiles where canopy snow model is used
  cube_land = cube_from_array(surftiles_to_gbm(progs%snow_grnd_surft,         &
  ainfo,canSnowTile))
  cube = map_from_land(cube_land)

CASE ( 'snow_ice_gb' )
  ! Calculate sum of sice along the snow layers dimension
  DO n = 1,nsurft
    DO j = 1,surft_pts(n)
      i = ainfo%surft_index(j,n)
      workspace_surft(i,n) = SUM(progs%sice_surft(i,n,1:progs%nsnow_surft(i,n)))
    END DO
  END DO
  cube_land = cube_from_array(surftiles_to_gbm(workspace_surft, ainfo))
  cube = map_from_land(cube_land)

CASE ( 'snow_liq_gb' )
  ! Calculate sum of sliq along the snow layers dimension
  DO n = 1,nsurft
    DO j = 1,surft_pts(n)
      i = ainfo%surft_index(j,n)
      workspace_surft(i,n) = SUM(progs%sliq_surft(i,n,1:progs%nsnow_surft(i,n)))
    END DO
  END DO
  cube_land = cube_from_array(surftiles_to_gbm(workspace_surft, ainfo))
  cube = map_from_land(cube_land)

CASE ( 'snow_melt_gb' )
  cube_land = cube_from_array(snow_melt_gb)
  cube = map_from_land(cube_land)

CASE ( 'soil_index' )
  DO l = 1,soil_pts
    workspace_land(ainfo%soil_index(l)) = REAL(ainfo%soil_index(l))
  END DO
  cube_land = cube_from_array(workspace_land)
  cube = map_from_land(cube_land)

CASE ( 'sthzw' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(toppdm%sthzw_soilt,ainfo))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for sthzw with nsoilt>1: " // &
                   "calculation is unphysical. Use toppdm%sthzw_soilt instead.")
  END IF

CASE ( 'sthzw_soilt' )
  cube_land = cube_from_array(toppdm%sthzw_soilt)
  cube = map_from_land(cube_land)

CASE ( 'sub_surf_roff' )
  cube_land = cube_from_array(sub_surf_roff_gb)
  cube = map_from_land(cube_land)

CASE ( 'surf_roff' )
  cube_land = cube_from_array(surf_roff_gb)
  cube = map_from_land(cube_land)

CASE ( 'swet_liq_tot' )
  IF ( nsoilt == 1) THEN
    ! Divide column unfrozen moisture content by saturated moisture content
    workspace_land(:) = 0.0
    DO j = 1,soil_pts
      i = ainfo%soil_index(j)
      workspace_land(i) = SUM(dzsoil(:) * MAX(0.0,                            &
                        psparms%sthu_soilt(i,1,:) *                           &
                        psparms%smvcst_soilt(i,1,:)))                         &
                        / SUM(dzsoil(:) * MAX(0.0, psparms%smvcst_soilt(i,1,:)))
    END DO
    ! Map the result onto the full grid
    cube_land = cube_from_array(workspace_land)
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for"                   //     &
                   "swet_liq_tot with nsoilt>1: "                      //     &
                   "soil tiled diagnostic not yet implemented.")
  END IF

CASE ( 'swet_tot' )
  IF ( nsoilt == 1) THEN
    ! Divide column total moisture content by saturated moisture content
    workspace_land(:) = 0.0
    DO j = 1,soil_pts
      i = ainfo%soil_index(j)
      workspace_land(i) = SUM(dzsoil(:) *                                     &
                            MAX(0.0,                                          &
                            (psparms%sthu_soilt(i,1,:) +                      &
                            psparms%sthf_soilt(i,1,:)) *                      &
                            psparms%smvcst_soilt(i,1,:)))                     &
                            /SUM(dzsoil(:) * MAX(0.0,                         &
                            psparms%smvcst_soilt(i,1,:)))
    END DO
    ! Map the result onto the full grid
    cube_land = cube_from_array(workspace_land)
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for "                   //    &
                   "swet_tot wth nsoilt>1: "                            //    &
                   "soil tiled diagnostic not yet implemented.")
  END IF

CASE ( 'sw_net' )
  ! Calculate the albedo as used in subroutine control when calculating
  ! the net shortwave on tiles
  ! Here we take the average of diffuse albedos in VIS and NIR.
  DO l = 1,land_pts
    j = (ainfo%land_index(l) - 1) / t_i_length + 1
    i = ainfo%land_index(l) - (j-1) * t_i_length
    workspace_land(l) = ( 1.0 - (                                             &
                          wght_alb(1) * land_albedo_ij(i,j,1) +               &
                          wght_alb(2) * land_albedo_ij(i,j,2) +               &
                          wght_alb(3) * land_albedo_ij(i,j,3) +               &
                          wght_alb(4) * land_albedo_ij(i,j,4)                 &
                ) ) * sw_down_ij(i,j)

  END DO
  cube_land = cube_from_array(workspace_land)
  cube = map_from_land(cube_land)

CASE ( 'NDVI_land' )
  ! Calculate NDVI on tiles
  ! Using the average of direct & diffuse albedos in VIS and NIR.
  DO l = 1,land_pts
    j = (ainfo%land_index(l) - 1) / t_i_length + 1
    i = ainfo%land_index(l) - (j-1) * t_i_length
    workspace_land(l) = (((land_albedo_ij(i,j,3) + land_albedo_ij(i,j,4)) -   &
                         (land_albedo_ij(i,j,1) + land_albedo_ij(i,j,2))) /   &
                        ((land_albedo_ij(i,j,3) + land_albedo_ij(i,j,4)) +    &
                         (land_albedo_ij(i,j,1) + land_albedo_ij(i,j,2))))

  END DO
  cube_land = cube_from_array(workspace_land)
  cube = map_from_land(cube_land)


CASE ( 'tfall' )
  cube_land = cube_from_array(tot_tfall_gb)
  cube = map_from_land(cube_land)

CASE ( 'trad' )
  ! Assuming emissivity=1.
  cube_land = cube_from_array((                                               &
                        surftiles_to_gbm(progs%tstar_surft**4, ainfo) )**0.25)
  cube = map_from_land(cube_land)

CASE ( 'zw' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(toppdm%zw_soilt,ainfo))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for zw with nsoilt>1. " //    &
                   "Use toppdm%zw_soilt instead.")
  END IF

CASE ( 'zw_soilt' )
  cube_land = cube_from_array(toppdm%zw_soilt)
  cube = map_from_land(cube_land)

CASE ( 'c_veg' )
  cube_land = cube_from_array(trifctltype%c_veg_pft)
  cube = map_from_land(cube_land)

CASE ( 'fapar' )
  cube_land = cube_from_array(trif_vars%fapar_diag_pft)
  cube = map_from_land(cube_land)

CASE ( 'apar' )
  cube_land = cube_from_array(trif_vars%apar_diag_pft)
  cube = map_from_land(cube_land)

CASE ( 'apar_gb' )
  cube_land = cube_from_array(trif_vars%apar_diag_gb)
  cube = map_from_land(cube_land)

CASE ( 'fao_et0' )
  cube_land = cube_from_array(trif_vars%fao_et0)
  cube = map_from_land(cube_land)

CASE ( 'canht' )
  cube_land = cube_from_array(progs%canht_pft)
  cube = map_from_land(cube_land)

CASE ( 'flux_o3_stom' )
  cube_land = cube_from_array(flux_o3_pft)
  cube = map_from_land(cube_land)

CASE ( 'fsmc' )
  cube_land = cube_from_array(fsmc_pft)
  cube = map_from_land(cube_land)

CASE ( 'g_leaf' )
  cube_land = cube_from_array(trifctltype%g_leaf_pft)
  cube = map_from_land(cube_land)

CASE ( 'g_leaf_day' )
  cube_land = cube_from_array(trifctltype%g_leaf_day_pft)
  cube = map_from_land(cube_land)

CASE ( 'g_leaf_dr_out' )
  cube_land = cube_from_array(trifctltype%g_leaf_dr_out_pft)
  cube = map_from_land(cube_land)

CASE ( 'g_leaf_phen' )
  cube_land = cube_from_array(trifctltype%g_leaf_phen_pft)
  cube = map_from_land(cube_land)

CASE ( 'gpp' )
  cube_land = cube_from_array(trifctltype%gpp_pft)
  cube = map_from_land(cube_land)

CASE ( 'lai' )
  cube_land = cube_from_array(progs%lai_pft)
  cube = map_from_land(cube_land)

CASE ( 'lai_phen' )
  cube_land = cube_from_array(trifctltype%lai_phen_pft)
  cube = map_from_land(cube_land)

CASE ( 'lit_c' )
  cube_land = cube_from_array(trifctltype%lit_c_pft)
  cube = map_from_land(cube_land)

CASE ( 'lit_c_ag' )
  cube_land = cube_from_array(trif_vars%lit_c_ag_pft)
  cube = map_from_land(cube_land)

CASE ( 'lit_c_orig' )
  cube_land = cube_from_array(trif_vars%lit_c_orig_pft)
  cube = map_from_land(cube_land)

CASE ( 'harvest' )
  cube_land = cube_from_array(trif_vars%harvest_pft)
  cube = map_from_land(cube_land)

CASE ( 'harvest_gb' )
  cube_land = cube_from_array(trif_vars%harvest_gb)
  cube = map_from_land(cube_land)

CASE ( 'harvest_n' )
  cube_land = cube_from_array(trif_vars%harvest_n_pft)
  cube = map_from_land(cube_land)

CASE ( 'harvest_n_gb' )
  cube_land = cube_from_array(trif_vars%harvest_n_gb)
  cube = map_from_land(cube_land)

CASE ( 'n_fertiliser' )
  cube_land = cube_from_array(trif_vars%n_fertiliser_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_fertiliser_gb' )
  cube_land = cube_from_array(trif_vars%n_fertiliser_gb)
  cube = map_from_land(cube_land)

CASE ( 'root_abandon' )
  cube_land = cube_from_array(trif_vars%root_abandon_pft)
  cube = map_from_land(cube_land)

CASE ( 'root_abandon_gb' )
  cube_land = cube_from_array(trif_vars%root_abandon_gb)
  cube = map_from_land(cube_land)

CASE ( 'root_abandon_n' )
  cube_land = cube_from_array(trif_vars%root_abandon_n_pft)
  cube = map_from_land(cube_land)

CASE ( 'root_abandon_n_gb' )
  cube_land = cube_from_array(trif_vars%root_abandon_n_gb)
  cube = map_from_land(cube_land)

CASE ( 'lit_n_ag' )
  cube_land = cube_from_array(trif_vars%lit_n_ag_pft)
  cube = map_from_land(cube_land)

CASE ( 'lit_n_orig' )
  cube_land = cube_from_array(trif_vars%lit_n_orig_pft)
  cube = map_from_land(cube_land)

CASE ( 'npp_dr_out' )
  cube_land = cube_from_array(trifctltype%npp_dr_out_pft)
  cube = map_from_land(cube_land)

CASE ( 'npp' )
  cube_land = cube_from_array(trifctltype%npp_pft)
  cube = map_from_land(cube_land)

CASE ( 'o3_exp_fac' )
  cube_land = cube_from_array(fo3_pft)
  cube = map_from_land(cube_land)

CASE ( 'co2_mmr' )
  co2_tmp(:)=co2_mmr
  cube_land = cube_from_array(co2_tmp)
  cube = map_from_land(cube_land)

CASE ( 'frac_agr' )
  cube_land = cube_from_array(trifctltype%frac_agr_gb)
  cube = map_from_land(cube_land)

CASE ( 'frac_past' )
  cube_land = cube_from_array(trif_vars%frac_past_gb)
  cube = map_from_land(cube_land)

CASE ( 'pc_s' )
  cube_land = cube_from_array(trif_vars%pc_s_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_leaf' )
  cube_land = cube_from_array(trif_vars%n_leaf_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_root' )
  cube_land = cube_from_array(trif_vars%n_root_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_stem' )
  cube_land = cube_from_array(trif_vars%n_stem_pft)
  cube = map_from_land(cube_land)

CASE ( 'lai_bal' )
  cube_land = cube_from_array(trif_vars%lai_bal_pft)
  cube = map_from_land(cube_land)

CASE ( 'resp_p' )
  cube_land = cube_from_array(trifctltype%resp_p_pft)
  cube = map_from_land(cube_land)

CASE ( 'resp_r' )
  cube_land = cube_from_array(trif_vars%resp_r_pft)
  cube = map_from_land(cube_land)

CASE ( 'resp_l' )
  cube_land = cube_from_array(trif_vars%resp_l_pft)
  cube = map_from_land(cube_land)

CASE ( 'resp_w_dr_out' )
  cube_land = cube_from_array(trifctltype%resp_w_dr_out_pft)
  cube = map_from_land(cube_land)

CASE ( 'resp_w' )
  cube_land = cube_from_array(trifctltype%resp_w_pft)
  cube = map_from_land(cube_land)

CASE ( 'cs' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array(progs%cs_pool_soilt(:,1,:,:))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for cs when nsoilt > 1.")
  END IF

CASE ( 'resp_s_to_atmos_gb' )
  SELECT CASE ( soil_bgc_model )
  CASE ( soil_model_rothc )
    cube_land = cube_from_array( SUM(trif_vars%resp_s_to_atmos_gb(:,:),2) )
  CASE ( soil_model_ecosse )
    ! Convert units from kg m-2 s-1 to kg m-2 (360 days)-1.
    cube_land = cube_from_array( soilecosse%co2_soil_gb(:) * secs_per_360days &
                                 / dt_soilc )
  CASE DEFAULT
    CALL log_fatal("extract_var",                                             &
                   "Cannot output resps_to_atmos_gb - no code for " //        &
                   "this soil model.")
  END SELECT
  cube = map_from_land(cube_land)

CASE ( 'resp_s' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array(trifctltype%resp_s_soilt(:,1,:,:))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for resp_s when nsoilt > 1.")
  END IF

CASE ( 'wood_prod_fast' )
  cube_land = cube_from_array(progs%wood_prod_fast_gb)
  cube = map_from_land(cube_land)

CASE ( 'wood_prod_med' )
  cube_land = cube_from_array(progs%wood_prod_med_gb)
  cube = map_from_land(cube_land)

CASE ( 'wood_prod_slow' )
  cube_land = cube_from_array(progs%wood_prod_slow_gb)
  cube = map_from_land(cube_land)

CASE ( 'WP_fast_in' )
  cube_land = cube_from_array(trif_vars%WP_fast_in_gb)
  cube = map_from_land(cube_land)

CASE ( 'WP_med_in' )
  cube_land = cube_from_array(trif_vars%WP_med_in_gb)
  cube = map_from_land(cube_land)

CASE ( 'WP_slow_in' )
  cube_land = cube_from_array(trif_vars%WP_slow_in_gb)
  cube = map_from_land(cube_land)

CASE ( 'WP_fast_out' )
  cube_land = cube_from_array(trif_vars%WP_fast_out_gb)
  cube = map_from_land(cube_land)

CASE ( 'WP_med_out' )
  cube_land = cube_from_array(trif_vars%WP_med_out_gb)
  cube = map_from_land(cube_land)

CASE ( 'WP_slow_out' )
  cube_land = cube_from_array(trif_vars%WP_slow_out_gb)
  cube = map_from_land(cube_land)

CASE ( 'cnsrv_carbon_veg2' )
  cube_land = cube_from_array(REAL(trif_vars%cnsrv_carbon_veg2_gb))
  cube = map_from_land(cube_land)

CASE ( 'cnsrv_carbon_triffid' )
  cube_land = cube_from_array(REAL(trif_vars%cnsrv_carbon_triffid_gb))
  cube = map_from_land(cube_land)

CASE ( 'cnsrv_veg_triffid' )
  cube_land = cube_from_array(REAL(trif_vars%cnsrv_veg_triffid_gb))
  cube = map_from_land(cube_land)

CASE ( 'cnsrv_soil_triffid' )
  cube_land = cube_from_array(REAL(trif_vars%cnsrv_soil_triffid_gb))
  cube = map_from_land(cube_land)

CASE ( 'cnsrv_prod_triffid' )
  cube_land = cube_from_array(REAL(trif_vars%cnsrv_prod_triffid_gb))
  cube = map_from_land(cube_land)

CASE ( 'cnsrv_nitrogen_triffid' )
  cube_land = cube_from_array(REAL(trif_vars%cnsrv_nitrogen_triffid_gb))
  cube = map_from_land(cube_land)

CASE ( 'cnsrv_vegN_triffid' )
  cube_land = cube_from_array(REAL(trif_vars%cnsrv_vegN_triffid_gb))
  cube = map_from_land(cube_land)

CASE ( 'cnsrv_soilN_triffid' )
  cube_land = cube_from_array(REAL(trif_vars%cnsrv_soilN_triffid_gb))
  cube = map_from_land(cube_land)

CASE ( 'cnsrv_n_inorg_triffid' )
  cube_land = cube_from_array(REAL(trif_vars%cnsrv_n_inorg_triffid_gb))
  cube = map_from_land(cube_land)

CASE ( 'lit_c_fire' )
  cube_land = cube_from_array(REAL(trif_vars%lit_c_fire_pft))
  cube = map_from_land(cube_land)

CASE ( 'burnt_carbon_dpm' )
  cube_land = cube_from_array(REAL(trif_vars%burnt_carbon_dpm))
  cube = map_from_land(cube_land)

CASE ( 'lit_n_fire' )
  cube_land = cube_from_array(REAL(trif_vars%lit_n_fire_pft))
  cube = map_from_land(cube_land)

CASE ( 'burnt_carbon_rpm' )
  cube_land = cube_from_array(REAL(trif_vars%burnt_carbon_rpm))
  cube = map_from_land(cube_land)

CASE ( 'veg_c_fire_emission_gb' )
  cube_land = cube_from_array(REAL(trif_vars%veg_c_fire_emission_gb))
  cube = map_from_land(cube_land)

CASE ( 'veg_c_fire_emission_pft' )
  cube_land = cube_from_array(REAL(trif_vars%veg_c_fire_emission_pft))
  cube = map_from_land(cube_land)

CASE ( 'n_inorg' )
  cube_land = cube_from_array(progs%n_inorg_soilt_lyrs(:,1,:))
  cube = map_from_land(cube_land)

CASE ( 'n_inorg_avail_pft' )
  cube_land = cube_from_array(progs%n_inorg_avail_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_inorg_gb' )
  IF ( nsoilt == 1 ) THEN
    IF ( soil_bgc_model == soil_model_rothc ) THEN
      cube_land = cube_from_array( SUM( progs%n_inorg_soilt_lyrs(:,1,:), 2 ) )
    ELSE IF (  soil_bgc_model == soil_model_ecosse ) THEN
      ! Sum over levels.
      cube_land = cube_from_array( SUM(                                       &
                                soilecosse%n_soil_pool_soilt(:,1,:,i_amm) +   &
                                soilecosse%n_soil_pool_soilt(:,1,:,i_nit), 2 ) )
    END IF
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for n_inorg_gb with " //      &
                   "nsoilt>1.")
  END IF

CASE ( 'deposition_n' )
  cube_land = cube_from_array(trif_vars%deposition_n_gb)
  cube = map_from_land(cube_land)

CASE ( 'substr_ch4' )
  cube_land = cube_from_array(progs%substr_ch4)
  cube = map_from_land(cube_land)

CASE ( 'mic_ch4' )
  cube_land = cube_from_array(progs%mic_ch4)
  cube = map_from_land(cube_land)

CASE ( 'mic_act_ch4' )
  cube_land = cube_from_array(progs%mic_act_ch4)
  cube = map_from_land(cube_land)

CASE ( 'acclim_ch4' )
  cube_land = cube_from_array(progs%acclim_ch4)
  cube = map_from_land(cube_land)

CASE ( 'g_burn_pft' )
  cube_land = cube_from_array(REAL(trif_vars%g_burn_pft))
  cube = map_from_land(cube_land)

CASE ( 'g_burn_gb' )
  cube_land = cube_from_array(REAL(trif_vars%g_burn_gb))
  cube = map_from_land(cube_land)

CASE ( 'leafC' )
  cube_land = cube_from_array(trif_vars%leafc_pft)
  cube = map_from_land(cube_land)

CASE ( 'rootC' )
  cube_land = cube_from_array(trif_vars%rootc_pft)
  cube = map_from_land(cube_land)

CASE ( 'stemC' )
  cube_land = cube_from_array(trif_vars%stemc_pft)
  cube = map_from_land(cube_land)

CASE ( 'woodC' )
  cube_land = cube_from_array(trif_vars%woodc_pft)
  cube = map_from_land(cube_land)

CASE ( 'dleaf' )
  cube_land = cube_from_array(trif_vars%dleaf_pft)
  cube = map_from_land(cube_land)

CASE ( 'droot' )
  cube_land = cube_from_array(trif_vars%droot_pft)
  cube = map_from_land(cube_land)

CASE ( 'dwood' )
  cube_land = cube_from_array(trif_vars%dwood_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_uptake' )
  cube_land = cube_from_array(trif_vars%n_uptake_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_uptake_gb' )
  cube_land = cube_from_array(trif_vars%n_uptake_gb)
  cube = map_from_land(cube_land)

CASE ( 'leaf_litC' )
  cube_land = cube_from_array(trif_vars%leaf_litc_pft)
  cube = map_from_land(cube_land)

CASE ( 'root_litC' )
  cube_land = cube_from_array(trif_vars%root_litc_pft)
  cube = map_from_land(cube_land)

CASE ( 'wood_litC' )
  cube_land = cube_from_array(trif_vars%wood_litc_pft)
  cube = map_from_land(cube_land)

CASE ( 'leaf_litN' )
  cube_land = cube_from_array(trif_vars%leaf_litn_pft)
  cube = map_from_land(cube_land)

CASE ( 'root_litN' )
  cube_land = cube_from_array(trif_vars%root_litn_pft)
  cube = map_from_land(cube_land)

CASE ( 'wood_litN' )
  cube_land = cube_from_array(trif_vars%wood_litn_pft)
  cube = map_from_land(cube_land)

CASE ( 'litterC' )
  cube_land = cube_from_array(trif_vars%litterc_pft)
  cube = map_from_land(cube_land)

CASE ( 'litterN' )
  cube_land = cube_from_array(trif_vars%littern_pft)
  cube = map_from_land(cube_land)

CASE ( 'exudates' )
  cube_land = cube_from_array(trif_vars%exudates_pft)
  cube = map_from_land(cube_land)

CASE ( 'exudates_gb' )
  cube_land = cube_from_array(trif_vars%exudates_gb)
  cube = map_from_land(cube_land)

CASE ( 'npp_n' )
  cube_land = cube_from_array(trif_vars%npp_n)
  cube = map_from_land(cube_land)

CASE ( 'npp_n_gb' )
  cube_land = cube_from_array(trif_vars%npp_n_gb)
  cube = map_from_land(cube_land)

CASE ( 'n_demand' )
  cube_land = cube_from_array(trif_vars%n_demand_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_fix' )
  cube_land = cube_from_array(trif_vars%n_fix_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_fix_gb' )
  cube_land = cube_from_array(trif_vars%n_fix_gb)
  cube = map_from_land(cube_land)

CASE ( 'n_gas' )
  cube_land = cube_from_array(trif_vars%n_gas_gb(:,:))
  cube = map_from_land(cube_land)

CASE ( 'n_gas_gb' )
  SELECT CASE ( soil_bgc_model )
  CASE ( soil_model_rothc )
    cube_land = cube_from_array(SUM(trif_vars%n_gas_gb(:,:),2))
  CASE ( soil_model_ecosse )
    cube_land = cube_from_array( soilecosse%n2_denitrif_gb                    &
                                                + soilecosse%no_soil_gb       &
                                                + soilecosse%n2o_soil_gb )
  END SELECT
  cube = map_from_land(cube_land)

CASE ( 'n_leach' )
  SELECT CASE ( soil_bgc_model )
  CASE ( soil_model_rothc )
    IF (nsoilt == 1) THEN
      cube_land = cube_from_array(trif_vars%n_leach_soilt(:,1))
    ELSE
      CALL log_fatal("extract_var",                                           &
                     "Cannot output GB mean value for n_leach when " //       &
                     "nsoilt > 1.")
    END IF
  CASE ( soil_model_ecosse )
    cube_land = cube_from_array( soilecosse%n_leach_amm_gb(:) +               &
                                 soilecosse%n_leach_nit_gb(:) )
  END SELECT
  cube = map_from_land(cube_land)

CASE ( 'n_demand_gb' )
  cube_land = cube_from_array(trif_vars%n_demand_gb)
  cube = map_from_land(cube_land)

CASE ( 'n_demand_growth' )
  cube_land = cube_from_array(trif_vars%n_demand_growth_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_uptake_growth' )
  cube_land = cube_from_array(trif_vars%n_uptake_growth_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_demand_lit' )
  cube_land = cube_from_array(trif_vars%n_demand_lit_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_demand_spread' )
  cube_land = cube_from_array(trif_vars%n_demand_spread_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_veg' )
  cube_land = cube_from_array(trif_vars%n_veg_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_veg_gb' )
  cube_land = cube_from_array(trif_vars%n_veg_gb)
  cube = map_from_land(cube_land)

CASE ( 'n_loss' )
  cube_land = cube_from_array(trif_vars%n_loss_gb)
  cube = map_from_land(cube_land)

CASE ( 'dpm_ratio' )
  cube_land = cube_from_array(trif_vars%dpm_ratio_gb)
  cube = map_from_land(cube_land)

CASE ( 'dnveg' )
  cube_land = cube_from_array(trif_vars%dnveg_pft)
  cube = map_from_land(cube_land)

CASE ( 'dnveg_gb' )
  cube_land = cube_from_array(trif_vars%dnveg_gb)
  cube = map_from_land(cube_land)

CASE ( 'dcveg' )
  cube_land = cube_from_array(trif_vars%dcveg_pft)
  cube = map_from_land(cube_land)

CASE ( 'dcveg_gb' )
  cube_land = cube_from_array(trif_vars%dcveg_gb)
  cube = map_from_land(cube_land)

CASE ( 'lit_N' )
  cube_land = cube_from_array(trif_vars%lit_n_pft)
  cube = map_from_land(cube_land)

CASE ( 'lit_N_t' )
  cube_land = cube_from_array(trif_vars%lit_n_t_gb)
  cube = map_from_land(cube_land)

CASE ( 'lit_n' )
  cube_land = cube_from_array(trif_vars%lit_n_pft)
  cube = map_from_land(cube_land)

CASE ( 'lit_n_t' )
  cube_land = cube_from_array(trif_vars%lit_n_t_gb)
  cube = map_from_land(cube_land)

CASE ( 'con_rain' )
  cube = cube_from_array(con_rain_ij)

CASE ( 'con_snow' )
  cube = cube_from_array(con_snow_ij)

CASE ( 'cosz' )
  cube = cube_from_array(psparms%cosz_ij)

CASE ( 'diff_frac' )
  cube = cube_from_array(RESHAPE(jules_vars%diff_frac, (/ t_i_length, t_j_length /)))

CASE ( 'ecan_gb' )
  cube = cube_from_array(ecan_ij)

CASE ( 'ei_gb' )
  cube = cube_from_array(ei_ij)

CASE ( 'esoil_gb' )
  IF (nsoilt == 1) THEN
    cube = cube_from_array(esoil_ij_soilt(:,:,1))
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for esoil_gb when nsoilt>1.")
  END IF

CASE ( 'fqw_gb' )
  cube = cube_from_array(fqw_1_ij)

CASE ( 'ftl_gb' )
  cube = cube_from_array(ftl_1_ij)

CASE ( 'tau_gb' )
  cube = cube_from_array(sf_diag%tau_1)

CASE ( 'land_albedo_1' )
  cube = cube_from_array(land_albedo_ij(:,:,1))

CASE ( 'land_albedo_2' )
  cube = cube_from_array(land_albedo_ij(:,:,2))

CASE ( 'land_albedo_3' )
  cube = cube_from_array(land_albedo_ij(:,:,3))

CASE ( 'land_albedo_4' )
  cube = cube_from_array(land_albedo_ij(:,:,4))

CASE ( 'upward_sw_1' )
  cube = cube_from_array(land_albedo_ij(:,:,1) * sw_down_ij(:,:) * wght_alb(1))

CASE ( 'upward_sw_2' )
  cube = cube_from_array(land_albedo_ij(:,:,2) * sw_down_ij(:,:) * wght_alb(2))

CASE ( 'upward_sw_3' )
  cube = cube_from_array(land_albedo_ij(:,:,3) * sw_down_ij(:,:) * wght_alb(3))

CASE ( 'upward_sw_4' )
  cube = cube_from_array(land_albedo_ij(:,:,4) * sw_down_ij(:,:) * wght_alb(4))

CASE ( 'latent_heat' )
  cube = cube_from_array(sf_diag%latent_heat)

CASE ( 'ls_rain' )
  cube = cube_from_array(ls_rain_ij)

CASE ( 'ls_snow' )
  cube = cube_from_array(ls_snow_ij)

CASE ( 'lw_down' )
  cube = cube_from_array(lw_down_ij)

CASE ( 'precip' )
  cube = cube_from_array(ls_rain_ij + con_rain_ij +                           &
                         ls_snow_ij + con_snow_ij)

CASE ( 'pstar' )
  cube = cube_from_array(pstar_ij)

CASE ( 'q1p5m_gb' )
  cube = cube_from_array(sf_diag%q1p5m)

CASE ( 'qw1' )
  cube = cube_from_array(qw_1_ij)

CASE ( 'rainfall' )
  cube = cube_from_array(ls_rain_ij + con_rain_ij)

CASE ( 'snomlt_surf_htf' )
  cube = cube_from_array(sf_diag%snomlt_surf_htf)

CASE ( 'snowfall' )
  cube = cube_from_array(ls_snow_ij + con_snow_ij)

CASE ( 'snow_mass_gb' )
  !     Don't use the snow_mass_ij variable as that is calculated under control.
  cube_land = cube_from_array(surftiles_to_gbm(progs%snow_grnd_surft +        &
                                               progs%snow_surft, ainfo))
  cube = map_from_land(cube_land)

CASE ( 'surf_ht_flux_gb' )
  cube = cube_from_array(surf_ht_flux_ij)

CASE ( 'sw_down' )
  cube = cube_from_array(sw_down_ij)

CASE ( 't1p5m_gb' )
  cube = cube_from_array(sf_diag%t1p5m)

CASE ( 'taux1' )
  cube = cube_from_array(taux_1_ij)

CASE ( 'tauy1' )
  cube = cube_from_array(tauy_1_ij)

CASE ( 'tl1' )
  cube = cube_from_array(tl_1_ij)

CASE ( 'tstar_gb' )
  !     Don't use tstar as that is calculated part-way through a timestep.
  cube_land = cube_from_array(surftiles_to_gbm( progs%tstar_surft, ainfo ))
  cube = map_from_land(cube_land)

CASE ( 'u1' )
  cube = cube_from_array(u_1_ij)

CASE ( 'u10m' )
  cube = cube_from_array(sf_diag%u10m)

CASE ( 'v1' )
  cube = cube_from_array(v_1_ij)

CASE ( 'v10m' )
  cube = cube_from_array(sf_diag%v10m)

CASE ( 'wind' )
  cube = cube_from_array(SQRT(u_1_ij**2 + v_1_ij**2))

CASE ( 'dt_range' )
  cube = cube_from_array(diurnal_temperature_range_ij)

CASE ( 'ext' )
  IF ( nsoilt == 1 ) THEN
    cube_land = cube_from_array(ext_soilt(:,1,:))
    cube      = map_from_land(cube_land)
  ELSE
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(ext_soilt(:,:,k),ainfo)
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube      = map_from_land(cube_land)
  END IF

CASE ( 'ext_soilt' )
  cube_land = cube_from_array(ext_soilt)
  cube = map_from_land(cube_land)

CASE ( 'smcl' )
  IF ( nsoilt == 1 ) THEN
    cube_land = cube_from_array(progs%smcl_soilt(:,1,:))
    cube      = map_from_land(cube_land)
  ELSE
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(progs%smcl_soilt(:,:,k),ainfo)
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube      = map_from_land(cube_land)
  END IF

CASE ( 'smcl_soilt' )
  cube_land = cube_from_array(progs%smcl_soilt)
  cube = map_from_land(cube_land)

CASE ( 'soil_wet' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(psparms%sthu_soilt(:,:,k),ainfo) &
                          + soiltiles_to_gbm(psparms%sthf_soilt(:,:,k),ainfo)
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for soil_wet with " //        &
                   "with nsoilt>1: diagnostic not yet implemented.")
  END IF

CASE ( 'sthf' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(psparms%sthf_soilt(:,:,k),ainfo)
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for sthf with nsoilt>1: " //  &
                   "calculation is unphysical. Use psparms%sthf_soilt instead.")
  END IF

CASE ( 'sthf_soilt' )
  cube_land = cube_from_array(psparms%sthf_soilt)
  cube = map_from_land(cube_land)

CASE ( 'sthu' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(psparms%sthu_soilt(:,:,k),ainfo)
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for sthu with nsoilt>1: " //  &
                   "calculation is unphysical. Use psparms%sthu_soilt instead.")
  END IF

CASE ( 'sthu_soilt' )
  cube_land = cube_from_array(psparms%sthu_soilt)
  cube = map_from_land(cube_land)

CASE ( 't_soil' )
  DO k = 1, sm_levels
    workspace_levs(:,k) = soiltiles_to_gbm(progs%t_soil_soilt(:,:,k),ainfo)
  END DO
  cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
  cube = map_from_land(cube_land)

CASE ( 't_soil_soilt' )
  cube_land = cube_from_array(progs%t_soil_soilt)
  cube = map_from_land(cube_land)

CASE ( 'tsoil_deep' )
  cube_land = cube_from_array(progs%tsoil_deep_gb(:,1:ns_deep))
  cube = map_from_land(cube_land)

CASE ( 'alb_tile_1' )
  cube_land = cube_from_array(alb_surft(:,:,1))
  cube = map_from_land(cube_land)

CASE ( 'alb_tile_2' )
  cube_land = cube_from_array(alb_surft(:,:,2))
  cube = map_from_land(cube_land)

CASE ( 'alb_tile_3' )
  cube_land = cube_from_array(alb_surft(:,:,3))
  cube = map_from_land(cube_land)

CASE ( 'alb_tile_4' )
  cube_land = cube_from_array(alb_surft(:,:,4))
  cube = map_from_land(cube_land)

CASE ( 'anthrop_heat' )
  cube_land = cube_from_array(anthrop_heat_surft)
  cube = map_from_land(cube_land)

CASE ( 'canopy' )
  cube_land = cube_from_array(progs%canopy_surft)
  cube = map_from_land(cube_land)

CASE ( 'catch' )
  cube_land = cube_from_array(psparms%catch_surft)
  cube = map_from_land(cube_land)

CASE ( 'ecan' )
  cube_land = cube_from_array(ecan_surft)
  cube = map_from_land(cube_land)

CASE ( 'ei' )
  cube_land = cube_from_array(ei_surft)
  cube = map_from_land(cube_land)

CASE ( 'emis' )
  cube_land = cube_from_array(emis_surft)
  cube = map_from_land(cube_land)

CASE ( 'esoil' )
  cube_land = cube_from_array(esoil_surft)
  cube = map_from_land(cube_land)

CASE ( 'fqw' )
  ! Note that fqw_surft does not (always) do this job!
  DO n = 1,nsurft
    DO j = 1,surft_pts(n)
      l = ainfo%surft_index(j,n)
      workspace_surft(l,n) = ecan_surft(l,n) + ei_surft(l,n) +                &
                             esoil_surft(l,n)
      ! Add lake evaporation
      IF ( n == lake )                                                        &
        workspace_surft(l,n) = workspace_surft(l,n) + fqw_surft(l,n)
    END DO
  END DO
  cube_land = cube_from_array(workspace_surft)
  cube = map_from_land(cube_land)

CASE ( 'ftl' )
  cube_land = cube_from_array(ftl_surft)
  cube = map_from_land(cube_land)

CASE ( 'tau' )
  cube_land = cube_from_array(sf_diag%tau_surft)
  cube = map_from_land(cube_land)

CASE ( 'gc' )
  cube_land = cube_from_array(progs%gc_surft)
  cube = map_from_land(cube_land)

CASE ( 'le' )
  cube_land = cube_from_array(le_surft)
  cube = map_from_land(cube_land)

CASE ( 'nsnow' )
  cube_land = cube_from_array(REAL(progs%nsnow_surft))
  cube = map_from_land(cube_land)

CASE ( 'q1p5m' )
  cube_land = cube_from_array(sf_diag%q1p5m_surft)
  cube = map_from_land(cube_land)

CASE ( 'rad_net_tile' )
  cube_land = cube_from_array(radnet_surft)
  cube = map_from_land(cube_land)

CASE ( 'sw_surft' )
  cube_land = cube_from_array(sw_surft)
  cube = map_from_land(cube_land)

CASE ( 'lw_up_surft' )
  cube_land = cube_from_array(sf_diag%lw_up_surft)
  cube = map_from_land(cube_land)

CASE ( 'lw_down_surft' )
  cube_land = cube_from_array(sf_diag%lw_down_surft)
  cube = map_from_land(cube_land)

CASE ( 'rgrain' )
  cube_land = cube_from_array(progs%rgrain_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_can_melt' )
  ! Only include tiles where canopy snow model is used
  DO i = 1,nsurft
    IF ( canSnowTile(i) ) workspace_surft(:,i) = melt_surft(:,i)
  END DO
  cube_land = cube_from_array(workspace_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_can' )
  ! Only include tiles where canopy snow model is used
  DO i = 1,nsurft
    IF ( canSnowTile(i) ) workspace_surft(:,i) = progs%snow_surft(:,i)
  END DO
  cube_land = cube_from_array(workspace_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_depth' )
  cube_land = cube_from_array(progs%snowdepth_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_grnd_rho' )
  cube_land = cube_from_array(progs%rho_snow_grnd_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_grnd' )
  ! Only include tiles where canopy snow model is used
  DO i = 1,nsurft
    IF ( canSnowTile(i) ) workspace_surft(:,i) = progs%snow_grnd_surft(:,i)
  END DO
  cube_land = cube_from_array(workspace_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_ground' )
  DO i = 1,nsurft
    IF ( canSnowTile(i) ) THEN
      workspace_surft(:,i) = progs%snow_grnd_surft(:,i)
    ELSE
      workspace_surft(:,i) = progs%snow_surft(:,i)
    END IF
  END DO
  cube_land = cube_from_array(workspace_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_ice_tile' )
  DO n = 1,nsurft
    DO j = 1,surft_pts(n)
      i = ainfo%surft_index(j,n)
      workspace_surft(i,n) = SUM(progs%sice_surft(i,n,1:progs%nsnow_surft(i,n)))
    END DO
  END DO
  cube_land = cube_from_array(workspace_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_liq_tile' )
  DO n = 1,nsurft
    DO j = 1,surft_pts(n)
      i = ainfo%surft_index(j,n)
      workspace_surft(i,n) = SUM(progs%sliq_surft(i,n,1:progs%nsnow_surft(i,n)))
    END DO
  END DO
  cube_land = cube_from_array(workspace_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_mass' )
  workspace_surft(:,:) = progs%snow_surft(:,:)
  ! Add snow below canopy
  DO n = 1,nsurft
    IF ( canSnowTile(n) )                                                     &
      workspace_surft(:,n) = workspace_surft(:,n) + progs%snow_grnd_surft(:,n)
  END DO
  cube_land = cube_from_array(workspace_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_melt' )
  workspace_surft(:,:) = melt_surft(:,:)
  ! Add melting of snow below canopy.
  !     IF ( can_model == 4 )
  !workspace_surft(:,:) = workspace_surft(:,:) + snowGMeltDiag(:,:)
  cube_land = cube_from_array(workspace_surft)
  cube = map_from_land(cube_land)

CASE ( 'surf_ht_flux' )
  cube_land = cube_from_array(surf_htf_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_soil_htf' )
  cube_land = cube_from_array(snow_soil_htf)
  cube = map_from_land(cube_land)

CASE ( 'snice_smb_surft' )
  cube_land = cube_from_array(sf_diag%snice_smb_surft)
  cube = map_from_land(cube_land)

CASE ( 'snice_m_surft' )
  cube_land = cube_from_array(sf_diag%snice_m_surft)
  cube = map_from_land(cube_land)

CASE ( 'snice_freez_surft' )
  cube_land = cube_from_array(sf_diag%snice_freez_surft)
  cube = map_from_land(cube_land)

CASE ( 'snice_sicerate_surft' )
  cube_land = cube_from_array(sf_diag%snice_sicerate_surft)
  cube = map_from_land(cube_land)

CASE ( 'snice_sliqrate_surft' )
  cube_land = cube_from_array(sf_diag%snice_sliqrate_surft)
  cube = map_from_land(cube_land)

CASE ( 'snice_runoff_surft' )
  cube_land = cube_from_array(sf_diag%snice_runoff_surft)
  cube = map_from_land(cube_land)

CASE ( 'surf_ht_store' )
  cube_land = cube_from_array(surf_ht_store_surft)
  cube = map_from_land(cube_land)

CASE ( 't1p5m' )
  cube_land = cube_from_array(sf_diag%t1p5m_surft)
  cube = map_from_land(cube_land)

CASE ( 'tstar' )
  cube_land = cube_from_array(progs%tstar_surft)
  cube = map_from_land(cube_land)

CASE ( 'tsurf_elev_surft' )
  cube_land = cube_from_array(progs%tsurf_elev_surft)
  cube = map_from_land(cube_land)

CASE ( 'z0' )
  IF ( ANY(l_vegdrag_surft) ) THEN
    cube_land = cube_from_array(z0m_surft)
  ELSE
    cube_land = cube_from_array(psparms%z0_surft)
  END IF
  cube = map_from_land(cube_land)

CASE ( 'z0h' )
  !     This diagnostic is set from psparms%z0_surft unless separately
  !     aggregated.
  IF ( l_aggregate .AND. (i_aggregate_opt == 1) ) THEN
    cube_land = cube_from_array(psparms%z0h_bare_surft)
  ELSE
    IF ( ANY(l_vegdrag_surft) ) THEN
      cube_land = cube_from_array(z0h_surft)
    ELSE
      DO n = 1,nsurft
        workspace_surft(:,n) =  z0h_z0m(n) * psparms%z0_surft(:,n)
      END DO
      cube_land = cube_from_array(workspace_surft)
    END IF
  END IF
  cube = map_from_land(cube_land)

CASE ( 'tile_index' )
  cube_land = cube_from_array(REAL(ainfo%surft_index))
  cube = map_from_land(cube_land)

CASE ( 'isoprene_gb' )
  cube_land = cube_from_array(isoprene_gb)
  cube = map_from_land(cube_land)

CASE ( 'isoprene' )
  cube_land = cube_from_array(isoprene_pft)
  cube = map_from_land(cube_land)

CASE ( 'terpene_gb' )
  cube_land = cube_from_array(terpene_gb)
  cube = map_from_land(cube_land)

CASE ( 'terpene' )
  cube_land = cube_from_array(terpene_pft)
  cube = map_from_land(cube_land)

CASE ( 'methanol_gb' )
  cube_land = cube_from_array(methanol_gb)
  cube = map_from_land(cube_land)

CASE ( 'methanol' )
  cube_land = cube_from_array(methanol_pft)
  cube = map_from_land(cube_land)

CASE ( 'acetone_gb' )
  cube_land = cube_from_array(acetone_gb)
  cube = map_from_land(cube_land)

CASE ( 'acetone' )
  cube_land = cube_from_array(acetone_pft)
  cube = map_from_land(cube_land)

CASE ( 'croprootc' )
  cube_land = cube_from_array(crop_vars%rootc_cpft)
  cube = map_from_land(cube_land)

CASE ( 'cropharvc' )
  cube_land = cube_from_array(crop_vars%harvc_cpft)
  cube = map_from_land(cube_land)

CASE ( 'cropreservec' )
  cube_land = cube_from_array(crop_vars%reservec_cpft)
  cube = map_from_land(cube_land)

CASE ( 'cropdvi' )
  cube_land = cube_from_array(crop_vars%dvi_cpft)
  cube = map_from_land(cube_land)

CASE ( 'cropyield' )
  cube_land = cube_from_array(crop_vars%yield_diag_cpft)
  cube = map_from_land(cube_land)

CASE ( 'harvest_trigger' )
  workspace_cpft(:,:) = REAL(crop_vars%harvest_trigger_cpft(:,:))
  cube_land = cube_from_array(workspace_cpft)
  cube = map_from_land(cube_land)

CASE ( 'harvest_counter' )
  workspace_cpft(:,:) = REAL(crop_vars%harvest_counter_cpft(:,:))
  cube_land = cube_from_array(workspace_cpft)
  cube = map_from_land(cube_land)

CASE ( 'cropstemc' )
  cube_land = cube_from_array(crop_vars%stemc_diag_cpft)
  cube = map_from_land(cube_land)

CASE ( 'cropleafc' )
  cube_land = cube_from_array(crop_vars%leafc_diag_cpft)
  cube = map_from_land(cube_land)

CASE ( 'croplai' )
  cube_land = cube_from_array(crop_vars%croplai_cpft)
  cube = map_from_land(cube_land)

CASE ( 'cropcanht' )
  cube_land = cube_from_array(crop_vars%cropcanht_cpft)
  cube = map_from_land(cube_land)

CASE ( 'cropsowdate' )
  cube_land = cube_from_array(crop_vars%sow_date_cpft)
  cube = map_from_land(cube_land)

CASE ( 'croplatestharvdate' )
  cube_land = cube_from_array(crop_vars%latestharv_date_cpft)
  cube = map_from_land(cube_land)

CASE ( 'cropttveg' )
  cube_land = cube_from_array(crop_vars%tt_veg_cpft)
  cube = map_from_land(cube_land)

CASE ( 'cropttrep' )
  cube_land = cube_from_array(crop_vars%tt_rep_cpft)
  cube = map_from_land(cube_land)

CASE ( 'frac_irrig' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array(crop_vars%frac_irr_soilt(:,1))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output frac_irrig when nsoilt > 1. " //            &
                   "Use frac_irrig_soilt instead.")
  END IF

CASE ( 'frac_irrig_soilt' )
  cube_land = cube_from_array(crop_vars%frac_irr_soilt)
  cube = map_from_land(cube_land)

CASE ( 'irrfrac_irrtiles' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array(crop_vars%irrfrac_irrtiles(:,1))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
         "Cannot output irrfrac_irrtiles when nsoilt > 1.")
  END IF

CASE ( 'sthu_irr' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array(crop_vars%sthu_irr_soilt(:,1,:))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output sthu_irr when nsoilt > 1. " //              &
                   "Use sthu_irr_soilt instead")
  END IF

CASE ( 'sthu_irr_soilt' )
  cube_land = cube_from_array(crop_vars%sthu_irr_soilt)
  cube = map_from_land(cube_land)

CASE ( 'irrDaysDiag' )
  cube_land = cube_from_array(crop_vars%irrDaysDiag_gb)
  cube = map_from_land(cube_land)

CASE ( 'irrig_water' )
  cube_land = cube_from_array(crop_vars%irrig_water_gb)
  cube = map_from_land(cube_land)

CASE ( 'rgrainl' )
  cube_land = cube_from_array(progs%rgrainl_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_ds' )
  cube_land = cube_from_array(progs%ds_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_ice' )
  cube_land = cube_from_array(progs%sice_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_liq' )
  cube_land = cube_from_array(progs%sliq_surft)
  cube = map_from_land(cube_land)

CASE ( 'tsnow' )
  cube_land = cube_from_array(progs%tsnow_surft)
  cube = map_from_land(cube_land)

CASE ( 'rflow' )
  cube_land = cube_from_array(rflow_gb)
  cube = map_from_land(cube_land)

CASE ( 'rrun' )
  cube_land = cube_from_array(rrun_gb)
  cube = map_from_land(cube_land)

CASE ( 'frac_fplain_lp' )
  cube_land = cube_from_array(frac_fplain_lp)
  cube = map_from_land(cube_land)

CASE ( 'lake_depth' )
  cube_land = cube_from_array(lake_depth_gb)
  cube = map_from_land(cube_land)

CASE ( 'lake_fetch_gb' )
  cube_land = cube_from_array(lake_fetch_gb)
  cube = map_from_land(cube_land)

CASE ( 'lake_t_mean_gb' )
  cube_land = cube_from_array(lake_t_mean_gb)
  cube = map_from_land(cube_land)

CASE ( 'lake_t_mxl_gb' )
  cube_land = cube_from_array(lake_t_mxl_gb)
  cube = map_from_land(cube_land)

CASE ( 'lake_t_ice_gb' )
  cube_land = cube_from_array(lake_t_ice_gb)
  cube = map_from_land(cube_land)

CASE ( 'lake_h_mxl_gb' )
  cube_land = cube_from_array(lake_h_mxl_gb)
  cube = map_from_land(cube_land)

CASE ( 'lake_h_ice_gb' )
  cube_land = cube_from_array(lake_h_ice_gb)
  cube = map_from_land(cube_land)

CASE ( 'lake_shape_factor_gb' )
  cube_land = cube_from_array(lake_shape_factor_gb)
  cube = map_from_land(cube_land)

CASE ( 'g_dt_gb' )
  cube_land = cube_from_array(g_dt_gb)
  cube = map_from_land(cube_land)

CASE ( 'lake_t_sfc_gb' )
  cube_land = cube_from_array(lake_t_sfc_gb)
  cube = map_from_land(cube_land)

CASE ( 'lake_t_snow_gb' )
  cube_land = cube_from_array(lake_t_snow_gb)
  cube = map_from_land(cube_land)

CASE ( 'lake_h_snow_gb' )
  cube_land = cube_from_array(lake_h_snow_gb)
  cube = map_from_land(cube_land)

CASE ( 'lake_albedo_gb' )
  cube_land = cube_from_array(lake_albedo_gb)
  cube = map_from_land(cube_land)

  !
  !-----------------------------------------------------------------------------
  ! Equivalent neutral winds
  !-----------------------------------------------------------------------------
CASE ( 'u10m_n' )
  cube = cube_from_array(sf_diag%u10m_n)

CASE ( 'v10m_n' )
  cube = cube_from_array(sf_diag%v10m_n)

CASE ( 'mu10m_n' )
  cube = cube_from_array(sf_diag%mu10m_n)

CASE ( 'mv10m_n' )
  cube = cube_from_array(sf_diag%mv10m_n)

  ! Fire module variables
CASE ('fire_mcarthur')
  cube_land = cube_from_array(fire_diag(:)%mcarthur%ffdi)
  cube      = map_from_land(cube_land)

CASE ('fire_canadian_ffmc')
  cube_land = cube_from_array(fire_prog(:)%canadian%ffmc)
  cube      = map_from_land(cube_land)

CASE ('fire_canadian_dmc')
  cube_land = cube_from_array(fire_prog(:)%canadian%dmc)
  cube      = map_from_land(cube_land)

CASE ('fire_canadian_dc')
  cube_land = cube_from_array(fire_prog(:)%canadian%dc)
  cube      = map_from_land(cube_land)

CASE ('fire_canadian_isi')
  cube_land = cube_from_array(fire_diag(:)%canadian%isi)
  cube      = map_from_land(cube_land)

CASE ('fire_canadian_bui')
  cube_land = cube_from_array(fire_diag(:)%canadian%bui)
  cube      = map_from_land(cube_land)

CASE ('fire_canadian')
  cube_land = cube_from_array(fire_diag(:)%canadian%fwi)
  cube      = map_from_land(cube_land)

CASE ('fire_nesterov')
  cube_land = cube_from_array(fire_prog(:)%nesterov%findex)
  cube      = map_from_land(cube_land)

  ! INFERNO variables
CASE ( 'flammability' )
  cube_land = cube_from_array(fire_vars%flammability_ft)
  cube = map_from_land(cube_land)

CASE ( 'burnt_area_gb' )
  cube_land = cube_from_array(fire_vars%burnt_area)
  cube = map_from_land(cube_land)

CASE ( 'burnt_area' )
  cube_land = cube_from_array(fire_vars%burnt_area_ft)
  cube = map_from_land(cube_land)

CASE ( 'emitted_carbon_gb' )
  cube_land = cube_from_array(fire_vars%emitted_carbon)
  cube = map_from_land(cube_land)
CASE ( 'emitted_carbon' )
  cube_land = cube_from_array(fire_vars%emitted_carbon_ft)
  cube = map_from_land(cube_land)
CASE ( 'emitted_carbon_DPM' )
  cube_land = cube_from_array(fire_vars%emitted_carbon_DPM)
  cube = map_from_land(cube_land)
CASE ( 'emitted_carbon_RPM' )
  cube_land = cube_from_array(fire_vars%emitted_carbon_RPM)
  cube = map_from_land(cube_land)

CASE ( 'fire_em_CO2_gb' )
  cube_land = cube_from_array(fire_vars%fire_em_CO2)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_CO2' )
  cube_land = cube_from_array(fire_vars%fire_em_CO2_ft)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_CO2_DPM' )
  cube_land = cube_from_array(fire_vars%fire_em_CO2_DPM)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_CO2_RPM' )
  cube_land = cube_from_array(fire_vars%fire_em_CO2_RPM)
  cube = map_from_land(cube_land)

CASE ( 'fire_em_CO_gb' )
  cube_land = cube_from_array(fire_vars%fire_em_CO)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_CO' )
  cube_land = cube_from_array(fire_vars%fire_em_CO_ft)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_CO_DPM' )
  cube_land = cube_from_array(fire_vars%fire_em_CO_DPM)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_CO_RPM' )
  cube_land = cube_from_array(fire_vars%fire_em_CO_RPM)
  cube = map_from_land(cube_land)

CASE ( 'fire_em_CH4_gb' )
  cube_land = cube_from_array(fire_vars%fire_em_CH4)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_CH4' )
  cube_land = cube_from_array(fire_vars%fire_em_CH4_ft)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_CH4_DPM' )
  cube_land = cube_from_array(fire_vars%fire_em_CH4_DPM)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_CH4_RPM' )
  cube_land = cube_from_array(fire_vars%fire_em_CH4_RPM)
  cube = map_from_land(cube_land)

CASE ( 'fire_em_NOx_gb' )
  cube_land = cube_from_array(fire_vars%fire_em_NOx)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_NOx' )
  cube_land = cube_from_array(fire_vars%fire_em_NOx_ft)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_NOx_DPM' )
  cube_land = cube_from_array(fire_vars%fire_em_NOx_DPM)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_NOx_RPM' )
  cube_land = cube_from_array(fire_vars%fire_em_NOx_RPM)
  cube = map_from_land(cube_land)

CASE ( 'fire_em_SO2_gb' )
  cube_land = cube_from_array(fire_vars%fire_em_SO2)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_SO2' )
  cube_land = cube_from_array(fire_vars%fire_em_SO2_ft)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_SO2_DPM' )
  cube_land = cube_from_array(fire_vars%fire_em_SO2_DPM)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_SO2_RPM' )
  cube_land = cube_from_array(fire_vars%fire_em_SO2_RPM)
  cube = map_from_land(cube_land)

CASE ( 'fire_em_OC_gb' )
  cube_land = cube_from_array(fire_vars%fire_em_OC)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_OC' )
  cube_land = cube_from_array(fire_vars%fire_em_OC_ft)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_OC_DPM' )
  cube_land = cube_from_array(fire_vars%fire_em_OC_DPM)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_OC_RPM' )
  cube_land = cube_from_array(fire_vars%fire_em_OC_RPM)
  cube = map_from_land(cube_land)

CASE ( 'fire_em_BC_gb' )
  cube_land = cube_from_array(fire_vars%fire_em_BC)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_BC' )
  cube_land = cube_from_array(fire_vars%fire_em_BC_ft)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_BC_DPM' )
  cube_land = cube_from_array(fire_vars%fire_em_BC_DPM)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_BC_RPM' )
  cube_land = cube_from_array(fire_vars%fire_em_BC_RPM)
  cube = map_from_land(cube_land)

  !-----------------------------------------------------------------------------
  ! Variables available with any multi-pool soil C model (e.g. RothC, ECOSSE).
  !-----------------------------------------------------------------------------
  ! Soil prognostic variables and single pools from multi-pool prognostic
  ! variables. These are all coded for a single soil tile only.
  !-----------------------------------------------------------------------------
CASE ( 'c_bio' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array(progs%cs_pool_soilt(:,1,:,3))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output c_bio when nsoilt > 1.")
  END IF

CASE ( 'c_dpm' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array(progs%cs_pool_soilt(:,1,:,1))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output c_dpm when nsoilt > 1.")
  END IF

CASE ( 'c_hum' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array(progs%cs_pool_soilt(:,1,:,4))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output c_hum when nsoilt > 1.")
  END IF

CASE ( 'c_rpm' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array(progs%cs_pool_soilt(:,1,:,2))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output c_rpm when nsoilt > 1.")
  END IF

CASE ( 'n_bio' )
  cube_land = cube_from_array(progs%ns_pool_gb(:,:,3))
  cube = map_from_land(cube_land)

CASE ( 'n_dpm' )
  cube_land = cube_from_array(progs%ns_pool_gb(:,:,1))
  cube = map_from_land(cube_land)

CASE ( 'n_hum' )
  cube_land = cube_from_array(progs%ns_pool_gb(:,:,4))
  cube = map_from_land(cube_land)

CASE ( 'n_rpm' )
  cube_land = cube_from_array(progs%ns_pool_gb(:,:,2))
  cube = map_from_land(cube_land)

  ! Gridbox totals of multi-pool soil prognostic variables.
  ! All are coded for a single soil tile only.
CASE ( 'c_bio_gb' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array( SUM(progs%cs_pool_soilt(:,1,:,3),2) )
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output c_bio_gb when nsoilt > 1.")
  END IF

CASE ( 'c_dpm_gb' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array( SUM(progs%cs_pool_soilt(:,1,:,1),2) )
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output c_dpm_gb when nsoilt > 1.")
  END IF

CASE ( 'c_hum_gb' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array( SUM(progs%cs_pool_soilt(:,1,:,4),2) )
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output c_hum_gb when nsoilt > 1.")
  END IF

CASE ( 'c_rpm_gb' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array( SUM(progs%cs_pool_soilt(:,1,:,2),2) )
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output c_rpm_gb when nsoilt > 1.")
  END IF

CASE ( 'n_bio_gb' )
  cube_land = cube_from_array( SUM(progs%ns_pool_gb(:,:,3),2) )
  cube = map_from_land(cube_land)

CASE ( 'n_dpm_gb' )
  cube_land = cube_from_array( SUM(progs%ns_pool_gb(:,:,1),2) )
  cube = map_from_land(cube_land)

CASE ( 'n_hum_gb' )
  cube_land = cube_from_array( SUM(progs%ns_pool_gb(:,:,4),2) )
  cube = map_from_land(cube_land)

CASE ( 'n_rpm_gb' )
  cube_land = cube_from_array( SUM(progs%ns_pool_gb(:,:,2),2) )
  cube = map_from_land(cube_land)

  ! Gridbox totals of derived variables for multi-pool soil models.
CASE ( 'n_soil_gb' )
  ! Coded for a single soil tile only.
  IF (nsoilt == 1) THEN
    IF ( soil_bgc_model == soil_model_rothc ) THEN
      ! Sum organic pools.
      workspace_land(:) = SUM( SUM(progs%ns_pool_gb(:,:,:),3), 2 )
      ! Add inorganic pools.
      workspace_land(:) = workspace_land(:) +                                 &
                          SUM( progs%n_inorg_soilt_lyrs(:,1,:), 2 )
    ELSE IF ( soil_bgc_model == soil_model_ecosse ) THEN
      ! Sum over all pools and levels.
      workspace_land(:) = SUM( SUM(                                           &
                        soilecosse%n_soil_pool_soilt(:,1,:,:), 3), 2 )
    END IF
    ! Map the result onto the full grid.
    cube_land = cube_from_array(workspace_land)
    cube      = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output n_soil_gb when nsoilt > 1.")
  END IF

  ! ECOSSE prognostic variables (and components).
  ! All are coded for a single soil tile only.
CASE ( 'n_soil' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soilecosse%n_soil_pool_soilt(:,1,:,:))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output n_soil when nsoilt > 1")
  END IF

CASE ( 'n_amm' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array(soilecosse%n_soil_pool_soilt(:,1,:,i_amm))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output n_amm when nsoilt > 1.")
  END IF

CASE ( 'n_nit' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array(soilecosse%n_soil_pool_soilt(:,1,:,i_nit))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output n_nit when nsoilt > 1.")
  END IF

  ! Gridbox totals of components of ECOSSE prognostic variables.
CASE ( 'n_amm_gb' )
  IF (nsoilt == 1) THEN
    ! Sum over levels.
    cube_land = cube_from_array( SUM(                                         &
                                soilecosse%n_soil_pool_soilt(:,1,:,i_amm), 2 ))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output n_amm_gb when nsoilt > 1.")
  END IF

CASE ( 'n_nit_gb' )
  IF (nsoilt == 1) THEN
    ! Sum over levels.
    cube_land = cube_from_array( SUM(                                         &
                                soilecosse%n_soil_pool_soilt(:,1,:,i_nit), 2 ))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output n_nit_gb when nsoilt > 1.")
  END IF

  !   ECOSSE fluxes.
CASE ( 'plant_input_c_gb' )
  cube_land = cube_from_array( soilecosse%plant_input_c_gb(:) )
  cube = map_from_land(cube_land)

CASE ( 'plant_input_n_gb' )
  cube_land = cube_from_array( soilecosse%plant_input_n_gb(:) )
  cube = map_from_land(cube_land)

  ! Imogen prognostics
CASE ( 'd_land_atmos_co2' )
  workspace_land(:) = d_land_atmos_co2
  cube_land = cube_from_array(workspace_land)
  cube = map_from_land(cube_land)

CASE ( 'd_ocean_atmos' )
  workspace_land(:) = d_ocean_atmos
  cube_land = cube_from_array(workspace_land)
  cube = map_from_land(cube_land)

CASE ( 'c_emiss_out' )
  workspace_land(:) = c_emiss_out
  cube_land = cube_from_array(workspace_land)
  cube = map_from_land(cube_land)
  ! Imogen end

CASE ( 'n_denitrif_gb' )
  cube_land = cube_from_array( soilecosse%n_denitrification_gb )
  cube = map_from_land(cube_land)

CASE ( 'n_nitrif_gb' )
  cube_land = cube_from_array( soilecosse%n_nitrification_gb )
  cube = map_from_land(cube_land)

CASE ( 'no_soil_gb' )
  cube_land = cube_from_array( soilecosse%no_soil_gb )
  cube = map_from_land(cube_land)

CASE ( 'n2_denitrif_gb' )
  cube_land = cube_from_array( soilecosse%n2_denitrif_gb )
  cube = map_from_land(cube_land)

CASE ( 'n2o_denitrif_gb' )
  cube_land = cube_from_array( soilecosse%n2o_denitrif_gb )
  cube = map_from_land(cube_land)

CASE ( 'n2o_nitrif_gb' )
  cube_land = cube_from_array( soilecosse%n2o_nitrif_gb )
  cube = map_from_land(cube_land)

CASE ( 'n2o_part_nitrif_gb' )
  cube_land = cube_from_array( soilecosse%n2o_partial_nitrif_gb )
  cube = map_from_land(cube_land)

CASE ( 'n2o_soil_gb' )
  cube_land = cube_from_array( soilecosse%n2o_soil_gb )
  cube = map_from_land(cube_land)

  ! ECOSSE leaching diagnostics.
CASE ( 'n_leach_amm_gb' )
  cube_land = cube_from_array( soilecosse%n_leach_amm_gb )
  cube = map_from_land(cube_land)

CASE ( 'n_leach_nit_gb' )
  cube_land = cube_from_array( soilecosse%n_leach_nit_gb )
  cube = map_from_land(cube_land)

  ! Atmospheric deposition variables.
CASE ( 'tracer_field' )
  cube_land = cube_from_array( tracer_field )
  cube = map_from_land(cube_land)

CASE ( 't_growth_gb' )
  cube_land = cube_from_array( t_growth_gb )
  cube = map_from_land(cube_land)

  ! Water resource variables.
CASE ( 'conveyance_loss' )
  cube_land = cube_from_array( conveyance_loss )
  cube = map_from_land(cube_land)

CASE ( 'demand_domestic' )
  cube_land = cube_from_array(demand_accum(:,use_domestic))
  cube = map_from_land(cube_land)

CASE ( 'demand_environment' )
  cube_land = cube_from_array(demand_accum(:,use_environment))
  cube = map_from_land(cube_land)

CASE ( 'demand_industry' )
  cube_land = cube_from_array(demand_accum(:,use_industry))
  cube = map_from_land(cube_land)

CASE ( 'demand_irrigation' )
  cube_land = cube_from_array(demand_accum(:,use_irrigation))
  cube = map_from_land(cube_land)

CASE ( 'demand_livestock' )
  cube_land = cube_from_array(demand_accum(:,use_livestock))
  cube = map_from_land(cube_land)

CASE ( 'demand_transfers' )
  cube_land = cube_from_array(demand_accum(:,use_transfers))
  cube = map_from_land(cube_land)

CASE ( 'irrig_eff' )
  cube_land = cube_from_array( irrig_eff )
  cube = map_from_land(cube_land)

CASE ( 'grid_area' )
  cube = cube_from_array( grid_area_ij )

CASE ( 'sfc_water_frac' )
  cube_land = cube_from_array( sfc_water_frac )
  cube = map_from_land(cube_land)

CASE ( 'water_demand' )
  ! Sum over active sectors.
  workspace_land(:) = 0.0
  IF ( l_water_domestic ) THEN
    workspace_land(:) =  workspace_land(:) + demand_accum(:,use_domestic)
  END IF
  IF ( l_water_environment ) THEN
    workspace_land(:) =  workspace_land(:) + demand_accum(:,use_environment)
  END IF
  IF ( l_water_industry ) THEN
    workspace_land(:) =  workspace_land(:) + demand_accum(:,use_industry)
  END IF
  IF ( l_water_irrigation ) THEN
    workspace_land(:) =  workspace_land(:) + demand_accum(:,use_irrigation)
  END IF
  IF ( l_water_livestock ) THEN
    workspace_land(:) =  workspace_land(:) + demand_accum(:,use_livestock)
  END IF
  IF ( l_water_transfers ) THEN
    workspace_land(:) =  workspace_land(:) + demand_accum(:,use_transfers)
  END IF
  cube_land = cube_from_array( workspace_land(:) )
  cube = map_from_land(cube_land)

CASE DEFAULT
  CALL log_fatal("extract_var",                                               &
                 "Unrecognised variable for output - '" //                    &
                 TRIM(get_string_identifier(var_id)) // "'. " //              &
                 "See docs for available variables")
END SELECT

!-----------------------------------------------------------------------------
! Free the land cube
!
! Note that is safe to call this routine even if the cube has not been
! allocated - it just won't do anything
!-----------------------------------------------------------------------------
CALL cube_free(cube_land)


RETURN

END FUNCTION extract_var
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

SUBROUTINE populate_var(var_id, cube, const_val)

USE jules_fields_mod, ONLY: crop_vars, psparms, toppdm, fire_vars, ainfo,     &
                            trif_vars, soilecosse, urban_param,progs, trifctltype, &
                            jules_vars

  !Science variables
! TYPE Definitions
USE jules_fields_mod, ONLY: toppdm

USE ancil_info, ONLY:                                                         &
  dim_cslayer, dim_soil_n_pool, land_pts, nsurft, nsoilt

USE model_grid_mod, ONLY:                                                     &
  grid_area_ij, latitude, longitude

USE coastal, ONLY:                                                            &
  flandg

USE fluxes, ONLY:                                                             &
  t_growth_gb

USE lake_mod, ONLY:                                                           &
  lake_depth_gb, lake_fetch_gb, lake_h_mxl_gb, lake_shape_factor_gb,          &
  lake_t_mean_gb, lake_t_mxl_gb, lake_t_ice_gb, lake_h_ice_gb,                &
  lake_t_sfc_gb

USE prognostics, ONLY:                                                        &
  l_broadcast_soilt

USE forcing, ONLY:                                                            &
  pstar_ij, qw_1_ij, tl_1_ij, lw_down_ij, sw_down_ij,                         &
  diff_rad_ij, ls_rain_ij, ls_snow_ij, con_rain_ij,                           &
  con_snow_ij, u_1_ij, v_1_ij, diurnal_temperature_range_ij

USE ozone_vars, ONLY:                                                         &
  o3_gb

USE jules_deposition_mod, ONLY:                                               &
  tracer_field

USE jules_soil_mod, ONLY:                                                     &
  sm_levels, l_tile_soil, l_broadcast_ancils

USE jules_rivers_mod, ONLY:                                                   &
  rivers_dir, rivers_seq, rivers_dra,                                         &
  rivers_xgrid, rivers_ygrid,                                                 &
  rivers_lat2d, rivers_lon2d,                                                 &
  rivers_sto_rp, rfm_surfstore_rp,                                            &
  rfm_substore_rp, rfm_flowin_rp, rfm_bflowin_rp

USE jules_water_resources_mod, ONLY:                                          &
  conveyance_loss, demand_rate_domestic, demand_rate_industry,                &
  demand_rate_livestock, demand_rate_transfers, irrig_eff, sfc_water_frac

USE overbank_inundation_mod, ONLY:                                            &
  logn_mean, logn_stdev

USE aero, ONLY:                                                               &
  co2_mmr

USE fire_mod, ONLY:                                                           &
  fire_prog

USE metstats_mod, ONLY:                                                       &
  metstats_prog

!Others
USE data_cube_mod, ONLY:                                                      &
  data_cube, cube_get_data, cube_free, cube_create

USE jules_surface_types_mod, ONLY: ncpft, nnpft

USE cable_prognostic_info_mod, ONLY: SoilTemp_CABLE,  SoilMoisture_CABLE,     &
                          FrozenSoilFrac_CABLE, SnowDepth_CABLE,              &
                          SnowMass_CABLE,  SnowTemp_CABLE,  SnowDensity_CABLE,&
                          ThreeLayerSnowFlag_CABLE, OneLyrSnowDensity_CABLE,  &
                          SnowAge_CABLE

IMPLICIT NONE

! Argument types
INTEGER, INTENT(IN) :: var_id  ! Identifies the variable to fill
TYPE(data_cube), INTENT(IN), OPTIONAL :: cube
                     ! The data to put in to the variable as a cube
REAL, INTENT(IN), OPTIONAL :: const_val
                     ! A constant value to fill all elements of the variable
                     ! with


! Work variables
TYPE(data_cube) :: cube_land  ! Workspace cube for land and river data
                              ! This is required so that it can be deallocated
                              ! to avoid memory leaks

REAL :: nsnow_real(land_pts, nsurft)  ! Real version of nsnow

REAL :: co2_tmp(land_pts)

REAL, ALLOCATABLE :: frac_cpft(:, :)

!-----------------------------------------------------------------------------


IF ( .NOT. PRESENT(cube) .AND. .NOT. PRESENT(const_val) )                     &
  CALL log_fatal("populate_var",                                              &
                 "Neither data cube or const_val have been provided for " //  &
                 "variable '" // TRIM(get_string_identifier(var_id)) // "'")

IF ( PRESENT(cube) .AND. PRESENT(const_val) )                                 &
  CALL log_warn("populate_var",                                               &
                "data cube and const_val both provided for variable '" //     &
                TRIM(get_string_identifier(var_id)) //                        &
                "' - using data in preference")

! We use the string identifier to select a CASE from the SELECT statement to
! avoid being dependent on how integer variable ids are implemented
SELECT CASE ( get_string_identifier(var_id) )

CASE ( 'latitude' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, latitude)
  ELSE
    latitude(:,:) = const_val
  END IF

CASE ( 'longitude' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, longitude)
  ELSE
    longitude(:,:) = const_val
  END IF

CASE ( 'land_fraction' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, flandg)
  ELSE
    flandg(:,:) = const_val
  END IF

CASE ( 'surf_hgt_surft' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, jules_vars%surf_hgt_surft)
  ELSE
    jules_vars%surf_hgt_surft(:,:) = const_val
  END IF

CASE ( 'z_land_land' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, jules_vars%z_land_land)
  ELSE
    jules_vars%z_land_land(:) = const_val
  END IF

CASE ( 'frac' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, ainfo%frac_surft)
  ELSE
    ainfo%frac_surft(:,:) = const_val
  END IF

CASE ( 'cropfrac' )
  ALLOCATE(frac_cpft(land_pts, ncpft))
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, frac_cpft)
  ELSE
    frac_cpft(:,:) = const_val
  END IF
  ainfo%frac_surft(:,1 + nnpft:ncpft + nnpft) = frac_cpft(:,:)
  DEALLOCATE(frac_cpft)

  !-----------------------------------------------------------------------------
  ! Soil properties have a corresponding 0-level version for input of constant
  ! z values. albsoil_gb, psparms%albobs_sw_gb, etc have no levels anyway
  !
  ! Each variable can be read in with _soilt if the soil tile dimension is
  ! present. Where appropriate, this can be combined with _const_z as well.
  ! There are a few possible outcomes
  !
  ! No suffix:
  !  l_tile_soil = T and l_broadcast_ancils = T
  !  --Spread around soilt dim
  !  l_tile_soil = F
  !  --Do nothing
  !
  ! _const_z suffix:
  !  l_tile_soil = T and l_broadcast_ancils = T (and l_const_z = T)
  !  --Spread around soilt and layer dims
  !  l_tile_soil = F
  !  --Spread around layer dim
  !
  ! _soilt suffix:
  ! l_tile_soil = T and l_broadcast_ancils = F
  ! --Do nothing
  !
  ! _const_z_soilt:
  ! l_tile_soil = T and l_broadcast_ancils = F (and l_const_z = T)
  ! --Spread around layer dim
  !
  !-----------------------------------------------------------------------------

      !=== Start cases for albsoil ===
CASE ( 'albsoil' )  !no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%albsoil_soilt(:,1))
      psparms%albsoil_soilt(:,:) = SPREAD(psparms%albsoil_soilt(:,1), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%albsoil_soilt(:,1))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate albsoil- failure in logic")
    END IF
  ELSE
    psparms%albsoil_soilt(:,:) = const_val
  END IF

CASE ( 'albsoil_soilt' ) !with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%albsoil_soilt(:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate psparms%albsoil_soilt- failure in logic")
    END IF
  ELSE
    psparms%albsoil_soilt(:,:) = const_val
  END IF
  !=== End cases for albsoil ===

  !=== Start cases for clay ===
CASE ( 'clay' )  ! With levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      ! Spread around soilt dim.
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%clay_soilt(:,1,1:dim_cslayer))
      psparms%clay_soilt(:,:,:) = SPREAD(psparms%clay_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      ! Do nothing special.
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%clay_soilt(:,1,1:dim_cslayer))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate clay - failure in logic")
    END IF
  ELSE
    psparms%clay_soilt(:,:,:) = const_val
  END IF

CASE ( 'clay_const_z' )  ! No levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      ! Spread around soilt and layer dims
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%clay_soilt(:,1,1))
      psparms%clay_soilt(:,1,1:dim_cslayer) =                                 &
                        SPREAD(psparms%clay_soilt(:,1,1), 2, dim_cslayer)
      psparms%clay_soilt(:,:,:) = SPREAD(psparms%clay_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      ! Copy the values from the first vertical level to all others
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%clay_soilt(:,1,1))
      psparms%clay_soilt(:,1,1:dim_cslayer) =                                 &
                        SPREAD(psparms%clay_soilt(:,1,1), 2, dim_cslayer)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate clay_const_z - failure in logic")
    END IF
  ELSE
    psparms%clay_soilt(:,:,:) = const_val
  END IF

CASE ( 'clay_soilt' )  ! With levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       ! Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%clay_soilt(:,:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate psparms%clay_soilt- failure in logic")
    END IF
  ELSE
    psparms%clay_soilt(:,:,:) = const_val
  END IF

CASE ( 'clay_const_z_soilt' ) ! No levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Spread around layer dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%clay_soilt(:,1:nsoilt,1))
      psparms%clay_soilt(:,:,:) =                                             &
                        SPREAD(psparms%clay_soilt(:,:,1), 3, dim_cslayer)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate clay_const_z_soilt- failure in logic")
    END IF
  ELSE
    psparms%clay_soilt(:,:,:) = const_val
  END IF
  !=== End cases for clay ===

  !=== Start cases for b ===
CASE ( 'b' )  !With levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%bexp_soilt(:,1,1:sm_levels))
      psparms%bexp_soilt(:,:,:) = SPREAD(psparms%bexp_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%bexp_soilt(:,1,1:sm_levels))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate b- failure in logic")
    END IF
  ELSE
    psparms%bexp_soilt(:,:,:) = const_val
  END IF

CASE ( 'b_const_z' ) !No levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      ! Spread around soilt and layer dims
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%bexp_soilt(:,1,1))
      psparms%bexp_soilt(:,1,1:sm_levels) =                                   &
                            SPREAD(psparms%bexp_soilt(:,1,1), 2, sm_levels)
      psparms%bexp_soilt(:,:,:)           =                                   &
                            SPREAD(psparms%bexp_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      ! Copy the values from the first vertical level to all others
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%bexp_soilt(:,1,1))
      psparms%bexp_soilt(:,1,1:sm_levels) =                                   &
                            SPREAD(psparms%bexp_soilt(:,1,1), 2, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate b_const_z- failure in logic")
    END IF
  ELSE
    psparms%bexp_soilt(:,:,:) = const_val
  END IF

CASE ( 'b_soilt' ) !With levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%bexp_soilt(:,:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate b_soilt- failure in logic")
    END IF
  ELSE
    psparms%bexp_soilt(:,:,:) = const_val
  END IF

CASE ( 'b_const_z_soilt' ) !No levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Spread around layer dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%bexp_soilt(:,1:nsoilt,1))
      psparms%bexp_soilt(:,:,:) =                                             &
                        SPREAD(psparms%bexp_soilt(:,:,1), 3, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate b_const_z_soilt- failure in logic")
    END IF
  ELSE
    psparms%bexp_soilt(:,:,:) = const_val
  END IF
  !=== End cases for b ===

  !=== Start cases for sathh ===
CASE ( 'sathh' )  !With levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%sathh_soilt(:,1,1:sm_levels))
      psparms%sathh_soilt(:,:,:) = SPREAD(psparms%sathh_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%sathh_soilt(:,1,1:sm_levels))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sathh- failure in logic")
    END IF
  ELSE
    psparms%sathh_soilt(:,:,:) = const_val
  END IF

CASE ( 'sathh_const_z' ) !No levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      ! Spread around soilt and layer dims
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%sathh_soilt(:,1,1))
      psparms%sathh_soilt(:,1,1:sm_levels) =                                  &
                            SPREAD(psparms%sathh_soilt(:,1,1), 2, sm_levels)
      psparms%sathh_soilt(:,:,:)           =                                  &
                            SPREAD(psparms%sathh_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      ! Copy the values from the first vertical level to all others
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%sathh_soilt(:,1,1))
      psparms%sathh_soilt(:,1,1:sm_levels) =                                  &
                            SPREAD(psparms%sathh_soilt(:,1,1), 2, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sathh_const_z- failure in logic")
    END IF
  ELSE
    psparms%sathh_soilt(:,:,:) = const_val
  END IF

CASE ( 'sathh_soilt' ) !With levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%sathh_soilt(:,:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate psparms%sathh_soilt- failure in logic")
    END IF
  ELSE
    psparms%sathh_soilt(:,:,:) = const_val
  END IF

CASE ( 'sathh_const_z_soilt' ) !No levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Spread around layer dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%sathh_soilt(:,1:nsoilt,1))
      psparms%sathh_soilt(:,:,:) =                                            &
                        SPREAD(psparms%sathh_soilt(:,:,1), 3, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sathh_const_z_soilt- failure in logic")
    END IF
  ELSE
    psparms%sathh_soilt(:,:,:) = const_val
  END IF
  !=== End cases for sathh ===

  !=== Start cases for sathh ===
CASE ( 'satcon' )  !With levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%satcon_soilt(:,1,1:sm_levels))
      psparms%satcon_soilt(:,:,:) =                                           &
                        SPREAD(psparms%satcon_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%satcon_soilt(:,1,1:sm_levels))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate satcon- failure in logic")
    END IF
  ELSE
    psparms%satcon_soilt(:,:,:) = const_val
  END IF

CASE ( 'satcon_const_z' ) !No levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      ! Spread around soilt and layer dims
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%satcon_soilt(:,1,1))
      psparms%satcon_soilt(:,1,1:sm_levels) =                                 &
                        SPREAD(psparms%satcon_soilt(:,1,1), 2, sm_levels)
      psparms%satcon_soilt(:,:,:)           =                                 &
                        SPREAD(psparms%satcon_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      ! Copy the values from the first vertical level to all others
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%satcon_soilt(:,1,1))
      psparms%satcon_soilt(:,1,1:sm_levels) =                                 &
                        SPREAD(psparms%satcon_soilt(:,1,1), 2, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate satcon_const_z- failure in logic")
    END IF
  ELSE
    psparms%satcon_soilt(:,:,:) = const_val
  END IF

CASE ( 'satcon_soilt' ) !With levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%satcon_soilt(:,:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate psparms%satcon_soilt- failure in logic")
    END IF
  ELSE
    psparms%satcon_soilt(:,:,:) = const_val
  END IF

CASE ( 'satcon_const_z_soilt' ) !No levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Spread around layer dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%satcon_soilt(:,1:nsoilt,1))
      psparms%satcon_soilt(:,:,:) =                                           &
                        SPREAD(psparms%satcon_soilt(:,:,1), 3, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate satcon_const_z_soilt- failure in logic")
    END IF
  ELSE
    psparms%satcon_soilt(:,:,:) = const_val
  END IF
  !=== End cases for satcon ===

  !=== Start cases for sm_sat ===
CASE ( 'sm_sat' )  !With levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%smvcst_soilt(:,1,1:sm_levels))
      psparms%smvcst_soilt(:,:,:) =                                           &
                        SPREAD(psparms%smvcst_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%smvcst_soilt(:,1,1:sm_levels))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_sat- failure in logic")
    END IF
  ELSE
    psparms%smvcst_soilt(:,:,:) = const_val
  END IF

CASE ( 'sm_sat_const_z' ) !No levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      ! Spread around soilt and layer dims
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%smvcst_soilt(:,1,1))
      psparms%smvcst_soilt(:,1,1:sm_levels) =                                 &
                        SPREAD(psparms%smvcst_soilt(:,1,1), 2, sm_levels)
      psparms%smvcst_soilt(:,:,:)           =                                 &
                        SPREAD(psparms%smvcst_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      ! Copy the values from the first vertical level to all others
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%smvcst_soilt(:,1,1))
      psparms%smvcst_soilt(:,1,1:sm_levels) =                                 &
                        SPREAD(psparms%smvcst_soilt(:,1,1), 2, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_sat_const_z- failure in logic")
    END IF
  ELSE
    psparms%smvcst_soilt(:,:,:) = const_val
  END IF

CASE ( 'sm_sat_soilt' ) !With levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%smvcst_soilt(:,:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_sat_soilt- failure in logic")
    END IF
  ELSE
    psparms%smvcst_soilt(:,:,:) = const_val
  END IF

CASE ( 'sm_sat_const_z_soilt' ) !No levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Spread around layer dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%smvcst_soilt(:,1:nsoilt,1))
      psparms%smvcst_soilt(:,:,:) =                                           &
                            SPREAD(psparms%smvcst_soilt(:,:,1), 3, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_sat_const_z_soilt- failure in logic")
    END IF
  ELSE
    psparms%smvcst_soilt(:,:,:) = const_val
  END IF
  !=== End cases for sm_sat ===

  !=== Start cases for sm_crit ===
CASE ( 'sm_crit' )  !With levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%smvccl_soilt(:,1,1:sm_levels))
      psparms%smvccl_soilt(:,:,:) = SPREAD(psparms%smvccl_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%smvccl_soilt(:,1,1:sm_levels))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_crit- failure in logic")
    END IF
  ELSE
    psparms%smvccl_soilt(:,:,:) = const_val
  END IF

CASE ( 'sm_crit_const_z' ) !No levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      ! Spread around soilt and layer dims
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%smvccl_soilt(:,1,1))
      psparms%smvccl_soilt(:,1,1:sm_levels) =                                 &
                            SPREAD(psparms%smvccl_soilt(:,1,1), 2, sm_levels)
      psparms%smvccl_soilt(:,:,:)           =                                 &
                            SPREAD(psparms%smvccl_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      ! Copy the values from the first vertical level to all others
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%smvccl_soilt(:,1,1))
      psparms%smvccl_soilt(:,1,1:sm_levels) =                                 &
                            SPREAD(psparms%smvccl_soilt(:,1,1), 2, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_crit_const_z- failure in logic")
    END IF
  ELSE
    psparms%smvccl_soilt(:,:,:) = const_val
  END IF

CASE ( 'sm_crit_soilt' ) !With levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%smvccl_soilt(:,:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_crit_soilt- failure in logic")
    END IF
  ELSE
    psparms%smvccl_soilt(:,:,:) = const_val
  END IF

CASE ( 'sm_crit_const_z_soilt' ) !No levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Spread around layer dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%smvccl_soilt(:,1:nsoilt,1))
      psparms%smvccl_soilt(:,:,:) =                                           &
                            SPREAD(psparms%smvccl_soilt(:,:,1), 3, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_crit_const_z_soilt- failure in logic")
    END IF
  ELSE
    psparms%smvccl_soilt(:,:,:) = const_val
  END IF
  !=== End cases for sm_crit ===

  !=== Start cases for sm_wilt ===
CASE ( 'sm_wilt' )  !With levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%smvcwt_soilt(:,1,1:sm_levels))
      psparms%smvcwt_soilt(:,:,:) =                                           &
                            SPREAD(psparms%smvcwt_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%smvcwt_soilt(:,1,1:sm_levels))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_wilt- failure in logic")
    END IF
  ELSE
    psparms%smvcwt_soilt(:,:,:) = const_val
  END IF

CASE ( 'sm_wilt_const_z' ) !No levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      ! Spread around soilt and layer dims
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%smvcwt_soilt(:,1,1))
      psparms%smvcwt_soilt(:,1,1:sm_levels) =                                 &
                            SPREAD(psparms%smvcwt_soilt(:,1,1), 2, sm_levels)
      psparms%smvcwt_soilt(:,:,:)           =                                 &
                            SPREAD(psparms%smvcwt_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      ! Copy the values from the first vertical level to all others
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%smvcwt_soilt(:,1,1))
      psparms%smvcwt_soilt(:,1,1:sm_levels) =                                 &
                            SPREAD(psparms%smvcwt_soilt(:,1,1), 2, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_wilt_const_z- failure in logic")
    END IF
  ELSE
    psparms%smvcwt_soilt(:,:,:) = const_val
  END IF

CASE ( 'sm_wilt_soilt' ) !With levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%smvcwt_soilt(:,:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_wilt_soilt- failure in logic")
    END IF
  ELSE
    psparms%smvcwt_soilt(:,:,:) = const_val
  END IF

CASE ( 'sm_wilt_const_z_soilt' ) !No levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Spread around layer dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%smvcwt_soilt(:,1:nsoilt,1))
      psparms%smvcwt_soilt(:,:,:) =                                           &
                            SPREAD(psparms%smvcwt_soilt(:,:,1), 3, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_wilt_const_z_soilt- failure in logic")
    END IF
  ELSE
    psparms%smvcwt_soilt(:,:,:) = const_val
  END IF
  !=== End cases for sm_wilt ===

  !=== Start cases for hcap ===
CASE ( 'hcap' )  !With levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%hcap_soilt(:,1,1:sm_levels))
      psparms%hcap_soilt(:,:,:) = SPREAD(psparms%hcap_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%hcap_soilt(:,1,1:sm_levels))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate hcap- failure in logic")
    END IF
  ELSE
    psparms%hcap_soilt(:,:,:) = const_val
  END IF

CASE ( 'hcap_const_z' ) !No levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      ! Spread around soilt and layer dims
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%hcap_soilt(:,1,1))
      psparms%hcap_soilt(:,1,1:sm_levels) =                                   &
                            SPREAD(psparms%hcap_soilt(:,1,1), 2, sm_levels)
      psparms%hcap_soilt(:,:,:)           =                                   &
                            SPREAD(psparms%hcap_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      ! Copy the values from the first vertical level to all others
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%hcap_soilt(:,1,1))
      psparms%hcap_soilt(:,1,1:sm_levels) =                                   &
                            SPREAD(psparms%hcap_soilt(:,1,1), 2, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate hcap_const_z- failure in logic")
    END IF
  ELSE
    psparms%hcap_soilt(:,:,:) = const_val
  END IF

CASE ( 'hcap_soilt' ) !With levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%hcap_soilt(:,:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate psparms%hcap_soilt- failure in logic")
    END IF
  ELSE
    psparms%hcap_soilt(:,:,:) = const_val
  END IF

CASE ( 'hcap_const_z_soilt' ) !No levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Spread around layer dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%hcap_soilt(:,1:nsoilt,1))
      psparms%hcap_soilt(:,:,:) =                                             &
                          SPREAD(psparms%hcap_soilt(:,:,1), 3, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate hcap_const_z_soilt- failure in logic")
    END IF
  ELSE
    psparms%hcap_soilt(:,:,:) = const_val
  END IF
  !=== End cases for hcap ===

  !=== Start cases for hcon ===
CASE ( 'hcon' )  !With levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%hcon_soilt(:,1,1:sm_levels))
      psparms%hcon_soilt(:,:,:) = SPREAD(psparms%hcon_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%hcon_soilt(:,1,1:sm_levels))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate hcon- failure in logic")
    END IF
  ELSE
    psparms%hcon_soilt(:,:,:) = const_val
  END IF

CASE ( 'hcon_const_z' ) !No levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      ! Spread around soilt and layer dims
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%hcon_soilt(:,1,1))
      psparms%hcon_soilt(:,1,1:sm_levels) =                                   &
                        SPREAD(psparms%hcon_soilt(:,1,1), 2, sm_levels)
      psparms%hcon_soilt(:,:,:)           =                                   &
                        SPREAD(psparms%hcon_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      ! Copy the values from the first vertical level to all others
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%hcon_soilt(:,1,1))
      psparms%hcon_soilt(:,1,1:sm_levels) =                                   &
                        SPREAD(psparms%hcon_soilt(:,1,1), 2, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate hcon_const_z- failure in logic")
    END IF
  ELSE
    psparms%hcon_soilt(:,:,:) = const_val
  END IF

CASE ( 'hcon_soilt' ) !With levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%hcon_soilt(:,:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate psparms%hcon_soilt- failure in logic")
    END IF
  ELSE
    psparms%hcon_soilt(:,:,:) = const_val
  END IF

CASE ( 'hcon_const_z_soilt' ) !No levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Spread around layer dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%hcon_soilt(:,1:nsoilt,1))
      psparms%hcon_soilt(:,:,:) =                                             &
                        SPREAD(psparms%hcon_soilt(:,:,1), 3, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate hcon_const_z_soilt- failure in logic")
    END IF
  ELSE
    psparms%hcon_soilt(:,:,:) = const_val
  END IF
  !=== End cases for hcon ===

  !-----------------------------------------------------------------------------
  ! Soil ancillaries on soil carbon layers.
  !-----------------------------------------------------------------------------

      !=== Start cases for soil_ph ===
CASE ( 'soil_ph' )  !With levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%soil_ph_soilt(:,1,1:dim_cslayer))
      psparms%soil_ph_soilt(:,:,:) =                                          &
                        SPREAD(psparms%soil_ph_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%soil_ph_soilt(:,1,1:dim_cslayer))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate soil_ph- failure in logic")
    END IF
  ELSE
    psparms%soil_ph_soilt(:,:,:) = const_val
  END IF

CASE ( 'soil_ph_const_z' ) !No levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      ! Spread around soilt and layer dims
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%soil_ph_soilt(:,1,1))
      psparms%soil_ph_soilt(:,1,1:dim_cslayer) =                              &
                        SPREAD(psparms%soil_ph_soilt(:,1,1), 2, dim_cslayer)
      psparms%soil_ph_soilt(:,:,:)             =                              &
                        SPREAD(psparms%soil_ph_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      ! Copy the values from the first vertical level to all others
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%soil_ph_soilt(:,1,1))
      psparms%soil_ph_soilt(:,1,1:dim_cslayer) =                              &
                        SPREAD(psparms%soil_ph_soilt(:,1,1), 2, dim_cslayer)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate soil_ph_const_z- failure in logic")
    END IF
  ELSE
    psparms%soil_ph_soilt(:,:,:) = const_val
  END IF

CASE ( 'soil_ph_soilt' ) !With levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%soil_ph_soilt(:,:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate psparms%soil_ph_soilt- failure in logic")
    END IF
  ELSE
    psparms%soil_ph_soilt(:,:,:) = const_val
  END IF

CASE ( 'soil_ph_const_z_soilt' ) !No levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Spread around layer dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, psparms%soil_ph_soilt(:,1:nsoilt,1))
      psparms%soil_ph_soilt(:,:,:) = SPREAD(psparms%soil_ph_soilt(:,:,1), 3,  &
                                      dim_cslayer)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate soil_ph_const_z_soilt- " //             &
                     "failure in logic")
    END IF
  ELSE
    psparms%soil_ph_soilt(:,:,:) = const_val
  END IF
  !=== End cases for soil_ph ===

  !-----------------------------------------------------------------------------
  ! TOPMODEL variables are just on land points
  !
  ! Each variable can be read in with _soilt if the soil tile dimension is
  ! present.
  ! There are a few possible outcomes:
  !
  ! No suffix:
  !  l_tile_soil = T and l_broadcast_ancils = T
  !  --Spread around soilt dim
  !  l_tile_soil = F
  !  --Do nothing
  !
  ! _soilt suffix:
  ! l_tile_soil = T and l_broadcast_ancils = F
  ! --Do nothing
  !-----------------------------------------------------------------------------
      !=== Start cases for fexp ===
CASE ( 'fexp' )  !no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, toppdm%fexp_soilt(:,1))
      toppdm%fexp_soilt(:,:) = SPREAD(toppdm%fexp_soilt(:,1), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, toppdm%fexp_soilt(:,1))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate fexp- failure in logic")
    END IF
  ELSE
    toppdm%fexp_soilt(:,:) = const_val
  END IF

CASE ( 'toppdm%fexp_soilt' ) !with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, toppdm%fexp_soilt(:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate toppdm%fexp_soilt- failure in logic")
    END IF
  ELSE
    toppdm%fexp_soilt(:,:) = const_val
  END IF
  !=== End cases for fexp ===

  !=== Start cases for ti_mean ===
CASE ( 'ti_mean' )  !no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, toppdm%ti_mean_soilt(:,1))
      toppdm%ti_mean_soilt(:,:) = SPREAD(toppdm%ti_mean_soilt(:,1), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, toppdm%ti_mean_soilt(:,1))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate ti_mean- failure in logic")
    END IF
  ELSE
    toppdm%ti_mean_soilt(:,:) = const_val
  END IF

CASE ( 'toppdm%ti_mean_soilt' ) !with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, toppdm%ti_mean_soilt(:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate toppdm%ti_mean_soilt- failure in logic")
    END IF
  ELSE
    toppdm%ti_mean_soilt(:,:) = const_val
  END IF
  !=== End cases for ti_mean ===

  !=== Start cases for ti_sig ===
CASE ( 'ti_sig' )  !no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, toppdm%ti_sig_soilt(:,1))
      toppdm%ti_sig_soilt(:,:) = SPREAD(toppdm%ti_sig_soilt(:,1), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, toppdm%ti_sig_soilt(:,1))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate ti_sig- failure in logic")
    END IF
  ELSE
    toppdm%ti_sig_soilt(:,:) = const_val
  END IF

CASE ( 'toppdm%ti_sig_soilt' ) !with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, toppdm%ti_sig_soilt(:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate toppdm%ti_sig_soilt- failure in logic")
    END IF
  ELSE
    toppdm%ti_sig_soilt(:,:) = const_val
  END IF
  !=== End cases for ti_sig ===

  !-----------------------------------------------------------------------------
  ! Agricultural fraction is on land points only
  !-----------------------------------------------------------------------------
CASE ( 'frac_agr' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, trifctltype%frac_agr_gb)
  ELSE
    trifctltype%frac_agr_gb(:) = const_val
  END IF

CASE ( 'frac_past' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, trif_vars%frac_past_gb)
  ELSE
    trif_vars%frac_past_gb(:) = const_val
  END IF

  !-----------------------------------------------------------------------------
  ! Irrigation fraction is on land points only
  !-----------------------------------------------------------------------------
CASE ( 'frac_irrig' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, crop_vars%frac_irr_all(:,1))
  ELSE
    crop_vars%frac_irr_all(:,1) = const_val
  END IF

CASE ( 'irrfrac_irrtiles' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, crop_vars%irrfrac_irrtiles(:,1))
  ELSE
    crop_vars%irrfrac_irrtiles(:,1) = const_val
  END IF

CASE ( 'frac_agr_prev' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%frac_agr_prev_gb)
  ELSE
    progs%frac_agr_prev_gb(:) = const_val
  END IF

CASE ( 'frac_past_prev' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%frac_past_prev_gb)
  ELSE
    progs%frac_past_prev_gb(:) = const_val
  END IF

CASE ( 'wood_prod_fast' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%wood_prod_fast_gb)
  ELSE
    progs%wood_prod_fast_gb(:) = const_val
  END IF

CASE ( 'wood_prod_med' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%wood_prod_med_gb)
  ELSE
    progs%wood_prod_med_gb(:) = const_val
  END IF

CASE ( 'wood_prod_slow' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%wood_prod_slow_gb)
  ELSE
    progs%wood_prod_slow_gb(:) = const_val
  END IF

CASE ( 'n_inorg_avail_pft' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%n_inorg_avail_pft)
  ELSE
    progs%n_inorg_avail_pft(:,:,:) = const_val
  END IF

CASE ( 'n_inorg' )
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_soilt) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, progs%n_inorg_soilt_lyrs(:,1,:))
      progs%n_inorg_soilt_lyrs(:,:,:) =                                       &
        SPREAD(progs%n_inorg_soilt_lyrs(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, progs%n_inorg_soilt_lyrs(:,1,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate n_inorg- failure in logic")
    END IF
  ELSE
    progs%n_inorg_soilt_lyrs(:,:,:) = const_val
  END IF

CASE ( 'substr_ch4' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%substr_ch4)
  ELSE
    progs%substr_ch4(:,:) = const_val
  END IF

CASE ( 'mic_ch4' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%mic_ch4)
  ELSE
    progs%mic_ch4(:,:) = const_val
  END IF

CASE ( 'mic_act_ch4' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%mic_act_ch4)
  ELSE
    progs%mic_act_ch4(:,:) = const_val
  END IF

CASE ( 'acclim_ch4' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%acclim_ch4)
  ELSE
    progs%acclim_ch4(:,:) = const_val
  END IF

  !-----------------------------------------------------------------------------
  ! Urban variables are on land points only
  !-----------------------------------------------------------------------------
CASE ( 'wrr' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, urban_param%wrr_gb)
  ELSE
    urban_param%wrr_gb(:) = const_val
  END IF

CASE ( 'hwr' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, urban_param%hwr_gb)
  ELSE
    urban_param%hwr_gb(:) = const_val
  END IF

CASE ( 'hgt' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, urban_param%hgt_gb)
  ELSE
    urban_param%hgt_gb(:) = const_val
  END IF

CASE ( 'ztm' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, urban_param%ztm_gb)
  ELSE
    urban_param%ztm_gb(:) = const_val
  END IF

CASE ( 'disp' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, urban_param%disp_gb)
  ELSE
    urban_param%disp_gb(:) = const_val
  END IF

CASE ( 'albwl' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, urban_param%albwl_gb)
  ELSE
    urban_param%albwl_gb(:) = const_val
  END IF

CASE ( 'albrd' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, urban_param%albrd_gb)
  ELSE
    urban_param%albrd_gb(:) = const_val
  END IF

CASE ( 'emisw' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, urban_param%emisw_gb)
  ELSE
    urban_param%emisw_gb(:) = const_val
  END IF

CASE ( 'emisr' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, urban_param%emisr_gb)
  ELSE
    urban_param%emisr_gb(:) = const_val
  END IF

  !-----------------------------------------------------------------------------
  ! FLake variables are on land points only
  !-----------------------------------------------------------------------------
CASE ( 'lake_depth' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, lake_depth_gb(:))
  ELSE
    lake_depth_gb(:) = const_val
  END IF

CASE ( 'lake_fetch_gb' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, lake_fetch_gb(:))
  ELSE
    lake_fetch_gb(:) = const_val
  END IF

CASE ( 'lake_h_mxl_gb' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, lake_h_mxl_gb(:))
  ELSE
    lake_h_mxl_gb(:) = const_val
  END IF

CASE ( 'lake_shape_factor_gb' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, lake_shape_factor_gb(:))
  ELSE
    lake_shape_factor_gb(:) = const_val
  END IF

CASE ( 'lake_t_mean_gb' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, lake_t_mean_gb(:))
  ELSE
    lake_t_mean_gb(:) = const_val
  END IF

CASE ( 'lake_t_mxl_gb' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, lake_t_mxl_gb(:))
  ELSE
    lake_t_mxl_gb(:) = const_val
  END IF

CASE ( 'lake_t_ice_gb' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, lake_t_ice_gb(:))
  ELSE
    lake_t_ice_gb(:) = const_val
  END IF

CASE ( 'lake_h_ice_gb' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, lake_h_ice_gb(:))
  ELSE
    lake_h_ice_gb(:) = const_val
  END IF

CASE ( 'lake_t_sfc_gb' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, lake_t_sfc_gb(:))
  ELSE
    lake_t_sfc_gb(:) = const_val
  END IF

  !-----------------------------------------------------------------------------
  ! Variables that are set as initial conditions
  ! These are all land points only, but with varying numbers of levels

  ! Some variables can be read in with _soilt if the soil tile dimension is
  ! present. There are a few possible outcomes
  !
  ! No suffix:
  !  l_tile_soil = T and l_broadcast_soilt = T
  !  --Spread around soilt dim
  !  l_tile_soil = F
  !  --Do nothing
  !
  ! _soilt suffix:
  ! l_tile_soil = T and l_broadcast_soilt = F
  ! --Do nothing
  !-----------------------------------------------------------------------------

CASE ( 'canopy' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%canopy_surft)
  ELSE
    progs%canopy_surft(:,:) = const_val
  END IF

CASE ( 'cs' )
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_soilt) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, progs%cs_pool_soilt(:,1,:,:))
      progs%cs_pool_soilt(:,:,:,:) =                                          &
                        SPREAD(progs%cs_pool_soilt(:,1,:,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, progs%cs_pool_soilt(:,1,:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate cs- failure in logic")
    END IF
  ELSE
    progs%cs_pool_soilt(:,:,:,:) = const_val
  END IF

CASE ( 'cs_soilt' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%cs_pool_soilt(:,:,:,:))
  ELSE
    progs%cs_pool_soilt(:,:,:,:) = const_val
  END IF

CASE ( 'ns' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%ns_pool_gb)
  ELSE
    progs%ns_pool_gb(:,:,:) = const_val
  END IF

CASE ( 'deposition_n' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, trif_vars%deposition_n_gb)
  ELSE
    trif_vars%deposition_n_gb(:) = const_val
  END IF

CASE ( 'trif_vars%g_burn_pft' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, trif_vars%g_burn_pft)
  ELSE
    trif_vars%g_burn_pft(:,:) = const_val
  END IF

CASE ( 'gs' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%gs_gb)
  ELSE
    progs%gs_gb(:) = const_val
  END IF

CASE ( 'snow_tile' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%snow_surft)
  ELSE
    progs%snow_surft(:,:) = const_val
  END IF

CASE ( 'sthuf' )
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_soilt) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, jules_vars%sthuf_soilt(:,1,:))
      jules_vars%sthuf_soilt(:,:,:) = SPREAD(jules_vars%sthuf_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, jules_vars%sthuf_soilt(:,1,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sthuf- failure in logic")
    END IF
  ELSE
    jules_vars%sthuf_soilt(:,:,:) = const_val
  END IF

CASE ( 'sthuf_soilt' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, jules_vars%sthuf_soilt)
  ELSE
    jules_vars%sthuf_soilt(:,:,:) = const_val
  END IF

CASE ( 't_soil' )
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_soilt) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, progs%t_soil_soilt(:,1,:))
      progs%t_soil_soilt(:,:,:) = SPREAD(progs%t_soil_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, progs%t_soil_soilt(:,1,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate t_soil- failure in logic")
    END IF
  ELSE
    progs%t_soil_soilt(:,:,:) = const_val
  END IF

CASE ( 't_soil_soilt' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%t_soil_soilt)
  ELSE
    progs%t_soil_soilt(:,:,:) = const_val
  END IF

CASE ( 'tsoil_deep' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%tsoil_deep_gb)
  ELSE
    progs%tsoil_deep_gb(:,:) = const_val
  END IF

CASE ( 'sthu_irr' )
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_soilt) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, crop_vars%sthu_irr_soilt(:,1,:))
      crop_vars%sthu_irr_soilt(:,:,:) =                                       &
        SPREAD(crop_vars%sthu_irr_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, crop_vars%sthu_irr_soilt(:,1,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sthu_irr- failure in logic")
    END IF
  ELSE
    crop_vars%sthu_irr_soilt(:,:,:) = const_val
  END IF

CASE ( 'tstar_tile' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%tstar_surft)
  ELSE
    progs%tstar_surft(:,:) = const_val
  END IF

CASE ( 'tsurf_elev_surft' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%tsurf_elev_surft)
  ELSE
    progs%tsurf_elev_surft(:,:) = const_val
  END IF

CASE ( 'cropdvi' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, crop_vars%dvi_cpft)
  ELSE
    crop_vars%dvi_cpft(:,:) = const_val
  END IF

CASE ( 'croprootc' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, crop_vars%rootc_cpft)
  ELSE
    crop_vars%rootc_cpft(:,:) = const_val
  END IF

CASE ( 'cropharvc' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, crop_vars%harvc_cpft)
  ELSE
    crop_vars%harvc_cpft(:,:) = const_val
  END IF

CASE ( 'cropreservec' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, crop_vars%reservec_cpft)
  ELSE
    crop_vars%reservec_cpft(:,:) = const_val
  END IF

CASE ( 'croplai' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, crop_vars%croplai_cpft)
  ELSE
    crop_vars%croplai_cpft(:,:) = const_val
  END IF

CASE ( 'cropcanht' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, crop_vars%cropcanht_cpft)
  ELSE
    crop_vars%cropcanht_cpft(:,:) = const_val
  END IF

CASE ( 'cropsowdate' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, crop_vars%sow_date_cpft)
  ELSE
    crop_vars%sow_date_cpft(:,:) = const_val
  END IF

CASE ( 'croplatestharvdate' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, crop_vars%latestharv_date_cpft)
  ELSE
    crop_vars%latestharv_date_cpft(:,:) = const_val
  END IF

CASE ( 'cropttveg' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, crop_vars%tt_veg_cpft)
  ELSE
    crop_vars%tt_veg_cpft(:,:) = const_val
  END IF

CASE ( 'cropttrep' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, crop_vars%tt_rep_cpft)
  ELSE
    crop_vars%tt_rep_cpft(:,:) = const_val
  END IF

CASE ( 'lai' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%lai_pft)
  ELSE
    progs%lai_pft(:,:) = const_val
  END IF

CASE ( 'canht' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%canht_pft)
  ELSE
    progs%canht_pft(:,:) = const_val
  END IF

CASE ( 'sthzw' )
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_soilt) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, toppdm%sthzw_soilt(:,1))
      toppdm%sthzw_soilt(:,:) = SPREAD(toppdm%sthzw_soilt(:,1), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, toppdm%sthzw_soilt(:,1))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sthzw- failure in logic")
    END IF
  ELSE
    toppdm%sthzw_soilt(:,:) = const_val
  END IF

CASE ( 'toppdm%sthzw_soilt' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, toppdm%sthzw_soilt)
  ELSE
    toppdm%sthzw_soilt(:,:) = const_val
  END IF

CASE ( 'zw' )
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_soilt) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, toppdm%zw_soilt(:,1))
      toppdm%zw_soilt(:,:) = SPREAD(toppdm%zw_soilt(:,1), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, toppdm%zw_soilt(:,1))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate zw- failure in logic")
    END IF
  ELSE
    toppdm%zw_soilt(:,:) = const_val
  END IF

CASE ( 'toppdm%zw_soilt' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, toppdm%zw_soilt)
  ELSE
    toppdm%zw_soilt(:,:) = const_val
  END IF

CASE ( 'slope' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, toppdm%slope_gb)
  ELSE
    toppdm%slope_gb(:) = const_val
  END IF

CASE ( 'rgrain' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%rgrain_surft)
  ELSE
    progs%rgrain_surft(:,:) = const_val
  END IF

CASE ( 'cv' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, trifctltype%cv_gb)
  ELSE
    trifctltype%cv_gb(:) = const_val
  END IF

CASE ( 'rho_snow' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%rho_snow_grnd_surft)
  ELSE
    progs%rho_snow_grnd_surft(:,:) = const_val
  END IF

CASE ( 'snow_depth' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%snowdepth_surft)
  ELSE
    progs%snowdepth_surft(:,:) = const_val
  END IF

CASE ( 'snow_grnd' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%snow_grnd_surft)
  ELSE
    progs%snow_grnd_surft(:,:) = const_val
  END IF

CASE ( 'nsnow' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, nsnow_real)
    progs%nsnow_surft(:,:) = NINT(nsnow_real)
  ELSE
    progs%nsnow_surft(:,:) = NINT(const_val)
  END IF

CASE ( 'snow_ds' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%ds_surft)
  ELSE
    progs%ds_surft(:,:,:) = const_val
  END IF

CASE ( 'snow_ice' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%sice_surft)
  ELSE
    progs%sice_surft(:,:,:) = const_val
  END IF

CASE ( 'snow_liq' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%sliq_surft)
  ELSE
    progs%sliq_surft(:,:,:) = const_val
  END IF

CASE ( 'tsnow' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%tsnow_surft)
  ELSE
    progs%tsnow_surft(:,:,:) = const_val
  END IF

CASE ( 'rgrainl' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, progs%rgrainl_surft)
  ELSE
    progs%rgrainl_surft(:,:,:) = const_val
  END IF

  !-----------------------------------------------------------------------------
  ! Forcing variables
  !-----------------------------------------------------------------------------
CASE ( 'pstar' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, pstar_ij)
  ELSE
    pstar_ij(:,:) = const_val
  END IF

CASE ( 'q' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, qw_1_ij)
  ELSE
    qw_1_ij(:,:) = const_val
  END IF
  qw_1_ij(:,:) = MAX(qw_1_ij, 0.0)

CASE ( 't' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, tl_1_ij)
  ELSE
    tl_1_ij(:,:) = const_val
  END IF

CASE ( 'rad_net' )
  ! Net downward radiation is stored in lw_down until it is processed
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, lw_down_ij)
  ELSE
    lw_down_ij(:,:) = const_val
  END IF

CASE ( 'lw_net' )
  ! Net LW downward radiation is stored in lw_down until it is processed
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, lw_down_ij)
  ELSE
    lw_down_ij(:,:) = const_val
  END IF

CASE ( 'sw_net' )
  ! Net SW downward radiation is stored in sw_down until it is processed
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, sw_down_ij)
  ELSE
    sw_down_ij(:,:) = const_val
  END IF

CASE ( 'lw_down' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, lw_down_ij)
  ELSE
    lw_down_ij(:,:) = const_val
  END IF

CASE ( 'sw_down' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, sw_down_ij)
  ELSE
    sw_down_ij(:,:) = const_val
  END IF

CASE ( 'diff_rad' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, diff_rad_ij)
  ELSE
    diff_rad_ij(:,:) = const_val
  END IF

CASE ( 'dt_range' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, diurnal_temperature_range_ij)
  ELSE
    diurnal_temperature_range_ij(:,:) = const_val
  END IF

CASE ( 'precip' )
  ! Store total precip as large-scale rainfall until it is partitioned
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, ls_rain_ij)
  ELSE
    ls_rain_ij(:,:) = const_val
  END IF

CASE ( 'tot_rain' )
  ! Store total rainfall as large-scale until it is partitioned
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, ls_rain_ij)
  ELSE
    ls_rain_ij(:,:) = const_val
  END IF

CASE ( 'tot_snow' )
  ! If given total snow, we assume it is all large-scale
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, ls_snow_ij)
  ELSE
    ls_snow_ij(:,:) = const_val
  END IF

CASE ( 'con_rain' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, con_rain_ij)
  ELSE
    con_rain_ij(:,:) = const_val
  END IF

CASE ( 'ls_rain' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, ls_rain_ij)
  ELSE
    ls_rain_ij(:,:) = const_val
  END IF

CASE ( 'con_snow' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, con_snow_ij)
  ELSE
    con_snow_ij(:,:) = const_val
  END IF

CASE ( 'ls_snow' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, ls_snow_ij)
  ELSE
    ls_snow_ij(:,:) = const_val
  END IF

CASE ( 'wind' )
  ! Wind speed just goes directly into u component
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, u_1_ij)
  ELSE
    u_1_ij(:,:) = const_val
  END IF

CASE ( 'u' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, u_1_ij)
  ELSE
    u_1_ij(:,:) = const_val
  END IF

CASE ( 'v' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, v_1_ij)
  ELSE
    v_1_ij(:,:) = const_val
  END IF

CASE ( 'z1_tq_in' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, ainfo%z1_tq_ij)
  ELSE
    ainfo%z1_tq_ij(:,:) = const_val
  END IF

  !-----------------------------------------------------------------------------
  ! River routing variables
  !-----------------------------------------------------------------------------

CASE ( 'latitude_2d' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, rivers_lat2d)
  ELSE
    rivers_lat2d(:,:) = const_val
  END IF

CASE ( 'longitude_2d' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, rivers_lon2d)
  ELSE
    rivers_lon2d(:,:) = const_val
  END IF

CASE ( 'direction' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, rivers_dir)
  ELSE
    rivers_dir(:,:) = const_val
  END IF

CASE ( 'sequence' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, rivers_seq)
  ELSE
    rivers_seq(:,:) = const_val
  END IF

CASE ( 'area' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, rivers_dra)
  ELSE
    rivers_dra(:,:) = const_val
  END IF

CASE ( 'rivers_ygrid' )
  IF ( PRESENT(cube) ) THEN
    cube_land = cube_create((/ SIZE(cube%values) /))
    cube_land%values(:) = RESHAPE(cube%values, (/ SIZE(cube%values) /))
    CALL cube_get_data(cube_land, rivers_ygrid)
  ELSE
    rivers_ygrid(:) = const_val
  END IF

CASE ( 'rivers_xgrid' )
  IF ( PRESENT(cube) ) THEN
    cube_land = cube_create((/ SIZE(cube%values) /))
    cube_land%values(:) = RESHAPE(cube%values, (/ SIZE(cube%values) /))
    CALL cube_get_data(cube_land, rivers_xgrid)
  ELSE
    rivers_xgrid(:) = const_val
  END IF

CASE ( 'rivers_sto_rp' )
  IF ( PRESENT(cube) ) THEN
    cube_land = cube_create((/ SIZE(cube%values) /))
    CALL cube_get_data(cube_land, rivers_sto_rp)
  ELSE
    rivers_sto_rp(:) = const_val
  END IF

CASE ( 'rfm_surfstore_rp' )
  IF ( PRESENT(cube) ) THEN
    cube_land = cube_create((/ SIZE(cube%values) /))
    CALL cube_get_data(cube_land, rfm_surfstore_rp)
  ELSE
    rfm_surfstore_rp(:) = const_val
  END IF

CASE ( 'rfm_substore_rp' )
  IF ( PRESENT(cube) ) THEN
    cube_land = cube_create((/ SIZE(cube%values) /))
    CALL cube_get_data(cube_land, rfm_substore_rp)
  ELSE
    rfm_substore_rp(:) = const_val
  END IF

CASE ( 'rfm_flowin_rp' )
  IF ( PRESENT(cube) ) THEN
    cube_land = cube_create((/ SIZE(cube%values) /))
    CALL cube_get_data(cube_land, rfm_flowin_rp)
  ELSE
    rfm_flowin_rp(:) = const_val
  END IF

CASE ( 'rfm_bflowin_rp' )
  IF ( PRESENT(cube) ) THEN
    cube_land = cube_create((/ SIZE(cube%values) /))
    CALL cube_get_data(cube_land, rfm_bflowin_rp)
  ELSE
    rfm_bflowin_rp(:) = const_val
  END IF

  !-----------------------------------------------------------------------------
  ! Overbank inundation variables
  !-----------------------------------------------------------------------------

CASE ( 'logn_mean' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, logn_mean)
  ELSE
    logn_mean(:,:) = const_val
  END IF

CASE ( 'logn_stdev' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, logn_stdev)
  ELSE
    logn_stdev(:,:) = const_val
  END IF

  !-----------------------------------------------------------------------------
  ! Fire and metstats variables- land points only
  !-----------------------------------------------------------------------------
  ! metstats module cases- land points only
CASE ( 'temp_max_00h_r' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%temp_max_00h%run)
  ELSE
    metstats_prog(:)%temp_max_00h%run = const_val
  END IF

CASE ( 'temp_ave_00h_r' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%temp_ave_00h%run)
  ELSE
    metstats_prog(:)%temp_ave_00h%run = const_val
  END IF

CASE ( 'prec_tot_00h_r' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%prec_tot_00h%run)
  ELSE
    metstats_prog(:)%prec_tot_00h%run = const_val
  END IF

CASE ( 'prec_tot_12h_r' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%prec_tot_12h%run)
  ELSE
    metstats_prog(:)%prec_tot_12h%run = const_val
  END IF

CASE ( 'rhum_min_00h_r' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%rhum_min_00h%run)
  ELSE
    metstats_prog(:)%rhum_min_00h%run = const_val
  END IF

CASE ( 'dewp_ave_00h_r' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%dewp_ave_00h%run)
  ELSE
    metstats_prog(:)%dewp_ave_00h%run = const_val
  END IF

CASE ( 'wind_ave_00h_r' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%wind_ave_00h%run)
  ELSE
    metstats_prog(:)%wind_ave_00h%run = const_val
  END IF

CASE ( 'temp_max_00h' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%temp_max_00h%fin)
  ELSE
    metstats_prog(:)%temp_max_00h%fin = const_val
  END IF

CASE ( 'temp_ave_00h' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%temp_ave_00h%fin)
  ELSE
    metstats_prog(:)%temp_ave_00h%fin = const_val
  END IF

CASE ( 'temp_ave_nday' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%temp_ave_nday%fin)
  ELSE
    metstats_prog(:)%temp_ave_nday%fin = const_val
  END IF

CASE ( 'temp_pnt_12h' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%temp_pnt_12h%fin)
  ELSE
    metstats_prog(:)%temp_pnt_12h%fin = const_val
  END IF

CASE ( 'prec_tot_00h' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%prec_tot_00h%fin)
  ELSE
    metstats_prog(:)%prec_tot_00h%fin = const_val
  END IF

CASE ( 'prec_tot_12h' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%prec_tot_12h%fin)
  ELSE
    metstats_prog(:)%prec_tot_12h%fin = const_val
  END IF

CASE ( 'rhum_min_00h' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%rhum_min_00h%fin)
  ELSE
    metstats_prog(:)%rhum_min_00h%fin = const_val
  END IF

CASE ( 'rhum_pnt_12h' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%rhum_pnt_12h%fin)
  ELSE
    metstats_prog(:)%rhum_pnt_12h%fin = const_val
  END IF

CASE ( 'dewp_ave_00h' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%dewp_ave_00h%fin)
  ELSE
    metstats_prog(:)%dewp_ave_00h%fin = const_val
  END IF

CASE ( 'wind_ave_00h' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%wind_ave_00h%fin)
  ELSE
    metstats_prog(:)%wind_ave_00h%fin = const_val
  END IF

CASE ( 'wind_pnt_12h' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%wind_pnt_12h%fin)
  ELSE
    metstats_prog(:)%wind_pnt_12h%fin = const_val
  END IF

  ! Fire module variables- land points only
CASE ( 'fire_mcarthur_r_dr' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, fire_prog(:)%mcarthur%r_dr)
  ELSE
    fire_prog(:)%mcarthur%r_dr = const_val
  END IF

CASE ( 'fire_mcarthur_n_dr' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, fire_prog(:)%mcarthur%n_dr)
  ELSE
    fire_prog(:)%mcarthur%n_dr = const_val
  END IF

CASE ( 'fire_canadian_ffmc' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, fire_prog(:)%canadian%ffmc)
  ELSE
    fire_prog(:)%canadian%ffmc = const_val
  END IF

CASE ( 'fire_canadian_ffmc_mois' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, fire_prog(:)%canadian%ffmc_mois)
  ELSE
    fire_prog(:)%canadian%ffmc_mois = const_val
  END IF

CASE ( 'fire_canadian_dmc' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, fire_prog(:)%canadian%dmc)
  ELSE
    fire_prog(:)%canadian%dmc = const_val
  END IF

CASE ( 'fire_canadian_dc' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, fire_prog(:)%canadian%dc)
  ELSE
    fire_prog(:)%canadian%dc = const_val
  END IF

CASE ( 'fire_nesterov' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, fire_prog(:)%nesterov%findex)
  ELSE
    fire_prog(:)%nesterov%findex = const_val
  END IF

  !-----------------------------------------------------------------------------
  ! Other variables that might be prescribed
  !-----------------------------------------------------------------------------
CASE ( 'bl_height' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, jules_vars%zh)
  ELSE
    jules_vars%zh(:,:) = const_val
  END IF

CASE ( 'level_separation' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, jules_vars%dzl)
  ELSE
    jules_vars%dzl(:,:,:) = const_val
  END IF

CASE ( 'ozone' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, o3_gb)
  ELSE
    o3_gb(:) = const_val
  END IF

CASE ( 'tracer_field' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, tracer_field)
  ELSE
    tracer_field(:,:) = const_val
  END IF

CASE ( 'albobs_sw' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, psparms%albobs_sw_gb)
  ELSE
    psparms%albobs_sw_gb(:) = const_val
  END IF

CASE ( 'albobs_vis' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, psparms%albobs_vis_gb)
  ELSE
    psparms%albobs_vis_gb(:) = const_val
  END IF

CASE ( 'albobs_nir' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, psparms%albobs_nir_gb)
  ELSE
    psparms%albobs_nir_gb(:) = const_val
  END IF

CASE ( 'co2_mmr' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, co2_tmp)

    IF ((MAXVAL(co2_tmp) /= co2_tmp(1)) .OR.                                  &
        (MINVAL(co2_tmp) /= co2_tmp(1))) THEN
      CALL log_fatal("populate_var",                                          &
                     "All land points must have same CO2 concetration")
    END IF
    co2_mmr = co2_tmp(1)
  ELSE
    co2_mmr = const_val
  END IF

CASE ( 'pop_den' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, fire_vars%pop_den)
  ELSE
    fire_vars%pop_den(:) = const_val
  END IF

CASE ( 'flash_rate' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, fire_vars%flash_rate)
  ELSE
    fire_vars%flash_rate(:) = const_val
  END IF

CASE ( 't_growth_gb' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, t_growth_gb)
  ELSE
    t_growth_gb(:) = const_val
  END IF

  !-----------------------------------------------------------------------------
  !   ECOSSE variables.
  !-----------------------------------------------------------------------------
  !   ECOSSE prognostic variables.

CASE ( 'n_soil' )
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_soilt) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, soilecosse%n_soil_pool_soilt(:,1,:,:))
      soilecosse%n_soil_pool_soilt(:,:,:,:) =                                 &
                    SPREAD(soilecosse%n_soil_pool_soilt(:,1,:,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, soilecosse%n_soil_pool_soilt(:,1,:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate n_soil- failure in logic")
    END IF
  ELSE
    soilecosse%n_soil_pool_soilt(:,:,:,:) = const_val
  END IF

  !-----------------------------------------------------------------------------
  ! Variables for water resource management.
  !-----------------------------------------------------------------------------
CASE ( 'conveyance_loss' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, conveyance_loss)
  ELSE
    conveyance_loss(:) = const_val
  END IF

CASE ( 'demand_rate_domestic' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, demand_rate_domestic)
  ELSE
    demand_rate_domestic(:) = const_val
  END IF

CASE ( 'demand_rate_industry' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, demand_rate_industry)
  ELSE
    demand_rate_industry(:) = const_val
  END IF

CASE ( 'demand_rate_livestock' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, demand_rate_livestock)
  ELSE
    demand_rate_livestock(:) = const_val
  END IF

CASE ( 'demand_rate_transfers' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, demand_rate_transfers)
  ELSE
    demand_rate_transfers(:) = const_val
  END IF
 
CASE ( 'irrig_eff' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, irrig_eff)
  ELSE
    irrig_eff(:) = const_val
  END IF
 
CASE ( 'grid_area' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, grid_area_ij)
  ELSE
    grid_area_ij(:,:) = const_val
  END IF
 
CASE ( 'sfc_water_frac' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, sfc_water_frac)
  ELSE
    sfc_water_frac(:) = const_val
  END IF

CASE ( 'SoilTemp_CABLE' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, SoilTemp_CABLE )
  ELSE
    SoilTemp_CABLE(:,:,:) = const_val
  END IF

CASE ( 'FrozenSoilFrac_CABLE' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, FrozenSoilFrac_CABLE )
  ELSE
    FrozenSoilFrac_CABLE(:,:,:) = const_val
  END IF

CASE ( 'SoilMoisture_CABLE' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, SoilMoisture_CABLE )
  ELSE
    SoilMoisture_CABLE(:,:,:) = const_val
  END IF

CASE ( 'SnowDepth_CABLE' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, SnowDepth_CABLE )
  ELSE
    SnowDepth_CABLE(:,:,:) = const_val
  END IF

CASE ( 'SnowAge_CABLE' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, SnowAge_CABLE )
  ELSE
    SnowAge_CABLE(:,:) = const_val
  END IF

CASE ( 'SnowMass_CABLE' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, SnowMass_CABLE )
  ELSE
    SnowMass_CABLE(:,:,:) = const_val
  END IF

CASE ( 'ThreeLayerSnowFlag_CABLE' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, ThreeLayerSnowFlag_CABLE)
  ELSE
    ThreeLayerSnowFlag_CABLE(:,:) = const_val
  END IF

CASE ( 'SnowTemp_CABLE' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, SnowTemp_CABLE )
  ELSE
    SnowTemp_CABLE(:,:,:) = const_val
  END IF

CASE ( 'SnowDensity_CABLE' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, SnowDensity_CABLE)
  ELSE
    SnowDensity_CABLE(:,:,:) = const_val
  END IF

CASE ( 'OneLyrSnowDensity_CABLE' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, OneLyrSnowDensity_CABLE)
  ELSE
    OneLyrSnowDensity_CABLE(:,:) = const_val
  END IF

CASE DEFAULT
  CALL log_fatal("populate_var",                                              &
                 "Unrecognised variable for input - '" //                     &
                 TRIM(get_string_identifier(var_id)) // "'. " //              &
                 "See docs for available variables")
END SELECT


! Free the land cube
CALL cube_free(cube_land)

RETURN

END SUBROUTINE populate_var

! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION map_to_land(cube_grid) RESULT(cube_land)

USE data_cube_mod, ONLY: data_cube, cube_create

USE ancil_info, ONLY: land_pts

USE jules_fields_mod, ONLY: ainfo

USE theta_field_sizes, ONLY: t_i_length, t_j_length

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes a cube of data on the full model grid and maps it to land points only
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(data_cube), INTENT(IN) :: cube_grid  ! The data on the full model grid
                                          ! This cube should have two grid
                                          ! dimensions with all other dimensions
                                          ! being levels dimensions
                                          ! Vertical levels are preserved in
                                          ! the mapping


! Return type
TYPE(data_cube) :: cube_land  ! The data mapped onto land points
                              ! This cube will have a 'land points' dimension
                              ! plus all the levels dimensions from the gridded
                              ! data


! Work variables
REAL, ALLOCATABLE :: data_land(:,:), data_grid(:,:,:)
                              ! Variables to hold the data that use a
                              ! 'collapsed levels dimension', i.e. all the
                              ! vertical levels dimensions are represented
                              ! by one dimension

INTEGER :: nlevs  ! The size of the combined levels dimension

! Work variables
INTEGER :: i, j, l  ! Indexing variables


!-----------------------------------------------------------------------------


! Check the data is on the full model grid
IF ( cube_grid%SHAPE(1) /= t_i_length .OR. cube_grid%SHAPE(2) /= t_j_length ) &
  CALL log_fatal("map_to_land",                                               &
                 "Input data must be on the full model grid")

! Allocate the data_land and data_grid arrays
nlevs = PRODUCT(cube_grid%SHAPE(3:))
ALLOCATE(data_land(land_pts, nlevs))
ALLOCATE(data_grid(t_i_length, t_j_length, nlevs))

! Copy the data from the cube into the data array
! Note that we can't use cube_get_data as our array has the combined z dimension
data_grid(:,:,:) = RESHAPE(cube_grid%values, SHAPE(data_grid))

! Do the mapping
DO l = 1,land_pts
  j = (ainfo%land_index(l) - 1) / t_i_length + 1
  i = ainfo%land_index(l) - (j-1) * t_i_length

  data_land(l,:) = data_grid(i,j,:)
END DO

! Copy the land data into the return cube
cube_land = cube_create((/ land_pts, cube_grid%SHAPE(3:) /))
cube_land%values(:) = RESHAPE(data_land, (/ SIZE(data_land) /))

! Deallocate the arrays
DEALLOCATE(data_land)
DEALLOCATE(data_grid)

RETURN

END FUNCTION map_to_land
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION map_from_land(cube_land) RESULT(cube_grid)

USE io_constants, ONLY: mdi

USE data_cube_mod, ONLY: data_cube, cube_create

USE ancil_info, ONLY: land_pts

USE jules_fields_mod, ONLY: ainfo

USE theta_field_sizes, ONLY: t_i_length, t_j_length

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes a data cube on land points only and maps it onto the full model grid
!   Only sets values at the land points - other values are left untouched
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(data_cube), INTENT(IN) :: cube_land  ! The land point data to map to the
                                          ! full model grid
                                          ! It is assumed that the cube will
                                          ! have a 'land points' dimension
                                          ! with all other dimensions being
                                          ! levels dimensions
                                          ! Vertical levels are preserved in
                                          ! the mapping


! Return type
TYPE(data_cube) :: cube_grid  ! The land data mapped onto the full model
                              ! grid
                              ! This cube will have two grid dimensions plus
                              ! the levels dimensions from data_land


! Work variables
REAL, ALLOCATABLE :: data_land(:,:), data_grid(:,:,:)
                              ! Variables to hold the data that use a
                              ! 'collapsed levels dimension', i.e. all the
                              ! vertical levels dimensions are represented
                              ! by one dimension

INTEGER :: nlevs  ! The size of the combined levels dimension

INTEGER :: i, j, l  ! Indexing variables


!-----------------------------------------------------------------------------


IF ( cube_land%SHAPE(1) /= land_pts )                                         &
  CALL log_fatal("map_from_land", "Input data must be on land points")

! Allocate the data_land and data_grid arrays
nlevs = PRODUCT(cube_land%SHAPE(2:))
ALLOCATE(data_land(land_pts, nlevs))
ALLOCATE(data_grid(t_i_length, t_j_length, nlevs))

! Copy the data from the cube into the data array
! Note that we can't use cube_get_data as our array has the combined z dimension
data_land(:,:) = RESHAPE(cube_land%values, SHAPE(data_land))

! Initialise data_grid to missing data - this is the value it will take at
! non-land points
data_grid(:,:,:) = mdi

! Do the mapping
DO l = 1,land_pts
  j = (ainfo%land_index(l) - 1) / t_i_length + 1
  i = ainfo%land_index(l) - (j-1) * t_i_length

  data_grid(i,j,:) = data_land(l,:)
END DO

! Copy the gridded data into the return cube
cube_grid = cube_create((/ t_i_length, t_j_length, cube_land%SHAPE(2:) /))
cube_grid%values(:) = RESHAPE(data_grid, (/ SIZE(data_grid) /))

! Deallocate the arrays
DEALLOCATE(data_land)
DEALLOCATE(data_grid)

RETURN

END FUNCTION map_from_land

END MODULE model_interface_mod
