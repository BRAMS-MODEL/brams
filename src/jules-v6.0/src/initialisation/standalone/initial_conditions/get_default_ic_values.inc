#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

FUNCTION get_default_ic_values(total_snow) RESULT(defaults_dict)

! USE statments
USE dictionary_mod, ONLY: dict, dict_set, dict_create
USE logging_mod,    ONLY: log_info

USE dump_mod,       ONLY: max_var_dump, required_vars_for_configuration

USE model_interface_mod, ONLY: identifier_len

USE fire_mod,        ONLY: fire_inis

USE lake_mod,        ONLY: lake_fetch_0, lake_h_mxl_0, lake_shape_0,          &
                           lake_T_mean_0, lake_T_mxl_0, lake_T_ice_0,         &
                           lake_H_ice_0

USE metstats_mod,    ONLY: metstats_inis

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Gets default values for initialisation
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Arguments
LOGICAL, INTENT(IN) :: total_snow

CHARACTER(LEN=identifier_len) :: identifiers(max_var_dump)
                             ! The variable identifiers of the required
                             ! variables

CHARACTER(LEN=identifier_len) :: vars_from_ancil(max_var_dump)
                             ! The variable identifiers of the ancil
                             ! variables (not used in this subroutine)

! Work variables
INTEGER :: i !Counter
INTEGER :: nvars_required
INTEGER :: nvars_from_ancil

! Return type
TYPE(dict) :: defaults_dict

LOGICAL, PARAMETER :: l_output_mode = .FALSE.

!End of header
!-----------------------------------------------------------------------------

  !Get the list of variables we need to have default variables for
CALL required_vars_for_configuration(nvars_required, identifiers,             &
                                     nvars_from_ancil, vars_from_ancil,       &
                                     l_output_mode,                           &
                                     total_snow, .FALSE.)

!Create the dictionary containing the default values
defaults_dict = dict_create(nvars_required, 1.0)

!Set the default values one by one
DO i = 1, nvars_required
  SELECT CASE ( identifiers(i) )

    !      CASE ( 'canopy' )
    !      CASE ( 'cs' )
    !      CASE ( 'gs' )
    !      CASE ( 'snow_tile' )
    !      CASE ( 'sthuf' )
    !      CASE ( 't_soil' )
    !      CASE ( 'tstar_tile' )
    !      CASE ( 'lai' )
    !      CASE ( 'canht' )
    !      CASE ( 'frac' )
    !      CASE ( 'sthzw' )
    !      CASE ( 'zw' )
    !      CASE ( 'rgrain' )

          !Metstat module variables
  CASE ( 'temp_max_00h_r' )
    CALL dict_set(defaults_dict, 'temp_max_00h_r', metstats_inis%temp_max_00h)
  CASE ( 'temp_ave_00h_r' )
    CALL dict_set(defaults_dict, 'temp_ave_00h_r', metstats_inis%temp_ave_00h)
  CASE ( 'prec_tot_00h_r' )
    CALL dict_set(defaults_dict, 'prec_tot_00h_r', metstats_inis%prec_tot_00h)
  CASE ( 'prec_tot_12h_r' )
    CALL dict_set(defaults_dict, 'prec_tot_12h_r', metstats_inis%prec_tot_12h)
  CASE ( 'rhum_min_00h_r' )
    CALL dict_set(defaults_dict, 'rhum_min_00h_r', metstats_inis%rhum_min_00h)
  CASE ( 'dewp_ave_00h_r' )
    CALL dict_set(defaults_dict, 'dewp_ave_00h_r', metstats_inis%dewp_ave_00h)
  CASE ( 'wind_ave_00h_r' )
    CALL dict_set(defaults_dict, 'wind_ave_00h_r', metstats_inis%wind_ave_00h)
  CASE ( 'temp_max_00h' )
    CALL dict_set(defaults_dict, 'temp_max_00h',   metstats_inis%temp_max_00h)
  CASE ( 'temp_ave_00h' )
    CALL dict_set(defaults_dict, 'temp_ave_00h',   metstats_inis%temp_ave_00h)
  CASE ( 'temp_ave_nday' )
    CALL dict_set(defaults_dict, 'temp_ave_nday',  metstats_inis%temp_ave_nday)
  CASE ( 'temp_pnt_12h' )
    CALL dict_set(defaults_dict, 'temp_pnt_12h',   metstats_inis%temp_pnt_12h)
  CASE ( 'prec_tot_00h' )
    CALL dict_set(defaults_dict, 'prec_tot_00h',   metstats_inis%prec_tot_00h)
  CASE ( 'prec_tot_12h' )
    CALL dict_set(defaults_dict, 'prec_tot_12h',   metstats_inis%prec_tot_12h)
  CASE ( 'rhum_min_00h' )
    CALL dict_set(defaults_dict, 'rhum_min_00h',   metstats_inis%rhum_min_00h)
  CASE ( 'rhum_pnt_12h' )
    CALL dict_set(defaults_dict, 'rhum_pnt_12h',   metstats_inis%rhum_pnt_12h)
  CASE ( 'dewp_ave_00h' )
    CALL dict_set(defaults_dict, 'dewp_ave_00h',   metstats_inis%dewp_ave_00h)
  CASE ( 'wind_ave_00h' )
    CALL dict_set(defaults_dict, 'wind_ave_00h',   metstats_inis%wind_ave_00h)
  CASE ( 'wind_pnt_12h' )
    CALL dict_set(defaults_dict, 'wind_pnt_12h',   metstats_inis%wind_pnt_12h)

    !Fire module variables
  CASE ( 'fire_canadian_ffmc' )
    CALL dict_set(defaults_dict, 'fire_canadian_ffmc'     , fire_inis%canadian_ffmc )
  CASE ( 'fire_canadian_ffmc_mois' )
    CALL dict_set(defaults_dict, 'fire_canadian_ffmc_mois', fire_inis%canadian_ffmc_mois)
  CASE ( 'fire_canadian_dmc' )
    CALL dict_set(defaults_dict, 'fire_canadian_dmc'      , fire_inis%canadian_dmc )
  CASE ( 'fire_canadian_dc' )
    CALL dict_set(defaults_dict, 'fire_canadian_dc'       , fire_inis%canadian_dc )
  CASE ( 'fire_mcarthur_r_dr' )
    CALL dict_set(defaults_dict, 'fire_mcarthur_r_dr'     , fire_inis%mcarthur_r_dr )
  CASE ( 'fire_mcarthur_n_dr' )
    CALL dict_set(defaults_dict, 'fire_mcarthur_n_dr'     , fire_inis%mcarthur_n_dr)
  CASE ( 'fire_nesterov' )
    CALL dict_set(defaults_dict, 'fire_nesterov'          , fire_inis%nesterov_index )

    !FLake variables
  CASE ( 'lake_fetch_gb' )
    CALL dict_set(defaults_dict, 'lake_fetch_gb'       , lake_fetch_0)
  CASE ( 'lake_h_mxl_gb' )
    CALL dict_set(defaults_dict, 'lake_h_mxl_gb'       , lake_h_mxl_0)
  CASE ( 'lake_shape_factor_gb' )
    CALL dict_set(defaults_dict, 'lake_shape_factor_gb', lake_shape_0)
  CASE ( 'lake_t_mean_gb' )
    CALL dict_set(defaults_dict, 'lake_t_mean_gb'      , lake_T_mean_0)
  CASE ( 'lake_t_mxl_gb' )
    CALL dict_set(defaults_dict, 'lake_t_mxl_gb'       , lake_T_mxl_0)
  CASE ( 'lake_t_ice_gb' )
    CALL dict_set(defaults_dict, 'lake_t_ice_gb'       , lake_T_ice_0)
  CASE ( 'lake_h_ice_gb' )
    CALL dict_set(defaults_dict, 'lake_h_ice_gb'       , lake_H_ice_0)

    !Snow variables
    !      CASE ( 'rho_snow' )
    !      CASE ( 'snow_depth' )
    !      CASE ( 'snow_grnd' )
    !      CASE ( 'nsnow' )
    !      CASE ( 'snow_ds' )
    !      CASE ( 'snow_ice' )
    !      CASE ( 'snow_liq' )
    !      CASE ( 'tsnow' )

  CASE DEFAULT
    CALL log_info("get_default_ic_values",                                    &
                   "No default value defined. Value must be defined in namelist - " &
                   // TRIM(identifiers(i)))
  END SELECT
END DO

RETURN

END FUNCTION get_default_ic_values
#endif
