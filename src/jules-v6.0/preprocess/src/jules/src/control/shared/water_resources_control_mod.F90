!******************************COPYRIGHT**************************************
! (c) UK Centre for Ecology & Hydrology.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT**************************************

MODULE water_resources_control_mod

!-----------------------------------------------------------------------------
! Description:
!   Control-level code for water resource management.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in HYDROLOGY
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE  !  Private scope by default.
PUBLIC water_resources_control

! Module parameters.

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                       &
  ModuleName = 'WATER_RESOURCES_CONTROL_MOD'

CONTAINS

!#############################################################################

SUBROUTINE water_resources_control(                                           &
             ntype, land_index, con_rain_ij, con_snow_ij,                     &
             conveyance_loss, demand_rate_domestic, demand_rate_industry,     &
             demand_rate_livestock, demand_rate_transfers, dvi_cpft,          &
             flandg, frac_irr_soilt, frac_soilt,                              &
             frac_surft, irrig_eff, grid_area_ij,                             &
             ls_rain_ij, ls_snow_ij, lw_down,                                 &
             smvccl_soilt, smvcst_soilt, sthf_soilt,                          &
             sw_surft, tl_1_ij, tstar_surft,                                  &
             icntmax_gb, plant_n_gb, demand_accum,                            &
             prec_1_day_av_gb, prec_1_day_av_use_gb,                          &
             rn_1_day_av_gb, rn_1_day_av_use_gb,                              &
             smcl_soilt, sthu_irr_soilt, sthu_soilt,                          &
             tl_1_day_av_gb, tl_1_day_av_use_gb,                              &
             priority_order, irrig_water_gb )

USE ancil_info, ONLY: land_pts, nsoilt, nsurft

USE atm_fields_bounds_mod, ONLY: pdims_s

USE crop_date_mod, ONLY: calc_crop_date

USE crop_vars_mod, ONLY: ndpy, nyav

USE jules_irrig_mod, ONLY: irr_crop_doell, irr_crop

USE jules_soil_mod, ONLY: sm_levels

USE jules_surface_types_mod, ONLY: ncpft

USE jules_water_resources_mod, ONLY:                                          &
  l_water_irrigation, nstep_water_res, nwater_use, water_res_count

USE theta_field_sizes, ONLY: row_length=>t_i_length, rows=>t_j_length

USE model_time_mod, ONLY: timestep_number=>timestep

USE water_resources_drive_mod, ONLY: water_resources_drive

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Top-level control routine for water resource management.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  ntype
    ! Number of surface types.

!-----------------------------------------------------------------------------
! Array arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  land_index(land_pts)
    ! Index of land points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  con_rain_ij(row_length,rows),                                               &
    ! Convective rainfall rate (kg m-2 s-1).
  con_snow_ij(row_length,rows),                                               &
    ! Convective snowfall rate (kg m-2 s-1).
  conveyance_loss(land_pts),                                                  &
    ! Fraction of water that is lost during conveyance from source to user.
  demand_rate_domestic(land_pts),                                             &
    ! Demand for water for domestic use (kg s-1).
  demand_rate_industry(land_pts),                                             &
    ! Demand for water for industrial use (kg s-1).
  demand_rate_livestock(land_pts),                                            &
    ! Demand for water for livestock (kg s-1).
  demand_rate_transfers(land_pts),                                            &
    ! Demand for water for (explicit) transfers (kg s-1).
  dvi_cpft(land_pts,ncpft),                                                   &
    ! Development index for crop tiles. 
  flandg(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),        &
    ! Land fraction.
  frac_irr_soilt(land_pts,nsoilt),                                            &
    ! Irrigation fraction.
  frac_soilt(land_pts,nsoilt),                                                &
    !  Fraction of gridbox for each soil tile.
  frac_surft(land_pts,ntype),                                                 &
    ! Fractions of surface types.
  irrig_eff(land_pts),                                                        &
    ! Irrigation efficiency i.e. the fraction of the water withdrawn for
    ! irrigation that is demanded by the crop scheme.
  grid_area_ij(row_length,rows),                                              &
    ! Area of each gridbox (m2).
  ls_rain_ij(row_length,rows),                                                &
    ! Large-scale rainfall rate (kg m-2 s-1).
  ls_snow_ij(row_length,rows),                                                &
    ! Large-scale snowfall rate (kg m-2 s-1).
  lw_down(row_length,rows),                                                   &
    ! Surface downward longwave radiation (W m-2).
  smvccl_soilt(land_pts,nsoilt,sm_levels),                                    &
    ! Critical volumetric SMC (cubic m per cubic m of soil).
  smvcst_soilt(land_pts,nsoilt,sm_levels),                                    &
    ! Volumetric saturation point (m3/m3 of soil).
  sthf_soilt(land_pts,nsoilt,sm_levels),                                      &
    ! Frozen soil moisture content of each layer as a fraction of saturation.
  sw_surft(land_pts,nsurft),                                                  &
    ! Surface net shortwave radiation on land tiles (W m-2).
  tl_1_ij(row_length,rows),                                                   &
    ! Ice/liquid water temperature (K).
  tstar_surft(land_pts,nsurft)
    ! Surface temperature (K).

!-----------------------------------------------------------------------------
! Array arguments with intent(in out)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN OUT) ::                                                    &
  icntmax_gb(land_pts),                                                       &
    ! Counter for start date for non-rice crops.
  plant_n_gb(land_pts)
    ! Best planting date for non-rice crops.

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                     &
  demand_accum(land_pts,nwater_use),                                          &
    ! Demands for water accumulated over the water resource timestep (kg).
  prec_1_day_av_gb(land_pts),                                                 &
    ! Average precipitation rate for the current day (kg m-2 s-1).
  prec_1_day_av_use_gb(land_pts,ndpy,nyav),                                   &
    ! Daily average precipitation rate (kg m-2 s-1).
  rn_1_day_av_gb(land_pts),                                                   &
    ! Average net radiation for the current day (W m-2).
  rn_1_day_av_use_gb(land_pts,ndpy,nyav),                                     &
    ! Daily average net radiation (W m-2).
  smcl_soilt(land_pts,nsoilt,sm_levels),                                      &
    ! Soil moisture content of each layer (kg/m2).
  sthu_irr_soilt(land_pts,nsoilt,sm_levels),                                  &
    ! Unfrozen soil moisture content of each layer as a fraction of
    ! saturation in irrigated fraction.
  sthu_soilt(land_pts,nsoilt,sm_levels),                                      &
    ! Unfrozen soil moisture content of each layer as a fraction of
    ! saturation.
  tl_1_day_av_gb(land_pts),                                                   &
    ! Average air temperature for the current day (K).
  tl_1_day_av_use_gb(land_pts,ndpy,nyav)
    ! Daily average air temperature (K).

!-----------------------------------------------------------------------------
! Array arguments with intent(out)
!-----------------------------------------------------------------------------
INTEGER, INTENT(OUT) ::                                                       &
  priority_order(land_pts,nwater_use)
    ! Priorities of water demands at each gridpoint, in order of decreasing
    ! priority. Values are the index in multi-sector arrays.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  irrig_water_gb(land_pts)
    ! Water added to soil via irrigation (kg m-2 s-1).

!-----------------------------------------------------------------------------
! Local parameters
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'WATER_RESOURCES_CONTROL'

!-----------------------------------------------------------------------------
! Local scalar variables
!-----------------------------------------------------------------------------
LOGICAL ::                                                                    &
  l_water_res_call
    ! TRUE on timesteps when the water resource model is to be run,
    ! FALSE otherwise.

! Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------------
!end of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! On the first call to this routine, initialise further variables.
!-----------------------------------------------------------------------------
IF ( timestep_number == 1 ) THEN
  CALL initialise_water_resources( priority_order )
END IF

!-----------------------------------------------------------------------------
! Work out if the water resource model is to be called this timestep.
!-----------------------------------------------------------------------------
! Increment timestep counter.
water_res_count = water_res_count + 1

IF ( water_res_count == nstep_water_res ) THEN
  l_water_res_call = .TRUE.
ELSE
  l_water_res_call = .FALSE.
END IF

!-----------------------------------------------------------------------------
! Initialise accumulations to zero.
! We do this at the start of the water resource timestep (rather than at the
! end) so the accumulations are available for potential use as diagnostics.
! Don't reinitialise rate diagnostics - so they retain the same value until
! next updated.
!-----------------------------------------------------------------------------
IF ( water_res_count == 1 ) THEN
  demand_accum(:,:) = 0.0
END IF

!-----------------------------------------------------------------------------
! Add to the accumulated demands.
!-----------------------------------------------------------------------------
CALL accumulate_demand( conveyance_loss, demand_rate_domestic,                &
                        demand_rate_industry, demand_rate_livestock,          &
                        demand_rate_transfers, demand_accum)

!-----------------------------------------------------------------------------
! Call the top-level routine for the chosen model.
!-----------------------------------------------------------------------------
IF ( l_water_res_call ) THEN

  ! Initially only one model is available.
  CALL water_resources_drive( ntype, land_index, plant_n_gb,                  &
             dvi_cpft, flandg, frac_irr_soilt, frac_soilt,                    &
             irrig_eff, grid_area_ij, smvccl_soilt, smvcst_soilt, sthf_soilt, &
             demand_accum, smcl_soilt, sthu_irr_soilt, sthu_soilt,            &
             irrig_water_gb )

  ! Reset timestep counter.
  water_res_count = 0

END IF  !  l_water_res_call

!-----------------------------------------------------------------------------
! Update crop/irrigation-related information.
! This is called every timestep.
!-----------------------------------------------------------------------------
IF ( l_water_irrigation .AND. irr_crop == irr_crop_doell ) THEN
  CALL calc_crop_date(land_index, land_pts, row_length, rows,                 &
                      nsurft, frac_surft,                                     &
                      sw_surft, tstar_surft, lw_down, tl_1_ij,                &
                      con_rain_ij, ls_rain_ij, con_snow_ij, ls_snow_ij,       &
                      prec_1_day_av_gb, prec_1_day_av_use_gb,                 &
                      rn_1_day_av_gb, rn_1_day_av_use_gb,                     &
                      tl_1_day_av_gb, tl_1_day_av_use_gb,                     &
                      icntmax_gb, plant_n_gb)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE water_resources_control

!#############################################################################
!#############################################################################

SUBROUTINE initialise_water_resources( priority_order )

!-----------------------------------------------------------------------------
! Description:
!   Initialise further aspects of the water resource management code.
!-----------------------------------------------------------------------------

USE ancil_info, ONLY: land_pts

USE ereport_mod, ONLY: ereport

USE jules_water_resources_mod, ONLY:                                          &
  l_prioritise, name_domestic, name_environment, name_industry,               &
  name_irrigation, name_livestock, name_transfers,                            &
  nwater_use,  priority, use_domestic, use_environment, use_industry,         &
  use_irrigation, use_livestock, use_transfers

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Array arguments with intent(out)
!-----------------------------------------------------------------------------
INTEGER, INTENT(OUT) ::                                                       &
  priority_order(land_pts,nwater_use)
    ! Water demands at each gridpoint, in order of decreasing priority.
    ! Values are the index in multi-sector arrays.

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'INITIALISE_WATER_RESOURCES'

!-----------------------------------------------------------------------------
! Local scalar variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  error_status,                                                               &
    ! Error status.
  i
    ! Loop counter.

! Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------------
!end of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( l_prioritise ) THEN
  ! Set sector priorities at each location.
  ! At present these are the same at all locations and it is simply a case
  ! of setting values based on the priority variable.
  ! In future this information might come from an ancillary file.
  ! The ancillary could list the sector names and use grids of numerical
  ! values [indicating the index in the name variable]. The names can be
  ! checked against those known to this code, to ensure the ancil uses a
  ! scheme that is consistent with this code.
  DO i = 1,nwater_use
    SELECT CASE ( priority(i) )
    CASE ( name_domestic )
      priority_order(:,i) = use_domestic
    CASE ( name_environment )
      priority_order(:,i) = use_environment
    CASE ( name_industry )
      priority_order(:,i) = use_industry
    CASE ( name_irrigation )
      priority_order(:,i) = use_irrigation
    CASE ( name_livestock )
      priority_order(:,i) = use_livestock
    CASE ( name_transfers )
      priority_order(:,i) = use_transfers
    CASE DEFAULT
      ! Set error status to show a fatal error.
      error_status = 101
      CALL ereport ( RoutineName, error_status,                               &
                     "Priority name not valid: " // TRIM(priority(i)) )
    END SELECT
  END DO

END IF  !  l_prioritise

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE initialise_water_resources

!#############################################################################
!#############################################################################

SUBROUTINE accumulate_demand( conveyance_loss, demand_rate_domestic,          &
             demand_rate_industry, demand_rate_livestock,                     &
             demand_rate_transfers, demand_accum )

!-----------------------------------------------------------------------------
! Add to the accumulated demand in each sector over the water resource
! timestep, converting from kg s-1 to kg, and increasing by conveyance loss
! for relevant demands.
! We only do this for demands that are prescribed rather than calculated
! within the model.
!-----------------------------------------------------------------------------

USE ancil_info, ONLY: land_pts

USE jules_water_resources_mod, ONLY:                                          &
  l_water_domestic, l_water_industry, l_water_livestock, l_water_transfers,   &
  nwater_use, use_domestic, use_industry, use_livestock, use_transfers

USE timestep_mod, ONLY: timestep_len=>timestep

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Array arguments with intent(in)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  conveyance_loss(land_pts),                                                  &
    ! Fraction of water that is lost during conveyance from source to user.
  demand_rate_domestic(land_pts),                                             &
    ! Demand for water for domestic use (kg s-1).
  demand_rate_industry(land_pts),                                             &
    ! Demand for water for industrial use (kg s-1).
  demand_rate_livestock(land_pts),                                            &
    ! Demand for water for livestock (kg s-1).
  demand_rate_transfers(land_pts)
    ! Demand for water for (explicit) transfers (kg s-1).

!-----------------------------------------------------------------------------
! Array arguments with intent(in out)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                     &
  demand_accum(land_pts,nwater_use)
    ! Demands for water accumulated over the water resource timestep (kg).

!-----------------------------------------------------------------------------
! Local parameters
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ACCUMULATE_DEMAND'

! Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------------
!end of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( l_water_domestic ) THEN
  demand_accum(:,use_domestic) = demand_accum(:,use_domestic)                 &
                                 + demand_rate_domestic(:) * timestep_len     &
                                 * ( 1.0 + conveyance_loss(:) )
END IF

IF ( l_water_industry ) THEN
  demand_accum(:,use_industry) = demand_accum(:,use_industry)                 &
                                 + demand_rate_industry(:) * timestep_len     &
                                 * ( 1.0 + conveyance_loss(:) )
END IF

IF ( l_water_livestock ) THEN
  demand_accum(:,use_livestock) = demand_accum(:,use_livestock)               &
                                  + demand_rate_livestock(:) * timestep_len   &
                                  * ( 1.0 + conveyance_loss(:) )
END IF

IF ( l_water_transfers ) THEN
  demand_accum(:,use_transfers) = demand_accum(:,use_transfers)               &
                                  + demand_rate_transfers(:) * timestep_len
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE accumulate_demand

!#############################################################################
!#############################################################################

END MODULE water_resources_control_mod
