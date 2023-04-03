MODULE surf_couple_extra_mod
! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT****************************************
USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE

PUBLIC :: surf_couple_extra

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SURF_COUPLE_EXTRA_MOD'

CONTAINS

!===============================================================================
!Note that the commented intents may not be correct. The declarations section
!will be more correct

!Routine Layout **PLEASE READ**

!-Header (please follow ordering guidance)
!-Type correction for timestep
!-SELECT for LSM
!--JULES
!---Redimensioning of arrays (Note 1)
!---Compression to land points (Note 2)
!---Science calls (Note 3)
!---Expansion to full grid
!---UM Diagnostics calls
!--CABLE
!---Details to follow
!-END

!Note 1- This is for any variable that might have differnt numbers of
!        dimensions in standalone vs UM mode. These generally arise due to the
!        differing science payloads.

!Note 2- This is where all gridded variables (with _ij suffix) are compressed
!        onto land points. This allows efficient use of OpenMP

!Note 3- New additions here should ideally consist ONLY of calls to other
!        routines. Please encapsulate longer code within a subroutine.

!===============================================================================
SUBROUTINE surf_couple_extra(                                                  &
   ! Arguments used by JULES-standalone
   !Driving data and associated INTENT(IN)
   ls_rain_ij, con_rain_ij, ls_snow_ij, con_snow_ij, tl_1_ij, lw_down,         &
   qw_1_ij, u_1_ij, v_1_ij,                                                    &
   pstar_ij,                                                                   &
   !Fluxes INTENT(IN)
   ei_surft, surf_htf_surft, ecan_surft, ext_soilt, sw_surft,                  &
   !Misc INTENT(IN)
   a_step, smlt, tile_frac, hcons_soilt, rhostar,                              &
   !Fluxes INTENT(INOUT)
   melt_surft,                                                                 &
   !Misc INTENT(INOUT)
   lake_albedo_gb, lake_t_sfc_gb, lake_h_snow_gb, lake_t_snow_gb,              &
   !Fluxes INTENT(OUT)
   snomlt_surf_htf, snowmelt_ij, snomlt_sub_htf, sub_surf_roff, surf_roff,     &
   tot_tfall, snowmelt_gb, rrun, rflow, snow_soil_htf,                         &
   !Arguments for the UM-----------------------------------------
   !IN
   land_pts, row_length, rows, river_row_length, river_rows,                   &
   ls_graup_ij,                                                                &
   cca_2d, nsurft, surft_pts,                                                  &
   lice_pts, soil_pts,                                                         &
   stf_sub_surf_roff,                                                          &
   fexp_soilt, gamtot_soilt, ti_mean_soilt, ti_sig_soilt,                      &
   cs_ch4_soilt, flash_rate_ancil, pop_den_ancil,                              &
   a_fsat_soilt, c_fsat_soilt, a_fwet_soilt, c_fwet_soilt,                     &
   ntype, fqw_surft,                                                           &
   delta_lambda, delta_phi, xx_cos_theta_latitude,                             &
   aocpl_row_length, aocpl_p_rows, xpa, xua, xva, ypa, yua, yva,               &
   g_p_field, g_r_field, n_proc, global_row_length, global_rows,               &
   global_river_row_length, global_river_rows, flandg,                         &
   trivdir, trivseq, r_area, slope, flowobs1, r_inext, r_jnext, r_land,        &
   frac_agr_gb, soil_clay_ij, resp_s_soilt, npp_gb,                            &
   u_s_std_surft,                                                              &
   !INOUT
   a_steps_since_riv,                                                          &
   fsat_soilt, fwetl_soilt,                                                    &
   zw_soilt, sthzw_soilt,                                                      &
   ls_rainfrac_gb,                                                             &
   substore, surfstore, flowin, bflowin,                                       &
   tot_surf_runoff_gb, tot_sub_runoff_gb, acc_lake_evap_gb, twatstor,          &
   asteps_since_triffid,                                                       &
   resp_s_acc_gb_um, cs_pool_gb_um,                                            &
   inlandout_atm_gb,                                                           &
   !OUT
   dhf_surf_minus_soil,                                                        &
   land_sea_mask,                                                              &
   !TYPES containing field data (IN OUT)
   crop_vars,psparms,toppdm,fire_vars,ainfo,trif_vars,soilecosse,urban_param,  &
   progs,trifctltype,jules_vars)

!Module imports

!TYPE definitions
USE crop_vars_mod, ONLY: crop_vars_type
USE p_s_parms, ONLY: psparms_type
USE top_pdm, ONLY: top_pdm_type
USE fire_vars_mod, ONLY: fire_vars_type
USE ancil_info,    ONLY: ainfo_type
USE trif_vars_mod, ONLY: trif_vars_type
USE soil_ecosse_vars_mod, ONLY: soil_ecosse_vars_type
USE urban_param_mod, ONLY:urban_param_type
USE prognostics, ONLY: progs_type
USE trifctl, ONLY: trifctl_type
USE jules_vars_mod, ONLY: jules_vars_type

!Import interfaces to subroutines called
USE hydrol_mod,               ONLY: hydrol
USE snow_mod,                 ONLY: snow
USE jules_rivers_mod,         ONLY: l_rivers, l_inland, rivers_call

! Code which isn't currently suitable for building into LFRic
USE river_control_mod,        ONLY: river_control
USE inferno_mod,              ONLY: calc_soil_carbon_pools
USE inferno_io_mod,           ONLY: inferno_io
USE irrigation_mod,           ONLY: irrigation_control
USE veg_control_mod,          ONLY: veg_control
USE next_gen_biogeochem_mod,  ONLY: next_gen_biogeochem

!ifdef required to manage the different science payloads available for
!coupled and standalone use
USE crop_mod,                 ONLY: crop
USE fire_timestep_mod,        ONLY: fire_timestep
USE metstats_timestep_mod,    ONLY: metstats_timestep
USE gridbox_mean_mod,         ONLY: soiltiles_to_gbm, surftiles_to_gbm
USE fao_evapotranspiration,   ONLY: fao_ref_evapotranspiration
USE water_resources_control_mod, ONLY: water_resources_control
USE zenith_mod,               ONLY: photoperiod
! End of standalone Jules code

! End of code excluded from LFRic builds

!Variables- The rules here are:
!Each module is only USED once, using an ifdef to control the variables if req.
!Whole modules are grouped into a single ifdef/else
!Alphabetical order

!Modules common to JULES and UM

USE ancil_info,               ONLY:                                            &
                                    dim_cs1,                                   &
                                    dim_cslayer, nsoilt, nmasst,               &
                                  ! Replaces USE statement
                                    dim_soil_n_pool

USE atm_fields_bounds_mod,    ONLY: tdims, tdims_s, pdims, pdims_s

USE jules_hydrology_mod,      ONLY: l_hydrology, l_pdm, l_top, l_var_rainfrac

USE lake_mod, ONLY:  u_s_lake_gb                                               &
                     , surf_ht_flux_lake_ij                                    &
                     , surf_ht_flux_lk_gb                                      &
                     , sw_down_gb                                              &
                     , lake_depth_gb                                           &
                     , lake_fetch_gb                                           &
                     , coriolis_param_gb                                       &
                     , lake_t_ice_gb                                           &
                     , lake_t_mean_gb                                          &
                     , lake_t_mxl_gb                                           &
                     , lake_shape_factor_gb                                    &
                     , lake_h_ice_gb                                           &
                     , lake_h_mxl_gb                                           &
                     , ts1_lake_gb                                             &
                     , g_dt_gb                                                 &
                     , h_snow_sw_att                                           &
                     , trap_frozen                                             &
                     , trap_unfrozen

USE jules_snow_mod,      ONLY: nsmax

USE jules_soil_biogeochem_mod, ONLY:                                           &
  soil_model_1pool, soil_model_ecosse, soil_model_rothc, soil_bgc_model

USE jules_soil_mod,           ONLY: sm_levels, l_soil_sat_down, confrac

USE jules_surface_mod,        ONLY: l_aggregate, l_flake_model

USE jules_surface_types_mod,  ONLY: npft, ncpft, nnpft, soil, lake

USE jules_vegetation_mod,     ONLY:                                            &
  l_crop, l_triffid, l_trif_eq, l_phenol, phenol_period, triffid_period,       &
  l_inferno, ignition_method, l_fao_ref_evapotranspiration, frac_min

USE jules_irrig_mod, ONLY: l_irrig_dmd

USE jules_water_resources_mod, ONLY: irrig_eff, l_water_resources

USE sf_diags_mod,             ONLY: sf_diag

USE theta_field_sizes,        ONLY: t_i_length, t_j_length

USE veg3_parm_mod,            ONLY:                                            &
  l_veg3,veg3_ctrl,litter_parms,red_parms

USE veg3_field_mod,           ONLY: veg_state,red_state

USE water_constants_mod,      ONLY: rho_water

USE ereport_mod,              ONLY: ereport

!Modules specific to the UM
!Modules specific to JULES
USE conversions_mod,          ONLY: isec_per_day, rsec_per_day

USE fire_mod,                 ONLY: fire_prog, fire_diag, l_fire

USE fluxes,                   ONLY: surf_ht_flux_ij

USE forcing,                  ONLY: sw_down_ij, lw_down_ij

USE metstats_mod,             ONLY: metstats_prog, l_metstats

USE model_grid_mod,           ONLY: grid_area_ij

USE model_time_mod,           ONLY: timestep_len, current_time

USE jules_rivers_mod, ONLY: &
  rivers_sto_per_m2_on_landpts, rivers_adj_on_landpts

USE jules_water_resources_mod, ONLY:                                          &
  conveyance_loss, demand_accum, demand_rate_domestic, demand_rate_industry,  &
  demand_rate_livestock, demand_rate_transfers, priority_order

!Technical modules
USE parkind1,                    ONLY: jprb, jpim
USE yomhook,                     ONLY: lhook, dr_hook
USE ereport_mod,                 ONLY: ereport
USE jules_print_mgr,             ONLY: jules_message, jules_print
USE jules_model_environment_mod, ONLY: lsm_id, jules, cable

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Private subroutine for accessing the JULES science routines that are called
!   after the implicit code
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

!Subroutine Arguments, ordered by intent and type

INTEGER, INTENT(IN) ::                                                         &
  land_pts,                                                                    &
  row_length,                                                                  &
  rows,                                                                        &
  river_row_length,                                                            &
  river_rows,                                                                  &
  nsurft,                                                                      &
  a_step,                                                                      &
  surft_pts(nsurft),                                                           &
  lice_pts,                                                                    &
  soil_pts,                                                                    &
  ntype,                                                                       &
  aocpl_row_length,                                                            &
  aocpl_p_rows,                                                                &
  g_p_field,                                                                   &
  g_r_field,                                                                   &
  n_proc,                                                                      &
  global_row_length,                                                           &
  global_rows,                                                                 &
  global_river_row_length,                                                     &
  global_river_rows

LOGICAL, INTENT(IN) ::                                                         &
  smlt,                                                                        &
  stf_sub_surf_roff,                                                           &
  land_sea_mask(row_length, rows)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  !Driving data
  ls_rain_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
  con_rain_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
  ls_snow_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
  con_snow_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
  pstar_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
  ls_graup_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
  tl_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
  lw_down(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
  qw_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
  u_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
  v_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &

  !Fluxes
  ei_surft(land_pts,nsurft),                                                   &
       !Sublimation of snow (kg/m2/s)
  surf_htf_surft(land_pts,nsurft),                                             &
       !Surface heat flux (W/m2)
  ecan_surft(land_pts,nsurft),                                                 &
       !Canopy evaporation from land tiles (kg/m2/s).
  ext_soilt(land_pts,nsoilt,sm_levels),                                        &
       !Extraction of water from each soil layer (kg/m2/s).
  tile_frac(land_pts,nsurft),                                                  &
  rhostar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
       !Surface air density
  sw_surft(land_pts,nsurft),                                                   &
       !Surface net SW radiation on land tiles (W/m2)
  cca_2d(row_length,rows),                                                     &
  fexp_soilt(land_pts,nsoilt),                                                 &
  gamtot_soilt(land_pts,nsoilt),                                               &
  ti_mean_soilt(land_pts,nsoilt),                                              &
  ti_sig_soilt(land_pts,nsoilt),                                               &
  npp_gb(land_pts),                                                            &
  a_fsat_soilt(land_pts,nsoilt),                                               &
  c_fsat_soilt(land_pts,nsoilt),                                               &
  a_fwet_soilt(land_pts,nsoilt),                                               &
  c_fwet_soilt(land_pts,nsoilt),                                               &
  fqw_surft(land_pts,nsurft),                                                  &
  frac_agr_gb(land_pts),                                                       &
  soil_clay_ij(row_length,rows),                                               &
  flash_rate_ancil(row_length,rows),                                           &
  pop_den_ancil(row_length,rows),                                              &
  u_s_std_surft(land_pts, nsurft),                                             &

  !River routing
  delta_lambda,                                                                &
  delta_phi,                                                                   &
  xx_cos_theta_latitude(tdims_s%i_start:tdims_s%i_end,                         &
                        tdims_s%j_start:tdims_s%j_end),                        &
  xpa(aocpl_row_length+1),                                                     &
  xua(0:aocpl_row_length),                                                     &
  xva(aocpl_row_length+1),                                                     &
  ypa(aocpl_p_rows),                                                           &
  yua(aocpl_p_rows),                                                           &
  yva(0:aocpl_p_rows),                                                         &
  flandg(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),         &
  trivdir(river_row_length, river_rows),                                       &
  trivseq(river_row_length, river_rows),                                       &
  r_area(row_length, rows),                                                    &
  slope(row_length, rows),                                                     &
  flowobs1(row_length, rows),                                                  &
  r_inext(row_length, rows),                                                   &
  r_jnext(row_length, rows),                                                   &
  r_land(row_length, rows)

INTEGER, INTENT(IN OUT)  ::                                                    &
  a_steps_since_riv,                                                           &
  asteps_since_triffid

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  melt_surft(land_pts,nsurft),                                                 &
        !Surface or canopy snowmelt rate (kg/m2/s)
        !On output, this is the total melt rate for the tile
        !(i.e. sum of  melt on canopy and ground).
  snomlt_sub_htf(land_pts),                                                    &
        !Sub-canopy snowmelt heat flux (W/m2)
  hcons_soilt(land_pts,nsoilt),                                                &
  lake_albedo_gb(land_pts),                                                    &
        ! lake albedo
  lake_t_sfc_gb(land_pts),                                                     &
        !temperature (of water, ice or snow) at surface of lake (K)
  lake_h_snow_gb(land_pts),                                                    &
        !snow thickness (m)
  lake_t_snow_gb(land_pts),                                                    &
        !temperature at the air-snow interface (K)
  dhf_surf_minus_soil(land_pts),                                               &
        ! Heat flux difference across the FLake snowpack (W/m2)
  resp_s_soilt(land_pts,nsoilt,1,dim_cs1),                                     &
  ls_rainfrac_gb(land_pts),                                                    &
  fsat_soilt(land_pts,nsoilt),                                                 &
  fwetl_soilt(land_pts,nsoilt),                                                &
  zw_soilt(land_pts,nsoilt),                                                   &
  sthzw_soilt(land_pts,nsoilt),                                                &
  substore(row_length, rows),                                                  &
  surfstore(row_length, rows),                                                 &
  flowin(row_length, rows),                                                    &
  bflowin(row_length, rows),                                                   &
  tot_surf_runoff_gb(land_pts),                                                &
  tot_sub_runoff_gb(land_pts),                                                 &
  acc_lake_evap_gb(row_length,rows),                                           &
  twatstor(river_row_length, river_rows),                                      &
  resp_s_acc_gb_um(land_pts,dim_cs1),                                          &
  cs_ch4_soilt(land_pts,nsoilt),                                               &
         ! soil carbon used in wetland CH4 emissions model if TRIFFID
         ! is switched off
  cs_pool_gb_um(land_pts,dim_cs1),                                             &
  inlandout_atm_gb(land_pts)

REAL(KIND=real_jlslsm), INTENT(OUT)::                                          &
  sub_surf_roff(land_pts),                                                     &
        !Sub-surface runoff (kg/m2/s).
  surf_roff(land_pts),                                                         &
        !Surface runoff (kg/m2/s).
  tot_tfall(land_pts),                                                         &
        !Total throughfall (kg/m2/s).
  snomlt_surf_htf(row_length,rows),                                            &
        !Gridbox snowmelt heat flux (W/m2)
  snow_soil_htf(land_pts,nsurft),                                              &
        !Tiled snow->soil heat flux (W/m2)
  snowmelt_gb(land_pts),                                                       &
  snowmelt_ij(row_length,rows),                                                &
  rrun(land_pts),                                                              &
         !Surface runoff after river routing (kg/m2/s)
  rflow(land_pts)
         !River runoff (kg/m2/s)

!TYPES containing field data (IN OUT)
TYPE(crop_vars_type), INTENT(IN OUT) :: crop_vars
TYPE(psparms_type), INTENT(IN OUT) :: psparms
TYPE(top_pdm_type), INTENT(IN OUT) :: toppdm
TYPE(fire_vars_type), INTENT(IN OUT) :: fire_vars
TYPE(ainfo_type), INTENT(IN OUT) :: ainfo
TYPE(trif_vars_type), INTENT(IN OUT) :: trif_vars
TYPE(soil_ecosse_vars_type), INTENT(IN OUT) :: soilecosse
TYPE(urban_param_type), INTENT(IN OUT) :: urban_param
TYPE(progs_type), INTENT(IN OUT) :: progs
TYPE(trifctl_type), INTENT(IN OUT) :: trifctltype
TYPE(jules_vars_type), INTENT(IN OUT) :: jules_vars

!==============================================================================
!Local variables

INTEGER ::                                                                     &
   i,j,k,l,n,m,                                                                &
        !Various counters
  errcode,                                                                     &
        ! error code to pass to ereport.
  crop_call
        !indicates whether crop model is to be called

INTEGER, PARAMETER :: crop_period = 1
        ! Crop code hard wired to run daily : crop_period = 1

REAL(KIND=real_jlslsm) :: timestep_real ! Model timestep (s) in REAL

REAL(KIND=real_jlslsm) ::                                                      &
  !Gridbox versions of forcing data (ls_rainfrac_gb is an argument)
  con_rainfrac_gb(land_pts),                                                   &
  ls_rain_gb(land_pts),                                                        &
  con_rain_gb(land_pts),                                                       &
  ls_snow_gb(land_pts),                                                        &
  ls_graup_gb(land_pts),                                                       &
  con_snow_gb(land_pts),                                                       &
  pstar_gb(land_pts),                                                          &
  tl_1_gb(land_pts),                                                           &
  qw_1_gb(land_pts),                                                           &
  u_1_gb(land_pts),                                                            &
  v_1_gb(land_pts),                                                            &
  !Passed between different science schemes (various)
  qbase_l_soilt(land_pts,nsoilt,sm_levels+1),                                  &
        ! Base flow from each soil layer (kg m-2 s-1).
  w_flux_soilt(land_pts,nsoilt,0:sm_levels),                                   &
        ! Fluxes of water between layers (kg m-2 s-1).
  surf_ht_flux_ld(land_pts),                                                   &
  dwsw_sub_snow,                                                               &
        !remaining SW flux under the snowpack (W/m2), for FLake with ML-snow.

  !Passed between fire & INFERNO routines only
  smc_gb(land_pts),                                                            &
         !To allow GBM soil moisture to be passed down to fire
  c_soil_dpm_gb(land_pts),                                                     &
         ! Gridbox soil C in the Decomposable Plant Material pool (kg m-2).
  c_soil_rpm_gb(land_pts),                                                     &
         ! Gridbox soil C in the Resistant Plant Material pool (kg m-2).

  !Passed between evapotranspiration routines only
  trad(land_pts),                                                              &
         ! gridbox effective radiative temperature (assuming emissivity=1)

  !Passed from snow or hydrol to diagnostics_hyd only
  snow_mass_gb(land_pts),                                                      &
  dun_roff_soilt(land_pts,nsoilt),                                             &
  drain_soilt(land_pts,nsoilt),                                                &

  !Passed between veg routines and diagnostics_veg only
  resp_s_dr_out_gb_um(land_pts,dim_cs1+1),                                     &

  !Passed between river routing and diagnostics_riv only
  riverout(row_length, rows),                                                  &
  riverout_rgrid(river_row_length, river_rows),                                &
  box_outflow(river_row_length, river_rows),                                   &
  box_inflow(river_row_length, river_rows),                                    &
  inlandout_riv(river_row_length,river_rows)

INTEGER                      :: errorstatus
CHARACTER(LEN=*), PARAMETER  :: RoutineName = 'SURF_COUPLE_EXTRA'

REAL :: diff

!Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!===============================================================================
!End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!Care needed as the timestep is used in arithmetic to control the
!timing of when various science schemes are called.
timestep_real = REAL(timestep_len)

!CABLE_LSM: implement switching based on lsm_id
SELECT CASE( lsm_id )

  !=============================================================================
CASE ( jules )

  !=============================================================================
  ! Redimensioning of arrays


  !Encompassing IF statement for land_pts > 0. We exit/re-enter this IF about
  !half way down to allow for river routing

  IF (l_hydrology .AND. land_pts > 0 ) THEN

    !===========================================================================
    ! Compression to land points

    DO l = 1, land_pts
      j = (ainfo%land_index(l) - 1) / row_length + 1
      i = ainfo%land_index(l) - (j-1) * row_length
      ls_rain_gb(l)    = ls_rain_ij(i,j)
      con_rain_gb(l)   = con_rain_ij(i,j)
      con_snow_gb(l)   = con_snow_ij(i,j)
      ls_snow_gb(l)    = ls_snow_ij(i,j)
      ls_graup_gb(l)   = ls_graup_ij(i,j)
      pstar_gb(l)      = pstar_ij(i,j)
      tl_1_gb(l)       = tl_1_ij(i,j)
      qw_1_gb(l)       = qw_1_ij(i,j)
      u_1_gb(l)        = u_1_ij(i,j)
      v_1_gb(l)        = v_1_ij(i,j)
    END DO

    !Pass jules the modelled rain fractions
    !In standalone mode, both rainfracs are passed in as zeroed arrays
    IF (l_var_rainfrac) THEN
      DO l = 1, land_pts
        j = (ainfo%land_index(l) - 1) / row_length + 1
        i = ainfo%land_index(l) - (j-1) * row_length

        con_rainfrac_gb(l) = MIN(cca_2d(i,j),0.5) !As in diagnostics_conv

        !provide some safety checking for convective rain with no CCA
        IF (con_rain_gb(l) > 0.0 ) THEN
          IF (con_rainfrac_gb(l) == 0.0) THEN
            con_rainfrac_gb(l) = confrac
            !and for very small CCA amounts
          ELSE IF (con_rainfrac_gb(l) < 0.01) THEN
            con_rainfrac_gb(l) = 0.01
          END IF
        END IF

        !provide some safety checking for ls rain with no rainfrac
        IF (ls_rain_gb(l) > 0.0) THEN
          IF (ls_rainfrac_gb(l) == 0.0) THEN
            ls_rainfrac_gb(l) = 0.5
            !and for very small rainfrac amounts
          ELSE IF (ls_rainfrac_gb(l) < 0.01) THEN
            ls_rainfrac_gb(l) = 0.01
          END IF
        END IF
      END DO

    ELSE !use original default values
      DO l = 1, land_pts
        con_rainfrac_gb(l) = confrac
        ls_rainfrac_gb(l)  = 1.0
      END DO
    END IF !l_var_rainfrac

    !===========================================================================
    ! Science calls

    !Snow (standalone and UM)
    CALL snow ( land_pts,timestep_real,smlt,nsurft,surft_pts,                 &
                ainfo%surft_index,psparms%catch_snow_surft,con_snow_gb,       &
                con_rain_gb, tile_frac,ls_snow_gb,ls_graup_gb,ls_rain_gb,     &
                ei_surft,psparms%hcap_soilt(:,:,1),hcons_soilt,melt_surft,    &
                progs%smcl_soilt(:,:,1),psparms%sthf_soilt(:,:,1),            &
                surf_htf_surft, progs%t_soil_soilt(:,:,1),                    &
                progs%tsurf_elev_surft, progs%tstar_surft,                    &
                psparms%smvcst_soilt(:,:,1),progs%rgrain_surft,               &
                progs%rgrainl_surft, progs%rho_snow_grnd_surft,               &
                progs%sice_surft,progs%sliq_surft,progs%snow_grnd_surft,      &
                progs%snow_surft, progs%snowdepth_surft, progs%tsnow_surft,   &
                progs%nsnow_surft,progs%ds_surft,                             &
                snomlt_surf_htf,snow_mass_gb,progs%rho_snow_surft,            &
                snomlt_sub_htf, snowmelt_gb,snow_soil_htf,surf_ht_flux_ld,    &
                sf_diag, dhf_surf_minus_soil,                                 &
                ! New Arguments to replace USE statements
                ! jules_internal
                jules_vars%unload_backgrnd_pft, npft,                         &
                !Ancil info (IN)
                ainfo%l_lice_point)

    !Hydrology (standalone and UM)
    CALL hydrol (                                                             &
      lice_pts,ainfo%lice_index,soil_pts,ainfo%soil_index, progs%nsnow_surft, &
      land_pts,sm_levels,psparms%bexp_soilt,psparms%catch_surft,con_rain_gb,  &
      ecan_surft,ext_soilt,psparms%hcap_soilt,psparms%hcon_soilt,ls_rain_gb,  &
      con_rainfrac_gb, ls_rainfrac_gb,                                        &
      psparms%satcon_soilt,psparms%sathh_soilt,progs%snowdepth_surft,         &
      snow_soil_htf, surf_ht_flux_ld,timestep_real,                           &
      psparms%smvcst_soilt,psparms%smvcwt_soilt,progs%canopy_surft,           &
      stf_sub_surf_roff,progs%smcl_soilt,psparms%sthf_soilt,                  &
      psparms%sthu_soilt, progs%t_soil_soilt,progs%tsurf_elev_surft,          &
      progs%canopy_gb,progs%smc_soilt,snowmelt_gb,                            &
      sub_surf_roff,surf_roff,tot_tfall,                                      &
      ! add new inland basin variable
      inlandout_atm_gb,l_inland,                                              &
      ! Additional variables for MOSES II
      nsurft,surft_pts,ainfo%surft_index,                                     &
      psparms%infil_surft, melt_surft,tile_frac,                              &
      ! Additional variables required for large-scale hydrology:
      l_top,l_pdm,fexp_soilt,ti_mean_soilt,cs_ch4_soilt,progs%cs_pool_soilt,  &
      dun_roff_soilt,drain_soilt,fsat_soilt,fwetl_soilt,toppdm%qbase_soilt,   &
      qbase_l_soilt, toppdm%qbase_zw_soilt, w_flux_soilt,                     &
      zw_soilt,sthzw_soilt,a_fsat_soilt,c_fsat_soilt,a_fwet_soilt,            &
      c_fwet_soilt,                                                           &
      resp_s_soilt,npp_gb,toppdm%fch4_wetl_soilt,                             &
      toppdm%fch4_wetl_cs_soilt,toppdm%fch4_wetl_npp_soilt,                   &
      toppdm%fch4_wetl_resps_soilt, crop_vars%sthu_irr_soilt,                 &
      crop_vars%frac_irr_soilt, crop_vars%ext_irr_soilt,                      &
      ! New arguments replacing USE statements
      ! trif_vars_mod
      trif_vars%n_leach_soilt, trif_vars%n_leach_gb_acc,                      &
      ! prognostics
      progs%n_inorg_soilt_lyrs, progs%n_inorg_avail_pft, npft, &
      progs%t_soil_soilt_acc, progs%tsoil_deep_gb,                            &
      ! top_pdm_alloc
       toppdm%fch4_wetl_acc_soilt,                                            &
      ! microbial methane scheme
      progs%substr_ch4, progs%mic_ch4, progs%mic_act_ch4, progs%acclim_ch4,   &
      ! pdm_vars
       toppdm%slope_gb, dim_cs1, l_soil_sat_down, asteps_since_triffid)

  END IF ! ( l_hydrology .AND. land_pts /= 0 )

  ! Code not yet ported to LFRic

    !Here we need to exit the land_pts IF to allow river routing to be called
    !on all MPI ranks. This is because it is possible that a rank with
    !land_pts = 0 still has a river routing point.

    ! Water resources (standalone; not yet allowed in UM).
    IF ( l_water_resources ) THEN
      CALL water_resources_control(                                           &
             ntype, ainfo%land_index, con_rain_ij, con_snow_ij,               &
             conveyance_loss, demand_rate_domestic, demand_rate_industry,     &
             demand_rate_livestock, demand_rate_transfers, crop_vars%dvi_cpft,&
             flandg, crop_vars%frac_irr_soilt, ainfo%frac_soilt,              &
             ainfo%frac_surft, irrig_eff, grid_area_ij,                       &
             ls_rain_ij, ls_snow_ij, lw_down,                                 &
             psparms%smvccl_soilt, psparms%smvcst_soilt, psparms%sthf_soilt,  &
             sw_surft, tl_1_ij, progs%tstar_surft,                            &
             crop_vars%icntmax_gb, crop_vars%plant_n_gb, demand_accum,        &
             crop_vars%prec_1_day_av_gb, crop_vars%prec_1_day_av_use_gb,      &
             crop_vars%rn_1_day_av_gb, crop_vars%rn_1_day_av_use_gb,          &
             progs%smcl_soilt, crop_vars%sthu_irr_soilt, psparms%sthu_soilt,  &
             crop_vars%tl_1_day_av_gb, crop_vars%tl_1_day_av_use_gb,          &
             priority_order, crop_vars%irrig_water_gb )
    END IF

    ! River Routing (standalone and UM, but differing options)
  IF ( l_rivers ) THEN
    CALL river_control(                                                       &
      !INTEGER, INTENT(IN)
      land_pts,                                                               &
      !REAL, INTENT(IN)
      sub_surf_roff, surf_roff,                                               &
      !INTEGER, INTENT(INOUT)
      a_steps_since_riv,                                                      &
      !REAL, INTENT (INOUT)
      tot_surf_runoff_gb, tot_sub_runoff_gb, acc_lake_evap_gb,                &
      !REAL, INTENT (OUT)
      rivers_sto_per_m2_on_landpts, rflow, rrun,                              &
      !Arguments for the UM-----------------------------------------
      !INTEGER, INTENT(IN)
      n_proc, row_length, rows, river_row_length, river_rows,                 &
      ainfo%land_index, ntype, aocpl_row_length, aocpl_p_rows, g_p_field,     &
      g_r_field, global_row_length, global_rows, global_river_row_length,     &
       global_river_rows, nsurft,                                             &
      !REAL, INTENT(IN)
      fqw_surft, delta_lambda, delta_phi, xx_cos_theta_latitude,              &
      xpa, xua, xva, ypa, yua, yva, flandg, trivdir,                          &
      trivseq, r_area, slope, flowobs1, r_inext, r_jnext, r_land,             &
      psparms%smvcst_soilt, psparms%smvcwt_soilt, ainfo%frac_surft,           &
      !REAL, INTENT(INOUT)
      substore, surfstore, flowin, bflowin, twatstor,                         &
      progs%smcl_soilt, psparms%sthu_soilt,                                   &
      !REAL, INTENT(OUT)
      inlandout_atm_gb, inlandout_riv, riverout, box_outflow, box_inflow,     &
      riverout_rgrid                                                          &
      )
  END IF ! l_rivers (ATMOS)

  !Irrigation (standalone and UM for some options).
  IF ( l_irrig_dmd .AND. .NOT. l_water_resources ) THEN
    CALL irrigation_control( a_step, land_pts, ntype,                         &
                       ainfo%land_index,                                      &
                       con_rain_ij, con_snow_ij, crop_vars%dvi_cpft,          &
                       crop_vars%frac_irr_soilt, ainfo%frac_surft,            &
                       ls_rain_ij, ls_snow_ij, lw_down,                       &
                       psparms%smvccl_soilt, psparms%smvcst_soilt,            &
                       psparms%smvcwt_soilt, psparms%sthf_soilt,              &
                       sw_surft, tl_1_ij, progs%tstar_surft,                  &
                       crop_vars%icntmax_gb, crop_vars%plant_n_gb,            &
                       crop_vars%irrDaysDiag_gb, crop_vars%prec_1_day_av_gb,  &
                       crop_vars%prec_1_day_av_use_gb,                        &
                       crop_vars%rn_1_day_av_gb, crop_vars%rn_1_day_av_use_gb,&
                       crop_vars%tl_1_day_av_gb, crop_vars%tl_1_day_av_use_gb,&
                       progs%smcl_soilt,  crop_vars%sthu_irr_soilt,           &
                       psparms%sthu_soilt, sthzw_soilt,                       &
                       crop_vars%irrig_water_gb,                              &
                       !New arguments replacing USE statements
                       !jules_rivers_mod
                       rivers_sto_per_m2_on_landpts, rivers_adj_on_landpts )
  ELSE IF ( .NOT. l_water_resources ) THEN
    ! Set sthu_irr_soilt to 0.0 in case it is still reported
    ! It would be better to not allocate this array when it is not being used.
    ! This led to a memory leak in the D1 array (UM).
    crop_vars%sthu_irr_soilt(:,:,:) = 0.0
  END IF ! l_irrig_dmd

  !Restart the IF for land_pts > 0 now that river routing is done with
  IF (land_pts > 0) THEN

    !Crops (standalone only)
    IF ( l_crop ) THEN
      crop_call = MOD ( REAL(a_step),                                          &
                        REAL(crop_period) * rsec_per_day / timestep_real )

      DO n = 1,ncpft
        DO l = 1,land_pts
          trifctltype%npp_acc_pft(l,nnpft + n) =                              &
                        trifctltype%npp_acc_pft(l,nnpft + n)                  &
                        + (trifctltype%npp_pft(l,nnpft + n) * timestep_real)
        END DO
      END DO

      CALL photoperiod(t_i_length * t_j_length, crop_vars%phot, &
                       crop_vars%dphotdt)

      CALL crop(t_i_length * t_j_length, land_pts, ainfo%land_index, a_step,  &
                crop_call, sm_levels, ainfo%frac_surft, crop_vars%phot,       &
                crop_vars%dphotdt,                                            &
                sf_diag%t1p5m_surft, progs%t_soil_soilt, psparms%sthu_soilt,  &
                psparms%smvccl_soilt,                                         &
                psparms%smvcst_soilt, trifctltype%npp_acc_pft,                &
                progs%canht_pft, progs%lai_pft, crop_vars%dvi_cpft,           &
                crop_vars%rootc_cpft,                                         &
                crop_vars%harvc_cpft,                                         &
                crop_vars%reservec_cpft, crop_vars%croplai_cpft,              &
                crop_vars%cropcanht_cpft,                                     &
                psparms%catch_surft, psparms%z0_surft,                        &
                !New arguments replacing USE statements
                !crop_vars_mod
                crop_vars%sow_date_cpft, crop_vars%tt_veg_cpft,               &
                crop_vars%tt_rep_cpft,                                        &
                crop_vars%latestharv_date_cpft,                               &
                crop_vars%yield_diag_cpft, crop_vars%stemc_diag_cpft,         &
                crop_vars%leafc_diag_cpft, crop_vars%nonyield_diag_cpft,      &
                crop_vars%harvest_trigger_cpft, crop_vars%harvest_counter_cpft,&
                !ancil_info (IN)
                ainfo%l_lice_point)

    END IF  ! l_crop

    !Metstats (standalone only)
        !Beware- does not currently account for graupel
    IF ( l_metstats ) THEN
      CALL metstats_timestep(tl_1_gb, qw_1_gb, u_1_gb, v_1_gb, ls_rain_gb,     &
                            con_rain_gb, ls_snow_gb, con_snow_gb, pstar_gb,    &
                            metstats_prog,                                     &
                            !Vars that should be USED but can't due to
                            !UM/standalone differences
                            current_time%time, timestep_real,land_pts)
    END IF

    !Fire (standalone only)
    IF ( l_fire ) THEN
      !Calculate the gridbox mean soil moisture
      smc_gb = soiltiles_to_gbm(progs%smc_soilt, ainfo)
      CALL fire_timestep(metstats_prog, smc_gb, fire_prog, fire_diag,          &
                        !Vars that should be USED but can't due to
                        !UM/standalone differences
                        current_time%time, current_time%month, timestep_real,  &
                        land_pts)
    END IF

    !INFERNO (standalone and UM)
    IF ( l_inferno ) THEN
      CALL calc_soil_carbon_pools(land_pts, soil_pts, ainfo%soil_index, &
                                  dim_cs1, progs%cs_pool_soilt,               &
                                  c_soil_dpm_gb, c_soil_rpm_gb)

      CALL inferno_io( sf_diag%t1p5m_surft, sf_diag%q1p5m_surft, pstar_gb,    &
                      psparms%sthu_soilt, sm_levels,                          &
                      ainfo%frac_surft, c_soil_dpm_gb, c_soil_rpm_gb,         &
                      progs%canht_pft,                                        &
                      ls_rain_gb, con_rain_gb,                                &
                      fire_vars,                                              &
                      land_pts, ignition_method,                              &
                      nsurft, asteps_since_triffid,                           &
                    ! New Arguments to replace USE statements
                    ! TRIF_VARS_MOD
                      trif_vars%g_burn_pft_acc)
    END IF  !  l_inferno

    IF (l_phenol .OR. l_triffid) THEN
      IF (l_veg3) THEN

        !Update veg3 structures
        !Structures can then be passed straight through
        ! Copy prognostics and forcings in and then out
        ! Same calling when coulpled to atmos
        DO n = 1, nnpft
          DO k = 1, land_pts
            l = ainfo%land_index(k)
            veg_state%npp_acc(l,n) = trifctltype%npp_acc_pft(l,n)
            veg_state%frac(l,n)    = ainfo%frac_surft(l,n)
            veg_state%canht(l,n)   = progs%canht_pft(l,n)
          END DO
        END DO

        CALL next_gen_biogeochem(                                              &

          !IN control vars
            asteps_since_triffid,land_pts,nnpft,nmasst,veg3_ctrl,              &
          !IN parms
            litter_parms,red_parms,                                            &
          !INOUT data structures
            veg_state,red_state                                                &
          !OUT diagnostics
          )

        DO n = 1, nnpft
          DO k = 1, land_pts
            l = ainfo%land_index(k)
            progs%canht_pft(l,n)   = veg_state%canht(l,n)
            ainfo%frac_surft(l,n)  = veg_state%frac(l,n)
            trifctltype%npp_acc_pft(l,n) = veg_state%npp_acc(l,n)
          END DO
        END DO

      ELSE ! Triffid based model

        !Vegetation (standalone and UM- differing functionality)

        CALL veg_control(                                                      &
          land_pts, nsurft, dim_cs1,                                           &
          a_step, asteps_since_triffid,                                        &
          phenol_period, triffid_period,                                       &
          l_phenol, l_triffid, l_trif_eq,                                      &
          timestep_real, frac_agr_gb, trif_vars%frac_past_gb,                  &
          psparms%satcon_soilt,                                                &
          trifctltype%g_leaf_acc_pft, trifctltype%g_leaf_phen_acc_pft,         &
          trifctltype%npp_acc_pft,                                             &
          resp_s_acc_gb_um, trifctltype%resp_w_acc_pft,                        &
          cs_pool_gb_um, ainfo%frac_surft, progs%lai_pft,                      &
          psparms%clay_soilt, psparms%z0m_soil_gb,                             &
          progs%canht_pft,                                                     &
          psparms%catch_snow_surft, psparms%catch_surft,                       &
          psparms%infil_surft, psparms%z0_surft, psparms%z0h_bare_surft,       &
          trifctltype%c_veg_pft, trifctltype%cv_gb, trifctltype%lit_c_pft,     &
          trifctltype%lit_c_mn_gb, trifctltype%g_leaf_day_pft,                 &
          trifctltype%g_leaf_phen_pft,trifctltype%lai_phen_pft,                &
          trifctltype%g_leaf_dr_out_pft, trifctltype%npp_dr_out_pft,           &
          trifctltype%resp_w_dr_out_pft,                                       &
          resp_s_dr_out_gb_um, qbase_l_soilt,                                  &
          psparms%sthf_soilt, psparms%sthu_soilt,                              &
          w_flux_soilt, progs%t_soil_soilt,progs%cs_pool_soilt,                &
          !New arguments replacing USE statements
          !trif_vars_mod (IN OUT)
          trif_vars,                                                           &
          !crop_vars_mod (IN)
          crop_vars%rootc_cpft, crop_vars%harvc_cpft, crop_vars%reservec_cpft, &
          crop_vars%stemc_diag_cpft,                                           &
          crop_vars%leafc_diag_cpft, crop_vars%dvi_cpft,                       &
          ! prognostics (IN)
          progs%wood_prod_fast_gb, progs%wood_prod_med_gb,                     &
          progs%wood_prod_slow_gb, progs%frac_agr_prev_gb,                     &
          progs%frac_past_prev_gb, progs%n_inorg_gb, progs%n_inorg_soilt_lyrs, &
          progs%n_inorg_avail_pft, progs%ns_pool_gb,                           &
          progs%triffid_co2_gb, progs%t_soil_soilt_acc,                        &
          ! p_s_parms (IN)
          psparms%bexp_soilt, psparms%sathh_soilt,                             &
          psparms%smvcst_soilt, psparms%smvcwt_soilt,                          &
          psparms%soil_ph_soilt,                                               &
          ! soil_ecosse_vars_mod (OUT)
          soilecosse%n_soil_pool_soilt, dim_soil_n_pool,                       &
          soilecosse%co2_soil_gb, soilecosse%n2o_soil_gb,                      &
          soilecosse%n2o_denitrif_gb, soilecosse%n2o_nitrif_gb,                &
          soilecosse%n2o_partial_nitrif_gb, soilecosse%n2_denitrif_gb,         &
          soilecosse%n_denitrification_gb, soilecosse%n_leach_amm_gb,          &
          soilecosse%n_leach_nit_gb, soilecosse%n_nitrification_gb,            &
          soilecosse%no_soil_gb, soilecosse%soil_c_add, soilecosse%soil_n_add, &
          ! soil_ecosse_vars_mod (IN)
          soilecosse%qbase_l_driver, soilecosse%sthf_driver,                   &
          soilecosse%sthu_driver, soilecosse%tsoil_driver,                     &
          soilecosse%wflux_driver,                                             &
          !ancil_info (IN)
          ainfo%l_lice_point, ainfo%soil_index,                                &
          !TYPES
          soilecosse, urban_param,trifctltype)

      END IF
    END IF

    !Reference evapotranspiration (standalone only)
    IF (l_fao_ref_evapotranspiration) THEN
      trad = ( surftiles_to_gbm(progs%tstar_surft**4, ainfo) )**0.25
      CALL fao_ref_evapotranspiration(soil_pts, ainfo%soil_index,              &
        land_pts, ainfo%land_index, sf_diag%t1p5m,                             &
        sw_down_ij, lw_down_ij, surf_ht_flux_ij, sf_diag%u10m,                 &
        sf_diag%v10m, sf_diag%q1p5m, pstar_ij, trad, trif_vars%fao_et0)
    END IF
    ! End of standalone Jules code

    !-----------------------------------------------------------------------
    !     call to the FLake interface
    !-----------------------------------------------------------------------
    IF (       l_flake_model                                                   &
         .AND. ( .NOT. l_aggregate)                                            &
         .AND. (land_pts > 0 ) ) THEN

      DO k = 1,surft_pts(lake)
        l = ainfo%surft_index(k,lake)
        j=(ainfo%land_index(l) - 1) / row_length + 1
        i = ainfo%land_index(l) - (j-1) * row_length

        !       U* of lake obtained by matching surface stress
        u_s_lake_gb(l) =   u_s_std_surft(l, lake)                              &
                      * SQRT( rhostar(i,j) / rho_water )

        !       Downwelling SW on lake tile
        sw_down_gb(l) = sw_surft(l,lake) / (1.0 - lake_albedo_gb(l))

        !       Take the net SW flux out of the surface heat flux
        !       since this is done separately within FLake.
        surf_ht_flux_lk_gb(l) =   surf_ht_flux_lake_ij(i,j) - sw_surft(l,lake)

        IF ( (nsmax > 0) .AND. (lake_h_snow_gb(l) > EPSILON(1.0)) ) THEN
          ! For the new snow scheme, FLake is forced with zero snow
          ! and the forcing fluxes are those at the bottom of the snowpack.
          ! The order of the following calculations is important.

          ! attenuated DWSW below snow
          dwsw_sub_snow = sw_down_gb(l) *                                      &
                    EXP( -lake_h_snow_gb(l) / h_snow_sw_att )

          ! following last calculation, set snow depth to zero
          lake_h_snow_gb(l) = 0.0

          ! heat flux into the lake becomes
          ! the sub-snow value minus remaining DWSW

          surf_ht_flux_lk_gb(l) =   surf_ht_flux_lake_ij(i,j)                  &
                                  - dhf_surf_minus_soil(l)                     &
                                  - dwsw_sub_snow

          ! now overwrite the DWSW with the remaining sub-snow amount
          sw_down_gb(l) = dwsw_sub_snow

        END IF

      END DO

      trap_frozen   = 0
      trap_unfrozen = 0

        ! DEPENDS ON: flake_interface
      CALL flake_interface( land_pts                                           &
                           ,surft_pts(lake)                                    &
                           ,ainfo%surft_index(:,lake)                          &
                           ,u_s_lake_gb                                        &
                           ,surf_ht_flux_lk_gb                                 &
                           ,sw_down_gb                                         &
                           ,lake_depth_gb                                      &
                           ,lake_fetch_gb                                      &
                           ,coriolis_param_gb                                  &
                           ,timestep_real                                      &
                           ,lake_albedo_gb                                     &
                           ,lake_t_snow_gb                                     &
                           ,lake_t_ice_gb                                      &
                           ,lake_t_mean_gb                                     &
                           ,lake_t_mxl_gb                                      &
                           ,lake_shape_factor_gb                               &
                           ,lake_h_snow_gb                                     &
                           ,lake_h_ice_gb                                      &
                           ,lake_h_mxl_gb                                      &
                           ,lake_t_sfc_gb                                      &
                           ,ts1_lake_gb                                        &
                           ,g_dt_gb                                            &
                           ,trap_frozen                                        &
                           ,trap_unfrozen )

      IF ( trap_frozen > 0 ) THEN
        errcode = -1
        CALL ereport('control', errcode,                                       &
                     'surf_couple_extra-FLake: # zero-divide (frozen) avoided =')
      END IF
      IF ( trap_unfrozen > 0 ) THEN
        errcode = -1
        CALL ereport('control', errcode,                                       &
                     'surf_couple_extra-FLake: # zero-divide (unfrozen) avoided =')
      END IF

    END IF ! Flake

! End of code excluded from LFRic builds

!End of science calls
!===========================================================================
!Expansion to full grid

    DO l = 1, land_pts
      j=(ainfo%land_index(l) - 1) / row_length + 1
      i = ainfo%land_index(l) - (j-1) * row_length
      progs%snow_mass_ij(i,j) = snow_mass_gb(l)
      snowmelt_ij(i,j)  = snowmelt_gb(l)
    END DO

  END IF ! land_pts > 0

  !=============================================================================
  !UM Diagnostics calls

    !For the UM, soil tiling has not been implemented, ie nsoilt = 1, so we can
    !hard-code _soilt variables with index 1 using m = 1

  !=============================================================================
CASE ( cable )
  ! for testing LSM switch
  WRITE(jules_message,'(A)') "CABLE not yet implemented"
  CALL jules_print(RoutineName, jules_message)

  ! initialise all INTENT(OUT) for now until CABLE is implemented
  melt_surft(:,:) = 0.0
  snomlt_surf_htf(:,:) = 0.0
  snowmelt_ij(:,:) = 0.0
  snomlt_sub_htf(:) = 0.0
  sub_surf_roff(:) = 0.0
  surf_roff(:) = 0.0
  tot_tfall(:) = 0.0
  snowmelt_gb(:) = 0.0
  rrun(:) = 0.0
  rflow(:) = 0.0
  snow_soil_htf(:,:) = 0.0

CASE DEFAULT
  errorstatus = 101
  WRITE(jules_message,'(A,I0)') 'Unrecognised surface scheme. lsm_id = ',     &
     lsm_id
  CALL ereport(RoutineName, errorstatus, jules_message)

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE surf_couple_extra
END MODULE surf_couple_extra_mod
