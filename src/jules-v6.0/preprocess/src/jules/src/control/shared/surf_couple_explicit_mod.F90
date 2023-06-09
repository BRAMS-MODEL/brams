! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT****************************************

MODULE surf_couple_explicit_mod

USE jules_gridinit_sf_explicit_mod, ONLY: jules_gridinit_sf_explicit
USE jules_land_sf_explicit_mod,     ONLY: jules_land_sf_explicit
USE jules_ssi_sf_explicit_mod,      ONLY: jules_ssi_sf_explicit
USE jules_griddiag_sf_explicit_mod, ONLY: jules_griddiag_sf_explicit

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE
PUBLIC :: surf_couple_explicit

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SURF_COUPLE_EXPLICIT_MOD'

CONTAINS

!===============================================================================
! Public subroutine
!===============================================================================
SUBROUTINE surf_couple_explicit(                                              &
    !Arguments used by JULES-standalone
    !Misc INTENT(IN) DONE
    bq_1, bt_1, photosynth_act_rad,                                           &
    curr_day_number,                                                          &
    !Forcing INTENT(IN) DONE
    qw, tl, pstar, lw_down,                                                   &
    !Fluxes INTENT(IN) DONE
    sw_surft, tstar,                                                          &
    !Diagnostics, INTENT(INOUT)
    sf_diag,                                                                  &
    !Fluxes INTENT(OUT) DONE
    fqw_1,ftl_1,ftl_surft, fqw_surft, fqw_ice, ftl_ice, fsmc, emis_surft,     &
    !Misc INTENT(OUT)
    radnet_sice, rhokm_1, rhokm_land, rhokm_ssi,                              &
    !Out of explicit and into implicit only INTENT(OUT)
    cdr10m,                                                                   &
    alpha1, alpha1_sea, alpha1_sice, ashtf_prime, ashtf_prime_sea,            &
    ashtf_prime_surft, epot_surft,                                            &
    fraca, resfs, resft, rhokh, rhokh_surft, rhokh_sice, rhokh_sea,           &
    dtstar_surft, dtstar_sea, dtstar_sice,                                    &
    z0hssi, z0h_surft, z0mssi, z0m_surft, chr1p5m, chr1p5m_sice, canhc_surft, &
    wt_ext_surft, flake,                                                      &
    !Out of explicit and into extra only INTENT(OUT)
    hcons_soilt,                                                              &
    !Out of explicit and into implicit and extra INTENT(OUT)
    tile_frac,                                                                &
    !Additional arguments for the UM-----------------------------------------
    !JULES prognostics module
    cs_pool_gb_um,                                                            &
    ti_cat_sicat,                                                             &
    !JULES ancil_info module
    land_pts, nsurft, surft_pts,                                              &
    ! IN input data from the wave model
    charnock_w,                                                               &
    !JULES coastal module
    flandg,                                                                   &
    !JULES aero module
    co2_3d_ij,                                                                &
    !JULES trifctl module
    asteps_since_triffid, resp_s_acc_gb_um, resp_s_gb_um,                     &
    !JULES p_s_parms module
    soil_clay_ij,                                                             &
    !JULES switches module
    l_spec_z0,                                                                &
    !Not in a JULES module
    numcycles, cycleno, z1_uv_top, z1_tq_top, ddmfx,                          &
    l_aero_classic, z0m_scm, z0h_scm, recip_l_mo_sea, rib,                    &
    flandfac, fseafac, fb_surf, u_s, t1_sd, q1_sd, rhostar,                   &
    vshr, resp_s_tot_soilt, emis_soil,                                        &
    !TYPES containing field data (IN OUT)
    crop_vars,psparms,ainfo,trif_vars,aerotype,urban_param,progs,             &
    trifctltype,coast, jules_vars)

!Module imports

!TYPE definitions
USE crop_vars_mod, ONLY: crop_vars_type
USE p_s_parms, ONLY: psparms_type
USE ancil_info,    ONLY: ainfo_type
USE trif_vars_mod, ONLY: trif_vars_type
USE aero,          ONLY: aero_type
USE urban_param_mod, ONLY: urban_param_type
USE prognostics, ONLY: progs_type
USE trifctl, ONLY: trifctl_type
USE coastal, ONLY: coastal_type
USE jules_vars_mod, ONLY: jules_vars_type

!Common modules
USE ereport_mod,              ONLY:                                           &
  ereport
USE jules_sea_seaice_mod,     ONLY:                                           &
  nice, nice_use
USE jules_soil_mod,           ONLY:                                           &
  sm_levels
USE jules_vegetation_mod,     ONLY:                                           &
  l_phenol
USE jules_surface_types_mod,  ONLY:                                           &
  npft, ntype
USE sf_diags_mod, ONLY: strnewsfdiag
USE ozone_vars,               ONLY:                                           &
  o3_gb
USE atm_fields_bounds_mod,    ONLY:                                           &
  pdims_s, pdims, tdims
USE dust_param,               ONLY:                                           &
  ndiv
USE ancil_info,               ONLY:                                           &
  nsoilt

  !New arguments replacing USE statements
USE Fluxes, ONLY:                                                             &
  t_growth_gb, anthrop_heat_surft
USE metstats_mod, ONLY:                                                       &
  metstats_prog
USE bvoc_vars, ONLY:                                                          &
  isoprene_gb, isoprene_pft, terpene_gb , terpene_pft,                        &
  methanol_gb, methanol_pft, acetone_gb, acetone_pft
USE lake_mod, ONLY:                                                           &
  lake_t_ice_gb,lake_t_snow_gb,lake_t_mxl_gb,lake_t_mean_gb,                  &
  lake_h_snow_gb,lake_h_ice_gb,lake_h_mxl_gb,lake_depth_gb,                   &
  ts1_lake_gb,nusselt_gb,g_dt_gb
USE fluxes, ONLY:                                                             &
  sw_sicat, sw_sea


! Modules used by UM only

!Modules that change name between JULES and UM
USE ancil_info,             ONLY:                                             &
  dim_cs1, dim_cs2, co2_dim_len,co2_dim_row
USE switches,               ONLY:                                             &
  l_co2_interactive, l_dust, l_mr_physics
USE aero,                   ONLY:                                             &
  co2_mmr


! reads the lsm switch set in a namelist
USE jules_model_environment_mod, ONLY: lsm_id, jules, cable

! for testing lsm switch
USE jules_print_mgr,             ONLY: jules_message, jules_print

!Dr Hook
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Coupling routine between the UM or JULES system code and land surface
!   explicit science routines. Calls the appropriate LSM-specific code.
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Subroutine arguments

!UM-only arguments
!JULES ancil_info module
INTEGER, INTENT(IN) ::                                                        &
  land_pts,                                                                   &
  nsurft
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  charnock_w(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
INTEGER, INTENT(OUT)  ::                                                      &
  surft_pts(nsurft)

!JULES prognostics module
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  cs_pool_gb_um(land_pts,dim_cs1),                                            &
  ti_cat_sicat(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice)
!JULES coastal module
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  flandg(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)
!JULES aero module
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  co2_3d_ij(co2_dim_len, co2_dim_row)
!JULES trifctl module
INTEGER, INTENT(IN)::                                                         &
  asteps_since_triffid
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                     &
  resp_s_acc_gb_um(land_pts,dim_cs1)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  resp_s_gb_um(land_pts,dim_cs1)
!JULES p_s_parms module
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  soil_clay_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)

!JULES switches module
LOGICAL, INTENT(IN) ::                                                        &
  l_spec_z0
!Not in a JULES module
INTEGER, INTENT(IN) :: numcycles, cycleno
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  !These variables are INTENT(IN) to sf_expl, but not used with the
  !current configuration of standalone JULES (initialised to 0 below)
  z1_uv_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    ! Height of top of lowest uv-layer
  z1_tq_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    ! Height of top of lowest Tq-layer
  ddmfx(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
    ! Convective downdraught mass-flux at cloud base
LOGICAL, INTENT(IN) :: l_aero_classic
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  !These variables are required for prescribed roughness lengths in
  !SCM mode in UM - not used standalone
  z0m_scm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Fixed Sea-surface roughness length for momentum (m).(SCM)
  z0h_scm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
    !Fixed Sea-surface roughness length for heat (m). (SCM)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  recip_l_mo_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
    !Reciprocal of the surface Obukhov  length at sea points. (m-1).
  rib(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
    !Mean bulk Richardson number for lowest layer.
  flandfac(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),      &
  fseafac(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),       &
  fb_surf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Surface flux buoyancy over density (m^2/s^3)
  u_s(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
    !Surface friction velocity (m/s)
  t1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !Standard deviation of turbulent fluctuations of layer 1 temp; used in
    !initiating convection.
  q1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !Standard deviation of turbulent flux of layer 1 humidity; used in
    !initiating convection.
  rhostar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Surface air density
  vshr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    !Magnitude of surface-to-lowest atm level wind shear (m per s).
  resp_s_tot_soilt(dim_cs2,nsoilt),                                           &
    !Total soil respiration (kg C/m2/s).
  emis_soil(land_pts)
    !Emissivity of underlying soil

!Misc INTENT(IN)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  bq_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    !A buoyancy parameter (beta q tilde).
  bt_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    !A buoyancy parameter (beta T tilde).
  photosynth_act_rad(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
    !Net downward shortwave radiation in band 1 (w/m2).
INTEGER, INTENT(IN) ::                                                        &
  curr_day_number

!Forcing INTENT(IN)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                    &
    !Total water content
  tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                    &
    !Ice/liquid water temperature
  pstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
    !Surface pressure (Pascals).
  lw_down(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
    !Surface downward LW radiation (W/m2).

!Fluxes INTENT(IN)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  sw_surft(land_pts,nsurft),                                                  &
    !Surface net SW radiation on land tiles (W/m2).
  tstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
    !GBM surface temperature (K).

!diagnostic array
TYPE (strnewsfdiag), INTENT(IN OUT) :: sf_diag

!Fluxes INTENT(OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  fqw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !Moisture flux between layers (kg per square metre per sec).
  ftl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !FTL(,K) contains net turbulent sensible heat flux into layer K from below;
    !so FTL(,1) is the surface sensible heat, H.(W/m2)
  ftl_surft(land_pts,nsurft),                                                 &
    !Surface FTL for land tiles
  fqw_surft(land_pts,nsurft),                                                 &
    !Surface FQW for land tiles
  fqw_ice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),      &
    !Surface FQW for sea-ice
  ftl_ice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),      &
    !Surface FTL for sea-ice
  fsmc(land_pts,npft),                                                        &
    !Moisture availability factor.
  emis_surft(land_pts,nsurft)
    !Emissivity for land tiles

!Misc INTENT(OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  radnet_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),  &
    !Surface net radiation on sea-ice (W/m2)
  rhokm_1(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),       &
    !Exchange coefficients for momentum on P-grid
  rhokm_land(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),    &
  rhokm_ssi(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)

!Out of explicit and into implicit only INTENT(OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  cdr10m(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),        &
!   Interpolation coefficient for diagnosis of 10 m winds
  alpha1(land_pts,nsurft),                                                    &
    !Mean gradient of saturated specific humidity with respect to temperature
    !between the bottom model layer and tile surfaces
  alpha1_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),  &
    !ALPHA1 for sea-ice.
  alpha1_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
    !ALPHA1 for sea.
  ashtf_prime(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),  &
    !Coefficient to calculate surface heat flux into sea-ice.
  ashtf_prime_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),       &
    !Coefficient to calculate surface heat flux into sea.
  ashtf_prime_surft(land_pts,nsurft),                                         &
    !Coefficient to calculate surface heat flux into land tiles.
  epot_surft(land_pts,nsurft),                                                &
    !Local EPOT for land tiles.
  fraca(land_pts,nsurft),                                                     &
    !Fraction of surface moisture flux with only aerodynamic resistance for
    !snow-free land tiles.
  resfs(land_pts,nsurft),                                                     &
    !Combined soil, stomatal and aerodynamic resistance factor for fraction
    !(1-FRACA) of snow-free land tiles.
  resft(land_pts,nsurft),                                                     &
    !Total resistance factor. FRACA+(1-FRACA)*RESFS for snow-free land, 1 for
    !snow.
  rhokh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !Grid-box surface exchange coefficients
  rhokh_surft(land_pts,nsurft),                                               &
    !Surface exchange coefficients for land tiles
  rhokh_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),   &
    !Surface exchange coefficients for sea-ice
  rhokh_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    !Surface exchange coefficients for sea
  dtstar_surft(land_pts,nsurft),                                              &
    !Change in TSTAR over timestep for land tiles
  dtstar_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
    !Change is TSTAR over timestep for open sea
  dtstar_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),  &
    !Change is TSTAR over timestep for sea-ice
  z0hssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
    !Roughness length for heat and moisture over sea (m).
  z0h_surft(land_pts,nsurft),                                                 &
    !Tile roughness lengths for heat and moisture (m).
  z0mssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
    !Roughness length for momentum over sea (m).
  z0m_surft(land_pts,nsurft),                                                 &
    !Tile roughness lengths for momentum.
  chr1p5m(land_pts,nsurft),                                                   &
    !Ratio of coefffs for calculation of 1.5m temp for land tiles.
  chr1p5m_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
    !CHR1P5M for sea and sea-ice (leads ignored).
  canhc_surft(land_pts,nsurft),                                               &
    !Areal heat capacity of canopy for land tiles (J/K/m2).
  wt_ext_surft(land_pts,sm_levels,nsurft),                                    &
    !Fraction of evapotranspiration which is extracted from each soil layer
    !by each tile.
  flake(land_pts,nsurft)
    !Lake fraction.

!Out of explicit and into extra only INTENT(OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  hcons_soilt(land_pts,nsoilt)
    !Soil thermal conductivity including water and ice

!Out of explicit and into implicit and extra INTENT(OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  tile_frac(land_pts,nsurft)
    !Tile fractions including snow cover in the ice tile.

LOGICAL :: l_dust_diag
  !In standalone, this switch essentially does the same job as l_dust

!TYPES containing field data (IN OUT)
TYPE(crop_vars_type), INTENT(IN OUT) :: crop_vars
TYPE(psparms_type), INTENT(IN OUT) :: psparms
TYPE(ainfo_type), INTENT(IN OUT) :: ainfo
TYPE(trif_vars_type), INTENT(IN OUT) :: trif_vars
TYPE(aero_type), INTENT(IN OUT) :: aerotype
TYPE(urban_param_type), INTENT(IN OUT) :: urban_param
TYPE(progs_type), INTENT(IN OUT) :: progs
TYPE(trifctl_type), INTENT(IN OUT) :: trifctltype
TYPE(coastal_type), INTENT(IN OUT) :: coast
TYPE(jules_vars_type), INTENT(IN OUT) :: jules_vars

INTEGER ::                                                                    &
   i,j,l,n         !Various counters
INTEGER :: errorstatus

!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------

!  Workspace for sea and sea-ice leads
REAL(KIND=real_jlslsm) ::                                                     &
 cd_ssi_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)               &
                             ! Bulk transfer coefficient for
                             ! momentum over sea mean.
,ch_ssi_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)               &
                             ! Bulk transfer coefficient for heat
                             ! and/or moisture over sea mean.
,rhokh_1_sice_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)         &
                             ! Surface exchange coefficient for
                             ! sea-ice.
,rib_sea_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)              &
                             ! Bulk Richardson number
,z0h_sea_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! Roughness length for heat and
                             ! moisture transport

!  Workspace for sea-ice and marginal ice zone
REAL(KIND=real_jlslsm) ::                                                     &
 cd_land_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)              &
                             ! Bulk transfer coefficient for
                             ! momentum over land.
,rib_ice_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,              &
                                                        nice_use)             &
                             ! Bulk Richardson number
,rib_surft(land_pts,nsurft)                                                   &
                             ! RIB for land tiles.
,z0m_ice_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,              &
                                                        nice_use)             &
                             ! Momentum Roughness length.
,z0h_ice_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,              &
                                                        nice_use)             &
                             ! Thermal Roughness length.
,ice_fract_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! Sea ice fraction summed over all cats

!  Workspace for land tiles
REAL(KIND=real_jlslsm) ::                                                     &
 ch_surft_classic(land_pts,nsurft)                                            &
                             ! Bulk transfer coefficient for
                             ! heat for aerosol deposition.
,cd_std_classic(land_pts,nsurft)
                             ! Bulk transfer coefficient for
                             ! momentum for aerosol deposition.
! Workspace for air density calculation
REAL(KIND=real_jlslsm) ::                                                     &
 rhostar_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),           &
 rhostar_ssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)

!  Workspace for variables passed between land and sea and sea-ice routines
LOGICAL ::                                                                    &
 l_cdr10m_snow
                             ! IN Flag indicating if cdr10m
                             !    (an interpolation coefficient) is
                             !    to be calculated for use with
                             !    snow unloading.

! Temp until the model switching is implemented for coupled mode
! Having a parameter until then should hopefully help the compiler eliminate
! dead code

!Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SURF_COUPLE_EXPLICIT'

!-----------------------------------------------------------------------------
!End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

SELECT CASE( lsm_id )
CASE ( jules )

  ! Change 2d UM clay to 1d jules clay content for soil respiration
  ! Soil tiling not currently in the UM, so broadcast ij value to all tiles.
  ! Multi-layer clay not currently in UM so set all layers to same value.
  l_dust_diag = l_dust

  CALL jules_gridinit_sf_explicit (                                           &
    !IN soil/vegetation/land surface data :
    flandg,                                                                   &
    !IN everything not covered so far :
    pstar,tstar,                                                              &
    !IN variables for message passing
    jules_vars%u_1_p_ij, jules_vars%v_1_p_ij,                                 &
    jules_vars%u_0_p_ij, jules_vars%v_0_p_ij,                                 &
    ! INOUT
    sf_diag,                                                                  &
    !OUT Diagnostic not requiring STASH flags :
    fqw_1,ftl_1,                                                              &
    !OUT variables for message passing
    flandfac, fseafac, cdr10m,                                                &
    !OUT data required for mineral dust scheme
    t1_sd,q1_sd,                                                              &
    !OUT data required elsewhere in boundary layer or surface code
    rhostar,vshr,coast%vshr_land_ij,coast%vshr_ssi_ij                         &
    )

  DO i = tdims%i_start,tdims%i_end
    DO j = tdims%j_start,tdims%j_end
      rhostar_land(i,j) = rhostar(i,j)
      rhostar_ssi(i,j) = rhostar(i,j)
    END DO
  END DO


  CALL jules_land_sf_explicit (                                               &
    !IN date-related values
    curr_day_number,                                                          &
    !IN values defining field dimensions and subset to be processed :
    land_pts,                                                                 &
    !IN  parameters for iterative SISL scheme
    numcycles, cycleno,                                                       &
    !IN parameters required from boundary-layer scheme :
    bq_1,bt_1,ainfo%z1_uv_ij,z1_uv_top,ainfo%z1_tq_ij,z1_tq_top,qw,tl,        &
    !IN soil/vegetation/land surface data :
    ainfo%land_index,nsurft,sm_levels,progs%canopy_surft,psparms%catch_surft, &
    psparms%catch_snow_surft, psparms%hcon_soilt,jules_vars%ho2r2_orog_gb,    &
    flandg(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
    progs%snow_surft,jules_vars%sil_orog_land_gb,psparms%smvccl_soilt,        &
    psparms%smvcst_soilt,psparms%smvcwt_soilt,                                &
    psparms%sthf_soilt, psparms%sthu_soilt,psparms%z0_surft,                  &
    psparms%z0h_bare_surft, psparms%z0m_soil_gb,                              &
    !IN input data from the wave model
    charnock_w,                                                               &
    !IN everything not covered so far :
    pstar,lw_down,sw_surft,jules_vars%zh,ddmfx,                               &
    co2_mmr,co2_3d_ij,l_co2_interactive,l_phenol,                             &
    ! asteps_since_triffid,progs%cs_pool_soilt,frac_surft,canht_pft,           &
    asteps_since_triffid,progs%cs_pool_soilt,ainfo%frac_surft,progs%canht_pft,&
    photosynth_act_rad, progs%lai_pft,                                        &
    l_mr_physics,progs%t_soil_soilt,progs%tsurf_elev_surft,                   &
    tstar,progs%tstar_surft,jules_vars%z_land_ij,psparms%albsoil_soilt,       &
    psparms%cosz_ij,                                                          &
    l_aero_classic,l_dust,l_dust_diag,psparms%clay_soilt,o3_gb,               &
    !INOUT diagnostics
    sf_diag,                                                                  &
    !INOUT data :
    progs%gs_gb,trifctltype%g_leaf_acc_pft,trifctltype%npp_acc_pft,           &
    trifctltype%resp_w_acc_pft,                                               &
    trifctltype%resp_s_acc_soilt,rhostar_land,                                &
    !INOUT Diagnostic not requiring STASH flags :
    fqw_1,ftl_1,t1_sd,q1_sd,vshr,coast%vshr_land_ij,                          &
    !OUT Diagnostic not requiring STASH flags :
    ftl_surft,                                                                &
    !OUT variables for message passing
    rhokm_land, cdr10m,                                                       &
    !OUT data required for mineral dust scheme
    aerotype%u_s_std_surft,                                                   &
    !OUT data required elsewhere in boundary layer or surface code
    alpha1,ashtf_prime_surft,fqw_surft,epot_surft,fraca,                      &
    resfs,resft,rhokh_surft,dtstar_surft,z0h_surft,z0m_surft,                 &
    chr1p5m,progs%smc_soilt,hcons_soilt,trifctltype%gpp_gb,trifctltype%npp_gb,&
    trifctltype%resp_p_gb,trifctltype%g_leaf_pft,                             &
    trifctltype%gpp_pft,trifctltype%npp_pft,trifctltype%resp_p_pft,           &
    trifctltype%resp_s_soilt,resp_s_tot_soilt,                                &
    trif_vars%resp_l_pft,trif_vars%resp_r_pft,trifctltype%resp_w_pft,         &
    trif_vars%n_leaf_pft,trif_vars%n_root_pft,trif_vars%n_stem_pft,           &
    trif_vars%lai_bal_pft,                                                    &
    progs%gc_surft,canhc_surft,wt_ext_surft,flake,                            &
    ainfo%surft_index,surft_pts,tile_frac,fsmc,emis_surft,emis_soil,          &
    ! OUT required for classic aerosols
    cd_land_ij,rib_surft,ch_surft_classic,cd_std_classic,                     &
    ! OUT required for sea and sea-ice calculations
    l_cdr10m_snow,                                                            &
    !New arguments replacing USE statements
    !Fluxes (OUT)
    t_growth_gb,                                                              &
    !urban_param (IN)
    urban_param%emisr_gb, urban_param%emisw_gb, urban_param%hwr_gb,           &
    !jules_vars_mod (IN OUT)
    jules_vars%albobs_scaling_surft,                                          &
    !metstats (IN)
    metstats_prog,                                                            &
    !bvoc_vars (OUT)
    isoprene_gb, isoprene_pft, terpene_gb , terpene_pft,                      &
    methanol_gb, methanol_pft, acetone_gb, acetone_pft,                       &
    !trif_vars_mod (OUT)
    trif_vars%fapar_diag_pft, trif_vars%apar_diag_pft, trif_vars%apar_diag_gb,&
    trif_vars%gpp_gb_acc, trif_vars%gpp_pft_acc,                              &
    !crop_vars_mod (IN)
    crop_vars%rootc_cpft, crop_vars%sthu_irr_soilt,                           &
    crop_vars%frac_irr_soilt, crop_vars%frac_irr_surft, crop_vars%dvi_cpft,   &
    !crop_vars_mod (IN OUT)
    crop_vars%resfs_irr_surft,                                                &
    !crop_vars_mod (OUT)
    crop_vars%gs_irr_surft, crop_vars%smc_irr_soilt,                          &
    crop_vars%wt_ext_irr_surft, crop_vars%gc_irr_surft,                       &
    !p_s_parms (IN)
    psparms%bexp_soilt, psparms%sathh_soilt, psparms%v_close_pft,             &
    psparms%v_open_pft,                                                       &
    !urban_param (IN)
    urban_param%wrr_gb,                                                       &
    !Fluxes (IN OUT)
    anthrop_heat_surft,                                                       &
    !prognostics (IN)
    progs%nsnow_surft, progs%sice_surft, progs%sliq_surft,                    &
    progs%snowdepth_surft, progs%tsnow_surft, progs%ds_surft,                 &
    !c_elevate (OUT)
    jules_vars%surf_hgt_surft, jules_vars%lw_down_elevcorr_surft,             &
    !jules_vars_mod (OUT)
    jules_vars%snowdep_surft,                                                 &
    !urban_param (IN)
    urban_param%hgt_gb, urban_param%ztm_gb, urban_param%disp_gb,              &
    !lake_mod (IN)
    lake_t_ice_gb,lake_t_snow_gb,lake_t_mxl_gb,lake_t_mean_gb,                &
    lake_h_snow_gb,lake_h_ice_gb,lake_h_mxl_gb,lake_depth_gb,                 &
    g_dt_gb,                                                                  &
    !lake_mod (OUT)
    nusselt_gb, ts1_lake_gb,                                                  &
    !ancil_info (IN)
    ainfo%l_lice_point, ainfo%l_soil_point,                                   &
    !jules_surface_types (IN)
    jules_vars%diff_frac)


  CALL jules_ssi_sf_explicit (                                                &
    !IN values defining field dimensions and subset to be processed :
    nice, nice_use,                                                           &
    !IN parameters required from boundary-layer scheme :
    bq_1,bt_1,ainfo%z1_uv_ij,z1_uv_top,ainfo%z1_tq_ij,z1_tq_top,qw,tl,        &
    !IN soil/vegetation/land surface data :
    flandg(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
    !IN sea/sea-ice data :
    ainfo%ice_fract_ncat_sicat,progs%k_sice_sicat,                            &
    progs%z0m_sice_fmd, progs%z0m_sice_skin,                                  &
    !IN input data from the wave model
    charnock_w,                                                               &
    !IN everything not covered so far :
    pstar,lw_down,jules_vars%zh,ddmfx,l_mr_physics,progs%ti_sicat,            &
    ainfo%ti_cat_sicat,                                                       &
    coast%tstar_sea_ij,coast%tstar_sice_sicat,jules_vars%z_land_ij,           &
    !IN idealised and SCM things
    l_spec_z0, z0m_scm, z0h_scm,                                              &
    ! IN calcualted in land code
    l_cdr10m_snow,                                                            &
    !INOUT diagnostics
    sf_diag,                                                                  &
    !INOUT data :
    progs%z0msea_ij,rhostar_ssi,fqw_1,ftl_1,t1_sd,q1_sd,cdr10m,               &
    !OUT Diagnostic not requiring STASH flags :
    recip_l_mo_sea,radnet_sice,                                               &
    !OUT variables for message passing
    rhokm_ssi,                                                                &
    !OUT data required elsewhere in boundary layer or surface code
    alpha1_sea,alpha1_sice,ashtf_prime,ashtf_prime_sea,fqw_ice,ftl_ice,       &
    rhokh_sice,rhokh_sea,dtstar_sea,dtstar_sice,chr1p5m_sice,                 &
    coast%vshr_ssi_ij,                                                        &
    ! OUT required for classic aerosols
    cd_ssi_ij,ch_ssi_ij,rhokh_1_sice_ij,rib_sea_ij,z0h_sea_ij,                &
    rib_ice_ij,z0m_ice_ij,z0h_ice_ij,ice_fract_ij,                            &
    !New arguments replacing USE statements
    !ancil_info (IN)
    ainfo%ssi_index, ainfo%sea_index, ainfo%fssi_ij, ainfo%sea_frac,          &
    ainfo%sice_index_ncat, ainfo%sice_frac_ncat,                              &
    !fluxes (IN)
    sw_sicat, sw_sea                                                          &
    )

  CALL jules_griddiag_sf_explicit (                                           &
    !IN values defining field dimensions and subset to be processed :
    land_pts, nice_use,                                                       &
    !IN parameters required from boundary-layer scheme :
    bq_1,bt_1,ainfo%z1_uv_ij,                                                 &
    !IN soil/vegetation/land surface data :
    ainfo%land_index,nsurft,jules_vars%ho2r2_orog_gb,                         &
    flandg(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
    jules_vars%sil_orog_land_gb,                                              &
    !IN everything not covered so far :
    pstar,rhostar_land,rhostar_ssi,jules_vars%zh,tstar,l_aero_classic,        &
    !INOUT data :
    progs%z0msea_ij,cdr10m,                                                   &
    !IN Diagnostic not requiring STASH flags :
    fqw_1,ftl_1,                                                              &
    !IN variables for message passing
    rhokm_land, rhokm_ssi,                                                    &
    !IN data required elsewhere in boundary layer or surface code
    rhokh_surft,rhokh_sice,rhokh_sea,z0h_surft,z0m_surft,                     &
    coast%vshr_land_ij,coast%vshr_ssi_ij,ainfo%surft_index,surft_pts,tile_frac, &
    ! IN required for classic aerosols
    cd_ssi_ij,ch_ssi_ij,rhokh_1_sice_ij,rib_sea_ij,z0h_sea_ij,                &
    cd_land_ij,rib_ice_ij,rib_surft,z0m_ice_ij,z0h_ice_ij,ice_fract_ij,       &
    ch_surft_classic,cd_std_classic,                                          &
    ! INOUT data
    rhostar,                                                                  &
    !INOUT diagnostics
    sf_diag,                                                                  &
    !OUT Diagnostic not requiring STASH flags :
    rhokm_1,rib,                                                              &
    !OUT data required for tracer mixing :
    aerotype%rho_aresist_ij,aerotype%aresist_ij,aerotype%resist_b_ij,         &
    aerotype%rho_aresist_surft,aerotype%aresist_surft,aerotype%resist_b_surft,&
    !OUT data required for mineral dust scheme
    aerotype%r_b_dust_ij,aerotype%cd_std_dust_ij,                             &
    !OUT data required elsewhere in UM system :
    fb_surf,u_s,                                                              &
    !OUT data required elsewhere in boundary layer or surface code
    rhokh,jules_vars%h_blend_orog_ij,z0hssi,z0mssi,jules_vars%z0m_eff_ij,     &
   !ancil_info (IN)
    ainfo%ssi_index, ainfo%sice_index_ncat, ainfo%sice_frac_ncat,             &
   !jules_internal
    jules_vars%unload_backgrnd_pft                                            &
    )



CASE ( cable )
  ! for testing LSM switch
  WRITE(jules_message,'(A)') "CABLE not yet implemented"
  CALL jules_print(RoutineName, jules_message)

  ! initialise all INTENT(OUT) for now until CABLE is implemented
  fqw_1(:,:) = 0.0
  ftl_1(:,:) = 0.0
  ftl_surft(:,:) = 0.0
  fqw_surft(:,:) = 0.0
  fqw_ice(:,:,:) = 0.0
  ftl_ice(:,:,:) = 0.0
  fsmc(:,:) = 0.0
  emis_surft(:,:) = 0.0
  radnet_sice(:,:,:) = 0.0
  rhokm_1(:,:) = 0.0
  rhokm_land(:,:) = 0.0
  rhokm_ssi(:,:) = 0.0
  cdr10m(:,:) = 0.0
  alpha1(:,:) = 0.0
  alpha1_sea(:,:) = 0.0
  alpha1_sice(:,:,:) = 0.0
  ashtf_prime(:,:,:) = 0.0
  ashtf_prime_sea(:,:) = 0.0
  ashtf_prime_surft(:,:) = 0.0
  epot_surft(:,:) = 0.0
  fraca(:,:) = 0.0
  resfs(:,:) = 0.0
  resft(:,:) = 0.0
  rhokh(:,:) = 0.0
  rhokh_surft(:,:) = 0.0
  rhokh_sice(:,:,:) = 0.0
  rhokh_sea(:,:) = 0.0
  dtstar_surft(:,:) = 0.0
  dtstar_sea(:,:) = 0.0
  dtstar_sice(:,:,:) = 0.0
  z0hssi(:,:) = 0.0
  z0h_surft(:,:) = 0.0
  z0mssi(:,:) = 0.0
  z0m_surft(:,:) = 0.0
  chr1p5m(:,:) = 0.0
  chr1p5m_sice(:,:) = 0.0
  canhc_surft(:,:) = 0.0
  wt_ext_surft(:,:,:) = 0.0
  flake(:,:) = 0.0
  hcons_soilt(:,:) = 0.0
  tile_frac(:,:) = 0.0

CASE DEFAULT
  errorstatus = 101
  WRITE(jules_message,'(A,I0)') 'Unrecognised surface scheme. lsm_id = ',     &
     lsm_id
  CALL ereport(RoutineName, errorstatus, jules_message)

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE surf_couple_explicit

END MODULE surf_couple_explicit_mod
