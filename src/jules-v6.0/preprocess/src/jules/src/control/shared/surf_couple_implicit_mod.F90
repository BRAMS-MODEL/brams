! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT****************************************

MODULE surf_couple_implicit_mod

USE jules_grid_update_implicit_mod, ONLY: jules_grid_update_implicit
USE jules_land_sf_implicit_mod,     ONLY: jules_land_sf_implicit
USE jules_ssi_sf_implicit_mod,      ONLY: jules_ssi_sf_implicit
USE jules_griddiag_sf_implicit_mod, ONLY: jules_griddiag_sf_implicit

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE
PUBLIC :: surf_couple_implicit

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SURF_COUPLE_IMPLICIT_MOD'

CONTAINS

!===============================================================================
! Public subroutine
!===============================================================================
SUBROUTINE surf_couple_implicit(                                              &
  !Important switch
  l_correct,                                                                  &
  !Forcing INTENT(IN)
  pstar, lw_down, qw_1, tl_1, u_1, v_1, u_0, v_0,                             &
  !Fluxes INTENT(IN)
  sw_surft, emis_surft,                                                       &
  !Misc INTENT(IN) Many of these simply come out of explicit and into here.
  rhokm_u, rhokm_v, gamma1, gamma2, alpha1, alpha1_sea, alpha1_sice,          &
  ashtf, ashtf_sea, ashtf_surft, du, dv, fraca, resfs, resft, rhokh,          &
  rhokh_surft, rhokh_sice, rhokh_sea,                                         &
  z0hssi, z0mssi, z0h_surft, z0m_surft, chr1p5m,                              &
  chr1p5m_sice, canhc_surft, flake, tile_frac, wt_ext_surft,                  &
  cdr10m_u, cdr10m_v, r_gamma,                                                &
  !Diagnostics, INTENT(INOUT)
  sf_diag,                                                                    &
  !Fluxes INTENT(INOUT)
  fqw_ice, ftl_ice, fqw_surft, fqw_1, ftl_1, ftl_surft,                       &
  !Misc INTENT(INOUT)
  epot_surft, dtstar_surft, dtstar_sea, dtstar_sice, radnet_sice, olr,        &
  !Fluxes INTENT(OUT)
  tstar, le_surft, radnet_surft, e_sea, h_sea, taux_1, tauy_1, ecan_surft, ei,&
  esoil_soilt, ext_soilt, snowmelt, melt_surft,                               &
  ecan, ei_surft, esoil_surft, sea_ice_htf, surf_ht_flux, surf_htf_surft,     &
  !Misc INTENT(OUT)
  error,                                                                      &
  !UM-only arguments
  !JULES ancil_info module
  nsurft, land_pts, surft_pts,                                                &
  !JULES coastal module
  flandg,                                                                     &
  ! Coastal OUT - do this here as needs to be OUT
  tstar_land_ij, tstar_sice_ij,                                               &
  !JULES switches module
  l_co2_interactive, l_mr_physics,                                            &
  co2_3d_ij,                                                                  &
  !Arguments without a JULES module
  ctctq1,dqw1_1,dtl1_1,du_star1,dv_star1,cq_cm_u_1,cq_cm_v_1,flandg_u,        &
  flandg_v, rho1, f3_at_p, uStarGBM,tscrndcl_ssi,tscrndcl_surft,tStbTrans,    &
  ei_sice,rhokh_mix, ti_gb,                                                   &
  !TYPES containing field data (IN OUT)
  crop_vars,ainfo,aerotype,progs,coast, jules_vars)

!Module Imports

!TYPE definitions
USE crop_vars_mod, ONLY: crop_vars_type
USE ancil_info,    ONLY: ainfo_type
USE aero,          ONLY: aero_type
USE prognostics,   ONLY: progs_type
USE coastal,       ONLY: coastal_type
USE jules_vars_mod, ONLY: jules_vars_type

!Common modules
USE ereport_mod,           ONLY: ereport

USE jules_soil_mod,        ONLY: sm_levels
USE jules_sea_seaice_mod,  ONLY: nice_use, nice, l_sice_multilayers
USE sf_diags_mod,          ONLY: strnewsfdiag
USE ancil_info,            ONLY: nsoilt
USE jules_surface_types_mod, ONLY: ntype
!Potential troublemakers
USE atm_fields_bounds_mod, ONLY: tdims,   udims,   vdims,   pdims,            &
                                 tdims_s, udims_s, vdims_s

!New arguments replacing USE statements
USE lake_mod, ONLY:                                                           &
lake_h_ice_gb
USE lake_mod, ONLY:                                                           &
surf_ht_flux_lake_ij
USE fluxes, ONLY:                                                             &
anthrop_heat_surft
USE fluxes, ONLY:                                                             &
surf_ht_store_surft

! Module switches name between UM and JULES-standalone
USE aero, ONLY:                                                               &
  co2_mmr

USE jules_model_environment_mod, ONLY:                                        &
  lsm_id, jules, cable

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

! for testing LSM switch
USE jules_print_mgr,          ONLY: jules_message, jules_print

IMPLICIT NONE


!-----------------------------------------------------------------------------
! Description:
!   Coupling routine between the UM or JULES system code and land surface
!   implicit science routines. Calls the appropriate LSM-specific code.
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

!Arguments

!UM-only arguments
!JULES ancil_info module
INTEGER, INTENT(IN) ::                                                        &
  nsurft, land_pts,                                                           &
  ! ainfo%surft_index(nsurft),                                                        &
  surft_pts(ntype)
!JULES coastal module
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  flandg(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!JULES switches module
LOGICAL, INTENT(IN) ::                                                        &
  l_co2_interactive, l_mr_physics

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  co2_3d_ij(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end)

!Arguments without a JULES module
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  ctctq1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
  dqw1_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
  dtl1_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
  du_star1(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end),      &
  dv_star1(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end),      &
  cq_cm_u_1(udims%i_start:udims%i_end,udims%j_start:udims%j_end),             &
    ! Coefficient in U tri-diagonal implicit matrix
  cq_cm_v_1(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),             &
    ! Coefficient in V tri-diagonal implicit matrix
  flandg_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),              &
    !Land frac (on U-grid, with 1st and last rows undefined or, at present,
    !set to "missing data")
  flandg_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),              &
    !Land frac (on V-grid, with 1st and last rows undefined or, at present,
    !set to "missing data")
  rho1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    !Density on lowest level
  f3_at_p(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Coriolis parameter
  ustargbm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)

REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
    ! BM surface friction velocity
  tscrndcl_ssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
    !Decoupled screen-level temperature over sea or sea-ice
  tscrndcl_surft(land_pts,nsurft),                                            &
    !Decoupled screen-level temperature over land tiles
  tstbtrans(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  tstar_land_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
  tstar_sice_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
  ei_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),      &
    !Sea ice sublimation
  rhokh_mix(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    !Exchange coeffs for moisture.
  ti_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)

! Important switch
LOGICAL, INTENT(IN) ::                                                        &
  l_correct                             ! flag used by the new BL solver

!Forcing INTENT(IN)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  pstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
    !Surface pressure (Pascals).
  lw_down(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Surface downward LW radiation (W/m2).
  qw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    !Total water content
  tl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    !Ice/liquid water temperature
  u_1(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end),           &
    !W'ly wind component (m/s)
  v_1(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end),           &
    !S'ly wind component (m/s)
  u_0(udims%i_start:udims%i_end,udims%j_start:udims%j_end),                   &
    !W'ly component of surface current (m/s).
  v_0(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
    !S'ly component of surface current (m/s).

!Fluxes INTENT(IN)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  sw_surft(land_pts,nsurft),                                                  &
    !Surface net SW radiation on land tiles (W/m2).
  emis_surft(land_pts,nsurft)
    !Emissivity for land tiles

!Misc INTENT(IN) Many of these simply come out of explicit and come
!back into here. Need a module for them.

!Arrays
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  rhokm_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),               &
    !Exchange coefficients for momentum (on U-grid, with 1st and last rows
    !undefined or, at present, set to "missing data")
  rhokm_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),               &
    !Exchange coefficients for momentum (on V-grid, with 1st and last rows
    !undefined or, at present, set to "missing data")
  gamma1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
    !weights for new BL solver
  gamma2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
  alpha1(land_pts,nsurft),                                                    &
    !Mean gradient of saturated specific humidity with respect to temperature
    !between the bottom model layer and tile surfaces
  alpha1_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
    !ALPHA1 for sea.
  alpha1_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),  &
    !ALPHA1 for sea-ice.
  ashtf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),        &
    !Adjusted SEB coefficient for sea-ice
  ashtf_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    !Adjusted SEB coefficient for sea
  ashtf_surft(land_pts,nsurft),                                               &
    !Adjusted SEB coefficient for land tiles.
  du(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end),            &
    !Level 1 increment to u wind field
  dv(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end),            &
    !Level 1 increment to v wind field
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
    !Grid-box surface exchange coefficients (not used for JULES)
  rhokh_surft(land_pts,nsurft),                                               &
    !Surface exchange coefficients for land tiles
  rhokh_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),   &
    !Surface exchange coefficients for sea sea-ice
  rhokh_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    !Surface exchange coefficients for sea
  z0hssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
  z0mssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
    !Roughness lengths over sea (m)
  z0h_surft(land_pts,nsurft),                                                 &
    !Tile roughness lengths for heat and moisture (m).
  z0m_surft(land_pts,nsurft),                                                 &
    !Tile roughness lengths for momentum.
  chr1p5m(land_pts,nsurft),                                                   &
    !Ratio of coefffs for calculation of 1.5m temp for land tiles.
  chr1p5m_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
    !CHR1P5M for sea and sea-ice (leads ignored).
  canhc_surft(land_pts,nsurft),                                               &
    !Areal heat capacity of canopy for land tiles (J/K/m2).
  flake(land_pts,nsurft),                                                     &
    !Lake fraction.
  tile_frac(land_pts,nsurft),                                                 &
    !Tile fractions including snow cover in the ice tile.
  wt_ext_surft(land_pts,sm_levels,nsurft),                                    &
    !Fraction of evapotranspiration extracted from each soil layer by each tile.
  cdr10m_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),              &
    !Ratio of CD's reqd for calculation of 10 m wind. On U-grid; comments as
    !per RHOKM.
  cdr10m_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
    !Ratio of CD's reqd for calculation of 10 m wind. On V-grid; comments as

! Implicit weighting for B.L.
REAL(KIND=real_jlslsm), INTENT(IN)  :: r_gamma

!diagnostic array
TYPE (strnewsfdiag), INTENT(INOUT) :: sf_diag

!Fluxes INTENT(INOUT)
REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
  fqw_ice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),      &
    !Surface FQW for sea-ice
  ftl_ice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),      &
    !Surface FTL for sea-ice
  fqw_surft(land_pts,nsurft),                                                 &
    !Surface FQW for land tiles
  ftl_surft(land_pts,nsurft),                                                 &
    !Surface FTL for land tiles
  fqw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !Moisture flux between layers (kg per square metre per sec) FQW(,1) is
    !total water flux from surface, 'E'.
  ftl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
    !FTL(,K) contains net turbulent sensible heat flux into layer K from below;
    !so FTL(,1) is the surface sensible heat, H.(W/m2)

!Misc INTENT(INOUT)
REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
  epot_surft(land_pts,nsurft),                                                &
    !surface tile potential evaporation
  dtstar_surft(land_pts,nsurft),                                              &
    !Change in TSTAR over timestep for land tiles
  dtstar_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
    !Change is TSTAR over timestep for open sea
  dtstar_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),  &
    !Change is TSTAR over timestep for sea-ice
  radnet_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),  &
    !Surface net radiation on sea-ice (W/m2)
  olr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
    ! IN    TOA - surface upward LW on last radiation timestep
    ! OUT   Corrected TOA outward LW

!Fluxes INTENT(OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  tstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !GBM surface temperature (K).
  le_surft(land_pts,nsurft),                                                  &
    !Surface latent heat flux for land tiles
  radnet_surft(land_pts,nsurft),                                              &
    !Surface net radiation on land tiles (W/m2)
  e_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !Evaporation from sea times leads fraction. Zero over land.
    !(kg per square metre per sec).
  h_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !Surface sensible heat flux over sea times leads fraction (W/m2)
  taux_1(udims%i_start:udims%i_end,udims%j_start:udims%j_end),                &
    !W'ly component of surface wind stress (N/sq m). (On UV-grid with first
    !and last rows undefined or, at present, set to missing data
  tauy_1(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),                &
    !S'ly component of surface wind stress (N/sq m).  On UV-grid; comments as
    !per TAUX
  ecan_surft(land_pts,nsurft),                                                &
    !ECAN for snow-free land tiles
  ei(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                    &
    !Sublimation from lying snow or sea-ice (kg/m2/s).
  ei_surft(land_pts,nsurft),                                                  &
    !EI for land tiles.
  esoil_soilt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nsoilt),    &
    !Surface evapotranspiration from soil moisture store (kg/m2/s).
  ext_soilt(land_pts,nsoilt,sm_levels),                                       &
    !Extraction of water from each soil layer (kg/m2/s).
  snowmelt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
    !Snowmelt (kg/m2/s).
  melt_surft(land_pts,nsurft),                                                &
    !Snowmelt on land tiles (kg/m2/s)
  ecan(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    !Gridbox mean evaporation from canopy/surface store (kg/m2/s). Zero over sea
  esoil_surft(land_pts,nsurft),                                               &
    !ESOIL for snow-free land tiles
  sea_ice_htf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice),      &
    !Heat flux through sea-ice (W/m2, positive downwards). (Not used for JULES)
  surf_ht_flux(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
    !Net downward heat flux at surface over land and sea-ice fraction of gridbox
    !(W/m2).
  surf_htf_surft(land_pts,nsurft)
    !Net downward surface heat flux on tiles (W/m2)

!Misc INTENT(OUT)
INTEGER, INTENT(OUT) ::                                                       &
  error          !0 - AOK; 1 to 7  - bad grid definition detected

!TYPES containing field data (IN OUT)
TYPE(crop_vars_type), INTENT(IN OUT) :: crop_vars
TYPE(ainfo_type), INTENT(IN OUT) :: ainfo
TYPE(aero_type), INTENT(IN OUT) :: aerotype
TYPE(progs_type), INTENT(IN OUT) :: progs
TYPE(coastal_type), INTENT(IN OUT) :: coast
TYPE(jules_vars_type), INTENT(IN OUT) :: jules_vars

!-----------------------------------------------------------------------------
! Workspace
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                     &
 ice_fract_cat_use(tdims%i_start:tdims%i_end,                                 &
                   tdims%j_start:tdims%j_end,nice_use)                        &
                             ! Sea ice category fractions
                             ! If nice_use=1, this is the total ice
                             ! fraction
,surf_ht_flux_sice_sm_ij(tdims%i_start:tdims%i_end,                           &
                      tdims%j_start:tdims%j_end)                              &
                             ! Sea area mean seaice surface heat flux
,ei_sice_sm_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)           &
                             ! Sea area mean sea ice sublimation
,ei_land_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)              &
                             ! Sublimation from lying snow
                             ! (kg/m2/s).
,tstar_ssi_old_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)        &
                             ! Sea and sea-ice surface temperature
                             ! at beginning of timestep --
                             ! Required only for decoupled diagnosis,
                             ! so allocatable, and local since it is
                             ! used only on the predictor step
,tstar_surft_old(land_pts,nsurft)
                             ! Tile surface temperatures at
                             ! beginning of timestep.


!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
INTEGER :: i,j,n

!Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SURF_COUPLE_IMPLICIT'

!-----------------------------------------------------------------------------
!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

SELECT CASE( lsm_id )
CASE ( jules )

  CALL  jules_grid_update_implicit(                                           &
  ! IN values defining field dimensions and subset to be processed :
          land_pts,ainfo%land_index,nice,nice_use,nsurft,ainfo%surft_index,   &
          surft_pts, tile_frac,flandg,                                        &
  ! IN sea/sea-ice data :
          ainfo%ice_fract_ij,ainfo%ice_fract_ncat_sicat,                      &
  ! IN everything not covered so far :
          rhokm_u,rhokm_v,r_gamma,                                            &
          gamma1,gamma2,alpha1,alpha1_sea,alpha1_sice,                        &
          ashtf,ashtf_sea,ashtf_surft,                                        &
          du,dv,resft,rhokh_surft,rhokh_sice,rhokh_sea,ctctq1,                &
          dqw1_1,dtl1_1,du_star1,dv_star1,cq_cm_u_1,cq_cm_v_1,                &
          l_correct,flandg_u,flandg_v,progs%snow_surft,                       &
  ! INOUT data :
          epot_surft,fqw_ice,ftl_ice,dtstar_surft,dtstar_sea,dtstar_sice,     &
          fqw_surft, fqw_1,ftl_1,ftl_surft,                                   &
          coast%taux_land_ij,coast%taux_ssi_ij,coast%tauy_land_ij,            &
          coast%tauy_ssi_ij, coast%taux_land_star,coast%tauy_land_star,       &
          coast%taux_ssi_star, coast%tauy_ssi_star,                           &
  ! OUT data required elsewhere in UM system :
          e_sea,h_sea,taux_1,tauy_1,ice_fract_cat_use                         &
        )

  IF ( .NOT. l_correct ) THEN
    CALL jules_land_sf_implicit (                                             &
    ! IN values defining field dimensions and subset to be processed :
            land_pts,ainfo%land_index,nsurft,ainfo%surft_index,surft_pts,     &
            sm_levels,                                                        &
            canhc_surft,progs%canopy_surft,flake,progs%smc_soilt,tile_frac,   &
            wt_ext_surft,coast%fland,flandg,                                  &
    ! IN everything not covered so far :
            lw_down,sw_surft,progs%t_soil_soilt,r_gamma,alpha1,ashtf_surft,   &
            jules_vars%dtrdz_charney_grid_1_ij,fraca,resfs,resft,rhokh_surft, &
            emis_surft,progs%snow_surft,dtstar_surft,                         &
    ! INOUT data :
            progs%tstar_surft,fqw_surft,fqw_1,ftl_1,ftl_surft,sf_diag,        &
    ! OUT Diagnostic not requiring STASH flags :
            ecan,ei_surft,esoil_surft,coast%surf_ht_flux_land_ij,             &
            ei_land_ij,surf_htf_surft,                                        &
    ! OUT data required elsewhere in UM system :
            tstar_land_ij,le_surft,radnet_surft,ecan_surft,esoil_soilt,       &
            ext_soilt,snowmelt,melt_surft,tstar_surft_old,error,              &
            !New arguments replacing USE statements
            ! lake_mod (IN)
            lake_h_ice_gb,                                                    &
            ! lake_mod (OUT)
            surf_ht_flux_lake_ij,                                             &
            ! fluxes (IN)
            anthrop_heat_surft,                                               &
            ! fluxes (OUT)
            surf_ht_store_surft,                                              &
            ! c_elevate (IN)
            jules_vars%lw_down_elevcorr_surft,                                &
            ! prognostics (IN)
            progs%nsnow_surft,                                                &
            ! jules_vars_mod (IN)
            jules_vars%snowdep_surft,                                         &
            !TYPES containing field data (IN OUT)
            crop_vars)
  END IF ! IF .NOT. L_correct

  CALL jules_ssi_sf_implicit (                                                &
  ! IN values defining field dimensions and subset to be processed :
          nice,nice_use,flandg,                                               &
  ! IN sea/sea-ice data :
          ainfo%ice_fract_ij,ainfo%ice_fract_ncat_sicat,ice_fract_cat_use,    &
          progs%k_sice_sicat,progs%di_ncat_sicat,ainfo%sstfrz_ij,             &
  ! IN everything not covered so far :
          lw_down,r_gamma,alpha1_sice,ashtf,                                  &
          jules_vars%dtrdz_charney_grid_1_ij,rhokh_sice,l_correct,            &
  ! INOUT data :
          fqw_ice,ftl_ice,coast%tstar_sice_sicat,coast%tstar_ssi_ij,          &
          coast%tstar_sea_ij,                                                 &
          radnet_sice,fqw_1,ftl_1,progs%ti_sicat,sf_diag,                     &
  ! OUT Diagnostic not requiring STASH flags :
          ti_gb,sea_ice_htf,coast%surf_ht_flux_sice_sicat,                    &
  ! OUT data required elsewhere in UM system :
          tstar_sice_ij,e_sea,h_sea,ei_sice,dtstar_sea,dtstar_sice,           &
          surf_ht_flux_sice_sm_ij,ei_sice_sm_ij,tstar_ssi_old_ij,             &
  ! ancil_info (IN)
          ainfo%ssi_index, ainfo%sice_index, ainfo%sice_index_ncat,           &
          ainfo%fssi_ij, ainfo%sice_frac, ainfo%sice_frac_ncat                &
        )

  CALL jules_griddiag_sf_implicit (                                           &
  ! IN values defining field dimensions and subset to be processed :
          land_pts,ainfo%land_index,nice_use,nsurft,ainfo%surft_index,        &
          surft_pts, tile_frac,flandg,l_mr_physics,                           &
  ! IN sea/sea-ice data :
          ainfo%ice_fract_ij,ice_fract_cat_use,u_0,v_0,ei_sice,               &
  ! IN everything not covered so far :
          pstar,qw_1,tl_1,u_1,v_1,du,dv,                                      &
          resft,rhokh,ainfo%z1_tq_ij,                                         &
          z0hssi,z0mssi,z0h_surft,z0m_surft,                                  &
          cdr10m_u,cdr10m_v,                                                  &
          chr1p5m,chr1p5m_sice,ctctq1,                                        &
          dqw1_1,dtl1_1,du_star1,dv_star1,cq_cm_u_1,cq_cm_v_1,                &
          l_correct,emis_surft,                                               &
          coast%tstar_sice_sicat,coast%tstar_ssi_ij,progs%tstar_surft,        &
          coast%tstar_sea_ij,fqw_1,                                           &
          coast%surf_ht_flux_land_ij,surf_ht_flux_sice_sm_ij,ei_land_ij,      &
          ei_sice_sm_ij,tstar_land_ij,tstar_ssi_old_ij,tstar_surft_old,       &
          taux_1,tauy_1,                                                      &
  ! IN variables used to calculate cooling at the screen level
          l_co2_interactive, co2_mmr, co2_3d_ij,rho1, f3_at_p,                &
          ustargbm,                                                           &
  ! INOUT data :
          ftl_1,olr,                                                          &
          tscrndcl_ssi,tscrndcl_surft,tStbTrans,sf_diag,                      &
  ! OUT Diagnostic not requiring STASH flags :
          surf_ht_flux,                                                       &
  ! OUT data required elsewhere in UM system :
          tstar,ei,rhokh_mix                                                  &
        )


CASE ( cable )
  ! for testing LSM switch
  WRITE(jules_message,'(A)') "CABLE not yet implemented"
  CALL jules_print(RoutineName, jules_message)

  ! initialise all INTENT(OUT) for now until CABLE is implemented
  tstar(:,:) = 0.0
  le_surft(:,:) = 0.0
  radnet_surft(:,:) = 0.0
  e_sea(:,:) = 0.0
  h_sea(:,:) = 0.0
  taux_1(:,:) = 0.0
  tauy_1(:,:) = 0.0
  ecan_surft(:,:) = 0.0
  ei(:,:) = 0.0
  esoil_soilt(:,:,:) = 0.0
  ext_soilt(:,:,:) = 0.0
  snowmelt(:,:) = 0.0
  melt_surft(:,:) = 0.0
  ecan(:,:) = 0.0
  ei_surft(:,:) = 0.0
  esoil_surft(:,:) = 0.0
  sea_ice_htf(:,:,:) = 0.0
  surf_ht_flux(:,:) = 0.0
  surf_htf_surft(:,:) = 0.0
  error = 0

CASE DEFAULT
  error = 101
  WRITE(jules_message,'(A,I0)') 'Unrecognised surface scheme. lsm_id = ',     &
     lsm_id
  CALL ereport(RoutineName, error, jules_message)

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE surf_couple_implicit

END MODULE surf_couple_implicit_mod
