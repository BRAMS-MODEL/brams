#if !defined(UM_JULES)
SUBROUTINE control (                                                          &
!   Scalar arguments (INTENT IN)
    timestep,                                                                 &
!   Atmospheric forcing (INTENT IN)
    qw_1_ij, tl_1_ij, u_0_ij, v_0_ij, u_1_ij, v_1_ij, pstar_ij,               &
    ls_rain_ij, con_rain_ij, ls_snow_ij, con_snow_ij, sw_down_ij, lw_down_ij, &
!   Gridbox mean surface fluxes (INTENT OUT)
    fqw_1_ij, ftl_1_ij, taux_1_ij, tauy_1_ij,                                 &
    !TYPES containing field data (IN OUT)
    crop_vars,psparms,toppdm,fire_vars,ainfo,trif_vars,soilecosse, aerotype,  &
    urban_param,progs,trifctltype,coast,jules_vars)


!-------------------------------------------------------------------------------
! Control level routine, to call the main parts of the model.
!-------------------------------------------------------------------------------

!TYPE definitions
USE crop_vars_mod, ONLY: crop_vars_type
USE p_s_parms, ONLY: psparms_type
USE top_pdm, ONLY: top_pdm_type
USE fire_vars_mod, ONLY: fire_vars_type
USE ancil_info,    ONLY: ainfo_type
USE trif_vars_mod, ONLY: trif_vars_type
USE soil_ecosse_vars_mod,     ONLY: soil_ecosse_vars_type
USE aero, ONLY: aero_type
USE urban_param_mod, ONLY: urban_param_type
USE prognostics, ONLY: progs_type
USE coastal, ONLY: coastal_type
USE top_pdm, ONLY: top_pdm_type
USE trifctl, ONLY: trifctl_type
USE jules_vars_mod, ONLY: jules_vars_type

!Import the main subroutines we will use
USE surf_couple_radiation_mod, ONLY: surf_couple_radiation
USE surf_couple_explicit_mod,  ONLY: surf_couple_explicit
USE surf_couple_implicit_mod,  ONLY: surf_couple_implicit
USE surf_couple_extra_mod,     ONLY: surf_couple_extra
USE update_mod,                ONLY: calc_downward_rad

!Variables- modules in alphabetical order

USE ancil_info,               ONLY:                                           &
  land_pts, nsurft, surft_pts, row_length, rows, dim_cs1, dim_cs2,            &
  lice_pts, soil_pts, nsoilt

USE atm_fields_bounds_mod,    ONLY:                                           &
  tdims, pdims, udims, vdims, tdims_s, pdims_s, udims_s, vdims_s

USE bl_option_mod,            ONLY: on

USE coastal,                  ONLY: flandg

USE datetime_mod,             ONLY: l_360, l_leap

USE datetime_utils_mod,       ONLY: day_of_year

USE diag_swchs,               ONLY: stf_sub_surf_roff

USE fluxes,                   ONLY:                                           &
  alb_surft, e_sea_ij, ecan_ij, ecan_surft, ei_ij, ei_surft, esoil_ij_soilt,  &
  esoil_surft, ext_soilt, fqw_surft, fqw_sicat, fsmc_pft,                     &
  ftl_sicat, ftl_surft, h_sea_ij, hf_snow_melt_gb, land_albedo_ij,            &
  le_surft, melt_surft, sea_ice_htf_sicat, snomlt_sub_htf_gb, snow_melt_gb,   &
  snowmelt_ij, sub_surf_roff_gb, surf_ht_flux_ij, surf_htf_surft,             &
  surf_roff_gb, radnet_surft, tot_tfall_gb, tstar_ij,                         &
  sw_surft, emis_surft, rflow_gb, rrun_gb, snow_soil_htf, z0m_surft, z0h_surft

USE jules_radiation_mod,      ONLY: i_sea_alb_method, l_cosz

USE jules_sea_seaice_mod,     ONLY: l_ctile, buddy_sea, nice_use

USE jules_snow_mod,           ONLY: nsmax

USE jules_soil_mod,           ONLY: sm_levels

USE jules_surface_mod,        ONLY: l_aggregate, l_flake_model

USE jules_surface_types_mod,  ONLY: ntype, lake

USE lake_mod, ONLY: lake_albedo_gb                                            &
                   , lake_t_snow_gb                                           &
                   , lake_h_snow_gb                                           &
                   , lake_h_ice_gb                                            &
                   , lake_t_sfc_gb

USE model_time_mod,           ONLY: current_time

USE planet_constants_mod,     ONLY: c_virtual

USE sf_diags_mod,             ONLY: sf_diag

USE switches,                 ONLY:                                           &
  l_mr_physics, l_co2_interactive, l_spec_z0

USE theta_field_sizes,        ONLY: t_i_length, t_j_length

USE trifctl,                  ONLY: asteps_since_triffid

USE water_constants_mod,      ONLY: rho_water

USE zenith_mod,               ONLY: zenith

USE jules_rivers_mod,         ONLY:                                           &
   tot_surf_runoff_gb, tot_sub_runoff_gb, acc_lake_evap_gb


!-------------------------------------------------------------------------------

IMPLICIT NONE
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-------------------------------------------------------------------------------
INTEGER, INTENT(IN) :: timestep    ! IN Atmospheric timestep number.

!-------------------------------------------------------------------------------
! Arguments with intent(in)
!-------------------------------------------------------------------------------
!Forcing INTENT(IN)
REAL, INTENT(IN) ::                                                           &
  qw_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
!    Total water content (Kg/Kg)
  tl_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
!    Ice/liquid water temperature (k)
  u_0_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
!    W'ly component of surface current (m/s)
  v_0_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
!    S'ly component of surface current (m/s)
  u_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
!    W'ly wind component (m/s)
  v_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
!    S'ly wind component (m/s)
  pstar_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
!    Surface pressure (Pascals)
  ls_rain_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
!    Large-scale rain (kg/m2/s)
  con_rain_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),           &
!     Convective rain (kg/m2/s)
  ls_snow_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
!    Large-scale snowfall (kg/m2/s)
  con_snow_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),           &
!    Convective snowfall (kg/m2/s)
  sw_down_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
!    Surface downward SW radiation (W/m2)
  lw_down_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!    Surface downward LW radiation (W/m2)

!-------------------------------------------------------------------------------
! Arguments with intent(out)
!-------------------------------------------------------------------------------
TYPE(crop_vars_type), INTENT(IN OUT) :: crop_vars
TYPE(psparms_type), INTENT(IN OUT) :: psparms
TYPE(top_pdm_type), INTENT(IN OUT) :: toppdm
TYPE(fire_vars_type), INTENT(IN OUT) :: fire_vars
TYPE(ainfo_type), INTENT(IN OUT) :: ainfo
TYPE(trif_vars_type), INTENT(IN OUT) :: trif_vars
TYPE(soil_ecosse_vars_type), INTENT(IN OUT) :: soilecosse
TYPE(aero_type), INTENT(IN OUT) :: aerotype
TYPE(urban_param_type), INTENT(IN OUT) :: urban_param
TYPE(progs_type), INTENT(IN OUT) :: progs
TYPE(trifctl_type), INTENT(IN OUT) :: trifctltype
TYPE(coastal_type), INTENT(IN OUT) :: coast
TYPE(jules_vars_type), INTENT(IN OUT) :: jules_vars

!-------------------------------------------------------------------------------
! Arguments with intent(out)
!-------------------------------------------------------------------------------
!Fluxes INTENT(OUT)
REAL, INTENT(OUT) ::                                                          &
  fqw_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
!    Moisture flux between layers (kg per square metre per sec).
  ftl_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
!    FTL(,K) contains net turbulent sensible heat flux into layer K from below;
!    so FTL(,1) is the surface sensible heat, H.(W/m2)
  taux_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
!    W'ly component of surface wind stress (N/sq m)
  tauy_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!    S'ly component of surface wind stress (N/sq m)
!    On V-grid; comments as per TAUX

!-------------------------------------------------------------------------------
!Parameters
!-------------------------------------------------------------------------------

LOGICAL, PARAMETER ::                                                         &
  l_aero_classic = .FALSE.
    !switch for CLASSIC aerosol- NEVER USED!

INTEGER, PARAMETER ::                                                         &
  n_swbands = 1,                                                              &
    !  Nnumber of SW bands
  numcycles = 1,                                                              &
    !Number of cycles (iterations) for iterative SISL.
  cycleno = 1,                                                                &
    !Iteration no
  n_proc = 2,                                                                 &
  river_row_length = 1,                                                       &
  river_rows = 1,                                                             &
  aocpl_row_length = 1,                                                       &
  aocpl_p_rows = 1

REAL, PARAMETER ::                                                            &
  r_gamma = 2.0,                                                              &
    ! Implicit weighting coefficient
  spec_band_bb(n_swbands) = -1.0
    !spectral band boundary (bb=-1)

!TYPES containing field data (IN OUT)

!-------------------------------------------------------------------------------
! Local variables
! In the vast majority of cases, these are arguments only required when
! coupled to the UM.
!-------------------------------------------------------------------------------

LOGICAL ::                                                                    &
  l_correct,                                                                  &
  land_sea_mask(row_length,rows)

INTEGER ::                                                                    &
  a_steps_since_riv = 0,                                                      &
  g_p_field,                                                                  &
  g_r_field,                                                                  &
  global_row_length,                                                          &
  global_rows,                                                                &
  global_river_row_length,                                                    &
  global_river_rows,                                                          &
  error,                                                                      &
    ! OUT 0 - AOK; 1 to 7  - bad grid definition detected
  curr_day_number

REAL ::                                                                       &
  olr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
    !TOA-surface upward LW on last radiation timestep Corrected TOA outward LW
  sea_ice_albedo(row_length,rows,4),                                          &
    ! Sea ice albedo
  ws10m(t_i_length,t_j_length),                                               &
    ! 10m wind speed (m s-1)
  photosynth_act_rad(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),    &
    ! Net downward shortwave radiation in band 1 (w/m2).
  bt_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    ! A buoyancy parameter (beta T tilde).
  bq_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    ! A buoyancy parameter (beta q tilde).
  radnet_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),  &
    ! Surface net radiation on sea-ice (W/m2)
  rhokm_land(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),    &
  rhokm_ssi(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),     &
  rhostar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    ! surface air density from ideal gas equation, used for FLake call
  cdr10m(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),        &
  alpha1(land_pts,nsurft),                                                    &
    !Mean gradient of saturated specific humidity with respect to temperature
    !between the bottom model layer and tile surfaces
  alpha1_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
    ! ALPHA1 for open sea.
  alpha1_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),  &
    ! ALPHA1 for sea-ice.
  ashtf_prime(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),  &
    ! Adjusted SEB coefficient for sea-ice
  ashtf_prime_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),       &
    ! Adjusted SEB coefficient for open sea
  ashtf_prime_surft(land_pts,nsurft),                                         &
    ! Adjusted SEB coefficient for land tiles
  rhokm_1(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),       &
    ! Exchange coefficients for momentum on P-grid
  epot_surft(land_pts,nsurft),                                                &
    ! Local EPOT for land tiles.
  fraca(land_pts,nsurft),                                                     &
    !Fraction of surface moisture flux with only aerodynamic resistance for
    !snow-free land tiles.
  resfs(land_pts,nsurft),                                                     &
    ! Combined soil, stomatal and aerodynamic resistance factor for fraction
    !(1-FRACA) of snow-free land tiles.
  resft(land_pts,nsurft),                                                     &
    !Total resistance factor. FRACA+(1-FRACA)*RESFS for snow-free land, 1
    !for snow.
  rhokh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    ! Grid-box surface exchange coefficients
  rhokh_surft(land_pts,nsurft),                                               &
    ! Surface exchange coefficients for land tiles
  rhokh_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),   &
    ! Surface exchange coefficients for sea-ice
  rhokh_sea_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
    ! Surface exchange coefficients for sea
  dtstar_ij_surft(land_pts,nsurft),                                           &
    ! Change in tstar_ij over timestep for land tiles
  dtstar_ij_sea(tdims%i_start:tdims%i_end,                                    &
                tdims%j_start:tdims%j_end),                                   &
    ! Change is tstar_ij over timestep for open sea
  dtstar_ij_sice(tdims%i_start:tdims%i_end,                                   &
                tdims%j_start:tdims%j_end,nice_use),                          &
    ! Change is tstar_ij over timestep for sea-ice
  z0hssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
    ! Roughness length for heat and moisture over sea (m).
  z0mssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
    ! Roughness length for momentum over sea (m).
  chr1p5m(land_pts,nsurft),                                                   &
    ! Ratio of coefffs for calculation of 1.5m temp for land tiles.
  chr1p5m_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
    ! CHR1P5M for sea and sea-ice (leads ignored).
  canhc_surft(land_pts,nsurft),                                               &
    ! Areal heat capacity of canopy for land tiles (J/K/m2).
  wt_ext_surft(land_pts,sm_levels,nsurft),                                    &
    !Fraction of evapotranspiration which is extracted from each soil layer
    !by each tile.
  flake(land_pts,nsurft),                                                     &
    !Lake fraction.
  tile_frac(land_pts,nsurft),                                                 &
    !Tile fractions including snow cover in the ice tile.
  hcons_soilt(land_pts,nsoilt),                                               &
    !Thermal conductivity of top soil layer, including water and ice (W/m/K)
  flandfac_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),            &
  flandfac_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),            &
  fseafac_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),             &
  fseafac_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),             &
  du(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end),            &
    !Level 1 increment to u wind field
  dv(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end),            &
    !Level 1 increment to v wind field
  r_gamma1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
    !weights for new BL solver
  r_gamma2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
  chloro(row_length,rows),                                                    &
    !nr surface chlorophyll content

!Radiation
  albobs_sc_ij(t_i_length,t_j_length,nsurft,2),                               &
    !albedo scaling factors to obs
  open_sea_albedo(row_length,rows,2,n_swbands),                               &
    !Surface albedo for Open Sea (direct and diffuse components, for each
    !band, with zeros for safety where no value applies)

!Implicit
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
  uStargbm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
    ! BM surface friction velocity
  tscrndcl_ssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
    !Decoupled screen-level temperature over sea or sea-ice
  tscrndcl_surft(land_pts,nsurft),                                            &
    !Decoupled screen-level temperature over land tiles
  tstbtrans(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    !Time since the transition
  ei_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),      &
    !Sea ice sublimation
  rhokh_mix(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    !Exchange coeffs for moisture.
  ti_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &

!Explicit
  !These variables are INTENT(IN) to sf_expl, but not used with the
  !current configuration of standalone JULES (initialised to 0 below)
  z1_uv_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    ! Height of top of lowest uv-layer
  z1_tq_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    ! Height of top of lowest Tq-layer
  ddmfx(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    ! Convective downdraught mass-flux at cloud base
  cs_pool_gb_um(land_pts,dim_cs1),                                            &
  resp_s_acc_gb_um(land_pts,dim_cs1),                                         &
  resp_s_gb_um(land_pts,dim_cs1),                                             &
  soil_clay_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
  ! Charnock parameter from the wave model
  charnock_w(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &

  !These variables are required for prescribed roughness lengths in
  !SCM mode in UM - not used standalone
  z0m_scm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Fixed Sea-surface roughness length for momentum (m).(SCM)
  z0h_scm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Fixed Sea-surface roughness length for heat (m). (SCM)
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
  vshr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    !Magnitude of surface-to-lowest atm level wind shear (m per s).
  resp_s_tot_soilt(dim_cs2,nsoilt),                                           &
    !Total soil respiration (kg C/m2/s).
  emis_soil(land_pts),                                                        &
  cca_2d(row_length,rows),                                                    &
  cs_ch4_soilt(land_pts,nsoilt),                                              &
  dhf_surf_minus_soil(land_pts),                                              &
  flash_rate_ancil(row_length,rows),                                          &
  pop_den_ancil(row_length,rows),                                             &
  inlandout_atm_gb(land_pts),                                                 &
  delta_lambda,                                                               &
  delta_phi,                                                                  &
  xx_cos_theta_latitude(tdims_s%i_start:tdims_s%i_end,                        &
                        tdims_s%j_start:tdims_s%j_end),                       &
  xpa(aocpl_row_length+1),                                                    &
  xua(0:aocpl_row_length),                                                    &
  xva(aocpl_row_length+1),                                                    &
  ypa(aocpl_p_rows),                                                          &
  yua(aocpl_p_rows),                                                          &
  yva(0:aocpl_p_rows),                                                        &
  trivdir(river_row_length, river_rows),                                      &
  trivseq(river_row_length, river_rows),                                      &
  r_area(row_length, rows),                                                   &
  slope(row_length, rows),                                                    &
  flowobs1(row_length, rows),                                                 &
  r_inext(row_length, rows),                                                  &
  r_jnext(row_length, rows),                                                  &
  r_land(row_length, rows),                                                   &
  substore(row_length, rows),                                                 &
  surfstore(row_length, rows),                                                &
  flowin(row_length, rows),                                                   &
  bflowin(row_length, rows),                                                  &
  twatstor(river_row_length, river_rows),                                     &
  ls_graup_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),           &
  ls_rainfrac_land(land_pts)

INTEGER :: l      !loop counter
INTEGER :: point  !loop counter
INTEGER :: i, j   !point indices

!------------------------------------------------------------------------------
!End of header

!Initialise the olr diagnostic
olr(:,:) = 0.0

!-----------------------------------------------------------------------
!   Conditional allocation for calculation of diagnostics
!-----------------------------------------------------------------------

IF (sf_diag%suv10m_n) THEN
  ALLOCATE(sf_diag%cdr10m_n(pdims_s%i_start:pdims_s%i_end,                    &
                            pdims_s%j_start:pdims_s%j_end))
  ALLOCATE(sf_diag%cdr10m_n_u(udims%i_start:udims%i_end,                      &
                              udims%j_start:udims%j_end))
  ALLOCATE(sf_diag%cdr10m_n_v(vdims%i_start:vdims%i_end,                      &
                              vdims%j_start:vdims%j_end))
  ALLOCATE(sf_diag%cd10m_n(pdims_s%i_start:pdims_s%i_end,                     &
                           pdims_s%j_start:pdims_s%j_end))
ELSE
  ALLOCATE(sf_diag%cdr10m_n(1,1))
  ALLOCATE(sf_diag%cdr10m_n_u(1,1))
  ALLOCATE(sf_diag%cdr10m_n_v(1,1))
  ALLOCATE(sf_diag%cd10m_n(1,1))
END IF

!Calculate the cosine of the zenith angle
IF ( l_cosz ) THEN
  CALL zenith(psparms%cosz_ij)
ELSE
  !     Set cosz to a default of 1.0
  psparms%cosz_ij(:,:) = 1.0
END IF

IF ( i_sea_alb_method == 3 ) THEN
  ! Calculate the 10m wind speed.
  IF (timestep == 1) THEN
    ! the u/v at 10m not calculated yet, so make a first guess:
    ws10m(:,:) = 2.0
  ELSE
    ws10m(:,:) = SQRT(sf_diag%u10m(:,:)**2 + sf_diag%v10m(:,:)**2)
  END IF
END IF

CALL surf_couple_radiation(                                                   &
  !Fluxes INTENT(IN)
  tstar_ij,                                                                   &
  !Misc INTENT(IN)
  ws10m, chloro,                                                              &
  n_swbands, n_swbands, spec_band_bb, spec_band_bb,                           &
  !Misc INTENT(OUT)
  sea_ice_albedo,                                                             &
  !Fluxes INTENT(OUT)
  alb_surft, land_albedo_ij,                                                  &
  !(ancil_info mod)
  nsurft, land_pts, surft_pts, row_length, rows,                              &
  !UM-only args: INTENT(IN)
  !(coastal mod)
  flandg,                                                                     &
  !(prognostics mod)
  !Warning snow_surft passed as array due to UM problems
  progs%snow_surft,                                                           &
  !UM-only args: INTENT(OUT)
  albobs_sc_ij, open_sea_albedo,                                              &
  !TYPES containing field data (IN OUT)
  psparms, ainfo, urban_param, progs, coast, jules_vars                       &
  )

!All the calculations for the radiation downward components have been put into
!a subroutine. For standalone only.
CALL calc_downward_rad(sea_ice_albedo, photosynth_act_rad, ainfo, progs,      &
                       coast, jules_vars)

!Calculate buoyancy parameters bt and bq. set ct_ctq_1, cq_cm_u_1,
!cq_cm_v_1, dtl_1, dqw_1, du , dv all to zero (explicit coupling)
!and boundary-layer depth.
bt_1(:,:) = 1.0 / tl_1_ij(:,:)
bq_1(:,:) = c_virtual / (1.0 + c_virtual * qw_1_ij(:,:))

!-----------------------------------------------------------------------
! Set date-related information for generate_anthropogenic_heat in sf_expl
!-----------------------------------------------------------------------

curr_day_number = day_of_year(                                                &
            current_time%year, current_time%month, current_time%day,          &
            l_360, l_leap)

z0m_scm(:,:)       = 0.0
z0h_scm(:,:)       = 0.0
z1_uv_top          = 0.0
z1_tq_top          = 0.0
ddmfx(:,:)         = 0.0
charnock_w(:,:)    = 0.0

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Initialise FLake variables required for melt calculations
!-----------------------------------------------------------------------

IF ( l_flake_model                                                            &
     .AND. ( .NOT. l_aggregate   )                                            &
     .AND. (land_pts > 0    ) ) THEN

  ! Calculate mean albedo of lake tile
  lake_albedo_gb(:) = 0.0
  DO point = 1, surft_pts(lake)
    l = ainfo%surft_index(point,lake)
    j = (ainfo%land_index(l) - 1) / row_length + 1
    i = ainfo%land_index(l) - (j-1) * row_length
    IF (sw_surft(l,lake) > EPSILON(sw_surft(l,lake))) THEN
      lake_albedo_gb(l) = 1.0 - (sw_surft(l,lake) / sw_down_ij(i,j))
    END IF
  END DO

  ! initialise FLake variables needed in surface exchange
  DO l = 1, land_pts
    ! set surface T based on T* of the lake tile
    lake_t_sfc_gb(l) = progs%tstar_surft(l, lake)

    ! set the FLake snow depth,
    ! BUT depending on the presence of ice
    lake_h_snow_gb(l) = 0.0
    IF (lake_h_ice_gb(l) > 0.0) THEN
      lake_h_snow_gb(l) = progs%snowdepth_surft(l, lake)
    END IF

    ! set the snow surface T based on T* of the lake tile
    lake_t_snow_gb(l) = progs%tstar_surft(l, lake)
  END DO

END IF ! l_flake_model

!-----------------------------------------------------------------------
!   Explicit calculations.
!-----------------------------------------------------------------------
CALL surf_couple_explicit(                                                    &
  !Misc INTENT(IN)
  bq_1, bt_1, photosynth_act_rad,                                             &
  curr_day_number,                                                            &
  !Forcing INTENT(IN)
  qw_1_ij, tl_1_ij, pstar_ij, lw_down_ij,                                     &
  !Fluxes INTENT(IN)
  sw_surft, tstar_ij,                                                         &
  !Diagnostics, INTENT(INOUT)
  sf_diag,                                                                    &
  !Fluxes INTENT(OUT)
  fqw_1_ij,ftl_1_ij, ftl_surft, fqw_surft, fqw_sicat, ftl_sicat, fsmc_pft,    &
  emis_surft,                                                                 &
  !Misc INTENT(OUT)
  radnet_sice, rhokm_1, rhokm_land, rhokm_ssi,                                &
  !Out of explicit and into implicit only INTENT(OUT)
  cdr10m,                                                                     &
  alpha1, alpha1_sea, alpha1_sice, ashtf_prime, ashtf_prime_sea,              &
  ashtf_prime_surft, epot_surft,                                              &
  fraca, resfs, resft, rhokh, rhokh_surft, rhokh_sice, rhokh_sea_ij,          &
  dtstar_ij_surft, dtstar_ij_sea, dtstar_ij_sice,                             &
  z0hssi, z0h_surft, z0mssi, z0m_surft, chr1p5m, chr1p5m_sice, canhc_surft,   &
  wt_ext_surft, flake,                                                        &
  !Out of explicit and into extra only INTENT(OUT)
  hcons_soilt,                                                                &
  !Out of explicit and into implicit and extra INTENT(OUT)
  tile_frac,                                                                  &
  !Additional arguments for the UM-----------------------------------------
  !JULES prognostics module
  cs_pool_gb_um,                                                              &
  ainfo%ti_cat_sicat,                                                         &
  !JULES ancil_info module
  land_pts, nsurft,                                                           &
  surft_pts,                                                                  &
  ! IN input data from the wave model
  charnock_w,                                                                 &
  !JULES coastal module
  flandg,                                                                     &
  !JULES aero module
  aerotype%co2_3d_ij,                                                         &
  !JULES trifctl module
  asteps_since_triffid, resp_s_acc_gb_um, resp_s_gb_um,                       &
  !JULES p_s_parms module
  soil_clay_ij,                                                               &
  !JULES switches module
  l_spec_z0,                                                                  &
  !Not in a JULES module
  numcycles, cycleno, z1_uv_top, z1_tq_top, ddmfx,                            &
  l_aero_classic, z0m_scm, z0h_scm, recip_l_mo_sea, rib,                      &
  flandfac, fseafac, fb_surf, u_s, t1_sd, q1_sd, rhostar,                     &
  vshr, resp_s_tot_soilt, emis_soil,                                          &
  !TYPES containing field data (IN OUT)
  crop_vars,psparms,ainfo,trif_vars,aerotype,urban_param,progs,trifctltype,   &
  coast, jules_vars)

!  Items that are only passed from explicit to implicit.
!  Cross-check for exlusive use in the UM too
!  rhokm_1, cdr10m, alpha1, alpha1_sice, ashtf_prime, ashtf_prime_surft,
!  epot_surft, fraca, rhokh, rhokh_surft, rhokh_sice, dtstar_ij_surft,
!  dtstar_ij_sea, dtstar_ij_sice, z0hssi, z0h_surft,
!  z0mssi, z0m_surft, chr1p5m, chr1p5m_sice, canhc_surft, wt_ext_surft, flake

!  Items that are passed from explicit to implicit and/or extra.
!  Cross-check in UM as above.
!  hcons_soilt, tile_frac

!Calculate variables required by sf_impl2
!In the UM, these are message passing variables and variables required by
!the implicit solver

!(i.e. no points are part land/part sea)
flandfac_u(:,:)   = 1.0
flandfac_v(:,:)   = 1.0
fseafac_u(:,:)    = 1.0
fseafac_v(:,:)    = 1.0

IF (l_ctile .AND. buddy_sea == on) THEN
  coast%taux_land_ij(:,:) =                                                   &
              rhokm_land(:,:) *  u_1_ij(:,:)                * flandfac_u(:,:)
  coast%taux_ssi_ij(:,:)  =                                                   &
              rhokm_ssi(:,:)  * (u_1_ij(:,:) - u_0_ij(:,:)) * fseafac_u(:,:)
  coast%tauy_land_ij(:,:) =                                                   &
              rhokm_land(:,:) *  v_1_ij(:,:)                * flandfac_v(:,:)
  coast%tauy_ssi_ij(:,:)  =                                                   &
              rhokm_ssi(:,:)  * (v_1_ij(:,:) - v_0_ij(:,:)) * fseafac_v(:,:)
ELSE   ! Standard code
  coast%taux_land_ij(:,:) = rhokm_land(:,:) *  u_1_ij(:,:)
  coast%taux_ssi_ij(:,:)  = rhokm_ssi(:,:)  * (u_1_ij(:,:) - u_0_ij(:,:))
  coast%tauy_land_ij(:,:) = rhokm_land(:,:) *  v_1_ij(:,:)
  coast%tauy_ssi_ij(:,:)  = rhokm_ssi(:,:)  * (v_1_ij(:,:) - v_0_ij(:,:))
END IF

taux_1_ij(:,:) = flandg(:,:) * coast%taux_land_ij(:,:) +                      &
                 (1.0 - flandg(:,:)) * coast%taux_ssi_ij(:,:)
tauy_1_ij(:,:) = flandg(:,:) * coast%tauy_land_ij(:,:) +                      &
                 (1.0 - flandg(:,:)) * coast%tauy_ssi_ij(:,:)

IF (sf_diag%suv10m_n) THEN
  sf_diag%cdr10m_n_u(:,:) = sf_diag%cdr10m_n(:,:)
  sf_diag%cdr10m_n_v(:,:) = sf_diag%cdr10m_n(:,:)
END IF

!-----------------------------------------------------------------------
! Implicit calculations.
!
! In the new boundary layer implicit solver in the UM, sf_impl2 is
! called twice - once with l_correct = .FALSE. and once with
! l_correct = .TRUE.
! The variables r_gamma1 and r_gamma2 that determine weights for this new
! solver are set so that the new scheme is the same as the old scheme
! (i.e. fully explicit coupling)
!-----------------------------------------------------------------------

! Set values of inputs for the implicit solver so that we get
! explicit coupling
du(:,:)   = 0.0
dv(:,:)   = 0.0
r_gamma1(:,:) = r_gamma

! To get explicit coupling with the new scheme, r_gamma2 is 0 everywhere
! for the 1st call, and 1 everywhere for the 2nd call
r_gamma2(:,:) = 0.0
! Call sf_impl2 with
l_correct = .FALSE.

ctctq1(:,:)     = 0.0
dqw1_1(:,:)     = 0.0
dtl1_1(:,:)     = 0.0
du_star1(:,:)   = 0.0
dv_star1(:,:)   = 0.0
cq_cm_u_1(:,:)  = 0.0
cq_cm_v_1(:,:)  = 0.0
rho1(:,:)       = 0.0
f3_at_p(:,:)    = 0.0
uStargbm(:,:)   = 0.0
flandg_u(:,:)   = flandg(:,:)
flandg_v(:,:)   = flandg(:,:)

CALL surf_couple_implicit(                                                    &
  !Important switch
  l_correct,                                                                  &
  !Forcing INTENT(IN)
  pstar_ij, lw_down_ij, qw_1_ij, tl_1_ij, u_1_ij, v_1_ij, u_0_ij, v_0_ij,     &
  !Fluxes INTENT(IN)
  sw_surft, emis_surft,                                                       &
  !Misc INTENT(IN) Many of these simply come out of explicit and into here.
  rhokm_1, rhokm_1, r_gamma1, r_gamma2, alpha1, alpha1_sea, alpha1_sice,      &
  ashtf_prime, ashtf_prime_sea, ashtf_prime_surft, du, dv, fraca, resfs,      &
  resft, rhokh, rhokh_surft, rhokh_sice, rhokh_sea_ij, z0hssi, z0mssi,        &
  z0h_surft, z0m_surft, chr1p5m, chr1p5m_sice, canhc_surft, flake, tile_frac, &
  wt_ext_surft, cdr10m, cdr10m, r_gamma,                                      &
  !Diagnostics, INTENT(INOUT)
  sf_diag,                                                                    &
  !Fluxes INTENT(INOUT)
  fqw_sicat, ftl_sicat, fqw_surft, fqw_1_ij, ftl_1_ij, ftl_surft,             &
  !Misc INTENT(INOUT)
  epot_surft, dtstar_ij_surft, dtstar_ij_sea, dtstar_ij_sice, radnet_sice,    &
  olr,                                                                        &
  !Fluxes INTENT(OUT)
  tstar_ij, le_surft, radnet_surft, e_sea_ij, h_sea_ij, taux_1_ij, tauy_1_ij, &
  ecan_surft, ei_ij,                                                          &
  esoil_ij_soilt, ext_soilt, snowmelt_ij, melt_surft,                         &
  ecan_ij, ei_surft, esoil_surft, sea_ice_htf_sicat, surf_ht_flux_ij,         &
  surf_htf_surft,                                                             &
  !Misc INTENT(OUT)
  error,                                                                      &
  !UM-only arguments
  !JULES ancil_info module
  nsurft, land_pts, surft_pts,                                                &
  !JULES prognostics module
  !JULES coastal module
  flandg, coast%tstar_land_ij, coast%tstar_sice_ij,                           &
  !JULES switches module
  l_co2_interactive, l_mr_physics,                                            &
  aerotype%co2_3d_ij,                                                         &
  !Arguments without a JULES module
  ctctq1,dqw1_1,dtl1_1,du_star1,dv_star1,cq_cm_u_1,cq_cm_v_1,flandg_u,flandg_v, &
  rho1, f3_at_p, uStarGBM,tscrndcl_ssi,tscrndcl_surft,tStbTrans,              &
  ei_sice,rhokh_mix, ti_gb,                                                   &
  !TYPES containing field data (IN OUT)
  crop_vars,ainfo,aerotype,progs, coast, jules_vars)


! Adjust the value of r_gamma2 to ensure explicit coupling
r_gamma2(:,:) = 1.0
! Call sf_impl2 again with
l_correct = .TRUE.

ctctq1(:,:)     = 0.0
dqw1_1(:,:)     = 0.0
dtl1_1(:,:)     = 0.0
du_star1(:,:)   = 0.0
dv_star1(:,:)   = 0.0
cq_cm_u_1(:,:)  = 0.0
cq_cm_v_1(:,:)  = 0.0
rho1(:,:)       = 0.0
f3_at_p(:,:)    = 0.0
uStargbm(:,:)   = 0.0
flandg_u(:,:)   = flandg(:,:)
flandg_v(:,:)   = flandg(:,:)


CALL surf_couple_implicit(                                                    &
  !Important switch
  l_correct,                                                                  &
  !Forcing INTENT(IN)
  pstar_ij, lw_down_ij, qw_1_ij, tl_1_ij, u_1_ij, v_1_ij, u_0_ij, v_0_ij,     &
  !Fluxes INTENT(IN)
  sw_surft, emis_surft,                                                       &
  !Misc INTENT(IN) Many of these simply come out of explicit and into here.
  rhokm_1, rhokm_1, r_gamma1, r_gamma2, alpha1, alpha1_sea, alpha1_sice,      &
  ashtf_prime, ashtf_prime_sea, ashtf_prime_surft, du, dv, fraca, resfs,      &
  resft, rhokh, rhokh_surft, rhokh_sice, rhokh_sea_ij, z0hssi, z0mssi,        &
  z0h_surft, z0m_surft, chr1p5m, chr1p5m_sice, canhc_surft, flake, tile_frac, &
  wt_ext_surft, cdr10m, cdr10m, r_gamma,                                      &
  !Diagnostics, INTENT(INOUT)
  sf_diag,                                                                    &
  !Fluxes INTENT(INOUT)
  fqw_sicat, ftl_sicat, fqw_surft, fqw_1_ij, ftl_1_ij, ftl_surft,             &
  !Misc INTENT(INOUT)
  epot_surft, dtstar_ij_surft, dtstar_ij_sea, dtstar_ij_sice, radnet_sice,    &
  olr,                                                                        &
  !Fluxes INTENT(OUT)
  tstar_ij, le_surft, radnet_surft, e_sea_ij, h_sea_ij, taux_1_ij, tauy_1_ij, &
  ecan_surft, ei_ij,                                                          &
  esoil_ij_soilt, ext_soilt, snowmelt_ij, melt_surft,                         &
  ecan_ij, ei_surft, esoil_surft, sea_ice_htf_sicat, surf_ht_flux_ij,         &
  surf_htf_surft,                                                             &
  !Misc INTENT(OUT)
  error,                                                                      &
  !UM-only arguments
  !JULES ancil_info module
  nsurft, land_pts, surft_pts,                                                &
  !JULES prognostics module
  !JULES coastal module
  flandg, coast%tstar_land_ij, coast%tstar_sice_ij,                           &
  !JULES switches module
  l_co2_interactive, l_mr_physics,                                            &
  aerotype%co2_3d_ij,                                                         &
  !Arguments without a JULES module
  ctctq1,dqw1_1,dtl1_1,du_star1,dv_star1,cq_cm_u_1,cq_cm_v_1,flandg_u,flandg_v, &
  rho1, f3_at_p, uStarGBM,tscrndcl_ssi,tscrndcl_surft,tStbTrans,              &
  ei_sice,rhokh_mix, ti_gb,                                                   &
  !TYPES containing field data (IN OUT)
  crop_vars,ainfo,aerotype,progs, coast, jules_vars)

!-------------------------------------------------------------------------------
! Calculation of pseudostresses.
!-------------------------------------------------------------------------------
IF (sf_diag%suv10m_n) THEN
  sf_diag%mu10m_n(:,:) = (taux_1_ij(:,:)) / sf_diag%cd10m_n(:,:)
  sf_diag%mv10m_n(:,:) = (tauy_1_ij(:,:)) / sf_diag%cd10m_n(:,:)
END IF
DEALLOCATE(sf_diag%cdr10m_n)
DEALLOCATE(sf_diag%cdr10m_n_u)
DEALLOCATE(sf_diag%cdr10m_n_v)
DEALLOCATE(sf_diag%cd10m_n)


!Arguments to surf_couple_extra that are not actively used in standalone JULES
!However, we expect a sensible value to be stored in them
ls_graup_ij(:,:) = 0.0
cca_2d(:,:) = 0.0
ls_rainfrac_land(:) = 0.0

  !Note that the commented intents may not be correct. The declarations in
  !surf_couple_extra will be more correct

CALL surf_couple_extra(                                                       &
  !Driving data and associated INTENT(IN)
  ls_rain_ij, con_rain_ij, ls_snow_ij, con_snow_ij, tl_1_ij, lw_down_ij,      &
  qw_1_ij, u_1_ij, v_1_ij, pstar_ij,                                          &
  !Fluxes INTENT(IN)
  ei_surft, surf_htf_surft, ecan_surft, ext_soilt, sw_surft,                  &
  !Misc INTENT(IN)
  timestep, sf_diag%smlt, tile_frac, hcons_soilt, rhostar,                    &
  !Fluxes INTENT(INOUT)
  melt_surft,                                                                 &
  !Misc INTENT(INOUT)
  lake_albedo_gb, lake_t_sfc_gb, lake_h_snow_gb, lake_t_snow_gb,              &
  !Fluxes INTENT(OUT)
  hf_snow_melt_gb, snowmelt_ij, snomlt_sub_htf_gb, sub_surf_roff_gb,          &
  surf_roff_gb, tot_tfall_gb, snow_melt_gb, rrun_gb, rflow_gb,                &
  snow_soil_htf,                                                              &
  !Arguments for the UM-----------------------------------------
  !IN
  land_pts, row_length, rows, river_row_length, river_rows,                   &
  ls_graup_ij,                                                                &
  cca_2d, nsurft, surft_pts,                                                  &
  lice_pts, soil_pts,                                                         &
  stf_sub_surf_roff,                                                          &
  toppdm%fexp_soilt, toppdm%gamtot_soilt, toppdm%ti_mean_soilt,               &
  toppdm%ti_sig_soilt, cs_ch4_soilt, flash_rate_ancil, pop_den_ancil,         &
  toppdm%a_fsat_soilt, toppdm%c_fsat_soilt, toppdm%a_fwet_soilt,              &
  toppdm%c_fwet_soilt,ntype, fqw_surft,                                       &
  delta_lambda, delta_phi, xx_cos_theta_latitude,                             &
  aocpl_row_length, aocpl_p_rows, xpa, xua, xva, ypa, yua, yva,               &
  g_p_field, g_r_field, n_proc, global_row_length, global_rows,               &
  global_river_row_length, global_river_rows, flandg,                         &
  trivdir, trivseq, r_area, slope, flowobs1, r_inext, r_jnext, r_land,        &
  trifctltype%frac_agr_gb, soil_clay_ij, trifctltype%resp_s_soilt,            &
  trifctltype%npp_gb, aerotype%u_s_std_surft,                                 &
  !INOUT
  a_steps_since_riv,                                                          &
  toppdm%fsat_soilt,                                                          &
  toppdm%fwetl_soilt,  toppdm%zw_soilt, toppdm%sthzw_soilt,                   &
  ls_rainfrac_land,                                                           &
  substore, surfstore, flowin, bflowin,                                       &
  tot_surf_runoff_gb, tot_sub_runoff_gb, acc_lake_evap_gb, twatstor,          &
  asteps_since_triffid, resp_s_acc_gb_um,                                     &
  cs_pool_gb_um,                                                              &
  inlandout_atm_gb,                                                           &
  !OUT
  dhf_surf_minus_soil, land_sea_mask,                                         &
  !TYPES containing field data (IN OUT)
  crop_vars,psparms,toppdm,fire_vars,ainfo,trif_vars,soilecosse, urban_param, &
  progs,trifctltype, jules_vars                                               &
  )

END SUBROUTINE control
#endif
