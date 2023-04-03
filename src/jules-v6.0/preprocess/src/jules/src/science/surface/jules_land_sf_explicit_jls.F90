! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE jules_land_sf_explicit_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                       &
                  ModuleName='JULES_LAND_SF_EXPLICIT_MOD'

CONTAINS
!  SUBROUTINE JULES_LAND_SF_EXPLICIT ---------------------------------
!
!  Purpose: Calculate explicit surface fluxes of heat, moisture and
!           momentum over land. Also calculates surface exchange 
!           coefficients required for implicit update of surface
!           fluxes and surface information required by the
!           explicit boundary layer routine
!
!
!  Documentation: UMDP 24.
!
!---------------------------------------------------------------------
!    Arguments :-
SUBROUTINE jules_land_sf_explicit (                                           &
! IN date-related values
 curr_day_number,                                                             &
! IN values defining field dimensions and subset to be processed :
 land_pts,                                                                    &
! IN  parameters for iterative SISL scheme
 numcycles, cycleno,                                                          &
! IN parameters required from boundary-layer scheme :
 bq_1,bt_1,z1_uv,z1_uv_top,z1_tq,z1_tq_top,qw_1,tl_1,                         &
! IN soil/vegetation/land surface data :
 land_index,nsurft,sm_levels,canopy,catch,catch_snow,hcon_soilt,              &
 ho2r2_orog, flandg,                                                          &
 snow_surft,sil_orog_land,smvccl_soilt,smvcst_soilt,smvcwt_soilt,sthf_soilt,  &
 sthu_soilt,z0_surft,z0h_surft_bare, z0m_soil_in,                             &
! IN input data from the wave model
 charnock_w,                                                                  &
! IN everything not covered so far :
 pstar,lw_down,sw_surft,zh,ddmfx,                                             &
 co2_mmr,co2_3d,l_co2_interactive,l_phenol,                                   &
 asteps_since_triffid,cs_pool_soilt,frac,canht_pft,photosynth_act_rad,lai_pft,&
 l_mr_physics,t_soil_soilt,tsurf_elev_surft,tstar,tstar_surft,z_land,         &
 albsoil_soilt,cos_zenith_angle,l_aero_classic,l_dust,l_dust_diag,            &
 clay_soilt,o3,                                                               &
! INOUT diagnostics
 sf_diag,                                                                     &
! INOUT data :
 gs,g_leaf_acc,npp_pft_acc,resp_w_pft_acc,resp_s_acc_soilt,                   &
 rhostar,fqw_1,ftl_1,t1_sd,q1_sd,vshr,vshr_land,                              &
! OUT Diagnostic not requiring STASH flags :
 ftl_surft,                                                                   &
! OUT variables for message passing
 rhokm_land, cdr10m,                                                          &
! OUT data required for mineral dust scheme
 u_s_std_surft,                                                               &
! OUT data required elsewhere in boundary layer or surface code
 alpha1,ashtf_prime_surft,fqw_surft,epot_surft,fraca,                         &
 resfs,resft,rhokh_surft,dtstar_surft,z0h_surft, z0m_surft,                   &
 chr1p5m,smc_soilt,hcons_soilt,gpp,npp,resp_p,g_leaf,gpp_pft,npp_pft,         &
 resp_p_pft,resp_s_soilt,resp_s_tot_soilt,resp_l_pft,resp_r_pft,              &
 resp_w_pft,n_leaf,n_root,n_stem,lai_bal,gc_surft,canhc_surft,wt_ext_surft,   &
 flake,surft_index,surft_pts,tile_frac,fsmc_pft,emis_surft,emis_soil,         &
! OUT required for classic aerosols
 cd_land,rib_surft,ch_surft_classic,cd_std_classic,                           &
! OUT required for sea and sea-ice calculations
 l_cdr10m_snow,                                                               &
 !New arguments replacing USE statements
 !Fluxes (OUT)
 t_growth_gb,                                                                 &
 !urban_param (IN)
 emisr_gb, emisw_gb, hwr_gb,                                                  &
 !jules_mod (IN OUT)
 albobs_scaling_surft,                                                        &
 !metstats (IN)
 metstats_prog,                                                               &
 !bvoc_vars (OUT)
 isoprene_gb, isoprene_pft, terpene_gb , terpene_pft,                         &
 methanol_gb, methanol_pft, acetone_gb, acetone_pft,                          &
 !trif_vars_mod (OUT)
 fapar_diag_pft, apar_diag_pft, apar_diag_gb, gpp_gb_acc, gpp_pft_acc,        &
 !crop_vars_mod (IN)
 rootc_cpft, sthu_irr_soilt, frac_irr_soilt, frac_irr_surft, dvi_cpft,        &
 !crop_vars_mod (IN OUT)
 resfs_irr_surft,                                                             &
 !crop_vars_mod (OUT)
 gs_irr_surft, smc_irr_soilt, wt_ext_irr_surft, gc_irr_surft,                 &
 !p_s_parms (IN) 
 bexp_soilt, sathh_soilt, v_close_pft, v_open_pft,                            &
 !urban_param (IN)
 wrr_gb,                                                                      &
 !Fluxes (IN OUT)
 anthrop_heat_surft,                                                          &
 !prognostics (IN)
 nsnow_surft, sice_surft, sliq_surft, snowdepth_surft,                        &
                        tsnow_surft, ds_surft,                                &
 !c_elevate (OUT)
 surf_hgt_surft, lw_down_elevcorr_surft,                                      &
 !jules_mod (OUT)
 snowdep_surft,                                                               &
 !urban_param (IN)
 hgt_gb, ztm_gb, disp_gb,                                                     &
 !lake_mod (IN)
 lake_t_ice_gb,lake_t_snow_gb,lake_t_mxl_gb,lake_t_mean_gb,                   &
 lake_h_snow_gb,lake_h_ice_gb,lake_h_mxl_gb,lake_depth_gb,                    &
 g_dt_gb,                                                                     &
 !lake_mod (OUT)
 nusselt_gb, ts1_lake_gb,                                                     &
 !ancil_info
 l_lice_point, l_soil_point,                                                  &
 !jules_surface_types (IN)
 diff_frac)

USE sf_flux_mod, ONLY: sf_flux
USE sfl_int_mod, ONLY: sfl_int
USE fcdch_mod, ONLY: fcdch
USE stdev1_mod, ONLY: stdev1
USE sf_rib_mod, ONLY: sf_rib
USE sf_resist_mod, ONLY: sf_resist
USE snowtherm_mod,         ONLY: snowtherm
USE gen_anthrop_heat_mod,  ONLY: generate_anthropogenic_heat
USE tilepts_mod,           ONLY: tilepts
USE heat_con_mod,          ONLY: heat_con
USE physiol_mod,           ONLY: physiol
USE qsat_mod, ONLY: qsat, qsat_mix
USE elevate_mod, ONLY: elevate
USE urbanz0_mod, ONLY: urbanz0
USE sf_orog_mod, ONLY: sf_orog
USE can_drag_mod, ONLY: can_drag_z0, can_drag_phi_m_h
USE calc_air_dens_mod, ONLY: calc_air_dens

!Use in relevant variables
USE theta_field_sizes, ONLY: t_i_length,t_j_length

USE dust_param, ONLY: ndiv, z0_soil
USE planet_constants_mod, ONLY: cp, vkman
USE csigma, ONLY: sbcon
USE jules_soil_biogeochem_mod, ONLY:                                          &
! imported scalar parameters
   soil_model_rothc,                                                          &
! imported scalar variables (IN)
   soil_bgc_model

USE jules_soil_mod, ONLY: dzsoil, dzsoil_elev, hcice, hcwat, hcondeep

USE jules_surface_types_mod, ONLY: npft, nnpft, ntype,                        &
                                   urban_canyon, urban_roof, soil, lake, ncpft

USE veg_param, ONLY: secs_per_360days

USE sf_diags_mod, ONLY: strnewsfdiag
USE timestep_mod, ONLY: timestep

USE ancil_info, ONLY:  dim_cs1, dim_cs2,                                      &
     co2_dim_len,co2_dim_row

USE ancil_info, ONLY: dim_cslayer, l_lice_surft, nsoilt, rad_nband

USE bl_option_mod, ONLY: l_quick_ap2
USE solinc_data, ONLY: sky, l_skyview

USE atm_fields_bounds_mod, ONLY:                                              &
   pdims_s, pdims, tdims

USE c_z0h_z0m, ONLY: z0h_z0m, z0h_z0m_classic
USE water_constants_mod, ONLY: lc, rho_ice, tm

USE jules_snow_mod, ONLY: cansnowtile                                         &
                          ,rho_snow_const                                     &
                          ,snow_hcon                                          &
                          ,l_snowdep_surf                                     &
                          ,l_snow_nocan_hc                                    &
                          ,nsmax                                              &
                          ,unload_rate_cnst                                   &
                          ,unload_rate_u


USE jules_surface_mod, ONLY: l_aggregate, formdrag, orog_drag_param,          &
                             fd_stab_dep, l_anthrop_heat_src,                 &
                             i_aggregate_opt, cor_mo_iter,                    &
                             use_correct_ustar, iscrntdiag,                   &
                             l_flake_model,l_elev_lw_down,                    &
                             effective_z0,                                    &
                             IP_ScrnDecpl2, IP_ScrnDecpl3,                    &
                             l_vary_z0m_soil, l_elev_land_ice

USE jules_vegetation_mod, ONLY: can_model, can_rad_mod, ilayers, l_triffid,   &
                                l_vegdrag_surft
                                
USE jules_irrig_mod, ONLY: l_irrig_dmd

USE jules_sea_seaice_mod, ONLY: l_ctile, charnock, ip_ss_solid

USE c_elevate, ONLY: l_elev_absolute_height

USE switches_urban, ONLY: l_moruses_rough, l_moruses_storage

USE jules_science_fixes_mod, ONLY: l_fix_ctile_orog, l_fix_wind_snow,         &
                                   l_accurate_rho

USE stochastic_physics_run_mod, ONLY: l_rp2, i_rp_scheme, i_rp2b, z0hm_pft_rp

USE urban_param_mod, ONLY: z0m_mat

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

!Use in TYPE definitions
USE metstats_mod, ONLY: metstats_prog_struct

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE ozone_vars, ONLY: flux_o3_pft, fo3_pft


IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs :-
!-----------------------------------------------------------------------
! (a) Defining horizontal grid and subset thereof to be processed.
!    Checked for consistency.
INTEGER, INTENT(IN) ::                                                        &
 curr_day_number,                                                             &
            ! IN current day of year
 land_pts,                                                                    &
            ! IN No of land points being processed.
 numcycles,                                                                   &
            ! Number of cycles (iterations) for iterative SISL.
 cycleno
            ! Iteration no

! Defining vertical grid of model atmosphere.
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
 bq_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                             ! IN A buoyancy parameter
                             !    (beta q tilde).
,bt_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                             ! IN A buoyancy parameter
                             !    (beta T tilde).
,z1_uv(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                             ! IN Height of lowest uv level (m).
,z1_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                             ! IN Height of lowest tq level (m).
                             !    Note, if the grid used is
                             !    staggered in the vertical,
                             !    Z1_UV and Z1_TQ can be
                             !    different.
,qw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                             ! IN Total water content
,tl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! IN Ice/liquid water temperature

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  charnock_w(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
! Charnock's coefficient from wave model

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
                    z1_uv_top(tdims%i_start:tdims%i_end,                      &
                              tdims%j_start:tdims%j_end)
                             ! Height of top of lowest uv-layer
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
                    z1_tq_top(tdims%i_start:tdims%i_end,                      &
                              tdims%j_start:tdims%j_end)
                             ! Height of top of lowest Tq-layer

! (c) Soil/vegetation/land surface parameters (mostly constant).
LOGICAL, INTENT(IN) ::                                                        &
 l_co2_interactive                                                            &
                             ! IN Switch for 3D CO2 field
,l_phenol 
                             ! IN Indicates whether phenology
                             !    in use

INTEGER, INTENT(IN) ::                                                        &
 land_index(land_pts)        ! IN LAND_INDEX(I)=J => the Jth
                             !    point in ROW_LENGTH,ROWS is the
                             !    land point.

INTEGER, INTENT(IN) ::                                                        &
 sm_levels                                                                    &
                             ! IN No. of soil moisture levels
,nsurft                                                                       &
                             ! IN No. of land-surface tiles
,asteps_since_triffid
                             ! IN Number of atmospheric
                             !    timesteps since last call
                             !    to TRIFFID.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
 canopy(land_pts,nsurft)                                                      &
                             ! IN Surface/canopy water for
                             !    snow-free land tiles (kg/m2)
,catch(land_pts,nsurft)                                                       &
                             ! IN Surface/canopy water capacity
                             !    of snow-free land tiles (kg/m2).
,catch_snow(land_pts,nsurft)                                                  &
                             ! IN Snow interception capacity of
                             !    tiles (kg/m2).
,hcon_soilt(land_pts,nsoilt)                                                  &
                             ! IN Soil thermal conductivity
                             !    (W/m/K).
,snow_surft(land_pts,nsurft)                                                  &
                             ! IN Lying snow on tiles (kg/m2)
,smvccl_soilt(land_pts,nsoilt,sm_levels)                                      &
                             ! IN Critical volumetric SMC
                             !    (cubic m per cubic m of soil).
,smvcst_soilt(land_pts,nsoilt,sm_levels)                                      &
                             ! IN Volumetric saturation point
                             !    (m3/m3 of soil).
,smvcwt_soilt(land_pts,nsoilt,sm_levels)                                      &
                             ! IN Volumetric wilting point
                             !    (cubic m per cubic m of soil).
,sthf_soilt(land_pts,nsoilt,sm_levels)                                        &
                             ! IN Frozen soil moisture content of
                             !    each layer as a fraction of
                             !    saturation.
,sthu_soilt(land_pts,nsoilt,sm_levels)                                        &
                             ! IN Unfrozen soil moisture content
                             !    of each layer as a fraction of
                             !    saturation.
,z0_surft(land_pts,nsurft)                                                    &
                             ! IN Tile roughness lengths (m).
,z0h_surft_bare(land_pts,nsurft)                                              &
                             ! IN Tile thermal roughness lengths
                             !    without snow cover(m).
,z0m_soil_in(land_pts)                                                        &
                             ! IN bare soil momentum z0 (m).
,sil_orog_land(land_pts)                                                      &
                             ! IN Silhouette area of unresolved
                             !    orography per unit horizontal
                             !    area on land points only.
,ho2r2_orog(land_pts)                                                         &
                             ! IN Standard Deviation of orography.
                             !    equivilent to peak to trough
                             !    height of unresolved orography
,flandg(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! IN Land fraction on all tiles.
                             !    divided by 2SQRT(2) on land
                             !    points only (m)

! (f) Atmospheric + any other data not covered so far, incl control.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
 pstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)                   &
                             ! IN Surface pressure (Pascals).
,lw_down(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                             ! IN Surface downward LW radiation
                             !    (W/m2).
,sw_surft(land_pts,nsurft)                                                    &
                             ! IN Surface net SW radiation on
                             !    land tiles (W/m2).
,zh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                      &
                             ! IN Height above surface of top of
                             !    boundary layer (metres).
,ddmfx(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                             ! IN Convective downdraught
                             !    mass-flux at cloud base
,co2_mmr                                                                      &
                             ! IN CO2 Mass Mixing Ratio
,co2_3d(co2_dim_len,co2_dim_row)                                              &
                             ! IN 3D CO2 field if required.
,cs_pool_soilt(land_pts,nsoilt,dim_cslayer,dim_cs1)                           &
                             ! IN Soil carbon (kg C/m2).
,frac(land_pts,ntype)                                                         &
                             ! IN Fractions of surface types.
,canht_pft(land_pts,npft)                                                     &
                             ! IN Canopy height (m)
,photosynth_act_rad(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)      &
                             ! IN Net downward shortwave radiation
                             !    in band 1 (w/m2).
,lai_pft(land_pts,npft)                                                       &
                             ! IN Leaf area index
,t_soil_soilt(land_pts,nsoilt,sm_levels)                                      &
                             ! IN Soil temperatures (K).
,tsurf_elev_surft(land_pts,nsurft)                                            &
                             ! IN Tiled ice sub-surface temperature (K)
,tstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                             ! IN GBM surface temperature (K).
,tstar_surft(land_pts,nsurft)                                                 &
                             ! IN Surface tile temperatures
,z_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                             ! IN Land height (m).
,albsoil_soilt(land_pts,nsoilt)                                               &
                             ! IN Soil albedo.
, cos_zenith_angle(tdims%i_start:tdims%i_end,                                 &
                   tdims%j_start:tdims%j_end)                                 &
                             ! IN Cosine of the zenith angle
,clay_soilt(land_pts,nsoilt,dim_cslayer)                                      &
                             ! IN Soil clay fraction.
,o3(land_pts)
                             ! IN Surface ozone concentration (ppb).

LOGICAL, INTENT(IN) ::                                                        &
 l_aero_classic                                                               &
                             ! IN switch for using CLASSIC aerosol
                             !    scheme
,l_dust                                                                       &
                             ! IN switch for mineral dust
,l_dust_diag                                                                  &
                             ! IN Switch for diagnostic mineral dust
                             !    lifting
,l_mr_physics
                             ! IN Switch for when mixing ratios are used
!-----------------------------------------------------------------------
!  In/outs :-
!-----------------------------------------------------------------------
!Diagnostics
TYPE (strnewsfdiag), INTENT(INOUT) :: sf_diag

REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
 gs(land_pts)                                                                 &
                              ! INOUT "Stomatal" conductance to
                                  !        evaporation (m/s).
,g_leaf_acc(land_pts,npft)                                                    &
                              ! INOUT Accumulated G_LEAF
,npp_pft_acc(land_pts,npft)                                                   &
                              ! INOUT Accumulated NPP_pft
,resp_w_pft_acc(land_pts,npft)                                                &
                              ! INOUT Accum RESP_W_pft
,resp_s_acc_soilt(land_pts,nsoilt,dim_cslayer,dim_cs1)                        &
                              ! INOUT Accumulated RESP_S
,rhostar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                              ! INOUT Surface air density
,fqw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                              ! INOUT Moisture flux between layers
                              !       (kg per square metre per sec).
                              !       FQW(,1) is total water flux
                              !       from surface, 'E'.
,ftl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                              ! INOUT FTL(,K) contains net turbulent
                              !       sensible heat flux into layer K
                              !       from below; so FTL(,1) is the
                              !       surface sensible heat, H.(W/m2)
,t1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                              ! OUT Standard deviation of turbulent
                              !     fluctuations of layer 1 temp;
                              !     used in initiating convection.
,q1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                              ! OUT Standard deviation of turbulent
                              !     flucs of layer 1 humidity;
                              !     used in initiating convection.
,vshr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                              ! OUT Magnitude of surface-to-lowest
                              !     atm level wind shear (m per s).
,vshr_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                              ! OUT Magnitude of surface-to-lowest
                              !     atm level wind shear (m per s).


!-----------------------------------------------------------------------
!  Outputs :-
!-----------------------------------------------------------------------
!-1 Diagnostic (or effectively so - includes coupled model requisites):-
INTEGER, INTENT(OUT) ::                                                       &
 surft_index(land_pts,ntype)                                                  &
                              ! OUT Index of tile points
,surft_pts(ntype)             ! OUT Number of tile points

!  (a) Calculated anyway (use STASH space from higher level) :-

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
 ftl_surft(land_pts,nsurft)                                                   &
                             ! OUT Surface FTL for land tiles
, u_s_std_surft(land_pts,nsurft)                                              &
                             ! OUT Surface friction velocity
                             !     (standard value)
                             !     for mineral dust
,emis_surft(land_pts,nsurft)                                                  &
                             ! OUT Emissivity for land tiles
,emis_soil(land_pts)
                             ! OUT Emissivity of underlying soil

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
 rhokm_land(pdims_s%i_start:pdims_s%i_end,                                    &
            pdims_s%j_start:pdims_s%j_end),                                   &
 cdr10m(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)

!  (b) Not passed between lower-level routines (not in workspace at this
!      level) :-

!-2 Genuinely output, needed by other atmospheric routines :-

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
 alpha1(land_pts,nsurft)                                                      &
                             ! OUT Mean gradient of saturated
                             !     specific humidity with respect
                             !     to temperature between the
                             !     bottom model layer and tile
                             !     surfaces
,ashtf_prime_surft(land_pts,nsurft)                                           &
                             ! OUT Coefficient to calculate
                             !     surface heat flux into land
                             !     tiles.
,fqw_surft(land_pts,nsurft)                                                   &
                             ! OUT Surface FQW for land tiles
,epot_surft(land_pts,nsurft)                                                  &
                             ! OUT Local EPOT for land tiles.
,fraca(land_pts,nsurft)                                                       &
                             ! OUT Fraction of surface moisture
                             !     flux with only aerodynamic
                             !     resistance for snow-free land
                             !     tiles.
,resfs(land_pts,nsurft)                                                       &
                             ! OUT Combined soil, stomatal
                             !     and aerodynamic resistance
                             !     factor for fraction (1-FRACA)
                             !     of snow-free land tiles.
,resft(land_pts,nsurft)                                                       &
                             ! OUT Total resistance factor.
                             !     FRACA+(1-FRACA)*RESFS for
                             !     snow-free land, 1 for snow.
,rhokh_surft(land_pts,nsurft)                                                 &
                             ! OUT Surface exchange coefficients
                             !     for land tiles
,dtstar_surft(land_pts,nsurft)                                                &
                             ! OUT Change in TSTAR over timestep
                             !     for land tiles
,z0h_surft(land_pts,nsurft)                                                   &
                             ! OUT Tile roughness lengths for heat
                             !     and moisture (m).
,z0m_surft(land_pts,nsurft)                                                   &
                             ! OUT Tile roughness lengths for
                             !     momentum.
,chr1p5m(land_pts,nsurft)                                                     &
                             ! OUT Ratio of coefffs for
                             !     calculation of 1.5m temp for
                             !     land tiles.
,smc_soilt(land_pts,nsoilt)                                                   &
                             ! OUT Available moisture in the
                             !     soil profile (mm).
,hcons_soilt(land_pts,nsoilt)                                                 &
                             ! OUT Soil thermal conductivity
                             !     including water and ice
,gpp(land_pts)                                                                &
                             ! OUT Gross primary productivity
                             !     (kg C/m2/s).
,npp(land_pts)                                                                &
                             ! OUT Net primary productivity
                             !     (kg C/m2/s).
,resp_p(land_pts)                                                             &
                             ! OUT Plant respiration (kg C/m2/s).
,g_leaf(land_pts,npft)                                                        &
                             ! OUT Leaf turnover rate (/360days).
,gpp_pft(land_pts,npft)                                                       &
                             ! OUT Gross primary productivity
                             !     on PFTs (kg C/m2/s).
,npp_pft(land_pts,npft)                                                       &
                             ! OUT Net primary productivity
                             !     (kg C/m2/s).
,resp_p_pft(land_pts,npft)                                                    &
                             ! OUT Plant respiration on PFTs
                             !     (kg C/m2/s).
,resp_s_soilt(land_pts,nsoilt,dim_cslayer,dim_cs1)                            &
                          ! OUT Soil respiration (kg C/m2/s).
,resp_s_tot_soilt(dim_cs2,nsoilt)                                             &
                             ! OUT Total soil respiration
                             ! (kg C/m2/s).
,resp_l_pft(land_pts,npft)                                                    &
                             ! OUT Leaf maintenance respiration
                             !     (kg C/m2/s).
,resp_r_pft(land_pts,npft)                                                    &
                             ! OUT Root maintenance respiration
                             !     (kg C/m2/s).
,resp_w_pft(land_pts,npft)                                                    &
                             ! OUT Wood maintenance respiration
                             !     (kg C/m2/s).
,n_leaf(land_pts,npft)                                                        &
                             ! OUT Leaf N content scaled by LAI
                             !     (kg N/m2).
,n_root(land_pts,npft)                                                        &
                             ! OUT Root N content scaled by LAI_bal
                             !     (kg N/m2).
,n_stem(land_pts,npft)                                                        &
                             ! OUT Stem N content scaled by LAI_bal
                             !     (kg N/m2).
,lai_bal(land_pts,npft)                                                       &
                             ! OUT LAI_bal
,gc_surft(land_pts,nsurft)                                                    &
                             ! OUT "Stomatal" conductance to
                             !      evaporation for land tiles
                             !      (m/s).
,canhc_surft(land_pts,nsurft)                                                 &
                             ! OUT Areal heat capacity of canopy
                             !    for land tiles (J/K/m2).
,wt_ext_surft(land_pts,sm_levels,nsurft)                                      &
                             ! OUT Fraction of evapotranspiration
                             !     which is extracted from each
                             !     soil layer by each tile.
,flake(land_pts,nsurft)                                                       &
                             ! OUT Lake fraction.
,tile_frac(land_pts,nsurft)                                                   &
                             ! OUT Tile fractions including
                             !     snow cover in the ice tile.
,fsmc_pft(land_pts,npft)     ! OUT Moisture availability factor.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
 cd_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                             ! OUT Bulk transfer coefficient for
                             !     momentum over land.
,rib_surft(land_pts,nsurft)                                                   &
                             ! OUT RIB for land tiles.
,ch_surft_classic(land_pts,nsurft)                                            &
                             ! OUT Bulk transfer coefficient for
                             !     heat for aerosol deposition.
,cd_std_classic(land_pts,nsurft)
                             ! OUT Bulk transfer coefficient for
                             !     momentum for aerosol deposition.
LOGICAL, INTENT(OUT) ::                                                       &
 l_cdr10m_snow
                             ! OUT Flag indicating if cdr10m 
                             !     (an interpolation coefficient) is
                             !     to be calculated for use with 
                             !     snow unloading.

!New arguments replacing USE statements
!urban_param
REAL(KIND=real_jlslsm), INTENT(IN) :: emisr_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: emisw_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: hwr_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: wrr_gb(land_pts)

!jules_mod
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                     &
          albobs_scaling_surft(land_pts,ntype,rad_nband)

!p_s_parms (IN) 
REAL(KIND=real_jlslsm), INTENT(IN) :: bexp_soilt(land_pts,nsoilt,sm_levels)
REAL(KIND=real_jlslsm), INTENT(IN) :: sathh_soilt(land_pts,nsoilt,sm_levels)
REAL(KIND=real_jlslsm), INTENT(IN) :: v_close_pft(land_pts,sm_levels,npft)
REAL(KIND=real_jlslsm), INTENT(IN) :: v_open_pft(land_pts,sm_levels,npft)

!crop_vars_mod (IN)
REAL(KIND=real_jlslsm), INTENT(IN) :: rootc_cpft(land_pts,ncpft)
REAL(KIND=real_jlslsm), INTENT(IN) :: sthu_irr_soilt(land_pts,nsoilt,sm_levels)
REAL(KIND=real_jlslsm), INTENT(IN) :: frac_irr_soilt(land_pts,nsoilt)
REAL(KIND=real_jlslsm), INTENT(IN) :: frac_irr_surft(land_pts,nsurft)
REAL(KIND=real_jlslsm), INTENT(IN) :: dvi_cpft(land_pts,ncpft)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: resfs_irr_surft(land_pts,nsurft)

!metstats (IN)
TYPE(metstats_prog_struct), INTENT(IN) :: metstats_prog(land_pts)

!Fluxes (IN OUT)
REAL(KIND=real_jlslsm), INTENT(IN OUT)    :: anthrop_heat_surft(land_pts,nsurft)

!Fluxes
REAL(KIND=real_jlslsm), INTENT(OUT) :: t_growth_gb(land_pts)

!bvoc_vars
REAL(KIND=real_jlslsm), INTENT(OUT) :: isoprene_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: terpene_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: methanol_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: acetone_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: isoprene_pft(land_pts,npft)
REAL(KIND=real_jlslsm), INTENT(OUT) :: terpene_pft(land_pts,npft)
REAL(KIND=real_jlslsm), INTENT(OUT) :: methanol_pft(land_pts,npft)
REAL(KIND=real_jlslsm), INTENT(OUT) :: acetone_pft(land_pts,npft)

!trif_vars_mod
REAL(KIND=real_jlslsm), INTENT(OUT) :: fapar_diag_pft(land_pts,npft)
REAL(KIND=real_jlslsm), INTENT(OUT) :: apar_diag_pft(land_pts,npft)
REAL(KIND=real_jlslsm), INTENT(OUT) :: apar_diag_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: gpp_gb_acc(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: gpp_pft_acc(land_pts,npft)

!crop_vars_mod (OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) :: gs_irr_surft(land_pts,nsurft)
REAL(KIND=real_jlslsm), INTENT(OUT) :: smc_irr_soilt(land_pts,nsoilt)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
        wt_ext_irr_surft(land_pts,sm_levels,nsurft)
REAL(KIND=real_jlslsm), INTENT(OUT) :: gc_irr_surft(land_pts,nsurft)

!prognostics (IN)
INTEGER, INTENT(IN) :: nsnow_surft(land_pts,nsurft)
REAL(KIND=real_jlslsm), INTENT(IN) :: sice_surft(land_pts,nsurft,nsmax),      &
                                      sliq_surft(land_pts,nsurft,nsmax),      &
                                      snowdepth_surft(land_pts,nsurft),       &
                                      tsnow_surft(land_pts,nsurft,nsmax),     &
                                      ds_surft(land_pts,nsurft,nsmax)

!c_elevate (OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) :: surf_hgt_surft(land_pts,nsurft),       &
                                       lw_down_elevcorr_surft(land_pts,nsurft)

!jules_mod (OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) :: snowdep_surft(land_pts,nsurft)

!urban_param (IN)
REAL(KIND=real_jlslsm), INTENT(IN) :: hgt_gb(land_pts),                       &
                                      ztm_gb(land_pts),                       &
                                      disp_gb(land_pts)

!lake_mod (IN)
REAL(KIND=real_jlslsm), INTENT(IN) :: lake_t_ice_gb(land_pts),                &
                                      lake_t_snow_gb(land_pts),               &
                                      lake_t_mxl_gb(land_pts),                &
                                      lake_t_mean_gb(land_pts),               &
                                      lake_h_snow_gb(land_pts),               &
                                      lake_h_ice_gb(land_pts),                &
                                      lake_h_mxl_gb(land_pts),                &
                                      lake_depth_gb(land_pts),                &
                                      g_dt_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: ts1_lake_gb(land_pts),                 &
                                      nusselt_gb(land_pts)

!ancil_info (IN)                                      
LOGICAL, INTENT(IN) :: l_lice_point(land_pts)
LOGICAL, INTENT(IN) :: l_soil_point(land_pts)

!JULES surface_types_mod (IN)
REAL(KIND=real_jlslsm), INTENT(IN) :: diff_frac(t_i_length * t_j_length)

!-----------------------------------------------------------------------
! LOCAL variables
!-----------------------------------------------------------------------
!  Workspace :-
REAL(KIND=real_jlslsm) :: work_clay         ! working variable

REAL(KIND=real_jlslsm) ::                                                     &
 vfrac_surft(land_pts,nsurft)                                                 &
                             ! Fractional canopy coverage for
                             ! land tiles.
,radnet_surft(land_pts,nsurft)                                                &
                             ! Surface net radiation on tiles
,csnow(land_pts,nsmax)                                                        &
                             ! Areal heat capacity of snow (J/K/m2)
,ksnow(land_pts,nsmax)                                                        &
                             ! Thermal conductivity of snow (W/m/K)
,hcons_snow(land_pts,nsurft)                                                  &
                             ! Snow thermal conductivity
,resp_frac(dim_cs2,dim_cslayer)                                               &
                             ! respired fraction of RESP_S
,gc_stom_surft(land_pts,nsurft)
                             ! canopy conductance

REAL(KIND=real_jlslsm) ::                                                     &
 lh0                         ! Latent heat for snow free surface
                             !   =LS for sea-ice, =LC otherwise


!  Workspace for sea-ice and marginal ice zone
REAL(KIND=real_jlslsm) ::                                                     &
 ch_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! Bulk transfer coefficient for
                             !      het and moisture over land.

!  Workspace for land tiles
REAL(KIND=real_jlslsm) ::                                                     &
 cd_std(land_pts,nsurft)                                                      &
                             ! Local drag coefficient for calc
                             ! of interpolation coefficient
,cd_surft(land_pts,nsurft)                                                    &
                             ! Drag coefficient
,ch_surft(land_pts,nsurft)                                                    &
                             ! Transfer coefficient for heat and
                             ! moisture
,chn(land_pts,nsurft)                                                         &
                             ! Neutral value of CH.
,dq(land_pts)                                                                 &
                             ! Sp humidity difference between
                             ! surface and lowest atmospheric lev
,epdt(land_pts)                                                               &
                             ! "Potential" Evaporation * Timestep
,pstar_land(land_pts)                                                         &
                             ! Surface pressure for land points.
,qstar_surft(land_pts,nsurft)                                                 &
                             !Surface saturated sp humidity.
,rhokh_can(land_pts,nsurft)                                                   &
                             ! Exchange coefficient for canopy
                             ! air to surface
,rhokm_1_surft(land_pts,nsurft)                                               &
                             ! Momentum exchange coefficient.
,tsurf(land_pts)                                                              &
                             ! Surface layer temp (snow or soil) (
,dzsurf(land_pts)                                                             &
                             ! Surface layer thickness
                             ! (snow or soil) (m)
,canhc_surf(land_pts)                                                         &
                             ! Surface layer thickness
                             ! (snow or soil) (m)
,hcons_surf(land_pts)                                                         &
                             ! Thermal conductivity
                             ! (snow or soil)  (W/m/K)
,wind_profile_factor(land_pts,nsurft)                                         &
                             ! For transforming effective surface
                             ! transfer coefficients to those
                             ! excluding form drag.
,z0m_eff_surft(land_pts,nsurft)                                               &
                             ! Effective momentum roughness length
,db_surft(land_pts,nsurft)                                                    &
                             ! Buoyancy difference for surface
                             ! tile
,v_s_surft(land_pts,nsurft)                                                   &
                             ! Surface layer scaling velocity
                             ! for tiles (m/s).
,v_s_std(land_pts,nsurft)                                                     &
                             ! Surface layer scaling velocity
                             ! for tiles excluding orographic
                             ! form drag (m/s).
,u_s_iter_surft(land_pts,nsurft)                                              &
                             ! Scaling velocity from middle of
                             ! MO scheme - picked up in error by
                             ! dust code!
,recip_l_mo_surft(land_pts,nsurft)                                            &
                             ! Reciprocal of the Monin-Obukhov
                             ! length for tiles (m^-1).
,z0m_soil(land_pts,nsurft)                                                    &
                             ! Bare soil momentum roughness length
                             ! for use in 1 tile dust scheme
,z0h_soil(land_pts,nsurft)                                                    &
                             ! Bare soil roughness length for heat
                             ! for use in 1 tile dust scheme
,wind_profile_fac_soil(land_pts,nsurft)                                       &
                             ! Equivalent of wind_profile_factor
                             ! for use in 1 tile dust scheme
,cd_surft_soil(land_pts,nsurft), ch_surft_soil(land_pts,nsurft)               &
,cd_std_soil(land_pts,nsurft), v_s_surft_soil(land_pts,nsurft)                &
,recip_l_mo_surft_soil(land_pts,nsurft)                                       &
                             ! Dummy output variables from extra
                             ! call to fcdch needed for
                             ! 1 tile dust scheme
,v_s_std_soil(land_pts,nsurft)                                                &
                             ! Bare soil surface layer scaling
                             ! velocity for tiles excluding
                             ! orographic form drag (m/s)
                             ! for use in 1 tile dust scheme
,u_s_iter_soil(land_pts,nsurft)                                               &
                             ! Bare soil scaling velocity from
                             ! middle of MO scheme for use in
                             ! 1 tile dust scheme - picked up in
                             ! error by dust code
,z0h_surft_classic(land_pts,nsurft)                                           &
                             ! z0h to be used in calculation for
                             ! CLASSIC aerosol deposition
,cd_surft_classic(land_pts,nsurft)                                            &
,v_s_surft_classic(land_pts,nsurft)                                           &
,recip_l_mo_surft_classic(land_pts,nsurft)                                    &
,v_s_std_classic(land_pts,nsurft)                                             &
,u_s_iter_classic(land_pts,nsurft)                                            &
                             ! Dummy output variables from extra
                             ! call to fcdch needed for aerosol
                             ! deposition with different z0h
,t_elev(land_pts,nsurft)                                                      &
                             ! Temperature at elevated height (k)
,q_elev(land_pts,nsurft)                                                      &
                             ! Specific humidity at elevated
                             !     height (kg per kg air)
,qs1_elev(land_pts,nsurft)                                                    &
                             ! Saturated specific humidity at elev
                             !     height (kg per kg air)
,scaling_urban(land_pts,nsurft)                                               &
                             ! MORUSES: ground heat flux scaling;
                             ! canyon tile only coupled to soil
,hcon_lake(land_pts)                                                          &
                             ! "average" thermal conductivity of
                             ! lake-ice, lake and soil sandwich (W/m/K).
,zdt_surft(land_pts,nsurft)                                                   &
                             ! Difference between the canopy height and
                             ! displacement height (m)
,phi_m(land_pts)                                                              &
                             ! Monin-Obukhov stability function for momentum
                             ! integrated to the model's lowest wind level.
,phi_h(land_pts)                                                              &
                             ! Monin-Obukhov stability function for scalars
                             ! integrated to the model's lowest temperature
                             ! and humidity level.
,rhostar_mom(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! Surface air density for momentum

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
 tau_surft(:,:)
                             ! Local tau_1 for land tiles.

! dummy arrays required for sea and se-ice to create universal
! routines for all surfaces
REAL(KIND=real_jlslsm) ::                                                     &
 array_zero(t_i_length * t_j_length)                                          &
                                ! Array of zeros
,zdt_dummy(t_i_length * t_j_length)
                                ! Dummy array for zdt

!Gridbox mean values calculated from soil tiled versions for FLAKE
REAL(KIND=real_jlslsm) ::                                                     &
  hcons_mean_soil(land_pts),                                                  &
  tsoil_mean_soil(land_pts)


!  Local scalars :-

INTEGER ::                                                                    &
 i,j                                                                          &
             ! Loop counter (horizontal field index).
,k                                                                            &
             ! Loop counter (tile field index).
,l                                                                            &
             ! Loop counter (land point field index).
,n                                                                            &
             ! Loop counter (tile index).
,nn                                                                           &
             ! Loop counter (soil carbon layers)
,m                                                                            &
             ! Index for soil tile
,n_veg
             ! Actual or dummy pointer to array
             ! defined only on PFTs

REAL(KIND=real_jlslsm) ::                                                     &
 zetam                                                                        &
             ! Temporary in calculation of CHN.
,zetah                                                                        &
             ! Temporary in calculation of CHN.
,zeta1                                                                        &
             ! Work space
,z0                                                                           &
             ! yet more workspace
,ws10                                                                         &
             ! 10-m wind speed, used in calculation of unloading of snow from
             ! vegetation (m s-1)

             ! Temporary variables for adjustment of downwelling
             ! logwave to elevation tiles and correction back to
             ! conserve gridbox mean
,t_rad                                                                        &
,lw_down_surftsum                                                             &
,lw_down_surftabs

LOGICAL :: l_vegdrag_active_here
             ! Logical to indicate whether the vegetative drag scheme
             ! is active on the current surface tile in cases where
             ! l_vegdrag_surft itself may not be applicable

REAL(KIND=real_jlslsm) :: sea_point

INTEGER :: n_diag

INTEGER :: errcode
CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='JULES_LAND_SF_EXPLICIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


!-----------------------------------------------------------------------
!  0. Initialisations
!-----------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(I,j)                                     &
!$OMP SHARED(t_i_length,t_j_length,array_zero,tdims,rhokm_land,cd_land,       &
!$OMP ch_land,zdt_dummy)
!$OMP DO SCHEDULE(STATIC)
DO i = 1,t_i_length * t_j_length
  array_zero(i)     = 0.0
  zdt_dummy(i)      = 0.0
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    rhokm_land(i,j) = 0.0
    cd_land(i,j) = 0.0
    ch_land(i,j) = 0.0
  END DO
END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

IF (sf_diag%l_tau_surft .OR. sf_diag%l_tau_1) THEN
  ALLOCATE(tau_surft(land_pts,nsurft))
ELSE
  ALLOCATE(tau_surft(1,nsurft))
END IF

!-----------------------------------------------------------------------
!  1. Initialise FTL_SURFT and RIB_SURFT on all tiles at all points,
!     to allow STASH to process these as diagnostics.
!-----------------------------------------------------------------------
!$OMP PARALLEL PRIVATE(l,n) DEFAULT(NONE)                                 &
!$OMP SHARED(land_pts,nsurft,ftl_surft,fqw_surft,tau_surft,               &
!$OMP rib_surft,z0m_surft,u_s_std_surft,chr1p5m,resfs,radnet_surft,       &
!$OMP snowdep_surft,snowdepth_surft,sf_diag)
DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
  DO l = 1,land_pts
    ! MORUSES Initialise urban roughness array
    ftl_surft(l,n) = 0.0
    fqw_surft(l,n) = 0.0
    rib_surft(l,n) = 0.0
    z0m_surft(l,n) = 0.0
    u_s_std_surft(l,n) = 0.0
    chr1p5m(l,n) = 0.0
    resfs(l,n) = 0.0
    radnet_surft(l,n)  = 0.0
    ! Equivalent snowdepth for surface calculations.
    snowdep_surft(l,n) = snowdepth_surft(l,n)
    IF (sf_diag%l_et_stom .OR. sf_diag%l_et_stom_surft) THEN
      sf_diag%resfs_stom(l,n) = 0.0
    END IF
  END DO
!$OMP END DO NOWAIT
END DO
IF (sf_diag%l_tau_surft) THEN
  DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
    DO l = 1,land_pts
      tau_surft(l,n) = 0.0
    END DO
!$OMP END DO NOWAIT
  END DO
END IF
!$OMP END PARALLEL


!-----------------------------------------------------------------------
! Call TILEPTS to calculate surft_pts and surft_index for surface types
!-----------------------------------------------------------------------
CALL tilepts(land_pts,frac,surft_pts,surft_index, l_lice_point)


!-----------------------------------------------------------------------
!   Generate the anthropogenic heat for surface calculations
!-----------------------------------------------------------------------
CALL generate_anthropogenic_heat( curr_day_number, land_pts, frac,            &
                                  l_anthrop_heat_src,                         &
                                  !New arguments replacing USE statements
                                  !urban_param (IN)
                                  wrr_gb,                                     &
                                  !Fluxes (IN OUT)
                                  anthrop_heat_surft,                         &
                                  !Ancil_info (IN)
                                  l_lice_point)

!-----------------------------------------------------------------------
! Call physiology routine to calculate surface conductances and carbon
! fluxes.
!-----------------------------------------------------------------------
CALL physiol (                                                                &
  land_pts,land_index,                                                        &
  sm_levels,nsurft,surft_pts,surft_index,                                     &
  dim_cs1,dim_cs2,                                                            &
  co2_mmr,co2_3d,co2_dim_len, co2_dim_row,l_co2_interactive,                  &
  can_model,cs_pool_soilt,frac,canht_pft,photosynth_act_rad,                  &
  lai_pft,pstar,qw_1,sthu_soilt,sthf_soilt,t_soil_soilt,tstar_surft,          &
  smvccl_soilt,smvcst_soilt,smvcwt_soilt,vshr,z0_surft,z1_uv,o3,              &
  canhc_surft,vfrac_surft,emis_surft,emis_soil,flake,                         &
  g_leaf,gs,gc_surft,gc_stom_surft,gpp,gpp_pft,npp,npp_pft,                   &
  resp_p,resp_p_pft,resp_s_soilt,resp_l_pft,                                  &
  resp_r_pft,resp_w_pft,n_leaf,                                               &
  n_root,n_stem,lai_bal,                                                      &
  smc_soilt,wt_ext_surft,fsmc_pft,                                            &
  albsoil_soilt,cos_zenith_angle,                                             &
  can_rad_mod,ilayers,flux_o3_pft,fo3_pft,sf_diag,asteps_since_triffid,       &
  !New arguments replacing USE statements
  !Fluxes (OUT)
  t_growth_gb,                                                                &
  !urban_param (IN)
  emisr_gb, emisw_gb, hwr_gb,                                                 &
  !jules_mod (IN OUT)
  albobs_scaling_surft,                                                       &
  !metstats (IN)
  metstats_prog,                                                              &
  !bvoc_vars (OUT)
  isoprene_gb, isoprene_pft, terpene_gb , terpene_pft,                        &
  methanol_gb, methanol_pft, acetone_gb, acetone_pft,                         &
  !trif_vars_mod (OUT)
  fapar_diag_pft, apar_diag_pft, apar_diag_gb, gpp_gb_acc, gpp_pft_acc,       &
  !crop_vars_mod (IN)
  rootc_cpft, sthu_irr_soilt, frac_irr_soilt, frac_irr_surft, dvi_cpft,       &
  !crop_vars_mod (OUT)
  gs_irr_surft, smc_irr_soilt, wt_ext_irr_surft, gc_irr_surft,                &
  !p_s_parms (IN) 
  bexp_soilt, sathh_soilt, v_close_pft, v_open_pft,                           &
  !ancil_info
  l_soil_point,                                                               &
  !jules_surface_types (IN)
  diff_frac)

!----------------------------------------------------------------------
! If TRIFFID is being used apply any correction to the land-atmosphere
! fluxes on the first timestep after the last TRIFFID call. Such a
! correction will typically be associated with a total depletion of
! carbon or with maintanence of the seed fraction. The corrections
! are stored in the accumulation variables after the call to TRIFFID.
! The correction is added to the instantaneous land-atmosphere fluxes
! (so that the atmospheric carbon budget is corrected) but is not
! included in the accumulation variables which drive TRIFFID, since
! this has already been dealt with during the last TRIFFID call.
!----------------------------------------------------------------------
IF (l_triffid .AND. (asteps_since_triffid == 1)                               &
    .AND. ( cycleno == numcycles .OR. l_quick_ap2) ) THEN

  DO n = 1,nnpft
    DO l = 1,land_pts
      npp_pft(l,n) = npp_pft(l,n) + npp_pft_acc(l,n) / timestep
      resp_p_pft(l,n) = resp_p_pft(l,n) - npp_pft_acc(l,n) / timestep
      npp_pft_acc(l,n)=-npp_pft_acc(l,n)
    END DO
  END DO

  ! Here we have assumed that RothC must be used with TRIFFID, and is called
  ! on the same timestep.
  IF ( soil_bgc_model == soil_model_rothc ) THEN
    DO n = 1,dim_cs1
      DO l = 1,land_pts
        DO nn = 1,dim_cslayer
          !soil tiling is not compatible with triffid. OK to hard-code soil
          !tile index to 1 here
          resp_s_soilt(l,1,nn,n) = resp_s_soilt(l,1,nn,n)                     &
                                   + (resp_s_acc_soilt(l,1,nn,n) / timestep)
          resp_s_acc_soilt(l,1,nn,n) = -resp_s_acc_soilt(l,1,nn,n)
        END DO
      END DO
    END DO
  END IF

END IF

!----------------------------------------------------------------------
! Increment accumulation of leaf turnover rate.
! This is required for leaf phenology and/or TRIFFID, either of
! which can be enabled independently of the other.
!----------------------------------------------------------------------
IF ( cycleno == numcycles .OR. l_quick_ap2 ) THEN

  IF (l_phenol .OR. l_triffid) THEN
    DO n = 1,nnpft
      DO l = 1,land_pts
        g_leaf_acc(l,n) = g_leaf_acc(l,n) +                                   &
                          g_leaf(l,n) * ( timestep / secs_per_360days )
      END DO
    END DO
  END IF

  !----------------------------------------------------------------------
  ! Increment accumulation prognostics for TRIFFID
  !----------------------------------------------------------------------
  IF (l_triffid) THEN
    DO n = 1,nnpft
      DO l = 1,land_pts
        npp_pft_acc(l,n) = npp_pft_acc(l,n) + npp_pft(l,n) * timestep
        resp_w_pft_acc(l,n) = resp_w_pft_acc(l,n)                             &
                              + resp_w_pft(l,n) * timestep
      END DO
    END DO
  END IF

  IF ( soil_bgc_model == soil_model_rothc ) THEN
    DO n = 1,dim_cs1
      DO l = 1,land_pts
        DO nn = 1,dim_cslayer
          !soil tiling is not compatible with triffid. OK to hard-code soil
          !tile index to 1 here
          resp_s_acc_soilt(l,1,nn,n) = resp_s_acc_soilt(l,1,nn,n)             &
                                       + resp_s_soilt(l,1,nn,n) * timestep
        END DO
      END DO
    END DO
  END IF

END IF ! CycleNo == NumCycles

!-----------------------------------------------------------------------
! calculate CO2:(BIO+HUM) ratio, dependent on soil clay content, and
! sum soil respiration components
! (RESP_FRAC here then contains the fraction of soil respiration which
! is respired to the atmos. the rest is re-partitioned into BIO+HUM)

! resp_s_acc_soilt contains the full amount, and this is carried forward to
! VEG_CTL for use in updating soil carbon pools. RESP_S_TOT calculated
! here is passed to BL_TRMIX as the fraction which is respired as CO2
! to the atmosphere. RESP_S_TOT, and RESP_S are also passed out for
! storage in diagnostics 3293, and 3467-470.

!-----------------------------------------------------------------------
IF ( soil_bgc_model == soil_model_rothc ) THEN
  DO i = 1,land_pts
    !soil tiling is not compatible with triffid. OK to hard-code soil tile
    !index to 1 here by setting m = 1
    m = 1
    resp_s_tot_soilt(i,m) = 0.0
    DO nn = 1,dim_cslayer
      work_clay = EXP(-0.0786 * 100.0 * clay_soilt(i,m,nn))
      resp_frac(i,nn) = (3.0895+2.672 * work_clay) /                          &
                        (4.0895+2.672 * work_clay)
      resp_s_soilt(i,m,nn,1)  = resp_s_soilt(i,m,nn,1) * resp_frac(i,nn)
      resp_s_soilt(i,m,nn,2)  = resp_s_soilt(i,m,nn,2) * resp_frac(i,nn)
      resp_s_soilt(i,m,nn,3)  = resp_s_soilt(i,m,nn,3) * resp_frac(i,nn)
      resp_s_soilt(i,m,nn,4)  = resp_s_soilt(i,m,nn,4) * resp_frac(i,nn)
      resp_s_tot_soilt(i,m)   = resp_s_tot_soilt(i,m)                         &
                                + resp_s_soilt(i,m,nn,1)                      &
                                + resp_s_soilt(i,m,nn,2)                      &
                                + resp_s_soilt(i,m,nn,3)                      &
                                + resp_s_soilt(i,m,nn,4)
    END DO  !  layers
  END DO  !  points
END IF


!-----------------------------------------------------------------------
! Reset surft_pts and surft_index and set tile fractions to 1 if aggregate
! tiles are used (L_AGGREGATE=.T.).
! Otherwise, set tile fractions to surface type fractions.
!-----------------------------------------------------------------------
IF (l_aggregate) THEN
  surft_pts(1) = land_pts
  DO l = 1,land_pts
    tile_frac(l,1) = 1.0
    surft_index(l,1) = l
  END DO
ELSE
  DO n = 1,ntype
    DO l = 1, land_pts
      tile_frac(l,n) = frac(l,n)
    END DO
  END DO
END IF

IF (land_pts >  0) THEN    ! Omit if no land points

  !-----------------------------------------------------------------------
  ! Calculate the thermal conductivity of the top soil layer.
  !-----------------------------------------------------------------------
  DO m = 1, nsoilt
    CALL heat_con (land_pts,hcon_soilt,sthu_soilt(:,m,1),                     &
                   sthf_soilt(:,m,1),smvcst_soilt(:,m,1),hcons_soilt(:,m))
  END DO

  ! Thermal conductvity of top snow layer if nsmax > 0
  IF (nsmax > 0) THEN
    DO n = 1,nsurft
      CALL snowtherm(land_pts,surft_pts(n),nsnow_surft(:,n),                  &
                     surft_index(:,n),ds_surft(:,n,:),sice_surft(:,n,:),      &
                     sliq_surft(:,n,:),csnow,ksnow)
      DO l = 1,land_pts
        hcons_snow(l,n) = ksnow(l,1)
      END DO
    END DO
  END IF

END IF                     ! End test on land points

!-----------------------------------------------------------------------
! Calculate net radiation on land tiles
!-----------------------------------------------------------------------
!$OMP PARALLEL                                                                &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(l,k,j,i,n)                                                      &
!$OMP SHARED(surft_pts,surft_index,land_index,pdims,radnet_surft,sw_surft,    &
!$OMP        emis_surft,sky,lw_down,tstar_surft,nsurft,l_skyview)
IF (l_skyview) THEN
  DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
    DO k = 1,surft_pts(n)
      l = surft_index(k,n)
      j=(land_index(l) - 1) / pdims%i_end + 1
      i = land_index(l) - (j-1) * pdims%i_end
      radnet_surft(l,n) = sw_surft(l,n) + emis_surft(l,n) *                   &
        sky(i,j) * ( lw_down(i,j) - sbcon * tstar_surft(l,n)**4 )
    END DO
!$OMP END DO NOWAIT
  END DO
ELSE
  DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
    DO k = 1,surft_pts(n)
      l = surft_index(k,n)
      j=(land_index(l) - 1) / pdims%i_end + 1
      i = land_index(l) - (j-1) * pdims%i_end
      radnet_surft(l,n) = sw_surft(l,n) + emis_surft(l,n) *                   &
                 ( lw_down(i,j) - sbcon * tstar_surft(l,n)**4 )
    END DO
!$OMP END DO NOWAIT
  END DO
END IF
!$OMP END PARALLEL


!-----------------------------------------------------------------------
! 4.  Surface turbulent exchange coefficients and "explicit" fluxes
!     (P243a, routine SF_EXCH).
!     Wind mixing "power" and some values required for other, later,
!     diagnostic calculations, are also evaluated if requested.
!-----------------------------------------------------------------------

IF ( l_rp2 .AND. i_rp_scheme == i_rp2b) THEN
  DO n = 1,npft
    z0h_z0m(n) = z0hm_pft_rp(n)
  END DO
END IF

!$OMP PARALLEL PRIVATE(i,j,l,n) DEFAULT(NONE)                             &
!$OMP SHARED(land_pts,nsurft,land_index, flandg,z_land,t_i_length,        &
!$OMP can_model,cansnowtile,l_snowdep_surf,snow_surft,rho_snow_const,     &
!$OMP snowdep_surft,surf_hgt_surft,l_fix_ctile_orog,l_ctile)    
DO n = 1,nsurft
  IF ( (can_model == 4) .AND. cansnowtile(n) .AND. l_snowdep_surf ) THEN
!$OMP DO SCHEDULE(STATIC)
    DO l = 1,land_pts
      snowdep_surft(l,n) = snow_surft(l,n) / rho_snow_const
    END DO
!$OMP END DO NOWAIT
  END IF
END DO

IF (l_fix_ctile_orog .AND. l_ctile) THEN
!$OMP DO SCHEDULE(STATIC)
  DO l = 1,land_pts
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    IF (flandg(i,j) > 0.0 .AND. flandg(i,j) < 1.0) THEN
      ! calculate height of orography relative to grid-box mean
      ! limit this to 1000m. z_land already limited to be > 0
      ! in ni_bl_ctl
      surf_hgt_surft(l,:) = MIN(z_land(i,j) * (1.0 / flandg(i,j) - 1.0),1000.0)
    END IF
  END DO
!$OMP END DO NOWAIT
END IF
!$OMP END PARALLEL


! Calculate temeprature and specific humidity for elevation bands
CALL elevate(                                                                 &
 land_pts,nsurft,surft_pts,land_index,surft_index,                            &
 tl_1,qw_1,pstar,surf_hgt_surft,l_elev_absolute_height,z_land,                &
 t_elev,q_elev)


!-----------------------------------------------------------------------
!  2.  Calculate QSAT values required later.
!-----------------------------------------------------------------------

!$OMP PARALLEL DO IF(land_pts > 1) DEFAULT(NONE) PRIVATE(i, j, l)             &
!$OMP SHARED(land_index, land_pts, pstar, pstar_land, t_i_length)             &
!$OMP SCHEDULE(STATIC)
DO l = 1,land_pts
  j=(land_index(l) - 1) / t_i_length + 1
  i = land_index(l) - (j-1) * t_i_length
  pstar_land(l) = pstar(i,j)
END DO
!$OMP END PARALLEL DO

IF (l_mr_physics) THEN
!$OMP PARALLEL DO IF(nsurft > 1) DEFAULT(NONE) PRIVATE(n)                     &
!$OMP SHARED(nsurft,qstar_surft,tstar_surft,pstar_land,land_pts,              &
!$OMP qs1_elev,t_elev)  SCHEDULE(STATIC)
  DO n = 1,nsurft
    CALL qsat_mix(qstar_surft(:,n),tstar_surft(:,n),pstar_land,land_pts)
    CALL qsat_mix(qs1_elev(:,n),t_elev(:,n),pstar_land,land_pts)
  END DO
!$OMP END PARALLEL DO
ELSE
!$OMP PARALLEL DO IF(nsurft > 1) DEFAULT(NONE) PRIVATE(n)                     &
!$OMP SHARED(nsurft,qstar_surft,tstar_surft,pstar_land,land_pts,              &
!$OMP qs1_elev,t_elev)  SCHEDULE(STATIC)
  DO n = 1,nsurft
    CALL qsat(qstar_surft(:,n),tstar_surft(:,n),pstar_land,land_pts)
    CALL qsat(qs1_elev(:,n),t_elev(:,n),pstar_land,land_pts)
  END DO
!$OMP END PARALLEL DO
END IF

!-----------------------------------------------------------------------
!  3. Calculation of transfer coefficients and surface layer stability
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  3.1 Calculate neutral roughness lengths
!-----------------------------------------------------------------------

! Land tiles
! Z0_SURFT contains the appropriate value for land-ice points, but has to
! be modified for snow-cover on non-land-ice points. Thermal roughness
! lengths are set to be proportional to the tiled roughness length in
! the case of multiple tiles, but if tiled properties have already
! been aggregated, an adjustment for snow cover is required. In this case
! a ratio of 0.1 between the thermal and momentum roughness lengths over
! snow is assumed; to do otherwise would require reaggregation. In the
! case of multiple tiles, the assignment is delayed until the urban
! options have been considered.
!
!$OMP PARALLEL DO IF(nsurft > 1) DEFAULT(NONE) PRIVATE(k, l, n, z0, zeta1)    &
!$OMP          SHARED(nsurft, surft_pts, surft_index, snow_surft,             &
!$OMP                 l_soil_point, z0_surft, snowdep_surft, l_aggregate,     &
!$OMP                 i_aggregate_opt, z0h_surft_bare, z0m_surft,             &
!$OMP                 l_moruses_rough, urban_canyon, urban_roof,              &
!$OMP                 ztm_gb, z0h_z0m, db_surft, z0h_surft,                   &
!$OMP                 wind_profile_factor, z0m_eff_surft, z0h_surft_classic,  &
!$OMP                 z0h_z0m_classic)   SCHEDULE(STATIC)
DO n = 1,nsurft
  DO k = 1,surft_pts(n)
    l = surft_index(k,n)
    IF ( snow_surft(l,n) > 0.0 .AND. l_soil_point(l) ) THEN
      z0 = z0_surft(l,n) - 0.1 * snowdep_surft(l,n)
      zeta1 = MIN( 5.0e-4 , z0_surft(l,n)  )
      z0m_surft(l,n) = MAX( zeta1 , z0 )
      !     Set z0h_surft explicitly if this option is selected,
      !     otherwise, it will be set for the first tile below.
      IF (l_aggregate .AND. i_aggregate_opt == 1) THEN
        z0 = z0h_surft_bare(l,n) - 0.1 * 0.1 * snowdep_surft(l,n)
        zeta1 = MIN( 5.0e-5 , z0h_surft_bare(l,n)  )
        z0h_surft(l,n) = MAX( zeta1 , z0 )
      END IF
    ELSE
      z0m_surft(l,n) = z0_surft(l,n)
      !     Set z0h_surft explicitly if this option is selected,
      !     otherwise, it will be set for the first tile below.
      IF (l_aggregate .AND. i_aggregate_opt == 1)                             &
        z0h_surft(l,n) = z0h_surft_bare(l,n)
    END IF

    ! MORUSES: z0m_surft (z0_surft previously set in sparm) reset with ztm_gb
    ! (calculated in init_urban) to remove effects of snow. MORUSES does not
    ! yet contain a parametrisation for snow although snow should not affect
    ! the behaviour of bluff bodies rather it affects the  material roughness
    ! of the road & roof
    ! Note: The test on smvcst (not replaced with l_soil_point) is used to alter
    !       just the roof tile and not the
    !       land-ice tiles as the tile is shared at the moment. This will not
    !       be required when flexible tiles are introduced.

    IF ( .NOT. l_aggregate                                                    &
       .AND. l_moruses_rough                                                  &
       .AND. l_soil_point(l)                                                  &
       .AND. ( n == urban_canyon .OR. n == urban_roof ) ) THEN
      ! Snow could be added here to the material roughness length for
      ! momentum before passing to urbanz0. Only the road & roof should be
      ! affected as the walls will be essentially snow-free
      z0m_surft(l,n) = ztm_gb(l)
    END IF

    !
    !   Set the thermal roughness length if aggregation is not being
    !   carried out, or if the original scheme is being used.
    !   It must be done here for consistency with the urban options.
    IF ( ( .NOT. l_aggregate) .OR.                                            &
         (l_aggregate .AND. i_aggregate_opt == 0) )                           &
      z0h_surft(l,n) = z0h_z0m(n) * z0m_surft(l,n)
    db_surft(l,n) = 0.0
    wind_profile_factor(l,n) = 1.0
    z0m_eff_surft(l,n) = z0m_surft(l,n)

    ! Also set additional roughness length for use in CLASSIC aerosol
    ! deposition scheme
    z0h_surft_classic(l,n) = z0h_z0m_classic(n) * z0m_surft(l,n)
  END DO
END DO
!$OMP END PARALLEL DO

! MORUSES: urbanz0 calculates r.l. for heat and updates z0h_surft, which is
! previously set above. Two calls are required; one for urban_canyon and one
! for urban_roof. It is done this way to remove the need for do loop over tile
! type and to avoid over-writing any land-ice tiles. If l_aggregate roughness
! lengths are already calculated: z0m in sparm & z0h above.
IF ( .NOT. l_aggregate ) THEN
  IF ( l_moruses_rough ) THEN
    n = urban_canyon
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,l,j,i)             &
!$OMP SHARED(surft_pts,n,surft_index,land_index,t_i_length,z1_uv,z1_tq,hgt_gb,&
!$OMP hwr_gb,disp_gb,ztm_gb,z0h_surft,urban_roof,                             &
!$OMP z0h_surft_classic)
    DO k = 1,surft_pts(n)
      l = surft_index(k,n)
      j = ( land_index(l) - 1 ) / t_i_length + 1
      i = land_index(l) - ( j - 1 ) * t_i_length
      CALL urbanz0(                                                           &
         n, z1_uv(i,j), z1_tq(i,j), hgt_gb(l), hwr_gb(l), disp_gb(l),         &
         z0m_mat, ztm_gb(l),                                                  &
         z0h_surft(l,n) )
      CALL urbanz0(                                                           &
         urban_roof, z1_uv(i,j), z1_tq(i,j), hgt_gb(l), hwr_gb(l), disp_gb(l),&
         z0m_mat, ztm_gb(l),                                                  &
         z0h_surft(l,urban_roof) )
      ! Make CLASSIC aerosol roughness length for urban tiles consistent
      ! with those for heat and momentum
      z0h_surft_classic(l,n) = z0h_surft(l,n)
      z0h_surft_classic(l,urban_roof) = z0h_surft(l,urban_roof)
    END DO
!$OMP END PARALLEL DO

  END IF

END IF


! Calculate roughness length affected by roughness sublayer in neutral
! conditions.
zdt_surft(:,:) = 0.0
DO n = 1,nsurft
  IF (l_vegdrag_surft(n)) THEN
    CALL can_drag_z0(                                                         &
      land_pts, surft_pts(n), surft_index(:,n),                               &
      array_zero, canht_pft(:,n), lai_pft(:,n),                               &
      z0m_surft(:,n), z0h_surft(:,n), zdt_surft(:,n))
  END IF
END DO

! Calculate orographic effective parameter for neutral conditions
! if using orographic roughness scheme
IF (formdrag ==  effective_z0) THEN
  DO n = 1,nsurft
    CALL sf_orog (                                                            &
     land_pts,surft_pts(n),land_index,surft_index(:,n),                       &
     fd_stab_dep,orog_drag_param,                                             &
     ho2r2_orog,rib_surft(:,n),sil_orog_land,z0m_surft(:,n),z1_uv,            &
     wind_profile_factor(:,n),z0m_eff_surft(:,n)                              &
     )
  END DO
END IF

!-----------------------------------------------------------------------
! Calculate RESFT with neutral CH and EPDT = 0 for use in calculation
! of Richardson number. RESFT=1 for snow and land-ice.
!-----------------------------------------------------------------------

DO n = 1,nsurft
  IF (l_vegdrag_surft(n)) THEN
    CALL can_drag_phi_m_h(                                                    &
      land_pts, surft_pts(n), surft_index(:,n), land_index,                   &
      array_zero, z1_uv, z1_tq, canht_pft(:,n), lai_pft(:,n),                 &
      z0m_surft(:,n), z0h_surft(:,n),                                         &
      phi_m, phi_h)

!$OMP PARALLEL DO IF(surft_pts(n)>1) DEFAULT(NONE)                            &
!$OMP PRIVATE(i, j, k, l)                                                     &
!$OMP SHARED(surft_pts, surft_index, t_i_length, land_index,                  &
!$OMP        phi_m, phi_h, chn, wind_profile_factor, dq, qw_1,                &
!$OMP        qstar_surft, epdt, n)    SCHEDULE(STATIC)
    DO k = 1,surft_pts(n)
      l = surft_index(k,n)
      j=(land_index(l) - 1) / t_i_length + 1
      i = land_index(l) - (j-1) * t_i_length
      chn(l,n) = vkman**2 / (phi_m(l) * phi_h(l)) * wind_profile_factor(l,n)
      dq(l) = qw_1(i,j) - qstar_surft(l,n)
      epdt(l) = 0.0
    END DO
!$OMP END PARALLEL DO

  ELSE ! l_vegdrag_surft = F

!$OMP PARALLEL DO IF(surft_pts(n)>1) DEFAULT(NONE)                            &
!$OMP PRIVATE(i, j, k, l, zetah, zetam)                                       &
!$OMP SHARED(surft_pts, surft_index, t_i_length, land_index, z1_uv,           &
!$OMP        z0m_surft, z1_tq, z0h_surft, chn, wind_profile_factor, dq, qw_1, &
!$OMP        qstar_surft, epdt, n)    SCHEDULE(STATIC)
    DO k = 1,surft_pts(n)
      l = surft_index(k,n)
      j=(land_index(l) - 1) / t_i_length + 1
      i = land_index(l) - (j-1) * t_i_length
      zetam = LOG ( (z1_uv(i,j) + z0m_surft(l,n)) / z0m_surft(l,n) )
      zetah = LOG ( (z1_tq(i,j) + z0m_surft(l,n)) / z0h_surft(l,n) )
      chn(l,n) = (vkman / zetah) * (vkman / zetam) *                          &
        wind_profile_factor(l,n)
      dq(l) = qw_1(i,j) - qstar_surft(l,n)
      epdt(l) = 0.0
    END DO
!$OMP END PARALLEL DO
  END IF

  ! We should only attempt to access sf_diag%resfs_stom(:,n) if it has
  ! been fully allocated.
  IF (sf_diag%l_et_stom .OR. sf_diag%l_et_stom_surft) THEN
    n_diag = n
  ELSE
    n_diag = 1
  END IF

  CALL sf_resist (                                                            &
   land_pts,surft_pts(n),land_index,surft_index(:,n),                         &
   canopy(:,n),catch(:,n),chn(:,n),dq,epdt,flake(:,n),gc_surft(:,n),          &
   gc_stom_surft(:,n),snowdep_surft(:,n),vshr_land,fraca(:,n),resfs(:,n),     &
   resft(:,n),gc_irr_surft(:,n),resfs_irr_surft(:,n),                         &
   sf_diag%resfs_stom(:,n_diag),sf_diag%l_et_stom, sf_diag%l_et_stom_surft)

END DO

! RESFT < 1 for snow on canopy if canopy snow model used
! N.B. chn is calculated in the preceding loop over tiles and
! used in the next loop over npft. This works only if the
! first npft tiles are the vegetated ones.
IF ( .NOT. l_aggregate .AND. can_model == 4) THEN
  DO n = 1,npft
    IF ( cansnowtile(n) ) THEN
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(i, j, k, l)      &
!$OMP          SHARED(surft_pts, surft_index, land_index, t_i_length,         &
!$OMP                 snow_surft, gc_surft, catch_snow, tstar_surft,          &
!$OMP                 vshr_land, fraca, resfs, chn, resft, l_irrig_dmd,       &
!$OMP                 resfs_irr_surft, n)   SCHEDULE(STATIC)
      DO k = 1,surft_pts(n)
        l = surft_index(k,n)
        j=(land_index(l) - 1) / t_i_length + 1
        i = land_index(l) - (j-1) * t_i_length
        IF (snow_surft(l,n) >  0.0) THEN
          gc_surft(l,n) = 0.06 * snow_surft(l,n)**0.6 * catch_snow(l,n)**0.4  &
                       * 2.06e-5 * (tm / tstar_surft(l,n))**1.75              &
                       * (1.79+3 * SQRT(vshr_land(i,j)))                      &
                       / (2 * rho_ice * 5.0e-4**2)
          fraca(l,n) = 0.0
          resfs(l,n) = gc_surft(l,n) /                                        &
                  (gc_surft(l,n) + chn(l,n) * vshr_land(i,j))
          resft(l,n) = resfs(l,n)
          IF (l_irrig_dmd) THEN
            resfs_irr_surft(l,n) = resfs(l,n)
            ! print*,'changing resfs cos of canopy snow?'
          END IF
        END IF
      END DO
!$OMP END PARALLEL DO
    END IF
  END DO
END IF

!-----------------------------------------------------------------------
!  3.2 Calculate bulk Richardson number for the lowest model level.
!-----------------------------------------------------------------------

! Land tiles
DO n = 1,nsurft
  CALL sf_rib (                                                               &
   land_pts,surft_pts(n),land_index,surft_index(:,n),                         &
   bq_1,bt_1,qstar_surft(:,n),q_elev(:,n),resft(:,n),t_elev(:,n),             &
   tstar_surft(:,n),vshr_land,z0h_surft(:,n),z0m_surft(:,n),zdt_surft(:,n),   &
   z1_tq,z1_uv,l_vegdrag_surft(n),                                            &
   rib_surft(:,n),db_surft(:,n)                                               &
   )
END DO

!-----------------------------------------------------------------------
!  3.3 Calculate stability corrected effective roughness length.
!  Stability correction only applies to land points.
!-----------------------------------------------------------------------

IF (formdrag ==  effective_z0) THEN
  DO n = 1,nsurft
    CALL sf_orog (                                                            &
     land_pts,surft_pts(n),land_index,surft_index(:,n),                       &
     fd_stab_dep,orog_drag_param,                                             &
     ho2r2_orog,rib_surft(:,n),sil_orog_land,z0m_surft(:,n),z1_uv,            &
     wind_profile_factor(:,n),z0m_eff_surft(:,n)                              &
     )
  END DO
END IF

!-----------------------------------------------------------------------
!  3.4 Calculate CD, CH via routine FCDCH.
!-----------------------------------------------------------------------

! Land tiles
DO n = 1,nsurft
  IF (l_vegdrag_surft(n)) THEN
    n_veg = n
    z0m_eff_surft(:,n) = z0m_surft(:,n)
  ELSE
    n_veg = 1
  END IF
  CALL fcdch (                                                                &
   cor_mo_iter,land_pts,surft_pts(n),                                         &
   surft_index(:,n),land_index,                                               &
   db_surft(:,n),vshr_land,                                                   &
   z0m_eff_surft(:,n),z0h_surft(:,n),zdt_surft(:,n),zh,                       &
   z1_uv,z1_uv_top,z1_tq,z1_tq_top,wind_profile_factor(:,n),                  &
   ddmfx,ip_ss_solid,charnock,                                                &
   charnock_w,                                                                &
   l_vegdrag_surft(n),canht_pft(:,n_veg),lai_pft(:,n_veg),                    &
   cd_surft(:,n),ch_surft(:,n),cd_std(:,n),                                   &
   v_s_surft(:,n),v_s_std(:,n),recip_l_mo_surft(:,n),                         &
   u_s_iter_surft(:,n)                                                        &
   )
END DO

! As roughness length have been changed by vegetation drag effect, effective
! roughness length should be updated.
DO n = 1,nsurft
  IF (l_vegdrag_surft(n)) THEN
    z0m_surft(:,n) = z0m_eff_surft(:,n)
    IF (formdrag == effective_z0) THEN
      CALL sf_orog (                                                          &
       land_pts,surft_pts(n),land_index,surft_index(:,n),                     &
       fd_stab_dep,orog_drag_param,                                           &
       ho2r2_orog,rib_surft(:,n),sil_orog_land,z0m_surft(:,n),z1_uv,          &
       wind_profile_factor(:,n),z0m_eff_surft(:,n)                            &
       )
    END IF
  END IF
END DO

!$OMP PARALLEL IF(nsurft > 1) DEFAULT(NONE) PRIVATE(k, l, n)                   &
!$OMP          SHARED(cor_mo_iter, nsurft, surft_pts, surft_index,             &
!$OMP                 u_s_iter_surft, u_s_std_surft, v_s_std)
IF ( cor_mo_iter >= use_correct_ustar ) THEN
  !       Use correct "standard" ustar
!$OMP DO SCHEDULE(STATIC)
  DO n = 1,nsurft
    DO k = 1,surft_pts(n)
      l = surft_index(k,n)
      u_s_std_surft(l,n) = v_s_std(l,n)
    END DO
  END DO
!$OMP END DO
ELSE
  !       Use ustar from mid-iteration
!$OMP DO SCHEDULE(STATIC)
  DO n = 1,nsurft
    DO k = 1,surft_pts(n)
      l = surft_index(k,n)
      u_s_std_surft(l,n) = u_s_iter_surft(l,n)
    END DO
  END DO
!$OMP END DO
END IF
!$OMP END PARALLEL

!-----------------------------------------------------------------------
!  3.5 Recalculate friction velocity for the dust scheme using the
!      bare soil roughness length if using only 1 aggregated tile
!-----------------------------------------------------------------------

IF ((l_dust .OR. l_dust_diag) .AND. l_aggregate) THEN

  n = 1
  ! Calculate z0m and z0h for bare soil
  IF (l_vary_z0m_soil) THEN
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(k, l)            &
!$OMP             SHARED(n, soil, surft_pts, surft_index, z0m_soil_in,        &
!$OMP                   z0h_soil,z0h_z0m, z0m_soil, wind_profile_fac_soil)    &
!$OMP             SCHEDULE(STATIC)
    DO k = 1, surft_pts(n)
      l = surft_index(k,n)
      z0m_soil(l,n) = z0m_soil_in(l)
      z0h_soil(l,n) = z0h_z0m(soil) * z0m_soil(l,n)
      !     Set wind profile factor to 1 as not using orog term in z0m
      wind_profile_fac_soil(l,n) = 1.0
    END DO
!$OMP END PARALLEL DO
  ELSE
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(k, l)            &
!$OMP             SHARED(n, soil, surft_pts, surft_index, z0_soil, z0h_soil,  &
!$OMP                   z0h_z0m, z0m_soil, wind_profile_fac_soil)             &
!$OMP             SCHEDULE(STATIC)
    DO k = 1, surft_pts(n)
      l = surft_index(k,n)
      z0m_soil(l,n) = z0_soil
      z0h_soil(l,n) = z0h_z0m(soil) * z0m_soil(l,n)
      !     Set wind profile factor to 1 as not using orog term in z0m
      wind_profile_fac_soil(l,n) = 1.0
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Call fcdch again to calculate dust friction velocity on bare soil.
  ! The canopy drag scheme is not available on the aggregated tile
  ! and is disabled.
  l_vegdrag_active_here = .FALSE.
  CALL fcdch (                                                                &
   cor_mo_iter,land_pts,surft_pts(n),                                         &
   surft_index(:,n),land_index,                                               &
   db_surft(:,n),vshr_land,                                                   &
   z0m_soil(:,n),z0h_soil(:,n),zdt_dummy,zh,                                  &
   z1_uv,z1_uv_top,z1_tq,z1_tq_top,wind_profile_fac_soil(:,n),                &
   ddmfx,ip_ss_solid,charnock,                                                &
   charnock_w,                                                                &
   l_vegdrag_active_here,array_zero,array_zero,                               &
  !  Following tiled outputs (except v_s_std_soil and u_s_iter_soil)
  !  are dummy variables not needed from this call
     cd_surft_soil(:,n),ch_surft_soil(:,n),cd_std_soil(:,n),                  &
     v_s_surft_soil(:,n),v_s_std_soil(:,n),recip_l_mo_surft_soil(:,n),        &
     u_s_iter_soil(:,n)                                                       &
     )

!$OMP PARALLEL IF(nsurft > 1) DEFAULT(NONE) PRIVATE(k, l, n)                   &
!$OMP          SHARED(cor_mo_iter, nsurft, surft_index, surft_pts,             &
!$OMP                 u_s_iter_soil, u_s_std_surft, v_s_std_soil)
  IF ( cor_mo_iter >= use_correct_ustar ) THEN
    !       Use correct "standard" ustar
!$OMP DO SCHEDULE(STATIC)
    DO n = 1,nsurft
      DO k = 1,surft_pts(n)
        l = surft_index(k,n)
        u_s_std_surft(l,n) = v_s_std_soil(l,n)
      END DO
    END DO
!$OMP END DO
  ELSE
    !       Use ustar from mid-iteration
!$OMP DO SCHEDULE(STATIC)
    DO n = 1,nsurft
      DO k = 1,surft_pts(n)
        l = surft_index(k,n)
        u_s_std_surft(l,n) = u_s_iter_soil(l,n)
      END DO
    END DO
!$OMP END DO
  END IF
!$OMP END PARALLEL

END IF

!-----------------------------------------------------------------------
!  3.6 Recalculate cd, ch etc. using z0h=z0h_classic. The parameters
!       calculated using this additional roughness length are for
!       CLASSIC aerosol deposition only.
!-----------------------------------------------------------------------

IF (l_aero_classic) THEN
  ! The canopy drag scheme is not supported for this aerosol scheme and
  ! is turned off.
  l_vegdrag_active_here = .FALSE.
  ! Land tiles
  DO n = 1,nsurft
    CALL fcdch (                                                              &
    !  Input variables identical to main call except using different z0h
         cor_mo_iter,land_pts,surft_pts(n),                                   &
         surft_index(:,n),land_index,                                         &
         db_surft(:,n),vshr_land,                                             &
         z0m_eff_surft(:,n),z0h_surft_classic(:,n),zdt_dummy,zh,              &
         z1_uv,z1_uv_top,z1_tq,z1_tq_top,                                     &
         wind_profile_factor(:,n),                                            &
         ddmfx,ip_ss_solid,charnock,                                          &
         charnock_w,                                                          &
         l_vegdrag_active_here,array_zero,array_zero,                         &
    !  Following tiled outputs (except cd_std_classic and ch_surft_classic)
    !  are dummy variables not needed from this call
         cd_surft_classic(:,n),ch_surft_classic(:,n),                         &
         cd_std_classic(:,n),v_s_surft_classic(:,n),                          &
         v_s_std_classic(:,n),recip_l_mo_surft_classic(:,n),                  &
         u_s_iter_classic(:,n)                                                &
         )
  END DO
END IF

!-----------------------------------------------------------------------
! 4.0 If requested, improve accuracy of air density, rhostar
!     On input rhostar = pstar/R*Tstar
!-----------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j) SHARED(tdims,rhostar_mom,rhostar)
!$OMP DO SCHEDULE(STATIC)
DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    ! original approximation for surface air density
    rhostar_mom(i,j) = rhostar(i,j)
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

IF (l_accurate_rho) THEN
  ! More accurate expressions for surface air density.
  ! Use bottom level vapour as a better approximation over land 
  ! than qsat(tstar)
  CALL calc_air_dens(l_mr_physics,qw_1,rhostar,rhostar_mom)

END IF ! l_accurate_rho

!-----------------------------------------------------------------------
!  4.1 Recalculate RESFT using "true" CH and EPDT for land tiles
!-----------------------------------------------------------------------

DO n = 1,nsurft
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(i, j, k, l) &
!$OMP SHARED(ch_surft, land_index, qstar_surft, qw_1, dq, epdt,n,        &
!$OMP rhostar, surft_pts, surft_index, timestep, t_i_length,             &
!$OMP vshr_land) SCHEDULE(STATIC)
  DO k = 1,surft_pts(n)
    l = surft_index(k,n)
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    dq(l) = qw_1(i,j) - qstar_surft(l,n)
    epdt(l) = - rhostar(i,j) * ch_surft(l,n) * vshr_land(i,j)                 &
      *dq(l) * timestep
  END DO
!$OMP END PARALLEL DO

  ! We should only attempt to access sf_diag%resfs_stom(:,n) if it has
  ! been fully allocated.
  IF (sf_diag%l_et_stom .OR. sf_diag%l_et_stom_surft) THEN
    n_diag = n
  ELSE
    n_diag = 1
  END IF

  CALL sf_resist (                                                            &
   land_pts,surft_pts(n),land_index,surft_index(:,n),                         &
   canopy(:,n),catch(:,n),ch_surft(:,n),dq,epdt,flake(:,n),                   &
   gc_surft(:,n),gc_stom_surft(:,n),snowdep_surft(:,n),vshr_land,fraca(:,n),  &
   resfs(:,n),resft(:,n),                                                     &
   gc_irr_surft(:,n),resfs_irr_surft(:,n),sf_diag%resfs_stom(:,n_diag),       &
   sf_diag%l_et_stom,sf_diag%l_et_stom_surft)

END DO

IF ( .NOT. l_aggregate .AND. can_model == 4) THEN
  DO n = 1,npft
    IF ( cansnowtile(n) ) THEN
!$OMP PARALLEL DO IF(surft_pts(n) > 1) DEFAULT(NONE) PRIVATE(i, j, k, l)       &
!$OMP             SHARED(ch_surft, fraca, gc_surft, land_index, l_irrig_dmd, n,&
!$OMP                    snow_surft, resfs, resfs_irr_surft, resft, surft_pts, &
!$OMP                    surft_index, t_i_length, vshr_land) SCHEDULE(STATIC)
      DO k = 1,surft_pts(n)
        l = surft_index(k,n)
        j=(land_index(l) - 1) / t_i_length + 1
        i = land_index(l) - (j-1) * t_i_length
        IF (snow_surft(l,n) >  0.0) THEN
          fraca(l,n) = 0.0
          resfs(l,n) = gc_surft(l,n) /                                        &
                  (gc_surft(l,n) + ch_surft(l,n) * vshr_land(i,j))
          resft(l,n) = resfs(l,n)
          IF (l_irrig_dmd) THEN
            resfs_irr_surft(l,n) = resfs(l,n)
          END IF
        END IF
      END DO
!$OMP END PARALLEL DO
    END IF
  END DO
END IF

!-----------------------------------------------------------------------
! Calculate gridbox-means of transfer coefficients.
!-----------------------------------------------------------------------

! Land tiles
DO n = 1,nsurft
  DO k = 1,surft_pts(n)
    l = surft_index(k,n)
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    cd_land(i,j) = cd_land(i,j) + tile_frac(l,n) * cd_surft(l,n)
    ch_land(i,j) = ch_land(i,j) + tile_frac(l,n) * ch_surft(l,n)
  END DO
END DO

! aerodynamic resistance diagnostic
IF (sf_diag%l_ra) THEN
  DO l = 1,land_pts
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    sf_diag%ra(l) = 1.0 / (ch_land(i,j) * vshr_land(i,j))
  END DO
END IF

!-----------------------------------------------------------------------
!  4.3 Calculate the surface exchange coefficients RHOK(*) and
!       resistances for use in CLASSIC aerosol scheme
!       (Note that CD_STD, CH and VSHR should never = 0)
!     RHOSTAR * CD * VSHR stored for diagnostic output before
!     horizontal interpolation.
!-----------------------------------------------------------------------

! Land tiles
DO n = 1,nsurft
  DO k = 1,surft_pts(n)
    l = surft_index(k,n)
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    rhokm_1_surft(l,n) = rhostar_mom(i,j) * cd_surft(l,n) * vshr_land(i,j)
    !                                                         ! P243.124
    rhokm_land(i,j) = rhokm_land(i,j) +                                       &
           tile_frac(l,n) * rhokm_1_surft(l,n)
    rhokh_surft(l,n) = rhostar_mom(i,j) * ch_surft(l,n) * vshr_land(i,j)
    !                                                         ! P243.125
  END DO
END DO

!-----------------------------------------------------------------------
!  Calculate local and gridbox-average surface fluxes of heat and
!  moisture.
!-----------------------------------------------------------------------

! Adjust ASHTF for sens. heat flux to ground beneath coniferous canopy
rhokh_can(:,:) = 0.0

IF ( .NOT. l_aggregate .AND. can_model == 4) THEN
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(n, k, l, j, i) SHARED(cansnowtile,     &
!$OMP             cd_surft, land_index, npft, rhokh_can, rhostar, surft_pts,   &
!$OMP             surft_index, t_i_length, vshr_land, cp) SCHEDULE(STATIC)
  DO n = 1,npft
    IF ( cansnowtile(n) ) THEN
      DO k = 1,surft_pts(n)
        l = surft_index(k,n)
        j=(land_index(l) - 1) / t_i_length + 1
        i = land_index(l) - (j-1) * t_i_length
        rhokh_can(l,n) = rhostar(i,j) * cp /                                  &
                     (43.0 / (SQRT(cd_surft(l,n)) * vshr_land(i,j)))
      END DO
    END IF
  END DO
!$OMP END PARALLEL DO
END IF

! Initialise scaling_urban to 1.0 so that it only affects urban tiles when
! MORUSES used with no aggregation.
scaling_urban(:,:) = 1.0
IF ( .NOT. l_aggregate .AND. l_moruses_storage ) THEN
  n = urban_canyon
!$OMP PARALLEL DO IF(land_pts > 1) DEFAULT(NONE) PRIVATE(l) SHARED(land_pts,   &
!$OMP             n, scaling_urban, tile_frac, urban_roof) SCHEDULE(STATIC)
  DO l = 1, land_pts
    IF ( tile_frac(l,n) > 0.0 ) THEN
      scaling_urban(l,n) =                                                    &
         ( tile_frac(l,n) + tile_frac(l,urban_roof) ) /                       &
         tile_frac(l,n)
    END IF
  END DO
!$OMP END PARALLEL DO
END IF

! Calculate average layer temperature and conductivity for lakes.
! This is a fudge - a layer with average properties won't
! really behave like a stack of layers with different properties.
!
IF (     l_flake_model                                                        &
    .AND. ( .NOT. l_aggregate)) THEN

  !==============================================================================
  ! *NOTICE REGARDING SOIL TILING**
  !
  !The following section facilitates the use of soil tiling. As implemented,
  !there are two soil tiling options:
  !
  !nsoilt == 1
  !Operate as with a single soil tile, functionally identical to JULES upto
  ! at least vn4.7 (Oct 2016)
  ! This means that a soilt variable being passed 'up' to the surface is
  ! broadcast to the surft variable (with weighting by frac if requred)
  !
  !nsoilt > 1
  !Operate with nsoilt = nsurft, with a direct mapping between them
  ! This means that a soilt variable being passed 'up' to the surface is simply
  ! copied into the surft variable
  !
  ! This will need to be refactored for other tiling approaches. This note
  ! will be replicated elsewhere in the code as required
  !
  !These comments apply until **END NOTICE REGARDING SOIL TILING**
  !==============================================================================


    !Calculate GBM values for soilt tiled variables
  IF (nsoilt == 1) THEN
    !Just 1 soil tile
    m = 1
    hcons_mean_soil(:) = hcons_soilt(:,m)
    tsoil_mean_soil(:) = t_soil_soilt(:,m,1)
  ELSE
    !Surface tiles map directly on to soil tiles
    hcons_mean_soil(:) = 0.0
    tsoil_mean_soil(:) = 0.0
    DO m = 1,nsoilt
      hcons_mean_soil(:) = hcons_mean_soil(:)                                 &
                          + (tile_frac(:,m) * hcons_soilt(:,m))
      tsoil_mean_soil(:) = tsoil_mean_soil(:)                                 &
                          + (tile_frac(:,m) * t_soil_soilt(:,m,1))
    END DO
  END IF


  !==============================================================================
  ! *END NOTICE REGARDING SOIL TILING**
  !==============================================================================

!$OMP PARALLEL DO IF(land_pts > 1) DEFAULT(NONE)                              &
!$OMP PRIVATE(l, errcode, cmessage)                                           &
!$OMP SHARED(land_pts, hcon_lake, dzsoil, hcons_mean_soil, lake_h_snow_gb,    &
!$OMP        lake_h_ice_gb, snow_hcon, lake_depth_gb, nusselt_gb,             &
!$OMP        g_dt_gb, ts1_lake_gb, lake_t_snow_gb, lake_t_mxl_gb,             &
!$OMP        lake_t_ice_gb, lake_h_mxl_gb, lake_t_mean_gb, tsoil_mean_soil)   &
!$OMP SCHEDULE(STATIC)
  DO l = 1, land_pts

    !
    ! first calculate HCON_LAKE
    !---------------------------
    hcon_lake(l) = 0.0

    IF (dzsoil(1) <= 0.0) THEN
      !
      ! catch-all for sillies - just use land value
      !
      errcode = -1
      WRITE(cmessage, '(A,F16.4)') 'Negative value of dzsoil = ', dzsoil(1)
      CALL ereport(routinename, errcode, cmessage)
      hcon_lake(l) = hcons_mean_soil(l)

    ELSE IF (     (dzsoil(1)  > 0.0)                                          &
             .AND. (dzsoil(1) <= lake_h_snow_gb(l)) ) THEN

      hcon_lake(l) = snow_hcon

    ELSE IF (     (dzsoil(1)  >  lake_h_snow_gb(l))                           &
             .AND. (dzsoil(1) <= (lake_h_snow_gb(l)                           &
                             +lake_h_ice_gb( l)))) THEN

      hcon_lake(l) = ( snow_hcon * lake_h_snow_gb(l)                          &
                      +hcice * (dzsoil(1) - lake_h_snow_gb(l)))               &
                    / dzsoil(1)

    ELSE IF (     (dzsoil(1)  > (lake_h_snow_gb(l)                            &
                             +lake_h_ice_gb( l)))                             &
             .AND. (dzsoil(1) <= (lake_h_snow_gb(l)                           &
                             +lake_h_ice_gb( l)                               &
                             +lake_depth_gb( l)))) THEN

      nusselt_gb(l) = g_dt_gb(l) * (dzsoil(1) - lake_h_snow_gb(l)             &
                                     - lake_h_ice_gb(l) )                     &
                   / ( 2.0 * hcwat )
      nusselt_gb(l) = MAX( nusselt_gb(l), 1.0 )
      hcon_lake(l) = ( snow_hcon * lake_h_snow_gb(l)                          &
                      +hcice * lake_h_ice_gb(l)                               &
                      +hcwat * (dzsoil(1) - lake_h_snow_gb(l)                 &
                                       - lake_h_ice_gb( l)))                  &
                             * nusselt_gb(l)                                  &
                    / dzsoil(1)

    ELSE

      errcode = -1
      WRITE(cmessage, '(A,F16.4)') 'Unusual value of dzsoil found when '//    &
                                   'computing hcon_lake. Found depth '  //    &
                                   'of first soil layer comparable to ' //    &
                                   'lake depth: dzsoil= ', dzsoil(1)
      CALL ereport(routinename, errcode, cmessage)
      nusselt_gb(l) = g_dt_gb(l) * lake_depth_gb(l)                           &
                   / ( 2.0 * hcwat )
      nusselt_gb(l) = MAX( nusselt_gb(l), 1.0 )
      hcon_lake(l) = ( snow_hcon * lake_h_snow_gb(l)                          &
                      +hcice * lake_h_ice_gb(l)                               &
                      +hcwat * lake_depth_gb(l)                               &
                             * nusselt_gb(l)                                  &
                      +hcons_mean_soil(l) * (dzsoil(1) - lake_h_snow_gb(l)    &
                                             - lake_h_ice_gb( l)              &
                                             - lake_depth_gb( l)))            &
                    / dzsoil(1)

    END IF
    !
    ! second calculate TS1_LAKE
    !----------------------------
    ts1_lake_gb(l) = 0.0

    IF (dzsoil(1) <= 0.0) THEN
      !
      ! catch-all for sillies - just use land value
      !

      ts1_lake_gb(l) = tsoil_mean_soil(l)

    ELSE IF (     ((dzsoil(1) / 2.0)  > 0.0)                                  &
             .AND. ((dzsoil(1) / 2.0) <= lake_h_snow_gb(l)) ) THEN

      ts1_lake_gb(l) = lake_t_snow_gb(l) +                                    &
                    (lake_t_ice_gb(l) - lake_t_snow_gb(l))                    &
                   *(dzsoil(1) / 2.0)                                         &
                   / lake_h_snow_gb(l)

    ELSE IF (     ((dzsoil(1) / 2.0)  >  lake_h_snow_gb(l))                   &
             .AND. ((dzsoil(1) / 2.0) <= (lake_h_snow_gb(l)                   &
                                   +lake_h_ice_gb( l)))) THEN

      ts1_lake_gb(l) = lake_t_ice_gb(l) +                                     &
                    (tm - lake_t_ice_gb(l))                                   &
                   *((dzsoil(1) / 2.0) - lake_h_snow_gb(l))                   &
                   / lake_h_ice_gb(l)

    ELSE IF (     ((dzsoil(1) / 2.0)  > (lake_h_snow_gb(l)                    &
                                   +lake_h_ice_gb( l)))                       &
             .AND. ((dzsoil(1) / 2.0) <= (lake_h_snow_gb(l)                   &
                                   +lake_h_ice_gb( l)                         &
                                   +lake_h_mxl_gb( l)))) THEN

      ts1_lake_gb(l) = lake_t_mxl_gb(l)

    ELSE IF (     ((dzsoil(1) / 2.0)  > (lake_h_snow_gb(l)                    &
                                   +lake_h_snow_gb(l)                         &
                                   +lake_h_mxl_gb( l)))                       &
             .AND. ((dzsoil(1) / 2.0) <= (lake_h_snow_gb(l)                   &
                                   +lake_h_ice_gb( l)                         &
                                   +lake_depth_gb( l)))) THEN

      ts1_lake_gb(l) = lake_t_mean_gb(l)

    ELSE
      errcode = -1
      WRITE(cmessage, '(A,F16.4)') 'Unusual value of dzsoil found when '//    &
                                   'computing ts1_lake_gb. Found depth '//    &
                                   'of first soil layer comparable to ' //    &
                                   'lake depth: dzsoil= ', dzsoil(1)
      CALL ereport(routinename, errcode, cmessage)
      ts1_lake_gb(l) = tsoil_mean_soil(l)

    END IF
  END DO
!$OMP END PARALLEL DO
END IF


lw_down_elevcorr_surft(:,:) = 0.0

IF (l_elev_lw_down) THEN
  ! for tiles at different elevations, adjust downwelling longwave
  ! according to the anomaly in surface temperature that has been
  ! calculated with LW ~ T^4
  DO l = 1,land_pts
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length

    IF (lw_down(i,j) > 0.0) THEN

      lw_down_surftsum = 0.0
      lw_down_surftabs = 0.0

      DO n = 1,nsurft
        t_rad = 0.0
        IF (t_elev(l,n) > 0.0 ) THEN
          t_rad = (lw_down(i,j) / sbcon)**(1.0 / 4.0)
          t_rad = t_rad + t_elev(l,n) - tl_1(i,j)
          lw_down_elevcorr_surft(l,n) =  sbcon * (t_rad**4) - lw_down(i,j)
        END IF
        lw_down_surftsum = lw_down_surftsum +                                 &
                           lw_down_elevcorr_surft(l,n) * tile_frac(l,n)
        lw_down_surftabs = lw_down_surftabs +                                 &
                           ABS(lw_down_elevcorr_surft(l,n)) * tile_frac(l,n)
      END DO

      IF (lw_down_surftabs > EPSILON(0.0)) THEN
        ! correct each adjustment to preserve the gridbox mean.
        ! size of correction in proportion to the size of adjustment
        ! so that unadjusted tiles remain unaffected
        DO n = 1,nsurft
          lw_down_elevcorr_surft(l,n) = lw_down_elevcorr_surft(l,n)           &
                                        - lw_down_surftsum                    &
                                          * ABS(lw_down_elevcorr_surft(l,n))  &
                                            / lw_down_surftabs
        END DO
      END IF

      DO n = 1,nsurft
        IF (l_skyview) THEN
          radnet_surft(l,n) = radnet_surft(l,n)                               &
                              + sky(i,j) * emis_surft(l,n)                    &
                                * lw_down_elevcorr_surft(l,n)
        ELSE
          radnet_surft(l,n) = radnet_surft(l,n)                               &
                              + emis_surft(l,n) * lw_down_elevcorr_surft(l,n)
        END IF
      END DO

    END IF

  END DO
END IF

!==============================================================================
! *NOTICE REGARDING SOIL TILING**
!
!The following section facilitates the use of soil tiling. As implemented,
!there are two soil tiling options:
!
!nsoilt == 1
!Operate as with a single soil tile, functionally identical to JULES upto
! at least vn4.7 (Oct 2016)
! This means that a soilt variable being passed 'up' to the surface is
! broadcast to the surft variable (with weighting by frac if requred)
!
!nsoilt > 1
!Operate with nsoilt = nsurft, with a direct mapping between them
! This means that a soilt variable being passed 'up' to the surface is simply
! copied into the surft variable
!
! This will need to be refactored for other tiling approaches. This note
! will be replicated elsewhere in the code as required
!
!These comments apply until **END NOTICE REGARDING SOIL TILING**
!==============================================================================

DO n = 1,nsurft
  lh0 = lc

  !Set the current soil tile (see notice above)
  IF (nsoilt == 1) THEN
    !There is only 1 soil tile
    m = 1
  ELSE ! nsoilt == nsurft
    !Soil tiles map directly on to surface tiles
    m = n
  END IF !nsoilt

!$OMP PARALLEL DO IF(land_pts > 1) DEFAULT(NONE) PRIVATE(i, j, l, t_rad)      &
!$OMP SHARED(land_pts, land_index, t_i_length, tsurf, t_soil_soilt, t_elev,   &
!$OMP        tl_1,                                                            &
!$OMP        dzsurf, dzsoil, hcons_surf, hcons_soilt, canhc_surf, canhc_surft,&
!$OMP        nsmax, nsnow_surft, tsnow_surft, ds_surft, hcons_snow,           &
!$OMP        cansnowtile, snowdepth_surft, snow_hcon,                         &
!$OMP        l_snow_nocan_hc, l_flake_model, l_aggregate, lake,               &
!$OMP        n, m, ts1_lake_gb, hcon_lake, l_elev_land_ice, l_lice_point,     &
!$OMP        tsurf_elev_surft, dzsoil_elev, l_lice_surft, hcondeep)           &
!$OMP SCHEDULE(STATIC)
  DO l = 1,land_pts
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    IF (l_elev_land_ice .AND. l_lice_point(l)) THEN
      tsurf(l) = tsurf_elev_surft(l,n)
      dzsurf(l) = dzsoil_elev
      canhc_surf(l) = 0.0
      IF (l_lice_surft(n)) THEN
        hcons_surf(l) = snow_hcon
      ELSE
        hcons_surf(l) = hcondeep
      END IF
    ELSE
      tsurf(l) = t_soil_soilt(l,m,1) + t_elev(l,n) - tl_1(i,j)
      dzsurf(l) = dzsoil(1)
      hcons_surf(l) = hcons_soilt(l,m)
      canhc_surf(l) = canhc_surft(l,n)
    END IF
    IF ((nsmax > 0) .AND. (nsnow_surft(l,n) > 0)) THEN
      tsurf(l) = tsnow_surft(l,n,1)
      ! change the effective surface layer thickness for snow
      dzsurf(l) = ds_surft(l,n,1)
      hcons_surf(l) = hcons_snow(l,n)
      IF ( ( .NOT. cansnowtile(n)) .AND. l_snow_nocan_hc .AND.                &
           (nsmax > 0) .AND. (nsnow_surft(l,n) > 0) ) canhc_surf(l) = 0.0
    END IF

    IF (     (l_flake_model   )                                               &
        .AND. ( .NOT. l_aggregate)                                            &
        .AND. (n == lake       )                                              &
        ! single-layer snow scheme OR negligible snow
        .AND. ( (nsmax == 0) .OR. (snowdepth_surft(l,n) < EPSILON(1.0) ) )    &
       ) THEN
      tsurf(l) = ts1_lake_gb(l)
      hcons_surf(l) = hcon_lake(l)
      dzsurf(l) = dzsoil(1)
    END IF

  END DO
!$OMP END PARALLEL DO

  !==============================================================================
  ! *END NOTICE REGARDING SOIL TILING**
  !==============================================================================

  sea_point = 0.0
  CALL sf_flux (                                                              &
   land_pts,surft_pts(n),flandg,                                              &
   land_index,surft_index(:,n),                                               &
   nsnow_surft(:,n),n,canhc_surf(:),dzsurf,hcons_surf,                        &
   qs1_elev(:,n),qstar_surft(:,n),q_elev(:,n),                                &
   radnet_surft(:,n),resft(:,n),rhokh_surft(:,n),l_soil_point,                &
   snowdepth_surft(:,n),tile_frac(:,n),timestep,t_elev(:,n),tsurf,            &
   tstar_surft(:,n),vfrac_surft(:,n),rhokh_can(:,n),z0h_surft(:,n),           &
   z0m_eff_surft(:,n),zdt_surft(:,n),z1_tq,lh0,emis_surft(:,n),emis_soil,     &
   1.0,anthrop_heat_surft(:,n),scaling_urban(:,n),l_vegdrag_surft(n),         &
   fqw_1,ftl_1,                                                               &
   alpha1(:,n),ashtf_prime_surft(:,n),fqw_surft(:,n),                         &
   epot_surft(:,n),ftl_surft(:,n),                                            &
   rhokm_1_surft(:,n),vshr_land,tau_surft(:,n),                               &
   dtstar_surft(:,n),sea_point, sf_diag                                       &
 )
END DO


! Set surface stress on tiles diagnostic
IF (sf_diag%l_tau_surft) THEN
  DO n = 1,nsurft
    DO l = 1,land_pts
      sf_diag%tau_surft(l,n) = tau_surft(l,n)
    END DO
  END DO
END IF

! tau_surft no longer required
DEALLOCATE(tau_surft)


!-----------------------------------------------------------------------
!  4.4   Calculate the standard deviations of layer 1 turbulent
!        fluctuations of temperature and humidity using approximate
!        formulae from first order closure.
!-----------------------------------------------------------------------

! Land tiles
DO n = 1,nsurft
  CALL stdev1 (                                                               &
   land_pts,surft_pts(n),land_index,surft_index(:,n),flandg,                  &
   bq_1,bt_1,fqw_surft(:,n),ftl_surft(:,n),rhokm_1_surft(:,n),                &
   rhostar,vshr_land,z0m_surft(:,n),z1_tq,tile_frac(:,n),                     &
   q1_sd,t1_sd                                                                &
   )
END DO

!-----------------------------------------------------------------------
! Call SFL_INT to calculate CDR10M and CHR1P5M - interpolation coeffs
! used to calculate screen temperature, humidity and 10m winds.
!-----------------------------------------------------------------------

! Set flag if cdr10m is to be calculated for use with snow unloading.
! If l_fix_wind_snow=F this calculation will depend on the su10 and sv10
! switches in sf_diag.
l_cdr10m_snow = .FALSE.
IF ( l_fix_wind_snow ) THEN
  DO n = 1,nsurft
    IF ( canSnowTile(n) .AND. unload_rate_u(n) /= 0.0 ) l_cdr10m_snow = .TRUE.
  END DO
END IF

IF (sf_diag%suv10m_n) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j)                 &
!$OMP SHARED(pdims,sf_diag)
  DO j = pdims%j_start,pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      sf_diag%cdr10m_n(i,j) = 0.0
      sf_diag%cd10m_n(i,j)  = 0.0
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! Land tiles
IF (sf_diag%su10 .OR. sf_diag%sv10 .OR. sf_diag%sq1p5 .OR.                    &
    sf_diag%st1p5 .OR. sf_diag%suv10m_n .OR.                                  &
    l_cdr10m_snow .OR.                                                        &
    (IScrnTDiag == IP_ScrnDecpl2) .OR.                                        &
    (IScrnTDiag == IP_ScrnDecpl3) ) THEN
  DO n = 1,nsurft
    CALL sfl_int (                                                            &
     land_pts,surft_pts(n),l_cdr10m_snow,surft_index(:,n),land_index,flandg,  &
     vshr_land,cd_std(:,n),cd_surft(:,n),ch_surft(:,n),                       &
     tile_frac(:,n),                                                          &
     z0m_eff_surft(:,n),z0m_surft(:,n),z0h_surft(:,n),                        &
     recip_l_mo_surft(:,n),                                                   &
     v_s_surft(:,n),v_s_std(:,n),                                             &
     z1_uv,z1_tq,db_surft(:,n),                                               &
     sf_diag,                                                                 &
     cdr10m,sf_diag%cdr10m_n,sf_diag%cd10m_n,chr1p5m(:,n)                     &
     )
  END DO

END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE jules_land_sf_explicit
END MODULE jules_land_sf_explicit_mod





