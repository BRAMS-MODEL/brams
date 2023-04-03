! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE HYDROL--------------------------------------------------------

! Description:
!     Surface hydrology module which also updates the
!     sub-surface temperatures. Includes soil water phase
!     changes and the effect of soil water and ice on the
!     thermal and hydraulic characteristics of the soil.
!     This version is for use with MOSES II land surface
!     scheme.

! Documentation : UM Documentation Paper 25

MODULE hydrol_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='HYDROL_MOD'

CONTAINS

SUBROUTINE hydrol (                                                           &
    lice_pts, lice_index, soil_pts, soil_index,  nsnow_surft,                 &
    land_pts, sm_levels, bexp_soilt, catch_surft, con_rain_land,              &
    ecan_surft, ext_soilt, hcap_soilt, hcon_soilt, ls_rain_land,              &
    con_rainfrac_land,  ls_rainfrac_land,                                     &
    satcon_soilt, sathh_soilt, snowdepth_surft,  snow_soil_htf,               &
    surf_ht_flux_ld, timestep,                                                &
    smvcst_soilt, smvcwt_soilt, canopy_surft,                                 &
    stf_sub_surf_roff, smcl_soilt, sthf_soilt, sthu_soilt,                    &
    t_soil_soilt, tsurf_elev_surft, canopy_gb, smc_soilt, snow_melt,          &
    sub_surf_roff_gb, surf_roff_gb, tot_tfall_gb,                             &
    ! add new inland basin variable
    inlandout_atm_gb, l_inland,                                               &
    ! Additional variables for MOSES II
    nsurft, surft_pts, surft_index,                                           &
    infil_surft,  melt_surft, tile_frac,                                      &
    ! Additional variables required for large-scale hydrology:
    l_top, l_pdm, fexp_soilt, ti_mean_soilt, cs_ch4_soilt, cs_pool_soilt,     &
    dun_roff_soilt, drain_soilt, fsat_soilt, fwetl_soilt, qbase_soilt,        &
    qbase_l_soilt, qbase_zw_soilt, w_flux_soilt,                              &
    zw_soilt, sthzw_soilt, a_fsat_soilt, c_fsat_soilt, a_fwet_soilt,          &
    c_fwet_soilt,                                                             &
    resp_s_soilt, npp_soilt, fch4_wetl_soilt,                                 &
    fch4_wetl_cs_soilt, fch4_wetl_npp_soilt, fch4_wetl_resps_soilt,           &
    sthu_irr_soilt, frac_irr_soilt, ext_irr_soilt,                            &
    ! New arguments replacing USE statements
    ! trif_vars_mod
     n_leach_soilt, n_leach_gb_acc,                                           &
    ! prognostics
     n_inorg_soilt_lyrs, n_inorg_avail_pft, npft, t_soil_soilt_acc,           &
     tsoil_deep_gb,                                                           &
    ! top_pdm_alloc
     fch4_wetl_acc_soilt, substr_ch4, mic_ch4, mic_act_ch4, acclim_ch4,       &
    ! pdm_vars
     slope_gb, dim_cs1, l_soil_sat_down, asteps_since_triffid)

!Use in relevant subroutines
USE calc_zw_inund_mod,        ONLY: calc_zw_inund
USE ch4_wetl_mod,             ONLY: ch4_wetl
USE surf_hyd_mod,             ONLY: surf_hyd
USE soilt_mod,                ONLY: soilt
USE soilmc_mod,               ONLY: soilmc
USE soil_hyd_wt_mod,          ONLY: soil_hyd_wt
USE soil_hyd_update_mod,      ONLY: soil_hyd_update
USE soil_hyd_mod,             ONLY: soil_hyd
USE soil_htc_mod,             ONLY: soil_htc
USE ice_htc_mod,              ONLY: ice_htc
USE calc_baseflow_jules_mod,  ONLY: calc_baseflow_jules
USE n_leach_mod,              ONLY: n_leach

!Use in relevant variables
USE ancil_info,               ONLY:                                           &
  nsoilt

USE jules_hydrology_mod, ONLY:                                                &
  l_wetland_unfrozen, ti_max, zw_max

USE ancil_info, ONLY: dim_cslayer

USE elev_htc_mod, ONLY: elev_htc

USE jules_soil_biogeochem_mod, ONLY:                                          &
  ! imported scalar parameters
  soil_model_ecosse, soil_model_rothc, soil_model_1pool,                      &
  ! imported scalar variables
  soil_bgc_model, l_ch4_tlayered, l_layeredc, dim_ch4layer

USE jules_soil_mod, ONLY:                                                     &
  dzsoil,                                                                     &
     ! Thicknesses of the soil layers (m).
  dzsoil_elev, ns_deep

USE jules_vegetation_mod, ONLY: l_nitrogen

USE jules_irrig_mod, ONLY: l_irrig_dmd

USE jules_surface_mod, ONLY: l_elev_land_ice

USE ereport_mod, ONLY: ereport

USE water_constants_mod, ONLY:                                                &
  rho_water  ! Density of pure water (kg/m3).

USE update_mod,     ONLY: l_imogen
USE model_time_mod, ONLY: timestep_len

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Scalar arguments with INTENT(IN):
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  lice_pts,                                                                   &
    ! Number of land ice points.
  land_pts,                                                                   &
    ! Number of gridpoints.
  sm_levels,                                                                  &
    ! Number of soil moisture levels.
  soil_pts,                                                                   &
    ! Number of soil points.
  nsurft,                                                                     &
    ! Number of tiles
  dim_cs1,                                                                    &
    ! Number of soil carbon pools
  npft
    !

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  timestep
    ! Model timestep (s).

LOGICAL , INTENT(IN) ::                                                       &
  stf_sub_surf_roff,                                                          &
    ! Stash flag for sub-surface runoff.
  l_top,                                                                      &
    ! Flag for TOPMODEL-based hydrology.
  l_pdm,                                                                      &
    ! Flag for PDM hydrology.
  l_soil_sat_down
    ! Switch controlling direction of movement of
    ! soil moisture in excess of saturation.

!-----------------------------------------------------------------------------
! Array arguments with INTENT(IN):
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  lice_index(land_pts),                                                       &
    ! Array of land ice points.
  soil_index(land_pts),                                                       &
    ! Array of soil points.
  nsnow_surft(land_pts,nsurft)
    ! Number of snow layers

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  bexp_soilt(land_pts,nsoilt,sm_levels),                                      &
    ! Clapp-Hornberger exponent.
  catch_surft(land_pts,nsurft),                                               &
    ! Canopy/surface capacity of land tiles (kg/m2).
  con_rain_land(land_pts),                                                    &
    ! Convective rain (kg/m2/s).
  ecan_surft(land_pts,nsurft),                                                &
    ! Canopy evaporation from land tiles (kg/m2/s).
  ext_soilt(land_pts,nsoilt,sm_levels),                                       &
    ! Extraction of water from each soil layer (kg/m2/s).
  hcap_soilt(land_pts,nsoilt,sm_levels),                                      &
    ! Soil heat capacity (J/K/m3).
  hcon_soilt(land_pts,nsoilt,0:sm_levels),                                    &
    ! Soil thermal conductivity (W/m/K).
  ls_rain_land(land_pts),                                                     &
    ! Large-scale rain (kg/m2/s).
  con_rainfrac_land(land_pts),                                                &
    ! Convective rain fraction
  ls_rainfrac_land(land_pts),                                                 &
    ! large scale rain fraction
  satcon_soilt(land_pts,nsoilt,0:sm_levels),                                  &
    ! Saturated hydraulic conductivity (kg/m2/s).
  sathh_soilt(land_pts,nsoilt,sm_levels),                                     &
    ! Saturated soil water pressure (m).
  snow_melt(land_pts),                                                        &
    ! Snowmelt (kg/m2/s).
  snowdepth_surft(land_pts,nsurft),                                           &
    ! Snow depth (on ground) (m)
  snow_soil_htf(land_pts,nsurft),                                             &
    ! Tiled snowpack-> soil heat flux.
  surf_ht_flux_ld(land_pts),                                                  &
    ! Net downward surface heat flux (W/m2).
  smvcst_soilt(land_pts,nsoilt,sm_levels),                                    &
    ! Volumetric soil moisture concentration at saturation (m3 H2O/m3 soil).
  smvcwt_soilt(land_pts,nsoilt,sm_levels),                                    &
    ! Volumetric soil moisture concentration below which
    ! stomata close (m3 H2O/m3 soil).
  fexp_soilt(land_pts,nsoilt),                                                &
    ! Decay factor in Sat. Conductivity in deep LSH/TOPMODEL layer.
  ti_mean_soilt(land_pts,nsoilt),                                             &
    ! Mean topographic index.
  npp_soilt(land_pts,nsoilt),                                                 &
    ! Gridbox mean net primary productivity (kg C/m2/s).
  a_fsat_soilt(land_pts,nsoilt),                                              &
    ! Fitting parameter for Fsat in LSH model.
  c_fsat_soilt(land_pts,nsoilt),                                              &
    ! Fitting parameter for Fsat in LSH model.
  a_fwet_soilt(land_pts,nsoilt),                                              &
    ! Fitting parameter for Fwet in LSH model.
  c_fwet_soilt(land_pts,nsoilt)
    ! Fitting parameter for Fwet in LSH model.

!-----------------------------------------------------------------------------
! Array arguments with INTENT(INOUT):
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                     &
  cs_ch4_soilt(land_pts,nsoilt),                                              &
    ! Soil carbon used in CH4 wetlands if TRIFFID is switched off (kg C/m2).
  canopy_surft(land_pts,nsurft),                                              &
    ! Canopy water content for land tiles (kg/m2).
  smcl_soilt(land_pts,nsoilt,sm_levels),                                      &
    ! Soil moisture content of each
!                          !       layer (kg/m2).
  sthf_soilt(land_pts,nsoilt,sm_levels),                                      &
    ! Frozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
  sthu_soilt(land_pts,nsoilt,sm_levels),                                      &
    ! Unfrozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
  t_soil_soilt(land_pts,nsoilt,sm_levels),                                    &
    ! Sub-surface temperatures (K).
  tsurf_elev_surft(land_pts,nsurft),                                          &
    ! Tiled sub-surface temperatures (K).
  cs_pool_soilt(land_pts,nsoilt,dim_cslayer,dim_cs1),                         &
    ! Soil carbon (kg C/m2).
!                          !   For RothC (dim_cs1=4), the pools
!                          !   are DPM, RPM, biomass and humus.
  resp_s_soilt(land_pts,nsoilt,dim_cslayer,dim_cs1),                          &
    ! Soil respiration in pools (kg C/m2/s).
  fsat_soilt(land_pts,nsoilt),                                                &
     ! Surface saturation fraction.
  fwetl_soilt(land_pts,nsoilt),                                               &
     ! Wetland fraction.
  zw_soilt(land_pts,nsoilt),                                                  &
     ! Water table depth (m).
  sthzw_soilt(land_pts,nsoilt),                                               &
     ! Soil moisture fraction in deep LSH/TOPMODEL layer.
  substr_ch4(land_pts,dim_ch4layer),                                          &
    ! Dissolved substrate that methaogens consume (kg C/m2)
  mic_ch4(land_pts,dim_ch4layer),                                             &
    ! Methanogenic biomass (kg C/m2)
  mic_act_ch4(land_pts,dim_ch4layer),                                         &
    ! Activity level of methanogenic biomass (fraction)
  acclim_ch4(land_pts,dim_ch4layer),                                          &
    ! Acclimation factor for microbial trait adaptation
  sthu_irr_soilt(land_pts, nsoilt, sm_levels),                                &
  frac_irr_soilt(land_pts, nsoilt),                                           &
  ext_irr_soilt(land_pts, nsoilt, sm_levels)

!-----------------------------------------------------------------------------
! Array arguments with INTENT(OUT):
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  canopy_gb(land_pts),                                                        &
    ! Gridbox canopy water content (kg/m2).
  smc_soilt(land_pts,nsoilt),                                                 &
    ! Available soil moisture in a layer at the surface (kg/m2)
  sub_surf_roff_gb(land_pts),                                                 &
    ! Sub-surface runoff (kg/m2/s).
  surf_roff_gb(land_pts),                                                     &
    ! Surface runoff (kg/m2/s).
  tot_tfall_gb(land_pts),                                                     &
    ! Total throughfall (kg/m2/s).
  dun_roff_soilt(land_pts,nsoilt),                                            &
    ! Dunne part of sfc runoff (kg/m2/s).
  qbase_soilt(land_pts,nsoilt),                                               &
    ! Base flow (kg/m2/s).
  qbase_l_soilt(land_pts,nsoilt,sm_levels+1),                                 &
    ! Base flow from each level (kg/m2/s).
  qbase_zw_soilt(land_pts,nsoilt),                                            &
    ! Base flow from deep LSH/TOPMODEL layer (kg/m2/s).
  w_flux_soilt(land_pts,nsoilt,0:sm_levels),                                  &
    ! Fluxes of water between layers (kg/m2/s).
  drain_soilt(land_pts,nsoilt),                                               &
    ! Drainage out of sm_levels'th level (kg/m2/s).
  fch4_wetl_soilt(land_pts,nsoilt),                                           &
    ! Scaled wetland methane flux (default substrate) for use in
    ! atmos chemistry model (10^-9 kg C/m2/s).
  fch4_wetl_cs_soilt(land_pts,nsoilt),                                        &
    ! Scaled methane flux (soil carbon substrate) (kg C/m2/s).
  fch4_wetl_npp_soilt(land_pts,nsoilt),                                       &
    ! Scaled methane flux (npp substrate) (kg C/m2/s).
  fch4_wetl_resps_soilt(land_pts,nsoilt)
    ! Scaled methane flux (soil respiration substrate) (kg C/m2/s).

! Additional variables for MOSES II
INTEGER, INTENT(IN) ::                                                        &
  surft_pts(nsurft),                                                          &
    ! Number of tile points.
  surft_index(land_pts,nsurft)
    ! Index of tile points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  infil_surft(land_pts,nsurft),                                               &
    ! Maximum surface infiltration
  melt_surft(land_pts,nsurft),                                                &
    ! Snowmelt on tiles (kg/m2/s).
  tile_frac(land_pts,nsurft),                                                 &
    ! Tile fractions.
! Declare variable for inland basin outflow
  inlandout_atm_gb(land_pts)
    ! IN TRIP INLAND BASIN OUTFLOW FOR LAND POINTS ONLY,kg/m2/s=mm.

LOGICAL, INTENT(IN) ::                                                        &
  l_inland
    ! True if re-routing inland basin flow to soil moisture.

!-----------------------------------------------------------------------------
! New arguments replacing USE statements:
!-----------------------------------------------------------------------------
!trif_vars_mod
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
     n_leach_soilt(land_pts,nsoilt),                                          &
     n_leach_gb_acc(land_pts)
!prognostics
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
     n_inorg_soilt_lyrs(land_pts,nsoilt,dim_cslayer),                         &
     n_inorg_avail_pft(land_pts,npft,dim_cslayer),                            &
     t_soil_soilt_acc(land_pts,nsoilt,sm_levels)
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                     &
     tsoil_deep_gb(land_pts,ns_deep)
       ! Deep soil temperature (K).
!top_pdm
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: fch4_wetl_acc_soilt(land_pts,nsoilt)
!pdm_vars
REAL(KIND=real_jlslsm), INTENT(IN) :: slope_gb(land_pts)

!-----------------------------------------------------------------------------
! Local scalars:
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  i, j,                                                                       &
  n,                                                                          &
    ! Counter for soil level.
  m,                                                                          &
    ! Counter for soil tile.
  errorstatus

!-----------------------------------------------------------------------------
! Local arrays:
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                     &
  dsmc_dt_soilt(land_pts,nsoilt),                                             &
    ! Rate of change of soil moisture due to water falling onto the
    ! surface after surface runoff (kg/m2/s).
  ksz_soilt(land_pts,nsoilt,0:sm_levels),                                     &
    ! Saturated hydraulic conductivity in layer (kg/m2/s).
  qbase_unfr_soilt(land_pts,nsoilt),                                          &
    ! Base flow in unfrozen soil (kg/m2/s).
  qbase_l_unfr_soilt(land_pts,nsoilt,sm_levels+1),                            &
    ! As qbase_l but for unfrozen soil (kg/m2/s).
  top_crit_soilt(land_pts,nsoilt),                                            &
    ! Critical TI when ZW <=0.0
  dumtop_crit_soilt(land_pts,nsoilt),                                         &
    ! Dummy for top_crit_soilt
  dumsthf_soilt(land_pts,nsoilt,sm_levels),                                   &
    ! Dummy Frozen soil moisture content of each layer as a fraction of
    ! saturation (always set to 0).
  zdepth(0:sm_levels),                                                        &
    ! Lower soil layer boundary depth (m).
  tsoil_d_soilt(land_pts,nsoilt),                                             &
    ! Soil temperature in the top metre
  zw_inund_soilt(land_pts,nsoilt),                                            &
    ! Water table depth used
  wutot_soilt(land_pts,nsoilt),                                               &
    ! Ratio of unfrozen to total soil moisture at ZW.
  surf_roff_inc_soilt(land_pts,nsoilt),                                       &
    ! Increment to tiled surface runoff (kg m-2 s-1).
  surf_roff_soilt(land_pts,nsoilt),                                           &
    ! Soil-tiled contributions to surface runoff (kg m-2 s-1).
  sub_surf_roff_soilt(land_pts,nsoilt)
    ! Soil-tiled contributions to subsurface runoff (kg m-2 s-1).

REAL(KIND=real_jlslsm), PARAMETER :: to_kg_conversion = 1.0e-9
    ! multiplier for converting to kgC for wetland CH4 and IMOGEN

! Variables required for irrigation code
REAL(KIND=real_jlslsm) ::                                                     &
  w_flux_irr_soilt(land_pts,nsoilt,0:sm_levels),                              &
    ! The fluxes of water between layers in irrigated fraction (kg/m2/s).
  w_flux_nir_soilt(land_pts,nsoilt,0:sm_levels),                              &
    ! The fluxes of water between layers in non-irrigated fraction (kg/m2/s).
  smcl_irr_soilt(land_pts,nsoilt,sm_levels),                                  &
    ! Total soil moisture contents of each layer in irrigated
    ! fraction (kg/m2).
  smcl_nir_soilt(land_pts,nsoilt,sm_levels),                                  &
    ! Total soil moisture contents of each layer in non-irrigated
    ! fraction (kg/m2).
  sthu_nir_soilt(land_pts,nsoilt,sm_levels),                                  &
    ! Unfrozen soil moisture content of each layer as a fraction of
    ! saturation in irrigated fraction.
  ext_nir_soilt(land_pts,nsoilt,sm_levels),                                   &
    ! Extraction of water from each soil layer in non-irrigated fraction
    ! (kg/m2/s).
  smclsat_soilt(land_pts,nsoilt,sm_levels),                                   &
    ! The saturation moisture content of each layer (kg/m2).
  smclzw_soilt(land_pts,nsoilt),                                              &
    ! moisture content in deep layer(kg/m2).
  smclsatzw_soilt(land_pts,nsoilt)
    ! moisture content in deep layer (kg/m2).

INTEGER :: asteps_since_triffid

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='HYDROL'

! End of header---------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate soil carbon for use in the wetland CH4 scheme only
! (only used if single-pool C model is used):
IF ( soil_bgc_model == soil_model_1pool ) THEN
  DO m = 1,nsoilt
    DO j = 1,soil_pts
      i = soil_index(j)
      cs_ch4_soilt(i,m) = 0.0
      DO n = 1,dim_cslayer
        cs_ch4_soilt(i,m) = cs_ch4_soilt(i,m) + cs_pool_soilt(i,m,n,1)
      END DO
    END DO
  END DO
END IF

! Initialise w_flux variables that are used in irrigation code
w_flux_soilt(:,:,:)     = 0.0
w_flux_irr_soilt(:,:,:) = 0.0 ! to prevent random values reported over areas
                    ! that are not included as soil points (i.e. ice points)

!-----------------------------------------------------------------------------
! Set up variables required for LSH scheme:
!-----------------------------------------------------------------------------
zdepth(:) = 0.0

DO n = 1,sm_levels
  zdepth(n) = zdepth(n-1) + dzsoil(n)
END DO

! Initialise runoff increment.
surf_roff_inc_soilt(:,:) = 0.0

!-----------------------------------------------------------------------------
! Calculate throughfall and surface runoff, and update the canopy water
! content
!-----------------------------------------------------------------------------
CALL surf_hyd (land_pts, nsurft, surft_pts, surft_index,                      &
               catch_surft, ecan_surft, tile_frac, infil_surft, con_rain_land,&
               ls_rain_land,  con_rainfrac_land,  ls_rainfrac_land,           &
               melt_surft, snow_melt, timestep,                               &
               canopy_surft, canopy_gb, dsmc_dt_soilt,                        &
               l_top, l_pdm, sm_levels, soil_pts, soil_index,                 &
               surf_roff_gb, tot_tfall_gb,                                    &
               dun_roff_soilt, fsat_soilt, smvcst_soilt, sthu_soilt,          &
               sthf_soilt, surf_roff_soilt,                                   &
               ! New arguments to replace USE statements
               ! pdm_vars
               slope_gb)

!-----------------------------------------------------------------------------
! Specify the reduction of hydraulic conductivity with depth.
! Initialiase base flow to zero.
!-----------------------------------------------------------------------------
DO m = 1, nsoilt
  DO n = 0,sm_levels
    !CDIR NODEP
    DO j = 1,soil_pts
      i = soil_index(j)
      ksz_soilt(i,m,n) = satcon_soilt(i,m,n)
    END DO
  END DO
END DO

DO m = 1, nsoilt
  DO n = 1,sm_levels
    !CDIR NODEP
    DO j = 1,soil_pts
      smclsat_soilt(soil_index(j),m,n)      = rho_water * dzsoil(n) *         &
                                              smvcst_soilt(soil_index(j),m,n)
      qbase_l_soilt(soil_index(j),m,n)      = 0.0
      qbase_l_unfr_soilt(soil_index(j),m,n) = 0.0
      dumsthf_soilt(soil_index(j),m,n)      = 0.0
    END DO
  END DO
END DO

DO m = 1, nsoilt
  DO i = 1,land_pts
    qbase_soilt(i,m)      = 0.0
    qbase_zw_soilt(i,m)   = 0.0
    wutot_soilt(i,m)      = 0.0
    drain_soilt(i,m)      = 0.0
    qbase_unfr_soilt(i,m) = 0.0
    zw_inund_soilt(i,m)   = 0.0
  END DO
END DO

IF (l_top) THEN
  IF (soil_pts /= 0) THEN
    DO m = 1, nsoilt
      CALL calc_baseflow_jules(                                               &
        soil_pts, soil_index, land_pts, sm_levels,                            &
        zdepth, ksz_soilt(:,m,:),                                             &
        bexp_soilt(:,m,:), fexp_soilt(:,m), ti_mean_soilt(:,m), zw_soilt(:,m),&
        sthf_soilt(:,m,:),                                                    &
        dumtop_crit_soilt(:,m), qbase_soilt(:,m), qbase_l_soilt(:,m,:))
    END DO

    IF (l_wetland_unfrozen) THEN
      DO m = 1,nsoilt
        CALL calc_zw_inund(land_pts, sm_levels, soil_pts, soil_index, zdepth, &
          bexp_soilt(:,m,:), sathh_soilt(:,m,:), smclsat_soilt(:,m,:),        &
          smcl_soilt(:,m,:), sthu_soilt(:,m,:), sthzw_soilt(:,m),             &
          zw_soilt(:,m), zw_inund_soilt(:,m), wutot_soilt(:,m))

        ! Now call again to get the unfrozen equivalents to calculate fsat and
        ! fwet:
        CALL calc_baseflow_jules(                                             &
          soil_pts, soil_index, land_pts, sm_levels,                          &
          zdepth, ksz_soilt(:,m,:),                                           &
          bexp_soilt(:,m,:), fexp_soilt(:,m), ti_mean_soilt(:,m),             &
          zw_inund_soilt(:,m), dumsthf_soilt(:,m,:),                          &
          top_crit_soilt(:,m), qbase_unfr_soilt(:,m),                         &
          qbase_l_unfr_soilt(:,m,:))
      END DO
    ELSE
      top_crit_soilt(:,:) = dumtop_crit_soilt(:,:)
    END IF

  END IF
END IF  !  l_top

IF (l_inland) THEN
  DO i = 1,land_pts

    ! Add inland basin outflow to change in soil moisture store.
    ! Note for soil tiling- this is only used by the riv_intctl_1a, which is
    ! not compatible with nsoilt > 1.
    dsmc_dt_soilt(i,1) = dsmc_dt_soilt(i,1) + inlandout_atm_gb(i)
  END DO
END IF

!-----------------------------------------------------------------------------
! Update the layer soil moisture contents and calculate the
! gravitational drainage.
!-----------------------------------------------------------------------------
IF (soil_pts /= 0) THEN

  IF (l_irrig_dmd) THEN

    !-------------------------------------------------------------------------
    ! If l_irrig_dmd = TRUE, call soil_hyd separately for irrigated and
    ! non-irrigated fraction
    ! afterwards, call soil_hyd ONLY to update water table/drainage with
    ! gridbox total w_flux_soilt, smcl_soilt
    !-------------------------------------------------------------------------

    ! Split into irrigated / non-irrigated fraction.
    DO n = 1,sm_levels
      DO m = 1, nsoilt
        DO j = 1,soil_pts
          i = soil_index(j)
          smclsat_soilt(i,m,n) = rho_water * dzsoil(n) * smvcst_soilt(i,m,n)

          ! hadrd - gridbox sthu_soilt is assumed to be combination of
          ! sthu_soilt of non-irrigated fraction and sthu_irr_soilt, i.e.
          ! sthu_soilt = frac_irr_soilt*sthu_irr_soilt + (1-frac_irr_soilt)
          !              *sthu_nir_soilt
          IF ( frac_irr_soilt(i,m) < 1.0 ) THEN
            sthu_nir_soilt(i,m,n) = (sthu_soilt(i,m,n) - frac_irr_soilt(i,m)  &
                                    * sthu_irr_soilt(i,m,n))                  &
                                    / (1.0 - frac_irr_soilt(i,m))
            ext_nir_soilt(i,m,n)  = (ext_soilt(i,m,n)  - frac_irr_soilt(i,m)  &
                                    * ext_irr_soilt(i,m,n))                   &
                                    / (1.0 - frac_irr_soilt(i,m))
          ELSE
            sthu_nir_soilt(i,m,n) = sthu_soilt(i,m,n)
            ext_nir_soilt(i,m,n)  = ext_soilt(i,m,n)
          END IF

          smcl_irr_soilt(i,m,n) = smcl_soilt(i,m,n)                           &
                                  + (sthu_irr_soilt(i,m,n)                    &
                                     - sthu_soilt(i,m,n))                     &
                                  * smclsat_soilt(i,m,n)
          smcl_nir_soilt(i,m,n) = smcl_soilt(i,m,n)                           &
                                  + (sthu_nir_soilt(i,m,n)                    &
                                     - sthu_soilt(i,m,n))                     &
                                  * smclsat_soilt(i,m,n)
        END DO
      END DO
    END DO

    DO m = 1, nsoilt
      ! First call soil_hyd for non-irrigated fraction.
      CALL soil_hyd (                                                         &
        land_pts, sm_levels, soil_pts, soil_index, bexp_soilt(:,m,:), dzsoil, &
        ext_nir_soilt(:,m,:), dsmc_dt_soilt(:,m), ksz_soilt(:,m,:),           &
        sathh_soilt(:,m,:), timestep, smvcst_soilt(:,m,:),                    &
        smcl_nir_soilt(:,m,:), sthu_nir_soilt(:,m,:), w_flux_nir_soilt(:,m,:),&
        sthzw_soilt(:,m), zdepth, qbase_l_soilt(:,m,:), l_top,                &
        l_soil_sat_down, smclzw_soilt(:,m), smclsatzw_soilt(:,m),             &
        smclsat_soilt(:,m,:))

      ! Next call soil_hyd for irrigated fraction.
      CALL soil_hyd (                                                         &
        land_pts, sm_levels, soil_pts, soil_index, bexp_soilt(:,m,:), dzsoil, &
        ext_irr_soilt(:,m,:), dsmc_dt_soilt(:,m), ksz_soilt(:,m,:),           &
        sathh_soilt(:,m,:), timestep, smvcst_soilt(:,m,:),                    &
        smcl_irr_soilt(:,m,:), sthu_irr_soilt(:,m,:), w_flux_irr_soilt(:,m,:),&
        sthzw_soilt(:,m), zdepth, qbase_l_soilt(:,m,:), l_top,                &
        l_soil_sat_down, smclzw_soilt(:,m), smclsatzw_soilt(:,m),             &
        smclsat_soilt(:,m,:))
    END DO

    ! Re-calculate total grid box soil moisture.
    ! hadrd - perhaps this could be done more efficiently with WHERE
    DO m = 1,nsoilt
      DO n = 0,sm_levels
        DO j = 1,soil_pts
          i = soil_index(j)

          ! Ensure sensible values if irrigation fraction is very small:
          IF ( frac_irr_soilt(i,m) <= EPSILON(1.0) ) THEN
            w_flux_irr_soilt(i,m,n) = 0.0
            w_flux_soilt(i,m,n)     = w_flux_nir_soilt(i,m,n)
            IF ( n >= 1 ) THEN
              sthu_irr_soilt(i,m,n) = sthu_nir_soilt(i,m,n)
              smcl_soilt(i,m,n)     = smcl_nir_soilt(i,m,n)
              sthu_soilt(i,m,n)     = sthu_nir_soilt(i,m,n)
            END IF
          ELSE
            w_flux_soilt(i,m,n) = frac_irr_soilt(i,m)                         &
                                  * w_flux_irr_soilt(i,m,n)                   &
                                  + ( 1.0 - frac_irr_soilt(i,m) )             &
                                  * w_flux_nir_soilt(i,m,n)
            IF ( n >= 1 ) THEN
              smcl_soilt(i,m,n) = frac_irr_soilt(i,m) * smcl_irr_soilt(i,m,n) &
                                  + ( 1.0 - frac_irr_soilt(i,m) )             &
                                  * smcl_nir_soilt(i,m,n)
              sthu_soilt(i,m,n) = frac_irr_soilt(i,m) * sthu_irr_soilt(i,m,n) &
                                  + ( 1.0 - frac_irr_soilt(i,m) )             &
                                  * sthu_nir_soilt(i,m,n)
            END IF
          END IF
        END DO  !  soil points
      END DO  !  layers
    END DO  !  tiles

    DO m = 1, nsoilt
      CALL soil_hyd_update(land_pts, sm_levels, soil_pts, soil_index, dzsoil, &
                           smvcst_soilt(:,m,:), zdepth, smclzw_soilt(:,m),    &
                           sthzw_soilt(:,m), smclsat_soilt(:,m,:),            &
                           smclsatzw_soilt(:,m))
    END DO

  ELSE
    ! .NOT. l_irrig_dmd

    DO m = 1, nsoilt
      CALL soil_hyd (                                                         &
        land_pts, sm_levels, soil_pts, soil_index, bexp_soilt(:,m,:), dzsoil, &
        ext_soilt(:,m,:), dsmc_dt_soilt(:,m), ksz_soilt(:,m,:),               &
        sathh_soilt(:,m,:), timestep, smvcst_soilt(:,m,:),                    &
        smcl_soilt(:,m,:), sthu_soilt(:,m,:), w_flux_soilt(:,m,:),            &
        sthzw_soilt(:,m), zdepth, qbase_l_soilt(:,m,:), l_top,                &
        l_soil_sat_down, smclzw_soilt(:,m), smclsatzw_soilt(:,m),             &
        smclsat_soilt(:,m,:))
    END DO

  END IF  !  l_irrig_dmd

  IF (nsoilt == 1) THEN
    !To maintain bit-comparability, we need to call with the _gb version of
    !surface runoff.
    m = 1
    CALL soil_hyd_wt (                                                        &
      land_pts, sm_levels, soil_pts, soil_index,                              &
      bexp_soilt(:,m,:), dsmc_dt_soilt(:,m), sathh_soilt(:,m,:), timestep,    &
      smvcst_soilt(:,m,:), sub_surf_roff_soilt(:,m), smcl_soilt(:,m,:),       &
      surf_roff_gb,                                                           &
      w_flux_soilt(:,m,:), stf_sub_surf_roff, zw_soilt(:,m),                  &
      sthzw_soilt(:,m), qbase_soilt(:,m), qbase_l_soilt(:,m,:),               &
      drain_soilt(:,m), l_top, smclzw_soilt(:,m), smclsatzw_soilt(:,m),       &
      smclsat_soilt(:,m,:), surf_roff_inc_soilt(:,m) )

    ! For a single soil tile, simply copy across to the output variable.
    sub_surf_roff_gb(:) = sub_surf_roff_soilt(:,m)
    ! Update the tiled surface runoff.
    surf_roff_soilt(:,m) = surf_roff_soilt(:,m) + surf_roff_inc_soilt(:,m)
  ELSE

    ! Initialise output variable
    sub_surf_roff_gb(:) = 0.0
    DO m = 1, nsoilt
      CALL soil_hyd_wt (                                                      &
        land_pts, sm_levels, soil_pts, soil_index,                            &
        bexp_soilt(:,m,:), dsmc_dt_soilt(:,m), sathh_soilt(:,m,:), timestep,  &
        smvcst_soilt(:,m,:), sub_surf_roff_soilt(:,m), smcl_soilt(:,m,:),     &
        surf_roff_soilt(:,m),                                                 &
        w_flux_soilt(:,m,:), stf_sub_surf_roff, zw_soilt(:,m),                &
        sthzw_soilt(:,m), qbase_soilt(:,m), qbase_l_soilt(:,m,:),             &
        drain_soilt(:,m), l_top, smclzw_soilt(:,m), smclsatzw_soilt(:,m),     &
        smclsat_soilt(:,m,:), surf_roff_inc_soilt(:,m) )

      ! For multiple soil tiles, add up the contributions, allowing for frac
      sub_surf_roff_gb(:) = sub_surf_roff_gb(:)                               &
                            + ( tile_frac(:,m) * sub_surf_roff_soilt(:,m))
      surf_roff_gb(:) = surf_roff_gb(:) + tile_frac(:,m)                      &
                                          * surf_roff_inc_soilt(:,m)

    END DO
  END IF !nsoilt == 1

  !---------------------------------------------------------------------------
  ! Calculate surface saturation and wetland fractions:
  !---------------------------------------------------------------------------
  IF (l_top) THEN
    DO m = 1, nsoilt

      DO i = 1,land_pts
        fsat_soilt(i,m)  = 0.0
        fwetl_soilt(i,m) = 0.0

        ! Zero soil porosity over land ice:
        IF (smvcst_soilt(i,m,sm_levels) <= 0.0) THEN
          zw_soilt(i,m) = zw_max
        END IF
      END DO

      DO j = 1,soil_pts
        i = soil_index(j)
        qbase_zw_soilt(i,m) = qbase_l_soilt(i,m,sm_levels+1)

        !Now use fit for fsat_soilt and fwet:
        IF (l_wetland_unfrozen) THEN
          fsat_soilt(i,m)  = wutot_soilt(i,m) * a_fsat_soilt(i,m)             &
                             * EXP(-c_fsat_soilt(i,m) * top_crit_soilt(i,m))
          fwetl_soilt(i,m) = wutot_soilt(i,m) * a_fwet_soilt(i,m)             &
                             * EXP(-c_fwet_soilt(i,m) * top_crit_soilt(i,m))
        ELSE
          fsat_soilt(i,m)  = a_fsat_soilt(i,m)                                &
                             * EXP(-c_fsat_soilt(i,m) * top_crit_soilt(i,m))
          fwetl_soilt(i,m) = a_fwet_soilt(i,m)                                &
                             * EXP(-c_fwet_soilt(i,m) * top_crit_soilt(i,m))
        END IF

        IF (top_crit_soilt(i,m) >= ti_max) THEN
          fsat_soilt(i,m)  = 0.0
          fwetl_soilt(i,m) = 0.0
        END IF

      END DO
    END DO
  END IF  !  l_top

ELSE ! soil pts

  !---------------------------------------------------------------------------
  ! If required by STASH flag and there are no soil points,
  ! set sub-surface runoff to zero.
  !---------------------------------------------------------------------------
  IF (stf_sub_surf_roff) THEN
    DO i = 1,land_pts
      sub_surf_roff_gb(i) = 0.0
    END DO
  END IF

END IF  !  soil_pts

!-----------------------------------------------------------------------------
! Update the soil temperatures and the frozen moisture fractions
!-----------------------------------------------------------------------------

!=============================================================================
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
!=============================================================================
IF (soil_pts /= 0) THEN
  ! When using soil tiling, we can use the surface tiled version of
  ! surf_ht_flux_ld, snow_soil_htf. The _ld version is a gridbox mean
  ! calculated in snow and passed through.
  IF (nsoilt == 1) THEN
    m = 1
    CALL soil_htc (                                                           &
      land_pts, sm_levels, nsurft, soil_pts, soil_index,                      &
      surft_pts, surft_index, nsnow_surft,                                    &
      bexp_soilt(:,m,:), dzsoil, tile_frac, hcap_soilt(:,m,:),                &
      hcon_soilt(:,m,:),                                                      &
      sathh_soilt(:,m,:), surf_ht_flux_ld, timestep, smvcst_soilt(:,m,:),     &
      w_flux_soilt(:,m,:), sthu_irr_soilt(:,m,:), smcl_soilt(:,m,:),          &
      snowdepth_surft, sthu_soilt(:,m,:), sthf_soilt(:,m,:),                  &
      t_soil_soilt(:,m,:),                                                    &
      !New arguments to replace USE statements
      ! prognostics
      tsoil_deep_gb )
  ELSE
    ! Surface and soil tiles map directly on to each other.
    DO m = 1, nsoilt
      n = m
      CALL soil_htc (                                                         &
        land_pts, sm_levels, nsurft, soil_pts, soil_index,                    &
        surft_pts, surft_index, nsnow_surft,                                  &
        bexp_soilt(:,m,:), dzsoil, tile_frac, hcap_soilt(:,m,:),              &
        hcon_soilt(:,m,:),                                                    &
        sathh_soilt(:,m,:), snow_soil_htf(:,n), timestep, smvcst_soilt(:,m,:),&
        w_flux_soilt(:,m,:), sthu_irr_soilt(:,m,:), smcl_soilt(:,m,:),        &
        snowdepth_surft, sthu_soilt(:,m,:), sthf_soilt(:,m,:),                &
        t_soil_soilt(:,m,:),                                                  &
        !New arguments to replace USE statements
        ! prognostics
        tsoil_deep_gb )
    END DO
  END IF
END IF
!=============================================================================
! *END NOTICE REGARDING SOIL TILING**
!=============================================================================

!-----------------------------------------------------------------------------
! Update the sub-surface temperatures for land ice.
!-----------------------------------------------------------------------------
IF (lice_pts /= 0) THEN
  IF ( .NOT. l_elev_land_ice) THEN
    DO m = 1, nsoilt
      CALL ice_htc (land_pts, sm_levels, lice_pts, lice_index, dzsoil,        &
                    surf_ht_flux_ld, timestep,                                &
                    t_soil_soilt(:,m,:))
    END DO
  ELSE
    CALL elev_htc (land_pts, lice_pts, lice_index, nsurft,                    &
                   dzsoil_elev, snow_soil_htf, timestep,                      &
                   tsurf_elev_surft)
  END IF
END IF

!-----------------------------------------------------------------------------
! Diagnose the available soil moisture in a layer at the surface.
!-----------------------------------------------------------------------------
DO m = 1, nsoilt
  CALL soilmc ( land_pts,sm_levels,soil_pts,soil_index,                       &
                dzsoil,sthu_soilt(:,m,:),smvcst_soilt(:,m,:),                 &
                smvcwt_soilt(:,m,:),smc_soilt(:,m) )
END DO

!-----------------------------------------------------------------------------
! Calculate mean soil temperature and scaled CH4 flux:
!-----------------------------------------------------------------------------
DO m = 1, nsoilt
  DO i = 1,land_pts
    fch4_wetl_soilt(i,m)       = 0.0
    fch4_wetl_cs_soilt(i,m)    = 0.0
    fch4_wetl_npp_soilt(i,m)   = 0.0
    fch4_wetl_resps_soilt(i,m) = 0.0
  END DO
  IF ( l_top .AND. soil_pts /= 0 ) THEN
    SELECT CASE ( soil_bgc_model )
    CASE ( soil_model_1pool, soil_model_rothc )
      IF ( l_ch4_tlayered ) THEN
        ! This variable is not used with layered CH4 calc.
        tsoil_d_soilt(:,m) = 0.0
      ELSE
        CALL soilt(land_pts, sm_levels, soil_pts, soil_index,                 &
                   dzsoil, t_soil_soilt(:,m,:), tsoil_d_soilt(:,m))
      END IF
      CALL ch4_wetl(land_pts, soil_pts, dim_cs1, soil_index, sm_levels,       &
                    tsoil_d_soilt(:,m), cs_ch4_soilt(:,m),                    &
                    t_soil_soilt(:,m,:), cs_pool_soilt(:,m,:,:),              &
                    resp_s_soilt(:,m,:,:), npp_soilt(:,m), fwetl_soilt(:,m),  &
                    fch4_wetl_soilt(:,m), fch4_wetl_cs_soilt(:,m),            &
                    fch4_wetl_npp_soilt(:,m), fch4_wetl_resps_soilt(:,m),     &
                    substr_ch4, mic_ch4, mic_act_ch4, acclim_ch4,             &
                    sthu_soilt(:,m,:), bexp_soilt(:,m,:), timestep,           &
                    l_ch4_tlayered)
      IF (l_imogen) THEN
        DO i = 1,land_pts
          ! fch4_wetl_acc_soilt in kg/m2 and fch4_wetl_soilt in 10e9kg/m2/s
          fch4_wetl_acc_soilt(i,m) = fch4_wetl_acc_soilt(i,m) +               &
                (fch4_wetl_soilt(i,m) * to_kg_conversion * timestep_len)
        END DO
      END IF
    END SELECT
  END IF
END DO  !  tiles

IF ( (soil_bgc_model == soil_model_rothc .AND. l_layeredc )                   &
     .OR. soil_bgc_model == soil_model_ecosse ) THEN
  !-----------------------------------------------------------------------
  !accumulate soil temperature for layered soil carbon and nitrogen
  !-----------------------------------------------------------------------
  DO m = 1, nsoilt
    IF (asteps_since_triffid == 1) THEN
      t_soil_soilt_acc(:,m,:) = 0.0
    END IF
    DO j = 1,soil_pts
      i = soil_index(j)
      t_soil_soilt_acc(i,m,:) = t_soil_soilt_acc(i,m,:) + t_soil_soilt(i,m,:)
    END DO
  END DO
END IF

!-----------------------------------------------------------------------
! Calculate Nitrogen Leaching
!-----------------------------------------------------------------------
IF (soil_bgc_model == soil_model_rothc .AND. l_nitrogen) THEN
  IF (nsoilt > 1) THEN
    errorstatus = 101
    CALL ereport("check hydrol_jls", errorstatus,                             &
                 "nsoilt>1 and l_nitrogen - n_leach not currently coded")
  END IF
  CALL n_leach(land_pts, timestep, smcl_soilt, w_flux_soilt,                  &
               sub_surf_roff_gb, qbase_l_soilt, asteps_since_triffid,         &
               ! New arguments replacing USE statements
               ! trif_vars_mod
                n_leach_soilt, n_leach_gb_acc,                                &
               ! prognostics
                n_inorg_soilt_lyrs, n_inorg_avail_pft, dim_cslayer)
END IF

!-----------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE hydrol
END MODULE hydrol_mod
