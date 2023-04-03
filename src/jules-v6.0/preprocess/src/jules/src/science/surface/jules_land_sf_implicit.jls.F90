! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE jules_land_sf_implicit_mod

USE sf_melt_mod, ONLY: sf_melt
USE screen_tq_mod, ONLY: screen_tq
USE sf_evap_mod, ONLY: sf_evap
USE im_sf_pt2_mod, ONLY: im_sf_pt2
USE sice_htf_mod, ONLY: sice_htf
 
USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                       &
                  ModuleName='JULES_LAND_SF_IMPLICIT_MOD'

CONTAINS
!  SUBROUTINE JULES_LAND_SF_IMPLICIT --------------------------------
!
!  Purpose: Calculate implicit correction for land point to surface
!           fluxes of heat,moisture and momentum, to be used by 
!           the unconditionally stable and non-oscillatory BL 
!           numerical solver.
!
!--------------------------------------------------------------------
!    Arguments :-
SUBROUTINE jules_land_sf_implicit (                                           &
! IN values defining field dimensions and subset to be processed :
 land_pts,land_index,nsurft,surft_index,surft_pts,sm_levels,                  &
 canhc_surft,canopy,flake,smc_soilt,tile_frac,wt_ext_surft,fland,flandg,      &
! IN everything not covered so far :
 lw_down,sw_surft,t_soil_soilt,r_gamma,alpha1,ashtf_prime_surft,              &
 dtrdz_charney_grid_1,fraca,resfs,resft,rhokh_surft,                          &
 emis_surft,snow_surft,dtstar_surft,                                          &
! INOUT data :
 tstar_surft,fqw_surft,fqw_1,ftl_1,ftl_surft,sf_diag,                         &
! OUT Diagnostic not requiring STASH flags :
 ecan,ei_surft,esoil_surft,surf_ht_flux_land,ei_land,surf_htf_surft,          &
! OUT data required elsewhere in UM system :
 tstar_land,le_surft,radnet_surft,ecan_surft,esoil_soilt,                     &
 ext_soilt,snowmelt,melt_surft,tstar_surft_old,error,                         &
 !New arguments replacing USE statements
 ! lake_mod (IN)
 lake_h_ice_gb,                                                               &
 ! lake_mod (OUT)
 surf_ht_flux_lake_ij,                                                        &
 ! fluxes (IN)
 anthrop_heat_surft,                                                          &
 ! fluxes (OUT)
 surf_ht_store_surft,                                                         &
 ! c_elevate (IN)
 lw_down_elevcorr_surft,                                                      &
 ! prognostics (IN)
 nsnow_surft,                                                                 &
 ! jules_mod (IN)
 snowdep_surft,                                                               &
 ! JULES Types containing field data (IN OUT)
 crop_vars                                                                    &
)

!TYPE definitions
USE crop_vars_mod, ONLY: crop_vars_type

USE csigma,                   ONLY: sbcon

USE planet_constants_mod,     ONLY: cp

USE atm_fields_bounds_mod,    ONLY: tdims, pdims

USE theta_field_sizes,        ONLY: t_i_length, t_j_length

USE jules_surface_mod,        ONLY: l_aggregate, l_flake_model, ls

USE jules_snow_mod,           ONLY:                                           &
  nsmax, rho_snow_const, cansnowtile, l_snow_nocan_hc

USE jules_surface_types_mod,  ONLY: lake

USE sf_diags_mod,             ONLY: strnewsfdiag

USE timestep_mod,             ONLY: timestep

USE ancil_info,               ONLY: nsoilt

USE jules_surface_mod,        ONLY: l_neg_tstar

USE water_constants_mod,      ONLY: lc, lf, rho_ice

USE solinc_data,              ONLY: sky, l_skyview

USE parkind1,                 ONLY: jprb, jpim
USE yomhook,                  ONLY: lhook, dr_hook

USE jules_print_mgr,          ONLY: jules_message, jules_print


IMPLICIT NONE
!--------------------------------------------------------------------
!  Inputs :-
! (a) Defining horizontal grid and subset thereof to be processed.
!    Checked for consistency.
!--------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
 land_pts    ! IN No of land points

! (c) Soil/vegetation/land surface parameters (mostly constant).
INTEGER, INTENT(IN) ::                                                        &
 land_index(land_pts)        ! IN LAND_INDEX(I)=J => the Jth
                             !    point in ROW_LENGTH,ROWS is the
                             !    Ith land point.

INTEGER, INTENT(IN) ::                                                        &
 sm_levels                                                                    &
                             ! IN No. of soil moisture levels
,nsurft                                                                       &
                             ! IN No. of land tiles
,surft_index(land_pts,nsurft)                                                 &
                             ! IN Index of tile points
,surft_pts(nsurft)
                             ! IN Number of tile points

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
 canhc_surft(land_pts,nsurft)                                                 &
                             ! IN Areal heat capacity of canopy
                             !    for land tiles (J/K/m2).
,canopy(land_pts,nsurft)                                                      &
                             ! IN Surface/canopy water for
                             !    snow-free land tiles (kg/m2)
,flake(land_pts,nsurft)                                                       &
                             ! IN Lake fraction.
,smc_soilt(land_pts,nsoilt)                                                   &
                             ! IN Available soil moisture (kg/m2).
,tile_frac(land_pts,nsurft)                                                   &
                             ! IN Tile fractions including
                             !    snow cover in the ice tile.
,wt_ext_surft(land_pts,sm_levels,nsurft)                                      &
                             ! IN Fraction of evapotranspiration
                             !    extracted from each soil layer
                             !    by each tile.
,fland(land_pts)                                                              &
                             ! IN Land fraction on land pts.
,flandg(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                  &
                             ! IN Land fraction on all pts.
,emis_surft(land_pts,nsurft)                                                  &
                             ! IN Emissivity for land tiles
,snow_surft(land_pts,nsurft)                                                  &
                             ! IN Lying snow on tiles (kg/m2)
,dtstar_surft(land_pts,nsurft)
                             ! IN Change in TSTAR over timestep
                             !    for land tiles

! (f) Atmospheric + any other data not covered so far, incl control.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
 lw_down(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                             ! IN Surface downward LW radiation
                                  !    (W/m2).
,sw_surft(land_pts,nsurft)                                                    &
                             ! IN Surface net SW radiation on
                                  !    land tiles (W/m2).
,t_soil_soilt(land_pts,nsoilt,sm_levels)
                             ! IN Soil temperatures (K).

REAL(KIND=real_jlslsm), INTENT(IN) :: r_gamma
                             ! IN implicit weight in level 1

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
 alpha1(land_pts,nsurft)                                                      &
                             ! IN Mean gradient of saturated
                             !    specific humidity with respect
                             !    to temperature between the
                             !    bottom model layer and tile
                             !    surfaces
,ashtf_prime_surft(land_pts,nsurft)                                           &
                             ! IN Adjusted SEB coefficient for
                             !    land tiles.
,dtrdz_charney_grid_1(pdims%i_start:pdims%i_end,                              &
                      pdims%j_start:pdims%j_end)                              &
                             ! IN -g.dt/dp for model layers.
,fraca(land_pts,nsurft)                                                       &
                             ! IN Fraction of surface moisture
                             !    flux with only aerodynamic
                             !    resistance for snow-free land
                             !    tiles.
,resfs(land_pts,nsurft)                                                       &
                             ! IN Combined soil, stomatal
                             !    and aerodynamic resistance
                             !    factor for fraction (1-FRACA) of
                             !    snow-free land tiles.
,resft(land_pts,nsurft)                                                       &
                             ! IN Total resistance factor.
                             !    FRACA+(1-FRACA)*RESFS for
                             !    snow-free land, 1 for snow.
,rhokh_surft(land_pts,nsurft)
                             ! IN Surface exchange coefficients
                             !    for land tiles


!--------------------------------------------------------------------
!  In/outs :-
!--------------------------------------------------------------------
TYPE (strnewsfdiag), INTENT(INOUT) :: sf_diag
REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
 tstar_surft(land_pts,nsurft)                                                 &
                             ! INOUT Surface tile temperatures
,fqw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                             ! INOUT Moisture flux between layers
                             !       (kg per square metre per sec)
                             !       FQW(,1) is total water flux
                             !       from surface, 'E'.
,ftl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                   &
                             ! INOUT FTL(,K) contains net
                             !       turbulent sensible heat flux
                             !       into layer K from below; so
                             !       FTL(,1) is the surface
                             !       sensible heat, H.(W/m2)
,ftl_surft(land_pts,nsurft)                                                   &
                             ! INOUT Surface FTL for land tiles
,fqw_surft(land_pts,nsurft) 
                             ! INOUT Surface FQW for land tiles
  
!--------------------------------------------------------------------
!  Outputs :-
!-1 Diagnostic (or effectively so - includes coupled model requisites):-

!  (a) Calculated anyway (use STASH space from higher level) :-
!--------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
 ecan(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                    &
                             ! OUT Gridbox mean evaporation from
                             !     canopy/surface store (kg/m2/s).
                             !     Zero over sea.
,esoil_surft(land_pts,nsurft)                                                 &
                             ! OUT ESOIL for snow-free land tiles
,surf_ht_flux_land(tdims%i_start:tdims%i_end,                                 &
                   tdims%j_start:tdims%j_end)                                 &
                             ! OUT Net downward heat flux at
                             !     surface over land
                             !     fraction of gridbox (W/m2).
,ei_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                 &
                             ! OUT Sublimation from lying snow
                             !     (kg/m2/s).
,surf_htf_surft(land_pts,nsurft)
                             ! OUT Net downward surface heat flux
                             !     on tiles (W/m2)

!-2 Genuinely output, needed by other atmospheric routines :-

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
 tstar_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)              &
                             ! OUT   Land mean sfc temperature (K)
,le_surft(land_pts,nsurft)                                                    &
                             ! OUT Surface latent heat flux for
                             !     land tiles
,radnet_surft(land_pts,nsurft)                                                &
                             ! OUT Surface net radiation on
                             !     land tiles (W/m2)
,ei_surft(land_pts,nsurft)                                                    &
                             ! OUT EI for land tiles.
,ecan_surft(land_pts,nsurft)                                                  &
                             ! OUT ECAN for snow-free land tiles
,esoil_soilt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nsoilt)      &
                             ! OUT Surface evapotranspiration
                             !     from soil moisture store
                             !     (kg/m2/s).
,ext_soilt(land_pts,nsoilt,sm_levels)                                         &
                             ! OUT Extraction of water from each
                             !     soil layer (kg/m2/s).
,snowmelt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                &
                             ! OUT Snowmelt (kg/m2/s).
,melt_surft(land_pts,nsurft)                                                  &
                             ! OUT Snowmelt on land tiles (kg/m2/s
,tstar_surft_old(land_pts,nsurft)
                             ! OUT Tile surface temperatures at
                             !     beginning of timestep.

INTEGER, INTENT(OUT) ::                                                       &
 error                       ! OUT 0 - AOK;
                             !     1 to 7  - bad grid definition detected;

!New arguments replacing USE statements
! lake_mod (IN)
REAL(KIND=real_jlslsm), INTENT(IN) :: lake_h_ice_gb(land_pts)
! lake_mod (OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) :: surf_ht_flux_lake_ij(t_i_length,t_j_length)
! fluxes (IN)
REAL(KIND=real_jlslsm), INTENT(IN) :: anthrop_heat_surft(land_pts,nsurft)
! fluxes (OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) :: surf_ht_store_surft(land_pts,nsurft)
! c_elevate (IN)
REAL(KIND=real_jlslsm), INTENT(IN) :: lw_down_elevcorr_surft(land_pts,nsurft)
! prognostics (IN)
INTEGER, INTENT(IN) :: nsnow_surft(land_pts,nsurft)
! jules_mod (IN)
REAL(KIND=real_jlslsm), INTENT(IN) :: snowdep_surft(land_pts,nsurft)

!TYPES containing field data (IN OUT)
TYPE(crop_vars_type), INTENT(IN OUT) :: crop_vars

!--------------------------------------------------------------------
!  Workspace :-
!--------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                     &
 elake_surft(land_pts,nsurft)                                                 &
                             ! Lake evaporation.
,melt_ice_surft(land_pts,nsurft)                                              &
                             ! Ice melt on FLake lake tile (kg/m2/s)
,lake_ice_mass(land_pts)                                                      &
                             ! areal density equivalent to
                             ! lake ice of a given depth (kg/m2)
,non_lake_frac(tdims%i_start:tdims%i_end,                                     &
               tdims%j_start:tdims%j_end)
                             ! total tile fraction for surface types
                             ! other than inland water

REAL(KIND=real_jlslsm) ::                                                     &
 canhc_surf(land_pts)
                             ! Areal heat capacity of canopy
                             ! for land tiles (J/K/m2).

!  Local scalars :-

INTEGER ::                                                                    &
 i,j                                                                          &
            ! LOCAL Loop counter (horizontal field index).
,k                                                                            &
            ! LOCAL Tile pointer
,l                                                                            &
            ! LOCAL Land pointer
,n                                                                            &
            ! LOCAL Loop counter (tile index).
,m
            ! Loop counter for soil tiles

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='JULES_LAND_SF_IMPLICIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


!-----------------------------------------------------------------------

! Calculate surface scalar fluxes, temperatures only at the 1st call
! of the subroutine (first stage of the new BL solver) using standard
! MOSES2 physics and equations. These are the final values for this
! timestep and there is no need to repeat the calculation.
!-----------------------------------------------------------------------

error = 0


!$OMP PARALLEL                                                                &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(l,n,j,i,k)                                                      &
!$OMP SHARED(tdims,nsurft,surft_pts,surft_index,                              &
!$OMP        ftl_surft,nsoilt,land_pts,t_soil_soilt,                          &
!$OMP        tstar_surft_old,tstar_surft,dtstar_surft,cp,error)

!-----------------------------------------------------------------------
! 6.1 Convert FTL to sensible heat flux in Watts per square metre.
!-----------------------------------------------------------------------

DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
  DO k = 1,surft_pts(n)
    l = surft_index(k,n)
    ftl_surft(l,n) = cp * ftl_surft(l,n)
  END DO
!$OMP END DO NOWAIT
END DO

!-----------------------------------------------------------------------
! Land surface calculations
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Optional error check : test for negative top soil layer temperature
!-----------------------------------------------------------------------
IF (l_neg_tstar) THEN
  DO m = 1,nsoilt
!$OMP DO SCHEDULE(STATIC)
    DO l = 1,land_pts
      IF (t_soil_soilt(l,m,1) < 0) THEN
        error = 1
        WRITE(jules_message,*)                                                &
              '*** ERROR DETECTED BY ROUTINE JULES_LAND_SF_IMPLICIT ***'
        CALL jules_print('jules_land_sf_implicit_jls',jules_message)
        WRITE(jules_message,*) 'NEGATIVE TEMPERATURE IN TOP SOIL LAYER AT '
        CALL jules_print('jules_land_sf_implicit_jls',jules_message)
        WRITE(jules_message,*) 'LAND POINT ',l
        CALL jules_print('jules_land_sf_implicit_jls',jules_message)
      END IF
    END DO
!$OMP END DO NOWAIT
  END DO
END IF

!-----------------------------------------------------------------------
!   Diagnose the land surface temperature
!-----------------------------------------------------------------------

DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
  DO k = 1,surft_pts(n)
    l = surft_index(k,n)
    tstar_surft_old(l,n) = tstar_surft(l,n)
    tstar_surft(l,n) = tstar_surft_old(l,n) + dtstar_surft(l,n)
  END DO
!$OMP END DO NOWAIT
END DO

!$OMP END PARALLEL

!-----------------------------------------------------------------------
! 7.  Surface evaporation components and updating of surface
!     temperature (P245, routine SF_EVAP).
!-----------------------------------------------------------------------
CALL sf_evap (                                                                &
  land_pts,nsurft,                                                            &
  land_index,surft_index,surft_pts,sm_levels,fland,                           &
  ashtf_prime_surft,canopy,dtrdz_charney_grid_1,flake,fraca,                  &
  snow_surft,resfs,resft,rhokh_surft,tile_frac,smc_soilt,wt_ext_surft,        &
  timestep,r_gamma,fqw_1,fqw_surft,ftl_1,ftl_surft,tstar_surft,               &
  ecan,ecan_surft,elake_surft,esoil_soilt,esoil_surft,ei_surft,ext_soilt,     &
  sf_diag,                                                                    &
  ! crop_vars_mod (IN)
  crop_vars%frac_irr_soilt, crop_vars%frac_irr_surft,                         &
  crop_vars%wt_ext_irr_surft, crop_vars%resfs_irr_surft,                      &
  ! crop_vars_mod (IN OUT)
  crop_vars%smc_irr_soilt,                                                    &
  ! crop_vars_mod (OUT)
  crop_vars%ext_irr_soilt)

!-----------------------------------------------------------------------
!     Surface melting of sea-ice and snow on land tiles.
!-----------------------------------------------------------------------

!$OMP PARALLEL                                                                &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(l,n,j,i)                                                        &
!$OMP SHARED(tdims,ei_land,snowmelt,nsurft,land_pts,melt_ice_surft)

!$OMP DO SCHEDULE(STATIC)
DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    ei_land(i,j)  = 0.0
    snowmelt(i,j) = 0.0
  END DO
END DO
!$OMP END DO NOWAIT

! Lake initialisation
DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
  DO l = 1,land_pts
    melt_ice_surft(l,n) = 0.0
  END DO
!$OMP END DO NOWAIT
END DO

!$OMP END PARALLEL

DO n = 1,nsurft
  CALL sf_melt (                                                              &
    land_pts,land_index,                                                      &
    surft_index(:,n),surft_pts(n),flandg,                                     &
    alpha1(:,n),ashtf_prime_surft(:,n),dtrdz_charney_grid_1,                  &
    resft(:,n),rhokh_surft(:,n),tile_frac(:,n),timestep,r_gamma,              &
    ei_surft(:,n),fqw_1,ftl_1,fqw_surft(:,n),ftl_surft(:,n),                  &
    tstar_surft(:,n),snow_surft(:,n),snowdep_surft(:,n),                      &
    melt_surft(:,n)                                                           &
    )

  !-----------------------------------------------------------------------
  ! thermodynamic, flux contribution of melting ice on the FLake lake tile
  !-----------------------------------------------------------------------
  IF (     (l_flake_model   )                                                 &
    .AND. ( .NOT. l_aggregate)                                                &
    .AND. (n == lake       ) ) THEN

    ! lake_h_ice_gb is only initialised if FLake is on.

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(l)                                                              &
!$OMP SHARED(land_pts,lake_ice_mass,lake_h_ice_gb)
    DO l = 1, land_pts
      lake_ice_mass(l) = lake_h_ice_gb(l) * rho_ice
    END DO
!$OMP END PARALLEL DO

    CALL sf_melt (                                                            &
      land_pts,land_index,                                                    &
      surft_index(:,n),surft_pts(n),flandg,                                   &
      alpha1(:,n),ashtf_prime_surft(:,n),dtrdz_charney_grid_1,                &
      resft(:,n),rhokh_surft(:,n),tile_frac(:,n),timestep,r_gamma,            &
      ei_surft(:,n),fqw_1,ftl_1,fqw_surft(:,n),ftl_surft(:,n),                &
      tstar_surft(:,n),lake_ice_mass,lake_ice_mass / rho_snow_const,          &
      melt_ice_surft(:,n)                                                     &
        )
  END IF

  !-----------------------------------------------------------------------
  !  Increment snow by sublimation and melt
  !-----------------------------------------------------------------------

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,l,j,i)                                                        &
!$OMP SHARED(surft_pts,surft_index,land_index,t_i_length,ei_land,tile_frac,   &
!$OMP        ei_surft,snowmelt,melt_surft,n)
  DO k = 1,surft_pts(n)
    l = surft_index(k,n)
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    ei_land(i,j) = ei_land(i,j) + tile_frac(l,n) * ei_surft(l,n)
    snowmelt(i,j) = snowmelt(i,j) +                                           &
                    tile_frac(l,n) * melt_surft(l,n)
  END DO
!$OMP END PARALLEL DO

END DO

!$OMP PARALLEL                                                                &
!$OMP DEFAULT(SHARED)                                                         &
!$OMP PRIVATE(l,n,j,i,k)

IF (sf_diag%smlt) THEN
!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      sf_diag%snomlt_surf_htf(i,j) = lf * snowmelt(i,j)
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

!$OMP DO SCHEDULE(STATIC)
DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    surf_ht_flux_land(i,j) = 0.0
  END DO
END DO
!$OMP END DO NOWAIT

IF (     (l_flake_model   )                                                   &
    .AND. ( .NOT. l_aggregate) ) THEN
!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      surf_ht_flux_lake_ij(i,j) = 0.0
      ! initialise the non-lake fraction to one, not zero,
      ! in case there should ever be more than one lake tile, see below
      non_lake_frac(    i,j) = 1.0
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

!$OMP DO SCHEDULE(STATIC)
DO l = 1,land_pts
  j=(land_index(l) - 1) / t_i_length + 1
  i = land_index(l) - (j-1) * t_i_length
  tstar_land(i,j) = 0.0
END DO
!$OMP END DO NOWAIT

! initialise diagnostics to 0 to avoid packing problems
DO n = 1, nsurft
!$OMP DO SCHEDULE(STATIC)
  DO l = 1, land_pts
    radnet_surft(l,n) = 0.0
    le_surft(l,n) = 0.0
  END DO
!$OMP END DO NOWAIT
END DO

IF (sf_diag%l_lw_surft) THEN
  DO n = 1, nsurft
!$OMP DO SCHEDULE(STATIC)
    DO l = 1, land_pts
      sf_diag%lw_up_surft(l,n) = 0.0
      sf_diag%lw_down_surft(l,n) = 0.0
    END DO
!$OMP END DO NOWAIT
  END DO
END IF

!$OMP BARRIER

IF (l_skyview) THEN
  DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
    DO k = 1,surft_pts(n)
      l = surft_index(k,n)
      j=(land_index(l) - 1) / tdims%i_end + 1
      i = land_index(l) - (j-1) * tdims%i_end
      radnet_surft(l,n) = sw_surft(l,n) +   emis_surft(l,n) *                 &
        sky(i,j) * ( lw_down(i,j) + lw_down_elevcorr_surft(l,n)               &
                                - sbcon * tstar_surft(l,n)**4 )
    END DO
!$OMP END DO
  END DO
  IF (sf_diag%l_lw_surft) THEN
    DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
      DO k = 1,surft_pts(n)
        l = surft_index(k,n)
        j=(land_index(l) - 1) / tdims%i_end + 1
        i = land_index(l) - (j-1) * tdims%i_end
        sf_diag%lw_up_surft(l,n)   = emis_surft(l,n) * sky(i,j) *             &
                                     sbcon * tstar_surft(l,n)**4
        sf_diag%lw_down_surft(l,n) = emis_surft(l,n) * sky(i,j) *             &
                                     (lw_down(i,j) +                          &
                                     lw_down_elevcorr_surft(l,n))
      END DO
!$OMP END DO
    END DO
  END IF
ELSE
  DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
    DO k = 1,surft_pts(n)
      l = surft_index(k,n)
      j=(land_index(l) - 1) / tdims%i_end + 1
      i = land_index(l) - (j-1) * tdims%i_end
      radnet_surft(l,n) = sw_surft(l,n) +   emis_surft(l,n) *                 &
                 ( lw_down(i,j) + lw_down_elevcorr_surft(l,n)                 &
                                - sbcon * tstar_surft(l,n)**4 )
    END DO
!$OMP END DO
  END DO
  IF (sf_diag%l_lw_surft) THEN
    DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
      DO k = 1,surft_pts(n)
        l = surft_index(k,n)
        j=(land_index(l) - 1) / tdims%i_end + 1
        i = land_index(l) - (j-1) * tdims%i_end
        sf_diag%lw_up_surft(l,n)   = emis_surft(l,n) * sbcon *                &
                                     tstar_surft(l,n)**4
        sf_diag%lw_down_surft(l,n) = emis_surft(l,n) * (lw_down(i,j) +        &
                                     lw_down_elevcorr_surft(l,n))
      END DO
!$OMP END DO
    END DO
  END IF
END IF

DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
  DO k = 1,surft_pts(n)
    l = surft_index(k,n)
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    canhc_surf(l) = canhc_surft(l,n)
    IF ( ( .NOT. cansnowtile(n)) .AND. l_snow_nocan_hc .AND.                  &
         (nsmax > 0) .AND. (nsnow_surft(l,n) > 0) ) canhc_surf(l) = 0.0
    le_surft(l,n) = lc * ecan_surft(l,n) + lc * esoil_surft(l,n) +            &
                   lc * elake_surft(l,n) + ls * ei_surft(l,n)
    surf_ht_store_surft(l,n) = (canhc_surf(l) / timestep) *                   &
                         (tstar_surft(l,n) - tstar_surft_old(l,n))
    surf_htf_surft(l,n) = radnet_surft(l,n) + anthrop_heat_surft(l,n) -       &
                        ftl_surft(l,n) -                                      &
                        le_surft(l,n) -                                       &
                        lf * (melt_surft(l,n) + melt_ice_surft(l,n)) -        &
                        surf_ht_store_surft(l,n)
    ! separate out the lake heat flux for FLake
    ! and replace the snow-melt (NSMAX=0 only) and ice-melt heat fluxes
    ! so Flake can do its melting
    IF (     (l_flake_model   )                                               &
        .AND. ( .NOT. l_aggregate)                                            &
        .AND. (n == lake       ) ) THEN
      IF (nsmax == 0) THEN
        surf_ht_flux_lake_ij(i,j) = surf_htf_surft(l,n)                       &
                      + lf * (melt_surft(l,n) + melt_ice_surft(l,n))
        non_lake_frac(    i,j) = non_lake_frac(i,j) - tile_frac(l,n)
      ELSE
        surf_ht_flux_lake_ij(i,j) = surf_htf_surft(l,n)                       &
                      + lf * melt_ice_surft(l,n)
        non_lake_frac(    i,j) = non_lake_frac(i,j) - tile_frac(l,n)
      END IF
    ELSE
      surf_ht_flux_land(i,j) = surf_ht_flux_land(i,j)                         &
                        + tile_frac(l,n) * surf_htf_surft(l,n)
    END IF
    tstar_land(i,j) = tstar_land(i,j)                                         &
               + tile_frac(l,n) * tstar_surft(l,n)
  END DO
!$OMP END DO
END DO

  ! normalise the non-lake surface heat flux
IF (     (l_flake_model   )                                                   &
    .AND. ( .NOT. l_aggregate) ) THEN
!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      ! be careful about gridboxes that are all lake
      IF (non_lake_frac(i,j) > EPSILON(0.0)) THEN
        surf_ht_flux_land(i,j) =   surf_ht_flux_land(i,j)                     &
                                 / non_lake_frac(i,j)
      END IF
    END DO
  END DO
!$OMP END DO
END IF

IF (sf_diag%l_lh_land) THEN
  DO l = 1,land_pts
    sf_diag%lh_land(l) = SUM((tile_frac(l,:) * le_surft(l,:)))
  END DO
END IF

!-----------------------------------------------------------------------
! Optional error check : test for negative surface temperature
!-----------------------------------------------------------------------
IF (l_neg_tstar) THEN
!$OMP DO SCHEDULE(STATIC)
  DO l = 1,land_pts
    j=(land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    IF (tstar_land(i,j) < 0) THEN
      error = 1
      WRITE(jules_message,*)                                                  &
           '*** ERROR DETECTED BY ROUTINE JULES_LAND_SF_IMPLICIT ***'
      CALL jules_print('jules_land_sf_implicit_jls',jules_message)
      WRITE(jules_message,*) 'NEGATIVE SURFACE TEMPERATURE AT LAND POINT ',l
      CALL jules_print('jules_land_sf_implicit_jls',jules_message)
    END IF
  END DO
!$OMP END DO
END IF

!$OMP END PARALLEL


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE jules_land_sf_implicit
END MODULE jules_land_sf_implicit_mod
