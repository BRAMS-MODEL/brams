! *****************************COPYRIGHT*******************************

! (c) [University of Edinburgh]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
!  SUBROUTINE SNOW-----------------------------------------------------

! Description:
!     Calling routine for snow module

! Subroutine Interface:
MODULE snow_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SNOW_MOD'

CONTAINS

SUBROUTINE snow ( land_pts, timestep, stf_hf_snow_melt, nsurft,               &
                  surft_pts, surft_index, catch_snow, con_snow,               &
                  con_rain, tile_frac, ls_snow, ls_graup, ls_rain, ei_surft,  &
                  hcaps1_soilt, hcons, melt_surft,                            &
                  smcl1_soilt, sthf1_soilt, surf_htf_surft,                   &
                  t_soil1_soilt, tsurf_elev_surft, tstar_surft,               &
                  smvcst1_soilt, rgrain, rgrainl, rho_snow_grnd, sice,        &
                  sliq, snow_grnd, snow_surft, snowdepth,                     &
                  tsnow, nsnow, ds, hf_snow_melt, lying_snow,                 &
                  rho_snow, snomlt_sub_htf, snow_melt,                        &
                  snow_soil_htf,  surf_ht_flux_ld,  sf_diag,                  &
                  dhf_surf_minus_soil,                                        &
                  ! New Arguments to replace USE statements
                  ! jules_internal
                  unload_backgrnd_pft, npft,                                  &
                  !Ancil info (IN)
                  l_lice_point)

USE canopysnow_mod,  ONLY: canopysnow
USE compactsnow_mod, ONLY: compactsnow
USE layersnow_mod,   ONLY: layersnow
USE relayersnow_mod, ONLY: relayersnow
USE snowgrain_mod,   ONLY: snowgrain
USE snowpack_mod,    ONLY: snowpack
USE snowtherm_mod,   ONLY: snowtherm

USE water_constants_mod, ONLY:                                                &
  ! imported scalar parameters
   lf
     ! Latent heat of fusion of water at 0degc (J kg-1).

USE jules_snow_mod, ONLY:                                                     &
  nsmax,                                                                      &
    ! Maximum possible number of snow layers.
  l_snow_infilt,                                                              &
    ! Include infiltration of rain into the snowpack.
  graupel_options,                                                            &
    ! Switch for treatment of graupel in the snow scheme.
  graupel_as_snow,                                                            &
    ! Include graupel in the surface snowfall.
  ignore_graupel,                                                             &
    ! Ignore graupel in the surface snowfall.
  r0,                                                                         &
    ! Grain size for fresh snow (microns).
  cansnowtile
    ! Switch for canopy snow model.

USE jules_surface_types_mod, ONLY: lake

USE jules_radiation_mod, ONLY: l_snow_albedo, l_embedded_snow

USE sf_diags_mod, ONLY: strnewsfdiag

USE ancil_info, ONLY: nsoilt

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  land_pts              ! Total number of land points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  timestep              ! Timestep length (s).

LOGICAL, INTENT(IN) ::                                                        &
  stf_hf_snow_melt      ! Stash flag for snowmelt heat flux.

!-----------------------------------------------------------------------------
! Array arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  nsurft,                                                                     &
    ! Number of land tiles.
  surft_pts(nsurft),                                                          &
    ! Number of tile points.
  surft_index(land_pts,nsurft)
    ! Index of tile points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  catch_snow(land_pts,nsurft),                                                &
    ! Canopy snow capacity (kg/m2).
  con_snow(land_pts),                                                         &
    ! Convective snowfall rate (kg/m2/s).
  tile_frac(land_pts,nsurft),                                                 &
    ! Tile fractions.
  ei_surft(land_pts,nsurft),                                                  &
    ! Sublimation of snow (kg/m2/s).
  hcaps1_soilt(land_pts,nsoilt),                                              &
    ! Soil heat capacity of top layer(J/K/m3).
  hcons(land_pts),                                                            &
    ! Thermal conductivity of top soil layer,
    ! including water and ice (W/m/K).
  smcl1_soilt(land_pts,nsoilt),                                               &
    ! Moisture content of surface soil layer (kg/m2).
  sthf1_soilt(land_pts,nsoilt),                                               &
    ! Frozen soil moisture content of surface layer as a fraction
    ! of saturation.
  surf_htf_surft(land_pts,nsurft),                                            &
    ! Surface heat flux (W/m2).
  tstar_surft(land_pts,nsurft),                                               &
    ! Tile surface temperature (K).
  smvcst1_soilt(land_pts,nsoilt)
    ! Surface soil layer volumetric moisture concentration at saturation.

!-----------------------------------------------------------------------------
! Array arguments with intent(inout)
!-----------------------------------------------------------------------------
TYPE (strnewsfdiag), INTENT(INOUT) :: sf_diag

INTEGER, INTENT(INOUT) ::                                                     &
  nsnow(land_pts,nsurft)   ! Number of snow layers.

REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
  con_rain(land_pts),                                                         &
    ! Convective rainfall rate (kg/m2/s).
  ls_snow(land_pts),                                                          &
    ! Large-scale frozen precip fall rate (kg/m2/s).
  ls_graup(land_pts),                                                         &
    ! Large-scale graupel fall rate (kg/m2/s).
  ls_rain(land_pts),                                                          &
    ! Large-scale rainfall rate (kg/m2/s).
  melt_surft(land_pts,nsurft),                                                &
    ! Surface or canopy snowmelt rate (kg/m2/s).
    ! On output, this is the total melt rate for the tile
    ! (i.e. sum of  melt on canopy and ground).
  t_soil1_soilt(land_pts,nsoilt),                                             &
    ! Soil surface layer temperature (K).
  tsurf_elev_surft(land_pts,nsurft),                                          &
    ! Temperature of elevated subsurface tiles (K).
  rgrain(land_pts,nsurft),                                                    &
    ! Snow surface grain size (microns).
  rgrainl(land_pts,nsurft,nsmax),                                             &
    ! Snow layer grain size (microns).
  rho_snow_grnd(land_pts,nsurft),                                             &
    ! Snowpack bulk density (kg/m3).
  sice(land_pts,nsurft,nsmax),                                                &
    ! Ice content of snow layers (kg/m2).
  sliq(land_pts,nsurft,nsmax),                                                &
    ! Liquid content of snow layers (kg/m2).
  snow_grnd(land_pts,nsurft),                                                 &
    ! Snow beneath canopy (kg/m2).
  snow_surft(land_pts,nsurft),                                                &
    ! Snow mass (kg/m2).
  snowdepth(land_pts,nsurft),                                                 &
    ! Snow depth (m).
  tsnow(land_pts,nsurft,nsmax)
    ! Snow layer temperatures (K).

!-----------------------------------------------------------------------------
! Array arguments with intent(out)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  ds(land_pts,nsurft,nsmax),                                                  &
    ! Snow layer thicknesses (m).
  hf_snow_melt(land_pts),                                                     &
    ! Gridbox snowmelt heat flux (W/m2).
  lying_snow(land_pts),                                                       &
    ! Gridbox snow mass (kg m-2).
  rho_snow(land_pts,nsurft,nsmax),                                            &
    ! Snow layer densities (kg/m3).
  snomlt_sub_htf(land_pts),                                                   &
    ! Sub-canopy snowmelt heat flux (W/m2).
  snow_melt(land_pts),                                                        &
    ! Gridbox snowmelt (kg/m2/s).
  surf_ht_flux_ld(land_pts),                                                  &
    ! Surface heat flux on land (W/m2).
  snow_soil_htf(land_pts,nsurft),                                             &
    ! Heat flux into the uppermost subsurface layer (W/m2).
    ! i.e. snow to ground, or into snow/soil composite layer.
  dhf_surf_minus_soil(land_pts)
    ! Heat flux difference across the FLake snowpack (W/m2).

!-----------------------------------------------------------------------------
! New arguments to replace USE statements
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) :: npft
REAL(KIND=real_jlslsm), INTENT(IN) :: unload_backgrnd_pft(land_pts,npft)

!ancil_info (IN)
LOGICAL, INTENT(IN) :: l_lice_point(land_pts)

!-----------------------------------------------------------------------------
! Local scalars
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  i,                                                                          &
    ! Land point index and loop counter.
  j,                                                                          &
    ! Tile pts loop counter.
  k,                                                                          &
    ! Tile number.
  n,                                                                          &
    ! Tile loop counter.
  m,                                                                          &
    ! Soil tile index.
  ns,                                                                         &
    ! Snow layer index
  n_can
    ! Actual or dummy pointer to array defined only on PFTs.

!-----------------------------------------------------------------------------
! Local arrays
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                     &
  csnow(land_pts,nsmax),                                                      &
    ! Areal heat capacity of layers (J/K/m2).
  ksnow(land_pts,nsmax),                                                      &
    ! Thermal conductivity of layers (W/m/K).
  rho0(land_pts),                                                             &
    ! Density of fresh snow (kg/m3).
    ! Where nsnow=0, rho0 is the density of the snowpack.
  snowfall(land_pts),                                                         &
    ! Total frozen precip reaching the ground in timestep
    ! (kg/m2) - includes any canopy unloading.
  graupfall(land_pts),                                                        &
    ! Graupel reaching the ground in timestep(kg/m2).
  infiltration(land_pts),                                                     &
    ! Infiltration of rainfall into snow pack
    ! on the current tile (kg/m2).
  infil_rate_con_gbm(land_pts),                                               &
    ! Grid-box mean rate of infiltration of convective rainfall into
    ! snowpack (kg/m2/s).
  infil_rate_ls_gbm(land_pts),                                                &
    ! Grid-box mean rate of infiltration of large-scale rainfall into
    ! snowpack (kg/m2/s).
  snowmass(land_pts),                                                         &
    ! Snow mass on the ground (kg/m2).
  rgrain0(land_pts),                                                          &
    ! Fresh snow grain size (microns).
  sice0(land_pts),                                                            &
   ! Ice content of fresh snow (kg/m2).
   ! Where nsnow=0, sice0 is the mass of the snowpack.
  snow_can(land_pts,nsurft),                                                  &
    ! Canopy snow load (kg/m2).
  tsnow0(land_pts)
    ! Temperature of fresh snow (K).

! Snow quantities on a single surft, i.e. over snow layers (sl)
REAL(KIND=real_jlslsm) ::                                                     &
  ds_sl(land_pts,nsmax),                                                      &
    ! Snow layer thicknesses (m).
  sice_sl(land_pts,nsmax),                                                    &
    ! Ice content of snow layers (kg/m2).
  sliq_sl(land_pts,nsmax),                                                    &
    ! Liquid content of snow layers (kg/m2).
  tsnow_sl(land_pts,nsmax),                                                   &
    ! Temperature of fresh snow (K).
  rho_snow_sl(land_pts,nsmax),                                                &
    ! Snow layer densities (kg/m3).
  rgrainl_sl(land_pts,nsmax)
    ! Snow layer grain size (microns).

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  snow_surft_old(:,:),                                                        &
  snow_grnd_old(:,:),                                                         &
  nsnow_old(:,:),                                                             &
  sice_old(:,:),                                                              &
  sliq_old(:,:)
    ! Temporary copies of things that may be needed for diagnostics.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SNOW'

!-----------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Initialise gridbox variables.
!-----------------------------------------------------------------------------

IF (sf_diag%l_snice) THEN
  ALLOCATE( snow_surft_old(land_pts,nsurft) )
  ALLOCATE( snow_grnd_old(land_pts,nsurft)  )
  ALLOCATE( nsnow_old(land_pts,nsurft)      )
  ALLOCATE( sice_old(land_pts,nsurft)       )
  ALLOCATE( sliq_old(land_pts,nsurft)       )

!$OMP PARALLEL                                                                &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(n,k,i)                                                          &
!$OMP SHARED(nsurft,surft_pts,land_pts,sf_diag,snow_surft_old,snow_surft,     &
!$OMP        snow_grnd_old,snow_grnd,sice_old,sliq_old,nsmax,nsnow,           &
!$OMP        sice,sliq,surft_index)

  ! Zeroing diagnostics
  DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
    DO i = 1, land_pts
      sf_diag%snice_m_surft(i,n)        = 0.0
      sf_diag%snice_freez_surft(i,n)    = 0.0
      sf_diag%snice_sicerate_surft(i,n) = 0.0
      sf_diag%snice_sliqrate_surft(i,n) = 0.0
      sf_diag%snice_runoff_surft(i,n)   = 0.0
      sf_diag%snice_smb_surft(i,n)      = 0.0
    END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
    DO k = 1,surft_pts(n)
      i = surft_index(k,n)
      snow_surft_old(i,n)               = snow_surft(i,n)
      snow_grnd_old(i,n)                = snow_grnd(i,n)
      sice_old(i,n)                     = 0.0
      sliq_old(i,n)                     = 0.0
    END DO
!$OMP END DO NOWAIT
  END DO

  IF (nsmax > 0) THEN
    DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
      DO k = 1,surft_pts(n)
        i = surft_index(k,n)
        IF (nsnow(i,n) > 0) THEN
          sice_old(i,n) = SUM(sice(i,n,1:nsnow(i,n)))
          sliq_old(i,n) = SUM(sliq(i,n,1:nsnow(i,n)))
        END IF
      END DO
!$OMP END DO
    END DO
  END IF

!$OMP END PARALLEL

END IF

!$OMP PARALLEL                                                                &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(i)                                                              &
!$OMP SHARED(land_pts,lying_snow,snomlt_sub_htf,snow_melt,infil_rate_con_gbm, &
!$OMP        infil_rate_ls_gbm,dhf_surf_minus_soil,graupel_options,           &
!$OMP        ls_graup,ls_snow)

!$OMP DO SCHEDULE(STATIC)
DO i = 1, land_pts
  lying_snow(i)         = 0.0
  snomlt_sub_htf(i)     = 0.0
  snow_melt(i)          = 0.0
  infil_rate_con_gbm(i) = 0.0
  infil_rate_ls_gbm(i)  = 0.0

  ! initialise FLake tile flux divergence
  dhf_surf_minus_soil(i) = 0.0
END DO
!$OMP END DO

! Set required snow fall
IF (graupel_options == graupel_as_snow) THEN
  ! leave graupel to be treated as snow (ls_snow is total frozen precip)
!$OMP DO SCHEDULE(STATIC)
  DO i = 1, land_pts
    ls_graup(i) = 0.0
  END DO
!$OMP END DO
ELSE IF (graupel_options == ignore_graupel) THEN
  ! get rid of graupel from ls_snow and throw it away completely
!$OMP DO SCHEDULE(STATIC)
  DO i = 1, land_pts
    ls_snow(i)  = MAX(0.0, ls_snow(i) - ls_graup(i) )
    ls_graup(i) = 0.0
  END DO
!$OMP END DO
END IF

!$OMP END PARALLEL

DO n = 1,nsurft
  !---------------------------------------------------------------------------
  ! Set snow mass variables
  !---------------------------------------------------------------------------

  IF ( cansnowtile(n) ) THEN

    !-------------------------------------------------------------------------
    ! With the canopy snow model, snow_surft is held on the canopy,
    ! while the mass of the snowpack (on the ground) is snow_grnd.
    !-------------------------------------------------------------------------
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i)                                                            &
!$OMP SHARED(surft_pts,surft_index,snow_can,snow_surft,snowmass,snow_grnd,    &
!$OMP        ei_surft,melt_surft,timestep,n)
    DO k = 1,surft_pts(n)
      i = surft_index(k,n)
      snow_can(i,n) = snow_surft(i,n)
      snowmass(i)   = snow_grnd(i,n)

      !-----------------------------------------------------------------------
      ! Subtract sublimation and melt of canopy snow.
      !-----------------------------------------------------------------------
      snow_can(i,n) = snow_can(i,n) -                                         &
                     ( ei_surft(i,n) + melt_surft(i,n) ) * timestep

    END DO
!$OMP END PARALLEL DO

  ELSE

    !-------------------------------------------------------------------------
    ! Without the snow canopy model, all the snow is in the snowpack
    !-------------------------------------------------------------------------
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i)                                                            &
!$OMP SHARED(surft_pts,surft_index,snowmass,snow_surft,n)
    DO k = 1,surft_pts(n)
      i = surft_index(k,n)
      snowmass(i) = snow_surft(i,n)
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Copy data for this surft to the snow layer (sl) arrays
!$OMP PARALLEL DEFAULT(NONE)                                                  &
!$OMP PRIVATE(ns,i)                                                           &
!$OMP SHARED(nsmax,land_pts,sice,sice_sl,sliq,sliq_sl,tsnow,tsnow_sl,         &
!$OMP        l_snow_albedo,l_embedded_snow,rgrainl,rgrainl_sl,n)
  DO ns = 1, nsmax
!$OMP DO SCHEDULE(STATIC)
    DO i = 1, land_pts
      sice_sl(i,ns) = sice(i,n,ns)
      sliq_sl(i,ns) = sliq(i,n,ns)
      tsnow_sl(i,ns) = tsnow(i,n,ns)
    END DO
!$OMP END DO NOWAIT
  END DO

  ! Only need rgrainl when certain science options are turned on
  IF (l_snow_albedo .OR. l_embedded_snow) THEN
    DO ns = 1, nsmax
!$OMP DO SCHEDULE(STATIC)
      DO i = 1, land_pts
        rgrainl_sl(i,ns) = rgrainl(i,n,ns)
      END DO
!$OMP END DO NOWAIT
    END DO
  END IF

!$OMP END PARALLEL

  !---------------------------------------------------------------------------
  ! Canopy interception, throughfall and unloading of snow
  !---------------------------------------------------------------------------
  ! The canopy model applies only on PFTs, so set a dummy pointer
  ! in the array on unused tiles.
  IF (cansnowtile(n)) THEN
    n_can = n
  ELSE
    n_can = 1
  END IF

  CALL canopysnow ( land_pts, surft_pts(n), timestep, cansnowtile(n),         &
                    surft_index(:,n), catch_snow(:,n), con_snow,              &
                    ls_snow, ls_graup, unload_backgrnd_pft(:,n_can),          &
                    melt_surft(:,n), snow_can(:,n), snowfall, graupfall )

  !---------------------------------------------------------------------------
  ! Divide snow pack into layers
  !---------------------------------------------------------------------------

  CALL layersnow ( land_pts, surft_pts(n), surft_index(:, n),                 &
                   snowdepth(:,n), nsnow(:,n), ds_sl )

  !---------------------------------------------------------------------------
  ! Thermal properties of snow layers
  !---------------------------------------------------------------------------
  IF ( nsmax > 0 ) THEN
    CALL snowtherm ( land_pts, surft_pts(n), nsnow(:,n),                      &
                     surft_index(:,n), ds_sl, sice_sl,                        &
                     sliq_sl, csnow, ksnow )
  END IF

  !---------------------------------------------------------------------------
  ! Snow thermodynamics and hydrology
  !---------------------------------------------------------------------------
  IF (l_snow_infilt) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i)                                                            &
!$OMP SHARED(surft_pts,surft_index,nsnow,infiltration,ls_rain,con_rain,       &
!$OMP        timestep,infil_rate_con_gbm,tile_frac,infil_rate_ls_gbm,n)
    DO k = 1,surft_pts(n)
      i = surft_index(k,n)

      ! Where there is snow on the ground direct rainfall into infiltration.
      IF (nsnow(i,n) > 0) THEN
        infiltration(i)       = ( ls_rain(i) + con_rain(i) ) * timestep
        infil_rate_con_gbm(i) = infil_rate_con_gbm(i) +                       &
                                tile_frac(i,n) * con_rain(i)
        infil_rate_ls_gbm(i)  = infil_rate_ls_gbm(i) +                        &
                                tile_frac(i,n) * ls_rain(i)
      ELSE
        infiltration(i) = 0.0
      END IF
    END DO
!$OMP END PARALLEL DO
  ELSE
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i)                                                            &
!$OMP SHARED(surft_pts,surft_index,infiltration,n)
    DO k = 1,surft_pts(n)
      i = surft_index(k,n)
      infiltration(i) = 0.0
    END DO
!$OMP END PARALLEL DO
  END IF

  !==========================================================================
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
  ! This means that a soilt variable being passed 'up' to the surface is
  ! simply copied into the surft variable
  !
  ! This will need to be refactored for other tiling approaches. This note
  ! will be replicated elsewhere in the code as required
  !
  !These comments apply until **END NOTICE REGARDING SOIL TILING**
  !===========================================================================

  !Set the current soil tile (see notice above).
  IF (nsoilt == 1) THEN
    !There is only 1 soil tile
    m = 1
  ELSE ! nsoilt == nsurft
    !Soil tiles map directly on to surface tiles
    m = n
  END IF !nsoilt

  CALL snowpack ( n, land_pts, surft_pts(n), timestep, cansnowtile(n),        &
                  nsnow(:,n), surft_index(:,n), csnow, ei_surft(:,n),         &
                  hcaps1_soilt(:,m), hcons, infiltration,                     &
                  ksnow, rho_snow_grnd(:,n), smcl1_soilt(:,m),                &
                  snowfall, sthf1_soilt(:,m), surf_htf_surft(:,n),            &
                  tile_frac(:,n), smvcst1_soilt(:,m), ds_sl,                  &
                  melt_surft(:,n), sice_sl, sliq_sl, snomlt_sub_htf,          &
                  snowdepth(:,n), snowmass, tsnow_sl,                         &
                  t_soil1_soilt(:,m), tsurf_elev_surft(:,n),                  &
                  snow_soil_htf(:,n), rho_snow_sl, rho0, sice0, tsnow0,       &
                  sf_diag,                                                    &
                  !Ancil info (IN)
                  l_lice_point)

  !===========================================================================
  ! *END NOTICE REGARDING SOIL TILING**
  !===========================================================================

  !---------------------------------------------------------------------------
  ! Growth of snow grains
  !---------------------------------------------------------------------------
  IF ( l_snow_albedo .OR. l_embedded_snow ) THEN
    CALL snowgrain ( land_pts, surft_pts(n), timestep, nsnow(:,n),            &
                     surft_index(:,n), sice_sl, snowfall,                     &
                     snowmass, tsnow_sl, tstar_surft(:,n),                    &
                     rgrain(:,n), rgrainl_sl, rgrain0 )
  ELSE
    ! Default initialization required for bit-comparison in the UM.
    rgrain0(:) = r0
  END IF

  IF ( nsmax > 0 ) THEN
    !-------------------------------------------------------------------------
    ! Mechanical compaction of snow
    !-------------------------------------------------------------------------
    CALL compactsnow ( land_pts, surft_pts(n), timestep, nsnow(:,n),          &
                       surft_index(:,n), sice_sl, sliq_sl,                    &
                       tsnow_sl, rho_snow_sl, ds_sl )

    !-------------------------------------------------------------------------
    ! Redivide snowpack after changes in depth, conserving mass and energy
    !-------------------------------------------------------------------------
    CALL relayersnow ( land_pts, surft_pts(n), surft_index(:,n),              &
                       rgrain0, rho0, sice0, snowfall,                        &
                       snowmass, tsnow0, nsnow(:,n), ds_sl,                   &
                       rgrain(:,n), rgrainl_sl, sice_sl,                      &
                       rho_snow_grnd(:,n), sliq_sl,                           &
                       tsnow_sl, rho_snow_sl,                                 &
                       snowdepth(:,n) )
  END IF  !  NSMAX>0

  !---------------------------------------------------------------------------
  ! Copy into final snow mass variables
  !---------------------------------------------------------------------------

!$OMP PARALLEL                                                                &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,ns)                                                         &
!$OMP SHARED(n,cansnowtile,surft_pts,surft_index,snow_grnd,snowmass,          &
!$OMP        snow_surft,snow_can,lying_snow,tile_frac,snow_melt,melt_surft,   &
!$OMP        nsmax,land_pts,ds,ds_sl,sice,sice_sl,sliq,sliq_sl,tsnow,         &
!$OMP        tsnow_sl,rho_snow,rho_snow_sl,l_snow_albedo,                     &
!$OMP        l_embedded_snow,rgrainl,rgrainl_sl)

  IF ( cansnowtile(n) ) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = 1,surft_pts(n)
      i = surft_index(k,n)
      snow_grnd(i,n)  = snowmass(i)
      snow_surft(i,n) = snow_can(i,n)
    END DO
!$OMP END DO
  ELSE
!$OMP DO SCHEDULE(STATIC)
    DO k = 1,surft_pts(n)
      i = surft_index(k,n)
      snow_surft(i,n) = snowmass(i)
    END DO
!$OMP END DO
  END IF

  ! Copy data for this surft from the snow layer (sl) arrays
  DO ns = 1, nsmax
!$OMP DO SCHEDULE(STATIC)
    DO i = 1, land_pts
      ds(i,n,ns) = ds_sl(i,ns)
      sice(i,n,ns) = sice_sl(i,ns)
      sliq(i,n,ns) = sliq_sl(i,ns)
      tsnow(i,n,ns) = tsnow_sl(i,ns)
      rho_snow(i,n,ns) = rho_snow_sl(i,ns)
    END DO
!$OMP END DO NOWAIT
  END DO

  IF (l_snow_albedo .OR. l_embedded_snow) THEN
    DO ns = 1, nsmax
!$OMP DO SCHEDULE(STATIC)
      DO i = 1, land_pts
        rgrainl(i,n,ns) = rgrainl_sl(i,ns)
      END DO
!$OMP END DO NOWAIT
    END DO
  END IF

  !---------------------------------------------------------------------------
  ! Increment gridbox lying snow and snow melt.
  !---------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
  DO k = 1,surft_pts(n)
    i = surft_index(k,n)
    lying_snow(i) = lying_snow(i) +                                           &
                    tile_frac(i,n) * snow_surft(i,n)
    ! Add snow beneath canopy.
    IF ( cansnowtile(n) )                                                     &
      lying_snow(i) = lying_snow(i) +                                         &
                      tile_frac(i,n) * snow_grnd(i,n)

    ! Snow melt.
    snow_melt(i) = snow_melt(i) + tile_frac(i,n) * melt_surft(i,n)
  END DO
!$OMP END DO
!$OMP END PARALLEL

END DO  !  tiles

!$OMP PARALLEL                                                                &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,n,j)                                                        &
!$OMP SHARED(stf_hf_snow_melt,land_pts,hf_snow_melt,snow_melt,surf_ht_flux_ld,&
!$OMP        nsurft,surft_pts,surft_index,tile_frac,snow_soil_htf,lake,       &
!$OMP        dhf_surf_minus_soil,surf_htf_surft)

!-----------------------------------------------------------------------------
! Calculate the total snowmelt heat flux.
!-----------------------------------------------------------------------------
IF ( stf_hf_snow_melt ) THEN
!$OMP DO SCHEDULE(STATIC)
  DO i = 1,land_pts
    hf_snow_melt(i) = lf * snow_melt(i)
  END DO
!$OMP END DO
END IF

!-----------------------------------------------------------------------------
! Calculate gridbox surface heat flux over land.
!-----------------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
DO i = 1,land_pts
  surf_ht_flux_ld(i) = 0.0
END DO
!$OMP END DO

DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
  DO j = 1,surft_pts(n)
    i = surft_index(j,n)
    surf_ht_flux_ld(i) = surf_ht_flux_ld(i) +                                 &
                         tile_frac(i,n) * snow_soil_htf(i,n)
    IF (n == lake) THEN
      dhf_surf_minus_soil(i) = surf_htf_surft(i,n) - snow_soil_htf(i,n)
    END IF
  END DO
!$OMP END DO
END DO

!$OMP END PARALLEL

!-----------------------------------------------------------------------------
! Calculate surface mass balance diagnostics
! Large packs can have a high mass (kg) with small rates of change (/s)
! - beware floating-point precision rounding errors
!-----------------------------------------------------------------------------

IF (sf_diag%l_snice) THEN

!$OMP PARALLEL                                                                &
!$OMP DEFAULT(SHARED)                                                         &
!$OMP PRIVATE(n,k,i)

  DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
    DO k = 1,surft_pts(n)
      i = surft_index(k,n)
      sf_diag%snice_smb_surft(i,n) = (snow_surft(i,n) - snow_surft_old(i,n))  &
                                     / timestep
    END DO
!$OMP END DO
  END DO

  DO n = 1,nsurft
    IF (cansnowtile(n)) THEN
!$OMP DO SCHEDULE(STATIC)
      DO k = 1,surft_pts(n)
        i = surft_index(k,n)
        sf_diag%snice_smb_surft(i,n) = sf_diag%snice_smb_surft(i,n)           &
                                   + (snow_grnd(i,n) - snow_grnd_old(i,n))    &
                                   / timestep
      END DO
!$OMP END DO
    END IF
  END DO

  DO n = 1,nsurft
!$OMP DO SCHEDULE(STATIC)
    DO k = 1,surft_pts(n)
      i = surft_index(k,n)
      IF ( snow_surft(i,n)     > EPSILON(snow_surft) .OR.                     &
           snow_surft_old(i,n) > EPSILON(snow_surft) ) THEN
        IF (nsnow(i,n) > 0) THEN
          sf_diag%snice_sicerate_surft(i,n) = ( SUM(sice(i,n,1:nsnow(i,n))) - &
                                              sice_old(i,n) ) / timestep
          sf_diag%snice_sliqrate_surft(i,n) = ( SUM(sliq(i,n,1:nsnow(i,n))) - &
                                              sliq_old(i,n) ) / timestep
          sf_diag%snice_runoff_surft(i,n) = con_rain(i) + ls_rain(i)          &
                                            + sf_diag%snice_m_surft(i,n)      &
                                            - sf_diag%snice_freez_surft(i,n)  &
                                            - sf_diag%snice_sliqrate_surft(i,n)
        ELSE
          sf_diag%snice_runoff_surft(i,n)   = con_rain(i) + ls_rain(i)        &
                                              + melt_surft(i,n)
          sf_diag%snice_sicerate_surft(i,n) = sf_diag%snice_smb_surft(i,n)
          sf_diag%snice_sliqrate_surft(i,n) = (0.0 - sliq_old(i,n)) / timestep
        END IF
      END IF
    END DO
!$OMP END DO
  END DO
!$OMP END PARALLEL

  DEALLOCATE(sliq_old)
  DEALLOCATE(sice_old)
  DEALLOCATE(nsnow_old)
  DEALLOCATE(snow_grnd_old)
  DEALLOCATE(snow_surft_old)

END IF

!-----------------------------------------------------------------------------
! Remove the rainfall redirected into precipitation.
!-----------------------------------------------------------------------------
IF (l_snow_infilt) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(i)                                                              &
!$OMP SHARED(land_pts,con_rain,infil_rate_con_gbm,ls_rain,infil_rate_ls_gbm)
  DO i = 1,land_pts
    con_rain(i) = con_rain(i) - infil_rate_con_gbm(i)
    ls_rain(i)  = ls_rain(i)  - infil_rate_ls_gbm(i)
  END DO
!$OMP END PARALLEL DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE snow
END MODULE snow_mod
