! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE SURF_HYD-----------------------------------------------

!  PURPOSE : TO CARRY OUT CANOPY AND SURFACE HYDROLOGY CALCULATIONS

!            CANOPY WATER CONTENT IS DEPRECIATED BY EVAPORATION

!            SNOWMELT IS RUNOFF THE SURFACE WITHOUT INTERACTING
!            WITH THE CANOPY

!            THE CANOPY INTERCEPTION AND SURFACE RUNOFF OF
!            LARGE-SCALE RAIN IS CALCUALTED

!            THE CANOPY INTERCEPTION AND SURFACE RUNOFF OF
!            CONVECTIVE RAIN IS CALCUALTED


!  SUITABLE FOR SINGLE COLUMN MODEL USE

!  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER NO 25

!-----------------------------------------------------------------------------

MODULE surf_hyd_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SURF_HYD_MOD'

CONTAINS

SUBROUTINE surf_hyd (                                                         &
            land_pts, nsurft, surft_pts, surft_index,                         &
            catch_surft, ecan_surft, tile_frac, infil_surft, con_rain_land,   &
            ls_rain_land, con_rainfrac_land, ls_rainfrac_land, melt_surft,    &
            snow_melt, timestep,                                              &
            canopy_surft, canopy_gb, dsmc_dt_soilt,                           &
            l_top, l_pdm, sm_levels, soil_pts, soil_index,                    &
            surf_roff_gb, tot_tfall_gb,                                       &
            dun_roff_soilt, fsat_soilt, smvcst_soilt, sthu_soilt, sthf_soilt, &
            surf_roff_soilt,                                                  &
            ! New arguments to replace USE statements
            ! pdm_vars
            slope_gb)


! Use in relevant subroutines
USE frunoff_mod,  ONLY: frunoff
USE sieve_mod,   ONLY: sieve
USE pdm_mod,     ONLY: pdm

! Use in relevant variables
USE jules_surface_mod, ONLY: l_point_data
! Switch for using point rainfall data

USE ancil_info,        ONLY: nsoilt

USE parkind1, ONLY: jprb, jpim
USE yomhook,  ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(IN):
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  land_pts,                                                                   &
    ! Total number of land points.
  nsurft,                                                                     &
    ! Number of tiles.
  sm_levels,                                                                  &
    ! Number of soil moisture levels.
  soil_pts,                                                                   &
    ! Number of soil points.
  surft_pts(nsurft),                                                          &
    ! Number of tile points.
  surft_index(land_pts,nsurft),                                               &
    ! Index of tile points.
  soil_index(land_pts)
    ! Array of soil points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  catch_surft(land_pts,nsurft),                                               &
    ! Canopy capacity for land tiles (kg/m2).
  ecan_surft(land_pts,nsurft),                                                &
    ! Canopy evaporation (kg/m2/s).
  tile_frac(land_pts,nsurft),                                                 &
    ! Tile fractions.
  infil_surft(land_pts,nsurft),                                               &
    ! Tile infiltration rate (kg/m2/s).
  con_rain_land(land_pts),                                                    &
    ! Convective rain (kg/m2/s).
  ls_rain_land(land_pts),                                                     &
    ! Large-scale rain (kg/m2/s).
  con_rainfrac_land(land_pts),                                                &
    ! Convective rain fraction
  ls_rainfrac_land(land_pts),                                                 &
    ! large scale rain fraction
  melt_surft(land_pts,nsurft),                                                &
    ! Snow melt on tiles (kg/m2/s).
  snow_melt(land_pts),                                                        &
    ! GBM snow melt (kg/m2/s).
  fsat_soilt(land_pts,nsoilt),                                                &
    ! Surface saturation fraction.
  smvcst_soilt(land_pts,nsoilt,sm_levels),                                    &
    ! Volumetric soil moisture concentration at saturation (m3 H2O/m3 soil).
  timestep,                                                                   &
    ! Timestep (s).
  sthf_soilt(land_pts,nsoilt,sm_levels),                                      &
    ! Frozen soil moisture content of each layer as a fraction of saturation.
  sthu_soilt(land_pts,nsoilt,sm_levels)
    ! Unfrozen soil moisture content of each layer as a fraction of
    ! saturation.


LOGICAL, INTENT(IN) ::                                                        &
  l_top,                                                                      &
    ! TOPMODEL-based hydrology logical.
  l_pdm
    ! PDM logical.

!-----------------------------------------------------------------------------
! Arguments with INTENT(INOUT):
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                     &
  canopy_surft(land_pts,nsurft)
    ! Tile canopy water contents (kg/m2).

!-----------------------------------------------------------------------------
! Arguments with INTENT(OUT):
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  canopy_gb(land_pts),                                                        &
    ! Gridbox canopy water content (kg/m2).
  dsmc_dt_soilt(land_pts,nsoilt),                                             &
    ! Rate of change of soil moisture content (kg/m2/s).
  surf_roff_gb(land_pts),                                                     &
    ! Cumulative surface runoff (kg/m2/s).
  tot_tfall_gb(land_pts),                                                     &
    ! Cumulative canopy throughfall (kg/m2/s).
  dun_roff_soilt(land_pts,nsoilt),                                            &
    ! Cumulative Dunne sfc runoff (kg/m2/s).
  surf_roff_soilt(land_pts,nsoilt)
    ! Soil-tiled contributions to surface runoff (kg/m2/s).

!-----------------------------------------------------------------------------
! New arguments to replace USE statements:
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN) :: slope_gb(land_pts)

!-----------------------------------------------------------------------------
! Local variables:
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  i,j,n,m
    ! Loop counters:
    ! i for land point
    ! j for tile point
    ! n for surface tile
    ! m for soil tile

REAL(KIND=real_jlslsm) ::                                                     &
  r,                                                                          &
    ! Total downward water flux (i.e. rain + condensation + snowmelt)
    ! (kg/m2/s).
  tfall,                                                                      &
    ! Cumulative canopy throughfall on a tile (kg/m2/s).
  s_roff,                                                                     &
    ! Cumulative surface runoff on a tile (kg/m2/s).
  can_cond(land_pts),                                                         &
    ! Canopy condensation (kg/m2/s).
  frac_cov(land_pts),                                                         &
    ! Fraction of gridbox with snowmelt or evap.
  surf_roff_surft(land_pts,nsurft),                                           &
    ! Surface-tiled contributions to surface runoff.
  tot_tfall_surft(land_pts,nsurft),                                           &
    ! Surface-tiled contributions to throughfall.
  tot_tfall_soilt(land_pts,nsoilt)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SURF_HYD'

!-----------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Assume snowmelt and evaporation cover 100% of the gridbox
frac_cov(:) = 1.0

! Zero cumulative stores
DO i = 1,land_pts
  tot_tfall_gb(i)     = 0.0
  surf_roff_gb(i)     = 0.0
  dsmc_dt_soilt(i,:)  = 0.0
  dun_roff_soilt(i,:) = 0.0
END DO
surf_roff_surft(:,:) = 0.0
tot_tfall_surft(:,:) = 0.0

! Reduce canopy water content by evaporation
DO n = 1,nsurft
  DO j = 1,surft_pts(n)
    i = surft_index(j,n)
    IF (ecan_surft(i,n) > 0.0) THEN
      canopy_surft(i,n) =  MAX(canopy_surft(i,n) - ecan_surft(i,n) * timestep,&
                               0.0 )
    END IF
  END DO
END DO

IF (l_point_data) THEN
  !-----------------------------------------------------------------------
  ! Using point precipitation data.
  !-----------------------------------------------------------------------
  DO n = 1,nsurft
    DO j = 1,surft_pts(n)
      i = surft_index(j,n)
      ! Calculate total downward flux of liquid water.
      r = ls_rain_land(i) + con_rain_land(i)
      ! Add any condensation
      IF ( ecan_surft(i,n) < 0.0 ) r = r - ecan_surft(i,n)
      ! Calculate throughfall.
      IF ( r <= catch_surft(i,n) / timestep .AND.                             &
           catch_surft(i,n) > 0 ) THEN
        tfall = r * canopy_surft(i,n) / catch_surft(i,n)
      ELSE
        tfall = r - ( catch_surft(i,n) - canopy_surft(i,n) ) / timestep
      END IF
      ! Update canopy water content.
      canopy_surft(i,n) = canopy_surft(i,n) + ( r - tfall ) * timestep
      ! Add melt to throughfall before calculation of runoff.
      tfall = tfall + melt_surft(i,n)
      ! Calculate surface runoff.
      IF ( tfall > infil_surft(i,n) ) THEN
        s_roff = tfall - infil_surft(i,n)
      ELSE
        s_roff = 0.0
      END IF

      ! Add to gridbox accumulations.
      ! Don't include melt in throughfall.
      tot_tfall_gb(i) = tot_tfall_gb(i)                                       &
                        + ( tfall - melt_surft(i,n) ) * tile_frac(i,n)
      surf_roff_gb(i) = surf_roff_gb(i) + s_roff * tile_frac(i,n)

      ! Tiled runoff and throughfall.
      tot_tfall_surft(i,n) = tfall - melt_surft(i,n)
      surf_roff_surft(i,n) = s_roff

    END DO
  END DO

ELSE

  !-----------------------------------------------------------------------
  ! Using area-average precipitation data. Assume spatial distribution.
  !-----------------------------------------------------------------------

  DO n = 1,nsurft

    ! Surface runoff of snowmelt, assumed to cover 100% of tile
    CALL frunoff (land_pts, surft_pts(n), surft_index(:,n), frac_cov,         &
                  catch_surft(:,n), catch_surft(:,n), infil_surft(:,n),       &
                  melt_surft(:,n), tile_frac(:,n), timestep,                  &
                  surf_roff_gb, surf_roff_surft(:,n))

    ! Define canopy condensation when evaporation is negative
    DO j = 1,surft_pts(n)
      i = surft_index(j,n)
      IF ( ecan_surft(i,n)  <  0.0 ) THEN
        can_cond(i) = - ecan_surft(i,n)
      ELSE
        can_cond(i) = 0.0
      END IF
    END DO

    ! Canopy interception, throughfall and surface runoff for condensation,
    ! assumed to cover 100% of gridbox
    CALL sieve (land_pts, surft_pts(n), surft_index(:,n),frac_cov,            &
                catch_surft(:,n), can_cond, tile_frac(:,n), timestep,         &
                canopy_surft(:,n), tot_tfall_gb, tot_tfall_surft(:,n))

    CALL frunoff (land_pts, surft_pts(n), surft_index(:,n), frac_cov,         &
                  catch_surft(:,n), canopy_surft(:,n), infil_surft(:,n),      &
                  can_cond, tile_frac(:,n), timestep,                         &
                  surf_roff_gb, surf_roff_surft(:,n))

    ! Canopy interception, throughfall and surface runoff for large-scale
    ! rain, assumed to cover "ls_rainfrac_land" of gridbox
    CALL sieve (land_pts, surft_pts(n), surft_index(:,n), ls_rainfrac_land,   &
                catch_surft(:,n), ls_rain_land, tile_frac(:,n), timestep,     &
                canopy_surft(:,n), tot_tfall_gb, tot_tfall_surft(:,n))

    CALL frunoff (land_pts, surft_pts(n), surft_index(:,n), ls_rainfrac_land, &
                  catch_surft(:,n), canopy_surft(:,n), infil_surft(:,n),      &
                  ls_rain_land, tile_frac(:,n), timestep,                     &
                  surf_roff_gb, surf_roff_surft(:,n))

    ! Canopy interception, throughfall and surface runoff for convective
    ! rain, assumed to cover fraction "con_rainfrac_land" of gridbox
    CALL sieve (land_pts, surft_pts(n), surft_index(:,n), con_rainfrac_land,  &
                catch_surft(:,n), con_rain_land, tile_frac(:,n), timestep,    &
                canopy_surft(:,n), tot_tfall_gb, tot_tfall_surft(:,n))

    CALL frunoff (land_pts, surft_pts(n), surft_index(:,n), con_rainfrac_land,&
                  catch_surft(:,n), canopy_surft(:,n), infil_surft(:,n),      &
                  con_rain_land, tile_frac(:,n), timestep,                    &
                  surf_roff_gb, surf_roff_surft(:,n))

  END DO
END IF   !   l_point_data

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

IF ( nsoilt == 1) THEN
  ! There is only one soil tile
  m = 1

  ! To maintain bit-comparability, we need to run pdm with the _gb versions.
  !---------------------------------------------------------------------------
  ! Calculate Saturation excess runoff through PDM:
  !---------------------------------------------------------------------------
  IF (soil_pts > 0 .AND. l_pdm) THEN
    CALL pdm(                                                                 &
            land_pts, soil_pts, soil_index, sm_levels,                        &
            tot_tfall_gb, snow_melt, surf_roff_gb, timestep,                  &
            smvcst_soilt(:,m,:), dun_roff_soilt(:,m),                         &
            sthu_soilt(:,m,:), sthf_soilt(:,m,:),                             &
            ! New arguments to replace USE statements
            ! pdm_vars
            slope_gb)
  END IF

  IF (l_top) THEN
    DO i = 1,land_pts
      dun_roff_soilt(i,m) = fsat_soilt(i,m) * (tot_tfall_gb(i)                &
                            + snow_melt(i) - surf_roff_gb(i))
    END DO
  END IF

  DO i = 1,land_pts
    IF (l_top .OR. l_pdm) THEN
      surf_roff_gb(i) = surf_roff_gb(i) + dun_roff_soilt(i,m)
    END IF
    dsmc_dt_soilt(i,m) = tot_tfall_gb(i) + snow_melt(i) - surf_roff_gb(i)
  END DO

  ! Copy the final answers across to the _soilt version so they can be used
  ! elsewhere.
  surf_roff_soilt(:,m) = surf_roff_gb(:)

ELSE
  ! nsoilt = nsurft
  surf_roff_soilt(:,:)   = surf_roff_surft(:,:)
  tot_tfall_soilt(:,:)   = tot_tfall_surft(:,:)

  !-----------------------------------------------------------------------
  ! Calculate saturation excess runoff through PDM:
  !-----------------------------------------------------------------------
  IF (soil_pts > 0 .AND. l_pdm) THEN
    DO m = 1, nsoilt
      n = m
      !Note we use the surface-tiled melt, rather than the GBM
      CALL pdm(                                                               &
              land_pts, soil_pts, soil_index, sm_levels,                      &
              tot_tfall_soilt(:,m), melt_surft(:,n), surf_roff_soilt(:,m),    &
              timestep, smvcst_soilt(:,m,:), dun_roff_soilt(:,m),             &
              sthu_soilt(:,m,:), sthf_soilt(:,m,:),                           &
              ! New arguments to replace USE statements
              ! pdm_vars
              slope_gb)
    END DO
  END IF

  IF (l_top) THEN
    DO m = 1, nsoilt
      n = m
      DO i = 1,land_pts
        dun_roff_soilt(i,m) = fsat_soilt(i,m) * (tot_tfall_soilt(i,m)         &
                              + melt_surft(i,n) - surf_roff_soilt(i,m))
      END DO
    END DO
  END IF

  DO m = 1, nsoilt
    n = m
    DO i = 1,land_pts
      IF (l_top .OR. l_pdm) THEN
        surf_roff_soilt(i,m) = surf_roff_soilt(i,m) + dun_roff_soilt(i,m)
        surf_roff_gb(i) = surf_roff_gb(i)                                     &
                          + tile_frac(i,n) * dun_roff_soilt(i,m)
      END IF
      dsmc_dt_soilt(i,m) = tot_tfall_soilt(i,m) + melt_surft(i,n)             &
                           - surf_roff_soilt(i,m)
    END DO
  END DO

END IF  !  nsoilt

!=============================================================================
! *END NOTICE REGARDING SOIL TILING**
!=============================================================================

!Calulcate GBM canopy water content
canopy_gb(:) = 0.0
DO n = 1, nsurft
  DO j = 1,surft_pts(n)
    i = surft_index(j,n)
    canopy_gb(i) = canopy_gb(i) + (tile_frac(i,n) * canopy_surft(i,n))
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE surf_hyd
END MODULE surf_hyd_mod
