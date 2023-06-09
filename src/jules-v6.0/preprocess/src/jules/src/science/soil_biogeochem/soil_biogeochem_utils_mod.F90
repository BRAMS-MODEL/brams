!******************************COPYRIGHT**************************************
! (c) Centre for Ecology and Hydrology.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT**************************************

MODULE soil_biogeochem_utils_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Utility routines for soil biogeochemistry.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in BIOGEOCHEMISTRY
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Private scope by default.
PRIVATE
PUBLIC soil_n_precalc

CONTAINS

!#############################################################################

SUBROUTINE soil_n_precalc( land_pts, trif_pts, nstep_trif, trif_index,        &
                           n_inorg_avail_pft, n_inorg_soilt_lyrs,             &
                           fert_layer, have_layers, dz, f_root_pft,           &
                           f_root_pft_dz, isunfrozen,                         &
                        ! These arguments replace USE statements
                           n_soil_pool_soilt, dim_soil_n_pool,                &
                           ! prognostics
                           t_soil_soilt_acc )

!-----------------------------------------------------------------------------
! Description:
!   Calculates soil nitrogen-related variables.
!   This is coded only for a single soil tile.
!-----------------------------------------------------------------------------

USE ancil_info, ONLY:                                                         &
  ! imported scalars
  dim_cslayer, nsoilt

USE conversions_mod, ONLY:                                                    &
  ! imported scalar parameters
  zerodegc

USE ecosse_utils_mod, ONLY:                                                   &
  ! imported procedures
  get_residual_n, transfer_layers

USE ereport_mod, ONLY:                                                        &
  ! imported procedures
  ereport

USE jules_soil_biogeochem_mod, ONLY:                                          &
  ! imported scalar parameters
  soil_model_rothc, soil_model_ecosse,                                        &
  ! imported scalars
  soil_bgc_model

USE jules_soil_ecosse_mod, ONLY:                                              &
  ! imported arrays
  dz_soilc

USE jules_soil_mod, ONLY:                                                     &
  ! imported scalars
  sm_levels,                                                                  &
  ! imported arrays
  dzsoil

USE jules_surface_types_mod, ONLY:                                            &
  ! imported scalars
  npft

USE pftparm, ONLY:                                                            &
  ! imported arrays
  rootd_ft

USE root_frac_mod, ONLY:                                                      &
  ! imported procedures
  root_frac

USE soil_ecosse_vars_mod, ONLY:                                               &
  ! imported parameters
  i_amm, i_nit

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  land_pts,                                                                   &
    ! Number of land points.
  trif_pts,                                                                   &
    ! Number of points on which TRIFFID may operate.
  nstep_trif
    ! Number of atmospheric timesteps between calls to TRIFFID.

!-----------------------------------------------------------------------------
! Array arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  trif_index(land_pts)
    ! Indices of land points on which TRIFFID may operate.

!-----------------------------------------------------------------------------
! Array arguments with intent(inout)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
  n_inorg_avail_pft(land_pts,npft,dim_cslayer),                               &
    ! Available inorganic N for PFTs (kg m-2).
  n_inorg_soilt_lyrs(land_pts,nsoilt,dim_cslayer)
    ! Gridbox inorganic N pool on soil levels (kg m-2).

!-----------------------------------------------------------------------------
! Scalar arguments with intent(out)
!-----------------------------------------------------------------------------
INTEGER, INTENT(OUT) ::                                                       &
  fert_layer
    ! Deepest soil level to which fertiliser is applied.

LOGICAL, INTENT(OUT) ::                                                       &
  have_layers
    !  T if we have multiple soil layers, else F.

!-----------------------------------------------------------------------------
! Array arguments with intent(out)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  dz(dim_cslayer),                                                            &
    ! Thicknesses of soil layers (m).
  f_root_pft(npft,dim_cslayer),                                               &
    ! Root fraction in each soil layer per PFT.
  f_root_pft_dz(npft,dim_cslayer),                                            &
    ! Normalised roots in each soil layer per PFT.
  isunfrozen(land_pts,dim_cslayer)
    ! Unfrozen soil flag. 1 in unfrozen layers, else 0.

!-----------------------------------------------------------------------------
! These arguments replace USE statements
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) :: dim_soil_n_pool
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  n_soil_pool_soilt(land_pts,nsoilt,dim_cslayer,dim_soil_n_pool)
REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
  t_soil_soilt_acc(land_pts,nsoilt,sm_levels)
    ! Sub-surface temperature on layers and soil tiles accumulated over
    ! TRIFFID timestep (K). On entry this is the accumulation, on exit this
    ! is the average temperature.

!-----------------------------------------------------------------------------
! Local scalar parameters.
!-----------------------------------------------------------------------------
INTEGER, PARAMETER ::                                                         &
  s = 1
    ! Soil tile number.

LOGICAL, PARAMETER ::                                                         &
  l_interp_true = .TRUE.,                                                     &
    ! Flag indicating that variables shouuld be interpolated for ECOSSE.
  l_midpoint_true = .TRUE.
    ! Flag indicating that variables apply at mid-point of each layer,
    ! used with interpolation for ECOSSE.

!-----------------------------------------------------------------------------
! Local scalar variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  errorstatus,                                                                &
    ! Value used for error.
  iz, l, n, t  ! Indices.

!-----------------------------------------------------------------------------
! Local array variables.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                     &
  residual_n(land_pts,dim_cslayer),                                           &
    ! Minimum-allowed (residual) inorganic N amount (kg m-2).
  t_soil(land_pts,dim_cslayer)
    ! Soil temperature on soil biogeochemistry layers (K).

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'SOIL_N_PRECALC'

!-----------------------------------------------------------------------------
!end of header

!-----------------------------------------------------------------------------
! Some of this subroutine assumes we have a single soil tile.
!-----------------------------------------------------------------------------
IF ( nsoilt > 1 ) THEN
  errorstatus = 101
  CALL ereport( TRIM(RoutineName), errorstatus,                               &
                "Only coded for a single soil tile." )
END IF

!-----------------------------------------------------------------------------
! Initialisation.
!-----------------------------------------------------------------------------
isunfrozen(:,:) = 1.0

!-----------------------------------------------------------------------------
! Set flag for multiple layers.
!-----------------------------------------------------------------------------
have_layers = dim_cslayer > 1

!-----------------------------------------------------------------------------
! Get soil layer thicknesses.
!-----------------------------------------------------------------------------
SELECT CASE ( soil_bgc_model )
CASE ( soil_model_rothc )
  dz(1:dim_cslayer) = dzsoil(1:dim_cslayer)
CASE ( soil_model_ecosse )
  dz(:) = dz_soilc(:)
END SELECT

IF ( have_layers ) THEN

  !---------------------------------------------------------------------------
  ! Find the layer that contains the 30cm depth, for fertiliser addition.
  !---------------------------------------------------------------------------
  fert_layer = 0
  DO iz = dim_cslayer,1,-1
    IF ( SUM(dz(1:iz)) >= 0.3 ) fert_layer = iz
  END DO
  ! If soil depth < 30cm, add fertiliser to all layers.
  IF ( fert_layer == 0 ) fert_layer = dim_cslayer

  !-------------------------------------------------------------------------
  ! Calculate the fraction and normalised roots in each soil layer.
  !-------------------------------------------------------------------------
  DO n = 1,npft
    CALL root_frac( n, dim_cslayer, dz, rootd_ft(n), f_root_pft(n,:) )
    f_root_pft_dz(n,:) = f_root_pft(n,:) / dzsoil(:) /                        &
                        (f_root_pft(n,1) / dzsoil(1))
  END DO

  !-------------------------------------------------------------------------
  ! Set flag in frozen layers.
  !-------------------------------------------------------------------------
  SELECT CASE( soil_bgc_model )

  CASE ( soil_model_rothc )

    ! Calculate time-averaged soil temperature.
    t_soil_soilt_acc(:,:,:) = t_soil_soilt_acc(:,:,:) / REAL(nstep_trif)

    DO t = 1,trif_pts
      l = trif_index(t)
      DO iz = 1,dim_cslayer
        IF ( t_soil_soilt_acc(l,1,iz) < zerodegc ) isunfrozen(l,iz) = 0.0
      END DO
    END DO

  CASE ( soil_model_ecosse )

    ! Calculate time-averaged soil temperature.
    t_soil_soilt_acc(:,:,:) = t_soil_soilt_acc(:,:,:) / REAL(nstep_trif)

    ! Get time-averaged soil temperature on soil N layers.
    t_soil(:,:) = transfer_layers(                                            &
                              land_pts, trif_pts, sm_levels, dim_cslayer,     &
                              l_midpoint_true, l_interp_true, trif_index,     &
                              dzsoil, dz_soilc,                               &
                              t_soil_soilt_acc(:,1,:) )
    ! Look for unfrozen layers.
    DO t = 1,trif_pts
      l = trif_index(t)
      DO iz = 1,dim_cslayer
        IF ( t_soil(l,iz) < zerodegc ) isunfrozen(l,iz) = 0.0
      END DO
    END DO

  END SELECT

ELSE

  ! .NOT. have_layers
  ! Indicate that fertiliser should be added to the only layer.
  fert_layer = 1

END IF  !  have_layers

!-----------------------------------------------------------------------------
! Load ECOSSE inorganic nitrogen variables.
!-----------------------------------------------------------------------------
IF ( soil_bgc_model == soil_model_ecosse ) THEN

  ! Calculate the minimum-allowed amount of inorganic N.
  CALL get_residual_n( land_pts, trif_pts, trif_index, residual_n )

  ! All PFTs can access the same soil N.
  DO n = 1,npft
    n_inorg_avail_pft(:,n,:) = ( n_soil_pool_soilt(:,s,:,i_amm)               &
                                 - residual_n(:,:) )                          &
                               + ( n_soil_pool_soilt(:,s,:,i_nit)             &
                                  - residual_n(:,:) )
  END DO

  n_inorg_soilt_lyrs(:,s,:)  = n_soil_pool_soilt(:,s,:,i_amm)                 &
                               + n_soil_pool_soilt(:,s,:,i_nit)

END IF  !  soil_model_ecosse

END SUBROUTINE soil_n_precalc

!#############################################################################

END MODULE soil_biogeochem_utils_mod
