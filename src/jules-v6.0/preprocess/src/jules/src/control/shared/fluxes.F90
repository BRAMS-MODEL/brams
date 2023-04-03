! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module containing surface fluxes.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.
!
! Code Owner: Please refer to ModuleLeaders.txt and UM file CodeOwners.txt
! This file belongs in section: Land
!
! Changes to variable names to enable soil tiling
! Anything named _tile is now ambiguous, so the following convention is
! adopted:
! _gb for variables on land points
! _ij for variables with i and j indices
! Surface heterogeneity...
! _surft for surface tiled variables (generally size n_surft)
! _pft for plant functional type surface tiled variables (generally sized n_pft)
! _sicat for sea ice catergories (generally size nice_use)
! Sub-surface heterogeneity...
! _soilt for soil tiled variables

MODULE fluxes

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! anthrop_heat is required by both the UM and standalone configurations
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE :: anthrop_heat_surft(:,:)
! Additional heat source on surface tiles used for anthropgenic urban heat
! source (W/m2)

REAL(KIND=real_jlslsm), ALLOCATABLE :: surf_ht_store_surft(:,:)
! Diagnostic to store values of C*(dT/dt) during calculation of energy
! balance

REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: sw_sicat(:,:)
! Net SW on sea ice categories
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: sw_rts_sicat(:,:)
! Net SW on sea ice categories on the radiative timestep
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: swup_rts_sicat(:,:)
! Upward SW on sea ice categories on the radiative timestep
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: swdn_rts_sicat(:,:)
! Downward SW on sea ice categories on the radiative timestep
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: alb_sicat(:,:,:)
! Albedo of sea ice categories (ij point, sicat, band- see below)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: penabs_rad_frac(:,:,:)
! Fraction of downward solar that penetrates the sea ice and is absorbed
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: penabs_rad_rts(:,:)
! Penetrated-absorbed radiation on sea ice categories on radiation timestep
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: sw_sea(:)
! Net SW on open sea
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: sw_rts_sea(:)
! Net SW on open sea on the radiative timestep
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: t_growth_gb(:)
! Average temperature (growth temperature) for thermal adaptation or
! acclimation of photosynthetic capacity (degrees Celsius).

!-----------------------------------------------------------------------------
! Everything else is required only in standalone configuration
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE :: alb_surft(:,:,:)
!   Albedo for surface tiles
!     (:,:,1) direct beam visible
!     (:,:,2) diffuse visible
!     (:,:,3) direct beam near-IR
!     (:,:,4) diffuse near-IR
REAL(KIND=real_jlslsm), ALLOCATABLE :: e_sea_ij(:,:)
!   Evaporation from sea times leads fraction. Zero over land
!                                (kg per square metre per sec)
REAL(KIND=real_jlslsm), ALLOCATABLE :: ecan_ij(:,:)
!   Gridbox mean evaporation from canopy/surface store (kg/m2/s)
!     Zero over sea
REAL(KIND=real_jlslsm), ALLOCATABLE :: ecan_surft(:,:)
!   Canopy evaporation from for snow-free land tiles
REAL(KIND=real_jlslsm), ALLOCATABLE :: ei_ij(:,:)
!   Sublimation from lying snow or sea-ice (kg/m2/s)
REAL(KIND=real_jlslsm), ALLOCATABLE :: ei_surft(:,:)
!   EI for land tiles
REAL(KIND=real_jlslsm), ALLOCATABLE :: esoil_ij_soilt(:,:,:)
!   Surface evapotranspiration from soil moisture store (kg/m2/s)
REAL(KIND=real_jlslsm), ALLOCATABLE :: esoil_surft(:,:)
!   ESOIL for snow-free land tiles
REAL(KIND=real_jlslsm), ALLOCATABLE :: ext_soilt(:,:,:)
!   Extraction of water from each soil layer (kg/m2/s)
REAL(KIND=real_jlslsm), ALLOCATABLE :: fqw_surft(:,:)
!   Surface FQW for land tiles
REAL(KIND=real_jlslsm), ALLOCATABLE :: fqw_sicat(:,:,:)
!   Surface FQW for sea-ice
REAL(KIND=real_jlslsm), ALLOCATABLE :: fsmc_pft(:,:)
!   Moisture availability factor.
REAL(KIND=real_jlslsm), ALLOCATABLE :: ftl_sicat(:,:,:)
!   Surface FTL for sea-ice
REAL(KIND=real_jlslsm), ALLOCATABLE :: ftl_surft(:,:)
!   Surface FTL for land tiles
REAL(KIND=real_jlslsm), ALLOCATABLE :: h_sea_ij(:,:)
!   Surface sensible heat flux over sea times leads fraction (W/m2)
REAL(KIND=real_jlslsm), ALLOCATABLE :: hf_snow_melt_gb(:)
!   Gridbox snowmelt heat flux (W/m2)
REAL(KIND=real_jlslsm), ALLOCATABLE :: land_albedo_ij(:,:,:)
!   GBM albedo
!     (:,:,1) direct beam visible
!     (:,:,2) diffuse visible
!     (:,:,3) direct beam near-IR
!     (:,:,4) diffuse near-IR
REAL(KIND=real_jlslsm), ALLOCATABLE :: le_surft(:,:)
!   Surface latent heat flux for land tiles
REAL(KIND=real_jlslsm), ALLOCATABLE :: melt_surft(:,:)
!   Snowmelt on land tiles (kg/m2/s)
REAL(KIND=real_jlslsm), ALLOCATABLE :: sea_ice_htf_sicat(:,:,:)
!   Heat flux through sea-ice (W/m2, positive downwards)
REAL(KIND=real_jlslsm), ALLOCATABLE :: snomlt_sub_htf_gb(:)
!   Sub-canopy snowmelt heat flux (W/m2)
REAL(KIND=real_jlslsm), ALLOCATABLE :: snow_melt_gb(:)
!   snowmelt on land points (kg/m2/s)
REAL(KIND=real_jlslsm), ALLOCATABLE :: snowmelt_ij(:,:)
!   Snowmelt (kg/m2/s)
REAL(KIND=real_jlslsm), ALLOCATABLE :: sub_surf_roff_gb(:)
!   Sub-surface runoff (kg/m2/s)
REAL(KIND=real_jlslsm), ALLOCATABLE :: surf_ht_flux_ij(:,:)
!   Net downward heat flux at surface over land and sea-ice fraction of
!gridbox (W/m2)
REAL(KIND=real_jlslsm), ALLOCATABLE ::  snow_soil_htf(:,:)
!   Heat flux under snow to subsurface on tiles (W/m2)
REAL(KIND=real_jlslsm), ALLOCATABLE :: surf_htf_surft(:,:)
!   Surface heat flux on land tiles (W/m2)
REAL(KIND=real_jlslsm), ALLOCATABLE :: surf_roff_gb(:)
!   Surface runoff (kg/m2/s)
REAL(KIND=real_jlslsm), ALLOCATABLE :: radnet_surft(:,:)
!   Surface net radiation on tiles ( W/m2)
REAL(KIND=real_jlslsm), ALLOCATABLE :: tot_tfall_gb(:)
!   Total throughfall (kg/m2/s)
REAL(KIND=real_jlslsm), ALLOCATABLE :: tstar_ij(:,:)
!   GBM surface temperature (K)
REAL(KIND=real_jlslsm), ALLOCATABLE :: emis_surft(:,:)
!   Tile emissivity
REAL(KIND=real_jlslsm), ALLOCATABLE :: sw_surft(:,:)
!   Surface net shortwave on tiles (W/m2)
REAL(KIND=real_jlslsm), ALLOCATABLE :: rflow_gb(:)
!   River outflow on model grid (kg/m2/s)
REAL(KIND=real_jlslsm), ALLOCATABLE :: rrun_gb(:)
!   Runoff after river routing on model grid (kg/m2/s)
REAL(KIND=real_jlslsm), ALLOCATABLE :: z0m_surft(:,:)
!   Tile roughness lengths for momentum.
REAL(KIND=real_jlslsm), ALLOCATABLE :: z0h_surft(:,:)
!   Tile roughness lengths for heat and moisture (m).

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='FLUXES'

CONTAINS
  
SUBROUTINE fluxes_alloc(land_pts, t_i_length, t_j_length,                     &
                  nsurft, npft, nsoilt, sm_levels,                            &
                  nice, nice_use)

USE jules_vegetation_mod, ONLY: photo_adapt, photo_acclim, photo_acclim_model

!No USE statements other than error reporting and Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
INTEGER, INTENT(IN) :: land_pts, t_i_length, t_j_length,                      &
                       nsurft, npft, nsoilt, sm_levels,                       &
                       nice, nice_use
!Local variables

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FLUXES_ALLOC'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!  ====fluxes module, common
ALLOCATE( surf_ht_store_surft(land_pts,nsurft))
ALLOCATE( anthrop_heat_surft(land_pts,nsurft))
ALLOCATE( sw_rts_sicat(t_i_length * t_j_length, nice_use))
ALLOCATE( swup_rts_sicat(t_i_length * t_j_length, nice_use))
ALLOCATE( swdn_rts_sicat(t_i_length * t_j_length, nice_use))
ALLOCATE( sw_sicat(t_i_length * t_j_length, nice_use))
ALLOCATE( alb_sicat(t_i_length * t_j_length, nice_use, 4))
ALLOCATE( penabs_rad_frac(t_i_length * t_j_length, nice_use, 4))
ALLOCATE( penabs_rad_rts(t_i_length * t_j_length, nice_use))
ALLOCATE( sw_rts_sea(t_i_length * t_j_length))
ALLOCATE( sw_sea(t_i_length * t_j_length))
SELECT CASE ( photo_acclim_model )
CASE ( photo_adapt, photo_acclim )
  ALLOCATE( t_growth_gb(land_pts))
  t_growth_gb(:)         = 0.0
CASE DEFAULT
  ALLOCATE( t_growth_gb(1))
  t_growth_gb(1)         = 0.0
END SELECT

surf_ht_store_surft(:,:) = 0.0
anthrop_heat_surft(:,:)  = 0.0
sw_rts_sicat(:,:)        = 0.0
swup_rts_sicat(:,:)      = 0.0
swdn_rts_sicat(:,:)      = 0.0
sw_sicat(:,:)            = 0.0
alb_sicat(:,:,:)         = 0.0
penabs_rad_frac(:,:,:)   = 0.0
penabs_rad_rts(:,:)      = 0.0
sw_rts_sea(:)            = 0.0
sw_sea(:)                = 0.0

!  ====fluxes module JULES====
! Runoff components.
ALLOCATE( sub_surf_roff_gb(land_pts))
ALLOCATE( surf_roff_gb(land_pts))

! (Forcing) fluxes
ALLOCATE( alb_surft(land_pts,nsurft,4) )
ALLOCATE( tstar_ij(t_i_length,t_j_length))
ALLOCATE( e_sea_ij(t_i_length,t_j_length))
ALLOCATE( fsmc_pft(land_pts,npft))
ALLOCATE( ftl_surft(land_pts,nsurft))
ALLOCATE( le_surft(land_pts,nsurft))
ALLOCATE( h_sea_ij(t_i_length,t_j_length))
ALLOCATE( fqw_surft(land_pts,nsurft))
ALLOCATE( fqw_sicat(t_i_length,t_j_length,nice_use))
ALLOCATE( ftl_sicat(t_i_length,t_j_length,nice_use))
ALLOCATE( ecan_ij(t_i_length,t_j_length))
ALLOCATE( esoil_surft(land_pts,nsurft))
ALLOCATE( sea_ice_htf_sicat(t_i_length,t_j_length,nice))
ALLOCATE( surf_ht_flux_ij(t_i_length,t_j_length))
ALLOCATE( surf_htf_surft(land_pts,nsurft))
ALLOCATE( snow_soil_htf(land_pts,nsurft))
ALLOCATE( land_albedo_ij(t_i_length,t_j_length,4) )
ALLOCATE( ei_ij(t_i_length,t_j_length))
ALLOCATE( ei_surft(land_pts,nsurft))
ALLOCATE( ecan_surft(land_pts,nsurft))
ALLOCATE( esoil_ij_soilt(t_i_length,t_j_length,nsoilt))
ALLOCATE( ext_soilt(land_pts,nsoilt,sm_levels))
ALLOCATE( snowmelt_ij(t_i_length,t_j_length))
ALLOCATE( hf_snow_melt_gb(land_pts))
ALLOCATE( radnet_surft(land_pts,nsurft) )
ALLOCATE( sw_surft(land_pts,nsurft) )
ALLOCATE( emis_surft(land_pts,nsurft) )
ALLOCATE( snow_melt_gb(land_pts))
ALLOCATE( snomlt_sub_htf_gb(land_pts))
ALLOCATE( tot_tfall_gb(land_pts))
ALLOCATE( melt_surft(land_pts,nsurft))
ALLOCATE( z0m_surft(land_pts,nsurft))
ALLOCATE( z0h_surft(land_pts,nsurft))
ALLOCATE( rflow_gb(land_pts))
ALLOCATE( rrun_gb(land_pts))

sub_surf_roff_gb(:)          = 0.0
surf_roff_gb(:)              = 0.0
alb_surft(:,:,:)             = 0.0
tstar_ij(:,:)                = 0.0
e_sea_ij(:,:)                = 0.0
fsmc_pft(:,:)                = 0.0
ftl_surft(:,:)               = 0.0
le_surft(:,:)                = 0.0
h_sea_ij(:,:)                = 0.0
fqw_surft(:,:)               = 0.0
fqw_sicat(:,:,:)             = 0.0
ftl_sicat(:,:,:)             = 0.0
ecan_ij(:,:)                 = 0.0
esoil_surft(:,:)             = 273.15
sea_ice_htf_sicat(:,:,:)     = 0.0
surf_ht_flux_ij(:,:)         = 0.0
surf_htf_surft(:,:)          = 0.0
snow_soil_htf(:,:)           = 0.0
land_albedo_ij(:,:,:)        = 0.0
ei_ij(:,:)                   = 0.0
ei_surft(:,:)                = 0.0
ecan_surft(:,:)              = 0.0
esoil_ij_soilt(:,:,:)        = 0.0
ext_soilt(:,:,:)             = 0.0
snowmelt_ij(:,:)             = 0.0
hf_snow_melt_gb(:)           = 0.0
radnet_surft(:,:)            = 0.0
sw_surft(:,:)                = 0.0
emis_surft(:,:)              = 0.0
snow_melt_gb(:)              = 0.0
snomlt_sub_htf_gb(:)         = 0.0
tot_tfall_gb(:)              = 0.0
melt_surft(:,:)              = 0.0
z0m_surft(:,:)               = 0.0
z0h_surft(:,:)               = 0.0
rflow_gb(:)                  = 0.0
rrun_gb(:)                   = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE fluxes_alloc

END MODULE fluxes
