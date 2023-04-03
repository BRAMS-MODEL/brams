! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description: Module containing the variables for the FLake lake scheme
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.
!
!   Code Owner: Please refer to ModuleLeaders.txt and UM file CodeOwners.txt
!   This file belongs in section: Land
!
MODULE lake_mod
USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

! parameters following values in FLake routine flake_parameters
REAL(KIND=real_jlslsm), PARAMETER  :: h_snow_min_flk = 1.0e-5  & ! (m)
                   ,h_ice_min_flk  = 1.0e-9  & ! (m)
                   ,h_ice_max      = 3.0       ! (m)

! parameters following values in FLake routine flake_albedo_ref
REAL(KIND=real_jlslsm), PARAMETER  :: albedo_whiteice_ref =  0.60 & ! White ice
                   ,albedo_blueice_ref  =  0.10 & ! Blue ice
                   ,c_albice_MR         = 95.6
                      ! Constant in the interpolation formula for
                      ! the ice albedo (Mironov and Ritter 2004)

! parameters for empirical initialisation of lake-ice thickness
REAL(KIND=real_jlslsm), PARAMETER :: lake_h_deltaT =  25.0 ! (K)
REAL(KIND=real_jlslsm), PARAMETER :: lake_h_scale  = 100.0 ! dimensionless

! approximate lengthscale for exponential DWSW attenuation in snow
REAL(KIND=real_jlslsm), PARAMETER :: h_snow_sw_att = 0.25 ! (m)

REAL(KIND=real_jlslsm), DIMENSION(:,:), ALLOCATABLE :: surf_ht_flux_lake_ij
                            ! Net downward heat flux at surface over
                            ! lake fraction of gridbox, all points (W/m2)

REAL(KIND=real_jlslsm), DIMENSION(:), ALLOCATABLE :: SURF_HT_FLUX_LK_gb
                            ! Net downward heat flux at surface over
                            ! lake fraction of gridbox, land points (W/m2)
REAL(KIND=real_jlslsm), DIMENSION(:), ALLOCATABLE :: U_S_LAKE_gb
                            ! lake subsurface friction velocity (m/s)
REAL(KIND=real_jlslsm), DIMENSION(:), ALLOCATABLE :: SW_DOWN_gb
                            ! downwelling shortwave irradiance [W m^{-2}]
REAL(KIND=real_jlslsm), DIMENSION(:), ALLOCATABLE :: CORIOLIS_PARAM_gb
                            ! Coriolis parameter (s^-1)
REAL(KIND=real_jlslsm), DIMENSION(:), ALLOCATABLE :: LAKE_DEPTH_gb
                            ! lake depth (m)
REAL(KIND=real_jlslsm), DIMENSION(:), ALLOCATABLE :: LAKE_FETCH_gb
                            ! Typical wind fetch (m)
REAL(KIND=real_jlslsm), DIMENSION(:), ALLOCATABLE, SAVE :: LAKE_ALBEDO_gb
                            ! lake albedo
REAL(KIND=real_jlslsm), DIMENSION(:), ALLOCATABLE :: LAKE_T_SNOW_gb
                            ! temperature at the air-snow interface (K)
REAL(KIND=real_jlslsm), DIMENSION(:), ALLOCATABLE :: LAKE_T_ICE_gb
                            ! temperature at upper boundary of lake ice(K)
REAL(KIND=real_jlslsm), DIMENSION(:), ALLOCATABLE :: LAKE_T_MEAN_gb
                            ! lake mean temperature (K)
REAL(KIND=real_jlslsm), DIMENSION(:), ALLOCATABLE :: LAKE_T_MXL_gb
                            ! lake mixed-layer temperature (K)
REAL(KIND=real_jlslsm), DIMENSION(:), ALLOCATABLE :: LAKE_SHAPE_FACTOR_gb
                            ! thermocline shape factor (dimensionless?)
REAL(KIND=real_jlslsm), DIMENSION(:), ALLOCATABLE :: LAKE_H_SNOW_gb
                            ! snow thickness (m)
REAL(KIND=real_jlslsm), DIMENSION(:), ALLOCATABLE :: LAKE_H_ICE_gb
                            ! lake ice thickness (m)
REAL(KIND=real_jlslsm), DIMENSION(:), ALLOCATABLE :: LAKE_H_MXL_gb
                            ! lake mixed-layer thickness (m)
REAL(KIND=real_jlslsm), DIMENSION(:), ALLOCATABLE :: LAKE_T_SFC_gb
                            ! temperature (of water, ice or snow) at
                            !surface of lake (K)
REAL(KIND=real_jlslsm), DIMENSION(:), ALLOCATABLE :: TS1_LAKE_gb
                            ! "average" temperature of
                            ! lake-ice, lake and soil sandwich
                            ! (K).
REAL(KIND=real_jlslsm), DIMENSION(:), ALLOCATABLE :: NUSSELT_gb
                                  ! Nusselt number
REAL(KIND=real_jlslsm), DIMENSION(:), ALLOCATABLE :: G_DT_gb
                            ! ground heat flux over delta T [W m-2 K-1]

! initial values
REAL(KIND=real_jlslsm) :: lake_depth_0  =    5.0  ! (m)
REAL(KIND=real_jlslsm) :: lake_fetch_0  =   25.0  ! (m)
REAL(KIND=real_jlslsm) :: lake_h_mxl_0  =    2.0  ! (m)
REAL(KIND=real_jlslsm) :: lake_shape_0  =    0.5  ! shape factor, dimensionless
REAL(KIND=real_jlslsm) :: nusselt_0     = 1000.0
                                ! Nusselt number, dimensionless
REAL(KIND=real_jlslsm) :: g_dt_0        = 1.0e-10 ! very small [W m-2 K-1]

REAL :: lake_T_mean_0 = 280.0   ! (K) ---|      Vaguely physically
REAL :: lake_T_mxl_0  = 275.0   ! (K)    |_____ plausible initial values,
REAL :: lake_T_ice_0  = 270.0   ! (K)    |      will be overwritten by 
REAL :: lake_H_ice_0  = 3.0     ! (mm)---|      flake_init_mod. 

! trap counters
INTEGER :: trap_frozen   = 0
INTEGER :: trap_unfrozen = 0

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LAKE_MOD'

CONTAINS
  
SUBROUTINE lake_alloc(land_pts,                                               &
                      l_flake_model)

!No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

!Except this one for now
USE atm_fields_bounds_mod,    ONLY: tdims

IMPLICIT NONE

!Arguments
INTEGER, INTENT(IN) :: land_pts

LOGICAL, INTENT(IN) :: l_flake_model

!Local variables
INTEGER :: temp_size

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LAKE_ALLOC'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!  ====lake_mod module UM====
IF ( l_flake_model ) THEN
  ALLOCATE( surf_ht_flux_lake_ij(tdims%i_start:tdims%i_end                    &
                             ,tdims%j_start:tdims%j_end))
  temp_size = land_pts
ELSE
  ALLOCATE( surf_ht_flux_lake_ij(1,1))
  temp_size = 1
END IF

surf_ht_flux_lake_ij(:,:) = 0.0

ALLOCATE( surf_ht_flux_lk_gb(temp_size))
ALLOCATE( sw_down_gb(temp_size))
ALLOCATE( coriolis_param_gb(temp_size))
ALLOCATE( u_s_lake_gb(temp_size))
ALLOCATE( lake_depth_gb(temp_size))
ALLOCATE( lake_fetch_gb(temp_size))
ALLOCATE( lake_albedo_gb(temp_size))
ALLOCATE( lake_t_snow_gb(temp_size))
ALLOCATE( lake_t_ice_gb(temp_size))
ALLOCATE( lake_t_mean_gb(temp_size))
ALLOCATE( lake_t_mxl_gb(temp_size))
ALLOCATE( lake_shape_factor_gb(temp_size))
ALLOCATE( lake_h_snow_gb(temp_size))
ALLOCATE( lake_h_ice_gb(temp_size))
ALLOCATE( lake_h_mxl_gb(temp_size))
ALLOCATE( lake_t_sfc_gb(temp_size))
ALLOCATE( ts1_lake_gb(temp_size))
ALLOCATE( nusselt_gb(temp_size))
ALLOCATE( g_dt_gb(temp_size))

surf_ht_flux_lk_gb(:)   = 0.0
sw_down_gb(:)           = 0.0
coriolis_param_gb(:)    = 0.0
u_s_lake_gb(:)          = 0.0
lake_depth_gb(:)        = 0.0
lake_fetch_gb(:)        = 0.0
lake_albedo_gb(:)       = 0.0
lake_t_snow_gb(:)       = 0.0
lake_t_ice_gb(:)        = 0.0
lake_t_mean_gb(:)       = 0.0
lake_t_mxl_gb(:)        = 0.0
lake_shape_factor_gb(:) = 0.0
lake_h_snow_gb(:)       = 0.0
lake_h_ice_gb(:)        = 0.0
lake_h_mxl_gb(:)        = 0.0
lake_t_sfc_gb(:)        = 0.0
ts1_lake_gb(:)          = 0.0
nusselt_gb(:)           = 0.0
g_dt_gb(:)              = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lake_alloc

END MODULE
