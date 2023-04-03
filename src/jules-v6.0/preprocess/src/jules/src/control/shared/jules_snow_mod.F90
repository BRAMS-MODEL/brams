! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE jules_snow_mod

USE max_dimensions, ONLY: snow_layers_max, npft_max, nsurft_max

!-----------------------------------------------------------------------------
! Description:
!   Contains snow options and a namelist for setting them
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

USE missing_data_mod, ONLY: rmdi
USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

INTEGER ::                                                                    &
  nsmax = 0        !  Maximum number of snow layers

!-----------------------------------------------------------------------------
! Switches
!-----------------------------------------------------------------------------
LOGICAL ::                                                                    &
  l_snowdep_surf = .FALSE.,                                                   &
                   ! use equivalent canopy snow depth for surface
                   ! calculations on tiles with a snow canopy
  l_rho_snow_corr  = .TRUE.,                                                  &
                   ! Switch for using a correction to the density of the
                   ! snow pack when nsnow=0 when relayering in the new snow
                   ! scheme
                   ! Has no effect for nsmax < 1
  l_et_metamorph = .FALSE.,                                                   &
                   ! Switch to include equitemperature metamorphism
  l_snow_infilt  = .FALSE.,                                                   &
                   ! Include infiltration of rain water into snow
                   ! (Only applies to the multilayer scheme)
  l_snow_nocan_hc = .FALSE.
                   ! Flag to negelect canopy heat capacity in the presence
                   ! of snow if not using an explicit canopy model

INTEGER ::                                                                    &
  frac_snow_subl_melt = 0
                   ! Switch for use of snow-cover fraction in the calculation
                   ! of sublimation and melting
                   !   0 = off
                   !   1 = on
! Parametrizations of thermal conductivity

INTEGER, PARAMETER ::                                                         &
  ip_snow_cond_yen81 = 0,                                                     &
                   ! Parametrization of snow conductivity following
                   ! Yen (1981), "Review of the thermal properties of snow,
                   ! ice and sea ice", Technical Report 81-10
                   ! Cold Regions Research and Engineering Laboratory,
                   ! Hanover, NH
  ip_snow_cond_calonne11 = 1
                   ! Parametrization of snow conductivity following
                   ! Calonne et al. (2011), GRL, 38, L23501

INTEGER ::                                                                    &
  i_snow_cond_parm = 0
                   ! Parametrization scheme for snow conductivity

!-----------------------------------------------------------------------------
! Parametrization of the rate of growth of snow grains
INTEGER, PARAMETER ::                                                         &
  ip_grain_growth_marshall89 = 0,                                             &
                   ! Original scheme following Marshall (1989)
  ip_grain_growth_taillandier_et = 1
                   ! Rate of growth of grains at temperatures below
                   ! freezing follows the equitemperature scheme of
                   ! Taillandier et al. (2007), J. Geophys. Res., 112,
                   ! F03003

INTEGER ::                                                                    &
  i_grain_growth_opt = 0
                   ! Option for rate of growth of snow grains

!-----------------------------------------------------------------------------
! Options for relayering the snow pack.
INTEGER, PARAMETER ::                                                         &
  ip_relayer_linear = 0,                                                      &
                   ! Relayer based on thickness
  ip_relayer_rgrain_inv = 1
                   ! As the linear scheme, but relayering the inverse of
                   ! the grain size

INTEGER ::                                                                    &
  i_relayer_opt = 0
                    ! Option for relayering the snow pack.

!-----------------------------------------------------------------------------
! Options for melting at the base of thin snow packs.
INTEGER, PARAMETER ::                                                         &
  ip_basal_melting_off = 0,                                                   &
                   ! Basal melting not considered in the zero-layer scheme
  ip_basal_melting_diag_0 = 1
                   ! Basal melting is diagnosed for thin snow packs
                   ! at points where the zero-layer scheme is active.
                   ! The multilayer scheme is unaffected.

INTEGER ::                                                                    &
  i_basal_melting_opt = 0
                    ! Option for basal melting

! ----------------------------------------------------------------------------
! Parametrisations for treatment of graupel
INTEGER ::                                                                    &
 graupel_options = 0
                   ! Switch for treatment of graupel in the snow scheme

INTEGER, PARAMETER ::                                                         &
  graupel_as_snow = 0,                                                        &
                   ! Include graupel in the surface snowfall
  ignore_graupel = 1
                   ! Ignore graupel in the surface snowfall
!-----------------------------------------------------------------------------
! Parameters for equitemperature metamorphism
!
REAL(KIND=real_jlslsm) ::                                                     &
  a_snow_et = 0.0,                                                            &
  b_snow_et = 0.0,                                                            &
  c_snow_et = 0.0
                   ! Parameters in rate equation for ET metamorphism
REAL(KIND=real_jlslsm) ::                                                     &
  rho_snow_et_crit = 0.0
                   ! Critical density

!-----------------------------------------------------------------------------
! Radiation parameters for namelist
!
! Defaults taken from JULES examples
!-----------------------------------------------------------------------------
! Parameters for snow grain size.
REAL(KIND=real_jlslsm) ::                                                     &
  r0 = 50.0,                                                                  &
                !  Grain size for fresh snow (microns)
  rmax = 2000.0
                !  Maximum snow grain size (microns)


! Parameters for snow grain size growth (HCTN30.16)
! Values are for melting snow, cold fresh snow and cold aged snow respectively
REAL(KIND=real_jlslsm) ::                                                     &
  snow_ggr(3)    !  snow grain area growth rates (microns**2 s-1)
DATA snow_ggr / 0.6, 0.06, 0.23e6 /


! Parameters for prognostic, spectral snow albedo.
REAL(KIND=real_jlslsm) ::                                                     &
 amax(2)        !  Maximum albedo for fresh snow (values for VIS and NIR)
DATA amax / 0.98, 0.7 /

REAL(KIND=real_jlslsm) ::                                                     &
 aicemax(2)     !  Maximum albedo for bare land ice (values for VIS and NIR)
                !  Values hand tuned by Sarah Shannon to match a range of
                !  glaciers in standalone JULES (pers comm.)
DATA aicemax / 0.78, 0.36 /

!-----------------------------------------------------------------------
! Parameters for prognostic, spectral snow albedo using the embedded
! canopy scheme.
!-----------------------------------------------------------------------
REAL(KIND=real_jlslsm) :: dce = 2.4848e-7
               !  Effective dimension of black carbon
REAL(KIND=real_jlslsm) :: cnr_g(2, 2) = RESHAPE(                              &
          (/ 0.900396, 0.00143416, 0.947363, 0.00570307 /),                   &
          (/ 2, 2 /) )
               !  Parameters fitting the asymmetry
REAL(KIND=real_jlslsm) :: cnr_om(3, 2, 2) = RESHAPE(                          &
          (/ 0.000869870, 0.000194185, 1.09046e-05,                           &
             0.0, 0.2709248, 0.0,                                             &
             0.103851, 0.0201497, 0.00101678,                                 &
             0.0, 0.4306217, 0.0                       /),                    &
          (/ 3, 2, 2 /) )
               !  Parameters fitting the single scattering albedo

! Parameters for (diagnostic) all-band snow albedo.
REAL(KIND=real_jlslsm) ::                                                     &
  dtland = 2.0,                                                               &
                ! Degrees C below freezing point at which snow albedo equals
                ! cold deep snow albedo. This is 2 in HCTN30.4.
                ! Must not be zero!
  kland_numerator = 0.3,                                                      &
                ! KLAND is calculated as KLAND_NUMERATOR / DTLAND once the
                ! namelist has been read
  maskd = 50.0
              ! Used in the exponent of equation weighting snow and snow-free
              ! albedo to get tile albedo
              ! This is 0.2 in HCTN30.5, where it is used as maskd*snowMass,
              ! assuming snow density=250kg/m3
              ! It is now used as maskd*snowDepth, so maskd=50 gives the same
              ! relationship as HCTN30.5

!-----------------------------------------------------------------------------
! Calculated radiation parameters
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) :: kland
                ! Used in snow-ageing effect on snow albedo
                ! This is 0.3 in HCTN30.4, although note that the last term
                ! in that eqn should be divided by dtland
                ! KLAND is calculated as KLAND_NUMERATOR / DTLAND once the
                ! namelist has been read
REAL(KIND=real_jlslsm) :: tcland
                ! Temperature below which snow albedo equals cold deep snow
                ! albedo (HCTN30.4)
                ! This is initialised as TM - DTLAND once the namelist has
                ! been read

!-----------------------------------------------------------------------------
! Other snow parameters for namelist
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                     &
  rho_snow_const = 250.0,                                                     &
                ! constant density of lying snow (kg per m**3)
                ! This is used:
                !   (a) as the snow density when nsmax=0
                !   (b) with nsmax>0 and l_snowdep_surf=.TRUE. this
                !       is the density of snow on the canopy
                !       (not on the ground)
  rho_snow_fresh = 100.0,                                                     &
                ! Density of fresh snow (kg per m**3)
                ! Only used with nsmax>0.
  snow_hcon = 0.265,                                                          &
                ! Thermal conductivity of lying snow (Watts per m per K)
  snow_hcap = 0.63e6,                                                         &
                ! Thermal capacity of lying snow (J/K/m3)
  snowliqcap = 0.05,                                                          &
                ! Liquid water holding capacity of lying snow as a fraction
                ! of snow mass
  snowinterceptfact = 0.7,                                                    &
                ! Constant in relationship between mass of intercepted snow
                ! and snowfall rate
  snowloadlai = 4.4,                                                          &
                ! Ratio of maximum canopy snow load to leaf area index(kg m-2)
  snowunloadfact = 0.4,                                                       &
                ! Constant in relationship between canopy snow unloading and
                ! canopy snow melt rate
  rho_firn_albedo = 550.0,                                                    &
                ! Threshold surface density where albedo parameterisation
                ! switches to scaling with surface density not grain size
                ! (kg per m**3)
  rho_firn_pore_restrict = 450.0,                                             &
                ! Threshold snowpack density where holding capacity/ability
                ! to percolate meltwater starts to be reduced
                ! (kg per m**3)
  rho_firn_pore_closure = 850.0
                ! Threshold snowpack density where holding capacity/ability
                ! to percolate meltwater becomes 0
                ! (kg per m**3)

!-----------------------------------------------------------------------------
! Allocatable versions of namelist arrays
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) :: dzsnow(snow_layers_max) = rmdi
                ! Prescribed thickness of snow layers (m)
                ! This is the thickness of each snow layer when it is not
                ! the bottom layer (note that dzSnow(nsMax) is not used
                ! because that is always the bottom layer)

LOGICAL :: cansnowtile(nsurft_max)
                ! Switch for canopy snow model on each tile
                ! Must be false for non-PFT tiles
DATA cansnowtile / nsurft_max * .FALSE. /

LOGICAL :: cansnowpft(npft_max)
                ! This is what appears in the namelist, to ensure that true
                ! values can only given for pfts
DATA cansnowpft / npft_max * .FALSE. /

REAL(KIND=real_jlslsm) :: unload_rate_cnst(npft_max)
!                 !  Constant canopy unloading rate for the tile
REAL(KIND=real_jlslsm) :: unload_rate_u(npft_max)
!                   !  Wind-speed dependent canopy unloading rate
!                   !  for the tile

DATA unload_rate_cnst / npft_max * 0.0 /
DATA unload_rate_u / npft_max * 0.0 /

REAL(KIND=real_jlslsm) :: can_clump(npft_max)
                 ! Clumping factor for snow in the canopy
REAL(KIND=real_jlslsm) :: lai_alb_lim_sn(npft_max)
                 ! Lower limit on permitted LAI in albedo with snow
REAL(KIND=real_jlslsm) :: n_lai_exposed(npft_max)
                 ! Shape parameter for exposed canopy with
                 ! embedded snow

DATA can_clump / npft_max * 0.0 /
DATA lai_alb_lim_sn / npft_max * 0.0 /
DATA n_lai_exposed / npft_max * 0.0 /


!-----------------------------------------------------------------------------
! Single namelist definition for UM and standalone
!-----------------------------------------------------------------------------
NAMELIST  / jules_snow/                                                       &
! Switches
    nsmax, l_snowdep_surf, l_rho_snow_corr, frac_snow_subl_melt,              &
    graupel_options,                                                          &
! Equitemperature metamorphism
    l_et_metamorph, a_snow_et, b_snow_et, c_snow_et, rho_snow_et_crit,        &
! Thermal conductivity of snow
    i_snow_cond_parm,                                                         &
! Radiation parameters
    r0, rmax, snow_ggr, amax, aicemax, maskd, dtland, kland_numerator,        &
! Other snow parameters
    rho_snow_const, rho_snow_fresh, rho_firn_albedo,                          &
    snow_hcon, snow_hcap, snowliqcap,                                         &
    snowinterceptfact, snowloadlai, snowunloadfact, dzsnow, cansnowpft,       &
    unload_rate_cnst, unload_rate_u, can_clump, lai_alb_lim_sn,               &
    n_lai_exposed, l_snow_infilt, l_snow_nocan_hc,                            &
    i_grain_growth_opt, i_relayer_opt, i_basal_melting_opt

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='JULES_SNOW_MOD'

CONTAINS

SUBROUTINE check_jules_snow()

USE ereport_mod, ONLY: ereport

USE water_constants_mod, ONLY: tm

USE jules_surface_types_mod, ONLY: npft
USE jules_surface_mod, ONLY: l_aggregate
USE jules_vegetation_mod, ONLY: can_model

!-----------------------------------------------------------------------------
! Description:
!   Checks JULES_SNOW namelist for consistency
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IMPLICIT NONE

INTEGER :: errorstatus

!-----------------------------------------------------------------------------


! Check that we have a suitable number of snow layers
IF ( nsmax > snow_layers_max ) THEN
  errorstatus = 101
  CALL ereport("check_jules_snow", errorstatus,                               &
  "Too many snow layers specified - increase snow_layers_max and recompile")
END IF

! Copy the values for the given number of pfts across
IF ( .NOT. l_aggregate .AND. can_model == 4 ) THEN
  canSnowTile(1:npft) = cansnowpft(1:npft)
END IF

! Derive kland and dtland
kland = kland_numerator / dtland
tcland = tm - dtland

END SUBROUTINE check_jules_snow


SUBROUTINE print_nlist_jules_snow()

USE jules_print_mgr, ONLY: jules_print

IMPLICIT NONE

CHARACTER(LEN=50000) :: lineBuffer


!-----------------------------------------------------------------------------


CALL jules_print('jules_snow',                                                &
                 'Contents of namelist jules_snow')

WRITE(lineBuffer, *) '  nsmax = ', nsmax
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  l_snowdep_surf = ', l_snowdep_surf
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  l_rho_snow_corr = ', l_rho_snow_corr
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  graupel_options = ', graupel_options
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  frac_snow_subl_melt = ', frac_snow_subl_melt
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  l_snow_infilt = ', l_snow_infilt
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  l_snow_nocan_hc = ', l_snow_nocan_hc
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  l_et_metamorph = ', l_et_metamorph
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  a_snow_et = ', a_snow_et
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  b_snow_et = ', b_snow_et
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  c_snow_et = ', c_snow_et
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  rho_snow_et_crit = ', rho_snow_et_crit
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  i_snow_cond_parm = ', i_snow_cond_parm
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  r0 = ', r0
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  rmax = ', rmax
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  snow_ggr = ', snow_ggr
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  amax = ', amax
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  aicemax = ', aicemax
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  maskd = ', maskd
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  dtland = ', dtland
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  kland_numerator = ', kland_numerator
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  rho_snow_const = ', rho_snow_const
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  rho_firn_albedo = ', rho_firn_albedo
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  rho_snow_fresh = ', rho_snow_fresh
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  snow_hcon = ', snow_hcon
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  snow_hcap = ', snow_hcap
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  snowliqcap = ', snowliqcap
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  snowinterceptfact = ', snowinterceptfact
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  snowloadlai = ', snowloadlai
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  snowunloadfact = ', snowunloadfact
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  dzsnow = ', dzsnow
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer, *) '  cansnowpft = ', cansnowpft
CALL jules_print('jules_snow', lineBuffer)
WRITE(lineBuffer,*)' unload_rate_cnst = ', unload_rate_cnst
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer,*)' unload_rate_u = ', unload_rate_u
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer,*)' can_clump = ', can_clump
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer,*)' lai_alb_lim_sn = ', lai_alb_lim_sn
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer,*)' n_lai_exposed = ', n_lai_exposed
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer,'(a,i3)')' i_grain_growth_opt = ', i_grain_growth_opt
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer,'(a,i3)')' i_relayer_opt = ', i_relayer_opt
CALL jules_print('jules_snow', lineBuffer)

WRITE(lineBuffer,'(a,i3)')' i_basal_melting_opt = ', i_basal_melting_opt
CALL jules_print('jules_snow', lineBuffer)

CALL jules_print('jules_snow',                                                &
    '- - - - - - end of namelist - - - - - -')

END SUBROUTINE print_nlist_jules_snow


END MODULE jules_snow_mod
