! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE jules_soil_mod

USE max_dimensions, ONLY: sm_levels_max

USE jules_irrig_mod, ONLY: l_irrig_dmd

USE jules_surface_mod, ONLY: l_elev_land_ice

!-----------------------------------------------------------------------------
! Description:
!   Contains soil options and a namelist for setting them
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Switches
!-----------------------------------------------------------------------------
INTEGER, PARAMETER ::                                                         &
  sm_levels_default = 4

INTEGER ::                                                                    &
  sm_levels = sm_levels_default
      ! Number of soil layers (standalone only - this variable is not used
      ! when run in the UM)
      ! When running in the UM, its value is checked against the default and
      ! a warning is issued if it has changed

LOGICAL ::                                                                    &
  l_vg_soil      = .FALSE.,                                                   &
      ! Switch for using Van Genuchten soil scheme
  l_dpsids_dsdz = .FALSE.,                                                    &
      ! Switch to calculate vertical gradient of soil suction with the
      ! assumption of linearity only for fractional saturation
      ! (consistent with the calculation of hydraulic conductivity)
  l_soil_sat_down = .FALSE.,                                                  &
      ! Switch for direction of super_saturated soil moisture
      !   TRUE excess water is pushed down
      !   FALSE excess water is pushed up (as in JULES2.0)
  l_bedrock = .FALSE.,                                                        &
      ! Switch to control the presence of bedrock
      ! This is an additional layer below the hydrologically active
      ! soil where only thermal diffusion occurs
  l_tile_soil = .FALSE.,                                                      &
      ! Switch to control soil tiling.
      ! When F nsoilt = 1, when T nsoilt = nsurft
  soil_props_const_z = .FALSE.,                                               &
      ! Switch for whether soil ancils has the same values on each layer.
      ! Set in the JULES_SOIL_PROPS namelist.
  l_holdwater = .FALSE.
      ! Switch to control how supersaturated and negative soil moisture is
      ! handled in the implicit calculation. FALSE: excess/required moisture
      ! is pushed out/in from the base of the soil. TRUE: water is added/
      ! taken from an adjacent layer.

LOGICAL ::                                                                    &
! Switches that are only present in standalone JULES.
  l_broadcast_ancils = .FALSE.
      ! Switch to allow soil ancillaries to be read in as gridbox values and
      ! broadcast to all soil tiles

INTEGER ::                                                                    &
  soilhc_method = 1
      ! Switch for the calculation method of soil thermal conductivity
      !   SOILHC_METHOD=1: Method of Cox et al (1999)
      !   SOILHC_METHOD=2: Simplified Johansen (1975)
      !   SOILHC_METHOD=3: Chadburn et al. (2015) best for organic soils


!-----------------------------------------------------------------------------
! Fixed parameters
!-----------------------------------------------------------------------------
! Thermal conductivities (Source: "The Frozen Earth" p.90)
REAL(KIND=real_jlslsm), PARAMETER ::                                          &
  hcair = 0.025,                                                              &
      ! Thermal conductivity of air (W/m/K)
  hcice = 2.24,                                                               &
      ! Thermal conductivity of ice (W/m/K)
  hcwat = 0.56
      ! Thermal conductivity of liquid water (W/m/K)

! Parameters for soil moisture update
REAL(KIND=real_jlslsm), PARAMETER ::                                          &
  gamma_w = 1.0
      ! Forward timestep weighting

! Parameters for soil temperature update
INTEGER, PARAMETER ::                                                         &
  mmax = 3
      ! Maximum number of iterations on temperature

REAL(KIND=real_jlslsm), PARAMETER ::                                          &
  facur = 0.01,                                                               &
      ! Required flux conservation accuracy (W/m2)
  gamma_t = 1.0,                                                              &
      ! Forward timestep weighting
  tacur = 0.00000
      ! Required accuracy of temperature calculation (celsius)


!-----------------------------------------------------------------------------
! Parameters that can be set by the namelist
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                     &
  cs_min = 1.0e-6,                                                            &
      ! Minimum soil carbon (kg C/m2)
  zsmc = 1.0,                                                                 &
      ! Depth of layer over which soil moisture diagnostic is averaged (m)
  zst = 1.0,                                                                  &
      ! Depth of layer over which soil temperature diagnostic is averaged (m)
  confrac = 0.3
      ! Fraction of the gridbox over which convective precipitation is
      ! assumed to fall

!-----------------------------------------------------------------------------
! Bedrock parameters : these can be set by the namelist JULES_BEDROCK
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  ns_deep = 100
      ! Number of bedrock layers

REAL(KIND=real_jlslsm) ::                                                     &
  hcapdeep = 2100000.0,                                                       &
      ! Heat capacity of bedrock
  hcondeep = 8.6,                                                             &
      ! Thermal conductivity of bedrock
  dzdeep = 0.5
      ! Thickness of bedrock layers


!-----------------------------------------------------------------------------
! Variable length arrays that can be set using the namelist
!-----------------------------------------------------------------------------
! We provide a fixed-length version that is read by the namelist, then point
! the actual version to the portion of it we will use
REAL(KIND=real_jlslsm), POINTER ::                                            &
  dzsoil(:)
      ! Thicknesses of the soil layers (m)

REAL(KIND=real_jlslsm), TARGET ::                                             &
  dzsoil_io(sm_levels_max)
      ! Fixed length equivalent of dzsoil for namelist IO

! Initialise dzsoil_io to negative values
! In check_jules_soil, we verify that all values of dzsoil we will use are >= 0
DATA dzsoil_io / sm_levels_max * -1.0 /

REAL(KIND=real_jlslsm)  ::                                                    &
  dzsoil_elev = -1.0
      ! Depth of tiled bedrock subsurfaces under elevated tiles (m)

!-----------------------------------------------------------------------------
! Namelist definition for UM and standalone
!-----------------------------------------------------------------------------
NAMELIST  / jules_soil/                                                       &
! Additional parameters for standalone JULES only (for soil tiling).
    l_broadcast_ancils,                                                       &
! Levels - not used in the UM
    sm_levels,                                                                &
! Switches
    l_vg_soil, l_dpsids_dsdz, l_soil_sat_down, soilhc_method, l_bedrock,      &
    l_holdwater, l_tile_soil,                                                 &
! Parameters
    cs_min, zsmc, zst, confrac, ns_deep, hcapdeep, hcondeep,                  &
    dzdeep, dzsoil_io, dzsoil_elev



CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='JULES_SOIL_MOD'

CONTAINS

SUBROUTINE check_jules_soil(sm_levels_in)

USE ereport_mod, ONLY: ereport

!-----------------------------------------------------------------------------
! Description:
!   Checks JULES_SOIL namelist for consistency
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IMPLICIT NONE

INTEGER, INTENT(IN) :: sm_levels_in
    ! The number of soil levels
    ! This is passed in because the value in the namelist is not used when
    ! running in the UM
    ! This is due to complications with reading it from dump etc.

INTEGER :: errorstatus

!-----------------------------------------------------------------------------
! Verify that a suitable sm_levels was given in the namelist
!
!   * If we are running in the UM and the user tried to specify sm_levels
!     using the namelist, emit an error as it is probably something that
!     needs looking at
!     We use the passed in value instead
!   * If we are running standalone, we ignore the passed in value and use the
!     value from the namelist
!-----------------------------------------------------------------------------

! Check that there are not too many levels
IF ( sm_levels > sm_levels_max ) THEN
  errorstatus = 101
  CALL ereport("check_jules_soil", errorstatus,                               &
    "Too many soil layers specified - increase sm_levels_max and recompile")
END IF

! Check a suitable soilhc_method was given
IF ( soilhc_method < 1 .OR. soilhc_method > 3 ) THEN
  errorstatus = 101
  CALL ereport("check_jules_soil", errorstatus,                               &
    "soilhc_method must be 1, 2 or 3")
END IF

! Associate the dzsoil pointer with the appropriate section of dzsoil_io
dzsoil => dzsoil_io(1:sm_levels)

! Check we have sensible values for dzsoil
IF ( ANY(dzsoil(:) < 0.0) ) THEN
  errorstatus = 101
  CALL ereport("check_jules_soil", errorstatus,                               &
               "dzsoil < 0 - check namelist jules_soil")
END IF

IF (l_elev_land_ice .AND. dzsoil_elev <= 0.0) THEN
  errorstatus = 101
  CALL ereport("check_jules_soil", errorstatus,                               &
               "dzsoil_elev <= 0 - check namelist jules_soil")
END IF

! Check that zsmc and zst are within soil depth.
IF ( SUM(dzsoil(:)) < zsmc ) THEN
  errorstatus = 101
  CALL ereport("check_jules_soil", errorstatus,                               &
               "zsmc is below bottom of soil column")
END IF

IF ( SUM(dzsoil(:)) < zst ) THEN
  errorstatus = 101
  CALL ereport("check_jules_soil", errorstatus,                               &
               "zst is below bottom of soil column")
END IF

IF ( l_irrig_dmd .AND. sm_levels < 2 ) THEN
  errorstatus = 101
  CALL ereport("check_jules_soil",errorstatus,                                &
               'Irrigation demand (l_irrig_dmd=T) can only be used with'//    &
               '2 or more soil levels (sm_levels >= 2)')
END IF

IF ( l_bedrock ) THEN

  ! Check that at least one layer is given
  IF ( ns_deep < 1 ) THEN
    errorstatus = 101
    CALL ereport("check_jules_soil", errorstatus,                             &
                 'Bedrock must have at least one layer')
  END IF
END IF

END SUBROUTINE check_jules_soil


SUBROUTINE print_nlist_jules_soil()

USE jules_print_mgr, ONLY: jules_print

IMPLICIT NONE

CHARACTER(LEN=50000) :: lineBuffer


!-----------------------------------------------------------------------------


CALL jules_print('jules_soil', 'Contents of namelist jules_soil')

WRITE(lineBuffer, *) '  l_vg_soil = ', l_vg_soil
CALL jules_print('jules_soil', lineBuffer)

WRITE(lineBuffer, *) '  l_dpsids_dsdz = ', l_dpsids_dsdz
CALL jules_print('jules_soil', lineBuffer)

WRITE(lineBuffer, *) '  l_soil_sat_down = ', l_soil_sat_down
CALL jules_print('jules_soil', lineBuffer)

WRITE(lineBuffer, *) '  l_holdwater = ', l_holdwater
CALL jules_print('jules_soil', lineBuffer)

WRITE(lineBuffer, *) '  l_bedrock = ', l_bedrock
CALL jules_print('jules_soil', lineBuffer)

WRITE(lineBuffer, *) '  l_tile_soil = ', l_tile_soil
CALL jules_print('jules_soil', lineBuffer)

WRITE(lineBuffer, *) '  soilhc_method = ', soilhc_method
CALL jules_print('jules_soil', lineBuffer)

WRITE(lineBuffer, *) '  cs_min = ', cs_min
CALL jules_print('jules_soil', lineBuffer)

WRITE(lineBuffer, *) '  zsmc = ', zsmc
CALL jules_print('jules_soil', lineBuffer)

WRITE(lineBuffer, *) '  zst = ', zst
CALL jules_print('jules_soil', lineBuffer)

WRITE(lineBuffer, *) '  confrac = ', confrac
CALL jules_print('jules_soil', lineBuffer)

WRITE(lineBuffer, *) '  ns_deep = ', ns_deep
CALL jules_print('jules_soil', lineBuffer)

WRITE(lineBuffer, *) '  hcapdeep = ', hcapdeep
CALL jules_print('jules_soil', lineBuffer)

WRITE(lineBuffer, *) '  hcondeep = ', hcondeep
CALL jules_print('jules_soil', lineBuffer)

WRITE(lineBuffer, *) '  dzdeep = ', dzdeep
CALL jules_print('jules_soil', lineBuffer)

WRITE(lineBuffer, *) '  dzsoil_io = ', dzsoil_io
CALL jules_print('jules_soil', lineBuffer)

WRITE(lineBuffer, *) '  dzsoil_elev = ', dzsoil_elev
CALL jules_print('jules_soil', lineBuffer)

CALL jules_print('jules_soil',                                                &
    '- - - - - - end of namelist - - - - - -')

END SUBROUTINE print_nlist_jules_soil



END MODULE jules_soil_mod
