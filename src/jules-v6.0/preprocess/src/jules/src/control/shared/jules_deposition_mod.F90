!******************************COPYRIGHT**************************************
! (c) Centre for Ecology and Hydrology.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
! ******************************** COPYRIGHT *********************************

!-----------------------------------------------------------------------------
! Description:
!   Contains dry deposition options and a namelist for setting them, and
!   variables for the schemes.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in BIOGENIC FLUXES.
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

MODULE jules_deposition_mod

USE max_dimensions,          ONLY: ndep_species_max

USE missing_data_mod,        ONLY: imdi, rmdi

USE ereport_mod,             ONLY: ereport

USE jules_print_mgr,         ONLY: jules_print

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

! Public scope by default.

!-----------------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='JULES_DEPOSITION_MOD'

! Parameters identifying alternative dry deposition schemes.
! These should be unique.
INTEGER, PARAMETER ::                                                         &
  dry_dep_ukca_jules = 1
    ! Deposition is calculated in JULES using code modelled on UKCA.

!-----------------------------------------------------------------------------
! Module variables - these are included in the namelist.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Variables that can apply to more than one deposition model.
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  dry_dep_model = imdi,                                                       &
    ! Indicates the chosen model for dry deposition.
    ! Valid values are given by the dry_dep_* parameters.
  ndry_dep_species = imdi
    ! Number of species for which dry deposition is calculated.

REAL(KIND=real_jlslsm) :: dzl_const = rmdi
  ! Constant value for separation of boundary layer levels (m).
  ! All layer thicknesses are set to this value. This is used as a simple way
  ! to prescribe the layer thicknesses in standalone mode.
  ! When/if this module is added to the UM, this variable should only be
  ! present in standalone JULES.

LOGICAL ::                                                                    &
  l_deposition = .FALSE.,                                                     &
    ! Switch for calculation of atmospheric dry deposition within JULES.
  l_deposition_flux = .FALSE.
    ! Switch for calculation of deposition fluxes.
    ! Only used if deposition is requested.
    ! T means the deposition flux is calculated.
    ! F means only the deposition velocity is calculated.

!-----------------------------------------------------------------------------
! Variables for the UKCA scheme.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                     &
  tundra_s_limit = rmdi
    ! Latitude of southern limit of tundra (degrees).

LOGICAL ::                                                                    &
  l_ukca_ddep_lev1  = .FALSE.
    ! Switch controlling which atmospheric levels are used for dry deposition.
    ! T means only the lowest layer is used.
    ! F means all layers in the boundary layer are used.
    ! **** NOTE ****
    ! In future (when the code is added) if this is FALSE we will be looping
    ! through boundary layer levels (held in dzl) until we reach the top of
    ! the boundary layer (at height zh). At present I am not 100% sure that
    ! the UM always ensures that zh<=SUM(dzl), though I think it does. This
    ! should be checked with UM code owners and if necessary the standalone
    ! code should check that zh<=SUM(dzl) on every timestep. Until then do
    ! not assume that the boundary layer top will be found!
    ! **** END OF NOTE ****

!-----------------------------------------------------------------------------
! Namelist definitions.
!-----------------------------------------------------------------------------
NAMELIST  / jules_deposition /                                                &
  dzl_const, dry_dep_model, l_deposition, l_deposition_flux,                  &
  l_ukca_ddep_lev1, ndry_dep_species, tundra_s_limit

!-----------------------------------------------------------------------------
! Variables for the dry deposition scheme(s).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  tracer_field(:,:)
    ! Concentration of chemical tracers in the atmosphere, for the
    ! calculation of deposition, as mass mixing ratio (kg kg-1).

CONTAINS

!#############################################################################

SUBROUTINE jules_deposition_alloc(land_pts)

!No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

USE ereport_mod, ONLY: ereport

IMPLICIT NONE

!Arguments
INTEGER, INTENT(IN) :: land_pts

!Local variables
INTEGER :: temp_size, temp_tiles, temp_layers

INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='JULES_DEPOSITION_ALLOC'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

errcode = 101

IF ( l_deposition ) THEN

  ! Model-dependent variables.
  SELECT CASE ( dry_dep_model )
  CASE ( dry_dep_ukca_jules )
    ! Deposition variables.
    IF ( l_deposition_flux ) THEN
      ALLOCATE( tracer_field(land_pts,ndry_dep_species))
      
      tracer_field(:,:) = 0.0
      
    END IF  !  l_deposition_flux
  CASE DEFAULT
    CALL ereport ("Invalid value for dry_dep_model. ", errcode,               &
                   "Please check jules_deposition_alloc")
  END SELECT
END IF  !  l_deposition


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE jules_deposition_alloc


SUBROUTINE check_jules_deposition()

!-----------------------------------------------------------------------------
! Description:
!   Checks values from the JULES_DEPOSITION namelist.
!-----------------------------------------------------------------------------

USE jules_surface_mod, ONLY:                                                  &
  l_aggregate

USE ereport_mod, ONLY: ereport

IMPLICIT NONE

INTEGER :: errorstatus

CHARACTER(LEN=*), PARAMETER ::                                                &
  RoutineName = 'CHECK_JULES_DEPOSITION'   ! Name of this procedure.

!-----------------------------------------------------------------------------
!end of header

IF ( .NOT. l_deposition ) THEN
  ! Nothing more to do here.
  RETURN
END IF

!-----------------------------------------------------------------------------
! Check that a valid dry deposition model is selected.
!-----------------------------------------------------------------------------
SELECT CASE ( dry_dep_model )
CASE ( dry_dep_ukca_jules )
  !  Acceptable values.
CASE DEFAULT
  errorstatus = 100
  CALL ereport( TRIM(routineName), errorstatus,                               &
                "Invalid value for dry_dep_model." )
END SELECT

!-----------------------------------------------------------------------------
! Deposition requires a tiled model.
!-----------------------------------------------------------------------------
IF ( l_aggregate ) THEN
  errorstatus = 100
  CALL ereport( TRIM(routineName), errorstatus,                               &
                "Deposition cannot use aggregate tiles." )
END IF

!-----------------------------------------------------------------------------
! Ensure number of species is reasonable.
!-----------------------------------------------------------------------------
IF ( ndry_dep_species < 1 .OR. ndry_dep_species > ndep_species_max ) THEN
  errorstatus = 100  !  a hard error
  CALL ereport( TRIM(routineName), errorstatus,                               &
                "Number of species must be in range 1 to " //                 &
                "ndep_species_max." )
END IF

!-----------------------------------------------------------------------------
! Ensure layer thickness is reasonable.
!-----------------------------------------------------------------------------
IF ( dzl_const <= 0.0 ) THEN
  errorstatus = 100  !  a hard error
  CALL ereport( TRIM(routineName), errorstatus,                               &
                "dzl_const must be greater than zero." )
END IF

!-----------------------------------------------------------------------------
! Ensure tundra_s_limit was provided - any non-missing value is allowed.
!-----------------------------------------------------------------------------
IF ( ABS( tundra_s_limit - rmdi ) < EPSILON(tundra_s_limit) ) THEN
  errorstatus = 100  !  a hard error
  CALL ereport( TRIM(routineName), errorstatus,                               &
                "tundra_s_limit must be provided." )
END IF

END SUBROUTINE check_jules_deposition

!-----------------------------------------------------------------------------

SUBROUTINE print_nlist_jules_deposition()

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER ::                                                &
  RoutineName = 'PRINT_NLIST_JULES_DEPOSITION'   ! Name of this procedure.

CHARACTER(LEN=50000) :: lineBuffer

CALL jules_print(RoutineName,                                                 &
                 'Contents of namelist jules_deposition')

WRITE(lineBuffer,"(A,L1)") '  l_deposition      = ', l_deposition
CALL jules_print(RoutineName, lineBuffer)

WRITE(lineBuffer,"(A,L1)") '  l_deposition_flux = ', l_deposition_flux
CALL jules_print(RoutineName, lineBuffer)

WRITE(lineBuffer,"(A,I0)") '  dry_dep_model     = ', dry_dep_model
CALL jules_print(RoutineName, lineBuffer)

WRITE(lineBuffer,"(A,L1)") '  l_ukca_ddep_lev1  = ', l_ukca_ddep_lev1
CALL jules_print(RoutineName, lineBuffer)

WRITE(lineBuffer,"(A,I0)") '  ndry_dep_species  = ', ndry_dep_species
CALL jules_print(RoutineName, lineBuffer)

WRITE(lineBuffer,"(A,G11.4E2)") '  tundra_s_limit    = ', tundra_s_limit
CALL jules_print(RoutineName, lineBuffer)

CALL jules_print(RoutineName, '- - - - - - end of namelist - - - - - -')

END SUBROUTINE print_nlist_jules_deposition

!-----------------------------------------------------------------------------

END MODULE jules_deposition_mod
