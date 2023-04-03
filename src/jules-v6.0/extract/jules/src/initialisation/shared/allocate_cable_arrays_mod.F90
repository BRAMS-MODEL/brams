MODULE allocate_cable_arrays_mod

!Common Non-science modules

USE yomhook,                  ONLY: lhook, dr_hook
USE ereport_mod,              ONLY: ereport
USE parkind1,                 ONLY: jprb, jpim

IMPLICIT NONE

PRIVATE

PUBLIC :: allocate_cable_arrays, allocate_cable_progs

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ALLOCATE_CABLE_ARRAYS_MOD'

CONTAINS

SUBROUTINE allocate_cable_arrays()

USE cable_types_mod,      ONLY: mp, soil, veg
USE ancil_info,           ONLY: surft_pts

!-----------------------------------------------------------------------------
! Description:
!   Allocates the CABLE model arrays using sizes determined during
!   initialisation
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

! Determine the number of active tiles
mp = SUM(surft_pts)

CALL allocate_soil_parameter_type(soil, mp)
CALL allocate_veg_parameter_type(veg, mp)

END SUBROUTINE allocate_cable_arrays

SUBROUTINE allocate_veg_parameter_type(var, mp)

USE cable_types_mod,            ONLY: veg_parameter_type
USE cable_other_constants_mod,  ONLY: nrb
USE jules_soil_mod,             ONLY: sm_levels   ! number of soil levels

!-----------------------------------------------------------------------------
! Description:
!   Allocates variable arrays for vegetation parameters.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

TYPE(veg_parameter_type), INTENT(INOUT) :: var
INTEGER, INTENT(IN) :: mp

INTEGER ::                                                                    &
  error        = 0,                                                           &
                       ! Variable for trapping the error from each
                       ! individual call to allocate
  error_sum    = 0,                                                           &

                       ! Variable to track the sum of all errors
                       ! resulting from calls to allocate. Hence we
                       ! know that everything was successful if and
                       ! only if this is zero at the end
  errcode      = 101
                       ! Variable to use in error report

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_VEG_PARAMETER_TYPE'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE( var% iveg(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var% hc(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var% vlai(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var%refl(mp,nrb), stat = error )
error_sum = error_sum + error
ALLOCATE( var%taul(mp,nrb), stat = error )
error_sum = error_sum + error

IF (error_sum == 0) THEN
  var%iveg(:) = 0
  var%hc(:) = 0.0
  var%vlai(:) = 0.0
  var%refl(:,:) = 0.0
  var%taul(:,:) = 0.0
ELSE
  IF ( error_sum /= 0 )                                                       &
    CALL ereport("allocate_cable_arrays", errcode,                            &
                 "Error allocating CABLE model veg parameter arrays")
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE allocate_veg_parameter_type


SUBROUTINE allocate_soil_parameter_type(var, mp)

USE cable_types_mod,            ONLY: soil_parameter_type
USE cable_other_constants_mod,  ONLY: nrb
USE jules_soil_mod,             ONLY: sm_levels   ! number of soil levels

!-----------------------------------------------------------------------------
! Description:
!   Allocates soil parameter arrays.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

TYPE(soil_parameter_type), INTENT(INOUT) :: var
INTEGER, INTENT(IN) :: mp

INTEGER ::                                                                    &
  error        = 0,                                                           &
                       ! Variable for trapping the error from each
                       ! individual call to allocate
  error_sum    = 0,                                                           &

                       ! Variable to track the sum of all errors
                       ! resulting from calls to allocate. Hence we
                       ! know that everything was successful if and
                       ! only if this is zero at the end
  errcode      = 101
                       ! Variable to use in error report

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_SOIL_PARAMETER_TYPE'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE( var% isoilm(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var% albsoil(mp, nrb), stat = error )
error_sum = error_sum + error
ALLOCATE( var% albsoilf(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var% soilcol(mp), stat = error )
error_sum = error_sum + error

IF (error_sum == 0) THEN
  var%isoilm(:) = 0
  var%albsoil(:,:) = 0.0
  var%albsoilf(:) = 0.0
  var%soilcol(:) = 0.0
ELSE
  IF ( error_sum /= 0 )                                                       &
    CALL ereport("allocate_cable_arrays", errcode,                            &
                 "Error allocating CABLE model soil parameter arrays")
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE allocate_soil_parameter_type


SUBROUTINE allocate_cable_progs()

!Replacements for the argument list
USE ancil_info,               ONLY:                                           &
  land_pts
USE jules_soil_mod,           ONLY:                                           &
  sm_levels
USE cable_prognostic_info_mod, ONLY: soiltemp_cable, soilmoisture_cable,      &
  frozensoilfrac_cable, snowdepth_cable, snowmass_cable, snowtemp_cable,      &
  snowdensity_cable, threelayersnowflag_cable, onelyrsnowdensity_cable,       &
  snowage_cable
USE cable_surface_types_mod,  ONLY:                                           &
  ntype_cable

!Common Non-science modules
USE jules_print_mgr,          ONLY:                                           &
  jules_message, jules_print, PrNorm

USE jules_snow_mod,           ONLY:                                           &
  nsmax

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Allocates CABLE prognostic variable arrays.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------
! Local variables for error trapping
!-----------------------------------------------------------------------
INTEGER ::                                                                    &
  error        = 0,                                                           &
                       ! Variable for trapping the error from each
                       ! individual call to allocate
  error_sum    = 0,                                                           &

                       ! Variable to track the sum of all errors
                       ! resulting from calls to allocate. Hence we
                       ! know that everything was successful if and
                       ! only if this is zero at the end

  temp_size,                                                                  &
  temp_tiles,                                                                 &
  temp_layers,                                                                &
                       ! For storing the size of array to allocate for variables
                       ! that are sometimes set to size 1.
                       ! Removes some duplicate allocate statements
  errcode
                       ! Variable to use in error report

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_CABLE_PROGS'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
errcode = 101

! CABLE vars to be initialized via JULES i/o
ALLOCATE( SoilTemp_CABLE(land_pts, ntype_cable, sm_levels), stat = error )
error_sum = error_sum + error

ALLOCATE( SoilMoisture_CABLE(land_pts, ntype_cable, sm_levels), stat = error )
error_sum = error_sum + error

ALLOCATE( FrozenSoilFrac_CABLE(land_pts, ntype_cable, sm_levels),             &
          stat = error )
error_sum = error_sum + error

ALLOCATE( SnowDepth_CABLE(land_pts, ntype_cable, nsmax), stat = error )
error_sum = error_sum + error

ALLOCATE( SnowMass_CABLE(land_pts, ntype_cable, nsmax), stat = error )
error_sum = error_sum + error

ALLOCATE( SnowTemp_CABLE(land_pts, ntype_cable, nsmax), stat = error )
error_sum = error_sum + error

ALLOCATE( SnowDensity_CABLE(land_pts, ntype_cable, nsmax), stat = error)
error_sum = error_sum + error

ALLOCATE( ThreeLayerSnowFlag_CABLE(land_pts, ntype_cable), stat = error)
error_sum = error_sum + error

ALLOCATE( OneLyrSnowDensity_CABLE(land_pts, ntype_cable), stat = error)
error_sum = error_sum + error

ALLOCATE( SnowAge_CABLE(land_pts, ntype_cable), stat = error )
error_sum = error_sum + error

!-----------------------------------------------------------------------
! Write out an error if there was one
!-----------------------------------------------------------------------
! Check for error.
IF ( error_sum /= 0 )                                                         &
  CALL ereport("allocate_cable_progs", errcode,                               &
               "Error allocating CABLE prognostic variable arrays")


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE allocate_cable_progs

END MODULE allocate_cable_arrays_mod
