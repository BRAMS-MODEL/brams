#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
MODULE init_check_compatibility_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE init_check_compatibility()
!-----------------------------------------------------------------------------
! Description:
!   Checks that the enabled science schemes are compatible. Refer to the JULES
!   user manual for more information.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

USE jules_vegetation_mod,    ONLY: l_triffid, l_inferno
USE jules_irrig_mod,         ONLY: l_irrig_dmd
USE jules_soil_mod,          ONLY: l_tile_soil, l_holdwater
USE jules_radiation_mod,     ONLY: l_albedo_obs
USE switches_urban,          ONLY: l_moruses_storage
USE nvegparm,                ONLY: vf_nvg
USE jules_surface_types_mod, ONLY: urban_roof, npft

USE logging_mod,                   ONLY: log_error, log_fatal
USE errormessagelength_mod,        ONLY: errormessagelength
USE check_unavailable_options_mod, ONLY: check_unavailable_options

IMPLICIT NONE

! Work variables
INTEGER :: error  ! Error indicator
CHARACTER(LEN=errormessagelength) :: iomessage

CHARACTER(LEN=*), PARAMETER :: routinename='INIT_CHECK_COMPATIBILITY'

! Not all JULES options are available to standalone. Check that no unavailable
! options are being attempted to be used. A CASE statement using lsm_id can be
! added here and the routine split up if it transpires that some are required
! by child models.
CALL check_unavailable_options()

error = 0

! Check for compatibility with soil tiling. At present it is not allowed for
! TRIFFID, INFERNO, and l_albedo_obs

IF ( l_triffid .AND. l_tile_soil ) THEN
  error = 1
  CALL log_error(routinename,                                                 &
                 "TRIFFID not presently compatible with soil tiling")
END IF

IF ( l_inferno .AND. l_tile_soil ) THEN
  error = 1
  CALL log_error(routinename,                                                 &
                 "INFERNO not presently compatible with soil tiling")
END IF

IF ( l_albedo_obs .AND. l_tile_soil ) THEN
  error = 1
  CALL log_error(routinename,                                                 &
                 "l_albedo_obs not presently compatible with soil tiling")
END IF

IF (l_irrig_dmd .AND. l_holdwater) THEN
  error = 1
  CALL log_error(routinename,                                                 &
                 "l_holdwater=T not presently compatible with l_irrig_dmd=T")
END IF

IF ( l_moruses_storage ) THEN
  ! The roof surface needs to be either 0 (meaning "uncoupled" in this case) or
  ! 1 (radiatively coupled).
  IF ( vf_nvg(urban_roof - npft) /= 1.0 .AND.                                 &
     vf_nvg(urban_roof - npft) /= 0.0 )                                       &
     THEN
    error = 1
    WRITE(iomessage,'(A,F3.1)')                                               &
       "MORUSES roof coupling needs vf_nvg = 0 or 1: vf_nvg = ",              &
       vf_nvg(urban_roof - npft)
    CALL log_error(routinename, iomessage)
  END IF
END IF

IF ( error /= 0 )                                                             &
   CALL log_fatal(routinename,                                                &
                  "Error(s) found in enabled science options and/or " //      &
                  "parameter choice(s) - see earlier error message(s)")

RETURN
END SUBROUTINE init_check_compatibility

END MODULE init_check_compatibility_mod
#endif
