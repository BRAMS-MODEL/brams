! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE jules_model_environment_mod

!-----------------------------------------------------------------------------
! Description:
!   Contains model interface options. Not all JULES options are available to all
!   environments in which JULES is run e.g. standalone, UM, CABLE, LIS, MONC.
!-----------------------------------------------------------------------------

USE missing_data_mod, ONLY: imdi

IMPLICIT NONE

PRIVATE
! Please do not make this namelist public. The only place this namelist
! should be used is in this module. Please see the JULES Working Practices.
! If required new items added can be made public, but as a default it should be
! private; l_jules_parent should not have science options attached.
PUBLIC :: print_nlist_jules_model_environment,                                &
          read_nml_jules_model_environment,                                   &
          check_jules_model_environment


INTEGER :: l_jules_parent = imdi ! Switch to identify UM-JULES environment
INTEGER, PUBLIC :: lsm_id = imdi ! Switch to identify land surface model
INTEGER, PARAMETER, PUBLIC :: jules = 1
INTEGER, PARAMETER, PUBLIC :: cable = 2
CHARACTER(LEN=5), DIMENSION(2) :: lsm_name = (/ 'JULES', 'CABLE' /)

NAMELIST  / jules_model_environment / l_jules_parent, lsm_id

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='JULES_MODEL_ENVIRONMENT_MOD'

CONTAINS

SUBROUTINE check_jules_model_environment()

!-----------------------------------------------------------------------------
! Description:
!   Checks JULES_MODEL_ENVIRONMENT namelist for consistency
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

USE ereport_mod, ONLY: ereport
USE jules_print_mgr, ONLY: jules_message, jules_print

IMPLICIT NONE

! Options for defining JULES parent models
INTEGER, PARAMETER ::                                                         &
  jules_standalone = 0,                                                       &
  um_jules         = 1

INTEGER :: errcode   ! error code to pass to ereport.
CHARACTER(LEN=*), PARAMETER :: RoutineName='CHECK_JULES_MODEL_ENVIRONMENT'

!-----------------------------------------------------------------------------

! Check that l_jules_parent is consistent with UM_JULES ifdef.
! MONC_JULES (jules:#347) should be added here when available.
IF ( l_jules_parent /= jules_standalone ) THEN
  errcode = 20
  WRITE(jules_message,'(A,I0)')                                               &
     "l_jules_parent should be 0 for standalone. l_jules_parent = ",          &
     l_jules_parent
  CALL ereport(RoutineName, errcode, jules_message )
END IF

IF ( l_jules_parent == imdi ) THEN
  errcode = 30
  WRITE(jules_message,'(A,I0)') "l_jules_parent is currently unset."
  CALL ereport(RoutineName, errcode, jules_message )
END IF

!-----------------------------------------------------------------------------
! Check if land surface model id is valid
!-----------------------------------------------------------------------------
SELECT CASE ( lsm_id )
CASE ( cable )
  ! Allowed LSM
CASE ( jules )
  ! Allowed LSM
CASE DEFAULT
  errcode = 40
  WRITE(jules_message,'(A,I0)') "Invalid LSM ID provided. lsm_id = ", lsm_id
  CALL ereport(RoutineName, errcode, jules_message )
END SELECT
WRITE(jules_message,'(A)') "Land surface model selected is " // lsm_name(lsm_id)
CALL jules_print(RoutineName, jules_message)


END SUBROUTINE check_jules_model_environment

!-----------------------------------------------------------------------------

SUBROUTINE print_nlist_jules_model_environment()

USE jules_print_mgr, ONLY: jules_print

IMPLICIT NONE

CHARACTER(LEN=50000)        :: lineBuffer
CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_JULES_MODEL_ENVIRONMENT'

CALL jules_print(RoutineName, 'Contents of namelist jules_model_environment')

WRITE(lineBuffer, '(A,I0)') 'l_jules_parent = ', l_jules_parent
CALL jules_print(RoutineName, lineBuffer)

WRITE(lineBuffer, '(A,I0)') 'lsm_id = ', lsm_id
CALL jules_print(RoutineName, lineBuffer)

CALL jules_print(RoutineName, '- - - - - - end of namelist - - - - - -')

END SUBROUTINE print_nlist_jules_model_environment

!-----------------------------------------------------------------------------


SUBROUTINE read_nml_jules_model_environment(nml_dir)

USE io_constants, ONLY: namelist_unit

USE string_utils_mod, ONLY: to_string

USE logging_mod, ONLY: log_info, log_fatal

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Arguments
CHARACTER(LEN=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                         ! namelists
! Work variables
INTEGER :: error  ! Error indicator
CHARACTER(LEN=errormessagelength) :: iomessage
CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_JULES_MODEL_ENVIRONMENT'

CALL log_info(RoutineName, "Reading JULES_MODEL_ENVIRONMENT namelist...")

OPEN(namelist_unit,                                                           &
     FILE=(TRIM(nml_dir) // '/' // 'model_environment.nml'),                  &
     STATUS='old', POSITION='rewind', ACTION='read', IOSTAT = error,          &
     IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal(RoutineName,                                                 &
                 "Error opening namelist file model_environment.nml " //      &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

READ(namelist_unit, NML = jules_model_environment, IOSTAT = error,            &
   IOMSG = iomessage)
IF ( error /= 0 )                                                             &
   CALL log_fatal(RoutineName,                                                &
                 "Error reading namelist JULES_MODEL_ENVIRONMENT " //         &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

CLOSE(namelist_unit, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
   CALL log_fatal(RoutineName,                                                &
                 "Error closing namelist file model_environment.nml " //      &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

RETURN
END SUBROUTINE read_nml_jules_model_environment

END MODULE jules_model_environment_mod
