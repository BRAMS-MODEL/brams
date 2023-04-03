!******************************COPYRIGHT**************************************
! (c) UK Centre for Ecology & Hydrology.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237] 
!******************************COPYRIGHT**************************************
MODULE init_water_resources_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE init_water_resources(nml_dir)
!-----------------------------------------------------------------------------
! Description:
!   Reads in water resource namelist items and checks them for consistency.
!   Calls a further routine to set further values.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in HYDROLOGY
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

USE io_constants, ONLY: namelist_unit

USE string_utils_mod, ONLY: to_string

USE jules_water_resources_mod, ONLY:                                          &
! imported procedures
  check_jules_water_resources, set_jules_water_resources,                     &
! imported scalar variables
  l_prioritise, l_water_domestic, l_water_environment, l_water_industry,      &
  l_water_irrigation,l_water_livestock, l_water_resources,                    &
  l_water_transfers, nwater_use,                                              &
! imported array variables
  priority,                                                                   &
! imported namelist
  jules_water_resources

USE logging_mod, ONLY: log_info, log_fatal
  
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                         ! namelists

!-----------------------------------------------------------------------------
! Local scalar parameters.
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER ::                                                &
  RoutineName = 'INIT_WATER_RESOURCES'   ! Name of this routine.

!-----------------------------------------------------------------------------
! Local variables.
!-----------------------------------------------------------------------------
INTEGER :: error  ! Error indicator
INTEGER :: i      ! Loop counter.
CHARACTER(LEN=errormessagelength) :: iomessage

!-----------------------------------------------------------------------------
! Open the namelist file.
!-----------------------------------------------------------------------------
CALL log_info( RoutineName, "Reading JULES_WATER_RESOURCES namelist..." )

OPEN(namelist_unit,                                                           &
     FILE=(TRIM(nml_dir) // '/' // 'jules_water_resources.nml'),              &
     STATUS='old', POSITION='rewind', ACTION='read', IOSTAT = error,          &
     IOMSG = iomessage)

IF ( error /= 0 ) THEN
  CALL log_fatal( RoutineName,                                                &
                  "Error opening namelist file jules_water_resources.nml " // &
                  "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG="    //     &
                  TRIM(iomessage) // ")" )
END IF

!-----------------------------------------------------------------------------
! Read the namelist.
!-----------------------------------------------------------------------------
READ(namelist_unit, NML = jules_water_resources, IOSTAT = error,              &
     IOMSG = iomessage)
IF ( error /= 0 ) THEN
  CALL log_fatal( RoutineName,                                                &
                  "Error reading namelist JULES_WATER_RESOURCES "   //        &
                  "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //        &
                  TRIM(iomessage) // ")" )
END IF

!-----------------------------------------------------------------------------
! Close namelist file.
!-----------------------------------------------------------------------------
CLOSE(namelist_unit, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 ) THEN
  CALL log_fatal( RoutineName,                                                &
                  "Error closing namelist file jules_water_resources.nml " // &
                  "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG="   //      &
                  TRIM(iomessage) // ")" )
END IF

!-----------------------------------------------------------------------------
! Check that settings are reasonable.
!-----------------------------------------------------------------------------
CALL check_jules_water_resources

!-----------------------------------------------------------------------------
! Set further values.
!-----------------------------------------------------------------------------
CALL set_jules_water_resources

!-----------------------------------------------------------------------------
! Print some human friendly summary information about the selected options.
!-----------------------------------------------------------------------------

IF ( l_water_resources ) THEN

  CALL log_info( RoutineName, "Water resources are modelled." )

  IF ( l_water_domestic ) THEN
    CALL log_info( RoutineName, "Domestic demand is included." )
  END IF
  IF ( l_water_environment ) THEN
    CALL log_info( RoutineName, "Environmental demand is included." )
  END IF
  IF ( l_water_industry ) THEN
    CALL log_info( RoutineName, "Industrial demand is included." )
  END IF
  IF ( l_water_livestock ) THEN
    CALL log_info( RoutineName, "Livestock demand is included." )
  END IF
  IF ( l_water_transfers ) THEN
    CALL log_info( RoutineName, "(Explicit) water transfers are included." )
  END IF

  IF ( l_prioritise ) THEN
    CALL log_info( RoutineName, "Demands are prioritised:" )
    DO i = 1,nwater_use
      CALL log_info( RoutineName, "  Priority " // TRIM(to_string(i)) // ": " &
                                  // priority(i) )
    END DO
  ELSE
    CALL log_info( RoutineName, "Demands are not prioritised." )
  END IF

ELSE

  CALL log_info( RoutineName, "Water resources are NOT modelled." )

END IF  !  l_water_resources

RETURN
END SUBROUTINE init_water_resources

END MODULE init_water_resources_mod
