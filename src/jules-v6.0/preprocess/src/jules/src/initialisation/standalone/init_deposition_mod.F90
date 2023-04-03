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
MODULE init_deposition_mod

!-----------------------------------------------------------------------------
! Description:
!   Reads the jules_deposition namelist items and checks them for consistency.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IMPLICIT NONE

! Private scope by default.
PRIVATE
PUBLIC init_deposition

CONTAINS

SUBROUTINE init_deposition(nml_dir)

USE io_constants, ONLY: namelist_unit

USE string_utils_mod, ONLY: to_string

USE jules_deposition_mod, ONLY:                                               &
! imported parameters
  dry_dep_ukca_jules,                                                         &
! imported scalars
  dry_dep_model, l_deposition, l_deposition_flux, l_ukca_ddep_lev1,           &
  ndry_dep_species, tundra_s_limit,                                           &
! imported procedure
  check_jules_deposition,                                                     &
! imported namelist
  jules_deposition

USE logging_mod, ONLY: log_info, log_fatal
  
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), INTENT(IN) ::                                               &
  nml_dir ! The directory containing the namelist files.

!-----------------------------------------------------------------------------
! Local parameters.
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER ::                                                &
  RoutineName = 'init_deposition' ! Name of this routine.

!-----------------------------------------------------------------------------
! Local variables.
!-----------------------------------------------------------------------------
INTEGER :: error

CHARACTER(LEN=errormessagelength) :: iomessage ! Error message.

!-----------------------------------------------------------------------------
! Open the namelist file.
!-----------------------------------------------------------------------------
OPEN(namelist_unit,                                                           &
     FILE=(TRIM(nml_dir) // '/' // 'jules_deposition.nml'),                   &
     STATUS='old', POSITION='rewind', ACTION='read', IOSTAT = error,          &
     IOMSG = iomessage)

IF ( error /= 0 ) THEN
  CALL log_fatal( RoutineName,                                                &
                  "Error opening namelist file "                    //        &
                  "jules_deposition.nml "                           //        &
                  "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //        &
                  TRIM(iomessage) // ")" )
END IF

!-----------------------------------------------------------------------------
! Read the namelist.
!-----------------------------------------------------------------------------
CALL log_info(RoutineName, "Reading JULES_DEPOSITION namelist...")

READ(namelist_unit, NML = jules_deposition, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 ) THEN
  CALL log_fatal( RoutineName,                                                &
                  "Error reading namelist JULES_DEPOSITION "        //        &
                  "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //        &
                  TRIM(iomessage) // ")" )
END IF

!-----------------------------------------------------------------------------
! Close namelist file.
!-----------------------------------------------------------------------------
CLOSE(namelist_unit, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 ) THEN
  CALL log_fatal( RoutineName,                                                &
                  "Error closing namelist file "                    //        &
                  "jules_deposition.nml "                           //        &
                  "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //        &
                  TRIM(iomessage) // ")" )
END IF

!-----------------------------------------------------------------------------
! Check that settings are reasonable.
!-----------------------------------------------------------------------------
CALL check_jules_deposition()

!-----------------------------------------------------------------------------
! Print summary information about the selected options.
!-----------------------------------------------------------------------------
IF ( l_deposition ) THEN
  CALL log_info( RoutineName, "Dry deposition is modelled." )

  ! Report which model is selected.
  SELECT CASE ( dry_dep_model )

  CASE ( dry_dep_ukca_jules )
    CALL log_info( RoutineName, "Dry deposition is modelled by JULES, " //    &
                                "using the approach of UKCA." )
    CALL log_info( RoutineName, "This code is still being developed and " //  &
                                "is not recommended for general use!" )

  CASE DEFAULT
    CALL log_fatal( RoutineName, "Invalid value for dry_dep_model." )

  END SELECT

  ! Report number of species to be modelled.
  CALL log_info( RoutineName, "Number of species for dry deposition = " //    &
                              to_string(ndry_dep_species) )
 
  !---------------------------------------------------------------------------
  ! Report scheme-specific values.
  !---------------------------------------------------------------------------
  SELECT CASE ( dry_dep_model )

  CASE ( dry_dep_ukca_jules )

    IF ( l_deposition_flux ) THEN
      IF ( l_ukca_ddep_lev1 ) THEN
        CALL log_info( RoutineName,                                           &
                       "Deposition only occurs from the lowest " //           &
                       "layer.")
      ELSE
        CALL log_info( RoutineName,                                           &
                       "Deposition occurs throughout the boundary layer." )
      END IF
    ELSE
      CALL log_info( RoutineName,                                             &
                     "Dry deposition velocities are modelled (not fluxes)." )
    END IF  !  l_deposition_flux

    CALL log_info( RoutineName, "tundra_s_limit = " //                        &
                   to_string(tundra_s_limit) )

  END SELECT  !  dry_dep_model 

END IF  !  l_deposition

RETURN
END SUBROUTINE init_deposition

END MODULE init_deposition_mod
