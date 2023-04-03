#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE internal_next_output_file(profile)

USE datetime_mod, ONLY: period_year, period_month, datetime_advance

USE templating_mod, ONLY: tpl_yr_4digit, tpl_mon_2digit

USE model_time_mod, ONLY: is_spinup, spinup_start, spinup_end, spinup_cycle

USE file_ts_mod, ONLY: file_ts_close

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Opens the next output file for the profile
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(output_profile), INTENT(INOUT) :: profile  ! The profile to open the next
                                                ! file for

! Work variables
TYPE(datetime) :: file_start, file_end  ! The start and end time for the next
                                        ! file

LOGICAL :: use_template  ! Indicates if we are using a file name template or not
CHARACTER(LEN=max_file_name_len) :: file_name
                                        ! The name/time template for the file



!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Close the current file if there is one
!-----------------------------------------------------------------------------
IF ( profile%has_open_file ) CALL file_ts_close(profile%fh)


!-----------------------------------------------------------------------------
! Work out when the data in the file will start and end
!
! We tell the file that data will end just after output_end, since this will
! allow us to output data exactly at output_end, but no later
!-----------------------------------------------------------------------------
IF ( is_spinup ) THEN
  file_start = spinup_start
  file_end   = datetime_advance(spinup_end, 1)
ELSE
  file_start = profile%output_start
  file_end   = datetime_advance(profile%output_end, 1)
END IF


!-----------------------------------------------------------------------------
! Work out what file name/template we want to use
!-----------------------------------------------------------------------------
! All file names start with the run id followed by the profile name
file_name = TRIM(run_id) // "." // TRIM(profile%profile_name)

IF ( is_spinup ) THEN
  ! If we are in a spinup cycle, add that to the file name
  file_name = TRIM(file_name) // "." // "spin" // TRIM(to_string(spinup_cycle))
END IF

! Indicate whether we will be using time templating or not and add the required
! time template to the file name
use_template = .FALSE.

IF ( profile%file_period == period_year .OR.                                  &
     profile%file_period == period_month ) THEN
  use_template = .TRUE.

  file_name = TRIM(file_name) // "." // tpl_yr_4digit
  IF ( profile%file_period == period_month )                                  &
    file_name = TRIM(file_name) // tpl_mon_2digit
END IF

! Add an extension based on the output format
SELECT CASE ( output_format )
CASE ( format_ascii )
  file_name = TRIM(file_name) // ".asc"

CASE ( format_ncdf )
  file_name = TRIM(file_name) // ".nc"

CASE DEFAULT
  CALL log_fatal("internal_next_output_file",                                 &
                 "Unrecognised file format - " // TRIM(output_format))
END SELECT

! Prepend the output directory
file_name = TRIM(output_dir) // "/" // TRIM(file_name)


!-----------------------------------------------------------------------------
! Open a new file using the properties we have gathered
!-----------------------------------------------------------------------------
profile%fh = internal_open_output_file(                                       &
  file_start, file_end, profile%output_period, use_template, file_name,       &
  profile%fields(1:profile%nfields)                                           &
)


!-----------------------------------------------------------------------------
! Indicate that the profile now has an open file
!-----------------------------------------------------------------------------
profile%has_open_file = .TRUE.

RETURN

END SUBROUTINE internal_next_output_file
#endif
