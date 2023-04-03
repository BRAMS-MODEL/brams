! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE output_mod

USE io_constants, ONLY: max_sdf_name_len, max_file_name_len,                  &
                         format_len, format_ascii, format_ncdf

USE string_utils_mod, ONLY: to_string

USE datetime_mod, ONLY: datetime,                                             &
! Also import the comparison operators on the datetime type
                           operator( == ), operator( /= ), operator( < ),     &
                           operator( > ), operator( <= ), operator( >= )

USE grid_utils_mod, ONLY: grid_info, subgrid_info

USE data_cube_mod, ONLY: data_cube

USE file_ts_mod, ONLY: file_ts

USE logging_mod, ONLY: log_info, log_warn, log_fatal

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------------
INTEGER, PARAMETER :: nprofiles_max = 50

! Constants for dump period unit
CHARACTER(LEN=1), PARAMETER ::                                                &
  ! Calendar year, default
  dump_period_year  = 'Y',                                                    &
  ! Calendar second of day
  dump_period_time  = 'T'

! Constants for types of output available
CHARACTER(LEN=1), PARAMETER ::                                                &
  output_snapshot = 'S',                                                      &
      ! Indicates output should be a snapshot every output period
  output_accum    = 'A',                                                      &
      ! Indicates output should be an accumulation over the output period
  output_mean     = 'M',                                                      &
      ! Indicates output should be a mean over the output period
  output_max      = 'X',                                                      &
      ! Indicates output should be a max over the output period    
  output_min      = 'N'
      ! Indicates output should be a min over the output period    

CHARACTER(LEN=1), PARAMETER, DIMENSION(5) :: allowed_output_types = (/        &
  output_snapshot, output_accum, output_mean, output_max, output_min          &
/)
 
! Dimension names for the output grid
! USED ONLY FOR 1D GRID
CHARACTER(LEN=max_sdf_name_len), PARAMETER :: grid_dim_name = "points"
    ! The dimension name of the single grid dimension

! USED ONLY FOR 2D GRID
CHARACTER(LEN=max_sdf_name_len), PARAMETER ::                                 &
  x_dim_name = "x",                                                           &
  y_dim_name = "y"

CHARACTER(LEN=max_sdf_name_len) :: time_dim_name = "time"
                       ! The dimension name for the time dimension

CHARACTER(LEN=format_len), PARAMETER :: output_format = format_ncdf
                                         ! The file format to use for output


!-----------------------------------------------------------------------------
! Module types
!-----------------------------------------------------------------------------
TYPE output_profile

  CHARACTER(LEN=max_sdf_name_len) :: profile_name  ! The name of the profile

  LOGICAL :: has_open_file = .FALSE.  ! Indicates if the profile is currently
                                      ! open
  TYPE(file_ts) :: fh  ! Handle to the underlying file

  LOGICAL :: output_initial ! T - this profile should output initial data
                            !     for each section it is outputting
                            ! F - this profile should not output initial data

  LOGICAL :: output_spinup  ! T - this profile is outputting during spinup
                            ! F - this profile is not outputting during spinup

  LOGICAL :: output_main_run  ! T - this profile is outputting (part of) the
                              !     main run
                              ! F - this profile is not outputting any of
                              !     the main run
  TYPE(datetime) :: output_start, output_end
                            ! The start and end times for output during the
                            ! main run

  INTEGER :: file_period    ! The period to use for files for this profile
  INTEGER :: output_period  ! The output period for this profile
  INTEGER :: sample_period  ! The sampling period for this profile (s).
                            ! Note that "special" periods, such as monthly,
                            ! are not allowed because time accumulations
                            ! need to know how many timesteps went into the
                            ! sample, which is less straightforward for
                            ! special periods (which generally don't have
                            ! constant length; a further counter would
                            ! likely be required to deal with that case).

  TYPE(datetime) :: next_sample_time  ! The time that the data will next be
                                      ! sampled
  TYPE(datetime) :: current_output_time  ! The time that the current output
                                         ! period started
  TYPE(datetime) :: next_output_time  ! The time that the next output period
                                      ! starts

  INTEGER :: samples_in_period  ! The number of times that data has been
                                ! sampled in this output period

  INTEGER :: nfields = 0  ! The number of fields in the profile
  TYPE(output_field), POINTER :: fields(:) => NULL()  ! The fields in the profile

END TYPE output_profile


TYPE output_field

  INTEGER :: var_id  ! The variable id of the model variable that
                     ! this field outputs - see model_interface_mod

  CHARACTER(LEN=max_sdf_name_len) :: output_name
                     ! The name of the variable in output files

  CHARACTER(LEN=1) :: field_type  ! The type of output to use

  INTEGER :: file_id  ! The id of the variable in the underlying file

  TYPE(data_cube) :: field_data ! The data collected for output as a cube

END TYPE output_field


!-----------------------------------------------------------------------------
! Module variables
!-----------------------------------------------------------------------------
CHARACTER(LEN=max_file_name_len) :: output_dir = "./output"
                                         ! The directory to use for output

CHARACTER(LEN=max_sdf_name_len) :: run_id = ""  ! Identifier for the run

TYPE(grid_info), SAVE :: grid  ! The output grid
                               ! gfortran requires the SAVE attribute
                               ! Other compilers don't care either way

LOGICAL :: use_subgrid = .FALSE.
    ! T => the model grid is a subset of the output grid
    ! F => the model grid is the output grid

TYPE(subgrid_info), SAVE :: subgrid  ! Used if use_subgrid=T
                                     ! Information about the subgrid to write
                                     ! gfortran requires the SAVE attribute
                                     ! Other compilers don't care either way


INTEGER :: nprofiles = 0  ! The number of output profiles requested
TYPE(output_profile), SAVE :: profiles(nprofiles_max)  ! The requested profiles

INTEGER :: dump_period = 1
!DSM CHARACTER(LEN=1) :: dump_period_unit = dump_period_year
CHARACTER(LEN=1) :: dump_period_unit = dump_period_time
!   The dump period, and its unit (calendar year or second of calendar day)
!   In calendar year mode, the number of years between model dumps.
!   Note that the calendar year (date) is used to determine whether a dump
!   is written, not the number of years simulated so far.
!   In second of calendar day mode, the number of seconds between model dumps.
!   Note that this is calculated using the number of seconds into the day,
!   not the number seconds simulated so far.
!   Dumps are also written at the start and end of the main run, and at the
!   start of each cycle of any spin up.

!-----------------------------------------------------------------------------
! Visibility declarations
!-----------------------------------------------------------------------------
PRIVATE
PUBLIC                                                                        &
! Constants
    dump_period_year, dump_period_time,                                       &
! Variables
    output_dir, run_id, grid, use_subgrid, subgrid, dump_period,              &
    dump_period_unit,                                                         &
! Routines
    register_output_profile, sample_data, output_data, close_all,             &
    output_initial_data


CONTAINS


! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE output_initial_data()

USE datetime_mod, ONLY: datetime_subtract

USE data_cube_mod, ONLY: cube_free

USE model_time_mod, ONLY: is_spinup, spinup_start, main_run_start,            &
                           spinup_cycle, current_time

USE file_ts_mod, ONLY: file_ts_write_var, file_ts_close

USE model_interface_mod, ONLY: extract_var

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   For each profile, decide if we need to output a file containing initial
!   data for a section of the run
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
TYPE(file_ts) :: FILE  ! The initial data file

CHARACTER(LEN=max_file_name_len) :: file_name
                                        ! The name for the initial data file

TYPE(datetime) :: file_start  ! The start time for the file

TYPE(data_cube) :: data_for_output  ! Cube to hold data for output

INTEGER :: i,j  ! Loop counter


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! If we are not in the first timestep of a section, there is nothing to do
! for any profile
!-----------------------------------------------------------------------------
IF ( is_spinup .AND. current_time /= spinup_start )                           &
  RETURN

IF ( .NOT. is_spinup .AND. current_time /= main_run_start )                   &
  RETURN


!-----------------------------------------------------------------------------
! Check what we need to do for each profile
!-----------------------------------------------------------------------------
DO i = 1,nprofiles

  !-----------------------------------------------------------------------------
  ! Check if we need to output initial data for this profile
  !-----------------------------------------------------------------------------
  IF ( .NOT. profiles(i)%output_initial )                                     &
  ! If initial state was not asked for, there is nothing to do
        CYCLE

  IF ( .NOT. ANY( profiles(i)%fields(1:profiles(i)%nfields)%field_type        &
                  == output_snapshot ) ) THEN
    ! If there is no snapshot data in the profile, there is nothing to output here
    CYCLE
  END IF

  IF ( is_spinup .AND. .NOT. profiles(i)%output_spinup )                      &
  ! If we are at the start of a spinup cycle but profile is not outputting during
  ! spinup, we don't need to do anything
        CYCLE

  IF ( .NOT. is_spinup ) THEN
    ! We are at the start of the main run
    IF ( .NOT. profiles(i)%output_main_run )                                  &
    ! If the profile is not outputting for the main run, there is nothing to do
            CYCLE

    IF ( profiles(i)%output_start /= main_run_start )                         &
    ! If the profile does not start outputting until later in the run, there is nothing to do
            CYCLE
  END IF


  !-----------------------------------------------------------------------------
  ! Open an initial data file for the profile
  !-----------------------------------------------------------------------------
  ! All file names start with the run id followed by the profile name
  file_name = TRIM(run_id) // "." // TRIM(profiles(i)%profile_name)

  IF ( is_spinup ) THEN
    ! If we are in a spinup cycle, add that to the file name
    file_name = TRIM(file_name) // "." // "spin" // TRIM(to_string(spinup_cycle))
  END IF

  ! Indicate that it is an initial data file
  file_name = TRIM(file_name) // ".initial"

  ! Add an extension based on the output format
  SELECT CASE ( output_format )
  CASE ( format_ascii )
    file_name = TRIM(file_name) // ".asc"

  CASE ( format_ncdf )
    file_name = TRIM(file_name) // ".nc"

  CASE DEFAULT
    CALL log_fatal("output_initial_data",                                     &
                   "Unrecognised file format - " // TRIM(output_format))
  END SELECT

  ! Prepend the output directory
  file_name = TRIM(output_dir) // "/" // TRIM(file_name)

  ! The timestamp in files is the end of the output period, so we create a file
  ! with exactly one output period ending at the current time
  ! We don't much care about time_bounds since the only real valid output is
  ! snapshot variables, for which the timestamp is the important thing
  ! So we create a file that ends at the current time, starts 1 second before
  ! and has an output period of 1
  file_start = datetime_subtract(current_time, 1)

  FILE=internal_open_output_file(                                             &
    file_start, current_time, 1, .FALSE., file_name,                          &
    profiles(i)%fields(1:profiles(i)%nfields)                                 &
  )


  !-----------------------------------------------------------------------------
  ! Populate any snapshot variables in the file
  !-----------------------------------------------------------------------------
  DO j = 1,profiles(i)%nfields
    ! Ignore variables that are not snapshot variables
    IF ( profiles(i)%fields(j)%field_type /= output_snapshot ) CYCLE

    ! Extract the data to output
    data_for_output = extract_var(profiles(i)%fields(j)%var_id)

    ! Output the data to file
    CALL file_ts_write_var(                                                   &
      FILE, profiles(i)%fields(j)%file_id, data_for_output,                   &
    ! Subgrid information (for writing a slab of the output grid in parallel mode)
            use_subgrid, subgrid                                              &
          )

    ! Free the temporary cube
    CALL cube_free(data_for_output)
  END DO


  !-----------------------------------------------------------------------------
  ! Close the file
  !-----------------------------------------------------------------------------
  CALL file_ts_close(FILE)

END DO

RETURN

END SUBROUTINE output_initial_data
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE register_output_profile(profile_name, output_initial,              &
                                   output_spinup, output_main_run,            &
                                   output_start, output_end,                  &
                                   output_period, sample_period, file_period, &
                                   identifiers, var_names, output_types)

USE datetime_mod, ONLY: period_month, period_year, secs_in_day,               &
                         datetime_to_string

USE model_time_mod, ONLY: max_spinup_cycles, main_run_start, main_run_end,    &
                           timestep_len

USE model_interface_mod, ONLY: get_var_id

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Register an output profile to output given variables over the given
!   period
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
CHARACTER(LEN=*), INTENT(IN) :: profile_name
    ! The name of this profile - only used in filenames

LOGICAL, INTENT(IN) :: output_initial
    ! T - this profile should output initial data
    ! F - this profile should not output initial data

LOGICAL, INTENT(IN) :: output_spinup
    ! T - profile will provide output during spinup
    ! F - profile will not provide output during spinup

LOGICAL, INTENT(IN) :: output_main_run
    ! T - profile will provide output for the specified portion of the
    !     main run only
    ! F - profile will not provide any output during the main run

TYPE(datetime), INTENT(IN) :: output_start, output_end
    ! USED ONLY IF output_main_run=T
    ! Start and end dates for output if providing output for main run

INTEGER, INTENT(IN) ::                                                        &
  output_period,                                                              &
    ! The output period for this profile - this should be a multiple of
    ! the model timestep or a 'special' period (monthly/yearly)
  sample_period,                                                              &
    ! The sampling period for this profile (s). This should be a multiple of
    ! the model timestep length and a factor of output_period.
  file_period
    ! This should be a special period for monthly or yearly files
    ! Any other value will result in a single file for the whole of the
    ! output period

CHARACTER(LEN=*), INTENT(IN) :: identifiers(:)
    ! The model identifiers of the variables this profile provides output for

CHARACTER(LEN=*), INTENT(IN) :: var_names(:)
    ! The name to use in output files for each variable
    ! Each value can be the empty string - if it is, the model identifier will
    ! be used

CHARACTER(LEN=*), INTENT(IN) :: output_types(:)
    ! The type of output to use for each variable


! Work variables
INTEGER :: nvars  ! The number of output fields in the profile
INTEGER :: test_length  ! A length of time used in testing (s).

! Local versions of INTENT(IN) arguments that we can modify
LOGICAL :: output_spinup_local
TYPE(datetime) :: output_start_local, output_end_local
CHARACTER(LEN=len(var_names)) :: var_names_local(SIZE(var_names))

INTEGER :: i  ! Loop counters


!-----------------------------------------------------------------------------

! Check that we have space to register another profile
IF ( nprofiles >= nprofiles_max )                                             &
  CALL log_fatal("register_output_profile",                                   &
                 "Too many profiles registered - try increasing nprofiles_max")

! Check that we have an output type for each variable
nvars = SIZE(identifiers)

IF ( nvars /= SIZE(var_names) )                                               &
  CALL log_fatal("register_output_profile",                                   &
                 "identifiers and var_names must have the same " //           &
                 "number of elements")

IF ( nvars /= SIZE(output_types) )                                            &
  CALL log_fatal("register_output_profile",                                   &
                 "identifiers and output_types must have the same " //        &
                 "number of elements")

!-----------------------------------------------------------------------------
! Check that arguments make sense
!-----------------------------------------------------------------------------
! First copy INTENT(IN) arguments that we might want to modify into local
! versions
output_spinup_local = output_spinup
output_start_local  = output_start
output_end_local    = output_end

! Check that we have been given a profile name
IF ( LEN_TRIM(profile_name) == 0 )                                            &
  CALL log_fatal("register_output_profile", "No profile name given")

! Check that this profile name has not already been used.
IF ( ANY( profiles(1:nprofiles)%profile_name == profile_name ) ) THEN
  CALL log_fatal("register_output_profile",                                   &
                 "Duplicate profile name: " // TRIM(profile_name))
END IF

IF ( output_spinup_local .AND. max_spinup_cycles <= 0 ) THEN
  ! If spinup output has been requested but there is no spinup, issue a warning
  ! and ignore it
  CALL log_warn("register_output_profile",                                    &
                 "Model has no spinup - ignoring request for output " //      &
                 "during spinup from profile " // TRIM(profile_name))
  output_spinup_local = .FALSE.
END IF

! Check if the above means that this profile has requested no output at all
IF ( .NOT. output_spinup_local .AND. .NOT. output_main_run ) THEN
  CALL log_warn("register_output_profile",                                    &
                 "Profile " // TRIM(profile_name) // " will provide no " //   &
                 "output with current model setup - ignoring")
  RETURN
END IF

IF ( output_main_run ) THEN
  ! If output has been requested for times that are outside the main run, we
  ! issue a warning and truncate
  IF ( output_start_local < main_run_start ) THEN
    CALL log_warn("register_output_profile",                                  &
                   "Output has been requested for times before the " //       &
                   "start of the main run - output will start at the " //     &
                   "start of the main run")
    output_start_local = main_run_start
  END IF

  IF ( output_end_local > main_run_end ) THEN
    CALL log_warn("register_output_profile",                                  &
                   "Output has been requested for times after the " //        &
                   "end of the main run - output will end at the " //         &
                   "end of the main run")
    output_end_local = main_run_end
  END IF

  ! Check that the times make sense
  IF ( output_end_local <= output_start_local )                               &
    CALL log_fatal("register_output_profile",                                 &
                   "Output cannot end before it has started." //              &
                   " Profile: " // TRIM(profile_name) )
END IF

! Check that the given periods make sense
IF ( file_period /= period_month .AND. file_period /= period_year )           &
! Warn that one file will be used for all output, since file_period is not
! a special period
    CALL log_info("register_output_profile",                                  &
                  "Since file_period is not a 'special' period, all " //      &
                  "output from each section will go into one file")

IF ( output_period /= period_month .AND. output_period /= period_year .AND.   &
     MOD(output_period, timestep_len) /= 0 )                                  &
! If output period is not a special period, then it must be a multiple of
! the model timestep
    CALL log_fatal("register_output_profile",                                 &
                   "Output period must be a 'special' period or a " //        &
                   "multiple of model timestep." //                           &
                   " Profile: " // TRIM(profile_name) )

IF ( sample_period < 1 .OR. MOD(sample_period, timestep_len) /= 0 )           &
  CALL log_fatal("register_output_profile",                                   &
                 "Sample period must be a multiple of model timestep." //     &
                 " Profile: " // TRIM(profile_name) )
! Insist that output_period is a multiple of sample_period. This test is
! designed to check cases when sample_period > timestep_len, but can always
! be applied. It is also more important for snapshot outputs which otherwise
! can be timestamped with a time other than the actual time of the data,
! which is confusing (whereas for other output types we anyway accept that
! any intermittent sampling only gives approximate results), but we impose
! this requirement for all output types for convenience and because the
! outputs will generally be better behaved with this restriction.
IF ( output_period < 1 ) THEN
  ! This is a "special" period, such as monthly, which is a whole number of
  ! days.
  ! We will check that one day is a multiple of sample_period.
  test_length = secs_in_day
ELSE
  ! We will check that output_period is a multiple of sample_period.
  test_length = output_period
END IF
IF ( MOD( test_length, sample_period ) /= 0 ) THEN
  CALL log_fatal("register_output_profile",                                   &
                 "Output period must be a multiple of sampling period. " //   &
                 "Profile: " // TRIM(profile_name) )
END IF

!-----------------------------------------------------------------------------
! Indicate what output we are providing
!-----------------------------------------------------------------------------
IF ( output_initial )                                                         &
  CALL log_info("register_output_profile",                                    &
                "Profile with name " // TRIM(profile_name) // " " //          &
                "will provide initial data for each section it is " //        &
                "outputting for")

IF ( output_spinup_local )                                                    &
  CALL log_info("register_output_profile",                                    &
                "Profile with name " // TRIM(profile_name) // " " //          &
                "registered to provide output during spinup")

IF ( output_main_run )                                                        &
  CALL log_info("register_output_profile",                                    &
                "Profile with name " // TRIM(profile_name) // " " //          &
                "registered to provide output for main run from " //          &
                datetime_to_string(output_start_local) // " to " //           &
                datetime_to_string(output_end_local))

!-----------------------------------------------------------------------------
! Set up the variable names that will be used for output
! Use user-supplied output name if given, otherwise use the model identifier
!
! We also check for duplicates in the output variable names as we go
! Note that duplicate identifiers are allowed, as long different names have
! been specified for them to use in output files
!-----------------------------------------------------------------------------
var_names_local(:) = var_names(:)
DO i = 1,nvars

  ! Check we don't have an empty identifier.
  IF ( LEN_TRIM(identifiers(i)) == 0 ) THEN
    CALL log_fatal("register_output_profile",                                 &
                   "Missing variable name for variable #" // to_string(i) //  &
                   " detected in output profile " // TRIM(profile_name))
  END IF

  IF ( LEN_TRIM(var_names(i)) <= 0 ) var_names_local(i) = identifiers(i)

  ! If the output name being processed matches any of those previously processed,
  ! issue a fatal error
  ! We issue a fatal error to force the user to resolve the problem and avoid any
  ! confusion over what has been output
  IF ( ANY(var_names_local(1:i-1) == var_names_local(i)) )                    &
    CALL log_fatal("register_output_profile",                                 &
                   "Duplicate variable name=" // TRIM(var_names_local(i)) //  &
                   " for output files " //                                    &
                   "detected in profile " // TRIM(profile_name))

  ! Check that only the allowed types have been given as output types
  IF ( .NOT. ANY(output_types(i) == allowed_output_types) ) THEN
    CALL log_fatal("register_output_profile",                                 &
                   "Unrecognised output type '" // TRIM(output_types(i)) //   &
                   "' in profile " // TRIM(profile_name) )
  END IF
END DO

!-----------------------------------------------------------------------------
! Set up the output_profile object and its corresponding output_field objects
!-----------------------------------------------------------------------------
! Store profile constant data
nprofiles = nprofiles + 1
profiles(nprofiles)%profile_name    = profile_name
profiles(nprofiles)%output_initial  = output_initial
profiles(nprofiles)%output_spinup   = output_spinup_local
profiles(nprofiles)%output_main_run = output_main_run
profiles(nprofiles)%file_period     = file_period
profiles(nprofiles)%output_period   = output_period
profiles(nprofiles)%sample_period   = sample_period
profiles(nprofiles)%output_start    = output_start_local
profiles(nprofiles)%output_end      = output_end_local

! Allocate space on the output_profile object for the fields
profiles(nprofiles)%nfields = nvars
ALLOCATE(profiles(nprofiles)%fields(nvars))

DO i = 1,nvars
  ! Store as much info about fields as we currently know for later
  profiles(nprofiles)%fields(i)%var_id      = get_var_id(identifiers(i))
  profiles(nprofiles)%fields(i)%output_name = var_names_local(i)
  profiles(nprofiles)%fields(i)%field_type  = output_types(i)
END DO

RETURN

END SUBROUTINE register_output_profile
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE sample_data()

USE datetime_mod, ONLY: datetime_advance

USE data_cube_mod, ONLY: data_cube, operator (+), cube_safe_copy,             &
                          cube_free, cube_min, cube_max

USE model_time_mod, ONLY: current_time, is_spinup, spinup_start

USE model_interface_mod, ONLY: extract_var

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   For each profile, decide if we need to sample data
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
TYPE(data_cube) :: var_data, accum
               ! Workspace cubes that can be deallocated to avoid memory leaks

INTEGER :: i,j  ! Loop counter

!-----------------------------------------------------------------------------

DO i = 1,nprofiles
  !-----------------------------------------------------------------------------
  ! First check if the profile is currently active
  !-----------------------------------------------------------------------------
  IF ( is_spinup .AND. .NOT. profiles(i)%output_spinup )                      &
  ! If we are in spinup but profile is not outputting, we don't need to do anything
        CYCLE

  IF ( .NOT. is_spinup .AND. ( .NOT. profiles(i)%output_main_run .OR.         &
       ( current_time < profiles(i)%output_start .OR.                         &
         profiles(i)%output_end < current_time ) ) )                          &
  ! If we are in the main run but either not outputting any of the main run or
  ! the time is out of the range of times the profile is outputting, we don't
  ! need to do anything
        CYCLE


  !-----------------------------------------------------------------------------
  ! Initialise the profile fields, if required
  !
  ! This is the case if:
  !
  !   - The profile is outputting during spinup and we are at in the first
  !     timestep of a spinup cycle
  !
  !   - The profile is outputting part of the main run and we are in the
  !     timestep in the main run when output starts
  !
  ! We do this check in two stages for readability
  !-----------------------------------------------------------------------------
  IF ( is_spinup .AND. profiles(i)%output_spinup .AND.                        &
       current_time == spinup_start )                                         &
    CALL internal_init_profile_vars(profiles(i))

  IF ( .NOT. is_spinup .AND. profiles(i)%output_main_run .AND.                &
       current_time == profiles(i)%output_start )                             &
    CALL internal_init_profile_vars(profiles(i))


  !-----------------------------------------------------------------------------
  ! If we are in the first timestep of a sample period (i.e. current time is
  ! now equal to (or greater than) the next sample time), then sample the data
  !-----------------------------------------------------------------------------
  IF ( profiles(i)%next_sample_time <= current_time ) THEN
    ! If we are in the first timestep of an output period, capture snapshot
    ! values
    DO j = 1,profiles(i)%nfields
      SELECT CASE ( profiles(i)%fields(j)%field_type )
      CASE ( output_snapshot )
        ! Free the previous cube
        CALL cube_free(profiles(i)%fields(j)%field_data)
        ! Extract the data from the variable and set it as the current cube
        profiles(i)%fields(j)%field_data = extract_var(profiles(i)%fields(j)%var_id)

      CASE ( output_accum, output_mean )
        ! Collect the accumulated value in the work cube
        var_data = extract_var(profiles(i)%fields(j)%var_id)
        accum = profiles(i)%fields(j)%field_data + var_data
        ! Copy the accumulated value safely into the data cube
        CALL cube_safe_copy(profiles(i)%fields(j)%field_data, accum)
        ! Free the work cubes
        CALL cube_free(var_data)
        CALL cube_free(accum)
            
      CASE ( output_min, output_max )
        ! If this is the first sample of a new output period, just use the data
        IF ( profiles(i)%samples_in_period == 0 ) THEN
          CALL cube_free(profiles(i)%fields(j)%field_data)
          profiles(i)%fields(j)%field_data =                                  &
                                 extract_var(profiles(i)%fields(j)%var_id)
        ELSE
          ! Otherwise, compare the data for this sample to the previous data and store
          ! appropriately
          var_data = extract_var(profiles(i)%fields(j)%var_id)
          SELECT CASE ( profiles(i)%fields(j)%field_type )            
          CASE ( output_min )
            accum = cube_min(profiles(i)%fields(j)%field_data, var_data)
            CALL cube_safe_copy(profiles(i)%fields(j)%field_data, accum)
          CASE ( output_max )
            accum = cube_max(profiles(i)%fields(j)%field_data, var_data)
            CALL cube_safe_copy(profiles(i)%fields(j)%field_data, accum)
          END SELECT
          ! Free the work cubes
          CALL cube_free(var_data)
          CALL cube_free(accum)
        END IF
                   
      END SELECT
    END DO

    ! Increment the sample count
    profiles(i)%samples_in_period = profiles(i)%samples_in_period+1

    ! Work out when we will next take a sample
    profiles(i)%next_sample_time = datetime_advance(                          &
      profiles(i)%next_sample_time, profiles(i)%sample_period                 &
    )
  END IF

  ! Nothing more to do!!
END DO

RETURN

END SUBROUTINE sample_data
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE output_data()

USE io_constants, ONLY: mdi

USE datetime_mod, ONLY: datetime_advance

USE data_cube_mod, ONLY: data_cube, cube_safe_copy, cube_free,                &
                          operator (+), operator (*), operator (/)

USE model_time_mod, ONLY: current_time, is_spinup, spinup_start, timestep_len

USE file_ts_mod, ONLY: file_ts_write_var, file_ts_advance

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   For each profile, decide if we need to output data at this timestep, and
!   output it if we do
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
TYPE(datetime) :: next_time  ! The time at the next timestep

TYPE(data_cube) :: data_for_output  ! The data cube to output

INTEGER :: i,j  ! Loop counter

!-----------------------------------------------------------------------------


next_time = datetime_advance(current_time, timestep_len)

DO i = 1,nprofiles
  !-----------------------------------------------------------------------------
  ! First check if the profile is currently active
  !-----------------------------------------------------------------------------
  IF ( is_spinup .AND. .NOT. profiles(i)%output_spinup )                      &
  ! If we are in spinup but profile is not outputting, we don't need to do anything
        CYCLE

  IF ( .NOT. is_spinup .AND. ( .NOT. profiles(i)%output_main_run .OR.         &
       ( current_time < profiles(i)%output_start .OR.                         &
         profiles(i)%output_end < current_time ) ) )                          &
  ! If we are in the main run but either not outputting any of the main run or
  ! the time is out of the range of times the profile is outputting, we don't
  ! need to do anything
        CYCLE


  !-----------------------------------------------------------------------------
  ! Open the next output file for the profile, if required
  !
  ! This is the case if:
  !
  !   - The profile is outputting during spinup and we are at in the first
  !     timestep of a spinup cycle
  !
  !   - The profile is outputting part of the main run and we are in the
  !     timestep in the main run when output starts
  !
  ! We do this check in two stages for readability
  !-----------------------------------------------------------------------------
  IF ( is_spinup .AND. profiles(i)%output_spinup .AND.                        &
       current_time == spinup_start )                                         &
    CALL internal_next_output_file(profiles(i))

  IF ( .NOT. is_spinup .AND. profiles(i)%output_main_run .AND.                &
       current_time == profiles(i)%output_start )                             &
    CALL internal_next_output_file(profiles(i))


  !-----------------------------------------------------------------------------
  ! If we are in the last timestep of an output period (i.e. the time at the
  ! start of the next model timestep is equal to (or greater than) the start of
  ! the next output interval), output data
  !-----------------------------------------------------------------------------
  IF ( next_time >= profiles(i)%next_output_time ) THEN

    DO j = 1,profiles(i)%nfields

      SELECT CASE ( profiles(i)%fields(j)%field_type )
      CASE ( output_snapshot, output_min, output_max )
        ! For snapshot, min or max variables, just copy the data array
        CALL cube_safe_copy(data_for_output,                                  &
                            profiles(i)%fields(j)%field_data)
        ! Once data has been output, reset it to mdi
        profiles(i)%fields(j)%field_data%values(:) = mdi

      CASE ( output_accum )
        ! In the case of an accumulation, we output the accumulated value multiplied
        ! by the number of timesteps in a sample period. This adjusts the accumulation
        ! for any intermittent sampling and is designed so that the time accumulation
        ! of a flux (e.g. kg s-1) can easily be converted to a total (e.g. kg) by
        ! subsequently multiplying by the model timestep length during post
        ! processing.
        ! Note that we don't have to worry about using cube_safe_copy since the
        ! multiplication creates a new cube
        data_for_output = profiles(i)%fields(j)%field_data                    &
                * ( REAL(profiles(i)%sample_period) / REAL(timestep_len) )

      CASE ( output_mean )
        ! In the case of a mean, divide the accumulated value by the number of samples
        ! taken
        ! Note that we don't have to worry about using cube_safe_copy since the
        ! division creates a new cube
        data_for_output = profiles(i)%fields(j)%field_data                    &
                        / REAL(profiles(i)%samples_in_period)
        ! Zero the field's data property and start gathering again
        profiles(i)%fields(j)%field_data%values(:) = 0.0
      END SELECT

      ! Write the data to file
      CALL file_ts_write_var(                                                 &
        profiles(i)%fh, profiles(i)%fields(j)%file_id, data_for_output,       &
      ! Subgrid information (for writing a slab of the output grid in parallel mode)
                use_subgrid, subgrid                                          &
              )

      ! Free the temporary cube
      CALL cube_free(data_for_output)
    END DO  ! fields

    !-----------------------------------------------------------------------------
    ! Update tracking variables and, if this is not the last output for the
    ! profile, advance the file
    !-----------------------------------------------------------------------------
    profiles(i)%samples_in_period = 0

    profiles(i)%current_output_time = profiles(i)%next_output_time
    profiles(i)%next_output_time = datetime_advance(                          &
      profiles(i)%next_output_time, profiles(i)%output_period                 &
    )

    ! Check if this is the last output time for the currently open output file
    ! If it is not, then advance the file
    IF ( profiles(i)%current_output_time                                      &
         < profiles(i)%fh%data_end ) THEN
      CALL file_ts_advance(profiles(i)%fh)
    END IF

  END IF  ! Last timestep of output period

  ! Nothing more to do!
END DO

RETURN

END SUBROUTINE output_data
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE close_all()

USE data_cube_mod, ONLY: cube_free

USE file_ts_mod, ONLY: file_ts_close

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Closes all the output files and frees all resources consumed by them
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
INTEGER :: i,j  ! Loop counters


!-----------------------------------------------------------------------------

DO i = 1,nprofiles
  ! Free all the data cubes for the output_fields
  DO j = 1,profiles(i)%nfields
    CALL cube_free(profiles(i)%fields(j)%field_data)
  END DO

  DEALLOCATE(profiles(i)%fields)
  NULLIFY(profiles(i)%fields)

  IF ( profiles(i)%has_open_file ) CALL file_ts_close(profiles(i)%fh)
END DO


END SUBROUTINE close_all

! Internal procedures
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE internal_init_profile_vars(profile)

USE io_constants, ONLY: max_dim_var

USE data_cube_mod, ONLY: cube_create, cube_free

USE datetime_mod, ONLY: datetime_advance

USE model_interface_mod, ONLY: get_var_levs_dims

USE model_time_mod, ONLY: is_spinup, spinup_start, timestep_len

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises fields required to start sampling data for the current section
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(output_profile), INTENT(INOUT) :: profile  ! The profile to initialise
                                                ! fields for


! Work variables
INTEGER :: grid_shape(2)  ! The size of the two grid dimensions (depending
                          ! on whether a subgrid is used)

INTEGER :: var_ndims   ! The number of levels dimensions the variable has
INTEGER :: lev_dim_sizes(max_dim_var)  ! The sizes of the levels dimensions

TYPE(datetime) :: section_start  ! The start time of this 'section' of output

INTEGER :: i  ! Index variable


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Allocate the data cubes for the output fields
!-----------------------------------------------------------------------------
! Work out what shape the grid is
IF ( use_subgrid ) THEN
  grid_shape(:) = (/ subgrid%nx, subgrid%ny /)
ELSE
  grid_shape(:) = (/ grid%nx, grid%ny /)
END IF

DO i = 1,profile%nfields
  ! Get the number and sizes of the levels dimensions for the variable
  CALL get_var_levs_dims(                                                     &
    profile%fields(i)%var_id, ndims = var_ndims, dim_sizes = lev_dim_sizes    &
  )

  ! Deallocate any existing cube first
  CALL cube_free(profile%fields(i)%field_data)

  ! Allocate a cube of the correct shape for the variable
  profile%fields(i)%field_data =                                              &
                   cube_create((/ grid_shape, lev_dim_sizes(1:var_ndims) /))
END DO


!-----------------------------------------------------------------------------
! Initialise the profile so that the sampling and output periods start at
! the correct time
!-----------------------------------------------------------------------------
! Work out the start time for the current section of output
IF ( is_spinup ) THEN
  section_start = spinup_start
ELSE
  section_start = profile%output_start
END IF

! First sample will be at time=sample_period after start of section. We
! reduce by one timestep so as to match the time at the start of the
! timestep.
profile%next_sample_time    = datetime_advance( section_start,                &
                                                profile%sample_period         &
                                                - timestep_len )
profile%current_output_time = section_start
profile%next_output_time    = datetime_advance(section_start, profile%output_period)
profile%samples_in_period   = 0

RETURN

END SUBROUTINE internal_init_profile_vars
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
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION internal_open_output_file(                                           &
  file_start, file_end, output_period, is_time_template, file_name, fields    &
) RESULT(out_file)

USE mpi, ONLY: mpi_comm_world, mpi_info_null

USE io_constants, ONLY: mode_write, max_dim_file, attr_global

USE dictionary_mod, ONLY: dict, dict_create, dict_free

USE data_cube_mod, ONLY: cube_get_data, cube_free

USE file_ts_mod, ONLY: file_ts_open, file_ts_def_grid, file_ts_def_time_dim,  &
                        file_ts_def_attr, file_ts_enddef, file_ts_write_var

USE model_interface_mod, ONLY: get_var_id, extract_var

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Opens a file_ts object for a new output file using the supplied information
!   and returns it
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(datetime), INTENT(IN) :: file_start  ! The start time of the file
TYPE(datetime), INTENT(IN) :: file_end    ! The end time of the file

INTEGER, INTENT(IN) :: output_period  ! The output period of the file

LOGICAL, INTENT(IN) :: is_time_template  ! Indicates if file_name includes
                                         ! time templating strings, or
                                         ! whether it is a single file name
CHARACTER(LEN=max_file_name_len), INTENT(IN) :: file_name
                                         ! The name or time template to use
                                         ! for the file

TYPE(output_field), INTENT(INOUT) :: fields(:)
                                         ! The fields we are defining in the
                                         ! file
                                         ! These are INOUT because the
                                         ! file_id field will be set to the
                                         ! id of the variable in the file

! Return type
TYPE(file_ts) :: out_file


! Work variables
INTEGER :: dummy  ! Dummy variable to receive dimension id of time dim

TYPE(dict) :: defined_dims  ! Dictionary of dimensions defined in the file
                            ! so far

TYPE(data_cube) :: lat_cube, lon_cube  ! Cubes for lat/lon data
REAL :: point_lat, point_lon  ! Real values for lat/lon of single point
                              ! for ASCII files
INTEGER :: lat_id, lon_id  ! Id of lat/lon variable in file for NetCDF files

CHARACTER(LEN=15) :: cell_methods_val  ! Value of cell_methods attribute
                                       ! for the current variable

INTEGER :: i  ! Loop variable


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Open a new file using the properties we were given
!-----------------------------------------------------------------------------
out_file = file_ts_open(mode_write, file_start, file_end,                     &
                                output_period, .FALSE.,                       &
! If is_time_template = T, file_name will be used as a time template
                                  is_time_template, file_name,                &
! If is_time_template = F, file_names and file_times will be used
! So we give a single file name starting at the start of the file (i.e. one
! file for all output)
                                  (/ file_name /), (/ file_start /),          &
! Pass MPI variables to file_ts_open so that parallel I/O will be used
                                  mpi_comm_world, mpi_info_null)


!-----------------------------------------------------------------------------
! Define the grid
!
! Before passing the grid_info object from output_mod.F90, overide the
! dimension names that will be used
!-----------------------------------------------------------------------------
grid%dim_name = grid_dim_name
grid%x_name   = x_dim_name
grid%y_name   = y_dim_name
CALL file_ts_def_grid(out_file, grid)


!-----------------------------------------------------------------------------
! Define the time dimension
!-----------------------------------------------------------------------------
dummy = file_ts_def_time_dim(out_file, time_dim_name)


!-----------------------------------------------------------------------------
! Create a dictionary to store the levels dimensions defined so far
!    dim_name => dim_id
!-----------------------------------------------------------------------------
defined_dims = dict_create(max_dim_file, INT(1))


!-----------------------------------------------------------------------------
! Get the latitude and longitude data cubes
!-----------------------------------------------------------------------------
lat_cube = extract_var(get_var_id('latitude'))
lon_cube = extract_var(get_var_id('longitude'))


!-----------------------------------------------------------------------------
! Do the things that need to be done differently for ASCII and NetCDF
!
! In particular, if we are using ASCII files every variable must have a time
! dimension. However, the grid is restricted to 1 x 1 - so we use global
! attributes for latitude and longitude
!-----------------------------------------------------------------------------
SELECT CASE ( output_format )
CASE ( format_ascii )
  ! For ASCII, populate the lat/lon attributes
  CALL cube_get_data(lat_cube, point_lat)
  CALL cube_get_data(lon_cube, point_lon)

  CALL file_ts_def_attr(out_file, attr_global, 'latitude', point_lat)
  CALL file_ts_def_attr(out_file, attr_global, 'longitude', point_lon)

CASE ( format_ncdf )
  ! For NetCDF, create non-time-varying lat/lon variables
  CALL internal_define_var(                                                   &
    out_file, defined_dims, get_var_id('latitude'), 'latitude', .FALSE., lat_id &
  )

  CALL internal_define_var(                                                   &
    out_file, defined_dims, get_var_id('longitude'), 'longitude', .FALSE., lon_id &
  )

  ! No default case, so that using a format other than those defined is a definite
  ! error
END SELECT


!-----------------------------------------------------------------------------
! Set up the output variables in the file
!-----------------------------------------------------------------------------
DO i = 1,SIZE(fields)

  ! Define the variable in the file (inc. dimensions and attributes)
  CALL internal_define_var(                                                   &
    out_file, defined_dims,                                                   &
    fields(i)%var_id, fields(i)%output_name, .TRUE., fields(i)%file_id        &
  )

  ! Add the CF convention "cell_methods" attribute to indicate whether the field
  ! type of output
  SELECT CASE ( fields(i)%field_type )
  CASE ( output_snapshot )
    cell_methods_val = "time : point"

  CASE ( output_accum )
    cell_methods_val = "time : sum"

  CASE ( output_mean )
    cell_methods_val = "time : mean"
        
  CASE ( output_min )
    cell_methods_val = "time : minimum"
        
  CASE ( output_max )
    cell_methods_val = "time : maximum"

  END SELECT

  CALL file_ts_def_attr(                                                      &
    out_file, fields(i)%file_id, "cell_methods", cell_methods_val             &
  )

END DO

! We have finished defining things on the file handle
CALL file_ts_enddef(out_file)

! We no longer need the defined_dims dictionary
CALL dict_free(defined_dims)


!-----------------------------------------------------------------------------
! Fill the latitude and longitude if we need to
!-----------------------------------------------------------------------------
SELECT CASE ( output_format )
CASE ( format_ascii )
  ! For ASCII files, we do nothing - we use a select statement so we get an
  ! error if someone has added a file type but not considered it's consequences
  ! here

CASE ( format_ncdf )
  CALL file_ts_write_var(                                                     &
    out_file, lat_id, lat_cube,                                               &
  ! Subgrid information (for writing a slab of the output grid in parallel mode)
          use_subgrid, subgrid                                                &
        )
  CALL file_ts_write_var(                                                     &
    out_file, lon_id, lon_cube,                                               &
  ! Subgrid information (for writing a slab of the output grid in parallel mode)
          use_subgrid, subgrid                                                &
        )

  ! No default case, so that using a format other than those defined is a definite
  ! error
END SELECT


!-----------------------------------------------------------------------------
! Deallocate the lat/lon data cubes
!-----------------------------------------------------------------------------
CALL cube_free(lat_cube)
CALL cube_free(lon_cube)


RETURN

END FUNCTION internal_open_output_file
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

SUBROUTINE internal_define_var(FILE, defined_dims, var_id, var_name,          &
                               use_time, var_file_id)

USE io_constants, ONLY: max_dim_var

USE dictionary_mod, ONLY: dict_key_len, dict, dict_char_val_len,              &
                           dict_has_key, dict_get, dict_set, dict_free

USE model_interface_mod, ONLY: get_var_levs_dims, get_var_attrs


USE file_ts_mod, ONLY: file_ts_def_dim, file_ts_def_var, file_ts_def_attr

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   INTERNAL PROCEDURE TO output_mod
!   Defines the variable identified by var_id on the given file_ts object
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ts), INTENT(INOUT) :: FILE  ! The file to define variables on
TYPE(dict), INTENT(INOUT) :: defined_dims
                                      ! Mapping of dim_name => dim_id for
                                      ! dimensions already defined on the
                                      ! file
INTEGER, INTENT(IN) :: var_id
                              ! The variable id of the model variable to
                              ! define on the file
                              ! See model_interface_mod
CHARACTER(LEN=*), INTENT(IN) :: var_name
                              ! The name to use for the variable in
                              ! output files
LOGICAL, INTENT(IN) :: use_time
                              !   T - create the variable as using time
                              !       dimension
                              !   F - create the variable as not using time
                              !       dimension
INTEGER, INTENT(OUT) :: var_file_id
                             ! The id of the variable in the file


! Work variables
INTEGER :: var_ndims   ! The number of levels dimensions the variable has
INTEGER :: lev_dim_sizes(max_dim_var)  ! The sizes of the levels dimensions
CHARACTER(LEN=max_sdf_name_len) :: dim_names(max_dim_var)
                                    ! The names of the levels dimensions for
                                    ! the variable
INTEGER :: dim_ids(max_dim_var)    ! The ids in file of the levels dimensions

TYPE(dict) :: int_attrs, real_attrs, char_attrs  ! Dictionaries containing
                                                 ! attribute values

CHARACTER(LEN=dict_key_len) :: key  ! Used when iterating over attribute
INTEGER :: int_val                  ! dictionaries
REAL :: real_val
CHARACTER(LEN=dict_char_val_len) :: char_val

INTEGER :: i  ! Index variables


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Get the ids of the levels dimensions for the variable, defining them if
! they have not been defined
!-----------------------------------------------------------------------------
CALL get_var_levs_dims(var_id, ndims = var_ndims, dim_names_out = dim_names,  &
                                                dim_sizes = lev_dim_sizes)

DO i = 1,var_ndims
  ! If it has not yet been defined, define the dimension, storing its id
  IF ( .NOT. dict_has_key(defined_dims, dim_names(i)) )                       &
    CALL dict_set(                                                            &
      defined_dims, dim_names(i),                                             &
      file_ts_def_dim(FILE, dim_names(i), lev_dim_sizes(i))                   &
    )

  ! Get the dimension id from the dict and add it to the list for this variable
  CALL dict_get(defined_dims, dim_names(i), dim_ids(i))
END DO

!-----------------------------------------------------------------------------
! Create the variable and store its id
!-----------------------------------------------------------------------------
var_file_id = file_ts_def_var(                                                &
  FILE, var_name, dim_ids(1:var_ndims), use_time                              &
)

!-----------------------------------------------------------------------------
! Define attributes
!-----------------------------------------------------------------------------
CALL get_var_attrs(var_id, int_attrs, real_attrs, char_attrs)

! First the integer valued attributes
DO i = 1,int_attrs%length
  key = int_attrs%keys(i)
  CALL dict_get(int_attrs, key, int_val)

  CALL file_ts_def_attr(FILE, var_file_id, key, int_val)
END DO

! Next, real valued attributes
DO i = 1,real_attrs%length
  key = real_attrs%keys(i)
  CALL dict_get(real_attrs, key, real_val)

  CALL file_ts_def_attr(FILE, var_file_id, key, real_val)
END DO

! Lastly, character valued attributes
DO i = 1,char_attrs%length
  key = char_attrs%keys(i)
  CALL dict_get(char_attrs, key, char_val)

  CALL file_ts_def_attr(FILE, var_file_id, key, char_val)
END DO

! Free the attribute dictionaries
CALL dict_free(int_attrs)
CALL dict_free(real_attrs)
CALL dict_free(char_attrs)

RETURN

END SUBROUTINE internal_define_var

END MODULE output_mod
