! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE time_varying_input_mod

USE io_constants, ONLY: max_sdf_name_len

USE data_cube_mod, ONLY: data_cube

USE file_ts_mod, ONLY: file_ts

USE datetime_mod, ONLY: datetime,                                             &
! Also import the comparison operators on the datetime type
                           operator( == ), operator( /= ), operator( < ),     &
                           operator( > ), operator( <= ), operator( >= )

USE input_mod, ONLY: grid, use_subgrid, subgrid

USE logging_mod, ONLY: log_fatal

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------------
INTEGER, PARAMETER :: nfiles_max = 50  ! The maximum number of file_ts objects
                                       ! that we can read data from

!-----------------------------------------------------------------------------
! Module types
!-----------------------------------------------------------------------------
TYPE input_file

  TYPE(file_ts) :: fh  ! The underlying timeseries file

  INTEGER :: times_lbound, times_ubound  ! The bounds on the times required
                                         ! to fulfill required interpolation
                                         !  -1 - previous data
                                         !   0 - current data
                                         !   1 - next data
                                         !   2 - 2 data steps forward

  TYPE(datetime) :: data_times(-1:2)  ! The times that data apply from
                                      !   -1 - the time that the previous
                                      !        data in file apply from
                                      !    0 - the time that the current
                                      !        data in file apply from
                                      !    1 - the time that the next data
                                      !        in file apply from
                                      !    2 - the time that the data after
                                      !        the next data in file apply
                                      !        from
                                      ! It is guaranteed that data_times(0)
                                      ! and data_times(1) will always be
                                      ! populated, even if not required
                                      ! Other times will only be populated
                                      ! if required

  INTEGER :: tsteps_in_data_period(-1:1)
                                      ! The number of model timesteps in
                                      ! each data period
                                      !   -1 - the number of model timesteps
                                      !        between data_times(-1) and
                                      !        data_times(0)
                                      !    0 - the number of model timesteps
                                      !        between data_times(0) and
                                      !        data_times(1)
                                      !    1 - the number of model timesteps
                                      !        between data_times(1) and
                                      !        data_times(2)

  INTEGER :: current_tstep = 0        ! The number of model timesteps since
                                      ! data_times(0) (i.e. how far into the
                                      ! current interpolation period we are)

  INTEGER :: nfields  ! The number of fields in the file
  TYPE(input_field), POINTER :: fields(:) => NULL()  ! The fields in this file

END TYPE input_file


TYPE input_field

  INTEGER :: var_id  ! The id of the model variable that this field
                     ! populates - see model_interface_mod

  INTEGER :: file_id  ! The id of the variable in the underlying file

  CHARACTER(LEN=2) :: interp_flag  ! The type of interpolation to use

  TYPE(data_cube), POINTER :: DATA(:) => NULL()
                      ! The current data cubes for the field
                      ! One cube per time required for interpolation

END TYPE input_field

!-----------------------------------------------------------------------------
! Module variables
!-----------------------------------------------------------------------------
INTEGER :: nfiles = 0  ! The number of files registered as containing time
                       ! varying input
TYPE(input_file), SAVE :: files(nfiles_max)  ! The registered files

CHARACTER(LEN=max_sdf_name_len) :: time_dim_name = "time"
                       ! The dimension name for the time dimension in input
                       ! files


!-----------------------------------------------------------------------------
! Visibility declarations
!-----------------------------------------------------------------------------
PRIVATE
PUBLIC time_dim_name,                                                         &
       register_input_file, advance_all, seek_all_to_current_datetime,        &
       update_model_variables, close_all


CONTAINS


!-----------------------------------------------------------------------------
! Module procedures
!-----------------------------------------------------------------------------
! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE register_input_file(data_start, data_end, data_period,             &
                               is_climatology, use_time_template, template,   &
                               file_names, file_times, identifiers,           &
                               sdf_names, interp_flags)

USE io_constants, ONLY: max_file_name_len, mode_read, max_sdf_name_len,       &
                         max_dim_file, max_dim_var

USE model_time_mod, ONLY: run_min_time, run_max_time, timestep_len

USE datetime_mod, ONLY: period_month, period_year, datetime_advance

USE dictionary_mod, ONLY: dict, dict_create, dict_get, dict_set,              &
                           dict_has_key, dict_free

USE file_ts_mod, ONLY: file_ts_open, file_ts_def_grid, file_ts_def_dim,       &
                        file_ts_def_time_dim, file_ts_def_var,                &
                        file_ts_enddef

USE model_interface_mod, ONLY: get_var_id, get_var_levs_dims

USE interpolation_mod, ONLY: get_required_time_bounds

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Registers a file or group of files as providing time-varying data for the
!   model variables specified by the given identifiers
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Argument types
!-----------------------------------------------------------------------------
TYPE(datetime), INTENT(IN) :: data_start
                            ! The date and time of the first data
TYPE(datetime), INTENT(IN) :: data_end
                            ! The date and time of the last data
INTEGER, INTENT(IN) :: data_period
                            ! The period of the data
                            ! (in seconds or a 'special' period)
LOGICAL, INTENT(IN) :: is_climatology
                            ! .TRUE. - the data is a climatology
                            ! .FALSE. - the data is not a climatology

LOGICAL, INTENT(IN) :: use_time_template
                            ! .TRUE. - use time templating
                            ! .FALSE. - use lists of file names and times
                            !           of first data in each file

! With time templating
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: template
                            ! The time template to use

! With a file list
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: file_names(:)
                            ! List of file names
TYPE(datetime), INTENT(IN), OPTIONAL :: file_times(:)
                            ! Time of first data for each file

CHARACTER(LEN=*), INTENT(IN) :: identifiers(:)
                            ! The model identifiers of the variables
                            ! this file(s) provides data for
CHARACTER(LEN=*), INTENT(IN) :: sdf_names(:)
                            ! The name of each variable in the file(s)
CHARACTER(LEN=*), INTENT(IN) :: interp_flags(:)
                            ! The type of interpolation to use for each
                            ! variable


!-----------------------------------------------------------------------------
! Work variables
!-----------------------------------------------------------------------------
! Local variables passed to file_ts_open
! In the case where optional variables are not given, these are set to values
! that will cause errors if they are required, but will be ignored if they
! are not required
CHARACTER(LEN=max_file_name_len) :: template_local
CHARACTER(LEN=max_file_name_len), ALLOCATABLE :: file_names_local(:)
TYPE(datetime), ALLOCATABLE :: file_times_local(:)

INTEGER :: nvars  ! The number of variables in the file

TYPE(dict) :: file_dim_ids  ! Dictionary containing the dimension ids defined
                            ! so far
                            ! Maps dim_name => dim_id

! Used to define dimensions and variables in file
INTEGER :: ndims  ! The number of levels dimensions the current variable has
CHARACTER(LEN=max_sdf_name_len) :: dim_names(max_dim_var)
                  ! The levels dimension names for the current variable
INTEGER :: dim_sizes(max_dim_var)
                  ! The sizes of the levels dimensions for the current variable
INTEGER :: dim_ids(max_dim_var)
                  ! The dimension ids in file of the levels dimensions for the
                  ! current variable

INTEGER :: dummy  ! This is used to store the result from defining the time
                  ! dimension - we don't need to keep it
CHARACTER(LEN=len(sdf_names)) :: sdf_name_local
                  ! A local copy of a value from sdf_names, for alteration.

INTEGER :: i, j  ! Index variables


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check arguments
!-----------------------------------------------------------------------------
! Check that we have space to register another file
IF ( nfiles >= nfiles_max )                                                   &
  CALL log_fatal("register_input_file",                                       &
                 "Too many files registered - try increasing NFILES_MAX")

! Check that we have an sdf_name and an interpolation flag for each variable
nvars = SIZE(identifiers)
IF ( nvars /= SIZE(sdf_names) )                                               &
  CALL log_fatal("register_input_file",                                       &
                 "identifiers and sdf_names must have the same number " //    &
                 "of elements")
IF ( nvars /= SIZE(interp_flags) )                                            &
  CALL log_fatal("register_input_file",                                       &
                 "identifiers and interp_flags must have the same " //        &
                 "number of elements")


!-----------------------------------------------------------------------------
! Check that data will be provided for the whole run at an appropriate
! timestep
!-----------------------------------------------------------------------------
! If not using a climatology, data must be provided for the whole run
! If using a climatology, the data must cover a whole year, but this is checked
! by file_ts_open
IF ( .NOT. is_climatology .AND.                                               &
     ( data_start > run_min_time .OR. data_end < run_max_time ) )             &
  CALL log_fatal("register_input_file",                                       &
                 "Each input file must provide data for the entire run")

! Unless data is given on a 'special' period, the data period must be a whole
! number of model timesteps
IF ( data_period /= period_month .AND. data_period /= period_year .AND.       &
     MOD(data_period, timestep_len) /= 0 )                                    &
  CALL log_fatal("register_input_file",                                       &
                 "Data period must be a special period or a multiple " //     &
                 "of timestep length")


!-----------------------------------------------------------------------------
! Set up the optional arguments depending on what is available
!-----------------------------------------------------------------------------
IF ( PRESENT(template) ) THEN
  template_local = template
ELSE
  ! If template is not provided, we provide an empty string
  ! This will cause an error unless a file list is specified
  template_local = ""
END IF

IF ( PRESENT(file_names) ) THEN
  ! If file_names is present, copy its values into the local counterpart
  ! that has fixed length strings
  ALLOCATE(file_names_local(SIZE(file_names)))
  file_names_local(:) = file_names(:)
ELSE
  ! If file_names is not provided, we provide one file that doesn't exist
  ! This will cause an error unless time templating is specified
  ALLOCATE(file_names_local(1))
  file_names_local(1) = ""
END IF

IF ( PRESENT(file_times) ) THEN
  ! If file_times is present, copy its values into the local counterpart
  ALLOCATE(file_times_local(SIZE(file_times)))
  file_times_local(:) = file_times(:)
ELSE
  ! If file_times is not provided, we provide one time that is not equal to
  ! data_start
  ! This will cause an error unless time templating is specified
  ALLOCATE(file_times_local(1))
  file_times_local(1) = datetime_advance(data_start, 1)
END IF


!-----------------------------------------------------------------------------
! Open the file handle
!-----------------------------------------------------------------------------
nfiles = nfiles + 1

files(nfiles)%fh = file_ts_open(mode_read, data_start, data_end,              &
                                data_period, is_climatology,                  &
                                use_time_template, template_local,            &
                                file_names_local, file_times_local)

! We have finished with the local pointers
DEALLOCATE(file_names_local)
DEALLOCATE(file_times_local)


!-----------------------------------------------------------------------------
! Define the grid
!-----------------------------------------------------------------------------
CALL file_ts_def_grid(files(nfiles)%fh, grid)


!-----------------------------------------------------------------------------
! Define the time dimension
!-----------------------------------------------------------------------------
dummy = file_ts_def_time_dim(files(nfiles)%fh, time_dim_name)


!-----------------------------------------------------------------------------
! Allocate space for the input fields
!-----------------------------------------------------------------------------
files(nfiles)%nfields = nvars
ALLOCATE(files(nfiles)%fields(nvars))


!-----------------------------------------------------------------------------
! Get the upper and lower bounds to use for the time dimension based on what
! interpolation each variable is using
!-----------------------------------------------------------------------------
CALL get_required_time_bounds(                                                &
  interp_flags, files(nfiles)%times_lbound, files(nfiles)%times_ubound        &
)


!-----------------------------------------------------------------------------
! Define the required variables
!-----------------------------------------------------------------------------
! Define a dictionary to gather the dimension ids
file_dim_ids = dict_create(max_dim_file, INT(1))

DO i = 1,nvars

  ! Get the integer id from model_interface_mod for the identifier
  files(nfiles)%fields(i)%var_id = get_var_id(identifiers(i))

  ! Get the levels dims used by this variable - we only care about the names
  ! used in input files
  CALL get_var_levs_dims(files(nfiles)%fields(i)%var_id, ndims = ndims,       &
                         dim_names_in = dim_names, dim_sizes = dim_sizes)

  DO j = 1,ndims
    ! If it has not yet been defined, define the dimension, storing its id
    IF ( .NOT. dict_has_key(file_dim_ids, dim_names(j)) )                     &
      CALL dict_set(                                                          &
        file_dim_ids, dim_names(j),                                           &
        file_ts_def_dim(files(nfiles)%fh, dim_names(j), dim_sizes(j))         &
      )

    ! Get the dimension id from the dict and add it to the list for this variable
    CALL dict_get(file_dim_ids, dim_names(j), dim_ids(j))
  END DO

  ! If sdf_names is empty, use the identifier.
  sdf_name_local = sdf_names(i)
  IF ( LEN_TRIM(sdf_name_local) == 0 ) THEN
    !     Check variable is long enough.
    !     Note that a fatal error will not occur as long as the code declares
    !     sdf_names with lengths >= those of identifiers.
    IF ( LEN_TRIM(identifiers(i)) > LEN(sdf_name_local) )                     &
      CALL log_fatal("register_input_file",                                   &
                     "identifier too long for sdf_name. " //                  &
                     "identifier: " // TRIM(identifiers(i)) )
    sdf_name_local = identifiers(i)
  END IF

  ! Create the variable and store its id
  files(nfiles)%fields(i)%file_id = file_ts_def_var(                          &
    files(nfiles)%fh, sdf_name_local, dim_ids(1:ndims), .TRUE.                &
  )

  ! Allocate space for the data cubes for the field
  ALLOCATE(files(nfiles)%fields(i)%DATA(                                      &
    files(nfiles)%times_lbound:files(nfiles)%times_ubound                     &
  ))

  ! Store the interpolation flag for the field
  files(nfiles)%fields(i)%interp_flag = interp_flags(i)

END DO


!-----------------------------------------------------------------------------
! We are done - clean up
!-----------------------------------------------------------------------------
! Free the dimensions dictionary
CALL dict_free(file_dim_ids)

! Take the file out of define mode
CALL file_ts_enddef(files(nfiles)%fh)

RETURN

END SUBROUTINE register_input_file
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE advance_all()

USE precision_mod, ONLY: int64

USE datetime_mod, ONLY: period_month, period_year,                            &
                         datetime_advance, datetime_diff, datetime_clone

USE data_cube_mod, ONLY: cube_safe_copy, cube_free

USE model_time_mod, ONLY: is_spinup, spinup_start, spinup_end, timestep_len

USE file_ts_mod, ONLY: file_ts_advance, file_ts_seek_to_datetime,             &
                        file_ts_read_var

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Advances all the input files so that they are able to provide data for
!   the current model time when the timestep has been advanced normally (i.e.
!   by one model timestep)
!   If the current model time has been set abnormally (i.e. resetting to the
!   start of spinup period), use seek_to_current_datetime instead
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
TYPE(datetime) :: corrected_dt  ! Next data time corrected for end of spinup
                                ! if required
INTEGER :: secs_into_period  ! The number of seconds into a data_period
                             ! that spinup_start is

INTEGER(KIND=int64) :: diff_secs
                      ! Used in calculation of tsteps_in_data_period
                      ! The number of seconds between two times

! Work variables
INTEGER :: i, j, k  ! Loop counters

!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Overview of algorithm (for each file) in psuedocode:
!
!   IF ( read is not required ) exit
!
!   advance underlying file
!
!   FOREACH variable IN file
!     shift data back one timestep
!
!     read next slab of data into last time slot (vacated by the shift above)
!   END FOREACH
!
!   adjust data times
!
!-----------------------------------------------------------------------------


DO i = 1,nfiles
  ! Every time we advance, we need to increase the counter for number of model
  ! timesteps into a data period we are
  files(i)%current_tstep = files(i)%current_tstep+1

  ! Check if the file is required to read more data (i.e. we have read all the
  ! timesteps in this data period)
  ! If not, then there is nothing more to do
  IF ( files(i)%current_tstep < files(i)%tsteps_in_data_period(0) ) CYCLE


  !-----------------------------------------------------------------------------
  ! If we get to here, then we need to read the next data for this file
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! First, update the metadata for data times etc.
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! Update data_times as if we are proceeding normally
  !-----------------------------------------------------------------------------
  ! We know that on entering this routine, data_times(-1:2) were all correctly
  ! populated (since a seek will have been done before)
  ! This must be the case for the new timestep when this routine exits

  ! The shift means that all times in the range -1:1 are correctly populated
  ! for the new timestep
  files(i)%data_times = EOSHIFT(files(i)%data_times, shift = 1,               &
  ! Since data_times is an array of derived type objects, Intel insists on
  ! boundary being provided (makes sense since it doesn't know what to put
  ! there!). We don't actually care what it is since it will be overwritten later
                                boundary = datetime_clone(files(i)%data_times(1)) &
      )
  ! Populate data_times(2) from data_times(1)
  files(i)%data_times(2) = datetime_advance(                                  &
    files(i)%data_times(1), files(i)%fh%data_period                           &
  )

  !-----------------------------------------------------------------------------
  ! Update tsteps_in_data_period in a similar fashion using the calculated
  ! data times
  !-----------------------------------------------------------------------------
  ! We know that on entering this routine, tsteps_in_data_period(-1:1) were
  ! all correctly populated, so after the shift, -1:0 will be correct for the
  ! new timestep
  files(i)%tsteps_in_data_period = EOSHIFT(files(i)%tsteps_in_data_period, shift = 1)

  ! All that remains (in a normal situation) is to calculate the value for
  ! index 1 from the new data_times
  diff_secs = datetime_diff(files(i)%data_times(1), files(i)%data_times(2))

  IF ( MOD(diff_secs, INT(timestep_len, int64)) /= 0 )                        &
    CALL log_fatal("advance_all",                                             &
                   "Data should be a whole number of model timesteps apart")

  files(i)%tsteps_in_data_period(1) = REAL(diff_secs) / REAL(timestep_len)


  !-----------------------------------------------------------------------------
  ! If we are coming to the end of a spinup cycle, we want to correct the
  ! data_times (and tsteps_in_data_period) so that the data ramp nicely across
  ! the break in time, rather than changing instantaneously
  !
  ! Note that this corrects the data times for times in the future, even
  ! when they are not required. Hence we check below to see if the required
  ! data times are out of sync with the underlying file - if they are, we need
  ! to seek the file
  !
  ! This assumes that we will be starting a new cycle of spinup
  ! In the case of the start of the main run, the files are seeked to the
  ! start of the main run with no ramping in next_time
  ! This means that running a spinup and starting from the end dump of that
  ! is equivalent to running a spinup and then carrying on with the main run
  !-----------------------------------------------------------------------------
  IF ( is_spinup .AND. files(i)%data_times(2) >= spinup_end ) THEN
    ! Adjust data_times(2) so that the time wraps around to the first data time
    ! after spinup start
    ! To do this, we have to find the first data time after spinup start
    SELECT CASE ( files(i)%fh%data_period )
    CASE ( 1: )
      ! In the case of a period of seconds >= 1, we start by assuming that the
      ! corrected datetime will be exactly spinup start
      corrected_dt = spinup_start

      ! Work out how many seconds into a data period that datetime is

      ! The IBM XL Fortran compiler requires explicit conversion between
      ! the default integer kind and int64:
      !    datetime_diff returns int64
      !    files(i)%fh%data_period => int64
      !    MOD in 64-bit
      !    result => default real kind for secs_into_period
      ! Some compilers apply implicit conversion:
      ! For gfortran, we can confirm this using -Wconversion
      ! For intel, results confirm that this is what is happening
      secs_into_period = MOD(                                                 &
        datetime_diff(files(i)%fh%data_start, corrected_dt),                  &
        INT(files(i)%fh%data_period, int64)                                   &
      )

      ! If we guessed the wrong time, make a correction by adding on enough
      ! seconds to take us to the next data time
      IF ( secs_into_period > 0 )                                             &
        corrected_dt = datetime_advance(                                      &
          corrected_dt, files(i)%fh%data_period - secs_into_period            &
        )

    CASE ( period_month )
      ! In the case of a monthly period, we know that we can get the start of the
      ! data period containing spinup_start just by setting the day and time to the
      ! beginning of the month
      corrected_dt = datetime_clone(spinup_start)
      corrected_dt%day  = 1
      corrected_dt%time = 0

      ! Only if that time is strictly less than spinup_start do we need to correct by
      ! adding on a data period
      IF ( corrected_dt < spinup_start )                                      &
        corrected_dt = datetime_advance(corrected_dt, period_month)

    CASE ( period_year )
      ! In the case of a yearly period, we know that we can get the start of the
      ! data period containing spinup_start just by setting the month, day and time
      ! to the beginning of the year
      corrected_dt = datetime_clone(spinup_start)
      corrected_dt%month = 1
      corrected_dt%day   = 1
      corrected_dt%time  = 0

      ! Only if that time is strictly less than spinup_start do we need to correct by
      ! adding on a data period
      IF ( corrected_dt < spinup_start )                                      &
        corrected_dt = datetime_advance(corrected_dt, period_year)
    END SELECT

    files(i)%data_times(2) = corrected_dt

    ! Adjust tsteps_in_data_period so that the number of timesteps in the affected
    ! period is the number until the end of spinup + the number from start of
    ! spinup until the next data
    diff_secs = datetime_diff(files(i)%data_times(1), spinup_end)             &
              + datetime_diff(spinup_start, files(i)%data_times(2))

    IF ( MOD(diff_secs, INT(timestep_len, int64)) /= 0 )                      &
      CALL log_fatal("advance_all",                                           &
                     "Each spinup cycle must contain a whole number of " //   &
                     "timesteps")

    files(i)%tsteps_in_data_period(1) = REAL(diff_secs) / REAL(timestep_len)
  END IF

  ! We are now in the first timestep of a new data period
  files(i)%current_tstep = 0


  !-----------------------------------------------------------------------------
  ! We need to make sure that the file is at times_ubound before we update the
  ! data
  !
  ! In the majority of cases this will be a regular advance, but in the case
  ! of rolling over a spinup cycle, it will require a seek
  !-----------------------------------------------------------------------------
  ! To determine if we need a seek, we compare data_times(times_ubound) to
  ! data_times(times_ubound - 1)
  IF ( files(i)%data_times(files(i)%times_ubound) >                           &
       files(i)%data_times(files(i)%times_ubound-1) ) THEN
    ! If we are moving forward in time normally, then assume that we just require
    ! an advance
    CALL file_ts_advance(files(i)%fh)
  ELSE
    ! Otherwise seek
    CALL file_ts_seek_to_datetime(                                            &
      files(i)%fh, files(i)%data_times(files(i)%times_ubound)                 &
    )
  END IF


  !-----------------------------------------------------------------------------
  ! Update the actual data
  !-----------------------------------------------------------------------------
  DO j = 1,files(i)%nfields
    ! Shift data back one timestep
    DO k = files(i)%times_lbound,(files(i)%times_ubound-1)
      ! Copy the state from the next timestep into the current timestep
      CALL cube_safe_copy(                                                    &
        files(i)%fields(j)%DATA(k), files(i)%fields(j)%DATA(k+1)              &
      )
    END DO

    ! Read next slab of data into the last timeslot
    CALL cube_free(files(i)%fields(j)%DATA(files(i)%times_ubound))
    files(i)%fields(j)%DATA(files(i)%times_ubound) = file_ts_read_var(        &
      files(i)%fh, files(i)%fields(j)%file_id,                                &
    ! Subgrid information from input_mod
            use_subgrid, subgrid                                              &
          )
  END DO  ! fields

END DO  ! files

RETURN

END SUBROUTINE advance_all
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE seek_all_to_current_datetime()

USE precision_mod, ONLY: int64

USE datetime_mod, ONLY: period_month, period_year,                            &
                         datetime_clone, datetime_diff,                       &
                         datetime_advance, datetime_subtract

USE model_time_mod, ONLY: current_time, timestep_len

USE file_ts_mod, ONLY: file_ts_seek_to_datetime, file_ts_read_var,            &
                        file_ts_advance

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Sets up the input files so that they are able to provide data for the
!   current timestep next time they are asked for data
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
TYPE(datetime) :: data_dt  ! The datetime for data that applies for
                           ! current_dt
                           ! See step 1 in algorithm description below
TYPE(datetime) :: seek_dt  ! The datetime to seek the underlying file to
                           ! See step 2 of algorithm below

INTEGER(KIND=int64) :: diff_secs
    ! Used as return value from calls to datetime_diff (hence why an int64)

INTEGER :: time_into_period
    ! For the 1: case, this is the number of seconds into a data period

INTEGER :: i,j,k  ! Loop counters

!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Overview of algorithm (for each file):
!
!   1. Find the data time for the current datetime (i.e. the closest
!      time going backwards from current_dt that data exists for)
!   2. Depending what times are needed for interpolation, seek the underlying
!      file to the appropriate time
!   3. Read the timesteps required for interpolation
!   4. Calculate the data times for all times from -1 to 2, even if they are
!      not required (and even if they are not actually part of the data)
!   5. Calculate the number of model timesteps in each data period and the
!      number of timesteps into the current interpolation period we are
!
!-----------------------------------------------------------------------------


DO i = 1,nfiles
  !-----------------------------------------------------------------------------
  ! Step 1 - find the data time for the requested datetime for this file
  !-----------------------------------------------------------------------------
  SELECT CASE ( files(i)%fh%data_period )
  CASE ( 1: )
    ! In the case of a period of seconds >= 1, we need to know how many seconds
    ! past the closest data time current_dt is, and take that off current_dt

    ! The IBM XL Fortran compiler requires explicit conversion between
    ! the default integer kind and int64:
    !    diff_secs is int64
    !    files(i)%fh%data_period => int64
    !    MOD in 64-bit
    !    result => default real kind for time_into_period
    ! Some compilers apply implicit conversion:
    ! For gfortran, we can confirm this using -Wconversion
    ! For intel, results confirm that this is what is happening
    diff_secs = datetime_diff(files(i)%fh%data_start, current_time)
    time_into_period = MOD(diff_secs, INT(files(i)%fh%data_period, int64))
    data_dt = datetime_subtract(current_time, time_into_period)

  CASE ( period_month )
    ! In the case of a monthly period, we know that the data applies from the
    ! start of the month
    data_dt = datetime_clone(current_time)
    data_dt%day  = 1
    data_dt%time = 0

  CASE ( period_year )
    ! In the case of a yearly period, we know that the data will apply from the
    ! start of the year
    data_dt = datetime_clone(current_time)
    data_dt%month = 1
    data_dt%day   = 1
    data_dt%time  = 0
  END SELECT

  ! Update the time of the current data
  files(i)%data_times(0) = data_dt


  !-----------------------------------------------------------------------------
  ! Step 2 - depending on the interpolation in use for this file, seek the
  !          underlying timeseries file to the correct place
  !-----------------------------------------------------------------------------
  seek_dt = datetime_clone(data_dt)
  ! Because data_period could be a special period, which we can't meaningfully
  ! multiply, we do this by looping
  DO j = 1,ABS(files(i)%times_lbound)
    IF ( files(i)%times_lbound < 0 ) THEN
      seek_dt = datetime_subtract(seek_dt, files(i)%fh%data_period)
    ELSE
      seek_dt = datetime_advance(seek_dt, files(i)%fh%data_period)
    END IF
  END DO

  CALL file_ts_seek_to_datetime(files(i)%fh, seek_dt)


  !-----------------------------------------------------------------------------
  ! Step 3 - read the number of timesteps required for interpolation for
  !          each field
  !-----------------------------------------------------------------------------
  DO j = files(i)%times_lbound,files(i)%times_ubound
    DO k = 1,files(i)%nfields
      ! Read the slab of data for this time
      files(i)%fields(k)%DATA(j) = file_ts_read_var(                          &
        files(i)%fh, files(i)%fields(k)%file_id,                              &
      ! Subgrid information from input_mod
                use_subgrid, subgrid                                          &
              )
    END DO  ! fields

    ! If we need to read more data now, advance the file
    IF ( j /= files(i)%times_ubound ) CALL file_ts_advance(files(i)%fh)
  END DO  ! times


  !-----------------------------------------------------------------------------
  ! Step 4 - Populate data_times for all times -1 to 2
  !-----------------------------------------------------------------------------
  ! We know data_times(0), so populate the other times from that
  files(i)%data_times(-1) = datetime_subtract(                                &
    files(i)%data_times(0), files(i)%fh%data_period                           &
  )

  DO j = 1,2
    files(i)%data_times(j) = datetime_advance(                                &
      files(i)%data_times(j-1), files(i)%fh%data_period                       &
    )
  END DO

  !-----------------------------------------------------------------------------
  ! Step 5 - Populate tsteps_in_data_period and current_tstep
  !-----------------------------------------------------------------------------
  ! Populate tsteps_in_data_period from the data times
  DO j = -1,1
    diff_secs = datetime_diff(                                                &
      files(i)%data_times(j), files(i)%data_times(j+1)                        &
    )

    IF ( MOD(diff_secs, INT(timestep_len, int64)) /= 0 )                      &
      CALL log_fatal("seek_all_to_current_datetime",                          &
                     "Data should be a whole number of model timesteps apart")

    files(i)%tsteps_in_data_period(j) = REAL(diff_secs) / REAL(timestep_len)
  END DO

  ! Work out how many model timesteps into the interpolation period we are
  diff_secs = datetime_diff(files(i)%data_times(0), current_time)

  IF ( MOD(diff_secs, INT(timestep_len, int64)) /= 0 )                        &
    CALL log_fatal("seek_all_to_current_datetime",                            &
                   "Data should be a whole number of model timesteps apart")

  files(i)%current_tstep = REAL(diff_secs) / REAL(timestep_len)

END DO  ! files

RETURN

END SUBROUTINE seek_all_to_current_datetime
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE update_model_variables()

USE data_cube_mod, ONLY: data_cube, cube_free

USE interpolation_mod, ONLY: interpolate

USE model_interface_mod, ONLY: populate_var

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Updates model variables with data for the current model time
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
TYPE(data_cube) :: interp_data  ! Cube for the interpolated data that can be freed

INTEGER :: i,j  ! Loop counter


!-----------------------------------------------------------------------------

DO i = 1,nfiles
  DO j = 1,files(i)%nfields
    ! We want to fill the model variable with interpolated data
    interp_data = interpolate(                                                &
      files(i)%fields(j)%DATA,                                                &
      files(i)%fields(j)%interp_flag,                                         &
      files(i)%tsteps_in_data_period,                                         &
      files(i)%current_tstep                                                  &
    )

    CALL populate_var(files(i)%fields(j)%var_id, interp_data)

    CALL cube_free(interp_data)
  END DO  ! fields
END DO  ! files

RETURN

END SUBROUTINE update_model_variables
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
!   Closes all the time varying input files and frees all resources
!   consumed by them
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
INTEGER :: i, j, k  ! Loop counters


!-----------------------------------------------------------------------------

DO i = 1,nfiles
  DO j = 1,files(i)%nfields
    ! Deallocate all the data cubes for the input field
    DO k = files(i)%times_lbound,files(i)%times_ubound
      CALL cube_free(files(i)%fields(j)%DATA(k))
    END DO

    ! Deallocate the cube store for the field
    DEALLOCATE(files(i)%fields(j)%DATA)
    NULLIFY(files(i)%fields(j)%DATA)
  END DO

  ! Deallocate the fields array for the file
  DEALLOCATE(files(i)%fields)
  NULLIFY(files(i)%fields)

  ! Close the underlying file
  CALL file_ts_close(files(i)%fh)
END DO

RETURN

END SUBROUTINE close_all

END MODULE time_varying_input_mod
