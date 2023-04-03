! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

MODULE file_ts_mod

USE io_constants, ONLY: max_file_name_len, max_sdf_name_len,                  &
                         max_attr_val_len, max_dim_file, max_var_file,        &
                         max_attr_file

USE datetime_mod, ONLY: datetime,                                             &
! Also import the comparison operators on the datetime type
                           operator( == ), operator( /= ), operator( < ),     &
                           operator( > ), operator( <= ), operator( >= )

USE grid_utils_mod, ONLY: grid_info

USE string_utils_mod, ONLY: to_string

USE file_gridded_mod, ONLY: file_gridded

USE dictionary_mod, ONLY: dict

USE logging_mod, ONLY: log_info, log_warn, log_fatal

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Constants
!-----------------------------------------------------------------------------
CHARACTER(LEN=max_sdf_name_len) :: time_index_var_name = 'time'

!-----------------------------------------------------------------------------
! Type definitions
!-----------------------------------------------------------------------------
! Internal type containing information about a dimension
TYPE dim_ts

  INTEGER :: current_id  ! The id of the dimension in the currently open file
  CHARACTER(LEN=max_sdf_name_len) :: NAME ! The name of the dimension
  INTEGER :: length  ! The length of the dimension

END TYPE dim_ts


! Internal type containing information about a variable
TYPE var_ts

  INTEGER :: current_id  ! The id of the variable in the currently open file

  CHARACTER(LEN=max_sdf_name_len) :: NAME  ! The name of the variable in the
                                           ! file(s)

  INTEGER, POINTER :: levels_dims(:)  ! The ids of the vertical levels dimensions
                                      ! (as indices in the dims array on the
                                      ! file_ts object)

  LOGICAL :: use_time   ! Indicates whether the variable uses the time dimension

  ! Dictionaries to hold attribute values
  TYPE(dict) :: attrs_real  ! The values for real valued attributes
                            ! Indexed by attribute name
  TYPE(dict) :: attrs_int  ! The values for integer valued attributes
                           ! Indexed by attribute name
  TYPE(dict) :: attrs_char  ! The values for char valued attributes
                            ! Indexed by attribute name

END TYPE var_ts


TYPE file_ts
  !-----------------------------------------------------------------------------
  ! This type and associated functions and subroutines are intended to provide
  ! access to timeseries data in a unified way, whether that data is contained
  ! in a single file or multiple files
  !-----------------------------------------------------------------------------
  INTEGER :: mode ! The mode in which to open files
                  ! This can be one of mode_read or mode_write
  LOGICAL :: use_mpiio = .FALSE.  ! T - Pass values for MPI_Comm and MPI_Info
                                  !     down the chain to open files for
                                  !     parallel access
                                  ! F - Open files for serial access
  INTEGER :: comm = -1  ! MPI communicator to use
  INTEGER :: info = -1  ! MPI info type to use


  !-----------------------------------------------------------------------------
  ! Properties defining the characteristics of the data
  !-----------------------------------------------------------------------------
  TYPE(datetime) :: data_start ! The date and time of the first data
  TYPE(datetime) :: data_end   ! The date and time of the last data
  INTEGER :: data_period ! The period of the data in the files
                         ! (in seconds or a 'special' period)
                         ! NOTE: that this means all data in a file must
                         ! have the same period
  LOGICAL :: is_climatology = .FALSE. ! Indicates whether the data represents
                                      ! a climatology.
                                      ! In this case, the same data are
                                      ! reused for every year


  !-----------------------------------------------------------------------------
  ! Properties containing information about the dimensions, variables and
  ! attributes contained in each file - this will need to be redefined for
  ! each file that is opened
  !-----------------------------------------------------------------------------
  LOGICAL :: define_mode = .TRUE.  ! Indicates whether object is in
                                   ! 'define mode', i.e. if we can define
                                   ! dimensions and variables

  TYPE(grid_info) :: grid  ! The grid that variables in the file(s) are on

  LOGICAL :: has_time_dim = .FALSE.  ! Indicates whether the time dimension
                                     ! has been defined
  TYPE(dim_ts) :: time_dim ! The time dimension for the file(s)

  INTEGER :: ndims = 0  ! The number of non-grid, non-time dimensions
                        ! that have been defined
  TYPE(dim_ts) :: dims(max_dim_file) ! The non-grid, non-time dimensions
                                     ! in the file(s)

  INTEGER :: nvars = 0  ! The number of variables that have been defined
  TYPE(var_ts) :: vars(max_var_file) ! The variables in the file(s)

  ! Dictionaries to hold global attribute values
  TYPE(dict) :: attrs_real  ! The values for real valued attributes
                                ! Indexed by attribute name
  TYPE(dict) :: attrs_int  ! The values for integer valued attributes
                               ! Indexed by attribute name
  TYPE(dict) :: attrs_char  ! The values for char valued attributes
                                ! Indexed by attribute name

  INTEGER :: time_index_var_id  ! The variable id for the time index
                                ! in the currently open file
  INTEGER :: time_bounds_var_id ! The variable id for the time bounds
                                ! in the currently open file
                                ! Note that since the time index and time
                                ! bounds do not have a grid, we have to
                                ! use the file_handle underlying
                                ! the file_gridded object to manipulate
                                ! it, hence these are variable ids in the
                                ! file_handle object, i.e. file_ts%open_file%fh


  !-----------------------------------------------------------------------------
  ! Properties used to determine what file to open
  !-----------------------------------------------------------------------------
  LOGICAL :: use_time_template

  ! With time templating
  CHARACTER(LEN=max_file_name_len) :: template
  INTEGER :: tpl_period

  ! With a file list
  INTEGER :: nfiles ! The number of files
  CHARACTER(LEN=max_file_name_len), POINTER :: files(:) => NULL()
                     ! List of file names
  TYPE(datetime), POINTER :: file_times(:) => NULL()
                     ! Time of first data for each file


  !-----------------------------------------------------------------------------
  ! Properties related to the currently open file
  !-----------------------------------------------------------------------------
  LOGICAL :: has_open_file = .FALSE.
                            ! T - the file in open_file represents an
                            !     actual open file
                            ! F - the file object in open_file is uninitialised

  INTEGER :: open_file_index ! The index in the files/file_times arrays
                             ! of the currently open file

  TYPE(file_gridded) :: open_file ! The currently open file

  TYPE(datetime) :: next_file_start ! The start time of the next file

  TYPE(datetime) :: current_datetime ! The datetime that this file_ts
                                     ! object is currently pointing to
                                     ! NOTE: it is possible for this datetime
                                     ! to contain a different value to
                                     ! the current model datetime - it is
                                     ! maintained to ensure the correct file
                                     ! is kept open

  LOGICAL :: timestamp_written = .FALSE.  ! Used in write mode only
                                          ! Indicates if the timestamp has
                                          ! been written for the current
                                          ! timestep
                                          ! It is written on first write

END TYPE file_ts


! Overloads for file_ts_def_attr
INTERFACE file_ts_def_attr
MODULE PROCEDURE file_ts_def_attr_real, file_ts_def_attr_int,                 &
                 file_ts_def_attr_char
END INTERFACE file_ts_def_attr


!-----------------------------------------------------------------------------
! Visibility declarations
!-----------------------------------------------------------------------------
PRIVATE
PUBLIC                                                                        &
! Types
         file_ts,                                                             &
! Routines for opening and closing
         file_ts_open, file_ts_close,                                         &
! Routines for definitions
         file_ts_def_grid, file_ts_def_dim,                                   &
         file_ts_def_time_dim, file_ts_def_var,                               &
         file_ts_def_attr, file_ts_enddef,                                    &
! Routines for seeking
         file_ts_seek_to_datetime, file_ts_advance,                           &
! Routines for reading and writing
         file_ts_read_var, file_ts_write_var


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


FUNCTION file_ts_open(mode, data_start, data_end, data_period, is_climatology,&
                      use_time_template, template, files, file_times,         &
                      comm, info)                                             &
                              RESULT(FILE)

USE io_constants, ONLY: mode_read, mode_write

USE datetime_mod, ONLY: period_month, period_year, datetime_advance

USE templating_mod, ONLY: tpl_detect_period

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Opens a timeseries file and returns a file_ts object representing it
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
INTEGER, INTENT(IN) :: mode ! The mode to open the file
                            ! One of mode_read or mode_write

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
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: template
                            ! The time template to use

! With a file list
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: files(:)
                            ! List of file names
TYPE(datetime), OPTIONAL, INTENT(IN) :: file_times(:)
                            ! Time of first data for each file

! MPI variables
INTEGER, OPTIONAL, INTENT(IN) :: comm  ! MPI communicator to use for parallel
                                       ! I/O
                                       ! If not given, serial I/O is used
INTEGER, OPTIONAL, INTENT(IN) :: info  ! MPI info object to use for parallel
                                       ! I/O
                                       ! If not given, serial I/O is used

! Return type
TYPE(file_ts) :: FILE


! Work variables
INTEGER :: nfiles ! The number of files in the files/file_times lists
INTEGER :: tpl_period ! The detected templating period

TYPE(datetime) :: prev_file_start  ! The start time of the previous file
                                   ! in the list
                                   ! Used when checking that files
                                   ! given in a list are in chronological
                                   ! order
LOGICAL :: file_exists  ! Indicates if the file exists
                        ! Used when checking that files given in a list
                        ! exist

INTEGER :: i ! Loop counter

CHARACTER(LEN=max_file_name_len) :: fname   !  Name of file or template, for
                                            !  error message.

!-----------------------------------------------------------------------------


! Initialise nfiles - this is not used when time-templating is enabled
nfiles = -1

! Load name of file or template, for error messages.
fname = ''  !  initialise, in case optional argument not provided
IF ( use_time_template ) THEN
  IF ( PRESENT(template) ) fname = 'template name: ' // TRIM(template)
ELSE
  IF ( PRESENT(files) ) THEN
    IF ( SIZE(files) > 0 ) fname = files(1)    ! Give name of first file.
    IF ( SIZE(files) > 1 ) fname = TRIM(fname) // ' and related files'
  END IF
END IF

!*****************************************************************************
! Check that a valid combination of options has been supplied
!*****************************************************************************
!-----------------------------------------------------------------------------
! Check that the start and end times for the data make sense
!-----------------------------------------------------------------------------
! Obviously, the data must start before they end, so bail if not true
IF ( data_start >= data_end )                                                 &
  CALL log_fatal("file_ts_open",                                              &
                 TRIM(fname) // ": " //                                       &
                 "Data start must be strictly before data end.")

! Check that the data period is allowed
SELECT CASE ( data_period )
CASE ( 1:, period_month, period_year )
  ! These are the allowed data periods
  CALL log_info("file_ts_open",                                               &
                "Opening time series with data_period=" //                    &
                TRIM(to_string(data_period)))

CASE DEFAULT
  CALL log_fatal("file_ts_open",                                              &
                 TRIM(fname) // ": " //                                       &
                 "Data period must be > 0 or a 'special' period " //          &
                 "(supplied " // TRIM(to_string(data_period)) // ")")
END SELECT
! Check that the start time is appropriate if using a special period
IF ( data_period == period_month .OR. data_period == period_year ) THEN
  IF ( data_start%time /= 0 .OR. data_start%day /= 1 )                        &
    CALL log_fatal("file_ts_open",                                            &
                   TRIM(fname) // ": " //                                     &
                   "When using data_period=" //                               &
                   TRIM(to_string(data_period)) //                            &
                   ", data must start at 00:00:00 on 1st of month")

  IF ( data_period == period_year .AND. data_start%month /= 1 )               &
    CALL log_fatal("file_ts_open",                                            &
                   TRIM(fname) // ": " //                                     &
                   "When using data_period=" //                               &
                   TRIM(to_string(data_period)) //                            &
                   ", data must start at 00:00:00 on 1st of January")
END IF

!-----------------------------------------------------------------------------
! Check that we have appropriate arguments to be able to locate the
! correct files to use
!-----------------------------------------------------------------------------
IF ( use_time_template ) THEN
  ! If using time templating, we need a template to use
  IF ( .NOT. PRESENT(template) )                                              &
    CALL log_fatal("file_ts_open",                                            &
                   "Time templating selected but no template given")

  ! Check that template is not an empty string
  IF ( LEN_TRIM(template) <= 0 )                                              &
    CALL log_fatal("file_ts_open",                                            &
                   "Time templating selected but template is empty string")

  tpl_period = tpl_detect_period(template)

  ! If the period is not one that has code written for it, abort
  SELECT CASE ( tpl_period )
  CASE ( period_month, period_year )
    ! Fine - just log some information
    CALL log_info("file_ts_open",                                             &
                  "Detected period=" // TRIM(to_string(tpl_period)) //        &
                  " for template " // TRIM(template))

  CASE DEFAULT
    CALL log_fatal("file_ts_open",                                            &
                   "Could not detect supported templating period " //         &
                   "for template " // TRIM(template))
  END SELECT
ELSE
  ! If using lists of files and start times, check that they are present and
  ! consistent
  IF ( .NOT. PRESENT(files) )                                                 &
    CALL log_fatal("file_ts_open",                                            &
                   "Time templating is not selected - a list of files " //    &
                   "must be given")

  IF ( .NOT. PRESENT(file_times) )                                            &
    CALL log_fatal("file_ts_open",                                            &
                   TRIM(fname) // ": " //                                     &
                   "Time templating is not selected - a list of file " //     &
                   "times must be given")

  nfiles = SIZE(files)

  IF ( nfiles <= 0 )                                                          &
    CALL log_fatal("file_ts_open",                                            &
                   "Time templating is not selected - list must contain" //   &
                   " at least one file")

  IF ( SIZE(file_times) /= nfiles )                                           &
    CALL log_fatal("file_ts_open",                                            &
                   TRIM(fname) // ": " //                                     &
                   "'files' and 'file_times' must have the same " //          &
                   "number of entries")

  ! If using yearly data, there must be one file containing all the data
  ! We know from previous checks on file_times(1) and data_start that this file
  ! starts on the 1st Jan for some year
  IF ( data_period == period_year .AND. nfiles > 1 )                          &
    CALL log_fatal("file_ts_open",                                            &
                   TRIM(fname) // ": " //                                     &
                   "Yearly data must be contained in a single file for " //   &
                   "the entirety of the data")

  ! Do checks on the individual files
  DO i = 1,nfiles
    IF ( mode == mode_read ) THEN
      ! By this point, we know that we don't have variable name or time templating,
      ! so we should check if the given files exist
      INQUIRE(FILE = files(i), EXIST = file_exists)
!JWA      IF ( .NOT. file_exists )                                                &
!JWA        CALL log_fatal("file_ts_open",                                        &
!JWA                       "Given file '" // TRIM(files(i)) // "' does not exist")
    ELSE
      ! If we are in write mode, just check that a non-empty file name has been
      ! given
      IF ( LEN_TRIM(files(i)) <= 0 )                                          &
        CALL log_fatal("file_ts_open",                                        &
                       "List of files given, but one or more file names" //   &
                       " are the empty string")
    END IF

    ! Check that the start times of the files are suitable
!     IF ( i == 1 ) THEN
!       ! Check that the first file starts at data start
!       IF ( file_times(i) /= data_start )                                      &
!         CALL log_fatal("file_ts_open",                                        &
!                        TRIM(files(i)) // ": " //                              &
!                        "Start time for the first file must match data_start")
! 
!       prev_file_start = file_times(i)
!     ELSE
!       ! Check that all subsequent files start after the previous file
!       IF ( file_times(i) <= prev_file_start )                                 &
!         CALL log_fatal("file_ts_open",                                        &
!                        TRIM(files(i)) // ": " //                              &
!                        "Files must be given in chronological order")
! 
!       prev_file_start = file_times(i)
!     END IF

    ! If using monthly data, files must start at midnight on the 1st of some
    ! month (but not necessarily only contain data for one month)
    IF ( data_period == period_month ) THEN
      IF ( file_times(i)%time /= 0 .OR. file_times(i)%day /= 1 )              &
        CALL log_fatal("file_ts_open",                                        &
                       TRIM(files(i)) // ": " //                              &
                       "When using monthly data with a list of files, " //    &
                       "all files must start at 00:00:00 on 1st of " //       &
                       "some month")
    END IF
  END DO
END IF  !  use_time_template

!-----------------------------------------------------------------------------
! If we have a climatology, check that the given arguments and detected
! properties make sense
!-----------------------------------------------------------------------------
IF ( is_climatology ) THEN
  ! It makes no sense to specify is_climatology if in write mode
  IF ( mode == mode_write )                                                   &
    CALL log_fatal("file_ts_open",                                            &
                   "Cannot open a file in write mode as a climatology.")

  ! The data must start at the beginning of the year
  IF ( data_start%time /= 0 .OR.                                              &
       data_start%day /= 1 .OR. data_start%month /= 1 )                       &
    CALL log_fatal("file_ts_open",                                            &
                   TRIM(fname) // ": " //                                     &
                   "When using data as a climatology, data must start " //    &
                   "at 00:00:00 on 1st of January")

  ! The data must apply for exactly one year
  IF ( datetime_advance(data_start, period_year) /= data_end )                &
    CALL log_fatal("file_ts_open",                                            &
                   TRIM(fname) // ": " //                                     &
                   "When using data as a climatology, exactly one year " //   &
                   "of data must be given.")

  ! Warn the user if they have supplied yearly data - this means the same data
  ! will be used every year, and is the same as having fixed data
  IF ( data_period == period_year )                                           &
    CALL log_warn("file_ts_open",                                             &
                  "When using a climatology, a data period of a year " //     &
                  "is equivalent to having fixed data")
END IF

!-----------------------------------------------------------------------------
! Check that the MPI variables are either specified together or not at all
!-----------------------------------------------------------------------------
IF ( PRESENT(comm) .NEQV. PRESENT(info) )                                     &
  CALL log_fatal("file_ts_open",                                              &
                 "Only one of comm and info is present - either give a " //   &
                 "value for both MPI variables for parallel access or " //    &
                 "omit both for serial access")


!-----------------------------------------------------------------------------
! Now we are happy that the arguments are consistent, we can set up the
! file_ts object
!-----------------------------------------------------------------------------
FILE%mode = mode

! Store info about MPI-IO settings
IF ( PRESENT(comm) ) THEN
  FILE%use_mpiio = .TRUE.
  FILE%comm = comm
  FILE%info = info
END IF

FILE%data_start  = data_start
FILE%data_end    = data_end
FILE%data_period = data_period
FILE%is_climatology = is_climatology

FILE%use_time_template = use_time_template

IF ( use_time_template ) THEN
  ! With time templating
  FILE%template   = template
  FILE%tpl_period = tpl_period
ELSE
  ! With a file list
  FILE%nfiles = nfiles

  ! Allocate space for the files and file_times arrays
  ALLOCATE(FILE%files(nfiles))
  ALLOCATE(FILE%file_times(nfiles))

  FILE%files(:)      = files(:)
  FILE%file_times(:) = file_times(:)
END IF

RETURN

END FUNCTION file_ts_open
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ts_def_grid(FILE, grid)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Defines the grid used by variables in the file
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ts), INTENT(INOUT) :: FILE
    ! The file to define the grid on

TYPE(grid_info), INTENT(IN) :: grid  ! The grid to define


!-----------------------------------------------------------------------------


! If we are not in define mode, error out
IF ( .NOT. FILE%define_mode )                                                 &
  CALL log_fatal("file_ts_def_grid",                                          &
                 "Cannot define grid - file is not in define mode")

FILE%grid = grid


RETURN

END SUBROUTINE file_ts_def_grid
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_ts_def_dim(FILE, dim_name, dim_len) RESULT(dim_id)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Defines a dimension on the given timeseries file, returning the dimension id
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ts), INTENT(INOUT) :: FILE
                                ! The file to define the dimension on
CHARACTER(LEN=*), INTENT(IN) :: dim_name
                                ! The name of the dimension
INTEGER, INTENT(IN) :: dim_len  ! The length of the dimension

! Return type
INTEGER :: dim_id               ! The dimension id


!-----------------------------------------------------------------------------

! If we are not in define mode, error out
IF ( .NOT. FILE%define_mode )                                                 &
  CALL log_fatal("file_ts_def_dim",                                           &
                 "Cannot define dimension - file is not in define mode")

! If adding another dimension will cause us to have too many dimensions,
! error out
IF ( FILE%ndims >= max_dim_file )                                             &
  CALL log_fatal("file_ts_def_dim",                                           &
                 "Too many dimensions in file - try increasing max_dim_file")

!-----------------------------------------------------------------------------
! Store the dimension attributes so that they can be used later to define
! dimensions on actual file(s)
!-----------------------------------------------------------------------------
FILE%ndims = FILE%ndims + 1

! The returned dimension id is just the index in the dims array on the file object
dim_id = FILE%ndims

FILE%dims(dim_id)%NAME   = dim_name
FILE%dims(dim_id)%length = dim_len

RETURN

END FUNCTION file_ts_def_dim
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_ts_def_time_dim(FILE, dim_name) RESULT(dim_id)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Defines a dimension on the given file, returning the dimension id
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ts), INTENT(INOUT) :: FILE
                                ! The file to define the dimension on
CHARACTER(LEN=*), INTENT(IN) :: dim_name
                                ! The name of the dimension

! Return type
INTEGER :: dim_id               ! The dimension id


!-----------------------------------------------------------------------------

! If we are not in define mode, error out
IF ( .NOT. FILE%define_mode )                                                 &
  CALL log_fatal("file_ts_def_time_dim",                                      &
                 "Cannot define time dimension - file is not in define mode")

! A file can only have one time dimension
IF ( FILE%has_time_dim )                                                      &
  CALL log_fatal("file_ts_def_time_dim",                                      &
                 "Time dimension has already been defined")

!-----------------------------------------------------------------------------
! Store the dimension attributes so that they can be used later to define
! the dimension on actual file(s)
!-----------------------------------------------------------------------------
FILE%has_time_dim = .TRUE.

FILE%time_dim%NAME   = dim_name
FILE%time_dim%length = -1  ! The time dimension does not have a length

! Return -1 as the dimension id, since the return value should never be used
! and -1 can never be a valid dimension id for a non-time dimension
dim_id = -1

RETURN

END FUNCTION file_ts_def_time_dim
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_ts_def_var(FILE, var_name, levels_dims, use_time) RESULT(var_id)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Defines a variable in the given file, returning the variable id
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ts), INTENT(INOUT) :: FILE
                                ! The file to define the variable in
CHARACTER(LEN=*), INTENT(IN) :: var_name
                                  ! The name of the variable
INTEGER, INTENT(IN), OPTIONAL :: levels_dims(:)
                                  ! The ids of the dimensions to use for the
                                  ! vertical levels of the variable, if
                                  ! required
                                  ! If not given or an array of 0 length is
                                  ! given, no vertical levels are used
LOGICAL, INTENT(IN) :: use_time   ! Indicates whether the variable uses the
                                  ! time dimension

! Return type
INTEGER :: var_id               ! The variable id


! Work variables


!-----------------------------------------------------------------------------

! If we are not in define mode, error out
IF ( .NOT. FILE%define_mode )                                                 &
  CALL log_fatal("file_ts_def_var",                                           &
                 "Cannot define variable - file is not in define mode")

! If adding another variable will cause us to have too many variables,
! error out
IF ( FILE%nvars >= max_var_file )                                             &
  CALL log_fatal("file_ts_def_var",                                           &
                 "Too many variables in file - try increasing max_var_file")


!-----------------------------------------------------------------------------
! Store the variable attributes so that they can be used later to define
! variables on actual file(s)
!-----------------------------------------------------------------------------

FILE%nvars = FILE%nvars + 1

! The returned variable id is just the index in the vars array on the file object
var_id = FILE%nvars

FILE%vars(var_id)%NAME       = var_name
FILE%vars(var_id)%use_time   = use_time

! This exploits the fact that providing a levels_dims array of size 0 to
! file_gridded_def_var is the same as not providing it
IF ( PRESENT(levels_dims) ) THEN
  ALLOCATE(FILE%vars(var_id)%levels_dims(SIZE(levels_dims)))
  FILE%vars(var_id)%levels_dims(:) = levels_dims(:)
ELSE
  ALLOCATE(FILE%vars(var_id)%levels_dims(0))
END IF

RETURN

END FUNCTION file_ts_def_var
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ts_def_attr_real(FILE, var_id, NAME, VALUE)

USE io_constants, ONLY: attr_global

USE dictionary_mod, ONLY: dict_create, dict_set

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Define a real valued attribute on the given variable with the given name
!   and value
!   To define a global attribute, specify attr_global as var_id
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ts), INTENT(INOUT) :: FILE
                                  ! The file to define the attribute in
INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                  ! attribute on
CHARACTER(LEN=*), INTENT(IN) :: NAME
                                  ! The name of the attribute
REAL, INTENT(IN) :: VALUE         ! The value of the attribute


!-----------------------------------------------------------------------------


! If we are not in define mode, error out
IF ( .NOT. FILE%define_mode )                                                 &
  CALL log_fatal("file_ts_def_attr_real",                                     &
                 "Cannot define attribute - file is not in define mode")

!-----------------------------------------------------------------------------
! Work out what dictionary we want to populate
!-----------------------------------------------------------------------------
IF ( var_id == attr_global ) THEN
  ! If it is a global attribute, use the global dictionary
  ! Create the required dictionary if it has not been used yet
  IF ( FILE%attrs_real%length == 0 )                                          &
    FILE%attrs_real = dict_create(max_attr_file, 1.0)

  CALL dict_set(FILE%attrs_real, NAME, VALUE)
ELSE
  ! Otherwise, use the dictionary on the specified variable
  ! Create the required dictionary if it has not been used yet
  IF ( FILE%vars(var_id)%attrs_real%length == 0 )                             &
    FILE%vars(var_id)%attrs_real = dict_create(max_attr_file, 1.0)

  CALL dict_set(FILE%vars(var_id)%attrs_real, NAME, VALUE)
END IF

RETURN

END SUBROUTINE file_ts_def_attr_real


SUBROUTINE file_ts_def_attr_int(FILE, var_id, NAME, VALUE)

USE io_constants, ONLY: attr_global

USE dictionary_mod, ONLY: dict_create, dict_set

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Define an integer valued attribute on the given variable with the given
!   name and value
!   To define a global attribute, specify attr_global as var_id
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ts), INTENT(INOUT) :: FILE
                                  ! The file to define the attribute in
INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                  ! attribute on
CHARACTER(LEN=*), INTENT(IN) :: NAME
                                  ! The name of the attribute
INTEGER, INTENT(IN) :: VALUE      ! The value of the attribute


!-----------------------------------------------------------------------------

! If we are not in define mode, error out
IF ( .NOT. FILE%define_mode )                                                 &
  CALL log_fatal("file_ts_def_attr_int",                                      &
                 "Cannot define attribute - file is not in define mode")

!-----------------------------------------------------------------------------
! Work out what dictionary we want to populate
!-----------------------------------------------------------------------------
IF ( var_id == attr_global ) THEN
  ! If it is a global attribute, use the global dictionary
  ! Create the required dictionary if it has not been used yet
  IF ( FILE%attrs_int%length == 0 )                                           &
    FILE%attrs_int = dict_create(max_attr_file, INT(1))

  CALL dict_set(FILE%attrs_int, NAME, VALUE)
ELSE
  ! Otherwise, use the dictionary on the specified variable
  ! Create the required dictionary if it has not been used yet
  IF ( FILE%vars(var_id)%attrs_int%length == 0 )                              &
    FILE%vars(var_id)%attrs_int = dict_create(max_attr_file, INT(1))

  CALL dict_set(FILE%vars(var_id)%attrs_int, NAME, VALUE)
END IF

RETURN

END SUBROUTINE file_ts_def_attr_int


SUBROUTINE file_ts_def_attr_char(FILE, var_id, NAME, VALUE)

USE io_constants, ONLY: attr_global

USE dictionary_mod, ONLY: dict_create, dict_set

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Define a character valued attribute on the given variable with the given
!   name and value
!   To define a global attribute, specify attr_global as var_id
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ts), INTENT(INOUT) :: FILE
                                  ! The file to define the attribute in
INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                  ! attribute on
CHARACTER(LEN=*), INTENT(IN) :: NAME
                                  ! The name of the attribute
CHARACTER(LEN=*), INTENT(IN) :: VALUE
                                  ! The value of the attribute


!-----------------------------------------------------------------------------

! If we are not in define mode, error out
IF ( .NOT. FILE%define_mode )                                                 &
  CALL log_fatal("file_ts_def_attr_char",                                     &
                 "Cannot define attribute - file is not in define mode")

!-----------------------------------------------------------------------------
! Work out what dictionary we want to populate
!-----------------------------------------------------------------------------
IF ( var_id == attr_global ) THEN
  ! If it is a global attribute, use the global dictionary
  ! Create the required dictionary if it has not been used yet
  IF ( FILE%attrs_char%length == 0 )                                          &
    FILE%attrs_char = dict_create(max_attr_file, NAME)

  CALL dict_set(FILE%attrs_char, NAME, VALUE)
ELSE
  ! Otherwise, use the dictionary on the specified variable
  ! Create the required dictionary if it has not been used yet
  IF ( FILE%vars(var_id)%attrs_char%length == 0 )                             &
    FILE%vars(var_id)%attrs_char = dict_create(max_attr_file, NAME)

  CALL dict_set(FILE%vars(var_id)%attrs_char, NAME, VALUE)
END IF

RETURN

END SUBROUTINE file_ts_def_attr_char
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ts_enddef(FILE)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes the file out of definition mode - no more dimensions or variables
!   may be defined after this
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ts), INTENT(INOUT) :: FILE


!-----------------------------------------------------------------------------

FILE%define_mode = .FALSE.


! Verify that a time dimension has been defined - this is necessary for a
! timeseries file
IF ( .NOT. FILE%has_time_dim )                                                &
  CALL log_fatal("file_ts_enddef",                                            &
                 "Time dimension has not been defined")

! Seek the file so that the first timestep will be read/written when requested
CALL file_ts_seek_to_datetime(FILE, FILE%data_start)

RETURN

END SUBROUTINE file_ts_enddef
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ts_seek_to_datetime(FILE, dt)

USE io_constants, ONLY: mode_read, mode_write

USE precision_mod, ONLY: int64

USE datetime_utils_mod, ONLY: days_in_month

USE datetime_mod, ONLY: period_month, period_year, l_360, l_leap,             &
                         datetime_create, datetime_clone, datetime_advance,   &
                         datetime_diff, datetime_to_string

USE templating_mod, ONLY: tpl_substitute_datetime

USE file_gridded_mod, ONLY: file_gridded_seek

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Seeks the file_ts to a particular datetime (i.e. the next time values are
!   read from the file_ts using file_ts_read_var, they will be read at
!   the requested datetime)
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ts), INTENT(INOUT) :: FILE  ! The file to seek
TYPE(datetime), INTENT(IN) :: dt        ! The datetime to seek to


! Work variables
TYPE(datetime) :: dt_local  ! Local version of dt that is possibly modified
                            ! to deal with climatology

CHARACTER(LEN=max_file_name_len) :: file_name ! The name of the file
                                              ! containing the requested
                                              ! timestep
TYPE(datetime) :: file_start, next_file_start ! The start times of the file
                                              ! containing the requested
                                              ! timestep and the next file
INTEGER :: file_index                         ! The index of the file
                                              ! containing the requested
                                              ! timestep in the files
                                              ! array on the file object
INTEGER :: tsteps_from_start                  ! How far into the file is
INTEGER(KIND=int64) :: secs_from_start      ! the requested timestep?

INTEGER :: i  ! Loop counters


!-----------------------------------------------------------------------------

! If file is still in define mode, error out
IF ( FILE%define_mode )                                                       &
  CALL log_fatal("file_ts_seek_to_datetime",                                  &
                 "Cannot seek while file is in define mode")

!-----------------------------------------------------------------------------
! Arbitrary seeks are only allowed in read mode
! Only seeking to the beginning of the data is permitted in write mode
!-----------------------------------------------------------------------------
IF ( FILE%mode /= mode_read .AND. dt /= FILE%data_start )                     &
  CALL log_fatal("file_ts_seek_to_datetime",                                  &
                 "Seeking to arbitrary datetime is only allowed in " //       &
                 "read mode")

!-----------------------------------------------------------------------------
! If we are using a climatology, create a new datetime to seek to where
! the year is the same as the year that the data runs for
! Otherwise we seek to the supplied datetime
!-----------------------------------------------------------------------------
dt_local = datetime_clone(dt)
IF ( FILE%is_climatology ) THEN
  dt_local%year = FILE%data_start%year
  ! We might need to adjust the day - if dt is in a leap year but the year we
  ! have data for is not a leap year and we are seekeing to somewhere in Feb 29
  dt_local%day = MIN(                                                         &
    dt_local%day, days_in_month(dt_local%year, dt_local%month, l_360, l_leap) &
  )
END IF

!-----------------------------------------------------------------------------
! Check that the datetime to seek to is appropriate for special periods
!-----------------------------------------------------------------------------
! Monthly data is assumed to apply at 00:00:00 on the first of the month
! Annual data adds the additional constraint that the month must be January
IF ( FILE%data_period == period_month .OR.                                    &
     FILE%data_period == period_year       ) THEN

  IF ( dt_local%time /= 0 .OR. dt_local%day /= 1 )                            &
    CALL log_fatal("file_ts_seek_to_datetime",                                &
                   "Cannot seek to datetime - data_period=" //                &
                   TRIM(to_string(FILE%data_period)) //                       &
                   " assumes data apply at 00:00:00 on 1st of month " //      &
                   "(given " // datetime_to_string(dt) // ")")

  IF ( FILE%data_period == period_year .AND. dt_local%month /= 1 )            &
    CALL log_fatal("file_ts_seek_to_datetime",                                &
                   "Cannot seek to datetime - data_period=" //                &
                   TRIM(to_string(FILE%data_period)) //                       &
                   " assumes data apply at 00:00:00 on 1st of January " //    &
                   "(given " // datetime_to_string(dt) // ")")

END IF

!-----------------------------------------------------------------------------
! Check that it is actually possible to seek to the requested datetime
!-----------------------------------------------------------------------------
IF ( dt_local < FILE%data_start .OR. FILE%data_end <= dt_local )              &
  CALL log_fatal("file_ts_seek_to_datetime",                                  &
                 "No data for datetime - out of range")

!-----------------------------------------------------------------------------
! We need to get the name and start time of the file containing the datetime
! that has been asked for
! This is done differently depending on whether we are using time templating
! or not
!-----------------------------------------------------------------------------
! Initialise file_index - note that this is never used when time templating
! is in use
file_index = 0

IF ( FILE%use_time_template ) THEN
  ! Time templating is in use
  ! We already know that the requested datetime is within the data range, so
  ! get the name of the file to open by substituting values from the given
  ! datetime into the template
  file_name = tpl_substitute_datetime(FILE%template, dt_local)

  ! Get the start and end time for the file based on the constraints imposed
  ! by the templating period being used
  SELECT CASE ( FILE%tpl_period )
  CASE ( period_month )
    ! The first data in a monthly file must apply at 00:00:00 on the 1st of the
    ! month
    file_start = datetime_create(                                             &
      dt_local%year, dt_local%month, 1, 0, 0, 0                               &
    )

  CASE ( period_year )
    ! The first data in a yearly file must apply at 00:00:00 on the 1st of Jan
    file_start = datetime_create(dt_local%year, 1, 1, 0, 0, 0)

  CASE DEFAULT
    CALL log_fatal("file_ts_seek_to_datetime",                                &
                   "No code for tpl_period=" //                               &
                   TRIM(to_string(FILE%tpl_period)))
  END SELECT

  ! In all cases, the start of the next file is one templating period on from
  ! the start of the file we are opening
  next_file_start = datetime_advance(file_start, FILE%tpl_period)
ELSE
  ! No time templating - use the list of files instead

  ! Assume we are in the last file unless we discover otherwise while scanning
  ! below
  file_index = FILE%nfiles

  ! Get the index of the file containing the correct date and time
  DO i = 1,FILE%nfiles-1
    IF ( FILE%file_times(i) <= dt_local .AND.                                 &
         dt_local < FILE%file_times(i+1)     ) THEN
      file_index = i
      EXIT
    END IF
  END DO

  file_name  = FILE%files(file_index)
  file_start = FILE%file_times(file_index)

  ! Get the start time for the next file
  ! If we are in the last file (i.e. there is no next file), then set this to
  ! 1 second past the end of the data so it never gets hit
  IF ( file_index == FILE%nfiles) THEN
    next_file_start = datetime_advance(FILE%data_end, 1)
  ELSE
    next_file_start = FILE%file_times(file_index+1)
  END IF
END IF

!-----------------------------------------------------------------------------
! Calculate where in the new file we need to seek to
!-----------------------------------------------------------------------------
SELECT CASE ( FILE%data_period )
CASE ( 1: )
  ! Any positive number is interpreted as a period in seconds

  !-----------------------------------------------------------------------------
  ! Get the timestep in the file that we want to seek to
  !-----------------------------------------------------------------------------
  secs_from_start = datetime_diff(file_start, dt_local)

  ! Check that the requested time is a whole number of timesteps into the file
  IF ( MOD(secs_from_start, INT(FILE%data_period, int64)) /= 0 )              &
    CALL log_fatal("file_ts_seek_to_datetime",                                &
                   "Cannot seek to datetime " //                              &
                   datetime_to_string(dt) //                                  &
                   " - it is not a valid timestep in file " //                &
                   TRIM(file_name))

  ! Integer division gives us the number of timesteps into the file that we
  ! need for the requested time
  ! To avoid cluttering the code with INT calls, we rely on implicit conversion
  ! betweeen int64 and the default integer kind to do the right thing, i.e.
  !    file%data_period => int64
  !    division in 64-bit,
  !    result => default real kind
  ! For gfortran, we can confirm this using -Wconversion
  ! For intel, results confirm that this is what is happening
  tsteps_from_start = ( secs_from_start / FILE%data_period ) + 1

CASE ( period_month )
  ! Monthly data is only allowed when every file starts at 00:00:00 on the 1st
  ! of some month - when using a list of files this is verified in
  ! file_ts_open, and when using time templating it is true by construction
  ! This means we can calculate the timestep in the file to use just by
  ! comparing the year and month of file_start and dt
  tsteps_from_start = ( ( dt_local%year - file_start%year ) * 12 )            &
                    + ( dt_local%month - file_start%month )                   &
                    + 1

CASE ( period_year )
  ! Yearly data is only allowed when there is one file for the whole span of
  ! the data and the data starts at 00:00:00 on 1st January for some year
  ! - this is verified in file_ts_open
  ! This means we can calculate the timestep in the file to use just by comparing
  ! the year of file_start and dt
  tsteps_from_start = ( dt_local%year - file_start%year ) + 1

CASE DEFAULT
  CALL log_fatal("file_ts_seek_to_datetime",                                  &
                 "No code for data_period = " //                              &
                 TRIM(to_string(FILE%data_period)))
END SELECT

!-----------------------------------------------------------------------------
! Open the file that we have found
!-----------------------------------------------------------------------------
CALL file_ts_internal_open_file(FILE, file_name)


! If we are in read mode, then seek the file to the correct record
! If we are in write mode, then we know (due to the check at the top of the
! file) that we are seeking to the start of the file(s) and so no further
! seeking is necessary
IF ( FILE%mode == mode_read)                                                  &
  CALL file_gridded_seek(FILE%open_file, tsteps_from_start)

! Update properties on the file_ts object
FILE%current_datetime = dt_local
FILE%open_file_index  = file_index
FILE%next_file_start  = next_file_start

! Indicate that the timestamp for this file timestep has not been written yet
FILE%timestamp_written = .FALSE.

RETURN

END SUBROUTINE file_ts_seek_to_datetime
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ts_advance(FILE)

USE io_constants, ONLY: mode_write

USE datetime_mod, ONLY: datetime_advance, datetime_to_string

USE templating_mod, ONLY: tpl_substitute_datetime

USE file_gridded_mod, ONLY: file_gridded_advance

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Advances the file by one timestep
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument type
TYPE(file_ts), INTENT(INOUT) :: FILE


! Work variables
TYPE(datetime) :: next_time  ! The datetime that we are trying to advance to

! These are used only if advancing requires opening a new file
CHARACTER(LEN=max_file_name_len) :: file_name  ! The name of the file
INTEGER :: file_index  ! The index of the file in the files array of the
                       ! file_ts object
TYPE(datetime) :: next_file_start  ! The start time of the next file after
                                   ! the one we need to open


!-----------------------------------------------------------------------------


! If file is still in define mode, error out
IF ( FILE%define_mode )                                                       &
  CALL log_fatal("file_ts_advance",                                           &
                 "Cannot seek while file is in define mode")


next_time = datetime_advance(                                                 &
  FILE%current_datetime, FILE%data_period                                     &
)

IF ( next_time >= FILE%data_end ) THEN
  IF ( FILE%is_climatology ) THEN
    ! If we are using a climatology, then advancing past the end of the data
    ! just means we go back to the start
    CALL file_ts_seek_to_datetime(FILE, FILE%data_start)
    RETURN
  ELSE
    ! If we are not using a climatology, stop with an error as we cannot advance
    ! beyond the end of the data
    CALL log_fatal("file_ts_advance",                                         &
                   "Cannot advance to " // datetime_to_string(next_time) //   &
                   " - data_end is " // datetime_to_string(FILE%data_end))
  END IF
END IF

! Check if we need to open a new file
IF ( next_time >= FILE%next_file_start ) THEN
  !-----------------------------------------------------------------------------
  ! We need to get the name of the file to open and the start time of the next
  ! file so we can update the file_ts object
  !
  ! We know that we will be using the first timestep in the file, so we don't
  ! need to calculate where in the file to start
  !
  ! This is done differently depending on whether we are using time templating
  ! or not
  !-----------------------------------------------------------------------------
  ! Initialise file_index - not used when time templating is in use
  file_index = 0

  IF ( FILE%use_time_template ) THEN
    ! Time templating is in use
    ! We already know that the next time is within the data range, so we get the
    ! file to open by just substituting values into the template
    file_name = tpl_substitute_datetime(FILE%template, next_time)
    ! We get the start time for the next file by skipping one templating period
    ! on from the previous value
    next_file_start = datetime_advance(                                       &
      FILE%next_file_start, FILE%tpl_period                                   &
    )
  ELSE
    ! No time templating - use the list of files instead

    ! We know that the file we are moving out of is not the last file, otherwise
    ! the test on file%data_end earlier would fail
    file_index = FILE%open_file_index + 1

    ! So we can easily get the name of the file now
    file_name = FILE%files(file_index)

    ! Get the start time for the next file
    ! If we are in the moving in to the last file (i.e. there is no next file),
    ! then set this to 1 second past the end of the data so it never gets hit
    IF ( file_index == FILE%nfiles ) THEN
      next_file_start = datetime_advance(FILE%data_end, 1)
    ELSE
      next_file_start = FILE%file_times(file_index+1)
    END IF
  END IF

  ! Make sure we have the correct file open
  CALL file_ts_internal_open_file(FILE, file_name)

  ! Update open file related properties on the file_ts object
  FILE%open_file_index  = file_index
  FILE%next_file_start  = next_file_start

ELSE
  !-----------------------------------------------------------------------------
  ! No new file required - all we have to do is advance the current file
  !-----------------------------------------------------------------------------
  CALL file_gridded_advance(FILE%open_file)
END IF

! We always update the current datetime
FILE%current_datetime = next_time


! Indicate that the timestamp for this file timestep has not been written yet
FILE%timestamp_written = .FALSE.

RETURN

END SUBROUTINE file_ts_advance
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_ts_read_var(FILE, var_id, extract_subgrid, subgrid) RESULT(cube)

USE grid_utils_mod, ONLY: subgrid_info

USE file_gridded_mod, ONLY: file_gridded_read_var

USE data_cube_mod, ONLY: data_cube

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads data from the given variable in the given file as a data cube
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ts), INTENT(INOUT) :: FILE
                               ! The file to read data from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
LOGICAL, INTENT(IN) :: extract_subgrid
                               ! T - extract a subgrid to return from the
                               !     full grid of the file
                               ! F - return the full grid
TYPE(subgrid_info), OPTIONAL :: subgrid  ! The subgrid to extract


! Return type
TYPE(data_cube) :: cube  ! The data cube read from file


! Work variables
INTEGER :: actual_var_id  ! The id of the variable in the underlying file


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Defer the actual reading to the file_gridded routine
!-----------------------------------------------------------------------------
actual_var_id = FILE%vars(var_id)%current_id

IF ( PRESENT(subgrid) ) THEN
  cube = file_gridded_read_var(                                               &
    FILE%open_file, actual_var_id, extract_subgrid, subgrid                   &
  )
ELSE
  cube = file_gridded_read_var(                                               &
    FILE%open_file, actual_var_id, extract_subgrid                            &
  )
END IF

RETURN

END FUNCTION file_ts_read_var
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ts_write_var(FILE, var_id, cube, write_subgrid, subgrid)

USE precision_mod, ONLY: int64

USE grid_utils_mod, ONLY: subgrid_info

USE file_gridded_mod, ONLY: file_gridded_write_var

USE datetime_mod, ONLY: datetime_advance, datetime_diff

USE file_mod, ONLY: file_write_var

USE data_cube_mod, ONLY: data_cube

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Writes data given as a data cube to the given variable in the given file
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ts), INTENT(INOUT) :: FILE
                               ! The file to write data to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
TYPE(data_cube), INTENT(IN) :: cube
                               ! The values to write to the file, contained
                               ! in a data cube
                               ! The first 2 dimensions of the cube should
                               ! be the grid dimensions
                               ! The rest of the dimensions will be taken to
                               ! be levels dimensions
LOGICAL, INTENT(IN) :: write_subgrid
                               ! T - write to a subgrid of the full
                               !     grid of the file
                               ! F - write the full grid
TYPE(subgrid_info), OPTIONAL :: subgrid  ! The subgrid to write


! Work variables
INTEGER :: actual_var_id  ! The id of the variable in the underlying file

INTEGER(KIND=int64) :: interval_start, interval_end
                          ! The start and end time of the current interval
                          ! in seconds since the start of the data file


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! If this is the first write since an advance, then write the timestamp
! information
!
! If there are no writes between advancing a file and closing it, nothing
! will be written for the final time
!-----------------------------------------------------------------------------
IF ( .NOT. FILE%timestamp_written ) THEN
  ! Work out the start and end of the current interval
  interval_start = datetime_diff(FILE%data_start, FILE%current_datetime)
  interval_end = datetime_diff(                                               &
    FILE%data_start,                                                          &
    datetime_advance(FILE%current_datetime, FILE%data_period)                 &
  )

  ! Write the index variable to be the end of the current interval
  CALL file_write_var(                                                        &
    FILE%open_file%fh, FILE%time_index_var_id, REAL(interval_end)             &
  )

  ! Write the bounds variable to indicate how long the interval lasts for
  CALL file_write_var(                                                        &
    FILE%open_file%fh, FILE%time_bounds_var_id,                               &
    REAL((/ interval_start, interval_end /))                                  &
  )

  FILE%timestamp_written = .TRUE.
END IF

!-----------------------------------------------------------------------------
! Defer the actual writing to the file_gridded routine
!-----------------------------------------------------------------------------
actual_var_id = FILE%vars(var_id)%current_id

IF ( PRESENT(subgrid) ) THEN
  CALL file_gridded_write_var(                                                &
    FILE%open_file, actual_var_id, cube, write_subgrid, subgrid               &
  )
ELSE
  CALL file_gridded_write_var(                                                &
    FILE%open_file, actual_var_id, cube, write_subgrid                        &
  )
END IF

RETURN

END SUBROUTINE file_ts_write_var
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ts_close(FILE)

USE file_gridded_mod, ONLY: file_gridded_close

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Closes and frees any resources consumed by the given file object
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ts), INTENT(INOUT) :: FILE ! The file to close


! Work variables
INTEGER :: i  ! Index variable


!-----------------------------------------------------------------------------

! Make sure the currently open file is closed
CALL file_gridded_close(FILE%open_file)

! Free the memory associated with the files and file_times arrays if it was
! allocated (i.e. if not using time templating)
IF ( .NOT. FILE%use_time_template ) THEN
  DEALLOCATE(FILE%files)
  NULLIFY(FILE%files)

  DEALLOCATE(FILE%file_times)
  NULLIFY(FILE%file_times)
END IF

! Free the memory associated with the levels dimensions for variables
DO i = 1,FILE%nvars
  IF ( ASSOCIATED(FILE%vars(i)%levels_dims) ) THEN
    DEALLOCATE(FILE%vars(i)%levels_dims)
    NULLIFY(FILE%vars(i)%levels_dims)
  END IF
END DO

RETURN

END SUBROUTINE file_ts_close

! Utility routines
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ts_internal_open_file(FILE, file_name)

USE io_constants, ONLY: mode_write, attr_global

USE datetime_mod, ONLY: l_360, l_leap, datetime_to_string

USE dictionary_mod, ONLY: dict_key_len, dict_char_val_len, dict_get

USE data_cube_mod, ONLY: data_cube

USE file_gridded_mod, ONLY: file_gridded_open, file_gridded_def_grid,         &
                             file_gridded_def_dim,                            &
                             file_gridded_def_record_dim,                     &
                             file_gridded_def_var, file_gridded_def_attr,     &
                             file_gridded_enddef, file_gridded_close

USE file_mod, ONLY: file_def_dim, file_def_var, file_def_attr

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   INTERNAL PROCEDURE TO file_ts_mod
!   Sets the currently open file to that specified by file_name, including
!   defining dimensions, variables etc.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ts), INTENT(INOUT) :: FILE  ! The file to update
CHARACTER(LEN=*), INTENT(IN) :: file_name  ! The file to open and set as
                                           ! current file


! Work variables
INTEGER :: attr_var_id    ! Variable ID to use when defining attributes

CHARACTER(LEN=dict_key_len) :: key  ! Used when iterating over attribute
INTEGER :: int_val                  ! dictionaries
REAL :: real_val
CHARACTER(LEN=dict_char_val_len) :: char_val

INTEGER, ALLOCATABLE :: zero_element_array(:)  ! This array is allocated
                                               ! to have size 0 - used to
                                               ! indicate passing no
                                               ! dimensions to file_def_var
                                               ! for the time index


INTEGER :: n_non_time_vars  ! The number of non-time-varying variables
INTEGER :: nx, ny  ! Number of x and y points in the grid - used to
                   ! dimension array below
TYPE(data_cube), ALLOCATABLE :: non_time_varying_data(:)
                         ! Array to hold the data cubes for non-time-varying
                         ! variables from the previous file until they can
                         ! be copied into the next file
INTEGER :: var_nlevs  ! The size of the vertical levels dimension for the
                      ! variable currently being processed
                      ! Only use when copying non-time-varying variables
                      ! from one file to another

INTEGER :: dim_id  ! Id of dimension for time_bounds

INTEGER :: i,j  ! Loop counters


!-----------------------------------------------------------------------------


IF ( FILE%has_open_file ) THEN
  IF ( FILE%mode == mode_write ) THEN
    !-----------------------------------------------------------------------------
    ! If we are in write mode, we want to copy the values of all the variables
    ! not using the time dimensions into the new file
    ! To do that, we have to grab their values here
    !-----------------------------------------------------------------------------
    ! Only bother doing this if we have non-time-varying variables
    n_non_time_vars = COUNT( .NOT. FILE%vars(1:FILE%nvars)%use_time)
    IF ( n_non_time_vars > 0 ) THEN
      ! We allocate (potentially) too much data, since we only need it for this routine
      ALLOCATE(non_time_varying_data(n_non_time_vars))

      ! Gather the data for the non-time varying variables
      j = 1
      DO i = 1,FILE%nvars
        IF ( .NOT. FILE%vars(i)%use_time ) THEN
          ! We are not using a subgrid, so give extract_subgrid=F
          non_time_varying_data(j) = file_ts_read_var(FILE, i, .FALSE.)
          ! Advance the non-time variable counter
          j = j + 1
        END IF
      END DO
    END IF
  END IF

  !-----------------------------------------------------------------------------
  ! Close the currently open file
  !-----------------------------------------------------------------------------
  CALL file_gridded_close(FILE%open_file)
END IF

!-----------------------------------------------------------------------------
! Open the new file
!-----------------------------------------------------------------------------
IF ( FILE%use_mpiio ) THEN
  ! If we are using parallel I/O, then pass MPI variables to file_gridded_open
  FILE%open_file = file_gridded_open(                                         &
    file_name, FILE%mode, FILE%comm, FILE%info                                &
  )
ELSE
  ! Otherwise omit the MPI variables
  FILE%open_file = file_gridded_open(file_name, FILE%mode)
END IF
FILE%has_open_file = .TRUE.

!-----------------------------------------------------------------------------
! Define the grid on the newly opened file
!-----------------------------------------------------------------------------
CALL file_gridded_def_grid(FILE%open_file, FILE%grid)

!-----------------------------------------------------------------------------
! Define the vertical level dimensions and update the dimension ids to the
! ones used by the newly opened file
!-----------------------------------------------------------------------------
DO i = 1,FILE%ndims
  FILE%dims(i)%current_id = file_gridded_def_dim(                             &
    FILE%open_file, FILE%dims(i)%NAME, FILE%dims(i)%length                    &
  )
END DO

!-----------------------------------------------------------------------------
! Define the time dimension as a record dimension on the open file and update
! the stored id to the one used by the newly opened file
!-----------------------------------------------------------------------------
FILE%time_dim%current_id = file_gridded_def_record_dim(                       &
  FILE%open_file, FILE%time_dim%NAME                                          &
)

!-----------------------------------------------------------------------------
! Define the time indexes if we are in write mode
!
! We define a variable time_bounds that defines the start and end of each
! "time cell". This is associated with the time coordinate using the
! "bounds" attribute.
!-----------------------------------------------------------------------------
! Since the time bounds and index do not want grid semantics, we have to go
! down to the raw file routines
IF ( FILE%mode == mode_write ) THEN
  !-----------------------------------------------------------------------------
  ! Create the time_bounds variable
  !-----------------------------------------------------------------------------
  ! Before creating the variable, we need to create a dimension of size 2 to use
  dim_id = file_def_dim(FILE%open_file%fh, "nt", 2)
  FILE%time_bounds_var_id = file_def_var(                                     &
    FILE%open_file%fh, "time_bounds", (/ dim_id /), .TRUE.                    &
  )

  ! Define it's attributes
  CALL file_def_attr(                                                         &
    FILE%open_file%fh, FILE%time_bounds_var_id,                               &
    "long_name", "Time bounds for each time stamp"                            &
  )
  CALL file_def_attr(                                                         &
    FILE%open_file%fh, FILE%time_bounds_var_id, "units",                      &
    "seconds since " // datetime_to_string(FILE%data_start)                   &
  )

  !-----------------------------------------------------------------------------
  ! Create the time index
  !-----------------------------------------------------------------------------
  ! We pass an array with 0 elements as the non-record dimensions, but indicate
  ! we want to use the record dim
  ALLOCATE(zero_element_array(0))
  FILE%time_index_var_id = file_def_var(                                      &
    FILE%open_file%fh, time_index_var_name, zero_element_array, .TRUE.        &
  )
  DEALLOCATE(zero_element_array)

  ! Define it's name and units
  CALL file_def_attr(                                                         &
    FILE%open_file%fh, FILE%time_index_var_id, "standard_name", "time"        &
  )
  CALL file_def_attr(                                                         &
    FILE%open_file%fh, FILE%time_index_var_id, "long_name", "Time of data"    &
  )
  CALL file_def_attr(                                                         &
    FILE%open_file%fh, FILE%time_index_var_id, "units",                       &
    "seconds since " // datetime_to_string(FILE%data_start)                   &
  )
  CALL file_def_attr(                                                         &
    FILE%open_file%fh, FILE%time_index_var_id, "bounds", "time_bounds"        &
  )

  ! Add the calendar attribute 
  IF ( l_360 ) THEN
    CALL file_def_attr(                                                       &
      FILE%open_file%fh, FILE%time_index_var_id, "calendar", "360_day"        &
    )
  ELSE IF ( l_leap ) THEN
    CALL file_def_attr(                                                       &
      FILE%open_file%fh, FILE%time_index_var_id, "calendar", "standard"       &
    )
  ELSE
    CALL file_def_attr(                                                       &
      FILE%open_file%fh, FILE%time_index_var_id, "calendar", "365_day"        &
    )
  END IF
END IF

!-----------------------------------------------------------------------------
! Define the variables and update the variable ids to the ones used by the
! newly opened file
!-----------------------------------------------------------------------------
DO i = 1,FILE%nvars
  ! Define the variable, indicating whether to use the record dimension
  FILE%vars(i)%current_id = file_gridded_def_var(                             &
    FILE%open_file, FILE%vars(i)%NAME,                                        &
  ! We have to convert the dimension ids to ids in the underlying file
        FILE%dims( FILE%vars(i)%levels_dims )%current_id,                     &
        FILE%vars(i)%use_time                                                 &
      )

  ! Define the attributes for the variable
  ! First the real valued attributes
  DO j = 1,FILE%vars(i)%attrs_real%length
    key = FILE%vars(i)%attrs_real%keys(j)
    CALL dict_get(FILE%vars(i)%attrs_real, key, real_val)

    CALL file_gridded_def_attr(                                               &
      FILE%open_file, FILE%vars(i)%current_id, key, real_val                  &
    )
  END DO

  ! Next, integer valued attributes
  DO j = 1,FILE%vars(i)%attrs_int%length
    key = FILE%vars(i)%attrs_int%keys(j)
    CALL dict_get(FILE%vars(i)%attrs_int, key, int_val)

    CALL file_gridded_def_attr(                                               &
      FILE%open_file, FILE%vars(i)%current_id, key, int_val                   &
    )
  END DO

  ! Lastly, character valued attributes
  DO j = 1,FILE%vars(i)%attrs_char%length
    key = FILE%vars(i)%attrs_char%keys(j)
    CALL dict_get(FILE%vars(i)%attrs_char, key, char_val)

    CALL file_gridded_def_attr(                                               &
      FILE%open_file, FILE%vars(i)%current_id, key, char_val                  &
    )
  END DO
END DO

!-----------------------------------------------------------------------------
! Lastly, define global attributes
!-----------------------------------------------------------------------------
! First the real valued attributes
DO i = 1,FILE%attrs_real%length
  key = FILE%attrs_real%keys(i)
  CALL dict_get(FILE%attrs_real, key, real_val)

  CALL file_gridded_def_attr(FILE%open_file, attr_global, key, real_val)
END DO

! Next, integer valued attributes
DO i = 1,FILE%attrs_int%length
  key = FILE%attrs_int%keys(i)
  CALL dict_get(FILE%attrs_int, key, int_val)

  CALL file_gridded_def_attr(FILE%open_file, attr_global, key, int_val)
END DO

! Lastly, character valued attributes
DO i = 1,FILE%attrs_char%length
  key = FILE%attrs_char%keys(i)
  CALL dict_get(FILE%attrs_char, key, char_val)

  CALL file_gridded_def_attr(FILE%open_file, attr_global, key, char_val)
END DO

!-----------------------------------------------------------------------------
! Take the file out of define mode, so that it is ready to read from/write to
!-----------------------------------------------------------------------------
CALL file_gridded_enddef(FILE%open_file)

!-----------------------------------------------------------------------------
! Now everything is defined on the new file, we write the data that we
! gathered for non-time-varying fields before we closed the old file
!-----------------------------------------------------------------------------
IF ( ALLOCATED(non_time_varying_data) ) THEN
  ! Write the data for the non-time varying variables
  j = 1
  DO i = 1,FILE%nvars
    IF ( .NOT. FILE%vars(i)%use_time ) THEN
      ! We are not using a subgrid, so give write_subgrid=F
      CALL file_ts_write_var(FILE, i, non_time_varying_data(j), .FALSE.)
      ! Advance the non-time variable counter
      j = j + 1
    END IF
  END DO

  DEALLOCATE(non_time_varying_data)
END IF

RETURN

END SUBROUTINE file_ts_internal_open_file

END MODULE file_ts_mod
