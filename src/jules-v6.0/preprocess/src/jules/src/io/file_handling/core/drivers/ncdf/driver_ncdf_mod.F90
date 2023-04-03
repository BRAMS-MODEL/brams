! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE driver_ncdf_mod

USE io_constants, ONLY: max_file_name_len, max_var_file

USE string_utils_mod, ONLY: to_string

! Import this globally as it is universally used
USE netcdf, ONLY: nf90_noerr

USE logging_mod, ONLY: log_info, log_fatal

IMPLICIT NONE


!-----------------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------------
CHARACTER(LEN=15) :: extensions(2)  ! The file name extensions recognised
DATA extensions / 'nc', 'cdf' /     ! for NetCDF files

!-----------------------------------------------------------------------------
! Type definitions
!-----------------------------------------------------------------------------
TYPE var_ncdf

  INTEGER :: id  ! The NetCDF id of the variable
  INTEGER :: ndims  ! The number of NON-RECORD dimensions that the
                     ! variable has
  LOGICAL :: is_record  ! Does the variable use the record dimension?

END TYPE var_ncdf


TYPE file_ncdf
  !-----------------------------------------------------------------------------
  ! This type contains information required to process a NetCDF file
  !-----------------------------------------------------------------------------
  INTEGER :: id  ! NetCDF id of file

  CHARACTER(LEN=max_file_name_len) :: NAME  ! The name of the file
                                            ! This exists purely for
                                            ! debugging purposes

  INTEGER :: mode ! The mode the file is open in
                  ! This can be one of mode_read or mode_write

  LOGICAL :: parallel = .FALSE.  ! Indicates if this file is being accessed
                                 ! in parallel
                                 ! When 1 is set, this will always
                                 ! be false, otherwise it will be true
                                 ! if the file was opened with comm and info

  INTEGER :: record_dim = -1  ! The id of the record dimension
                              ! -1 means no record dimension

  INTEGER :: current_record = 1  ! By default, we read from the 1st record

  INTEGER :: nvars = 0  ! The number of variables in the file
  TYPE(var_ncdf) :: vars(max_var_file)  ! The variables in the file

END TYPE file_ncdf


!-----------------------------------------------------------------------------
! Interface definitions for the overloads of read/write var for variables
! of different ranks
!-----------------------------------------------------------------------------
INTERFACE file_ncdf_read_var
MODULE PROCEDURE file_ncdf_read_var_scalar, file_ncdf_read_var_1d,            &
                 file_ncdf_read_var_2d,     file_ncdf_read_var_3d,            &
                 file_ncdf_read_var_4d,     file_ncdf_read_var_5d,            &
                 file_ncdf_read_var_6d,     file_ncdf_read_var_7d
END INTERFACE file_ncdf_read_var

INTERFACE file_ncdf_write_var
MODULE PROCEDURE file_ncdf_write_var_scalar, file_ncdf_write_var_1d,          &
                 file_ncdf_write_var_2d,     file_ncdf_write_var_3d,          &
                 file_ncdf_write_var_4d,     file_ncdf_write_var_5d,          &
                 file_ncdf_write_var_6d,     file_ncdf_write_var_7d
END INTERFACE file_ncdf_write_var


!-----------------------------------------------------------------------------
! Visibility declarations
!-----------------------------------------------------------------------------
PRIVATE
PUBLIC extensions,                                                            &
! Types
         file_ncdf,                                                           &
! Routines for opening and closing
         file_ncdf_open, file_ncdf_close,                                     &
! Routines for definitions
         file_ncdf_def_dim, file_ncdf_def_record_dim, file_ncdf_def_var,      &
         file_ncdf_def_attr_real, file_ncdf_def_attr_int,                     &
         file_ncdf_def_attr_char, file_ncdf_enddef,                           &
! Routines for seeking
         file_ncdf_seek, file_ncdf_advance,                                   &
! Routines for reading and writing
         file_ncdf_read_var, file_ncdf_write_var, file_ncdf_sync,             &
! Routines for inquiry
         file_ncdf_introspect, file_ncdf_inquire_dim, file_ncdf_inquire_var


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


FUNCTION file_ncdf_open(NAME, mode, comm, info) RESULT(FILE)

USE io_constants, ONLY: mode_read, mode_write

USE netcdf, ONLY:                                                             &
! Constants
    nf90_nowrite, nf90_clobber,                                               &
! Procedures
    nf90_open, nf90_create

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Opens a NetCDF file and returns a file_ncdf object representing it
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
CHARACTER(LEN=*), INTENT(IN) :: NAME
                            ! The name of the file
INTEGER, INTENT(IN) :: mode ! The mode to open the file
                            ! One of mode_read or mode_write
INTEGER, OPTIONAL, INTENT(IN) :: comm  ! MPI communicator to use for parallel
                                       ! I/O
                                       ! If not given, serial I/O is used
INTEGER, OPTIONAL, INTENT(IN) :: info  ! MPI info object to use for parallel
                                       ! I/O
                                       ! If not given, serial I/O is used

!-----------------------------------------------------------------------------
! We know from the conditions imposed by file_open that comm and info are
! either both present or both not present
!-----------------------------------------------------------------------------

! Return type
TYPE(file_ncdf) :: FILE


! Work variables
INTEGER :: ncid ! The NetCDF id of the opened file
INTEGER :: error ! Error code for any errors that occur

!-----------------------------------------------------------------------------


SELECT CASE ( mode )
CASE ( mode_read )
  CALL log_info("file_ncdf_open",                                             &
                "Opening file " // TRIM(NAME) // " for reading")
  ! Open file for reading only - file must already exist
  ! We don't need to specify comm or info to get parallel access here - any
  ! number of parallel readers are possible anyway...
  error = nf90_open(NAME, nf90_nowrite, ncid)

CASE ( mode_write )
  CALL log_info("file_ncdf_open",                                             &
                "Opening file " // TRIM(NAME) // " for writing")
  ! Create an empty file for (reading and) writing - if a file with the
  ! same name exists, overwrite it
  IF ( PRESENT(comm) ) THEN
    ! Don't try to specify MPI variables to NetCDF routines if we are using the
    ! dummy MPI library
    error = nf90_create(NAME, nf90_clobber, ncid)
  ELSE
    error = nf90_create(NAME, nf90_clobber, ncid)
  END IF

CASE DEFAULT
  ! Read and write are the only supported modes
  CALL log_fatal("file_ncdf_open",                                            &
                 "Unsupported mode - " // TRIM(to_string(mode)))

END SELECT

IF ( error /= nf90_noerr )                                                    &
  CALL log_fatal_ncdf("file_ncdf_open",                                       &
                      "Error opening file " // TRIM(NAME), error)


! Initialise the file_ncdf object
FILE%NAME = NAME
FILE%mode = mode
FILE%id   = ncid

RETURN

END FUNCTION file_ncdf_open
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_ncdf_def_dim(FILE, dim_name, dim_len) RESULT(dim_id)

USE io_constants, ONLY: max_sdf_name_len, mode_read, mode_write

USE netcdf, ONLY: nf90_inq_dimid, nf90_inquire_dimension, nf90_def_dim

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
TYPE(file_ncdf), INTENT(INOUT) :: FILE
                                ! The file to define the dimension on
CHARACTER(LEN=*), INTENT(IN) :: dim_name
                                ! The name of the dimension
INTEGER, INTENT(IN) :: dim_len  ! The length of the dimension

! Return type
INTEGER :: dim_id               ! The dimension id

! Work variables
INTEGER :: error  ! The current error code (if any)
INTEGER :: len_in_file  ! The length of the dimension in the file (only
                        ! used in read mode to validate dimension size)


!-----------------------------------------------------------------------------


! Check the dimension has a sensible size (i.e. > 0)
IF ( dim_len < 1 )                                                            &
  CALL log_fatal("file_ncdf_def_dim",                                         &
                 TRIM(FILE%NAME) // ": " //                                   &
                 "Dimension size should be >= 1 - size of " //                &
                 TRIM(to_string(dim_len)) // " given for dimension " //       &
                 TRIM(dim_name))

! We do different things depending on whether we are open for reading
! or writing
SELECT CASE ( FILE%mode )
CASE ( mode_read )
  ! In read mode, we just get the id of the dimension and validate its size
  error = nf90_inq_dimid(FILE%id, dim_name, dim_id)
  IF ( error /= nf90_noerr )                                                  &
    CALL log_fatal_ncdf("file_ncdf_def_dim",                                  &
                        TRIM(FILE%NAME) // ": " //                            &
                        "Error getting dimension id for dimension: " //       &
                        TRIM(dim_name), error)
      
  error = nf90_inquire_dimension(FILE%id, dim_id, LEN = len_in_file)
  IF ( error /= nf90_noerr )                                                  &
    CALL log_fatal_ncdf("file_ncdf_def_dim",                                  &
                        TRIM(FILE%NAME) // ": " //                            &
                        "dim_name =" // TRIM(dim_name) // ": " //             &
                        "Error getting dimension size", error)
  IF ( dim_len /= len_in_file )                                               &
    CALL log_fatal("file_ncdf_def_dim",                                       &
                   TRIM(FILE%NAME) // ": " //                                 &
                   "dim_name =" // TRIM(dim_name) // ": " //                  &
                   "Dimension length mismatch - length in file: " //          &
                   TRIM(to_string(len_in_file)) // ", expected: " //          &
                   TRIM(to_string(dim_len)))

CASE ( mode_write )
  ! In write mode, we define the dimension
  error = nf90_def_dim(FILE%id, dim_name, dim_len, dim_id)
  IF ( error /= nf90_noerr )                                                  &
    CALL log_fatal_ncdf("file_ncdf_def_dim",                                  &
                        TRIM(FILE%NAME) // ": " //                            &
                        "dim_name =" // TRIM(dim_name) // ": " //             &
                        "Error defining dimension", error)

  ! No default case as we already know that mode_read and mode_write are
  ! the only options
END SELECT

RETURN

END FUNCTION file_ncdf_def_dim
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_ncdf_def_record_dim(FILE, dim_name) RESULT(dim_id)

USE io_constants, ONLY: max_sdf_name_len, mode_read, mode_write

USE netcdf, ONLY: nf90_inq_dimid, nf90_def_dim, nf90_unlimited

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Defines the record dimension on the given file, returning the dimension id
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ncdf), INTENT(INOUT) :: FILE
                                ! The file to define the dimension on
CHARACTER(LEN=*), INTENT(IN) :: dim_name
                                ! The name of the dimension

! Return type
INTEGER :: dim_id               ! The dimension id

! Work variables
INTEGER :: error  ! The current error code (if any)


!-----------------------------------------------------------------------------

! We do different things depending on whether we are open for reading
! or writing
SELECT CASE ( FILE%mode )
CASE ( mode_read )
  ! In read mode, we just get the id of the named dimension
  ! We don't actually care whether it is a 'proper' unlimited dimension or not,
  ! since we won't be writing to it
  error = nf90_inq_dimid(FILE%id, dim_name, dim_id)
  IF ( error /= nf90_noerr )                                                  &
    CALL log_fatal_ncdf("file_ncdf_def_record_dim",                           &
                        TRIM(FILE%NAME) // ": " //                            &
                        "dim_name =" // TRIM(dim_name) // ": " //             &
                        "Error getting dimension id", error)

CASE ( mode_write )
  ! In write mode, we define an unlimited dimension
  error = nf90_def_dim(FILE%id, dim_name, nf90_unlimited, dim_id)
  IF ( error /= nf90_noerr )                                                  &
    CALL log_fatal_ncdf("file_ncdf_def_record_dim",                           &
                        TRIM(FILE%NAME) // ": " //                            &
                        "dim_name =" // TRIM(dim_name) // ": " //             &
                        "Error defining unlimited dimension", error)

  ! No default case as we already know that mode_read and mode_write are
  ! the only options
END SELECT

! Store the id of the record dimension for future use
FILE%record_dim = dim_id

RETURN

END FUNCTION file_ncdf_def_record_dim
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_ncdf_def_var(FILE, var_name, dims_in, is_record) RESULT(var_id)

USE io_constants, ONLY: max_sdf_name_len, max_dim_var, mode_read, mode_write

USE netcdf, ONLY: nf90_inq_varid, nf90_inquire_variable, nf90_def_var,        &
                   nf90_float
                     

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
TYPE(file_ncdf), INTENT(INOUT) :: FILE
                                  ! The file to define the variable in
CHARACTER(LEN=*), INTENT(IN) :: var_name
                                  ! The name of the variable
INTEGER, INTENT(IN) :: dims_in(:) ! The ids of the NON-RECORD dimensions of
                                  ! the variable
LOGICAL, INTENT(IN) :: is_record  ! Indicates whether the variable uses the
                                  ! record dimension

! Return type
INTEGER :: var_id                 ! The variable id - this is the index
                                  ! in the vars array of the file_ncdf object
                                  ! of the variable


! Work variables
INTEGER :: var_ncid  ! The NetCDF id of the variable

INTEGER :: ndims_in_file  ! The number of dimensions that the variable
                           ! has in the NetCDF file
                           ! Only used when in read mode for verification
                           ! purposes

INTEGER :: ndims     ! The number of dimensions supplied INCLUDING RECORD DIM

INTEGER, ALLOCATABLE :: dims_in_file(:) ! The dimension ids of the
                                        ! variable in the file

INTEGER :: n          ! Loop counters
INTEGER :: error      ! Error indicator


!-----------------------------------------------------------------------------
! Check that a variable name has been provided.
IF ( LEN_TRIM(var_name) < 1 )                                                 &
  CALL log_fatal( "file_ncdf_def_var",                                        &
                  TRIM(FILE%NAME) // ": No variable name provided." )

! The record dimension needs to be defined before any variables can use it
IF ( is_record .AND. FILE%record_dim < 0 )                                    &
  CALL log_fatal("file_ncdf_def_var",                                         &
                 TRIM(FILE%NAME) // ": " //                                   &
                 "var_name =" // TRIM(var_name) // ": " //                    &
                 "Record dimension must be defined before defining a " //     &
                 "variable as using it")

! Check that none of the given dimensions is the record dimension
! (note that we don't have to check first that a record dim has been defined,
! since file%record_dim takes the value -1 if not defined and no valid dimension
! should take this value anyway)
IF ( ANY(dims_in == FILE%record_dim) )                                        &
  CALL log_fatal("file_ncdf_def_var",                                         &
                 TRIM(FILE%NAME) // ": " //                                   &
                 "var_name =" // TRIM(var_name) // ": " //                    &
                 "Record dimension is in dimension list - use of the " //     &
                 "record dimension should be specified using 'is_record'")

! Work out how many dimensions the variable has in total
ndims = SIZE(dims_in)
IF ( ndims > max_dim_var )                                                    &
  CALL log_fatal("file_ncdf_def_var",                                         &
                 TRIM(FILE%NAME) // ": " //                                   &
                 "var_name =" // TRIM(var_name) // ": " //                    &
                 "Variable has too many dimensions - code only exists " //    &
                 "for variables with up to " //                               &
                 TRIM(to_string(max_dim_var)) // " dimensions")

IF ( is_record ) ndims = ndims + 1


SELECT CASE ( FILE%mode )
CASE ( mode_read )
  !-----------------------------------------------------------------------------
  ! In read mode, we just get the id of the variable and validate that it has
  ! the correct dimensions
  !-----------------------------------------------------------------------------

  ! Retrieve the variable id by its name
  error = nf90_inq_varid(FILE%id, var_name, var_ncid)
  IF ( error /= nf90_noerr )                                                  &
    CALL log_fatal_ncdf("file_ncdf_def_var",                                  &
                        TRIM(FILE%NAME) // ": " //                            &
                        "Error getting variable id for variable: " //         &
                        TRIM(var_name), error)

  ! Get the number of dimensions that the variable has in the file
  error = nf90_inquire_variable(FILE%id, var_ncid, ndims = ndims_in_file)
  IF ( error /= nf90_noerr )                                                  &
    CALL log_fatal_ncdf("file_ncdf_def_var",                                  &
                        TRIM(FILE%NAME) // ": " //                            &
                        "var_name =" // TRIM(var_name) // ": " //             &
                        "Error getting number of dimensions for variable",    &
                        error)

  ! Check that the number of dimensions match
  IF ( ndims_in_file /= ndims )                                               &
    CALL log_fatal("file_ncdf_def_var",                                       &
                   TRIM(FILE%NAME) // ": " //                                 &
                   "var_name =" // TRIM(var_name) // ": " //                  &
                   "Dimensions mismatch - number in file: " //                &
                   TRIM(to_string(ndims_in_file)) // ", expected: " //        &
                   TRIM(to_string(ndims)))

  ! Retrieve the dimension ids
  ALLOCATE(dims_in_file(ndims_in_file))
  error = nf90_inquire_variable(FILE%id, var_ncid, dimids = dims_in_file)
  IF ( error /= nf90_noerr )                                                  &
    CALL log_fatal_ncdf("file_ncdf_def_var",                                  &
                        TRIM(FILE%NAME) // ": " //                            &
                        "var_name =" // TRIM(var_name) // ": " //             &
                        "Error getting dimensions for variable", error)

  ! Check that the given dimension ids match
  DO n = 1,SIZE(dims_in)
    IF ( .NOT. ANY(dims_in_file == dims_in(n)) )                              &
      CALL log_fatal("file_ncdf_def_var",                                     &
                     TRIM(FILE%NAME) // ": " //                               &
                     "var_name =" // TRIM(var_name) // ": " //                &
                     "Given dimension ids do not match file")
  END DO

  ! Check that the variable is defined with the record dimension if it is
  ! supposed to be
  IF ( is_record .AND. ( .NOT. ANY(dims_in_file == FILE%record_dim) ) )       &
    CALL log_fatal("file_ncdf_def_var",                                       &
                   TRIM(FILE%NAME) // ": " //                                 &
                   "var_name =" // TRIM(var_name) // ": " //                  &
                   "Variable not defined with record dimension in file")

  ! Verification is over - the variable is defined as we were told it was
  DEALLOCATE(dims_in_file)

CASE ( mode_write )
  !-----------------------------------------------------------------------------
  ! In write mode, we need to define the variable
  !-----------------------------------------------------------------------------

  ! Construct a new dimensions array that has the record dimension appended
  ! if required
  ALLOCATE(dims_in_file(ndims))
  dims_in_file(1:SIZE(dims_in)) = dims_in(:)
  IF ( is_record ) dims_in_file(ndims) = FILE%record_dim

  ! Define the variable
  error = nf90_def_var(                                                       &
    FILE%id, var_name, nf90_float, dims_in_file, var_ncid                     &
  )
  IF ( error /= nf90_noerr )                                                  &
    CALL log_fatal_ncdf("file_ncdf_def_var",                                  &
                        TRIM(FILE%NAME) // ": " //                            &
                        "var_name =" // TRIM(var_name) // ": " //             &
                        "Error defining variable", error)

  DEALLOCATE(dims_in_file)
      

  ! No default case as we already know that mode_read and mode_write are
  ! the only options
END SELECT

!-----------------------------------------------------------------------------
! Set up the var_ncdf object corresponding to this variable and store it
!-----------------------------------------------------------------------------
IF ( FILE%nvars >= max_var_file )                                             &
  CALL log_fatal("file_ncdf_def_var",                                         &
                 TRIM(FILE%NAME) // ": " //                                   &
                 "Too many variables in file - try increasing max_var_file")

FILE%nvars = FILE%nvars + 1

! The return value is the index of the variable in the vars array on the
! file_ncdf object
var_id = FILE%nvars
FILE%vars(var_id)%id        = var_ncid
FILE%vars(var_id)%ndims     = SIZE(dims_in)
FILE%vars(var_id)%is_record = is_record

RETURN

END FUNCTION file_ncdf_def_var
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ncdf_def_attr_real(FILE, var_id, NAME, VALUE)

USE io_constants, ONLY: mode_read, mode_write, attr_global

USE netcdf, ONLY: nf90_put_att, nf90_global

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
TYPE(file_ncdf), INTENT(INOUT) :: FILE
                                  ! The file to define the attribute in
INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                  ! attribute on
CHARACTER(LEN=*), INTENT(IN) :: NAME
                                  ! The name of the attribute
REAL, INTENT(IN) :: VALUE         ! The value of the attribute


! Work variables
INTEGER :: var_id_local  ! Local version of var_id that we can modify
INTEGER :: error         ! Error indicator

!-----------------------------------------------------------------------------

SELECT CASE ( FILE%mode )
CASE ( mode_read )
  !-----------------------------------------------------------------------------
  ! In read mode, we just ignore the request, as we don't care about attributes
  !-----------------------------------------------------------------------------

CASE ( mode_write )
  !-----------------------------------------------------------------------------
  ! In write mode, we need to define the attribute
  !-----------------------------------------------------------------------------
  IF ( var_id == attr_global ) THEN
    ! If a global attribute has been requested, set var id to the NetCDF
    ! identifier for global attributes
    var_id_local = nf90_global
  ELSE
    ! Otherwise get the NetCDF id of the variable from the given index in vars
    ! array
    var_id_local = FILE%vars(var_id)%id
  END IF

  ! Try to actually create the attribute
  error = nf90_put_att(FILE%id, var_id_local, NAME, VALUE)
  IF ( error /= nf90_noerr )                                                  &
    CALL log_fatal_ncdf("file_ncdf_def_attr_real",                            &
                        "Error defining attribute", error)

  ! No default case as we already know that mode_read and mode_write are
  ! the only options
END SELECT

RETURN

END SUBROUTINE file_ncdf_def_attr_real


SUBROUTINE file_ncdf_def_attr_int(FILE, var_id, NAME, VALUE)

USE io_constants, ONLY: mode_read, mode_write, attr_global

USE netcdf, ONLY: nf90_put_att, nf90_global

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
TYPE(file_ncdf), INTENT(INOUT) :: FILE
                                  ! The file to define the attribute in
INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                  ! attribute on
CHARACTER(LEN=*), INTENT(IN) :: NAME
                                  ! The name of the attribute
INTEGER, INTENT(IN) :: VALUE      ! The value of the attribute


! Work variables
INTEGER :: var_id_local  ! Local version of var_id that we can modify
INTEGER :: error         ! Error indicator

!-----------------------------------------------------------------------------

SELECT CASE ( FILE%mode )
CASE ( mode_read )
  !-----------------------------------------------------------------------------
  ! In read mode, we just ignore the request, as we don't care about attributes
  !-----------------------------------------------------------------------------

CASE ( mode_write )
  !-----------------------------------------------------------------------------
  ! In write mode, we need to define the attribute
  !-----------------------------------------------------------------------------
  IF ( var_id == attr_global ) THEN
    ! If a global attribute has been requested, set var id to the NetCDF
    ! identifier for global attributes
    var_id_local = nf90_global
  ELSE
    ! Otherwise get the NetCDF id of the variable from the given index in vars
    ! array
    var_id_local = FILE%vars(var_id)%id
  END IF

  ! Try to actually create the attribute
  error = nf90_put_att(FILE%id, var_id_local, NAME, VALUE)
  IF ( error /= nf90_noerr )                                                  &
    CALL log_fatal_ncdf("file_ncdf_def_attr_real",                            &
                        "Error defining attribute", error)

  ! No default case as we already know that mode_read and mode_write are
  ! the only options
END SELECT

RETURN

END SUBROUTINE file_ncdf_def_attr_int


SUBROUTINE file_ncdf_def_attr_char(FILE, var_id, NAME, VALUE)

USE io_constants, ONLY: mode_read, mode_write, attr_global

USE netcdf, ONLY: nf90_put_att, nf90_global

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
TYPE(file_ncdf), INTENT(INOUT) :: FILE
                                  ! The file to define the attribute in
INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                  ! attribute on
CHARACTER(LEN=*), INTENT(IN) :: NAME
                                  ! The name of the attribute
CHARACTER(LEN=*), INTENT(IN) :: VALUE
                                  ! The value of the attribute


! Work variables
INTEGER :: var_id_local  ! Local version of var_id that we can modify
INTEGER :: error         ! Error indicator

!-----------------------------------------------------------------------------

SELECT CASE ( FILE%mode )
CASE ( mode_read )
  !-----------------------------------------------------------------------------
  ! In read mode, we just ignore the request, as we don't care about attributes
  !-----------------------------------------------------------------------------

CASE ( mode_write )
  !-----------------------------------------------------------------------------
  ! In write mode, we need to define the attribute
  !-----------------------------------------------------------------------------
  IF ( var_id == attr_global ) THEN
    ! If a global attribute has been requested, set var id to the NetCDF
    ! identifier for global attributes
    var_id_local = nf90_global
  ELSE
    ! Otherwise get the NetCDF id of the variable from the given index in vars
    ! array
    var_id_local = FILE%vars(var_id)%id
  END IF

  ! Try to actually create the attribute
  error = nf90_put_att(FILE%id, var_id_local, NAME, VALUE)
  IF ( error /= nf90_noerr )                                                  &
    CALL log_fatal_ncdf("file_ncdf_def_attr_real",                            &
                        "Error defining attribute", error)

  ! No default case as we already know that mode_read and mode_write are
  ! the only options
END SELECT

RETURN

END SUBROUTINE file_ncdf_def_attr_char
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ncdf_enddef(FILE)

USE io_constants, ONLY: mode_write

USE netcdf, ONLY: nf90_enddef

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
TYPE(file_ncdf), INTENT(INOUT) :: FILE

INTEGER :: error      ! Error indicator


!-----------------------------------------------------------------------------

! We are only genuinely in define mode if the file is open in write mode
IF ( FILE%mode == mode_write ) THEN
  error = nf90_enddef(FILE%id)
  IF ( error /= nf90_noerr )                                                  &
    CALL log_fatal_ncdf("file_ncdf_enddef",                                   &
                        "Error exiting define mode", error)
END IF

RETURN

END SUBROUTINE file_ncdf_enddef
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ncdf_seek(FILE, record)

USE io_constants, ONLY: mode_read

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Seeks the file to a particular record (i.e. the next time values are read
!   from the file using file_read_var, they will be read from the requested
!   record)
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ncdf), INTENT(INOUT) :: FILE  ! The file to seek
INTEGER, INTENT(IN) :: record           ! The record number to seek to


!-----------------------------------------------------------------------------

IF ( FILE%mode /= mode_read )                                                 &
  CALL log_fatal("file_ncdf_seek",                                            &
                 TRIM(FILE%NAME) // ": " //                                   &
                 "Arbitrary seeking is only allowed in read mode")

IF ( FILE%record_dim < 0 )                                                    &
  CALL log_fatal("file_ncdf_seek",                                            &
                 TRIM(FILE%NAME) // ": " //                                   &
                 "Cannot seek - no record dimension defined")

IF ( record < 1 )                                                             &
  CALL log_fatal("file_ncdf_seek",                                            &
                 TRIM(FILE%NAME) // ": " //                                   &
                 "Record number must be > 0")

FILE%current_record = record

RETURN

END SUBROUTINE file_ncdf_seek
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ncdf_advance(FILE)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Advances the file by one record
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ncdf), INTENT(INOUT) :: FILE  ! The file to advance


!-----------------------------------------------------------------------------

IF ( FILE%record_dim < 0 )                                                    &
  CALL log_fatal("file_ncdf_advance",                                         &
                 TRIM(FILE%NAME) // ": " //                                   &
                 "Cannot advance - no record dimension defined")

FILE%current_record = FILE%current_record + 1

RETURN

END SUBROUTINE file_ncdf_advance
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ncdf_read_var_scalar(FILE, var_id, VALUE, start)

USE io_constants, ONLY: max_dim_var, max_sdf_name_len

USE netcdf, ONLY: nf90_get_var, nf90_inquire_variable

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads a scalar value from a variable in the given file
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ncdf), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, INTENT(OUT) :: VALUE  ! The value read from file
INTEGER, INTENT(IN) :: start(max_dim_var)
                            ! The point to start reading from
                            ! This should contain one value for each
                            ! NON-RECORD dimension of the variable, and
                            ! any unused slots should be one

! Work variables
TYPE(var_ncdf) :: var  ! Structure containing information about the
                       ! variable

INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                       ! has, including record dimension

INTEGER, ALLOCATABLE :: local_start(:)

INTEGER :: error, error_name  ! Error indicators

CHARACTER(LEN=max_sdf_name_len) :: var_name  ! Name of variable.

!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
var = FILE%vars(var_id)

ndims_tot = var%ndims
IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we read from the correct record, if required
ALLOCATE(local_start(ndims_tot))
local_start(1:var%ndims) = start(1:var%ndims)
IF ( var%is_record ) local_start(ndims_tot) = FILE%current_record

! Try to read the variable
error = nf90_get_var(FILE%id, var%id, VALUE, local_start)
IF ( error /= nf90_noerr ) THEN
  !   Get variable name.
  error_name = nf90_inquire_variable(FILE%id, var%id, NAME = var_name )
  IF ( error_name /= nf90_noerr ) var_name = 'unknown' 
  CALL log_fatal_ncdf("file_ncdf_read_var_scalar",                            &
                      TRIM(FILE%NAME) // ": " //                              &
                      "Error reading variable '" // TRIM(var_name) // "'",    &
                       error)
END IF

DEALLOCATE(local_start)

RETURN

END SUBROUTINE file_ncdf_read_var_scalar

!#############################################################################

SUBROUTINE file_ncdf_read_var_1d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var, max_sdf_name_len

USE netcdf, ONLY: nf90_get_var, nf90_inquire_variable

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads a 1d array from a variable in the given file
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ncdf), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, DIMENSION(:), INTENT(OUT) :: values  ! The values read from file
INTEGER, INTENT(IN) :: start(max_dim_var)  ! The point to start reading from
                                 ! (one value for each NON-RECORD
                                 ! dimension of the variable)
INTEGER, INTENT(IN) :: counter(max_dim_var)  ! The number of points to read
                                 ! in each dimension of the variable


! Work variables
TYPE(var_ncdf) :: var  ! Structure containing information about the
                       ! variable

INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                       ! has, including record dimension

INTEGER, ALLOCATABLE :: local_start(:)
INTEGER, ALLOCATABLE :: local_count(:)

INTEGER :: error, error_name  ! Error indicators

CHARACTER(LEN=max_sdf_name_len) :: var_name  ! Name of variable.


!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
var = FILE%vars(var_id)

ndims_tot = var%ndims
IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we read from the correct record, if required
ALLOCATE(local_start(ndims_tot))
local_start(1:var%ndims) = start(1:var%ndims)
IF ( var%is_record ) local_start(ndims_tot) = FILE%current_record

! Set up the local_count array - this consists of one value for each
! dimension of the variable from the given counter array, plus a value so that
! we read from only the correct record, if required
ALLOCATE(local_count(ndims_tot))
local_count(1:var%ndims) = counter(1:var%ndims)
IF ( var%is_record ) local_count(ndims_tot) = 1

! Try to read the variable
error = nf90_get_var(FILE%id, var%id, values, local_start, local_count)
IF ( error /= nf90_noerr ) THEN
  !   Get variable name.
  error_name = nf90_inquire_variable(FILE%id, var%id, NAME = var_name )
  IF ( error_name /= nf90_noerr ) var_name = 'unknown' 
  CALL log_fatal_ncdf("file_ncdf_read_var_1d",                                &
                      TRIM(FILE%NAME) // ": " //                              &
                      "Error reading variable '" // TRIM(var_name) // "'",    &
                       error)
END IF

DEALLOCATE(local_start)
DEALLOCATE(local_count)

RETURN

END SUBROUTINE file_ncdf_read_var_1d

!#############################################################################

SUBROUTINE file_ncdf_read_var_2d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var, max_sdf_name_len

USE netcdf, ONLY: nf90_get_var, nf90_inquire_variable

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads a 2d array from a variable in the given file
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ncdf), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, DIMENSION(:,:), INTENT(OUT) :: values  ! The values read from file
INTEGER, INTENT(IN) :: start(max_dim_var)  ! The point to start reading from
                                 ! (one value for each NON-RECORD
                                 ! dimension of the variable)
INTEGER, INTENT(IN) :: counter(max_dim_var)  ! The number of points to read
                                 ! in each dimension of the variable


! Work variables
TYPE(var_ncdf) :: var  ! Structure containing information about the
                       ! variable

INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                       ! has, including record dimension

INTEGER, ALLOCATABLE :: local_start(:)
INTEGER, ALLOCATABLE :: local_count(:)

INTEGER :: error, error_name  ! Error indicators

CHARACTER(LEN=max_sdf_name_len) :: var_name  ! Name of variable.


!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
var = FILE%vars(var_id)

ndims_tot = var%ndims
IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we read from the correct record, if required
ALLOCATE(local_start(ndims_tot))
local_start(1:var%ndims) = start(1:var%ndims)
IF ( var%is_record ) local_start(ndims_tot) = FILE%current_record

! Set up the local_count array - this consists of one value for each
! dimension of the variable from the given counter array, plus a value so that
! we read from only the correct record, if required
ALLOCATE(local_count(ndims_tot))
local_count(1:var%ndims) = counter(1:var%ndims)
IF ( var%is_record ) local_count(ndims_tot) = 1

! Try to read the variable
error = nf90_get_var(FILE%id, var%id, values, local_start, local_count)
IF ( error /= nf90_noerr ) THEN
  !   Get variable name.
  error_name = nf90_inquire_variable(FILE%id, var%id, NAME = var_name )
  IF ( error_name /= nf90_noerr ) var_name = 'unknown' 
  CALL log_fatal_ncdf("file_ncdf_read_var_2d",                                &
                      TRIM(FILE%NAME) // ": " //                              &
                      "Error reading variable '" // TRIM(var_name) // "'",    &
                       error)
END IF

DEALLOCATE(local_start)
DEALLOCATE(local_count)

RETURN

END SUBROUTINE file_ncdf_read_var_2d

!#############################################################################

SUBROUTINE file_ncdf_read_var_3d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var, max_sdf_name_len

USE netcdf, ONLY: nf90_get_var, nf90_inquire_variable

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads a 3d array from a variable in the given file
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ncdf), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, DIMENSION(:,:,:), INTENT(OUT) :: values  ! The values read from file
INTEGER, INTENT(IN) :: start(max_dim_var)  ! The point to start reading from
                                 ! (one value for each NON-RECORD
                                 ! dimension of the variable)
INTEGER, INTENT(IN) :: counter(max_dim_var)  ! The number of points to read
                                 ! in each dimension of the variable


! Work variables
TYPE(var_ncdf) :: var  ! Structure containing information about the
                       ! variable

INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                       ! has, including record dimension

INTEGER, ALLOCATABLE :: local_start(:)
INTEGER, ALLOCATABLE :: local_count(:)

INTEGER :: error, error_name  ! Error indicators

CHARACTER(LEN=max_sdf_name_len) :: var_name  ! Name of variable.


!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
var = FILE%vars(var_id)

ndims_tot = var%ndims
IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we read from the correct record, if required
ALLOCATE(local_start(ndims_tot))
local_start(1:var%ndims) = start(1:var%ndims)
IF ( var%is_record ) local_start(ndims_tot) = FILE%current_record

! Set up the local_count array - this consists of one value for each
! dimension of the variable from the given counter array, plus a value so that
! we read from only the correct record, if required
ALLOCATE(local_count(ndims_tot))
local_count(1:var%ndims) = counter(1:var%ndims)
IF ( var%is_record ) local_count(ndims_tot) = 1

! Try to read the variable
error = nf90_get_var(FILE%id, var%id, values, local_start, local_count)
IF ( error /= nf90_noerr ) THEN
  !   Get variable name.
  error_name = nf90_inquire_variable(FILE%id, var%id, NAME = var_name )
  IF ( error_name /= nf90_noerr ) var_name = 'unknown' 
  CALL log_fatal_ncdf("file_ncdf_read_var_3d",                                &
                      TRIM(FILE%NAME) // ": " //                              &
                      "Error reading variable '" // TRIM(var_name) // "'",    &
                       error)
END IF

DEALLOCATE(local_start)
DEALLOCATE(local_count)

RETURN

END SUBROUTINE file_ncdf_read_var_3d

!#############################################################################

SUBROUTINE file_ncdf_read_var_4d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var, max_sdf_name_len

USE netcdf, ONLY: nf90_get_var, nf90_inquire_variable

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads a 4d array from a variable in the given file
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ncdf), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: values  ! The values read from file
INTEGER, INTENT(IN) :: start(max_dim_var)  ! The point to start reading from
                                 ! (one value for each NON-RECORD
                                 ! dimension of the variable)
INTEGER, INTENT(IN) :: counter(max_dim_var)  ! The number of points to read
                                 ! in each dimension of the variable


! Work variables
TYPE(var_ncdf) :: var  ! Structure containing information about the
                       ! variable

INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                       ! has, including record dimension

INTEGER, ALLOCATABLE :: local_start(:)
INTEGER, ALLOCATABLE :: local_count(:)

INTEGER :: error, error_name  ! Error indicators

CHARACTER(LEN=max_sdf_name_len) :: var_name  ! Name of variable.


!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
var = FILE%vars(var_id)

ndims_tot = var%ndims
IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we read from the correct record, if required
ALLOCATE(local_start(ndims_tot))
local_start(1:var%ndims) = start(1:var%ndims)
IF ( var%is_record ) local_start(ndims_tot) = FILE%current_record

! Set up the local_count array - this consists of one value for each
! dimension of the variable from the given counter array, plus a value so that
! we read from only the correct record, if required
ALLOCATE(local_count(ndims_tot))
local_count(1:var%ndims) = counter(1:var%ndims)
IF ( var%is_record ) local_count(ndims_tot) = 1

! Try to read the variable
error = nf90_get_var(FILE%id, var%id, values, local_start, local_count)
IF ( error /= nf90_noerr ) THEN
  !   Get variable name.
  error_name = nf90_inquire_variable(FILE%id, var%id, NAME = var_name )
  IF ( error_name /= nf90_noerr ) var_name = 'unknown' 
  CALL log_fatal_ncdf("file_ncdf_read_var_4d",                                &
                      TRIM(FILE%NAME) // ": " //                              &
                      "Error reading variable '" // TRIM(var_name) // "'",    &
                       error)
END IF

DEALLOCATE(local_start)
DEALLOCATE(local_count)

RETURN

END SUBROUTINE file_ncdf_read_var_4d

!#############################################################################

SUBROUTINE file_ncdf_read_var_5d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var, max_sdf_name_len

USE netcdf, ONLY: nf90_get_var, nf90_inquire_variable

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads a 5d array from a variable in the given file
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ncdf), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: values  ! The values read from file
INTEGER, INTENT(IN) :: start(max_dim_var)  ! The point to start reading from
                                 ! (one value for each NON-RECORD
                                 ! dimension of the variable)
INTEGER, INTENT(IN) :: counter(max_dim_var)  ! The number of points to read
                                 ! in each dimension of the variable


! Work variables
TYPE(var_ncdf) :: var  ! Structure containing information about the
                       ! variable

INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                       ! has, including record dimension

INTEGER, ALLOCATABLE :: local_start(:)
INTEGER, ALLOCATABLE :: local_count(:)

INTEGER :: error, error_name  ! Error indicators

CHARACTER(LEN=max_sdf_name_len) :: var_name  ! Name of variable.


!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
var = FILE%vars(var_id)

ndims_tot = var%ndims
IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we read from the correct record, if required
ALLOCATE(local_start(ndims_tot))
local_start(1:var%ndims) = start(1:var%ndims)
IF ( var%is_record ) local_start(ndims_tot) = FILE%current_record

! Set up the local_count array - this consists of one value for each
! dimension of the variable from the given counter array, plus a value so that
! we read from only the correct record, if required
ALLOCATE(local_count(ndims_tot))
local_count(1:var%ndims) = counter(1:var%ndims)
IF ( var%is_record ) local_count(ndims_tot) = 1

! Try to read the variable
error = nf90_get_var(FILE%id, var%id, values, local_start, local_count)
IF ( error /= nf90_noerr ) THEN
  !   Get variable name.
  error_name = nf90_inquire_variable(FILE%id, var%id, NAME = var_name )
  IF ( error_name /= nf90_noerr ) var_name = 'unknown' 
  CALL log_fatal_ncdf("file_ncdf_read_var_5d",                                &
                      TRIM(FILE%NAME) // ": " //                              &
                      "Error reading variable '" // TRIM(var_name) // "'",    &
                       error)
END IF

DEALLOCATE(local_start)
DEALLOCATE(local_count)

RETURN

END SUBROUTINE file_ncdf_read_var_5d

!#############################################################################

SUBROUTINE file_ncdf_read_var_6d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var, max_sdf_name_len

USE netcdf, ONLY: nf90_get_var, nf90_inquire_variable

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads a 6d array from a variable in the given file
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ncdf), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, DIMENSION(:,:,:,:,:,:), INTENT(OUT) :: values  ! The values read from file
INTEGER, INTENT(IN) :: start(max_dim_var)  ! The point to start reading from
                                 ! (one value for each NON-RECORD
                                 ! dimension of the variable)
INTEGER, INTENT(IN) :: counter(max_dim_var)  ! The number of points to read
                                 ! in each dimension of the variable


! Work variables
TYPE(var_ncdf) :: var  ! Structure containing information about the
                       ! variable

INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                       ! has, including record dimension

INTEGER, ALLOCATABLE :: local_start(:)
INTEGER, ALLOCATABLE :: local_count(:)

INTEGER :: error, error_name  ! Error indicators

CHARACTER(LEN=max_sdf_name_len) :: var_name  ! Name of variable.


!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
var = FILE%vars(var_id)

ndims_tot = var%ndims
IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we read from the correct record, if required
ALLOCATE(local_start(ndims_tot))
local_start(1:var%ndims) = start(1:var%ndims)
IF ( var%is_record ) local_start(ndims_tot) = FILE%current_record

! Set up the local_count array - this consists of one value for each
! dimension of the variable from the given counter array, plus a value so that
! we read from only the correct record, if required
ALLOCATE(local_count(ndims_tot))
local_count(1:var%ndims) = counter(1:var%ndims)
IF ( var%is_record ) local_count(ndims_tot) = 1

! Try to read the variable
error = nf90_get_var(FILE%id, var%id, values, local_start, local_count)
IF ( error /= nf90_noerr ) THEN
  !   Get variable name.
  error_name = nf90_inquire_variable(FILE%id, var%id, NAME = var_name )
  IF ( error_name /= nf90_noerr ) var_name = 'unknown' 
  CALL log_fatal_ncdf("file_ncdf_read_var_6d",                                &
                      TRIM(FILE%NAME) // ": " //                              &
                      "Error reading variable '" // TRIM(var_name) // "'",    &
                       error)
END IF

DEALLOCATE(local_start)
DEALLOCATE(local_count)

RETURN

END SUBROUTINE file_ncdf_read_var_6d

!#############################################################################

SUBROUTINE file_ncdf_read_var_7d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var, max_sdf_name_len

USE netcdf, ONLY: nf90_get_var, nf90_inquire_variable

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads a 7d array from a variable in the given file
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ncdf), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, DIMENSION(:,:,:,:,:,:,:), INTENT(OUT) :: values  ! The values read from file
INTEGER, INTENT(IN) :: start(max_dim_var)  ! The point to start reading from
                                 ! (one value for each NON-RECORD
                                 ! dimension of the variable)
INTEGER, INTENT(IN) :: counter(max_dim_var)  ! The number of points to read
                                 ! in each dimension of the variable


! Work variables
TYPE(var_ncdf) :: var  ! Structure containing information about the
                       ! variable

INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                       ! has, including record dimension

INTEGER, ALLOCATABLE :: local_start(:)
INTEGER, ALLOCATABLE :: local_count(:)

INTEGER :: error, error_name  ! Error indicators

CHARACTER(LEN=max_sdf_name_len) :: var_name  ! Name of variable.


!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
var = FILE%vars(var_id)

ndims_tot = var%ndims
IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we read from the correct record, if required
ALLOCATE(local_start(ndims_tot))
local_start(1:var%ndims) = start(1:var%ndims)
IF ( var%is_record ) local_start(ndims_tot) = FILE%current_record

! Set up the local_count array - this consists of one value for each
! dimension of the variable from the given counter array, plus a value so that
! we read from only the correct record, if required
ALLOCATE(local_count(ndims_tot))
local_count(1:var%ndims) = counter(1:var%ndims)
IF ( var%is_record ) local_count(ndims_tot) = 1

! Try to read the variable
error = nf90_get_var(FILE%id, var%id, values, local_start, local_count)
IF ( error /= nf90_noerr ) THEN
  !   Get variable name.
  error_name = nf90_inquire_variable(FILE%id, var%id, NAME = var_name )
  IF ( error_name /= nf90_noerr ) var_name = 'unknown' 
  CALL log_fatal_ncdf("file_ncdf_read_var_7d",                                &
                      TRIM(FILE%NAME) // ": " //                              &
                      "Error reading variable '" // TRIM(var_name) // "'",    &
                       error)
END IF

DEALLOCATE(local_start)
DEALLOCATE(local_count)

RETURN

END SUBROUTINE file_ncdf_read_var_7d
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ncdf_write_var_scalar(FILE, var_id, VALUE, start)

USE io_constants, ONLY: max_dim_var, max_sdf_name_len

USE netcdf, ONLY: nf90_inquire_variable, nf90_put_var

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Writes a scalar value to a variable in the given file
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ncdf), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, INTENT(IN) :: VALUE  ! The value to write
INTEGER, INTENT(IN) :: start(max_dim_var)
                            ! The point to start writing at
                            ! This should contain one value for each
                            ! NON-RECORD dimension of the variable, and
                            ! any unused slots should be one

! Work variables
TYPE(var_ncdf) :: var  ! Structure containing information about the
                       ! variable

INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                       ! has, including record dimension

INTEGER, ALLOCATABLE :: local_start(:)

INTEGER :: error, error_name  ! Error indicators

CHARACTER(LEN=max_sdf_name_len) :: var_name  ! Name of variable.

!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
var = FILE%vars(var_id)

ndims_tot = var%ndims
IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we write to the correct record, if required
ALLOCATE(local_start(ndims_tot))
local_start(1:var%ndims) = start(1:var%ndims)
IF ( var%is_record ) local_start(ndims_tot) = FILE%current_record

! Try to write the variable
error = nf90_put_var(FILE%id, var%id, VALUE, local_start)
IF ( error /= nf90_noerr ) THEN
  !   Get variable name.
  error_name = nf90_inquire_variable(FILE%id, var%id, NAME = var_name )
  IF ( error_name /= nf90_noerr ) var_name = 'unknown' 
  CALL log_fatal_ncdf("file_ncdf_write_var_scalar",                           &
                      TRIM(FILE%NAME) // ": " //                              &
                      "Error writing variable '" // TRIM(var_name) // "'",    &
                       error)
END IF

DEALLOCATE(local_start)

RETURN

END SUBROUTINE file_ncdf_write_var_scalar

!#############################################################################

SUBROUTINE file_ncdf_write_var_1d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var, max_sdf_name_len

USE netcdf, ONLY: nf90_inquire_variable, nf90_put_var

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Writes a 1d array to a variable in the given file
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ncdf), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, DIMENSION(:), INTENT(IN) :: values  ! The values to write
INTEGER, INTENT(IN) :: start(max_dim_var)  ! The point to start writing
                                 ! (one value for each NON-RECORD
                                 ! dimension of the variable)
INTEGER, INTENT(IN) :: counter(max_dim_var)  ! The number of points to write
                                 ! in each dimension of the variable


! Work variables
TYPE(var_ncdf) :: var  ! Structure containing information about the
                       ! variable

INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                       ! has, including record dimension

INTEGER, ALLOCATABLE :: local_start(:)
INTEGER, ALLOCATABLE :: local_count(:)

INTEGER :: error, error_name  ! Error indicators

CHARACTER(LEN=max_sdf_name_len) :: var_name  ! Name of variable.


!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
var = FILE%vars(var_id)

ndims_tot = var%ndims
IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we write to the correct record, if required
ALLOCATE(local_start(ndims_tot))
local_start(1:var%ndims) = start(1:var%ndims)
IF ( var%is_record ) local_start(ndims_tot) = FILE%current_record

! Set up the local_count array - this consists of one value for each
! dimension of the variable from the given counter array, plus a value so that
! we only write to the correct record, if required
ALLOCATE(local_count(ndims_tot))
local_count(1:var%ndims) = counter(1:var%ndims)
IF ( var%is_record ) local_count(ndims_tot) = 1

! Try to write the variable
error = nf90_put_var(FILE%id, var%id, values, local_start, local_count)
IF ( error /= nf90_noerr ) THEN
  !   Get variable name.
  error_name = nf90_inquire_variable(FILE%id, var%id, NAME = var_name )
  IF ( error_name /= nf90_noerr ) var_name = 'unknown' 
  CALL log_fatal_ncdf("file_ncdf_write_var_1d",                               &
                      TRIM(FILE%NAME) // ": " //                              &
                      "Error writing variable '" // TRIM(var_name) // "'",    &
                       error)
END IF

DEALLOCATE(local_start)
DEALLOCATE(local_count)

RETURN

END SUBROUTINE file_ncdf_write_var_1d

!#############################################################################

SUBROUTINE file_ncdf_write_var_2d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var, max_sdf_name_len

USE netcdf, ONLY: nf90_inquire_variable, nf90_put_var

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Writes a 2d array to a variable in the given file
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ncdf), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, DIMENSION(:,:), INTENT(IN) :: values  ! The values to write
INTEGER, INTENT(IN) :: start(max_dim_var)  ! The point to start writing
                                 ! (one value for each NON-RECORD
                                 ! dimension of the variable)
INTEGER, INTENT(IN) :: counter(max_dim_var)  ! The number of points to write
                                 ! in each dimension of the variable


! Work variables
TYPE(var_ncdf) :: var  ! Structure containing information about the
                       ! variable

INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                       ! has, including record dimension

INTEGER, ALLOCATABLE :: local_start(:)
INTEGER, ALLOCATABLE :: local_count(:)

INTEGER :: error, error_name  ! Error indicators

CHARACTER(LEN=max_sdf_name_len) :: var_name  ! Name of variable.


!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
var = FILE%vars(var_id)

ndims_tot = var%ndims
IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we write to the correct record, if required
ALLOCATE(local_start(ndims_tot))
local_start(1:var%ndims) = start(1:var%ndims)
IF ( var%is_record ) local_start(ndims_tot) = FILE%current_record

! Set up the local_count array - this consists of one value for each
! dimension of the variable from the given counter array, plus a value so that
! we only write to the correct record, if required
ALLOCATE(local_count(ndims_tot))
local_count(1:var%ndims) = counter(1:var%ndims)
IF ( var%is_record ) local_count(ndims_tot) = 1

! Try to write the variable
error = nf90_put_var(FILE%id, var%id, values, local_start, local_count)
IF ( error /= nf90_noerr ) THEN
  !   Get variable name.
  error_name = nf90_inquire_variable(FILE%id, var%id, NAME = var_name )
  IF ( error_name /= nf90_noerr ) var_name = 'unknown' 
  CALL log_fatal_ncdf("file_ncdf_write_var_2d",                               &
                      TRIM(FILE%NAME) // ": " //                              &
                      "Error writing variable '" // TRIM(var_name) // "'",    &
                       error)
END IF

DEALLOCATE(local_start)
DEALLOCATE(local_count)

RETURN

END SUBROUTINE file_ncdf_write_var_2d

!#############################################################################

SUBROUTINE file_ncdf_write_var_3d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var, max_sdf_name_len

USE netcdf, ONLY: nf90_inquire_variable, nf90_put_var

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Writes a 3d array to a variable in the given file
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ncdf), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, DIMENSION(:,:,:), INTENT(IN) :: values  ! The values to write
INTEGER, INTENT(IN) :: start(max_dim_var)  ! The point to start writing
                                 ! (one value for each NON-RECORD
                                 ! dimension of the variable)
INTEGER, INTENT(IN) :: counter(max_dim_var)  ! The number of points to write
                                 ! in each dimension of the variable


! Work variables
TYPE(var_ncdf) :: var  ! Structure containing information about the
                       ! variable

INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                       ! has, including record dimension

INTEGER, ALLOCATABLE :: local_start(:)
INTEGER, ALLOCATABLE :: local_count(:)

INTEGER :: error, error_name  ! Error indicators

CHARACTER(LEN=max_sdf_name_len) :: var_name  ! Name of variable.


!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
var = FILE%vars(var_id)

ndims_tot = var%ndims
IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we write to the correct record, if required
ALLOCATE(local_start(ndims_tot))
local_start(1:var%ndims) = start(1:var%ndims)
IF ( var%is_record ) local_start(ndims_tot) = FILE%current_record

! Set up the local_count array - this consists of one value for each
! dimension of the variable from the given counter array, plus a value so that
! we only write to the correct record, if required
ALLOCATE(local_count(ndims_tot))
local_count(1:var%ndims) = counter(1:var%ndims)
IF ( var%is_record ) local_count(ndims_tot) = 1

! Try to write the variable
error = nf90_put_var(FILE%id, var%id, values, local_start, local_count)
IF ( error /= nf90_noerr ) THEN
  !   Get variable name.
  error_name = nf90_inquire_variable(FILE%id, var%id, NAME = var_name )
  IF ( error_name /= nf90_noerr ) var_name = 'unknown' 
  CALL log_fatal_ncdf("file_ncdf_write_var_3d",                               &
                      TRIM(FILE%NAME) // ": " //                              &
                      "Error writing variable '" // TRIM(var_name) // "'",    &
                       error)
END IF

DEALLOCATE(local_start)
DEALLOCATE(local_count)

RETURN

END SUBROUTINE file_ncdf_write_var_3d

!#############################################################################

SUBROUTINE file_ncdf_write_var_4d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var, max_sdf_name_len

USE netcdf, ONLY: nf90_inquire_variable, nf90_put_var

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Writes a 4d array to a variable in the given file
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ncdf), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: values  ! The values to write
INTEGER, INTENT(IN) :: start(max_dim_var)  ! The point to start writing
                                 ! (one value for each NON-RECORD
                                 ! dimension of the variable)
INTEGER, INTENT(IN) :: counter(max_dim_var)  ! The number of points to write
                                 ! in each dimension of the variable


! Work variables
TYPE(var_ncdf) :: var  ! Structure containing information about the
                       ! variable

INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                       ! has, including record dimension

INTEGER, ALLOCATABLE :: local_start(:)
INTEGER, ALLOCATABLE :: local_count(:)

INTEGER :: error, error_name  ! Error indicators

CHARACTER(LEN=max_sdf_name_len) :: var_name  ! Name of variable.


!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
var = FILE%vars(var_id)

ndims_tot = var%ndims
IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we write to the correct record, if required
ALLOCATE(local_start(ndims_tot))
local_start(1:var%ndims) = start(1:var%ndims)
IF ( var%is_record ) local_start(ndims_tot) = FILE%current_record

! Set up the local_count array - this consists of one value for each
! dimension of the variable from the given counter array, plus a value so that
! we only write to the correct record, if required
ALLOCATE(local_count(ndims_tot))
local_count(1:var%ndims) = counter(1:var%ndims)
IF ( var%is_record ) local_count(ndims_tot) = 1

! Try to write the variable
error = nf90_put_var(FILE%id, var%id, values, local_start, local_count)
IF ( error /= nf90_noerr ) THEN
  !   Get variable name.
  error_name = nf90_inquire_variable(FILE%id, var%id, NAME = var_name )
  IF ( error_name /= nf90_noerr ) var_name = 'unknown' 
  CALL log_fatal_ncdf("file_ncdf_write_var_4d",                               &
                      TRIM(FILE%NAME) // ": " //                              &
                      "Error writing variable '" // TRIM(var_name) // "'",    &
                       error)
END IF

DEALLOCATE(local_start)
DEALLOCATE(local_count)

RETURN

END SUBROUTINE file_ncdf_write_var_4d

!#############################################################################

SUBROUTINE file_ncdf_write_var_5d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var, max_sdf_name_len

USE netcdf, ONLY: nf90_inquire_variable, nf90_put_var

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Writes a 5d array to a variable in the given file
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ncdf), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, DIMENSION(:,:,:,:,:), INTENT(IN) :: values  ! The values to write
INTEGER, INTENT(IN) :: start(max_dim_var)  ! The point to start writing
                                 ! (one value for each NON-RECORD
                                 ! dimension of the variable)
INTEGER, INTENT(IN) :: counter(max_dim_var)  ! The number of points to write
                                 ! in each dimension of the variable


! Work variables
TYPE(var_ncdf) :: var  ! Structure containing information about the
                       ! variable

INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                       ! has, including record dimension

INTEGER, ALLOCATABLE :: local_start(:)
INTEGER, ALLOCATABLE :: local_count(:)

INTEGER :: error, error_name  ! Error indicators

CHARACTER(LEN=max_sdf_name_len) :: var_name  ! Name of variable.


!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
var = FILE%vars(var_id)

ndims_tot = var%ndims
IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we write to the correct record, if required
ALLOCATE(local_start(ndims_tot))
local_start(1:var%ndims) = start(1:var%ndims)
IF ( var%is_record ) local_start(ndims_tot) = FILE%current_record

! Set up the local_count array - this consists of one value for each
! dimension of the variable from the given counter array, plus a value so that
! we only write to the correct record, if required
ALLOCATE(local_count(ndims_tot))
local_count(1:var%ndims) = counter(1:var%ndims)
IF ( var%is_record ) local_count(ndims_tot) = 1

! Try to write the variable
error = nf90_put_var(FILE%id, var%id, values, local_start, local_count)
IF ( error /= nf90_noerr ) THEN
  !   Get variable name.
  error_name = nf90_inquire_variable(FILE%id, var%id, NAME = var_name )
  IF ( error_name /= nf90_noerr ) var_name = 'unknown' 
  CALL log_fatal_ncdf("file_ncdf_write_var_5d",                               &
                      TRIM(FILE%NAME) // ": " //                              &
                      "Error writing variable '" // TRIM(var_name) // "'",    &
                       error)
END IF

DEALLOCATE(local_start)
DEALLOCATE(local_count)

RETURN

END SUBROUTINE file_ncdf_write_var_5d

!#############################################################################

SUBROUTINE file_ncdf_write_var_6d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var, max_sdf_name_len

USE netcdf, ONLY: nf90_inquire_variable, nf90_put_var

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Writes a 6d array to a variable in the given file
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ncdf), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, DIMENSION(:,:,:,:,:,:), INTENT(IN) :: values  ! The values to write
INTEGER, INTENT(IN) :: start(max_dim_var)  ! The point to start writing
                                 ! (one value for each NON-RECORD
                                 ! dimension of the variable)
INTEGER, INTENT(IN) :: counter(max_dim_var)  ! The number of points to write
                                 ! in each dimension of the variable


! Work variables
TYPE(var_ncdf) :: var  ! Structure containing information about the
                       ! variable

INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                       ! has, including record dimension

INTEGER, ALLOCATABLE :: local_start(:)
INTEGER, ALLOCATABLE :: local_count(:)

INTEGER :: error, error_name  ! Error indicators

CHARACTER(LEN=max_sdf_name_len) :: var_name  ! Name of variable.


!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
var = FILE%vars(var_id)

ndims_tot = var%ndims
IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we write to the correct record, if required
ALLOCATE(local_start(ndims_tot))
local_start(1:var%ndims) = start(1:var%ndims)
IF ( var%is_record ) local_start(ndims_tot) = FILE%current_record

! Set up the local_count array - this consists of one value for each
! dimension of the variable from the given counter array, plus a value so that
! we only write to the correct record, if required
ALLOCATE(local_count(ndims_tot))
local_count(1:var%ndims) = counter(1:var%ndims)
IF ( var%is_record ) local_count(ndims_tot) = 1

! Try to write the variable
error = nf90_put_var(FILE%id, var%id, values, local_start, local_count)
IF ( error /= nf90_noerr ) THEN
  !   Get variable name.
  error_name = nf90_inquire_variable(FILE%id, var%id, NAME = var_name )
  IF ( error_name /= nf90_noerr ) var_name = 'unknown' 
  CALL log_fatal_ncdf("file_ncdf_write_var_6d",                               &
                      TRIM(FILE%NAME) // ": " //                              &
                      "Error writing variable '" // TRIM(var_name) // "'",    &
                       error)
END IF

DEALLOCATE(local_start)
DEALLOCATE(local_count)

RETURN

END SUBROUTINE file_ncdf_write_var_6d

!#############################################################################

SUBROUTINE file_ncdf_write_var_7d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var, max_sdf_name_len

USE netcdf, ONLY: nf90_inquire_variable, nf90_put_var

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Writes a 7d array to a variable in the given file
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ncdf), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, DIMENSION(:,:,:,:,:,:,:), INTENT(IN) :: values  ! The values to write
INTEGER, INTENT(IN) :: start(max_dim_var)  ! The point to start writing
                                 ! (one value for each NON-RECORD
                                 ! dimension of the variable)
INTEGER, INTENT(IN) :: counter(max_dim_var)  ! The number of points to write
                                 ! in each dimension of the variable


! Work variables
TYPE(var_ncdf) :: var  ! Structure containing information about the
                       ! variable

INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                       ! has, including record dimension

INTEGER, ALLOCATABLE :: local_start(:)
INTEGER, ALLOCATABLE :: local_count(:)

INTEGER :: error, error_name  ! Error indicators

CHARACTER(LEN=max_sdf_name_len) :: var_name  ! Name of variable.

!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
var = FILE%vars(var_id)

ndims_tot = var%ndims
IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we write to the correct record, if required
ALLOCATE(local_start(ndims_tot))
local_start(1:var%ndims) = start(1:var%ndims)
IF ( var%is_record ) local_start(ndims_tot) = FILE%current_record

! Set up the local_count array - this consists of one value for each
! dimension of the variable from the given counter array, plus a value so that
! we only write to the correct record, if required
ALLOCATE(local_count(ndims_tot))
local_count(1:var%ndims) = counter(1:var%ndims)
IF ( var%is_record ) local_count(ndims_tot) = 1

! Try to write the variable
error = nf90_put_var(FILE%id, var%id, values, local_start, local_count)
IF ( error /= nf90_noerr ) THEN
  !   Get variable name.
  error_name = nf90_inquire_variable(FILE%id, var%id, NAME = var_name )
  IF ( error_name /= nf90_noerr ) var_name = 'unknown' 
  CALL log_fatal_ncdf("file_ncdf_write_var_7d",                               &
                      TRIM(FILE%NAME) // ": " //                              &
                      "Error writing variable '" // TRIM(var_name) // "'",    &
                       error)
END IF

DEALLOCATE(local_start)
DEALLOCATE(local_count)

RETURN

END SUBROUTINE file_ncdf_write_var_7d
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ncdf_close(FILE)

USE netcdf, ONLY: nf90_close
USE io_constants, ONLY: mode_write

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
TYPE(file_ncdf), INTENT(INOUT) :: FILE ! The file to close


! Work variables
INTEGER :: error ! Error code for any errors that occur

!-----------------------------------------------------------------------------

CALL log_info("file_ncdf_close",                                              &
              "Closing file " // TRIM(FILE%NAME))

IF ( FILE%mode == mode_write )                                                &
  CALL file_ncdf_sync(FILE)

error = nf90_close(FILE%id)
IF ( error /= nf90_noerr )                                                    &
  CALL log_fatal_ncdf("file_ncdf_close",                                      &
                      "Error closing file " // TRIM(FILE%NAME), error)

RETURN

END SUBROUTINE file_ncdf_close
!******************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT****************************************

SUBROUTINE file_ncdf_sync(FILE)

USE netcdf, ONLY: nf90_sync

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Syncs the file to a disk
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ncdf), INTENT(IN) :: FILE     ! The file to sync

! Local variables
INTEGER :: error

!-----------------------------------------------------------------------------

error = nf90_sync(FILE%id)

IF ( error /= nf90_noerr )                                                    &
  CALL log_fatal_ncdf("file_ncdf_sync",                                       &
                      "Error syncing file " // TRIM(FILE%NAME), error)

RETURN

END SUBROUTINE file_ncdf_sync

! Inquiry routines
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ncdf_introspect(FILE)

USE io_constants, ONLY: max_sdf_name_len, max_dim_file, mode_read

USE netcdf, ONLY: nf90_inquire, nf90_inquire_dimension, nf90_inquire_variable

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Given an open file in read mode, try to detect the dimensions and
!   variables in the file and define them on the file
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ncdf), INTENT(INOUT) :: FILE  ! The file to detect items in

! Work variables

! Information about the dimensions found in the file
INTEGER :: ndims
INTEGER :: record_dim_id
CHARACTER(LEN=max_sdf_name_len) :: dim_name
INTEGER :: dim_len

! Information about the variables found in the file
INTEGER :: nvars
CHARACTER(LEN=max_sdf_name_len) :: var_name
INTEGER :: var_ndims
INTEGER, ALLOCATABLE :: var_dim_ids(:)

INTEGER :: i  ! Index variables

INTEGER :: dummy  ! Throwaway value

INTEGER :: error  ! Error indicator


!-----------------------------------------------------------------------------

! We can only do this in read mode
IF ( FILE%mode /= mode_read )                                                 &
  CALL log_fatal("file_ncdf_introspect",                                      &
                 TRIM(FILE%NAME) // ": " //                                   &
                 "Can only introspect files in read mode")

!-----------------------------------------------------------------------------
! For NetCDF files, we can use netcdf library functions to recover the
! information that we need
!-----------------------------------------------------------------------------
! First find out general information about the file
error = nf90_inquire(FILE%id, nDimensions = ndims, nVariables = nvars,        &
                              unlimitedDimId = record_dim_id)
IF ( error /= nf90_noerr )                                                    &
  CALL log_fatal_ncdf("file_ncdf_introspect",                                 &
                      TRIM(FILE%NAME) // ": " //                              &
                      "Error reading numbers of dimensions and variables",    &
                      error)

! Check if we have code to handle the number of dimensions
IF ( ndims > max_dim_file )                                                   &
  CALL log_fatal("file_ncdf_introspect",                                      &
                 TRIM(FILE%NAME) // ": " //                                   &
                 "Too many dimensions in file - try increasing max_dim_file")

! Check if we have code to handle the number of variables
IF ( nvars > max_var_file )                                                   &
  CALL log_fatal("file_ncdf_introspect",                                      &
                 TRIM(FILE%NAME) // ": " //                                   &
                 "Too many variables in file - try increasing max_var_file")


!-----------------------------------------------------------------------------
! For each dimension, gather information about it and then define it
!-----------------------------------------------------------------------------
DO i = 1,ndims
  error = nf90_inquire_dimension(FILE%id, i, dim_name, dim_len)
  IF ( error /= nf90_noerr )                                                  &
    CALL log_fatal_ncdf("file_ncdf_introspect",                               &
                        TRIM(FILE%NAME) // ": " //                            &
                        "Error reading dimension info", error)

  ! If we hit the record dimension, then define it as such
  IF ( i == record_dim_id ) THEN
    dummy = file_ncdf_def_record_dim(FILE, dim_name)
    CYCLE
  END IF

  ! Otherwise, define a regular dimension. We don't need to collect the ids
  ! since this driver uses the NetCDF ids for dimensions to define variables
  dummy = file_ncdf_def_dim(FILE, dim_name, dim_len)
END DO


!-----------------------------------------------------------------------------
! For each variable, gather information about it and then define it
!-----------------------------------------------------------------------------
DO i = 1,nvars
  ! Retrieve the variable name and number of dimensions
  error = nf90_inquire_variable(FILE%id, i, NAME = var_name, ndims = var_ndims)
  IF ( error /= nf90_noerr )                                                  &
    CALL log_fatal_ncdf("file_ncdf_introspect",                               &
                        TRIM(FILE%NAME) // ": " //                            &
                        "Error reading variable name", error)

  ! Now we know how many dimensions there are, we can retrieve their ids
  ALLOCATE(var_dim_ids(var_ndims))

  error = nf90_inquire_variable(FILE%id, i, dimids = var_dim_ids)
  IF ( error /= nf90_noerr )                                                  &
    CALL log_fatal_ncdf("file_ncdf_introspect",                               &
                      TRIM(FILE%NAME) // ": " //                              &
                      "Error reading variable dimensions", error)

  ! Now we can define the variable on the file - we don't need to collect
  ! variable ids
  dummy = file_ncdf_def_var(                                                  &
    FILE, var_name,                                                           &
  ! The dimensions that we give are the dimensions we found with the record dim
  ! removed if present
        PACK(var_dim_ids, var_dim_ids /= record_dim_id),                      &
  ! is_record is .TRUE. if the record dimension is present in the list
        ANY(var_dim_ids == record_dim_id)                                     &
      )

  DEALLOCATE(var_dim_ids)
END DO

RETURN

END SUBROUTINE file_ncdf_introspect
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ncdf_inquire_dim(FILE, dim_name, dim_id, dim_len, is_record_dim)

USE netcdf, ONLY: nf90_inq_dimid, nf90_inquire_dimension

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Given a dimension name, returns its id and its size
!   If the returned id < 0, the dimension doesn't exist
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ncdf), INTENT(IN) :: FILE  ! The file to inspect
CHARACTER(LEN=*), INTENT(IN) :: dim_name  ! The name of the dimension to
                                          ! inquire about
INTEGER, INTENT(OUT) :: dim_id  ! The id of the dimension in the file
INTEGER, INTENT(OUT) :: dim_len  ! The length of the dimension in the
LOGICAL, INTENT(OUT) :: is_record_dim  ! Indicates if the named dimension
                                       ! is the record dimension for the file

! Work variables
INTEGER :: error  ! The current error code (if any)


!-----------------------------------------------------------------------------


is_record_dim = .FALSE.

! To recover all the required information about a dimension in a NetCDF file
! needs two steps

! First, recover the dimension id from the name
error = nf90_inq_dimid(FILE%id, dim_name, dim_id)
IF ( error /= nf90_noerr ) THEN
  ! If there was an error retrieving the dimension id, return a non-existent
  ! value
  dim_id = -1
  RETURN
END IF

! If the dimension is the record dimension, we can return early
IF ( dim_id == FILE%record_dim ) THEN
  is_record_dim = .TRUE.
  RETURN
END IF

! If we have a non-record dim, then recover the dim size from the id
error = nf90_inquire_dimension(FILE%id, dim_id, LEN = dim_len)
IF ( error /= nf90_noerr )                                                    &
  CALL log_fatal_ncdf("file_ncdf_inquire_dim",                                &
                      TRIM(FILE%NAME) // ": "     //                          &
                      "dim_name =" // TRIM(dim_name) // ": " //               &
                      "Error recovering dimension size", error)

RETURN

END SUBROUTINE file_ncdf_inquire_dim
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ncdf_inquire_var(FILE, var_name, var_id, ndims, dim_ids, is_record)

USE netcdf, ONLY: nf90_inq_varid, nf90_inquire_variable

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Given a variable name, returns its id and information about dimensions
!   If the returned id < 0, the variable doesn't exist
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ncdf), INTENT(IN) :: FILE  ! The file to inspect
CHARACTER(LEN=*), INTENT(IN) :: var_name  ! The name of the variable to
                                          ! inquire about
INTEGER, INTENT(OUT) :: var_id  ! The id of the variable in the file
INTEGER, INTENT(OUT) :: ndims  ! The number of dimensions that
                               ! the variable has
INTEGER, INTENT(OUT) :: dim_ids(:)  ! The dimension ids
LOGICAL, INTENT(OUT) :: is_record  ! Indicates if the variable
                                   ! uses the record dimension

! Work variables
INTEGER :: var_id_nc  ! The variable id in the NetCDF file

INTEGER :: i, j  ! Index variables

INTEGER :: error  ! The current error code (if any)


!-----------------------------------------------------------------------------

! The default is to assume that we won't find the variable
var_id = -1

! To recover all the required information about a variable in a NetCDF file
! needs two steps (assuming we have been provided with enough space in the
! dim_ids array)

! First, recover the variable id in the file from the name
error = nf90_inq_varid(FILE%id, var_name, var_id_nc)
! If there was an error retrieving the variable id, then the variable doesn't
! exist in the underlying file and we have nothing more to do
IF ( error /= nf90_noerr ) RETURN

! The actual returned var_id is the index in the file%vars array
! So convert it now, in case the variable we just found is one that we are
! not 'technically' aware of (i.e. it has not been specified using file_def_var)
DO i = 1,FILE%nvars
  IF ( FILE%vars(i)%id == var_id_nc ) THEN
    var_id = i
    EXIT
  END IF
END DO

! If we are not aware of the variable, we are done
IF ( var_id < 1 ) RETURN

! Now we can get the information about dimensions
! We assume we have been allocated enough space to return the dimension ids
error = nf90_inquire_variable(FILE%id, var_id_nc, ndims = ndims, dimids = dim_ids)
IF ( error /= nf90_noerr )                                                    &
  CALL log_fatal_ncdf("file_ncdf_inquire_var",                                &
                      TRIM(FILE%NAME) // ": " //                              &
                      "var_name =" // TRIM(var_name) // ": " //               &
                      "Error recovering variable dimensions", error)

! If the record dimension appears in the list of dimension ids, then we
! need to remove it
j = 0
DO i = 1,ndims
  IF ( dim_ids(i) == FILE%record_dim ) CYCLE

  j = j + 1
  dim_ids(j) = dim_ids(i)
END DO
! We will have j dimensions left
ndims = j

! Lastly, indicate if the variable uses the record dimension
is_record = FILE%vars(var_id)%is_record

RETURN

END SUBROUTINE file_ncdf_inquire_var

! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE log_fatal_ncdf(originating_proc, message, ncdf_err)

USE netcdf, ONLY: nf90_strerror

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Logs a NetCDF error
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
CHARACTER(LEN=*), INTENT(IN) :: originating_proc  ! The procedure where
                                                  ! the error occured
CHARACTER(LEN=*), INTENT(IN) :: message  ! Message to accompany error
INTEGER, INTENT(IN) :: ncdf_err  ! The NetCDF error code


!-----------------------------------------------------------------------------

CALL log_fatal(originating_proc, TRIM(message) //                             &
               " (NetCDF error - " // TRIM(nf90_strerror(ncdf_err)) // ")")

RETURN

END SUBROUTINE log_fatal_ncdf

END MODULE driver_ncdf_mod
