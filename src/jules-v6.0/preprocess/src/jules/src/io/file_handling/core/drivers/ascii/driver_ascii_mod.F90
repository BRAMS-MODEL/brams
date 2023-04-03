! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE driver_ascii_mod

USE io_constants, ONLY: max_dim_file, max_var_file, max_attr_file,            &
                         max_file_name_len, max_sdf_name_len,                 &
                         max_attr_val_len, max_dim_var

USE string_utils_mod, ONLY: to_string

USE logging_mod, ONLY: log_info, log_fatal

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------------
CHARACTER(LEN=15) :: extensions(3)       ! The file name extensions recognised
DATA extensions / 'asc', 'txt', 'dat' /  ! for ASCII files

INTEGER, PARAMETER :: max_unit_number = 80

CHARACTER(LEN=1) :: comment_chars(2)
DATA comment_chars / '!', '#' /


!-----------------------------------------------------------------------------
! Type definitions
!-----------------------------------------------------------------------------
TYPE file_ascii
  !-----------------------------------------------------------------------------
  ! This type contains information required to process an ASCII file
  !-----------------------------------------------------------------------------
  INTEGER :: UNIT  ! The unit number that the file is open on

  CHARACTER(LEN=max_file_name_len) :: NAME  ! The name of the file
                                            ! This exists purely for
                                            ! debugging purposes

  INTEGER :: mode ! The mode the file is open in
                  ! This can be one of mode_read or mode_write


  LOGICAL :: define_mode = .TRUE.  ! Determines if file is in define mode

  INTEGER :: ndims = 0  ! The number of dimensions in the file
  CHARACTER(LEN=max_sdf_name_len) :: dim_names(max_dim_file)
                                      ! The names of the dimensions
  INTEGER :: dim_sizes(max_dim_file)  ! The sizes of the dimensions

  LOGICAL :: has_record_dim = .FALSE.  ! Indicates if the file has a
                                       ! record dimension
  CHARACTER(LEN=max_sdf_name_len) :: record_dim_name
                                       ! The name of the record dimension

  INTEGER :: nvars = 0  ! The number of variables in the file
  CHARACTER(LEN=max_sdf_name_len) :: var_names(max_var_file)
                        ! The names of the variables
  INTEGER :: var_sizes(max_var_file)  ! The sizes of the variables
  INTEGER :: var_offsets(max_var_file)  ! The offset of each variable
                                        ! in the buffer
  INTEGER :: var_ndims(max_var_file)  ! The number of dimensions that
                                      ! each variable has
  INTEGER :: var_dims(max_var_file, max_dim_var)
                                      ! The dimension ids for each variable

  INTEGER :: nattrs = 0
  INTEGER :: attr_var_ids(max_attr_file) ! The id of the variable that each
                                         ! attribute is associated with
  CHARACTER(LEN=max_sdf_name_len + max_attr_val_len+3) :: attr_values(max_attr_file)
                                         ! The representation of each attribute
                                         ! that will be printed to file
                                         ! This is of the form
                                         !   name = value

  INTEGER :: record_len  ! The record length to use when reading lines from
                         ! this file
  REAL, POINTER :: buffer(:)  ! Stores data for all variables for the
                              ! current record (i.e. line)

  LOGICAL :: buffer_is_dirty = .FALSE.  ! Only used in write mode
                                        ! Indicates if the buffer has been written
                                        ! to

END TYPE file_ascii


!-----------------------------------------------------------------------------
! Module variables
!-----------------------------------------------------------------------------
CHARACTER(LEN=20) :: out_format_str  ! The format string to use to write
                                     ! output. This is calculated once all
                                     ! variables have been defined


!-----------------------------------------------------------------------------
! Interface definitions for the overloads of read/write var for variables
! of different ranks
!-----------------------------------------------------------------------------
INTERFACE file_ascii_read_var
MODULE PROCEDURE file_ascii_read_var_scalar, file_ascii_read_var_1d,          &
                 file_ascii_read_var_2d,     file_ascii_read_var_3d,          &
                 file_ascii_read_var_4d,     file_ascii_read_var_5d,          &
                 file_ascii_read_var_6d,     file_ascii_read_var_7d
END INTERFACE file_ascii_read_var

INTERFACE file_ascii_write_var
MODULE PROCEDURE file_ascii_write_var_scalar, file_ascii_write_var_1d,        &
                 file_ascii_write_var_2d,     file_ascii_write_var_3d,        &
                 file_ascii_write_var_4d,     file_ascii_write_var_5d,        &
                 file_ascii_write_var_6d,     file_ascii_write_var_7d
END INTERFACE file_ascii_write_var


!-----------------------------------------------------------------------------
! Visibility declarations
!-----------------------------------------------------------------------------
PRIVATE
PUBLIC extensions,                                                            &
! Types
         file_ascii,                                                          &
! Routines for opening and closing
         file_ascii_open, file_ascii_close,                                   &
! Routines for definitions
         file_ascii_def_dim, file_ascii_def_record_dim, file_ascii_def_var,   &
         file_ascii_def_attr_real, file_ascii_def_attr_int,                   &
         file_ascii_def_attr_char, file_ascii_enddef,                         &
! Routines for seeking
         file_ascii_seek, file_ascii_advance,                                 &
! Routines for reading and writing
         file_ascii_read_var, file_ascii_write_var,                           &
! Routines for inquiry
         file_ascii_introspect, file_ascii_inquire_dim, file_ascii_inquire_var


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


FUNCTION file_ascii_open(NAME, mode, comm, info) RESULT(FILE)

USE io_constants, ONLY: mode_read, mode_write, unit_stdin, unit_stdout
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Opens a file and returns a file_ascii object representing it
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
                                       ! If this is given, an error is raised
                                       ! as it is not possible to use parallel
                                       ! access with this ASCII driver
INTEGER, OPTIONAL, INTENT(IN) :: info  ! MPI info object to use for parallel
                                       ! I/O
                                       ! If this is given, an error is raised
                                       ! as it is not possible to use parallel
                                       ! access with this ASCII driver

!-----------------------------------------------------------------------------
! We know from the conditions imposed by file_open that comm and info are
! either both present or both not present
!-----------------------------------------------------------------------------

! Return type
TYPE(file_ascii) :: FILE


! Work variables
INTEGER :: UNIT ! The unit number for the opened file
INTEGER :: error ! Error code for any errors that occur
INTEGER :: ntasks  ! The number of MPI tasks

LOGICAL :: unit_in_use ! Used in call to INQUIRE to detect if a unit is in use
INTEGER :: i  ! Loop counter
  
CHARACTER(LEN=errormessagelength) :: iomessage

!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! If comm and info are specified and the number of available tasks is > 1,
! raise an error since parallel access is not possible with this ASCII driver
!-----------------------------------------------------------------------------
IF ( PRESENT(comm) ) THEN
  CALL mpi_comm_size(comm, ntasks, error)

  IF ( ntasks > 1 )                                                           &
    CALL log_fatal("file_ascii_open",                                         &
                   "Parallel access is not available for ASCII files. " //    &
                   "File: " // TRIM(NAME) )
END IF

!-----------------------------------------------------------------------------
! Get a spare unit to open the file
!-----------------------------------------------------------------------------
UNIT = 0
! Search the allowed unit numbers until we find one with nothing connected
DO i = 1,max_unit_number
  ! Avoid units for standard i/o
  IF ( i == unit_stdin .OR. i == unit_stdout ) CYCLE

  ! Find out if anything is connected to this unit.
  INQUIRE(UNIT = i, OPENED = unit_in_use )
  IF ( .NOT. unit_in_use ) THEN
    ! This unit is free, so exit the loop
    UNIT = i
    EXIT
  END IF
END DO

IF ( UNIT < 1 ) THEN
  CALL log_fatal("file_ascii_open",                                           &
                 "All allowed units are in use - try increasing " //          &
                 "MAX_UNIT_NUMBER")
END IF


!-----------------------------------------------------------------------------
! Open the file in the requested mode
!-----------------------------------------------------------------------------
SELECT CASE ( mode )
CASE ( mode_read )
  CALL log_info("file_ascii_open",                                            &
                "Opening file " // TRIM(NAME) // " for reading")
  ! Open file for reading only - file must already exist
  OPEN(UNIT, FILE=NAME, STATUS='old', POSITION='rewind', ACTION='read',       &
             IOSTAT = error, IOMSG = iomessage)

CASE ( mode_write )
  CALL log_info("file_ascii_open",                                            &
                "Opening file " // TRIM(NAME) // " for writing")
  ! Create an empty file for writing only - if a file with the
  ! same name exists, overwrite it
  OPEN(UNIT, FILE=NAME, STATUS='replace', POSITION='rewind',                  &
             ACTION='write', IOSTAT = error)

CASE DEFAULT
  ! Read and write are the only supported modes
  CALL log_fatal("file_ascii_open",                                           &
                 "Unsupported mode - " // TRIM(to_string(mode)))

END SELECT

!-----------------------------------------------------------------------------
! Report any errors
!-----------------------------------------------------------------------------
IF ( error /= 0 )                                                             &
  CALL log_fatal("file_ascii_open",                                           &
                 "Error opening file " // TRIM(NAME) //                       &
                 " (IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //        &
                 TRIM(iomessage) // ")")


! Initialise the file_ascii object
FILE%NAME = NAME
FILE%mode = mode
FILE%UNIT = UNIT

RETURN

END FUNCTION file_ascii_open
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_ascii_def_dim(FILE, dim_name, dim_len) RESULT(dim_id)

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
TYPE(file_ascii), INTENT(INOUT) :: FILE
                                ! The file to define the dimension on
CHARACTER(LEN=*), INTENT(IN) :: dim_name
                                ! The name of the dimension
INTEGER, INTENT(IN) :: dim_len  ! The length of the dimension

! Return type
INTEGER :: dim_id               ! The dimension id


!-----------------------------------------------------------------------------

IF ( .NOT. FILE%define_mode )                                                 &
  CALL log_fatal("file_ascii_def_dim",                                        &
                 "Cannot define dimension - file is not in define mode")

! If adding another dimension will cause us to have too many dimensions,
! error out
IF ( FILE%ndims >= max_dim_file )                                             &
  CALL log_fatal("file_ascii_def_dim",                                        &
                 "Too many dimensions in file - try increasing max_dim_file")

! Check the dimension has a sensible size (i.e. > 0)
IF ( dim_len < 1 )                                                            &
  CALL log_fatal("file_ascii_def_dim",                                        &
                 "Dimension size should be >= 1 - size of " //                &
                 TRIM(to_string(dim_len)) // " given for dimension " //       &
                 TRIM(dim_name))

!-----------------------------------------------------------------------------
! Store the dimension info to use later
!-----------------------------------------------------------------------------
FILE%ndims = FILE%ndims + 1

! The returned dimension id is just the index in the dims array on the file object
dim_id = FILE%ndims

FILE%dim_names(dim_id) = dim_name
FILE%dim_sizes(dim_id) = dim_len

RETURN

END FUNCTION file_ascii_def_dim
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_ascii_def_record_dim(FILE, dim_name) RESULT(dim_id)

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
TYPE(file_ascii), INTENT(INOUT) :: FILE
                                ! The file to define the dimension on
CHARACTER(LEN=*), INTENT(IN) :: dim_name
                                ! The name of the record dimension

! Return type
INTEGER :: dim_id               ! The dimension id


!-----------------------------------------------------------------------------

IF ( .NOT. FILE%define_mode )                                                 &
  CALL log_fatal("file_ascii_def_record_dim",                                 &
                 "Cannot define record dimension - file is not in define mode")

! Just indicate that the file has a record dimension
FILE%has_record_dim = .TRUE.
FILE%record_dim_name = dim_name

! Return -1 as the dimension id, since the return value should never be used
! and -1 can never be a valid dimension id for a non-record dimension
dim_id = -1

RETURN

END FUNCTION file_ascii_def_record_dim
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_ascii_def_var(FILE, var_name, dims, is_record) RESULT(var_id)

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
TYPE(file_ascii), INTENT(INOUT) :: FILE
                                  ! The file to define the variable in
CHARACTER(LEN=*), INTENT(IN) :: var_name
                                  ! The name of the variable
INTEGER, INTENT(IN) :: dims(:)    ! The ids of the NON-RECORD dimensions of
                                  ! the variable
LOGICAL, INTENT(IN) :: is_record  ! Indicates whether the variable uses the
                                  ! record dimension

! Return type
INTEGER :: var_id                 ! The variable id - this is the index
                                  ! in the vars array of the file_ascii object
                                  ! of the variable


! Work variables
INTEGER :: dim_sizes(SIZE(dims))  ! The size of each specified dimension

INTEGER :: i  ! Loop variable


!-----------------------------------------------------------------------------


! If we are not in define mode, error out
IF ( .NOT. FILE%define_mode )                                                 &
  CALL log_fatal("file_ascii_def_var",                                        &
                 "Cannot define variable - file is not in define mode")

! If adding another variable will cause us to have too many variables,
! error out
IF ( FILE%nvars >= max_var_file )                                             &
  CALL log_fatal("file_ascii_def_var",                                        &
                 "Too many variables in file - try increasing max_var_file")

! The record dimension needs to be defined before any variables can use it
IF ( is_record .AND. .NOT. FILE%has_record_dim )                              &
  CALL log_fatal("file_ascii_def_var",                                        &
                 "Record dimension must be defined before defining a " //     &
                 "variable as using it")

! If the record dimension is defined, all variables must use it
IF ( FILE%has_record_dim .AND. .NOT. is_record )                              &
  CALL log_fatal("file_ascii_def_var",                                        &
                 "Record dimension has been defined, so all variables " //    &
                 "must use it")

! Check that all the given dimensions are valid
IF ( ANY(dims < 1) )                                                          &
  CALL log_fatal("file_ascii_def_var",                                        &
                 "Negative dimension id given - note that use of the " //     &
                 "record dimension should be specified using 'is_record'")

! Get the dimension sizes from the given ids
DO i = 1,SIZE(dims)
  dim_sizes(i) = FILE%dim_sizes( dims(i) )
END DO


!-----------------------------------------------------------------------------
! Save information about the variable for later use
!-----------------------------------------------------------------------------
FILE%nvars = FILE%nvars + 1

! The return value is the index of the variable in the vars array on the
! file_ascii object
var_id = FILE%nvars

FILE%var_names(var_id) = var_name
FILE%var_sizes(var_id) = PRODUCT(dim_sizes)

FILE%var_ndims(var_id) = SIZE(dims)
FILE%var_dims(var_id,1:SIZE(dims)) = dims(:)

RETURN

END FUNCTION file_ascii_def_var
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ascii_def_attr_real(FILE, var_id, NAME, VALUE)

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
TYPE(file_ascii), INTENT(INOUT) :: FILE
                                  ! The file to define the attribute in
INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                  ! attribute on
CHARACTER(LEN=*), INTENT(IN) :: NAME
                                  ! The name of the attribute
REAL, INTENT(IN) :: VALUE         ! The value of the attribute


!-----------------------------------------------------------------------------

! Since all attributes on ASCII files will eventually be char when written
! to file, we just delegate to def_attr_char with a converted value
CALL file_ascii_def_attr_char(FILE, var_id, NAME, TRIM(to_string(VALUE)))

RETURN

END SUBROUTINE file_ascii_def_attr_real


SUBROUTINE file_ascii_def_attr_int(FILE, var_id, NAME, VALUE)

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
TYPE(file_ascii), INTENT(INOUT) :: FILE
                                  ! The file to define the attribute in
INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                  ! attribute on
CHARACTER(LEN=*), INTENT(IN) :: NAME
                                  ! The name of the attribute
INTEGER, INTENT(IN) :: VALUE      ! The value of the attribute


!-----------------------------------------------------------------------------

! Since all attributes on ASCII files will eventually be char when written
! to file, we just delegate to def_attr_char with a converted value
CALL file_ascii_def_attr_char(FILE, var_id, NAME, TRIM(to_string(VALUE)))

RETURN

END SUBROUTINE file_ascii_def_attr_int


SUBROUTINE file_ascii_def_attr_char(FILE, var_id, NAME, VALUE)

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
TYPE(file_ascii), INTENT(INOUT) :: FILE
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
  CALL log_fatal("file_ascii_def_attr_char",                                  &
                 "Cannot define attribute - file is not in define mode")

! If adding another attribute will cause us to have too many attributes,
! error out
IF ( FILE%nattrs >= max_attr_file )                                           &
  CALL log_fatal("file_ascii_def_attr_char",                                  &
                 "Too many attributes in file - try increasing max_attr_file")

!-----------------------------------------------------------------------------
! Store information about the attribute so that it can be used later
!-----------------------------------------------------------------------------
FILE%nattrs = FILE%nattrs + 1

FILE%attr_var_ids(FILE%nattrs) = var_id
! We trim trailing spaces from name and leading spaces from value so that
! they will be the correct distance apart. Then left align the whole lot
FILE%attr_values(FILE%nattrs)  = ADJUSTL(TRIM(NAME) // " = " // ADJUSTL(VALUE))

RETURN

END SUBROUTINE file_ascii_def_attr_char
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ascii_enddef(FILE)

USE io_constants, ONLY: mode_read, mode_write, attr_global

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
TYPE(file_ascii), INTENT(INOUT) :: FILE

! Work variables
CHARACTER(LEN=100) :: var_dims_str  ! Used when writing headers to store
                                    ! variables dimension names as a
                                    ! string for writing

INTEGER :: buffer_size  ! The size of the buffer
INTEGER :: offset  ! Used in calculation of offsets
INTEGER :: error  ! Error indicator

INTEGER :: i,j  ! Loop counters


!-----------------------------------------------------------------------------

FILE%define_mode = .FALSE.

!-----------------------------------------------------------------------------
! Allocate the buffer - it must be big enough to hold one records worth of
! data for each variable
!-----------------------------------------------------------------------------
! To get the size of the buffer, just sum the sizes of the defined variables
buffer_size = SUM(FILE%var_sizes(1:FILE%nvars))
ALLOCATE(FILE%buffer(buffer_size), stat = error)
IF ( error /= 0 )                                                             &
  CALL log_fatal("file_ascii_enddef",                                         &
                 "Error allocating memory - " //                              &
                 "(STAT=" // TRIM(to_string(error)) // ")")

! Now we know the size of the buffer, we can make a better guess at the
! record length we will need. For now, we'll assume no more than 30 characters
! per buffer element - this should be very conservative
FILE%record_len = buffer_size * 30

!-----------------------------------------------------------------------------
! Calculate offsets for each variable
! It is assumed that they appear in the file in the order in which they
! were defined
!-----------------------------------------------------------------------------
offset = 1
DO i = 1,FILE%nvars
  FILE%var_offsets(i) = offset
  offset = offset + FILE%var_sizes(i)
END DO


! Now we do different things depending on mode
SELECT CASE ( FILE%mode )
CASE ( mode_read )
  ! In read mode, we fill the buffer with the first data so we are ready to read
  CALL file_ascii_fill_buffer(FILE)

CASE ( mode_write )
  !-----------------------------------------------------------------------------
  ! In write mode, we write the headers describing the variables in the file
  !
  ! The header is of the form:
  !
  ! #
  ! # CREATED BY JULES LSM
  ! #
  ! # Global attributes:
  ! #     global_attr_1 = foobar ;
  ! #     global_attr_2 = baz ;
  ! #
  ! # Dimensions:
  ! #     dim1 = size1 ;
  ! #     dim2 = size2 ;
  ! #     dim3 = size3 ;
  ! #     dim4 = UNLIMITED ;
  ! #
  ! # Variables:
  ! #     variable1(dim1,dim2) ;
  ! #         variable1_attr_1 = foobaz ;
  ! #         variable1_attr_2 = blah ;
  ! #
  ! #     variable2(dim1) ;
  ! #
  ! #     variable3(dim1,dim2,dim3) ;
  ! #
  !-----------------------------------------------------------------------------
  ! First write the created by header
  WRITE(FILE%UNIT, "(A)") "#"
  WRITE(FILE%UNIT, "(A)") "# CREATED BY JULES LSM"
  WRITE(FILE%UNIT, "(A)") "#"

  !-----------------------------------------------------------------------------
  ! Write global attributes
  !-----------------------------------------------------------------------------
  WRITE(FILE%UNIT, "(A)") "# Global attributes:"
  DO i = 1,FILE%nattrs
    IF ( FILE%attr_var_ids(i) == attr_global )                                &
      WRITE(FILE%UNIT, "(A)") "#     " // TRIM(FILE%attr_values(i)) // " ;"
  END DO

  WRITE(FILE%UNIT, "(A)") "#"

  !-----------------------------------------------------------------------------
  ! Write information about the dimensions in the file
  !-----------------------------------------------------------------------------
  WRITE(FILE%UNIT, "(A)") "# Dimensions:"
  DO i = 1,FILE%ndims
    WRITE(FILE%UNIT, "(A)") "#     " // TRIM(FILE%dim_names(i)) //            &
                            " = " // TRIM(to_string(FILE%dim_sizes(i))) // " ;"

  END DO
  ! Write information about the record dimension
  IF ( FILE%has_record_dim )                                                  &
    WRITE(FILE%UNIT, "(A)") "#     " // TRIM(FILE%record_dim_name) //         &
                            " = UNLIMITED ;"

  WRITE(FILE%UNIT, "(A)") "#"

  !-----------------------------------------------------------------------------
  ! Write information about the variables in the file
  !-----------------------------------------------------------------------------
  WRITE(FILE%UNIT, "(A)") "# Variables:"
  DO i = 1,FILE%nvars
    ! Build the dimension names string for the variable
    var_dims_str = ""
    DO j = 1,FILE%var_ndims(i)
      var_dims_str = TRIM(var_dims_str) // FILE%dim_names( FILE%var_dims(i,j) )

      ! Add a comma unless this is the last dimension that the variable has, including
      ! the addition of the record dim below
      IF ( j < FILE%var_ndims(i) .OR. FILE%has_record_dim )                   &
        var_dims_str = TRIM(var_dims_str) // ","
    END DO
    ! Add the record dimension if it is in use
    IF ( FILE%has_record_dim ) THEN
      var_dims_str = TRIM(var_dims_str) // FILE%record_dim_name
    END IF

    WRITE(FILE%UNIT, "(A)") "#     " // TRIM(FILE%var_names(i)) //            &
                            "(" // TRIM(var_dims_str) // ") ;"

    ! Write info about each of the variables attributes
    DO j = 1,FILE%nattrs
      IF ( FILE%attr_var_ids(j) == i )                                        &
        WRITE(FILE%UNIT, "(A)") "#         " // TRIM(FILE%attr_values(j)) // " ;"
    END DO

    WRITE(FILE%UNIT, "(A)") "#"
  END DO

  ! Set the format string to use for output, now we know how many elements are in
  ! the buffer
  out_format_str = "(" // TRIM(to_string(buffer_size)) // "G15.7e2)"

END SELECT

RETURN

END SUBROUTINE file_ascii_enddef
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ascii_seek(FILE, record)

USE io_constants, ONLY: mode_read
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Seeks the file to a particular record (i.e. the next time
!   values are read from the file using file_read_var, they will be read from
!   the requested record)
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ascii), INTENT(INOUT) :: FILE  ! The file to seek
INTEGER, INTENT(IN) :: record            ! The record number to seek to


! Work variables
INTEGER :: error  ! Error indicator

INTEGER :: i  ! Loop counter
  
CHARACTER(LEN=errormessagelength) :: iomessage

!-----------------------------------------------------------------------------

! We can't seek the file if it's in define mode
IF ( FILE%define_mode )                                                       &
  CALL log_fatal("file_ascii_seek",                                           &
                 "Cannot advance - file is still in define mode")

! Arbitrary seeks are only allowed in read mode
IF ( FILE%mode /= mode_read )                                                 &
  CALL log_fatal("file_ascii_seek",                                           &
                 "Arbitrary seeks are only allowed in read mode")

! We can only seek if there is a record dimension
IF ( .NOT. FILE%has_record_dim )                                              &
  CALL log_fatal("file_ascii_seek",                                           &
                 "Cannot advance file - no record dimension defined")

!-----------------------------------------------------------------------------
! Actually seek the file
! We do this by rewinding the file to the start, then using fill_buffer to
! skip comment lines
!-----------------------------------------------------------------------------
REWIND(FILE%UNIT, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("file_ascii_seek",                                           &
                 "Error rewinding file " //                                   &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

DO i = 1,record
  CALL file_ascii_fill_buffer(FILE)
END DO

! By the time we get to here, the buffer is filled with the data from record
! #record, which is what we want...!

RETURN

END SUBROUTINE file_ascii_seek
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ascii_advance(FILE)

USE io_constants, ONLY: mode_read, mode_write

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
TYPE(file_ascii), INTENT(INOUT) :: FILE  ! The file to advance


!-----------------------------------------------------------------------------

! We can't advance the file if it's in define mode
IF ( FILE%define_mode )                                                       &
  CALL log_fatal("file_ascii_advance",                                        &
                 "Cannot advance - file is still in define mode")

! We can only advance if there is a record dimension
IF ( .NOT. FILE%has_record_dim )                                              &
  CALL log_fatal("file_ascii_advance",                                        &
                 "Cannot advance file - no record dimension defined")


! We do different things depending on mode
SELECT CASE ( FILE%mode )
CASE ( mode_read )
  ! In read mode, we just fill the buffer with the next record of data, so that
  ! the next set of reads will use that data
  CALL file_ascii_fill_buffer(FILE)

CASE ( mode_write )
  ! In write mode, we just flush the buffer, so we can start populating it
  ! with the next records worth of data
  CALL file_ascii_flush_buffer(FILE)

END SELECT

RETURN

END SUBROUTINE file_ascii_advance
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ascii_read_var_scalar(FILE, var_id, VALUE, start)

USE io_constants, ONLY: max_dim_var

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
TYPE(file_ascii), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, INTENT(OUT) :: VALUE  ! The value read from file
INTEGER, INTENT(IN) :: start(:)  ! The point to start reading from
                                 ! (one value for each dimension
                                 ! of the variable in the file)

! Work variables
REAL :: local_values(1)  ! An array version of value
INTEGER :: local_count(max_dim_var)  ! A version of count, which will be
                                     ! all 1s

!-----------------------------------------------------------------------------

! To implement this routine, we just delegate to the 1d routine and reshape
! the result
local_count(:) = 1

CALL file_ascii_read_var_1d(FILE, var_id, local_values, start, local_count)

VALUE = local_values(1)

RETURN

END SUBROUTINE file_ascii_read_var_scalar


SUBROUTINE file_ascii_read_var_1d(FILE, var_id, values, start, count_in)

USE io_constants, ONLY: mode_read

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
TYPE(file_ascii), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, DIMENSION(:), INTENT(OUT) :: values  ! The values read from file
INTEGER, INTENT(IN) :: start(:)  ! The point to start reading from
                                 ! (one value for each dimension
                                 ! of the variable in the file)
INTEGER, INTENT(IN) :: count_in(:)  ! The number of points to read
                                    ! in each dimension of the variable

! Work variables
INTEGER :: ndims  ! The number of dimensions that the variable has in file

INTEGER :: offset  ! The offset to use in the buffer

!-----------------------------------------------------------------------------

! We can't read data if the file is in define mode
IF ( FILE%define_mode )                                                       &
  CALL log_fatal("file_ascii_read_var_1d",                                    &
                 "Cannot read data - file is still in define mode")

! We can't read data unless we are in read mode...
IF ( FILE%mode /= mode_read )                                                 &
  CALL log_fatal("file_ascii_read_var_1d",                                    &
                 "Can only read data if file is opened in read mode")

ndims = FILE%var_ndims(var_id)

! ASCII files have the restriction that we can only read the whole variable
IF ( ANY(start(1:ndims) /= 1) )                                               &
  CALL log_fatal("file_ascii_read_var_1d",                                    &
                 "start must be 1 for all dimensions - reading part of " //   &
                 "a variable from an ASCII file is not supported")

IF ( SIZE(values) /= PRODUCT(count_in(1:ndims)) .AND.                         &
     SIZE(values) /= FILE%var_sizes(var_id) )                                 &
  CALL log_fatal("file_ascii_read_var_1d",                                    &
                 "values must have the same number of elements as the " //    &
                 "variable in file - reading part of a variable from an " //  &
                 "ASCII file is not supported")

! Now we know we can read correctly
offset = FILE%var_offsets(var_id)
values(:) = FILE%buffer(offset:(offset + FILE%var_sizes(var_id) - 1))

RETURN

END SUBROUTINE file_ascii_read_var_1d


SUBROUTINE file_ascii_read_var_2d(FILE, var_id, values, start, count_in)

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
TYPE(file_ascii), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, DIMENSION(:,:), INTENT(OUT) :: values  ! The values read from file
INTEGER, INTENT(IN) :: start(:)  ! The point to start reading from
                                 ! (one value for each dimension
                                 ! of the variable in the file)
INTEGER, INTENT(IN) :: count_in(:)  ! The number of points to read
                                    ! in each dimension of the variable


! Work variables
REAL, ALLOCATABLE :: local_values(:)  ! A 1d version of values

!-----------------------------------------------------------------------------

! To implement this routine, we just delegate to the 1d routine and reshape
! the result

ALLOCATE(local_values(PRODUCT(SHAPE(values))))

CALL file_ascii_read_var_1d(FILE, var_id, local_values, start, count_in)

values = RESHAPE(local_values, SHAPE(values))

DEALLOCATE(local_values)

RETURN

END SUBROUTINE file_ascii_read_var_2d


SUBROUTINE file_ascii_read_var_3d(FILE, var_id, values, start, count_in)

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
TYPE(file_ascii), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, DIMENSION(:,:,:), INTENT(OUT) :: values  ! The values read from file
INTEGER, INTENT(IN) :: start(:)  ! The point to start reading from
                                 ! (one value for each dimension
                                 ! of the variable in the file)
INTEGER, INTENT(IN) :: count_in(:)  ! The number of points to read
                                    ! in each dimension of the variable


! Work variables
REAL, ALLOCATABLE :: local_values(:)  ! A 1d version of values

!-----------------------------------------------------------------------------


! To implement this routine, we just delegate to the 1d routine and reshape
! the result

ALLOCATE(local_values(PRODUCT(SHAPE(values))))

CALL file_ascii_read_var_1d(FILE, var_id, local_values, start, count_in)

values = RESHAPE(local_values, SHAPE(values))

DEALLOCATE(local_values)

END SUBROUTINE file_ascii_read_var_3d


SUBROUTINE file_ascii_read_var_4d(FILE, var_id, values, start, count_in)

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
TYPE(file_ascii), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: values  ! The values read from file
INTEGER, INTENT(IN) :: start(:)  ! The point to start reading from
                                 ! (one value for each dimension
                                 ! of the variable in the file)
INTEGER, INTENT(IN) :: count_in(:)  ! The number of points to read
                                    ! in each dimension of the variable


! Work variables
REAL, ALLOCATABLE :: local_values(:)  ! A 1d version of values

!-----------------------------------------------------------------------------


! To implement this routine, we just delegate to the 1d routine and reshape
! the result

ALLOCATE(local_values(PRODUCT(SHAPE(values))))

CALL file_ascii_read_var_1d(FILE, var_id, local_values, start, count_in)

values = RESHAPE(local_values, SHAPE(values))

DEALLOCATE(local_values)

END SUBROUTINE file_ascii_read_var_4d


SUBROUTINE file_ascii_read_var_5d(FILE, var_id, values, start, count_in)

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
TYPE(file_ascii), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: values  ! The values read from file
INTEGER, INTENT(IN) :: start(:)  ! The point to start reading from
                                 ! (one value for each dimension
                                 ! of the variable in the file)
INTEGER, INTENT(IN) :: count_in(:)  ! The number of points to read
                                    ! in each dimension of the variable


! Work variables
REAL, ALLOCATABLE :: local_values(:)  ! A 1d version of values

!-----------------------------------------------------------------------------


! To implement this routine, we just delegate to the 1d routine and reshape
! the result

ALLOCATE(local_values(PRODUCT(SHAPE(values))))

CALL file_ascii_read_var_1d(FILE, var_id, local_values, start, count_in)

values = RESHAPE(local_values, SHAPE(values))

DEALLOCATE(local_values)

END SUBROUTINE file_ascii_read_var_5d


SUBROUTINE file_ascii_read_var_6d(FILE, var_id, values, start, count_in)

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
TYPE(file_ascii), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, DIMENSION(:,:,:,:,:,:), INTENT(OUT) :: values  ! The values read from file
INTEGER, INTENT(IN) :: start(:)  ! The point to start reading from
                                 ! (one value for each dimension
                                 ! of the variable in the file)
INTEGER, INTENT(IN) :: count_in(:)  ! The number of points to read
                                    ! in each dimension of the variable


! Work variables
REAL, ALLOCATABLE :: local_values(:)  ! A 1d version of values

!-----------------------------------------------------------------------------


! To implement this routine, we just delegate to the 1d routine and reshape
! the result

ALLOCATE(local_values(PRODUCT(SHAPE(values))))

CALL file_ascii_read_var_1d(FILE, var_id, local_values, start, count_in)

values = RESHAPE(local_values, SHAPE(values))

DEALLOCATE(local_values)

END SUBROUTINE file_ascii_read_var_6d


SUBROUTINE file_ascii_read_var_7d(FILE, var_id, values, start, count_in)

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
TYPE(file_ascii), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, DIMENSION(:,:,:,:,:,:,:), INTENT(OUT) :: values  ! The values read from file
INTEGER, INTENT(IN) :: start(:)  ! The point to start reading from
                                 ! (one value for each dimension
                                 ! of the variable in the file)
INTEGER, INTENT(IN) :: count_in(:)  ! The number of points to read
                                    ! in each dimension of the variable


! Work variables
REAL, ALLOCATABLE :: local_values(:)  ! A 1d version of values

!-----------------------------------------------------------------------------


! To implement this routine, we just delegate to the 1d routine and reshape
! the result

ALLOCATE(local_values(PRODUCT(SHAPE(values))))

CALL file_ascii_read_var_1d(FILE, var_id, local_values, start, count_in)

values = RESHAPE(local_values, SHAPE(values))

DEALLOCATE(local_values)

END SUBROUTINE file_ascii_read_var_7d
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ascii_write_var_scalar(FILE, var_id, VALUE, start)

USE io_constants, ONLY: max_dim_var

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
TYPE(file_ascii), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, INTENT(IN) :: VALUE  ! The value to write
INTEGER, INTENT(IN) :: start(:)  ! The point to start writing at
                                 ! (one value for each dimension
                                 ! of the variable in the file)

! Work variables
REAL :: local_values(1)  ! An array version of value
INTEGER :: local_count(max_dim_var)  ! A version of count, which will be
                                     ! all 1s

!-----------------------------------------------------------------------------

! To implement this routine, we just delegate to the 1d routine and reshape
! the result
local_count(:) = 1
local_values(1) = VALUE

CALL file_ascii_write_var_1d(FILE, var_id, local_values, start, local_count)

RETURN

END SUBROUTINE file_ascii_write_var_scalar


SUBROUTINE file_ascii_write_var_1d(FILE, var_id, values, start, count_in)

USE io_constants, ONLY: mode_write

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
TYPE(file_ascii), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, DIMENSION(:), INTENT(IN) :: values  ! The values to write
INTEGER, INTENT(IN) :: start(:)  ! The point to start writing at
                                 ! (one value for each dimension
                                 ! of the variable in the file)
INTEGER, INTENT(IN) :: count_in(:)  ! The number of points to write
                                    ! in each dimension of the variable

! Work variables
INTEGER :: ndims  ! The number of dimensions that the variable has in file

INTEGER :: offset  ! The offset to use in the buffer


!-----------------------------------------------------------------------------

! We can't write data if the file is in define mode
IF ( FILE%define_mode )                                                       &
  CALL log_fatal("file_ascii_write_var_1d",                                   &
                 "Cannot write data - file is still in define mode")

! We can't write data unless we are in write mode...
IF ( FILE%mode /= mode_write )                                                &
  CALL log_fatal("file_ascii_write_var_1d",                                   &
                 "Can only write data if file is opened in write mode")

ndims = FILE%var_ndims(var_id)

! ASCII files have the restriction that we can only write the whole variable
IF ( ANY(start(1:ndims) /= 1) )                                               &
  CALL log_fatal("file_ascii_write_var_1d",                                   &
                 "start must be 1 for all dimensions - writing part of " //   &
                 "a variable in an ASCII file is not supported")

IF ( SIZE(values) /= PRODUCT(count_in(1:ndims)) .AND.                         &
     SIZE(values) /= FILE%var_sizes(var_id) )                                 &
  CALL log_fatal("file_ascii_write_var_1d",                                   &
                 "values must have the same number of elements as the " //    &
                 "variable in file - writing part of a variable in an " //    &
                 "ASCII file is not supported")

! Now we know we can write to the buffer correctly
offset = FILE%var_offsets(var_id)
FILE%buffer(offset:(offset + FILE%var_sizes(var_id) - 1)) = values(:)

! Mark the buffer as dirty (i.e. has writes that have not been committed to file
FILE%buffer_is_dirty = .TRUE.

RETURN

END SUBROUTINE file_ascii_write_var_1d


SUBROUTINE file_ascii_write_var_2d(FILE, var_id, values, start, count_in)

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
TYPE(file_ascii), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, DIMENSION(:,:), INTENT(IN) :: values  ! The values to write
INTEGER, INTENT(IN) :: start(:)  ! The point to start writing at
                                 ! (one value for each dimension
                                 ! of the variable in the file)
INTEGER, INTENT(IN) :: count_in(:)  ! The number of points to write
                                    ! in each dimension of the variable


! Work variables
REAL, ALLOCATABLE :: local_values(:)  ! A 1d version of values

!-----------------------------------------------------------------------------

! To implement this routine, we just delegate to the 1d routine and reshape
! the result

ALLOCATE(local_values(PRODUCT(SHAPE(values))))

local_values(:) = RESHAPE(values, (/ PRODUCT(SHAPE(values)) /))

CALL file_ascii_write_var_1d(FILE, var_id, local_values, start, count_in)

DEALLOCATE(local_values)

RETURN

END SUBROUTINE file_ascii_write_var_2d


SUBROUTINE file_ascii_write_var_3d(FILE, var_id, values, start, count_in)

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
TYPE(file_ascii), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, DIMENSION(:,:,:), INTENT(IN) :: values  ! The values to write
INTEGER, INTENT(IN) :: start(:)  ! The point to start writing at
                                 ! (one value for each dimension
                                 ! of the variable in the file)
INTEGER, INTENT(IN) :: count_in(:)  ! The number of points to write
                                    ! in each dimension of the variable


! Work variables
REAL, ALLOCATABLE :: local_values(:)  ! A 1d version of values

!-----------------------------------------------------------------------------

! To implement this routine, we just delegate to the 1d routine and reshape
! the result

ALLOCATE(local_values(PRODUCT(SHAPE(values))))

local_values(:) = RESHAPE(values, (/ PRODUCT(SHAPE(values)) /))

CALL file_ascii_write_var_1d(FILE, var_id, local_values, start, count_in)

DEALLOCATE(local_values)

RETURN

END SUBROUTINE file_ascii_write_var_3d


SUBROUTINE file_ascii_write_var_4d(FILE, var_id, values, start, count_in)

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
TYPE(file_ascii), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: values  ! The values to write
INTEGER, INTENT(IN) :: start(:)  ! The point to start writing at
                                 ! (one value for each dimension
                                 ! of the variable in the file)
INTEGER, INTENT(IN) :: count_in(:)  ! The number of points to write
                                    ! in each dimension of the variable


! Work variables
REAL, ALLOCATABLE :: local_values(:)  ! A 1d version of values

!-----------------------------------------------------------------------------

! To implement this routine, we just delegate to the 1d routine and reshape
! the result

ALLOCATE(local_values(PRODUCT(SHAPE(values))))

local_values(:) = RESHAPE(values, (/ PRODUCT(SHAPE(values)) /))

CALL file_ascii_write_var_1d(FILE, var_id, local_values, start, count_in)

DEALLOCATE(local_values)

RETURN

END SUBROUTINE file_ascii_write_var_4d


SUBROUTINE file_ascii_write_var_5d(FILE, var_id, values, start, count_in)

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
TYPE(file_ascii), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, DIMENSION(:,:,:,:,:), INTENT(IN) :: values  ! The values to write
INTEGER, INTENT(IN) :: start(:)  ! The point to start writing at
                                 ! (one value for each dimension
                                 ! of the variable in the file)
INTEGER, INTENT(IN) :: count_in(:)  ! The number of points to write
                                    ! in each dimension of the variable


! Work variables
REAL, ALLOCATABLE :: local_values(:)  ! A 1d version of values

!-----------------------------------------------------------------------------

! To implement this routine, we just delegate to the 1d routine and reshape
! the result

ALLOCATE(local_values(PRODUCT(SHAPE(values))))

local_values(:) = RESHAPE(values, (/ PRODUCT(SHAPE(values)) /))

CALL file_ascii_write_var_1d(FILE, var_id, local_values, start, count_in)

DEALLOCATE(local_values)

RETURN

END SUBROUTINE file_ascii_write_var_5d


SUBROUTINE file_ascii_write_var_6d(FILE, var_id, values, start, count_in)

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
TYPE(file_ascii), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, DIMENSION(:,:,:,:,:,:), INTENT(IN) :: values  ! The values to write
INTEGER, INTENT(IN) :: start(:)  ! The point to start writing at
                                 ! (one value for each dimension
                                 ! of the variable in the file)
INTEGER, INTENT(IN) :: count_in(:)  ! The number of points to write
                                    ! in each dimension of the variable


! Work variables
REAL, ALLOCATABLE :: local_values(:)  ! A 1d version of values

!-----------------------------------------------------------------------------

! To implement this routine, we just delegate to the 1d routine and reshape
! the result

ALLOCATE(local_values(PRODUCT(SHAPE(values))))

local_values(:) = RESHAPE(values, (/ PRODUCT(SHAPE(values)) /))

CALL file_ascii_write_var_1d(FILE, var_id, local_values, start, count_in)

DEALLOCATE(local_values)

RETURN

END SUBROUTINE file_ascii_write_var_6d


SUBROUTINE file_ascii_write_var_7d(FILE, var_id, values, start, count_in)

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
TYPE(file_ascii), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, DIMENSION(:,:,:,:,:,:,:), INTENT(IN) :: values  ! The values to write
INTEGER, INTENT(IN) :: start(:)  ! The point to start writing at
                                 ! (one value for each dimension
                                 ! of the variable in the file)
INTEGER, INTENT(IN) :: count_in(:)  ! The number of points to write
                                    ! in each dimension of the variable


! Work variables
REAL, ALLOCATABLE :: local_values(:)  ! A 1d version of values

!-----------------------------------------------------------------------------

! To implement this routine, we just delegate to the 1d routine and reshape
! the result

ALLOCATE(local_values(PRODUCT(SHAPE(values))))

local_values(:) = RESHAPE(values, (/ PRODUCT(SHAPE(values)) /))

CALL file_ascii_write_var_1d(FILE, var_id, local_values, start, count_in)

DEALLOCATE(local_values)

RETURN

END SUBROUTINE file_ascii_write_var_7d
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ascii_close(FILE)

USE io_constants, ONLY: mode_write
USE errormessagelength_mod, ONLY: errormessagelength
USE file_ascii_generic_sync_mod, ONLY: file_sync
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
TYPE(file_ascii), INTENT(INOUT) :: FILE ! The file to close


! Work variables
INTEGER :: error ! Error code for any errors that occur
CHARACTER(LEN=errormessagelength) :: iomessage

!-----------------------------------------------------------------------------

CALL log_info("file_ascii_close",                                             &
              "Closing file " // TRIM(FILE%NAME))

! Flush the buffer to disk only if it has uncommitted writes
IF ( FILE%mode == mode_write ) THEN

  IF ( FILE%buffer_is_dirty ) THEN
    CALL file_ascii_flush_buffer(FILE)
  END IF

  CALL file_sync(FILE%UNIT)

END IF

! Deallocate the memory associated with the buffer
DEALLOCATE(FILE%buffer)
NULLIFY(FILE%buffer)

CLOSE(FILE%UNIT, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("file_ascii_close",                                          &
                 "Error closing file " // TRIM(FILE%NAME) //                  &
                 " (IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //        &
                 TRIM(iomessage) // ")")

RETURN

END SUBROUTINE file_ascii_close

! Inquiry routines
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ascii_introspect(FILE)

USE io_constants, ONLY: mode_read

USE string_utils_mod, ONLY: str_starts_with, str_replace, str_split
USE errormessagelength_mod, ONLY: errormessagelength

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
TYPE(file_ascii), INTENT(INOUT) :: FILE  ! The file to detect items in

! Work variables
CHARACTER(LEN=5000) :: line  ! The current line of the file
CHARACTER(LEN=errormessagelength) :: iomessage
! Information about the dimensions found in the file
INTEGER :: ndims
CHARACTER(LEN=max_sdf_name_len) :: dim_names(max_dim_file)
INTEGER :: dim_sizes(max_dim_file)
INTEGER :: dim_ids(max_dim_file)  ! The ids of the dimensions

! Information about the record dimension
LOGICAL :: has_record_dim
CHARACTER(LEN=max_sdf_name_len) :: record_dim_name

! Information about the variables found in the file
INTEGER :: nvars
CHARACTER(LEN=max_sdf_name_len) :: var_names(max_var_file)
INTEGER :: var_ndims(max_var_file)
CHARACTER(LEN=max_sdf_name_len) :: var_dim_names(max_var_file, max_dim_var)

! Used in calls to str_split
INTEGER :: nparts
CHARACTER(LEN=200) :: parts(2)
CHARACTER(LEN=200) :: dim_str

INTEGER :: dummy  ! Throwaway value

INTEGER :: var_dim_ids(max_dim_var)  ! Used when defining variables
                                     ! The ids of the dimensions found


INTEGER :: i,j,k  ! Index variables

INTEGER :: error  ! Error indicator


!-----------------------------------------------------------------------------

! We can only do this in read mode
IF ( FILE%mode /= mode_read )                                                 &
  CALL log_fatal("file_ascii_introspect",                                     &
                 "Can only introspect files in read mode")

!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
ndims          = 0
has_record_dim = .FALSE.
nvars          = 0

!-----------------------------------------------------------------------------
! Try to parse the header
! We can only do this for files that were created using this IO library and
! hence we know the format of the header...
!-----------------------------------------------------------------------------
! Start by rewinding the file to the start
REWIND(FILE%UNIT, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("file_ascii_introspect",                                     &
                 "Error rewinding file " //                                   &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

!-----------------------------------------------------------------------------
! Find the line indicating start of dimensions
!-----------------------------------------------------------------------------
DO
  READ(FILE%UNIT, "(A)", IOSTAT = error, IOMSG = iomessage) line
  ! Check for end of file
  IF ( error < 0 )                                                            &
    CALL log_fatal("file_ascii_introspect",                                   &
                   "File does not have expected format" //                    &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                    TRIM(iomessage) // ")")
  ! Report a general read error
  IF ( error > 0 )                                                            &
    CALL log_fatal("file_ascii_introspect",                                   &
                   "Error reading from file - likely bad file format " //     &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")

  ! We either keep reading until we get an error, or we exit when we find the
  ! dimensions line
  IF ( str_starts_with(line, "# Dimensions") ) EXIT
END DO

!-----------------------------------------------------------------------------
! Read information about the dimensions
!-----------------------------------------------------------------------------
DO
  READ(FILE%UNIT, "(A)", IOSTAT = error, IOMSG = iomessage) line
  ! Check for end of file
  IF ( error < 0 )                                                            &
    CALL log_fatal("file_ascii_introspect",                                   &
                   "File does not have expected format")
  ! Report a general read error
  IF ( error > 0 )                                                            &
    CALL log_fatal("file_ascii_introspect",                                   &
                   "Error reading from file - likely bad file format " //     &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")

  ! Chop the comment character off
  line = line(2:)

  ! If we reach an empty line, then we are done with dimensions
  IF ( LEN_TRIM(line) == 0 ) EXIT

  ! If we have a dimension that we can interpret, we should be able to split
  ! it into 2 parts about the '='
  CALL str_split(line, "=", nparts, parts)
  IF ( nparts /= 2 )                                                          &
  ! If we can't get two parts, then we have an unexpected format
        CALL log_fatal("file_ascii_introspect",                               &
                       "File does not have expected format")

  ! If we get to here, we have successfully found a dimension
  ! The dimension name is the bit before the '=' (parts(1))
  ! The dimension length is the bit after the '=' (parts(2))

  ! Check if it is the record dimension by looking for UNLIMITED as length
  IF ( INDEX(parts(2), "UNLIMITED") > 0 ) THEN
    has_record_dim = .TRUE.
    record_dim_name = ADJUSTL(parts(1))
    CYCLE
  END IF

  ! If we get to here then we have a proper dimension
  ndims = ndims + 1
  dim_names(ndims) = ADJUSTL(parts(1))
  ! Try to read an integer out of the size part using list directed IO
  READ(parts(2), *, IOSTAT = error) dim_sizes(ndims)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("file_ascii_introspect",                                   &
                   "File does not have expected format")
END DO


!-----------------------------------------------------------------------------
! Find the line indicating start of variables
!-----------------------------------------------------------------------------
DO
  READ(FILE%UNIT, "(A)", IOSTAT = error, IOMSG = iomessage) line
  ! Check for end of file
  IF ( error < 0 )                                                            &
    CALL log_fatal("file_ascii_introspect",                                   &
                   "File does not have expected format" //                    &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")
  ! Report a general read error
  IF ( error > 0 )                                                            &
    CALL log_fatal("file_ascii_introspect",                                   &
                   "Error reading from file - likely bad file format " //     &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")

  ! We either keep reading until we get an error, or we exit when we find the
  ! variables line
  IF ( str_starts_with(line, "# Variables") ) EXIT
END DO

!-----------------------------------------------------------------------------
! Read information about the variables
!-----------------------------------------------------------------------------
DO
  READ(FILE%UNIT, "(A)", IOSTAT = error, IOMSG = iomessage) line
  ! Check for end of file
  IF ( error < 0 )                                                            &
    CALL log_fatal("file_ascii_introspect",                                   &
                   "File does not have expected format" //                    &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")
  ! Report a general read error
  IF ( error > 0 )                                                            &
    CALL log_fatal("file_ascii_introspect",                                   &
                   "Error reading from file - likely bad file format " //     &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")

  ! If we've gone past the end of the header, we are done
  IF ( line(1:1) /= "#" ) EXIT

  ! Chop the comment character off
  line = line(2:)

  ! Ignore empty lines and lines with an '=' character in (attribute lines)
  IF ( LEN_TRIM(line) == 0 .OR. INDEX(line, "=") > 0 ) CYCLE

  ! If we have a variable that we can interpret, we should be able to split
  ! it into 2 parts about the '('
  CALL str_split(line, "(", nparts, parts)
  IF ( nparts /= 2 )                                                          &
  ! If we can't get two parts, then we have an unexpected format
        CALL log_fatal("file_ascii_introspect",                               &
                       "File does not have expected format")

  ! If we get to here, we have successfully found a variable
  nvars = nvars + 1
  !   parts(1) contains the variable name - all we need to do is left align it
  var_names(nvars) = ADJUSTL(parts(1))

  ! We want to get rid of the trailing ') ;' from parts(2) to get just the
  ! dimension names separated by commas
  dim_str = str_replace(parts(2), ") ;", " ")

  ! If there is no string left, then we have a variable with no non-record
  ! dimensions
  IF ( LEN_TRIM(dim_str) == 0 ) THEN
    var_ndims(nvars) = 0
    CYCLE
  END IF

  ! Now we can extract the dimensions by just splitting about commas
  CALL str_split(dim_str, ",", var_ndims(nvars), var_dim_names(nvars,:))
END DO

!-----------------------------------------------------------------------------
! We have finished reading the file, so rewind it for consistency
!-----------------------------------------------------------------------------
REWIND(FILE%UNIT, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("file_ascii_introspect",                                     &
                 "Error rewinding file " //                                   &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")


!-----------------------------------------------------------------------------
! Now we know what dimensions and variables are in the file, we can define
! them
!-----------------------------------------------------------------------------
! First, dimensions
DO i = 1,ndims
  dim_ids(i) = file_ascii_def_dim(FILE, dim_names(i), dim_sizes(i))
END DO
! And record dimension
IF ( has_record_dim )                                                         &
  dummy = file_ascii_def_record_dim(FILE, record_dim_name)

! Then define variables
DO i = 1,nvars
  ! First, gather up the dimension ids
  DO j = 1,var_ndims(i)
    ! If the variable has a dimension that isn't in the file, error out
    IF ( .NOT. ANY(dim_names == var_dim_names(i,j)) )                         &
      CALL log_fatal("file_ascii_introspect",                                 &
                     "Variable has a dimension not defined in file")

    ! Otherwise, get the corresponding dimension id
    DO k = 1,ndims
      IF ( var_dim_names(i,j) == dim_names(k) ) THEN
        var_dim_ids(j) = dim_ids(k)
        EXIT
      END IF
    END DO
  END DO

  ! Now we have the dimension ids for the variable, we can define it
  ! We throw away the id as we don't care what it is
  dummy = file_ascii_def_var(                                                 &
    FILE, var_names(i), var_dim_ids(1:var_ndims(i)), has_record_dim           &
  )
END DO

! Finally, we take the file out of define mode
CALL file_ascii_enddef(FILE)

RETURN

END SUBROUTINE file_ascii_introspect
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ascii_inquire_dim(FILE, dim_name, dim_id, dim_len, is_record_dim)

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
TYPE(file_ascii), INTENT(IN) :: FILE  ! The file to inspect
CHARACTER(LEN=*), INTENT(IN) :: dim_name  ! The name of the dimension to
                                          ! inquire about
INTEGER, INTENT(OUT) :: dim_id  ! The id of the dimension in the file
INTEGER, INTENT(OUT) :: dim_len  ! The length of the dimension in the file
LOGICAL, INTENT(OUT) :: is_record_dim  ! Indicates if the named dimension
                                       ! is the record dimension for the file

! Work variables
INTEGER :: i  ! Loop counter


!-----------------------------------------------------------------------------


! If we are still in define mode, we can't inspect
IF ( FILE%define_mode )                                                       &
  CALL log_fatal("file_ascii_inquire_dim",                                    &
                 "Cannot inquire file in define mode")

! Deal with the case that the record dimension has the specified name first
is_record_dim = FILE%has_record_dim .AND. ( FILE%record_dim_name == dim_name )
IF ( is_record_dim ) THEN
  RETURN
END IF

! Find the index in the dims arrays of the dimension with the given name
! This will be the returned dimension id
dim_id = -1
DO i = 1,FILE%ndims
  IF ( FILE%dim_names(i) == dim_name ) THEN
    dim_id = i
    EXIT
  END IF
END DO

! If we didn't find the dimension name, we have nothing more to do
IF ( dim_id < 1 ) RETURN

! Now we can easily get the size
dim_len = FILE%dim_sizes(dim_id)

RETURN

END SUBROUTINE file_ascii_inquire_dim
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ascii_inquire_var(FILE, var_name, var_id, ndims, dim_ids, is_record)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Given a variable name, returns its id and information about dimensions
!   If the returned id < 0, the varension doesn't exist
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ascii), INTENT(IN) :: FILE  ! The file to inspect
CHARACTER(LEN=*), INTENT(IN) :: var_name  ! The name of the variable to
                                          ! inquire about
INTEGER, INTENT(OUT) :: var_id  ! The id of the variable in the file
INTEGER, INTENT(OUT) :: ndims  ! The number of dimensions that
                                         ! the variable has
INTEGER, INTENT(OUT) :: dim_ids(:)  ! The dimension ids
LOGICAL, INTENT(OUT) :: is_record   ! Indicates if the variable
                                    ! uses the record dimension

! Work variables
INTEGER :: i  ! Loop counter


!-----------------------------------------------------------------------------


! If we are still in define mode, we can't inspect
IF ( FILE%define_mode )                                                       &
  CALL log_fatal("file_ascii_inquire_var",                                    &
                 "Cannot inquire file in define mode")

! Find the index in the vars arrays of the variable with the given name
! This will be the returned variable id
var_id = -1
DO i = 1,FILE%nvars
  IF ( FILE%var_names(i) == var_name ) THEN
    var_id = i
    EXIT
  END IF
END DO

! If we didn't find the variable name, we are done
IF ( var_id < 1 ) RETURN

! Otherwise, the var_id is the location in the array, and we can easily get
! the number of dimensions and their ids
ndims = FILE%var_ndims(var_id)
dim_ids(1:ndims) = FILE%var_dims(var_id,1:ndims)

! Indicate if the variable uses the record dimension
! In ASCII files, either all variables do or none do
is_record = FILE%has_record_dim

RETURN

END SUBROUTINE file_ascii_inquire_var

! Other utility routines
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ascii_fill_buffer(FILE)

USE io_constants, ONLY: mode_read
USE errormessagelength_mod, ONLY: errormessagelength
USE string_utils_mod, ONLY: str_starts_with

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   INTERNAL PROCEDURE - reads the next data from file, skipping any comment
!   lines, and fills the buffer with it
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ascii), INTENT(INOUT) :: FILE


! Work variables
CHARACTER(LEN=FILE%record_len) :: line
CHARACTER(LEN=errormessagelength) :: iomessage

LOGICAL :: is_comment  ! Indicates if the current line is a comment line
INTEGER :: i  ! Loop counter

INTEGER :: error  ! Error indicator


!-----------------------------------------------------------------------------

! We can't fill the buffer until enddef has been called
IF ( FILE%define_mode )                                                       &
  CALL log_fatal("file_ascii_fill_buffer",                                    &
                 "Cannot fill buffer - file is still in define mode")

IF ( FILE%mode /= mode_read )                                                 &
  CALL log_fatal("file_ascii_fill_buffer",                                    &
                 "Cannot read from file - file is not in read mode")

! Read lines from the file until we get a non-comment line
DO
  READ(FILE%UNIT, "(A)", IOSTAT = error, IOMSG = iomessage) line
  IF ( error /= 0 )                                                           &
    CALL log_fatal("file_ascii_fill_buffer",                                  &
                   "Error reading from file " //                              &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")

  ! Remove any leading spaces
  line = ADJUSTL(line)

  ! We skip over empty lines
  IF ( LEN_TRIM(line) == 0 ) CYCLE

  ! If the line does not start with a comment char, we have found a line of data
  ! and can exit the loop
  is_comment = .FALSE.
  DO i = 1,SIZE(comment_chars)
    IF ( str_starts_with(line, comment_chars(i)) ) THEN
      is_comment = .TRUE.
      EXIT
    END IF
  END DO
  IF ( .NOT. is_comment ) EXIT
END DO

! Fill the buffer from the line of data
READ(line, *) FILE%buffer(:)

RETURN

END SUBROUTINE file_ascii_fill_buffer
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ascii_flush_buffer(FILE)

USE io_constants, ONLY: mode_write, mdi
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   INTERNAL PROCEDURE - flushes the current buffer to file
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_ascii), INTENT(INOUT) :: FILE


! Work variables
INTEGER :: error  ! Error indicator
CHARACTER(LEN=errormessagelength) :: iomessage

!-----------------------------------------------------------------------------

! We can't flush the buffer until enddef has been called
IF ( FILE%define_mode )                                                       &
  CALL log_fatal("file_ascii_flush_buffer",                                   &
                 "Cannot flush buffer - file is still in define mode")

IF ( FILE%mode /= mode_write )                                                &
  CALL log_fatal("file_ascii_flush_buffer",                                   &
                 "Cannot write to file - file is not in write mode")

! Just write the current buffer to the file
WRITE(FILE%UNIT, out_format_str, IOSTAT = error, IOMSG = iomessage)           &
     FILE%buffer(:)
IF ( error /= 0 )                                                             &
  CALL log_fatal("file_ascii_flush_buffer",                                   &
                 "Error writing to file " //                                  &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

! Reset buffer values to missing data value
FILE%buffer(:) = mdi
! Buffer is clean again (i.e. has no writes that have not been committed to
! file)
FILE%buffer_is_dirty = .FALSE.

RETURN

END SUBROUTINE file_ascii_flush_buffer

END MODULE driver_ascii_mod
