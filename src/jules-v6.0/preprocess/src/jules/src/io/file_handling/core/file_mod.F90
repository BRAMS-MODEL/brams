! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE file_mod

USE driver_ascii_mod, ONLY: file_ascii
USE driver_ncdf_mod, ONLY: file_ncdf

USE logging_mod, ONLY: log_fatal

IMPLICIT NONE

! Constants for available drivers
INTEGER, PARAMETER :: driver_ascii = 1
INTEGER, PARAMETER :: driver_ncdf  = 2


! Type definitions
TYPE file_handle
  !-----------------------------------------------------------------------------
  ! This type and associated functions and subroutines are intended to provide
  ! a standard interface for different types of file
  !-----------------------------------------------------------------------------
  INTEGER :: driver ! The driver to use to process this file_handle object
                    ! This should be one of the drivers with a constant
                    ! defined above

  ! Only one of these representations is used for any given file_handle object
  TYPE(file_ascii) :: ascii ! The ASCII representation of the file
  TYPE(file_ncdf) :: ncdf   ! The NetCDF representation of the file

END TYPE file_handle


! Interface definition for overloads of def_attr for different types of
! variables
INTERFACE file_def_attr
MODULE PROCEDURE file_def_attr_real, file_def_attr_int, file_def_attr_char
END INTERFACE file_def_attr

! Interface definitions for the overloads of read/write var for variables
! of different ranks
INTERFACE file_read_var
MODULE PROCEDURE file_read_var_scalar, file_read_var_1d,                      &
                 file_read_var_2d,     file_read_var_3d,                      &
                 file_read_var_4d,     file_read_var_5d,                      &
                 file_read_var_6d,     file_read_var_7d
END INTERFACE file_read_var

INTERFACE file_write_var
MODULE PROCEDURE file_write_var_scalar, file_write_var_1d,                    &
                 file_write_var_2d,     file_write_var_3d,                    &
                 file_write_var_4d,     file_write_var_5d,                    &
                 file_write_var_6d,     file_write_var_7d
END INTERFACE file_write_var


!-----------------------------------------------------------------------------
! Visibility declarations
!-----------------------------------------------------------------------------
PRIVATE
PUBLIC                                                                        &
! Types
         file_handle,                                                         &
! Routines for opening and closing
         file_open, file_close,                                               &
! Routines for definitions
         file_def_dim, file_def_record_dim, file_def_var, file_def_attr,      &
         file_enddef,                                                         &
! Routines for seeking
         file_seek, file_advance,                                             &
! Routines for reading and writing
         file_read_var, file_write_var,                                       &
! Routines for inquiry
         file_introspect, file_inquire_dim, file_inquire_var


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


FUNCTION file_open(NAME, mode, comm, info) RESULT(FILE)

USE driver_ascii_mod, ONLY: extensions_ascii => extensions, file_ascii_open
USE driver_ncdf_mod, ONLY: extensions_ncdf => extensions, file_ncdf_open

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Opens a file and returns a file_handle object representing it
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

! Return type
TYPE(file_handle) :: FILE


! Work variables
INTEGER :: last_dot  ! The index of the last . (dot) in the file name
                     ! Used to extract the file extension
CHARACTER(LEN=15) :: extension  ! The extension for the file name


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that the MPI variables are either specified together or not at all
!-----------------------------------------------------------------------------
IF ( PRESENT(comm) .NEQV. PRESENT(info) )                                     &
  CALL log_fatal("file_open",                                                 &
                 "Only one of comm and info is present - either give a " //   &
                 "value for both MPI variables for parallel access or " //    &
                 "omit both for serial access")

!-----------------------------------------------------------------------------
! Get the extension for the file name (i.e. everything after the last dot)
!-----------------------------------------------------------------------------
last_dot = INDEX(NAME, '.', back = .TRUE.)
extension = NAME(last_dot+1:)

! Select a driver based on the file extension
IF ( ANY(extensions_ascii == extension) ) THEN
  ! Assign the correct driver
  FILE%driver = driver_ascii
  ! Initialise the ASCII representation of the file
  IF ( PRESENT(comm) ) THEN
    FILE%ascii = file_ascii_open(NAME, mode, comm, info)
  !DSM ELSE
    !DSM FILE%ascii = file_ascii_open(NAME, mode)
  END IF

ELSE IF ( ANY(extensions_ncdf == extension) ) THEN
  ! Assign the correct driver
  FILE%driver = driver_ncdf
  ! Initialise the NetCDF representation of the file
  IF ( PRESENT(comm) ) THEN
    FILE%ncdf = file_ncdf_open(NAME, mode, comm, info)
  ELSE
    if (NAME/='xxx.nc') FILE%ncdf = file_ncdf_open(NAME, mode)
  END IF

ELSE
  ! File type not recognised
  CALL log_fatal("file_open",                                                 &
                 "Unrecognised file extension for file " // TRIM(NAME) //     &
                 " - see docs for supported file types")
END IF

RETURN

END FUNCTION file_open
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_def_dim(FILE, dim_name, dim_len) RESULT(dim_id)

USE io_constants, ONLY: max_sdf_name_len

USE driver_ascii_mod, ONLY: file_ascii_def_dim
USE driver_ncdf_mod, ONLY: file_ncdf_def_dim

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
TYPE(file_handle), INTENT(INOUT) :: FILE
                                ! The file to define the dimension on
CHARACTER(LEN=*), INTENT(IN) :: dim_name
                                ! The name of the dimension
INTEGER, INTENT(IN) :: dim_len  ! The length of the dimension

! Return type
INTEGER :: dim_id               ! The dimension id


!-----------------------------------------------------------------------------

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  dim_id = file_ascii_def_dim(FILE%ascii, dim_name, dim_len)

CASE ( driver_ncdf )
  dim_id = file_ncdf_def_dim(FILE%ncdf, dim_name, dim_len)

CASE DEFAULT
  CALL log_fatal("file_def_dim",                                              &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END FUNCTION file_def_dim
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_def_record_dim(FILE, dim_name) RESULT(dim_id)

USE io_constants, ONLY: max_sdf_name_len

USE driver_ascii_mod, ONLY: file_ascii_def_record_dim
USE driver_ncdf_mod, ONLY: file_ncdf_def_record_dim

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
TYPE(file_handle), INTENT(INOUT) :: FILE
                                ! The file to define the dimension on
CHARACTER(LEN=*), INTENT(IN) :: dim_name
                                ! The name of the record dimension

! Return type
INTEGER :: dim_id               ! The dimension id


!-----------------------------------------------------------------------------

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  dim_id = file_ascii_def_record_dim(FILE%ascii, dim_name)

CASE ( driver_ncdf )
  dim_id = file_ncdf_def_record_dim(FILE%ncdf, dim_name)

CASE DEFAULT
  CALL log_fatal("file_def_record_dim",                                       &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END FUNCTION file_def_record_dim
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_def_var(FILE, var_name, dims, is_record) RESULT(var_id)

USE io_constants, ONLY: max_sdf_name_len

USE driver_ascii_mod, ONLY: file_ascii_def_var
USE driver_ncdf_mod, ONLY: file_ncdf_def_var

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
TYPE(file_handle), INTENT(INOUT) :: FILE
                                ! The file to define the variable in
CHARACTER(LEN=*), INTENT(IN) :: var_name
                                  ! The name of the variable
INTEGER, INTENT(IN) :: dims(:)    ! The ids of the NON-RECORD dimensions of
                                  ! the variable
LOGICAL, INTENT(IN) :: is_record  ! Indicates whether the variable uses the
                                  ! record dimension

! Return type
INTEGER :: var_id               ! The variable id


!-----------------------------------------------------------------------------

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  var_id = file_ascii_def_var(FILE%ascii, var_name, dims, is_record)

CASE ( driver_ncdf )
  var_id = file_ncdf_def_var(FILE%ncdf, var_name, dims, is_record)

CASE DEFAULT
  CALL log_fatal("file_def_var",                                              &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END FUNCTION file_def_var
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_def_attr_real(FILE, var_id, NAME, VALUE)

USE io_constants, ONLY: max_sdf_name_len

USE driver_ascii_mod, ONLY: file_ascii_def_attr_real
USE driver_ncdf_mod, ONLY: file_ncdf_def_attr_real

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
TYPE(file_handle), INTENT(INOUT) :: FILE
                                  ! The file to define the attribute in
INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                  ! attribute on
CHARACTER(LEN=*), INTENT(IN) :: NAME
                                  ! The name of the attribute
REAL, INTENT(IN) :: VALUE         ! The value of the attribute


!-----------------------------------------------------------------------------

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_def_attr_real(FILE%ascii, var_id, NAME, VALUE)

CASE ( driver_ncdf )
  CALL file_ncdf_def_attr_real(FILE%ncdf, var_id, NAME, VALUE)

CASE DEFAULT
  CALL log_fatal("file_def_attr_real",                                        &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_def_attr_real


SUBROUTINE file_def_attr_int(FILE, var_id, NAME, VALUE)

USE io_constants, ONLY: max_sdf_name_len

USE driver_ascii_mod, ONLY: file_ascii_def_attr_int
USE driver_ncdf_mod, ONLY: file_ncdf_def_attr_int

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
TYPE(file_handle), INTENT(INOUT) :: FILE
                                  ! The file to define the attribute in
INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                  ! attribute on
CHARACTER(LEN=*), INTENT(IN) :: NAME
                                  ! The name of the attribute
INTEGER, INTENT(IN) :: VALUE      ! The value of the attribute


!-----------------------------------------------------------------------------

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_def_attr_int(FILE%ascii, var_id, NAME, VALUE)

CASE ( driver_ncdf )
  CALL file_ncdf_def_attr_int(FILE%ncdf, var_id, NAME, VALUE)

CASE DEFAULT
  CALL log_fatal("file_def_attr_int",                                         &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_def_attr_int


SUBROUTINE file_def_attr_char(FILE, var_id, NAME, VALUE)

USE io_constants, ONLY: max_sdf_name_len

USE driver_ascii_mod, ONLY: file_ascii_def_attr_char
USE driver_ncdf_mod, ONLY: file_ncdf_def_attr_char

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
TYPE(file_handle), INTENT(INOUT) :: FILE
                                  ! The file to define the attribute in
INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                  ! attribute on
CHARACTER(LEN=*), INTENT(IN) :: NAME
                                  ! The name of the attribute
CHARACTER(LEN=*), INTENT(IN) :: VALUE
                                  ! The value of the attribute


!-----------------------------------------------------------------------------

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_def_attr_char(FILE%ascii, var_id, NAME, VALUE)

CASE ( driver_ncdf )
  CALL file_ncdf_def_attr_char(FILE%ncdf, var_id, NAME, VALUE)

CASE DEFAULT
  CALL log_fatal("file_def_attr_char",                                        &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_def_attr_char
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_enddef(FILE)

USE driver_ascii_mod, ONLY: file_ascii_enddef
USE driver_ncdf_mod, ONLY: file_ncdf_enddef

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
TYPE(file_handle), INTENT(INOUT) :: FILE


!-----------------------------------------------------------------------------

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_enddef(FILE%ascii)

CASE ( driver_ncdf )
  CALL file_ncdf_enddef(FILE%ncdf)

CASE DEFAULT
  CALL log_fatal("file_enddef",                                               &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_enddef
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_seek(FILE, record)

USE driver_ascii_mod, ONLY: file_ascii_seek
USE driver_ncdf_mod, ONLY: file_ncdf_seek

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Seeks the file to before a particular record (i.e. the next time a
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
TYPE(file_handle), INTENT(INOUT) :: FILE  ! The file to seek
INTEGER, INTENT(IN) :: record           ! The record number to seek to


!-----------------------------------------------------------------------------

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_seek(FILE%ascii, record)

CASE ( driver_ncdf )
  !DSM CALL file_ncdf_seek(FILE%ncdf, record)

CASE DEFAULT
  CALL log_fatal("file_seek",                                                 &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_seek
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_advance(FILE)

USE driver_ascii_mod, ONLY: file_ascii_advance
USE driver_ncdf_mod, ONLY: file_ncdf_advance

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
TYPE(file_handle), INTENT(INOUT) :: FILE  ! The file to seek


!-----------------------------------------------------------------------------

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_advance(FILE%ascii)

CASE ( driver_ncdf )
  CALL file_ncdf_advance(FILE%ncdf)

CASE DEFAULT
  CALL log_fatal("file_advance",                                              &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_advance
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_read_var_scalar(FILE, var_id, VALUE, start)

USE io_constants, ONLY: max_dim_var

USE driver_ascii_mod, ONLY: file_ascii_read_var
USE driver_ncdf_mod, ONLY: file_ncdf_read_var

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
TYPE(file_handle), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, INTENT(OUT) :: VALUE  ! The value read from file

! Optional arguments
INTEGER, INTENT(IN), OPTIONAL :: start(:)  ! The point to start reading from
                                           ! (one value for each NON-RECORD
                                           ! dimension of the variable)


! Work variables
INTEGER :: local_start(max_dim_var)

!-----------------------------------------------------------------------------
local_start(:) = 1
IF ( PRESENT(start) ) local_start(:SIZE(start)) = start(:)

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_read_var(FILE%ascii, var_id, VALUE, local_start)

CASE ( driver_ncdf )
  CALL file_ncdf_read_var(FILE%ncdf, var_id, VALUE, local_start)

CASE DEFAULT
  CALL log_fatal("file_read_var_scalar",                                      &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_read_var_scalar


SUBROUTINE file_read_var_1d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var

USE driver_ascii_mod, ONLY: file_ascii_read_var
USE driver_ncdf_mod, ONLY: file_ncdf_read_var

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
TYPE(file_handle), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, DIMENSION(:), INTENT(OUT) :: values  ! The values read from file

! Optional arguments
INTEGER, INTENT(IN), OPTIONAL :: start(:)  ! The point to start reading from
                                           ! (one value for each NON-RECORD
                                           ! dimension of the variable)
INTEGER, INTENT(IN), OPTIONAL :: counter(:)  ! The number of points to read
                                           ! in each dimension of the variable


! Work variables
INTEGER :: local_start(max_dim_var)
INTEGER :: local_count(max_dim_var)

!-----------------------------------------------------------------------------

! The default is to start at the 1st element of each dimension
local_start(:) = 1
IF ( PRESENT(start) ) local_start(:SIZE(start)) = start(:)

! The default counter is the size of the array to fill, with the rest of the
! dimensions being one
local_count(:) = 1
local_count(1) = SIZE(values)
IF ( PRESENT(counter) ) local_count(:SIZE(counter)) = counter(:)

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_read_var(                                                   &
    FILE%ascii, var_id, values, local_start, local_count                      &
  )

CASE ( driver_ncdf )
  CALL file_ncdf_read_var(                                                    &
    FILE%ncdf, var_id, values, local_start, local_count                       &
  )

CASE DEFAULT
  CALL log_fatal("file_read_var_1d",                                          &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_read_var_1d


SUBROUTINE file_read_var_2d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var

USE driver_ascii_mod, ONLY: file_ascii_read_var
USE driver_ncdf_mod, ONLY: file_ncdf_read_var

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
TYPE(file_handle), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, DIMENSION(:,:), INTENT(OUT) :: values  ! The values read from file

! Optional arguments
INTEGER, INTENT(IN), OPTIONAL :: start(:)  ! The point to start reading from
                                           ! (one value for each NON-RECORD
                                           ! dimension of the variable)
INTEGER, INTENT(IN), OPTIONAL :: counter(:)  ! The number of points to read
                                           ! in each dimension of the variable


! Work variables
INTEGER :: local_start(max_dim_var)
INTEGER :: local_count(max_dim_var)

!-----------------------------------------------------------------------------

! The default is to start at the 1st element of each dimension
local_start(:) = 1
IF ( PRESENT(start) ) local_start(:SIZE(start)) = start(:)

! The default counter is the size of the array to fill from the first
! available dimensions
local_count(:) = 1
local_count(1:2) = SHAPE(values)
IF ( PRESENT(counter) ) local_count(:SIZE(counter)) = counter(:)

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_read_var(                                                   &
    FILE%ascii, var_id, values, local_start, local_count                      &
  )

CASE ( driver_ncdf )
  CALL file_ncdf_read_var(                                                    &
    FILE%ncdf, var_id, values, local_start, local_count                       &
  )

CASE DEFAULT
  CALL log_fatal("file_read_var_2d",                                          &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_read_var_2d


SUBROUTINE file_read_var_3d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var

USE driver_ascii_mod, ONLY: file_ascii_read_var
USE driver_ncdf_mod, ONLY: file_ncdf_read_var

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
TYPE(file_handle), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, DIMENSION(:,:,:), INTENT(OUT) :: values  ! The values read from file

! Optional arguments
INTEGER, INTENT(IN), OPTIONAL :: start(:)  ! The point to start reading from
                                           ! (one value for each NON-RECORD
                                           ! dimension of the variable)
INTEGER, INTENT(IN), OPTIONAL :: counter(:)  ! The number of points to read
                                           ! in each dimension of the variable


! Work variables
INTEGER :: local_start(max_dim_var)
INTEGER :: local_count(max_dim_var)

!-----------------------------------------------------------------------------

! The default is to start at the 1st element of each dimension
local_start(:) = 1
IF ( PRESENT(start) ) local_start(:SIZE(start)) = start(:)

! The default counter is the size of the array to fill from the first
! available dimensions
local_count(:) = 1
local_count(1:3) = SHAPE(values)
IF ( PRESENT(counter) ) local_count(:SIZE(counter)) = counter(:)

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_read_var(                                                   &
    FILE%ascii, var_id, values, local_start, local_count                      &
  )

CASE ( driver_ncdf )
  CALL file_ncdf_read_var(                                                    &
    FILE%ncdf, var_id, values, local_start, local_count                       &
  )

CASE DEFAULT
  CALL log_fatal("file_read_var_3d",                                          &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_read_var_3d


SUBROUTINE file_read_var_4d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var

USE driver_ascii_mod, ONLY: file_ascii_read_var
USE driver_ncdf_mod, ONLY: file_ncdf_read_var

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
TYPE(file_handle), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: values  ! The values read from file

! Optional arguments
INTEGER, INTENT(IN), OPTIONAL :: start(:)  ! The point to start reading from
                                           ! (one value for each NON-RECORD
                                           ! dimension of the variable)
INTEGER, INTENT(IN), OPTIONAL :: counter(:)  ! The number of points to read
                                           ! in each dimension of the variable


! Work variables
INTEGER :: local_start(max_dim_var)
INTEGER :: local_count(max_dim_var)

!-----------------------------------------------------------------------------

! The default is to start at the 1st element of each dimension
local_start(:) = 1
IF ( PRESENT(start) ) local_start(:SIZE(start)) = start(:)

! The default counter is the size of the array to fill from the first
! available dimensions
local_count(:) = 1
local_count(1:4) = SHAPE(values)
IF ( PRESENT(counter) ) local_count(:SIZE(counter)) = counter(:)

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_read_var(                                                   &
    FILE%ascii, var_id, values, local_start, local_count                      &
  )

CASE ( driver_ncdf )
  CALL file_ncdf_read_var(                                                    &
    FILE%ncdf, var_id, values, local_start, local_count                       &
  )

CASE DEFAULT
  CALL log_fatal("file_read_var_4d",                                          &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_read_var_4d


SUBROUTINE file_read_var_5d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var

USE driver_ascii_mod, ONLY: file_ascii_read_var
USE driver_ncdf_mod, ONLY: file_ncdf_read_var

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
TYPE(file_handle), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: values  ! The values read from file

! Optional arguments
INTEGER, INTENT(IN), OPTIONAL :: start(:)  ! The point to start reading from
                                           ! (one value for each NON-RECORD
                                           ! dimension of the variable)
INTEGER, INTENT(IN), OPTIONAL :: counter(:)  ! The number of points to read
                                           ! in each dimension of the variable


! Work variables
INTEGER :: local_start(max_dim_var)
INTEGER :: local_count(max_dim_var)

!-----------------------------------------------------------------------------

! The default is to start at the 1st element of each dimension
local_start(:) = 1
IF ( PRESENT(start) ) local_start(:SIZE(start)) = start(:)

! The default counter is the size of the array to fill from the first
! available dimensions
local_count(:) = 1
local_count(1:5) = SHAPE(values)
IF ( PRESENT(counter) ) local_count(:SIZE(counter)) = counter(:)

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_read_var(                                                   &
    FILE%ascii, var_id, values, local_start, local_count                      &
  )

CASE ( driver_ncdf )
  CALL file_ncdf_read_var(                                                    &
    FILE%ncdf, var_id, values, local_start, local_count                       &
  )

CASE DEFAULT
  CALL log_fatal("file_read_var_5d",                                          &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_read_var_5d


SUBROUTINE file_read_var_6d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var

USE driver_ascii_mod, ONLY: file_ascii_read_var
USE driver_ncdf_mod, ONLY: file_ncdf_read_var

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
TYPE(file_handle), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, DIMENSION(:,:,:,:,:,:), INTENT(OUT) :: values  ! The values read from file

! Optional arguments
INTEGER, INTENT(IN), OPTIONAL :: start(:)  ! The point to start reading from
                                           ! (one value for each NON-RECORD
                                           ! dimension of the variable)
INTEGER, INTENT(IN), OPTIONAL :: counter(:)  ! The number of points to read
                                           ! in each dimension of the variable


! Work variables
INTEGER :: local_start(max_dim_var)
INTEGER :: local_count(max_dim_var)

!-----------------------------------------------------------------------------

! The default is to start at the 1st element of each dimension
local_start(:) = 1
IF ( PRESENT(start) ) local_start(:SIZE(start)) = start(:)

! The default counter is the size of the array to fill from the first
! available dimensions
local_count(:) = 1
local_count(1:6) = SHAPE(values)
IF ( PRESENT(counter) ) local_count(:SIZE(counter)) = counter(:)

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_read_var(                                                   &
    FILE%ascii, var_id, values, local_start, local_count                      &
  )

CASE ( driver_ncdf )
  CALL file_ncdf_read_var(                                                    &
    FILE%ncdf, var_id, values, local_start, local_count                       &
  )

CASE DEFAULT
  CALL log_fatal("file_read_var_6d",                                          &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_read_var_6d


SUBROUTINE file_read_var_7d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var

USE driver_ascii_mod, ONLY: file_ascii_read_var
USE driver_ncdf_mod, ONLY: file_ncdf_read_var

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
TYPE(file_handle), INTENT(INOUT) :: FILE  ! The file to read from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
REAL, DIMENSION(:,:,:,:,:,:,:), INTENT(OUT) :: values  ! The values read from file

! Optional arguments
INTEGER, INTENT(IN), OPTIONAL :: start(:)  ! The point to start reading from
                                           ! (one value for each NON-RECORD
                                           ! dimension of the variable)
INTEGER, INTENT(IN), OPTIONAL :: counter(:)  ! The number of points to read
                                           ! in each dimension of the variable


! Work variables
INTEGER :: local_start(max_dim_var)
INTEGER :: local_count(max_dim_var)

!-----------------------------------------------------------------------------

! The default is to start at the 1st element of each dimension
local_start(:) = 1
IF ( PRESENT(start) ) local_start(:SIZE(start)) = start(:)

! The default counter is the size of the array to fill from the first
! available dimensions
local_count(:) = 1
local_count(1:7) = SHAPE(values)
IF ( PRESENT(counter) ) local_count(:SIZE(counter)) = counter(:)

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_read_var(                                                   &
    FILE%ascii, var_id, values, local_start, local_count                      &
  )

CASE ( driver_ncdf )
  CALL file_ncdf_read_var(                                                    &
    FILE%ncdf, var_id, values, local_start, local_count                       &
  )

CASE DEFAULT
  CALL log_fatal("file_read_var_7d",                                          &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_read_var_7d

! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_write_var_scalar(FILE, var_id, VALUE, start)

USE io_constants, ONLY: max_dim_var

USE driver_ascii_mod, ONLY: file_ascii_write_var
USE driver_ncdf_mod, ONLY: file_ncdf_write_var

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
TYPE(file_handle), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, INTENT(IN) :: VALUE  ! The value to write

! Optional arguments
INTEGER, INTENT(IN), OPTIONAL :: start(:)  ! The point to start writing at
                                           ! (one value for each NON-RECORD
                                           ! dimension of the variable)


! Work variables
INTEGER :: local_start(max_dim_var)

!-----------------------------------------------------------------------------
local_start(:) = 1
IF ( PRESENT(start) ) local_start(:SIZE(start)) = start(:)

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_write_var(FILE%ascii, var_id, VALUE, local_start)

CASE ( driver_ncdf )
  CALL file_ncdf_write_var(FILE%ncdf, var_id, VALUE, local_start)

CASE DEFAULT
  CALL log_fatal("file_write_var_scalar",                                     &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_write_var_scalar


SUBROUTINE file_write_var_1d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var

USE driver_ascii_mod, ONLY: file_ascii_write_var
USE driver_ncdf_mod, ONLY: file_ncdf_write_var

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
TYPE(file_handle), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, DIMENSION(:), INTENT(IN) :: values  ! The values to write

! Optional arguments
INTEGER, INTENT(IN), OPTIONAL :: start(:)  ! The point to start writing at
                                           ! (one value for each NON-RECORD
                                           ! dimension of the variable)
INTEGER, INTENT(IN), OPTIONAL :: counter(:)  ! The number of points to write
                                           ! in each dimension of the variable


! Work variables
INTEGER :: local_start(max_dim_var)
INTEGER :: local_count(max_dim_var)

!-----------------------------------------------------------------------------

! The default is to start at the 1st element of each dimension
local_start(:) = 1
IF ( PRESENT(start) ) local_start(:SIZE(start)) = start(:)

! The default counter is the size of the array to write
local_count(:) = 1
local_count(1) = SIZE(values)
IF ( PRESENT(counter) ) local_count(:SIZE(counter)) = counter(:)

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_write_var(                                                  &
    FILE%ascii, var_id, values, local_start, local_count                      &
  )

CASE ( driver_ncdf )
  CALL file_ncdf_write_var(                                                   &
    FILE%ncdf, var_id, values, local_start, local_count                       &
  )

CASE DEFAULT
  CALL log_fatal("file_write_var_1d",                                         &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_write_var_1d


SUBROUTINE file_write_var_2d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var

USE driver_ascii_mod, ONLY: file_ascii_write_var
USE driver_ncdf_mod, ONLY: file_ncdf_write_var

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
TYPE(file_handle), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, DIMENSION(:,:), INTENT(IN) :: values  ! The values to write

! Optional arguments
INTEGER, INTENT(IN), OPTIONAL :: start(:)  ! The point to start writing at
                                           ! (one value for each NON-RECORD
                                           ! dimension of the variable)
INTEGER, INTENT(IN), OPTIONAL :: counter(:)  ! The number of points to write
                                           ! in each dimension of the variable


! Work variables
INTEGER :: local_start(max_dim_var)
INTEGER :: local_count(max_dim_var)

!-----------------------------------------------------------------------------

! The default is to start at the 1st element of each dimension
local_start(:) = 1
IF ( PRESENT(start) ) local_start(:SIZE(start)) = start(:)

! The default counter is the size of the array to write
local_count(:) = 1
local_count(1:2) = SHAPE(values)
IF ( PRESENT(counter) ) local_count(:SIZE(counter)) = counter(:)

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_write_var(                                                  &
    FILE%ascii, var_id, values, local_start, local_count                      &
  )

CASE ( driver_ncdf )
  CALL file_ncdf_write_var(                                                   &
    FILE%ncdf, var_id, values, local_start, local_count                       &
  )

CASE DEFAULT
  CALL log_fatal("file_write_var_2d",                                         &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_write_var_2d


SUBROUTINE file_write_var_3d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var

USE driver_ascii_mod, ONLY: file_ascii_write_var
USE driver_ncdf_mod, ONLY: file_ncdf_write_var

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
TYPE(file_handle), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, DIMENSION(:,:,:), INTENT(IN) :: values  ! The values to write

! Optional arguments
INTEGER, INTENT(IN), OPTIONAL :: start(:)  ! The point to start writing at
                                           ! (one value for each NON-RECORD
                                           ! dimension of the variable)
INTEGER, INTENT(IN), OPTIONAL :: counter(:)  ! The number of points to write
                                           ! in each dimension of the variable


! Work variables
INTEGER :: local_start(max_dim_var)
INTEGER :: local_count(max_dim_var)

!-----------------------------------------------------------------------------

! The default is to start at the 1st element of each dimension
local_start(:) = 1
IF ( PRESENT(start) ) local_start(:SIZE(start)) = start(:)

! The default counter is the size of the array to write
local_count(:) = 1
local_count(1:3) = SHAPE(values)
IF ( PRESENT(counter) ) local_count(:SIZE(counter)) = counter(:)

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_write_var(                                                  &
    FILE%ascii, var_id, values, local_start, local_count                      &
  )

CASE ( driver_ncdf )
  CALL file_ncdf_write_var(                                                   &
    FILE%ncdf, var_id, values, local_start, local_count                       &
  )

CASE DEFAULT
  CALL log_fatal("file_write_var_3d",                                         &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_write_var_3d


SUBROUTINE file_write_var_4d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var

USE driver_ascii_mod, ONLY: file_ascii_write_var
USE driver_ncdf_mod, ONLY: file_ncdf_write_var

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
TYPE(file_handle), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: values  ! The values to write

! Optional arguments
INTEGER, INTENT(IN), OPTIONAL :: start(:)  ! The point to start writing at
                                           ! (one value for each NON-RECORD
                                           ! dimension of the variable)
INTEGER, INTENT(IN), OPTIONAL :: counter(:)  ! The number of points to write
                                           ! in each dimension of the variable


! Work variables
INTEGER :: local_start(max_dim_var)
INTEGER :: local_count(max_dim_var)

!-----------------------------------------------------------------------------

! The default is to start at the 1st element of each dimension
local_start(:) = 1
IF ( PRESENT(start) ) local_start(:SIZE(start)) = start(:)

! The default counter is the size of the array to write
local_count(:) = 1
local_count(1:4) = SHAPE(values)
IF ( PRESENT(counter) ) local_count(:SIZE(counter)) = counter(:)

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_write_var(                                                  &
    FILE%ascii, var_id, values, local_start, local_count                      &
  )

CASE ( driver_ncdf )
  CALL file_ncdf_write_var(                                                   &
    FILE%ncdf, var_id, values, local_start, local_count                       &
  )

CASE DEFAULT
  CALL log_fatal("file_write_var_4d",                                         &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_write_var_4d


SUBROUTINE file_write_var_5d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var

USE driver_ascii_mod, ONLY: file_ascii_write_var
USE driver_ncdf_mod, ONLY: file_ncdf_write_var

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
TYPE(file_handle), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, DIMENSION(:,:,:,:,:), INTENT(IN) :: values  ! The values to write

! Optional arguments
INTEGER, INTENT(IN), OPTIONAL :: start(:)  ! The point to start writing at
                                           ! (one value for each NON-RECORD
                                           ! dimension of the variable)
INTEGER, INTENT(IN), OPTIONAL :: counter(:)  ! The number of points to write
                                           ! in each dimension of the variable


! Work variables
INTEGER :: local_start(max_dim_var)
INTEGER :: local_count(max_dim_var)

!-----------------------------------------------------------------------------

! The default is to start at the 1st element of each dimension
local_start(:) = 1
IF ( PRESENT(start) ) local_start(:SIZE(start)) = start(:)

! The default counter is the size of the array to write
local_count(:) = 1
local_count(1:5) = SHAPE(values)
IF ( PRESENT(counter) ) local_count(:SIZE(counter)) = counter(:)

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_write_var(                                                  &
    FILE%ascii, var_id, values, local_start, local_count                      &
  )

CASE ( driver_ncdf )
  CALL file_ncdf_write_var(                                                   &
    FILE%ncdf, var_id, values, local_start, local_count                       &
  )

CASE DEFAULT
  CALL log_fatal("file_write_var_5d",                                         &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_write_var_5d


SUBROUTINE file_write_var_6d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var

USE driver_ascii_mod, ONLY: file_ascii_write_var
USE driver_ncdf_mod, ONLY: file_ncdf_write_var

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
TYPE(file_handle), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, DIMENSION(:,:,:,:,:,:), INTENT(IN) :: values  ! The values to write

! Optional arguments
INTEGER, INTENT(IN), OPTIONAL :: start(:)  ! The point to start writing at
                                           ! (one value for each NON-RECORD
                                           ! dimension of the variable)
INTEGER, INTENT(IN), OPTIONAL :: counter(:)  ! The number of points to write
                                           ! in each dimension of the variable


! Work variables
INTEGER :: local_start(max_dim_var)
INTEGER :: local_count(max_dim_var)

!-----------------------------------------------------------------------------

! The default is to start at the 1st element of each dimension
local_start(:) = 1
IF ( PRESENT(start) ) local_start(:SIZE(start)) = start(:)

! The default counter is the size of the array to write
local_count(:) = 1
local_count(1:6) = SHAPE(values)
IF ( PRESENT(counter) ) local_count(:SIZE(counter)) = counter(:)

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_write_var(                                                  &
    FILE%ascii, var_id, values, local_start, local_count                      &
  )

CASE ( driver_ncdf )
  CALL file_ncdf_write_var(                                                   &
    FILE%ncdf, var_id, values, local_start, local_count                       &
  )

CASE DEFAULT
  CALL log_fatal("file_write_var_6d",                                         &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_write_var_6d


SUBROUTINE file_write_var_7d(FILE, var_id, values, start, counter)

USE io_constants, ONLY: max_dim_var

USE driver_ascii_mod, ONLY: file_ascii_write_var
USE driver_ncdf_mod, ONLY: file_ncdf_write_var

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
TYPE(file_handle), INTENT(INOUT) :: FILE  ! The file to write to
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
REAL, DIMENSION(:,:,:,:,:,:,:), INTENT(IN) :: values  ! The values to write

! Optional arguments
INTEGER, INTENT(IN), OPTIONAL :: start(:)  ! The point to start writing at
                                           ! (one value for each NON-RECORD
                                           ! dimension of the variable)
INTEGER, INTENT(IN), OPTIONAL :: counter(:)  ! The number of points to write
                                           ! in each dimension of the variable


! Work variables
INTEGER :: local_start(max_dim_var)
INTEGER :: local_count(max_dim_var)

!-----------------------------------------------------------------------------

! The default is to start at the 1st element of each dimension
local_start(:) = 1
IF ( PRESENT(start) ) local_start(:SIZE(start)) = start(:)

! The default counter is the size of the array to write
local_count(:) = 1
local_count(1:7) = SHAPE(values)
IF ( PRESENT(counter) ) local_count(:SIZE(counter)) = counter(:)

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_write_var(                                                  &
    FILE%ascii, var_id, values, local_start, local_count                      &
  )

CASE ( driver_ncdf )
  CALL file_ncdf_write_var(                                                   &
    FILE%ncdf, var_id, values, local_start, local_count                       &
  )

CASE DEFAULT
  CALL log_fatal("file_write_var_7d",                                         &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_write_var_7d
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_close(FILE)

USE driver_ascii_mod, ONLY: file_ascii_close
USE driver_ncdf_mod, ONLY: file_ncdf_close

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
TYPE(file_handle), INTENT(INOUT) :: FILE ! The file to close


!-----------------------------------------------------------------------------

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_close(FILE%ascii)

CASE ( driver_ncdf )
  CALL file_ncdf_close(FILE%ncdf)

CASE DEFAULT
  CALL log_fatal("file_close",                                                &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_close

! Inquiry routines
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_inquire_dim(FILE, dim_name, dim_id, dim_len, is_record_dim)

USE driver_ascii_mod, ONLY: file_ascii_inquire_dim
USE driver_ncdf_mod, ONLY: file_ncdf_inquire_dim

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
TYPE(file_handle), INTENT(IN) :: FILE  ! The file to inspect
CHARACTER(LEN=*), INTENT(IN) :: dim_name  ! The name of the dimension to
                                          ! inquire about
INTEGER, INTENT(OUT) :: dim_id  ! The id of the dimension in the file
INTEGER, INTENT(OUT) :: dim_len  ! The length of the dimension in the file
                                 ! This value is unset if the dimension is
                                 ! a record dimension
LOGICAL, INTENT(OUT) :: is_record_dim  ! Indicates if the named dimension
                                       ! is the record dimension for the file


!-----------------------------------------------------------------------------

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_inquire_dim(                                                &
    FILE%ascii, dim_name, dim_id, dim_len, is_record_dim                      &
  )

CASE ( driver_ncdf )
  CALL file_ncdf_inquire_dim(                                                 &
    FILE%ncdf, dim_name, dim_id, dim_len, is_record_dim                       &
  )

CASE DEFAULT
  CALL log_fatal("file_inquire_dim",                                          &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_inquire_dim
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_inquire_var(FILE, var_name, var_id, ndims, dim_ids, is_record)

USE io_constants, ONLY: max_dim_file, max_file_name_len

USE driver_ascii_mod, ONLY: file_ascii_inquire_var
USE driver_ncdf_mod, ONLY: file_ncdf_inquire_var

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
TYPE(file_handle), INTENT(IN) :: FILE  ! The file to inspect
CHARACTER(LEN=*), INTENT(IN) :: var_name  ! The name of the variable to
                                          ! inquire about
INTEGER, INTENT(OUT) :: var_id  ! The id of the variable in the file
INTEGER, INTENT(OUT), OPTIONAL :: ndims  ! The number of dimensions that
                                         ! the variable has
INTEGER, INTENT(OUT), OPTIONAL :: dim_ids(:)  ! The dimension ids
LOGICAL, INTENT(OUT), OPTIONAL :: is_record  ! Indicates if the variable
                                             ! uses the record dimension


! Work variables
! Local versions of the optional arguments to pass to the underlying
! implementations
INTEGER :: ndims_local
INTEGER :: dim_ids_local(max_dim_file)
LOGICAL :: is_record_local
! Other local variables.
CHARACTER(LEN=max_file_name_len) :: file_name

!-----------------------------------------------------------------------------

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  file_name = FILE % ascii % NAME
  CALL file_ascii_inquire_var(                                                &
    FILE%ascii, var_name, var_id,                                             &
    ndims_local, dim_ids_local, is_record_local                               &
  )

CASE ( driver_ncdf )
  file_name = FILE % ncdf % NAME
  CALL file_ncdf_inquire_var(                                                 &
    FILE%ncdf, var_name, var_id,                                              &
    ndims_local, dim_ids_local, is_record_local                               &
  )

CASE DEFAULT
  CALL log_fatal("file_inquire_var",                                          &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

! Check that the variable was found.
IF ( var_id < 1 ) CALL log_fatal( "file_inquire_var",                         &
      "Variable " // TRIM(var_name) // " not in file " // TRIM(file_name) )

! Copy the returned values into the optional arguments if present
IF ( PRESENT(ndims) ) ndims = ndims_local
IF ( PRESENT(dim_ids) ) dim_ids(1:ndims) = dim_ids_local(1:ndims)
IF ( PRESENT(is_record) ) is_record = is_record_local

RETURN

END SUBROUTINE file_inquire_var
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_introspect(FILE)

USE driver_ascii_mod, ONLY: file_ascii_introspect
USE driver_ncdf_mod, ONLY: file_ncdf_introspect

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Given an open file_handle in read mode, try to detect the dimensions and
!   variables in the file and define them on the file_handle
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_handle), INTENT(INOUT) :: FILE  ! The file to detect items in


!-----------------------------------------------------------------------------

SELECT CASE ( FILE%driver )
CASE ( driver_ascii )
  CALL file_ascii_introspect(FILE%ascii)

CASE ( driver_ncdf )
  CALL file_ncdf_introspect(FILE%ncdf)

CASE DEFAULT
  CALL log_fatal("file_introspect",                                           &
                 "Unrecognised driver - see docs for available drivers")
END SELECT

RETURN

END SUBROUTINE file_introspect

END MODULE file_mod
