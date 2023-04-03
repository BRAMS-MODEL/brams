! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE file_gridded_mod

USE io_constants, ONLY: max_dim_file, max_var_file, max_dim_var

USE grid_utils_mod, ONLY: grid_info

USE file_mod, ONLY: file_handle

USE logging_mod, ONLY: log_fatal

IMPLICIT NONE


! Type definitions
TYPE var_gridded
  !-----------------------------------------------------------------------------
  ! INTERNAL TYPE
  !
  ! This type contains information about the variables in a gridded file
  !-----------------------------------------------------------------------------
  INTEGER :: id  ! The id of the variable in the underlying file

  INTEGER, POINTER :: lev_sizes(:) => NULL() ! The sizes of the levels dimensions
END TYPE var_gridded


TYPE file_gridded

  TYPE(file_handle) :: fh  ! Handle to the underlying file

  !-----------------------------------------------------------------------------
  ! Information about the grid to use
  !-----------------------------------------------------------------------------
  TYPE(grid_info) :: grid  ! The grid that we are using

  INTEGER :: grid_dim_id = -1  ! Dimension id of the single grid dimension
                               ! if using a 1d grid
  INTEGER :: grid_x_dim_id = -1  ! Dimension ids of the x and y dimensions
  INTEGER :: grid_y_dim_id = -1  ! if using a 2d grid

  !-----------------------------------------------------------------------------
  ! Information about the 'levels' dimensions that have been defined
  ! (i.e. dimensions that define the number of vertical levels a variable has)
  ! Each variable can have at most one 'levels' dimension
  !-----------------------------------------------------------------------------
  INTEGER :: ndims = 0  ! The number of 'levels' dimensions defined
  INTEGER :: dim_ids(max_dim_file)  ! The ids of the 'levels' dimensions
                                    ! in the file
  INTEGER :: dim_sizes(max_dim_file)  ! The sizes of the 'levels' dimensions

  !-----------------------------------------------------------------------------
  ! Information about the number of levels that variable in the file have
  !-----------------------------------------------------------------------------
  INTEGER :: nvars = 0  ! The number of variables in the file
  TYPE(var_gridded) :: vars(max_var_file)  ! The variables in the file

END TYPE file_gridded


! Overloads for file_gridded_def_attr
INTERFACE file_gridded_def_attr
MODULE PROCEDURE file_gridded_def_attr_real, file_gridded_def_attr_int,       &
                 file_gridded_def_attr_char
END INTERFACE file_gridded_def_attr


!-----------------------------------------------------------------------------
! Visibility declarations
!-----------------------------------------------------------------------------
PRIVATE
PUBLIC                                                                        &
! Types
         file_gridded,                                                        &
! Routines for opening and closing
         file_gridded_open, file_gridded_close,                               &
! Routines for definitions
         file_gridded_def_grid, file_gridded_def_dim,                         &
         file_gridded_def_record_dim, file_gridded_def_var,                   &
         file_gridded_def_attr, file_gridded_enddef,                          &
! Routines for seeking
         file_gridded_seek, file_gridded_advance,                             &
! Routines for reading and writing
         file_gridded_read_var, file_gridded_write_var


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


FUNCTION file_gridded_open(NAME, mode, comm, info) RESULT(FILE)

USE file_mod, ONLY: file_open

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Opens a gridded file and returns a file_gridded object representing it
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
TYPE(file_gridded) :: FILE


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that the MPI variables are either specified together or not at all
!-----------------------------------------------------------------------------
IF ( PRESENT(comm) .NEQV. PRESENT(info) )                                     &
  CALL log_fatal("file_gridded_open",                                         &
                 "Only one of comm and info is present - either give a " //   &
                 "value for both MPI variables for parallel access or " //    &
                 "omit both for serial access")

! All we need to do is open the underlying file
IF ( PRESENT(comm) ) THEN
  FILE%fh = file_open(NAME, mode, comm, info)
ELSE
  FILE%fh = file_open(NAME, mode)
END IF

RETURN

END FUNCTION file_gridded_open
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_gridded_def_grid(FILE, grid)

USE file_mod, ONLY: file_def_dim

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
TYPE(file_gridded), INTENT(INOUT) :: FILE
    ! The file to define the grid on

TYPE(grid_info), INTENT(IN) :: grid  ! The grid to define


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Define the dimensions for the grid on the underlying file and store the
! dimension ids for later use
!-----------------------------------------------------------------------------
FILE%grid = grid

IF ( grid%is_1d ) THEN
  FILE%grid_dim_id  = file_def_dim(FILE%fh, grid%dim_name, grid%nx)
ELSE
  FILE%grid_x_dim_id = file_def_dim(FILE%fh, grid%x_name, grid%nx)
  FILE%grid_y_dim_id = file_def_dim(FILE%fh, grid%y_name, grid%ny)
END IF

RETURN

END SUBROUTINE file_gridded_def_grid
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_gridded_def_dim(FILE, dim_name, dim_len) RESULT(dim_id)

USE file_mod, ONLY: file_def_dim

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Defines a non-grid (levels) dimension on the given file
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_gridded), INTENT(INOUT) :: FILE
    ! The file to define the dimension on

CHARACTER(LEN=*), INTENT(IN) :: dim_name
    ! The name of the dimension
INTEGER, INTENT(IN) :: dim_len
    ! The size of the dimension

! Return type
INTEGER :: dim_id  ! The dimension id


!-----------------------------------------------------------------------------


! Check if we already have the maximum number of dimensions
IF ( FILE%ndims >= max_dim_file )                                             &
  CALL log_fatal("file_gridded_def_dim",                                      &
                 "Too many dimensions in file - try increasing max_dim_file")

! Define the dimension on the underlying file and store its id for later use
FILE%ndims = FILE%ndims + 1

! The returned dimension id is just the index in the dim_ids array
dim_id = FILE%ndims

FILE%dim_ids(dim_id)   = file_def_dim(FILE%fh, dim_name, dim_len)
FILE%dim_sizes(dim_id) = dim_len

RETURN

END FUNCTION file_gridded_def_dim
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_gridded_def_record_dim(FILE, dim_name) RESULT(dim_id)

USE file_mod, ONLY: file_def_record_dim

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Defines a record dimension on the given file, returning the dimension id
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_gridded), INTENT(INOUT) :: FILE
                                ! The file to define the dimension on
CHARACTER(LEN=*), INTENT(IN) :: dim_name
                                ! The name of the record dimension

! Return type
INTEGER :: dim_id               ! The dimension id


!-----------------------------------------------------------------------------

! Just defer to the underlying file
dim_id = file_def_record_dim(FILE%fh, dim_name)

RETURN

END FUNCTION file_gridded_def_record_dim
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_gridded_def_var(FILE, var_name, levels_dims, is_record)         &
                                                                RESULT(var_id)

USE file_mod, ONLY: file_def_var

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
TYPE(file_gridded), INTENT(INOUT) :: FILE
                                ! The file to define the variable in
CHARACTER(LEN=*), INTENT(IN) :: var_name
                                  ! The name of the variable
INTEGER, INTENT(IN), OPTIONAL :: levels_dims(:)
                                  ! The ids of the dimensions to use for the
                                  ! vertical levels of the variable, if
                                  ! required
                                  ! If not given or an array of 0 length is
                                  ! given, no vertical levels are used
LOGICAL, INTENT(IN) :: is_record  ! Indicates whether the variable uses the
                                  ! record dimension

! Return type
INTEGER :: var_id               ! The variable id


! Work variables
INTEGER :: ndims_levs  ! The number of levels dimensions that the variable has
INTEGER :: ndims       ! The total number of dimensions that the variable has
INTEGER, ALLOCATABLE :: dim_ids(:)  ! The ids in the underlying file of the
                                    ! dimensions of the variable
                                    ! The variable has:
                                    !   * 1 or 2 dimensions for grid
                                    !   * 0 or more levels

INTEGER :: i  ! Loop index


!-----------------------------------------------------------------------------


! Check if we already have the maximum number of variables
IF ( FILE%nvars >= max_var_file )                                             &
  CALL log_fatal("file_gridded_def_var",                                      &
                 "Too many variables in file - try increasing max_var_file")

!-----------------------------------------------------------------------------
! Construct the array of dimension ids for the variable
!-----------------------------------------------------------------------------
! How many levels dimensions do we have
ndims_levs = 0
IF ( PRESENT(levels_dims) ) ndims_levs = SIZE(levels_dims)

! Construct the array of dimension ids
! Grid dimensions come first
IF ( FILE%grid%is_1d ) THEN
  ALLOCATE(dim_ids(1 + ndims_levs))
  ndims = 1
  dim_ids(1) = FILE%grid_dim_id
ELSE
  ALLOCATE(dim_ids(2 + ndims_levs))
  ndims = 2
  dim_ids(1:2) = (/ FILE%grid_x_dim_id, FILE%grid_y_dim_id /)
END IF

DO i = 1,ndims_levs
  ndims = ndims + 1
  ! Each levels dim is an index in file%dim_ids, as returned by file_gridded_def_dim
  dim_ids(ndims) = FILE%dim_ids(levels_dims(i))
END DO

!-----------------------------------------------------------------------------
! Define the variable in the underlying file and store details for future use
!-----------------------------------------------------------------------------
FILE%nvars = FILE%nvars + 1

   ! The returned variable id is just the index in the var_ids array of the file
   var_id = FILE%nvars
!print*,'aaaa>>>',var_id,trim(var_name)
       
!DSM{
   if (trim(var_name)=='sw_downB') then
      FILE%vars(var_id)%id=1
   elseif (trim(var_name)=='lw_downB') then
      FILE%vars(var_id)%id=2
   elseif (trim(var_name)=='diff_radB') then
      FILE%vars(var_id)%id=3
   elseif (trim(var_name)=='precipB') then
      FILE%vars(var_id)%id=4
   elseif (trim(var_name)=='tB') then
      FILE%vars(var_id)%id=5
   elseif (trim(var_name)=='uB') then
      FILE%vars(var_id)%id=6
   elseif (trim(var_name)=='vB') then
      FILE%vars(var_id)%id=7
   elseif (trim(var_name)=='pstarB') then
      FILE%vars(var_id)%id=8
   elseif (trim(var_name)=='qB') then
      FILE%vars(var_id)%id=9
   else
      FILE%vars(var_id)%id = file_def_var(                                          &
        FILE%fh, var_name, dim_ids(1:ndims), is_record                              &
      )
   endif
!DSM}

!print*,'AA=>>',trim(var_name),FILE%vars(var_id)%id

! Note that we can do this even when ndims_levs = 0 - we just get a zero-sized
! array
! This actually simplifies things in file_gridded_read/write_var as well
ALLOCATE(FILE%vars(var_id)%lev_sizes(ndims_levs))
DO i = 1,ndims_levs
  FILE%vars(var_id)%lev_sizes(i) = FILE%dim_sizes(levels_dims(i))
END DO

DEALLOCATE(dim_ids)

RETURN

END FUNCTION file_gridded_def_var
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_gridded_def_attr_real(FILE, var_id, NAME, VALUE)

USE io_constants, ONLY: max_sdf_name_len, attr_global

USE file_mod, ONLY: file_def_attr

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
TYPE(file_gridded), INTENT(INOUT) :: FILE
                                  ! The file to define the attribute in
INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                  ! attribute on
CHARACTER(LEN=*), INTENT(IN) :: NAME
                                  ! The name of the attribute
REAL, INTENT(IN) :: VALUE         ! The value of the attribute

! Work variables
INTEGER :: var_id_local


!-----------------------------------------------------------------------------

! Look up the variable id in the underlying file, unless a global attribute
! has been requested
! var_id is an index in the var_ids array of the file_gridded object, as
! returned by file_gridded_def_var
var_id_local = var_id
IF ( var_id /= attr_global ) var_id_local = FILE%vars(var_id)%id

! Now just defer to the underlying file
CALL file_def_attr(FILE%fh, var_id_local, NAME, VALUE)

RETURN

END SUBROUTINE file_gridded_def_attr_real


SUBROUTINE file_gridded_def_attr_int(FILE, var_id, NAME, VALUE)

USE io_constants, ONLY: max_sdf_name_len, attr_global

USE file_mod, ONLY: file_def_attr

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
TYPE(file_gridded), INTENT(INOUT) :: FILE
                                  ! The file to define the attribute in
INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                  ! attribute on
CHARACTER(LEN=*), INTENT(IN) :: NAME
                                  ! The name of the attribute
INTEGER, INTENT(IN) :: VALUE      ! The value of the attribute

! Work variables
INTEGER :: var_id_local


!-----------------------------------------------------------------------------

! Look up the variable id in the underlying file, unless a global attribute
! has been requested
! var_id is an index in the var_ids array of the file_gridded object, as
! returned by file_gridded_def_var
var_id_local = var_id
IF ( var_id /= attr_global ) var_id_local = FILE%vars(var_id)%id

! Now just defer to the underlying file
CALL file_def_attr(FILE%fh, var_id_local, NAME, VALUE)

RETURN

END SUBROUTINE file_gridded_def_attr_int


SUBROUTINE file_gridded_def_attr_char(FILE, var_id, NAME, VALUE)

USE io_constants, ONLY: max_sdf_name_len, attr_global

USE file_mod, ONLY: file_def_attr

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
TYPE(file_gridded), INTENT(INOUT) :: FILE
                                  ! The file to define the attribute in
INTEGER, INTENT(IN) :: var_id     ! The id of the variable to define
                                  ! attribute on
CHARACTER(LEN=*), INTENT(IN) :: NAME
                                  ! The name of the attribute
CHARACTER(LEN=*), INTENT(IN) :: VALUE
                                  ! The value of the attribute

! Work variables
INTEGER :: var_id_local


!-----------------------------------------------------------------------------

! Look up the variable id in the underlying file, unless a global attribute
! has been requested
! var_id is an index in the var_ids array of the file_gridded object, as
! returned by file_gridded_def_var
var_id_local = var_id
IF ( var_id /= attr_global ) var_id_local = FILE%vars(var_id)%id

! Now just defer to the underlying file
CALL file_def_attr(FILE%fh, var_id_local, NAME, VALUE)

RETURN

END SUBROUTINE file_gridded_def_attr_char
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_gridded_enddef(FILE)

USE file_mod, ONLY: file_enddef

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
TYPE(file_gridded), INTENT(INOUT) :: FILE


!-----------------------------------------------------------------------------

! We just need to defer to the underlying file
CALL file_enddef(FILE%fh)

RETURN

END SUBROUTINE file_gridded_enddef
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_gridded_seek(FILE, record)

USE file_mod, ONLY: file_seek

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
TYPE(file_gridded), INTENT(INOUT) :: FILE  ! The file to seek
INTEGER, INTENT(IN) :: record              ! The record number to seek to


!-----------------------------------------------------------------------------

! We just need to defer to the underlying file
CALL file_seek(FILE%fh, record)

RETURN

END SUBROUTINE file_gridded_seek
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_gridded_advance(FILE)

USE file_mod, ONLY: file_advance

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
TYPE(file_gridded), INTENT(INOUT) :: FILE  ! The file to seek


!-----------------------------------------------------------------------------

! Just defer to the underlying file
CALL file_advance(FILE%fh)

RETURN

END SUBROUTINE file_gridded_advance
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_gridded_read_var(FILE, var_id, extract_subgrid, subgrid) RESULT(cube)

USE grid_utils_mod, ONLY: subgrid_info, operator( /= ), subgrid_extract

USE file_mod, ONLY: file_read_var

USE data_cube_mod, ONLY: data_cube, cube_create

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
TYPE(file_gridded), INTENT(INOUT) :: FILE
                               ! The file to read data from
INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
LOGICAL, INTENT(IN) :: extract_subgrid
                               ! T - extract a subgrid to return from the
                               !     full grid of the file
                               ! F - return the full grid
TYPE(subgrid_info), OPTIONAL, INTENT(IN) :: subgrid  ! The subgrid to extract


! Return type
TYPE(data_cube) :: cube  ! The data cube read from file


! Work variables
TYPE(var_gridded) :: var  ! The variable we are reading from

REAL, ALLOCATABLE :: values_3d(:,:,:)  ! 3D version of the cube values to
                                       ! simplify calculations using subgrids

REAL, ALLOCATABLE :: DATA(:,:,:)  ! Used to hold the full grid of data if
                                  ! extracting a subgrid in memory
                                  ! Uses a combined levels dimension to
                                  ! simplify calculations

INTEGER :: i, j, k, x, y  ! Index variables


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check if a subgrid has been requested and that the parent grid is the same
! as the grid the file is on
!-----------------------------------------------------------------------------
IF ( extract_subgrid ) THEN
  IF ( .NOT. PRESENT(subgrid) )                                               &
    CALL log_fatal("file_gridded_read_var",                                   &
                   "Subgrid extraction has been requested but no " //         &
                   "subgrid has been given")

  ! Check that the subgrid's parent grid matches the file's grid
  IF ( subgrid%parent /= FILE%grid )                                          &
    CALL log_fatal("file_gridded_read_var",                                   &
                   "Subgrid is inconsistent with grid in file")
END IF  ! extract_subgrid

!-----------------------------------------------------------------------------
! Extract information on the variable being read into a local variable for
! convenience
!-----------------------------------------------------------------------------
var = FILE%vars(var_id)

!-----------------------------------------------------------------------------
! Create a cube of the correct size
!-----------------------------------------------------------------------------
IF ( extract_subgrid ) THEN
  cube = cube_create((/ subgrid%nx, subgrid%ny, var%lev_sizes /))
ELSE
  cube = cube_create((/ FILE%grid%nx, FILE%grid%ny, var%lev_sizes /))
END IF

!-----------------------------------------------------------------------------
! Actually read the requested data
!
! The decision process for what to actually read from file is as follows:
!   * If no subgrid has been specified, read everything (obviously!)
!
!   * If a subgrid has been specified using a region, then read only that
!     region from file (easy to specify with start and count)
!
!   * If a subgrid has been specified using points, then one of two things
!     happens:
!       1. If only a small proportion of the full grid is to be extracted, we
!          read the values from file one at a time
!       2. If a significant proportion of the full grid is to extracted, we
!          read all the values and perform the extraction in memory
!
! This process provides the best trade-off between memory usage and time
! spent in I/O
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! If we are not extracting a subgrid, just read all the data into the cube
! and return
!-----------------------------------------------------------------------------
IF ( .NOT. extract_subgrid ) THEN
  IF ( FILE%grid%is_1d ) THEN
    CALL file_read_var(                                                       &
      FILE%fh, var%id, cube%values,                                           &
    ! We want to read the whole of the grid and levels dimensions, so set start
    ! and count accordingly
    ! For a 1D grid, we only need to consider the x axis of the grid
            (/ 1,            (1, i = 1,SIZE(var%lev_sizes)) /),               &
            (/ FILE%grid%nx, var%lev_sizes /)                                 &
          )
  ELSE

!JWA    CALL file_read_var(                                                       &
!JWA      FILE%fh, var%id, cube%values,                                           &
   ! We want to read the whole of the grid and levels dimensions, so set start
   ! and count accordingly
   ! For a 2D grid, we only need to consider both axes of the grid
!JWA            (/ 1,            1,            (1, i = 1,SIZE(var%lev_sizes)) /), &
!JWA            (/ FILE%grid%nx, FILE%grid%ny, var%lev_sizes /)                   &
!JWA          )
!Print*,'FILE%grid%nx=',FILE%grid%nx
!Print*,'FILE%grid%ny=',FILE%grid%ny
!Print*,'var%lev_sizes=',var%lev_sizes
!print*, 'size lev_sizes=', SIZE(var%lev_sizes)
      CALL read_drive_BRAMS(var%id,cube%values) !DSM,&
         !DSM (/ 1,            1,            (1, i = 1,SIZE(var%lev_sizes)) /),&
         !DSM (/ FILE%grid%nx, FILE%grid%ny, var%lev_sizes /))
  END IF

  RETURN
END IF

!-----------------------------------------------------------------------------
! If we get to here, we know we have to perform a subgrid extraction
!-----------------------------------------------------------------------------
IF ( .NOT. ASSOCIATED(subgrid%points) ) THEN

  ! If the subgrid is specified using a region, just extract the region from file

  IF ( FILE%grid%is_1d ) THEN
    CALL file_read_var(                                                       &
      FILE%fh, var%id, cube%values,                                           &
    ! We want to read the whole of the levels dimensions, but only the correct
    ! chunk of the grid dimensions, so set start and count accordingly
    ! For a 1D grid, we only need to consider the 'x axis' of the grid
            (/ subgrid%x_start, (1, i = 1,SIZE(var%lev_sizes)) /),            &
            (/ subgrid%nx,      var%lev_sizes /)                              &
          )
  ELSE
    CALL file_read_var(                                                       &
      FILE%fh, var%id, cube%values,                                           &
    ! We want to read the whole of the levels dimensions, but only the correct
    ! chunk of the grid dimensions, so set start and count accordingly
    ! For a 1D grid, we only need to consider the 'x axis' of the grid
            (/ subgrid%x_start, subgrid%y_start, (1, i = 1,SIZE(var%lev_sizes)) /), &
            (/ subgrid%nx,      subgrid%ny,      var%lev_sizes /)             &
          )
  END IF

ELSE

  ! If the subgrid is specified using points, then check whether we will do the
  ! extraction directly from file or in memory
  ! For now, we only extract directly from file if reading <= 1% of points
  ! Testing has shown this to be around the point where extracting in memory becomes
  ! faster for large files

  ! For this section, it is easier to deal with the data as if it has two grid
  ! dimensions and a combined z dimension that represents all the levels combined
  ALLOCATE(values_3d(subgrid%nx, subgrid%ny, PRODUCT(var%lev_sizes)))

  IF ( ( 100 * subgrid%nx * subgrid%ny ) > ( FILE%grid%nx * FILE%grid%ny ) ) THEN
    ! If the extraction is to take place in memory, then allocate and read the
    ! full grid of data
    ! Again, we use a combined z dimension here
    ! See above for explaination of start and count
    ALLOCATE(DATA(FILE%grid%nx, FILE%grid%ny, PRODUCT(var%lev_sizes)))
    IF ( FILE%grid%is_1d ) THEN
      CALL file_read_var(                                                     &
        FILE%fh, var%id, DATA,                                                &
        (/ 1,            (1, i = 1,SIZE(var%lev_sizes)) /),                   &
        (/ FILE%grid%nx, var%lev_sizes /)                                     &
      )
    ELSE
      CALL file_read_var(                                                     &
        FILE%fh, var%id, DATA,                                                &
        (/ 1,            1,            (1, i = 1,SIZE(var%lev_sizes)) /),     &
        (/ FILE%grid%nx, FILE%grid%ny, var%lev_sizes /)                       &
      )
    END IF

    ! Extract the points from the data - we use a utility routine to do this
    ! one level at a time
    DO k = 1,SIZE(DATA, 3)
      values_3d(:,:,k) = subgrid_extract(subgrid, DATA(:,:,k))
    END DO

    DEALLOCATE(DATA)

  ELSE

    ! Read the points one at a time directly from file
    DO j = 1,subgrid%ny
      DO i = 1,subgrid%nx
        ! Translate the index in subgrid%points into x and y coordinates in the file grid
        y = (subgrid%points(i,j) - 1) / FILE%grid%nx + 1
        x = subgrid%points(i,j) - (y-1) * FILE%grid%nx

        IF ( FILE%grid%is_1d ) THEN
          CALL file_read_var(FILE%fh, var%id, values_3d(i,j,:),               &
          ! For a 1D grid, we only need to specify the x index
                                         (/ x, (1, i = 1,SIZE(var%lev_sizes)) /), &
                                         (/ 1, var%lev_sizes /))
        ELSE
          CALL file_read_var(FILE%fh, var%id, values_3d(i,j,:),               &
          ! For a 2D grid, we specify the x and y indices
                                         (/ x, y, (1, i = 1,SIZE(var%lev_sizes)) /), &
                                         (/ 1, 1, var%lev_sizes /))
        END IF
      END DO
    END DO
  END IF

  ! Convert the 3d values back into a 1d array
  cube%values(:) = RESHAPE(values_3d, (/ SIZE(values_3d) /))
  DEALLOCATE(values_3d)

END IF  ! subgrid specified using region or points

RETURN

END FUNCTION file_gridded_read_var



!DSM ---- Lendo as variaveis do BRAMS ---{
SUBROUTINE read_drive_BRAMS(varID,var)

USE mem_brams_jules, ONLY: nxB,nyB,swdownB,lwdownB,diff_radB,precipB,tempB,upsB,vpsB,pstarB,qB

IMPLICIT NONE

INTEGER, INTENT(IN) :: varID 
INTEGER :: i, j, k

REAL :: var(nxB*nyB)

if (varID == 1) then
   k=1
   do j = 1 , nyB
      do i = 1 , nxB
         var(k) = swdownB(i,j)
         k=k+1
      enddo
   enddo
   !print*,'swdownB=',var(1:50)

else if (varID == 2) then
   k=1
   do j = 1 , nyB
      do i = 1 , nxB
         var(k) = lwdownB(i,j)
         k=k+1
      enddo
   enddo
   !print*,'lwdownB=',var(1:50)

else if (varID == 3) then

   k=1
   do j = 1 , nyB
      do i = 1 , nxB
         var(k) = diff_radB(i,j)
         k=k+1
      enddo
   enddo
   !print*,'diff_radB=',var(1:50)

else if (varID == 4) then
   k=1
   do j = 1 , nyB
      do i = 1 , nxB
         var(k) = precipB(i,j)
         k=k+1
      enddo
   enddo
   !print*,'precipB=',var(1:50)

else if (varID == 5) then
   k=1
   do j = 1 , nyB
      do i = 1 , nxB
         var(k) = tempB(i,j)
         k=k+1
      enddo
   enddo
   !print*,'tempB=',var(1:50)

else if (varID == 6) then
   k=1
   do j = 1 , nyB
      do i = 1 , nxB
         var(k) = upsB(i,j)
         k=k+1
      enddo
   enddo
   !print*,'upsB=',var(1:50)

else if (varID == 7) then
   k=1
   do j = 1 , nyB
      do i = 1 , nxB
         var(k) = vpsB(i,j)
         k=k+1
      enddo
   enddo
   !print*,'vpsB=',var(1:50)

else if (varID == 8) then
   k=1
   do j = 1 , nyB
      do i = 1 , nxB
         var(k) = pstarB(i,j)
         k=k+1
      enddo
   enddo
   !print*,'pstarB=',var(1:50)

else if (varID == 9) then
   k=1
   do j = 1 , nyB
      do i = 1 , nxB
         var(k) = qB(i,j)
         k=k+1
      enddo
   enddo
   !print*,'qB=',var(1:50)

end if


END SUBROUTINE read_drive_BRAMS
!DSM }


! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_gridded_write_var(FILE, var_id, cube, write_subgrid, subgrid)

USE grid_utils_mod, ONLY: subgrid_info, operator( /= )

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
TYPE(file_gridded), INTENT(INOUT) :: FILE
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
TYPE(var_gridded) :: var  ! The variable we are writing to

REAL, ALLOCATABLE :: values_3d(:,:,:)  ! 3D version of the cube values to
                                       ! simplify calculations using subgrids

INTEGER :: nx, ny  ! The sizes of the grid dimensions

INTEGER :: i, j, x, y  ! Index variables


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check if a subgrid has been requested and that the parent grid is the same
! as the grid the file is on
!-----------------------------------------------------------------------------
IF ( write_subgrid ) THEN
  IF ( .NOT. PRESENT(subgrid) )                                               &
    CALL log_fatal("file_gridded_write_var",                                  &
                   "Writing a subgrid has been requested but no " //          &
                   "subgrid has been given")

  ! Check that the subgrid's parent grid matches the file's grid
  IF ( subgrid%parent /= FILE%grid )                                          &
    CALL log_fatal("file_gridded_write_var",                                  &
                   "Subgrid is inconsistent with grid in file")
END IF  ! extract_subgrid

!-----------------------------------------------------------------------------
! Extract information on the variable being read into a local variable for
! convenience
!-----------------------------------------------------------------------------
var = FILE%vars(var_id)

!-----------------------------------------------------------------------------
! Check that the given cube has the correct size
!-----------------------------------------------------------------------------
! The cube should have two grid dimensions plus the correct number of levels
! dimensions
IF ( SIZE(cube%SHAPE) /= (2 + SIZE(var%lev_sizes)) )                          &
  CALL log_fatal("file_gridded_write_var",                                    &
                 "Given cube has incorrect number of dimensions")

! Check that the cube has the shape we are expecting
IF ( write_subgrid ) THEN
  nx = subgrid%nx
  ny = subgrid%ny
ELSE
  nx = FILE%grid%nx
  ny = FILE%grid%ny
END IF
IF ( .NOT. ALL(cube%SHAPE == (/ nx, ny, var%lev_sizes /)) )                   &
  CALL log_fatal("file_gridded_write_var",                                    &
                 "At least one dimension of given cube has incorrect size")

!-----------------------------------------------------------------------------
! Actually write the data
!
! The decision process for how to write to file is as follows:
!   * If no subgrid has been specified, write everything (obviously!)
!
!   * If a subgrid has been specified using a region, then write only that
!     region (easy to specify with start and count)
!
!   * If a subgrid has been specified using points, then write the values
!     to file one at a time. Note that the optimisation used in read_var
!     when the number of points in the subgrid is small compared to the
!     full grid does not apply here.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! If we are not writing a subgrid, just write all the data and return
!-----------------------------------------------------------------------------
IF ( .NOT. write_subgrid ) THEN
  IF ( FILE%grid%is_1d ) THEN
    CALL file_write_var(                                                      &
      FILE%fh, var%id, cube%values,                                           &
    ! We want to write the whole of the grid and levels dimensions, so set start
    ! and count accordingly
    ! For a 1D grid, we only need to consider the x axis of the grid
            (/ 1,            (1, i = 1,SIZE(var%lev_sizes)) /),               &
            (/ FILE%grid%nx, var%lev_sizes /)                                 &
          )
  ELSE
    CALL file_write_var(                                                      &
      FILE%fh, var%id, cube%values,                                           &
    ! We want to write the whole of the grid and levels dimensions, so set start
    ! and count accordingly
    ! For a 2D grid, we only need to consider both axes of the grid
            (/ 1,            1,            (1, i = 1,SIZE(var%lev_sizes)) /), &
            (/ FILE%grid%nx, FILE%grid%ny, var%lev_sizes /)                   &
          )
  END IF

  RETURN
END IF

!-----------------------------------------------------------------------------
! If we get to here, we know we have to write a subgrid
!-----------------------------------------------------------------------------
IF ( .NOT. ASSOCIATED(subgrid%points) ) THEN

  ! The subgrid is specified using a region - we just write a single slab
  ! using start/count appropriately

  IF ( FILE%grid%is_1d ) THEN
    CALL file_write_var(                                                      &
      FILE%fh, var%id, cube%values,                                           &
    ! We want to write the whole of the levels dimensions, but only the correct
    ! chunk of the grid dimensions, so set start and count accordingly
    ! For a 1D grid, we only need to consider the 'x axis' of the grid
            (/ subgrid%x_start, (1, i = 1,SIZE(var%lev_sizes)) /),            &
            (/ subgrid%nx,      var%lev_sizes /)                              &
          )
  ELSE
    CALL file_write_var(                                                      &
      FILE%fh, var%id, cube%values,                                           &
    ! We want to read the whole of the levels dimensions, but only the correct
    ! chunk of the grid dimensions, so set start and count accordingly
    ! For a 1D grid, we only need to consider the 'x axis' of the grid
            (/ subgrid%x_start, subgrid%y_start, (1, i = 1,SIZE(var%lev_sizes)) /), &
            (/ subgrid%nx,      subgrid%ny,      var%lev_sizes /)             &
          )
  END IF

ELSE

  ! The subgrid is specified using points - just write the points one at a time
  ! directly to file

  ! For this section, it is easier to deal with the data as if it has two grid
  ! dimensions and a combined z dimension that represents all the levels combined
  ALLOCATE(values_3d(subgrid%nx, subgrid%ny, PRODUCT(var%lev_sizes)))
  values_3d(:,:,:) = RESHAPE(cube%values, SHAPE(values_3d))

  DO j = 1,subgrid%ny
    DO i = 1,subgrid%nx
      ! Translate the index in subgrid%points into x and y coordinates in the file grid
      y = (subgrid%points(i,j) - 1) / FILE%grid%nx + 1
      x = subgrid%points(i,j) - (y-1) * FILE%grid%nx

      IF ( FILE%grid%is_1d ) THEN
        CALL file_write_var(FILE%fh, var%id, values_3d(i,j,:),                &
        ! For a 1D grid, we only need to specify the x index
                                      (/ x, (1, i = 1,SIZE(var%lev_sizes)) /),&
                                      (/ 1, var%lev_sizes /))
      ELSE
        CALL file_write_var(FILE%fh, var%id, values_3d(i,j,:),                &
        ! For a 2D grid, we specify the x and y indices
                                      (/ x, y, (1, i = 1,SIZE(var%lev_sizes)) /), &
                                      (/ 1, 1, var%lev_sizes /))
      END IF
    END DO
  END DO

  DEALLOCATE(values_3d)

END IF  ! subgrid specified using region or points


RETURN

END SUBROUTINE file_gridded_write_var
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_gridded_close(FILE)

USE file_mod, ONLY: file_close

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Closes and frees any resources consumed by the given gridded file object
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(file_gridded), INTENT(INOUT) :: FILE ! The file to close

! Work variables
INTEGER :: i  ! Loop index


!-----------------------------------------------------------------------------


! Close the underlying file handle
!DSM CALL file_close(FILE%fh)
if (index(trim(FILE%fh%ncdf%name)//' ',' ') > 2 .and. index(trim(FILE%fh%ncdf%name)//' ',' ')<450) CALL file_close(FILE%fh)

! Deallocate the sizes of any variables
DO i = 1,FILE%nvars
  IF ( ASSOCIATED(FILE%vars(i)%lev_sizes) ) THEN
    DEALLOCATE(FILE%vars(i)%lev_sizes)
    NULLIFY(FILE%vars(i)%lev_sizes)
  END IF
END DO

RETURN

END SUBROUTINE file_gridded_close

END MODULE file_gridded_mod
