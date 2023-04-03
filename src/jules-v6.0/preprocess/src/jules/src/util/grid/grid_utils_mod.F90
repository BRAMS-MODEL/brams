! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE grid_utils_mod

USE io_constants, ONLY: max_sdf_name_len

USE logging_mod, ONLY: log_fatal

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Type definitions
!-----------------------------------------------------------------------------
TYPE grid_info
  !-----------------------------------------------------------------------------
  ! This type encapsulates information about a grid
  !-----------------------------------------------------------------------------
  LOGICAL :: is_1d = .FALSE.

  CHARACTER(LEN=max_sdf_name_len) :: dim_name = ""
                               ! Used if is_1d = T
                               ! The name of the single grid dimension
  CHARACTER(LEN=max_sdf_name_len) :: x_name = "", y_name = ""
                               ! Used if is_1d = F
                               ! The names of the x and y dimensions of
                               ! the grid respectively

  INTEGER :: nx = 0, ny = 0  ! The sizes of the x and y dimensions of the grid
                             ! For 1d grids, ny = 1

END TYPE grid_info


TYPE subgrid_info

  ! Always set
  TYPE(grid_info) :: parent  ! The grid that this is a subgrid of

  INTEGER :: nx = 0, ny = 0  ! The size of the subgrid
                             ! For 1D subgrids, ny = 1

  ! Used if subgrid is constructed by specifying a region
  INTEGER :: x_start = -1, y_start = -1  ! The indices of the corner of the
                                         ! region
                                         ! The region consists of the area
                                         ! enclosed by the points
                                         !   (x_start,          y_start)
                                         !   (x_start + nx - 1, y_start)
                                         !   (x_start,          y_start + ny - 1)
                                         !   (x_start + nx - 1, y_start + ny - 1)

  ! Used if subgrid is constructed from point indices
  INTEGER, POINTER :: points(:,:) => NULL()
                                   ! For each point in the subgrid, this
                                   ! specifies the corresponding index in
                                   ! the parent grid
                                   ! For 2D grids, the first row is represented
                                   ! by indices 1:nx, the second by nx+1:2nx,
                                   ! the third by 2nx+1:3nx etc.

END TYPE subgrid_info

!-----------------------------------------------------------------------------
! Operator declarations
!-----------------------------------------------------------------------------
INTERFACE operator ( == )
MODULE PROCEDURE grid_eq
END INTERFACE

INTERFACE operator ( /= )
MODULE PROCEDURE grid_ne
END INTERFACE

!-----------------------------------------------------------------------------
! Overloads
!-----------------------------------------------------------------------------
! There are several ways to create a subgrid - we want to use them all by the
! same name externally
INTERFACE subgrid_create
MODULE PROCEDURE subgrid_create_mask, subgrid_create_region,                  &
                 subgrid_create_points
END INTERFACE subgrid_create

! A subgrid can also be further restricted in several ways - we also want to
! use them all be the same name
INTERFACE subgrid_restrict
MODULE PROCEDURE subgrid_restrict_mask, subgrid_restrict_region,              &
                 subgrid_restrict_points, subgrid_restrict_subgrid
END INTERFACE subgrid_restrict

!-----------------------------------------------------------------------------
! Visibility declarations
!-----------------------------------------------------------------------------
PRIVATE
PUBLIC                                                                        &
! Public types
    grid_info, subgrid_info,                                                  &
! Public operators
    operator ( == ), operator ( /= ),                                         &
! Public procedures
    grid_create, subgrid_create, subgrid_restrict, subgrid_extract


CONTAINS


! Includes for grid functionality
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION grid_create(grid_is_1d, dim_name, npoints, x_name, nx, y_name, ny)   &
                                                                  RESULT(grid)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Creates a grid_info object encapsulating the specified information
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
LOGICAL, INTENT(IN) :: grid_is_1d
    ! T - define a 1d grid (i.e. a vector)
    ! F - define a 2d grid

CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_name
    ! ONLY USED IF is_1d=T
    ! The name of the single grid dimension
INTEGER, INTENT(IN), OPTIONAL :: npoints
    ! ONLY USED IF is_1d=T
    ! The size of the single grid dimension

CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: x_name, y_name
    ! ONLY USED IF is_1d=F
    ! The names of the x and y dimensions of the grid respectively
INTEGER, INTENT(IN), OPTIONAL :: nx, ny
    ! ONLY USED IF is_1d=F
    ! The sizes of the x and y dimensions respectively

! Return type
TYPE(grid_info) :: grid  ! The grid_info object


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check an appropriate combination of variables has been given
!-----------------------------------------------------------------------------
IF ( grid_is_1d ) THEN
  IF ( ( .NOT. PRESENT(dim_name) ) .OR. ( .NOT. PRESENT(npoints) ) )          &
    CALL log_fatal("grid_create",                                             &
                   "To create a 1D grid, dim_name and npoints are required")
ELSE
  IF ( ( .NOT. PRESENT(x_name) ) .OR. ( .NOT. PRESENT(nx) ) .OR.              &
       ( .NOT. PRESENT(y_name) ) .OR. ( .NOT. PRESENT(ny) ) )                 &
    CALL log_fatal("grid_create",                                             &
                   "To create a 2D grid, x_name, y_name, nx and ny are " //   &
                   "all required")
END IF

!-----------------------------------------------------------------------------
! Build the grid_info object
!-----------------------------------------------------------------------------
grid%is_1d = grid_is_1d

IF ( grid_is_1d ) THEN
  grid%dim_name = dim_name

  grid%nx = npoints
  grid%ny = 1
ELSE
  grid%x_name = x_name
  grid%y_name = y_name

  grid%nx = nx
  grid%ny = ny
END IF

RETURN

END FUNCTION grid_create
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


LOGICAL FUNCTION grid_eq(grid, other_grid)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Determines if two grids are equal. Currently, two grids are considered
!   equal if they have the same dimensions (i.e. dimension names are not
!   taken into account).
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(grid_info), INTENT(IN) :: grid, other_grid  ! The grids to compare


!-----------------------------------------------------------------------------


grid_eq = ( grid%nx == other_grid%nx ) .AND. ( grid%ny == other_grid%ny )

RETURN

END FUNCTION grid_eq



LOGICAL FUNCTION grid_ne(grid, other_grid)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Determines if two grids are not equal.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(grid_info), INTENT(IN) :: grid, other_grid  ! The grids to compare


!-----------------------------------------------------------------------------


grid_ne = .NOT. grid_eq(grid, other_grid)

RETURN

END FUNCTION grid_ne

! Includes for subgrid functionality
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION subgrid_create_region(parent, x_start, y_start, nx, ny) RESULT(subgrid)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Creates a subgrid_info object representing the specified region
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(grid_info), INTENT(IN) :: parent  ! The parent grid

INTEGER, INTENT(IN) :: x_start, y_start, nx, ny
                                       ! The start and extent of the region
                                       ! the subgrid will cover

! Return type
TYPE(subgrid_info) :: subgrid  ! The subgrid object

! Work variables
INTEGER :: vertices(4,2)  ! The x/y coordinates of each of the 4 vertices
                          ! of the region, as determined by the starts and
                          ! extents


!-----------------------------------------------------------------------------


! Set the subgrid parent in all cases
subgrid%parent = parent

!-----------------------------------------------------------------------------
! If negative or 0 nx/ny have been specified, this results in an empty
! subgrid (i.e. return the subgrid in its default state)
!-----------------------------------------------------------------------------
IF ( nx < 1 .OR. ny < 1 ) RETURN

!-----------------------------------------------------------------------------
! Check that the subgrid lies within the parent grid
!-----------------------------------------------------------------------------
! Construct the vertices that define the region
vertices(1,:) = (/ x_start,          y_start /)
vertices(2,:) = (/ x_start + nx - 1, y_start /)
vertices(3,:) = (/ x_start,          y_start + ny - 1 /)
vertices(4,:) = (/ x_start + nx - 1, y_start + ny - 1 /)

! Check that all the values are within the bounds of the parent grid
IF ( ANY(vertices < 1) .OR.                                                   &
     ANY(vertices(:,1) > parent%nx) .OR. ANY(vertices(:,2) > parent%ny) )     &
  CALL log_fatal("subgrid_create_region",                                     &
                 "Specified region contains points beyond the bounds " //     &
                 "of the parent grid")

!-----------------------------------------------------------------------------
! Build and return the subgrid
!-----------------------------------------------------------------------------
subgrid%nx      = nx
subgrid%ny      = ny

subgrid%x_start = x_start
subgrid%y_start = y_start

RETURN

END FUNCTION subgrid_create_region


!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
!*****************************************************************************


FUNCTION subgrid_create_points(parent, points) RESULT(subgrid)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Creates a subgrid_info object from the specified points
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(grid_info), INTENT(IN) :: parent  ! The parent grid

INTEGER, INTENT(IN) :: points(:,:)  ! The indices in the parent grid of
                                    ! the points in the subgrid
                                    ! The shape of this argument will
                                    ! determine the shape of the subgrid

! Return type
TYPE(subgrid_info) :: subgrid  ! The subgrid object


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that the subgrid lies within the parent grid
!-----------------------------------------------------------------------------
! Check that all the values are within the bounds of the parent grid
IF ( ANY(points < 1) .OR. ANY(points > (parent%nx * parent%ny)) )             &
  CALL log_fatal("subgrid_create_points",                                     &
                 "Points beyond the bounds of the parent grid have " //       &
                 "been specified")

!-----------------------------------------------------------------------------
! Build and return the subgrid
!-----------------------------------------------------------------------------
subgrid%parent  = parent

subgrid%nx      = SIZE(points, 1)
subgrid%ny      = SIZE(points, 2)

ALLOCATE(subgrid%points(subgrid%nx,subgrid%ny))
subgrid%points(:,:) = points(:,:)

RETURN

END FUNCTION subgrid_create_points


!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
!*****************************************************************************


FUNCTION subgrid_create_mask(parent, mask, force_1d_grid) RESULT(subgrid)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Creates a subgrid_info object from the specified mask
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(grid_info), INTENT(IN) :: parent  ! The parent grid

LOGICAL, INTENT(IN) :: mask(:,:)  ! Mask indicating which points in the
                                  ! parent grid to use
                                  ! This should be the same size as the
                                  ! parent grid
                                    
LOGICAL, INTENT(IN) :: force_1d_grid  ! T - model grid is 1D
                          ! F - model grid is 1D unless input grid is 2D 
                          !     and model grid is the whole input grid

! Return type
TYPE(subgrid_info) :: subgrid  ! The subgrid object


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! When using a mask to create a subgrid, we would like to use a region
! if possible
! So that the logic to detect this is contained in one place, we create a
! "null subgrid" - a subgrid specified as a region that covers the whole
! grid - and just pass that to subgrid_restrict with the same mask
!-----------------------------------------------------------------------------
! By specifying x/y_start = 1 and nx/ny the same as the parent grid, we create
! a "null subgrid"
subgrid = subgrid_restrict(                                                   &
  subgrid_create(parent, 1, 1, parent%nx, parent%ny), mask, force_1d_grid     &
)

RETURN

END FUNCTION subgrid_create_mask
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION subgrid_restrict_region(parent, x_start, y_start, nx, ny) RESULT(subgrid)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes a subgrid_info object and returns a subgrid_info object that is
!   equivalent to extracting the specified region from the original subgrid
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(subgrid_info), INTENT(IN) :: parent  ! The parent subgrid

INTEGER, INTENT(IN) :: x_start, y_start, nx, ny
                                       ! The start and extent of the region
                                       ! the subgrid will cover

! Return type
TYPE(subgrid_info) :: subgrid  ! The restricted subgrid object

! Work variables
INTEGER :: vertices(4,2)  ! The x/y coordinates of each of the 4 vertices
                          ! of the region, as determined by the starts and
                          ! extents


!-----------------------------------------------------------------------------


! Set the subgrid parent to be the actual grid that is the parent of the
! parent subgrid
subgrid%parent = parent%parent

!-----------------------------------------------------------------------------
! If negative or 0 nx/ny have been specified, this results in an empty
! subgrid (i.e. return the subgrid in its default state)
!-----------------------------------------------------------------------------
IF ( nx < 1 .OR. ny < 1 ) RETURN

!-----------------------------------------------------------------------------
! Check that the subgrid lies within the parent grid
!-----------------------------------------------------------------------------
! Construct the vertices that define the region
vertices(1,:) = (/ x_start,          y_start /)
vertices(2,:) = (/ x_start + nx - 1, y_start /)
vertices(3,:) = (/ x_start,          y_start + ny - 1 /)
vertices(4,:) = (/ x_start + nx - 1, y_start + ny - 1 /)

! Check that all the values are within the bounds of the parent grid
IF ( ANY(vertices < 1) .OR.                                                   &
     ANY(vertices(:,1) > parent%nx) .OR. ANY(vertices(:,2) > parent%ny) )     &
  CALL log_fatal("subgrid_restrict_region",                                   &
                 "Specified region contains points beyond the bounds " //     &
                 "of the parent grid")

!-----------------------------------------------------------------------------
! Build and return the subgrid
!-----------------------------------------------------------------------------
subgrid%nx = nx
subgrid%ny = ny

IF ( ASSOCIATED(parent%points) ) THEN
  ! If the parent subgrid is specified using points, then our new subgrid can
  ! only be specified using points
  ! We extract the points from the parent subgrid associated with the specified
  ! region
  ALLOCATE(subgrid%points(subgrid%nx,subgrid%ny))
  subgrid%points(:,:) = parent%points(vertices(1,1):vertices(4,1),            &
                                      vertices(1,2):vertices(4,2))
ELSE
  ! If the parent subgrid is specified using a region, then we can just restrict
  ! that region by offsetting the given x/y_start by the parent x/y_start
  subgrid%x_start = parent%x_start + x_start - 1
  subgrid%y_start = parent%y_start + y_start - 1
END IF

RETURN

END FUNCTION subgrid_restrict_region


!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
!*****************************************************************************


FUNCTION subgrid_restrict_points(parent, points) RESULT(subgrid)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes a subgrid_info object and returns a subgrid_info object that is
!   equivalent to extracting the specified points from the original subgrid
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(subgrid_info), INTENT(IN) :: parent  ! The parent subgrid

INTEGER, INTENT(IN) :: points(:,:)  ! The indices in the parent subgrid of
                                    ! the points in the restricted subgrid
                                    ! The shape of this argument will
                                    ! determine the shape of the restricted
                                    ! subgrid

! Return type
TYPE(subgrid_info) :: subgrid  ! The restricted subgrid object

! Work variables
INTEGER :: i, j, x, y  ! Index variables


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that the subgrid lies within the parent subgrid
!-----------------------------------------------------------------------------
! Check that all the values are within the bounds of the parent grid
IF ( ANY(points < 1) .OR. ANY(points > (parent%nx * parent%ny)) )             &
  CALL log_fatal("subgrid_restrict_points",                                   &
                 "Points beyond the bounds of the parent grid have " //       &
                 "been specified")

!-----------------------------------------------------------------------------
! Build and return the subgrid
!-----------------------------------------------------------------------------
! The parent of the new subgrid is the same as the old subgrid
subgrid%parent  = parent%parent

! The size of the new subgrid is determined by the size of the points array
subgrid%nx = SIZE(points, 1)
subgrid%ny = SIZE(points, 2)

ALLOCATE(subgrid%points(subgrid%nx,subgrid%ny))

IF ( ASSOCIATED(parent%points) ) THEN
  ! If the parent subgrid is specified using points, then we construct our new
  ! subgrid by extracting those points from the parent subgrid
  DO j = 1,subgrid%ny
    DO i = 1,subgrid%nx
      ! Convert the 1D index in points into x/y indices in the parent subgrid
      y = (points(i,j) - 1) / parent%nx + 1
      x = points(i,j) - (y-1) * parent%nx

      subgrid%points(i,j) = parent%points(x,y)
    END DO
  END DO
ELSE
  ! If the parent subgrid is specified using a region, we construct our subgrid
  ! by offsetting each specified point using the parent subgrid's x_start and
  ! y_start
  DO j = 1,subgrid%ny
    DO i = 1,subgrid%nx
      ! Convert the 1D index in points into x/y indices in the parent subgrid
      y = (points(i,j) - 1) / parent%nx + 1
      x = points(i,j) - (y-1) * parent%nx

      ! Offset those indices by the parent subgrid's x/y_start to get the x/y indices
      ! in the full parent grid
      x = x + parent%x_start - 1
      y = y + parent%y_start - 1

      ! Convert back to a 1D index in the full parent grid
      subgrid%points(i,j) = (y-1) * subgrid%parent%nx + x
    END DO
  END DO
END IF

RETURN

END FUNCTION subgrid_restrict_points


!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
!*****************************************************************************


FUNCTION subgrid_restrict_mask(parent, mask, force_1d_grid) RESULT(subgrid)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes a subgrid_info object and returns a subgrid_info object that is
!   equivalent to extracting the points from the original subgrid where
!   mask is .TRUE.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(subgrid_info), INTENT(IN) :: parent  ! The parent subgrid

LOGICAL, INTENT(IN) :: mask(:,:)  ! Mask indicating which points in the
                                  ! parent grid to use
                                  ! This should be the same size as the
                                  ! parent grid
                                    
LOGICAL, INTENT(IN) :: force_1d_grid  ! T - model grid is 1D
                          ! F - model grid is 1D unless input grid is 2D 
                          !     and model grid is the whole input grid

! Return type
TYPE(subgrid_info) :: subgrid  ! The subgrid object

! Work variables
INTEGER :: nx, ny  ! The dimensions of the subgrid

INTEGER :: x_start, y_start  ! The index of the first column/row to contain
                             ! a .TRUE. value

INTEGER :: x_end, y_end  ! The index of the last column/row to contain
                         ! a .TRUE. value

INTEGER, ALLOCATABLE :: points(:,:)  ! The indices of the .TRUE. points
                                     ! in mask

INTEGER :: i, j, n  ! Index variables


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! First, verify that the mask is applicable to the parent subgrid
!-----------------------------------------------------------------------------
IF ( SIZE(mask,1) /= parent%nx .OR. SIZE(mask,2) /= parent%ny )               &
  CALL log_fatal("subgrid_restrict_mask",                                     &
                 "Mask must have the same dimensions as the grid it is " //   &
                 "being applied to")

!-----------------------------------------------------------------------------
! If mask contains no .TRUE. points, we return a 0 size subgrid by specifying
! a region with nx = ny = 0
!-----------------------------------------------------------------------------
IF ( .NOT. ANY(mask) ) THEN
  subgrid = subgrid_restrict(parent, 1, 1, 0, 0)
  RETURN
END IF

!-----------------------------------------------------------------------------
! Try to detect if the mask specifies a region - we do this because using a
! subgrid is more efficient (particularly when reading/writing subgrids
! during I/O)
!
! The fallback is to use a list of points
!-----------------------------------------------------------------------------
! Find the first column that contains a .TRUE. in mask
DO n = 1,parent%nx
  IF ( ANY(mask(n,:)) ) THEN
    x_start = n
    EXIT
  END IF
END DO

! Find the last column that contains a .TRUE.
DO n = parent%nx,1,-1
  IF ( ANY(mask(n,:)) ) THEN
    x_end = n
    EXIT
  END IF
END DO

! Find the first row that contains a .TRUE. in mask
DO n = 1,parent%ny
  IF ( ANY(mask(:,n)) ) THEN
    y_start = n
    EXIT
  END IF
END DO

! Find the last row that contains a .TRUE.
DO n = parent%ny,1,-1
  IF ( ANY(mask(:,n)) ) THEN
    y_end = n
    EXIT
  END IF
END DO

!-----------------------------------------------------------------------------
! We have found the extents of the .TRUE. points in mask (i.e. mask has no
! .TRUE. points outside of the bounds we have found
!
! If all the points within the bounds we have found are also .TRUE., we have
! a region
! Otherwise, we fall back to a list of points
!-----------------------------------------------------------------------------
IF ( ALL(mask(x_start:x_end,y_start:y_end)) .AND.                             &
    (parent%ny == 1 .OR. ( .NOT. force_1d_grid )) ) THEN
  ! Calculate the nx/ny we will pass to the region overload of subgrid_restrict
  ! to create the subgrid
  ! This is inclusive at both ends
  nx = x_end - x_start + 1
  ny = y_end - y_start + 1

  subgrid = subgrid_restrict(parent, x_start, y_start, nx, ny)
ELSE
  ! Build a list of points to pass to the points overload of subgrid_restrict
  nx = COUNT(mask)
  ny = 1

  ALLOCATE(points(nx,ny))

  n = 0
  DO j = 1,parent%ny
    DO i = 1,parent%nx
      IF ( mask(i,j) ) THEN
        n = n + 1
        points(n,1) = (j-1) * parent%nx + i
      END IF
    END DO
  END DO

  subgrid = subgrid_restrict(parent, points)

  DEALLOCATE(points)
END IF

RETURN

END FUNCTION subgrid_restrict_mask


!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
!*****************************************************************************


FUNCTION subgrid_restrict_subgrid(parent, subgrid) RESULT(new_subgrid)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes two subgrid objects and combines them into a new subgrid. The size
!   of the parent grid of the second subgrid must match the size of the first
!   subgrid
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(subgrid_info), INTENT(IN) :: parent  ! The parent subgrid

TYPE(subgrid_info), INTENT(IN) :: subgrid  ! The subgrid to use to restrict
                                           ! the parent subgrid

! Return type
TYPE(subgrid_info) :: new_subgrid  ! The combined subgrid object

! Work variables


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! First, verify that the combination of subgrids can be done
!-----------------------------------------------------------------------------
IF ( subgrid%parent%nx /= parent%nx .OR. subgrid%parent%ny /= parent%ny )     &
  CALL log_fatal("subgrid_restrict_subgrid",                                  &
                 "Grid dimensions are incompatible")

!-----------------------------------------------------------------------------
! From this point, we can defer to one of the other subgrid_restrict_*
! variants, depending on whether the second subgrid is restricting by
! region or points
!-----------------------------------------------------------------------------
IF ( ASSOCIATED(subgrid%points) ) THEN
  new_subgrid = subgrid_restrict(parent, subgrid%points)
ELSE
  new_subgrid = subgrid_restrict(                                             &
    parent, subgrid%x_start, subgrid%y_start, subgrid%nx, subgrid%ny          &
  )
END IF

RETURN

END FUNCTION subgrid_restrict_subgrid
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION subgrid_extract(subgrid, grid_data) RESULT(subgrid_data)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes data on the parent grid and a subgrid_info object, extracts the
!   data on the subgrid and returns it
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(subgrid_info), INTENT(IN) :: subgrid
                   ! The subgrid to extract
REAL, INTENT(IN) :: grid_data(subgrid%parent%nx,subgrid%parent%ny)
                   ! The data on the parent grid to extract from

! Return type
REAL :: subgrid_data(subgrid%nx,subgrid%ny)

! Work variables
INTEGER :: i, j, x, y  ! Index variables

!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! How we extract data depends on whether the subgrid is defined as a region
! or a set of points
!-----------------------------------------------------------------------------
IF ( ASSOCIATED(subgrid%points) ) THEN
  DO j = 1,subgrid%ny
    DO i = 1,subgrid%nx
      ! Convert the point indices in subgrid%points into x/y coords in the parent grid
      y = (subgrid%points(i,j) - 1) / subgrid%parent%nx + 1
      x = subgrid%points(i,j) - (y-1) * subgrid%parent%nx

      ! Set the subgrid point to the value of the corresponding point in the parent grid
      subgrid_data(i,j) = grid_data(x,y)
    END DO
  END DO
ELSE
  ! If the subgrid is specified using a region, just extract that region from the
  ! parent grid
  subgrid_data(:,:) = grid_data(                                              &
    subgrid%x_start:(subgrid%x_start + subgrid%nx-1),                         &
    subgrid%y_start:(subgrid%y_start + subgrid%ny-1)                          &
  )
END IF

RETURN

END FUNCTION subgrid_extract

END MODULE grid_utils_mod
