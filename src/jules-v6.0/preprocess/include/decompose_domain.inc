#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION decompose_domain(grid) RESULT(task_subgrid)

USE mpi, ONLY: mpi_comm_world, mpi_address_kind, mpi_real

USE grid_utils_mod, ONLY: grid_info, subgrid_info, subgrid_create

USE string_utils_mod, ONLY: to_string

IMPLICIT NONE
  
!-----------------------------------------------------------------------------
! Description:
!   Decomposes the given grid across the available MPI tasks
!   Returns a subgrid object representing the part of the grid that the
!   current task will be responsible for
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(grid_info), INTENT(IN) :: grid  ! The grid to decompose

! Return type
TYPE(subgrid_info) :: task_subgrid  ! The subgrid this task is responsible for

! Work variables
INTEGER :: ntasks_x, ntasks_y  ! The size of the "task grid" (i.e. the grid
                               ! will be split into ntasks_x by ntasks_y
                               ! blocks

INTEGER :: task_nx, task_ny    ! The size of the block for the current task

INTEGER :: task_x, task_y  ! The x and y coordinates in the "task grid"
                           ! of the current task

INTEGER :: x_start, y_start  ! The x/y coordinates in the grid of the start
                             ! of the task subgrid

LOGICAL :: found_decomposition  ! T - we found a usable decomposition
                                ! F - we did not find a usable decomposition

INTEGER :: leftover  ! The remainder when distributing columns along a task row

INTEGER :: mpi_type  ! Holds intermediate MPI datatype before the extent
                     ! is adjusted
INTEGER(KIND=mpi_address_kind) :: mpi_real_lb, mpi_real_extent
                     ! The lower bound and extent for the mpi_real type
                       
INTEGER(KIND=mpi_address_kind), PARAMETER :: mpi_zero = 0
                     ! A 'zero' of the correct kind to be used as an MPI address

INTEGER, ALLOCATABLE :: counts_2d(:,:), offsets_2d(:,:)
                     ! Used when calculating the counts and offsets on the
                     ! task grid

INTEGER :: error  ! Error indicator for MPI calls
                  ! This is ignored as (most) MPI implementations fail
                  ! rather than returning actual error codes

INTEGER :: i, j, n  ! Index variables


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! This routine currently implements a very naive decomposition
!
! The main concern when decomposing the grid is I/O, not MPI communication
! I.E. the grid must be split into contiguous regions that can be written to
! file in one write statement by specifying appropriate start and count
!-----------------------------------------------------------------------------

! First get the number of available tasks and the id of this task
!DSM CALL mpi_comm_size(mpi_comm_world, ntasks, error)
  ntasks  = 1  !DSM
  error = 0 !DSM
!DSM CALL mpi_comm_rank(mpi_comm_world, task_id, error)
  task_id  = 0 !DSM
  error = 0 !DSM

CALL log_info("decompose_domain",                                             &
              "Decomposing domain across " // TRIM(to_string(ntasks)) //      &
              " available tasks")

! We can only utilise at most 1 task per point
IF ( ntasks > grid%nx * grid%ny )                                             &
  CALL log_fatal("decompose_domain",                                          &
                 "More tasks are available than points in the model grid")

!-----------------------------------------------------------------------------
! Work out the decomposition, subject to the following rules:
!
!   * Each task must have the same number of grid rows, but can have a
!     varying number of columns
!
!   * Each row of the task grid must have at most grid%nx tasks
!
!   * Each row of the task grid must have the same number of tasks
!
!-----------------------------------------------------------------------------
! This is the minimum number of rows we need in the task grid to ensure that
! each row has <= grid%nx tasks
ntasks_y = (ntasks-1) / grid%nx + 1

! Loop until we find a suitable number of rows for the task grid
! Limited testing found that using as many rows of tasks as possible resulted
! in the most efficient decomposition more of the time (I realise that sounds
! a bit wooly!)
DO n = grid%ny,ntasks_y,-1
  found_decomposition =     ( MOD(grid%ny, n) == 0 ) & ! Each task gets the same number of grid rows
                      .AND. ( MOD(ntasks, n) == 0 )     ! Each row of the task grid has the same number of tasks

  IF ( found_decomposition ) THEN
    ntasks_x = ntasks / n
    ntasks_y = n
    EXIT
  END IF
END DO

! If we could not find a suitable decomposition, suggest changing the number
! of available processes
IF ( .NOT. found_decomposition )                                              &
  CALL log_fatal("decompose_domain",                                          &
                 "Unable to find a suitable decomposition - try " //          &
                 "using a different number of tasks")

CALL log_info("decompose_domain",                                             &
              "Tasks are arranged as a grid of size " //                      &
              TRIM(to_string(ntasks_x)) // " x " // TRIM(to_string(ntasks_y)))

! Each task has the same number of rows, which we can now calculate
task_ny = grid%ny / ntasks_y

!-----------------------------------------------------------------------------
! Build the MPI datatype that allows us to scatter to and gather from
! global arrays in blocks of size 1 x task_ny
!-----------------------------------------------------------------------------
! Get the lower bound and extent for the mpi_real type
!DSM CALL mpi_type_get_extent(mpi_real, mpi_real_lb, mpi_real_extent, error)
  mpi_real_lb = 0 !DSM
  mpi_real_extent = 4 !DSM
  error  = 0 !DSM

! Define a MPI type that selects columns from the full grid
!DSM CALL mpi_type_vector(task_ny, 1, grid%nx, mpi_real, mpi_type, error)
  mpi_type = 1 !DSM

! Restrict the extent of the datatype to 1 real value for use in offset
! calculations
!DSM CALL mpi_type_create_resized(                                                 &
!DSM  mpi_type, mpi_zero, mpi_real_extent, mpi_type_global_col, error             &
!DSM )
  mpi_type_global_col = 1 !DSM

! Commit the datatype
!DSM CALL mpi_type_commit(mpi_type_global_col, error)

!-----------------------------------------------------------------------------
! Work out how many grid columns each task in the task grid will have
!-----------------------------------------------------------------------------
ALLOCATE(counts_2d(ntasks_x,ntasks_y))
! Work out how many columns (most of) the tasks will get
counts_2d(:,:) = grid%nx / ntasks_x
! Distribute any leftover columns
leftover = MOD(grid%nx, ntasks_x)
IF ( leftover > 0 ) THEN
  DO n = 1,leftover
    counts_2d(n,:) = counts_2d(n,:) + 1
  END DO
END IF

! The counts used in MPI calls are a flattened version of this
ALLOCATE(counts(ntasks))
counts(:) = RESHAPE(counts_2d, (/ ntasks /))

!-----------------------------------------------------------------------------
! Calculate the offsets for each task
!
! Because we adjusted the extent of the column type, our offsets are
! calculated in actual grid cells
!-----------------------------------------------------------------------------
ALLOCATE(offsets_2d(ntasks_x,ntasks_y))
! Note that MPI offsets must start from 0, not 1!
DO j = 1,ntasks_y
  DO i = 1,ntasks_x
    ! First, sum along each row in the task grid to get the offsets within the row
    offsets_2d(i,j) = SUM(counts_2d(1:i-1,j))
  END DO
  ! Then for each row of tasks, add the offset to the start of that row
  offsets_2d(:,j) = offsets_2d(:,j) + (j-1) * grid%nx * task_ny
END DO

! The offsets used in MPI calls are a flattened version of this
ALLOCATE(offsets(ntasks))
offsets(:) = RESHAPE(offsets_2d, (/ ntasks /))

!-----------------------------------------------------------------------------
! Construct the subgrid that this task is responsible for
!-----------------------------------------------------------------------------
! Work out where the current task sits in the task grid
! Remember that task ids start from 0, not 1!
task_y = task_id / ntasks_x + 1
task_x = (task_id+1) - (task_y-1) * ntasks_x

! From that, we can retrieve the number of columns the current task is
! reponsible for
! We already know how many rows the task is responsible for
task_nx = counts_2d(task_x,task_y)

! Now we work out the position of the start of our grid in the full grid
! Don't forget that these offsets start from 1, not 0!
! We can get our x_start by summing the counts for the previous columns
x_start = SUM(counts_2d(1:task_x-1,task_y)) + 1
! Since each task has the same number of rows, getting y_start is easier
y_start = (task_y-1) * task_ny + 1
  
task_subgrid = subgrid_create(grid, x_start, y_start, task_nx, task_ny)

!-----------------------------------------------------------------------------
! Now we know the size of the subgrid, we define a MPI datatype that selects
! columns in the subgrid for the current task
!-----------------------------------------------------------------------------
!DSM CALL mpi_type_vector(                                                         &
!DSM  task_subgrid%ny, 1, task_subgrid%nx, mpi_real, mpi_type, error              &
!DSM )
  mpi_type = 1 !DSM

! Again, we set the extent of the type to one real value for offset
! calculations
!DSMCALL mpi_type_create_resized(                                                 &
!DSM  mpi_type, mpi_zero, mpi_real_extent, mpi_type_local_col, error              &
!DSM )
  mpi_type_global_col = 1 !DSM

!DSM CALL mpi_type_commit(mpi_type_local_col, error)

DEALLOCATE(counts_2d)
DEALLOCATE(offsets_2d)

RETURN

END FUNCTION decompose_domain
#endif
