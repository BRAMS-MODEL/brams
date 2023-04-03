#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE gather_land_field(field_local_land, field_global_land)

USE mpi, ONLY: mpi_comm_world

USE missing_data_mod, ONLY: rmdi

USE ancil_info, ONLY: land_pts

USE jules_fields_mod, ONLY: ainfo
  
USE theta_field_sizes, ONLY: t_i_length, t_j_length
  
USE model_grid_mod, ONLY: model_grid, global_land_pts, global_land_mask

IMPLICIT NONE
  
! Interface definition is only required if using the dummy MPI library - it
! allows the mpi_gatherv implementation to use assumed-shape arrays (without
! being in a module), which it needs to as it doesn't track MPI type information
! and relies purely on the shape of the array
#if defined(MPI_DUMMY)
INTERFACE
SUBROUTINE mpi_gatherv(sendbuf, sendcnt, sendtype,                            &
                       recvbuf, recvcnts, displs, recvtype,                   &
                       root, comm, error)

REAL, INTENT(IN) :: sendbuf(:,:)
INTEGER, INTENT(IN) :: sendcnt
INTEGER, INTENT(IN) :: sendtype

REAL, INTENT(OUT) :: recvbuf(:,:)
INTEGER, INTENT(IN) :: recvcnts(:)
INTEGER, INTENT(IN) :: displs(:)
INTEGER, INTENT(IN) :: recvtype

INTEGER, INTENT(IN) :: root
INTEGER, INTENT(IN) :: comm
INTEGER, INTENT(OUT) :: error

END SUBROUTINE mpi_gatherv
END INTERFACE
#endif

!-----------------------------------------------------------------------------
! Description:
!   Takes a field defined on the land points of each task and gathers it
!   onto the land points of the full model grid in the master task, in
!   accordance with the decomposition performed by decompose_domain
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
REAL, INTENT(IN) :: field_local_land(land_pts)  
                                       ! The field on the land points for the
                                       ! current task
                                       ! This must be defined for all tasks
REAL, INTENT(OUT) :: field_global_land(:)
                                     ! The field on the global land points
                                     ! This only has to be properly defined
                                     ! for the master task
                                         
! Work variables
REAL, ALLOCATABLE :: field_global_2d(:,:)
                                ! The global field on the full 2d model grid
                                ! Only allocated in master task
REAL :: field_local_2d(t_i_length,t_j_length)
                                ! The local field on the model grid for the
                                ! current task
                                ! Used in all tasks
                                  
INTEGER :: i,j,l  ! Indexing variables

INTEGER :: error  ! Error indicator for MPI calls
                  ! This is ignored as (most) MPI implementations fail
                  ! rather than returning actual error codes


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Verify that the global land points field is the right size
! We do this at runtime to avoid forcing it to be defined at that size for
! non-master tasks
! We also allocate the global 2d array at this point, for the same reason
!-----------------------------------------------------------------------------
IF ( is_master_task() ) THEN
  IF ( SIZE(field_global_land) /= global_land_pts )                           &
    CALL log_fatal("gather_land_field",                                       &
                   "Output field should be on global land points")
                                          
  ALLOCATE(field_global_2d(model_grid%nx,model_grid%ny))
ELSE
  ALLOCATE(field_global_2d(1,1))
END IF

!-----------------------------------------------------------------------------
! Convert the local land point array into a 2d array on the model grid for
! the task
!-----------------------------------------------------------------------------
field_local_2d(:,:) = rmdi
DO l = 1,land_pts
  j = (ainfo%land_index(l) - 1) / t_i_length + 1
  i = ainfo%land_index(l) - (j-1) * t_i_length

  field_local_2d(i,j) = field_local_land(l)
END DO

!-----------------------------------------------------------------------------
! Gather the values from each task into the global field in the main task
! using the block datatype and calculated counts and offsets
!-----------------------------------------------------------------------------
CALL mpi_gatherv(field_local_2d,  t_i_length,      mpi_type_local_col,        &
                 field_global_2d, counts, offsets, mpi_type_global_col,       &
                 master_task_id, mpi_comm_world, error)

!-----------------------------------------------------------------------------
! In the master task, convert the 2d array on the full model grid into
! an array on the global land points
!-----------------------------------------------------------------------------
IF ( is_master_task() )                                                       &
  field_global_land(:) = PACK(field_global_2d, global_land_mask)
  
!-----------------------------------------------------------------------------
! Deallocate the global field at the end
!-----------------------------------------------------------------------------
IF ( ALLOCATED(field_global_2d) ) DEALLOCATE(field_global_2d)

RETURN

END SUBROUTINE gather_land_field
#endif
