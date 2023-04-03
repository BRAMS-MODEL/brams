#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE scatter_land_field(field_global_land, field_local_land)

USE mpi, ONLY: mpi_comm_world

USE missing_data_mod, ONLY: rmdi

USE ancil_info, ONLY: land_pts

USE jules_fields_mod, ONLY: ainfo
  
USE theta_field_sizes, ONLY: t_i_length, t_j_length
  
USE model_grid_mod, ONLY: model_grid, global_land_pts, global_land_mask

IMPLICIT NONE
  
! Interface definition is only required if using the dummy MPI library - it
! allows the mpi_scatterv implementation to use assumed-shape arrays (without
! being in a module), which it needs to as it doesn't track MPI type information
! and relies purely on the shape of the array
#if defined(MPI_DUMMY)
INTERFACE
SUBROUTINE mpi_scatterv(sendbuf, sendcnts, displs, sendtype,                  &
                        recvbuf, recvcnt, recvtype,                           &
                        root, comm, error)

REAL, INTENT(IN) :: sendbuf(:,:)
INTEGER, INTENT(IN) :: sendcnts(:)
INTEGER, INTENT(IN) :: displs(:)
INTEGER, INTENT(IN) :: sendtype

REAL, INTENT(OUT) :: recvbuf(:,:)
INTEGER, INTENT(IN) :: recvcnt
INTEGER, INTENT(IN) :: recvtype

INTEGER, INTENT(IN) :: root
INTEGER, INTENT(IN) :: comm
INTEGER, INTENT(OUT) :: error

END SUBROUTINE mpi_scatterv
END INTERFACE
#endif

!-----------------------------------------------------------------------------
! Description:
!   Takes a field defined on the land points of the full model grid in the
!   master task and scatters it onto the land points for each task, in
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
REAL, INTENT(IN) :: field_global_land(:)
                                     ! The field on the global land points
                                     ! This only has to be properly defined
                                     ! for the master task
REAL, INTENT(OUT) :: field_local_land(land_pts)  
                                       ! The field on the land points for the
                                       ! current task
                                       ! This must be defined for all tasks
                                         
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
! In the master task, convert the full land points array into a 2d array
! on the full model grid
!-----------------------------------------------------------------------------
IF ( is_master_task() ) THEN
  ! Check that the global field is properly defined in the master task
  IF ( SIZE(field_global_land) /= global_land_pts )                           &
    CALL log_fatal("scatter_land_field",                                      &
                   "Input field should be on global land points")
                     
  ! We then need to map back onto the full model grid before scattering
  ALLOCATE(field_global_2d(model_grid%nx,model_grid%ny))
  field_global_2d(:,:) = rmdi
  field_global_2d(:,:) = UNPACK(                                              &
    field_global_land, global_land_mask, field_global_2d                      &
  )
ELSE
  ALLOCATE(field_global_2d(1,1))
END IF
  
!-----------------------------------------------------------------------------
! Scatter the field into local 2d parts using the block datatype and
! calculated counts and offsets
!-----------------------------------------------------------------------------
CALL mpi_scatterv(field_global_2d, counts, offsets, mpi_type_global_col,      &
                  field_local_2d,  t_i_length,      mpi_type_local_col,       &
                  master_task_id, mpi_comm_world, error)

!-----------------------------------------------------------------------------
! Convert the local 2d parts into land point arrays
!-----------------------------------------------------------------------------
DO l = 1,land_pts
  j = (ainfo%land_index(l) - 1) / t_i_length + 1
  i = ainfo%land_index(l) - (j-1) * t_i_length

  field_local_land(l) = field_local_2d(i,j)
END DO
  
!-----------------------------------------------------------------------------
! Deallocate the global field at the end
!-----------------------------------------------------------------------------
IF ( ALLOCATED(field_global_2d) ) DEALLOCATE(field_global_2d)

RETURN

END SUBROUTINE scatter_land_field
#endif
