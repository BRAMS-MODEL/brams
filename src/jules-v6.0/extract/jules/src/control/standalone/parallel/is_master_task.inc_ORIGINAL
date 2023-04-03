#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


LOGICAL FUNCTION is_master_task()

USE mpi, ONLY: mpi_comm_world

IMPLICIT NONE
  
!-----------------------------------------------------------------------------
! Description:
!   Returns .TRUE. if the current task is the master task, .FALSE. otherwise
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
INTEGER :: task_id  ! The id of this task
  
INTEGER :: error  ! Error indicator for MPI calls
                  ! This is ignored as (most) MPI implementations fail
                  ! rather than returning actual error codes


!-----------------------------------------------------------------------------

CALL mpi_comm_rank(mpi_comm_world, task_id, error)
  
is_master_task = (task_id == master_task_id)

RETURN

END FUNCTION is_master_task
#endif
