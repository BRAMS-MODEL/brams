#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE parallel_mod

USE logging_mod, ONLY: log_info, log_fatal

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   This module provides utilities used for running in parallel using MPI
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Module parameters
!-----------------------------------------------------------------------------
INTEGER, PARAMETER :: master_task_id = 0

!-----------------------------------------------------------------------------
! Module variables
!-----------------------------------------------------------------------------
! These are set by calling decompose_domain
INTEGER :: ntasks = 0   ! The number of MPI tasks available
INTEGER :: task_id = -1  ! The ID of this task

INTEGER :: mpi_type_global_col  ! MPI type for columns in the full grid
INTEGER :: mpi_type_local_col   ! MPI type for columns in the task local grid

! Used in scatter/gather calls
INTEGER, ALLOCATABLE :: counts(:)  ! The number of blocks sent to each
                                   ! task. This is 1 for every task
INTEGER, ALLOCATABLE :: offsets(:)  ! The offset, in numbers of blocks,
                                    ! for each task's block in the global
                                    ! array

!-----------------------------------------------------------------------------
! Visibilities
!-----------------------------------------------------------------------------
PRIVATE
PUBLIC master_task_id, is_master_task, decompose_domain,                      &
                       scatter_land_field, gather_land_field


CONTAINS


#include "is_master_task.inc"
#include "decompose_domain.inc"
#include "scatter_land_field.inc"
#include "gather_land_field.inc"

END MODULE parallel_mod
#endif
