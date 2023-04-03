!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine reduce_max_cfl_to_master(master_num, mynum, ngrids, cflxy, cflz)
  implicit none
  integer, intent(in) :: master_num
  integer, intent(in) :: mynum
  integer, intent(in) :: ngrids
  real,    intent(inout) :: cflxy(ngrids)
  real,    intent(inout) :: cflz(ngrids)

  include "mpif.h"

  integer :: ierr
  real    :: send_buff(2*ngrids)
  real    :: recv_buff(2*ngrids)

  ! build send buffer

  if (mynum == master_num) then
     send_buff(:) = 0.0
  else
     send_buff(1:ngrids) = cflxy(1:ngrids)
     send_buff(ngrids+1:2*ngrids) = cflz(1:ngrids)
  end if

  ! build receiving buffer

  recv_buff(:)  = 0.0

  ! reduce send buffer by MAX operation, 
  ! storing value at master's receiving buffer

  call MPI_REDUCE(send_buff, recv_buff, 2*ngrids, MPI_REAL, &
       MPI_MAX, master_num, MPI_COMM_WORLD, ierr)

  ! master unpacks receiving buffer

  if (mynum == master_num) then
     cflxy(1:ngrids) = recv_buff(1:ngrids)
     cflz (1:ngrids) = recv_buff(ngrids+1:2*ngrids)
  end if
end subroutine reduce_max_cfl_to_master

!-------------------------------------------------------------------------

subroutine reduce_max_cfl_and_broadcast(master_num, mynum, ngrids, cflxy, cflz)
  implicit none
  integer, intent(in) :: master_num
  integer, intent(in) :: mynum
  integer, intent(in) :: ngrids
  real,    intent(inout) :: cflxy(ngrids)
  real,    intent(inout) :: cflz(ngrids)

  include "mpif.h"

  integer :: ierr
  real    :: send_buff(2*ngrids)
  real    :: recv_buff(2*ngrids)

  ! build send buffer

  if (mynum == master_num) then
     send_buff(:) = 0.0
  else
     send_buff(1:ngrids) = cflxy(1:ngrids)
     send_buff(ngrids+1:2*ngrids) = cflz(1:ngrids)
  end if

  ! reduce send buffer by MAX operation, 
  ! storing result at receiving buffer

  call MPI_ALLREDUCE(send_buff, recv_buff, 2*ngrids, MPI_REAL, &
       MPI_MAX, MPI_COMM_WORLD, ierr)

  ! unpack receiving buffer

  cflxy(1:ngrids) = recv_buff(1:ngrids)
  cflz (1:ngrids) = recv_buff(ngrids+1:2*ngrids)
end subroutine reduce_max_cfl_and_broadcast

!-------------------------------------------------------------------------

subroutine gather_cpu_time_master_print(master_num, mynum, nmachs, tcpu)
  implicit none
  integer, intent(in) :: master_num
  integer, intent(in) :: mynum
  integer, intent(in) :: nmachs
  real,    intent(in) :: tcpu

  include "mpif.h"

  integer :: ierr
  integer :: first, last
  real    :: recv_buff(0:nmachs)

  call MPI_GATHER (&
       tcpu,      1, MPI_REAL, &
       recv_buff, 1, MPI_REAL, &
       master_num, MPI_COMM_WORLD, ierr)
  
  if (mynum == master_num) then
     write(*,"(a)") "; CPU per slave (first slave : last slave):"
     do first = 1, nmachs, 10
        last=min(nmachs, first+9)
        write(*,"(' (',i4.4,':',i4.4,')',10f8.2)") first, last, recv_buff(first:last)
     end do
  end if
end subroutine gather_cpu_time_master_print

 subroutine AllReduceMeteogramNpts(npts, totalPts)

  integer, intent(in) :: npts
  integer :: totalPts

  include "mpif.h"
    
   totalPts = 0  
   call MPI_Allreduce(npts, totalPts, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
 
 end subroutine AllReduceMeteogramNpts

 subroutine GatherMeteogram(masterNum, myNum, nCities, nVars, localSums, &
            localMins, localMaxs, fUnit, totalSums, totalMins, totalMaxs)
   implicit none
   integer, intent(in)                          :: masterNum
   integer, intent(in)                          :: myNum
   integer, intent(in)                          :: nCities
   integer, intent(in)                          :: nVars
   real, dimension(nCities, nVars), intent(in)  :: localSums
   real, dimension(nCities, nVars), intent(in)  :: localMaxs
   real, dimension(nCities, nVars), intent(in)  :: localMins
   integer, intent(in)                          :: fUnit
   real, dimension(nCities, nVars), intent(inout) :: totalSums
   real, dimension(nCities, nVars), intent(inout) :: totalMins
   real, dimension(nCities, nVars), intent(inout) :: totalMaxs
   
    include "mpif.h"
  
    integer :: nv
    integer :: nc
    integer :: ierr

       call MPI_Reduce(localSums, totalSums, nCities*nVars, MPI_REAL, MPI_SUM, masterNum, MPI_COMM_WORLD, ierr)
       call MPI_Reduce(localMaxs, totalMaxs, nCities*nVars, MPI_REAL, MPI_MAX, masterNum, MPI_COMM_WORLD, ierr)
       call MPI_Reduce(localMins, totalMins, nCities*nVars, MPI_REAL, MPI_MIN, masterNum, MPI_COMM_WORLD, ierr)   
       
 end subroutine GatherMeteogram
 
  subroutine MeteogramBarrier()
    
    implicit none
   
    include "mpif.h"
  
    integer :: ierr
  
    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
       
 end subroutine MeteogramBarrier
