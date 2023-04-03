!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


!---------------------------------------------------------------------------

!For Carma - CATT
subroutine mk_buff_carma(a,b,n1,n2,n3,m1,m2,m3,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,n3,m1,m2,m3,i1,i2,j1,j2
  real :: a(n1,n2,n3),b(m1,m2,m3)

  b(1:m1,1:m2,1:m3)=a(i1:i2,j1:j2,1:n3)

  return
end subroutine mk_buff_carma

!---------------------------------------------------------------------------

subroutine mk_2_buff(a,b,n1,n2,m1,m2,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,m1,m2,i1,i2,j1,j2
  real :: a(n1,n2),b(m1,m2)

  b(1:m1,1:m2)=a(i1:i2,j1:j2)

  return
end subroutine mk_2_buff

!---------------------------------------------------------------------------

subroutine mk_2p_buff(a,b,n1,n2,n3,m1,m2,m3,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,n3,m1,m2,m3,i1,i2,j1,j2
  real :: a(n1,n2,n3),b(m1,m2,m3)

  b(1:m1,1:m2,1:m3)=a(i1:i2,j1:j2,1:n3)

  return
end subroutine mk_2p_buff

!---------------------------------------------------------------------------

subroutine mk_3_buff(a,b,n1,n2,n3,m1,m2,m3,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,n3,m1,m2,m3,i1,i2,j1,j2
  real :: a(n1,n2,n3),b(m1,m2,m3)

  b(1:m1,1:m2,1:m3)=a(1:n1,i1:i2,j1:j2)

  return
end subroutine mk_3_buff

!---------------------------------------------------------------------------

subroutine mk_4_buff(a,b,n1,n2,n3,n4,m1,m2,m3,m4,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,n3,n4,m1,m2,m3,m4,i1,i2,j1,j2
  real :: a(n1,n2,n3,m4),b(m1,m2,m3,m4)

  b(1:m1,1:m2,1:m3,1:m4)=a(1:n1,i1:i2,j1:j2,1:n4)

  return
end subroutine mk_4_buff

!---------------------------------------------------------------------------

subroutine ex_buff_carma(a,b,n1,n2,n3,m1,m2,m3,i0,j0,i1,i2,j1,j2)
  !i0 = ixoff(nm,ng)
  !j0 = iyoff(nm,ng)
  !i1 = il1
  !i2 = ir2
  !j1 = jb1
  !j2 = jt2
  implicit none
  integer :: n1,n2,n3,m1,m2,m3,i0,j0,i1,i2,j1,j2
  real :: a(n1,n2,n3),b(m1,m2,m3)

  a(i1+i0:i2+i0,j1+j0:j2+j0,1:n3) = b(i1:i2,j1:j2,1:m3)

  return
end subroutine ex_buff_carma

!---------------------------------------------------------------------------

subroutine ex_2_buff(a,b,n1,n2,m1,m2,i0,j0,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,m1,m2,i0,j0,i1,i2,j1,j2
  real :: a(n1,n2),b(m1,m2)

  a(i1+i0:i2+i0,j1+j0:j2+j0) = b(i1:i2,j1:j2)

  return
end subroutine ex_2_buff

!---------------------------------------------------------------------------

subroutine ex_3_buff(a,b,n1,n2,n3,m1,m2,m3,i0,j0,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,n3,m1,m2,m3,i0,j0,i1,i2,j1,j2
  real,intent(inout) :: a(n1,n2,n3),b(m1,m2,m3)

  a(1:n1,i1+i0:i2+i0,j1+j0:j2+j0) = b(1:m1,i1:i2,j1:j2)

  return
end subroutine ex_3_buff

!---------------------------------------------------------------------------

subroutine ex_4_buff(a,b,n1,n2,n3,n4,m1,m2,m3,m4,i0,j0,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,n3,n4,m1,m2,m3,m4,i0,j0,i1,i2,j1,j2
  real :: a(n1,n2,n3,n4),b(m1,m2,m3,m4)

  a(1:n1,i1+i0:i2+i0,j1+j0:j2+j0,1:n4) = b(1:m1,i1:i2,j1:j2,1:m4)

  return
end subroutine ex_4_buff

!---------------------------------------------------------------------------

subroutine ex_2p_buff(a,b,n1,n2,n3,m1,m2,m3,i0,j0,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,n3,m1,m2,m3,i0,j0,i1,i2,j1,j2
  real :: a(n1,n2,n3),b(m1,m2,m3)

  a(i1+i0:i2+i0,j1+j0:j2+j0,1:n3) = b(i1:i2,j1:j2,1:m3)

  return
end subroutine ex_2p_buff

!---------------------------------------------------------------------------

subroutine GatherAllChunks (nmachs, master_num, toGather, sizeToGather, &
     gathered, sizeGathered, eachSizeGathered, displacement)
  implicit none
  integer, intent(in) :: nmachs
  integer, intent(in) :: master_num
  integer, intent(in) :: sizeToGather
  integer, intent(in) :: sizeGathered
  real, intent(in) :: toGather(sizeToGather)
  real, intent(out) :: gathered(sizeGathered)
  integer, intent(in) :: eachSizeGathered(nmachs)
  integer, intent(in) :: displacement(nmachs)

  include "mpif.h"
  integer :: ierr
  character(len=8) :: c0, c1, c2, c3
  character(len=*), parameter :: h="**(GatherAllChunks)**" 
  logical, parameter :: dumpLocal=.false.

  if (dumpLocal) then
     write(c0,"(i8)") sizeToGather
     write(c1,"(i8)") sizeGathered
     write(c2,"(i8)") eachSizeGathered(1)
     write(c3,"(i8)") displacement(1)
     write(*,"(a)") h//" sizeToGather="//trim(adjustl(c0))//&
          "; sizeGathered="//trim(adjustl(c1))//&
          "; eachSizeGathered(1)="//trim(adjustl(c2))//&
          "; displacement(1)="//trim(adjustl(c3))
     call flush(6)
     write(c0,"(e8.1)") sum(toGather)
     write(*,"(a)") h//" sum(toGather)="//trim(adjustl(c0))
     call flush(6)
     write(c0,"(e8.1)") sum(gathered)
     write(*,"(a)") h//" sum(gathered)="//trim(adjustl(c0))
     call flush(6)
  end if
  ! gather a field

  call MPI_GATHERV(toGather, sizeToGather, MPI_REAL, &
       gathered, eachSizeGathered, displacement, MPI_REAL, &
       master_num, MPI_COMM_WORLD, ierr)

  if (ierr /= MPI_SUCCESS) THEN
     write(c0,"(i8)") ierr
     call fatal_error(h//" gatherv fails with ierr="//trim(adjustl(c0)))
  else if (dumpLocal) then
     write(*,"(a)") h//" done"
     call flush(6)
  end if
end subroutine GatherAllChunks

!---------------------------------------------------------------------------

subroutine StoreGathered (idim_type, nmachs, sizeGathered, gathered, &
     nnxp, nnyp, nnzp, nzs, nzg, npatch, nwave, &
     nodemxp, nodemyp, nodemzp, nodeia, nodeiz, nodeja, nodejz, &
     nodei0, nodej0, il1, ir2, jb1, jt2, displacement, eachSizeGathered, &
     sizeStored, stored)
  implicit none
  integer, intent(in) :: idim_type
  integer, intent(in) :: nmachs
  integer, intent(in) :: sizeGathered
  real,    intent(in) :: gathered(sizeGathered)
  integer, intent(in) :: nnxp
  integer, intent(in) :: nnyp
  integer, intent(in) :: nnzp
  integer, intent(in) :: nzs
  integer, intent(in) :: nzg
  integer, intent(in) :: npatch
  integer, intent(in) :: nwave
  integer, intent(in) :: nodemxp(nmachs)
  integer, intent(in) :: nodemyp(nmachs)
  integer, intent(in) :: nodemzp(nmachs)
  integer, intent(in) :: nodeia(nmachs)
  integer, intent(in) :: nodeiz(nmachs)
  integer, intent(in) :: nodeja(nmachs)
  integer, intent(in) :: nodejz(nmachs)
  integer, intent(in) :: nodei0(nmachs)
  integer, intent(in) :: nodej0(nmachs)
  integer, intent(in) :: il1(nmachs)
  integer, intent(in) :: ir2(nmachs)
  integer, intent(in) :: jb1(nmachs)
  integer, intent(in) :: jt2(nmachs)
  integer, intent(in) :: displacement(nmachs)
  integer, intent(in) :: eachSizeGathered(nmachs)
  integer, intent(in) :: sizeStored
  real,    intent(out) :: stored(sizeStored)

  integer :: proc
  character(len=8) :: c0
  character(len=*), parameter :: h="**(StoreGathered)**"

  ! unpack gathered field, removing unnecessary ghost zones and
  ! placing entries at correct field positions

  select case (idim_type)
  case (2)
     do proc = 1, nmachs
        if (eachSizeGathered(proc)/=0) then
           call ex_2_buff(stored, gathered(displacement(proc)+1), &
                nnxp, nnyp, nodemxp(proc), nodemyp(proc), &
                nodei0(proc), nodej0(proc), &
                il1(proc), ir2(proc), jb1(proc), jt2(proc))
        end if
     end do
  case (3)
     do proc = 1, nmachs
        if (eachSizeGathered(proc)/=0) then
           call ex_3_buff(stored, gathered(displacement(proc)+1), &
                nnzp, nnxp, nnyp, nnzp, nodemxp(proc), nodemyp(proc), &
                nodei0(proc), nodej0(proc), &
                il1(proc), ir2(proc), jb1(proc), jt2(proc))
        end if
     end do
  case (4)
     do proc = 1, nmachs
        if (eachSizeGathered(proc)/=0) then
           call ex_4_buff(stored, gathered(displacement(proc)+1), &
                nzg, nnxp, nnyp, npatch, nzg, nodemxp(proc), nodemyp(proc), npatch, &
                nodei0(proc), nodej0(proc), &
                il1(proc), ir2(proc), jb1(proc), jt2(proc))
        end if
     end do
  case (5)
     do proc = 1, nmachs
        if (eachSizeGathered(proc)/=0) then
           call ex_4_buff(stored, gathered(displacement(proc)+1), &
                nzs, nnxp, nnyp, npatch, nzs, nodemxp(proc), nodemyp(proc), npatch, &
                nodei0(proc), nodej0(proc), &
                il1(proc), ir2(proc), jb1(proc), jt2(proc))
        end if
     end do
  case (6)
     do proc = 1, nmachs
        if (eachSizeGathered(proc)/=0) then
           call ex_2p_buff(stored, gathered(displacement(proc)+1), &
                nnxp, nnyp, npatch, nodemxp(proc), nodemyp(proc), npatch, &
                nodei0(proc), nodej0(proc), &
                il1(proc), ir2(proc), jb1(proc), jt2(proc))
        end if
     end do
  case (7)
     do proc = 1, nmachs
        if (eachSizeGathered(proc)/=0) then
           call ex_buff_carma(stored, gathered(displacement(proc)+1), &
                nnxp, nnyp, nwave, nodemxp(proc), nodemyp(proc), nwave, &
                nodei0(proc), nodej0(proc), &
                il1(proc), ir2(proc), jb1(proc), jt2(proc))
        end if
     end do
  case default
     write(c0,"(i8)") idim_type
     call fatal_error(h//" invoked with unknown idim_type ("//&
          trim(adjustl(c0))//")")
  end select
end subroutine StoreGathered

!---------------------------------------------------------------------------

subroutine SizesDisplacements (idim_type, nmachs, &
     nnxp, nnyp, nnzp, nzs, nzg, npatch, nwave, &
     nodemxp, nodemyp, nodemzp, nodeia, nodeiz, nodeja, nodejz, nodeibcon, &
     il1, ir2, jb1, jt2, eachSizeGathered, displacement)
  implicit none
  integer, intent(in) :: idim_type 
  integer, intent(in) :: nmachs
  integer, intent(in) :: nnxp
  integer, intent(in) :: nnyp
  integer, intent(in) :: nnzp
  integer, intent(in) :: nzs
  integer, intent(in) :: nzg
  integer, intent(in) :: npatch
  integer, intent(in) :: nwave
  integer, intent(in) :: nodemxp(nmachs)
  integer, intent(in) :: nodemyp(nmachs)
  integer, intent(in) :: nodemzp(nmachs)
  integer, intent(in) :: nodeia(nmachs)
  integer, intent(in) :: nodeiz(nmachs)
  integer, intent(in) :: nodeja(nmachs)
  integer, intent(in) :: nodejz(nmachs)
  integer, intent(in) :: nodeibcon(nmachs)
  integer, intent(out) :: il1(nmachs)
  integer, intent(out) :: ir2(nmachs)
  integer, intent(out) :: jb1(nmachs)
  integer, intent(out) :: jt2(nmachs)
  integer, intent(out) :: eachSizeGathered(nmachs)
  integer, intent(out) :: displacement(nmachs)

  integer :: proc
  character(len=8) :: c0
  character(len=*), parameter :: h="**(SizesDisplacements)**" 

  ! field portion to be gathered at each process

  select case (idim_type)
  case (2)
     do proc = 1, nmachs
        eachSizeGathered(proc) = nodemxp(proc)*nodemyp(proc)
     end do
  case (3)
     do proc = 1, nmachs
        eachSizeGathered(proc) = nodemzp(proc)*nodemxp(proc)*nodemyp(proc)
     end do
  case (4)
     do proc = 1, nmachs
        eachSizeGathered(proc) = nzg*nodemxp(proc)*nodemyp(proc)*npatch
     end do
  case (5)
     do proc = 1, nmachs
        eachSizeGathered(proc) = nzs*nodemxp(proc)*nodemyp(proc)*npatch
     end do
  case (6)
     do proc = 1, nmachs
        eachSizeGathered(proc) = nodemxp(proc)*nodemyp(proc)*npatch
     end do
  case (7)
     do proc = 1, nmachs
        eachSizeGathered(proc) = nodemxp(proc)*nodemyp(proc)*nwave
     end do
  case default
     write(c0,"(i8)") idim_type
     call fatal_error(h//" invoked with unknown idim_type ("//&
          trim(adjustl(c0))//")")
  end select

  ! where to place field portion at gathering result

  displacement(1)=0
  do proc = 2, nmachs
     displacement(proc)=displacement(proc-1)+eachSizeGathered(proc-1)
  end do

  ! eliminate unnecessary ghost zones while unpacking gathered result

  do proc = 1, nmachs
     if (btest(nodeibcon(proc),0)) then
        il1(proc) = nodeia(proc) - 1
     else
        il1(proc) = nodeia(proc)
     end if
     if (btest(nodeibcon(proc),1)) then
        ir2(proc) = nodeiz(proc) + 1
     else
        ir2(proc) = nodeiz(proc)
     end if
     if (btest(nodeibcon(proc),2)) then
        jb1(proc) = nodeja(proc) - 1
     else
        jb1(proc) = nodeja(proc)
     end if
     if (btest(nodeibcon(proc),3)) then
        jt2(proc) = nodejz(proc) + 1
     else
        jt2(proc) = nodejz(proc)
     end if
  end do
end subroutine SizesDisplacements

!---------------------------------------------------------------------------

subroutine OneGatherStoreAllChunks(mchnum, mynum, nmachs, master_num, &
     toGather, idim_type, eachSizeGathered, displacement, &
     nnxp, nnyp, nnzp, nzs, nzg, npatch, nwave, &
     nodemxp, nodemyp, nodemzp, nodeia, nodeiz, nodeja, nodejz, &
     nodei0, nodej0, il1, ir2, jb1, jt2, &
     sizeStored, stored)
  implicit none
  integer, intent(in) :: mchnum
  integer, intent(in) :: mynum
  integer, intent(in) :: nmachs
  integer, intent(in) :: master_num
  real,    intent(in) :: toGather(*)
  integer, intent(in) :: idim_type
  integer, intent(in) :: eachSizeGathered(nmachs)
  integer, intent(in) :: displacement(nmachs)
  integer, intent(in) :: nnxp
  integer, intent(in) :: nnyp
  integer, intent(in) :: nnzp
  integer, intent(in) :: nzs
  integer, intent(in) :: nzg
  integer, intent(in) :: npatch
  integer, intent(in) :: nwave
  integer, intent(in) :: nodemxp(nmachs)
  integer, intent(in) :: nodemyp(nmachs)
  integer, intent(in) :: nodemzp(nmachs)
  integer, intent(in) :: nodeia(nmachs)
  integer, intent(in) :: nodeiz(nmachs)
  integer, intent(in) :: nodeja(nmachs)
  integer, intent(in) :: nodejz(nmachs)
  integer, intent(in) :: il1(nmachs)
  integer, intent(in) :: ir2(nmachs)
  integer, intent(in) :: jb1(nmachs)
  integer, intent(in) :: jt2(nmachs)
  integer, intent(in) :: nodei0(nmachs)
  integer, intent(in) :: nodej0(nmachs)
  integer, intent(in) :: sizeStored
  real,    intent(out) :: stored(sizeStored)

  integer :: ierr
  integer :: sizeGathered
  real, allocatable :: gathered(:)
  character(len=8) :: c0, c1
  character(len=*), parameter :: h="**(OneGatherStoreAllChunks)**"
  logical, parameter :: dumpLocal=.false.

  ! size of gathered field (contains all ghost zones of all processes)

  sizeGathered = displacement(nmachs)+eachSizeGathered(nmachs)
  allocate(gathered(sizeGathered), stat=ierr)
  if (ierr /= 0) then
     write(c0,"(i8)") sizeGathered
     write(c1,"(i8)") ierr
     call fatal_error(h//" allocate gathered("//trim(adjustl(c0))//&
          ") failed with stat="//trim(adjustl(c1)))
  else if (dumpLocal) then
     write(c0,"(i8)") sizeGathered
     write(*,"(a)") h//" allocated gathered("//trim(adjustl(c0))//")"
     call flush(6)
  end if

  ! gathered field at master_num

  if (dumpLocal) then
     write(c0,"(i8)") eachSizeGathered(mynum)
     write(*,"(a)") h//" will gather with local size "//trim(adjustl(c0))
     call flush(6)
  end if
  call GatherAllChunks (nmachs, master_num, toGather, eachSizeGathered(mynum), &
     gathered(1), sizeGathered, eachSizeGathered, displacement)
  if (dumpLocal) then
     write(*,"(a)") h//" done gathering"
  end if

  ! master_num unpacks fields (removes unnecessary ghost zones and positions entries)
  
  if (mchnum == master_num) then
     if (dumpLocal) then
        write(*,"(a)") h//" master will StoreGathered"
     end if
     call StoreGathered (idim_type, nmachs, sizeGathered, gathered, &
          nnxp, nnyp, nnzp, nzs, nzg, npatch, nwave, &
          nodemxp, nodemyp, nodemzp, nodeia, nodeiz, nodeja, nodejz, &
          nodei0, nodej0, il1, ir2, jb1, jt2, displacement, eachSizeGathered, &
          sizeStored, stored)
     if (dumpLocal) then
        write(*,"(a)") h//" done StoreGathered"
     end if
  end if

  deallocate(gathered, stat=ierr)
  if (ierr /= 0) then
     write(c1,"(i8)") ierr
     call fatal_error(h//" deallocate gathered failed with stat="//trim(adjustl(c1)))
  end if
end subroutine OneGatherStoreAllChunks

!---------------------------------------------------------------------------

subroutine AllGatherStoreAllChunks(varType)

  use node_mod, only: &
       nmachs, master_num, mchnum, mynum,&
       nodemxp, nodemyp, &
       nodeia, nodeiz, nodeja, nodejz, &
       nodei0, nodej0, nodeibcon
      
  use mem_grid, only: &
       ngrids, nnxp, nnyp, nnzp, nzs, nzg, npatch, &
       time

  use mem_aerad, only: &
       nwave

  use var_tables, only: &
       num_var, &
       vtab_r


  use an_header, only: &
       head_table

  use mem_basic, only:  &
       basic_g

  use mem_turb,  only:   &
       turb_g, idiffk,  &
       xkhkm

  use io_params, only: &
       ioutput, &
       iclobber


  implicit none
  include "i8.h"
  character(len=*), intent(in) :: varType

  integer :: idim_type
  integer :: ng
  integer :: nv
  integer :: il1(nmachs,2:7)
  integer :: ir2(nmachs,2:7)
  integer :: jb1(nmachs,2:7)
  integer :: jt2(nmachs,2:7)
  integer :: eachSizeGathered(nmachs,2:7)
  integer :: displacement(nmachs,2:7)

  include "files.h"

  character(len=8)              :: c0, c1, c2
  character(len=3)              :: cProc
  character(len=*),   parameter :: h="**(AllGatherStoreAllChunks)**" 

  integer                       :: ierr
  integer                       :: sizeRes
  real,             allocatable :: Res(:)
  real, pointer                 :: field
  integer                       :: flag
  character(len=16)             :: varn
  character(len=f_name_length)  :: anamel(ngrids)
  character(len=f_name_length)  :: anamelh
  integer                       :: nvtot
  integer                       :: nvcnt
  integer                       :: nvMax
  integer                       :: ioaunt
  real                          :: timeold
  type(head_table), allocatable :: aw_table(:)
  integer(kind=i8)              :: npointer
  real,             pointer     :: scr1(:)
  integer                       :: nodemxyzp
  logical,          allocatable :: WillGather(:,:)

  ! return if no output is selected

  if (ioutput == 0) return
  ioaunt=10
  write(cProc,"(i3.3)") mchnum

  ! field independent constants for gather and unpacking

  do ng = 1, ngrids
     do idim_type=2,7
        call SizesDisplacements (idim_type, nmachs, &
             nnxp(ng), nnyp(ng), nnzp(ng), nzs, nzg, npatch, nwave, &
             nodemxp(:,ng), nodemyp(:,ng), nnzp(ng), &
             nodeia(:,ng), nodeiz(:,ng), nodeja(:,ng), nodejz(:,ng), &
             nodeibcon(:,ng),il1(:,idim_type), ir2(:,idim_type), &
             jb1(:,idim_type), jt2(:,idim_type),eachSizeGathered(:,idim_type),&
             displacement(:,idim_type))
     end do
  end do

  ! gather fields just for output
  ! find which fields should be gathered for output (value of varType)
  ! and store findings at WillGather

  nvMax = maxval(num_var(1:ngrids))
  allocate(WillGather(nvMax,ngrids), stat=ierr)
  if (ierr /= 0) then
     write(c0,"(i8)") nvMax
     write(c1,"(i8)") ngrids
     write(c2,"(i8)") ierr
     call fatal_error(h//" allocate WillGather("//trim(adjustl(c0))//&
          ","//trim(adjustl(c1))//") failed with stat="//trim(adjustl(c2)))
  end if

  do ng = 1, ngrids
     do nv = 1, num_var(ng)
        select case (varType)
        case ("INST")
           WillGather(nv,ng)= vtab_r(nv,ng)%ianal==1
        case ("LITE")
           WillGather(nv,ng)= vtab_r(nv,ng)%ilite==1
        case ("MEAN")
           WillGather(nv,ng)= vtab_r(nv,ng)%imean==1
        case ("BOTH")
           WillGather(nv,ng)= vtab_r(nv,ng)%ilite==1
        case default
           call fatal_error(h//" Proc "//cProc//" unknown varType=**"// &
                trim(varType)//"**")
        end select
     end do
     do nv = num_var(ng)+1, nvMax
        WillGather(nv,ng)=.false.
     end do
  end do

  ! how many fields to gather and output

  nvtot = count(WillGather)

  ! master builds output filename and create aw_table data structure

  if (mchnum == master_num) then

!**(JP)** nao vejo necessidade de descobrir se vai imprimir agora ou nao, pois
!**(JP)** AllGatherStoreAllChunks soh deve ser invocada quando for tempo de impressao,
!**(JP)** definido pelos valores de isendflg, isendlite, isendmean e isendboth em par_model,
!**(JP)** model, rams_node e OneProc. Se este calculo for diferente do que eh feito nessas
!**(JP)** variaveis, o calculo deve ser transferido para lah.

!!$     writeAnlFile = frqWrt(iflag,varType,time,avgtim,timmax,frqanl,frqlite, &
!!$          frqmean,frqboth,dtlongn)
!!$     if (writeAnlFile) then

     call AnlwrtInitialize(varType,anamel,anamelh,timeold)
     allocate(aw_table(nvtot), stat=ierr)
     if (ierr /= 0) then
        write(c0,"(i8)") nvtot
        write(c1,"(i8)") ierr
        call fatal_error(h//" allocate aw_table("//trim(adjustl(c0))//&
             ") failed with stat="//trim(adjustl(c1)))
     end if
  end if

  ! for all grids

  do ng = 1, ngrids
     
     ! master opens output file and initialize output counters
     
     if (mchnum == master_num) then
        call OpenAnlwrt(anamel)
        npointer = 0_i8 
        nvcnt = 0
     end if

     ! for all fields selected for output

     do nv = 1, num_var(ng)

        ! field to gather

        if (WillGather(nv,ng)) then
           select case (varType)
           case ("INST")
              field => vtab_r(nv,ng)%var_p
           case ("LITE")
              field => vtab_r(nv,ng)%var_p
           case ("MEAN")
              field => vtab_r(nv,ng)%var_m
           case ("BOTH")
              field => vtab_r(nv,ng)%var_m
           case default
              call fatal_error(h//" Proc "//cProc//" unknown varType=**"// &
                   trim(varType)//"**")
           end select
           
           ! allocate scratch area to store pre-processed field
           
           nodemxyzp=nodemxp(mynum,ng)*nodemyp(mynum,ng)* &
                nnzp(ng)
           allocate(scr1(nodemxyzp), stat=ierr)
           if (ierr /= 0) then
              write(c0,"(i8)") nodemxyzp
              write(c1,"(i8)") ierr
              call fatal_error(h//" allocate scr1("// &
                   trim(adjustl(c0))//") failed with stat="// &
                   trim(adjustl(c1)))
           end if
           varn= vtab_r(nv,ng)%name
           
           ! pre-process field before gathering
           
           if (varn == 'PP') then
              
              ! Output total Exner function
              
              call RAMS_aprep_p (nodemxyzp,field &
                   ,basic_g(ng)%pi0,scr1 )
              varn='PI'
              field =>scr1(1)
           elseif(varn == 'HKM') then
              
              ! Convert to HKM to HKH (note that VKH is HKH for Deardorff)
              
              call RAMS_aprep_hkh (nodemxyzp,field &
                   ,turb_g(ng)%vkh,basic_g(ng)%dn0 &
                   ,scr1,idiffk(ng),xkhkm(ng))
              varn='HKH'
              field =>scr1(1)
           elseif(varn == 'VKH') then
              
              ! Un-density weight VKH
              
              call RAMS_aprep_vkh (nodemxyzp,field &
                   ,basic_g(ng)%dn0,scr1)
              field =>scr1(1)
           end if

           ! classify field type

           idim_type=vtab_r(nv,ng)%idim_type
           if (idim_type < 2 .or. idim_type > 7) then
              write(c0,"(i8)") idim_type
              call fatal_error(h//" Proc "//cProc//"  unknown idim_type="// &
                   trim(adjustl(c0)))
           end if

           ! unpacked field size

           if (mchnum == master_num) then
              select case (idim_type)
              case (2)
                 sizeRes=nnxp(ng)*nnyp(ng)
              case (3)
                 sizeRes=nnzp(ng)*nnxp(ng)*nnyp(ng)
              case (4)
                 sizeRes=nzg*nnxp(ng)*nnyp(ng)*npatch
              case (5)
                 sizeRes=nzs*nnxp(ng)*nnyp(ng)*npatch
              case (6)
                 sizeRes=nnxp(ng)*nnyp(ng)*npatch
              case (7)
                 sizeRes=nnxp(ng)*nnyp(ng)*nwave
              case default
                 write(c0,"(i8)") vtab_r(nv,ng)%idim_type
                 call fatal_error(h//" Proc "//cProc// &
                      " unknown idim_type="//trim(adjustl(c0)))
              end select
           else
              sizeRes=1
           end if

           ! space for unpacked field

           allocate(Res(sizeRes), stat=ierr)
           if (ierr /= 0) then
              write(c0,"(i8)") sizeRes
              write(c1,"(i8)") ierr
              call fatal_error(h//" Proc "//cProc//" allocate Res("// &
                   trim(adjustl(c0))//") failed with stat="// &
                   trim(adjustl(c1)))
           end if

           ! gather field scaterred over slaves; master gets full unpacked field

           call OneGatherStoreAllChunks(mchnum, mynum, nmachs, &
                master_num, field, idim_type, &
                eachSizeGathered(:,idim_type), displacement(:,idim_type), &
                nnxp(ng), nnyp(ng), nnzp(ng), nzs, nzg, npatch, nwave, &
                nodemxp(:,ng), nodemyp(:,ng), nnzp(ng), &
                nodeia(:,ng), nodeiz(:,ng), nodeja(:,ng), &
                nodejz(:,ng), nodei0(:,ng), nodej0(:,ng), &
                il1(:,idim_type), ir2(:,idim_type), &
                jb1(:,idim_type), jt2(:,idim_type), &
                sizeRes, Res(:))
           
           ! master outputs field

           if (mchnum == master_num) then
              call OneFieldAnlwrt (ioaunt, sizeRes, Res, varn, idim_type, &
                   nnzp(ng), nnxp(ng), nnyp(ng), nzs, nzg, npatch, ng, &
                   npointer, nvtot, nvcnt, aw_table)
           end if

           ! dealocate space

           deallocate(scr1, stat=ierr)
           if (ierr /= 0) then
              write(c1,"(i8)") ierr
              call fatal_error(h//" Proc "//cProc//&
                   " deallocate scr1 failed with stat="//trim(adjustl(c1)))
           end if
           
           deallocate(Res, stat=ierr)
           if (ierr /= 0) then
              write(c1,"(i8)") ierr
              call fatal_error(h//" Proc "//cProc//&
                   " deallocate Res failed with stat="//trim(adjustl(c1)))
           end if

        end if

     end do

     ! master closes output file for this grid
     
     if (mchnum == master_num) then
        call CloseAnlwrt()
     end if

  end do

  ! master dumps one header file for all grids

  if (mchnum == master_num) then
     call AnlwrtFinalize(varType,anamelh,anamel,ioaunt,iclobber,nvcnt, &
          aw_table,nvtot,time,timeold)
  end if

  ! master deallocates aw_table

  if (mchnum == master_num) then
     deallocate(aw_table, stat=ierr)
     if (ierr /= 0) then
        write(c1,"(i8)") ierr
        call fatal_error(h//" Proc "//cProc//&
             " deallocate aw_table failed with stat="//trim(adjustl(c1)))
     end if
  end if

  ! deallocate WillGather

  deallocate(WillGather, stat=ierr)
  if (ierr /= 0) then
     write(c2,"(i8)") ierr
     call fatal_error(h//" deallocate WillGather failed with stat="//trim(adjustl(c2)))
  end if
end subroutine AllGatherStoreAllChunks

