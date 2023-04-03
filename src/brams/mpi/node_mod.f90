!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module node_mod

  use ModNamelistFile, only: namelistFile
  use ModParallelEnvironment, only: parallelEnvironment
  use grid_dims, only : maxgrds, maxmach
  use mem_grid,  only : ngrids,runtype
  use ModGridTree, only: GridTree, GridTreeRoot, NextOnGridTree
  use ModDomainDecomp, only: DomainDecomp

  implicit none
  private

  public :: StoreNamelistFileAtNode_mod
  public :: StoreParallelEnvironmentAtnode_mod
  public :: StoreDomainDecompAtNode_mod
  public :: DumpNodePaths

  ! MPI number of this process is mchnum
  ! MPI number of the process that does io is master_num

  integer, public :: nmachs             ! number of processes (MPI size)
  integer, public :: mchnum             ! MPI number of this process
  integer, public :: master_num         ! MPI number of the master process

  ! RAMS processes are mapped into MPI processes;
  ! RAMS processes are numbered from 1 to nmachs;
  ! RAMS process number i is mapped into MPI process machs(i);
  ! RAMS number of this process is mynum

  integer, public :: mynum              ! this process RAMS/BRAMS number; set by ProcessOrder
  integer, public, allocatable :: machs(:)     ! MPI process number, indexed by RAMS process number; set by ProcessOrder

  ! dump unit
  integer, parameter, public :: dumpUnit=21

  ! data structures for current grid (current = grid being
  ! time advanced at routine timestep) sub-domain at this process

  integer, public :: mxp
  !#x points (with ghost zone) on sub-domain
  integer, public :: myp
  !#y points (with ghost zone) on sub-domain
  integer, public :: mzp
  !#z points (with ghost zone) on sub-domain
  integer, public :: ia                 ! first x index of a sub-domain interior point
  integer, public :: iz                 ! last x index of a sub-domain interior point
  integer, public :: izu                ! last x index of a sub-domain point with interior wind x-component (u)
  integer, public :: ja                 ! first y index of a sub-domain interior point
  integer, public :: jz                 ! last y index of a sub-domain interior point
  integer, public :: jzv                ! last x index of a sub-domain point with interior wind y-component (v)
  integer, public :: ia_1               ! ia-1
  integer, public :: ia1                ! ia+1
  integer, public :: ja_1               ! ja-1
  integer, public :: ja1                ! ja+1
  integer, public :: iz_1               ! iz-1
  integer, public :: iz1                ! iz+1
  integer, public :: jz_1               ! jz-1
  integer, public :: jz1                ! jz+1
  integer, public :: i0                 ! x offset of this sub-domain in the global domain (global_x_index = local_x_index + i0)
  integer, public :: j0                 ! y offset of this sub-domain in the global domain (global_y_index = local_y_index + j0)
  integer, public :: load_bal           ! dynamic load balance flag; sould be 0 (unfinished); from RAMSIN
  integer, public :: ibcon              ! if any sub-domain boundary is a full domain boundary (bit coded)
  !                                     !  bit 1=west, bit 2=east, bit 3=south, bit 4=north


  integer, public :: newbuff_feed       ! set by StoreNodepathsBuffAlloc
  integer, public :: nbuff_feed
  integer, public :: nbuff_nest


  ! Domain Decomposition Data Structure - Part I
  !
  ! Each process stores domain decomposition info of all processes and grids
  ! on four arrays.
  ! All arrays are indexed by (rams process number, grid number).
  ! At this first part of the data structure, array values are global indices
  ! (not local indices), i.e., their range should be [1:nnxp,1:nnyp],
  ! but global indices do not include domain boundaries. Consequently,
  ! their range is limited to [2:nnxp-1,2:nnyp-1].
  ! The four arrays define sub-domain borders without ghost zone.
  ! For each grid, the set of rectangles [ixb:ixe,iyb:iye] defines
  ! a partition of grid sub-domain [2:nnxp-1,2:nnyp-1]

  ! Arrays are computed by procedure decomp_par (at para_init.f90)

  integer, allocatable, public, target :: ixb(:,:)   ! first x
  integer, allocatable, public, target :: ixe(:,:)   ! last x
  integer, allocatable, public, target :: iyb(:,:)   ! first y
  integer, allocatable, public, target :: iye(:,:)   ! last y

  ! The next set of four arrays are similar to the previous four,
  ! except that they define sub-domain borders with ghost zone.
  ! This set of four arrays do not compose a partition of the
  ! full domain, since the rectangles defined by
  ! [nxbeg:nxend, nybeg:nyend] intersect at the ghost zones and
  ! encompass domain boundaries.
  ! Again, the four arrays are indexed by (rams process number, grid number)
  ! and their values are global indices (not local indices).

  ! Arrays are computed by procedure decomp_bounds_par (at para_init.f90).

  integer, allocatable, public :: nxbeg(:,:) ! first x
  integer, allocatable, public :: nxend(:,:) ! last x
  integer, allocatable, public :: nybeg(:,:) ! first y
  integer, allocatable, public :: nyend(:,:) ! last x

  ! The auxiliar array npxy contains the total number of points
  ! at each sub-domain, including ghost zone.
  ! The array is computed by procedure decomp_bounds_par (at para_init.f90).

  integer, allocatable, public :: npxy(:,:)  ! # points

  ! Domain Decomposition Data Structure - Part II
  !
  ! These three arrays store sub-domain sizes with ghost zone.
  ! All arrays indexed by (rams process number, grid number)
  ! Array values are sizes of each sub-domain with ghost zone:
  !   (nodemxp, nodemyp, nodemzp);

  ! Arrays are computed by procedure decomp_bounds_par (at para_init.f90).

  integer, public, allocatable, target :: nodemxp(:,:) ! # x points
  integer, public, allocatable, target :: nodemyp(:,:) ! # y points
  integer, public, allocatable, target :: nodemzp(:,:) ! # z points

  ! Domain Decomposition Data Structure - Part III
  !
  ! These four arrays store first and last sub-domain position
  ! without ghost zone, that is, interior points.
  ! Each process stores domain decomposition info of all processes and grids.
  ! All arrays are indexed by (rams process number, grid number).
  ! Array values are local indices (not global indices).
  ! Sub-domain interior points are (nodeia:nodeiz, nodeja:nodejz);

  ! Arrays are computed by procedure decomp_bounds_par (at para_init.f90).

  integer, public, allocatable, target :: nodeia(:,:)  ! first x interior point
  integer, public, allocatable, target :: nodeiz(:,:)  ! last x interior point
  integer, public, allocatable, target :: nodeja(:,:)  ! first y interior point
  integer, public, allocatable, target :: nodejz(:,:)  ! last y interior point

  ! Domain Decomposition Data Structure - Part IV
  !
  ! These three arrays are auxiliaries for two situations:
  ! (1) translations of global indices to local indices and vice-versa;
  ! (2) to determine if sub-domain boundary is a global boundary or not.
  !
  ! Arrays nodei0 and nodej0 store offsets of local indices into
  ! global indices. That is, to compute a global index from a local index,
  ! add (nodei0, nodej0) to the local index.

  ! Array nodeibcon stores if sub-domain boundaries are global boundaries.
  ! Bit 1 is set iff west sub-domain boundary is a global boundary.
  ! Bit 2 is set iff east sub-domain boundary is a global boundary.
  ! Bit 3 is set iff south sub-domain boundary is a global boundary.
  ! Bit 4 is set iff north sub-domain boundary is a global boundary.

  ! Arrays are computed by procedure decomp_bounds_par (at para_init.f90).

  integer, public, allocatable, target :: nodei0(:,:)    ! x offset of a sub-domain
  integer, public, allocatable, target :: nodej0(:,:)    ! y offset of a sub-domain
  integer, public, allocatable, target :: nodeibcon(:,:) ! full domain boundary flag


  integer, allocatable, public :: nodebounds(:,:) ! Reprod.-Saulo Barros;

  ! types of messages are coded by:
  ! <type message> = 1 means communication on large timesteps to update ghost zones
  ! <type message> = 2 means communication for u on small timesteps to update ghost zones
  ! <type message> = 3 means communication for v on small timesteps to update ghost zones
  ! <type message> = 4 means communication for pi on small timesteps to update ghost zones
  ! <type message> = 5 means communication for coarser to finner grid to interpolate
  ! <type message> = 6 means communication for finner to coarser grid to interpolate
  ! <type message> = 7 means ???; not used for communication; see comment "A second index..."@node_paths_par

  ! ipaths(<info index>, <type message>, <grid>, <destination process array index>):
  !   = 0 if no message will be send to destination (for all <info index>)
  !   = global index of the first x to send to destination if <info index>=1;
  !   = global index of the last  x to send to destination if <info index>=2;
  !   = global index of the first y to send to destination if <info index>=3;
  !   = global index of the last  y to send to destination if <info index>=4;
  !   = <MPI destination process number>                   if <info index>=5; entry >= 0 on unified


  ! Allocated in NodePathsBuffAlloc, defined in node_paths_par
  public :: ipaths
  integer, allocatable :: ipaths(:,:,:,:)


  ! iget_paths(<type message>, <grid>, <source process array index>):
  !   = 0 if source process will not send message of type <type message> to this process
  !   = 1 if source process will     send message of type <type message> to this process
  ! Allocated in NodePathsBuffAlloc, defined in node_paths_par
  public :: iget_paths
  integer, allocatable :: iget_paths(:,:,:)


  ! area to store pending request ids from/to any other machine

  integer, allocatable, public :: irecv_req(:)
  integer, allocatable, public :: isend_req(:)

  ! lateral boundary conditions communication buffers

  public :: lbc_buffs
  type lbc_buffs
     real, pointer :: lbc_send_buff(:)  ! send buffer
     real, pointer :: lbc_recv_buff(:)  ! receive buffer
     integer       :: nsend             ! send buffer size
     integer       :: nrecv             ! receive buffer size
  end type lbc_buffs

  ! array of lateral boundary condition communication buffers, indexed by
  ! process array index

  type(lbc_buffs), allocatable, public :: node_buffs(:)
  type(lbc_buffs), allocatable, public :: node_buffs_st(:)

  public :: ProcessOrder

  integer, public :: f_ndmd_size   ! Stores MPI information to pack/unpack data
  integer, public :: mpi_int_size  ! Stores MPI information to pack/unpack data
  integer, public :: mpi_real_size ! Stores MPI information to pack/unpack data

  public :: alloc_paths
  public :: alloc_bounds
  public :: dealloc_bounds

  logical, parameter :: dumpLocal=.false.
contains


  !*** StoreParallelEnvironment: copy ParallelEnvironment fields to global variables

  subroutine StoreParallelEnvironmentAtNode_mod(oneParallelEnvironment)
    type(ParallelEnvironment), pointer :: oneParallelEnvironment
    character(len=16) :: fName="Dump.XXXXX.YYYYY"

    nmachs     = oneParallelEnvironment%nmachs
    mchnum     = oneParallelEnvironment%mchnum
    master_num = oneParallelEnvironment%master_num

    !**(JP)** this file open should be somewhere else
    if (dumpLocal) then
       write(fName(6:10),"(i5.5)") mchnum+1
       write(fName(12:16),"(i5.5)") nmachs
       open(dumpUnit, file=fName, action="write", status="replace")
    end if
  end subroutine StoreParallelEnvironmentAtNode_mod



  !*** ProcessOrder: order of RAMS processes and map into MPI processes ***



  subroutine ProcessOrder ()

    integer :: i, ierr
    character(len=8) :: c0, c1, c2
    character(len=*), parameter :: h="**(ProcessOrder)**"

    ! protection: bounds of array machs

    if (nmachs > maxmach) then
       write(c0,"(i8)") nmachs
       write(c1,"(i8)") maxmach
       call fatal_error(h//" selected number of processes ("//trim(adjustl(c0))//&
            ") exceeds maxmach ("//trim(adjustl(c1))//") defined at module grid_dims")
       stop
    end if

    ! map RAMS process number into MPI rank (arbitrary map)

    do i = 1, nmachs
       machs(i) = i-1
    end do

    ! this process RAMS number
    ! search is unnecessary, but allows keeping arbitrary map

    do i = 1, nmachs
       if (machs(i) == mchnum) then
          mynum = i
          exit
       end if
    end do

    ! protection: programming error

    if (i > nmachs) then
       write(c0,"(i8)") mchnum
       call fatal_error(h//" Searching for mchnum "//trim(adjustl(c0))//&
            &" was unsuccessfull")
       stop
    end if
  end subroutine ProcessOrder




  subroutine DumpOneNodePaths(grid, itype, msg)
    integer, intent(in):: grid
    integer, intent(in) :: itype
    character(len=*), intent(in) :: msg

    character(len=10) :: cProc
    character(len=7) :: cGrid
    character(len=8) :: c0, c1, c2, c3, c4, c5
    character(len=*), parameter :: h="**(DumpNodePaths)**"
    integer :: proc

    if (any(ipaths(1,itype,grid,1:nmachs) /= 0)) then
       write(dumpUnit,"(a)") h//" sends msgs with "//trim(msg)
       do proc = 1, nmachs
          if (ipaths(1,itype,grid,proc) /= 0) then
             write(c0,"(i8)") proc
             write(c1,"(i8)") ipaths(1,itype,grid,proc)
             write(c2,"(i8)") ipaths(2,itype,grid,proc)
             write(c3,"(i8)") ipaths(3,itype,grid,proc)
             write(c4,"(i8)") ipaths(4,itype,grid,proc)
             write(c5,"(i8)") ipaths(5,itype,grid,proc)
             write(dumpUnit,"(a)") h//" sends to proc "//trim(adjustl(c0))//&
                  " field horizontal portion ["//&
                  trim(adjustl(c1))//":"//trim(adjustl(c2))//","//&
                  trim(adjustl(c3))//":"//trim(adjustl(c4))//"]"
          end if
       end do
    else
       write(dumpUnit,"(a)") h//" sends no msgs with "//trim(msg)
    end if

    if (any(iget_paths(itype,grid,1:nmachs) /= 0)) then
       write(dumpUnit,"(a)") h//" recvs msgs with "//trim(msg)
       do proc = 1, nmachs
          write(c0,"(i8)") proc
          if (iget_paths(itype,grid,proc) /= 0) then
             write (dumpUnit,"(a)") h//&
                  " will receive fields from proc "//trim(adjustl(c0))
          end if
       end do
    else
       write(dumpUnit,"(a)") h//" recvs no msg with "//trim(msg)
    end if
  end subroutine DumpOneNodePaths



  subroutine DumpNodePaths(grid)
    integer, intent(in):: grid
    character(len=10) :: cProc
    character(len=7) :: cGrid
    character(len=*), parameter :: h="**(DumpNodePaths)**"
    integer :: proc, itype

    write(cProc,"(a5,i5.5)") " Proc",mynum
    write(cGrid,"(a6,i1)") " grid ",grid
    write(dumpUnit,"(a)") h//cProc//" at"//cGrid//":"

    call DumpOneNodePaths(grid, 1, "wp and pp on small ts or ghost zone on large ts")
    call DumpOneNodePaths(grid, 2, "up on small ts")
    call DumpOneNodePaths(grid, 3, "vp on small ts")
    call DumpOneNodePaths(grid, 4, "pp on small ts")
  end subroutine DumpNodePaths



  subroutine DumpUnifiedVars(ngrids)
    integer, intent(in) :: ngrids

    integer :: grid, proc
    character(len=8) :: c0, c1, c2
    character(len=256) :: line
    integer :: last
    character(len=*), parameter :: h="**(DumpUnifiedVars)**"

    write(c0,"(i8)") mchnum
    write(c1,"(i8)") nmachs
    write(c2,"(i8)") master_num
    write(*,"(a)") h//" Proc "//trim(adjustl(c0))//&
         " out of "//trim(adjustl(c1))//&
         "; master is "//trim(adjustl(c2))

    do grid = 1, ngrids
       do proc = 1, nmachs
          line = h//" Proc "; last=len_trim(line)+1
          write(line(last+1:last+4),"(i4.4)") mchnum; last=last+4
          write(line(last+1:last+10),"(a5,i4.4,a1)") "(ind ",mynum,")";last=last+10
          line(last+1:last+11) = " sees proc "; last=last+11
          write(line(last+1:last+4),"(i4.4)") proc-1; last=last+4
          line(last+1:last+9) = " at grid "; last=last+9
          write(line(last+1:last+1),"(i1)") grid; last=last+1
          line(last+1:last+7) = ": size "; last=last+7
          write(line(last+1:last+12),"(a1,i3.3,a1,i3.3,a1,i2.2,a1)")&
               "(",nodemxp(proc,grid),",",nodemyp(proc,grid),",",&
               nodemzp(proc,grid),")"; last=last+12
          line(last+1:last+8) = ", inner "; last=last+8
          write(line(last+1:last+17),"(4(a1,i3.3),a1)")&
               "(",nodeia(proc,grid),&
               ":",nodeiz(proc,grid),&
               ",",nodeja(proc,grid),&
               ":",nodejz(proc,grid),")"; last=last+17
          line(last+1:last+8) = ", start "; last=last+8
          write(line(last+1:last+9),"(2(a1,i3.3),a1)")&
               "(",nodei0(proc,grid),&
               ",",nodej0(proc,grid),")"; last=last+9
          line(last+1:last+9) = ", bounds "; last=last+9
          write(line(last+1:last+9),"(4(a1,l1),a1)") &
               "(",btest(nodeibcon(proc,grid),0),&
               ",",btest(nodeibcon(proc,grid),1),&
               ",",btest(nodeibcon(proc,grid),2),&
               ",",btest(nodeibcon(proc,grid),3),")"
          write(*,"(a)") line(1:last+9)
       end do
    end do
  end subroutine DumpUnifiedVars

  ! *******************************************************

  subroutine alloc_paths(ngrids, nmachs)

    implicit none
    ! Arguments
    integer, intent(in) :: ngrids, nmachs
    ! Local variable
    integer :: ierr

    ! Allocating ipaths and iget_paths
    allocate(ipaths(5, 7, ngrids, nmachs), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating ipaths (alloc_paths)")

    allocate(iget_paths(6, ngrids, nmachs), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating iget_paths (aloc_paths)")

  end subroutine alloc_paths

  !*************************************************************************

  subroutine alloc_bounds(ngrids, nmachs)

    implicit none
    ! Arguments
    integer, intent(in) :: ngrids, nmachs
    ! Local variable
    integer :: ierr

    ! Allocating data

    ! Allocating "machs"
    allocate(machs(nmachs), stat=ierr)
    if (ierr /= 0) call fatal_error("ERROR allocating machs (alloc_bounds)")

    allocate(nxbeg(nmachs, ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating nxbeg (alloc_bounds)")
    allocate(nxend(nmachs, ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating nxend (alloc_bounds)")
    allocate(nybeg(nmachs, ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating nybeg (alloc_bounds)")
    allocate(nyend(nmachs, ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating nyend (alloc_bounds)")
    allocate(ixb(nmachs, ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating ixb (alloc_bounds)")
    allocate(ixe(nmachs, ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating ixe (alloc_bounds)")
    allocate(iyb(nmachs, ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating iyb (alloc_bounds)")
    allocate(iye(nmachs, ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating iye (alloc_bounds)")
    allocate(npxy(nmachs, ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating npxy (alloc_bounds)")

    allocate(nodemxp(nmachs, ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating nodemxp (alloc_bounds)")
    allocate(nodemyp(nmachs, ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating nodemyp (alloc_bounds)")
    allocate(nodemzp(nmachs, ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating nodemzp (alloc_bounds)")
    allocate(nodeia(nmachs, ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating nodeia (alloc_bounds)")
    allocate(nodeiz(nmachs, ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating nodeiz (alloc_bounds)")
    allocate(nodeja(nmachs, ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating nodeja (alloc_bounds)")
    allocate(nodejz(nmachs, ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating nodejz (alloc_bounds)")
    allocate(nodei0(nmachs, ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating nodei0 (alloc_bounds)")
    allocate(nodej0(nmachs, ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating nodej0 (alloc_bounds)")
    allocate(nodeibcon(nmachs, ngrids), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating nodeibcon (alloc_bounds)")

    allocate(irecv_req(nmachs), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating irecv_req (alloc_bounds)")
    allocate(isend_req(nmachs), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating isend_req (alloc_bounds)")

    allocate(node_buffs(nmachs), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating node_buffs (alloc_bounds)")
    allocate(node_buffs_st(nmachs), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating node_buffs_st (alloc_bounds)")

    ! maxmach,maxgrds
    allocate(nodebounds(maxgrds,8), STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR allocating nodebounds (alloc_bounds)")

  end subroutine alloc_bounds

  ! *************************************************************************

  subroutine dealloc_bounds()

    implicit none
    ! Local variable
    integer :: ierr

!--(inspxe)--------------------------------------------------------------------------------------------------
! Correcao bug: Memory leak
    integer :: nm
!--(inspxe)--------------------------------------------------------------------------------------------------

    ! Allocating data
    deallocate(machs, stat=ierr)
    if (ierr /= 0) call fatal_error("ERROR deallocating machs (dealloc_bounds)")

    deallocate(nxbeg, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating nxbeg (dealloc_bounds)")
    deallocate(nxend, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating nxend (dealloc_bounds)")
    deallocate(nybeg, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating nybeg (dealloc_bounds)")
    deallocate(nyend, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating nyend (dealloc_bounds)")
    deallocate(ixb, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating ixb (dealloc_bounds)")
    deallocate(ixe, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating ixe (dealloc_bounds)")
    deallocate(iyb, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating iyb (dealloc_bounds)")
    deallocate(iye, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating iye (dealloc_bounds)")
    deallocate(npxy, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating npxy (dealloc_bounds)")

    deallocate(nodemxp, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating nodemxp (dealloc_bounds)")
    deallocate(nodemyp, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating nodemyp (dealloc_bounds)")
    deallocate(nodemzp, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating nodemzp (dealloc_bounds)")
    deallocate(nodeia, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating nodeia (dealloc_bounds)")
    deallocate(nodeiz, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating nodeiz (dealloc_bounds)")
    deallocate(nodeja, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating nodeja (dealloc_bounds)")
    deallocate(nodejz, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating nodejz (dealloc_bounds)")
    deallocate(nodei0, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating nodei0 (dealloc_bounds)")
    deallocate(nodej0, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating nodej0 (dealloc_bounds)")
    deallocate(nodeibcon, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating nodeibcon (dealloc_bounds)")

    deallocate(nodebounds, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating nodebounds (dealloc_bounds)")

    deallocate(irecv_req, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating irecv_req (dealloc_bounds)")

    deallocate(isend_req, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating isend_req (dealloc_bounds)")

!--(inspxe)--------------------------------------------------------------------------------------------------
! Correcao bug: Memory leak
    if (.not. (runtype(1:7)=='MAKESFC' .or. runtype(1:9)=='MAKEVFILE')) then
      do nm=1,nmachs
        deallocate(node_buffs(nm)%lbc_send_buff, STAT=ierr)
        if (ierr/=0) call fatal_error("ERROR deallocating node_buffs%lbc_send_buff (dealloc_bounds)")
        deallocate(node_buffs(nm)%lbc_recv_buff, STAT=ierr)
        if (ierr/=0) call fatal_error("ERROR deallocating node_buffs%lbc_recv_buff (dealloc_bounds)")
        deallocate(node_buffs_st(nm)%lbc_send_buff, STAT=ierr)
        if (ierr/=0) call fatal_error("ERROR deallocating node_buffs_st%lbc_send_buff (dealloc_bounds)")
        deallocate(node_buffs_st(nm)%lbc_recv_buff, STAT=ierr)
        if (ierr/=0) call fatal_error("ERROR deallocating node_buffs_st%lbc_recv_buff (dealloc_bounds)")
      enddo
    endif
!--(inspxe)--------------------------------------------------------------------------------------------------

    deallocate(node_buffs, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating node_buffs (dealloc_bounds)")
    deallocate(node_buffs_st, STAT=ierr)
    if (ierr/=0) call fatal_error("ERROR deallocating node_buffs_st (dealloc_bounds)")

  end subroutine dealloc_bounds

  subroutine StoreNamelistFileAtNode_mod(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile

    load_bal = oneNamelistFile%load_bal
  end subroutine StoreNamelistFileAtNode_mod




  subroutine StoreDomainDecompAtNode_mod(AllGrids)
    type(GridTree), pointer :: AllGrids

    integer :: gridID, node
    type(GridTree), pointer :: OneGridTree => null()
    type(DomainDecomp), pointer :: GlobalNoGhost => null()

    OneGridTree => GridTreeRoot(AllGrids)
    do while (associated(OneGridTree))
       gridId = OneGridTree%curr%Id
       GlobalNoGhost => OneGridTree%curr%GlobalNoGhost
       do node = 1, OneGridTree%curr%ParEnv%nmachs
          ixb(node,gridId) = GlobalNoGhost%xb(node)
          ixe(node,gridId) = GlobalNoGhost%xe(node)
          iyb(node,gridId) = GlobalNoGhost%yb(node)
          iye(node,gridId) = GlobalNoGhost%ye(node)
       end do
       OneGridTree => NextOnGridTree(OneGridTree)
    end do
  end subroutine StoreDomainDecompAtNode_mod
end module node_mod
