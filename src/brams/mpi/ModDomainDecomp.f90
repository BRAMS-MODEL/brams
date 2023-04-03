module ModDomainDecomp
  use ModParallelEnvironment, only: &
       ParallelEnvironment, MsgDump
  use ModGridDims, only: GridDims

  implicit none
  private
  public :: DomainDecomp
  public :: CreateGlobalNoGhost
  public :: CreateGlobalWithGhost
  public :: CreateLocalInterior
  public :: DestroyDomainDecomp
  public :: DumpDomainDecomp
  public :: DumpDomainDecompHistogram


  ! DomainDecomp: stores indices of a domain decomposed grid.
  !               at node i, indices of the sub-domain stored at
  !               this node are are [xb(i):xe(i), yb(i),ye(i)] and
  !               ibcon(i) stores if any sub-domain boundary is
  !               a full domain boundary.
  !               GhostZoneLength is the length of the ghost zone, if any.
  !
  !               There are three variables of this type:
  !               GlobalNoGhost: partition of domain [2:nnxp,2:nnyp];
  !                              indices are global indices;
  !                              no ghost zone.
  !               GlobalWithGhost: GlobalNoGhost with ghost zone;
  !                                not a real partition, since sub-domains instersect
  !                                indices are global indices;
  !                                ghost zone is valid and included.
  !               LocalInterior: Local indices of interior columns of GlobalNoGhost;
  !                              Stores required indices for dynamics and physics;
  !                              xb and xe are the local indices converted from
  !                              GlobalNoGhost%xb and GlobalNoGhost%xe; same for yb, ye.
  !                              nx and ny are total sub-domain size (same as GlobalWithGhost)

  type DomainDecomp

     ! GhostZoneLength: 0 at GlobalNoGhost; 
     !                  exact value at GlobalWithGhost;
     !                  exact value at LocalInterior

     integer :: GhostZoneLength

     ! x axis indexed from xb to xe with nx indices;
     ! y axis indexed from yb to ye with ny indices;
     ! global indices at GlobalNoGhost and GlobalWithGhost;
     ! local indices at LocalInterior;
     ! local x index = global x index - (GlobalWithGhost%xb-1)
     ! local y index = global y index - (GlobalWithGhost%yb-1)

     integer, allocatable :: xb(:)   ! first x
     integer, allocatable :: xe(:)   ! last x 
     integer, allocatable :: nx(:)   ! # points x
     integer, allocatable :: yb(:)   ! first y
     integer, allocatable :: ye(:)   ! last y 
     integer, allocatable :: ny(:)   ! # points y

     ! Array ibcon (boundary condition), also indexed by rams process number, stores
     ! if sub-domain boundary is a global boundary or not. Info is coded
     ! on a bit structure:
     ! Bit 1 is set iff west sub-domain boundary (low x value) is a global boundary
     ! Bit 2 is set iff east sub-domain boundary (high x value) is a global boundary
     ! Bit 3 is set iff south sub-domain boundary (low y value) is a global boundary
     ! Bit 4 is set iff north sub-domain boundary (high y value) is a global boundary

     ! To check if a sub-domain boundary is a global domain boundary, use:
     !  if btest(ibcon(mach),4) is true, sub-domain north boundary is a global boundary
     !  if btest(ibcon(mach),3) is true, sub-domain south boundary is a global boundary
     !  if btest(ibcon(mach),2) is true, sub-domain east boundary is a global boundary
     !  if btest(ibcon(mach),1) is true, sub-domain west boundary is a global boundary

     integer, allocatable :: ibcon(:) ! full domain boundary flag

  end type DomainDecomp

  logical, parameter :: dumpLocal = .false.
contains



  ! CreateGlobalNoGhost: Creates a variable of type DomainDecomp for
  !                            given grid and parallel environment. Performs domain
  !                            decomposition, filling all components of the created
  !                            variable.



  subroutine CreateGlobalNoGhost(GridSize, ParEnv, GlobalNoGhost)
    type(ParallelEnvironment), pointer :: ParEnv 
    type(GridDims), pointer :: GridSize
    type(DomainDecomp), pointer :: GlobalNoGhost

    character(len=8) :: c0, c1, c2
    character(len=*), parameter :: h="**(CreateGlobalNoGhost)**"
    logical, parameter :: dumpLocal=.false.

    integer :: nxp
    integer :: nyp
    integer :: nmachs

    if (.not. associated(GridSize)) then
       call fatal_error(h//" invoked with null GridSize")
    else if (.not. associated(ParEnv)) then
       call fatal_error(h//" invoked with null ParEnv")
    else if (associated(GlobalNoGhost)) then
       call fatal_error(h//" invoked with GlobalNoGhost "//&
            "already created; invoke DestroyDomainDecomp")
    end if

    nxp = GridSize%nnxp
    nyp = GridSize%nnyp
    nmachs = ParEnv%nmachs

    if (dumpLocal) then
       write(c0,"(i8)") nxp
       write(c1,"(i8)") nyp
       write(c2,"(i8)") nmachs
       call MsgDump (h//" starts with domain ["//&
            trim(adjustl(c0))//" x "//trim(adjustl(c1))//&
            "] to decompose into "//trim(adjustl(c2))//" sub-domains")
    end if

    allocate(GlobalNoGhost)
    allocate(GlobalNoGhost%xb(nmachs))
    allocate(GlobalNoGhost%xe(nmachs))
    allocate(GlobalNoGhost%nx(nmachs))
    allocate(GlobalNoGhost%yb(nmachs))
    allocate(GlobalNoGhost%ye(nmachs))
    allocate(GlobalNoGhost%ny(nmachs))
    allocate(GlobalNoGhost%ibcon(nmachs))

    ! no ghost zone info

    GlobalNoGhost%GhostZoneLength = 0

    ! find Rams domain decomposition xb, xe, yb, ye

    call DomainDecompRams(nxp, nyp, nmachs, &
         GlobalNoGhost%xb, GlobalNoGhost%xe, &
         GlobalNoGhost%yb, GlobalNoGhost%ye)

    ! fill boundary condition ibcon

    call MarkBoundary(nxp, nyp, nmachs, &
         GlobalNoGhost%xb, GlobalNoGhost%xe, &
         GlobalNoGhost%yb, GlobalNoGhost%ye, &
         GlobalNoGhost%ibcon)

    ! fill number of points at each sub-domain 

    GlobalNoGhost%nx = GlobalNoGhost%xe - GlobalNoGhost%xb + 1
    GlobalNoGhost%ny = GlobalNoGhost%ye - GlobalNoGhost%yb + 1

    ! Verify partition correctness

    call CheckPartition(nxp, nyp, nmachs, &
         GlobalNoGhost%xb, GlobalNoGhost%xe, &
         GlobalNoGhost%yb, GlobalNoGhost%ye, &
         GlobalNoGhost%ibcon)

  end subroutine CreateGlobalNoGhost



  ! DomainDecompRams: original RAMS domain decomposition




  subroutine DomainDecompRams(nxp,nyp,nodes,&
       xb,xe,yb,ye)
    ! Arguments:
    integer, intent(in)  :: nxp
    integer, intent(in)  :: nyp
    integer, intent(in)  :: nodes
    integer, intent(out) :: xb(nodes)
    integer, intent(out) :: xe(nodes)
    integer, intent(out) :: yb(nodes)
    integer, intent(out) :: ye(nodes)

    ! Local variables:
    real  :: work(nxp,nyp)
    real :: workrow(nyp)
    real :: workcol(nxp)
    real :: workload(nodes)
    real :: workblock(nodes)
    integer :: jrows(nodes)
    integer :: jrow(nodes)
    integer :: nblocks(nodes)
    real :: relspeed(nodes)
    real :: bfact
    character(len=*), parameter :: h="**(DomainDecompRams)**"

    integer :: inode,i,j,islab,jnodes,nslabs,min_blocks,nbigslabs,iblock &
         ,jnode,knode
    real :: anodes,aslabs,totspeed,workdom,workaccum,worksofar &
         ,slabspeed,workslab

    ! default relspeed = 1.0 for nodes of uniform speed.

    relspeed = 1.0

    ! work factor of 1.0 at interior points and .2 at boundary points

    work = 1.0
    bfact=.2
    do j = 1,nyp
       work(1,j) = bfact
       work(nxp,j) = bfact
    enddo
    do i = 1,nxp
       work(i,1) = bfact
       work(i,nyp) = bfact
    enddo

    ! This routine decomposes grid domains of size (nxp,nyp) into a number,
    ! specified by nodes, of rectangular subdomains.  The convention is followed
    ! that any internal boundaries (between subdomains) that are parallel to
    ! the x-axis run continuously across the full domain, while boundaries
    ! parallel to the y-axis may or may not run the full distance across the
    ! domain.  For convenience, regions of the domain bounded by adjacent
    ! east-west internal boundaries are termed "slabs", while smaller divisions
    ! within each slab are termed "blocks".  Each block is required to have
    ! a minimum dimension of 6 by 6 grid cells.  If this cannot be satisfied
    ! with the given input parameters, the subroutine stops.


    ! Estimate the number of slabs to be used (aslabs), and compute a final
    ! nearest integer value (nslabs) which is limited to allowable values.
    ! Zero out array for accumulating number of columns for each node.

    anodes = float(nodes)
    aslabs = sqrt(anodes * float(nyp) / float(nxp))
    nslabs = min(nodes,max(1,nint(aslabs)))

    totspeed = 0.
    do inode = 1,nodes
       xe(inode) = 0
       totspeed = totspeed + relspeed(inode)
    enddo

    ! Compute total work load over each row and over entire domain.

    workdom = 0.
    do j = 1,nyp
       workrow(j) = 0.
       do i = 1,nxp
          workrow(j) = workrow(j) + work(i,j)
       enddo
       workdom = workdom + workrow(j)
    enddo
    workrow(2) = workrow(2) + workrow(1)
    workrow(nyp-1) = workrow(nyp-1) + workrow(nyp)

    ! Determine number of blocks and the average workload for each slab.

    min_blocks = nodes / nslabs
    nbigslabs = nodes - min_blocks * nslabs
    inode = 0
    do islab = 1,nslabs
       workload(islab) = 0.
       nblocks(islab) = min_blocks
       if (islab .le. nbigslabs) nblocks(islab) = min_blocks + 1
       do iblock = 1,nblocks(islab)
          inode = inode + 1
          workload(islab) = workload(islab)  &
               + workdom * relspeed(inode) / totspeed
       enddo
    enddo

    ! Assign all j-rows to their respective slabs in a way that balances the work
    ! load among slabs according to their respective numbers of nodes (blocks).
    ! The array jrows counts the number of rows in each slab, and the array
    ! jrow is the index of the southernmost row in each slab.

    do islab = 1,nslabs
       jrows(islab) = 0
    enddo

    workaccum = 0.
    worksofar = 0.
    islab = 0

    do j = 2,nyp-1
       workaccum = workaccum + workrow(j)
       if (workaccum - .5 * workrow(j) .gt. worksofar .and.  &
            islab .lt. nslabs) then
          islab = islab + 1
          jrow(islab) = j
          worksofar = worksofar + workload(islab)
       endif
       jrows(islab) = jrows(islab) + 1
    enddo

    inode = 0
    jnode = 0
    knode = 0
    do islab = 1,nslabs

       ! Compute the total work load for each slab and for each i-column in the
       ! slab.

       slabspeed = 0.
       workslab = 0.
       do i = 1,nxp
          workcol(i) = 0.
          do j = jrow(islab),jrow(islab)+jrows(islab)-1
             workcol(i) = workcol(i) + work(i,j)
          enddo
          workslab = workslab + workcol(i)
       enddo
       workcol(2) = workcol(2) + workcol(1)
       workcol(nxp-1) = workcol(nxp-1) + workcol(nxp)

       ! Determine average workload for each block.

       do iblock = 1,nblocks(islab)
          jnode = jnode + 1
          slabspeed = slabspeed + relspeed(jnode)
       enddo
       do iblock = 1,nblocks(islab)
          knode = knode + 1
          workblock(iblock) = workslab  &
               * relspeed(knode) / slabspeed
       enddo

       ! Assign the i-columns of each slab to their respective blocks in a way that
       ! balances the work load among the blocks.  The array ncols counts the number
       ! of i-columns on each node, and the array ncol is the index of the
       ! westernmost i-column on each node.

       workaccum = 0.
       worksofar = 0.

       iblock = 0
       do i = 2,nxp-1
          workaccum = workaccum + workcol(i)
          if (workaccum - .5 * workcol(i) .gt. worksofar .and.  &
               iblock .lt. nblocks(islab)) then
             iblock = iblock + 1
             inode = inode + 1
             yb(inode) = jrow(islab)
             xb(inode) = i
             ye(inode) = yb(inode) + jrows(islab) - 1
             worksofar = worksofar + workblock(iblock)
          endif
          xe(inode) = xe(inode) + 1
       enddo
    enddo

    do jnode = 1,nodes
       xe(jnode) = xb(jnode) + xe(jnode) - 1
    enddo

    ! Check to make sure that each subdomain has at least 2 interior 
    ! rows and columns.

    do jnode = 1,nodes
       if (ye(jnode) - yb(jnode) .lt. 1 .or.  &
            xe(jnode) - xb(jnode) .lt. 1) then
          call fatal_error(h//" too many processors")
          return
       endif
    enddo
  end subroutine DomainDecompRams




  ! CheckPartition: verifies if xb, xe, yb, ye
  !                 partitions domain [2:nxp-1,2:nyp-1] across nodes



  subroutine CheckPartition(nxp, nyp, nodes, &
       xb, xe, yb, ye, ibcon)
    integer, intent(in)  :: nxp
    integer, intent(in)  :: nyp
    integer, intent(in)  :: nodes
    integer, intent(in)  :: xb(:)
    integer, intent(in)  :: xe(:)
    integer, intent(in)  :: yb(:)
    integer, intent(in)  :: ye(:)
    integer, intent(in) :: ibcon(:)

    integer :: node
    integer :: ix, xstart, xend
    integer :: iy, ystart, yend
    integer, parameter :: UNASSIGNED=-1
    integer :: owner(nxp,nyp)
    character(len=8) :: c0, c1, c2, c3
    character(len=*), parameter :: h="**(CheckPartition)**"

    owner = UNASSIGNED
    do node = 1, nodes
       if (btest(ibcon(node),1)) then
          xstart = xb(node)-1
       else
          xstart = xb(node)
       end if

       if (btest(ibcon(node),2)) then
          xend = xe(node)+1
       else
          xend = xe(node)
       end if

       if (btest(ibcon(node),3)) then
          ystart = yb(node)-1
       else
          ystart = yb(node)
       end if

       if (btest(ibcon(node),4)) then
          yend = ye(node)+1
       else
          yend = ye(node)
       end if

       do iy = ystart, yend
          do ix = xstart, xend
             if (owner(ix,iy) == UNASSIGNED) then
                owner(ix,iy) = node
             else
                write(c0,"(i8)") ix
                write(c1,"(i8)") iy
                write(c2,"(i8)") node
                write(c3,"(i8)") owner(ix,iy)
                call fatal_error(h//"  point "//&
                     "("//trim(adjustl(c1))//","//trim(adjustl(c2))//")"//&
                     " is assigned to node "//trim(adjustl(c2))//&
                     " and to node "//trim(adjustl(c3)))
             end if
          end do
       end do
    end do
    if (any(owner == UNASSIGNED)) then
       call fatal_error(h//" there are unassigned domain points")
    end if
  end subroutine CheckPartition



  ! MarkBoundary: verify if sub-domain boundaries are global domain boundaries
  !               or not; fills array ibcon



  subroutine MarkBoundary(nxp, nyp, nodes, &
       xb, xe, yb, ye, ibcon)
    integer, intent(in)  :: nxp
    integer, intent(in)  :: nyp
    integer, intent(in)  :: nodes
    integer, intent(in)  :: xb(:)
    integer, intent(in)  :: xe(:)
    integer, intent(in)  :: yb(:)
    integer, intent(in)  :: ye(:)
    integer, intent(out) :: ibcon(:)

    integer :: node
    character(len=*), parameter :: h="**(MarkBoundary)**"

    ! position nodes into domain points [1:nxp,1:nyp] making sure
    ! that nodes do not intersect

    do node = 1, nodes
       ibcon(node) = 0
       if (xb(node) == 2) then
          ibcon(node)=ibset(ibcon(node),1)
       end if
       if (xe(node) == nxp-1) then
          ibcon(node)=ibset(ibcon(node),2)
       end if
       if (yb(node) == 2) then
          ibcon(node)=ibset(ibcon(node),3)
       end if
       if (ye(node) == nyp-1) then
          ibcon(node)=ibset(ibcon(node),4)
       end if
    end do
  end subroutine MarkBoundary




  ! CreateGlobalWithGhost: Creates a variable of type DomainDecomp that extends
  !                        the sub-domains of GlobalNoGhost to include ghost zone
  !                        of length GhostZoneLength. Sub-domain [xb:xe] is limited
  !                        by [1:nxp], same with [yb:ye]. Number of points nx, ny
  !                        include ghost zone. Boundary conditions ibcon are copied
  !                        from GlobalWithGhost.



  subroutine CreateGlobalWithGhost(GridSize, ParEnv, GhostZoneLength, &
       GlobalNoGhost, GlobalWithGhost)
    type(GridDims), pointer :: GridSize              ! intent(in)
    type(ParallelEnvironment), pointer :: ParEnv     ! intent(in)
    integer, intent(in) :: GhostZoneLength
    type(domainDecomp), pointer :: GlobalNoGhost     ! intent(in)
    type(domainDecomp), pointer :: GlobalWithGhost   ! intent(out)

    integer :: nxp
    integer :: nyp
    integer :: nmachs
    integer :: cell
    character(len=8) :: c0, c1, c2
    character(len=*), parameter :: h="**(CreateGlobalWithGhost)**"

    if (.not. associated(GridSize)) then
       call fatal_error(h//" invoked with null GridSize")
    else if (.not. associated(ParEnv)) then
       call fatal_error(h//" invoked with null ParEnv")
    else if (.not. associated(GlobalNoGhost)) then
       call fatal_error(h//" invoked with null GlobalNoGhost ")
    else if (associated(GlobalWithGhost)) then
       call fatal_error(h//" invoked with GlobalWithGhost "//&
            "already created; invoke DestroyDomainDecomp")
    else if (GhostZoneLength < 1) then
       write (c0,"(i8)") GhostZoneLength
       call fatal_error(h//" GhostZoneLength sould be >= 1 but it is "//&
            trim(adjustl(c0)))
    end if


    nmachs = ParEnv%nmachs
    nxp = GridSize%nnxp
    nyp = GridSize%nnyp

    if (dumpLocal) then
       write(c0,"(i8)") nxp
       write(c1,"(i8)") nyp
       write(c2,"(i8)") GhostZoneLength
       call MsgDump (h//" extends the domain decomposed GlobalNoGhost ["//&
            trim(adjustl(c0))//" x "//trim(adjustl(c1))//&
            "] with ghost zone of length "//trim(adjustl(c2)))
    end if

    allocate(GlobalWithGhost)
    allocate(GlobalWithGhost%xb(nmachs))
    allocate(GlobalWithGhost%xe(nmachs))
    allocate(GlobalWithGhost%nx(nmachs))
    allocate(GlobalWithGhost%yb(nmachs))
    allocate(GlobalWithGhost%ye(nmachs))
    allocate(GlobalWithGhost%ny(nmachs))
    allocate(GlobalWithGhost%ibcon(nmachs))

    ! 

    GlobalWithGhost%GhostZoneLength = GhostZoneLength

    do cell = 1, nmachs
       GlobalWithGhost%ibcon(cell) = GlobalNoGhost%ibcon(cell)
       if (btest(GlobalWithGhost%ibcon(cell),1)) then
          GlobalWithGhost%xb(cell) = 1
       else
          GlobalWithGhost%xb(cell) = max(1,GlobalNoGhost%xb(cell) - GhostZoneLength)
       end if
       if (btest(GlobalWithGhost%ibcon(cell),2)) then
          GlobalWithGhost%xe(cell) = nxp
       else
          GlobalWithGhost%xe(cell) = min(nxp,GlobalNoGhost%xe(cell) + GhostZoneLength)
       end if
       if (btest(GlobalWithGhost%ibcon(cell),3)) then
          GlobalWithGhost%yb(cell) = 1
       else
          GlobalWithGhost%yb(cell) = max(1,GlobalNoGhost%yb(cell) - GhostZoneLength)
       end if
       if (btest(GlobalWithGhost%ibcon(cell),4)) then
          GlobalWithGhost%ye(cell) = nyp
       else
          GlobalWithGhost%ye(cell) = min(nyp,GlobalNoGhost%ye(cell) + GhostZoneLength)
       end if
       GlobalWithGhost%nx(cell) = GlobalWithGhost%xe(cell)-GlobalWithGhost%xb(cell)+1
       GlobalWithGhost%ny(cell) = GlobalWithGhost%ye(cell)-GlobalWithGhost%yb(cell)+1
    end do

  end subroutine CreateGlobalWithGhost



  ! CreateLocalInterior: Computes the interior points of GlobalWithGhost.
  !                      This variable stores at [xb:xe,yb:ye] the local indices
  !                      of interior points that should be computed at most physics
  !                      and dynamics routines. Variables nx, ny should be used to
  !                      dimension fields, since they include ghost zones. 



  subroutine CreateLocalInterior(ParEnv, GlobalWithGhost, &
       GlobalNoGhost, LocalInterior)
    type(ParallelEnvironment), pointer :: ParEnv   ! intent(in)
    type(domainDecomp), pointer :: GlobalWithGhost ! intent(in)
    type(domainDecomp), pointer :: GlobalNoGhost   ! intent(in)
    type(domainDecomp), pointer :: LocalInterior   ! intent(out)

    integer :: nmachs
    integer :: cell
    integer :: x0
    integer :: y0
    character(len=*), parameter :: h="**(CreateLocalInterior)**"

    if (.not. associated(ParEnv)) then
       call fatal_error(h//" invoked with null ParEnv")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" invoked with null GlobalWithGhost ")
    else if (.not. associated(GlobalNoGhost)) then
       call fatal_error(h//" invoked with null GlobalNoGhost ")
    else if (associated(LocalInterior)) then
       call fatal_error(h//" invoked with LocalInterior "//&
            "already created; invoke DestroyDomainDecomp")
    end if

    if (dumpLocal) then
       call MsgDump (h//" generates local indices from GlobalWithGhost")
    end if

    nmachs=ParEnv%nmachs
    allocate(LocalInterior)
    allocate(LocalInterior%xb(nmachs))
    allocate(LocalInterior%xe(nmachs))
    allocate(LocalInterior%nx(nmachs))
    allocate(LocalInterior%yb(nmachs))
    allocate(LocalInterior%ye(nmachs))
    allocate(LocalInterior%ny(nmachs))
    allocate(LocalInterior%ibcon(nmachs))


    LocalInterior%GhostZoneLength = GlobalWithGhost%GhostZoneLength 
    do cell = 1, nmachs
       x0 = GlobalWithGhost%xb(cell)-1
       LocalInterior%xb(cell) = GlobalNoGhost%xb(cell) - x0
       LocalInterior%xe(cell) = GlobalNoGhost%xe(cell) - x0
       LocalInterior%nx(cell) = GlobalWithGhost%nx(cell)
       y0 = GlobalWithGhost%yb(cell)-1
       LocalInterior%yb(cell) = GlobalNoGhost%yb(cell) - y0
       LocalInterior%ye(cell) = GlobalNoGhost%ye(cell) - y0
       LocalInterior%ny(cell) = GlobalWithGhost%ny(cell)
       LocalInterior%ibcon(cell) = GlobalNoGhost%ibcon(cell)
    end do

  end subroutine CreateLocalInterior



  ! DestroyDomainDecomp: Removes all storage allocated to a variable 
  !                             of type DomainDecomp, if any.



  subroutine DestroyDomainDecomp(OneDomainDecomp)
    type(domainDecomp), pointer :: OneDomainDecomp

    character(len=*), parameter :: h="**(DestroyDomainDecomp)**"

    if (associated(oneDomainDecomp)) then
       deallocate(oneDomainDecomp%xb)
       deallocate(oneDomainDecomp%xe)
       deallocate(oneDomainDecomp%nx)
       deallocate(oneDomainDecomp%yb)
       deallocate(oneDomainDecomp%ye)
       deallocate(oneDomainDecomp%ny)
       deallocate(oneDomainDecomp%ibcon)
       deallocate(oneDomainDecomp)
    end if
    nullify(oneDomainDecomp)
  end subroutine DestroyDomainDecomp



  ! DumpDomainDecomp: Dumps info stored at a variable 
  !                          of type DomainDecomp, if any.



  subroutine DumpDomainDecomp(oneDomainDecomp, name)
    type(DomainDecomp), pointer :: oneDomainDecomp
    character(len=*), intent(in) :: name

    integer :: nmachs
    integer :: mach
    character(len=56) :: msg
    character(len=8) :: bnd
    character(len=*), parameter :: h="**(DomainDecompDump)**"
    character(len=8) :: c0, c1, c2

    ! dumps at selected unit

    if (.not. associated(oneDomainDecomp)) then
       call MsgDump (h//" empty DomainDecomp named "//trim(name))
       return
    else if (.not. allocated(oneDomainDecomp%xb)) then
       call fatal_error(h//" oneDomainDecomp%xb not allocated")
    else if (.not. allocated(oneDomainDecomp%xe)) then
       call fatal_error(h//" oneDomainDecomp%xe not allocated")
    else if (.not. allocated(oneDomainDecomp%nx)) then
       call fatal_error(h//" oneDomainDecomp%nx not allocated")
    else if (.not. allocated(oneDomainDecomp%yb)) then
       call fatal_error(h//" oneDomainDecomp%yb not allocated")
    else if (.not. allocated(oneDomainDecomp%ye)) then
       call fatal_error(h//" oneDomainDecomp%ye not allocated")
    else if (.not. allocated(oneDomainDecomp%ny)) then
       call fatal_error(h//" oneDomainDecomp%ny not allocated")
    else if (.not. allocated(oneDomainDecomp%ibcon)) then
       call fatal_error(h//" oneDomainDecomp%ibcon not allocated")
    end if

    nmachs = size(oneDomainDecomp%xb)

    write(c0,"(i8)") nmachs
    write(c1,"(i8)") oneDomainDecomp%GhostZoneLength
    call MsgDump (h//" named "//trim(name)//" has "//&
         trim(adjustl(c0))//" sub-domains and ghost zone of length "//&
         trim(adjustl(c1)))

    call MsgDump('    node   x-beg   x-end   y-beg   y-end    cols   ibcon')
    do mach = 1,nmachs
       bnd=""
       if (btest(oneDomainDecomp%ibcon(mach),1)) bnd=trim(bnd)//"X-"
       if (btest(oneDomainDecomp%ibcon(mach),2)) bnd=trim(bnd)//"X+"
       if (btest(oneDomainDecomp%ibcon(mach),3)) bnd=trim(bnd)//"Y-"
       if (btest(oneDomainDecomp%ibcon(mach),4)) bnd=trim(bnd)//"Y+"
       write(c0,"(i8)") mach
       msg=c0
       write(c0,"(i8)") oneDomainDecomp%xb(mach)
       msg=trim(msg)//c0
       write(c0,"(i8)") oneDomainDecomp%xe(mach)
       msg=trim(msg)//c0
       write(c0,"(i8)") oneDomainDecomp%yb(mach)
       msg=trim(msg)//c0
       write(c0,"(i8)") oneDomainDecomp%ye(mach)
       msg=trim(msg)//c0
       write(c0,"(i8)") oneDomainDecomp%nx(mach)*oneDomainDecomp%ny(mach)
       msg=trim(msg)//c0
       msg=trim(msg)//adjustr(bnd)
       call MsgDump (msg)
    end do
  end subroutine DumpDomainDecomp



  ! DumpDomainDecompHistogram: Dumps an histogram of sub-domains
  !                                   at a variable of type DomainDecomp,
  !                                   if any.



  subroutine DumpDomainDecompHistogram(oneDomainDecomp, name)
    type(DomainDecomp), pointer :: oneDomainDecomp
    character(len=*), intent(in) :: name

    integer :: nmachs
    integer :: mach
    integer, allocatable :: ncols(:)
    integer, allocatable :: bins(:) 
    integer :: low, high, ibin
    character(len=8) :: c0, c1, c2, c3
    character(len=*), parameter :: h="**(DumpDomainDecompHistogram)**"
    logical, parameter :: dumpLocal=.false.

    ! dumps at selected unit

    if (.not. associated(oneDomainDecomp)) then
       call fatal_error(h//" null oneDomainDecomp")
    else if (.not. allocated(oneDomainDecomp%xb)) then
       call fatal_error(h//" oneDomainDecomp%xb not allocated")
    else if (.not. allocated(oneDomainDecomp%xe)) then
       call fatal_error(h//" oneDomainDecomp%xe not allocated")
    else if (.not. allocated(oneDomainDecomp%nx)) then
       call fatal_error(h//" oneDomainDecomp%nx not allocated")
    else if (.not. allocated(oneDomainDecomp%yb)) then
       call fatal_error(h//" oneDomainDecomp%yb not allocated")
    else if (.not. allocated(oneDomainDecomp%ye)) then
       call fatal_error(h//" oneDomainDecomp%ye not allocated")
    else if (.not. allocated(oneDomainDecomp%ny)) then
       call fatal_error(h//" oneDomainDecomp%ny not allocated")
    else if (.not. allocated(oneDomainDecomp%ibcon)) then
       call fatal_error(h//" oneDomainDecomp%ibcon not allocated")
    end if

    ! number of ranks

    nmachs = size(oneDomainDecomp%xb)

    ! dumps

    write(c0,"(i8)") nmachs
    write(c1,"(i8)") oneDomainDecomp%GhostZoneLength
    call MsgDump (h//" named "//trim(name)//" has "//&
         trim(adjustl(c0))//" sub-domains and ghost zone of length "//&
         trim(adjustl(c1)))

    ! grid points per rank

    allocate(ncols(nmachs))
    do mach = 1, nmachs
       ncols(mach)= &
            (1+oneDomainDecomp%xe(mach)-oneDomainDecomp%xb(mach))  &
            *(1+oneDomainDecomp%ye(mach)-oneDomainDecomp%yb(mach))
    end do

    ! build histogram

    low=minval(ncols(:))
    high=maxval(ncols(:))

    allocate(bins(low:high))
    bins = 0
    do mach=1, nmachs
       bins(ncols(mach)) = bins(ncols(mach)) + 1
    end do

    ! dump histogram

    do ibin = low, high
       if (bins(ibin) == 0) cycle
       write(c0,"(i8)") bins(ibin)
       write(c1,"(i8)") ibin
       write(c2,"(i8)") (100*bins(ibin))/nmachs
       call MsgDump(trim(adjustl(c0))//" ranks with "//&
            trim(adjustl(c1))//" surface points; "//trim(adjustl(c2))//"% of ranks")
    end do
    
    deallocate(ncols)
    deallocate(bins)
  end subroutine DumpDomainDecompHistogram
end module ModDomainDecomp
