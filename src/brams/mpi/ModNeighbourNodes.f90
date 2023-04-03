module ModNeighbourNodes

  ! ModNeighbourNodes: get all neighbours of this node.
  !                    A neighbour is a node which interior points with ghost zone 
  !                    intersects this node inner points (no ghost zone)
  !                    or is a node which interior points intersects
  !                    this node interior points with ghost zone, except that
  !                    this node is never a neighbour of itself.


  use ModParallelEnvironment, only: &
       ParallelEnvironment, &
       MsgDump

  use ModGridDims, only: &
       GridDims

  use ModDomainDecomp, only: &
       DomainDecomp

  implicit none
  private
  public :: NeighbourNodes
  public :: CreateNeighbourNodes
  public :: DumpNeighbourNodes
  public :: DestroyNeighbourNodes
  public :: Inter
  public :: NodesToSendRecvMessages
  public :: IncludeDomainBoundaries
  public :: GetNumberOfNeighbours

  type NeighbourNodes
     integer :: nNeigh ! how many neighbours
     integer, allocatable :: neigh(:) ! neighbours list
  end type NeighbourNodes

  logical, parameter :: dumpLocal=.false.

contains



  ! Inter: intersection of areas [xs1:xe1,ys1:ye1] and [xs2:xe2,ys2:ye2].
  !        Returns in hasInter if there is (or not) intersection and the
  !        intersection itself at [xsInter:xeInter,ysInter:yeInter].



  subroutine Inter(&
       xs1, xe1, ys1, ye1, &
       xs2, xe2, ys2, ye2, &
       xsInter, xeInter, ysInter, yeInter, &
       hasInter)
    integer, intent(in) :: xs1
    integer, intent(in) :: xe1
    integer, intent(in) :: ys1
    integer, intent(in) :: ye1
    integer, intent(in) :: xs2
    integer, intent(in) :: xe2
    integer, intent(in) :: ys2
    integer, intent(in) :: ye2
    integer, intent(out) :: xsInter
    integer, intent(out) :: xeInter
    integer, intent(out) :: ysInter
    integer, intent(out) :: yeInter
    logical, intent(out) :: hasInter

    xsInter = max(xs1,xs2)
    xeInter = min(xe1,xe2)
    ysInter = max(ys1,ys2)
    yeInter = min(ye1,ye2)
    hasInter = &
         xsInter <= xeInter .and. &
         ysInter <= yeInter
  end subroutine Inter



  ! CreateNeighbourNodes: Finds all neighbours of this node
  !                       Returns a null pointer if no neighbours.



  subroutine CreateNeighbourNodes(ParEnv, &
       GlobalNoGhost, GlobalWithGhost, &
       OneNeighbourNodes)
    type(ParallelEnvironment), pointer :: ParEnv
    type(DomainDecomp), pointer :: GlobalNoGhost
    type(DomainDecomp), pointer :: GlobalWithGhost
    type(NeighbourNodes), pointer :: OneNeighbourNodes

    integer :: myNum
    integer :: nmachs
    integer :: node
    integer :: nNeigh
    integer :: cnt
    integer :: xsInter, xeInter, ysInter, yeInter
    logical :: myNumSend, myNumRecv
    logical, allocatable :: isNeighbour(:)
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(CreateNeighbourNodes)**"

    if (.not. associated(ParEnv)) then
       call fatal_error(h//" invoked with null ParEnv")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" invoked with null GlobalWithGhost ")
    else if (.not. associated(GlobalNoGhost)) then
       call fatal_error(h//" invoked with null GlobalNoGhost ")
    else if (associated(oneNeighbourNodes)) then
       call fatal_error(h//" starts with already associated oneNeighbourNodes")
    end if

    nmachs = size(GlobalNoGhost%xb)
    myNum = ParEnv%myNum

    ! scratch area storing which nodes are neighbour of this node:
    ! these are nodes with extended domain (with ghost zone) intersect
    ! with this node original domain, or 
    ! nodes where the original domain intersects with this node
    ! extended domain (with ghost zone)

    allocate(isNeighbour(nmachs))
    do node = 1, nmachs
       call Inter(&
            GlobalNoGhost%xb(myNum), GlobalNoGhost%xe(myNum), &
            GlobalNoGhost%yb(myNum), GlobalNoGhost%ye(myNum), &
            GlobalWithGhost%xb(node), GlobalWithGhost%xe(node), & 
            GlobalWithGhost%yb(node), GlobalWithGhost%ye(node), &
            xsInter, xeInter, ysInter, yeInter, myNumSend)
       call Inter(&
            GlobalNoGhost%xb(node), GlobalNoGhost%xe(node), &
            GlobalNoGhost%yb(node), GlobalNoGhost%ye(node), &
            GlobalWithGhost%xb(myNum), GlobalWithGhost%xe(myNum), & 
            GlobalWithGhost%yb(myNum), GlobalWithGhost%ye(myNum), &
            xsInter, xeInter, ysInter, yeInter, myNumRecv)
       isNeighbour(node) = myNumSend .or. myNumRecv
    end do

    ! exclude own node

    isNeighbour(myNum) = .false.

    ! how many neighbours; default is no neighbour

    nNeigh = count(isNeighbour)
    OneNeighbourNodes => null()
    if (nNeigh /= 0) then

       ! finishes building OneNeighbourNodes

       allocate(OneNeighbourNodes)
       OneNeighbourNodes%nNeigh = nNeigh
       allocate(OneNeighbourNodes%neigh(OneNeighbourNodes%nNeigh))

       cnt = 0
       do node =  1, nmachs
          if (isNeighbour(node)) then
             cnt = cnt + 1
             OneNeighbourNodes%neigh(cnt) = node
          end if
       end do
       if (cnt /= OneNeighbourNodes%nNeigh) then
          write(c0,"(i8)") cnt
          write(c1,"(i8)") OneNeighbourNodes%nNeigh
          call fatal_error(h//" inconsistence: cnt ("//trim(adjustl(c0))//&
               ") differs from nNeigh ("//trim(adjustl(c1))//")")
       end if
    end if
    deallocate(isNeighbour)

    if (dumpLocal) then
       call DumpNeighbourNodes(OneNeighbourNodes)
    end if
  end subroutine CreateNeighbourNodes



  ! DumpNeighbourNodes: at this process dump file



  subroutine DumpNeighbourNodes(OneNeighbourNodes)
    type(NeighbourNodes), pointer :: OneNeighbourNodes
    integer :: neigh
    character(len=8) :: c0
    character(len=*), parameter :: h="**(DumpNeighbourNodes)**"

    if (.not. associated(OneNeighbourNodes)) then
       call MsgDump(h//" empty")
    else
       write(c0,"(i8)") OneNeighbourNodes%nNeigh
       call MsgDump(h//" there are "//trim(adjustl(c0))//&
            " neighbours: nodes", .true.)
       do neigh = 1, OneNeighbourNodes%nNeigh-1
          write(c0,"(i8)") OneNeighbourNodes%neigh(neigh)
          call MsgDump(" "//trim(adjustl(c0))//",",.true.)
       end do
       write(c0,"(i8)") OneNeighbourNodes%neigh(neigh)
       call MsgDump(" "//trim(adjustl(c0)))
    end if
  end subroutine DumpNeighbourNodes



  ! DestroyNeighbourNodes: returns allocated area, if any



  subroutine DestroyNeighbourNodes(OneNeighbourNodes)
    type(NeighbourNodes), pointer :: OneNeighbourNodes
    if (associated(OneNeighbourNodes)) then
       deallocate(OneNeighbourNodes%neigh)
!--(inspxe)-----------------------------------------------
! Correcao de bug: Memory leak
       deallocate(OneNeighbourNodes)
!--(inspxe)-----------------------------------------------
    end if
    nullify(OneNeighbourNodes)
  end subroutine DestroyNeighbourNodes



  ! NodesToSendRecvMessages: given a region to be updated at each process,
  !                          find which processes this node will send messages
  !                          to and receive messages from to update the area.
  !                          Also returns the region of each process to be
  !                          sent or received.



  subroutine NodesToSendRecvMessages(thisNode, Neigh, GlobalNoGhost, &
       xbToBeUpdated, xeToBeUpdated, ybToBeUpdated, yeToBeUpdated, &
       xbSend, xeSend, ybSend, yeSend, willSend, &
       xbRecv, xeRecv, ybRecv, yeRecv, willRecv)

    integer, intent(in) :: thisNode              ! BRAMS process number
    type(NeighbourNodes), pointer :: Neigh       ! intent(in)
    type(DomainDecomp), pointer :: GlobalNoGhost ! intent(in)

    ! region to be updated at each node (included at GlobalWithGhost)

    integer, intent(in) :: xbToBeUpdated(:)      ! global index, size of ParEnv%nmachs
    integer, intent(in) :: xeToBeUpdated(:)      ! global index, size of ParEnv%nmachs
    integer, intent(in) :: ybToBeUpdated(:)      ! global index, size of ParEnv%nmachs
    integer, intent(in) :: yeToBeUpdated(:)      ! global index, size of ParEnv%nmachs

    ! region to be sent to each neighbour

    integer, intent(out) :: xbSend(:)            ! global index, size of Neigh%nNeigh
    integer, intent(out) :: xeSend(:)            ! global index, size of Neigh%nNeigh
    integer, intent(out) :: ybSend(:)            ! global index, size of Neigh%nNeigh
    integer, intent(out) :: yeSend(:)            ! global index, size of Neigh%nNeigh

    ! will send msgs to which neighbours

    logical, intent(out) :: willSend(:)          ! global index, size of Neigh%nNeigh

    ! region to be received from each neighbour

    integer, intent(out) :: xbRecv(:)            ! global index, size of Neigh%nNeigh
    integer, intent(out) :: xeRecv(:)            ! global index, size of Neigh%nNeigh
    integer, intent(out) :: ybRecv(:)            ! global index, size of Neigh%nNeigh
    integer, intent(out) :: yeRecv(:)            ! global index, size of Neigh%nNeigh

    ! will receive msgs from which neighbours

    logical, intent(out) :: willRecv(:)          ! global index, size of Neigh%nNeigh

    integer :: otherNode
    integer :: nNeigh
    integer :: i, indMsg
    character(len=8) :: c0, c1, c2, c3, c4
    character(len=*), parameter :: h="**(NodesToSendRecvMessages)**"

    ! check arguments

    if (.not. associated(GlobalNoGhost)) then
       call fatal_error(h//" starts with null GlobalNoGhost")
    end if

    ! default output

    xbSend=0
    xeSend=0
    ybSend=0
    yeSend=0
    willSend=.false.

    xbRecv=0
    xeRecv=0
    ybRecv=0
    yeRecv=0
    willRecv=.false.

    ! no messages if no neighbours

    if (.not. associated(Neigh)) then
       if (dumpLocal) then
          call MsgDump(h//" no neighbour")
       end if
       return
    end if

    ! auxiliar variables

    nNeigh = Neigh%nNeigh

    ! this node will send messages to other node
    ! if the region of the other node to be updated
    ! intersects with the region that this node owns.
    ! find which neighbours this node will send messages to
    ! and the regions (global indices) to be received.

    do i = 1, nNeigh
       otherNode = Neigh%neigh(i)
       call Inter(&
            xbToBeUpdated(otherNode), &
            xeToBeUpdated(otherNode), &
            ybToBeUpdated(otherNode), &
            yeToBeUpdated(otherNode), &
            GlobalNoGhost%xb(thisNode), &
            GlobalNoGhost%xe(thisNode), &
            GlobalNoGhost%yb(thisNode), &
            GlobalNoGhost%ye(thisNode), &
            xbSend(i), xeSend(i), &
            ybSend(i), yeSend(i), &
            willSend(i))

    !!!!!! Modificar        
       if (willSend(i) .and. dumpLocal) then
          write(c0,"(i8)") otherNode
          write(c1,"(i8)") xbSend(i)
          write(c2,"(i8)") xeSend(i)
          write(c3,"(i8)") ybSend(i)
          write(c4,"(i8)") yeSend(i)
          call MsgDump(h//&
               " will send to BRAMS process "//trim(adjustl(c0))//&
               " the region ["//&
               trim(adjustl(c1))//":"//trim(adjustl(c2))//","//&
               trim(adjustl(c3))//":"//trim(adjustl(c4))//"]")
       end if
    end do

    ! this node will receive messages from other nodes
    ! if the region of this node to be updated
    ! intersects with the region owned by the other node.
    ! find which neighbours will receive messages from this node
    ! and the regions (global indices) to be received.

    do i = 1, nNeigh
       otherNode = Neigh%neigh(i)
       call Inter(&
            xbToBeUpdated(thisNode), &
            xeToBeUpdated(thisNode), &
            ybToBeUpdated(thisNode), &
            yeToBeUpdated(thisNode), &
            GlobalNoGhost%xb(otherNode), &
            GlobalNoGhost%xe(otherNode), &
            GlobalNoGhost%yb(otherNode), &
            GlobalNoGhost%ye(otherNode), &
            xbRecv(i), xeRecv(i), &
            ybRecv(i), yeRecv(i), &
            willRecv(i))

            
       if (willRecv(i) .and. dumpLocal) then
          write(c0,"(i8)") otherNode
          write(c1,"(i8)") xbRecv(i)
          write(c2,"(i8)") xeRecv(i)
          write(c3,"(i8)") ybRecv(i)
          write(c4,"(i8)") yeRecv(i)
          call MsgDump(h//" will recv from RAMS node "//trim(adjustl(c0))//" ["//&
               trim(adjustl(c1))//":"//trim(adjustl(c2))//","//&
               trim(adjustl(c3))//":"//trim(adjustl(c4))//"]")
       end if
    end do
  end subroutine NodesToSendRecvMessages






  subroutine IncludeDomainBoundaries(Neigh, GridSize, GlobalNoGhost, &
       xbComm, xeComm, ybComm, yeComm, willComm)

    type(NeighbourNodes), pointer :: Neigh       ! intent(in)
    type(GridDims), pointer :: GridSize       ! intent(in)
    type(DomainDecomp), pointer :: GlobalNoGhost ! intent(in)

    integer, intent(inout) :: xbComm(:)          ! global index, size of Neigh%nNeigh
    integer, intent(inout) :: xeComm(:)          ! global index, size of Neigh%nNeigh
    integer, intent(inout) :: ybComm(:)          ! global index, size of Neigh%nNeigh
    integer, intent(inout) :: yeComm(:)          ! global index, size of Neigh%nNeigh
    logical, intent(inout) :: willComm(:)        ! global index, size of Neigh%nNeigh

    integer :: indNeigh
    integer :: proc
    character(len=8) :: c0, c1, c2, c3, c4
    character(len=*), parameter :: h="**(IncludeDomainBoundaries)**"

    ! check arguments

    if (.not. associated(GlobalNoGhost)) then
       call fatal_error(h//" starts with null GlobalNoGhost")
    end if

    if (.not. associated(Neigh)) then
       return
    end if

    do indNeigh = 1, Neigh%nNeigh
       if (willComm(indNeigh)) then
          proc = Neigh%neigh(indNeigh)
          if (xbComm(indNeigh) == 2) then
             xbComm(indNeigh) = 1
             if (dumpLocal) then
                write(c0,"(i8)") proc
                write(c1,"(i8)") indNeigh
                write(c2,"(i8)") xbComm(indNeigh)
                call MsgDump(h//" xbComm("//trim(adjustl(c1))//&
                     ") = "//trim(adjustl(c2))//&
                     " (proc="//trim(adjustl(c0))//")")
             end if
          end if
          if (xeComm(indNeigh) == GridSize%nnxp-1) then
             xeComm(indNeigh) = GridSize%nnxp
             if (dumpLocal) then
                write(c0,"(i8)") proc
                write(c1,"(i8)") indNeigh
                write(c2,"(i8)") xeComm(indNeigh)
                call MsgDump(h//" xeComm("//trim(adjustl(c1))//&
                     ") = "//trim(adjustl(c2))//&
                     " (proc="//trim(adjustl(c0))//")")
             end if
          end if
          if (ybComm(indNeigh) == 2) then
             ybComm(indNeigh) = 1
             if (dumpLocal) then
                write(c0,"(i8)") proc
                write(c1,"(i8)") indNeigh
                write(c2,"(i8)") ybComm(indNeigh)
                call MsgDump(h//" ybComm("//trim(adjustl(c1))//&
                     ") = "//trim(adjustl(c2))//&
                     " (proc="//trim(adjustl(c0))//")")
             end if
          end if
          if (yeComm(indNeigh) == GridSize%nnyp-1) then
             yeComm(indNeigh) = GridSize%nnyp
             if (dumpLocal) then
                write(c0,"(i8)") proc
                write(c1,"(i8)") indNeigh
                write(c2,"(i8)") yeComm(indNeigh)
                call MsgDump(h//" yeComm("//trim(adjustl(c1))//&
                     ") = "//trim(adjustl(c2))//&
                     " (proc="//trim(adjustl(c0))//")")
             end if
          end if
       end if
    end do
  end subroutine IncludeDomainBoundaries






  integer function GetNumberOfNeighbours(OneNeighbourNodes)
    type(NeighbourNodes), pointer :: OneNeighbourNodes
    integer :: neigh
    character(len=8) :: c0
    character(len=*), parameter :: h="**(GetNumberOfNeighbours)**"

    if (.not. associated(OneNeighbourNodes)) then
       call MsgDump(h//" empty")
    else
       GetNumberOfNeighbours = OneNeighbourNodes%nNeigh
    end if
  end function GetNumberOfNeighbours
end module ModNeighbourNodes

