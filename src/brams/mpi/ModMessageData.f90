module ModMessageData
  use ModParallelEnvironment, only: &
       ParallelEnvironment, &
       MsgDump

  use ModNeighbourNodes, only: &
       NeighbourNodes

  use ModDomainDecomp, only: &
       DomainDecomp

  use ModFieldSectionList, only: &
       FieldSection, &
       CreateFieldSection, &
       FieldSectionList, &
       CreateFieldSectionList, &
       DestroyFieldSectionList, &
       InsertAtFieldSectionList, &
       DumpFieldSectionList

  use var_tables, only: &
       var_tables_r
  

  implicit none
  private
  public :: MessageData
  public :: InitializeMessageData
  public :: TransferMessageData
  public :: InsertFieldSectionAtSendRecvMessageData
  public :: CleanMessageData
  public :: DestroyMessageData

  ! data to send/receive to/from one node in one message

  type MessageData
     real, allocatable :: buf(:)  ! message buffer
     integer :: bufSize=0         ! message buffer size
     type (FieldSectionList), pointer :: fieldList=> null() ! field sections to communicate
  end type MessageData

contains



  ! InitializeMessageData: Create MessageData components of each entry
  !                        of a previously allocated MesageData array 



  subroutine InitializeMessageData (MsgData)
    type(MessageData), intent(inout) :: MsgData(:)

    integer :: i
    character(len=*), parameter :: h="**(InitializeMessageData)**"

    do i = 1, size(MsgData)
       MsgData(i)%bufSize = 0
       MsgData(i)%fieldList => CreateFieldSectionList()
    end do
  end subroutine InitializeMessageData



  ! TransferMessageData: assignment of MessageData variables



  subroutine TransferMessageData(left, right)
    type(MessageData) :: left
    type(MessageData) :: right
    left%bufSize = right%bufSize
    left%fieldList => right%fieldList
  end subroutine TransferMessageData



  ! InsertFieldSectionAtSendRecvMessageData: 
  !    Insert, at previosly existing message data variables,
  !    field sections to communicate. 
  !    The field sections are arrays indexed by number of neighbours.
  !    Logical arrays indicate if each neighbour has field
  !    sections to communicate.
  !    Insertion is performed only for neighbours that have
  !    field sections to communicate.



  subroutine InsertFieldSectionAtSendRecvMessageData(&
       vTabPtr, ParEnv, Neigh, GlobalWithGhost, &
       xbSend, xeSend, ybSend, yeSend, willSend, SendMsgData, &
       xbRecv, xeRecv, ybRecv, yeRecv, willRecv, RecvMsgData)

    type(var_tables_r), pointer :: vTabPtr         ! intent(in)
    type(ParallelEnvironment), pointer :: ParEnv   ! intent(in)
    type(NeighbourNodes), pointer :: Neigh         ! intent(in)
    type(DomainDecomp), pointer :: GlobalWithGhost ! intent(in)

    ! all remaining arguments are dimensioned by number of neighbours 
    ! and indexed by neighbour number

    ! region to be sent to each neighbour (global indices)

    integer, intent(in) :: xbSend(:)
    integer, intent(in) :: xeSend(:)
    integer, intent(in) :: ybSend(:)
    integer, intent(in) :: yeSend(:)

    ! which neighbours will receive msgs from this node

    logical, intent(in) :: willSend(:)

    ! potential sending message data

    type(MessageData), intent(inout) :: SendMsgData(:)

    ! region to be received from each neighbour (global indices)

    integer, intent(in) :: xbRecv(:)
    integer, intent(in) :: xeRecv(:)
    integer, intent(in) :: ybRecv(:)
    integer, intent(in) :: yeRecv(:)

    ! which neighbours will send msgs to this node

    logical, intent(in) :: willRecv(:)

    ! potential receiving message data

    type(MessageData), intent(inout) :: RecvMsgData(:)

    integer :: i
    integer :: thisNode
    integer :: x0, y0
    type(FieldSection), pointer :: oneFieldSection
    character(len=*), parameter :: h="**(InsertFieldSectionAtSendRecvMessageData)**"

    ! check arguments

    if (.not. associated(vTabPtr)) then
       call fatal_error(h//" vTabPtr not associated")
    else if (.not. associated(ParEnv)) then
       call fatal_error(h//" ParEnv not associated")
    else if (.not. associated(Neigh)) then
       call fatal_error(h//" Neigh not associated")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" GlobalWithGhost not associated")
    end if

    ! offsets to convert global indices to local indices at this proc

    thisNode = ParEnv%myNum
    x0 = GlobalWithGhost%xb(thisNode) - 1
    y0 = GlobalWithGhost%yb(thisNode) - 1

    ! create list of Field Sections to send and to receive

    do i = 1, Neigh%nNeigh
       if (willSend(i)) then
          oneFieldSection =>  CreateFieldSection(&
               vTabPtr, &
               xbSend(i)-x0, xeSend(i)-x0, &
               ybSend(i)-y0, yeSend(i)-y0, &
               GlobalWithGhost)
          call InsertAtFieldSectionList(oneFieldSection, &
               SendMsgData(i)%fieldList)
          SendMsgData(i)%bufSize = SendMsgData(i)%bufSize+oneFieldSection%fieldSectionSize
       end if
       if (willRecv(i)) then
          oneFieldSection =>  CreateFieldSection(&
               vTabPtr, &
               xbRecv(i)-x0, xeRecv(i)-x0, &
               ybRecv(i)-y0, yeRecv(i)-y0, &
               GlobalWithGhost)
          call InsertAtFieldSectionList(oneFieldSection, &
               RecvMsgData(i)%fieldList)
          RecvMsgData(i)%bufSize = RecvMsgData(i)%bufSize+oneFieldSection%fieldSectionSize
       end if
    end do
  end subroutine InsertFieldSectionAtSendRecvMessageData



  ! CleanMessageData: To clean a message data array,
  !                   destroy field section lists of null entries,
  !                   create fresh areas for non-empty field section lists,
  !                   since these non-empty entries are pointed by someone else



  subroutine CleanMessageData (MsgData)
    type(MessageData), intent(inout) :: MsgData(:)

    integer :: i
    character(len=*), parameter :: h="**(InitializeMessageData)**"

    do i = 1, size(MsgData)
       if (MsgData(i)%bufSize == 0) then
          call DestroyFieldSectionList(MsgData(i)%fieldList)
       else
          MsgData(i)%bufSize = 0
          MsgData(i)%fieldList => null()
       end if
    end do
  end subroutine CleanMessageData



  ! DestroyMessageData: destroy a variable of this type



  subroutine DestroyMessageData(msgData)
    type(MessageData) :: msgData
    character(len=*), parameter :: h="**(DestroyMessageData)**"

!!$    call MsgDump(h//" starts")
    msgData%bufSize = 0
    if (allocated(msgData%buf)) then
       call MsgDump(h//" msgData%buf is allocated")
       deallocate(msgData%buf)
!!$       call MsgDump(h//" msgData%buf was deallocated")
!!$    else
!!$       call MsgDump(h//" msgData%buf is not allocated")
    end if
!!$    call MsgDump(h//" will destroy field section list:")
!!$    call DumpFieldSectionList(msgData%fieldList)
    call DestroyFieldSectionList(msgData%fieldList)
!!$    call MsgDump(h//" done ")
  end subroutine DestroyMessageData
end module ModMessageData
