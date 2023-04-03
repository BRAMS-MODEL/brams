module ModMessageSet
  use ModParallelEnvironment, only: &
       ParallelEnvironment, &
       Brams2MpiProcNbr, &
       Mpi2BramsProcNbr, &
       MsgDump

  use ModMessageData, only: &
       MessageData, &
       TransferMessageData, &
       DestroyMessageData

  use ModNeighbourNodes, only: &
       NeighbourNodes

  use ModFieldSectionList, only: &
       FieldSection, &
       CreateFieldSection, &
       NextFieldSection, &
       DumpFieldSection, &
       CreateFieldSectionList, &
       InsertAtFieldSectionList, &
       DumpFieldSectionList

  use ParLib, only: &
       parf_get_noblock_real, &
       parf_send_noblock_real, &
       parf_wait_any_nostatus, &
       parf_wait_all_nostatus

  use ModBuffering, only: &
       FieldSection2Buffer, &
       Buffer2FieldSection

  use ModDomainDecomp, only: &
       DomainDecomp

  use var_tables, only: &
       var_tables_r

!  use CUPARM_GRELL3, only: g3d_g

  implicit none
  private
  public :: MessageSet
  public :: CreateMessageSet
  public :: InsertFieldSectionAtMessageSet
  public :: DumpMessageSet
  public :: DestroyMessageSet
  public :: PostRecvSendMsgs
  public :: WaitRecvMsgs

  include "mpif.h"

  integer, parameter :: UNDEFINED=-1

  ! all messages to send/receive from one process
  ! to update a field.
  ! arrays are indexed 1:nMsgs

  type MessageSet
     character(len=64) :: name                   ! msg name
     integer :: nMsgs=UNDEFINED                  ! # procs on communication
     integer :: tag                              ! same tag for all messages in the set
     type(MessageData), allocatable :: oneMsg(:) ! data on communication
     integer, allocatable :: request(:)          ! communication request
     integer, allocatable :: otherProc(:)        ! MPI procs to communicate
  end type MessageSet

  logical, parameter :: dumpLocal=.false.

contains


  ! CreateMessageSet: Generates a variable of type MessageSet containing
  !                   all message envelopes and no message data
  !                   to be sent by this node to neighbour nodes or
  !                   to be received by this node from neighbour nodes
  !                   during the communication denoted by "name".
  !                   Input includes "hasMsg", a logical array indexed by
  !                   neighbour number that stores if a neighbour will or
  !                   will not communicate with this node in this
  !                   communication.
  !                   Same tag should be used on send and recieve.
  !                   Returning variable has as many messages as the number
  !                   of true values in hasMsg. Arrays in returning variable
  !                   are indexed by true value count, at range 1:nMsgs.



  function CreateMessageSet (name, tag, hasMsg, Neigh) result(Msgs)
    character(len=*), intent(in) :: name
    integer, intent(in) :: tag
    logical, intent(in) :: hasMsg(:)
    type(NeighbourNodes), pointer :: Neigh
    type(MessageSet), pointer :: Msgs

    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(CreateMessageSet)**"

    integer :: nNeigh
    integer :: nMsgs
    integer :: iNeigh
    integer :: cntMsg

    ! no message if no neighbours or no node has messages

    nNeigh = Neigh%nNeigh
    nMsgs = count(hasMsg)
    if (nNeigh == 0 .or. nMsgs == 0) then
       Msgs => null()
       return
    end if

    ! there are messages: allocate area

    allocate(Msgs)
    allocate(Msgs%oneMsg(nMsgs))
    allocate(Msgs%request(nMsgs))
    allocate(Msgs%otherProc(nMsgs))

    ! store basic info

    Msgs%name = name
    Msgs%nMsgs = nMsgs
    Msgs%tag = tag

    ! for each neighbour node that will communicate,
    ! build message envelop to/from neighbour and
    ! create empty message data.

    cntMsg=0
    do iNeigh = 1, nNeigh
       if (hasMsg(iNeigh)) then
          cntMsg = cntMsg + 1
          Msgs%oneMsg(cntMsg)%bufSize=0
          Msgs%oneMsg(cntMsg)%fieldList => CreateFieldSectionList()
          Msgs%request(cntMsg) = MPI_REQUEST_NULL
          Msgs%otherProc(cntMsg)= Brams2MpiProcNbr(Neigh%neigh(iNeigh))
       end if
    end do

    if (cntMsg /= nMsgs) then
       write(c0,"(i8)") cntMsg
       write(c1,"(i8)") nMsgs
       call fatal_error(h//" inconsistency: cntMsg="//trim(adjustl(c0))//&
            " while nMsgs="//trim(adjustl(c1)))
    end if
  end function CreateMessageSet



  ! InsertFieldSectionAtMessageSet:
  ! CreateMessageSet: Generates a variable of type MessageSet containing
  !                   all message envelopes and no message data
  !                   to be sent by this node to neighbour nodes or
  !                   to be received by this node from neighbour nodes
  !                   during the communication denoted by "name".
  !                   Input includes "hasMsg", a logical array indexed by
  !                   neighbour number that stores if a neighbour will or
  !                   will not communicate with this node in this
  !                   communication.
  !                   Same tag should be used on send and recieve.
  !                   Returning variable has as many messages as the number
  !                   of true values in hasMsg. Arrays in returning variable
  !                   are indexed by true value count, at range 1:nMsgs.



  subroutine InsertFieldSectionAtMessageSet(&
       myNum, vTabPtr, Neigh, GlobalWithGhost, &
       xbComm, xeComm, ybComm, yeComm, willComm, Msgs)

    integer, intent(in) :: myNum
    type(var_tables_r), pointer :: vTabPtr         ! intent(in)
    type(ParallelEnvironment), pointer :: ParEnv   ! intent(in)
    type(NeighbourNodes), pointer :: Neigh         ! intent(in)
    type(DomainDecomp), pointer :: GlobalWithGhost ! intent(in)

    ! all remaining arguments are dimensioned by number of neighbours
    ! and indexed by neighbour number

    ! region to be sent to each neighbour (global indices)

    integer, intent(in) :: xbComm(:)
    integer, intent(in) :: xeComm(:)
    integer, intent(in) :: ybComm(:)
    integer, intent(in) :: yeComm(:)

    ! which neighbours will receive msgs from this node

    logical, intent(in) :: willComm(:)

    ! potential sending message data

    type(MessageSet), pointer :: Msgs

    integer :: nMsgs
    integer :: x0, y0
    integer :: cntMsg
    integer :: iNeigh
    type(FieldSection), pointer :: oneFieldSection
    character(len=8) :: c0
    character(len=*), parameter :: h="**(InsertFieldSectionAtMessageSet)**"

    ! check arguments

    if (.not. associated(vTabPtr)) then
       call fatal_error(h//" vTabPtr not associated")
    else if (.not. associated(Neigh)) then
       call fatal_error(h//" Neigh not associated")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" GlobalWithGhost not associated")
    end if

    ! return if no messages to send

    if (.not. associated(Msgs)) then
       !print *,'Msgs not associated, return!'; call flush(6)
       return
    end if
    nMsgs = Msgs%nMsgs

    ! offsets to convert global indices to local indices at this proc

    x0 = GlobalWithGhost%xb(myNum) - 1
    y0 = GlobalWithGhost%yb(myNum) - 1

    ! create list of Field Sections to communicate

    cntMsg = 0
    do iNeigh = 1, Neigh%nNeigh
       if (willComm(iNeigh)) then
          cntMsg = cntMsg + 1
          if (cntMsg > nMsgs) then
             write(c0,"(i8)") nMsgs
             call fatal_error(h//" nMsgs ("//&
                  trim(adjustl(c0))//") exceeded while inserting field "//&
                  trim(adjustl(vTabPtr%name))//&
                  " at message "//trim(adjustl(Msgs%name)))
          end if
          oneFieldSection =>  CreateFieldSection(&
               vTabPtr, &
               xbComm(iNeigh)-x0, xeComm(iNeigh)-x0, &
               ybComm(iNeigh)-y0, yeComm(iNeigh)-y0, &
               GlobalWithGhost)
          call InsertAtFieldSectionList(&
               oneFieldSection, &
               Msgs%oneMsg(cntMsg)%fieldList)
          Msgs%oneMsg(cntMsg)%bufSize = Msgs%oneMsg(cntMsg)%bufSize + &
               oneFieldSection%fieldSectionSize
       end if
    end do
  end subroutine InsertFieldSectionAtMessageSet




  subroutine DumpMessageSet(Msgs)
    type(MessageSet), pointer :: Msgs

    character(len=8) :: c0, c1, c2, c3, c4
    character(len=*), parameter :: h="**(DumpMessageSet)**"
    integer :: i

    if (.not. associated(Msgs)) then
       call MsgDump(h//" empty messages")
    else
       write(c0,"(i8)") Msgs%nMsgs
       if (Msgs%tag == UNDEFINED) then
          c2="UNDEF"
       else
          write(c2,"(i8)") Msgs%tag
       end if
       call MsgDump(h//" named "//trim(adjustl(Msgs%name))//&
            " with "//trim(adjustl(c0))//" messages and tag "//&
            trim(adjustl(c2)))
       do i = 1, Msgs%nMsgs
          write(c0,"(i8)") i
          if (Msgs%otherProc(i) == UNDEFINED) then
             c1="UNDEF"
          else
             write(c1,"(i8)") Mpi2BramsProcNbr(Msgs%otherProc(i))
          end if
          if (Msgs%oneMsg(i)%bufSize == UNDEFINED) then
             c3="UNDEF"
          else
             write(c3,"(i8)") Msgs%oneMsg(i)%bufSize
          end if
          if (Msgs%request(i) == MPI_REQUEST_NULL) then
             c4="NULL"
          else
             write(c4,"(Z8)") Msgs%request(i)
          end if
          call MsgDump(h//" message "//trim(adjustl(c0))//&
               " to/from node "//trim(adjustl(c1))//&
               ", request "//trim(adjustl(c4))//&
               ", size "//trim(adjustl(c3))//&
               " and field sections:")
          call DumpFieldSectionList(Msgs%oneMsg(i)%fieldList)
       end do
    end if
  end subroutine DumpMessageSet





  subroutine DestroyMessageSet (Msgs)
    type(MessageSet), pointer :: Msgs

    integer :: msg
    character(len=8) :: c0
    character(len=*), parameter :: h="**(DestroyMessageSet)**"

    if (associated(Msgs)) then
       do msg = 1, Msgs%nMsgs
          call DestroyMessageData(Msgs%oneMsg(msg))
       end do
       deallocate(Msgs%oneMsg)
       deallocate(Msgs%request)
       deallocate(Msgs%otherProc)
       deallocate(Msgs)
    end if
    Msgs => null()
  end subroutine DestroyMessageSet




  subroutine PostRecvSendMsgs(SendMsg, RecvMsg)
    type(MessageSet), pointer :: SendMsg
    type(MessageSet), pointer :: RecvMsg

    integer :: iSend
    integer :: iRecv
    integer :: firstBuffer
    integer :: lastBuffer
    type(MessageData), pointer :: msgData => null()
    type(FieldSection), pointer :: node => null()
    character(len=8) :: c0, c1, c2, c3, c4, c5
    character(len=*), parameter :: h="**(PostRecvSendMsgs)**"

    if (dumpLocal) then
       call MsgDump(h//" enter with SendMsg:")
       call DumpMessageSet(SendMsg)
       call MsgDump(h//" enter with RecvMsg:")
       call DumpMessageSet(RecvMsg)
       call MsgDump(h//" finished dumping input arguments")
    end if

    ! post non-blocking receive for each receiving message

    if (associated(RecvMsg)) then
       do iRecv= 1,RecvMsg%nMsgs
          msgData => RecvMsg%oneMsg(iRecv)

          ! allocate receive buffer

          allocate(msgData%buf(msgData%bufSize))

          ! post receive

          call parf_get_noblock_real(msgData%buf, msgData%bufSize, &
               RecvMsg%otherProc(iRecv), &
               RecvMsg%tag, RecvMsg%request(iRecv))

          if (dumpLocal) then
             write(c1,"(i8)") RecvMsg%otherProc(iRecv)
             write(c2,"(i8)") size(msgData%buf)
             write(c3,"(i8)") RecvMsg%tag
             if (RecvMsg%request(iRecv) == MPI_REQUEST_NULL) then
                c4="NULL"
             else
                write(c4,"(Z8)") RecvMsg%request(iRecv)
             end if
             call MsgDump(h//" for "//trim(adjustl(RecvMsg%name))//&
                  " post recv from MPI node "//trim(adjustl(c1))//&
                  " with buffer size "//trim(adjustl(c2))//&
                  " tag "//trim(adjustl(c3))//" and request "//trim(adjustl(c4)))
          end if
       end do
    else
       if (dumpLocal) then
          call MsgDump(h//" no message to receive")
       end if
    end if

    ! for each sending message,
    ! build sending buffer copying field sections to the buffer and
    ! send non-blocking message

    if (associated(SendMsg)) then
       do iSend = 1,SendMsg%nMsgs
          if (dumpLocal) then
             write(c0,"(i8)") iSend
             call MsgDump(h//" sending message "//trim(adjustl(c0)))
          end if
          msgData => SendMsg%oneMsg(iSend)

          allocate(msgData%buf(msgData%bufSize))

          node => null()
          lastBuffer=0
          if (dumpLocal) then
             call MsgDump(h//" starts building sending buffer")
          end if
          do
             node => NextFieldSection(node, msgData%fieldList)
             if (dumpLocal) then
                call MsgDump(h//" next FieldSection to insert at sending buffer:")
                call DumpFieldSection(node)
             end if
             if (associated(node)) then
                if (associated(node%vTabPtr%var_p_2D)) then
                   firstBuffer = lastBuffer+1
                   call FieldSection2Buffer(node%vTabPtr%var_p_2D, &
                        node%vTabPtr%idim_type, &
                        node%xStart, node%xEnd,&
                        node%yStart, node%yEnd,&
                        msgData%buf, lastBuffer)

                   if (dumpLocal) then
                      write(c0,"(i8)") firstBuffer
                      write(c1,"(i8)") lastBuffer
                      write(c2,"(i8)") node%xStart
                      write(c3,"(i8)") node%xEnd
                      write(c4,"(i8)") node%yStart
                      write(c5,"(i8)") node%yEnd
                      call MsgDump(h//" filled buf["//trim(adjustl(c0))//&
                           ":"//trim(adjustl(c1))//&
                           "] with 2D field "//trim(adjustl(node%vTabPtr%name))//"["//&
                           trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
                           trim(adjustl(c4))//":"//trim(adjustl(c5))//"]")
                   end if

                else if (associated(node%vTabPtr%var_p_3D)) then
                   firstBuffer = lastBuffer+1
                   call FieldSection2Buffer(node%vTabPtr%var_p_3D, &
                        node%vTabPtr%idim_type, &
                        node%xStart, node%xEnd,&
                        node%yStart, node%yEnd,&
                        msgData%buf, lastBuffer)

                   if (dumpLocal) then
                      write(c0,"(i8)") firstBuffer
                      write(c1,"(i8)") lastBuffer
                      write(c2,"(i8)") node%xStart
                      write(c3,"(i8)") node%xEnd
                      write(c4,"(i8)") node%yStart
                      write(c5,"(i8)") node%yEnd
                      select case (node%vTabPtr%idim_type)
                      case(3)
                         call MsgDump(h//" filled buf["//trim(adjustl(c0))//&
                              ":"//trim(adjustl(c1))//&
                              "] with 3D field "//trim(adjustl(node%vTabPtr%name))//"[:,"//&
                              trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
                              trim(adjustl(c4))//":"//trim(adjustl(c5))//"]")
                      case(6:7)
                         call MsgDump(h//" filled buf["//trim(adjustl(c0))//&
                              ":"//trim(adjustl(c1))//&
                              "] with 3D field "//trim(adjustl(node%vTabPtr%name))//"["//&
                              trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
                              trim(adjustl(c4))//":"//trim(adjustl(c5))//",:]")
                      end select
                   end if

                else if (associated(node%vTabPtr%var_p_4D)) then
                   firstBuffer = lastBuffer+1
                   call FieldSection2Buffer(node%vTabPtr%var_p_4D, &
                        node%vTabPtr%idim_type, &
                        node%xStart, node%xEnd,&
                        node%yStart, node%yEnd,&
                        msgData%buf, lastBuffer)

                   if (dumpLocal) then
                      write(c0,"(i8)") firstBuffer
                      write(c1,"(i8)") lastBuffer
                      write(c2,"(i8)") node%xStart
                      write(c3,"(i8)") node%xEnd
                      write(c4,"(i8)") node%yStart
                      write(c5,"(i8)") node%yEnd
                      call MsgDump(h//" filled buf["//trim(adjustl(c0))//&
                           ":"//trim(adjustl(c1))//&
                           "] with 4D field "//trim(adjustl(node%vTabPtr%name))//"[:,"//&
                           trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
                           trim(adjustl(c4))//":"//trim(adjustl(c5))//",:]")
                   end if
                else
                   call fatal_error(h//" inconsistent var_tables entry named "//&
                        trim(adjustl(node%vTabPtr%name)))
                end if
             else
                if (dumpLocal) then
                   call MsgDump(h//" finishes building sending buffer")
                end if
                exit
             end if
          end do
          if (lastBuffer /= size(msgData%buf)) then
             write(c0,"(i8)") lastBuffer
             write(c1,"(i8)") size(msgData%buf)
             call fatal_error(h//" the send buffer of size "//trim(adjustl(c1))//&
                  " was filled with "//trim(adjustl(c0))//" entries")
          end if
          call parf_send_noblock_real(msgData%buf, msgData%bufSize, &
               SendMsg%otherProc(iSend), &
               SendMsg%tag, SendMsg%request(iSend))
          if (dumpLocal) then
             write(c1,"(i8)") SendMsg%otherProc(iSend)
             write(c2,"(i8)") size(msgData%buf)
             write(c3,"(i8)") SendMsg%tag
             if (SendMsg%request(iSend) == MPI_REQUEST_NULL) then
                c4="NULL"
             else
                write(c4,"(Z8)") SendMsg%request(iSend)
             end if
             call MsgDump(h//" sends to MPI node "//trim(adjustl(c1))//&
                  " buffer of size "//trim(adjustl(c2))//&
                  " tag "//trim(adjustl(c3))//" and request "//trim(adjustl(c4)))
          end if
       end do
    else
       if (dumpLocal) then
          call MsgDump(h//" no message to send")
       end if
    end if
  end subroutine PostRecvSendMsgs




  subroutine WaitRecvMsgs(SendMsg, RecvMsg)
    type(MessageSet), pointer :: SendMsg
    type(MessageSet), pointer :: RecvMsg

    integer :: i
    integer :: iSend
    integer :: iRecv
    integer :: firstBuffer
    integer :: lastBuffer
    integer :: recvNbr
    integer :: sendNbr
    type(MessageData), pointer :: msgData => null()
    type(FieldSection), pointer :: node => null()
    character(len=8) :: c0, c1, c2, c3, c4, c5
    character(len=*), parameter :: h="**(WaitRecvMsgs)**"

    ! for each receive message:
!print *,'LFR-DBG: inside WMS 1: ',size(g3d_g(1)%cugd_ttens,1), &
!                 size(g3d_g(1)%cugd_ttens,2),size(g3d_g(1)%cugd_ttens,3); call flush(6)
    if (associated(RecvMsg)) then
       if (dumpLocal) then
          write(c0,"(i8)") RecvMsg%nMsgs
          call MsgDump(h//" for "//trim(adjustl(RecvMsg%name))//&
               " waits on "//trim(adjustl(c0))//" recvs with requests: ",.true.)
          do iRecv= 1,RecvMsg%nMsgs-1
             if (RecvMsg%request(iRecv) == MPI_REQUEST_NULL) then
                c0="NULL"
             else
                write(c0,"(Z8)") RecvMsg%request(iRecv)
             end if
             call MsgDump(" "//trim(adjustl(c0))//",",.true.)
          end do
          if (RecvMsg%request(RecvMsg%nMsgs) == MPI_REQUEST_NULL) then
             c0="NULL"
          else
             write(c0,"(Z8)") RecvMsg%request(RecvMsg%nMsgs)
          end if
          call MsgDump(" "//trim(adjustl(c0)))
       end if
!print *,'LFR-DBG: inside WMS 2: ',size(g3d_g(1)%cugd_ttens,1), &
!                 size(g3d_g(1)%cugd_ttens,2),size(g3d_g(1)%cugd_ttens,3); call flush(6)
       !

       do iRecv= 1,RecvMsg%nMsgs

          ! wait on any arrived message

          if (dumpLocal) then
             call MsgDump(h//" Dump requests before wait: ",.true.)
             do i= 1,RecvMsg%nMsgs-1
                if (RecvMsg%request(i) == MPI_REQUEST_NULL) then
                   c0="NULL"
                else
                   write(c0,"(Z8)") RecvMsg%request(i)
                end if
                call MsgDump(" "//trim(adjustl(c0))//",",.true.)
             end do
             if (RecvMsg%request(RecvMsg%nMsgs) == MPI_REQUEST_NULL) then
                c0="NULL"
             else
                write(c0,"(Z8)") RecvMsg%request(RecvMsg%nMsgs)
             end if
             call MsgDump(" "//trim(adjustl(c0)))
          end if
          call parf_wait_any_nostatus(RecvMsg%nMsgs, &
               RecvMsg%request, recvNbr)
          if (dumpLocal) then
             write(c0,"(i8)") recvNbr
             call MsgDump(h//" received message #"//trim(adjustl(c0)))
             call MsgDump(h//" Dump requests after wait: ",.true.)
             do i= 1,RecvMsg%nMsgs-1
                if (RecvMsg%request(i) == MPI_REQUEST_NULL) then
                   c0="NULL"
                else
                   write(c0,"(Z8)") RecvMsg%request(i)
                end if
                call MsgDump(" "//trim(adjustl(c0))//",",.true.)
             end do
             if (RecvMsg%request(RecvMsg%nMsgs) == MPI_REQUEST_NULL) then
                c0="NULL"
             else
                write(c0,"(Z8)") RecvMsg%request(RecvMsg%nMsgs)
             end if
             call MsgDump(" "//trim(adjustl(c0)))
          end if
          msgData => RecvMsg%oneMsg(recvNbr)
!print *,'LFR-DBG: inside WMS 3: ',size(g3d_g(1)%cugd_ttens,1), &
!                 size(g3d_g(1)%cugd_ttens,2),size(g3d_g(1)%cugd_ttens,3); call flush(6)
          ! extract field sections from incoming buffer

          node => null()
          lastBuffer=0
          do
!            print *,'LFR-DBG: inside WMS 3.1: ',size(g3d_g(1)%cugd_ttens,1), &
!                 size(g3d_g(1)%cugd_ttens,2),size(g3d_g(1)%cugd_ttens,3); call flush(6)
             node => NextFieldSection(node, msgData%fieldList)
             if (associated(node)) then
                if (associated(node%vTabPtr%var_p_2D)) then
                   firstBuffer = lastBuffer+1


!            print *,'LFR-DBG: inside WMS 3.1.0: ', &
!                        node%vTabPtr%idim_type, &
!                        node%xStart, node%xEnd,&
!                        node%yStart, node%yEnd,node%vTabPtr%name; call flush(6)


!                   print *,'LFR-DBG: 2D: ',trim(adjustl(node%vTabPtr%name)),node%vTabPtr%idim_type; call flush(6)
                   call Buffer2FieldSection(node%vTabPtr%var_p_2D, &
                        node%vTabPtr%idim_type, &
                        node%xStart, node%xEnd,&
                        node%yStart, node%yEnd,&
                        msgData%buf, lastBuffer)

!            print *,'LFR-DBG: inside WMS 3.1.0.1:',size(node%vTabPtr%var_p_2D,1),size(node%vTabPtr%var_p_2D,2),lastBuffer


!            print *,'LFR-DBG: inside WMS 3.1.1: ',size(g3d_g(1)%cugd_ttens,1), &
!                 size(g3d_g(1)%cugd_ttens,2),size(g3d_g(1)%cugd_ttens,3); call flush(6)





                   if (dumpLocal) then
                      write(c0,"(i8)") firstBuffer
                      write(c1,"(i8)") lastBuffer
                      write(c2,"(i8)") node%xStart
                      write(c3,"(i8)") node%xEnd
                      write(c4,"(i8)") node%yStart
                      write(c5,"(i8)") node%yEnd
                      call MsgDump(h//&
                           "filled 2D field "//trim(adjustl(node%vTabPtr%name))//"["//&
                           trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
                           trim(adjustl(c4))//":"//trim(adjustl(c5))//"]"//&
                           " with buf["//trim(adjustl(c0))//&
                           ":"//trim(adjustl(c1))//"]")
                   end if
!print *,'LFR-DBG: inside WMS 3.2: ',size(g3d_g(1)%cugd_ttens,1), &
!                 size(g3d_g(1)%cugd_ttens,2),size(g3d_g(1)%cugd_ttens,3); call flush(6)
                else if (associated(node%vTabPtr%var_p_3D)) then
                   firstBuffer = lastBuffer+1
                   !print *,'LFR-DBG: 3D: ',trim(adjustl(node%vTabPtr%name)),node%vTabPtr%idim_type; call flush(6)
                   call Buffer2FieldSection(node%vTabPtr%var_p_3D, &
                        node%vTabPtr%idim_type, &
                        node%xStart, node%xEnd,&
                        node%yStart, node%yEnd,&
                        msgData%buf, lastBuffer)

                   if (dumpLocal) then
                      write(c0,"(i8)") firstBuffer
                      write(c1,"(i8)") lastBuffer
                      write(c2,"(i8)") node%xStart
                      write(c3,"(i8)") node%xEnd
                      write(c4,"(i8)") node%yStart
                      write(c5,"(i8)") node%yEnd
                      select case (node%vTabPtr%idim_type)
                      case(3)
                         call MsgDump(h//&
                              " filled 3D field "//trim(adjustl(node%vTabPtr%name))//"[:,"//&
                              trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
                              trim(adjustl(c4))//":"//trim(adjustl(c5))//"]"//&
                              " with buf["//trim(adjustl(c0))//&
                              ":"//trim(adjustl(c1))//"]")
                      case(6:7)
                         call MsgDump(h//&
                              " filled 3D field "//trim(adjustl(node%vTabPtr%name))//"["//&
                              trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
                              trim(adjustl(c4))//":"//trim(adjustl(c5))//",:]"//&
                              " with buf["//trim(adjustl(c0))//&
                              ":"//trim(adjustl(c1))//"]")
                      end select
                   end if
!print *,'LFR-DBG: inside WMS 3.3: ',size(g3d_g(1)%cugd_ttens,1), &
!                 size(g3d_g(1)%cugd_ttens,2),size(g3d_g(1)%cugd_ttens,3); call flush(6)
                else if (associated(node%vTabPtr%var_p_4D)) then
                   firstBuffer = lastBuffer+1
                   !print *,'LFR-DBG: 4D: ',trim(adjustl(node%vTabPtr%name)),node%vTabPtr%idim_type; call flush(6)
                   call Buffer2FieldSection(node%vTabPtr%var_p_4D, &
                        node%vTabPtr%idim_type, &
                        node%xStart, node%xEnd,&
                        node%yStart, node%yEnd,&
                        msgData%buf, lastBuffer)

                   if (dumpLocal) then
                      write(c0,"(i8)") firstBuffer
                      write(c1,"(i8)") lastBuffer
                      write(c2,"(i8)") node%xStart
                      write(c3,"(i8)") node%xEnd
                      write(c4,"(i8)") node%yStart
                      write(c5,"(i8)") node%yEnd
                      call MsgDump(h//&
                           " filled 4D field "//trim(adjustl(node%vTabPtr%name))//"[:,"//&
                           trim(adjustl(c2))//":"//trim(adjustl(c3))//","//&
                           trim(adjustl(c4))//":"//trim(adjustl(c5))//",:]"//&
                           " with buf["//trim(adjustl(c0))//&
                           ":"//trim(adjustl(c1))//"]")
                   end if
!print *,'LFR-DBG: inside WMS 3.4: ',size(g3d_g(1)%cugd_ttens,1), &
!                 size(g3d_g(1)%cugd_ttens,2),size(g3d_g(1)%cugd_ttens,3); call flush(6)
                else
                   call fatal_error(h//" inconsistent var_tables entry")
                end if
             else
                exit
             end if
          end do

          ! done with this message; deallocate buffer

          deallocate(msgData%buf)
       end do
    end if
!print *,'LFR-DBG: inside WMS 4: ',size(g3d_g(1)%cugd_ttens,1), &
!                 size(g3d_g(1)%cugd_ttens,2),size(g3d_g(1)%cugd_ttens,3); call flush(6)
    ! for all posted send messages, wait on pending request,
    ! deallocate buffer and empty request

    if (associated(SendMsg)) then
!CDIR$ NOVECTOR
       do iSend = 1,SendMsg%nMsgs
          call parf_wait_any_nostatus(SendMsg%nMsgs, &
               SendMsg%request, sendNbr)
          msgData => SendMsg%oneMsg(sendNbr)
          deallocate(msgData%buf)
       end do
    end if

    !print *,'LFR-DBG: inside WMS 5: ',size(g3d_g(1)%cugd_ttens,1), &
!                 size(g3d_g(1)%cugd_ttens,2),size(g3d_g(1)%cugd_ttens,3); call flush(6)
  end subroutine WaitRecvMsgs
end module ModMessageSet
