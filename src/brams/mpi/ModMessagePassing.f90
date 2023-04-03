module ModMessagePassing

  use ModGridDims, only: &
       GridDims

  use ModParallelEnvironment, only: &
       ParallelEnvironment, &
       MsgDump

  use ModNeighbourNodes, only: &
       NeighbourNodes, &
       NodesToSendRecvMessages, &
       IncludeDomainBoundaries

  use ModDomainDecomp, only: &
       DomainDecomp

  use ModMessageSet, only: &
       MessageSet, &
       CreateMessageSet, &
       InsertFieldSectionAtMessageSet, &
       DestroyMessageSet

  use var_tables, only: &
       var_tables_r, &
       GetVTabEntry

  use ModNamelistFile, only: &
       NamelistFile

  use mem_grid, only : &
   	 dyncore_flag

  implicit none
  private
  public :: CreateAcousticMessagePassing
  public :: DestroyAcousticMessagePassing

  public :: CreateDn0MessagePassing
  public :: DestroyDn0MessagePassing

  public :: CreateG3DMessagePassing
  public :: DestroyG3DMessagePassing

  public :: CreateSelectedGhostZoneMessagePassing
  public :: DestroySelectedGhostZoneMessagePassing

  public :: CreateAllGhostZoneMessagePassing
  public :: DestroyAllGhostZoneMessagePassing
contains




  subroutine CreateAcousticMessagePassing(&
       gridId, GridSize, ParEnv, Neigh, &
       GlobalNoGhost, GlobalWithGhost, &
       AcouSendU, AcouRecvU, &
       AcouSendV, AcouRecvV, &
       AcouSendP, AcouRecvP, &
       AcouSendUV, AcouRecvUV, &
       AcouSendWP, AcouRecvWP)

    integer, intent(in) :: gridId
    type(GridDims), pointer :: GridSize
    type(ParallelEnvironment), pointer :: ParEnv    ! intent(in)
    type(NeighbourNodes), pointer :: Neigh          ! intent(in)
    type(DomainDecomp), pointer :: GlobalNoGhost    ! intent(in)
    type(DomainDecomp), pointer :: GlobalWithGhost  ! intent(in)
    type(MessageSet), pointer :: AcouSendU          ! intent(out)
    type(MessageSet), pointer :: AcouRecvU          ! intent(out)
    type(MessageSet), pointer :: AcouSendV          ! intent(out)
    type(MessageSet), pointer :: AcouRecvV          ! intent(out)
    type(MessageSet), pointer :: AcouSendP          ! intent(out)
    type(MessageSet), pointer :: AcouRecvP          ! intent(out)
    type(MessageSet), pointer :: AcouSendUV         ! intent(out)
    type(MessageSet), pointer :: AcouRecvUV         ! intent(out)
    type(MessageSet), pointer :: AcouSendWP         ! intent(out)
    type(MessageSet), pointer :: AcouRecvWP         ! intent(out)

    integer :: nMachs
    integer :: myNum
    integer :: nNeigh
    character(len=*), parameter :: h="**(CreateAcousticMessagePassing)**"

    ! verify input arguments

    if (.not. associated(ParEnv)) then
       call fatal_error(h//" starts with null ParEnv")
    else if (.not. associated(GlobalNoGhost)) then
       call fatal_error(h//" starts with null GlobalNoGhost")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" starts with null GlobalWithGhost")
    end if

    ! default output (case no neighbours)

    AcouSendU => null()
    AcouRecvU => null()
    AcouSendV => null()
    AcouRecvV => null()
    AcouSendP => null()
    AcouRecvP => null()
    AcouSendUV => null()
    AcouRecvUV => null()
    AcouSendWP => null()
    AcouRecvWP => null()

    if (associated(Neigh)) then

       myNum  = ParEnv%myNum
       nMachs = ParEnv%nMachs
       nNeigh = Neigh%nNeigh

       call CreateAcousticSendRecvU(&
            gridId, nMachs, nNeigh, myNum, &
            Neigh, GlobalNoGhost, GlobalWithGhost, &
            AcouSendU, AcouRecvU)

       call CreateAcousticSendRecvV(&
            gridId, nMachs, nNeigh, myNum, &
            Neigh, GlobalNoGhost, GlobalWithGhost, &
            AcouSendV, AcouRecvV)

       call CreateAcousticSendRecvP(&
            gridId, nMachs, nNeigh, myNum, &
            Neigh, GlobalNoGhost, GlobalWithGhost, &
            AcouSendP, AcouRecvP)

       call CreateAcousticSendRecvUV(&
            gridId, nMachs, nNeigh, myNum, &
            GridSize, Neigh, GlobalNoGhost, GlobalWithGhost, &
            AcouSendUV, AcouRecvUV)

       call CreateAcousticSendRecvWP(&
            gridId, nMachs, nNeigh, myNum, &
            GridSize, Neigh, GlobalNoGhost, GlobalWithGhost, &
            AcouSendWP, AcouRecvWP)
    end if
  end subroutine CreateAcousticMessagePassing




  subroutine CreateDn0MessagePassing(&
       gridId, GridSize, ParEnv, Neigh, &
       GlobalNoGhost, GlobalWithGhost, &
       SendDn0u, RecvDn0u, &
       SendDn0v, RecvDn0v)

    integer, intent(in) :: gridId
    type(GridDims), pointer :: GridSize
    type(ParallelEnvironment), pointer :: ParEnv    ! intent(in)
    type(NeighbourNodes), pointer :: Neigh          ! intent(in)
    type(DomainDecomp), pointer :: GlobalNoGhost    ! intent(in)
    type(DomainDecomp), pointer :: GlobalWithGhost  ! intent(in)
    type(MessageSet), pointer :: SendDn0u          ! intent(out)
    type(MessageSet), pointer :: RecvDn0u          ! intent(out)
    type(MessageSet), pointer :: SendDn0v          ! intent(out)
    type(MessageSet), pointer :: RecvDn0v          ! intent(out)

    integer :: nMachs
    integer :: myNum
    integer :: nNeigh
    character(len=*), parameter :: h="**(CreateDn0MessagePassing)**"

    ! verify input arguments

    if (.not. associated(ParEnv)) then
       call fatal_error(h//" starts with null ParEnv")
    else if (.not. associated(GlobalNoGhost)) then
       call fatal_error(h//" starts with null GlobalNoGhost")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" starts with null GlobalWithGhost")
    end if

    ! default output (case no neighbours)

    SendDn0u => null()
    RecvDn0u => null()
    SendDn0v => null()
    RecvDn0v => null()

    if (associated(Neigh)) then

       myNum  = ParEnv%myNum
       nMachs = ParEnv%nMachs
       nNeigh = Neigh%nNeigh

       call CreateSendRecvDn0u(&
            gridId, nMachs, nNeigh, myNum, &
            Neigh, GlobalNoGhost, GlobalWithGhost, &
            SendDn0u, RecvDn0u)

       call CreateSendRecvDn0v(&
            gridId, nMachs, nNeigh, myNum, &
            Neigh, GlobalNoGhost, GlobalWithGhost, &
            SendDn0v, RecvDn0v)
    end if
  end subroutine CreateDn0MessagePassing




  subroutine CreateG3DMessagePassing(&
       gridId, GridSize, ParEnv, Neigh, &
       GlobalNoGhost, GlobalWithGhost, &
       Ramsin, &
       SendG3D, RecvG3D)

    integer, intent(in) :: gridId
    type(GridDims), pointer :: GridSize
    type(ParallelEnvironment), pointer :: ParEnv    ! intent(in)
    type(NeighbourNodes), pointer :: Neigh          ! intent(in)
    type(DomainDecomp), pointer :: GlobalNoGhost    ! intent(in)
    type(DomainDecomp), pointer :: GlobalWithGhost  ! intent(in)
    type(NamelistFile), pointer :: Ramsin           ! intent(in)
    type(MessageSet), pointer :: SendG3D          ! intent(out)
    type(MessageSet), pointer :: RecvG3D          ! intent(out)

    integer :: nMachs
    integer :: myNum
    integer :: nNeigh
    integer :: g3d_spread
    integer :: g3d_smoothh
    character(len=*), parameter :: h="**(CreateG3DMessagePassing)**"

    ! verify input arguments

    if (.not. associated(ParEnv)) then
       call fatal_error(h//" starts with null ParEnv")
    else if (.not. associated(GlobalNoGhost)) then
       call fatal_error(h//" starts with null GlobalNoGhost")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" starts with null GlobalWithGhost")
    end if

    ! default output (case no neighbours or
    ! selected namelist variables were not set)

    SendG3D => null()
    RecvG3D => null()

    ! there will be messages if there are neighbour nodes
    ! and any of the namelist variables
    ! g3d_spread or g3d_smoothh were set

    g3d_spread = Ramsin%g3d_spread
    g3d_smoothh = Ramsin%g3d_smoothh

    if (associated(Neigh) .and. &
         (g3d_spread /= 0 .or. g3d_smoothh /= 0)) then

       myNum  = ParEnv%myNum
       nMachs = ParEnv%nMachs
       nNeigh = Neigh%nNeigh

       call CreateG3DSendRecv(&
            gridId, nMachs, nNeigh, myNum, &
            g3d_spread, g3d_smoothh, &
            GridSize, Neigh, GlobalNoGhost, GlobalWithGhost, &
            SendG3D, RecvG3D)
    end if
  end subroutine CreateG3DMessagePassing





  subroutine CreateSelectedGhostZoneMessagePassing(&
       gridId, num_var, vtab_r, &
       GridSize, ParEnv, Neigh, &
       GlobalNoGhost, GlobalWithGhost, &
       SelectedGhostZoneSend, SelectedGhostZoneRecv)

    integer, intent(in) :: gridId
    integer, intent(in) :: num_var(:)
    type(var_tables_r), target, intent(in) ::  vtab_r(:,:)
    type(GridDims), pointer :: GridSize
    type(ParallelEnvironment), pointer :: ParEnv    ! intent(in)
    type(NeighbourNodes), pointer :: Neigh          ! intent(in)
    type(DomainDecomp), pointer :: GlobalNoGhost    ! intent(in)
    type(DomainDecomp), pointer :: GlobalWithGhost  ! intent(in)
    type(MessageSet), pointer :: SelectedGhostZoneSend          ! intent(out)
    type(MessageSet), pointer :: SelectedGhostZoneRecv          ! intent(out)


    integer :: nMachs
    integer :: myNum
    integer :: nNeigh
    character(len=*), parameter :: h="**(CreateAcousticMessagePassing)**"

    ! verify input arguments

    if (.not. associated(ParEnv)) then
       call fatal_error(h//" starts with null ParEnv")
    else if (.not. associated(GlobalNoGhost)) then
       call fatal_error(h//" starts with null GlobalNoGhost")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" starts with null GlobalWithGhost")
    end if

    ! default output (case no neighbours)

    SelectedGhostZoneSend => null()
    SelectedGhostZoneRecv => null()

    if (associated(Neigh)) then

       myNum  = ParEnv%myNum
       nMachs = ParEnv%nMachs
       nNeigh = Neigh%nNeigh

       call CreateSelectedGhostZoneSendRecv(&
            gridId, nMachs, nNeigh, myNum, num_var, vtab_r, &
            GridSize, Neigh, GlobalNoGhost, GlobalWithGhost, &
            SelectedGhostZoneSend, SelectedGhostZoneRecv)
    end if
  end subroutine CreateSelectedGhostZoneMessagePassing





  subroutine CreateAllGhostZoneMessagePassing(&
       gridId, num_var, vtab_r, &
       GridSize, ParEnv, Neigh, &
       GlobalNoGhost, GlobalWithGhost, &
       AllGhostZoneSend, AllGhostZoneRecv)

    integer, intent(in) :: gridId
    integer, intent(in) :: num_var(:)
    type(var_tables_r), target, intent(in) ::  vtab_r(:,:)
    type(GridDims), pointer :: GridSize
    type(ParallelEnvironment), pointer :: ParEnv    ! intent(in)
    type(NeighbourNodes), pointer :: Neigh          ! intent(in)
    type(DomainDecomp), pointer :: GlobalNoGhost    ! intent(in)
    type(DomainDecomp), pointer :: GlobalWithGhost  ! intent(in)
    type(MessageSet), pointer :: AllGhostZoneSend          ! intent(out)
    type(MessageSet), pointer :: AllGhostZoneRecv          ! intent(out)


    integer :: nMachs
    integer :: myNum
    integer :: nNeigh
    character(len=*), parameter :: h="**(CreateAcousticMessagePassing)**"

    ! verify input arguments

    if (.not. associated(ParEnv)) then
       call fatal_error(h//" starts with null ParEnv")
    else if (.not. associated(GlobalNoGhost)) then
       call fatal_error(h//" starts with null GlobalNoGhost")
    else if (.not. associated(GlobalWithGhost)) then
       call fatal_error(h//" starts with null GlobalWithGhost")
    end if

    ! default output (case no neighbours)

    AllGhostZoneSend => null()
    AllGhostZoneRecv => null()

    if (associated(Neigh)) then

       myNum  = ParEnv%myNum
       nMachs = ParEnv%nMachs
       nNeigh = Neigh%nNeigh

       call CreateAllGhostZoneSendRecv(&
            gridId, nMachs, nNeigh, myNum, num_var, vtab_r, &
            GridSize, Neigh, GlobalNoGhost, GlobalWithGhost, &
            AllGhostZoneSend, AllGhostZoneRecv)
    end if
  end subroutine CreateAllGhostZoneMessagePassing





  subroutine DestroyAcousticMessagePassing( &
       AcouSendU, AcouRecvU, &
       AcouSendV, AcouRecvV, &
       AcouSendP, AcouRecvP, &
       AcouSendUV, AcouRecvUV, &
       AcouSendWP, AcouRecvWP)

    type(MessageSet), pointer :: AcouSendU            ! intent(out)
    type(MessageSet), pointer :: AcouRecvU            ! intent(out)
    type(MessageSet), pointer :: AcouSendV            ! intent(out)
    type(MessageSet), pointer :: AcouRecvV            ! intent(out)
    type(MessageSet), pointer :: AcouSendP            ! intent(out)
    type(MessageSet), pointer :: AcouRecvP            ! intent(out)
    type(MessageSet), pointer :: AcouSendUV           ! intent(out)
    type(MessageSet), pointer :: AcouRecvUV           ! intent(out)
    type(MessageSet), pointer :: AcouSendWP           ! intent(out)
    type(MessageSet), pointer :: AcouRecvWP           ! intent(out)

    call DestroyMessageSet(AcouSendU)
    call DestroyMessageSet(AcouRecvU)

    call DestroyMessageSet(AcouSendV)
    call DestroyMessageSet(AcouRecvV)

    call DestroyMessageSet(AcouSendP)
    call DestroyMessageSet(AcouRecvP)

    call DestroyMessageSet(AcouSendUV)
    call DestroyMessageSet(AcouRecvUV)

    call DestroyMessageSet(AcouSendWP)
    call DestroyMessageSet(AcouRecvWP)
  end subroutine DestroyAcousticMessagePassing




  subroutine DestroyDn0MessagePassing( &
       SendDn0u, RecvDn0u, SendDn0v, RecvDn0v)

    type(MessageSet), pointer :: SendDn0u
    type(MessageSet), pointer :: RecvDn0u
    type(MessageSet), pointer :: SendDn0v
    type(MessageSet), pointer :: RecvDn0v

    call DestroyMessageSet(SendDn0u)
    call DestroyMessageSet(RecvDn0u)
    call DestroyMessageSet(SendDn0v)
    call DestroyMessageSet(RecvDn0v)

  end subroutine DestroyDn0MessagePassing




  subroutine DestroyG3DMessagePassing( &
       SendG3D, RecvG3D)

    type(MessageSet), pointer :: SendG3D
    type(MessageSet), pointer :: RecvG3D

    call DestroyMessageSet(SendG3D)
    call DestroyMessageSet(RecvG3D)

  end subroutine DestroyG3DMessagePassing




  subroutine DestroySelectedGhostZoneMessagePassing( &
       SelectedGhostZoneSend, SelectedGhostZoneRecv)

    type(MessageSet), pointer :: SelectedGhostZoneSend            ! intent(out)
    type(MessageSet), pointer :: SelectedGhostZoneRecv            ! intent(out)

    call DestroyMessageSet(SelectedGhostZoneSend)
    call DestroyMessageSet(SelectedGhostZoneRecv)

  end subroutine DestroySelectedGhostZoneMessagePassing






  subroutine DestroyAllGhostZoneMessagePassing( &
       AllGhostZoneSend, AllGhostZoneRecv)

    type(MessageSet), pointer :: AllGhostZoneSend            ! intent(out)
    type(MessageSet), pointer :: AllGhostZoneRecv            ! intent(out)

    call DestroyMessageSet(AllGhostZoneSend)
    call DestroyMessageSet(AllGhostZoneRecv)

  end subroutine DestroyAllGhostZoneMessagePassing




  subroutine CreateAcousticSendRecvU(&
       gridId, nMachs, nNeigh, myNum, &
       Neigh, GlobalNoGhost, GlobalWithGhost, &
       AcouSendU, AcouRecvU)

    integer, intent(in) :: gridId
    integer, intent(in) :: nMachs
    integer, intent(in) :: nNeigh
    integer, intent(in) :: myNum

    type(NeighbourNodes), pointer :: Neigh          ! intent(in)
    type(DomainDecomp), pointer :: GlobalNoGhost    ! intent(in)
    type(DomainDecomp), pointer :: GlobalWithGhost  ! intent(in)
    type(MessageSet), pointer :: AcouSendU          ! intent(out)
    type(MessageSet), pointer :: AcouRecvU          ! intent(out)

    integer, parameter :: TagU=22
    character(len=*), parameter :: NameSendU="AcouSendU"
    character(len=*), parameter :: NameRecvU="AcouRecvU"

    ! scratch arrays of size number of BRAMS processes
    ! containing global indices

    integer :: xbToBeUpdated(nMachs)
    integer :: xeToBeUpdated(nMachs)
    integer :: ybToBeUpdated(nMachs)
    integer :: yeToBeUpdated(nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(nNeigh)
    integer :: xeSend(nNeigh)
    integer :: ybSend(nNeigh)
    integer :: yeSend(nNeigh)
    integer :: xbRecv(nNeigh)
    integer :: xeRecv(nNeigh)
    integer :: ybRecv(nNeigh)
    integer :: yeRecv(nNeigh)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(nNeigh)
    logical :: willRecv(nNeigh)

    type(var_tables_r), pointer   :: vtabPtr => null()
    character(len=*), parameter :: h="**(CreateAcousticSendRecvU)**"
    character(len=30) :: tmp_name

    ! AcouSendU, AcouRecvU:
    ! messages update GlobalNoGhost [xb-1:xb-1,yb:ye]

    xbToBeUpdated = GlobalNoGhost%xb - 1
    xeToBeUpdated = xbToBeUpdated
    ybToBeUpdated = GlobalNoGhost%yb
    yeToBeUpdated = GlobalNoGhost%ye

    ! which neighbour nodes will send and receive

    call NodesToSendRecvMessages(myNum, Neigh, GlobalNoGhost, &
         xbToBeUpdated, xeToBeUpdated, ybToBeUpdated, yeToBeUpdated, &
         xbSend, xeSend, ybSend, yeSend, willSend, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv)

    ! build message set

    AcouSendU => CreateMessageSet(NameSendU, TagU, willSend, Neigh)
    AcouRecvU => CreateMessageSet(NameRecvU, TagU, willRecv, Neigh)

    ! get field

    vTabPtr => null()
    if(dyncore_flag==2) then
      tmp_name='UC'
    else
      tmp_name='UP'
    endif
    call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

    ! build field sections to be sent and received

    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcouSendU)
    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcouRecvU)
  end subroutine CreateAcousticSendRecvU





  subroutine CreateAcousticSendRecvV(&
       gridId, nMachs, nNeigh, myNum, &
       Neigh, GlobalNoGhost, GlobalWithGhost, &
       AcouSendV, AcouRecvV)

    integer, intent(in) :: gridId
    integer, intent(in) :: nMachs
    integer, intent(in) :: nNeigh
    integer, intent(in) :: myNum

    type(NeighbourNodes), pointer :: Neigh          ! intent(in)
    type(DomainDecomp), pointer :: GlobalNoGhost    ! intent(in)
    type(DomainDecomp), pointer :: GlobalWithGhost  ! intent(in)
    type(MessageSet), pointer :: AcouSendV          ! intent(out)
    type(MessageSet), pointer :: AcouRecvV          ! intent(out)

    integer, parameter :: TagV=23
    character(len=*), parameter :: NameSendV="AcouSendV"
    character(len=*), parameter :: NameRecvV="AcouRecvV"
    character(len=30) :: tmp_name
    ! scratch arrays of size number of BRAMS processes
    ! containing global indices

    integer :: xbToBeUpdated(nMachs)
    integer :: xeToBeUpdated(nMachs)
    integer :: ybToBeUpdated(nMachs)
    integer :: yeToBeUpdated(nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(nNeigh)
    integer :: xeSend(nNeigh)
    integer :: ybSend(nNeigh)
    integer :: yeSend(nNeigh)
    integer :: xbRecv(nNeigh)
    integer :: xeRecv(nNeigh)
    integer :: ybRecv(nNeigh)
    integer :: yeRecv(nNeigh)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(nNeigh)
    logical :: willRecv(nNeigh)


    type(var_tables_r), pointer   :: vtabPtr => null()
    character(len=*), parameter :: h="**(CreateAcousticSendRecvV)**"

    ! AcouSendV, AcouRecvV:
    ! messages update GlobalNoGhost [xb:xe,yb-1:yb-1]

    xbToBeUpdated = GlobalNoGhost%xb
    xeToBeUpdated = GlobalNoGhost%xe
    ybToBeUpdated = GlobalNoGhost%yb - 1
    yeToBeUpdated = ybToBeUpdated

    ! which neighbour nodes will send and receive

    call NodesToSendRecvMessages(myNum, Neigh, GlobalNoGhost, &
         xbToBeUpdated, xeToBeUpdated, ybToBeUpdated, yeToBeUpdated, &
         xbSend, xeSend, ybSend, yeSend, willSend, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv)

    ! build message set

    AcouSendV => CreateMessageSet(NameSendV, TagV, willSend, Neigh)
    AcouRecvV => CreateMessageSet(NameRecvV, TagV, willRecv, Neigh)

    ! get field

    vTabPtr => null()
    if(dyncore_flag==2) then
      tmp_name='VC'
    else
      tmp_name='VP'
    endif
    call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

    ! build field sections to be sent and received

    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcouSendV)
    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcouRecvV)
  end subroutine CreateAcousticSendRecvV





  subroutine CreateAcousticSendRecvP(&
       gridId, nMachs, nNeigh, myNum, &
       Neigh, GlobalNoGhost, GlobalWithGhost, &
       AcouSendP, AcouRecvP)

    integer, intent(in) :: gridId
    integer, intent(in) :: nMachs
    integer, intent(in) :: nNeigh
    integer, intent(in) :: myNum

    type(NeighbourNodes), pointer :: Neigh          ! intent(in)
    type(DomainDecomp), pointer :: GlobalNoGhost    ! intent(in)
    type(DomainDecomp), pointer :: GlobalWithGhost  ! intent(in)
    type(MessageSet), pointer :: AcouSendP          ! intent(out)
    type(MessageSet), pointer :: AcouRecvP          ! intent(out)

    integer, parameter :: TagP=24
    character(len=*), parameter :: NameSendP="AcouSendP"
    character(len=*), parameter :: NameRecvP="AcouRecvP"

    ! scratch arrays of size number of BRAMS processes
    ! containing global indices

    integer :: xbToBeUpdated(nMachs)
    integer :: xeToBeUpdated(nMachs)
    integer :: ybToBeUpdated(nMachs)
    integer :: yeToBeUpdated(nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(nNeigh)
    integer :: xeSend(nNeigh)
    integer :: ybSend(nNeigh)
    integer :: yeSend(nNeigh)
    integer :: xbRecv(nNeigh)
    integer :: xeRecv(nNeigh)
    integer :: ybRecv(nNeigh)
    integer :: yeRecv(nNeigh)

    ! second region to build union of regions

    integer :: xbSend_2(nNeigh)
    integer :: xeSend_2(nNeigh)
    integer :: ybSend_2(nNeigh)
    integer :: yeSend_2(nNeigh)
    integer :: xbRecv_2(nNeigh)
    integer :: xeRecv_2(nNeigh)
    integer :: ybRecv_2(nNeigh)
    integer :: yeRecv_2(nNeigh)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(nNeigh)
    logical :: willRecv(nNeigh)

    logical :: willSend_2(nNeigh)
    logical :: willRecv_2(nNeigh)

    type(var_tables_r), pointer   :: vtabPtr => null()
    character(len=*), parameter :: h="**(CreateAcousticSendRecvP)**"
    character(len=30) :: tmp_name

    ! AcouSendP, AcouRecvP: union of
    !               GlobalNoGhost [xe+1:xe+1,yb:ye] with
    !               GlobalNoGhost [xb:xe,ye+1:ye+1]

    ! first part of the union ([xe+1:xe+1,yb:ye])

    xbToBeUpdated = GlobalNoGhost%xe+1
    xeToBeUpdated = xbToBeUpdated
    ybToBeUpdated = GlobalNoGhost%yb
    yeToBeUpdated = GlobalNoGhost%ye

    ! which neighbour nodes will send and receive

    call NodesToSendRecvMessages(myNum, Neigh, GlobalNoGhost, &
         xbToBeUpdated, xeToBeUpdated, ybToBeUpdated, yeToBeUpdated, &
         xbSend, xeSend, ybSend, yeSend, willSend, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv)

    ! second part of the union ([xb:xe,ye+1:ye+1])

    xbToBeUpdated = GlobalNoGhost%xb
    xeToBeUpdated = GlobalNoGhost%xe
    ybToBeUpdated = GlobalNoGhost%ye + 1
    yeToBeUpdated = ybToBeUpdated


    ! which neighbour nodes will send and receive

    call NodesToSendRecvMessages(myNum, Neigh, GlobalNoGhost, &
         xbToBeUpdated, xeToBeUpdated, ybToBeUpdated, yeToBeUpdated, &
         xbSend_2, xeSend_2, ybSend_2, yeSend_2, willSend_2, &
         xbRecv_2, xeRecv_2, ybRecv_2, yeRecv_2, willRecv_2)

    ! make union

    call BuildUnion(nNeigh, &
         xbSend_2, xeSend_2, ybSend_2, yeSend_2, willSend_2, &
         xbSend, xeSend, ybSend, yeSend, willSend)
    call BuildUnion(nNeigh, &
         xbRecv_2, xeRecv_2, ybRecv_2, yeRecv_2, willRecv_2, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv)

    ! build message set

    AcouSendP => CreateMessageSet(NameSendP, TagP, willSend, Neigh)
    AcouRecvP => CreateMessageSet(NameRecvP, TagP, willRecv, Neigh)

    ! get field

    vTabPtr => null()
    if(dyncore_flag==2) then
      tmp_name='PC'
    else
      tmp_name='PP'
    endif
    call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

    ! build field sections to be sent and received

    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcouSendP)
    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcouRecvP)
  end subroutine CreateAcousticSendRecvP





  subroutine BuildUnion(nNeigh, &
       xbComm1, xeComm1, ybComm1, yeComm1, willComm1, &
       xbComm2, xeComm2, ybComm2, yeComm2, willComm2)
    integer, intent(in) :: nNeigh
    integer, intent(in) :: xbComm1(nNeigh)
    integer, intent(in) :: xeComm1(nNeigh)
    integer, intent(in) :: ybComm1(nNeigh)
    integer, intent(in) :: yeComm1(nNeigh)
    logical, intent(in) :: willComm1(nNeigh)

    integer, intent(inout) :: xbComm2(nNeigh)
    integer, intent(inout) :: xeComm2(nNeigh)
    integer, intent(inout) :: ybComm2(nNeigh)
    integer, intent(inout) :: yeComm2(nNeigh)
    logical, intent(inout) :: willComm2(nNeigh)

    integer :: iNeigh
    character(len=8) :: c0, c1, c2
    character(len=128) :: inter1, inter2
    character(len=*), parameter :: h="**(BuildUnion)**"

    do iNeigh = 1, nNeigh
       if (willComm1(iNeigh) .and. willComm2(iNeigh)) then
          write(c0,"(i8)") iNeigh
          inter1="inter1"
          write(c1,"(i8)") xbComm1
          write(c2,"(i8)") xeComm1
          inter1=trim(inter1)//"("//trim(adjustl(c1))//":"//trim(adjustl(c2))
          write(c1,"(i8)") ybComm1
          write(c2,"(i8)") yeComm1
          inter1=trim(inter1)//","//trim(adjustl(c1))//":"//trim(adjustl(c2))//")"
          inter2="inter2"
          write(c1,"(i8)") xbComm2
          write(c2,"(i8)") xeComm2
          inter2=trim(inter2)//"("//trim(adjustl(c1))//":"//trim(adjustl(c2))
          write(c1,"(i8)") ybComm2
          write(c2,"(i8)") yeComm2
          inter2=trim(inter2)//","//trim(adjustl(c1))//":"//trim(adjustl(c2))//")"
          call fatal_error(h//" both willComm("//trim(adjustl(c0))//&
               "); intervals="//trim(adjustl(inter1))//", "//trim(adjustl(inter2)))
       else if (willComm1(iNeigh)) then
          xbComm2(iNeigh) = xbComm1(iNeigh)
          xeComm2(iNeigh) = xeComm1(iNeigh)
          ybComm2(iNeigh) = ybComm1(iNeigh)
          yeComm2(iNeigh) = yeComm1(iNeigh)
          willComm2(iNeigh) = willComm1(iNeigh)
       end if
    end do
  end subroutine BuildUnion




  subroutine CreateAcousticSendRecvUV(&
       gridId, nMachs, nNeigh, myNum, &
       GridSize, Neigh, GlobalNoGhost, GlobalWithGhost, &
       AcouSendUV, AcouRecvUV)

    integer, intent(in) :: gridId
    integer, intent(in) :: nMachs
    integer, intent(in) :: nNeigh
    integer, intent(in) :: myNum

    type(GridDims), pointer :: GridSize
    type(NeighbourNodes), pointer :: Neigh          ! intent(in)
    type(DomainDecomp), pointer :: GlobalNoGhost    ! intent(in)
    type(DomainDecomp), pointer :: GlobalWithGhost  ! intent(in)
    type(MessageSet), pointer :: AcouSendUV          ! intent(out)
    type(MessageSet), pointer :: AcouRecvUV          ! intent(out)

    integer, parameter :: TagUV=25
    character(len=*), parameter :: NameSendUV="AcouSendUV"
    character(len=*), parameter :: NameRecvUV="AcouRecvUV"

    ! scratch arrays of size number of BRAMS processes
    ! containing global indices

    integer :: xbToBeUpdated(nMachs)
    integer :: xeToBeUpdated(nMachs)
    integer :: ybToBeUpdated(nMachs)
    integer :: yeToBeUpdated(nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(nNeigh)
    integer :: xeSend(nNeigh)
    integer :: ybSend(nNeigh)
    integer :: yeSend(nNeigh)
    integer :: xbRecv(nNeigh)
    integer :: xeRecv(nNeigh)
    integer :: ybRecv(nNeigh)
    integer :: yeRecv(nNeigh)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(nNeigh)
    logical :: willRecv(nNeigh)

    type(var_tables_r), pointer   :: vtabPtr => null()
    character(len=*), parameter :: h="**(CreateAcousticSendRecvUV)**"
    character(len=30) :: tmp_name

    ! AcouSendUV, AcouRecvUV:
    ! messages update entire GhostZone

    xbToBeUpdated = GlobalWithGhost%xb
    xeToBeUpdated = GlobalWithGhost%xe
    ybToBeUpdated = GlobalWithGhost%yb
    yeToBeUpdated = GlobalWithGhost%ye

    ! which neighbour nodes will send and receive

    call NodesToSendRecvMessages(myNum, Neigh, GlobalNoGhost, &
         xbToBeUpdated, xeToBeUpdated, ybToBeUpdated, yeToBeUpdated, &
         xbSend, xeSend, ybSend, yeSend, willSend, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv)

    ! include real domain boundaries

    call IncludeDomainBoundaries(Neigh, GridSize, GlobalNoGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend)
    call IncludeDomainBoundaries(Neigh, GridSize, GlobalNoGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv)

    ! build message set

    AcouSendUV => CreateMessageSet(NameSendUV, TagUV, willSend, Neigh)
    AcouRecvUV => CreateMessageSet(NameRecvUV, TagUV, willRecv, Neigh)

    ! get field UP

    vTabPtr => null()
    if(dyncore_flag==2) then
      tmp_name='UC'
    else
      tmp_name='UP'
    endif
    call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

    ! build field sections to be sent and received

    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcouSendUV)
    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcouRecvUV)

    ! get field VP

    vTabPtr => null()
    if(dyncore_flag==2) then
      tmp_name='VC'
    else
      tmp_name='VP'
    endif
    call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

    ! build field sections to be sent and received

    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcouSendUV)
    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcouRecvUV)
  end subroutine CreateAcousticSendRecvUV










  subroutine CreateAcousticSendRecvWP(&
       gridId, nMachs, nNeigh, myNum, &
       GridSize, Neigh, GlobalNoGhost, GlobalWithGhost, &
       AcouSendWP, AcouRecvWP)

    integer, intent(in) :: gridId
    integer, intent(in) :: nMachs
    integer, intent(in) :: nNeigh
    integer, intent(in) :: myNum

    type(GridDims), pointer :: GridSize
    type(NeighbourNodes), pointer :: Neigh          ! intent(in)
    type(DomainDecomp), pointer :: GlobalNoGhost    ! intent(in)
    type(DomainDecomp), pointer :: GlobalWithGhost  ! intent(in)
    type(MessageSet), pointer :: AcouSendWP          ! intent(out)
    type(MessageSet), pointer :: AcouRecvWP          ! intent(out)

    integer, parameter :: TagWP=26
    character(len=*), parameter :: NameSendWP="AcouSendWP"
    character(len=*), parameter :: NameRecvWP="AcouRecvWP"

    ! scratch arrays of size number of BRAMS processes
    ! containing global indices

    integer :: xbToBeUpdated(nMachs)
    integer :: xeToBeUpdated(nMachs)
    integer :: ybToBeUpdated(nMachs)
    integer :: yeToBeUpdated(nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(nNeigh)
    integer :: xeSend(nNeigh)
    integer :: ybSend(nNeigh)
    integer :: yeSend(nNeigh)
    integer :: xbRecv(nNeigh)
    integer :: xeRecv(nNeigh)
    integer :: ybRecv(nNeigh)
    integer :: yeRecv(nNeigh)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(nNeigh)
    logical :: willRecv(nNeigh)

    type(var_tables_r), pointer   :: vtabPtr => null()
    character(len=*), parameter :: h="**(CreateAcousticSendRecvWP)**"
    character(len=30) :: tmp_name

    ! AcouSendWP, AcouRecvWP:
    ! messages update entire GhostZone

    xbToBeUpdated = GlobalWithGhost%xb
    xeToBeUpdated = GlobalWithGhost%xe
    ybToBeUpdated = GlobalWithGhost%yb
    yeToBeUpdated = GlobalWithGhost%ye

    ! which neighbour nodes will send and receive

    call NodesToSendRecvMessages(myNum, Neigh, GlobalNoGhost, &
         xbToBeUpdated, xeToBeUpdated, ybToBeUpdated, yeToBeUpdated, &
         xbSend, xeSend, ybSend, yeSend, willSend, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv)

    ! include real domain boundaries

    call IncludeDomainBoundaries(Neigh, GridSize, GlobalNoGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend)
    call IncludeDomainBoundaries(Neigh, GridSize, GlobalNoGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv)

    ! build message set

    AcouSendWP => CreateMessageSet(NameSendWP, TagWP, willSend, Neigh)
    AcouRecvWP => CreateMessageSet(NameRecvWP, TagWP, willRecv, Neigh)

    ! get field UP

    vTabPtr => null()
    if(dyncore_flag==2) then
      tmp_name='WC'
    else
      tmp_name='WP'
    endif
    call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

    ! build field sections to be sent and received

    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcouSendWP)
    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcouRecvWP)

    ! get field VP

    vTabPtr => null()
    if(dyncore_flag==2) then
      tmp_name='PC'
    else
      tmp_name='PP'
    endif
    call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

    ! build field sections to be sent and received

    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, AcouSendWP)
    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AcouRecvWP)
  end subroutine CreateAcousticSendRecvWP






  subroutine CreateSelectedGhostZoneSendRecv(&
       gridId, nMachs, nNeigh, myNum, num_var, vtab_r, &
       GridSize, Neigh, GlobalNoGhost, GlobalWithGhost, &
       SelectedGhostZoneSend, SelectedGhostZoneRecv)

    integer, intent(in) :: gridId
    integer, intent(in) :: nMachs
    integer, intent(in) :: nNeigh
    integer, intent(in) :: myNum
    integer, intent(in) :: num_var(:)
    type(var_tables_r), target, intent(in) ::  vtab_r(:,:)

    type(GridDims), pointer :: GridSize
    type(NeighbourNodes), pointer :: Neigh          ! intent(in)
    type(DomainDecomp), pointer :: GlobalNoGhost    ! intent(in)
    type(DomainDecomp), pointer :: GlobalWithGhost  ! intent(in)
    type(MessageSet), pointer :: SelectedGhostZoneSend          ! intent(out)
    type(MessageSet), pointer :: SelectedGhostZoneRecv          ! intent(out)

    integer :: vTabNbr
    integer, parameter :: TagSelectedGhostZone=27
    character(len=*), parameter :: NameSendSelectedGhostZone="SelectedGhostZoneSend"
    character(len=*), parameter :: NameRecvSelectedGhostZone="SelectedGhostZoneRecv"

    ! scratch arrays of size number of BRAMS processes
    ! containing global indices

    integer :: xbToBeUpdated(nMachs)
    integer :: xeToBeUpdated(nMachs)
    integer :: ybToBeUpdated(nMachs)
    integer :: yeToBeUpdated(nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(nNeigh)
    integer :: xeSend(nNeigh)
    integer :: ybSend(nNeigh)
    integer :: yeSend(nNeigh)
    integer :: xbRecv(nNeigh)
    integer :: xeRecv(nNeigh)
    integer :: ybRecv(nNeigh)
    integer :: yeRecv(nNeigh)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(nNeigh)
    logical :: willRecv(nNeigh)

    type(var_tables_r), pointer   :: vtabPtr => null()
    character(len=*), parameter :: h="**(CreateSelectedGhostZoneSendRecv)**"

    ! SelectedGhostZoneSend, SelectedGhostZoneRecv:
    ! messages update entire GhostZone

    xbToBeUpdated = GlobalWithGhost%xb
    xeToBeUpdated = GlobalWithGhost%xe
    ybToBeUpdated = GlobalWithGhost%yb
    yeToBeUpdated = GlobalWithGhost%ye

    ! which neighbour nodes will send and receive

    call NodesToSendRecvMessages(myNum, Neigh, GlobalNoGhost, &
         xbToBeUpdated, xeToBeUpdated, ybToBeUpdated, yeToBeUpdated, &
         xbSend, xeSend, ybSend, yeSend, willSend, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv)

    ! include real domain boundaries

    call IncludeDomainBoundaries(Neigh, GridSize, GlobalNoGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend)
    call IncludeDomainBoundaries(Neigh, GridSize, GlobalNoGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv)

    ! build message set

    SelectedGhostZoneSend => CreateMessageSet(&
         NameSendSelectedGhostZone, TagSelectedGhostZone, &
         willSend, Neigh)
    SelectedGhostZoneRecv => CreateMessageSet(&
         NameRecvSelectedGhostZone, TagSelectedGhostZone, &
         willRecv, Neigh)

    ! take all var_tables field that should be communicated

    do vTabNbr = 1, num_var(gridId)

       if (vtab_r(vTabNbr,gridId)%impt1 == 1) then
          vTabPtr => vtab_r(vTabNbr,gridId)

          ! build field sections to be sent and received

          call InsertFieldSectionAtMessageSet(&
               myNum, vTabPtr, Neigh, GlobalWithGhost, &
               xbSend, xeSend, ybSend, yeSend, willSend, &
               SelectedGhostZoneSend)
          call InsertFieldSectionAtMessageSet(&
               myNum, vTabPtr, Neigh, GlobalWithGhost, &
               xbRecv, xeRecv, ybRecv, yeRecv, willRecv, &
               SelectedGhostZoneRecv)
       end if
    end do
  end subroutine CreateSelectedGhostZoneSendRecv






  subroutine CreateAllGhostZoneSendRecv(&
       gridId, nMachs, nNeigh, myNum, num_var, vtab_r, &
       GridSize, Neigh, GlobalNoGhost, GlobalWithGhost, &
       AllGhostZoneSend, AllGhostZoneRecv)

    integer, intent(in) :: gridId
    integer, intent(in) :: nMachs
    integer, intent(in) :: nNeigh
    integer, intent(in) :: myNum
    integer, intent(in) :: num_var(:)
    type(var_tables_r), target, intent(in) ::  vtab_r(:,:)

    type(GridDims), pointer :: GridSize
    type(NeighbourNodes), pointer :: Neigh          ! intent(in)
    type(DomainDecomp), pointer :: GlobalNoGhost    ! intent(in)
    type(DomainDecomp), pointer :: GlobalWithGhost  ! intent(in)
    type(MessageSet), pointer :: AllGhostZoneSend          ! intent(out)
    type(MessageSet), pointer :: AllGhostZoneRecv          ! intent(out)

    integer :: vTabNbr
    integer, parameter :: TagAllGhostZone=28
    character(len=*), parameter :: NameSendAllGhostZone="AllGhostZoneSend"
    character(len=*), parameter :: NameRecvAllGhostZone="AllGhostZoneRecv"

    ! scratch arrays of size number of BRAMS processes
    ! containing global indices

    integer :: xbToBeUpdated(nMachs)
    integer :: xeToBeUpdated(nMachs)
    integer :: ybToBeUpdated(nMachs)
    integer :: yeToBeUpdated(nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(nNeigh)
    integer :: xeSend(nNeigh)
    integer :: ybSend(nNeigh)
    integer :: yeSend(nNeigh)
    integer :: xbRecv(nNeigh)
    integer :: xeRecv(nNeigh)
    integer :: ybRecv(nNeigh)
    integer :: yeRecv(nNeigh)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(nNeigh)
    logical :: willRecv(nNeigh)

    type(var_tables_r), pointer   :: vtabPtr => null()
    character(len=*), parameter :: h="**(CreateAllGhostZoneSendRecv)**"

    ! AllGhostZoneSend, AllGhostZoneRecv:
    ! messages update entire GhostZone

    xbToBeUpdated = GlobalWithGhost%xb
    xeToBeUpdated = GlobalWithGhost%xe
    ybToBeUpdated = GlobalWithGhost%yb
    yeToBeUpdated = GlobalWithGhost%ye

    ! which neighbour nodes will send and receive

    call NodesToSendRecvMessages(myNum, Neigh, GlobalNoGhost, &
         xbToBeUpdated, xeToBeUpdated, ybToBeUpdated, yeToBeUpdated, &
         xbSend, xeSend, ybSend, yeSend, willSend, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv)

    ! include real domain boundaries

    call IncludeDomainBoundaries(Neigh, GridSize, GlobalNoGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend)
    call IncludeDomainBoundaries(Neigh, GridSize, GlobalNoGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv)

    ! build message set

    AllGhostZoneSend => CreateMessageSet(NameSendAllGhostZone, TagAllGhostZone, willSend, Neigh)
    AllGhostZoneRecv => CreateMessageSet(NameRecvAllGhostZone, TagAllGhostZone, willRecv, Neigh)

    ! take all var_tables field that should be communicated

    do vTabNbr = 1, num_var(gridId)

       vTabPtr => vtab_r(vTabNbr,gridId)

       if (&
            trim(adjustl(vTabPtr%name)) /= "LPU" .and. &
            trim(adjustl(vTabPtr%name)) /= "LPV" .and. &
            trim(adjustl(vTabPtr%name)) /= "LPW" ) then

          ! build field sections to be sent and received

          call InsertFieldSectionAtMessageSet(&
               myNum, vTabPtr, Neigh, GlobalWithGhost, &
               xbSend, xeSend, ybSend, yeSend, willSend, AllGhostZoneSend)
          call InsertFieldSectionAtMessageSet(&
               myNum, vTabPtr, Neigh, GlobalWithGhost, &
               xbRecv, xeRecv, ybRecv, yeRecv, willRecv, AllGhostZoneRecv)
       end if
    end do
  end subroutine CreateAllGhostZoneSendRecv





  subroutine CreateSendRecvDn0u(&
       gridId, nMachs, nNeigh, myNum, &
       Neigh, GlobalNoGhost, GlobalWithGhost, &
       SendDn0u, RecvDn0u)

    integer, intent(in) :: gridId
    integer, intent(in) :: nMachs
    integer, intent(in) :: nNeigh
    integer, intent(in) :: myNum

    type(NeighbourNodes), pointer :: Neigh          ! intent(in)
    type(DomainDecomp), pointer :: GlobalNoGhost    ! intent(in)
    type(DomainDecomp), pointer :: GlobalWithGhost  ! intent(in)
    type(MessageSet), pointer :: SendDn0u          ! intent(out)
    type(MessageSet), pointer :: RecvDn0u          ! intent(out)

    integer, parameter :: TagDn0u=29
    character(len=*), parameter :: NameSendDn0u="SendDn0u"
    character(len=*), parameter :: NameRecvDn0u="RecvDn0u"

    ! scratch arrays of size number of BRAMS processes
    ! containing global indices

    integer :: xbToBeUpdated(nMachs)
    integer :: xeToBeUpdated(nMachs)
    integer :: ybToBeUpdated(nMachs)
    integer :: yeToBeUpdated(nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(nNeigh)
    integer :: xeSend(nNeigh)
    integer :: ybSend(nNeigh)
    integer :: yeSend(nNeigh)
    integer :: xbRecv(nNeigh)
    integer :: xeRecv(nNeigh)
    integer :: ybRecv(nNeigh)
    integer :: yeRecv(nNeigh)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(nNeigh)
    logical :: willRecv(nNeigh)

    type(var_tables_r), pointer   :: vtabPtr => null()
    character(len=*), parameter :: h="**(CreateSendRecvDn0u)**"
    character(len=30) :: tmp_name

    ! SendDn0u, RecvDn0u:
    ! messages update GlobalNoGhost [xe+1:xe+1,yb:ye]

    xbToBeUpdated = GlobalNoGhost%xe + 1
    xeToBeUpdated = xbToBeUpdated
    ybToBeUpdated = GlobalNoGhost%yb
    yeToBeUpdated = GlobalNoGhost%ye

    ! which neighbour nodes will send and receive

    call NodesToSendRecvMessages(myNum, Neigh, GlobalNoGhost, &
         xbToBeUpdated, xeToBeUpdated, ybToBeUpdated, yeToBeUpdated, &
         xbSend, xeSend, ybSend, yeSend, willSend, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv)

    ! build message set

    SendDn0u => CreateMessageSet(NameSendDn0u, TagDn0u, willSend, Neigh)
    RecvDn0u => CreateMessageSet(NameRecvDn0u, TagDn0u, willRecv, Neigh)

    ! get field

    vTabPtr => null()
    tmp_name='DN0U'
    call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

    ! build field sections to be sent and received

    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, SendDn0u)
    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, RecvDn0u)
  end subroutine CreateSendRecvDn0u




  subroutine CreateSendRecvDn0v(&
       gridId, nMachs, nNeigh, myNum, &
       Neigh, GlobalNoGhost, GlobalWithGhost, &
       SendDn0v, RecvDn0v)

    integer, intent(in) :: gridId
    integer, intent(in) :: nMachs
    integer, intent(in) :: nNeigh
    integer, intent(in) :: myNum

    type(NeighbourNodes), pointer :: Neigh          ! intent(in)
    type(DomainDecomp), pointer :: GlobalNoGhost    ! intent(in)
    type(DomainDecomp), pointer :: GlobalWithGhost  ! intent(in)
    type(MessageSet), pointer :: SendDn0v          ! intent(out)
    type(MessageSet), pointer :: RecvDn0v          ! intent(out)

    integer, parameter :: TagDn0v=30
    character(len=*), parameter :: NameSendDn0v="SendDn0v"
    character(len=*), parameter :: NameRecvDn0v="RecvDn0v"

    ! scratch arrays of size number of BRAMS processes
    ! containing global indices

    integer :: xbToBeUpdated(nMachs)
    integer :: xeToBeUpdated(nMachs)
    integer :: ybToBeUpdated(nMachs)
    integer :: yeToBeUpdated(nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(nNeigh)
    integer :: xeSend(nNeigh)
    integer :: ybSend(nNeigh)
    integer :: yeSend(nNeigh)
    integer :: xbRecv(nNeigh)
    integer :: xeRecv(nNeigh)
    integer :: ybRecv(nNeigh)
    integer :: yeRecv(nNeigh)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(nNeigh)
    logical :: willRecv(nNeigh)

    type(var_tables_r), pointer   :: vtabPtr => null()
    character(len=*), parameter :: h="**(CreateSendRecvDn0v)**"
    character(len=30) :: tmp_name

    ! SendDn0v, RecvDn0v:
    ! messages update GlobalNoGhost [xb:xe,ye+1:ye+1]

    xbToBeUpdated = GlobalNoGhost%xb
    xeToBeUpdated = GlobalNoGhost%xe
    ybToBeUpdated = GlobalNoGhost%ye+1
    yeToBeUpdated = ybToBeUpdated

    ! which neighbour nodes will send and receive

    call NodesToSendRecvMessages(myNum, Neigh, GlobalNoGhost, &
         xbToBeUpdated, xeToBeUpdated, ybToBeUpdated, yeToBeUpdated, &
         xbSend, xeSend, ybSend, yeSend, willSend, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv)

    ! build message set

    SendDn0v => CreateMessageSet(NameSendDn0v, TagDn0v, willSend, Neigh)
    RecvDn0v => CreateMessageSet(NameRecvDn0v, TagDn0v, willRecv, Neigh)

    ! get field

    vTabPtr => null()
    tmp_name='DN0V'
    call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

    ! build field sections to be sent and received

    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend, SendDn0v)
    call InsertFieldSectionAtMessageSet(&
         myNum, vTabPtr, Neigh, GlobalWithGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv, RecvDn0v)
  end subroutine CreateSendRecvDn0v






  subroutine CreateG3DSendRecv(&
       gridId, nMachs, nNeigh, myNum, &
       g3d_spread, g3d_smoothh, &
       GridSize, Neigh, GlobalNoGhost, GlobalWithGhost, &
       SendG3D, RecvG3D)

    integer, intent(in) :: gridId
    integer, intent(in) :: nMachs
    integer, intent(in) :: nNeigh
    integer, intent(in) :: myNum
    integer, intent(in) :: g3d_spread
    integer, intent(in) :: g3d_smoothh

    type(GridDims), pointer :: GridSize
    type(NeighbourNodes), pointer :: Neigh          ! intent(in)
    type(DomainDecomp), pointer :: GlobalNoGhost    ! intent(in)
    type(DomainDecomp), pointer :: GlobalWithGhost  ! intent(in)
    type(MessageSet), pointer :: SendG3D          ! intent(out)
    type(MessageSet), pointer :: RecvG3D          ! intent(out)

    integer :: vTabNbr
    integer, parameter :: TagG3D=31
    character(len=*), parameter :: NameSendG3D="SendG3D"
    character(len=*), parameter :: NameRecvG3D="RecvG3D"

    ! scratch arrays of size number of BRAMS processes
    ! containing global indices

    integer :: xbToBeUpdated(nMachs)
    integer :: xeToBeUpdated(nMachs)
    integer :: ybToBeUpdated(nMachs)
    integer :: yeToBeUpdated(nMachs)

    ! scratch arrays of size number of neighbour nodes
    ! containing global indices of regions for send and receive

    integer :: xbSend(nNeigh)
    integer :: xeSend(nNeigh)
    integer :: ybSend(nNeigh)
    integer :: yeSend(nNeigh)
    integer :: xbRecv(nNeigh)
    integer :: xeRecv(nNeigh)
    integer :: ybRecv(nNeigh)
    integer :: yeRecv(nNeigh)

    ! scratch arrays of size number of neighbour nodes
    ! containing which neighbour nodes will send of receive

    logical :: willSend(nNeigh)
    logical :: willRecv(nNeigh)

    type(var_tables_r), pointer   :: vtabPtr => null()
    character(len=*), parameter :: h="**(CreateG3DSendRecv)**"
    character(len=30) :: tmp_name

    ! SendG3D, RecvG3D:
    ! messages update entire GhostZone

    xbToBeUpdated = GlobalWithGhost%xb
    xeToBeUpdated = GlobalWithGhost%xe
    ybToBeUpdated = GlobalWithGhost%yb
    yeToBeUpdated = GlobalWithGhost%ye

    ! which neighbour nodes will send and receive

    call NodesToSendRecvMessages(myNum, Neigh, GlobalNoGhost, &
         xbToBeUpdated, xeToBeUpdated, ybToBeUpdated, yeToBeUpdated, &
         xbSend, xeSend, ybSend, yeSend, willSend, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv)

    ! include real domain boundaries

    call IncludeDomainBoundaries(Neigh, GridSize, GlobalNoGhost, &
         xbSend, xeSend, ybSend, yeSend, willSend)
    call IncludeDomainBoundaries(Neigh, GridSize, GlobalNoGhost, &
         xbRecv, xeRecv, ybRecv, yeRecv, willRecv)

    ! build message set

    SendG3D => CreateMessageSet(&
         NameSendG3D, TagG3D, &
         willSend, Neigh)
    RecvG3D => CreateMessageSet(&
         NameRecvG3D, TagG3D, &
         willRecv, Neigh)

    ! when g3d_spread is selected, send and receive fields TTENS and QVTTENS

    if (g3d_spread == 1) then

       vTabPtr => null()
       tmp_name='TTENS'
       call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

       call InsertFieldSectionAtMessageSet(&
            myNum, vTabPtr, Neigh, GlobalWithGhost, &
            xbSend, xeSend, ybSend, yeSend, willSend, &
            SendG3D)
       call InsertFieldSectionAtMessageSet(&
            myNum, vTabPtr, Neigh, GlobalWithGhost, &
            xbRecv, xeRecv, ybRecv, yeRecv, willRecv, &
            RecvG3D)

       vTabPtr => null()
       tmp_name='QVTTENS'
       call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

       call InsertFieldSectionAtMessageSet(&
            myNum, vTabPtr, Neigh, GlobalWithGhost, &
            xbSend, xeSend, ybSend, yeSend, willSend, &
            SendG3D)
       call InsertFieldSectionAtMessageSet(&
            myNum, vTabPtr, Neigh, GlobalWithGhost, &
            xbRecv, xeRecv, ybRecv, yeRecv, willRecv, &
            RecvG3D)
    end if

    ! when g3d_smoothh is selected, send and receive fields THSRC and RTSRC

    if (g3d_smoothh == 1) then

       vTabPtr => null()
       tmp_name='THSRC'
       call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

       call InsertFieldSectionAtMessageSet(&
            myNum, vTabPtr, Neigh, GlobalWithGhost, &
            xbSend, xeSend, ybSend, yeSend, willSend, &
            SendG3D)
       call InsertFieldSectionAtMessageSet(&
            myNum, vTabPtr, Neigh, GlobalWithGhost, &
            xbRecv, xeRecv, ybRecv, yeRecv, willRecv, &
            RecvG3D)

       vTabPtr => null()
       tmp_name='RTSRC'
       call GetVTabEntry(trim(tmp_name), gridId, vTabPtr)

       call InsertFieldSectionAtMessageSet(&
            myNum, vTabPtr, Neigh, GlobalWithGhost, &
            xbSend, xeSend, ybSend, yeSend, willSend, &
            SendG3D)
       call InsertFieldSectionAtMessageSet(&
            myNum, vTabPtr, Neigh, GlobalWithGhost, &
            xbRecv, xeRecv, ybRecv, yeRecv, willRecv, &
            RecvG3D)
    end if
  end subroutine CreateG3DSendRecv
end module ModMessagePassing
