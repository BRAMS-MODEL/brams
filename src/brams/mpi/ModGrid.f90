module ModGrid

  ! ModGrid: 

  use ModNamelistFile, only: &
       NamelistFile

  use ModParallelEnvironment, only: &
       ParallelEnvironment, &
       MsgDump

  use ModGridDims, only: &
       GridDims, &
       CreateGridDims, &
       DumpGridDims, &
       DestroyGridDims

  use ModDomainDecomp, only: &
       DomainDecomp, &
       CreateGlobalNoGhost, &
       CreateGlobalWithGhost, &
       CreateLocalInterior, &
       DumpDomainDecomp, &
       DestroyDomainDecomp

  use ModNeighbourNodes, only: &
       NeighbourNodes, &
       CreateNeighbourNodes, &
       DumpNeighbourNodes, &
       DestroyNeighbourNodes

  use ModMessageSet, only: &
       MessageSet, &
       DumpMessageSet

  use ModMessagePassing, only: &
       CreateAcousticMessagePassing, &
       DestroyAcousticMessagePassing, &
       CreateDn0MessagePassing, &
       DestroyDn0MessagePassing, &
       CreateG3DMessagePassing, &
       DestroyG3DMessagePassing, &
       CreateSelectedGhostZoneMessagePassing, &
       DestroySelectedGhostZoneMessagePassing, &
       CreateAllGhostZoneMessagePassing, &
       DestroyAllGhostZoneMessagePassing


  ! JP: temporariamente usa variaveis globais enquanto
  !     var_tables nao for inclusa no tipo Grid

  use var_tables, only: &
    num_var, &
    vtab_r
    
  use meteogramType, only: &
      PolygonContainer

  implicit none

  private
  public :: Grid
  public :: CreateGrid
  public :: InsertMessagePassingAtOneGrid
  public :: DestroyGrid
  public :: DumpGrid


  type Grid
     integer :: Id    ! grid number on Namelist
     type(NamelistFile), pointer :: Ramsin => null()
     type(ParallelEnvironment), pointer :: ParEnv => null()
     type(GridDims), pointer :: GridSize => null()
     type(DomainDecomp), pointer :: GlobalNoGhost => null()
     type(DomainDecomp), pointer :: GlobalWithGhost => null()
     type(DomainDecomp), pointer :: LocalInterior => null()
     type(NeighbourNodes), pointer :: Neigh => null()
     type(MessageSet), pointer :: AcouSendU
     type(MessageSet), pointer :: AcouRecvU
     type(MessageSet), pointer :: AcouSendV
     type(MessageSet), pointer :: AcouRecvV
     type(MessageSet), pointer :: AcouSendP
     type(MessageSet), pointer :: AcouRecvP
     type(MessageSet), pointer :: AcouSendUV
     type(MessageSet), pointer :: AcouRecvUV
     type(MessageSet), pointer :: AcouSendWP
     type(MessageSet), pointer :: AcouRecvWP
     type(MessageSet), pointer :: SendDn0u
     type(MessageSet), pointer :: RecvDn0u
     type(MessageSet), pointer :: SendDn0v
     type(MessageSet), pointer :: RecvDn0v
     type(MessageSet), pointer :: SendG3D
     type(MessageSet), pointer :: RecvG3D
     type(MessageSet), pointer :: SelectedGhostZoneSend
     type(MessageSet), pointer :: SelectedGhostZoneRecv
     type(MessageSet), pointer :: AllGhostZoneSend
     type(MessageSet), pointer :: AllGhostZoneRecv
     
     type(PolygonContainer), pointer :: meteoPolygons

  end type Grid


  logical, parameter :: dumpLocal=.false.

contains



  ! CreateGrid: create and fill variable of this type,
  !             extracting info from the Namelist File.



  subroutine CreateGrid(gridId, GhostZoneLength, &
       oneNamelistFile, oneParallelEnvironment, oneGrid)
    integer, intent(in) :: gridId
    integer, intent(in) :: GhostZoneLength
    type(NamelistFile), pointer :: oneNamelistFile
    type(ParallelEnvironment), pointer :: oneParallelEnvironment
    type(Grid), pointer :: oneGrid

    character(len=16) :: c0, c1
    character(len=*), parameter :: h="**(CreateGrid)**"

    ! correctness of input arguments

    if (.not. associated(oneNamelistFile)) then
       call fatal_error(h//" invoked with null oneNamelistFile")
    else if (.not. associated(oneParallelEnvironment)) then
       call fatal_error(h//" invoked with null oneParallelEnvironment")
    else if (associated(oneGrid)) then
       call fatal_error(h//" invoked with already associated oneGrid")
    end if

    ! create a variable of type grid and fill entries

    allocate(oneGrid)

    oneGrid%id = gridId
    oneGrid%Ramsin => oneNamelistFile
    oneGrid%ParEnv => oneParallelEnvironment
    call CreateGridDims(gridId, &
         oneNamelistFile, &
         oneGrid%GridSize)
    call CreateGlobalNoGhost(oneGrid%GridSize, &
         oneGrid%ParEnv, &
         oneGrid%GlobalNoGhost)
    call CreateGlobalWithGhost(oneGrid%GridSize, &
         oneGrid%ParEnv, &
         GhostZoneLength, &
         oneGrid%GlobalNoGhost, &
         oneGrid%GlobalWithGhost)
    call CreateLocalInterior(oneGrid%ParEnv, &
         oneGrid%GlobalWithGhost, &
         oneGrid%GlobalNoGhost, &
         oneGrid%LocalInterior)
    call CreateNeighbourNodes(oneGrid%ParEnv, &
         oneGrid%GlobalNoGhost, &
         oneGrid%GlobalWithGhost, &
         oneGrid%Neigh)
	 
    oneGrid%AcouSendU => null()
    oneGrid%AcouRecvU => null()
    oneGrid%AcouSendV => null()
    oneGrid%AcouRecvV => null()
    oneGrid%AcouSendP => null()
    oneGrid%AcouRecvP => null()
    oneGrid%AcouSendUV => null()
    oneGrid%AcouRecvUV => null()
    oneGrid%AcouSendWP => null()
    oneGrid%AcouRecvWP => null()
    oneGrid%SendDn0u => null()
    oneGrid%RecvDn0u => null()
    oneGrid%SendDn0v => null()
    oneGrid%RecvDn0v => null()
    oneGrid%SendG3D => null()
    oneGrid%RecvG3D => null()
    oneGrid%SelectedGhostZoneSend => null()
    oneGrid%SelectedGhostZoneRecv => null()
    oneGrid%AllGhostZoneSend => null()
    oneGrid%AllGhostZoneRecv => null()
    oneGrid%meteoPolygons => null()
    
  end subroutine CreateGrid





  subroutine InsertMessagePassingAtOneGrid(oneGrid)
    type(Grid), pointer :: oneGrid

    character(len=16) :: c0, c1
    character(len=*), parameter :: h="**(InsertMessagePassingAtOneGrid)**"

    if (.not. associated(oneGrid)) then
       call fatal_error(h//" invoked with null grid")
    end if

    call CreateAcousticMessagePassing(oneGrid%Id, &
         oneGrid%GridSize, oneGrid%ParEnv, oneGrid%Neigh, &
         oneGrid%GlobalNoGhost, &
         oneGrid%GlobalWithGhost, &
         oneGrid%AcouSendU, oneGrid%AcouRecvU, &
         oneGrid%AcouSendV, oneGrid%AcouRecvV,&
         oneGrid%AcouSendP, oneGrid%AcouRecvP, &
         oneGrid%AcouSendUV, oneGrid%AcouRecvUV, &
         oneGrid%AcouSendWP, oneGrid%AcouRecvWP)

    call CreateDn0MessagePassing(oneGrid%Id, &
         oneGrid%GridSize, oneGrid%ParEnv, oneGrid%Neigh, &
         oneGrid%GlobalNoGhost, oneGrid%GlobalWithGhost, &
         oneGrid%SendDn0u, oneGrid%RecvDn0u, &
         oneGrid%SendDn0v, oneGrid%RecvDn0v)

    call CreateG3DMessagePassing(oneGrid%Id, &
         oneGrid%GridSize, oneGrid%ParEnv, oneGrid%Neigh, &
         oneGrid%GlobalNoGhost, oneGrid%GlobalWithGhost, &
         oneGrid%Ramsin, &
         oneGrid%SendG3D, oneGrid%RecvG3D)

    ! temporariamente, num_var e vtab_r sao variaveis globais,
    ! enquanto nao inclusas no tipo Grid

    call CreateSelectedGhostZoneMessagePassing(&
       oneGrid%Id, num_var, vtab_r, &
       oneGrid%GridSize, oneGrid%ParEnv, oneGrid%Neigh, &
       oneGrid%GlobalNoGhost, oneGrid%GlobalWithGhost, &
       oneGrid%SelectedGhostZoneSend, oneGrid%SelectedGhostZoneRecv)

    call CreateAllGhostZoneMessagePassing(&
       oneGrid%Id, num_var, vtab_r, &
       oneGrid%GridSize, oneGrid%ParEnv, oneGrid%Neigh, &
       oneGrid%GlobalNoGhost, oneGrid%GlobalWithGhost, &
       oneGrid%AllGhostZoneSend, oneGrid%AllGhostZoneRecv)

    if (dumpLocal) then
       call MsgDump(h//" dumping oneGrid")
       call DumpGrid(OneGrid)
    end if
  end subroutine InsertMessagePassingAtOneGrid



  ! DestroyGrid: deallocate area of a variable of type grid



  subroutine DestroyGrid(oneGrid)
    type(Grid), pointer :: oneGrid

    if (associated(oneGrid)) then
       call DestroyGridDims(oneGrid%GridSize)
       call DestroyDomainDecomp(oneGrid%GlobalNoGhost)
       call DestroyDomainDecomp(oneGrid%GlobalWithGhost)
       call DestroyDomainDecomp(oneGrid%LocalInterior)
       call DestroyNeighbourNodes(oneGrid%Neigh)
       call DestroyAcousticMessagePassing(&
            oneGrid%AcouSendU, oneGrid%AcouRecvU, &
            oneGrid%AcouSendV, oneGrid%AcouRecvV, &
            oneGrid%AcouSendP, oneGrid%AcouRecvP, &
            oneGrid%AcouSendUV, oneGrid%AcouRecvUV, &
            oneGrid%AcouSendWP, oneGrid%AcouRecvWP)
       call DestroyDn0MessagePassing( &
            oneGrid%SendDn0u, oneGrid%RecvDn0u, &
            oneGrid%SendDn0v, oneGrid%RecvDn0v)
       call DestroyG3DMessagePassing( &
            oneGrid%SendG3D, oneGrid%RecvG3D)
       call DestroySelectedGhostZoneMessagePassing( &
            oneGrid%SelectedGhostZoneSend, &
            oneGrid%SelectedGhostZoneRecv)
       call DestroyAllGhostZoneMessagePassing( &
            oneGrid%AllGhostZoneSend, &
            oneGrid%AllGhostZoneRecv)

       deallocate(oneGrid)
    end if
    nullify(oneGrid)
  end subroutine DestroyGrid



  ! DumpGrid:


  subroutine DumpGrid(oneGrid)
    type(Grid), pointer :: oneGrid

    character(len=8) :: c0
    character(len=*), parameter :: h="**(DumpGrid)**"

    if (.not. associated(oneGrid)) then
       call fatal_error(h//" invoked with null oneGrid")
    else if (.not. associated(oneGrid%Ramsin)) then
       call fatal_error(h//" invoked with null oneGrid%Ramsin")
    else if (.not. associated(oneGrid%ParEnv)) then
       call fatal_error(h//" invoked with null oneGrid%ParEnv")
    end if

    write(c0,"(i8)") oneGrid%Id
    call MsgDump(h//" for grid "//trim(adjustl(c0)))

    call MsgDump(h//" dumping component GridSize")
    call DumpGridDims(oneGrid%GridSize)

    call MsgDump(h//" dumping domain decomposed components")
    call DumpDomainDecomp(oneGrid%GlobalNoGhost, "GlobalNoGhost")
    call DumpDomainDecomp(oneGrid%GlobalWithGhost, "GlobalWithGhost")
    call DumpDomainDecomp(oneGrid%LocalInterior, "LocalInterior")

    call MsgDump(h//" dumping neighborhood")
    call DumpNeighbourNodes(oneGrid%Neigh)

    call MsgDump(h//" dumping AcouSendU")
    call DumpMessageSet(oneGrid%AcouSendU)
    call MsgDump(h//" dumping AcouRecvU")
    call DumpMessageSet(oneGrid%AcouRecvU)
    call MsgDump(h//" dumping AcouSendV")
    call DumpMessageSet(oneGrid%AcouSendV)
    call MsgDump(h//" dumping AcouRecvV")
    call DumpMessageSet(oneGrid%AcouRecvV)
    call MsgDump(h//" dumping AcouSendP")
    call DumpMessageSet(oneGrid%AcouSendP)
    call MsgDump(h//" dumping AcouRecvP")
    call DumpMessageSet(oneGrid%AcouRecvP)
    call MsgDump(h//" dumping AcouSendUV")
    call DumpMessageSet(oneGrid%AcouSendUV)
    call MsgDump(h//" dumping AcouRecvUV")
    call DumpMessageSet(oneGrid%AcouRecvUV)
    call MsgDump(h//" dumping AcouSendWP")
    call DumpMessageSet(oneGrid%AcouSendWP)
    call MsgDump(h//" dumping AcouRecvWP")
    call DumpMessageSet(oneGrid%AcouRecvWP)
    call MsgDump(h//" dumping SendDn0u")
    call DumpMessageSet(oneGrid%SendDn0u)
    call MsgDump(h//" dumping RecvDn0u")
    call DumpMessageSet(oneGrid%RecvDn0u)
    call MsgDump(h//" dumping SendDn0v")
    call DumpMessageSet(oneGrid%SendDn0v)
    call MsgDump(h//" dumping RecvDn0v")
    call DumpMessageSet(oneGrid%RecvDn0v)
    call MsgDump(h//" dumping SendG3D")
    call DumpMessageSet(oneGrid%SendG3D)
    call MsgDump(h//" dumping RecvG3D")
    call DumpMessageSet(oneGrid%RecvG3D)
    call MsgDump(h//" dumping SelectedGhostZoneSend")
    call DumpMessageSet(oneGrid%SelectedGhostZoneSend)
    call MsgDump(h//" dumping SelectedGhostZoneRecv")
    call DumpMessageSet(oneGrid%SelectedGhostZoneRecv)
    call MsgDump(h//" dumping AllGhostZoneSend")
    call DumpMessageSet(oneGrid%AllGhostZoneSend)
    call MsgDump(h//" dumping AllGhostZoneRecv")
    call DumpMessageSet(oneGrid%AllGhostZoneRecv)
  end subroutine DumpGrid
end module ModGrid
