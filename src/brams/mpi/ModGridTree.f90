module ModGridTree



  ! ModGridTree: The set of nesting grids froms a tree, with
  !              root at the outermost grid. A variable of type
  !              "GridTree" stores grid info on the entire tree.
  !              Each node is a variable of type "Grid".



  use ModNamelistFile, only: &
       NamelistFile

  use ModGrid, only: &
       Grid, &
       CreateGrid, &
       DestroyGrid, &
       DumpGrid

  use ModParallelEnvironment, only: &
       ParallelEnvironment, &
       MsgDump

  implicit none
  private
  public :: GridTree
  public :: CreateGridTree
  public :: GridTreeRoot
  public :: NextOnGridTree
  public :: FetchGrid
  public :: DumpGridTree
  public :: DestroyGridTree

  type GridTree
     type(Grid), pointer :: curr => null()           ! current grid
     type(GridTree), pointer :: ancestor => null()   ! one level up
     type(GridTree), pointer :: sibling => null()    ! at this level
     type(GridTree), pointer :: descendent => null() ! first node one level down
     type(GridTree), pointer :: PreOrder => null()   ! next in PreOrder
  end type GridTree

contains



  ! CreateGridTree: build all grids tree from the namelist file



  subroutine CreateGridTree(oneNamelistFile, oneParallelEnvironment, allGrids)
    type(namelistFile), pointer :: oneNamelistFile
    type(ParallelEnvironment), pointer :: oneParallelEnvironment
    type(GridTree), pointer :: allGrids

    ! AllGridNodes is an array of pointers to GridTree nodes,
    ! used as temporary area to build the entire grid.

    type GridTreePointer
       type(GridTree), pointer :: this
    end type GridTreePointer
    type(GridTreePointer), allocatable :: AllGridNodes(:)

    integer, parameter :: GhostZoneLength=1
    integer :: gridId
    type(GridTree), pointer :: oneGridNode, ancestor, brother, root
    type(Grid), pointer :: oneGrid
    character(len=16) :: c0, c1
    character(len=*), parameter :: h="**(CreateGridTree)**"

    ! input arguments consistency

    if (.not. associated(oneNamelistFile)) then
       call fatal_error(h//" invoked with null oneNamelistFile")
    else if (.not. associated(oneParallelEnvironment)) then
       call fatal_error(h//" invoked with null oneParallelEnvironment")
    else if (oneNamelistFile%ngrids <= 0) then
       write(c0,"(i16)") oneNamelistFile%ngrids 
       call fatal_error(h//" ngrids on namelist file "//&
            "should be >= 1 but is"//trim(adjustl(c0)))
    end if

    ! copy namelist file values into array of GridTrees pointers

    allocate(AllGridNodes(oneNamelistFile%ngrids))
    do gridId = 1, oneNamelistFile%ngrids

       ! create and fill a grid node

       allocate(AllGridNodes(gridId)%this)
       AllGridNodes(gridId)%this%curr => null()
       AllGridNodes(gridId)%this%ancestor => null()
       AllGridNodes(gridId)%this%sibling => null()
       AllGridNodes(gridId)%this%descendent => null()
       AllGridNodes(gridId)%this%PreOrder => null()
       call CreateGrid(gridId, GhostZoneLength, &
            oneNamelistFile, oneParallelEnvironment, &
            AllGridNodes(gridId)%this%curr)
    end do

    ! verify nesting correction

    if (oneNamelistFile%nxtnest(1) /= 0) then
       write(c0,"(i16)") oneNamelistFile%nxtnest(1)
       call fatal_error(h//" nxtnext(1) on namelist file "//&
            "should be 0 but is"//trim(adjustl(c0)))
    end if
    do gridId = 2, oneNamelistFile%ngrids
       if (oneNamelistFile%nxtnest(gridId) >= gridId) then
          write(c0,"(i16)") gridId
          write(c1,"(i16)") oneNamelistFile%nxtnest(gridId)
          call fatal_error(h//" nxtnext("//trim(adjustl(c0))//&
               ") on namelist file should be <= "//trim(adjustl(c0))//&
               " but is"//trim(adjustl(c1)))
       end if
    end do

    ! fill ancestors

    do gridId = 2, oneNamelistFile%ngrids
       oneGridNode => AllGridNodes(gridId)%this
       ancestor => AllGridNodes(oneNamelistFile%nxtnest(gridId))%this
       oneGridNode%ancestor => ancestor
       if (.not. associated(ancestor%descendent)) then
          ancestor%descendent => oneGridNode
       else
          brother => ancestor%descendent
          do 
             if (associated(brother%sibling)) then
                brother => brother%sibling
             else
                brother%sibling => oneGridNode
                exit
             end if
          end do
       end if
    end do

    ! build PreOrder

    root => AllGridNodes(1)%this
    call PreOrder(root, root%descendent)

    ! return coarser grid

    allGrids => root

    ! destroy scratch area

    deallocate(AllGridNodes)
  end subroutine CreateGridTree



  recursive subroutine PreOrder(previous, this)
    type(GridTree), pointer :: previous
    type(GridTree), pointer :: this
    character(len=*), parameter :: h="**(PreOrder)**"

    ! invariant: previous is always associated

    if (.not. associated(previous)) then
       call fatal_error(h//" invoked with null previous")
    end if
    if (associated(this)) then
       previous%PreOrder => this
       previous => this
       call PreOrder(previous, this%descendent)
       call PreOrder(previous, this%sibling)
    end if
  end subroutine PreOrder



  function GridTreeRoot(OneGridTree) result(this)
    type(GridTree), pointer :: OneGridTree
    type(GridTree), pointer :: this
    character(len=*), parameter :: h="**(GridTreeRoot)**"

    if (.not. associated(OneGridTree)) then
       call fatal_error(h//" invoked with null tree")
    end if
    this => OneGridTree
    do
       if (.not. associated(this%ancestor)) then
          return
       else
          this => this%ancestor
       end if
    end do
  end function GridTreeRoot
       


  function NextOnGridTree(OneGridTree) result(this)
    type(GridTree), pointer :: OneGridTree
    type(GridTree), pointer :: this
    if (associated(OneGridTree)) then
       this => OneGridTree%PreOrder
    else
       nullify(this)
    end if
  end function NextOnGridTree


  ! DumpGridTree: dumps the tree in pre-order


  subroutine DumpGridTree(OneGridTree)
    type(GridTree), pointer :: OneGridTree
    type(GridTree), pointer :: this

    character(len=*), parameter :: h="**(DumpGridTree)**"

    if (.not. associated(OneGridTree)) then
       call MsgDump(h//" empty tree")
       return
    else
       this => GridTreeRoot(OneGridTree)
       do
          if (.not. associated(this)) then
             exit
          else
             call DumpOneGridTreeNode(this)
             this => NextOnGridTree(this)
          end if
       end do
    end if
  end subroutine DumpGridTree


  ! DumpOneGridTreeNode: dumps a single grid



  subroutine DumpOneGridTreeNode(this)
    type(GridTree), pointer :: this

    character(len=8) :: c0, c1, c2, c3
    character(len=*), parameter :: h="**(DumpOneGridTreeNode)**"

    if (.not. associated(this)) then
       return
    else if (.not. associated(this%curr)) then
       call fatal_error(h//" inconsistent grid tree")
    end if

    ! grid id

    write(c0,"(i8)") this%curr%Id
    call MsgDump(h//" grid "//trim(adjustl(c0)), .true.)

    ! ancestors

    if (associated(this%ancestor)) then
       write(c0,"(i8)") this%ancestor%curr%Id
       call MsgDump(" is nested on grid "//trim(adjustl(c0)),.true.)
    else if (this%curr%id == 1) then         
       call MsgDump(" is outermost grid ",.true.)
    else
       call fatal_error(h//" inconsistent grid tree")
    end if

    ! sibling

    if (associated(this%sibling)) then
       write(c0,"(i8)") this%sibling%curr%Id
       call MsgDump(", has grid "//trim(adjustl(c0))//&
            " as first sibling",.true.)
    else
       call MsgDump(", has no sibling",.true.)
    end if

    ! descendent

    if (associated(this%descendent)) then
       write(c0,"(i8)") this%descendent%curr%Id
       call MsgDump(" and has grid "//trim(adjustl(c0))//&
            " as first descendent")
    else
       call MsgDump(" and has no descendent")
    end if
    call DumpGrid(this%curr)
  end subroutine DumpOneGridTreeNode


  ! DestroyGridTree: destroys the tree in depth first order


  recursive subroutine DestroyGridTree(OneGridTree)
    type(GridTree), pointer :: OneGridTree

    if (.not. associated(OneGridTree)) then
       return
    else
       call DestroyGridTree(OneGridTree%descendent)
       call DestroyGridTree(OneGridTree%sibling)
       call DestroyOneGridTreeNode(OneGridTree)
    end if
  end subroutine DestroyGridTree


  ! DestroyOneGridTreeNode: destroys one node



  subroutine DestroyOneGridTreeNode(this)
    type(GridTree), pointer :: this

    character(len=*), parameter :: h="**(DestroyOneGridTreeNode)**"

    if (.not. associated(this)) then
       return
    else if (.not. associated(this%curr)) then
       call fatal_error(h//" not conforming grid tree")
    end if

    call DestroyGrid(this%curr)
    deallocate(this)
    nullify(this)
  end subroutine DestroyOneGridTreeNode



  ! FetchGrid: returns pointer to grid number "nbr"



  function FetchGrid(AllGrids, nbr) result(OneGrid)
    type(GridTree), pointer :: AllGrids
    integer, intent(in) :: nbr
    type(Grid), pointer :: OneGrid

    type(GridTree), pointer :: OneGridTreeNode
    logical :: found
    character(len=8) :: c0
    character(len=*), parameter :: h="**(FetchGrid)**"

    if (.not. associated(AllGrids)) then
       call fatal_error(h//" null AllGrids")
    end if
    OneGridTreeNode => GridTreeRoot(AllGrids)
    found = .false.
    do while (associated(OneGridTreeNode))
       found = OneGridTreeNode%curr%Id == nbr
       if (found) then
          exit
       else
          OneGridTreeNode => NextOnGridTree(OneGridTreeNode)
       end if
    end do
    if (.not. found) then
       write(c0,"(i8)") nbr
       call fatal_error(h//" grid number "//trim(adjustl(c0))//" not found")
    end if
    OneGrid => OneGridTreeNode%curr
  end function FetchGrid
end module ModGridTree
