module ModFieldSectionList


  ! ModFieldSectionList: List of field sections to be communicated
  !                      among processes in a single message


  use var_tables, only: &
       var_tables_r, &
       GetVTabSectionSize, &
       VerifyVTabEntry, &
       StringIndexing

  use ModParallelEnvironment, only: &
       MsgDump


  use ModDomainDecomp, only: &
       DomainDecomp

  implicit none

  private
  public :: FieldSection
  public :: CreateFieldSection
  public :: StringFieldSection
  public :: DumpFieldSection
  public :: DestroyFieldSection
  public :: FieldSectionList
  public :: CreateFieldSectionList
  public :: InsertAtFieldSectionList
  public :: RemoveFromFieldSectionList
  public :: DestroyFieldSectionList
  public :: DumpFieldSectionList
  public :: NextFieldSection


  ! FieldSection: one entry of a list of field sections to be sent/received 
  !               in a single message.
  !               It represents the section [xStart:xEnd,yStart:yEnd]
  !               of the field pointed by vTabPtr (2D, 3D or 4D).
  !               Number of points (reals) to be communicated is fieldSectionSize.


  type FieldSection
     type(var_tables_r), pointer :: vTabPtr ! field to communicate
     integer :: xStart, xEnd                ! local indices
     integer :: yStart, yEnd                ! local indices
     integer :: fieldSectionSize            ! # reals to communicate
     type(FieldSection), pointer :: next=>null()
     type(FieldSection), pointer :: previous=>null()
  end type FieldSection


  ! FieldSectionList: list of FieldSection


  type FieldSectionList
     type(FieldSection), pointer :: head=>null()
     type(FieldSection), pointer :: tail=>null()
  end type FieldSectionList

  logical, parameter :: dumpLocal=.false.

contains


  ! CreateFieldSection: returns pointer to a newly created 
  !                     FieldSection variable


  function CreateFieldSection(vTabPtr, &
       xStart, xEnd, yStart, yEnd, GlobalWithGhost) result(oneEntry)
    type(var_tables_r), pointer :: vTabPtr
    integer, intent(in) :: xStart
    integer, intent(in) :: xEnd
    integer, intent(in) :: yStart
    integer, intent(in) :: yEnd
    type(DomainDecomp), pointer :: GlobalWithGhost
    type(FieldSection), pointer :: oneEntry

    character(len=*), parameter :: h="**(CreateFieldSection)**"

    ! verify input arguments

    call VerifyVTabEntry(vTabPtr)

    ! allocate and fill entry

    allocate(oneEntry)
    oneEntry%vTabPtr => vTabPtr
    oneEntry%xStart = xStart
    oneEntry%xEnd = xEnd
    oneEntry%yStart = yStart
    oneEntry%yEnd = yEnd
    oneEntry%fieldSectionSize = GetVTabSectionSize(vTabPtr, &
         xStart, xEnd, yStart, yEnd)
    oneEntry%next=>null()
    oneEntry%previous=>null()
    if (dumpLocal) then
       call MsgDump(h//" with entries:")
       call DumpFieldSection(oneEntry)
    end if
  end function CreateFieldSection


  ! StringFieldSection: Returns a string with the fields of 
  !                     a variable of type FieldSection


  function StringFieldSection(oneEntry) result(res)
    type(FieldSection), pointer :: oneEntry
    character(len=256) :: res

    character(len=128) :: indexing
    character(len=8) :: c0
    character(len=*), parameter :: h="**(StringFieldSection)**"

    if (.not. associated(oneEntry)) then
       res = " null FieldSection"
    else if (.not. associated(oneEntry%vTabPtr)) then
       call fatal_error(h//&
            " FieldSection associated but vTabPtr not associated")
    else
       call StringIndexing(oneEntry%vTabPtr, &
            oneEntry%xStart, oneEntry%xEnd, &
            oneEntry%yStart, oneEntry%yEnd, indexing)
       write(c0,"(i8)") oneEntry%FieldSectionSize
       res = "field section (local indices) "//trim(oneEntry%vTabPtr%name)//&
            trim(indexing)//" of size "//trim(adjustl(c0))
    end if
  end function StringFieldSection


  ! DumpFieldSection: Dumps a variable of type FieldSection at 
  !                   this processor dump file


  subroutine DumpFieldSection(oneEntry)
    type(FieldSection), pointer :: oneEntry

    character(len=256) :: res
    character(len=*), parameter :: h="**(DumpFieldSection)**"

    res = StringFieldSection(oneEntry)
    call MsgDump(h//trim(adjustl(res)))
  end subroutine DumpFieldSection
            

  ! DestroyFieldSection: reclaims memory area and returns
  !                      null pointer


  subroutine DestroyFieldSection(oneEntry)
    type(FieldSection), pointer :: oneEntry

    character(len=128) :: name
    character(len=*), parameter :: h="**(DestroyFieldSection)**"

    if (associated(oneEntry)) then
       name = oneEntry%vTabPtr%name
       deallocate(oneEntry)
       if (dumpLocal) then
          call MsgDump(h//" named "//trim(adjustl(name)))
       end if
    end if
    nullify(oneEntry)
  end subroutine DestroyFieldSection


  ! CreateFieldSectionList: Creates empty list from a null pointer


  function CreateFieldSectionList() result (OneList)
    type (FieldSectionList), pointer :: OneList

    character(len=*), parameter :: h="**(CreateFieldSectionList)**"

    allocate(OneList)
    OneList%head=>null()
    OneList%tail=>null()
  end function CreateFieldSectionList


  ! InsertAtFieldSectionList: Given a FieldSectionList variable "list",
  !                           append FieldSection variable "node" to the
  !                           list


  subroutine InsertAtFieldSectionList(node, list)
    type(FieldSection), pointer :: node
    type(FieldSectionList), pointer :: list

    character(len=*), parameter :: h="**(InsertAtFieldSectionList)**"

    if (.not. associated(node)) then
       call fatal_error(h//" node not associated")
    else if (.not. associated(list)) then
       call fatal_error(h//" list not associated")
    end if

    if (.not. associated(list%head)) then

       ! if empty list, just point head and tail to the node

       list%head => node
       list%tail => node
    else if (.not. associated(list%tail)) then

       ! badly formed list if head is associated and tail is not associated

       call fatal_error(h//" bad list: head associated but null tail ")
    else

       ! append node to the end of the list

       list%tail%next => node
       node%previous => list%tail
       node%next => null()
       list%tail => node
    end if
    if (dumpLocal) then
       call MsgDump(h//" inserted the node "//&
            trim(adjustl(StringFieldSection(node))))
    end if
  end subroutine InsertAtFieldSectionList


  ! RemoveFromFieldSectionList: removes a node from the list,


  subroutine RemoveFromFieldSectionList(node, list)
    type(FieldSection), pointer :: node
    type(FieldSectionList), pointer :: list

    logical :: found
    type(FieldSection), pointer :: isThis
    character(len=*), parameter :: h="**(RemoveFromFieldSectionList)**"

    if (.not. associated(node)) then
       call fatal_error(h//" node not associated")
    else if (.not. associated(list)) then
       call fatal_error(h//" list not associated")
    end if

    found = .false.
    isThis => list%head
    do while (associated(isThis))
       found = associated(isThis, node)
       if (found) exit
       isThis => isThis%next
    end do
    if (.not. found) then
       call fatal_error(h//" node "//&
            trim(adjustl(StringFieldSection(node)))//&
            " is not in list")
    end if

    if (associated(list%head, node)) then
       list%head => node%next
       if (associated(list%head)) then
          list%head%previous => null()
       end if
    end if
    if (associated(list%tail, node)) then
       list%tail => node%previous
       if (associated(list%tail)) then
          list%tail%next => null()
       end if
    end if
    if (associated(node%previous)) then
       node%previous%next => node%next
    end if
    if (associated(node%next)) then
       node%next%previous => node%previous
    end if
    node%next => null()
    node%previous => null()
    if (dumpLocal) then
       call MsgDump(h//" removed the node "//&
            trim(adjustl(StringFieldSection(node))))
    end if
    call DestroyFieldSection(node)
  end subroutine RemoveFromFieldSectionList


  ! DestroyFieldSectionList: deallocates all nodes and the list


  subroutine DestroyFieldSectionList(list)
    type(FieldSectionList), pointer :: list

    type(FieldSection), pointer :: node
    character(len=*), parameter :: h="**(DestroyFieldSectionList)**"

    if (.not. associated(list)) then
       call MsgDump(h//" list not associated")
       return
    end if

!!$    call MsgDump(h//" list ")
!!$    call DumpFieldSectionList(list)
    do
       if (.not. associated(list%head) .and. &
            .not. associated(list%tail)) then
!!$          call MsgDump(h//" null tail, null head, will dealocate list")
          deallocate(list)
          list => null()
!!$          call MsgDump(h//" done deallocating list")
          exit
       else if (associated(list%head) .and. &
            associated(list%tail)) then
!!$          call MsgDump(h//" both tail and head associated, will dealocate node:")
          node => list%head
!!$          call DumpFieldSection(node)
          call RemoveFromFieldSectionList(node, list)
!!$          call MsgDump(h//" done deallocating node")
       else if (associated(list%head)) then
          call fatal_error(h//&
               " head associated but tail not associated")
       else if (associated(list%tail)) then
          call fatal_error(h//&
               " head not associated but tail associated")
       end if
    end do
  end subroutine DestroyFieldSectionList


  ! DumpFieldSectionList: Dumps a variable of type FieldSectionList at 
  !                       this processor dump file


  subroutine DumpFieldSectionList(list)
    type(FieldSectionList), pointer :: list

    integer :: cnt
    character(len=8) :: c0
    type(FieldSection), pointer :: node
    character(len=*), parameter :: h="**(DumpFieldSectionList)**"

    if (.not. associated(list)) then
       call MsgDump(h//" list is not associated")
    else

       ! list length

       cnt = 0
       node=>list%head
       do
          if (.not. associated(node)) exit
          cnt = cnt+1
          node => node%next
       end do
       select case (cnt)
       case (0)
          call MsgDump(h//" list with no entries")
       case (1)
          call MsgDump(h//" list has a single entry: "//&
               trim(adjustl(StringFieldSection(list%head))))
       case default
          write(c0,"(i8)") cnt
          call MsgDump(h//" list has the following "//&
               trim(adjustl(c0))//" entries:")
          
          ! dump each node
          
          node=>list%head
          do
             if (.not. associated(node)) then
                exit
             end if
             call DumpFieldSection(node)
             node => node%next
          end do
       end select
    end if
  end subroutine DumpFieldSectionList


  ! NextFieldSection: returns node following "node" at the list;
  !                   if input "node" is empty, returns list head;
  !                   if no more nodes in the list, returns null


  function NextFieldSection(node, list) result(next)
    type(FieldSection), pointer :: node
    type(FieldSectionList), pointer :: list
    type(FieldSection), pointer :: next

    character(len=*), parameter :: h="**(NextFieldSection)**"

    next => null()
    if (.not. associated(list)) then
       call fatal_error(h//" invoked with not associated list")
    else
       if (associated(list%head)) then
          if (.not. associated(node)) then
             next => list%head
          else
             next => node%next
          end if
       end if
    end if

    if (dumpLocal) then
       call MsgDump(h//" returns "//&
            trim(adjustl(StringFieldSection(next))))
    end if
  end function NextFieldSection
end module ModFieldSectionList
