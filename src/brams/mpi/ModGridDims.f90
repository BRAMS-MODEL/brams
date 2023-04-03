module ModGridDims

  ! ModGridDims: stores grid dimensions of a single grid
  !          at a variable of type "GridDims"

  use ModNamelistFile, only: NamelistFile
  use ModParallelEnvironment, only: ParallelEnvironment, MsgDump
  implicit none
  private
  public :: GridDims
  public :: CreateGridDims
  public :: DumpGridDims
  public :: StringGridDims
  public :: DestroyGridDims

  type GridDims
     integer :: nnxp  ! # x points
     integer :: nnyp  ! # y points
     integer :: nnzp  ! # z points
  end type GridDims

contains



  ! CreateGridDims: create and fill variable of this type,
  !             extracting info from the Namelist File.



  subroutine CreateGridDims(gridId, oneNamelistFile, oneGridDims)
    integer, intent(in) :: gridId
    type(NamelistFile), pointer :: oneNamelistFile
    type(GridDims), pointer :: oneGridDims
    
    character(len=16) :: c0, c1
    character(len=512) :: str
    character(len=*), parameter :: h="**(CreateGridDims)**"
    logical, parameter :: dumpLocal=.false.

    ! correctness of input arguments

    if (.not. associated(oneNamelistFile)) then
       call fatal_error(h//" invoked with null oneNamelistFile")
    else if (oneNamelistFile%ngrids <= 0) then
       write(c0,"(i16)") oneNamelistFile%ngrids 
       call fatal_error(h//" ngrids on namelist file "//&
            "should be >= 1 but is"//trim(adjustl(c0)))
    else if (gridId > oneNamelistFile%ngrids) then
       write(c0,"(i16)") oneNamelistFile%ngrids 
       write(c1,"(i16)") gridId 
       call fatal_error(h//" invoked with gridId ("//&
            trim(adjustl(c1))//") > number of grids on Namelist File ("//&
            trim(adjustl(c0))//")")
    else if (gridId < 1) then
       write(c1,"(i16)") gridId 
       call fatal_error(h//" invoked with gridId ("//&
            trim(adjustl(c1))//") < 1 ")
    else if (associated(oneGridDims)) then
       call fatal_error(h//" invoked with oneGridDims already asociated")
    end if

    ! create a variable of type grid and fill entries

    allocate(oneGridDims)
    
    oneGridDims%nnxp = oneNamelistFile%nnxp(gridId)
    oneGridDims%nnyp = oneNamelistFile%nnyp(gridId)
    oneGridDims%nnzp = oneNamelistFile%nnzp(gridId)
    if (dumpLocal) then
       call StringGridDims(oneGridDims, str)
       call MsgDump(h//" ends producing grid "//trim(str))
    end if
  end subroutine CreateGridDims



  ! StringGridDims: creates string with info of a variable of type GridDims



  subroutine StringGridDims(oneGridDims, str)
    type(GridDims), pointer :: oneGridDims
    character(len=*), intent(out) :: str

    integer :: lenIn, lenOut
    character(len=8) :: c0, c1, c2
    character(len=*), parameter :: h="**(StringGridDims)**"

    if (.not. associated(oneGridDims)) then
       lenIn = len(str)
       lenOut = len(" empty grid")
       if (lenIn >= lenOut) then
          str = " empty grid"
       else
          write(c0,"(i8)") lenIn
          write(c1,"(i8)") lenOut
          call fatal_error(h//" input string has size "//trim(adjustl(c0))//&
               " but should be at least "//trim(adjustl(c1)))
       end if
    else
       write(c0,"(i8)") oneGridDims%nnxp
       write(c1,"(i8)") oneGridDims%nnyp
       write(c2,"(i8)") oneGridDims%nnzp
       lenOut = len(" with dimensions "//&
            "["//trim(adjustl(c0))//","//trim(adjustl(c1))//&
            ","//trim(adjustl(c2))//"]")
       lenIn = len(str)
       if (lenIn >= lenOut) then
          str = " with dimensions "//&
               "["//trim(adjustl(c0))//","//trim(adjustl(c1))//&
               ","//trim(adjustl(c2))//"]"
       else
          write(c0,"(i8)") lenIn
          write(c1,"(i8)") lenOut
          call fatal_error(h//" input string has size "//trim(adjustl(c0))//&
               " but should be at least "//trim(adjustl(c1)))
       end if
    end if
  end subroutine StringGridDims



  ! DumpGridDims: prints info of a variable of type grid on Dump file



  subroutine DumpGridDims(oneGridDims)
    type(GridDims), pointer :: oneGridDims

    character(len=512) :: str
    character(len=*), parameter :: h="**(DumpGridDims)**"

    if (.not. associated(oneGridDims)) then
       call MsgDump(h//" null GridDims")
    else
       call StringGridDims(oneGridDims, str)
       call MsgDump(h//trim(str))
    end if
  end subroutine DumpGridDims



  ! DestroyGridDims: deallocate area of a variable of type grid



  subroutine DestroyGridDims(oneGridDims)
    type(GridDims), pointer :: oneGridDims

    if (associated(oneGridDims)) then
       deallocate(oneGridDims)
    end if
    nullify(oneGridDims)
  end subroutine DestroyGridDims
end module ModGridDims
