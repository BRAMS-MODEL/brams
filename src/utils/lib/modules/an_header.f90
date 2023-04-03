!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

module an_header

  implicit none

  include "i8.h"
  include "files.h"

  ! head_table: single field on i/o file descriptor data structure

  type head_table
     character(len=16) :: string     ! field name
     integer           :: idim_type  ! field dimensionality (coded)
     integer           :: ngrid      ! grid number
     integer(kind=i8)  :: npointer   ! on ASCII coded files, field starting position
     integer(kind=i8)  :: nvalues    ! field size
  end type head_table

  ! IOHeadTable: i/o file descriptor data structure, with all fields

  type IOHeadTable
     integer                   :: maxSize    ! maximum number of field entries
     integer                   :: lastUsed   ! last used field entry
     type(head_table), pointer :: f(:)       ! field entries
  end type IOHeadTable

  ! IOFileDS: i/o file data structure

  type IOFileDS
     logical                     :: enable      ! if this i/o is enable
     character(len=4)            :: fId         ! file ID (one of "HIST", "INST", "LITE", "MEAN")
     character(len=f_name_length)          :: fHeadName   ! header file name
     logical                     :: collapsed   ! one file collapsing all grids (true) or one file per grid (false)
     logical                     :: binary      ! binary or ascii coded data files
     character(len=f_name_length), pointer :: fName(:)    ! file names (indexed by grid number if not collapsed)
     integer(kind=i8)            :: fPos        ! last position on current file
     integer                     :: unit        ! Fortran i/o unit or C FILE handler for current opened file
     integer                     :: iclobber    ! 
     type(IOHeadTable)           :: ht          ! all fields on all data files
  end type IOFileDS

  type (head_table), allocatable,save :: anal_table(:)  !**(JP)** old
  integer, save:: nvbtab                                !**(JP)** old

  interface FieldWriteASCIIStoreInfo
     module procedure ArrayWriteASCIIStoreInfo, PointerWriteASCIIStoreInfo
  end interface

  interface FieldWriteBinStoreInfo
     module procedure ArrayWriteBinStoreInfo, PointerWriteBinStoreInfo
  end interface

  interface FieldWriteStoreInfo
     module procedure ArrayWriteStoreInfo, PointerWriteStoreInfo
  end interface

  logical, parameter :: dumpLocal=.false.
contains



  subroutine CreateIOHeadTable (maxSize, oneIOHeadTable)
    integer,       intent(in)    :: maxSize
    type(IOHeadTable), intent(inout) :: oneIOHeadTable

    integer :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(CreateIOHeadTable)**"

    ! protection: procedure arguments

!--(DMK-CCATT-INI)-------------------------------------------------------------------
!    if (associated(oneIOHeadTable%f)) then
!       deallocate(oneIOHeadTable%f, stat=ierr)
!       if (ierr /= 0) then
!          write(c0,"(i8)") ierr
!          call fatal_error(h//" deallocate f fails with stat="//trim(adjustl(c0)))
!       end if
!       nullify(oneIOHeadTable%f)
!    end if
!--(DMK-CCATT-FIM)-------------------------------------------------------------------

    if (maxSize < 1) then
       write(c0,"(i8)") maxSize
       call fatal_error(h//" invoked with non-positive maxSize="//trim(adjustl(c0)))
    end if

    ! create and initialize variable

    allocate(oneIOHeadTable%f(maxSize), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       write(c1,"(i8)") maxSize
       call fatal_error(h//" allocate f("//trim(adjustl(c1))//&
            ") fails with stat="//trim(adjustl(c0)))
    end if
    oneIOHeadTable%maxSize = maxSize
    oneIOHeadTable%lastUsed = 0

    if (dumpLocal) then
       write(c0,"(i8)") maxSize
       write(*,"(a)") h//" with maxSize="//trim(adjustl(c0))
    end if
  end subroutine CreateIOHeadTable




  subroutine DestroyIOHeadTable(oneIOHeadTable)
    type(IOHeadTable), intent(inout) :: oneIOHeadTable

    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(DestroyIOHeadTable)**"

    ! deallocate fields

    if (associated(oneIOHeadTable%f)) then
       deallocate(oneIOHeadTable%f, stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") ierr
          call fatal_error(h//" deallocate f fails with stat="//trim(adjustl(c0)))
       end if
       nullify(oneIOHeadTable%f)
    end if
    oneIOHeadTable%maxSize=0
    oneIOHeadTable%lastUsed=0
    if (dumpLocal) then
       write(*,"(a)") h//" done"
    end if
  end subroutine DestroyIOHeadTable





  subroutine InsertFieldOnIOHeadTable(oneIOHeadTable, &
       string, idim_type, ngrid, npointer, nvalues)
    type(IOHeadTable), intent(inout) :: oneIOHeadTable
    character(len=*), intent(in) :: string
    integer,          intent(in) :: idim_type
    integer,          intent(in) :: ngrid
    integer(kind=i8), intent(in) :: npointer
    integer(kind=i8), intent(in) :: nvalues

    integer :: ierr
    character(len=8) :: c0, c1, c2, c3
    character(len=*), parameter :: h="**(InsertFieldOnIOHeadTable)**"

    if (.not. associated(oneIOHeadTable%f)) then
       call fatal_error(h//" invoked before creating IOHeadTable variable")
    end if

    oneIOHeadTable%lastUsed = oneIOHeadTable%lastUsed + 1
    if (oneIOHeadTable%lastUsed > oneIOHeadTable%maxSize) then
       write(c0,"(i8)") oneIOHeadTable%maxSize
       call fatal_error(h//" exceeded maximum table size of "//&
            " ("//trim(adjustl(c0))//")")
    end if

    oneIOHeadTable%f(oneIOHeadTable%lastUsed)%string = string
    oneIOHeadTable%f(oneIOHeadTable%lastUsed)%idim_type = idim_type
    oneIOHeadTable%f(oneIOHeadTable%lastUsed)%ngrid = ngrid
    oneIOHeadTable%f(oneIOHeadTable%lastUsed)%npointer = npointer
    oneIOHeadTable%f(oneIOHeadTable%lastUsed)%nvalues = nvalues

    if (dumpLocal) then
       write(c0,"(i8)") oneIOHeadTable%lastUsed
       write(c1,"(i8)") idim_type
       write(c2,"(i8)") ngrid
       write(c3,"(i8)") nvalues
       write(*,"(a)") h//" field "//trim(string)//" at pos "//trim(adjustl(c0))//&
            "; idim_type="//trim(adjustl(c1))//"; grid="//trim(adjustl(c2))//&
            "; nvalues="//trim(adjustl(c3))
    end if
  end subroutine InsertFieldOnIOHeadTable




  subroutine CreateDisabledIOFileDS(cfId, oneIOFileDS)
    character(len=4), intent(in ) :: cfId          ! file ID (one of "HIST", "INST", "LITE", "MEAN")
    type(IOFileDS),   intent(out) :: oneIOFileDS

    character(len=*), parameter :: h="**(CreateDisabledIOFileDS)**"
    oneIOFileDS%fId = cfId
    oneIOFileDS%enable = .false.
    if (dumpLocal) then
       write(*,"(a)") h//" for fId "//trim(cfId)
    end if
  end subroutine CreateDisabledIOFileDS





  subroutine CreateEnabledIOFileDS (cfId, cPrefix, cType, iclobber, &
       tinc, iyr, imn, idy, itm, collapsed, binary, ngrids, maxNFields, &
       oneIOFileDS)
    character(len=4), intent(in ) :: cfId          ! file ID (one of "HIST", "INST", "LITE", "MEAN")
    character(len=f_name_length), intent(in ) :: cPrefix       ! file name prefix
    character,        intent(in ) :: cType
    integer,          intent(in ) :: iclobber
    real,             intent(in ) :: tinc
    integer,          intent(in ) :: iyr
    integer,          intent(in ) :: imn
    integer,          intent(in ) :: idy
    integer,          intent(in ) :: itm
    logical,          intent(in ) :: collapsed
    logical,          intent(in ) :: binary
    integer,          intent(in ) :: ngrids
    integer,          intent(in ) :: maxNFields
    type(IOFileDS),   intent(out) :: oneIOFileDS

    integer :: ng
    character(len=2) :: cGrid
    character(len=3) :: cFmt
    integer :: ierr
    character(len=8) :: c0, c1
    character(len=10) :: c2
    character :: l0, l1
    character(len=*), parameter :: h="**(CreateEnabledIOFileDS)**"

    ! register input values

    oneIOFileDS%fId = cfId
    oneIOFileDS%iclobber = iclobber
    oneIOFileDS%collapsed = collapsed
    oneIOFileDS%binary = binary
    oneIOFileDS%enable = .true.

    if (dumpLocal) then
       write(c0,"(i8)") iclobber
       write(l0,"(l1)") collapsed
       write(l1,"(l1)") binary
       write(c1,"(f8.3)") tinc
       write(c2,"(i4.4,3i2.2)") iyr, imn, idy, itm
       write(*,"(a)") h//" for "//cfId//", iclobber="//trim(adjustl(c0))//&
            ", collapsed="//l0//", binary="//l1//", tinc="//c1//&
            ", iyr,imn,idy,itm="//c2
       call flush(6)
    end if
       
    ! build header file name

    call makefnam (oneIOFileDS%fHeadName, cPrefix, & !trim(cPrefix), &
         tinc, iyr, imn, idy, itm*100, cType, "head", "txt")

    if (dumpLocal) then
       write(*,"(a)") h//" named header file "//trim(oneIOFileDS%fHeadName)
       call flush(6)
    end if

    ! build data file names sufix: "bin" for binary files, "vfm" for ascii files

    if (oneIOFileDS%binary) then
       cFmt = "bin"
    else
       cFmt = "vfm"
    end if

    ! allocate and build data file names

    if (oneIOFileDS%collapsed) then

       ! collapse all grids into a single file

       allocate(oneIOFileDS%fName(1), stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") ierr
          call fatal_error(h//" allocate fName(1) failed with stat="//trim(adjustl(c0)))
       end if

       call makefnam (oneIOFileDS%fName(1), cPrefix, tinc, iyr, imn, idy, itm*100, &
            cType, "$", cFmt)

       if (dumpLocal) then
          write(*,"(a)") h//" named data file "//trim(oneIOFileDS%fName(1))
       end if

    else

       ! one file per grid

       allocate(oneIOFileDS%fName(ngrids), stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") ierr
          write(c1,"(i8)") ngrids
          call fatal_error(h//" allocate fName("//trim(adjustl(c1))//&
               ") failed with stat="//trim(adjustl(c0)))
       end if

       do ng = 1, ngrids
          write(cGrid, "(a1,i1)") "g",ng
          call makefnam (oneIOFileDS%fName(ng), cPrefix, tinc, iyr, imn, idy, itm*100, &
               cType, cGrid, cFmt)
          if (dumpLocal) then
             write(*,"(a)") h//" named data file "//trim(oneIOFileDS%fName(ng))
          end if
       end do

    end if

    ! create header table

    call CreateIOHeadTable (maxNFields, oneIOFileDS%ht)

    ! mark data file not opened

    oneIOFileDS%unit=-1
  end subroutine CreateEnabledIOFileDS




  logical function IsIOFileDSEnabled(oneIOFileDS)
    type(IOFileDS), intent(in) :: oneIOFileDS
    IsIOFileDSEnabled = oneIOFileDS%enable
  end function IsIOFileDSEnabled




  subroutine OpenDataFile (oneIOFileDS, ngrid)
    type(IOFileDS), intent(inout) :: oneIOFileDS
    integer,        intent(in   ) :: ngrid

    integer, parameter :: unitLow=10
    integer, parameter :: unitHigh=99
    integer :: iunit
    integer :: index
    logical :: op
    character(len=8) :: c0
    character(len=*), parameter :: h="**(OpenDataFile)**"

    ! return if disabled 

    if (.not. oneIOFileDS%enable) return

    if (oneIOFileDS%unit == -1) then

       ! open one data file for all grids (if collapsed) 
       ! or one data file for this grid (if not collapsed)

       if (oneIOFileDS%collapsed) then
          index = 1
       else
          index = ngrid
       end if

       ! open data file: Fortran if binary, C if ascii

       if (oneIOFileDS%binary) then

          ! find available Fortran unit

          do iunit = unitLow, unitHigh
             inquire (unit=iunit, opened=op)
             if (.not. op) exit
          end do
          if (iunit <= unitHigh) then
             oneIOFileDS%unit = iunit
          else
             call fatal_error(h//" Fortran i/o units exausted")
          end if

          ! append to existing file (if collapsed and not first grid)
          ! or create a new one; reset file position if new file

          if (oneIOFileDS%collapsed .and. ngrid /= 1) then
             call rams_f_open_u (oneIOFileDS%unit, trim(oneIOFileDS%fName(index)), &
                  "UNFORMATTED", "OLD", "WRITE", "APPEND", &
                  oneIOFileDS%iclobber)
             if (dumpLocal) then
                write(c0,"(i8)") oneIOFileDS%unit
                write(*,"(a)") h//" open to append old Fortran file "//&
                     trim(oneIOFileDS%fName(index))//&
                     " at unit "//trim(adjustl(c0))
             end if
          else
             call rams_f_open_u (oneIOFileDS%unit, trim(oneIOFileDS%fName(index)), &
                  "UNFORMATTED", "REPLACE", "WRITE", "ASIS", &
                  oneIOFileDS%iclobber)
             oneIOFileDS%fPos = 0_i8
             if (dumpLocal) then
                write(c0,"(i8)") oneIOFileDS%unit
                write(*,"(a)") h//" open new Fortran file "//trim(oneIOFileDS%fName(index))//&
                     " at unit "//trim(adjustl(c0))
             end if
          end if
       else

          ! open C file; store returned "C unit" handle;
          ! reset file position

          call rams_c_open_u(trim(oneIOFileDS%fName(index))//char(0), "w"//char(0), &
               oneIOFileDS%unit)
          oneIOFileDS%fPos = 0_i8
          if (dumpLocal) then
             write(c0,"(i8)") oneIOFileDS%unit
             write(*,"(a)") h//" open new C file "//trim(oneIOFileDS%fName(index))//&
                  " at unit "//trim(adjustl(c0))
          end if
       end if

       ! file opened OK?

       if (oneIOFileDS%unit == -1) then
          call fatal_error(h//" open file for "//trim(oneIOFileDS%fId)//" fails")
       end if

    else

       ! previous opened file was not closed

       call fatal_error(h//" previously opened file for "//oneIOFileDS%fId//" not closed")
    end if
  end subroutine OpenDataFile




  subroutine CloseDataFile (oneIOFileDS)
    type(IOFileDS), intent(inout) :: oneIOFileDS

    character(len=*), parameter :: h="**(CloseDataFile)**"


    ! return if disabled 

    if (.not. oneIOFileDS%enable) return

    if (oneIOFileDS%unit /= -1) then

       ! close data file: Fortran if binary, C if ascii

       if (oneIOFileDS%binary) then
          close(oneIOFileDS%unit)
          if (dumpLocal) then
             write (*,"(a)") h//" close Fortran file for "//oneIOFileDS%fId
          end if
       else
          call rams_c_close_u(oneIOFileDS%unit)
          if (dumpLocal) then
             write (*,"(a)") h//" close C file for "//oneIOFileDS%fId
          end if
       end if
       oneIOFileDS%unit=-1

    else
       call fatal_error(h//" trying to close an unopened file for "//oneIOFileDS%fId)
    end if
  end subroutine CloseDataFile





  subroutine ArrayWriteASCIIStoreInfo (field, fieldSize, &
       varn, idim_type, ngr, oneIOFileDS)

    integer,           intent(in)    :: fieldSize         ! size of field to write
    real,              intent(in)    :: field(fieldSize)  ! field to write
    character(len=16), intent(in)    :: varn              ! field name
    integer,           intent(in)    :: idim_type         ! field dimensionality (coded)
    integer,           intent(in)    :: ngr               ! grid number
    type(IOFileDS),    intent(inout) :: oneIOFileDS

    ! local variables

    integer          :: ierr
    integer(kind=i8) :: fieldSize_i8
    character(len=8) :: c0
    character(len=8) :: c1
    character(len=*), parameter :: h="**(ArrayWriteASCIIStoreInfo)**"

    integer :: tam
!!$    character(len=128) :: line

    real, allocatable :: scr(:)
    real, allocatable :: cscr(:)


    ! return if disabled 

    if (.not. oneIOFileDS%enable) return

    ! scratch areas

    allocate(scr(fieldSize), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") fieldSize
       write(c1,"(i8)") ierr
       call fatal_error(h//" allocate scr("// &
            trim(adjustl(c0))//") failed with stat="// &
            trim(adjustl(c1)))
    end if

    allocate(cscr(fieldSize), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") fieldSize
       write(c1,"(i8)") ierr
       call fatal_error(h//" allocate cscr("// &
            trim(adjustl(c0))//") failed with stat="// &
            trim(adjustl(c1)))
    end if

    ! store info on field to be written info at next available table entry

    if (dumpLocal) then
       write(*,"(a)") h//" will insert field "//varn//" at "//oneIOFileDS%fId
    end if
    fieldSize_i8 = int(fieldSize,i8)
    call InsertFieldOnIOHeadTable(oneIOFileDS%ht, &
         varn, idim_type, ngr, oneIOFileDS%fPos, fieldSize_i8)

    ! code field and dump coded field on file

    if (dumpLocal) then
       write(c0,"(i8)") oneIOFileDS%unit
       write(*,"(a)") h//" will write field "//varn//" at unit "//trim(adjustl(c0))
    end if
    call vforecr_u(oneIOFileDS%unit, field, fieldSize_i8, 18,  scr, cscr, &
         'LIN', oneIOFileDS%fPos)

    ! deallocate scratch area

    deallocate(scr, stat=ierr)
    if (ierr /= 0) then
       write(c1,"(i8)") ierr
       call fatal_error(h//" deallocate scr failed with stat="//trim(adjustl(c1)))
    end if

    deallocate(cscr, stat=ierr)
    if (ierr /= 0) then
       write(c1,"(i8)") ierr
       call fatal_error(h//" deallocate cscr failed with stat="//trim(adjustl(c1)))
    end if

  end subroutine ArrayWriteASCIIStoreInfo




  subroutine PointerWriteASCIIStoreInfo (field, fieldSize, &
       varn, idim_type, ngr, oneIOFileDS)

    integer,           intent(in)    :: fieldSize         ! size of field to write
    real,              pointer       :: field             ! field to write (points to first position of array)
    character(len=16), intent(in)    :: varn              ! field name
    integer,           intent(in)    :: idim_type         ! field dimensionality (coded)
    integer,           intent(in)    :: ngr               ! grid number
    type(IOFileDS),    intent(inout) :: oneIOFileDS

    ! local variables

    integer          :: ierr
    integer(kind=i8) :: fieldSize_i8
    character(len=8) :: c0
    character(len=8) :: c1
    character(len=*), parameter :: h="**(PointerWriteASCIIStoreInfo)**"

    integer :: tam
!!$    character(len=128) :: line

    real, allocatable :: scr(:)
    real, allocatable :: cscr(:)

    ! return if disabled 

    if (.not. oneIOFileDS%enable) return

    ! scratch areas

    allocate(scr(fieldSize), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") fieldSize
       write(c1,"(i8)") ierr
       call fatal_error(h//" allocate scr("// &
            trim(adjustl(c0))//") failed with stat="// &
            trim(adjustl(c1)))
    end if

    allocate(cscr(fieldSize), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") fieldSize
       write(c1,"(i8)") ierr
       call fatal_error(h//" allocate cscr("// &
            trim(adjustl(c0))//") failed with stat="// &
            trim(adjustl(c1)))
    end if

    ! store info on field to be written info at next available table entry

    if (dumpLocal) then
       write(*,"(a)") h//" will insert field "//varn//" at "//oneIOFileDS%fId
    end if
    fieldSize_i8 = int(fieldSize,i8)
    call InsertFieldOnIOHeadTable(oneIOFileDS%ht, &
         varn, idim_type, ngr, oneIOFileDS%fPos, fieldSize_i8)

    ! code field and dump coded field on file

    if (dumpLocal) then
       write(c0,"(i8)") oneIOFileDS%unit
       write(*,"(a)") h//" will write field "//varn//" at unit "//trim(adjustl(c0))
    end if
    call vforecr_u(oneIOFileDS%unit, field, fieldSize_i8, 18,  scr, cscr, &
         'LIN', oneIOFileDS%fPos)

    ! deallocate scratch area

    deallocate(scr, stat=ierr)
    if (ierr /= 0) then
       write(c1,"(i8)") ierr
       call fatal_error(h//" deallocate scr failed with stat="//trim(adjustl(c1)))
    end if

    deallocate(cscr, stat=ierr)
    if (ierr /= 0) then
       write(c1,"(i8)") ierr
       call fatal_error(h//" deallocate cscr failed with stat="//trim(adjustl(c1)))
    end if

  end subroutine PointerWriteASCIIStoreInfo







  subroutine ArrayWriteBinStoreInfo (field, fieldSize, &
       varn, idim_type, ngr, oneIOFileDS)

    integer,           intent(in)    :: fieldSize         ! size of field to write
    real,              intent(in)    :: field(fieldSize)  ! field to write
    character(len=16), intent(in)    :: varn              ! field name
    integer,           intent(in)    :: idim_type         ! field dimensionality (coded)
    integer,           intent(in)    :: ngr               ! grid number
    type(IOFileDS),    intent(inout) :: oneIOFileDS

    ! local variables

    integer          :: ierr
    integer(kind=i8) :: fieldSize_i8
    character(len=8) :: c0
    character(len=8) :: c1
    character(len=*), parameter :: h="**(ArrayWriteBinStoreInfo)**"

    integer :: tam
!!$    character(len=128) :: line

    ! return if disabled 

    if (.not. oneIOFileDS%enable) return

    ! store info on field to be written info at next available table entry

    if (dumpLocal) then
       write(*,"(a)") h//" will insert field "//varn//" at "//oneIOFileDS%fId
    end if
    fieldSize_i8 = int(fieldSize,i8)
    call InsertFieldOnIOHeadTable(oneIOFileDS%ht, &
         varn, idim_type, ngr, oneIOFileDS%fPos, fieldSize_i8)

    ! dump binary field on file

    if (dumpLocal) then
       write(c0,"(i8)") oneIOFileDS%unit
       write(*,"(a)") h//" will write field "//varn//" at unit "//trim(adjustl(c0))
    end if
    call writebin(oneIOFileDS%unit, field, fieldSize_i8)
    oneIOFileDS%fPos = oneIOFileDS%fPos + fieldSize_i8

  end subroutine ArrayWriteBinStoreInfo







  subroutine PointerWriteBinStoreInfo (field, fieldSize, &
       varn, idim_type, ngr, oneIOFileDS)

    integer,           intent(in)    :: fieldSize         ! size of field to write
    real,              pointer       :: field             ! points to field first position
    character(len=16), intent(in)    :: varn              ! field name
    integer,           intent(in)    :: idim_type         ! field dimensionality (coded)
    integer,           intent(in)    :: ngr               ! grid number
    type(IOFileDS),    intent(inout) :: oneIOFileDS

    ! local variables

    integer          :: ierr
    integer(kind=i8) :: fieldSize_i8
    character(len=8) :: c0
    character(len=8) :: c1
    character(len=*), parameter :: h="**(PointerWriteBinStoreInfo)**"

    integer :: tam
!!$    character(len=128) :: line

    ! return if disabled 

    if (.not. oneIOFileDS%enable) return

    ! store info on field to be written info at next available table entry

    if (dumpLocal) then
       write(*,"(a)") h//" will insert field "//varn//" at "//oneIOFileDS%fId
    end if
    fieldSize_i8 = int(fieldSize,i8)
    call InsertFieldOnIOHeadTable(oneIOFileDS%ht, &
         varn, idim_type, ngr, oneIOFileDS%fPos, fieldSize_i8)

    ! dump binary field on file

    if (dumpLocal) then
       write(c0,"(i8)") oneIOFileDS%unit
       write(*,"(a)") h//" will write field "//varn//" at unit "//trim(adjustl(c0))
    end if
    call writebin(oneIOFileDS%unit, field, fieldSize)
    oneIOFileDS%fPos = oneIOFileDS%fPos + fieldSize

  end subroutine PointerWriteBinStoreInfo




  subroutine PointerWriteStoreInfo (field, fieldSize, &
       varn, idim_type, ngr, oneIOFileDS)

    integer,           intent(in)    :: fieldSize         ! size of field to write
    real,              pointer       :: field             ! points to field first position
    character(len=16), intent(in)    :: varn              ! field name
    integer,           intent(in)    :: idim_type         ! field dimensionality (coded)
    integer,           intent(in)    :: ngr               ! grid number
    type(IOFileDS),    intent(inout) :: oneIOFileDS

    ! return if disabled 

    if (.not. oneIOFileDS%enable) return

    if (oneIOFileDS%binary) then
       call PointerWriteBinStoreInfo (field, fieldSize, &
            varn, idim_type, ngr, oneIOFileDS)
    else
       call PointerWriteASCIIStoreInfo (field, fieldSize, &
            varn, idim_type, ngr, oneIOFileDS)
    end if
  end subroutine PointerWriteStoreInfo




  subroutine ArrayWriteStoreInfo (field, fieldSize, &
       varn, idim_type, ngr, oneIOFileDS)

    integer,           intent(in)    :: fieldSize         ! size of field to write
    real,              intent(in)    :: field(fieldSize)  ! points to field first position
    character(len=16), intent(in)    :: varn              ! field name
    integer,           intent(in)    :: idim_type         ! field dimensionality (coded)
    integer,           intent(in)    :: ngr               ! grid number
    type(IOFileDS),    intent(inout) :: oneIOFileDS

    ! return if disabled 

    if (.not. oneIOFileDS%enable) return

    if (oneIOFileDS%binary) then
       call ArrayWriteBinStoreInfo (field, fieldSize, &
            varn, idim_type, ngr, oneIOFileDS)
    else
       call ArrayWriteASCIIStoreInfo (field, fieldSize, &
            varn, idim_type, ngr, oneIOFileDS)
    end if
  end subroutine ArrayWriteStoreInfo




  subroutine DumpIOHeadTable(oneIOFileDS)

    type(IOFileDS), intent(inout) :: oneIOFileDS

    integer, parameter :: unitLow=10
    integer, parameter :: unitHigh=99
    integer :: iunit
    integer :: nv
    logical :: op
    character(len=8) :: c0
    character(len=*), parameter :: h="**(DumpIOHeadTable)**"

    ! return if disabled 

    if (.not. oneIOFileDS%enable) return

    ! find available Fortran unit
    
    do iunit = unitLow, unitHigh
       inquire (unit=iunit, opened=op)
       if (.not. op) exit
    end do
    if (iunit <= unitHigh) then
       oneIOFileDS%unit = iunit
    else
       call fatal_error(h//" Fortran i/o units exausted")
    end if

    ! open file, write header information, close file

    call rams_f_open_u(oneIOFileDS%unit, oneIOFileDS%fHeadName,  &
         'FORMATTED','REPLACE','WRITE', "ASIS", oneIOFileDS%iclobber)
    if (dumpLocal) then
       write(c0,"(i8)") oneIOFileDS%unit
       write(*, "(a)") h//" opened new file "//trim(oneIOFileDS%fHeadName)//&
            " at unit "//trim(adjustl(c0))
    end if

    write(oneIOFileDS%unit,'(i6)') oneIOFileDS%ht%lastUsed
    do nv = 1, oneIOFileDS%ht%lastUsed
       write(oneIOFileDS%unit,fmt='(a16,1x,i12,i3,i3,1x,i9)') &
            oneIOFileDS%ht%f(nv)%string,                      &
            oneIOFileDS%ht%f(nv)%npointer,                    &
            oneIOFileDS%ht%f(nv)%idim_type,                   &
            oneIOFileDS%ht%f(nv)%ngrid,                       &
            oneIOFileDS%ht%f(nv)%nvalues
    end do

    call commio(oneIOFileDS%fId,'WRITE',oneIOFileDS%unit)
    close(oneIOFileDS%unit)
    if (dumpLocal) then
       write(c0,"(i8)") oneIOFileDS%unit
       write(*, "(a)") h//" dumped and close file "//trim(oneIOFileDS%fHeadName)//&
            " at unit "//trim(adjustl(c0))
    end if

    oneIOFileDS%unit = -1
  end subroutine DumpIOHeadTable




  subroutine DestroyIOFileDS(oneIOFileDS)
    type(IOFileDS), intent(inout) :: oneIOFileDS

    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(DestroyIOFileDS)**"

    if (dumpLocal) then
       write(*,"(a)") h//" named "//oneIOFileDS%fId
    end if

    if (oneIOFileDS%enable) then
       deallocate (oneIOFileDS%fName, stat=ierr)
       if (ierr /= 0) then
          write(c0,"(i8)") ierr
          call fatal_error(h//" deallocate fName fails with stat="//trim(adjustl(c0)))
       end if
       call DestroyIOHeadTable(oneIOFileDS%ht)
    end if

    oneIOFileDS%enable = .false.
    oneIOFileDS%fId = "    "
  end subroutine DestroyIOFileDS
end module an_header
