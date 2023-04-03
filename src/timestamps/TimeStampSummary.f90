module InputTimeStamp

  implicit none

  private

  include "i8.h"

  type directory
     character(len=128)         :: fName
     integer                    :: MPISize
     integer                    :: count_rate
     integer                    :: count_max
     integer                    :: lastEventName
     character(len=12), pointer :: EventName(:)
     integer                    :: lastTS
  end type directory

  type TSFile
     character(len=128)  :: fName
     integer             :: MPIRank
     integer             :: lastTS
     integer(kind=i8), pointer     :: TS(:)
     integer, pointer    :: EventId(:)
     integer             :: eventMin
     integer             :: eventMax
     integer(kind=i8)    :: tsMax
  end type TSFile


  type Report
     character(len=256) :: fName
     integer            :: sizeNames
     integer, pointer   :: mapEventToPos(:)
  end type Report

  integer, parameter :: unitDir=59

  public :: directory
  public :: ReadDirectory
  public :: DestroyDirectory
  public :: TSFile
  public :: ReadTSFile
  public :: DumpTSFile
  public :: DestroyTSFile
  public :: Report
  public :: CreateReport
  public :: ReportOneTSFile
  public :: DestroyReport
contains


  subroutine ReadDirectory (machsize, dirName, dir)
    integer, intent(in) :: machsize
    character(len=*), intent(in) :: dirName
    type(directory), intent(out) :: dir

    integer :: name
    integer :: lastChar

    dir%fName="ts.XXXXX.Directory"
    write(dir%fName(4:8),"(i5.5)") machsize
    lastChar = len_trim(dirName)
    if (dirName(lastChar:lastChar) == "/") then
       dir%fName= dirName(1:lastChar)//trim(dir%fName)
    else
       dir%fName= dirName(1:lastChar)//"/"//trim(dir%fName)
    end if

    dir%MPISize = machsize

    open (unitDir, file=trim(dir%fName), status="old", action="read")
    read (unitDir, "(i12)") dir%count_rate
    read (unitDir, "(i12)") dir%count_max
    read (unitDir, "(i12)") dir%lastEventName
    allocate(dir%EventName(dir%lastEventName))
    do name = 1, dir%lastEventName
       read (unitDir, "(a)") dir%EventName(name)
    end do
    read (unitDir, "(i12)") dir%lastTS
    print *,dir%count_rate
    print *,dir%count_max
    print *,dir%lastEventName
    close (unitDir)
  end subroutine ReadDirectory


  subroutine DestroyDirectory(dir)
    type(directory), intent(inout) :: dir
    deallocate(dir%EventName)
  end subroutine DestroyDirectory



  subroutine ReadTSFile (machsize, machrank, dirName, ts)
    integer, intent(in) :: machsize
    integer, intent(in) :: machrank
    character(len=*), intent(in) :: dirName
    type(TSFile), intent(out) :: ts
    integer :: i
    integer :: lastChar


    ts%fName = "tsXXXXXYYYYY.out"
    write(ts%fName(3: 7),"(i5.5)") machsize
    write(ts%fName(8:12),"(i5.5)") machrank
    lastChar = len_trim(dirName)
    if (dirName(lastChar:lastChar) == "/") then
       ts%fName= dirName(1:lastChar)//trim(ts%fName)
    else
       ts%fName= dirName(1:lastChar)//"/"//trim(ts%fName)
    end if
    print *,ts%fName
    ts%MPIRank=machrank

    ! File Containts

    open (unitDir, file=trim(ts%fName), status="old", form="unformatted", action="read")
    read (unitDir) ts%lastTS
    allocate(ts%TS(0:ts%lastTS))
    allocate(ts%EventId(0:ts%lastTS))
    do i = 0, ts%lastTS
       read (unitDir) ts%TS(i), ts%EventId(i)
       print *,i,ts%TS(i),ts%EventId(i)
    end do
    close (unitDir)

    ! extremes

    ts%eventMin = minval(ts%EventId(:))
    ts%eventMax = maxval(ts%EventId(:))
    ts%tsMax    = ts%TS(ts%lastTS)
  end subroutine ReadTSFile




  subroutine DumpTSFile(ts)
    type(TSFile), intent(in) :: ts

    integer :: i
    integer :: event
    integer, allocatable :: eventDist(:)
    character(len=*), parameter :: h="**(DumpTSFile)**"
    character(len=16) :: c0, c1

    write(c0,"(i16)") ts%lastTS
    write(c1,"(i16)") ts%tsMax
    write(*,"(a)") h//" file "//trim(ts%fName)//" has "//trim(adjustl(c0))//&
         " time stamps up to time "//trim(adjustl(c1))

    allocate (eventDist(ts%eventMin:ts%eventMax))
    eventDist(:) = 0
    do i = 0, ts%lastTS
       event = ts%EventId(i)
       eventDist(event) = eventDist(event) + 1
    end do

    write(*,"(a)",advance="no") h//" event:      "
    do event = ts%eventMin, ts%eventMax
       write(*, "(i6)", advance="no") event
    end do
    write(*,"(1x)")
    write(*,"(a)",advance="no") h//" occurences: "
    do event = ts%eventMin, ts%eventMax
       write(*, "(i6)", advance="no") eventDist(event)
    end do
    write(*,"(1x)")

    deallocate(eventDist)
  end subroutine DumpTSFile




  subroutine DestroyTSFile(ts)
    type(TSFile), intent(inout) :: ts
    deallocate(ts%TS)
    deallocate(ts%EventId)
  end subroutine DestroyTSFile




  subroutine CreateReport (dir, dirName, ts, rep)
    type(directory),  intent(in) :: dir
    character(len=*), intent(in) :: dirName
    type(TSFile),     intent(in) :: ts
    type(Report),     intent(out) :: rep

    integer :: iTS
    integer :: iTSFile
    integer :: lastChar
    integer :: i
    integer :: pos
    integer :: eventNumber
    integer :: ierr
    integer :: posNames
    integer :: posEventNames
    integer :: sizeNames
    character(len=9), allocatable :: names(:)
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(CreateReport)**"

    ! report file name

    rep%fName="Summary.XXXXX"
    write(rep%fName(9:13), "(i5.5)") dir%MPISize
    lastChar = len_trim(dirName)
    if (dirName(lastChar:lastChar) == "/") then
       rep%fName = dirName(1:lastChar)//trim(rep%fName)
    else
       rep%fName = dirName(1:lastChar)//"/"//trim(rep%fName)
    end if

    ! open report file

    open(unit=unitDir, file=trim(rep%fName), status="replace", action="write", iostat=ierr)
    if (ierr == 0) then
       write(*,"(a)") h//" Building file "//trim(rep%fName)
    else
       write(c0,"(i8)") ierr
       write(*,"(a)") h//"**(ERROR)** open file "//trim(rep%fName)//&
            &" returns iostat="//trim(adjustl(c0))
       stop
    end if

    ! event names on report's first line printing position

    rep%sizeNames = dir%lastEventName + 1  ! events, total
    allocate(names(rep%sizeNames), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       write(c1,"(i8)") rep%sizeNames
       write(*,"(a)") h//"**(ERROR)** allocate names("//trim(adjustl(c1))//&
            &") returns iostat="//trim(adjustl(c0))
       stop
    end if

    ! map event id to printing position

    allocate(rep%mapEventToPos(-dir%lastEventName:dir%lastEventName), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       write(c1,"(i8)") dir%lastEventName
       write(*,"(a)") h//"**(ERROR)** allocate mapEventToPos(-"//trim(adjustl(c1))//&
            &":"//trim(adjustl(c1))//") returns iostat="//trim(adjustl(c0))
       stop
    end if
    rep%mapEventToPos = 0

    ! position events (and possible synchronizations) for printing

    posNames=0
    do posEventNames = 1, dir%lastEventName  ! skip event 0 (baseline)
       posNames = posNames+1
       names(posNames) = dir%EventName(posEventNames)
       rep%mapEventToPos(posEventNames) = posNames
       rep%mapEventToPos(-posEventNames) = posNames
    end do
    posNames = posNames + 1
    if (posNames == rep%sizeNames) then
       names(posNames) = "Total"
    else
       write(c0,"(i8)") posNames
       write(c1,"(i8)") rep%sizeNames
       write(*,"(a)") h//"**(ERROR)** logical error: posNames ("//trim(adjustl(c0))//&
            &") and sizeNames ("//trim(adjustl(c1))//") differ"
       stop
    end if

    ! header line

    write(unitDir,"(a4)", advance="no") "proc"
    do i = 1, size(names)
       write(unitDir,"(1x,a9)", advance="no") adjustr(names(i))
    end do
    write(unitDir,"(1x)")

    close(unitDir)
    deallocate(names)
  end subroutine CreateReport




  subroutine ReportOneTSFile (dir, ts, rep)
    type(directory), intent(in) :: dir
    type(TSFile),    intent(in) :: ts
    type(Report),    intent(in) :: rep
    
    integer :: i
    integer :: iTS
    integer :: pos
    integer :: ierr
    integer(kind=i8), allocatable :: counters(:)
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(ReportOneTSFile)**"

    ! open report file
    write(*, "(a,i6.6)") h//" for proc ", ts%MPIrank

    open(unit=unitDir, file=trim(rep%fName), status="old", action="write", position="append", iostat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       write(*,"(a)") h//"**(ERROR)** open file "//trim(rep%fName)//&
            &" to append returns iostat="//trim(adjustl(c0))
       stop
    end if

    ! process TS file

    allocate(counters(1:rep%sizeNames))
    counters = 0_i8
    do iTS = 1, ts%lastTS
       pos = rep%mapEventToPos(ts%EventId(iTS))
       counters(pos) = counters(pos) + ts%TS(iTS) - ts%TS(iTS-1)
    end do
    counters(rep%sizeNames) = ts%TS(iTS-1) - ts%TS(0)

    write(unitDir,"(i4.4)",advance="no") dir%MPISize
    do i = 1, size(counters)
       write(unitDir,"(1x,f9.2)", advance="no") real(counters(i))/real(dir%count_rate)
    end do
    write(unitDir,"(1x)")

    ! finalize

    deallocate(counters)
    close(unitDir)
  end subroutine ReportOneTSFile




  subroutine DestroyReport(rep)
    type(Report),     intent(out) :: rep
    deallocate(rep%mapEventToPos)
  end subroutine DestroyReport
end module InputTimeStamp




program TimeStampReport

  use InputTimeStamp
  implicit none

  ! namelist 

  integer            :: machsize ! MPI size
  character(len=128) :: pathIn   ! data directory
  character(len=128) :: pathOut  ! output directory
  namelist /TSInput/ machsize, pathIn, pathOut

  integer :: mach
  integer :: ierr
  character(len=8)  :: c0, c1
  character(len=*), parameter :: h="**(LeTimeStamp)**" 

  type(directory) :: dir
  type(TSFile) :: ts
  type(Report) :: rep

  ! read namelist and echo 

  read (*, NML=TSInput, iostat=ierr)
  if (ierr == 0) then
     write(*, "(a,i4.4)") " MPI rank is ",machsize
     write(*, "(a)") " Input Directory: "//trim(pathIn)
     write(*, "(a)") " Output Directory: "//trim(pathOut)
  else
     write(c0,"(i8)") ierr
     write(*,"(a)") "**(ERROR)** read namelist returns istat="&
          &//trim(adjustl(c0))
     stop
  end if

  ! read directory

  call ReadDirectory(machsize, pathIn, dir)

  ! read proc 0 trace file; write first line and proc 0 line; destroy trace file

  call ReadTSFile (machsize, 0, pathIn, ts)
  call CreateReport (dir, pathOut, ts, rep)
  call ReportOneTSFile (dir, ts, rep)
  call DestroyTSFile(ts)

  ! read one trace file, report and destroy

! do mach = 1, machsize-1
!    call ReadTSFile (machsize, mach, pathIn, ts)
!    call ReportOneTSFile (dir, ts, rep)
!    call DestroyTSFile(ts)
! end do

  ! finalize

  call DestroyDirectory(dir)
  call DestroyReport(rep)
  stop
end program TimeStampReport
