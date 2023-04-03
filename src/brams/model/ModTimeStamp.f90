module ModTimeStamp

  ! PURPOSE:
  !
  ! Stores wall clock execution time stamps at selected events of an MPI run.
  ! Recorded time stamps are dump to a file, for off-line post-processing.
  !
  ! USAGE:
  !
  ! Select source code regions to be measured. Enumerate regions.
  ! Assign a name to each region. 
  !
  ! Invoke CreateTimeStamp at the beginning of the run 
  ! (just after MPI_INIT), passing region names as input argument.
  !
  ! Invoke TimeStamp (or SynchronizedTimeStamp) at the end of 
  ! each region, passing region enumeration as input argument. 
  ! 
  ! Invoke DestroyTimeStamp at the end of the run (just before
  ! MPI_FINALIZE) to dump measured execution times to a file.
  !
  ! Run produces one time stamp file for each MPI process and
  ! a directory file, to be post-processed.
  !
  ! It is central that all processes execute the same set of time stamp
  ! procedures, since each SynchronizedTimeStamp invocation has a build-in
  ! global barrier. 
  !
  ! METHOD:
  !
  ! Each call to TimeStamp records wall clock at the point of measurement.
  ! Each call to SynchronizedTimeStamp records two wall clocks: first one
  ! just before a build-in MPI_BARRIER and the second one just after the barrier.
  ! CreateTimeStamp records execution time at the beginning of the run.
  ! DestroyTimeStamp records execution time at the end of the run and
  ! dumps execution time measurements to a file.
  !
  ! Time difference among consecutive measurements is assigned to region id.
  ! Post-processing adds up all time differences assigned to a region id
  ! and produces a report, one line per process.


  implicit none

  private

  include "i8.h"
  integer              :: count_rate                  ! ticks per second
  integer              :: count_max                   ! measurement overflow
  integer              :: nodes                       ! MPI size
  integer              :: pId                         ! MPI rank
  integer, parameter   :: sizeTS = 2000000             ! array size
  integer              :: lastTS                      ! last entry used
  integer(kind=i8), allocatable :: TS(:)                       ! time stamps indexed by event
  integer, allocatable :: EventId(:)                  ! time stamp event identifier
  integer, parameter   :: sizeEventName = 100         ! maximum number of events
  integer              :: lastEventName               ! number of distinct events
  character(len=12), allocatable :: EventName(:)      ! name of distinct events
  character(len=4)     :: cpId                        ! MPI rank in characters
  integer, parameter   :: unitDir=59                  ! unit directory file name
  character(len=18)    :: fNameDir                    ! directory file name 
  logical, parameter   :: dumpLocal=.false.           ! module debug (set to .true. for dumping exec info)
  logical, parameter   :: noInstrumentation=.false.    ! disable instrumentation (if set to true)

  public :: CreateTimeStamp
  public :: TimeStamp
  public :: SynchronizedTimeStamp
  public :: DestroyTimeStamp

  

contains



  ! CreateTimeStamp: module initialization
  !                  number of events is the size of input array "names"




  subroutine CreateTimeStamp(machsize, machnum, names)
    integer,          intent(in) :: machsize  ! MPI size
    integer,          intent(in) :: machnum   ! MPI rank
    character(len=*), intent(in) :: names(:)  ! names of events

    integer :: ierr
    integer :: name
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(CreateTimeStamp)**"

    ! null action if disabled instrumentation

    if (noInstrumentation) return

    ! module initialization

    allocate(TS(0:sizeTS), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" allocate TS fails with stat="//trim(adjustl(c0)))
       stop
    end if
    allocate(EventId(0:sizeTS), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" allocate EventId fails with stat="//trim(adjustl(c0)))
       stop
    end if
    allocate(EventName(sizeEventName), stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" allocate EventName fails with stat="//trim(adjustl(c0)))
       stop
    end if
    lastTS=-1
    lastEventName = size(names)
    nodes = machsize
    pId = machnum
    write(cpId,"(i4.4)") machnum

    ! names

    if (lastEventName > sizeEventName) then
       write(c0, "(i8)") sizeEventName
       write(c1, "(i8)") lastEventName
       call fatal_error (h//" size of names ("//trim(adjustl(c1))//&
            &") exceeds maximum ("//trim(adjustl(c1))//")")
       stop
    end if
    EventName(1:lastEventName) = names(1:lastEventName)
    
    ! base measurement

    call timeBlocked(0)

    ! directory file header

    if (pId == 0) then
       fNameDir="ts.XXXXX.Directory"
       write(fNameDir(4:8),"(i5.5)") machsize
       open  (unitDir, file=fNameDir, status="replace", action="write")
       write (unitDir, "(i12)") count_rate
       write (unitDir, "(i12)") count_max
       write (unitDir, "(i12)") lastEventName
       do name = 1, lastEventName
          write (unitDir, "(a)") EventName(name)
       end do
       close (unitDir)
    end if

    ! dump

    if (dumpLocal) then
       write(c0, "(i8)") machsize
       write(c1, "(i8)") lastEventName
       write(*,*) h//" done at proc "//cpId//" of "//trim(adjustl(c0))//" with "//trim(adjustl(c1))//" events"
    end if
  end subroutine CreateTimeStamp



  ! TimeStamp: store a time stamp to event "event"



  subroutine TimeStamp(event)
    integer, intent(in) :: event  ! event id

    integer :: ii
    integer, save :: previous=0
    integer(i8), save :: adjustOverflow=0_i8
    character(len=16) :: c0, c1, c2
    character(len=*), parameter :: h="**(TimeStamp)**"

    ! null action if disabled instrumentation

    if (noInstrumentation) return

    ! measure time accountig for overflow

    call system_clock(ii, count_rate, count_max)
    if (ii < previous) then
       adjustOverflow = adjustOverflow + count_max
    end if
    previous = ii

    ! store time stamp at next available position

    lastTS = lastTS + 1
    if (lastTS > sizeTS) then
       write(c0,"(i8)") lastTS
       write(c1,"(i8)") sizeTS
       call fatal_error(h//" more TimeStamps to store ("//&
            trim(adjustl(c0))//") than antecipated ("//trim(adjustl(c1))//")")
       stop
    else
       TS(lastTS) = int(ii,i8) + adjustOverflow
       EventId(lastTS) = event
    end if

    ! dump

    if (dumpLocal) then
       write(c0,"(i8)") event
       write(c1,"(i16)") TS(lastTS)
       write(c2,"(i8)") lastTS
       write(*,*) h//" at Proc "//cpId//" time stamp #"//trim(adjustl(c2))//&
            " with event "//trim(adjustl(c0))//" at time "//trim(adjustl(c1))
    end if
  end subroutine TimeStamp


  ! timeBlocked: module private routine
  !              records a time stamp after synchronization


  subroutine timeBlocked(event)
    integer, intent(in) :: event
    integer :: ierr
    
    character(len=8) :: c0
    character(len=*), parameter :: h="**(timeBlocked)**"
    include "mpif.h"

    ! barrier followed by time stamp

    if (dumpLocal) then
       write(c0,"(i8)") event
       write(*,*) h//" proc "//cpId//" will wait on barrier of event "//trim(adjustl(c0))
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call TimeStamp(event)
  end subroutine timeBlocked



  ! SynchronizedTimeStamp: records two time stamps, first on procedure entrance,
  !                        second after a global barrier. Time stamp difference
  !                        on the two measures represent synchronization time.
  !                        Input argument is event Id.



  subroutine SynchronizedTimeStamp(event)
    integer, intent(in) :: event

    ! null action if disabled instrumentation

    if (noInstrumentation) return

    call TimeStamp(event)
    call timeBlocked(-1*event)
  end subroutine SynchronizedTimeStamp



  ! DestroyTimeStamp: Dumps measured time stamps to a file and
  !                   re-initialize module



  subroutine DestroyTimeStamp()
    character(len=16)   :: fName
    integer             :: i
    integer             :: ierr
    character(len=8) :: c0, c1
    character(len=*), parameter :: h="**(DestroyTimeStamp)**"

    ! null action if disabled instrumentation

    if (noInstrumentation) return

    ! finish directory file 

    if (pId == 0) then
       open  (unitDir, file=fNameDir, status="old", action="write", position="append")
       write (unitDir, "(i12)") lastTS
       close (unitDir)
    end if

    ! File Name

    fName = "tsXXXXXYYYYY.out"
    write(fName(3: 7),"(i5.5)") nodes
    write(fName(8:12),"(i5.5)") pId

    ! File Containts

    open  (unitDir, file=fName, status="replace", form="unformatted", action="write")
    write (unitDir) lastTS
    do i = 0, lastTS
       write (unitDir) TS(i), EventId(i)
    end do
    close (unitDir)

    ! dump

    if (dumpLocal) then
       write(*,*) h//" Proc "//cpId//" finishes "
    end if
    deallocate(TS, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" deallocate TS fails with stat="//trim(adjustl(c0)))
       stop
    end if
    deallocate(EventId, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" deallocate EventId fails with stat="//trim(adjustl(c0)))
       stop
    end if
    deallocate(EventName, stat=ierr)
    if (ierr /= 0) then
       write(c0,"(i8)") ierr
       call fatal_error(h//" deallocate EventName fails with stat="//trim(adjustl(c0)))
       stop
    end if
  end subroutine DestroyTimeStamp
end module ModTimeStamp
