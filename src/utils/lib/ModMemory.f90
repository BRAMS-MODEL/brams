module ModMemory
!  use mpi
  implicit none
include 'mpif.h' !if using intel at tupa swap "use mpi" for this.
  private
  public :: CreateMemory
  public :: DumpMemory
  public :: DestroyMemory

  integer, parameter :: unitDump=77
  logical, parameter :: dumpLocal=.false.
  logical, parameter :: noInstrumentation=.false.   ! disable instrumentation (if set to true)
  logical :: thisProcessDumps=.false.
contains





  subroutine CreateMemory(mchnum, nmachs, master_num)
    integer, intent(in) :: nmachs
    integer, intent(in) :: mchnum
    integer, intent(in) :: master_num
    character(len=8) :: c0

    if (noInstrumentation) return

    thisProcessDumps = mchnum == master_num
    if (thisProcessDumps) then
       write(c0,"(i8)") nmachs
       open(unitDump, file="DumpMemory."//trim(adjustl(c0)), action="write", status="replace")
    end if
  end subroutine CreateMemory





  subroutine DumpMemory(str)
    character(LEN=*), intent(IN) :: str

    integer, parameter :: mem_unit=44
    integer :: myId, ierr
    integer :: memHWM, memRSS
    integer :: memTotal, memFree, memUsed
    character(len=*), parameter :: h="**(DumpMemory)**"
    character(len=*), parameter :: file1 = '/proc/self/status'
    character(len=*), parameter :: file2 = '/proc/meminfo'
    integer, parameter :: lenString=256
    character(len=lenString) :: string
    character(LEN=8) :: c0, c1, c2, c3, c4

    if (noInstrumentation) return

    if (thisProcessDumps) then

       call MPI_Comm_rank(MPI_COMM_WORLD, myId, ierr)

       open (unit=mem_unit, file=file1, status="old", action="read")
       do
          read (mem_unit,'(a)', end=10) string
          if (index(string,'VmHWM:') == 1) then
             read (string(7:lenString-2),*) memHWM
             if (dumpLocal) write(*,*) h//" VmHWM=",memHWM,trim(string)
          else if (index(string,'VmRSS:') == 1) then
             read (string(7:lenString-2),*) memRSS
             if (dumpLocal) write(*,*) h//" VmRSS=",memRSS,trim(string)
             exit
          end if
       end do
10     continue
       close(mem_unit)

       open (unit=mem_unit, file=file2, status="old", action="read")
       do
          read (mem_unit,'(a)', end=20) string
          if (index(string,'MemTotal:') == 1) then
             read (string(10:lenString-2),*) memTotal
             if (dumpLocal) write(*,*) h//" memTotal=",memTotal
          else if (index(string,'MemFree:') == 1) then
             read (string(9:lenString-2),*) memFree
             if (dumpLocal) write(*,*) h//" memFree=",memFree
             exit
          end if
       end do
20     continue
       close(mem_unit)

       memUsed = memTotal-memFree
       write(c0,"(i8)") myId
       write(c1,"(i8)") memRSS/1024
       write(c2,"(i8)") memHWM/1024
       write(c3,"(i8)") memUsed/1024
       write(c4,"(i8)") memFree/1024

       write(*,"(a)") h//trim(adjustl(str))//" proc "//trim(adjustl(c0))//&
            "; RSS(MB): cur="//trim(adjustl(c1))//&
            " hwm="//trim(adjustl(c2))// &
            "; mem(MB): used="//trim(adjustl(c3))// &
            " free="//trim(adjustl(c4))

       write(unitDump,"(4(i5.5,1x),a)") &
            memRSS/1024, &
            memHWM/1024, &
            memUsed/1024, &
            memFree/1024, &
            trim(adjustl(str))
    end if
  end subroutine DumpMemory





  subroutine DestroyMemory()
    if (noInstrumentation) return
    close(unitDump)
  end subroutine DestroyMemory
end module ModMemory
