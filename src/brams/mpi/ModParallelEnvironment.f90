module ModParallelEnvironment

  ! ModParallelEnvironment: creates a type that stores info about the 
  !                         MPI parallel environment. Creates a dump file
  !                         for each process. Provides procedures for
  !                         dumping a msg at the dump file or at stdout.

  use ISO_FORTRAN_ENV
  private
  public :: ParallelEnvironment
  public :: CreateParallelEnvironment
  public :: DestroyParallelEnvironment
  public :: MsgDump
  public :: MsgOutput
  public :: Brams2MpiProcNbr
  public :: Mpi2BramsProcNbr
  public :: GetNumberOfProcesses
  public :: GetThisBramsProcessNumber

  type ParallelEnvironment
     integer :: communicator 
     integer :: nmachs             ! number of processes (MPI size)
     integer :: mchnum             ! MPI number of this process
     integer :: master_num         ! MPI number of the master process
     integer :: myNum              ! Brams number of this process
  end type ParallelEnvironment

  include "mpif.h"
  integer, parameter :: DumpUnit=21

contains



  !*** CreateParallelEnvironment: create and fill variable of this type.
  !                               first call opens dump file for one process



  subroutine CreateParallelEnvironment(nmachs, mchnum, master_num, &
       comm, oneParallelEnvironment)
    integer, intent(in) :: nmachs     ! number of processes (0 iff sequential run)
    integer, intent(in) :: master_num ! this process rank (0:nmachs-1); 0 on sequential runs
    integer, intent(in) :: mchnum     ! this process rank (0:nmachs-1); 0 on sequential runs
    integer, intent(in) :: comm       ! MPI communicator

    logical :: op
    type(ParallelEnvironment), pointer :: oneParallelEnvironment

    ! creates variable and fill components

    allocate(oneParallelEnvironment)
    oneParallelEnvironment%nmachs = nmachs
    oneParallelEnvironment%mchnum = mchnum
    oneParallelEnvironment%master_num = master_num
    oneParallelEnvironment%communicator = comm
    oneParallelEnvironment%myNum = Mpi2BramsProcNbr(mchnum)

  end subroutine CreateParallelEnvironment



  ! DestroyParallelEnvironment: returns storage, destroy variable
  !                             closes dump file



  subroutine DestroyParallelEnvironment(oneParallelEnvironment)
    type(parallelEnvironment), pointer :: oneParallelEnvironment

    logical :: op
    if (associated(oneParallelEnvironment)) then
       deallocate(oneParallelEnvironment)

       ! only the first execution of this procedure closes dump file

       inquire(unit=DumpUnit, opened=op)
       if (op) then
          close(DumpUnit)
       end if
    end if
    nullify(oneParallelEnvironment)
  end subroutine DestroyParallelEnvironment



  ! MsgDump: Writes string at dump file



  subroutine MsgDump(str, noAdvance)
    character(len=*), intent(in) :: str
    logical, optional, intent(in) :: noAdvance

    integer :: mchnum
    integer :: nmachs
    integer :: ierr
    logical :: op
    character(len=16) :: dumpFName="Dump.XXXXX.YYYYY"

    ! only the first execution of this procedure opens dump file

    inquire(unit=DumpUnit, opened=op)
    if (.not. op) then
       call MPI_Comm_rank(MPI_COMM_WORLD, mchnum, ierr)
       if(ierr /= MPI_SUCCESS) then
          call fatal_error("Error in MPI_Comm_rank")
       endif

       call MPI_Comm_size(MPI_COMM_WORLD, nmachs, ierr)
       if(ierr /= MPI_SUCCESS) then
          call fatal_error("Error in MPI_Comm_size")
       endif

       write(dumpFName(6:10),"(i5.5)") mchnum
       write(dumpFName(12:16),"(i5.5)") nmachs
       open(DumpUnit, file=dumpFName)
    end if

    if (present(noAdvance)) then
       if (noAdvance) then
          write(DumpUnit,"(a)",advance="no") trim(str)
          return
       end if
    end if
    write(DumpUnit,"(a)") trim(str)
    flush(DumpUnit)
  end subroutine MsgDump



  ! MsgOutput: Writes string at output file



  subroutine MsgOutput(str, noAdvance)
    character(len=*), intent(in) :: str
    logical, optional, intent(in) :: noAdvance

    if (present(noAdvance)) then
       if (noAdvance) then
          write(OUTPUT_UNIT,"(a)",advance="no") trim(str)
          return
       end if
    end if
    write(OUTPUT_UNIT,"(a)") trim(str)
  end subroutine MsgOutput




  integer function Brams2MpiProcNbr(procNbr)
    integer, intent(in) :: procNbr
    Brams2MpiProcNbr = procNbr-1
  end function Brams2MpiProcNbr




  integer function Mpi2BramsProcNbr(procNbr)
    integer, intent(in) :: procNbr
    Mpi2BramsProcNbr = procNbr+1
  end function Mpi2BramsProcNbr




  integer function GetNumberOfProcesses(oneParallelEnvironment)
    type(parallelEnvironment), pointer :: oneParallelEnvironment

    character(len=*), parameter :: h="**(GetNumberOfProcesses)**"

    if (.not. associated(oneParallelEnvironment)) then
       call fatal_error(h//" oneParallelEnvironment not associated")
    else
       GetNumberOfProcesses = oneParallelEnvironment%nmachs
    end if
  end function GetNumberOfProcesses




  integer function GetThisBramsProcessNumber(oneParallelEnvironment)
    type(parallelEnvironment), pointer :: oneParallelEnvironment

    character(len=*), parameter :: h="**(GetThisBramsProcessNumber)**"

    if (.not. associated(oneParallelEnvironment)) then
       call fatal_error(h//" oneParallelEnvironment not associated")
    else
       GetThisBramsProcessNumber = oneParallelEnvironment%myNum
    end if
  end function GetThisBramsProcessNumber
end module ModParallelEnvironment

