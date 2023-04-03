module modUtils
    use modMemory

    use dump, only: dumpMessage
    implicit none

    include "constants.h"
    character(len=*),parameter :: fname='./namelist'
    character(len=*),parameter :: srcFile='modUtils.f90'

contains

subroutine readNamelist()
    
    character(len=*),parameter :: header='**(readnamelist)**'
    integer, parameter :: iunit=33
    logical :: ex
    integer :: err

    namelist /MODEL/     &
        prefix,outfolder,imonth1,idate1,iyear1,itime1,ntimes,tincrem,source

    inquire(file=trim(fname), exist=ex)
    if(.not. ex) iErrNumber=dumpMessage(c_tty,c_yes,header,srcFile &
              ,c_fatal,trim(fname)//' not found. Please, check it!')

    open(iunit, file=trim(fname), status="old", action="read",&
         iostat=err)
    write(*,fmt='("Reading Namelist from ",A)') trim(fname)
    read (iunit, iostat=err, NML=MODEL)
    close(iunit)

    write(*,fmt='(A,A)')  'prefix    : ',trim(prefix)
    write(*,fmt='(A,I2)') 'imonth1   : ',imonth1
    write(*,fmt='(A,I2)') 'idate1    : ',idate1 
    write(*,fmt='(A,I4)') 'iyear1    : ',iyear1 
    write(*,fmt='(A,I4)') 'itime1    : ',itime1 
    write(*,fmt='(A,I3)') 'ntimes    : ',ntimes 
    write(*,fmt='(A,I2)') 'tincrem   : ',tincrem
    write(*,fmt='(A,A4)') 'source    : ',source
    write(*,fmt='(A,A)')  'outFolder : ',trim(outFolder)
    print *,''
    lastTime=(nTimes-1)*tincrem
    write(*,fmt='(A,I4)') 'lastTime  : ',lastTime
    print *,''

end subroutine readNamelist


subroutine init(timeDivision,lastTime)
    use mpi

    integer, intent(in)  :: timeDivision
    integer, intent(in)  :: lastTime

    integer :: i
    real :: fpp,totf,np,count

    call MPI_Init ( ierr )

    call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )

    call MPI_Comm_size ( MPI_COMM_WORLD, p, ierr )

    allocate(nFilesPerProc(0:p-1))
    allocate(ini(0:p-1))
    allocate(fim(0:p-1))

    if(id==0) then
        totF=lastTime/timeDivision+1
        fpp=(real(lastTime/timeDivision)+1)/real(p)
        !print *,'fpp=',fpp
        do i=0,p-1
            nFilesPerProc(i)=int(fpp)
        enddo
    
        do i=1,ceiling((fpp-int(fpp))*p)
            nFilesPerProc(i-1)=nFilesPerProc(i-1)+1
        enddo
    
        ini(0)=0
        fim(0)=ini(0)+timeDivision*(nFilesPerProc(0)-1)
        do i=1,p-1
            ini(i)=fim(i-1)+timeDivision
            fim(i)=ini(i)+timeDivision*(nFilesPerProc(i)-1)
        enddo
    
        do i=0,p-1
            write(*,fmt='("nFilesPerProc(",I2.2,")=",3(I3.3,1X))') i,nFilesPerProc(i),ini(i),fim(i)
        enddo
    endif
    
    call MPI_BCAST(nFilesPerProc, p, MPI_INTEGER, 0, &
           MPI_COMM_WORLD, ierr) 
    call MPI_BCAST(ini, p, MPI_INTEGER, 0, &
           MPI_COMM_WORLD, ierr) 
    call MPI_BCAST(fim, p, MPI_INTEGER, 0, &
           MPI_COMM_WORLD, ierr) 

end subroutine init

integer function encerra()
    call MPI_Finalize ( encerra )
end function encerra

subroutine str2int(str,int,stat)
  implicit none
  ! Arguments
  character(len=*),intent(in) :: str
  integer,intent(out)         :: int
  integer,intent(out)         :: stat
  read(str,*,iostat=stat)  int
end subroutine str2int

character(len=3) function monthName(cMonth)
    character(len=*),parameter :: header='**(monthName)**'
    character(len=2), intent(in) :: cMonth
    select case (cMOnth)
        case("01")
            monthName='jan'
        case("02")
            monthName='feb'    
        case("03")
            monthName='mar'  
        case("04")
            monthName='apr'  
        case("05")
            monthName='may'  
        case("06")
            monthName='jun'  
        case("07")
            monthName='jul'  
        case("08")
            monthName='aug'  
        case("09")
            monthName='sep'  
        case("10")
            monthName='oct'  
        case("11")
            monthName='nov'  
        case("12")
            monthName='dec' 
        case default
            iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion &
                ,c_fatal,'month '//cMOnth//' no exist!')
    end select

end function monthName

end module modUtils