!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

   integer function outRealSize()
       !# get the output byte size accordingly the machine
       !#
       !# @note
       !# ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
       !#
       !# **Brief**: Return the output real size (1,4,8,etc) accordingly the machine
       !#
       !# **Documentation**: <http://brams.cptec.inpe.br/documentation/>
       !#
       !# **Author(s)**: Luiz Flavio Rodrigues **&#9993;**<luiz.rodrigues@inpe.br>
       !#
       !# **Date**: 28 July 2020 (Tuesday)
       !# @endnote
       !#
       !# @changes
       !# &#9744; <br/>
       !# @endchanges
       !# @bug
       !#
       !#@endbug
       !#
       !#@todo
       !#  &#9744; <br/>
       !# @endtodo
       !#
       !# @warning
       !# Now is under CC-GPL License, please see
       !# &copy; <https://creativecommons.org/licenses/GPL/2.0/legalcode.pt>
       !# @endwarning
       !#
       
       !Use area
       use dump
   
       implicit none
   
       include "constants.f90"
       character(len=*),parameter :: sourceName='generic.f90' !Name of source code
       character(len=*),parameter :: procedureName='**getOutputByteSize**' !Name of this procedure
       !
       !Local Parameters
   
       !Input/Output variables
   
       !Local variables
       integer :: output_byte_size
       !# output_byte_size
       real :: bytes_in_float
       !# bytes_in_float
   
       !Code
       inquire(iolength=output_byte_size) bytes_in_float
   
       outRealSize=output_byte_size
   
   end function outRealSize 

!***************************************************************************

real function walltime()
  implicit none
  integer :: ii,ir,im
  integer, save :: previous=0
  real, save :: adjustOverflow=0.0
  call system_clock(count=ii,count_rate=ir,count_max=im)
  if (ii < previous) then
    adjustOverflow=adjustOverflow+real(im)/real(ir)
  end if
  previous=ii
  walltime=adjustOverflow + real(ii)/real(ir) ! real(ii)/float(ir)
  return
end function walltime

real function cputime(w1)
  implicit none
  real :: w1
  real :: cc,fsecs
  real, external :: walltime

  call timing(2,cc)
  cputime=cc
  fsecs=72559200.
  w1=walltime()
  return
end function cputime

!***************************************************************************

subroutine rearrange(nzp, nxp, nyp, a, b)
  implicit none
  ! Arguments:
  integer, intent(IN) :: nzp, nxp, nyp
  real, intent(IN)    :: a(nzp,nxp,nyp)
  real, intent(OUT)   :: b(nxp,nyp,nzp)
  ! Local variables:
  integer :: k, i, j

  do i=1,nxp
     do j=1,nyp
        do k=1,nzp
           b(i,j,k) = a(k,i,j)
        enddo
     enddo
  enddo
end subroutine rearrange

!***************************************************************************

subroutine unarrange(nzp, nxp, nyp, a, b)
  implicit none
  ! Arguments:
  integer, intent(IN) :: nzp, nxp, nyp
  real, intent(IN)    :: a(nxp,nyp,nzp)
  real, intent(OUT)   :: b(nzp,nxp,nyp)
  ! Local variables:
  integer :: k, i, j

  do i=1,nxp
     do j=1,nyp
        do k=1,nzp
           b(k,i,j) = a(i,j,k)
        enddo
     enddo
  enddo
end subroutine unarrange

!***************************************************************************

subroutine makefnam (fname,prefix,tinc,iyr,imn,idy,itm,type,post,fmt)
  use ModDateUtils

  ! creates standard timestamped filename

  implicit none

  include "files.h"

  character(len=*), intent(out) :: fname
  character(len=*), intent(in ) :: prefix
  real,             intent(in ) :: tinc
  integer,          intent(in ) :: iyr
  integer,          intent(in ) :: imn
  integer,          intent(in ) :: idy
  integer,          intent(in ) :: itm
  character,        intent(in ) :: type
  character(len=*), intent(in ) :: post
  character(len=3), intent(in ) :: fmt

  integer :: oyr, omn, ody, otm
  integer :: ib1,ib2
  character(len=20) :: dstring
  character(len=8) :: c0, c1, c2, c3, c4
  character(len=16) :: d0
  character(len=*), parameter :: h="**(makefnam)**"
  logical, parameter :: dumpLocal=.false.

  if (dumpLocal) then
     write(d0,"(e16.7)") tinc
     write(c1,"(i8)") iyr
     write(c2,"(i8)") imn
     write(c3,"(i8)") idy
     write(c4,"(i8)") itm
     write(*, "(a)") h//" prefix="//trim(prefix)//"; tinc="//trim(adjustl(d0))//&
          "; iyr="//trim(adjustl(c1))//"; imn="//trim(adjustl(c2))//"; idy="//&
          trim(adjustl(c3))//"; itm="//trim(adjustl(c4))//"; type="//type//&
          "; post="//trim(post)//"; fmt ="//fmt
  end if

  if(tinc == 0.) then
     oyr=iyr ; omn=imn ; ody=idy ; otm=itm
  else
     call date_add_to(iyr,imn,idy,itm,tinc,'s',oyr,omn,ody,otm)
  endif

  write(dstring,100) '-',type,'-',oyr,'-',omn,'-',ody,'-',otm
100 format(3a1,i4.4,a1,i2.2,a1,i2.2,a1,i6.6)

  if (dumpLocal) then
     write(c1,"(i8)") oyr
     write(c2,"(i8)") omn
     write(c3,"(i8)") ody
     write(c4,"(i8)") otm
     write(*, "(a)") h//" oyr="//trim(adjustl(c1))//"; omn="//trim(adjustl(c2))//&
          "; ody="//trim(adjustl(c3))//"; otm="//trim(adjustl(c4))//&
          "; dstring="//dstring
  end if

  ib1=len_trim(prefix)
  fname=prefix(1:ib1)//dstring(1:20)
  if (post(1:1) /= '$') then
     ib1=len_trim(fname)
     ib2=len_trim(post)
     fname=fname(1:ib1)//'-'//post(1:ib2)
  endif
  ib1=len_trim(fname)
  fname=fname(1:ib1)//'.'//fmt(1:3)

  if (dumpLocal) then
     write(*,"(a)") h//" fname="//trim(fname)
  end if
end subroutine makefnam


!***************************************************************************

subroutine rams_f_open(iunit, filenm, formt, stat, act, iclob)

  ! replaces old jclopen and jclget
  ! files are overwritten unless iclob (ICLOBBER) set to 1

  implicit none

  include "files.h"

  integer :: iunit, iclob
  character(len=*) :: formt, stat, act
  character(len=*) :: filenm
  logical :: exans,opnd

  inquire(FILE=filenm(1:len_trim(filenm)),EXIST=exans)

  if(exans.and.iclob.eq.0.and.  &
       (act(1:4).eq.'WRIT'.or.act(1:4).eq.'writ')) then
     print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     print*,'!!!   trying to open file name :'
     print*,'!!!       ',filenm
     print*,'!!!   but it already exists. run is ended.'
     print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     stop 'rams_f_open - exists'
  endif

  open(iunit,STATUS=stat,FILE=filenm(1:len_trim(filenm)),FORM=formt)
!!$  print*,'F_open - ',filenm(1:len_trim(filenm))

  return
end subroutine rams_f_open


subroutine rams_f_open_u(iunit, filenm, formt, stat, act, pos, iclob)

  ! replaces old jclopen and jclget
  ! files are overwritten unless iclob (ICLOBBER) set to 1
  use dump, only: &
    dumpMessage

  implicit none

  include "constants.f90"
  integer,          intent(in) :: iunit
  integer,          intent(in) :: iclob
  character(len=*), intent(in) :: filenm
  character(len=*), intent(in) :: formt
  character(len=*), intent(in) :: stat
  character(len=*), intent(in) :: act
  character(len=*), intent(in) :: pos

  logical :: exans
  integer :: ierr
  character(len=8) :: c0
  character(len=*), parameter :: h="**(rams_f_open_u)**"
  character(len=*), parameter :: header="**(rams_f_open_u)**"

  inquire(FILE=filenm(1:len_trim(filenm)),EXIST=exans)

  if(exans.and.iclob.eq.0.and.  &
       (act(1:4).eq.'WRIT'.or.act(1:4).eq.'writ')) then
     print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     print*,'!!!   trying to open file name :'
     print*,'!!!       ',filenm
     print*,'!!!   but it already exists. run is ended.'
     print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     stop 'rams_f_open - exists'
  endif

  open(iunit,FILE=filenm(1:len_trim(filenm)), &
       FORM=formt, STATUS=stat, ACTION=act, POSITION=pos, IOSTAT=ierr)
  if (ierr /= 0) then
     write(c0, "(i8)") ierr
     !call fatal_error(h//" open file "//trim(filenm)//", FORM="//formt//&
    !      ", STATUS="//stat//", ACTION="//act//", POSITION="//pos//&
    !      " fails with IOSTAT="//trim(adjustl(c0)))
     iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion,c_fatal, &
          " open file "//trim(filenm)//", FORM="//formt//&
          ", STATUS="//stat//", ACTION="//act//", POSITION="//pos//&
          " fails with IOSTAT="//trim(adjustl(c0)))
  end if
end subroutine rams_f_open_u

integer function AvailableFileUnit()
  use dump, only: &
    dumpMessage

  implicit none

  include "constants.f90"
  integer, parameter :: firstUnit=20  ! lowest io unit number available
  integer, parameter :: lastUnit=99   ! highest io unit number available

  integer :: iunit
  logical :: op
  character(len=*), parameter :: h="**(AvailableFileUnit)**"
  character(len=*), parameter :: header="**(AvailableFileUnit)**"

  ! select unused i/o unit

  do iunit = firstUnit, lastUnit
     inquire(iunit,opened=op)
     if (.not. op) exit
  end do

  if (iunit > lastUnit) then
     !call fatal_error(h//" all i/o units in use")
     iErrNumber=dumpMessage(c_tty,c_yes,header,modelVersion,c_fatal, &
          "all i/o units in use")
  else
     AvailableFileUnit = iunit
  end if
end function AvailableFileUnit
