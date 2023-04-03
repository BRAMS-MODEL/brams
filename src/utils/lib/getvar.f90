!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

integer function RAMS_getvar(string, ngrd, a, b, flnm)

  use an_header,only:&
       anal_table,    &
       nvbtab

  implicit none

  include "files.h"
  include "i8.h"

  ! Arguments:
  character(len=*), intent(in) :: string
  integer, intent(in)          :: ngrd
  real, intent(inout)          :: a(*), b(*)
  character(len=*), intent(in) :: flnm
  ! Local variables:
  integer                      :: itype !!, rams_c_pos
  character(LEN=1)             :: cgrid
  character(LEN=f_name_length) :: flng
  character(LEN=120)           :: errmsg
  logical                      :: there
  integer                      :: ni
  integer(kind=i8)             :: npts, iword

  print*,'getvar:',string

  do ni=1,nvbtab

     if(string==anal_table(ni)%string .and. ngrd==anal_table(ni)%ngrid) then

        write(cgrid,'(i1)') ngrd
        flng=trim(flnm)//'-g'//cgrid//'.vfm'

        inquire(file=flng(1:len_trim(flng)),exist=there)
        if(.not.there) then
           errmsg='File not found - '//trim(flng)
           RAMS_getvar = 1
           call error_mess(errmsg)
           return
        endif

        npts=anal_table(ni)%nvalues
        itype=anal_table(ni)%idim_type
        iword=anal_table(ni)%npointer

        !  print*,'gv:opening'
        call RAMS_c_open(flng(1:len_trim(flng))//char(0),'r'//char(0))
        !  print*,'gv:opened'
        call vfirecr(10,a,npts,'LIN',b,iword)
        !  print*,'gv:vfirecr'
        call RAMS_c_close()
        !  print*,'gv:closed'

        RAMS_getvar=0
        print*,'getvar good:',string
        return

     endif
  enddo

  errmsg='Variable not available in this run - '//string
  call error_mess(errmsg)
  RAMS_getvar=1

  return
end function RAMS_getvar




subroutine FindFieldInAnalysisFile(string, ngrd, flnm, flng, npts, fPosition)
  use dump, only: &
   dumpMessage

  use an_header,only:&
       anal_table,    &
       nvbtab

  implicit none

  include "i8.h"
  include "constants.f90"
  ! Arguments:
  character(len=*), intent(in ) :: string
  integer,          intent(in ) :: ngrd
  character(len=*), intent(in ) :: flnm
  character(len=*), intent(out) :: flng
  integer(kind=i8), intent(out) :: npts
  integer(kind=i8), intent(out) :: fPosition

  ! Local variables:
  character(LEN=1)  :: cgrid
  logical           :: there
  integer           :: ni
  character(len=*), parameter  :: h="**(FindFieldInAnalysisFile)**"

  do ni=1,nvbtab

     if(string==anal_table(ni)%string .and. ngrd==anal_table(ni)%ngrid) then

        write(cgrid,'(i1)') ngrd
        flng=trim(flnm)//'-g'//cgrid//'.vfm'

        inquire(file=flng(1:len_trim(flng)),exist=there)
        if(.not.there) then
           !call fatal_error(h//'File not found - '//trim(flng)//&
            !    ' while searching for field '//trim(string)//&
            !    ' at grid '//cgrid)
           iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
                'File not found - '//trim(flng)//&
                ' while searching for field '//trim(string)//&
                ' at grid '//cgrid)
        else
           npts=anal_table(ni)%nvalues
           fPosition=anal_table(ni)%npointer
           return
        end if
     end if
  end do

  !call fatal_error(h//'Variable '//trim(string)//&
  !     &' not available in the analysis file')
  iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
       'Variable '//trim(string)//&
       ' not available in the analysis file')
       
end subroutine FindFieldInAnalysisFile



subroutine GetFieldInAnalysisFile(flng, npts, fPosition, a, b)
  implicit none
  include "i8.h"
  character(len=*),             intent(in   ) :: flng
  integer(kind=i8),             intent(in   ) :: npts
  integer(kind=i8),             intent(in   ) :: fPosition
  real,                         intent(inout) :: a(*)
  real,                         intent(inout) :: b(*)

  call RAMS_c_open(trim(flng)//char(0),'r'//char(0))

  call vfirecr(10,a,npts,'LIN',b,fPosition)

  call RAMS_c_close()
end subroutine GetFieldInAnalysisFile
