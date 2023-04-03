!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################




!***************************************************************************

subroutine RAMS_dintsort(ni, chnums, cstr)
  implicit none

  include "files.h"

  integer, intent(IN)                         :: ni
!!$character(len=14) :: chnums(*)
!!$character(len=*) :: cstr(*)
  character(len=14), intent(INOUT)            :: chnums(ni)
  character(len=f_name_length), intent(INOUT) :: cstr(ni)
  
  ! sort an array of character strings by an associated character field
  
  character(len=f_name_length) :: cscr
  character(len=14)            :: mini, nscr
  integer :: n, nm, nmm
  
  do n=1,ni
     mini='99999999999999'
     do nm=n,ni
        if(chnums(nm).lt.mini) then
           nmm=nm
           mini=chnums(nm)
        endif
     enddo
     nscr=chnums(n)
     chnums(n)=chnums(nmm)
     chnums(nmm)=nscr
     cscr=cstr(n)
     cstr(n)=cstr(nmm)
     cstr(nmm)=cscr
  enddo
  
end subroutine RAMS_dintsort

!***************************************************************************

subroutine RAMS_sort_dint3 (n1,ia1,n2,ia2,n3,ia3,nt,iall)
implicit none
integer :: n1,n2,n3,nt
character(len=14) :: ia1(*),ia2(*),ia3(*),iall(*)

!     sort 3 arrays of char's, put back in 1 array
!     copy all to output array

character(len=14) :: mini,nscr
integer :: n,nm,nmm

nt=0
do n=1,n1
   nt=nt+1
   iall(nt)=ia1(n)
enddo
do n=1,n2
   nt=nt+1
   iall(nt)=ia2(n)
enddo
do n=1,n3
   nt=nt+1
   iall(nt)=ia3(n)
enddo

do n=1,nt
   mini='99999999999999'
   do nm=n,nt
      if(iall(nm).lt.mini) then
         nmm=nm
         mini=iall(nm)
      endif
   enddo
   nscr=iall(n)
   iall(n)=iall(nmm)
   iall(nmm)=nscr
enddo

return
end

!***************************************************************************

subroutine RAMS_unique_dint (n1,ia1)
implicit none
integer :: n1
character(len=14) :: ia1(*)

integer :: n,nt,nn

! reduce an array to get rid of duplicate entries


nt=n1
10 continue
do n=2,nt
   if(ia1(n).eq.ia1(n-1)) then
      do nn=n,nt
         ia1(nn-1)=ia1(nn)
      enddo
      nt=nt-1
      goto 10
   endif
enddo
n1=nt

return
end

!***************************************************************************

!--(DMK-CCATT-INI)----------------------------------------------------------------
INTEGER FUNCTION julday (imonth,iday,iyear)
  IMPLICIT NONE
  INTEGER :: imonth,iday,iyear

  ! compute the julian day from a normal date

  julday= iday  &
       + MIN(1,MAX(0,imonth-1))*31  &
       + MIN(1,MAX(0,imonth-2))*(28+(1-MIN(1,MOD(iyear,4))))  &
       + MIN(1,MAX(0,imonth-3))*31  &
       + MIN(1,MAX(0,imonth-4))*30  &
       + MIN(1,MAX(0,imonth-5))*31  &
       + MIN(1,MAX(0,imonth-6))*30  &
       + MIN(1,MAX(0,imonth-7))*31  &
       + MIN(1,MAX(0,imonth-8))*31  &
       + MIN(1,MAX(0,imonth-9))*30  &
       + MIN(1,MAX(0,imonth-10))*31  &
       + MIN(1,MAX(0,imonth-11))*30  &
       + MIN(1,MAX(0,imonth-12))*31

  RETURN
END FUNCTION julday
!--(DMK-CCATT-FIM)----------------------------------------------------------------
