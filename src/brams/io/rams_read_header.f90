!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine rams_read_header(flnm)

use an_header,only: &
    anal_table,     &
    nvbtab

implicit none

include "files.h"

!--(DMK-CCATT-INI)-----------------------------------------------------
character(len=*) :: flnm
!--(DMK-CCATT-OLD)-----------------------------------------------------
!character*(f_name_length) flnm
!--(DMK-CCATT-FIM)-----------------------------------------------------

character(len=f_name_length) :: flnm2
integer lenf,nv


if(allocated(anal_table)) deallocate(anal_table)

! open analysis file and read in commons

!--(DMK-CCATT-INI)-----------------------------------------------------
open(10,file=flnm(1:len_trim(flnm))//'-head.txt',status='old')
!--(DMK-CCATT-OLD)-----------------------------------------------------
!flnm2=flnm(1:len_trim(flnm))//'-head.txt'
!open(10,file=flnm2(1:len_trim(flnm2)),status='old')
!--(DMK-CCATT-FIM)-----------------------------------------------------

read(10,*) nvbtab
allocate (anal_table(nvbtab))
do nv=1,nvbtab
   read(10,*)  anal_table(nv)%string   &
              ,anal_table(nv)%npointer  &
              ,anal_table(nv)%idim_type  &
              ,anal_table(nv)%ngrid  &
              ,anal_table(nv)%nvalues
enddo

!call commio('ANAL','READ',10)
close(10)


return
end
