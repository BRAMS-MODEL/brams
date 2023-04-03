!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990,1995,1999,2000,2002,2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


subroutine nud_read(initflag)
use ModDateUtils
use mem_grid
use mem_varinit

!--(DMK-CCATT-INI)-----------------------------------------------------
use mem_chem1, only: &
     chem_assim
!--(DMK-CCATT-FIM)-----------------------------------------------------

implicit none

integer :: initflag

character(len=14)  :: itotdate_start
integer :: iyears,imonths,idates,ihours,nf,ifm

if (initflag == 1) then   ! Initialization

   ! Inventory all history files. 

   call nud_file_inv (nud_hfile,iyear1,imonth1,idate1,itime1)


   ! Find past time nudging file

   if (runtype == 'HISTORY') then
      call date_add_to(iyear1,imonth1,idate1,itime1*100  &
                   ,time,'s',iyears,imonths,idates,ihours)
      call date_make_big(iyears,imonths,idates,ihours  &
                   ,itotdate_start)
   elseif (runtype == 'INITIAL') then
      call date_make_big(iyear1,imonth1,idate1,itime1*100  &
                   ,itotdate_start)
   endif

   do nf=1,nnudfiles
      if(itotdate_start >= itotdate_nud(nf) .and.  &
            itotdate_start <  itotdate_nud(nf+1) ) then
         nnudfl=nf
         exit
      endif
   enddo

   print*,'nud starting at history file:',nnudfl

   ! Calculate varweights just like var init
   print *,'LFR-DEB->nud_read.f90'
   call varweight(nnzp(1),nnxp(1),nnyp(1),varinit_g(1)%varwts(1,1,1)  &
       ,grid_g(1)%topt(1,1),grid_g(1)%rtgt(1,1))

!--(DMK-CCATT-INI)---------------------------------------------------------
   if(chem_assim == 1) &
     call varweight_chem(nnzp(1),nnxp(1),nnyp(1),varinit_g(1)%varwts_chem(1,1,1)  &
                   ,grid_g(1)%topt(1,1),grid_g(1)%rtgt(1,1))
!--(DMK-CCATT-FIM)---------------------------------------------------------
   
   ! Read and interpolate files to new grid 1. Put stuff in varinit arrays.
   
   call nud_update(0,nnudfl)
   
   do ifm = 2,ngrids     
      call newgrid(ifm)
      call vfintrpf(ifm,1)
   enddo

elseif (initflag == 2 .or. initflag == 4) then   ! Runtime file increment

   if ( time >= nud_times(nnudfl+1) ) nnudfl = nnudfl + 1

endif


! Read and interpolate files to new grid 1.

call nud_update(1,nnudfl+1)

! Fill nested grid nudging arrays for grids we didn't fill directly

!do ifm = 2,ngrids
!   if(igrid_match(ifm) == -1) then
!      call newgrid(ifm)
!      call vfintrpf(ifm,2)
!   endif
!enddo

htime1=nud_times(nnudfl)
htime2=nud_times(nnudfl+1)

return
end



subroutine nud_file_inv (hfilin,iyear1,imonth1,idate1,itime1)
use ModDateUtils
use isan_coms, only: &
       ISAN_INC
use mem_varinit, only: &
    fnames_nud,        &
    itotdate_nud,      &
    maxnudfiles,       &
    nnudfiles,         &
    nud_times
use mem_grid, only: &
    time, &
    timmax 

implicit none

include "files.h"

character(len=f_name_length) :: hfilin
integer :: iyear1,imonth1,idate1,itime1
logical there 
integer :: localTime
character(len=f_name_length) :: sVarName 
integer :: iyears,imonths,idates,ihours
integer :: indice

integer :: nc,nf,lnf,nhftot
integer :: inyear,inmonth,indate,inhour


character(len=f_name_length), dimension(maxnudfiles) :: fnames
character(len=f_name_length) :: hpref
character(len=f_name_length) :: rams_filelist_arg
character(len=14)  :: itotdate
real(kind=8) :: secs_init,secs_nud

! Get abs seconds of run start

call date_abs_secs2(iyear1,imonth1,idate1,itime1*100,secs_init)

! Go through history files and make inventory

nc=len_trim(hfilin)
nhftot=-1
hpref=hfilin(1:nc-26)

!!$ rams_filelist_arg = trim(hpref)//'????-??-??-??????-head.txt'

!!$ call RAMS_filelist(fnames, rams_filelist_arg, nhftot)

nhftot = ((timmax/3600) / (isan_inc/100)) + 1    
           
call date_add_to(iyear1,imonth1,idate1,itime1*100  &
                ,time,'s',iyears,imonths,idates,ihours)

localTime = itime1
indice = 1  

do nf=1,nhftot
   call date_add_to(iyears, imonths, idates, localTime*100,  &
                    0., 's', iyears, imonths, idates, ihours)

   call makefnam(sVarName,hpref,time,iyears,imonths,idates,localTime*100  &
          ,'R','head','txt')

   inquire(file=sVarName(1:len_trim(sVarName)),exist=there)

   if (there) then
      fnames(indice) = trim(sVarName)//"-head.txt"
      indice = indice + 1
   endif

!--(DMK-CCATT-INI)---------------------------------------------------------
   if(localTime .LE. 2359)then
!--(DMK-CCATT-OLD)---------------------------------------------------------
!   if(localTime .LE. 1800)then
!--(DMK-CCATT-FIM)---------------------------------------------------------

      localTime = localTime + isan_inc
   else
      localTime = 000000
      localTime = localTime + isan_inc
   endif 
enddo

nhftot = indice - 1

if(nhftot > maxnudfiles) then
   print*,'too many nud history files'
   stop 'lots_of_nud_history'
endif

nnudfiles=0
do nf=1,nhftot
   lnf=len_trim(fnames(nf))
   read(fnames(nf)(lnf-25:lnf-9),20) inyear,inmonth,indate,inhour
   20 format(i4,1x,i2,1x,i2,1x,i6)

   call date_make_big(inyear,inmonth,indate,inhour,itotdate)

   nnudfiles=nnudfiles+1
   fnames_nud(nnudfiles)=fnames(nf)
   itotdate_nud(nnudfiles)=itotdate
   
   call date_abs_secs2(inyear,inmonth,indate,inhour,secs_nud)
   nud_times(nnudfiles)=secs_nud - secs_init

enddo

call RAMS_dintsort(nnudfiles,itotdate_nud,fnames_nud)

!  start printing section
!--------------------------------------------------------------

print*,' '
print*,' '
print*,' '
print*,'-------------------------------------------------------------'
print*,'-----------  History Nudging Input File Inventory -----------'
print*,'-------------------------------------------------------------'
do nf=1,nnudfiles
   print*,  itotdate_nud(nf),'   ',nud_times(nf)  &
           ,trim(fnames_nud(nf))
enddo
print*,'------------------------------------------------------'

!--------------------------------------------------------------


return
end


