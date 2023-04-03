!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine oda_read()
  use ModDateUtils
  use mem_oda,only:  &
       fnames_sfc,    &  
       fnames_upa,    &
       itotdate_sfc,  &
       itotdate_upa,  &
       maxodafiles,   &
       nsfcfiles,     &
       nupafiles,     &
       oda_sfcprefix, &
       oda_upaprefix
  use mem_grid

  implicit none

  integer :: iyears,imonths,idates,ihours,iyearf,imonthf,idatef,ihourf

  ! Inventory all observation files. Due to possible mismatches with 
  !   start/end time and obs file times, we are going to start 12 hours
  !   before the actual start and end 12 hours after timmax

  call date_add_to(iyear1,imonth1,idate1,itime1*100  &
       ,-43200.,'s',iyears,imonths,idates,ihours)
  call date_add_to(iyear1,imonth1,idate1,itime1*100  &
       ,timmax+43200.,'s',iyearf,imonthf,idatef,ihourf)

  call oda_file_inv (iyears,imonths,idates,ihours  &
       ,iyearf,imonthf,idatef,ihourf)
  print*,'++++++++++++++++: ',nupafiles

  ! First pass through the files: find number of unique station ids

  call oda_sta_count(platn(1),plonn(1),ngrids)
  print*,'++++++++++++++++: ',nupafiles

  ! Allocate obs structures

  call oda_obs_alloc()

  ! Fill data structures

  call oda_sta_input(platn(1),plonn(1),ngrids)

  ! Do this later!!!!!!!!!
!!!!!!     ! Sort obs according to time
!!! call sort_me()

!!! stop

  return
end subroutine oda_read

!--------------------------------------------------------------------

subroutine oda_obs_alloc ()

  use mem_oda

  implicit none

  integer :: ns

  allocate (oda_sfc_info(num_oda_sfc), oda_sfc_obs(num_oda_sfc))
  allocate (oda_upa_info(num_oda_upa), oda_upa_obs(num_oda_upa))

  !  Assuming one data time per ralph file for the allocations
  !    This wasn't enough!!!!! Actual count of times have been done

  do ns=1,num_oda_sfc
     allocate(oda_sfc_obs(ns)%time (maxtimes_sfc))
     allocate(oda_sfc_obs(ns)%temp (maxtimes_sfc))
     allocate(oda_sfc_obs(ns)%dewpt(maxtimes_sfc))
     allocate(oda_sfc_obs(ns)%us(maxtimes_sfc))
     allocate(oda_sfc_obs(ns)%vs(maxtimes_sfc))
     allocate(oda_sfc_obs(ns)%u (maxtimes_sfc))
     allocate(oda_sfc_obs(ns)%v (maxtimes_sfc))
     allocate(oda_sfc_obs(ns)%ps(maxtimes_sfc))
  enddo

  ! Upper air
  do ns=1,num_oda_upa
     allocate(oda_upa_obs(ns)%time(maxtimes_upa))
     allocate(oda_upa_obs(ns)%lp  (maxtimes_upa))
     allocate(oda_upa_obs(ns)%lz  (maxtimes_upa))
     allocate(oda_upa_obs(ns)%theta(maxupalevs,maxtimes_upa))
     allocate(oda_upa_obs(ns)%rv(maxupalevs,maxtimes_upa))
     allocate(oda_upa_obs(ns)%us(maxupalevs,maxtimes_upa))
     allocate(oda_upa_obs(ns)%vs(maxupalevs,maxtimes_upa))
     allocate(oda_upa_obs(ns)%zz(maxupalevs,maxtimes_upa))
     allocate(oda_upa_obs(ns)%u (maxupalevs,maxtimes_upa))
     allocate(oda_upa_obs(ns)%v (maxupalevs,maxtimes_upa))
     allocate(oda_upa_obs(ns)%pi(maxupalevs,maxtimes_upa))
     allocate(oda_upa_obs(ns)%zgeo (maxupalevs,maxtimes_upa))
  enddo

  return
end subroutine oda_obs_alloc

!--------------------------------------------------------------------

subroutine oda_file_inv (iyear1,imonth1,idate1,itime1  &
     ,iyear2,imonth2,idate2,itime2)
  use ModDateUtils
  use mem_oda, only: &
       ODA_UPAPREFIX, &
       NUPAFILES, &
       FNAMES_UPA, &
       ITOTDATE_UPA, &
       ODA_SFCPREFIX, &
       NSFCFILES, &
       FNAMES_SFC, &
       ITOTDATE_SFC, &
       MAXODAFILES
  use mem_grid, only: &
       time, &
       timmax

  use isan_coms, only: &
       ISAN_INC


  implicit none

  include "files.h"

  integer :: iyear1,imonth1,idate1,itime1,iyear2,imonth2,idate2,itime2


  integer :: nc,nf,lnf,nsfctot,nupatot
  integer :: inyear,inmonth,indate,inhour

  logical there 
  integer :: localTime
  character(len=f_name_length) :: sVarName 
  integer :: iyears,imonths,idates,ihours
  integer :: indice
  character(len=4) :: cdummy(2)
  character(len=f_name_length), dimension(maxodafiles) :: fnames
  character(len=f_name_length) :: rams_filelist_arg
  character(len=14)  :: itotdate,itotdate_start,itotdate_end

  !          Go through upper air and surface input files
  !            and make inventory

  print*,'st:',iyear1,imonth1,idate1,itime1
  print*,'en:',iyear2,imonth2,idate2,itime2
  call date_make_big(iyear1,imonth1,idate1,itime1,itotdate_start)
  call date_make_big(iyear2,imonth2,idate2,itime2,itotdate_end)

  if(oda_upaprefix(1:1) /= ' ' .and. oda_upaprefix(1:1) /= char(0) ) then

     nc=len_trim(oda_upaprefix)
     nupatot=-1

!!$ rams_filelist_arg = oda_upaprefix(1:nc)//'????-??-??-????'

!!$ call RAMS_filelist(fnames, rams_filelist_arg, nupatot)

     nupatot = ((timmax/3600) / (isan_inc/100)) + 1    

     call date_add_to(iyear1,imonth1,idate1,itime1*100  &
          ,time,'s',iyears,imonths,idates,ihours)

     localTime = itime1
     indice = 1  

     do nf=1,nupatot

        call date_add_to(iyears, imonths, idates, localTime*100,  &
             0., 's', iyears, imonths, idates, ihours)

        write(sVarName,100) oda_upaprefix(1:len_trim(oda_upaprefix)),iyears,'-',imonths,'-',idates,'-',ihours/100   
100     format(a,i4.4,a1,i2.2,a1,i2.2,a1,i4.4)

        inquire(file=sVarName(1:len_trim(sVarName)),exist=there)

        if (there) then
           fnames(indice) = trim(sVarName)
           indice = indice + 1
        endif

        if(localTime .LE. 1800)then
           localTime = localTime + isan_inc
        else
           localTime = 000000
           localTime = localTime + isan_inc
        endif
     enddo

     nupatot = indice - 1

     if(nupatot > maxodafiles) then
        print*,'too many oda upper air files'
        stop 'lots_of_oda_upper_air'
     endif

     nupafiles=0
     do nf=1,nupatot
        lnf=len_trim(fnames(nf))
        !   print*,lnf ,nupatot,fnames(nf)
        !   print*,lnf ,fnames(nf)(lnf-14:lnf)
        read(fnames(nf)(lnf-14:lnf),20) inyear,inmonth,indate,inhour
20      format(i4,1x,i2,1x,i2,1x,i4)

        call date_make_big(inyear,inmonth,indate,inhour*100,itotdate)
        ! print*, inyear,inmonth,indate,inhour
        !print*,itotdate,itotdate_start,itotdate_end
        if(itotdate >= itotdate_start .and. itotdate <= itotdate_end) then
           nupafiles=nupafiles+1
           fnames_upa(nupafiles)=fnames(nf)
           itotdate_upa(nupafiles)=itotdate
        endif

     enddo

     call RAMS_dintsort(nupafiles,itotdate_upa,fnames_upa)

     ! do nf=1,nupafiles
     !    print*,'up files:',nf,itotdate_upa(nf),fnames_upa(nf)
     ! enddo

  endif


  if(oda_sfcprefix(1:1) /= ' '.and. oda_sfcprefix(1:1) /= char(0) ) then

     nc=len_trim(oda_sfcprefix)
     nsfctot=-1

!!$    rams_filelist_arg = oda_sfcprefix(1:nc)//'????-??-??-????'

!!$    call RAMS_filelist(fnames, rams_filelist_arg, nsfctot)

     nsfctot = ((timmax/3600) / (isan_inc/100)) + 1    

     call date_add_to(iyear1,imonth1,idate1,itime1*100  &
          ,time,'s',iyears,imonths,idates,ihours)

     localTime = itime1
     indice = 1  

     do nf=1,nsfctot

        call date_add_to(iyears, imonths, idates, localTime*100,  &
             0., 's', iyears, imonths, idates, ihours)

        write(sVarName,200) oda_sfcprefix(1:len_trim(oda_sfcprefix)),iyears,'-',imonths,'-',idates,'-',ihours/100   
200     format(a,i4.4,a1,i2.2,a1,i2.2,a1,i4.4)

        inquire(file=sVarName(1:len_trim(sVarName)),exist=there)

        if (there) then
           fnames(indice) = trim(sVarName)
           indice = indice + 1
        endif

        if(localTime .LE. 1800)then
           localTime = localTime + isan_inc
        else
           localTime = 000000
           localTime = localTime + isan_inc
        endif
     enddo

     nsfctot = indice - 1

     if(nsfctot > maxodafiles) then
        print*,'too many oda surface air files'
        stop 'lots_of_oda_surface'
     endif

     nsfcfiles=0
     do nf=1,nsfctot
        lnf=len_trim(fnames(nf))
        !  print*,lnf ,nsfctot,fnames(nf)
        !  print*,lnf ,fnames(nf)(lnf-14:lnf)
        read(fnames(nf)(lnf-14:lnf),20) inyear,inmonth,indate,inhour

        call date_make_big(inyear,inmonth,indate,inhour*100,itotdate)

        if(itotdate >= itotdate_start .and. itotdate <= itotdate_end) then
           nsfcfiles=nsfcfiles+1
           fnames_sfc(nsfcfiles)=fnames(nf)
           itotdate_sfc(nsfcfiles)=itotdate
        endif

     enddo

     call RAMS_dintsort(nsfcfiles,itotdate_sfc,fnames_sfc)

     ! do nf=1,nsfcfiles
     !    print*,'sf files:',nf,itotdate_sfc(nf),fnames_sfc(nf)
     ! enddo

  endif


  !  start printing section
  !--------------------------------------------------------------

  print*,' '
  print*,' '
  print*,' '
  print*,'-------------------------------------------------------------'
  print*,'-----------  Obs 4DDA Input File Date Inventory -------------'
  print*,'-------------------------------------------------------------'
  print*,'---- Upper air   files:'
  do nf=1,nupafiles
     print*,  itotdate_upa(nf),'   ',trim(fnames_upa(nf))
  enddo
  print*,'---- Surface obs files:'
  do nf=1,nsfcfiles
     print*,  itotdate_sfc(nf),'   ',trim(fnames_sfc(nf))
  enddo
  print*,'------------------------------------------------------'

  !--------------------------------------------------------------

  return
end subroutine oda_file_inv
