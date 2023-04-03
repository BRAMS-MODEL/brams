!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

SUBROUTINE isan_driver (name_name)
  use ModDateUtils
  use isan_coms,only: &
        dnref,        &
        guess1st,     &
        hybbot,       &
        hybtop,       &
        idate,        &
        ihour,        &
        iproc_dates,  &
        iproc_flag,   &
        iproc_names,  &
        imonth,       &
        innpr,        &
        ioflgisz,     &
        ioflgvar,     &
        iszstage,     &
        is_grids,     &
        ivrstage,     &
        iyear,        &
        maxisfiles,   &
        maxiy,        &
        maxix,        & 
        maxiz,        &
        maxsigz,      &
        natime,       &
        nfeedvar,     &
        nigrids,      &
        nisn,         &
        npry,         &
        nprx,         &
        nprz,         &
        npdates,      &
        nsigz,        &
        piref,        &
        pi_p,         &
        pi_r,         &
        pi_s,         &
        pi_scrb,      &
        pi_scra,      &
        pi_u,         &
        pi_v,         &
        ps_p,         &
        ps_r,         &
        ps_scrb,      &
        ps_scra,      &
        ps_t,         &
        ps_u,         &
        ps_v,         &
        p_u,          &
        rr_scr1,      &
        rr_scr2,      &
        rr_vt2da,     &
        rs_p,         &
        rs_r,         &
        rs_s,         &
        rs_t,         &
        rs_top,       &
        rs_u,         &
        rs_v,         &
        rs_qual,      &
        rs_sfp,       &
        rs_sft,       &
        rs_slp,       &
        rs_snow,      &
        rs_sst,       &
        rtref,        &
        sigz,         &
        thref,        &
        topsigz,      &
        varpfx
  use mem_grid
  use io_params

  implicit none

  include "files.h"

  character(len=*), intent(IN) :: name_name

  character(len=3) :: csuff
  character(len=f_name_length) :: locfn, locfna

!!$  ! fnames dimensioned by a max number of files possible
!!$  character(len=128) :: fnames(maxisfiles)

  integer, dimension(maxgrds) :: itoptn,iglatn,iglonn
  integer :: ifm,icm,ng,i,k,ifileok

  ! Not necessary because all information is read by Namelist

!!$  ! Read input ISAN namelists
!!$
!!$  call rams_f_open(1,name_name(1:len_trim(name_name)),'FORMATTED','OLD','READ',0)
!!$  call namein_isan(1,'$ISAN_CONTROL',1)
!!$  call namein_isan(1,'$ISAN_ISENTROPIC',1)
!!$  close(1)
  call opspec4

  !  Allocate grid data type since almost all model memory is not used for
  !     ISAN.

  print*,'start ISAN grid alloc'
  allocate(grid_g(ngrids))
  do ng=1,ngrids
     call nullify_grid(grid_g(ng))
     call alloc_grid(grid_g(ng),nnzp(ng),nnxp(ng),nnyp(ng),ng,if_adap) 
  enddo

  ! Allocate global grid variables data type.
  allocate(oneGlobalGridData(ngrids))
  do ng=1,ngrids
     call nullify_GlobalGridData(oneGlobalGridData(ng))
     call alloc_GlobalGridData(oneGlobalGridData(ng), nnxp(ng), nnyp(ng))
  end do

  ! ISAN runs all the time in sigma-z vertical coordinate. Reset ADAP flag...
  if_adap=0 

  ! Topo is read from surface files.

  do ng=1,ngrids
     call newgrid(ng)
     call TopReadStoreOwnChunk(ng)
  enddo

  !  Setup RAMS hroizontal and vertical grid structure.
  call gridSetup(3) !(2)

  !  Allocate RAMS grid arrays where data analysis will be put 
  !     for output and feedback.  Old way was to use "A"

  maxix=maxval(nnxp(1:nigrids))   
  maxiy=maxval(nnyp(1:nigrids))   
  maxiz=maxval(nnzp(1:nigrids))   

  do ngrid=1,nigrids  
     allocate(is_grids(ngrid)%rr_u (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_v (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_t (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_r (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_p (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_ug (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_vg (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_tg (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_rg (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_pg (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_pi0 (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_th0 (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_dn0 (nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_dn0u(nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_dn0v(nnzp(ngrid),nnxp(ngrid),nnyp(ngrid)))

     allocate(is_grids(ngrid)%rr_slp (nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_sfp (nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_sft (nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_snow(nnxp(ngrid),nnyp(ngrid)))
     allocate(is_grids(ngrid)%rr_sst (nnxp(ngrid),nnyp(ngrid)))
  enddo
  allocate(rr_scr1(maxix*maxiy*maxiz))
  allocate(rr_scr2(maxix*maxiy*maxiz))
  allocate(rr_vt2da(maxix*maxiy))

  ! Do inventory of input file names, determine which times to process

  call ISAN_file_inv(iyear1,imonth1,idate1,itime1,timmax)

  ! Start main processing loop for each data time

  do natime=1,npdates
     if (iproc_flag(natime,1) == 0) cycle
     call date_unmake_big(iyear,imonth,idate,ihour,iproc_dates(natime))
     print*
     print*,'========================================================'
     print*,'ISAN processing time: ',natime,' ',iyear,imonth,idate,ihour
     print*
     ihour=ihour/100

     innpr=iproc_names(natime,1)

     if(iproc_flag(natime,2).eq.1.and.guess1st.eq.'PRESS') then

        if(nhemgrd2 == 0) then
           call pressure_stage(nnxp(1),nnyp(1),nhemgrd2  &
                ,grid_g(1)%glat(1,1),grid_g(1)%glon(1,1)  &
                ,grid_g(1)%glat(1,1),grid_g(1)%glon(1,1))
        else
           call pressure_stage(nnxp(1),nnyp(1),nhemgrd2  &
                ,grid_g(1)%glat(1,1),grid_g(1)%glon(1,1)  &
                ,grid_g(nhemgrd2)%glat(1,1)  &
                ,grid_g(nhemgrd2)%glon(1,1))
        endif
        print*,'after pressure stage',nprx,npry,p_u(nprx/2,npry/2,1:nprz)

     endif

     ! Isentropic/sigma-z analysis to all requested RAMS grids

     do ngrid=1,nigrids

        ! Find number of sigma-z levels

        if(guess1st == 'RAMS') then
           topsigz=50.e3 ;   hybtop =50.e3;   hybbot =50.e3
           print*,'************FIRST GUESS RAMS*************'
           print*,'*resetting topsigz,hybbot,hybtop to 50km*'
           print*,'****************************************'
        endif
        do k=1,nnzp(ngrid)
           if(ztn(k,ngrid) > topsigz) goto 100
           sigz(k)=ztn(k,ngrid)
        enddo
        k=nnzp(ngrid)+1
100     continue
        nsigz=k-1


        !         Allocate memory for isentropic analysis
        !         --------------------------------------------------------
        print*,'Allocating RAMS polar/isentropic grid-'  &
             ,ngrid,nnxp(ngrid),nnyp(ngrid),nisn

        allocate(pi_u(nnxp(ngrid),nnyp(ngrid),nisn))
        allocate(pi_v(nnxp(ngrid),nnyp(ngrid),nisn))
        allocate(pi_p(nnxp(ngrid),nnyp(ngrid),nisn))
        allocate(pi_s(nnxp(ngrid),nnyp(ngrid),nisn))
        allocate(pi_r(nnxp(ngrid),nnyp(ngrid),nisn))
        allocate(pi_scra(nnxp(ngrid),nnyp(ngrid),nisn))
        allocate(pi_scrb(nnxp(ngrid),nnyp(ngrid),nisn))

        !         Allocate memory for sigma-z analysis 
        !         --------------------------------------------------------

        print*,'Allocating RAMS polar/sigmaz grid-'  &
             ,ngrid,nnxp(ngrid),nnyp(ngrid),nsigz

        allocate(ps_u(nnxp(ngrid),nnyp(ngrid),nsigz))
        allocate(ps_v(nnxp(ngrid),nnyp(ngrid),nsigz))
        allocate(ps_p(nnxp(ngrid),nnyp(ngrid),nsigz))
        allocate(ps_t(nnxp(ngrid),nnyp(ngrid),nsigz))
        allocate(ps_r(nnxp(ngrid),nnyp(ngrid),nsigz))
        allocate(ps_scra(nnxp(ngrid),nnyp(ngrid),nsigz))
        allocate(ps_scrb(nnxp(ngrid),nnyp(ngrid),nsigz))

        !         Allocate memory for surface analysis 
        !         --------------------------------------------------------

        print*,'Allocating RAMS polar/surface grid-'  &
             ,ngrid,nnxp(ngrid),nnyp(ngrid)

        allocate(rs_u(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_v(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_p(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_t(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_r(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_s(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_top(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_qual(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_slp(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_sfp(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_sft(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_snow(nnxp(ngrid),nnyp(ngrid)))
        allocate(rs_sst(nnxp(ngrid),nnyp(ngrid)))


        ! Do isentropic and sigma-z analysis

        if(iszstage == 1) then
           call isnstage ()

           ! Output isentropic file if desired

           if(ioflgisz == 1)then
              write(csuff,'(a1,i1)') 'g',ngrid
              call makefnam (locfn,varpfx,0,iyear,imonth,idate,  &
                   ihour*100,'I',csuff,'vfm')
              call rams_f_open(13,locfn(1:len_trim(locfn)),'FORMATTED','REPLACE','READ',iclobber)
              call isenio ('OUT',13,nnxp(ngrid),nnyp(ngrid))
              call sigzio ('OUT',13,nnxp(ngrid),nnyp(ngrid))
              close(13)
           endif

        elseif(ivrstage == 1) then
           write(csuff,'(a1,i1)') 'g',ngrid
           call makefnam(locfn,varpfx,0,iyear,imonth,idate,  &
                ihour*100,'I',csuff,'vfm')
           call rams_f_open(13,locfn(1:len_trim(locfn)),'FORMATTED','OLD','READ',iclobber)
           call isenio('IN',13,nnxp(ngrid),nnyp(ngrid))
           call sigzio('IN',13,nnxp(ngrid),nnyp(ngrid))
           close(13)
        endif

        ! Prepare "varfiles"

        if(ivrstage == 1) then
           call makevarf(ngrid)
        endif

        deallocate(pi_u,pi_v,pi_p,pi_s,pi_r,pi_scra,pi_scrb)
        deallocate(ps_u,ps_v,ps_p,ps_t,ps_r,ps_scra,ps_scrb)
        deallocate(rs_u,rs_v,rs_p,rs_t,rs_r,rs_s,rs_top,rs_qual)
        deallocate(rs_slp,rs_sfp,rs_sft,rs_snow,rs_sst)

     enddo

     ! Do the nesting feedback and write out the "varfiles"

     if(ivrstage == 1) then

        if(nfeedvar == 1 .and. nigrids > 1) then

           ! fill reference states for all grids

           call varfile_refstate(nnzp(1),nnxp(1),nnyp(1)  &
                ,is_grids(1)%rr_t(1,1,1),is_grids(1)%rr_p(1,1,1)  &
                ,is_grids(1)%rr_pi0(1,1,1),is_grids(1)%rr_th0(1,1,1)  &
                ,is_grids(1)%rr_r(1,1,1),is_grids(1)%rr_dn0(1,1,1)  &
                ,is_grids(1)%rr_dn0u(1,1,1),is_grids(1)%rr_dn0v(1,1,1)  &
                ,grid_g(1)%topt(1,1),grid_g(1)%rtgt(1,1)  &
                ,ztn(1,1),ztop,piref(1,1),thref(1,1),dnref(1,1),rtref(1,1))
           is_grids(1)%rr_p(1:nnzp(1),1:nnxp(1),1:nnyp(1)) =  &
                is_grids(1)%rr_p  (1:nnzp(1),1:nnxp(1),1:nnyp(1)) &
                -is_grids(1)%rr_pi0(1:nnzp(1),1:nnxp(1),1:nnyp(1))

           do ifm=2,nigrids
              icm=nxtnest(ifm)
              call fmrefs1d_isan(ifm,icm,maxsigz,nnzp(ifm)  &
                   ,piref,thref,dnref,rtref)
              call fmrefs3d_isan(ifm,icm,nnzp(ifm),nnxp(ifm),nnyp(ifm) &
                   ,nnzp(icm),nnxp(icm),nnyp(icm),maxiz,maxix,maxiy  &
                   ,nnstbot(ifm),nnsttop(ifm),jdim  &
                   ,rr_scr1(1),rr_scr2(1),rr_vt2da(1)  &
                   ,grid_g(ifm)%topt(1,1),grid_g(icm)%topt(1,1) &
                   ,is_grids(icm)%rr_dn0(1,1,1),is_grids(ifm)%rr_dn0(1,1,1) &
                   ,is_grids(icm)%rr_th0(1,1,1),is_grids(ifm)%rr_th0(1,1,1) &
                   ,is_grids(ifm)%rr_pi0(1,1,1),is_grids(ifm)%rr_dn0u(1,1,1) &
                   ,is_grids(ifm)%rr_dn0v(1,1,1),ztn(1,ifm),ztop )

              call fmdn0_isan(ifm,icm,nnzp(ifm),nnxp(ifm),nnyp(ifm) &
                   ,nnzp(icm),nnxp(icm),nnyp(icm),maxiz,maxix,maxiy &
                   ,rr_scr1(1),rr_scr2(1)  &
                   ,grid_g(ifm)%topt(1,1),grid_g(icm)%topt(1,1) &
                   ,is_grids(ifm)%rr_dn0(1,1,1),is_grids(ifm)%rr_dn0u(1,1,1) &
                   ,is_grids(ifm)%rr_dn0v(1,1,1),ztn(1,ifm),ztop )

              is_grids(ifm)%rr_p(1:nnzp(ifm),1:nnxp(ifm),1:nnyp(ifm)) =  &
                   is_grids(ifm)%rr_p  (1:nnzp(ifm),1:nnxp(ifm),1:nnyp(ifm)) &
                   -is_grids(ifm)%rr_pi0(1:nnzp(ifm),1:nnxp(ifm),1:nnyp(ifm))

           enddo

           ! Feed back u,v,pi,t,and rt fields

           do ifm=nigrids,2,-1
              icm=nxtnest(ifm)
              if (icm == 0) cycle
              call varfile_nstfeed(ifm,icm,nnzp(ifm),nnxp(ifm),nnyp(ifm) &
                   ,nnzp(icm),nnxp(icm),nnyp(icm)  &
                   ,nnstbot(icm),nnsttop(icm))
           enddo

           do ifm=1,nigrids
              is_grids(ifm)%rr_p(1:nnzp(ifm),1:nnxp(ifm),1:nnyp(ifm)) =  &
                   is_grids(ifm)%rr_p  (1:nnzp(ifm),1:nnxp(ifm),1:nnyp(ifm)) &
                   +is_grids(ifm)%rr_pi0(1:nnzp(ifm),1:nnxp(ifm),1:nnyp(ifm))
           enddo

        endif


        if(ioflgvar == 1) then
           do ng=1,nigrids
              nxyzp=nnxp(ng)*nnyp(ng)*nnzp(ng)
              nxyp =nnxp(ng)*nnyp(ng)
              write(csuff,'(a1,i1)') 'g',ng
              call makefnam (locfn,varpfx,0,iyear,imonth,idate  &
                   ,ihour*100,'V',csuff,'vfm')
              call rams_f_open (2,locfn(1:len_trim(locfn)),'FORMATTED','REPLACE','WRITE',iclobber)
              write(2,11) 999999,2
11            format(i7,i3)
              write(2,10) iyear,imonth,idate,ihour  &
                   ,nnxp(ng),nnyp(ng),nnzp(ng)  &
                   ,platn(ng),plonn(ng),deltaxn(ng)  &
                   ,deltayn(ng),deltazn(ng),dzrat,dzmax
10            format(7i5,7f14.5)
              call vforec(2,is_grids(ng)%rr_u(1,1,1),nxyzp,18  &
                   ,rr_scr1(1),'LIN')
              call vforec(2,is_grids(ng)%rr_v(1,1,1),nxyzp,18  &
                   ,rr_scr1(1),'LIN')
              call vforec(2,is_grids(ng)%rr_p(1,1,1),nxyzp,18  &
                   ,rr_scr1(1),'LIN')
              call vforec(2,is_grids(ng)%rr_t(1,1,1),nxyzp,18  &
                   ,rr_scr1(1),'LIN')
              call vforec(2,is_grids(ng)%rr_r(1,1,1),nxyzp,18  &
                   ,rr_scr1(1),'LIN')

              call vmissw(is_grids(ng)%rr_slp(1,1),nxyp,rr_vt2da(1),1E30,-1.)
              call vforec(2,rr_vt2da(1),nxyp,18,rr_scr1(1),'LIN')
              call vmissw(is_grids(ng)%rr_sfp(1,1),nxyp,rr_vt2da(1),1E30,-1.)
              call vforec(2,rr_vt2da(1),nxyp,18,rr_scr1(1),'LIN')
              call vmissw(is_grids(ng)%rr_sft(1,1),nxyp,rr_vt2da(1),1E30,-1.)
              call vforec(2,rr_vt2da(1),nxyp,18,rr_scr1(1),'LIN')
              call vmissw(is_grids(ng)%rr_snow(1,1),nxyp,rr_vt2da(1),1E30,-1.)
              call vforec(2,rr_vt2da(1),nxyp,18,rr_scr1(1),'LIN')
              call vmissw(is_grids(ng)%rr_sst(1,1),nxyp,rr_vt2da(1),1E30,-1.)
              call vforec(2,rr_vt2da(1),nxyp,18,rr_scr1(1),'LIN')
              close(2)
           enddo
           call makefnam (locfn,varpfx,0,iyear,imonth,idate  &
                ,ihour*100,'V','$','tag')
           call rams_f_open (2,locfn(1:len_trim(locfn)),'FORMATTED','REPLACE','WRITE',iclobber)
           write(2,*) nigrids
           close(2)

        endif

     endif

  enddo

  return
end SUBROUTINE isan_driver

!******************************************************************************

SUBROUTINE OPSPEC4

  use mem_grid
  use isan_coms

  implicit none

  integer :: ifaterr,iwarerr,infoerr

  ! This routine checks the option specifications in the $MODEL_GRIDS
  !    $ISAN_ISENTROPIC namelists for consistency.

  IFATERR=0
  IWARERR=0
  INFOERR=0

  ! Don't allow NIGRIDS <= NGRIDS

  if(nigrids.gt.ngrids)then
     print*,' FATAL - NIGRIDS must be <= NGRIDS'
     IFATERR=IFATERR+1
  endif

  ! Make sure that TIMMAX >= ISAN_INC

  !fisan_inc=int(ISAN_INC/100)*3600+float(mod(ISAN_INC,100))
  !if(fisan_inc.gt.TIMMAX)then
  !  print*,' WARNING - TIMMAX must be <= ISAN_INC'
  !  print*,'           resetting to TIMMAX to (s) ',fisan_inc
  !  timmax=fisan_inc
  !  timstr=fisan_inc
  !  IWARERR=IWARERR+1
  !endif

  ! Stop the run if there are any fatal errors.  List how many
  !   warning and informative errors.

  PRINT*,' -----------opspec4--------------------------'
  PRINT*,' FATAL     errors - ',IFATERR
  PRINT*,' WARNING   errors - ',IWARERR
  PRINT*,' INFORM  messages - ',INFOERR
  PRINT*,' -----------------------------------------------'

  IF(IFATERR.GT.0) STOP 'OPSPEC4'

  RETURN
END SUBROUTINE OPSPEC4
