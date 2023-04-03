!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine opspec1(MPISize, mchnum, master_num)

  ! this routine checks the option specifications in the $model_grids
  !   namelist for consistency, and overrides settings of icloud,
  !   irain, ipris, isnow, iaggr, igraup, and ihail, setting them
  !   all to zero if level is less than 3.

  use grid_dims, only: &
       maxgrds,        &
       nxpmax,         &
       nypmax,         &
       nzpmax,         &
       nzgmax

  use mem_grid, only: &
       runtype,       &
       ngrids,        &
       ngrid,         &
       nnxp,          &
       nnyp,          &
       nnzp,          &
       nzg,           &
       nzs,           &
       nstratx,       &
       nstraty,       &
       nxtnest,       &
       ninest,        &
       njnest,        &
       polelat,       &
       polelon,       &
       centlat,       &
       centlon,       &
       nhemgrd2,      &
       ibnd,          &
       jbnd,          &
       nnsttop,       &
       nnstbot,       &
       naddsc,        &
       dyncore_flag,  &
       pd_or_mnt_constraint, &
       order_h, &
       order_v


  use micphys, only: &
       mcphys_type,  &
       level,        &
       icloud,       &
       idriz,        &
       irain,        &
       ipris,        &
       isnow,        &
       iaggr,        &
       igraup,       &
       ihail,        &
       irime,iplaws,idust,isalt,imbudget,imbudtot,        &
       iccnlev,        &
       imd1flg,        &
       imd2flg,        &
       epsil


!--(DMK-CCATT-INI)------------------------------------------------------------
  use ccatt_start, only: &
       ccatt
  use mem_chem1, only: &
       chemistry
!--(DMK-CCATT-FIM)------------------------------------------------------------

  implicit none
  ! Arguments:
  integer, intent(in) :: MPISize
  integer, intent(in) :: mchnum
  integer, intent(in) :: master_num
  ! local Variables:
  integer :: ierr
  integer :: ifm
  integer :: icm
  integer :: ng
  integer :: ifaterr
  integer :: lev4bins
  integer :: lev5bins
  integer :: nhemgrds
  logical :: twod
  character(len=10) :: c0, c1, c2
  character(len=*), parameter :: h="**(opspec1)**"

  ifaterr=0

!--(DMK-CCATT-INI)------------------------------------------------------------
  if(ccatt == 0 .and. chemistry > -1) then
     print *, "FATAL - It is not allowed running chemistry greater than -1 without CCATT"
     ifaterr=ifaterr+1
  end if
!--(DMK-CCATT-FIM)------------------------------------------------------------

  ! Only serial run is allowed for MAKEVFILE (ISAN) or MAKESFC phase
  if (runtype(1:7)=='MAKESFC' .or. runtype(1:9)=='MAKEVFILE') then
     if (MPIsize>1) then
        if (mchnum==master_num) then
!!$           call fatal_error("Only serial run is allowed for "//&
!!$                trim(runtype)//" phase")
           print *, "*** For "//trim(runtype)//&
                " phase parallel run is not necessary! ***"
        endif
     endif
  endif

  ! check if number of grids is within bounds;
  ! if not, stop analyzing input data since it may be corrupted

  if (ngrids < 1 .or. ngrids > maxgrds) then
     write(c0,"(i10)") ngrids
     write(c1,"(i10)") maxgrds
     call fatal_error(h//" desired number of grids (ngrids="//&
          &trim(adjustl(c0))//") is outside allowed interval [1:"//&
          &trim(adjustl(c1))//"]")
  endif

  ! flag 2d simulation

  twod = all(nnyp(1:ngrids) == 1)

  ! check if number of grid points is within bounds

  do ngrid=1,ngrids
     if (nnxp(ngrid) < 4 .or. nnxp(ngrid) > nxpmax) then
        write(c0,"(i10)") ngrid
        write(c1,"(i10)") nnxp(ngrid)
        write(c2,"(i10)") nxpmax
        write(*,"(a)") h//" at grid "//trim(adjustl(c0))//&
             &" desired number of x intervals (nnxp="//trim(adjustl(c1))//&
             &") is outside allowed interval [4:"//trim(adjustl(c2))//"]"
        ifaterr = ifaterr + 1
     end if
     if (.not. twod .and. (nnyp(ngrid) < 4 .or. nnyp(ngrid) > nypmax)) then
        write(c0,"(i10)") ngrid
        write(c1,"(i10)") nnyp(ngrid)
        write(c2,"(i10)") nypmax
        write(*,"(a)") h//" at grid "//trim(adjustl(c0))//&
             &" desired number of y intervals (nnyp="//trim(adjustl(c1))//&
             &") is outside allowed interval [4:"//trim(adjustl(c2))//"]"
        ifaterr = ifaterr + 1
     end if
     if (nnzp(ngrid) < 11 .or. nnzp(ngrid) > nzpmax) then
        write(c0,"(i10)") ngrid
        write(c1,"(i10)") nnzp(ngrid)
        write(c2,"(i10)") nzpmax
        write(*,"(a)") h//" at grid "//trim(adjustl(c0))//&
             &" desired number of z intervals (nnzp="//trim(adjustl(c1))//&
             &") is outside allowed interval [11:"//trim(adjustl(c2))//"]"
        ifaterr = ifaterr + 1
     end if
  enddo

  ! check if number of soil levels is bellow maximum
  ! **(jp)**:  what is the minumum number of soil levels? i imposed 0.

  if (nzg < 0 .or. nzg > nzgmax) then
     write(c0,"(i10)") nzg
     write(c1,"(i10)") nzgmax
     write(*,"(a)") h//" desired number of soil levels "//&
          &"(nzg="//trim(adjustl(c0))//&
          &") is outside allowed interval [0:"//trim(adjustl(c1))//"]"
     ifaterr = ifaterr + 1
  end if

  ! atmospheric levels must exceed soil levels

  do ngrid=1,ngrids
     if (nnzp(ngrid) <= nzg+nzs) then
        write(c0,"(i10)") ngrid
        write(c1,"(i10)") nnzp(ngrid)
        write(c2,"(i10)") nzg+nzs
        write(*,"(a)") h//" at grid "//trim(adjustl(c0))//&
             &" number of z intervals (nnzp="//trim(adjustl(c1))//&
             &") should exceed number of soil levels (nzg+nzs="//&
             &trim(adjustl(c2))//")"
        ifaterr = ifaterr + 1
     endif
  enddo

  ! consistency of nstratx, nstraty

  do ifm=1,ngrids
     if (nstratx(ifm) < 1) then
        write(c0,"(i10)") ifm
        write(c1,"(i10)") nstratx(ifm)
        write(*,"(a)") h//" at grid "//trim(adjustl(c0))//&
             &" nest x ratio should be at least 1 (nstratx="//trim(adjustl(c1))//&
             &")"
        ifaterr = ifaterr + 1
     endif
     if (nstraty(ifm) < 1) then
        write(c0,"(i10)") ifm
        write(c1,"(i10)") nstraty(ifm)
        write(*,"(a)") h//" at grid "//trim(adjustl(c0))//&
             &" nest y ratio should be at least 1 (nstraty="//trim(adjustl(c1))//&
             &")"
        ifaterr = ifaterr + 1
     endif
  end do

  if (twod .and. any(nstraty(1:ngrids) /= 1)) then
     write(*,"(a)") h//" for 2d simulations, all nstraty should be 1; at least one is not"
     ifaterr = ifaterr + 1
  end if

  do ifm=1,ngrids
     icm = nxtnest(ifm)
     if (icm .ge. 1 .and. nnyp(ifm) .eq. 1 .and.  &
          (ninest(ifm) .lt. 3 .or. njnest(ifm) .lt. 1)) then
        print*, ' fatal - nested 2d grid must have ninest > 2 '  &
             ,'and njnest = 1 in namelist.'
        ifaterr=ifaterr+1
     endif
  enddo

  ! allowable values of centlat, centlon, polelat, polelon
  !   (severity - f)

  if(polelat.lt.-90..or.polelat.gt.90.) then
     print*,' fatal - polelat outside of legal bounds.'
     ifaterr=ifaterr+1
  endif

  if(polelon.lt.-180..or.polelon.gt.180.) then
     print*,' fatal - polelon outside of legal bounds.'
     ifaterr=ifaterr+1
  endif

  do ng=1,ngrids
     if(centlat(ng).lt.-90..or.centlat(ng).gt.90.) then
        print*,' fatal - centlat outside of legal bounds.'
        ifaterr=ifaterr+1
     endif

     if(centlon(ng).lt.-180..or.centlon(ng).gt.180.) then
        print*,' fatal - centlon outside of legal bounds.'
        ifaterr=ifaterr+1
     endif
  enddo

  ! check nxtnest values for validity and whether this is a global simulation
  !   (severity - f)

  if (nxtnest(1) .ne. 0) then
     print*, ' fatal - grid # 1 must have its parent mesh'  &
          ,' designated as grid # 0 (nxtnest(1) = 0)'
     ifaterr=ifaterr+1
  endif

  nhemgrds = 0
  do ifm = 1,ngrids
     icm = nxtnest(ifm)

     if (icm .ge. ifm) then
        print 1, ifm
1       format (' fatal - nest #',i3,' has specified parent'  &
             ,' mesh of equal or higher number')
        ifaterr=ifaterr+1
     endif

     if (icm .lt. 0) then
        print 2, ifm
2       format (' fatal - nest #',i3,' has specified parent'  &
             ,' mesh of number less than 0')
        ifaterr=ifaterr+1
     endif

     if (icm .eq. 0) then
        nhemgrds = nhemgrds + 1
        nhemgrd2 = ifm
     endif

  enddo

  if (nhemgrds .gt. 2) then
     print*, ' fatal - more than two grids have grid # 0 specified'  &
          ,' as their parent'
     ifaterr=ifaterr+1
  endif

  ! if this is a global simulation, print message that deltax and deltaz will be
  ! redefined by nnxp(1) and nnyp(1).  check that nnxp, nnyp, and nnzp of the
  ! two top grids are identical, and that nnxp and nnyp are equal to each
  ! other.  check to make sure that cyclic lateral boundary conditions are
  ! not specified.

  !print*, 'nhemgrd2,nhemgrds',nhemgrd2,nhemgrds
  !print*, 'nnxp(1),nnxp(nhemgrd2)',nnxp(1),nnxp(nhemgrd2)
  !print*, 'nnyp(1),nnyp(nhemgrd2)',nnyp(1),nnyp(nhemgrd2)
  !print*, 'nnzp(1),nnzp(nhemgrd2)',nnzp(1),nnzp(nhemgrd2)

  if (nhemgrds .eq. 2) then
     if (nnxp(1) .ne. nnxp(nhemgrd2) .or.  &
          nnyp(1) .ne. nnyp(nhemgrd2) .or.  &
          nnzp(1) .ne. nnzp(nhemgrd2)) then
        print*, ' fatal - for a global simulation, nnxp, nnyp, and nnzp'
        print*, ' must be identical between both hemispheric grids'
        ifaterr=ifaterr+1
     endif

     if (nnxp(1) .ne. nnyp(1)) then
        print*, ' '
        print*, 'fatal - for a global simulation, nnxp must equal'
        print*, ' nnyp in both hemispheric grids'
        ifaterr = ifaterr + 1
     endif

     if (ibnd .eq. 4 .or. jbnd .eq. 4) then
        print*, ' '
        print*, 'fatal - for a global simulation, ibnd and jbnd'
        print*, ' must not be set to 4'
        ifaterr = ifaterr + 1
     endif

     print*, ' '
     print*, 'because two values of nxtnest are set to 0, this'
     print*, ' is configured as a global simulation.'
     print*, 'consequently, ihtran will automatically be set'
     print*, ' to 1, and deltax and deltay will be redefined'
     print*, ' in terms of nnxp(1) and nnyp(1), ignoring'
     print*, ' the values specified in the namelist.'

  endif

  ! check to make sure that top grids have nsttop and nstbot set to 1
  !   (severity - f)

  if (nnsttop(1) .ne. 1 .or. nnsttop(nhemgrd2) .ne. 1) then
     print*,'nsttop not set to 1 for a hemispheric grid'
     ifaterr=ifaterr+1
  endif

  if (nnstbot(1) .ne. 1 .or. nnstbot(nhemgrd2) .ne. 1) then
     print*,'nstbot not set to 1 for a hemispheric grid'
     ifaterr=ifaterr+1
  endif

  IF(mcphys_type .le. 1) then
  ! if level is less than 3, set microphysics parameters to zero.
  ! if level equals 3, check for values of microphysics parameters
  ! that are out of bounds.  if level is equal to 4, set microphysics
  ! parameters other than icloud to zero.

    if (level .le. 2) then

     icloud = 0
     irain = 0
     ipris = 0
     isnow = 0
     iaggr = 0
     igraup = 0
     ihail = 0
     !Glauber 2014
     idriz = 0
     iccnlev = 0
     imbudget = 0
     imbudtot = 0
     !Glauber 2014

    elseif (level .eq. 3) then

     if (icloud .lt. 0 .or. icloud .gt. 7) then
        print*,'fatal - icloud out of range'
        ifaterr = ifaterr + 1
     endif
     if (irain .lt. 0 .or. irain .gt. 5) then
        print*,'fatal - irain out of range'
        ifaterr = ifaterr + 1
     endif
     if (ipris .lt. 0 .or. ipris .gt. 7) then
        print*,'fatal - ipris out of range'
        ifaterr = ifaterr + 1
     endif
     if (isnow .lt. 0 .or. isnow .gt. 5) then
        print*,'fatal - isnow out of range'
        ifaterr = ifaterr + 1
     endif
     if (iaggr .lt. 0 .or. iaggr .gt. 5) then
        print*,'fatal - iaggr out of range'
        ifaterr = ifaterr + 1
     endif
     if (igraup .lt. 0 .or. igraup .gt. 5) then
        print*,'fatal - igraup out of range'
        ifaterr = ifaterr + 1
     endif
     if (ihail .lt. 0 .or. ihail .gt. 5) then
        print*,'fatal - ihail out of range'
        ifaterr = ifaterr + 1
     endif
     !Glauber 2014
     if (idriz .lt. 0 .or. idriz .gt. 7) THEN
        print*,'FATAL - IDRIZ OUT OF RANGE'
        IFATERR = IFATERR + 1
     endif
     if (iccnlev .lt. 0 .or. iccnlev .gt. 2) THEN
        print*,'FATAL - ICCNLEV OUT OF RANGE: MUST BE 0-2'
        IFATERR = IFATERR + 1
     endif
     if (imbudget .lt. 0 .or. imbudget .gt. 2) THEN
        print*,'FATAL - IMBUDGET OUT OF RANGE'
        IFATERR = IFATERR + 1
     endif
     if (imbudtot .lt. 0 .or. imbudtot .gt. 2) THEN
        print*,'FATAL - IMBUDTOT OUT OF RANGE'
        IFATERR = IFATERR + 1
     endif
     if (epsil .lt. 0.05 .or. epsil .gt. 1.0) THEN
        print*,'FATAL - EPSIL OUT OF RANGE (0.05 to 1.0)'
        IFATERR = IFATERR + 1
     endif
     !Glauber 2014

    elseif (level .eq. 4) then

     ipris = 0
     isnow = 0
     iaggr = 0
     igraup = 0
     ihail = 0
     !Glauber 2014
     idriz = 0
     !Glauber 2014

    endif

    ! if level is 4, make sure that naddsc is large enough for number
    !   of bins specified in icloud.

    if (level .eq. 4) then
     if (irain .eq. 0) then
        lev4bins = 2 * icloud + 1
        if (naddsc .lt. lev4bins) then
           print*, 'fatal - naddsc is not large enough for icloud'
           print*, 'value with level = 4.'
           print*, 'naddsc must be at least ',lev4bins
           ifaterr = ifaterr + 1
        endif
     else
        lev4bins = 2 * icloud + 1 + 67
        if (naddsc .lt. lev4bins) then
           print*, 'fatal - naddsc is not large enough for icloud'
           print*, 'value with level = 4 and irain = 1.'
           print*, 'naddsc must be at least ',lev4bins
           ifaterr = ifaterr + 1
        endif
     endif
   endif

  ! if level is 5, make sure that naddsc is large enough for number
  !   of bins specified in icloud.

   if (level .eq. 5) then
     if (irain .eq. 0) then
        lev5bins = 2 * (icloud + ipris + iaggr + igraup) + 5
        if (naddsc .lt. lev5bins) then
           print*, 'fatal - naddsc is not large enough for icloud,'
           print*, 'ipris, iaggr, and igraup values with level = 5.'
           print*, 'naddsc must be at least ',lev5bins
           ifaterr = ifaterr + 1
        endif
     else
        lev5bins = 2 * (icloud + ipris + iaggr + igraup) + 5 + 67
        if (naddsc .lt. lev5bins) then
           print*, 'fatal - naddsc is not large enough for icloud,'
           print*, 'ipris, iaggr, and igraup values with level = 5'
           print*, 'and irain = 1.'
           print*, 'naddsc must be at least ',lev5bins
           ifaterr = ifaterr + 1
        endif
     endif
   endif

  ENDIF

  ! stop the run if there are any fatal errors.  list how many
  !   warning and informative errors.

  if (ifaterr > 0) then
     print*,' -----------opspec1--------------------------'
     print*,' fatal     errors - ',ifaterr
     print*,' -----------------------------------------------'
     call fatal_error(h//" fatal errors at namelist")
  end if
end subroutine opspec1

! ************************************************************************

subroutine opspec2

  ! check that fine mesh is a valid subset of its coarser mesh.
  !   and that top and bottom boundary flags are set correctly.
  !   (severity - f)

  use mem_grid, only: &
       ngrids,  &
       nxtnest, &
       ninest,  &
       njnest,  &
       nknest,  &
       nnxp,    &
       nnyp,    &
       nstratx, &
       nstraty, &
       nnsttop, &
       nnstbot, &
       nestz1,  &
       nnzp,    &
       nrz

  use mem_varinit, only: &
       vwaittot, &
       vwait1

  use mem_radiate, only: ISWRTYP, ILWRTYP ! Intent(in)
  use mem_globrad, only: raddatfn ! Intent(in)

  implicit none

  integer :: icm,ifm,ifaterr,ncx,ncy,nfx,nfxp,nfy,nfyp
  integer :: ng,nesta,nfz,kc
  character(len=*), parameter :: h="**(opspec2)**"
  logical :: ex

  ifaterr=0

  do ifm=1,ngrids

     icm = nxtnest(ifm)
     if (icm .ge. 1) then

        ncx=(nnxp(ifm)-2)/nstratx(ifm)
        ncy=(nnyp(ifm)-2)/nstraty(ifm)

        if ((nnyp(ifm).eq.1.and.nnyp(icm).ne.1).or.  &
             (nnyp(ifm).ne.1.and.nnyp(icm).eq.1)) then
           print*,' fatal - grids must be either all 3-d or all 2-d'
           ifaterr=ifaterr+1
        endif

        if (ninest(ifm).lt.3) then
           print 11, ifm
11         format (' fatal - nest #',i3,' too close to western'  &
                ,' boundary of coarser mesh')
           ifaterr=ifaterr+1
        endif

        if (njnest(ifm).lt.3.and.nnyp(ifm).gt.1) then
           print 12, ifm
12         format (' fatal - nest #',i3,' too close to southern'  &
                ,' boundary of coarser mesh')
           ifaterr=ifaterr+1
        endif

        if (nknest(ifm).lt.3.and.nnstbot(ifm).eq.0) then
           print 13, ifm
13         format (' fatal - nest #',i3,' too close to lower'  &
                ,' boundary of coarser mesh or nnstbot incorrect')
           ifaterr=ifaterr+1
        endif

        if (nknest(ifm).ne.1.and.nnstbot(ifm).eq.1) then
           print 14, ifm
14         format (' fatal - nest #',i3,' not to lower boundary of'  &
                ,' coarser mesh or nnstbot flag set incorrectly')
           ifaterr=ifaterr+1
        endif

        if (ninest(ifm)+ncx.gt.nnxp(icm)-3) then
           print 15, ifm
15         format (' fatal - nest #',i3,' too close to eastern'  &
                ,' boundary of coarser mesh')
           ifaterr=ifaterr+1
        endif

        if (njnest(ifm)+ncy.gt.nnyp(icm)-3.and.nnyp(ifm).gt.1) then
           print 16, ifm
16         format (' fatal - nest #',i3,' too close to northern'  &
                ,' boundary of coarser mesh')
           ifaterr=ifaterr+1
        endif

        if (ncx.lt.2.or.(ncy.lt.2.and.nnyp(ifm).ne.1)) then
           print 17, ifm
17         format (' fatal - nest #',i3,' dimensioned too small'  &
                ,' in at least one horizontal direction')
           ifaterr=ifaterr+1
        endif

        nfx=ncx*nstratx(ifm)+2
        if (nnxp(ifm).ne.nfx) then
           nfxp=nfx+nstratx(ifm)
           print 18, ifm,nfxp,nfx
18         format (' fatal - nest #',i3,' nnxp incompatible with'  &
                ,' nstratx:  may increase nnxp to',i5  &
                ,' or decrease it to',i5)
           ifaterr=ifaterr+1
        endif

        nfy=ncy*nstraty(ifm)+2
        if (nnyp(ifm).ne.nfy.and.nnyp(ifm).gt.1) then
           nfyp=nfy+nstraty(ifm)
           print 19, ng,nfyp,nfy
19         format (' fatal - nest #',i3,' nnyp incompatible with'  &
                ,' nstraty:  may increase nnyp to',i5  &
                ,' or decrease it to',i5)
           ifaterr=ifaterr+1
        endif

        if (nnstbot(ifm).eq.1.and.nnstbot(icm).ne.1) then
           print 20, ifm
20         format (' fatal - nest #',i3,' nnstbot flag'  &
                ,' incompatible with coarser mesh')
           ifaterr=ifaterr+1
        endif

        if (nnsttop(ifm).eq.1.and.nnsttop(icm).ne.1) then
           print 21, ifm
21         format (' fatal - nest #',i3,' nnsttop flag incompatible'  &
                ,' with coarser mesh')
           ifaterr=ifaterr+1
        endif

     endif
  enddo

  nesta=abs(nestz1)
  if(nestz1.ne.0.and.nesta.le.ngrids)then
     nfz=nnzp(nesta)-2
     kc=nknest(nesta)
1002 continue
     kc=kc+1
     nfz=nfz-nrz(kc,nesta)
     if(nfz.lt.0)then
        print 195,nesta
195     format(' fatal - vertically nested grid #',i3,  &
             ' has illegal number of levels for given nstratz values')
        ifaterr=ifaterr+1
     endif
     if(nfz.gt.0)go to 1002
     if(nfz.eq.0)then
        if(kc.gt.nnzp(nxtnest(nesta))-3.and.nnsttop(nesta).ne.1)then
           print 22, nesta
22         format (' fatal - nest #',i3,' too high'  &
                ,' or nnsttop flag set incorrectly')
           ifaterr=ifaterr+1
        endif

        if(kc.ne.nnzp(nxtnest(nesta))-1.and.nnsttop(nesta).eq.1)then
           print 23, nesta
23         format (' fatal - nest #',i3,' not to upper boundary of'  &
                ,' coarser mesh or nnsttop flag set incorrectly')
           ifaterr=ifaterr+1
        endif
     endif
  endif

  do ifm = 1,ngrids
     icm = nxtnest(ifm)
     if (ifm .ne. nesta .and. icm .ge. 1)then
        if(nnzp(ifm).gt.nnzp(icm)-nknest(ifm)-1.and.nnsttop(ifm).eq.0)then
           print 24, ifm
24         format (' fatal - nest #',i3,' nnzp incompatible with'  &
                ,' parent grid or nnsttop flag set incorrectly')
           ifaterr=ifaterr+1
        endif

        if(nnzp(ifm).ne.nnzp(icm)-nknest(ifm)+1.and.nnsttop(ifm).eq.1)then
           print 25, ifm
25         format (' fatal - nest #',i3,' not to upper boundary of'  &
                ,' coarser mesh or nnsttop flag set incorrectly')
           ifaterr=ifaterr+1
        endif
     endif
  enddo

  ! Checking problems with CARMA Radiation - BRAMS 4
  if (ISWRTYP==4 .or. ILWRTYP==4) then
     if(trim(raddatfn)=='') then
        ifaterr=ifaterr+1
        print *,'FATAL ERROR: rad data file name not set (empty).'
        print *,'Please, check the RADDATFN variable in RAMSIN.'
        print *, "Program will stop."
     else
        inquire(FILE=raddatfn(1:len_trim(raddatfn)),exist=ex)
        if (.not.ex) then
           ifaterr=ifaterr+1
           print *, "FATAL ERROR: file=", trim(raddatfn), &
                " does not exist."
           print *, "Program will stop."
        endif
     endif
  endif
  !

  ! this need to be done here since varfiles are filled before opspec3
  if (vwaittot.lt.vwait1) then
     print*,'total wait time must be <= individual varfile wait'
     print*,'      resetting vwaittot to ',vwait1
     vwaittot=vwait1
  endif

  ! stop the run if there are any fatal errors.  list how many
  !   warning and informative errors.


  if (ifaterr > 0) then
     print*,' -----------opspec2--------------------------'
     print*,' fatal     errors - ',ifaterr
     print*,' -----------------------------------------------'
     call fatal_error(h//" fatal errors at namelist")
  end if
end subroutine opspec2

! ********************************************************************

subroutine opspec3

  use mem_varinit
  use mem_grid
  use micphys
  use io_params
  use mem_radiate
  use mem_cuparm
  use mem_turb
  use mem_leaf

  use mem_grell_param, only : CLOSURE_TYPE ! INTENT(IN)

  ! TEB_SPM
  use teb_spm_start, only: TEB_SPM ! INTENT(IN)
  use mem_emiss, only: ichemi, isource ! INTENT(IN)

!--(DMK-CCATT-INI)-----------------------------------------------------
  use ccatt_start, only: &
       ccatt
  use mem_stilt, only: &
       iexev,          &
       imassflx
  use chem1_list, only: &
       nspecies,        &
       PhotojMethod   ! intent(in)
  use chem1aq_list, only: &
       nspeciesaq            ! intent(in)
  use mem_chem1, only: &
       chemistry,      &
       chem_assim,     &
       split_method,   &
       chem_timestep,  &
       nsrc,           &
       diur_cycle,     &
       bburn,          &
       src_name  ! intent(in)
  use mem_chem1aq, only: &
       chemistry_aq
  use mem_aer1, ONLY:  &
                        aerosol&
		       ,aer_timestep

  use aer1_list, ONLY:  &
                        aerosol_mechanism
  use chem_sources, only: &
       def_proc_src
  use shcu_vars_const, only: &
       nnshcu
  use modIau

!  use catt_start, only: CATT ! INTENT(IN)
!--(DMK-CCATT-END)-----------------------------------------------------

  implicit none

  integer :: ip,k,ifaterr,iwarerr,infoerr,ng,ngr
  character(len=*), parameter :: h="**(opspec3)**"
  character(len=*), parameter :: header="**(opspec3)**"
  character(len=*), parameter :: version="5.4"

!--(DMK-CCATT-INI)-----------------------------------------------------
  integer :: isrc
!--(DMK-CCATT-FIM)-----------------------------------------------------

  ifaterr=0
  iwarerr=0
  infoerr=0

  ! TEB_SPM
  !##########################################################################
  !EDF - Adition tho check if isource is activated if ichemi is
  !##########################################################################
  if (TEB_SPM==1) then
     if (ichemi==1) then
        if (isource==0) then
           print*, 'FATAL - The SPM can not be activated without sources.'
           print*, '        ISOURCE must be equal to 1.'
           IFATERR = IFATERR + 1
        endif
     endif
  endif
  !##########################################################################

  ! CCATT

!--(DMK-CCATT-INI)-----------------------------------------------------
  if (ccatt == 1) then
     ! Consistency in CCATT
     ! Checking the tracers
!--(DMK-CCATT-OLD)-----------------------------------------------------
!  if (CATT==1) then
!     if (naddsc < 4) then
!        print*, 'FATAL - If using CCATT, the variable NADDSC must be >= 4.'
!        IFATERR = IFATERR + 1
!     endif
!--(DMK-CCATT-FIM)-----------------------------------------------------
     ! Fatal error if using Shaved-ETA
     if (IF_ADAP==1) then
        print*, 'FATAL - It is not allowed running CCATT with Shaved-ETA vertical coordinate system.'
        IFATERR = IFATERR + 1
     endif
  endif

  ! check that moisture is turned on if radiation is used.
  !   (severity - f)

  if(ilwrtyp+iswrtyp.gt.0.and.level.eq.0)then
     print*,' fatal  - radiation scheme must be run with moisture.'
     ifaterr=ifaterr+1
  endif

  if(ilwrtyp==1 .or. ilwrtyp==2 .or. ilwrtyp==3 .or. ilwrtyp==5 .or. &
     iswrtyp==1 .or. iswrtyp==2 .or. iswrtyp==3 .or. iswrtyp==5 ) then
     print*,' fatal  - radiation schemes 1,2,3 and 5 are not allowed'
     print*,' fatal  - BRAMS 5.2+ allows only 4 (CARMA) or 6 (RRTM)'
     ifaterr=ifaterr+1
  endif

  if(isfcl ==3  .or. isfcl == 4) then
     print*,' fatal  - surface schemes 3 or 4 are not allowed'
     print*,' fatal  - BRAMS 5.2+ allows only 1,2 (LEAF) or 5 (JULES)'
     ifaterr=ifaterr+1
  endif

  ! microphysics flags and parameter settings

!  if((irain.ge.2.and.irain.le.4.and.rparm.le.0.)  &
!       .or.(icloud.ge.2.and.icloud.le.5.and.cparm.le.0.)  &
!       .or.(ipris.ge.2.and.ipris.le.4.and.pparm.le.0.)  &
!       .or.(isnow.ge.2.and.isnow.le.4.and.sparm.le.0.)  &
!       .or.(igraup.ge.2.and.igraup.le.4.and.gparm.le.0.)  &
!       .or.(iaggr.ge.2.and.iaggr.le.4.and.aparm.le.0.)  &
!       .or.(ihail.ge.2.and.ihail.le.4.and.hparm.le.0.)) then
!     print 26,ng,rparm,pparm,sparm,gparm,aparm,hparm
!26   format (' fatal - microphysics - xparm must be positive'  &
!          ,' if micro flags are set to 2, 3, or 4,'  &
!          ,' or up to 5 for icloud. ',i3,5f10.7)
!     ifaterr=ifaterr+1
!  endif
!Glauber 2014 ---------------------------
!  if((irain.ge.2.and.irain.le.4.and.rparm.le.0.)  &
!       .or.(icloud.ge.2.and.icloud.le.5.and.cparm.le.0.)  &
!       .or.(ipris.ge.2.and.ipris.le.4.and.pparm.le.0.)  &
!       .or.(isnow.ge.2.and.isnow.le.4.and.sparm.le.0.)  &
!       .or.(igraup.ge.2.and.igraup.le.4.and.gparm.le.0.)  &
!       .or.(iaggr.ge.2.and.iaggr.le.4.and.aparm.le.0.)  &
!       .or.(ihail.ge.2.and.ihail.le.4.and.hparm.le.0.)) then
!     print 26,ng,rparm,pparm,sparm,gparm,aparm,hparm
!26   format (' fatal - microphysics - xparm must be positive'  &
!          ,' if micro flags are set to 2, 3, or 4,'  &
!          ,' or up to 5 for icloud. ',i3,5f10.7)
!     ifaterr=ifaterr+1
!  endif
!print*,cparm,dparm,rparm,pparm,Sparm,gparm,Aparm,hparm
IF(  (icloud .GE. 2 .AND. icloud .LE. 7 .AND. cparm .LE. 0.)  &
 .OR.(idriz  .GE. 2 .AND. idriz  .LE. 7 .AND. dparm .LE. 0.)  &
 .OR.(irain  .GE. 2 .AND. irain  .LE. 4 .AND. rparm .LE. 0.)  &
 .OR.(ipris  .GE. 2 .AND. ipris  .LE. 7 .AND. pparm .LE. 0.)  &
 .OR.(isnow  .GE. 2 .AND. isnow  .LE. 4 .AND. Sparm .LE. 0.)  &
 .OR.(igraup .GE. 2 .AND. igraup .LE. 4 .AND. gparm .LE. 0.)  &
 .OR.(iaggr  .GE. 2 .AND. iaggr  .LE. 4 .AND. Aparm .LE. 0.)  &
 .OR.(ihail  .GE. 2 .AND. ihail  .LE. 4 .AND. hparm .LE. 0.)) THEN
   PRINT 26,ng,cparm,dparm,rparm,pparm,sparm,gparm,aparm,hparm
   26 FORMAT (' FATAL - Microphysics - xPARM must be positive'  &
             ,' if micro flags are set to 2, 3, or 4,'  &
             ,' or up to 5 for icloud. ',i3,5f10.7)
   IFATERR=IFATERR+1
ENDIF
IF(idriz.GE.5 .and. icloud.lt.5) then
   PRINT*,' FATAL - Microphysics - ICLOUD must be >= 5 if IDRIZ >= 5'
   IFATERR=IFATERR+1
ENDIF
IF(iccnlev.GT.0 .and. icloud.lt.5) then
   PRINT*,' FATAL - Microphysics - ICCNLEV must be 0 if ICLOUD < 5'
   IFATERR=IFATERR+1
ENDIF
!IF(level.lt.3 .and. IDUST.EQ.1) then
!   PRINT*,' FATAL - DustModel - Dust source model requires Level=3 Micro'
!   IFATERR=IFATERR+1
!ENDIF
!IF(level.lt.3 .and. ISALT.EQ.1) then
!   PRINT*,' FATAL - Sea-Salt Model - Salt source model requires Level=3 Micro'
!   IFATERR=IFATERR+1
!ENDIF
!IF(level.lt.3 .and. (IAMSFLG.EQ.1 .or. ISS1FLG.EQ.1 .or. ISS2FLG.EQ.1 &
!    .or. IMD1FLG.EQ.1 .or. IMD2FLG.EQ.1)) then
!   PRINT*,' FATAL - Aerosol - Radiative aerosol requires Level=3 Micro'
!   IFATERR=IFATERR+1
!ENDIF

!IF(ibubble.lt.0 .and. ibubble.gt.2) then
!   PRINT*,' FATAL - IBUBBLE must be 0, 1, or 2'
!   PRINT*,'         0 = off, 1 = RAMSIN set bubble, 2 = Random bubble in ruser'
!   IFATERR=IFATERR+1
!ENDIF
!Glauber 2014 ---------------------------

  ! convective parameterization flags and parameter settings

  do ng=1,ngrids
     if (nnqparm(ng)/=0 .and. level==0) then
        print 27
27      format (' fatal - level must be at least'  &
             ,' 1 for the cumulus parameterization')
        ifaterr = ifaterr + 1
     endif
  enddo

! time integration schemes
  IF (dyncore_flag .gt. 3) THEN
     PRINT *, 'FATAL - Only the hybrid time integration schemes are allowed in 5.2'
     call fatal_error('set dyncore_flag == 0 to 3')
     IFATERR=IFATERR+1
  ENDIF

  IF (pd_or_mnt_constraint .gt. 1 .or. pd_or_mnt_constraint .lt. 0) THEN
     call fatal_error('set pd_or_mnt_constraint == 0 to 1')
     IFATERR=IFATERR+1
  ENDIF
!print *, 'Ordens:', order_h, order_v
  IF (order_h .gt. 6 .or. order_h .lt. 1) THEN
     call fatal_error('set order_h == 1 to 6')
     IFATERR=IFATERR+1
  ENDIF

  IF (order_v .gt. 6 .or. order_v .lt. 1) THEN
     call fatal_error('set order_v == 1 to 6')
     IFATERR=IFATERR+1
  ENDIF

!--(DMK-CCATT-INI)-----------------------------------------------------
![ML
! Complete Exner tendency and vertical coordinate.
  IF (IEXEV .EQ. 2 .AND. IF_ADAP .NE. 0) THEN
     call fatal_error('IEXEV cannot be set to 2 with shaved-eta coordinate')
     IFATERR=IFATERR+1
  ENDIF
!srf
  IF (IEXEV .EQ. 2 .AND. level .EQ. 0) THEN
     call fatal_error('IEXEV cannot be set to 2 with microphyics level = 0 ')
     IFATERR=IFATERR+1
  ENDIF

  !Mass flux  cannot be output with Kuo parameterization
  DO NG=1, NGRIDS
     IF (IMASSFLX == 1 .AND. NNQPARM(NG) == 1) THEN
        PRINT *, 'FATAL - Convective mass flux cannot be used with Kuo convective parameterization (NNQPARM=1)'
        IFATERR=IFATERR+1
     ENDIF

     IF (IMASSFLX == 1 .AND. IDIFFK(NG) /= 7) THEN
        PRINT *, 'FATAL - Mass flux output can be used only with Nakanishi parameterization (IDIFFK=7)'
        IFATERR=IFATERR+1
     ENDIF
  ENDDO
!ML]
!--(DMK-CCATT-END)-----------------------------------------------------

  ! moving grids and topography interpolation

  do ng=1,ngrids
     if((abs(gridu(ng)) .gt. 1.e-20 .or. abs(gridv(ng)) .gt. 1.e-20)  &
          .and. itoptflg(ng) .ne. 0) then
        print 28
28      format (' fatal - nested grid topography must be interpolated'  &
             ,' from parent grid if nested grid moves')
        ifaterr=ifaterr+1
     endif
  enddo

  ! check horizontal and vertical grid spacings.

  if(dzmax.lt.deltaz)then
     print*,' warning - deltaz is being reduced by a low value',  &
          ' of dzmax.'
     iwarerr=iwarerr+1
  endif

  if(dzrat.gt.1.2)then
     print*,' warning - large vertical stretch ratios sacrifice',  &
          ' second order accuracy in the vertical differencing.'
     iwarerr=iwarerr+1
  endif

  ! check numerical schemes.

  if((sspct.lt.0.2.or.sspct.gt.1.0).and.sspct.ne.0.0)then
     print*,' warning - sspct should normally range from 0.2 to 1.0'
     iwarerr=iwarerr+1
  endif

  if(nfpt.gt.nnzp(1))then
     print*,' fatal - nfpt must be less than nnzp(1).'
     ifaterr=ifaterr+1
  endif

  if(iadvl.ne.2.and.iadvl.ne.4)then
     print*,' fatal - iadvl must be 2 or 4'
     ifaterr=ifaterr+1
  endif

  if(iadvf.ne.2.and.iadvf.ne.6)then
     print*,' fatal - iadvf must be 2 or 6'
     ifaterr=ifaterr+1
  endif

  ! check turbulence parameterization.

  do ngr=1,ngrids

     !_stc
     !_stc.............................................
     !_stc correction to account for the new options
     !_stc for the e-l and e-eps closures: idiffk <= 6
     !_stc  (s. trini castelli)
     !_stc.............................................
     !_stc  if(idiffk(ngr).lt.1.or.idiffk(ngr).gt.4)then
     !_stc    print*,' fatal - idiffk must be 1, 2, 3,or 4.'
     !_stc
     if(idiffk(ngr).lt.1.or.idiffk(ngr).gt.8)then
        print*,' fatal - idiffk must be 1, 2, 3, 4, 5, 6 , 7 or 8'
        ifaterr=ifaterr+1
     endif
     !srf-opt
     if(level.lt.1.AND.idiffk(ngr).EQ.7)then
        print*,' fatal - idiffk 7 cannot be used with microphysics level 0'
        ifaterr=ifaterr+1
     endif
  enddo
  ! check that diffusion flags are compatible if using ihorgrad=1

  if(ihorgrad.eq.2)then
     if(idiffk(ngr) >= 3)then
        print*,' fatal - cant use ihorgrad=2 if idiffk >= 3'
        ifaterr=ifaterr+1
     endif
  endif

  ! check whether the soil model will be run and make sure that the
  !   number of soil levels are correct.(severity - f,i )

  if(isfcl.eq.0.and.nzg.gt.1)then
     print*,' info  - more soil levels specified than needed.'
     infoerr=infoerr+1
  endif

  if (isfcl == 0 .and. npatch /= 2) then
     print*, ' fatal  - when isfcl = 0, npatch must be 2. '
     ifaterr = ifaterr + 1
  endif

  if(isfcl.gt.0.and.nzg.le.2)then
     print*,  &
          ' fatal  - at least 2 soil levels are needed for soil'  &
          ,' model.'
     ifaterr=ifaterr+1
  endif

  do k=1,nzg
     if (slz(k) .gt. -.001) then
        print*, 'fatal - soil level',k,' not (enough) below ground'  &
             ,' level'
        ifaterr=ifaterr+1
     endif
  enddo

  do k=1,nzg-1
     if (slz(k)-slz(k+1) .gt. .001) then
        print*, 'fatal - soil level',k,' not (enough) deeper than'  &
             ,' soil level',k+1
        ifaterr=ifaterr+1
     endif
  enddo

  ! if the soil model will be run with no radiation, make a suggestion
  !   that the radiation be turned on. (severity - f )

  do ngr=1,ngrids
     if(isfcl.gt.0.and.ilwrtyp+iswrtyp.eq.0)then
        print*,' fatal  - radiation scheme must be run with soil',  &
             ' model.'
        ifaterr=ifaterr+1
     endif
  enddo


  ! make sure that if nudging, nudging time scales are greater than
  ! the model coarse grid timestep, and that rayleigh friction nudging
  ! is not done with variable initialization.

  if (initial .eq. 1) then
     if (nfpt .gt. 0 .and. distim .gt. 0. .and.  &
          distim .lt. dtlongn(1)) then
        print*, 'rayleigh friction nudging is being done'
        print*, 'and distim is less than dtlongn(1).'
        print*, 'this nudging is too strong.'
        ifaterr=ifaterr+1
     endif
  endif

  if (initial .eq. 2) then

     if (nfpt .gt. 0 .and. distim .gt. 0.) then
        print*, 'rayleigh friction nudging may not be used when'
        print*, 'variable initialization is used.'
        ifaterr=ifaterr+1
     endif

     if (nudlat .ge. 1 .and. tnudlat .gt. 0. .and.  &
          tnudlat .lt. dtlongn(1)) then
        print*, 'lateral boundary nudging is being done'
        print*, 'and tnudlat is less than dtlongn(1).'
        print*, 'this nudging is too strong.'
        ifaterr=ifaterr+1
     endif

     if (tnudcent .gt. 0. .and. tnudcent .lt. dtlongn(1)) then
        print*, 'center nudging is being done'
        print*, 'and tnudcent is less than dtlongn(1).'
        print*, 'this nudging is too strong.'
        ifaterr=ifaterr+1
     endif
     if (tnudtop .gt. 0. .and. tnudtop .lt. dtlongn(1)) then
        print*, 'top boundary nudging is being done'
        print*, 'and tnudtop is less than dtlongn(1).'
        print*, 'this nudging is too strong.'
        ifaterr=ifaterr+1
     endif
   
     if (applyIAU > 0 .and. tnudcent .gt. 0.  ) then
         print*, 'The IAU procedure requires no nudging in the center domain,'
	 print*, 'set tnudcent equal to zero.'
         ifaterr=ifaterr+1
     endif

  endif


  !     check the averaging and analysis frequencies for consistency.

  if (abs(avgtim).gt.0.0.and.frqmean.le.0.0.and.frqboth.le.0.) then
     print*,'have frqmean=0 & frqboth=0 even though avgtim=',avgtim
     print*,'respecifying avgtim=0.'
     avgtim=0.
     iwarerr=iwarerr+1
  endif
  if (frqlite.gt.0.) then
     !   if ( nl3d(1)+nl2d(1)+nl3ds(1).eq.0)then
     !      print*,'have no lite variables even though frqlite=',frqlite
     !      print*,'respecify in vtables or set frqlite=0 in namelist'
     !      ifaterr=ifaterr+1
     !   endif
  endif
  if (frqmean.gt.0.0.and.abs(avgtim).gt.0.) then
     if ( abs(avgtim).gt.frqmean ) then
        print*,'avgtim must be <= frqmean'
        ifaterr=ifaterr+1
     endif
     !   if ( nm3d(1)+nm2d(1)+nm3ds(1).eq.0)then
     !      print*,'have no mean variables even though frqmean=',frqmean
     !      print*,'respecify in vtables or set frqmean=0 in namelist'
     !      ifaterr=ifaterr+1
     !   endif
  endif
  if (frqboth.gt.0.0.and.abs(avgtim).gt.0.) then
     if ( abs(avgtim).gt.frqboth ) then
        print*,'avgtim must be <= frqboth'
        ifaterr=ifaterr+1
     endif
     !   if ( nb3d(1)+nb2d(1)+nb3ds(1).eq.0)then
     !      print*,'have no both variables even though frqboth=',frqboth
     !      print*,'respecify in vtables or set frqboth=0 in namelist'
     !      ifaterr=ifaterr+1
     !   endif
  endif
  if (frqmean.gt.0.0.and.frqboth.gt.0.0.and.abs(avgtim).gt.0.) then
     if ( (frqmean.gt.frqboth.and.mod(frqmean,frqboth).ne.0.).or.  &
          (frqmean.lt.frqboth.and.mod(frqboth,frqmean).ne.0.) ) then
        print*,'frqmean must be a multiple of frqboth or vice versa'
        ifaterr=ifaterr+1
     endif
  endif

  ! check printout parameters.

  do ngr=1,ngrids
     do ip=1,nplt
        if(ixsctn(ip).eq.1.and.isbval(ip).gt.nnyp(ngr))then
           print 1,ip,ngr
1          format (' fatal - isbval(',i2,') is out of bounds in'  &
                ,' y-direction for grid number ',i2,'.')
           ifaterr=ifaterr+1
        elseif(ixsctn(ip).eq.2.and.isbval(ip).gt.nnxp(ngr))then
           print 2,ip,ngr
2          format (' fatal - isbval(',i2,') is out of bounds in'  &
                ,' x-direction for grid number ',i2,'.')
           ifaterr=ifaterr+1
        elseif(ixsctn(ip).eq.3.and.isbval(ip).gt.nnzp(ngr))then
           print 3,ip,ngr
3          format (' fatal - isbval(',i2,') is out of bounds in'  &
                ,' z-direction for grid number ',i2,'.')
           ifaterr=ifaterr+1
        endif
     enddo
  enddo

!if( (mcphys_type == 2 .or. mcphys_type == 3) .and. (isfcl<=2 .and. isfcl >=1)) then
!    print*,' FATAL - cannot use LEAF-3 scheme with GT microphysics'
!    print*," FATAL - arrays qpcpg and dpcpg are not avalaible in this scheme"
!    IFATERR=IFATERR+1
!endif

do ngr=1,ngrids
!   if( (mcphys_type == 2 .or. mcphys_type == 3) .and. nnshcu(NGR) ==2) then
!    print*,' FATAL - cannot use nnhscu = 2 scheme with GT microphysics'
!    IFATERR=IFATERR+1
!   endif
   if( (mcphys_type >= 2 ) .and. nnqparm(NGR) ==2) then
    print*,' FATAL - cannot use nnqparm = 2 scheme with  microphysics >= 2'
    IFATERR=IFATERR+1
   endif
enddo


if(ccatt == 0 ) then
  IF( CHEMISTRY >= 0) then
    print*,' FATAL - Cant set CCATT = 0 and  CHEMISTRY >= 0 '
    IFATERR=IFATERR+1
  endif
endif

if(ccatt == 1 ) then
  IF( CHEMISTRY >= 0 .and. nspecies == 0) then
    print*,' FATAL - Cant use CHEMISTRY with nspecies = 0 '
    IFATERR=IFATERR+1
  endif
  IF( CHEMISTRY < 0 .and. CHEM_ASSIM == 1) then
    print*,' FATAL - Cant use  4DDA CHEM ASSIM  with CHEMISTRY < 0'
    !print*,'  setting  zero to parameter CHEM_ASSIM'
    !CHEM_ASSIM = 0
    IFATERR=IFATERR+1
  endif
  IF( CHEMISTRY < 0 .and. AEROSOL > 0) then
    print*,' FATAL - For AEROSOL levels 1 or 2, Chemistry must be >= 0'
    IFATERR=IFATERR+1
  endif
  IF(CHEMISTRY >= 0) then
   IF( trim(aerosol_mechanism) /= 'MATRIX' .and. AEROSOL == 2 ) then
     print*,' FATAL - AEROSOL = 2 reguires Aerosol_mechanism MATRIX'
     IFATERR=IFATERR+1
   endif
 ENDIF
 IF( AEROSOL <= 1 .and. aer_timestep .ne. dtlong) then
    print*,'AER_TIMESTEP for AEROSOL <= 1 must be equal to dtlong'
    print*,"AER_TIMESTEP will be defined as DTLONG"
    AER_TIMESTEP = DTLONG
  endif

  IF( CHEMISTRY < 1 .and. CHEMISTRY_AQ > 0) then
    print*,' FATAL - Cant use  CHEMISTRY_AQ  with CHEMISTRY < 1'
    IFATERR=IFATERR+1
  endif

  IF( CHEMISTRY_AQ > 0 .and. nspeciesaq == 0 ) then
    print*,' FATAL - Cant use  CHEMISTRY_AQ  with nspeciesaq = 0'
    print*,' Check your chemical mechanism used to create spack files'
    IFATERR=IFATERR+1
  endif

!  IF( CHEMISTRY < 0 .and. trim(PhotojMethod) == 'FAST-JX') then
!     print*,' FATAL - Cant use  FAST=JX without chemistry '
!     IFATERR=IFATERR+1
!  endif
  IF(CHEMISTRY > 0) then
   IF( (ilwrtyp .ne. 4 .or. iswrtyp .ne. 4) .AND. trim(PhotojMethod) == 'FAST-JX')THEN
    PRINT*,' FATAL  -  CARMA radiation scheme must be run with FAST-JX WHEN CHEMISTRY  IS ON.'
    IFATERR=IFATERR+1
   ENDIF
   IF (ilwrtyp == 6 .or. iswrtyp == 6) THEN
      IF(trim(PhotojMethod) /= 'FAST-TUV' .and. trim(PhotojMethod) /= 'LUT' ) THEN
         PRINT*,' FATAL  -  CARMA/RRTM radiation schemes must be run with FAST-TUV WHEN CHEMISTRY  IS ON.'
         IFATERR=IFATERR+1
      END IF
   ENDIF
 ENDIF

  IF( CHEMISTRY >= 0) then
     IF( trim(adjustl(def_proc_src)) == 'last_sources') def_proc_src='LAST_SOURCES'
     IF( trim(adjustl(def_proc_src)) == 'stop'        ) def_proc_src='STOP'

     IF( trim(adjustl(def_proc_src)) .ne. 'STOP'        .and. &
         trim(adjustl(def_proc_src)) .ne. 'LAST_SOURCES') THEN
         PRINT*,'FATAL - unknow def_proc_src parameter: ',def_proc_src
         PRINT*,'Allowed def_proc_src parameters are : LAST_SOURCES or STOP'
         IFATERR=IFATERR+1
     endif
     do isrc=1,nsrc
       if(diur_cycle(isrc) .ne. 0 .and. diur_cycle(isrc) .ne. 1) then
          print*, 'Wrong value for diur_cycle=',src_name(isrc),' :',diur_cycle(isrc)
          print*, 'Allowed vaues are 0 or 1'
	  IFATERR=IFATERR+1
       endif
     enddo
     if(diur_cycle(bburn) .ne. 1) then
          print*, 'Diur_cycle for bburn must be =1'
	  IFATERR=IFATERR+1
     endif

  endif


  IF(CHEMISTRY > 0) then
     IF( trim(adjustl(SPLIT_METHOD)) == 'sequential') SPLIT_METHOD='SEQUENTIAL'
     IF( trim(adjustl(SPLIT_METHOD)) == 'symmetric' ) SPLIT_METHOD='SYMMETRIC'
     IF( trim(adjustl(SPLIT_METHOD)) == 'parallel'  ) SPLIT_METHOD='PARALLEL'

     IF(trim(adjustl(SPLIT_METHOD)) .ne. 'SEQUENTIAL' .and. &
        trim(adjustl(SPLIT_METHOD)) .ne. 'SYMMETRIC'  .and. &
        trim(adjustl(SPLIT_METHOD)) .ne. 'PARALLEL'	    ) THEN
	PRINT*,'FATAL - unknow splitting method for chemistry int: ',SPLIT_METHOD
	PRINT*,'Allowed methods are : PARALLEL , SYMMETRIC or SEQUENTIAL'
        IFATERR=IFATERR+1
     ENDIF
     IF(CHEMISTRY == 1) then
         PRINT*,'Chemistry with solver 1 - QSSA - is now obsolet'
	 PRINT*,'Only solver 2 - 3 or 4 is recommended'
	 IFATERR=IFATERR+1
        if(trim(adjustl(SPLIT_METHOD)) .ne. 'PARALLEL') then
	  PRINT*,'FATAL - chemistry 1 must be used with SPLIT_METHOD = PARALLEL'
          IFATERR=IFATERR+1
        ENDIF
        if(CHEM_TIMESTEP .ne. DTLONG) then
	  PRINT*,'FATAL - chemistry 1 => CHEM_TIMESTEP = DTLONG'
	  CHEM_TIMESTEP = DTLONG
        ENDIF
     ENDIF
  ELSEIF(CHEMISTRY == 0) then
	!PRINT*,'WARNING - chemistry 0 requires SPLIT_METHOD = PARALLEL and CHEM_TIMESTEP = DTLONG'
	!PRINT*,'Warning - chemistry 0 must be used CHEM_TIMESTEP = DTLONG'
	!PRINT*,'Warning - model will run with these settings'
        SPLIT_METHOD  = 'PARALLEL'
	CHEM_TIMESTEP = DTLONG
  ENDIF

endif
  do ng=1,ngrids
  ![MLO - Blocking Grell deep/shallow convection without TKE
     if (nnshcu(ng) == 2 .or. NNQPARM(NG) == 2  .or. NNQPARM(NG) == 5 .or. NNQPARM(NG) == 6) then
        if (idiffk(ng) == 2 .or. idiffk(ng) == 3) then
           print *, 'FATAL - deep (nnqparm 2 or  5 or 6) and shallow (nnshcu 2) requires turbulence scheme with TKE (1,4,5,6,7)'
           print *, 'Please change your setup for grid ',ng,'...'
           IFATERR=IFATERR+1
        endif
     endif
  end do
  !MLO]

!--(DMK-CCATT-END)-----------------------------------------------------
!  if (ISWRTYP==4 .or. ILWRTYP==4) then
!       if( (nnqparm(ng) /= 5 .and.&
!	    nnqparm(ng) /= 6 .and.&
!	    nnqparm(ng) /= 2)	 )    THEN
!          ifaterr=ifaterr+1
!          print *,'FATAL ERROR: CARMA radiation requires NNQPARM= 2, 5 or 6 *'
!          print *, "Program will stop."
!    endif
!  endif
  do ng=1,ngrids
   if (ISWRTYP==6 .or. ILWRTYP==6) then
     if( mcphys_type <= 1 .and. icloud < 5 )    THEN
       ifaterr=ifaterr+1
       print *,'FATAL ERROR: RRTM radiation requires ICLOUD >=5'
       !print*,"values=",icloud
       print *, "Program will stop."
     endif
   endif
  end do


  ! stop the run if there are any fatal errors.  list how many
  !   warning and informative errors.

  ! Checking closure type for New Grell Param.
  if  (CLOSURE_TYPE/='EN' .and. &
       CLOSURE_TYPE/='GR' .and. &
       CLOSURE_TYPE/='LO' .and. &
       CLOSURE_TYPE/='MC' .and. &
       CLOSURE_TYPE/='SC' .and. &
       CLOSURE_TYPE/='AS' .and. &
       CLOSURE_TYPE/='PB') then
     print *, "FATAL - Grell Closure type ERROR!"
     print *, "        Program will stop!"
     ifaterr=ifaterr+1
  endif


  if (ifaterr > 0) then
     print*,' -----------opspec3--------------------------'
     print*,' fatal     errors - ',ifaterr
     print*,' warning   errors - ',iwarerr
     print*,' inform  messages - ',infoerr
     print*,' -----------------------------------------------'
     call fatal_error(h//" fatal errors at namelist")
  end if
end subroutine opspec3
