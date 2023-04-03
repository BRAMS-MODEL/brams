!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

!
! gridloc_prt: prints position on Earth of each grid,
!              both in lat/lon as in km
!

subroutine gridloc_prt()

  use mem_grid, only: &
       oneGlobalGridData, ngrids, nnxp, nnyp, nnzp, &
       xtn, xmn, ytn, ymn, zmn, platn, plonn, deltaxn, deltayn

  use mem_varinit, only: &
    nudlat

  implicit none

  real :: centx
  real :: centy
  real :: latCenter
  real :: lonCenter
  integer :: nx
  integer :: ny
  integer :: ngr
  integer :: midx
  integer :: midy
  real :: NudLatMax
  real :: NudLonMax

  include "constants.f90"

  character(len=*), parameter :: h="**(gridloc_prt)**"

  do ngr=1,ngrids

     ! local nx, ny just for convinience

     nx=nnxp(ngr); ny=nnyp(ngr)

     ! NW and NE LAT/LON

     write(logUnit,"(/,a,i4,a,/)") &
          "---------------------------LOCATION AND DIMENSIONS OF GRID", ngr, &
          "---------------------------------------"

      write(*,*) 'Model geographics: NW & NE: '
     write(logUnit,101) &
          oneGlobalGridData(ngr)%global_glat(1,ny), &
          oneGlobalGridData(ngr)%global_glon(1,ny), &
          oneGlobalGridData(ngr)%global_glat(nx,ny), &
          oneGlobalGridData(ngr)%global_glon(nx,ny)
101  format('  NW LAT/LON      NE LAT/LON (DEG)= '  &
          ,2X,'(',F7.3,',',F8.3,')',4X,'(',F7.3,',',F8.3,')')

     NudLatMax=(oneGlobalGridData(ngr)%global_glat(1,2)-oneGlobalGridData(ngr)%global_glat(1,1))*nudlat
     NudLonMax=(oneGlobalGridData(ngr)%global_glon(2,1)-oneGlobalGridData(ngr)%global_glon(1,1))*nudlat
     write(*,*)
      write(*,*) 'Minimum Boundaries IC geographics: NW & NE: '
     write(logUnit,115) &
          oneGlobalGridData(ngr)%global_glat(1,ny)+NudLatMax, &
          oneGlobalGridData(ngr)%global_glon(1,ny)-NudLonMax, &
          oneGlobalGridData(ngr)%global_glat(nx,ny)+NudLatMax, &
          oneGlobalGridData(ngr)%global_glon(nx,ny)+NudLonMax

115  format('  NW LAT/LON      NE LAT/LON (DEG)= '  &
          ,2X,'(',F7.3,',',F8.3,')',4X,'(',F7.3,',',F8.3,')')

     midx = (nx+1)/2
     midy = (ny+1)/2
     if (nx/2 == midx) then
        centx = xmn(nx/2,ngr)
     else
        centx = xtn(midx,ngr)
     endif
     if (ny/2 == midy) then
        centy = ymn(ny/2,ngr)
     else
        centy = ytn(midy,ngr)
     endif
     call xy_ll(latCenter,lonCenter,platn(ngr),plonn(ngr)  &
          ,centx,centy)
     write(*,*)
     write(logUnit,102)latCenter,lonCenter
     write(*,*)
102  format('        CENT LAT/LON (DEG)        = '  &
          ,13X,'(',F7.3,',',F8.3,')')

    write(*,*) 'Model geographics: SW & SE: '
     write(logUnit,105) &
          oneGlobalGridData(ngr)%global_glat(1,1), &
          oneGlobalGridData(ngr)%global_glon(1,1), &
          oneGlobalGridData(ngr)%global_glat(nx,1), &
          oneGlobalGridData(ngr)%global_glon(nx,1)
105  format('  SW LAT/LON      SE LAT/LON (DEG)= '  &
       ,2X,'(',F7.3,',',F8.3,')',4X,'(',F7.3,',',F8.3,')',/)

    write(*,*) 'Minimum Boundaries IC geographics: SW & SE: '
     write(logUnit,118) &
          oneGlobalGridData(ngr)%global_glat(1,1)-NudLatMax, &
          oneGlobalGridData(ngr)%global_glon(1,1)-NudLonMax, &
          oneGlobalGridData(ngr)%global_glat(nx,1)-NudLatMax, &
          oneGlobalGridData(ngr)%global_glon(nx,1)+NudLonMax
118  format('  SW LAT/LON      SE LAT/LON (DEG)= '  &
       ,2X,'(',F7.3,',',F8.3,')',4X,'(',F7.3,',',F8.3,')',/)


     write(logUnit,106)xtn(1,ngr),xtn(nx,ngr)
106  format('     WEST COORD (KM)=',-3PF10.3  &
          ,'               EAST COORD (KM)=',-3PF10.3)

     write(logUnit,107)ytn(1,ngr),ytn(ny,ngr)
107  format('    SOUTH COORD (KM)=',-3PF10.3  &
          ,'              NORTH COORD (KM)=',-3PF10.3)

     write(logUnit,108)zmn(1,ngr),zmn(nnzp(ngr)-1,ngr)
108  format('    BOTTOM COORD (M)=',F10.1  &
          ,'            TOP COORDINATE (M)=',F10.1)

     write(logUnit,109)deltaxn(ngr),deltayn(ngr)
109  format('         DELTA-X (M)=',F10.1  &
          ,'                   DELTA-Y (M)=',F10.1)

     write(logUnit,110)zmn(2,ngr)-zmn(1,ngr)
110  format('  BOTTOM DELTA-Z (M)=',F10.1)

     write(logUnit,111) &
          maxval(oneGlobalGridData(ngr)%global_glat(:,ny)), &
          minval(oneGlobalGridData(ngr)%global_glat(:,1)), &
          minval(oneGlobalGridData(ngr)%global_glon(1,:)), &
          maxval(oneGlobalGridData(ngr)%global_glon(nx,:))
111  format(/,' MAX N LAT=',f9.3,'  MIN S LAT=',f9.3  &
          ,/,' MIN W LON=',f9.3,' MAX E LON=',f9.3 )

  enddo

  write(logUnit,"(104('-'),//)")

  return
end subroutine gridloc_prt
!
!     ******************************************************************
!
subroutine refs3d(n1,n2,n3,pi0,dn0,dn0u,dn0v,th0,topt,rtgt)

  use mem_grid, only: &
       nxp, nyp, if_adap, ngrid, zt, nzp, ztop, dzm

  use mem_scratch, only: &
       vctr2

  use ref_sounding, only: &
       pi01dn, th01dn

  use rconstants, only: &
       g, cpor, cp, p00, rgas

  implicit none
  integer, intent(in)  :: n1
  integer, intent(in)  :: n2
  integer, intent(in)  :: n3
  real,    intent(out) :: pi0(n1,n2,n3)
  real,    intent(out) :: dn0(n1,n2,n3)
  real,    intent(out) :: dn0u(n1,n2,n3)
  real,    intent(out) :: dn0v(n1,n2,n3)
  real,    intent(out) :: th0(n1,n2,n3)
  real,    intent(in)  :: topt(n2,n3)
  real,    intent(in)  :: rtgt(n2,n3)

  integer :: i,j,k
  real :: c1,c2,c3

  ! +---------------------------------------------------------------------
  ! _    This routine initializes the 3-D reference state arrays from the
  !        1-D reference state.
  ! +---------------------------------------------------------------------

  do j=1,n3
     do i=1,n2

        if (if_adap == 1) then

           do k = 1,n1
              pi0(k,i,j) = pi01dn(k,ngrid)
              th0(k,i,j) = th01dn(k,ngrid)
           end do
           c1 = g * 2.

        else

           do k = 1,n1
              vctr2(k) = zt(k) * rtgt(i,j) + topt(i,j)
           end do
           call htint(nzp,pi01dn(1,ngrid),zt,nzp,pi0(1,i,j),vctr2)
           call htint(nzp,th01dn(1,ngrid),zt,nzp,th0(1,i,j),vctr2)
           c1 = g * 2. * (1. - topt(i,j) / ztop)

        end if

        c2 = 1. - cpor
        c3 = cp ** c2
        do k = n1-1,1,-1
           pi0(k,i,j) = pi0(k+1,i,j)  &
                + c1 / ((th0(k,i,j) + th0(k+1,i,j)) * dzm(k))
        end do

        do k = 1,n1
           dn0(k,i,j) = (c3 * p00) / (rgas * th0(k,i,j) * pi0(k,i,j) ** c2)
        end do

     end do
  end do

  call fill_dn0uv(n1,n2,n3,dn0,dn0u,dn0v)
end subroutine refs3d

!     *****************************************************************

subroutine fill_dn0uv(n1,n2,n3,dn0,dn0u,dn0v)
  implicit none
  integer, intent(in)  :: n1
  integer, intent(in)  :: n2
  integer, intent(in)  :: n3
  real,    intent(in)  :: dn0(n1,n2,n3)
  real,    intent(out) :: dn0u(n1,n2,n3)
  real,    intent(out) :: dn0v(n1,n2,n3)
  integer :: i,j,k,i1,j1

  do j = 1,n3
     j1 = min(j+1,n3)
     do i = 1,n2
        i1 = min(i+1,n2)
        do k = 1,n1
           dn0u(k,i,j) = .5 * (dn0(k,i,j) + dn0(k,i1,j))
           dn0v(k,i,j) = .5 * (dn0(k,i,j) + dn0(k,i,j1))
        end do
     end do
  end do
end subroutine fill_dn0uv

!****************************************************************************

subroutine toptsmth(n2,n3,topt,vt2da,dxt,dyt,dxm,dym)

  use mem_grid, only: &
       jdim, nxp, nyp, nx, ny1

  use io_params, only: &
       ntopsmth, izflat

  implicit none
  integer :: n2,n3
  real :: topt(n2,n3),vt2da(n2,n3),dxt(n2,n3),dyt(n2,n3)  &
       ,dxm(n2,n3),dym(n2,n3)

  integer :: i,j,iter

  !          Smooth the topography if desired.

  do iter=1,ntopsmth
     do j=jdim+1,ny1*jdim+1
        do i=2,nx
           vt2da(i,j)=.25/(dxt(i,j)*dyt(i,j))*(  &
                ( (topt(i+1,j)-topt(i,j))*dxm(i,j)  &
                -(topt(i,j)-topt(i-1,j)  )*dxm(i-1,j))*dxt(i,j)  &
                +( (topt(i,j+jdim)-topt(i,j))*dym(i,j)  &
                -(topt(i,j)-topt(i,j-jdim))*dym(i,j-jdim))*dyt(i,j))
        enddo
     enddo
     do j=1,nyp
        do i=1,nxp
           topt(i,j)=topt(i,j)+vt2da(i,j)
        enddo
     enddo
  enddo

  !     FLATTEN TOPOGRAPHY AT BOUNDARIES

  if(izflat.gt.0)then
     do j=1,nyp
        do i=1,izflat
           topt(i,j)=topt(izflat+1,j)
           topt(nxp+1-i,j)=topt(nxp-izflat,j)
        enddo
     enddo

     if(jdim==1)then
        do j=1,izflat
           do i=1,nxp
              topt(i,j)=topt(i,izflat+1)
              topt(i,nyp+1-j)=topt(i,nyp-izflat)
           enddo
        enddo
     endif

  endif
  return
end subroutine toptsmth


!=================================================================


!     *****************************************************************







subroutine FieldInit(initflg)

  use mem_grid, only: &
       ngrid, initial

  use node_mod, only: &
       mxp, myp, mzp, ibcon

  use mem_basic, only: &
       basic_g

  use mem_micro, only: &
       micro_g

  use micphys, only: &
       icloud, irain, ipris, isnow, iaggr, igraup, ihail

  implicit none

  integer, intent(in) :: initflg

  ! finish initializing past time level variables

  if (initial /= 3) then

     basic_g(ngrid)%up = basic_g(ngrid)%uc
     basic_g(ngrid)%vp = basic_g(ngrid)%vc
     basic_g(ngrid)%wp = basic_g(ngrid)%wc
     basic_g(ngrid)%pp = basic_g(ngrid)%pc

     if(initflg==1) then
        if(icloud>=1) micro_g(ngrid)%rcp = 0.
        if(irain>=1)  micro_g(ngrid)%rrp = 0.
        if(ipris>=1)  micro_g(ngrid)%rpp = 0.
        if(isnow>=1)  micro_g(ngrid)%rsp = 0.
        if(iaggr>=1)  micro_g(ngrid)%rap = 0.
        if(igraup>=1) micro_g(ngrid)%rgp = 0.
        if(ihail>=1)  micro_g(ngrid)%rhp = 0.
        basic_g(ngrid)%wp = 0.
     endif

     call tkeinit(mzp,mxp,myp)

  endif

  call dumset(mzp,mxp,myp,1,mxp,1,myp,ibcon,basic_g(ngrid)%uc(1,1,1),'U')
  call dumset(mzp,mxp,myp,1,mxp,1,myp,ibcon,basic_g(ngrid)%vc(1,1,1),'V')
  call dumset(mzp,mxp,myp,1,mxp,1,myp,ibcon,basic_g(ngrid)%wc(1,1,1),'W')
  call dumset(mzp,mxp,myp,1,mxp,1,myp,ibcon,basic_g(ngrid)%up(1,1,1),'U')
  call dumset(mzp,mxp,myp,1,mxp,1,myp,ibcon,basic_g(ngrid)%vp(1,1,1),'V')
  call dumset(mzp,mxp,myp,1,mxp,1,myp,ibcon,basic_g(ngrid)%wp(1,1,1),'W')
end subroutine FieldInit
