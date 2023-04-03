!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine advectc(varn,mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
!> @brief: advectc
!! @date:  18/Nov/2015
!! @version:  5.2
!! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
!! 
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
  use grid_dims, only: maxgrds
  use mem_tend, only: tend
  use var_tables, only: num_scalar, scalar_tab
  use mem_scratch, only: scratch, vctr1, vctr2
  use mem_grid, only: ngrid, nzpmax, grid_g, dtlt, if_adap, jdim, time, &
       zt, zm, dzm, dzt, hw4, dyncore_flag

  use mem_basic, only: basic_g

  !--(DMK-CCATT-INI)-----------------------------------------------------
  !-srf for aerosols sedimentation
  use ccatt_start, only: &
       ccatt      ! intent(in)

  use mem_aer1, only: &
       aerosol,    &
       num_scalar_aer_1st

  use mem_chem1, only: &
       nspecies_transported

  use module_dry_dep, only: &
       dd_sedim,         &
       fa_preptc_with_sedim

  use monotonic_adv, ONLY: advmnt 
  !-srf
  !--(DMK-CCATT-FIM)-----------------------------------------------------

  use ModNamelistFile, only: &
         namelistfile

  implicit none
  include "i8.h"
  integer :: mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum,n
  integer(kind=i8) :: mxyzp
  character(len=*) :: varn

  real :: dtlto2
  real, dimension(nzpmax*7) :: advscr
  real, dimension(maxgrds), save :: save_dtlt
  integer :: i,j,k,ind

  real, pointer :: scalarp, scalart

  !--(DMK-CCATT-INI)-----------------------------------------------------
  integer :: i_scl
  !--(DMK-CCATT-FIM)-----------------------------------------------------

  integer, dimension(maxgrds), save :: ncall
  data ncall/maxgrds*0/

  !--(DMK-LFR NEC-SX6)----------------------------------------------
  !real, external :: valugp

  !--(DMK-LFR NEC-SX6)----------------------------------------------

  if (ncall(ngrid) == 0 .or. dtlt .ne. save_dtlt(ngrid) ) then
     ncall(ngrid) = 1
     save_dtlt(ngrid) = dtlt
  endif
  mxyzp = mxp * myp * mzp
  !  basic_g(ngrid)%up    (1:nzp,1:nxp,1:nyp) = 10.
  !  basic_g(ngrid)%uc    (1:nzp,1:nxp,1:nyp) = 10.
  !  basic_g(ngrid)%vp    (1:nzp,1:nxp,1:nyp) = 0.
  !  basic_g(ngrid)%vc    (1:nzp,1:nxp,1:nyp) = 0.
  !  basic_g(ngrid)%wp    (1:nzp,1:nxp,1:nyp) = 0.
  !  basic_g(ngrid)%wc    (1:nzp,1:nxp,1:nyp) = 0.

  if (trim(varn) .eq. 'V' .or. trim(varn) .eq. 'ALL') then
     ! Advect  U, V, and W

     if (if_adap == 0) then

        if ( (dyncore_flag == 0) .or. (dyncore_flag == 1) ) then 

          call vel_advectc(mzp,mxp,myp,ia,iz,ja,jz,izu,jzv  &
             ,tend%ut              (1)      &
             ,tend%vt              (1)     ,tend%wt              (1)      &
             ,scratch%vt3da        (1)      &
             ,scratch%vt3db        (1)     ,scratch%vt3dc        (1))

        else if (dyncore_flag == 2) then 
          !MB: this branch must probably be transferred also to the other if-statements in this subroutine!
     	  print*, "velocity advection"
          call vel_advectc(mzp,mxp,myp,ia,iz,ja,jz,izu,jzv  &
             ,tend%ut_rk           (1)      &
             ,tend%vt_rk           (1)     ,tend%wt_rk           (1)      &
             ,scratch%vt3da        (1)      &
             ,scratch%vt3db        (1)     ,scratch%vt3dc        (1))
        end if
	
     else

        call vel_advectc_adap(mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,jdim    &
             ,grid_g(ngrid)%lpu   (1,1)   ,grid_g(ngrid)%lpv   (1,1)    &
             ,grid_g(ngrid)%lpw   (1,1)   ,basic_g(ngrid)%uc   (1,1,1)  &
             ,basic_g(ngrid)%vc   (1,1,1) ,basic_g(ngrid)%wc   (1,1,1)  &
             ,tend%ut             (1)     ,tend%vt             (1)      &
             ,tend%wt             (1)     ,basic_g(ngrid)%dn0  (1,1,1)  &
             ,basic_g(ngrid)%dn0u (1,1,1) ,basic_g(ngrid)%dn0v (1,1,1)  &
             ,grid_g(ngrid)%aru   (1,1,1) ,grid_g(ngrid)%arv   (1,1,1)  &
             ,grid_g(ngrid)%arw   (1,1,1) ,grid_g(ngrid)%volu  (1,1,1)  &
             ,grid_g(ngrid)%volv  (1,1,1) ,grid_g(ngrid)%volw  (1,1,1)  &
             ,scratch%vt3da       (1)     ,scratch%vt3db       (1)      &
             ,scratch%vt3dc       (1)     ,scratch%vt3dd       (1)      &
             ,scratch%vt3de       (1)     ,scratch%vt3df       (1) ,time)

     endif

  endif

  
  if (trim(varn) .eq. 'T' .or. trim(varn) .eq. 'SCALAR') THEN

     ! Advect  scalars

     dtlto2 = .5 * dtlt
     ind = 0
     do j = 1,myp
        do i = 1,mxp
           do k = 1,mzp
              ind = ind + 1
              scratch%vt3da(ind) = (basic_g(ngrid)%up(k,i,j)  &
                   + basic_g(ngrid)%uc(k,i,j)) * dtlto2
              scratch%vt3db(ind) = (basic_g(ngrid)%vp(k,i,j)  &
                   + basic_g(ngrid)%vc(k,i,j)) * dtlto2
              scratch%vt3dc(ind) = (basic_g(ngrid)%wp(k,i,j)  &
                   + basic_g(ngrid)%wc(k,i,j)) * dtlto2
              ! if(i==4.and.j==6.and.k==2) print*,'adv:',dtlt  &
              ! ,basic_g(ngrid)%up(k,i,j),basic_g(ngrid)%vp(k,i,j),basic_g(ngrid)%wp(k,i,j)&
              ! ,basic_g(ngrid)%uc(k,i,j),basic_g(ngrid)%vc(k,i,j),basic_g(ngrid)%wc(k,i,j)
           enddo
        enddo
     enddo

     if (if_adap == 0) then

        call fa_preptc(mzp,mxp,myp        &
             ,scratch%vt3da        (1)     ,scratch%vt3db        (1)      &
             ,scratch%vt3dc        (1)     ,scratch%vt3dd        (1)      &
             ,scratch%vt3de        (1)     ,scratch%vt3df        (1)      &
             ,scratch%vt3dh        (1)     ,scratch%vt3di        (1)      &
             ,scratch%vt3dj        (1)     ,scratch%vt3dk        (1)      &
             ,mynum                         )

     else

        call fa_preptc_adap(mzp,mxp,myp                               &
             ,scratch%vt3da       (1)     ,scratch%vt3db       (1)      &
             ,scratch%vt3dc       (1)     ,scratch%vt3dd       (1)      &
             ,scratch%vt3de       (1)     ,scratch%vt3df       (1)      &
             ,scratch%vt3dh       (1)     ,basic_g(ngrid)%dn0  (1,1,1)  &
             ,basic_g(ngrid)%dn0u (1,1,1) ,basic_g(ngrid)%dn0v (1,1,1)  &
             ,grid_g(ngrid)%aru   (1,1,1) ,grid_g(ngrid)%arv   (1,1,1)  &
             ,grid_g(ngrid)%arw   (1,1,1) ,grid_g(ngrid)%volt  (1,1,1)  &
             ,grid_g(ngrid)%dxu   (1,1)   ,grid_g(ngrid)%dyv   (1,1)    &
             ,grid_g(ngrid)%dxt   (1,1)   ,grid_g(ngrid)%dyt   (1,1)    &
             ,zt,zm,dzm,vctr1,vctr2,jdim,mynum                          )

     endif

     !--(DMK-CCATT-INI)-----------------------------------------------------
     !- combine old adv for thermo + micro with mnt adv for tracers
     if( advmnt== 0) then
        i_scl=num_scalar(ngrid)  !- all scalars 
     elseif(advmnt == 2) then
        i_scl= num_scalar(ngrid) - NSPECIES_TRANSPORTED !- only theta_il+water+tke
     elseif(advmnt == 3) then
        i_scl=1  !- only theta_il
     endif

     do n=1,i_scl
        
        !- if RK or ABM3 schemes, THP/THC are not transported here
        if (dyncore_flag == 2 .or. dyncore_flag == 3) then
          if (scalar_tab(n,ngrid)%name == 'THC' .or. &
              scalar_tab(n,ngrid)%name == 'THP') cycle
        endif
        
	!tmp
        !srf - somente para TEMP, RV, RTP, TKE ...
        !      do n=1,num_scalar(ngrid) - NSPECIES_TRANSPORTED
        !      if (scalar_tab(n,ngrid)%name .ne. 'COP') cycle
        !      print*,' scal orig=',n,scalar_tab(n,ngrid)%name; call flush(6)     
        !      do n=1,1! for only thp
        !      if (scalar_tab(n,ngrid)%name /= 'THP') stop 999
        !srf -  somente para TEMP, RV, RTP, TKE ...
        !tmp

        if(ccatt == 1 .and. aerosol == 1) then
           if( n >= num_scalar_aer_1st ) then 
              !cycle! para testar ADVMNT
              !srf-  We are going to include sedimentation of aerosols at
              !      vertical advection tendency. It is supposed that scalars 
              !      with  N >= num_scalar_aer_1st are _all_ aerosols with 
              !      sedimentation speeds previously calculated.
              if (if_adap == 0) then
	         !print*,'applying sedim for aerosol=', n; call flush(6)
                 call fa_preptc_with_sedim(mzp,mxp,myp	    &
                      ,scratch%vt3da             ,scratch%vt3db	      &
                      ,scratch%vt3dc              ,scratch%vt3df &
                      ,scratch%vt3dk	      &
                      ,basic_g(ngrid)%dn0        ,grid_g(ngrid)%rtgt       &
                      ,grid_g(ngrid)%f13t        ,grid_g(ngrid)%f23t       &
                                !srf- aerosol section 
                      ,dtlt  		    &
                      ,N			    & ! current scalar being transported
                      ,num_scalar_aer_1st	    & ! 1st aerosol at scalar table
                      ,basic_g(ngrid)%wp         & ! air vertical velocity (P time)
                      ,basic_g(ngrid)%wc         & ! air vertical velocity (C time)
                      ,scratch%vt3dp             & ! ) ! to save horizontal contribution on the sigmaz velocity
                      ,nzpmax,hw4,dzm,dzt,dd_sedim(:,ngrid))   ! (DMK) deposicao seca (cod. limpo)
              else
                 print*,'sedim not yet prepared for shaved eta'
                 stop 3333
              endif
           endif
        endif
        !--(DMK-CCATT-FIM)-----------------------------------------------------

        scalarp => scalar_tab(n,ngrid)%var_p
        scalart => scalar_tab(n,ngrid)%var_t
        call atob_long(mxyzp, scalarp, scratch%scr1(1))

        if (if_adap == 0) then

           call fa_xc(mzp,mxp,myp,ia,iz,1,myp        &
                ,scalarp           ,scratch%scr1  (1)  &
                ,scratch%vt3da (1) ,scratch%vt3dd (1)  &
                ,scratch%vt3dg (1) ,scratch%vt3dh (1)  &
                ,scratch%vt3di (1) ,mynum              )

           if (jdim .eq. 1)  &
                call fa_yc(mzp,mxp,myp,ia,iz,ja,jz        &
                ,scalarp           ,scratch%scr1  (1)  &
                ,scratch%vt3db (1) ,scratch%vt3de (1)  &
                ,scratch%vt3dg (1) ,scratch%vt3dj (1)  &
                ,scratch%vt3di (1) ,jdim,mynum         )

           call fa_zc(mzp,mxp,myp,ia,iz,ja,jz        &
                ,scalarp           ,scratch%scr1  (1)  &
                ,scratch%vt3dc (1) ,scratch%vt3df (1)  &
                ,scratch%vt3dg (1) ,scratch%vt3dk (1)  &
                ,vctr1,vctr2,mynum                     )

           call advtndc(mzp,mxp,myp,ia,iz,ja,jz    &
                ,scalarp          ,scratch%scr1 (1)  &
                ,scalart          ,dtlt,mynum        )

        else

           call fa_xc_adap(mzp,mxp,myp,ia,iz,1,myp         &
                ,grid_g(ngrid)%lpw (1,1) ,scalarp            &
                ,scratch%scr1      (1)   ,scratch%vt3da (1)  &
                ,scratch%vt3dd     (1)   ,scratch%vt3dg (1)  &
                ,scratch%vt3dh     (1)   ,mynum              )

           if (jdim .eq. 1)                                &
                call fa_yc_adap(mzp,mxp,myp,ia,iz,ja,jz         &
                ,grid_g(ngrid)%lpw (1,1) ,scalarp            &
                ,scratch%scr1      (1)   ,scratch%vt3db (1)  &
                ,scratch%vt3de     (1)   ,scratch%vt3dg (1)  &
                ,scratch%vt3dh (1)  &
                ,jdim                    ,mynum              )

           call fa_zc_adap(mzp,mxp,myp,ia,iz,ja,jz         &
                ,grid_g(ngrid)%lpw (1,1) ,scalarp            &
                ,scratch%scr1      (1)   ,scratch%vt3dc (1)  &
                ,scratch%vt3df     (1)   ,scratch%vt3dg (1)  &
                ,scratch%vt3dh     (1)   ,vctr1              &
                ,vctr2                   ,mynum              )

           call advtndc_adap(mzp,mxp,myp,ia,iz,ja,jz  &
                ,grid_g(ngrid)%lpw (1,1) ,scalarp       &
                ,scratch%scr1      (1)   ,scalart       &
                ,dtlt                    ,mynum         )

        endif

     enddo

  endif

  return
end subroutine advectc

!     *********************************************************************



subroutine vel_advectc(m1,m2,m3,ia,iz,ja,jz,izu,jzv  &
     ,ut,vt,wt,flxu,flxv,flxw)
!> @brief: vel_advectc
!! @date:  18/Nov/2015
!! @version:  5.2
!! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
!! 
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
  use mem_grid
  use mem_basic
  use node_mod, only :  mynum

  implicit none

  integer,intent(in) :: m1,m2,m3,ia,iz,ja,jz,izu,jzv
  real, dimension(m1,m2,m3) :: ut,vt,wt,flxu,flxv,flxw


  integer :: j,i,k,jm,im
  real :: c1z,c1x,c1y

  ! Compute momentum fluxes flxu, flxv, flxw

  do j = 1,m3
     do i = 1,m2
        do k = 1,m1
           flxu(k,i,j) = basic_g(ngrid)%uc(k,i,j) * basic_g(ngrid)%dn0u(k,i,j) * grid_g(ngrid)%rtgu(i,j)  &
                * grid_g(ngrid)%fmapui(i,j)
           flxv(k,i,j) = basic_g(ngrid)%vc(k,i,j) * basic_g(ngrid)%dn0v(k,i,j) * grid_g(ngrid)%rtgv(i,j)  &
                * grid_g(ngrid)%fmapvi(i,j)
        enddo
     enddo
  enddo

  if(itopo.eq.0) then
     do j = 1,m3
        do i = 1,m2
           do k = 1,m1-1
              flxw(k,i,j) = basic_g(ngrid)%wc(k,i,j)  &
                   * .5 * (basic_g(ngrid)%dn0(k,i,j) + basic_g(ngrid)%dn0(k+1,i,j))
           enddo
        enddo
     enddo
  else
     do j = 1,m3
        jm = max(j-1,1)
        do i = 1,m2
           im = max(i-1,1)
           do k = 1,m1-1
              flxw(k,i,j) = basic_g(ngrid)%wc(k,i,j)  &
                   * .5 * (basic_g(ngrid)%dn0(k,i,j) + basic_g(ngrid)%dn0(k+1,i,j))  &
                   + hw4(k) * ((flxu(k,i,j) + flxu(k+1,i,j)  &
                   + flxu(k,im,j) + flxu(k+1,im,j)) * grid_g(ngrid)%f13t(i,j)  &
                   + (flxv(k,i,j) + flxv(k+1,i,j)  &
                   + flxv(k,i,jm) + flxv(k+1,i,jm)) * grid_g(ngrid)%f23t(i,j))
           enddo
        enddo
     enddo
  endif

  ! Compute advection contribution to U tendency

  do j = ja,jz
     do i = ia,izu
        c1z = .25 / grid_g(ngrid)%rtgu(i,j)
        c1x = c1z * grid_g(ngrid)%fmapu(i,j) * grid_g(ngrid)%dxu(i,j)
        c1y = c1z * grid_g(ngrid)%fmapu(i,j) * grid_g(ngrid)%dyu(i,j)

        do k = 2,m1-1
           ut(k,i,j) = ut(k,i,j) + c1x / basic_g(ngrid)%dn0u(k,i,j) * (  &
                (flxu(k,i,j) + flxu(k,i-1,j))  &
                * (basic_g(ngrid)%uc(k,i,j) + basic_g(ngrid)%uc(k,i-1,j))  &
                - (flxu(k,i,j) + flxu(k,i+1,j))  &
                * (basic_g(ngrid)%uc(k,i,j) + basic_g(ngrid)%uc(k,i+1,j))  &
                + (flxu(k,i+1,j) - flxu(k,i-1,j)) * 2.* basic_g(ngrid)%uc(k,i,j) )
        enddo

        do k = 2,m1-1
           ut(k,i,j) = ut(k,i,j) + c1y / basic_g(ngrid)%dn0u(k,i,j) * (  &
                (flxv(k,i,j-jdim) + flxv(k,i+1,j-jdim))  &
                * (basic_g(ngrid)%uc(k,i,j) + basic_g(ngrid)%uc(k,i,j-jdim))  &
                - (flxv(k,i,j) + flxv(k,i+1,j))  &
                * (basic_g(ngrid)%uc(k,i,j) + basic_g(ngrid)%uc(k,i,j+jdim))&
                + (flxv(k,i,j) + flxv(k,i+1,j) - flxv(k,i,j-jdim)  &
                - flxv(k,i+1,j-jdim)) * 2.* basic_g(ngrid)%uc(k,i,j) )
        enddo

        do k = 2,m1-1
           ut(k,i,j) = ut(k,i,j) + c1z * dzt(k) / basic_g(ngrid)%dn0u(k,i,j) * (  &
                (flxw(k-1,i,j) + flxw(k-1,i+1,j))  &
                * (basic_g(ngrid)%uc(k,i,j) + basic_g(ngrid)%uc(k-1,i,j))  &
                - (flxw(k,i,j) + flxw(k,i+1,j))  &
                * (basic_g(ngrid)%uc(k,i,j) + basic_g(ngrid)%uc(k+1,i,j))   &
                + (flxw(k,i,j) + flxw(k,i+1,j) - flxw(k-1,i,j)  &
                - flxw(k-1,i+1,j)) * 2.* basic_g(ngrid)%uc(k,i,j) )
        enddo
     enddo
  enddo

  ! Compute advection contribution to V tendency

  do j = ja,jzv
     do i = ia,iz
        c1z = .25 / grid_g(ngrid)%rtgv(i,j)
        c1x = c1z * grid_g(ngrid)%fmapv(i,j) * grid_g(ngrid)%dxv(i,j)
        c1y = c1z * grid_g(ngrid)%fmapv(i,j) * grid_g(ngrid)%dyv(i,j)

        do k = 2,m1-1
           vt(k,i,j) = vt(k,i,j) + c1x / basic_g(ngrid)%dn0v(k,i,j) * (  &
                (flxu(k,i-1,j) + flxu(k,i-1,j+jdim))  &
                * (basic_g(ngrid)%vc(k,i,j) + basic_g(ngrid)%vc(k,i-1,j))  &
                - (flxu(k,i,j) + flxu(k,i,j+jdim))  &
                * (basic_g(ngrid)%vc(k,i,j) + basic_g(ngrid)%vc(k,i+1,j))  &
                + (flxu(k,i,j) + flxu(k,i,j+jdim) - flxu(k,i-1,j)  &
                - flxu(k,i-1,j+jdim)) * 2.* basic_g(ngrid)%vc(k,i,j) )
        enddo

        do k = 2,m1-1
           vt(k,i,j) = vt(k,i,j) + c1y / basic_g(ngrid)%dn0v(k,i,j) * (  &
                (flxv(k,i,j) + flxv(k,i,j-jdim))  &
                * (basic_g(ngrid)%vc(k,i,j) + basic_g(ngrid)%vc(k,i,j-jdim))  &
                - (flxv(k,i,j) + flxv(k,i,j+jdim))  &
                * (basic_g(ngrid)%vc(k,i,j) + basic_g(ngrid)%vc(k,i,j+jdim))  &
                + (flxv(k,i,j+jdim) - flxv(k,i,j-jdim))  &
                * 2.* basic_g(ngrid)%vc(k,i,j) )
        enddo

        do k = 2,m1-1
           vt(k,i,j) = vt(k,i,j) + c1z * dzt(k) / basic_g(ngrid)%dn0v(k,i,j) * (  &
                (flxw(k-1,i,j) + flxw(k-1,i,j+jdim))  &
                * (basic_g(ngrid)%vc(k,i,j) + basic_g(ngrid)%vc(k-1,i,j))  &
                - (flxw(k,i,j) + flxw(k,i,j+jdim))  &
                * (basic_g(ngrid)%vc(k,i,j) + basic_g(ngrid)%vc(k+1,i,j))  &
                + (flxw(k,i,j) + flxw(k,i,j+jdim) - flxw(k-1,i,j)  &
                - flxw(k-1,i,j+jdim)) * 2.* basic_g(ngrid)%vc(k,i,j) )
        enddo
     enddo
  enddo

  ! Compute advection contribution to W tendency
  !write (*,fmt='("rbasic_g(ngrid)%vc1-wt: ",I2.2,1X,I1.1,1X,4(E9.3,1X),"/",4(E9.3,1X))') mynum,1,(wt(k,2,2),k=1,4),(wt(k,3,3),k=1,4)
  !call flush(6)
  do j = ja,jz
     do i = ia,iz
        c1z = .5 / grid_g(ngrid)%rtgt(i,j)
        c1x = c1z * grid_g(ngrid)%fmapt(i,j) * grid_g(ngrid)%dxt(i,j)
        c1y = c1z * grid_g(ngrid)%fmapt(i,j) * grid_g(ngrid)%dyt(i,j)

        do k = 2,m1-2
           wt(k,i,j) = wt(k,i,j)  &
                + c1x / (basic_g(ngrid)%dn0(k,i,j) + basic_g(ngrid)%dn0(k+1,i,j)) * (  &
                (flxu(k,i-1,j) + flxu(k+1,i-1,j))  &
                * (basic_g(ngrid)%wc(k,i,j) + basic_g(ngrid)%wc(k,i-1,j))  &
                - (flxu(k,i,j) + flxu(k+1,i,j))  &
                * (basic_g(ngrid)%wc(k,i,j) + basic_g(ngrid)%wc(k,i+1,j))  &
                + (flxu(k,i,j) + flxu(k+1,i,j) - flxu(k,i-1,j)  &
                - flxu(k+1,i-1,j)) * 2.* basic_g(ngrid)%wc(k,i,j) )
        enddo

        do k = 2,m1-2
           wt(k,i,j) = wt(k,i,j)  &
                + c1y / (basic_g(ngrid)%dn0(k,i,j) + basic_g(ngrid)%dn0(k+1,i,j)) * (  &
                (flxv(k,i,j-jdim) + flxv(k+1,i,j-jdim))  &
                * (basic_g(ngrid)%wc(k,i,j) + basic_g(ngrid)%wc(k,i,j-jdim))  &
                - (flxv(k,i,j) + flxv(k+1,i,j))  &
                * (basic_g(ngrid)%wc(k,i,j) + basic_g(ngrid)%wc(k,i,j+jdim))  &
                + (flxv(k,i,j) + flxv(k+1,i,j) - flxv(k,i,j-jdim)  &
                - flxv(k+1,i,j-jdim)) * 2.* basic_g(ngrid)%wc(k,i,j) )
        enddo

        do k = 2,m1-2
           wt(k,i,j) = wt(k,i,j)  &
                + c1z * dzm(k) / (basic_g(ngrid)%dn0(k,i,j) + basic_g(ngrid)%dn0(k+1,i,j)) * (  &
                (flxw(k,i,j) + flxw(k-1,i,j))  &
                * (basic_g(ngrid)%wc(k,i,j) + basic_g(ngrid)%wc(k-1,i,j))  &
                - (flxw(k,i,j) + flxw(k+1,i,j))  &
                * (basic_g(ngrid)%wc(k,i,j) + basic_g(ngrid)%wc(k+1,i,j))   &
                + (flxw(k+1,i,j) - flxw(k-1,i,j)) * 2.* basic_g(ngrid)%wc(k,i,j) )
        enddo
     enddo
  enddo
  !write (*,fmt='("rvc2-wt: ",I2.2,1X,I1.1,1X,4(E9.3,1X),"/",4(E9.3,1X))') mynum,1,(wt(k,2,2),k=1,4),(wt(k,3,3),k=1,4)
  !call flush(6)

  return
end subroutine vel_advectc




!     *********************************************************************

subroutine fa_preptc(m1,m2,m3,vt3da,vt3db,vt3dc,vt3dd,vt3de,vt3df  &
     ,vt3dh,vt3di,vt3dj,vt3dk,mynum)
!> @brief:  fa_preptc                                      
!! @author:  unknow
!! @date:  18/Nov/2015
!! @version:  5.2
!! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
!! 
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt

  use mem_grid
  use mem_basic
  use mem_scratch
                       
  implicit none

  integer,intent(in) :: m1,m2,m3,mynum
  integer :: j,i,k,im,ip,jm,jp

  real :: c1,c2,c3,c4,rtgti

  real, dimension(m1,m2,m3) :: vt3da,vt3db,vt3dc,vt3dd,vt3de,vt3df  &
       ,vt3dh,vt3di,vt3dj,vt3dk

  ! VT3DA, VT3DB, and VT3DC are input as the velocity components (averaged
  ! between past and current time levels) times dtlt.

  ! Add contribution to VT3DC from horiz winds crossing sloping sigma surfaces,
  !    and include 1/rtgt factor in VT3DC
  ! Compute half Courant numbers: VT3DD, VT3DE, and VT3DF
  ! Compute weight at scalar point: VT3DH
  ! Compute advective weights for the linear term: VT3DI, VCTR1, and VCTR2

  do j = 1,m3
     jm = max(1,j-1)
     jp = min(m3,j+1)
     do i = 1,m2
        im = max(1,i-1)
        ip = min(m2,i+1)
        rtgti = 1. / grid_g(ngrid)%rtgt(i,j)
        c1 = .5 * grid_g(ngrid)%dxu(i,j)
        c2 = .5 * grid_g(ngrid)%dyv(i,j)
        c3 = grid_g(ngrid)%dxt(i,j) * grid_g(ngrid)%fmapt(i,j) * rtgti
        c4 = grid_g(ngrid)%dyt(i,j) * grid_g(ngrid)%fmapt(i,j) * rtgti

        do k = 1,m1-1
           vt3dc(k,i,j) = ((vt3da(k,i,j) + vt3da(k+1,i,j)  &
                + vt3da(k,im,j) + vt3da(k+1,im,j)) * grid_g(ngrid)%f13t(i,j)  &
                + (vt3db(k,i,j) + vt3db(k+1,i,j) + vt3db(k,i,jm)  &
                + vt3db(k+1,i,jm)) * grid_g(ngrid)%f23t(i,j)) * hw4(k)  &
                + vt3dc(k,i,j) * rtgti
           vt3dd(k,i,j) = c1 * vt3da(k,i,j)
           vt3de(k,i,j) = c2 * vt3db(k,i,j)
           vt3df(k,i,j) = .5 * vt3dc(k,i,j) * dzm(k)
           vctr3(k) = 1. / basic_g(ngrid)%dn0(k,i,j)
           vt3dh(k,i,j) = c3 * vctr3(k)
           vt3dj(k,i,j) = c4 * vctr3(k)
           vt3dk(k,i,j) = dzt(k) * vctr3(k)
        enddo

        !            vt3di(1,i,j) = grid_g(ngrid)%dxu(i,j) / (grid_g(ngrid)%dxu(i,j) + grid_g(ngrid)%dxt(ip,j))
        !            vt3di(2,i,j) = grid_g(ngrid)%dxu(i,j) / (grid_g(ngrid)%dxu(i,j) + grid_g(ngrid)%dxt(i,j))
        !            vt3di(3,i,j) = grid_g(ngrid)%dyv(i,j) / (grid_g(ngrid)%dyv(i,j) + grid_g(ngrid)%dyt(i,jp))
        !            vt3di(4,i,j) = grid_g(ngrid)%dyv(i,j) / (grid_g(ngrid)%dyv(i,j) + grid_g(ngrid)%dyt(i,j))
        ! temporary override
        vt3di(1,i,j) = .5
        vt3di(2,i,j) = .5
        vt3di(3,i,j) = .5
        vt3di(4,i,j) = .5
     enddo
  enddo

  do k = 1,m1-1
     vctr1(k) = (zt(k+1) - zm(k)) * dzm(k)
     vctr2(k) =  (zm(k) - zt(k)) * dzm(k)
  enddo

  ! Convert velocity components * dtlt (VT3DA, VT3DB, VT3DC)
  ! into mass fluxes times dtlt.

  do j = 1,m3
     do i = 1,m2
        c1 = grid_g(ngrid)%fmapui(i,j) * grid_g(ngrid)%rtgu(i,j)
        c2 = grid_g(ngrid)%fmapvi(i,j) * grid_g(ngrid)%rtgv(i,j)
        do k = 1,m1-1
           vt3da(k,i,j) = vt3da(k,i,j) * c1 * basic_g(ngrid)%dn0u(k,i,j)
           vt3db(k,i,j) = vt3db(k,i,j) * c2 * basic_g(ngrid)%dn0v(k,i,j)
           vt3dc(k,i,j) = vt3dc(k,i,j) * .5  &
                * (basic_g(ngrid)%dn0(k,i,j) + basic_g(ngrid)%dn0(k+1,i,j))
        enddo
     enddo
  enddo

  return
end subroutine fa_preptc

!     *********************************************************************

subroutine fa_xc(m1,m2,m3,ia,iz,ja,jz  &
     ,scp,scr1,vt3da,vt3dd,vt3dg,vt3dh,vt3di,mynum)
!> @brief: Compute scalar flux times dtlt
!! @date:  18/Nov/2015
!! @version:  5.2
!! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
!! 
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
  implicit none
  integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k,mynum

  real :: dfact
  real, dimension(m1,m2,m3) :: scp,scr1,vt3da,vt3dd,vt3dg,vt3dh,vt3di

  dfact = .5
  do j = 1,m3
     do i = 1,iz

        ! Compute scalar flux times dtlt [VT3DG]

        do k = 2,m1-1
           vt3dg(k,i,j) = vt3da(k,i,j)  &
                * (vt3di(1,i,j) * scr1(k,i,j)  &
                +  vt3di(2,i,j) * scr1(k,i+1,j)  &
                +  vt3dd(k,i,j) * (scr1(k,i,j) - scr1(k,i+1,j)))
        enddo

        ! Modify fluxes to retain positive-definiteness on scalar quantities.
        !    If a flux will remove 1/2 quantity during a timestep,
        !    reduce to first order flux. This will remain positive-definite
        !    under the assumption that ABS(CFL(i)) + ABS(CFL(i-1)) < 1.0 if
        !    both fluxes are evacuating the box.

        do k = 2,m1-1
           if (vt3da(k,i,j) .gt. 0.) then
              if (vt3dg(k,i,j) * vt3dh(k,i,j) .gt.  &
                   dfact * scr1(k,i,j)) then
                 vt3dg(k,i,j) = vt3da(k,i,j) * scr1(k,i,j)
              endif
           elseif (vt3da(k,i,j) .lt. 0.) then
              if (-vt3dg(k,i,j) * vt3dh(k,i+1,j) .gt.  &
                   dfact * scr1(k,i+1,j)) then
                 vt3dg(k,i,j) = vt3da(k,i,j) * scr1(k,i+1,j)
              endif
           endif
        enddo
     enddo
  enddo

  ! Compute flux divergence

  do j = 1,m3
     do i = ia,iz
        do k = 2,m1-1
           !if(k==2.and.j==6) print*,'flux:',i,scr1(k,i,j),vt3dg(k,i,j)  &
           !            , vt3dh(k,i,j) * (vt3dg(k,i-1,j) - vt3dg(k,i,j))  &
           !            , vt3dh(k,i,j) * scp(k,i,j) * (vt3da(k,i,j) - vt3da(k,i-1,j))
           scr1(k,i,j) = scr1(k,i,j)  &
                + vt3dh(k,i,j) * (vt3dg(k,i-1,j) - vt3dg(k,i,j)  &
                + scp(k,i,j) * (vt3da(k,i,j) - vt3da(k,i-1,j)))
        enddo
     enddo
  enddo

  return
end subroutine fa_xc

!     *********************************************************************

subroutine fa_yc(m1,m2,m3,ia,iz,ja,jz  &
     ,scp,scr1,vt3db,vt3de,vt3dg,vt3dj,vt3di,jdim,mynum)
!> @brief: Compute scalar flux 
!! @date:  18/Nov/2015
!! @version:  5.2
!! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
!! 
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
  implicit none
  integer :: m1,m2,m3,ia,iz,ja,jz,jdim,mynum,i,j,k

  real :: dfact
  real, dimension(m1,m2,m3) :: scp,scr1,vt3db,vt3de,vt3dg,vt3dj,vt3di

  dfact = .5
  do j = 1,jz
     do i = ia,iz

        ! Compute scalar flux VT3DG

        do k = 2,m1-1
           vt3dg(k,i,j) = vt3db(k,i,j)  &
                * (vt3di(3,i,j) * scr1(k,i,j)  &
                +  vt3di(4,i,j) * scr1(k,i,j+jdim)  &
                +  vt3de(k,i,j) * (scr1(k,i,j) - scr1(k,i,j+jdim)))
        enddo

        !      Modify fluxes to retain positive-definiteness on scalar quantities.
        !         If a flux will remove 1/2 quantity during a timestep,
        !         reduce to first order flux. This will remain positive-definite
        !         under the assumption that ABS(CFL(i)) + ABS(CFL(i-1)) < 1.0 if
        !         both fluxes are evacuating the box.

        do k = 2,m1-1
           if (vt3db(k,i,j) .gt. 0.) then
              if (vt3dg(k,i,j) * vt3dj(k,i,j) .gt.  &
                   dfact * scr1(k,i,j)) then
                 vt3dg(k,i,j) = vt3db(k,i,j) * scr1(k,i,j)
              endif
           elseif (vt3db(k,i,j) .lt. 0.) then
              if (-vt3dg(k,i,j) * vt3dj(k,i,j+jdim) .gt.  &
                   dfact * scr1(k,i,j+jdim)) then
                 vt3dg(k,i,j) = vt3db(k,i,j) * scr1(k,i,j+jdim)
              endif
           endif
        enddo
     enddo
  enddo

  ! Compute flux divergence

  do j = ja,jz
     do i = ia,iz
        do k = 2,m1-1
           scr1(k,i,j) = scr1(k,i,j)  &
                + vt3dj(k,i,j) * (vt3dg(k,i,j-jdim) - vt3dg(k,i,j)  &
                + scp(k,i,j) * (vt3db(k,i,j) - vt3db(k,i,j-jdim)))
        enddo
     enddo
  enddo

  return
end subroutine fa_yc

!     *********************************************************************

subroutine fa_zc(m1,m2,m3,ia,iz,ja,jz  &
     ,scp,scr1,vt3dc,vt3df,vt3dg,vt3dk,vctr1,vctr2,mynum)
!> @brief: Compute scalar flux 
!! @date:  18/Nov/2015
!! @version:  5.2
!! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
!! 
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
  implicit none
  integer :: m1,m2,m3,ia,iz,ja,jz,mynum,i,j,k

  real :: dfact
  real, dimension(m1,m2,m3) :: scp,scr1,vt3dc,vt3df,vt3dg,vt3dk
  real, dimension(*) :: vctr1,vctr2

  dfact = .5
  do j = ja,jz
     do i = ia,iz

        ! Compute scalar flux VT3DG

        do k = 1,m1-1
           vt3dg(k,i,j) = vt3dc(k,i,j)  &
                * (vctr1(k) * scr1(k,i,j)  &
                +  vctr2(k) * scr1(k+1,i,j)  &
                +  vt3df(k,i,j) * (scr1(k,i,j) - scr1(k+1,i,j)))
        enddo

        ! Modify fluxes to retain positive-definiteness on scalar quantities.
        !    If a flux will remove 1/2 quantity during a timestep,
        !    reduce to first order flux. This will remain positive-definite
        !    under the assumption that ABS(CFL(i)) + ABS(CFL(i-1)) < 1.0 if
        !    both fluxes are evacuating the box.

        do k = 1,m1-1
           if (vt3dc(k,i,j) .gt. 0.) then
              if (vt3dg(k,i,j) * vt3dk(k,i,j) .gt.  &
                   dfact * scr1(k,i,j)) then
                 vt3dg(k,i,j) = vt3dc(k,i,j) * scr1(k,i,j)
              endif
           elseif (vt3dc(k,i,j) .lt. 0.) then
              if (-vt3dg(k,i,j) * vt3dk(k+1,i,j) .gt.  &
                   dfact * scr1(k+1,i,j)) then
                 vt3dg(k,i,j) = vt3dc(k,i,j) * scr1(k+1,i,j)
              endif
           endif
        enddo
     enddo
  enddo

  ! Compute flux divergence

  do j = ja,jz
     do i = ia,iz
        do k = 2,m1-1
           scr1(k,i,j) = scr1(k,i,j)  &
                + vt3dk(k,i,j) * (vt3dg(k-1,i,j) - vt3dg(k,i,j)  &
                + scp(k,i,j) * (vt3dc(k,i,j) - vt3dc(k-1,i,j)))
        enddo
     enddo
  enddo

  return
end subroutine fa_zc

!     ****************************************************************

subroutine advtndc(m1,m2,m3,ia,iz,ja,jz,scp,sca,sct,dtl,mynum)
!> @brief: advtndc
!! @date:  18/Nov/2015
!! @version:  5.2
!! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
!! 
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
  implicit none
  integer :: m1,m2,m3,ia,iz,ja,jz,mynum,i,j,k

  real :: dtl,dtli
  real, dimension(m1,m2,m3) :: scp,sca,sct

  dtli = 1. / dtl
  do j = ja,jz
     do i = ia,iz
        do k = 2,m1-1
           sct(k,i,j) = sct(k,i,j) + (sca(k,i,j)-scp(k,i,j)) * dtli
        enddo
     enddo
  enddo

  return
end subroutine advtndc
