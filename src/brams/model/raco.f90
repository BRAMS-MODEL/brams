!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
!  Modifications done for Runge-Kutta dynamical core
!  Dec/2015 by Saulo Freitas (INPE), Michael Baldauf (DWD)
!###########################################################################
module ModAcoust

  !- divergence damping coefficient [m^2/s];
  !- defined as a 3d field to optionally reduce it over steep orography
  real, private, allocatable :: alpha_div(:,:,:)


  real, private :: div_damp_strength    ! dim.less coeff for div-damping
  real :: div_damp_slope

  !- .TRUE. is necessary at least for dyncore_flags 2 and 3
  !- currently div damping is not applied to leapfrog solver
  logical,parameter :: apply_div_damping =.TRUE.
  integer,parameter :: damp_formulation = 1 !1=orig, 2=wrf
contains



  subroutine coefz(mzp,mxp,myp,ia,iz,ja,jz,  &
       acoc,acof,acog,dn0,pi0,th0,rtgt,a1da2,amoe,amof,acoaa)
!> @brief: Calculate coefficients for the vertical pressure gradient
!!  and divergence terms.  These will be combined later for the
!!  implicit computations.
!! @author:  unknow
!! @date:  18/Nov/2015
!! @version:  5.2
!! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
!!
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
    use mem_grid, only: impl, dts, sspct, dzm, dzt, nzp, nz
    use mem_scratch, only: vctr11, vctr12
    use rconstants, only: rgas, cv

    implicit none
    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    real, intent(out) :: acoc(mzp,mxp,myp)
    real, intent(out) :: acof(mzp,mxp,myp)
    real, intent(out) :: acog(mzp,mxp,myp)
    real, intent(in) :: dn0(mzp,mxp,myp)
    real, intent(in) :: pi0(mzp,mxp,myp)
    real, intent(in) :: th0(mzp,mxp,myp)
    real, intent(in) :: rtgt(mxp,myp)
    real, intent(out) :: a1da2
    real, intent(out) :: amoe(mzp,mxp,myp)
    real, intent(out) :: amof(mzp,mxp,myp)
    real, intent(out) :: acoaa(mzp,mxp,myp)



    integer :: i,j,k
    real :: dt2al2,rdto2cv,dt2al2r,rdtr
    real :: acobb(mzp)
    real :: acocc(mzp)

    !acof=0.0 !LFR because acof(k+1,:,:) isn't made
    if (impl .eq. 1) then
       dt2al2 = dts * .75
       a1da2 = 1. / 3.
    else
       dt2al2 = dts
       a1da2 = 1.
    endif
    rdto2cv = sspct ** 2 * rgas * dts / (2.0 * cv)

    do j = ja,jz
       do i = ia,iz

          !         Coefficient for the vertical pressure gradient term

          dt2al2r = .5 * dt2al2 / rtgt(i,j)
          do k = 1,mzp-1
             acoc(k,i,j) = dt2al2r * dzm(k) * (th0(k,i,j) + th0(k+1,i,j))
          enddo

          !         Coefficients for the vertical divergence term

          rdtr = rdto2cv / rtgt(i,j)
          do k = 2,mzp
             vctr12(k) = dn0(k,i,j) * th0(k,i,j)
             vctr11(k) = rdtr * pi0(k,i,j) * dzt(k) / vctr12(k)
          enddo
          vctr12(1) = dn0(1,i,j) * th0(1,i,j)
          do k = 2,mzp-1
             acof(k,i,j) = -vctr11(k) * (vctr12(k) + vctr12(k+1))
             acog(k,i,j) = vctr11(k) * (vctr12(k) + vctr12(k-1))
          enddo
          acog(mzp,i,j) = vctr11(nzp) * (vctr12(nzp) + vctr12(nz))

!srf - bug-fix 04092012
	  acof(mzp,i,j) = -vctr11(nzp) * (vctr12(nzp) + vctr12(nz))

          do k = 2,mzp-1
             acoaa(k,i,j) = acoc(k,i,j) * acog(k,i,j)
             acobb(k) = acoc(k,i,j) * (acof(k,i,j) - acog(k+1,i,j)) - 1.
             acocc(k) = -acoc(k,i,j) * acof(k+1,i,j)
          enddo
          acobb(1) = -1.
          acocc(1) = 0.
          acoaa(mzp,i,j) = 0.
          acobb(mzp) = -1.
          acocc(mzp) = 0.

          amof(1,i,j) = acobb(1)
          amoe(1,i,j) = acocc(1) / amof(1,i,j)
          do k = 2,mzp
             amof(k,i,j) = acobb(k) - acoaa(k,i,j) * amoe(k-1,i,j)
             amoe(k,i,j) = acocc(k) / amof(k,i,j)
          enddo

       enddo
    enddo
    return
  end subroutine coefz


  subroutine rayf_u(mzp,mxp,myp,ia,izu,ja,jz,up,rtgx,topx)
!> @brief: This routine calculates rayleigh friction terms velocity and theta_il
!! @author:  unknow
!! @date:  18/Nov/2015
!! @version:  5.2
!! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
!!
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
    use mem_grid, only : nfpt, distim, nnz, zmn, ztop, dts, nzp, zt, ngrid, nz
    use mem_scratch, only : vctr2, vctr5
    use ref_sounding, only : u01dn

    implicit none
    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: ia
    integer, intent(in) :: izu
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    real, intent(inout) :: up(mzp,mxp,myp)
    real, intent(in) :: rtgx(mxp,myp)
    real, intent(in) :: topx(mxp,myp)

    real :: zmkf,c1,c2
    integer :: kf,i,j,k

    !

    if (nfpt /= 0 .and. distim > 0) then
       kf = nnz(1) - nfpt
       zmkf = zmn(kf,1)
       c1 = 1. / (distim * (ztop - zmkf))
       c2 = dts * c1

       ! u friction

       do j = ja,jz
          do i = ia,izu
             do k = 1,nzp
                vctr2(k) = zt(k) * rtgx(i,j) + topx(i,j)
             enddo
             call htint(nzp,u01dn(1,ngrid),zt,nzp,vctr5,vctr2)
             do k = nz,2,-1
                if (vctr2(k) .le. zmkf) exit
                up(k,i,j) = up(k,i,j) + c2 * (vctr2(k) - zmkf)  &
                     * (vctr5(k) - up(k,i,j))
             enddo
          enddo
       enddo
    end if
  end subroutine rayf_u



  subroutine prdctu(mzp,mxp,myp,ia,izu,ja,jz,  &
       up,ut,pp,th0,f13u,rtgu,rtgt,dxu,topu)
!> @brief: U prediction
!! @author:  unknow
!! @date:  18/Nov/2015
!! @version:  5.2
!! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
!!
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt

    use mem_grid, only : hw4, dzt, distim, dts, nstbot, itopo

    implicit none
    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: ia
    integer, intent(in) :: izu
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    real, intent(inout) :: up(mzp,mxp,myp)
    real, intent(in) :: ut(mzp,mxp,myp)
    real, intent(in) :: pp(mzp,mxp,myp)
    real, intent(in) :: th0(mzp,mxp,myp)
    real, intent(in) :: f13u(mxp,myp)
    real, intent(in) :: rtgu(mxp,myp)
    real, intent(in) :: rtgt(mxp,myp)
    real, intent(in) :: dxu(mxp,myp)
    real, intent(in) :: topu(mxp,myp)

    integer :: i,j,k
    real :: dxl
    real :: vt3da(mzp,mxp,myp)
    real :: dpdx(mzp,mxp,myp)

    !

    dpdx = 0.

    !     Calculate acoustic tendency (horizontal pressure gradient)

    do j = ja,jz
       do i = ia,izu
          do k = 1,mzp-1
             vt3da(k,i,j) = (pp(k,i,j) + pp(k+1,i,j)  &
                  + pp(k,i+1,j) + pp(k+1,i+1,j)) * hw4(k)
          enddo
       enddo
    enddo

    do j = ja,jz
       do i = ia,izu
          dxl = dxu(i,j) / rtgu(i,j)
          do k = 2,mzp-1
             dpdx(k,i,j) = -(th0(k,i,j) + th0(k,i+1,j)) * .5  &
                  * ((pp(k,i+1,j) * rtgt(i+1,j) - pp(k,i,j) * rtgt(i,j)) * dxl  &
                  + (vt3da(k,i,j) - vt3da(k-1,i,j)) * dzt(k) * f13u(i,j))
          enddo
       enddo
    enddo

    if (distim .ne. 0.) then
       call rayf_u(mzp,mxp,myp,ia,izu,ja,jz,up,rtgu,topu)
    endif

    do j = 1,myp
       do i = 1,mxp
          do k = 1,mzp
             up(k,i,j) = up(k,i,j) + dts * (dpdx(k,i,j) + ut(k,i,j))
          enddo
       enddo
    enddo

    if (nstbot .eq. 1 .and. itopo .eq. 1) then
       do i = 1,mxp
          do j = 1,myp
             up(1,i,j) = up(2,i,j)
          enddo
       enddo
    end if
  end subroutine prdctu




  subroutine rayf_v(mzp,mxp,myp,ia,iz,ja,jz,vp,rtgx,topx)
!> @brief:  This routine calculates rayleigh friction terms velocity and theta_il
!! @author:  unknow
!! @date:  18/Nov/2015
!! @version:  5.2
!! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
!!
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
    use mem_grid, only : nfpt, distim, nnz, zmn, ztop, &
         dts, nzp, zt, ngrid, nz, jdim, icorflg
    use mem_scratch, only : vctr2, vctr5
    use ref_sounding, only: v01dn

    implicit none
    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    real, intent(inout) :: vp(mzp,mxp,myp)
    real, intent(in) :: rtgx(mxp,myp)
    real, intent(in) :: topx(mxp,myp)

    real :: zmkf,c1,c2
    integer :: kf,i,j,k

    if (nfpt /= 0 .and. distim > 0 .and. (jdim /= 0 .or. icorflg /= 0)) then
       kf = nnz(1) - nfpt
       zmkf = zmn(kf,1)
       c1 = 1. / (distim * (ztop - zmkf))
       c2 = dts * c1

       ! v friction

       do j = ja,jz
          do i = ia,iz
             do k = 1,nzp
                vctr2(k) = zt(k) * rtgx(i,j) + topx(i,j)
             enddo
             call htint(nzp,v01dn(1,ngrid),zt,nzp,vctr5,vctr2)
             do k = nz,2,-1
                if (vctr2(k) .le. zmkf) exit
                vp(k,i,j) = vp(k,i,j) + c2 * (vctr2(k) - zmkf)  &
                     * (vctr5(k) - vp(k,i,j))
             enddo
          enddo
       enddo
    end if
  end subroutine rayf_v





  subroutine prdctv(mzp,mxp,myp,ia,iz,ja,jzv,  &
       vp,vt,pp,th0,f23v,rtgv,rtgt,dyv,topv)
!> @brief: V prediction
!! @author:  unknow
!! @date:  18/Nov/2015
!! @version:  5.2
!! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
!!
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
    use mem_grid, only : hw4, dzt, distim, dts, nstbot, itopo, jdim

    implicit none
    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jzv
    real, intent(inout) :: vp(mzp,mxp,myp)
    real, intent(in) :: vt(mzp,mxp,myp)
    real, intent(in) :: pp(mzp,mxp,myp)
    real, intent(in) :: th0(mzp,mxp,myp)
    real, intent(in) :: f23v(mxp,myp)
    real, intent(in) :: rtgv(mxp,myp)
    real, intent(in) :: rtgt(mxp,myp)
    real, intent(in) :: dyv(mxp,myp)
    real, intent(in) :: topv(mxp,myp)

    integer :: i,j,k
    real :: dyl
    real :: vt3da(mzp,mxp,myp)
    real :: dpdy(mzp,mxp,myp)

    !     V prediction

    dpdy = 0.

    if (jdim .eq. 1) then

       ! calculate acoustic tendency (horizontal pressure gradient)

       do j = ja,jzv
          do i = ia,iz
             do k = 1,mzp-1
                vt3da(k,i,j) =(pp(k,i,j) + pp(k+1,i,j)  &
                     + pp(k,i,j+1) + pp(k+1,i,j+1)) * hw4(k)
             enddo
          enddo
       enddo

       do j = ja,jzv
          do i = ia,iz
             dyl = dyv(i,j) / rtgv(i,j)
             do k = 2,mzp-1
                dpdy(k,i,j) = -(th0(k,i,j) + th0(k,i,j+1)) * .5  &
                     * ((pp(k,i,j+1)*rtgt(i,j+1) - pp(k,i,j)*rtgt(i,j)) * dyl  &
                     + (vt3da(k,i,j) - vt3da(k-1,i,j)) * dzt(k) * f23v(i,j))
             enddo
          enddo
       enddo
    endif

    if (distim .ne. 0.) then
       call rayf_v(mzp,mxp,myp,ia,iz,ja,jzv,vp,rtgv,topv)
    endif

    do j = 1,myp
       do i = 1,mxp
          do k = 1,mzp
             vp(k,i,j) = vp(k,i,j) + dts * (dpdy(k,i,j) + vt(k,i,j))
          enddo
       enddo
    enddo

    if (nstbot .eq. 1 .and. itopo .eq. 1) then
       do i = 1,mxp
          do j = 1,myp
             vp(1,i,j) = vp(2,i,j)
          enddo
       enddo
    end if
  end subroutine prdctv




  subroutine rayf_w(mzp,mxp,myp,ia,iz,ja,jz,wp,rtgx,topx)
!> @brief:    This routine calculates rayleigh friction terms velocity and theta_il
!! @author:  unknow
!! @date:  18/Nov/2015
!! @version:  5.2
!! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
!!
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
    use mem_grid, only : nfpt, distim, nnz, zmn, ztop, dts, zt, nz
    use mem_scratch, only : vctr2

    implicit none
    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    real, intent(inout) :: wp(mzp,mxp,myp)
    real, intent(in) :: rtgx(mxp,myp)
    real, intent(in) :: topx(mxp,myp)

    real :: zmkf,c1,c2
    integer :: kf,i,j,k

    if (nfpt /= 0 .and. distim > 0) then
       kf = nnz(1) - nfpt
       zmkf = zmn(kf,1)
       c1 = 1. / (distim * (ztop - zmkf))
       c2 = dts * c1

       !     W friction

       do j = ja,jz
          do i = ia,iz
             do k = nz,2,-1
                vctr2(k) = zt(k) * rtgx(i,j) + topx(i,j)
                if (vctr2(k) .le. zmkf) exit
                wp(k,i,j) = wp(k,i,j) - c2 * (vctr2(k) - zmkf) * wp(k,i,j)
             enddo
          enddo
       enddo
    end if
  end subroutine rayf_w





  subroutine prdctw1(mzp,mxp,myp,ia,iz,ja,jz,  &
       wp,wt,pp,acoc,a1da2,rtgt,topt)
!> @brief: W prediction
!! @author:  unknow
!! @date:  18/Nov/2015
!! @version:  5.2
!! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
!!
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
    use mem_grid, only : distim, dts
    use node_mod, only :  mynum

    implicit none
    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    real, intent(inout) :: wp(mzp,mxp,myp)
    real, intent(in) :: wt(mzp,mxp,myp)
    real, intent(in) :: pp(mzp,mxp,myp)
    real, intent(in) :: acoc(mzp,mxp,myp)
    real, intent(in) :: a1da2
    real, intent(in) :: rtgt(mxp,myp)
    real, intent(in) :: topt(mxp,myp)

    integer :: i,j,k

    !     First part of prediction at I,J point

    !     Compute forward part of Crank-Nickelson scheme. This will be total
    !     W prediction for explicit case.

    if (distim .ne. 0.) then
       call rayf_w(mzp,mxp,myp,ia,iz,ja,jz,wp,rtgt,topt)
    endif

    do j = 1,myp
       do i = 1,mxp
          do k = 1,mzp-2
             wp(k,i,j) = wp(k,i,j) + dts * wt(k,i,j)
          enddo
       enddo
    enddo

    do j = ja,jz
       do i = ia,iz
          do k = 1,mzp-2
             wp(k,i,j) = wp(k,i,j) + &
                  a1da2 * acoc(k,i,j) * (pp(k,i,j)-pp(k+1,i,j))
          enddo
       enddo
    enddo
!   write (*,fmt='("wp-pr4: ",I2.2,1X,I1.1,1X,4(E9.3,1X),"/",4(E9.3,1X))') mynum,1,(wp(k,2,2),k=1,4),(wp(k,3,3),k=1,4)
!  call flush(6)

  end subroutine prdctw1


  subroutine prdctp1(mzp,mxp,myp,ia,iz,ja,jz,  &
       pp,up,vp,pi0,dn0,th0,pt,f13t,f23t,rtgt,rtgu,rtgv,  &
       heatfx1,fmapui,fmapvi,dxt,dyt,fmapt)
!> @brief: Divergence calculations for topographical transformation
!! @author:  unknow
!! @date:  18/Nov/2015
!! @version:  5.2
!! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
!!
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
    use mem_grid, only : sspct, & !intent(in)
         jdim,                  & !intent(in)
         hw4,                   & !intent(in)
         dzt,                   & !intent(in)
         dts                      !intent(in)
    use rconstants, only : rocv   !intent(in)

    implicit none
    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    real, intent(inout) :: pp(mzp,mxp,myp)
    real, intent(in) :: up(mzp,mxp,myp)
    real, intent(in) :: vp(mzp,mxp,myp)
    real, intent(in) :: pi0(mzp,mxp,myp)
    real, intent(in) :: dn0(mzp,mxp,myp)
    real, intent(in) :: th0(mzp,mxp,myp)
    real, intent(in) :: pt(mzp,mxp,myp)
    real, intent(in) :: f13t(mxp,myp)
    real, intent(in) :: f23t(mxp,myp)
    real, intent(in) :: rtgt(mxp,myp)
    real, intent(in) :: rtgu(mxp,myp)
    real, intent(in) :: rtgv(mxp,myp)
    real, intent(out) :: heatfx1(mxp,myp)
    real, intent(in) :: fmapui(mxp,myp)
    real, intent(in) :: fmapvi(mxp,myp)
    real, intent(in) :: dxt(mxp,myp)
    real, intent(in) :: dyt(mxp,myp)
    real, intent(in) :: fmapt(mxp,myp)

    ! Local Variables
    integer :: i, j, k
    real :: rocvpct
    real :: heatdv(mzp,mxp,myp)
    real :: heatfx(mzp,mxp,myp)

    heatdv = 0.
    rocvpct =rocv *sspct ** 2

    !     Divergence calculations for topographical transformation

    !     First calculate vertically transformed heat flux

    do j = ja,jz
       do i = ia,iz
          do k = 1,mzp
             heatfx(k,i,j) = ((up(k,i,j) + up(k,i-1,j)) * f13t(i,j)  &
                  + (vp(k,i,j) + vp(k,i,j-jdim)) * f23t(i,j)  &
                  ) * dn0(k,i,j) * th0(k,i,j)
          enddo
       enddo
    enddo
    do j = ja,jz
       do i = ia,iz
          do k = 1,mzp-1
             heatfx(k,i,j) = (heatfx(k,i,j) + heatfx(k+1,i,j)) * hw4(k)
          enddo
          heatfx1(i,j) = heatfx(1,i,j) / (.5 * (dn0(1,i,j) * th0(1,i,j)  &
               + dn0(2,i,j) * th0(2,i,j)))
       enddo
    enddo

    do j = ja,jz
       do i = ia,iz
          do k = 2,mzp-1
             heatdv(k,i,j) = (heatfx(k,i,j) - heatfx(k-1,i,j)) * dzt(k)
          enddo
       enddo
    enddo
    do j = 1,myp
       do i = 1,mxp
          do k = 1,mzp
             heatfx(k,i,j) = dn0(k,i,j) * th0(k,i,j)
          enddo
       enddo
    enddo
    do j = ja,jz
       do i = ia,iz
          do k = 2,mzp-1

             heatdv(k,i,j) = -rocvpct * pi0(k,i,j) / heatfx(k,i,j)  &
                  * (heatdv(k,i,j) + fmapt(i,j)  &
                  * ((up(k,i,j) * rtgu(i,j) * fmapui(i,j)  &
                  * (heatfx(k,i,j) + heatfx(k,i+1,j))  &
                  - up(k,i-1,j) * rtgu(i-1,j) * fmapui(i-1,j)  &
                  * (heatfx(k,i,j) + heatfx(k,i-1,j))) * dxt(i,j) * .5  &
                  + (vp(k,i,j) * rtgv(i,j) * fmapvi(i,j)  &
                  * (heatfx(k,i,j) + heatfx(k,i,j+jdim))  &
                  - vp(k,i,j-jdim) * rtgv(i,j-jdim)  &
                  * fmapvi(i,j-jdim)  &
                  * (heatfx(k,i,j) + heatfx(k,i,j-jdim)))  &
                  * dyt(i,j) * .5) / rtgt(i,j))

          enddo
       enddo
    enddo


    do j = ja,jz
       do i = ia,iz
          do k = 1,mzp
             pp(k,i,j) = pp(k,i,j) + (pt(k,i,j) + heatdv(k,i,j)) * dts
          enddo
       enddo
    enddo
  end subroutine prdctp1




  subroutine  prdctw2(mzp,mxp,myp,ia,iz,ja,jz,  &
       wp,pp,acoc,amof,amog,acoaa,rtgt,heatfx1)
!> @brief: W2 prediction
!! @author:  unknow
!! @date:  18/Nov/2015
!! @version:  5.2
!! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
!!
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
    use mem_grid, only: nstbot, nsttop, impl, nz

    implicit none
    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    real, intent(inout) :: wp(mzp,mxp,myp)
    real, intent(in) :: pp(mzp,mxp,myp)
    real, intent(in) :: acoc(mzp,mxp,myp)
    real, intent(in) :: amof(mzp,mxp,myp)
    real, intent(out) :: amog(mzp,mxp,myp)
    real, intent(in) :: acoaa(mzp,mxp,myp)
    real, intent(in) :: rtgt(mxp,myp)
    real, intent(in) :: heatfx1(mxp,myp)

    integer :: i,j,k

    if (nstbot .eq. 1) then
       do j = ja,jz
          do i = ia,iz
             wp(1,i,j) = -heatfx1(i,j) * rtgt(i,j)
          enddo
       enddo
    endif

    if (nsttop .eq. 1) then
       do j = ja,jz
          do i = ia,iz
             wp(nz,i,j) = 0.
          enddo
       enddo
    endif

    if (impl .eq. 1) then

       !  First implicit part of the w prediction

       do j = ja,jz
          do i = ia,iz
             do k = 2,mzp-2
                wp(k,i,j) = wp(k,i,j) - (pp(k+1,i,j) - pp(k,i,j)) * acoc(k,i,j)
             enddo
          enddo
       enddo

       do j = ja,jz
          do i = ia,iz
             amog(1,i,j) = -wp(1,i,j) / amof(1,i,j)
          enddo
          do k = 2,mzp-2
             do i = ia,iz
                amog(k,i,j) = (-wp(k,i,j) - acoaa(k,i,j) * amog(k-1,i,j))  &
                     / amof(k,i,j)
             enddo
          enddo
       enddo
    endif
  end subroutine prdctw2



  subroutine prdctw3(mzp,mxp,myp,ia,iz,ja,jz,wp,amog,amoe,impl)
!> @brief: Conclusion of implicit w prediction
!! @author:  unknow
!! @date:  18/Nov/2015
!! @version:  5.2
!! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
!!
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
    implicit none
    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    real, intent(inout) :: wp(mzp,mxp,myp)
    real, intent(in) :: amog(mzp,mxp,myp)
    real, intent(in) :: amoe(mzp,mxp,myp)
    integer, intent(in) :: impl

    integer :: i,j,k

    if (impl .eq. 1) then
       do k = mzp-2,2,-1
          do j = ja,jz
             do i = ia,iz
                wp(k,i,j) = amog(k,i,j) - amoe(k,i,j) * wp(k+1,i,j)
             enddo
          enddo
       enddo
    endif
  end subroutine prdctw3




  subroutine prdctp2(mzp,mxp,myp,ia,iz,ja,jz,pp,wp,acof,acog)
!> @brief: Finish pressure prediction
!! @author:  unknow
!! @date:  18/Nov/2015
!! @version:  5.2
!! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
!!
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
    use mem_grid, only: nstbot, dzm
    use node_mod, only :  mynum

    implicit none
    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    real, intent(inout) :: pp(mzp,mxp,myp)
    real, intent(in) :: wp(mzp,mxp,myp)
    real, intent(in) :: acof(mzp,mxp,myp)
    real, intent(in) :: acog(mzp,mxp,myp)

    integer :: i,j,k
    real :: dzmr

    !  Finish pressure prediction
!write (*,fmt='("Pr-pre1: ",I2.2,1X,I1.1,1X,4(E9.3,1X),"/",4(E9.3,1X))') mynum,1,(pp(k,2,2),k=1,4),(pp(k,3,3),k=1,4)
!write (*,fmt='("Pr-wp  : ",4(E9.3,1X))') (wp(k,2,2),k=1,4)
!write (*,fmt='("Pr-acof: ",4(E9.3,1X))') (acof(k,2,2),k=1,4)
!write (*,fmt='("Pr-acog: ",4(E9.3,1X))') (acog(k,2,2),k=1,4)

!call flush(6)

    do j = ja,jz
       do i = ia,iz
          do k = 2,mzp-1
             pp(k,i,j) = pp(k,i,j)  &
                  + (wp(k,i,j) * acof(k,i,j) + wp(k-1,i,j) * acog(k,i,j))
          enddo
       enddo
    enddo
!write (*,fmt='("Pr-mid1: ",I2.2,1X,I1.1,1X,4(E9.3,1X),"/",4(E9.3,1X))') mynum,1,(pp(k,2,2),k=1,4),(pp(k,3,3),k=1,4)

    if (nstbot .eq. 1) then
       dzmr = dzm(2) / dzm(1)
       do i = 1,mxp
          do j = 1,myp
             pp(1,i,j) = pp(2,i,j) + (pp(2,i,j) - pp(3,i,j)) * dzmr
          enddo
       enddo
    end if

!write (*,fmt='("Pr-pos1: ",I2.2,1X,I1.1,1X,4(E9.3,1X),"/",4(E9.3,1X))') mynum,1,(pp(k,2,2),k=1,4),(pp(k,3,3),k=1,4)
!call flush(6)

  end subroutine prdctp2




  subroutine acoustic_new(OneGrid, nnacoust_loc,lrk)
    !> @brief: Acoustic terms small time-step driver
    !! @details: This routine calls all the necessary routines to march the model
    !!     through the small timesteps.
    !! @author:  unknow
    !! @date:  18/Nov/2015
    !! @version:  5.2
    !! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
    !!
    !! @copyright Under CC-GPL License by INPE/CPTEC
    !! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
    use dump, only: &
      dumpMessage
    use mem_tend
    use mem_grid, only: dts, grid_g, ngrid, if_adap, dyncore_flag
    use mem_basic
    use mem_scratch
    use node_mod
    use ModGrid, only: &
         Grid
    use ModAcoust_adap, only: acoust_adap

    implicit none
    include "constants.f90"
    character(len=*),parameter :: h='**(acoustic_new)**'
    type(Grid), pointer :: OneGrid
    integer, intent(in) :: nnacoust_loc  ! number of small time steps
    integer, intent(in) :: lrk
    real :: scr2(mzp*mxp*myp)

    integer :: nmbr_gpts  !MB: only for testing

    !print*, "subr. [acoustic_new] ..."

    scr2=0.0

    if ( ( dyncore_flag == 0 ) .or. ( dyncore_flag == 1 ) ) then

       if ( if_adap == 0) then

          call acoust_new(OneGrid, &
               mzp,mxp,myp, &
               nnacoust_loc,  &
               basic_g(ngrid)%dn0,basic_g(ngrid)%pi0,  &
               basic_g(ngrid)%th0,basic_g(ngrid)%up,  &
               basic_g(ngrid)%vp,basic_g(ngrid)%wp,  &
               basic_g(ngrid)%pp,  &
               tend%ut,tend%vt,tend%wt,tend%pt,  &
               grid_g(ngrid)%topt,grid_g(ngrid)%topu,  &
               grid_g(ngrid)%topv,grid_g(ngrid)%rtgt,  &
               grid_g(ngrid)%rtgu,grid_g(ngrid)%f13u,  &
               grid_g(ngrid)%dxu,grid_g(ngrid)%rtgv,  &
               grid_g(ngrid)%dyv,  &
               grid_g(ngrid)%f23v,grid_g(ngrid)%f13t,  &
               grid_g(ngrid)%f23t,grid_g(ngrid)%fmapui,  &
               grid_g(ngrid)%fmapvi,grid_g(ngrid)%dxt,  &
               grid_g(ngrid)%dyt,grid_g(ngrid)%fmapt,0)

       else

          call acoust_adap(OneGrid, &
               mzp,mxp,myp   &
               ,grid_g(ngrid)%lpu      ,grid_g(ngrid)%lpv       &
               ,grid_g(ngrid)%lpw      ,scratch%scr1              &
               ,scr2             ,scratch%vt3da       &
               ,scratch%vt3db        ,scratch%vt3dc       &
               ,scratch%vt3dd        ,scratch%vt3de       &
               ,scratch%vt3df        ,scratch%vt3dg       &
               ,scratch%vt3dh        ,scratch%vt2da       &
               ,basic_g(ngrid)%dn0   ,basic_g(ngrid)%pi0  &
               ,basic_g(ngrid)%th0   ,basic_g(ngrid)%up   &
               ,basic_g(ngrid)%vp    ,basic_g(ngrid)%wp   &
               ,basic_g(ngrid)%pp    ,tend%ut	       &
               ,tend%vt              ,tend%wt	       &
               ,tend%pt              ,grid_g(ngrid)%dxu   &
               ,grid_g(ngrid)%dyv    ,grid_g(ngrid)%fmapu &
               ,grid_g(ngrid)%fmapvi ,grid_g(ngrid)%dxt   &
               ,grid_g(ngrid)%dyt    ,grid_g(ngrid)%fmapt &
               ,grid_g(ngrid)%aru    ,grid_g(ngrid)%arv   &
               ,grid_g(ngrid)%arw    ,grid_g(ngrid)%volt  &
               ,grid_g(ngrid)%volu   ,grid_g(ngrid)%volv  &
               ,grid_g(ngrid)%volw  		       )

       endif

    else if ( dyncore_flag == 2 .or. dyncore_flag == 3) then
       ! Runge-Kutta based solver
       ! remark: the philosphy of the time levels is slightly different
       ! than in leapfrog: here always fields uc, vc, .. are updated

       !call acoust_new_rk(OneGrid, &
       call acoust_new(OneGrid, &
            mzp,mxp,myp, &
            nnacoust_loc,   &
            basic_g(ngrid)%dn0, basic_g(ngrid)%pi0,  &
            basic_g(ngrid)%th0,                      &
	    basic_g(ngrid)%uc, basic_g(ngrid)%vc,    &
	    basic_g(ngrid)%wc, basic_g(ngrid)%pc,    &
            tend%ut_rk, tend%vt_rk,            &
	    tend%wt_rk, tend%pt_rk,            &
            grid_g(ngrid)%topt, grid_g(ngrid)%topu,  &
            grid_g(ngrid)%topv, grid_g(ngrid)%rtgt,  &
            grid_g(ngrid)%rtgu, grid_g(ngrid)%f13u,  &
            grid_g(ngrid)%dxu,  grid_g(ngrid)%rtgv,  &
            grid_g(ngrid)%dyv,  &
            grid_g(ngrid)%f23v, grid_g(ngrid)%f13t,  &
            grid_g(ngrid)%f23t,grid_g(ngrid)%fmapui,  &
            grid_g(ngrid)%fmapvi,grid_g(ngrid)%dxt,  &
            grid_g(ngrid)%dyt,grid_g(ngrid)%fmapt,lrk)

     !MB:
     !nmbr_gpts = mxp * myp * mzp    !MB: only for testing!!!
     !write(*,"(A,I2,4F15.8)") "u    : ", 78, minval( basic_g(ngrid)%uc    ), maxval( basic_g(ngrid)%uc    ), &
     !    sum( basic_g(ngrid)%uc    )/nmbr_gpts, basic_g(ngrid)%uc(7,8,9)
     !write(*,"(A,I2,4F15.8)") "v    : ", 78, minval( basic_g(ngrid)%vc    ), maxval( basic_g(ngrid)%vc    ), &
     !    sum( basic_g(ngrid)%vc    )/nmbr_gpts, basic_g(ngrid)%vc(7,8,9)
     !write(*,"(A,I2,4F15.8)") "w    : ", 78, minval( basic_g(ngrid)%wc    ), maxval( basic_g(ngrid)%wc    ), &
     !    sum( basic_g(ngrid)%wc    )/nmbr_gpts, basic_g(ngrid)%wc(7,8,9)
     !write(*,"(A,I2,4F15.8)") "pi   : ", 78, minval( basic_g(ngrid)%pc    ), maxval( basic_g(ngrid)%pc    ), &
     !    sum( basic_g(ngrid)%pc    )/nmbr_gpts, basic_g(ngrid)%pc(7,8,9)

    else
       !call fatal_error("ERROR in acoustic_new: value of dyncore_flag is unknown!")
       iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
         "ERROR in acoustic_new: value of dyncore_flag is unknown!")
    end if

    return
  end subroutine acoustic_new



  subroutine acoust_new(OneGrid, mzp,mxp,myp,  &
       nnacoust_loc,                           &
       dn0,pi0,th0,                            &
       up,vp,wp,pp,                            &
       ut,vt,wt,pt,                            &
       topt,topu,topv,rtgt,rtgu,f13u,dxu,rtgv, &
       dyv,f23v,f13t,f23t,fmapui,fmapvi,dxt,dyt,fmapt,lrk)
!> @brief: Acoustic terms small time-step driver
!! @details: This routine calls all the necessary routines to march the model
!!     through the nnacoust_loc small timesteps.
!!     The integration starts form the state up, vp, wp, pp; these fields
!!     are overwritten with the updated fields at output.
!! @author:  unknow
!! @date:  18/Nov/2015
!! @version:  5.2
!! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
!!
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt


    use ModParallelEnvironment, only: MsgDump

    use ModNamelistFile, only: namelistfile

    use ModGrid, only: Grid, DumpGrid

    use ModMessageSet, only: &
         PostRecvSendMsgs, &
         WaitRecvMsgs

    use mem_grid, only : nnacoust, & ! intent(in)
         ngrid,                    & ! intent(in)
         dts,                      & ! intent(out)
         dtlt,                     & ! intent(in)
         nxtnest,                  & ! intent(in)
         nzp, nxp, nyp,            & ! intent(in)
         impl,                     & ! intent(in)
         time,                     &
	 dyncore_flag,             &
         grid_g                      ! for get_wind_div

    use node_mod, only :  nmachs,  & ! intent(in)
         ia,                       & ! intent(in)
         iz,                       & ! intent(in)
         ja,                       & ! intent(in)
         jz,                       & ! intent(in)
         izu,                      & ! intent(in)
         jzv,                      & ! intent(in)
         mynum

    use ModNamelistFile, only: &
         namelistfile

    use mem_scratch, only: scratch

    use ModComm, only: commHaloAcou

    implicit none

    type(Grid), pointer :: OneGrid
    integer, intent(in) :: lrk

    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp

    integer, intent(in) :: nnacoust_loc  ! number of small time steps

    real, intent(in) :: dn0(mzp,mxp,myp)
    real, intent(in) :: pi0(mzp,mxp,myp)
    real, intent(in) :: th0(mzp,mxp,myp)

    real, intent(inout) :: up(mzp,mxp,myp)
    real, intent(inout) :: vp(mzp,mxp,myp)
    real, intent(inout) :: wp(mzp,mxp,myp)
    real, intent(inout) :: pp(mzp,mxp,myp)

    real, intent(in) :: topt(mxp,myp)
    real, intent(in) :: topu(mxp,myp)
    real, intent(in) :: topv(mxp,myp)
    real, intent(in) :: rtgt(mxp,myp)
    real, intent(in) :: rtgu(mxp,myp)
    real, intent(in) :: f13u(mxp,myp)
    real, intent(in) :: dxu(mxp,myp)
    real, intent(in) :: rtgv(mxp,myp)
    real, intent(in) :: dyv(mxp,myp)
    real, intent(in) :: f23v(mxp,myp)
    real, intent(in) :: f13t(mxp,myp)
    real, intent(in) :: f23t(mxp,myp)
    real, intent(in) :: fmapui(mxp,myp)
    real, intent(in) :: fmapvi(mxp,myp)
    real, intent(in) :: dxt(mxp,myp)
    real, intent(in) :: dyt(mxp,myp)
    real, intent(in) :: fmapt(mxp,myp)

    real, intent(in) :: ut(mzp,mxp,myp)
    real, intent(in) :: vt(mzp,mxp,myp)
    real, intent(in) :: wt(mzp,mxp,myp)
    real, intent(in) :: pt(mzp,mxp,myp)

    include "tsNames.h"

    integer :: iter,k
    integer :: lastIter
    logical :: singleProcRun
    logical :: outermostGrid
    real :: a1da2
    real :: acoaa(mzp,mxp,myp)
    real :: acoc(mzp,mxp,myp)
    real :: acof(mzp,mxp,myp)
    real :: acog(mzp,mxp,myp)
    real :: amoe(mzp,mxp,myp)
    real :: amof(mzp,mxp,myp)
    real :: amog(mzp,mxp,myp)
    real :: heatfx1(mxp,myp)

    integer :: i, j
    character(len=1) :: citer

    real, dimension(:,:,:), allocatable :: &
       div,        &
       pp_minus_div&
      ,pp_t_minus_dt
    character(len=*), parameter :: h="**(acoust_new)**"
    logical, parameter :: dumpLocal=.false.
    character(LEN=5) :: ctime
    integer :: nmbr_gpts
    real :: dtacum

    lastIter = nnacoust_loc
    singleProcRun = nmachs == 1
    outermostGrid = nxtnest(ngrid) == 0

    if ( apply_div_damping .and. (dyncore_flag == 2 .or. dyncore_flag == 3) ) then
      allocate( pp_minus_div ( mzp, mxp, myp) )
      if(damp_formulation==1) allocate( div          ( mzp, mxp, myp) )
      if(damp_formulation==2) allocate( pp_t_minus_dt( mzp, mxp, myp) )
    end if

    ! computes acoc, acof, acog, a1da2, amoe, amof, acoaa
    ! uses dn0, pi0, th0, rtgt
    call coefz(mzp,mxp,myp,ia,iz,ja,jz,  &
         acoc,acof,acog,dn0,pi0,th0,rtgt,a1da2,amoe,amof,acoaa)

    !dtacum=0.
    !
    ! run small time steps
    do iter=1,lastIter

       !dtacum=dtacum+dts
       !print*,"NAC=",iter,dts
       !write(citer,fmt='(I1)') iter
       ! if not first loop iteration,
       !   receives pp(i+1,j) and pp(i,j+1)
       !   if outermost grid,
       !     receives pp on full grid boundaries
      if (iter > 1 .and. .not. singleProcRun) then
        call WaitRecvMsgs(OneGrid%AcouSendP, OneGrid%AcouRecvP)
      end if

      if ( apply_div_damping .and. ( dyncore_flag == 2 .or. dyncore_flag == 3 )&
        .and. ( iter > 1 ) ) then
        if(damp_formulation==1) then
          ! application of divergence damping in the 1st timestep is detrimental to stability
          call get_wind_div_v2 (mzp,mxp,myp,ia,iz,ja,jz,izu,jzv    &
            ,up(1,1,1),vp(1,1,1),wp(1,1,1)                   &
            ,scratch%vt3da(1), scratch%vt3db(1), scratch%vt3dh(1) &
            ,div(1,1,1)          &
            ,rtgt(1,1),rtgu(1,1),dxu(1,1),rtgv(1,1),dyv(1,1),f13t(1,1)&
            ,f23t(1,1),fmapui(1,1),fmapvi(1,1)   &
            ,grid_g(ngrid)%fmapu(1,1), grid_g(ngrid)%fmapv(1,1)&
            ,grid_g(ngrid)%dxt,grid_g(ngrid)%dyt,fmapt)
            if (.not. singleProcRun) then
              call commHaloAcou(div,mxp,myp,mzp,myNum,'div')
            endif
            ! as proposed in Wicker, Skamarock (2002) (?) divergence damping is
            ! used in an approximated form by adding the following term to the pressure:
            pp_minus_div(:,:,:) = pp(:,:,:) - alpha_div(:,:,:) / th0(:,:,:) * div(:,:,:)
!      call dumpVarAllLatLonk(div,'aDiv'  //crk//citer,1,mxp,1,myp,1,mzp,0.0,0.0)
        elseif(damp_formulation==2) then
          !- alternative formulation
          if (.not. singleProcRun) then
            call commHaloAcou(pp_t_minus_dt,mxp,myp,mzp,myNum,'pp_t_minus_dt')
          endif
          pp_minus_div(:,:,:) = pp(:,:,:) + div_damp_strength*(pp(:,:,:) - pp_t_minus_dt(:,:,:))
        endif
      end if

      if(damp_formulation==2) pp_t_minus_dt(:,:,:)=pp(:,:,:)

      ! up = f(up, pp); uses pp(i+1,j)
      ! remaining arguments are loop invariant
      if ( apply_div_damping .AND. ( dyncore_flag == 2 .or. dyncore_flag == 3 ) &
        .AND. (iter > 1) )then
        call prdctu(mzp,mxp,myp,ia,izu,ja,jz,  &
            up,ut,pp_minus_div,th0,f13u,rtgu,rtgt,dxu,topu)
      else
        call prdctu(mzp,mxp,myp,ia,izu,ja,jz,  &
              up,ut,pp,th0,f13u,rtgu,rtgt,dxu,topu)
      end if
 !call dumpVarAllLatLonk(up, 'UP'  ,1155,lrk,iter,1,mxp,1,myp,1,mzp,0.0,0.0)
      ! if not last loop iteration,  sends up(i-1,j)
      if (iter < lastIter) then
        call PostRecvSendMsgs(OneGrid%AcouSendU, OneGrid%AcouRecvU)
      end if

      ! vp = f(vp, pp); uses pp(i,j+1)
      ! remaining arguments are loop invariant
      if ( apply_div_damping .AND. ( dyncore_flag == 2 .or. dyncore_flag == 3 ) &
          .AND. (iter > 1) )then
        call prdctv(mzp,mxp,myp,ia,iz,ja,jzv,&
               vp,vt, pp_minus_div, th0,f23v,rtgv,rtgt,dyv,topv)
      else
        call prdctv(mzp,mxp,myp,ia,iz,ja,jzv,&
               vp,vt, pp,th0,f23v,rtgv,rtgt,dyv,topv)
      end if


      ! on parallel runs,
      !   if not last loop iteration
      !      sends vp(i,j-1)
      !   else if last loop iteration,
      !      sends up, vp to update all neighbour processes ghost zone
      if (.not. singleProcRun) then
!        call commHaloAcou(up,mxp,myp,mzp,myNum,'up')
        if (iter < lastIter) then
          call PostRecvSendMsgs(OneGrid%AcouSendV, OneGrid%AcouRecvV)
        else
          call PostRecvSendMsgs(OneGrid%AcouSendUV, OneGrid%AcouRecvUV)
        end if
      end if

      ! wp = f(wp, pp); uses node inner cells only
      ! remaining arguments are loop invariant
      call prdctw1(mzp,mxp,myp,ia,iz,ja,jz, &
            wp,wt,pp,acoc,a1da2,rtgt,topt)

      ! on parallel runs,
      !   if not last loop iteration
      !      receives up(i-1,j)
      !      if outermost grid,
      !        receives up on full grid boundaries
      !      receives vp(i,j-1)
      !      if outermost grid,
      !        receives vp on full grid boundaries
      !   else if last loop iteration,
      !      receives up, vp and updates this process ghost zone
      !      if outermost grid,
      !        receives up, vp on full grid boundaries

      if (.not. singleProcRun) then
        if (iter < lastIter) then
          call WaitRecvMsgs(OneGrid%AcouSendU, OneGrid%AcouRecvU)
          call WaitRecvMsgs(OneGrid%AcouSendV, OneGrid%AcouRecvV)
        else
          call WaitRecvMsgs(OneGrid%AcouSendUV, OneGrid%AcouRecvUV)
        end if
      end if

      ! pp = f(pp, up, vp); uses up(i-1,j), vp(i,j-1)
      ! also computes heatfx1
      ! remaining arguments are loop invariant
      call prdctp1(mzp,mxp,myp,ia,iz,ja,jz,  &
            pp,up,vp,pi0,dn0,th0,pt,f13t,f23t,rtgt,rtgu,rtgv, &
            heatfx1,fmapui,fmapvi,dxt,dyt,fmapt)

      ! wp = f(wp, pp); uses node inner cells only
      ! uses heatfx1
      ! also computes amog
      ! remaining arguments are loop invariant
      call prdctw2(mzp,mxp,myp,ia,iz,ja,jz, &
            wp,pp,acoc,amof,amog,acoaa,rtgt,heatfx1)

      ! finishes updating wp=f(wp); uses node inner cells only
      ! uses amog
      ! remaining arguments are loop invariant
      call prdctw3(mzp,mxp,myp,ia,iz,ja,jz,wp,amog,amoe,impl)
!      call dumpVarAllLatLonk(wp,'aWPg'  //crk//citer,1,mxp,1,myp,1,mzp,0.0,0.0)
      ! finishes updating pp=f(pp,wp); uses node inner cells only
      ! remaining arguments are loop invariant
      call prdctp2(mzp,mxp,myp,ia,iz,ja,jz,pp,wp,acof,acog)

      ! on parallel runs,
      !   if not last loop iteration
      !      sends pp(i,j+1) and pp(i+1,j)
      !   else if last loop iteration,
      !      sends wp, pp to update neighbour processes ghost zone
      !      receives wp and pp and updates this process ghost zone
      if (.not. singleProcRun) then
        if (iter < lastIter) then
          call PostRecvSendMsgs(OneGrid%AcouSendP, OneGrid%AcouRecvP)
        else
          call PostRecvSendMsgs(OneGrid%AcouSendWP, OneGrid%AcouRecvWP)
          call WaitRecvMsgs(OneGrid%AcouSendWP, OneGrid%AcouRecvWP)
        end if
      end if

    enddo

    if ( apply_div_damping .AND. ( dyncore_flag == 2 .or. dyncore_flag == 3) ) then
      if (allocated(div          )) deallocate(div          )
      if (allocated(pp_minus_div )) deallocate(pp_minus_div )
      if (allocated(pp_t_minus_dt)) deallocate(pp_t_minus_dt)
    end if

  end subroutine acoust_new


  subroutine init_div_damping_coeff( dts )

    ! Reduction of the divergence damping coefficient over steep terrain
    ! to increase numerical stability
    ! see sec. 7 in M. Baldauf (2013) COSMO Techn. Rep., nr. 21, www.cosmo-model.org
    use dump, only: &
      dumpMessage

    use mem_grid, only: &
       grid_g,    &
       ngrid,     &
       dyncore_flag, &
       zt, zm

    use node_mod, only:  &
       mxp,  &
       myp,  &
       mzp,  &
       nmachs, &
       mynum

    use ModComm, only: commHaloAcou

    implicit none

    include "constants.f90"
    character(len=*),parameter :: h='**(init_div_damping_coeff)**'
    real, intent(in) :: dts  ! small time step [s]

    real :: c_sound   ! rough estimate for maximum value of sound velocity [m/s]
    real :: alpha_div_limit
    logical :: limit_alpha_div_by_slope_stability

    real :: z_kij, z_mij, z_kpj, z_mpj, z_kip, z_mip
    real :: ax, ay
    real :: work
    real :: delta_h_x_at_u(mxp, myp)
    real :: delta_h_y_at_v(mxp, myp)
    real :: delta_h_x, delta_h_y
    real :: delta_z

    integer :: i, j, k
    integer :: istat
    logical :: singleProcRun

    !print*, "Subr. init_div_damping_coeff ..."
    singleProcRun = nmachs == 1
    limit_alpha_div_by_slope_stability = .true.

    div_damp_strength =  0.1   !MB: might be a namelist-parameter (?)

    div_damp_slope    =  1.0   !MB: this was until now a namelist parameter in COSMO.
                               ! However, with the bug fix done here, it can be
                               ! probably set equal to 1 in accordance with stability theory.

    !MB:
    !if(mynum==1) &
    !write(logUnit,fmt="(A,F12.6,A,F12.6,1X,A,L1)") "div_damp_strength=", div_damp_strength, &
    !        "  div_damp_slope=", div_damp_slope,  &
    !        "limit_alpha_div_by_slope_stability=",limit_alpha_div_by_slope_stability

    allocate( alpha_div( mzp, mxp, myp), STAT=istat )
    if ( istat /= 0 )  &! call fatal_error("ERROR allocating alpha_div")
      iErrNumber=dumpMessage(c_tty,c_yes,h,modelVersion,c_fatal, &
      "ERROR allocating alpha_div")
    c_sound = 330.0   ! this rough estimation of sound velocity is sufficient here
                      ! see Wicker, Skamarock (2002) MWR:
    alpha_div_limit = div_damp_strength * c_sound**2 * dts

    if ( .NOT. limit_alpha_div_by_slope_stability ) then
       alpha_div(:,:,:) = alpha_div_limit
    else
      do k=2, mzp
        do j=1, myp-1
          do i=1, mxp-1
            ! it is assumed here that z(k,i,j) is the height of the grid point
            ! in the position of w(k,i,j)   !!???
            z_kij = zt(k)   * grid_g(ngrid)%rtgt(i  ,j  ) + grid_g(ngrid)%topt(i  ,j  )
            z_mij = zt(k-1) * grid_g(ngrid)%rtgt(i  ,j  ) + grid_g(ngrid)%topt(i  ,j  )
            z_kpj = zt(k)   * grid_g(ngrid)%rtgt(i+1,j  ) + grid_g(ngrid)%topt(i+1,j  )
            z_mpj = zt(k-1) * grid_g(ngrid)%rtgt(i+1,j  ) + grid_g(ngrid)%topt(i+1,j  )
            z_kip = zt(k)   * grid_g(ngrid)%rtgt(i,  j+1) + grid_g(ngrid)%topt(i  ,j+1)
            z_mip = zt(k-1) * grid_g(ngrid)%rtgt(i,  j+1) + grid_g(ngrid)%topt(i  ,j+1)
            !orography steps in x- and y-direction at the grid position of u and v, respectively:
            delta_h_x_at_u(i,j) = 0.5 * ( ( z_kpj + z_mpj )   &
                                            - ( z_kij + z_mij ) )
            delta_h_y_at_v(i,j) = 0.5 * ( ( z_kip + z_mip )   &
                                            - ( z_kij + z_mij ) )
	        end do
        end do
        do j=2, myp-1
          do i=2, mxp-1
            delta_z = zm(k) - zm(k-1)
            delta_h_x = 0.5 * ( abs(delta_h_x_at_u(i,j)) + abs(delta_h_x_at_u(i-1,j  )) )
            delta_h_y = 0.5 * ( abs(delta_h_y_at_v(i,j)) + abs(delta_h_y_at_v(i  ,j-1)) )
            ax = grid_g(ngrid)%dxt(i,j) * (  2.0 + abs( delta_h_x / delta_z ) )
            ay = grid_g(ngrid)%dyt(i,j) * (  2.0 + abs( delta_h_y / delta_z ) )
            work = div_damp_slope * 2.0 / ( dts * ( ax**2 + ay**2 ) )
            alpha_div(k,i,j) = MIN( work, alpha_div_limit )
          end do
        end do
      end do
      ! fill up boundary values (although these shoudn't be needed)
      alpha_div(:,1  ,:  ) = alpha_div(:,2    ,:    )
      alpha_div(:,mxp,:  ) = alpha_div(:,mxp-1,:    )
      alpha_div(:,:  ,1  ) = alpha_div(:,:    ,2    )
      alpha_div(:,:  ,myp) = alpha_div(:,:    ,myp-1)
      alpha_div(:,1  ,1  ) = alpha_div(:,2    ,2    )
      alpha_div(:,mxp,1  ) = alpha_div(:,mxp-1,2    )
      alpha_div(:,1  ,myp) = alpha_div(:,2    ,myp-1)
      alpha_div(:,mxp,myp) = alpha_div(:,mxp-1,myp-1)
      alpha_div(1,:  ,:  ) = alpha_div(2,:    ,:    )
    end if

!--- mpi paralelization :
    !MB: here an exchange of alpha_div for MPI-parallelization is necessary!
    if (.not. singleProcRun) then
      call commHaloAcou(alpha_div,mxp,myp,mzp,myNum,'alpha_div')
    endif
!--- mpi paralelization :

     !print*, "alpha_div: ", minval(alpha_div), maxval(alpha_div) ;call flush(6)

  end subroutine init_div_damping_coeff

  !-------------------------------------------------------

  subroutine deallocate_alpha_div

     deallocate( alpha_div )

  end subroutine deallocate_alpha_div


    !.. subroutine acoust_new_rk(OneGrid, mzp,mxp,myp,  &
       !.. nnacoust_loc,                           &
       !.. dn0,pi0,th0,                            &
       !.. up,vp,wp,pp,                            &
       !.. ut,vt,wt,pt,                            &
       !.. topt,topu,topv,rtgt,rtgu,f13u,dxu,rtgv, &
       !.. dyv,f23v,f13t,f23t,fmapui,fmapvi,dxt,dyt,fmapt,lrk)
!.. !> @brief: Acoustic terms small time-step driver
!.. !! @details: This routine calls all the necessary routines to march the model
!.. !!     through the nnacoust_loc small timesteps.
!.. !!     The integration starts form the state up, vp, wp, pp; these fields
!.. !!     are overwritten with the updated fields at output.
!.. !! @author:  unknow
!.. !! @date:  18/Nov/2015
!.. !! @version:  5.2
!.. !! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv
!.. !!
!.. !! @copyright Under CC-GPL License by INPE/CPTEC
!.. !! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt


    !.. use ModParallelEnvironment, only: MsgDump

    !.. use ModNamelistFile, only: namelistfile

    !.. use ModGrid, only: Grid, DumpGrid

    !.. use ModMessageSet, only: &
         !.. PostRecvSendMsgs, &
         !.. WaitRecvMsgs

    !.. use mem_grid, only : nnacoust, & ! intent(in)
         !.. ngrid,                    & ! intent(in)
         !.. dts,                      & ! intent(out)
         !.. dtlt,                     & ! intent(in)
         !.. nxtnest,                  & ! intent(in)
         !.. nzp, nxp, nyp,            & ! intent(in)
         !.. impl,                     & ! intent(in)
         !.. time,                     &
	 !.. dyncore_flag,             &
         !.. grid_g                      ! for get_wind_div

    !.. use node_mod, only :  nmachs,  & ! intent(in)
         !.. ia,                       & ! intent(in)
         !.. iz,                       & ! intent(in)
         !.. ja,                       & ! intent(in)
         !.. jz,                       & ! intent(in)
         !.. izu,                      & ! intent(in)
         !.. jzv,                      & ! intent(in)
         !.. mynum

    !.. use ModNamelistFile, only: &
         !.. namelistfile

    !.. use mem_scratch, only: scratch

    !.. use ModComm, only: commHaloAcou

    !.. implicit none

    !.. type(Grid), pointer :: OneGrid
    !.. integer, intent(in) :: lrk

    !.. integer, intent(in) :: mzp
    !.. integer, intent(in) :: mxp
    !.. integer, intent(in) :: myp

    !.. integer, intent(in) :: nnacoust_loc  ! number of small time steps

    !.. real, intent(in) :: dn0(mzp,mxp,myp)
    !.. real, intent(in) :: pi0(mzp,mxp,myp)
    !.. real, intent(in) :: th0(mzp,mxp,myp)

    !.. real, intent(inout) :: up(mzp,mxp,myp)
    !.. real, intent(inout) :: vp(mzp,mxp,myp)
    !.. real, intent(inout) :: wp(mzp,mxp,myp)
    !.. real, intent(inout) :: pp(mzp,mxp,myp)

    !.. real, intent(in) :: topt(mxp,myp)
    !.. real, intent(in) :: topu(mxp,myp)
    !.. real, intent(in) :: topv(mxp,myp)
    !.. real, intent(in) :: rtgt(mxp,myp)
    !.. real, intent(in) :: rtgu(mxp,myp)
    !.. real, intent(in) :: f13u(mxp,myp)
    !.. real, intent(in) :: dxu(mxp,myp)
    !.. real, intent(in) :: rtgv(mxp,myp)
    !.. real, intent(in) :: dyv(mxp,myp)
    !.. real, intent(in) :: f23v(mxp,myp)
    !.. real, intent(in) :: f13t(mxp,myp)
    !.. real, intent(in) :: f23t(mxp,myp)
    !.. real, intent(in) :: fmapui(mxp,myp)
    !.. real, intent(in) :: fmapvi(mxp,myp)
    !.. real, intent(in) :: dxt(mxp,myp)
    !.. real, intent(in) :: dyt(mxp,myp)
    !.. real, intent(in) :: fmapt(mxp,myp)

    !.. real, intent(in) :: ut(mzp,mxp,myp)
    !.. real, intent(in) :: vt(mzp,mxp,myp)
    !.. real, intent(in) :: wt(mzp,mxp,myp)
    !.. real, intent(in) :: pt(mzp,mxp,myp)

    !.. include "tsNames.h"

    !.. integer :: iter,k
    !.. integer :: lastIter
    !.. logical :: singleProcRun
    !.. logical :: outermostGrid
    !.. real :: a1da2
    !.. real :: acoaa(mzp,mxp,myp)
    !.. real :: acoc(mzp,mxp,myp)
    !.. real :: acof(mzp,mxp,myp)
    !.. real :: acog(mzp,mxp,myp)
    !.. real :: amoe(mzp,mxp,myp)
    !.. real :: amof(mzp,mxp,myp)
    !.. real :: amog(mzp,mxp,myp)
    !.. real :: heatfx1(mxp,myp)

    !.. integer :: i, j
    !.. character(len=1) :: citer

    !.. real, dimension(:,:,:), allocatable :: &
       !.. div,        &
       !.. pp_minus_div&
      !.. ,pp_t_minus_dt
    !.. character(len=*), parameter :: h="**(acoust_new)**"
    !.. logical, parameter :: dumpLocal=.false.
    !.. character(LEN=5) :: ctime
    !.. integer :: nmbr_gpts

    !.. lastIter = nnacoust_loc
    !.. singleProcRun = nmachs == 1
    !.. outermostGrid = nxtnest(ngrid) == 0

    !.. if ( apply_div_damping .and. (dyncore_flag == 2 .or. dyncore_flag == 3) ) then
      !.. allocate( pp_minus_div ( mzp, mxp, myp) )
      !.. if(damp_formulation==1) allocate( div          ( mzp, mxp, myp) )
      !.. if(damp_formulation==2) allocate( pp_t_minus_dt( mzp, mxp, myp) )
    !.. end if

    !.. ! computes acoc, acof, acog, a1da2, amoe, amof, acoaa
    !.. ! uses dn0, pi0, th0, rtgt
    !.. call coefz(mzp,mxp,myp,ia,iz,ja,jz,  &
         !.. acoc,acof,acog,dn0,pi0,th0,rtgt,a1da2,amoe,amof,acoaa)

    !.. ! run small time steps
    !.. do iter=1,lastIter
      !.. write(citer,fmt='(I1)') iter

      !.. if ( apply_div_damping .and. ( dyncore_flag == 2 .or. dyncore_flag == 3 )&
        !.. .and. ( iter > 1 ) ) then
        !.. if(damp_formulation==1) then
          !.. ! application of divergence damping in the 1st timestep is detrimental to stability
          !.. call get_wind_div_v2 (mzp,mxp,myp,ia,iz,ja,jz,izu,jzv    &
            !.. ,up(1,1,1),vp(1,1,1),wp(1,1,1)                   &
            !.. ,scratch%vt3da(1), scratch%vt3db(1), scratch%vt3dh(1) &
            !.. ,div(1,1,1)          &
            !.. ,rtgt(1,1),rtgu(1,1),dxu(1,1),rtgv(1,1),dyv(1,1),f13t(1,1)&
            !.. ,f23t(1,1),fmapui(1,1),fmapvi(1,1)   &
            !.. ,grid_g(ngrid)%fmapu(1,1), grid_g(ngrid)%fmapv(1,1)&
            !.. ,grid_g(ngrid)%dxt,grid_g(ngrid)%dyt,fmapt)
            !.. if (.not. singleProcRun) then
              !.. call commHaloAcou(div,mxp,myp,mzp,myNum,'div')
            !.. endif
            !.. ! as proposed in Wicker, Skamarock (2002) (?) divergence damping is
            !.. ! used in an approximated form by adding the following term to the pressure:
            !.. pp_minus_div(:,:,:) = pp(:,:,:) - alpha_div(:,:,:) / th0(:,:,:) * div(:,:,:)
        !.. elseif(damp_formulation==2) then
          !.. !- alternative formulation
          !.. pp_minus_div(:,:,:) = pp(:,:,:) + div_damp_strength*(pp(:,:,:) - pp_t_minus_dt(:,:,:))
        !.. endif
      !.. end if

      !.. if(damp_formulation==2) pp_t_minus_dt(:,:,:)=pp(:,:,:)

!.. call dumpVarAllLatLonk3P(pp , 'PP'  ,1612,lrk,iter,1,mxp,1,myp,1,mzp,0.0,0.0,h)
!.. call dumpVarAllLatLonk3P(up , 'UP'  ,1612,lrk,iter,1,mxp,1,myp,1,mzp,0.0,0.0,h)
!.. call dumpVarAllLatLonk3P(ut , 'UT'  ,1612,lrk,iter,1,mxp,1,myp,1,mzp,0.0,0.0,h)
!.. call dumpVarAllLatLonk3P(pp_minus_div, 'PPMD'  ,1612,lrk,iter,1,mxp,1,myp,1,mzp,0.0,0.0,h)

      !.. ! up = f(up, pp); uses pp(i+1,j)
      !.. ! remaining arguments are loop invariant
      !.. if ( apply_div_damping .AND. ( dyncore_flag == 2 .or. dyncore_flag == 3 ) &
        !.. .AND. (iter > 1) )then
        !.. call prdctu(mzp,mxp,myp,ia,izu,ja,jz,  &
            !.. up,ut,pp_minus_div,th0,f13u,rtgu,rtgt,dxu,topu)
      !.. else
        !.. call prdctu(mzp,mxp,myp,ia,izu,ja,jz,  &
              !.. up,ut,pp,th0,f13u,rtgu,rtgt,dxu,topu)
      !.. end if

!.. call dumpVarAllLatLonk3P(up , 'UP'  ,1631,lrk,iter,1,mxp,1,myp,1,mzp,0.0,0.0,h)
!.. call dumpVarAllLatLonk3P(vp , 'VP'  ,1631,lrk,iter,1,mxp,1,myp,1,mzp,0.0,0.0,h)
!.. call dumpVarAllLatLonk3P(vt , 'VT'  ,1631,lrk,iter,1,mxp,1,myp,1,mzp,0.0,0.0,h)

      ! vp = f(vp, pp); uses pp(i,j+1)
      !.. ! remaining arguments are loop invariant
      !.. if ( apply_div_damping .AND. ( dyncore_flag == 2 .or. dyncore_flag == 3 ) &
          !.. .AND. (iter > 1) )then
        !.. call prdctv(mzp,mxp,myp,ia,iz,ja,jzv,&
               !.. vp,vt, pp_minus_div, th0,f23v,rtgv,rtgt,dyv,topv)
      !.. else
        !.. call prdctv(mzp,mxp,myp,ia,iz,ja,jzv,&
               !.. vp,vt, pp,th0,f23v,rtgv,rtgt,dyv,topv)
      !.. end if

!.. call dumpVarAllLatLonk3P(vp , 'VP'  ,1662,lrk,iter,1,mxp,1,myp,1,mzp,0.0,0.0,h)

! wp = f(wp, pp); uses node inner cells only
      !.. ! remaining arguments are loop invariant
      !.. call prdctw1(mzp,mxp,myp,ia,iz,ja,jz, &
            !.. wp,wt,pp,acoc,a1da2,rtgt,topt)

      !.. if (.not. singleProcRun) then
        !.. call commHaloAcou(up,mxp,myp,mzp,myNum,'up')
        !.. call commHaloAcou(vp,mxp,myp,mzp,myNum,'vp')
      !.. end if

!.. call dumpVarAllLatLonk3P(up , 'UP'  ,1685,lrk,iter,1,mxp,1,myp,1,mzp,0.0,0.0,h)
!.. call dumpVarAllLatLonk3P(vp , 'VP'  ,1685,lrk,iter,1,mxp,1,myp,1,mzp,0.0,0.0,h)


      !.. ! pp = f(pp, up, vp); uses up(i-1,j), vp(i,j-1)
      !.. ! also computes heatfx1
      !.. ! remaining arguments are loop invariant
      !.. call prdctp1(mzp,mxp,myp,ia,iz,ja,jz,  &
            !.. pp,up,vp,pi0,dn0,th0,pt,f13t,f23t,rtgt,rtgu,rtgv, &
            !.. heatfx1,fmapui,fmapvi,dxt,dyt,fmapt)

      !.. if (.not. singleProcRun) then
        !.. call commHaloAcou(wp,mxp,myp,mzp,myNum,'wp')
        !.. call commHaloAcou(pp,mxp,myp,mzp,myNum,'pp')
      !.. end if

!.. call dumpVarAllLatLonk3P(wp , 'WP'  ,1704,lrk,iter,1,mxp,1,myp,1,mzp,0.0,0.0,h)
!.. call dumpVarAllLatLonk3P(pp , 'PP'  ,1704,lrk,iter,1,mxp,1,myp,1,mzp,0.0,0.0,h)

! wp = f(wp, pp); uses node inner cells only
      !.. ! uses heatfx1
      !.. ! also computes amog
      !.. ! remaining arguments are loop invariant
      !.. call prdctw2(mzp,mxp,myp,ia,iz,ja,jz, &
            !.. wp,pp,acoc,amof,amog,acoaa,rtgt,heatfx1)

      !.. ! finishes updating wp=f(wp); uses node inner cells only
      !.. ! uses amog
      !.. ! remaining arguments are loop invariant
      !.. call prdctw3(mzp,mxp,myp,ia,iz,ja,jz,wp,amog,amoe,impl)

      !.. ! finishes updating pp=f(pp,wp); uses node inner cells only
      !.. ! remaining arguments are loop invariant
      !.. call prdctp2(mzp,mxp,myp,ia,iz,ja,jz,pp,wp,acof,acog)

      !.. if (.not. singleProcRun) then
        !.. call commHaloAcou(pp,mxp,myp,mzp,myNum,'pp')
        !.. call commHaloAcou(wp,mxp,myp,mzp,myNum,'wp')
      !.. end if

!.. call dumpVarAllLatLonk3P(pp , 'PP'  ,1732,lrk,iter,1,mxp,1,myp,1,mzp,0.0,0.0,h)
!.. call dumpVarAllLatLonk3P(wp , 'WP'  ,1732,lrk,iter,1,mxp,1,myp,1,mzp,0.0,0.0,h)
!.. call dumpVarAllLatLonk3P(up , 'UP'  ,1732,lrk,iter,1,mxp,1,myp,1,mzp,0.0,0.0,h)
!.. call dumpVarAllLatLonk3P(vp , 'VP'  ,1732,lrk,iter,1,mxp,1,myp,1,mzp,0.0,0.0,h)

    !.. enddo

    !.. if ( apply_div_damping .AND. ( dyncore_flag == 2 .or. dyncore_flag == 3) ) then
      !.. if (allocated(div          )) deallocate(div          )
      !.. if (allocated(pp_minus_div )) deallocate(pp_minus_div )
      !.. if (allocated(pp_t_minus_dt)) deallocate(pp_t_minus_dt)
    !.. end if

  !.. end subroutine acoust_new_rk



end module ModAcoust

!-------------------------------------------------------

subroutine get_wind_div_v2(m1,m2,m3,ia,iz,ja,jz,izu,jzv &
           ,up,vp,wp,flxu,flxv,flxw,div  &
	   ,rtgt,rtgu,dxu,rtgv          &
           ,dyv,f13t,f23t,fmapui,fmapvi,fmapu,fmapv,dxt,dyt,fmapt)
   use mem_grid, only :       & !intent(in)
       jdim,                  & !intent(in)
       hw4,                   & !intent(in)
       dzm, itopo,dzt         !intent(in)

 ! author: Saulo Freitas

 implicit none
 integer, intent(in):: m1,m2,m3,ia,iz,ja,jz,izu,jzv
 real, dimension(m1,m2,m3),intent(in):: up,vp,wp
 real, dimension(m2,m3)   ,intent(in):: rtgt,rtgu,dxu,rtgv  &
                           ,dyv,f13t,f23t,fmapui,fmapvi&
			   ,fmapu,fmapv,dxt,dyt,fmapt
 real, dimension(m1,m2,m3),intent(inout):: flxu,flxv,flxw
 real, dimension(m1,m2,m3),intent(out):: div

 integer i,j,k,im,jm
 real dummy,c1x,c1z,c1y

 div=0
   !print*, "subr. [get_wind_div_v2] ..."

   !- get div_X--------------------------
   do j = 1,m3
      do i = 1,m2
    	dummy = rtgu(i,j)  * fmapui(i,j)

    	do k = 1,m1
    	   flxu(k,i,j) = up(k,i,j)  * dummy
    	enddo
      enddo
    enddo


    do j = ja,jz
     do i = ia,izu
        c1z = 1.0 / rtgt(i,j)
        c1x = c1z * fmapt(i,j) * dxt(i,j)
        do k = 2,m1-1
          div(k,i,j) = c1x * ( flxu(k,i,j) - flxu(k,i-1,j) )
    	enddo
      enddo
    enddo

    !- get div_Y --------------------------
    do j = 1,m3
      do i = 1,m2
    	dummy = rtgv(i,j)  * fmapvi(i,j)

    	do k = 1,m1
    	   flxv(k,i,j) = vp(k,i,j)  * dummy
    	enddo
      enddo
    enddo
    do j = ja,jzv
     do i = ia,iz
        c1z = 1.0 / rtgt(i,j)
        c1y = c1z * fmapt(i,j) * dyt(i,j)
        do k = 2,m1-1
           div(k,i,j) = c1y * ( flxv(k,i,j) - flxv(k,i,j-1) ) + div(k,i,j)
    	enddo; enddo;enddo

    !- get div_Z --------------------------
  if(itopo.eq.0) then
    do j = 1,m3
      do i = 1,m2
         do k = 1,m1-1
            flxw(k,i,j) = wp(k,i,j)
         enddo
      enddo
    enddo
  else
    do j = 1,m3
      jm = max(j-1,1)
      do i = 1,m2
         im = max(i-1,1)
         do k = 1,m1-1
            flxw(k,i,j) = wp(k,i,j)  &
               + hw4(k) * ((up(k,i ,j) + up(k+1,i ,j)  &
               +            up(k,im,j) + up(k+1,im,j)) * f13t(i,j)  &
               +           (vp(k,i ,j) + vp(k+1,i ,j)  &
               +            vp(k,i,jm) + vp(k+1,i,jm)) * f23t(i,j)  &
	                  )
         enddo
      enddo
    enddo
  endif

  do j = ja,jz
   do i = ia,iz
      c1z = 1.0 / rtgt(i,j)

      do k = 2,m1-2
!srf check
        div(k,i,j) = c1z * dzt(k) * ( flxw(k,i,j) - flxw(k-1,i,j) ) + div(k,i,j)
!       div(k,i,j) =       dzt(k) * ( flxw(k,i,j) - flxw(k-1,i,j) ) + div(k,i,j)
      enddo
    enddo
  enddo

  if(jzv /= jz) div(:  , :,jz)=div(   :,: ,jzv)
  if(izu /= iz) div(:  ,iz,: )=div(   :,izu,: )
  div(:   , 1,: )=div(   :,ia,: )
  div(:   ,m2,: )=div(   :,iz,: )
  div(:   ,: ,1 )=div(   :,: ,ja)
  div(:   ,: ,m3)=div(   :,: ,jz)
  div(1   ,: ,: )=div(   2,: ,: )
  div(m1  ,: ,: )=div(m1-2,: ,: )
  div(m1-1,: ,: )=div(m1-2,: ,: )

end subroutine get_wind_div_v2


!******************************************************************************

subroutine buoyancy( wt )
!> @brief: buoyancy
!! @author:  unknow
!! @date:  18/Nov/2015
!! @version:  5.2
!! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv, wt_ptr
!!
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt
  use mem_basic,only: basic_g
  use mem_grid, only: ngrid
  use node_mod, only: mzp,mxp,myp,ia,iz,ja,jz,mynum
  use micphys , only: level

  implicit none
  real, intent(inout) :: wt( mzp, mxp, myp )

  call boyanc(mzp,mxp,myp,ia,iz,ja,jz,level   &
   ,wt                                       &
   ,basic_g(ngrid)%theta ,basic_g(ngrid)%rtp  &
   ,basic_g(ngrid)%rv    ,basic_g(ngrid)%th0  &
   ,mynum                      )


end subroutine buoyancy

!******************************************************************************


subroutine boyanc(m1,m2,m3,ia,iz,ja,jz,level,wt,theta,rtp,rv,th0,mynum)
!> @brief: boyancy
!! Some remarks about the meaning of the timelevels of the fields:
!! in the leapfrog dynamical core, theta and rv are defined at timelevel t.
!! In RK they have the same "timelevel-meaning" as uc, vc, wc, pc.
!! @author:  unknow
!! @date:  18/Nov/2015
!! @version:  5.2
!! @param: mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv, wt
!!
!! @copyright Under CC-GPL License by INPE/CPTEC
!! Please, read @link https://creativecommons.org/licenses/GPL/2.0/legalcode.pt

use rconstants, only : gg

implicit none

integer, intent(in) :: m1,m2,m3,ia,iz,ja,jz,level,mynum
real, dimension(m1,m2,m3),intent(in   ) :: theta,rtp,rv,th0
real, dimension(m1,m2,m3),intent(inout) :: wt

integer, dimension(m2,m3) :: lpw
real, dimension(m1,m2,m3) :: vtemp
integer :: i,j,k

lpw(:,:)= 2

if (level .ge. 1) then
   do j = ja,jz
      do i = ia,iz
         do k = lpw(i,j),m1-1
            vtemp(k,i,j) = gg * ((theta(k,i,j) * (1. + .61 * rv(k,i,j))  &
               - th0(k,i,j)) / th0(k,i,j) - (rtp(k,i,j) - rv(k,i,j)) )
         enddo
      enddo
   enddo
else
   do j = ja,jz
      do i = ia,iz
         do k = lpw(i,j),m1-1
            vtemp(k,i,j) = gg * (theta(k,i,j) / th0(k,i,j) - 1.)
         enddo
      enddo
   enddo
endif

do j = ja,jz
   do i = ia,iz
      do k = lpw(i,j),m1-2
         wt(k,i,j) = wt(k,i,j) + vtemp(k,i,j) + vtemp(k+1,i,j)
      enddo
   enddo
enddo

end subroutine boyanc
