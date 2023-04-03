!###########################################################################
!  Advection scheme for the Runge-Kutta dynamical core
!    (after Wicker, Skamarock, 2002, MWR)
!  Dec/2015 by Saulo Freitas (INPE), Michael Baldauf (DWD)
!  Dec/2017 paralelized and several bug fixes by Luiz F. Rodrigues (INPE)
!###########################################################################
!
!-this routine is called by timestep_rk
!call advectc_rk('V'      ,mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum,l_rk)
!call advectc_rk('THETAIL',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum,l_rk)
!call advectc_rk('PI'     ,mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum,l_rk)
!call advectc_rk('SCALAR' ,mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum,l_rk)
!
!###########################################################################
  MODULE advRkParam
    !- parameters for various advection
    real, parameter :: f30 = 7./12.
    real, parameter :: f31 = 1./12.
    real, parameter :: f40 = 7./12.
    real, parameter :: f41 = 1./12.
    real, parameter :: f50 = 37./60.
    real, parameter :: f51 = 2./15.
    real, parameter :: f52 = 1./60.
    real, parameter :: f60 = 37./60.
    real, parameter :: f61 = 2./15.
    real, parameter :: f62 = 1./60.
    real, parameter :: eps = 1.e-20
    real, parameter :: real_init=-0.0*huge(1.)
    real :: fifth_order
  end module advRkParam
!
  subroutine advectc_rk(varn,mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum,l_rk)
    !use var_tables, only: scalar_table
    use grid_dims, only: maxgrds
    use mem_tend, only: tend
    use var_tables, only: num_scalar, scalar_tab
    use mem_grid, only: ngrid, nzpmax, grid_g, dtlt, if_adap, jdim, time, &
       zt, zm, dzm, dzt, hw4,itopo,pd_or_mnt_constraint,order_h,order_v

    use mem_basic, only: basic_g
    use mem_chem1, only: nspecies_transported
    use mem_stilt, only: stilt_g,iexev

    implicit none
    include "i8.h"
    integer, intent(in) :: mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum,l_rk
    character(len=*), intent(in) :: varn
    !
    ! Local Variables
    integer(kind=i8) :: mxyzp
    integer :: i,j,k,n
    real, pointer :: scalarp, scalart
    integer :: i_scl,is,js,ks
    !- scratchs (local arrays)
    real :: vt3da(mzp,mxp,myp)
    real :: vt3db(mzp,mxp,myp)
    real :: vt3dc(mzp,mxp,myp)
    real :: vt3dh(mzp,mxp,myp)
    real :: vt3dj(mzp,mxp,myp)
    real :: vt3dk(mzp,mxp,myp)
    real :: vctr1(mzp)
    real :: vctr2(mzp)

    real :: mfx_wind(mzp,mxp,myp)
    real :: mfy_wind(mzp,mxp,myp)
    real :: mfz_wind(mzp,mxp,myp)

    mxyzp = mxp * myp * mzp

    vt3da=0.0
    vt3db=0.0
    vt3dc=0.0
    vt3dh=0.0
    vt3dj=0.0
    vt3dk=0.0
    vctr1=0.0
    vctr2=0.0

    if (trim(varn) .eq. 'V') then

      mfx_wind=0.0 ;mfy_wind=0.0 ;mfz_wind =0.0

      ! Advect  U, V, and W
      ! input: mzp,mxp,myp,ia,iz,ja,jz,izu,jzv
      !        basic_g%uc,%vc,%wc,%dn0,%dn0u,%dn0v
      !        grid_g%dxt,%dxu,%dxv,%dyt,%dyu,%dyv,%rtgt,%rtgu,%rtgv
      !        grid_g%f13t,%f23t,%fmapt,%fmapu,%fmapv,%fmapui,%fmapvi
      ! output: tend%ut,%vt,%wt
      !
      !--------------- U-advect
      is=1
      js=0
      ks=0

      call mf_wind(mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,itopo,hw4,jdim,dzt,dzm  &
          ,basic_g(ngrid)%uc,basic_g(ngrid)%vc,basic_g(ngrid)%wc              &
          ,basic_g(ngrid)%dn0,basic_g(ngrid)%dn0u,basic_g(ngrid)%dn0v&
          ,grid_g(ngrid)%dxt,grid_g(ngrid)%dxu,grid_g(ngrid)%dxv     &
          ,grid_g(ngrid)%dyt,grid_g(ngrid)%dyu,grid_g(ngrid)%dyv     &
          !
          ,grid_g(ngrid)%rtgt,grid_g(ngrid)%rtgu,grid_g(ngrid)%rtgv  &
          ,grid_g(ngrid)%f13t,grid_g(ngrid)%f23t,grid_g(ngrid)%fmapt &
          ,grid_g(ngrid)%fmapu ,grid_g(ngrid)%fmapv                  &
          ,grid_g(ngrid)%fmapui,grid_g(ngrid)%fmapvi                 &
          ,vt3da,vt3db,vt3dc,mfx_wind,mfy_wind,mfz_wind,is,js,ks)

      call advect_ws(mzp,mxp,myp,ia,iz,ja,jz &
                    ,basic_g(ngrid)%uc &! field being advected
                    ,vt3da    & ! uc*dn0u*fmapui*rtgu = rhou*U
                    ,vt3db    & ! similar for v
                    ,vt3dc    & ! similar for sigma_dot
                    ,mfx_wind & ! fmapt*rtgti*dxt/dn0 = 1(rho dx)
                    ,mfy_wind & ! similar for v
                    ,mfz_wind & ! similar for sigma_dot
                    !
                    ,tend%ut_rk   &
                    ,is,js,ks     &
                    ,pd_or_mnt_constraint &
                    ,order_h,order_v      &
                    ,dtlt,                 &
                    'uc' &
                    )

      !--------------- V-advect
      is=0
      js=1
      ks=0

      call mf_wind(mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,itopo,hw4,jdim,dzt,dzm  &
           ,basic_g(ngrid)%uc,basic_g(ngrid)%vc,basic_g(ngrid)%wc              &
           ,basic_g(ngrid)%dn0,basic_g(ngrid)%dn0u,basic_g(ngrid)%dn0v&
           ,grid_g(ngrid)%dxt,grid_g(ngrid)%dxu,grid_g(ngrid)%dxv     &
           ,grid_g(ngrid)%dyt,grid_g(ngrid)%dyu,grid_g(ngrid)%dyv     &
           !
           ,grid_g(ngrid)%rtgt,grid_g(ngrid)%rtgu,grid_g(ngrid)%rtgv  &
           ,grid_g(ngrid)%f13t,grid_g(ngrid)%f23t,grid_g(ngrid)%fmapt &
           ,grid_g(ngrid)%fmapu ,grid_g(ngrid)%fmapv                  &
           ,grid_g(ngrid)%fmapui,grid_g(ngrid)%fmapvi                 &
           ,vt3da,vt3db,vt3dc,mfx_wind,mfy_wind,mfz_wind,is,js,ks)

      call advect_ws(mzp,mxp,myp,ia,iz,ja,jz,basic_g(ngrid)%vc &
                     ,vt3da    & ! uc*dn0u*fmapui*rtgu = rhou*V
                     ,vt3db    & ! similar for v
                     ,vt3dc    & ! similar for sigma_dot
                     ,mfx_wind & ! fmapt*rtgti*dxt/dn0 = 1(rho dx)
                     ,mfy_wind & ! similar for v
                     ,mfz_wind & ! similar for sigma_dot
                     !
                     ,tend%vt_rk   &
                     ,is,js,ks     &
                     ,pd_or_mnt_constraint &
                     ,order_h,order_v      &
                     ,dtlt,                 &
                     'vc' &
                     )


      !--------------- W-advect
      is=0
      js=0
      ks=1

      call mf_wind(mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,itopo,hw4,jdim,dzt,dzm  &
           ,basic_g(ngrid)%uc,basic_g(ngrid)%vc,basic_g(ngrid)%wc              &
           ,basic_g(ngrid)%dn0,basic_g(ngrid)%dn0u,basic_g(ngrid)%dn0v&
           ,grid_g(ngrid)%dxt,grid_g(ngrid)%dxu,grid_g(ngrid)%dxv     &
           ,grid_g(ngrid)%dyt,grid_g(ngrid)%dyu,grid_g(ngrid)%dyv     &
           !
           ,grid_g(ngrid)%rtgt,grid_g(ngrid)%rtgu,grid_g(ngrid)%rtgv  &
           ,grid_g(ngrid)%f13t,grid_g(ngrid)%f23t,grid_g(ngrid)%fmapt &
           ,grid_g(ngrid)%fmapu ,grid_g(ngrid)%fmapv                  &
           ,grid_g(ngrid)%fmapui,grid_g(ngrid)%fmapvi                 &
           ,vt3da,vt3db,vt3dc,mfx_wind,mfy_wind,mfz_wind,is,js,ks)

      call advect_ws(mzp,mxp,myp,ia,iz,ja,jz,basic_g(ngrid)%wc &
                    ,vt3da    & ! uc*dn0u*fmapui*rtgu = rhou*W
                    ,vt3db    & ! similar for v
                    ,vt3dc    & ! similar for sigma_dot
                    ,mfx_wind & ! fmapt*rtgti*dxt/dn0 = 1(rho dx)
                    ,mfy_wind & ! similar for v
                    ,mfz_wind & ! similar for sigma_dot
                    !
                    ,tend%wt_rk   &
                    ,is,js,ks     &
                    ,pd_or_mnt_constraint&
                    ,order_h,order_v &
                    ,dtlt,                 &
                    'wc' &
                    )
    endif  ! endif of varn .eq. 'V'

    if (trim(varn) .eq. 'T'   .or. trim(varn) .eq. "PI"      .or. &
        trim(varn) .eq. "THA" .or. trim(varn) .eq. "THETAIL" .or. &
        trim(varn) .eq. 'SCALAR'  ) THEN

      !-------------- scalars advect
      is=0
      js=0
      ks=0

      ! input: basic_g%up,%uc,%vp,%vc,%wp,%wc,%dn0,%dn0u,%dn0v
      !        grid_g%rtgt,%rtgu,%rtgv,%fmapt,%fmapui,%fmapvi,%f13t,%f23t,%dxu,%dyv,%dxt,%dyt
      !        scalar_tab%var_p, %var_t
      ! output: scalar_tab%var_t

      if (trim(varn) .eq. 'T' .or. trim(varn) .eq. 'SCALAR' ) then
        !- combine the 2-time levels wind fields for tracers
        do j = 1,myp
          do i = 1,mxp
             do k = 1,mzp
                vt3da(k,i,j) = (basic_g(ngrid)%up(k,i,j)  &
                              + basic_g(ngrid)%uc(k,i,j)) * 0.5
                vt3db(k,i,j) = (basic_g(ngrid)%vp(k,i,j)  &
                              + basic_g(ngrid)%vc(k,i,j)) * 0.5
                vt3dc(k,i,j) = (basic_g(ngrid)%wp(k,i,j)  &
                              + basic_g(ngrid)%wc(k,i,j)) * 0.5
             end do
          end do
        end do
      else
        do j = 1,myp
          do i = 1,mxp
             do k = 1,mzp
                vt3da(k,i,j) = basic_g(ngrid)%uc(k,i,j)
                vt3db(k,i,j) = basic_g(ngrid)%vc(k,i,j)
                vt3dc(k,i,j) = basic_g(ngrid)%wc(k,i,j)
             end do
          end do
        end do
      endif

      ! input: vt3da,vt3db,vt3dc
      !        basic_g%dn0,%dn0u,%dn0v
      !        grid_g%rtgt,%rtgu,%rtgv,%fmapt,%fmapui,%fmapvi,%f13t,%f23t,%dxu,%dyv,%dxt,%dyt
      ! output:vt3da,vt3db,vt3dc,vt3dh,vt3dj,vt3dk
      !
      call fa_preptc_rk(mzp,mxp,myp    &
                   ,vt3da,vt3db,vt3dc,vt3dh,vt3dj,vt3dk,vctr1,vctr2              &
                   ,basic_g(ngrid)%dn0,basic_g(ngrid)%dn0u,basic_g(ngrid)%dn0v   &
                   ,grid_g(ngrid)%rtgt,grid_g(ngrid)%rtgu,grid_g(ngrid)%rtgv     &
                   ,grid_g(ngrid)%fmapt,grid_g(ngrid)%fmapui,grid_g(ngrid)%fmapvi&
                   ,grid_g(ngrid)%f13t,grid_g(ngrid)%f23t                        &
                   ,grid_g(ngrid)%dxu,grid_g(ngrid)%dyv                          &
                   ,grid_g(ngrid)%dxt,grid_g(ngrid)%dyt,hw4,dzm,dzt,zm,zt)

!      if(trim(varn)=='THETAIL') call dumpVarAllLatLonk(basic_g(ngrid)%thc,'THC',274,l_rk,0,1,mxp,1,myp,1,mzp,0.0,0.0) !ok

      if ( trim(varn) .eq. "THETAIL" ) THEN

        !.. call dumpVarAllLatLonk(basic_g(ngrid)%thc,'THC',278,l_rk,0,1,mxp,1,myp,1,mzp,0.0,0.0) !
        !.. call dumpVarAllLatLonk(tend%tht_rk,'THT',279,l_rk,0,1,mxp,1,myp,1,mzp,0.0,0.0)

        call advect_ws(mzp,mxp,myp,ia,iz,ja,jz,basic_g(ngrid)%thc &
                      ,vt3da & ! uc*dn0u*fmapui*rtgu = rhou*U
                      ,vt3db & ! similar for v
                      ,vt3dc & ! similar for sigma_dot
                      ,vt3dh & ! fmapt*rtgti*dxt/dn0 = 1(rho dx)
                      ,vt3dj & ! similar for v
                      ,vt3dk & ! similar for sigma_dot
                      !
                      ,tend%tht_rk         &
                      ,is,js,ks            &
                      ,pd_or_mnt_constraint&
                      ,order_h,order_v     &
                      ,dtlt,                &
                      'thc'  &
                      )

!        call dumpVarAllLatLonk(tend%tht_rk,'THT',297,l_rk,0,1,mxp,1,myp,1,mzp,0.0,0.0)

        return
      endif !endif of varn .eq. 'THETAIL'

      if ( trim(varn) .eq. "THA" .and. iexev == 2) THEN
           !-get log(thetav)
        call prep_lnthetv(mzp,mxp,myp,ia,iz,ja,jz&
                      ,basic_g(ngrid)%theta &
                      ,basic_g(ngrid)%rtp   &
                      ,basic_g(ngrid)%rv    &
                      ,stilt_g(ngrid)%lnthetav)

        call advect_ws(mzp,mxp,myp,ia,iz,ja,jz &
                      ,stilt_g(ngrid)%lnthetav   &! advected field
                      ,vt3da & ! uc*dn0u*fmapui*rtgu = rhou*U
                      ,vt3db & ! similar for v
                      ,vt3dc & ! similar for sigma_dot
                      ,vt3dh & ! fmapt*rtgti*dxt/dn0 = 1(rho dx)
                      ,vt3dj & ! similar for v
                      ,vt3dk & ! similar for sigma_dot
                      !
                      ,stilt_g(ngrid)%lnthvadv & !tendency field including advection
                      ,is,js,ks,pd_or_mnt_constraint&
                      ,order_h,order_v               &
                      ,dtlt,                &
                      'lnthetav' &
                      )

        return
      endif !endif og varn .eq. 'THA'
      !
      if ( trim(varn) .eq. "PI" .and. iexev == 2) THEN
        call advect_ws(mzp,mxp,myp,ia,iz,ja,jz&
                      ,basic_g(ngrid)%pc & !advected field
                      ,vt3da & ! uc*dn0u*fmapui*rtgu = rhou*U
                      ,vt3db & ! similar for v
                      ,vt3dc & ! similar for sigma_dot
                      ,vt3dh & ! fmapt*rtgti*dxt/dn0 = 1(rho dx)
                      ,vt3dj & ! similar for v
                      ,vt3dk & ! similar for sigma_dot
                      !
                      ,tend%pt_rk          & !tendency from advection
                      ,is,js,ks            &
                      ,pd_or_mnt_constraint&
                      ,order_h,order_v               &
                      ,dtlt,                &
                      'pc' &
                      )
        return
      endif !endif og varn .eq. 'PI'

      if (trim(varn) .eq. 'T' .or. trim(varn) .eq. 'SCALAR') THEN

        i_scl=num_scalar(ngrid)  !- all scalars
        !i_scl= num_scalar(ngrid) - NSPECIES_TRANSPORTED !- only theta_il+water+tke

        do n=1,i_scl
          !
	  !- if RK or ABM3 schemes, THP/THC are not transported here
          if (scalar_tab(n,ngrid)%name == 'THC' .or. &
              scalar_tab(n,ngrid)%name == 'THP') cycle

         !if(mynum==1)print*,"scalars=",scalar_tab(n,ngrid)%name,order_h,order_v;call flush(6)
          scalarp => scalar_tab(n,ngrid)%var_p
          scalart => scalar_tab(n,ngrid)%var_t

          ! input: scalarp, scalart, dtlt
          ! output: scalart

          call advect_ws(mzp,mxp,myp,ia,iz,ja,jz &
	              ,scalarp & !scalar being advected 
                      ,vt3da   & ! 0.5(up+uc)*dn0u*fmapui*rtgu = rhou*U
                      ,vt3db   & ! similar for v
                      ,vt3dc   & ! similar for sigma_dot
                      ,vt3dh   & ! fmapt*rtgti*dxt/dn0 = 1(rho dx)
                      ,vt3dj   & ! similar for v
                      ,vt3dk   & ! similar for sigma_dot
                      !
                      ,scalart,is,js,ks    & !scalar tendency
                      ,pd_or_mnt_constraint& ! 
                      ,order_h,order_v     & !order horiz/vert 
                      ,dtlt                & !timestep
                      ,scalar_tab(n,ngrid)%name & ! scalar name
                      )

        end do

      endif  !endif og varn .eq. 'T'

    endif !endif of varn .eq. 'T' .or. varn .eq. "PI"

  end subroutine advectc_rk

!---------------------------------------------------------------------

  subroutine fa_preptc_rk(m1,m2,m3,vt3da,vt3db,vt3dc, &
       vt3dh,vt3dj,vt3dk,vctr1, vctr2,dn0,dn0u,dn0v,  &
       rtgt,rtgu,rtgv,fmapt,fmapui,fmapvi,f13t,f23t,  &
       dxu,dyv,dxt,dyt,hw4, dzm, dzt, zm, zt)


    implicit none
    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    real, intent(in) :: hw4(m1), dzm(m1), dzt(m1), zm(m1), zt(m1)
    real, intent(inout) :: vt3da(m1,m2,m3)
    real, intent(inout) :: vt3db(m1,m2,m3)
    real, intent(inout) :: vt3dc(m1,m2,m3)
    real, intent(out) :: vt3dh(m1,m2,m3)
    real, intent(out) :: vt3dj(m1,m2,m3)
    real, intent(out) :: vt3dk(m1,m2,m3)
    real, intent(in) :: dn0(m1,m2,m3)
    real, intent(in) :: dn0u(m1,m2,m3)
    real, intent(in) :: dn0v(m1,m2,m3)
    real, intent(in) :: rtgt(m2,m3)
    real, intent(in) :: rtgu(m2,m3)
    real, intent(in) :: rtgv(m2,m3)
    real, intent(in) :: fmapt(m2,m3)
    real, intent(in) :: fmapui(m2,m3)
    real, intent(in) :: fmapvi(m2,m3)
    real, intent(in) :: f13t(m2,m3)
    real, intent(in) :: f23t(m2,m3)
    real, intent(in) :: dxu(m2,m3)
    real, intent(in) :: dyv(m2,m3)
    real, intent(in) :: dxt(m2,m3)
    real, intent(in) :: dyt(m2,m3)
    real, intent(out) :: vctr1(m1)
    real, intent(out) :: vctr2(m1)

    ! Local Variables
    integer :: j,i,k,im,ip,jm,jp
    real :: c1,c2,c3,c4,rtgti
    real :: vctr3(m1)


    ! VT3DA, VT3DB, and VT3DC are input as the velocity components (averaged
    ! between past and current time levels) times dtlt.

    ! Add contribution to VT3DC from horiz winds crossing sloping sigma surfaces,
    !    and include 1/rtgt factor in VT3DC

    do j = 1,m3
       jm = max(1,j-1)
       jp = min(m3,j+1)
       do i = 1,m2
          im = max(1,i-1)
          ip = min(m2,i+1)
          rtgti = 1. / rtgt(i,j)
          c1 = .5 * dxu(i,j)
          c2 = .5 * dyv(i,j)
          c3 = dxt(i,j) * fmapt(i,j) * rtgti
          c4 = dyt(i,j) * fmapt(i,j) * rtgti

          do k = 1,m1-1
             vt3dc(k,i,j) = ((vt3da(k,i,j) + vt3da(k+1,i,j)  &
                  + vt3da(k,im,j) + vt3da(k+1,im,j)) * f13t(i,j)  &
                  + (vt3db(k,i,j) + vt3db(k+1,i,j) + vt3db(k,i,jm)  &
                  + vt3db(k+1,i,jm)) * f23t(i,j)) * hw4(k)  &
                  + vt3dc(k,i,j) * rtgti
             !vt3dd(k,i,j) = c1 * vt3da(k,i,j)
             !vt3de(k,i,j) = c2 * vt3db(k,i,j)
             !vt3df(k,i,j) = .5 * vt3dc(k,i,j) * dzm(k)
             vctr3(k) = 1. / dn0(k,i,j)
             !-- these are the metric factors for u,v and z directions
             !-- divided by ( air_dens times grid spacing).
             vt3dh(k,i,j) = c3 * vctr3(k)
             vt3dj(k,i,j) = c4 * vctr3(k)
             vt3dk(k,i,j) = dzt(k) * vctr3(k)
          end do

          !            vt3di(1,i,j) = dxu(i,j) / (dxu(i,j) + dxt(ip,j))
          !            vt3di(2,i,j) = dxu(i,j) / (dxu(i,j) + dxt(i,j))
          !            vt3di(3,i,j) = dyv(i,j) / (dyv(i,j) + dyt(i,jp))
          !            vt3di(4,i,j) = dyv(i,j) / (dyv(i,j) + dyt(i,j))
          ! temporary override
          !vt3di(1,i,j) = .5
          !vt3di(2,i,j) = .5
          !vt3di(3,i,j) = .5
          !vt3di(4,i,j) = .5
       end do
    end do

    do k = 1,m1-1
       vctr1(k) = (zt(k+1) - zm(k)) * dzm(k)
       vctr2(k) =  (zm(k) - zt(k)) * dzm(k)
    end do

    ! Convert velocity components * dtlt (VT3DA, VT3DB, VT3DC)
    ! into mass fluxes times dtlt.

    do j = 1,m3
       do i = 1,m2
          c1 = fmapui(i,j) * rtgu(i,j)
          c2 = fmapvi(i,j) * rtgv(i,j)
          do k = 1,m1-1
             vt3da(k,i,j) = vt3da(k,i,j) * c1 * dn0u(k,i,j)
             vt3db(k,i,j) = vt3db(k,i,j) * c2 * dn0v(k,i,j)
             vt3dc(k,i,j) = vt3dc(k,i,j) * .5  &
                  * (dn0(k,i,j) + dn0(k+1,i,j))
          end do
       end do
    end do
  end subroutine fa_preptc_rk

  !---------------------------------------------------------------------
  subroutine advect_ws(mzp,mxp,myp,ia,iz,ja,jz,scp,ufx,vfx,wfx &
                      ,vt3dh,vt3dj,vt3dk,sct,is,js,ks          &
                      ,pd_or_mnt_constraint,order_h,order_v,dt,vname)
    use advRkParam, only: fifth_order, eps,real_init
    use ModComm, only: copyMyPart, commHalo, border, &
        north, south, west, east,expandBorder
    use node_mod, only:  nmachs, myNum,nodei0,nodej0
    use mem_grid, only: time

    implicit none
    integer, intent(in) :: mzp !- z
    integer, intent(in) :: mxp !- x
    integer, intent(in) :: myp !- y
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: is,js,ks,pd_or_mnt_constraint,order_h,order_v
    real, intent(in)    :: dt
    real, intent(in)    :: scp(mzp,mxp,myp)
    real, intent(inout) :: sct(mzp,mxp,myp)
    real, dimension(mzp,mxp,myp), intent(in) :: ufx, vfx,wfx,vt3dh,vt3dj,vt3dk
    character(len=*),intent(in) :: vname

    logical :: scalar
    real,allocatable :: qz(:,:,:)
    real,allocatable :: qx(:,:,:)
    real,allocatable :: qy(:,:,:)
    real,allocatable :: scr      (:,:,:)
    real,allocatable :: ufx_local(:,:,:)
    real,allocatable :: vfx_local(:,:,:)
    real,allocatable :: wfx_local(:,:,:)

!---tmp 
!real   :: scp_new(mzp,mxp,myp)
!---tmp 
    !external functions
    real, external :: flux_upwind,fq2, fq3, fq4, fq5, fq6, fq

    logical, parameter :: IsToDump=.false.
    !< Do a dump os communications?
    logical, parameter :: variable=.true.
    !<
    integer, parameter :: mzi=-2, myi=-2, mxi=-2

    integer :: mzpp3,mxpp3,mypp3
    integer :: mzppks,mxppis,myppjs

    mzpp3=mzp+3; mxpp3=mxp+3; mypp3=myp+3
    mzppks=mzp+ks; mxppis=mxp+is; myppjs=myp+js

    allocate(qx(mzppks,mxppis,myppjs));qx=real_init
    allocate(qy(mzppks,mxppis,myppjs));qy=real_init
    allocate(qz(mzppks,mxppis,myppjs));qz=real_init
    allocate(scr      (mzi:mzpp3,mxi:mxpp3,myi:mypp3));scr=real_init
    allocate(ufx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3));ufx_local=real_init
    allocate(vfx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3));vfx_local=real_init
    allocate(wfx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3));wfx_local=real_init

    !- flag to determine if a scalar is being advected
    scalar = .false.
    IF(is==0 .and. js==0 .and. ks==0) scalar = .true.

    if(IsToDump) &
       call dumpXYvar(scp,vname,'a',1,mxp,1,myp,1,mzp,600.0,660.0)

    !- copy input arrays to extended arrays
    call copyMyPart(scp,scr,ufx_local,vfx_local,wfx_local, &
          ufx,vfx,wfx,mxp,myp,mzp,is,js,ks, &
          mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
          vname)

    if(.not. variable) &
      call expandBorder(mxp,myp,mzp,is,js,ks, &
             mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
             scr,ufx_local,vfx_local,wfx_local)

    if(IsToDump) &
       call dumpXYvar(wfx_local,vname,'b',-2,mxp+3,-2,myp+3,-2,mzp+3,0.0,0.0)

    ! Set x & y boundary values in halo zones
    if(nmachs>1) call commHalo(scr,ufx_local,vfx_local,wfx_local, &
                            mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
                            mxp,myp,mzp,myNum,nmachs,nodei0,nodej0,vname)

    if(IsToDump) &
       call dumpXYvar(wfx_local,vname,'c',-2,mxp+3,-2,myp+3,-2,mzp+3,0.0,0.0)

    select case (order_h)
      case (1)
        call compXYInterface_or1(mxp,myp,mzp,ks,is,js,&
                                   mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
                                   mzppks,mxppis,myppjs, &
                                   scr,ufx_local,vfx_local,&
                                   border(mynum,:),qx,qy, &
                                   variable, vname)
      case (2)
        call compXYInterface_or2(mxp,myp,mzp,ks,is,js,&
                                   mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
                                   mzppks,mxppis,myppjs, &
                                   scr,ufx_local,vfx_local,&
                                   border(mynum,:),qx,qy, &
                                   variable, vname)
      case (3)
        call compXYInterface_or3(mxp,myp,mzp,ks,is,js,&
                                   mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
                                   mzppks,mxppis,myppjs, &
                                   scr,ufx_local,vfx_local,&
                                   border(mynum,:),qx,qy, &
                                   variable, vname)
      case (4)
        call compXYInterface_or4(mxp,myp,mzp,ks,is,js,&
                                   mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
                                   mzppks,mxppis,myppjs, &
                                   scr,ufx_local,vfx_local,&
                                   border(mynum,:),qx,qy, &
                                   variable, vname)
      case(5,6)
        call compXYInterface_or56(mxp,myp,mzp,ks,is,js,&
                                   mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
                                   mzppks,mxppis,myppjs, &
                                   scr,ufx_local,vfx_local,&
                                   border(mynum,:),qx,qy, &
                                   variable, vname, order_h)
      case default
        write (*,fmt='(A)') 'Advect Error : the order_h must be from 1 to 6'
        stop 'ERROR!'
    end select

    select case (order_v)
      case (1)
        call compZInterface_or1(mxp,myp,mzp,ks,is,js,&
                                  mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
                                  mzppks,mxppis,myppjs, &
                                  scr,wfx_local,qz, &
                                  variable, vname)
      case (2)
        call compZInterface_or2(mxp,myp,mzp,ks,is,js,&
                                  mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
                                  mzppks,mxppis,myppjs, &
                                  scr,wfx_local,qz, &
                                  variable, vname)
      case (3)
        call compZInterface_or3(mxp,myp,mzp,ks,is,js,&
                                  mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
                                  mzppks,mxppis,myppjs, &
                                  scr,wfx_local,qz, &
                                  variable, vname)
      case (4)
        call compZInterface_or4(mxp,myp,mzp,ks,is,js,&
                                  mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
                                  mzppks,mxppis,myppjs, &
                                  scr,wfx_local,qz, &
                                  variable, vname)
      case(5,6)
        call compZInterface_or56(mxp,myp,mzp,ks,is,js,&
                                  mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
                                  mzppks,mxppis,myppjs, &
                                  scr,wfx_local,qz, &
                                  variable, vname,order_v)
      case default
        write (*,fmt='(A)') 'Advect Error : the order_v must be from 1 to 6'
        stop 'ERROR!'
    end select

    !print *, 'Is scalar? ',scalar,is,js,ks;call flush(6)
    !iErrNumber=dumpMessage(c_tty,c_yes,header,version,c_notice,vname//' is scalar?: ',(/logical2Int(scalar),is,js,ks/),'I1')

    IF(pd_or_mnt_constraint > 0 .and. scalar .and. vname .ne. 'thc' .and. vname .ne. 'pc') THEN  
   
    !-- positivity/monotonicity constraints
        call PosMonConstraints(mxp,myp,mzp,is,js,ks,ia,iz,ja,jz,&
                               pd_or_mnt_constraint, &
                               dt,ufx,vfx,wfx,vt3dh,vt3dj,vt3dk,scp, &
                               scr,ufx_local,vfx_local,wfx_local, &
                               qx,qy,qz,mzi,mzpp3,mxi,mxpp3,myi, &
                               mypp3,mzppks,mxppis,myppjs,mynum,vname,sct)

    ENDIF

    call CreateTendency(mxp,myp,mzp,is,js,ks,ia,iz,ja,jz, &
                                mzppks,mxppis,myppjs, &
                               dt,ufx,vfx,wfx, &
                               vt3dh,vt3dj,vt3dk,scp, &
                               qx,qy,qz,sct,vname,mynum)
!-----check neg mass
!IF(pd_or_mnt_constraint > 0 .and. scalar .and. vname .ne. 'thc' .and. vname .ne. 'pc') THEN  
!IF( scalar .and. vname .ne. 'thc' .and. vname .ne. 'pc') THEN  
!scp_new=scp+sct*dt
!where(scp_new<0.) sct=-scp/dt
!print*,"maxmin=",trim(vname),1.-abs(minval(scp_new)-maxval(scp_new))/(1.e-20+maxval(scp_new))
!call flush(6)
!ENDIF
!-----check neg mass


    deallocate(qx ,qy,    qz,    scr,          ufx_local,    vfx_local ,   wfx_local)

  end subroutine advect_ws
  !---------------------------------------------------------------------

  subroutine mf_wind(m1,m2,m3,ia,iz,ja,jz,izu,jzv,itopo, hw4, jdim, dzt, dzm,&
                     uc,vc,wc,dn0,dn0u,dn0v,                                 &
		     dxt,dxu,dxv,dyt,dyu,dyv,rtgt,rtgu,rtgv,f13t,f23t,       &
                     fmapt,fmapu,fmapv,fmapui,fmapvi,                        &
                     !
                     flxu,flxv,flxw, mfx_wind,mfy_wind,mfz_wind,is,js,ks)

    implicit none
    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: izu
    integer, intent(in) :: jzv
    integer, intent(in) :: itopo, jdim,is,js,ks
    real, intent(in) :: dzt(m1), dzm(m1), hw4(m1)
    real, intent(in) :: uc  (m1,m2,m3)
    real, intent(in) :: vc  (m1,m2,m3)
    real, intent(in) :: wc  (m1,m2,m3)
    real, intent(in) :: dn0 (m1,m2,m3)
    real, intent(in) :: dn0u(m1,m2,m3)
    real, intent(in) :: dn0v(m1,m2,m3)
    real, intent(in) :: dxt(m2,m3)
    real, intent(in) :: dxu(m2,m3)
    real, intent(in) :: dxv(m2,m3)
    real, intent(in) :: dyt(m2,m3)
    real, intent(in) :: dyu(m2,m3)
    real, intent(in) :: dyv(m2,m3)
    real, intent(in) :: rtgt(m2,m3)
    real, intent(in) :: rtgu(m2,m3)
    real, intent(in) :: rtgv(m2,m3)
    real, intent(in) :: f13t(m2,m3)
    real, intent(in) :: f23t(m2,m3)
    real, intent(in) :: fmapt(m2,m3)
    real, intent(in) :: fmapu(m2,m3)
    real, intent(in) :: fmapv(m2,m3)
    real, intent(in) :: fmapui(m2,m3)
    real, intent(in) :: fmapvi(m2,m3)

    real, dimension(m1,m2,m3),intent(inout) :: mfx_wind,mfy_wind,mfz_wind,&
                                               flxu,flxv,flxw

    ! Local Variables
    integer :: j,i,k,jm,im
    real :: c1z,c1x,c1y

    ! Compute momentum fluxes flxu, flxv, flxw
    !- only necessary one time per timestep
    !
    if(is == 1 .and. js == 0 .and. ks == 0) then

     do j = 1,m3
       do i = 1,m2
          do k = 1,m1
             flxu(k,i,j) = uc(k,i,j) * dn0u(k,i,j) * rtgu(i,j)  &
                         * fmapui(i,j)
             flxv(k,i,j) = vc(k,i,j) * dn0v(k,i,j) * rtgv(i,j)  &
                         * fmapvi(i,j)
          end do
       end do
     end do

     if(itopo.eq.0) then
       do j = 1,m3
          do i = 1,m2
             do k = 1,m1-1
                flxw(k,i,j) = wc(k,i,j)  &
                            * .5 * (dn0(k,i,j) + dn0(k+1,i,j))
             end do
          end do
       end do
     else
       do j = 1,m3
          jm = max(j-1,1)
          do i = 1,m2
             im = max(i-1,1)
             do k = 1,m1-1
                flxw(k,i,j) = wc(k,i,j)  &
                     * .5 * (dn0(k,i,j) + dn0(k+1,i,j))  &
                     + hw4(k) * ((flxu(k,i,j) + flxu(k+1,i,j)  &
                     + flxu(k,im,j) + flxu(k+1,im,j)) * f13t(i,j)  &
                     + (flxv(k,i,j) + flxv(k+1,i,j)  &
                     + flxv(k,i,jm) + flxv(k+1,i,jm)) * f23t(i,j))
             end do
          end do
       end do
     end if
    endif

    if(is == 1 .and. js == 0 .and. ks == 0) then
    !- Compute metric factors to compute U tendency
     do j = ja,jz
       do i = ia,izu
          c1z = 1. / rtgu(i,j)
          c1x = c1z * fmapu(i,j) * dxu(i,j)
          c1y = c1z * fmapu(i,j) * dyu(i,j)

          do k = 2,m1-1
              mfx_wind(k,i,j)= c1x / dn0u(k,i,j)
          end do

          do k = 2,m1-1
            mfy_wind(k,i,j)= c1y / dn0u(k,i,j)
          end do

          do k = 2,m1-1
            mfz_wind(k,i,j)= c1z * dzt(k) / dn0u(k,i,j)
          end do
       end do
     end do


    elseif(is == 0 .and. js == 1 .and. ks == 0) then
     !- Compute metric factors to compute V tendency
     do j = ja,jzv
       do i = ia,iz
          c1z = 1. / rtgv(i,j)
          c1x = c1z * fmapv(i,j) * dxv(i,j)
          c1y = c1z * fmapv(i,j) * dyv(i,j)

          do k = 2,m1-1
             mfx_wind(k,i,j) = c1x / dn0v(k,i,j)
          end do

          do k = 2,m1-1
            mfy_wind(k,i,j) = c1y / dn0v(k,i,j)
          end do

          do k = 2,m1-1
            mfz_wind(k,i,j)= c1z * dzt(k) / dn0v(k,i,j)
          end do
       end do
     end do

    elseif(is == 0 .and. js == 0 .and. ks == 1) then
     !- Compute metric factors to compute W tendency
     do j = ja,jz
       do i = ia,iz
          c1z = 1. / rtgt(i,j)
          c1x = c1z * fmapt(i,j) * dxt(i,j)
          c1y = c1z * fmapt(i,j) * dyt(i,j)

          do k = 2,m1-1
            mfx_wind(k,i,j) = c1x / (dn0(k,i,j) + dn0(k+1,i,j))
          end do

          do k = 2,m1-1
            mfy_wind(k,i,j) = c1y / (dn0(k,i,j) + dn0(k+1,i,j))
          end do

          do k = 2,m1-1
            mfz_wind(k,i,j) = c1z * dzm(k) / (dn0(k,i,j) + dn0(k+1,i,j))
          end do
       end do
     end do
    endif
  end subroutine mf_wind

  !!- compute interface values of scalars/wind
  subroutine compXYInterface_or1(mxp,myp,mzp,ks,isi,js,&
                                 mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
                                 mzppks,mxppis,myppjs, &
                                 scr,ufx_local,vfx_local,&
                                 border,qx,qy,variable, vname)
    implicit none
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: mzp
    integer,intent(in)  :: ks
    integer,intent(in)  :: isi
    integer,intent(in)  :: js
    integer, intent(in) :: mzi
    integer, intent(in) :: mzpp3
    integer, intent(in) :: mxi
    integer, intent(in) :: mxpp3
    integer, intent(in) :: myi
    integer, intent(in) :: mypp3
    integer, intent(in) :: mzppks
    integer, intent(in) :: mxppis
    integer, intent(in) :: myppjs
    real,intent(in)     :: scr      (mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    real,intent(in)     :: ufx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    real,intent(in)     :: vfx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    logical, intent(in) :: border(4)
    character(len=*), intent(in) :: vname
    logical, intent(in) :: variable
    real,intent(out)    :: qx(mzppks,mxppis,myppjs)
    real,intent(out)    :: qy(mzppks,mxppis,myppjs)

    integer, parameter :: west=1, east=2, north=3, south=4
    real, external :: flux_upwind
    integer :: i,j,k
    real :: dir

    !- compute x-interface values upwind order
    do j = 1,myp
      do i = 1,mxp-1
        do k = 1,mzp
          dir = SIGN(1.0,ufx_local(k,i,j)+ufx_local(k+ks,i+isi,j+js))
          qx(k,i,j)=flux_upwind(scr(k,i,j),scr(k,i+1,j),dir)
        enddo
      enddo
    enddo
    !- compute y-interface values upwind order
    do j = 1,myp-1
      do i = 1,mxp
        do k = 1,mzp
          dir = SIGN(1.0,vfx_local(k,i,j)+vfx_local(k+ks,i+isi,j+js))
          qy(k,i,j)=flux_upwind(scr(k,i,j),scr(k,i,j+1),dir)
        enddo
      enddo
    enddo

  end subroutine compXYInterface_or1

  subroutine compXYInterface_or2(mxp,myp,mzp,ks,isi,js,&
                                 mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
                                 mzppks,mxppis,myppjs, &
                                 scr,ufx_local,vfx_local,&
                                 border,qx,qy,variable, vname)
    implicit none
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: mzp
    integer,intent(in)  :: ks
    integer,intent(in)  :: isi
    integer,intent(in)  :: js
    integer, intent(in) :: mzi
    integer, intent(in) :: mzpp3
    integer, intent(in) :: mxi
    integer, intent(in) :: mxpp3
    integer, intent(in) :: myi
    integer, intent(in) :: mypp3
    integer, intent(in) :: mzppks
    integer, intent(in) :: mxppis
    integer, intent(in) :: myppjs
    real,intent(in)     :: scr      (mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    real,intent(in)     :: ufx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    real,intent(in)     :: vfx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    logical, intent(in) :: border(4)
    character(len=*), intent(in) :: vname
    logical, intent(in) :: variable
    real,intent(out)    :: qx(mzppks,mxppis,myppjs)
    real,intent(out)    :: qy(mzppks,mxppis,myppjs)

    integer, parameter :: west=1, east=2, north=3, south=4
    real, external :: fq2
    integer :: i,j,k

    !- compute x-interface values
    do j = 1,myp
      do i = 1,mxp-1
        do k = 1,mzp
          qx(k,i,j) = fq2(scr(k,i,j),scr(k,i+1,j))
        enddo
      enddo
    enddo
    !- compute y-interface values
    do j = 1,myp-1
      do i = 1,mxp
        do k = 1,mzp
          qy(k,i,j) = fq2(scr(k,i,j),scr(k,i,j+1))
        enddo
      enddo
    enddo

  end subroutine compXYInterface_or2

  !- compute x-interface values
  subroutine compXYInterface_or3(mxp,myp,mzp,ks,isi,js,&
                                 mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
                                 mzppks,mxppis,myppjs, &
                                 scr,ufx_local,vfx_local,&
                                 border,qx,qy,variable, vname)
    implicit none
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: mzp
    integer,intent(in)  :: ks
    integer,intent(in)  :: isi
    integer,intent(in)  :: js
    integer, intent(in) :: mzi
    integer, intent(in) :: mzpp3
    integer, intent(in) :: mxi
    integer, intent(in) :: mxpp3
    integer, intent(in) :: myi
    integer, intent(in) :: mypp3
    integer, intent(in) :: mzppks
    integer, intent(in) :: mxppis
    integer, intent(in) :: myppjs
    real,intent(in)     :: scr      (mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    real,intent(in)     :: ufx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    real,intent(in)     :: vfx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    logical, intent(in) :: border(4)
    character(len=*), intent(in) :: vname
    logical, intent(in) :: variable
    real,intent(out)    :: qx(mzppks,mxppis,myppjs)
    real,intent(out)    :: qy(mzppks,mxppis,myppjs)

    integer, parameter :: west=1, east=2, north=3, south=4
    real,external :: flux_upwind
    real,external :: fq2,fq3
    real :: dir
    integer :: i,j,k

    qx=0.0
    qy=0.0

    do j = 1,myp
      do i = 1,mxp-1
        do k = 1,mzp
          if(((border(west) .and. i==1) .or. (border(east) .and. i>mxp-2)) &
            .and. variable) then
            !Use order 1
            dir = SIGN(1.0,ufx_local(k,i,j)+ufx_local(k+ks,i+isi,j+js))
            qx(k,i,j)=flux_upwind(scr(k,i,j),scr(k,i+1,j),dir)
          else
            !use order 3
            dir = SIGN(1.0,ufx_local(k,i,j)+ufx_local(k+ks,i+isi,j+js))
            qx(k,i,j) = fq3(scr(k,i-1,j),scr(k,i,j),scr(k,i+1,j),scr(k,i+2,j),dir)
          endif
        enddo
      enddo
    enddo
    !- compute y-interface values
    do j = 1,myp-1
      do i = 1,mxp
        do k = 1,mzp
          if(((border(north) .and. j==1) .or. (border(south) .and. j>myp-2)) &
            .and. variable) then
            dir = SIGN(1.0,vfx_local(k,i,j)+vfx_local(k+ks,i+isi,j+js))
            qy(k,i,j)=flux_upwind(scr(k,i,j),scr(k,i,j+1),dir)
          else
            dir = SIGN(1.0,vfx_local(k,i,j)+vfx_local(k+ks,i+isi,j+js))
            qy(k,i,j) = fq3(scr(k,i,j-1),scr(k,i,j),scr(k,i,j+1), scr(k,i,j+2),dir)
          endif
        enddo
      enddo
    enddo

    !.. call dumpVarOnly1Lat(qx,vname//   '-x3-',25,4,1,mxp+1,1,myp+1,1,mzp+1,0.0,60.0)
    !.. call dumpVarOnly1Lat(qy,vname//   '-y3-',25,4,1,mxp+1,1,myp+1,1,mzp+1,0.0,60.0)

  end subroutine compXYInterface_or3

  !- compute x-interface values
  subroutine compXYInterface_or4(mxp,myp,mzp,ks,isi,js,&
                                 mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
                                 mzppks,mxppis,myppjs, &
                                 scr,ufx_local,vfx_local,&
                                 border,qx,qy,variable, vname)
    implicit none
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: mzp
    integer,intent(in)  :: ks
    integer,intent(in)  :: isi
    integer,intent(in)  :: js
    integer, intent(in) :: mzi
    integer, intent(in) :: mzpp3
    integer, intent(in) :: mxi
    integer, intent(in) :: mxpp3
    integer, intent(in) :: myi
    integer, intent(in) :: mypp3
    integer, intent(in) :: mzppks
    integer, intent(in) :: mxppis
    integer, intent(in) :: myppjs
    real,intent(in)     :: scr      (mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    real,intent(in)     :: ufx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    real,intent(in)     :: vfx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    logical, intent(in) :: border(4)
    character(len=*), intent(in) :: vname
    logical, intent(in) :: variable
    real,intent(out)    :: qx(mzppks,mxppis,myppjs)
    real,intent(out)    :: qy(mzppks,mxppis,myppjs)

    integer, parameter :: west=1, east=2, north=3, south=4
    real, external :: fq4
    integer :: i,j,k

    do j = 1,myp
      do i = 1,mxp-1
        do k = 1,mzp
          qx(k,i,j) = fq4(scr(k,i-1,j),scr(k,i,j),scr(k,i+1,j),scr(k,i+2,j))
        enddo
      enddo
    enddo
    !- compute y-interface values
    do j = 1,myp-1
      do i = 1,mxp
        do k = 1,mzp
          qy(k,i,j) = fq4(scr(k,i,j-1),scr(k,i,j),scr(k,i,j+1),scr(k,i,j+2))
        enddo
      enddo
    enddo

  end subroutine compXYInterface_or4

  !- compute x-interface values
  subroutine compXYInterface_or56(mxp,myp,mzp,ks,isi,js,&
                                 mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
                                 mzppks,mxppis,myppjs, &
                                 scr,ufx_local,vfx_local,&
                                 border,qx,qy,variable, vname, order_h)
    use advRkParam, only: fifth_order
    implicit none
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: mzp
    integer,intent(in)  :: ks
    integer,intent(in)  :: isi
    integer,intent(in)  :: js
    integer, intent(in) :: mzi
    integer, intent(in) :: mzpp3
    integer, intent(in) :: mxi
    integer, intent(in) :: mxpp3
    integer, intent(in) :: myi
    integer, intent(in) :: mypp3
    integer, intent(in) :: mzppks
    integer, intent(in) :: mxppis
    integer, intent(in) :: myppjs
    integer, intent(in) :: order_h
    real,intent(in)     :: scr      (mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    real,intent(in)     :: ufx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    real,intent(in)     :: vfx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    logical, intent(in) :: border(4)
    character(len=*), intent(in) :: vname
    logical, intent(in) :: variable
    real,intent(out)    :: qx(mzppks,mxppis,myppjs)
    real,intent(out)    :: qy(mzppks,mxppis,myppjs)

    integer, parameter :: west=1, east=2, north=3, south=4
    real, external :: fq,fq3,flux_upwind
    real :: dir
    integer :: i,j,k

    fifth_order = 1.0                ! 5th order
    if(order_h == 6)fifth_order = 0.0  ! 6th order

    do j = 1,myp
      do i = 1,mxp-1
        do k = 1,mzp
          if(((border(west) .and. i==1) .or. (border(east) .and. i==mxp-1)) &
            .and. variable) then
            !Order=1
            dir = SIGN(1.0,ufx_local(k,i,j)+ufx_local(k+ks,i+isi,j+js))
            qx(k,i,j)=flux_upwind(scr(k,i,j),scr(k,i+1,j),dir)
          elseif((border(west) .and. i==2) .or. (border(east) .and. i==mxp-2) .and. variable) then
            !use order 3
            dir = SIGN(1.0,ufx_local(k,i,j)+ufx_local(k+ks,i+isi,j+js))
            qx(k,i,j) = fq3(scr(k,i-1,j),scr(k,i,j),scr(k,i+1,j),scr(k,i+2,j),dir)
          else
            !Use order 5 or 6
            dir = SIGN(1.0,ufx_local(k,i,j)+ufx_local(k+ks,i+isi,j+js))
            qx(k,i,j) = fq(scr(k,i-2,j),scr(k,i-1,j),scr(k,i,j),scr(k,i+1,j),scr(k,i+2,j),scr(k,i+3,j),dir)
          endif
        enddo
      enddo
    enddo
    !- compute y-interface values
    do j = 1,myp-1
      do i = 1,mxp
        do k = 1,mzp
          if(((border(north) .and. j==1) .or. (border(south) .and. j==myp-1)) &
            .and. variable) then
            ! Order 1
            dir = SIGN(1.0,vfx_local(k,i,j)+vfx_local(k+ks,i+isi,j+js))
            qy(k,i,j)=flux_upwind(scr(k,i,j),scr(k,i,j+1),dir)
          elseif( (border(north) .and. j==2)  .or. (border(south) .and. j==myp-2) .and. variable) then
            !Order 3
            dir = SIGN(1.0,vfx_local(k,i,j)+vfx_local(k+ks,i+isi,j+js))
            qy(k,i,j) = fq3(scr(k,i,j-1),scr(k,i,j),scr(k,i,j+1), scr(k,i,j+2),dir)
          else
            !Use order 6
            dir = SIGN(1.0,vfx_local(k,i,j)+vfx_local(k+ks,i+isi,j+js))
            qy(k,i,j) = fq(scr(k,i,j-2),scr(k,i,j-1),scr(k,i,j),scr(k,i,j+1),scr(k,i,j+2),scr(k,i,j+3),dir)
          endif
        enddo
      enddo
    enddo

  end subroutine compXYInterface_or56

  !- compute z-interface values upwind order
  subroutine compZInterface_or1(mxp,myp,mzp,ks,isi,js,&
                                 mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
                                 mzppks,mxppis,myppjs, &
                                 scr,wfx_local,qz,variable,vname)
    implicit none
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: mzp
    integer,intent(in)  :: ks
    integer,intent(in)  :: isi
    integer,intent(in)  :: js
    integer, intent(in) :: mzi
    integer, intent(in) :: mzpp3
    integer, intent(in) :: mxi
    integer, intent(in) :: mxpp3
    integer, intent(in) :: myi
    integer, intent(in) :: mypp3
    integer, intent(in) :: mzppks
    integer, intent(in) :: mxppis
    integer, intent(in) :: myppjs
    real,intent(in)     :: scr      (mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    real,intent(in)     :: wfx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    character(len=*),intent(in)    :: vname
    logical, intent(in) :: variable
    real,intent(out)    :: qz(mzppks,mxppis,myppjs)

    integer, parameter :: west=1, east=2, north=3, south=4
    real, external :: flux_upwind
    integer :: i,j,k
    real :: dir

    do j = 1,myp
      do i = 1,mxp
        do k = 1,mzp-1
          dir = SIGN(1.0,wfx_local(k,i,j)+wfx_local(k+ks,i+isi,j+js))
          qz(k,i,j)=flux_upwind(scr(k,i,j),scr(k+1,i,j),dir)
        enddo
      enddo
    enddo

  end subroutine compZInterface_or1

  !- compute z-interface values
  subroutine compZInterface_or2(mxp,myp,mzp,ks,isi,js,&
                                 mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
                                 mzppks,mxppis,myppjs, &
                                 scr,wfx_local,qz,variable,vname)
    implicit none
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: mzp
    integer,intent(in)  :: ks
    integer,intent(in)  :: isi
    integer,intent(in)  :: js
    integer, intent(in) :: mzi
    integer, intent(in) :: mzpp3
    integer, intent(in) :: mxi
    integer, intent(in) :: mxpp3
    integer, intent(in) :: myi
    integer, intent(in) :: mypp3
    integer, intent(in) :: mzppks
    integer, intent(in) :: mxppis
    integer, intent(in) :: myppjs
    real,intent(in)     :: scr      (mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    real,intent(in)     :: wfx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    character(len=*),intent(in)    :: vname
    logical, intent(in) :: variable
    real,intent(out)    :: qz(mzppks,mxppis,myppjs)

    integer, parameter :: west=1, east=2, north=3, south=4
    real, external :: fq2
    integer :: i,j,k

    do j = 1,myp
      do i = 1,mxp
        do k = 1,mzp-1
          qz(k,i,j) = fq2(scr(k,i,j),scr(k+1,i,j))
        enddo
      enddo
    enddo

  end subroutine compZInterface_or2

  !- compute z-interface values
  subroutine compZInterface_or3(mxp,myp,mzp,ks,isi,js,&
                                 mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
                                 mzppks,mxppis,myppjs, &
                                 scr,wfx_local,qz,variable,vname)
    implicit none
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: mzp
    integer,intent(in)  :: ks
    integer,intent(in)  :: isi
    integer,intent(in)  :: js
    integer, intent(in) :: mzi
    integer, intent(in) :: mzpp3
    integer, intent(in) :: mxi
    integer, intent(in) :: mxpp3
    integer, intent(in) :: myi
    integer, intent(in) :: mypp3
    integer, intent(in) :: mzppks
    integer, intent(in) :: mxppis
    integer, intent(in) :: myppjs
    real,intent(in)     :: scr      (mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    real,intent(in)     :: wfx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    character(len=*),intent(in)    :: vname
    logical, intent(in) :: variable
    real,intent(out)    :: qz(mzppks,mxppis,myppjs)

    integer, parameter :: west=1, east=2, north=3, south=4
    real, external :: fq2,fq3,flux_upwind
    real :: dir
    integer :: i,j,k

    !qz=0.0

    do j = 1,myp
      do i = 1,mxp
        do k = 1,mzp-1
          if((k==1 .or. k>mzp-2) .and. variable) then
            ! order 1
            dir = SIGN(1.0,wfx_local(k,i,j)+wfx_local(k+ks,i+isi,j+js))
            qz(k,i,j)=flux_upwind(scr(k,i,j),scr(k+1,i,j),dir)
          else
            !Order 3
            dir = SIGN(1.0,wfx_local(k,i,j)+wfx_local(k+ks,i+isi,j+js))
            qz(k,i,j) = fq3(scr(k-1,i,j),scr(k,i,j),scr(k+1,i,j), scr(k+2,i,j),dir)
          ENDIF
        enddo
      enddo
    enddo

  end subroutine compZInterface_or3

  !- compute z-interface values
  subroutine compZInterface_or4(mxp,myp,mzp,ks,isi,js,&
                                 mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
                                 mzppks,mxppis,myppjs, &
                                 scr,wfx_local,qz,variable,vname)
    implicit none
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: mzp
    integer,intent(in)  :: ks
    integer,intent(in)  :: isi
    integer,intent(in)  :: js
    integer, intent(in) :: mzi
    integer, intent(in) :: mzpp3
    integer, intent(in) :: mxi
    integer, intent(in) :: mxpp3
    integer, intent(in) :: myi
    integer, intent(in) :: mypp3
    integer, intent(in) :: mzppks
    integer, intent(in) :: mxppis
    integer, intent(in) :: myppjs
    real,intent(in)     :: scr      (mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    real,intent(in)     :: wfx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    character(len=*),intent(in)    :: vname
    logical, intent(in) :: variable
    real,intent(out)    :: qz(mzppks,mxppis,myppjs)

    integer, parameter :: west=1, east=2, north=3, south=4
    real, external :: fq4,fq2
    integer :: i,j,k

    do j = 1,myp
      do i = 1,mxp
        do k = 1,mzp-1
          if(k==1 .or. k>mzp-2) then !use order 2
            qz(k,i,j) = fq2(scr(k,i,j),scr(k+1,i,j))
          else
            qz(k,i,j) = fq4(scr(k-1,i,j),scr(k,i,j),scr(k+1,i,j),scr(k+2,i,j))
          endif
        enddo
      enddo
    enddo

  end subroutine compZInterface_or4

  subroutine compZInterface_or56(mxp,myp,mzp,ks,isi,js,&
                                 mzi,mzpp3,mxi,mxpp3,myi,mypp3, &
                                 mzppks,mxppis,myppjs, &
                                 scr,wfx_local,qz,variable,vname,order_v)
    use advRkParam, only: fifth_order
    implicit none
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: mzp
    integer,intent(in)  :: ks
    integer,intent(in)  :: isi
    integer,intent(in)  :: js
    integer, intent(in) :: mzi
    integer, intent(in) :: mzpp3
    integer, intent(in) :: mxi
    integer, intent(in) :: mxpp3
    integer, intent(in) :: myi
    integer, intent(in) :: mypp3
    integer, intent(in) :: mzppks
    integer, intent(in) :: mxppis
    integer, intent(in) :: myppjs
    integer, intent(in) :: order_v
    real,intent(in)     :: scr      (mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    real,intent(in)     :: wfx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    character(len=*),intent(in)    :: vname
    logical, intent(in) :: variable
    real,intent(out)    :: qz(mzppks,mxppis,myppjs)

    integer, parameter :: west=1, east=2, north=3, south=4
    real, external :: fq,fq3,flux_upwind
    real :: dir
    integer :: i,j,k

    fifth_order = 1.0                   ! 5th order
    if(order_v == 6) fifth_order = 0.0  ! 6th order
    !- compute z-interface values
    do j = 1,myp
      do i = 1,mxp
        do k = 1,mzp-1
          if(k==1 .or. k==mzp-1) then
            ! order 1
            dir = SIGN(1.0,wfx_local(k,i,j)+wfx_local(k+ks,i+isi,j+js))
            qz(k,i,j)=flux_upwind(scr(k,i,j),scr(k+1,i,j),dir)
          elseif(k==2 .or. k==mzp-2) then
            !Order 3
            dir = SIGN(1.0,wfx_local(k,i,j)+wfx_local(k+ks,i+isi,j+js))
            qz(k,i,j) = fq3(scr(k-1,i,j),scr(k,i,j),scr(k+1,i,j), scr(k+2,i,j),dir)
          else
            !Use order 5
            dir = SIGN(1.0,wfx_local(k,i,j)+wfx_local(k+ks,i+isi,j+js))
            qz(k,i,j) = fq(scr(k-2,i,j),scr(k-1,i,j),scr(k,i,j),scr(k+1,i,j), &
                         scr(k+2,i,j),scr(k+3,i,j),dir)
          endif
        enddo
      enddo
    enddo

  end subroutine compZInterface_or56


  !positivity/monotonicity constraints
  subroutine posMonConstraints(mxp,myp,mzp,is,js,ks,ia,iz,ja,jz, &
                               pd_or_mnt_constraint, &
                               dt,ufx,vfx,wfx,vt3dh,vt3dj,vt3dk, scp, &
                               scr, ufx_local,vfx_local,wfx_local, &
                               qx,qy,qz,mzi,mzpp3,mxi,mxpp3,myi, &
                               mypp3,mzppks,mxppis,myppjs,mynum,vname,sct)
    use advRkParam, only: eps
    implicit none
    character(len=*), intent(in) :: vname
    integer, intent(in) :: mynum
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: mzp
    integer, intent(in) :: mzi
    integer, intent(in) :: mzpp3
    integer, intent(in) :: mxi
    integer, intent(in) :: mxpp3
    integer, intent(in) :: myi
    integer, intent(in) :: mypp3
    integer, intent(in) :: mzppks
    integer, intent(in) :: mxppis
    integer, intent(in) :: myppjs
    integer, intent(in) :: is
    integer, intent(in) :: js
    integer, intent(in) :: ks
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: pd_or_mnt_constraint
    real,    intent(in) :: dt
    real,    intent(in) :: ufx(mzp,mxp,myp)
    real,    intent(in) :: vfx(mzp,mxp,myp)
    real,    intent(in) :: wfx(mzp,mxp,myp)
    real,    intent(in) :: vt3dh(mzp,mxp,myp)
    real,    intent(in) :: vt3dj(mzp,mxp,myp)
    real,    intent(in) :: vt3dk(mzp,mxp,myp)
    real,    intent(in) :: scp(mzp,mxp,myp)
    real,    intent(in) :: sct(mzp,mxp,myp)
    real,    intent(in) :: scr      (mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    real,    intent(in) :: ufx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    real,    intent(in) :: vfx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    real,    intent(in) :: wfx_local(mzi:mzpp3,mxi:mxpp3,myi:mypp3)
    real,    intent(inout) :: qx(mzppks,mxppis,myppjs)
    real,    intent(inout) :: qy(mzppks,mxppis,myppjs)
    real,    intent(inout) :: qz(mzppks,mxppis,myppjs)

    real, external :: flux_upwind
    real,dimension(mzp+ks,mxp+is,myp+js) :: qxl,qyl,qzl
    integer :: i,j,k
    real :: dir,wtop,wbot,div_term,scale
    real :: ue,uw,vn,vs,scl_low,flux_out
    real :: flux_out_x,flux_out_y,flux_out_z

    !-compute x,y,z-interfaces values with upwind scheme
    do j = 1,myp
      do i = 1,mxp-1
        do k = 1,mzp
          dir = SIGN(1.0,ufx_local(k,i,j))
          !- upwind flux (pd and monotonic)
          qxl(k,i,j)=flux_upwind(scr(k,i,j),scr(k,i+1,j),dir)
          !- difference of high order and upwind fluxes
          qx (k,i,j)=qx (k,i,j)-qxl(k,i,j)
        enddo
      enddo
    enddo
    qxl(:,mxp,:)=0.
    !
    do j = 1,myp-1
      do i = 1,mxp
        do k = 1,mzp
          dir = SIGN(1.0,vfx_local(k,i,j))
          !- upwind flux (pd and monotonic)
          qyl(k,i,j)=flux_upwind(scr(k,i,j),scr(k,i,j+1),dir)
          !- difference of high order and upwind fluxes
          qy (k,i,j)=qy (k,i,j)-qyl(k,i,j)
        enddo
      enddo
    enddo
    qyl(:,:,myp)=0.
    !
    do j = 1,myp
      do i = 1,mxp
        do k = 1,mzp-1
          dir = SIGN(1.0,wfx_local(k,i,j))
          !- upwind flux (pd and monotonic)
          qzl(k,i,j)=flux_upwind(scr(k,i,j),scr(k+1,i,j),dir)
          !- difference of high order and upwind fluxes
          qz (k,i,j)=qz (k,i,j)-qzl(k,i,j)
        enddo
      enddo
    enddo
    qzl(mzp,:,:)=0.

    !-- this section only imposes positivity constraint on scalars
    IF(pd_or_mnt_constraint == 1) then
      do i = ia,iz
        do j = ja,jz
          do k = 2,mzp
            ue  =ufx(k,i  ,j)
            uw  =ufx(k,i-1,j)
            vn  =vfx(k,i,j  )
            vs  =vfx(k,i,j-1)
            wtop=wfx(k  ,i,j)
            wbot=wfx(k-1,i,j)
            div_term = scp(k,i,j)*(vt3dh(k,i,j)*(ue   - uw  ) + &
                                   vt3dj(k,i,j)*(vn   - vs  ) + &
                                   vt3dk(k,i,j)*(wtop - wbot))
            !- 1st order update
            scl_low = scp(k,i,j) +   &
                 dt*(- vt3dh(k,i,j)*(ue  *qxl(k,i,j) - uw  *qxl(k,i-1,j)) &
                     - vt3dj(k,i,j)*(vn  *qyl(k,i,j) - vs  *qyl(k,i,j-1)) &
                     - vt3dk(k,i,j)*(wtop*qzl(k,i,j) - wbot*qzl(k-1,i,j)) &
                     + div_term)
            !- net flux out
            flux_out_x = vt3dh(k,i,j)*(max(0.,ue  *qx(k,i,j)) - min(0.,uw  *qx(k,i-1,j)))
            flux_out_y = vt3dj(k,i,j)*(max(0.,vn  *qy(k,i,j)) - min(0.,vs  *qy(k,i,j-1)))
            flux_out_z = vt3dk(k,i,j)*(max(0.,wtop*qz(k,i,j)) - min(0.,wbot*qz(k-1,i,j)))
            flux_out = flux_out_x + flux_out_y + flux_out_z
            !- include divergence term
            flux_out = (flux_out - div_term)*dt
            !- re-scale the fluxes to keep positivity


            IF( flux_out > scl_low ) THEN
              !- scale factor (=< 1.)
              scale = max(0.,scl_low/(flux_out+eps))
              !-faces - x
              if( ue*qx(k,i  ,j) > 0.  ) qx(k,i  ,j)=scale*qx(k,i  ,j)
              if( uw*qx(k,i-1,j) < 0.  ) qx(k,i-1,j)=scale*qx(k,i-1,j)
              !-faces - y
              if( vn*qy(k,i  ,j) > 0.  ) qy(k,i  ,j)=scale*qy(k,i  ,j)
              if( vs*qy(k,i,j-1) < 0.  ) qy(k,i,j-1)=scale*qy(k,i,j-1)
              !-faces - z
              if( wtop*qz(k  ,i,j) > 0.) qz(k  ,i,j)=scale*qz(k  ,i,j)
              if( wbot*qz(k-1,i,j) < 0.) qz(k-1,i,j)=scale*qz(k-1,i,j)
            ENDIF
          enddo
        enddo
      enddo
    ENDIF ! endif for pd_or_mnt_constraint == 1

    !-- this section imposes monotonicity constraint on scalars
    !-- flux renormalization for monotonicity constraint following Durran (1998)
    IF(pd_or_mnt_constraint == 2) then
      stop "pd_or_mnt_constraint == 2 is not ready yet"
      !--- needs to:
      ! a) expand to 3 dimensions (z)
      ! b) change (i,j,k) +1 to (i,j,k) 
      ! c) change (i,j,k)    to (i,j,k) -1
      !
      !!- flux renormalization for monotonicity constraint following Durran (1998)
      !  DO j = 1+js,ny-1
      !    DO i = 1+is,nx-1
      !
      !      uw = u(i,j)   
      !      ue = u(i+1,j) 
      !      vs = v(i,j)   
      !      vn = v(i,j+1) 
      !      !====	
      !      !div = rdx * (ue - uw) + rdy * (vn - vs)
      !      !divplus = 1.0 - dt*div
      !      
      !      !- upwind/monotonic (1st order) update
      !      sn_td(i,j) =sn(i,j) - dt*( rdx * (ue*qxl(i+1,j) - uw*qxl(i,j))  &
      ! 			      + rdy * (vn*qyl(i,j+1) - vs*qyl(i,j))  )
      !      !
      !      !- in the future include the change associated to the wind divergence
      !      !scl_td=scl_td/divplus
      !     
      !      !-this is Pj+ antidiffusive fluxes into grid point (i,j)
      !      xflux_in(i,j) =  dt*( rdx*( max(0.,u(i  ,j) *qx (i  ,j))  &
      ! 				-min(0.,u(i+1,j) *qx (i+1,j)) )&
      ! 			 + rdy*( max(0.,v(i  ,j) *qy (i,j  ))  &
      ! 				-min(0.,v(i,j+1) *qy (i,j+1)) )) 
      !     
      !     !flux_out = flux_out-sn(i,j)*dt*div
      !      
      !      !-this is Pj-  antidiffusive fluxes out grid point (i,j)
      !      xflux_out(i,j)  =  dt*( rdx*( max(0.,u(i+1,j)*qx (i+1,j)) &
      ! 				  -min(0.,u(i  ,j)*qx (i  ,j)))&  
      ! 			   + rdy*( max(0.,v(i,j+1)*qy (i,j+1)) &
      ! 				  -min(0.,v(i,  j)*qy (i,j  )))) 
      !
      !    ENDDO
      !  ENDDO
      ! !-Durran 1998, pg 278
      !  DO j = 1+js,ny+js
      !    DO i = 1+is,nx+is
      !     !-permissible values of tracer concentration at time n+1
      !     smaxa(i,j)=max(s(i,j),sn_td(i,j)) 
      !     sminb(i,j)=min(s(i,j),sn_td(i,j)) 
      !    ENDDO
      !  ENDDO
      !
      !  DO j = 1+js,ny+js
      !    jm1=max(1 ,j-1)
      !    jp1=min(ny,j+1)
      !    DO i = 1+is,nx+is
      !     
      !     im1=max(1 ,i-1)
      !     ip1=min(nx,i+1)
      !     !-permissible values of tracer concentration at time n+1
      !     smax(i,j)=max(smaxa(i,j), smaxa(i,jm1), smaxa(i,jp1), &
      ! 			      smaxa(im1,j), smaxa(ip1,j)  )
      !     smin(i,j)=min(sminb(i,j), sminb(i,jm1), sminb(i,jp1), &
      ! 			      sminb(im1,j), sminb(ip1,j)  )
      !    ENDDO
      !  ENDDO
      !
      !  DO j = 1+js,ny-1
      !   DO i = 1+is,nx-1
      !    
      !    !- this is Rj+ the required limitation on the net antidiffusive flux _into_ 
      !    !- grid point (i,j). Attention: from here xflux_in becames Rj+
      !     IF( xflux_in (i,j) > 0.) &
      !     xflux_in (i,j)=min(1.,(smax(i,j)-sn_td(i,j))/(xflux_in(i,j)+eps))
      !
      !    !- this is Rj- the required limitation on the net antidiffusive flux _out of_
      !    !- grid point (i,j). Attention: from here xflux_out becames Rj-
      !     IF( xflux_out(i,j) > 0.) &
      !     xflux_out(i,j)=min(1.,(sn_td(i,j)-smin(i,j))/(xflux_out(i,j)+eps))
      !!-srf cm     
      !!**** by definition xflux_in and xflux_out are both >= 0. In the future, think 
      !!**** about eliminating the two IF statements above to speed up the code execution
      !    ENDDO
      !  ENDDO
      !
      !!--- optional preliminary step (not using,see page 260)
      !  !DO j = 1+js,ny-2
      !  !  DO i = 1+is,nx-2
      !    !if(u(i+1,j)*qx(i+1,j) * (sn_td(i+1,j) - sn_td(i,j)) < 0.  ) qx(i+1,j)=0.
      !    !if(u(i+1,j)*qx(i+1,j) * (sn_td(i+2,j) - sn_td(i+1,j)) < 0.) qx(i+1,j)=0.
      !    !if(u(i+1,j)*qx(i+1,j) * (sn_td(i,j) - sn_td(i-1,j)) < 0.  ) qx(i+1,j)=0.
      !    !if(v(i,j+1)*qy( ...
      !  !  ENDDO
      !  !ENDDO
      !!-------
      !  
      !  DO j = 1+js,ny
      !   DO i = 1+is+1,nx
      !     if(u(i,j) *qx (i,j) >= 0.) then
      ! 	 rmon_i(i,j) = min(xflux_in(i  ,j),xflux_out(i-1,j))
      !     else
      ! 	 rmon_i(i,j) = min(xflux_in(i-1,j),xflux_out(i  ,j))
      !     endif
      !   ENDDO 
      !  ENDDO
      !  rmon_i(1,:) = 0. ! b.c.
      !    
      !  DO j = 1+js+1,ny
      !   DO i = 1+is,nx
      !     if(v(i,j) *qy (i,j) >= 0.) then
      ! 	 rmon_j(i,j) = min(xflux_in(i,j  ),xflux_out(i,j-1))
      !     else
      ! 	 rmon_j(i,j) = min(xflux_in(i,j-1),xflux_out(i,j  ))
      !     endif
      !   ENDDO
      !  ENDDO
      !  rmon_j(:,1) = 0.0 ! b.c.
      !  
      !  !- setting rmon_ij = 0. => full upwind scheme
      !  !- setting rmon_ij = 1. => full high order scheme
      !  !-
      !  !- implements the monotonicity selectivity
      !  !- cm: no futuro reescrever toda esta se��o de modo que  somente seja 
      !  !-	feito todos estes calculos nos casos em que smoothx ou smoothy < smooth_max
      !  IF(mnt_selectivity == 1)THEN 
      !    DO j = 1+js,ny+js!-1
      !      DO i = 1+is,nx+is!-1
      !        if( smoothx(i,j) < smooth_max) rmon_i(i,j)=1.
      !        if( smoothy(i,j) < smooth_max) rmon_j(i,j)=1.  
      !!	if( log(smoothx(i,j)+1) < 0.02) rmon_i(i+1,j)=1.
      !!	if( log(smoothy(i,j)+1) < 0.02) rmon_j(i,j+1)=1.  
      !     ENDDO
      !    ENDDO
      !  ENDIF 
      !  !- get back the "qx fluxes", but now renormalized 
      !  DO j = 1+js,ny+js-1
      !    DO i = 1+is,nx+is
      !       qx(i,j) = rmon_i(i,j)* qx(i,j)
      !    ENDDO
      !  ENDDO
      !  DO j = 1+js,ny+js
      !    DO i = 1+is,nx+is-1
      !       qy(i,j) = rmon_j(i,j)* qy(i,j)
      !    ENDDO
      !  ENDDO
    ENDIF

    !-- this section imposes monotonicity constraint on scalars
    !-- flux renormalization of Blossey&Durran 2008 + Skamarock 2006
    IF(pd_or_mnt_constraint == 3) then
    !- reserved
    ENDIF ! endif for pd_or_mnt_constraint == 3
    !
    !-- get back the total "flux" , but now renormalized to guarantee
    !-- positivity and/or monotonicity
    do j = 1,myp
      do i = 1,mxp
        do k = 1,mzp-1
          qx(k,i,j) = qx(k,i,j) + qxl(k,i,j)
          qy(k,i,j) = qy(k,i,j) + qyl(k,i,j)
          qz(k,i,j) = qz(k,i,j) + qzl(k,i,j)
        enddo
      enddo
    enddo

  end subroutine posMonConstraints

  subroutine createTendency(mxp,myp,mzp,is,js,ks,ia,iz,ja,jz, &
                               mzppks,mxppis,myppjs, &
                               dt,ufx,vfx,wfx, &
                               vt3dh,vt3dj,vt3dk,scp, &
                               qx,qy,qz,sct,vname,mynum)
    use advRkParam, only: real_init
    implicit none
    integer, intent(in) :: mynum
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: mzp
    integer, intent(in) :: is
    integer, intent(in) :: js
    integer, intent(in) :: ks
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: mzppks
    integer, intent(in) :: mxppis
    integer, intent(in) :: myppjs
    real,    intent(in) :: dt
    real,    intent(in) :: vt3dh(mzp,mxp,myp)
    real,    intent(in) :: vt3dj(mzp,mxp,myp)
    real,    intent(in) :: vt3dk(mzp,mxp,myp)
    real,    intent(in) :: scp(mzp,mxp,myp)
    real,    intent(in) :: ufx(mzp,mxp,myp)!(-2:mzp+3,-2:mxp+3,-2:myp+3)
    real,    intent(in) :: vfx(mzp,mxp,myp)!(-2:mzp+3,-2:mxp+3,-2:myp+3)
    real,    intent(in) :: wfx(mzp,mxp,myp)!(-2:mzp+3,-2:mxp+3,-2:myp+3)
    real,    intent(in) :: qx(mzppks,mxppis,myppjs)
    real,    intent(in) :: qy(mzppks,mxppis,myppjs)
    real,    intent(in) :: qz(mzppks,mxppis,myppjs)
    character(len=*), intent(in) :: vname
    real,    intent(inout) :: sct(mzp,mxp,myp)

    integer :: i,j,k
    real :: dir,wtop,wbot,div_term
    real :: ue,uw,vn,vs
 
!call dumpVarAllLatLonk(qx,vname//'xS10',1,mxppis,1,myppjs,1,mzppks,0.0,60.0)
!call dumpVarAllLatLonk(qy,vname//'yS10',1,mxppis,1,myppjs,1,mzppks,0.0,60.0)
!call dumpVarAllLatLonk(qz,vname//'zS10',1,mxppis,1,myppjs,1,mzppks,0.0,60.0)

    do j = ja,jz
      do i = ia,iz
        do k = 2,mzp-1
          ue  =0.5*( ufx(k,i  ,j) + ufx(k+ks,i  +is,j+js) )
          uw  =0.5*( ufx(k,i-1,j) + ufx(k+ks,i-1+is,j+js) )
          vn  =0.5*( vfx(k,i,j  ) + vfx(k+ks,i+is,j  +js) )
          vs  =0.5*( vfx(k,i,j-1) + vfx(k+ks,i+is,j-1+js) )
          wtop=0.5*( wfx(k  ,i,j) + wfx(k+ks  ,i+is,j+js) )
          wbot=0.5*( wfx(k-1,i,j) + wfx(k-1+ks,i+is,j+js) )
          div_term = scp(k,i,j)*(vt3dh(k,i,j)*(ue   - uw  ) + &
                                 vt3dj(k,i,j)*(vn   - vs  ) + &
                                 vt3dk(k,i,j)*(wtop - wbot))
          sct(k,i,j) = sct(k,i,j) &
                    - vt3dh(k,i,j)*(ue  *qx(k,i,j) - uw  *qx(k,i-1,j)) &
                    - vt3dj(k,i,j)*(vn  *qy(k,i,j) - vs  *qy(k,i,j-1)) &
                    - vt3dk(k,i,j)*(wtop*qz(k,i,j) - wbot*qz(k-1,i,j)) &
                    + div_term
        end do
      end do
    end do

!call dumpVarAllLatLonk(sct,vname//'tS20',1,mxp,1,myp,1,mzp,0.0,60.0)

  end subroutine createTendency

  !--- 1st order or upwind scheme
  REAL FUNCTION flux_upwind(qi, qip1, dir )
    implicit none
    real, intent(in) :: dir
    real, intent(in) :: qi
    real, intent(in) :: qip1

    flux_upwind = 0.5*(1.+dir)*qi - 0.5*(dir-1.)*qip1

  END FUNCTION flux_upwind

  !--- 2nd order interpolation operator
  real function fq2(qi,qip1)
    implicit none
    real, intent(in) :: qi
    real, intent(in) :: qip1

    fq2= 0.5*(qip1+qi)

  end function fq2

  !--- 3rd order interpolation operator
  real function fq3(qim1,qi,qip1,qip2,dir)
    use advRkParam, only: f40,f41
    implicit none
    real,intent(in) :: qim1
    real,intent(in) :: qi
    real,intent(in) :: qip1
    real,intent(in) :: qip2
    real,intent(in) :: dir

    fq3= f40*(qip1+qi)-f41*(qip2+qim1) &
        - f41*(3.*(qip1-qi)-(qip2-qim1))*dir

  end function fq3

  !--- 4th order interpolation operator
  real function fq4(qim1,qi,qip1,qip2)
    use advRkParam, only: f40,f41
    implicit none
    real, intent(in) :: qim1
    real, intent(in) :: qi
    real, intent(in) :: qip1
    real, intent(in) :: qip2

    fq4= f40*(qip1+qi)-f41*(qip2+qim1)

  end function fq4

  !- polynomial interpolation operators
  ! for 5th and 6th orders
  real function fq(qim2,qim1,qi,qip1,qip2,qip3,dir)
    use advRkParam, only: f50,f51,f52,fifth_order
    implicit none
    real, intent(in) :: qim2
    real, intent(in) :: qim1
    real, intent(in) :: qi
    real, intent(in) :: qip1
    real, intent(in) :: qip2
    real, intent(in) :: qip3
    real, intent(in) :: dir

    fq = f50*(qip1 + qi) - f51*(qip2 + qim1) + f52*(qip3 + qim2) &
     - fifth_order* f52*(qip3-qim2-5.*(qip2-qim1)+10.*(qip1-qi))*dir

  end function fq

  subroutine dumpXYvar(var,vname,fase,ii,ie,ji,je,ki,ke,begtime,endtime)
    use mem_grid, only: time
    use node_mod, only:  &
       nodei0, &
       nodej0, mynum
    implicit none

    real,intent(in) :: var(ki:ke,ii:ie,ji:je)
    character(len=*), intent(in) :: vname
    character(len=*), intent(in) :: fase
    integer, intent(in) :: ii,ie,ji,je,ki,ke
    real, intent(in) :: begtime
    real, intent(in) :: endtime

    character(len=6) :: ctim
    character(len=2) :: cmzp,cmxp
    character :: cnum,cmyn
    character(len=256) :: fname
    integer :: i,j,ia,ja,k,ios

    if(time>endtime .or. time<begtime) return

    write(ctim,fmt='(I3.3)') int(time)
    write(cmyn,fmt='(I1)') mynum

    write(cmxp,fmt='(I2.2)') ie-ii+1
    fname='dumpDir/'//trim(vname)//'_'//trim(fase)//'_'//cmyn//'_'//trim(ctim)//'.dat'
    !print *,'file=',trim(fname)

    open(unit=65,file=trim(fname))
    !if(ios/=0) then
    !  print *,'Err: ',ios
    !  stop 'Error open file'
    !end if
    do k=ki,ke
      write(65,fmt='(A,I2.2,A)') '--------------- k=',k,' -------------------'
      write(65,fmt='(1(A6,2X),'//cmxp//'(I16,1X))') &
        '',(i,i=ii,ie)
      write(65,fmt='(A)') ''
      do j=ji,je
          write(65,fmt='(I6,2X,'//cmxp//'(E16.8,1X),I4)') &
            j,(var(k,i,j),i=ii,ie),j+nodej0(mynum,1)
      enddo
      write(65,fmt='(A)') ''
      write(65,fmt='(1(A6,2X),'//cmxp//'(I16,1X))') &
        '',(i+nodei0(mynum,1),i=ii,ie)
      call flush(65)
    end do
    close(65)

  end subroutine dumpXYvar

  subroutine dumpVarOnly1Col(var,varn,ipos,jpos,xi,xe,yi,ye,zi,ze,begtime,endtime)
    use mem_grid, only: time
    use node_mod, only:  &
       nodei0, &
       nodej0, mynum, nmachs
    implicit none

    integer         , intent(in) :: ipos,jpos,xi,xe,yi,ye,zi,ze
    real            , intent(in) :: begtime,endtime
    real            , intent(in) :: var(zi:ze,xi:xe,yi:ye)
    character(len=*), intent(in) :: varn

    integer :: i,j
    character :: cnm,cmyp
    character(len=3) :: ctim
    character(len=2) :: cmzp

    write(cnm,fmt='(I1)') nmachs
    write(cmyp,fmt='(I1)') mynum
    write(ctim,fmt='(I3.3)') int(time)
    write(cmzp,fmt='(I2.2)') ze-zi+1

    if(time<begtime .or. time>endtime) return

    i= ipos-nodei0(mynum,1)
    j= jpos-nodej0(mynum,1)

    open(unit=88,file='dumpDir/'//trim(varn)//trim(cnm)//trim(cmyp)//'-'//trim(ctim)//'.dat')
    if(i<xi .or. i>xe) then
      write (88,fmt='(A,I2,A,I2,A,I1,A,I2,1X,I2)') 'Ipos=',ipos,', local=',i,' is out of P',mynum,' area. Min,Max=',xi,xe
      return
    end if
    if(j<yi .or. j>ye) then
      write (88,fmt='(A,I2,A,I2,A,I1,A,I2,1X,I2)') 'Jpos=',jpos,', local=',j,' is out of P',mynum,' area. Min,Max=',yi,ye
      return
    end if
    write(88,fmt='(4(A,I2.2,1X),A,'//cmzp//'(E16.6,1X))') 'iPos=',ipos,',jPos=',jpos,',i=',i,',j=',j, &
             ',Var(zi:ze)=',var(zi:ze,i,j)
    close(88)

  end subroutine dumpVarOnly1Col

  subroutine dumpVarOnly1Lat(var,varn,jpos,kpos,xi,xe,yi,ye,zi,ze,begtime,endtime)
    use mem_grid, only: time
    use node_mod, only:  &
       nodei0, &
       nodej0, mynum, nmachs
    implicit none

    integer         , intent(in) :: jpos,kpos,xi,xe,yi,ye,zi,ze
    real            , intent(in) :: begtime,endtime
    real            , intent(in) :: var(zi:ze,xi:xe,yi:ye)
    character(len=*), intent(in) :: varn

    integer :: i,j,iaux
    character :: cnm,cmyp
    character(len=3) :: ctim
    character(len=2) :: cmzp,cjp,ckp

    write(cnm,fmt='(I1)') nmachs
    write(cmyp,fmt='(I1)') mynum
    write(ctim,fmt='(I3.3)') int(time)
		write(cjp,fmt='(I2.2)') jpos
    write(ckp,fmt='(I2.2)') kpos

    if(time<begtime .or. time>endtime) return

    j= jpos-nodej0(mynum,1)
    if(j<yi .or. j>ye) then
      write (88,fmt='(A,I2,A,I2,A,I1,A,I2,1X,I2)') 'Jpos=',jpos,', local=',j,' is out of P',mynum,' area. Min,Max=',yi,ye
      return
    end if
    if(nmachs==1) then
      iaux=(xe-xi+1)/2+1
      write(cmzp,fmt='(I2.2)') iaux
      open(unit=88,file='dumpDir/'//trim(varn)//'J'//cjp//'K'//ckp//trim(cnm)//'1'//'-'//trim(ctim)//'.dat')
      write(88,fmt='(3(A,I3,1X),A,I3,A,I3,A,'//cmzp//'(E16.6,1X))') 'kPos=',kpos,',jPos=',jpos,'j=',j, &
             'var(',xi+nodei0(mynum,1),':',iaux+nodei0(mynum,1),')=',var(kpos,xi:iaux,j)
      close(88)
      open(unit=88,file='dumpDir/'//trim(varn)//'J'//cjp//'K'//ckp//trim(cnm)//'2'//'-'//trim(ctim)//'.dat')
      write(88,fmt='(3(A,I3,1X),A,I3,A,I3,A,'//cmzp//'(E16.6,1X))') 'kPos=',kpos,',jPos=',jpos,'j=',j, &
             'var(',(iaux-1)+nodei0(mynum,1),':',xe+nodei0(mynum,1),')=',var(kpos,iaux-1:xe,j)
      close(88)

    else
      write(cmzp,fmt='(I2.2)') xe-xi+1
      open(unit=88,file='dumpDir/'//trim(varn)//'J'//cjp//'K'//ckp//trim(cnm)//trim(cmyp)//'-'//trim(ctim)//'.dat')
      if(j<yi .or. j>ye) then
        write (88,fmt='(A,I2,A,I2,A,I1,A,I2,1X,I2)') 'Jpos=',jpos,', local=',j,' is out of P',mynum,' area. Min,Max=',yi,ye
        return
      end if
      write(88,fmt='(3(A,I3,1X),A,I3,A,I3,A,'//cmzp//'(E16.6,1X))') 'kPos=',kpos,',jPos=',jpos,'j=',j, &
             'var(',xi+nodei0(mynum,1),':',xe+nodei0(mynum,1),')=',var(kpos,xi:xe,j)
      close(88)
    endif


  end subroutine dumpVarOnly1Lat

  subroutine dumpVarOnly1Lon(var,varn,ipos,kpos,xi,xe,yi,ye,zi,ze,begtime,endtime)
    use mem_grid, only: time
    use node_mod, only:  &
       nodei0, &
       nodej0, mynum, nmachs
    implicit none

    integer         , intent(in) :: ipos,kpos,xi,xe,yi,ye,zi,ze
    real            , intent(in) :: begtime,endtime
    real            , intent(in) :: var(zi:ze,xi:xe,yi:ye)
    character(len=*), intent(in) :: varn

    integer :: i,j,iaux
    character :: cnm,cmyp
    character(len=3) :: ctim
    character(len=2) :: cmzp,cip,ckp

    write(cnm,fmt='(I1)') nmachs
    write(cmyp,fmt='(I1)') mynum
    write(ctim,fmt='(I3.3)') int(time)
		write(cip,fmt='(I2.2)') ipos
    write(ckp,fmt='(I2.2)') kpos

    if(time<begtime .or. time>endtime) return

    i= ipos-nodei0(mynum,1)
    if(i<xi .or. i>xe) then
      write (88,fmt='(A,I2,A,I2,A,I1,A,I2,1X,I2)') 'Ipos=',ipos,', local=',i,' is out of P',mynum,' area. Min,Max=',xi,xe
      return
    end if
    if(nmachs==1) then
      iaux=(ye-yi+1)/2+1
      write(cmzp,fmt='(I2.2)') iaux
      open(unit=88,file='dumpDir/'//trim(varn)//'I'//cip//'K'//ckp//trim(cnm)//'1'//'-'//trim(ctim)//'.dat')
      write(88,fmt='(3(A,I3,1X),A,I3,A,I3,A,'//cmzp//'(E16.6,1X))') 'kPos=',kpos,',iPos=',ipos,'i=',i, &
             'var(',yi+nodej0(mynum,1),':',iaux+nodej0(mynum,1),')=',var(kpos,i,yi:iaux)
      close(88)
      open(unit=88,file='dumpDir/'//trim(varn)//'I'//cip//'K'//ckp//trim(cnm)//'2'//'-'//trim(ctim)//'.dat')
      write(88,fmt='(3(A,I3,1X),A,I3,A,I3,A,'//cmzp//'(E16.6,1X))') 'kPos=',kpos,',iPos=',ipos,'i=',i, &
             'var(',(iaux-1)+nodej0(mynum,1),':',ye+nodej0(mynum,1),')=',var(kpos,i,iaux-1:ye)
      close(88)

    else
      write(cmzp,fmt='(I2.2)') ye-yi+1
      open(unit=88,file='dumpDir/'//trim(varn)//'J'//cip//'K'//ckp//trim(cnm)//trim(cmyp)//'-'//trim(ctim)//'.dat')
      if(j<yi .or. j>ye) then
        write (88,fmt='(A,I2,A,I2,A,I1,A,I2,1X,I2)') 'Ipos=',ipos,', local=',i,' is out of P',mynum,' area. Min,Max=',xi,xe
        return
      end if
			close(88)
			open(unit=88,file='dumpDir/'//trim(varn)//'I'//cip//'K'//ckp//trim(cnm)//trim(cmyp)//'-'//trim(ctim)//'.dat')
      write(88,fmt='(3(A,I3,1X),A,I3,A,I3,A,'//cmzp//'(E16.6,1X))') 'kPos=',kpos,',iPos=',ipos,'i=',i, &
             'var(',yi+nodej0(mynum,1),':',ye+nodej0(mynum,1),')=',var(kpos,i,yi:ye)
      close(88)
    endif


  end subroutine dumpVarOnly1Lon

  subroutine dumpVarAllLatLonk(var,varn,ilin,iter,ipos,xi,xe,yi,ye,zi,ze,begtime,endtime,rout)
    use mem_grid, only: time
    use node_mod, only:  &
       nodei0, &
       nodej0, mynum, nmachs
    implicit none

    integer         , intent(in) :: xi,xe,yi,ye,zi,ze
    integer         , intent(in) :: iter,ipos,ilin
    real            , intent(in) :: begtime,endtime
    real            , intent(in) :: var(zi:ze,xi:xe,yi:ye)
    character(len=*), intent(in) :: varn
    character(len=*), intent(in) :: rout

    integer :: i,j,iaux
    character :: cnm,cmyp,citer,cipos
    character(len=4) :: ctim
    character(len=2) :: cmzp,cip,ckp
    character(len=4) :: cilin
    character(len=48) :: header

    if(time<begtime .or. time>endtime) return

    write(cnm,fmt='(I1)') nmachs
    write(cmyp,fmt='(I1)') mynum
    write(ctim,fmt='(I4.4)') int(time)
    write(cmzp,fmt='(I2.2)') ze-zi+1

    write(citer,fmt='(I1)') iter
    write(cipos,fmt='(I1)') ipos
    write(cilin,fmt='(I4.4)') ilin

    iaux=len(rout)-3

    header='dumpDir/'//rout(4:iaux)//'.'

    if(nmachs==1) then
      iaux=(ye-yi+1)/2+1
      open(unit=88,file=trim(header)//trim(varn)//'.L'//cilin//'.I'//citer//'.N'//cipos//'.P'&
      //trim(cnm)//'1'//'.T'//trim(ctim)//'.dat')
      do i=1,iaux
        do j=yi,ye
          write(88,fmt='(2(I2.2,1X),'//cmzp//'(E16.6))') i+nodei0(mynum,1),j+nodej0(mynum,1),var(zi:ze,i,j)
        enddo
      enddo
      close(88)
      open(unit=88,file=trim(header)//trim(varn)//'.L'//cilin//'.I'//citer//'.N'//cipos//'.P'&
      //trim(cnm)//'2'//'.T'//trim(ctim)//'.dat')
      do i=iaux-1,xe
        do j=yi,ye
          write(88,fmt='(2(I2.2,1X),'//cmzp//'(E16.6))') i+nodei0(mynum,1),j+nodej0(mynum,1),var(zi:ze,i,j)
        enddo
      enddo
      close(88)
    else
      open(unit=88,file=trim(header)//trim(varn)//'.L'//cilin//'.I'//citer//'.N'//cipos//'.P'&
      //trim(cnm)//trim(cmyp)//'.T'//trim(ctim)//'.dat')
      do i=xi,xe
        do j=yi,ye
          write(88,fmt='(2(I2.2,1X),'//cmzp//'(E16.6))') i+nodei0(mynum,1),j+nodej0(mynum,1),var(zi:ze,i,j)
        enddo
      enddo
      close(88)
    endif

  end subroutine dumpVarAllLatLonk

  subroutine dumpVarAllLatLonkSurf(var,varn,ilin,iter,ipos,xi,xe,yi,ye,zi,ze,begtime,endtime,rout)
    use mem_grid, only: time
    use node_mod, only:  &
       nodei0, &
       nodej0, mynum, nmachs
    implicit none

    integer         , intent(in) :: xi,xe,yi,ye,zi,ze
    integer         , intent(in) :: iter,ipos,ilin
    real            , intent(in) :: begtime,endtime
    real            , intent(in) :: var(xi:xe,yi:ye,zi:ze)
    character(len=*), intent(in) :: varn
    character(len=*), intent(in) :: rout

    integer :: i,j,iaux
    character :: cnm,cmyp,citer,cipos
    character(len=4) :: ctim
    character(len=2) :: cmzp,cip,ckp
    character(len=4) :: cilin
    character(len=48) :: header

    if(time<begtime .or. time>endtime) return

    write(cnm,fmt='(I1)') nmachs
    write(cmyp,fmt='(I1)') mynum
    write(ctim,fmt='(I4.4)') int(time)
    write(cmzp,fmt='(I2.2)') ze-zi+1

    write(citer,fmt='(I1)') iter
    write(cipos,fmt='(I1)') ipos
    write(cilin,fmt='(I4.4)') ilin

    iaux=len(rout)-3

    header='dumpDir/'//rout(4:iaux)//'.'

    if(nmachs==1) then
      iaux=(ye-yi+1)/2+1
      open(unit=88,file=trim(header)//trim(varn)//'.L'//cilin//'.I'//citer//'.N'//cipos//'.P'&
      //trim(cnm)//'1'//'.T'//trim(ctim)//'.dat')
      do i=1,iaux
        do j=yi,ye
          write(88,fmt='(2(I2.2,1X),'//cmzp//'(E16.6))') i+nodei0(mynum,1),j+nodej0(mynum,1),var(i,j,zi:ze)
        enddo
      enddo
      close(88)
      open(unit=88,file=trim(header)//trim(varn)//'.L'//cilin//'.I'//citer//'.N'//cipos//'.P' &
      //trim(cnm)//'2'//'.T'//trim(ctim)//'.dat')
      do i=iaux-1,xe
        do j=yi,ye
          write(88,fmt='(2(I2.2,1X),'//cmzp//'(E16.6))') i+nodei0(mynum,1),j+nodej0(mynum,1),var(i,j,zi:ze)
        enddo
      enddo
      close(88)
    else
      open(unit=88,file=trim(header)//trim(varn)//'.L'//cilin//'.I'//citer//'.N'//cipos//'.P' &
      //trim(cnm)//trim(cmyp)//'.T'//trim(ctim)//'.dat')
      do i=xi,xe
        do j=yi,ye
          write(88,fmt='(2(I2.2,1X),'//cmzp//'(E16.6))') i+nodei0(mynum,1),j+nodej0(mynum,1),var(i,j,zi:ze)
        enddo
      enddo
      close(88)
    endif

   end subroutine dumpVarAllLatLonkSurf
!     *********************************************************************

   subroutine dumpVarAllLatLonk3P(var,varn,ilin,iter,ipos,xi,xe,yi,ye,zi,ze,begtime,endtime,rout)
    use mem_grid, only: time
    use node_mod, only:  &
       nodei0, &
       nodej0, mynum, nmachs
    implicit none

    integer         , intent(in) :: xi,xe,yi,ye,zi,ze
    integer         , intent(in) :: iter,ipos,ilin
    real            , intent(in) :: begtime,endtime
    real            , intent(in) :: var(zi:ze,xi:xe,yi:ye)
    character(len=*), intent(in) :: varn
    character(len=*), intent(in) :: rout

    integer :: i,j,iaux
    character :: cnm,cmyp,citer,cipos
    character(len=4) :: ctim
    character(len=2) :: cmzp,cip,ckp
    character(len=4) :: cilin
    character(len=48) :: header

    if(time<begtime .or. time>endtime) return

    write(cnm,fmt='(I1)') nmachs
    write(cmyp,fmt='(I1)') mynum
    write(ctim,fmt='(I4.4)') int(time)
    write(cmzp,fmt='(I2.2)') ze-zi+1

    write(citer,fmt='(I1)') iter
    write(cipos,fmt='(I1)') ipos
    write(cilin,fmt='(I4.4)') ilin

    iaux=len(rout)-3

    header='dumpDir/'//rout(4:iaux)//'.'

    if(nmachs==1) then
      iaux=(ye-yi+1)/2+1
      open(unit=88,file=trim(header)//trim(varn)//'.L'//cilin//'.I'//citer//'.N'//cipos//'.P' &
      //trim(cnm)//'1'//'.T'//trim(ctim)//'.dat')
      do i=1,31
        do j=1,41
          write(88,fmt='(2(I2.2,1X),'//cmzp//'(E16.6))') i+nodei0(mynum,1),j+nodej0(mynum,1),var(zi:ze,i,j)
        enddo
      enddo
      close(88)

      open(unit=88,file=trim(header)//trim(varn)//'.L'//cilin//'.I'//citer//'.N'//cipos//'.P' &
      //trim(cnm)//'2'//'.T'//trim(ctim)//'.dat')
      do i=30,60
        do j=1,41
          write(88,fmt='(2(I2.2,1X),'//cmzp//'(E16.6))') i+nodei0(mynum,1),j+nodej0(mynum,1),var(zi:ze,i,j)
        enddo
      enddo
      close(88)

      open(unit=88,file=trim(header)//trim(varn)//'.L'//cilin//'.I'//citer//'.N'//cipos//'.P' &
      //trim(cnm)//'3'//'.T'//trim(ctim)//'.dat')
      do i=1,60
        do j=40,60
          write(88,fmt='(2(I2.2,1X),'//cmzp//'(E16.6))') i+nodei0(mynum,1),j+nodej0(mynum,1),var(zi:ze,i,j)
        enddo
      enddo
      close(88)

    else

      open(unit=88,file=trim(header)//trim(varn)//'.L'//cilin//'.I'//citer//'.N'//cipos//'.P' &
      //trim(cnm)//trim(cmyp)//'.T'//trim(ctim)//'.dat')
      do i=xi,xe
        do j=yi,ye
          write(88,fmt='(2(I2.2,1X),'//cmzp//'(E16.6))') i+nodei0(mynum,1),j+nodej0(mynum,1),var(zi:ze,i,j)
        enddo
      enddo
      close(88)
    endif

  end subroutine dumpVarAllLatLonk3P




!end module ModAdvectc_rk
