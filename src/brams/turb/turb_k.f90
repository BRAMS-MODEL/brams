!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine diffuse()

  ! +-----------------------------------------------------------------+
  ! \	this routine is the subdriver to compute tendencies due to    \
  ! \	  subgrid-scale turbulence.				      \
  ! +-----------------------------------------------------------------+

  use mem_tend, only:    &
       tend                   ! %tket, %epst, %ut, %vt, %wt

  use mem_basic, only:   &
       basic_g                !  %up(IN), %vp(IN)

  use var_tables, only:  &
       num_scalar,       &    ! INTENT(IN)
       scalar_tab             ! %var_p, %var_t

  use mem_turb, only:    &
       idiffk,           &    !INTENT(IN)
       turb_g,           &    ! %tkep, %hkm, %vkh
       xkhkm

  use mem_grid, only:     &
       if_adap,           &  !           INTENT(IN)
       jdim,              &  !           INTENT(IN)
       ngrid,             &  !           INTENT(IN)
       grid_g,            &  ! %rtgt     INTENT(IN)
       dzm,               &  !           INTENT(IN)
       dzt,               &  !           INTENT(IN)
       npatch,            &  !           INTENT(IN)
       nstbot,            &  !           INTENT(IN)
       nscl,              &  !           INTENT(IN)
       naddsc,            &  !           INTENT(IN)
       dtlt,            &  !           INTENT(IN)
       nnzp

  use mem_leaf, only:     &
       leaf_g                  !INTENT(IN)

  use mem_micro, only:    &
       micro_g                 !%rcp

  use mem_scratch, only:  &
       scratch,           &    ! %vt3da, %vt3db, %vt3dc, %vt3dd, %vt3de,
                                ! %vt3df, %vt3dg, %vt3dh, %vt3di, %vt3dj,
                                ! %vt3dn, %scr2
       vctr34                  !

  use node_mod, only:     &
       mxp,               &  !INTENT(IN)
       myp,               &  !INTENT(IN)
       mzp,               &  !INTENT(IN)
       ia,                &  !INTENT(IN)
       iz,                &  !INTENT(IN)
       ja,                &  !INTENT(IN)
       jz,                &  !INTENT(IN)
       ia_1,              &  !INTENT(IN)
       ja_1,              &  !INTENT(IN)
       iz1,               &  !INTENT(IN)
       jz1,               &  !INTENT(IN)
       ibcon,             &  !INTENT(IN)
       mynum,             &  !INTENT(IN)
       nodei0,            &  !INTENT(IN)
       nodej0,            &  !INTENT(IN)
       nodemyp,           &  !INTENT(IN)
       nodemxp,           &  !INTENT(IN)
       izu,               &  !INTENT(IN)
       jzv,               &  !INTENT(IN)
       ia1,               &  !INTENT(IN)
       ja1,               &  !INTENT(IN)
       iz_1,              &  !INTENT(IN)
       jz_1,              &  !INTENT(IN)
       mynum

  use ke_coms, only:      &
       alf_eps,           &
       alf_tke

  use micphys, only:      &
       level

  use mem_turb_scalar, only:   &
       turb_s

  use mem_grell,  only:        &
       cuforc_sh_g

  use mem_cuparm, only : &
       nnqparm

  use mem_chem1, only: &
       chemistry, &
       nspecies_transported

  !ml/srf- for stilt/new turb scheme
  use mem_stilt, only: &
  	stilt_g,       &
	imassflx

  use ccatt_start, only:        &
       ccatt           ! intent(in)


  implicit none
  include "i8.h"
  !local variables:
  integer(kind=i8) :: mxyzp, ind
  integer :: n
  real :: s1,s2,s3
  real, pointer :: scalarp(:), scalart(:), vkh_p(:), hkh_p(:)
  integer :: i,j,k,ksf

  ! srf - Large Scale Forcing for shallow and deep cumulus
  !real, pointer :: lsfcupar_p(:,:,:)

  real, target :: scr3(mxp*myp*mzp)
  real, target :: scr2(mxp*myp*mzp)


  !interface
  !   subroutine PBLforcing(ngrid, m1, m2, m3, ia, iz, ja, jz, &
  !        vt3df, scp, lsfcupar, nsc)
  !     ! Arguments:
  !    integer, INTENT(IN) :: ngrid, m1, m2, m3, ia, iz, ja, jz
  !     real, INTENT(IN)    :: vt3df(m1,m2,m3), scp(m1,m2,m3)
  !     real, pointer       :: lsfcupar(:,:,:)
  !     integer, intent(IN) :: nsc
  !   end subroutine PBLforcing
  !end interface

  ! CATT
  ! Nullifing pointer to Large Scale Forcing for GRELL CUPAR - Not used
  !nullify(lsfcupar_p)
  !

  mxyzp = mxp*myp*mzp

  scr2 = 0.0
  scr3 = 0.0
  scratch%vt3dg = 0. !LFR - Avoiding overflow in strain
  !-srf - opt turb_k com nakanishi
  !-srf passando rv e rtp em vez de vt3dp e vt3dq
  !if(level == 0) then
  !   scratch%vt3dp= 0.
  !   scratch%vt3dq= 0.
  !elseif (level > 0) then
  !   call atob(mxyzp,basic_g(ngrid)%rv,scratch%vt3dp)
  !   call atob(mxyzp,basic_g(ngrid)%rtp,scratch%vt3dq)
  !end if

  if (if_adap==0) then

     call strain(mzp,mxp,myp,ia,iz,ja,jz                       &
          ,ia_1,ja_1,iz1,jz1,jdim                                &
          ,basic_g(ngrid)%up (:,:,:) ,basic_g(ngrid)%vp (:,:,:)  &
          ,basic_g(ngrid)%wp (:,:,:) ,scratch%vt3da     (:)      &
          ,scratch%vt3db     (:)     ,scratch%vt3dc     (:)      &
          ,scratch%vt3dd     (:)     ,scratch%vt3de     (:)      &
          ,scratch%vt3df     (:)     ,scratch%vt3dg     (:)      &
          ,scratch%vt3dh     (:)     ,scratch%vt3di     (:)      &
          ,scratch%vt3dn     (:)     ,scr2      (:)      &
          ,idiffk(ngrid))

  else

     call strain_adap(mzp,mxp,myp,ia,iz,ja,jz                  &
          ,ia_1,ja_1,iz1,jz1,jdim                                &
          ,grid_g(ngrid)%lpu (:,:)   ,grid_g(ngrid)%lpv (:,:)    &
          ,grid_g(ngrid)%lpw (:,:)   ,basic_g(ngrid)%up (:,:,:)  &
          ,basic_g(ngrid)%vp (:,:,:) ,basic_g(ngrid)%wp (:,:,:)  &
          ,scratch%vt3da     (:)     ,scratch%vt3db     (:)      &
          ,scratch%vt3dc     (:)     ,scratch%vt3dd     (:)      &
          ,scratch%vt3de     (:)     ,scratch%vt3df     (:)      &
          ,scratch%vt3dg     (:)     ,scratch%vt3dh     (:)      &
          ,scratch%vt3di     (:)     ,scratch%vt3dn     (:)      &
          ,scr2      (:)     ,idiffk(ngrid)              &
          ,grid_g(ngrid)%dxm (:,:)   ,grid_g(ngrid)%dxt (:,:)    &
          ,grid_g(ngrid)%dxu (:,:)   ,grid_g(ngrid)%dxv (:,:)    &
          ,grid_g(ngrid)%dym (:,:)   ,grid_g(ngrid)%dyt (:,:)    &
          ,grid_g(ngrid)%dyu (:,:)   ,grid_g(ngrid)%dyv (:,:)    &
          ,dzm,dzt)

  endif

  if (level<=1) &
       scratch%vt3dp = 0.

  if (level>=2) &
       call ae1_l(int(mxyzp,i8), scratch%vt3dp(:), micro_g(ngrid)%rcp(:,:,:))

  !-srf 29/12/2008 adapted from OLAM
  call bruvais_OLAM(mzp, mxp, myp, ia, iz, ja, jz,                  &
       basic_g(ngrid)%theta(1,1,1), basic_g(ngrid)%rtp(1,1,1), &
       basic_g(ngrid)%rv(1,1,1),    scratch%vt3dp(1),          &
       basic_g(ngrid)%pp(1,1,1),    basic_g(ngrid)%pi0(1,1,1), &
       scratch%vt3dj(1),            grid_g(ngrid)%rtgt(1,1),   &
       grid_g(ngrid)%lpw(1,1))

  !ml/srf- for new turn scheme
  if (idiffk(ngrid) <= 3 .or. idiffk(ngrid) == 7.or. idiffk(ngrid) == 8) then
     call mxdefm(mzp,mxp,myp,ia,iz,ja,jz,ibcon,jdim            &
          ,scratch%vt3dh      (:)     ,scratch%vt3di      (:)    &
          ,scratch%vt3dj      (:)     ,scratch%vt3dk      (:)    &
          ,scratch%scr1       (:)     ,scr2       (:)    &
          ,basic_g(ngrid)%dn0 (:,:,:) ,grid_g(ngrid)%rtgt (:,:)  &
          ,grid_g(ngrid)%dxt  (:,:)   ,grid_g(ngrid)%dyt  (:,:)  &
          ,grid_g(ngrid)%lpw  (:,:)   ,mynum  )

     if (CCATT==1 .and. chemistry >= 0) then
        !srf------
        !coef de difusao horizontal diferente  para tracers
        call mxdefm_tracer(mzp,mxp,myp,ia,iz,ja,jz  &
             ,ibcon,jdim,scratch%vt3dh(:),scr3(:) &
             ,basic_g(ngrid)%dn0(:,:,:),grid_g(ngrid)%dxt(:,:),&
             grid_g(ngrid)%dyt(:,:),grid_g(ngrid)%lpw(:,:),mynum)

        !srf------
     endif

  endif

  !ML -> Nananishi and Niino (2004) scheme based on Mellor-Yamada Level 2.5
  if (idiffk(ngrid) == 7) then

!    if(IMASSFLX==1)then
      call nakanishi(mzp, mxp, myp, npatch, ia, iz, ja, jz, jdim                           &
             ,turb_g(ngrid)%tkep       ,tend%tket                ,scratch%vt3dd            &
             ,scratch%vt3de            ,scratch%vt3dh            ,scratch%vt3di            &
             ,scratch%vt3dj            ,scratch%scr1             ,grid_g(ngrid)%rtgt       &
!srf         ,basic_g(ngrid)%theta     ,scratch%vt3dp            ,scratch%vt3dq            &
             ,basic_g(ngrid)%theta     ,basic_g(ngrid)%rv        ,basic_g(ngrid)%rtp     &
!
             ,basic_g(ngrid)%dn0       ,basic_g(ngrid)%up        ,basic_g(ngrid)%vp        &
             ,leaf_g(ngrid)%veg_rough  ,leaf_g(ngrid)%patch_rough,leaf_g(ngrid)%tstar      &
             ,leaf_g(ngrid)%ustar      ,leaf_g(ngrid)%patch_area ,turb_g(ngrid)%sflux_u    &
             ,turb_g(ngrid)%sflux_v    ,turb_g(ngrid)%sflux_t    ,grid_g(ngrid)%lpu        &
             ,grid_g(ngrid)%lpv        ,grid_g(ngrid)%lpw         ,turb_g(ngrid)%kpbl      &
             ,stilt_g(ngrid)%pblhgt    ,stilt_g(ngrid)%lmo       ,stilt_g(ngrid)%ltscale   &
             ,stilt_g(ngrid)%sigw                                                          )


     !- from Marcos Longo
     !---------------------------------------------------------------------------------------!
     !    Integrating average turbulence parameters for mass flux.                           !
     !---------------------------------------------------------------------------------------!
        if(IMASSFLX==1)  call prepare_timeavg_driver(mzp,mxp,myp,ia,iz,ja,jz,dtlt,ngrid,idiffk(ngrid))
     !---------------------------------------------------------------------------------------!
!    else
!      call nakanishi_light(mzp, mxp, myp, npatch, ia, iz, ja, jz, jdim                     &
!             ,turb_g(ngrid)%tkep       ,tend%tket                ,scratch%vt3dd            &
!             ,scratch%vt3de            ,scratch%vt3dh            ,scratch%vt3di            &
!             ,scratch%vt3dj            ,scratch%scr1             ,grid_g(ngrid)%rtgt       &
!srf-opt     ,basic_g(ngrid)%theta     ,scratch%vt3dp            ,scratch%vt3dq            &
!             ,basic_g(ngrid)%theta     ,basic_g(ngrid)%rv        ,basic_g(ngrid)%rtp     &
!
!             ,basic_g(ngrid)%dn0       ,basic_g(ngrid)%up        ,basic_g(ngrid)%vp        &
!             ,leaf_g(ngrid)%veg_rough  ,leaf_g(ngrid)%patch_rough,leaf_g(ngrid)%tstar      &
!             ,leaf_g(ngrid)%ustar      ,leaf_g(ngrid)%patch_area ,turb_g(ngrid)%sflux_u    &
!             ,turb_g(ngrid)%sflux_v    ,turb_g(ngrid)%sflux_t    ,grid_g(ngrid)%lpu        &
!             ,grid_g(ngrid)%lpv        ,grid_g(ngrid)%lpw                                  &
!	     !,turb_g(ngrid)%kpbl      & !srf
!             !,stilt_g(ngrid)%pblhgt   & !srf
!	     ,stilt_g(ngrid)%lmo      &
!	     !,stilt_g(ngrid)%ltscale  & !srf
!             !,stilt_g(ngrid)%sigw     & !srf
!	     )
!    endif
  endif
  !ML

  if (idiffk(ngrid)==1) then
     call tkemy(mzp,mxp,myp,ia,iz,ja,jz,ibcon,jdim,nodei0(mynum,ngrid),nodej0(mynum,ngrid)  &
          ,turb_g(ngrid)%tkep   (:,:,:) ,tend%tket            (:)      &
          ,scratch%vt3dh        (:)     ,scratch%vt3di        (:)      &
          ,scratch%vt3dj        (:)     ,scratch%scr1         (:)      &
          ,grid_g(ngrid)%rtgt   (:,:)   ,basic_g(ngrid)%theta (:,:,:)  &
          ,basic_g(ngrid)%dn0   (:,:,:) ,basic_g(ngrid)%up    (:,:,:)  &
          ,basic_g(ngrid)%vp    (:,:,:) ,basic_g(ngrid)%wp    (:,:,:)  &
          ,turb_g(ngrid)%sflux_u(:,:)   ,turb_g(ngrid)%sflux_v(:,:)    &
          ,turb_g(ngrid)%sflux_w(:,:)   ,turb_g(ngrid)%sflux_t(:,:),vctr34 &
          ,grid_g(ngrid)%lpw    (:,:)   ,grid_g(ngrid)%lpu    (:,:)   &
          ,grid_g(ngrid)%lpv    (:,:))
  endif

  if (idiffk(ngrid)==4) then
     call mxtked(mzp,mxp,myp,ia,iz,ja,jz  &
          ,ibcon,jdim  &
          ,turb_g(ngrid)%tkep   (:,:,:) ,tend%tket            (:)      &
          ,basic_g(ngrid)%up    (:,:,:) ,basic_g(ngrid)%vp    (:,:,:)  &
          ,basic_g(ngrid)%wp    (:,:,:) ,basic_g(ngrid)%rtp   (:,:,:)  &
          ,basic_g(ngrid)%rv    (:,:,:) ,basic_g(ngrid)%theta (:,:,:)  &
          ,scratch%vt3da        (:)     ,scratch%vt3dc        (:)      &
          ,scratch%vt3dh        (:)     ,scratch%vt3dj        (:)      &
          ,scratch%scr1         (:)     ,scr2         (:)      &
          ,turb_g(ngrid)%sflux_u(:,:)   ,turb_g(ngrid)%sflux_v(:,:)    &
          ,turb_g(ngrid)%sflux_w(:,:)   ,turb_g(ngrid)%sflux_t(:,:)    &
          ,grid_g(ngrid)%dxt    (:,:)   ,grid_g(ngrid)%rtgt   (:,:)    &
          ,grid_g(ngrid)%lpw    (:,:)   )
  endif

  !_STC............................................................
  !_STC Call to subroutine tkescl for E-l closure
  !_STC (S. Trini Castelli)
  !_STC............................................................
  if (idiffk(ngrid)==5) then
     call tkescl(mzp,mxp,myp,npatch,ia,iz,ja,jz  &
          ,turb_g(ngrid)%tkep(:,:,:),tend%tket(:)  &
          ,turb_g(ngrid)%epsp(:,:,:),tend%epst(:)  &
          ,scratch%vt3da(:),scratch%vt3dc(:)  &
          ,scratch%vt3dh(:),scratch%vt3di(:)  &
          ,scratch%vt3dj(:),scratch%scr1(:)  &
          ,scr2(:) ,grid_g(ngrid)%rtgt(:,:)  &
          ,scratch%vt3dd(:),scratch%vt3de(:),grid_g(ngrid)%dxt(:,:)  &
          ,leaf_g(ngrid)%ustar(:,:,:),leaf_g(ngrid)%patch_area(:,:,:) &
          ,grid_g(ngrid)%lpw(:,:),basic_g(ngrid)%dn0(:,:,:)  )
  endif
  !_STC............................................................
  !_STC Call to subroutine tkeeps for E-eps closure
  !_STC (S. Trini Castelli)
  !_STC............................................................
  if (idiffk(ngrid)==6) then
     call tkeeps(mzp,mxp,myp,npatch,ia,iz,ja,jz  &
          ,turb_g(ngrid)%tkep(:,:,:),tend%tket(:)  &
          ,turb_g(ngrid)%epsp(:,:,:),tend%epst(:)  &
          ,scratch%vt3da(:),scratch%vt3dc(:)  &
          ,scratch%vt3dh(:),scratch%vt3di(:)  &
          ,scratch%vt3dj(:),scratch%scr1(:)  &
          ,scr2(:) ,grid_g(ngrid)%rtgt(:,:)  &
          ,leaf_g(ngrid)%ustar(:,:,:),leaf_g(ngrid)%patch_area(:,:,:) &
          ,grid_g(ngrid)%lpw(:,:),basic_g(ngrid)%dn0(:,:,:)  )
  endif
  !_STC..................................................
  !_STC    Note: from subroutines TKESCL, TKEEPS :
  !_STC           VT3DI=Ke
  !_STC           SCR1=Km
  !_STC           VT3DH = Kh
  !_STC           SCR2 = SCR1 = Km
  !_STC..................................................
  !_STC............................................................

  if (idiffk(ngrid) == 8) then
  !-----------------------------------
  !Call subrotines HCV_driver (Campos Velho - KZZ)
  !-----------------------------------
      call HCV_driver(mzp,mxp,myp,ia,iz,ja,jz,mynum   &
             ,basic_g(ngrid)%up    (:,:,:)        &
	     ,basic_g(ngrid)%vp    (:,:,:)        &
             ,basic_g(ngrid)%theta (:,:,:)        &
	     ,basic_g(ngrid)%rv    (:,:,:)	  &
             ,basic_g(ngrid)%dn0   (:,:,:)        &
	     ,turb_g(ngrid)%sflux_t(:,:)          &
	     ,turb_g(ngrid)%sflux_r(:,:)          &
	     ,turb_g(ngrid)%sflux_u(:,:)          &
	     ,turb_g(ngrid)%sflux_v(:,:)          &
	     ,grid_g(ngrid)%rtgt   (:,:)          &
	     ,grid_g(ngrid)%lpw    (:,:)          &
	     ,scratch%scr1	   (:)            &
	     ,scratch%vt3dh	   (:)            &
             ,scratch%vt3di        (:)            &
             ,scratch%vt3dj        (:)            &
	     ,scratch%vt3dk        (:)            )
   endif

  call klbnd(mzp,mxp,myp,ibcon,jdim  &
       ,scratch%scr1 (:),basic_g(ngrid)%dn0(:,:,:),grid_g(ngrid)%lpw(:,:))
  call klbnd(mzp,mxp,myp,ibcon,jdim  &
       ,scr2 (:),basic_g(ngrid)%dn0(:,:,:),grid_g(ngrid)%lpw(:,:))
  call klbnd(mzp,mxp,myp,ibcon,jdim  &
       ,scratch%vt3dh(:),basic_g(ngrid)%dn0(:,:,:),grid_g(ngrid)%lpw(:,:))


  ! CATT
  if (CCATT==1 .and. chemistry >= 0) then
     !srf----
     call klbnd(mzp,mxp,myp,ibcon,jdim,scr3(:) &
          ,basic_g(ngrid)%dn0(:,:,:),grid_g(ngrid)%lpw(:,:))
  endif

  !_STC ....... boundary conditions even on Ke diffusion coefficient
  if (idiffk(ngrid)==5 .or. idiffk(ngrid)==6) &
       call klbnd(mzp,mxp,myp,ibcon,jdim  &
       ,scratch%vt3di(:),basic_g(ngrid)%dn0(:,:,:),grid_g(ngrid)%lpw(:,:))

  !bob  swap new hkm, vkm, and vkh with past time level:  lagged K's have
  !bob  internal lateral boundary values from neighboring nodes

  ind = 0
  do j = 1,nodemyp(mynum,ngrid)
     do i = 1,nodemxp(mynum,ngrid)
        do k = 1,nnzp(ngrid)
           ind = ind + 1
           s1 = scr2(ind)
           s2 = scratch%scr1(ind)
           s3 = scratch%vt3dh(ind)
           scr2(ind) = turb_g(ngrid)%hkm(k,i,j)
           scratch%scr1(ind) = turb_g(ngrid)%vkm(k,i,j)
           scratch%vt3dh(ind) = turb_g(ngrid)%vkh(k,i,j)
           !! also for vt3di = K(tke) ?????    22 March 02
           turb_g(ngrid)%hkm(k,i,j) = s1
           turb_g(ngrid)%vkm(k,i,j) = s2
           turb_g(ngrid)%vkh(k,i,j) = s3
        enddo
     enddo
  enddo

  if (if_adap==0) then
!call dumpVarAllLatLonk(basic_g(ngrid)%up, 'tUP'  ,415,0,0,1,mxp,1,myp,1,mzp,0.0,0.0) !ok
!call dumpVarAllLatLonk(basic_g(ngrid)%Vp, 'tVP'  ,415,0,0,1,mxp,1,myp,1,mzp,0.0,0.0) !ok
!call dumpVarAllLatLonk(basic_g(ngrid)%Wp, 'tWP'  ,415,0,0,1,mxp,1,myp,1,mzp,0.0,0.0) !ok
!call dumpVarAllLatLonk(basic_g(ngrid)%dn0, 'DN0'  ,415,0,0,1,mxp,1,myp,1,mzp,0.0,0.0) !ok
!call dumpVarAllLatLonk(basic_g(ngrid)%dn0v, 'DN0V'  ,415,0,0,1,mxp,1,myp,1,mzp,0.0,0.0) !ok
!call dumpVarAllLatLonk(basic_g(ngrid)%dn0u, 'DN0U'  ,415,0,0,1,mxp,1,myp,1,mzp,0.0,0.0) !ok
!call dumpVarAllLatLonk(grid_g(ngrid)%rtgu   ,'rtgu',415,0,0,1,mxp,1,myp,1,1,0.0,0.0) !ok
!call dumpVarAllLatLonk(grid_g(ngrid)%rtgv   ,'rtgv',415,0,0,1,mxp,1,myp,1,1,0.0,0.0) !ok
!call dumpVarAllLatLonk(grid_g(ngrid)%rtgt   ,'rtgt',415,0,0,1,mxp,1,myp,1,1,0.0,0.0) !ok
!call dumpVarAllLatLonk(turb_g(ngrid)%sflux_u,'sflux_u',415,0,0,1,mxp,1,myp,1,1,0.0,0.0) !ok
!call dumpVarAllLatLonk(turb_g(ngrid)%sflux_v,'sflux_v',415,0,0,1,mxp,1,myp,1,1,0.0,0.0) !ok
!call dumpVarAllLatLonk(turb_g(ngrid)%sflux_w,'sflux_w',415,0,0,1,mxp,1,myp,1,1,0.0,0.0) !***BAD***
     call diffvel(mzp,mxp,myp,ia,iz,ja,jz,jdim,ia_1,ja_1             &
          ,ia1,ja1,iz_1,jz_1,iz1,jz1,izu,jzv,idiffk(ngrid)             &
          ,basic_g(ngrid)%up    (:,:,:) ,basic_g(ngrid)%vp    (:,:,:)  &
          ,basic_g(ngrid)%wp    (:,:,:) ,tend%ut              (:)      &
          ,tend%vt              (:)     ,tend%wt              (:)      &
          ,scratch%vt3da        (:)     ,scratch%vt3db        (:)      &
          ,scratch%vt3dc        (:)     ,scratch%vt3dd        (:)      &
          ,scratch%vt3de        (:)     ,scratch%vt3df        (:)      &
          ,scratch%vt3dg        (:)     ,scratch%vt3dj        (:)      &
          ,scratch%vt3dk        (:)     ,scratch%vt3dl        (:)      &
          ,scratch%vt3dm        (:)     ,scratch%vt3dn        (:)      &
          ,scratch%vt3do        (:)     ,grid_g(ngrid)%rtgu   (:,:)    &
          ,grid_g(ngrid)%rtgv   (:,:)   ,grid_g(ngrid)%rtgt   (:,:)    &
          ,turb_g(ngrid)%sflux_u(:,:)   ,turb_g(ngrid)%sflux_v(:,:)    &
          ,turb_g(ngrid)%sflux_w(:,:)   ,basic_g(ngrid)%dn0   (:,:,:)  &
          ,basic_g(ngrid)%dn0u  (:,:,:) ,basic_g(ngrid)%dn0v  (:,:,:)  &
          ,scratch%scr1         (:)     ,scr2         (:),ibcon,mynum)
!call dumpVarAllLatLonk(tend%wt, 'tWT'  ,433,0,0,1,mxp,1,myp,1,mzp,0.0,0.0)
  else

     call diffvel_adap(mzp,mxp,myp,ia,iz,ja,jz,jdim                  &
          ,iz1,jz1,izu,jzv,idiffk(ngrid)                               &
          ,basic_g(ngrid)%up    (:,:,:) ,basic_g(ngrid)%vp    (:,:,:)  &
          ,basic_g(ngrid)%wp    (:,:,:) ,tend%ut              (:)      &
          ,tend%vt              (:)     ,tend%wt              (:)      &
          ,scratch%vt3da        (:)     ,scratch%vt3db        (:)      &
          ,scratch%vt3dc        (:)     ,scratch%vt3dd        (:)      &
          ,scratch%vt3de        (:)     ,scratch%vt3df        (:)      &
          ,scratch%vt3dg        (:)     ,scratch%vt3dj        (:)      &
          ,scratch%vt3dk        (:)     ,scratch%vt3dl        (:)      &
          ,scratch%vt3dm        (:)     ,scratch%vt3dn        (:)      &
          ,scratch%vt3do        (:)     ,grid_g(ngrid)%aru    (:,:,:)  &
          ,grid_g(ngrid)%arv    (:,:,:) ,grid_g(ngrid)%arw    (:,:,:)  &
          ,grid_g(ngrid)%volu   (:,:,:) ,grid_g(ngrid)%volv   (:,:,:)  &
          ,grid_g(ngrid)%volw   (:,:,:) ,grid_g(ngrid)%lpu    (:,:)    &
          ,grid_g(ngrid)%lpv    (:,:)   ,grid_g(ngrid)%lpw    (:,:)    &
          ,turb_g(ngrid)%sflux_u(:,:)   ,turb_g(ngrid)%sflux_v(:,:)    &
          ,turb_g(ngrid)%sflux_w(:,:)   ,basic_g(ngrid)%dn0   (:,:,:)  &
          ,basic_g(ngrid)%dn0u  (:,:,:) ,basic_g(ngrid)%dn0v  (:,:,:)  &
          ,scratch%scr1         (:)     ,scr2         (:)      &
          ,grid_g(ngrid)%topma  (:,:)   ,ibcon,mynum)

  endif

  ! Convert momentum K's to scalar K's, if necessary
  !-ml/srf - for new turb scheme
  if (idiffk(ngrid) <= 3 .or. idiffk(ngrid) == 7 .or. idiffk(ngrid) == 8) then
     do ind = 1,mxyzp
        scr2(ind) = scr2(ind) * xkhkm(ngrid)
     enddo
  elseif (idiffk(ngrid) == 4) then
     do ind = 1,mxyzp
        scratch%vt3di(ind) = 2. * scratch%scr1(ind)
     enddo
  endif

  !- CATT
  if (CCATT==1.and. CHEMISTRY >=0) then

     ind = 0
     do j = 1,nodemyp(mynum,ngrid)
        do i = 1,nodemxp(mynum,ngrid)
           do k = 1,nnzp(ngrid)
              ind = ind + 1

              s1 = scr3(ind)
              scr3(ind)  = turb_s(ngrid)%hksc(k,i,j)
              turb_s(ngrid)%hksc(k,i,j) = s1

              ! salva o atual coef para o proximo passo no tempo.
           enddo
        enddo
     enddo

     if (idiffk(ngrid)<=3 .or. idiffk(ngrid) == 7.or. idiffk(ngrid) == 8) then

        ind = 0
        do j = 1,nodemyp(mynum,ngrid)
           do i = 1,nodemxp(mynum,ngrid)
              do k = 1,nnzp(ngrid)
                 ind = ind + 1
                 scr3(ind)  =  scr3(ind)  * xkhkm(ngrid)
              enddo
           enddo
        enddo

     endif
  endif

  do n = 1,num_scalar(ngrid)

     scalarp => scalar_tab(n,ngrid)%a_var_p
     scalart => scalar_tab(n,ngrid)%a_var_t


     scratch%vt2da = 0.

     if (nstbot==1) then
        if (scalar_tab(n,ngrid)%name=='THP' .or. &
            scalar_tab(n,ngrid)%name=='THC') then
           call atob(mxp*myp, turb_g(ngrid)%sflux_t(:,:), scratch%vt2da(:))

           ! Large Scale Forcing for GRELL CUPAR
           !if (nnqparm(ngrid)>=2) then
           !   !---------------------------------------
           !   !srf- Large Scale Forcing for GRELL CUPAR
           !   lsfcupar_p => cuforc_sh_g(ngrid)%lsfth
              !---------------------------------------
          ! endif
        elseif (scalar_tab(n,ngrid)%name=='RTP') then
          call atob(mxp*myp, turb_g(ngrid)%sflux_r(:,:), scratch%vt2da(:))

           ! Large Scale Forcing for GRELL CUPAR not used
           ! CATT
           !if (nnqparm(ngrid)>=2) then
           !   !---------------------------------------
           !   !srf- Large Scale Forcing for GRELL CUPAR
           !   lsfcupar_p => cuforc_sh_g(ngrid)%lsfrt
           !   !---------------------------------------
           !endif

        endif
     endif

     ! 3/10/01 - Define ksf below, the "K scalar flag", to let subroutine diffsclr
     ! know which vertical K is being passed to it.  If diffsclr sees that it's
     ! a different K from the previous one, diffsclr will re-compute the tridiff
     ! matrix coefficients.  In order to use vertical scalar K's other than
     ! vt3dh and vt3di, use ksf = 3, ksf = 4, etc. for each different K.

     !_STC..................................................
     !_STC Corrections to account for the new idiffk options
     !_STC for E-l and E-eps closure. Isotropy hypothesis.
     !_STC (S. Trini Castelli)
     !_STC..................................................

     if (scalar_tab(n,ngrid)%name=='TKEP') then
        vkh_p => scratch%vt3di
        hkh_p => scr2
        !-ml/srf - for new turb scheme
        if (idiffk(ngrid) >= 4 .and. idiffk(ngrid) /= 7 .and. idiffk(ngrid) /= 8) &
	 hkh_p => scratch%vt3di
!       if (idiffk(ngrid)>=4) hkh_p => scratch%vt3di
        ksf = 1
     elseif (scalar_tab(n,ngrid)%name=='EPSP') then
        vkh_p => scratch%vt3di
        hkh_p => scr2
        !-ml/srf - for new turb scheme
        if (idiffk(ngrid) >= 4 .and. idiffk(ngrid) /= 7.and. idiffk(ngrid) /= 8) &
	 hkh_p => scratch%vt3di
!        if (idiffk(ngrid)>=4)  hkh_p => scratch%vt3di
        ksf = 3
        ! Convert Ktke to Keps; it will be converted back after use below
        call ae1t0_l(mxyzp, vkh_p, vkh_p, (ALF_EPS/ALF_TKE))
        call ae1t0_l(mxyzp, hkh_p, hkh_p, (ALF_EPS/ALF_TKE))
     else
        vkh_p => scratch%vt3dh
        hkh_p => scr2
        !-ml/srf - for new turb scheme
        if (idiffk(ngrid)>=4 .and. idiffk(ngrid) /= 7.and. idiffk(ngrid) /= 8) &
	 hkh_p => scratch%vt3dh
!       if (idiffk(ngrid)>=4)  hkh_p => scratch%vt3di
        ksf = 2
     endif


     ! CATT
     if (CCATT==1 .and. CHEMISTRY >=0) then
        !srf----------------- Hor. Diffusion Coef for tracers
        if(n > (num_scalar(ngrid) - (NADDSC + NSPECIES_TRANSPORTED))) then
           if (idiffk(ngrid) < 4 .or.  idiffk(ngrid) == 7.or.  idiffk(ngrid) == 8) then

              hkh_p => scr3
           endif
        endif
        !srf--------------------
     endif


     if (if_adap==0) then

        ! NEW WAY
        call diffsclr(mzp, mxp, myp, ia, iz, ja, jz, jdim,         &
             ia_1, ja_1, ia1, ja1, iz_1, jz_1, iz1, jz1, n, ksf,   &
             scalarp(:), scalart(:),                               &
             scratch%vt3da(:), scratch%vt3db(:), scratch%vt3df(:), &
             scratch%vt3dg(:), scratch%vt3dj(:), scratch%vt3dk(:), &
             scratch%vt3do(:), scratch%vt3dc(:), scratch%vt3dd(:), &
             scratch%vt3dl(:), scratch%vt3dm(:), scratch%vt2db(:), &
             grid_g(ngrid)%rtgt(:,:), scratch%vt2da(:),            &
             basic_g(ngrid)%dn0(:,:,:),                            &
             vkh_p(:)                 , hkh_p(:)                   )
        !
        IF (nnqparm(ngrid)>=2) then
          ! SGScale Forcing for GRELL CUPAR
	  if (scalar_tab(n,ngrid)%name=='THP' .or. scalar_tab(n,ngrid)%name=='THC')     &
	       call PBLforcing(ngrid, mzp, mxp, myp, ia, iz, ja, jz, &
                             scratch%vt3df, scalarp(:), cuforc_sh_g(ngrid)%lsfth, n)

          if (scalar_tab(n,ngrid)%name=='RTP')                                          &
	       call PBLforcing(ngrid, mzp, mxp, myp, ia, iz, ja, jz, &
                             scratch%vt3df, scalarp(:), cuforc_sh_g(ngrid)%lsfrt, n)

       ENDIF
     else

        call diffsclr_adap(mzp,mxp,myp,ia,iz,ja,jz,jdim,n,ksf      &
             ,grid_g(ngrid)%lpw(:,:)     ,scalarp                    &
             ,scalart                    ,scratch%vt3da     (:)      &
             ,scratch%vt3dc      (:)     ,scratch%vt3df     (:)      &
             ,scratch%vt3dg      (:)     ,scratch%vt3dj     (:)      &
             ,scratch%vt3dk      (:)     ,scratch%vt3dl     (:)      &
             ,scratch%vt3dm      (:)     ,scratch%vt3do     (:)      &
             ,scratch%vt2da      (:)     ,scratch%vt2db     (:)      &
             ,basic_g(ngrid)%dn0 (:,:,:) ,vkh_p                      &
             ,hkh_p                      ,grid_g(ngrid)%aru (:,:,:)  &
             ,grid_g(ngrid)%arv  (:,:,:) ,grid_g(ngrid)%arw (:,:,:)  &
             ,grid_g(ngrid)%volt (:,:,:) ,scratch%vt3db     (:)      &
             ,grid_g(ngrid)%dxu  (:,:)   ,grid_g(ngrid)%dyv (:,:)    &
             ,grid_g(ngrid)%topma(:,:)                               )
     endif

     if (scalar_tab(n,ngrid)%name == 'EPSP') then
        call ae1t0_l(mxyzp, vkh_p, vkh_p, (ALF_TKE/ALF_EPS))
        call ae1t0_l(mxyzp, hkh_p, hkh_p, (ALF_TKE/ALF_EPS))
     endif

  enddo

end subroutine diffuse


!     ******************************************************************

subroutine strain(m1,m2,m3,ia,iz,ja,jz,ia_1,ja_1,iz1,jz1  &
     ,jd,up,vp,wp,vt3da,vt3db,vt3dc,vt3dd,vt3de  &
     ,vt3df,vt3dg,vt3dh,vt3di,vt3dn,scr2,idiffk)

  implicit none

  integer, intent(in)  :: m1   &
                        , m2   &
                        , m3   &
                        , ia   &
                        , iz   &
                        , ja   &
                        , jz   &
                        , ia_1 &
                        , ja_1 &
                        , iz1  &
                        , jz1  &
                        , jd   &
                        , idiffk

  real, dimension(m1,m2,m3), intent(in) :: up,vp,wp

  real, dimension(m1,m2,m3), intent(inout):: vt3da,      &
                                             vt3db,      &
                                             vt3dc,      &
                                             vt3dd,      &
                                             vt3de,      &
                                             vt3df,      &
                                             vt3dg,      &
                                             vt3dh,      &
                                             vt3di,      &
                                             vt3dn,      &
                                             scr2

  !local variables:
  integer              :: i,j,k


  call grad(m1, m2, m3, ia  , iz1, ja  , jz , up, vt3da, 'XDIR', 'UPNT')
  call grad(m1, m2, m3, ia_1, iz , ja_1, jz , vp, vt3db, 'XDIR', 'VPNT')
  call grad(m1, m2, m3, ia_1, iz , ja  , jz , wp, vt3df, 'XDIR', 'WPNT')

  call grad(m1, m2, m3, ia_1, iz , ja_1, jz , up, vt3dn, 'YDIR', 'UPNT')
  call grad(m1, m2, m3, ia  , iz , ja  , jz1, vp, vt3dc, 'YDIR', 'VPNT')
  call grad(m1, m2, m3, ia  , iz , ja_1, jz , wp, vt3dg, 'YDIR', 'WPNT')

  call grad(m1, m2, m3, ia_1, iz , ja  , jz , up, vt3dd, 'ZDIR', 'UPNT')
  call grad(m1, m2, m3, ia  , iz , ja_1, jz , vp, vt3de, 'ZDIR', 'VPNT')
  if(idiffk.ge.3 .and. idiffk /= 7)then
     call grad(m1,m2,m3,ia,iz,ja,jz,wp,scr2,'ZDIR','WPNT')
  endif

  if (idiffk .le. 2 .or. idiffk == 7 .or. idiffk == 8) then
     do j = ja,jz
        do i = ia,iz
           do k = 2,m1-1
              vt3dh(k,i,j) =2. * (vt3da(k,i,j) * vt3da(k,i,j)  &
                   + vt3dc(k,i,j) * vt3dc(k,i,j))  &
                   + .0625 * (vt3db(k,i,j) + vt3db(k,i-1,j)  &
                   + vt3db(k,i,j-jd) + vt3db(k,i-1,j-jd)  &
                   + vt3dn(k,i,j) + vt3dn(k,i-1,j)  &
                   + vt3dn(k,i,j-jd) + vt3dn(k,i-1,j-jd)) ** 2
              vt3di(k,i,j) = .0625 * ((vt3dd(k,i,j) + vt3dd(k-1,i,j)  &
                   + vt3dd(k,i-1,j) + vt3dd(k-1,i-1,j)) ** 2  &
                   + (vt3de(k,i,j) + vt3de(k-1,i,j)  &
                   + vt3de(k,i,j-jd) + vt3de(k-1,i,j-jd)) ** 2)
           enddo
        enddo
     enddo
  else
     do j = ja,jz
        do i = ia,iz
           do k = 2,m1-1
              vt3da(k,i,j) = 2. * vt3da(k,i,j)
              vt3dc(k,i,j) = 2. * vt3dc(k,i,j)
              scr2(k,i,j) = 2. * scr2(k,i,j)
              vt3db(k,i,j) = vt3db(k,i,j) + vt3dn(k,i,j)
              vt3dn(k,i,j) = vt3db(k,i,j)
              vt3dd(k,i,j) = vt3dd(k,i,j) + vt3df(k,i,j)
              vt3de(k,i,j) = vt3de(k,i,j) + vt3dg(k,i,j)
              vt3di(k,i,j) = 0.333333  &
                   * (vt3da(k,i,j) + vt3dc(k,i,j) + scr2(k,i,j))
           enddo
        enddo

        do k = 2,m1-1
           vt3da(k,iz1,j) = 2. * vt3da(k,iz1,j)
           vt3db(k,ia_1,j) = vt3db(k,ia_1,j) + vt3dn(k,ia_1,j)
           vt3dn(k,ia_1,j) = vt3db(k,ia_1,j)
           vt3dd(k,ia_1,j) = vt3dd(k,ia_1,j) + vt3df(k,ia_1,j)
        enddo
     enddo

     do i = ia_1,iz
        do k = 2,m1-1
           vt3dc(k,i,jz1) = 2. * vt3dc(k,i,jz1)
           vt3db(k,i,ja_1) = vt3db(k,i,ja_1) + vt3dn(k,i,ja_1)
           vt3dn(k,i,ja_1) = vt3db(k,i,ja_1)
           vt3de(k,i,ja_1) = vt3de(k,i,ja_1) + vt3dg(k,i,ja_1)
        enddo
     enddo

     do j = ja,jz
        do i = ia,iz
           do k = 2,m1-1
              vt3dh(k,i,j) = .5 * (  &
                   (vt3da(k,i,j) - vt3di(k,i,j)) ** 2  &
                   + (vt3dc(k,i,j) - vt3di(k,i,j)) ** 2  &
                   + ( scr2(k,i,j) - vt3di(k,i,j)) ** 2)  &
                   + .0625 * ((vt3db(k,i,j) + vt3db(k,i-1,j)  &
                   + vt3db(k,i,j-jd) + vt3db(k,i-1,j-jd)) ** 2  &
                   + (vt3dd(k,i,j) + vt3dd(k,i-1,j)  &
                   + vt3dd(k-1,i,j) + vt3dd(k-1,i-1,j)) ** 2  &
                   + (vt3de(k,i,j) + vt3de(k-1,i,j)  &
                   + vt3de(k,i,j-jd) + vt3de(k-1,i,j-jd)) ** 2)
              vt3di(k,i,j) = vt3dh(k,i,j)
           enddo
        enddo
     enddo
  endif

  return
end subroutine strain

!     *****************************************************************

subroutine bruvais(m1, m2, m3, ia, iz, ja, jz, theta, rtp, rv, rcp, &
     pp, pi0, en2, rtgt, lpw_R)

  use mem_scratch, only: &
       vctr11,    &     !INTENT(INOUT)
       vctr12,    &     !INTENT(INOUT)
       vctr32,    &     !INTENT(INOUT)
       vctr1,     &     !INTENT(INOUT)
       vctr2,     &     !INTENT(INOUT)
       vctr3,     &     !INTENT(INOUT)
       vctr4            !INTENT(INOUT)

  use micphys, only: &
       level,     &     !INTENT(in)
       ipris,     &     !INTENT(in)
       isnow,     &     !INTENT(in)
       igraup,    &     !INTENT(in)
       iaggr,     &     !INTENT(in)
       ihail            !INTENT(in)

  use mem_grid, only: &
       zt,        &     !INTENT(in)
       nzp,       &     !INTENT(in)
       nz,        &     !INTENT(in)
       nzpmax           !INTENT(in)

  use rconstants, only: &
       alvl,      &     !INTENT(in)
       rgas,      &     !INTENT(in)
       ep,        &     !INTENT(in)
       cp,        &     !INTENT(in)
       alvi,      &     !INTENT(in)
       p00,       &     !INTENT(in)
       cpor,      &     !INTENT(in)
       g                !INTENT(in)

  implicit none
  ! Arguments:
  integer, INTENT(IN) :: m1, m2, m3, ia, iz, ja, jz
  real, INTENT(IN)    :: theta(m1,m2,m3), pp(m1,m2,m3), pi0(m1,m2,m3), &
       rtp(m1,m2,m3), rcp(m1,m2,m3), rv(m1,m2,m3)
  real, INTENT(INOUT) :: en2(m1,m2,m3)
  real, INTENT(IN)    :: rtgt(m2,m3)
  real, INTENT(IN) :: lpw_R(m2,m3)

  !local variables
  integer :: lpw(m2,m3)
  integer :: i, j, k, iweten, iwdiffk, ki, k2, k1
  real :: c1, c2, c3, ci1, ci2, ci3, rvlsi, rvii
  real :: pi(nzpmax), temp(nzpmax), prt(nzpmax), rvls(nzpmax), rc(nzpmax)
  ! **(JP)** fatora expressoes logicas para fora dos lacos
  logical :: log1, log2, log3, log4

  !lfr: Solving a  problem with integer inside vtables
  lpw=int(lpw_R)

  !     calculate brunt-vaisalla frequency squared (en2)

  iweten = 1

  iwdiffk = 0
  c1  = alvl/rgas
  c2  = ep*alvl**2/(cp*rgas)
  c3  = alvl/cp
  ci1 = alvi/rgas
  ci2 = ep*alvi**2/(cp*rgas)
  ci3 = alvi/cp

  ! **(JP)** fatora expressoes logicas para fora dos lacos
  log1 = level>=1
  log2 = level>=2 .and. iweten==1
  log3 = (ipris>=1 .or. isnow>=1 .or. igraup>=1 .or. iaggr>=1 .or. ihail>=1) &
       .and. level==3
  log4 = level==3
  ! **(JP)** fim modificacao

  !     calculate potential temperature profile

  do j=ja,jz
     do i=ia,iz

        k2 = lpw(i,j)
        k1 = k2 - 1

        do k=k1,m1
           vctr11(k) = theta(k,i,j)
           vctr12(k) = theta(k,i,j)
           vctr32(k) = 0.
        enddo
        ! **(JP)** fatora expressoes logicas para fora dos lacos
        if (log1) then
        ! **(JP)** fim modificacao
           do k=k1,m1
              vctr12(k) = vctr11(k)*(1. + 0.61*rv(k,i,j))
              vctr32(k) = (rtp(k,i,j) - rv(k,i,j))
           enddo
        endif

        !     check for saturation if level is 2 or greater.
        ! **(JP)** fatora expressoes logicas para fora dos lacos
        if (log2) then
        ! **(JP)** fim modificacao
           do k=k1,m1
              pi(k)   = (pp(k,i,j) + pi0(k,i,j))/cp
              temp(k) = theta(k,i,j)*pi(k)
              prt(k)  = p00*pi(k)**cpor
           enddo
           call mrsl(m1,prt(:), temp(:), rvls(:))
           do k=k2-1,m1
              vctr2(k) = c1
              vctr3(k) = c2
              vctr4(k) = c3
           enddo
           ki = m1 + 1

           !     if any ice phase microphysics are activated ....
           ! **(JP)** fatora expressoes logicas para fora dos lacos
           if (log3) then
           ! **(JP)** fim modificacao
              !  find level of -20 c.  assume ice saturation above this
              !   level.
              do k=k1,m1
                 if (temp(k)<=253.16) then
                    ki = k
                    go to 10
                 endif
              enddo
              ki = m1 + 1
10            continue
              call mrsi(m1-ki+1, prt(ki), temp(ki), rvls(ki))
              !-srf 19/03/2005
              !bug: Subscript out of range for array prt
              !    subscript=0, lower bound=1, upper bound=132, dimension=1
              ! tempc < 253 sobre parte da Antartica.
              !              ki = max(ki,k2)                        !solucao 1
              !              call mrsi(1,prt(ki-1),temp(ki-1),rvlsi)
              if (ki>1) call mrsi(1, prt(ki-1), temp(ki-1), rvlsi) ! solucao 2
              !srf

              do k=ki,m1
                 vctr2(k) = ci1
                 vctr3(k) = ci2
                 vctr4(k) = ci3
              enddo
           endif

           ! **(JP)** fatora expressoes logicas para fora dos lacos
           if (log4) then
           ! **(JP)** fim modificacao
              do k=k1,m1
                 rc(k) = rcp(k,i,j)
              enddo
           else
              do k=k1,m1
                 rc(k) = max(rv(k,i,j)/rvls(k) - 0.999, 0.)
              enddo
           endif

        endif

        do k=k2,m1-1
           vctr1(k) = g/((zt(k+1) - zt(k-1))*rtgt(i,j))
        enddo

        ! **(JP)** fatora expressoes logicas para fora dos lacos
        if (log2) then
        ! **(JP)** fim modificacao
           do k=k2,m1-1
              if (rc(k)>0.) then
                 rvii = rvls(k-1)
                 if (k==ki) rvii = rvlsi
                 en2(k,i,j) = vctr1(k)*                        &
                      ((1. + vctr2(k)*rvls(k)/temp(k))/        &
                      (1. + vctr3(k)*rvls(k)/temp(k)**2)*      &
                      ((vctr11(k+1) - vctr11(k-1))/vctr11(k) + &
                      vctr4(k)/temp(k)*(rvls(k+1) - rvii)) -   &
                      (rtp(k+1,i,j) - rtp(k-1,i,j)))
              else
                 en2(k,i,j) = vctr1(k)*((vctr12(k+1)-vctr12(k-1))/vctr12(k) - &
                      (vctr32(k+1) - vctr32(k-1)))
              endif

           enddo
        else
           do k=k2,m1-1
              en2(k,i,j) = vctr1(k)*((vctr12(k+1)-vctr12(k-1))/vctr12(k) - &
                   (vctr32(k+1) - vctr32(k-1)))
           enddo
        endif
        ! **(JP)** remove para fora do laco, permitindo vetorizacao
        !en2(k1,i,j) = en2(k2,i,j)
        !en2(nzp,i,j)=en2(nz,i,j)
        ! **(JP)** fim da modificacao

     enddo
  enddo

  ! **(JP)** removido de dentro do laco
  do j=ja,jz
     do i=ia,iz
        en2(lpw(i,j)-1,i,j) = en2(lpw(i,j),i,j)
     enddo
  enddo
  do j=ja,jz
     do i=ia,iz
        en2(nzp,i,j) = en2(nz,i,j)
     enddo
  enddo
  ! **(JP)** fim da modificacao

end subroutine bruvais

!     *****************************************************************

subroutine mxdefm(m1, m2, m3, ia, iz, ja, jz, ibcon, jd,  &
     vt3dh, vt3di, vt3dj, vt3dk, scr1, scr2, dn0, rtgt, dxt, dyt, lpw_R, mynum)

  !     +-------------------------------------------------------------+
  !     \   this routine calculates the mixing coefficients with a    \
  !     \     smagorinsky-type deformational based k with an optional \
  !     \     unstable brunt-vaisala enhancement and an optional      \
  !     \     richardson number modification.                         \
  !     +-------------------------------------------------------------+

  use mem_scratch, only: &
       vctr1,            &     !INTENT(INOUT)
       vctr2                   !INTENT(INOUT)

  use mem_grid, only:    &
       ngrid,            &     !INTENT(in)
       zm,               &     !INTENT(in)
       zt,               &     !INTENT(in)
       akminvar                !INTENT(in) ! For new turb scheme

  use mem_turb, only:    &
       csx,              &     !INTENT(in)
       idiffk,           &     !INTENT(in)
       rmin,             &     !INTENT(out)
       rmax,             &     !INTENT(out)
       zkhkm,            &     !INTENT(in)
       csz,              &     !INTENT(in)
       akmin                   !INTENT(in)

  use rconstants, only:  &
       vonk

  implicit none
  ! Arguments:
  integer, intent(IN) :: m1, m2, m3, ia, iz, ja, jz
  integer, intent(IN) :: ibcon, jd, mynum   !- EHE -> nao sao usadas!!!
  real, INTENT(INOUT) :: vt3dh(m1,m2,m3), vt3dk(m1,m2,m3), &
       scr1(m1,m2,m3), scr2(m1,m2,m3)
  real, INTENT(IN)    :: dn0(m1,m2,m3), vt3di(m1,m2,m3), vt3dj(m1,m2,m3)
  real, INTENT(IN)    :: rtgt(m2,m3), dxt(m2,m3)
  real, INTENT(IN)    :: dyt(m2,m3)  !- EHE -> nao e' usada!!!
  real, INTENT(IN) :: lpw_R(m2,m3)
  ! Local variables:
  integer :: lpw(m2,m3)
  integer :: i, j, k, irich, ienfl
  real :: csx2, sq300, enfl, rchmax, c1, c2, c3, c4, akm, ambda, vkz2

  lpw=int(lpw_R)

  irich = 1
  ienfl = 1

  csx2 = csx(ngrid)*csx(ngrid)
  sq300 = 90000.
  if (idiffk(ngrid)==2 .or. idiffk(ngrid)==3) then
     rmin = -100.
     rmax = 1./zkhkm(ngrid)
     do j=ja,jz
        do i=ia,iz
           do k=lpw(i,j),m1-1
              vt3dk(k,i,j) = max(min(vt3dj(k,i,j)/&
                   max(vt3di(k,i,j), 1.e-15), rmax), rmin)
           enddo
        enddo
     enddo
     enfl   = float(ienfl)
     rchmax = 1.0 + 9.0*float(irich)
     !--(DMK-CCATT)---------------------------------------------------------
     !mudanca rams60 nova versao
     do k = 2, m1
        !--(DMK-original)------------------------------------------------------
        !     do k = lpw(i,j),m1
        !--(DMK-CCATT)---------------------------------------------------------
        vctr1(k) = csz(ngrid)*(zm(k) - zm(k-1))
        vctr2(k) = vctr1(k)*vctr1(k)
     enddo
  endif

  if (idiffk(ngrid)==1 .or. idiffk(ngrid)==7.or. idiffk(ngrid)==8) then
     do j=ja,jz
        do i=ia,iz
           c2 = 1.0/(dxt(i,j)*dxt(i,j))
           c3 = csx2*c2
           akm = abs(akmin(ngrid))*0.075*c2**(0.666667)
           if (AKMIN(ngrid)<0.) then !ml/srf - for new turb scheme
              !----
              !srf-define diferentes AKMINs para melhorar estabilidade
              !srf-sobre os Andes
              !!akm=extra2d(5,ngrid)%d2(i,j)*akm
              akm = akminvar(ngrid)%akmin2d(i,j)*akm
           endif
           do k=lpw(i,j),m1-1
              scr2(k,i,j) = dn0(k,i,j)*max(akm, c3*sqrt(vt3dh(k,i,j)))
           enddo
        enddo
     enddo
  elseif (idiffk(ngrid)==2) then
     do j=ja,jz
        do i=ia,iz
           c1  = rtgt(i,j)*rtgt(i,j)
           c2  = 1.0/(dxt(i,j)*dxt(i,j))
           c3  = csx2*c2
           akm = abs(akmin(ngrid))*0.075*c2**(0.666667)
           c4  = vonk*vonk*c1
           do k=lpw(i,j),m1-1
              ! old csz*dz len  scr1(k,i,j) = dn0(k,i,j)*c1*vctr2(k)

              ! asymptotic vertical scale length from bjorn with modifications:
              ! c3 is (csx*dx)^2, c1*vctr2(k) is (csz*dz)^2, sq300 is the square
              ! of 300 meters (used as a limit for horizontal grid spacing influence
              ! on vertical scale length), ambda is (asymptotic_vertical_length_scale)^2,
              ! and vkz2 is (vonk*height_above_surface)^2.

              ambda        = max(c1*vctr2(k), min(sq300,c3))
              vkz2         = c4*zt(k)*zt(k)
              scr1(k,i,j)  = dn0(k,i,j)*vkz2/(vkz2/ambda + 1)*&
                   (sqrt(vt3di(k,i,j)) +                      &
                   enfl*sqrt(max(0., -vt3dj(k,i,j))))*        &
                   min(rchmax, sqrt(max(0.,(1.-zkhkm(ngrid)*vt3dk(k,i,j)))))

              scr2(k,i,j)  = dn0(k,i,j)*max(akm, c3*sqrt(vt3dh(k,i,j)))
              vt3dh(k,i,j) = scr1(k,i,j)*zkhkm(ngrid)
           enddo
        enddo
     enddo
     !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     !     call friclyr(nzp,nxp,nyp,a(iscr1),a(iustarl),a(itstarl)
     !    +    ,a(iustarw),a(itstarw),a(ipctlnd),a(itheta),a(irtgt))
     !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  elseif (idiffk(ngrid)==3) then
     do j=ja,jz
        do i=ia,iz
           c1 = rtgt(i,j)*rtgt(i,j)
           do k=lpw(i,j),m1-1
              scr1(k,i,j)  = dn0(k,i,j)*c1*vctr2(k)*(sqrt(vt3dh(k,i,j)) + &
                   enfl*sqrt(max(0.,-vt3dj(k,i,j))))*                     &
                   min(rchmax, sqrt(max(0.,(1.-zkhkm(ngrid)*vt3dk(k,i,j)))))
              scr2(k,i,j)  = scr1(k,i,j)
              vt3dh(k,i,j) = scr1(k,i,j)*zkhkm(ngrid)
           enddo
        enddo
     enddo
     !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     !     call friclyr(nzp,nxp,nyp,a(iscr1),a(iustarl),a(itstarl)
     !    +    ,a(iustarw),a(itstarw),a(ipctlnd),a(itheta),a(irtgt))
     !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  endif

end subroutine mxdefm

!     *****************************************************************

subroutine klbnd(m1,m2,m3,ibcon,jd,akay,dn0,lpw_R)

  implicit none

  integer, INTENT(IN)       :: m1     &
                             , m2     &
                             , m3     &
                             , ibcon  &
                             , jd

  real, dimension(m1,m2,m3), INTENT(INOUT) :: akay

  real, dimension(m1,m2,m3), INTENT(IN)    :: dn0

  real, dimension(m2,m3), INTENT(IN)    :: lpw_R

  !local variables:
  integer ::i,j,k,k2
  integer,dimension(m2,m3)  :: lpw
  !     boundary conditions on a mixing coefficient

  lpw=int(lpw_R)

  do j = 1,m3
     do i = 1,m2
        k2=lpw(i,j)
        do k=1,lpw(i,j)-1
           akay(k,i,j) = akay(k2,i,j) * dn0(k,i,j) / dn0(k2,i,j)
        enddo
        akay(m1,i,j) = akay(m1-1,i,j) * dn0(m1,i,j)  &
             / dn0(m1-1,i,j)
     enddo
  enddo
  if (iand(ibcon,1) .ne. 0) then

     do j = 1,m3
        do k = 1,m1
           akay(k,1,j) = akay(k,2,j)
        enddo
     enddo
  endif

  if (iand(ibcon,2) .ne. 0) then
     do j = 1,m3
        do k = 1,m1
           akay(k,m2,j) = akay(k,m2-1,j)
        enddo
     enddo
  endif

  if (jd .eq. 1) then
     if (iand(ibcon,4) .ne. 0) then
        do i = 1,m2
           do k = 1,m1
              akay(k,i,1) = akay(k,i,2)
           enddo
        enddo
     endif

     if (iand(ibcon,8) .ne. 0) then
        do i = 1,m2
           do k = 1,m1
              akay(k,i,m3) = akay(k,i,m3-1)
           enddo
        enddo
     endif
  endif

  return
end subroutine klbnd

!     *****************************************************************
subroutine mxdefm_tracer(m1, m2, m3, ia, iz, ja, jz, ibcon, jd,  &
     vt3dh, khtr, dn0, dxt, dyt, lpw_R, mynum)  !ALF

  !     +-------------------------------------------------------------+
  !     \   this routine calculates the mixing coefficients with a    \
  !     \     smagorinsky-type deformational based k                  \
  !     +-------------------------------------------------------------+
  !       khtr = coef. dif. horizontal para tracers
  !       frtr = fator de reducao do Akmin dos campos meteorologicos
  !
  use mem_grid, only:  &
       ngrid,          &    !INTENT(IN)

!--(DMK-CCATT-INI)-------------------------------------------------------
       ngrids
!--(DMK-CCATT-FIM)-------------------------------------------------------

  use mem_turb, only:  &
       csx,            &    !INTENT(IN)
       akmin                !INTENT(IN)

!--(DMK-CCATT-INI)-------------------------------------------------------
  use monotonic_adv, ONLY: advmnt !INTENT(IN)
!--(DMK-CCATT-FIM)-------------------------------------------------------

  implicit none
  ! Arguments:
  integer, intent(in) :: m1, m2, m3, ia, iz, ja, jz, ibcon, jd, mynum
  real, intent(in) :: lpw_R(m2,m3)
  real, intent(in)    :: dxt(m2,m3), dyt(m2,m3)    !dyt nao eh usada
  real, intent(in)    :: vt3dh(m1,m2,m3), dn0(m1,m2,m3)
  real, intent(inout) :: khtr(m1,m2,m3)
  ! local variables:
  integer :: lpw(m2,m3)
  integer :: i, j, k
  real    :: csx2, c2, c3, akm!, frtr

!--(DMK-CCATT-INI)-------------------------------------------------------
  real    :: frtr(ngrids)

  lpw=int(lpw_R)
  frtr = 0.1
  if(advmnt > 0) frtr = 0.01

!--(DMK-CCATT-OLD)-------------------------------------------------------
!  frtr = 0.05
!--(DMK-CCATT-FIM)-------------------------------------------------------

  csx2 = csx(ngrid)*csx(ngrid)

  do j=min(1,ja),max(jz,m3)
     do i=min(1,ia),max(iz,m2)

        c2  = 1.0/(dxt(i,j)*dxt(i,j))
        c3  = csx2*c2
!--(DMK-CCATT-INI)-------------------------------------------------------
        akm = frtr(ngrid) * abs(akmin(ngrid)) * 0.075 * c2 ** (0.666667)
!--(DMK-CCATT-OLD)-------------------------------------------------------
!        akm = frtr*abs(akmin(ngrid))*0.075*c2**(0.666667)
!--(DMK-CCATT-FIM)-------------------------------------------------------

        do k=lpw(i,j),m1-1
           khtr(k,i,j) = dn0(k,i,j)*max(akm, c3*sqrt(vt3dh(k,i,j)))
        enddo

     enddo
  enddo

end subroutine mxdefm_tracer


!     *****************************************************************
!srf - OLAM

subroutine bruvais_OLAM(m1, m2, m3, ia, iz, ja, jz, theta, rtp, rv, rcp, &
     pp, pi0, en2, rtgt, lpw_R)

  use mem_scratch, only: &
       vctr11,    &     !INTENT(INOUT)
       vctr12,    &     !INTENT(INOUT)
       vctr32,    &     !INTENT(INOUT)
       vctr1,     &     !INTENT(INOUT)
       vctr2,     &     !INTENT(INOUT)
       vctr3,     &     !INTENT(INOUT)
       vctr4            !INTENT(INOUT)

  use micphys, only: &
       level,     &     !INTENT(in)
       ipris,     &     !INTENT(in)
       isnow,     &     !INTENT(in)
       igraup,    &     !INTENT(in)
       iaggr,     &     !INTENT(in)
       ihail            !INTENT(in)

  use mem_grid, only: &
       zt,        &     !INTENT(in)
       nzp,       &     !INTENT(in)
       nz,        &     !INTENT(in)
       nzpmax           !INTENT(in)

  use rconstants, only: &
       alvl,      &     !INTENT(in)
       rgas,      &     !INTENT(in)
       ep,        &     !INTENT(in)
       cp,        &     !INTENT(in)
       alvi,      &     !INTENT(in)
       p00,       &     !INTENT(in)
       cpor,      &     !INTENT(in)
       g                !INTENT(in)

  implicit none
  ! Arguments:
  integer, INTENT(IN) :: m1, m2, m3, ia, iz, ja, jz
  real, INTENT(IN)    :: theta(m1,m2,m3), pp(m1,m2,m3), pi0(m1,m2,m3), &
       rtp(m1,m2,m3), rcp(m1,m2,m3), rv(m1,m2,m3)
  real, INTENT(INOUT) :: en2(m1,m2,m3)
  real, INTENT(IN)    :: rtgt(m2,m3)
  real, INTENT(IN) :: lpw_R(m2,m3)
  !local variables
  integer :: i, j, k, iweten, iwdiffk, ki, k2, k1
  real :: c1, c2, c3, ci1, ci2, ci3, rvlsi, rvii
  real :: pi(nzpmax), temp(nzpmax), prt(nzpmax), rvls(nzpmax), rc(nzpmax)
  ! **(JP)** fatora expressoes logicas para fora dos lacos
  logical :: log1, log2, log3, log4
  integer :: lpw(m2,m3)

  lpw=int(lpw_R)
  !-------------------------------------------------
  !srf- 29/12/2008:  adaptado da versao OLAM3.0
  !-------------------------------------------------
  !     calculate brunt-vaisalla frequency squared (en2)

  !iweten = 1
  !
  !iwdiffk = 0
  !c1 = alvl / rgas
  !c2 = ep * alvl ** 2 / (cp * rgas)
  !c3 = alvl / cp
  !ci1 = alvi / rgas
  !ci2 = ep * alvi ** 2 / (cp * rgas)
  !ci3 = alvi / cp

  ! **(JP)** fatora expressoes logicas para fora dos lacos
  log1 = level>=1
!  log2 = level .ge. 2 .and. iweten .eq. 1
!  log3 = (ipris  .ge. 1 .or. isnow .ge. 1 .or.  &
!       igraup .ge. 1 .or. iaggr .ge. 1 .or.  &
!       ihail  .ge. 1) .and. level .eq. 3
!  log4 = level .eq. 3
! **(JP)** fim modificacao

  !     calculate potential temperature profile

  do j=ja,jz
     do i=ia,iz

        k2 = lpw(i,j)
        k1 = k2 - 1

        do k=k1,m1
           vctr11(k) = theta(k,i,j)
           vctr12(k) = theta(k,i,j)
           !vctr32(k) = 0.
        enddo
        ! **(JP)** fatora expressoes logicas para fora dos lacos
        if (log1) then
        ! **(JP)** fim modificacao
           do k=k1,m1
              vctr12(k) = vctr11(k)*(1. + 0.61*rv(k,i,j))
              !vctr32(k) = (rtp(k,i,j) - rv(k,i,j))
           enddo
        endif

        !     check for saturation if level is 2 or greater.
        ! **(JP)** fatora expressoes logicas para fora dos lacos
        !if (log2) then
        ! **(JP)** fim modificacao
        !   do k = k1,m1
        !      pi(k) = (pp(k,i,j) + pi0(k,i,j)) / cp
        !      temp(k) = theta(k,i,j) * pi(k)
        !      prt(k) = p00 * pi(k) ** cpor
        !   enddo
        !   call mrsl(m1,prt(:),temp(:),rvls(:))
        !   do k = k2-1,m1
        !      vctr2(k) = c1
        !      vctr3(k) = c2
        !      vctr4(k) = c3
        !   enddo
        !   ki = m1 + 1
        !
        !   !	  if any ice phase microphysics are activated ....
        !   ! **(JP)** fatora expressoes logicas para fora dos lacos
        !   if (log3) then
        !   ! **(JP)** fim modificacao
        !      !  find level of -20 c.  assume ice saturation above this
        !      !   level.
	!
        !      do k = k1,m1
        !	  if (temp(k) .le. 253.16) then
        !	     ki = k
        !	     go to 10
        !	  endif
        !      enddo
        !      ki = m1 + 1
        !10            continue
        !      call mrsi(m1-ki+1,prt(ki),temp(ki),rvls(ki))
        !-srf 19/03/2005
        !bug: Subscript out of range for array prt
        !    subscript=0, lower bound=1, upper bound=132, dimension=1
        ! tempc < 253 sobre parte da Antartica.
        !              ki = max(ki,k2)                             !solucao 1
        !              call mrsi(1,prt(ki-1),temp(ki-1),rvlsi)
        !       if(ki > 1) call mrsi(1,prt(ki-1),temp(ki-1),rvlsi) ! solucao 2
        !srf
        !
        !
        !      do k = ki,m1
        !         vctr2(k) = ci1
        !         vctr3(k) = ci2
        !         vctr4(k) = ci3
        !      enddo
        !   endif
        !
        !   ! **(JP)** fatora expressoes logicas para fora dos lacos
        !   if (log4) then
        !   ! **(JP)** fim modificacao
        !      do k = k1,m1
        !         rc(k) = rcp(k,i,j)
        !      enddo
        !   else
        !      do k = k1,m1
        !         rc(k) = max(rv(k,i,j) / rvls(k) - .999,0.)
        !      enddo
        !   endif
        !
        !endif

        do k = k2,m1-1
           vctr1(k) = g/((zt(k+1) - zt(k-1))*rtgt(i,j))
        enddo

        ! **(JP)** fatora expressoes logicas para fora dos lacos
        !if (log2) then
        ! **(JP)** fim modificacao
        !   do k = k2,m1-1
        !      if (rc(k) .gt. 0.) then
        !         rvii = rvls(k-1)
        !         if (k .eq. ki) rvii = rvlsi
        !         en2(k,i,j) = vctr1(k) * (  &
        !              (1. + vctr2(k) * rvls(k) / temp(k))  &
        !              / (1. + vctr3(k) * rvls(k) / temp(k) ** 2)  &
        !              * ((vctr11(k+1) - vctr11(k-1)) / vctr11(k)  &
        !              + vctr4(k) / temp(k) * (rvls(k+1) - rvii))  &
        !              - (rtp(k+1,i,j) - rtp(k-1,i,j)))
        !      else
        !         en2(k,i,j) = vctr1(k)*((vctr12(k+1)-vctr12(k-1))  &
        !              / vctr12(k) - (vctr32(k+1) - vctr32(k-1)))
        !      endif
        !
        !   enddo
        !else
        do k=k2,m1-1
           en2(k,i,j) = vctr1(k)*((vctr12(k+1)-vctr12(k-1))/vctr12(k))
           !          / vctr12(k) - (vctr32(k+1) - vctr32(k-1)))
        enddo
        !endif
        ! **(JP)** remove para fora do laco, permitindo vetorizacao
        !en2(k1,i,j) = en2(k2,i,j)
        !en2(nzp,i,j)=en2(nz,i,j)
        ! **(JP)** fim da modificacao

     enddo
  enddo

  ! **(JP)** removido de dentro do laco
  do j=ja,jz
     do i=ia,iz
        en2(lpw(i,j)-1,i,j) = en2(lpw(i,j),i,j)
     enddo
  enddo
  do j=ja,jz
     do i=ia,iz
        en2(nzp,i,j) = en2(nz,i,j)
     enddo
  enddo
  ! **(JP)** fim da modificacao

end subroutine bruvais_OLAM

!     *****************************************************************
!-----------------------------------------------------------------------

! Begin subroutine HCV_driver (IDIFF==7)

!-----------------------------------------------------------------------
subroutine HCV_driver(m1,m2,m3,ia,iz,ja,jz,mynum  &
             ,up          &
	     ,vp          &
             ,theta       &
	     ,rv          &
	     ,dn0         &
             ,sflux_t     &
	     ,sflux_r     &
	     ,sflux_u     &
	     ,sflux_v     &
	     ,rtgt        &
	     ,lpw_R       &
	     ,scr1   	  &
	     ,vt3dh       &
	     ,vt3di       &
	     ,vt3dj       &
	     ,vt3dk       )


use mem_grid, only: zt, dzt, dzm, ngrid     !INTENT(IN)
use mem_turb, only: zkhkm
use rconstants, only:  &
                tkmin, &   !INTENT(IN)
                g,     &   !INTENT(IN)
		vonk       !INTENT(IN)


use mem_scratch, only: vctr30    !INTENT(OUT)

implicit none
real, intent(out)  , dimension(m1,m2,m3) :: scr1,vt3dh,vt3di,vt3dj,vt3dk

real, intent(in)   , dimension(m2,m3)    :: rtgt

integer, intent(in) :: m1,m2,m3,ia,iz,ja,jz,mynum

real, intent(in)   , dimension(m1,m2,m3) :: theta, &
					    rv,   &
					    up,   &
					    vp,   &
					    dn0
real, intent(in), dimension(m2,m3) :: lpw_R

real, intent(in), dimension(m2,m3) :: sflux_t, &
                                      sflux_r, &
                                      sflux_u, &
                                      sflux_v

integer :: i,j,k,k2,kzi,lpw
integer,dimension(m2,m3) :: kzi_2d
real :: sbf, zl, wstar, ustar, h


do j=ja,jz
   do i=ia,iz
      lpw=int(lpw_R(i,j))
      k2=lpw


! ---------------
! Fill column vector with virtual potential temperature
! ---------------
!      do k=k2,m1-1
      do k=1,m1-1
        vctr30(k) = theta(k,i,j) * (1. + .61 * rv(k,i,j))
      enddo
!
! ---------------
! Compute velocity scale for neutra/stable layer [ustar]
! ---------------

      ustar = max(0.1,sqrt(sqrt(sflux_u(i,j)**2 + sflux_v(i,j)**2)))

! ---------------
! Compute surface buoyancy flux [sbf] and Monin Obukhov height [zl]
! ---------------

      sbf = (g/vctr30(k2)) &

            * (sflux_t(i,j) * (1. + .61 * rv(k2,i,j)) &

	     + sflux_r(i,j) * .61 * theta(k2,i,j))


      if(sbf .ne. 0.) then
          zl = (-ustar** 3)/(vonk*sbf)
      else
          zl = 500.
      endif

!
! ---------------
! Compute boundary layer depth [h]


! Compute boundary layer depth [h] for CLE (Zilitinkevich's, 1972)

      if (zl <= 500.0 .and. zl > 0) then

           h = 2.4E3 * (ustar**(3.0/2.0))

      else

!Compute boundary layer depth [h] for CLC - using gradient Richardon number [Ri]

           call get_richgrad(h,k2,m1,zt,rtgt(i,j),vt3dj(1,i,j),vt3di(1,i,j)&
	                  ,vt3dk(1,i,j),zkhkm(ngrid),i,j)
     endif


! ---------------
! Compute convective velocity scale [wstar] / boundary layer depth [vt2da]
! ---------------

      wstar = (max(0., sbf * h)) ** (1./3.)

! ---------------
! Find the highest model level zt(k) that is at or below boundary layer
! height h.
! ---------------
      do k=k2,m1-1
         kzi = k
         if (zt(k+1)*rtgt(i,j) .gt. h) exit
      enddo
      kzi_2d(i,j)=kzi
! ---------------
! get KZZ
! ---------------

      vt3dh(:,i,j)=1.e-3
      scr1 (:,i,j)=1.e-3

      do k=k2,kzi

!******************************************************************************
!
! Call subroutine kcgcv2 for to calculate KZZ - [scr1(k,i,j]
!
!******************************************************************************

	 call kcgcv2(scr1(k,i,j), ustar, wstar, zl, zt(k)*rtgt(i,j),h, i,j)


! Kmom = scr1(k,i,j)
	 scr1(k,i,j) = max(1.e-3,scr1(k,i,j) * dn0(k,i,j))
! For Ksclr = Kmom * 3.0
	 vt3dh(k,i,j)=scr1(k,i,j)!*3.0 - done outside of this subroutine

      enddo
   enddo
enddo
scr1 (1,:,:)=scr1 (2,:,:)
vt3dh(1,:,:)=vt3dh(2,:,:)

!write(mynum,*),"kzz=",maxval(scr1),minval(scr1),maxval(kzi_2d),maxloc(kzi_2d)&
!;call flush(mynum)
!i=7
!j=15
!if(mynum==3) then
!   do k=1,kzi_2d(i,j)+1
!     write(mynum,*)"k kzz=",k,scr1(k,i,j)
!     call flush(3)
!   enddo
!endif

return
end subroutine HCV_driver

!******************************************************************
!
! Begin subroutine get_richgrad for calculate boundary layer depth [h] for CLC
!
!******************************************************************************
subroutine get_richgrad(h,k2,m1,zt,rtgt,vt3dj,vt3di,vt3dk,zkhkm,i,j)

implicit none
integer, intent (in) :: k2,m1,i,j
real, intent(out) :: h
real, intent(in)  :: zkhkm,rtgt
real, dimension(m1), intent(in)     :: zt,vt3dj,vt3di
real, dimension(m1), intent(inout)  :: vt3dk
real :: rmin,rmax,agrich,rchmax
integer k,irich
   irich = 1
   rmin = -100.
   rmax = 1. / zkhkm
   do k = k2,m1-1
      vt3dk(k) = max(min(vt3dj(k)  &
   	            / max(vt3di(k),1.e-15),rmax),rmin)
   enddo
   rchmax = 1.0 + 9.0 * float(irich)

   h=zt(k2)*rtgt

   do k = k2,m1-1
      agrich=min(rchmax,sqrt(max(0.,(1.-zkhkm*vt3dk(k)))))
      if(agrich .le. 1.e-5) then
         h=max(rtgt*zt(k2), rtgt*(zt(k)+zt(k-1))*.5)
	 exit
      endif
   enddo


return
end subroutine get_richgrad

!     *****************************************************************
! Begin subroutine kcgcv2
! The goal is to provide exchange turbulent coefficients (KZZ) for the B-RAMS system.
!     *****************************************************************

      SUBROUTINE kcgcv2 (KZZ, ustar, wstar, zl, z, h, i, j)

!################################################################
!#                                                              #
!# PURPOUSE: This routine was designed to run in  RMAS, see     #
!#           reference 4.                                       #
!#                                                              #
!#                                                              #
!# Programer: Haroldo Fraga de Campos Velho                     #
!#                                                              #
!#            Permanet Address:                                 #
!#                                                              #
!#                LAC-INPE                                      #
!#                P.O. Box 515                                  #
!#                12.201-970 - Sao Jose dos Campos (SP)         #
!#                Brasil                                        #
!#                                                              #
!#                Fax:   +55  012  345-6375                     #
!#                Phone: +55  012  345-6356                     #
!#                E-mail: haroldo@lac.inpe.br                   #
!#                Home-page: http://www.lac.inpe.br/~haroldo/   #
!#                                                              #
!#                                                              #
!# Date: Abril 15, 2000                                         #
!#       Last alteration: June 19, 2006                        #
!#                                                              #
!#                                                              #
!# REFERENCES:                                                  #
!#                                                              #
!# 1. G.A. Degrazia,  H.F. Campos Velho,  J.C. Carvalho (1997): #
!#      "Nonlocal  Exchange  Coefficients  for the  Convective  #
!#      Boundary  Layer  Derived  from  Spectral  Properties",  #
!#      Beitrage  zur  Physik der Atmosphare, Vol. 70,  No. 1,  #
!#      pp. 57-64.                                              #
!#                                                              #
!# 2. G.A. Degrazia,  O.L.L. Moraes  (1992): "A Model for Eddy  #
!#      Diffusivity  in  a  Stable Boundary Layer", Boundary-   #
!#      Layer Meteorology, Vol. 58, pp. 205-124.                #
!#                                                              #
!# 3. G.A. Degrazia, D. Anfossi,  J.C. Carvalho,  T. Tirabassi, #
!#      H.F. Campos Velho (2000): "Turbulence Parameterization  #
!#      for PBL Dispersion Models in All Stability Conditions", #
!#      Atmospheric Environment (to appear)                     #
!#                                                              #
!# 4. R.A. Pielke, W.R. Cotton, R.L. Walko, C.J. Tremback, W.A. #
!#      Lyons,  L.D. Grasso,  M.E. Nicholls,  M.D. Moran,  D.A. #
!#      Wesley, T.J. Lee, J. Copeland (1992): "A Comprehensive  #
!#      Meteorological Modeling System - RAMS", Meteorologycal  #
!#      and Atmospheric Physics, Vol. 49, pp. 69-91.            #
!#                                                              #
!################################################################
!#                                                              #
!# SUBROUTINES, FUNCTIONS REQUIRED: none                        #
!#                                                              #
!# INPUTS PARAMETERS:                                           #
!#  Z: vertical coordinate                                      #
!#  H: boundary layer height                                    #
!#  ZL: Obukhov length                                          #
!#  USTAR: velocity scale for neutral/stable layer              #
!#  WSTAR: velocity scale for convective layer                  #
!#                                                              #
!# OUTPUT PARAMETERS:                                           #
!#  KZZ: Vertical eddy diffusivity                              #
!#                                                              #
!#                                                              #
!# LOCAL VARIABLES:                                             #
!#  LAMB: local Monin-Obukhov length, see reference 2.          #
!#     Q: stability function, Eq. 2.17 in reference 1.          #
!#   PSI: Dissipation function, Eq. 2.19 in reference 1.        #
!#    FC: Coriolis parameter.                                   #
!#  CORR: Correction factor.                                    #
!#   ZLM: Average Monin-Obukhov length.                         #
!#    KN: Vertical eddy diffusivity from wind shear.            #
!#                                                              #
!################################################################


      IMPLICIT NONE

      REAL, INTENT(OUT) :: KZZ !scr1(k,i,j)
      REAL, INTENT(IN)  :: zl, wstar, h, z, ustar

      INTEGER :: i, j

      REAL LAMB, Q, PSI, FC, CORR, KN, ZIL, ZH
      REAL AUX1, AUX2, AUXL
      REAL ALP1, ALP2
      REAL RNU, DEN

! Nondimensional vertical coordinate
      ZH = min(1.,Z/H)

! Coriolis parameter (see ref. 3)
      FC = 1.0E-04

! Contribution from wind shear to turbulence (see Eq(28) in ref. 3)

      AUX1 = 0.4*(1-ZH)**0.85 * (USTAR*Z)
      AUX2 = (1.0 + 15.0*(FC*Z/USTAR))**(4.0/3.0)
      KN = AUX1/AUX2

      IF(ABS(ZL).GT.500) THEN

	 KZZ = KN

         RETURN

      ENDIF

      IF(ZL.GT.0.0) THEN

! 1 - Stable parameterization (see ref. 2)
! 1.1 - Set experimental parameters

	 IF(ZL.le.500.0) THEN

            ALP1 = 1.5
            ALP2 = 1.0

	 ENDIF

! 1.2 - Calculate the local Monin-Obukhov lenght (see Eq(3) in ref. 2)

     AUXL = (3.0*ALP1)/2.0 - ALP2
	 LAMB = ABS(ZL)*(1.0 - ZH)**AUXL

! 1.3 - Stable eddy diffusivity (see Eq(29) in ref. 3)

         RNU = 0.4 * (1.-ZH)**(3.0/4.0)
         DEN = 1. + 3.7* (Z/LAMB)
         KZZ = RNU/DEN * (USTAR*Z)

	 RETURN

      ELSE

! 2 - Convective parameterization

! The stability parameter ZIL is the ratio between the
! convective boundary layer (CBL) height and Monin-Obukhov
! length. For a full convective CBL: 50 < ZIL < 360.
! See comment in page 60 in reference 1.

! Correction factor (see Eq. (14) in ref. 3)


         CORR = SQRT(0.01*(H /(-ZL))) !See Eq. 14 4

! Here will be used the same value in KAMM model.
         ZIL = 136.0

! 2.1 - Stability function: Eq. (2.17) in reference 1.

         Q = 1.0 - EXP(-4*ZH) - 0.0003*EXP(8*ZH)


! 2.2 - Eddy diffusivity: Eq. (2.18) in reference 1.

         AUX2 = 4.0/3.0


! 2.3 - Eddy diffusivity: Eq. (26) in reference 3.
	 KZZ = 0.16* (WSTAR*H) * CORR * Q**AUX2


! Shear effect contribution (see Eq. (15), ref. 3)


          KZZ = KZZ + KN


      ENDIF
 END SUBROUTINE kcgcv2
!-----------------------------------------------------------
