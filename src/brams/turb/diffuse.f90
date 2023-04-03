subroutine diffuse_brams31()

  ! +-----------------------------------------------------------------+
  ! \     this routine is the subdriver to compute tendencies due to  \
  ! \       subgrid-scale turbulence.                                 \
  ! +-----------------------------------------------------------------+
  
  use mem_tend,  only:    &
       tend                   ! %tket, %epst, %ut, %vt, %wt

  use mem_basic, only:    &
       basic_g                ! %up (IN), %vp(IN)
 
  use var_tables, only:   &
       num_scalar,        &
       scalar_tab             ! %var_p, %var_t

  use mem_turb, only:     &
       idiffk,            &   ! INTENT(IN)
       xkhkm,             &   ! INTENT(IN)
       turb_g                 ! %tkep, %hkm, %vkh


  use mem_grid, only:     &
       jdim,              &   !        INTENT(IN)
       ngrid,             &   !        INTENT(IN)
       grid_g,            &   ! %rtgt  INTENT(IN)
       npatch,            &   !        INTENT(IN)
       zt,                &   !        INTENT(IN)
       nstbot,            &   !        INTENT(IN)
       nnzp
       
  use mem_leaf, only:     &
       leaf_g                 ! %ustar, %patch_area

       
  use mem_micro, only:    &
       micro_g                ! %rcp


  use mem_scratch, only:  &
       scratch,           &   ! %vt3da, %vt3db, %vt3dc, %vt3dd, %vt3de, 
                              ! %vt3df, %vt3dg, %vt3dh, %vt3di, %vt3dj,
                              ! %vt3dn, %scr2
       vctr34,            &   !
       vctr1,             &   !
       vctr2,             &   !
       vctr3,             &   !
       vctr4,             &   !
       vctr11,            &   !
       vctr12,            &   !
       vctr32                 !

  use node_mod, only:     &
       mxp,     &     !INTENT(IN)
       myp,     &     !INTENT(IN)
       mzp,     &     !INTENT(IN)
       ia,      &     !INTENT(IN)
       iz,      &     !INTENT(IN)
       ja,      &     !INTENT(IN)
       jz,      &     !INTENT(IN)
       ia_1,    &     !INTENT(IN)
       ja_1,    &     !INTENT(IN)
       iz1,     &     !INTENT(IN)
       jz1,     &     !INTENT(IN)
       ibcon,   &     !INTENT(IN)
       mynum,   &     !INTENT(IN)
       nodei0,  &     !INTENT(IN)
       nodej0,  &     !INTENT(IN)
       nodemyp, &     !INTENT(IN)
       nodemxp, &     !INTENT(IN)
       ia1,     &     !INTENT(IN)
       ja1,     &     !INTENT(IN)
       iz_1,    &     !INTENT(IN)
       jz_1,    &     !INTENT(IN)
       izu,     &     !INTENT(IN)
       jzv,     &     !INTENT(IN)
       mynum


  use ke_coms, only:      &
       alf_eps,   &
       alf_tke

  use micphys, only:      &
       level                !INTENT(IN)

  use mem_opt, only:      &      ! For optmization
       jstep,   &
       istep,   &
       opt                       ! %ind1_x_a, 
                                 ! %ind2_x_a, %weight_x_a, %ind1_x_b  , %ind2_x_b, %weight_x_b,
                                 ! %ind1_y_a, %ind2_y_a  , %weight_y_a, %ind2_y_b, %weight_y_b
  implicit none
  include "i8.h"
  integer(kind=i8) :: mxyzp, ind
  integer :: n
  real :: s1,s2,s3
  real, pointer :: scalarp, scalart,vkh_p,hkh_p
  integer :: i,j,k,ksf

  !! For Optimization
  integer      :: htint_i, htint_j, iia, jja, iiz, jjz, iistep, jjstep
  !! ALF
  !
  ! (JP) removing scr2 from scratch
  real, target :: scr2(mxp*myp*mzp)
  
  mxyzp = mxp*myp*mzp

  ! Shaved-Eta Coordinate not available in this optimization method
  
  
  call strain(mzp,mxp,myp,ia,iz,ja,jz                       &
       ,ia_1,ja_1,iz1,jz1,jdim                                &
       ,basic_g(ngrid)%up (1,1,1) ,basic_g(ngrid)%vp (1,1,1)  &
       ,basic_g(ngrid)%wp (1,1,1) ,scratch%vt3da     (1)      &
       ,scratch%vt3db     (1)     ,scratch%vt3dc     (1)      &
       ,scratch%vt3dd     (1)     ,scratch%vt3de     (1)      &
       ,scratch%vt3df     (1)     ,scratch%vt3dg     (1)      &
       ,scratch%vt3dh     (1)     ,scratch%vt3di     (1)      &
       ,scratch%vt3dn     (1)     ,scr2      (1)      &
       ,idiffk(ngrid))
     
  if (level<=1) &
!!$       call azero(mxyzp,scratch%vt3dp(1))
       scratch%vt3dp = 0.

  if (level>=2) &
!!$       call ae1(mxyzp,scratch%vt3dp(1),micro_g(ngrid)%rcp(1,1,1))
       call ae1_l(mxyzp, scratch%vt3dp(1), micro_g(ngrid)%rcp(1,1,1))
  
  
  call bruvais(mzp,mxp,myp,ia,iz,ja,jz                          &
       ,basic_g(ngrid)%theta (1,1,1) ,basic_g(ngrid)%rtp (1,1,1)  &
       ,basic_g(ngrid)%rv    (1,1,1) ,scratch%vt3dp(1)            &
       ,basic_g(ngrid)%pp    (1,1,1) ,basic_g(ngrid)%pi0 (1,1,1)  &
       ,scratch%vt3dj        (1)     ,grid_g(ngrid)%rtgt (1,1)    &
       ,grid_g(ngrid)%lpw    (1,1)   )
  
  if (idiffk(ngrid) <= 3) then
     call mxdefm(mzp,mxp,myp,ia,iz,ja,jz,ibcon,jdim            &
          ,scratch%vt3dh      (1)     ,scratch%vt3di      (1)    &
          ,scratch%vt3dj      (1)     ,scratch%vt3dk      (1)    &
          ,scratch%scr1       (1)     ,scr2       (1)    &
          ,basic_g(ngrid)%dn0 (1,1,1) ,grid_g(ngrid)%rtgt (1,1)  &
          ,grid_g(ngrid)%dxt  (1,1)   ,grid_g(ngrid)%dyt  (1,1)  &
          ,grid_g(ngrid)%lpw  (1,1)   ,mynum  )
  endif
  
  if (idiffk(ngrid) == 1) then
     call tkemy(mzp,mxp,myp,ia,iz,ja,jz,ibcon,jdim,nodei0(mynum,ngrid),nodej0(mynum,ngrid)  &
          ,turb_g(ngrid)%tkep   (1,1,1) ,tend%tket            (1)      &
          ,scratch%vt3dh        (1)     ,scratch%vt3di        (1)      &
          ,scratch%vt3dj        (1)     ,scratch%scr1         (1)      &
          ,grid_g(ngrid)%rtgt   (1,1)   ,basic_g(ngrid)%theta (1,1,1)  &
          ,basic_g(ngrid)%dn0   (1,1,1) ,basic_g(ngrid)%up    (1,1,1)  &
          ,basic_g(ngrid)%vp    (1,1,1) ,basic_g(ngrid)%wp    (1,1,1)  &
          ,turb_g(ngrid)%sflux_u(1,1)   ,turb_g(ngrid)%sflux_v(1,1)    &
          ,turb_g(ngrid)%sflux_w(1,1)   ,turb_g(ngrid)%sflux_t(1,1),vctr34 &
          ,grid_g(ngrid)%lpw    (1,1)   ,grid_g(ngrid)%lpu    (1,1)   &
          ,grid_g(ngrid)%lpv    (1,1))
  endif
  
  if (idiffk(ngrid) == 4) then
     call mxtked(mzp,mxp,myp,ia,iz,ja,jz  &
          ,ibcon,jdim  &
          ,turb_g(ngrid)%tkep   (1,1,1) ,tend%tket            (1)      &
          ,basic_g(ngrid)%up    (1,1,1) ,basic_g(ngrid)%vp    (1,1,1)  &
          ,basic_g(ngrid)%wp    (1,1,1) ,basic_g(ngrid)%rtp   (1,1,1)  &
          ,basic_g(ngrid)%rv    (1,1,1) ,basic_g(ngrid)%theta (1,1,1)  &
          ,scratch%vt3da        (1)     ,scratch%vt3dc        (1)      &
          ,scratch%vt3dh        (1)     ,scratch%vt3dj        (1)      &
          ,scratch%scr1         (1)     ,scr2         (1)      &
          ,turb_g(ngrid)%sflux_u(1,1)   ,turb_g(ngrid)%sflux_v(1,1)    &
          ,turb_g(ngrid)%sflux_w(1,1)   ,turb_g(ngrid)%sflux_t(1,1)    &
          ,grid_g(ngrid)%dxt    (1,1)   ,grid_g(ngrid)%rtgt   (1,1)    &
          ,grid_g(ngrid)%lpw    (1,1)   )
  endif
  
  !_STC............................................................
  !_STC Call to subroutine tkescl for E-l closure
  !_STC (S. Trini Castelli)
  !_STC............................................................
  if (idiffk(ngrid) == 5) then
     call tkescl(mzp,mxp,myp,npatch,ia,iz,ja,jz  &
          ,turb_g(ngrid)%tkep(1,1,1),tend%tket(1)  &
          ,turb_g(ngrid)%epsp(1,1,1),tend%epst(1)  &
          ,scratch%vt3da(1),scratch%vt3dc(1)  &
          ,scratch%vt3dh(1),scratch%vt3di(1)  &
          ,scratch%vt3dj(1),scratch%scr1(1)  &
          ,scr2(1) ,grid_g(ngrid)%rtgt(1,1)  &
          ,scratch%vt3dd(1),scratch%vt3de(1),grid_g(ngrid)%dxt(1,1)  &
          ,leaf_g(ngrid)%ustar(1,1,1),leaf_g(ngrid)%patch_area(1,1,1) &
          ,grid_g(ngrid)%lpw(1,1),basic_g(ngrid)%dn0(1,1,1)  )
  endif
  !_STC............................................................
  !_STC Call to subroutine tkeeps for E-eps closure
  !_STC (S. Trini Castelli)
  !_STC............................................................
  if (idiffk(ngrid) == 6) then
     call tkeeps(mzp,mxp,myp,npatch,ia,iz,ja,jz  &
          ,turb_g(ngrid)%tkep(1,1,1),tend%tket(1)  &
          ,turb_g(ngrid)%epsp(1,1,1),tend%epst(1)  &
          ,scratch%vt3da(1),scratch%vt3dc(1)  &
          ,scratch%vt3dh(1),scratch%vt3di(1)  &
          ,scratch%vt3dj(1),scratch%scr1(1)  &
          ,scr2(1) ,grid_g(ngrid)%rtgt(1,1)  &
          ,leaf_g(ngrid)%ustar(1,1,1),leaf_g(ngrid)%patch_area(1,1,1) &
          ,grid_g(ngrid)%lpw(1,1),basic_g(ngrid)%dn0(1,1,1)  )
  endif
  !_STC..................................................
  !_STC    Note: from subroutines TKESCL, TKEEPS :
  !_STC           VT3DI=Ke
  !_STC           SCR1=Km
  !_STC           VT3DH = Kh
  !_STC           SCR2 = SCR1 = Km
  !_STC..................................................
  !_STC............................................................
  
  call klbnd(mzp,mxp,myp,ibcon,jdim  &
       ,scratch%scr1 (1),basic_g(ngrid)%dn0(1,1,1),grid_g(ngrid)%lpw(1,1))
  call klbnd(mzp,mxp,myp,ibcon,jdim  &
       ,scr2 (1),basic_g(ngrid)%dn0(1,1,1),grid_g(ngrid)%lpw(1,1))
  call klbnd(mzp,mxp,myp,ibcon,jdim  &
       ,scratch%vt3dh(1),basic_g(ngrid)%dn0(1,1,1),grid_g(ngrid)%lpw(1,1))
  !_STC ....... boundary conditions even on Ke diffusion coefficient
  if(idiffk(ngrid) ==  5 .or. idiffk(ngrid) == 6) &
       call klbnd(mzp,mxp,myp,ibcon,jdim  &
       ,scratch%vt3di(1),basic_g(ngrid)%dn0(1,1,1),grid_g(ngrid)%lpw(1,1))
  
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
           !!         scratch%vt3di(ind) = turb_g(ngrid)%vke(k,i,j)
           turb_g(ngrid)%hkm(k,i,j) = s1
           turb_g(ngrid)%vkm(k,i,j) = s2
           turb_g(ngrid)%vkh(k,i,j) = s3
        enddo
     enddo
  enddo
  
  call diffvel(mzp,mxp,myp,ia,iz,ja,jz,jdim,ia_1,ja_1             &
       ,ia1,ja1,iz_1,jz_1,iz1,jz1,izu,jzv,idiffk(ngrid)             &
       ,basic_g(ngrid)%up    (1,1,1) ,basic_g(ngrid)%vp    (1,1,1)  &
       ,basic_g(ngrid)%wp    (1,1,1) ,tend%ut              (1)      &
       ,tend%vt              (1)     ,tend%wt              (1)      &
       ,scratch%vt3da        (1)     ,scratch%vt3db        (1)      &
       ,scratch%vt3dc        (1)     ,scratch%vt3dd        (1)      &
       ,scratch%vt3de        (1)     ,scratch%vt3df        (1)      &
       ,scratch%vt3dg        (1)     ,scratch%vt3dj        (1)      &
       ,scratch%vt3dk        (1)     ,scratch%vt3dl        (1)      &
       ,scratch%vt3dm        (1)     ,scratch%vt3dn        (1)      &
       ,scratch%vt3do        (1)     ,grid_g(ngrid)%rtgu   (1,1)    &
       ,grid_g(ngrid)%rtgv   (1,1)   ,grid_g(ngrid)%rtgt   (1,1)    &
       ,turb_g(ngrid)%sflux_u(1,1)   ,turb_g(ngrid)%sflux_v(1,1)    &
       ,turb_g(ngrid)%sflux_w(1,1)   ,basic_g(ngrid)%dn0   (1,1,1)  &
       ,basic_g(ngrid)%dn0u  (1,1,1) ,basic_g(ngrid)%dn0v  (1,1,1)  &
       ,scratch%scr1         (1)     ,scr2         (1),ibcon,mynum)
     
  ! Convert momentum K's to scalar K's, if necessary
  
  if (idiffk(ngrid) <= 3) then
     do ind = 1,mxyzp
        scr2(ind) = scr2(ind) * xkhkm(ngrid)
     enddo
  elseif (idiffk(ngrid) == 4) then
     do ind = 1,mxyzp
        scratch%vt3di(ind) = 2. * scratch%scr1(ind)
     enddo
  endif


  ! Calculating Index and Weights for interpolation (htint)

  ! Setting indexes and Weights for HTINT with Optmization

  ! Loop for a block of spatial local region

  jjstep = min(jstep,(jz-ja))
  iistep = min(istep,(iz-ia))

  do jja = ja, jz, jjstep

     jjz = MIN(jja+jjstep-1, jz) ! Calculating last element

     do iia = ia, iz, iistep

        iiz = MIN(iia+iistep-1, iz) ! Calculating last element

        ! Compute the Indexes and Weights for x_dir - ALF
        do j=jja,jjz
           do i=iia,iiz
              do k=1,mzp !m1
                 vctr1(k)=grid_g(ngrid)%topt(i-1,j)+  &
                      zt(k)*grid_g(ngrid)%rtgt(i-1,j)

                 vctr2(k)=grid_g(ngrid)%topt(i  ,j)+  &
                      zt(k)*grid_g(ngrid)%rtgt(i ,j)

                 vctr3(k)=grid_g(ngrid)%topt(i+1,j)+  &
                      zt(k)*grid_g(ngrid)%rtgt(i+1,j)
              enddo

              htint_i = i-iia+1
              htint_j = j-jja+1

              call htint_index(mzp, vctr1, mzp, vctr2,  &
                   opt%ind1_x_a(1,htint_i,htint_j),   &
                   opt%ind2_x_a(1,htint_i,htint_j),   &
                   opt%weight_x_a(1,htint_i,htint_j))

              call htint_index(mzp, vctr3, mzp, vctr2,  &
                   opt%ind1_x_b(1,htint_i,htint_j),   &
                   opt%ind2_x_b(1,htint_i,htint_j),   &
                   opt%weight_x_b(1,htint_i,htint_j))


              ! Compute the Indexes and Weights for y_dir - ALF

              do k=1,mzp !m1
                 vctr1(k)=grid_g(ngrid)%topt(i,j-jdim)+  &
                      zt(k)*grid_g(ngrid)%rtgt(i,j-jdim)

                 vctr3(k)=grid_g(ngrid)%topt(i,j+jdim)+  &
                      zt(k)*grid_g(ngrid)%rtgt(i,j+jdim)
              enddo

              call htint_index(mzp, vctr1, mzp, vctr2,  &
                   opt%ind1_y_a(1,htint_i,htint_j),   &
                   opt%ind2_y_a(1,htint_i,htint_j),   &
                   opt%weight_y_a(1,htint_i,htint_j))

              call htint_index(mzp, vctr3, mzp, vctr2,  &
                   opt%ind1_y_b(1,htint_i,htint_j),   &
                   opt%ind2_y_b(1,htint_i,htint_j),   &
                   opt%weight_y_b(1,htint_i,htint_j))

           enddo
        enddo
        ! ALF

        !! End of calculating Index and Weights for interpolation (htint)

  
        do n = 1,num_scalar(ngrid)

           scalarp => scalar_tab(n,ngrid)%var_p
           scalart => scalar_tab(n,ngrid)%var_t
     
!!$           call azero(mxp*myp,scratch%vt2da(1))
           scratch%vt2da = 0.
           if (nstbot == 1) then
              if (scalar_tab(n,ngrid)%name == 'THP' .or. &
	          scalar_tab(n,ngrid)%name == 'THC' ) then
                 call atob(mxp*myp, turb_g(ngrid)%sflux_t(1,1), scratch%vt2da(1))
              elseif (scalar_tab(n,ngrid)%name == 'RTP') then
                 call atob(mxp*myp, turb_g(ngrid)%sflux_r(1,1), scratch%vt2da(1))
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
     
           if (scalar_tab(n,ngrid)%name == 'TKEP') then
              vkh_p => scratch%vt3di(1)
              hkh_p => scr2(1)
              if (idiffk(ngrid) >= 4) hkh_p => scratch%vt3di(1)
              ksf = 1   
           elseif (scalar_tab(n,ngrid)%name == 'EPSP') then
              vkh_p => scratch%vt3di(1)
              hkh_p => scr2(1)
              if (idiffk(ngrid) >= 4)  hkh_p => scratch%vt3di(1)
              ksf = 3
              ! Convert Ktke to Keps; it will be converted back after use below
!!$              call ae1t0 (mxyzp,vkh_p,vkh_p,ALF_EPS/ALF_TKE)
!!$              call ae1t0 (mxyzp,hkh_p,hkh_p,ALF_EPS/ALF_TKE)
              call ae1t0_l(mxyzp, vkh_p, vkh_p, (ALF_EPS/ALF_TKE))
              call ae1t0_l(mxyzp, hkh_p, hkh_p, (ALF_EPS/ALF_TKE))
           else
              vkh_p => scratch%vt3dh(1)
              hkh_p => scr2(1)
              if (idiffk(ngrid) >= 4) hkh_p => scratch%vt3dh(1)
              ksf = 2
           endif
     
     
           call diffsclr_brams31(mzp,mxp,myp,iia,iiz,jja,jjz,jdim        &
                ,ia_1,ja_1,ia1,ja1,iz_1,jz_1,iz1,jz1,n,ksf               &
                ,scalarp,scalart            ,scratch%vt3da(1)            &
                ,scratch%vt3db      (1)     ,scratch%vt3df      (1)      &
                ,scratch%vt3dg      (1)     ,scratch%vt3dj      (1)      &
                ,scratch%vt3dk      (1)     ,scratch%vt3do      (1)      &
                ,scratch%vt3dc      (1)     ,scratch%vt3dd      (1)      &
                ,scratch%vt3dl      (1)     ,scratch%vt3dm      (1)      &
                ,scratch%vt2db      (1)     ,grid_g(ngrid)%rtgt (1,1)    &
                ,scratch%vt2da      (1)     ,basic_g(ngrid)%dn0 (1,1,1)  &
                ,vkh_p                      ,hkh_p                       )

           if (scalar_tab(n,ngrid)%name == 'EPSP') then
!!$              call ae1t0 (mxyzp,vkh_p,vkh_p,ALF_TKE/ALF_EPS)
!!$              call ae1t0 (mxyzp,hkh_p,hkh_p,ALF_TKE/ALF_EPS)
              call ae1t0_l(mxyzp, vkh_p, vkh_p, (ALF_TKE/ALF_EPS))
              call ae1t0_l(mxyzp, hkh_p, hkh_p, (ALF_TKE/ALF_EPS))
           endif
     
        enddo

     enddo

  enddo
  
  return
end subroutine diffuse_brams31
