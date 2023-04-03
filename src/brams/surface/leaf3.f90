!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine sfclyr(mzp,mxp,myp,ia,iz,ja,jz,ibcon)

  use mem_all

  ! TEB_SPM
  use teb_spm_start, only: TEB_SPM    ! INTENT(IN)
  USE mem_teb, only: &
       teb_g,        &        !Data type
       teb_vars               !Type
  USE mem_teb_common, only: &
       tebc_g,              & !Data Type
       teb_common             !Type

  implicit none

  !Arguments:
  integer, intent(in) :: mzp,mxp,myp,ia,iz,ja,jz,ibcon

  !Local Variables
  real :: rslif
  integer :: ng
  !integer, save :: ncall=0
  !Pointer to TEB data. Can be associated or not
  type(teb_vars), pointer   :: p_teb_g
  type(teb_common), pointer :: p_tebc_g
  ! New automatic arrays for use as scratch instead of:
  ! scratch%vt2da, scratch%vt2db, scratch%vt2dc, scratch%vt2dd,
  ! scratch%vt2de, scratch%vt2df, scratch%vt3da
  real, dimension(mxp,myp) :: l_ths2, l_rvs2, l_pis2, l_dens2, &
       l_ups2, l_vps2, l_zts2

  ! Interface necessary to use pointer as argument - TEB
  interface
     subroutine leaf3(m1,m2,m3,mzg,mzs,np,ia,iz,ja,jz    &
          ,leaf,basic,turb,radiate,grid,cuparm,micro     &
          ,ths2,rvs2,pis2,dens2,ups2,vps2,zts2           &
          ! For TEB
          ,pteb,ptebc                                    &
          !
          )
       ! for interface
       use mem_leaf, only: leaf_vars  ! type
       use mem_basic, only: basic_vars! type
       use mem_turb, only: turb_vars! type
       use mem_radiate, only: radiate_vars! type
       use mem_grid, only: grid_vars! type
       use mem_cuparm, only: cuparm_vars! type
       use mem_micro, only: micro_vars! type
       use mem_teb, only: teb_vars               !Type
       USE mem_teb_common, only: teb_common             !Type
       !
       integer, intent(in) :: m1,m2,m3,mzg,mzs,np,ia,iz,ja,jz
       type (leaf_vars)    :: leaf
       type (basic_vars)   :: basic
       type (turb_vars)    :: turb
       type (radiate_vars) :: radiate
       type (grid_vars)    :: grid
       type (cuparm_vars)  :: cuparm
       type (micro_vars)   :: micro
       ! For TEB
       TYPE (teb_vars), pointer   :: pteb
       TYPE (teb_common), pointer :: ptebc
       real, dimension(m2,m3), intent(out) :: ths2,rvs2,pis2,dens2,ups2,vps2,zts2
     end subroutine leaf3
  end interface

  if (nstbot == 0) return

  !print*,'ncall=',ncall
  !if(ncall == 0) then
  !   ncall=1
  !   open(11,file='leaf.lis', status='unknown')
  !   rewind 11

  !print*,'calling alloc=',ncall
  !   call alloc_leafcol(nzg,nzs)

  !endif

  ng=ngrid

  ! TEB
  if (TEB_SPM==1) then
     p_teb_g  => teb_g(ngrid)
     p_tebc_g => tebc_g(ngrid)
  else
     nullify(p_teb_g)
     nullify(p_tebc_g)
  endif

  call leaf3(mzp,mxp,myp,nzg,nzs,npatch,ia,iz,ja,jz             &
       ,leaf_g (ng), basic_g (ng), turb_g (ng), radiate_g(ng)   &
       ,grid_g (ng), cuparm_g(ng), micro_g(ng)                  &
       ,l_ths2(1,1), l_rvs2(1,1), l_pis2(1,1)                   &
       ,l_dens2(1,1),l_ups2(1,1), l_vps2(1,1)                   &
       ,l_zts2(1,1)                                             &
       ! For TEB
       ,p_teb_g, p_tebc_g                                       &
       !
       )

  if (isfcl == 2) then
     call hydro(mxp,myp,nzg,nzs,npatch         &
          ,leaf_g(ng)%soil_water      (1,1,1,1)  &
          ,leaf_g(ng)%soil_energy     (1,1,1,1)  &
          ,leaf_g(ng)%soil_text       (1,1,1,1)  &
          ,leaf_g(ng)%sfcwater_mass   (1,1,1,1)  &
          ,leaf_g(ng)%sfcwater_energy (1,1,1,1)  &
          ,leaf_g(ng)%patch_area      (1,1,1)    &
          ,leaf_g(ng)%patch_wetind    (1,1,1)    )
  endif

  ! Apply lateral boundary conditions to leaf3 arrays

  call leaf_bcond(mxp,myp,nzg,nzs,npatch,jdim     &
       ,leaf_g(ng)%soil_water (1,1,1,1) ,leaf_g(ng)%sfcwater_mass  (1,1,1,1)  &
       ,leaf_g(ng)%soil_energy(1,1,1,1) ,leaf_g(ng)%sfcwater_energy(1,1,1,1)  &
       ,leaf_g(ng)%soil_text  (1,1,1,1) ,leaf_g(ng)%sfcwater_depth (1,1,1,1)  &
       ,leaf_g(ng)%ustar        (1,1,1) ,leaf_g(ng)%tstar            (1,1,1)  &
       ,leaf_g(ng)%rstar        (1,1,1) ,leaf_g(ng)%veg_albedo       (1,1,1)  &
       ,leaf_g(ng)%veg_fracarea (1,1,1) ,leaf_g(ng)%veg_lai          (1,1,1)  &
       ,leaf_g(ng)%veg_tai      (1,1,1)                                       &
       ,leaf_g(ng)%veg_rough    (1,1,1) ,leaf_g(ng)%veg_height       (1,1,1)  &
       ,leaf_g(ng)%patch_area   (1,1,1) ,leaf_g(ng)%patch_rough      (1,1,1)  &
       ,leaf_g(ng)%patch_wetind (1,1,1) ,leaf_g(ng)%leaf_class       (1,1,1)  &
       ,leaf_g(ng)%soil_rough   (1,1,1) ,leaf_g(ng)%sfcwater_nlev    (1,1,1)  &
       ,leaf_g(ng)%stom_resist  (1,1,1) ,leaf_g(ng)%ground_rsat      (1,1,1)  &
       ,leaf_g(ng)%ground_rvap  (1,1,1) ,leaf_g(ng)%veg_water        (1,1,1)  &
       ,leaf_g(ng)%veg_temp     (1,1,1) ,leaf_g(ng)%can_rvap         (1,1,1)  &
       ,leaf_g(ng)%can_temp     (1,1,1) ,leaf_g(ng)%veg_ndvip        (1,1,1)  &
       ,leaf_g(ng)%veg_ndvic    (1,1,1) ,leaf_g(ng)%veg_ndvif        (1,1,1)  )

  return
end subroutine sfclyr

!*****************************************************************************

subroutine leaf3(m1,m2,m3,mzg,mzs,np,ia,iz,ja,jz  &
     ,leaf,basic,turb,radiate,grid,cuparm,micro     &
     ,ths2,rvs2,pis2,dens2,ups2,vps2,zts2           &
                                ! For TEB
     ,pteb,ptebc                                    &
                                !
     )

  use mem_leaf, only: leaf_vars, isfcl, dthcon, drtcon, pctlcon
  use mem_basic, only: basic_vars
  use mem_turb, only: turb_vars
  use mem_radiate, only: radiate_vars, ilwrtyp, iswrtyp
  use mem_grid, only: grid_vars, zt, zm, ngrid, time, dtlt, if_adap, jdim, dzt, istp, itime1
  use mem_micro, only: micro_vars
  use mem_cuparm, only: cuparm_vars, nnqparm
  use micphys, only: level
  use io_params, only: iupdsst, ssttime1, ssttime2
  use leaf_coms, only: timefac_sst, niter_leaf, niter_can, dtll, dtll_factor, &
       dtlc_factor, dtlc, hcapcan, wcapcan, hcapveg, dtllohcc, dtllowcc, &
       dtlcohcc, dtlcowcc, dtlcohcv, z0fac_water, snowrough, ups, vps, ths, rvs, &
       zts, pis, dens, prss, vels, gzotheta, tempk, snowfac, thetacan, &
       vels_pat, rdi, pcpgl
  use rconstants, only: g, cpor, p00, vonk, cp, cpi, rocp, alvl
  use node_mod, only : mynum  ! INTENT(IN)

!--(DMK-CCATT-INI)-----------------------------------------------------------
  use ccatt_start, only: &
       ccatt
!--(DMK-CCATT-OLD)-----------------------------------------------------------
!  use catt_start, only: catt  ! INTENT(IN)
!--(DMK-CCATT-FIM)-----------------------------------------------------------

  use teb_spm_start, only: teb_spm    ! INTENT(IN)
  use mem_teb, only: teb_vars               !Type
  use mem_teb_common, only: teb_common             !Type




  implicit none

  ! Arguments:
  integer, intent(in) :: m1,m2,m3,mzg,mzs,np,ia,iz,ja,jz
  type (leaf_vars)    :: leaf
  type (basic_vars)   :: basic
  type (turb_vars)    :: turb
  type (radiate_vars) :: radiate
  type (grid_vars)    :: grid
  type (cuparm_vars)  :: cuparm
  type (micro_vars)   :: micro
  ! For TEB
  TYPE (teb_vars), pointer   :: pteb
  TYPE (teb_common), pointer :: ptebc
  real, dimension(m2,m3), intent(out) :: ths2,rvs2,pis2,dens2,ups2,vps2,zts2

  ! Local variables:
  integer :: i,j,ip,iter_leaf
  real :: rslif
  integer :: k2
  real :: dvelu,dvelv,velnew,sflux_uv,cosine1,sine1
  ! Local variables used for TEB
  real :: PSUP1,PSUP2,depe,alt2,deze,dpdz,exn1ST,airt
  real :: ZH_TOWN,ZLE_TOWN,ZSFU_TOWN,ZSFV_TOWN
  real, pointer :: EMIS_TOWN, ALB_TOWN, TS_TOWN
  real, pointer :: G_URBAN

  INTERFACE

     SUBROUTINE sfcrad(mzg, mzs, ip,                  &
          soil_energy, soil_water, soil_text,         &
          sfcwater_energy, sfcwater_depth,            &
          patch_area, can_temp, veg_temp, leaf_class, &
          veg_height, veg_fracarea, veg_albedo,       &
          sfcwater_nlev, rshort, rlong, albedt,       &
          rlongup, cosz,                              &
          G_URBAN, ETOWN, ALBTOWN, TSTOWN             )
       integer, intent(IN) :: mzg, mzs, ip
       real, intent(IN)    :: soil_energy(mzg)
       real, intent(IN)    :: soil_water(mzg)
       real, intent(IN)    :: soil_text(mzg)
       real, intent(IN)    :: sfcwater_energy(mzs)
       real, intent(IN)    :: sfcwater_depth(mzs)
       real, intent(IN)    :: patch_area
       real, intent(IN)    :: can_temp
       real, intent(IN)    :: veg_temp
       real, intent(IN)    :: leaf_class
       real, intent(IN)    :: veg_height
       real, intent(IN)    :: veg_fracarea
       real, intent(IN)    :: veg_albedo
       real, intent(IN)    :: sfcwater_nlev
       real, intent(IN)    :: rshort
       real, intent(IN)    :: rlong
       real, intent(INOUT) :: albedt
       real, intent(INOUT) :: rlongup
       real, intent(IN)    :: cosz
       real, pointer, optional :: G_URBAN, ETOWN, ALBTOWN, TSTOWN
     END SUBROUTINE sfcrad

     subroutine sfclmcv(i,j,ustar, tstar, rstar, vels, vels_pat, ups, vps, &
          gzotheta, patch_area, sflux_u, sflux_v, sflux_w, sflux_t, sflux_r)
       real, intent(in)    :: ustar, tstar, rstar, vels,vels_pat, ups, vps, &
            gzotheta, patch_area
       real, intent(inout) :: sflux_u, sflux_v, sflux_w, sflux_t, sflux_r
       !real, pointer, optional :: G_URBAN
       integer, intent(in) :: i,j
     end subroutine sfclmcv

  END INTERFACE

  dtll=0.0
  !TEB
  nullify(G_URBAN)
  nullify(EMIS_TOWN)
  nullify(ALB_TOWN)
  nullify(TS_TOWN)

  ! Time interpolation factor for updating SST

  if (iupdsst == 0) then
     timefac_sst = 0.
  else
     timefac_sst = (time - ssttime1(ngrid)) / (ssttime2(ngrid) - ssttime1(ngrid))
  endif

  ! Define leaf3 and canopy time-split timesteps here.  This ensures that leaf3
  ! will not use a timestep longer than about 40 seconds, and canopy will not
  ! use a timestep longer than about 15 seconds.  This allows values of
  ! hcapcan = 2.e4, wcapcan = 2.e1, and hcapveg = 3.e4 as are now defined below.

  niter_leaf  = max(1,nint(dtlt/40.+.4))
  niter_can   = max(1,nint(dtll/15.+.4))
  dtll_factor = 1. / float(niter_leaf)
  dtll        = dtlt * dtll_factor
  dtlc_factor = 1. / float(niter_can)
  dtlc        = dtll * dtlc_factor

  hcapcan = 2.0e4
  wcapcan = 2.0e1
  hcapveg = 3.e4

  dtllohcc = dtll / hcapcan
  dtllowcc = dtll / wcapcan
  dtlcohcc = dtlc / hcapcan
  dtlcowcc = dtlc / wcapcan
  dtlcohcv = dtlc / hcapveg

  z0fac_water = .016 / g
  snowrough = .001


  ! Copy surface atmospheric variables into 2d arrays for input to leaf

  if (if_adap == 1) then
     call sfc_fields_adap(m1,m2,m3,ia,iz,ja,jz,jdim             &
          ,grid%lpu   (1,1)   ,grid%lpv (1,1)   ,grid%lpw(1,1)    &
          ,grid%topma (1,1)   ,grid%aru (1,1,1) ,grid%arv(1,1,1)  &
          ,basic%theta(1,1,1) ,basic%rv (1,1,1) ,basic%up(1,1,1)  &
          ,basic%vp   (1,1,1) ,basic%dn0(1,1,1) ,basic%pp(1,1,1)  &
          ,basic%pi0  (1,1,1) ,zt,zm,dzt                          &
          ,ths2,rvs2,ups2,vps2,pis2,dens2,zts2                    )
  else
     call sfc_fields(m1,m2,m3,ia,iz,ja,jz,jdim                  &
          ,basic%theta(1,1,1) ,basic%rv (1,1,1) ,basic%up(1,1,1)  &
          ,basic%vp   (1,1,1) ,basic%dn0(1,1,1) ,basic%pp(1,1,1)  &
          ,basic%pi0  (1,1,1) ,grid%rtgt(1,1)   ,zt               &
          ,ths2,rvs2,ups2,vps2,pis2,dens2,zts2                    )
  endif


  do j = ja,jz
     do i = ia,iz

        ! Copy surface variables to single-column values

        ups = ups2(i,j)
        vps = vps2(i,j)
        ths = ths2(i,j)
        rvs = rvs2(i,j)
        zts = zts2(i,j)
        pis = pis2(i,j)
        dens = dens2(i,j)

        prss = pis ** cpor * p00
        vels = sqrt(ups ** 2 + vps ** 2)
        gzotheta = g * zts / ths

        ! Update water internal energy from time-dependent SST

        leaf%soil_energy(mzg,i,j,1) = 334000.  &
             + 4186. * (leaf%seatp(i,j) + (leaf%seatf(i,j) - leaf%seatp(i,j))  &
             * timefac_sst - 273.15)

        ! Fill surface precipitation arrays for input to leaf

        call sfc_pcp(nnqparm(ngrid),level,i,j,cuparm,micro)

        ! Zero out albedo, upward surface longwave, and momentum, heat, and moisture
        ! flux arrays before summing over patches

        if (ilwrtyp > 0 .or. iswrtyp > 0) then
           radiate%albedt(i,j) = 0.
           radiate%rlongup(i,j) = 0.
        endif

        turb%sflux_u(i,j) = 0.
        turb%sflux_v(i,j) = 0.
        turb%sflux_w(i,j) = 0.
        turb%sflux_t(i,j) = 0.
        turb%sflux_r(i,j) = 0.

        ! For no soil model (patch 2) fill "canopy" temperature and moisture

        if (isfcl == 0) then
           leaf%can_temp  (i,j,2) = (ths - dthcon) * pis
           leaf%can_rvap  (i,j,2) = rvs - drtcon
           leaf%patch_area(i,j,1) = 1. - pctlcon
           leaf%patch_area(i,j,2) = pctlcon
        endif

        ! Begin patch loop
        do ip = 1,np

           ! Update time-dependent vegetation LAI and fractional coverage
           if (ip >= 2 .and. leaf%patch_area(i,j,ip) >= .009) then

              if (ip >= 2 .and. isfcl >= 1) call vegndvi(ngrid        &
                   ,leaf%patch_area  (i,j,ip) ,leaf%leaf_class(i,j,ip)  &
                   ,leaf%veg_fracarea(i,j,ip) ,leaf%veg_lai   (i,j,ip)  &
                   ,leaf%veg_tai     (i,j,ip) ,leaf%veg_rough (i,j,ip)  &
                   ,leaf%veg_height  (i,j,ip) ,leaf%veg_albedo(i,j,ip)  &
                   ,leaf%veg_ndvip   (i,j,ip) ,leaf%veg_ndvic (i,j,ip)  &
                   ,leaf%veg_ndvif   (i,j,ip)                           )

           endif

           ! Begin leaf small timestep here.
           do iter_leaf = 1,niter_leaf

              ! Calculate radiative fluxes between atmosphere, vegetation, and ground/snow
              ! based on already-computed downward shortwave and longwave fluxes from
              ! the atmosphere.  Fill tempk array with soil and snow temperature (C) and
              ! fracliq array with liquid fraction of water content in soil and snow.
              ! Other snowcover properties are also computed here.

              if (iswrtyp > 0 .or. ilwrtyp > 0) then

                 if (ip == 1 .or. leaf%patch_area(i,j,ip) >= .009) then

                    ! TEB - defining pointers
                    if (TEB_SPM==1) then
                       G_URBAN   => leaf%G_URBAN(i,j,ip)
                       EMIS_TOWN => ptebc%EMIS_TOWN(i,j)
                       ALB_TOWN  => ptebc%ALB_TOWN(i,j)
                       TS_TOWN   => ptebc%TS_TOWN(i,j)
                       call sfcrad(mzg, mzs, ip,                              &
                            leaf%soil_energy(1:mzg,i,j,ip),                   &
                            leaf%soil_water(1:mzg,i,j,ip),                    &
                            leaf%soil_text(1:mzg,i,j,ip),                     &
                            leaf%sfcwater_energy(1:mzs,i,j,ip),               &
                            leaf%sfcwater_depth(1:mzs,i,j,ip),                &
                            leaf%patch_area(i,j,ip),                          &
                            leaf%can_temp(i,j,ip), leaf%veg_temp(i,j,ip),     &
                            leaf%leaf_class(i,j,ip), leaf%veg_height(i,j,ip), &
                            leaf%veg_fracarea(i,j,ip),                        &
                            leaf%veg_albedo(i,j,ip),                          &
                            leaf%sfcwater_nlev(i,j,ip),                       &
                            radiate%rshort(i,j), radiate%rlong(i,j),          &
                            radiate%albedt(i,j), radiate%rlongup(i,j),        &
                            radiate%cosz(i,j),                                &
                            ! For TEB
                            G_URBAN, EMIS_TOWN, ALB_TOWN, TS_TOWN             &
                            !
                            )
                    else
                       call sfcrad(mzg, mzs, ip,                              &
                            leaf%soil_energy(1:mzg,i,j,ip),                   &
                            leaf%soil_water(1:mzg,i,j,ip),                    &
                            leaf%soil_text(1:mzg,i,j,ip),                     &
                            leaf%sfcwater_energy(1:mzs,i,j,ip),               &
                            leaf%sfcwater_depth(1:mzs,i,j,ip),                &
                            leaf%patch_area(i,j,ip),                          &
                            leaf%can_temp(i,j,ip), leaf%veg_temp(i,j,ip),     &
                            leaf%leaf_class(i,j,ip), leaf%veg_height(i,j,ip), &
                            leaf%veg_fracarea(i,j,ip),                        &
                            leaf%veg_albedo(i,j,ip),                          &
                            leaf%sfcwater_nlev(i,j,ip),                       &
                            radiate%rshort(i,j), radiate%rlong(i,j),          &
                            radiate%albedt(i,j), radiate%rlongup(i,j),        &
                            radiate%cosz(i,j)                                 &
                            )
                    endif

!!$                    call sfcrad(mzg, mzs, ip,                                 &
!!$                         leaf%soil_energy(1:mzg,i,j,ip),                      &
!!$                         leaf%soil_water(1:mzg,i,j,ip),                       &
!!$                         leaf%soil_text(1:mzg,i,j,ip),                        &
!!$                         leaf%sfcwater_energy(1:mzs,i,j,ip),                  &
!!$                         leaf%sfcwater_depth(1:mzs,i,j,ip),                   &
!!$                         leaf%patch_area(i,j,ip),                             &
!!$                         leaf%can_temp(i,j,ip), leaf%veg_temp(i,j,ip),        &
!!$                         leaf%leaf_class(i,j,ip), leaf%veg_height(i,j,ip),    &
!!$                         leaf%veg_fracarea(i,j,ip), leaf%veg_albedo(i,j,ip),  &
!!$                         leaf%sfcwater_nlev(i,j,ip),                          &
!!$                         radiate%rshort(i,j), radiate%rlong(i,j),             &
!!$                         radiate%albedt(i,j), radiate%rlongup(i,j),           &
!!$                         radiate%cosz(i,j),                                   &
!!$                         ! For TEB
!!$                         G_URBAN, EMIS_TOWN, ALB_TOWN, TS_TOWN                &
!!$                         !
!!$                         )

                 endif

              endif

              ! For water surface (patch 1), compute surface saturation mixing ratio
              ! and roughness length based on previous ustar.
              ! For soil patches, compute roughness length based on vegetation and snow.

              if (ip == 1) then

                 leaf%ground_rsat(i,j,ip) = rslif(prss,tempk(mzg))
                 leaf%patch_rough(i,j,ip)  &
                      = max(z0fac_water * leaf%ustar(i,j,ip) ** 2,.0001)

              elseif (isfcl >= 1) then

                 if (leaf%patch_area(i,j,ip) >= .009) then
                    leaf%patch_rough(i,j,ip)   &
                         = max(grid%topzo(i,j),leaf%soil_rough(i,j,ip)  &
                         ,leaf%veg_rough(i,j,ip)) * (1. - snowfac)   &
                         + snowrough * snowfac

                 endif

              endif

              ! Calculate turbulent fluxes between atmosphere and canopy (or "canopy")

              if (leaf%patch_area(i,j,ip) < .009 .and. isfcl == 1 .and. ip >= 2) then
                 thetacan = ths
              else
                 thetacan = leaf%can_temp(i,j,ip) / pis
              endif

              call stars(leaf%ustar(i,j,ip),leaf%tstar(i,j,ip)  &
                   ,leaf%rstar(i,j,ip),ths,rvs,thetacan,leaf%can_rvap(i,j,ip)  &
                   ,zts,leaf%patch_rough(i,j,ip),leaf%patch_area(i,j,ip)  &
                                !kml drydep
                   ,vels,vels_pat,vonk,dtllohcc,dens,dtll             &
                   ,leaf%R_aer(i,j,ip))
              !kml drydep

              ! TEB
              if (TEB_SPM==1) then
                 G_URBAN => leaf%G_URBAN(i,j,ip)
                !  call sfclmcv(i,j,leaf%ustar(i,j,ip), leaf%tstar(i,j,ip),         &
                !       leaf%rstar(i,j,ip), vels, vels_pat, ups, vps, gzotheta, &
                !       leaf%patch_area(i,j,ip), turb%sflux_u(i,j),             &
                !       turb%sflux_v(i,j), turb%sflux_w(i,j),                   &
                !       turb%sflux_t(i,j), turb%sflux_r(i,j),                   &
                !       ! For TEB
                !       G_URBAN                                                 &
                !       !
                !       )
              else
                 call sfclmcv(i,j,leaf%ustar(i,j,ip), leaf%tstar(i,j,ip),         &
                      leaf%rstar(i,j,ip), vels, vels_pat, ups, vps, gzotheta, &
                      leaf%patch_area(i,j,ip), turb%sflux_u(i,j),             &
                      turb%sflux_v(i,j), turb%sflux_w(i,j),                   &
                      turb%sflux_t(i,j), turb%sflux_r(i,j)                    &
                      )
              endif
!call dumpVarAllLatLonk(turb%sflux_w,'Lsflux_w',520,ip,iter_leaf,ia-1,iz+1,ja-1,jz+1,1,1,0.0,0.0) !
!!$              call sfclmcv(leaf%ustar(i,j,ip),leaf%tstar(i,j,ip)     &
!!$                   ,leaf%rstar(i,j,ip),vels,vels_pat,ups,vps,gzotheta  &
!!$                   ,leaf%patch_area(i,j,ip),turb%sflux_u(i,j)          &
!!$                   ,turb%sflux_v(i,j),turb%sflux_w(i,j)                &
!!$                   ,turb%sflux_t(i,j),turb%sflux_r(i,j)                &
!!$                                ! For TEB
!!$                   ,G_URBAN                                            &
!!$                                !
!!$                   )

              ! For water patches, update temperature and moisture of "canopy" from
              ! divergence of fluxes with water surface and atmosphere.  rdi = ustar/5
              ! is the viscous sublayer conductivity from Garratt (1992).

              if (ip == 1) then

                 rdi = .2 * leaf%ustar(i,j,1)

                 leaf%can_temp(i,j,1) = leaf%can_temp(i,j,1)        &
                      + dtllohcc * dens * cp                          &
                      * ((tempk(mzg) - leaf%can_temp(i,j,1)) * rdi    &
                      + leaf%ustar(i,j,1) * leaf%tstar(i,j,1) * pis)

                 leaf%can_rvap(i,j,1) = leaf%can_rvap(i,j,1) + dtllowcc * dens       &
                      * ((leaf%ground_rsat(i,j,1) - leaf%can_rvap(i,j,1)) * rdi  &
                      + leaf%ustar(i,j,1) * leaf%rstar(i,j,1))

                 if (CCATT==1) then
                    leaf%R_aer(i,j,ip) = rdi
                 endif

              endif

              if (TEB_SPM==1) then

                 !TEB Urban canopy parameterization starts here
                 if (ip >= 2) then

                    if (NINT(leaf%G_URBAN(i,j,ip))/=0.) then

                       !edmilson
                       PSUP1  = ((basic%PI0(1,I,J) + basic%pp(1,i,j))*cpi)** &
                            CPOR*P00
                       PSUP2  = ((basic%PI0(2,I,J) + basic%pp(2,i,j))*cpi)** &
                            CPOR*P00
                       depe   = psup2 - psup1
                       alt2   = zt(1)*grid%rtgt(i,j)
                       deze   = zts - alt2
                       dpdz   = depe/deze
                       exn1ST = (psup2/p00)**rocp

                       airt   = basic%theta(2,i,j)*exn1st

                       ! TEB - defining pointers
                       if (TEB_SPM==1) then
                          G_URBAN => leaf%G_URBAN(i,j,ip)
                       endif

                       CALL LEAF3_TEB_INTERFACE(ISTP,DTLT,DTLL,          &
                            radiate%COSZ(i,j), ZTS,                      &
                            radiate%rlong(i,j), radiate%rshort(i,j),     &
                            psup2, airt, ups, vps, basic%rv(2,i,j),      &
                            pcpgl/dtlt, pteb%fuso(i,j),                  &
                            pteb%T_CANYON(i,j), pteb%R_CANYON(i,j),      &
                            pteb%TS_ROOF(i,j), pteb%TS_ROAD(i,j),        &
                            pteb%TS_WALL(i,j), pteb%TI_ROAD(i,j),        &
                            pteb%TI_BLD(i,j), pteb%WS_ROOF(i,j),         &
                            pteb%WS_ROAD(i,j), pteb%T_ROOF(2:4,i,j),     &
                            pteb%T_ROAD(2:4,i,j), pteb%T_WALL(2:4,i,j),  &
                            ZH_TOWN, ZLE_TOWN, ptebc%EMIS_TOWN(i,j),     &
                            ZSFU_TOWN, ZSFV_TOWN, ptebc%TS_TOWN(i,j),    &
                            ptebc%ALB_TOWN(i,j), NINT(G_URBAN),          &
                            pteb%H_TRAFFIC(i,j), pteb%H_INDUSTRY(i,j),   &
                            pteb%LE_TRAFFIC(i,j), pteb%LE_INDUSTRY(i,j), &
                            pteb%T2M_TOWN(i,j), pteb%R2M_TOWN(i,j),      &
                            time, itime1, dpdz, dens                     )

                       turb%sflux_u(i,j) = &
                            turb%sflux_u(i,j) + leaf%patch_area(i,j,ip)*ZSFU_TOWN
                       turb%sflux_v(i,j) = &
                            turb%sflux_v(i,j) + leaf%patch_area(i,j,ip)*ZSFV_TOWN
                       turb%sflux_t(i,j) = &
                            turb%sflux_t(i,j) + &
                            leaf%patch_area(i,j,ip)*ZH_TOWN/(CP*DENS)
                       turb%sflux_r(i,j) = &
                            turb%sflux_r(i,j) + &
                            leaf%patch_area(i,j,ip)*ZLE_TOWN/(ALVL*DENS)

                    endif
                 endif
              endif

              ! For soil model patches, update temperature and moisture of soil,
              ! vegetation, and canopy

              if (isfcl >= 1 .and. ip >= 2) then

                 if (leaf%patch_area(i,j,ip) >= .009) then

                    call leaftw(mzg,mzs,np,              &
                         leaf%soil_water(1,i,j,ip),      &
                         leaf%soil_energy(1,i,j,ip),     &
                         leaf%soil_text(1,i,j,ip),       &
                         leaf%sfcwater_mass(1,i,j,ip),   &
                         leaf%sfcwater_energy(1,i,j,ip), &
                         leaf%sfcwater_depth (1,i,j,ip), &
                         leaf%ustar(i,j,ip),             &
                         leaf%tstar(i,j,ip),             &
                         leaf%rstar(i,j,ip),             &
                         leaf%veg_albedo(i,j,ip),        &
                         leaf%veg_fracarea(i,j,ip),      &
                         leaf%veg_lai(i,j,ip),           &
                         leaf%veg_tai(i,j,ip),           &
                         leaf%veg_rough(i,j,ip),         &
                         leaf%veg_height(i,j,ip),        &
                         leaf%patch_area(i,j,ip),        &
                         leaf%patch_rough(i,j,ip),       &
                         leaf%patch_wetind(i,j,ip),      &
                         leaf%leaf_class(i,j,ip),        &
                         leaf%soil_rough(i,j,ip),        &
                         leaf%sfcwater_nlev(i,j,ip),     &
                         leaf%stom_resist(i,j,ip),       &
                         leaf%ground_rsat(i,j,ip),       &
                         leaf%ground_rvap(i,j,ip),       &
                         leaf%veg_water(i,j,ip),         &
                         leaf%veg_temp(i,j,ip),          &
                         leaf%can_rvap(i,j,ip),          &
                         leaf%can_temp(i,j,ip),          &
                         leaf%veg_ndvip(i,j,ip),         &
                         leaf%veg_ndvic(i,j,ip),         &
                         leaf%veg_ndvif(i,j,ip),         &
                         radiate%rshort(i,j),            &
                         radiate%cosz(i,j),              &
                         ip,i,j                          )

                 endif
              endif
           enddo
        enddo

     enddo
  enddo

  ! Normalize accumulated fluxes and albedo seen by atmosphere over model
  ! timestep dtlt.


  do j = ja,jz
     do i = ia,iz
        turb%sflux_u(i,j) = turb%sflux_u(i,j) * dtll_factor * dens2(i,j)
        turb%sflux_v(i,j) = turb%sflux_v(i,j) * dtll_factor * dens2(i,j)
        turb%sflux_w(i,j) = turb%sflux_w(i,j) * dtll_factor * dens2(i,j)
        turb%sflux_t(i,j) = turb%sflux_t(i,j) * dtll_factor * dens2(i,j)
        turb%sflux_r(i,j) = turb%sflux_r(i,j) * dtll_factor * dens2(i,j)
     enddo
  enddo
!call dumpVarAllLatLonk(turb%sflux_w,'Lsflux_w',683,0,0,ia-1,iz+1,ja-1,jz+1,1,1,0.0,0.0) !
  if (ilwrtyp > 0 .or. iswrtyp > 0) then
     do j = ja,jz
        do i = ia,iz
           radiate%albedt (i,j) = radiate%albedt (i,j) * dtll_factor
           radiate%rlongup(i,j) = radiate%rlongup(i,j) * dtll_factor
        enddo
     enddo
  endif

  return
end subroutine leaf3

!***************************************************************************

subroutine leaftw(mzg,mzs,np  &
   ,soil_water     ,soil_energy      ,soil_text       &
   ,sfcwater_mass  ,sfcwater_energy  ,sfcwater_depth  &
   ,ustar          ,tstar            ,rstar           &
   ,veg_albedo     ,veg_fracarea     ,veg_lai         &
   ,veg_tai                                           &
   ,veg_rough      ,veg_height       ,patch_area      &
   ,patch_rough    ,patch_wetind     ,leaf_class      &
   ,soil_rough     ,sfcwater_nlev    ,stom_resist     &
   ,ground_rsat    ,ground_rvap      ,veg_water       &
   ,veg_temp       ,can_rvap         ,can_temp        &
   ,veg_ndvip      ,veg_ndvic        ,veg_ndvif       &
   ,rshort         ,cosz             ,ip,i,j              )

use leaf_coms
use mem_leaf
use rconstants
use mem_scratch

use ccatt_start, only: &
       ccatt

implicit none

integer :: mzg,mzs,np

real, dimension(mzg) :: soil_water,soil_energy,soil_text
real, dimension(mzs) :: sfcwater_mass,sfcwater_energy,sfcwater_depth

real                 :: ustar        ,tstar         ,rstar        &
                       ,veg_albedo   ,veg_fracarea  ,veg_lai      &
                       ,veg_tai                                   &
                       ,veg_rough    ,veg_height    ,patch_area   &
                       ,patch_rough  ,patch_wetind  ,leaf_class   &
                       ,soil_rough   ,sfcwater_nlev ,stom_resist  &
                       ,ground_rsat  ,ground_rvap   ,veg_water    &
                       ,veg_temp     ,can_rvap      ,can_temp     &
                       ,veg_ndvip    ,veg_ndvic     ,veg_ndvif    &
                       ,rshort       ,cosz

real, dimension(nzgmax,nstyp) :: slcons1
real, dimension(nstyp) :: slcons0,fhydraul
common/efold/slcons1,slcons0,fhydraul

integer :: ip,k,nveg,ksn,nsoil,ksnnew,newlayers,nlayers,kold,ktrans,nsl,kzs
integer, save :: ncall=0

integer :: i,j

!-----------------------------------------------------------------------------
! parameters for new soil heat conductivity (8/17/00):  Move to sfcdata later
real, dimension(12) :: soilcond0,soilcond1,soilcond2
data soilcond0/ 0.30, 0.30, 0.29, 0.27, 0.28, 0.28  &
               ,0.26, 0.27, 0.27, 0.25, 0.25, 0.06/
data soilcond1/ 4.80, 4.66, 4.27, 3.47, 3.63, 3.78  &
               ,2.73, 3.23, 3.32, 2.58, 2.40, 0.46/
data soilcond2/-2.70,-2.60,-2.31,-1.74,-1.85,-1.96  &
              ,-1.20,-1.56,-1.63,-1.09,-0.96, 0.00/
!-----------------------------------------------------------------------------

real :: stretch,thik,pf,snden,vegfracc,qwfree,wfree          &
   ,depthgain,totsnow,qw,w,qwt,wt,soilhcap,fac,wfreeb,depthloss,soilcap    &
   ,sndenmax,sndenmin,snowmin,hold,wtnew,wtold,wdiff,watermid,availwat,wg  &
   ,wloss,soilcond,waterfrac

real, save, dimension(20) :: thicknet
real, save, dimension(20,20) :: thick

!---------srf-05052006--------------------------- root profiles
integer, parameter :: iroot=1
real, dimension(nzgmax) ::fswpk
real fswp_equiv
! use mem_grid, only: ngrid
! use extras, only:   extra2d,extra3d, NA_EXTRA2D
!
!---------srf-05052006--------------------------- root profiles




do k = 1,mzg
   dslz   (k) = slz(k+1) - slz(k)
   dslzi  (k) = 1. / dslz(k)
   dslzidt(k) = dslzi(k) * dtll
   slzt   (k) = .5 * (slz(k) + slz(k+1))
enddo

do k = 2,mzg
   dslzt   (k) = slzt(k) - slzt(k-1)
   dslzti  (k) = 1. / dslzt(k)
   dslztidt(k) = dslzti(k) * dtll
enddo

! Initialize snow thickness scaling array

if (ncall /= 40) then
   ncall = 40
   stretch = 2.0
   do kzs = 1,mzs
      thik = 1.0
      thicknet(kzs) = 0.0
      do k = 1,(kzs+1)/2
         thick(k,kzs) = thik
         thick(kzs+1-k,kzs) = thik
         thicknet(kzs) = thicknet(kzs) + 2. * thik
         thik = thik * stretch
      enddo
      if ((kzs+1)/2 /= kzs/2) thicknet(kzs) = thicknet(kzs) - thik/stretch
      do k = 1,kzs
         thick(k,kzs) = thick(k,kzs) / thicknet(kzs)
      enddo
   enddo
endif

nveg = nint(leaf_class)
ksn = nint(sfcwater_nlev)

! Evaluate any exchanges of heat and moisture to or from vegetation, apply
! moisture and heat changes to vegetation, and evaluate the resistance
! parameter rd between canopy air and the top soil or snow surface.

call canopy(mzg,mzs,ksn,nveg  &
   ,soil_energy,soil_water,soil_text,sfcwater_mass  &
   ,ustar,tstar,rstar,soil_rough,veg_rough,veg_height  &
   ,veg_lai,veg_tai,veg_water,veg_temp,leaf_class,veg_fracarea  &
   ,stom_resist,can_temp,can_rvap,ground_rsat,ground_rvap,rshort  &
!   ,i,j,ip)
!---------srf-05052006---------------------------
   ,i,j,ip,iroot,fswpk,fswp_equiv)
!---------srf-05052006---------------------------

! Compute soil and effective snow heat resistance times layer depth (rfactor).

do k = 1,mzg
   nsoil = nint(soil_text(k))
   waterfrac = soil_water(k) / slmsts(nsoil)
   soilcond = soilcond0(nsoil)  &
      + waterfrac * (soilcond1(nsoil) + waterfrac * soilcond2(nsoil))
   rfactor(k) = dslz(k) / soilcond
enddo

do k = 1,ksn
   snden = sfcwater_mass(k) / sfcwater_depth(k)
   rfactor(k+mzg) = sfcwater_depth(k)  &
      / (1.093e-3 * exp(.028 * tempk(k+mzg)) * (.030 + snden  &
      * (.303e-3 + snden * (-.177e-6 + snden * 2.25e-9))))
enddo

! Find soil and snow internal sensible heat fluxes [W/m2]

hfluxgsc(1) = 0.

do k = 2,mzg+ksn
   hfluxgsc(k) = - (tempk(k) - tempk(k-1)) / ((rfactor(k) + rfactor(k-1)) * .5)
enddo

! Heat flux at soil or snow top from longwave, sensible, and
! upward latent heat fluxes [W/m^2]

hfluxgsc(ksn+1+mzg) = hflxgc + wflxgc * alvl - rlonga_gs  &
   - rlongv_gs + rlonggs_v + rlonggs_a

! Update soil Q values [J/m3] and snow Q values [J/kg] from sensible heat,
! upward water vapor (latent heat), longwave, and shortwave fluxes.
! This excludes effects of dew/frost formation, precipitation, shedding,
! and percolation.  Update top soil or snow moisture from evaporation only.

do k = 1,mzg
   soil_energy(k) = soil_energy(k) + dslzidt(k) * (hfluxgsc(k) - hfluxgsc(k+1))
enddo

soil_energy(mzg) = soil_energy(mzg) + dslzidt(mzg) * rshort_g

do k = 1,ksn
   sfcwater_energy(k) = sfcwater_energy(k) + dtll &
      * (hfluxgsc(k+mzg) - hfluxgsc(k+1+mzg) + rshort_s(k)) / sfcwater_mass(k)
enddo

if (ksn == 0) then
   !DSM soil_water(mzg) = soil_water(mzg) - 1.e-3 * wflxgc * dslzidt(mzg)
   soil_water(mzg) = max(0.07,soil_water(mzg) - 1.e-3 * wflxgc * dslzidt(mzg))
else
   sfcwater_mass(ksn) = max(0.,sfcwater_mass(ksn) - wflxgc * dtll)
endif

! New moisture, qw, and depth from dew/frost formation, precipitation,
! shedding, and percolation.  ksnnew is the layer that receives the new
! condensate that comes directly from the air above.  If there is no
! pre-existing snowcover, this is a temporary "snow" layer.

if (pcpgl + wshed + dewgnd > 1.e-9) then
   ksnnew = max(ksn,1)
   vegfracc = 1. - veg_fracarea
   qwfree = dewgnd * alvl + qpcpgl * vegfracc + qwshed * veg_fracarea
   wfree = dewgnd + pcpgl * vegfracc + wshed * veg_fracarea
   depthgain = dpcpgl * vegfracc + (dewgnd + wshed) * veg_fracarea * .001
else
   ksnnew = ksn
   qwfree = 0.
   wfree = 0.
   depthgain = 0.
endif

if (ksnnew > 0) then

! Transfer water downward through snow layers by percolation.
! Fracliq is the fraction of liquid in the snowcover or surface water.  wfree
! is the quantity of that liquid in kg/m2 which is free (not attached to
! snowcover) and therefore available to soak into the layer below).
! soilcap is the capacity of the top soil layer in kg/m2 to accept surface
! water.  wfree in the lowest snow layer is limited by this value.

   totsnow = 0.
   nsoil = nint(soil_text(mzg))

   do k = ksnnew,1,-1
      qw = sfcwater_energy(k) * sfcwater_mass(k) + qwfree
      w = sfcwater_mass(k) + wfree

! If (only) snow layer is too thin for computational stability, bring
! it to thermal equilibrium with the top soil layer by exchanging
! heat between them.

      if (ksnnew == 1 .and. sfcwater_mass(k) < 3.) then
         qwt = qw + soil_energy(mzg) * dslz(mzg)
         wt = w + soil_water(mzg) * 1.e3 * dslz(mzg)
         soilhcap = slcpd(nsoil) * dslz(mzg)
         call qwtk(qwt,wt,soilhcap,tempk(k+mzg),fracliq(k+mzg))
         fac = 4186.
         if (fracliq(k+mzg) <= .0001) fac = 2093.
         qw = (fac * (tempk(k+mzg) - 273.15) + fracliq(k+mzg) * 334000.) * w
         tempk(mzg) = tempk(k+mzg)
         fracliq(mzg) = fracliq(k+mzg)
         soil_energy(mzg) = (qwt - qw) * dslzi(mzg)
      else
         call qwtk(qw,w,100.,tempk(k+mzg),fracliq(k+mzg))
      endif

! Shed liquid in excess of a 1:9 liquid-to-ice ratio.  Limit this shed amount
! (wfreeb) in lowest snow layer to amount top soil layer can hold.

      wfreeb = max (0.,w * (fracliq(k+mzg) - .1) / 0.9)
      depthloss = wfreeb * 1.e-3
      if (k == 1) then
         soilcap = 1.e3 * max (0.,-slz(mzg) * (slmsts(nsoil) - soil_water(mzg)))
         wfreeb = min (wfreeb, soilcap)
         qwfree = wfreeb * 4186. * (tempk(k+mzg) - 193.36)
         soil_water(mzg) = soil_water(mzg) + 1.e-3 * wfreeb * dslzi(mzg)
         soil_energy(mzg) = soil_energy(mzg) + qwfree * dslzi(mzg)
      endif
      qwfree = wfreeb * 4186. * (tempk(k+mzg) - 193.36)

      !srf-25/02/2005
      sfcwater_mass(k) = w - wfreeb ! from RAMS 6.0
      !sfcwater_energy(k) = w - wfreeb ! old way


      sfcwater_depth(k) = sfcwater_depth(k) + depthgain - depthloss

      totsnow = totsnow + sfcwater_mass(k)
      sfcwater_energy(k) = (qw - qwfree) / (max(1.e-9,sfcwater_mass(k)))
      sfcwater_energy(k) = max (-1.6e5, min (4.8e5, sfcwater_energy(k)))

! Temporary simple evolution of snow layer depth and density

      sfcwater_depth(k) = sfcwater_depth(k) * (1. - dtll / 1.e5)
      snden = sfcwater_mass(k) / max(1.e-6,sfcwater_depth(k))
      sndenmax = 1000.
      sndenmin = max(30.,200. * (wfree + wfreeb) / max(1.e-9,sfcwater_mass(k)))
      snden = min (sndenmax, max (sndenmin, snden))
      sfcwater_depth(k) = sfcwater_mass(k) / snden
      wfree = wfreeb
      depthgain = depthloss
   enddo

! Re-distribute snow layers to maintain prescribed distribution of mass

   if (totsnow < 1.e-9) then
      sfcwater_nlev = 0.
      sfcwater_mass(1) = 0.
      sfcwater_energy(1) = 0.
      sfcwater_depth(1) = 0.
   else
      nlayers = ksnnew
      snowmin = 3.0
      newlayers = 1
      do k = 2,mzs
         if (snowmin * thicknet(k) <= totsnow .and.  &
            sfcwater_energy(k) < 334000.) newlayers = newlayers + 1
      enddo
      newlayers = min (newlayers, mzs, nlayers+1)
      sfcwater_nlev = float(newlayers)
      kold = 1
      wtnew = 1.
      wtold = 1.
      do k = 1,newlayers
         vctr14(k) = totsnow * thick(k,newlayers)
         vctr16(k) = 0.
         vctr18(k) = 0.
10       continue
         wdiff = wtnew * vctr14(k) - wtold * sfcwater_mass(kold)
         if (wdiff > 0.) then
            vctr16(k) = vctr16(k) + wtold * sfcwater_mass(kold)  &
               * sfcwater_energy(kold)
            vctr18(k) = vctr18(k) + wtold * sfcwater_depth(kold)
            wtnew = wtnew - wtold * sfcwater_mass(kold) / vctr14(k)
            kold = kold + 1
            wtold = 1.
            if (kold <= ksn) go to 10
         else
            vctr16(k) = vctr16(k) + wtnew * vctr14(k) * sfcwater_energy(kold)

            !DSM vctr18(k) = vctr18(k) + wtnew * vctr14(k) * sfcwater_depth(kold)  &
            !DSM    / max(1.e-9,sfcwater_mass(kold))

            vctr18(k) = vctr18(k) + wtnew * vctr14(k) * (sfcwater_depth(kold)  &
               / max(1.e-9,sfcwater_mass(kold)))
            wtold = wtold - wtnew * vctr14(k) / sfcwater_mass(kold)
            wtnew = 1.
         endif
      enddo

      do k = 1,newlayers
         sfcwater_mass(k) = vctr14(k)
         sfcwater_energy(k) = vctr16(k) / sfcwater_mass(k)
         sfcwater_depth(k) = vctr18(k)
      enddo

   endif
endif

! Compute gravitational potential plus moisture potential z + psi [m],
! liquid water content [m], and half the remaining water capacity [m].

do k = 1,mzg
   nsoil = nint(soil_text(k))
   psiplusz(k) = slzt(k) + slpots(nsoil)   &
               * (slmsts(nsoil) / soil_water(k)) ** slbs(nsoil)
   soil_liq(k) = dslz(k) * min(soil_water(k) - soilcp(nsoil)  &
      ,soil_water(k) * fracliq(k) )
   half_soilair(k) = (slmsts(nsoil) - soil_water(k)) * dslz(k) * .5
enddo

! Find amount of water transferred between soil layers (wflux) [m]
! modulated by the liquid water fraction

wflux(1) = 0.
wflux(mzg+1) = 0.
qwflux(1) = 0.
qwflux(mzg+1) = 0.

do k = 2,mzg
   nsoil = nint(soil_text(k))
   watermid = 0.5 * (soil_water(k) + soil_water(k-1))
   wflux(k) = dslztidt(k) * slcons1(k,nsoil)  &
      * (watermid / slmsts(nsoil)) ** (2. * slbs(nsoil) + 3.)  &
      * (psiplusz(k-1) - psiplusz(k)) * .5 * (fracliq(k) + fracliq(k-1))

! Limit water transfers to prevent over-saturation and over-depletion
! Compute q transfers between soil layers (qwflux) [J/m2]

   if (wflux(k) > 0.) then
      wflux(k) = min(wflux(k),soil_liq(k-1),half_soilair(k))
   else
      wflux(k) = - min(-wflux(k),soil_liq(k),half_soilair(k-1))
   endif

   qwflux(k) = wflux(k) * 4.186e6 * (tempk(k) - 193.36)

enddo

! Update soil moisture (impose minimum value of soilcp) and q value.

do k = 1,mzg
   nsoil = nint(soil_text(k))
   !
   !soil_water(k) = max(soilcp(nsoil),soil_water(k)  &
   !   - dslzi(k) * (wflux(k+1) - wflux(k)))
   !srf-rams60 mod
   soil_water(k) = min(slmsts(nsoil), max(soilcp(nsoil)  &
                  ,soil_water(k)- dslzi(k) * (wflux(k+1) - wflux(k))))

   soil_energy(k) = soil_energy(k) - dslzi(k) * (qwflux(k+1) - qwflux(k))
enddo

! Remove water only from moistest level in root zone for transpiration.
! Bottom k-index in root zone is kroot(nveg).  Limit soil moisture
! to values above soilcp.  More sophisticated use of root
! profile function and water extractibility function, and
! improved minimum value are being considered.
! Units of wloss are m3/m3, of transp are kg/m2/s.
!
!cccccccccccccccccccccccccccccc  TRANSPIRATION bypass  ccccccccccccccccc
!      go to 4044
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!---------srf-2006-2009--------------------------- root profiles
if( iroot == 1 .and. transp > 1.e-18) then

  if(maxval(fswpk(1:mzg)) == 0.) stop "wrong fswpk"

  wg=0.
  do k = 1,mzg
     nsoil = nint(soil_text(k))
     availwat = soil_water(k) * fracliq(k)
     wloss = min(&
             transp*dslzidt(k)*1.e-3*root(nveg,k)*fswp_equiv/fswpk(k) &
            ,availwat                                  &
            ,soil_water(k) - soilcp(nsoil))

     soil_water (k) = soil_water (k) - wloss


     wg=wg+wloss/(dslzidt(k) * 1.e-3)! total water lost by transpiration

     soil_energy(k) = soil_energy(k)  &
  		      - wloss * 4.186e6 * (tempk(k) - 193.36)

  enddo

else  ! RAMS original way

wg = 0.
nsl = nint(soil_text(mzg))
ktrans = mzg
do k = kroot(nveg),mzg
   nsoil = nint(soil_text(k))
   availwat = soil_water(k) * fracliq(k)
   if (wg < availwat) then
      wg = availwat
      ktrans = k
      nsl = nsoil
   endif
enddo

wloss = min(transp * dslzidt(ktrans) * 1.e-3,wg   &
       ,soil_water(ktrans) - soilcp(nsl))

soil_water(ktrans) = soil_water(ktrans) - wloss
soil_energy(ktrans) = soil_energy(ktrans)  &
   - wloss * 4.186e6 * (tempk(ktrans) - 193.36)


endif

!---------srf-2006-2009--------------------------- root profiles

4044 continue

! Compute ground vap mxrat for availability on next timestep; put into
! ground_rsat.

call grndvap(soil_energy(mzg),soil_water(mzg),soil_text(mzg)  &
   ,sfcwater_energy(mzs),sfcwater_nlev,ground_rsat            &
   ,ground_rvap,can_temp,can_rvap,prss                        )

return
end subroutine leaftw

!***************************************************************************

subroutine canopy(mzg,mzs,ksn,nveg  &
   ,soil_energy,soil_water,soil_text,sfcwater_mass  &
   ,ustar,tstar,rstar,soil_rough,veg_rough,veg_height &
   ,veg_lai,veg_tai,veg_water,veg_temp,leaf_class,veg_fracarea  &
   ,stom_resist,can_temp,can_rvap,ground_rsat,ground_rvap,rshort  &
!   ,i,j,ip)
   ,i,j,ip,iroot,fswpk,fswp_equiv)
!---------srf-05052006---------------------------

use leaf_coms
use rconstants

!--(DMK-CCATT-INI)-----------------------------------------------------------
  use ccatt_start, only: &
       ccatt
!--(DMK-CCATT-FIM)-----------------------------------------------------------

!!$!CATT
!!$use catt_start, only: CATT

implicit none

integer :: mzg,mzs,ksn,nveg

real, dimension(mzg) :: soil_energy,soil_water,soil_text
real, dimension(mzs) :: sfcwater_mass
real :: ustar,tstar,rstar,soil_rough,veg_rough,veg_height             &
       ,veg_lai,veg_tai,veg_water,veg_temp,leaf_class,veg_fracarea    &
       ,stom_resist,can_temp,can_rvap,ground_rsat,ground_rvap,rshort

integer :: k,kk,nsoil,iter_can

integer :: i,j,ip

real :: aux,brad,bswp,bthi,btlo,bvpd,c2,c3,c4,dsm,es,eslf,fac,factv         &
   ,fracliqv,fthi,ftlo,frad,fswp,fvpd,qwtot,rasgnd,rasveg,rleaf,rsatveg     &
   ,sigmaw ,slai,stai,slpotv,srad,sswp,sthi,stlo,svpd,swp,tveg,tvegc,tvegk  &
   ,vpd,wtemp,wtroot,x,zognd,zoveg,zdisp,zveg,wflx,dewgndflx,ustar0         &
   ,transp_test,rc
real*8 :: rc_inf, qwshed_tmp, aux1
! From RAMS 6.0
real :: wshed0, transp0

!     Note: c1=261.5*sqrt((1.-exp(-2.*exar))/(2.*exar))
!     from Lee's dissertation, Eq. 3.36.  The factor of 261.5 is
!     100 * ln((h-d)/zo) / vonk   where d = .63 * h and zo = .13 * h.
!     The factor of 100 is 1/L in Eq. 3.37.  Thus, c1 * ustar is the
!     total expression inside the radical in Eq. 3.37.
!bob      parameter(exar=3.5,covr=2.16,c1=98.8)
real, parameter :: exar=2.5,covr=2.16,c1=116.6

data   brad,  srad/    196.,   0.047/,  &
       btlo,  stlo/   281.5,    0.26/,  &
       bthi,  sthi/   310.1,  -0.124/,  &
       bvpd,  svpd/  4850.0, -0.0051/,  &
       bswp,  sswp/ -1.07e6, 7.42e-6/
!     save
!---------srf-05052006---------------------------
integer, intent(in) :: iroot
real, dimension(mzg), intent(out) ::fswpk
real, intent(out) :: fswp_equiv
!---------srf-05052006---------------------------

!if (CATT==1) then
!srf - define the heat capacity of the vegetation (hcapveg)
   hcapveg  = 3.e4                    !para cerrado, pastagem,...
   if (nveg<=7) hcapveg = 1.e5        !para floresta
   dtlcohcv = dtlc/hcapveg
!endif
!srf -------------------------------------------------------------

! Compute ground-canopy resistance rd.  Assume zognd not affected by snow.
! Assume (zoveg,zdisp) decrease linearly with snow depth, attaining
! the values (zognd,0) when veg covered.

zognd = soil_rough
zoveg = veg_rough * (1.-snowfac) + zognd * snowfac
zdisp = veg_height * (1.-snowfac)
zveg = zdisp / 0.63
!bob      rasgnd = log(zts / zognd) * log((zdisp + zoveg) / zognd)
!bob     +      / (vonk * vonk * vels)

! The following value of ground-canopy resistance for the
! nonvegetated (bare soil or water) surface is from John Garratt.
! It is 5/ustar and replaces the one from old leaf.

rasgnd = 5. / ustar

if (veg_tai >= .1 .and. snowfac < .9) then

! If vegetation is sufficiently abundant and not covered by snow, compute
! heat and moisture fluxes from vegetation to canopy, and flux resistance
! from soil or snow to canopy.

   factv = log(zts / zoveg) / (vonk * vonk * vels)
   aux = exp(exar * (1. - (zdisp + zoveg) / zveg))
   rasveg = factv * zveg / (exar * (zveg - zdisp)) * (exp(exar) - aux)
   c2 = max(0.,min(1., 1.1 * veg_tai / covr))
   rd = rasgnd * (1. - c2) + rasveg * c2
   ! From RAMS 6.0
   wshed = 0.
   qwshed = 0.
   transp = 0.

else

! If the TAI is very small or if snow mostly covers the vegetation, bypass
! vegetation computations.  Set heat and moisture flux resistance rd between
! the "canopy" and snow or soil surface to its bare soil value.  Set shed
! precipitation heat and moisture to unintercepted values.

   wshed = pcpgl
   qwshed = qpcpgl
   transp = 0.
   rd = rasgnd

endif

! Compute sensible heat and moisture fluxes between top soil or snow surface
! and canopy air.  wflxgc [kg/m2/s] is the upward vapor flux from soil or snow
! evaporation and dewgnd is the mass of dew that forms on the snow/soil
! surface this timestep; both are defined as always positive or zero.

hflxgc = cp * dens * (tempk(mzg+ksn) - can_temp) / rd
wflx = dens * (ground_rsat - can_rvap) / rd
dewgndflx = max(0.,-wflx)
dewgnd = dewgndflx * dtll
if (ksn == 0) then
   wflxgc = max(0.,dens * (ground_rvap - can_rvap) / rd)
else
   wflxgc = max(0.,min(wflx,sfcwater_mass(ksn)/dtll))
endif

if (veg_tai >= .1 .and. snowfac < .9) then

! If vegetation is sufficiently abundant and not covered by snow, compute
! heat and moisture fluxes from vegetation to canopy, and flux resistance
! from soil or snow to canopy.

! Compute veg-canopy resistance rb.  Note that rb and rc are both defined
! WITHOUT the LAI factor; the LAI factor is included later in computing
! fluxes that involve rb and/or rc.

   ustar0 = max(.1,ustar)
   c4 = .01 * sqrt(ustar0 * c1)

   rb  = (1. + .5 * veg_tai) / (.01 * sqrt(ustar0 * c1))

! Soil water potential factor for stomatal control
!---------srf-2006-2009--------------------------- root profiles
IF(  iroot == 1) THEN
   fswp_equiv= 0.
   swp =  -200.0
   nveg = nint(leaf_class)

   do k = mzg,1,-1
      !if(root(nveg,k) < 0.02) cycle

      nsoil = nint(soil_text(k))
      slpotv = slpots(nsoil) * (slmsts(nsoil) / soil_water(k)) ** slbs(nsoil)

      swp = max(-200., slpotv)
      !- water factor k-dependent
      fswpk(k)= ( 1. + exp(-sswp * (swp * 9810.- bswp)) )
      !- equivalent water
      fswp_equiv=fswp_equiv+root(nveg,k)/fswpk(k)

   enddo
   !- equivalent fswp (paralel array)
   fswp_equiv=1./fswp_equiv

ELSE  !RAMS original way


   swp = -200.0
   nveg = nint(leaf_class)

   do k = kroot(nveg),mzg
      nsoil = nint(soil_text(k))
      slpotv = slpots(nsoil) * (slmsts(nsoil) / soil_water(k)) ** slbs(nsoil)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Bypass this computation, and replace with wtroot = 1 for root zone
!      wtroot = float(int(root(k) + .95))
!      wtroot = 1.
!      swp = swp + wtroot * (max(swp,slpotv) - swp)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (slpotv > swp) swp = slpotv
   enddo
   swp = swp * 9810.

ENDIF
!---------srf-2006-2009--------------------------- root profiles


! Begin canopy time-split iterations

   do iter_can = 1,niter_can

! Calculate the saturation vapor pressure at TVEG

      tveg = veg_temp
      tvegc = tveg - 273.15
      es = eslf(tvegc)
      rsatveg = .622 * es / (prss - es)

! Compute mixing ratio at leaf surface using previous rc

      rc = stom_resist
      rleaf = (rb * rsatveg + rc * can_rvap) / (rb + rc)
      vpd = max((es - rleaf * prss / (.622 + rleaf)),0.)

! Evaluate 5 environmental factors and new rc

      ftlo = 1. + exp(-stlo * (tveg - btlo))
      fthi = 1. + exp(-sthi * (tveg - bthi))
      frad = 1. + exp(-srad * (rshort - brad))

!---------srf-2006-2009--------------------------- root profiles
      IF( iroot == 1) THEN
        fswp = fswp_equiv
      ELSE ! RAMS orig way
        fswp = 1. + exp(-sswp * (swp - bswp))
      ENDIF
!---------srf-05052006---------------------------

      fvpd = 1. + exp(-svpd * (vpd - bvpd))

! 15-minute response time for stomatal conductance (must have dtlc <= 900.)

      !DSM rc_inf = ftlo * fthi * frad * fvpd * fswp * rcmin(nveg)
      rc_inf = ftlo * fthi
      rc_inf = rc_inf * frad
      rc_inf = rc_inf * fvpd
      rc_inf = rc_inf * fswp
      rc_inf = rc_inf * rcmin(nveg)
      rc = 1. / (1. / rc + .0011 * dtlc * (1. / rc_inf - 1. / rc))

! Limit maximum transpiration to be <= 400 w m-2 by increasing rc if necessary.

      transp_test = alvl * dens * veg_lai * (rsatveg - can_rvap) / (rb + rc)
      if (transp_test > 400.) then
         rc = (rb + rc) * transp_test * .0025 - rb
      endif

      stom_resist = rc

! Flux of heat and moisture from vegetation to canopy

      stai = veg_tai * (1. - snowfac)
      slai = veg_lai * (1. - snowfac)

      hflxvc = 2. * stai * cp * dens * (tveg - can_temp) / rb
      c3 = dens * (rsatveg - can_rvap)
      if (c3 >= 0.) then
         sigmaw = min(1.,(veg_water / (.2 * stai)) ** .66667)
         wflxvc = min(c3 * 2. * stai * sigmaw / rb,veg_water/dtlc)
         ! From RAMS 6.0
         !transp = max(0.,c3 * slai * (1.-sigmaw) / (rb + rc))
         if((rb + rc)<1.0E-38 .and. (rb+rc)>-1.0E-38) then
            transp0=0.0
         else
            transp0 = max(0.,c3 * slai * (1.-sigmaw) / (rb + rc))
         endif
      else
         wflxvc = c3 * 2. * stai / rb
         ! From RAMS 6.0
         !transp = 0.
         transp0 = 0.
      endif

! Update vegetation moisture from precipitation, vapor flux with canopy,
! and shedding of excess moisture.

      wtemp = veg_water - wflxvc * dtlc + pcpgc * vf
      ! From RAMS 6.0
      !wshed = max(0.,wtemp - 0.2 * stai)
      !veg_water = max(0.,wtemp - wshed)
      wshed0 = max(0.,wtemp - 0.2 * stai)
      wshed = wshed + wshed0
      transp = transp + transp0
      veg_water = max(0.,wtemp - wshed0)

! Update vegetation temperature from radiative fluxes and sensible and
! latent heat transfer with canopy air.  Heat capacity of the vegetation
! is defined as hcapveg; it requires further testing and needs to be
! increased for dense vegetation.  A value of hcapveg of 3.e4 corresponds
! to the heat capacity of a layer of water 7.5 mm thick.

      ! From RAMS 6.0
      !tvegc = veg_temp - 273.15 + dtlcohcv  &
      !   * (rshort_v + rlonga_v + rlonggs_v - rlongv_gs  &
      !   -  rlongv_a - hflxvc - (wflxvc + transp) * alvl)
      tvegc = veg_temp - 273.15 + dtlcohcv  &
         * (rshort_v + rlonga_v + rlonggs_v - rlongv_gs  &
         -  rlongv_a - hflxvc - (wflxvc + transp0) * alvl)

! Exchange heat between vegetation and precipitation in implicit scheme

      qwtot = qpcpgc * vf + hcapveg * tvegc
      call qwtk(qwtot,pcpgc * vf,hcapveg,tvegk,fracliqv)
      veg_temp = tvegk
      fac = 4186.
      if (fracliqv <= .0001) fac = 2093.
      ! From RAMS 6.0
      !qwshed = (fac * (tvegk - 273.15) + fracliqv * 334000.) * wshed
      !DSM  qwshed = qwshed + (fac * (tvegk - 273.15) + fracliqv * 334000.) * wshed0

      !{DSM
      aux1=fac * (tvegk - 273.15) + fracliqv * 334000.
      qwshed_tmp=aux1 * wshed0

      if (qwshed_tmp > 3.0e38) then
         qwshed=3.0e38
      else
         qwshed=qwshed  + qwshed_tmp
      endif
      !DSM}


! Update temperature and moisture of canopy.  hcapcan [J/m2/K] and
! wcapcan [kg_air/m2] are
! the heat and moisture capacities of the canopy, while c2 and c3 are deltat
! times their inverses, respectively.  They require further testing.

      can_temp = can_temp + dtlcohcc * (hflxgc + hflxvc  &
               + dens * cp * ustar * tstar * pis)

      ! From RAMS 6.0
      !can_rvap = can_rvap + dtlcowcc * (wflxgc - dewgndflx + wflxvc  &
      !         + transp + dens * ustar * rstar)
      can_rvap = can_rvap + dtlcowcc * (wflxgc - dewgndflx + wflxvc  &
           + transp0 + dens * ustar * rstar)

   enddo

else

   can_temp = can_temp + dtllohcc * (hflxgc + dens * cp * ustar * tstar * pis)
   can_rvap = can_rvap + dtllowcc * (wflxgc - dewgndflx + dens * ustar * rstar)

endif

return
end subroutine canopy

!*****************************************************************************

!subroutine stars(t,ng,mynum,i0,j0,i,j,ip  &
subroutine stars(ustar,tstar,rstar,ths,rvs,thetacan,can_rvap,zts,patch_rough  &
!kml drydep
   ,patch_area,vels,vels_pat,vonk,dtllohcc,dens,dtll                        &
   ,R_aer)
!kml drydep

!CATT

!--(DMK-CCATT-INI)-----------------------------------------------------------
use ccatt_start, only: CCATT
!--(DMK-CCATT-OLD)-----------------------------------------------------------
!use catt_start, only: CATT
!--(DMK-CCATT-FIM)-----------------------------------------------------------

implicit none

real :: ustar,tstar,rstar,ths,rvs,thetacan,can_rvap,zts,patch_rough  &
       ,patch_area,vels,vels_pat,vonk,dtllohcc,dens,dtll,t

!kml drydep
real :: R_aer
!kml drydep

!integer i,j,ip,ng,mynum,i0,j0

real :: b,csm,csh,d,a2,c1,ri,fm,fh,c2,cm,ch,c3
!real :: tstaro,d_theta,d_thetan,theta_new
real :: ustaro,d_vel,d_veln,vel_new,delz
integer :: ifix,ifixu

real, parameter :: ubmin = .25   & ! should use ubmin=1.0 for convec case
                  ,ustmin = .1     !            ubmin=0.1 for stable case

! Routine to compute Louis (1981) surface layer parameterization.

vels_pat = max(vels,ubmin)
b = 5.
csm = 7.5
csh = 5.
d = 5.

! a2 is the drag coefficient in neutral conditions, here same for h/m
! ri is the bulk richardson numer, eq. 3.45 in Garratt

a2 = (vonk / log(zts / patch_rough)) ** 2
c1 = a2 * vels_pat
ri = 9.8 * zts * (ths - thetacan)  &
   / (.5 * (ths + thetacan) * vels_pat * vels_pat)

if (ths - thetacan > 0.) then   ! STABLE CASE

   fm = 1. / (1. + (2 * b * ri / sqrt(1 + d * ri)))
   fh = 1. / (1. + (3 * b * ri * sqrt(1 + d * ri)))

else                            ! UNSTABLE CASE

   c2 = b * a2 * sqrt(zts / patch_rough * (abs(ri)))
   cm = csm * c2
   ch = csh * c2
   fm = (1. - 2 * b * ri / (1. + 2 * cm))
   fh = (1. - 3 * b * ri / (1. + 3 * ch))

endif

ustar = max(ustmin,sqrt(c1 * vels_pat * fm))
c3 = c1 * fh / ustar
rstar = c3 * (rvs - can_rvap)
tstar = c3 * (ths - thetacan)

if (CCATT==1) then
   !kml drydep
   R_aer = 1./(a2*vels_pat*fh)
   !kml drydep
endif

! Limit ustar so that the flux cannot take more than 1/2 velocity in a timestep

ifixu=0
ustaro=ustar
delz = 2.*zts
d_vel =  - ustar * ustar *dtll / delz
vel_new = vels_pat + d_vel
if(vel_new < .5 * vels_pat) then
      ifixu=1
      d_veln = .5 * vels_pat
      ustar=sqrt(d_veln*delz/dtll)
endif

! Limit tstar and rstar so that the direction of the gradients cannot change
!   sign due to the fluxes
   !   can_thetan = thetacan + dtlcohcc * dens * cp * ustar * tstar * pis
!ifix=0
!tstaro=tstar
!delz = 2.*zts
!d_theta = dtllohcc * dens * cp * ustar * tstar
!theta_new = thetacan + d_theta
!if(ths - thetacan > 0.) then
!   if(theta_new > ths) then
!      ifix=1
!      d_thetan = ths - thetacan
!      tstar=d_thetan/(dtllohcc * dens * cp * ustar)
!   endif
!else
!   if(theta_new < ths) then
!      ifix=2
!      d_thetan = ths - thetacan
!      tstar=d_thetan/(dtllohcc * dens * cp * ustar)
!   endif
!endif
!if( (ng == 1 .and. i0+i == 53 .and. j0+j == 32) .or. ifix /= 0 .or. ifixu /= 0) then
!if( (ng == 1 .and. i0+i == 53 .and. j0+j == 32) ) then
!write(11,'(a7,f8.0,6i3,2f8.0,9f10.3)')  &
!    'leaf--f',t,ng,i0+i,j0+j,ip,ifixu,ifix,-1004.*1.15*ustar*tstar,-2.5e6*1.15*ustar*rstar  &
!   ,ustar,ustaro,d_vel,d_veln,vels_pat,vel_new,patch_area,patch_rough
   !,ustar,tstar,tstaro,ths,thetacan,theta_new,d_theta,d_thetan,vels_pat,patch_area,patch_rough
!endif

!if(ng == 1 .and. i0+i == 53 .and. j0+j == 32)  &
!if(abs(-1004.*1.15*ustar*tstar) > 1000. .or. abs(-2.5e6*1.15*ustar*rstar) > 1000.)  &
!write(11,'(a7,f8.0,4i3,2f8.0,2f8.3,f10.5,f9.4,2f7.1,2f9.4,2f7.2,2f8.4,f8.3)')  &
!    'leaf--3',t,ng,i0+i,j0+j,ip,-1004.*1.15*ustar*tstar,-2.5e6*1.15*ustar*rstar  &
!   ,ustar,tstar,rstar,ri,ths,thetacan,rvs,can_rvap,vels_pat,patch_rough,c1,fh,patch_area

return
end subroutine stars

!*****************************************************************************

subroutine sfclmcv(i,j,ustar, tstar, rstar, vels, vels_pat, ups, vps, gzotheta, &
     patch_area, sflux_u, sflux_v, sflux_w, sflux_t, sflux_r)
     !
  !  +---------------------------------------------------------------+
  !  \  This routine computes the turbulent fluxes of momentum,      \
  !  \  heat and moisture from the surface layer using the           \
  !  \  Manton-Cotton algebraic surface layer equations.             \
  !  +---------------------------------------------------------------+

  use rconstants

  ! TEB_SPM
  use teb_spm_start, only: TEB_SPM !INTENT(IN)

  implicit none
  ! Arguments:
  integer, intent(in) :: i,j
  real, intent(in)    :: ustar, tstar, rstar, vels,vels_pat, ups, vps, &
       gzotheta, patch_area
  real, intent(inout) :: sflux_u, sflux_v, sflux_w, sflux_t, sflux_r
  ! For TEB
  !real, pointer, optional :: G_URBAN
! Local Variables:
  real :: zoverl
  real :: wtol, cosine1, sine1, vtscr, cx, psin

  data wtol/1e-20/

  cosine1 = ups/vels_pat
  sine1   = vps/vels_pat

  vtscr   = ustar*patch_area


  sflux_u = sflux_u - ustar*vtscr*cosine1
  sflux_v = sflux_v - ustar*vtscr*sine1
  sflux_t = sflux_t - tstar*vtscr
  sflux_r = sflux_r - rstar*vtscr

  zoverl = gzotheta*vonk*tstar/(ustar*ustar)

  if (zoverl<0.)then
     cx = zoverl*sqrt(sqrt(1. - 15.*zoverl))
  else
     cx = zoverl/(1.0 + 4.7*zoverl)
  endif

  psin    = sqrt((1.-2.86*cx)/(1. + cx*(-5.39 + cx*6.998)))
  sflux_w = sflux_w + &
       (0.27*max(6.25*(1. - cx)*psin, wtol) - 1.18*cx*psin)*ustar*vtscr

end subroutine sfclmcv

!*****************************************************************************

subroutine grndvap(soil_energy,soil_water,soil_text,sfcwater_energy  &
   ,sfcwater_nlev,ground_rsat,ground_rvap,can_temp,can_rvap,prsg)

use leaf_coms
use rconstants

implicit none

real :: soil_energy,soil_water,soil_text,sfcwater_energy
real :: sfcwater_nlev,ground_rsat,ground_rvap,can_temp,can_rvap,prsg

integer :: ksn,nsoil

real :: gdrm,tempkk,fracliqq,slpotvn,alpha,beta
real :: rslif
real :: auxVar


!write (88,fmt='(A,7(E18.6,1X))') 'slbs(1:7)= ',slbs(1:7)

gdrm = g / rm

ksn = nint(sfcwater_nlev)

! ground_rsat(i,j) is the saturation mixing ratio of the top soil/snow surface
! and is used for dew formation and snow evaporation.

if (ksn > 0) then
   call qtk(sfcwater_energy,tempkk,fracliqq)
   ground_rsat = rslif(prsg,tempkk)
else

! Without snowcover, ground_rvap is the effective saturation mixing
! ratio of soil and is used for soil evaporation.  First, compute the
! "alpha" term or soil "relative humidity" and the "beta" term.

   nsoil = nint(soil_text)

   call qwtk(soil_energy,soil_water*1.e3,slcpd(nsoil),tempkk,fracliqq)
   ground_rsat = rslif(prsg,tempkk)

    !write(*,fmt='("nsoil: ",I3.3)') nsoil
    !write (*,fmt='(A,E18.6)') 'slpots:',slpots(nsoil)
    !write (*,fmt='(A,E18.6)') 'slmsts:',slmsts(nsoil)
    !write (*,fmt='(A,E18.6)') 'soil_water:',soil_water
    !write (*,fmt='(A,E18.6)') 'slbs(nsoil)=',slbs(nsoil)

   slpotvn = slpots(nsoil) * (slmsts(nsoil) / soil_water) ** slbs(nsoil)
   !write (*,fmt='(A,E18.6)') 'grdm:',gdrm
   !write (*,fmt='(A,E18.6)') 'slpotvn:',slpotvn
   !write (*,fmt='(A,E18.6)') 'tempkk:',tempkk
   auxVar=gdrm * slpotvn / tempkk
   if(auxVar<-87.0) auxVar=-87.0
   !write (*,fmt='(A,E18.6)') 'Produto:',auxVar
   call flush(6)
   alpha = exp(auxVar)
   beta = .25 * (1. - cos (min(1.,soil_water / sfldcap(nsoil)) * 3.14159)) ** 2
   ground_rvap = ground_rsat * alpha * beta + (1. - beta) * can_rvap

endif

return
end subroutine grndvap

!*****************************************************************************

subroutine leaf_bcond(m2,m3,mzg,mzs,npat,jdim  &

   ,soil_water       ,sfcwater_mass  ,soil_energy     &
   ,sfcwater_energy  ,soil_text      ,sfcwater_depth  &
   ,ustar            ,tstar          ,rstar           &
   ,veg_albedo       ,veg_fracarea   ,veg_lai         &
   ,veg_tai                                           &
   ,veg_rough        ,veg_height     ,patch_area      &
   ,patch_rough      ,patch_wetind   ,leaf_class      &
   ,soil_rough       ,sfcwater_nlev  ,stom_resist     &
   ,ground_rsat      ,ground_rvap    ,veg_water       &
   ,veg_temp         ,can_rvap       ,can_temp        &
   ,veg_ndvip        ,veg_ndvic      ,veg_ndvif)

implicit none

integer :: m2,m3,mzg,mzs,npat,jdim

real, dimension(mzg,m2,m3,npat) :: soil_water,soil_energy,soil_text
real, dimension(mzs,m2,m3,npat) :: sfcwater_mass,sfcwater_energy            &
                                  ,sfcwater_depth
real, dimension(m2,m3,npat)     :: ustar        ,tstar         ,rstar        &
                                  ,veg_albedo   ,veg_fracarea  ,veg_lai      &
                                  ,veg_tai                                   &
                                  ,veg_rough    ,veg_height    ,patch_area   &
                                  ,patch_rough  ,patch_wetind  ,leaf_class   &
                                  ,soil_rough   ,sfcwater_nlev ,stom_resist  &
                                  ,ground_rsat  ,ground_rvap   ,veg_water    &
                                  ,veg_temp     ,can_rvap      ,can_temp     &
                                  ,veg_ndvip    ,veg_ndvic     ,veg_ndvif

integer :: i,j,k,ipat

do ipat = 1,npat
   do j = 1,m3

      ustar          (1,j,ipat) = ustar            (2,j,ipat)
      tstar          (1,j,ipat) = tstar            (2,j,ipat)
      rstar          (1,j,ipat) = rstar            (2,j,ipat)
      veg_fracarea   (1,j,ipat) = veg_fracarea     (2,j,ipat)
      veg_lai        (1,j,ipat) = veg_lai          (2,j,ipat)
      veg_tai        (1,j,ipat) = veg_tai          (2,j,ipat)
      veg_rough      (1,j,ipat) = veg_rough        (2,j,ipat)
      veg_height     (1,j,ipat) = veg_height       (2,j,ipat)
      patch_area     (1,j,ipat) = patch_area       (2,j,ipat)
      patch_rough    (1,j,ipat) = patch_rough      (2,j,ipat)
      patch_wetind   (1,j,ipat) = patch_wetind     (2,j,ipat)
      leaf_class     (1,j,ipat) = leaf_class       (2,j,ipat)
      soil_rough     (1,j,ipat) = soil_rough       (2,j,ipat)
      sfcwater_nlev  (1,j,ipat) = sfcwater_nlev    (2,j,ipat)
      stom_resist    (1,j,ipat) = stom_resist      (2,j,ipat)
      ground_rsat    (1,j,ipat) = ground_rsat      (2,j,ipat)
      ground_rvap    (1,j,ipat) = ground_rvap      (2,j,ipat)
      veg_water      (1,j,ipat) = veg_water        (2,j,ipat)
      veg_temp       (1,j,ipat) = veg_temp         (2,j,ipat)
      can_rvap       (1,j,ipat) = can_rvap         (2,j,ipat)
      can_temp       (1,j,ipat) = can_temp         (2,j,ipat)
      veg_ndvip      (1,j,ipat) = veg_ndvip        (2,j,ipat)
      veg_ndvic      (1,j,ipat) = veg_ndvic        (2,j,ipat)
      veg_ndvif      (1,j,ipat) = veg_ndvif        (2,j,ipat)

      ustar         (m2,j,ipat) = ustar         (m2-1,j,ipat)
      tstar         (m2,j,ipat) = tstar         (m2-1,j,ipat)
      rstar         (m2,j,ipat) = rstar         (m2-1,j,ipat)
      veg_albedo    (m2,j,ipat) = veg_albedo    (m2-1,j,ipat)
      veg_fracarea  (m2,j,ipat) = veg_fracarea  (m2-1,j,ipat)
      veg_lai       (m2,j,ipat) = veg_lai       (m2-1,j,ipat)
      veg_tai       (m2,j,ipat) = veg_tai       (m2-1,j,ipat)
      veg_rough     (m2,j,ipat) = veg_rough     (m2-1,j,ipat)
      veg_height    (m2,j,ipat) = veg_height    (m2-1,j,ipat)
      patch_area    (m2,j,ipat) = patch_area    (m2-1,j,ipat)
      patch_rough   (m2,j,ipat) = patch_rough   (m2-1,j,ipat)
      patch_wetind  (m2,j,ipat) = patch_wetind  (m2-1,j,ipat)
      leaf_class    (m2,j,ipat) = leaf_class    (m2-1,j,ipat)
      soil_rough    (m2,j,ipat) = soil_rough    (m2-1,j,ipat)
      sfcwater_nlev (m2,j,ipat) = sfcwater_nlev (m2-1,j,ipat)
      stom_resist   (m2,j,ipat) = stom_resist   (m2-1,j,ipat)
      ground_rsat   (m2,j,ipat) = ground_rsat   (m2-1,j,ipat)
      ground_rvap   (m2,j,ipat) = ground_rvap   (m2-1,j,ipat)
      veg_water     (m2,j,ipat) = veg_water     (m2-1,j,ipat)
      veg_temp      (m2,j,ipat) = veg_temp      (m2-1,j,ipat)
      can_rvap      (m2,j,ipat) = can_rvap      (m2-1,j,ipat)
      can_temp      (m2,j,ipat) = can_temp      (m2-1,j,ipat)
      veg_ndvip     (m2,j,ipat) = veg_ndvip     (m2-1,j,ipat)
      veg_ndvic     (m2,j,ipat) = veg_ndvic     (m2-1,j,ipat)
      veg_ndvif     (m2,j,ipat) = veg_ndvif     (m2-1,j,ipat)

      do k = 1,mzg
         soil_water       (k,1,j,ipat) = soil_water         (k,2,j,ipat)
         soil_energy      (k,1,j,ipat) = soil_energy        (k,2,j,ipat)
         soil_text        (k,1,j,ipat) = soil_text          (k,2,j,ipat)

         soil_water      (k,m2,j,ipat) = soil_water      (k,m2-1,j,ipat)
         soil_energy     (k,m2,j,ipat) = soil_energy     (k,m2-1,j,ipat)
         soil_text       (k,m2,j,ipat) = soil_text       (k,m2-1,j,ipat)
      enddo

      do k = 1,mzs
         sfcwater_mass    (k,1,j,ipat) = sfcwater_mass      (k,2,j,ipat)
         sfcwater_energy  (k,1,j,ipat) = sfcwater_energy    (k,2,j,ipat)
         sfcwater_depth   (k,1,j,ipat) = sfcwater_depth     (k,2,j,ipat)

         sfcwater_mass   (k,m2,j,ipat) = sfcwater_mass   (k,m2-1,j,ipat)
         sfcwater_energy (k,m2,j,ipat) = sfcwater_energy (k,m2-1,j,ipat)
         sfcwater_depth  (k,m2,j,ipat) = sfcwater_depth  (k,m2-1,j,ipat)
      enddo

   enddo

   if (jdim == 1) then

      do i = 1,m2
         ustar          (i,1,ipat) = ustar            (i,2,ipat)
         tstar          (i,1,ipat) = tstar            (i,2,ipat)
         rstar          (i,1,ipat) = rstar            (i,2,ipat)
         veg_albedo     (i,1,ipat) = veg_albedo       (i,2,ipat)
         veg_fracarea   (i,1,ipat) = veg_fracarea     (i,2,ipat)
         veg_lai        (i,1,ipat) = veg_lai          (i,2,ipat)
         veg_tai        (i,1,ipat) = veg_tai          (i,2,ipat)
         veg_rough      (i,1,ipat) = veg_rough        (i,2,ipat)
         veg_height     (i,1,ipat) = veg_height       (i,2,ipat)
         patch_area     (i,1,ipat) = patch_area       (i,2,ipat)
         patch_rough    (i,1,ipat) = patch_rough      (i,2,ipat)
         patch_wetind   (i,1,ipat) = patch_wetind     (i,2,ipat)
         leaf_class     (i,1,ipat) = leaf_class       (i,2,ipat)
         soil_rough     (i,1,ipat) = soil_rough       (i,2,ipat)
         sfcwater_nlev  (i,1,ipat) = sfcwater_nlev    (i,2,ipat)
         stom_resist    (i,1,ipat) = stom_resist      (i,2,ipat)
         ground_rsat    (i,1,ipat) = ground_rsat      (i,2,ipat)
         ground_rvap    (i,1,ipat) = ground_rvap      (i,2,ipat)
         veg_water      (i,1,ipat) = veg_water        (i,2,ipat)
         veg_temp       (i,1,ipat) = veg_temp         (i,2,ipat)
         can_rvap       (i,1,ipat) = can_rvap         (i,2,ipat)
         can_temp       (i,1,ipat) = can_temp         (i,2,ipat)
         veg_ndvip      (i,1,ipat) = veg_ndvip        (i,2,ipat)
         veg_ndvic      (i,1,ipat) = veg_ndvic        (i,2,ipat)
         veg_ndvif      (i,1,ipat) = veg_ndvif        (i,2,ipat)

         ustar         (i,m3,ipat) = ustar         (i,m3-1,ipat)
         tstar         (i,m3,ipat) = tstar         (i,m3-1,ipat)
         rstar         (i,m3,ipat) = rstar         (i,m3-1,ipat)
         veg_albedo    (i,m3,ipat) = veg_albedo    (i,m3-1,ipat)
         veg_fracarea  (i,m3,ipat) = veg_fracarea  (i,m3-1,ipat)
         veg_lai       (i,m3,ipat) = veg_lai       (i,m3-1,ipat)
         veg_tai       (i,m3,ipat) = veg_tai       (i,m3-1,ipat)
         veg_rough     (i,m3,ipat) = veg_rough     (i,m3-1,ipat)
         veg_height    (i,m3,ipat) = veg_height    (i,m3-1,ipat)
         patch_area    (i,m3,ipat) = patch_area    (i,m3-1,ipat)
         patch_rough   (i,m3,ipat) = patch_rough   (i,m3-1,ipat)
         patch_wetind  (i,m3,ipat) = patch_wetind  (i,m3-1,ipat)
         leaf_class    (i,m3,ipat) = leaf_class    (i,m3-1,ipat)
         soil_rough    (i,m3,ipat) = soil_rough    (i,m3-1,ipat)
         sfcwater_nlev (i,m3,ipat) = sfcwater_nlev (i,m3-1,ipat)
         stom_resist   (i,m3,ipat) = stom_resist   (i,m3-1,ipat)
         ground_rsat   (i,m3,ipat) = ground_rsat   (i,m3-1,ipat)
         ground_rvap   (i,m3,ipat) = ground_rvap   (i,m3-1,ipat)
         veg_water     (i,m3,ipat) = veg_water     (i,m3-1,ipat)
         veg_temp      (i,m3,ipat) = veg_temp      (i,m3-1,ipat)
         can_rvap      (i,m3,ipat) = can_rvap      (i,m3-1,ipat)
         can_temp      (i,m3,ipat) = can_temp      (i,m3-1,ipat)
         veg_ndvip     (i,m3,ipat) = veg_ndvip     (i,m3-1,ipat)
         veg_ndvic     (i,m3,ipat) = veg_ndvic     (i,m3-1,ipat)
         veg_ndvif     (i,m3,ipat) = veg_ndvif     (i,m3-1,ipat)

         do k = 1,mzg
            soil_water       (k,i,1,ipat) = soil_water         (k,i,2,ipat)
            soil_energy      (k,i,1,ipat) = soil_energy        (k,i,2,ipat)
            soil_text        (k,i,1,ipat) = soil_text          (k,i,2,ipat)

            soil_water      (k,i,m3,ipat) = soil_water      (k,i,m3-1,ipat)
            soil_energy     (k,i,m3,ipat) = soil_energy     (k,i,m3-1,ipat)
            soil_text       (k,i,m3,ipat) = soil_text       (k,i,m3-1,ipat)
         enddo

         do k = 1,mzs
            sfcwater_mass    (k,i,1,ipat) = sfcwater_mass      (k,i,2,ipat)
            sfcwater_energy  (k,i,1,ipat) = sfcwater_energy    (k,i,2,ipat)
            sfcwater_depth   (k,i,1,ipat) = sfcwater_depth     (k,i,2,ipat)

            sfcwater_mass   (k,i,m3,ipat) = sfcwater_mass   (k,i,m3-1,ipat)
            sfcwater_energy (k,i,m3,ipat) = sfcwater_energy (k,i,m3-1,ipat)
            sfcwater_depth  (k,i,m3,ipat) = sfcwater_depth  (k,i,m3-1,ipat)
         enddo

      enddo

   endif

enddo
return
end subroutine leaf_bcond

!****************************************************************************

subroutine sfc_fields(m1,m2,m3,ia,iz,ja,jz,jd  &
   ,theta,rv,up,vp,dn0,pp,pi0,rtgt,zt,ths2,rvs2,ups2,vps2,pis2,dens2,zts2)

use leaf_coms
use rconstants

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,jd
real, dimension(m1,m2,m3) :: theta,rv,up,vp,dn0,pp,pi0
real, dimension(m2,m3) :: rtgt,ths2,rvs2,ups2,vps2,pis2,dens2,zts2
real, dimension(m1) :: zt

integer :: i,j
real :: hcpi

! Compute surface atmospheric conditions

hcpi = .5 * cpi

do j = ja,jz
   do i = ia,iz
      ths2(i,j) = theta(2,i,j)
      rvs2(i,j) = rv(2,i,j)
      ups2(i,j) = (up(2,i-1,j) + up(2,i,j)) * .5
      vps2(i,j) = (vp(2,i,j-jd) + vp(2,i,j)) * .5
      zts2(i,j) = zt(2) * rtgt(i,j)
      pis2(i,j) = (pp(1,i,j) + pi0(1,i,j) + pp(2,i,j) + pi0(2,i,j)) * hcpi
      dens2(i,j) = (dn0(1,i,j) + dn0(2,i,j)) * .5
   enddo
enddo

return
end subroutine sfc_fields

!****************************************************************************

subroutine sfc_fields_adap(m1,m2,m3,ia,iz,ja,jz,jd,lpu,lpv,lpw  &
   ,topma,aru,arv,theta,rv,up,vp,dn0,pp,pi0,zt,zm,dzt       &
   ,ths2,rvs2,ups2,vps2,pis2,dens2,zts2)

use leaf_coms
use rconstants

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,jd
integer, dimension(m2,m3) :: lpu,lpv,lpw
real, dimension(m1,m2,m3) :: aru,arv,theta,rv,up,vp,dn0,pp,pi0
real, dimension(m2,m3) :: topma,ths2,rvs2,ups2,vps2,pis2,dens2,zts2
real, dimension(m1) :: zt,zm,dzt

integer :: i,j,k1,k2,k3
real :: topma_t,wtw,wtu1,wtu2,wtv1,wtv2

! Compute surface atmospheric conditions

do j = ja,jz
   do i = ia,iz
      k2 = lpw(i,j)
      k1 = k2 - 1
      k3 = k2 + 1

      topma_t = .25 * (topma(i,j) + topma(i-1,j)  &
              + topma(i,j-jd) + topma(i-1,j-jd))

! weights for lowest predicted points, relative to points above them

      wtw = (zm(k2) - topma_t) * dzt(k2)
      wtu1 = aru(lpu(i-1,j),i-1,j)   / aru(lpu(i-1,j)+1,i-1,j)
      wtu2 = aru(lpu(i,j),i,j)       / aru(lpu(i,j)+1,i,j)
      wtv1 = arv(lpv(i,j-jd),i,j-jd) / arv(lpv(i,j-jd)+1,i,j-jd)
      wtv2 = arv(lpv(i,j),i,j)       / arv(lpv(i,j)+1,i,j)

      ths2(i,j) =  wtw * theta(k2,i,j) + (1. - wtw)  * theta(k3,i,j)

      rvs2(i,j) =  wtw * rv(k2,i,j)    + (1. - wtw)  * rv(k3,i,j)

      ups2(i,j) = (wtu1        * up(lpu(i-1,j),i-1,j)    &
                +  (1. - wtu1) * up(lpu(i-1,j)+1,i-1,j)  &
                +  wtu2        * up(lpu(i,j),i,j)        &
                +  (1. - wtu2) * up(lpu(i,j)+1,i,j)) * .5

      vps2(i,j) = (wtv1        * vp(lpv(i,j-jd),i,j-jd)    &
                +  (1. - wtv1) * vp(lpv(i,j-jd)+1,i,j-jd)  &
                +  wtv2        * vp(lpv(i,j),i,j)          &
                +  (1. - wtv2) * vp(lpv(i,j)+1,i,j)) * .5

      zts2(i,j) = (wtw * (zt(k2) - zm(k1))  &
                + (1. - wtw) * (zt(k3) - zm(k2)))

      if (wtw >= .5) then
         pis2(i,j)  = ((wtw - .5) * (pp(k1,i,j) + pi0(k1,i,j))  &
                    + (1.5 - wtw) * (pp(k2,i,j) + pi0(k2,i,j))) * cpi
         dens2(i,j) = (wtw - .5)  * dn0(k1,i,j)  &
                    + (1.5 - wtw) * dn0(k2,i,j)
      else
         pis2(i,j)  = ((wtw + .5) * (pp(k2,i,j) + pi0(k2,i,j))  &
                    + (.5 - wtw) * (pp(k3,i,j) + pi0(k3,i,j))) * cpi
         dens2(i,j) = (wtw + .5) * dn0(k2,i,j)  &
                    + (.5 - wtw) * dn0(k3,i,j)
      endif

   enddo
enddo

return
end subroutine sfc_fields_adap

!****************************************************************************

subroutine sfc_pcp(nqparm,level,i,j,cuparm,micro)

use mem_basic
use mem_micro
use mem_cuparm
use leaf_coms
use micphys, only: mcphys_type ! INTENT(IN)
use node_mod, only: mynum ! INTENT(IN)

implicit none

integer :: nqparm,level,i,j
type (cuparm_vars)  cuparm
type (micro_vars)   micro

if (nqparm>0) then

   pcpgl  = cuparm%conprr(i,j) * dtll
   qpcpgl = pcpgl  * 4186. * (ths * pis - 193.36)
   dpcpgl = pcpgl  * .001
   pcpgc  = dtlc_factor * pcpgl
   qpcpgc = dtlc_factor * qpcpgl

else

   pcpgl  = 0.
   qpcpgl = 0.
   dpcpgl = 0.
   pcpgc  = 0.
   qpcpgc = 0.

endif

if(mcphys_type == 2 .or. mcphys_type == 3 .or. mcphys_type == 4 .or. mcphys_type == 6 .or. mcphys_type == 7) then

   pcpgl  = pcpgl  + dtll_factor * micro%pcpg(i,j)
   qpcpgl = qpcpgl + dtll_factor * micro%pcpg(i,j)*4186. * (ths * pis - 193.36)
   dpcpgl = dpcpgl + dtll_factor * micro%pcpg(i,j) * .001
   pcpgc  = pcpgc  + dtlc_factor * dtll_factor * micro%pcpg(i,j)
   qpcpgc = qpcpgc + dtlc_factor * dtll_factor * micro%pcpg(i,j)*4186. * (ths * pis - 193.36)

else
  if (level >= 3) then
   pcpgl  = pcpgl  + dtll_factor * micro%pcpg(i,j)
   qpcpgl = qpcpgl + dtll_factor * micro%qpcpg(i,j)
   dpcpgl = dpcpgl + dtll_factor * micro%dpcpg(i,j)
   pcpgc  = pcpgc  + dtlc_factor * dtll_factor * micro%pcpg(i,j)
   qpcpgc = qpcpgc + dtlc_factor * dtll_factor * micro%qpcpg(i,j)
  endif
endif

return
end subroutine sfc_pcp

!****************************************************************************

subroutine vegndvi(ifm    &
   ,patch_area,leaf_class,veg_fracarea,veg_lai,veg_tai,veg_rough  &
   ,veg_height,veg_albedo,veg_ndvip,veg_ndvic,veg_ndvif)

use leaf_coms
use rconstants
use io_params
use mem_grid

!--(DMK-CCATT-INI)-----------------------------------------------------------
  use ccatt_start, only: &
       ccatt
!--(DMK-CCATT-FIM)-----------------------------------------------------------

!!$!CATT
!!$use catt_start, only: CATT  ! INTENT(IN)

implicit none

integer :: ifm
real :: patch_area,leaf_class,veg_fracarea,veg_lai,veg_tai,veg_rough  &
   ,veg_height,veg_albedo,veg_ndvip,veg_ndvic,veg_ndvif

integer :: nveg
integer, save :: nvcall = 0

real :: timefac_ndvi,sr,fpar,dead_lai,green_frac
real, save :: sr_min=1.081,fpar_min=.001,fpar_max=.950,fpcon=-.3338082
real, save :: ccc=-2.9657
real, save :: bz=.91,hz=.0075,extinc_veg=.5

real, dimension(nvtyp+nvtyp_teb), save :: dfpardsr

!  Initialize dfpardsr array

if (nvcall == 0) then
   nvcall = 1
   do nveg = 1,(nvtyp+nvtyp_teb)
      dfpardsr(nveg) = (fpar_max - fpar_min) / (sr_max(nveg) - sr_min)
   enddo
endif

!  Compute LAI, vegetation roughness, albedo, vegfrac from time-dependent NDVI

nveg = nint(leaf_class)

if (tai_max(nveg) < .1) then

   veg_lai = 0.
   veg_tai = 0.
   veg_rough = 0.
   veg_albedo = 0.
   veg_fracarea = 0.

else

   if (iupdndvi == 0) then
      timefac_ndvi = 0.
   else
      timefac_ndvi = (time - ndvitime1(ifm)) / (ndvitime2(ifm) - ndvitime1(ifm))
   endif

!  Time-interpolate ndvi to get current value veg_ndvic(i,j) for this patch

   veg_ndvic = veg_ndvip + (veg_ndvif - veg_ndvip) * timefac_ndvi

!  Limit ndvi to prevent values > .99 to prevent division by zero.

   if (veg_ndvic > .99) veg_ndvic = .99

! Compute "simple ratio" and limit between sr_min and sr_max(nveg).

   sr = (1. + veg_ndvic) / (1. - veg_ndvic)

   if (sr < sr_min) then
      sr = sr_min
   elseif (sr > sr_max(nveg)) then
      sr = sr_max(nveg)
   endif

! Compute fpar

   fpar = fpar_min + (sr - sr_min) * dfpardsr(nveg)

! Compute green leaf area index (veg_lai), dead leaf area index (dead_lai),
! total area index (tai), and green fraction

   veg_lai = glai_max(nveg) * (veg_clump(nveg) * fpar / fpar_max  &
           + (1. - veg_clump(nveg)) * alog(1. - fpar) * fpcon)

!--(DMK-CCATT-INI)-----------------------------------------------------------
!------srf------
    if (ccatt == 1) veg_lai = max(veg_lai,1.5)
!--(DMK-CCATT-FIM)-----------------------------------------------------------

   dead_lai = (glai_max(nveg) - veg_lai) * dead_frac(nveg)
   veg_tai = veg_lai + sai(nveg) + dead_lai
   green_frac = veg_lai / veg_tai

! Compute vegetation roughness height, albedo, and fractional area

   veg_rough = veg_height * (1. - bz * exp(-hz * veg_tai))
   veg_albedo = albv_green(nveg) * green_frac  &
              + albv_brown(nveg) * (1. - green_frac)
   veg_fracarea = veg_frac(nveg) * (1. - exp(-extinc_veg * veg_tai))

endif
return
end subroutine vegndvi

! Remaining issues:
!
! 1. Relationship between clumping, V, vegfrac
! 2. Impact of V on radiation
! 3. Put in 900 second response time for rc?
! 4. Build lookup tables, especially for things with exponentials?

!****************************************************************************

subroutine sfcrad(mzg, mzs, ip,                                            &
     soil_energy, soil_water, soil_text, sfcwater_energy, sfcwater_depth,  &
     patch_area, can_temp, veg_temp, leaf_class, veg_height, veg_fracarea, &
     veg_albedo, sfcwater_nlev, rshort, rlong, albedt, rlongup, cosz,      &
     ! For TEB
     G_URBAN, ETOWN, ALBTOWN, TSTOWN                                       &
     !
     )

  use mem_leaf
  use leaf_coms
  use rconstants
  use mem_scratch

  !TEB_SPM
  use teb_spm_start, only: TEB_SPM ! INTENT(IN)

  implicit none
  ! Arguments:
  integer, intent(IN) :: mzg, mzs, ip
  real, intent(IN)    :: soil_energy(mzg)
  real, intent(IN)    :: soil_water(mzg)
  real, intent(IN)    :: soil_text(mzg)
  real, intent(IN)    :: sfcwater_energy(mzs)
  real, intent(IN)    :: sfcwater_depth(mzs)
  real, intent(IN)    :: patch_area
  real, intent(IN)    :: can_temp
  real, intent(IN)    :: veg_temp
  real, intent(IN)    :: leaf_class
  real, intent(IN)    :: veg_height
  real, intent(IN)    :: veg_fracarea
  real, intent(IN)    :: veg_albedo
  real, intent(IN)    :: sfcwater_nlev
  real, intent(IN)    :: rshort
  real, intent(IN)    :: rlong
  real, intent(INOUT) :: albedt
  real, intent(INOUT) :: rlongup
  real, intent(IN)    :: cosz
  ! for TEB
  real, pointer, optional :: G_URBAN, ETOWN, ALBTOWN, TSTOWN
  !
  ! Local Variables:
  integer :: k, nsoil, nveg, ksn
  real :: alb, vfc, fcpct, alg, rad, als, fractrans, absg, algs, emv, emgs, &
       gslong, vlong, alv

  ! This routine is called by the radiation parameterization and by leaf.
  ! It computes net surface albedo plus radiative exchange between the
  ! atmosphere, vegetation, and the snow/ground given previously computed
  ! downward longwave and shortwave fluxes from the atmosphere.
  ! Also computed are functions of snowcover that are required for the above
  ! radiation calculations as well as other calculations in leaf.

  ! The shortwave parameterizations are only valid if the cosine of the zenith
  ! angle is greater than .03 .  Water albedo from Atwater and Bell (1981)

  ! alg, als, and alv are the albedos of the ground, snow, and vegetation
  ! (als needs a better formula based on age of the surface snow).

  ! absg and vctr32 are the actual fractions of shortwave incident on snow
  ! plus ground that get absorbed by the ground and each snow layer,
  ! respectively.  They currently use the variable fractrans, which is the
  ! fraction of light transmitted through each layer based on mass per square
  ! meter.  algs is the resultant albedo from snow plus ground.

  if (ip==1) then

     if (cosz>0.03) then
        alb    = min(max(-.0139 + .0467*tan(acos(cosz)), 0.03), 0.999)
        albedt = albedt + patch_area*alb
     endif

     call qtk(soil_energy(mzg), tempk(mzg), fracliq(mzg))
     rlongup = rlongup + patch_area*stefan*tempk(mzg)**4

  elseif (isfcl==0) then

     albedt  = albedt + patch_area*albedo
     rlongup = rlongup + patch_area*stefan*can_temp**4

  else

     ! Diagnose soil temperature and liquid fraction

     do k=1,mzg
        nsoil = nint(soil_text(k))
        call qwtk(soil_energy(k), soil_water(k)*1.e3,  &
             slcpd(nsoil), tempk(k), fracliq(k))
     enddo

     ! Diagnose snow temperature and the influence of snow covering veg.

     nveg    = nint(leaf_class)
     ksn     = nint(sfcwater_nlev)
     snowfac = 0.
     do k=1,ksn
        snowfac = snowfac + sfcwater_depth(k)
        call qtk(sfcwater_energy(k), tempk(k+mzg), fracliq(k+mzg))
     enddo
     snowfac = min(0.99, snowfac/max(0.001, veg_height))

     vf  = veg_fracarea*(1. - snowfac)
     vfc = 1. - vf

     ! Shortwave radiation calculations

     !srf-25-02-2005
     nsoil = nint(soil_text(mzg))
     !srf

     fcpct = soil_water(mzg)/slmsts(nsoil)
     if (fcpct>0.5) then
        alg = 0.14
     else
        alg = 0.31 - 0.34*fcpct
     endif
     alv = veg_albedo

     rad = 1.
     if (ksn>0) then

        ! als = .14 (the wet soil value) for all-liquid

        als = 0.5 - 0.36*fracliq(ksn + mzg)
        rad = 1.  - als
     endif
     do k=ksn,1,-1
        fractrans = exp(-20.*sfcwater_depth(k))
        vctr32(k) = rad*(1. - fractrans)
        rad = rad * fractrans
     enddo
     absg = (1. - alg)*rad
     algs = 1. - absg
     do k=ksn,1,-1
        algs        = algs - vctr32(k)
        rshort_s(k) = rshort*vfc*vctr32(k)
     enddo
     rshort_g = rshort*vfc*absg
     rshort_v = rshort*vf*(1. - alv + vfc*algs)
     !  rshort_a = rshort*(vf*alv + vfc*vfc*algs)

     alb = vf*alv + vfc*vfc*algs

     ! TEB
     if (TEB_SPM==1 .and. present(G_URBAN) .and. present(ALBTOWN)) then
        if (NINT(G_URBAN)==0) then
           albedt = albedt + patch_area*alb
        else
           albedt = albedt + patch_area*ALBTOWN
        endif
     else
        albedt = albedt + patch_area*alb
     endif

     ! Longwave radiation calculations

     emv  = emisv(nveg)
     emgs = emisg(nsoil)
     if (ksn>0) emgs = 1.0
     gslong = emgs*stefan*tempk(ksn+mzg)**4
     vlong  = emv*stefan*veg_temp**4

     rlonga_v  = rlong*vf*(emv + vfc*(1. - emgs))
     rlonga_gs = rlong*vfc*emgs
     rlongv_gs = vlong*vf*emgs
     rlongv_a  = vlong*vf*(2. - emgs - vf + emgs*vf)
     rlonggs_v = gslong*vf*emv
     rlonggs_a = gslong*vfc
     rlonga_a  = rlong*(vf*(1. - emv) + vfc*vfc*(1. - emgs))

     ! TEB
     if (TEB_SPM==1 .and. present(G_URBAN) .and. present(ETOWN) .and. &
          present(TSTOWN)) then
        IF (NINT(G_URBAN)==0) THEN
           rlongup = rlongup + patch_area*(rlongv_a + rlonggs_a + rlonga_a)
        else
           rlongup = rlongup + patch_area*ETOWN*STEFAN*TSTOWN**4
        ENDIF
     else
        rlongup = rlongup + patch_area*(rlongv_a + rlonggs_a + rlonga_a)
     endif

     !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     ! In case rlong is not computed, zero out all longwave fluxes other
     ! than rlongup.  [On the first timestep, radiative fluxes may not
     ! be available until microphysics is called, and zeroing these fluxes
     ! avoids the imbalance of having upward longwave without downward longwave
     ! in LEAF-2.  Also, this allows LEAF-2 to run without radiation for all
     ! timesteps, if desired for testing.]

     if (rlong<0.1) then
        rlonga_v  = 0.
        rlonga_gs = 0.
        rlongv_gs = 0.
        rlongv_a  = 0.
        rlonggs_v = 0.
        rlonggs_a = 0.
        rlonga_a  = 0.
     endif
     !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  endif

end subroutine sfcrad
