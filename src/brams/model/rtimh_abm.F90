
!############################# Change Log ##################################
!    New file for the Adams-Bashforth-Moulton dynamical core
!    (after Wicker, Skamarock, 2002 and Wicker 2009 MWR)
!    Saulo Freitas (INPE), Michael Baldauf (DWD)
!
!############################# Change Log ##################################

module ModTimestep_ABM

  !MB: for testing only
  logical :: flag_mb_adv_test=.false.

contains
subroutine timestep_abm(OneGrid,oneNamelistFile)

  use ModMessageSet, only: &
       PostRecvSendMsgs, &
       WaitRecvMsgs

  use ModAcoust, only:         &
      acoustic_new,            &
      init_div_damping_coeff,  &
      deallocate_alpha_div,    &
      apply_div_damping

  use ModGrid, only: &
       Grid

  use mem_basic, only: &
       basic_g  ! INTENT(INOUT)

  use node_mod, only: &
       mzp, mxp, myp,  & ! INTENT(IN)
       ia, iz, ja, jz, & ! INTENT(IN)
       i0, j0,         & ! INTENT(IN)
       izu, jzv,       & ! INTENT(IN)
       mynum,          & ! INTENT(IN)
       ibcon,          & ! INTENT(IN)
       nmachs            ! INTENT(IN)

  use mem_cuparm, only: &
       NNQPARM, & ! INTENT(IN)
       IF_CUINV   ! INTENT(IN)

  use mem_varinit, only: &
       NUD_TYPE ! INTENT(IN)

  use mem_turb, only: &
       IF_URBAN_CANOPY, & ! INTENT(IN)
       ihorgrad           ! INTENT(IN)

  use mem_oda,   only: &
       if_oda ! INTENT(IN)

  use micphys,   only: &
       mcphys_type,  &! INTENT(IN)
       level          ! INTENT(IN)

  use mem_grid, only: &
       ngrids,     & ! INTENT(IN)
       ngrid,      & ! INTENT(IN)
       npatch,     & ! INTENT(IN)
       time,       & ! INTENT(IN)
       dts,        &
       dtlong,     & ! INTENT(IN)
       dtlongn,    & ! INTENT(IN)
       ideltat,    &
       nnacoust,   &
       iyear1,     & ! INTENT(IN)
       imonth1,    & ! INTENT(IN)
       idate1,     & ! INTENT(IN)
       grid_g,     & ! INTENT(INOUT)
       nxtnest,    & ! INTENT(IN)
       if_adap,    & ! INTENT(IN)
       dtlt,       & ! INTENT(IN)
       istp,       & ! INTENT(IN)
       jdim,       & ! INTENT(IN)
       nzp,        & ! INTENT(IN)
       f_thermo_e, & ! INTENT(IN)
       f_thermo_w, & ! INTENT(IN)
       f_thermo_s, & ! INTENT(IN)
       f_thermo_n, &   ! INTENT(IN)
       zt,         &
       zm,         &
       dzt,        &
       nzpmax,     &
       itime1,     &
       vveldamp    ! INTENT(IN)

  use shcu_vars_const, only: & ! For Shallow Cumulus Paramet.
       NNSHCU ! INTENT(IN)

  use mem_scalar, only: & ! For SiB
       scalar_g ! INTENT(IN)

  use mem_leaf, only: & ! For SiB
       ISFCL ! INTENT(IN)

  ! TEB_SPM
  use teb_spm_start, only: &
       TEB_SPM ! INTENT(IN)
  use mem_emiss, only: &
       ichemi,         & ! INTENT(IN)
       isource           ! INTENT(IN)

  !ALF- Necessary in new advection scheme
  use advect_kit, only : &
       calc_advec          ! Subroutine

  ! For specific optimization depending the type of machine
  use machine_arq, only: &
       machine ! INTENT(IN)

  use rconstants, only:  &
       g,                & ! (IN)
       cp,               & ! (IN)
       cpor,             & ! (IN)
       p00,              & ! (IN)
       rgas,             & ! (IN)
       pi180               ! (IN)

  use ccatt_start, only: &
       ccatt               ! (IN)

  use mem_tend, only: &
       tend

  use mem_chem1, only: &
       nvert_src=>chem1_src_z_dim_g, & ! (IN)
       chem1_g,                      & ! (INOUT)
       nsrc,                         & ! (IN)
       chem1_src_g,                  & ! %sc_src(INOUT)
       chemistry,                    & ! (IN)
       split_method,                 & ! (IN)
       n_dyn_chem,                   &
       ntimes_src

  use mem_aer1, only:                  &
       aerosol,                        &! (IN)
       aer1_g,                         &! %sc_src(INOUT)
       aer2_g,                         &! %sc_src(INOUT)
       aer_nvert_src=>aer1_src_z_dim_g  ! (IN)

  use mem_plume_chem1, only:  plume_mean_g &   ! %flam_frac(IN), %fire_size(IN)
                             ,plume_fre_g

  use mem_stilt, only: &
       iexev,          &  ! (IN)
       imassflx,       &  ! (IN)
       stilt_g            ! %dnp (IN)

  use mem_radiate, only: radiate_g

  use chem_sources, only :     &
       alloc_emiss_cycle,      &  ! Subroutine
       init_actual_time_index, &  ! Subroutine
       emiss_cycle,            &  ! (INOUT)
       emiss_cycle_alloc,      &
       srcmapfn                   ! (IN)

  use ChemSourcesDriver, only:  sources_driver            ! Subroutine

  use ChemDryDepDriver , only:  drydep_driver             ! Subroutine

  use module_chemistry_driver, only: chemistry_driver ! Subroutine

  use radiation, only: radiate ! Subroutine

  use ModTimeStamp, only: SynchronizedTimeStamp, TimeStamp

  use cuparm_grell3, only: cuparm_grell3_catt &  ! subroutine
                          ,g3d_g

  use digitalFilter, only: 	        &
                    applyDigitalFilter, & ! subroutine
                    fileNameDF,		& ! intent(inout) - file control
                    dfVars,             &
                    applyDF

  USE monotonic_adv, only:                 &
                           advmnt_driver,  &        ! subroutine
                           advmnt

  USE DriverMatrix, ONLY: MatrixDriver  !Matrix Aerosol Model


  USE rams_microphysics_2M, only: micro_2M_rams60,negadj1_2M_rams60

  use rrtm_driv, only:  rrtm_driver

  USE mem_radiate, ONLY: &
         ilwrtyp, iswrtyp

  use CUPARM_GRELL3, only: g3d_g

  use wind_farm, only: wind_farm_driver,windfarm

  use optical, only: &
            aodDriver

  use ModNamelistFile, only : namelistFile

  implicit none

  type(Grid), pointer :: OneGrid
  type(namelistFile), pointer :: oneNamelistFile

  ! execution time instrumentation
  include "tsNames.h"
  include "i8.h"

  !srf- for ABM3 scheme
  LOGICAL, parameter :: flag_Coriolis_in_every_RK_step = .FALSE.

  integer :: l_abm
  real    :: abm(2,3)
  integer :: abm_nmbr_small_timesteps(3)

  integer :: i,j,k
  integer :: n

  integer :: nmbr_gpts
  integer(kind=i8) :: mzxyp
  real,dimension(mzp*mxp*myp) :: vt3da,vt3db,vt3dc,vt3de,vt3dd
  character(len=256) :: julesFile

  !        +-------------------------------------------------------------------------+
  !        |   Timestep driver for the ABM3(Wicker, 2009) non-hydrostatic time-split |
  !        |      model.                                                             |
  !        +-------------------------------------------------------------------------+

  julesFile=oneNamelistFile%julesin

  mzxyp =mzp*mxp*myp

  !  Zero out all tendency arrays.
  !--------------------------------
  call TEND0()

  !  Thermodynamic diagnosis
  !--------------------------------
  if (mcphys_type <= 1 .and. level/=3) then
     call THERMO(mzp,mxp,myp,ia,iz,ja,jz,'SUPSAT')
  endif

  if (iexev == 2) &
     ! evolution of the Exner pressure: compression term
     call exevolve(mzp,mxp,myp,ngrid,ia,iz,ja,jz,izu,jzv,jdim,mynum,dtlt,'ADV')

  if (CCATT==1 .and. chemistry >= 1) call aodDriver(mxp,myp,mzp,ia,iz,ja,jz,ngrids)

  !  Radiation parameterization
  !--------------------------------
  call RADIATE(mzp,mxp,myp,ia,iz,ja,jz,mynum)

  !  Surface layer, soil and veggie model
  !----------------------------------------
  if (isfcl<=2) then
     call SFCLYR(mzp,mxp,myp,ia,iz,ja,jz,ibcon)
#ifdef JULES
  elseif (isfcl == 5) then
     if (time==0.) &
     call SFCLYR      (mzp,mxp,myp,ia,iz,ja,jz,ibcon)
     call SFCLYR_JULES(mzp,mxp,myp,ia,iz,ja,jz,jdim,julesFile)
  !DSM}
#endif
  endif

     !
     !-LFR Sea salt Aerossol inline source
     call SeaSaltDriver(ia,iz,ja,jz,ngrid,mxp,myp)

  !-- emission/deposition for CCATT chemistry models
  if (CCATT==1 .and. chemistry >= 0) then

     !- emission
     !- allocation for diurnal cycle of emission arrays   and
     !- actual_time_index on nodes
     if( (.not. emiss_cycle_alloc) .and. (chemistry >= 0) .and. &
          trim(srcmapfn) .ne.  'NONE' .and. trim(srcmapfn) .ne.  'none') then

	  call alloc_emiss_cycle(mxp,myp,ngrids,nsrc)

          call init_actual_time_index(nsrc,ntimes_src)

     endif


     !plume_mean_g(:,:) instead of plume_mean_g(:,ngrid) to avoid memory errors.
     !emiss_cycle(:,:)  instead of emiss_cycle(:,ngrid)  to avoid memory errors.
     !the same for the others var
     call sources_driver(ngrid, mzp,mxp,myp,ia,iz,ja,jz,                          &
                         g,cp,cpor,p00,rgas,pi180,                                &
                         radiate_g(ngrid)%cosz,basic_g(ngrid)%theta,              &
                         basic_g(ngrid)%pp,basic_g(ngrid)%pi0,basic_g(ngrid)%rv,  &
                         basic_g(ngrid)%dn0,basic_g(ngrid)%up,basic_g(ngrid)%vp,  &
	   		 time,iyear1,imonth1,idate1,itime1,dtlt,                  &
                         grid_g(ngrid)%rtgt,grid_g(ngrid)%lpw,grid_g(ngrid)%glat, &
                         grid_g(ngrid)%glon,zt,zm,dzt,nzpmax,                     &
                         nvert_src    (:,:),                                  &
                         chem1_g      (:,:),                                  &
                         chem1_src_g  (:,:,:,:),                              &
                         aer1_g       (:,:,:),                                &
                         aer_nvert_src(:,:),                                  &
                         plume_mean_g (:,:),                                  &
                         stilt_g(ngrid)%dnp,                                  &
                         emiss_cycle  (:,:),                                  &
                         aer2_g       (:,:),                                  &
			 plume_fre_g  (:,:)                                   )


     !- call dry deposition and sedimentation routines
     call drydep_driver(mzp,mxp,myp,ia,iz,ja,jz)

     !- call Matrix Aerosol Model
     !----------------------------------------
     !if(AEROSOL==2) then
     !print*,"not doing matrix";call flush(6)
     !    CALL MatrixDriver(ia,iz,ja,jz,mzp,mxp,myp)
     !endif

  endif

!!$  call TimeStamp(TS_PHYSICS)

  !  Send boundaries to adjoining nodes
  !-------------------------------------------
  if (nmachs > 1) then
     call PostRecvSendMsgs(OneGrid%SelectedGhostZoneSend, OneGrid%SelectedGhostZoneRecv)
  endif

  !  Coriolis terms
  !  ----------------------------------------
  if ( .NOT. flag_Coriolis_in_every_RK_step ) then
    call CORLOS(mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv, tend%ut, tend%vt)
  end if

  !  Cumulus parameterization version 1
  !----------------------------------------
  if (NNQPARM(ngrid)==1 .or. IF_CUINV==1) then
     call cuparm()
  end if

  !  Urban canopy parameterization
  !----------------------------------------
  if (IF_URBAN_CANOPY==1) call urban_canopy()

  !  Analysis nudging and boundary condition
  !------------------------------------------
!!$  call TimeStamp(TS_PHYSICS)
  if (NUD_TYPE>0) call DATASSIM()
!!$  call TimeStamp(TS_DINAMICS)

  !  Observation data assimilation
  !----------------------------------------
  if (IF_ODA==1) call oda_nudge()
!!$  call TimeStamp(TS_PHYSICS)

  !  Nested grid boundaries
  !----------------------------------------
  if (nxtnest(ngrid)>=1) call nstbdriv()

  !  Rayleigh friction for theta
  !----------------------------------------
  call RAYFT()

  !  Get the overlap region between parallel nodes
  !---------------------------------------------------
  if (nmachs > 1) then
     call WaitRecvMsgs(OneGrid%SelectedGhostZoneSend, OneGrid%SelectedGhostZoneRecv)
  endif

  if (iexev == 2) &
     call exevolve(mzp,mxp,myp,ngrid,ia,iz,ja,jz,izu,jzv,jdim,mynum,dtlt,'THA')

  !- task 2:  NO production by "eclair"
  if (ccatt == 1) &
     call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,2,50)

  !- CATT & Chemistry == CCATT
  !----------------------------------------
  if (ccatt==1 .and. split_method== 'PARALLEL' .and. n_dyn_chem==1) then
     ! task 3 : production/loss by chemical processes and inclusion of the
     ! chemistry tendency at the total tendency
     CALL chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,3,50)
  endif
  if (ccatt==1 ) then
     ! task 4 : mass transfer between gas and liquid
     call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,4,50)
  endif

  !---------------------------------------------------
  ! Shallow  cumulus parameterization by Souza
  if (NNSHCU(ngrid)==1) call SHCUPA()
  !---------------------------------------------------

  if (TEB_SPM==1) then
     ! Update urban emissions
     !----------------------------------------
     if (isource==1) then
        call sources_teb(mzp, mxp, myp, ia, iz, ja, jz, ngrid, ngrids)
     endif
     !  Update chemistry
     !----------------------------------------
     if (ichemi==1) then
        call ozone(mzp, mxp, myp, ia, iz, ja, jz, ngrid, dtlt)
     endif
  endif

  if (iexev == 2) &
  call exevolve(mzp,mxp,myp,ngrid,ia,iz,ja,jz,izu,jzv,jdim,mynum,dtlt,'THS')

  !  Sub-grid diffusion terms
  !----------------------------------------
  if ((if_adap==0) .and. (ihorgrad==2)) then
     call diffuse_brams31() !call optimized subroutine
  else
     call diffuse()
  endif

  !- STILT-BRAMS coupling (ML)
  if (imassflx == 1) call prep_advflx_to_stilt(mzp,mxp,myp,ia,iz,ja,jz,ngrid)

  !- large and subgrid scale forcing for shallow and deep cumulus
  IF( NNQPARM(ngrid) >=2 .OR. NNSHCU(ngrid)>=2 ) CALL prepare_lsf(NNQPARM(ngrid), NNSHCU(ngrid),1)

  !- cumulus parameterizations options: G3d - GD-FIM and GF
  IF(NNQPARM(ngrid)>=3) CALL CUPARM_GRELL3_CATT(OneGrid,1,NNQPARM(ngrid),NNSHCU(ngrid))

  !------------------------------------------------------------------------------
  ! init preparations for ABM3 -loop
  !------------------------------------------------------------------------------
  ! Wicker, 2009 ABM3 time-stepping scheme:
  ! 1)  q*      = q(t) +  1/2 dt(3tend[q(t)] -  tend[q(t-dt)] )
  ! 2)  q(t+dt) = q(t) + 1/12 dt(5tend[q*  ] + 8tend[q(t)   ] - tend[q(t-dt)])
  ! where tend = slow mode evolution operator
  abm(1,2) =  3.0 / 2.0
  abm(1,3) = -1.0 / 2.0

  abm_nmbr_small_timesteps(1) = nnacoust(ngrid)
  abm_nmbr_small_timesteps(2) = nnacoust(ngrid)

  abm(2,1) =  5.0 / 12.0
  abm(2,2) =  8.0 / 12.0
  abm(2,3) = -1.0 / 12.0

  dts = dtlt / nnacoust(ngrid)

  if ( apply_div_damping ) then
     ! calculation of divergence damping coeff:
     !MB: ATTENTION: dts different for different nests?!!
     if ( ideltat == 0 ) then
        ! it is sufficient to calculate alpha_div only once:
        if ( istp == 1 ) then
           call init_div_damping_coeff( dts )
        end if
     else
        call init_div_damping_coeff( dts )
     end if
  end if
  !  Lateral velocity boundaries - radiative
  !-------------------------------------------
  call LATBND()

  ! begin of ABM3 loop
  !---------------------------------------------------
  do l_abm = 1, 2

     !- initialize the tendencies with the physics tendencies
     tend%ut_rk (:) =tend%ut (:)
     tend%vt_rk (:) =tend%vt (:)
     tend%wt_rk (:) =tend%wt (:)
     tend%pt_rk (:) =tend%pt (:)
     tend%tht_rk(:) =tend%tht(:)

     !- advection should give back tendencies
     !- ut_rk, vt_rk, wt_rk, pt_rk, tht_rk = physics tend + advection tendency

     !  Velocity advection
     !----------------------------------------
     call advectc_rk('V',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)

     !  advection of pi and theta_il
     !----------------------------------------
     call advectc_rk('THETAIL',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
     call advectc_rk('PI'     ,mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)

     !  coriolis term
     !-----------------------------------------
     if ( flag_Coriolis_in_every_RK_step ) then
       !if ( .not. flag_mb_adv_test)
          call CORLOS(mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv, tend%ut_rk, tend%vt_rk)
       !end if
     endif

     !  Buoyancy term for w equation
     !----------------------------------------
     !if ( .not. flag_mb_adv_test)
        call BUOYANCY( tend%wt_rk )
     !end if

     !-- slow mode tendencies
     if(l_abm == 1) then
       !- save tendencies at time T_n in the scratchs vt3d(a,...,e) arrays for the 2nd step
       call copy_long_rk(mzxyp,tend%ut_rk ,vt3da)
       call copy_long_rk(mzxyp,tend%vt_rk ,vt3db)
       call copy_long_rk(mzxyp,tend%wt_rk ,vt3dc)
       call copy_long_rk(mzxyp,tend%pt_rk ,vt3dd)
       call copy_long_rk(mzxyp,tend%tht_rk,vt3de)

       if(time < dtlongn(1)) then
         !- only for the 1st model timestep
	 !- means:  tend%ut_past (:)=tend%ut_rk (:)
         call copy_long_rk(mzxyp,tend%ut_rk ,tend%ut_past )
         call copy_long_rk(mzxyp,tend%vt_rk ,tend%vt_past )
         call copy_long_rk(mzxyp,tend%wt_rk ,tend%wt_past )
         call copy_long_rk(mzxyp,tend%pt_rk ,tend%pt_past )
         call copy_long_rk(mzxyp,tend%tht_rk,tend%tht_past)

       endif
       !- get the ab2 tendencies for time integration together with acoustic sub-steps
       !- will provide u*,v*,w*,pi* and th* after the acoustic integration
       !- means = tend%ut_rk (:) = abm(1,2)*tend%ut_rk (:) + abm(1,3)*tend%ut_past (:)
       call get_ab2_tend(mzxyp,abm(1,2),tend%ut_rk ,abm(1,3) ,tend%ut_past ,tend%ut_rk)
       call get_ab2_tend(mzxyp,abm(1,2),tend%vt_rk ,abm(1,3) ,tend%vt_past ,tend%vt_rk)
       call get_ab2_tend(mzxyp,abm(1,2),tend%wt_rk ,abm(1,3) ,tend%wt_past ,tend%wt_rk)
       call get_ab2_tend(mzxyp,abm(1,2),tend%pt_rk ,abm(1,3) ,tend%pt_past ,tend%pt_rk)
       call get_ab2_tend(mzxyp,abm(1,2),tend%tht_rk,abm(1,3) ,tend%tht_past,tend%tht_rk)
    else
       !- get the am3 (final) tendencies for time integration together with acoustic sub-steps
       !- will provide u ,v ,w ,pi  and th at time T_n+1 after the acoustic integration
       !- means = tend%Xt_rk (:) =abm(2,1)*tend%Xt_rk (:) +abm(2,2)*v3tdx(:,:,:)+ abm(2,3)*tend%Xt_past (:)
       call get_am3_tend(mzxyp,abm(2,1),tend%ut_rk , abm(2,2),vt3da, abm(2,3) ,tend%ut_past  ,tend%ut_rk)
       call get_am3_tend(mzxyp,abm(2,1),tend%vt_rk , abm(2,2),vt3db, abm(2,3) ,tend%vt_past  ,tend%vt_rk)
       call get_am3_tend(mzxyp,abm(2,1),tend%wt_rk , abm(2,2),vt3dc, abm(2,3) ,tend%wt_past  ,tend%wt_rk)
       call get_am3_tend(mzxyp,abm(2,1),tend%pt_rk , abm(2,2),vt3dd, abm(2,3) ,tend%pt_past  ,tend%pt_rk)
       call get_am3_tend(mzxyp,abm(2,1),tend%tht_rk, abm(2,2),vt3de, abm(2,3) ,tend%tht_past ,tend%tht_rk)

    endif

    if ( l_abm == 2 ) then
       ! (not necessary in the first ABM substep)
       basic_g(ngrid)%uc (:,:,:) = basic_g(ngrid)%up (:,:,:)
       basic_g(ngrid)%vc (:,:,:) = basic_g(ngrid)%vp (:,:,:)
       basic_g(ngrid)%wc (:,:,:) = basic_g(ngrid)%wp (:,:,:)
       basic_g(ngrid)%pc (:,:,:) = basic_g(ngrid)%pp (:,:,:)
       basic_g(ngrid)%thc(:,:,:) = basic_g(ngrid)%thp(:,:,:)

       !- save tendencies for the next timestep
       call copy_long_rk(mzxyp,vt3da,tend%ut_past )
       call copy_long_rk(mzxyp,vt3db,tend%vt_past )
       call copy_long_rk(mzxyp,vt3dc,tend%wt_past )
       call copy_long_rk(mzxyp,vt3dd,tend%pt_past )
       call copy_long_rk(mzxyp,vt3de,tend%tht_past)

    endif

    !-  Acoustic small timesteps
    !- (here u,v,w and pi are actually moved forward from time t to time t+dt
    if ( .not. flag_mb_adv_test) then
       call acoustic_new(OneGrid, abm_nmbr_small_timesteps(l_abm),0 )
    else
      ! - just for testing
      ! without the acoustic subroutine, tendencies must be added "manually"
      ! update pp -> pc (similar for u, v, w, if needed.)
      call update_long_rk(int(mxp*myp*mzp,i8),dtlt,1.0 &
                        ,basic_g(ngrid)%pc,basic_g(ngrid)%pp  &
                        ,tend%pt_rk)
    endif

    !- update thp -> thc (theta_il is not contained in acoustic_new)
    call update_long_rk(int(mxp*myp*mzp,i8),dtlt,1.0 &
                        ,basic_g(ngrid)%thc,basic_g(ngrid)%thp  &
                        ,tend%tht_rk)
    !- determine theta (dry potential temp.) for the buoyancy term:
    call theta_thp_rk(mzp,mxp,myp,ia,iz,ja,jz,"get_theta")

    !-damping on vertical velocity to keep stability
    !MB: does this act on wc???
    if(vveldamp == 1) call w_damping(mzp,mxp,myp,ia,iz,ja,jz,mynum)

  end do
  !------------------------------------------------------------------------------
  ! end of ABM3 loop
  ! i.e. the fields
  !    basic_g%uc, ..%vc, ..%wc, ..%pc, ..%thc
  ! contain the fields at the new time level n+1
  !------------------------------------------------------------------------------

  if ( apply_div_damping ) then
     if ( ideltat /= 0 ) then
        call deallocate_alpha_div
     end if
  end if
  !
  !
  !  water species, tke and tracers advection
  !----------------------------------------
  IF(advmnt == 1) THEN
      !- monotonic advection scheme
      call advmnt_driver('SCALAR',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
  ELSE
  !>>>>
      CALL ADVECTc   ('SCALAR',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
  !   CALL advectc_rk('SCALAR',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
  !>>>>
  ENDIF

  !  Update scalars (water species, tke and tracers)
  !----------------------------------------
  !
   call PREDTR()

  !-  copy current time into past time (u,v,w,pi,thetail)
  !---> thp   must be changed to THC for microphysics/bc/theta update
  !---> pp    must be changed to PC  for microphysics
  !---> wp    must be changed to WC  for microphysics
  !---> up,vp must be changed to UC,VC for output
   basic_g(ngrid)%up (:,:,:) = basic_g(ngrid)%uc (:,:,:)
   basic_g(ngrid)%vp (:,:,:) = basic_g(ngrid)%vc (:,:,:)
   basic_g(ngrid)%wp (:,:,:) = basic_g(ngrid)%wc (:,:,:)
   basic_g(ngrid)%pp (:,:,:) = basic_g(ngrid)%pc (:,:,:)
   basic_g(ngrid)%thp(:,:,:) = basic_g(ngrid)%thc(:,:,:)
  !---->
  !---->

  !
  !  Moisture variables positive definite
  !----------------------------------------
  if     (mcphys_type == 0) then
           call negadj1(mzp,mxp,myp)

   elseif(mcphys_type == 1) then
           call negadj1_2M_rams60(mzp,mxp,myp)
  endif

  !  Microphysics (applied on THP, just updated)
  !----------------------------------------
  if (mcphys_type == 0 .and. level==3) then
     if (machine==1 .and. TEB_SPM==0) then
        !- optimized version only for SX-6
        call micro_opt()
     else
        !- original Version used in a Generic IA32 machine
        call micro()
     endif
  endif
  if (mcphys_type == 1 .and. level==3) then
        !- 2M rams microphysics
        call micro_2M_rams60()
  endif
  if (mcphys_type == 2 .or. mcphys_type == 3 ) then
        !- G. Thompson microphysics
        call micro_thompson()
  endif
  !----------------------------------------

  !- Thermodynamic diagnosis
  if (mcphys_type <= 1 .and. level==3)  then
     call THERMO(mzp,mxp,myp,1,mxp,1,myp,'MICRO')
  endif

  !  Apply scalar b.c.'s (THP is changed here)
  !----------------------------------------
  call TRSETS()
  !---> THC must be changed to THP to include microphysics/trsets changes
  !---> for the next timestep
        basic_g(ngrid)%thc(:,:,:) = basic_g(ngrid)%thp(:,:,:)
  !--->

  !  Lateral velocity boundaries - radiative
  !-------------------------------------------
  !call LATBND()

  !  Velocity/pressure boundary conditions
  !----------------------------------------
  call VPSETS()

  !- call THERMO on the boundaries
  call thermo_boundary_driver((time+dtlongn(ngrid)), dtlong, &
       f_thermo_e, f_thermo_w, f_thermo_s, f_thermo_n, &
       nzp, mxp, myp, jdim)

  if (iexev == 2) call get_true_air_density(mzp,mxp,myp,ia,iz,ja,jz)

  !----------------------------------------
  !- chemistry - microphysics tranfers - sedimentation and tranfer from clouds to rain
  if (ccatt==1) then
     ! task 5 : sedimentation and mass transfer between clouds and rain
     call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,5,50)
  endif
  !
  !----------------------------------------
  !- chemistry/aerosol solvers
  if (ccatt==1) THEN
    if ( (split_method== 'PARALLEL' .AND. N_DYN_CHEM > 1) .OR.  &
         (split_method== 'SEQUENTIAL'                   ) .OR.  &
 	 (split_method== 'SYMMETRIC'                    )       )THEN

       ! task 3 : production/loss by chemical processes and final updated
       !  	of each specie
       call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,3,50)
    endif

    !- call Matrix Aerosol Model
    !- using symmetric/sequential spliting operator
    if(AEROSOL==2) then

       CALL MatrixDriver(ia,iz,ja,jz,mzp,mxp,myp)
    endif
  endif
  if (ccatt==1 .and. aerosol == 1) THEN
    call aer_background(ngrid,mzp,mxp,myp,ia,iz,ja,jz)
  endif
  !----------------------------------------

  if (TEB_SPM==1) then
     !EDF  emission module
     if (isource==1) then
        ! Apply only for last finner grid
        if (ngrid==ngrids) then
           CALL le_fontes(ngrid, mzp, mxp, myp, &
                npatch, ia, iz, ja, jz, (time+dtlongn(1)))
        endif
     endif
     !EDF
  endif

  !windfarm
  call wind_farm_driver(ngrid,mzp,mxp,myp,ia,iz,ja,jz)

! >>>>srf -> check  latter for ABM3
 !- apply digital filter
 if(applyDF) call applyDigitalFilter(fileNameDF, dfVars)

!!$  call TimeStamp(TS_PHYSICS)

end subroutine timestep_abm


end module ModTimestep_ABM
