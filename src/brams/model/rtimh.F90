!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

module ModTimestep
contains
subroutine timestep(OneGrid,oneNamelistFile)

  use ModMessageSet, only: &
       PostRecvSendMsgs, &
       WaitRecvMsgs

  use ModAcoust, only: acoustic_new

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
       vveldamp

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
                             ,plume_fre_g      !

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

  use cuparm_grell3, only: cuparm_grell3_catt 

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

  USE wind_Farm, ONLY: wind_farm_driver,windfarm

  use optical, only: &
            aodDriver

  use ModNamelistFile, only : namelistFile

 ! use ModAdvectc_rk, only: &
 !        advectc_rk

  implicit none

  type(Grid), pointer :: OneGrid
  type(namelistFile), pointer :: oneNamelistFile

  ! execution time instrumentation
  include "tsNames.h"

  INTEGER, PARAMETER :: acoshdp = 0
  character(len=256) :: julesFile

  julesFile=oneNamelistFile%julesin
  !        +-------------------------------------------------------------+
  !        |   Timestep driver for the hybrid non-hydrostatic time-split |
  !        |      model.                                                 |
  !        +-------------------------------------------------------------+

  !  Zero out all tendency arrays.
  !--------------------------------
  call TEND0()

!!!!  IF( NNQPARM(ngrid) >=2 .OR. NNSHCU(ngrid) >=2 )CALL prepare_lsf_OLD(NNQPARM(ngrid), NNSHCU(ngrid),1)


  !  Thermodynamic diagnosis
  !--------------------------------
  if (mcphys_type <= 1 .and. level/=3) then
     call THERMO(mzp,mxp,myp,ia,iz,ja,jz,'SUPSAT')
  endif

  if (iexev == 2) &
     call exevolve(mzp,mxp,myp,ngrid,ia,iz,ja,jz,izu,jzv,jdim,mynum,dtlt,'ADV')

!		Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
!  call SynchronizedTimeStamp(TS_DYNAMICS)
  if (CCATT==1 .and. chemistry >= 1) call aodDriver(mzp,mxp,myp,ia,iz,ja,jz,ngrids)

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
     !
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
     if(AEROSOL==2) then
     !print*,"not doing matrix";call flush(6)
         CALL MatrixDriver(ia,iz,ja,jz,mzp,mxp,myp)
     endif

  endif

!!!  !srf- large and subgrid scale forcing for shallow and deep cumulus
!!!  IF( NNQPARM(ngrid) >=2 .OR. NNSHCU(ngrid) >=2 )CALL prepare_lsf_OLD(NNQPARM(ngrid), NNSHCU(ngrid),2)

!		Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
!  call SynchronizedTimeStamp(TS_PHYSICS)

  !  Send boundaries to adjoining nodes
  !-------------------------------------------
  if (nmachs > 1) then
     call PostRecvSendMsgs(OneGrid%SelectedGhostZoneSend, OneGrid%SelectedGhostZoneRecv)
  endif

  !  Coriolis terms
  !  ----------------------------------------
  call CORLOS(mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv, tend%ut, tend%vt)

  !  Velocity advection
  !----------------------------------------
  ! Use Optmized advection only in SX-6, for the moment
  if (machine==0) then
     ! If Generic IA32 use old Advction Scheme

     CALL ADVECTc   ('V',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
!    CALL advectc_rk('V',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)

  elseif (machine==1) then
     ! Using optmized advection scheme only in SX-6
     call calc_advec('V',ngrid,mzp,mxp,myp)
  endif

!		Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
!  call SynchronizedTimeStamp(TS_DYNAMICS)

  !  Cumulus parameterization version 1
  !----------------------------------------
  if (NNQPARM(ngrid)==1 .or. IF_CUINV==1) then
     call cuparm()
  end if

  !  Urban canopy parameterization
  !----------------------------------------
  if (IF_URBAN_CANOPY==1) call urban_canopy()

!		Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
!  call SynchronizedTimeStamp(TS_PHYSICS)

  !  Analysis nudging and boundary condition
  !------------------------------------------

  if (NUD_TYPE>0) call DATASSIM()


  !  Observation data assimilation
  !----------------------------------------
  if (IF_ODA==1) call oda_nudge()

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

  !  Sub-grid diffusion terms
  !----------------------------------------
  if ((if_adap==0) .and. (ihorgrad==2)) then
     call diffuse_brams31() !call optimized subroutine
  else
     call diffuse()
  endif

!!!!!  IF( NNQPARM(ngrid) >=2 .OR. NNSHCU(ngrid)>=2 ) CALL prepare_lsf_OLD(NNQPARM(ngrid), NNSHCU(ngrid),3)

  !  Velocity advection
  !----------------------------------------
  ! Use Optmized advection only in SX-6, for the moment
  if (machine==0) then

     IF(advmnt >= 1) THEN
      !-srf monotonic advection scheme
         call advmnt_driver('T',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
         if(advmnt >= 2) &
            CALL ADVECTc('T',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
     ELSE

         CALL ADVECTc   ('T',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
!        CALL advectc_rk('T',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)

     ENDIF     ! If Generic IA32 use old Advction Scheme
!----- for R-Ks time integration ---

  elseif (machine==1) then
     ! Using optmized advection scheme only in SX-6
     if (ngrid<=2) then
        call calc_advec('T',ngrid,mzp,mxp,myp)
     else
        call ADVECTc('T',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
     endif
  endif

!		Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
!  call SynchronizedTimeStamp(TS_DYNAMICS)

  !- STILT-BRAMS coupling (ML)
  if (imassflx == 1) call prep_advflx_to_stilt(mzp,mxp,myp,ia,iz,ja,jz,ngrid)

  !- large and subgrid scale forcing for shallow and deep cumulus
!!1  IF(  NNQPARM(ngrid) >=2 .OR. NNSHCU(ngrid)>=2 ) CALL prepare_lsf_OLD(NNQPARM(ngrid), NNSHCU(ngrid),4)
  IF( NNQPARM(ngrid) >=2 .OR. NNSHCU(ngrid)>=2 ) CALL prepare_lsf(NNQPARM(ngrid), NNSHCU(ngrid),1)

  !-   Cumulus parameterization options 2->6:
  !                    Deep Convection scheme
     !- call deep first, if there is deep convection , turn off shallow.
  IF(NNQPARM(ngrid)==2) CALL CUPARM_GRELL_CATT(1)
     !
     !                    Shallow Convection scheme
  IF(NNSHCU(ngrid)==2 ) CALL CUPARM_GRELL_CATT(2)
     !
     !- G3d - GD-FIM and GF
  IF(NNQPARM(ngrid)>=3) CALL CUPARM_GRELL3_CATT(OneGrid,1,NNQPARM(ngrid),NNSHCU(ngrid))

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

!		Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
!  call SynchronizedTimeStamp(TS_PHYSICS)

  !  Update scalars
  !----------------------------------------
  !
  call PREDTR()
  !
  !  Moisture variables positive definite
  !----------------------------------------
  if    (mcphys_type == 0) then
           call negadj1(mzp,mxp,myp)

  elseif(mcphys_type == 1) then
           call negadj1_2M_rams60(mzp,mxp,myp)

  endif

!		Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
!  call SynchronizedTimeStamp(TS_DYNAMICS)

  !  Microphysics
  !----------------------------------------
  if (mcphys_type == 0 .and. level==3) then
     if (machine==1 .and. TEB_SPM==0) then
        ! Optimized version only for SX-6
        call micro_opt()
     else
        ! Original Version used in a Generic IA32 machine
        call micro()
     endif
  endif
  if (mcphys_type == 1 .and. level==3) then
        ! 2M rams microphysics
        call micro_2M_rams60()
  endif
  if (mcphys_type == 2 .or. mcphys_type == 3 ) then
        ! G. Thompson microphysics
        call micro_thompson()
  endif
  if (mcphys_type == 4 ) then
        call micro_gfdl()
        
  elseif(mcphys_type >= 5 ) then
        call micro_wsm()
  
  endif

  !----------------------------------------
  !- chemistry - microphysics tranfers - sedimentation and tranfer from clouds to rain
  if (ccatt==1) then
     ! task 5 : sedimentation and mass transfer between clouds and rain
     call chemistry_driver(mzp,mxp,myp,ia,iz,ja,jz,5,50)
  endif

!		Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
!  call SynchronizedTimeStamp(TS_PHYSICS)

  !
  !  Thermodynamic diagnosis
  !----------------------------------------
  if (mcphys_type <= 1 .and. level==3)  then
     call THERMO(mzp,mxp,myp,1,mxp,1,myp,'MICRO')
  endif

  if (iexev == 2) call exevolve(mzp,mxp,myp,ngrid,ia,iz,ja,jz,izu,jzv,jdim,mynum,dtlt,'THS')

  !-damping on vertical velocity to keep stability
  if(vveldamp == 1) call w_damping(mzp,mxp,myp,ia,iz,ja,jz,mynum)

  !  Apply scalar b.c.'s
  !----------------------------------------
  call TRSETS()

  !  Lateral velocity boundaries - radiative
  !-------------------------------------------
  call LATBND()

  !  Apply Asselin time filter
  !---------------------------------------------------
  !
  !  Af(t)=A(t)+gama*(Af(t-deltat)-2*A(t)+A(t+deltat))
  !
  !  Where A denotes quantities and,
  !  Af is the filtered quantities
  !---------------------------------------------------

  !  First stage Asselin filter
  !------------------------------------------
  !
  !  scratch=A(t)+gama*(Af(t-deltat)-2*A(t))
  !
  !       +--------------+--------------+
  !       | IN           | OUT          |
  !  +----+--------------+--------------|
  !  | AP | Af(t-deltat) | Af(t-deltat) |
  !  |----+--------------+--------------|
  !  | AC | A(t)         | scratch      |
  !  +----+--------------+--------------+
  !
  !------------------------------------------
  call HADVANCE(1)
  !  Buoyancy term for w equation
  !----------------------------------------
  call BUOYANCY( tend%wt )

  !  Acoustic small timesteps
  !----------------------------------------
  ! AF: OUT: A(t+deltat)

  dts = 2. * dtlt / nnacoust(ngrid)

  call acoustic_new(OneGrid, nnacoust(ngrid),0 )   !MB:

  !  Last stage of Asselin filter
  !------------------------------------------
  !
  !  Af(t)=scratch+gama*A(t+deltat)
  !
  !       +--------------+--------------+
  !       | IN           | OUT          |
  !  +----+--------------+--------------|
  !  | AP ! A(t+deltat)  | Af(t)        |
  !  |----+--------------+--------------|
  !  | AC ! scratch      | A(t+deltat)  |
  !  +----+--------------+--------------+
  !
  !------------------------------------------
  call HADVANCE(2)

  !  Velocity/pressure boundary conditions
  !----------------------------------------
  call VPSETS()

  if (iexev == 2) call get_true_air_density(mzp,mxp,myp,ia,iz,ja,jz)

  ! Call THERMO on the boundaries
  call thermo_boundary_driver((time+dtlongn(ngrid)), dtlong, &
       f_thermo_e, f_thermo_w, f_thermo_s, f_thermo_n, &
       nzp, mxp, myp, jdim)

!		Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
!  call SynchronizedTimeStamp(TS_DYNAMICS)

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
    !if(AEROSOL==2) then
    !
    !   CALL MatrixDriver(ia,iz,ja,jz,mzp,mxp,myp)
    !endif
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

 !- apply digital filter
 if(applyDF) call applyDigitalFilter(fileNameDF, dfVars)

!		Uncoment to calculate execution time and set noInstrumentation = false in ModTimestamp.f90
! call SynchronizedTimeStamp(TS_PHYSICS)

end subroutine timestep
end module ModTimestep

!*************************************************************************

subroutine mass_flux(n1,n2,n3,m1,m2,m3,up,vp,wp  &
     ,dn0,rtgu,rtgv,dyu,dxv,pp,pi0)

  use mem_grid
  use rconstants

  implicit none
  integer :: n1,n2,n3,m1,m2,m3
  real :: up(m1,m2,m3),vp(m1,m2,m3),wp(m1,m2,m3)  &
       ,dn0(n1,n2,n3),rtgu(n2,n3),dyu(n2,n3),dxv(n2,n3)  &
       ,rtgv(n2,n3),pp(m1,m2,m3),pi0(n1,n2,n3)

  real, save :: aintmass=0.

  integer :: i,j,k
  real :: wmass,emass,smass,nmass,prtot,tmass,ppp,area

  !cc      if (mod(time,300.).gt..1) return

  !  west/east bound
  wmass=0.
  emass=0.
  do j=2,nyp-1
     do k=2,nzp-1
        i=1
        wmass=wmass +  &
             up(k,i,j)*rtgu(i,j)/(dyu(i,j)*dzt(k))  &
             *(dn0(k,i,j)+dn0(k,i+1,j))*.5
        i=nxp-1
        emass=emass -  &
             up(k,i,j)*rtgu(i,j)/(dyu(i,j)*dzt(k))  &
             *(dn0(k,i,j)+dn0(k,i+1,j))*.5
     enddo
  enddo

  !  north/south bound
  smass=0.
  nmass=0.
  do i=2,nxp-1
     do k=2,nzp-1
        j=1
        smass=smass +  &
             vp(k,i,j)*rtgv(i,j)/(dxv(i,j)*dzt(k))  &
             *(dn0(k,i,j)+dn0(k,i,j+1))*.5
        j=nyp-1
        nmass=nmass -  &
             vp(k,i,j)*rtgv(i,j)/(dxv(i,j)*dzt(k))  &
             *(dn0(k,i,j)+dn0(k,i,j+1))*.5
     enddo
  enddo

  k=2
  prtot=0.
  do j=2,nyp-1
     do i=2,nxp-1
        ppp= ( (pp(k,i,j)+pi0(k,i,j))/cp )**cpor*p00
        prtot=prtot+ppp/(dyu(i,j)*dxv(i,j))
     enddo
  enddo


  tmass=wmass+emass+smass+nmass
  aintmass=aintmass+tmass*dtlong
  area=(nxp-2)*deltax*(nyp-2)*deltay


  print*,'==============================='
  print*,' Mass flux - W, E, S, N'
  print*,  wmass,emass,smass,nmass
  print*, 'total (kg/(m2 s):',tmass/area
  print*, 'total (kg/m2):',aintmass/area
  print*, 'total pr change (pa):',aintmass/area*9.8
  print*, 'computed mean press:',prtot/area
  print*,'==============================='

  return
end subroutine mass_flux
!     *****************************************************************

subroutine w_damping(mzp,mxp,myp,ia,iz,ja,jz,mynum)

  use mem_basic, only: &
       basic_g           ! intent(in)

  use mem_grid, only: &
       ngrid,         & ! intent(in)
       grid_g           ! intent(in)
  use mem_tend, only: tend

    implicit none
    integer, intent(in) :: mzp
    integer, intent(in) :: mxp
    integer, intent(in) :: myp
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: mynum

  call apply_wdamp(mzp,mxp,myp,ia,iz,ja,jz,mynum   &
       ,basic_g(ngrid)%up   ,basic_g(ngrid)%vp     &
       ,basic_g(ngrid)%wp   ,grid_g(ngrid)%rtgt    &
       ,grid_g(ngrid)%f13t  ,grid_g(ngrid)%f23t    &
       ,grid_g(ngrid)%dxt   ,grid_g(ngrid)%dyt     &
       ,tend%ut             ,tend%vt               &
       ,tend%wt                                            )

end subroutine w_damping

! **********************************************************************

subroutine apply_wdamp(m1,m2,m3,ia,iz,ja,jz,mynum,up,vp,wp,rtgt,f13t,f23t,dxt,dyt,ut,vt,wt)

  use mem_grid, only: &
       jdim,          & ! intent(in)
       ngrid,         & ! intent(in)
       dtlt,          & ! intent(in)
       ht,            & ! intent(in)
       dzt              ! intent(in)

  implicit none
  integer, intent(in) :: m1
  integer, intent(in) :: m2
  integer, intent(in) :: m3
  integer, intent(in) :: ia
  integer, intent(in) :: iz
  integer, intent(in) :: ja
  integer, intent(in) :: jz
  integer, intent(in) :: mynum
  real,    intent(in) :: up(m1,m2,m3)
  real,    intent(in) :: vp(m1,m2,m3)
  real,    intent(in) :: wp(m1,m2,m3)
  real,    intent(in) :: rtgt(m2,m3)
  real,    intent(in) :: f13t(m2,m3)
  real,    intent(in) :: f23t(m2,m3)
  real,    intent(in) :: dxt(m2,m3)
  real,    intent(in) :: dyt(m2,m3)
  real,    intent(inout) :: ut(m1,m2,m3)
  real,    intent(inout) :: vt(m1,m2,m3)
  real,    intent(inout) :: wt(m1,m2,m3)

  integer :: i,j,k,ifm,icm,innest
  real :: c1x,c1y,c1z,cflnumh,cflnumv,cflz
  real :: vctr1(m1)
  real :: vctr2(m1)
  real :: vctr3(m1)
  real , parameter ::gama_w=0.3 !m/s²
  !     This routine damps the vertical velocity when CFLZ is exceed
  !     (Actually check on 80% of CFL)

  cflnumh = .80
  cflnumv = .50
  cflz=0.0
  !
  !c1x=0.0
  !vctr3(:)=0.0

  do j = ja,jz
    do i = ia,iz
      do k = 2,m1-1

           !cflx
	   !vctr1(k) = .5*(up(k,i,j)+up(k,i-1,j))*dtlt*dxt(i,j)
           !cfly
	   !vctr2(k) = .5*(vp(k,i,j)+vp(k,i,j-jdim))*dtlt*dyt(i,j)
           !cflz
	   vctr3(k) = ((wp(k,i,j)+wp(k-1,i,j))  &
                +(up(k,i,j)+up(k,i-1,j))*f13t(i,j)*ht(k)*rtgt(i,j)  &
                +(vp(k,i,j)+vp(k,i,j-jdim))*f23t(i,j)*ht(k)*rtgt(i,j)  &
                )*.5*dtlt*dzt(k)
	   c1z=abs(vctr3(k))
	   if(c1z > cflnumv) then
	     wt(k,i,j) = wt(k,i,j) -gama_w*sign(1.,wp(k,i,j))*(c1z-cflnumv)
	     print*,'wdamp applied at=',k,i,j,mynum,c1z,wp(k,i,j)
	     call flush(6)
	   endif
       enddo
       !do k = 2,m1-1
       !   c1x = abs(vctr1(k))
       !   c1y = abs(vctr2(k))
       !   c1z = abs(vctr3(k))
       !
       !   if (c1x .gt. cflxy(ngrid)) cflxy(ngrid) = c1x
       !   if (c1y .gt. cflxy(ngrid)) cflxy(ngrid) = c1y
       !    if (c1z .gt. cflz) cflz = c1z
       !enddo
     enddo
  enddo
  !print*,'at wdamp2-max cflz',mynum,cflz
  !call flush(6)

end subroutine apply_wdamp

! ***************************************************************************
