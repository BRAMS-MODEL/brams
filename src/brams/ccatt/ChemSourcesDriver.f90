MODULE ChemSourcesDriver

!#ifdef SIMPLE
!  USE aer1_list, ONLY :             &
!       aer_nspecies=>nspecies,      & ! (IN)
!       spc_aer_alloc=>spc_alloc,    & ! (IN)
!       spc_aer_name=>spc_name,      & ! (IN)
!       nmodes,                      & ! (IN)
!       aer_bburn => bburn,          & ! (IN)
!       aer_sdust => sdust,          & ! (IN)
!       aer_urban => urban,          & ! (IN)
!       aer_bioge => bioge,          & ! (IN)
!       aer_marin => marin,          & ! (IN)
!       aer_v_ash => v_ash,          & ! (IN)
!       matrix_level                   ! (IN)

!#elif MATRIX
!---srf incluir instrucoes para o MATRIX
 
     
  IMPLICIT NONE


CONTAINS

!=================================================================
  SUBROUTINE sources_driver(ngrid, mzp,mxp,myp,ia,iz,ja,jz,            &
                            g,cp,cpor,p00,rgas,pi180,           &
                            cosz,theta,pp,pi0,                  &
                            rv,dn0,up,vp,                       &
                            time,iyear1,imonth1,idate1,itime1,  &
                            dtlt,rtgt,lpw_r,glat,glon,zt,zm,dzt,  &
                            nzpmax,                             &
                            nvert_src,                          &
                            chem1_g,                            &
                            chem1_src_g,                        &
                            aer1_g,                             &
                            aer_nvert_src,                      &
                            plume_mean_g,                       &
                            dnp,                                &
                            emiss_cycle,                        &
                            aer2_g,                             &
                            plume_fre_g                         )
  USE chem1_list, ONLY :            &
       chem_nspecies =>nspecies,    & ! (IN)
       spc_chem_alloc=>spc_alloc,   & ! (IN)
       spc_chem_name =>spc_name,    & ! (IN)
       src,                         & ! (IN)
       on,                          & ! (IN)
       off,                         & ! (IN)
       CO2,                         & ! (IN)  !DSM
       SO2,                         & ! (IN)  !DSM
       transport                      ! (IN)

   USE aer1_list, ONLY :            &
       aer_nspecies=>nspecies,      & ! (IN)
       spc_aer_alloc=>spc_alloc,    & ! (IN)
       spc_aer_name=>aer_name,      & ! (IN)
       nmodes,                      & ! (IN)
       aer_bburn => bburn,          & ! (IN)
       aer_sdust => sdust,          & ! (IN)
       aer_urban => urban,          & ! (IN)
       aer_bioge => bioge,          & ! (IN)
       aer_marin => marin,          & ! (IN)
       aer_v_ash => v_ash,          & ! (IN)
       matrix_level                   ! (IN)

  USE mem_chem1, ONLY:              &
       chem1_vars,                  & ! Type
       chem1_src_vars,              & ! Type
       nsrc,                        & ! (IN)
       max_ntimes_src,              & ! (IN)
       nsrc,                        & ! (IN)
       bburn,                       & ! (IN)
       antro,                       & ! (IN)
       bioge,                       & ! (IN)
       geoge,                       & ! (IN)
       CHEMISTRY,                   & ! (IN)
       ntimes_src,                  & ! (IN)
       diur_cycle                     ! (IN)
  USE mem_aer1, ONLY:               &
       aerosol, &
       aer1_vars                 

  USE mem_leaf , ONLY: isfcl ! (IN)   !DSM

  USE mem_plume_chem1, ONLY:        &
       plume_vars,                  & ! Type
       plume_mean_vars,             & ! Type
       plume_fre_vars,             & ! Type
       plumerise,                   & ! (IN)
       prfrq,                       & ! (IN)
       nveg_agreg,                  & ! (IN)
       tropical_forest,             & ! (IN)
       boreal_forest,               & ! (IN)
       savannah,                    & ! (IN)
       grassland                      ! (IN)

  USE mem_stilt, ONLY:              &
       iexev                          ! (IN)

  USE mem_volc_chem1, ONLY:         &
       volcanoes                      ! (IN)

  USE chem_sources, ONLY:           &
       cycle_emission,              & ! Type
       get_diurnal_cycle_normalized,& ! Subroutine
       sources,                     & ! Subroutine
       emiss_cycle_time,            & ! (IN)
       srcmapfn                       ! (IN)

  USE mod_chem_plumerise_scalar, ONLY: &
       plumerise_driver               ! Subroutine

  USE io_params, ONLY:    &
       srctime1,          &
       srctime2

    INTEGER, INTENT(IN)  :: ngrid
    ! original
    INTEGER , INTENT(IN) :: mzp
    INTEGER , INTENT(IN) :: mxp
    INTEGER , INTENT(IN) :: myp
    INTEGER , INTENT(IN) :: ia
    INTEGER , INTENT(IN) :: iz
    INTEGER , INTENT(IN) :: ja
    INTEGER , INTENT(IN) :: jz

    ! rconstants
    REAL,  INTENT(IN) :: g
    REAL,  INTENT(IN) :: cp
    REAL,  INTENT(IN) :: cpor
    REAL,  INTENT(IN) :: p00
    REAL,  INTENT(IN) :: rgas
    REAL , INTENT(IN) :: pi180

    ! mem_radiate
    REAL , INTENT(IN) :: cosz(mxp,myp)

    ! mem_basic
    REAL,    INTENT(IN) :: theta(mzp,mxp,myp)
    REAL,    INTENT(IN) :: pp(mzp,mxp,myp)
    REAL,    INTENT(IN) :: pi0(mzp,mxp,myp)
    REAL,    INTENT(IN) :: rv(mzp,mxp,myp)
    REAL,    INTENT(IN) :: up(mzp,mxp,myp)!srf-AWE
    REAL,    INTENT(IN) :: vp(mzp,mxp,myp)!srf-AWE
    REAL,    POINTER    :: dn0(:,:,:)

    ! mem_grid
    REAL    , INTENT(IN) :: time
    INTEGER , INTENT(IN) :: iyear1
    INTEGER , INTENT(IN) :: imonth1
    INTEGER , INTENT(IN) :: idate1
    INTEGER , INTENT(IN) :: itime1
    REAL    , INTENT(IN) :: dtlt
    REAL    , INTENT(IN) :: rtgt(mxp,myp)
    real , INTENT(IN) :: lpw_r(mxp,myp)
    REAL    , INTENT(IN) :: glat(mxp,myp)
    REAL    , INTENT(IN) :: glon(mxp,myp)
    REAL    , INTENT(IN) :: zt(nzpmax)
    REAL    , INTENT(IN) :: zm(nzpmax)
    REAL    , INTENT(IN) :: dzt(nzpmax)

    ! grid_dims
    INTEGER, INTENT(IN) :: nzpmax

    ! mem_chem1
    INTEGER              , INTENT(IN)    :: nvert_src(nsrc)
    TYPE(chem1_vars)     , INTENT(INOUT) :: chem1_g(chem_nspecies)
    TYPE(chem1_src_vars) , INTENT(INOUT) :: chem1_src_g(2,nsrc,chem_nspecies)

    ! mem_aer1
    TYPE (aer1_vars) , INTENT(INOUT) :: aer1_g(nmodes,aer_nspecies),aer2_g(nmodes)
    INTEGER          , INTENT(IN)    :: aer_nvert_src(aer_nspecies)

    ! mem_plume_chem1
    TYPE (plume_mean_vars), INTENT(INOUT) :: plume_mean_g(:,:)
    TYPE (plume_fre_vars ), INTENT(INOUT) :: plume_fre_g (:,:)

    ! mem_stilt
    REAL    , POINTER    :: dnp(:,:,:) ! in

    TYPE(cycle_emission), INTENT(INOUT) :: emiss_cycle(nsrc)

    INTEGER :: lpw(mxp,myp)

    lpw=int(lpw_R)
    

    IF( CHEMISTRY < 0 .OR. srcmapfn(1:LEN_TRIM(srcmapfn)) == 'NONE' &
                      .OR. srcmapfn(1:LEN_TRIM(srcmapfn)) == 'none') RETURN


    IF(MOD(time+0.001,emiss_cycle_time) .LT. dtlt .OR. time .LT. .01)     & 
         CALL get_diurnal_cycle_normalized(mxp,myp,ia,iz,ja,jz,dtlt,glat,       &
                                           glon,imonth1,idate1,iyear1,itime1,   &
                                           antro,bioge,pi180,nsrc,emiss_cycle(:))

    !-  call plumerise and sources  
    IF( plumerise > 0 ) &
         CALL plumerise_driver(mzp,mxp,myp,ia,iz,ja,jz,                        &
                               srctime1,chem_nspecies,spc_chem_alloc,src,      &
                               on,off,nmodes,aer_nspecies,spc_aer_alloc,       &
                               aer_bburn,nsrc,bburn,nveg_agreg,tropical_forest,&
                               boreal_forest,savannah,grassland,nzpmax,dtlt,   &
                               time,zt,zm,dzt,                                 &
                               g,cp,cpor,p00,rgas,theta,pp,pi0,                &
                               rv,up,vp,rtgt,lpw,prfrq,aerosol,                &
                               chem1_src_g(:,:,:),                             &
                               aer1_g(:,:),                                    &
                               plume_mean_g(:,ngrid),                          &
                               spc_aer_name,                                   &
                               plume_fre_g(:,ngrid),                           &
                               plumerise                                       )

    CALL sources(mzp,mxp,myp,ia,iz,ja,jz,itime1,time,imonth1,idate1,iyear1,         &
                 glon,chem_nspecies,spc_chem_alloc,src,on,off,transport,nsrc,       &
                 nvert_src(:),chem1_src_g(:,:,:),chem1_g(:),bburn,bioge,antro,geoge,&
                 diur_cycle,aer_nvert_src(:),aer1_g(:,:),                           &
                 aerosol,aer_nspecies,spc_aer_alloc,nmodes,aer_bburn,aer_sdust,     &
                 aer_urban,aer_bioge,aer_marin,dnp,iexev,dn0,cosz,spc_chem_name,    &
                 emiss_cycle,volcanoes,aer_v_ash,CO2,isfcl,spc_aer_name,matrix_level,&
                 aer2_g(:),SO2)
  
  END SUBROUTINE  sources_driver
!
END MODULE ChemSourcesDriver
