MODULE module_dry_dep


  USE chem1_list, ONLY:           &
       nspecies_chem =>nspecies,  & ! (IN)
       spc_alloc_chem=>spc_alloc, & ! (IN)
       spc_name_chem=>spc_name,   & ! (IN)
       ddp,                       & ! (IN)
       transport,                 & ! (IN)
       off,                       & ! (IN)
       O3,                        & ! (IN)
       SO2,                       & ! (IN)
       NH3,                       & ! (IN)
       SULF,                      & ! (IN)
       dvj,                       & ! (IN)
       hstar,                     & ! (IN)
       ak0,                       & ! (IN)
       dak,                       & ! (IN)
       dhr,                       & ! (IN)
       f0                           ! (IN)

  USE aer1_list, ONLY:            &
       nspecies_aer=> nspecies,   & ! (IN)
       spc_alloc_aer=>spc_alloc,  & ! (IN)
       part_radius,               & ! (IN)
       part_dens,                 & ! (IN)
       mode_alloc,                & ! (IN)
       on,                        & ! (IN)
       nmodes                       ! (IN)

  USE mem_chem1, ONLY:            &
       chem1_vars,                & ! Type
       chem1_g                      ! (INOUT)

  USE mem_aer1, ONLY:             &
       aer1_vars,                 & ! Type
       aer1_g                       ! (INOUT)

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
  USE ModDateUtils, ONLY:         &
       julday
!--(DMK-CCATT-BRAMS-5.0-OLD)------------------------------------------------------------------
!  USE mod_therm_lib, ONLY: &
!       rs                  ! Function
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------


  IMPLICIT NONE

  PRIVATE


  INTEGER, PARAMETER :: dep_seasons = 5
  INTEGER, PARAMETER :: nlu = 25
  INTEGER, PARAMETER :: nseason = 1
  INTEGER, PARAMETER :: nseasons = 2

  ! following currently hardwired to leaf-3 classes
  INTEGER, PARAMETER :: isice=2
  INTEGER, PARAMETER :: iswater=1
  !- aerosol section
  INTEGER, PARAMETER :: ntotal=nspecies_aer*nmodes
  INTEGER, PUBLIC :: naer_a,naer_z
  INTEGER, PUBLIC :: ind_aer(ntotal)
  INTEGER, PUBLIC :: ind_mode(ntotal)
  INTEGER, PUBLIC :: NAER_TRANSPORTED
  INTEGER, PUBLIC :: NAER_TRANSPORTED_numb

  TYPE, PUBLIC :: sedim_type

     REAL,POINTER :: v_sed_part(:,:,:)
     REAL,POINTER :: r_lsl_part(:,:,:)
     REAL,POINTER :: v_dep_part(:,:)

  END TYPE sedim_type

  TYPE(sedim_type), PUBLIC, ALLOCATABLE :: dd_sedim(:,:)
  TYPE(sedim_type), PUBLIC, ALLOCATABLE :: dd_sedim_numb(:)
  LOGICAL :: aer_alloc = .FALSE., aer_alloc_numb = .FALSE.

  INTEGER :: ixxxlu(nlu)
  REAL    :: kpart(nlu)
  REAL    :: rac(nlu,dep_seasons)
  REAL    :: rclo(nlu,dep_seasons)
  REAL    :: rcls(nlu,dep_seasons)
  REAL    :: rgso(nlu,dep_seasons)
  REAL    :: rgss(nlu,dep_seasons)
  REAL    :: ri(nlu,dep_seasons)
  REAL    :: rlu(nlu,dep_seasons)
  REAL    :: dratio(nspecies_chem)

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
  REAL,    EXTERNAL :: rs        ! Function
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

  PUBLIC :: dep_init             ! Subroutine
  PUBLIC :: alloc_aer_sedim
  PUBLIC :: alloc_aer_sedim_numb
  PUBLIC :: fa_preptc_with_sedim ! Subroutine
  PUBLIC :: dry_dep              ! Subroutine


CONTAINS

  !=========================================================================
  SUBROUTINE alloc_aer_sedim(npatch,ngrids &
                            ,nmodes,nspecies_aer    &
                            ,mode_alloc,on          &
                            ,mmzp,mmxp,mmyp)

    INTEGER , INTENT(IN) :: npatch
    INTEGER , INTENT(IN) :: ngrids

    ! aer1_list
    INTEGER , INTENT(IN) :: nmodes
    INTEGER , INTENT(IN) :: nspecies_aer
    INTEGER , INTENT(IN) :: mode_alloc(nmodes,nspecies_aer)
    INTEGER , INTENT(IN) :: on

    ! node_mod
    INTEGER , INTENT(IN) :: mmxp(ngrids)
    INTEGER , INTENT(IN) :: mmyp(ngrids)
    INTEGER , INTENT(IN) :: mmzp(ngrids)

    INTEGER :: ispc,imode,ng


    IF(aer_alloc) THEN
       PRINT *,'ERROR: aer_alloc already allocated'
       PRINT *,'Routine: aer_alloc File: chem_dry_dep.f90'
       STOP
    END IF

    NAER_TRANSPORTED = 0

    ! aerosol section : mapping
    naer_a = 1
    naer_z = 0

    DO ispc = 1,nspecies_aer
       DO imode=1,nmodes

          IF(mode_alloc   (imode,ispc) == on ) THEN

             naer_z = naer_z + 1
             ind_aer (naer_z) = ispc
             ind_mode(naer_z) = imode

          ENDIF
       ENDDO
    ENDDO

    ! total number of species (aer) to be transported
    NAER_TRANSPORTED  =  naer_z
    IF(NAER_TRANSPORTED == 0) RETURN

    ALLOCATE (dd_sedim(NAER_TRANSPORTED,ngrids))
    DO ng=1,ngrids
       DO ispc=1,NAER_TRANSPORTED
          ALLOCATE(dd_sedim(ispc,ng)%v_sed_part (mmzp(ng),mmxp(ng),mmyp(ng)))
          dd_sedim(ispc,ng)%v_sed_part = 0.
          ALLOCATE(dd_sedim(ispc,ng)%r_lsl_part (mmxp(ng),mmyp(ng),npatch))
          dd_sedim(ispc,ng)%r_lsl_part = 0.
          ALLOCATE(dd_sedim(ispc,ng)%v_dep_part (mmxp(ng),mmyp(ng)))
          dd_sedim(ispc,ng)%v_dep_part = 0.
       ENDDO
    ENDDO


    aer_alloc=.TRUE.

  END SUBROUTINE alloc_aer_sedim


    !=========================================================================
    SUBROUTINE alloc_aer_sedim_numb(npatch,ngrids &
                              ,nmodes  &
                              ,mode_alloc,on          &
                              ,mmzp,mmxp,mmyp)

      INTEGER , INTENT(IN) :: npatch
      INTEGER , INTENT(IN) :: ngrids

      ! aer1_list
      INTEGER , INTENT(IN) :: nmodes
      INTEGER , INTENT(IN) :: mode_alloc(nmodes)
      INTEGER , INTENT(IN) :: on

      ! node_mod
      INTEGER , INTENT(IN) :: mmxp(ngrids)
      INTEGER , INTENT(IN) :: mmyp(ngrids)
      INTEGER , INTENT(IN) :: mmzp(ngrids)

      INTEGER :: ispc,imode,ng


      IF(aer_alloc_numb) THEN
         PRINT *,'ERROR: aer_alloc for number already allocated'
         PRINT *,'Routine: aer_alloc File: chem_dry_dep.f90'
         STOP
      END IF

      NAER_TRANSPORTED = 0

      ! aerosol section : mapping
      naer_a = 1
      naer_z = 0

      !DO ispc = 1,nspecies_aer
         DO imode=1,nmodes

            IF(mode_alloc   (imode) == on ) THEN

               naer_z = naer_z + 1
               ind_aer (naer_z) = ispc
               ind_mode(naer_z) = imode

            ENDIF
         ENDDO
      !ENDDO

      ! total number of species (aer) to be transported
      NAER_TRANSPORTED_numb  =  naer_z
      IF(NAER_TRANSPORTED_numb == 0) RETURN

      ALLOCATE (dd_sedim_numb(ngrids))
      DO ng=1,ngrids
         DO ispc=1,NAER_TRANSPORTED
            ALLOCATE(dd_sedim_numb(ng)%v_sed_part (mmzp(ng),mmxp(ng),mmyp(ng)))
            dd_sedim_numb(ng)%v_sed_part = 0.
            ALLOCATE(dd_sedim_numb(ng)%r_lsl_part (mmxp(ng),mmyp(ng),npatch))
            dd_sedim_numb(ng)%r_lsl_part = 0.
            ALLOCATE(dd_sedim_numb(ng)%v_dep_part (mmxp(ng),mmyp(ng)))
            dd_sedim_numb(ng)%v_dep_part = 0.
         ENDDO
      ENDDO


      aer_alloc_numb=.TRUE.

    END SUBROUTINE alloc_aer_sedim_numb


  !========================================================================
  SUBROUTINE dry_dep(m1,m2,m3,ia,iz,ja,jz                 &
                    ,cpi,cpor,p00,g,vonk                  &
                    ,jdim,dzt,zt,nzpmax                   &
                    ,npatch,dt,level,nnqparm              &
                    ,imonth1,idate1,iyear1                &
                    ,chemistry,aerosol                    &
                    ,theta,rv,pp,dn0,pi0,up,vp,tke        &
                    ,sflux_t,sflux_r,sflux_u,sflux_v      &
                    ,r_aer,ustar,tstar,patch_area,veg,Z0m &
                    ,rcp,pcpg,rtgt,rshort,conprr          &
                    ,aer1_g                               &
                    ,chem1_g                              &
                    ,dd_sedim)

    USE Extras, ONLY : extra2d,NA_EXTRA2D
    INTEGER , INTENT(IN)    :: m1
    INTEGER , INTENT(IN)    :: m2
    INTEGER , INTENT(IN)    :: m3
    INTEGER , INTENT(IN)    :: ia
    INTEGER , INTENT(IN)    :: iz
    INTEGER , INTENT(IN)    :: ja
    INTEGER , INTENT(IN)    :: jz

    ! rconstants
    REAL    , INTENT(IN)    :: cpi
    REAL    , INTENT(IN)    :: cpor
    REAL    , INTENT(IN)    :: p00
    REAL    , INTENT(IN)    :: g
    REAL    , INTENT(IN)    :: vonk

    ! mem_grid
    INTEGER , INTENT(IN)    :: npatch
    INTEGER , INTENT(IN)    :: jdim
    REAL    , INTENT(IN)    :: dt
    REAL    , INTENT(IN)    :: dzt(nzpmax)
    REAL    , INTENT(IN)    :: zt(nzpmax)
    INTEGER , INTENT(IN)    :: nzpmax
    INTEGER , INTENT(IN)    :: imonth1
    INTEGER , INTENT(IN)    :: idate1
    INTEGER , INTENT(IN)    :: iyear1

    ! micphys
    INTEGER , INTENT(IN)    :: level

    ! mem_basic
    REAL    , INTENT(IN)    :: theta(m1,m2,m3)
    REAL    , INTENT(IN)    :: rv(m1,m2,m3)
    REAL    , INTENT(IN)    :: pp(m1,m2,m3)
    REAL    , INTENT(IN)    :: dn0(m1,m2,m3)
    REAL    , INTENT(IN)    :: pi0(m1,m2,m3)
    REAL    , INTENT(IN)    :: up(m1,m2,m3)
    REAL    , INTENT(IN)    :: vp(m1,m2,m3)

    ! mem_turb
    REAL    , INTENT(IN)    :: tke(m1,m2,m3)
    REAL    , INTENT(IN)    :: sflux_t(m2,m3)
    REAL    , INTENT(IN)    :: sflux_r(m2,m3)
    REAL    , INTENT(IN)    :: sflux_u(m2,m3)
    REAL    , INTENT(IN)    :: sflux_v(m2,m3)

    ! mem_radiate
    REAL    , INTENT(IN)    :: rshort(m2,m3)

    ! mem_micro
    REAL    , INTENT(IN)    :: rcp(m1,m2,m3)
    REAL    , INTENT(IN)    :: pcpg(m2,m3)

    ! mem_grid
    REAL    , INTENT(IN)    :: rtgt(m2,m3)

    ! mem_leaf
    REAL    , INTENT(IN)    :: r_aer(m2,m3,npatch)
    REAL    , INTENT(IN)    :: ustar(m2,m3,npatch)
    REAL    , INTENT(IN)    :: tstar(m2,m3,npatch)
    REAL    , INTENT(IN)    :: patch_area(m2,m3,npatch)
    REAL    , INTENT(IN)    :: veg(m2,m3,npatch)
    REAL    , INTENT(IN)    :: Z0m(m2,m3,npatch)

    ! mem_cuparm
    INTEGER , INTENT(IN)    :: nnqparm
    REAL    , INTENT(IN)    :: conprr(m2,m3)

    ! mem_aer1
    INTEGER         , INTENT(IN)    :: aerosol
    TYPE(aer1_vars) , INTENT(INOUT) :: aer1_g(nmodes,nspecies_aer)

    ! mem_chem1
    INTEGER         , INTENT(IN)    :: chemistry
    TYPE(chem1_vars), INTENT(INOUT) :: chem1_g(nspecies_chem)

    TYPE(sedim_type), INTENT(INOUT) :: dd_sedim(naer_transported)


    REAL, PARAMETER ::  ubmin = 0.25

    INTEGER :: i,j,ipatch,idry_part,ispc,k
    REAL    :: ups,vps,pis,sbf

    REAL    :: v_dep_gas(nspecies_chem,m2,m3)
    REAL    :: check_rain(m2,m3              )
    REAL    :: rmol(m2,m3              )!1./Monin-Obukhob
    REAL    :: rhchem(m2,m3)

    ! (DMK) old scratches vars
    REAL    :: temps(m2,m3), prss(m2,m3), dens(m2,m3), vels(m2,m3), rvs(m2,m3), Zi(m2,m3)
    REAL    :: temp3d(m1,m2,m3), air_dens3d(m1,m2,m3)

    !- no tracers, aerosols or chemistry, go back
    IF(CHEMISTRY < 0) RETURN

    v_dep_gas = 0.0
    check_rain = 0.
    rmol = 0.
    !-aux variables

    DO j = ja,jz
       DO i = ia,iz
          rvs  (i,j) = rv(2,i,j)
          pis        = ( pp(1,i,j) + pp(2,i,j) + pi0(1,i,j) + pi0(2,i,j) ) * .5 * cpi
          prss (i,j) = pis ** cpor * p00
          dens (i,j) = ( dn0(1,i,j) + dn0(2,i,j) ) * .5
          temps(i,j) = theta(2,i,j) * pis        ! temps=theta*Exner/CP
          rhchem(i,j)=100.*MIN(1.,MAX(0.05,rvs(i,j)/rs(prss(i,j),temps(i,j))))
          ups        = (up(2,i,j) + up(2,i-1,j   )) * .5
          vps        = (vp(2,i,j) + vp(2,i,j-jdim)) * .5
          vels (i,j) = MAX(ubmin,SQRT(ups** 2 + vps** 2))
          !- compute surface buoyancy flux [sbf] and 1 / Monin Obukhov height [1/zl]
          sbf = (g/theta(2,i,j) * (1. + .61 * rv(2,i,j))) &
               * ( sflux_t(i,j) * (1. + .61 * rv(2,i,j)) &
               +   sflux_r(i,j) * .61 * theta(2,i,j)   )
          !-  reciprocal of the Monin-Obukhov length (1/m)
          rmol (i,j) = - vonk*sbf / MAX(0.1,SQRT(SQRT(sflux_u(i,j)**2 + sflux_v(i,j)**2)))
          IF(ABS(rmol (i,j)) < 1.e-7) rmol (i,j) = 0.
          !    print*,'rmol=',rmol(i,j)
       ENDDO
    ENDDO
    IF (nnqparm > 0) THEN
       check_rain(ia:iz,ja:jz)  = conprr(ia:iz,ja:jz)
    ENDIF
    IF (level >= 3) THEN
       check_rain(ia:iz,ja:jz)  = check_rain(ia:iz,ja:jz)  +  pcpg(ia:iz,ja:jz)
    ENDIF

    CALL define_PBL_height(m1,m2,m3,npatch,ia,iz,ja,jz,zt     &
         ,tke,sflux_t,rcp,rtgt &
         ,zi)


    !---  dry deposition for gases
    CALL dry_dep_gases(m1,m2,m3,nspecies_chem,npatch,ia,iz,ja,jz  &
         ,v_dep_gas,r_aer  		           &
         ,prss,temps,dens,vels,rvs,rcp,Zi		   &
         ,ustar,tstar,patch_area,veg,Z0m		   &
         ,rshort,rtgt,dzt,check_rain,rmol,rhchem &
         ,O3,SULF,spc_alloc_chem,transport &
         ,on,dvj,hstar,ak0,dak,dhr,f0,imonth1,idate1,iyear1)

    !-apply dry deposition on the tracers concentration and get the deposited mass on surface
    DO ispc=1,nspecies_chem
       IF (spc_alloc_chem(ddp,ispc) == off .OR. &
           spc_alloc_chem(transport,ispc) == off) CYCLE

       !-temporary for save dep velocity of ozone
       !- to use it NA_EXTRA2D must be greater than zero
       !if(ispc==O3 .and. NA_EXTRA2D > 0)  then
       !  extra2d(1,1)%d2(1:m2,1:m3)=v_dep_gas(ispc,1:m2,1:m3)
       !endif

       CALL apply_drydep(m1,m2,m3,ia,iz,ja,jz		&
            ,v_dep_gas(ispc,1:m2,1:m3)	& ! don't change this
            ,chem1_g  (ispc)%sc_t     &! tendency array
            ,chem1_g  (ispc)%sc_dd    &! deposited mass
            ,chem1_g  (ispc)%sc_p     &! mixing ratio
            ,dens,rtgt,dzt,dt)
    ENDDO

    !--- aerosol dry deposition and sedimentation

    IF(AEROSOL == 1) THEN

       DO j = ja,jz
          DO i = ia,iz
             DO k = 1, m1-1
                temp3d(k,i,j) = 0.5*( theta(k+1,i,j) * ( pp(k+1,i,j)+pi0(k+1,i,j) ) * cpi + &
                     theta(k  ,i,j) * ( pp(k  ,i,j)+pi0(k  ,i,j) ) * cpi   )
                air_dens3d(k,i,j) = 0.5 *(  dn0(k+1,i,j) +  dn0(k  ,i,j) )
             ENDDO
             temp3d(m1,i,j) =	 temp3d(m1-1,i,j)

             air_dens3d(m1,i,j) =    air_dens3d(m1-1,i,j)
          ENDDO
       ENDDO
       !print*,'temp=',temp3d(1:2,int(iz/2),int(jz/2))

       !--  dry deposition/sedimentation for particles
       CALL dry_dep_sedim_particles(m1,m2,m3,nspecies_aer              &
                                   ,npatch,ia,iz,ja,jz                 &
                                   ,r_aer,temp3d,air_dens3d,temps,dens &
                                   ,vels,rvs,Zi,ustar,tstar,patch_area &
                                   ,veg,Z0m,nmodes,nspecies_aer        &
                                   ,part_radius,part_dens,mode_alloc   &
                                   ,on,g,dd_sedim)

       !- loop over the transported aerosols particles
       DO ispc=naer_a,naer_z


          IF(spc_alloc_aer(ddp      ,ind_mode(ispc),ind_aer(ispc)) == off .OR. &
               spc_alloc_aer(transport,ind_mode(ispc),ind_aer(ispc)) == off) CYCLE


          !print*,'dd=',dd_sedim(ispc,ngrid)%v_dep_part(int(iz/2),int(jz/2))*100.

          CALL apply_drydep(m1,m2,m3,ia,iz,ja,jz	       &
               ,dd_sedim(ispc)%v_dep_part   	               &
               ,aer1_g  (ind_mode(ispc),ind_aer(ispc))%sc_t    &! tendency array
               ,aer1_g  (ind_mode(ispc),ind_aer(ispc))%sc_dd   &! deposited mass
               ,aer1_g  (ind_mode(ispc),ind_aer(ispc))%sc_p    &! mixing ratio
               ,dens,rtgt,dzt,dt)

       ENDDO
    ENDIF  ! endif of aerosol

  END SUBROUTINE dry_dep
  !========================================================================

  SUBROUTINE dry_dep_gases(m1,m2,m3,nspecies_chem,npatch,ia,iz,ja,jz         &
                          ,V_dep,r_aer,prss,temps,dens,vels,rvs,rcp,Zi       &
                          ,ustar,tstar,patch_area,veg,Z0m,rshort,rtgt,dzt    &
                          ,check_rain,rmol,rhchem,O3,SULF,spc_alloc_chem,transport&
                          ,on,dvj,hstar,ak0,dak,dhr,f0,imonth1,idate1,iyear1)

    INTEGER , INTENT(IN)    :: m1
    INTEGER , INTENT(IN)    :: m2
    INTEGER , INTENT(IN)    :: m3
    INTEGER , INTENT(IN)    :: nspecies_chem
    INTEGER , INTENT(IN)    :: npatch
    INTEGER , INTENT(IN)    :: ia
    INTEGER , INTENT(IN)    :: iz
    INTEGER , INTENT(IN)    :: ja
    INTEGER , INTENT(IN)    :: jz
    REAL    , INTENT(INOUT) :: V_dep(nspecies_chem,m2,m3)
    REAL    , INTENT(IN)    :: r_aer(m2,m3,npatch)
    REAL    , INTENT(IN)    :: prss(m2,m3)
    REAL    , INTENT(IN)    :: temps(m2,m3)

    REAL    , INTENT(IN)    :: dens(m2,m3)         ! (DMK) not used
    REAL    , INTENT(IN)    :: vels(m2,m3)         ! (DMK) not used
    REAL    , INTENT(IN)    :: rvs(m2,m3)          ! (DMK) not used
    REAL    , INTENT(IN)    :: rcp(m1,m2,m3)
    REAL    , INTENT(IN)    :: Zi(m2,m3)           ! (DMK) not used
    REAL    , INTENT(IN)    :: ustar(m2,m3,npatch)
    REAL    , INTENT(IN)    :: tstar(m2,m3,npatch) ! (DMK) not used
    REAL    , INTENT(IN)    :: patch_area(m2,m3,npatch)
    REAL    , INTENT(IN)    :: veg(m2,m3,npatch)
    REAL    , INTENT(IN)    :: Z0m(m2,m3,npatch)   ! (DMK) not used
    REAL    , INTENT(IN)    :: rshort(m2,m3)
    REAL    , INTENT(IN)    :: rtgt(m2,m3)
    REAL    , INTENT(IN)    :: dzt(m1)
    REAL    , INTENT(IN)    :: check_rain(m2,m3)
    REAL    , INTENT(IN)    :: rmol(m2,m3)
    REAL    , INTENT(IN)    :: rhchem(m2,m3)

    ! chem1_list
    INTEGER , INTENT(IN)    :: O3
    INTEGER , INTENT(IN)    :: SULF
    INTEGER , INTENT(IN)    :: spc_alloc_chem(6,nspecies_chem)
    INTEGER , INTENT(IN)    :: transport
    INTEGER , INTENT(IN)    :: on
    REAL    , INTENT(IN)    :: dvj(nspecies_chem)
    REAL    , INTENT(IN)    :: hstar(nspecies_chem)
    REAL    , INTENT(IN)    :: ak0(nspecies_chem)
    REAL    , INTENT(IN)    :: dak(nspecies_chem)
    REAL    , INTENT(IN)    :: dhr(nspecies_chem)
    REAL    , INTENT(IN)    :: f0(nspecies_chem)

    ! mem_grid
    INTEGER , INTENT(IN)    :: imonth1
    INTEGER , INTENT(IN)    :: idate1
    INTEGER , INTENT(IN)    :: iyear1

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
!    INTEGER,EXTERNAL :: julday
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

    INTEGER :: i,j,ipatch

    !- local var
    REAL :: vgs0d(nspecies_chem)  !,srfres
    INTEGER iveg,iseason,jday
    LOGICAL :: highnh3, rainflag, vegflag, wetflag
    REAL dvpart,dvfog,clwchem,z1,ustarw,vgpart
    !real , dimension(m2,m3) :: aer_res
    !real   ::  r_lsl(m2,m3,npatch,nspecies_chem),rcx(npatch,nspecies_chem)
    REAL   ::  r_lsl(nspecies_chem),rcx(nspecies_chem)

    !   Set the reference height (10.0 m)
    REAL, PARAMETER ::  zr = 10.0


    jday=julday(imonth1,idate1,iyear1)
    iseason = 1
    IF( jday.LT.90 .OR. jday.GT.270)  iseason=2

    DO j = ja,jz
       DO i = ia,iz
          ustarw = 0.
          !
          !     Set logical default values
          rainflag = .FALSE.
          IF (check_rain(i,j) > 0.  ) rainflag = .TRUE.

          wetflag  = .FALSE.
          IF (rhchem(i,j) >= 95.) wetflag  = .TRUE.

          clwchem = rcp(2,i,j)

          !-not using
          highnh3 = .FALSE.
          !      if(chem(i,kts,j,p_nh3).gt.2.*chem(i,kts,j,p_so2))highnh3 = .true.
          !----


          !zntt = znt(i,j)  ! Surface roughness height (m)

          z1 =   rtgt(i,j)/dzt(2)! height of the fisrt cell

          DO ipatch = 1,npatch
             IF (patch_area(i,j,ipatch) .GE. .009) THEN

                IF(ipatch == 1) THEN
                   iveg = 1           ! rc routine needs water as 1
                ELSE
                   iveg = NINT(veg(i,j,ipatch))
                ENDIF

                ustarw=ustarw+ustar(i,j,ipatch)*patch_area(i,j,ipatch)
                !
                !- calculates surface resistance (rcx)
                !

                CALL rc(rcx,temps(i,j),rshort(i,j),rhchem(i,j),iveg,iseason   &
                     ,nspecies_chem,wetflag,rainflag,highnh3,spc_alloc_chem   &
                     ,hstar,ak0,dak,dhr,f0,transport,O3,on)

                !- get R_b
                CALL get_rb(nspecies_chem,r_lsl,temps(i,j) &
                     ,prss(i,j),ustar(i,j,ipatch),dvj)
                !
                !- deposition velocity for each patch
                vgs0d(1:nspecies_chem) = 1./( r_aer(i,j,ipatch     )  + &
                     r_lsl(1:nspecies_chem)  + &
                     rcx  (1:nspecies_chem)    )
                !
                !- special treatment for SULF
	    if( spc_alloc_chem(transport,SULF) == ON) then
             CALL deppart(rmol(i,j),ustar(i,j,ipatch),rhchem(i,j),clwchem,iveg,dvpart,dvfog)
             CALL depvel(nspecies_chem,rmol(i,j),zr,Z0m(i,j,ipatch)   &
	               ,ustar(i,j,ipatch),vgpart)
             vgs0d(sulf)=1.0/((1.0/vgpart)+(1.0/dvpart))
	    endif
                !- end SULF treatment
                !
                !
                !- effective deposition velocity for all patches
                V_dep(1:nspecies_chem,i,j) = V_dep(1:nspecies_chem,i,j) + &
                     patch_area(i,j,ipatch)*vgs0d(1:nspecies_chem)

             ENDIF
          ENDDO

          !	print*,'O3-1 cm/s=',V_dep(O3,i,j)*100; call flush(6)
          !- get the effective deposition velocity
          CALL cellvg(V_dep(1:nspecies_chem,i,j),ustarw,z1,zr,rmol(i,j),nspecies_chem)
          !	print*,'O3-2 cm/s=',V_dep(O3,i,j)*100; call flush(6)



       ENDDO
    ENDDO


    RETURN
  END SUBROUTINE dry_dep_gases
  !========================================================================

  SUBROUTINE get_rb(nspecies_chem,r_lsl,temps,prss,ustar,dvj)

    INTEGER , INTENT(IN)  :: nspecies_chem
    REAL    , INTENT(OUT) :: r_lsl(nspecies_chem)
    REAL    , INTENT(IN)  :: temps
    REAL    , INTENT(IN)  :: prss
    REAL    , INTENT(IN)  :: ustar

    ! chem1_list
    REAL    , INTENT(IN)  :: dvj(nspecies_chem)

    !local var
    REAL, PARAMETER :: i23=2./3., STDTEMP = 273.15, STDATMPA = 101325.0
    INTEGER ispc
    REAL dvj_tmp,sc,scpr23,eta


    !-  kinematic viscosity of air at surface (cm^2/s)
    eta = -1.1555E-10*temps**3 + 9.5728E-07*temps**2 + 3.7604E-04*temps - 3.4484E-02

    DO ispc=1,nspecies_chem

       !         dvj(ispc) = dvj(ispc)*(293.15/298.15)**1.75
       !         dvj_tmp = dvj(ispc)*(temps/(temps+5.))**1.75
       dvj_tmp = dvj(ispc)* (( temps / STDTEMP ) ** 1.75)&
            * (STDATMPA / prss ) ! cm^2/s
       !          dratio(ispc) = 0.242/dvj(ispc)

       !         sc = 0.15/dvj_tmp ! Schmidt Number at 20�C
       sc =  eta/dvj_tmp ! Schmidt Number at any temps

       scpr23 = (sc/0.72)**(i23) ! (Schmidt # / Prandtl #)**
       r_lsl(ispc) = 5.*scpr23/ustar
    ENDDO
  END SUBROUTINE get_rb

  ! **********************************************************************
  ! **************************  SUBROUTINE RC  ***************************
  ! **********************************************************************
  SUBROUTINE rc(rcx,t,rad,rh,iland,iseason,nspecies_chem,        &
                wetflag,rainflag,highnh3,spc_alloc_chem,hstar,   &
                ak0,dak,dhr,f0,transport,O3,on )
    !     THIS SUBROUTINE CALCULATES SURFACE RESISTENCES ACCORDING
    !     TO THE MODEL OF
    !     M. L. WESELY,
    !     ATMOSPHERIC ENVIRONMENT 23 (1989), 1293-1304
    !     WITH SOME ADDITIONS ACCORDING TO
    !     J. W. ERISMAN, A. VAN PUL, AND P. WYERS,
    !     ATMOSPHERIC ENVIRONMENT 28 (1994), 2595-2607
    !     WRITTEN BY  WINFRIED SEIDL, APRIL 1997
    !     MODYFIED BY WINFRIED SEIDL, MARCH 2000
    !                    FOR MM5 VERSION 3
    !----------------------------------------------------------------------

    REAL    , INTENT(IN)    :: t
    REAL    , INTENT(IN)    :: rad
    REAL    , INTENT(IN)    :: rh ! (DMK) not used
    INTEGER , INTENT(IN)    :: iland   ! (DMK) IN
    INTEGER , INTENT(IN)    :: iseason ! (DMK) IN
    INTEGER , INTENT(IN)    :: nspecies_chem
    LOGICAL , INTENT(IN)    :: wetflag
    LOGICAL , INTENT(IN)    :: rainflag
    LOGICAL , INTENT(IN)    :: highnh3 ! (DMK) not used
    REAL    , INTENT(INOUT) :: rcx(nspecies_chem)

    ! chem1_list
    INTEGER , INTENT(IN)    :: spc_alloc_chem(6,nspecies_chem)
    REAL    , INTENT(IN)    :: hstar(nspecies_chem)
    REAL    , INTENT(IN)    :: ak0(nspecies_chem)
    REAL    , INTENT(IN)    :: dak(nspecies_chem)
    REAL    , INTENT(IN)    :: dhr(nspecies_chem)
    REAL    , INTENT(IN)    :: f0(nspecies_chem)
    INTEGER , INTENT(IN)    :: transport
    INTEGER , INTENT(IN)    :: O3
    INTEGER , INTENT(IN)    :: on

!------------old--------------------------------------------
! .. Scalar Arguments ..
!        REAL :: rad, rh, t, temp_i
!        INTEGER :: iland, iseason, nspecies_chem
!        LOGICAL :: highnh3, rainflag, wetflag
!        real :: rcx(nspecies_chem)
!------------old--------------------------------------------

    REAL :: temp_i
    ! ..
    ! .. Array Arguments ..
    INTEGER :: iprt
    ! ..
    ! .. Local Scalars ..
    REAL :: rclx, rdc, resice, rgsx, rluo1, rluo2, rlux, rmx, rs, rsmx, &
         tc, rdtheta, z, hplus, corrh
    INTEGER :: n
    ! ..
    ! .. Local Arrays ..
    REAL :: hstary(nspecies_chem)
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC exp
    ! ..
    !INTEGER, PARAMETER :: SO2=10000

    DO n = 1, nspecies_chem
       rcx(n) = 1.
    END DO

    tc = t - 273.15
    rdtheta = 0.
    temp_i = 1./t-1./298.

    z = 200./(rad+0.1)

!!!  HARDWIRE VALUES FOR TESTING
    !       z=0.4727409
    !       tc=22.76083
    !       t=tc+273.15
    !       rad = 412.8426
    !       rainflag=.false.
    !       wetflag=.false.

    IF ((tc<=0.) .OR. (tc>=40.)) THEN
       rs = 9999.
    ELSE

!--(DMK-CCATT-INI)-------------------------------------------------------------
       rs = ri(ixxxlu(iland),iseason)*(1+z*z)*(400./(tc*(40.-tc)))
!--(DMK-CCATT-OLD)-------------------------------------------------------------
!       rs = ri(iland,iseason)*(1+z*z)*(400./(tc*(40.-tc)))
!--(DMK-CCATT-FIM)-------------------------------------------------------------

    END IF
    rdc = 100*(1.+1000./(rad+10))/(1+1000.*rdtheta)

!--(DMK-CCATT-INI)-------------------------------------------------------------
    rluo1 = 1./(1./3000.+1./3./rlu(ixxxlu(iland),iseason))
    rluo2 = 1./(1./1000.+1./3./rlu(ixxxlu(iland),iseason))
!--(DMK-CCATT-OLD)-------------------------------------------------------------
!    rluo1 = 1./(1./3000.+1./3./rlu(iland,iseason))
!    rluo2 = 1./(1./1000.+1./3./rlu(iland,iseason))
!--(DMK-CCATT-FIM)-------------------------------------------------------------

    resice = 1000.*EXP(-tc-4.)

    ! change MP 11/12/07 taking into account the acid dissociation constant
    ! and assuming pH=7 in the vegetation
    !srf-opt  hstary(n) = hstar(n)*exp(dhr(n)*(1./t-1./298.))
    !                    *(1+ak0(n)*exp(dak(n)*(1./t-1/298.))/hplus)
    !          pH=7  hplus=10**(-pH)
    hplus=1.E-7
    DO n = 1, nspecies_chem
       !srf      IF (hstar(n)==0.) GO TO 10
       IF (hstar(n)==0.) CYCLE
       !srf
       corrh=1+ak0(n)*EXP(dak(n)*(temp_i))/hplus
       hstary(n) = hstar(n)*EXP(dhr(n)*(temp_i))*corrh
       !
       !         if ((n == 8).or.(n == 52))then
       !           print*, n, hstary(n),t,temp_i,ak0(n), dak(n),hplus,corrh &
       !                , hstar(n), dhr(n)
       !         endif
       ! end change MP
       rmx = 1./(hstary(n)/3000.+100.*f0(n))
       rsmx = rs*dratio(n) + rmx

!--(DMK-CCATT-INI)-------------------------------------------------------------
          rclx = 1./(hstary(n)/1.E+5/rcls(ixxxlu(iland),iseason)+f0(n)/rclo(ixxxlu(iland), &
            iseason)) + resice
          rgsx = 1./(hstary(n)/1.E+5/rgss(ixxxlu(iland),iseason)+f0(n)/rgso(ixxxlu(iland), &
            iseason)) + resice
          rlux = rlu(ixxxlu(iland),iseason)/(1.E-5*hstary(n)+f0(n)) + resice
          IF (wetflag) THEN
            rlux = 1./(1./3./rlu(ixxxlu(iland),iseason)+1.E-7*hstary(n)+f0(n)/rluo1)
          END IF
        IF (rainflag) THEN
            rlux = 1./(1./3./rlu(ixxxlu(iland),iseason)+1.E-7*hstary(n)+f0(n)/rluo2)
          END IF
          rcx(n) = 1./(1./rsmx+1./rlux+1./(rdc+rclx)+1./(rac(ixxxlu(iland), &
            iseason)+rgsx))
!--(DMK-CCATT-OLD)-------------------------------------------------------------
!       rclx = 1./(hstary(n)/1.E+5/rcls(iland,iseason)+f0(n)/rclo(iland, &
!            iseason)) + resice
!       rgsx = 1./(hstary(n)/1.E+5/rgss(iland,iseason)+f0(n)/rgso(iland, &
!            iseason)) + resice
!       rlux = rlu(iland,iseason)/(1.E-5*hstary(n)+f0(n)) + resice
!       IF (wetflag) THEN
!          rlux = 1./(1./3./rlu(iland,iseason)+1.E-7*hstary(n)+f0(n)/rluo1)
!       END IF
!       IF (rainflag) THEN
!          rlux = 1./(1./3./rlu(iland,iseason)+1.E-7*hstary(n)+f0(n)/rluo2)
!       END IF
!       rcx(n) = 1./(1./rsmx+1./rlux+1./(rdc+rclx)+1./(rac(iland, &
!            iseason)+rgsx))
!--(DMK-CCATT-FIM)-------------------------------------------------------------

       IF (rcx(n)<1.) rcx(n) = 1.
10  END DO
    !

    IF( spc_alloc_chem(transport,O3) == ON) THEN

       !--------------------------------- O3
       !        SPECIAL TREATMENT FOR OZONE
       !srf       hstary(O3) = hstar(O3)*exp(dhr(O3)*(1./t-1./298.))
       hstary(O3) = hstar(O3)*EXP(dhr(O3)*(temp_i))
       rmx = 1./(hstary(O3)/3000.+100.*f0(O3))
       rsmx = rs*dratio(O3) + rmx

!--(DMK-CCATT-INI)-------------------------------------------------------------
       rlux = rlu(ixxxlu(iland),iseason)/(1.E-5*hstary(O3)+f0(O3)) + resice
       rclx = rclo(ixxxlu(iland),iseason) + resice
       rgsx = rgso(ixxxlu(iland),iseason) + resice
!--(DMK-CCATT-OLD)-------------------------------------------------------------
!       rlux = rlu(iland,iseason)/(1.E-5*hstary(O3)+f0(O3)) + resice
!       rclx = rclo(iland,iseason) + resice
!       rgsx = rgso(iland,iseason) + resice
!--(DMK-CCATT-FIM)-------------------------------------------------------------

       IF (wetflag) rlux = rluo1
       IF (rainflag) rlux = rluo2

!--(DMK-CCATT-INI)-------------------------------------------------------------
        rcx(O3) = 1./(1./rsmx+1./rlux+1./(rdc+rclx)+1./(rac(ixxxlu(iland), &
          iseason)+rgsx))
!--(DMK-CCATT-OLD)-------------------------------------------------------------
!       rcx(O3) = 1./(1./rsmx+1./rlux+1./(rdc+rclx)+1./(rac(iland, &
!            iseason)+rgsx))
!--(DMK-CCATT-FIM)-------------------------------------------------------------

       rcx(O3)=MAX(rcx(O3),1.) !- old: IF (rcx(O3)<1.) rcx(O3) = 1.
    ENDIF


      if( spc_alloc_chem(transport,SO2) == ON) then
!!$
!!$!    	      SPECIAL TREATMENT FOR SO2 (Wesely)
!!$!    		HSTARY(SO2)=HSTAR(SO2)*EXP(DHR(SO2)*(1./T-1./298.))
!!$!    		RMX=1./(HSTARY(SO2)/3000.+100.*F0(SO2))
!!$!    		RSMX=RS*DRATIO(SO2)+RMX
!!$!    		RLUX=RLU(ILAND,ISEASON)/(1.E-5*HSTARY(SO2)+F0(SO2))
!!$!    	     &       +RESICE
!!$!    		RCLX=RCLS(ILAND,ISEASON)+RESICE
!!$!    		RGSX=RGSS(ILAND,ISEASON)+RESICE
!!$!    		IF ((wetflag).OR.(RAINFLAG)) THEN
!!$!    		  IF (ILAND.EQ.1) THEN
!!$!    		    RLUX=50.
!!$!    		  ELSE
!!$!    		    RLUX=100.
!!$!    		  END IF
!!$!    		END IF
!!$!    		RCX(SO2)=1./(1./RSMX+1./RLUX+1./(RDC+RCLX)
!!$!    	     &  	      +1./(RAC(ILAND,ISEASON)+RGSX))
!!$!    		IF (RCX(SO2).LT.1.) RCX(SO2)=1.
!!$!
!!$!
!!$!----------------------------------------- SO2
!!$!    	      SO2 according to Erisman et al. 1994
!!$!    		R_STOM
    		rsmx = rs*dratio(SO2)
!!$!    		R_EXT
     		IF (tc>(-1.)) THEN
     		  IF (rh<81.3) THEN
     		    rlux = 25000.*exp(-0.0693*rh)
     		  ELSE
     		    rlux = 0.58E12*exp(-0.278*rh)
     		  END IF
     		END IF
     		IF (((wetflag) .OR. (rainflag)) .AND. (tc>(-1.))) THEN
     		  rlux = 1.
     		END IF
     		IF ((tc>=(-5.)) .AND. (tc<=(-1.))) THEN
     		  rlux = 200.
     		END IF
     		IF (tc<(-5.)) THEN
     		  rlux = 500.
     		END IF
!!$!    		INSTEAD OF R_INC R_CL and R_DC of Wesely are used

!--(DMK-CCATT-INI)-------------------------------------------------------------
                rclx = rcls(ixxxlu(iland),iseason)
!--(DMK-CCATT-OLD)-------------------------------------------------------------
!     		rclx = rcls(iland,iseason)
!--(DMK-CCATT-FIM)-------------------------------------------------------------

!!$!    		DRY SURFACE
     		rgsx = 1000.
!!$!    		WET SURFACE
     		IF ((wetflag) .OR. (rainflag)) THEN
     		  IF (highnh3) THEN
     		    rgsx = 0.
     		  ELSE
     		    rgsx = 500.
     		  END IF
     		END IF
!!$!    		WATER
     		IF (iland==iswater) THEN
     		  rgsx = 0.
     		END IF
!!$!    		SNOW
     		IF ((iseason==4) .OR. (iland==isice)) THEN
     		  IF (tc>2.) THEN
     		    rgsx = 0.
     		  END IF
     		  IF ((tc>=(-1.)) .AND. (tc<=2.)) THEN
     		    rgsx = 70.*(2.-tc)
     		  END IF
     		  IF (tc<(-1.)) THEN
     		    rgsx = 500.
     		  END IF
     		END IF
!!$!    		TOTAL SURFACE RESISTENCE
     		IF ((iseason/=4) .AND. (ixxxlu(iland)/=1) .AND. (iland/=iswater) .AND. &
     		    (iland/=isice)) THEN
     		  rcx(SO2) = 1./(1./rsmx+1./rlux+1./(rclx+rdc+rgsx))
     		ELSE
     		  rcx(SO2) = rgsx
     		END IF
     		IF (rcx(SO2)<1.) rcx(SO2) = 1.

      endif

      !KML DEBUG
 !      print*, 'KML, iland, ixxxlu, rcx(O3)',iland, ixxxlu(iland), rcx(O3)
    !------------ NOT IN USE BELOW -----------------------------------------------
    !srf- dry dep for NH3 NOT in use :
    !RETURN
    !srf ---

    !     NH3 according to Erisman et al. 1994
    !       R_STOM
    rsmx = rs*dratio(NH3)
    !       GRASSLAND (PASTURE DURING GRAZING)
!kml    IF (ixxxlu(iland)==3) THEN
    IF (ixxxlu(iland)==2 .OR. ixxxlu(iland)==3 .OR. ixxxlu(iland)==4) THEN
       IF (iseason==1) THEN
          !           SUMMER
          rcx(NH3) = 1000.
       END IF
       IF ((iseason==2) .OR. (iseason==3) .OR. (iseason==5)) THEN
          !           WINTER, NO SNOW
          IF (tc>-1.) THEN
             IF (rad/=0.) THEN
                rcx(NH3) = 50.
             ELSE
                rcx(NH3) = 100.
             END IF
             IF ((wetflag) .OR. (rainflag)) THEN
                rcx(NH3) = 20.
             END IF
          END IF
          IF ((tc>=(-5.)) .AND. (tc<=-1.)) THEN
             rcx(NH3) = 200.
          END IF
          IF (tc<(-5.)) THEN
             rcx(NH3) = 500.
          END IF
       END IF
    END IF
    !       AGRICULTURAL LAND (CROPS AND UNGRAZED PASTURE)
!kml    IF (ixxxlu(iland)==2) THEN
    IF (ixxxlu(iland)==5 .OR. ixxxlu(iland)==6) THEN

       IF (iseason==1) THEN
          !           SUMMER
          IF (rad/=0.) THEN
             rcx(NH3) = rsmx
          ELSE
             rcx(NH3) = 200.
          END IF
          IF ((wetflag) .OR. (rainflag)) THEN
             rcx(NH3) = 50.
          END IF
       END IF
       IF ((iseason==2) .OR. (iseason==3) .OR. (iseason==5)) THEN
          !           WINTER, NO SNOW
          IF (tc>-1.) THEN
             IF (rad/=0.) THEN
                rcx(NH3) = rsmx
             ELSE
                rcx(NH3) = 300.
             END IF
             IF ((wetflag) .OR. (rainflag)) THEN
                rcx(NH3) = 100.
             END IF
          END IF
          IF ((tc>=(-5.)) .AND. (tc<=-1.)) THEN
             rcx(NH3) = 200.
          END IF
          IF (tc<(-5.)) THEN
             rcx(NH3) = 500.
          END IF
       END IF
    END IF
    !       SEMI-NATURAL ECOSYSTEMS AND FORESTS
!kml    IF ((ixxxlu(iland)==4) .OR. (ixxxlu(iland)==5) .OR. (ixxxlu( &
!kml         iland)==6)) THEN

    IF ((ixxxlu(iland)==10) .OR. (ixxxlu(iland)==11) .OR.  &
        (ixxxlu(iland)==12) .OR. (ixxxlu(iland)==13) .OR.  &
        (ixxxlu(iland)==14) .OR. (ixxxlu(iland)==15) .OR.  &
        (ixxxlu(iland)==18)) THEN

       IF (rad/=0.) THEN
          rcx(NH3) = 500.
       ELSE
          rcx(NH3) = 1000.
       END IF
       IF ((wetflag) .OR. (rainflag)) THEN
          IF (highnh3) THEN
             rcx(NH3) = 100.
          ELSE
             rcx(NH3) = 0.
          END IF
       END IF
       IF ((iseason==2) .OR. (iseason==3) .OR. (iseason==5)) THEN
          !           WINTER, NO SNOW
          IF ((tc>=(-5.)) .AND. (tc<=-1.)) THEN
             rcx(NH3) = 200.
          END IF
          IF (tc<(-5.)) THEN
             rcx(NH3) = 500.
          END IF
       END IF
    END IF
    !       WATER
    IF (iland==iswater) THEN
       rcx(NH3) = 0.
    END IF
    !       URBAN AND DESERT (SOIL SURFACES)
!kml    IF (ixxxlu(iland)==1) THEN
    IF (ixxxlu(iland)==1 .OR. ixxxlu(iland)==19 ) THEN

       IF ( .NOT. wetflag) THEN
          rcx(NH3) = 50.
       ELSE
          rcx(NH3) = 0.
       END IF
    END IF
    !       SNOW COVERED SURFACES OR PERMANENT ICE
    IF ((iseason==4) .OR. (iland==isice)) THEN
       IF (tc>2.) THEN
          rcx(NH3) = 0.
       END IF
       IF ((tc>=(-1.)) .AND. (tc<=2.)) THEN
          rcx(NH3) = 70.*(2.-tc)
       END IF
       IF (tc<(-1.)) THEN
          rcx(NH3) = 500.
       END IF
    END IF
    IF (rcx(NH3)<1.) rcx(NH3) = 1.

  END SUBROUTINE rc
  ! **********************************************************************
  SUBROUTINE cellvg(vgtemp,ustar,dz,zr,rmol,nspec)
    !     THIS PROGRAM HAS BEEN DESIGNED TO CALCULATE THE CELL AVERAGE
    !     DEPOSITION VELOCITY GIVEN THE VALUE OF VG AT SOME REFERENCE
    !     HEIGHT ZR WHICH IS MUCH SMALLER THAN THE CELL HEIGHT DZ.
    !       PROGRAM WRITTEN BY GREGORY J.MCRAE (NOVEMBER 1977)
    !         Modified by Darrell A. Winner    (February 1991)
    !.....PROGRAM VARIABLES...
    !     VgTemp   - DEPOSITION VELOCITY AT THE REFERENCE HEIGHT
    !     USTAR    - FRICTION VELOCITY
    !     RMOL     - RECIPROCAL OF THE MONIN-OBUKHOV LENGTH
    !     ZR       - REFERENCE HEIGHT
    !     DZ       - CELL HEIGHT
    !     CELLVG   - CELL AVERAGE DEPOSITION VELOCITY
    !     VK       - VON KARMAN CONSTANT

    REAL    , INTENT(IN)  :: ustar
    REAL    , INTENT(IN)  :: dz
    REAL    , INTENT(IN)  :: zr
    REAL    , INTENT(IN)  :: rmol
    INTEGER , INTENT(IN)  :: nspec
    REAL    , INTENT(OUT) :: vgtemp(nspec)

    ! ..
    ! .. Local Scalars ..
    REAL :: a, fac, pdz, pzr,fx
    INTEGER :: nss
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC alog, sqrt
    ! ..
    !     Set the von Karman constant
    REAL, PARAMETER :: vk = 0.4

    !     DETERMINE THE STABILITY BASED ON THE CONDITIONS
    !             1/L < 0 UNSTABLE
    !             1/L = 0 NEUTRAL
    !             1/L > 0 STABLE


    IF (rmol<0) THEN
       pdz = SQRT(1.0-9.0*dz*rmol)
       pzr = SQRT(1.0-9.0*zr*rmol)
       fac = ((pdz-1.0)/(pzr-1.0))*((pzr+1.0)/(pdz+1.0))
       a = 0.74*dz*alog(fac) + (0.164/rmol)*(pdz-pzr)
    ELSE IF (rmol==0) THEN
       a = 0.74*(dz*alog(dz/zr)-dz+zr)
    ELSE
       a = 0.74*(dz*alog(dz/zr)-dz+zr) + (2.35*rmol)*(dz-zr)**2
    END IF

    fx = a/(vk*ustar*(dz-zr))

    !      CALCULATE THE DEPOSITION VELOCITIY
    DO nss = 1, nspec
       vgtemp(nss) = vgtemp(nss)/(1.0+vgtemp(nss)*fx)
    END DO

    RETURN
  END SUBROUTINE cellvg


  !========================================================================

  SUBROUTINE apply_drydep(m1,m2,m3,ia,iz,ja,jz,V_dep,sclt,M_dep,sclp &
                         ,dens,rtgt,dzt,dt)

    INTEGER , INTENT(IN)    :: m1
    INTEGER , INTENT(IN)    :: m2
    INTEGER , INTENT(IN)    :: m3
    INTEGER , INTENT(IN)    :: ia
    INTEGER , INTENT(IN)    :: iz
    INTEGER , INTENT(IN)    :: ja
    INTEGER , INTENT(IN)    :: jz

    REAL    , INTENT(IN)    :: V_dep(m2,m3)   ! dry deposition velocity (m/s)
    REAL    , INTENT(INOUT) :: sclt(m1,m2,m3) ! tendency (kg[tracer] kg[air]^-1 s^-1)
    REAL    , INTENT(INOUT) :: M_dep(m2,m3)   ! accumulated mass on surface due dry dep
                              ! process (kg m^-2)
    REAL    , INTENT(IN)    :: sclp(m1,m2,m3) ! tracer mixing ratio (kg/kg)
    REAL    , INTENT(IN)    :: dens(m2,m3)
    REAL    , INTENT(IN)    :: rtgt(m2,m3) ! (DMK) not used
    REAL    , INTENT(IN)    :: dzt(m1)

    REAL    , INTENT(IN)    :: dt

    INTEGER :: i,j
    REAL    :: dz,tend_drydep

    DO j=ja,jz
       DO i=ia,iz
          !-1st vertical layer thickness
          dz = rtgt(i,j)/dzt(2) ! dzt=1/(z(k)-z(k-1))
          !- tendency to dry deposition process
          tend_drydep  = - V_dep(i,j)*sclp(2,i,j)/(dz)
          !- update total tendency kg[tracer]/kg[air]/s
          sclt(2,i,j) = sclt(2,i,j) + tend_drydep
          !- accumulate the surface deposited mass of the tracer by this process
          !	M_dep(i,j) = M_dep(i,j) + dt*tend_drydep*dens(i,j) ! kg[tracer] m^-2
          M_dep(i,j) = M_dep(i,j) - dt*tend_drydep*dens(i,j)*dz ! kg[tracer] m^-2
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE apply_drydep
  !========================================================================

  SUBROUTINE define_PBL_height(m1,m2,m3,npatch,ia,iz,ja,jz,zt,tke,shf,rcp,rtgt,Zi)

    INTEGER , INTENT(IN)    :: m1
    INTEGER , INTENT(IN)    :: m2
    INTEGER , INTENT(IN)    :: m3
    INTEGER , INTENT(IN)    :: npatch ! (DMK) not used
    INTEGER , INTENT(IN)    :: ia
    INTEGER , INTENT(IN)    :: iz
    INTEGER , INTENT(IN)    :: ja
    INTEGER , INTENT(IN)    :: jz

    REAL    , INTENT(IN)    :: zt(m1)
    REAL    , INTENT(IN)    :: tke(m1,m2,m3)
    REAL    , INTENT(IN)    :: shf(m2,m3)
    REAL    , INTENT(IN)    :: rcp(m1,m2,m3)
    REAL    , INTENT(IN)    :: rtgt(m2,m3)
    REAL    , INTENT(INOUT) :: Zi(m2,m3)

    REAL,PARAMETER :: tkethrsh=0.001       !   tke threshold for PBL height in m2/s2
    !   tkmin    = 5.e-4   minimum TKE in RAMS
    REAL,PARAMETER :: rcmin=1.e-4          !   min liq water = 0.1 g/kg

    REAL pblht
    INTEGER :: i,j,k


    DO j=ja,jz
       DO i=ia,iz

          Zi(i,j) = 0.
          !- convective layer
          IF(shf(i,j) >= 1.e-8) THEN

             pblht=0.
             DO k=2,m1-1
                pblht=zt(k)*rtgt(i,j)
                !                  if(i.ge.10.and.i.le.25.and.j.ge.13.and.j.le.25)
                !     &               print*,'i,j,k,z,pbl=',i,j,k,ztn(k,ngrd),pblht
                IF( rcp(k,i,j) .GT. rcmin     )GOTO 10 ! dry convective layer  ! (DMK) EXIT
                IF( tke(k,i,j) .LE. tkethrsh  )GOTO 10                         ! (DMK) EXIT
             ENDDO
10           CONTINUE

             Zi(i,j)=pblht
          ENDIF
          !print*,'PBLh',Zi(i,j)

       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE define_PBL_height


  !========================================================================
  SUBROUTINE dep_init(nspecies_chem,dvj)

    INTEGER , INTENT(IN) :: nspecies_chem
    REAL    , INTENT(IN) :: dvj(nspecies_chem)

    ! ..
    ! ..
    ! .. Local Scalars ..
    REAL :: sc
    INTEGER :: iland, iseason, l
    INTEGER :: iprt
    ! ..
    ! .. Local Arrays ..
    REAL :: dat1(nlu,dep_seasons), dat2(nlu,dep_seasons),         &
         dat3(nlu,dep_seasons), dat4(nlu,dep_seasons),         &
         dat5(nlu,dep_seasons), dat6(nlu,dep_seasons),         &
         dat7(nlu,dep_seasons)
    ! ..
    ! .. Data Statements ..
    !     RI for stomatal resistance
    !      data ((ri(ILAND,ISEASON),ILAND=1,nlu),ISEASON=1,dep_seasons)/0.10E+11, &
    DATA ((dat1(iland,iseason),iland=1,nlu),iseason=1,dep_seasons)/0.10E+11, &
         0.60E+02, 0.60E+02, 0.60E+02, 0.60E+02, 0.70E+02, 0.12E+03, &
         0.12E+03, 0.12E+03, 0.12E+03, 0.70E+02, 0.13E+03, 0.70E+02, &
         0.13E+03, 0.10E+03, 0.10E+11, 0.80E+02, 0.10E+03, 0.10E+11, &
         0.80E+02, 0.10E+03, 0.10E+03, 0.10E+11, 0.10E+11, 0.10E+11, &
         0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
         0.10E+11, 0.10E+11, 0.10E+11, 0.12E+03, 0.10E+11, 0.10E+11, &
         0.70E+02, 0.25E+03, 0.50E+03, 0.10E+11, 0.10E+11, 0.50E+03, &
         0.10E+11, 0.10E+11, 0.50E+03, 0.50E+03, 0.10E+11, 0.10E+11, &
         0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
         0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.12E+03, 0.10E+11, &
         0.10E+11, 0.70E+02, 0.25E+03, 0.50E+03, 0.10E+11, 0.10E+11, &
         0.50E+03, 0.10E+11, 0.10E+11, 0.50E+03, 0.50E+03, 0.10E+11, &
         0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
         0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
         0.10E+11, 0.10E+11, 0.70E+02, 0.40E+03, 0.80E+03, 0.10E+11, &
         0.10E+11, 0.80E+03, 0.10E+11, 0.10E+11, 0.80E+03, 0.80E+03, &
         0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.12E+03, 0.12E+03, &
         0.12E+03, 0.12E+03, 0.14E+03, 0.24E+03, 0.24E+03, 0.24E+03, &
         0.12E+03, 0.14E+03, 0.25E+03, 0.70E+02, 0.25E+03, 0.19E+03, &
         0.10E+11, 0.16E+03, 0.19E+03, 0.10E+11, 0.16E+03, 0.19E+03, &
         0.19E+03, 0.10E+11, 0.10E+11, 0.10E+11/
    ! ..
    IF (nlu/=25) THEN
       PRINT *, 'number of land use classifications not correct '
       STOP
    END IF
    IF (dep_seasons/=5) THEN
       PRINT *, 'number of dep_seasons not correct '
       STOP
    END IF

    !     SURFACE RESISTANCE DATA FOR DEPOSITION MODEL OF
    !     M. L. WESELY, ATMOSPHERIC ENVIRONMENT 23 (1989) 1293-1304

    !     Seasonal categories:
    !     1: midsummer with lush vegetation
    !     2: autumn with unharvested cropland
    !     3: late autumn with frost, no snow
    !     4: winter, snow on ground and subfreezing
    !     5: transitional spring with partially green short annuals

    !     Land use types:
    !     USGS type                                Wesely type
    !      1: Urban and built-up land              1  Urban land
    !      2: Dryland cropland and pasture         2  agricultural land
    !      3: Irrigated cropland and pasture       2
    !      4: Mix. dry/irrg. cropland and pasture  2
    !      5: Cropland/grassland mosaic            2
    !      6: Cropland/woodland mosaic             4  deciduous forest
    !      7: Grassland                            3  range land
    !      8: Shrubland                            3
    !      9: Mixed shrubland/grassland            3
    !     10: Savanna                              3, range land:    always summer
    !     11: Deciduous broadleaf forest           4
    !     12: Deciduous needleleaf forest          5, coniferous forest:autumn and winter modi
    !     13: Evergreen broadleaf forest           4, always summer
    !     14: Evergreen needleleaf forest          5
    !     15: Mixed Forest                         6 Mixed forest
    !     16: Water Bodies                         7 Water
    !     17: Herbaceous wetland                   9
    !     18: Wooded wetland                       6
    !     19: Barren or sparsely vegetated         8 Barren land, most desert
    !     20: Herbaceous Tundra                    9 Nonforested wetland
    !     21: Wooded Tundra                        6
    !     22: Mixed Tundra                         6
    !     23: Bare Ground Tundra                   8
    !     24: Snow or Ice                          -, always winter
    !     25: No data                              8


    !     Order of data:
    !      |
    !      |   seasonal category
    !     \|/
    !     ---> landuse USGS type  !kml
    !     1       2       3       4       5       6      ..... 25    !kml
    !     RLU for outer surfaces in the upper canopy
    DO iseason = 1, dep_seasons
       DO iland = 1, nlu

!--(DMK-CCATT-INI)-------------------------------------------------------------
          ri((iland),iseason) = dat1(iland,iseason)
!--(DMK-CCATT-OLD)-------------------------------------------------------------
!          ri(iland,iseason) = dat1(iland,iseason)
!--(DMK-CCATT-FIM)-------------------------------------------------------------

       END DO
    END DO
    !      data ((rlu(ILAND,ISEASON),ILAND=1,25),ISEASON=1,5)/0.10E+11, &
    DATA ((dat2(iland,iseason),iland=1,nlu),iseason=1,dep_seasons)/0.10E+11, &
         0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, &
         0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, &
         0.20E+04, 0.20E+04, 0.10E+11, 0.25E+04, 0.20E+04, 0.10E+11, &
         0.25E+04, 0.20E+04, 0.20E+04, 0.10E+11, 0.10E+11, 0.10E+11, &
         0.10E+11, 0.90E+04, 0.90E+04, 0.90E+04, 0.90E+04, 0.90E+04, &
         0.90E+04, 0.90E+04, 0.90E+04, 0.20E+04, 0.90E+04, 0.90E+04, &
         0.20E+04, 0.40E+04, 0.80E+04, 0.10E+11, 0.90E+04, 0.80E+04, &
         0.10E+11, 0.90E+04, 0.80E+04, 0.80E+04, 0.10E+11, 0.10E+11, &
         0.10E+11, 0.10E+11, 0.90E+04, 0.90E+04, 0.90E+04, 0.90E+04, &
         0.90E+04, 0.90E+04, 0.90E+04, 0.90E+04, 0.20E+04, 0.90E+04, &
         0.90E+04, 0.20E+04, 0.40E+04, 0.80E+04, 0.10E+11, 0.90E+04, &
         0.80E+04, 0.10E+11, 0.90E+04, 0.80E+04, 0.80E+04, 0.10E+11, &
         0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
         0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
         0.10E+11, 0.10E+11, 0.20E+04, 0.60E+04, 0.90E+04, 0.10E+11, &
         0.90E+04, 0.90E+04, 0.10E+11, 0.90E+04, 0.90E+04, 0.90E+04, &
         0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.40E+04, 0.40E+04, &
         0.40E+04, 0.40E+04, 0.40E+04, 0.40E+04, 0.40E+04, 0.40E+04, &
         0.20E+04, 0.40E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.30E+04, &
         0.10E+11, 0.40E+04, 0.30E+04, 0.10E+11, 0.40E+04, 0.30E+04, &
         0.30E+04, 0.10E+11, 0.10E+11, 0.10E+11/
    DO iseason = 1, dep_seasons
       DO iland = 1, nlu

!--(DMK-CCATT-INI)-------------------------------------------------------------
          rlu((iland),iseason) = dat2(iland,iseason)
!--(DMK-CCATT-OLD)-------------------------------------------------------------
!          rlu(iland,iseason) = dat2(iland,iseason)
!--(DMK-CCATT-FIM)-------------------------------------------------------------

       END DO
    END DO
    !     RAC for transfer that depends on canopy height and density
    !      data ((rac(ILAND,ISEASON),ILAND=1,25),ISEASON=1,5)/0.10E+03, &
    DATA ((dat3(iland,iseason),iland=1,nlu),iseason=1,dep_seasons)/0.10E+03, &
         0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+04, 0.10E+03, &
         0.10E+03, 0.10E+03, 0.10E+03, 0.20E+04, 0.20E+04, 0.20E+04, &
         0.20E+04, 0.20E+04, 0.00E+00, 0.30E+03, 0.20E+04, 0.00E+00, &
         0.30E+03, 0.20E+04, 0.20E+04, 0.00E+00, 0.00E+00, 0.00E+00, &
         0.10E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.15E+04, &
         0.10E+03, 0.10E+03, 0.10E+03, 0.10E+03, 0.15E+04, 0.20E+04, &
         0.20E+04, 0.20E+04, 0.17E+04, 0.00E+00, 0.20E+03, 0.17E+04, &
         0.00E+00, 0.20E+03, 0.17E+04, 0.17E+04, 0.00E+00, 0.00E+00, &
         0.00E+00, 0.10E+03, 0.10E+02, 0.10E+02, 0.10E+02, 0.10E+02, &
         0.10E+04, 0.10E+03, 0.10E+03, 0.10E+03, 0.10E+03, 0.10E+04, &
         0.20E+04, 0.20E+04, 0.20E+04, 0.15E+04, 0.00E+00, 0.10E+03, &
         0.15E+04, 0.00E+00, 0.10E+03, 0.15E+04, 0.15E+04, 0.00E+00, &
         0.00E+00, 0.00E+00, 0.10E+03, 0.10E+02, 0.10E+02, 0.10E+02, &
         0.10E+02, 0.10E+04, 0.10E+02, 0.10E+02, 0.10E+02, 0.10E+02, &
         0.10E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.15E+04, 0.00E+00, &
         0.50E+02, 0.15E+04, 0.00E+00, 0.50E+02, 0.15E+04, 0.15E+04, &
         0.00E+00, 0.00E+00, 0.00E+00, 0.10E+03, 0.50E+02, 0.50E+02, &
         0.50E+02, 0.50E+02, 0.12E+04, 0.80E+02, 0.80E+02, 0.80E+02, &
         0.10E+03, 0.12E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.15E+04, &
         0.00E+00, 0.20E+03, 0.15E+04, 0.00E+00, 0.20E+03, 0.15E+04, &
         0.15E+04, 0.00E+00, 0.00E+00, 0.00E+00/
    DO iseason = 1, dep_seasons
       DO iland = 1, nlu

!--(DMK-CCATT-INI)-------------------------------------------------------------
         rac((iland),iseason) = dat3(iland,iseason)
!--(DMK-CCATT-OLD)-------------------------------------------------------------
!          rac(iland,iseason) = dat3(iland,iseason)
!--(DMK-CCATT-FIM)-------------------------------------------------------------

       END DO
    END DO
    !     RGSS for ground surface  SO2
    !      data ((rgss(ILAND,ISEASON),ILAND=1,25),ISEASON=1,5)/0.40E+03, &
    DATA ((dat4(iland,iseason),iland=1,nlu),iseason=1,dep_seasons)/0.40E+03, &
         0.15E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.50E+03, 0.35E+03, &
         0.35E+03, 0.35E+03, 0.35E+03, 0.50E+03, 0.50E+03, 0.50E+03, &
         0.50E+03, 0.10E+03, 0.10E+01, 0.10E+01, 0.10E+03, 0.10E+04, &
         0.10E+01, 0.10E+03, 0.10E+03, 0.10E+04, 0.10E+03, 0.10E+04, &
         0.40E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.50E+03, &
         0.35E+03, 0.35E+03, 0.35E+03, 0.35E+03, 0.50E+03, 0.50E+03, &
         0.50E+03, 0.50E+03, 0.10E+03, 0.10E+01, 0.10E+01, 0.10E+03, &
         0.10E+04, 0.10E+01, 0.10E+03, 0.10E+03, 0.10E+04, 0.10E+03, &
         0.10E+04, 0.40E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.15E+03, &
         0.50E+03, 0.35E+03, 0.35E+03, 0.35E+03, 0.35E+03, 0.50E+03, &
         0.50E+03, 0.50E+03, 0.50E+03, 0.20E+03, 0.10E+01, 0.10E+01, &
         0.20E+03, 0.10E+04, 0.10E+01, 0.20E+03, 0.20E+03, 0.10E+04, &
         0.10E+03, 0.10E+04, 0.10E+03, 0.10E+03, 0.10E+03, 0.10E+03, &
         0.10E+03, 0.10E+03, 0.10E+03, 0.10E+03, 0.10E+03, 0.10E+03, &
         0.10E+03, 0.10E+03, 0.50E+03, 0.10E+03, 0.10E+03, 0.10E+01, &
         0.10E+03, 0.10E+03, 0.10E+04, 0.10E+03, 0.10E+03, 0.10E+03, &
         0.10E+04, 0.10E+03, 0.10E+04, 0.50E+03, 0.15E+03, 0.15E+03, &
         0.15E+03, 0.15E+03, 0.50E+03, 0.35E+03, 0.35E+03, 0.35E+03, &
         0.35E+03, 0.50E+03, 0.50E+03, 0.50E+03, 0.50E+03, 0.20E+03, &
         0.10E+01, 0.10E+01, 0.20E+03, 0.10E+04, 0.10E+01, 0.20E+03, &
         0.20E+03, 0.10E+04, 0.10E+03, 0.10E+04/
    DO iseason = 1, dep_seasons
       DO iland = 1, nlu

!--(DMK-CCATT-INI)-------------------------------------------------------------
         rgss((iland),iseason) = dat4(iland,iseason)
!--(DMK-CCATT-OLD)-------------------------------------------------------------
!          rgss(iland,iseason) = dat4(iland,iseason)
!--(DMK-CCATT-FIM)-------------------------------------------------------------

       END DO
    END DO
    !     RGSO for ground surface  O3
    !      data ((rgso(ILAND,ISEASON),ILAND=1,25),ISEASON=1,5)/0.30E+03, &
    DATA ((dat5(iland,iseason),iland=1,nlu),iseason=1,dep_seasons)/0.30E+03, &
         0.15E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.20E+03, 0.20E+03, &
         0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, &
         0.20E+03, 0.30E+03, 0.20E+04, 0.10E+04, 0.30E+03, 0.40E+03, &
         0.10E+04, 0.30E+03, 0.30E+03, 0.40E+03, 0.35E+04, 0.40E+03, &
         0.30E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.20E+03, &
         0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, &
         0.20E+03, 0.20E+03, 0.30E+03, 0.20E+04, 0.80E+03, 0.30E+03, &
         0.40E+03, 0.80E+03, 0.30E+03, 0.30E+03, 0.40E+03, 0.35E+04, &
         0.40E+03, 0.30E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.15E+03, &
         0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, &
         0.20E+03, 0.20E+03, 0.20E+03, 0.30E+03, 0.20E+04, 0.10E+04, &
         0.30E+03, 0.40E+03, 0.10E+04, 0.30E+03, 0.30E+03, 0.40E+03, &
         0.35E+04, 0.40E+03, 0.60E+03, 0.35E+04, 0.35E+04, 0.35E+04, &
         0.35E+04, 0.35E+04, 0.35E+04, 0.35E+04, 0.35E+04, 0.35E+04, &
         0.35E+04, 0.35E+04, 0.20E+03, 0.35E+04, 0.35E+04, 0.20E+04, &
         0.35E+04, 0.35E+04, 0.40E+03, 0.35E+04, 0.35E+04, 0.35E+04, &
         0.40E+03, 0.35E+04, 0.40E+03, 0.30E+03, 0.15E+03, 0.15E+03, &
         0.15E+03, 0.15E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, &
         0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.30E+03, &
         0.20E+04, 0.10E+04, 0.30E+03, 0.40E+03, 0.10E+04, 0.30E+03, &
         0.30E+03, 0.40E+03, 0.35E+04, 0.40E+03/
    DO iseason = 1, dep_seasons
       DO iland = 1, nlu

!--(DMK-CCATT-INI)-------------------------------------------------------------
          rgso((iland),iseason) = dat5(iland,iseason)
!--(DMK-CCATT-OLD)-------------------------------------------------------------
!          rgso(iland,iseason) = dat5(iland,iseason)
!--(DMK-CCATT-FIM)-------------------------------------------------------------

       END DO
    END DO
    !     RCLS for exposed surfaces in the lower canopy  SO2
    !      data ((rcls(ILAND,ISEASON),ILAND=1,25),ISEASON=1,5)/0.10E+11, &
    DATA ((dat6(iland,iseason),iland=1,nlu),iseason=1,dep_seasons)/0.10E+11, &
         0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, &
         0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, &
         0.20E+04, 0.20E+04, 0.10E+11, 0.25E+04, 0.20E+04, 0.10E+11, &
         0.25E+04, 0.20E+04, 0.20E+04, 0.10E+11, 0.10E+11, 0.10E+11, &
         0.10E+11, 0.90E+04, 0.90E+04, 0.90E+04, 0.90E+04, 0.90E+04, &
         0.90E+04, 0.90E+04, 0.90E+04, 0.20E+04, 0.90E+04, 0.90E+04, &
         0.20E+04, 0.20E+04, 0.40E+04, 0.10E+11, 0.90E+04, 0.40E+04, &
         0.10E+11, 0.90E+04, 0.40E+04, 0.40E+04, 0.10E+11, 0.10E+11, &
         0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
         0.90E+04, 0.90E+04, 0.90E+04, 0.90E+04, 0.20E+04, 0.90E+04, &
         0.90E+04, 0.20E+04, 0.30E+04, 0.60E+04, 0.10E+11, 0.90E+04, &
         0.60E+04, 0.10E+11, 0.90E+04, 0.60E+04, 0.60E+04, 0.10E+11, &
         0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
         0.10E+11, 0.90E+04, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
         0.90E+04, 0.90E+04, 0.20E+04, 0.20E+03, 0.40E+03, 0.10E+11, &
         0.90E+04, 0.40E+03, 0.10E+11, 0.90E+04, 0.40E+03, 0.40E+03, &
         0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.40E+04, 0.40E+04, &
         0.40E+04, 0.40E+04, 0.40E+04, 0.40E+04, 0.40E+04, 0.40E+04, &
         0.20E+04, 0.40E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.30E+04, &
         0.10E+11, 0.40E+04, 0.30E+04, 0.10E+11, 0.40E+04, 0.30E+04, &
         0.30E+04, 0.10E+11, 0.10E+11, 0.10E+11/
    DO iseason = 1, dep_seasons
       DO iland = 1, nlu
          rcls(iland,iseason) = dat6(iland,iseason)
       END DO
    END DO
    !     RCLO for exposed surfaces in the lower canopy  O3
    !      data ((rclo(ILAND,ISEASON),ILAND=1,25),ISEASON=1,5)/0.10E+11, &
    DATA ((dat7(iland,iseason),iland=1,nlu),iseason=1,dep_seasons)/0.10E+11, &
         0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, &
         0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, &
         0.10E+04, 0.10E+04, 0.10E+11, 0.10E+04, 0.10E+04, 0.10E+11, &
         0.10E+04, 0.10E+04, 0.10E+04, 0.10E+11, 0.10E+11, 0.10E+11, &
         0.10E+11, 0.40E+03, 0.40E+03, 0.40E+03, 0.40E+03, 0.40E+03, &
         0.40E+03, 0.40E+03, 0.40E+03, 0.10E+04, 0.40E+03, 0.40E+03, &
         0.10E+04, 0.10E+04, 0.60E+03, 0.10E+11, 0.40E+03, 0.60E+03, &
         0.10E+11, 0.40E+03, 0.60E+03, 0.60E+03, 0.10E+11, 0.10E+11, &
         0.10E+11, 0.10E+11, 0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, &
         0.40E+03, 0.40E+03, 0.40E+03, 0.40E+03, 0.10E+04, 0.40E+03, &
         0.40E+03, 0.10E+04, 0.10E+04, 0.60E+03, 0.10E+11, 0.80E+03, &
         0.60E+03, 0.10E+11, 0.80E+03, 0.60E+03, 0.60E+03, 0.10E+11, &
         0.10E+11, 0.10E+11, 0.10E+11, 0.10E+04, 0.10E+04, 0.10E+04, &
         0.10E+04, 0.40E+03, 0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, &
         0.40E+03, 0.40E+03, 0.10E+04, 0.15E+04, 0.60E+03, 0.10E+11, &
         0.80E+03, 0.60E+03, 0.10E+11, 0.80E+03, 0.60E+03, 0.60E+03, &
         0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+04, 0.10E+04, &
         0.10E+04, 0.10E+04, 0.50E+03, 0.50E+03, 0.50E+03, 0.50E+03, &
         0.10E+04, 0.50E+03, 0.15E+04, 0.10E+04, 0.15E+04, 0.70E+03, &
         0.10E+11, 0.60E+03, 0.70E+03, 0.10E+11, 0.60E+03, 0.70E+03, &
         0.70E+03, 0.10E+11, 0.10E+11, 0.10E+11/
    DO iseason = 1, dep_seasons
       DO iland = 1, nlu

!--(DMK-CCATT-INI)-------------------------------------------------------------
          rclo((iland),iseason) = dat7(iland,iseason)
!--(DMK-CCATT-OLD)-------------------------------------------------------------
!          rclo(iland,iseason) = dat7(iland,iseason)
!--(DMK-CCATT-FIM)-------------------------------------------------------------

       END DO
    END DO

    DO l = 1, nspecies_chem
       ! if(dvj(l) < 0.) dvj(l) = 1./(sqrt(weight(l)))
       !          hstar4(l) = hstar(l) ! preliminary
       ! Correction of diff. coeff
       !          dvj(l) = dvj(l)*(293.15/298.15)**1.75
       !          sc = 0.15/dvj(l) ! Schmidt Number at 20�C
       dratio(l) = 0.242/dvj(l) !    ! of water vapor and gas at
       ! Ratio of diffusion coeffi
       !          scpr23(l) = (sc/0.72)**(2./3.) ! (Schmidt # / Prandtl #)**
    END DO


    !     DATA FOR AEROSOL PARTICLE DEPOSITION FOR THE MODEL OF
    !     J. W. ERISMAN, A. VAN PUL AND P. WYERS
    !     ATMOSPHERIC ENVIRONMENT 28 (1994), 2595-2607

    !     vd = (u* / k) * CORRECTION FACTORS

    !     CONSTANT K FOR LANDUSE TYPES:
    !srf- general initialization
    kpart(:) = 500.
    ! urban and built-up land
    kpart(19) = 500.
    ! dryland cropland and pasture
    kpart(2) = 500.
    ! irrigated cropland and pasture
    kpart(16) = 500.
    ! mixed dryland/irrigated cropland and past
    kpart(16) = 500.
    ! cropland/grassland mosaic
    kpart(15) = 500.
    ! cropland/woodland mosaic
    kpart(14) = 100.
    ! grassland
    kpart(8) = 500.
    ! shrubland
    kpart(13) = 500.
    ! mixed shrubland/grassland
    kpart(9) = 500.
    ! savanna
    kpart(9) = 500.
    ! deciduous broadleaf forest
    kpart(5) = 100.
    ! deciduous needleleaf forest
    kpart(6) = 100.
    ! evergreen broadleaf forest
    kpart(7) = 100.
    ! evergreen needleleaf forest
    kpart(4) = 100.
    ! mixed forest
    kpart(14) = 100.
    ! water bodies
    kpart(1) = 500.
    ! herbaceous wetland
    kpart(17) = 500.
    ! wooded wetland
    kpart(18) = 500.
    ! barren or sparsely vegetated
    kpart(10) = 500.
    ! herbaceous tundra
    kpart(11) = 500.
    ! wooded tundra
    kpart(11) = 100.
    ! mixed tundra
    kpart(11) = 500.
    ! bare ground tundra
    kpart(3) = 500.
    ! snow or ice
    kpart(2) = 500.
    !     Comments:
    kpart(25) = 500.
    !     Erisman et al. (1994) give
    !     k = 500 for low vegetation and k = 100 for forests.

    !     For desert k = 500 is taken according to measurements
    !     on bare soil by
    !     J. Fontan, A. Lopez, E. Lamaud and A. Druilhet (1997)
    !     "Vertical Flux Measurements of the Submicronic Aerosol Particles
    !     and Parametrisation of the Dry Deposition Velocity"
    !     in: "Biosphere-Atmosphere Exchange of Pollutants
    !     and Trace Substances"
    !     Editor: S. Slanina. Springer-Verlag Berlin, Heidelberg, 1997
    !     pp. 381-390

    !     For coniferous forest the Erisman value of  k = 100 is taken.
    !     Measurements of Erisman et al. (1997) in a coniferous forest
    !     in the Netherlands, lead to values of k between 20 and 38
    !     (Atmospheric Environment 31 (1997), 321-332).
    !     However, these high values of vd may be reached during
    !     instable cases. The eddy correlation measurements
    !     of Gallagher et al. (1997) made during the same experiment
    !     show for stable cases (L>0) values of k between 200 and 250
    !     at minimum (Atmospheric Environment 31 (1997), 359-373).
    !     Fontan et al. (1997) found k = 250 in a forest
    !     of maritime pine in southwestern France.

    !     For gras, model calculations of Davidson et al. support
    !     the value of 500.
    !     C. I. Davidson, J. M. Miller and M. A. Pleskov
    !     "The Influence of Surface Structure on Predicted Particles
    !     Dry Deposition to Natural Gras Canopies"
    !     Water, Air, and Soil Pollution 18 (1982) 25-43

    !     Snow covered surface: The experiment of Ibrahim et al. (1983)
    !     gives k = 436 for 0.7 um diameter particles.
    !     The deposition velocity of Milford and Davidson (1987)
    !     gives k = 154 for continental sulfate aerosol.
    !     M. Ibrahim, L. A. Barrie and F. Fanaki
    !     Atmospheric Environment 17 (1983), 781-788

    !     J. B. Milford and C. I. Davidson
    !     "The Sizes of Particulate Sulfate and Nitrate in the Atmosphere
    !     - A Review"
    !     JAPCA 37 (1987), 125-134
    !
    !------------!-------------!----------------------------------------------------------------
    !  LEAF-3 CLASS (20)                      !  Wesely type !   Land use types:
    !---AND DESCRIPTION 			  !		 !   USGS type
    !  0  Ocean                               !   7 	 !    1: Urban and built-up land
    !  1  Lakes, rivers, streams              !   7 	 !    2: Dryland cropland and pasture
    !  2  Ice cap/glacier                     !   - 	 !    3: Irrigated cropland and pasture
    !  3  Desert, bare soil                   !   8 	 !    4: Mix. dry/irrg. cropland and pasture
    !  4  Evergreen needleleaf tree           !   5 	 !    5: Cropland/grassland mosaic
    !  5  Deciduous needleleaf tree           !   5 	 !    6: Cropland/woodland mosaic
    !  6  Deciduous broadleaf tree            !   4 	 !    7: Grassland
    !  7  Evergreen broadleaf tree            !   4 	 !    8: Shrubland
    !  8  Short grass                         !   2 	 !    9: Mixed shrubland/grassland
    !  9  Tall grass                          !   3 	 !   10: Savanna
    ! 10  Semi-desert                         !   8 	 !   11: Deciduous broadleaf forest
    ! 11  Tundra                              !   8 	 !   12: Deciduous needleleaf forest
    ! 12  Evergreen shrub                     !   3 	 !   13: Evergreen broadleaf forest
    ! 13  Deciduous shrub                     !   3 	 !   14: Evergreen needleleaf forest
    ! 14  Mixed woodland                      !   6 	 !   15: Mixed Forest
    ! 15  Crop/mixed farming, C3 grassland    !   2 	 !   16: Water Bodies
    ! 16  Irrigated crop                      !   2 	 !   17: Herbaceous wetland
    ! 17  Bog or marsh                        !   9 	 !   18: Wooded wetland
    ! 18  Wooded grassland                    !   6 	 !   19: Barren or sparsely vegetated
    ! 19  Urban and built up                  !   1 	 !   20: Herbaceous Tundra
    ! 20  Wetland evergreen broadleaf tree    !   6 	 !   21: Wooded Tundra
    !		                                         !   22: Mixed Tundra
    !		                                         !   23: Bare Ground Tundra
    !	 	                                         !   24: Snow or Ice
    !--------------------------------------------------------------------------------------------
!KML: Os valores de ixxxlu foram trocados para fazer a equivalencia entre as classes do
! do LEAF3 e do USGS porque as tabelas do Wesely estao em funcao das classes do USGS (e nao das classes do
!Wesely).

!kml    ixxxlu(1)  =7
!kml    ixxxlu(2)  =7
!kml    ixxxlu(3)  =8
!kml    ixxxlu(4)  =5
!kml    ixxxlu(5)  =5
!kml    ixxxlu(6)  =4
!kml    ixxxlu(7)  =4
!kml    ixxxlu(8)  =2
!kml    ixxxlu(9)  =3
!kml    ixxxlu(10) =8
!kml    ixxxlu(11) =8
!kml    ixxxlu(12) =3
!kml    ixxxlu(13) =3
!kml    ixxxlu(14) =6
!kml    ixxxlu(15) =2
!kml    ixxxlu(16) =2
!kml    ixxxlu(17) =9
!kml    ixxxlu(18) =6
!kml    ixxxlu(19) =1
!kml    ixxxlu(20) =6



    ixxxlu(1)  =16
    ixxxlu(2)  =24
    ixxxlu(3)  =19
    ixxxlu(4)  =14
    ixxxlu(5)  =12
    ixxxlu(6)  =11
    ixxxlu(7)  =13
    ixxxlu(8)  =7
    ixxxlu(9)  =8
    ixxxlu(10) =19
    ixxxlu(11) =22
    ixxxlu(12) =8
    ixxxlu(13) =8
    ixxxlu(14) =15
    ixxxlu(15) =5
    ixxxlu(16) =3
    ixxxlu(17) =17
    ixxxlu(18) =10
    ixxxlu(19) =1
    ixxxlu(20) =18


  END SUBROUTINE dep_init

  ! **********************************************************************
  SUBROUTINE deppart(rmol,ustar,rh,clw,iland,dvpart,dvfog)
    !     THIS SUBROUTINE CALCULATES SURFACE DEPOSITION VELOCITIES
    !     FOR FINE AEROSOL PARTICLES ACCORDING TO THE MODEL OF
    !     J. W. ERISMAN, A. VAN PUL, AND P. WYERS,
    !     ATMOSPHERIC ENVIRONMENT 28 (1994), 2595-2607
    !     WRITTEN BY WINFRIED SEIDL, APRIL 1997
    !     MODIFIED BY WINFRIED SEIDL, MARCH 2000
    !            FOR MM5 VERSION 3
    ! ----------------------------------------------------------------------

    REAL    , INTENT(IN)  :: rmol
    REAL    , INTENT(IN)  :: ustar
    REAL    , INTENT(IN)  :: rh
    REAL    , INTENT(IN)  :: clw
    INTEGER , INTENT(IN)  :: iland
    REAL    , INTENT(OUT) :: dvpart
    REAL    , INTENT(OUT) :: dvfog

    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC exp
    ! ..
    dvpart = ustar/kpart(iland)
    !print*,ustar,iland,kpart(iland)

    IF (rmol<0.) THEN
       !         INSTABLE LAYERING CORRECTION
       dvpart = dvpart*(1.+(-300.*rmol)**0.66667)
    END IF
    IF (rh>80.) THEN
       !         HIGH RELATIVE HUMIDITY CORRECTION
       !         ACCORDING TO J. W. ERISMAN ET AL.
       !         ATMOSPHERIC ENVIRONMENT 31 (1997), 321-332
       dvpart = dvpart*(1.+0.37*EXP((rh-80.)/20.))
    END IF
    !print*,'dvpart=', dvpart
    !srf - not using fog deposition (in case on, set clw,iland=iveg)
    RETURN

    !       SEDIMENTATION VELOCITY OF FOG WATER ACCORDING TO
    !       R. FORKEL, W. SEIDL, R. DLUGI AND E. DEIGELE
    !       J. GEOPHYS. RES. 95D (1990), 18501-18515
    dvfog = 0.06*clw
!kml    IF (ixxxlu(iland)==5) THEN
    IF (ixxxlu(iland)==12 .OR.ixxxlu(iland)==14 ) THEN

       !         TURBULENT DEPOSITION OF FOG WATER IN CONIFEROUS FOREST ACCORDI
       !         A. T. VERMEULEN ET AL.
       !         ATMOSPHERIC ENVIRONMENT 31 (1997), 375-386
       dvfog = dvfog + 0.195*ustar*ustar
    END IF

  END SUBROUTINE deppart
  ! **********************************************************************
  SUBROUTINE depvel(numchem,rmol,zr,z0,ustar,vgpart)!,polint)
    !     THIS FUNCTION HAS BEEN DESIGNED TO EVALUATE AN UPPER LIMIT
    !     FOR THE POLLUTANT DEPOSITION VELOCITY AS A FUNCTION OF THE
    !     SURFACE ROUGHNESS AND METEOROLOGICAL CONDITIONS.
    !     PROGRAM WRITTEN BY GREGORY J.MCRAE (NOVEMBER 1977)
    !         Modified by Darrell A. Winner  (Feb. 1991)
    !                  by Winfried Seidl     (Aug. 1997)
    !.....PROGRAM VARIABLES...
    !     RMOL     - RECIPROCAL OF THE MONIN-OBUKHOV LENGTH
    !     ZR       - REFERENCE HEIGHT
    !     Z0       - SURFACE ROUGHNESS HEIGHT
    !     SCPR23   - (Schmidt #/Prandtl #)**(2/3) Diffusion correction fact
    !     UBAR     - ABSOLUTE VALUE OF SURFACE WIND SPEED
    !     DEPVEL   - POLLUTANT DEPOSITION VELOCITY
    !     Vk       - VON KARMAN CONSTANT
    !     USTAR    - FRICTION VELOCITY U*
    !     POLINT   - POLLUTANT INTEGRAL
    !.....REFERENCES...
    !     MCRAE, G.J. ET AL. (1983) MATHEMATICAL MODELING OF PHOTOCHEMICAL
    !       AIR POLLUTION, ENVIRONMENTAL QUALITY LABORATORY REPORT 18,
    !       CALIFORNIA INSTITUTE OF TECHNOLOGY, PASADENA, CALIFORNIA.
    !.....RESTRICTIONS...
    !     1. THE MODEL EDDY DIFFUSIVITIES ARE BASED ON MONIN-OBUKHOV
    !        SIMILARITY THEORY AND SO ARE ONLY APPLICABLE IN THE
    !        SURFACE LAYER, A HEIGHT OF O(30M).
    !     2. ALL INPUT UNITS MUST BE CONSISTENT
    !     3. THE PHI FUNCTIONS USED TO CALCULATE THE FRICTION
    !        VELOCITY U* AND THE POLLUTANT INTEGRALS ARE BASED
    !        ON THE WORK OF BUSINGER ET AL.(1971).
    !     4. THE MOMENTUM AND POLLUTANT DIFFUSIVITIES ARE NOT
    !        THE SAME FOR THE CASES L<0 AND L>0.

    INTEGER , INTENT(IN)  :: numchem ! (DMK) not used
    REAL    , INTENT(IN)  :: rmol
    REAL    , INTENT(IN)  :: zr
    REAL    , INTENT(IN)  :: z0
    REAL    , INTENT(IN)  :: ustar
    REAL    , INTENT(OUT) :: vgpart

    ! ..
    ! .. Array Arguments ..
    !        REAL :: depv(numchem)
    ! ..
    ! .. Local Scalars ..
    REAL :: ao, ar, polint
    INTEGER :: l
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC alog
    ! ..
    !     Set the von Karman constant
    REAL, PARAMETER :: vk = 0.4

    !     DETERMINE THE STABILITY BASED ON THE CONDITIONS
    !             1/L < 0 UNSTABLE
    !             1/L = 0 NEUTRAL
    !             1/L > 0 STABLE

    IF (rmol<0) THEN
       ar = ((1.0-9.0*zr*rmol)**(0.25)+0.001)**2
       ao = ((1.0-9.0*z0*rmol)**(0.25)+0.001)**2
       polint = 0.74*(alog((ar-1.0)/(ar+1.0))-alog((ao-1.0)/(ao+1.0)))
    ELSE IF (rmol==0) THEN
       polint = 0.74*alog(zr/z0)
    ELSE
       polint = 0.74*alog(zr/z0) + 4.7*rmol*(zr-z0)
    END IF

    !     CALCULATE THE Maximum DEPOSITION VELOCITY

    !DO l = 1, numchem
    !  depv(l) = ustar*vk/(2.0*scpr23(l)+polint)
    !END DO
    vgpart = ustar*vk/polint
    !print*,polint,ustar,vgpart,z0,ustar*vk/polint
    RETURN
  END SUBROUTINE depvel

  !========================================================================

  ! -----------------   AEROSOL DRY DEP AND SEDIMENTATION -----------------

  !========================================================================

  SUBROUTINE dry_dep_sedim_particles(m1,m2,m3,naddsc,npatch              &
                                    ,ia,iz,ja,jz,r_aer,temp3d,air_dens3d &
                                    ,temps,dens,vels,rvs,Zi,ustar,tstar  &
                                    ,patch_area,veg,Z0m,nmodes,nspecies  &
                                    ,part_radius,part_dens,mode_alloc,on &
                                    ,g,dd_sedim)

    INTEGER , INTENT(IN) :: m1
    INTEGER , INTENT(IN) :: m2
    INTEGER , INTENT(IN) :: m3
    INTEGER , INTENT(IN) :: naddsc ! (DMK) not used
    INTEGER , INTENT(IN) :: npatch
    INTEGER , INTENT(IN) :: ia
    INTEGER , INTENT(IN) :: iz
    INTEGER , INTENT(IN) :: ja
    INTEGER , INTENT(IN) :: jz

    REAL    , INTENT(IN) :: r_aer(m2,m3,npatch)
    REAL    , INTENT(IN) :: temp3d(m1,m2,m3)
    REAL    , INTENT(IN) :: air_dens3d(m1,m2,m3)
    REAL    , INTENT(IN) :: temps(m2,m3)
    REAL    , INTENT(IN) :: dens(m2,m3)
    REAL    , INTENT(IN) :: vels(m2,m3)
    REAL    , INTENT(IN) :: rvs(m2,m3) ! (DMK) not used
    REAL    , INTENT(IN) :: Zi(m2,m3)
    REAL    , INTENT(IN) :: ustar(m2,m3,npatch)
    REAL    , INTENT(IN) :: tstar(m2,m3,npatch)
    REAL    , INTENT(IN) :: patch_area(m2,m3,npatch)
    REAL    , INTENT(IN) :: veg(m2,m3,npatch)
    REAL    , INTENT(IN) :: Z0m(m2,m3,npatch)

    ! aer1_list
    INTEGER , INTENT(IN) :: nmodes
    INTEGER , INTENT(IN) :: nspecies
    REAL    , INTENT(IN) :: part_radius(nmodes,nspecies)
    REAL    , INTENT(IN) :: part_dens(nmodes,nspecies)
    INTEGER , INTENT(IN) :: mode_alloc(nmodes,nspecies_aer)
    INTEGER , INTENT(IN) :: on

    ! rconstants
    REAL    , INTENT(IN) :: g

    TYPE(sedim_type), INTENT(INOUT) :: dd_sedim(naer_transported)

    REAL    :: vdtmp
    INTEGER :: i,j,ipatch,ispc

!    IF(.NOT. aer_alloc) THEN
!       CALL alloc_aer_sedim(m1,m2,m3,npatch,ngrids, &
!                            nmodes,nspecies,mode_alloc, &
!                            on,mmzp,mmxp,mmyp)
!    END IF

    IF( NAER_TRANSPORTED == 0 ) RETURN

    !- sedimentation  parameterization
    CALL sedim_particles_3d(m1,m2,m3,npatch,ia,iz,ja,jz,temp3d,air_dens3d, &
                            nmodes,nspecies,part_radius,part_dens,g,dd_sedim)


    !- dry deposition parameterization

    !- laminar sub-layer resistance
    CALL lsl_particles(m2,m3,npatch,ia,iz,ja,jz &
         ,temps,dens,vels,rvs,Zi,ustar,tstar,patch_area,veg,Z0m &
         ,nmodes,nspecies,part_radius,part_dens &
         ,dd_sedim)

    !- particles deposition velocity (m/s)

    DO j = ja,jz
       DO i = ia,iz
          DO ispc=naer_a,naer_z

             dd_sedim(ispc)%v_dep_part(i,j) = 0.

             DO ipatch = 1,npatch

                IF (patch_area(i,j,ipatch) .GE. .009) THEN


                   vdtmp = dd_sedim(ispc)%v_sed_part(1,i,j) + &
                        1./( r_aer(i,j,ipatch) + dd_sedim(ispc)%r_lsl_part(i,j,ipatch) + &
                        r_aer(i,j,ipatch) * dd_sedim(ispc)%r_lsl_part(i,j,ipatch) * &
                        dd_sedim(ispc)%v_sed_part(1,i,j)				  )

                   dd_sedim(ispc)%v_dep_part(i,j) = dd_sedim(ispc)%v_dep_part(i,j)  + &
                        patch_area(i,j,ipatch)*vdtmp

                ENDIF
             ENDDO

             !print*,'V_dep (cm/s)',i,j,dd_sedim(ispc,ngrid)%v_dep_part(i,j)*100.
             !print*,'r_lsl-v_sed =', dd_sedim(ispc,ngrid)%r_lsl_part(i,j,1:npatch), &
             !	dd_sedim(ispc,ngrid)%v_sed_part(1,i,j) !,r_aer(i,j,ipatch)

          ENDDO
       ENDDO
    ENDDO

    RETURN

    !- replace the sedimentation velocity at surface with the deposition velocity
    !- (which includes v_sed_part plus the turbulent flux in laminar sub-layer)

    DO j = ja,jz
       DO i = ia,iz
          DO ispc=naer_a,naer_z
             dd_sedim(ispc)%v_sed_part(1,i,j) = dd_sedim(ispc)%v_dep_part(i,j)
          ENDDO
       ENDDO
    ENDDO


  END SUBROUTINE dry_dep_sedim_particles
  !========================================================================

  SUBROUTINE lsl_particles(m2,m3,npatch,ia,iz,ja,jz  &
                          ,temps,dens,vels,rvs,Zi,ustar &
                          ,tstar,patch_area,veg,Z0m     &
                          ,nmodes,nspecies,part_radius  &
                          ,part_dens                    &
                          ,dd_sedim)

    INTEGER , INTENT(IN) :: m2
    INTEGER , INTENT(IN) :: m3
    INTEGER , INTENT(IN) :: npatch
    INTEGER , INTENT(IN) :: ia
    INTEGER , INTENT(IN) :: iz
    INTEGER , INTENT(IN) :: ja
    INTEGER , INTENT(IN) :: jz

    REAL    , INTENT(IN) :: temps(m2,m3)
    REAL    , INTENT(IN) :: dens(m2,m3)
    REAL    , INTENT(IN) :: vels(m2,m3)
    REAL    , INTENT(IN) :: rvs(m2,m3) ! (DMK) not used
    REAL    , INTENT(IN) :: Zi(m2,m3)
    REAL    , INTENT(IN) :: ustar(m2,m3,npatch)
    REAL    , INTENT(IN) :: tstar(m2,m3,npatch)
    REAL    , INTENT(IN) :: patch_area(m2,m3,npatch)
    REAL    , INTENT(IN) :: veg(m2,m3,npatch)
    REAL    , INTENT(IN) :: Z0m(m2,m3,npatch) ! (DMK) not used

    ! aer1_list
    INTEGER , INTENT(IN) :: nmodes
    INTEGER , INTENT(IN) :: nspecies
    REAL    , INTENT(IN) :: part_radius(nmodes,nspecies)
    REAL    , INTENT(IN) :: part_dens(nmodes,nspecies)

    TYPE(sedim_type), INTENT(INOUT) :: dd_sedim(naer_transported)

    REAL,PARAMETER :: kB = 1.3807e-23      ! const Boltzmann - kg m^2 s^-2 K^-1 molecule^-1
    REAL,PARAMETER :: ASP = 1.257          ! 1.249
    REAL,PARAMETER :: BSP = 0.4            ! 0.42
    REAL,PARAMETER :: CSP = 1.1            ! 0.87

    REAL,PARAMETER :: M_AVEG = 4.8096e-26  ! average mass of one molecure - kg  molecule^-1
    REAL,PARAMETER :: vonK = 0.40          ! von Karman constant
    REAL,PARAMETER :: Cpd = 1004.67        ! specific heat of dry air [J/kg/K]
    REAL,PARAMETER :: em23 = -2./3., ep23 = +2./3., ep13=1./3.    ! exponents 2./3. 1/3
    REAL,PARAMETER :: pi = 3.1415927, g = 9.80

    INTEGER :: i,j,ipatch,ispc
    REAL    :: wptp,wstar
    !
    ! dd_sedim(ispc,ng)%v_sed      == particle sedimentation velocity   (m/s)
    ! dd_sedim(ispc,ng)%r_lsl_part == resistance to molecular diffusion (s/m)
    !
    REAL Kn,n_air,Gi,v_air  ,mfp,D,nu,Sc,St,Kd,Dh,Z0h,Pr

    REAL,PARAMETER :: limite = -30.

    DO j = ja,jz
       DO i = ia,iz

          !- several particle/environment properties

          !- mean speed of air molecules (m/s)
          !  v_air = sqrt( 8. * kB   * temps(i,j) / (pi * M_AVEG) )
          v_air = SQRT( 7.3102e+2 * temps(i,j)			  )

          !-dynamic viscosity of air (kg m^-1 s^-1)
          !  n_air = 1.8325e-5*(416.16/(temps(i,j)+120.))*(temps(i,j)/296.16)**1.5
          !optimized version
          n_air = 1.8325e-5*(416.16/(temps(i,j)+120.))*(temps(i,j)/296.16) * &
               SQRT(temps(i,j)/296.16)

          !- mean free path of an air molecule (m)
          mfp = 2.* n_air /(dens(i,j)*v_air)

          !-  kinematic viscosity of air ()
          nu = n_air/dens(i,j)

          DO ispc=naer_a,naer_z

             !- Knudsen number
             Kn = mfp/part_radius(ind_mode(ispc),ind_aer(ispc))


             !- Slip correction factor (Gi)
             Gi = 1. + Kn*( ASP + BSP*EXP(-CSP/Kn) )

             !- Schmidt number determination (Jacobson)
             !-- molecular diffusivity (Brownian diffusivity coeficient)
             !    D =  (kB/M_AVEG) * temps(i,j) * Gi / (6.*pi*part_radius*n_air)

             D =  kB * temps(i,j) * Gi / (6.*pi*part_radius(ind_mode(ispc),ind_aer(ispc))*n_air)


             !-  Schmidt number
             Sc = nu/D

             !- laminar sub-layer resistance for particles (s m^-1)
             DO ipatch=1,npatch

                !-  Stokes number determination (Slinn, 1980)
                !-  St =     ustar^2 * V_sedim /( g * kinematic viscosity of air)
                St = MAX(ustar(i,j,ipatch)**2. * dd_sedim(ispc)%v_sed_part(1,i,j) &
                     / (g * nu), 0.01)

                !-  laminar sub-layer resistance for particles (s m^-1)

                IF(ipatch == 1) THEN !- water

                   !-from Slinn 1980 (Atmos. Env.)
                   !dd_sedim(ispc,ng)%r_lsl_part(i,j,ipatch) = ( vonK * vels(i,j)/ustar(i,j,ipatch)**2. ) / &
                   !	                  ( 1./sqrt(Sc) + 10.**(-3./St) )
                   dd_sedim(ispc)%r_lsl_part(i,j,ipatch) = ( vonK * vels(i,j)/ustar(i,j,ipatch)**2. ) / &
                        ( 1./SQRT(Sc) + 10.**MAX((-3./St),limite) )


                ELSE
                   IF(patch_area(i,j,ipatch) > 0.009) THEN

                      !-- for smooth land surfaces !- bare ground (3), desert(3) or ice (2)
                      !-- Seinfeld & Pandis (1998)
                      IF(NINT(veg(i,j,ipatch)) == 3 ) THEN

                         ! dd_sedim(ispc,ng)%r_lsl_part(i,j,ipatch) = 1./(ustar(i,j,ipatch) * (Sc**em23 + 10.**(-3./St)))
                         dd_sedim(ispc)%r_lsl_part(i,j,ipatch) = 1./(ustar(i,j,ipatch) * &
                              (Sc**em23 + 10.**MAX((-3./St),limite)))


                      ELSE
                         !- lsl resistance according  Jacobson(1999)
                         !
                         !- thermal conductivity of dry air (Kd)
                         !Kd = 0.023807 + 7.1128e-5*(temps(i,j) - 273.15) !- Eq.(2.3)
                         !- Prandt number
                         !Pr =  n_air*Cpd*(1.+0.859*rvs(i,j))/Kd           !- Eq.(17.32)

                         !- energy moisture roughness lengths (Z0h)                 !- Eq.(8.10)
                         !-- molecular thermal diffusion coeff. (m^2 s^-1)
                         !Dh=Kd/(dens(i,j)*Cpd)
                         !- Z0h
                         !Z0h=Dh/(vonK*ustar(i,j,ipatch))
                         !- lsl resistance according Jacobson Eq. (20.14)
                         ! dd_sedim(ispc,ng)%r_lsl_part(i,j,ipatch) =  log( Z0m(i,j,ipatch)/Z0h )* ( (Sc/Pr)**ep23 ) &
                         !	                     /  ( vonK*ustar(i,j,ipatch) )

                         !---------
                         !- lsl resistance according :
                         !- from Wesely et al. (1985) following Slinn (1982)
                         !- also Binkowski & Shankar, JGR, 100,D12,26191-26209, 1995
                         !- Rb= (u* (1+ (w*/u*)^2)) (Sc^2./3 + 10^(-3/St))

                         wptp  = - ustar(i,j,ipatch) * tstar(i,j,ipatch) ! sensible heat flux
                         wstar = ( MAX (0., g* Zi(i,j)*  wptp/temps(i,j) ) )**ep13
                         !
                         !dd_sedim(ispc,ng)%r_lsl_part(i,j,ipatch) = 1./(  					  &
                         !				     ustar(i,j,ipatch)* 			 &
                         !				     (1. + 0.24*(wstar/ustar(i,j,ipatch))**2.)*  &
                         !				     ( (Sc**em23 + 10.**(-3./St)) )		 &
                         !				    )
                         dd_sedim(ispc)%r_lsl_part(i,j,ipatch) = 1./(						&
                              ustar(i,j,ipatch)*  			&
                              (1. + 0.24*(wstar/ustar(i,j,ipatch))**2.)*  &
                              ( (Sc**em23 + 10.**MAX((-3./St),limite) ) ) )

                         !       print*,'================'
                         !       print*,'LSL_WE',i,j,ipatch,nint(veg(i,j,ipatch)),patch_area(i,j,ipatch)
                         !	PRINT*,wstar,ustar(i,j,ipatch),r_lsl(i,j,ipatch)

                      ENDIF
                   ENDIF
                ENDIF


                !print*,'laminar sub-layer resistance'
                !print*,'LSL_J',i,j,ipatch,dd_sedim(ispc,ng)%r_lsl_part(i,j,ipatch)
                !print*,'ZOM',Z0m(i,j,ipatch)
                !print*,'Z0H',Z0h
                !print*,'SC PR',Sc,Pr
                !print*,'u*',ustar(i,j,ipatch)

             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE lsl_particles


  !------------------------------------------------------------------------
  SUBROUTINE sedim_particles_3d(m1,m2,m3,npatch,ia,iz,ja,jz,temp3d       &
                               ,air_dens3d,nmodes,nspecies,part_radius   &
                               ,part_dens,g                              &
                               ,dd_sedim)

    INTEGER , INTENT(IN) :: m1
    INTEGER , INTENT(IN) :: m2
    INTEGER , INTENT(IN) :: m3
    INTEGER , INTENT(IN) :: npatch
    INTEGER , INTENT(IN) :: ia
    INTEGER , INTENT(IN) :: iz
    INTEGER , INTENT(IN) :: ja
    INTEGER , INTENT(IN) :: jz
    REAL    , INTENT(IN) :: temp3d(m1,m2,m3)
    REAL    , INTENT(IN) :: air_dens3d(m1,m2,m3)

    ! aer1_list
    INTEGER , INTENT(IN) :: nmodes
    INTEGER , INTENT(IN) :: nspecies
    REAL    , INTENT(IN) :: part_radius(nmodes,nspecies)
    REAL    , INTENT(IN) :: part_dens(nmodes,nspecies)

    ! rconstants
    REAL    , INTENT(IN) :: g

    TYPE(sedim_type), INTENT(INOUT) :: dd_sedim(naer_transported)

    REAL,PARAMETER :: ASP = 1.257          ! 1.249
    REAL,PARAMETER :: BSP = 0.4            ! 0.42
    REAL,PARAMETER :: CSP = 1.1            ! 0.87
    REAL,PARAMETER :: onesix = 1./6.

    INTEGER :: i,j,ipatch,k,ispc

    !
    !For small particles diameter (<20 micrometer) (low Reynolds number)
    !the fall velocity is given by (Seinfeld & Pandis, Jacobson)
    ! V_s = 2 r**2 (rho_p - rho_air) * g * Gi / 9 n_air
    ! where:
    ! r == radius of particle (m)
    ! rho_p,a = density of particle, air (kg/m^3)
    ! g = gravity accel (9.8 m/s^2)
    ! Gi = Kn(A' + B' + C' exp(-C'/Kn) , Kn = Knudsen number
    ! n_air = dynamic viscosity of air
    !- constantes
    !parameter (ASP=1.257,BSP=0.4,CSP=1.1) ! Seinfeld & Pandis
    !parameter (kB = 1.3807e-23)     ! const Boltzmann - kg m^2 s^-2 K^-1 molecule^-1
    !parameter (M_AVEG = 4.8096e-26) ! average mass of one molecure - kg  molecule^-1
    !
    !real, dimension(m1,m2,m3,nmodes,aer_nspecies)    :: v_sed

    REAL A,B,C,Kn,n_air,Gi,v_air,mfp,nu,re,z,b0,x,y
    !- for FALL 3d formulation
    REAL Cd
    REAL, PARAMETER :: eps = 1.e-6
    REAL, PARAMETER :: psi = 0.95
    REAL, PARAMETER :: k1 = 3./(1. + 2.*psi**(-0.5))    ,&
                       k2 = 10**(1.84148*(-log(psi))**0.5743)

    real :: d1,d2,d3,d4,d5


    DO j = ja,jz
       DO i = ia,iz
          DO k=1,m1

             !- several particle/environment properties

             !- mean speed of air molecules (m/s)
             !  v_air = sqrt( 8. * kB	* temp3d(k,i,j) / (pi * M_AVEG) )
             v_air = SQRT( 7.3102e+2 * temp3d(k,i,j)  		      )

             !-dynamic viscosity of air (kg m^-1 s^-1)
             !  n_air = 1.8325e-5*(416.16/(temp3d(k,i,j)+120.))*(temp3d(k,i,j)/296.16)**1.5
             !optimized version
             n_air = 1.8325e-5*(416.16/(temp3d(k,i,j)+120.))*(temp3d(k,i,j)/296.16) * &
                  SQRT(temp3d(k,i,j)/296.16)
             !- kinematic viscosity of air ()
             nu = n_air/air_dens3d(k,i,j)
             !- mean free path of an air molecule (m)
             mfp = 2.* n_air /(air_dens3d(k,i,j)*v_air)

             DO ispc=naer_a,naer_z

                !print*,'particle=',ispc


                !- Knudsen number
                Kn = mfp/part_radius(ind_mode(ispc),ind_aer(ispc))

                !- Slip correction factor (Gi)
                Gi = 1. + Kn*( ASP + BSP*EXP(-CSP/Kn) )


                !- This is regime 1 (Pruppacher and Klett (chap. 10)) , first guess for terminal velocity
                !- 0.5 < part_radius < 10 micrometers  / 1.e-6 < Re < 1e-2
                !
                !- particle sedimentation velocity (m/s)
                !  v_sed(k,i,j) = (2./9.)*g*part_radius**2  *(part_dens-air_dens3d(k,i,j))*Gi/n_air
                !- opt
                dd_sedim(ispc)%v_sed_part(k,i,j) = 2.18*part_radius(ind_mode(ispc),ind_aer(ispc))**2.  &
                     *(part_dens(ind_mode(ispc),ind_aer(ispc)))*Gi/n_air !part_dens >>dens_air


                ! Reynolds number
                !   re = 2*part_radius*v_sed_part(i,j)/nu
                re = 2. * air_dens3d(k,i,j) * part_radius(ind_mode(ispc),ind_aer(ispc)) * &
                     dd_sedim(ispc)%v_sed_part(k,i,j)/n_air

                go to 343
                !-----
		!- srf addded new formulation
		!- Approach of FALL 3d (Costa et al., 2005)
		!
		!- avoid division by zero in case  re = zero.
		re = re + eps
		!- drag coefficient
		Cd = 24./(re * k1 ) * (1. + 0.1118 *( re *( k1*k2 )**0.6567 )) &
		     + 0.4305 * k2 / ( 1. + 3305. / (re * k1 * k2 ))

                !-- particle sedimentation velocity (m/s)
		dd_sedim(ispc)%v_sed_part(k,i,j) = (4.* g *part_dens(ind_mode(ispc),ind_aer(ispc)) * &
		                                    2.  *part_radius(ind_mode(ispc),ind_aer(ispc))) /&
		                                   (3.* Cd * air_dens3d(k,i,j))

                dd_sedim(ispc)%v_sed_part(k,i,j) = ( dd_sedim(ispc)%v_sed_part(k,i,j) ) ** 0.5

		cycle ! not using the corrections below

		!------
 343            continue
                IF( re .GE. 0.01 .AND. re .LE. 300.) THEN
                   !  This is "regime 2" in Pruppacher and Klett (chap. 10).
                   x = LOG(24.*re/Gi)
                   y = -3.18657 + x*(0.992696    - x*(.00153193 &
                        + x*(0.000987059 + x*(.000578878&
                        - x*(8.55176E-05 - x* 3.27815E-06 )))))
                   !if( y .lt. -675. ) y = -675.
                   !if( y .ge.  741. ) y =  741.
                   y=MIN(MAX(y,-675.),741.)

                   re = EXP(y)*Gi

                   dd_sedim(ispc)%v_sed_part(k,i,j) = re * n_air / &
                        (2.*part_radius(ind_mode(ispc),ind_aer(ispc))* air_dens3d(k,i,j))
                ENDIF

                IF( re > 300. )THEN
                   !  This is "regime 3" in Pruppacher and Klett (chap. 10).
		   ! - units mus be CGS
		   d1= air_dens3d(k,i,j) *1.e-3 ! g/cm�
                   d2= part_dens(ind_mode(ispc),ind_aer(ispc)) *1.e-3 ! g/cm�
                   d3= g *100. ! cm/s2
                   d4=n_air*10. ! g/cm/s
                   ! - old
		   !z  = ((1.e6*air_dens3d(k,i,j)**2) &
                   !     /  (g * part_dens(ind_mode(ispc),ind_aer(ispc)) * n_air**4) )**(onesix)

                   z  = ((1.e6*d1**2)  /  (d3 * d2* d4**4) )**(onesix)

                   d5=dd_sedim(ispc)%v_sed_part(k,i,j)*100. !cm/s

                   !b0 = (24.*dd_sedim(ispc)%v_sed_part(k,i,j) *n_air)/100.
                   b0  = (24.*d5*d4)/100.

		   x  = LOG(z*b0)
                   y  = -5.00015 + x*(5.23778   - x*(2.04914 - x*(0.475294 &
                        - x*(0.0542819 - x* 0.00238449 ))))
                   !if( y .lt. -675. )  y = -675.0
                   !if( y .ge.  741. )  y =  741.0
                   y=MIN(MAX(y,-675.),741.)
                   re = z*EXP(y)*Gi

                   dd_sedim(ispc)%v_sed_part(k,i,j) = re * n_air / &
                        (2.*part_radius(ind_mode(ispc),ind_aer(ispc)) * air_dens3d(k,i,j))

                ENDIF
                !tmp
                !for testing
                !dd_sedim(ispc,ngrid)%v_sed_part(k,i,j) = float(ispc)/500.
                !if(k==1 .or. K==2)print*,'v-sed=',k,dd_sedim(ispc,ngrid)%v_sed_part(k,i,j)
                !tmp


             ENDDO ! ispc loop
          ENDDO  !   k loop
       ENDDO
    ENDDO


  !print*,'max val min val '
  !do  ispc=naer_a,naer_z
  !
  ! if(ispc>7)   print*,1000.*part_radius(ind_mode(ispc),ind_aer(ispc)),&
  !	     maxval(dd_sedim(ispc)%v_sed_part),minval(dd_sedim(ispc)%v_sed_part)
  !
  !enddo
  !stop 333

      RETURN
  END SUBROUTINE sedim_particles_3d


  !------------------------------------------------------------------------
  SUBROUTINE fa_preptc_with_sedim(m1,m2,m3,vt3da,vt3db,vt3dc,vt3df    &
                                 ,vt3dk,dn0,rtgt,f13t,f23t,dtlt       &
                                 ,N_current_scalar,num_scalar_aer_1st &
                                 ,wp,wc,vt3dp,nzpmax,hw4,dzm,dzt      &
                                 ,dd_sedim)

    INTEGER , INTENT(IN)    :: m1
    INTEGER , INTENT(IN)    :: m2
    INTEGER , INTENT(IN)    :: m3

    REAL    , INTENT(IN)    :: vt3da(m1,m2,m3)
    REAL    , INTENT(IN)    :: vt3db(m1,m2,m3)
    REAL    , INTENT(OUT)   :: vt3dc(m1,m2,m3)

    REAL    , INTENT(INOUT) :: vt3df(m1,m2,m3)

    REAL    , INTENT(INOUT) :: vt3dk(m1,m2,m3)
    REAL    , INTENT(IN)    :: dn0(m1,m2,m3)

    REAL    , INTENT(IN)    :: rtgt(m2,m3)

    REAL    , INTENT(IN)    :: f13t(m2,m3)
    REAL    , INTENT(IN)    :: f23t(m2,m3)

    !srf- aerosol section
    REAL    , INTENT(IN)    :: dtlt
    INTEGER , INTENT(IN)    :: N_current_scalar
    INTEGER , INTENT(IN)    :: num_scalar_aer_1st

    REAL    , INTENT(IN)    :: wp(m1,m2,m3)
    REAL    , INTENT(IN)    :: wc(m1,m2,m3)
    REAL    , INTENT(INOUT) :: vt3dp(m1,m2,m3)

    ! grid_dims
    INTEGER , INTENT(IN)    :: nzpmax

    ! mem_grid
    REAL    , INTENT(IN)    :: hw4(nzpmax)
    REAL    , INTENT(IN)    :: dzm(nzpmax)
    REAL    , INTENT(IN)    :: dzt(nzpmax)

    TYPE(sedim_type), INTENT(INOUT) :: dd_sedim(naer_transported)

    INTEGER :: j,i,k,im,ip,jm,jp
    REAL    :: c1,c2,c3,c4,rtgti

    ! VT3DA, VT3DB, and VT3DC are input as the velocity components (averaged
    ! between past and current time levels) times dtlt.

    ! Add contribution to VT3DC from horiz winds crossing sloping sigma surfaces,
    !    and include 1/rtgt factor in VT3DC
    ! Compute half Courant numbers: VT3DD, VT3DE, and VT3DF
    ! Compute weight at scalar point: VT3DH
    ! Compute advective weights for the linear term: VT3DI, VCTR1, and VCTR2

    INTEGER ISPC

    !------------------ aerosol mapping
    ! num_scalar_aer_1st                  corresponds do aerosol  1 (=naer_a)
    ! num_scalar_aer_1st +1               corresponds do aerosol  2
    ! ...
    ! num_scalar_aer_1st +NAER_TRANSPORTED corresponds do aerosol naer_z
    ! in this case, to access the correct V_SED use "ISPC"
    ISPC = N_current_scalar - num_scalar_aer_1st + 1

    ! sedim is included at vertical direction (vt3dc)
    DO j = 1,m3
       DO i = 1,m2
          DO k = 1,m1
             !if(j==20 .and. i==20) print*,'particle',ISPC,sedim(ispc,ngrid)%v_sed_part(k,i,j),wp(k,i,j)
             vt3dc(k,i,j) = ( ( wp(k,i,j) + wc(k,i,j) )* 0.5 - dd_sedim(ispc)%v_sed_part(k,i,j) )* dtlt
             !-test	 vt3dc(k,i,j) = ( - sedim(ispc,ngrid)%v_sed_part(k,i,j) )* dtlt
          ENDDO
       ENDDO
    ENDDO
    !print*,'W=',wp(1:2,int(m2/2),int(m3/2)),wc(1:2,int(m2/2),int(m3/2))
!!!! (be sure the lines below are executed each timestep)!!!!
    !- only necessary one time
    IF(ISPC==1) THEN
       DO j = 1,m3
          jm = MAX(1,j-1)
          !jp = min(m3,j+1)
          DO i = 1,m2
             im = MAX(1,i-1)
             !ip = min(m2,i+1)
             !rtgti = 1. / rtgt(i,j)
             DO k = 1,m1-1
                vt3dp(k,i,j) = ((vt3da(k,i,j) + vt3da(k+1,i,j)  &
                     + vt3da(k,im,j) + vt3da(k+1,im,j)) * f13t(i,j)  &
                     + (vt3db(k,i,j) + vt3db(k+1,i,j) + vt3db(k,i,jm)  &
                     + vt3db(k+1,i,jm)) * f23t(i,j)) * hw4(k)
                vt3dk(k,i,j) = dzt(k) / dn0(k,i,j)
             ENDDO

          ENDDO
       ENDDO
    ENDIF


    !-  then the fluxes are re-evaluated
    DO j = 1,m3
       jm = MAX(1,j-1)
       !jp = min(m3,j+1)
       DO i = 1,m2
          im = MAX(1,i-1)
          !ip = min(m2,i+1)
          rtgti = 1. / rtgt(i,j)
          DO k = 1,m1-1
             !-orig way
             !        vt3dc(k,i,j) = ((vt3da(k,i,j) + vt3da(k+1,i,j)  &
             !           + vt3da(k,im,j) + vt3da(k+1,im,j)) * f13t(i,j)  &
             !           + (vt3db(k,i,j) + vt3db(k+1,i,j) + vt3db(k,i,jm)  &
             !           + vt3db(k+1,i,jm)) * f23t(i,j)) * hw4(k)  &
             !           + vt3dc(k,i,j) * rtgti
             !opt way
             vt3dc(k,i,j) = vt3dp(k,i,j) &
                  + vt3dc(k,i,j) * rtgti


             vt3df(k,i,j) = .5 * vt3dc(k,i,j) * dzm(k)
             !vt3dk(k,i,j) = dzt(k) / dn0(k,i,j)

          ENDDO

       ENDDO
    ENDDO
    !do k = 1,m1-1
    !   vctr1(k) = (zt(k+1) - zm(k)) * dzm(k)
    !   vctr2(k) =  (zm(k) - zt(k)) * dzm(k)
    !enddo
    !print*,'W=',wp(1:2,int(m2/2),int(m3/2)),wc(1:2,int(m2/2),int(m3/2))

    !convert velocity components * dtlt (VT3DA, VT3DB, VT3DC)
    ! into mass fluxes times dtlt.

    DO j = 1,m3
       DO i = 1,m2
          !c1 = fmapui(i,j) * rtgu(i,j)
          !c2 = fmapvi(i,j) * rtgv(i,j)
          DO k = 1,m1-1
             vt3dc(k,i,j) = vt3dc(k,i,j) * .5  &
                  * (dn0(k,i,j) + dn0(k+1,i,j))      !air_dens3d ????
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE fa_preptc_with_sedim
  !--------------------------------------------------------------------------



END MODULE module_dry_dep
