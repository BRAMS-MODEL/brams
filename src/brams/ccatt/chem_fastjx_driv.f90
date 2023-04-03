MODULE FastJX

  USE chem1_list, ONLY: &
       maxJcomb,        & ! PARAMETER
       nfactors,        & ! PARAMETER
       factor,          & ! PARAMETER
       nr_photo,        & ! PARAMETER
       JReactionComp      ! PARAMETER

  USE rconstants , ONLY : &
       cp,                & ! PARAMETER
       cpor,              & ! PARAMETER
       p00,               & ! PARAMETER
       stefan               ! PARAMETER

  USE jx_data, ONLY: &
       jvn_,         & ! PARAMETER
       x_,           & ! PARAMETER
       w_,           & ! PARAMETER
       m_,           & ! PARAMETER
       n_,           & ! PARAMETER
       rad,          & ! PARAMETER
       szamax,       & ! PARAMETER
       zzht,         & ! PARAMETER
       a_,           & ! PARAMETER
       njval,        &
       jlabel,       &
       jfacta,       &
       atau0,        &
       atau,         &
       jtaumx,       &
       naa,          &
       nw2,          &
       nw1,          &
       wl,           &
       tqq,          &
       qrayl,        &
       fl,           &
       qo3,          &
       qo2,          &
       qaa,          &
       paa,          &
       ssa,          &
       qqq,          &
       titlej,       &
       q1d,          &
       title0

!  USE mod_dateutils, ONLY: &
!       julday                ! Function
  
  USE FastJX57, ONLY: &
       rd_xxx,        & ! Subroutine
       rd_mie,        & ! Subroutine
       rd_js,         & ! Subroutine
       photoj           ! Subroutine

  IMPLICIT NONE

  PRIVATE

  INTEGER,DIMENSION(maxJcomb,nr_photo)    :: jx2spack
  INTEGER, PUBLIC                         :: Fast_JX_initialized=0
  INTEGER,PARAMETER                       :: yes=1,no=0
  DOUBLE PRECISION  , ALLOCATABLE         :: valjl(:,:) !2-D array of J_s 

  TYPE, PUBLIC :: fast_JX_vars
     DOUBLE PRECISION, POINTER,DIMENSION  (:,:,:,:)  :: jphoto 
  END TYPE fast_JX_vars
  TYPE(fast_JX_vars), PUBLIC, ALLOCATABLE,DIMENSION(:)  :: fast_JX_g

  PUBLIC :: FastJX_driver, &
            Initialize_Fast_JX ! Subroutine

CONTAINS


  !------------------------------------------------------------------------------------------ 
  SUBROUTINE FastJX_driver(m1,m2,m3,ia,iz,ja,jz,nzpmax,ngrid,ngrids,lpw_R,rtgt,zm,imonth1, &
                           idate1,iyear1,time,dzt,mmzp,mmxp,mmyp,i0,j0,maxgrds,pp,pi0, &
                           theta,dn0,rlongup,cosz,albedt,raddatfn,do3,daer,na)
  !------------------------------------------------------------------------------------------  

    !Input
    INTEGER , INTENT(IN)  :: ia			   !i (left) model domain
    INTEGER , INTENT(IN)  :: iz			   !i (right) model domain
    INTEGER , INTENT(IN)  :: ja			   !j (top) model domain
    INTEGER , INTENT(IN)  :: jz			   !j (bottom) model domain
    INTEGER , INTENT(IN)  :: m1,m2,m3		   ! vertical levels    

    ! grid_dims
    INTEGER , INTENT(IN) :: nzpmax

    ! mem_all
    INTEGER , INTENT(IN) :: ngrid
    INTEGER , INTENT(IN) :: ngrids

    ! mem_grid
    real , INTENT(IN) :: lpw_R(m2,m3)
    REAL    , INTENT(IN) :: rtgt(m2,m3)
    REAL    , INTENT(IN) :: zm(nzpmax)
    INTEGER , INTENT(IN) :: imonth1
    INTEGER , INTENT(IN) :: idate1
    INTEGER , INTENT(IN) :: iyear1
    REAL    , INTENT(IN) :: time
    REAL    , INTENT(IN) :: dzt(nzpmax)

    ! nodemod
    INTEGER , INTENT(IN) :: mmzp(maxgrds)
    INTEGER , INTENT(IN) :: mmxp(maxgrds)
    INTEGER , INTENT(IN) :: mmyp(maxgrds)
    INTEGER , INTENT(IN) :: i0
    INTEGER , INTENT(IN) :: j0
    INTEGER , INTENT(IN) :: maxgrds

    ! mem_basic
    REAL    , INTENT(IN) :: pp(m1,m2,m3)
    REAL    , INTENT(IN) :: pi0(m1,m2,m3)
    REAL    , INTENT(IN) :: theta(m1,m2,m3)
    REAL    , POINTER    :: dn0(:,:,:)

    ! mem_radiate
    REAL    , INTENT(IN) :: rlongup(m2,m3)
    REAL    , INTENT(IN) :: cosz(m2,m3)
    REAL    , INTENT(IN) :: albedt(m2,m3)

    ! mem_globrad
    CHARACTER(LEN=256), INTENT(IN) :: raddatfn

    ! carma_fastjx
    REAL    , INTENT(IN) :: do3(m1,m2,m3)
    REAL    , INTENT(IN) :: daer(na,m1,m2,m3)
    INTEGER , INTENT(IN) :: na

    INTEGER  :: ncldx(m1,ia:iz,ja:jz)   !Index for cloud at layer L
    INTEGER  :: naer(na,m1,ia:iz,ja:jz) !Index for aerosol 1 at layer L
    !------
    INTEGER :: lpw(m2,m3)

    DOUBLE PRECISION  :: frefl        !Reflected flux?????

    !local var
    INTEGER :: ilng,jlat
    INTEGER :: l_,l1_,l2_
    INTEGER,PARAMETER :: iph=27       !File number to read
    DOUBLE PRECISION, DIMENSION(nzpmax)   :: dair,pird,prd,zhl,dztr
    DOUBLE PRECISION, DIMENSION(nzpmax+1) :: temprd
    REAL solfac
    DOUBLE PRECISION :: solf,zangle
    INTEGER :: j,k1,k2,k,l,jj,infactors,jvl_
    DOUBLE PRECISION  :: do3_(m1)
    CHARACTER(len=256)  :: dirfastjxdata 
    LOGICAL :: there

    l1_=m1
    l_=l1_-1
    l2_=2*l_+2
    jvl_=l1_

    lpw=int(lpw_R)
    !For the first time the code must read parameters files
    IF(Fast_JX_initialized==no) THEN

       !- determine dir where needed dat files are (must be same of carma tables)
       dirfastjxdata=TRIM(raddatfn(1:LEN_TRIM(raddatfn)-14))

       ! Read in fast-J X-sections (spectral data) <<<<<<<<<<<<<< new fast-JX
       INQUIRE(file=dirfastjxdata(1:LEN_TRIM(dirfastjxdata))//'FJX_spec.dat',exist=there)
       IF(.NOT.there) CALL error_fastjx(dirfastjxdata(1:LEN_TRIM(dirfastjxdata))//'FJX_spec.dat')

       CALL rd_xxx(iph,dirfastjxdata(1:LEN_TRIM(dirfastjxdata))//'FJX_spec.dat',njval,nw1,nw2, &
            title0,x_,titlej,tqq,w_,wl,fl,qrayl,qo2,qo3,q1d,qqq)


       ! Read in aerosol/cloud scattering data <<<<<<<<<<<<<<<<<< new fast-JX
       INQUIRE(file=dirfastjxdata(1:LEN_TRIM(dirfastjxdata))//'FJX_scat.dat',exist=there)
       IF(.NOT.there) CALL error_fastjx(dirfastjxdata(1:LEN_TRIM(dirfastjxdata))//'FJX_scat.dat')

       CALL rd_mie(iph,dirfastjxdata(1:LEN_TRIM(dirfastjxdata))//'FJX_scat.dat',atau, &
            atau0,jtaumx,naa,title0,a_,qaa,ssa,paa)

       ! Read in labels of photolysis rates required   >>>>> keyed to chem code
       INQUIRE(file=dirfastjxdata(1:LEN_TRIM(dirfastjxdata))//'ratj.dat',exist=there)
       IF(.NOT.there) CALL error_fastjx(dirfastjxdata(1:LEN_TRIM(dirfastjxdata))//'ratj.dat')
       CALL rd_js(iph,dirfastjxdata(1:LEN_TRIM(dirfastjxdata))//'ratj.dat',jvn_, &
            x_,njval,titlej,jlabel,jfacta)


       CALL initialize_fast_JX(mmxp,mmyp,mmzp,ngrids)

       !Crossreference between FastJX and Spack
       ALLOCATE(valjl(jvl_,njval))

       Fast_JX_initialized = yes

    END IF

    !-----------tmp
    ncldx  = 1 !int
    naer   = 1 !int
    !-----------tmp
    ! in:
    !        press(m1) = pressure profile at edges
    !        temp(m1) = = temperatures at mid-level
    ! out:

    CALL get_solfac(imonth1,idate1,iyear1,time,solfac)


    DO jlat=ja,jz
       DO ilng=ia,iz
          !-  get O3 profile from CARMA RAD model and transform to the FAST-JX axis 
          !-  orientation (height increase with k)
          !- srf: it is not being used at CARMA - fix this later
	  !DO k=1,m1
          !   do3_(k)=DBLE(do3(reverse(k,m1),ilng,jlat))
          !END DO

          k2=lpw(ilng,jlat)
          k1=k2-1

          DO k = 1,m1
             pird  (k-k1+1) = (pp(k,ilng,jlat) + pi0(k,ilng,jlat)) / cp ! Exner function
             temprd(k-k1+1) =  theta(k,ilng,jlat) * pird(k-k1+1)   ! temp (K)
             prd   (k-k1+1) = pird(k-k1+1) ** cpor * p00  * 0.01             ! units (mbars)

             dztr(k-k1+1) =  rtgt(ilng,jlat)/ dzt(k) ! meters (thickness of each layer, dzt=1/(z(k)-z(k-1)))

             !- air density in molecules.cm-2
             dair(k-k1+1) = dn0(k,ilng,jlat) * 6.022D+23 /(28.96  * 1.e-3)!(from kg.m-3 to molecules.m-3)
             !	 dair(k-k1+1) = basic_g(ngrid)%dn0(k,ilng,jlat) * 6.022D+23 /(28.96  * 1.e-3)!(from kg.m-3 to molecules.m-3)
             dair(k-k1+1) = dair(k-k1+1) * dztr(k-k1+1) *1.e-4  !(from molecules.m-3, molecules.m-2 and molecules.cm-2)

             !print*,prd(k-k1+1),dair(k-k1+1), basic_g(ngrid)% dn0(k,ilng,jlat),6.022D+23 /(28.96  * 1.e-3),dztr(k-k1+1)
             zhl(k-k1+1) =zm(k) * rtgt(ilng,jlat) * 1.e+2 !unit (cm)

             ! do3(k-k1+1,ilng,jlat) = do3(k-k1+1,ilng,jlat)*6.02D23*28.96* prd(k-k1+1)*100.D0* dztr(k-k1+1) *1.e-4 &
             !                        /(8.314D0*temprd(k-k1+1)*48.0D0)*1000.D0
             !- v1
             do3_(k-k1+1) = do3_(k-k1+1)*6.023D+23*28.96D0* prd(k-k1+1)*100.D0 *&
                  dztr(k-k1+1) *1.e-4 &
                  /(8.314D0*temprd(k-k1+1)*48.0D0)

             !- v2
             !do3_(k-k1+1) = do3_(k-k1+1) * basic_g(ngrid)% dn0(k,ilng,jlat) & ! kg[O3]/m^3
             !                            * 6.023D+23 / ( 48.0D0 * 1.e-3)    & ! molec[O3]/m^3
             !                            * dztr(k-k1+1) *1.e-4                ! molec[O3]/cm^2

             !print*,'o3=', do3_(k-k1+1),k,dair(k-k1+1); call flush(6)
          ENDDO

          temprd(1) = (rlongup (ilng,jlat) / stefan) ** 0.25

          temprd(m1+1-k1+1) = temprd(m1-k1+1)

          solf = DBLE(solfac)

          zangle = DBLE( ACOS(cosz(ilng,jlat)) )

          frefl = 0. ! set zero for now	
          !-----------------------------------------------------------------------
          CALL photoj(zangle,             &! Zenital angle
               DBLE(cosz(ilng,jlat)),     &! Zenital Cos's angle
               frefl,                     &! ???
               solf,                      &! Solar Factor
               .FALSE.,                   &! Debug ON or OFF (.true.=ON)
               l_,l1_,l2_,jvl_,           &! Size of arrays
               prd,                       &! Pressure on each level
               temprd,                    &! Temperature on each level
               dair,                      &! Air density
               do3_,                      &! Ozone column  # molecules/cm^2
               !
               DBLE(daer(1,:,ilng,jlat)), &! Aerosol Cloud optical depth 1
               DBLE(daer(2,:,ilng,jlat)), &!                             2
               DBLE(daer(3,:,ilng,jlat)), &!                             3
               naer(1,:,ilng,jlat),       &! Index for Aerosol 1
               naer(2,:,ilng,jlat),       &!                   2
               naer(3,:,ilng,jlat),       &!                   3
               DBLE(daer(4,:,ilng,jlat)), &! Cloud optical depth
               ncldx(:,ilng,jlat),        &! Index for Rayleigh phase
               zhl,                       &! Effective altitudes
               DBLE(albedt(ilng,jlat)),   &! Surface Albedo
               valjl,                     &!
               njval,                     &!
               x_,atau0,atau,w_,jtaumx,naa,m_,n_,nw2, nw1, rad, &
               szamax, zzht, wl, tqq, qrayl, fl, qo3, qo2, a_,qaa, &
               paa, ssa, qqq, titlej, q1d &
               )
          !-----------------------------------------------------------------------    
          !- converting reactions from FastJX to CATT-BRAMS/Spack
          DO l = 2,m1-1  ! loop at vertical grid
             DO j = 1,nr_photo  ! loop at photolysis reactions 

                !- zerout for accumulation
                fast_JX_g(ngrid)%jphoto(j,l,ilng,jlat) = 0.0D0

                DO infactors=1,nfactors(j)

                   jj = jx2spack(infactors,j)

                   !>>>>>>>>   ! solucao temporaria para i-np e n-pn => zerados
                   IF(jj==0) CYCLE
                   !>>>>>>>>  

                   !-srf: factor(infactors,j) must be indexed by "j" and not by "jj"	     
                   fast_JX_g(ngrid)%jphoto(j,l,ilng,jlat) = fast_JX_g(ngrid)%jphoto(j,l,ilng,jlat) +     &
                        factor(infactors,j) * valjl(l,jj) * jfacta(jj)

                END DO
             END DO
          END DO

          !---Printout J's:
          IF(ilng+i0==100 .AND. jlat+j0==120) THEN      
             WRITE(4,'(a)')
             WRITE(4,'(a)') ' fast-JX (5.7)----J-values----'
             !   	 write(4,'(a,i3,a,i3,a,f4.1,a,f7.2,a,2f7.2)')  &
             !   	   ' Step=',i,'  Day=',iday,' UT(hr)=',gmtau,'  SZA=',sza,  &
             !   	   '  LAT x LONG=',ydgrd(jlat),xdgrd(ilng)
             WRITE(4,'(1x,a,64(a7,2x))') 'L=  ',(jlabel(jx2spack(1,k)), k=1,16)
             DO l=jvl_-1,2,-1
                WRITE(4,'(i3,1p, 64e9.2)') l,(fast_JX_g(ngrid)%jphoto(j,l,ilng,jlat),j=1,16)
             ENDDO
             CALL flush(4)
          ENDIF

       END DO
    END DO
  END SUBROUTINE FastJX_driver


  !---------------------------------------------
  INTEGER FUNCTION reverse(pos,pmax)

    INTEGER, INTENT(IN) :: pos
    INTEGER, INTENT(IN) :: pmax

    reverse=pmax-pos+1

  END FUNCTION reverse


  !---------------------------------------------  
  SUBROUTINE get_solfac(imonth1,idate1,iyear1,time,solfac)

    INTEGER , INTENT(IN)  :: imonth1
    INTEGER , INTENT(IN)  :: idate1
    INTEGER , INTENT(IN)  :: iyear1
    REAL    , INTENT(IN)  :: time
    REAL    , INTENT(OUT) :: solfac

    INTEGER :: jday

!    INTEGER, EXTERNAL :: julday
    INTEGER :: julday

    REAL d0,d02
    ! solfac 
    jday = julday(imonth1,idate1,iyear1)
    jday = jday + NINT(time/86400.)
    d0 = 6.2831853 * float(jday-1) / 365.
    d02 = d0 * 2.
    solfac = 1.000110 + 0.034221 * COS (d0) + 0.001280 * SIN(d0)  &
         + 0.000719 * COS(d02) + 0.000077 * SIN(d02)
  END SUBROUTINE get_solfac


  !---------------------------------------------
  SUBROUTINE initialize_fast_JX(mmxp, mmyp, mmzp, ngrids)

    INTEGER, INTENT(IN) :: mmxp(ngrids)
    INTEGER, INTENT(IN) :: mmyp(ngrids) 
    INTEGER, INTENT(IN) :: mmzp(ngrids) 
    INTEGER, INTENT(IN) :: ngrids

    INTEGER :: ireaction,infactors,j,ng

    DO ireaction=1, nr_photo
       ! print*,'==============',ireaction,nfactors(ireaction),nr_photo

       DO infactors=1,nfactors(ireaction)
          !print*,'reac=',ireaction, ' factor=',infactors ;call flush(6)

          DO j=1,njval  !!! nratj
             ! print*,'1 ',trim(JReactionComp(infactors,ireaction)),' ',trim(jlabel(j));call flush(6)
             IF(TRIM(JReactionComp(infactors,ireaction))==TRIM(jlabel(j))) THEN
                !  print*,'2 ',trim(JReactionComp(infactors,ireaction)),' ',trim(jlabel(j)), j ;call flush(6)
	        jx2spack(infactors,ireaction)=j !; exit
             END IF
          END DO
       END DO
    END DO
    !do ireaction=1,nr_photo
    !   do infactors=1,nfactors(ireaction)
    !   if(infactors==1) print*,'reaction=',ireaction
    !   print*,'comp=',real(factor(infactors,ireaction)),trim(JReactionComp(infactors,ireaction))&
    !            ,jx2spack(infactors,ireaction)
    !enddo;enddo
    !          do j=1,njval  ; print*,'jlab=',trim(jlabel(j)),j,njval; enddo; call flush(6)


    ALLOCATE (fast_JX_g(ngrids))
    DO ng=1,ngrids
       ALLOCATE(fast_JX_g(ng)%jphoto (nr_photo,mmzp(ng),mmxp(ng),mmyp(ng)) )
       fast_JX_g(ng)%jphoto= 0.d0
    ENDDO

  END SUBROUTINE initialize_fast_JX
  !---------------------------------------------
  SUBROUTINE error_fastjx(fname)

    CHARACTER(len=*), INTENT(IN) :: fname

    PRINT*,' FAST JX ERROR : data file does not exist',fname
    STOP 3333
  END SUBROUTINE error_fastjx

  !------------------------------------------------------------------------
END MODULE FastJX
