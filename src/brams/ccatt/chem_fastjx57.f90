MODULE FastJX57

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rd_xxx, & ! Subroutine
            rd_mie, & ! Subroutine
            rd_js,  & ! Subroutine
            flint,  & ! Subroutine
            photoj    ! Subroutine


CONTAINS

  !  >>>>>>>>>>>>>>>>current code revised to JX ver 5.7 (3/07)<<<<<<<<<<<<
  ! version 5.7
  !     adds the new flux diagnostics (including heating rates)
  !        accurate fluxes for spherical atmos and SZA > 90 !
  !     recommend geometric delta-tau factor from 1.18 to 1.12 for more accurate
  !        heating rates (but more layers!)
  !     tuned and corrected to be almost flux conserving (1.e-5), except
  !        deep clouds, where diffusive flux is created (1.e-4)
  !     still needs to return to the original 1970-code for the block-tri solution
  !        after extensive profiling with F95 and 'modern' versions
  !        it was found that they are much more expensive!!!
  !     corrects typo in JAC(2000) fast-J paper on I+ (reflected from l.b.):
  !        I+(lb) = refl/(1+refl) * (4*Integ[j(lb)*mu*dmu] + mu0*Fdirect(lb))
  ! version 5.6 adds
  !      clean up problems with thick clouds does correct solar attenuation
  !        into cloud sub-layers and into the mid-point of the CTM level
  !      New calculated upward and downward FLUXES at each wavelength at TOP/BOT
  !      Correct deposition of solar flux in each CTM layer (spherical)
  !        awaits new diagnostics of the h's for heating rates.
  !      back to old matrix solver (UCI blocksolver and matinv-4)
  ! version 5.5 adds
  !      new code for generating and solving the block tri-diagonal scattering
  !           problem.  Uses single call to GEM and general 4x4 block-tri solver.
  ! version 5.3c adds
  !      calculates reflected UV-vis solar energy (relative to 1.0)
  !      new solar spectrum (J-O2 increases in strat by 10%, J-NO by 15+%)

  ! version 5.3b changes include:
  !      new data files for specral Xsection and mie-scattering.
  !      add sub-layers (JXTRA) to thick cloud/aerosol layers,
  !           sets up log-spaced sub-layers of increasing thickness ATAU
  !      correction 'b' does massive clean up of the linking code,
  !           now the only subroutine that has access to CTM arrays is PHOTOJ
  !           Also, the access to the cmn_JVdat.f is 'read-only' after init.
  !           This should enable safe openMP/MPI coding.

  ! common files and what they mean:
  !   parm_CTM.f  dimensions & params for code (CTM and fast-JX)
  !   parm_MIE.f  dimensions for mie code variables.
  !   cmn_metdat.f  CTM 3-D arrays, time of day, grid,  etc.
  !   cmn_JVdat.f   Xsects, Mie, etc., (initialized and then read-only)


  !     RD_JS(NJ1,NAMFIL):  Read labels of photo. rates, called once by INPHOT.
  !              COMMON BLOCKS: cmn_metdat.f, cmn_JVdat.f
  !              Input files: ratj.dat



  !<<<<<<<<<<<<<<<<<<<<<begin CTM-fastJX linking subroutines<<<<<<<<<<<<<<

  !     PHOTOJ(UTIME,IDAY,ILNG,JLAT, SZA,ZPJ)
  !              Gateway to fast-JX, Update the photolysis rates
  !              COMMON BLOCKS: cmn_metdat.f, cmn_JVdat.f

  !<<<<<<<<<<<<<<<<<<<<<begin core fast-J subroutines<<<<<<<<<<<<<<<<<<<<<
  !  N.B. all these need access to cmn_JVdat.f, but do NOT write into it.
  !           also have no need to access cmn_metdat.f

  !     JRATET(PPJ,TTJ,FFF, VALJL):  Calculate J-value, called by PTOTOJ.
  !              COMMON BLOCKS: cmn_JVdat.f

  !     JP_ATM(PPJ,TTJ,DDJ,ZZJ,ZHL,ZZHT,AER,ABX,ADX,JXTRA)
  !              print out atmosphere used in J-value calc.
  !              COMMON BLOCKS: cmn_JVdat.f


  !     RD_XXX(NJ1,NAMFIL):  Read wavelength bins, solar fluxes, Rayleigh
  !             parameters, TEM-dependent X-sections, called once by INPHOT.
  !              COMMON BLOCKS: cmn_JVdat.f
  !              Input files: FJX_spec.dat

  !     RD_MIE(NJ1,NAMFIL):  Set aerosols/cloud scattering, called once by INPHOT
  !              COMMON BLOCKS: cmn_JVdat.f
  !              Input files: FJX_scat.dat

  !     FUNCTION FLINT(TINT,T1,T2,T3,F1,F2,F3)

  !     SOLARZ(GMTIME,NDAY,YGRDJ,XGRDI, SZA,COSSZA,SOLFX)
  !              calc SZA and Solar Flux factor for given lat/lon/UT

  !     SPHERE(U0,RAD,ZHL,ZZHT,AMF,L1_):
  !              calculate spherical geometry, air-mass factors

  !     EXTRAL(AER,ADX,L1X,L2X,NX,JTAUMX,ATAU,ATAU0, JXTRA)
  !              add sub-layers (JXTRA) to thick cloud/aerosol layers

  !     OPMIE (KW,KM,WAVEL,ABX,AER,ADX,U0,RFLECT,AMF,
  !    &    JXTRA,FJACT,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0)
  !              calculate mean intensity (actinic) at each CTM levels
  !              calculate fluxes and deposition (heating rates)
  !              COMMON BLOCKS: cmn_JVdat.f

  !<<<<<<<<<<<<<<<<<<<<<<<begin core scattering subroutines<<<<<<<<<<<<<<<

  !      MIESCT (FJ,POMEGA,FZ,ZTAU,ZFLUX,ZREFL,ZU0,MFIT,ND)
  !            include 'parm_MIE.f' = dimension parameters

  !      BLKSLV (FJ,POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0,M,N,MFIT,ND)
  !              PARAMETER FILE: parm_MIE.f

  !      GEN (POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0,B,CC,AA,A,H,C1
  !            ,M,N,MFIT,ND,ID)
  !              PARAMETER FILE: parm_MIE.f

  !      LEGND0 (X,PL,N)

  !      MATIN4 (A)

  !      GAUSSP (N,XPT,XWT)

  !      EFOLD  (F0, F1, N, F)

  !
  !-----------------------------------------------------------------------
  SUBROUTINE photoj(sza,u0,frefl,solf,wrt,l_,l1_,l2_,jvl_,pp,tj,dm,do3, &
                    daer1,daer2,daer3,naer1,naer2,naer3,odcld,ncldx,zh, &
                    sa,valjl,njval,x_,atau0,atau,w_,jtaumx,naa,m_,n_,nw2, &
                    nw1,rad,szamax,zzht,wl,tqq,qrayl,fl,qo3,qo2,a_,qaa,paa, &
                    ssa,qqq,titlej,q1d)

    !-----------------------------------------------------------------------

    !  PHOTOJ is the gateway to fast-JX calculations:
    !	   only access to CTM 3-D GLOBAL arrays
    !	   sets up the 1-D column arrays for calculating J's


    !-----------------------------------------------------------------------
    !	AVGF   Attenuation of beam at each level for each wavelength
    !	FFF    Actinic flux at each desired level
    !	XQO2   Absorption cross-section of O2
    !	XQO3   Absorption cross-section of O3
    !-----------------------------------------------------------------------

    INTEGER          , INTENT(IN)    :: l_
    INTEGER          , INTENT(IN)    :: l1_
    INTEGER          , INTENT(IN)    :: l2_
    INTEGER          , INTENT(IN)    :: jvl_
    INTEGER          , INTENT(IN)    :: njval
    DOUBLE PRECISION , INTENT(IN)    :: pp(l1_)     !Pressure
    DOUBLE PRECISION , INTENT(IN)    :: sa          !Surface Albedo
    DOUBLE PRECISION , INTENT(IN)    :: tj(l1_)     !Temperature
    DOUBLE PRECISION , INTENT(IN)    :: dm(l1_)     !Air column for each model level (molecules.cm-2)
    DOUBLE PRECISION , INTENT(IN)    :: do3(l1_)    !Ozone column for each model level (molecules.cm-2)
    DOUBLE PRECISION , INTENT(IN)    :: zh(l1_)    
    DOUBLE PRECISION , INTENT(IN)    :: daer1(l1_)
    DOUBLE PRECISION , INTENT(IN)    :: daer2(l1_)
    DOUBLE PRECISION , INTENT(IN)    :: daer3(l1_)    
    INTEGER          , INTENT(IN)    :: naer1(l1_)
    INTEGER          , INTENT(IN)    :: naer2(l1_)
    INTEGER          , INTENT(IN)    :: naer3(l1_)
    INTEGER          , INTENT(IN)    :: ncldx(l1_)    
    DOUBLE PRECISION , INTENT(IN)    :: odcld(l1_)   
    DOUBLE PRECISION , INTENT(IN)    :: sza
    DOUBLE PRECISION , INTENT(IN)    :: solf
    LOGICAL          , INTENT(IN)    :: wrt
    DOUBLE PRECISION , INTENT(IN)    :: u0
    DOUBLE PRECISION , INTENT(INOUT) :: frefl  ! (DMK) alterado (OUT) para (INOUT)
    !DOUBLE PRECISION, INTENT(OUT) :: zpj(jvl_,jvn_) !2-D array of J's indexed to CTM chemistry!
    DOUBLE PRECISION , INTENT(OUT)   :: valjl(jvl_,njval) !2-D array of J_s returned by JRATET
    INTEGER          , INTENT(IN)    :: x_
    DOUBLE PRECISION , INTENT(IN)    :: atau0
    DOUBLE PRECISION , INTENT(IN)    :: atau
    INTEGER          , INTENT(IN)    :: w_
    INTEGER          , INTENT(IN)    :: jtaumx
    INTEGER          , INTENT(IN)    :: naa
    INTEGER          , INTENT(IN)    :: m_
    INTEGER          , INTENT(IN)    :: n_
    INTEGER          , INTENT(IN)    :: nw2
    INTEGER          , INTENT(IN)    :: nw1
    DOUBLE PRECISION , INTENT(IN)    :: rad
    DOUBLE PRECISION , INTENT(IN)    :: szamax
    DOUBLE PRECISION , INTENT(IN)    :: zzht
    DOUBLE PRECISION , INTENT(IN)    :: wl(w_)
    DOUBLE PRECISION , INTENT(INOUT) :: tqq(3,x_)
    DOUBLE PRECISION , INTENT(IN)    :: qrayl(w_+1)
    DOUBLE PRECISION , INTENT(IN)    :: fl(w_)
    DOUBLE PRECISION , INTENT(IN)    :: qo3(w_,3)
    DOUBLE PRECISION , INTENT(IN)    :: qo2(w_,3)
    INTEGER          , INTENT(IN)    :: a_
    DOUBLE PRECISION , INTENT(IN)    :: qaa(5,a_)
    DOUBLE PRECISION , INTENT(IN)    :: paa(8,5,a_)
    DOUBLE PRECISION , INTENT(IN)    :: ssa(5,a_)
    DOUBLE PRECISION , INTENT(IN)    :: qqq(w_,2,x_)
    CHARACTER(LEN=7) , INTENT(IN)    :: titlej(x_)  
    DOUBLE PRECISION , INTENT(IN)    :: q1d(w_,3)


    !-----------------------------------------------------------------------

    !--------key amtospheric data needed to solve plane-parallel J---------
    DOUBLE PRECISION, DIMENSION(5,l1_) :: aer
    INTEGER, DIMENSION(5,l1_) :: adx
    DOUBLE PRECISION, DIMENSION(l1_+1) :: abx, ttj,ddj,zzj,zhl
    DOUBLE PRECISION, DIMENSION(l1_+1) :: ppj
    INTEGER, DIMENSION(l2_+1) :: jxtra

    DOUBLE PRECISION	     :: rflect,frefs,frefi,fjtop,fjbot,fsbot
    DOUBLE PRECISION	     :: fjflx(l_),flxd(l1_),flxd0
    DOUBLE PRECISION	     :: amf(l1_+1,l1_+1)
    !------------key arrays AFTER solving for J's---------------------------
    DOUBLE PRECISION  :: fff(w_,jvl_)
    DOUBLE PRECISION  :: flxup(w_),flxdn(w_),dirup(w_),dirdn(w_)

    INTEGER :: i,k,l,km, ratio(w_)
    !TMP : temp fix
    !srf  DOUBLE PRECISION  :: avgf(l_),xqo3,xqo2	 ,wave, flint, ttt
    DOUBLE PRECISION  :: avgf(l1_),xqo3,xqo2	 ,wave, ttt

!    DOUBLE PRECISION, EXTERNAL :: flint

    !---flux/heating arrays (along with FJFLX,FLXD,FLXD0)
    DOUBLE PRECISION  :: flxj(l1_),ffx(w_,l1_),ffxnet(w_,8),ffx0,fxbot,fabot
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------

    fff(:,:) = 0.d0
    frefl = 0.d0
    frefs = 0.d0
    frefi = 0.d0

    !---check for dark conditions SZA > 98.0 deg => tan ht = 63 km
    !			   or	      99.		   80 km
    !LFR>   IF (sza > szamax) GO TO 99
    IF (sza > szamax) RETURN

    !---load the amtospheric column data
    DO l = 1,l1_
       ppj(l) = pp(l)   !LFR_>Calculated outside of this routine 
       ttj(l) = tj(l)
       ddj(l) = dm(l)
       zzj(l) = do3(l)
       !print*,'do3=',real(do3(l)),l;call flush(6)
    END DO
    ppj(l1_+1) = 0.d0

    !---calculate spherical weighting functions (AMF: Air Mass Factor)
    DO l = 1,l1_
       zhl(l) = zh(l)
    END DO
    zhl(l1_+1) = zhl(l1_) + zzht
    !print*,'ZHL=',zzht
    !-----------------------------------------------------------------------
    CALL sphere(u0,rad,zhl,amf,l1_)
    !-----------------------------------------------------------------------

    !---load the profiles of aerosols & clouds: treated as the same from here
    DO l = 1,l1_
       aer(1,l)  = daer1(l)  ! Opt. Depth aerosol 1 in layer L
       aer(2,l)  = daer2(l)  ! Opt. Depth aerosol 2 in layer L
       aer(3,l)  = daer3(l)  ! Opt. Depth aerosol 3 in layer L
       aer(4,l)  = odcld(l)  ! cloud Opt. Depth in L
       aer(5,l)  = 0.d0      ! save space for Rayleigh
       adx(1,l)  = naer1(l)  ! index for aerosol 1 at layer L
       adx(2,l)  = naer2(l)  ! index for aerosol 2 at layer L
       adx(3,l)  = naer3(l)  ! index for aerosol 3 at layer L
       adx(4,l)  = ncldx(l)  ! index for cloud at layer L
       adx(5,l)  = 1	  ! index for Rayleigh phase in L
    END DO

    DO l = 1,l1_
       DO i=1,4
          adx(i,l) = MIN(naa, MAX(0, adx(i,l)))
       END DO
    END DO

    !---Now given the aerosol+cloud OD/layer in visible (600 nm) can calculate
    !	   how to add additonal levels at top of clouds (now uses log spacing)
    !-----------------------------------------------------------------------
    CALL extral(aer,adx,l1_,l2_,n_,jtaumx,atau,atau0, jxtra)
    !-----------------------------------------------------------------------

    !---set surface reflectance
    rflect = sa
    rflect = MAX(0.d0,MIN(1.d0,rflect))

    !---Loop over all wavelength bins to calc mean actinic flux AVGF(L)
    DO k = nw1,nw2
       wave = wl(k)
       !---Pick nearest Mie wavelength, no interpolation--------------
       km=1  ! use 300 nm aerosol properties for <355 nm
       IF( wave > 355.d0 ) km=2  ! use 400 nm prop for 355-500 nm
       IF( wave > 500.d0 ) km=3
       IF( wave > 800.d0 ) km=4

       !---Loop over CTM layers L=1:L1_ = 1:L_+1,
       !	values at L1_=L_+1 are a pseudo layer above the top CTM layer (L_)
       DO l = 1,l1_
          ttt     = ttj(l)
          xqo3 = flint(ttt,tqq(1,2),tqq(2,2),tqq(3,2) ,qo3(k,1),qo3(k,2),qo3(k,3))

          xqo2 = flint(ttt,tqq(1,1),tqq(2,1),tqq(3,1) ,qo2(k,1),qo2(k,2),qo2(k,3))

          abx(l) = xqo3*zzj(l) + xqo2*ddj(l)*0.20948D0

          aer(5,l) = ddj(l)*qrayl(k)

          !  print*,'before opmie',ttt,zzj(l),ddj(l),xqo3


       END DO
       ! stop 333
       !-----------------------------------------------------------------------

       CALL opmie (km,abx,aer,adx,u0,rflect,amf,jxtra,avgf,fjtop,fjbot,fsbot, &
                   fjflx,flxd,flxd0,l1_,l2_,l_,m_,n_,atau0,qaa,a_,paa,ssa)

       !-----------------------------------------------------------------------

       !----direct(DIR) and diffuse(FLX) fluxes at top(UP) (solar = negative by convention)
       !----     also at bottom (DN), does not include diffuse reflected flux.
       flxup(k) =  fjtop
       dirup(k) = -flxd0
       flxdn(k) = -fjbot
       dirdn(k) = -fsbot

       DO l = 1,jvl_
          fff(k,l) = fff(k,l) + solf*fl(k)*avgf(l)
          !      print*,'XX',real(fff(k,l)),real(solf),real(fl(k)),real(avgf(l)),l,jvl_
       END DO
       frefi = frefi + solf*fl(k)*flxd0/wave
       frefl = frefl + solf*fl(k)*fjtop/wave
       frefs = frefs + solf*fl(k)/wave

       !---for each wavelength calculate the flux budget/heating rates:
       !  FLXD(L) = direct flux deposited in layer L  [approx = MU0*(F(L+1) -F(L)]
       !	       but for spherical atmosphere!
       !  FJFLX(L) = diffuse flux across top of layer L

       !---calculate divergence of diffuse flux in each CTM layer (& t-o-a)
       !---     need special fix at top and bottom:
       !---FABOT = total abs at L.B. &  FXBOT = net diffusive flux at L.B.
       fabot = (1.d0-rflect)*(fjbot+fsbot)
       fxbot = -fjbot + rflect*(fjbot+fsbot)
       flxj(1) = fjflx(1) - fxbot
       DO l=2,l_
          flxj(l) = fjflx(l) - fjflx(l-1)
       END DO
       flxj(l_+1) = fjtop - fjflx(l_)
       !---calculate net flux deposited in each CTM layer (direct & diffuse):
       ffx0 = 0.d0
       DO l=1,l1_
          ffx(k,l) = flxd(l) - flxj(l)
          ffx0 = ffx0 + ffx(k,l)
       END DO

       !  NB: the radiation level ABOVE the top CTM level is included in these budgets
       !	 these are the flux budget/heating terms for the column:
       !  FFXNET(K,1) = FLXD0	direct(solar) flux dep into atmos (spherical)
       !  FFXNET(K,2) = FSBOT	direct(solar) flux dep onto LB (surface)
       !  FFXNET(K,3) = FLXD0+FSBOT  TOTAL solar into atmopshere+surface
       !  FFXNET(K,4) = FJTOP	diffuse flux leaving top-of-atmos
       !  FFXNET(K,5) = FFX0 	diffuse flux absorbed in atmos
       !  FFXNET(K,6) = FABOT	total (dir+dif) absorbed at LB (surface)
       !	  these are surface fluxes to compare direct vs. diffuse:
       !  FFXNET(K,7) = FSBOT	direct flux dep onto LB (surface) - for srf diags
       !  FFXNET(K,8) = FJBOT	diffuse flux dep onto LB (surface)

       ffxnet(k,1) = flxd0
       ffxnet(k,2) = fsbot
       ffxnet(k,3) = flxd0+fsbot
       ffxnet(k,4) = fjtop
       ffxnet(k,5) = ffx0
       ffxnet(k,6) = fabot
       ffxnet(k,7) = fsbot
       ffxnet(k,8) = fjbot



    END DO       ! end loop over wavelength K

    frefl = frefl/frefs	   !calculate reflected flux (energy weighted)
    frefi = frefi/frefs

    !---NB UVB = 280-320 = bins 12:15, UVA = 320-400 = bins 16:17, VIS = bin 18 (++)

    !-----------------------------------------------------------------------
    CALL jratet(ppj,ttj,w_,fff, valjl, l1_,jvl_,nw2, nw1, njval, x_, tqq, &
         qqq,titlej,qo2,qo3,q1d)
    !-----------------------------------------------------------------------

    !---map the J-values from fast-JX onto ASAD ones (use JIND & JFACTA)
    !DO l = 1,jvl_
    !  DO j = 1,nratj
    !    IF (jind(j) > 0) THEN
    ! 	zpj(l,j) = valjl(l,jind(j))*jfacta(j)
    !    ELSE
    !	zpj(l,j) = 0.d0
    !    END IF
    !  END DO
    !END DO

    !---diagnostics that are NOT returned to the CTM code

    IF(wrt) THEN !If is to diagnostic print
       !-----------------------------------------------------------------------
       WRITE(6,*)'fast-JX-(5.7)----PHOTOJ internal print: Atmosphere----'

       CALL jp_atm(ppj,ttj,ddj,zzj,zhl,zzht,aer,adx,jxtra,l1_,l2_)

       !---PRINT SUMMARY of mean intensity, flux, heating rates:
       WRITE(6,*)
       WRITE(6,*)'fast-JX(5.7)----PHOTOJ internal print: Mean Intens----'
       WRITE(6,'(a,5f10.4)') ' SUMMARY fast-JX: albedo/SZA/u0/F-incd/F-refl/',  &
            rflect,sza,u0,frefi,frefl

       WRITE(6,'(a5,18i8)')   ' bin:',(k, k=nw2,nw1,-1)
       WRITE(6,'(a5,18f8.1)') ' wvl:',(wl(k), k=nw2,nw1,-1)
       WRITE(6,'(a)') ' ----  100000=Fsolar   MEAN INTENSITY per wvl bin'
       DO l = jvl_,1,-1
          DO k=nw1,nw2
             ratio(k) = (1.d5*fff(k,l)/fl(k))
          END DO
          WRITE(6,'(i3,2x,18i8)') l,(ratio(k),k=nw2,nw1,-1)
       END DO

       WRITE(6,*)
       WRITE(6,*)'fast-JX(5.7)----PHOTOJ internal print: Net Fluxes----'
       WRITE(6,'(a11,18i8)')   ' bin:',(k, k=nw2,nw1,-1)
       WRITE(6,'(a11,18f8.1)') ' wvl:',(wl(k), k=nw2,nw1,-1)
       !	 write(6,'(a11,18f8.4)') ' sol in atm',(FFXNET(K,1), K=NW2,NW1,-1)
       !	 write(6,'(a11,18f8.4)') ' sol at srf',(FFXNET(K,2), K=NW2,NW1,-1)
       WRITE(6,*) ' ---NET FLUXES--- '
       WRITE(6,'(a11,18f8.4)') ' sol TOTAL ',(ffxnet(k,3), k=nw2,nw1,-1)
       WRITE(6,'(a11,18f8.4)') ' dif outtop',(ffxnet(k,4), k=nw2,nw1,-1)
       WRITE(6,'(a11,18f8.4)') ' abs in atm',(ffxnet(k,5), k=nw2,nw1,-1)
       WRITE(6,'(a11,18f8.4)') ' abs at srf',(ffxnet(k,6), k=nw2,nw1,-1)
       WRITE(6,*) ' ---SRF FLUXES--- '
       WRITE(6,'(a11,18f8.4)') ' srf direct',(ffxnet(k,7), k=nw2,nw1,-1)
       WRITE(6,'(a11,18f8.4)') ' srf diffus',(ffxnet(k,8), k=nw2,nw1,-1)
       WRITE(6,'(2a)') '  ---NET ABS per layer:	 10000=Fsolar',  &
            '  [NB: values <0 = numerical error w/clouds or SZA>90, colm OK]'
       DO l = jvl_,1,-1
          DO k=nw1,nw2
             ratio(k) = 1.d5*ffx(k,l)
          END DO
          WRITE(6,'(i9,2x,18i8)') l,(ratio(k),k=nw2,nw1,-1)
       END DO
    END IF
    !-----------------------------------------------------------------------
  END SUBROUTINE photoj

  !<<<<<<<<<<<<<<<<<<<<<begin core fast-J subroutines<<<<<<<<<<<<<<<<<<<<<


  !-----------------------------------------------------------------------  
  SUBROUTINE jratet(ppj,ttj,w_,fff, valjl,l1_,jvl_,nw2,nw1,njval,x_,tqq, &
                    qqq,titlej,qo2,qo3,q1d)
    !-----------------------------------------------------------------------
    ! in:
    !	   PPJ(L1_+1) = pressure profile at edges
    !	   TTJ(L1_) = = temperatures at mid-level
    !	   FFF(K=1:NW, L=1:JVL_) = mean actinic flux
    ! out:
    !	   VALJ(JVL_,NJVAL)  JVL_ = no of levels
    !-----------------------------------------------------------------------

    !LFR>       include 'parm_CTM.f'
    !LFR>       include 'cmn_JVdat.f'
    INTEGER          , INTENT(IN)    :: l1_
    INTEGER          , INTENT(IN)    :: jvl_
    DOUBLE PRECISION , INTENT(IN)    :: ppj(l1_+1)
    DOUBLE PRECISION , INTENT(IN)    :: ttj(l1_)
    INTEGER          , INTENT(IN)    :: w_
    DOUBLE PRECISION , INTENT(IN)    :: fff(w_,jvl_)
    DOUBLE PRECISION , INTENT(OUT)   :: valjl(jvl_,njval)
    INTEGER          , INTENT(IN)    :: nw2
    INTEGER          , INTENT(IN)    :: nw1
    INTEGER          , INTENT(IN)    :: njval
    INTEGER          , INTENT(IN)    :: x_
    DOUBLE PRECISION , INTENT(INOUT) :: tqq(3,x_)
    DOUBLE PRECISION , INTENT(IN)    :: qqq(w_,2,x_)
    CHARACTER(LEN=7) , INTENT(IN)    :: titlej(x_)  
    DOUBLE PRECISION , INTENT(IN)    :: qo2(w_,3)
    DOUBLE PRECISION , INTENT(IN)    :: qo3(w_,3)
    DOUBLE PRECISION , INTENT(IN)    :: q1d(w_,3)

!    DOUBLE PRECISION, EXTERNAL  :: flint_ 	    ! external function for X-sections
    DOUBLE PRECISION  :: valj(x_)	    ! temp for call J's at one L
    DOUBLE PRECISION  :: qo2tot, qo3tot, qo31dy, qo31d, qqqt, tfact
    DOUBLE PRECISION  :: tt,pp,dd,tt200,tfaca,tfac0,tfac1,tfac2,qqqa,qq2,qq1a,qq1b
    INTEGER :: j,k,l, iv

    DO l = 1,jvl_    ! master loop over layer = L

       !---need temperature and density (for some quantum yields):
       !---in this case the Pressures PPJ are defined at the boundaries,
       !---  	      Temperatures in the middle of each layer
       tt   = ttj(l)
       pp  = (ppj(l)+ppj(l+1))*0.5D0
       IF (l == 1) pp = ppj(1)
       dd = 7.24E18*pp/tt

       DO j = 1,njval
          valj(j) = 0.d0
       END DO

       DO k = nw1,nw2		      ! Using model 'T's here
          qo3tot = flint(tt,tqq(1,2),tqq(2,2),tqq(3,2) ,qo3(k,1),qo3(k,2),qo3(k,3))
          qo2tot = flint(tt,tqq(1,1),tqq(2,1),tqq(3,1) ,qo2(k,1),qo2(k,2),qo2(k,3))
          qo31dy = flint(tt,tqq(1,3),tqq(2,3),tqq(3,3) ,q1d(k,1),q1d(k,2),q1d(k,3))
          qo31d  = qo31dy*qo3tot
          valj(1) = valj(1) + qo2tot*fff(k,l)
          valj(2) = valj(2) + qo3tot*fff(k,l)
          valj(3) = valj(3) + qo31d*fff(k,l)
       END DO

       DO j = 4,njval

          IF (tqq(2,j) > tqq(1,j)) THEN
             tfact = MAX(0.d0,MIN(1.d0,(tt-tqq(1,j))/(tqq(2,j)-tqq(1,j))))
          ELSE
             tfact = 0.d0
          END IF

          DO k = nw1,nw2
             qqqt	= qqq(k,1,j) + (qqq(k,2,j) - qqq(k,1,j))*tfact
             valj(j) = valj(j) + qqqt*fff(k,l)
          END DO

          ! #52 Methylvinyl ketone   'MeVK  '	  q(M) = 1/(1 + 1.67e-19*[M])
          IF (titlej(j) == 'MeVK  ') THEN
             valj(j) = valj(j)/(1.0 + 1.67E-19*dd)
          END IF
          ! #55 Methylethyl ketone   MEKeto	q(M) = 1/(1 + 2.0*[M/2.5e19])
          IF (titlej(j) == 'MEKeto') THEN
             valj(j) = valj(j)/(1.0 + 0.80E-19*dd)
          END IF
          ! #57 Methyl glyoxal       MGlyxl	q(M) = 1/(1 + 4.15*[M/2.5E19])
          IF (titlej(j) == 'MGlyxl') THEN
             valj(j) = valj(j)/(1.0 + 1.66E-19*dd)
          END IF

       END DO

       IF (titlej(njval-1) == 'Acet-a') THEN
          !--------------J-ref v8.3 includes Blitz ACETONE q-yields--------------
          !---Acetone is a special case:   (as per Blitz et al GRL, 2004)
          !---     61 = NJVAL-1 = J1(acetone-a) ==> CH3CO + CH3
          !---     62 = NJVAL	= J2(acetone-b) ==> CH3 + CO + CH3
          valj(njval-1) = 0.d0
          valj(njval)   = 0.d0
          !---IV=NJVAL-1 = Xsect (total abs) for Acetone - pre-calc Temp interp factors
          iv    = njval-1
          tfaca = (tt-tqq(1,iv))/(tqq(2,iv)-tqq(1,iv))
          tfaca = MAX(0.d0, MIN(1.d0, tfaca))
          !---IV=NJVAL = Q2 for Acetone=>(2), specifically designed for quadratic interp.
          !---      but force to Q2=0 by 210K
          iv    = njval
          tfac0 = ( (tt-tqq(1,iv))/(tqq(2,iv)-tqq(1,iv)) )**2
          IF (tt < tqq(1,iv)) THEN
             tfac0 = (tt - 210.d0)/(tqq(1,iv)-210.d0)
          END IF
          tfac0 = MAX(0.d0, MIN(1.d0, tfac0))
          !---IV=NJVAL+1 = Q1A for Acetone => (1), allow full range of T = 200K-300K
          iv    = njval+1
          tt200 = MIN(300.d0, MAX(200.d0, tt))
          tfac1 = (tt200-tqq(1,iv))/(tqq(2,iv)-tqq(1,iv))
          !---IV=NJVAL+2 = Q1B for Acetone => (1)
          iv    = njval+2
          tfac2 = (tt200-tqq(1,iv))/(tqq(2,iv)-tqq(1,iv))

          !---now integrate over wavelengths
          DO k = nw1,nw2
             !---NJVAL-1 = Xsect (total abs) for Acetone
             iv   = njval-1
             qqqa = qqq(k,1,iv) + (qqq(k,2,iv)-qqq(k,1,iv))*tfaca
             !---NJVAL   = Q2 for Acetone=>(2), specifically designed for quadratic interp.
             iv   = njval
             qq2  = qqq(k,1,iv) + (qqq(k,2,iv)-qqq(k,1,iv))*tfac0
             IF (tt < tqq(1,iv)) THEN
                qq2 = qqq(k,1,iv)*tfac0
             END IF
             !---NJVAL+1 = Q1A for Acetone => (1), allow full range of T = 200K-300K
             iv   = njval+1
             qq1a = qqq(k,1,iv) + (qqq(k,2,iv)-qqq(k,1,iv))*tfac1
             !---NJVAL+2 = Q1B for Acetone => (1)   ! scaled to [M]=2.5e19
             iv   = njval+2
             qq1b = qqq(k,1,iv) + (qqq(k,2,iv)-qqq(k,1,iv))*tfac2
             qq1b = qq1b*4.d-20
             !---J(61)
             valj(njval-1) = valj(njval-1)  &
                  + fff(k,l)*qqqa*(1.d0-qq2)/(qq1a + qq1b*dd)
             !---J(62)
             valj(njval) = valj(njval) + fff(k,l)*qqqa*qq2

          END DO	!K
          !-----------end v-8.3 includes Blitz ACETONE q-yields--------------
       END IF

       !----Load array of J-values in native order, need to be indexed/scaled
       !    by ASAD-related code later: ZPJ(L,JJ) = VALJL(L,JIND(JJ))*JFACTA(JJ)
       DO j=1,njval
          valjl(l,j) = valj(j)
       END DO

    END DO    ! master loop over L=1,JVL_

  END SUBROUTINE jratet


  !-----------------------------------------------------------------------
  SUBROUTINE jp_atm(ppj,ttj,ddj,zzj,zhl,zzht1,aer,adx,jxtra,l1_,l2_)
  !-----------------------------------------------------------------------

    !LFR>       include 'parm_CTM.f'

    INTEGER                           , INTENT(IN)    :: l1_
    INTEGER                           , INTENT(IN)    :: l2_
    DOUBLE PRECISION, DIMENSION(l1_+1), INTENT(INOUT) :: ppj
    DOUBLE PRECISION, DIMENSION(l1_+1), INTENT(INOUT) :: ttj
    DOUBLE PRECISION, DIMENSION(l1_+1), INTENT(IN)    :: ddj
    DOUBLE PRECISION, DIMENSION(l1_+1), INTENT(IN)    :: zzj
    DOUBLE PRECISION, DIMENSION(l1_+1), INTENT(IN)    :: zhl
    DOUBLE PRECISION                  , INTENT(IN)    :: zzht1
    DOUBLE PRECISION, DIMENSION(5,l1_), INTENT(IN)    :: aer
    INTEGER, DIMENSION(5,l1_)         , INTENT(INOUT) :: adx
    INTEGER, DIMENSION(l2_+1)         , INTENT(IN)    :: jxtra
    !-----------------------------------------------------------------------
    !--------key amtospheric data needed to solve plane-parallel J---------  
    !-----------------------------------------------------------------------
    INTEGER :: i,l
    DOUBLE PRECISION  :: col(4),colo2,colo3,zkm,delz,ztop

    WRITE(6,'(4a)') '   L z(km)	  p	 T   ',  &
         '    d(air)   d(O3)','  col(O2)  col(O3)  ndx colm(aer/cld)',  &
         ' added cld lvls: top/bot CTM lyr=>'

    colo2 = 0.d0
    colo3 = 0.d0
    DO i=1,4
       col(i) = 0.d0
    END DO
    ztop = zhl(l1_) + zzht1

    DO l = l1_,1,-1

       DO i=1,4
          col(i) = col(i) + aer(i,l)
       END DO
       colo2 = colo2 + ddj(l)*0.20948D0
       colo3 = colo3 + zzj(l)
       delz = ztop-zhl(l)
       ztop = zhl(l)
       zkm = zhl(l)*1.d-5

       WRITE(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,4(i3,e9.2),2i3)')  &
            l,zkm,ppj(l),ttj(l),ddj(l)/delz,zzj(l)/delz, colo2,colo3,  &
            (adx(i,l),col(i), i=1,4), jxtra(l+l),jxtra(l+l-1)
    END DO

  END SUBROUTINE jp_atm


  !-----------------------------------------------------------------------
  SUBROUTINE rd_xxx(nj1,namfil,njval,nw1,nw2,title0,x_,titlej,tqq,w_,wl,fl, &
                    qrayl,qo2,qo3,q1d,qqq)

    !-----------------------------------------------------------------------
    !  Read in wavelength bins, solar fluxes, Rayleigh parameters,
    !	 T-dependent X-sections.

    !  >>>current code revised to JPL-02 ver 8.5 (5/05)<<<<<

    !-----------------------------------------------------------------------
    !	NAMFIL   Name of spectral data file (j2_spec.dat) >> j2 for fast-J2
    !	NJ1	 Channel number for reading data file

    !	NJVAL	 Number of species to calculate J-values for
    !	NWWW	 Number of wavelength bins, from 1:NWWW
    !	WBIN	 Boundaries of wavelength bins
    !	WL	 Centres of wavelength bins - 'effective wavelength'
    !	FL	 Solar flux incident on top of atmosphere (cm-2.s-1)
    !	QRAYL	 Rayleigh parameters (effective cross-section) (cm2)
    !	QO2	 O2 cross-sections
    !	QO3	 O3 cross-sections
    !	Q1D	 O3 => O(1D) quantum yield
    !	TQQ	 Temperature for supplied cross sections
    !	QQQ	 Supplied cross sections in each wavelength bin (cm2)
    !-----------------------------------------------------------------------

    !LFR>       include 'parm_CTM.f'
    !LFR>       include 'cmn_JVdat.f'

    INTEGER          , INTENT(IN)  :: nj1
    CHARACTER(LEN=*) , INTENT(IN)  ::  namfil

    INTEGER          , INTENT(OUT) :: njval
    INTEGER          , INTENT(OUT) :: nw1
    INTEGER          , INTENT(OUT) :: nw2
    CHARACTER(LEN=78), INTENT(OUT) :: title0
    INTEGER          , INTENT(IN)  :: x_
    CHARACTER(LEN=7) , INTENT(OUT) :: titlej(x_)  
    DOUBLE PRECISION , INTENT(OUT) :: tqq(3,x_)
    INTEGER          , INTENT(IN)  :: w_
    DOUBLE PRECISION , INTENT(OUT) :: wl(w_)
    DOUBLE PRECISION , INTENT(OUT) :: fl(w_)
    DOUBLE PRECISION , INTENT(OUT) :: qrayl(w_+1)
    DOUBLE PRECISION , INTENT(OUT) :: qo2(w_,3)
    DOUBLE PRECISION , INTENT(OUT) :: qo3(w_,3)
    DOUBLE PRECISION , INTENT(OUT) :: q1d(w_,3)
    DOUBLE PRECISION , INTENT(OUT) :: qqq(w_,2,x_)


    CHARACTER(LEN=7) :: titlej2 ! (DMK) scratch
    CHARACTER(LEN=7) :: titlej3 ! (DMK) scratch


    INTEGER :: i, j, iw, nqqq, nwww

    tqq(:,:) = 0.d0

    !----------spectral data----set for new format data J-ver8.3------------------
    !	    note that NJVAL = # J-values, but NQQQ (>NJVAL) = # Xsects read in
    !	    for 2005a data, NJVAL = 62 (including a spare XXXX) and
    !		 NQQQ = 64 so that 4 wavelength datasets read in for acetone
    !	    note NQQQ is not used outside this subroutine!

    OPEN (nj1,FILE=namfil,STATUS='old',FORM='formatted')
    READ (nj1,100) title0
    READ (nj1,101) njval,nqqq, nwww,nw1,nw2
    IF (njval > x_ .OR. nqqq > x_) THEN
       WRITE(6,201) njval,x_
       STOP
    END IF
    WRITE(6,'(1X,A)') title0
    !----J-values:  1=O2, 2=O3P,3=O3D 4=readin Xsects
    READ (nj1,102) (wl(iw),iw=1,nwww)
    READ (nj1,102) (fl(iw),iw=1,nwww)
    READ (nj1,102) (qrayl(iw),iw=1,nwww)

    !---Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps)
    READ (nj1,103) titlej(1),tqq(1,1), (qo2(iw,1),iw=1,nwww)
    READ (nj1,103) titlej2,  tqq(2,1), (qo2(iw,2),iw=1,nwww)
    READ (nj1,103) titlej3,  tqq(3,1), (qo2(iw,3),iw=1,nwww)

    READ (nj1,103) titlej(2),tqq(1,2), (qo3(iw,1),iw=1,nwww)
    READ (nj1,103) titlej2,  tqq(2,2), (qo3(iw,2),iw=1,nwww)
    READ (nj1,103) titlej3,  tqq(3,2), (qo3(iw,3),iw=1,nwww)

    READ (nj1,103) titlej(3),tqq(1,3), (q1d(iw,1),iw=1,nwww)
    READ (nj1,103) titlej2,  tqq(2,3), (q1d(iw,2),iw=1,nwww)
    READ (nj1,103) titlej3,  tqq(3,3), (q1d(iw,3),iw=1,nwww)

    DO j = 1,3
       WRITE(6,200) titlej(j),(tqq(i,j),i=1,3)
    END DO

    !---Read remaining species:  X-sections at 2 T_s
    DO j = 4,nqqq
       READ (nj1,103) titlej(j),tqq(1,j),(qqq(iw,1,j),iw=1,nwww)
       READ (nj1,103) titlej2,  tqq(2,j),(qqq(iw,2,j),iw=1,nwww)
       WRITE(6,200) titlej(j),(tqq(i,j),i=1,2)
    END DO

    !  Reset the titles for NJVAL-1 & NJVAL to be the two acetone J_s
    !   61: C3H6O  = Acet-a     (CH3CO + CH3)
    !   62: Q2-Ac  = Acet-b     (CH3 + CO + CH3)

    titlej(njval-1) = 'Acet-a'
    titlej(njval)   = 'Acet-b'

    CLOSE(nj1)

100 FORMAT(a)
101 FORMAT(10X,5I5)
102 FORMAT(10X,    6E10.3/(10X,6E10.3)/(10X,6E10.3))
103 FORMAT(a7,f3.0,6E10.3/(10X,6E10.3)/(10X,6E10.3))
200 FORMAT(1X,' x-sect:',a10,3(3X,f6.2))
201 FORMAT(' Number of x-sections supplied to Fast-J2: ',i3,/,  &
         ' Maximum number allowed (X_) only set to: ',i3, ' - increase in cmn_jv.f')

  END SUBROUTINE rd_xxx

  !-----------------------------------------------------------------------  
  SUBROUTINE rd_mie(nj1,namfil,atau,atau0,jtaumx,naa,title0,a_,qaa,ssa,paa)
    !-----------------------------------------------------------------------
    !-------aerosols/cloud scattering data set for fast-JX (ver 5.3+)
    !  >>>>>>>>>>>>>>>>spectral data rev to J-ref ver8.5 (5/05)<<<<<<<<<<<<
    !-----------------------------------------------------------------------
    !	NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
    !	NJ1	 Channel number for reading data file
    !	NAA	 Number of categories for scattering phase functions
    !	QAA	 Aerosol scattering phase functions
    !	NK	 Number of wavelengths at which functions supplied (set as 4)
    !	WAA	 Wavelengths for the NK supplied phase functions
    !	PAA	 Phase function: first 8 terms of expansion
    !	RAA	 Effective radius associated with aerosol type
    !	SSA	 Single scattering albedo
    !-----------------------------------------------------------------------

    !LFR>       include 'parm_CTM.f'
    !LFR>       include 'cmn_JVdat.f'

    INTEGER          , INTENT(IN)  :: nj1
    CHARACTER(LEN=*) , INTENT(IN)  :: namfil

    DOUBLE PRECISION , INTENT(OUT) :: atau
    DOUBLE PRECISION , INTENT(OUT) :: atau0
    INTEGER          , INTENT(OUT) :: jtaumx
    INTEGER          , INTENT(OUT) :: naa
    CHARACTER(LEN=78), INTENT(OUT) :: title0
    INTEGER          , INTENT(IN)  :: a_
    DOUBLE PRECISION , INTENT(OUT) :: qaa(5,a_)
    DOUBLE PRECISION , INTENT(OUT) :: ssa(5,a_)
    DOUBLE PRECISION , INTENT(OUT) :: paa(8,5,a_)


    INTEGER :: i, j, k

    CHARACTER(LEN=20) :: titlea(a_) ! (DMK) scratch
    DOUBLE PRECISION  :: waa(5,a_)  ! (DMK) scratch
    DOUBLE PRECISION  :: raa(5,a_)  ! (DMK) scratch
    

    OPEN (nj1,FILE=namfil,STATUS='old',FORM='formatted')

    READ (nj1,'(i2,a78)') naa,title0
    READ (nj1,'(5x,i5,2f10.5)') jtaumx,atau,atau0
    WRITE(6,'(a,2f9.5,i5)') ' ATAU/ATAU0/JMX',atau,atau0,jtaumx
    READ (nj1,*)

    DO j = 1,naa
       READ (nj1,'(3x,a20)') titlea(j)
       DO k = 1,4     ! Fix number of aerosol wavelengths at 4
          READ (nj1,'(f5.0,f8.1,f7.3,f8.4,f7.3,7f6.3)')  &
               waa(k,j),qaa(k,j),raa(k,j),ssa(k,j),(paa(i,k,j),i=1,8)
       END DO
    END DO

    CLOSE(nj1)

    WRITE(6,*) 'Aerosol phase functions & wavelengths'
    WRITE(6,*) title0
    DO j=1,naa
       WRITE(6,'(1x,A8,I2,A,9F8.1)') titlea(j),j,'  wavel=',(waa(k,j),k=1,4)
       WRITE(6,'(9x,I2,A,9F8.4)') j,'  Qext =',(qaa(k,j),k=1,4)
    END DO


  END SUBROUTINE rd_mie


  !-----------------------------------------------------------------------  
  DOUBLE PRECISION FUNCTION flint (tint,t1,t2,t3,f1,f2,f3)
    !-----------------------------------------------------------------------
    !  Three-point linear interpolation function
    !-----------------------------------------------------------------------

    DOUBLE PRECISION, INTENT(INOUT) :: tint
    DOUBLE PRECISION, INTENT(INOUT) :: t1
    DOUBLE PRECISION, INTENT(INOUT) :: t2
    DOUBLE PRECISION, INTENT(INOUT) :: t3
    DOUBLE PRECISION, INTENT(IN)    :: f1
    DOUBLE PRECISION, INTENT(IN)    :: f2
    DOUBLE PRECISION, INTENT(IN)    :: f3

    IF (tint <= t2)  THEN
       IF (tint <= t1)  THEN
          flint = f1
       ELSE
          flint = f1 + (f2 - f1)*(tint -t1)/(t2 -t1)
       END IF
    ELSE
       IF (tint >= t3)  THEN
          flint = f3
       ELSE
          flint = f2 + (f3 - f2)*(tint -t2)/(t3 -t2)
       END IF
    END IF

  END FUNCTION flint


  !-----------------------------------------------------------------------  
  SUBROUTINE sphere(gmu,rad,zhl,amf,l1_)
    !-----------------------------------------------------------------------
    !  Calculation of spherical geometry; derive tangent heights, slant path
    !  lengths and air mass factor for each layer. Not called when
    !  SZA > 98 degrees.  Beyond 90 degrees, include treatment of emergent
    !  beam (where tangent height is below altitude J-value desired at).
    !-----------------------------------------------------------------------
    ! in:
    !	GMU	= MU0 = cos(solar zenith angle)
    !	RAD	radius of Earth mean sea level (cm)
    !	ZHL(L)  height (cm) of the bottome edge of CTM level L
    !	ZZHT	scale height (cm) used above top of CTM (ZHL(L_+1)
    !	L1_	dimension of CTM = levels +1
    ! out:
    !	AMF(I,J) = air mass factor for CTM level I for sunlight reaching J
    !-----------------------------------------------------------------------

    DOUBLE PRECISION , INTENT(IN)  :: gmu
    DOUBLE PRECISION , INTENT(IN)  :: rad
    DOUBLE PRECISION , INTENT(IN)  :: zhl(l1_+1)
    DOUBLE PRECISION , INTENT(OUT) :: amf(l1_+1,l1_+1)
    INTEGER          , INTENT(IN)  :: l1_

    !	RZ	Distance from centre of Earth to each point (cm)
    !	RQ	Square of radius ratios
    !	TANHT	Tangent height for the current SZA
    !	XL	Slant path between points

    INTEGER :: i, j, ii
    DOUBLE PRECISION  :: xmu1,xmu2,xl,diff,tanht,rz(l1_+1),rq(l1_+1)

    !  Inlined air mass factor function for top of atmosphere
    !	 AIRMAS(Ux,H) = (1.0d0+H)/SQRT(Ux*Ux+2.0d0*H*(1.0d0-
    !	&	  0.6817d0*EXP(-57.3d0*abs(Ux)/SQRT(1.0d0+5500.d0*H))/
    !	&					      (1.0d0+0.625d0*H)))
    !  Use function and scale height to provide AMF above top of model
    !	   ZBYR  = ZZHT/RAD
    !	   AMF(L_+1,J)  = AIRMAS(XMU1,ZBYR)

    !--- must have top-of-atmos (NOT top-of-CTM) defined
    !	 ZHL(L1_+1) = ZHL(L1_) + ZZHT

    rz(1) = rad + zhl(1)

    DO ii = 2,l1_+1
       rz(ii)   = rad + zhl(ii)
       rq(ii-1) = (rz(ii-1)/rz(ii))**2
    END DO
    IF (gmu < 0.0D0)  THEN
       tanht = rz(1)/DSQRT(1.0D0-gmu**2)
    ELSE
       tanht = 0.d0
    END IF

    !  Go up from the surface calculating the slant paths between each level
    !  and the level above, and deriving the appropriate Air Mass Factor
    DO  j = 1,l1_+1

       DO i = 1,l1_+1
          amf(i,j) = 0.d0
       END DO
       !  Air Mass Factors all zero if below the tangent height
       IF (rz(j) < tanht) CYCLE
       !  Ascend from layer J calculating AMFs
       xmu1 = ABS(gmu)
       DO i = j,l1_
          xmu2     = DSQRT(1.0D0 - rq(i)*(1.0D0-xmu1**2))
          xl       = rz(i+1)*xmu2 - rz(i)*xmu1
          amf(i,j) = xl / (rz(i+1)-rz(i))
          !print*,'AMF=',real(amf(i,j)),xl,rz(i+1),rz(i),i
          xmu1     = xmu2
       END DO
       !--fix above top-of-atmos (L=L1_+1), must set DTAUX(L1_+1)=0
       amf(l1_+1,j) = 1.d0

       !  Twilight case - Emergent Beam, calc air mass factors below layer
       IF (gmu >= 0.0D0) CYCLE

       !  Descend from layer J
       xmu1       = ABS(gmu)
       DO ii = j-1,1,-1
          diff	  = rz(ii+1)*SQRT(1.0D0-xmu1**2)-rz(ii)
          IF (ii == 1)  diff = MAX(diff,0.d0)   ! filter
          !  Tangent height below current level - beam passes through twice
          IF (diff < 0.0D0)  THEN
             xmu2	  = SQRT(1.0D0 - (1.0D0-xmu1**2)/rq(ii))
             xl	  = ABS(rz(ii+1)*xmu1-rz(ii)*xmu2)
             amf(ii,j) = 2.d0*xl/(rz(ii+1)-rz(ii))
             xmu1	  = xmu2
             !  Lowest level intersected by emergent beam
          ELSE
             xl	  = rz(ii+1)*xmu1*2.0D0
             amf(ii,j) = xl/(rz(ii+1)-rz(ii))
             CYCLE
          END IF
       END DO

    END DO
  END SUBROUTINE sphere

  !-----------------------------------------------------------------------
  SUBROUTINE extral(aer,adx,l1x,l2x,nx,jtaumx,atau,atau0,jxtra)
    !-----------------------------------------------------------------------  
    !    new version 5.3, add sub-layers (JXTRA) to thick cloud/aerosol layers
    !    this version sets up log-spaced sub-layers of increasing thickness ATAU

    !	AER(5,L=1:L1X) = Optical Depth in layer L (general visible OD)
    !	   This is not used in the calculation of J's but in calculating
    !	   the number in levels to insert in each layer L
    !	   Set for log-spacing of tau levels, increasing top-down.
    !	AER(1:4,L) = aerosol+cloud OD (up to 4 types, index type = ADX(1:4,L)
    !	AER(5,L) = Rayleigh-not used here

    !	N.B. the TTAU, etc caluclate here are not really used elsewhere

    !---The log-spacing parameters have been tested for convergence and chosen
    !---  to be within 0.5% for ranges OD=1-500, rflect=0-100%, mu0=0.1-1.0
    !---  use of ATAU = 1.18 and min = 0.01, gives at most +135 pts for OD=100
    !-----------------------------------------------------------------------

    DOUBLE PRECISION , INTENT(IN)    :: aer(5,l1x) !cloud+3aerosol OD in each layer
    INTEGER          , INTENT(IN)    :: adx(5,l1x)
    INTEGER          , INTENT(IN)    :: l1x !index of cloud/aerosol
    INTEGER          , INTENT(IN)    :: l2x !index of cloud/aerosol
    INTEGER          , INTENT(IN)    :: nx
    INTEGER          , INTENT(IN)    :: jtaumx !index of cloud/aerosol
    DOUBLE PRECISION , INTENT(IN)    :: atau
    DOUBLE PRECISION , INTENT(IN)    :: atau0
    INTEGER          , INTENT(INOUT) :: jxtra(l2x+1) !number of sub-layers to be added  ! (DMK) alterado (OUT) para (INOUT)

    INTEGER :: jtotl,i,l,l2
    DOUBLE PRECISION  :: dtaux(l1x),ttau(l2x+1),dtauj, atau1,atauln,ataum,ataun1

    !---Reinitialize arrays
    ttau(:)  = 0.d0
    jxtra(:) = 0

    !---Set up total optical depth over each CTM level (L=1:L1X)
    !---   DTAUX(L) = bulk properties (OD) of each CTM layer

    DO l = 1,l1x
       !  Total optical depth from all elements I=1:4 = clouds  + 3 aerosol types
       !    Assume here that AER(1:4, 1:L1X) is 'visible' Optical Depth = 600 nm
       !    do NOT count if ADX = 0
       dtaux(l)   = 0.d0
       DO i = 1,4
          IF (adx(i,l) > 0) dtaux(l) = dtaux(l) + aer(i,l)
       END DO
    END DO

    !---combine these edge- and mid-layer points into grid of size:
    !---  	    L2X+1 = 2*L1X+1 = 2*L_+3
    !---calculate column optical depths above each level, TTAU(1:L2X+1)
    !---      note that TTAU(L2X+1)=0 and TTAU(1)=total OD

    !---Divide thick layers to achieve better accuracy in the scattering code
    !---In the original fast-J, equal sub-layers were chosen, this is wasteful
    !---and this new code (ver 5.3) uses log-scale:
    !---        Each succesive layer (down) increase thickness by ATAU > 1
    !---        e.g., if ATAU = 2, a layer with OD = 15 could be divided into
    !---        4 sub-layers with ODs = 1 - 2 - 4 - 8
    !---The key parameters are:
    !---        ATAU = factor increase from one layer to the next
    !---        ATAUMN = the smallest OD layer desired
    !---        JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)
    !---These are hardwired below, can be changed, but have been tested/optimized

    atau1  = atau - 1.d0
    atauln = LOG(atau)
    ttau(l2x+1)  = 0.0D0
    DO l2 = l2x,1,-1
       l	      = (l2+1)/2
       dtauj     = 0.5D0 * dtaux(l)
       ttau(l2)  = ttau(l2+1) + dtauj
       !---Now compute the number of log-spaced sub-layers to be added in
       !---   the interval TTAU(L2) > TTAU(L2+1)
       !---The objective is to have successive TAU-layers increasing by factor ATAU >1
       !---the number of sub-layers + 1
       IF (ttau(l2) < atau0) THEN
          jxtra(l2) = 0
       ELSE
          ataum    = MAX(atau0, ttau(l2+1))
          ataun1 = LOG(ttau(l2)/ataum) / atauln
          jxtra(l2) = MIN(jtaumx, MAX(0, INT(ataun1 - 0.5D0)))
       END IF
    END DO

    !---check on overflow of arrays, cut off JXTRA at lower L if too many levels
    jtotl    = l2x + 2
    DO l2 = l2x,1,-1
       jtotl  = jtotl + jxtra(l2)
       IF (jtotl > nx/2)  THEN
          WRITE(6,'(A,2I5,F9.2)') 'N_/L2_/L2-cutoff JXTRA:',nx,l2x,l2
          DO l = l2,1,-1
             jxtra(l) = 0
          END DO
          EXIT
       END IF
    END DO

  END SUBROUTINE extral


  !-----------------------------------------------------------------------  
  SUBROUTINE opmie (km,abx,aer,adx,u0,rflect,amf,jxtra,fjact,fjtop,fjbot, &
                    fsbot,fjflx,flxd,flxd0,l1_,l2_,l_,m_,n_,atau0,qaa,a_,paa, &
                    ssa)  
    !-----------------------------------------------------------------------
    ! in:
    !	KW = wavelength bin # (1:18)
    !	KM = wavelength index for Mie properties (1:4 = 300-400-600-999 nm)
    !	WAVEL = wavelength of bin (in nm) - not now used
    !	ABX(L1_) = vertical profiles of ABSORPTION Optical Depth in each layer
    !		       includes O2 and O3 for now (BC under aerosols)
    !	AER(1:5,1:L1_) = 5 vertical profiles of Optical Depth in each layer
    !	ADX(1:5,1:L1_) = integer index of the scattering properties of each AER
    !	       1:4 are reserved for aerosols and clouds
    !	       5 is meant only for Rayleigh scattering!
    !	JXTRA(1:L1_) = number 0:J = no. of additional levels to be inserted
    !	U0  = cos (SZA)
    !	RFLECT = Lambertian albedo of surface
    ! out:
    !	FJACT(1:L_) = mean actinic flux(diff+direct) at std CTM levels(mid-layer)
    !  (new ver 5.7 diagnostics for fluxes, deposition)  fluxes 'down' are <0
    !	FJTOP = diffuse flux out top-of-atmosphere (TAU=0, above top model lauer)
    !	FJBOT = diffuse flux onto surface (<0 by definition)
    !	FSBOT = direct/solar flux onto surface  (<0 by definition)
    !	FJFLX(1:L_) = diffuse flux across top of model layer L
    !	   this connects with FJBOT = FJFLX(0) & FJTOP = FJFLX(L_+1) (not dim!!)
    !	FLXD(1:L_+1) = solar flux deposited in layer L (includes layer above CTM)
    !	   this should take into account sphericity, and is not just = mu0
    !	FLXD0 = sum of solar flux deposited in atmos
    !	   does NOT include flux on lower surface, does NOT mean absorbed!
    !-----------------------------------------------------------------------

    !	DTAUX	 Local optical depth of each CTM level
    !	TTAU	 Optical depth of air vertically above each point (to top of atm)
    !	FTAU	 Attenuation of solar beam
    !	POMEGA   Scattering phase function

    !---new Ver 5.3 code adds sub-layers (# = JXTRA(L2)) using ATAU as the
    !   factor increase from sub-layer to sub-layer

    !  fast-J Mie code for J_s, only uses 8-term expansion, 4-Gauss pts
    !  Currently allow up to A_ aerosol phase functions (at all altitudes) to
    !  be associated with optical depth AER(1:L2_) = aerosol opt.depth @ 1000 nm

    !  Pick Mie-wavelength with phase function and Qext: e.g.,
    !  01 RAYLE = Rayleigh phase
    !  02 ISOTR = isotropic
    !  03 ABSRB = fully absorbing 'soot', wavelength indep.
    !  04 X_Bkg = backgrnd stratospheric sulfate (n=1.46,log-norm:r=.09um/sigma=.6)
    !  05 X_Vol = volcanic stratospheric sulfate (n=1.46,log-norm:r=.08um/sigma=.8)
    !. . .
    !  11 W_C13 = water cloud (C1/Deirm.) (n=1.335, gamma:  r-mode=13.3um /alpha=6)
    !. . .
    !  13 Ice-H = hexagonal ice cloud (Mishchenko)
    !  14 Ice-I = irregular ice cloud (Mishchenko)

    !  Choice of the 4 aerosol/cloud indices ADX is made earlier
    !  Optical depths for the 4 (aerosol+clouds) = AER

    !-----------------------------------------------------------------------

    !LFR>       include 'parm_CTM.f'
    !LFR>       include 'parm_MIE.f'
    !LFR>       include 'cmn_JVdat.f'

    INTEGER          , INTENT(IN)    :: l1_
    INTEGER          , INTENT(IN)    :: l2_
    INTEGER          , INTENT(IN)    :: l_
    INTEGER          , INTENT(IN)    :: km
    DOUBLE PRECISION , INTENT(IN)    :: abx(l1_)
    DOUBLE PRECISION , INTENT(IN)    :: aer(5,l1_)
    INTEGER          , INTENT(IN)    :: adx(5,l1_)
    DOUBLE PRECISION , INTENT(IN)    :: u0
    DOUBLE PRECISION , INTENT(IN)    :: rflect
    DOUBLE PRECISION , INTENT(IN)    :: amf(l1_+1,l1_+1)
    INTEGER          , INTENT(IN)    :: jxtra(l2_+1)
    DOUBLE PRECISION , INTENT(OUT)   :: fjact(l_)
    DOUBLE PRECISION , INTENT(INOUT) :: fjtop ! (DMK) alterado (OUT) para (INOUT)
    DOUBLE PRECISION , INTENT(INOUT) :: fjbot ! (DMK) alterado (OUT) para (INOUT)
    DOUBLE PRECISION , INTENT(INOUT) :: fsbot ! (DMK) alterado (OUT) para (INOUT)
    DOUBLE PRECISION , INTENT(OUT)   :: fjflx(l_)
    DOUBLE PRECISION , INTENT(INOUT) :: flxd(l1_) ! (DMK) alterado (OUT) para (INOUT)
    DOUBLE PRECISION , INTENT(INOUT) :: flxd0 ! (DMK) alterado (OUT) para (INOUT)
    INTEGER          , INTENT(IN)    :: m_
    INTEGER          , INTENT(IN)    :: n_
    DOUBLE PRECISION , INTENT(IN)    :: atau0
    INTEGER          , INTENT(IN)    :: a_
    DOUBLE PRECISION , INTENT(IN)    :: qaa(5,a_)
    DOUBLE PRECISION , INTENT(IN)    :: paa(8,5,a_)
    DOUBLE PRECISION , INTENT(IN)    :: ssa(5,a_)

    INTEGER :: jndlev(l_),jnelev(l1_)
    INTEGER :: jaddlv(l2_+1),jaddto(l2_+1),l2lev(l2_+1)
    INTEGER :: i,k,l,l2,l2l,l22,lz,lzz,ndz
    INTEGER :: lz0,lz1,lzmid
    DOUBLE PRECISION  :: sumt,sumj,fjact2(l_)

    DOUBLE PRECISION :: qxmie(5),xlaer(5),ssalb(5),dtaux(l1_+1),piaer(5,l1_)  &
         ,pomegaj(2*m_,l2_+1),ttau(l2_+1),ftau(l1_+1)  &
         ,ftau2(l2_+1) &
         ,xltau,taudn,tauup  &
         ,dtauj,taubtm,tautop,fbtm,ftop,pomegab(2*m_) ,ataua,atauz
    !--- variables used in mie code-----------------------------------------
    DOUBLE PRECISION  :: fj(n_),pomega(2*m_,n_),fz(n_),ztau(n_),zrefl,zu0,zflux
    DOUBLE PRECISION  ::  fjt,fjb, fjflx0
    INTEGER :: mfit, nd

    !---Reinitialize arrays
    ztau(:)     = 0.d0
    fz(:)       = 0.d0
    pomega(:,:) = 0.d0

    !---Set up optical depth DTAUX(L), and scattering fraction PIAER(1:5,L)
    !---    where L = 1:L1_ = bulk properties of each CTM layer.
    DO l = 1,l1_
       DO i = 1,4
          IF (adx(i,l) == 0)  THEN
             qxmie(i) = 0.d0
             ssalb(i) = 0.d0
          ELSE
             !---for Mie code scale extinction at 600 nm to wavelength WAVEL (QXMIE)
             qxmie(i) = qaa(km,adx(i,l))/qaa(3,adx(i,l))
             ssalb(i) = ssa(km,adx(i,l))
             !print*,'qxmie=', real(qaa(3,adx(i,l))), real(adx(i,l)),real(qxmie(i))
          END IF
       END DO
       !---special case for Rayleigh scattering
       qxmie(5) = 1.d0
       ssalb(5) = 1.d0
       DO i = 1,5
          xlaer(i) = aer(i,l)*qxmie(i)
          !print*,'xlaer=', real(xlaer(i)), real(aer(i,l)),real(qxmie(i))
       END DO

       dtaux(l) = abx(l)
       !print*,'dtaux=',abx(l)
       DO i = 1,5
          dtaux(l) = dtaux(l) + xlaer(i)
       END DO
       !---fractional extinction for Rayleigh scattering and each aerosol type
       DO i = 1,5
          piaer(i,l) = ssalb(i)*xlaer(i)/dtaux(l)
       END DO
    END DO
    dtaux(l1_+1) = 0.d0

    !---Calculate attenuated incident beam exp(-TTAU/U0 = DTAUX * AirMassFactor)
    !---      at the edges of the CTM layers L=1:L1_
    !---  L1_ is top-edge of CTM (ie, L=38 = 2 hPa) which has TAU > 0
    !---  note that DTAUX(L1_) is optical depth in the above-CTM layer
    !---  L1_+1 = top-of-atmos which has TAU=0 by definition,
    !---     but has spherical attenuation thru atmos below if U0 < 0
    DO l = 1,l1_+1
       ftau(l) = 0.d0
    END DO
    DO l = 1,l1_+1
       IF (amf(l,l) > 0.0D0) THEN
          xltau = 0.0D0
          DO i = 1,l1_+1
             xltau = xltau + dtaux(i)*amf(i,l)
             !print*,'FTAU XLTAU=',ftau(l) ,dtaux(i),amf(i,l)
          END DO
          IF (xltau < 76.d0) THEN	! zero out flux at 1e-33
             ftau(l) = EXP(-xltau)
          END IF

       END IF
    END DO

    !---calculate direct solar flux deposited in each CTM layer: L=1:L_
    !---     use FSBOT for surface flux, cannot do layer above CTM (L_+1)
    DO l = 1,l1_
       IF (amf(l,l) > 0.d0) THEN
          flxd(l) = (ftau(l+1) - ftau(l))/amf(l,l)
       ELSE
          flxd(l) = 0.d0
       END IF
    END DO
    IF (amf(1,1) > 0.d0) THEN
       fsbot = ftau(1)/amf(1,1)
    ELSE
       fsbot = 0.d0
    END IF
    !---integrate solar flux depositied in CTM layers L=1:L_, cannot do top layer
    !---        FLXD0 = (1.d0 - FTAU(L_+1))/AMF(L_+1,L_+1) not with spherical atmos
    flxd0 = 0.d0
    IF (amf(l1_,l1_) > 0.d0) THEN
       DO l=1,l1_
          flxd0 = flxd0 + flxd(l)
       END DO
    END IF

    !---Define the total scattering phase fn for each CTM layer L=1:L_+1
    !---   from a DTAUX-wt_d mix of aerosols, cloud & Rayleigh
    !---No. of quadrature pts fixed at 4(M_), expansion of phase fn @ 8
    mfit = 2*m_
    DO l = 1,l1_
       DO i = 1,mfit
          pomegaj(i,l) = 0.d0
          DO k = 1,5
             IF (adx(k,l) > 0)  THEN
                pomegaj(i,l)=pomegaj(i,l) + piaer(k,l)*paa(i,km,adx(k,l))
             END IF
          END DO
       END DO
    END DO


    !------------------------------------------------------------------------
    !  Take optical properties on CTM layers and convert to a photolysis
    !  level grid corresponding to layer centres and boundaries. This is
    !  required so that J-values can be calculated for the centre of CTM
    !  layers; the index of these layers is kept in the JNDLEV array.
    !------------------------------------------------------------------------

    !---Now combine the CTM layer edges (1:L_+2) with the CTM mid-layer
    !---    points (1:L_) plus 1 for the mid point of added top layer.

    !---combine these edge- and mid-layer points into grid of size:
    !---  	    L2_+1 = 2*L1_+1 = 2*L_+3
    !---calculate column optical depths above each level, TTAU(1:L2_+1)
    !---      note that TTAU(L2_+1)=0 and TTAU(1)=total OD
    ttau(l2_+1) = 0.0D0
    DO l2 = l2_,1,-1
       l	       = (l2+1)/2
       dtauj      = 0.5D0 * dtaux(l)
       ttau(l2)   = ttau(l2+1) + dtauj
       !print*,'ttau=',l2,ttau(l2),dtaux(l)
    END DO

    !----solar flux incident on lower boundary & Lambertian reflect factor:
    IF (fsbot > 0.d0) THEN
       !---        FSBOT = U0*FTAU(1)
       !---        ZFLUX = U0*FTAU(1)*RFLECT/(1.d0+RFLECT)
       zflux = fsbot*rflect/(1.d0+rflect)
    ELSE
       zflux = 0.d0
    END IF

    !---calculate attenuated beam FTAU2 on the new doubled-levels L2=1:L2_+1
    !---       calculate FTAU2 at CTM mid-layers from sqrt
    !---       L2_ = 2*L1_ = 2*L_+2

    !---version 5.6 fix	(mp, 3/2007)
    !---solar flux at interp-levels: use exp(-TAU/U0) if U0>0
    IF (u0 > 0.d0) THEN
       ftau2(l2_+1) = 1.0D0
       DO l2 = l2_,1,-2
          if(l2<1) exit
          l  = l2/2
          ftau2(l2  ) = ftau2(l2+1)*EXP((ttau(l2+1)-ttau(l2))/u0)
          ftau2(l2-1) = ftau(l)
       END DO
       !---old, pre-5.6 version dropped flux in deep clouds to zero (FBTM)
    ELSE
       ftau2(l2_+1) = 1.0D0
       ftau2(l2_)   = SQRT(1.0D0*ftau(l1_))
       ftau2(l2_-1) = ftau(l1_)
       DO l2 = l2_-3,1,-2
          if(l2<1) exit
          l 	  = (l2+1)/2
          ftau2(l2)   = ftau(l)
          ftau2(l2+1) = SQRT(ftau(l+1)*ftau(l))
       END DO
    END IF


    !  Calculate scattering properties, level centres then level boundaries
    ! ***be careful of order, we are shifting the 'POMEGAJ' upward in index***
    DO l2 = l2_,2,-2
       if(l2<2) exit
       l	= l2/2
       DO i = 1,mfit
          pomegaj(i,l2) = pomegaj(i,l)
       END DO
    END DO
    !---lower boundary value is set (POMEGAJ(I,1), but set upper:
    DO i = 1,mfit
       pomegaj(i,l2_+1) = pomegaj(i,l2_)
    END DO
    !---now have POMEGAJ filled at even points from L2=3:L2_-1
    !---use inverse interpolation for correct tau-weighted values at edges
    DO l2 = 3,l2_-1,2
       taudn = ttau(l2-1)-ttau(l2)
       tauup = ttau(l2)-ttau(l2+1)
       DO i = 1,mfit
          pomegaj(i,l2) = (pomegaj(i,l2-1)*taudn +  &
               pomegaj(i,l2+1)*tauup) / (taudn+tauup)
       END DO
    END DO
    !---at this point FTAU2(1:L2_+1) and POMEAGJ(1:8, 1:L2_+1)
    !---    where FTAU2(L2_+1) = 1.0 = top-of-atmos, FTAU2(1) = surface

    !------------------------------------------------------------------------
    !  Calculate cumulative total and define levels we want J-values at.
    !  Sum upwards for levels, and then downwards for Mie code readjustments.

    !	JXTRA(L2)  Number of new levels to add between (L2) and (L2+1)
    !	      ***JXTRA(1:L2_+1) is calculated based on the aerosol+cloud OD_s
    !	JADDLV(L2)  Number of new levels actually added at each wavelength
    !	       where JADDLV = 0 when there is effectively no FTAU2
    !	JADDTO(L2)   Total number of new levels to add to and above level (L2)
    !	JNDLEV(L) = L2 index that maps on CTM mid-layer L

    !------------------------------------------------------------------------

    !---JADDLV(L2=1:L2_) = number of levels to add between TTAU2(L2) and TTAU(L2+1)
    !---    JADDLV is taken from JXTRA, which is based on visible OD.
    !---    JADDTO(L2=1:L2_+1) is the cumulative number of levels to be added
    !---note that JADDLV and JADDTO will change with wavelength and solar zenith

    !--now try to insert additional levels for thick clouds, ONLY IF FTAU2 > 1.e-8
    !-- this will cut off additional levels where the solar beam is negligible.

    !---new v5.6-----keep all wavelengths the same for now
    !	 do L2 = 1,L2_,1
    !	   if (FTAU2(L2+1) .lt. 1.d-30) then
    !	     JADDLV(L2) = 0
    !	   else
    !	     JADDLV(L2) = JXTRA(L2)
    !	   endif
    !	 enddo
    DO l2 = 1,l2_,1
       jaddlv(l2) = jxtra(l2)
    END DO
    jaddto(l2_+1) = 0
    DO l2 = l2_,1,-1
       if(l2<1) exit
       jaddto(l2) = jaddto(l2+1) + jaddlv(l2)
    END DO

    !---expanded grid now included CTM edge and mid layers plus expanded
    !---    grid to allow for finer delta-tau at tops of clouds.
    !---    DIM of new grid = L2_ + JADDTO(1) + 1

    !---L2LEV(L2) = L2-index for old level L2 in expanded J-grid (w/JADDLV)
    !	in absence of JADDLV, L2LEV(L2) = L2
    l2lev(1)  = 1
    DO l2 = 2,l2_+1
       l2lev(l2) = l2lev(l2-1) + 1 + jaddlv(l2-1)
    END DO

    !---JNDLEV(L=1:L_) = L2-index in expanded grid for CTM mid-layer L
    !---JNELEV(L=1:L_) = L2-index for top of layer L
    DO l = 1,l_
       jndlev(l) = l2lev(2*l)
       jnelev(l) = l2lev(2*l+1)
    END DO
    jnelev(l_+1) = 0  !need to set this to top-of-atmosphere
    !---------------------SET UP FOR MIE CODE-------------------------------

    !  Transpose the ascending TTAU grid to a descending ZTAU grid.
    !  Double the resolution - TTAU points become the odd points on the
    !  ZTAU grid, even points needed for asymm phase fn soln, contain 'h'.
    !  Odd point added at top of grid for unattenuated beam   (Z='inf')

    !  The following mapping holds for JADDLV=0
    !	   Surface:   TTAU(1)	 ==> ZTAU(2*L2_+1)
    !	   Top:       TTAU(L2_)  ==> ZTAU(3)
    !	   Infinity:	 0.0	 ==> ZTAU(1)
    !	   index: 2*(L2_+1-L2)+1 ==> LZ

    !  Mie scattering code only used from surface to level L2_
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    !  Insert new levels, working downwards from the top of the atmosphere
    !  to the surface (down in 'LZ', up in 'L2'). This allows ztau and pomega
    !  to be incremented linearly (in a +ve sense), and the flux fz to be
    !  attenuated top-down (avoiding problems where lower level fluxes are
    !  zero).

    !    zk	 fractional increment in level
    !    dTTAU	 change in ttau per increment	 (linear, positive)
    !    dPOMEGA   change in pomega per increment  (linear)
    !    ftaulog   change in ftau per increment	 (exponential, normally < 1)

    !------------------------------------------------------------------------

    !  Ascend through atmosphere transposing grid and adding extra points
    !  remember L2=1 is surface of CTM, but last layer (LZ) in scattering code.
    !  there are twice the number of layers in the LZ arrays (2*L2_ + 2*JADDTO + 1)
    !    because we need to insert the intermediate layers (even LZ) for the
    !    asymmetric scattering code.


    !  Transfer the L2=1:L2_+1 values (TTAU,FTAU2,POMEGAJ) onto the reverse
    !    order, expanded, doubled-level scatter grid.
    !    Note that we need to deal with the expansion by JADD levels (L2L).
    !	 These JADDLV levels are skipped and need to be interpolated later.
    !    Note that only odd LZ levels are filled,

    ndz = 2*l2_ + 2*jaddto(1) + 1

    !   Note that the successive sub-layers have the ratio in OD of ATAU
    !	 ATAUA = (ATAU - 1.d0)/ATAU	! this is the limit for L22=>inf

    DO l2 = 1,l2_+1	   ! L2 = index of CTM edge- and mid-layers
       l2l = l2lev(l2)	   ! L2L = index for L2 in expanded scale(JADD)
       lz  = ndz + 2 - 2*l2l  ! LZ = index for L2 in scatt arrays
       ztau(lz) = ttau(l2)
       fz(lz)   = ftau2(l2)
       !print*,'ztau1=',ztau(lz)
       DO i=1,mfit
          pomega(i,lz) = pomegaj(i,l2)
       END DO
    END DO

    !   Now go thru the pairs of L2 levels to see if we need JADD levels
    DO l2 = 1,l2_ 	    ! L2 = index of CTM edge- and mid-layers
       l2l = l2lev(l2)	    ! L2L = index for L2 in expanded scale(JADD)
       lz  = ndz + 2 - 2*l2l   ! LZ = index for L2 in scatt arrays
       l22 = l2lev(l2+1) - l2lev(l2) - 1	! L22 = 0 if no added levels
       IF (l22 > 0) THEN
          taubtm = ttau(l2)
          tautop = ttau(l2+1)
          fbtm   = ftau2(l2)
          ftop   = ftau2(l2+1)
          DO i = 1,mfit
             pomegab(i) = pomegaj(i,l2)
          END DO

          !---to fit L22 new layers between TAUBOT > TAUTOP, calculate new 1/ATAU factor
          !---  such that TAU(just above TAU-btm) = ATUAZ * TAUBTM < TAUBTM

          atauz = EXP(-LOG(taubtm/MAX(tautop,atau0))/FLOAT(l22+1))

          DO l = 1,l22	     ! add odd levels between L2LEV(L2) & L2LEV(L2+1)
             lzz = lz - 2*l       ! LZZ = index(odd) of added level in scatt arrays
             ztau(lzz) = taubtm * atauz

             ataua=(taubtm-ztau(lzz))/(taubtm-tautop) !fraction from TAUBTM=>TAUTOP

             !---version 5.6 fix	(mp, 3/2007)
             !---solar flux at interp-levels: use exp(TAU/U0) if U0>0, else scale by TAU
             IF (u0 > 0.d0) THEN
                fz(lzz) = ftop * EXP((tautop-ztau(lzz))/u0)
                !print*,'fz 2=',fz(lzz),n_
             ELSE
                IF (fbtm < 1.d-32) THEN
                   fz(lzz) = 0.d0
                ELSE
                   fz(lzz) = fbtm * (ftop/fbtm)**ataua
                   !print*,'fz 3=',fz(lzz),n_
                END IF
             END IF

             DO i = 1,mfit
                pomega(i,lzz) = pomegab(i) + ataua*(pomegaj(i,l2+1)-pomegab(i))
             END DO
             taubtm       = ztau(lzz)
             fbtm	     = fz(lzz)
             DO i = 1,mfit
                pomegab(i) = pomega(i,lzz)
             END DO
          END DO
       END IF
    END DO

    !   Now fill in the even points with simple interpolation in scatter arrays:
    DO lz = 2,ndz-1,2
       ztau(lz) = 0.5D0*(ztau(lz-1)+ztau(lz+1))
       fz(lz)   = SQRT(fz(lz-1)*fz(lz+1))
       DO i=1,mfit
          pomega(i,lz) = 0.5D0*(pomega(i,lz-1)+pomega(i,lz+1))
       END DO
    END DO

    nd = ndz
    zu0 = u0
    zrefl = rflect

    !---PRINT diagnostics
    !----now check integral of FZ over each CTM layer to ensure that it equals FLXD
    !	 if (KW.eq.18) then
    !	  write(6,'(A,3I6)') 'OPMIE levels: L,TAU,Fs,pomega',KW, ND,N_
    !	 do L=1,ND
    !	  write(6,'(i5,f15.5,1p,e15.5,0p,10f10.5)') L,
    !	& ZTAU(L),FZ(L),(POMEGA(I,L),I=1,3)
    !	 enddo

    !	 write(6,'(a,2i6)') '  net flux in each CTM layer',KW,L_
    !	 write(6,'(2i5,a8,1p,e12.5)')(L,JNELEV(L),' flxd ',FLXD(L),L=1,L1_)

    !	 endif

    IF(nd > n_) THEN
       WRITE(6,'(a,2i9)') ' overflow of scatter arrays:',nd,n_
       STOP
    END IF
    ! print*,'XX1',fj,fjt
    ! print*,'XX2',fjb,pomega
    ! print*,'XX1',ztau
    ! print*,'XX2',zflux
    ! print*,'XX4',zrefl
    ! print*,'XX5',zu0
    ! print*,'XX6',fz
    ! print*,'XX7',mfit
    ! print*,'XX8',nd
    !-----------------------------------------------------------------------
    CALL miesct(fj,fjt,fjb,pomega,fz,ztau,zflux,zrefl,zu0,mfit,nd,m_,n_)
    !  print*,'XX2',fj(1:l_), fz(1:l_)
    !-----------------------------------------------------------------------

    !---Move mean intensity from scatter array FJ(LZ=1:ND)
    !---  	    to CTM mid-level array FJACT(L=1:L_)

    !---mean intensity:  4*<I> + solar at mid-layer
    DO l = 1,l_
       l2l = jndlev(l)
       lz  = nd+2 - 2*l2l
       fjact(l) = 4.d0*fj(lz) + fz(lz)
       !print*,'opmie:',l,fjact(l)
    END DO

    !---mean diffuse flux:  4<I*mu> (not solar) at top of layer L
    !---      average (tau-wtd) the h's just above and below the L-edge
    DO l = 1,l_
       l2l = jnelev(l)
       lz  = nd+2 - 2*l2l
       !---       FJFLX(L) = 2.0d0*(FJ(LZ-1) + FJ(LZ+1))
       fjflx0 = (ztau(lz+1)-ztau(lz))/(ztau(lz+1)-ztau(lz-1))
       fjflx(l) = 4.d0*(fj(lz-1)*fjflx0 + fj(lz+1)*(1.d0-fjflx0))
    END DO

    !---NB if one needs the mean intensity throughout layer L (instead of mid-pt)
    !---   then average (tau-weighted) the odd-points from: NELEV(L-1) to NELEV(L)
    !---NB This is NOT now used.  Errors in cloudy layers are 0-10%  (low sun)
    !---   the results are stored locally in FJACT2 (actinic)
    lz1 = nd
    DO l = 1,l_
       lzmid = nd+2-2*jndlev(l)
       lz0 = nd+2-2*jnelev(l)
       sumt = 0.d0
       sumj = 0.d0
       DO l2 = lz0,lz1-2,2
          sumt = sumt + ztau(l2+2)-ztau(l2)
          sumj = sumj + (ztau(l2+2)-ztau(l2))*(fj(l2)+fj(l2+2))*0.5
       END DO
       fjact2(l) = 4.d0*sumj/sumt + fz(lzmid)
       !-----LZ are indices to the 1:ND array     LZ1 > LZMID > LZ0
       !	  write(6,'(4i5,3f10.3)') L, LZ0,LZMID,LZ1,
       !	&      ZTAU(LZMID),FJACT(L),FJACT2(L)
       lz1 = lz0
    END DO

    !---diffuse fluxes reflected at top, incident at bottom
    fjtop = fjt
    fjbot = fjb

  END SUBROUTINE opmie

  !<<<<<<<<<<<<<<<<<<<<<<<end core fast-J subroutines<<<<<<<<<<<<<<<<<<<<<

  SUBROUTINE blkslv(fj,pomega,fz,ztau,zflux,zrefl,wt,emu,pm,pm0,fjtop,fjbot,m,n,mfit,nd,m_,n_)
    !-----------------------------------------------------------------------
    !  Sets up and solves the block tri-diagonal system:
    !		  A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
    !  Needs subroutine to generate the matrix (GEN) solve it (DECBT+SOLBT)
    !	      can also use non-LLNL solver (TRIDAG)
    !-----------------------------------------------------------------------

    !LFR>       include 'parm_MIE.f'
    !--- expect parameters M_, N_ in parm_MIE.f------------------------------  

    INTEGER          , INTENT(IN)    :: m
    INTEGER          , INTENT(IN)    :: n
    INTEGER          , INTENT(IN)    :: mfit
    INTEGER          , INTENT(IN)    :: nd
    DOUBLE PRECISION , INTENT(IN)    :: wt(m_)
    DOUBLE PRECISION , INTENT(IN)    :: zflux
    DOUBLE PRECISION , INTENT(IN)    :: pomega(2*m_,n_)
    DOUBLE PRECISION , INTENT(IN)    :: fz(n_)
    DOUBLE PRECISION , INTENT(IN)    :: ztau(n_)
    DOUBLE PRECISION , INTENT(IN)    :: zrefl
    DOUBLE PRECISION , INTENT(IN)    :: emu(m_)
    DOUBLE PRECISION , INTENT(IN)    :: pm(m_,2*m_)
    DOUBLE PRECISION , INTENT(IN)    :: pm0(2*m_)
    DOUBLE PRECISION , INTENT(INOUT) :: fjtop ! (DMK) alterado (OUT) para (INOUT)
    DOUBLE PRECISION , INTENT(INOUT) :: fjbot ! (DMK) alterado (OUT) para (INOUT)
    DOUBLE PRECISION , INTENT(INOUT) :: fj(n_) ! (DMK) alterado (OUT) para (INOUT)
    INTEGER          , INTENT(IN)    :: m_
    INTEGER          , INTENT(IN)    :: n_


    DOUBLE PRECISION, DIMENSION(m_,n_)    :: hh, rr
    DOUBLE PRECISION, DIMENSION(m_,m_,n_) :: aa, bb, cc, dd
    INTEGER :: i, id
    DOUBLE PRECISION  :: fiplus
    !-----------------------------------------------------------------------

    CALL gen (pomega,fz,ztau,zflux,zrefl,wt,emu,pm,pm0,aa,bb,cc,hh, m,n,mfit,nd,m_,n_)

    !----UCI solver (.le. 5.3c):
    CALL tridag (aa,bb,cc,hh, dd,rr, n,nd,m_,n_)

    !----------MEAN J & H
    DO id = 1,nd,2
       fj(id) = 0.0D0
       DO i = 1,n
          fj(id) = fj(id) + rr(i,id)*wt(i)
       END DO
    END DO
    DO id = 2,nd,2
       fj(id) = 0.0D0
       DO i = 1,n
          fj(id) = fj(id) + rr(i,id)*wt(i)*emu(i)
       END DO
    END DO

    fjtop = 0.0D0
    fjbot = 0.0D0
    DO i = 1,n
       fjtop = fjtop + rr(i,1)*wt(i)*emu(i)
       fjbot = fjbot + rr(i,nd)*wt(i)*emu(i)
    END DO
    !---FJTOP = scaled diffuse flux out top-of-atmosphere (limit = mu0)
    fjtop = 4.d0*fjtop
    !---FJBOT = scaled diffuse flux onto surface:
    !---         ZFLUX = reflect/(1 + reflect) * mu0 * Fsolar(lower boundary)
    fiplus = 4.d0*zrefl*fjbot/(1.0D0 + zrefl) + zflux
    fjbot = 4.d0*fjbot - fiplus

  END SUBROUTINE blkslv


  SUBROUTINE gen(pomega,fz,ztau,zflux,zrefl,wt,emu,pm,pm0,aa,bb,cc,hh, m,n,mfit,nd,m_,n_)
    !-----------------------------------------------------------------------
    !  Generates coefficient matrices for the block tri-diagonal system:
    !		  A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
    !-----------------------------------------------------------------------

    !LFR>       include 'parm_MIE.f'
    !--- expect parameters M_, N_ in parm_MIE.f------------------------------

    INTEGER          , INTENT(IN)    :: m
    INTEGER          , INTENT(IN)    :: n
    INTEGER          , INTENT(IN)    :: mfit
    INTEGER          , INTENT(IN)    :: nd
    DOUBLE PRECISION , INTENT(IN)    :: wt(m_)
    DOUBLE PRECISION , INTENT(IN)    :: zflux
    DOUBLE PRECISION , INTENT(IN)    :: pomega(2*m_,n_)
    DOUBLE PRECISION , INTENT(IN)    :: fz(n_)
    DOUBLE PRECISION , INTENT(IN)    :: ztau(n_)
    DOUBLE PRECISION , INTENT(IN)    :: zrefl
    DOUBLE PRECISION , INTENT(IN)    :: emu(m_)
    DOUBLE PRECISION , INTENT(IN)    :: pm(m_,2*m_)
    DOUBLE PRECISION , INTENT(IN)    :: pm0(2*m_)
    DOUBLE PRECISION , INTENT(INOUT) :: aa(m_,m_,n_) ! (DMK) alterado (OUT) para (INOUT)
    DOUBLE PRECISION , INTENT(INOUT) :: bb(m_,m_,n_) ! (DMK) alterado (OUT) para (INOUT)
    DOUBLE PRECISION , INTENT(INOUT) :: cc(m_,m_,n_) ! (DMK) alterado (OUT) para (INOUT)
    DOUBLE PRECISION , INTENT(INOUT) :: hh(m_,n_) ! (DMK) alterado (OUT) para (INOUT)
    INTEGER          , INTENT(IN)    :: m_
    INTEGER          , INTENT(IN)    :: n_

    ! LOCAL:
    INTEGER :: id, id0, id1, im, i, j, k, mstart
    DOUBLE PRECISION :: sum0, sum1, sum2, sum3
    DOUBLE PRECISION :: deltau, d1, d2, surfac

    DOUBLE PRECISION :: s(m_,m_), w(m_,m_), w1(m_,m_), u1(m_,m_), u(m_,m_)
    DOUBLE PRECISION :: b1(m_,m_),c(m_),c1(m_),h1(m_)
    !---------------------------------------------

    DO id = 1,nd
       DO i = 1,n
          hh(i,id) = 0.0D0
       END DO
       DO i = 1,n
          DO j = 1,n
             aa(i,j,id) = 0.d0
             bb(i,j,id) = 0.d0
             cc(i,j,id) = 0.d0
          END DO
       END DO
    END DO

    !-------------upper boundary (top of atmosphere) (1): 2nd-order b.c.
    id0 = 1
    id1 = 2

    DO i = 1,n
       sum0 = 0.0D0
       sum1 = 0.0D0
       sum2 = 0.0D0
       sum3 = 0.0D0
       DO im = m,mfit,2
          sum0 = sum0 + pomega(im,id0)*pm(i,im)*pm0(im)
          sum2 = sum2 + pomega(im,id1)*pm(i,im)*pm0(im)
       END DO
       DO im = m+1,mfit,2
          sum1 = sum1 + pomega(im,id0)*pm(i,im)*pm0(im)
          sum3 = sum3 + pomega(im,id1)*pm(i,im)*pm0(im)
       END DO
       h1(i) = 0.5D0*(sum0*fz(id0) + sum2*fz(id1))
       c(i) = 0.5D0*(sum1*fz(id0) + sum3*fz(id1))
       DO j = 1,i
          sum0 = 0.0D0
          sum1 = 0.0D0
          sum2 = 0.0D0
          sum3 = 0.0D0
          DO im = m,mfit,2
             sum0 = sum0 + pomega(im,id0)*pm(i,im)*pm(j,im)
             sum2 = sum2 + pomega(im,id1)*pm(i,im)*pm(j,im)
          END DO
          DO im = m+1,mfit,2
             sum1 = sum1 + pomega(im,id0)*pm(i,im)*pm(j,im)
             sum3 = sum3 + pomega(im,id1)*pm(i,im)*pm(j,im)
          END DO
          s(i,j) = - sum2*wt(j)
          s(j,i) = - sum2*wt(i)
          w(i,j) = - sum1*wt(j)
          w(j,i) = - sum1*wt(i)
          u(i,j) = - sum3*wt(j)
          u(j,i) = - sum3*wt(i)
          sum0 = 0.5D0*(sum0 + sum2)
          b1(i,j) = - sum0*wt(j)
          b1(j,i) = - sum0*wt(i)
       END DO
       s(i,i) = s(i,i) + 1.0D0
       w(i,i) = w(i,i) + 1.0D0
       u(i,i) = u(i,i) + 1.0D0
       b1(i,i) = b1(i,i) + 1.0D0
    END DO

    DO i = 1,n
       c1(i) = 0.0D0
       DO j = 1,n
          c1(i) = c1(i) + s(i,j)*c(j)/emu(j)
       END DO
    END DO

    DO i = 1,n
       DO j = 1,n
          w1(j,i) = 0.0D0
          u1(j,i) = 0.0D0
          DO k = 1,n
             w1(j,i) = w1(j,i) + s(j,k)*w(k,i)/emu(k)
             u1(j,i) = u1(j,i) + s(j,k)*u(k,i)/emu(k)
          END DO
       END DO
    END DO

    deltau = ztau(2) - ztau(1)
    d2 = 0.25D0*deltau

    DO i = 1,n
       d1 = emu(i)/deltau
       DO j = 1,n
          bb(i,j,1) = b1(i,j) + d2*w1(i,j)
          cc(i,j,1) = d2*u1(i,j)
       END DO
       bb(i,i,1) = bb(i,i,1) + d1
       cc(i,i,1) = cc(i,i,1) - d1
       hh(i,1) = h1(i) + 2.d0*d2*c1(i)
       !cccc     HH(I,1) = HH(I,1) + D1*SISOTP    ! w/isotropic flux at top
    END DO


    !-------------lower (surface) boundary (ND): 2nd-order b.c.
    id0 = nd
    id1 = nd-1

    DO i = 1,n
       sum0 = 0.0D0
       sum1 = 0.0D0
       sum2 = 0.0D0
       sum3 = 0.0D0
       DO im = m,mfit,2
          sum0 = sum0 + pomega(im,id0)*pm(i,im)*pm0(im)
          sum2 = sum2 + pomega(im,id1)*pm(i,im)*pm0(im)
       END DO
       DO im = m+1,mfit,2
          sum1 = sum1 + pomega(im,id0)*pm(i,im)*pm0(im)
          sum3 = sum3 + pomega(im,id1)*pm(i,im)*pm0(im)
       END DO
       h1(i) = 0.5D0*(sum0*fz(id0) + sum2*fz(id1))
       c(i) = 0.5D0*(sum1*fz(id0) + sum3*fz(id1))
       DO j = 1,i
          sum0 = 0.0D0
          sum1 = 0.0D0
          sum2 = 0.0D0
          sum3 = 0.0D0
          DO im = m,mfit,2
             sum0 = sum0 + pomega(im,id0)*pm(i,im)*pm(j,im)
             sum2 = sum2 + pomega(im,id1)*pm(i,im)*pm(j,im)
          END DO
          DO im = m+1,mfit,2
             sum1 = sum1 + pomega(im,id0)*pm(i,im)*pm(j,im)
             sum3 = sum3 + pomega(im,id1)*pm(i,im)*pm(j,im)
          END DO
          s(i,j) = - sum2*wt(j)
          s(j,i) = - sum2*wt(i)
          w(i,j) = - sum1*wt(j)
          w(j,i) = - sum1*wt(i)
          u(i,j) = - sum3*wt(j)
          u(j,i) = - sum3*wt(i)
          sum0 = 0.5D0*(sum0 + sum2)
          b1(i,j) = - sum0*wt(j)
          b1(j,i) = - sum0*wt(i)
       END DO
       s(i,i) = s(i,i) + 1.0D0
       w(i,i) = w(i,i) + 1.0D0
       u(i,i) = u(i,i) + 1.0D0
       b1(i,i) = b1(i,i) + 1.0D0
    END DO

    DO i = 1,n
       c1(i) = 0.0D0
       DO j = 1,n
          c1(i) = c1(i) + s(i,j)*c(j)/emu(j)
       END DO
    END DO

    DO i = 1,n
       DO j = 1,n
          w1(j,i) = 0.0D0
          u1(j,i) = 0.0D0
          DO k = 1,n
             w1(j,i) = w1(j,i) + s(j,k)*w(k,i)/emu(k)
             u1(j,i) = u1(j,i) + s(j,k)*u(k,i)/emu(k)
          END DO
       END DO
    END DO

    deltau = ztau(nd) - ztau(nd-1)
    d2 = 0.25D0*deltau
    surfac = 4.0D0*zrefl/(1.0D0 + zrefl)

    DO i = 1,n
       d1 = emu(i)/deltau
       sum0 = 0.0D0
       DO j = 1,n
          sum0 = sum0 + w1(i,j)
       END DO
       sum0 = d1 + d2*sum0
       sum1 = surfac*sum0
       DO j = 1,n
          bb(i,j,nd) = b1(i,j) + d2*w1(i,j) - sum1*emu(j)*wt(j)
          aa(i,j,nd) = -d2*u1(i,j)
       END DO
       bb(i,i,nd) = bb(i,i,nd) + d1
       aa(i,i,nd) = aa(i,i,nd) + d1
       hh(i,nd) = h1(i) - 2.0D0*d2*c1(i) + sum0*zflux
    END DO


    !------------intermediate points (2:ND-1), can be even or odd, A & C diagonal
    DO id = 2,nd-1
       deltau = ztau(id+1) - ztau(id-1)
       mstart = m + MOD(id+1,2)

       !---right-hand side:  HH
       DO i = 1,n
          DO im = mstart,mfit,2
             hh(i,id) = hh(i,id) + pomega(im,id)*pm(i,im)*pm0(im)*fz(id)
          END DO
       END DO

       !---off-diagonal blocks: AA & CC
       DO i = 1,n
          aa(i,i,id) =  emu(i)/deltau
          cc(i,i,id) = -emu(i)/deltau
       END DO

       !---diagonal blocks: BB
       DO i = 1,n
          DO j = 1,i
             sum0 = 0.0D0
             DO im = mstart,mfit,2
                sum0 = sum0 + pomega(im,id)*pm(i,im)*pm(j,im)
             END DO
             bb(i,j,id) =  -sum0*wt(j)
             bb(j,i,id) =  -sum0*wt(i)
          END DO
          bb(i,i,id) = bb(i,i,id) + 1.0D0
       END DO

    END DO    ! ID=2,ND-1
    !------------

  END SUBROUTINE gen

  !<<<<<<<<<<<<<<<<<<<<<<<begin core scattering subroutines<<<<<<<<<<<<<<<
  !-----------------------------------------------------------------------
  SUBROUTINE miesct(fj,fjt,fjb, pomega,fz,ztau,zflux,zrefl,zu0,mfit,nd,m_,n_)
    !-----------------------------------------------------------------------
    !   This is an adaption of the Prather radiative transfer code, (mjp, 10/95)
    !	Prather, 1974, Astrophys. J. 192, 787-792.
    !	    Sol_n of inhomogeneous Rayleigh scattering atmosphere.
    !	    (original Rayleigh w/ polarization)
    !	Cochran and Trafton, 1978, Ap.J., 219, 756-762.
    !	    Raman scattering in the atmospheres of the major planets.
    !	    (first use of anisotropic code)
    !	Jacob, Gottlieb and Prather, 1989, J.Geophys.Res., 94, 12975-13002.
    !	    Chemistry of a polluted cloudy boundary layer,
    !	    (documentation of extension to anisotropic scattering)

    !    takes atmospheric structure and source terms from std J-code
    !    ALSO limited to 4 Gauss points, only calculates mean field!

    !   mean rad. field ONLY (M=1)
    !-----------------------------------------------------------------------

    !LFR>       include 'parm_MIE.f'
    !--- expect parameters M_, N_ in parm_MIE.f------------------------------

    INTEGER          , INTENT(IN)    :: mfit
    INTEGER          , INTENT(IN)    :: nd
    DOUBLE PRECISION , INTENT(INOUT) :: fj(n_) ! (DMK) alterado (OUT) para (INOUT)
    DOUBLE PRECISION , INTENT(INOUT) :: fjt ! (DMK) alterado (OUT) para (INOUT)
    DOUBLE PRECISION , INTENT(INOUT) :: fjb ! (DMK) alterado (OUT) para (INOUT)
    DOUBLE PRECISION , INTENT(IN)    :: pomega(2*m_,n_)
    DOUBLE PRECISION , INTENT(IN)    :: fz(n_)
    DOUBLE PRECISION , INTENT(IN)    :: ztau(n_)
    DOUBLE PRECISION , INTENT(IN)    :: zflux
    DOUBLE PRECISION , INTENT(IN)    :: zrefl
    DOUBLE PRECISION , INTENT(IN)    :: zu0
    INTEGER          , INTENT(IN)    :: m_
    INTEGER          , INTENT(IN)    :: n_

    DOUBLE PRECISION :: wt(m_),emu(m_),pm(m_,2*m_),pm0(2*m_),cmeq1
    INTEGER :: i, im, m, n
    !-----------------------------------------------------------------------
    !---fix scattering to 4 Gauss pts = 8-stream
    CALL gaussp (n,m_,emu,wt)
    !---calc in OPMIE:  ZFLUX = (ZU0*FZ(ND)*ZREFL)/(1.0d0+ZREFL)
    m = 1
    DO i = 1,n
       CALL legnd0 (emu(i),pm0,mfit)
       DO im = m,mfit
          pm(i,im) = pm0(im)
       END DO
    END DO

    cmeq1 = 0.25D0
    CALL legnd0 (-zu0,pm0,mfit)
    DO im=m,mfit
       pm0(im) = cmeq1*pm0(im)
    END DO

    CALL blkslv(fj,pomega,fz,ztau,zflux,zrefl,wt,emu,pm,pm0  &
         ,fjt,fjb,  m,n,mfit,nd,m_,n_)

    !	 do ID=1,ND,2
    !	   FJ(ID) = 4.0d0*FJ(ID) + FZ(ID)
    !	 enddo


  END SUBROUTINE miesct



  SUBROUTINE tridag (aa,bb,cc,hh, dd,rr, n,nd,m_,n_)
    !-----------------------------------------------------------------------
    !  Solves the block tri-diagonal system for R(N,ND):
    !		  A(I)*R(I-1) + B(I)*R(I) + C(I)*R(I+1) = H(I)
    !-----------------------------------------------------------------------

    !LFR>       include 'parm_MIE.f'
    !--- expect parameters M_, N_ in parm_MIE.f------------------------------

    INTEGER          , INTENT(IN)  :: n
    INTEGER          , INTENT(IN)  :: nd
    DOUBLE PRECISION , INTENT(IN)  :: aa(m_,m_,n_)
    DOUBLE PRECISION , INTENT(IN)  :: bb(m_,m_,n_)
    DOUBLE PRECISION , INTENT(IN)  :: cc(m_,m_,n_)
    DOUBLE PRECISION , INTENT(IN)  :: hh(m_,n_)
    DOUBLE PRECISION , INTENT(OUT) :: dd(m_,m_,n_)
    DOUBLE PRECISION , INTENT(OUT) :: rr(m_,n_)
    INTEGER          , INTENT(IN)  :: m_
    INTEGER          , INTENT(IN)  :: n_


    DOUBLE PRECISION, DIMENSION(m_,m_)    :: b
    DOUBLE PRECISION, DIMENSION(m_)       :: r
    INTEGER :: i, j, k, id

    rr(:,:) = 0.d0
    dd(:,:,:) = 0.d0

    id = 1
    b(:,:) = bb(:,:,id)
    CALL matin4 (b)
    DO i = 1,n
       DO j = 1,n
          DO k = 1,n
             dd(i,j,id) = dd(i,j,id) - b(i,k)*cc(k,j,id)
          END DO
          rr(i,id) = rr(i,id) + b(i,j)*hh(j,id)
       END DO
    END DO

    DO id = 2,nd-1
       r(:) = hh(:,id)
       b(:,:) = bb(:,:,id)
       DO j = 1,n
          DO i = 1,n
             DO k = 1,n
                b(i,j) = b(i,j) + aa(i,k,id)*dd(k,j,id-1)
             END DO
             r(j) = r(j) - aa(j,i,id)*rr(i,id-1)
          END DO
       END DO
       CALL matin4 (b)

       DO i = 1,n
          DO j = 1,n
             DO k = 1,n
                dd(i,j,id) = dd(i,j,id) - b(i,k)*cc(k,j,id)
             END DO
             rr(i,id) = rr(i,id) + b(i,j)*r(j)
          END DO
       END DO
    END DO

    id = nd
    r(:) = hh(:,id)
    b(:,:) = bb(:,:,id)
    DO j = 1,n
       DO i = 1,n
          DO k = 1,n
             b(i,j) = b(i,j) + aa(i,k,id)*dd(k,j,id-1)
          END DO
          r(j) = r(j) - aa(j,i,id)*rr(i,id-1)
       END DO
    END DO
    CALL matin4 (b)
    DO i = 1,n
       DO j = 1,n
          rr(i,id) = rr(i,id) + b(i,j)*r(j)
       END DO
    END DO

    DO id = nd-1,1,-1
       DO i = 1,n
          DO j = 1,n
             rr(i,id) = rr(i,id) + dd(i,j,id)*rr(j,id+1)
          END DO
       END DO
    END DO


  END SUBROUTINE tridag


  SUBROUTINE legnd0 (x,pl,n)
    !---Calculates ORDINARY Legendre fns of X (real)
    !---   from P[0] = PL(1) = 1,  P[1] = X, .... P[N-1] = PL(N)

    INTEGER          , INTENT(IN)  :: n
    DOUBLE PRECISION , INTENT(IN)  :: x
    DOUBLE PRECISION , INTENT(OUT) :: pl(n)

    INTEGER :: i
    DOUBLE PRECISION  :: den
    !---Always does PL(2) = P[1]
    pl(1) = 1.d0
    pl(2) = x
    DO i = 3,n
       den = (i-1)
       pl(i) = pl(i-1)*x*(2.d0-1.0/den) - pl(i-2)*(1.d0-1.d0/den)
    END DO

  END SUBROUTINE legnd0


  SUBROUTINE matin4 (a)
    !-----------------------------------------------------------------------
    !  invert 4x4 matrix A(4,4) in place with L-U decomposition (mjp, old...)
    !-----------------------------------------------------------------------

    DOUBLE PRECISION, INTENT(INOUT)  :: a(4,4)

    !---SETUP L AND U
    a(2,1) = a(2,1)/a(1,1)
    a(2,2) = a(2,2)-a(2,1)*a(1,2)
    a(2,3) = a(2,3)-a(2,1)*a(1,3)
    a(2,4) = a(2,4)-a(2,1)*a(1,4)
    a(3,1) = a(3,1)/a(1,1)
    a(3,2) = (a(3,2)-a(3,1)*a(1,2))/a(2,2)
    a(3,3) = a(3,3)-a(3,1)*a(1,3)-a(3,2)*a(2,3)
    a(3,4) = a(3,4)-a(3,1)*a(1,4)-a(3,2)*a(2,4)
    a(4,1) = a(4,1)/a(1,1)
    a(4,2) = (a(4,2)-a(4,1)*a(1,2))/a(2,2)
    a(4,3) = (a(4,3)-a(4,1)*a(1,3)-a(4,2)*a(2,3))/a(3,3)
    a(4,4) = a(4,4)-a(4,1)*a(1,4)-a(4,2)*a(2,4)-a(4,3)*a(3,4)
    !---INVERT L
    a(4,3) = -a(4,3)
    a(4,2) = -a(4,2)-a(4,3)*a(3,2)
    a(4,1) = -a(4,1)-a(4,2)*a(2,1)-a(4,3)*a(3,1)
    a(3,2) = -a(3,2)
    a(3,1) = -a(3,1)-a(3,2)*a(2,1)
    a(2,1) = -a(2,1)
    !---INVERT U
    a(4,4) = 1.d0/a(4,4)
    a(3,4) = -a(3,4)*a(4,4)/a(3,3)
    a(3,3) = 1.d0/a(3,3)
    a(2,4) = -(a(2,3)*a(3,4)+a(2,4)*a(4,4))/a(2,2)
    a(2,3) = -a(2,3)*a(3,3)/a(2,2)
    a(2,2) = 1.d0/a(2,2)
    a(1,4) = -(a(1,2)*a(2,4)+a(1,3)*a(3,4)+a(1,4)*a(4,4))/a(1,1)
    a(1,3) = -(a(1,2)*a(2,3)+a(1,3)*a(3,3))/a(1,1)
    a(1,2) = -a(1,2)*a(2,2)/a(1,1)
    a(1,1) = 1.d0/a(1,1)
    !---MULTIPLY (U-INVERSE)*(L-INVERSE)
    a(1,1) = a(1,1)+a(1,2)*a(2,1)+a(1,3)*a(3,1)+a(1,4)*a(4,1)
    a(1,2) = a(1,2)+a(1,3)*a(3,2)+a(1,4)*a(4,2)
    a(1,3) = a(1,3)+a(1,4)*a(4,3)
    a(2,1) = a(2,2)*a(2,1)+a(2,3)*a(3,1)+a(2,4)*a(4,1)
    a(2,2) = a(2,2)+a(2,3)*a(3,2)+a(2,4)*a(4,2)
    a(2,3) = a(2,3)+a(2,4)*a(4,3)
    a(3,1) = a(3,3)*a(3,1)+a(3,4)*a(4,1)
    a(3,2) = a(3,3)*a(3,2)+a(3,4)*a(4,2)
    a(3,3) = a(3,3)+a(3,4)*a(4,3)
    a(4,1) = a(4,4)*a(4,1)
    a(4,2) = a(4,4)*a(4,2)
    a(4,3) = a(4,4)*a(4,3)

  END SUBROUTINE matin4


  SUBROUTINE gaussp (n,m,xpt,xwt)
    !-----------------------------------------------------------------------
    !  Loads in pre-set Gauss points for 4 angles from 0 to +1 in cos(theta)=mu
    !-----------------------------------------------------------------------

    INTEGER          , INTENT(OUT)   :: n
    INTEGER          , INTENT(IN)    :: m
    DOUBLE PRECISION , INTENT(INOUT) :: xpt(m)
    DOUBLE PRECISION , INTENT(INOUT) :: xwt(m)

    DOUBLE PRECISION  :: gpt4(4),gwt4(4)
    INTEGER :: i
    DATA gpt4/.06943184420297D0,.33000947820757D0,.66999052179243D0,  &
         .93056815579703D0/
    DATA gwt4/.17392742256873D0,.32607257743127D0,.32607257743127D0,  &
         .17392742256873D0/

    n = 4
    DO i = 1,n
       xpt(i) = gpt4(i)
       xwt(i) = gwt4(i)
    END DO

  END SUBROUTINE gaussp


  !SUBROUTINE efold (f0, f1, n, f)
  !  !-----------------------------------------------------------------------
  !  !	 ***not used in fast-J, part of original scattering code
  !  !---Speciality subroutine for calculating consistent exp(-tau/mu0)
  !  !---  values on the tau grid so that photons are conserved.
  !---  ***only works for plane-parallel, NOT psuedo-spherical atmos
  !  
  !  !---  calculate the e-fold between two boundaries, given the value
  !  !---     at both boundaries F0(x=0) = top, F1(x=1) = bottom.
  !  !---  presume that F(x) proportional to exp[-A*x] for x=0 to x=1
  !  !--- 	 d2F/dx2 = A*A*F  and thus expect F1 = F0 * exp[-A]
  !  !--- 	  alternatively, could define A = ln[F0/F1]
  !  !---  let X = A*x, d2F/dX2 = F
  !  !---  assume equal spacing (not necessary, but makes this easier)
  !  !---      with N-1 intermediate points (and N layers of thickness dX = A/N)
  !  !---
  !  !---  2nd-order finite difference:  (F(i-1) - 2F(i) + F(i+1)) / dX*dX = F(i)
  !---      let D = 1 / dX*dX:
  !  
  !  !  1  |   1        0        0	0	 0	  0   |    | F0 |
  !  !	 |						      |    | 0  |
  !  !  2  |  -D      2D+1      -D	0	 0	  0   |    | 0  |
  !  !	 |						      |    | 0  |
  !  !  3  |   0       -D      2D+1      -D	 0	  0   |    | 0  |
  !  !	 |						      |    | 0  |
  !  !	 |   0        0       -D      2D+1	-D	  0   |    | 0  |
  !  !	 |						      |    | 0  |
  !  !  N  |   0        0        0       -D      2D+1	 -D   |    | 0  |
  !  !	 |						      |    | 0  |
  ! N+1 |   0        0        0	0	 0	  1   |    | F1 |
  !
  !  !-----------------------------------------------------------------------
  !  !  Advantage of scheme over simple attenuation factor: conserves total
  !  !  number of photons - very useful when using scheme for heating rates.
  !  !  Disadvantage: although reproduces e-folds very well for small flux
  !  !  differences, starts to drift off when many orders of magnitude are
  !  !  involved.
  !  !-----------------------------------------------------------------------
  !  IMPLICIT NONE
  !    
  !  INTEGER	   , intent(in)    :: n
  !   DOUBLE PRECISION, intent(in)    :: f0
  !  DOUBLE PRECISION, intent(in)    :: f1
  !  DOUBLE PRECISION, intent(out)   :: f(101)  
  !  INTEGER :: i
  !  DOUBLE PRECISION :: a,dx,d,dsq,ddp1, b(101),r(101)
  !  
  !  IF (f0 == 0.d0) THEN
  !    DO i = 1,n
  !      f(i)=0.d0
  !    END DO
  !    RETURN
  !  ELSE IF (f1 == 0.d0) THEN
  !     a = LOG(f0/1.d-250)
  !  ELSE
  !    a = LOG(f0/f1)
  !  END IF
  !  dx = FLOAT(n)/a
  !  d = dx*dx
  !  dsq = d*d
  !  ddp1 = d+d+1.d0
  !  b(2) = ddp1
  !  r(2) = +d*f0
  !  DO i = 3,n
  !    b(i) = ddp1 - dsq/b(i-1)
  !     r(i) = +d*r(i-1)/b(i-1)
  !  END DO
  !  f(n+1) = f1
  !  DO i = n,2,-1
  !    f(i) = (r(i) + d*f(i+1))/b(i)
  !  END DO
  !  f(1) = f0
  !
  !END SUBROUTINE efold


  SUBROUTINE rd_js(nj1,namfil,jvn_,x_,njval,titlej,jlabel,jfacta)
    !-----------------------------------------------------------------------
    !  Reread the ratj.dat file to map photolysis rate to reaction
    !  Read in quantum yield 'jfacta' and fastj2 label 'jlabel'
    !-----------------------------------------------------------------------

    !	jfacta    Quantum yield (or multiplication factor) for photolysis
    !	jlabel    Reference label identifying appropriate J-value to use
    !	ipr	  Photolysis reaction counter - should total 'JVN_'

    !-----------------------------------------------------------------------

    !LFR>       include 'parm_CTM.f'
    !LFR>       include 'cmn_metdat.f'
    !LFR>       include 'cmn_JVdat.f'

    INTEGER          , INTENT(IN)    :: nj1
    INTEGER          , INTENT(IN)    :: jvn_
    CHARACTER(LEN=*) , INTENT(IN)    :: namfil  

    INTEGER          , INTENT(IN)    :: x_
    INTEGER          , INTENT(IN)    :: njval
    CHARACTER(LEN=7) , INTENT(IN)    :: titlej(x_)
    CHARACTER(LEN=7) , INTENT(INOUT) :: jlabel(jvn_)
    DOUBLE PRECISION , INTENT(INOUT) :: jfacta(jvn_)



    INTEGER :: ipr, j, k
    CHARACTER (LEN=120) :: cline

    INTEGER :: jind(jvn_) ! (DMK) scratch
    INTEGER :: nratj      ! (DMK) scratch


    ! Reread the ratj.dat file to map photolysis rate to reaction
    !			Read in quantum yield jfacta and fastj2 label jlabel
    ipr = 0
    OPEN (nj1,FILE=namfil,STATUS='old',FORM='formatted')
10  READ (nj1,'(A)',ERR=20)  cline
    IF (ipr == jvn_) GO TO 20

    IF (cline(2:5) == '9999') THEN
       GO TO 20
    ELSE IF (cline(1:1) == '#') THEN
       GO TO 10
    ELSE IF (cline(5:5) == '$') THEN
       GO TO 10
    ELSE
       ipr = ipr+1
       READ (cline(79:83),'(F5.1)') jfacta(ipr)
       READ (cline(86:92),'(A7)')   jlabel(ipr)
       jfacta(ipr) = jfacta(ipr)/100.d0
       GO TO 10
    END IF
20  CLOSE(nj1)

    nratj = ipr
    ! print*,'nrat=', ipr, jvn_!,cline(2:5),cline(1:1),cline(5:5)
    ! stop 31

    !-----------------------------------------------------------------------
    !  compare Xsections titles with J-values listed in chem code (jratd.dat)
    !  map the J-values needed for chemistry (ratj.dat) onto the fast-JX rates
    !  >>>>>>>>>>>>>>>>current code revised to JPL-02 ver 8.5 (5/05)<<<<<<<<<
    !	     >>>this must now follow the read in of Xsects, etc<<<
    !-----------------------------------------------------------------------

    !---Zero / Set index arrays that map Jvalue(j) onto rates
    DO j = 1,jvn_
       jind(j) = 0
    END DO
    DO j = 1,njval
       DO k = 1,nratj
          IF (jlabel(k) == titlej(j)) jind(k)=j
       END DO
    END DO

    WRITE(6,'(a,i4,a)') ' Photochemistry Scheme with ',ipr,' J-values'
    DO k=1,nratj
       j = jind(k)
       IF (j == 0) THEN
          WRITE(6,'(i5,a9,f6.2,a,i4,a9)') k,jlabel(k),jfacta(k),  &
               ' has no mapping onto onto fast-JX'
       ELSE
          WRITE(6,'(i5,a9,f6.2,a,i4,a9)') k,jlabel(k),jfacta(k),  &
               ' mapped onto fast-JX:',j,titlej(j)
       END IF
    END DO


  END SUBROUTINE rd_js

  !<<<<<<<<<<<<<<<<<<<<<<<<<end core scattering subroutines<<<<<<<<<<<<<<<

END MODULE FastJX57
