 MODULE kbcon_ecmwf

 CONTAINS 

 subroutine trigg_ecmwf(mentr_rate2,PO,DNO,KIDIA, KFDIA,mgmxp,mgmzp, KLEV, PTENH,  PRMISH, PRMISSH, PGEOH, PAPH,&
                       &  PQHFL,PAHFS, PTEN, PRMIS,PRMISS, PGEO,ierr,k22,kbcon,ktop,m1,m2,m3,j,wu)  

!  ROTINA CUBASEN - MECANISMO DE DISPARO DE CONVECCAO UTILIZADO NO MODELO DO ECMWF
!  Referencia: Jakob, K. and Siebesma, P. A new subcloud Model for Mass-Flux Convection 
!  Schemes: Influence on Triggering, updrafts properties, and Climate Model.
!  Monthly Weather Review, v. 131, p. 2765-2778, 2003.
!   
!  Adaptada para o modelo CATT-BRAMS por Claudio Silva, Saulo Freitas, Nelson Pereira Filho e Rafael Mello
!
!  claudio.moises@cptec.inpe.br
!
 IMPLICIT NONE
 INTEGER           ::  KLON, KLEV, KIDIA, KFDIA, KTDIA,mgmxp,mgmzp
 INTEGER           ::  JK, JL, JKE
 INTEGER 	   ::  ierr(mgmxp), kbcon(mgmxp), ktop(mgmxp), aux(mgmxp), k22(mgmxp)
 integer           ::  m1, m2, m3, j!intent in
 real              ::  wu(mgmxp,mgmzp) !intent out
 REAL              ::  PO(mgmxp), DNO(mgmxp) 
 REAL		   ::  PTEN(mgmxp,mgmzp)   ! TEMPERATURA DO MCGAS (oC)
 REAL		   ::  PQEN(mgmxp,mgmzp)   ! UMIDADE ESPECIFICA (KG/KG)
 REAL		   ::  PTENL(mgmxp,mgmzp)  ! TEMPERATURA DO MCGAS (oC)
 REAL		   ::  PQENL(mgmxp,mgmzp)  ! UMIDADE ESPECIFICA (KG/KG)
 REAL		   ::  PQSEN(mgmxp,mgmzp)  ! UMIDADE ESPECIFICA DE SATURACAO (KG/KG)
 REAL		   ::  PQSENL(mgmxp,mgmzp) ! UMIDADE ESPECIFICA DE SATURACAO (KG/KG)
 REAL		   ::  PGEO(mgmxp,mgmzp)   ! ALTURA GEOPOTENCIAL (M) 
 REAL		   ::  PGEOL(mgmxp,mgmzp)  ! ALTURA GEOPOTENCIAL (M) 

!srf REAL	   ::  PGEOH(mgmxp,mgmzp+1)  ! ALTURA GEOPOETNCIAL NA GRADE DO MODELO DE NUVEM (M)
 REAL		   ::  PGEOH(mgmxp,mgmzp  )  ! ALTURA GEOPOETNCIAL NA GRADE DO MODELO DE NUVEM (M)

 REAL		   ::  PGEOHL(mgmxp,mgmzp+1) ! ALTURA GEOPOETNCIAL NA GRADE DO MODELO DE NUVEM (M)

! REAL		   ::  PGEOH(mgmxp,mgmzp)  ! ALTURA GEOPOETNCIAL NA GRADE DO MODELO DE NUVEM (M)
! REAL		   ::  PGEOHL(mgmxp,mgmzp) ! ALTURA GEOPOETNCIAL NA GRADE DO MODELO DE NUVEM (M)

!srf REAL	   ::  PAPH(mgmxp,mgmzp+1)  ! PRESSAO (hPa)
 REAL		   ::  PAPH(mgmxp,mgmzp  )  ! PRESSAO (hPa)

 REAL		   ::  PAPHL(mgmxp,mgmzp+1) ! PRESSAO (hPa)
 REAL		   ::  PTENH(mgmxp,mgmzp)   ! TEMPERATURA (oC) NA GRADE DO MODELO DE NUVEM
 REAL		   ::  PTENHL(mgmxp,mgmzp)  ! TEMPERATURA (oC) NA GRADE DO MODELO DE NUVEM
 REAL		   ::  PQENH(mgmxp,mgmzp)   ! UMIDADE ESPECIFICA (KG/KG) NA GRADE DO MODELO DE NUVEM
 REAL		   ::  PQENHL(mgmxp,mgmzp)  ! UMIDADE ESPECIFICA (KG/KG) NA GRADE DO MODELO DE NUVEM
 REAL		   ::  PQSENH(mgmxp,mgmzp)  ! UMIDADE ESPECIFICA DE SATURACAO (KG/KG) NA GRADE DO MODELO DE NUVEM
 REAL		   ::  PQSENHL(mgmxp,mgmzp) ! UMIDADE ESPECIFICA DE SATURACAO (KG/KG) NA GRADE DO MODELO DE NUVEM
 REAL		   ::  PRMIS(mgmxp,mgmzp)   ! RAZAO DE MISTURA DE VAPOR (KG/KG)
 REAL		   ::  PRMISL(mgmxp,mgmzp)  ! RAZAO DE MISTURA DE VAPOR (KG/KG)
 REAL	           ::  PRMISS(mgmxp,mgmzp)  ! RAZAO DE MISTURA DE SATURACAO 
 REAL	           ::  PRMISSL(mgmxp,mgmzp) ! RAZAO DE MISTURA DE SATURACAO 
 REAL              ::  PSAT(mgmxp,mgmzp)    ! PRESSAO DE SATURACAO DO VAPOR (hPa)
 REAL              ::  PSATL(mgmxp,mgmzp)   ! PRESSAO DE SATURACAO DO VAPOR (hPa)
 LOGICAL	   ::  LLFLAG(mgmxp)       ! FLAG QUE INDICA A OCORRENCIA DE CONVECCAO CUMULUS
 REAL		   ::  LATI(mgmxp,mgmzp)    ! LATITUDE INICIAL EM GRAUS (PARA LER OS ARQUIVOS DE RADIOSSONDAGENS DA
 REAL		   ::  PUREL(mgmxp,mgmzp)   ! UMIDADE RELATIVA DO AR (%)
 REAL		   ::  PURELL(mgmxp,mgmzp)  ! UMIDADE RELATIVA DO AR (%)
 REAL  		   ::  PQHFL(mgmxp)        ! FLUXO DE UMIDADE (KG/M2S)   
 REAL		   ::  PAHFS(mgmxp)        ! FLUXO DE CALOR SENSIVEL (W/M2)
 LOGICAL	   ::  LDCUM(mgmxp)        ! FLAG PARA OCORRENCIA DE CONVECCAO DEEP-CUMULUS 
 LOGICAL	   ::  LDSC(mgmxp)         ! FLAG PARA OCORRENCIA DE CONVECAO SC
 INTEGER	   ::  KCBOT(mgmxp)        ! NIVEL DA BASE DE NUVENS CUMULUS
 INTEGER	   ::  KBOTSC(mgmxp)       ! NIVEL DA BASE DE NUVENS STRATO CUMULUS
 INTEGER	   ::  KCTOP(mgmxp)        ! NIVEL DO TOPO DE NUVEM CUMULUS
 INTEGER	   ::  KDPL(mgmxp)         ! EQUIVALENTE AO k22
 REAL		   ::  PCAPE(mgmxp)        ! CAPE (J/KG)
 REAL, PARAMETER   ::  RCPD = 1004.64
 REAL              ::  PRMISSH(mgmxp,mgmzp)
 REAL              ::  PRMISH(mgmxp,mgmzp)
 REAL              ::  PRMISSHL(mgmxp,mgmzp)
 REAL              ::  PRMISHL(mgmxp,mgmzp)
 REAL              ::  mentr_rate2L(mgmxp,mgmzp)
 REAL              ::  mentr_rate2(mgmxp,mgmzp)
 REAL              ::  ZWU2H(mgmxp,mgmzp)
  
   DO JL=KIDIA,KFDIA
       DO JK=KLEV, 1, -1  ! INVERTIDO PARA FICAR COMPATIVEL COM O ESQUEMA DE GRELL !
	 JKE = KLEV+1-JK
	 PTENHL   (JL,JK)     =  PTENH     (JL,JKE)
	 PRMISHL  (JL,JK)     =  PRMISH    (JL,JKE)
	 PRMISSHL (JL,JK)     =  PRMISSH   (JL,JKE)
	 PQENHL   (JL,JK)     =  PRMISHL   (JL,JK)/(1.+PRMISHL(JL,JK))   !CONVERTENDO DE "r" PARA "q" 
	 PQSENHL  (JL,JK)     =  PRMISSHL  (JL,JK)/(1.+PRMISSHL(JL,JK)) 
	 PGEOHL   (JL,JK)     =  9.8*PGEOH (JL,JKE)  !convertendo z para geopotencial *g
	 PAPHL    (JL,JK)     =  100*PAPH  (JL,JKE)   ! convertende de hPa para Pa
	 PTENL    (JL,JK)     =  PTEN      (JL,JKE)
	 PRMISL   (JL,JK)     =  PRMIS     (JL,JKE)
	 PRMISSL  (JL,JK)     =  PRMISS    (JL,JKE)
	 PQENL    (JL,JK)     =  PRMISL    (JL,JK)/(1.+PRMISL (JL,JK))  
         PQSENL   (JL,JK)     =  PRMISSL   (JL,JK)/(1.+PRMISSL(JL,JK))
         PGEOL    (JL,JK)     =  9.8*PGEO  (JL,JKE)
	 mentr_rate2L(JL,JK)  =  mentr_rate2(JL,JKE)
      !print*,'1',JL,JK, PTENHL   (JL,JK), PRMISHL  (JL,JK),PRMISSHL (JL,JK),PQENHL   (JL,JK)  
      
      
      ENDDO
         PAPHL(JL,KLEV+1)    =  PO(JL)
	 PGEOHL(JL,KLEV+1)   =  0.9*PGEOHL(JL,klev)
       	 PQHFL(JL)           =  PQHFL(JL) 
         PAHFS(JL)           =  RCPD*DNO(JL)*PAHFS(JL)
       !print*,'11',JL,JK, PAPHL(JL,KLEV+1),PGEOHL(JL,KLEV+1) , PQHFL(JL),  PAHFS(JL) 
  
    ENDDO  	
    
 
  CALL CUBASEN(mentr_rate2L, PO, DNO, KIDIA, KFDIA, mgmxp,mgmzp,  KLEV, PTENHL,PQENHL,PQSENHL,PGEOHL, PAPHL,&
                & PQHFL, PAHFS, PTENL, PQENL,PQSENL, PGEOL,ierr,k22,kbcon,ktop,m1,m2,m3,j,wu)  


	    
 END SUBROUTINE trigg_ecmwf
 
 !------------------------------------ subrotina cubasen (Jakob e Siebesma, 2003)

 
 SUBROUTINE CUBASEN(mentr_rate2,PO, DNO, KIDIA, KFDIA, mgmxp,mgmzp, KLEV, PTENH,  PQENH, PQSENH, PGEOH, PAPH,&
                  & PQHFL, PAHFS, PTEN, PQEN,PQSEN, PGEO,ierr,k22,kcbot,kctop,m1,m2,m3,j,wu)  

 IMPLICIT NONE
 

!          THIS ROUTINE CALCULATES CLOUD BASE FIELDS
!          CLOUD BASE HEIGHT AND CLOUD TOP HEIGHT

!          A. Pier Siebesma   KNMI ********      
!          modified C Jakob (ECMWF) (01/2001) 
!          modified P Bechtold (ECMWF) (08/2002) 
!          (include cycling over levels to find unstable departure/base level+
!           mixed layer properties +w Trigger)

!          PURPOSE.
!          --------
!          TO PRODUCE CLOUD BASE AND CLOUD TOP VALUES FOR CU-PARAMETRIZATION

!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
!          IT RETURNS CLOUD FIELDS VALUES AND FLAGS AS FOLLOWS;
!                 KLAB=0 FOR STABLE LAYERS
!                 KLAB=1 FOR SUBCLOUD LEVELS
!                 KLAB=2 FOR CLOUD LEVELS LEVEL

!          METHOD.
!          --------
!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD TOP
!          (ENTRAINING PLUME, WITH ENTRAINMENT PROPORTIONAL TO (1/Z))

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS

!    INPUT PARAMETERS (REAL):

!    not used at the moment because we want to use linear intepolation
!    for fields on the half levels.

!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS           K
!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS      KG/KG
!    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
!    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2
!    *PSSTRU*       KINEMATIC surface U-MOMENTUM FLUX             (M/S)^2
!    *PSSTRV*       KINEMATIC surface V-MOMENTUM FLUX             (M/S)^2
!    *PWN*          NORMALIZED LARGE-SCALE VERTICAL VELOCITY      (M/S)
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS             PA
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
!    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
!    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2

!    UPDATED PARAMETERS (REAL):

!    *PTU*          TEMPERATURE IN UPDRAFTS                         K
!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                    KG/KG
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG
!    *PUU*          U-VELOCITY IN UPDRAFTS                         M/S
!    *PVU*          V-VELOCITY IN UPDRAFTS                         M/S

!    UPDATED PARAMETERS (INTEGER):

!    *KLAB*         FLAG KLAB=1 FOR SUBCLOUD LEVELS
!                        KLAB=2 FOR CLOUD LEVELS

!    OUTPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 
!    *LDSC*         FLAG: .TRUE. IF BL-CLOUDS EXIST

!    OUTPUT PARAMETERS (INTEGER):

!    *KCBOT*       CLOUD BASE LEVEL !    
!    *KCTOP*       CLOUD TOP LEVEL = HEIGHEST HALF LEVEL 
!                  WITH A NON-ZERO CLOUD UPDRAFT.
!    *KBOTSC*      CLOUD BASE LEVEL OF BL-CLOUDS
!    *KDPL*        DEPARTURE LEVEL
!    *PCAPE*       PSEUDOADIABATIQUE max CAPE (J/KG)

!          EXTERNALS
!          ---------
!          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO CONDENSATION IN ASCENT

!          MODIFICATIONS
!          -------------
!             92-09-21 : Update to Cy44      J.-J. MORCRETTE
!             02-11-02 : Use fixed last possible departure level and 
!                        last updraft computation level for bit-reproducibility
!                                            D.Salmond &  J. Hague
!             03-07-03 : Tuning for p690     J. Hague
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!----------------------------------------------------------------------

  REAL, PARAMETER  :: RCPD   = 1004.64
  REAL, PARAMETER  :: RETV   = 0.608 
  REAL, PARAMETER  :: RD     = 287.04 
  REAL, PARAMETER  :: RV     = 461.50  
  REAL, PARAMETER  :: RG     = 9.8200
  REAL, PARAMETER  :: RKAP   = 0.41
  INTEGER          :: NJKT1  
  INTEGER          :: NJKT2 
  REAL, PARAMETER  :: RDEPTHS = 20000.
  REAL, PARAMETER  :: ENTRPEN = 1.0E-4
  REAL, PARAMETER  :: ENTRSCV = 1.E-3
  LOGICAL          :: LMFDUDV = .false. 
  REAL, PARAMETER  :: RLMIN = 1E-5 
  REAL, PARAMETER  :: RESTT = 611.14
  REAL, PARAMETER  :: RMD = 28.9644
  REAL, PARAMETER  :: RMV = 18.0153
  REAL, PARAMETER  :: RTT = 273.16
  REAL, PARAMETER  :: R2ES = RESTT*(RD/RV)
  REAL, PARAMETER  :: R3LES = 17.269 
  REAL, PARAMETER  :: R3IES = 21.875
  REAL, PARAMETER  :: R4LES = 35.86 
  REAL, PARAMETER  :: R4IES = 7.66
  REAL, PARAMETER  :: R5LES = R3LES*(273.16-R4LES)
  REAL, PARAMETER  :: R5IES = R3IES*(273.16-R4IES)
  REAL, PARAMETER  :: RKBOL = 1.380658E-23
  REAL, PARAMETER  :: RNAVO = 6.0221367E+23
  REAL, PARAMETER  :: R = RNAVO*RKBOL
  REAL, PARAMETER  :: RLVTT = 2.5008E+6
  REAL, PARAMETER  :: RLSTT = 2.8345E+6
  REAL, PARAMETER  :: R5ALVCP = R5LES*RLVTT/RCPD
  REAL, PARAMETER  :: R5ALSCP = R5IES*RLSTT/RCPD 
  REAL, PARAMETER  :: RALVDCP = RLVTT/RCPD
  REAL, PARAMETER  :: RLMLT=RLSTT-RLVTT
  REAL, PARAMETER  :: RALSDCP = RLSTT/RCPD
  REAL, PARAMETER  :: RALFDCP = RLMLT/RCPD
  REAL, PARAMETER  :: RTWAT = 273.16 
  REAL, PARAMETER  :: RTBER = RTWAT-5.
  REAL, PARAMETER  :: RTICE = RTT-0.1
  REAL, PARAMETER  :: RTICECU = RTT-23.0
  REAL, PARAMETER  :: RTWAT_RTICECU_R  = 1./(RTWAT-RTICECU)
  REAL, PARAMETER  :: RTWAT_RTICE_R = 1./(RTWAT-RTICE)
  
  
  INTEGER :: KLON, mgmxp,mgmzp
  INTEGER :: KLEV 
  INTEGER :: KIDIA 
  INTEGER :: KFDIA 
  INTEGER :: KTDIA ! Argument NOT used
  REAL    :: PTENH(mgmxp,mgmzp) 
  
  REAL    :: PQENH(mgmxp,mgmzp)
   
  REAL    :: PGEOH(mgmxp,mgmzp+1) 
  REAL    :: PAPH(mgmxp,mgmzp+1) 
  REAL    :: PQHFL(mgmxp) 
  REAL    :: PAHFS(mgmxp) 
  REAL    :: PSSTRU(mgmxp) ! Argument NOT used
  REAL    :: PSSTRV(mgmxp) ! Argument NOT used
  REAL    :: PWN(mgmxp,mgmzp) ! Argument NOT used 
  REAL    :: PTEN(mgmxp,mgmzp)  
  REAL    :: PQEN(mgmxp,mgmzp) 
  REAL    :: PGEO(mgmxp,mgmzp) 
  REAL    :: PUEN(mgmxp,mgmzp) 
  REAL    :: PVEN(mgmxp,mgmzp)  
  REAL    :: PTU(mgmxp,mgmzp) 
  REAL    :: PQU(mgmxp,mgmzp)  
  REAL    :: PLU(mgmxp,mgmzp) 
  REAL    :: PUU(mgmxp,mgmzp) 
  REAL    :: PVU(mgmxp,mgmzp) 
  REAL    :: PWUBASE(mgmxp) 
  INTEGER :: KLAB(mgmxp,mgmzp) 
  LOGICAL :: LDCUM(mgmxp) 
  LOGICAL :: LDSC(mgmxp) 
  INTEGER :: KCBOT(mgmxp) 
  INTEGER :: KBOTSC(mgmxp) 
  INTEGER :: KCTOP(mgmxp) 
  INTEGER :: KDPL(mgmxp) 
  REAL    :: PCAPE(mgmxp) 
  INTEGER :: ICTOP(mgmxp)
  INTEGER :: ICBOT(mgmxp)
  INTEGER :: IBOTSC(mgmxp)
  INTEGER :: ILAB(mgmxp,mgmzp)
  INTEGER :: IDPL(mgmxp)  
  INTEGER :: ierr(mgmxp)   ! flag que passa ldcum (logica) para o inteiro correspondente 
  	  		  ! correspondente do esquema de grell, if ierr = 0, LDCUM = .T. --
	  		  ! CLAUDIO SILVA
 INTEGER  :: AUX(mgmxp) ! para inverter os niveis kctop e kcbot
 INTEGER  :: kbcon(mgmxp), ktop(mgmxp), k22(mgmxp)
 
 REAL  :: PQSEN(mgmxp,mgmzp)
 REAL  :: PQSENH(mgmxp,mgmzp)
 REAL  :: PRMISS(mgmxp,mgmzp)
 REAL  :: PRMIS(mgmxp,mgmzp)
 REAL  :: PSAT(mgmxp,mgmzp)
 REAL  :: PVAP(mgmxp,mgmzp)
 REAL  :: PO(mgmxp), DNO(mgmxp) 
 REAL  :: PRMISSH(mgmxp,mgmzp)
 REAL  :: PRMISH(mgmxp,mgmzp)
 
 REAL  :: mentr_rate2(mgmxp,mgmzp) ! TAXA DE ENTRANHAMENTO DE GRELL

 integer           ::  m1, m2, m3, j, jke, jkes !intent in
 real              ::  wu(mgmxp,mgmzp)
 !real              ::  wu_C(m1,m2,m3) ! matriz auxiliar para calcular a 
 !                                             ! velocidade vertical - e a energia cinetica
					      
!intent out
!             LOCAL STORAGE
!             ----- -------

  LOGICAL ::         LL_LDBASE(mgmxp),&
 & LLGO_ON(mgmxp),&
 & LLDEEP(mgmxp),    LLDCUM(mgmxp), &!*UPG PB
 & LLDSC(mgmxp),     LLFIRST(mgmxp)  
 
  LOGICAL ::     LLRESET,        LLRESETJL(mgmxp)
  INTEGER :: ICALL, IK, IKB, IS, JK, JL, JKK, JKT, JKB !*UPG PB
  INTEGER :: JKT1
  INTEGER :: JKT2   
  REAL    :: ZS(mgmxp,mgmzp)
  REAL    :: ZSENH(mgmxp,mgmzp+1)
  REAL    :: ZQENH(mgmxp,mgmzp+1)
  REAL    :: ZSUH (mgmxp,mgmzp)
  REAL    :: ZWU2H(mgmxp,mgmzp)
  REAL    :: ZBUOH(mgmxp,mgmzp)  
  REAL    :: ZQOLD(mgmxp),ZPH(mgmxp)
  REAL    :: ZMIX(mgmxp)
  REAL    :: ZDZ(mgmxp)
  REAL    :: ZCBASE(mgmxp)
  REAL    :: ZLU(mgmxp,mgmzp)
  REAL    :: ZQU(mgmxp,mgmzp)
  REAL    :: ZTU(mgmxp,mgmzp)
  REAL    :: ZUU(mgmxp,mgmzp)
  REAL    :: ZVU(mgmxp,mgmzp)  
  REAL    :: WVEL(mgmxp,mgmzp)
  REAL    :: ZCAPE(mgmxp,mgmzp) ! local for CAPE at every departure level
  REAL    :: ZBUOF, ZZ, ZC2, ZEPSADD
  REAL    :: ZRHO      ! DENSITY AT SURFACE (KG/M^3) 
  REAL    :: ZKHVFL    ! SURFACE BUOYANCY FLUX (K M/S)
  REAL    :: ZWS       ! SIGMA_W AT LOWEST MODEL HALFLEVEL (M/S)
  REAL    :: ZQEXC,ZQEXCS(mgmxp)! HUMIDITY EXCESS AT LOWEST MODEL HALFLEVEL (KG/KG)
  REAL    :: ZTEXC,ZTEXCS(mgmxp)! TEMPERATURE EXCESS AT LOWEST MODEL HALFLEVEL (K)
  REAL    :: ZEPS      ! FRACTIONAL ENTRAINMENT RATE   [M^-1]
  REAL    :: ZTVENH    ! ENVIRONMENT VIRTUAL TEMPERATURE AT HALF LEVELS (K)  
  REAL    :: ZTVUH     ! UPDRAFT VIRTUAL TEMPERATURE AT HALF LEVELS     (K)
  REAL    :: ZLGLAC    ! UPDRAFT LIQUID WATER FROZEN IN ONE LAYER
  REAL    :: zqsu
  REAL    :: ZCOR
  REAL    :: zdq
  REAL    :: zalfaw
  REAL    :: zfacw
  REAL    :: zfaci
  REAL    :: zfac
  REAL    :: zesdp
  REAL    :: zdqsdt
  REAL    :: zdtdp
  REAL    :: zdp
  REAL    :: zpdifftop
  REAL    :: zpdiffbot
  REAL    :: ZSF
  REAL    :: ZQF
  REAL    :: zaw
  REAL    :: zbw  
  REAL    :: ZTVEN1
  REAL    :: ZTVEN2
  REAL    :: ZTVU1
  REAL    :: ZTVU2 ! pseudoadiabatique T_v
  REAL    :: ZDTVTRIG(mgmxp) ! virtual temperatures
  REAL    :: ZWORK1
  REAL    :: ZWORK2! work arrays for T and w perturbations
  REAL    :: ZRCPD
  REAL    :: ZRG
  REAL    :: ZTMP
  REAL    :: ZTINY
  REAL    :: TINY
  REAL    :: ZLUMIN
  INTEGER :: INDB, INDT ! PARA LIMITAR O NJKT1 E O NJKT2
   
  REAL :: ZWS1D(mgmxp),ZRHO1D(mgmxp)!srf 
    
!----------------------------------------------------------------------
!     0.           INITIALIZE CONSTANTS AND FIELDS
!                  -------------------------------
!----------------------------------------------------------------------
!      CONTEUDO MINIMO DE AGUA PARA ENCONTRAR A BASE DA NUVEM! CLAUDIO SILVA


 	ZLUMIN = 5.e-4
	ZC2    = 0.55
	ZAW    = 1.0
	ZBW    = 2.0
        ZEPSADD= 1.E-4

 	DO JL=KIDIA,KFDIA
 	    PWUBASE(JL)=0.0
  	    LLGO_ON(JL)=.TRUE.
  	    LLFIRST(JL)=.TRUE.
  	    KDPL(JL)=KLEV
  	    ZTEXCS(JL)=0.2
  	    ZQEXCS(JL)=1.E-4
	ENDDO

!!!!!!!!!!!!!!! Para determinar os limites de busca da base e topo ! 
    
      INDB = 0
      INDT = 0
      DO JL = KIDIA, KFDIA 
        DO JK = 1, KLEV
	  IF(PAPH(JL,JK).GT.60.AND.INDB.EQ.0) THEN
             JKT2 = JK   
             INDB = 1
	  ENDIF

          IF(PAPH(JL,JK).GT.350E2.AND.INDT.EQ.0) THEN  ! o limite inferior de busca da base
             JKT1 = JK                                 ! da nuvem esta muito profundo,
	                                               ! testar 600
             INDT = 1
          ENDIF
	
        ENDDO
      ENDDO
      
!!!!!!!!!!!!!!!!!!!
	
    ZRG=1.0/RG
    ZRCPD=1.0/RCPD
    ZTINY=TINY(ZRG)

	DO JK=1,KLEV
  	   DO JL=KIDIA,KFDIA
    		ZTU(JL,JK) = PTENH(JL,JK) !PTU(JL,JK)
    		ZQU(JL,JK) = PQENH(JL,JK) !PQU(JL,JK)
    		ZLU(JL,JK) = PLU(JL,JK)
    		ZUU(JL,JK) = PUU(JL,JK)
    		ZVU(JL,JK) = PVU(JL,JK)
    		ILAB(JL,JK)= KLAB(JL,JK)
    		ZCAPE(JL,JK)= 0.0
 	    ENDDO
	ENDDO
!----------------------------------------------------------------------
!       -----------------------------------------------------------
!       1.1  PREPARE FIELDS ON HALF LEVELS BY LINEAR INTERPOLATION
!             OF SPECIFIC HUMIDITY AND STATIC ENERGY
!       -----------------------------------------------------------

 	DO JK=1,KLEV
  	   DO JL=KIDIA,KFDIA
    		ZWU2H(JL,JK) = 0.0
   		ZS   (JL,JK) = RCPD*PTEN(JL,JK) + PGEO(JL,JK)
    		ZQENH(JL,JK) = PQENH(JL,JK)
    		ZSENH(JL,JK) = RCPD*PTENH(JL,JK)+PGEOH(JL,JK)
  	   ENDDO
	ENDDO
                
		          
	DO JKK = KLEV,JKT1,-1  ! Big external loop for level testing:
                               ! find first departure level that produces deepest cloud top
                               ! or take surface level for shallow convection and Sc
   ! 
   !        ---------------------------------------------------------
   !        1.2    INITIALISE FIELDS AT DEPARTURE HALF MODEL LEVEL
   !        ---------------------------------------------------------
   !
 		 IS=0
		 DO JL=KIDIA,KFDIA
    			IF (LLGO_ON(JL)) THEN
      			IS=IS+1
		        IDPL(JL)      = JKK      ! departure level
		        ICBOT  (JL)   = JKK      ! cloud base level for convection, (-1 if not found)
		        IBOTSC (JL)   = KLEV-1   ! sc    base level for sc-clouds , (-1 if not found)
		        ICTOP(JL)     = KLEV-1   ! cloud top for convection (-1 if not found) --->>> Mas o topo fica embaixo ??
		        LLDCUM(JL)    = .FALSE.  ! on exit: true if cloudbase=found
		        LLDSC (JL)    = .FALSE.  ! on exit: true if cloudbase=found
		        LL_LDBASE(JL) =.FALSE.   ! on exit: true if cloudbase=found
		        ZDTVTRIG(JL)  = 0.0
  		        ZUU(JL,JKK)   = PUEN(JL,JKK)*(PAPH(JL,JKK+1)-PAPH(JL,JKK))
		        ZVU(JL,JKK)   = PVEN(JL,JKK)*(PAPH(JL,JKK+1)-PAPH(JL,JKK))
		        ENDIF 
   	   	 ENDDO

 		 IF(IS /= 0) THEN

			IF(JKK == KLEV) THEN

			      DO JL=KIDIA,KFDIA
        			   IF (LLGO_ON(JL)) THEN
        			   ZRHO  = PAPH(JL,JKK+1)/(RD*(PTEN(JL,JKK)*(1.+RETV*PQEN(JL,JKK))))
			           !- buoyancy flux (H+LE)
				   ZKHVFL= (PAHFS(JL)/RCPD+RETV*PTEN(JL,JKK)*PQHFL(JL))/ZRHO
				   !-convective-scale velocity w*
				   ZWS=0.001-1.5*RKAP*ZKHVFL*PGEOH(JL,KLEV)/PTEN(JL,KLEV)
				   !srf
				   ZWS1D(JL)=1.2*ZWS**.3333
				   ZRHO1D(JL)=ZRHO
				   !srf
				   				   
				   IF(ZWS >= ZTINY) THEN
				   !-convective-scale velocity w*
				   ZWS=1.2*ZWS**.3333
			           ILAB(JL,JKK)= 1
			           !- temperature excess 
				   ZTEXC     = MAX(-1.5*PAHFS(JL)/(ZRHO*ZWS*RCPD),0.0)
     		          	   !- moisture  excess
				   ZQEXC     = MAX(-1.5*PQHFL(JL)/(ZRHO*ZWS),0.0)
					  
				   ZTEXCS(JL)=ZTEXC
  			           ZQEXCS(JL)=ZQEXC
                	           
				   !- initial values for updrafts
				   !- humidty
				   ZQU (JL,JKK) = ZQENH(JL,JKK) + ZQEXC
			           !- static energy (?)
				   ZSUH (JL,JKK) = ZSENH(JL,JKK) + RCPD*ZTEXC
			           !- temperature
				   ZTU (JL,JKK) = (ZSENH(JL,JKK)-PGEOH(JL,JKK))/RCPD + ZTEXC
				   !- kinetic energy
			           ZWU2H(JL,JKK) = ZWS**2				   
       		                   !-  determine buoyancy at lowest half level
				   ZTVENH            = (1.0+RETV*ZQENH(JL,JKK)) &
             				& *(ZSENH(JL,JKK)-PGEOH(JL,JKK))/RCPD  
            			   ZTVUH             = (1.0+RETV*ZQU(JL,JKK))*ZTU(JL,JKK)
			           ZBUOH(JL,JKK) = (ZTVUH-ZTVENH)*RG/ZTVENH
				  
				  ELSE
          			  LLGO_ON(JL)=.FALSE.      ! non-convective point
          			  
				  ENDIF
        			  ENDIF
      			      ENDDO
   
    		        ELSE ! se jkk /= klev

    			DO JL=KIDIA,KFDIA
        		      IF (LLGO_ON(JL)) THEN
         		      !- air density
			      ZRHO  = PAPH(JL,JKK+1)/(RD*(PTEN(JL,JKK)*(1.+RETV*PQEN(JL,JKK))))
			      ILAB(JL,JKK)= 1
			      !- temperature excess
			      ZTEXC=.2
			      !- moisture  excess
			      ZQEXC=1.E-4

			     !srf--------------------------------------------------- 
			      !- temperature excess 			     
			      if(jkk == klev -1) then
			       ZTEXC     = min(-1.5*PAHFS(JL)/(ZRHO1D(JL)*ZWS1d(JL)*RCPD),0.2)
  		              !- moisture  excess
			       ZQEXC     = min(-1.5*PQHFL(JL)/(ZRHO1D(JL)*ZWS1D(JL)),1.E-4)
			      endif
			     !srf--------------------------------------------------- 
			     
			      !- initial values for updrafts
			      !- humidty
			      ZQU (JL,JKK) = ZQENH(JL,JKK) + ZQEXC
			      !- static energy (?)
			      ZSUH (JL,JKK) = ZSENH(JL,JKK) + RCPD*ZTEXC
			      !- temperature
			      ZTU (JL,JKK) = (ZSENH(JL,JKK)-PGEOH(JL,JKK))*ZRCPD + ZTEXC
			 ! construct mixed layer for parcels emanating in lowest 60 hPa
			 ! Esabilizando ainda mais as parcelas, fazendo a Camada de mistura ateh 80 hPa
          		     IF (PAPH(JL,KLEV+1)-PAPH(JL,JKK-1)<80.E2) THEN
            		     ZQU(JL,JKK) =0.0
		             ZSUH(JL,JKK)=0.0
		             ZWORK1      =0.0
			            DO JK=JKK+1,JKK-1,-1
              				IF( ZWORK1 < 50.E2 ) then
			                ZWORK2=PAPH(JL,JK)-PAPH(JL,JK-1)
			                ZWORK1      =ZWORK1+ZWORK2
			                ZQU(JL,JKK) =ZQU(JL,JKK) +ZQENH(JL,JK)*ZWORK2
			                ZSUH(JL,JKK)=ZSUH(JL,JKK)+ZSENH(JL,JK)*ZWORK2
					ENDIF
            			    ENDDO
	
        		     ZQU(JL,JKK) =ZQU(JL,JKK) /ZWORK1+ZQEXC
		             ZSUH(JL,JKK)=ZSUH(JL,JKK)/ZWORK1+RCPD*ZTEXC
		             ZTU(JL,JKK) =(ZSUH(JL,JKK)-PGEOH(JL,JKK))/RCPD+ZTEXC
		             ENDIF
		             !- kinetic energy
			     ZWU2H(JL,JKK) = 1.0 
			      !
			      !  determine buoyancy at lowest half level
			      !
          		     ZTVENH            = (1.0 +RETV*ZQENH(JL,JKK)) &
			           & *(ZSENH(JL,JKK)-PGEOH(JL,JKK))/RCPD  
		             ZTVUH             = (1.0 +RETV*ZQU(JL,JKK))*ZTU(JL,JKK)
           	             ZBUOH(JL,JKK) = (ZTVUH-ZTVENH)*RG/ZTVENH
      			     ENDIF
     			ENDDO
    			ENDIF
		   ENDIF

   !----------------------------------------------------------------------
   !     2.0          DO ASCENT IN SUBCLOUD AND LAYER,
   !                  CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
   !                  ADJUST T,Q AND L ACCORDINGLY IN *CUADJTQ*,
   !                  CHECK FOR BUOYANCY AND SET FLAGS
   !                  -------------------------------------
   !       ------------------------------------------------------------
   !        1.2  DO THE VERTICAL ASCENT UNTIL VELOCITY BECOMES NEGATIVE
   !       ------------------------------------------------------------
  		 DO JK=JKK-1,JKT2,-1
    		 IS=0
			IF(JKK==KLEV) THEN ! 1/z mixing for shallow
			      DO JL=KIDIA,KFDIA
        			        IF (LLGO_ON(JL)) THEN
          				IS         = IS+1
				        ZDZ(JL)        = (PGEOH(JL,JK) - PGEOH(JL,JK+1))*ZRG
				        
					zeps = mentr_rate2(jl,jk)
					
					!ZEPS = ZC2/(PGEO(JL,JK)*ZRG + ZDZ(jl)) + ZEPSADD
			                !mentr_rate2(jl,jk) = ZC2/(PGEO(JL,JK)*ZRG + ZDZ(jl)) + ZEPSADD
			            	!mentr_rate2(jl,jk)=zeps
					
					ZMIX(JL)   = 0.5*ZDZ(JL)*ZEPS
					ZQF = (PQENH(JL,JK+1) + PQENH(JL,JK))*0.5
				        ZSF = (ZSENH(JL,JK+1) + ZSENH(JL,JK))*0.5
				        ZTMP = 1.0/(1.0+ZMIX(JL))
				        ZQU(JL,JK)= (ZQU(JL,JK+1)*(1.0-ZMIX(JL))&
				           & +2.0*ZMIX(jl)*ZQF) * ZTMP  
				        ZSUH (JL,JK)= (ZSUH(JL,JK+1)*(1.0-ZMIX(JL))&
				           & +2.0*ZMIX(jl)*ZSF) * ZTMP  
			                ZQOLD(JL)  = ZQU(JL,JK)
				        ZTU (JL,JK) = (ZSUH(JL,JK)-PGEOH(JL,JK))*ZRCPD
				        ZPH  (JL)    = PAPH(JL,JK)
				        ENDIF
      			      ENDDO
    			ELSE

		        DO JL=KIDIA,KFDIA
        		      IF (LLGO_ON(JL)) THEN
	         	      IS         = IS+1
			      ZDZ(JL)        = (PGEOH(JL,JK) - PGEOH(JL,JK+1))*ZRG
		              ZQF = (PQENH(JL,JK+1) + PQENH(JL,JK))*0.5
		              ZSF = (ZSENH(JL,JK+1) + ZSENH(JL,JK))*0.5

!srf			      ZMIX(JL)= ENTRPEN           *(PGEOH(JL,JK)-PGEOH(JL,JK+1))/RG
			      ZMIX(JL)= mentr_rate2(jl,jk)*(PGEOH(JL,JK)-PGEOH(JL,JK+1))/RG

			      ZQU(JL,JK)= ZQU(JL,JK+1)*(1.0 -ZMIX(JL))+ ZQF*ZMIX(JL)
			      ZSUH(JL,JK)= ZSUH(JL,JK+1)*(1.0 -ZMIX(JL))+ ZSF*ZMIX(JL)
		              ZQOLD(JL)  = ZQU(JL,JK)
		              ZTU (JL,JK) = (ZSUH(JL,JK)-PGEOH(JL,JK))*ZRCPD
			      ZPH  (JL)    = PAPH(JL,JK)
        		     ENDIF
      			ENDDO
    			ENDIF


               IF (IS == 0) EXIT
         	IK=JK
    		ICALL=1
                CALL CUADJTQ(KIDIA,KFDIA,mgmxp,mgmzp,KTDIA,KLEV,IK,ZPH,ZTU,ZQU,LLGO_ON,ICALL)  
    		 DO JL=KIDIA,KFDIA
      		 	IF(LLGO_ON(JL)) THEN
   
			! add condensation to water
          	        ZDQ=MAX(ZQOLD(JL)-ZQU(JL,JK),0.0)
			
			! O CALCULO DA BASE DA NUVEM DEPENDE EXCLUSIAMENTE DE ZLU !
		        ZLU(JL,JK)=ZLU(JL,JK+1)+ZDQ
        
			! freezing
          	        ZLGLAC=ZDQ*((1.0-FOEALFCU(ZTU(JL,JK)))-&
		         & (1.0-FOEALFCU(ZTU(JL,JK+1))))  

			! pseudo-microphysics
        			IF(JKK==KLEV) THEN  ! no precip for shallow
			        ZLU(JL,JK)=MIN(ZLU(JL,JK),5.E-3)
	  
   			!* chose a more pseudo-adiabatic formulation as original overestimates
		        !* water loading efect and therefore strongly underestimates cloud thickness
        
				ELSE
			        ZLU(JL,JK)=0.5*ZLU(JL,JK) 
        		ENDIF
   
		       ! update dry static energy after condensation + freezing
		        ZSUH(JL,JK) = RCPD*(ZTU(JL,JK)+RALFDCP*ZLGLAC)+PGEOH(JL,JK)
    
		      ! Buoyancy on half and full levels
        		ZTVUH           = (1.0+RETV*ZQU(JL,JK)-ZLU(JL,JK))*ZTU(JL,JK)&
		         & +RALFDCP*ZLGLAC  
		       
		        ZTVENH          = (1.0+RETV*ZQENH(JL,JK)) &
		         & *(ZSENH(JL,JK)-PGEOH(JL,JK))*ZRCPD  
		       
		        ZBUOH(JL,JK)   = (ZTVUH-ZTVENH)*RG/ZTVENH
		        ZBUOF          = (ZBUOH(JL,JK) + ZBUOH(JL,JK+1))*0.5
   
   		     ! solve kinetic energy equation
                        ZTMP=1.0/(1.0+2.0*ZBW*ZMIX(jl))
		        ZWU2H(JL,JK) = (ZWU2H(JL,JK+1)*(1.0-2.0*ZBW*ZMIX(jl))&
		         & +2.0*ZAW*ZBUOF*ZDZ(jl))*ZTMP 
	 
		     ! compute pseudoadiabatique CAPE for diagnostics
                        ZTVU2 = ZTU(JL,JK)  *(1.0 +RETV*ZQU(JL,JK))
		        ZTVEN2= PTENH(JL,JK)*(1.0 +RETV*PQENH(JL,JK))
	
        		IF (JK == JKK-1) THEN
         		ZTVU1  = ZTVU2
		        ZTVEN1 = ZTVEN2
		        ENDIF
		        ZBUOF = (ZTVU2+ZTVU1-ZTVEN1-ZTVEN2)/ZTVEN2
		        ZBUOF = ZBUOF*ZDZ(JL)*RG
		        ZCAPE(JL,JKK)  = ZCAPE(JL,JKK) + MAX(0.0,ZBUOF)
		        ZTVU1=ZTVU2
		        ZTVEN1=ZTVEN2
 
	
               IF(ZLU(JL,JK) >ZLUMIN.AND.ILAB(JL,JK+1)==1) THEN
                IK=JK+1
          	ZQSU=FOEEWM(ZTU(JL,IK))/PAPH(JL,IK)
          	ZQSU=MIN(0.5,ZQSU)
          	ZCOR=1.0/(1.0-RETV*ZQSU)
	        ZQSU=ZQSU*ZCOR
          	ZDQ=MIN(0.,ZQU(JL,IK)-ZQSU)
          	ZALFAW=FOEALFA(ZTU(JL,IK))
          	ZFACW=R5LES/((ZTU(JL,IK)-R4LES)**2)
          	ZFACI=R5IES/((ZTU(JL,IK)-R4IES)**2)
          	ZFAC=ZALFAW*ZFACW+(1.-ZALFAW)*ZFACI
          	ZESDP=FOEEWM(ZTU(JL,IK))/PAPH(JL,IK)
          	ZCOR=1.0/(1.0-RETV*ZESDP)
          	ZDQSDT=ZFAC*ZCOR*ZQSU
          	ZDTDP=RD*ZTU(JL,IK)/(RCPD*PAPH(JL,IK))
          	ZDP=ZDQ/(ZDQSDT*ZDTDP)
          	ZCBASE(JL)=PAPH(JL,IK)+ZDP
      	       
	       ! chose nearest half level as cloud base
                ZPDIFFTOP=ZCBASE(JL)-PAPH(JL,JK)
          	ZPDIFFBOT=PAPH(JL,JK+1)-ZCBASE(JL)
           
          	IF(ZPDIFFTOP > ZPDIFFBOT.AND.ZWU2H(JL,JK+1)>0.0) THEN
            	JKB=MIN(KLEV-1,JK+1)
            	ILAB(JL,JKB)=2 !*UPG
            	ILAB(JL,JK)=2
            	LL_LDBASE(JL) =.TRUE.
            	LLDSC(JL)   =.TRUE.
            	IBOTSC(JL) =JKB
            	ICBOT(JL)  =JKB
            	ZLU(JL,JK+1) = RLMIN

		ELSEIF(ZPDIFFTOP <= ZPDIFFBOT.AND.ZWU2H(JL,JK)>0.0) THEN
            	ILAB(JL,JK)=2
            	LL_LDBASE(JL) =.TRUE.
            	LLDSC(JL)   =.TRUE.
            	IBOTSC(JL) =JK
            	ICBOT(JL)  =JK

          	ENDIF
          	JKB=ICBOT(JL)
        	ENDIF
   
   ! decide on presence of convection, cloud base and cloud top based on
   ! kinetic energy
        	IF (ZWU2H(JL,JK) < 0.0) THEN
          	LLGO_ON(JL) = .FALSE. 
          		IF (ZLU(JL,JK+1)>0.0) THEN  ! AQUI ZLU PRECISA SER MAIOR QUE ZERO
	                ICTOP(JL)   = JK
            		LLDCUM(JL)   = .TRUE.
          		ELSE
            		LLDCUM(JL)   = .FALSE.
          		ENDIF
        	ELSE
        	IF (ZLU(JL,JK)>ZLUMIN) THEN
	    	ILAB(JL,JK) = 2
          	ELSE
            	ILAB(JL,JK) = 1
          	ENDIF
        	ENDIF
      ENDIF
    ENDDO
   
  IF(LMFDUDV.AND.JKK==KLEV) THEN
      DO JL=KIDIA,KFDIA
        IF(.NOT.LL_LDBASE(JL).AND.LLGO_ON(JL)) THEN
          ZUU(JL,JKK)=ZUU(JL,JKK)+PUEN(JL,JK)*(PAPH(JL,JK+1)-PAPH(JL,JK))
          ZVU(JL,JKK)=ZVU(JL,JKK)+PVEN(JL,JK)*(PAPH(JL,JK+1)-PAPH(JL,JK))
        ENDIF
      ENDDO
    ENDIF
   
  ENDDO
   
   
  IF( JKK==KLEV) THEN
      ! set values for departure level for PBL clouds = first model level
    DO JL=KIDIA,KFDIA
      LDSC(JL)  = LLDSC(JL)
      IF(LDSC(JL)) THEN
        KBOTSC(JL)= IBOTSC(JL)
      ELSE
        KBOTSC(JL)=-1
      ENDIF
    
      LLGO_ON(JL) = .FALSE.
      JKT=ICTOP(JL)
      JKB=ICBOT(JL)
      LLDEEP(JL)=PAPH(JL,JKB)-PAPH(JL,JKT)>RDEPTHS
     
      IF(LLDEEP(JL)) LLDCUM(JL)=.FALSE. ! no deep allowed for KLEV
      lldeep(jl)=.false.                ! for deep convection start only at level KLEV-1
                                        ! and form mixed layer, so go on
                                        ! test further for deep convective columns as not yet found
      
      IF ( LLDEEP(JL) ) LLFIRST(JL)=.FALSE.
      LLGO_ON(JL) = .NOT.LLDEEP(JL)
      
      IF(LLDCUM(JL)) THEN
        KCBOT(JL)= ICBOT(JL)
        KCTOP(JL)= ICTOP(JL)
        KDPL(JL)  = IDPL(JL)
        LDCUM(JL) = LLDCUM(JL)
        PWUBASE(JL)=SQRT(MAX(ZWU2H(JL,JKB),0.0))
      ELSE
        KCTOP(JL)=-1
        KCBOT(JL)=-1
        KDPL(JL) =KLEV-1
        LDCUM(JL)=.FALSE.
        PWUBASE(JL)=0.0
      ENDIF
    ENDDO
    DO JK=KLEV,1,-1
      DO JL=KIDIA,KFDIA
        JKT=ICTOP(JL)
        IF ( JK>=JKT ) THEN
          KLAB(JL,JK)=ILAB(JL,JK)
          PTU(JL,JK)=ZTU(JL,JK)
          PQU(JL,JK)=ZQU(JL,JK)
          PLU(JL,JK)=ZLU(JL,JK)
        ENDIF
      ENDDO
    ENDDO
    IF(LMFDUDV) THEN
      DO JL=KIDIA,KFDIA
        IF(LDCUM(JL)) THEN
          IKB=KCBOT(JL)
          ZZ=1.0 /(PAPH(JL,JKK+1)-PAPH(JL,IKB))
          PUU(JL,JKK)=ZUU(JL,JKK)*ZZ
          PVU(JL,JKK)=ZVU(JL,JKK)*ZZ
        ENDIF
      ENDDO
    ENDIF
  ENDIF
   
  IF( JKK < KLEV ) THEN
    LLRESET=.FALSE.
    DO JL=KIDIA,KFDIA
      IF ( .NOT.LLDEEP(JL) ) THEN
        JKT=ICTOP(JL)
        JKB=ICBOT(JL)
        LLDEEP(JL)=PAPH(JL,JKB)-PAPH(JL,JKT)>=RDEPTHS 
      ENDIF
      LLRESETJL(JL)=LLDEEP(JL).AND.LLFIRST(JL)
      LLRESET=LLRESET.OR.LLRESETJL(JL)
    ENDDO

    IF(LLRESET) THEN
      DO JK=KLEV,1,-1
        DO JL=KIDIA,KFDIA
         ! keep first departure level that produces deep cloud
          IF ( LLRESETJL(JL) ) THEN 
            JKT=ICTOP(JL)
            JKB=IDPL(JL)
            IF ( JK<=JKB .AND. JK>=JKT ) THEN
              KLAB(JL,JK)=ILAB(JL,JK)
              PTU(JL,JK)=ZTU(JL,JK)
              PQU(JL,JK)=ZQU(JL,JK)
              PLU(JL,JK)=ZLU(JL,JK)
            ELSE 
              KLAB(JL,JK)=1
              PTU(JL,JK)=PTENH(JL,JK)
              PQU(JL,JK)=PQENH(JL,JK)
              PLU(JL,JK)=0.0
            ENDIF
            IF ( JK<JKT ) KLAB(JL,JK)=0
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    DO JL=KIDIA,KFDIA
      IF (LLDEEP(JL) .AND. LLFIRST(JL)) THEN
        KDPL(JL)  = IDPL(JL)
        KCTOP(JL) = ICTOP(JL)
        KCBOT(JL) = ICBOT(JL)
        LDCUM(JL) = LLDCUM(JL)
        LDSC(JL)  = .FALSE.
        KBOTSC(JL)= -1
        JKB=KCBOT(JL)
        PWUBASE(JL)=SQRT(MAX(ZWU2H(JL,JKB),0.0))
	
!  no initialization of wind for deep here, this is done in
!  CUINI and CUASCN
        LLFIRST(JL)=.FALSE.
      ENDIF
      LLGO_ON(JL) = .NOT.LLDEEP(JL)
    ENDDO
  ENDIF
ENDDO ! end of big loop for search of departure level     



! chose maximum CAPE value
 DO JL=KIDIA,KFDIA
  PCAPE(JL) = MAXVAL(ZCAPE(JL,:))
 ENDDO
 
!!!! PASSANDO O LDCUM DO ESQUEMA "CUBASEN" PARA IERR DO GRELL
 
 DO JL = KIDIA, KFDIA
   IF (LDCUM(JL)) THEN
      ierr(jl) = 0
   ELSE
      ierr(jl) = 2
   ENDIF

 ENDDO
 
 
!!!!!! INVERTENDO A ORDEM DE KCBOT E KCTOP
   
    DO JL = KIDIA, KFDIA
     if(ierr(jl).eq.0) then 
       KCTOP(JL) = KLEV+1-KCTOP(JL)
       KCBOT(JL) = KLEV+1-KCBOT(JL)
       K22(JL)   = KLEV+1-KDPL(JL)
     ELSE
       KCTOP(JL) = -1
       KCBOT(JL) = -1  
       K22 (JL)  = -1    
     ENDIF
    ENDDO
    
  
  ! A variavel ZWU2H(JL,JK) corresponde a energia cinetica, 
  ! como pode ser visto pelo calculo da velocidade na base das nuvens
  ! PWUBASE(JL)=SQRT(MAX(ZWU2H(JL,JKB),0.0))
  ! portanto, faz-se necessario calcular esta velocidade para toda a coluna.
  ! mas e importante notar que nao podemos ter valores negativos da energia cinetica!
       
   DO JL=KIDIA,KFDIA 
      if(ierr(jl).eq.0) then 
        DO JK=KLEV,1,-1
           JKE = KLEV+1-JK  
	   
           !wu_C(JKE,JL,J) = ZWU2H(JL,JK)
	   
	   !if(jke.ge.kctop(jl)) then ! zero a matriz do topo para cima, com isso 
	   !wu_C(JKE,JL,J) = 0. ! excluimos a energia cinetica negativa (que o 
	   !endif                     ! algoritmo permite, embora nao seja possivel 
	                             ! na realidade.
          !wu(JKE,JL,J) = sqrt(max(2*wu_C(JKE,JL,J),0.0))
          wu(JL,JKE) = sqrt(max(2*ZWU2H(JL,JK),0.0))
	
	enddo
      endif
      if(ierr(jl).eq.2) then 
	DO JK=KLEV,1,-1
           JKE = KLEV+1-JK
           !wu_C(JKE,JL,J) = 0.0
	   !wu(JKE,JL,J) = sqrt(max(2*wu_C(JKE,JL,J),0.0))
	   wu(JL,JKE) = 0.
        ENDDO
      endif   
    enddo
    
      	   
    	   
   

END SUBROUTINE CUBASEN

!!----------------------------------------------------------------
!!-----------------SUBROTINA QUE AJUSTA OS CAMPOS DE UMIDADE E ETEMPERATURA
!!----------------- DEVIDO A CONDENSACAO
!!--------------------------------------------------------------------------------
    SUBROUTINE CUADJTQ(KIDIA,KFDIA,mgmxp,mgmzp,KTDIA,KLEV,KK,PSP,PT,PQ,LDFLAG,KCALL)
      
      !          M.TIEDTKE         E.C.M.W.F.     12/89 
      !          MODIFICATIONS
      !          -------------
      !          D.SALMOND         CRAY(UK))      12/8/91
      !          J.J. MORCRETTE    ECMWF          92-09-18   Update to Cy44
      !          J.F. MAHFOUF      ECMWF          96-06-11   Smoothing option
      !          PURPOSE.
      !          --------
      !          TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT
      
      !          INTERFACE
      !          ---------
      !          THIS ROUTINE IS CALLED FROM SUBROUTINES:
      !              *COND*     (T AND Q AT CONDENSATION LEVEL)
      !              *CUBASE*   (T AND Q AT CONDENSATION LEVEL)
      !              *CUASC*    (T AND Q AT CLOUD LEVELS)
      !              *CUINI*    (ENVIRONMENTAL T AND QS VALUES AT HALF LEVELS)
      !              *CUSTRAT*  (T AND Q AT CONDENSATION LEVEL)
      !          INPUT ARE UNADJUSTED T AND Q VALUES,
      !          IT RETURNS ADJUSTED VALUES OF T AND Q
      
      !     PARAMETER     DESCRIPTION                                   UNITS
      !     ---------     -----------                                   -----
      !     INPUT PARAMETERS (INTEGER):
      
      !    *KIDIA*        START POINT
      !    *KFDIA*        END POINT
      !    *KLON*         NUMBER OF GRID POINTS PER PACKET
      !    *KTDIA*        START OF THE VERTICAL LOOP
      !    *KLEV*         NUMBER OF LEVELS
      !    *KK*           LEVEL
      !    *KCALL*        DEFINES CALCULATION AS
      !                      KCALL=0  ENV. T AND QS IN*CUINI*
      !                      KCALL=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
      !                      KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)
      
      !     INPUT PARAMETERS (LOGICAL):
      
      !    *LDLAND*       LAND-SEA MASK (.TRUE. FOR LAND POINTS)
      
      !     INPUT PARAMETERS (REAL):
   
      !    *PSP*          PRESSURE                                        PA
      
      !     UPDATED PARAMETERS (REAL):
      
      !    *PT*           TEMPERATURE                                     K
      !    *PQ*           SPECIFIC HUMIDITY                             KG/KG
      
      
      !          EXTERNALS   
      !          ---------
      !          3 LOOKUP TABLES ( TLUCUA, TLUCUB, TLUCUC )
      !          FOR CONDENSATION CALCULATIONS.
      !          THE TABLES ARE INITIALISED IN *SUPHEC*.
      
      !----------------------------------------------------------------------
            IMPLICIT LOGICAL (L)
            REAL, PARAMETER :: RETV    = 0.608 
	    REAL, PARAMETER :: RTT     = 273.16 
	    REAL, PARAMETER :: RESTT   = 611.14
	    REAL, PARAMETER :: RMD     = 28.9644
	    REAL, PARAMETER :: RMV     = 18.0153
	    REAL, PARAMETER :: RV      = 461.50
	    REAL, PARAMETER :: RD      = 287.04
	    REAL, PARAMETER :: R2ES    = RESTT*(RD/RV)
            REAL, PARAMETER :: R3LES   = 17.269 
            REAL, PARAMETER :: R3IES   = 21.875
            REAL, PARAMETER :: R4LES   = 35.86 
            REAL, PARAMETER :: R4IES   = 7.66
            REAL, PARAMETER :: R5LES   = R3LES*(273.16-R4LES)
            REAL, PARAMETER :: R5IES   = R3IES*(273.16-R4IES)
	    REAL, PARAMETER :: RCPD    = 3.5*RD
	    REAL, PARAMETER :: RLVTT   = 2.5008E+6
	    REAL, PARAMETER :: RLSTT   = 2.8345E+6
            REAL, PARAMETER :: R5ALVCP = R5LES*RLVTT/RCPD
            REAL, PARAMETER :: R5ALSCP = R5IES*RLSTT/RCPD 
            REAL, PARAMETER :: RALVDCP = RLVTT/RCPD
	    REAL, PARAMETER :: RLMLT   = RLSTT-RLVTT
            REAL, PARAMETER :: RALSDCP = RLSTT/RCPD
            REAL, PARAMETER :: RALFDCP = RLMLT/RCPD
            REAL, PARAMETER :: RTWAT   = 273.16 
            REAL, PARAMETER :: RTBER   = RTWAT-5.
            REAL, PARAMETER :: RTICE   = RTT-0.1
            REAL, PARAMETER :: RTICECU = RTT-23.0
 	    LOGICAL         :: LPHYLIN = .true. ! if not linearized physics is used
            REAL, PARAMETER :: RLPTRC  = 266.425 
 	    REAL, PARAMETER :: RLPAL1  = 0.15 
            REAL, PARAMETER :: RLPAL2  = 20. 
	    REAL, PARAMETER :: ZQMAX   = 0.5
      	    LOGICAL    LDFLAG(mgmxp)
            REAL ::    PT(mgmxp,mgmzp)
	    REAL ::    PQ(mgmxp,mgmzp)
	    REAL ::    PSP(mgmxp)
            REAL ::    ZCOND(mgmxp)
	    REAL ::    ZQP(mgmxp)
	    REAL ::    ZQSAT1
	 
!*********************************************
            IF (LPHYLIN) THEN ! para usar a fIsica nAo linearizada (ClAudio Silva)
!*********************************************                 
!     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY, 
!                  -----------------------------------------------------
              IF (KCALL.EQ.1 ) THEN
                ISUM=0
                DO JL=KIDIA,KFDIA
                    IF(LDFLAG(JL)) THEN
		      ZQP(JL)=1./PSP(JL)
		      ZQSAT=FOEEWM(PT(JL,KK))*ZQP(JL)
                      ZQSAT=MIN(0.5,ZQSAT)
                      ZCOR=1./(1.-RETV *ZQSAT)
                      ZQSAT=ZQSAT*ZCOR
		      ZCOND(JL)=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
                      ZCOND(JL)=MAX(ZCOND(JL),0.)
                      PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND(JL)
                      PQ(JL,KK)=PQ(JL,KK)-ZCOND(JL)
	            IF(ZCOND(JL).NE.0.0) ISUM=ISUM+1
                    ELSE
                    ZCOND(JL)=0.0
                    ENDIF
                ENDDO

                IF(ISUM.EQ.0) GO TO 230
                DO JL=KIDIA,KFDIA
                  IF(LDFLAG(JL).AND.ZCOND(JL).NE.0.) THEN
                    ZQSAT=FOEEWM(PT(JL,KK))*ZQP(JL)
                    ZQSAT=MIN(0.5,ZQSAT)
                    ZCOR=1./(1.-RETV*ZQSAT)
                    ZQSAT=ZQSAT*ZCOR
                    ZCOND1=(PQ(JL,KK)-ZQSAT) /(1.+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
                    PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND1
                    PQ(JL,KK)=PQ(JL,KK)-ZCOND1
                  ENDIF
                ENDDO
      
       230     CONTINUE
      
	      ENDIF
      
              IF(KCALL.EQ.2) THEN
                ISUM=0
               DO JL=KIDIA,KFDIA
                  IF(LDFLAG(JL)) THEN
                    ZQP(JL)=1./PSP(JL)
                    ZQSAT=FOEEWM(PT(JL,KK))*ZQP(JL)
		    ZQSAT=MIN(0.5,ZQSAT)
                    ZCOR=1./(1.-RETV  *ZQSAT)
                    ZQSAT=ZQSAT*ZCOR
                    ZCOND(JL)=(PQ(JL,KK)-ZQSAT) /(1.+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
                    ZCOND(JL)=MIN(ZCOND(JL),0.)
                    PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND(JL)
                    PQ(JL,KK)=PQ(JL,KK)-ZCOND(JL)
                    IF(ZCOND(JL).NE.0.0) ISUM=ISUM+1
                  ELSE
                    ZCOND(JL)=0.0
                  ENDIF
                ENDDO
      
                IF(ISUM.EQ.0) GO TO 330
 
                DO JL=KIDIA,KFDIA
                 IF(LDFLAG(JL).AND.ZCOND(JL).NE.0.) THEN
                    ZQSAT=FOEEWM(PT(JL,KK))*ZQP(JL)
                    ZQSAT=MIN(0.5,ZQSAT)
                    ZCOR=1./(1.-RETV  *ZQSAT)
                    ZQSAT=ZQSAT*ZCOR
                    ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
                    PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND1
                    PQ(JL,KK)=PQ(JL,KK)-ZCOND1
                  ENDIF
                ENDDO
      
        330     CONTINUE
      
             ENDIF
      
              IF(KCALL.EQ.0) THEN
                DO JL=KIDIA,KFDIA
                  ZQP(JL)=1./PSP(JL)
                  ZQSAT=FOEEWM(PT(JL,KK))*ZQP(JL)
                  ZQSAT=MIN(0.5,ZQSAT)
		  ZCOR=1./(1.-RETV  *ZQSAT)
                  ZQSAT=ZQSAT*ZCOR
                  ZCOND1=(PQ(JL,KK)-ZQSAT) /(1.+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
		  ZCOND1=max(0.,ZCOND1)
		  PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND1
		  PQ(JL,KK)=PQ(JL,KK)-ZCOND1
       	        ENDDO
      
               DO JL=KIDIA,KFDIA
                  ZQSAT=FOEEWM(PT(JL,KK))*ZQP(JL)
                  ZQSAT=MIN(0.5,ZQSAT)
                  ZCOR=1./(1.-RETV  *ZQSAT)
                  ZQSAT=ZQSAT*ZCOR
                  ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
                  PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND1
                  PQ(JL,KK)=PQ(JL,KK)-ZCOND1
	       ENDDO
      
              ENDIF
              IF(KCALL.EQ.4) THEN
                DO JL=KIDIA,KFDIA
                  ZQP(JL)=1./PSP(JL)
                  ZQSAT=FOEEWM(PT(JL,KK))*ZQP(JL)
                  ZQSAT=MIN(0.5,ZQSAT)
                  ZCOR=1./(1.-RETV  *ZQSAT)
                  ZQSAT=ZQSAT*ZCOR
                  ZCOND(JL)=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
                  PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND(JL)
                  PQ(JL,KK)=PQ(JL,KK)-ZCOND(JL)
                ENDDO
 
               DO JL=KIDIA,KFDIA
                  ZQSAT=FOEEWM(PT(JL,KK))*ZQP(JL)
                  ZQSAT=MIN(0.5,ZQSAT)
                  ZCOR=1./(1.-RETV  *ZQSAT)
                  ZQSAT=ZQSAT*ZCOR
                  ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*FOEDEM(PT(JL,KK)))
                  PT(JL,KK)=PT(JL,KK)+FOELDCPM(PT(JL,KK))*ZCOND1
                  PQ(JL,KK)=PQ(JL,KK)-ZCOND1
                ENDDO
  
              ENDIF
      
      !*********************************************
            ELSE ! line 103
      !*********************************************                 
      !     2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
      !                  -----------------------------------------------------
 
              IF (KCALL.EQ.1 ) THEN
      
                ISUM=0
                DO JL=KIDIA,KFDIA
                  IF(LDFLAG(JL)) THEN
                    ZQP(JL)=1./PSP(JL)
                    ZTARG=PT(JL,KK)
                    ZOEALFA=0.5*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.)
	            ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
                    ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
                    ZQSAT=ZQP(JL)*(ZOEALFA*ZFOEEWL+(1.-ZOEALFA)*ZFOEEWI)
                    Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
                    ZQSAT=0.5*((1.-Z1S)*ZQSAT+(1.+Z1S)*ZQMAX) 
                    ZCOR=1./(1.-RETV  *ZQSAT)
                    ZQSAT=ZQSAT*ZCOR
                    Z2S=    ZOEALFA *R5ALVCP*(1./(ZTARG-R4LES)**2)+(1.-ZOEALFA)*R5ALSCP*(1./(ZTARG-R4IES)**2)
                    ZCOND(JL)=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*Z2S)
                    ZCOND(JL)=MAX(ZCOND(JL),0.)
                    PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.-ZOEALFA)*RALSDCP)*ZCOND(JL)
                    PQ(JL,KK)=PQ(JL,KK)-ZCOND(JL)
                    IF(ZCOND(JL).NE.0.0) ISUM=ISUM+1
                  ELSE
                    ZCOND(JL)=0.0
                  ENDIF
               ENDDO
      
                IF(ISUM.EQ.0) GO TO 231
      

                DO JL=KIDIA,KFDIA
                  IF(LDFLAG(JL).AND.ZCOND(JL).NE.0.) THEN
      
                    ZTARG=PT(JL,KK)
                    ZOEALFA=0.5*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.)
                    ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
                    ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
                    ZQSAT=ZQP(JL)*(ZOEALFA*ZFOEEWL+(1.-ZOEALFA)*ZFOEEWI)
                    Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
                    ZQSAT=0.5*((1.-Z1S)*ZQSAT+(1.+Z1S)*ZQMAX) 
                    ZCOR=1./(1.-RETV  *ZQSAT)
                    ZQSAT=ZQSAT*ZCOR
                    Z2S=    ZOEALFA *R5ALVCP*(1./(ZTARG-R4LES)**2)+ (1.-ZOEALFA)*R5ALSCP*(1./(ZTARG-R4IES)**2)
                    ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*Z2S)
                    PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.-ZOEALFA)*RALSDCP)*ZCOND1
                    PQ(JL,KK)=PQ(JL,KK)-ZCOND1
                  ENDIF
               ENDDO
      
        231     CONTINUE
      
              ENDIF
  
              IF(KCALL.EQ.2) THEN
      
                ISUM=0
                 DO JL=KIDIA,KFDIA
                  IF(LDFLAG(JL)) THEN
                    ZQP(JL)=1./PSP(JL)
                    ZTARG=PT(JL,KK)
                    ZOEALFA=0.5*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.)
                    ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
                    ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
                    ZQSAT=ZQP(JL)*(ZOEALFA*ZFOEEWL+(1.-ZOEALFA)*ZFOEEWI)
                    Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
                    ZQSAT=0.5*((1.-Z1S)*ZQSAT+(1.+Z1S)*ZQMAX) 
                    ZCOR=1./(1.-RETV  *ZQSAT)
                    ZQSAT=ZQSAT*ZCOR
                    Z2S=    ZOEALFA *R5ALVCP*(1./(ZTARG-R4LES)**2)+(1.-ZOEALFA)*R5ALSCP*(1./(ZTARG-R4IES)**2)
                    ZCOND(JL)=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*Z2S)
                    ZCOND(JL)=MIN(ZCOND(JL),0.)
                    PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.-ZOEALFA)*RALSDCP)*ZCOND(JL)
                    PQ(JL,KK)=PQ(JL,KK)-ZCOND(JL)
                    IF(ZCOND(JL).NE.0.0) ISUM=ISUM+1
                  ELSE
                    ZCOND(JL)=0.0
                  ENDIF
                ENDDO
      
                IF(ISUM.EQ.0) GO TO 331
      
                 DO JL=KIDIA,KFDIA
                  IF(LDFLAG(JL).AND.ZCOND(JL).NE.0.) THEN
                    ZTARG=PT(JL,KK)
                    ZOEALFA=0.5*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.)
                    ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
                    ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
                    ZQSAT=ZQP(JL)*(ZOEALFA*ZFOEEWL+(1.-ZOEALFA)*ZFOEEWI)
                    Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
                    ZQSAT=0.5*((1.-Z1S)*ZQSAT+(1.+Z1S)*ZQMAX) 
                    ZCOR=1./(1.-RETV  *ZQSAT)
                    ZQSAT=ZQSAT*ZCOR
                    Z2S=    ZOEALFA *R5ALVCP*(1./(ZTARG-R4LES)**2)+(1.-ZOEALFA)*R5ALSCP*(1./(ZTARG-R4IES)**2)
                    ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*Z2S)
                    PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.-ZOEALFA)*RALSDCP)*ZCOND1
                    PQ(JL,KK)=PQ(JL,KK)-ZCOND1
                 ENDIF
                ENDDO
      
        331     CONTINUE
      
              ENDIF
      
              IF(KCALL.EQ.0) THEN
      
                DO JL=KIDIA,KFDIA
                  ZQP(JL)=1./PSP(JL)
                  ZTARG=PT(JL,KK)
                  ZOEALFA=0.5*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.)
	          ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
                  ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
                  ZQSAT=ZQP(JL)*(ZOEALFA*ZFOEEWL+(1.-ZOEALFA)*ZFOEEWI)
		  Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
                  ZQSAT=0.5*((1.-Z1S)*ZQSAT+(1.+Z1S)*ZQMAX) 
		  ZCOR=1./(1.-RETV  *ZQSAT)
                  ZQSAT=ZQSAT*ZCOR
		  Z2S=    ZOEALFA *R5ALVCP*(1./(ZTARG-R4LES)**2)+(1.-ZOEALFA)*R5ALSCP*(1./(ZTARG-R4IES)**2)
                  ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*Z2S)
                  PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.-ZOEALFA)*RALSDCP)*ZCOND1
                  PQ(JL,KK)=PQ(JL,KK)-ZCOND1
                  PT(JL,KK)=PT(JL,KK)
                ENDDO
  

               DO JL=KIDIA,KFDIA
                  ZTARG=PT(JL,KK)
                  ZOEALFA=0.5*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.)
                  ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
                  ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
                  ZQSAT=ZQP(JL)*(ZOEALFA*ZFOEEWL+(1.-ZOEALFA)*ZFOEEWI)
                  Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
                  ZQSAT=0.5*((1.-Z1S)*ZQSAT+(1.+Z1S)*ZQMAX) 
                  ZCOR=1./(1.-RETV  *ZQSAT)
                  ZQSAT=ZQSAT*ZCOR
                  Z2S=    ZOEALFA *R5ALVCP*(1./(ZTARG-R4LES)**2)+(1.-ZOEALFA)*R5ALSCP*(1./(ZTARG-R4IES)**2)
                  ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*Z2S)
                  PT(JL,KK)=PT(JL,KK)+ (ZOEALFA*RALVDCP+(1.-ZOEALFA)*RALSDCP)*ZCOND1
                  PQ(JL,KK)=PQ(JL,KK)-ZCOND1
	       ENDDO
      
      
              ENDIF
      
            IF(KCALL.EQ.4) THEN
      
           DO JL=KIDIA,KFDIA
                  ZQP(JL)=1./PSP(JL)
                  ZTARG=PT(JL,KK)
                  ZOEALFA=0.5*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.)
	          ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
                  ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
                  ZQSAT=ZQP(JL)*(ZOEALFA*ZFOEEWL+(1.-ZOEALFA)*ZFOEEWI)
                  Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
                  ZQSAT=0.5*((1.-Z1S)*ZQSAT+(1.+Z1S)*ZQMAX) 
                  ZCOR=1./(1.-RETV  *ZQSAT)
                  ZQSAT=ZQSAT*ZCOR
                  Z2S=    ZOEALFA *R5ALVCP*(1./(ZTARG-R4LES)**2)+(1.-ZOEALFA)*R5ALSCP*(1./(ZTARG-R4IES)**2)
                  ZCOND(JL)=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*Z2S)
                  PQ(JL,KK)=PQ(JL,KK)-ZCOND(JL)
                ENDDO

                DO JL=KIDIA,KFDIA
                  ZTARG=PT(JL,KK)
                  ZOEALFA=0.5*(TANH(RLPAL1*(ZTARG-RLPTRC))+1.)
                  ZFOEEWL=R2ES*EXP(R3LES*(ZTARG-RTT)/(ZTARG-R4LES))
                  ZFOEEWI=R2ES*EXP(R3IES*(ZTARG-RTT)/(ZTARG-R4IES))
                  ZQSAT=ZQP(JL)*(ZOEALFA*ZFOEEWL+(1.-ZOEALFA)*ZFOEEWI)
                  Z1S=TANH(RLPAL2*(ZQSAT-ZQMAX))
                  ZQSAT=0.5*((1.-Z1S)*ZQSAT+(1.+Z1S)*ZQMAX) 
                  ZQSAT=MIN(ZQMAX,ZQSAT)
                  ZCOR=1./(1.-RETV  *ZQSAT)
                  ZQSAT=ZQSAT*ZCOR
                  Z2S=    ZOEALFA *R5ALVCP*(1./(ZTARG-R4LES)**2)+(1.-ZOEALFA)*R5ALSCP*(1./(ZTARG-R4IES)**2)
                  ZCOND1=(PQ(JL,KK)-ZQSAT)/(1.+ZQSAT*ZCOR*Z2S)
                  PT(JL,KK)=PT(JL,KK)+(ZOEALFA*RALVDCP+(1.-ZOEALFA)*RALSDCP)*ZCOND1
                  PQ(JL,KK)=PQ(JL,KK)-ZCOND1
                ENDDO
             ENDIF
      
!*********************************************
            ENDIF
!*********************************************          
      
  
 END SUBROUTINE CUADJTQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 REAL FUNCTION FOEALFA (PTARG) 
 REAL  PTARG
 REAL, PARAMETER :: RTICE = 273.16-0.1
 REAL, PARAMETER :: RTWAT = 273.16 
 FOEALFA = MIN(1.,((MAX(RTICE,MIN(RTWAT,PTARG))-RTICE)/(RTWAT-RTICE))**2) 
 RETURN
 END FUNCTION FOEALFA


 REAL FUNCTION FOEDELTA (PTARG) 
 REAL  PTARG
 REAL, PARAMETER :: RTT = 273.16
 FOEDELTA = MAX (0.,SIGN(1.,PTARG-RTT))
 RETURN
 END FUNCTION FOEDELTA

 
 REAL FUNCTION FOEEWM (PTARG)
 REAL  PTARG
 REAL, PARAMETER :: RTT     = 273.16
 REAL, PARAMETER :: RESTT   = 611.14
 REAL, PARAMETER :: RD      = 287.04
 REAL, PARAMETER :: RV      = 461.50
 REAL, PARAMETER :: R2ES    = RESTT*(RD/RV)
 REAL, PARAMETER :: R3LES   = 17.269 
 REAL, PARAMETER :: R3IES   = 21.875
 REAL, PARAMETER :: R4LES   = 35.86 
 REAL, PARAMETER :: R4IES   = 7.66
 FOEEWM =(R2ES*(FOEALFA(PTARG)*EXP(R3LES*(PTARG-RTT)/(PTARG-R4LES))+(1.-FOEALFA(PTARG))*EXP(R3IES*(PTARG-RTT)/(PTARG-R4IES))))
 RETURN
 END FUNCTION FOEEWM


 REAL FUNCTION FOEALFCU (PTARG) 
 REAL  PTARG
 REAL, PARAMETER :: RTWAT = 273.16 
 REAL, PARAMETER :: RTT = 273.16
 REAL, PARAMETER :: RTICECU = RTT-23.0
 REAL, PARAMETER :: RTWAT_RTICECU_R  = 1./(RTWAT-RTICECU)
 FOEALFCU = MIN(1.,((MAX(RTICECU,MIN(RTWAT,PTARG))-RTICECU)*RTWAT_RTICECU_R)**2)
 RETURN
 END FUNCTION FOEALFCU


 REAL FUNCTION FOEDEM (PTARG)
 REAL  PTARG
 REAL, PARAMETER :: R3LES   = 17.269 
 REAL, PARAMETER :: R3IES   = 21.875
 REAL, PARAMETER :: R4LES   = 35.86 
 REAL, PARAMETER :: R4IES   = 7.66
 REAL, PARAMETER :: R5LES   = R3LES*(273.16-R4LES)
 REAL, PARAMETER :: R5IES   = R3IES*(273.16-R4IES)
 REAL, PARAMETER :: RD      = 287.04
 REAL, PARAMETER :: RV      = 461.50
 REAL, PARAMETER :: RCPD    = 3.5*RD
 REAL, PARAMETER :: RLVTT   = 2.5008E+6
 REAL, PARAMETER :: RLSTT   = 2.8345E+6
 REAL, PARAMETER :: R5ALVCP = R5LES*RLVTT/RCPD
 REAL, PARAMETER :: R5ALSCP = R5IES*RLSTT/RCPD 
 REAL, PARAMETER :: RALVDCP = RLVTT/RCPD
 FOEDEM = FOEALFA(PTARG)*R5ALVCP*(1./(PTARG-R4LES)**2)+(1.-FOEALFA(PTARG))*R5ALSCP*(1./(PTARG-R4IES)**2)
 RETURN
 END FUNCTION FOEDEM


 REAL FUNCTION FOELDCPM(PTARG)
 REAL  PTARG
 REAL, PARAMETER :: RD      = 287.04
 REAL, PARAMETER :: RCPD    = 3.5*RD
 REAL, PARAMETER :: RLSTT   = 2.8345E+6
 REAL, PARAMETER :: RLVTT   = 2.5008E+6
 REAL, PARAMETER :: RALSDCP = RLSTT/RCPD
 REAL, PARAMETER :: RALVDCP = RLVTT/RCPD
 FOELDCPM = FOEALFA(PTARG)*RALVDCP+(1.-FOEALFA(PTARG))*RALSDCP
 RETURN
 END FUNCTION FOELDCPM

 END MODULE kbcon_ecmwf
 
 
 
