!############################# Change Log ##################################
! 5.0.2
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003, 2006 - All Rights Reserved
!  Brazilian Regional Atmospheric Modeling System - BRAMS
!###########################################################################

!Subroutine to link leaf3 to TEB for landclass = urban type 1 and 
!urban type 2
!Adapted by Edmilson Freitas 
!DATE : Jul 7th 2006
!Last modification:

! TEB
SUBROUTINE LEAF3_TEB_INTERFACE(             &
     ISTP, ZTSTEPFRC, ZTSTEP, COSZ, ZZREF,  &
     ZRAT, SOLAR,                           &
     ZPA, ZTA, ZUA, ZVA, RV,                &
     ZPCP_IN, fuso,                         &
     ZT_CANYON, ZR_CANYON,                  &
     ZTS_ROOF, ZTS_ROAD, ZTS_WALL,          &
     ZTI_ROAD, ZTI_BLD,                     &
     ZWS_ROOF, ZWS_ROAD,                    &
     ZT_ROOF, ZT_ROAD, ZT_WALL,             &
     ZH_TOWN, ZLE_TOWN, ZEMIS_TOWN,         &
     ZSFU_TOWN, ZSFV_TOWN, ZTS_TOWN,        &
     ZALB_TOWN, IG_URBAN,                   &
     ZH_TRAFFIC, ZH_INDUSTRY,               &
     ZLE_TRAFFIC, ZLE_INDUSTRY,             &
     T2M, R2M,                              &
     time, itime1, dpdz, dens               )

  USE teb_vars_const
  USE mem_emiss, ONLY: EFSAT, EFSUN, WEEKDAYIN ! INTENT(IN)

  IMPLICIT NONE
  ! constants:
  INTEGER, PARAMETER  :: INTEB = 3    ! TEB Vertical dimension
  ! Arguments:
  INTEGER, INTENT(in) :: ISTP ! Not used
  REAL, INTENT(in)    :: ZTSTEPFRC !time step of the input forcing ! Not used
  REAL, INTENT(in)    :: ZTSTEP    ! Time step of integration (sec.)
  REAL, INTENT(in)    :: COSZ      ! Cos(Zenith Angle)
  REAL, INTENT(in)    :: ZZREF     ! reference height of the first atm. level
  REAL, INTENT(in)    :: ZRAT      ! atmospheric infrared radiation
  REAL, INTENT(in)    :: SOLAR
  REAL, INTENT(in)    :: ZPA       ! atmospheric level pressure
  REAL, INTENT(in)    :: ZTA       ! atmospheric temperature at level za
  REAL, INTENT(in)    :: ZUA, ZVA  ! wind speeds in the x and y dirs at level za
  REAL, INTENT(in)    :: RV        ! atmospheric mixing ratio at level za
  REAL, INTENT(in)    :: ZPCP_IN   ! Precip. Rate (kg/m^2/s)
  REAL, INTENT(in)    :: fuso      ! local time correction
  REAL, INTENT(inout) :: ZT_CANYON ! canyon air temperature
  REAL, INTENT(inout) :: ZR_CANYON ! canyon air vapor mixing ratio
  REAL, INTENT(inout) :: ZTS_ROOF  ! roof surface temperature
  REAL, INTENT(inout) :: ZTS_ROAD  ! road surface temperature
  REAL, INTENT(inout) :: ZTS_WALL  ! wall surface temperature
  REAL, INTENT(inout) :: ZTI_ROAD  ! road deep temperature
  REAL, INTENT(inout) :: ZTI_BLD   ! INTERNAL BLDING TEMP (K)
  REAL, INTENT(inout) :: ZWS_ROOF  ! roof water reservoir
  REAL, INTENT(inout) :: ZWS_ROAD  ! road water reservoir
  REAL, INTENT(inout) :: ZT_ROOF(INTEB) ! roof layers temperatures
  REAL, INTENT(inout) :: ZT_ROAD(INTEB) ! road layers temperatures
  REAL, INTENT(inout) :: ZT_WALL(INTEB) ! wall layers temperatures
  REAL, INTENT(out)   :: ZH_TOWN        ! TOWN AVE. SENS. HEAT FLUX
  REAL, INTENT(out)   :: ZLE_TOWN       !TOWN AVE. LAT. HEAT FLUX
  REAL, INTENT(out)   :: ZEMIS_TOWN     ! town equivalent emissivity
  REAL, INTENT(out)   :: ZSFU_TOWN      !TOWN SCALE EDDY U-MOM FLUX
  REAL, INTENT(out)   :: ZSFV_TOWN      !TOWN SCALE EDDY V-MOM FLUX
  REAL, INTENT(out)   :: ZTS_TOWN       !TOWN SFC TEMP
  REAL, INTENT(out)   :: ZALB_TOWN      !TOWN EQV. ALBEDO
  INTEGER, INTENT(inout) :: IG_URBAN    ! Flag for urban(1) and suburban(2)
  REAL, INTENT(out)   :: ZH_TRAFFIC
  REAL, INTENT(out)   :: ZH_INDUSTRY
  REAL, INTENT(out)   :: ZLE_TRAFFIC
  REAL, INTENT(out)   :: ZLE_INDUSTRY
  REAL, INTENT(out)   :: T2M            ! Extrapolated 2 m temperature
  REAL, INTENT(out)   :: R2M            ! Extrapolated 2 m specif humidity
  INTEGER, INTENT(in) :: itime1
  REAL, INTENT(in)    :: time
  REAL, INTENT(in)    :: dpdz
  REAL, INTENT(in)    :: dens
  ! Local Variables:
  ! Declarations of variables
  !INPUT / OUTPUT VARIABLES
  REAL :: pfat
  
  REAL :: ZRG        ! incoming solar radiation
  REAL :: ZPS        ! surface pressure
  REAL :: ZQA        ! atmospheric specific humidity at level za
  REAL :: ZRR        ! rain rate
!!$  real :: ZPRECIP    ! precipitation rate ! not used
  !    0.9.D   TEB Diagnostics:
  !           ----------------
  REAL :: ZRN_ROOF     ! net radiation over roof
  REAL :: ZH_ROOF      ! sensible heat flux over roof
  REAL :: ZLE_ROOF     ! latent heat flux over roof
  REAL :: ZGFLUX_ROOF  ! flux through the roof
  REAL :: ZRUNOFF_ROOF ! runoff over the ground
  REAL :: ZRN_ROAD     ! net radiation over road
  REAL :: ZH_ROAD      ! sensible heat flux over road
  REAL :: ZLE_ROAD     ! latent heat flux over road
  REAL :: ZGFLUX_ROAD  ! flux through the road
  REAL :: ZRUNOFF_ROAD ! runoff over the ground
  REAL :: ZRN_WALL     ! net radiation over wall
  REAL :: ZH_WALL      ! sensible heat flux over wall
  REAL :: ZLE_WALL     ! latent heat flux over wall
  REAL :: ZGFLUX_WALL  ! flux through the wall
  !
  !only OUTPUT VARIABLES
  REAL :: ZU_CANYON    ! CANYON HOR. WIND
  REAL :: ZRN_TOWN     ! TOWN SCALE NET RAD
  REAL :: ZCH_TOWN     ! TOWN AVERAGED HEAT TRANSFER
  REAL :: ZGFLUX_TOWN  ! TOWN SCALE GROUND HEAT STORAGE
  REAL :: ZRUNOFF_TOWN ! TOWN SCALE RUNOFF
  !  MODEL FORCING VARIABLES
  REAL :: ZSVF_ROAD     ! ROAD SKY VIEW FACTOR
  REAL :: ZSVF_WALL     ! WALL SKY VIEW FACTOR
  REAL :: ZCAN_HW_RATIO ! CANYON HEIGHT TO WIDTH RATIO
  REAL :: ZDIR_SW_RAD   ! incoming direct solar radiation on an horiz.surface
  REAL :: ZSCA_SW_RAD   ! scattered incoming solar rad.

  INTEGER :: I  !!, JDT, JINDT
  REAL ::       ZZ0_TOWN	    &
       ,	ZBLD		    &
       ,	ZBLD_HEIGHT	    &
       ,	ZBLD_HL_RATIO	    &
       ,	ZALB_ROOF	    &
       ,	ZEMIS_ROOF	    &
       ,	ZALB_ROAD	    &
       ,	ZEMIS_ROAD	    &
       ,	ZALB_WALL	    &
       ,	ZEMIS_WALL	    &
       ,        ax1, ax2, tmp_teb, bx1, bx2 &
       ,        timeq1, timeq2, tign
  INTEGER :: idays
  CHARACTER(len=3) :: cday
  REAL :: ZDIRCOSZW, ZTANZEN
  !                  ZDIRCOSZW = Cosinus of the angle between the 
  !                  normal to the surface and the vertical
  !                  ZTANZEN   = tangent of solar zenith angle
  !
  REAL   :: ZRHOA,  ZEXNA, ZEXNS, ZVMOD, ZTVI, ZAZENIT, ZEXN2, ZP2
  !                  ZRHOA   = air density
  !                  ZEXNA   = Exner function
  !                  ZEXNS   = Exner function
  !                  ZVMOD   = modulus of the wind parallel to the orography
  !                  ZTVI    = virtual temperature 
  !                  ZAZENIT = solar zenith angle
  REAL, EXTERNAL ::  rslf ! Function in "therm_lib.f90"

  
  ZRG = SOLAR
  
  ZDIRCOSZW = 1.0
  
  ZQA = (1./((1000./rv) + 1.))*1000.
  
  
  ZZ0_TOWN        =  Z0_TOWN(IG_URBAN)    
  ZBLD            =  BLD (IG_URBAN) 
  ZBLD_HEIGHT     =  BLD_HEIGHT (IG_URBAN)    
  ZBLD_HL_RATIO   =  BLD_HL_RATIO(IG_URBAN)
  ZALB_ROOF       =  AROOF(IG_URBAN)
  ZEMIS_ROOF      =  EROOF(IG_URBAN)
  ZALB_ROAD       =  AROAD(IG_URBAN)
  ZEMIS_ROAD      =  EROAD(IG_URBAN)
  ZALB_WALL       =  AWALL(IG_URBAN)
  ZEMIS_WALL      =  EWALL(IG_URBAN)
  
  ZH_TRAFFIC      =  HTRAF(IG_URBAN)
  ZH_INDUSTRY     =  HINDU(IG_URBAN)	
  ZLE_TRAFFIC     =  PLETRAF(IG_URBAN)	 
  ZLE_INDUSTRY    =  PLEINDU(IG_URBAN)
  
  !
  !*      1.     Calculate the canyon shape
  !              --------------------------
  !
  ZCAN_HW_RATIO = ZBLD_HL_RATIO*ZBLD/(1.0 - ZBLD)
  !write(*,*)'ZCAN_HW_RATIO=',ZCAN_HW_RATIO
  !write(*,*)'ZBLD_HL_RATIO=',ZBLD_HL_RATIO
  !pause				
  !!*      2.     Calculate sky view factors
  !              --------------------------
  !
  ZSVF_ROAD  = SQRT(ZCAN_HW_RATIO*ZCAN_HW_RATIO + 1.0) - ZCAN_HW_RATIO
  ZSVF_WALL  = 0.5*(1.0 - ZSVF_ROAD)/ZCAN_HW_RATIO
  !*      3.     Calculate town albedo and emissivity
  !              ------------------------------------
  !
  ZALB_TOWN  = (1.0 - ZBLD)*ZALB_ROAD + ZBLD*ZALB_ROOF
  ZEMIS_TOWN = (1.0 - ZBLD)*ZEMIS_ROAD*ZSVF_ROAD + &
       (1.0 - ZBLD)*ZEMIS_WALL*(1.0 - ZSVF_ROAD) + ZBLD*ZEMIS_ROOF
  ZRR       = ZPCP_IN
  
  IF (ZRG < 0.0) ZRG = 0.0

  ! Initialize forcing values:
  ZAZENIT = ACOS(COSZ)
  ZTANZEN = TAN(ZAZENIT)
  IF (ZRG>0.0 .AND. ZTANZEN<0.0) THEN
     ZTANZEN = SQRT(9999.0)
  ENDIF
  
  ! Perform some conversions and final calculations using forcing variables:
  !
  ZVMOD   = SQRT(ZUA*ZUA + ZVA*ZVA)
  !
  ! Make sure wind magnitude is above a minimum threshold (original =1.0, 
  ! for test purpose = 0.5) :
  !
  ZVMOD   = MAX(0.5, ZVMOD)
  ZTVI    = ZTA*(1. + ((XRV/XRD) - 1.)*ZQA )
  ZRHOA   = ZPA/XRD/ZTVI
  !            ZRHOA   = dens
  ZEXNA   = (ZPA/100000.)**(XRD/XCPD)
  
  ZPS     = DPDZ*(ZBLD_HEIGHT-ZZREF) + ZPA
  ZEXNS   = (ZPS/100000.)**(XRD/XCPD)
  
  ! Calculate forcing needed by TEB from existing forcing
  ! for now assume all incoming radiation is direct
  !
  ZDIR_SW_RAD = ZRG
  ZSCA_SW_RAD = 0.0
  !
  
  tmp_teb = time + (itime1/100 + MOD(itime1,100)/60.)*3600
  idays   = INT((tmp_teb/3600.)/24.)  !number of days of simulation
  tign    = REAL(idays)*24.*3600.
  
  pfat=1.
  
  CALL EMFACTOR(WEEKDAYIN, idays, cday)
  IF (cday=='SAT') pfat = EFSAT
  IF (cday=='SUN') pfat = EFSUN
  
  !Fonte urbana (castanho, 1999 - tese de mestrado)
  bx1    = RUSHH1 - fuso + DAYLIGHT
  bx2    = RUSHH2 - fuso + DAYLIGHT
  timeq1 = (tmp_teb - tign)/3600. - bx1
  timeq2 = (tmp_teb - tign)/3600. - bx2
    
  !fonte de calor sensivel veicular
  ax2        = ZH_TRAFFIC - 5.
  ax1        = 0.63*ax2
  ZH_TRAFFIC = ((ax1*EXP(-(timeq1)**2/8.5) + &
       ax2*EXP(-(timeq2)**2/10.6))*pfat + 5.)
  
  !fonte de calor latente veicular
  ax2         = ZLE_TRAFFIC - 5.
  ax1         = 0.63*ax2
  ZLE_TRAFFIC = ((ax1*EXP(-(timeq1)**2/8.5) + &
       ax2*EXP(-(timeq2)**2/10.6))*pfat + 5.)
  
  !                              -------------------------------
  ZRN_ROOF     = 0.
  ZH_ROOF      = 0.
  ZLE_ROOF     = 0.
  ZGFLUX_ROOF  = 0.
  ZRUNOFF_ROOF = 0.
  ZRN_ROAD     = 0.
  ZH_ROAD      = 0.
  ZLE_ROAD     = 0.
  ZGFLUX_ROAD  = 0.
  ZRUNOFF_ROAD = 0.
  ZRN_WALL     = 0.
  ZH_WALL      = 0
  ZLE_WALL     = 0.
  ZGFLUX_WALL  = 0.
  ZRN_TOWN     = 0.
  ZH_TOWN      = 0.
  ZLE_TOWN     = 0.
  ZGFLUX_TOWN  = 0.
  ZRUNOFF_TOWN = 0.
  ZSFU_TOWN    = 0.
  ZSFV_TOWN    = 0.
  ZCH_TOWN     = 0.
  
  !*                        14.B   CALL THE MAIN SUBROUTINE OF TEB
  CALL URBAN(ZTS_TOWN, ZEMIS_TOWN, ZALB_TOWN,                   &
       ZT_CANYON, ZR_CANYON, ZU_CANYON,                         &
       ZTS_ROOF,ZTS_ROAD,ZTS_WALL,ZTI_ROAD,ZTI_BLD,             &
       ZT_ROOF, ZT_ROAD, ZT_WALL, ZWS_ROOF,ZWS_ROAD,            &
       ZPS, ZPA, ZEXNS, ZEXNA, ZTA, ZQA, ZRHOA,                 &
       ZRAT, ZDIR_SW_RAD, ZSCA_SW_RAD, ZTANZEN,                 &
       ZRR, ZZREF, ZDIRCOSZW, ZUA, ZVA, ZVMOD,                  &
       ZH_TRAFFIC, ZLE_TRAFFIC, ZH_INDUSTRY, ZLE_INDUSTRY,      &
       ZTSTEP,                                                  &
       ZZ0_TOWN,                                                &
       ZBLD, ZBLD_HEIGHT, (2.*ZBLD_HL_RATIO*ZBLD),              &
       ZCAN_HW_RATIO, ZALB_ROOF, ZEMIS_ROOF,                    &
       HC_ROOF(1:3), TC_ROOF(1:3), D_ROOF(1:3),                 &
       ZALB_ROAD, ZEMIS_ROAD, ZSVF_ROAD,                        &
       HC_ROAD(1:3), TC_ROAD(1:3), D_ROAD(1:3),                 &
       ZALB_WALL, ZEMIS_WALL, ZSVF_WALL,                        &
       HC_WALL(1:3), TC_WALL(1:3), D_WALL(1:3),                 &
       ZRN_ROOF, ZH_ROOF, ZLE_ROOF, ZGFLUX_ROOF,                &
       ZRUNOFF_ROOF,                                            &
       ZRN_ROAD, ZH_ROAD, ZLE_ROAD, ZGFLUX_ROAD,                &
       ZRUNOFF_ROAD,                                            &
       ZRN_WALL, ZH_WALL, ZLE_WALL, ZGFLUX_WALL,                &
       ZRN_TOWN, ZH_TOWN, ZLE_TOWN, ZGFLUX_TOWN,                &
       ZRUNOFF_TOWN, ZSFU_TOWN, ZSFV_TOWN, ZCH_TOWN)

  !Extrapolating temperature and specific humidity to 2 m height
  
  ZP2   = DPDZ*(2. - ZBLD_HEIGHT) + ZPS
  
  ZEXN2 = (ZP2/100000.)**(XRD/XCPD)
  
  T2M   = ZT_WALL(1)*ZEXN2/ZEXNS
  
  R2M   = ZR_CANYON*(rslf(ZP2,T2M)/rslf(ZPS, ZT_WALL(1)))
  
END SUBROUTINE LEAF3_TEB_INTERFACE

! TEB
!     ############################################
SUBROUTINE INI_TG_PROFILE(PTS, PTI, PTC, PD, PT)
  !   ############################################
  !
  !!****  *INI_TG_PROFILE*  
  !!
  !!    PURPOSE
  !!    -------
  !
  !     Computes the equilibrium of a temperature profile through a
  !     conductive material, given the two extreme surface temperatures.
  !         
  !!      
  !!    AUTHOR
  !!    ------
  !!
  !!	V. Masson           * Meteo-France *
  !!
  !!    MODIFICATIONS
  !!    -------------
  !!      Original 02/11/98 
  !!               17/05/2002 Edmilson Freitas (IAG-USP). Elimination of the
  !!                          the declarations that are not needed. Change from
  !!                          Module to normal subroutine.
  !----------------------------------------------------------------------------
  !
  !*       0.     DECLARATIONS
  !               ------------
  !
  !
  IMPLICIT NONE
  !
  !*      0.1    declarations of arguments
  !
  !INPUT VARIABLES:
  REAL, INTENT(in)  :: PTS      ! surface temperature
  REAL, INTENT(in)  :: PTI      ! internal temperature
  REAL, INTENT(in)  :: PTC(3)   ! thermal conductivity for roof layers
  REAL, INTENT(in)  :: PD(3)    ! depth of roof layers
  !
  !OUTPUT VARIABLES
  REAL, INTENT(out) :: PT(3)    ! layers temperatures
  !
  !*      0.2    declarations of local variables
  !
  REAL, DIMENSION(SIZE(PT))   :: ZA ! lower diag.
  REAL, DIMENSION(SIZE(PT))   :: ZB ! main  diag.
  REAL, DIMENSION(SIZE(PT))   :: ZC ! upper diag.
  REAL, DIMENSION(SIZE(PT))   :: ZY ! r.h.s.
  REAL, DIMENSION(SIZE(PT))   :: ZX ! solution
  !
  REAL, DIMENSION(0:SIZE(PT)) :: ZMTC_O_D
  ! mean thermal conductivity over distance between 2 layers
  !
  INTEGER                     :: ILAYER ! number of roof,road or wall layers
  INTEGER                     :: JLAYER ! loop counter
  !----------------------------------------------------------------------------
  !
  !*      1.     Layer thermal properties
  !              ------------------------
  !
  ILAYER      = SIZE(PT)
  ZA(:)       = 0.
  ZB(:)       = 0.
  ZC(:)       = 0.
  ZX(:)       = 0.
  ZY(:)       = 0.
  ZMTC_O_D(:) = 0.
  !
  ZMTC_O_D(0) = 2.*PTC(1)/PD(1)
  !
  DO JLAYER=1,ILAYER-1
     ZMTC_O_D(JLAYER) = 2./( PD(JLAYER)/PTC(JLAYER) + &
          PD(JLAYER+1)/PTC(JLAYER+1) )
  END DO
  !
  ZMTC_O_D(ILAYER) = 2.*PTC(ILAYER)/PD(ILAYER)
  !
  !----------------------------------------------------------------------------
  !
  !*      2.     Surface layer coefficients
  !              --++++++------------------
  !
  !
  ZA(1) =   0.
  
  ZB(1) =   ZMTC_O_D(0) + ZMTC_O_D(1)
  
  ZC(1) = - ZMTC_O_D(1)
  !
  ZY(1) =   ZMTC_O_D(0)*PTS
  !
  !
  !----------------------------------------------------------------------------
  !
  !*      3.     Other layers coefficients
  !              -------------------------
  !
  DO JLAYER=2,ILAYER-1
     ZA(JLAYER) = - ZMTC_O_D(JLAYER - 1)

     ZB(JLAYER) =   ZMTC_O_D(JLAYER - 1) + ZMTC_O_D(JLAYER)
     
     ZC(JLAYER) = - ZMTC_O_D(JLAYER)

     ZY(JLAYER) =   0.
  END DO
  !
  !----------------------------------------------------------------------------
  !
  !*      4.     Inside layer coefficients
  !              -------------------------
  !
  ZA(ILAYER) = - ZMTC_O_D(ILAYER - 1)

  ZB(ILAYER) =   ZMTC_O_D(ILAYER - 1) + ZMTC_O_D(ILAYER)
  
  ZC(ILAYER) =   0.

  ZY(ILAYER) =   ZMTC_O_D(ILAYER)*PTI
  !
  !
  !----------------------------------------------------------------------------
  !
  !*      5.     Tri-diagonal system resolution
  !              ------------------------------
  !
  CALL TRID(ZX, ZA, ZB, ZC, ZY, ILAYER)
  !
  PT(:) = ZX(:)
  !
  !----------------------------------------------------------------------------

END SUBROUTINE INI_TG_PROFILE

!-------------------------------------------------------------------------------
! TEB
SUBROUTINE TEB_INIT(n1, n2, n3, np, vegt, theta, rv, pi, pp,       &
     TROOF, TROAD, TWALL, TIBLD, TIROAD, TCANYON, RCANYON, TSROOF, &
     TSROAD, TSWALL, HT, LET, HIN, LEIN, WSROOF, WSROAD, EMISTOWN, &
     ALBTOWN, TSTOWN, G_URBAN)
  
  USE teb_vars_const, ONLY: &
       TMINBLD,             & ! INTENT(IN)
       D_ROAD,              & ! INTENT(IN)
       TC_ROAD,             & ! INTENT(IN)
       D_WALL,              & ! INTENT(IN)
       TC_WALL,             & ! INTENT(IN)
       D_ROOF,              & ! INTENT(IN)
       TC_ROOF,             & ! INTENT(IN)
       NURBTYPE,            & ! INTENT(IN)
       ILEAFCOD               ! INTENT(IN)

  IMPLICIT NONE
  ! Arguments:
  INTEGER, INTENT(IN) :: n1, n2, n3, np
  REAL, INTENT(IN)    :: vegt(n2,n3,np)
  REAL, INTENT(IN)    :: theta(n1,n2,n3)
  REAL, INTENT(IN)    :: rv(n1,n2,n3)
  REAL, INTENT(IN)    :: pi(n1,n2,n3)
  REAL, INTENT(IN)    :: pp(n1,n2,n3)
  REAL, INTENT(OUT)   :: TROOF(n1,n2,n3)
  REAL, INTENT(OUT)   :: TROAD(n1,n2,n3)
  REAL, INTENT(OUT)   :: TWALL(n1,n2,n3)
  REAL, INTENT(OUT)   :: TIBLD(n2,n3)
  REAL, INTENT(OUT)   :: TIROAD(n2,n3)
  REAL, INTENT(OUT)   :: TCANYON(n2,n3)
  REAL, INTENT(OUT)   :: RCANYON(n2,n3)
  REAL, INTENT(OUT)   :: TSROOF(n2,n3)
  REAL, INTENT(OUT)   :: TSROAD(n2,n3)
  REAL, INTENT(OUT)   :: TSWALL(n2,n3)
  REAL, INTENT(OUT)   :: HT(n2,n3)
  REAL, INTENT(OUT)   :: LET(n2,n3)
  REAL, INTENT(OUT)   :: HIN(n2,n3)
  REAL, INTENT(OUT)   :: LEIN(n2,n3)
  REAL, INTENT(OUT)   :: WSROOF(n2,n3)
  REAL, INTENT(OUT)   :: WSROAD(n2,n3)
  REAL, INTENT(INOUT) :: EMISTOWN(n2,n3)
  REAL, INTENT(INOUT) :: ALBTOWN(n2,n3)
  REAL, INTENT(INOUT) :: TSTOWN(n2,n3)
  REAL, INTENT(OUT)   :: G_URBAN(n2,n3,np)

  ! Local Variables:
  INTEGER :: i, j, k, ilf, inp 
  REAL    :: cpi, hcpi, pis, pl2

  ! Initiating data
  HT      = 0.
  LET     = 0.
  HIN     = 0.
  LEIN    = 0.
  WSROOF  = 0.
  WSROAD  = 0.
  TIBLD   = 0.
  TIROAD  = 0.
  TCANYON = 0.
  RCANYON = 0.
  TSROOF  = 0.
  TSROAD  = 0.
  TSWALL  = 0.
  TROOF   = 0.
  TROAD   = 0.
  TWALL   = 0.
  G_URBAN = 0.
  
  cpi = 1./1004.
  hcpi = 0.5*cpi
  
  DO i=1,n2
     
     DO j=1,n3
        
        pis = (pp(1,i,j) + pi(1,i,j) + pp(2,i,j) + pi(2,i,j)) * hcpi
        
        pl2 = (pp(2,i,j) + pi(2,i,j)) * cpi
        
        DO inp=2,np !patch looping
           
           DO ilf=1,NURBTYPE
              IF (NINT(vegt(i,j,inp))==ILEAFCOD(ILF)) &
                   G_URBAN(i,j,inp)=float(ilf)
           ENDDO
           
           IF (NINT(G_URBAN(i,j,inp))/=0.) THEN

              !internal temperature defined in RAMSIN
              TIBLD(i,j)    = 273.16 + TMINBLD
              !surface
              !TIROAD(i,j)   = (theta(2,i,j)+theta(1,i,j))*0.5*pis
              TIROAD(i,j)   = 281.16
              !model's first level
              TSROOF(i,j)   = (theta(2,i,j))*pl2
              !surface
              TSROAD(i,j)   = (theta(2,i,j)+theta(1,i,j))*0.5*pis
              !average surface and first level
              TSWALL(i,j)   = (TSROOF(i,j)+TSROAD(i,j))*0.5
              !average surface and first level
              TCANYON(i,j)  = TSWALL(i,j)
              !average surface and first level
              TSTOWN(i,j)   = TSWALL(i,j)
              EMISTOWN(i,j) =   0.8878604
              ALBTOWN(i,j)  =   0.129
              RCANYON(i,j)  = (1./((1000./rv(2,i,j))+1.))*1000.
              HT(i,j)       =   5.
              LET(i,j)      =   5.
              HIN(i,j)      =  10.
              LEIN(i,j)     =  40.
              
              CALL INI_TG_PROFILE(TSROOF(i,j), TIBLD(i,j), &
                   TC_ROOF(1:3),D_ROOF(1:3), TROOF(2:4,i,j))

              CALL INI_TG_PROFILE(TSROAD(i,j), TIROAD(i,j), &
                   TC_ROAD(1:3),D_ROAD(1:3), TROAD(2:4,i,j))

              CALL INI_TG_PROFILE(TSWALL(i,j), TIBLD(i,j), &
                   TC_WALL(1:3),D_WALL(1:3), TWALL(2:4,i,j))

           ENDIF
           
        ENDDO
        
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE TEB_INIT


! TEB
SUBROUTINE TEBC_INIT(n2, n3, np, G_URBAN, EMIS, ALB, TS)
  
  IMPLICIT NONE
  ! Arguments:
  INTEGER, INTENT(in) :: n2, n3, np
  REAL, INTENT(out)   :: G_URBAN(n2,n3,np)
  REAL, INTENT(out)   :: EMIS(n2,n3)
  REAL, INTENT(out)   :: ALB(n2,n3)
  REAL, INTENT(out)   :: TS(n2,n3)

  G_URBAN = 0.
  EMIS    = 0.
  ALB     = 0.
  TS      = 0.

END SUBROUTINE TEBC_INIT
