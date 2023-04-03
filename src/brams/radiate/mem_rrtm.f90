MODULE mem_rrtm
   use parkind, only : im => kind_im, rb => kind_rb
   USE rconstants  , ONLY : cp,cpor,p00,stefan,cpi
   USE rrtmg_lw_init, ONLY: rrtmg_lw_ini
   USE rrtmg_sw_init, ONLY: rrtmg_sw_ini
   use parrrsw, only : nbndsw
   use parrrtm, only : nbndlw
   
   REAL(kind=rb), ALLOCATABLE, dimension(:) :: sig,delsig,sigmid
   INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
   INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
   REAL(KIND=rb), PARAMETER :: pi=3.1415926e00_rb
   REAL(KIND=rb), PARAMETER :: adjust=0.01_rb !to convert from Pascal to mbar
   INTEGER :: nls
   
   !Pesos moleculares [g/mol] e razoes
   REAL(KIND=rb), PARAMETER :: pmAr =28.96_rb
   REAL(KIND=rb), PARAMETER :: pmO3 =47.998_rb  ,rO3Ar =1/(pmO3/pmAr)
   REAL(KIND=rb), PARAMETER :: pmCO2=44.01_rb   ,rCO2Ar=1/(pmCO2/pmAr)
   REAL(KIND=rb), PARAMETER :: pmCH4=16.04_rb   ,rCh4Ar=1/(pmCH4/pmAr)
   REAL(KIND=rb), PARAMETER :: pmNO2=46.0055_rb ,rNO2Ar=1/(pmNO2/pmAr)
   REAL(KIND=rb), PARAMETER :: pmO2 =64.0_rb    ,rO2AR =1/(pmO2/pmAr)
   REAL(KIND=rb), PARAMETER :: pmH2O=18.01528_rb,rH2OAr=1/(pmH2O/pmAr)   
   INTEGER,ALlocatable,DIMENSION(:) :: flip
   REAL(kind=rb) :: prsnz
   
   REAL(KIND=rb),ALlocatable :: aot_rrtm_lw(:,:,:) !m2,m3,ngptlw
   REAL(KIND=rb),ALlocatable :: aot_rrtm_sw(:,:,:) !m2,m3,ngptsw
   REAL(KIND=rb),ALlocatable :: co3FromTuv(:,:,:)
   
   LOGICAL :: o3isPresent,co2isPresent,ch4isPresent
   LOGICAL :: no2isPresent,o2isPresent,h2oispresent
   INTEGER :: o3pos,co2pos,ch4pos,no2pos,o2pos,h2opos
   
   !Number of band to convert from TUV-Carma to nbndsw of RRTM
   INTEGER, PARAMETER, DIMENSION(14) :: bndmin_sw=(/32,36,39,40,42,44,44,46,48,49,50,50,51,29/)
   INTEGER, PARAMETER, DIMENSION(14) :: bndmax_sw=(/33,36,39,41,42,44,45,47,48,50,50,51,51,29/)
   INTEGER, PARAMETER, DIMENSION(16) :: bndmin_lw=(/11,12,16,18,20,22,23,23,24,25,26,28,28,28,29,30/)
   INTEGER, PARAMETER, DIMENSION(16) :: bndmax_lw=(/12,13,16,18,20,23,23,24,25,25,27,28,29,29,30,31/)
   
   INTEGER, PARAMETER :: namax=10
   
!   real(kind=rb) :: adjes              ! Flux adjustment for Earth/Sun distance
!   real(kind=rb) :: scon               ! Solar constant (W/m2)
!   !    1) Input in-cloud optical depth, single scattering albedo, asymmetry parameter,
!   !       and forward scattering fraction  directly (inflgsw = 0)
!   !       Values of tauc, ssac, asmc, fsfc are set 0 in rrtm_drive.
!   !    2) Input cloud fraction and cloud physical properties: ice fracion,
!   !       ice and liquid particle sizes (inflgsw = 1 or 2);  
!   integer(kind=im) :: inflgsw         ! Flag for cloud optical properties   
!   integer(kind=im) :: iceflgsw        ! Flag for ice particle specification (see reicmcl above)
!   integer(kind=im) :: liqflgsw        ! Flag for liquid droplet specification
!   !
!   integer(kind=im) :: inflglw         ! Flag for cloud optical properties
!   integer(kind=im) :: iceflglw        ! Flag for ice particle specification
!   integer(kind=im) :: liqflglw        ! Flag for liquid droplet specification
!   integer(kind=im) :: idrv            ! Flag for calculation of dFdT, the change

   REAL(KIND=rb), PARAMETER ::  scon=1.3533e03_rb
   REAL(KIND=rb), PARAMETER ::  adjes=1.0_rb
   INTEGER(kind=im), PARAMETER :: inflgsw  = 2 
   INTEGER(kind=im), PARAMETER :: iceflgsw = 3
   INTEGER(kind=im), PARAMETER :: liqflgsw = 1
   INTEGER(kind=im), PARAMETER :: inflglw =2
   INTEGER(kind=im), PARAMETER :: iceflglw =3
   INTEGER(kind=im), PARAMETER :: liqflglw =1
   INTEGER(kind=im), PARAMETER :: idrv=1
   ! icld = 0, clear only
   ! icld = 1, with clouds using random cloud overlap (McICA only)
   ! icld = 2, with clouds using maximum/random cloud overlap (McICA only)
   ! icld = 3, with clouds using maximum cloud overlap (McICA only)  
!kml 19 ago - Changing icld to the default option (3)
!   INTEGER(kind=im) :: icld=1
   INTEGER(kind=im) :: icld=3

!kml - 19 ago
! If the cloud generator is called multiple times, permute the seed between each call.
! Between calls for LW and SW, recommended permuteseed differes by 'ngpt'. This is done in rtm_driver
   INTEGER(kind=im) :: permuteseed=1
! flag for random number generator
!  0 = kissvec
!  1 = Mersenne Twister
   INTEGER(kind=im) :: irng=1
!kml - 19 ago
   
   LOGICAL :: firstTime=.true.
   
   CONTAINS
   
   subroutine initRRTM()
      use mem_chem1 , only: chemistry 
      USE chem1_list !, ONLY: nspecies,spc_name,O3,CO2,NO2,H2O,CH4,O2      
      use node_mod  , only: mxp,myp,mzp
      USE mem_grid  , only: maxnzp

      IMPLICIT NONE
      
      INTEGER :: i,k,m
      real(kind=rb) :: cpdair     ! Specific heat capacity of dry air
                                              ! at constant pressure at 273 K
                                              ! (J kg-1 K-1)
      
      o3isPresent =.false.
      co2isPresent=.false.
      ch4isPresent=.false.
      no2isPresent=.false.
      o2isPresent =.false.
      
      allocate (aot_rrtm_lw(mxp,myp,nbndlw))     ; aot_rrtm_lw=0.0
      allocate (aot_rrtm_sw(mxp,myp,nbndsw))     ; aot_rrtm_sw=0.0
      allocate (co3FromTuv(maxnzp+namax,mxp,myp)); co3FromTuv =0.0
      
      cpdair=cp
      call rrtmg_lw_ini(cpdair)
      CALL rrtmg_sw_ini(cpdair)

      IF(chemistry >= 0) THEN
         DO i=1,nspecies
            IF(trim(spc_name(i))=='O3') THEN
               o3isPresent=.true.
               o3pos=i
            END IF
            IF(trim(spc_name(i))=='CO2') THEN
               co2isPresent=.true.
               co2pos=i
            END IF
            IF(trim(spc_name(i))=='CH4') THEN
               ch4isPresent=.true.
               ch4pos=i
            END IF
            IF(trim(spc_name(i))=='NO2') THEN
               nO2isPresent=.true.
               no2pos=i
            END IF
            IF(trim(spc_name(i))=='O2') THEN
               O2isPresent=.true.
               o2pos=i
            END IF
            IF(trim(spc_name(i))=='H2O') THEN
               h2oisPresent=.true.
               h2opos=i
            END IF
         END DO
      END IF

      allocate(flip(mzp),sig(mzp+1),delsig(mzp+1),sigmid(mzp+1))
      flip  =0
      sig   =0.0
      delsig=0.0
      sigmid=0.0
      
      m=mzp
      DO k=1,mzp
         flip(k)=m
         m=m-1
      END DO
      
   END SUBROUTINE initRRTM
   
END MODULE mem_rrtm
