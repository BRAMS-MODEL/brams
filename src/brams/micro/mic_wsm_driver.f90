!==========================================================================================
!- This code is a 5, 6 and 7-class phase microphyiscs scheme of the 
!  Single-Moment MicroPhyiscs (WSMMP).
!- Adapted to BRAMS 6.0+ by Saulo Freitas Apr/2022
!- version from WRF 4.3.3
!==========================================================================================

SUBROUTINE micro_wsm( )

  use mem_basic, only : &
       basic_g            ! INTENT(INOUT)

  use mem_micro, only:  &
       micro_g            ! INTENT(INOUT)

  use mem_grid, only:   &
       ngrids,          & ! INTENT(IN)
       ngrid,           & ! INTENT(IN)
       zm,              & ! INTENT(IN)
       dzt,             & ! INTENT(IN)
       dtlt,            & ! INTENT(IN)
       jdim,            & ! INTENT(IN)
       maxnzp,          & ! INTENT(IN)
       time,            & ! INTENT(IN)
       zt,              & ! INTENT(IN)
       itime1,          & ! INTENT(IN)
       if_adap,         & ! INTENT(IN)
       grid_g,          & ! INTENT(IN)
       nnzp,npatch,imonth1,dtlongn,timmax! INTENT(IN)
!
  use node_mod, only :  &
       mzp,             & ! INTENT(IN)
       mxp,             & ! INTENT(IN)
       myp,             & ! INTENT(IN)
       ja,              & ! INTENT(IN)
       jz,              & ! INTENT(IN)
       ia,              & ! INTENT(IN)
       iz,              & ! INTENT(IN)
       mynum,i0,j0        ! INTENT(IN)


  use io_params, only : frqanl !INTENT(IN)
  
  use micphys,   only:  &
       mcphys_type        ! INTENT(IN)

  use mem_radiate, ONLY: ilwrtyp, iswrtyp ! INTENT(IN)

  use mem_leaf, only: leaf_g
   IMPLICIT NONE

   INTEGER,PARAMETER :: &
         IDS=1, IDE=2, JDS=1, JDE=2, KDS=1, &
         IMS=1, IME=2, JMS=1, JME=2, KMS=1, &
         ITS=1, ITE=1, JTS=1, JTE=1, KTS=1  !- rams 1st level is below surface => kts=2

   INTEGER :: KDE, &
              KME, &
              KTE

   INTEGER :: i,j

   LOGICAL :: diagflag=.false. 
   REAL :: ocean_fraction

   !- converting WRF setting to BRAMS
   !ids=1   ;ide=mxp ;jds=1   ;jde=myp ;kds=1; kde=mzp                 
   !ims=1   ;ime=mxp ;jms=1   ;jme=myp ;kms=1; kme=mzp                           
   !its=ia  ;ite=iz  ;jts=ja  ;jte=jz  ;kts=2; kte=mzp-1  
   !- converting WRF setting to BRAMS
   kde=mzp             
   kme=mzp                       
   kte=mzp-1  
   !
   !- flag for diagnostic time
   diagflag = .false.
   if(mod(time,frqanl)<dtlongn(1).or.time>=timmax - 0.01*dtlongn(1) ) then 
      diagflag = .true. 
   endif

   do j = ja,jz
      do i = ia,iz
       ocean_fraction = leaf_g(ngrid)%patch_area(i,j,1)
              
       call brams_to_mic_wsm(ia,ja,iz,jz, &
             mcphys_type &
            ,ilwrtyp     &
            ,iswrtyp     &
            ,j           &
            ,i           &
            ,IDS, IDE, JDS, JDE, KDS, KDE   &
            ,IMS, IME, JMS, JME, KMS, KME   &
            ,ITS, ITE, JTS, JTE, KTS, KTE   &
            ,mzp      &
            ,ngrid    &
            ,mynum    &
            ,if_adap  &
            !
            ,diagflag &
            !
            ,dtlt     &
            ,time     &
            ,zm       &
            ,dzt      &            
            ,zt       &
            ,basic_g(ngrid) &
            ,grid_g (ngrid) &
            ,micro_g(ngrid) &
	         ,ocean_fraction &
            )
      enddo
  enddo
  !- for consistency with surface and radiation schemes, the total
  !- precip will be also stored in the pcpg array
  micro_g(ngrid)%pcpg(:,:)=micro_g(ngrid)%pcprr(:,:)
   
 
END SUBROUTINE micro_wsm
!=======================================================================================
!
  SUBROUTINE brams_to_mic_wsm(ia,ja,iz,jz, &
             mcphys_type &
            ,ilwrtyp     &
            ,iswrtyp     &
            ,j &
            ,i &
            ,IDS, IDE, JDS, JDE, KDS, KDE   &
            ,IMS, IME, JMS, JME, KMS, KME   &
            ,ITS, ITE, JTS, JTE, KTS, KTE   &
            ,m1    &
            ,ngrid &
            ,mynum &
            ,if_adap  &
            !
            ,diagflag &
            !
            ,dtlt  &
            ,time  &
            ,zm    &
            ,dzt   &
            ,zt    &
            !
            ,basic &
            !
            ,grd &
            !
            ,mic &
            !
            ,ocean_fraction &
            )
         
   USE module_mp_wsm5, only: wsm5init, wsm5
   USE module_mp_wsm6, only: wsm6init, wsm6
   USE module_mp_wsm7, only: wsm7init, wsm7
   USE rconstants, only: p00,cpor,alvl,alvi,cpi,cpi4,cp253i
   use mem_basic, only : &
       basic_vars            ! INTENT(INOUT)
   use mem_micro, only:  &
       micro_vars            ! INTENT(INOUT)

   use mem_grid, only:   &
       grid_vars             ! INTENT(IN)
    
   
   IMPLICIT NONE

   type(basic_vars) ::basic
   type(grid_vars ) ::grd
   type(micro_vars) ::mic

   INTEGER, INTENT(IN) ::  &          
             mcphys_type   &
            ,ilwrtyp       &
            ,iswrtyp       &
            ,j &
            ,i &
            ,IDS, IDE, JDS, JDE, KDS, KDE   &
            ,IMS, IME, JMS, JME, KMS, KME   &
            ,ITS, ITE, JTS, JTE, KTS, KTE   &
            ,m1    &
            ,ngrid &
            ,mynum &
            ,if_adap ,ia,ja,iz,jz

   REAL, INTENT(IN) ::   &
             dtlt  &
            ,time  &
	        ,ocean_fraction 

   REAL,  INTENT(IN)   ,DIMENSION(m1) :: &
             zm    &
            ,dzt   &
            ,zt    
              
   LOGICAL, INTENT(IN) ::  diagflag

!- in the context of BRAMS, the variables below are "local"
   REAL, DIMENSION( ims:ime , kms:kme , jms:jme )  ::   &
                  th       &
                 ,dz8w     &
                 ,pi_phy   &
                 ,p        &
                 ,air_dens 

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ) ::                &
                  w                                               &
                 ,qv_curr,qc_curr,qr_curr,qi_curr,qs_curr,qg_curr &
                 ,qh_curr &
                 ,re_cloud, re_ice, re_snow, orho
                 !,qnc_curr, qnr_curr,  qni_curr          &
                 !,qnwfa_curr,qnifa_curr                  &
                     
   REAL, DIMENSION( ims:ime , jms:jme )  ::  &
                   RAINNC     &
                  ,RAINNCV    &
                  ,SNOWNC     &
                  ,SNOWNCV    &
                  ,GRAUPELNC  &
                  ,GRAUPELNCV &
                  ,HAIL       & 
                  ,HAILNCV    &
                  ,SR
                                                                   
!----------------------------------------------------------------------
! qv              water vapor    mixing ratio (kg/kg)
! qc              cloud water    mixing ratio (kg/kg)
! qr              rain water     mixing ratio (kg/kg)
! qi              cloud ice      mixing ratio (kg/kg)
! qs              snow            mixing ratio (kg/kg)
! qg              graupel            mixing ratio (kg/kg)
!
! qnc             cloud water number concentration (#/kg)
! qni             cloud ice   number concentration (#/kg)
! qnr             rain        number concentration (#/kg)
! qnwfa      water friendly aerosol number concentration (#/kg) - CCN
! qnifa      ice   friendly aerosol number concentration (#/kg) - IFN
!
!-- th            potential temperature    (K)
!-- w             vertical velocity (cartesian) (m/s)
!-- rho           density of air           (kg/m^3)
!-- pi_phy        exner function           (dimensionless)
!-- p             pressure                 (Pa)
!-- RAINNC        grid scale precipitation (mm)
!-- RAINNCV       one time step grid scale precipitation (mm/step)
!-- SNOWNC        grid scale snow and ice (mm)
!-- SNOWNCV       one time step grid scale snow and ice (mm/step)
!-- GRAUPELNC     grid scale graupel (mm)
!-- GRAUPELNCV    one time step grid scale graupel (mm/step)
!-- HAILNC        grid scale hail (mm)
!-- HAILNCV       one time step grid scale hail (mm/step)
!-- SR            one time step mass ratio of snow to total precip
!-- z             Height above sea level   (m)
!-- dt            Time step              (s)
!-- G             acceleration due to gravity  (m/s^2)
!-- CP            heat capacity at constant pressure for dry air (J/kg/K)
!-- R_d           gas constant for dry air (J/kg/K)
!-- R_v           gas constant for water vapor (J/kg/K)
!-- XLS           latent heat of sublimation   (J/kg)
!-- XLV           latent heat of vaporization  (J/kg)
!-- XLF           latent heat of melting       (J/kg)
!-- id            grid id number
!-- ids           start index for i in domain
!-- ide           end index for i in domain
!-- jds           start index for j in domain
!-- jde           end index for j in domain
!-- kds           start index for k in domain
!-- kde           end index for k in domain
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-- i_start       start indices for i in tile
!-- i_end         end indices for i in tile
!-- j_start       start indices for j in tile
!-- j_end         end indices for j in tile
!-- its           start index for i in tile
!-- ite           end index for i in tile
!-- jts           start index for j in tile
!-- jte           end index for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!-- diagflag      Logical to tell us when to produce diagnostics for history or restart
!-- rainprod      total tendency of conversion of cloud water/ice and graupel to rain (kg kg-1 s-1)

     REAL, DIMENSION(ims:ime,kms:kme,jms:jme) :: HGT,refl_10cm 
     REAL, DIMENSION(ims:ime, kms:kme, jms:jme):: rainprod,evapprod
     
     REAL :: tempK,rliq,rice,til,qhydm,tairstr
    ! INTEGER :: itimestep = 1    ! not used in mp_wsm6
     INTEGER :: do_radar_ref = 0 ! flag to compute radar reflectivity  
     INTEGER :: has_reqc, has_reqi, has_reqs ! flags to calculate effec radius
                                                        ! for radiation (1=ON, 0=OFF
     REAL :: dx, dy              ! grid spacing (m)
     REAL :: dt                  ! model timestep (s)
     LOGICAL :: start_of_simulation =.true. 
     integer, save ::it=0
     INTEGER  :: ke_diag              ! check this latter
     LOGICAL  :: wetscav_on = .false. ! check this latter

     REAL,  DIMENSION(m1) :: &
      thp    &
     ,theta  &
     ,pp     &
     ,rtp    &
     ,rv     &
     ,wp     &
     ,dn0    &
     ,pi0

     REAL :: rtgt

     REAL,DIMENSION(m1) :: &
      rcp     &
     ,rrp     &
     ,rpp     &
     ,rsp     &
     ,rgp     &
     ,rhp     &
     ,rei     &
     ,rel     &
     ,crp     &
     ,cpp     &
     ,ccp     &
     ,cccnp   &
     ,cifnp   


     REAL ::  &
      accpr   &! kg/m2 - rain+ice+snow+graupel
     ,pcprr   &! kg/m2 - rain+ice+snow+graupel
     ,accps   &! kg/m2 - ice+snow
     ,pcprs   &! kg/m2 - ice+snow
     ,accpg   &! kg/m2 - graupel
     ,pcprg   & ! kg/m2 - graupel
     ,accph   &
     ,pcprh 

     REAL, PARAMETER :: nt_c_ocean=100.E6 &
                       ,nt_c_land =200.E6 

     REAL :: nt_c_var
     REAL    , PARAMETER :: epsilon         = 1.E-15
     REAL    , PARAMETER :: g = 9.81  ! acceleration due to gravity (m {s}^-2)
  
     REAL    , PARAMETER :: r_d          = 287.       ! gas constant of dry air (J deg^-1 kg^-1)
     REAL    , PARAMETER :: cp           = 7.*r_d/2.  ! 
     REAL    , PARAMETER :: rhowater     = 1000.      ! density of liquid water at 0^oC (kg m^-3)
     REAL    , PARAMETER :: rhosnow      = 100.       ! density of snow (kg m^-3)
     REAL    , PARAMETER :: rhoair0      = 1.28       ! density of dry air at 0^oC and 1000mb pressure (kg m^-3)
     REAL    , PARAMETER :: r_v          = 461.6      ! gas constant for water vapor (J deg^-1 kg^-1)
     REAL    , PARAMETER :: cv           = cp-r_d     ! Specific heat of air at contant volume (J deg^-1 kg^-1)
     REAL    , PARAMETER :: cpv          = 4.*r_v
     REAL    , PARAMETER :: cvv          = cpv-r_v    ! 
     REAL    , PARAMETER :: cvpm         = -cv/cp
     REAL    , PARAMETER :: cliq         = 4190.      ! specific heat of liquid water at 0^oC
     REAL    , PARAMETER :: cice         = 2106.      ! specific heat of ice at 0^oC
     REAL    , PARAMETER :: psat         = 610.78
     REAL    , PARAMETER :: XLV0         = 3.15E6     !  constant defined for calculation of latent heating 
     REAL    , PARAMETER :: XLV1         = 2370.      !  constant defined for calculation of latent heating 
     REAL    , PARAMETER :: XLS0         = 2.905E6    !  constant defined for calculation of latent heating
     REAL    , PARAMETER :: XLS1         = 259.532    !  constant defined for calculation of latent heating
     REAL    , PARAMETER :: SVP1         = 0.6112     !  constant for saturation vapor pressure calculation (dimensionless)
     REAL    , PARAMETER :: SVP2         = 17.67      ! constant for saturation vapor pressure calculation (dimensionless)
     REAL    , PARAMETER :: SVP3         = 29.65      ! constant for saturation vapor pressure calculation (K)
     REAL    , PARAMETER :: SVPT0        = 273.15     ! constant for saturation vapor pressure calculation (K)

     REAL    , PARAMETER :: XLS          = 2.85E6     ! latent heat of sublimation of water at 0^oC (J kg^-1) 
     REAL    , PARAMETER :: XLV          = 2.5E6      ! latent heat of vaporization of water at 0^oC (J kg^-1)
     REAL    , PARAMETER :: XLF          = 3.50E5     ! latent heat of fusion of water at 0^oC (J kg^-1)
     REAL    , PARAMETER :: EP_1         = R_v/R_d-1. !  constant for virtual temperature (r_v/r_d - 1) (dimensionless)
     REAL    , PARAMETER :: EP_2         = R_d/R_v    ! constant for specific humidity calculation (dimensionless)

     INTEGER :: k,kr

        nt_c_var =    ocean_fraction *nt_c_ocean + &
	               (1.-ocean_fraction)*nt_c_land

        !- column quantities
	     thp  (1:m1)= basic%thp  (1:m1,i,j)
        theta(1:m1)= basic%theta(1:m1,i,j)
        pp   (1:m1)= basic%pp   (1:m1,i,j)
        rtp  (1:m1)= basic%rtp  (1:m1,i,j)
        rv   (1:m1)= basic%rv   (1:m1,i,j)
        wp   (1:m1)= basic%wp   (1:m1,i,j)
        dn0  (1:m1)= basic%dn0  (1:m1,i,j)
        pi0  (1:m1)= basic%pi0  (1:m1,i,j)
        !--- mass mixing ratio
        rcp  (1:m1)= mic%rcp    (1:m1,i,j)
        rrp  (1:m1)= mic%rrp    (1:m1,i,j)
        rpp  (1:m1)= mic%rpp    (1:m1,i,j)
        rsp  (1:m1)= mic%rsp    (1:m1,i,j)
        
	!- surface quantities
	     rtgt = grd%rtgt (i,j)
	     accpr= mic%accpr(i,j)
        pcprr= mic%pcprr(i,j)
        accps= mic%accps(i,j)
        pcprs= mic%pcprs(i,j)      

        if(mcphys_type == 7 )  then 
          rhp  (1:m1)= mic%rhp  (1:m1,i,j)
          accph      = mic%accph(i,j)
          pcprh      = mic%pcprh(i,j)
        else
          rhp  (1:m1)= 0.
          accph      = 0.
          pcprh      = 0.
        endif

        if(mcphys_type /= 5) then
          rgp  (1:m1)= mic%rgp  (1:m1,i,j)
          accpg      = mic%accpg(i,j)
          pcprg      = mic%pcprg(i,j)
        else
          rgp  (1:m1)= 0.
          accpg      = 0.
          pcprg      = 0.
        endif

!
!- for coupling with brams
!       !- converting WRF setting to BRAMS
!       ids=1   ;ide=mxp ;jds=1   ;jde=myp ;kds=1; kde=mzp              
!       ims=1   ;ime=mxp ;jms=1   ;jme=myp ;kms=1; kme=mzp                        
!       its=ia  ;ite=iz  ;jts=ja  ;jte=jz  ;kts=1; kte=mzp-1  

        ! flags to calculate effec radius
        IF( (ilwrtyp==6 .or. iswrtyp==6)) then           
           has_reqc= 1 ; has_reqi= 1 ; has_reqs= 1 
        ELSE
           has_reqc= 0 ; has_reqi= 0 ; has_reqs= 0 
        ENDIF
        
        dt= dtlt        ! time step            (s)
        rainprod  =0.0  ! for scaveging   aerosols/gases
        evapprod  =0.0  ! for evaporation aerosols/gases
        SR        =0.0  ! fraction of snow of the total water
                        ! ( for land surface models)
        refl_10cm =0.0  ! 
        ke_diag   = kte
       
        !- surface precipitation (total accumulated)
        RAINNC    (1,1)=  accpr !- rain+ice+snow+graupel+hail
        SNOWNC    (1,1)=  accps !- ice+snow
        GRAUPELNC (1,1)=  accpg !- graupel
        HAIL      (1,1)=  accph !- hail
        
        DO k=1,kme-1
          kr = k + 1
          qv_curr (1,k,1)= max(1.e-12,rtp(kr) - &    ! QV
                               (rcp(kr)+rrp(kr)+rpp(kr)+rsp(kr)+rgp(kr)+rhp(kr)))                    
          qc_curr (1,k,1)= max(0.0,rcp(kr))          ! QC     
          qr_curr (1,k,1)= max(0.0,rrp(kr))          ! QR   
          qi_curr (1,k,1)= max(0.0,rpp(kr))          ! QI   
          qs_curr (1,k,1)= max(0.0,rsp(kr))          ! QS   
          qg_curr (1,k,1)= max(0.0,rgp(kr))          ! QG

          qh_curr (1,k,1)= max(0.0,rhp(kr))          ! QH   
         
          pi_phy  (1,k,1)= (pp(kr)+pi0(kr))*cpi ! Exner function/cp (dimensionless)
          
          P       (1,k,1)= ( (pp(kr)+pi0(kr))*cpi )** cpor * p00        ! pressure(Pa)
         !W       (1,k,1)= wp(kr)       ! vertical velocity (m/s) ! must be at center or face? ASK
          
          dz8w    (1,k,1)= rtgt/dzt(kr) ! layer thickness (m) 

        ENDDO

        !- get potential temperature (theta) from theta_il (thp) and condensates
        DO k=1,kme -1 
           kr = k + 1
           tempK    = theta(kr)* (pp(kr)+pi0(kr))*cpi 
              til   = thp  (kr)* (pp(kr)+pi0(kr))*cpi 
         
           rliq     =  qc_curr(1,k,1) + qr_curr(1,k,1)                
           rice     =  qi_curr(1,k,1) + qs_curr(1,k,1) + qg_curr(1,k,1) + qh_curr(1,k,1)
           qhydm    =  alvl * rliq + alvi * rice
           
           if (tempK .gt. 253.) then
              tairstr = 0.5 * (til + sqrt(til * (til + cpi4 * qhydm)))
           else
              tairstr = til * (1. + qhydm * cp253i)
           endif
           !- updated potential temperature TH in Kelvin (adv+dif+rad+conv+)
           TH (1,k,1) = tairstr / pi_phy(1,k,1)

           !- air density
           air_dens(1,k,1) = P(1,k,1)/(287.04*tempK*(1.+0.608*qv_curr(1,k,1)))
          !air_dens(1,k,1)= dn0(kr) 

        ENDDO
        
        !
        IF(start_of_simulation) THEN !.or.restart.)   
           IF(mcphys_type == 5 ) &   
              CALL wsm5init(rhoair0,rhowater,rhosnow,cliq,cpv,   .false. )
           IF(mcphys_type == 6 ) &   
              CALL wsm6init(rhoair0,rhowater,rhosnow,cliq,cpv, 0,.false. )
           IF(mcphys_type == 7 ) &   
              CALL wsm7init(rhoair0,rhowater,rhosnow,cliq,cpv,   .false. )

           start_of_simulation =.false.
        ENDIF
        
        IF(mcphys_type == 5)                 &   
             
             CALL wsm5(                         &
                     TH,                        &! potential temperature    (K)
                     qv_curr,                   &! QV=qv_curr,     
                     qc_curr,                   &! QC=qc_curr,     
                     qr_curr,                   &! QR=qr_curr,     
                     qi_curr,                   &! QI=qi_curr,     
                     qs_curr,                   &! QS=qs_curr,     
                     !
                     air_dens,                  &          
                     pi_phy,                    &! exner function (dimensionless)
                     P,                         &! pressure(Pa)
                     dz8w,                      &! deltaz
                     !
                     dt,      &                  ! time step              (s)
                     g,       &
                     cp,      &
                     cpv,     &
                     r_d,     &
                     r_v,     &
                     svpt0,   &
                     ep_1,    &
                     ep_2,    &
                     epsilon, &
                     xls,     &
                     xlv,     &
                     xlf,     &
                     rhoair0, &
                     rhowater,&  
                     cliq,    &
                     cice,    &
                     psat,    & 
                    
                     RAINNC,                    &
                     RAINNCV,                   &
                     SNOWNC,                    &
                     SNOWNCV,                   &
                     SR,                        &
                     !
                     refl_10cm,                 &
                     diagflag,                  &
                     do_radar_ref,              &
                     ! 
                     has_reqc,                  & 
                     has_reqi,                  &  
                     has_reqs,                  & 
                     !
                     re_cloud,                  & 
                     re_ice,                    &
                     re_snow,                   &
                     IDS,IDE, JDS,JDE, KDS,KDE, &
                     IMS,IME, JMS,JME, KMS,KME, &
                     ITS,ITE, JTS,JTE, KTS,KTE  &
                     )
        
        IF(mcphys_type == 6)                    &   
             
             CALL wsm6(                         &
                     TH,                        &! potential temperature    (K)
                     qv_curr,                   &! QV=qv_curr,     
                     qc_curr,                   &! QC=qc_curr,     
                     qr_curr,                   &! QR=qr_curr,     
                     qi_curr,                   &! QI=qi_curr,     
                     qs_curr,                   &! QS=qs_curr,     
                     qg_curr,                   &! QG=qg_curr,     
                     
                     air_dens,                  &          
                     pi_phy,                    &! exner function (dimensionless)
                     P,                         &! pressure(Pa)
                     dz8w,                      &! deltaz

                     dt,      &                  ! time step              (s)
                     g,       &
                     cp,      &
                     cpv,     &
                     r_d,     &
                     r_v,     &
                     svpt0,   &
                     ep_1,    &
                     ep_2,    &
                     epsilon, &
                     xls,     &
                     xlv,     &
                     xlf,     &
                     rhoair0, &
                     rhowater,&  
                     cliq,    &
                     cice,    &
                     psat,    & 
                    
                     RAINNC,                    &
                     RAINNCV,                   &
                     SNOWNC,                    &
                     SNOWNCV,                   &
                     SR,                        &
                     !
                     refl_10cm,                 &
                     diagflag,                  &
                     do_radar_ref,              &
                     GRAUPELNC,                 & 
                     GRAUPELNCV,                &
                     !
                     !ke_diag,                  &
                     ! 
                     has_reqc,                  & 
                     has_reqi,                  &  
                     has_reqs,                  & 
                     !
                     re_cloud,                  & 
                     re_ice,                    &
                     re_snow,                   &
                     IDS,IDE, JDS,JDE, KDS,KDE, &
                     IMS,IME, JMS,JME, KMS,KME, &
                     ITS,ITE, JTS,JTE, KTS,KTE, &
                     !
                     wetscav_on,                &
                     evapprod,                  &
                     rainprod                   &
                     )

           IF(mcphys_type == 7)                 &   
             
             CALL wsm7(                         &
                     TH,                        &! potential temperature    (K)
                     qv_curr,                   &! QV=qv_curr,     
                     qc_curr,                   &! QC=qc_curr,     
                     qr_curr,                   &! QR=qr_curr,     
                     qi_curr,                   &! QI=qi_curr,     
                     qs_curr,                   &! QS=qs_curr,     
                     qg_curr,                   &! QG=qg_curr,     
                     qh_curr,                   &! Qh=qg_curr,     
                     
                     air_dens,                  &          
                     pi_phy,                    &! exner function (dimensionless)
                     P,                         &! pressure(Pa)
                     dz8w,                      &! deltaz

                     dt,      &                  ! time step              (s)
                     g,       &
                     cp,      &
                     cpv,     &
                     r_d,     &
                     r_v,     &
                     svpt0,   &
                     ep_1,    &
                     ep_2,    &
                     epsilon, &
                     xls,     &
                     xlv,     &
                     xlf,     &
                     rhoair0, &
                     rhowater,&  
                     cliq,    &
                     cice,    &
                     psat,    & 
                    
                     RAINNC,                    &
                     RAINNCV,                   &
                     SNOWNC,                    &
                     SNOWNCV,                   &
                     SR,                        &
                     !
                     refl_10cm,                 &
                     diagflag,                  &
                     do_radar_ref,              &
                     GRAUPELNC,                 & 
                     GRAUPELNCV,                &
                     HAIL,                      & 
                     HAILNCV,                   &
                     !
                     !ke_diag,                  &
                     ! 
                     has_reqc,                  & 
                     has_reqi,                  &  
                     has_reqs,                  & 
                     !
                     re_cloud,                  & 
                     re_ice,                    &
                     re_snow,                   &
                     IDS,IDE, JDS,JDE, KDS,KDE, &
                     IMS,IME, JMS,JME, KMS,KME, &
                     ITS,ITE, JTS,JTE, KTS,KTE  &
                     )


                                                                    
        !- updated variables after microphysics processes (from WSM to BRAMS)
        DO k=1,kme-1
         kr=k+1
         rtp(kr)= qv_curr(1,k,1) + &
                  qc_curr(1,k,1) + &     
                  qr_curr(1,k,1) + &    
                  qi_curr(1,k,1) + &    
                  qs_curr(1,k,1) + &    
                  qg_curr(1,k,1) + &    
                  qh_curr(1,k,1)  

         rcp(kr)= qc_curr(1,k,1)
         rrp(kr)= qr_curr(1,k,1)
         rpp(kr)= qi_curr(1,k,1)
         rsp(kr)= qs_curr(1,k,1)
         rgp(kr)= qg_curr(1,k,1)
         rhp(kr)= qh_curr(1,k,1)  
         
         rv (kr)= max(1.0e-12, rtp(kr) -(rcp(kr)+rrp(kr)+rpp(kr)+rsp(kr)+rgp(kr)+rhp(kr)))
      
         theta(kr) =  TH(1,k,1)
         tempK     =  TH(1,k,1)*pi_phy(1,k,1)
         
         rliq     =  qc_curr(1,k,1) + qr_curr(1,k,1)                 
         rice     =  qi_curr(1,k,1) + qs_curr(1,k,1) + qg_curr(1,k,1) + qh_curr(1,k,1)
         
         !- update liq-ice potential temperature THP in Kelvin including microphysics processes
         thp(kr)  =   TH(1,k,1)*(1. + alvl * rliq/(cp * max(tempK,253.))  &
                                    + alvi * rice/(cp * max(tempK,253.)) ) **(-1.0)      
        ENDDO
        !- definition for k=1
          rtp(1)  = rtp(2)  
          rcp(1)  = rcp(2)  
          rrp(1)  = rrp(2)  
          rpp(1)  = rpp(2)  
          rsp(1)  = rsp(2)  
          rgp(1)  = rgp(2)  
          rhp(1)  = rhp(2)
          rv (1)  = rv (2)  
          theta(1)= theta(2)
          thp(1)  = thp(2)  
        !
        IF( (ilwrtyp==6 .or. iswrtyp==6)) then           
          DO k=1,kme-1
            kr=k+1
            rel (kr) = re_cloud (1,k,1) * 1.e+6 ! RRTM requires in micrometer
            rei (kr) = re_ice   (1,k,1) * 1.e+6 ! RRTM requires in micrometer
          ENDDO
          rel (1) =rel (2) !;  rel (kme) =rel (kme-1) 
          rei (1) =rei (2) !;  rei (kme) =rei (kme-1)
        ENDIF        

        !- surface precipitation (units are kg/m^2 = mm)
        !- RAINNC and RAINNCV constains all precipitation hidrometeors (rain, graupel, snow, ...)
        accpr = RAINNC    (1,1) ! a = accum
        pcprr = RAINNCV   (1,1) ! p = for each dt  (or per time step)
        accps = SNOWNC    (1,1) 
        pcprs = SNOWNCV   (1,1) 
        accpg = GRAUPELNC (1,1) 
        pcprg = GRAUPELNCV(1,1) 
        accph = HAIL      (1,1) 
        pcprh = HAILNCV   (1,1) 

        !- column quantities
        basic%thp  (1:m1,i,j) =thp  (1:m1) 
        basic%theta(1:m1,i,j) =theta(1:m1)
        basic%rtp  (1:m1,i,j) =rtp  (1:m1)
        basic%rv   (1:m1,i,j) =rv   (1:m1)  
              
        mic%rcp    (1:m1,i,j) =rcp  (1:m1)   
        mic%rrp    (1:m1,i,j) =rrp  (1:m1)   
        mic%rpp    (1:m1,i,j) =rpp  (1:m1)   
        mic%rsp    (1:m1,i,j) =rsp  (1:m1)   

        IF(mcphys_type /= 5) mic%rgp(1:m1,i,j) =rgp  (1:m1)  
        IF(mcphys_type == 7) mic%rhp(1:m1,i,j) =rhp  (1:m1)
        
        if( (ilwrtyp==6 .or. iswrtyp==6)) then
         mic%rei   (1:m1,i,j) =rei  (1:m1)
         mic%rel   (1:m1,i,j) =rel  (1:m1)
        endif 

	!- surface quantities
        mic%accpr(i,j) = accpr !constains all precipitation hidrometeors (rain, graupel, snow, ...)
        mic%pcprr(i,j) = pcprr !constains all precipitation hidrometeors (rain, graupel, snow, ...)
        mic%accps(i,j) = accps
        mic%pcprs(i,j) = pcprs
        if(mcphys_type /= 5) then
          mic%accpg(i,j) = accpg
          mic%pcprg(i,j) = pcprg
        endif
        if(mcphys_type == 7) then
          mic%accph(i,j) = accph
          mic%pcprh(i,j) = pcprh
        endif
!
END  SUBROUTINE brams_to_mic_wsm
