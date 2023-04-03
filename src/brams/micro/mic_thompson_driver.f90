!==========================================================================================
!- Thompson and Eidhammer microphysics (aerosol aware).
!- Adapted to BRAMS 5.0+ by Saulo Freitas/Karla Longo in july-2015/Oct-2018.
!- version from WRF 3.7/3.8/4.01
!- ref: Thompson, G. and T. Eidhammer, 2014:  A study of aerosol impacts on clouds and
!-      and precipitation development in a large winter cyclone. J. Atmos. Sci., 71, 3636-3658.  
!==========================================================================================

SUBROUTINE micro_thompson( )

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
         ITS=1, ITE=1, JTS=1, JTE=1, KTS=2  !- rams 1st level is below surface => kts=2

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

       !if(mcphys_type == 3) then
       !  ccp1d  (1:mzp)   = max(micro_g(ngrid)%ccp  (1:mzp,i,j) , 0.) 
       !  cccnp1d(1:mzp)   = max(micro_g(ngrid)%cccnp(1:mzp,i,j) , 0.)   
       !  cifnp1d(1:mzp)   = max(micro_g(ngrid)%cifnp(1:mzp,i,j) , 0.) 
       ! print*,"mic1:",maxval(  ccp1d)
       !endif
              
       call brams_to_mic_thompson(ia,ja,iz,jz, &
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
       
       !if(mcphys_type == 3) then
       !   micro_g(ngrid)%ccp  (1:mzp,i,j)  = ccp1d  (1:mzp)
       !  !- commented out for now
       !  !micro_g(ngrid)%cccnp(1:mzp,i,j)  = cccnp1d(1:mzp)!checar necessidade pois sera um array nao         
       !  !micro_g(ngrid)%cifnp(1:mzp,i,j)  = cifnp1d(1:mzp)!checar necessidade pois sera um array nao         
       !endif
      
      enddo
  enddo
  !- for consistency with surface and radiation schemes, the total
  !- precip will be also stored in the pcpg array
  micro_g(ngrid)%pcpg(:,:)=micro_g(ngrid)%pcprr(:,:)
   
 
END SUBROUTINE micro_thompson
!=======================================================================================
!
  SUBROUTINE brams_to_mic_thompson(ia,ja,iz,jz, &
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
         
   USE module_mp_thompson, only : thompson_init, mp_gt_driver
   USE rconstants, only: p00,cp,cpor,alvl,alvi,cpi,cpi4,cp253i
   use mem_basic, only : &
       basic_vars            ! INTENT(INOUT)
   use mem_micro, only:  &
       micro_vars          ! INTENT(INOUT)

   use mem_grid, only:   &
       grid_vars            ! INTENT(IN)
    
   
   IMPLICIT NONE

   type(basic_vars) ::basic
   type(grid_vars)  ::grd
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
                 ,p

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ) ::                &
                  w                                               &
                 ,qv_curr,qc_curr,qr_curr,qi_curr,qs_curr,qg_curr &
                 ,qnc_curr, qnr_curr,  qni_curr          &
                 ,qnwfa_curr,qnifa_curr                  &
                 ,re_cloud, re_ice, re_snow, orho
                     
   REAL, DIMENSION( ims:ime , jms:jme )  ::  &
                   RAINNC     &
                  ,RAINNCV    &
                  ,SNOWNC     &
                  ,SNOWNCV    &
                  ,GRAUPELNC  &
                  ,GRAUPELNCV &
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

     REAL, ALLOCATABLE,DIMENSION(:,:,:) ,save::  NWFA , NIFA
     REAL, ALLOCATABLE,DIMENSION(:,:)   ,save:: NWFA2D 

     REAL, DIMENSION(ims:ime, jms:jme):: qnwfa2d, qnifa2d
     REAL, DIMENSION(ims:ime,kms:kme,jms:jme) :: HGT,refl_10cm 
     REAL, DIMENSION(ims:ime, kms:kme, jms:jme):: rainprod,evapprod
     
     REAL :: tempK,rliq,rice,til,qhydm,tairstr
     INTEGER :: itimestep = 1    ! not used in mp_thompson
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
     ,pcprg    ! kg/m2 - graupel
     REAL, PARAMETER :: nt_c_ocean=100.E6 &
                       ,nt_c_land =200.E6 

     REAL :: nt_c_var

     !-tmp for output
     real, dimension(1,kme,100):: vargrads
     character (len=40),dimension(100,2) :: gradsname
     INTEGER :: k,jl,jk,nvar3d,nvar2d,nvar,nvx,nrec,gate
!_____________________________________________________________
!if(i==iz .and. j==jz .and. maxval( mic%rrp  )>1.e-6) then
!if(mynum == 72 .and. maxval( basic%wp)>10.) then
!	 print*,"1THP=" ,maxval(basic%thp),minval(basic%thp),i,j
!	 print*,"1theta",maxval( basic%theta),minval( basic%theta),mynum
!	 print*,"1pp "  ,maxval( basic%pp   ),minval( basic%pp   ),mynum
!	 print*,"1rtp"  ,maxval( basic%rtp  ),minval( basic%rtp  ),mynum
!	 print*,"1rv "  ,maxval( basic%rv   ),minval( basic%rv   ),mynum
!	 print*,"1wp "  ,maxval( basic%wp   ),minval( basic%wp   ),mynum
!	 !print*,"1dn0"  ,maxval( basic%dn0  ),minval( basic%dn0  ),mynum
!	 !print*,"1pi0"  ,maxval( basic%pi0  ),minval( basic%pi0  ),mynum
!
!	 print*,"1rcp " ,maxval( mic%rcp  )  ,minval( mic%rcp  ),mynum
!	 print*,"1rrp " ,maxval( mic%rrp  )  ,minval( mic%rrp  ),mynum
!	 print*,"1rpp " ,maxval( mic%rpp  )  ,minval( mic%rpp  ),mynum
!	 print*,"1rsp " ,maxval(mic%rsp  )   ,minval(mic%rsp  ),mynum
!	 print*,"1rgp " ,maxval(mic%rgp  )   ,minval(mic%rgp  ),mynum
!		     
!	 print*,"1crp"  ,maxval(mic%crp )    ,minval(mic%crp ),mynum
!	 print*,"1cpp"  ,maxval(mic%cpp )    ,minval(mic%cpp ),mynum
!	call flush(6)
!endif
!if( minval(mic%cpp(:,i,j)) < 0.) then
!     print*,"cpp1=",minval(mic%cpp(:,i,j)),i,j,mynum,it
!     call flush(6)
!endif    
!________________________________________________________

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
        rgp  (1:m1)= mic%rgp    (1:m1,i,j)
        
	!--- number concentration
        crp  (1:m1)= mic%crp    (1:m1,i,j) ! rain        number concentration (#/kg)
        cpp  (1:m1)= mic%cpp    (1:m1,i,j) ! cloud ice   number concentration (#/kg)
      
        if(mcphys_type == 3) then 
	     ccp  (1:m1) = max(mic%ccp  (1:m1,i,j) , 0.) ! cloud water        number concentration (#/kg)
             cccnp(1:m1) = max(mic%cccnp(1:m1,i,j) , 0.) ! water friendly aer number concentration (#/kg)  
             cifnp(1:m1) = max(mic%cifnp(1:m1,i,j) , 0.) ! ice   friendly aer number concentration (#/kg)  
        endif
        
	!- surface quantities
	     rtgt = grd%rtgt(i,j)
	     accpr= mic%accpr(i,j)
        pcprr= mic%pcprr(i,j)
        accps= mic%accps(i,j)
        pcprs= mic%pcprs(i,j)
        accpg= mic%accpg(i,j)
        pcprg= mic%pcprg(i,j)

!- for coupling with brams
!       !- converting WRF setting to BRAMS
!       ids=1   ;ide=mxp ;jds=1   ;jde=myp ;kds=1; kde=mzp              
!       ims=1   ;ime=mxp ;jms=1   ;jme=myp ;kms=1; kme=mzp                        
!       its=ia  ;ite=iz  ;jts=ja  ;jte=jz  ;kts=2; kte=mzp-1  

        ! flags to calculate effec radius
        IF( (ilwrtyp==6 .or. iswrtyp==6) .and. mcphys_type == 2 ) then           
           has_reqc= 1 ; has_reqi= 1 ; has_reqs= 1 
        ELSE
           has_reqc= 0 ; has_reqi= 0 ; has_reqs= 0 
        ENDIF
        
        dt= dtlt         ! time step            (s)
        
        rainprod  =0.0  ! for scaveging   aerosols/gases
        evapprod  =0.0  ! for evaporation aerosols/gases
        SR        =0.0  ! fraction of snow of the total water
                        ! ( for land surface models)
        
        refl_10cm =0.0  ! 

        ke_diag   = kte
        
        dx=10000. !- typical x- horizontal grid spacing (only for the local aerosol emission)
        dy=10000. !- typical y- horizontal grid spacing (only for the local aerosol emission)

        !- surface precipitation (total accumulated)
        RAINNC    (1,1)=  accpr !- rain+ice+snow+graupel
        SNOWNC    (1,1)=  accps !- ice+snow
        GRAUPELNC (1,1)=  accpg !- graupel
        
        
        DO k=1,kme
          qv_curr (1,k,1)= max(1.e-12,rtp(k) - &    ! QV
                               (rcp(k)+rrp(k)+rpp(k)+rsp(k)+rgp(k)))                    
          qc_curr (1,k,1)= max(0.0,rcp(k))          ! QC     
          qr_curr (1,k,1)= max(0.0,rrp(k))          ! QR   
          qi_curr (1,k,1)= max(0.0,rpp(k))          ! QI   
          qs_curr (1,k,1)= max(0.0,rsp(k))          ! QS   
          qg_curr (1,k,1)= max(0.0,rgp(k))          ! QG   
          
          qni_curr(1,k,1)= max(0.0,cpp(k))          ! NI    
          qnr_curr(1,k,1)= max(0.0,crp(k))          ! NR  
          
          pi_phy  (1,k,1)= (pp(k)+pi0(k))*cpi ! Exner function/cp (dimensionless)
          
          P       (1,k,1)= (pi_phy(1,k,1))** cpor * p00        ! pressure(Pa)
          W       (1,k,1)= wp(k)       ! vertical velocity (m/s) ! must be at center or face? ASK
          
          dz8w    (1,k,1)= rtgt/dzt(k) ! layer thickness (m) 
        ENDDO

        !- special setting for height - keep the range (:,k,:) 
        DO k=1,kme
          HGT     (:,k,:)=  zt(k)*rtgt  ! height above local surface (m)
        
        ENDDO
        !- get potential temperature (theta) from theta_il (thp) and condensates
        DO k=1,kme
           tempK    = theta(k)*pi_phy(1,k,1)
              til   = thp  (k)*pi_phy(1,k,1) 
         
           rliq     =  qc_curr(1,k,1) + qr_curr(1,k,1)                
           rice     =  qi_curr(1,k,1) + qs_curr(1,k,1)+qg_curr(1,k,1)
           qhydm    =  alvl * rliq + alvi * rice
           
           if (tempK .gt. 253.) then
              tairstr = 0.5 * (til + sqrt(til * (til + cpi4 * qhydm)))
           else
              tairstr = til * (1. + qhydm * cp253i)
           endif
           !- updated potential temperature TH in Kelvin (adv+dif+rad+conv+)
           TH (1,k,1) = tairstr / pi_phy(1,k,1)
        ENDDO
        
        !
        IF(start_of_simulation) THEN !.or.restart.or.config_flags%cycling)     &
          IF(mcphys_type == 2 ) then   
          
            CALL thompson_init(HGT,         &
                          DX, DY,           &
                          start_of_simulation,            &
                          IDS, IDE, JDS, JDE, KDS, KDE,   &
                          IMS, IME, JMS, JME, KMS, KME,   &
                          ITS, ITE, JTS, JTE, KTS, KTE)
          ENDIF
          
          IF(mcphys_type == 3 ) then  
            if( .not. allocated(NWFA)  ) allocate(NWFA  (ims:ime,kms:kme,jms:jme))
            if( .not. allocated(NIFA)  ) allocate(NIFA  (ims:ime,kms:kme,jms:jme))
            if( .not. allocated(NWFA2D)) allocate(NWFA2D(ims:ime,jms:jme))
            NWFA=0.0
            NIFA=0.0
            NWFA2D=0.0
	    orho(1,1:m1,1) =1./dn0(1:m1)
            CALL thompson_init(HGT,         &
                          DX, DY,           &
                          start_of_simulation,            &
                          IDS, IDE, JDS, JDE, KDS, KDE,   &
                          IMS, IME, JMS, JME, KMS, KME,   &
                          ITS, ITE, JTS, JTE, KTS, KTE,   &
                          orho,             &
			                 NWFA2D,           &
                          NWFA,             &
                          NIFA)

          ENDIF
          start_of_simulation =.false.
        ENDIF
                
        if(mcphys_type == 3) then
          ! - CCN and IN fields
          qnwfa_curr(1,:,1)=NWFA(1,:,1) !- CCN
          qnifa_curr(1,:,1)=NIFA(1,:,1) !- IN
          qnwfa2d   = 0.0 ! no emission for aerosols.
          qnifa2d   = 0.0 ! no emission for aerosols.
          !
          DO k=1,kme
           qnc_curr  (1,k,1) = ccp  (k) ! NC
           qnwfa_curr(1,k,1) = cccnp(k) !- CCN 
           qnifa_curr(1,k,1) = cifnp(k) !- IN  
          ENDDO
        endif

        !- this call cloud water 1-mom/ no aerosol aware scheme
        IF(mcphys_type == 2 ) & 
           CALL mp_gt_driver(                   &
                     qv_curr,                   &! QV=qv_curr,     
                     qc_curr,                   &! QC=qc_curr,     
                     qr_curr,                   &! QR=qr_curr,     
                     qi_curr,                   &! QI=qi_curr,     
                     qs_curr,                   &! QS=qs_curr,     
                     qg_curr,                   &! QG=qg_curr,     
                     qni_curr,                  &! NI=qni_curr,    
                     qnr_curr,                  &! NR=qnr_curr,    
!-these are optional arrays (only for mcphys_type == 3)
!                    qnc_curr,                  &! NC=qnc_curr,     
!                    qnwfa_curr,                &! NWFA=qnwfa_curr, 
!                    qnifa_curr,                &! NIFA=qnifa_curr, 
!                    qnwfa2d,                   &! NWFA2D=qnwfa2d,  
!-
                     TH,                        &! potential temperature    (K)
                     pi_phy,                    &! exner function (dimensionless)
                     P,                         &! pressure(Pa)
                     W,                         &
                     dz8w,                      &
                     dt,                        &! time step              (s)
                     itimestep,                 &
                     RAINNC,                    &
                     RAINNCV,                   &
                     SNOWNC,                    &
                     SNOWNCV,                   &
                     GRAUPELNC,                 & 
                     GRAUPELNCV,                & 
                     SR,                        &

                     wetscav_on,                &
                     
                     rainprod,                  &
                     evapprod,                  &
                     refl_10cm,                 &
                     diagflag,                  &

                     ke_diag,                   &

                     do_radar_ref,              &
                     re_cloud,                  & 
                     re_ice,                    &
                     re_snow,                   &
                     has_reqc,                  & ! G. Thompson
                     has_reqi,                  & ! G. Thompson
                     has_reqs,                  & ! G. Thompson
                     IDS,IDE, JDS,JDE, KDS,KDE, &
                     IMS,IME, JMS,JME, KMS,KME, &
                     ITS,ITE, JTS,JTE, KTS,KTE  &
!		               nt_c_var&
                     )


        !- this call cloud water 2-mom / aerosol aware scheme
        !- there are 4 additional arrays in the end of the argument list.
        IF(mcphys_type == 3 ) &        
           CALL mp_gt_driver(                   &
                     qv_curr,                   &! QV=qv_curr,     
                     qc_curr,                   &! QC=qc_curr,     
                     qr_curr,                   &! QR=qr_curr,     
                     qi_curr,                   &! QI=qi_curr,     
                     qs_curr,                   &! QS=qs_curr,     
                     qg_curr,                   &! QG=qg_curr,     
                     qni_curr,                  &! NI=qni_curr,    
                     qnr_curr,                  &! NR=qnr_curr,    
!-these are optional arrays and were moved to the end of this list
!                    qnc_curr,                  &! NC=qnc_curr,     
!                    qnwfa_curr,                 &! NWFA=qnwfa_curr, 
!                    qnifa_curr,                 &! NIFA=qnifa_curr, 
!                    qnwfa2d,                   &! NWFA2D=qnwfa2d,  
!-
                     TH,                        &! potential temperature    (K)
                     pi_phy,                    &! exner function (dimensionless)
                     P,                         &! pressure(Pa)
                     W,                         &
                     dz8w,                      &
                     dt,                        &! time step              (s)
                     itimestep,                 &
                     RAINNC,                    &
                     RAINNCV,                   &
                     SNOWNC,                    &
                     SNOWNCV,                   &
                     GRAUPELNC,                 & 
                     GRAUPELNCV,                & 
                     SR,                        &

                     wetscav_on,                &
                     
                     rainprod,                  &
                     evapprod,                  &
                     refl_10cm,                 &
                     diagflag,                  &

                     ke_diag,                   &

                     do_radar_ref,              &
                     re_cloud,                  & 
                     re_ice,                    &
                     re_snow,                   &
                     has_reqc,                  & ! G. Thompson
                     has_reqi,                  & ! G. Thompson
                     has_reqs,                  & ! G. Thompson
                     IDS,IDE, JDS,JDE, KDS,KDE, &
                     IMS,IME, JMS,JME, KMS,KME, &
                     ITS,ITE, JTS,JTE, KTS,KTE, &
!		               nt_c_var,                  &
!- moving optional arrays to the end of the argument list
                     qnc_curr,                  &! NC=qnc_curr,     
                     qnwfa_curr,                &! NWFA=qnwfa_curr, 
                     qnifa_curr,                &! NIFA=qnifa_curr, 
                     qnwfa2d,                   &! NWFA2D=qnwfa2d,  
                     qnifa2d                    &! NiFA2D=qnifa2d,  
                                                )
                                                                    
        !- updated variables after microphysics processes
        DO k=2,kme
         rtp(k)= qv_curr(1,k,1) + &
                 qc_curr(1,k,1) + &     
                 qr_curr(1,k,1) + &    
                 qi_curr(1,k,1) + &    
                 qs_curr(1,k,1) + &    
                 qg_curr(1,k,1)     
         
         rcp(k)= qc_curr(1,k,1)
         rrp(k)= qr_curr(1,k,1)
         rpp(k)= qi_curr(1,k,1)
         rsp(k)= qs_curr(1,k,1)
         rgp(k)= qg_curr(1,k,1)
         rv (k)= max(1.0e-12, rtp(k) -(rcp(k)+rrp(k)+rpp(k)+rsp(k)+rgp(k)))
          
         cpp(k)= qni_curr(1,k,1)
         crp(k)= qnr_curr(1,k,1)
          
         theta(k) =  TH(1,k,1)
         tempK    =  theta(k)*pi_phy(1,k,1)
         
         rliq     =  qc_curr(1,k,1) + qr_curr(1,k,1)                 
         rice     =  qi_curr(1,k,1) + qs_curr(1,k,1) + qg_curr(1,k,1)
         
         !- update liq-ice potential temperature THP in Kelvin including microphysics processes
         thp(k)   =  theta(k)*(1. + alvl * rliq/(cp * max(tempK,253.))  &
                                  + alvi * rice/(cp * max(tempK,253.)) ) **(-1.0)      
        ENDDO
        !- definition for k=1
          rtp(1)  = rtp(2)  
          rcp(1)  = rcp(2)  
          rrp(1)  = rrp(2)  
          rpp(1)  = rpp(2)  
          rsp(1)  = rsp(2)  
          rgp(1)  = rgp(2)  
          rv (1)  = rv (2)  
          cpp(1)  = cpp(2)  
          crp(1)  = crp(2)  
          theta(1)= theta(2)
          thp(1)  = thp(2)  
        !
        !
        IF(mcphys_type == 3) then
          DO k=2,kme
            ccp  (k) = qnc_curr  (1,k,1) 
            cccnp(k) = qnwfa_curr(1,k,1)                                          
            cifnp(k) = qnifa_curr(1,k,1) 
          ENDDO
          !- definition for k=1
          ccp  (1) = ccp  (2)
          cccnp(1) = cccnp(2)
          cifnp(1) = cifnp(2) 
        ENDIF        
        
        IF( (ilwrtyp==6 .or. iswrtyp==6) .and. mcphys_type == 2 ) then           
          DO k=2,kme-1
            rel (k) = re_cloud (1,k,1) * 1.e+6 ! RRTM requires in micrometer
            rei (k) = re_ice   (1,k,1) * 1.e+6 ! RRTM requires in micrometer
          ENDDO
          rel (1) =rel (2);  rel (kme) =rel (kme-1) 
          rei (1) =rei (2);  rei (kme) =rei (kme-1)
        ENDIF        

        !- surface precipitation (units are kg/m^2 = mm)
        accpr = RAINNC    (1,1) ! a = accum
        pcprr = RAINNCV   (1,1) ! p = for each dt  (or per time step)
        accps = SNOWNC    (1,1) 
        pcprs = SNOWNCV   (1,1) 
        accpg = GRAUPELNC (1,1) 
        pcprg = GRAUPELNCV(1,1) 

        !- column quantities
        basic%thp  (1:m1,i,j) =thp  (1:m1) 
        basic%theta(1:m1,i,j) =theta(1:m1)
        basic%rtp  (1:m1,i,j) =rtp  (1:m1)
        basic%rv   (1:m1,i,j) =rv   (1:m1)  
              
        mic%rcp    (1:m1,i,j) =rcp  (1:m1)   
        mic%rrp    (1:m1,i,j) =rrp  (1:m1)   
        mic%rpp    (1:m1,i,j) =rpp  (1:m1)   
        mic%rsp    (1:m1,i,j) =rsp  (1:m1)   
        mic%rgp    (1:m1,i,j) =rgp  (1:m1)  
        
        if( (ilwrtyp==6 .or. iswrtyp==6) .and. mcphys_type == 2 ) then  	 
         mic%rei   (1:m1,i,j) =rei  (1:m1)
         mic%rel   (1:m1,i,j) =rel  (1:m1)
        endif 
        
        mic%crp    (1:m1,i,j) =crp  (1:m1)
        mic%cpp    (1:m1,i,j) =cpp  (1:m1)
        
	     if(mcphys_type == 3) then
         mic%ccp  (1:m1,i,j) = ccp  (1:m1) 
         mic%cccnp(1:m1,i,j) = cccnp(1:m1) 
         mic%cifnp(1:m1,i,j) = cifnp(1:m1) 
        endif
	
	!- surface quantities
        mic%accpr(i,j) = accpr
        mic%pcprr(i,j) = pcprr
        mic%accps(i,j) = accps
        mic%pcprs(i,j) = pcprs
        mic%accpg(i,j) = accpg
        mic%pcprg(i,j) = pcprg

!if(i==iz .and. j==jz .and. maxval( mic%rrp  )>1.e-6) then
!if(mynum == 72 .and. maxval( basic%wp)>10.)then
!	 print*,"2THP=" ,maxval(basic%thp),minval(basic%thp),i,j
!	 print*,"2theta",maxval( basic%theta),minval( basic%theta),mynum
!	 print*,"2pp "  ,maxval( basic%pp   ),minval( basic%pp   ),mynum
!	 print*,"2rtp"  ,maxval( basic%rtp  ),minval( basic%rtp  ),mynum
!	 print*,"2rv "  ,maxval( basic%rv   ),minval( basic%rv   ),mynum
!	 print*,"2wp "  ,maxval( basic%wp   ),minval( basic%wp   ),mynum
!	 !print*,"2dn0"  ,maxval( basic%dn0  ),minval( basic%dn0  ),mynum
!	 !print*,"2pi0"  ,maxval( basic%pi0  ),minval( basic%pi0  ),mynum
!
!	 print*,"2rcp " ,maxval( mic%rcp  )  ,minval( mic%rcp  ),mynum
!	 print*,"2rrp " ,maxval( mic%rrp  )  ,minval( mic%rrp  ),mynum
!	 print*,"2rpp " ,maxval( mic%rpp  )  ,minval( mic%rpp  ),mynum
!	 print*,"2rsp " ,maxval(mic%rsp  )   ,minval(mic%rsp  ),mynum
!	 print*,"2rgp " ,maxval(mic%rgp  )   ,minval(mic%rgp  ),mynum
!		     
!	 print*,"2crp"  ,maxval(mic%crp )    ,minval(mic%crp ),mynum
!	 print*,"2cpp"  ,maxval(mic%cpp )    ,minval(mic%cpp ),mynum
!	call flush(6)
!endif
!if(  minval(mic%cpp(:,i,j)) < 0.) then
!     print*,"cpp2=",minval(mic%cpp(:,i,j)),i,j,mynum
!     call flush(6)
!endif    
!if(i==iz .and. j==jz ) it=it+1
!if(i==iz .and. j==jz .and. it==2 ) stop 33
  RETURN
!
!
END  SUBROUTINE brams_to_mic_thompson
