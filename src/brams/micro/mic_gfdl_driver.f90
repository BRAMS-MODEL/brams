!==========================================================================================
!- GFDL cloud microphysics 
!- Adapted to BRAMS 5.0+ by Saulo Freitas 2018/2021
!==========================================================================================

SUBROUTINE micro_gfdl( )

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
       deltaxn,         & ! INTENT(IN)
       deltayn,         & ! INTENT(IN)
       time,            & ! INTENT(IN)
       zt,              & ! INTENT(IN)
       itime1,          & ! INTENT(IN)
       if_adap,         & ! INTENT(IN)
       grid_g,          & ! INTENT(IN)
       nnzp,npatch,imonth1,dtlongn,timmax! INTENT(IN)

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
 
  use mem_radiate, ONLY: ilwrtyp, iswrtyp ! INTENT(IN)

   IMPLICIT NONE
   INTEGER :: ims, ime, jms, jme, kms, kme  
   
   !---bounds for the GFDL microphysics
   ims = 1; ime = mxp
   jms = 1; jme = myp 
   kms = 1; kme = mzp-2
   	  
   call brams_to_mic_gfdl(ia,ja,iz,jz	&
   	,ilwrtyp     &
   	,iswrtyp     &
   	,ims, ime, jms, jme, kms, kme	&
   	,mzp	    &
   	,ngrid    &
   	,mynum    &
   	!
   	,deltaxn(ngrid) &   ! INTENT(IN)
   	,deltayn(ngrid) &   ! INTENT(IN)
   	!
   	,dtlt	  &
   	,time	  &
   	,zm	  &
   	,dzt	  &	       
   	,zt	  &
   	,basic_g(ngrid) &
   	,grid_g (ngrid) &
   	,micro_g(ngrid) &
   	!
   	)

   !- for consistency with surface and radiation schemes, the total
   !- precip will be also stored in the pcpg array
   micro_g(ngrid)%pcpg(:,:)=micro_g(ngrid)%pcprr(:,:)
   
END SUBROUTINE micro_gfdl
!=======================================================================================
!
  SUBROUTINE brams_to_mic_gfdl(ia,ja,iz,jz  &
            ,ilwrtyp     &
            ,iswrtyp     &
            ,ims, ime, jms, jme, kms, kme   &
            ,m1    &
            ,ngrid &
            ,mynum &
            !
            ,deltaxn&   ! INTENT(IN)
            ,deltayn&   ! INTENT(IN)
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
            !
            )
         
   use gfdl_cloud_microphys_mod, only : gfdl_cloud_microphys_driver &
                                      , gfdl_cloud_microphys_init
  
   use rconstants, only: p00,cp,cpor,alvl,alvi,cpi,cpi4,cp253i
   use mem_basic, only : &
       basic_vars            ! INTENT(INOUT)
   use mem_micro, only:  &
       micro_vars          ! INTENT(INOUT)

   use mem_grid, only:   &
       grid_vars            ! INTENT(IN)
    
   use mem_leaf          , only: leaf_g, isfcl
   
   IMPLICIT NONE

   type(basic_vars) ::basic
   type(grid_vars)  ::grd
   type(micro_vars) ::mic

   INTEGER, INTENT(IN) ::  &          
             ilwrtyp       &
            ,iswrtyp       &
            ,ims, ime, jms, jme, kms, kme   &
            ,m1    &
            ,ngrid &
            ,mynum &
            ,ia,ja,iz,jz

   REAL, INTENT(IN) ::   &
             dtlt        &
            ,time        &
            ,deltaxn     &  
            ,deltayn   

   REAL,  INTENT(IN)   ,DIMENSION(m1) :: &
             zm    &
            ,dzt   &
            ,zt    
              
   REAL, DIMENSION( m1, ims:ime , jms:jme )  ::   &
                 exner


   REAL, DIMENSION( ims:ime , jms:jme, kms:kme )  ::   &
                  temp,dz,w,u,v, dp

   REAL, DIMENSION( ims:ime, jms:jme , kms:kme) ::                  &
                  qv_curr,qc_curr,qr_curr,qi_curr,qs_curr,qg_curr , &
                  re_cloud, re_ice ,re_snow, rad_cf , &
                  DQVDT_micro , & 
		            DQLDT_micro , & 
		            DQRDT_micro , & 
		            DQIDT_micro , & 
		            DQSDT_micro , & 
		            DQGDT_micro , & 
		            DQADT_micro , & 
		             DUDT_micro , & 
		             DVDT_micro , & 
		             DTDT_micro , &
		             NACTLI     , &
                   REV_MC_X, RSU_MC_X, EVAPC_X
		     
   REAL, DIMENSION( ims:ime , jms:jme )  ::  &
                   RAINNC     &
                  ,RAINNCV    &
                  ,SNOWNC     &
                  ,SNOWNCV    &
                  ,GRAUPELNC  &
                  ,GRAUPELNCV 

   REAL, DIMENSION( ims:ime , jms:jme )  ::  &
                   area    &
		            ,frland  &
		            ,cnv_fraction
                                                                   
   real, dimension(ims:ime , jms:jme ,0:kme) :: &
                  PFI_LS_X &
                 ,PFL_LS_X 

   real, dimension(ims:ime , jms:jme ) :: PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL
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
!-- press         pressure                 (Pa)
!-- RAINNC        grid scale precipitation (mm)
!-- RAINNCV       one time step grid scale precipitation (mm/step)
!-- SNOWNC        grid scale snow and ice (mm)
!-- SNOWNCV       one time step grid scale snow and ice (mm/step)
!-- GRAUPELNC     grid scale graupel (mm)
!-- GRAUPELNCV    one time step grid scale graupel (mm/step)
!-- HAILNC        grid scale hail (mm)
!-- HAILNCV       one time step grid scale hail (mm/step)
!-- z             Height above sea level   (m)
!-- dt            Time step              (s)


     real :: tempk,rliq,rice,til,qhydm,tairstr,exner_tmp,rcond
     real :: dx, dy		    ! grid spacing (m)
     real :: dt_moist		  ! model timestep (s)
     logical :: start_of_simulation =.true. 
     integer, allocatable, dimension(:),save:: flip
     real   , dimension(m1)	:: dp_tmp

     logical  :: lhydrostatic	   =.false.
     logical  :: lphys_hydrostatic =.false.
     integer :: i,j,k,m
     real, parameter :: ANV_ICEFALL=1.0, LS_ICEFALL=1.0
     

	 IF(start_of_simulation) THEN
	   !-initialization
	   call gfdl_cloud_microphys_init()
	   !- define the vector "flip" to invert the z-axis orientation
	   allocate(flip(m1)); flip(:)=-9999
     call flipz_minus1(flip,kme+1)
	   !-
	   start_of_simulation =.false.
	 ENDIF

   !---- zero-out microphysics tendencies
   DQVDT_micro = 0.
   DQLDT_micro = 0.
   DQRDT_micro = 0.
   DQIDT_micro = 0.
   DQSDT_micro = 0.
   DQGDT_micro = 0.
   DQADT_micro = 0.
    DUDT_micro = 0.
    DVDT_micro = 0.
    DTDT_micro = 0.
   !--- zero-out 3D Precipitation Fluxes 
   PFI_LS_X = 0. ! Ice
   PFL_LS_X = 0. ! Liquid
   !--- output rain re-evaporation and sublimation
   REV_MC_X= 0. ;RSU_MC_X= 0. ;EVAPC_X= 0. 

   dt_moist= dtlt  ! time step   (s)
                
	 DO j = jms,jme ; DO i = ims,ime 
            
            area	      (i,j) = deltaxn*deltayn
            frland	   (i,j) = 1.-leaf_g(ngrid)%patch_area(i,j,1)
            cnv_fraction(i,j) = 0.1
           
            DO k=1,kme
              !if(mynum==1 .and. j==jms .and. i==ims) print*,">>AA=",k,flip(k);call flush(6)
              qc_curr (i,j,k)= max(0.0,mic%rcp(flip(k),i,j))	      ! QC     
              qr_curr (i,j,k)= max(0.0,mic%rrp(flip(k),i,j))	      ! QR   
              qi_curr (i,j,k)= max(0.0,mic%rpp(flip(k),i,j))	      ! QI   
              qs_curr (i,j,k)= max(0.0,mic%rsp(flip(k),i,j))	      ! QS   
              qg_curr (i,j,k)= max(0.0,mic%rgp(flip(k),i,j))	      ! QG   

              rad_cf  (i,j,k)= max(0.0,mic%cldfr(flip(k),i,j))      ! cloud fraction
              
              rcond = qc_curr (i,j,k) + qr_curr (i,j,k) + &
	      	    qi_curr (i,j,k) + qs_curr (i,j,k) + qg_curr (i,j,k)
	     	       
              qv_curr (i,j,k)= max(1.e-12, basic%rtp(flip(k),i,j) - rcond)		      
              
              W       (i,j,k)= basic%wp(flip(k),i,j)	    ! vertical velocity (m/s) 
              U       (i,j,k)= basic%up(flip(k),i,j)	    ! U velocity (m/s) 
              V       (i,j,k)= basic%vp(flip(k),i,j)	    ! V velocity (m/s) 
              
            ! Delta-Z layer thickness (gfdl expects this to be negative)
              DZ      (i,j,k)= - grd%rtgt(i,j)/dzt(flip(k)) ! layer thickness (m) 
           
            ENDDO
            !
            !--special treatment for pressure layers
            DO k=1,kme+2     ! note that kme = mzp-2
              ! Exner function (brams data structure)
              exner  (k,i,j)= ((basic%pp(k,i,j)+basic%pi0(k,i,j))*cpi) 
            ENDDO
            DO k=1,kme      
              ! pressure thickness (Pa)
              dp_tmp(k) = 0.5*(exner(k+2,i,j)** cpor - exner(k,i,j)** cpor) * p00
              !if(mynum==1 .and. i==ims .and. j==jms) print*,"x=",dp_tmp(k) !, press(k+2,i,j) ,press(k,i,j),k
            ENDDO
            m=kme
            DO k=1,kme      
              DP(i,j,k) = - dp_tmp(m) ; m=m-1 !invert press
            ENDDO
            !
            !--special treatment for CCN/ICN numbers conce
            DO k=1,kme             
              NACTLI (i,j,k)= 3.e+2 !cm-2 !???????????????
              !qni_curr(i,j,k)= max(0.0,mic%cpp(flip(k),i,j))	      ! NI  !??????????????? 
              !qnr_curr(i,j,k)= max(0.0,mic%crp(flip(k),i,j))	      ! NR  !???????????????
            ENDDO
            !
            !- get  temperature (temp) from theta_il (thp) and condensates
            DO k=1,kme
               exner_tmp = exner (flip(k),i,j)
               tempK	   = basic%theta(flip(k),i,j)*exner_tmp
         	     til	     = basic%thp  (flip(k),i,j)*exner_tmp
             
               rliq	=  qc_curr(i,j,k) + qr_curr(i,j,k)		  
               rice	=  qi_curr(i,j,k) + qs_curr(i,j,k) + qg_curr(i,j,k)
               qhydm	=  alvl * rliq + alvi * rice
               
               if (tempK .gt. 253.) then
         	       tairstr = 0.5 * (til + sqrt(til * (til + cpi4 * qhydm)))
               else
         	       tairstr = til * (1. + qhydm * cp253i)
               endif
               !- updated  temperature TEMP in Kelvin (adv+dif+rad+conv+)
               TEMP (i,j,k) = tairstr
     	     ENDDO
	 ENDDO; ENDDO
	
  	
   !--- Execute GFDL microphysics
	 call gfdl_cloud_microphys_driver( &
		                 ! input water/cloud species and liquid+ice CCN [NACTL+NACTI]
                     ! RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QS, RAD_QG, RAD_CF, (NACTL+NACTI)/1.e6, &
                       qv_curr, 		  &! QV=qv_curr,     
                       qc_curr, 		  &! QC=qc_curr,     
                       qr_curr, 		  &! QR=qr_curr,     
                       qi_curr, 		  &! QI=qi_curr,     
                       qs_curr, 		  &! QS=qs_curr,     
                       qg_curr, 		  &! QG=qg_curr,     
                       rad_cf,  		  &! cloud fraction  
                       NACTLI,  		  &!@NR=qnr_curr, 
                     ! Output tendencies
                       DQVDT_micro, DQLDT_micro, DQRDT_micro, DQIDT_micro, &
                       DQSDT_micro, DQGDT_micro, DQADT_micro, DTDT_micro, &
                     ! Input fields
                       TEMP, W, U, V, DUDT_micro, DVDT_micro, DZ, DP, &
                     ! constant inputs
                       AREA, dt_moist, FRLAND, CNV_FRACTION, &
		       
                       ANV_ICEFALL, LS_ICEFALL, &
                     ! Output rain re-evaporation and sublimation
                       REV_MC_X, RSU_MC_X, &   ! EVAPC_X, & 
		       
                     ! Output precipitates
                       PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL, &
                     ! Output mass flux during sedimentation (Pa kg/kg)
                       PFL_LS_X(:,:,1:kme), PFI_LS_X(:,:,1:kme), &
                     ! constant grid/time information
                       LHYDROSTATIC, LPHYS_HYDROSTATIC, &
                       ims,ime, jms,jme, 1,kme, 1, kme, &
                       re_cloud, re_ice, re_snow)
                
   !-------------------- print ---------->
	 ! if(mynum==10000 ) then
   !         print*,"-----------------------------------------------------------------------"
   !         print*,"WX=",maxval(V), maxval(W),  maxval(U),mynum
   !         print*,"WN=",minval(V), minval(W),  minval(U),mynum
   !         print*,"TX=",maxval(TEMP), maxval(qv_curr),  maxval(qc_curr),mynum
   !         print*,"TN=",minval(TEMP), minval(qv_curr),  minval(qc_curr),mynum
   !         print*,"PX=",maxval(dp), maxval(dz),  maxval(qc_curr),mynum
   !         print*,"PN=",minval(dp), minval(dz),  minval(qc_curr),mynum
   !         call flush(6)	
	 ! endif
   !-------------------- print ----------<


   !--- update variables after cloud microphysics processes (only in the halo)
   DO j = ja,jz   ; DO i = ia,iz

         !- surface precipitation (units are kg/m^2 = mm and mm/s)
	       !- accpr,pcprr   kg/m2 - rain+ice+snow+graupel ! p = for each dt  (or per time step)
         mic%pcprr(i,j) = dt_moist*(PRCP_RAIN(i,j) + PRCP_SNOW   (i,j) + &
                                    PRCP_ICE (i,j) + PRCP_GRAUPEL(i,j)) / 86400. 

         mic%accpr(i,j) = mic%accpr(i,j) + mic%pcprr(i,j)	  
          
         !- accps,pcprs  kg/m2 - ice+snow
         mic%pcprs(i,j) = dt_moist*(PRCP_SNOW(i,j) + PRCP_ICE(i,j))/ 86400.    
         mic%accps(i,j) = mic%accps(i,j) + mic%pcprs(i,j)
          
         !- accpg, pcprg  &! kg/m2 - graupel
         mic%pcprg(i,j) = dt_moist* PRCP_GRAUPEL(i,j)/ 86400.
         mic%accpg(i,j) = mic%accpg(i,j) + mic%pcprg(i,j)


	      DO k=kme,1,-1
!-srf 	   
!if(i==89 .and. j==22 .and. minval(mic%rcp(:,i,j))< 0.)  then
!      print*,'rrp',k, mic%rcp(flip(k),i,j), qi_curr(i,j,k) , qi_curr(i,j,k) + DQIDT_micro(i,j,k) * DT_MOIST
!      print*,'rv ',k, basic%rv(flip(k),i,j), qv_curr(i,j,k), qv_curr(i,j,k) + DQVDT_micro(i,j,k) * DT_MOIST
!call flush(6)
!endif
!-srf 
	        mic%rcp(flip(k),i,j)= max(0., qc_curr(i,j,k) + DQLDT_micro(i,j,k) * DT_MOIST)
	        mic%rrp(flip(k),i,j)= max(0., qr_curr(i,j,k) + DQRDT_micro(i,j,k) * DT_MOIST)
	        mic%rpp(flip(k),i,j)= max(0., qi_curr(i,j,k) + DQIDT_micro(i,j,k) * DT_MOIST)
	        mic%rsp(flip(k),i,j)= max(0., qs_curr(i,j,k) + DQSDT_micro(i,j,k) * DT_MOIST)
	        mic%rgp(flip(k),i,j)= max(0., qg_curr(i,j,k) + DQGDT_micro(i,j,k) * DT_MOIST)
	   
          mic%cldfr(flip(k),i,j) = max(0., rad_cf(i,j,k) + DQADT_micro(i,j,k) * DT_MOIST)

          !- cloud ice and liq effective radius,  RRTMG requires in micrometer
          mic%rei (flip(k),i,j)= re_ice  (i,j,k)  
          mic%rel (flip(k),i,j)= re_cloud(i,j,k)
	       !mic%snow(flip(k),i,j)= re_snow (i,j,k)  !-- srf check this latter.
    

	        basic%rv(flip(k),i,j)= max(1.e-12, qv_curr(i,j,k) + DQVDT_micro(i,j,k) * DT_MOIST)
	 	 
 	        tempK                =   temp(i,j,k) + DTDT_micro (i,j,k) * DT_MOIST
          
	        basic%theta(flip(k),i,j) = tempK/ exner(flip(k),i,j)

         ENDDO

	      DO k=2,kme
	          !- update liq-ice potential temperature THP in Kelvin including microphysics processes
	          rliq     = mic%rcp(k,i,j) + mic%rrp(k,i,j)   
	          rice     = mic%rpp(k,i,j) + mic%rsp(k,i,j) + mic%rgp(k,i,j)	  
	  
	          tempK    = basic%theta(k,i,j) * exner(k,i,j)
	  
	          basic%thp(k,i,j) = basic%theta(k,i,j)*(1. + alvl * rliq/(cp * max(tempK,253.))  &
 			                       + alvi * rice/(cp * max(tempK,253.)) ) **(-1.0)	
            !- update total water 
	          basic%rtp(k,i,j) = basic%rv(k,i,j) + rliq + rice
         ENDDO
!reserved	   
!	  cpp(k)= qni_curr(1,k,1)
!	  crp(k)= qnr_curr(1,k,1)
!         U1   = U1   + DUDT_micro  * DT_MOIST
!         V1   = V1   + DVDT_micro  * DT_MOIST	   
!         W1 was updated in gfdl_microphsycis, fill WI tendency export
!         WI =  (W1 - W)/DT_MOIST
!         RAD_CF = RAD_CF + DQADT_micro * DT_MOIST
!reserved

	       !- definition for k=1
	       basic%rtp(1,i,j)  = basic%rtp(2,i,j)  
	       mic%rcp  (1,i,j)  = mic%rcp  (2,i,j)  
	       mic%rrp  (1,i,j)  = mic%rrp  (2,i,j)  
	       mic%rpp  (1,i,j)  = mic%rpp  (2,i,j)  
	       mic%rsp  (1,i,j)  = mic%rsp  (2,i,j)  
	       mic%rgp  (1,i,j)  = mic%rgp  (2,i,j)  
         mic%cldfr(1,i,j)  = mic%cldfr(2,i,j)

	       basic%rv   (1,i,j) = basic%rv   (2,i,j)  
	       basic%theta(1,i,j) = basic%theta(2,i,j)
	       basic%thp  (1,i,j) = basic%thp  (2,i,j)
         
         !- cloud ice and liq effective radius 
         mic%rei(1,i,j)= mic%rei(2,i,j) 
         mic%rel(1,i,j)= mic%rel(2,i,j)
        !mic%snow(1,i,j)=mic%snow(2,i,j)

        !temporary 
        !mic%rei   (1:m1,i,j) =5.01E-6 * 1.e+6 ! RRTM requires in micrometer
        !mic%rel   (1:m1,i,j) =2.51E-6 * 1.e+6 ! RRTM requires in micromete


!reserved
!	  mic%cpp (1,i,j)  = mic%cpp (2,i,j)  
!	  mic%crp (1,i,j)  = mic%crp (2,i,j)  
!	  mic%ccp (1,i,j)  = mic%ccp (2,i,j) 

ENDDO;ENDDO

END  SUBROUTINE brams_to_mic_gfdl
!------------------------------------------------------------------------------------

 SUBROUTINE FLIPZ_MINUS1(flip,mzp)
    implicit none
    integer, intent(In) :: mzp
    integer, dimension(mzp), INtent(inout) :: flip
    integer :: m,k
    m=mzp
    do k=1,mzp
     flip(k)=m
     m=m-1
    enddo
    flip(mzp)=1
 END SUBROUTINE FLIPZ_MINUS1
