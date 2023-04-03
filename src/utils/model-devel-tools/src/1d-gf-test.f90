program GF_1d_driver

  USE module_gate
  use module_cu_gf2   , only: GFDRV2
  use module_cu_gf_1d_oldversion  , only: GFDRV1d
  use module_cu_gd_fim, only:  GRELLDRV_FIM
  IMPLICIT NONE

    
  INTEGER,PARAMETER :: mzp=klev ,mxp=1   ,myp=1
  
  !-- ichoice_deep:     0 ensemble, 1 grell, 4 omega , 7 MCONV, 10 = KF, 13 = PB
  INTEGER, PARAMETER :: icoic=0
  
  !-- ichoice_shallow   0 ensemble, 1 Wstar, 4 heat-engine or 7 BLQE 
  INTEGER, PARAMETER :: icoic_sh=1

  INTEGER,PARAMETER :: &
         CCATT =1    , ngrid=1 ,ngrids_cp=1     & 
	 ,iens  =1   , mynum=1 ,npatch=1        &
	 ,i0    =1   , j0=1                     &	 
	 ,ia    =1   , ja=1 ,iz=1 ,jz=1         ! 1-d column
	 
  INTEGER,PARAMETER :: mgmxp=mxp &
                      ,mgmyp=myp &
		      ,mgmzp=mzp

	!- converting WRF setting to BRAMS
  INTEGER,PARAMETER :: &
         ids=1   ,ide=mxp ,jds=1   ,jde=myp ,kds=1, kde=mzp	 &	
	,ims=1   ,ime=mxp ,jms=1   ,jme=myp ,kms=1, kme=mzp	 &	 
	,ips=ia+1,ipe=iz-2,jps=ja+1,jpe=jz-2,kps=1, kpe=mzp	 &	   
	,its=ia  ,ite=iz  ,jts=ja  ,jte=jz  ,kts=1, kte=mzp!-1	 

  
   INTEGER ,PARAMETER ::  &
	   cugd_avedx  =1 &
          ,ishallow_g3 =1 &
	  ,imomentum   =0 &
	  ,training    =0 &
	  ,level=3     
  integer,parameter :: autoconv = 1!2 ! =1, Kessler
                                      ! =2, Berry 
  integer,parameter :: aerovap = 1!3  ! =1, orig
                                    ! =2, mix orig+new
				    ! =3, new 

  integer, parameter ::               &
       maxiens    = 2,  & !Cloud spectral size
       maxens     = 3,  & ! 3  ensemble one on cap_max
       maxens2    = 3,  & ! 3  ensemble two on precip efficiency
       maxens_sh  = 3,  & ! 3  ensemble one on mbdt
       maxens2_sh = 1,  & ! 1  ensemble two on precip efficiency
       maxens3_sh = 10, & !10 ensemble three done in cup_forcing_ens16
       maxens3    = 16, & !16 ensemble three done in cup_forcing_ens16

       maxens_g3d  = 1,  & ! 1  ensemble one on cap_max for G3d
       maxens2_g3d = 1,  & ! 1  ensemble two on precip efficiency for G3d
       maxens3_g3d = 16, & !16 ensemble three done in cup_forcing_ens16 for G3d
       ensdim      = 1*maxens    *maxens2    *maxens3    , &
       ensdim_g3d  = 1*maxens_g3d*maxens2_g3d*maxens3_g3d

   REAL, PARAMETER ::      &
       rgas    = 287.,    &
       cp      = 1004.,   &
       rm      = 461.,    &
       p00     = 1.e5,    &
       tcrit   = 273.15,  &
       g       = 9.80,    &
       cpor    = cp/rgas, &
       xl      = 2.5e6,   &
       akmin   = 1.0,     &
       tkmin   = 1.e-5
   
   LOGICAL,PARAMETER :: do_cupar_mcphys_coupling = .false. ! direct link cupar-microphysics
  
   REAL,  DIMENSION(kms:kme)  ::  zmn,ztn

   REAL,  DIMENSION(kms:kme ,  ims:ime , jms:jme ) :: &
                                                          UP,   &
                                                          VP,   &
                                                          WP,   &
                                                          rv,   &
                                                          rtp,  &
	           				       theta   ,& 
                   				       thp     ,& 
                   				       pp      ,& 
                   				       pi0     ,& 
                                                       dn0     ,&
                                                       tend_pt ,&
		                                       rcp,tkep           

   REAL, DIMENSION(ims:ime , jms:jme,npatch)  :: patch_area


   REAL, DIMENSION( ims:ime , jms:jme ) ::  aot500 , temp2m , rtgt &
                                           ,sflux_r, sflux_t, topt &
					   ,rshort 

   REAL :: DTLT,  time

   REAL, DIMENSION( ims:ime , jms:jme ) ::      &
                          conprr,xmb_deep,                     &
                          apr_gr,apr_w,apr_mc,apr_st,apr_as,    &
                          xmb_shallow,err_deep
   REAL, DIMENSION( ims:ime , jms:jme ) ::        &
                  weight_GR,weight_W,weight_MC,weight_ST,weight_AS
			  
  
! Optionals
!
   REAL, DIMENSION(kms:kme , ims:ime ,  jms:jme )   ::   &  
              THSRC	 & ! temp tendency	       
	     ,RTSRC	 & ! rv tendency	       
	     ,CLSRC	 & ! cloud/ice tendency        
	     ,USRC	 & ! U tendency 	       
	     ,VSRC	 & ! V tendency 	       
	     ,MUP	 & ! updraft mass flux         
	     ,lsfth	 & ! forcing for theta deep    
	     ,lsfrt	 & ! forcing for rv deep       
	     ,lsfth_sh	 & ! forcing for theta shallow 
	     ,lsfrt_sh	   ! forcing for rv shallow    
!
! Flags relating to the optional tendency arrays declared above
! Models that carry the optional tendencies will provdide the
! optional arguments at compile time; these flags all the model
! to determine at run-time whether a particular tracer is in
! use or not.
!
   LOGICAL           ::     F_QV      &
                           ,F_QC      &
                           ,F_QR      &
                           ,F_QI      &
                           ,F_QS
!- for convective transport-start
  integer, dimension(mgmxp,mgmyp,maxiens,ngrids_cp) ::         &
               ierr4d  		     &  	     
	      ,jmin4d  		     & 
	      ,kdet4d  		     & 
	      ,k224d	             & 
	      ,kbcon4d 		     & 
	      ,ktop4d  		     & 
	      ,kpbl4d  		     & 
	      ,kstabi4d		     & 
	      ,kstabm4d		   

   real,dimension(mgmxp,mgmyp,maxiens,ngrids_cp) :: &
 	       xmb4d		     & 
	      ,edt4d		     & 
	      ,pwav4d		     	      
   real,dimension(mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp) ::&
 	       pcup5d 		     & 
              ,up_massentr5d	     &        
	      ,up_massdetr5d	     &
	      ,dd_massentr5d	     &
	      ,dd_massdetr5d	     &
 	      ,zup5d		     &
	      ,zdn5d   		     & 
	      ,prup5d  		     & 
	      ,prdn5d  		     & 
	      ,clwup5d 		     & 
	      ,tup5d   		      
  
  real :: grid_length
!  
!- this for the namelist gf.inp
      !character (len=50) :: runname, runlabel
      namelist /run/ runname, runlabel, version

!- Here are the place for data related with the GATE soundings
!- soundings arrays
      integer ::jk, nruns, version

!- for grads output
      integer :: nrec,nvx,nvar,klevgrads,nvartotal,int_byte_size
      real    :: real_byte_size
!--- allocation      
      allocate(cupout(nvar_grads))
      do nvar=1,nvar_grads
        allocate(cupout(nvar)%varp(klon,klev))
        allocate(cupout(nvar)%varn(2))
        cupout(nvar)%varp(:,:)=0.0
        cupout(nvar)%varn(:)  ="xxxx"
      enddo
      if(.not. use_gate) then
       print*,"====================================================================="
       print*, "use_gate logical flag must be true to run in 1-d, model will stop"
       print*,"====================================================================="
       stop "use_gate flag"
      endif
!
!------------------- simulation begins  ------------------
!
!- reads namelist                
     open(15,file='gf.inp',status='old',form='formatted')	  
     read(15,nml=run)
     close(15)


!- reads gate soundings                
    open(7,file="gate.dat",form="formatted")
     read(7,*)
     do jl=1,klon
     	read(7,*)
     	!z(m)  p(hpa) t(c) q(g/kg) u  v (m/s) w(pa/s) q1 q2 !!!qr (k/d) advt(k/d) advq(1/s)
     	do jk=klev,1,-1
     	read(7,*)pgeo(jl,jk),ppres(jl,jk),ptemp(jl,jk),pq(jl,jk),        &
     		 pu(jl,jk),pv(jl,jk),pvervel(jl,jk), &
     		 zq1(jl,jk),zq2(jl,jk),zqr(jl,jk),zadvt(jl,jk),&
     		 zadvq(jl,jk)			    
     	!print*,"GATE=",jl,jk,pgeo(jl,jk),zadvq(jl,jk)
     	end do
     enddo
    close(7)
!-



!- general  initialization ---------------------------------------
   grid_length = 100000. ! meters
   dtlt=60. !seconds
   time=0.
   
   != nz
   zmn	       = 0.
   ztn	       = 0.
   != nx,ny
   RTGT    (:,:) = 1.     !don´t change this
   aot500  (:,:) = 0.1    ! #
   temp2m  (:,:) = 310.   ! Kelvin
   sflux_r (:,:) = 300./(1.1*xl) !(kg/kg/s)
   sflux_t (:,:) = 200./(1.1*cp) !(K/s)
   
   CONPRR  (:,:) = 0.
   apr_GR  (:,:) = 0.
   apr_W   (:,:) = 0.
   apr_MC  (:,:) = 0.
   apr_ST  (:,:) = 0.
   apr_AS  (:,:) = 0.
   xmb_deep(:,:) = 0.
   err_deep(:,:) = 0.
   xmb_shallow(:,:) = 0.
   weight_GR  (:,:) = 0.
   weight_W   (:,:) = 0.
   weight_MC  (:,:) = 0.
   weight_ST  (:,:) = 0.
   weight_AS  (:,:) = 0.
   topt       (:,:) = 0.
   rshort     (:,:) = 0.
   
   != nx,ny,npatch
   patch_area (:,:,:)= 1.!don´t change this

   != nz,nx,ny
   dn0    (:,:,:)= 1.
   up     (:,:,:)= 1.
   vp     (:,:,:)= 1.
   theta  (:,:,:)= 300. 
   thp    (:,:,:)= 300. 
   pp     (:,:,:)= 1000.
   pi0    (:,:,:)= 0.1
   wp     (:,:,:)= 0.
   rv     (:,:,:)= 0.001
   rtp    (:,:,:)= 0.001
   tend_PT(:,:,:)= 0.

   THSRC  (:,:,:)= 0.
   RTSRC  (:,:,:)= 0.
   CLSRC  (:,:,:)= 0.
   USRC   (:,:,:)= 0.
   VSRC   (:,:,:)= 0.
   MUP	  (:,:,:)= 0.
   lsfth  (:,:,:)= 0.
   lsfrt  (:,:,:)= 0.
   lsfth_sh  (:,:,:)= 0.
   lsfrt_sh  (:,:,:)= 0.
   
   rcp    (:,:,:)= 0.
   tkep   (:,:,:)= tkmin
!- end of  initialization ---------------------------------------


!- big loop on the gate soundings
      
   do jl=1,klon !klon=number of soundings
     write(0,*) "#######################################",jl

     IF(VERSION==1) then ! BRAMS 3-d version
        CALL GFDRV2(       &
               mgmxp,mgmyp,mgmzp,ngrid,ngrids_cp,iens &
              ,mynum,i0,j0,time,mzp,mxp,myp &
              ,dtlt        & !
              ,grid_length & !
              ,autoconv    & !
              ,aerovap     & !
              ,dn0	   & !basic_g(ngrid)%dn0     
	      ,CONPRR      & !cupout_g(ngrid)%CONPRR 
              ,up	   & !basic_g(ngrid)%up      
              ,vp	   & !basic_g(ngrid)%vp      
              ,theta	   & !basic_g(ngrid)%theta   
              ,thp	   & !basic_g(ngrid)%thp     
              ,pp	   & !basic_g(ngrid)%pp      
              ,pi0	   & !basic_g(ngrid)%pi0     
	      ,wp	   & !basic_g(ngrid)%wp      
	      ,rv	   & !basic_g(ngrid)%rv      
	      ,rtp	   & !basic_g(ngrid)%rtp     
              ,RTGT        & !grid_g(ngrid)%RTGT 
              ,tend_PT     & !tend%PT
	      ,XL	   & !
	      ,CP	   & !
	      ,G	   & !
	      ,rm          &
	      ,p00         &
	      ,cpor        & !
	      ,rgas        & !
!
	      ,zmn	   & !zmn(:,ngrid) 
	      ,ztn         & !ztn(:,ngrid) 
              ,apr_GR      & !g3d_ens_g(apr_gr,ngrid)%apr  
              ,apr_W       & !g3d_ens_g(apr_w ,ngrid)%apr  
              ,apr_MC      & !g3d_ens_g(apr_mc,ngrid)%apr  
              ,apr_ST      & !g3d_ens_g(apr_st,ngrid)%apr  
              ,apr_AS      & !g3d_ens_g(apr_as,ngrid)%apr  
              ,xmb_deep    & !g3d_g(ngrid)%xmb_deep   
              ,err_deep    & !g3d_g(ngrid)%err_deep   
              ,xmb_shallow & !g3d_g(ngrid)%xmb_shallow         
	      ,weight_GR   & !g3d_ens_g(apr_gr,ngrid)%weight
              ,weight_W    & !g3d_ens_g(apr_w ,ngrid)%weight
              ,weight_MC   & !g3d_ens_g(apr_mc,ngrid)%weight
              ,weight_ST   & !g3d_ens_g(apr_st,ngrid)%weight
              ,weight_AS   & !g3d_ens_g(apr_as,ngrid)%weight
              ,training    &
!
	      ,topt	   &!grid_g(ngrid)%topt     
              ,patch_area  &!leaf_g(ngrid)%patch_area
	      ,npatch      &
              ,rshort      &!radiate_g(ngrid)%rshort 
!
	      ,cugd_avedx					& 
	      ,imomentum          				& 
              ,ensdim_g3d,maxiens,maxens_g3d,maxens2_g3d,maxens3_g3d,icoic      & 
              ,ishallow_g3                                      & 
	      ,ids,ide, jds,jde, kds,kde                        & 
              ,ims,ime, jms,jme, kms,kme                        & 
              ,ips,ipe, jps,jpe, kps,kpe                        & 
              ,its,ite, jts,jte, kts,kte                        & 
!
!
              ,THSRC	  & ! temp tendency		   g3d_g(ngrid)%THSRC	   
              ,RTSRC	  & ! rv tendency		   g3d_g(ngrid)%RTSRC	   
              ,CLSRC	  & ! cloud/ice tendency	   g3d_g(ngrid)%CLSRC	   
              ,USRC	  & ! U tendency		   g3d_g(ngrid)%USRC	   
              ,VSRC	  & ! V tendency		   g3d_g(ngrid)%VSRC	   
              ,MUP	  & ! updraft mass flux 	   g3d_g(ngrid)%MUP	   
	      ,lsfth	  & ! forcing for theta deep    cuforc_g(ngrid)% lsfth  
	      ,lsfrt	  & ! forcing for rv deep       cuforc_g(ngrid)% lsfrt  
	      ,lsfth_sh   & ! forcing for theta shallow cuforc_sh_g(ngrid)%lsfth
	      ,lsfrt_sh   & ! forcing for rv shallow    cuforc_sh_g(ngrid)%lsfrt
!
              ,level      &
!
	      ,rcp	  & !micro_g(ngrid)%rcp 
	      ,aot500     & ! aot at 500nm
	      ,temp2m     & ! 2m-temp
  	      ,sflux_r    & !turb_g(ngrid)%sflux_r  
              ,sflux_t    & !turb_g(ngrid)%sflux_t  
              ,tkep	  & !turb_g(ngrid)%tkep     
!
              ,TKMIN      &
	      ,akmin      & !akmin(ngrid)
	      ,do_cupar_mcphys_coupling   &

     !- for convective transport-start
              ,ierr4d  		    &		    
	      ,jmin4d  		    & 
	      ,kdet4d  		    & 
	      ,k224d	            & 
	      ,kbcon4d 		    & 
	      ,ktop4d  		    & 
	      ,kpbl4d  		    & 
	      ,kstabi4d		    & 
	      ,kstabm4d		    & 
	      ,xmb4d		    & 
	      ,edt4d		    & 
	      ,pwav4d		    & 
	      ,pcup5d  		    & 
              ,up_massentr5d	    &	     
	      ,up_massdetr5d	    &
	      ,dd_massentr5d	    &
	      ,dd_massdetr5d	    &
	      ,zup5d		    &
	      ,zdn5d   		    & 
	      ,prup5d  		    & 
	      ,prdn5d  		    & 
	      ,clwup5d 		    & 
	      ,tup5d   		    & 
     !- for convective transport- end
          			    )

     ELSEIF(version == 2) then ! version 1-d
        CALL GFDRV1d(       &
               mgmxp,mgmyp,mgmzp,ngrid,ngrids_cp,iens &
              ,mynum,i0,j0,time,mzp,mxp,myp &
              ,dtlt        & !
              ,grid_length & !
              ,autoconv    & !
              ,aerovap     & !
              ,dn0	   & !basic_g(ngrid)%dn0     
	      ,CONPRR      & !cupout_g(ngrid)%CONPRR 
              ,up	   & !basic_g(ngrid)%up      
              ,vp	   & !basic_g(ngrid)%vp      
              ,theta	   & !basic_g(ngrid)%theta   
              ,thp	   & !basic_g(ngrid)%thp     
              ,pp	   & !basic_g(ngrid)%pp      
              ,pi0	   & !basic_g(ngrid)%pi0     
	      ,wp	   & !basic_g(ngrid)%wp      
	      ,rv	   & !basic_g(ngrid)%rv      
	      ,rtp	   & !basic_g(ngrid)%rtp     
              ,RTGT        & !grid_g(ngrid)%RTGT 
              ,tend_PT     & !tend%PT
	      ,XL	   & !
	      ,CP	   & !
	      ,G	   & !
	      ,rm          &
	      ,p00         &
	      ,cpor        & !
	      ,rgas        & !
!
	      ,zmn	   & !zmn(:,ngrid) 
	      ,ztn         & !ztn(:,ngrid) 
              ,apr_GR      & !g3d_ens_g(apr_gr,ngrid)%apr  
              ,apr_W       & !g3d_ens_g(apr_w ,ngrid)%apr  
              ,apr_MC      & !g3d_ens_g(apr_mc,ngrid)%apr  
              ,apr_ST      & !g3d_ens_g(apr_st,ngrid)%apr  
              ,apr_AS      & !g3d_ens_g(apr_as,ngrid)%apr  
              ,xmb_deep    & !g3d_g(ngrid)%xmb_deep   
              ,err_deep    & !g3d_g(ngrid)%err_deep   
              ,xmb_shallow & !g3d_g(ngrid)%xmb_shallow         
	      ,weight_GR   & !g3d_ens_g(apr_gr,ngrid)%weight
              ,weight_W    & !g3d_ens_g(apr_w ,ngrid)%weight
              ,weight_MC   & !g3d_ens_g(apr_mc,ngrid)%weight
              ,weight_ST   & !g3d_ens_g(apr_st,ngrid)%weight
              ,weight_AS   & !g3d_ens_g(apr_as,ngrid)%weight
              ,training    &
!
	      ,topt	   &!grid_g(ngrid)%topt     
              ,patch_area  &!leaf_g(ngrid)%patch_area
	      ,npatch      &
              ,rshort      &!radiate_g(ngrid)%rshort 
!
	      ,cugd_avedx					& 
	      ,imomentum          				& 
              ,ensdim_g3d,maxiens,maxens_g3d,maxens2_g3d,maxens3_g3d,icoic      & 
              ,ishallow_g3                                      & 
	      ,ids,ide, jds,jde, kds,kde                        & 
              ,ims,ime, jms,jme, kms,kme                        & 
              ,ips,ipe, jps,jpe, kps,kpe                        & 
              ,its,ite, jts,jte, kts,kte                        & 
!
!
              ,THSRC	  & ! temp tendency		   g3d_g(ngrid)%THSRC	   
              ,RTSRC	  & ! rv tendency		   g3d_g(ngrid)%RTSRC	   
              ,CLSRC	  & ! cloud/ice tendency	   g3d_g(ngrid)%CLSRC	   
              ,USRC	  & ! U tendency		   g3d_g(ngrid)%USRC	   
              ,VSRC	  & ! V tendency		   g3d_g(ngrid)%VSRC	   
              ,MUP	  & ! updraft mass flux 	   g3d_g(ngrid)%MUP	   
	      ,lsfth	  & ! forcing for theta deep    cuforc_g(ngrid)% lsfth  
	      ,lsfrt	  & ! forcing for rv deep       cuforc_g(ngrid)% lsfrt  
	      ,lsfth_sh      & ! forcing for theta shallow cuforc_sh_g(ngrid)%lsfth
	      ,lsfrt_sh      & ! forcing for rv shallow    cuforc_sh_g(ngrid)%lsfrt
!
              ,level      &
!
	      ,rcp	  & !micro_g(ngrid)%rcp 
	      ,aot500     & ! aot at 500nm
	      ,temp2m     & ! 2m-temp
  	      ,sflux_r    & !turb_g(ngrid)%sflux_r  
              ,sflux_t    & !turb_g(ngrid)%sflux_t  
              ,tkep	  & !turb_g(ngrid)%tkep     
!
              ,TKMIN      &
	      ,akmin      & !akmin(ngrid)

   !- for convective transport-start
              ,ierr4d  		    &		    
	      ,jmin4d  		    & 
	      ,kdet4d  		    & 
	      ,k224d	            & 
	      ,kbcon4d 		    & 
	      ,ktop4d  		    & 
	      ,kpbl4d  		    & 
	      ,kstabi4d		    & 
	      ,kstabm4d		    & 
	      ,xmb4d		    & 
	      ,edt4d		    & 
	      ,pwav4d		    & 
	      ,pcup5d  		    & 
              ,up_massentr5d	    &	     
	      ,up_massdetr5d	    &
	      ,dd_massentr5d	    &
	      ,dd_massdetr5d	    &
	      ,zup5d		    &
	      ,zdn5d   		    & 
	      ,prup5d  		    & 
	      ,prdn5d  		    & 
	      ,clwup5d 		    & 
	      ,tup5d   		    & 
    !- for convective transport- end
          			    )
     ELSEIF(version == 3) then ! version 3-d GF-FIM
        CALL GRELLDRV_FIM(                            &
               mgmxp,mgmyp,mgmzp,ngrid,ngrids_cp,iens &
              ,mynum,i0,j0,time,mzp,mxp,myp &
              ,dtlt         		    & !
              ,grid_length                  & !
              ,autoconv                     & !
              ,aerovap                      & !
              ,dn0	     & !
	      ,CONPRR        & !
              ,up	     & !
              ,vp	     & !
              ,theta	     & !
              ,thp	     & !
              ,pp	     & !
              ,pi0	     & !
	      ,wp	     & !
	      ,rv	     & !
	      ,rtp	     & !
              ,RTGT          & !
              ,tend_PT            & !
!	      ,XL		  & !
!	      ,CP		  & !
!	      ,G		  & !
!	      ,rm                 &
	      ,p00                &
	      ,cpor               & !
	      ,rgas               & !
	      ,zmn		  & !
	      ,ztn		  & !
              ,apr_GR   &
              ,apr_W    &
              ,apr_MC   &
              ,apr_ST   &
              ,apr_AS   &
              ,xmb_deep        &
              ,err_deep        &
              ,xmb_shallow     &
	      ,weight_gr &
              ,weight_w  &
              ,weight_mc &
              ,weight_st &
              ,weight_as &
              ,training  &
	      ,topt		 &
              ,patch_area	 &
	      ,npatch            &
              ,rshort            &
	      ,cugd_avedx					& 
	      ,imomentum          				& 
              ,ensdim_g3d,maxiens,maxens_g3d,maxens2_g3d,maxens3_g3d,icoic      & 
              ,ishallow_g3                                      & 
	      ,ids,ide, jds,jde, kds,kde                        & 
              ,ims,ime, jms,jme, kms,kme                        & 
              ,ips,ipe, jps,jpe, kps,kpe                        & 
              ,its,ite, jts,jte, kts,kte                        & 
              ,THSRC	  & ! temp tendency
              ,RTSRC	  & ! rv tendency
              ,CLSRC	  & ! cloud/ice tendency
              ,USRC	  & ! U tendency
              ,VSRC	  & ! V tendency
              ,MUP	  & ! updraft mass flux
	      ,lsfth	  & ! forcing for theta deep
	      ,lsfrt	  & ! forcing for rv deep 
	      ,lsfth_sh   & ! forcing for theta shallow
	      ,lsfrt_sh   & ! forcing for rv shallow 
              ,level      &
	      ,rcp	  & ! liquid water
	      ,aot500     &! aot at 500nm
	      ,temp2m     &! aot at 500nm
  	      ,sflux_r    &
              ,sflux_t    &
              ,tkep	  &
              ,TKMIN      &
	      ,akmin      &
!- for convective transport-start
              ,ierr4d  		     &  	     
	      ,jmin4d  		     & 
	      ,kdet4d  		     & 
	      ,k224d	             & 
	      ,kbcon4d 		     & 
	      ,ktop4d  		     & 
	      ,kpbl4d  		     & 
	      ,kstabi4d		     & 
	      ,kstabm4d		     & 
	      ,xmb4d		     & 
	      ,edt4d		     & 
	      ,pwav4d		     & 
	      ,pcup5d  		     & 
              ,up_massentr5d	     &        
	      ,up_massdetr5d	     &
	      ,dd_massentr5d	     &
	      ,dd_massdetr5d	     &
	      ,zup5d		     &
	      ,zdn5d   		     & 
	      ,prup5d  		     & 
	      ,prdn5d  		     & 
	      ,clwup5d 		     & 
	      ,tup5d   		     & 
!- for convective transport- end
                                     )


     ENDIF

     !-save for output
     do nvar=1,nvar_grads
        if(cupout(nvar)%varn(1) .ne. 'xxxx') &
	print*,"VAR= ",cupout(nvar)%varn(1)(1:len_trim(cupout(nvar)%varn(1)))&
	              ,"  ",cupout(nvar)%varn(2)(1:len_trim(cupout(nvar)%varn(2)))&
		      ," max-min=", maxval(cupout(nvar)%varp(jl,:)), minval(cupout(nvar)%varp(jl,:))
     enddo
				     
				     
   enddo ! loop over gate soundings				     
   !
   !
   !-- output
   print*,"writing grads control file:',trim(runname)//'.ctl"
   !
   !number of variables to be written
   nvartotal=0
   do nvar=1,nvar_grads
     if(cupout(nvar)%varn(1) .ne. "xxxx") nvartotal=nvartotal+1
   enddo

   !- ctl file
   open(20,file=trim(runname)//'.ctl',status='unknown')
   write(20,2001) '^'//trim(runname)//'.gra'
   write(20,2002) 'undef -9.99e33'
   write(20,2002) 'options' ! zrev'
   write(20,2002) 'title '//trim(runlabel)
   write(20,2003) 1,0.,1. ! units m/km
   write(20,2004) klon,1.,1.
   write(20,2005) klev,(ppres(1,jk),jk=1,klev)
   write(20,2006) 1,'00:00Z01JAN2000','1mn'
   write(20,2007) nvartotal
   do nvar=1,nvar_grads
    if(cupout(nvar)%varn(1) .ne. "xxxx") then
     klevgrads=klev
     write(20,2008) cupout(nvar)%varn(1)(1:len_trim(cupout(nvar)%varn(1)))&
                   ,klevgrads,cupout(nvar)%varn(2)(1:len_trim(cupout(nvar)%varn(2)))
    endif
   enddo
   !do nvx=21,20+nvar2d
   !  klevgrads=0
   !  write(20,2008) gradsname(nvx,1),klevgrads,gradsname(nvx,2)
   !enddo
   write(20,2002) 'endvars'
   close(20)
   !- binary file 
   inquire (iolength=int_byte_size) real_byte_size  ! inquire by output list
   print*, 'opening grads file:',trim(runname)//'.gra'
   open(19,file= trim(runname)//'.gra',form='unformatted',&
           access='direct',status='unknown', recl=int_byte_size*(klon))
   nrec=0
   do nvar=1,nvar_grads
       if(cupout(nvar)%varn(1) .ne. "xxxx") then
        klevgrads=klev
        do jk=1,klevgrads
          nrec=nrec+1
          write(19,REC=nrec) real((cupout(nvar)%varp(:,jk)),4)
        enddo
       endif
   enddo
   !do nvx=21,20+nvar2d
   !   klevgrads=1
   !   do jk=1,klevgrads
   !     nrec=nrec+1
   !     WRITE(19,REC=nrec) real(vargrads(:,jk,nvx),4)
   !   enddo
   !enddo


   close (19)

  2001 format('dset ',a)
  2002 format(a)
  2003 format('xdef ',i4,' linear ',2f15.3)
  2004 format('ydef ',i4,' linear ',2f15.3)
  2005 format('zdef ',i4,' levels ',60f6.0)
  2006 format('tdef ',i4,' linear ',2a15)
  2007 format('vars ',i4)
  2008 format(a10,i4,' 99 ',a40)!'[',a8,']')
  2055 format(60f7.0)
   133 format (1x,F7.0)

END PROGRAM GF_1d_driver

