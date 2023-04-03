!---------------------------------------------------------------------
!---------------------------------------------------------------------
!-srf: convective transport for GF scheme 
!
subroutine trans_conv_mflx_GF(iens,iinqparm)
!------------------------------------------------------------------
!-Convective transport for non-hygroscopic gases and aerosol 
!-Developed by Saulo Freitas (sfreitas@cptec.inpe.br)
!-ref: Freitas, et.al,: Monitoring the transport of biomass burning 
!      emissions in South America. Environmental Fluid Mechanics, 
!      Kluwer Academic Publishers, 2005.
!------------------------------------------------------------------
  use mem_tconv         ,only: stcum1d,dn01d,se,se_cup,sc_up,sc_dn,sc_up_c,sc_dn_c, &
                               henry_coef,pw_up,pw_dn,&
                               trans_conv_alloc, &    !Control alloc variable
                               alloc_trans_conv, &    !alloc subroutine
                               zero_tconv
  !use chem1_list, only: CO
  use mem_chem1         ,only: chem1_g, CHEMISTRY,nspecies=>NSPECIES_TRANSPORTED
  use node_mod          ,only: m1=>mzp,m2=>mxp,m3=>myp,ia,iz,ja,jz,i0,j0,mynum
  use mem_grid          ,only: dtlt,ngrid
 !use mem_scratch       ,only: scratch
 !use mem_basic         ,only: basic_g
 !use mem_cuparm        ,only: cuparm_g
  use mem_grell_param   ,only: mgmxp,mgmyp,mgmzp,maxiens,ngrids_cp
  use mem_scratch1_grell, only: ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d,kpbl4d,   &
                                kstabi4d,kstabm4d,xmb4d,edt4d,pwav4d,               &
                                zup5d, zdn5d, pcup5d,prup5d,prdn5d,clwup5d,tup5d,   &
                                up_massdetr5d, up_massentr5d,                       &
                                dd_massdetr5d ,dd_massentr5d



  implicit none  

  integer, intent(IN) :: iens,iinqparm
  integer :: i,j,k,kr,ipr,jpr,iscl
  !real,intent(INOUT) :: stcum(m1,m2,m3)


  !          se(:,:)  ! environment scalar profile z levels
  !      se_cup(:,:)  ! environment scalar profile z_cup levels
  !       sc_up(:,:)  ! updraft	gas-phase  scalar profile
  !       sc_dn(:,:)  ! downdraft gas-phase  scalar profile
  !     stcum1d(:,:)  ! 1d convective tendency
  !     sc_up_c(:,:)  ! updraft	aqueous-phase scalar profile
  !     sc_dn_c(:,:)  ! downdraft aqueous-phase scalar profile
  !       pw_up(:,:)  ! updraft precitable gas/aer
  !       pw_dn(:,:)  ! downdraft evaporated gas/aer
  !  henry_coef(:,:)  ! henry's constant for gases wet removal
  !       dn01d(:)    ! 1d air density
		      
  if(CHEMISTRY < 0) return 

  if(.not. trans_conv_alloc) then
     call alloc_trans_conv(nspecies,mgmzp) 
     call zero_tconv()
  end if

  if(nspecies == 0) return ! if there is not any specie to be trasnported
                           ! => return
 
  !---- Coordenadas para escrita ascii
  ipr=0 - i0
  jpr=0 - j0
  
  do j=ja,jz
    do i=ia,iz	   
	   
     !do iens=1,maxiens            !loop no espectro de nuvens

      !- check if at (i,j) there is convection or not 
       if(ierr4d(i,j,iens,ngrid) .eq. 0) then       
!tmp
!---------------------------
!edt4d  (i,j,iens,ngrid)=0.0!<<<<<<<<<<
!    if(iens==2) print*,"2shallow",j,i,ierr4d(i,j,iens,ngrid),xmb4d(i,j,iens,ngrid)&
!                ,up_massentr5d(kbcon4d(i,j,iens,ngrid),i,j,iens,ngrid)
!if(mynum==41 .and. i==4  .and. j==6 )then  
!       print*,'=in conv trans======================', i,j,i0,j0
!       print*,   'x=',         xmb4d(i,j,iens,ngrid),   edt4d(i,j,iens,ngrid), &
!                   	       jmin4d(i,j,iens,ngrid),  kdet4d(i,j,iens,ngrid), &
!                   	        k224d(i,j,iens,ngrid), kbcon4d(i,j,iens,ngrid), &
!                   	       ktop4d(i,j,iens,ngrid),  kpbl4d(i,j,iens,ngrid), &
!                   	     kstabi4d(i,j,iens,ngrid),kstabm4d(i,j,iens,ngrid)!, &
!       call flush(6)
!       do k=1,m1-1
!       WRITE (*, 400), k,        zup5d(k,i,j,iens,ngrid),  &
!			  up_massdetr5d(k,i,j,iens,ngrid),  &
!			  up_massentr5d(k,i,j,iens,ngrid),  &
!			          zdn5d(k,i,j,iens,ngrid),  &
!			  dd_massdetr5d(k,i,j,iens,ngrid),  & 
!			  dd_massentr5d(k,i,j,iens,ngrid)!,  &
!       
!       enddo
!       call flush(6)
!  400 FORMAT(I2,6F10.4)	  
!endif   
!---------------------------
! 
! if(mynum==1200 .and. i==8 .and. j == 7) then
!   print*,'xmb4d=',xmb4d(i,j,iens,ngrid), 'clwup5d ' , maxval(clwup5d      (:,i,j,iens,ngrid)),minval(clwup5d       (:,i,j,iens,ngrid)),&
!	     'prup5d '        ,maxval(prup5d	    (:,i,j,iens,ngrid)),minval(prup5d	     (:,i,j,iens,ngrid)),&
!	     'prdn5d '        ,maxval(prdn5d	    (:,i,j,iens,ngrid)),minval(prdn5d	     (:,i,j,iens,ngrid)),&
!	     'pwav4d '        ,pwav4d,      &
!	     'tup5d '	      ,maxval(tup5d	    (:,i,j,iens,ngrid)),minval(tup5d	     (:,i,j,iens,ngrid)),&
!	     'zup5d '	      ,maxval(zup5d	    (:,i,j,iens,ngrid)),minval(zup5d	     (:,i,j,iens,ngrid)),&
!	     'zdn5d '	      ,maxval(zdn5d	    (:,i,j,iens,ngrid)),minval(zdn5d	     (:,i,j,iens,ngrid)),&
!	     'up_massdetr5d ' ,maxval(up_massdetr5d (:,i,j,iens,ngrid)),minval(up_massdetr5d (:,i,j,iens,ngrid)),&
!	     'up_massentr5d ' ,maxval(up_massentr5d (:,i,j,iens,ngrid)),minval(up_massentr5d (:,i,j,iens,ngrid)),&
!	     'dd_massdetr5d ' ,maxval(dd_massdetr5d (:,i,j,iens,ngrid)),minval(dd_massdetr5d (:,i,j,iens,ngrid)),&
!	     'dd_massentr5d ' ,maxval(dd_massentr5d (:,i,j,iens,ngrid)),minval(dd_massentr5d (:,i,j,iens,ngrid)),'mynum', &
!	     mynum, i,j
!   call flush(6)
!
!    print*,	 'pcup5d     ',maxval(pcup5d	    (:,i,j,iens,ngrid)),minval(pcup5d	     (:,i,j,iens,ngrid)), &   
!	  'zup5d	     ',maxval(zup5d	    (:,i,j,iens,ngrid)),minval(zup5d	     (:,i,j,iens,ngrid)), &
!	  'up_massdetr5d     ',maxval(up_massdetr5d (:,i,j,iens,ngrid)),minval(up_massdetr5d (:,i,j,iens,ngrid)), &
!	  'up_massentr5d     ',maxval(up_massentr5d (:,i,j,iens,ngrid)),minval(up_massentr5d (:,i,j,iens,ngrid)), &
!	  'zdn5d	     ',maxval(zdn5d	    (:,i,j,iens,ngrid)),minval(zdn5d	     (:,i,j,iens,ngrid)), &
!	  'dd_massdetr5d     ',maxval(dd_massdetr5d (:,i,j,iens,ngrid)),minval(dd_massdetr5d (:,i,j,iens,ngrid)), & 
!	 'dd_massentr5d     ',maxval(dd_massentr5d (:,i,j,iens,ngrid)),minval(dd_massentr5d (:,i,j,iens,ngrid)),'mynum' ,&
!	    mynum, i,j
!  			 
!   call flush(6)			 
! endif  
!---------------------------

	  stcum1d = 0.
	  
	  call get_se_GF(nspecies,mgmzp,ngrid,m1,m2,m3,i,j,se,se_cup,mynum)

	  call get_incloud_sc_GF(nspecies,mgmzp,m1    ,&
				se,se_cup,sc_up,sc_dn,pw_up,pw_dn ,&
	    			henry_coef             ,  &   
!
	    			k224d  (i,j,iens,ngrid),  &
            			kbcon4d(i,j,iens,ngrid),  &
	    			ktop4d (i,j,iens,ngrid),  &
	    			jmin4d (i,j,iens,ngrid),  &
	    			kdet4d (i,j,iens,ngrid),  &
!
	    			edt4d  (i,j,iens,ngrid),  &
!
!	    		      zcup5d (1,i,j,iens,ngrid),  &
!
	    		      clwup5d      (1,i,j,iens,ngrid),  &
	    		      prup5d       (1,i,j,iens,ngrid),  &
	    		      prdn5d       (1,i,j,iens,ngrid),  &
	    		      pwav4d         (i,j,iens,ngrid),  &
            		      tup5d        (1,i,j,iens,ngrid),  &
	    		      zup5d        (1,i,j,iens,ngrid),  &
	    		      zdn5d        (1,i,j,iens,ngrid),  &
	    		      up_massdetr5d(1,i,j,iens,ngrid),  &
	    		      up_massentr5d(1,i,j,iens,ngrid),  &
	    		      dd_massdetr5d(1,i,j,iens,ngrid),  & 
	    		      dd_massentr5d(1,i,j,iens,ngrid)	&
	    		      )


          call get_stcum_GF(nspecies,mgmzp,m1,stcum1d,se,se_cup,sc_up,sc_dn,sc_up_c,	&
        		    xmb4d(i,j,iens,ngrid),   edt4d(i,j,iens,ngrid), &
        		   jmin4d(i,j,iens,ngrid),  kdet4d(i,j,iens,ngrid), &
        		    k224d(i,j,iens,ngrid), kbcon4d(i,j,iens,ngrid), &
        		   ktop4d(i,j,iens,ngrid),  kpbl4d(i,j,iens,ngrid), &
        		 kstabi4d(i,j,iens,ngrid),kstabm4d(i,j,iens,ngrid), &
!
!       		 zcup5d(1,i,j,iens,ngrid), &
!
        		 pcup5d        (1,i,j,iens,ngrid), &   
        		 zup5d         (1,i,j,iens,ngrid), &
	    		 up_massdetr5d (1,i,j,iens,ngrid), &
	    		 up_massentr5d (1,i,j,iens,ngrid), &
	    		 zdn5d         (1,i,j,iens,ngrid), &
	    		 dd_massdetr5d (1,i,j,iens,ngrid), & 
	    		 dd_massentr5d (1,i,j,iens,ngrid), &
			 mynum,i,j,iinqparm )

	 call apply_wet_removal_GF(nspecies,mgmzp,m1,dtlt,se,pw_up,pw_dn, &
	                             sc_up,sc_dn, stcum1d  , &
	    			     xmb4d (i,j,iens,ngrid), & 
	    			     edt4d (i,j,iens,ngrid), &
	    			     ktop4d(i,j,iens,ngrid), &
	    			   pcup5d(1,i,j,iens,ngrid), &
	    			   i,j,ngrid,mynum           )
     
         call apply_conv_tend_GF(nspecies,mgmzp,m1,m2,m3,stcum1d,i,j,ngrid,ktop4d(i,j,iens,ngrid))

!-srf- testing if negative mixing ratios are created by the conv transport
!        call update_gf(nspecies,mgmzp,ngrid,m1,m2,m3,i,j,se,se_cup,mynum,stcum1d,dtlt)
!-
       endif       
    enddo     !loop  i
  enddo        !loop  j
 
end subroutine trans_conv_mflx_GF
 
!--------------------------------------------------

subroutine get_se_GF(nspecies,mgmzp,ngrid,n1,n2,n3,i,j,se,se_cup,mynum) 
!  use chem1_list, only: spc_alloc_chem  =>spc_alloc,transport,on,off
  use mem_chem1,  only: chem1_g
  use mem_aer1,  only: aer1_g,aer2_g,aer1_inorg_g,aerosol
  use aer1_list, only: nmodes
 
  use mem_tconv, only : nchem_a,nchem_z, ind_chem &
                       ,naer_a,naer_z,ind_aer, ind_mode&
		       ,ind_mode_number,naer_a_number,naer_z_number&
		       ,ind_aer_inorg  ,naer_a_inorg ,naer_z_inorg
  
  implicit none
  integer :: nspecies,ngrid,mgmzp, n1, n2, n3, i, j,mynum
  real    :: se(nspecies,mgmzp), se_cup(nspecies,mgmzp)
  !local Variables
  integer :: k, kr,ispc

  do k=2,n1-1
     kr = k+1   ! level k of conv param corresponds to level K + 1 of rams.

     !-- chemistry section
     do ispc = nchem_a,nchem_z
       
       !- mixing ratio on Z (model) levels
       se    (ispc,k) =      chem1_g(ind_chem(ispc),ngrid)%sc_p(kr,i,j)         
       
       !- mixing ratio on Z_CUP (cloud) levels
       se_cup(ispc,k) = .5*( chem1_g(ind_chem(ispc),ngrid)%sc_p(kr-1,i,j) + &
                             chem1_g(ind_chem(ispc),ngrid)%sc_p(kr  ,i,j)   )  
     
       !if(se_cup(ispc,k) > 0)print*,'st=',ispc,k,se_cup(ispc,k),chem1_g(ind_chem(ispc),ngrid)%sc_p(kr,i,j)

!     do ispc=1,nspecies    
!       if(spc_alloc_chem(transport,ispc) == off) cycle
!       se    (ispc,k) =      chem1_g(ispc,ngrid)%sc_p(kr,i,j)  ! mixing ratio on Z levels
!       se_cup(ispc,k) = .5*( chem1_g(ispc,ngrid)%sc_p(kr-1,i,j) + chem1_g(ispc,ngrid)%sc_p(kr,i,j) )   ! mixing ratio on Z_CUP levels
!     enddo
  
     enddo
     
!     !aerosol section mass 
!     do ispc = naer_a,naer_z
!       se    (ispc,k) =      aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_p(kr,i,j)         ! mixing ratio on Z levels
!
!       se_cup(ispc,k) = .5*( aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_p(kr-1,i,j) + &
!                             aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_p(kr  ,i,j)   )   ! mixing ratio on Z_CUP levels
!     enddo
!     
  enddo

  do ispc = nchem_a,nchem_z
     
      se    (ispc,1) = chem1_g(ind_chem(ispc),ngrid)%sc_p(2,i,j)
      se_cup(ispc,1) = chem1_g(ind_chem(ispc),ngrid)%sc_p(2,i,j)

!  do ispc=1,nspecies 
!    if(spc_alloc_chem(transport,ispc) == off) cycle   
!    se    (ispc,1) = chem1_g(ispc,ngrid)%sc_p(2,i,j)
!    se_cup(ispc,1) = chem1_g(ispc,ngrid)%sc_p(2,i,j)
  enddo
!
!- aer section
!  do ispc = naer_a,naer_z
!     
!      se    (ispc,1) = aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_p(2,i,j)
!      se_cup(ispc,1) = aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_p(2,i,j)!
!
!  enddo

  if(aerosol >= 1) then
    
    do k=2,n1-1
     kr = k+1   ! level k of conv param corresponds to level K + 1 of rams.
     !- aerosol section mass concentration
     do ispc = naer_a,naer_z
            
         se    (ispc,k) =      aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_p(kr,i,j)         ! mixing ratio on Z levels

         se_cup(ispc,k) = .5*( aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_p(kr-1,i,j) + &
                               aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_p(kr  ,i,j)   )   ! mixing ratio on Z_CUP levels
     enddo
    enddo
    !- level 1
    do ispc = naer_a,naer_z
     
        se    (ispc,1) = aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_p(2,i,j)
        se_cup(ispc,1) = aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_p(2,i,j)
    enddo

  
    if(aerosol == 2) then
      !- aerosol section number concentration 
      do k=2,n1-1
       kr = k+1  
       
       do ispc = naer_a_number,naer_z_number
    	      
    	   se	 (ispc,k) =	 aer2_g(ind_mode_number(ispc),ngrid)%sc_p(kr,i,j)	 

    	   se_cup(ispc,k) = .5*( aer2_g(ind_mode_number(ispc),ngrid)%sc_p(kr-1,i,j) + &
    				 aer2_g(ind_mode_number(ispc),ngrid)%sc_p(kr  ,i,j)   )  
       enddo
      enddo
      do ispc = naer_a_number,naer_z_number
        se    (ispc,1) = aer2_g(ind_mode_number(ispc),ngrid)%sc_p(2,i,j)
        se_cup(ispc,1) = aer2_g(ind_mode_number(ispc),ngrid)%sc_p(2,i,j)
      enddo
  
  
       do k=2,n1-1
       kr = k+1  
       !- aerosol section mass of inorganics
       do ispc = naer_a_inorg,naer_z_inorg
    	      
    	   se	 (ispc,k) =	 aer1_inorg_g(ind_aer_inorg(ispc),ngrid)%sc_p(kr,i,j)	 

    	   se_cup(ispc,k) = .5*( aer1_inorg_g(ind_aer_inorg(ispc),ngrid)%sc_p(kr-1,i,j) + &
    				 aer1_inorg_g(ind_aer_inorg(ispc),ngrid)%sc_p(kr  ,i,j)   )  
       enddo
      enddo
      do ispc = naer_a_inorg,naer_z_inorg
        se    (ispc,1) = aer1_inorg_g(ind_aer_inorg(ispc),ngrid)%sc_p(2,i,j)
        se_cup(ispc,1) = aer1_inorg_g(ind_aer_inorg(ispc),ngrid)%sc_p(2,i,j)
      enddo
   endif
  endif
 end subroutine get_se_GF
!---------------------------------------------------------------------------

subroutine get_incloud_sc_GF(nspecies,mgmzp,n1 , &
                 se,se_cup,sc_up,sc_dn, pw_up,pw_dn        , &
		 henry_coef                            , &
		 k22,kbcon,ktop                            , &
		 jmin, kdet,edt                            , &
!
!		 z_cup                                     , &
!
                 cupclw,pw,pwd,pwav                        , &
		 tup,zu,zd                                 , &
                 up_massdetr, up_massentr                  , &
                 dd_massdetr, dd_massentr    		     &
			    )
  use chem1_list, only : CO!, H2O2
  use mem_tconv, only : nchem_a,nchem_z, ind_chem &
                       ,naer_a,naer_z,ind_aer, ind_mode

  implicit none
  integer, intent(in) :: nspecies, mgmzp,n1,k22,kbcon,ktop, jmin, kdet
  
 ! integer, intent(in)           :: spc_number
 ! character (len=*), intent(in) :: spc_name
  real, intent(in) ::  pwav, edt
  
  real, intent(in) , dimension(mgmzp) ::& !z_cup,  
         tup ,       & ! local temperature in cloud updraft [K]
       cupclw,       & ! up cloud liquid water [kg(water)/kg(air)] from cum. scheme
           pw,       & ! 
          pwd,       & ! u
           zu,       & ! norm mass flux updraft
           zd          ! norm mass flux downdraft
                   
  real, intent(in) , dimension(mgmzp) :: up_massdetr, up_massentr,&
                                         dd_massdetr, dd_massentr
                                     
  real, intent(in), dimension(nspecies,mgmzp) :: &
       se,         &
       se_cup      ! gas-phase in environment mixing ratio [kg(gas phase)/kg(air)]

  real, intent(out), dimension(nspecies,mgmzp) :: &
       henry_coef, & ! Henry's constant [(kg(aq)/kg(water))/(kg(gas phase)/kg(air))]
       sc_up ,     & ! gas-phase in updraft     mixing ratio [kg(gas phase)/kg(air)]
       sc_dn ,     & ! gas-phase in updraft     mixing ratio [kg(gas phase)/kg(air)]
!       sc_up_c,    & ! aqueous-phase in updraft mixing ratio [kg(aq)/kg(air)]
       pw_up,      & ! precitable gas/aer amount in updradt [kg[..]/kg[air]
       pw_dn         ! precitable gas/aer amount in downdradt[kg[..]/kg[air]

  !-local var
  integer ispc,k
  !real dz,evaporate,pwdper
  real evaporate,pwdper
  real, dimension(nspecies) ::  conc_equi,conc_mxr,qrch
  real, parameter :: scav_eff = 0.6  ! for smoke : Chuang et al. (1992) J. Atmos. Sci.


  !--- Set zero pw_up, henry_coef and sc_up_c arrays
  !sc_up_c    = 0.
  pw_up      = 0.
  pw_dn      = 0.
  henry_coef = 0.
  sc_up      = 0.
  sc_dn      = 0.
  
  !--- Get Henry's law constants - unit [(kg(aq)/kg(water))/(kg(gas phase)/kg(air))]
  CALL henry(nspecies,mgmzp,n1,kbcon,ktop,henry_coef,tup) 

  !--- updraft section -------------------------------------------------------------------------
  !
  !--- scalar concentration in-cloud updraft - gas phase 
  do k=1,k22-1
       do ispc=1,nspecies
              sc_up(ispc,k) = se_cup(ispc,k)
       enddo
  enddo

  do k=k22,kbcon
       do ispc=1,nspecies
              sc_up(ispc,k) = se_cup(ispc,k22)
       enddo
  enddo

  do k=kbcon+1,ktop
       !dz=z_cup(k)-z_cup(k-1)
       !
       !---- steady state plume equation, for what could
       !---- be in cloud before anything happens (kg/kg*1.e9 = ppbm)
       !
       do ispc=1,nspecies
        sc_up(ispc,k)=  (              zu(k-1)*sc_up(ispc,k-1)  - &
	                  0.5*up_massdetr(k-1)*sc_up(ispc,k-1)  + &
                              up_massentr(k-1)*   se(ispc,k-1)  )     /  &
                         (zu(k-1)-0.5*up_massdetr(k-1)+up_massentr(k-1)+1.e-8)
       enddo
       
       
      !- chemistry section
      !
       do ispc = nchem_a,nchem_z
            !--- equilibrium tracer concentration - Henry's law
	    !--- cloud liquid water tracer concentration
            conc_mxr(ispc) = (henry_coef(ispc,k)*cupclw(k) /(1.+henry_coef(ispc,k)*cupclw(k)) )* sc_up(ispc,k) 

            !---   aqueous-phase concentration in rain water            
	    pw_up(ispc,k) = conc_mxr(ispc)*pw(k)/(1.e-8+cupclw(k))
	    
	    !--- total mixing ratio in gas and aqueous (in cloud) phases
	    sc_up(ispc,k) = sc_up(ispc,k) - pw_up(ispc,k)
	  
       enddo

      !- aerosol section (including now, mass,number and inorg)
      !do ispc = naer_a,naer_z       
       do ispc = naer_a,nspecies
     
            conc_mxr(ispc)  =  scav_eff* sc_up(ispc,k) !unit [kg(aq)/kg(air)]  for aerosol/smoke 
            !
            !---   aqueous-phase concentration in rain water            
	    pw_up(ispc,k) = conc_mxr(ispc)*pw(k)/(1.e-8+cupclw(k))
	    
	    !--- total mixing ratio in gas and aqueous (in cloud) phases
	    sc_up(ispc,k) = sc_up(ispc,k) - pw_up(ispc,k)
	  
       enddo
       
  enddo
  !-----

  do k=ktop+1,n1-1
       do ispc = 1,nspecies
          sc_up(ispc,k) = se_cup(ispc,k)
       enddo
  enddo
  !
  !  !----- get back the in-cloud updraft gas-phase mixing ratio : sc_up(ispc,k)
  !  !      
  !  do k=kbcon+1,ktop
  !       do ispc = 1,nspecies
  !          sc_up(ispc,k) = sc_up(ispc,k) - sc_up_c(ispc,k)
  !       enddo
  !  enddo


  !-----------------------------------------------------------------------------------------------
  !--- downdraft section -------------------------------------------------------------------------
  !
  if(jmin == 0) return 
  
  !--- scalar concentration in-cloud - downdraft
  do k=jmin+1,n1-1
     do ispc=1,nspecies
       sc_dn(ispc,k) = se_cup(ispc,k)
     enddo
  enddo

  !--- at k=jmim
  do ispc=1,nspecies
       sc_dn(ispc,jmin) = se_cup(ispc,jmin)
  enddo

  !
  !--- calculate downdraft mass terms
  do k=jmin-1,1,-1
     !dz=z_cup(k+1)-z_cup(k)
     
     do ispc=1,nspecies
     

          sc_dn(ispc,k) = (            zd(k+1)*sc_dn(ispc,k+1)     -  &
                          0.5*dd_massdetr(k  )*sc_dn(ispc,k+1)     +  &
                              dd_massentr(k  )*se   (ispc,k  )   ) /  &
                        (zd(k+1)-0.5*dd_massdetr(k)+dd_massentr(k)    )

     enddo
  enddo
!--- calculate scalar/moisture properties of downdraft
  do k=jmin-1,1,-1
     !dz       = z_cup(k+1)-z_cup(k)
     pwdper   = -edt* pwd(k)/pwav
     do ispc=1,nspecies
        
        !- amount evaporated by the downdraft and returned to the atmosphere  
        !- in gas phase in downdraft 
        evaporate= pwdper *pw_up(ispc,k)
        !
        pw_dn(ispc,k) = pw_dn(ispc,k) + evaporate
        !
        sc_dn(ispc,k) = sc_dn(ispc,k) + pw_dn(ispc,k)
     enddo
  enddo


end subroutine get_incloud_sc_GF
!------------------------------------------------------------------------

subroutine get_stcum_GF(nspecies,mgmzp,n1,stcum1d,se,se_cup,sc_up,sc_dn,sc_up_c, &
                     xmb,edt                                     , &
     	 	     jmin,kdet,k22,kbcon,ktop,kpbl,kstabi,kstabm , &
!     	 	     z_cup,
		     p_cup                                       , &
		     zu, up_massdetr,up_massentr                 , &
		     zd, dd_massdetr,dd_massentr                 , &
		     mynum,i,j,iinqparm)

  use Phys_const, only: g

  implicit none

  ! Arguments
  integer, intent(in) :: nspecies,mgmzp, n1, jmin, kdet, k22, kbcon, ktop, kpbl, kstabi &
                       , kstabm,mynum,i,j,iinqparm
  real   , intent(in)                :: xmb,edt
  real, dimension(nspecies,mgmzp),intent(in) :: se,se_cup,sc_up,sc_dn,sc_up_c
  real, intent(in), dimension(mgmzp) ::  p_cup !z_cup

  real, dimension(mgmzp),intent(in) :: dd_massdetr,dd_massentr,zu,zd, &
                                       up_massdetr,up_massentr

  real, dimension(nspecies,mgmzp),intent(out) :: stcum1d

  ! Local variables
  integer                :: k,ispc
  real                   :: dz, dp, entup, detup, entdoj, entupk, detupk, detdo, entdo,   &
                            subin, subdown, detdo1, detdo2
  real :: mass

  real, dimension(nspecies) :: dummy

!
!------------  cloud level ktop
!
!- - - - - - - model level ktop-1
!    .   .
!    .   .
!    .   .
!
!------------  cloud level k+2
!
!- - - - - - - model level k+1
!
!------------  cloud level k+1
!
!- - - - - - - model level k
!
!------------  cloud level k
!
!    .   .
!    .   .
!
!------------  cloud level 3 - Zu,Zd, se_cup, p_cup, z_cup, sc_up, sc_dn
!
!- - - - - - - model level 2 - up_massentr, up_massdetr, se
!
!------------  cloud level 2 - Zu,Zd, se_cup, p_cup, z_cup, sc_up, sc_dn
!
!- - - - - - - model level 1 - up_massentr, up_massdetr , se


  do k=2,ktop
     !dz =   z_cup(k+1) - z_cup(k)
     !-- layer thickness in terms of pressure
      dp =  100.*( p_cup(k)   - p_cup(k+1) )
     !definicao em funcao de dn0 (densidade do ar do estado basico)
     ! o sinal - e' contabilizado na expressao para o ds/dt
     ! dp      = g*dn0(k)*dz  !dn0 aqui ja' esta' na grade DO GRELL.
     !
     !-- these three are only used at or near mass detrainment and/or entrainment levels
     entdoj  = 0.
     entupk  = 0.
     detupk  = 0.
     !-- detrainment and entrainment for downdraft
     detdo   = edt*dd_massdetr(k)
     entdo   = edt*dd_massentr(k)
     !-- entrainment/detrainment for updraft
     detup = up_massdetr(k)
     entup = up_massentr(k)
     !-- subsidence 
     subin   = zu(k+1) - edt*zd(k+1)
     subdown = zu(k  ) - edt*zd(k  )
     !-- special levels
     dummy = 0.
     if(k.eq.jmin .and. iinqparm .ne. 6)  then 
        entdoj    = edt*zd(k)
        dummy(:)  = se_cup(:,jmin) 
     endif
     
     if(k.eq.ktop) then
        detupk  = zu(ktop)
        subin   = 0.
        detdo   = 0.
        entdo   = 0.
        entup   = 0.
        detup   = 0.
     endif
     !-- tendency due cumulus transport ( k>=2) 
     do ispc=1,nspecies
          
        stcum1d(ispc,k) = xmb*(		                                  & 
                               detup*0.5*(sc_up(ispc,k+1)+ sc_up(ispc,k)) &
                         +     detdo*0.5*(sc_dn(ispc,k+1)+ sc_dn(ispc,k)) &
			 -     entup    *    se(ispc,k  )                 &
			 -     entdo    *    se(ispc,k  )                 &    
			 +     subin    *se_cup(ispc,k+1)  &
                    	 -     subdown  *se_cup(ispc,k  )  &
			 +     detupk   *sc_up (ispc,ktop) &
                         -     entupk   *se_cup(ispc,k22)  &				 			 
			 -     entdoj   * dummy(ispc)      &
                              )*g/dp	       
		       
         !if(stcum1d(ispc,k) > 0)print*,'st=',ispc,k,stcum1d(ispc,k)
      end do

  end do
  !
  !---tendency due cumulus transport (at bottom, k = 1)
  !dz        =       z_cup(2)-z_cup(1)
  !
  dp        = 100.*(p_cup(1)-p_cup(2))
  !definicao em funcao de dn0 (densidade do ar do estado basico)
  ! o sinal - e' contabilizado na expressao para o ds/dt
  ! dp =   g*dn0(2)*dz ! nivel k da grade do grell corresponde ao nivel k + 1 do rams
  !dp =   g*dn0(1)*dz ! ja' esta' na grade do rams

  !-- old way
  !  detdo1    = edt*zd(2)*cdd(1)*dz
  !  detdo2    = edt*zd(1)
  !  entdo     = edt*zd(2)*entrd(1)*dz
  !  subin     =-edt*zd(2)
  !  do ispc=1,nspecies
  !       stcum1d(ispc,1) = stcum1d(ispc,1) +		      &
  !   		    xmb*(			      &
  !		     detdo1*(sc_dn(ispc,1)+sc_dn(ispc,2))*.5	 &
  !		     + detdo2*    sc_dn(ispc,1)        &
  !		     + subin *    se_cup(ispc,2)       &
  !		     - entdo *    se(ispc,1)	       &
  !		     )*g/dp
  !
  !  enddo
  !- new formulation
   do ispc=1,nspecies
 	stcum1d(ispc,1) = xmb*			       &
 		          edt*zd(2)*( sc_dn(ispc,2) - se_cup(ispc,2) )*&
			  g/dp     
   enddo
  
  return
end subroutine get_stcum_GF

!---------------------------------------------------------------------
subroutine apply_wet_removal_GF(nspecies,mgmzp,m1,dt,se,pw_up,pw_dn, &
                               sc_up,sc_dn, &
                               stcum1d    , &
                               xmb,edt,ktop,p_cup,i,j,ngrid,mynum)
use mem_chem1 ,ONLY: chem1_g
use mem_aer1  ,ONLY: aer1_g,aer2_g,aer1_inorg_g,aerosol
use mem_tconv, only : nchem_a,nchem_z, ind_chem &
                     ,naer_a,naer_z,ind_aer, ind_mode&
		     ,ind_mode_number,naer_a_number,naer_z_number&
		     ,ind_aer_inorg  ,naer_a_inorg ,naer_z_inorg

use chem1_list, only : chem_spc_alloc =>spc_alloc, off, wdp
use aer1_list , only : aer_spc_alloc  =>spc_alloc     &
                      ,aer_numb_alloc =>numb_alloc&
		      ,aer_inorg_alloc=>inorg_alloc
use Phys_const, only: g

  implicit none
  integer, intent(in) :: nspecies,mgmzp,m1,ktop,i,j,ngrid,mynum
  real   , intent(in) :: dt,xmb,edt
  real, intent(in),dimension(mgmzp) :: p_cup
  real, intent(in),dimension(nspecies,mgmzp) :: pw_up,pw_dn,stcum1d,se,sc_up,sc_dn
  real :: mass(nspecies), massi(nspecies),dp
  integer :: k,ispc, nrec
  !---  wetdep units: kg m^-2
  !---  wet depositon mass is accumulated over the time
  do  k=1,ktop+1

    !- chemistry section
    do ispc = nchem_a,nchem_z

       if(chem_spc_alloc(wdp,ind_chem(ispc)) == off) cycle
       
       chem1_g(ind_chem(ispc),ngrid)%sc_wd(i,j) = &
       chem1_g(ind_chem(ispc),ngrid)%sc_wd(i,j) + dt*(pw_up(ispc,k)-pw_dn(ispc,k))*xmb       
    enddo
    
!    !- aerosol section
!    do ispc = naer_a,naer_z
!      
!      if(aer_spc_alloc(wdp,ind_mode(ispc),ind_aer(ispc)) == off) cycle
!
!       aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_wd(i,j) = &
!       aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_wd(i,j) + dt*(pw_up(ispc,k)-pw_dn(ispc,k))*xmb 
!
!      enddo

  !   print *,'pw_up,pw_dn,WM:',k,pw_up(ispc,k),pw_dn(ispc,k),wetdep
  enddo
  
  IF(AEROSOL >= 1) THEN
   do  k=1,ktop+1
    !- aerosol section mass concentration
    do ispc = naer_a,naer_z
      
      if(aer_spc_alloc(wdp,ind_mode(ispc),ind_aer(ispc)) == off) cycle

       aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_wd(i,j) = &
       aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_wd(i,j) + dt*(pw_up(ispc,k)-pw_dn(ispc,k))*xmb 

      enddo
      if(AEROSOL == 2) THEN
         !- aerosol section number concentration
         do ispc = naer_a_number,naer_z_number
           if(aer_numb_alloc(wdp,ind_mode_number(ispc)) == off) cycle
           
	   aer2_g(ind_mode_number(ispc),ngrid)%sc_wd(i,j) = &
           aer2_g(ind_mode_number(ispc),ngrid)%sc_wd(i,j) + dt*(pw_up(ispc,k)-pw_dn(ispc,k))*xmb 
         enddo
         
         !- aerosol section mass inorg
	 do ispc = naer_a_inorg,naer_z_inorg
           if(aer_inorg_alloc(wdp,ind_aer_inorg(ispc)) == off) cycle
           
	   aer1_inorg_g(ind_aer_inorg(ispc),ngrid)%sc_wd(i,j) = &
           aer1_inorg_g(ind_aer_inorg(ispc),ngrid)%sc_wd(i,j) + dt*(pw_up(ispc,k)-pw_dn(ispc,k))*xmb 
         enddo
 
      endif
   enddo
  ENDIF
  return
  
  !-- check mass conservation error
  mass =0.
  massi=0.
  DO k=ktop+1,1,-1
  
      dp=100.*(p_cup(k)-p_cup(k+1))
      do ispc = 1,nspecies
        !-- initial mass in this column
        massi(ispc)= massi(ispc)+ se(ispc,k) *dp/g
        !-- final mass in this column after transport and rain out
        mass (ispc) =mass(ispc) + se(ispc,k)*dp/g  &  
                        + dt*stcum1d(ispc,k)*dp/g  &	   ! transport term (stcum1d already includes 'xmb')
        		+ dt*(pw_up(ispc,k)-pw_dn(ispc,k))*xmb ! wet removal term
       if(mynum==41 .and. i==4 .and. j==6 ) then
         if(ispc==12)  print*,'co',k,se(ispc,k),stcum1d(ispc,k),pw_up(ispc,k),pw_dn(ispc,k)
         call flush(6)
       endif
     enddo
  ENDDO 
  !- print the error (%)
  do ispc = 1,nspecies
     if(mynum==41 .and. i==4 .and. j==6 .and. ispc == 12) then     
        print*,'error CO %',100.*(mass(ispc)-massi(ispc))/(massi(ispc)+1.e-8)
	call flush(6)
     endif	
  enddo
  
!  if(mynum==41 .and. i==4 .and. j==6) then     
!    nrec=0
!    OPEN(19,file='gf.gra',         &
!     form='unformatted',access='direct',status='unknown',  &
!     recl=4*m1)  !PC 
     
!     nrec=nrec+1
!     WRITE(19,rec=nrec) real(se(12,1:m1),4)
!     nrec=nrec+1
!     WRITE(19,rec=nrec) real(stcum1d(12,1:m1),4)
!     nrec=nrec+1
!     WRITE(19,rec=nrec) real(pw_up(12,1:m1),4)
!     nrec=nrec+1
!     WRITE(19,rec=nrec) real(pw_dn(12,1:m1),4)
!     nrec=nrec+1
!     WRITE(19,rec=nrec) real(sc_dn(12,1:m1),4)
!     nrec=nrec+1
!     WRITE(19,rec=nrec) real(sc_up(12,1:m1),4)
!    close(19)
!    stop 333
!   endif
  
  

end subroutine apply_wet_removal_GF
!---------------------------------------------------------------------

subroutine apply_conv_tend_GF(nspecies,mgmzp,m1,m2,m3,stcum1d,i,j,ngrid,ktop)
use mem_chem1 ,ONLY: chem1_g
use mem_aer1  ,ONLY: aer1_g,aer2_g,aer1_inorg_g,aerosol
use mem_tconv, only : nchem_a,nchem_z, ind_chem &
                     ,naer_a,naer_z,ind_aer, ind_mode&
		     ,ind_mode_number,naer_a_number,naer_z_number&
		     ,ind_aer_inorg  ,naer_a_inorg ,naer_z_inorg
implicit none
integer, intent(IN)::nspecies,mgmzp,m1,m2,m3,i,j,ngrid,ktop
real, dimension(nspecies,mgmzp) :: stcum1d

!  local variables
integer ij,kij,ispc,k,kr
! the expression: kij= k +  m1*(i-1) + (m1*(m2-1)+m1)*(j-1)
! provides the memory position of (k,i,j) at tendency array (1-d)

 ij=  m1*(i-1) + (m1*(m2-1)+m1)*(j-1)   
 
! do k=1,m1-1
 do k=1,ktop
    kr= k + 1   !nivel k da grade do grell corresponde ao nivel k + 1 do rams
    
    kij = kr + ij
    
    ! chemistry section
    do ispc = nchem_a,nchem_z

       ! old way : stcum(kr,i,j)=stcum1d(k)
       chem1_g(ind_chem(ispc),ngrid)%sc_t(kij)  = chem1_g(ind_chem(ispc),ngrid)%sc_t(kij) +&
                                                  stcum1d(ispc,k)
       !call xxxxxx(mgmzp,m1,m2,m3,stcum1d,i,j,k,kr,chem1_g(ind_chem(ispc),ngrid)%sc_t,stcum1d(ispc,k))
    enddo
    
   ! do ispc = naer_a,naer_z
   !         
   !   aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_t(kij) = &
   !   aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_t(kij) + stcum1d(ispc,k)
   !   
   ! 
   ! enddo
    
 enddo

IF(AEROSOL >= 1) THEN
 do k=1,ktop
    kr= k + 1   !nivel k da grade do grell corresponde ao nivel k + 1 do rams
    
    kij = kr + ij
    !-- aer mass concentration    
    do ispc = naer_a,naer_z
            
      aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_t(kij) = &
      aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_t(kij) + stcum1d(ispc,k)
    enddo
  enddo
   
    IF(AEROSOL == 2) THEN
     do k=1,ktop
       kr= k + 1  
    
       kij = kr + ij
       !-- aer number concentration    
       do ispc = naer_a_number,naer_z_number
            
          aer2_g(ind_mode_number(ispc),ngrid)%sc_t(kij) = &
          aer2_g(ind_mode_number(ispc),ngrid)%sc_t(kij) + stcum1d(ispc,k)
       enddo
     enddo
     !-- aer mass inorg 
     do k=1,ktop
       kr= k + 1   
       kij = kr + ij
       !-- aer mass concentration    inorg
       do ispc = naer_a_inorg,naer_z_inorg
           
          aer1_inorg_g(ind_aer_inorg(ispc),ngrid)%sc_t(kij) = &
          aer1_inorg_g(ind_aer_inorg(ispc),ngrid)%sc_t(kij) + stcum1d(ispc,k)
       enddo
     enddo

    ENDIF   
ENDIF
!  !Impose the Mass conservation
!  rmass=0.
!  rmass_p=0.
!  rmass_n=0.  
!  do k=1,n1
!     dp = 100.*(p_cup(k)-p_cup(k+1))
!     rmass   = rmass   + stcum1d(k)*dp
!     if(stcum1d(k) .gt. 0.) rmass_p = rmass_p + stcum1d(k)*dp
!     if(stcum1d(k) .lt. 0.) rmass_n = rmass_n + stcum1d(k)*dp
!  enddo
!  fr = abs(rmass_p/rmass_n)  
!  rmass=0.
!  do k=1,n1
!     dp = 100.*(p_cup(k)-p_cup(k+1))
!     if(stcum1d(k) .lt. 0) stcum1d(k) = fr  *  stcum1d(k)
!     rmass = rmass + stcum1d(k)*dp
!  enddo

end subroutine apply_conv_tend_GF

!----------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!-srf: convective transport for GD scheme 
!
subroutine trans_conv_mflx(iens,stcum)
!------------------------------------------------------------------
!-Convective transport for non-hygroscopic gases and aerosol 
!-Developed by Saulo Freitas (sfreitas@cptec.inpe.br)
!-ref: Freitas, et.al,: Monitoring the transport of biomass burning 
!      emissions in South America. Environmental Fluid Mechanics, 
!      Kluwer Academic Publishers, 2005.
!------------------------------------------------------------------
  use mem_tconv         ,only: stcum1d,dn01d,se,se_cup,sc_up,sc_dn,sc_up_c,sc_dn_c, &
                               henry_coef,pw_up,pw_dn,&
                               trans_conv_alloc, &    !Control alloc variable
                               alloc_trans_conv, &    !alloc subroutine
                               zero_tconv
  !use chem1_list, only: CO
  use mem_chem1         ,only: chem1_g, CHEMISTRY,nspecies=>NSPECIES_TRANSPORTED
  use node_mod          ,only: m1=>mzp,m2=>mxp,m3=>myp,ia,iz,ja,jz,i0,j0,mynum
  use mem_grid          ,only: dtlt,ngrid,naddsc
  use mem_scratch       ,only: scratch
  use mem_basic         ,only: basic_g
  use mem_cuparm        ,only: cuparm_g
  use mem_grell_param   ,only: mgmxp,mgmyp,mgmzp,maxiens,ngrids_cp
  use mem_scratch1_grell, only: ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d,kpbl4d,   &
                                kstabi4d,kstabm4d,xmb4d,edt4d,enup5d,endn5d,deup5d, &
                                dedn5d,zup5d,zdn5d,iruncon,zcup5d,pcup5d,prup5d,    &
				prdn5d,pwav4d,clwup5d,tup5d
  implicit none  

  integer, intent(IN) :: iens
  integer :: i,j,k,kr,ipr,jpr,iscl,iconv,iwet
  real,intent(INOUT) :: stcum(m1,m2,m3)


  !          se(:,:)  ! environment scalar profile z levels
  !      se_cup(:,:)  ! environment scalar profile z_cup levels
  !       sc_up(:,:)  ! updraft	gas-phase  scalar profile
  !       sc_dn(:,:)  ! downdraft gas-phase  scalar profile
  !     stcum1d(:,:)  ! 1d convective tendency
  !     sc_up_c(:,:)  ! updraft	aqueous-phase scalar profile
  !     sc_dn_c(:,:)  ! downdraft aqueous-phase scalar profile
  !       pw_up(:,:)  ! updraft precitable gas/aer
  !       pw_dn(:,:)  ! downdraft precitable gas/aer
  !  henry_coef(:,:)  ! henry's constant for gases wet removal
  !       dn01d(:)    ! 1d air density

  real, dimension(2) :: c0
  data (c0(i),i=1,2)  /0.002 &  ! deep    convection (iens=1), unit of c0 is m^-1
                      ,0.000 /  ! shallow convection (iens=2)
		      
  if(CHEMISTRY < 0) return 

  if(.not. trans_conv_alloc) then
     call alloc_trans_conv(nspecies,mgmzp) 
     call zero_tconv()
  end if

  if(nspecies == 0) return ! if there is not any specie to be trasnported
                           ! => return
  
  iwet = 1
  
  !---- Coordenadas para escrita ascii
  ipr=0 - i0
  jpr=0 - j0
  
  do j=ja,jz
        do i=ia,iz
	   !call azero(nspecies*mgmzp,stcum1d) ! ver necessidade <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
           
	   iconv = 0       !flag para definir se ha conveccao ou nao
           !srf- fev-2003: por enquanto nao usando o loop no espectro de nuvens
           !               mude posteriormente quando o codigo de shallow e deep
           !               estiverem fundidos  em um so.
           ! DO iens=1,maxiens            !loop no espectro de nuvens

           !     verifica se ha conveccao em i,j,iens
           if(ierr4d(i,j,iens,ngrid) .eq. 0) then       
	     
              stcum1d = 0.
!--(DMK-BRAMS-5.0-OLD)-------------------------------------------  
!	      call azero(nspecies*mgmzp,stcum1d)
!--(DMK-BRAMS-5.0-FIM)-------------------------------------------  


              iconv = 1      
              !srf -fev-2003: use a linha abaixo somente qdo o espectro de nuvens estiver funcionando
              !        if(iens.eq.1) call get_dn01d(mgmzp,m1,dn01d,dn0(1,i,j))
              
	      call get_dn01d(mgmzp,m1,dn01d,basic_g(ngrid)%dn0(1,i,j))
              
	      call get_se(nspecies,mgmzp,ngrid,m1,m2,m3,i,j,se,se_cup)

              if(iwet == 1) then

	       call get_sc_up_wet(nspecies,  mgmzp,m1,se,se_cup,sc_up,  &
				    k224d  (i,j,iens,ngrid),  &
        			    kbcon4d(i,j,iens,ngrid),  &
				    ktop4d (i,j,iens,ngrid),  &
        			  deup5d (1,i,j,iens,ngrid),  &
				  enup5d (1,i,j,iens,ngrid),  &
        			  zcup5d (1,i,j,iens,ngrid),  &
				  sc_up_c,henry_coef,pw_up ,  &
				  dn01d,c0(iens)           ,  &        
				  clwup5d(1,i,j,iens,ngrid),  &
        			  tup5d  (1,i,j,iens,ngrid),  &
				  zup5d  (1,i,j,iens,ngrid)   )
	      else 
	       call get_sc_up    (nspecies, mgmzp,m1,se,se_cup,sc_up,   &
		                       k224d(i,j,iens,ngrid),  &
                                     kbcon4d(i,j,iens,ngrid),  &
				    deup5d(1,i,j,iens,ngrid),  &
                                    enup5d(1,i,j,iens,ngrid),  &
				    zcup5d(1,i,j,iens,ngrid)   )
	      endif

              call get_sc_dn(nspecies,mgmzp,m1,se,se_cup,sc_dn,jmin4d(i,j,iens,ngrid),   &
                               kdet4d(i,j,iens,ngrid),dedn5d(1,i,j,iens,ngrid), &
                             endn5d(1,i,j,iens,ngrid),zcup5d(1,i,j,iens,ngrid)  )

              call get_stcum(nspecies,mgmzp,m1,dn01d,stcum1d,se,se_cup,sc_up,sc_dn,      &
                   	        xmb4d(i,j,iens,ngrid),   edt4d(i,j,iens,ngrid), &
                   	       jmin4d(i,j,iens,ngrid),  kdet4d(i,j,iens,ngrid), &
                   	        k224d(i,j,iens,ngrid), kbcon4d(i,j,iens,ngrid), &
                   	       ktop4d(i,j,iens,ngrid),  kpbl4d(i,j,iens,ngrid), &
                   	     kstabi4d(i,j,iens,ngrid),kstabm4d(i,j,iens,ngrid), &
                   	     zcup5d(1,i,j,iens,ngrid),pcup5d(1,i,j,iens,ngrid), &
                   	     deup5d(1,i,j,iens,ngrid),enup5d(1,i,j,iens,ngrid), &
                   	     dedn5d(1,i,j,iens,ngrid),endn5d(1,i,j,iens,ngrid), &
                   	      zup5d(1,i,j,iens,ngrid), zdn5d(1,i,j,iens,ngrid)  )

              if(iwet == 1) then    
             !---------------------- WET REMOVAL SCHEMES ------------------------------
                 !-- contribuicao devido `a concentracao aquosa desentranhada
                 !-- junto com a agua liquida
		 call get_stcum_detrain(nspecies,mgmzp,m1,dn01d,stcum1d,sc_up_c,sc_dn_c,    &
				   xmb4d(i,j,iens,ngrid),   edt4d(i,j,iens,ngrid), &
				  jmin4d(i,j,iens,ngrid),  kdet4d(i,j,iens,ngrid), &
				   k224d(i,j,iens,ngrid), kbcon4d(i,j,iens,ngrid), &
				  ktop4d(i,j,iens,ngrid),  kpbl4d(i,j,iens,ngrid), &
				kstabi4d(i,j,iens,ngrid),kstabm4d(i,j,iens,ngrid), &
				zcup5d(1,i,j,iens,ngrid),pcup5d(1,i,j,iens,ngrid), &
				deup5d(1,i,j,iens,ngrid),enup5d(1,i,j,iens,ngrid), &
				dedn5d(1,i,j,iens,ngrid),endn5d(1,i,j,iens,ngrid), &
				 zup5d(1,i,j,iens,ngrid), zdn5d(1,i,j,iens,ngrid)  )

		 !-- calcula a massa removida e depositada na superficie
		 call get_wet_deposited_mass(nspecies,mgmzp,m1,dtlt,dn01d,pw_up,pw_dn, &
					 xmb4d(i,j,iens,ngrid), edt4d(i,j,iens,ngrid), &
				       kbcon4d(i,j,iens,ngrid),ktop4d(i,j,iens,ngrid), &
				       i,j,ngrid  )


              endif

           call apply_conv_tend(nspecies,mgmzp,m1,m2,m3,stcum1d,i,j,ngrid,ktop4d(i,j,iens,ngrid))



 
           endif    ! endif of test if is there or not convection
           !      enddo     ! fim do loop no espectro de nuvens

       !    !Change to the 3d array tendency 
       !    if(iconv .eq. 1) then
       !       do k=1,m1-1
       !          kr= k + 1   !nivel k da grade do grell corresponde ao nivel k + 1 do rams
       !!srf          stcum(kr,i,j)=stcum1d(k)
       !       enddo
       !    endif  !endif do iconv      
    enddo     !loop em i
  enddo        !loop em j

     !     acumula a tendencia assoc. ao trans convectivo `a tendencia total
!!      call accum(m1*m2*m3,chem1_g(iscl,ngrid)%sc_t(1),stcum)

end subroutine trans_conv_mflx
 
!--------------------------------------------------

subroutine get_se(nspecies,mgmzp,ngrid,n1,n2,n3,i,j,se,se_cup) 
!  use chem1_list, only: spc_alloc_chem  =>spc_alloc,transport,on,off
  use mem_chem1,  only: chem1_g
  use mem_aer1,  only: aer1_g
  use aer1_list, only: nmodes
 
  use mem_tconv, only : nchem_a,nchem_z, ind_chem &
                       ,naer_a,naer_z,ind_aer, ind_mode
  
  implicit none
  integer :: nspecies,ngrid,mgmzp, n1, n2, n3, i, j
  real    :: se(nspecies,mgmzp), se_cup(nspecies,mgmzp)
  !local Variables
  integer :: k, kr,ispc

  do k=2,n1-1
     kr = k+1   ! nivel K da grade DO Grell corresponde ao nivel K + 1 DO RAMS

     ! chemistry section
     do ispc = nchem_a,nchem_z
       se    (ispc,k) =      chem1_g(ind_chem(ispc),ngrid)%sc_p(kr,i,j)         ! mixing ratio on Z levels
       
       se_cup(ispc,k) = .5*( chem1_g(ind_chem(ispc),ngrid)%sc_p(kr-1,i,j) + &
                             chem1_g(ind_chem(ispc),ngrid)%sc_p(kr  ,i,j)   )   ! mixing ratio on Z_CUP levels
     
       !if(se_cup(ispc,k) > 0)print*,'st=',ispc,k,se_cup(ispc,k),chem1_g(ind_chem(ispc),ngrid)%sc_p(kr,i,j)

!     do ispc=1,nspecies    
!       if(spc_alloc_chem(transport,ispc) == off) cycle
!       se    (ispc,k) =      chem1_g(ispc,ngrid)%sc_p(kr,i,j)  ! mixing ratio on Z levels
!       se_cup(ispc,k) = .5*( chem1_g(ispc,ngrid)%sc_p(kr-1,i,j) + chem1_g(ispc,ngrid)%sc_p(kr,i,j) )   ! mixing ratio on Z_CUP levels
!     enddo
  
     enddo
     
     !aerosol section
     
     do ispc = naer_a,naer_z
            
       se    (ispc,k) =      aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_p(kr,i,j)         ! mixing ratio on Z levels

       se_cup(ispc,k) = .5*( aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_p(kr-1,i,j) + &
                             aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_p(kr  ,i,j)   )   ! mixing ratio on Z_CUP levels
     enddo
     
  enddo

  do ispc = nchem_a,nchem_z
     
      se    (ispc,1) = chem1_g(ind_chem(ispc),ngrid)%sc_p(2,i,j)
      se_cup(ispc,1) = chem1_g(ind_chem(ispc),ngrid)%sc_p(2,i,j)

!  do ispc=1,nspecies 
!    if(spc_alloc_chem(transport,ispc) == off) cycle   
!    se    (ispc,1) = chem1_g(ispc,ngrid)%sc_p(2,i,j)
!    se_cup(ispc,1) = chem1_g(ispc,ngrid)%sc_p(2,i,j)
  enddo


  do ispc = naer_a,naer_z
     
      se    (ispc,1) = aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_p(2,i,j)
      se_cup(ispc,1) = aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_p(2,i,j)

  enddo



end subroutine get_se

!----------------------------------------------------------------------

subroutine get_sc_up(nspecies,mgmzp,n1,se,se_cup,sc_up,k22,kbcon,cd,entr,z_cup)

  implicit none
  integer :: nspecies,mgmzp,n1,k22,kbcon
  real, dimension(nspecies,mgmzp) ::se,se_cup,sc_up
  real, dimension(mgmzp) :: cd,entr,z_cup
  real ::dz
  integer :: k,ispc

  ! Scalar concentration in-cloud - updraft
  
    do k=1,k22-1
       do ispc=1,nspecies
              sc_up(ispc,k) = se_cup(ispc,k)
       enddo
    enddo

    do k=k22,kbcon
       do ispc=1,nspecies
              sc_up(ispc,k) = se_cup(ispc,k22)
       enddo
    enddo

    do k=kbcon+1,n1-1
 
       dz=z_cup(k)-z_cup(k-1)
       do ispc=1,nspecies
           sc_up(ispc,k)= ((1.-.5*cd(k)*dz)*sc_up(ispc,k-1)+entr(k)*dz*se(ispc,k-1)) / &
   		          (1. + entr(k)*dz - .5*cd(k)*dz )
          !expression II based on dH_c/dz = entr*(H_env - H_c)
          !	   sc_up(k)= ( (1.- .5*entr(k)*dz)*sc_up(k-1) + entr(k)*dz*se(k-1) ) &
          !		   / ( 1. + .5*entr(k)*dz )
       enddo
    enddo

end subroutine get_sc_up

!-------------------------------------------------------------------------

subroutine get_sc_dn(nspecies, mgmzp,n1,se,se_cup,sc_dn,jmin,kdet,cdd,entrd,z_cup )

  implicit none

  ! Arguments
  integer                :: nspecies,mgmzp, n1, jmin, kdet
  real, dimension(nspecies,mgmzp) ::  se, se_cup, sc_dn
  real, dimension(mgmzp) ::           cdd, entrd, z_cup

  ! Local variables
  integer                :: k,ispc
  real                   :: dz


  !Scalar concentration in-cloud - DOwndraft
  do k=jmin+1,n1-1
     do ispc=1,nspecies
       sc_dn(ispc,k) = se_cup(ispc,k)
     enddo
  enddo

  !k=jmim
  do ispc=1,nspecies
       sc_dn(ispc,jmin) = se_cup(ispc,jmin)
  enddo

  !added fev-2003 for shallow scheme
  if(jmin == 1 .or. jmin == 0 ) return

  do k=jmin-1,1,-1
     dz=z_cup(k+1)-z_cup(k)
     
     do ispc=1,nspecies
     
         sc_dn(ispc,k) = ((1.-.5*cdd(k)*dz)*sc_dn(ispc,K+1) + entrd(k)*dz*se(ispc,k)) /&
                          (1. + entrd(k)*dz - .5*cdd(k)*dz )
     enddo
  
  enddo

end subroutine get_sc_dn

!------------------------------------------------------------------------

subroutine get_stcum(nspecies,mgmzp,n1,dn0,stcum1d,se,se_cup,sc_up,sc_dn,xmb,edt,&
     	 	     jmin,kdet,k22,kbcon,ktop,kpbl,kstabi,kstabm,z_cup, &
     	 	     p_cup,cd,entr,cdd,entrd,zu,zd)

  use Phys_const, only: g

  implicit none

  ! Arguments
  integer                :: nspecies,mgmzp, n1, jmin, kdet, k22, kbcon, ktop, kpbl, kstabi, kstabm
  real                   :: xmb,edt
  real, dimension(nspecies,mgmzp) :: se,se_cup,sc_up,sc_dn,stcum1d
  real, dimension(mgmzp) :: z_cup,p_cup,cd,entr,cdd,entrd,zu,zd,dn0

  ! Local variables
  integer                :: k,ispc
  real                   :: dz, dp, entup, detup, entdoj, entupk, detupk, detdo, entdo,   &
                            subin, subdown, detdo1, detdo2

  !! edt = 0.  ! para nao considerar transporte por downdrafts

  do k=2,ktop
     dz =   z_cup(k+1) - z_cup(k)
     !definicao original de dp
     ! dp =  100.*( p_cup(k)   - p_cup(k+1) )
     !definicao em funcao de dn0 (densidade do ar do estado basico)
     ! o sinal - e' contabilizado na expressao para o ds/dt
     ! dp =   g*dn0(k+1)*dz ! nivel k da grade do grell corresponde ao nivel k + 1 do rams

     dp      = g*dn0(k)*dz  !dn0 aqui ja' esta' na grade DO GRELL.
     entup   = 0.
     detup   = 0.
     entdoj  = 0.
     entupk  = 0.
     detupk  = 0.
     detdo   = edt*  cdd(k)*dz*zd(k+1)
     entdo   = edt*entrd(k)*dz*zd(k+1)
     subin   = zu(k+1) - edt*zd(k+1)
     subdown = zu(k  ) - edt*zd(k  )

     if(k.ge.kbcon .and. k.lt.ktop)   then
        detup =   cd(k+1) *dz*zu(k)
        entup = entr(k)   *dz*zu(k)
     endif

     if(k.eq.jmin)          entdoj = edt*zd(k)
     if(k.eq.k22-1)         entupk = zu(kpbl)
     !if(k.eq.kpbl)         entupk = zu(kpbl)	 
     if(k.gt.kdet)          detdo  = 0.
     if(k.lt.kbcon)         detup  = 0.
     if(k.eq.ktop) then
        detupk  = zu(ktop)
        subin   = 0.
     endif

     do ispc=1,nspecies
          !tendency due cumulus transport ( k>=2)
          stcum1d(ispc,k) = stcum1d(ispc,k) +			    & !stcum1d comparece
        	       xmb*(				    & !no lado direito
        	       subin    *se_cup(ispc,k+1)		    & !para computar
        	       - subdown*se_cup(ispc,k  )		    & !a soma sobre o
        	       + detup*( sc_up(ispc,k+1)+ sc_up(ispc,k) )*.5  & !espectro de nuvens
        	       + detdo*( sc_dn(ispc,k+1)+ sc_dn(ispc,k) )*.5  & !(maxiens)
        	       - entup * se(ispc,k)			    &
        	       - entdo * se(ispc,k)			    &
        	       - entupk* se_cup(ispc,k22)		    &
        	       - entdoj* se_cup(ispc,jmin)		    &
        	       + detupk* sc_up(ispc,ktop)		    &
        	       )*g/dp

         !if(stcum1d(ispc,k) > 0)print*,'st=',ispc,k,stcum1d(ispc,k)
      end do

  end do
  !
  !tendency due cumulus transport (at bottom, k = 1)
  dz        =       z_cup(2)-z_cup(1)
  !definicao original de dp
  !      dp        = 100.*(p_cup(1)-p_cup(2))
  !definicao em funcao de dn0 (densidade do ar do estado basico)
  ! o sinal - e' contabilizado na expressao para o ds/dt
  ! dp =   g*dn0(2)*dz ! nivel k da grade do grell corresponde ao nivel k + 1 do rams
  dp =   g*dn0(1)*dz ! ja' esta' na grade do rams

  detdo1    = edt*zd(2)*cdd(1)*dz
  detdo2    = edt*zd(1)
  entdo     = edt*zd(2)*entrd(1)*dz
  subin     =-edt*zd(2)

  do ispc=1,nspecies
       stcum1d(ispc,1) = stcum1d(ispc,1) +		      &
   		    xmb*(			      &
   		    detdo1*(sc_dn(ispc,1)+sc_dn(ispc,2))*.5     &
   		    + detdo2*	 sc_dn(ispc,1)	      &
   		    + subin *	 se_cup(ispc,2)	      &
   		    - entdo *	 se(ispc,1)  	      &
   		    )*g/dp

  enddo
  return
end subroutine get_stcum

!---------------------------------------------------------------------------

subroutine get_sc_up_wet(nspecies,  mgmzp,n1,se,se_cup,sc_up,k22,kbcon,ktop,cd, &
     entr,z_cup,sc_up_c,henry_coef,pw_up,dn0,c0, &
     cupclw,tup,zu)
  !DSM use chem1_list, only : CO!, H2O2
  use mem_tconv, only : nchem_a,nchem_z, ind_chem &
                       ,naer_a,naer_z,ind_aer, ind_mode

  implicit none
  integer nspecies, mgmzp,n1,k,k22,kbcon,ktop, ispc
  
 ! integer, intent(in)           :: spc_number
 ! character (len=*), intent(in) :: spc_name
  real dz,c0
  real, dimension(nspecies) ::  conc_equi,conc_mxr,qrch
  real, dimension(mgmzp) :: cd,entr,z_cup, &
         tup ,       & ! local temperature in cloud updraft [K]
       cupclw,       & ! up cloud liquid water [kg(water)/kg(air)] from cum. scheme
           zu,       & ! norm mass flux updrat
          dn0         ! basic state air density [kg[air]/m^3]

  real, dimension(nspecies,mgmzp) :: &
       se,         &
       se_cup ,    & ! gas-phase in environment mixing ratio [kg(gas phase)/kg(air)]
       sc_up ,     & ! gas-phase in updraft     mixing ratio [kg(gas phase)/kg(air)]
       henry_coef, & ! Henry's constant [(kg(aq)/kg(water))/(kg(gas phase)/kg(air))]
       sc_up_c,    & ! aqueous-phase in updraft mixing ratio [kg(aq)/kg(air)]
       pw_up         ! precitable gas/aer amount [kg[..]/kg[air]

  real, parameter :: scav_eff = 0.6  ! for smoke : Chuang et al. (1992) J. Atmos. Sci.


  !--- Set zero pw_up, henry_coef and sc_up_c arrays
!--(DMK-BRAMS-5.0-INI)-------------------------------------------  
  sc_up_c    = 0.
  pw_up      = 0.
  henry_coef = 0.
!--(DMK-BRAMS-5.0-OLD)-------------------------------------------  
!  call azero3(nspecies*mgmzp,sc_up_c,pw_up,henry_coef)
!--(DMK-BRAMS-5.0-FIM)-------------------------------------------  


  !--- Get Henry's law constants - unit [(kg(aq)/kg(water))/(kg(gas phase)/kg(air))]
  CALL henry(nspecies,mgmzp,n1,kbcon,ktop,henry_coef,tup) 


  !--- Scalar concentration in-cloud updraft - gas phase 
  do k=1,k22-1
       do ispc=1,nspecies
              sc_up(ispc,k) = se_cup(ispc,k)
       enddo
  enddo

  do k=k22,kbcon
       do ispc=1,nspecies
              sc_up(ispc,k) = se_cup(ispc,k22)
       enddo
  enddo

  do k=kbcon+1,ktop
       dz=z_cup(k)-z_cup(k-1)
       !
       !---- steady state plume equation, for what could
       !---- be in cloud before anything happens (kg/kg*1.e9 = ppbm)
       !
       do ispc=1,nspecies
         sc_up(ispc,k)= ((1.-.5*cd(k)*dz)*sc_up(ispc,k-1)+entr(k)*dz*se(ispc,k-1))/ &
                        ( 1. + entr(k)*dz - .5*cd(k)*dz )
   
       !expression II based on dH_c/dz = entr*(H_env - H_c)
       !           sc_up(k)= ( (1.- .5*entr(k)*dz)*sc_up(k-1) + entr(k)*dz*se(k-1) ) &
       !                   / ( 1. + .5*entr(k)*dz )
       !
       enddo
       
       
       !- chemistry section
       !
       do ispc = nchem_a,nchem_z
            !--- equilibrium tracer concentration - Henry's law
	    !--- cloud liquid water tracer concentration
            conc_mxr(ispc) = (henry_coef(ispc,k)*cupclw(k) /(1.+henry_coef(ispc,k)*cupclw(k)) )* sc_up(ispc,k) 

            !---   gas-phase concentration in updraft (not incorporated into cloud liquid water)
            qrch(ispc) = sc_up(ispc,k) - conc_mxr(ispc)  
	    
            !if(ind_chem(ispc) == CO .or. ind_chem(ispc)== H2O2) then
        	!print*,'xx1=',   k,cupclw(k)*1000 , sc_up(ispc,k)
	        ! print*,'xx2=',  (henry_coef(ispc,k)/(1+henry_coef(ispc,k))) , conc_mxr(ispc),  qrch(ispc) 
	        !print*,'-----------------------------------------------'
	    ! endif
       enddo

       !- aerosol section
       
       do ispc = naer_a,naer_z
     
            conc_mxr(ispc)  =  scav_eff* sc_up(ispc,k) !unit [kg(aq)/kg(air)]  for aerosol/smoke 
            !
            !---   gas-phase concentration in updraft (not incorporated into cloud liquid water)
            qrch(ispc) = sc_up(ispc,k) - conc_mxr(ispc)                 
        
            !xxx               conc_mxr = sc_up(k)   ! for sulphur ???
       
       enddo
       
       
       !---   aqueous-phase concentration in updraft after rainout
       !---  'sc_up_c' here would be the part that is carried in cloud water
       do ispc = 1,nspecies
        !sc_up_c(k) =(sc_up(k)-qrch)/(1.+c0*dz      ) ! bug ??
        sc_up_c(ispc,k) = (sc_up(ispc,k)-qrch(ispc))& ! conc. dentro da gota, ja computando 
	                   /(1.+c0*dz*zu(k))           ! a parte removida pela precipitacao
        !
        !---  pw_up is the part that wil fall out in rain
        !
        pw_up(ispc,k) = c0*dz*sc_up_c(ispc,k)*zu(k)      
        !
        if(sc_up_c(ispc,k).lt.0.) sc_up_c(ispc,k) = 0.
       enddo
       
       
       do ispc = 1,nspecies
            !----- set next level
            sc_up(ispc,k) = sc_up_c(ispc,k) + qrch(ispc)  ! conc total = conc na agua liquida   (= sc_up_c ) 
                                                          ! + conc na corrente de ar (= qrch    )
       enddo

  enddo
  !-----

  do k=ktop+1,n1-1
       do ispc = 1,nspecies
          sc_up(ispc,k) = se_cup(ispc,k)
       enddo
  enddo
  !
  !----- get back the in-cloud updraft gas-phase mixing ratio : sc_up(ispc,k)
  !      
  do k=kbcon+1,ktop
       do ispc = 1,nspecies
          sc_up(ispc,k) = sc_up(ispc,k) - sc_up_c(ispc,k)
       enddo
  enddo

end subroutine get_sc_up_wet

!---------------------------------------------------------------------

subroutine henry(nspecies,mgmzp,n1,kbcon,ktop,henry_coef,temp)
  use mem_tconv, only : nchem_a,nchem_z, ind_chem
! change MP 11/12/07 
  use chem1_list, only : hstar,weight,dhr, ak0, dak! , CO, HNO3, ORA1 !, H2O2
! end change MP
implicit none
integer nspecies,mgmzp,n1,kbcon,ktop,k,ispc
real, dimension(mgmzp) :: temp!,dn0
real, dimension(nspecies,mgmzp) :: henry_coef
real pi ,rgas ,avogad ,rhoh2o, temp0 ,temp0i ,conv3,conv4 ,conv5 ,conv6,conv7 ,fct ,tcorr, &
     hplus,corrh
real XXXw,XXXh  !molecular weight
! define some constants!
parameter( pi    =	3.14159265)! pi number
parameter( rgas  =	  8.32e-2 )!atm M^-1 K^-1 !   8.314)! gas constant [J/(mol*K)]
parameter( avogad=        6.022e23)! Avogadro constant [1/mol]
parameter( rhoH2O=        999.9668)! density of water [kg/m3]
parameter( temp0 =          298.15)! standard temperature [K]
parameter( temp0i =      1./298.15)! inverse of standard temperature [K]

! miscellaneous
! PRES	   pressure [Pa]
! TCORR    temperature correction term [dimensionless]
! TEMP	   temperature [K]
! LWC      liquid water content [m3(water)/m3(air)]

!constants of convertion
parameter( conv3     = avogad / 1.0e6)!  [mol(g)/m3(air)]  to [molec(g)/cm3(air)]
parameter( conv4     = 100.          )!  [m]		to [cm]
parameter( conv5     = 1000.         )!  [m^3]  	to [l]
parameter( conv7     = 1/conv5       )!  [l]  	to [m^3]
parameter( conv6     = 1. / 101325.  )!  [Pa]              to [atm]


! calculate the molecular weight
!! call calc_mol_weight(XXXw)  !units [kg/mol]


!loop over the atmosph. column which may have liq. water content
do k=kbcon+1,ktop


! aqueous-phase concentrations XXXa [mol/m3(air)]!
! gas-phase concentrations XXXg [mol/m3(air)]!
! Henry constants XXXh for scavenging [mol/(l*atm)]!
! converted to [(mol(aq)/m3(aq))/(mol(g)/m3(air))], i.e. dimensionless!
! in equilibrium XXXa = XXXh * LWC * XXXg!
 tcorr = 1./temp(k) - temp0i
 !fct   = conv5 * rgas * temp(k) * conv6
 fct   = conv7 * rgas * temp(k)
 
! change MP 11/12/07 taking into account the acid dissociation constant
! ak=ak0*exp(dak*(1/t-1/298))
! for cloud water. pH is asuumed to be 3.93 as in RAMS-chimie
!   pH=3.93
!   hplus=10**(-pH)
    hplus=1.175E-4
    do ispc = nchem_a,nchem_z
       corrh=1+ak0(ind_chem(ispc))*exp(dak(ind_chem(ispc))*tcorr)/hplus     
   
       henry_coef(ispc,k) =  hstar(ind_chem(ispc))* exp( dhr(ind_chem(ispc))*tcorr) * fct &
                                      * corrh !* dn0(k) / rhoH2O
!    if(ind_chem(ispc) == HNO3 .or. ind_chem(ispc) == ORA1) then
!       print*,k,ind_chem(ispc), hstar(ind_chem(ispc)),tcorr,fct, henry_coef(ind_chem(ispc),k),corrh
!     endif
!
! end change MP 
!
      ! print*,'yy', tcorr ,XXXh ,hstar(ind_chem(ispc)),henry_coef(ind_chem(ispc),k)
      !if(ind_chem(ispc) == CO .or. ind_chem(ispc) == H2O2) then
      !endif
! select the appropiate gas 
! O3h	  = 1.2e-2	* exp( 2560.*tcorr) * fct!			    
! O2h	  = 1.3e-3	* exp( 1500.*tcorr) * fct!
! XXXh = O2h
!
! henry_coef(k) = XXXh      ! unit [(mol(aq)/m3(aq))/(mol(gas phase)/m3(air))]
!
! henry_coef(k) = XXXh * dn0(k) / rhoH2O   ! unit [(mol(aq)/kg(water))/(mol(gas phase)/kg(air))]
!observe that multiplying and dividing by XXXw [kg/mol] => the units are [(kg(aq)/kg(water))/(kg(gas phase)/kg(air))]
!so the units of the last calculation of henry_coef(k) is also [( kg(aq)/kg(water))/( kg(gas phase)/kg(air))]
! print*,k,temp(k),fct,XXXh,henry_coef(k)
 
   enddo
enddo

end subroutine henry

!---------------------------------------------------------------------

subroutine get_stcum_detrain(nspecies,mgmzp,n1,dn0,stcum1d,sc_up_c,sc_dn_c,  &
     	  		     xmb,edt,jmin,kdet,k22,kbcon,ktop,kpbl, &
     	  		     kstabi,kstabm,z_cup,p_cup,cd,entr,cdd, &
     	  		     entrd,zu,zd)

  use Phys_const, only: g  			

  implicit none

  ! Arguments
  integer                :: nspecies,mgmzp,n1,jmin,kdet,k22,kbcon,ktop,&
                            kpbl,kstabi,kstabm,ispc
  real                   :: xmb,edt
  real, dimension(nspecies,mgmzp) :: sc_up_c,sc_dn_c,stcum1d
  real, dimension(mgmzp) :: z_cup,p_cup,cd,entr,cdd,entrd,zu,zd,dn0

  ! Local Variables
  integer :: k
  real    :: dz, dp,detup, detupk

  do k=kbcon+1,ktop
     dz =   z_cup(k+1) - z_cup(k)
     dp =   g*dn0(k)*dz  
     detup  = 0.
     if(k.lt.ktop)  detup  = cd(k+1) *dz*zu(k)
     detupk = 0.
     if(k.eq.ktop)  detupk = zu(ktop)
     
     do ispc=1,nspecies
     
      !- tendency due cumulus transport ( k>=2)
       stcum1d(ispc,k) = stcum1d(ispc,k) +                 & ! stcum1d do lado direito corresponde 
                  xmb*( 				   & !         `a contribuicao da fase 
                  + detup*( sc_up_c(ispc,k+1)+ sc_up_c(ispc,k) )*.5  & !  	gasosa.
                  + detupk* sc_up_c(ispc,ktop)		   & ! stcum1d final = fase gasosa +
                  )*g/dp				     !  	       fase liquida
     end do

  end do
end subroutine get_stcum_detrain

!---------------------------------------------------------------------

subroutine get_wet_deposited_mass(nspecies,mgmzp,m1,dt,dn01d,pw_up,pw_dn, &
                                  xmb,edt,kbcon,ktop,i,j,ngrid)
use mem_chem1 ,ONLY: chem1_g
use mem_aer1  ,ONLY: aer1_g
use mem_tconv, only : nchem_a,nchem_z, ind_chem &
                     ,naer_a,naer_z,ind_aer, ind_mode
use chem1_list, only : chem_spc_alloc=>spc_alloc, off, wdp
use aer1_list , only : aer_spc_alloc =>spc_alloc

  implicit none
  integer k,nspecies,mgmzp,m1,kbcon,ktop,i,j,ngrid,ispc
  real dt,wetdep,xmb,edt
  real, dimension(mgmzp) :: dn01d
  real, dimension(nspecies,mgmzp) :: pw_up,pw_dn

  !--- wetdep units kg m^-2
  !!wetdep = 0.  ! use for instantaneous rate (not integrated over the time)
  do  k=1,ktop+1

    !- chemistry section
    do ispc = nchem_a,nchem_z

       if(chem_spc_alloc(wdp,ind_chem(ispc)) == off) cycle
       
       chem1_g(ind_chem(ispc),ngrid)%sc_wd(i,j) = chem1_g(ind_chem(ispc),ngrid)%sc_wd(i,j) & 
	                                  + dt*(pw_up(ispc,k)+edt*pw_dn(ispc,k))*xmb       
    enddo
    
    !- aerosol section
    do ispc = naer_a,naer_z
      
      if(aer_spc_alloc(wdp,ind_mode(ispc),ind_aer(ispc)) == off) cycle

       aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_wd(i,j) = &
       aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_wd(i,j) + dt*(pw_up(ispc,k)+edt*pw_dn(ispc,k))*xmb 

      enddo

  !   print *,'pw_up,pw_dn,WM:',k,pw_up(ispc,k),pw_dn(ispc,k),wetdep
  enddo

end subroutine get_wet_deposited_mass

!---------------------------------------------------------------------

subroutine get_dn01d(mgmzp,m1,dn01d,dn0)

  implicit none
  integer mgmzp,m1,k
  real, dimension(mgmzp) :: dn01d
  real, dimension(m1)    :: dn0

  do k=1,m1-1
     dn01d(k) = dn0(k+1)  ! nivel k da grade do grell corresponde ao nivel k + 1 do rams
  enddo
  dn01d(m1)=dn01d(m1-1)

end subroutine get_dn01d
!-------------------------------------------------------------------------------

subroutine apply_conv_tend(nspecies,mgmzp,m1,m2,m3,stcum1d,i,j,ngrid,ktop)
use mem_chem1 ,ONLY: chem1_g
use mem_aer1  ,ONLY: aer1_g
use mem_tconv, only : nchem_a,nchem_z, ind_chem &
                     ,naer_a,naer_z,ind_aer, ind_mode
implicit none
integer, intent(IN)::nspecies,mgmzp,m1,m2,m3,i,j,ngrid,ktop
real, dimension(nspecies,mgmzp) :: stcum1d

!  local variables
integer ij,kij,ispc,k,kr
! the expression: kij= k +  m1*(i-1) + (m1*(m2-1)+m1)*(j-1)
! provides the memory position of (k,i,j) at tendency array (1-d)

 ij=  m1*(i-1) + (m1*(m2-1)+m1)*(j-1)   
 
! do k=1,m1-1
 do k=1,ktop
    kr= k + 1   !nivel k da grade do grell corresponde ao nivel k + 1 do rams
    
    kij = kr + ij
    
    ! chemistry section
    do ispc = nchem_a,nchem_z

       ! old way : stcum(kr,i,j)=stcum1d(k)
       chem1_g(ind_chem(ispc),ngrid)%sc_t(kij)  = chem1_g(ind_chem(ispc),ngrid)%sc_t(kij) +&
                                                  stcum1d(ispc,k)
       !call xxxxxx(mgmzp,m1,m2,m3,stcum1d,i,j,k,kr,chem1_g(ind_chem(ispc),ngrid)%sc_t,stcum1d(ispc,k))
    enddo
    
    do ispc = naer_a,naer_z
            
      aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_t(kij) = &
      aer1_g(ind_mode(ispc),ind_aer(ispc),ngrid)%sc_t(kij) + stcum1d(ispc,k)
   
    enddo
    
enddo

end subroutine apply_conv_tend
!--------------------------------------------------------------

subroutine xxxxxx(mgmzp,m1,m2,m3,stcum1d,i,j,k,kr,st,s1t)
real st(m1,m2,m3)
real s1t
!stop 333
st(kr,i,j)=st(kr,i,j)+s1t
if(abs(s1t) > 1.e-6) stop 3355
end subroutine xxxxxx
!--------------------------------------------------------------

subroutine mass_conserv(mgmzp,n1,dt,dn0,se,stcum1d,pw_up,pw_dn,wetdep, &
                        z_cup,p_cup,xmb,edt)

  use Phys_const, only: g  			

  implicit none

  integer mgmzp,n1,k
  real, dimension(mgmzp) :: se,stcum1d,z_cup,p_cup,dn0,pw_up,pw_dn
  real :: wetdep,dt,rmi,rmf,rmass,dp,dz,xmb,edt
  real :: rmass_p,rmass_n,fr

  !- Mass conservation:
  ! INTEGRAL ( stcum1d * Air_dens * dZ ) = 0
  ! Air_dens = -1/g * Dp/Dz
  ! ---
  rmi  =0.
  rmf  =0.
  rmass=0.
  do k=1,n1-1
     dz =   z_cup(k+1) - z_cup(k)
     dp =   -g*dn0(k)*dz  !dn0 aqui ja' esta' na grade DO GRELL.
     ! dp = -100.*(p_cup(k)-p_cup(k+1))
     ! PRINT*,-g*dn0(k)*dz,-100.*(p_cup(k)-p_cup(k+1))
     rmass = rmass + stcum1d(k)*(-dp/g) + pw_up(k)*xmb
     rmi   = rmi   +   se(k)                  *(-dp/g)
     rmf   = rmf   + ( se(k) + stcum1d(k)*dt )*(-dp/g)
  enddo
  rmf = rmf
  print*,'rmass,rmi,rmf,(rmf-rmi)/rmi*100.'
  print*,rmass,rmi,rmf,(rmf-rmi)/rmi*100.
  rmf = rmf + wetdep
  print*,'com wetdep=',wetdep
  print*,rmass,rmi,rmf,(rmf-rmi)/rmi*100.
  print*,'-------------------------------------------------------'
  return

  !Impose the Mass conservation
  rmass=0.
  rmass_p=0.
  rmass_n=0.  
  do k=1,n1
     dp = 100.*(p_cup(k)-p_cup(k+1))
     rmass   = rmass   + stcum1d(k)*dp
     if(stcum1d(k) .gt. 0.) rmass_p = rmass_p + stcum1d(k)*dp
     if(stcum1d(k) .lt. 0.) rmass_n = rmass_n + stcum1d(k)*dp
  enddo
  fr = abs(rmass_p/rmass_n)  
  rmass=0.
  do k=1,n1
     dp = 100.*(p_cup(k)-p_cup(k+1))
     if(stcum1d(k) .lt. 0) stcum1d(k) = fr  *  stcum1d(k)
     rmass = rmass + stcum1d(k)*dp
  enddo

end subroutine mass_conserv

!--------------------------------------------------------------

subroutine wet_removal(mgmzp,m1,dt,dn0,wdepmass,stcum1d,z_cup)

  implicit none

  real dt,wdepmass,wr,dz
  integer mgmzp,m1,k
  real, dimension(mgmzp) :: stcum1d,z_cup,dn0
  !REAL, DIMENSION(m1)    :: dn0
  data wr/0.9/    !factor based on Andreae et al. GRL 2001 (80 - 95%)

  do k=1,m1-1
     if(stcum1d(k) .gt. 0.) then
        dz         = z_cup(k+1) - z_cup(k)
        wdepmass   = wdepmass + dt * (stcum1d(k)*wr)*dn0(k)*dz
        stcum1d(k) = stcum1d(k)*(1.-wr)
     endif
  enddo

end subroutine wet_removal

!----------------------------------------------------------------------

subroutine wet_removal_Berge(mgmzp,m1,dt,dn0,wdepmass,precip,stcum1d, &
     sc_up,k22,kbcon,ktop,z_cup,p_lw)

  implicit none
  integer k,mgmzp,m1,k22,kbcon,ktop
  real scav_eff,coll_eff,A,v,Mw,dz,dt,wdepmass,tEND_wd,precip
  real, dimension(mgmzp) :: sc_up,stcum1d,z_cup,p_lw,dn0
  !REAL, DIMENSION(m1)    :: dn0
  data scav_eff/0.6/    ! for smoke : Chuang et al. (1992) J. Atmos. Sci.
  data coll_eff/0.01/   ! the mean collection efficiency averaged over all raindrop water
  data A/5.2/           ! constant, unit : m^3/(kg s)
  data v/5./            ! raindrop fall speed (m/s)

  !------------ In-cloud wet deposition ' washout'----
  ! Following Berge, TELLUS, 1993, eq. 3.7a
  do  k=kbcon,ktop

     !Units:
     ! sc_up (in-cloud aerosol mixing ratio)   : kg[aer]   / kg[air]
     ! p_lw  (precip rate/ liq water content)  :(kg[water] / (m^2 s)) / (kg[water]/kg[air])
     !                                           = kg[air]/( m^2 s)
     ! dno   (basic state air density)         : kg[air]   / m^3
     ! dz    (thickness of vertical layer)     : m
     ! tEND_wd (wet dep tENDency)              : kg[aer]/ (kg[air] s)

     dz         = z_cup(k+1) - z_cup(k)
     tEND_wd    = scav_eff*sc_up(k)*p_lw(k) /( dn0(k)*dz )

     !Accumulate the deposited mass on surface.
     wdepmass   = wdepmass + dt * tEND_wd*dn0(k)*dz


     !------
     print*,'WET1',k,kbcon,ktop,wdepmass
     print*,'WET1',tEND_wd,stcum1d(k),-tEND_wd+stcum1d(k)

     !Final tENDency = convective transport - wet deposition removal
     stcum1d(k) = stcum1d(k) - tEND_wd
  enddo

  !------------ Below-cloud wet deposition 'rainout'----
  ! Following Berge, TELLUS, 1993, eq. 3.8
  !Units:
  ! sc_up (below cloud aerosol mixing ratio): kg[aer]   / kg[air]
  ! precip  (total precip rate)             :(kg[water] / (m^2 s)) 
  ! Mw (water precipitation concentration   : kg[water]/m^3 
  ! tEND_wd (wet dep tENDency)              : kg[aer]/ (kg[air] s)
  Mw = precip/v

  do  k=1,kbcon-1
     dz         = z_cup(k+1) - z_cup(k)
     tend_wd    = a * mw * coll_eff * sc_up(k)
     !accumulate the deposited mass on surface.
     wdepmass   = wdepmass + dt * tEND_wd*dn0(k)*dz
     !Final tENDency = convective transport - wet deposition removal
     stcum1d(k) = stcum1d(k) - tEND_wd
  enddo

end subroutine wet_removal_Berge

!--------------------------------------------------------------------

subroutine wet_removal_Berge_mod(mgmzp,m1,dt,dn0,wdepmass,precip, &
     stcum1d,sc_up,k22,kbcon,ktop,z_cup, &
     p_lw,cd,zu)

  implicit none
  integer k,mgmzp,m1,k22,kbcon,ktop
  real scav_eff,coll_eff,A,v,Mw,dz,dt,wdepmass,tEND_wd,precip
  real Maer,Idetr,M_I,tEND_wd_detr,detup
  real, dimension(mgmzp) :: sc_up,stcum1d,z_cup,p_lw,cd,zu,dn0
  !REAL, DIMENSION(m1)    :: dn0
  data scav_eff/0.6/    ! for smoke : Chuang et al. (1992) J. Atmos. Sci.
  data coll_eff/0.01/   ! the mean collection efficiency averaged over all raindrop water
  data A/5.2/           ! constant, unit : m^3/(kg s)
  data v/5./            ! raindrop fall speed (m/s)

  !------------ In-cloud wet deposition ' washout'----
  ! Following Berge, TELLUS, 1993, eq. 3.7a
  ! ModIFicaDO para colocar a remocao somente na regiao de desentranhamento
  Maer  = 0.
  Idetr = 0.

  do  k=kbcon,ktop

     !Units:
     ! sc_up (in-cloud aerosol mixing ratio)   : kg[aer]   / kg[air]
     ! p_lw  (precip rate/ liq water content)  :(kg[water] / (m^2 s)) / (kg[water]/kg[air])
     !                                           = kg[air]/( m^2 s)
     ! dno   (basic state air density)         : kg[air]   / m^3 (on the Grell's grid)
     ! dz    (thickness of vertical layer)     : m
     ! tend_wd (wet dep tendency) - original   : kg[aer]/ (kg[air] s)
     ! tend_wd_detr (wet dep tendency) - modif : kg[aer]/ (kg[air] s)
     ! Maer (in-cloud aerosol mass removed)    : kg[aer] / (m^2 s )
     ! Idetr (in-cloud detrainment integral)   : kg[air]/m^3
     ! cd (detrainment rate)                   : 1/m
     dz         = z_cup(k+1) - z_cup(k)
     tend_wd    = scav_eff*sc_up(k)*p_lw(k) /( dn0(k)*dz )

     !total mass removed:
     Maer       = Maer + tEND_wd*dn0(k)*dz
     !detrainment vertiCALLy integrated
     if(k.lt.ktop)   then
        detup =   cd(k+1) *dz*zu(k)
     else
        detup =   zu(k)
     endif
     ! Idetr      = Idetr + cd(k)*dn0(k)*dz
     Idetr      = Idetr + detup*dn0(k)*dz
  enddo

  M_I = Maer/(Idetr + 1.e-16)

  do  k=kbcon,ktop
     dz         = z_cup(k+1) - z_cup(k)
     !reallocated wet remotion tendency
     if(k.lt.ktop)   then
        detup =   cd(k+1) *dz*zu(k)
     else
        detup =   zu(k)
     endif
     ! tend_wd_detr = M_I*cd(k)
     tend_wd_detr = M_I*detup

     !Accumulate the deposited mass on surface.
     wdepmass   = wdepmass + dt * tend_wd_detr*dn0(k)*dz

     !Final tendency = convective transport - wet deposition removal
     stcum1d(k) = stcum1d(k) - tend_wd_detr
  enddo


  !--tmp------
  !PRINT*,'WET2',k,kbcon,ktop,wdepmass
  return
  !--tmp------


  !------------ Below-cloud wet deposition 'rainout'----
  ! Following Berge, TELLUS, 1993, eq. 3.8
  !Units:
  ! sc_up (below cloud aerosol mixing ratio): kg[aer]   / kg[air]
  ! precip  (total precip rate)             :(kg[water] / (m^2 s)) 
  ! Mw (water precipitation concentration   : kg[water]/m^3 
  ! tend_wd (wet dep tendency)              : kg[aer]/ (kg[air] s)

  Mw = precip/v

  do  k=1,kbcon-1

     dz         = z_cup(k+1) - z_cup(k)
     tend_wd    = A * Mw * coll_eff * sc_up(k)

     !Accumulate the deposited mass on surface.
     wdepmass   = wdepmass + dt * tend_wd*dn0(k)*dz

     !Final tendency = convective transport - wet deposition removal
     stcum1d(k) = stcum1d(k) - tend_wd
  enddo

end subroutine wet_removal_Berge_mod
!-------------------------------------------------------------------------------

subroutine print1d(mgmzp,n1,mynum,i,j,k22,kbcon,ktop,kpbl,xmb,edt, &
     jmin,kdet,stcum1d,se,se_cup,sc_up,sc_dn)

  implicit none

  ! Arguments
  integer                :: mgmzp, n1, mynum, i, j, k22, kbcon, ktop, kpbl, jmin, kdet
  real                   :: xmb, edt
  real, dimension(mgmzp) :: se, se_cup, sc_up, sc_dn, stcum1d

  ! Local variables
  integer :: k
  real    :: pmar, pmco, pmco2, fi, fcu, f_CO

  data PMAR/28.96/
  data PMCO/28./
  data PMCO2/44./

  fi=((PMAR/PMCO)*1.E+9)  ! Converte ppb <==> mg/kg
  fcu =1.e-6              ! mg [gas/part] ==> /kg [ar]
  f_CO=fcu*fi
  print*,'============================================================='
  print*,'MYNUM   I   J  =',mynum,i,j
  print*,'K22 KBCON KTOP =',K22,kbcon,ktop
  print*,'JMIN KDET KPBL =',JMIN,KDET,kpbl
  print*,'MASS FLUX - EDT=',xmb,edt
  do k=n1,1,-1
     print*,k,f_CO*se(k),f_CO*sc_up(k),f_CO*sc_dn(k),f_CO*86400.*stcum1d(k)
  enddo
  print*,'============================================================='

end subroutine print1d

!-------------------------------------------------------------------------------
subroutine update_gf(nspecies,mgmzp,ngrid,n1,n2,n3,i,j,se,se_cup,mynum,stcum1d,dt) 
!  use chem1_list, only: spc_alloc_chem  =>spc_alloc,transport,on,off
  use mem_chem1,  only: chem1_g
  use mem_aer1,  only: aer1_g
  use aer1_list, only: nmodes
 
  use mem_tconv, only : nchem_a,nchem_z, ind_chem &
                       ,naer_a,naer_z,ind_aer, ind_mode
  
  implicit none
  integer :: nspecies,ngrid,mgmzp, n1, n2, n3, i, j,mynum
  real    :: se(nspecies,mgmzp), se_cup(nspecies,mgmzp),dt
  !local Variables
  integer :: k, kr,ispc
  real, dimension(nspecies,mgmzp) :: stcum1d

  do k=2,n1-1
     kr = k+1   ! level k of conv param corresponds to level K + 1 of rams.

     !-- chemistry section
     do ispc = nchem_a,nchem_z
       
       !- mixing ratio on Z (model) levels
       se(ispc,k) = se(ispc,k) + dt *  stcum1d(ispc,k)           
       
       if(se(ispc,k) < -0.5) then
          write(mynum,fmt='(A,5(I5.5,1X))') 'CHEM i, j,mynum,k,ispc',i, j,mynum,k,ispc
          write(mynum,fmt='(2(E18.5,1X))') se(ispc,k),stcum1d(ispc,k)
          flush(mynum)
         ! stop 444
       endif  
	
       
       
!     enddo
  
     enddo
     
     !aerosol section
     
     do ispc = naer_a,naer_z
            
      
     ! if(ispc==34) print*,'SE=',se(ispc,k),se(ispc,k)+dt*stcum1d(ispc,k)  
     
       se    (ispc,k) =  se    (ispc,k)  + dt *  stcum1d(ispc,k)           ! mixing ratio on Z levels
    
        if(se(ispc,k) < -0.5) then
          write(mynum,fmt='(A,5(I5.5,1X))') 'AER i, j,mynum,k,ispc',i, j,mynum,k,ispc
          write(mynum,fmt='(2(E18.5,1X))') se(ispc,k),stcum1d(ispc,k)
          flush(mynum)
         ! stop 500
       endif  
      

     enddo
     
  enddo


end subroutine update_gf
!---------------------------------------------------------------------------
