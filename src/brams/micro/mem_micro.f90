!
! Copyright (C) 1991-2004  ; All Rights Reserved ; Colorado State University
! Colorado State University Research Foundation ; ATMET, LLC
! 
! This file is free software; you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software 
! Foundation; either version 2 of the License, or (at your option) any later version.
! 
! This software is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with this 
! program; if not, write to the Free Software Foundation, Inc., 
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!======================================================================================
!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################


Module mem_micro

   Type micro_vars
   
      ! Variables to be dimensioned by (nzp,nxp,nyp)
   real, pointer, dimension(:,:,:) :: &
                          rcp,rdp,rrp,rpp,rsp,rap,rgp,rhp &
                         ,ccp,cdp,crp,cpp,csp,cap,cgp,chp &
                         ,cccnp,gccnp,cifnp,q2,q6,q7 &
                         ,rei,rel,cldfr                 &
			                ,cccmp,gccmp,cnm1p,cnm2p,cnm3p,cnm8p &
                         ,md1np,md2np,salt_filmp,salt_jetp,salt_spmp  &
                         ,pcpvr,pcpvp,pcpvs,pcpva,pcpvg,pcpvh,pcpvd   &
           !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES
           ,nuccldr,nuccldc,nucicer,nucicec,inuchomr                  &
           ,inuchomc,inuccontr,inuccontc,inucifnr,inucifnc            &
           ,inuchazr,inuchazc,vapliq,vapice,vapcld                    &
           ,vaprain,vappris,vapsnow,vapaggr,vapgrau                   &
           ,vaphail,vapdriz,meltice,meltpris,meltsnow                 &
           ,meltaggr,meltgrau,melthail,cld2rain,rimecld               &
           ,rimecldsnow,rimecldaggr,rimecldgrau,rimecldhail,rain2ice  &
           ,rain2pr,rain2sn,rain2ag,rain2gr,rain2ha                   &
           ,rain2ha_xtra,ice2rain,aggregate,aggrselfpris,aggrselfsnow &
           ,aggrprissnow,latheatvap,latheatfrz                        &
           !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES (totals)
           ,nuccldrt,nuccldct,nucicert,nucicect,inuchomrt                  &
           ,inuchomct,inuccontrt,inuccontct,inucifnrt,inucifnct            &
           ,inuchazrt,inuchazct,vapliqt,vapicet,vapcldt                    &
           ,vapraint,vapprist,vapsnowt,vapaggrt,vapgraut                   &
           ,vaphailt,vapdrizt,melticet,meltprist,meltsnowt                 &
           ,meltaggrt,meltgraut,melthailt,cld2raint,rimecldt               &
           ,rimecldsnowt,rimecldaggrt,rimecldgraut,rimecldhailt,rain2icet  &
           ,rain2prt,rain2snt,rain2agt,rain2grt,rain2hat                   &
           ,rain2ha_xtrat,ice2raint,aggregatet,aggrselfprist,aggrselfsnowt &
           ,aggrprissnowt,latheatvapt,latheatfrzt


      ! Variables to be dimensioned by (nnxp,nyp)
   real, pointer, dimension(:,:) :: &
                          accpr,accpp,accps,accpa,accpg,accph,accpd &
                         ,pcprr,pcprp,pcprs,pcpra,pcprg,pcprh,pcprd &
                         ,pcpg,qpcpg,dpcpg 
                          
   End Type               
                          
   type (micro_vars), allocatable :: micro_g(:), microm_g(:)
                          
Contains                  
                          
   subroutine alloc_micro(micro,n1,n2,n3,ng)
   
   USE micphys, only : level,idriz,irain,ipris,isnow,igraup,ihail,jnmb,&
                       icloud,iccnlev,idust,imd1flg,imd2flg, isalt,&
		                 imbudget,imbudtot,iaggr,mcphys_type
   USE mem_radiate, ONLY: ilwrtyp, iswrtyp       ! INTENT(IN)
   USE mem_cuparm , ONLY: nnqparm                ! INTENT(IN)

   implicit none          
   type (micro_vars) :: micro
   integer, intent(in) :: n1,n2,n3,ng

! Allocate arrays based on options (if necessary)
   
   ! gthompson microphysics, GFDL and WSM
   IF(mcphys_type == 2 .or. mcphys_type == 3 .or. mcphys_type == 4 .or. &
      mcphys_type == 5 .or. mcphys_type == 6 .or. mcphys_type == 7) then 
     
          level = 3
          idriz = 0 ; icloud = 1
          irain = 1 ; ipris  = 1 
          isnow = 1 ; iaggr  = 0
          jnmb(1) = 1 !cloud
          jnmb(2) = 1 !rain
          jnmb(3) = 1 !pristine
          jnmb(4) = 1 !snow
          jnmb(5) = 0 !agg

          igraup = 1 !graupel
          IF(mcphys_type == 5) igraup = 0 
          jnmb(6) = igraup !graupel
          
          ihail = 0
          IF(mcphys_type == 7) ihail = 1
          jnmb(7) =  ihail
          
          jnmb(8) = 0 !idriz

          !- cloud liq water
          allocate (micro%rcp(n1,n2,n3))   ;micro%rcp  =0.0
          
          !- rain
          allocate (micro%rrp  (n1,n2,n3)) ;micro%rrp  =0.0
          !- for this scheme, the rain rate below will
	       !- account for rain+ice+snow+graupel
	       allocate (micro%accpr(n2,n3)) ;   micro%accpr=0.0
          allocate (micro%pcprr(n2,n3)) ;   micro%pcprr=0.0
      
          !- ice
          allocate (micro%rpp  (n1,n2,n3)) ;micro%rpp  =0.0
          !- don t need to be allocated, see coments above
	       !allocate (micro%accpp(n2,n3)) ;   micro%accpp=0.0
          !allocate (micro%pcprp(n2,n3)) ;   micro%pcprp=0.0
      
          !- snow
          allocate (micro%rsp  (n1,n2,n3)) ;micro%rsp	=0.0
          !- the rates bellow will account for snow and ice
	       allocate (micro%accps(n2,n3)) ;   micro%accps =0.0
          allocate (micro%pcprs(n2,n3)) ;   micro%pcprs =0.0
      
          !- graupel
          IF(mcphys_type /= 5)  then 
            allocate (micro%rgp  (n1,n2,n3)) ;micro%rgp  =0.0
            !- the rates bellow will account for only graupel
	         allocate (micro%accpg(n2,n3)) ;   micro%accpg=0.0
            allocate (micro%pcprg(n2,n3)) ;   micro%pcprg=0.0
          ENDIF
          !- hail
          IF(mcphys_type == 7) then
            allocate (micro%rhp  (n1,n2,n3)) ;micro%rhp  =0.0
            allocate (micro%accph(n2,n3)) ;   micro%accph=0.0
            allocate (micro%pcprh(n2,n3)) ;   micro%pcprh=0.0
            !allocate (micro%pcpvh(n1,n2,n3)) ;micro%pcpvh=0.0
            !allocate (micro%q7   (n1,n2,n3)) ;micro%q7   =0.0    
          ENDIF
      
	       IF(mcphys_type  == 2 .or. mcphys_type  == 3) then ! only for double-moment and 
            !- number concentration for cloud/rain/ice
            !- obs : ccp don t need to be allocated for the single-moment
	         !- cloud water scheme (the same for CCN and IFN).
            allocate(micro%crp  (n1,n2,n3)) ;micro%crp  =0.0 
            allocate(micro%cpp  (n1,n2,n3)) ;micro%cpp  =0.0 
!---these should not be allocated for mcphys_type  == 2 because
!---they are not used for this option
!            !ST
!	          allocate(micro%ccp  (n1,n2,n3)) ;micro%ccp  =0.0 
!            allocate(micro%cccnp(n1,n2,n3)) ;micro%cccnp=0.0 !;endif 
!            allocate(micro%cifnp(n1,n2,n3)) ;micro%cifnp=0.0 !;endif 
!            !ST
          ENDIF
         !- only for cloud water double-moment and aerosol aware microphysics         
	       IF(mcphys_type  == 3) then ! only for double-moment and 
	        allocate(micro%ccp  (n1,n2,n3)) ;micro%ccp  =0.0 
           allocate(micro%cccnp(n1,n2,n3)) ;micro%cccnp=0.0 !;endif 
           allocate(micro%cifnp(n1,n2,n3)) ;micro%cifnp=0.0 !;endif 
          ENDIF
          
          !- 3D cloud fraction from GFDL cloud microphysics and GF convection
          IF(mcphys_type  == 4 .or. nnqparm(ng) == 8) then 
           allocate(micro%cldfr  (n1,n2,n3)) ;micro%cldfr  =0.0 
          ENDIF
      
     !- for consistency with the other parts of BRAMS
	  !- pgcp will be the total precipitation rate
	        allocate (micro%pcpg (n2,n3));micro%pcpg =0.0
          !-the allocations below are tmp for leaf-3
	        allocate (micro%qpcpg(n2,n3));micro%qpcpg=0.0         
	        allocate (micro%dpcpg(n2,n3));micro%dpcpg=0.0
      
     !- allocation of memory for effective radius for RRTMG
	       IF(ilwrtyp==6 .or. iswrtyp==6 ) THEN
           allocate (micro%rei  (n1,n2,n3)) ;micro%rei  =0.0  
	        allocate (micro%rel  (n1,n2,n3)) ;micro%rel  =0.0
	       ENDIF
      
    ELSE  ! for the traditional RAMS microphysics
 
      if (level >= 2 ) then
            allocate (micro%rcp(n1,n2,n3))   ;micro%rcp    =0.0
      endif
      if (level >= 3) then         
	      if(irain >= 1)  then
            allocate (micro%rrp  (n1,n2,n3)) ;micro%rrp  =0.0
            allocate (micro%accpr(n2,n3)) ;   micro%accpr=0.0
            allocate (micro%pcprr(n2,n3)) ;   micro%pcprr=0.0
            allocate (micro%pcpvr(n1,n2,n3)) ;micro%pcpvr=0.0
            allocate (micro%q2   (n1,n2,n3)) ;micro%q2   =0.0
         endif
         if(ipris >= 1)  then
            allocate (micro%rpp  (n1,n2,n3)) ;micro%rpp  =0.0
            allocate (micro%accpp(n2,n3)) ;   micro%accpp=0.0
            allocate (micro%pcprp(n2,n3)) ;   micro%pcprp=0.0
            allocate (micro%pcpvp(n1,n2,n3)) ;micro%pcpvp=0.0
         endif
         if(isnow >= 1)  then
            allocate (micro%rsp  (n1,n2,n3)) ;micro%rsp   =0.0
            allocate (micro%accps(n2,n3)) ;   micro%accps =0.0
            allocate (micro%pcprs(n2,n3)) ;   micro%pcprs =0.0
            allocate (micro%pcpvs(n1,n2,n3)) ;micro%pcpvs =0.0
         endif
         if(iaggr >= 1)  then
            allocate (micro%rap  (n1,n2,n3)) ;micro%rap  =0.0
            allocate (micro%accpa(n2,n3)) ;   micro%accpa=0.0
            allocate (micro%pcpra(n2,n3)) ;   micro%pcpra=0.0
            allocate (micro%pcpva(n1,n2,n3)) ;micro%pcpva=0.0
         endif						 
         if(igraup >= 1) then
            allocate (micro%rgp  (n1,n2,n3)) ;micro%rgp  =0.0
            allocate (micro%accpg(n2,n3)) ;   micro%accpg=0.0
            allocate (micro%pcprg(n2,n3)) ;   micro%pcprg=0.0
            allocate (micro%pcpvg(n1,n2,n3)) ;micro%pcpvg=0.0
            allocate (micro%q6   (n1,n2,n3)) ;micro%q6   =0.0
         endif
         if(ihail >= 1)  then
            allocate (micro%rhp  (n1,n2,n3)) ;micro%rhp  =0.0
            allocate (micro%accph(n2,n3)) ;   micro%accph=0.0
            allocate (micro%pcprh(n2,n3)) ;   micro%pcprh=0.0
            allocate (micro%pcpvh(n1,n2,n3)) ;micro%pcpvh=0.0
            allocate (micro%q7   (n1,n2,n3)) ;micro%q7   =0.0
         endif
         if(jnmb(1) >= 5)  then; allocate(micro%ccp  (n1,n2,n3)) ;micro%ccp  =0.0 ;endif
         if(jnmb(2) == 5)  then; allocate(micro%crp  (n1,n2,n3)) ;micro%crp  =0.0 ;endif
         if(jnmb(3) >= 5)  then; allocate(micro%cpp  (n1,n2,n3)) ;micro%cpp  =0.0 ;endif
         if(jnmb(4) == 5)  then; allocate(micro%csp  (n1,n2,n3)) ;micro%csp  =0.0 ;endif
         if(jnmb(5) == 5)  then; allocate(micro%cap  (n1,n2,n3)) ;micro%cap  =0.0 ;endif
         if(jnmb(6) == 5)  then; allocate(micro%cgp  (n1,n2,n3)) ;micro%cgp  =0.0 ;endif
         if(jnmb(7) == 5)  then; allocate(micro%chp  (n1,n2,n3)) ;micro%chp  =0.0 ;endif
         if(icloud  >= 5)  then; allocate(micro%cccnp(n1,n2,n3)) ;micro%cccnp=0.0 ;endif 
         if(ipris   >= 5)  then; allocate(micro%cifnp(n1,n2,n3)) ;micro%cifnp=0.0 ;endif 
         
         if(icloud >= 5)   then; allocate(micro%cccmp(n1,n2,n3)) ;micro%cccmp=0.0 ;endif
  	 
	 allocate (micro%pcpg (n2,n3));micro%pcpg =0.0
  	 allocate (micro%qpcpg(n2,n3));micro%qpcpg=0.0
  	 allocate (micro%dpcpg(n2,n3));micro%dpcpg=0.0
 
         !- only for 2M microphysics
         if(mcphys_type == 1)  then
             if(idriz >= 1 )  then
               allocate (micro%rdp  (n1,n2,n3)) ;micro%rdp  =0.0
               allocate (micro%accpd(n2,n3))    ;micro%accpd=0.0
               allocate (micro%pcprd(n2,n3))    ;micro%pcprd=0.0
               allocate (micro%pcpvd(n1,n2,n3)) ;micro%pcpvd=0.0
             endif

  	     if(jnmb(8) >= 5)  then; allocate(micro%cdp  (n1,n2,n3)) ;micro%cdp  =0.0 ;endif
  	     if(idriz	>= 5)  then; allocate(micro%gccnp(n1,n2,n3)) ;micro%gccnp=0.0 ;endif 
  	     if(idriz   >= 5)  then; allocate(micro%gccmp(n1,n2,n3)) ;micro%gccmp=0.0 ;endif
  	     
  	     if(iccnlev >= 2 .and. jnmb(1) >= 5)  then; allocate(micro%cnm1p(n1,n2,n3)) ;micro%cnm1p=0.0 ;endif
  	     if(iccnlev >= 2 .and. jnmb(2) >= 1)  then; allocate(micro%cnm2p(n1,n2,n3)) ;micro%cnm2p=0.0 ;endif
  	     if(iccnlev >= 2 .and. jnmb(3) >= 1)  then; allocate(micro%cnm3p(n1,n2,n3)) ;micro%cnm3p=0.0 ;endif
  	     if(iccnlev >= 2 .and. jnmb(8) >= 1)  then; allocate(micro%cnm8p(n1,n2,n3)) ;micro%cnm8p=0.0 ;endif
  	     
  	     if(idust == 1 .or. imd1flg == 1)	  then; allocate(micro%md1np(n1,n2,n3)) ;micro%md1np=0.0 ;endif
  	     if(idust == 1 .or. imd2flg == 1)	  then; allocate(micro%md2np(n1,n2,n3)) ;micro%md2np=0.0 ;endif
  	     if(isalt == 1) then; allocate(micro%salt_filmp(n1,n2,n3))  ;micro%salt_filmp =0.0  	 ;endif
  	     if(isalt == 1) then; allocate(micro%salt_jetp (n1,n2,n3))  ;micro%salt_jetp  =0.0  	 ;endif
  	     if(isalt == 1) then; allocate(micro%salt_spmp (n1,n2,n3))  ;micro%salt_spmp  =0.0  	 ;endif

 		
  	     !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES
  	     if(imbudget>=1 .or. imbudtot>=1) then
  	       allocate (micro%latheatvap(n1,n2,n3));micro%latheatvap=0.0
  	       allocate (micro%latheatfrz(n1,n2,n3));micro%latheatfrz=0.0
  	     endif
  	     if(imbudget>=1) then
  	       allocate (micro%nuccldr  (n1,n2,n3));micro%nuccldr  =0.0
  	       allocate (micro%nuccldc  (n1,n2,n3));micro%nuccldc  =0.0
  	       allocate (micro%cld2rain (n1,n2,n3));micro%cld2rain =0.0
  	       allocate (micro%ice2rain (n1,n2,n3));micro%ice2rain =0.0
  	       allocate (micro%nucicer  (n1,n2,n3));micro%nucicer  =0.0
  	       allocate (micro%nucicec  (n1,n2,n3));micro%nucicec  =0.0
  	       allocate (micro%vapliq	(n1,n2,n3));micro%vapliq   =0.0 	  
  	       allocate (micro%vapice	(n1,n2,n3));micro%vapice   =0.0
  	       allocate (micro%meltice  (n1,n2,n3));micro%meltice  =0.0
  	       allocate (micro%rimecld  (n1,n2,n3));micro%rimecld  =0.0
  	       allocate (micro%rain2ice (n1,n2,n3));micro%rain2ice =0.0
  	       allocate (micro%aggregate(n1,n2,n3));micro%aggregate=0.0
  	     endif
  	     if(imbudget==2) then
  	       allocate (micro%inuchomr    (n1,n2,n3)) ;micro%inuchomr    =0.0 
  	       allocate (micro%inuchomc    (n1,n2,n3)) ;micro%inuchomc    =0.0
  	       allocate (micro%inuccontr   (n1,n2,n3)) ;micro%inuccontr   =0.0
  	       allocate (micro%inuccontc   (n1,n2,n3)) ;micro%inuccontc   =0.0
  	       allocate (micro%inucifnr    (n1,n2,n3)) ;micro%inucifnr    =0.0
  	       allocate (micro%inucifnc    (n1,n2,n3)) ;micro%inucifnc    =0.0
  	       allocate (micro%inuchazr    (n1,n2,n3)) ;micro%inuchazr    =0.0
  	       allocate (micro%inuchazc    (n1,n2,n3)) ;micro%inuchazc    =0.0
  	       allocate (micro%vapcld	   (n1,n2,n3)) ;micro%vapcld	  =0.0
  	       allocate (micro%vaprain     (n1,n2,n3)) ;micro%vaprain	  =0.0
  	       allocate (micro%vappris     (n1,n2,n3)) ;micro%vappris	  =0.0
  	       allocate (micro%vapsnow     (n1,n2,n3)) ;micro%vapsnow	  =0.0   
  	       allocate (micro%vapaggr     (n1,n2,n3)) ;micro%vapaggr	  =0.0
  	       allocate (micro%vapgrau     (n1,n2,n3)) ;micro%vapgrau	  =0.0
  	       allocate (micro%vaphail     (n1,n2,n3)) ;micro%vaphail	  =0.0
  	       allocate (micro%vapdriz     (n1,n2,n3)) ;micro%vapdriz	  =0.0
  	       allocate (micro%meltpris    (n1,n2,n3)) ;micro%meltpris    =0.0
  	       allocate (micro%meltsnow    (n1,n2,n3)) ;micro%meltsnow    =0.0
  	       allocate (micro%meltaggr    (n1,n2,n3)) ;micro%meltaggr    =0.0
  	       allocate (micro%meltgrau    (n1,n2,n3)) ;micro%meltgrau    =0.0
  	       allocate (micro%melthail    (n1,n2,n3)) ;micro%melthail    =0.0
  	       allocate (micro%rimecldsnow (n1,n2,n3)) ;micro%rimecldsnow =0.0
  	       allocate (micro%rimecldaggr (n1,n2,n3)) ;micro%rimecldaggr =0.0
  	       allocate (micro%rimecldgrau (n1,n2,n3)) ;micro%rimecldgrau =0.0
  	       allocate (micro%rimecldhail (n1,n2,n3)) ;micro%rimecldhail =0.0
  	       allocate (micro%rain2pr     (n1,n2,n3)) ;micro%rain2pr	  =0.0
  	       allocate (micro%rain2sn     (n1,n2,n3)) ;micro%rain2sn	  =0.0
  	       allocate (micro%rain2ag     (n1,n2,n3)) ;micro%rain2ag	  =0.0
  	       allocate (micro%rain2gr     (n1,n2,n3)) ;micro%rain2gr	  =0.0
  	       allocate (micro%rain2ha     (n1,n2,n3)) ;micro%rain2ha	  =0.0
  	       allocate (micro%rain2ha_xtra(n1,n2,n3)) ;micro%rain2ha_xtra=0.0
  	       allocate (micro%aggrselfpris(n1,n2,n3)) ;micro%aggrselfpris=0.0
  	       allocate (micro%aggrselfsnow(n1,n2,n3)) ;micro%aggrselfsnow=0.0
  	       allocate (micro%aggrprissnow(n1,n2,n3)) ;micro%aggrprissnow=0.0
  	     endif
  	     !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES (totals)
  	     if(imbudtot>=1) then
  	       allocate (micro%nuccldrt   (n1,n2,n3)) ;micro%nuccldrt	=0.0
  	       allocate (micro%nuccldct   (n1,n2,n3)) ;micro%nuccldct	=0.0
  	       allocate (micro%cld2raint  (n1,n2,n3)) ;micro%cld2raint  =0.0
  	       allocate (micro%ice2raint  (n1,n2,n3)) ;micro%ice2raint  =0.0
  	       allocate (micro%nucicert   (n1,n2,n3)) ;micro%nucicert	=0.0
  	       allocate (micro%nucicect   (n1,n2,n3)) ;micro%nucicect	=0.0
  	       allocate (micro%vapliqt    (n1,n2,n3)) ;micro%vapliqt	=0.0
  	       allocate (micro%vapicet    (n1,n2,n3)) ;micro%vapicet	=0.0
  	       allocate (micro%melticet   (n1,n2,n3)) ;micro%melticet	=0.0
  	       allocate (micro%rimecldt   (n1,n2,n3)) ;micro%rimecldt	=0.0
  	       allocate (micro%rain2icet  (n1,n2,n3)) ;micro%rain2icet  =0.0
  	       allocate (micro%aggregatet (n1,n2,n3)) ;micro%aggregatet =0.0
  	       allocate (micro%latheatvapt(n1,n2,n3)) ;micro%latheatvapt=0.0
  	       allocate (micro%latheatfrzt(n1,n2,n3)) ;micro%latheatfrzt=0.0
  	     endif							
  	     if(imbudtot==2) then
  	       allocate (micro%inuchomrt(n1,n2,n3)) ;	  micro%inuchomrt    =0.0
  	       allocate (micro%inuchomct(n1,n2,n3)) ;	  micro%inuchomct    =0.0
  	       allocate (micro%inuccontrt(n1,n2,n3)) ;    micro%inuccontrt   =0.0
  	       allocate (micro%inuccontct(n1,n2,n3)) ;    micro%inuccontct   =0.0
  	       allocate (micro%inucifnrt(n1,n2,n3)) ;	  micro%inucifnrt    =0.0
  	       allocate (micro%inucifnct(n1,n2,n3)) ;	  micro%inucifnct    =0.0
  	       allocate (micro%inuchazrt(n1,n2,n3)) ;	  micro%inuchazrt    =0.0
  	       allocate (micro%inuchazct(n1,n2,n3)) ;	  micro%inuchazct    =0.0
  	       allocate (micro%vapcldt(n1,n2,n3)) ;	  micro%vapcldt      =0.0
  	       allocate (micro%vapraint(n1,n2,n3)) ;	  micro%vapraint     =0.0
  	       allocate (micro%vapprist(n1,n2,n3)) ;	  micro%vapprist     =0.0
  	       allocate (micro%vapsnowt(n1,n2,n3)) ;	  micro%vapsnowt     =0.0
  	       allocate (micro%vapaggrt(n1,n2,n3)) ;	  micro%vapaggrt     =0.0
  	       allocate (micro%vapgraut(n1,n2,n3)) ;	  micro%vapgraut     =0.0
  	       allocate (micro%vaphailt(n1,n2,n3)) ;	  micro%vaphailt     =0.0
  	       allocate (micro%vapdrizt(n1,n2,n3)) ;	  micro%vapdrizt     =0.0
  	       allocate (micro%meltprist(n1,n2,n3)) ;	  micro%meltprist    =0.0
  	       allocate (micro%meltsnowt(n1,n2,n3)) ;	  micro%meltsnowt    =0.0
  	       allocate (micro%meltaggrt(n1,n2,n3)) ;	  micro%meltaggrt    =0.0
  	       allocate (micro%meltgraut(n1,n2,n3)) ;	  micro%meltgraut    =0.0
  	       allocate (micro%melthailt(n1,n2,n3)) ;	  micro%melthailt    =0.0
  	       allocate (micro%rimecldsnowt(n1,n2,n3)) ;  micro%rimecldsnowt =0.0
  	       allocate (micro%rimecldaggrt(n1,n2,n3)) ;  micro%rimecldaggrt =0.0
  	       allocate (micro%rimecldgraut(n1,n2,n3)) ;  micro%rimecldgraut =0.0
  	       allocate (micro%rimecldhailt(n1,n2,n3)) ;  micro%rimecldhailt =0.0
  	       allocate (micro%rain2prt(n1,n2,n3)) ;	  micro%rain2prt     =0.0
  	       allocate (micro%rain2snt(n1,n2,n3)) ;	  micro%rain2snt     =0.0
  	       allocate (micro%rain2agt(n1,n2,n3)) ;	  micro%rain2agt     =0.0
  	       allocate (micro%rain2grt(n1,n2,n3)) ;	  micro%rain2grt     =0.0
  	       allocate (micro%rain2hat(n1,n2,n3)) ;	  micro%rain2hat     =0.0
  	       allocate (micro%rain2ha_xtrat(n1,n2,n3)) ; micro%rain2ha_xtrat=0.0
  	       allocate (micro%aggrselfprist(n1,n2,n3)) ; micro%aggrselfprist=0.0
  	       allocate (micro%aggrselfsnowt(n1,n2,n3)) ; micro%aggrselfsnowt=0.0
  	       allocate (micro%aggrprissnowt(n1,n2,n3)) ; micro%aggrprissnowt=0.0
  	     endif							     
  	 endif	! mcphys_type=1							     
          !- allocation of memory for effective radius for RRTMG
	  IF(ilwrtyp==6 .or. iswrtyp==6 ) THEN
       allocate (micro%rei  (n1,n2,n3)) ;micro%rei  =0.0  
	    allocate (micro%rel  (n1,n2,n3)) ;micro%rel  =0.0
	  ENDIF
      endif ! level >=3							 
   ENDIF							 
   return								 
   end subroutine alloc_micro							 
									 
									 
   subroutine nullify_micro(micro)					 
									 
   implicit none
   type (micro_vars) :: micro
   
   if (associated(micro%rcp))     nullify (micro%rcp)
   if (associated(micro%rdp))     nullify (micro%rdp)
   if (associated(micro%rrp))     nullify (micro%rrp)
   if (associated(micro%rpp))     nullify (micro%rpp)
   if (associated(micro%rsp))     nullify (micro%rsp)
   if (associated(micro%rap))     nullify (micro%rap)
   if (associated(micro%rgp))     nullify (micro%rgp)
   if (associated(micro%rhp))     nullify (micro%rhp)
   if (associated(micro%ccp))     nullify (micro%ccp)
   if (associated(micro%cdp))     nullify (micro%cdp)
   if (associated(micro%crp))     nullify (micro%crp)
   if (associated(micro%cpp))     nullify (micro%cpp)
   if (associated(micro%csp))     nullify (micro%csp)
   if (associated(micro%cap))     nullify (micro%cap)
   if (associated(micro%cgp))     nullify (micro%cgp)
   if (associated(micro%chp))     nullify (micro%chp)
   if (associated(micro%cccnp))   nullify (micro%cccnp)
   if (associated(micro%gccnp))   nullify (micro%gccnp)
   if (associated(micro%cifnp))   nullify (micro%cifnp)
   if (associated(micro%q2))      nullify (micro%q2)
   if (associated(micro%q6))      nullify (micro%q6)
   if (associated(micro%q7))      nullify (micro%q7)

   if (associated(micro%cccmp))   nullify (micro%cccmp)
   if (associated(micro%gccmp))   nullify (micro%gccmp)
   if (associated(micro%cnm1p))   nullify (micro%cnm1p)
   if (associated(micro%cnm2p))   nullify (micro%cnm2p)
   if (associated(micro%cnm3p))   nullify (micro%cnm3p)
   if (associated(micro%cnm8p))   nullify (micro%cnm8p)
   if (associated(micro%md1np))   nullify (micro%md1np)
   if (associated(micro%md2np))   nullify (micro%md2np)
   if (associated(micro%salt_filmp)) nullify (micro%salt_filmp)
   if (associated(micro%salt_jetp))  nullify (micro%salt_jetp)
   if (associated(micro%salt_spmp))  nullify (micro%salt_spmp)
   if (associated(micro%pcpvr))   nullify (micro%pcpvr)
   if (associated(micro%pcpvp))   nullify (micro%pcpvp)
   if (associated(micro%pcpvs))   nullify (micro%pcpvs)
   if (associated(micro%pcpva))   nullify (micro%pcpva)
   if (associated(micro%pcpvg))   nullify (micro%pcpvg)
   if (associated(micro%pcpvh))   nullify (micro%pcpvh)
   if (associated(micro%pcpvd))   nullify (micro%pcpvd)

   if (associated(micro%accpr))   nullify (micro%accpr)
   if (associated(micro%accpp))   nullify (micro%accpp)
   if (associated(micro%accps))   nullify (micro%accps)
   if (associated(micro%accpa))   nullify (micro%accpa)
   if (associated(micro%accpg))   nullify (micro%accpg)
   if (associated(micro%accph))   nullify (micro%accph)
   if (associated(micro%accpd))   nullify (micro%accpd)
   if (associated(micro%pcprr))   nullify (micro%pcprr)
   if (associated(micro%pcprp))   nullify (micro%pcprp)
   if (associated(micro%pcprs))   nullify (micro%pcprs)
   if (associated(micro%pcpra))   nullify (micro%pcpra)
   if (associated(micro%pcprg))   nullify (micro%pcprg)
   if (associated(micro%pcprh))   nullify (micro%pcprh)
   if (associated(micro%pcprd))   nullify (micro%pcprd)
   if (associated(micro%pcpg))    nullify (micro%pcpg)
   if (associated(micro%qpcpg))   nullify (micro%qpcpg)
   if (associated(micro%dpcpg))   nullify (micro%dpcpg)

    !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES
    if (associated(micro%nuccldr))      nullify (micro%nuccldr)
    if (associated(micro%nuccldc))      nullify (micro%nuccldc)
    if (associated(micro%nucicer))      nullify (micro%nucicer)
    if (associated(micro%nucicec))      nullify (micro%nucicec)
    if (associated(micro%inuchomr))     nullify (micro%inuchomr)
    if (associated(micro%inuchomc))     nullify (micro%inuchomc)
    if (associated(micro%inuccontr))    nullify (micro%inuccontr)
    if (associated(micro%inuccontc))    nullify (micro%inuccontc)
    if (associated(micro%inucifnr))     nullify (micro%inucifnr)
    if (associated(micro%inucifnc))     nullify (micro%inucifnc)
    if (associated(micro%inuchazr))     nullify (micro%inuchazr)
    if (associated(micro%inuchazc))     nullify (micro%inuchazc)
    if (associated(micro%vapliq))       nullify (micro%vapliq)
    if (associated(micro%vapice))       nullify (micro%vapice)
    if (associated(micro%vapcld))       nullify (micro%vapcld)
    if (associated(micro%vaprain))      nullify (micro%vaprain)
    if (associated(micro%vappris))      nullify (micro%vappris)
    if (associated(micro%vapsnow))      nullify (micro%vapsnow)
    if (associated(micro%vapaggr))      nullify (micro%vapaggr)
    if (associated(micro%vapgrau))      nullify (micro%vapgrau)
    if (associated(micro%vaphail))      nullify (micro%vaphail)
    if (associated(micro%vapdriz))      nullify (micro%vapdriz)
    if (associated(micro%meltice))      nullify (micro%meltice)
    if (associated(micro%meltpris))     nullify (micro%meltpris)
    if (associated(micro%meltsnow))     nullify (micro%meltsnow)
    if (associated(micro%meltaggr))     nullify (micro%meltaggr)
    if (associated(micro%meltgrau))     nullify (micro%meltgrau)
    if (associated(micro%melthail))     nullify (micro%melthail)
    if (associated(micro%cld2rain))     nullify (micro%cld2rain)
    if (associated(micro%rimecld))      nullify (micro%rimecld)
    if (associated(micro%rimecldsnow))  nullify (micro%rimecldsnow)
    if (associated(micro%rimecldaggr))  nullify (micro%rimecldaggr)
    if (associated(micro%rimecldgrau))  nullify (micro%rimecldgrau)
    if (associated(micro%rimecldhail))  nullify (micro%rimecldhail)
    if (associated(micro%rain2ice))     nullify (micro%rain2ice)
    if (associated(micro%rain2pr))      nullify (micro%rain2pr)
    if (associated(micro%rain2sn))      nullify (micro%rain2sn)
    if (associated(micro%rain2ag))      nullify (micro%rain2ag)
    if (associated(micro%rain2gr))      nullify (micro%rain2gr)
    if (associated(micro%rain2ha))      nullify (micro%rain2ha)
    if (associated(micro%rain2ha_xtra)) nullify (micro%rain2ha_xtra)
    if (associated(micro%ice2rain))     nullify (micro%ice2rain)
    if (associated(micro%aggregate))    nullify (micro%aggregate)
    if (associated(micro%aggrselfpris)) nullify (micro%aggrselfpris)
    if (associated(micro%aggrselfsnow)) nullify (micro%aggrselfsnow)
    if (associated(micro%aggrprissnow)) nullify (micro%aggrprissnow)
    if (associated(micro%latheatvap))   nullify (micro%latheatvap)
    if (associated(micro%latheatfrz))   nullify (micro%latheatfrz)

    !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES (totals)
    if (associated(micro%nuccldrt))      nullify (micro%nuccldrt)
    if (associated(micro%nuccldct))      nullify (micro%nuccldct)
    if (associated(micro%nucicert))      nullify (micro%nucicert)
    if (associated(micro%nucicect))      nullify (micro%nucicect)
    if (associated(micro%inuchomrt))     nullify (micro%inuchomrt)
    if (associated(micro%inuchomct))     nullify (micro%inuchomct)
    if (associated(micro%inuccontrt))    nullify (micro%inuccontrt)
    if (associated(micro%inuccontct))    nullify (micro%inuccontct)
    if (associated(micro%inucifnrt))     nullify (micro%inucifnrt)
    if (associated(micro%inucifnct))     nullify (micro%inucifnct)
    if (associated(micro%inuchazrt))     nullify (micro%inuchazrt)
    if (associated(micro%inuchazct))     nullify (micro%inuchazct)
    if (associated(micro%vapliqt))       nullify (micro%vapliqt)
    if (associated(micro%vapicet))       nullify (micro%vapicet)
    if (associated(micro%vapcldt))       nullify (micro%vapcldt)
    if (associated(micro%vapraint))      nullify (micro%vapraint)
    if (associated(micro%vapprist))      nullify (micro%vapprist)
    if (associated(micro%vapsnowt))      nullify (micro%vapsnowt)
    if (associated(micro%vapaggrt))      nullify (micro%vapaggrt)
    if (associated(micro%vapgraut))      nullify (micro%vapgraut)
    if (associated(micro%vaphailt))      nullify (micro%vaphailt)
    if (associated(micro%vapdrizt))      nullify (micro%vapdrizt)
    if (associated(micro%melticet))      nullify (micro%melticet)
    if (associated(micro%meltprist))     nullify (micro%meltprist)
    if (associated(micro%meltsnowt))     nullify (micro%meltsnowt)
    if (associated(micro%meltaggrt))     nullify (micro%meltaggrt)
    if (associated(micro%meltgraut))     nullify (micro%meltgraut)
    if (associated(micro%melthailt))     nullify (micro%melthailt)
    if (associated(micro%cld2raint))     nullify (micro%cld2raint)
    if (associated(micro%rimecldt))      nullify (micro%rimecldt)
    if (associated(micro%rimecldsnowt))  nullify (micro%rimecldsnowt)
    if (associated(micro%rimecldaggrt))  nullify (micro%rimecldaggrt)
    if (associated(micro%rimecldgraut))  nullify (micro%rimecldgraut)
    if (associated(micro%rimecldhailt))  nullify (micro%rimecldhailt)
    if (associated(micro%rain2icet))     nullify (micro%rain2icet)
    if (associated(micro%rain2prt))      nullify (micro%rain2prt)
    if (associated(micro%rain2snt))      nullify (micro%rain2snt)
    if (associated(micro%rain2agt))      nullify (micro%rain2agt)
    if (associated(micro%rain2grt))      nullify (micro%rain2grt)
    if (associated(micro%rain2hat))      nullify (micro%rain2hat)
    if (associated(micro%rain2ha_xtrat)) nullify (micro%rain2ha_xtrat)
    if (associated(micro%ice2raint))     nullify (micro%ice2raint)
    if (associated(micro%aggregatet))    nullify (micro%aggregatet)
    if (associated(micro%aggrselfprist)) nullify (micro%aggrselfprist)
    if (associated(micro%aggrselfsnowt)) nullify (micro%aggrselfsnowt)
    if (associated(micro%aggrprissnowt)) nullify (micro%aggrprissnowt)
    if (associated(micro%latheatvapt))   nullify (micro%latheatvapt)
    if (associated(micro%latheatfrzt))   nullify (micro%latheatfrzt)

   if (associated(micro%rei))     nullify (micro%rei)
   if (associated(micro%rel))     nullify (micro%rel)
   if (associated(micro%cldfr))   nullify (micro%cldfr)
   return
   end subroutine nullify_micro

   subroutine dealloc_micro(micro)

   implicit none
   type (micro_vars) :: micro
   
   if (associated(micro%rcp))     deallocate (micro%rcp)
   if (associated(micro%rdp))     deallocate (micro%rdp)
   if (associated(micro%rrp))     deallocate (micro%rrp)
   if (associated(micro%rpp))     deallocate (micro%rpp)
   if (associated(micro%rsp))     deallocate (micro%rsp)
   if (associated(micro%rap))     deallocate (micro%rap)
   if (associated(micro%rgp))     deallocate (micro%rgp)
   if (associated(micro%rhp))     deallocate (micro%rhp)
   if (associated(micro%ccp))     deallocate (micro%ccp)
   if (associated(micro%cdp))     deallocate (micro%cdp)
   if (associated(micro%crp))     deallocate (micro%crp)
   if (associated(micro%cpp))     deallocate (micro%cpp)
   if (associated(micro%csp))     deallocate (micro%csp)
   if (associated(micro%cap))     deallocate (micro%cap)
   if (associated(micro%cgp))     deallocate (micro%cgp)
   if (associated(micro%chp))     deallocate (micro%chp)
   if (associated(micro%cccnp))   deallocate (micro%cccnp)
   if (associated(micro%gccnp))   deallocate (micro%gccnp)
   if (associated(micro%cifnp))   deallocate (micro%cifnp)
   if (associated(micro%q2))      deallocate (micro%q2)
   if (associated(micro%q6))      deallocate (micro%q6)
   if (associated(micro%q7))      deallocate (micro%q7)

   if (associated(micro%cccmp))   deallocate (micro%cccmp)
   if (associated(micro%gccmp))   deallocate (micro%gccmp)
   if (associated(micro%cnm1p))   deallocate (micro%cnm1p)
   if (associated(micro%cnm2p))   deallocate (micro%cnm2p)
   if (associated(micro%cnm3p))   deallocate (micro%cnm3p)
   if (associated(micro%cnm8p))   deallocate (micro%cnm8p)
   if (associated(micro%md1np))   deallocate (micro%md1np)
   if (associated(micro%md2np))   deallocate (micro%md2np)
   if (associated(micro%salt_filmp)) deallocate (micro%salt_filmp)
   if (associated(micro%salt_jetp))  deallocate (micro%salt_jetp)
   if (associated(micro%salt_spmp))  deallocate (micro%salt_spmp)
   if (associated(micro%pcpvr))   deallocate (micro%pcpvr)
   if (associated(micro%pcpvp))   deallocate (micro%pcpvp)
   if (associated(micro%pcpvs))   deallocate (micro%pcpvs)
   if (associated(micro%pcpva))   deallocate (micro%pcpva)
   if (associated(micro%pcpvg))   deallocate (micro%pcpvg)
   if (associated(micro%pcpvh))   deallocate (micro%pcpvh)
   if (associated(micro%pcpvd))   deallocate (micro%pcpvd)

   if (associated(micro%accpr))   deallocate (micro%accpr)
   if (associated(micro%accpp))   deallocate (micro%accpp)
   if (associated(micro%accps))   deallocate (micro%accps)
   if (associated(micro%accpa))   deallocate (micro%accpa)
   if (associated(micro%accpg))   deallocate (micro%accpg)
   if (associated(micro%accph))   deallocate (micro%accph)
   if (associated(micro%accpd))   deallocate (micro%accpd)
   if (associated(micro%pcprr))   deallocate (micro%pcprr)
   if (associated(micro%pcprp))   deallocate (micro%pcprp)
   if (associated(micro%pcprs))   deallocate (micro%pcprs)
   if (associated(micro%pcpra))   deallocate (micro%pcpra)
   if (associated(micro%pcprg))   deallocate (micro%pcprg)
   if (associated(micro%pcprh))   deallocate (micro%pcprh)
   if (associated(micro%pcprd))   deallocate (micro%pcprd)
   if (associated(micro%pcpg))    deallocate (micro%pcpg)
   if (associated(micro%qpcpg))   deallocate (micro%qpcpg)
   if (associated(micro%dpcpg))   deallocate (micro%dpcpg)

    !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES
    if (associated(micro%nuccldr))      deallocate (micro%nuccldr)
    if (associated(micro%nuccldc))      deallocate (micro%nuccldc)
    if (associated(micro%nucicer))      deallocate (micro%nucicer)
    if (associated(micro%nucicec))      deallocate (micro%nucicec)
    if (associated(micro%inuchomr))     deallocate (micro%inuchomr)     
    if (associated(micro%inuchomc))     deallocate (micro%inuchomc) 
    if (associated(micro%inuccontr))    deallocate (micro%inuccontr)     
    if (associated(micro%inuccontc))    deallocate (micro%inuccontc)
    if (associated(micro%inucifnr))     deallocate (micro%inucifnr)     
    if (associated(micro%inucifnc))     deallocate (micro%inucifnc) 
    if (associated(micro%inuchazr))     deallocate (micro%inuchazr)     
    if (associated(micro%inuchazc))     deallocate (micro%inuchazc)
    if (associated(micro%vapliq))       deallocate (micro%vapliq)
    if (associated(micro%vapice))       deallocate (micro%vapice)
    if (associated(micro%vapcld))       deallocate (micro%vapcld)
    if (associated(micro%vaprain))      deallocate (micro%vaprain)
    if (associated(micro%vappris))      deallocate (micro%vappris)
    if (associated(micro%vapsnow))      deallocate (micro%vapsnow)
    if (associated(micro%vapaggr))      deallocate (micro%vapaggr)
    if (associated(micro%vapgrau))      deallocate (micro%vapgrau)
    if (associated(micro%vaphail))      deallocate (micro%vaphail)
    if (associated(micro%vapdriz))      deallocate (micro%vapdriz)      
    if (associated(micro%meltice))      deallocate (micro%meltice) 
    if (associated(micro%meltpris))     deallocate (micro%meltpris) 
    if (associated(micro%meltsnow))     deallocate (micro%meltsnow) 
    if (associated(micro%meltaggr))     deallocate (micro%meltaggr)
    if (associated(micro%meltgrau))     deallocate (micro%meltgrau) 
    if (associated(micro%melthail))     deallocate (micro%melthail) 
    if (associated(micro%cld2rain))     deallocate (micro%cld2rain)
    if (associated(micro%rimecld))      deallocate (micro%rimecld)
    if (associated(micro%rimecldsnow))  deallocate (micro%rimecldsnow)
    if (associated(micro%rimecldaggr))  deallocate (micro%rimecldaggr)
    if (associated(micro%rimecldgrau))  deallocate (micro%rimecldgrau)
    if (associated(micro%rimecldhail))  deallocate (micro%rimecldhail)
    if (associated(micro%rain2ice))     deallocate (micro%rain2ice)
    if (associated(micro%rain2pr))      deallocate (micro%rain2pr) 
    if (associated(micro%rain2sn))      deallocate (micro%rain2sn) 
    if (associated(micro%rain2ag))      deallocate (micro%rain2ag) 
    if (associated(micro%rain2gr))      deallocate (micro%rain2gr) 
    if (associated(micro%rain2ha))      deallocate (micro%rain2ha)
    if (associated(micro%rain2ha_xtra)) deallocate (micro%rain2ha_xtra)
    if (associated(micro%ice2rain))     deallocate (micro%ice2rain)
    if (associated(micro%aggregate))    deallocate (micro%aggregate) 
    if (associated(micro%aggrselfpris)) deallocate (micro%aggrselfpris)
    if (associated(micro%aggrselfsnow)) deallocate (micro%aggrselfsnow)
    if (associated(micro%aggrprissnow)) deallocate (micro%aggrprissnow)
    if (associated(micro%latheatvap))   deallocate (micro%latheatvap)
    if (associated(micro%latheatfrz))   deallocate (micro%latheatfrz)

    !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES (totals)
    if (associated(micro%nuccldrt))      deallocate (micro%nuccldrt)
    if (associated(micro%nuccldct))      deallocate (micro%nuccldct)
    if (associated(micro%nucicert))      deallocate (micro%nucicert)
    if (associated(micro%nucicect))      deallocate (micro%nucicect)
    if (associated(micro%inuchomrt))     deallocate (micro%inuchomrt)     
    if (associated(micro%inuchomct))     deallocate (micro%inuchomct) 
    if (associated(micro%inuccontrt))    deallocate (micro%inuccontrt)     
    if (associated(micro%inuccontct))    deallocate (micro%inuccontct)
    if (associated(micro%inucifnrt))     deallocate (micro%inucifnrt)     
    if (associated(micro%inucifnct))     deallocate (micro%inucifnct) 
    if (associated(micro%inuchazrt))     deallocate (micro%inuchazrt)     
    if (associated(micro%inuchazct))     deallocate (micro%inuchazct)
    if (associated(micro%vapliqt))       deallocate (micro%vapliqt)
    if (associated(micro%vapicet))       deallocate (micro%vapicet)
    if (associated(micro%vapcldt))       deallocate (micro%vapcldt)
    if (associated(micro%vapraint))      deallocate (micro%vapraint)
    if (associated(micro%vapprist))      deallocate (micro%vapprist)
    if (associated(micro%vapsnowt))      deallocate (micro%vapsnowt)
    if (associated(micro%vapaggrt))      deallocate (micro%vapaggrt)
    if (associated(micro%vapgraut))      deallocate (micro%vapgraut)
    if (associated(micro%vaphailt))      deallocate (micro%vaphailt)
    if (associated(micro%vapdrizt))      deallocate (micro%vapdrizt)      
    if (associated(micro%melticet))      deallocate (micro%melticet) 
    if (associated(micro%meltprist))     deallocate (micro%meltprist) 
    if (associated(micro%meltsnowt))     deallocate (micro%meltsnowt) 
    if (associated(micro%meltaggrt))     deallocate (micro%meltaggrt)
    if (associated(micro%meltgraut))     deallocate (micro%meltgraut) 
    if (associated(micro%melthailt))     deallocate (micro%melthailt) 
    if (associated(micro%cld2raint))     deallocate (micro%cld2raint)
    if (associated(micro%rimecldt))      deallocate (micro%rimecldt)
    if (associated(micro%rimecldsnowt))  deallocate (micro%rimecldsnowt)
    if (associated(micro%rimecldaggrt))  deallocate (micro%rimecldaggrt)
    if (associated(micro%rimecldgraut))  deallocate (micro%rimecldgraut)
    if (associated(micro%rimecldhailt))  deallocate (micro%rimecldhailt)
    if (associated(micro%rain2icet))     deallocate (micro%rain2icet)
    if (associated(micro%rain2prt))      deallocate (micro%rain2prt) 
    if (associated(micro%rain2snt))      deallocate (micro%rain2snt) 
    if (associated(micro%rain2agt))      deallocate (micro%rain2agt) 
    if (associated(micro%rain2grt))      deallocate (micro%rain2grt) 
    if (associated(micro%rain2hat))      deallocate (micro%rain2hat)
    if (associated(micro%rain2ha_xtrat)) deallocate (micro%rain2ha_xtrat)
    if (associated(micro%ice2raint))     deallocate (micro%ice2raint)
    if (associated(micro%aggregatet))    deallocate (micro%aggregatet) 
    if (associated(micro%aggrselfprist)) deallocate (micro%aggrselfprist)
    if (associated(micro%aggrselfsnowt)) deallocate (micro%aggrselfsnowt)
    if (associated(micro%aggrprissnowt)) deallocate (micro%aggrprissnowt)
    if (associated(micro%latheatvapt))   deallocate (micro%latheatvapt)
    if (associated(micro%latheatfrzt))   deallocate (micro%latheatfrzt)

   if (associated(micro%rei))     deallocate (micro%rei)
   if (associated(micro%rel))     deallocate (micro%rel)
   if (associated(micro%cldfr))   deallocate (micro%cldfr)
   
   return
   end subroutine dealloc_micro


subroutine filltab_micro(micro,microm,imean,n1,n2,n3,ng)

   use var_tables, only: InsertVTab

   implicit none
   include "i8.h"
   type (micro_vars) :: micro,microm
   integer, intent(in) :: imean,n1,n2,n3,ng
   integer(kind=i8) :: npts
   real, pointer :: var,varm

! Fill pointers to arrays into variable tables

   npts=n1*n2*n3
   if (associated(micro%rcp))   &
      call InsertVTab (micro%rcp,microm%rcp  &
                 ,ng, npts, imean,  &
                 'RCP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%rdp))   &
      call InsertVTab (micro%rdp,microm%rdp  &
                 ,ng, npts, imean,  &
                 'RDP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%rrp))   &
      call InsertVTab (micro%rrp,microm%rrp  &
                 ,ng, npts, imean,  &
                 'RRP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%rpp))   &
      call InsertVTab (micro%rpp,microm%rpp  &
                 ,ng, npts, imean,  &
                 'RPP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%rsp))   &
      call InsertVTab (micro%rsp,microm%rsp  &
                 ,ng, npts, imean,  &
                 'RSP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%rap))   &
      call InsertVTab (micro%rap,microm%rap  &
                 ,ng, npts, imean,  &
                 'RAP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%rgp))   &
      call InsertVTab (micro%rgp,microm%rgp  &
                 ,ng, npts, imean,  &
                 'RGP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%rhp))   &
      call InsertVTab (micro%rhp,microm%rhp  &
                 ,ng, npts, imean,  &
                 'RHP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%ccp))   &
      call InsertVTab (micro%ccp,microm%ccp  &
                 ,ng, npts, imean,  &
                 'CCP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%cdp))   &
      call InsertVTab (micro%cdp,microm%cdp  &
                 ,ng, npts, imean,  &
                 'CDP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%crp))   &
      call InsertVTab (micro%crp,microm%crp  &
                 ,ng, npts, imean,  &
                 'CRP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%cpp))   &
      call InsertVTab (micro%cpp,microm%cpp  &
                 ,ng, npts, imean,  &
                 'CPP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%csp))   &
      call InsertVTab (micro%csp,microm%csp  &
                 ,ng, npts, imean,  &
                 'CSP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%cap))   &
      call InsertVTab (micro%cap,microm%cap  &
                 ,ng, npts, imean,  &
                 'CAP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%cgp))   &
      call InsertVTab (micro%cgp,microm%cgp  &
                 ,ng, npts, imean,  &
                 'CGP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%chp))   &
      call InsertVTab (micro%chp,microm%chp  &
                 ,ng, npts, imean,  &
                 'CHP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%cccnp)) &
      call InsertVTab (micro%cccnp,microm%cccnp  &
                 ,ng, npts, imean,  &
                 'CCCNP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%gccnp)) &
      call InsertVTab (micro%gccnp,microm%gccnp  &
                 ,ng, npts, imean,  &
                 'GCCNP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%cifnp)) &
      call InsertVTab (micro%cifnp,microm%cifnp  &
                 ,ng, npts, imean,  &
                 'CIFNP :3:hist:anal:mpti:mpt3:mpt1')

   if (associated(micro%q2))   &
      call InsertVTab (micro%q2,microm%q2  &
                 ,ng, npts, imean,  &
                 'Q2 :3:hist:anal:mpti:mpt3')
   if (associated(micro%q6)) &
      call InsertVTab (micro%q6,microm%q6  &
                 ,ng, npts, imean,  &
                 'Q6 :3:hist:anal:mpti:mpt3')
   if (associated(micro%q7)) &
      call InsertVTab (micro%q7,microm%q7  &
                 ,ng, npts, imean,  &
                 'Q7 :3:hist:anal:mpti:mpt3')

   if (associated(micro%cccmp)) &
      call InsertVTab (micro%cccmp,microm%cccmp  &
                 ,ng, npts, imean,  &
                 'CCCMP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%gccmp)) &
      call InsertVTab (micro%gccmp,microm%gccmp  &
                 ,ng, npts, imean,  &
                 'GCCMP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%cnm1p)) &
      call InsertVTab (micro%cnm1p,microm%cnm1p  &
                 ,ng, npts, imean,  &
                 'CNM1P :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%cnm2p)) &
      call InsertVTab (micro%cnm2p,microm%cnm2p  &
                 ,ng, npts, imean,  &
                 'CNM2P :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%cnm3p)) &
      call InsertVTab (micro%cnm3p,microm%cnm3p  &
                 ,ng, npts, imean,  &
                 'CNM3P :3:hist:anal:mpti:mpt3:mpt1')
  if (associated(micro%cnm8p)) &
      call InsertVTab (micro%cnm8p,microm%cnm8p  &
                 ,ng, npts, imean,  &
                 'CNM8P :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%md1np)) &
      call InsertVTab (micro%md1np,microm%md1np  &
                 ,ng, npts, imean,  &
                 'MD1NP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%md2np)) &
      call InsertVTab (micro%md2np,microm%md2np  &
                 ,ng, npts, imean,  &
                 'MD2NP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%salt_filmp)) &
      call InsertVTab (micro%salt_filmp,microm%salt_filmp  &
                 ,ng, npts, imean,  &
                 'SALT_FILMP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%salt_jetp)) &
      call InsertVTab (micro%salt_jetp,microm%salt_jetp  &
                 ,ng, npts, imean,  &
                 'SALT_JETP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%salt_spmp)) &
      call InsertVTab (micro%salt_spmp,microm%salt_spmp  &
                 ,ng, npts, imean,  &
                 'SALT_SPMP :3:hist:anal:mpti:mpt3:mpt1')
   if (associated(micro%rei))   &
      call InsertVTab (micro%rei,microm%rei  &
                 ,ng, npts, imean,  &
                 'REI :3:hist:anal:mpti:mpt3')
   if (associated(micro%rel))   &
      call InsertVTab (micro%rel,microm%rel  &
                 ,ng, npts, imean,  &
                 'REL :3:hist:anal:mpti:mpt3')
   if (associated(micro%cldfr))   &
      call InsertVTab (micro%cldfr,microm%cldfr  &
                 ,ng, npts, imean,  &
                 'CLDFR :3:hist:anal:mpti:mpt3')
  
  
   !VERTICAL PRECIPITATION RATES
   if (associated(micro%pcpvr)) &
      call InsertVTab (micro%pcpvr,microm%pcpvr  &
                 ,ng, npts, imean,  &
                 'PCPVR :3:hist:anal:mpti:mpt3')
   if (associated(micro%pcpvp)) &
      call InsertVTab (micro%pcpvp,microm%pcpvp  &
                 ,ng, npts, imean,  &
                 'PCPVP :3:hist:anal:mpti:mpt3')
   if (associated(micro%pcpvs)) &
      call InsertVTab (micro%pcpvs,microm%pcpvs  &
                 ,ng, npts, imean,  &
                 'PCPVS :3:hist:anal:mpti:mpt3')
   if (associated(micro%pcpva)) &
      call InsertVTab (micro%pcpva,microm%pcpva  &
                 ,ng, npts, imean,  &
                 'PCPVA :3:hist:anal:mpti:mpt3')
   if (associated(micro%pcpvg)) &
      call InsertVTab (micro%pcpvg,microm%pcpvg  &
                 ,ng, npts, imean,  &
                 'PCPVG :3:hist:anal:mpti:mpt3')
   if (associated(micro%pcpvh)) &
      call InsertVTab (micro%pcpvh,microm%pcpvh  &
                 ,ng, npts, imean,  &
                 'PCPVH :3:hist:anal:mpti:mpt3')
   if (associated(micro%pcpvd)) &
      call InsertVTab (micro%pcpvd,microm%pcpvd  &
                 ,ng, npts, imean,  &
                 'PCPVD :3:hist:anal:mpti:mpt3')

    !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES (instantaneous)
    if (associated(micro%nuccldr)) &
      call InsertVTab (micro%nuccldr,microm%nuccldr  &
                 ,ng, npts, imean,  &
                 'NUCCLDR :3:hist:anal:mpt3')
    if (associated(micro%nuccldc)) &
      call InsertVTab (micro%nuccldc,microm%nuccldc  &
                 ,ng, npts, imean,  &
                 'NUCCLDC :3:hist:anal:mpt3')

    if (associated(micro%nucicer)) &
      call InsertVTab (micro%nucicer,microm%nucicer  &
                 ,ng, npts, imean,  &
                 'NUCICER :3:hist:anal:mpt3')
    if (associated(micro%nucicec)) &
      call InsertVTab (micro%nucicec,microm%nucicec  &
                 ,ng, npts, imean,  &
                 'NUCICEC :3:hist:anal:mpt3')
    if (associated(micro%inuchomr)) &
      call InsertVTab (micro%inuchomr,microm%inuchomr  &
                 ,ng, npts, imean,  &
                 'INUCHOMR :3:hist:anal:mpt3')
    if (associated(micro%inuchomc)) &
      call InsertVTab (micro%inuchomc,microm%inuchomc  &
                 ,ng, npts, imean,  &
                 'INUCHOMC :3:hist:anal:mpt3')
    if (associated(micro%inuccontr)) &
      call InsertVTab (micro%inuccontr,microm%inuccontr  &
                 ,ng, npts, imean,  &
                 'INUCCONTR :3:hist:anal:mpt3')
    if (associated(micro%inuccontc)) &
      call InsertVTab (micro%inuccontc,microm%inuccontc  &
                 ,ng, npts, imean,  &
                 'INUCCONTC :3:hist:anal:mpt3')
   if (associated(micro%inucifnr)) &
      call InsertVTab (micro%inucifnr,microm%inucifnr  &
                 ,ng, npts, imean,  &
                 'INUCIFNR :3:hist:anal:mpt3')
   if (associated(micro%inucifnc)) &
      call InsertVTab (micro%inucifnc,microm%inucifnc  &
                 ,ng, npts, imean,  &
                 'INUCIFNC :3:hist:anal:mpt3')
   if (associated(micro%inuchazr)) &
      call InsertVTab (micro%inuchazr,microm%inuchazr  &
                 ,ng, npts, imean,  &
                 'INUCHAZR :3:hist:anal:mpt3')
   if (associated(micro%inuchazc)) &
      call InsertVTab (micro%inuchazc,microm%inuchazc  &
                 ,ng, npts, imean,  &
                 'INUCHAZC :3:hist:anal:mpt3')

    if (associated(micro%vapliq)) &
      call InsertVTab (micro%vapliq,microm%vapliq  &
                 ,ng, npts, imean,  &
                 'VAPLIQ :3:hist:anal:mpt3')
    if (associated(micro%vapice)) &
      call InsertVTab (micro%vapice,microm%vapice  &
                 ,ng, npts, imean,  &
                 'VAPICE :3:hist:anal:mpt3')
    if (associated(micro%vapcld)) &
      call InsertVTab (micro%vapcld,microm%vapcld  &
                 ,ng, npts, imean,  &
                 'VAPCLD :3:hist:anal:mpt3')
    if (associated(micro%vaprain)) &
      call InsertVTab (micro%vaprain,microm%vaprain  &
                 ,ng, npts, imean,  &
                 'VAPRAIN :3:hist:anal:mpt3')
    if (associated(micro%vappris)) &
      call InsertVTab (micro%vappris,microm%vappris  &
                 ,ng, npts, imean,  &
                 'VAPPRIS :3:hist:anal:mpt3')
    if (associated(micro%vapsnow)) &
      call InsertVTab (micro%vapsnow,microm%vapsnow  &
                 ,ng, npts, imean,  &
                 'VAPSNOW :3:hist:anal:mpt3')
    if (associated(micro%vapaggr)) &
      call InsertVTab (micro%vapaggr,microm%vapaggr  &
                 ,ng, npts, imean,  &
                 'VAPAGGR :3:hist:anal:mpt3')
    if (associated(micro%vapgrau)) &
      call InsertVTab (micro%vapgrau,microm%vapgrau  &
                 ,ng, npts, imean,  &
                 'VAPGRAU :3:hist:anal:mpt3')
    if (associated(micro%vaphail)) &
      call InsertVTab (micro%vaphail,microm%vaphail  &
                 ,ng, npts, imean,  &
                 'VAPHAIL :3:hist:anal:mpt3')
    if (associated(micro%vapdriz)) &
      call InsertVTab (micro%vapdriz,microm%vapdriz  &
                 ,ng, npts, imean,  &
                 'VAPDRIZ :3:hist:anal:mpt3')
 
   if (associated(micro%meltice)) &
      call InsertVTab (micro%meltice,microm%meltice  &
                 ,ng, npts, imean,  &
                 'MELTICE :3:hist:anal:mpt3')
   if (associated(micro%meltpris)) &
      call InsertVTab (micro%meltpris,microm%meltpris  &
                 ,ng, npts, imean,  &
                 'MELTPRIS :3:hist:anal:mpt3')
    if (associated(micro%meltsnow)) &
      call InsertVTab (micro%meltsnow,microm%meltsnow  &
                 ,ng, npts, imean,  &
                 'MELTSNOW :3:hist:anal:mpt3')
    if (associated(micro%meltaggr)) &
      call InsertVTab (micro%meltaggr,microm%meltaggr  &
                 ,ng, npts, imean,  &
                 'MELTAGGR :3:hist:anal:mpt3')
    if (associated(micro%meltgrau)) &
      call InsertVTab (micro%meltgrau,microm%meltgrau  &
                 ,ng, npts, imean,  &
                 'MELTGRAU :3:hist:anal:mpt3')
    if (associated(micro%melthail)) &
      call InsertVTab (micro%melthail,microm%melthail  &
                 ,ng, npts, imean,  &
                 'MELTHAIL :3:hist:anal:mpt3')

    if (associated(micro%cld2rain)) &
      call InsertVTab (micro%cld2rain,microm%cld2rain  &
                 ,ng, npts, imean,  &
                 'CLD2RAIN :3:hist:anal:mpt3')
    if (associated(micro%rimecld)) &
      call InsertVTab (micro%rimecld,microm%rimecld  &
                 ,ng, npts, imean,  &
                 'RIMECLD :3:hist:anal:mpt3')
    if (associated(micro%rimecldsnow)) &
      call InsertVTab (micro%rimecldsnow,microm%rimecldsnow  &
                 ,ng, npts, imean,  &
                 'RIMECLDSNOW :3:hist:anal:mpt3')
    if (associated(micro%rimecldaggr)) &
      call InsertVTab (micro%rimecldaggr,microm%rimecldaggr  &
                 ,ng, npts, imean,  &
                 'RIMECLDAGGR :3:hist:anal:mpt3')
    if (associated(micro%rimecldgrau)) &
      call InsertVTab (micro%rimecldgrau,microm%rimecldgrau  &
                 ,ng, npts, imean,  &
                 'RIMECLDGRAU :3:hist:anal:mpt3')
    if (associated(micro%rimecldhail)) &
      call InsertVTab (micro%rimecldhail,microm%rimecldhail  &
                 ,ng, npts, imean,  &
                 'RIMECLDHAIL :3:hist:anal:mpt3')

    if (associated(micro%rain2ice)) &
      call InsertVTab (micro%rain2ice,microm%rain2ice  &
                 ,ng, npts, imean,  &
                 'RAIN2ICE :3:hist:anal:mpt3')
    if (associated(micro%rain2pr)) &
      call InsertVTab (micro%rain2pr,microm%rain2pr  &
                 ,ng, npts, imean,  &
                 'RAIN2PR :3:hist:anal:mpt3')
    if (associated(micro%rain2sn)) &
      call InsertVTab (micro%rain2sn,microm%rain2sn  &
                 ,ng, npts, imean,  &
                 'RAIN2SN :3:hist:anal:mpt3')
    if (associated(micro%rain2ag)) &
      call InsertVTab (micro%rain2ag,microm%rain2ag  &
                 ,ng, npts, imean,  &
                 'RAIN2AG :3:hist:anal:mpt3')
    if (associated(micro%rain2gr)) &
      call InsertVTab (micro%rain2gr,microm%rain2gr  &
                 ,ng, npts, imean,  &
                 'RAIN2GR :3:hist:anal:mpt3')
    if (associated(micro%rain2ha)) &
      call InsertVTab (micro%rain2ha,microm%rain2ha  &
                 ,ng, npts, imean,  &
                 'RAIN2HA :3:hist:anal:mpt3')
    if (associated(micro%rain2ha_xtra)) &
      call InsertVTab (micro%rain2ha_xtra,microm%rain2ha_xtra  &
                 ,ng, npts, imean,  &
                 'RAIN2HA_XTRA :3:hist:anal:mpt3')
    if (associated(micro%ice2rain)) &
      call InsertVTab (micro%ice2rain,microm%ice2rain  &
                 ,ng, npts, imean,  &
                 'ICE2RAIN :3:hist:anal:mpt3')

   if (associated(micro%aggregate)) &
      call InsertVTab (micro%aggregate,microm%aggregate  &
                 ,ng, npts, imean,  &
                 'AGGREGATE :3:hist:anal:mpt3')
   if (associated(micro%aggrselfpris)) &
      call InsertVTab (micro%aggrselfpris,microm%aggrselfpris  &
                 ,ng, npts, imean,  &
                 'AGGRSELFPRIS :3:hist:anal:mpt3')
   if (associated(micro%aggrselfsnow)) &
      call InsertVTab (micro%aggrselfsnow,microm%aggrselfsnow  &
                 ,ng, npts, imean,  &
                 'AGGRSELFSNOW :3:hist:anal:mpt3')
   if (associated(micro%aggrprissnow)) &
      call InsertVTab (micro%aggrprissnow,microm%aggrprissnow  &
                 ,ng, npts, imean,  &
                 'AGGRPRISSNOW :3:hist:anal:mpt3')

   if (associated(micro%latheatvap)) &
      call InsertVTab (micro%latheatvap,microm%latheatvap  &
                 ,ng, npts, imean,  &
                 'LATHEATVAP :3:hist:anal:mpt3')
   if (associated(micro%latheatfrz)) &
      call InsertVTab (micro%latheatfrz,microm%latheatfrz  &
                 ,ng, npts, imean,  &
                 'LATHEATFRZ :3:hist:anal:mpt3')
   !END MICRO BUDGET PROCESSES (instantaneous)

    !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES (totals)
    if (associated(micro%nuccldrt)) &
      call InsertVTab (micro%nuccldrt,microm%nuccldrt  &
                 ,ng, npts, imean,  &
                 'NUCCLDRT :3:hist:anal:mpti:mpt3')
    if (associated(micro%nuccldct)) &
      call InsertVTab (micro%nuccldct,microm%nuccldct  &
                 ,ng, npts, imean,  &
                 'NUCCLDCT :3:hist:anal:mpti:mpt3')

    if (associated(micro%nucicert)) &
      call InsertVTab (micro%nucicert,microm%nucicert  &
                 ,ng, npts, imean,  &
                 'NUCICERT :3:hist:anal:mpti:mpt3')
    if (associated(micro%nucicect)) &
      call InsertVTab (micro%nucicect,microm%nucicect  &
                 ,ng, npts, imean,  &
                 'NUCICECT :3:hist:anal:mpti:mpt3')
    if (associated(micro%inuchomrt)) &
      call InsertVTab (micro%inuchomrt,microm%inuchomrt  &
                 ,ng, npts, imean,  &
                 'INUCHOMRT :3:hist:anal:mpti:mpt3')
    if (associated(micro%inuchomct)) &
      call InsertVTab (micro%inuchomct,microm%inuchomct  &
                 ,ng, npts, imean,  &
                 'INUCHOMCT :3:hist:anal:mpti:mpt3')
    if (associated(micro%inuccontrt)) &
      call InsertVTab (micro%inuccontrt,microm%inuccontrt  &
                 ,ng, npts, imean,  &
                 'INUCCONTRT :3:hist:anal:mpti:mpt3')
    if (associated(micro%inuccontct)) &
      call InsertVTab (micro%inuccontct,microm%inuccontct  &
                 ,ng, npts, imean,  &
                 'INUCCONTCT :3:hist:anal:mpti:mpt3')
   if (associated(micro%inucifnrt)) &
      call InsertVTab (micro%inucifnrt,microm%inucifnrt  &
                 ,ng, npts, imean,  &
                 'INUCIFNRT :3:hist:anal:mpti:mpt3')
   if (associated(micro%inucifnct)) &
      call InsertVTab (micro%inucifnct,microm%inucifnct  &
                 ,ng, npts, imean,  &
                 'INUCIFNCT :3:hist:anal:mpti:mpt3')
   if (associated(micro%inuchazrt)) &
      call InsertVTab (micro%inuchazrt,microm%inuchazrt  &
                 ,ng, npts, imean,  &
                 'INUCHAZRT :3:hist:anal:mpti:mpt3')
   if (associated(micro%inuchazct)) &
      call InsertVTab (micro%inuchazct,microm%inuchazct  &
                 ,ng, npts, imean,  &
                 'INUCHAZCT :3:hist:anal:mpti:mpt3')

    if (associated(micro%vapliqt)) &
      call InsertVTab (micro%vapliqt,microm%vapliqt  &
                 ,ng, npts, imean,  &
                 'VAPLIQT :3:hist:anal:mpti:mpt3')
    if (associated(micro%vapicet)) &
      call InsertVTab (micro%vapicet,microm%vapicet  &
                 ,ng, npts, imean,  &
                 'VAPICET :3:hist:anal:mpti:mpt3')
    if (associated(micro%vapcldt)) &
      call InsertVTab (micro%vapcldt,microm%vapcldt  &
                 ,ng, npts, imean,  &
                 'VAPCLDT :3:hist:anal:mpti:mpt3')
    if (associated(micro%vapraint)) &
      call InsertVTab (micro%vapraint,microm%vapraint  &
                 ,ng, npts, imean,  &
                 'VAPRAINT :3:hist:anal:mpti:mpt3')
    if (associated(micro%vapprist)) &
      call InsertVTab (micro%vapprist,microm%vapprist  &
                 ,ng, npts, imean,  &
                 'VAPPRIST :3:hist:anal:mpti:mpt3')
    if (associated(micro%vapsnowt)) &
      call InsertVTab (micro%vapsnowt,microm%vapsnowt  &
                 ,ng, npts, imean,  &
                 'VAPSNOWT :3:hist:anal:mpti:mpt3')
    if (associated(micro%vapaggrt)) &
      call InsertVTab (micro%vapaggrt,microm%vapaggrt  &
                 ,ng, npts, imean,  &
                 'VAPAGGRT :3:hist:anal:mpti:mpt3')
    if (associated(micro%vapgraut)) &
      call InsertVTab (micro%vapgraut,microm%vapgraut  &
                 ,ng, npts, imean,  &
                 'VAPGRAUT :3:hist:anal:mpti:mpt3')
    if (associated(micro%vaphailt)) &
      call InsertVTab (micro%vaphailt,microm%vaphailt  &
                 ,ng, npts, imean,  &
                 'VAPHAILT :3:hist:anal:mpti:mpt3')
    if (associated(micro%vapdrizt)) &
      call InsertVTab (micro%vapdrizt,microm%vapdrizt  &
                 ,ng, npts, imean,  &
                 'VAPDRIZT :3:hist:anal:mpti:mpt3')
 
   if (associated(micro%melticet)) &
      call InsertVTab (micro%melticet,microm%melticet  &
                 ,ng, npts, imean,  &
                 'MELTICET :3:hist:anal:mpti:mpt3')
   if (associated(micro%meltprist)) &
      call InsertVTab (micro%meltprist,microm%meltprist  &
                 ,ng, npts, imean,  &
                 'MELTPRIST :3:hist:anal:mpti:mpt3')
    if (associated(micro%meltsnowt)) &
      call InsertVTab (micro%meltsnowt,microm%meltsnowt  &
                 ,ng, npts, imean,  &
                 'MELTSNOWT :3:hist:anal:mpti:mpt3')
    if (associated(micro%meltaggrt)) &
      call InsertVTab (micro%meltaggrt,microm%meltaggrt  &
                 ,ng, npts, imean,  &
                 'MELTAGGRT :3:hist:anal:mpti:mpt3')
    if (associated(micro%meltgraut)) &
      call InsertVTab (micro%meltgraut,microm%meltgraut  &
                 ,ng, npts, imean,  &
                 'MELTGRAUT :3:hist:anal:mpti:mpt3')
    if (associated(micro%melthailt)) &
      call InsertVTab (micro%melthailt,microm%melthailt  &
                 ,ng, npts, imean,  &
                 'MELTHAILT :3:hist:anal:mpti:mpt3')

    if (associated(micro%cld2raint)) &
      call InsertVTab (micro%cld2raint,microm%cld2raint  &
                 ,ng, npts, imean,  &
                 'CLD2RAINT :3:hist:anal:mpti:mpt3')
    if (associated(micro%rimecldt)) &
      call InsertVTab (micro%rimecldt,microm%rimecldt  &
                 ,ng, npts, imean,  &
                 'RIMECLDT :3:hist:anal:mpti:mpt3')
    if (associated(micro%rimecldsnowt)) &
      call InsertVTab (micro%rimecldsnowt,microm%rimecldsnowt  &
                 ,ng, npts, imean,  &
                 'RIMECLDSNOWT :3:hist:anal:mpti:mpt3')
    if (associated(micro%rimecldaggrt)) &
      call InsertVTab (micro%rimecldaggrt,microm%rimecldaggrt  &
                 ,ng, npts, imean,  &
                 'RIMECLDAGGRT :3:hist:anal:mpti:mpt3')
    if (associated(micro%rimecldgraut)) &
      call InsertVTab (micro%rimecldgraut,microm%rimecldgraut  &
                 ,ng, npts, imean,  &
                 'RIMECLDGRAUT :3:hist:anal:mpti:mpt3')
    if (associated(micro%rimecldhailt)) &
      call InsertVTab (micro%rimecldhailt,microm%rimecldhailt  &
                 ,ng, npts, imean,  &
                 'RIMECLDHAILT :3:hist:anal:mpti:mpt3')

    if (associated(micro%rain2icet)) &
      call InsertVTab (micro%rain2icet,microm%rain2icet  &
                 ,ng, npts, imean,  &
                 'RAIN2ICET :3:hist:anal:mpti:mpt3')
    if (associated(micro%rain2prt)) &
      call InsertVTab (micro%rain2prt,microm%rain2prt  &
                 ,ng, npts, imean,  &
                 'RAIN2PRT :3:hist:anal:mpti:mpt3')
    if (associated(micro%rain2snt)) &
      call InsertVTab (micro%rain2snt,microm%rain2snt  &
                 ,ng, npts, imean,  &
                 'RAIN2SNT :3:hist:anal:mpti:mpt3')
    if (associated(micro%rain2agt)) &
      call InsertVTab (micro%rain2agt,microm%rain2agt  &
                 ,ng, npts, imean,  &
                 'RAIN2AGT :3:hist:anal:mpti:mpt3')
    if (associated(micro%rain2grt)) &
      call InsertVTab (micro%rain2grt,microm%rain2grt  &
                 ,ng, npts, imean,  &
                 'RAIN2GRT :3:hist:anal:mpti:mpt3')
    if (associated(micro%rain2hat)) &
      call InsertVTab (micro%rain2hat,microm%rain2hat  &
                 ,ng, npts, imean,  &
                 'RAIN2HAT :3:hist:anal:mpti:mpt3')
    if (associated(micro%rain2ha_xtrat)) &
      call InsertVTab (micro%rain2ha_xtrat,microm%rain2ha_xtrat  &
                 ,ng, npts, imean,  &
                 'RAIN2HA_XTRAT :3:hist:anal:mpti:mpt3')
    if (associated(micro%ice2raint)) &
      call InsertVTab (micro%ice2raint,microm%ice2raint  &
                 ,ng, npts, imean,  &
                 'ICE2RAINT :3:hist:anal:mpti:mpt3')

   if (associated(micro%aggregatet)) &
      call InsertVTab (micro%aggregatet,microm%aggregatet  &
                 ,ng, npts, imean,  &
                 'AGGREGATET :3:hist:anal:mpti:mpt3')
   if (associated(micro%aggrselfprist)) &
      call InsertVTab (micro%aggrselfprist,microm%aggrselfprist  &
                 ,ng, npts, imean,  &
                 'AGGRSELFPRIST :3:hist:anal:mpti:mpt3')
   if (associated(micro%aggrselfsnowt)) &
      call InsertVTab (micro%aggrselfsnowt,microm%aggrselfsnowt  &
                 ,ng, npts, imean,  &
                 'AGGRSELFSNOWT :3:hist:anal:mpti:mpt3')
   if (associated(micro%aggrprissnowt)) &
      call InsertVTab (micro%aggrprissnowt,microm%aggrprissnowt  &
                 ,ng, npts, imean,  &
                 'AGGRPRISSNOWT :3:hist:anal:mpti:mpt3')

   if (associated(micro%latheatvapt)) &
      call InsertVTab (micro%latheatvapt,microm%latheatvapt  &
                 ,ng, npts, imean,  &
                 'LATHEATVAPT :3:hist:anal:mpti:mpt3')
   if (associated(micro%latheatfrzt)) &
      call InsertVTab (micro%latheatfrzt,microm%latheatfrzt  &
                 ,ng, npts, imean,  &
                 'LATHEATFRZT :3:hist:anal:mpti:mpt3')
   !END MICRO BUDGET PROCECCES (totals)
                 
   npts=n2*n3
   if (associated(micro%accpr)) &
      call InsertVTab (micro%accpr,microm%accpr  &
                 ,ng, npts, imean,  &
                 'ACCPR :2:hist:anal:mpti:mpt3')
   if (associated(micro%accpp)) &
      call InsertVTab (micro%accpp,microm%accpp  &
                 ,ng, npts, imean,  &
                 'ACCPP :2:hist:anal:mpti:mpt3')
   if (associated(micro%accps)) &
      call InsertVTab (micro%accps,microm%accps  &
                 ,ng, npts, imean,  &
                 'ACCPS :2:hist:anal:mpti:mpt3')
   if (associated(micro%accpa)) &
      call InsertVTab (micro%accpa,microm%accpa  &
                 ,ng, npts, imean,  &
                 'ACCPA :2:hist:anal:mpti:mpt3')
   if (associated(micro%accpg)) &
      call InsertVTab (micro%accpg,microm%accpg  &
                 ,ng, npts, imean,  &
                 'ACCPG :2:hist:anal:mpti:mpt3')
   if (associated(micro%accph)) &
      call InsertVTab (micro%accph,microm%accph  &
                 ,ng, npts, imean,  &
                 'ACCPH :2:hist:anal:mpti:mpt3')
   if (associated(micro%accpd)) &
      call InsertVTab (micro%accpd,microm%accpd  &
                 ,ng, npts, imean,  &
                 'ACCPD :2:hist:anal:mpti:mpt3')
   if (associated(micro%pcprr)) &
      call InsertVTab (micro%pcprr,microm%pcprr  &
                 ,ng, npts, imean,  &
                 'PCPRR :2:hist:anal:mpt3')
   if (associated(micro%pcprp)) &
      call InsertVTab (micro%pcprp,microm%pcprp  &
                 ,ng, npts, imean,  &
                 'PCPRP :2:hist:anal:mpt3')
   if (associated(micro%pcprs)) &
      call InsertVTab (micro%pcprs,microm%pcprs  &
                 ,ng, npts, imean,  &
                 'PCPRS :2:hist:anal:mpt3')
   if (associated(micro%pcpra)) &
      call InsertVTab (micro%pcpra,microm%pcpra  &
                 ,ng, npts, imean,  &
                 'PCPRA :2:hist:anal:mpt3')
   if (associated(micro%pcprg)) &
      call InsertVTab (micro%pcprg,microm%pcprg  &
                 ,ng, npts, imean,  &
                 'PCPRG :2:hist:anal:mpt3')
   if (associated(micro%pcprh)) &
      call InsertVTab (micro%pcprh,microm%pcprh  &
                 ,ng, npts, imean,  &
                 'PCPRH :2:hist:anal:mpt3')
   if (associated(micro%pcprd)) &
      call InsertVTab (micro%pcprd,microm%pcprd  &
                 ,ng, npts, imean,  &
                 'PCPRD :2:hist:anal:mpt3')
   if (associated(micro%pcpg)) &
      call InsertVTab (micro%pcpg,microm%pcpg  &
                 ,ng, npts, imean,  &
                 'PCPG :2:hist:mpti:mpt3')
   if (associated(micro%qpcpg)) &
      call InsertVTab (micro%qpcpg,microm%qpcpg  &
                 ,ng, npts, imean,  &
                 'QPCPG :2:hist:mpti:mpt3')
   if (associated(micro%dpcpg)) &
      call InsertVTab (micro%dpcpg,microm%dpcpg  &
                 ,ng, npts, imean,  &
                 'DPCPG :2:hist:mpti:mpt3')

   return
   end subroutine filltab_micro

End Module mem_micro
