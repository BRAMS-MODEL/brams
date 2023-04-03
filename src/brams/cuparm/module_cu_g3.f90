!---------------------------------------------------------------------!
! Grell cumulus version 3d - implemented in CCATT-BRAMS/BRAMS feb 2012!
! Saulo R. Freitas: 1st version                                       !
! Rafael Mello: included parallelization for spread/smoothh arrays    !
!---------------------------------------------------------------------!

MODULE module_cu_g3

CONTAINS

!-------------------------------------------------------------
   SUBROUTINE G3DRV(mynum,i0,j0,time                            &
              ,DT                                               &
	      ,DX                            			&
              ,autoconv                                         & !
              ,aeroevap                                         & !
              ,rho						&
!	      ,RAINCV						&
	      ,PRATEC                        			&
              ,U                                                &
	      ,V						&
	      ,theta						&
              ,thetail                                          &
              ,pp                                               &
              ,pi0   	                                        &
	      ,W						&
	      ,rv						&
              ,rtgt                                             &
              ,pt                                               &
	      ,XLV						&
	      ,CP						&
	      ,G						&
	      ,r_v                           			&
	      ,p00                           			&
	      ,cpor                           			&
              ,APR_GR						&
	      ,APR_W						&
	      ,APR_MC						&
	      ,APR_ST						&
	      ,APR_AS             				&
!             ,APR_CAPMA					&
!	      ,APR_CAPME					&
!	      ,APR_CAPMI          				&
              ,MASS_FLUX					&
	      ,xmb_shallow        				&
!
              ,weight_GR  					&
              ,weight_W   					&
              ,weight_MC  					&
              ,weight_ST  					&
              ,weight_AS  					&
!
              ,training   					&
!
!	      ,XF_ENS						&
!	      ,PR_ENS						&
	      ,HT						&
	      ,patch_area					&
              ,npat                                             &
	      ,gsw						&
!	      ,edt_out   					&
!             ,GDC						&
!	      ,GDC2						&
!	      ,kpbl						&
!	      ,k22_shallow					&
!	      ,kbcon_shallow      				&
!             ,ktop_shallow					&
	      ,cugd_avedx					&
	      ,imomentum          				&
              ,ensdim,maxiens,maxens,maxens2,maxens3,ichoice    &
              ,ishallow_g3                                      &
	      ,ids,ide, jds,jde, kds,kde                        &
              ,ims,ime, jms,jme, kms,kme                        &
              ,ips,ipe, jps,jpe, kps,kpe                        &
              ,its,ite, jts,jte, kts,kte                        &
	      ,RTHCUTEN                     		        &
              ,RQVCUTEN                                         &
	      ,RQCCUTEN 				        &
!	      ,RQICUTEN                     		        &
              ,cugd_ttens					&
	      ,cugd_qvtens					&
!--
!             ,cugd_tten					&
!	      ,cugd_qvten					&
!	      ,cugd_qcten         				&
!--
! forcings - for deep/shallow
	      ,RTHFTEN					        &
              ,RQVFTEN					        &
	      ,rthblten                     		        &
              ,rqvblten 				        &
              ,level              				&
              ,rcp	          				&
              ,rrp	          				&
              ,rpp	          				&
              ,rsp	          				&
              ,rap	          				&
              ,rgp	          				&
              ,rhp	          				&
	      ,aot500					        &
!srf-tmp
!              ,extra2d &
!              ,extra3d ,extra3d2,extra3d3&
!srf-tmp
	      ,F_QV					        &
	      ,F_QC					        &
	      ,F_QR					        &
	      ,F_QI					        &
	      ,F_QS					       )!&
!#if ( WRF_DFI_RADAR == 1 )
!                 ! Optional CAP suppress option
!              ,do_capsuppress,cap_suppress_loc                  &
!#endif
!
!
!--------------------------------------------------------------------
!- variables :
!DT           : model timestep (s)                INTENT(in)
!itimestep    : not used
!DX	      : horizontal grid spacing	(m)	  INTENT(in)
!rho          : air density (kg/m^3)              INTENT(in)
!RAINCV       : accumulated rain over one timestep (kg/m^2 ) INTENT(out)
!PRATEC	      :	tendency of rain at that timestep (kg/m^2 s) INTENT(out)
!U            : wind x-dir (m/s)                  INTENT(in)
!V	      : wind y-dir (m/s)                  INTENT(in)
!t	      : air temperature (K)               iNTENT(in)
!W	      : vertical wind velocity (m/s)      iNTENT(in)
!q	      : water vapor mixing ratio (kg/kg)  INTENT(in)
!p	      : air pressure (hPa)                INTENT(in)
!pi 	      : (temperature air tendency to convert to project temp)
!dz8w	      : vertical grid spacing
!p8w	      : pressure at layer interfaces.
!XLV	      : cte
!CP	      : dry air specific heat (J/kg K)
!G	      : CTE GRAVITY (m/s^2)
!r_v  	      :  cte
!STEPCU       : not used
!htop	      : not used
!hbot	      :  not used
!CU_ACT_FLAG  :  not used
!warm_rain    :  not used
!APR_GR       : rainfall tendency for each closure (take out latter)
!APR_W	      : rainfall tendency for each closure (take out latter)
!APR_MC       : rainfall tendency for each closure (take out latter)
!APR_ST       : rainfall tendency for each closure (take out latter)
!APR_AS	      : rainfall tendency for each closure (take out latter)
!APR_CAPMA    : rainfall tendency for each closure (take out latter)
!APR_CAPME    : rainfall tendency for each closure (take out latter)
!APR_CAPMI    : rainfall tendency for each closure (take out latter)
!MASS_FLUX    : not used
!XF_ENS       : output
!PR_ENS       : output
!HT	      : terrain height
!XLAND        : land use (1 or 2)
!gsw	      : short wave radiation (only for change cap_max)
!edt_out      : output
!GDC	      : cloud water mixing ratio (output) - time average over the radiation
!                timestep
!GDC2	      : ice water mixing ratio (output) - time average over the radiation
!                timestep
!kpbl	      : the level of pbl height for shallow scheme
!k22_shallow  : output for tracer transport  for shallow convection
!kbcon_shallow: output for tracer transport  for shallow convection
!ktop_shallow : output for tracer transport  for shallow convection
!xmb_shallow  : output for tracer transport  for shallow convection
!
!     only for of the convective column
!cugd_tten    : tendencies for temp
!cugd_qvten   : tendencies for qv
!cugd_qcten   : tendencies for qc
!
!     only for the neighboors of the convective column (might include the central column
!cugd_ttens   : tendencies for temp due to subsidence/detrainment at top
!cugd_qvtens  : tendencies for qv due to subsidence/detrainment at top
!cugd_avedx   : 1 for the not-spreading - 3 for  spreading
!               use it only for grid spacing less than 10 km
!imomentum    : not used
!ensdim       :
!maxiens      :
!maxens       :
!maxens2      :
!maxens3      :
!ichoice      :
!ishallow_g3  : to turnon-off shallow convection
!,ids,ide, jds,jde, kds,kde  :   not in use
!,ims,ime, jms,jme, kms,kme  :	 not in use
!,ips,ipe, jps,jpe, kps,kpe  :	 not in use
!,its,ite, jts,jte, kts,kte  :	only this is in use
!periodic_x   : for global -wrf applications (use false for limited area domain)
!periodic_y   : for global -wrf applications (use false for limited area domain

!RQVCUTEN     : output tendencies for water vapor
!RQCCUTEN     : output tendencies for cloud liq
!RQICUTEN     : output tendencies for ice
!RTHCUTEN     : output tendencies for temp

!RQVFTEN      : input forcing for water vapor
!RTHFTEN      : input forcing for temp
!rqvblten     : forcing for only PBL (shallow) - water vapor
!rthblten     : forcing for only PBL (shallow) - temp

!F_QV	      : logical for existence of this variable - microphysics
!F_QC         : logical for existence of this variable - microphysics
!F_QR         : logical for existence of this variable - microphysics
!F_QI         : logical for existence of this variable - microphysics
!F_QS	      : logical for existence of this variable - microphysics
!-------------------------------------------------------------
!-------------------------------------------------------------
   IMPLICIT NONE
! beta is the massflux percentage of downdraft massflx that will reach the ground
      real, parameter :: beta=0.05
!-------------------------------------------------------------
   INTEGER,      INTENT(IN   ) ::                               &
                                  ids,ide, jds,jde, kds,kde,    &
                                  ims,ime, jms,jme, kms,kme,    &
                                  ips,ipe, jps,jpe, kps,kpe,    &
                                  its,ite, jts,jte, kts,kte,    &
				  mynum,i0,j0
   LOGICAL :: periodic_x=.false. ,periodic_y=.false.
               integer, parameter  :: ens4_spread = 3 ! max(3,cugd_avedx)
               integer, parameter  :: ens4=ens4_spread*ens4_spread

   integer, intent (in   )              ::                      &
                       ensdim,maxiens,maxens,maxens2,maxens3,ichoice,NPAT,level&
		       ,training

!   INTEGER,      INTENT(IN   ) :: STEPCU, ITIMESTEP,cugd_avedx, &
!                                  ishallow_g3,imomentum
   INTEGER,      INTENT(IN   ) ::  cugd_avedx,          &
                                  ishallow_g3,imomentum,&
				  autoconv,             & !
                                  aeroevap

   INTEGER ::  ITIMESTEP=0,STEPCU=0 !SRF

!   LOGICAL,      INTENT(IN   ) :: warm_rain !SRF

   REAL,         INTENT(IN   ) :: XLV, R_v
   REAL,         INTENT(IN   ) :: CP,G, cpor, p00

!srf
!  REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,    &
   REAL,  DIMENSION(kms:kme ,  ims:ime , jms:jme )         ,    &
          INTENT(IN   ) ::                                      &
                                                          U,    &
                                                          V,    &
                                                          W,    &
                                                          rv,   &
!                                                       dz8w,    &
!                                                       p8w,     &
	           				       theta   ,&
                   				       thetail ,&
                   				       pp      ,&
                   				       pi0     ,&
                                                       rho     ,&
                                                       pt      ,&
		                  rcp,rrp,rpp,rsp,rap,rgp,rhp






  real, dimension(ims:ime , jms:jme,npat), intent(in)  :: patch_area
!-srf
!   REAL,  DIMENSION( kms:kme, ims:ime , jms:jme )         ,    &
!          OPTIONAL                                         ,    &
!          INTENT(INOUT   ) ::                                   &
!               GDC,GDC2

   REAL, DIMENSION( ims:ime , jms:jme ),INTENT(IN) :: GSW,HT,rtgt,aot500
   REAL, DIMENSION( ims:ime , jms:jme ) :: XLAND
!srf-tmp
!   REAL, DIMENSION( ims:ime , jms:jme ),INTENT(INOUT) :: extra2d
!   REAL, DIMENSION( kms:kme ,ims:ime , jms:jme ),INTENT(INOUT) :: extra3d, extra3d2, extra3d3
!-srf
!   INTEGER, DIMENSION( ims:ime , jms:jme ),INTENT(IN) :: KPBL
!   INTEGER, DIMENSION( ims:ime , jms:jme ),INTENT(INOUT) :: k22_shallow, &
!                 kbcon_shallow,ktop_shallow
   INTEGER, DIMENSION( ims:ime , jms:jme )  :: KPBL
   INTEGER, DIMENSION( ims:ime , jms:jme )  :: k22_shallow &
                                              ,kbcon_shallow&
					      ,ktop_shallow
!-srf
!
   REAL, INTENT(IN   ) :: DT, DX,time
!
!-srf
!   REAL, DIMENSION( ims:ime , jms:jme ),                        &
!         INTENT(INOUT) ::           pratec,RAINCV, MASS_FLUX,   &
!                          APR_GR,APR_W,APR_MC,APR_ST,APR_AS,    &
!                         edt_out,APR_CAPMA,APR_CAPME,APR_CAPMI, &
!                         xmb_shallow
   REAL, DIMENSION( ims:ime , jms:jme ),                         &
         INTENT(INOUT) :: pratec,MASS_FLUX,                     &
                          APR_GR,APR_W,APR_MC,APR_ST,APR_AS,    &
                          xmb_shallow
   REAL, DIMENSION( ims:ime , jms:jme ),                      &
         INTENT(IN) :: weight_GR,weight_W,weight_MC,weight_ST,weight_AS

   REAL, DIMENSION( its:ite , jts:jte ) ::RAINCV,               &
                         edt_out,APR_CAPMA,APR_CAPME,APR_CAPMI
!-srf


!                         htop,hbot,xmb_shallow !SRF
!+lxz
 REAL, DIMENSION( ims:ime , jms:jme ) :: & !, INTENT(INOUT) ::       &
        HTOP,     &! highest model layer penetrated by cumulus since last reset in radiation_driver
        HBOT       ! lowest  model layer penetrated by cumulus since last reset in radiation_driver
!                  ! HBOT>HTOP follow physics leveling convention

!SRF
!   LOGICAL, DIMENSION( ims:ime , jms:jme ),                     &
!         INTENT(INOUT) ::                       CU_ACT_FLAG
   LOGICAL, DIMENSION( ims:ime , jms:jme )                      &
                        ::                       CU_ACT_FLAG
!SRF
!
! Optionals
!
   REAL, DIMENSION(kms:kme , ims:ime ,  jms:jme ),              &
         OPTIONAL,                                              &
         INTENT(IN) ::      RTHFTEN,  RQVFTEN,                  &
                            RTHBLTEN,RQVBLTEN

   REAL, DIMENSION(kms:kme , ims:ime ,  jms:jme ),              &
         OPTIONAL,                                              &
         INTENT(INOUT) ::                                       &
!                        cugd_tten,cugd_qvten,cugd_qcten,    &
                            cugd_ttens,cugd_qvtens


   REAL, DIMENSION(kms:kme , ims:ime ,  jms:jme ),              &
        OPTIONAL,                                               &
         INTENT(INOUT) ::                                       &
                                                   RTHCUTEN,    &
                                                   RQVCUTEN,    &
                                                   RQCCUTEN!srf,    &
                                                   !srf RQICUTEN
!
! Flags relating to the optional tendency arrays declared above
! Models that carry the optional tendencies will provdide the
! optional arguments at compile time; these flags all the model
! to determine at run-time whether a particular tracer is in
! use or not.
!
   LOGICAL, OPTIONAL ::                                      &
                                                   F_QV      &
                                                  ,F_QC      &
                                                  ,F_QR      &
                                                  ,F_QI      &
                                                  ,F_QS


!#if ( WRF_DFI_RADAR == 1 )
!
!  option of cap suppress:
!        do_capsuppress = 1   do
!        do_capsuppress = other   don't
!
!
!   INTEGER,      INTENT(IN   ) ,OPTIONAL   :: do_capsuppress
!   REAL, DIMENSION( ims:ime, jms:jme ),INTENT(IN   ),OPTIONAL  :: cap_suppress_loc
!   REAL, DIMENSION( its:ite ) :: cap_suppress_j
!#endif

! LOCAL VARS
!-srf
!    real,    dimension(ims:ime,jms:jme,1:ensdim),intent(inout) ::      &
!        xf_ens,pr_ens
     real,    dimension(ims:ime,jms:jme,1:ensdim) ::      &
        xf_ens,pr_ens
!-srf
     real,    dimension ( its:ite , jts:jte , 1:ensdim) ::      &
        massflni,xfi_ens,pri_ens
     REAL, DIMENSION( its:ite , jts:jte ) ::            MASSI_FLX,    &
                         APRi_GR,APRi_W,APRi_MC,APRi_ST,APRi_AS,    &
                         edti_out,APRi_CAPMA,APRi_CAPME,APRi_CAPMI,gswi
     real,    dimension (its:ite,kts:kte) ::                    &
        SUBT,SUBQ,OUTT,OUTQ,OUTQC,phh,subm,cupclw,dhdt,         &
        outts,outqs
     real,    dimension (its:ite)         ::                    &
        ccn,pret, ter11, aa0, fp,xlandi
!+lxz
     integer, dimension (its:ite) ::                            &
        kbcon, ktop,kpbli,k22s,kbcons,ktops,ierr
!.lxz
     integer, dimension (its:ite,jts:jte) ::                    &
        iact_old_gr
     integer :: iens,ibeg,iend,jbeg,jend,n,nn,ens4n
     integer :: ibegh,iendh,jbegh,jendh
     integer :: ibegc,iendc,jbegc,jendc
     integer,parameter :: iprint=0

!
! basic environmental input includes moisture convergence (mconv)
! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
! convection for this call only and at that particular gridpoint
!
     real,    dimension (its:ite,kts:kte) ::                    &
        T2d,q2d,PO,P2d,US,VS,tn,qo,tshall,qshall,rho2d
     real,    dimension (ips-2:ipe+2,kps:kpe,jps-2:jpe+2) ::    &
        ave_f_t,ave_f_q
     real,    dimension (its:ite,kts:kte,1:ens4) ::                    &
        omeg,tx,qx
     real, dimension (its:ite)            ::                    &
        Z1,PSUR,AAEQ,direction,cuten,umean,vmean,pmean,xmbs
     real, dimension (its:ite,1:ens4)     ::                    &
        mconv

   INTEGER :: i,j,k,ICLDCK,ipr,jpr,kr
   REAL    :: tcrit,tscl_KF,dp,dq,sub_spread,subcenter,exner,cpdTdt
   INTEGER :: itf,jtf,ktf,iss,jss,nbegin,nend
   INTEGER :: high_resolution
   REAL    :: rkbcon,rktop,r_liq,r_sol,tempk,fxc,dfxcdt,outxx,ccnclean    !-lxz
! ruc variable
   real, dimension (its:ite)            ::  tkm
! - relative humidity
   real,    dimension (its:ite,kts:kte) ::                    &
        rh2d,rs2d
   character*50 :: ierrc(its:ite)

  ! A. Betts for shallow convection: suggestion for the KF timescale < DELTAX  / 25 m/s
   tscl_kf=dx/25.
   ccn(:)=50.
  !
   high_resolution=0
   if(cugd_avedx.gt.1) high_resolution=1
   subcenter=0.
!  subcenter=1./float(cugd_avedx)
   sub_spread=max(1.,float(cugd_avedx*cugd_avedx-1))
   sub_spread=(1.-subcenter)/sub_spread
   !print*,'spread=',cugd_avedx,sub_spread;call flush(6)


   iens=1
   ipr=0
   jpr=0
   ipr=0
   jpr=0
!  if(itimestep.eq.8)then
!   ipr=37
!   jpr=16
!  endif
   IF ( periodic_x ) THEN ! ONLY FOR GLOBAL
      ibeg=max(its,ids)
      iend=min(ite,ide-1)
      ibegc=max(its,ids)
      iendc=min(ite,ide-1)
   ELSE
!srf      ibeg=max(its,ids)
!srf      iend=min(ite,ide-1)
!srf      ibegc=max(its,ids+4)
!srf      iendc=min(ite,ide-5)
      ibeg=its
      iend=ite
      ibegc=its
      iendc=ite
   END IF
   IF ( periodic_y ) THEN ! ONLY FOR GLOBAL
      jbeg=max(jts,jds)
      jend=min(jte,jde-1)
      jbegc=max(jts,jds)
      jendc=min(jte,jde-1)
   ELSE
!srf	  jbeg=max(jts,jds)
!srf	  jend=min(jte,jde-1)
!srf	  jbegc=max(jts,jds+4)
!srf	  jendc=min(jte,jde-5)
    jbeg=JTS
    jend=JTE
    jbegc=JTS
    jendc=JTE

   END IF
   do j=jts,jte
   do i=its,ite
 !  print*,i,j; call flush(6)
     k22_shallow(i,j)=0
     kbcon_shallow(i,j)=0
     ktop_shallow(i,j)=0
     xmb_shallow(i,j)=0
   enddo
   enddo
   tcrit=258.
   ave_f_t=0.
   ave_f_q=0.

   itf=MIN(ite,ide-1)
   ktf=MIN(kte,kde-1)
   jtf=MIN(jte,jde-1)


!print*,'g3d start 2 ',cugd_avedx,high_resolution;call flush(6)
!print*,'jts jtf its itf kts ktf', jts, jtf ,its ,itf, kts ,ktf;call flush(6)
!
!#if ( EM_CORE == 1 )
     if(high_resolution.eq.1)then
!
! calculate these on the halo...the incominh tendencies have been exchanged on a 24pt halo
! only neede for high resolution run
!
     ibegh=its
     jbegh=jts
     iendh=ite
     jendh=jte
     if(its.eq.ips)ibegh=max(its-1,ids)
     if(jts.eq.jps)jbegh=max(jts-1,jds)
     if(jte.eq.jpe)jendh=min(jte+1,jde-1)
     if(ite.eq.ipe)iendh=min(ite+1,ide-1)
        DO J = jbegh,jendh
        DO k= kts,ktf
	kr=k+1
        DO I= ibegh,iendh
          !ave_f_t(i,k,j)=(rthften(i-1,k,j-1)+rthften(i-1,k,j) + rthften(i-1,k,j+1)+ &
          !               rthften(i,k,j-1)   +rthften(i,k,j)   +rthften(i,k,j+1)+         &
          !               rthften(i+1,k,j-1) +rthften(i+1,k,j) +rthften(i+1,k,j+1))/9.
          !ave_f_q(i,k,j)=(rqvften(i-1,k,j-1)+rqvften(i-1,k,j) + rqvften(i-1,k,j+1)+ &
          !               rqvften(i,k,j-1)   +rqvften(i,k,j)   +rqvften(i,k,j+1)+         &
          !               rqvften(i+1,k,j-1) +rqvften(i+1,k,j) +rqvften(i+1,k,j+1))/9.
         ave_f_t(i,k,j)=rthften(kr,i,j)
         ave_f_q(i,k,j)=rqvften(kr,i,j)
        ENDDO
        ENDDO
        ENDDO
     endif  ! endif of high_resolution.eq.1
!#endif
     DO 100 J = jts,jtf  ! J LOOP
     DO n= 1,ensdim
     DO I= its,itf
       xfi_ens(i,j,n)=0.
       pri_ens(i,j,n)=0.
!      xfi_ens(i,j,n)=xf_ens(i,j,n)
!      pri_ens(i,j,n)=pr_ens(i,j,n)
     ENDDO
     ENDDO
     DO I= its,itf
        kbcon(i)=0
        ktop(i)=0
        tkm(i)=0.
        HBOT(I,J)  =REAL(KTE)
        HTOP(I,J)  =REAL(KTS)
        iact_old_gr(i,j)=0
        mass_flux(i,j)=0.
        massi_flx(i,j)=0.
        raincv(i,j)=0.
        pratec (i,j)=0.
        edt_out(i,j)=0.
        edti_out(i,j)=0.
        gswi(i,j)=gsw(i,j)
        xland(i,j)       = patch_area(i,j,1) !flag < 1 para land
	                                     !flag  =1 para water
        xlandi(i)=xland(i,j)
!-srf tmp
!        extra2d(i,j)=-9999.
!        extra3d(:,i,j)=0.
!        extra3d2(:,i,j)=0.
!        extra3d3(:,i,j)=0.
!-srf tmp

!---srf
!        APRi_GR(i,j)=apr_gr(i,j)
!        APRi_w(i,j)=apr_w(i,j)
!        APRi_mc(i,j)=apr_mc(i,j)
!        APRi_st(i,j)=apr_st(i,j)
!        APRi_as(i,j)=apr_as(i,j)
!        APRi_capma(i,j)=apr_capma(i,j)
!        APRi_capme(i,j)=apr_capme(i,j)
!        APRi_capmi(i,j)=apr_capmi(i,j)
!----srf
        CU_ACT_FLAG(i,j) = .true.
     ENDDO
     if(autoconv == 2) then
       DO I= its,itf
         ccn(i) = max( 100., ( 370.37*(0.01+MAX(0.,aot500(i,j))))**1.555 )
	 !testando
	 !ccn(i)=170.

       ENDDO
     ELSE
       DO I= its,itf
         ccn(i) = 100.
       ENDDO
     ENDIF
!srf
!     do k=kts,kte
!     DO I= its,itf
!       cugd_tten(i,k,j)=0.
!       cugd_ttens(i,k,j)=0.
!       cugd_qvten(i,k,j)=0.
!       cugd_qvtens(i,k,j)=0.
!       cugd_qcten(i,k,j)=0.
!     ENDDO
!     ENDDO

     DO n=1,ens4
     DO I= its,itf
        mconv(i,n)=0.
     ENDDO
     do k=kts,kte
     DO I= its,itf
         omeg(i,k,n)=0.
         tx(i,k,n)=0.
         qx(i,k,n)=0.
     ENDDO
     ENDDO
     ENDDO
     DO k=1,ensdim
     DO I= its,itf
        massflni(i,j,k)=0.
     ENDDO
     ENDDO
     !  put hydrostatic pressure on half levels
     DO K=kts,ktf
     DO I=ITS,ITF
!srf     phh(i,k) = p(i,k,j)
         kr=k+1
         phh(i,k) =((pp(kr,i,j)+pi0(kr,i,j))/cp)**cpor*p00  !*1.e-2
     ENDDO
     ENDDO

     DO I=ITS,ITF
!srf     PSUR(I)=p8w(I,1,J)*.01
         PSUR(I) = .5*( ((pp(1,i,j)+pi0(1,i,j))/cp)**cpor*p00 +  &
                        ((pp(2,i,j)+pi0(2,i,j))/cp)**cpor*p00 )*1.e-2

!        PSUR(I)=p(I,1,J)*.01
         TER11(I)=HT(i,j)
         aaeq(i)=0.
         direction(i)=0.
         pret(i)=0.
         umean(i)=0.
         vmean(i)=0.
         pmean(i)=0.
         kpbli(i)=kpbl(i,j)
     ENDDO
     if(j.eq.jpr)write(0,*)'psur(ipr),ter11(ipr),kpbli(ipr)'
     if(j.eq.jpr)write(4,*)psur(ipr),ter11(ipr)!,kpbli(ipr),r_v
     DO K=kts,ktf
     DO I=ITS,ITF
         po(i,k)=phh(i,k)*.01
         subm(i,k)=0.
         P2d(I,K)=PO(i,k)
         !srf -----
	 kr=k+1
	 US(I,K) =.5*( u(kr,i,j) + u(kr,i-1,j) )
         VS(I,K) =.5*( v(kr,i,j) + v(kr,i,j-1) )
         T2d(I,K)  = theta(kr,i,j)*(pp(kr,i,j)+pi0(kr,i,j))/cp
         q2d(I,K)  = rv(kr,i,j)
	 rs2d(i,k) = rs(phh(i,k),T2d(i,k))
	 rh2d(i,k) = min(1.,max(0.,q2d(i,k)/rs2d(i,k)))
	 rho2d(i,k)= 0.001*rho(kr,i,j) ! unidades: g/cm�


!---lixo
!	 extra3d2(k+1,i,j)=rs2d  (i,k)
!---lixo
	 !print*,'rs,rv,rh=',k,rs2d(i,k),q2d(i,k),rh2d(i,k),T2d(i,k);call flush(6)

         !srf -----

         IF(Q2d(I,K).LT.1.E-08)Q2d(I,K)=1.E-08
         SUBT(I,K)=0.
         SUBQ(I,K)=0.
         OUTT(I,K)=0.
         OUTQ(I,K)=0.
         OUTQC(I,K)=0.
         OUTTS(I,K)=0.
         OUTQS(I,K)=0.

!srf     TN(I,K)=t2d(i,k)+RTHFTEN(i,k,j)*dt
         exner= pp(kr,i,j)+pi0(kr,i,j)
         cpdTdt= exner*RTHFTEN(kr,i,j) + theta(kr,i,j)*pt(kr,i,j)
         TN(I,K)= T2d(I,K) + ( cpdTdt/cp )*dt


!srf	 QO(I,K)=q2d(i,k)+RQVFTEN(i,k,j)*dt
	 QO(I,K)=q2d(i,k)+RQVFTEN(kr,i,j)*dt

!srf - ver no futuro quando ligar o shallow
!srf - for shallow convection only
!         TSHALL(I,K)=t2d(i,k)+RTHBLTEN(i,k,j)*pi(i,k,j)*dt
!         DHDT(I,K)=cp*RTHBLTEN(i,k,j)*pi(i,k,j)+ XLV*RQVBLTEN(i,k,j)
!	 QSHALL(I,K)=q2d(i,k)+RQVBLTEN(i,k,j)*dt
!srf - ver no futuro


!-srf- ver no futuro
!	 if(high_resolution.eq.1)then
!            TN(I,K)=t2d(i,k)+ave_f_t(i,k,j)*dt
!            QO(I,K)=q2d(i,k)+ave_f_q(i,k,j)*dt
!         endif
!-srf- ver no futuro

	 IF(TN(I,K).LT.200.)    TN(I,K)=T2d(I,K)
         IF(QO(I,K).LT.1.E-08)  QO(I,K)=1.E-08
     ENDDO
     ENDDO
     ens4n=0
     nbegin=0
     nend=0
     if(ens4_spread.gt.1)then
     nbegin=-ens4_spread/2
     nend=ens4_spread/2
     endif
     !-----------  ENSEMBLE FORCING
     do nn=nbegin,nend,1
       jss=j !max(j+nn,jds+0)
       jss=j !min(jss,jde-1)
       do n=nbegin,nend,1
         ens4n=ens4n+1
         DO K=kts,ktf
         DO I=ITS,ITF
          iss=i !max(i+n,ids+0)
          iss=i !min(iss,ide-1)
          kr=k+1
!srf      omeg(I,K,ens4n)= -g*rho(i,k,j)*w(iss,k,jss)
          omeg(I,K,ens4n)= -g*rho(kr,I,j)*w(kr,i,j)
!         omeg(I,K,ens4n)= -g*rho(i,k,j)*w(i,k,j)
!
!srf     Tx(I,K,ens4n)=t2d(i,k)+RTHFTEN(iss,k,jss)*dt
         Tx(I,K,ens4n)=TN(I,K) ! for now t2d(i,k)+RTHFTEN(kr,iss,jss)*dt

!        Tx(I,K,ens4n)=t2d(i,k)+RTHFTEN(i,k,j)*dt

!srf -ver no futuro
!         if(high_resolution.eq.1)Tx(I,K,ens4n)=t2d(i,k)+ave_f_t(iss,k,jss)*dt
!        IF(Tx(I,K,ens4n).LT.200.)Tx(I,K,ens4n)=T2d(I,K)
!srf -ver no futuro
!
!srf     Qx(I,K,ens4n)=q2d(i,k)+RQVFTEN(iss,k,jss)*dt
         Qx(I,K,ens4n)=QO(i,k)
!        Qx(I,K,ens4n)=q2d(i,k)+RQVFTEN(kr,iss,jss)*dt
!        Qx(I,K,ens4n)=q2d(i,k)+RQVFTEN(i,k,j)*dt

!srf -ver no futuro
!         if(high_resolution.eq.1)qx(I,K,ens4n)=q2d(i,k)+ave_f_q(iss,k,jss)*dt
!         IF(Qx(I,K,ens4n).LT.1.E-08)Qx(I,K,ens4n)=1.E-08
!srf -ver no futuro

        enddo
        enddo
      enddo !n
      enddo !nn
      do k=  kts+1,ktf-1
      DO I = its,itf
         if((p2d(i,1)-p2d(i,k)).gt.150.and.p2d(i,k).gt.300)then
            dp=-.5*(p2d(i,k+1)-p2d(i,k-1))
            umean(i)=umean(i)+us(i,k)*dp
            vmean(i)=vmean(i)+vs(i,k)*dp
            pmean(i)=pmean(i)+dp
         endif
      enddo
      enddo
      DO I = its,itf
         umean(i)=umean(i)/pmean(i)
         vmean(i)=vmean(i)/pmean(i)
         direction(i)=(atan2(umean(i),vmean(i))+3.1415926)*57.29578
         if(direction(i).gt.360.)direction(i)=direction(i)-360.
      ENDDO
      do n=1,ens4
      DO K=kts,ktf-1
      DO I = its,itf
        dq=(q2d(i,k+1)-q2d(i,k))
        mconv(i,n)=mconv(i,n)+omeg(i,k,n)*dq/g
      enddo
      ENDDO
      ENDDO
      do n=1,ens4
      DO I = its,itf
       if(mconv(i,n).lt.0.)mconv(i,n)=0.
      ENDDO
      ENDDO
!
!---- CALL CUMULUS PARAMETERIZATION
!
!#if ( WRF_DFI_RADAR == 1 )
!      if(do_capsuppress == 1 ) then
!        DO I= its,itf
!            cap_suppress_j(i)=cap_suppress_loc(i,j)
!        ENDDO
!      endif
!#endif
      CALL CUP_enss_3d(&
           outqc,   & !1
	   j,   &
	   AAEQ,   &
	   T2d,   &
	   Q2d,   &
	   TER11,   &
	   subm,   &
	   TN,   &
	   QO,   &
	   PO,   &
	   PRET,   &
           P2d,   &
	   OUTT,   &
	   OUTQ,   &
	   DT,   &
	   itimestep,   &
	   tkm,   &
	   PSUR,   &
	   US,   &
	   VS,   &
	   tcrit,   &
	   iens,   &
	   tx,   &
	   qx,   &
           tshall,   &
	   qshall,   &
	   kpbli,   &
	   DHDT,   &
	   outts,   &
	   outqs,   &
	   tscl_kf,   &

	   k22s,   &
	   kbcons,   &
	   ktops,   &
	   xmbs,   &
	   ccn,   &
	   ccnclean,   &
	   rho2d,   &
	   dx,   &

	   mconv,   &
	   massflni,   &
	   iact_old_gr,   &
	   omeg,   &
	   direction,   &
	   MASSi_FLX,   &
           maxiens,   &
	   maxens,   &
	   maxens2,   &
	   maxens3,   &
	   ensdim,   &

	   APR_GR,APR_W,APR_MC,APR_ST,APR_AS,   &
           APR_CAPMA,APR_CAPME,APR_CAPMI,   &
	   kbcon,   &
	   ktop,   &
	   cupclw,   &
!
           xfi_ens,   &
	   pri_ens,   &
	   XLANDi,   &
	   gswi,   &
	   subt,   &
	   subq,   &
           xlv,   r_v,   cp,   g,   &
	   ichoice,   &
	   ipr,   &
	   jpr,   &
	   ierrc, &
	   ierr,  &
	   ens4,  &
	   high_resolution,   &
           beta,   &
	   ishallow_g3,   &
	   autoconv,   &
	   aeroevap,   &
	   itf,jtf,ktf,its,ite, jts,jte, kts,kte,         &
	   weight_GR,weight_W,weight_MC,weight_ST,weight_AS,   &
	   training)

!      CALL CUP_enss_3d(&
!           outqc,j,AAEQ,T2d,Q2d,TER11,subm,TN,QO,PO,PRET,   &
!           P2d,OUTT,OUTQ,DT,itimestep,tkm,PSUR,US,VS,tcrit, &
!   iens,tx,qx,tshall,qshall,kpbli,DHDT,outts,outqs, &
!	   tscl_kf, k22s,kbcons,ktops,xmbs, ccn,ccnclean,rho2d,dx,               &
!
!	   mconv,massflni,iact_old_gr,omeg,direction,MASSi_FLX,         &
!           maxiens,maxens,maxens2,maxens3,ensdim,                       &
!
!!srf       APRi_GR,APRi_W,APRi_MC,APRi_ST,APRi_AS,		  &
!!srf       APRi_CAPMA,APRi_CAPME,APRi_CAPMI,kbcon,ktop,cupclw,    &
!
!	   APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                     &
!           APR_CAPMA,APR_CAPME,APR_CAPMI,kbcon,ktop,cupclw,       &
!
!           xfi_ens,pri_ens,XLANDi,gswi,subt,subq,                 &
!           xlv,r_v,cp,g,ichoice,ipr,jpr,ierrc,ierr,ens4,high_resolution,     &
!           beta,ishallow_g3,autoconv,aeroevap,itf,jtf,ktf,                               &
!           its,ite, jts,jte, kts,kte,           &
!	   weight_GR,weight_W,weight_MC,weight_ST,weight_AS,training&
!
!	   !,extra2d ,extra3d ,extra3d2,extra3d3&
!#if ( WRF_DFI_RADAR == 1 )
!           ,do_capsuppress,cap_suppress_j                        &
!#endif
!                                                                  )
!
            if(j.lt.jbegc.or.j.gt.jendc)go to 100
            DO I=ibegc,iendc
              xmb_shallow(i,j)=xmbs(i)

              k22_shallow(i,j)=k22s(i)
              kbcon_shallow(i,j)=kbcons(i)
              ktop_shallow(i,j)=ktops(i)
              cuten(i)=0.
              if(pret(i).gt.0.)then
                 cuten(i)=1.
!                raincv(i,j)=pret(i)*dt
              endif
            ENDDO
            !if(j.eq.jpr)write(0,*)'precip,ktop,kbcon = ',pret(ipr),ktop(ipr),kbcon(ipr)
            !
	    DO I=ibegc,iendc
             DO K=kts,ktf
	       kr=k+1
	       !- subsidence tendencies
               cugd_ttens (Kr,i,J)=subt(i,k)*cuten(i)*sub_spread
               cugd_qvtens(Kr,i,J)=subq(i,k)*cuten(i)*sub_spread
!!             cugd_tten(I,K,J)=outt(i,k)*cuten(i)
!!             cugd_qvten(I,K,J)=outq(i,k)*cuten(i)
!!             cugd_tten (Kr,i,J)=outts(i,k)+outt(i,k)*cuten(i)
!!             cugd_qvten(Kr,i,J)=outqs(i,k)+outq(i,k)*cuten(i)
!!             cugd_qcten(Kr,i,J)=outqc(i,k)*cuten(i)
               !-
	       RTHCUTEN(Kr,i,J)=outts(i,k)+outt(i,k)*cuten(i)
               RQVCUTEN(Kr,i,J)=outqs(i,k)+outq(i,k)*cuten(i)
               RQCCUTEN(Kr,i,J)=outqc(i,k)*cuten(i)

               !- original version
!              RTHCUTEN(Kr,i,J)=(subt(i,k)*sub_spread+outt(i,k))*cuten(i)+outts(i,k)
!              RQVCUTEN(Kr,i,J)=(subq(i,k)*sub_spread+outq(i,k))*cuten(i)+outqs(i,k)
!              RQCCUTEN(Kr,i,J)=outqc(i,k)*cuten(i)

             ENDDO
!00000000000000000000000000000000000000000000000000000
    	    if(iprint==1 .and. pret(i) > 1.e-6) then
	      write(mynum+1,*) "--------------",i,j,'=>',i+i0,j+j0
	      write(mynum+1,*)psur(i),ter11(i),pret(i),time
	     do K=kts,ktf

	      write(mynum+1,123)k,p2d(i,k),t2d(i,k),tn(i,k),q2d(i,k),QO(i,k),US(i,k),Vs(i,k)&
	                    ,-g*rho(k+1,I,j)*w(k+1,i,j)
	      call flush(mynum+1)
	     enddo
	     write(mynum+1,*) "-----------------------------"
	     do K=kts,ktf
	       write(mynum+1,125)k,86400.*(subt(i,k)*cuten(i)+outts(i,k)+outt(i,k)*cuten(i))&
	                        ,86400.*(subq(i,k)*cuten(i)+outqs(i,k)+outq(i,k)*cuten(i))&
	                        ,86400.*outqc(i,k)*cuten(i)&
				,86400.*subt(i,k)*cuten(i)
	      call flush(mynum+1)
	     enddo
	     !stop 3333
 125  format(1x,i2,4(1x,e12.4))
 123  format(1x,i2,f8.0,1x,2(1x,f8.3),5(1x,e12.4))

	    endif
!00000000000000000000000000000000000000000000000000000
            ENDDO
!          ! Converte tend da temperatura (OUTT) em tend de theta (OUTTEM)
!          ! cp*T=Pi*Theta => cp dT/dt = Theta*dPi/dt + Pi*dTheta/dt,
!
!	   if(level <=2) then
!
!	    DO I=ibegc,iendc
!	     DO K=kts,ktf
!	       kr=k+1
!
!		 ! Converte tend da temperatura (OUTT) em tend de theta (OUTTEM)
!		 ! cp*T=Pi*Theta => cp dT/dt = Theta*dPi/dt + Pi*dTheta/dt,
!		 ! Exner's function = pp(kr,i,j)+pi0(kr,i,j)
!		  exner 	 = pp(kr,i,j) + pi0(kr,i,j)
!		 ! tendencia do theta devida a conv profunda
!		 RTHCUTEN(Kr,i,J) = cp/exner * RTHCUTEN(Kr,i,J) - theta(kr,i,j)*pt(kr,i,j)/exner
!		 !RQVCUTEN(Kr,i,J) = RQVCUTEN(Kr,i,J)+ RQCCUTEN(Kr,i,J)
!		 !RQCCUTEN(Kr,i,J) = 0.
!
!	    ENDDO
!	    RTHCUTEN(1,i,J)=RTHCUTEN(2,i,J)
!	   ENDDO
!	    elseif(level > 2) then
!
!
!	    DO I=ibegc,iendc
!	     DO K=kts,ktf
!	       kr=k+1
!	       outxx=(subt(i,k)*sub_spread+outt(i,k))*cuten(i)+outts(i,k)
!		   ! converte tend da temperatura (outt) em tend de theta (outtem)
!		    ! cp*T=Pi*Theta => cp dT/dt = Theta*dPi/dt + Pi*dTheta/dt,
!		    ! Exner's function = pp(kr,i,j)+pi0(kr,i,j)
!		    exner= pp(kr,i,j) + pi0(kr,i,j)
!		    ! tendencia do theta  devida a conv profunda
!		    RTHCUTEN (kr,i,j) = cp/exner * RTHCUTEN(kr,i,j) - theta(kr,i,j)*pt(kr,i,j)/exner
!
!		    ! tendencia do theta_il devida a conv profunda
!		    r_liq= max(0.,rcp(kr,i,j) + rrp(kr,i,j))
!
!		    r_sol= max(0.,rsp(kr,i,j)+rpp(kr,i,j)+	&
!				 rap(kr,i,j)+rgp(kr,i,j)+  &
!				 rhp(kr,i,j))
!
!		   tempk = theta(kr,i,j)*(exner)/cp ! air temp (Kelvin)
!
!		   if(tempk.le.253) then
!		     fxc =   (2.5e6*r_liq+2.83e6*r_sol)/(cp*amax1(tempk,253.))
!!		     dfxcdt = 2.83e6*OUTQC(I,K)*cuten(i)/(cp*amax1(tempk,253.))
!		     dfxcdt = 2.83e6*RQCCUTEN(Kr,i,J)/(cp*amax1(tempk,253.))
!		     RTHCUTEN (kr,i,j) = (1./(1+fxc))*( RTHCUTEN (kr,i,j) - thetail(kr,i,j)*dfxcdt )
!
!		   else
!
!		     fxc =   (2.5e6*r_liq+2.83e6*r_sol)/(cp*amax1(tempk,253.))
!		     dfxcdt = 2.5e6*RQCCUTEN(Kr,i,J)/(cp*amax1(tempk,253.)) - &
!!		     dfxcdt = 2.5e6*OUTQC(I,K)*cuten(i)/(cp*amax1(tempk,253.)) - &
!!			      fxc/(cp*amax1(tempk,253.)) * cp * OUTT(I,K)
!			      fxc/(cp*amax1(tempk,253.)) * cp * OUTxx
!
!		      RTHCUTEN (kr,i,j) = (1./(1+fxc))*( RTHCUTEN (kr,i,j) - thetail(kr,i,j)*dfxcdt )
!
!		   endif
!
!		     ! tendencia da vapor d'agua devida a conv profunda
!		     ! outrt	(kr,i,j) = OUTQ(I,K) !+ OUTQC(I,K)
!		     ! tendencia da agua condensada devida a conv profunda
!		     !if (CATT==1)sgrell3_3d(kr,i,j)= OUTQC(I,K)
!
!		      !RQVCUTEN(Kr,i,J) = RQVCUTEN(Kr,i,J)
!		      !RQCCUTEN(Kr,i,J) = outqc(I,K)*cuten(i)	 !<<<<<< ???????????
!
!	     ENDDO
!	     ENDDO
!	    endif
!40 continue

           DO I=ibegc,iendc
             if(pret(i).gt.0.)then
		 !raincv(i,j)=pret(i)*dt

!--              rainfall output
        	 pratec(i,j)=pret(i)
!--
                 rkbcon = kte+kts - kbcon(i)
                 rktop  = kte+kts -  ktop(i)
                 if (ktop(i)  > HTOP(i,j)) HTOP(i,j) = ktop(i)+.001
                 if (kbcon(i) < HBOT(i,j)) HBOT(i,j) = kbcon(i)+.001
              else
	         RTHCUTEN(:,i,j)=0.
	         RQVCUTEN(:,i,j)=0.
                 RQCCUTEN(:,i,j)=0.
	      endif
            ENDDO
            DO n= 1,ensdim
            DO I= ibegc,iendc
              xf_ens(i,j,n)=xfi_ens(i,j,n)
              pr_ens(i,j,n)=pri_ens(i,j,n)
            ENDDO
            ENDDO
            DO I= ibegc,iendc
!srf
!               APR_GR(i,j)=apri_gr(i,j)
!               APR_w(i,j)=apri_w(i,j)
!               APR_mc(i,j)=apri_mc(i,j)
!               APR_st(i,j)=apri_st(i,j)
!               APR_as(i,j)=apri_as(i,j)
!               APR_capma(i,j)=apri_capma(i,j)
!               APR_capme(i,j)=apri_capme(i,j)
!               APR_capmi(i,j)=apri_capmi(i,j)
!srf
                mass_flux(i,j)=massi_flx(i,j)
!srf            edt_out(i,j)=edti_out(i,j)
            ENDDO
!
!
!	    IF(PRESENT(RQCCUTEN)) THEN
!              IF ( F_QC ) THEN
!                DO K=kts,ktf
!		kr=k+1
!                DO I=ibegc,iendc
!                   RQCCUTEN(I,K,J)=outqc(I,K)*cuten(i)
!                   IF ( PRESENT( GDC ) ) GDC(Kr,i,J)=CUPCLW(I,K)*cuten(i)
!                   IF ( PRESENT( GDC2 ) ) GDC2(Kr,i,J)=0.
!                ENDDO
!                ENDDO
!              ENDIF
!            ENDIF
!
!......     QSTEN STORES GRAUPEL TENDENCY IF IT EXISTS, OTHERISE SNOW (V2)
!
!            IF(PRESENT(RQICUTEN).AND.PRESENT(RQCCUTEN))THEN
!              IF (F_QI) THEN
!                DO K=kts,ktf
!		kr=k+1
!                  DO I=ibegc,iendc
!                   if(t2d(i,k).lt.258.)then
!                      RQICUTEN(I,K,J)=outqc(I,K)*cuten(i)
!                      cugd_qcten(i,k,j)=0.
!                      RQCCUTEN(I,K,J)=0.
!                      IF ( PRESENT( GDC2 ) ) GDC2(kr,i,J)=CUPCLW(I,K)*cuten(i)
!                   else
!                      RQICUTEN(I,K,J)=0.
!                      RQCCUTEN(I,K,J)=outqc(I,K)*cuten(i)
!                      IF ( PRESENT( GDC ) ) GDC(kr,i,J)=CUPCLW(I,K)*cuten(i)
!                   endif
!                ENDDO
!                ENDDO
!              ENDIF
!            ENDIF
!
!-srf
 100    continue

   END SUBROUTINE G3DRV
!----------------------------------------------------------------------------
   SUBROUTINE CUP_enss_3d(OUTQC,J,AAEQ,T,Q,Z1,sub_mas,                    &
              TN,QO,PO,PRE,P,OUTT,OUTQ,DTIME,ktau,tkmax,PSUR,US,VS,    &
              TCRIT,iens,tx,qx,                                        &
              tshall,qshall,kpbl,dhdt,outts,outqs,tscl_kf,             &
              k23,kbcon3,ktop3,xmb3,ccn,ccnclean,rho,dx,               &
              mconv,massfln,iact,                                      &
              omeg,direction,massflx,maxiens,                          &
              maxens,maxens2,maxens3,ensdim,                           &
              APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                       &
              APR_CAPMA,APR_CAPME,APR_CAPMI,kbcon,ktop,cupclw,         &   !-lxz
              xf_ens,pr_ens,xland,gsw,subt,subq,                       &
              xl,rv,cp,g,ichoice,ipr,jpr,ierrc,ierr,ens4,high_resolution,  &
	      beta,ishallow_g3,autoconv,aeroevap,itf,jtf,ktf,              &
              its,ite, jts,jte, kts,kte,                                   &
              weight_GR,weight_W,weight_MC,weight_ST,weight_AS,training    &
!#if ( WRF_DFI_RADAR == 1 )
!                 ! Optional CAP suppress option
!                     ,do_capsuppress,cap_suppress_j                  &
!#endif
                                                )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        autoconv,aeroevap,itf,jtf,ktf,ktau,training,                   &
        its,ite, jts,jte, kts,kte,ipr,jpr,ens4,high_resolution
     integer, intent (in   )              ::                           &
        j,ensdim,maxiens,ishallow_g3,maxens,maxens2,maxens3,ichoice,iens
  !
  !
  !
     real,    dimension (its:ite,jts:jte,1:ensdim)                     &
        ,intent (inout)                   ::                           &
        massfln,xf_ens,pr_ens
     real,    dimension (its:ite,jts:jte)                              &
        ,intent (inout )                  ::                           &
               APR_GR,APR_W,APR_MC,APR_ST,APR_AS,APR_CAPMA,     &
               APR_CAPME,APR_CAPMI,massflx
    real, dimension( its:ite , jts:jte ),                              &
         intent(in) :: weight_GR,weight_W,weight_MC,weight_ST,weight_AS
     real,    dimension (its:ite,jts:jte)                              &
        ,intent (in   )                   ::                           &
               gsw
     integer, dimension (its:ite,jts:jte)                              &
        ,intent (in   )                   ::                           &
        iact
  ! outtem = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout  )                   ::                           &
        DHDT,OUTT,OUTQ,OUTQC,subt,subq,sub_mas,cupclw,outts,outqs
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pre,xmb3
     integer,    dimension (its:ite)                                   &
        ,intent (out  )                   ::                           &
        kbcon,ktop,k23,kbcon3,ktop3
     integer,    dimension (its:ite)                                   &
        ,intent (in  )                   ::                           &
        kpbl
  !
  ! basic environmental input includes moisture convergence (mconv)
  ! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
  ! convection for this call only and at that particular gridpoint
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        rho,T,PO,P,US,VS,tn,tshall,qshall
     real,    dimension (its:ite,kts:kte,1:ens4)                       &
        ,intent (inout   )                   ::                           &
        omeg,tx,qx
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
         Q,QO
     real, dimension (its:ite)                                         &
        ,intent (in   )                   ::                           &
        ccn,Z1,PSUR,AAEQ,direction,tkmax,xland
     real, dimension (its:ite,1:ens4)                                         &
        ,intent (in   )                   ::                           &
        mconv


       real                                                            &
        ,intent (in   )                   ::                           &
        beta,dx,ccnclean,dtime,tcrit,xl,cp,rv,g,tscl_kf

!#if ( WRF_DFI_RADAR == 1 )
!
!  option of cap suppress:
!        do_capsuppress = 1   do
!        do_capsuppress = other   don't
!
!
!   INTEGER,      INTENT(IN   ) ,OPTIONAL   :: do_capsuppress
!   REAL, DIMENSION( its:ite ),INTENT(IN   ) ,OPTIONAL   :: cap_suppress_j
!#endif

!
!  local ensemble dependent variables in this routine
!
     real,    dimension (its:ite,1:maxens)  ::                         &
        xaa0_ens
     real,    dimension (1:maxens)  ::                                 &
        mbdt_ens
     real,    dimension (1:maxens2) ::                                 &
        edt_ens
     real,    dimension (its:ite,1:maxens2) ::                         &
        edtc
     real,    dimension (its:ite,kts:kte,1:maxens2) ::                 &
        dellat_ens,dellaqc_ens,dellaq_ens,pwo_ens,subt_ens,subq_ens
!
!
!
!***************** the following are your basic environmental
!                  variables. They carry a "_cup" if they are
!                  on model cloud levels (staggered). They carry
!                  an "o"-ending (z becomes zo), if they are the forced
!                  variables. They are preceded by x (z becomes xz)
!                  to indicate modification by some typ of cloud
!
  ! z           = heights of model levels
  ! q           = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! p           = environmental pressure
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  ! z_cup       = heights of model cloud levels
  ! q_cup       = environmental q on model cloud levels
  ! qes_cup     = saturation q on model cloud levels
  ! t_cup       = temperature (Kelvin) on model cloud levels
  ! p_cup       = environmental pressure
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! gamma_cup = gamma on model cloud levels
!
!
  ! hcd = moist static energy in downdraft
  ! zd normalized downdraft mass flux
  ! dby = buoancy term
  ! entr = entrainment rate
  ! zd   = downdraft normalized mass flux
  ! entr= entrainment rate
  ! hcd = h in model cloud
  ! bu = buoancy term
  ! zd = normalized downdraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! mentr_rate = entrainment rate
  ! qcd = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! pwd = evaporate at that level
  ! pwev = total normalized integrated evaoprate (I2)
  ! entr= entrainment rate
  ! z1 = terrain elevation
  ! entr = downdraft entrainment rate
  ! jmin = downdraft originating level
  ! kdet = level above ground where downdraft start detraining
  ! psur        = surface pressure
  ! z1          = terrain elevation
  ! pr_ens = precipitation ensemble
  ! xf_ens = mass flux ensembles
  ! massfln = downdraft mass flux ensembles used in next timestep
  ! omeg = omega from large scale model
  ! mconv = moisture convergence from large scale model
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! dir     = "storm motion"
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! iact_gr_old = flag to tell where convection was active
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! dby = buoancy term
  ! ktop = cloud top (output)
  ! xmb    = total base mass flux
  ! hc = cloud moist static energy
  ! hkb = moist static energy at originating level
  ! mentr_rate = entrainment rate
     real,    dimension (its:ite,kts:kte) ::                           &
        he3,hes3,qes3,z3,zdo3,zu3_0,hc3_0,dby3_0,                      &
        qes3_cup,q3_cup,he3_cup,hes3_cup,z3_cup,gamma3_cup,t3_cup,     &
        xhe3,xhes3,xqes3,xz3,xt3,xq3,                                  &
        xqes3_cup,xq3_cup,xhe3_cup,xhes3_cup,xz3_cup,xgamma3_cup,      &
        xt3_cup,                                                       &
        xdby3,xqc3,xhc3,xqrc3,xzu3,                                    &
        dby3,qc3,pw3,hc3,qrc3,zu3,cd3,DELLAH3,DELLAQ3,                 &
        dsubt3,dsubq3,DELLAT3,DELLAQC3

     real,    dimension (its:ite,kts:kte) ::                           &
        entr_rate_2d,mentrd_rate_2d,he,hes,qes,z,                                     &
        heo,heso,qeso,zo,                                              &
        xhe,xhes,xqes,xz,xt,xq,                                        &

        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,      &
        qeso_cup,qo_cup,heo_cup,heso_cup,zo_cup,po_cup,gammao_cup,     &
        tn_cup,                                                        &
        xqes_cup,xq_cup,xhe_cup,xhes_cup,xz_cup,xp_cup,xgamma_cup,     &
        xt_cup,                                                        &

        xlamue,dby,qc,qrcd,pwd,pw,hcd,qcd,dbyd,hc,qrc,zu,zd,clw_all,   &
        dbyo,qco,qrcdo,pwdo,pwo,hcdo,qcdo,dbydo,hco,qrco,zuo,zdo,      &
        xdby,xqc,xqrcd,xpwd,xpw,xhcd,xqcd,xhc,xqrc,xzu,xzd,            &

  ! cd  = detrainment function for updraft
  ! cdd = detrainment function for downdraft
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble

        cd,cdd,DELLAH,DELLAQ,DELLAT,DELLAQC,dsubt,dsubq

  ! aa0 cloud work function for downdraft
  ! edt = epsilon
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects (ensemble dependent)
  ! edt     = epsilon

     real,    dimension (its:ite) ::                                   &
       aa3_0,aa3,hkb3,qkb3,pwav3,bu3,xaa3,xhkb3,                       &
       hkb3_0,edt,edto,edtx,AA1,AA0,XAA0,HKB,                          &
       HKBO,aad,XHKB,QKB,QKBO,edt3,                                    &
       XMB,XPWAV,XPWEV,PWAV,PWEV,PWAVO,                                &
       PWEVO,BU,BUO,cap_max,xland1,                                    &
       cap_max_increment,closure_n,cap_max3,psum,psumh,sig,zuhe
     real,    dimension (its:ite,1:ens4) ::                                   &
        axx
     integer,    dimension (its:ite) ::                                &
       kzdown,KDET,K22,KB,JMIN,kstabi,kstabm,K22x,jmin3,kdet3,         &   !-lxz
       KBCONx,KBx,KTOPx,ierr,ierr2,ierr3,KBMAX,ierr5,ierr5_0

     integer                              ::                           &
       nall,iedt,nens,nens3,ki,I,K,KK,iresult
     real                                 ::                           &
      day,dz,dzo,mbdt,mbdt_s,entr_rate,radius,entrd_rate,mentr_rate,mentrd_rate,  &
      zcutdown,edtmax,edtmin,depth_min,zkbmax,z_detr,zktop,      &
      massfld,dh,cap_maxs,trash,entr_rate3,mentr_rate3,frh,xlamdd
      real detdo1,detdo2,entdo,dp,subin,detdo,entup,                &
      detup,subdown,entdoj,entupk,detupk,totmas
      real :: power_entr

     integer :: jmini,levadj
     logical :: keep_going
     real xff_shal(9),blqe,xkshal
     character*50 :: ierrc(its:ite)
     real,    dimension (its:ite,kts:kte) ::                           &
       up_massentr,up_massdetr,dd_massentr,dd_massdetr                 &
      ,up_massentro,up_massdetro,dd_massentro,dd_massdetro
     real,    dimension (kts:kte) :: smth


     levadj=5
     power_entr=1.2
      day=86400.
      do i=its,itf
        xmb3(i)=0.
        closure_n(i)=16.
        xland1(i)=1.
        if(xland(i).gt.1.5)xland1(i)=0.
!       cap_max_increment(i)=50.
        cap_max_increment(i)=25.
      enddo
!
!--- specify entrainmentrate and detrainmentrate
!
      if(iens.le.4)then
      radius=14000.-float(iens)*2000.
      else
      radius=12000.
      endif
!
!--- gross entrainment rate (these may be changed later on in the
!--- program, depending what your detrainment is!!)
!
      entr_rate =.2/radius
      entr_rate=1.e-4
      radius=.2/entr_rate
      if(radius/dx .gt. 0.5)then
         radius=.5*dx
         entr_rate=.2/radius
      endif
      do i=its,itf
      !sig(i)=(1.-radius/dx)**2
      sig(i)=1.
      !write(11,*)'this run with dx,radius,entr_rate = ',dx,radius,entr_rate,sig(i)
      !write(0,*)'this run with dx,radius,entr_rate = ',i,its,itf,dx,radius,entr_rate,sig(i)
      enddo

      entr_rate3=.2/200.
!
!--- entrainment of mass
!
      mentrd_rate=0.
      xlamdd=mentrd_rate
      mentr_rate=entr_rate
      mentr_rate3=entr_rate3
!
!--- initial detrainmentrates
!
      do k=kts,ktf
      do i=its,itf
        cupclw(i,k)=0.
        cd(i,k)=1.*entr_rate
        cd3(i,k)=entr_rate3
        cdd(i,k)=0.
        zdo3(i,k)=0.
        hcdo(i,k)=0.
        qrcdo(i,k)=0.
        dellaqc(i,k)=0.
      enddo
      enddo
!
!--- max/min allowed value for epsilon (ratio downdraft base mass flux/updraft
!    base mass flux
!
      edtmax=1.
      edtmin=.1
!
!--- minimum depth (m), clouds must have
!
      depth_min=500.
!
!--- maximum depth (mb) of capping
!--- inversion (larger cap = no convection)
!
!      cap_maxs=125.
     cap_maxs=75.
      DO i=its,itf
        kbmax(i)=1
        jmin3(i)=0
        kdet3(i)=0
        aa0(i)=0.
        aa3_0(i)=0.
        aa1(i)=0.
        aa3(i)=0.
        aad(i)=0.
        edt(i)=0.
        edt3(i)=0.
        kstabm(i)=ktf-1
        IERR(i)=0
        IERR2(i)=0
        IERR3(i)=0
        IERR5(i)=0
        IERR5_0(i)=0
 enddo
!
!--- first check for upstream convection
!
!#if ( WRF_DFI_RADAR == 1 )
!  write(0,*)"WRFDFI_RADAR"
!  if(do_capsuppress == 1) then
!      do i=its,itf
!          cap_max(i)=cap_maxs
!          cap_max3(i)=25.
!          if(gsw(i,j).lt.1.or.high_resolution.eq.1)cap_max(i)=25.
!          if (abs(cap_suppress_j(i) - 1.0 ) < 0.1 ) then
!             cap_max(i)=cap_maxs+75.
!          elseif (abs(cap_suppress_j(i) - 0.0 ) < 0.1 ) then
!             cap_max(i)=10.0
!          endif
!          iresult=0
!      enddo
!  else
!     do i=its,itf
!         cap_max(i)=cap_maxs
!          cap_max3(i)=25.
!         if(gsw(i,j).lt.1.or.high_resolution.eq.1)cap_max(i)=25.
!       iresult=0
!     enddo
!  endif

!#else
!  write(0,*)"NO WRFDFI_RADAR",its,itf,jts,jtf,kts,ktf,ite,jte,kte
      do i=its,itf
          cap_max(i)=cap_maxs
          cap_max3(i)=25.
          if(gsw(i,j).lt.1.or.high_resolution.eq.1)cap_max(i)=25.
        iresult=0

      enddo
!#endif
!
!--- max height(m) above ground where updraft air can originate
!
      zkbmax=4000.
!
!--- height(m) above which no downdrafts are allowed to originate
!
      zcutdown=3000.
!
!--- depth(m) over which downdraft detrains all its mass
!
      z_detr=1250.
!
      do nens=1,maxens
         mbdt_ens(nens)=(float(nens)-3.)*dtime*1.e-3+dtime*5.E-03
      enddo
      do nens=1,maxens2
         edt_ens(nens)=.95-float(nens)*.01
      enddo
!
!--- environmental conditions, FIRST HEIGHTS
!
      do i=its,itf
         if(ierr(i).ne.20)then
            do k=1,maxens*maxens2*maxens3
               xf_ens(i,j,(iens-1)*maxens*maxens2*maxens3+k)=0.
               pr_ens(i,j,(iens-1)*maxens*maxens2*maxens3+k)=0.
            enddo
         endif
      enddo
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(z,qes,he,hes,t,q,p,z1, &
           psur,ierr,tcrit,0,xl,cp,   &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_env(zo,qeso,heo,heso,tn,qo,po,z1, &
           psur,ierr,tcrit,0,xl,cp,   &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)

!
!--- environmental values on cloud levels
!
      call cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,he_cup, &
           hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
           ierr,z1,xl,rv,cp,          &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_env_clev(tn,qeso,qo,heo,heso,zo,po,qeso_cup,qo_cup, &
           heo_cup,heso_cup,zo_cup,po_cup,gammao_cup,tn_cup,psur,  &
           ierr,z1,xl,rv,cp,          &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      do i=its,itf
        if(aaeq(i).lt.-0.1)then
           ierr(i)=20
        endif
!     if(ierr(i).eq.0)then
!
      do k=kts,ktf
        if(zo_cup(i,k).gt.zkbmax+z1(i))then
          kbmax(i)=k
          go to 25
        endif
      enddo
 25   continue
!
!--- level where detrainment for downdraft starts
!
      do k=kts,ktf
        if(zo_cup(i,k).gt.z_detr+z1(i))then
          kdet(i)=k
          go to 26
        endif
      enddo
 26   continue
!
!     endif
      enddo
      if(ipr.eq.1)then
        !write(0,*)'k,z(i,k),he(i,k),hes(i,k)'
      i=1
      do k=kts,ktf
        !write(0,*)k,z(i,k),he(i,k),hes(i,k)
      enddo
      endif

!
!
!
!------- DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22
!
      CALL cup_MAXIMI(HEO_CUP,3,KBMAX,K22,ierr, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
       !write(0,*)'initial k22 = ',k22(1)
       DO 36 i=its,itf
         IF(ierr(I).eq.0)THEN
         IF(K22(I).GE.KBMAX(i))ierr(i)=2
           ierrc(i)="could not find k22"
         endif
 36   CONTINUE
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!
      call cup_kbcon(ierrc,cap_max_increment,1,k22,kbcon,heo_cup,heso_cup, &
           ierr,kbmax,po_cup,cap_max, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!--- increase detrainment in stable layers
!
      CALL cup_minimi(HEso_cup,Kbcon,kstabm,kstabi,ierr,  &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      DO i=its,ite
         IF(ierr(I).eq.0)THEN
            do k=kts,ktf
   ! Saulo's approach
               frh = min(qo_cup(i,k)/qeso_cup(i,k),1.)
               entr_rate_2d(i,k)=1.e-4*(1.3-frh)
            enddo
            zuhe(i)=0.1
            frh=(1.-zuhe(i))/((float(kbcon(i))**power_entr)-float(k22(i))**power_entr)
            dh=zuhe(i)-frh*float(k22(i))**power_entr
            do k=k22(i),kbcon(i)-1
             dz=z_cup(i,k+1)-z_cup(i,k)
             entr_rate_2d(i,k)=((frh*(float((k+1))**power_entr)+dh)/zuhe(i)-1.)/dz
             zuhe(i)=zuhe(i)+entr_rate_2d(i,k)*dz*zuhe(i)
!             entr_rate_2d(i,k)=1./(0.2*dz)
!             entr_rate_2d(i,k)=(1.+.2**frh)/(dz)
            enddo
               entr_rate=entr_rate_2d(i,kbcon(i))
               mentr_rate=entr_rate
            do k=kts,ktf
               cd(i,k)=entr_rate_2d(i,kbcon(i))
               !write(0,*)'frh ',k,entr_rate_2d(i,k),entr_rate,cd(i,k)
            enddo
            if(kstabm(i)-1.gt.kstabi(i))then
               do k=kstabi(i),kstabm(i)-1
                 cd(i,k)=cd(i,k-1)+1.0*entr_rate
                 if(cd(i,k).gt.10.0*entr_rate)cd(i,k)=10.0*entr_rate
               !write(0,*)'frh2 ',k,entr_rate_2d(i,k),entr_rate,cd(i,k)
               enddo
               do k=kts+1,ktf-1
                 smth(k)=.25*(cd(i,k-1)+2.*cd(i,k)+cd(i,k+1))
               enddo
               do k=kts+1,ktf-1
                 cd(i,k)=smth(k)
               enddo
	     endif
               smth(:)=0.
               do k=k22(i)+0,ktf-1
                 smth(k)=.25*(entr_rate_2d(i,k-1)+2.*entr_rate_2d(i,k)+entr_rate_2d(i,k+1))
               enddo
               do k=k22(i)+0,kbcon(i)
                 entr_rate_2d(i,k)=smth(k)
               enddo
               zuhe(i)=.1
               do k=k22(i)+0,kbcon(i)+2
                 dz=z_cup(i,k+1)-z_cup(i,k)
                 frh=zuhe(i)
                 zuhe(i)=zuhe(i)+entr_rate_2d(i,k)*dz*zuhe(i)
                 !write(11,*)'zuhe,dz,frh = ',zuhe(i),frh,dz,entr_rate_2d(i,k)
                   if(zuhe(i).gt.1.)then
                     entr_rate_2d(i,k)=(1.-frh)/dz/frh
                     zuhe(i)=frh+entr_rate_2d(i,k)*dz*frh
                   endif
                 !write(11,*)'zuhe,dz,frh = ',entr_rate_2d(i,k)
               enddo
         ENDIF
       enddo
!     do i=its,itf
!     IF(ierr(I).eq.0)THEN
!       if(kstabm(i)-1.gt.kstabi(i))then
!          do k=kstabi(i),kstabm(i)-1
!            cd(i,k)=cd(i,k-1)+1.0*entr_rate
!            if(cd(i,k).gt.10.0*entr_rate)cd(i,k)=10.0*entr_rate
!          enddo
!       ENDIF
!     ENDIF
!     ENDDO
!
! calculate mass entrainment and detrainment
!
      do k=kts,ktf
      do i=its,itf
         hc(i,k)=0.
         DBY(I,K)=0.
         hco(i,k)=0.
         DBYo(I,K)=0.
      enddo
      enddo
      do i=its,itf
         hkb(i)=he_cup(i,k22(i))
         hkbo(i)=heo_cup(i,k22(i))
         do k=1,k22(i)
            hc(i,k)=he_cup(i,k)
            hco(i,k)=heo_cup(i,k)
         enddo
         do k=k22(i),kbcon(i)-1
            hc(i,k)=hkb(i)
            hco(i,k)=hkbo(i)
         enddo
          k=kbcon(i)
          hc(i,k)=hkb(i)
         DBY(I,Kbcon(i))=Hkb(I)-HES_cup(I,K)
          hco(i,k)=hkbo(i)
         DBYo(I,Kbcon(i))=Hkbo(I)-HESo_cup(I,K)
      enddo
!
!
      do i=its,itf
      if(ierr(i).eq.0)then
      do k=1,k22(i)-1
       zu(i,k)=0.
       zuo(i,k)=0.
      enddo
      zu(i,k22(i))=0.1
      zuo(i,k22(i))=0.1
!     zu(i,k22(i))=entr_rate_2d(i,k-1)*dz*zuo(i,k-1)
!     zuo(i,k22(i))=entr_rate_2d(i,k-1)*dz*zuo(i,k-1)
      do k=k22(i)+1,ktf
       dz=zo_cup(i,k)-zo_cup(i,k-1)
       up_massentro(i,k)=entr_rate_2d(i,k-1)*dz*zuo(i,k-1)
       up_massdetro(i,k)=cd(i,k-1)*dz*zuo(i,k-1)
       zuo(i,k)=zuo(i,k-1)+up_massentro(i,k)-up_massdetro(i,k)
       dz=z_cup(i,k)-z_cup(i,k-1)
       up_massentr(i,k)=entr_rate_2d(i,k-1)*dz*zu(i,k-1)
       up_massdetr(i,k)=cd(i,k-1)*dz*zu(i,k-1)
       zu(i,k)=zu(i,k-1)+up_massentr(i,k)-up_massdetr(i,k)
      enddo
      do k=kbcon(i)+1,ktf
       hc(i,k)=(hc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k)*hc(i,k-1)+ &
                         up_massentr(i,k)*he(i,k-1))   /            &
                         (zu(i,k-1)-.5*up_massdetr(i,k)+up_massentr(i,k))
       dby(i,k)=hc(i,k)-hes_cup(i,k)
       hco(i,k)=(hco(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k)*hco(i,k-1)+ &
                         up_massentro(i,k)*heo(i,k-1))   /            &
                         (zuo(i,k-1)-.5*up_massdetro(i,k)+up_massentro(i,k))
       dbyo(i,k)=hco(i,k)-heso_cup(i,k)
      enddo
      do k=kbcon(i)+1,ktf
       if(dbyo(i,k).lt.0)then
           ktop(i)=k-1
           go to 41
       endif
      enddo
41    continue
      if(ktop(i).lt.kbcon(i)+2)ierr(i)=5
      do k=1,k22(i)
           up_massentr(i,k)=0.
           up_massdetr(i,k)=0.
           up_massentro(i,k)=0.
           up_massdetro(i,k)=0.
      enddo
      do k=ktop(i)+1,ktf
           HC(i,K)=hes_cup(i,k)
           HCo(i,K)=heso_cup(i,k)
           DBY(I,K)=0.
           DBYo(I,K)=0.
           zu(i,k)=0.
           zuo(i,k)=0.
           up_massentr(i,k)=0.
           up_massdetr(i,k)=0.
           up_massentro(i,k)=0.
           up_massdetro(i,k)=0.
      enddo
        !write(0,*)'hcnew = '
      do k=1,ktf
        !write(0,*)k,hco(i,k),dbyo(i,k)
!       write(0,*)k,hc(i,k),.5*up_massdetr(i,k-1)*hc(i,k-1)/zu(i,k-1), &
!                  up_massentr(i,k-1)*he(i,k-1)/zu(i,k-1),(zu(i,k)+.5*up_massdetr(i,k-1))/zu(i,k-1)
      enddo
      endif
      enddo
!
!--- calculate incloud moist static energy
!
!     call cup_up_he(k22,hkb,z_cup,cd,mentr_rate,he_cup,hc, &
!          kbcon,ierr,dby,he,hes_cup,'deep', &
!          itf,jtf,ktf, &
!          its,ite, jts,jte, kts,kte)
!     call cup_up_he(k22,hkbo,zo_cup,cd,mentr_rate,heo_cup,hco, &
!          kbcon,ierr,dbyo,heo,heso_cup,'deep', &
!          itf,jtf,ktf, &
!          its,ite, jts,jte, kts,kte)
!     if(ipr.eq.1)then
!       write(0,*)"heso_cup(i,k),hco(i,k),dbyo(i,k)"
!      i=1
!      do k=kts,ktf
!       write(0,*)k,hco(i,k),dbyo(i,k)
!      enddo
!     endif

!--- DETERMINE CLOUD TOP - KTOP
!
!     call cup_ktop(ierrc,1,dbyo,kbcon,ktop,ierr, &
!          itf,jtf,ktf, &
!          its,ite, jts,jte, kts,kte)
      DO 37 i=its,itf
         kzdown(i)=0
         if(ierr(i).eq.0)then
            zktop=(zo_cup(i,ktop(i))-z1(i))*.6
            zktop=min(zktop+z1(i),zcutdown+z1(i))
            do k=kts,kte
              if(zo_cup(i,k).gt.zktop)then
                 kzdown(i)=k
                 go to 37
              endif
              enddo
         endif
 37   CONTINUE
!     do i=its,itf
!     do k=ktop(i),ktf
!          HC(i,K)=hes_cup(i,k)
!          HCo(i,K)=heso_cup(i,k)
!          DBY(I,K)=0.
!          DBYo(I,K)=0.
!     enddo
!     enddo
      do k=kts,ktf
        i=1
        !write(9,*)p_cup(i,k),heo_cup(i,k),heso_cup(i,k),hco(i,k)
      enddo
!
!--- DOWNDRAFT ORIGINATING LEVEL - JMIN
!
      call cup_minimi(HEso_cup,K22,kzdown,JMIN,ierr, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      DO 100 i=its,ite
      IF(ierr(I).eq.0)THEN
!
!--- check whether it would have buoyancy, if there where
!--- no entrainment/detrainment
!
      jmini = jmin(i)
      keep_going = .TRUE.
      do while ( keep_going )
        keep_going = .FALSE.
        if ( jmini - 1 .lt. kdet(i)   ) kdet(i) = jmini-1
        if ( jmini     .ge. ktop(i)-1 ) jmini = ktop(i) - 2
        ki = jmini
        hcdo(i,ki)=heso_cup(i,ki)
        DZ=Zo_cup(i,Ki+1)-Zo_cup(i,Ki)
        dh=0.
        do k=ki-1,1,-1
          hcdo(i,k)=heso_cup(i,jmini)
          DZ=Zo_cup(i,K+1)-Zo_cup(i,K)
          dh=dh+dz*(HCDo(i,K)-heso_cup(i,k))
          if(dh.gt.0.)then
            jmini=jmini-1
            if ( jmini .gt. 5 ) then
              keep_going = .TRUE.
            else
              ierr(i) = 9
              ierrc(i) = "could not find jmini9"
              exit
            endif
          endif
        enddo
      enddo
      jmin(i) = jmini
      if ( jmini .le. 5 ) then
        ierr(i)=4
        ierrc(i) = "could not find jmini4"
      endif
      ENDIF
100   continue
!
! - Must have at least depth_min m between cloud convective base
!     and cloud top.
!
      do i=its,itf
      IF(ierr(I).eq.0)THEN
      IF(-zo_cup(I,KBCON(I))+zo_cup(I,KTOP(I)).LT.depth_min)then
            ierr(i)=6
            ierrc(i)="cloud depth very shallow"
      endif
      endif
      enddo

!
!c--- normalized updraft mass flux profile
!
!
!
!
!c--- normalized downdraft mass flux profile,also work on bottom detrainment
!--- in this routine
!
      do k=kts,ktf
      do i=its,itf
       zd(i,k)=0.
       zdo(i,k)=0.
       cdd(i,k)=0.
       dd_massentr(i,k)=0.
       dd_massdetr(i,k)=0.
       dd_massentro(i,k)=0.
       dd_massdetro(i,k)=0.
       hcdo(i,k)=heso_cup(i,k)
       dbydo(i,k)=0.
      enddo
      enddo
      do i=its,itf
          IF(ierr(I).eq.0)then
            dz = z_cup(i,kbcon(i))/float(kbcon(i))
            frh = 1./float(kbcon(i))
            dzo = (1.-beta**frh)/dz
            do ki=jmin(i),1,-1
              if (ki.ge.kbcon(i)) then
                cdd(i,ki) = xlamdd
              else
                cdd(i,ki) = xlamdd + dzo
              endif
            end do
            smth(:)=0.
            do k=2,jmin(i)-1
              smth(k)=(cdd(i,k-1)+2.*cdd(i,k)+cdd(i,k+1))*.25
            enddo
            do k=2,jmin(i)-1
              cdd(i,k)=smth(k)
            enddo
            zd(i,jmin(i))=0.2
            zdo(i,jmin(i))=0.2
            frh=.8/(float(levadj)**power_entr-1.)
            dh=.2-frh
            zuhe(i)=.2
            mentrd_rate_2d(i,:)=mentrd_rate
            k=1
            do ki=jmin(i),jmin(i)-(levadj-2),-1
             k=k+1
             dz=z_cup(i,ki)-z_cup(i,ki-1)
             mentrd_rate_2d(i,ki)=((frh*(float((k))**power_entr)+dh)/zuhe(i)-1.)/dz
             zuhe(i)=zuhe(i)+mentrd_rate_2d(i,ki)*dz*zuhe(i)
            enddo
            smth(:)=0.
            do k=2,jmin(i)-1
              smth(k)=(mentrd_rate_2d(i,k-1)+2.*mentrd_rate_2d(i,k) &
                       +mentrd_rate_2d(i,k+1))*.25
            enddo
            do k=2,jmin(i)-1
              mentrd_rate_2d(i,k)=smth(k)
            enddo
            if(i.eq.ipr .and. j.eq.jpr)write(11,*)'k22,kbcon,jmin = ',k22(i),kbcon(i),jmin(i)
!
! now that we have entrainment and detrainment rates,
! calculate downdraft mass terms
           do ki=jmin(i)-1,1,-1
               mentrd_rate=mentrd_rate_2d(i,ki+1)
!               dz=z_cup(i,ki+1)-z_cup(i,ki)
!              dd_massentr(i,ki)=mentrd_rate*dz*zd(i,ki+1)
!              dd_massdetr(i,ki)=cdd(i,ki+1)*dz*zd(i,ki+1)
!              zd(i,ki)=zd(i,ki+1)+dd_massentr(i,ki)-dd_massdetr(i,ki)

               dzo=zo_cup(i,ki+1)-zo_cup(i,ki)
               dd_massentro(i,ki)=mentrd_rate*dzo*zdo(i,ki+1)
               dd_massdetro(i,ki)=cdd(i,ki+1)*dzo*zdo(i,ki+1)
               zdo(i,ki)=zdo(i,ki+1)+dd_massentro(i,ki)-dd_massdetro(i,ki)
            enddo
! downdraft moist static energy + moisture budget
!           hcd(i,jmin(i))=hes_cup(i,jmin(i))
!           dbyd(i,jmin(i))=hcd(i,jmin(i))-hes_cup(i,jmin(i))
            do ki=jmin(i)-1,1,-1
!            hcd(i,ki)=(hcd(i,ki+1)*zd(i,ki+1)                          &
!                       -.5*dd_massdetr(i,ki)*hcd(i,ki+1)+            &
!                        dd_massentr(i,ki)*he(i,ki+1))   /            &
!                      (zd(i,ki+1)-.5*dd_massdetr(i,ki)+dd_massentr(i,ki))
!            dbyd(i,ki)=hcd(i,ki)-hes_cup(i,ki)
             hcdo(i,ki)=(hcdo(i,ki+1)*zdo(i,ki+1)                       &
                         -.5*dd_massdetro(i,ki)*hcdo(i,ki+1)+ &
                        dd_massentro(i,ki)*heo(i,ki+1))   /            &
                        (zdo(i,ki+1)-.5*dd_massdetro(i,ki)+dd_massentro(i,ki))
             dbydo(i,ki)=hcdo(i,ki)-heso_cup(i,ki)
             enddo
          endif

      enddo
!        do k=kts,ktop(1)
!            write(0,*)'zdnew = ',k,zd(i,k),zdo(i,k),hcd(i,k)
!        enddo
!     call cup_dd_nms(zd,z_cup,cdd,mentrd_rate,jmin,ierr, &
!          0,kdet,z1,                 &
!             kbcon,beta,xlamdd,                      &
!          itf,jtf,ktf, &
!          its,ite, jts,jte, kts,kte)
!     call cup_dd_nms(zdo,zo_cup,cdd,mentrd_rate,jmin,ierr, &
!          1,kdet,z1,                 &
!             kbcon,beta,xlamdd,                      &
!          itf,jtf,ktf, &
!          its,ite, jts,jte, kts,kte)
!        do k=kts,ktop(1)
!            write(0,*)'zdold = ',k,zd(i,k),zdo(i,k),hcd(i,k)
!        enddo
      !i=1
      !do k=kts+1,ktop(1)
         !write(13,*)po_cup(1,k),zuo(1,k),up_massentro(i,k),up_massdetro(i,k)
      !enddo
      !do k=jmin(1)-1,kts,-1
         !write(14,*)po_cup(1,k),zdo(1,k),dd_massentro(i,k),dd_massdetro(i,k)
      !enddo
!
!--- downdraft moist static energy
!
!     call cup_dd_he(hes_cup,zd,hcd,z_cup,cdd,mentrd_rate, &
!          jmin,ierr,he,dbyd,he_cup,  &
!          itf,jtf,ktf, &
!          its,ite, jts,jte, kts,kte)
!     call cup_dd_he(heso_cup,zdo,hcdo,zo_cup,cdd,mentrd_rate, &
!          jmin,ierr,heo,dbydo,he_cup,&
!          itf,jtf,ktf, &
!          its,ite, jts,jte, kts,kte)
!
!--- calculate moisture properties of downdraft
!
!      call cup_dd_moisture_new(ierrc,zd,hcd,hes_cup,qcd,qes_cup, &
!           pwd,q_cup,z_cup,dd_massentr,dd_massdetr,jmin,ierr,gamma_cup, &
!           pwev,bu,qrcd,q,he,t_cup,2,xl,high_resolution, &
!           itf,jtf,ktf, &
!           its,ite, jts,jte, kts,kte)
      call cup_dd_moisture_new(ierrc,zdo,hcdo,heso_cup,qcdo,qeso_cup, &
           pwdo,qo_cup,zo_cup,dd_massentro,dd_massdetro,jmin,ierr,gammao_cup, &
           pwevo,bu,qrcdo,qo,heo,tn_cup,1,xl,high_resolution, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!--- calculate moisture properties of updraft
!
!     call cup_up_moisture('deep',ierr,z_cup,qc,qrc,pw,pwav, &
!          ccnclean,p_cup,kbcon,ktop,cd,dby,mentr_rate,clw_all,      &
!          t_cup,q,GAMMA_cup,zu,qes_cup,k22,q_cup,xl,ccn,rho,psum,psumh, &
!          autoconv,aeroevap,1,itf,jtf,ktf, &
!          its,ite, jts,jte, kts,kte)
     call cup_up_moisture('deep',ierr,zo_cup,qco,qrco,pwo,pwavo, &
           ccnclean,p_cup,kbcon,ktop,cd,dbyo,mentr_rate,clw_all, &
           t_cup,qo,GAMMAo_cup,zuo,qeso_cup,k22,qo_cup,xl,        &
           ccn,rho,up_massentr,up_massdetr,psum,psumh,&
           autoconv,aeroevap,1,itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)



      do k=kts,ktf
      do i=its,itf
         cupclw(i,k)=qrco(i,k)
      enddo
      enddo
!
!--- calculate workfunctions for updrafts
!
      call cup_up_aa0(aa0,z,zu,dby,GAMMA_CUP,t_cup, &
           kbcon,ktop,ierr,           &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_up_aa0(aa1,zo,zuo,dbyo,GAMMAo_CUP,tn_cup, &
           kbcon,ktop,ierr,           &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      do i=its,itf
         if(ierr(i).eq.0)then
           if(aa1(i).eq.0.)then
               ierr(i)=17
               ierrc(i)="cloud work function zero"
           endif
         endif
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    NEXT section for shallow convection
!
      if(ishallow_g3.eq.1)then
!     write(0,*)'now do shallow for j = ',j
      call cup_env(z3,qes3,he3,hes3,tshall,qshall,po,z1, &
           psur,ierr5,tcrit,0,xl,cp,   &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_env_clev(tshall,qes3,qshall,he3,hes3,z3,po,qes3_cup,q3_cup, &
           he3_cup,hes3_cup,z3_cup,po_cup,gamma3_cup,t3_cup,psur,  &
           ierr5,z1,xl,rv,cp,          &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      CALL cup_MAXIMI(HE3_CUP,1,kbmax,K23,ierr5, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
       DO i=its,itf
         if(kpbl(i).gt.5)cap_max3(i)=po_cup(i,kpbl(i))
         IF(ierr5(I).eq.0.)THEN
         IF(K23(I).Gt.Kbmax(i))ierr5(i)=2
         if(kpbl(i).gt.5)k23(i)=kpbl(i)
         endif
         ierr5_0(i)=ierr5(i)
       ENDDO
      call cup_kbcon(ierrc,cap_max_increment,5,k23,kbcon3,he3_cup,hes3_cup, &
           ierr5,kbmax,po_cup,cap_max3, &
!          ierr5,kpbl,po_cup,cap_max3, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_up_he(k23,hkb3,z3_cup,cd3,mentr_rate3,he3_cup,hc3, &
           kbcon3,ierr5,dby3,he3,hes3_cup,'shallow', &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_up_he(k23,hkb3_0,z_cup,cd3,mentr_rate3,he_cup,hc3_0, &
           kbcon3,ierr5,dby3_0,he,hes_cup,'shallow', &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_ktop(ierrc,1,dby3,kbcon3,ktop3,ierr5, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_up_nms(0,zu3,z3_cup,mentr_rate3,cd3,kbcon3,ktop3,    &
           ierr5,k23, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_up_nms(0,zu3_0,z_cup,mentr_rate3,cd3,kbcon3,ktop3,    &
           ierr5,k23, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
! first calculate aa3_0_cup
!
      call cup_up_aa0(aa3_0,z,zu3_0,dby3_0,GAMMA3_CUP,t_cup, &
           kbcon3,ktop3,ierr5,           &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!  now what is necessary for aa3 and feedbacks
!
      call cup3_up_moisture('shallow',ierr5,z3_cup,qc3,qrc3,pw3,pwav3, &
           kbcon3,ktop3,cd3,dby3,mentr_rate3,clw_all, &
           qshall,GAMMA3_cup,zu3,qes3_cup,k23,q3_cup,xl,&
           itf,jtf,ktf, &
	   its,ite, jts,jte, kts,kte)
      call cup_up_aa0(aa3,z3,zu3,dby3,GAMMA3_CUP,t3_cup, &
           kbcon3,ktop3,ierr5,           &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!     do i=its,itf
!        if(ierr5(i).eq.0)then
!          if(aa3(i).eq.0.)then
!              ierr5(i)=17
!          endif
!        endif
!     enddo
!     call cup_dellabot('shallow',ipr,jpr,q3_cup,ierr5,z3_cup,po,qrcdo,edto, &
!          zdo,cdd,q3,dellaq3,dsubq,j,mentrd_rate,z3,g,&
!          itf,jtf,ktf, &
!          its,ite, jts,jte, kts,kte)
      call cup_dellas_3d(ierr5,z3_cup,po_cup,hcdo,edt3,zdo3,cdd,    &
           he3,dellah3,dsubt3,j,mentrd_rate,zu3,g,                     &
           cd3,hc3,ktop3,k23,kbcon3,mentr_rate3,jmin,he3_cup,kdet, &
           k23,ipr,jpr,'shallow',0,                                 &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_dellas_3d(ierr5,z3_cup,po_cup,qrcdo,edt3,zdo3,cdd, &
           qshall,dellaq3,dsubq3,j,mentrd_rate,zu3,g, &
           cd3,qc3,ktop3,k23,kbcon3,mentr_rate3,jmin,q3_cup,kdet, &
           k23,ipr,jpr,'shallow',0,               &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte    )
              mbdt_s=1.e-1*mbdt_ens(1)
              do k=kts,ktf
              do i=its,itf
                 dellat3(i,k)=0.
                 if(ierr5(i).eq.0)then
                    trash=dsubt3(i,k)
                    XHE3(I,K)=(dsubt3(i,k)+DELLAH3(I,K))*MBDT_S+HE3(I,K)
                    XQ3(I,K)=(dsubq3(i,k)+DELLAQ3(I,K))*MBDT_S+QSHALL(I,K)
                    DELLAT3(I,K)=(1./cp)*(DELLAH3(I,K)-xl*DELLAQ3(I,K))
                    dSUBT3(I,K)=(1./cp)*(dsubt3(i,k)-xl*dsubq3(i,k))
                    XT3(I,K)= (DELLAT3(I,K)+dsubt3(i,k))*MBDT_S+TSHALL(I,K)
                    IF(XQ3(I,K).LE.0.)XQ3(I,K)=1.E-08
!                    if(i.eq.ipr.and.j.eq.jpr)then
!                      write(0,*)k,trash,DELLAQ3(I,K),dsubq3(I,K),dsubt3(i,k)
!                    endif
                 ENDIF
              enddo
              enddo
      do i=its,itf
      if(ierr5(i).eq.0)then
      XHE3(I,ktf)=HE3(I,ktf)
      XQ3(I,ktf)=QSHALL(I,ktf)
      XT3(I,ktf)=TSHALL(I,ktf)
      IF(XQ3(I,ktf).LE.0.)XQ3(I,ktf)=1.E-08
      endif
      enddo
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(xz3,xqes3,xhe3,xhes3,xt3,xq3,po,z1, &
           psur,ierr5,tcrit,2,xl,cp,   &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!--- environmental values on cloud levels
!
      call cup_env_clev(xt3,xqes3,xq3,xhe3,xhes3,xz3,po,xqes3_cup,xq3_cup, &
           xhe3_cup,xhes3_cup,xz3_cup,po_cup,gamma3_cup,xt3_cup,psur,   &
           ierr5,z1,xl,rv,cp,          &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!
!**************************** static control
!
!--- moist static energy inside cloud
!
      do i=its,itf
        if(ierr5(i).eq.0)then
          xhkb3(i)=xhe3(i,k23(i))
        endif
      enddo
      call cup_up_he(k23,xhkb3,xz3_cup,cd3,mentr_rate3,xhe3_cup,xhc3, &
           kbcon3,ierr5,xdby3,xhe3,xhes3_cup,'shallow', &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!c--- normalized mass flux profile and CWF
!
      call cup_up_nms(0,xzu3,xz3_cup,mentr_rate3,cd3,kbcon3,ktop3,ierr5,k23, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_up_aa0(xaa3,xz3,xzu3,xdby3,GAMMA3_CUP,xt3_cup, &
           kbcon3,ktop3,ierr5,           &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
! now for shallow forcing
!
       do i=its,itf
        xmb3(i)=0.
        xff_shal(1:9)=0.
        if(ierr5(i).eq.0)then
          xkshal=(xaa3(i)-aa3(i))/mbdt_s
          if(xkshal.ge.0.)xkshal=+1.e6
          if(xkshal.gt.-1.e-4 .and. xkshal.lt.0.)xkshal=-1.e-4
          xff_shal(1)=max(0.,-(aa3(i)-aa3_0(i))/(xkshal*dtime))
          xff_shal(2)=max(0.,-(aa3(i)-aa3_0(i))/(xkshal*dtime))
          xff_shal(3)=max(0.,-(aa3(i)-aa3_0(i))/(xkshal*dtime))
          if(aa3_0(i).le.0)then
           xff_shal(1)=0.
           xff_shal(2)=0.
           xff_shal(3)=0.
          endif
          if(aa3(i)-aa3_0(i).le.0.)then
           xff_shal(1)=0.
           xff_shal(2)=0.
           xff_shal(3)=0.
          endif
! boundary layer QE (from Saulo Freitas)
          blqe=0.
          trash=0.
          if(k23(i).lt.kpbl(i)+1)then
             do k=1,kbcon3(i)-1
                blqe=blqe+100.*dhdt(i,k)*(p_cup(i,k)-p_cup(i,k+1))/g
             enddo
             trash=max((hc3(i,kbcon3(i))-he_cup(i,kbcon3(i))),1.e1)
             xff_shal(7)=max(0.,blqe/trash)
             xff_shal(7)=min(0.1,xff_shal(7))
          else
             xff_shal(7)=0.
          endif
          if((xkshal.lt.-1.1e-04) .and.  &
             ((aa3(i)-aa3_0(i).gt.0.) .or. (xff_shal(7).gt.0)))then
          xff_shal(4)=max(0.,-aa3(i)/(xkshal*tscl_KF))
          xff_shal(4)=min(0.1,xff_shal(4))
          xff_shal(5)=xff_shal(4)
          xff_shal(6)=xff_shal(4)
          else
           xff_shal(4)=0.
           xff_shal(5)=0.
           xff_shal(6)=0.
          endif
!         write(0,888)'i0=',i,j,kpbl(i),blqe,xff_shal(7)
888       format(a3,3(1x,i3),2e12.4)
          xff_shal(8)= xff_shal(7)
          xff_shal(9)= xff_shal(7)
          do k=1,9
           xmb3(i)=xmb3(i)+xff_shal(k)
          enddo
          xmb3(i)=min(.1,xmb3(i)/9.)
!         if(xmb3(i).eq.10.1 )then
!           write(0,*)'i0,xmb3,blqe,xkshal = ',i,j,xmb3(i),blqe,xkshal
!           if(xff_shal(7).ge.0.1)then
!             write(0,*)'i1,blqe,trash = ',blqe,trash
!           endif
!           if(xff_shal(7).eq.0 .and. xff_shal(1).ge.0.1)then
!              write(0,*)'i2,aa3_0(i),aa3(i),xaa3(i) = ',aa3_0(i),aa3(i),xaa3(i)
!           endif
!           if(xff_shal(5).ge.0.1)then
!              write(0,*)'i3,aa3(i),a0,xkshal= ',aa3(i),aa3_0(i),xkshal
!           endif
!           write(0,*)'i0, xff_shallow = ',xff_shal
!         endif
!!         if(xff_shal(7).eq.0 .and. xff_shal(4).gt.0 .and. xmb3(i).eq.0.5)then
!!           write(0,*)'i4,xmb3 = ',i,j,xmb3(i),xkshal
!!           write(0,*)'xff_shallow = ',xff_shal
!!           write(0,*)aa3(i),xaa3(i),blqe
!!         endif
          if(xmb3(i).eq.0.)ierr5(i)=22
          if(xmb3(i).lt.0.)then
             ierr5(i)=21
!            write(0,*)'neg xmb,i,j,xmb3 for shallow = ',i,j,k23(i),ktop3(i),kbcon3(i),kpbl(i)
          endif
        endif
!         if(ierr5(i).eq.0)write(0,*)'i,j,xmb3 for shallow = ',i,j,xmb3(i),k23(i),ktop3(i)
!         if(ierr5(i).eq.0.and.i.eq.12.and.j.eq.25)write(0,*)'i,j,xmb3 for shallow = ',k23(i),ktop3(i),kbcon3(i),kpbl(i)
!         if(ierr5(i).eq.0)write(0,*)'i,j,xmb3 for shallow = ',i,j,k23(i),ktop3(i),kbcon3(i),kpbl(i)
        if(ierr5(i).ne.0)then
           k23(i)=0
           kbcon3(i)=0
           ktop3(i)=0
           xmb3(i)=0
           do k=kts,ktf
              outts(i,k)=0.
              outqs(i,k)=0.
           enddo
        else if(ierr5(i).eq.0)then
!
! got the mass flux, sanity check, first for heating rates
!
          trash=0.
          do k=2,ktop3(i)
           trash=max(trash,86400.*(dsubt3(i,k)+dellat3(i,k))*xmb3(i))
          enddo
          if(trash.gt.150.)xmb3(i)=xmb3(i)*150./trash
!
! sanity check on moisture tendencies: do not allow anything that may allow neg tendencies
!
          do k=2,ktop3(i)
           trash=q(i,k)+(dsubq3(i,k)+dellaq3(i,k))*xmb3(i)*dtime
          if(trash.lt.1.e-12)then
! max allowable tendency over tendency that would lead to too small mix ratios
!
            trash=((1.e-12-q(i,k))/dtime)                   &
                  /((dsubq3(i,k)+dellaq3(i,k))*xmb3(i))
            trash=max(0.,trash)
            trash=min(1.,trash)
            xmb3(i)=trash*xmb3(i)
          endif
          enddo
!
! final tendencies
!
          do k=2,ktop3(i)
           outts(i,k)=(dsubt3(i,k)+dellat3(i,k))*xmb3(i)
           outqs(i,k)=(dsubq3(i,k)+dellaq3(i,k))*xmb3(i)
          enddo
        endif
       enddo
!       if(j.eq.-25)then
!!        write(0,*)'!!!!!!!! j = ',j,' !!!!!!!!!!!!!!!!!!!!'
        i=12
!        write(0,*)k23(i),kbcon3(i),ktop3(i)
!        write(0,*)kpbl(i),ierr5(i),ierr(i)
!        write(0,*)xmb3(i),xff_shal(1:9)
!        write(0,*)xaa3(i),aa1(i),aa0(i),aa3(i)
!        do k=1,ktf
!          write(0,*)po(i,k),he3(i,k),hes3(i,k),dellah3(i,k)
!        enddo
!        do k=1,ktf
!          write(0,*)zu3(i,k),hc3(i,k),dsubt3(i,k),dellat3(i,k)
!        enddo
!        do k=1,ktop3(i)+1
!          blqe=cp*outts(i,k)+xl*outqs(i,k)
!          write(0,*)outts(i,k),outqs(i,k),blqe
!        enddo
!       endif
!
! done shallow
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!     call cup_axx(ierrc,tcrit,kbmax,z1,p,psur,xl,rv,cp,tx,qx,axx,ierr,    &
!	  cap_max,cap_max_increment,entr_rate,mentr_rate,&
!	  j,itf,jtf,ktf, &
!	  its,ite, jts,jte, kts,kte,ens4)

!2012 - tmp
      do i=1,ens4
       axx(:,i)=aa1(:)
      enddo
!2012 - tmp

!
!--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
!
      call cup_dd_edt(ierr,us,vs,zo,ktop,kbcon,edt,po,pwavo, &
           pwo,ccn,pwevo,edtmax,edtmin,maxens2,edtc,psum,psumh, &
           ccnclean,rho,aeroevap,itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      do 250 iedt=1,maxens2
        do i=its,itf
         if(ierr(i).eq.0)then
         edt(i)=edtc(i,iedt)
         edto(i)=edtc(i,iedt)
         edtx(i)=edtc(i,iedt)
         if(high_resolution.eq.1)then
            edt(i)=edtc(i,3)
            edto(i)=edtc(i,3)
            edtx(i)=edtc(i,3)
         endif
         endif
        enddo
        do k=kts,ktf
        do i=its,itf
           subt_ens(i,k,iedt)=0.
           subq_ens(i,k,iedt)=0.
           dellat_ens(i,k,iedt)=0.
           dellaq_ens(i,k,iedt)=0.
           dellaqc_ens(i,k,iedt)=0.
           pwo_ens(i,k,iedt)=0.
        enddo
        enddo
!
!      if(j.eq.jpr.and.iedt.eq.1.and.ipr.gt.its.and.ipr.lt.ite)then
!!      if(j.eq.jpr)then
!         i=ipr
!!        write(0,*)'in 250 loop ',iedt,edt(ipr),ierr(ipr)
!!       if(ierr(i).eq.0.or.ierr(i).eq.3)then
!         write(0,*)'250',k22(I),kbcon(i),ktop(i),jmin(i)
!         write(0,*)edt(i),aa0(i),aa1(i)
!         do k=kts,ktf
!           write(0,*)k,z(i,k),he(i,k),hes(i,k)
!         enddo
!         write(0,*)'end 250 loop ',iedt,edt(ipr),ierr(ipr)
          !i=1
          !do k=1,ktop(i)+1
            !write(0,*)zu(i,k),zd(i,k),pwo(i,k),pwdo(i,k)
          !enddo
!!        endif
!      endif
      do i=its,itf
        aad(i)=0.
      enddo
!
!--- change per unit mass that a model cloud would modify the environment
!
!--- 1. in bottom layer
!
      do k=kts,ktf
      do i=its,itf
        dellah(i,k)=0.
        dsubt(i,k)=0.
        dellaq(i,k)=0.
        dsubq(i,k)=0.
      enddo
      enddo

      do i=its,itf
        if(ierr(i).eq.0)then
         dp=100.*(po_cup(i,1)-po_cup(i,2))
         dellah(i,1)=(edto(i)*zdo(i,2)*hcdo(i,2)   &
                     -edto(i)*zdo(i,2)*heo_cup(i,2))*g/dp
         dellaq(i,1)=(edto(i)*zdo(i,2)*qrcdo(i,2)   &
                     -edto(i)*zdo(i,2)*qo_cup(i,2))*g/dp
         do k=kts+1,ktop(i)
! these three are only used at or near mass detrainment and/or entrainment levels
            entupk=0.
            detupk=0.
            entdoj=0.
! detrainment and entrainment for fowndrafts
            detdo=edto(i)*dd_massdetro(i,k)
            entdo=edto(i)*dd_massentro(i,k)
! entrainment/detrainment for updraft
            entup=up_massentro(i,k+1)
            detup=up_massdetro(i,k+1)
! subsidence by downdrafts only
            subin=-zdo(i,k+1)*edto(i)
            subdown=-zdo(i,k)*edto(i)
!
!         SPECIAL LEVELS
!
! updraft originates at k22, only updraft term at k22-1 is (zu-entupk)*he
            if(k.eq.k22(i)-1)then
               entupk=zuo(i,k+1)
!              entup=0.
!              detup=0.
            endif
! downdraft originating level, similiar to k22-1 for updraft
            if(k.eq.jmin(i))then
               entdoj=edto(i)*zdo(i,k)
!              subin=0.
!              detdo=0.
!              entdo=0.
            endif
            if(k.eq.ktop(i))then
               detupk=zuo(i,ktop(i))
               subin=0.
               subdown=0.
               detdo=0.
               entdo=0.
               entup=0.
               detup=0.
            endif
!           if(k.lt.kbcon(i))then
!              detup=0.
!              entup=0.
!           endif
            totmas=subin-subdown+detup-entup-entdo+ &
             detdo-entupk-entdoj+detupk+zuo(i,k+1)-zuo(i,k)
!               print *,'*********************',k,totmas
!              write(0,123)k,subin+zuo(i,k+1),subdown-zuo(i,k),detup,entup, &
!                          detdo,entdo,entupk,detupk
!             write(8,*)'totmas = ',k,totmas
            if(abs(totmas).gt.1.e-6)then
               !write(0,*)'*********************',i,j,k,totmas
!              print *,jmin(i),k22(i),kbcon(i),ktop(i)
               !write(0,123)k,subin,subdown,detup,entup, &
               !            detdo,entdo,entupk,detupk
123     formAT(1X,i2,8E12.4)
!        call wrf_error_fatal ( 'totmas .gt.1.e-6' )
            endif
            dp=100.*(po_cup(i,k-1)-po_cup(i,k))
            dellah(i,k)=(detup*.5*(HCo(i,K+1)+HCo(i,K)) &
                    +detdo*.5*(HCDo(i,K+1)+HCDo(i,K)) &
                    -entup*heo(i,k) &
                    -entdo*heo(i,k) &
                    +subin*heo_cup(i,k+1) &
                    -subdown*heo_cup(i,k) &
                    +detupk*(hco(i,ktop(i))-heo_cup(i,ktop(i)))    &
                    -entupk*heo_cup(i,k22(i)) &
                    -entdoj*heo_cup(i,jmin(i)) &
                     )*g/dp
!-srf changed on 27/feb/2013
            dellaq(i,k)=(detup*.5*(qco(i,K+1)+qco(i,K)-qrco(i,k+1)-qrco(i,k)) &
                    +detdo*.5*(qrcdo(i,K+1)+qrcdo(i,K)) &
                    -entup*qo(i,k) &
                    -entdo*qo(i,k) &
                    +subin*qo_cup(i,k+1) &
                    -subdown*qo_cup(i,k) &
                    +detupk*(qco(i,ktop(i))-qrco(i,ktop(i))-qo_cup(i,ktop(i)))    &
                    -entupk*qo_cup(i,k22(i)) &
                    -entdoj*qo_cup(i,jmin(i)) &
                     )*g/dp
!
!            dellaq(i,k)=(detup*.5*(qco(i,K+1)+qco(i,K)) &
!                    +detdo*.5*(qrcdo(i,K+1)+qrcdo(i,K)) &
!                    -entup*qo(i,k) &
!                    -entdo*qo(i,k) &
!                    +subin*qo_cup(i,k+1) &
!                    -subdown*qo_cup(i,k) &
!                    +detupk*(qco(i,ktop(i))-qo_cup(i,ktop(i)))    &
!                    -entupk*qo_cup(i,k22(i)) &
!                    -entdoj*qo_cup(i,jmin(i)) &
!                     )*g/dp
!
            if(high_resolution.eq.1)then
! the first term includes entr and detr into/from updraft as well as (entup-detup)*he(i,k) from
!  neighbouring point, to make things mass consistent....
!            if(k.ge.k22(i))then
                dellah(i,k)=(                          &
                    detup*.5*(HCo(i,K+1)+HCo(i,K))-entup*heo(i,k)+(entup-detup)*heo(i,k) &
                    +detdo*.5*(HCDo(i,K+1)+HCDo(i,K)) &
                    -entdo*heo(i,k) &
                    +subin*heo_cup(i,k+1) &
                    -subdown*heo_cup(i,k) &
                    +detupk*(hco(i,ktop(i))-heo(i,ktop(i)))    &
                    -entdoj*heo_cup(i,jmin(i)) &
                    -entupk*heo_cup(i,k22(i))+entupk*heo(i,k) &
                     )*g/dp
                dellaq(i,k)=(                          &
                    detup*.5*(qCo(i,K+1)+qCo(i,K))-entup*qo(i,k)+(entup-detup)*qo(i,k) &
                    +detdo*.5*(qrcdo(i,K+1)+qrcdo(i,K)) &
                    -entdo*qo(i,k) &
                    +subin*qo_cup(i,k+1) &
                    -subdown*qo_cup(i,k) &
                    +detupk*(qco(i,ktop(i))-qo(i,ktop(i)))    &
                    -entdoj*qo_cup(i,jmin(i)) &
                    -entupk*qo_cup(i,k22(i))+entupk*qo(i,k) &
                     )*g/dp
           endif
!
! updraft subsidence only
!
           if(k.ge.k22(i)-1.and.k.lt.ktop(i))then
             dsubt(i,k)=(zuo(i,k+1)*heo_cup(i,k+1) &
                    -zuo(i,k)*heo_cup(i,k))*g/dp
             dsubq(i,k)=(zuo(i,k+1)*qo_cup(i,k+1) &
                    -zuo(i,k)*qo_cup(i,k))*g/dp
           endif
!
! in igh res case, subsidence terms are for meighbouring points only. This has to be
! done mass consistent with the della term
         if(high_resolution.eq.1)then
            if(k.ge.k22(i).and.k.lt.ktop(i))then
               dsubt(i,k)=(zuo(i,k+1)*heo_cup(i,k+1)-zuo(i,k)*heo_cup(i,k)-(entup-detup)*heo(i,k))*g/dp
               dsubq(i,k)=(zuo(i,k+1)*qo_cup(i,k+1)-zuo(i,k)*qo_cup(i,k)-(entup-detup)*qo(i,k))*g/dp
            else if(k.eq.ktop(i))then
               dsubt(i,k)=detupk*(heo(i,ktop(i))-heo_cup(i,ktop(i)))*g/dp
               dsubq(i,k)=detupk*(qo(i,ktop(i))-qo_cup(i,ktop(i)))*g/dp
            else if(k.eq.k22(i)-1)then
               dsubt(i,k)=(entupk*heo(i,k)-entupk*heo_cup(i,k))*g/dp
               dsubq(i,k)=(entupk*qo(i,k)-entupk*qo_cup(i,k))*g/dp
            endif
         endif
       enddo   ! k

        endif
      enddo
!
!-- take out cloud liquid water for detrainment
!
      do k=kts,ktf-1
      do i=its,itf
       dellaqc(i,k)=0.
       if(ierr(i).eq.0)then
         if(k.eq.ktop(i)-0)dellaqc(i,k)= &
                      .01*zuo(i,ktop(i))*qrco(i,ktop(i))* &
                      9.81/(po_cup(i,k-1)-po_cup(i,k))
         if(k.lt.ktop(i).and.k.gt.kbcon(i))then
           dz=zo_cup(i,k+1)-zo_cup(i,k)
           dellaqc(i,k)=.01*9.81*cd(i,k+1)*dz*zuo(i,k) &
                        *.5*(qrco(i,k)+qrco(i,k+1))/ &
                        (po_cup(i,k-1)-po_cup(i,k))
         endif
       endif
       !write(11,*)'dellaqc = ',k,po_cup(i,k),dellaqc(i,k)
      enddo
      enddo
!
!--- using dellas, calculate changed environmental profiles
!
!     do 200 nens=1,maxens
      mbdt=mbdt_ens(2)
      do i=its,itf
      xaa0_ens(i,1)=0.
      xaa0_ens(i,2)=0.
      xaa0_ens(i,3)=0.
      enddo

!      if(j.eq.jpr)then
!               write(0,*)'xt',xl,'DELLAH(I,K),DELLAQ(I,K),dsubq(I,K),dsubt(i,k)'
!      endif
      !write(0,*)'DELLAS'
      do k=kts,ktf
      do i=its,itf
         dellat(i,k)=0.
         if(ierr(i).eq.0)then
            trash=dsubt(i,k)
            XHE(I,K)=(dsubt(i,k)+DELLAH(I,K))*MBDT+HEO(I,K)
            XQ(I,K)=(dsubq(i,k)+DELLAQ(I,K))*MBDT+QO(I,K)
            DELLAT(I,K)=(1./cp)*(DELLAH(I,K)-xl*DELLAQ(I,K))
            dSUBT(I,K)=(1./cp)*(dsubt(i,k)-xl*dsubq(i,k))
            XT(I,K)= (DELLAT(I,K)+dsubt(i,k))*MBDT+TN(I,K)
            IF(XQ(I,K).LE.0.)XQ(I,K)=1.E-08
!             if(i.eq.ipr.and.j.eq.jpr)then
                !write(0,*)k,trash+dellah(i,k),DELLAQ(I,K)+dsubq(I,K),dellat(i,k)+dsubt(i,k)
!             endif
         ENDIF
      enddo
      enddo
      do i=its,itf
      if(ierr(i).eq.0)then
      XHE(I,ktf)=HEO(I,ktf)
      XQ(I,ktf)=QO(I,ktf)
      XT(I,ktf)=TN(I,ktf)
      IF(XQ(I,ktf).LE.0.)XQ(I,ktf)=1.E-08
      endif
      enddo
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(xz,xqes,xhe,xhes,xt,xq,po,z1, &
           psur,ierr,tcrit,2,xl,cp,   &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!--- environmental values on cloud levels
!
      call cup_env_clev(xt,xqes,xq,xhe,xhes,xz,po,xqes_cup,xq_cup, &
           xhe_cup,xhes_cup,xz_cup,po_cup,gamma_cup,xt_cup,psur,   &
           ierr,z1,xl,rv,cp,          &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!
!**************************** static control
!
!--- moist static energy inside cloud
!
      do i=its,itf
        if(ierr(i).eq.0)then
          xhkb(i)=xhe(i,k22(i))
        endif
      enddo
      do k=kts,ktf
      do i=its,itf
         xhc(i,k)=0.
         xDBY(I,K)=0.
      enddo
      enddo
      do i=its,itf
        if(ierr(i).eq.0)then
         xhkb(i)=xhe_cup(i,k22(i))
         do k=1,k22(i)
            xhc(i,k)=xhe_cup(i,k)
         enddo
         do k=k22(i),kbcon(i)-1
            xhc(i,k)=xhkb(i)
         enddo
          k=kbcon(i)
          xhc(i,k)=xhkb(i)
          xDBY(I,Kbcon(i))=xHkb(I)-xHES_cup(I,K)
        endif
      enddo
!
!
      do i=its,itf
      if(ierr(i).eq.0)then
      do k=1,k22(i)-1
       xzu(i,k)=0.
      enddo
      xzu(i,k22(i))=1.
      do k=k22(i)+1,ktf
       dz=xz_cup(i,k)-xz_cup(i,k-1)
       up_massentr(i,k)=entr_rate_2d(i,k-1)*dz*xzu(i,k-1)
       up_massdetr(i,k)=cd(i,k-1)*dz*xzu(i,k-1)
       xzu(i,k)=xzu(i,k-1)+up_massentr(i,k)-up_massdetr(i,k)
      enddo
      do k=kbcon(i)+1,ktf
       xhc(i,k)=(xhc(i,k-1)*xzu(i,k-1)-.5*up_massdetr(i,k)*xhc(i,k-1)+ &
                         up_massentr(i,k)*xhe(i,k-1))   /            &
                         (xzu(i,k-1)-.5*up_massdetr(i,k)+up_massentr(i,k))
       xdby(i,k)=xhc(i,k)-xhes_cup(i,k)
      enddo
      do k=ktop(i)+1,ktf
           xHC(i,K)=xhes_cup(i,k)
           xDBY(I,K)=0.
           xzu(i,k)=0.
      enddo
      endif
      enddo

!     call cup_up_he(k22,xhkb,xz_cup,cd,mentr_rate,xhe_cup,xhc, &
!          kbcon,ierr,xdby,xhe,xhes_cup,'deep', &
!          itf,jtf,ktf, &
!          its,ite, jts,jte, kts,kte)
!
!c--- normalized mass flux profile
!
!     call cup_up_nms(0,xzu,xz_cup,mentr_rate,cd,kbcon,ktop,ierr,k22, &
!          itf,jtf,ktf, &
!          its,ite, jts,jte, kts,kte)
!
!--- moisture downdraft
!
!     call cup_dd_nms(xzd,xz_cup,cdd,mentrd_rate,jmin,ierr, &
!          1,kdet,z1,                 &
!             kbcon,beta,xlamdd,                      &
!          itf,jtf,ktf, &
!          its,ite, jts,jte, kts,kte)
!     call cup_dd_he(xhes_cup,xzd,xhcd,xz_cup,cdd,mentrd_rate, &
!          jmin,ierr,xhe,dbyd,xhe_cup,&
!          itf,jtf,ktf, &
!          its,ite, jts,jte, kts,kte)
!     call cup_dd_moisture_3d(xzd,xhcd,xhes_cup,xqcd,xqes_cup, &
!          xpwd,xq_cup,xz_cup,cdd,mentrd_rate,jmin,ierr,gamma_cup, &
!          xpwev,bu,xqrcd,xq,xhe,xt_cup,3,xl,high_resolution, &
!          itf,jtf,ktf, &
!          its,ite, jts,jte, kts,kte)

!
!------- MOISTURE updraft
!
!     call cup_up_moisture('deep',ierr,xz_cup,xqc,xqrc,xpw,xpwav, &
!          ccnclean,p_cup,kbcon,ktop,cd,xdby,mentr_rate,clw_all, &
!          t_cup,xq,GAMMA_cup,xzu,xqes_cup,k22,xq_cup,xl,ccn,rho,psum,psumh, &
!          autoconv,aeroevap,0,itf,jtf,ktf, &
!          its,ite, jts,jte, kts,kte)
!
!--- workfunctions for updraft
!
      call cup_up_aa0(xaa0,xz,xzu,xdby,GAMMA_CUP,xt_cup, &
           kbcon,ktop,ierr,           &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      do 200 nens=1,maxens
      do i=its,itf
         if(ierr(i).eq.0)then
           xaa0_ens(i,nens)=xaa0(i)
           nall=(iens-1)*maxens3*maxens*maxens2 &
                +(iedt-1)*maxens*maxens3 &
                +(nens-1)*maxens3
           do k=kts,ktf
              if(k.le.ktop(i))then
                 do nens3=1,maxens3
                 if(nens3.eq.7)then
!--- b=0
                 pr_ens(i,j,nall+nens3)=pr_ens(i,j,nall+nens3)  &
!                                +edto(i)*pwdo(i,k)             &
                                    +pwo(i,k)
!--- b=beta
                 else if(nens3.eq.8)then
                 pr_ens(i,j,nall+nens3)=pr_ens(i,j,nall+nens3)+ &
                                    pwo(i,k)
!--- b=beta/2
                 else if(nens3.eq.9)then
                 pr_ens(i,j,nall+nens3)=pr_ens(i,j,nall+nens3)  &
!                                +.5*edto(i)*pwdo(i,k)          &
                                 +  pwo(i,k)
                 else
                 pr_ens(i,j,nall+nens3)=pr_ens(i,j,nall+nens3)+ &
                                    pwo(i,k) ! +edto(i)*pwdo(i,k)
                 endif
                 enddo
              endif
           enddo
         if(pr_ens(i,j,nall+7).lt.1.e-6)then
            ierr(i)=18
            ierrc(i)="total normalized condensate too small"
            do nens3=1,maxens3
               pr_ens(i,j,nall+nens3)=0.
            enddo
         endif
         do nens3=1,maxens3
           if(pr_ens(i,j,nall+nens3).lt.1.e-4)then
            pr_ens(i,j,nall+nens3)=0.
           endif
         enddo
         endif
      enddo
 200  continue
!
!--- LARGE SCALE FORCING
!
!
!------- CHECK wether aa0 should have been zero
!
!
      CALL cup_MAXIMI(HEO_CUP,3,KBMAX,K22x,ierr, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      do i=its,itf
         ierr2(i)=ierr(i)
         ierr3(i)=ierr(i)
      enddo
      call cup_kbcon(ierrc,cap_max_increment,2,k22x,kbconx,heo_cup, &
           heso_cup,ierr2,kbmax,po_cup,cap_max, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      call cup_kbcon(ierrc,cap_max_increment,3,k22x,kbconx,heo_cup, &
           heso_cup,ierr3,kbmax,po_cup,cap_max, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!

      call cup_forcing_ens_3d(closure_n,xland1,aa0,aa1,xaa0_ens,mbdt_ens,dtime,   &
           ierr,ierr2,ierr3,xf_ens,j,'deeps',axx,                 &
           maxens,iens,iedt,maxens2,maxens3,mconv,            &
           po_cup,ktop,omeg,zdo,k22,zuo,pr_ens,edto,kbcon,    &
           massflx,iact,direction,ensdim,massfln,ichoice,     &
           high_resolution,itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte,ens4,ktau)
!
      do k=kts,ktf
      do i=its,itf
        if(ierr(i).eq.0)then
           subt_ens(i,k,iedt)=dsubt(i,k)
           subq_ens(i,k,iedt)=dsubq(i,k)
           dellat_ens(i,k,iedt)=dellat(i,k)
           dellaq_ens(i,k,iedt)=dellaq(i,k)
           dellaqc_ens(i,k,iedt)=dellaqc(i,k)
           pwo_ens(i,k,iedt)=pwo(i,k)+edt(i)*pwdo(i,k)
        else
           subt_ens(i,k,iedt)=0.
           subq_ens(i,k,iedt)=0.
           dellat_ens(i,k,iedt)=0.
           dellaq_ens(i,k,iedt)=0.
           dellaqc_ens(i,k,iedt)=0.
           pwo_ens(i,k,iedt)=0.
        endif
!       if(i.eq.ipr.and.j.eq.jpr)then
!         write(0,*)'1',iens,iedt,dellat(i,k),dellat_ens(i,k,iedt), &
!           dellaq(i,k), dellaqc(i,k)
          !write(0,*)'2',k,subt_ens(i,k,iedt),subq_ens(i,k,iedt)
!       endif
      enddo
      enddo
 250  continue
!
!--- FEEDBACK
!

      call cup_output_ens_3d(xf_ens,ierr,dellat_ens,dellaq_ens, &
            dellaqc_ens,subt_ens,subq_ens,subt,subq,outt,     &
            outq,outqc,zuo,sub_mas,pre,pwo_ens,xmb,ktop,      &
            j,'deep',maxens2,maxens,iens,ierr2,ierr3,         &
            pr_ens,maxens3,ensdim,massfln,                    &
            sig,APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                &
            APR_CAPMA,APR_CAPME,APR_CAPMI,closure_n,xland1,   &
            weight_GR,weight_W,weight_MC,weight_ST,weight_AS,training, &
 	   itf,jtf,ktf,                        &
            its,ite, jts,jte, kts,kte  )

      k=1
      do i=its,itf
          if(ierr(i).eq.0.and.ierr5(i).eq.0.and.kbcon(i).lt.ktop3(i)+1)then
!            write(0,*)'both ier and ier5=0 at i,j=',i,j,kbcon(i),ktop3(i)
             if(high_resolution.eq.1)then
                outts(i,kts:kte)=0.
                outqs(i,kts:kte)=0.
             endif
          elseif (ierr5(i).eq.0)then
!            write(0,*)'ier5=0 at i,j=',i,j,k23(i),ktop3(i)
          endif

           PRE(I)=MAX(PRE(I),0.)
!           if(i.eq.ipr.and.j.eq.jpr)then
!             write(0,*)'i,j,pre(i),aa0(i),aa1(i)'
!             write(0,*)i,j,pre(i),aa0(i)
!           endif
      enddo
!
!---------------------------done------------------------------
!
!      do i=its,itf
!        if(ierr(i).eq.0)then
!       if(i.eq.ipr.and.j.eq.jpr)then
!         write(0,*)'on output, pre =',pre(i),its,itf,kts,ktf
!         do k=kts,ktf
!           write(0,*)z(i,k),outt(i,k)*86400.,subt(i,k)*86400.
!         enddo
!         write(0,*)i,j,(axx(i,k),k=1,ens4)
!       endif
!       endif
!      enddo
!     print *,'ierr(i) = ',ierr(i),pre(i)

   END SUBROUTINE CUP_enss_3d


   SUBROUTINE cup_dd_aa0(edt,ierr,aa0,jmin,gamma_cup,t_cup, &
              hcd,hes_cup,z,zd,                             &
              itf,jtf,ktf,                    &
              its,ite, jts,jte, kts,kte                    )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,                                     &
        its,ite, jts,jte, kts,kte
  ! aa0 cloud work function for downdraft
  ! gamma_cup = gamma on model cloud levels
  ! t_cup = temperature (Kelvin) on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! hcd = moist static energy in downdraft
  ! edt = epsilon
  ! zd normalized downdraft mass flux
  ! z = heights of model levels
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        z,zd,gamma_cup,t_cup,hes_cup,hcd
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        edt
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        jmin
!
! input and output
!


     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        aa0
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,kk
     real                                 ::                           &
        dz
!
       do i=its,itf
        aa0(i)=0.
       enddo
!
!??    DO k=kts,kte-1
       DO k=kts,ktf-1
       do i=its,itf
         IF(ierr(I).eq.0.and.k.lt.jmin(i))then
         KK=JMIN(I)-K
!
!--- ORIGINAL
!
         DZ=(Z(I,KK)-Z(I,KK+1))
         AA0(I)=AA0(I)+zd(i,kk)*EDT(I)*DZ*(9.81/(1004.*T_cup(I,KK))) &
            *((hcd(i,kk)-hes_cup(i,kk))/(1.+GAMMA_cup(i,kk)))
         endif
      enddo
      enddo

   END SUBROUTINE CUP_dd_aa0


   SUBROUTINE cup_dd_edt(ierr,us,vs,z,ktop,kbcon,edt,p,pwav, &
              pw,ccn,pwev,edtmax,edtmin,maxens2,edtc,psum2,psumh, &
              ccnclean,rho,aeroevap,itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        aeroevap,itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
     integer, intent (in   )              ::                           &
        maxens2
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        rho,us,vs,z,p,pw
     real,    dimension (its:ite,1:maxens2)                            &
        ,intent (out  )                   ::                           &
        edtc
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        edt
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        pwav,pwev,ccn,psum2,psumh
     real                                                              &
        ,intent (in   )                   ::                           &
        ccnclean,edtmax,edtmin
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        ktop,kbcon
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!

     integer i,k,kk
     real    einc,pef,pefb,prezk,zkbc
     real,    dimension (its:ite)         ::                           &
      vshear,sdp,vws
     real :: prop_c,pefc,aeroadd,alpha3,beta3,rhoc
     prop_c=8. !10.386
     alpha3 = 1.9
     beta3  = -1.13
     pefc=0.

!
!--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
!
! */ calculate an average wind shear over the depth of the cloud
!
       do i=its,itf
        edt(i)=0.
        vws(i)=0.
        sdp(i)=0.
        vshear(i)=0.
       enddo
       do k=1,maxens2
       do i=its,itf
        edtc(i,k)=0.
       enddo
       enddo
       do kk = kts,ktf-1
         do 62 i=its,itf
          IF(ierr(i).ne.0)GO TO 62
          if (kk .le. min0(ktop(i),ktf) .and. kk .ge. kbcon(i)) then
             vws(i) = vws(i)+ &
              (abs((us(i,kk+1)-us(i,kk))/(z(i,kk+1)-z(i,kk))) &
          +   abs((vs(i,kk+1)-vs(i,kk))/(z(i,kk+1)-z(i,kk)))) * &
              (p(i,kk) - p(i,kk+1))
            sdp(i) = sdp(i) + p(i,kk) - p(i,kk+1)
          endif
          if (kk .eq. ktf-1)vshear(i) = 1.e3 * vws(i) / sdp(i)
   62   continue
       end do
      do i=its,itf
         IF(ierr(i).eq.0)then
            pef=(1.591-.639*VSHEAR(I)+.0953*(VSHEAR(I)**2) &
               -.00496*(VSHEAR(I)**3))
            if(pef.gt.0.9)pef=0.9
            if(pef.lt.0.1)pef=0.1
!
!--- cloud base precip efficiency
!
            zkbc=z(i,kbcon(i))*3.281e-3
            prezk=.02
            if(zkbc.gt.3.)then
               prezk=.96729352+zkbc*(-.70034167+zkbc*(.162179896+zkbc &
               *(- 1.2569798E-2+zkbc*(4.2772E-4-zkbc*5.44E-6))))
            endif
            if(zkbc.gt.25)then
               prezk=2.4
            endif
            pefb=1./(1.+prezk)
            !write(11,*)pefb,prezk,zkbc
            if(pefb.gt.0.9)pefb=0.9
            if(pefb.lt.0.1)pefb=0.1
            EDT(I)=1.-.5*(pefb+pef)
            if(aeroevap.gt.1)then
               aeroadd=(ccnclean**beta3)*((psumh(i))**(alpha3-1)) !*1.e6
               prop_c=.9/aeroadd
               !write(11,*) "aeroadd clean condition (100 cm^-3) = ",aeroadd,prop_c
               prop_c=.5*(pefb+pef)/aeroadd
               !write(11,*) "aeroadd clean condition (100 cm^-3) = ",aeroadd,prop_c
               aeroadd=(ccn(i)**beta3)*((psum2(i))**(alpha3-1)) !*1.e6
               aeroadd=prop_c*aeroadd
               !write(11,*) "aeroadd final, prop constant= ",aeroadd,prop_c
               pefc=aeroadd
               if(pefc.gt.0.9)pefc=0.9
               if(pefc.lt.0.1)pefc=0.1
               EDT(I)=1.-pefc
               if(aeroevap.eq.2)EDT(I)=1.-.25*(pefb+pef+2.*pefc)
            endif
            !write(11,*)aeroevap
            !write(11,*)'edt,pef,pefb,pefc = ',edt(i),pef,pefb,pefc


!--- edt here is 1-precipeff!
            einc=.2*edt(i)
            do k=1,maxens2
                edtc(i,k)=edt(i)+float(k-2)*einc
            enddo
         endif
      enddo
      do i=its,itf
         IF(ierr(i).eq.0)then
            do k=1,maxens2
               EDTC(I,K)=-EDTC(I,K)*PWAV(I)/PWEV(I)
               !write(11,*)'k,edtc,pw,pe = ',edtc(i,k),PWAV(I),PWEV(I)
               IF(EDTC(I,K).GT.edtmax)EDTC(I,K)=edtmax
               IF(EDTC(I,K).LT.edtmin)EDTC(I,K)=edtmin
            enddo
         endif
      enddo

   END SUBROUTINE cup_dd_edt


   SUBROUTINE cup_dd_he(hes_cup,zd,hcd,z_cup,cdd,entr,       &
              jmin,ierr,he,dby,he_cup,                       &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  itf,jtf,ktf,           &
                                  its,ite, jts,jte, kts,kte
  ! hcd = downdraft moist static energy
  ! he = moist static energy on model levels
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! dby = buoancy term
  ! cdd= detrainment function
  ! z_cup = heights of model cloud levels
  ! entr = entrainment rate
  ! zd   = downdraft normalized mass flux
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        he,he_cup,hes_cup,z_cup,cdd,zd
  ! entr= entrainment rate
     real                                                              &
        ,intent (in   )                   ::                           &
        entr
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        jmin
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        hcd,dby
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,ki
     real                                 ::                           &
        dz


      do k=kts+1,ktf
      do i=its,itf
      dby(i,k)=0.
      IF(ierr(I).eq.0)then
         hcd(i,k)=hes_cup(i,k)
      endif
      enddo
      enddo
!
      do 100 i=its,itf
      IF(ierr(I).eq.0)then
      k=jmin(i)
      hcd(i,k)=hes_cup(i,k)
      dby(i,k)=hcd(i,jmin(i))-hes_cup(i,k)
!
      do ki=jmin(i)-1,1,-1
         DZ=Z_cup(i,Ki+1)-Z_cup(i,Ki)
         HCD(i,Ki)=(HCD(i,Ki+1)*(1.-.5*CDD(i,Ki+1)*DZ) &
                  +entr*DZ*HE(i,Ki) &
                  )/(1.+entr*DZ-.5*CDD(i,Ki+1)*DZ)
         dby(i,ki)=HCD(i,Ki)-hes_cup(i,ki)
      enddo
!
      endif
!--- end loop over i
100    continue


   END SUBROUTINE cup_dd_he


   SUBROUTINE cup_dd_moisture_new(ierrc,zd,hcd,hes_cup,qcd,qes_cup,    &
              pwd,q_cup,z_cup,dd_massentr,dd_massdetr,jmin,ierr,            &
              gamma_cup,pwev,bu,qrcd,                        &
              q,he,t_cup,iloop,xl,high_resolution,           &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  itf,jtf,ktf,           &
                                  its,ite, jts,jte, kts,kte,high_resolution
  ! cdd= detrainment function
  ! q = environmental q on model levels
  ! q_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! hes_cup = saturation h on model cloud levels
  ! hcd = h in model cloud
  ! bu = buoancy term
  ! zd = normalized downdraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! mentr_rate = entrainment rate
  ! qcd = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! pwd = evaporate at that level
  ! pwev = total normalized integrated evaoprate (I2)
  ! entr= entrainment rate
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        zd,t_cup,hes_cup,hcd,qes_cup,q_cup,z_cup,                      &
        dd_massentr,dd_massdetr,gamma_cup,q,he
     real                                                              &
        ,intent (in   )                   ::                           &
        xl
     integer                                                           &
        ,intent (in   )                   ::                           &
        iloop
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        jmin
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qcd,qrcd,pwd
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pwev,bu
     character*50 :: ierrc(its:ite)
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,ki
     real                                 ::                           &
        dh,dz,dqeva

      do i=its,itf
         bu(i)=0.
         pwev(i)=0.
      enddo
      do k=kts,ktf
      do i=its,itf
         qcd(i,k)=0.
         qrcd(i,k)=0.
         pwd(i,k)=0.
      enddo
      enddo
!
!
!
      do 100 i=its,itf
      IF(ierr(I).eq.0)then
      k=jmin(i)
      DZ=Z_cup(i,K+1)-Z_cup(i,K)
      qcd(i,k)=q_cup(i,k)
      if(high_resolution.eq.1)qcd(i,k)=.5*(qes_cup(i,k)+q_cup(i,k))
      qrcd(i,k)=qes_cup(i,k)
      pwd(i,jmin(i))=min(0.,qcd(i,k)-qrcd(i,k))
      pwev(i)=pwev(i)+pwd(i,jmin(i))
      qcd(i,k)=qes_cup(i,k)
!
      DH=HCD(I,k)-HES_cup(I,K)
      bu(i)=dz*dh
      do ki=jmin(i)-1,1,-1
         DZ=Z_cup(i,Ki+1)-Z_cup(i,Ki)
!        QCD(i,Ki)=(qCD(i,Ki+1)*(1.-.5*CDD(i,Ki+1)*DZ) &
!                 +entr*DZ*q(i,Ki) &
!                )/(1.+entr*DZ-.5*CDD(i,Ki+1)*DZ)
!        dz=qcd(i,ki)
         qcd(i,ki)=(qcd(i,ki+1)*zd(i,ki+1)                          &
                  -.5*dd_massdetr(i,ki+1)*qcd(i,ki+1)+ &
                  dd_massentr(i,ki+1)*q(i,ki+1))   /            &
                  (zd(i,ki+1)-.5*dd_massdetr(i,ki+1)+dd_massentr(i,ki+1))
         !write(0,*)'qcd in dd_moi = ',qcd(i,ki)

!
!--- to be negatively buoyant, hcd should be smaller than hes!
!
         DH=HCD(I,ki)-HES_cup(I,Ki)
         bu(i)=bu(i)+dz*dh
         QRCD(I,Ki)=qes_cup(i,ki)+(1./XL)*(GAMMA_cup(i,ki) &
                  /(1.+GAMMA_cup(i,ki)))*DH
         dqeva=qcd(i,ki)-qrcd(i,ki)
         if(dqeva.gt.0.)dqeva=0.
         pwd(i,ki)=zd(i,ki)*dqeva
         qcd(i,ki)=qrcd(i,ki)
         pwev(i)=pwev(i)+pwd(i,ki)
!        if(iloop.eq.1.and.i.eq.102.and.j.eq.62)then
!         print *,'in cup_dd_moi ', hcd(i,ki),HES_cup(I,Ki),dh,dqeva
!        endif
      enddo
!
!--- end loop over i
       if(pwev(I).eq.0.and.iloop.eq.1)then
!        print *,'problem with buoy in cup_dd_moisture',i
         ierr(i)=7
         ierrc(i)="problem with buoy in cup_dd_moisture"
       endif
       if(BU(I).GE.0.and.iloop.eq.1)then
!        print *,'problem with buoy in cup_dd_moisture',i
         ierr(i)=7
         ierrc(i)="problem2 with buoy in cup_dd_moisture"
       endif
      endif
100    continue

   END SUBROUTINE cup_dd_moisture_new

   SUBROUTINE cup_dd_moisture_3d(ierrc,zd,hcd,hes_cup,qcd,qes_cup,    &
              pwd,q_cup,z_cup,cdd,entr,jmin,ierr,            &
              gamma_cup,pwev,bu,qrcd,                        &
              q,he,t_cup,iloop,xl,high_resolution,           &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  itf,jtf,ktf,           &
                                  its,ite, jts,jte, kts,kte,high_resolution
  ! cdd= detrainment function
  ! q = environmental q on model levels
  ! q_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! hes_cup = saturation h on model cloud levels
  ! hcd = h in model cloud
  ! bu = buoancy term
  ! zd = normalized downdraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! mentr_rate = entrainment rate
  ! qcd = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! pwd = evaporate at that level
  ! pwev = total normalized integrated evaoprate (I2)
  ! entr= entrainment rate
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        zd,t_cup,hes_cup,hcd,qes_cup,q_cup,z_cup,cdd,gamma_cup,q,he
     real                                                              &
        ,intent (in   )                   ::                           &
        entr,xl
     integer                                                           &
        ,intent (in   )                   ::                           &
        iloop
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        jmin
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qcd,qrcd,pwd
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pwev,bu
     character*50 :: ierrc(its:ite)
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,ki
     real                                 ::                           &
        dh,dz,dqeva

      do i=its,itf
         bu(i)=0.
         pwev(i)=0.
      enddo
      do k=kts,ktf
      do i=its,itf
         qcd(i,k)=0.
         qrcd(i,k)=0.
         pwd(i,k)=0.
      enddo
      enddo
!
!
!
      do 100 i=its,itf
      IF(ierr(I).eq.0)then
      k=jmin(i)
      DZ=Z_cup(i,K+1)-Z_cup(i,K)
      qcd(i,k)=q_cup(i,k)
      if(high_resolution.eq.1)qcd(i,k)=.5*(qes_cup(i,k)+q_cup(i,k))
      qrcd(i,k)=qes_cup(i,k)
      pwd(i,jmin(i))=min(0.,qcd(i,k)-qrcd(i,k))
      pwev(i)=pwev(i)+pwd(i,jmin(i))
      qcd(i,k)=qes_cup(i,k)
!
      DH=HCD(I,k)-HES_cup(I,K)
      bu(i)=dz*dh
      do ki=jmin(i)-1,1,-1
         DZ=Z_cup(i,Ki+1)-Z_cup(i,Ki)
         QCD(i,Ki)=(qCD(i,Ki+1)*(1.-.5*CDD(i,Ki+1)*DZ) &
                  +entr*DZ*q(i,Ki) &
                  )/(1.+entr*DZ-.5*CDD(i,Ki+1)*DZ)
         !write(0,*)'qcd in dd_moiold = ',qcd(i,ki)
!
!--- to be negatively buoyant, hcd should be smaller than hes!
!
         DH=HCD(I,ki)-HES_cup(I,Ki)
         bu(i)=bu(i)+dz*dh
         QRCD(I,Ki)=qes_cup(i,ki)+(1./XL)*(GAMMA_cup(i,ki) &
                  /(1.+GAMMA_cup(i,ki)))*DH
         dqeva=qcd(i,ki)-qrcd(i,ki)
         if(dqeva.gt.0.)dqeva=0.
         pwd(i,ki)=zd(i,ki)*dqeva
         qcd(i,ki)=qrcd(i,ki)
         pwev(i)=pwev(i)+pwd(i,ki)
!        if(iloop.eq.1.and.i.eq.102.and.j.eq.62)then
!         print *,'in cup_dd_moi ', hcd(i,ki),HES_cup(I,Ki),dh,dqeva
!        endif
      enddo
!
!--- end loop over i
       if(pwev(I).eq.0.and.iloop.eq.1)then
!        print *,'problem with buoy in cup_dd_moisture',i
         ierr(i)=7
         ierrc(i)="problem with buoy in cup_dd_moisture"
       endif
       if(BU(I).GE.0.and.iloop.eq.1)then
!        print *,'problem with buoy in cup_dd_moisture',i
         ierr(i)=7
         ierrc(i)="problem2 with buoy in cup_dd_moisture"
       endif
      endif
100    continue

   END SUBROUTINE cup_dd_moisture_3d


   SUBROUTINE cup_dd_nms(zd,z_cup,cdd,entr,jmin,ierr,        &
              itest,kdet,z1,                                 &
              kbcon,beta,xlamdd,                      &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  itf,jtf,ktf,           &
                                  its,ite, jts,jte, kts,kte
  ! z_cup = height of cloud model level
  ! z1 = terrain elevation
  ! entr = downdraft entrainment rate
  ! jmin = downdraft originating level
  ! kdet = level above ground where downdraft start detraining
  ! itest = flag to whether to calculate cdd

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        z_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        z1

     real                                                              &
        ,intent (in   )                   ::                           &
        entr,xlamdd,beta
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        jmin,kdet,kbcon
     integer                                                           &
        ,intent (in   )                   ::                           &
        itest
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
                                                                 ierr
   ! zd is the normalized downdraft mass flux
   ! cdd is the downdraft detrainmen function

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
                                                             zd
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
                                                             cdd
!
!  local variables in this routine
!

     integer                              ::                           &
                                                  i,k,ki
     real                                 ::                           &
                                            a,perc,dz,xlamd,tem

!
!--- perc is the percentage of mass left when hitting the ground
!
      perc=.03

      do k=kts,ktf
      do i=its,itf
         zd(i,k)=0.
!        if(itest.eq.0)cdd(i,k)=0.
      enddo
      enddo
      a=1.-perc
!
!
!
      if (itest.eq.0) then

        !write(0,*)'itest = 0 here '
        do k=kts,ktf
        do i=its,itf
          cdd(i,k)=0.
        enddo
        enddo

        do i=its,itf
          IF(ierr(I).eq.0)then
            dz = z_cup(i,kbcon(i))/float(kbcon(i))
            tem = 1./float(kbcon(i))
            xlamd = (1.-beta**tem)/dz
            !write(0,*)'b',beta**tem,tem,beta,xlamd,xlamdd

            do ki=jmin(i)-0,1,-1
              if (ki.ge.kbcon(i)) then
                cdd(i,ki) = xlamdd
              else
                cdd(i,ki) = xlamdd + xlamd
              endif
!           if(jpr.eq.1)!write(6,*)ki,cdd(i,ki),xlamdd,xlamd
            end do

          ENDIF
        enddo
      endif ! itest=0

      do 100 i=its,itf
      IF(ierr(I).eq.0)then
      zd(i,jmin(i))=1.
!
!--- integrate downward, specify detrainment(cdd)!
!
      do ki=jmin(i)-1,1,-1
         DZ=Z_cup(i,Ki+1)-Z_cup(i,Ki)
         zd(i,ki)=zd(i,ki+1)*(1.+(entr-cdd(i,ki+1))*dz)
      enddo
!
      endif
!--- end loop over i
100    continue

   END SUBROUTINE cup_dd_nms


   SUBROUTINE cup_dellabot(name,ipr,jpr,he_cup,ierr,z_cup,p_cup,  &
              hcd,edt,zd,cdd,he,della,subs,j,mentrd_rate,z,g,     &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
     integer, intent (in   )              ::                           &
        j,ipr,jpr
      character *(*), intent (in)        ::                           &
       name
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        della,subs
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in  )                   ::                           &
        z_cup,p_cup,hcd,zd,cdd,he,z,he_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        edt
     real                                                              &
        ,intent (in   )                   ::                           &
        g,mentrd_rate
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!

      integer i
      real detdo,detdo1,detdo2,entdo,dp,dz,subin,                      &
      totmas
!
!
!     if(name.eq.'shallow')then
!        edt(:)=0.
!        cdd(:,:)=0.
!     endif
      do 100 i=its,itf
      della(i,1)=0.
      subs(i,1)=0.
      if(ierr(i).ne.0)go to 100
      dz=z_cup(i,2)-z_cup(i,1)
      DP=100.*(p_cup(i,1)-P_cup(i,2))
      detdo1=edt(i)*zd(i,2)*CDD(i,1)*DZ
      detdo2=edt(i)*zd(i,1)
      entdo=edt(i)*zd(i,2)*mentrd_rate*dz
      subin=-EDT(I)*zd(i,2)
      detdo=detdo1+detdo2-entdo+subin
      DELLA(I,1)=(detdo1*.5*(HCD(i,1)+HCD(i,2)) &
                 +detdo2*hcd(i,1) &
                 +subin*he_cup(i,2) &
                 -entdo*he(i,1))*g/dp
       della(i,1)=(edt(i)*zd(i,2)*hcd(i,2) &
                  - edt(i)*zd(i,2)*he_cup(i,2))*g/dp
      SUBS(I,1)=0.
!      if(i.eq.ipr.and.j.eq.jpr)then
!       write(0,*)'db1',della(i,1),subs(i,1),subin,entdo
!       write(0,*)'db2',detdo1,detdo2,detdo1+detdo2-entdo+subin
!      endif
 100  CONTINUE

   END SUBROUTINE cup_dellabot


   SUBROUTINE cup_dellas_3d(ierr,z_cup,p_cup,hcd,edt,zd,cdd,              &
              he,della,subs,j,mentrd_rate,zu,g,                             &
              cd,hc,ktop,k22,kbcon,mentr_rate,jmin,he_cup,kdet,kpbl,   &
              ipr,jpr,name,high_res,                                            &
              itf,jtf,ktf,                               &
              its,ite, jts,jte, kts,kte                               )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
     integer, intent (in   )              ::                           &
        j,ipr,jpr,high_res
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        della,subs
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in  )                   ::                           &
        z_cup,p_cup,hcd,zd,cdd,he,hc,cd,zu,he_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        edt
     real                                                              &
        ,intent (in   )                   ::                           &
        g,mentrd_rate,mentr_rate
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop,k22,jmin,kdet,kpbl
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
      character *(*), intent (in)        ::                           &
       name
!
!  local variables in this routine
!

      integer i,k,kstart
      real detdo1,detdo2,entdo,dp,dz,subin,detdo,entup,                &
      detup,subdown,entdoj,entupk,detupk,totmas
!
      i=ipr
      kstart=kts+1
      if(name.eq.'shallow')kstart=kts
       DO K=kstart,ktf
       do i=its,itf
          della(i,k)=0.
          subs(i,k)=0.
       enddo
       enddo
!
! no downdrafts for shallow convection
!
       DO 100 k=kts+1,ktf-1
       DO 100 i=its,ite
         IF(ierr(i).ne.0)GO TO 100
         IF(K.Gt.KTOP(I))GO TO 100
         if(k.lt.k22(i)-1.and.name.eq.'shallow')GO TO 100
         DZ=Z_cup(I,K+1)-Z_cup(I,K)
! these three are only used at or near mass detrainment and/or entrainment levels
         entupk=0.
         detupk=0.
         entdoj=0.
! detrainment and entrainment for fowndrafts
         detdo=edt(i)*CDD(i,K+1)*DZ*ZD(i,k+1)
         entdo=edt(i)*mentrd_rate*dz*zd(i,k+1)
! entrainment/detrainment for updraft
         entup=mentr_rate*dz*zu(i,k)
         detup=CD(i,K)*DZ*ZU(i,k)
! subsidence by downdrafts only
         subin=-zd(i,k+1)*edt(i)
         subdown=-zd(i,k)*edt(i)
!
!         SPECIAL LEVELS
!
! updraft originates at k22, only updraft term at k22-1 is (zu-entupk)*he
         if(k.eq.k22(i)-1)then
            entupk=zu(i,k+1)
            entup=0.
            detup=0.
         endif
! downdraft originating level, similiar to k22-1 for updraft
         if(k.eq.jmin(i))then
            entdoj=edt(i)*zd(i,k)
            subin=0.
            detdo=0.
            entdo=0.
         endif
         if(k.eq.ktop(i))then
            detupk=zu(i,ktop(i))
            subin=0.
            subdown=0.
            detdo=0.
            entdo=0.
            entup=0.
            detup=0.
         endif
         if(k.lt.kbcon(i))then
            detup=0.
            entup=0.
         endif
!
!--- SPECIFY DETRAINMENT OF DOWNDRAFT, HAS TO BE CONSISTENT
!--- WITH ZD CALCULATIONS IN SOUNDD.
!
!              write(9,111),k,mentr_rate*dz*zu(i,k),cd(i,k)*dz*zu(i,k), &
!                mentr_rate,cd(i,k),dz,zu(i,k)
111 format(1x,i3,4e13.6,f8.1,f7.4)

!C
!C--- CHANGED DUE TO SUBSIDENCE AND ENTRAINMENT
!C
         totmas=subin-subdown+detup-entup-entdo+ &
             detdo-entupk-entdoj+detupk+zu(i,k+1)-zu(i,k)
!         if(j.eq.jpr.and.i.eq.ipr)print *,'k,totmas,sui,sud = ',k,
!     1   totmas,subin,subdown
!         if(j.eq.jpr.and.i.eq.ipr)print *,'updr stuff = ',detup,
!     1      entup,entupk,detupk
!         if(j.eq.jpr.and.i.eq.ipr)print *,'dddr stuff = ',entdo,
!     1      detdo,entdoj
         if(abs(totmas).gt.1.e-6)then
           print *,'*********************',i,j,k,totmas,name
           print *,jmin(i),k22(i),kbcon(i),ktop(i)
          print *,'updr stuff = ',subin, &
           subdown,detup,entup,entupk,detupk
!        call wrf_error_fatal ( 'totmas .gt.1.e-6' )
         endif
         dp=100.*(p_cup(i,k-1)-p_cup(i,k))
         della(i,k)=(detup*.5*(HC(i,K+1)+HC(i,K)) &
                    +detdo*.5*(HCD(i,K+1)+HCD(i,K)) &
                    -entup*he(i,k) &
                    -entdo*he(i,k) &
                    +subin*he_cup(i,k+1) &
                    -subdown*he_cup(i,k) &
                    +detupk*(hc(i,ktop(i))-he_cup(i,ktop(i)))    &
!                   +detupk*hc(i,ktop(i))    &
                    -entupk*he_cup(i,k22(i)) &
                    -entdoj*he_cup(i,jmin(i)) &
                     )*g/dp
           if(high_res.eq.1)then
! the first term includes entr and detr into/from updraft as well as (entup-detup)*he(i,k) from
!  neighbouring point, to make things mass consistent....
!            if(k.ge.k22(i))then
                della(i,k)=(                          &
                    detup*.5*(HC(i,K+1)+HC(i,K))-entup*he(i,k)+(entup-detup)*he(i,k) &
                    +detdo*.5*(HCD(i,K+1)+HCD(i,K)) &
                    -entdo*he(i,k) &
                    +subin*he_cup(i,k+1) &
                    -subdown*he_cup(i,k) &
                    +detupk*(hc(i,ktop(i))-he(i,ktop(i)))    &
                    -entdoj*he_cup(i,jmin(i)) &
                    -entupk*he_cup(i,k22(i))+entupk*he(i,k) &
                     )*g/dp
!             else if(k.eq.k22(i)-1)then
!                  della(i,k)=(-entupk*he_cup(i,k22(i))+entupk*he(i,k))*g/dp
           endif
!3d        subin=zu(i,k+1)-zd(i,k+1)*edt(i)
!
! updraft subsidence only
!
         if(k.ge.k22(i)-1.and.k.lt.ktop(i))then
         subs(i,k)=(zu(i,k+1)*he_cup(i,k+1) &
                    -zu(i,k)*he_cup(i,k))*g/dp
         !write(0,*)'subs = ',zu(i,k+1)*he_cup(i,k+1),zu(i,k)*he_cup(i,k),subs(i,k)
!        else if(k.eq.ktop(i))then
!        subs(i,k)=-detupk*he_cup(i,k)*g/dp
         endif
           if(k.le.kbcon(i))then
             !write(8,*)k,detup*.5*(HC(i,K+1)+HC(i,K)), &
             !            detdo*.5*(HCD(i,K+1)+HCD(i,K)), &
             !            entup*he(i,k),                            &
             !            detupk*(hc(i,ktop(i))-he_cup(i,ktop(i))), &
             !            entupk*he_cup(i,k22(i)),                  &
             !            entdo*he(i,k),                            &
             !            entdoj*he_cup(i,jmin(i)), &
             !            zu(i,k+1)*he_cup(i,k+1),  &
             !            zu(i,k)*he_cup(i,k),      &
             !            subin*he_cup(i,k+1),      &
             !            subdown*he_cup(i,k)
           endif
!
! in igh res case, subsidence terms are for meighbouring points only. This has to be
! done mass consistent with the della term
         if(high_res.eq.1)then
            if(k.ge.k22(i).and.k.lt.ktop(i))then
               subs(i,k)=(zu(i,k+1)*he_cup(i,k+1)-zu(i,k)*he_cup(i,k)-(entup-detup)*he(i,k))*g/dp
            else if(k.eq.ktop(i))then
               subs(i,k)=detupk*(he(i,ktop(i))-he_cup(i,ktop(i)))*g/dp
            else if(k.eq.k22(i)-1)then
               subs(i,k)=(entupk*he(i,k)-entupk*he_cup(i,k))*g/dp
         endif
         endif
!       if(i.eq.ipr.and.j.eq.jpr)then
!         write(0,*)'d',k,della(i,k),subs(i,k),subin,subdown
!!        write(0,*)'d',detup,entup,entdo,entupk,entdoj
!!        print *,k,della(i,k),subin*he_cup(i,k+1),subdown*he_cup(i,k),
!!     1            detdo*.5*(HCD(i,K+1)+HCD(i,K))
!!        print *,k,detup*.5*(HC(i,K+1)+HC(i,K)),detupk*hc(i,ktop(i)),
!!     1         entup*he(i,k),entdo*he(i,k)
!!        print *,k,he_cup(i,k+1),he_cup(i,k),entupk*he_cup(i,k)
!       endif

 100  CONTINUE

   END SUBROUTINE cup_dellas_3d


   SUBROUTINE cup_direction2(i,j,dir,id,massflx,             &
              iresult,imass,massfld,                         &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
     integer, intent (in   )              ::                           &
        i,j,imass
     integer, intent (out  )              ::                           &
        iresult
  !
  ! ierr error value, maybe modified in this routine
  !
     integer,    dimension (its:ite,jts:jte)                           &
        ,intent (in   )                   ::                           &
        id
     real,    dimension (its:ite,jts:jte)                              &
        ,intent (in   )                   ::                           &
        massflx
     real,    dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        dir
     real                                                              &
        ,intent (out  )                   ::                           &
        massfld
!
!  local variables in this routine
!

       integer k,ia,ja,ib,jb
       real diff
!
!
!
       if(imass.eq.1)then
           massfld=massflx(i,j)
       endif
       iresult=0
!      return
       diff=22.5
       if(dir(i).lt.22.5)dir(i)=360.+dir(i)
       if(id(i,j).eq.1)iresult=1
!      ja=max(2,j-1)
!      ia=max(2,i-1)
!      jb=min(mjx-1,j+1)
!      ib=min(mix-1,i+1)
       ja=j-1
       ia=i-1
       jb=j+1
       ib=i+1
        if(dir(i).gt.90.-diff.and.dir(i).le.90.+diff)then
!--- steering flow from the east
          if(id(ib,j).eq.1)then
            iresult=1
            if(imass.eq.1)then
               massfld=max(massflx(ib,j),massflx(i,j))
            endif
            return
          endif
        else if(dir(i).gt.135.-diff.and.dir(i).le.135.+diff)then
!--- steering flow from the south-east
          if(id(ib,ja).eq.1)then
            iresult=1
            if(imass.eq.1)then
               massfld=max(massflx(ib,ja),massflx(i,j))
            endif
            return
          endif
!--- steering flow from the south
        else if(dir(i).gt.180.-diff.and.dir(i).le.180.+diff)then
          if(id(i,ja).eq.1)then
            iresult=1
            if(imass.eq.1)then
               massfld=max(massflx(i,ja),massflx(i,j))
            endif
            return
          endif
!--- steering flow from the south west
        else if(dir(i).gt.225.-diff.and.dir(i).le.225.+diff)then
          if(id(ia,ja).eq.1)then
            iresult=1
            if(imass.eq.1)then
               massfld=max(massflx(ia,ja),massflx(i,j))
            endif
            return
          endif
!--- steering flow from the west
        else if(dir(i).gt.270.-diff.and.dir(i).le.270.+diff)then
          if(id(ia,j).eq.1)then
            iresult=1
            if(imass.eq.1)then
               massfld=max(massflx(ia,j),massflx(i,j))
            endif
            return
          endif
!--- steering flow from the north-west
        else if(dir(i).gt.305.-diff.and.dir(i).le.305.+diff)then
          if(id(ia,jb).eq.1)then
            iresult=1
            if(imass.eq.1)then
               massfld=max(massflx(ia,jb),massflx(i,j))
            endif
            return
          endif
!--- steering flow from the north
        else if(dir(i).gt.360.-diff.and.dir(i).le.360.+diff)then
          if(id(i,jb).eq.1)then
            iresult=1
            if(imass.eq.1)then
               massfld=max(massflx(i,jb),massflx(i,j))
            endif
            return
          endif
!--- steering flow from the north-east
        else if(dir(i).gt.45.-diff.and.dir(i).le.45.+diff)then
          if(id(ib,jb).eq.1)then
            iresult=1
            if(imass.eq.1)then
               massfld=max(massflx(ib,jb),massflx(i,j))
            endif
            return
          endif
        endif

   END SUBROUTINE cup_direction2


   SUBROUTINE cup_env(z,qes,he,hes,t,q,p,z1,                 &
              psur,ierr,tcrit,itest,xl,cp,                   &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
  !
  ! ierr error value, maybe modified in this routine
  ! q           = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! tv          = environmental virtual temp
  ! p           = environmental pressure
  ! z           = environmental heights
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  ! psur        = surface pressure
  ! z1          = terrain elevation
  !
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        p,t,q
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        he,hes,qes
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
        z
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        psur,z1
     real                                                              &
        ,intent (in   )                   ::                           &
        xl,cp
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     integer                                                           &
        ,intent (in   )                   ::                           &
        itest
!
!  local variables in this routine
!

     integer                              ::                           &
       i,k,iph
      real, dimension (1:2) :: AE,BE,HT
      real, dimension (its:ite,kts:kte) :: tv
      real :: tcrit,e,tvbar


      HT(1)=XL/CP
      HT(2)=2.834E6/CP
      BE(1)=.622*HT(1)/.286
      AE(1)=BE(1)/273.+ALOG(610.71)
      BE(2)=.622*HT(2)/.286
      AE(2)=BE(2)/273.+ALOG(610.71)
!      print *, 'TCRIT = ', tcrit,its,ite
      DO k=kts,ktf
      do i=its,itf
        if(ierr(i).eq.0)then
!Csgb - IPH is for phase, dependent on TCRIT (water or ice)
        IPH=1
        IF(T(I,K).LE.TCRIT)IPH=2
!       print *, 'AE(IPH),BE(IPH) = ',AE(IPH),BE(IPH),AE(IPH)-BE(IPH),T(i,k),i,k
        E=EXP(AE(IPH)-BE(IPH)/T(I,K))
!       print *, 'P, E = ', P(I,K), E
        QES(I,K)=.622*E/(100.*P(I,K)-E)
        IF(QES(I,K).LE.1.E-08)QES(I,K)=1.E-08
        IF(QES(I,K).LT.Q(I,K))QES(I,K)=Q(I,K)
!       IF(Q(I,K).GT.QES(I,K))Q(I,K)=QES(I,K)
        TV(I,K)=T(I,K)+.608*Q(I,K)*T(I,K)
        endif
      enddo
      enddo
!
!--- z's are calculated with changed h's and q's and t's
!--- if itest=2
!
      if(itest.ne.2)then
         do i=its,itf
           if(ierr(i).eq.0)then
             Z(I,1)=max(0.,Z1(I))-(ALOG(P(I,1))- &
                 ALOG(PSUR(I)))*287.*TV(I,1)/9.81
           endif
         enddo

! --- calculate heights
         DO K=kts+1,ktf
         do i=its,itf
           if(ierr(i).eq.0)then
              TVBAR=.5*TV(I,K)+.5*TV(I,K-1)
              Z(I,K)=Z(I,K-1)-(ALOG(P(I,K))- &
               ALOG(P(I,K-1)))*287.*TVBAR/9.81
           endif
         enddo
         enddo
      else
         do k=kts,ktf
         do i=its,itf
           if(ierr(i).eq.0)then
             z(i,k)=(he(i,k)-1004.*t(i,k)-2.5e6*q(i,k))/9.81
             z(i,k)=max(1.e-3,z(i,k))
           endif
         enddo
         enddo
      endif
!
!--- calculate moist static energy - HE
!    saturated moist static energy - HES
!
       DO k=kts,ktf
       do i=its,itf
         if(ierr(i).eq.0)then
         if(itest.eq.0)HE(I,K)=9.81*Z(I,K)+1004.*T(I,K)+2.5E06*Q(I,K)
         HES(I,K)=9.81*Z(I,K)+1004.*T(I,K)+2.5E06*QES(I,K)
         IF(HE(I,K).GE.HES(I,K))HE(I,K)=HES(I,K)
         endif
      enddo
      enddo

   END SUBROUTINE cup_env


   SUBROUTINE cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,   &
              he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
              ierr,z1,xl,rv,cp,                                &
              itf,jtf,ktf,                       &
              its,ite, jts,jte, kts,kte                       )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
  !
  ! ierr error value, maybe modified in this routine
  ! q           = environmental mixing ratio
  ! q_cup       = environmental mixing ratio on cloud levels
  ! qes         = environmental saturation mixing ratio
  ! qes_cup     = environmental saturation mixing ratio on cloud levels
  ! t           = environmental temp
  ! t_cup       = environmental temp on cloud levels
  ! p           = environmental pressure
  ! p_cup       = environmental pressure on cloud levels
  ! z           = environmental heights
  ! z_cup       = environmental heights on cloud levels
  ! he          = environmental moist static energy
  ! he_cup      = environmental moist static energy on cloud levels
  ! hes         = environmental saturation moist static energy
  ! hes_cup     = environmental saturation moist static energy on cloud levels
  ! gamma_cup   = gamma on cloud levels
  ! psur        = surface pressure
  ! z1          = terrain elevation
  !
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        qes,q,he,hes,z,p,t
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        psur,z1
     real                                                              &
        ,intent (in   )                   ::                           &
        xl,rv,cp
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!

     integer                              ::                           &
       i,k


      do k=kts+1,ktf
      do i=its,itf
        if(ierr(i).eq.0)then
        qes_cup(i,k)=.5*(qes(i,k-1)+qes(i,k))
        q_cup(i,k)=.5*(q(i,k-1)+q(i,k))
        hes_cup(i,k)=.5*(hes(i,k-1)+hes(i,k))
        he_cup(i,k)=.5*(he(i,k-1)+he(i,k))
        if(he_cup(i,k).gt.hes_cup(i,k))he_cup(i,k)=hes_cup(i,k)
        z_cup(i,k)=.5*(z(i,k-1)+z(i,k))
        p_cup(i,k)=.5*(p(i,k-1)+p(i,k))
        t_cup(i,k)=.5*(t(i,k-1)+t(i,k))
        gamma_cup(i,k)=(xl/cp)*(xl/(rv*t_cup(i,k) &
                       *t_cup(i,k)))*qes_cup(i,k)
        endif
      enddo
      enddo
      do i=its,itf
        if(ierr(i).eq.0)then
        qes_cup(i,1)=qes(i,1)
        q_cup(i,1)=q(i,1)
        hes_cup(i,1)=hes(i,1)
        he_cup(i,1)=he(i,1)
        z_cup(i,1)=.5*(z(i,1)+z1(i))
        p_cup(i,1)=.5*(p(i,1)+psur(i))
        t_cup(i,1)=t(i,1)
        gamma_cup(i,1)=xl/cp*(xl/(rv*t_cup(i,1) &
                       *t_cup(i,1)))*qes_cup(i,1)
        endif
      enddo

   END SUBROUTINE cup_env_clev


   SUBROUTINE cup_forcing_ens_3d(closure_n,xland,aa0,aa1,xaa0,mbdt,dtime,ierr,ierr2,ierr3,&
              xf_ens,j,name,axx,maxens,iens,iedt,maxens2,maxens3,mconv,    &
              p_cup,ktop,omeg,zd,k22,zu,pr_ens,edt,kbcon,massflx,      &
              iact_old_gr,dir,ensdim,massfln,icoic,            &
              high_resolution,itf,jtf,ktf,               &
              its,ite, jts,jte, kts,kte,ens4,ktau                )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte,ens4,high_resolution,ktau
     integer, intent (in   )              ::                           &
        j,ensdim,maxens,iens,iedt,maxens2,maxens3
  !
  ! ierr error value, maybe modified in this routine
  ! pr_ens = precipitation ensemble
  ! xf_ens = mass flux ensembles
  ! massfln = downdraft mass flux ensembles used in next timestep
  ! omeg = omega from large scale model
  ! mconv = moisture convergence from large scale model
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects (ensemble dependent)
  ! edt     = epsilon
  ! dir     = "storm motion"
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! iact_gr_old = flag to tell where convection was active
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! name        = deep or shallow convection flag
  !
     real,    dimension (its:ite,jts:jte,1:ensdim)                     &
        ,intent (inout)                   ::                           &
        pr_ens
     real,    dimension (its:ite,jts:jte,1:ensdim)                     &
        ,intent (out  )                   ::                           &
        xf_ens,massfln
     real,    dimension (its:ite,jts:jte)                              &
        ,intent (in   )                   ::                           &
        massflx
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        zd,zu,p_cup
     real,    dimension (its:ite,kts:kte,1:ens4)                              &
        ,intent (in   )                   ::                           &
        omeg
     real,    dimension (its:ite,1:maxens)                             &
        ,intent (in   )                   ::                           &
        xaa0
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        aa1,edt,dir,xland
     real,    dimension (its:ite,1:ens4)                                      &
        ,intent (in   )                   ::                           &
        mconv,axx
     real,    dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        aa0,closure_n
     real,    dimension (1:maxens)                                     &
        ,intent (in   )                   ::                           &
        mbdt
     real                                                              &
        ,intent (in   )                   ::                           &
        dtime
     integer, dimension (its:ite,jts:jte)                              &
        ,intent (in   )                   ::                           &
        iact_old_gr
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        k22,kbcon,ktop
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr,ierr2,ierr3
     integer                                                           &
        ,intent (in   )                   ::                           &
        icoic
      character *(*), intent (in)         ::                           &
       name
!
!  local variables in this routine
!

     real,    dimension (1:maxens3)       ::                           &
       xff_ens3
     real,    dimension (1:maxens)        ::                           &
       xk
     integer                              ::                           &
       i,k,nall,n,ne,nens,nens3,iresult,iresultd,iresulte,mkxcrt,kclim
     parameter (mkxcrt=15)
     real                                 ::                           &
       fens4,a1,massfld,a_ave,xff0,xff00,xxx,xomg,aclim1,aclim2,aclim3,aclim4
     real,    dimension(1:mkxcrt)         ::                           &
       pcrit,acrit,acritt

     integer :: nall2,ixxx,irandom
     !integer,  dimension (12) :: seed !LFR
     integer,  allocatable, dimension (:) :: seed !lFR

      DATA PCRIT/850.,800.,750.,700.,650.,600.,550.,500.,450.,400.,    &
                 350.,300.,250.,200.,150./
      DATA ACRIT/.0633,.0445,.0553,.0664,.075,.1082,.1521,.2216,       &
                 .3151,.3677,.41,.5255,.7663,1.1686,1.6851/
!  GDAS DERIVED ACRIT
      DATA ACRITT/.203,.515,.521,.566,.625,.665,.659,.688,             &
                  .743,.813,.886,.947,1.138,1.377,1.896/
!
    !LFR adapting to new version of fortran
    call random_seed(size=i)
    allocate(seed(i))

       seed=0
       seed(2)=j
       seed(3)=ktau
       nens=0
       irandom=1
       if(high_resolution.eq.1)irandom=0
       irandom=0
       fens4=float(ens4)

!--- LARGE SCALE FORCING
!
       DO 100 i=its,itf
          if(name.eq.'deeps'.and.ierr(i).gt.995)then
           aa0(i)=0.
           ierr(i)=0
          endif
          IF(ierr(i).eq.0)then
!
!---
!
             if(name.eq.'deeps')then
!
                a_ave=0.
                do ne=1,ens4
                  a_ave=a_ave+axx(i,ne)
                  !write(0,*)'in forcing, a_ave,axx(i,ne) = ',a_ave,axx(i,ne)
                enddo
                a_ave=max(0.,a_ave/fens4)
                a_ave=min(a_ave,aa1(i))
                a_ave=max(0.,a_ave)
                do ne=1,16
                  xff_ens3(ne)=0.
                enddo
                xff0= (AA1(I)-AA0(I))/DTIME
                if(high_resolution.eq.1)xff0= (a_ave-AA0(I))/DTIME
                xff_ens3(1)=(AA1(I)-AA0(I))/dtime
                xff_ens3(2)=(a_ave-AA0(I))/dtime
                if(irandom.eq.1)then
                   seed(1)=i
                   call random_seed (PUT=seed)
                   call random_number (xxx)
                   ixxx=min(ens4,max(1,int(fens4*xxx+1.e-8)))
                   xff_ens3(3)=(axx(i,ixxx)-AA0(I))/dtime
                   call random_number (xxx)
                   ixxx=min(ens4,max(1,int(fens4*xxx+1.e-8)))
                   xff_ens3(13)=(axx(i,ixxx)-AA0(I))/dtime
                else
                   xff_ens3(3)=(AA1(I)-AA0(I))/dtime
                   xff_ens3(13)=(AA1(I)-AA0(I))/dtime
                endif
                if(high_resolution.eq.1)then
                   xff_ens3(1)=(a_ave-AA0(I))/dtime
                   xff_ens3(2)=(a_ave-AA0(I))/dtime
                   xff_ens3(3)=(a_ave-AA0(I))/dtime
                   xff_ens3(13)=(a_ave-AA0(I))/dtime
                endif
!
!--- more original Arakawa-Schubert (climatologic value of aa0)
!
!
!--- omeg is in bar/s, mconv done with omeg in Pa/s
!     more like Brown (1979), or Frank-Cohen (199?)
!
                xff_ens3(14)=0.
                do ne=1,ens4
                  xff_ens3(14)=xff_ens3(14)-omeg(i,k22(i),ne)/(fens4*9.81)
                enddo
                if(xff_ens3(14).lt.0.)xff_ens3(14)=0.
                xff_ens3(5)=0.
                do ne=1,ens4
                  xff_ens3(5)=xff_ens3(5)-omeg(i,kbcon(i),ne)/(fens4*9.81)
                enddo
                if(xff_ens3(5).lt.0.)xff_ens3(5)=0.
!
! minimum below kbcon
!
                if(high_resolution.eq.0)then
                   xff_ens3(4)=-omeg(i,2,1)/9.81
                   do k=2,kbcon(i)-1
                   do ne=1,ens4
                     xomg=-omeg(i,k,ne)/9.81
                     if(xomg.lt.xff_ens3(4))xff_ens3(4)=xomg
                   enddo
                   enddo
                   if(xff_ens3(4).lt.0.)xff_ens3(4)=0.
!
! max below kbcon
                   xff_ens3(6)=-omeg(i,2,1)/9.81
                   do k=2,kbcon(i)-1
                   do ne=1,ens4
                     xomg=-omeg(i,k,ne)/9.81
                     if(xomg.gt.xff_ens3(6))xff_ens3(6)=xomg
                   enddo
                   enddo
                   if(xff_ens3(6).lt.0.)xff_ens3(6)=0.
                endif
                if(high_resolution.eq.1)then
                   xff_ens3(5)=min(xff_ens3(5),xff_ens3(14))
                   xff_ens3(4)=xff_ens3(5)
                   xff_ens3(6)=xff_ens3(5)
                endif
!
!--- more like Krishnamurti et al.; pick max and average values
!
                xff_ens3(7)=mconv(i,1)
                xff_ens3(8)=mconv(i,1)
                xff_ens3(9)=mconv(i,1)
                if(ens4.gt.1)then
                   do ne=2,ens4
                      if (mconv(i,ne).gt.xff_ens3(7))xff_ens3(7)=mconv(i,ne)
                   enddo
                   do ne=2,ens4
                      if (mconv(i,ne).lt.xff_ens3(8))xff_ens3(8)=mconv(i,ne)
                   enddo
                   do ne=2,ens4
                      xff_ens3(9)=xff_ens3(9)+mconv(i,ne)
                   enddo
                   xff_ens3(9)=xff_ens3(9)/fens4
                endif
                if(high_resolution.eq.1)then
                   xff_ens3(7)=xff_ens3(9)
                   xff_ens3(8)=xff_ens3(9)
                   xff_ens3(15)=xff_ens3(9)
                endif
!
                if(high_resolution.eq.0)then
                if(irandom.eq.1)then
                   seed(1)=i
                   call random_seed (PUT=seed)
                   call random_number (xxx)
                   ixxx=min(ens4,max(1,int(fens4*xxx+1.e-8)))
                   xff_ens3(15)=mconv(i,ixxx)
                else
                   xff_ens3(15)=mconv(i,1)
                endif
                endif
!
!--- more like Fritsch Chappel or Kain Fritsch (plus triggers)
!
                xff_ens3(10)=A_AVE/(60.*40.)
                xff_ens3(11)=AA1(I)/(60.*40.)
                if(irandom.eq.1)then
                   seed(1)=i
                   call random_seed (PUT=seed)
                   call random_number (xxx)
                   ixxx=min(ens4,max(1,int(fens4*xxx+1.e-8)))
                   xff_ens3(12)=AXX(I,ixxx)/(60.*40.)
                else
                   xff_ens3(12)=AA1(I)/(60.*40.)
                endif
                if(high_resolution.eq.1)then
                   xff_ens3(11)=xff_ens3(10)
                   xff_ens3(12)=xff_ens3(10)
                endif
                !write(0,*)'1',xff_ens3(1:16)
!
!--- more original Arakawa-Schubert (climatologic value of aa0)
!
                if(icoic.eq.0)then
                if(xff0.lt.0.)then
                     xff_ens3(1)=0.
                     xff_ens3(2)=0.
                     xff_ens3(3)=0.
                     xff_ens3(13)=0.
                     xff_ens3(10)=0.
                     xff_ens3(11)=0.
                     xff_ens3(12)=0.
                endif
                endif


                !write(0,*)'2',xff_ens3(1:16)

                do nens=1,maxens
                   XK(nens)=(XAA0(I,nens)-AA1(I))/MBDT(2)
                   if(xk(nens).le.0.and.xk(nens).gt.-1.e-6) &
                           xk(nens)=-1.e-6
                   if(xk(nens).gt.0.and.xk(nens).lt.1.e-6) &
                           xk(nens)=1.e-6
                !write(0,*)'3',xk(nens),XAA0(I,nens)
                enddo
!
!--- add up all ensembles
!
                do 350 ne=1,maxens
!
!--- for every xk, we have maxens3 xffs
!--- iens is from outermost ensemble (most expensive!
!
!--- iedt (maxens2 belongs to it)
!--- is from second, next outermost, not so expensive
!
!--- so, for every outermost loop, we have maxens*maxens2*3
!--- ensembles!!! nall would be 0, if everything is on first
!--- loop index, then ne would start counting, then iedt, then iens....
!
                   iresult=0
                   iresultd=0
                   iresulte=0
                   nall=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3 &
                        +(ne-1)*maxens3
!
! over water, enfor!e small cap for some of the closures
!
                if(xland(i).lt.0.1)then
                 if(ierr2(i).gt.0.or.ierr3(i).gt.0)then
                      xff_ens3(1) =0.
                      massfln(i,j,nall+1)=0.
                      xff_ens3(2) =0.
                      massfln(i,j,nall+2)=0.
                      xff_ens3(3) =0.
                      massfln(i,j,nall+3)=0.
                      xff_ens3(10) =0.
                      massfln(i,j,nall+10)=0.
                      xff_ens3(11) =0.
                      massfln(i,j,nall+11)=0.
                      xff_ens3(12) =0.
                      massfln(i,j,nall+12)=0.
                      xff_ens3(7) =0.
                      massfln(i,j,nall+7)=0.
                      xff_ens3(8) =0.
                      massfln(i,j,nall+8)=0.
                      xff_ens3(9) =0.
                      massfln(i,j,nall+9)=0.
                      xff_ens3(13) =0.
                      massfln(i,j,nall+13)=0.
                      xff_ens3(15) =0.
                      massfln(i,j,nall+15)=0.
                endif
                endif
!
! end water treatment
!
!
!--- check for upwind convection
!                  iresult=0
                   massfld=0.

!                  call cup_direction2(i,j,dir,iact_old_gr, &
!                       massflx,iresult,1,                  &
!                       massfld,                            &
!                       itf,jtf,ktf,          &
!                       ims,ime, jms,jme, kms,kme,          &
!                       its,ite, jts,jte, kts,kte          )
!                  if(i.eq.ipr.and.j.eq.jpr.and.iedt.eq.1.and.ne.eq.1)then
!                  if(iedt.eq.1.and.ne.eq.1)then
!                   print *,massfld,ne,iedt,iens
!                   print *,xk(ne),xff_ens3(1),xff_ens3(2),xff_ens3(3)
!                  endif
!                  print *,i,j,massfld,aa0(i),aa1(i)
                   IF(XK(ne).lt.0.and.xff0.gt.0.)iresultd=1
                   iresulte=max(iresult,iresultd)
                   iresulte=1
                   if(iresulte.eq.1)then
!
!--- special treatment for stability closures
!

                      if(xff0.ge.0.)then
                         xf_ens(i,j,nall+1)=massfld
                         xf_ens(i,j,nall+2)=massfld
                         xf_ens(i,j,nall+3)=massfld
                         xf_ens(i,j,nall+13)=massfld
                         if(xff_ens3(1).gt.0)xf_ens(i,j,nall+1)=max(0.,-xff_ens3(1)/xk(ne)) &
                                        +massfld
                         if(xff_ens3(2).gt.0)xf_ens(i,j,nall+2)=max(0.,-xff_ens3(2)/xk(ne)) &
                                        +massfld
                         if(xff_ens3(3).gt.0)xf_ens(i,j,nall+3)=max(0.,-xff_ens3(3)/xk(ne)) &
                                        +massfld
                         if(xff_ens3(13).gt.0)xf_ens(i,j,nall+13)=max(0.,-xff_ens3(13)/xk(ne)) &
                                        +massfld
!                       endif
                      else
                         xf_ens(i,j,nall+1)=massfld
                         xf_ens(i,j,nall+2)=massfld
                         xf_ens(i,j,nall+3)=massfld
                         xf_ens(i,j,nall+13)=massfld
                      endif
!
!--- if iresult.eq.1, following independent of xff0
!
                         xf_ens(i,j,nall+4)=max(0.,xff_ens3(4) &
                            +massfld)
                         xf_ens(i,j,nall+5)=max(0.,xff_ens3(5) &
                                        +massfld)
                         xf_ens(i,j,nall+6)=max(0.,xff_ens3(6) &
                                        +massfld)
                         xf_ens(i,j,nall+14)=max(0.,xff_ens3(14) &
                                        +massfld)
                         a1=max(1.e-3,pr_ens(i,j,nall+7))
                         xf_ens(i,j,nall+7)=max(0.,xff_ens3(7) &
                                     /a1)
                         a1=max(1.e-3,pr_ens(i,j,nall+8))
                         xf_ens(i,j,nall+8)=max(0.,xff_ens3(8) &
                                     /a1)
                         a1=max(1.e-3,pr_ens(i,j,nall+9))
                         xf_ens(i,j,nall+9)=max(0.,xff_ens3(9) &
                                     /a1)
                         a1=max(1.e-3,pr_ens(i,j,nall+15))
                         xf_ens(i,j,nall+15)=max(0.,xff_ens3(15) &
                                     /a1)
                         if(XK(ne).lt.0.)then
                            xf_ens(i,j,nall+10)=max(0., &
                                        -xff_ens3(10)/xk(ne)) &
                                        +massfld
                            xf_ens(i,j,nall+11)=max(0., &
                                        -xff_ens3(11)/xk(ne)) &
                                        +massfld
                            xf_ens(i,j,nall+12)=max(0., &
                                        -xff_ens3(12)/xk(ne)) &
                                        +massfld
                         else
                            xf_ens(i,j,nall+10)=massfld
                            xf_ens(i,j,nall+11)=massfld
                            xf_ens(i,j,nall+12)=massfld
                         endif
                      !if(nall.eq.48)write(0,*)nall,xf_ens(i,j,nall+1:nall+16)
                      if(icoic.ge.1)then
                      closure_n(i)=0.
                      xf_ens(i,j,nall+1)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+2)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+3)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+4)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+5)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+6)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+7)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+8)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+9)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+10)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+11)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+12)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+13)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+14)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+15)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+16)=xf_ens(i,j,nall+icoic)
                      endif
!
! 16 is a randon pick from the oher 15
!
                if(irandom.eq.1)then
                   call random_number (xxx)
                   ixxx=min(15,max(1,int(15.*xxx+1.e-8)))
                   xf_ens(i,j,nall+16)=xf_ens(i,j,nall+ixxx)
                else
                   xf_ens(i,j,nall+16)=xf_ens(i,j,nall+1)
                endif
!
!
!--- store new for next time step
!
                      do nens3=1,maxens3
                        massfln(i,j,nall+nens3)=edt(i) &
                                                *xf_ens(i,j,nall+nens3)
                        massfln(i,j,nall+nens3)=max(0., &
                                              massfln(i,j,nall+nens3))
                      enddo
!
!
!--- do some more on the caps!!! ne=1 for 175, ne=2 for 100,....
!
!     do not care for caps here for closure groups 1 and 5,
!     they are fine, do not turn them off here
!
!
                if(ne.eq.2.and.ierr2(i).gt.0)then
                      !write(0,*)'in forcing, setting xf_ens to zero for ierr2 '
                      xf_ens(i,j,nall+1) =0.
                      xf_ens(i,j,nall+2) =0.
                      xf_ens(i,j,nall+3) =0.
                      xf_ens(i,j,nall+4) =0.
                      xf_ens(i,j,nall+5) =0.
                      xf_ens(i,j,nall+6) =0.
                      xf_ens(i,j,nall+7) =0.
                      xf_ens(i,j,nall+8) =0.
                      xf_ens(i,j,nall+9) =0.
                      xf_ens(i,j,nall+10)=0.
                      xf_ens(i,j,nall+11)=0.
                      xf_ens(i,j,nall+12)=0.
                      xf_ens(i,j,nall+13)=0.
                      xf_ens(i,j,nall+14)=0.
                      xf_ens(i,j,nall+15)=0.
                      xf_ens(i,j,nall+16)=0.
                      massfln(i,j,nall+1)=0.
                      massfln(i,j,nall+2)=0.
                      massfln(i,j,nall+3)=0.
                      massfln(i,j,nall+4)=0.
                      massfln(i,j,nall+5)=0.
                      massfln(i,j,nall+6)=0.
                      massfln(i,j,nall+7)=0.
                      massfln(i,j,nall+8)=0.
                      massfln(i,j,nall+9)=0.
                      massfln(i,j,nall+10)=0.
                      massfln(i,j,nall+11)=0.
                      massfln(i,j,nall+12)=0.
                      massfln(i,j,nall+13)=0.
                      massfln(i,j,nall+14)=0.
                      massfln(i,j,nall+15)=0.
                      massfln(i,j,nall+16)=0.
                endif
                if(ne.eq.3.and.ierr3(i).gt.0)then
                      !write(0,*)'in forcing, setting xf_ens to zero for ierr3 '
                      xf_ens(i,j,nall+1) =0.
                      xf_ens(i,j,nall+2) =0.
                      xf_ens(i,j,nall+3) =0.
                      xf_ens(i,j,nall+4) =0.
                      xf_ens(i,j,nall+5) =0.
                      xf_ens(i,j,nall+6) =0.
                      xf_ens(i,j,nall+7) =0.
                      xf_ens(i,j,nall+8) =0.
                      xf_ens(i,j,nall+9) =0.
                      xf_ens(i,j,nall+10)=0.
                      xf_ens(i,j,nall+11)=0.
                      xf_ens(i,j,nall+12)=0.
                      xf_ens(i,j,nall+13)=0.
                      xf_ens(i,j,nall+14)=0.
                      xf_ens(i,j,nall+15)=0.
                      xf_ens(i,j,nall+16)=0.
                      massfln(i,j,nall+1)=0.
                      massfln(i,j,nall+2)=0.
                      massfln(i,j,nall+3)=0.
                      massfln(i,j,nall+4)=0.
                      massfln(i,j,nall+5)=0.
                      massfln(i,j,nall+6)=0.
                      massfln(i,j,nall+7)=0.
                      massfln(i,j,nall+8)=0.
                      massfln(i,j,nall+9)=0.
                      massfln(i,j,nall+10)=0.
                      massfln(i,j,nall+11)=0.
                      massfln(i,j,nall+12)=0.
                      massfln(i,j,nall+13)=0.
                      massfln(i,j,nall+14)=0.
                      massfln(i,j,nall+15)=0.
                      massfln(i,j,nall+16)=0.
                endif

                   endif
 350            continue
! ne=1, cap=175
!
                   nall=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3
! ne=2, cap=100
!
                   nall2=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3 &
                        +(2-1)*maxens3
                      xf_ens(i,j,nall+4) = xf_ens(i,j,nall2+4)
                      xf_ens(i,j,nall+5) =xf_ens(i,j,nall2+5)
                      xf_ens(i,j,nall+6) =xf_ens(i,j,nall2+6)
                      xf_ens(i,j,nall+14) =xf_ens(i,j,nall2+14)
                      xf_ens(i,j,nall+7) =xf_ens(i,j,nall2+7)
                      xf_ens(i,j,nall+8) =xf_ens(i,j,nall2+8)
                      xf_ens(i,j,nall+9) =xf_ens(i,j,nall2+9)
                      xf_ens(i,j,nall+15) =xf_ens(i,j,nall2+15)
                      xf_ens(i,j,nall+10)=xf_ens(i,j,nall2+10)
                      xf_ens(i,j,nall+11)=xf_ens(i,j,nall2+11)
                      xf_ens(i,j,nall+12)=xf_ens(i,j,nall2+12)
                go to 100
             endif
          elseif(ierr(i).ne.20.and.ierr(i).ne.0)then
             do n=1,ensdim
               xf_ens(i,j,n)=0.
               massfln(i,j,n)=0.
             enddo
          endif
 100   continue
      nall=48
      !write(0,*)nall,xf_ens(1,1,nall+1:nall+16)

   END SUBROUTINE cup_forcing_ens_3d


   SUBROUTINE cup_kbcon(ierrc,cap_inc,iloop,k22,kbcon,he_cup,hes_cup, &
              ierr,kbmax,p_cup,cap_max,                         &
              itf,jtf,ktf,                        &
              its,ite, jts,jte, kts,kte                        )

   IMPLICIT NONE
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
  !
  !
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        he_cup,hes_cup,p_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        cap_max,cap_inc
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbmax
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        kbcon,k22,ierr
     integer                                                           &
        ,intent (in   )                   ::                           &
        iloop
     character*50 :: ierrc(its:ite)
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
     real                                 ::                           &
        pbcdif,plus,hetest
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!
       DO 27 i=its,itf
      kbcon(i)=1
      IF(ierr(I).ne.0)GO TO 27
      KBCON(I)=K22(I)
      GO TO 32
 31   CONTINUE
      KBCON(I)=KBCON(I)+1
      IF(KBCON(I).GT.KBMAX(i)+2)THEN
         if(iloop.ne.4)ierr(i)=3
         if(iloop.lt.4)ierrc(i)="could not find reasonable kbcon in cup_kbcon"
!        if(iloop.lt.4)ierr(i)=997
         !write(0,*)'kbcon,kbmax,k22 = ',KBCON(I),KBMAX(i),K22(I)
         !write(0,*)HE_cup(I,K22(I)),HES_cup(I,KBCON(I)),-P_cup(I,KBCON(I))+P_cup(I,K22(I))
        GO TO 27
      ENDIF
 32   CONTINUE
      hetest=HE_cup(I,K22(I))
      if(iloop.eq.5)then
       do k=1,k22(i)
         hetest=max(hetest,he_cup(i,k))
       enddo
      endif
      IF(HETEST.LT.HES_cup(I,KBCON(I)))then
        !write(0,*)'htest',k22(i),kbcon(i),HETEST,-P_cup(I,KBCON(I))+P_cup(I,K22(I))
        GO TO 31
      endif

!     cloud base pressure and max moist static energy pressure
!     i.e., the depth (in mb) of the layer of negative buoyancy
      if(KBCON(I)-K22(I).eq.1)go to 27
      if(iloop.eq.5 .and. (KBCON(I)-K22(I)).eq.0)go to 27
      PBCDIF=-P_cup(I,KBCON(I))+P_cup(I,K22(I))
      plus=max(25.,cap_max(i)-float(iloop-1)*cap_inc(i))
      if(iloop.eq.4)plus=cap_max(i)
!
! for shallow convection, if cap_max is greater than 25, it is the pressure at pbltop
      if(iloop.eq.5)plus=25.
      if(iloop.eq.5.and.cap_max(i).gt.25)pbcdif=-P_cup(I,KBCON(I))+cap_max(i)
      IF(PBCDIF.GT.plus)THEN
        !write(0,*)'htest',k22(i),kbcon(i),plus,-P_cup(I,KBCON(I))+P_cup(I,K22(I))
        K22(I)=K22(I)+1
        KBCON(I)=K22(I)
        GO TO 32
      ENDIF
 27   CONTINUE

   END SUBROUTINE cup_kbcon


   SUBROUTINE cup_ktop(ierrc,ilo,dby,kbcon,ktop,ierr,              &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
  ! dby = buoancy term
  ! ktop = cloud top (output)
  ! ilo  = flag
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
        dby
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon
     integer                                                           &
        ,intent (in   )                   ::                           &
        ilo
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        ktop
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     character*50 :: ierrc(its:ite)
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
!
        DO 42 i=its,itf
        ktop(i)=1
         IF(ierr(I).EQ.0)then
          DO 40 K=KBCON(I)+1,ktf-1
            IF(DBY(I,K).LE.0.)THEN
                KTOP(I)=K-1
                GO TO 41
             ENDIF
  40      CONTINUE
          if(ilo.eq.1)ierr(i)=5
          if(ilo.eq.1)ierrc(i)="problem with defining ktop"
!         if(ilo.eq.2)ierr(i)=998
          GO TO 42
  41     CONTINUE
         do k=ktop(i)+1,ktf
           dby(i,k)=0.
         enddo
         if(kbcon(i).eq.ktop(i))then
            ierr(i)=55
            ierrc(i)="kbcon == ktop "
         endif
         endif
  42     CONTINUE

   END SUBROUTINE cup_ktop


   SUBROUTINE cup_MAXIMI(ARRAY,KS,KE,MAXX,ierr,              &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         itf,jtf,ktf,                                    &
         its,ite, jts,jte, kts,kte
  ! array input array
  ! x output array with return values
  ! kt output array of levels
  ! ks,kend  check-range
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
         array
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
         ierr,ke
     integer                                                           &
        ,intent (in   )                   ::                           &
         ks
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
         maxx
     real,    dimension (its:ite)         ::                           &
         x
     real                                 ::                           &
         xar
     integer                              ::                           &
         i,k

       DO 200 i=its,itf
       MAXX(I)=KS
       if(ierr(i).eq.0)then
      X(I)=ARRAY(I,KS)
!
       DO 100 K=KS,KE(i)
         XAR=ARRAY(I,K)
         IF(XAR.GE.X(I)) THEN
            X(I)=XAR
            MAXX(I)=K
         ENDIF
 100  CONTINUE
      endif
 200  CONTINUE

   END SUBROUTINE cup_MAXIMI


   SUBROUTINE cup_minimi(ARRAY,KS,KEND,KT,ierr,              &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         itf,jtf,ktf,                                    &
         its,ite, jts,jte, kts,kte
  ! array input array
  ! x output array with return values
  ! kt output array of levels
  ! ks,kend  check-range
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
         array
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
         ierr,ks,kend
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
         kt
     real,    dimension (its:ite)         ::                           &
         x
     integer                              ::                           &
         i,k,kstop

       DO 200 i=its,itf
      KT(I)=KS(I)
      if(ierr(i).eq.0)then
      X(I)=ARRAY(I,KS(I))
       KSTOP=MAX(KS(I)+1,KEND(I))
!
       DO 100 K=KS(I)+1,KSTOP
         IF(ARRAY(I,K).LT.X(I)) THEN
              X(I)=ARRAY(I,K)
              KT(I)=K
         ENDIF
 100  CONTINUE
      endif
 200  CONTINUE

   END SUBROUTINE cup_MINIMI



   SUBROUTINE cup_up_aa0(aa0,z,zu,dby,GAMMA_CUP,t_cup,       &
              kbcon,ktop,ierr,                               &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,                                     &
        its,ite, jts,jte, kts,kte
  ! aa0 cloud work function
  ! gamma_cup = gamma on model cloud levels
  ! t_cup = temperature (Kelvin) on model cloud levels
  ! dby = buoancy term
  ! zu= normalized updraft mass flux
  ! z = heights of model levels
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        z,zu,gamma_cup,t_cup,dby
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop
!
! input and output
!


     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        aa0
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
     real                                 ::                           &
        dz,da
!
        do i=its,itf
         aa0(i)=0.
        enddo
        DO 100 k=kts+1,ktf
        DO 100 i=its,itf
         IF(ierr(i).ne.0)GO TO 100
         IF(K.LE.KBCON(I))GO TO 100
         IF(K.Gt.KTOP(I))GO TO 100
         DZ=Z(I,K)-Z(I,K-1)
         da=zu(i,k)*DZ*(9.81/(1004.*( &
                (T_cup(I,K)))))*DBY(I,K-1)/ &
             (1.+GAMMA_CUP(I,K))
         IF(K.eq.KTOP(I).and.da.le.0.)go to 100
         AA0(I)=AA0(I)+da
         if(aa0(i).lt.0.)aa0(i)=0.
100     continue

   END SUBROUTINE cup_up_aa0


   SUBROUTINE cup_up_he(k22,hkb,z_cup,cd,entr,he_cup,hc,     &
              kbcon,ierr,dby,he,hes_cup,name,                &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  itf,jtf,ktf,           &
                                  its,ite, jts,jte, kts,kte
      character *(*), intent (in)        ::                           &
       name
  ! hc = cloud moist static energy
  ! hkb = moist static energy at originating level
  ! he = moist static energy on model levels
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function
  ! z_cup = heights of model cloud levels
  ! entr = entrainment rate
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        he,he_cup,hes_cup,z_cup,cd
  ! entr= entrainment rate
     real                                                              &
        ,intent (in   )                   ::                           &
        entr
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,k22
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        hc,dby
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        hkb
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
     real                                 ::                           &
        dz
!
!--- moist static energy inside cloud
!
      do k=kts,ktf
      do i=its,itf
         hc(i,k)=0.
         DBY(I,K)=0.
      enddo
      enddo
      do i=its,itf
         hkb(i)=0.
      enddo
      do i=its,itf
        if(ierr(I).eq.0)then
          hkb(i)=he_cup(i,k22(i))
          if(name.eq.'shallow')then
             do k=1,k22(i)
               hkb(i)=max(hkb(i),he_cup(i,k))
             enddo
          endif
          do k=1,k22(i)
              hc(i,k)=he_cup(i,k)
          enddo
          do k=k22(i),kbcon(i)-1
              hc(i,k)=hkb(i)
          enddo
          k=kbcon(i)
          hc(i,k)=hkb(i)
          DBY(I,Kbcon(i))=Hkb(I)-HES_cup(I,K)
        endif
      enddo
      do k=kts+1,ktf
      do i=its,itf
        if(k.gt.kbcon(i).and.ierr(I).eq.0)then
           DZ=Z_cup(i,K)-Z_cup(i,K-1)
           HC(i,K)=(HC(i,K-1)*(1.-.5*CD(i,K-1)*DZ)+entr* &
                DZ*HE(i,K-1))/(1.+entr*DZ-.5*cd(i,k-1)*dz)
           DBY(I,K)=HC(I,K)-HES_cup(I,K)
        endif
      enddo

      enddo
   END SUBROUTINE cup_up_he


   SUBROUTINE cup_up_moisture(name,ierr,z_cup,qc,qrc,pw,pwav,     &
              ccnclean,p_cup,kbcon,ktop,cd,dby,mentr_rate,clw_all,&
              t_cup,q,GAMMA_cup,zu,qes_cup,k22,qe_cup,xl,ccn,rho, &
              up_massentr,up_massdetr,psum,psumh,                 &
              autoconv,aeroevap,itest,itf,jtf,ktf,                &
              its,ite, jts,jte, kts,kte                     )

  IMPLICIT NONE
  real, parameter :: BDISPM = 0.366       !Berry--size dispersion (maritime)
  REAL, PARAMETER :: BDISPC = 0.146       !Berry--size dispersion (continental)
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  itest,autoconv,aeroevap,itf,jtf,ktf,           &
                                  its,ite, jts,jte, kts,kte
  ! cd= detrainment function
  ! q = environmental q on model levels
  ! qe_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function
  ! zu = normalized updraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! mentr_rate = entrainment rate
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        t_cup,p_cup,rho,q,zu,gamma_cup,qe_cup,                         &
        up_massentr,up_massdetr,dby,qes_cup,z_cup,cd
  ! entr= entrainment rate
     real                                                              &
        ,intent (in   )                   ::                           &
        ccnclean,mentr_rate,xl
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop,k22
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
      character *(*), intent (in)        ::                           &
       name
   ! qc = cloud q (including liquid water) after entrainment
   ! qrch = saturation q in cloud
   ! qrc = liquid water content in cloud after rainout
   ! pw = condensate that will fall out at that level
   ! pwav = totan normalized integrated condensate (I1)
   ! c0 = conversion rate (cloud to rain)

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qc,qrc,pw,clw_all
     real,    dimension (its:ite,kts:kte) ::                           &
        qch,qrcb,pwh,clw_allh
     real,    dimension (its:ite)         ::                           &
        pwavh
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pwav,psum,psumh
     real,    dimension (its:ite)                                      &
        ,intent (in  )                   ::                           &
        ccn
!
!  local variables in this routine
!

     integer                              ::                           &
        iounit,iprop,iall,i,k
     real                                 ::                           &
        prop_ave,qrcb_h,bdsp,dp,g,rhoc,dh,qrch,c0,dz,radius,berryc0,q1,berryc
     real,    dimension (kts:kte)         ::                           &
        prop_b
!
        prop_b(kts:kte)=0
        iall=0
        c0=.002
        g=9.81
        bdsp=BDISPM
!
!--- no precip for small clouds
!
        if(name.eq.'shallow')c0=0.
        do i=its,itf
          pwav(i)=0.
          pwavh(i)=0.
          psum(i)=0.
          psumh(i)=0.
        enddo
        do k=kts,ktf
        do i=its,itf
          pw(i,k)=0.
          pwh(i,k)=0.
          qc(i,k)=0.
          if(ierr(i).eq.0)qc(i,k)=qes_cup(i,k)
          if(ierr(i).eq.0)qch(i,k)=qes_cup(i,k)
          clw_all(i,k)=0.
          clw_allh(i,k)=0.
          qrc(i,k)=0.
          qrcb(i,k)=0.
        enddo
        enddo
      do i=its,itf
      if(ierr(I).eq.0)then
      do k=k22(i),kbcon(i)-1
        qc(i,k)=qe_cup(i,k22(i))
        qch(i,k)=qe_cup(i,k22(i))
        qc(i,k)=qe_cup(i,k)
        qch(i,k)=qe_cup(i,k)
      enddo
      endif
      enddo

        DO 100 k=kts+1,ktf
        DO 100 i=its,itf
         IF(ierr(i).ne.0)GO TO 100
         IF(K.Lt.KBCON(I))GO TO 100
         IF(K.Gt.KTOP(I))GO TO 100
         rhoc=.5*(rho(i,k)+rho(i,k+1))
         DZ=Z_cup(i,K)-Z_cup(i,K-1)
         DP=p_cup(i,K)-p_cup(i,K-1)
!
!--- saturation  in cloud, this is what is allowed to be in it
!
         QRCH=QES_cup(I,K)+(1./XL)*(GAMMA_cup(i,k) &
              /(1.+GAMMA_cup(i,k)))*DBY(I,K)
!
!------    1. steady state plume equation, for what could
!------       be in cloud without condensation
!
!
       qc(i,k)=(qc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k)*qc(i,k-1)+ &
                         up_massentr(i,k)*q(i,k-1))   /            &
                         (zu(i,k-1)-.5*up_massdetr(i,k)+up_massentr(i,k))
       qch(i,k)=(qch(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k)*qch(i,k-1)+ &
                         up_massentr(i,k)*q(i,k-1))   /            &
                         (zu(i,k-1)-.5*up_massdetr(i,k)+up_massentr(i,k))

        if(qc(i,k).le.qrch)go to 100
        if(qch(i,k).le.qrch)go to 100
!
!--- saturation  in cloud, this is what is allowed to be in it
!
!        QRCH=QES_cup(I,K)+(1./XL)*(GAMMA_cup(i,k) &
!             /(1.+GAMMA_cup(i,k)))*DBY(I,K)
!
!------- LIQUID WATER CONTENT IN cloud after rainout
!
        clw_all(i,k)=QC(I,K)-QRCH
        QRC(I,K)=(QC(I,K)-QRCH) ! /(1.+C0*DZ*zu(i,k))
        clw_allh(i,k)=QCH(I,K)-QRCH
        QRCB(I,K)=(QCH(I,K)-QRCH) ! /(1.+C0*DZ*zu(i,k))
    IF(autoconv.eq.2) then


!
! normalized berry
!
! first calculate for average conditions, used in cup_dd_edt!
! this will also determine proportionality constant prop_b, which, if applied,
! would give the same results as c0 under these conditions
!
!  do not dot his for ice only conversion, in that case use c0
!
!     if(t_cup(i,k).gt.263.)then
         q1=1.e3*rhoc*qrcb(i,k)  ! g/m^3 ! g[h2o]/cm^3
         berryc0=q1*q1/(60.0*(5.0 + 0.0366*CCNclean/ &
                ( q1 * BDSP)  ) ) !/(
         qrcb_h=qrcb(i,k)/(1.+c0*dz*zu(i,k))
         prop_b(k)=c0*qrcb_h*zu(i,k)/(1.e-3*berryc0)
!        berryc=berryc0*1.e-3/c0/qrcb(i,k)/zu(i,k)
!        prop_b(k)=1./berryc
!        pwh(i,k)=zu(i,k)*1.e-3*berryc0*dz*prop_b(k) ! 2.
!        if(t_cup(i,k).lt.263.)prop_b(k)=prop_b(k-1)
         pwh(i,k)=1.e-3*berryc0*dz*prop_b(k) ! 2.
         !write(11,*)'prop_b,prop_c = ',k,prop_b(k),berryc0,pwh(i,k)
         QRCb(I,K) = qrcb(i,k) - pwh(i,k)
         if(qrcb(i,k).lt.0.)then
           pwh(i,k)=qrcb(i,k)
           qrcb(i,k)=0.
         endif
!     else
!        if(t_cup(i,k).lt.263.)prop_b(k)=prop_b(k-1)
!        QRCB(I,K)=(QCH(I,K)-QRCH)/(1.+C0*DZ*zu(i,k))
!        PWH(i,k)=c0*dz*QRCB(I,K)*zu(i,k)
!        QRCb(I,K) = qrcb(i,k) - pwh(i,k)
!         if(qrcb(i,k).lt.0.)then
!           qrcb(i,k)=0.
!         endif
!     endif
      QCh(I,K)=QRCb(I,K)+qrch
      PWAVH(I)=PWAVH(I)+pwh(I,K)
      Psumh(I)=Psumh(I)+clw_allh(I,K)*zu(i,k) *dz
!
! then the real berry
!
!      if(t_cup(i,k).gt.263.)then
          q1=1.e3*rhoc*qrc(i,k)  ! g/m^3 ! g[h2o]/cm^3
          berryc0=q1*q1/(60.0*(5.0 + 0.0366*CCN(i)/ &
                ( q1 * BDSP)  ) ) !/(
!         berryc0=zu(i,k)*1.e-3*berryc0*dz*prop_b(k) ! 2.
          berryc0=1.e-3*berryc0*dz*prop_b(k) ! 2.
         !write(11,*)'prop_b = ',k,prop_b(k),qrc(i,k),berryc0
        !- condensed water remained in cloud
          berryc=qrc(i,k)
          QRC(I,K) = qrc(i,k) - berryc0 ! *zu(i,k)
          if(qrc(i,k).lt.0.)then
            !write(11,*) 'pw was very big ! '
            berryc0=berryc
            qrc(i,k)=0.
          endif
          pw(i,k)=berryc0
          if(qrc(i,k).lt.0.)then
            qrc(i,k)=0.
          endif
!       else
!          QRC(I,K)=(QC(I,K)-QRCH)/(1.+C0*DZ*zu(i,k))
!          PW(i,k)=c0*dz*QRC(I,K)*zu(i,k)
!       endif

!
!
!  if not running with berry at all, do the following
!
       ELSE       !c0=.002
!        QRC(I,K)=(QC(I,K)-QRCH)/(1.+C0*DZ*zu(i,k))
         QRC(I,K)=(QC(I,K)-QRCH)/(1.+C0*DZ)
         PW(i,k)=c0*dz*QRC(I,K)*zu(i,k)


!
!-------   3.Condensation
!
        if(iall.eq.1)then
          qrc(i,k)=0.
          pw(i,k)=(QC(I,K)-QRCH)*zu(i,k)
          if(pw(i,k).lt.0.)pw(i,k)=0.
        endif
      endif !autoconv
      iounit=22+(autoconv-1)*10+(aeroevap-1)*10
         !if(itest.eq.1)write(iounit,*)p_cup(i,k),pw(i,k)*rhoc

!
!----- set next level
!
         QC(I,K)=QRC(I,K)+qrch
!
!--- integrated normalized ondensate
!
         PWAV(I)=PWAV(I)+PW(I,K)
         Psum(I)=Psum(I)+clw_all(I,K)*zu(i,k) *dz
 100     CONTINUE
       prop_ave=0.
       iprop=0
       do k=kts,kte
        prop_ave=prop_ave+prop_b(k)
        if(prop_b(k).gt.0)iprop=iprop+1
       enddo
       iprop=max(iprop,1)
       !write(11,*)'prop_ave = ',prop_ave/float(iprop)
       !print *,'pwav = ',pwav(1)

   END SUBROUTINE cup_up_moisture


   SUBROUTINE cup_up_nms(itest,zu,z_cup,entr,cd,kbcon,ktop,ierr,k22,  &
              itf,jtf,ktf,                        &
              its,ite, jts,jte, kts,kte                        )

   IMPLICIT NONE

!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         itf,jtf,ktf,itest,                                    &
         its,ite, jts,jte, kts,kte
  ! cd= detrainment function
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
         z_cup,cd
  ! entr= entrainment rate
     real                                                              &
        ,intent (in   )                   ::                           &
         entr
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
         kbcon,ktop,k22
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
         ierr
   ! zu is the normalized mass flux

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
         zu
!
!  local variables in this routine
!

     integer                              ::                           &
         i,k
     real                                 ::                           &
         dz
!
!   initialize for this go around
!
       do k=kts,ktf
       do i=its,itf
         zu(i,k)=0.
       enddo
       enddo
!
! do normalized mass budget
!
       do i=its,itf
          IF(ierr(I).eq.0)then
             do k=k22(i),kbcon(i)
               zu(i,k)=1.
             enddo
             DO K=KBcon(i)+1,KTOP(i)
               DZ=Z_cup(i,K)-Z_cup(i,K-1)
               ZU(i,K)=ZU(i,K-1)*(1.+(entr-cd(i,k-1))*DZ)
!              if(itest.eq.1)write(9,111)'entr, detr at level ',k,' = ', &
!                             entr*dz*zu(i,k-1),cd(i,k-1)*dz*zu(i,k-1), &
!                             cd(i,k-1),dz,zu(i,k-1)
!              if(itest.eq.1)write(9,111),k, &
!                             entr*dz*zu(i,k-1),cd(i,k-1)*dz*zu(i,k-1), &
!                             entr,cd(i,k-1),dz,zu(i,k-1)
             enddo
          endif
       enddo
!111 format(1x,20A,i2,3A,3E13.6,1F8.1,1x,1F7.4)
111 format(1x,i3,4E13.6,F8.1,1x,F7.4)
!111 format(1x,i2,3E13.6,1F8.1,1x,1F7.4)

   END SUBROUTINE cup_up_nms

!====================================================================
   SUBROUTINE g3init(RTHCUTEN,RQVCUTEN,RQCCUTEN,RQICUTEN,           &
                        MASS_FLUX,cp,restart,                       &
                        P_QC,P_QI,P_FIRST_SCALAR,                   &
                        RTHFTEN, RQVFTEN,                           &
                        APR_GR,APR_W,APR_MC,APR_ST,APR_AS,          &
                        APR_CAPMA,APR_CAPME,APR_CAPMI,              &
                        cugd_tten,cugd_ttens,cugd_qvten,            &
                        cugd_qvtens,cugd_qcten,                     &
                        allowed_to_read,                            &
                        ids, ide, jds, jde, kds, kde,               &
                        ims, ime, jms, jme, kms, kme,               &
                        its, ite, jts, jte, kts, kte               )
!--------------------------------------------------------------------
   IMPLICIT NONE
!--------------------------------------------------------------------
   LOGICAL , INTENT(IN)           ::  restart,allowed_to_read
   INTEGER , INTENT(IN)           ::  ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      its, ite, jts, jte, kts, kte
   INTEGER , INTENT(IN)           ::  P_FIRST_SCALAR, P_QI, P_QC
   REAL,     INTENT(IN)           ::  cp

   REAL,     DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) ::       &
                                                          CUGD_TTEN,         &
                                                          CUGD_TTENS,        &
                                                          CUGD_QVTEN,        &
                                                          CUGD_QVTENS,       &
                                                          CUGD_QCTEN
   REAL,     DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) ::       &
                                                          RTHCUTEN, &
                                                          RQVCUTEN, &
                                                          RQCCUTEN, &
                                                          RQICUTEN

   REAL,     DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) ::       &
                                                          RTHFTEN,  &
                                                          RQVFTEN

   REAL,     DIMENSION( ims:ime , jms:jme ) , INTENT(OUT) ::        &
                                APR_GR,APR_W,APR_MC,APR_ST,APR_AS,  &
                                APR_CAPMA,APR_CAPME,APR_CAPMI,      &
                                MASS_FLUX

   INTEGER :: i, j, k, itf, jtf, ktf

   jtf=min0(jte,jde-1)
   ktf=min0(kte,kde-1)
   itf=min0(ite,ide-1)

   IF(.not.restart)THEN
     DO j=jts,jte
     DO k=kts,kte
     DO i=its,ite
        RTHCUTEN(i,k,j)=0.
        RQVCUTEN(i,k,j)=0.
     ENDDO
     ENDDO
     ENDDO
     DO j=jts,jte
     DO k=kts,kte
     DO i=its,ite
       cugd_tten(i,k,j)=0.
       cugd_ttens(i,k,j)=0.
       cugd_qvten(i,k,j)=0.
       cugd_qvtens(i,k,j)=0.
     ENDDO
     ENDDO
     ENDDO

     DO j=jts,jtf
     DO k=kts,ktf
     DO i=its,itf
        RTHFTEN(i,k,j)=0.
        RQVFTEN(i,k,j)=0.
     ENDDO
     ENDDO
     ENDDO

     IF (P_QC .ge. P_FIRST_SCALAR) THEN
        DO j=jts,jtf
        DO k=kts,ktf
        DO i=its,itf
           RQCCUTEN(i,k,j)=0.
           cugd_qcten(i,k,j)=0.
        ENDDO
        ENDDO
        ENDDO
     ENDIF

     IF (P_QI .ge. P_FIRST_SCALAR) THEN
        DO j=jts,jtf
        DO k=kts,ktf
        DO i=its,itf
           RQICUTEN(i,k,j)=0.
        ENDDO
        ENDDO
        ENDDO
     ENDIF

     DO j=jts,jtf
     DO i=its,itf
        mass_flux(i,j)=0.
     ENDDO
     ENDDO

     DO j=jts,jtf
     DO i=its,itf
        APR_GR(i,j)=0.
        APR_ST(i,j)=0.
        APR_W(i,j)=0.
        APR_MC(i,j)=0.
        APR_AS(i,j)=0.
        APR_CAPMA(i,j)=0.
        APR_CAPME(i,j)=0.
        APR_CAPMI(i,j)=0.
     ENDDO
     ENDDO

   ENDIF

   END SUBROUTINE g3init


   SUBROUTINE cup_axx(ierrc,tcrit,kbmax,z1,p,psur,xl,rv,cp,tx,qx,axx,ierr,    &
           cap_max,cap_max_increment,entr_rate,mentr_rate,&
           j,itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte,ens4)
   IMPLICIT NONE
   INTEGER,      INTENT(IN   ) ::                                             &
                                  j,itf,jtf,ktf,                &
                                  its,ite, jts,jte, kts,kte,ens4
     real, dimension (its:ite,kts:kte,1:ens4)                                 &
         , intent(inout) ::                                                   &
           tx,qx
     real, dimension (its:ite,kts:kte)                                 &
         , intent(in) ::                                                   &
           p
     real, dimension (its:ite)                                 &
         , intent(in) ::                                                   &
           z1,psur,cap_max,cap_max_increment
     real, intent(in) ::                                                   &
           tcrit,xl,rv,cp,mentr_rate,entr_rate
     real, dimension (its:ite,1:ens4)                                 &
         , intent(out) ::                                                   &
           axx
     integer, dimension (its:ite), intent (in) ::                             &
           ierr,kbmax
     integer, dimension (its:ite) ::                             &
           ierrxx,k22xx,kbconxx,ktopxx,kstabm,kstabi
      real, dimension (1:2) :: AE,BE,HT
      real, dimension (its:ite,kts:kte) :: tv
      real :: e,tvbar
     integer n,i,k,iph
     real,    dimension (its:ite,kts:kte) ::                           &
        he,hes,qes,z,                                                  &
        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,      &
        tn_cup,                                                        &
        dby,qc,qrcd,pwd,pw,hcd,qcd,dbyd,hc,qrc,zu,zd,cd

     real,    dimension (its:ite) ::                                   &
       AA0,HKB,QKB,          &
       PWAV,BU
     character*50 :: ierrc(its:ite)
      do n=1,ens4
      do i=its,ite
       axx(i,n)=0.
      enddo
      enddo
     HT(1)=XL/CP
     HT(2)=2.834E6/CP
     BE(1)=.622*HT(1)/.286
     AE(1)=BE(1)/273.+ALOG(610.71)
     BE(2)=.622*HT(2)/.286
     AE(2)=BE(2)/273.+ALOG(610.71)
!
!
     do 100 n=1,ens4

      do k=kts,ktf
      do i=its,itf
        cd(i,k)=0.1*entr_rate
      enddo
      enddo


      do i=its,itf
        ierrxx(i)=ierr(i)
        k22xx(i)=1
        kbconxx(i)=1
        ktopxx(i)=1
        kstabm(i)=ktf-1
      enddo
      DO k=kts,ktf
      do i=its,itf
        if(ierrxx(i).eq.0)then
        IPH=1
        IF(Tx(I,K,n).LE.TCRIT)IPH=2
        E=EXP(AE(IPH)-BE(IPH)/TX(I,K,N))
        QES(I,K)=.622*E/(100.*P(I,K)-E)
        IF(QES(I,K).LE.1.E-08)QES(I,K)=1.E-08
        IF(Qx(I,K,N).GT.QES(I,K))Qx(I,K,N)=QES(I,K)
        TV(I,K)=Tx(I,K,N)+.608*Qx(I,K,N)*Tx(I,K,N)
        endif
      enddo
      enddo
!
         do i=its,itf
           if(ierrxx(i).eq.0)then
             Z(I,KTS)=max(0.,Z1(I))-(ALOG(P(I,KTS))- &
                 ALOG(PSUR(I)))*287.*TV(I,KTS)/9.81
           endif
         enddo

! --- calculate heights
         DO K=kts+1,ktf
         do i=its,itf
           if(ierrxx(i).eq.0)then
              TVBAR=.5*TV(I,K)+.5*TV(I,K-1)
              Z(I,K)=Z(I,K-1)-(ALOG(P(I,K))- &
               ALOG(P(I,K-1)))*287.*TVBAR/9.81
           endif
         enddo
         enddo
!
!--- calculate moist static energy - HE
!    saturated moist static energy - HES
!
       DO k=kts,ktf
       do i=its,itf
         if(ierrxx(i).eq.0)then
         HE(I,K)=9.81*Z(I,K)+1004.*Tx(I,K,n)+2.5E06*Qx(I,K,n)
         HES(I,K)=9.81*Z(I,K)+1004.*Tx(I,K,n)+2.5E06*QES(I,K)
         IF(HE(I,K).GE.HES(I,K))HE(I,K)=HES(I,K)
         endif
      enddo
      enddo

! cup levels
!
      do k=kts+1,ktf
      do i=its,itf
        if(ierrxx(i).eq.0)then
        qes_cup(i,k)=.5*(qes(i,k-1)+qes(i,k))
        q_cup(i,k)=.5*(qx(i,k-1,n)+qx(i,k,n))
        hes_cup(i,k)=.5*(hes(i,k-1)+hes(i,k))
        he_cup(i,k)=.5*(he(i,k-1)+he(i,k))
        if(he_cup(i,k).gt.hes_cup(i,k))he_cup(i,k)=hes_cup(i,k)
        z_cup(i,k)=.5*(z(i,k-1)+z(i,k))
        p_cup(i,k)=.5*(p(i,k-1)+p(i,k))
        t_cup(i,k)=.5*(tx(i,k-1,n)+tx(i,k,n))
        gamma_cup(i,k)=(xl/cp)*(xl/(rv*t_cup(i,k) &
                       *t_cup(i,k)))*qes_cup(i,k)
        endif
      enddo
      enddo
      do i=its,itf
        if(ierrxx(i).eq.0)then
        qes_cup(i,1)=qes(i,1)
        q_cup(i,1)=qx(i,1,n)
        hes_cup(i,1)=hes(i,1)
        he_cup(i,1)=he(i,1)
        z_cup(i,1)=.5*(z(i,1)+z1(i))
        p_cup(i,1)=.5*(p(i,1)+psur(i))
        t_cup(i,1)=tx(i,1,n)
        gamma_cup(i,1)=xl/cp*(xl/(rv*t_cup(i,1) &
                       *t_cup(i,1)))*qes_cup(i,1)
        endif
      enddo
!
!
!------- DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22
!
      CALL cup_MAXIMI(HE_CUP,3,KBMAX,K22XX,ierrxx, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
       DO 36 i=its,itf
         IF(ierrxx(I).eq.0.)THEN
         IF(K22xx(I).GE.KBMAX(i))ierrxx(i)=2
         endif
 36   CONTINUE
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!
      call cup_kbcon(ierrc,cap_max_increment,1,k22xx,kbconxx,he_cup,hes_cup, &
           ierrxx,kbmax,p_cup,cap_max, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!--- increase detrainment in stable layers
!
      CALL cup_minimi(HEs_cup,Kbconxx,kstabm,kstabi,ierrxx,  &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      do i=its,itf
      IF(ierrxx(I).eq.0.)THEN
        if(kstabm(i)-1.gt.kstabi(i))then
           do k=kstabi(i),kstabm(i)-1
             cd(i,k)=cd(i,k-1)+1.5*entr_rate
             if(cd(i,k).gt.10.0*entr_rate)cd(i,k)=10.0*entr_rate
           enddo
        ENDIF
      ENDIF
      ENDDO
!
!--- calculate incloud moist static energy
!
      call cup_up_he(k22xx,hkb,z_cup,cd,mentr_rate,he_cup,hc, &
           kbconxx,ierrxx,dby,he,hes_cup,'deep', &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)

!--- DETERMINE CLOUD TOP - KTOP
!
      call cup_ktop(ierrc,1,dby,kbconxx,ktopxx,ierrxx, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!c--- normalized updraft mass flux profile
!
      call cup_up_nms(0,zu,z_cup,mentr_rate,cd,kbconxx,ktopxx,ierrxx,k22xx, &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
!
!--- calculate workfunctions for updrafts
!
      call cup_up_aa0(aa0,z,zu,dby,GAMMA_CUP,t_cup, &
           kbconxx,ktopxx,ierrxx,           &
           itf,jtf,ktf, &
           its,ite, jts,jte, kts,kte)
      do i=its,itf
       if(ierrxx(i).eq.0)axx(i,n)=aa0(i)
      enddo
100   continue
     END SUBROUTINE cup_axx


   SUBROUTINE cup_forcing_ens(closure_n,xland,aa0,aa1,xaa0,mbdt,dtime,ierr,ierr2,ierr3,&
              xf_ens,j,name,maxens,iens,iedt,maxens2,maxens3,mconv,    &
              p_cup,ktop,omeg,zd,k22,zu,pr_ens,edt,kbcon,massflx,      &
              iact_old_gr,dir,ensdim,massfln,icoic,                    &
              ids,ide, jds,jde, kds,kde,                               &
              ims,ime, jms,jme, kms,kme,                               &
              its,ite, jts,jte, kts,kte                               )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ids,ide, jds,jde, kds,kde,           &
        ims,ime, jms,jme, kms,kme,           &
        its,ite, jts,jte, kts,kte
     integer, intent (in   )              ::                           &
        j,ensdim,maxens,iens,iedt,maxens2,maxens3
  !
  ! ierr error value, maybe modified in this routine
  ! pr_ens = precipitation ensemble
  ! xf_ens = mass flux ensembles
  ! massfln = downdraft mass flux ensembles used in next timestep
  ! omeg = omega from large scale model
  ! mconv = moisture convergence from large scale model
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects (ensemble dependent)
  ! edt     = epsilon
  ! dir     = "storm motion"
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! iact_gr_old = flag to tell where convection was active
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! name        = deep or shallow convection flag
  !
     real,    dimension (ims:ime,jms:jme,1:ensdim)                     &
        ,intent (inout)                   ::                           &
        pr_ens
     real,    dimension (ims:ime,jms:jme,1:ensdim)                     &
        ,intent (out  )                   ::                           &
        xf_ens,massfln
     real,    dimension (ims:ime,jms:jme)                              &
        ,intent (in   )                   ::                           &
        massflx
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        omeg,zd,zu,p_cup
     real,    dimension (its:ite,1:maxens)                             &
        ,intent (in   )                   ::                           &
        xaa0
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        aa1,edt,dir,mconv,xland
     real,    dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        aa0,closure_n
     real,    dimension (1:maxens)                                     &
        ,intent (in   )                   ::                           &
        mbdt
     real                                                              &
        ,intent (in   )                   ::                           &
        dtime
     integer, dimension (its:ite,jts:jte)                              &
        ,intent (in   )                   ::                           &
        iact_old_gr
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        k22,kbcon,ktop
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr,ierr2,ierr3
     integer                                                           &
        ,intent (in   )                   ::                           &
        icoic
      character *(*), intent (in)         ::                           &
       name
!
!  local variables in this routine
!

     real,    dimension (1:maxens3)       ::                           &
       xff_ens3
     real,    dimension (1:maxens)        ::                           &
       xk
     integer                              ::                           &
       i,k,nall,n,ne,nens,nens3,iresult,iresultd,iresulte,mkxcrt,kclim
     parameter (mkxcrt=15)
     real                                 ::                           &
       a1,massfld,xff0,xomg,aclim1,aclim2,aclim3,aclim4
     real,    dimension(1:mkxcrt)         ::                           &
       pcrit,acrit,acritt

     integer :: itf,nall2
     integer :: ixxx,irandom
     real :: xxx(3)
     integer,  dimension (2) :: seed

     itf=ite
       seed=0
       irandom=1

      DATA PCRIT/850.,800.,750.,700.,650.,600.,550.,500.,450.,400.,    &
                 350.,300.,250.,200.,150./
      DATA ACRIT/.0633,.0445,.0553,.0664,.075,.1082,.1521,.2216,       &
                 .3151,.3677,.41,.5255,.7663,1.1686,1.6851/
!  GDAS DERIVED ACRIT
      DATA ACRITT/.203,.515,.521,.566,.625,.665,.659,.688,             &
                  .743,.813,.886,.947,1.138,1.377,1.896/
!
       nens=0

!--- LARGE SCALE FORCING
!
       do i=its,itf
        if(ierr(i).eq.0)then
          seed(1)=int(aa0(i))
          seed(2)=int(aa1(i))
          exit
        endif
       enddo
       DO 100 i=its,itf
!       if(i.eq.ipr.and.j.eq.jpr)print *,'ierr = ',ierr(i)
          if(name.eq.'deeps'.and.ierr(i).gt.995)then
!          print *,i,j,ierr(i),aa0(i)
           aa0(i)=0.
           ierr(i)=0
          endif
          IF(ierr(i).eq.0)then
!           kclim=0
           do k=mkxcrt,1,-1
             if(p_cup(i,ktop(i)).lt.pcrit(k))then
               kclim=k
               go to 9
             endif
           enddo
           if(p_cup(i,ktop(i)).ge.pcrit(1))kclim=1
 9         continue
           kclim=max(kclim,1)
           k=max(kclim-1,1)
           aclim1=acrit(kclim)*1.e3
           aclim2=acrit(k)*1.e3
           aclim3=acritt(kclim)*1.e3
           aclim4=acritt(k)*1.e3
!           print *,'p_cup(ktop(i)),kclim,pcrit(kclim)'
!           print *,p_cup(i,ktop(i)),kclim,pcrit(kclim)
!           print *,'aclim1,aclim2,aclim3,aclim4'
!           print *,aclim1,aclim2,aclim3,aclim4
!           print *,dtime,name,ierr(i),aa1(i),aa0(i)
!          print *,dtime,name,ierr(i),aa1(i),aa0(i)
!
!--- treatment different for this closure
!
             if(name.eq.'deeps')then
!
                xff0= (AA1(I)-AA0(I))/DTIME
                xff_ens3(1)=(AA1(I)-AA0(I))/dtime
                xff_ens3(2)=.9*xff_ens3(1)
                xff_ens3(3)=1.1*xff_ens3(1)
                call random_number (xxx)
                xxx(2)=min(1.,xxx(2))
                if(xxx(3).gt.xxx(1))xxx(2)=-xxx(2)
                xff_ens3(13)=xff_ens3(1)*(1.-xxx(2))
!
!--- more original Arakawa-Schubert (climatologic value of aa0)
!
!
!--- omeg is in bar/s, mconv done with omeg in Pa/s
!     more like Brown (1979), or Frank-Cohen (199?)
!
                xff_ens3(4)=-omeg(i,k22(i))/9.81
                xff_ens3(5)=-omeg(i,kbcon(i))/9.81
                xff_ens3(6)=-omeg(i,1)/9.81
                do k=2,kbcon(i)-1
                  xomg=-omeg(i,k)/9.81
                  if(xomg.gt.xff_ens3(6))xff_ens3(6)=xomg
                enddo
                call random_number (xxx)
                xxx(2)=min(1.,xxx(2))
                if(xxx(3).gt.xxx(1))xxx(2)=-xxx(2)
                xff_ens3(14)=xff_ens3(6)*(1.-xxx(2))
!
!--- more like Krishnamurti et al.
!
                xff_ens3(7)=mconv(i)
                xff_ens3(8)=mconv(i)
                xff_ens3(9)=mconv(i)
                call random_number (xxx)
                xxx(2)=min(1.,xxx(2))
                if(xxx(3).gt.xxx(1))xxx(2)=-xxx(2)
                xff_ens3(15)=xff_ens3(8)*(1.-xxx(2))
!
!--- more like Fritsch Chappel or Kain Fritsch (plus triggers)
!
                xff_ens3(10)=AA1(I)/(60.*20.)
                xff_ens3(11)=AA1(I)/(60.*30.)
                xff_ens3(12)=AA1(I)/(60.*40.)
                call random_number (xxx)
                xxx(2)=min(1.,xxx(2))
                if(xxx(3).gt.xxx(1))xxx(2)=-xxx(2)
                xff_ens3(16)=xff_ens3(11)*(1.-xxx(2))
                !write(0,*)'mid_forc ',xff_ens3(1:16)
!
!--- more original Arakawa-Schubert (climatologic value of aa0)
!
!               xff_ens3(13)=max(0.,(AA1(I)-aclim1)/dtime)
!               xff_ens3(14)=max(0.,(AA1(I)-aclim2)/dtime)
!               xff_ens3(15)=max(0.,(AA1(I)-aclim3)/dtime)
!               xff_ens3(16)=max(0.,(AA1(I)-aclim4)/dtime)
!               if(ierr2(i).gt.0.or.ierr3(i).gt.0)then
!                 xff_ens3(10)=0.
!                 xff_ens3(11)=0.
!                 xff_ens3(12)=0.
!                 xff_ens3(13)=0.
!                 xff_ens3(14)=0.
!                 xff_ens3(15)=0.
!                 xff_ens3(16)=0.
!               endif

                do nens=1,maxens
                   XK(nens)=(XAA0(I,nens)-AA1(I))/MBDT(2)
                   if(xk(nens).le.0.and.xk(nens).gt.-1.e-6) &
                           xk(nens)=-1.e-6
                   if(xk(nens).gt.0.and.xk(nens).lt.1.e-6) &
                           xk(nens)=1.e-6
                enddo
!
!--- add up all ensembles
!
                do 350 ne=1,maxens
!
!--- for every xk, we have maxens3 xffs
!--- iens is from outermost ensemble (most expensive!
!
!--- iedt (maxens2 belongs to it)
!--- is from second, next outermost, not so expensive
!
!--- so, for every outermost loop, we have maxens*maxens2*3
!--- ensembles!!! nall would be 0, if everything is on first
!--- loop index, then ne would start counting, then iedt, then iens....
!
                   iresult=0
                   iresultd=0
                   iresulte=0
                   nall=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3 &
                        +(ne-1)*maxens3
!
!--- check for upwind convection
!                  iresult=0
                   massfld=0.

!                  call cup_direction2(i,j,dir,iact_old_gr, &
!                       massflx,iresult,1,                  &
!                       massfld,                            &
!                       ids,ide, jds,jde, kds,kde,          &
!                       ims,ime, jms,jme, kms,kme,          &
!                       its,ite, jts,jte, kts,kte          )
!                  if(i.eq.ipr.and.j.eq.jpr.and.iedt.eq.1.and.ne.eq.1)then
!                  if(iedt.eq.1.and.ne.eq.1)then
!                   print *,massfld,ne,iedt,iens
!                   print *,xk(ne),xff_ens3(1),xff_ens3(2),xff_ens3(3)
!                  endif
!                  print *,i,j,massfld,aa0(i),aa1(i)
                   IF(XK(ne).lt.0.and.xff0.gt.0.)iresultd=1
                   iresulte=max(iresult,iresultd)
                   iresulte=1
                   if(iresulte.eq.1)then
!
!--- special treatment for stability closures
!

                      if(xff0.gt.0.)then
                         xf_ens(i,j,nall+1)=max(0.,-xff_ens3(1)/xk(ne)) &
                                        +massfld
                         xf_ens(i,j,nall+2)=max(0.,-xff_ens3(2)/xk(ne)) &
                                        +massfld
                         xf_ens(i,j,nall+3)=max(0.,-xff_ens3(3)/xk(ne)) &
                                        +massfld
                         xf_ens(i,j,nall+13)=max(0.,-xff_ens3(13)/xk(ne)) &
                                        +massfld
                      else
                         xf_ens(i,j,nall+1)=massfld
                         xf_ens(i,j,nall+2)=massfld
                         xf_ens(i,j,nall+3)=massfld
                         xf_ens(i,j,nall+13)=massfld
                      endif
!
!--- if iresult.eq.1, following independent of xff0
!
                         xf_ens(i,j,nall+4)=max(0.,xff_ens3(4) &
                            +massfld)
                         xf_ens(i,j,nall+5)=max(0.,xff_ens3(5) &
                                        +massfld)
                         xf_ens(i,j,nall+6)=max(0.,xff_ens3(6) &
                                        +massfld)
                         xf_ens(i,j,nall+14)=max(0.,xff_ens3(14) &
                                        +massfld)
                         a1=max(1.e-3,pr_ens(i,j,nall+7))
                         xf_ens(i,j,nall+7)=max(0.,xff_ens3(7) &
                                     /a1)
                         a1=max(1.e-3,pr_ens(i,j,nall+8))
                         xf_ens(i,j,nall+8)=max(0.,xff_ens3(8) &
                                     /a1)
                         a1=max(1.e-3,pr_ens(i,j,nall+9))
                         xf_ens(i,j,nall+9)=max(0.,xff_ens3(9) &
                                     /a1)
                         a1=max(1.e-3,pr_ens(i,j,nall+8))
                         xf_ens(i,j,nall+15)=max(0.,xff_ens3(15) &
                                     /a1)
                         if(XK(ne).lt.0.)then
                            xf_ens(i,j,nall+10)=max(0., &
                                        -xff_ens3(10)/xk(ne)) &
                                        +massfld
                            xf_ens(i,j,nall+11)=max(0., &
                                        -xff_ens3(11)/xk(ne)) &
                                        +massfld
                            xf_ens(i,j,nall+12)=max(0., &
                                        -xff_ens3(12)/xk(ne)) &
                                        +massfld
                            xf_ens(i,j,nall+16)=max(0., &
                                        -xff_ens3(16)/xk(ne)) &
                                        +massfld
                         else
                            xf_ens(i,j,nall+10)=massfld
                            xf_ens(i,j,nall+11)=massfld
                            xf_ens(i,j,nall+12)=massfld
                            xf_ens(i,j,nall+16)=massfld
                         endif
                      if(icoic.ge.1)then
                      closure_n(i)=0.
                      xf_ens(i,j,nall+1)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+2)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+3)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+4)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+5)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+6)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+7)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+8)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+9)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+10)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+11)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+12)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+13)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+14)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+15)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+16)=xf_ens(i,j,nall+icoic)
                      endif
!
! replace 13-16 for now with other stab closures
! (13 gave problems for mass model)
!
!                     xf_ens(i,j,nall+13)=xf_ens(i,j,nall+1)
!                     if(icoic.eq.0)xf_ens(i,j,nall+14)=xf_ens(i,j,nall+13)
!                     xf_ens(i,j,nall+15)=xf_ens(i,j,nall+11)
!                     xf_ens(i,j,nall+16)=xf_ens(i,j,nall+11)
!                     xf_ens(i,j,nall+7)=xf_ens(i,j,nall+4)
!                     xf_ens(i,j,nall+8)=xf_ens(i,j,nall+5)
!                     xf_ens(i,j,nall+9)=xf_ens(i,j,nall+6)
!
!--- store new for next time step
!
                      do nens3=1,maxens3
                        massfln(i,j,nall+nens3)=edt(i) &
                                                *xf_ens(i,j,nall+nens3)
                        massfln(i,j,nall+nens3)=max(0., &
                                              massfln(i,j,nall+nens3))
                      enddo
!
!
!--- do some more on the caps!!! ne=1 for 175, ne=2 for 100,....
!
!     do not care for caps here for closure groups 1 and 5,
!     they are fine, do not turn them off here
!
!
                if(ne.eq.2.and.ierr2(i).gt.0)then
                      xf_ens(i,j,nall+1) =0.
                      xf_ens(i,j,nall+2) =0.
                      xf_ens(i,j,nall+3) =0.
                      xf_ens(i,j,nall+4) =0.
                      xf_ens(i,j,nall+5) =0.
                      xf_ens(i,j,nall+6) =0.
                      xf_ens(i,j,nall+7) =0.
                      xf_ens(i,j,nall+8) =0.
                      xf_ens(i,j,nall+9) =0.
                      xf_ens(i,j,nall+10)=0.
                      xf_ens(i,j,nall+11)=0.
                      xf_ens(i,j,nall+12)=0.
                      xf_ens(i,j,nall+13)=0.
                      xf_ens(i,j,nall+14)=0.
                      xf_ens(i,j,nall+15)=0.
                      xf_ens(i,j,nall+16)=0.
                      massfln(i,j,nall+1)=0.
                      massfln(i,j,nall+2)=0.
                      massfln(i,j,nall+3)=0.
                      massfln(i,j,nall+4)=0.
                      massfln(i,j,nall+5)=0.
                      massfln(i,j,nall+6)=0.
                      massfln(i,j,nall+7)=0.
                      massfln(i,j,nall+8)=0.
                      massfln(i,j,nall+9)=0.
                      massfln(i,j,nall+10)=0.
                      massfln(i,j,nall+11)=0.
                      massfln(i,j,nall+12)=0.
                      massfln(i,j,nall+13)=0.
                      massfln(i,j,nall+14)=0.
                      massfln(i,j,nall+15)=0.
                      massfln(i,j,nall+16)=0.
                endif
                if(ne.eq.3.and.ierr3(i).gt.0)then
                      xf_ens(i,j,nall+1) =0.
                      xf_ens(i,j,nall+2) =0.
                      xf_ens(i,j,nall+3) =0.
                      xf_ens(i,j,nall+4) =0.
                      xf_ens(i,j,nall+5) =0.
                      xf_ens(i,j,nall+6) =0.
                      xf_ens(i,j,nall+7) =0.
                      xf_ens(i,j,nall+8) =0.
                      xf_ens(i,j,nall+9) =0.
                      xf_ens(i,j,nall+10)=0.
                      xf_ens(i,j,nall+11)=0.
                      xf_ens(i,j,nall+12)=0.
                      xf_ens(i,j,nall+13)=0.
                      xf_ens(i,j,nall+14)=0.
                      xf_ens(i,j,nall+15)=0.
                      xf_ens(i,j,nall+16)=0.
                      massfln(i,j,nall+1)=0.
                      massfln(i,j,nall+2)=0.
                      massfln(i,j,nall+3)=0.
                      massfln(i,j,nall+4)=0.
                      massfln(i,j,nall+5)=0.
                      massfln(i,j,nall+6)=0.
                      massfln(i,j,nall+7)=0.
                      massfln(i,j,nall+8)=0.
                      massfln(i,j,nall+9)=0.
                      massfln(i,j,nall+10)=0.
                      massfln(i,j,nall+11)=0.
                      massfln(i,j,nall+12)=0.
                      massfln(i,j,nall+13)=0.
                      massfln(i,j,nall+14)=0.
                      massfln(i,j,nall+15)=0.
                      massfln(i,j,nall+16)=0.
                endif

                   endif
 350            continue
! ne=1, cap=175
!
                   nall=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3
! ne=2, cap=100
!
                   nall2=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3 &
                        +(2-1)*maxens3
                      xf_ens(i,j,nall+4) = xf_ens(i,j,nall2+4)
                      xf_ens(i,j,nall+5) =xf_ens(i,j,nall2+5)
                      xf_ens(i,j,nall+6) =xf_ens(i,j,nall2+6)
                      xf_ens(i,j,nall+14) =xf_ens(i,j,nall2+14)
                      xf_ens(i,j,nall+7) =xf_ens(i,j,nall2+7)
                      xf_ens(i,j,nall+8) =xf_ens(i,j,nall2+8)
                      xf_ens(i,j,nall+9) =xf_ens(i,j,nall2+9)
                      xf_ens(i,j,nall+15) =xf_ens(i,j,nall2+15)
                      xf_ens(i,j,nall+10)=xf_ens(i,j,nall2+10)
                      xf_ens(i,j,nall+11)=xf_ens(i,j,nall2+11)
                      xf_ens(i,j,nall+12)=xf_ens(i,j,nall2+12)
                      xf_ens(i,j,nall+16)=xf_ens(i,j,nall2+16)
                go to 100
             endif
          elseif(ierr(i).ne.20.and.ierr(i).ne.0)then
             do n=1,ensdim
               xf_ens(i,j,n)=0.
               massfln(i,j,n)=0.
             enddo
          endif
 100   continue
      nall=48
      !write(0,*)'end of forcing',nall,xf_ens(1,1,nall+1:nall+16)


   END SUBROUTINE cup_forcing_ens



!====================================================================
   SUBROUTINE gdinit(RTHCUTEN,RQVCUTEN,RQCCUTEN,RQICUTEN,           &
                        MASS_FLUX,cp,restart,                       &
                        P_QC,P_QI,P_FIRST_SCALAR,                   &
                        RTHFTEN, RQVFTEN,                           &
                        APR_GR,APR_W,APR_MC,APR_ST,APR_AS,          &
                        APR_CAPMA,APR_CAPME,APR_CAPMI,              &
                        allowed_to_read,                            &
                        ids, ide, jds, jde, kds, kde,               &
                        ims, ime, jms, jme, kms, kme,               &
                        its, ite, jts, jte, kts, kte               )
!--------------------------------------------------------------------
   IMPLICIT NONE
!--------------------------------------------------------------------
   LOGICAL , INTENT(IN)           ::  restart,allowed_to_read
   INTEGER , INTENT(IN)           ::  ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      its, ite, jts, jte, kts, kte
   INTEGER , INTENT(IN)           ::  P_FIRST_SCALAR, P_QI, P_QC
   REAL,     INTENT(IN)           ::  cp

   REAL,     DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) ::       &
                                                          RTHCUTEN, &
                                                          RQVCUTEN, &
                                                          RQCCUTEN, &
                                                          RQICUTEN

   REAL,     DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) ::       &
                                                          RTHFTEN,  &
                                                          RQVFTEN

   REAL,     DIMENSION( ims:ime , jms:jme ) , INTENT(OUT) ::        &
                                APR_GR,APR_W,APR_MC,APR_ST,APR_AS,  &
                                APR_CAPMA,APR_CAPME,APR_CAPMI,      &
                                MASS_FLUX

   IF(.not.restart)THEN
        RTHCUTEN=0.
        RQVCUTEN=0.
        RTHFTEN=0.
        RQVFTEN=0.

     IF (P_QC .ge. P_FIRST_SCALAR) THEN
           RQCCUTEN=0.
     ENDIF

     IF (P_QI .ge. P_FIRST_SCALAR) THEN
           RQICUTEN=0.
     ENDIF

        mass_flux=0.

   ENDIF
        APR_GR=0.
        APR_ST=0.
        APR_W=0.
        APR_MC=0.
        APR_AS=0.
        APR_CAPMA=0.
        APR_CAPME=0.
        APR_CAPMI=0.

   END SUBROUTINE gdinit



   SUBROUTINE neg_check(dt,q,outq,outt,outqc,pret,its,ite,kts,kte,itf,ktf)

   INTEGER,      INTENT(IN   ) ::            its,ite,kts,kte,itf,ktf

     real, dimension (its:ite,kts:kte+1  )                    ,                 &
      intent(inout   ) ::                                                     &
       outq,outt,outqc
     real, dimension (its:ite,kts:kte  )                    ,                 &
      intent(inout   ) ::                                                     &
       q
     real, dimension (its:ite  )                            ,                 &
      intent(inout   ) ::                                                     &
       pret
     real                                                                     &
        ,intent (in  )                   ::                                   &
        dt
     real :: thresh,qmem,qmemf,qmem2,qtest,qmem1
!
! first do check on vertical heating rate
!
      thresh=200.01
      do i=its,itf
      qmemf=1.
      qmem=0.
      do k=kts,ktf
         qmem=outt(i,k)*86400.
         if(qmem.gt.2.*thresh)then
           qmem2=2.*thresh/qmem
           qmemf=min(qmemf,qmem2)
!
!
           print *,'1',' adjusted massflux by factor ',i,k,qmem,qmem2,qmemf
         endif
         if(qmem.lt.-thresh)then
           qmem2=-thresh/qmem
           qmemf=min(qmemf,qmem2)
!
!
           print *,'2',' adjusted massflux by factor ',i,k,qmem,qmem2,qmemf
         endif
      enddo
      do k=kts,ktf
         outq(i,k)=outq(i,k)*qmemf
         outt(i,k)=outt(i,k)*qmemf
         outqc(i,k)=outqc(i,k)*qmemf
      enddo
      pret(i)=pret(i)*qmemf
      enddo
!
! check whether routine produces negative q's. This can happen, since
! tendencies are calculated based on forced q's. This should have no
! influence on conservation properties, it scales linear through all
! tendencies
!
      thresh=1.e-10
      do i=its,itf
      qmemf=1.
      do k=kts,ktf
         qmem=outq(i,k)
         if(abs(qmem).gt.0.)then
         qtest=q(i,k)+outq(i,k)*dt
         if(qtest.lt.thresh)then
!
! qmem2 would be the maximum allowable tendency
!
           qmem1=outq(i,k)
           qmem2=(thresh-q(i,k))/dt
           qmemf=min(qmemf,qmem2/qmem1)
           qmemf=max(0.,qmemf)
           print *,'4 adjusted tendencies ',i,k,qmem,qmem2,qmemf
         endif
         endif
      enddo
      do k=kts,ktf
         outq(i,k)=outq(i,k)*qmemf
         outt(i,k)=outt(i,k)*qmemf
         outqc(i,k)=outqc(i,k)*qmemf
      enddo
      pret(i)=pret(i)*qmemf
      enddo

   END SUBROUTINE neg_check


!-------------------------------------------------------
   SUBROUTINE cup3_up_aa0(aa0,z,zu,dby,GAMMA_CUP,t_cup,       &
              kbcon,ktop,ierr,                               &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,                                     &
        its,ite, jts,jte, kts,kte
  ! aa0 cloud work function
  ! gamma_cup = gamma on model cloud levels
  ! t_cup = temperature (Kelvin) on model cloud levels
  ! dby = buoancy term
  ! zu= normalized updraft mass flux
  ! z = heights of model levels
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        z,zu,gamma_cup,t_cup,dby
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop
!
! input and output
!


     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        aa0
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
     real                                 ::                           &
        dz,da
!
        do i=its,itf
         aa0(i)=0.
        enddo
        DO 100 k=kts+1,ktf
        DO 100 i=its,itf
         IF(ierr(i).ne.0)GO TO 100
         IF(K.LE.KBCON(I))GO TO 100
         IF(K.Gt.KTOP(I))GO TO 100
         DZ=Z(I,K)-Z(I,K-1)
         da=zu(i,k)*DZ*(9.81/(1004.*( &
                (T_cup(I,K)))))*DBY(I,K-1)/ &
             (1.+GAMMA_CUP(I,K))
         IF(K.eq.KTOP(I).and.da.le.0.)go to 100
         AA0(I)=AA0(I)+da
         if(aa0(i).lt.0.)aa0(i)=0.
100     continue

   END SUBROUTINE cup3_up_aa0


   SUBROUTINE cup3_up_he(k22,hkb,z_cup,cd,entr,he_cup,hc,     &
              kbcon,ierr,dby,he,hes_cup,name,                     &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  itf,jtf,ktf,           &
                                  its,ite, jts,jte, kts,kte
      character *(*), intent (in)        ::                           &
       name

  ! hc = cloud moist static energy
  ! hkb = moist static energy at originating level
  ! he = moist static energy on model levels
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function
  ! z_cup = heights of model cloud levels
  ! entr = entrainment rate
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        he,he_cup,hes_cup,z_cup,cd
  ! entr= entrainment rate
     real                                                              &
        ,intent (in   )                   ::                           &
        entr
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,k22
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        hc,dby
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        hkb
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
     real                                 ::                           &
        dz
!
!--- moist static energy inside cloud
!
      do k=kts,ktf
      do i=its,itf
       hc(i,k)=0.
       DBY(I,K)=0.
      enddo
      enddo
      do i=its,itf
       hkb(i)=0.
      enddo
      do i=its,itf
      if(ierr(I).eq.0)then
      hkb(i)=he_cup(i,k22(i))
          if(name.eq.'shallow')then
             do k=1,k22(i)
               hkb(i)=max(hkb(i),he_cup(i,k))
             enddo
          endif
      do k=1,k22(i)
        hc(i,k)=he_cup(i,k)
!       DBY(I,K)=0.
      enddo
      do k=k22(i),kbcon(i)-1
        hc(i,k)=hkb(i)
!       DBY(I,K)=0.
      enddo
        k=kbcon(i)
        hc(i,k)=hkb(i)
        DBY(I,Kbcon(i))=Hkb(I)-HES_cup(I,K)
      endif
      enddo
      do k=kts+1,ktf
      do i=its,itf
        if(k.gt.kbcon(i).and.ierr(I).eq.0)then
           DZ=Z_cup(i,K)-Z_cup(i,K-1)
           HC(i,K)=(HC(i,K-1)*(1.-.5*CD(i,K)*DZ)+entr* &
                DZ*HE(i,K-1))/(1.+entr*DZ-.5*cd(i,k)*dz)
           DBY(I,K)=HC(I,K)-HES_cup(I,K)
        endif
      enddo

      enddo

   END SUBROUTINE cup3_up_he


   SUBROUTINE cup3_up_moisture(name,ierr,z_cup,qc,qrc,pw,pwav,     &
              kbcon,ktop,cd,dby,mentr_rate,clw_all,                  &
              q,GAMMA_cup,zu,qes_cup,k22,qe_cup,xl,          &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  itf,jtf,ktf,           &
                                  its,ite, jts,jte, kts,kte
      character *(*), intent (in)        ::                           &
       name

  ! cd= detrainment function
  ! q = environmental q on model levels
  ! qe_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function
  ! zu = normalized updraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! mentr_rate = entrainment rate
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        q,zu,gamma_cup,qe_cup,dby,qes_cup,z_cup,cd
  ! entr= entrainment rate
     real                                                              &
        ,intent (in   )                   ::                           &
        mentr_rate,xl
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop,k22
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
   ! qc = cloud q (including liquid water) after entrainment
   ! qrch = saturation q in cloud
   ! qrc = liquid water content in cloud after rainout
   ! pw = condensate that will fall out at that level
   ! pwav = totan normalized integrated condensate (I1)
   ! c0 = conversion rate (cloud to rain)

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qc,qrc,pw,clw_all
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pwav
!
!  local variables in this routine
!

     integer                              ::                           &
        iall,i,k
     real                                 ::                           &
        dh,qrch,c0,dz,radius
!
        iall=0
        c0=.002
!
!--- no precip for small clouds
!
        if(name.eq.'shallow')c0=0.

!       if(mentr_rate.gt.0.)then
!         radius=.2/mentr_rate
!         if(radius.lt.900.)c0=0.
!         if(radius.lt.900.)iall=0
!       endif
        do i=its,itf
          pwav(i)=0.
        enddo
        do k=kts,ktf
        do i=its,itf
          pw(i,k)=0.
          qc(i,k)=0.
          if(ierr(i).eq.0)qc(i,k)=qes_cup(i,k)
          clw_all(i,k)=0.
          qrc(i,k)=0.
        enddo
        enddo
      do i=its,itf
      if(ierr(I).eq.0)then
      do k=k22(i),kbcon(i)-1
        qc(i,k)=qe_cup(i,k22(i))
      enddo
      endif
      enddo

        DO 100 k=kts+1,ktf
        DO 100 i=its,itf
         IF(ierr(i).ne.0)GO TO 100
         IF(K.Lt.KBCON(I))GO TO 100
         IF(K.Gt.KTOP(I))GO TO 100
         DZ=Z_cup(i,K)-Z_cup(i,K-1)
!
!------    1. steady state plume equation, for what could
!------       be in cloud without condensation
!
!
        QC(i,K)=(QC(i,K-1)*(1.-.5*CD(i,K)*DZ)+mentr_rate* &
                DZ*Q(i,K-1))/(1.+mentr_rate*DZ-.5*cd(i,k)*dz)
!
!--- saturation  in cloud, this is what is allowed to be in it
!
         QRCH=QES_cup(I,K)+(1./XL)*(GAMMA_cup(i,k) &
              /(1.+GAMMA_cup(i,k)))*DBY(I,K)
!
!------- LIQUID WATER CONTENT IN cloud after rainout
!
        clw_all(i,k)=QC(I,K)-QRCH
        QRC(I,K)=(QC(I,K)-QRCH)/(1.+C0*DZ)
        if(qrc(i,k).lt.0.)then
          qrc(i,k)=0.
        endif
!
!-------   3.Condensation
!
         PW(i,k)=c0*dz*QRC(I,K)*zu(i,k)
        if(iall.eq.1)then
          qrc(i,k)=0.
          pw(i,k)=(QC(I,K)-QRCH)*zu(i,k)
          if(pw(i,k).lt.0.)pw(i,k)=0.
        endif
!
!----- set next level
!
         QC(I,K)=QRC(I,K)+qrch
!
!--- integrated normalized ondensate
!
         PWAV(I)=PWAV(I)+PW(I,K)
 100     CONTINUE

   END SUBROUTINE cup3_up_moisture


   SUBROUTINE cup3_up_nms(zu,z_cup,entr,cd,kbcon,ktop,ierr,k22,  &
              itf,jtf,ktf,                        &
              its,ite, jts,jte, kts,kte                        )

   IMPLICIT NONE

!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         itf,jtf,ktf,                                    &
         its,ite, jts,jte, kts,kte
  ! cd= detrainment function
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
         z_cup,cd
  ! entr= entrainment rate
     real                                                              &
        ,intent (in   )                   ::                           &
         entr
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
         kbcon,ktop,k22
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
         ierr
   ! zu is the normalized mass flux

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
         zu
!
!  local variables in this routine
!

     integer                              ::                           &
         i,k
     real                                 ::                           &
         dz
!
!   initialize for this go around
!
       do k=kts,ktf
       do i=its,itf
         zu(i,k)=0.
       enddo
       enddo
!
! do normalized mass budget
!
       do i=its,itf
          IF(ierr(I).eq.0)then
             do k=k22(i),kbcon(i)
               zu(i,k)=1.
             enddo
             DO K=KBcon(i)+1,KTOP(i)
               DZ=Z_cup(i,K)-Z_cup(i,K-1)
               ZU(i,K)=ZU(i,K-1)*(1.+(entr-cd(i,k))*DZ)
             enddo
          endif
       enddo

   END SUBROUTINE cup3_up_nms

   SUBROUTINE cup3_kbcon(cap_inc,iloop,k22,kbcon,he_cup,hes_cup, &
              ierr,kbmax,p_cup,cap_max,                         &
              itf,jtf,ktf,                        &
              its,ite, jts,jte, kts,kte                        )

   IMPLICIT NONE
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
  !
  !
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        he_cup,hes_cup,p_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        cap_max,cap_inc
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbmax
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        kbcon,k22,ierr
     integer                                                           &
        ,intent (in   )                   ::                           &
        iloop
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
     real                                 ::                           &
        pbcdif,plus,hetest
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!
       DO 27 i=its,itf
      kbcon(i)=1
      IF(ierr(I).ne.0)GO TO 27
      KBCON(I)=K22(I)
      GO TO 32
 31   CONTINUE
      KBCON(I)=KBCON(I)+1
 32   CONTINUE
      IF(KBCON(I).GT.KBMAX(i)+2)THEN
         if(iloop.ne.4)ierr(i)=3
!        if(iloop.lt.4)ierr(i)=997
        GO TO 27
      ENDIF
! 32   CONTINUE  changed July2012
      hetest=HE_cup(I,K22(I))
      if(iloop.eq.5)then
       do k=1,k22(i)
         hetest=max(hetest,he_cup(i,k))
       enddo
      endif
      IF(HETEST.LT.HES_cup(I,KBCON(I)))GO TO 31

!     cloud base pressure and max moist static energy pressure
!     i.e., the depth (in mb) of the layer of negative buoyancy
      if(KBCON(I)-K22(I).eq.0)go to 27
      PBCDIF=-P_cup(I,KBCON(I))+P_cup(I,K22(I))
      plus=max(25.,cap_max(i)-float(iloop-1)*cap_inc(i))
      if(iloop.eq.4)plus=cap_max(i)
!
! for shallow convection, if cap_max is greater than 25, it is the pressure at
! pbltop
      if(iloop.eq.5)plus=25.
      if(iloop.eq.5.and.cap_max(i).gt.25)pbcdif=-P_cup(I,KBCON(I))+cap_max(i)
      IF(PBCDIF.GT.plus)THEN
        K22(I)=K22(I)+1
        KBCON(I)=K22(I)+1 ! changed to +1, July2012
        GO TO 32
      ENDIF
 27   CONTINUE

   END SUBROUTINE cup3_kbcon
   SUBROUTINE cup3_ktop(ilo,dby,kbcon,ktop,ierr,              &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
  ! dby = buoancy term
  ! ktop = cloud top (output)
  ! ilo  = flag
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
        dby
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon
     integer                                                           &
        ,intent (in   )                   ::                           &
        ilo
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        ktop
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
!
        DO 42 i=its,itf
        ktop(i)=1
         IF(ierr(I).EQ.0)then
          DO 40 K=KBCON(I)+1,ktf-1
            IF(DBY(I,K).LE.0.)THEN
                KTOP(I)=K-1
                GO TO 41
             ENDIF
  40      CONTINUE
          if(ilo.eq.1)ierr(i)=5
!         if(ilo.eq.2)ierr(i)=998
          GO TO 42
  41     CONTINUE
         do k=ktop(i)+1,ktf
           dby(i,k)=0.
         enddo
         if(kbcon(i).eq.ktop(i))then
            ierr(i)=55
         endif
         endif
  42     CONTINUE

   END SUBROUTINE cup3_ktop


   SUBROUTINE cup3_MAXIMI(ARRAY,KS,KE,MAXX,ierr,              &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         itf,jtf,ktf,                                    &
         its,ite, jts,jte, kts,kte
  ! array input array
  ! x output array with return values
  ! kt output array of levels
  ! ks,kend  check-range
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
         array
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
         ierr,ke
     integer                                                           &
        ,intent (in   )                   ::                           &
         ks
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
         maxx
     real,    dimension (its:ite)         ::                           &
         x
     real                                 ::                           &
         xar
     integer                              ::                           &
         i,k

       DO 200 i=its,itf
       MAXX(I)=KS
       if(ierr(i).eq.0)then
      X(I)=ARRAY(I,KS)
!
       DO 100 K=KS,KE(i)
         XAR=ARRAY(I,K)
         IF(XAR.GE.X(I)) THEN
            X(I)=XAR
            MAXX(I)=K
         ENDIF
 100  CONTINUE
      endif
 200  CONTINUE

   END SUBROUTINE cup3_MAXIMI

   SUBROUTINE cup3_env(z,qes,he,hes,t,q,p,z1,                 &
              psur,ierr,tcrit,itest,xl,cp,                   &
              itf,jtf,ktf,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
  !
  ! ierr error value, maybe modified in this routine
  ! q           = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! tv          = environmental virtual temp
  ! p           = environmental pressure
  ! z           = environmental heights
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  ! psur        = surface pressure
  ! z1          = terrain elevation
  !
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        p,t,q
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        he,hes,qes
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
        z
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        psur,z1
     real                                                              &
        ,intent (in   )                   ::                           &
        xl,cp
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     integer                                                           &
        ,intent (in   )                   ::                           &
        itest
!
!  local variables in this routine
!

     integer                              ::                           &
       i,k,iph
      real, dimension (1:2) :: AE,BE,HT
      real, dimension (its:ite,kts:kte) :: tv
      real :: tcrit,e,tvbar


      HT(1)=XL/CP
      HT(2)=2.834E6/CP
      BE(1)=.622*HT(1)/.286
      AE(1)=BE(1)/273.+ALOG(610.71)
      BE(2)=.622*HT(2)/.286
      AE(2)=BE(2)/273.+ALOG(610.71)
!      print *, 'TCRIT = ', tcrit,its,ite
      DO k=kts,ktf
      do i=its,itf
        if(ierr(i).eq.0)then
!Csgb - IPH is for phase, dependent on TCRIT (water or ice)
        IPH=1
        IF(T(I,K).LE.TCRIT)IPH=2
!       print *, 'AE(IPH),BE(IPH) = ',AE(IPH),BE(IPH),AE(IPH)-BE(IPH),T(i,k),i,k
        E=EXP(AE(IPH)-BE(IPH)/T(I,K))
!       print *, 'P, E = ', P(I,K), E
        QES(I,K)=.622*E/(100.*P(I,K)-E)
        IF(QES(I,K).LE.1.E-08)QES(I,K)=1.E-08
        IF(QES(I,K).LT.Q(I,K))QES(I,K)=Q(I,K)
!       IF(Q(I,K).GT.QES(I,K))Q(I,K)=QES(I,K)
        TV(I,K)=T(I,K)+.608*Q(I,K)*T(I,K)
        endif
      enddo
      enddo
!
!--- z's are calculated with changed h's and q's and t's
!--- if itest=2
!
      if(itest.ne.2)then
         do i=its,itf
           if(ierr(i).eq.0)then
             Z(I,1)=max(0.,Z1(I))-(ALOG(P(I,1))- &
                 ALOG(PSUR(I)))*287.*TV(I,1)/9.81
           endif
         enddo

! --- calculate heights
         DO K=kts+1,ktf
         do i=its,itf
           if(ierr(i).eq.0)then
              TVBAR=.5*TV(I,K)+.5*TV(I,K-1)
              Z(I,K)=Z(I,K-1)-(ALOG(P(I,K))- &
               ALOG(P(I,K-1)))*287.*TVBAR/9.81
           endif
         enddo
         enddo
      else
         do k=kts,ktf
         do i=its,itf
           if(ierr(i).eq.0)then
             z(i,k)=(he(i,k)-1004.*t(i,k)-2.5e6*q(i,k))/9.81
             z(i,k)=max(1.e-3,z(i,k))
           endif
         enddo
         enddo
      endif
!
!--- calculate moist static energy - HE
!    saturated moist static energy - HES
!
       DO k=kts,ktf
       do i=its,itf
         if(ierr(i).eq.0)then
         if(itest.eq.0)HE(I,K)=9.81*Z(I,K)+1004.*T(I,K)+2.5E06*Q(I,K)
         HES(I,K)=9.81*Z(I,K)+1004.*T(I,K)+2.5E06*QES(I,K)
         IF(HE(I,K).GE.HES(I,K))HE(I,K)=HES(I,K)
         endif
      enddo
      enddo

   END SUBROUTINE cup3_env


   SUBROUTINE cup3_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,   &
              he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
              ierr,z1,xl,rv,cp,                                &
              itf,jtf,ktf,                       &
              its,ite, jts,jte, kts,kte                       )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,           &
        its,ite, jts,jte, kts,kte
  !
  ! ierr error value, maybe modified in this routine
  ! q           = environmental mixing ratio
  ! q_cup       = environmental mixing ratio on cloud levels
  ! qes         = environmental saturation mixing ratio
  ! qes_cup     = environmental saturation mixing ratio on cloud levels
  ! t           = environmental temp
  ! t_cup       = environmental temp on cloud levels
  ! p           = environmental pressure
  ! p_cup       = environmental pressure on cloud levels
  ! z           = environmental heights
  ! z_cup       = environmental heights on cloud levels
  ! he          = environmental moist static energy
  ! he_cup      = environmental moist static energy on cloud levels
  ! hes         = environmental saturation moist static energy
  ! hes_cup     = environmental saturation moist static energy on cloud levels
  ! gamma_cup   = gamma on cloud levels
  ! psur        = surface pressure
  ! z1          = terrain elevation
  !
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        qes,q,he,hes,z,p,t
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        psur,z1
     real                                                              &
        ,intent (in   )                   ::                           &
        xl,rv,cp
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!

     integer                              ::                           &
       i,k


      do k=kts+1,ktf
      do i=its,itf
        if(ierr(i).eq.0)then
        qes_cup(i,k)=.5*(qes(i,k-1)+qes(i,k))
        q_cup(i,k)=.5*(q(i,k-1)+q(i,k))
        hes_cup(i,k)=.5*(hes(i,k-1)+hes(i,k))
        he_cup(i,k)=.5*(he(i,k-1)+he(i,k))
        if(he_cup(i,k).gt.hes_cup(i,k))he_cup(i,k)=hes_cup(i,k)
        z_cup(i,k)=.5*(z(i,k-1)+z(i,k))
        p_cup(i,k)=.5*(p(i,k-1)+p(i,k))
        t_cup(i,k)=.5*(t(i,k-1)+t(i,k))
        gamma_cup(i,k)=(xl/cp)*(xl/(rv*t_cup(i,k) &
                       *t_cup(i,k)))*qes_cup(i,k)
        endif
      enddo
      enddo
      do i=its,itf
        if(ierr(i).eq.0)then
        qes_cup(i,1)=qes(i,1)
        q_cup(i,1)=q(i,1)
        hes_cup(i,1)=hes(i,1)
        he_cup(i,1)=he(i,1)
        z_cup(i,1)=.5*(z(i,1)+z1(i))
        p_cup(i,1)=.5*(p(i,1)+psur(i))
        t_cup(i,1)=t(i,1)
        gamma_cup(i,1)=xl/cp*(xl/(rv*t_cup(i,1) &
                       *t_cup(i,1)))*qes_cup(i,1)
        endif
      enddo

   END SUBROUTINE cup3_env_clev

   SUBROUTINE cup_output_ens_3d(xf_ens,ierr,dellat,dellaq,dellaqc,  &
              subt_ens,subq_ens,subt,subq,outtem,outq,outqc,     &
              zu,sub_mas,pre,pw,xmb,ktop,                 &
              j,name,nx,nx2,iens,ierr2,ierr3,pr_ens,             &
              maxens3,ensdim,massfln,                            &
              sig,APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                 &
              APR_CAPMA,APR_CAPME,APR_CAPMI,closure_n,xland1,    &
              weight_GR,weight_W,weight_MC,weight_ST,weight_AS,training,  &
	      itf,jtf,ktf, &
              its,ite, jts,jte, kts,kte)

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,jtf,ktf,     &
        its,ite, jts,jte, kts,kte
     integer, intent (in   )              ::                           &
        j,ensdim,nx,nx2,iens,maxens3,training
  ! xf_ens = ensemble mass fluxes
  ! pr_ens = precipitation ensembles
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble
  ! outtem = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
  ! xmb    = total base mass flux
  ! xfac1  = correction factor
  ! pw = pw -epsilon*pd (ensemble dependent)
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,jts:jte,1:ensdim)                     &
        ,intent (inout)                   ::                           &
       xf_ens,pr_ens,massfln
!srf ------
!    real,    dimension (its:ite,jts:jte)                              &
     real,    dimension (its:ite,jts:jte)                              &
         ,intent (inout)                   ::                          &
               APR_GR,APR_W,APR_MC,APR_ST,APR_AS,APR_CAPMA,            &
               APR_CAPME,APR_CAPMI
     real, dimension( its:ite , jts:jte )                      &
         ,intent(in) :: weight_gr,weight_w,weight_mc,weight_st,weight_as
!-srf---
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        outtem,outq,outqc,subt,subq,sub_mas
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in  )                   ::                           &
        zu
    real,   dimension (its:ite)                                      &
         ,intent (in  )                   ::                           &
        sig
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pre,xmb
     real,    dimension (its:ite)                                      &
        ,intent (inout  )                   ::                           &
        closure_n,xland1
     real,    dimension (its:ite,kts:kte,1:nx)                     &
        ,intent (in   )                   ::                           &
       subt_ens,subq_ens,dellat,dellaqc,dellaq,pw
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        ktop
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr,ierr2,ierr3
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,n,ncount
     real                                 ::                           &
        outtes,ddtes,dtt,dtq,dtqc,dtpw,tuning,prerate,clos_wei,xmbhelp
     real                                 ::                           &
        dtts,dtqs
     real,    dimension (its:ite)         ::                           &
       xfac1,xfac2
     real,    dimension (its:ite)::                           &
       xmb_ske,xmb_ave,xmb_std,xmb_cur,xmbweight
     real,    dimension (its:ite)::                           &
       pr_ske,pr_ave,pr_std,pr_cur
     real,    dimension (its:ite,jts:jte)::                           &
               pr_gr,pr_w,pr_mc,pr_st,pr_as,pr_capma,     &
               pr_capme,pr_capmi
     real, dimension (5) :: weight,wm,wm1,wm2,wm3
     real, dimension (its:ite,5) :: xmb_w

!
      character *(*), intent (in)        ::                           &
       name

!
     weight(1) = -999.  !this will turn off weights
     wm(1)=-999.

     tuning=0.
!
!
      DO k=kts,ktf
      do i=its,itf
        outtem(i,k)=0.
        outq(i,k)=0.
        outqc(i,k)=0.
        subt(i,k)=0.
        subq(i,k)=0.
        sub_mas(i,k)=0.
      enddo
      enddo
      do i=its,itf
        pre(i)=0.
        xmb(i)=0.
         xfac1(i)=0.
         xfac2(i)=0.
        xmbweight(i)=1.
      enddo
      do i=its,itf
        IF(ierr(i).eq.0)then
        do n=(iens-1)*nx*nx2*maxens3+1,iens*nx*nx2*maxens3
           if(pr_ens(i,j,n).le.0.)then
             xf_ens(i,j,n)=0.
           endif
        enddo
        endif
      enddo
!
!--- calculate ensemble average mass fluxes
!

       call massflx_stats(xf_ens,ensdim,nx2,nx,maxens3,  &
            xmb_ave,xmb_std,xmb_cur,xmb_ske,j,ierr,1,    &
            APR_GR,APR_W,APR_MC,APR_ST,APR_AS,           &
            APR_CAPMA,APR_CAPME,APR_CAPMI,               &
            pr_gr,pr_w,pr_mc,pr_st,pr_as,                &
            pr_capma,pr_capme,pr_capmi,                  &
            itf,jtf,ktf,                                 &
            its,ite, jts,jte, kts,kte, &
	    weight_GR,weight_W,weight_MC,weight_ST,weight_AS, training)
       xmb_w=0.
       call massflx_stats(pr_ens,ensdim,nx2,nx,maxens3,  &
            pr_ave,pr_std,pr_cur,pr_ske,j,ierr,2,        &
            APR_GR,APR_W,APR_MC,APR_ST,APR_AS,           &
            APR_CAPMA,APR_CAPME,APR_CAPMI,               &
            pr_gr,pr_w,pr_mc,pr_st,pr_as,                &
            pr_capma,pr_capme,pr_capmi,                  &
            itf,jtf,ktf,                                 &
            its,ite, jts,jte, kts,kte,  &
	    weight_GR,weight_W,weight_MC,weight_ST,weight_AS,training)

!
!-- now do feedback
!
      ddtes=100.
      do i=its,itf
        if(ierr(i).eq.0)then
         if(xmb_ave(i).le.0.)then
              ierr(i)=13
              xmb_ave(i)=0.
         endif
!--srf
         xmb(i)=xmb_ave(i)*sig(i)
!        xmb(i)=max(.1*xmb_ave(i),xmb_ave(i)-tuning*xmb_std(i))
!--srf
! --- Now use proper count of how many closures were actually
!       used in cup_forcing_ens (including screening of some
!       closures over water) to properly normalize xmb
           clos_wei=16./max(1.,closure_n(i))
!!if(i==19 .and. j==27)  write(0,*)'out1',xmb_ave(i),xmb_std(i),closure_n(i)

!srf -tmp          if (xland1(i).lt.0.5)xmb(i)=xmb(i)*clos_wei
!

           if(xmb(i).eq.0.)then
              ierr(i)=19
           endif
           if(xmb(i).gt.100.)then
              ierr(i)=19
           endif
           xfac1(i)=xmb(i)
           xfac2(i)=xmb(i)

        endif
!       if(weight(1).lt.-100.)xfac1(i)=xmb_ave(i)
!       if(weight(1).lt.-100.)xfac2(i)=xmb_ave(i)
      ENDDO
      DO k=kts,ktf
      do i=its,itf
            dtt =0.
            dtts=0.
            dtq =0.
            dtqs=0.
            dtqc=0.
            dtpw=0.
        IF(ierr(i).eq.0.and.k.le.ktop(i))then
           do n=1,nx
              dtt =dtt  + dellat  (i,k,n)
              dtts=dtts + subt_ens(i,k,n)
              dtq =dtq  + dellaq  (i,k,n)
              dtqs=dtqs + subq_ens(i,k,n)
              dtqc=dtqc + dellaqc (i,k,n)
              dtpw=dtpw + pw      (i,k,n)
           enddo
           OUTTEM(I,K)= XMB(I)* dtt /float(nx)
           SUBT  (I,K)= XMB(I)* dtts/float(nx)
           OUTQ  (I,K)= XMB(I)* dtq /float(nx)
           SUBQ  (I,K)= XMB(I)* dtqs/float(nx)
           OUTQC (I,K)= XMB(I)* dtqc/float(nx)

	   PRE(I)=PRE(I)+XMB(I)*dtpw/float(nx)

	   !if(training == 1)
	   !if(k==kts) PRE(I)=pr_ave(i)
	   !if(k==kts) PRE(I)=(apr_w(i,j)+apr_st(i,j)+apr_gr(i,j)+apr_mc(i,j)+apr_as(i,j))*0.2

           sub_mas(i,k)=zu(i,k)*xmb(i)
	   !if(i==19 .and. j==27)  write(0,*)'out2', xmb(i),dtpw,OUTQ(I,K)+SUBQ(I,K),pre(i)
        endif
       enddo
      enddo

      do i=its,itf
        if(ierr(i).eq.0)then
        do k=(iens-1)*nx*nx2*maxens3+1,iens*nx*nx2*maxens3
          massfln(i,j,k)=massfln(i,j,k)*xfac1(i)
          xf_ens(i,j,k)=xf_ens(i,j,k)*xfac1(i)
        enddo
        endif
      ENDDO

!srf-fix for preci
      do i=its,itf
        if(ierr(i).ne. 0)then
            apr_w (i,j)=0.0
	    apr_st(i,j)=0.0
	    apr_gr(i,j)=0.0
	    apr_mc(i,j)=0.0
	    apr_as(i,j)=0.0
        endif
      ENDDO
!srf
   END SUBROUTINE cup_output_ens_3d

   SUBROUTINE massflx_stats(xf_ens,ensdim,maxens,maxens2,maxens3, &
              xt_ave,xt_std,xt_cur,xt_ske,j,ierr,itest,           &
              APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                  &
              APR_CAPMA,APR_CAPME,APR_CAPMI,                      &
              pr_gr,pr_w,pr_mc,pr_st,pr_as,                       &
              pr_capma,pr_capme,pr_capmi,                         &
              itf,jtf,ktf,                                        &
              its,ite, jts,jte, kts,kte,         &
              weight_GR,weight_W,weight_MC,weight_ST,weight_AS, training)
   IMPLICIT NONE

    integer, intent (in   )              ::                                    &
                     j,ensdim,maxens3,maxens,maxens2,itest
     INTEGER,  INTENT(IN   ) :: 					    &
                	      itf,jtf,ktf,		    &
                	      its,ite, jts,jte, kts,kte

     INTEGER,  INTENT(IN   ) ::  training

     real, dimension (its:ite)                                                &
         , intent(inout) ::                                                   &
           xt_ave,xt_cur,xt_std,xt_ske
     integer, dimension (its:ite), intent (in) ::                             &
           ierr
     real, dimension (its:ite,jts:jte,1:ensdim)                               &
         , intent(in   ) ::                                                   &
           xf_ens
!srf-----
!    real, dimension (its:ite,jts:jte)
     real, dimension (its:ite,jts:jte)                                     &
         , intent(inout) ::                                                   &
           APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                                 &
           APR_CAPMA,APR_CAPME,APR_CAPMI
!srf----

     real, dimension (its:ite,jts:jte)                                        &
         , intent(inout) ::                                                   &
           pr_gr,pr_w,pr_mc,pr_st,pr_as,                                      &
           pr_capma,pr_capme,pr_capmi

!--------------  firefly weights
   REAL, DIMENSION( its:ite , jts:jte ),                      &
         INTENT(IN) :: weight_GR,weight_W,weight_MC,weight_ST,weight_AS
!--------------

!
! local stuff
!
     real, dimension (its:ite , 1:maxens3 )       ::                          &
           x_ave,x_cur,x_std,x_ske
     real, dimension (its:ite , 1:maxens  )       ::                          &
           x_ave_cap


      integer, dimension (1:maxens3) :: nc1
      integer :: i,k
      integer :: num,kk,num2,iedt
      real :: a3,a4

      num=ensdim/maxens3
      num2=ensdim/maxens
      if(itest.eq.1)then
      do i=its,ite
       pr_gr(i,j) =  0.
       pr_w(i,j) =  0.
       pr_mc(i,j) = 0.
       pr_st(i,j) = 0.
       pr_as(i,j) = 0.
       pr_capma(i,j) =  0.
       pr_capme(i,j) = 0.
       pr_capmi(i,j) = 0.
      enddo
      endif

      do k=1,maxens
      do i=its,ite
        x_ave_cap(i,k)=0.
      enddo
      enddo
      do k=1,maxens3
      do i=its,ite
        x_ave(i,k)=0.
        x_std(i,k)=0.
        x_ske(i,k)=0.
        x_cur(i,k)=0.
      enddo
      enddo
      do i=its,ite
        xt_ave(i)=0.
        xt_std(i)=0.
        xt_ske(i)=0.
        xt_cur(i)=0.
      enddo
      do kk=1,num
      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0)then
        x_ave(i,k)=x_ave(i,k)+xf_ens(i,j,maxens3*(kk-1)+k)
        endif
      enddo
      enddo
      enddo
      do iedt=1,maxens2
      do k=1,maxens
      do kk=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0)then
        x_ave_cap(i,k)=x_ave_cap(i,k)                               &
            +xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+kk)
        endif
      enddo
      enddo
      enddo
      enddo
      do k=1,maxens
      do i=its,ite
        if(ierr(i).eq.0)then
        x_ave_cap(i,k)=x_ave_cap(i,k)/float(num2)
        endif
      enddo
      enddo

      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0)then
        x_ave(i,k)=x_ave(i,k)/float(num)
        endif
      enddo
      enddo
      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0)then
        xt_ave(i)=xt_ave(i)+x_ave(i,k)
        endif
      enddo
      enddo
      do i=its,ite
        if(ierr(i).eq.0)then
        xt_ave(i)=xt_ave(i)/float(maxens3)
        endif
      enddo
!
!--- now do std, skewness,curtosis
!
      do kk=1,num
      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0.and.x_ave(i,k).gt.0.)then
!       print *,i,j,k,kk,x_std(i,k),xf_ens(i,j,maxens3*(kk-1)+k),x_ave(i,k)
        x_std(i,k)=x_std(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**2
        x_ske(i,k)=x_ske(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**3
        x_cur(i,k)=x_cur(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**4
        endif
      enddo
      enddo
      enddo
      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0.and.xt_ave(i).gt.0.)then
        xt_std(i)=xt_std(i)+(x_ave(i,k)-xt_ave(i))**2
        xt_ske(i)=xt_ske(i)+(x_ave(i,k)-xt_ave(i))**3
        xt_cur(i)=xt_cur(i)+(x_ave(i,k)-xt_ave(i))**4
        endif
      enddo
      enddo
      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0.and.x_std(i,k).gt.0.)then
           x_std(i,k)=x_std(i,k)/float(num)
           a3=max(1.e-6,x_std(i,k))
           x_std(i,k)=sqrt(a3)
           a3=max(1.e-6,x_std(i,k)**3)
           a4=max(1.e-6,x_std(i,k)**4)
           x_ske(i,k)=x_ske(i,k)/float(num)/a3
           x_cur(i,k)=x_cur(i,k)/float(num)/a4
        endif
!       print*,'                               '
!       print*,'Some statistics at gridpoint i,j, ierr',i,j,ierr(i)
!       print*,'statistics for closure number ',k
!       print*,'Average= ',x_ave(i,k),'  Std= ',x_std(i,k)
!       print*,'Skewness= ',x_ske(i,k),' Curtosis= ',x_cur(i,k)
!       print*,'                               '

      enddo
      enddo
      do i=its,ite
!srf-tmp        if(ierr(i).eq.0.and.xt_std(i).gt.0.)then
        if(ierr(i).eq.0)then
           xt_std(i)=xt_std(i)/float(maxens3)
           a3=max(1.e-6,xt_std(i))
           xt_std(i)=sqrt(a3)
           a3=max(1.e-6,xt_std(i)**3)
           a4=max(1.e-6,xt_std(i)**4)
           xt_ske(i)=xt_ske(i)/float(maxens3)/a3
           xt_cur(i)=xt_cur(i)/float(maxens3)/a4
!       print*,'                               '
!       print*,'Total ensemble independent statistics at i =',i
!       print*,'Average= ',xt_ave(i),'  Std= ',xt_std(i)
!       print*,'Skewness= ',xt_ske(i),' Curtosis= ',xt_cur(i)
!       print*,'                               '
!
!  first go around: store massflx for different closures/caps
!
      if(itest.eq.1)then
             pr_gr(i,j) = .25*(x_ave(i,1)+x_ave(i,2)+x_ave(i,3)+x_ave(i,13))
             pr_w (i,j) = .25*(x_ave(i,4)+x_ave(i,5)+x_ave(i,6)+x_ave(i,14))
             pr_mc(i,j) = .25*(x_ave(i,7)+x_ave(i,8)+x_ave(i,9)+x_ave(i,15))
             pr_st(i,j) = .333*(x_ave(i,10)+x_ave(i,11)+x_ave(i,12))
             pr_as(i,j) = x_ave(i,16)
             pr_capma(i,j) = x_ave_cap(i,1)
             pr_capme(i,j) = x_ave_cap(i,2)
             pr_capmi(i,j) = x_ave_cap(i,3)
!-----
 	  !- training on closures
             if(training==1) then
          !  print*,'xt 1 = ', xt_ave(i),ffly_weights(i,j,1)*pr_gr(i,j)+ &
          !  		ffly_weights(i,j,2)*pr_w (i,j)+ &
	  !   		ffly_weights(i,j,3)*pr_mc(i,j)+ &
          !  		ffly_weights(i,j,4)*pr_st(i,j)+ &
	  !   		ffly_weights(i,j,5)*pr_as(i,j); call flush(6)

               xt_ave(i) =   weight_GR(i,j)*pr_gr(i,j)+ &
            	   	     weight_W (i,j)*pr_w (i,j)+ &
	     		     weight_MC(i,j)*pr_mc(i,j)+ &
            		     weight_ST(i,j)*pr_st(i,j)+ &
	     		     weight_AS(i,j)*pr_as(i,j)

	  !- training on CAPS
	     elseif(training==2) then

               xt_ave(i) =  weight_GR(i,j)*pr_capma(i,j)+ &
            		    weight_W (i,j)*pr_capme(i,j)+ &
	     		    weight_MC(i,j)*pr_capmi(i,j)

	     endif
!-----
!
!  second go around: store preciprates (mm/hour) for different closures/caps
!
        else if (itest.eq.2)then

 	  !- training on closures
         if(training==1) then

     	     APR_GR(i,j)=.25*(x_ave(i,1)+x_ave(i,2)+x_ave(i,3)+x_ave(i,13))*	  &
     			pr_gr(i,j) !+APR_GR(i,j)
     	     APR_W(i,j)=.25*(x_ave(i,4)+x_ave(i,5)+x_ave(i,6)+x_ave(i,14))*	  &
     		        pr_w(i,j)  !+APR_W(i,j)
     	     APR_MC(i,j)=.25*(x_ave(i,7)+x_ave(i,8)+x_ave(i,9)+x_ave(i,15))*	  &
     			pr_mc(i,j) !+APR_MC(i,j)
     	     APR_ST(i,j)=.333*(x_ave(i,10)+x_ave(i,11)+x_ave(i,12))*   &
     			pr_st(i,j) !+APR_ST(i,j)
     	     APR_AS(i,j)=x_ave(i,16)*			    &
     			pr_as(i,j) !+APR_AS(i,j)

	     xt_ave(i) =  weight_GR(i,j)*apr_gr(i,j)+ &
			  weight_W (i,j)*apr_w (i,j)+ &
			  weight_MC(i,j)*apr_mc(i,j)+ &
			  weight_ST(i,j)*apr_st(i,j)+ &
			  weight_AS(i,j)*apr_as(i,j)

	  !- training on CAPS
	  elseif(training==2) then

     	     apr_capma(i,j) = x_ave_cap(i,1)*pr_capma(i,j) !+apr_capma(i,j)
     	     apr_capme(i,j) = x_ave_cap(i,2)*pr_capme(i,j) !+apr_capme(i,j)
     	     apr_capmi(i,j) = x_ave_cap(i,3)*pr_capmi(i,j) !+apr_capmi(i,j)

	     !- the weights must add up 1
	     xt_ave(i) =  weight_GR(i,j)*pr_capma(i,j)+ &
                	  weight_W (i,j)*pr_capme(i,j)+ &
	        	  weight_MC(i,j)*pr_capmi(i,j)

	     !- only for output purposes
     	     APR_GR(i,j)=apr_capma(i,j)
     	     APR_W (i,j)=apr_capme(i,j)
     	     APR_MC(i,j)=apr_capmi(i,j)
     	     APR_ST(i,j)=0.
     	     APR_AS(i,j)=0.


	     apr_w (i,j) = max(0.,apr_w (i,j))
	     xt_ave(i)   =	  apr_w (i,j)

          else

	     APR_GR(i,j)=.25*(x_ave(i,1)+x_ave(i,2)+x_ave(i,3)+x_ave(i,13))*	  &
     			3600.*pr_gr(i,j) +APR_GR(i,j)
     	     APR_W(i,j)=.25*(x_ave(i,4)+x_ave(i,5)+x_ave(i,6)+x_ave(i,14))*	  &
     			3600.*pr_w(i,j) +APR_W(i,j)
     	     APR_MC(i,j)=.25*(x_ave(i,7)+x_ave(i,8)+x_ave(i,9)+x_ave(i,15))*	  &
     			3600.*pr_mc(i,j) +APR_MC(i,j)
     	     APR_ST(i,j)=.333*(x_ave(i,10)+x_ave(i,11)+x_ave(i,12))*   &
     			3600.*pr_st(i,j) +APR_ST(i,j)
     	     APR_AS(i,j)=x_ave(i,16)*			    &
     			3600.*pr_as(i,j) +APR_AS(i,j)

     	     APR_CAPMA(i,j) = x_ave_cap(i,1)*			       &
     			3600.*pr_capma(i,j) +APR_CAPMA(i,j)
     	     APR_CAPME(i,j) = x_ave_cap(i,2)*			       &
     			3600.*pr_capme(i,j) +APR_CAPME(i,j)
     	     APR_CAPMI(i,j) = x_ave_cap(i,3)*			       &
     			3600.*pr_capmi(i,j) +APR_CAPMI(i,j)

          endif


        endif
        endif
      enddo

   END SUBROUTINE massflx_stats

!-------------------------------------------------------
FUNCTION RS(P,T)

!     This function calculates the liquid saturation vapor mixing ratio as
!     a function of pressure and Kelvin temperature

implicit none
real esl,rs,x,t,p,c0,c1,c2,c3,c4,c5,c6,c7,c8,es
parameter (c0= .6105851e+03,c1= .4440316e+02,c2= .1430341e+01)
parameter (c3= .2641412e-01,c4= .2995057e-03,c5= .2031998e-05)
parameter (c6= .6936113e-08,c7= .2564861e-11,c8=-.3704404e-13)


!ES=610.78*EXP(17.269*(T-273.16)/(T-35.86))
!RS=.622*ES/(P-ES)
!return

x=max(-80.,t-273.16)

esl=c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
rs=.622*esl/(p-esl)

END FUNCTION RS
!-------------------------------------------------------

!-------------------------------------------------------
END MODULE module_cu_g3
