!---------------------------GRELL SHALLOW CUMULUS SCHEME---------------------------
subroutine cuparth_shal(CCATT,     &
     	  		mynum,     &	!1
     	  		mgmxp,     &	!2
     	  		mgmyp,     &	!3
     	  		mgmzp,     &	!4
     	  		m1,	   &	!5
     	  		m2,	   &	!6
     	  		m3,	   &	!7
     	  		ia,	   &	!8
     	  		iz,	   &	!9
     	  		!
     	  		ja,	   &	!10
     	  		jz,	   &	!11
     	  		i0,	   &	!12
     	  		j0,	   &	!13
     	  		maxiens,   &	!15
     	  		iens,	   &	!16
     	  		ngrid,     &	!17
     	  		ngrids_cp, &	!18
     	  		DTIME,     &	!19
     	  		!
     	  		time,	   &	!20
     	  		UA,	   &	!21
     	  		VA,	   &	!22
     	  		WA,	   &	!23
     	  		THETA,     &	!24
                        thetail,   &	!24
     	  		PP,	   &	!25
     	  		PI0,	   &	!26
     	  		DN0,	   &	!27
     	  		RV,	   &	!28
     	  		TKE,	   &	!29
     	  		!
     	  		TKMIN,     &	!30
     	  		rcp,       &
     	  		topt,	   &	!32
     	  		RTGT,	   &	!33
     	  		THT,	   &	!34
     	  		RTT,	   &	!35
     	  		PT,	   &	!36
     	  		OUTTEM,    &	!37
     	  		OUTRT,     &	!38
!     	  		SGRELL,    &	!39
!     	  		sgrell2_2d,&	!37
     	  		ierr4d,    &	!40
     	  		jmin4d,    &	!41
     	  		kdet4d,    &	!42
     	  		k224d,     &	!43
     	  		kbcon4d,   &	!44
     	  		ktop4d,    &	!45
     	  		kpbl4d,    &	!46
     	  		kstabi4d,  &	!47
     	  		kstabm4d,  &	!48
     	  		xmb4d,     &	!49
     	  		edt4d,     &	!50
     	  		zcup5d,    &	!51
     	  		pcup5d,    &	!52
     	  		enup5d,    &	!53
     	  		endn5d,    &	!54
     	  		deup5d,    &	!55
     	  		dedn5d,    &	!56
     	  		zup5d,     &	!57
     	  		zdn5d,     &	!58
     	  		prup5d,    &	!59
     	  		clwup5d,   &	!60
     	  		tup5d,     &	!61
     	  		upmf,	   &	!62
     	  		xierr,     &	!63
     	  		xktop,     &	!64
     	  		xkbcon,    &	!65
     	  		xk22,	   &	!66
     	  		xierr_dp,  &	!67
     	  		confrq,    &
     	  		frqanl,    &
     	  		deltaxy2,  &
     	  		patch_area, npat,level,trigg,autoconv    )

  use mem_grell_param, only : &
           		     maxens =>maxens_sh,    & !INTENT(IN)
           		     maxens2=>maxens2_sh,   & !INTENT(IN)
           		     maxens3=>maxens3_sh,   & !INTENT(IN)
           		     icoic  =>icoic_sh	  !,& !INTENT(IN)
           		     !ensdim,		      !INTENT(IN)

  use Phys_const, only: rgas,cp,rm,p00,tcrit,g,cpor,pkdcut,xl

  use mem_scratch2_grell_sh
  use mem_varinit, only: nudlat
  use mem_grid, only: nxtnest, initial, jdim
  use node_mod, only: ibcon

  implicit none
  integer mgmxp,mgmyp,mgmzp,ngrid,ngrids_cp
  integer maxiens,ensdim,iens
  !
  ! maxens =3    ensemble one on mbdt
  ! maxens2=1    ensemble two on precip efficiency
  ! maxens3=10   ensemble three done in cup_forcing
  ! ensdim = maxiens*maxens*maxens2*maxens3 : ensemble dimension
  !
  !--- RAMS arrays
  integer, intent(IN):: CCATT,m1,m2,m3,ia,iz,ja,jz,i0,j0,mynum,npat,level,trigg,autoconv
  real :: time,tkmin,confrq,frqanl,dti,deltaxy2,tscl_KF
  real, dimension(m1,m2,m3) :: ua,va,wa
  real, dimension(m1,m2,m3) :: theta,pp,pi0,dn0,rv,tht,rtt, &
                               pt,outtem,outrt,tke,rcp,thetail
  real, dimension(m2,m3) :: topt,rtgt,upmf,xierr,xktop,xkbcon,xk22, &
                            xierr_dp

  real, dimension (mgmxp,mgmzp) :: dn0_2


  real, dimension(m2,m3,npat)  :: patch_area




  !-------salva parametros da param. para uso no transporte convectivo:
  integer, dimension(mgmxp,mgmyp,maxiens,ngrids_cp) :: ierr4d,jmin4d, &
                                  kdet4d,k224d,kbcon4d,ktop4d,kpbl4d, &
                                  kstabi4d,kstabm4d
  real,dimension(mgmxp,mgmyp,maxiens,ngrids_cp) :: xmb4d,edt4d
  real,dimension(mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp) :: enup5d,endn5d, &
                                                         deup5d,dedn5d, &
							 zup5d,zdn5d,   &
							 prup5d,clwup5d,&
							 tup5d,zcup5d,  &
							 pcup5d

  !------------------------------------- variaveis locais:
  integer istart,iend,i,j,k,mix,mjx,mkx,kr,ipr,jpr,j1,j2     ! kk,m
  real dtime,cpdTdt,exner
  real :: r_liq, r_sol,tempk,fxc,dfxcdt

  !Ensemble dimension
  !srf- 04-fev-2003 - mudanca no calculo DO 'ensdim'
  integer, parameter :: iens_tmp =1
  !integer, save :: iens_tmp
  !data iens_tmp /1/
  !----------------------------------------------------------------------
  ! iens_tmp substitui o parametro iens ate que se resolva a unIFicacao DO shallow e deep
  ! numa so rotina.
  !(apesar de haver espectro de nuvens (shallow e deep) estas sao tratadas separadamente em
  ! rotinas dIFerentes isto eh maxiens seria igual a 1. Porem no transporte convectivo
  ! o tratamento eh com maxiens =2.
  !ensdim=maxiens*maxens*maxens2*maxens3
   ensdim=      1*maxens*maxens2*maxens3


  !--srf: por  enquanto fixe aqui o fechamento
  !       modifique depois na rotina rams_mem_alloc.f90
  ! icoic is used for choice a specific closure for shallow cumulus
  ! icoic = 0 -> ensemble (all closures)
  ! icoic = 1 -> Grell
  ! icoic = 4 -> boundary layer Quasi-Equilibrium (former Arakawa-Schubert)
  ! icoic = 8 -> like Fritsch Chappel or Kain Fritsch
  icoic = 4
  !

  ISTART = ia
  IEND   = iz
  j1     = ja
  j2     = jz
  MKX    = m1 - 1    !MKX nao deve ser igual a m1
  MIX    = m2
  MJX    = m3

!  If variable initialization, on a coarse grid, not in a global simulation,
!  do not run convective parameterization in the lateral boundary region.

  IF (INITIAL == 2 .AND. nxtnest(ngrid) == 0 ) THEN

     if (iand(ibcon,1) > 0) ISTART = 1  + (nudlat)
     if (iand(ibcon,2) > 0) IEND   = m2 - (nudlat)
     if (iand(ibcon,4) > 0) j1 = 1 + (nudlat) * jdim
     if (iand(ibcon,8) > 0) j2 = max (1,m3 - (nudlat))

  ENDIF

  !---- Coordenadas para escrita ascii
  !      ipr=19 - i0
  !      jpr=19 - j0
  ipr=0 - i0
  jpr=0 - j0
  !srf- out/2004
  !- for ensemble average
  if(frqanl>0.) then
   dti = confrq/frqanl
  else
   dti = 0.
  endif
  ! A. Betts: suggestion for the KF timescale < DELTAX  / 25 m/s
  tscl_KF =  sqrt(deltaxy2) / 25.  !units: sec

  ! Loop externo : j
  do J=j1,j2

     do I = ISTART,IEND
        aa0(i)     = 0.
        xland(i,j) = patch_area(i,j,1) ! land < 1 /water=1 flag
        if(xland(i,j) < 0.95) xland(i,j) = 0.
     enddo
     !------- Transfere valores do RAMS para o esquema
     do K=1,MKX
        kr = K + 1 ! nivel K da grade do Grell corresponde ao nivel K + 1 do RAMS
        do I = ISTART,IEND
           z1(I)= topt(i,j)
           PSUR(I) = .5*( ((pp(1,i,j)+pi0(1,i,j))/cp)**cpor*p00 + &
                          ((pp(2,i,j)+pi0(2,i,j))/cp)**cpor*p00 )*1.e-2 ! Pressure in mbar
           PO(I,K) = ((pp(kr,i,j)+pi0(kr,i,j))/cp)**cpor*p00*1.e-2      ! Pressure in mbar
           US(I,K) = .5*( ua(kr,i,j) + ua(kr,i-1,j) )
           VS(I,K) = .5*( va(kr,i,j) + va(kr,i,j-1) )
           OMEG(I,K)   = -g*dn0(kr,i,j)*.5*( wa(kr,i,j)+wa(kr-1,i,j) )
           T(I,K)  = theta(kr,i,j)*(pp(kr,i,j)+pi0(kr,i,j))/cp
           Q(I,K)  = rv(kr,i,j)

	   ! fazendo a densidade em duas dimensoes
	   dn0_2(i,k) = dn0(kr,i,j)


           !variables for PBL top height
           TKEG(I,K) = TKE(kr,i,j)
           RCPG(I,K) = RCP(kr,i,j)
           !Calcula tendencia projetada na temperatura em funcao
           !das tendencias de theta e PI : cp*T=Pi*Theta
           exner= pp(kr,i,j)+pi0(kr,i,j)
!          cpdTdt  = exner*tht(kr,i,j)    !assuminDO PT(KR,I,J) << exner*THT(KR,I,J)/theta
           cpdTdt= exner*tht(kr,i,j) + theta(kr,i,j)*pt(kr,i,j)

           !Temperatura projetada se a conveccao nao ocorrer
           TN(I,K) = T(I,K) + ( cpdTdt/cp )*dtime
           !Umidade projetada se a conveccao nao ocorrer
           QO(I,K) = Q(I,K) +   rtt(kr,i,j)*dtime
           !Atribuicoes do esquema
           P(I,K)  = PO(I,K)
           !srf PSUR(I) = 0.5*(PO(I,1)+PO(I,2))
           if(TN(I,K).lt.200.)    TN(I,K) = T(I,K)
           if(QO(I,K).lt.1.E-08)  QO(I,K) = 1.E-08

	   !- srf-25062007- for bound-layer-quasi-equil closure (BLQE)
	   dhdt(i,k) = cpdTdt+xl*rtt(kr,i,j) ! large-scale/rad/turb moist static energy change

           OUTT(I,K) = 0.   !Tendencia no campo de temperatura associada aos cumulus
           OUTQ(I,K) = 0.   !Tendencia na razao de mist. de vapor d'agua assoc. aos cumulus
        enddo
     enddo
     !---  CUMULUS PARAMETERIZATION
     !      iens = 2              !--- shallow convection
     !      iens_tmp=1 substitui iens=2, ate a unificacao
     call CUP_enss_shal(ccatt,ngrid,mynum,m1,m2,m3,i0,j0,ipr,jpr, &
          mgmxp,mgmyp,mgmzp,maxiens,maxens,maxens2,maxens3,ensdim, &
          icoic,j,iens_tmp,istart,iend,mix,mjx,mkx,xland,z1, &
          aa0, t, q, tn, qo, po, p, outt, outq, dtime, &
          psur,tcrit,time,omeg, &
          TKEG,RCPG,tkmin, &
            ierr4d(1,j,iens,ngrid),   jmin4d(1,j,iens,ngrid), &
            kdet4d(1,j,iens,ngrid),    k224d(1,j,iens,ngrid), &
           kbcon4d(1,j,iens,ngrid),   ktop4d(1,j,iens,ngrid), &
            kpbl4d(1,j,iens,ngrid),                           &
          kstabi4d(1,j,iens,ngrid), kstabm4d(1,j,iens,ngrid), &
             xmb4d(1,j,iens,ngrid),    edt4d(1,j,iens,ngrid), &
          zcup5d(1,1,j,iens,ngrid), pcup5d(1,1,j,iens,ngrid), &
          enup5d(1,1,j,iens,ngrid), endn5d(1,1,j,iens,ngrid), &
          deup5d(1,1,j,iens,ngrid), dedn5d(1,1,j,iens,ngrid), &
           zup5d(1,1,j,iens,ngrid),  zdn5d(1,1,j,iens,ngrid), &
          prup5d(1,1,j,iens,ngrid),clwup5d(1,1,j,iens,ngrid), &
           tup5d(1,1,j,iens,ngrid),                           &
          upmf,xierr,xktop,xkbcon,xk22,xierr_dp, dti,tscl_KF,dhdt,trigg, iens,autoconv    )


     do I = ISTART,IEND
	if(int(xierr(i,j)) /= 0) cycle
         do K=1,MKX-1
            kr = K + 1
           !      Converte tendencia da temperatura (OUTT) em tendencia de theta (OUTTEM)
           !      cp*T=Pi*Theta => cp dT/dt = Theta*dPi/dt + Pi*dTheta/dt,
           !      Exner's function = pp(kr,i,j)+pi0(kr,i,j)
           !if(OUTT(I,K) /= 0.) then!-srf 26062012
            exner          =   pp(kr,i,j) + pi0(kr,i,j)
            ! tendencia do theta devida a conv profunda
            !-srf 26jun2012 -ate resolver a correta convers�o de dTemp/dt para dTheta/dt
	    outtem(kr,i,j) = cp/exner * OUTT(I,K) !- theta(kr,i,j)*pt(kr,i,j)/exner
           !endif!-srf 26062012

           outrt(kr,i,j) =   OUTQ(I,K)   ! tendencia DO Rtotal devida aos cumulus
           !
        enddo
     enddo
  enddo     ! loop externo - j -
           !
end subroutine CUPARTH_shal

!----------------------------------------------------------------------

subroutine CUP_enss_shal(ccatt,ngrid,mynum,m1,m2,m3,i0,j0,ipr,jpr, &
     	  		 mgmxp,mgmyp,mgmzp,maxiens,maxens,maxens2,maxens3,ensdim, &
     	  		 icoic,j,iens,ISTART,IEND,mix,mjx,mkx,xland,Z1, &
     	  		 AAEQ, T, Q, TN, QO, PO, P, OUTT, OUTQ, DTIME, &
     	  		 PSUR,TCRIT,time,omeg, &
     	  		 TKEG,RCPG,tkmin, &
     	  		 ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d,kpbl4d, &
     	  		 kstabi4d,kstabm4d,xmb4d,edt4d, &
     	  		 zcup5d,pcup5d,enup5d,endn5d,deup5d,dedn5d,zup5d,zdn5d, &
     	  		 prup5d,clwup5d,tup5d, &
     	  		 upmf,xierr,xktop,xkbcon,xk22,xierr_dp, &
			 	 dti,tscl_KF,dhdt,trigg,actual_iens,autoconv)

 use mem_scratch3_grell_sh, only: &
       CD, &
       AA0, &
       AA1, &
       IERR, &
       CAP_MAX_INCREMENT, &
       CAP_MAX, &
       MBDT_ENS, &
       XF_ENS, &
       Z, &
       QES, &
       HE, &
       HES, &
       ZO, &
       QESO, &
       HEO, &
       HESO, &
       QES_CUP, &
       Q_CUP, &
       HE_CUP, &
       HES_CUP, &
       Z_CUP, &
       P_CUP, &
       GAMMA_CUP, &
       T_CUP, &
       QESO_CUP, &
       QO_CUP, &
       HEO_CUP, &
       HESO_CUP, &
       ZO_CUP, &
       PO_CUP, &
       GAMMAO_CUP, &
       TN_CUP, &
       KBMAX, &
       KZI, &
       K22, &
       KBCON, &
       HKBO, &
       HKB, &
       HC, &
       DBY, &
       HCO, &
       DBYO, &
       KTOP, &
       ZU, &
       ZUO, &
       QC, &
       QRC, &
       PW, &
       PWAV, &
       QCO, &
       QRCO, &
       PWO, &
       PWAVO, &
       DELLAT_ENS, &
       DELLAQ_ENS, &
       DELLAH, &
       DELLAQ, &
       XAA0_ENS, &
       DELLAT, &
       XHE, &
       XQ, &
       DELLAT, &
       XT, &
       XQ, &
       XHE, &
       XQ, &
       XT, &
       XQ, &
       XZ, &
       XQES, &
       XHES, &
       XQES_CUP, &
       XQ_CUP, &
       XHE_CUP, &
       XHES_CUP, &
       XZ_CUP, &
       XT_CUP, &
       XHKB, &
       XHC, &
       XDBY, &
       XZU, &
       XQC, &
       XQRC, &
       XPW, &
       XPWAV, &
       XAA0, &
       XMB, &
       XFAC1, &
       KSTABI, &
       KSTABM
  use Phys_const, only: g

  implicit none
  integer :: ccatt,ngrid,mynum,i0,j0,ipr,jpr,m1,m2,m3
  integer :: maxiens,maxens,maxens2,maxens3,ensdim,actual_iens,autoconv
  integer :: mix,mjx,mkx,mgmxp, mgmyp, mgmzp
  integer :: nens,iens,iedt,icoic,trigg       !,nall,nens3,izero,ktau
  real ::    time,tkmin,dti,tscl_KF
  !--- Input variables -----------------------------
  real :: mconv(mgmxp),Z1(mgmxp),AAEQ(mgmxp),PSUR(mgmxp)
  real :: T(mgmxp,mgmzp), Q(mgmxp,mgmzp),  TN(mgmxp,mgmzp),  QO(mgmxp,mgmzp), &
       P(mgmxp,mgmzp),PO(mgmxp,mgmzp),omeg(mgmxp,mgmzp),TKEG(mgmxp,mgmzp), &
       RCPG(mgmxp,mgmzp),dhdt(mgmxp,mgmzp)

  real, dimension(m1,m2,m3) :: dn0

  real xland(mgmxp,mgmyp)

  !-------Salva parametros da CUP para uso no transporte convectivo:
  integer,dimension(mgmxp) :: ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d
  integer,dimension(mgmxp) :: kstabi4d,kstabm4d,kpbl4d
  real,   dimension(mgmxp) :: xmb4d,edt4d
  real,   dimension(mgmzp,mgmxp) :: zcup5d,pcup5d,enup5d,endn5d,deup5d,dedn5d,&
                                    zup5d,zdn5d,prup5d,clwup5d,tup5d

  real,   dimension(mgmxp,mgmzp) :: dn0_2

  !------Variables saved in RAMS Analisys
  ! use (m1,m2,m3) para dimensionar os vetores que sao escritos nas analises DO RAMS
  real, dimension(m2,m3) :: upmf,xierr,xktop,xkbcon,xk22,xierr_dp

  integer ::fquasi
  integer ::kstart,kstart_way
  integer ::I,J,K,ISTART,IEND
  real ::day,tcrit,dtime
  real ::dellaqsum,dp,mbdt
!  real cap_max_increment,cap_maxs,fdiur,tke_start
  real ::cap_maxs,fdiur,tke_start,tkeminnot
  !
  !--- New entrainment/detrainment related stuff --------------------
  !real mentr_rate,entr_rate,radius,depth_min,zkbmax

  real :: entr_rate,radius,depth_min,zkbmax
  real :: mentr_rate
  real :: mentr_rate2(mgmxp,mgmzp),wu(mgmxp,mgmzp),ccn(mgmxp)  !fazendo o entranhamento como matriz (C. Silva)

  !--- Output variables ----------------------------
  real :: outt(mgmxp,mgmzp),outq(mgmxp,mgmzp)

  day=86400.
  !--- specify entrainmentrate and detrainmentrate
  !radius=666. !same used at ECMWF model
  radius=200.  !same used at WRF/G3d model
  !--- gross entrainment rate (these may be changed later on in the
  !--- program, depending what your detrainment is!!)
  entr_rate=.2/radius
  !--- entrainment of mass

  mentr_rate =entr_rate  ! utilizar a matriz (C. Silva)

  !--- initial detrainment rates
  ! strong lateral mixing, entrtainment=dedtrainment,
  ! constant mass flux with height
  do k=1,mkx
     do i=istart,iend
         mentr_rate2(i,k) = mentr_rate
	 wu(i,k) = 0.   ! only for shallow conv
     enddo
  enddo
  ccn(istart:iend)= 300. ! only for shallow conv


  do k=1,mkx
     do i=istart,iend
        cd(i,k)   = entr_rate
     enddo
  enddo

  !--- minimum depth (m), clouds must have
  depth_min=50.


  do I=ISTART,IEND
     aa0(i)=0.
     aa1(i)=0.

     !srf -fev 2003
     !don't permite shallow if there is deep convection
     !if(int(xierr_dp(i,j)) .eq. 0 ) then
     !if(int(xierr_dp(i,j)) .gt. 10000 ) then
         !print*,'deep on => shallow off'
     !    ierr(i)  = 20
     !else
         ierr(i)  = 0
        xierr(i,j)= 0.
     !endif
  enddo

  !--- cap_maxs trigger parameter
  cap_max_increment(:)=0.
  cap_maxs=120.

  !---  cap_max initialization
  tkeminnot = 0.08
  tke_start = 0.05
  do i=istart,iend
  ! ciclo diurno off
    !cap_max(i)=cap_maxs
  !--- diurnal cycle on cap_maxs:
     fdiur = 1.0
     if( tkeg(i,1) - tkeminnot <  tke_start .and. tkeg(i,1) - tkeminnot > tkmin) fdiur = 0.5
     cap_max(i)= 70. + fdiur * max(0.,cap_maxs - 70.) &
     	      *( tkeg(i,1) - tkmin - tkeminnot) / max(tkeg(i,1), tkmin + tkeminnot)
  !
  !---srf out/2004
  !-  xland alterando o capmax
  ! if(xland(i,j) > 0.95) cap_max(i)=cap_maxs+50.
  end do


  !--- max height(m) above ground where shallow clouds can originate
  zkbmax=4000.

  do nens=1,maxens
     mbdt_ens(nens)=(float(nens)-3. )*dtime*1.e-3+dtime*5.E-03
    !mbdt_ens(nens)=(float(nens)-1.5)*dtime*2.e-3+dtime*5.E-03
  enddo

  !--- environmental conditions, FIRST HEIGHTS
  do i=istart,iend
     if(ierr(i).ne.20)then
        do k=1,maxens*maxens2*maxens3
           xf_ens(  i,j,(iens-1)*maxens*maxens2*maxens3+k)= 0.
        enddo
     endif
  enddo

  !--- calculate moist static energy, heights, qes
  call cup_env(j,ipr,jpr,z,qes,he,hes,t,q,p,z1,mix,mgmxp, &
       mkx,mgmzp,istart,iend,psur,ierr,tcrit,0)
  call cup_env(j,ipr,jpr,zo,qeso,heo,heso,tn,qo,po,z1,mix,mgmxp, &
       mkx,mgmzp,istart,iend,psur,ierr,tcrit,0)

  !--- environmental values on cloud levels
  call cup_env_clev(j,ipr,jpr,t,qes,q,he,hes,z,p,qes_cup, &
       q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
       mix,mgmxp,mkx,mgmzp,istart,iend,ierr,z1)
  call cup_env_clev(j,ipr,jpr,tn,qeso,qo,heo,heso,zo,po, &
       qeso_cup,qo_cup,heo_cup,heso_cup,zo_cup,po_cup, &
       gammao_cup,tn_cup,psur,mix,mgmxp,mkx,mgmzp,istart, &
       iend,ierr,z1)

  do i=istart,iend
     if(ierr(i).eq.0)then
        do k=1,mkx
           if(zo_cup(i,k).gt.zkbmax+z1(i))then
              kbmax(i)=k
              ! PRINT*,'kbmax', zkbmax,z1(i),zo_cup(i,k),kbmax(i)
              exit
           endif
        enddo
     endif
  enddo

  !--- DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22
  !
  !Gr-dec2002
  !   here, ,should just take top of PBL for shallow clouds, also maybe for
  !   deep clouds in tropics: k22=level(pbltop)
  !
  !srf-fev2003
  kstart_way = 1

  if (kstart_way == 1) then
     !	  New way to define K22
     !----   Determine PBL top using TKE (TKEG) and liquid water mixing ratio (RCPG)
     call get_zi(mix,mgmxp,mkx,mgmzp,istart,iend,j,ierr,kzi,&
          TKEG,RCPG,zo,z1,tkmin)

     do i=istart,iend
        !srf-14-fev-2003
        !A segunda forma produz uma cobertura de shallow mais representativa
        !Em caso de alta resolucao vertical tente a versao original (k22=kzi)
        !	     IF(ierr(i).eq.0) k22(i) = kzi(i)
        if(ierr(i).eq.0) k22(i) = max(2, kzi(i) - 1)
        !	    IF(ierr(i).eq.0) PRINT*,k22(i)
     enddo
!srf - out 2004
!srf- TESTAR NO FUTURO -
!      new way to define cap_max diurnal cycle
!	DO  I=ISTART,IEND
!	  IF(ierr(I).eq.0.)THEN
!	  tkemax(i) = 0.
!	   do k=1,kzi(i)
!		   tkemax(i) = max(tkemax(i),tkeg(i,k))
!	  enddo
!	  if(tkemax(i) < 1.5 )  cap_max(i) = 75.
!	  if(tkemax(i) < 0.5 )  cap_max(i) = 25
!	 endif
!      enddo
  else
     !---- old way to define k22
     kstart=2
     call maximi(heo_cup,mix,mgmxp,mkx,mgmzp,kstart,kbmax,k22,istart,iend,ierr)
  endif

  do I=ISTART,IEND
     if (ierr(I).eq.0) then
        if(K22(I).ge.KBMAX(i))ierr(i)=2
     endif
  enddo
  !if (j.eq.jpr) print *,ierr(ipr),k22(ipr),kbmax(ipr)

  !--- Determine the level of convective cloud base  - kbcon
  call cup_kbcon_catt(cap_max_increment,1,k22,kbcon,heo_cup,heso_cup,&
       hkbo,kzi,mix,mgmxp,mkx,mgmzp,istart,iend,ierr,kbmax,po_cup,cap_max,&
       j)

   !- srf: out2004
   !- TESTAR NO FUTURO
   !do i=istart,iend
   ! if(ierr(i).eq.0) hkb(i)=hkbo(i)
   !enddo

  !--PRINT-----srf
  !	 DO I=ISTART,IEND
  !	  IF(i.eq.ipr.and.j.eq.jpr) THEN
  !	  PRINT*,'-------------------------------'
  !	  PRINT*,'MYNUM I J=',mynum,i,j
  !	  PRINT*,'k22-kbcon-ierr=',k22(i),kbcon(i),ierr(i)
  !	  PRINT*,'-------------------------------'
  !	  ENDIF
  !	 ENDDO
  !--PRINT-----srf

  !--- Calculate incloud moist static energy
  call cup_up_he(k22,hkb,z_cup,cd,mentr_rate2,he_cup,hc,mix,mgmxp,mkx, &
       mgmzp,kbcon,ierr,istart,iend,dby,he,hes_cup)
  call cup_up_he(k22,hkbo,zo_cup,cd,mentr_rate2,heo_cup,hco, &
       mix,mgmxp,mkx,mgmzp,kbcon,ierr,istart,iend,dbyo, &
       heo,heso_cup)

  !srf--PRINT--------------------
  !	DO I=ISTART,IEND
  !	 IF(i.eq.ipr.and.j.eq.jpr) THEN
  !	  PRINT*,'-----------K T O P-1-----------'
  !	  PRINT*,'MYNUM I J=',mynum,i,j
  !	  PRINT*,'k22-kbcon-ierr-Ktop=',k22(i),kbcon(i),ierr(i),ktop(i)
  !	   DO K=KBCON(I)+1,MKX-2
  !	    PRINT*,'k-heo_cup-dbyo=',k,heo_cup(i,k),dbyo(i,k)
  !	   ENDDO
  !	  PRINT*,'-------------------------------'
  !	 ENDIF
  !	ENDDO
  !srf--PRINT--------------------

  !--- Determine cloud top - KTOP
  call cup_ktop(1,dbyo,kbcon,ktop,mix,mgmxp,mkx,mgmzp,istart,iend,ierr)

  !--- Normalized updraft mass flux profile
  call cup_up_nms(zu,z_cup,mentr_rate2,cd,kbcon,ktop,mix,mgmxp, &
       mkx,mgmzp,istart,iend,ierr,k22)
  call cup_up_nms(zuo,zo_cup,mentr_rate2,cd,kbcon,ktop,mix,mgmxp, &
       mkx,mgmzp,istart,iend,ierr,k22)

  !--- Calculate moisture properties of updraft
  !call cup_up_moisture(ierr,z_cup, qc, qrc, pw, pwav,kbcon,ktop,mix,&
  !     mgmxp,mkx,mgmzp,istart,iend,cd,dby,mentr_rate2, &
  !     q,GAMMA_cup, zu, qes_cup, k22,q_cup)
  !call cup_up_moisture(ierr,zo_cup,qco,qrco,pwo,pwavo,kbcon,ktop,mix, &
   !    mgmxp,mkx,mgmzp,istart,iend,cd,dbyo,mentr_rate2, &
   !    qo,GAMMAo_cup,zuo,qeso_cup,k22,qo_cup)


  !--- Calculate moisture properties of updraft

  call cup_up_moisture(ierr,z_cup, qc, qrc, pw, pwav,kbcon,ktop,mix,&
       mgmxp,mkx,mgmzp,istart,iend,cd,dby,mentr_rate2, &
       q,GAMMA_cup, zu, qes_cup, k22,q_cup, wu, dn0_2,ccn,trigg,actual_iens,autoconv )


  call cup_up_moisture(ierr,zo_cup,qco,qrco,pwo,pwavo,kbcon,ktop,mix, &
       mgmxp,mkx,mgmzp,istart,iend,cd,dbyo,mentr_rate2, &
       qo,GAMMAo_cup,zuo,qeso_cup,k22,qo_cup,  wu, dn0_2,ccn,trigg,actual_iens,autoconv )

  !--- Calculate workfunctions for updrafts
  call cup_up_aa0(aa0,z, zu, dby, GAMMA_CUP, t_cup,kbcon,ktop,mix, &
       mgmxp,mkx,mgmzp,istart,iend,ierr)
  call cup_up_aa0(aa1,zo,zuo,dbyo,GAMMAo_CUP,tn_cup,kbcon,ktop,mix,&
       mgmxp,mkx,mgmzp,istart,iend,ierr)

  do i=istart,iend
!srf - print------------------
!	   if(j==7 .and. i==48) then
!	     !if( k==1)
!	     print*,'BRAMS: AA1'
!	     print*,AA1(I),kbcon(i),ktop(i)
!	   endif
!srf - print------------------

     if(ierr(i).eq.0) then
        if(aa1(i).eq.0.) then
           ierr(i)=17
        endif
     endif
  enddo

  !srf - Big loop starts here!
  iedt=1
  do k=1,mkx
     do i=istart,iend
        dellat_ens(i,k,iedt)=0.
        dellaq_ens(i,k,iedt)=0.
     enddo
  enddo

  !PRINT-----------------------
  !	IF(j.eq.jpr)THEN
  !	  i=ipr
  !	  IF(ierr(i).eq.0.or.ierr(i).eq.3)THEN
  !	  PRINT *,k22(I),kbcon(i),ktop(i)
  !	  PRINT *,aa0(i),aa1(i)
  !	  DO k=1,mkx
  !	    PRINT *,z(i,k),he(i,k),hes(i,k)
  !	  ENDDO
  !	  DO k=1,ktop(i)+1
  !	    PRINT *,zu(i,k),pw(i,k)
  !	  ENDDO
  !	  ENDIF
  !	ENDIF
  !PRINT-----------------------

  !--- changes per unit mass from shallow convection
  call cup_dellas_shallow(ierr,zo_cup,po_cup,heo,mix,mgmxp,mkx,mgmzp, &
                          istart,iend,dellah,1,j,zuo,cd,hco,ktop,k22, &
                          kbcon,mentr_rate2,heo_cup,ipr,jpr,'shallow')

  call cup_dellas_shallow(ierr,zo_cup,po_cup,qo ,mix,mgmxp,mkx,mgmzp, &
                          istart,iend,dellaq,2,j,zuo,cd,qco,ktop,k22, &
                          kbcon,mentr_rate2, qo_cup,ipr,jpr,'shallow')

  !srf-----PRINT-------
  !do i=istart,iend
  !   if(ierr(i).eq.0)then
  !	 if(j.eq.jpr.and.i.eq.ipr) then
  !	    dellaqsum  =0.
  !	    do k=1,mkx
  !	       dp=-100.*(p_cup(i,k+1)-p_cup(i,k))
  !	       dellaqsum  =  dellaqsum +dellaq(i,k) *dp/g
  !	       !   IF(k.eq.1) THEN
  !	       !	   WRITE(6,'(a1,78a1)') ' ',('-',m=1,78)
  !	       !    PRINT*,'i k Z_cup P_cup dp dellaq dellaq dellaqsum'
  !	       !   ENDIF
  !	       !	  WRITE(6,'(2i4,5F12.4,3f16.4)') i,k,Z_cup(i,k),p_cup(i,k)
  !	       !    &	       ,dp,
  !	       !		1.e+3*86400.*dellaq(i,k)
  !	       !    &	       ,1.e+3*86400.*dellaqsum,
  !	       !		1.e+3*86400.*dellaqsum
  !	    enddo
  !	 endif
  !   endif
  !enddo
  !srf-----PRINT-------

  !--- Using dellas, calculate changed environmental profiles
  do nens=1,maxens
     mbdt=mbdt_ens(nens)
     do i=istart,iend
        xaa0_ens(i,nens)=0.
     enddo
     do k=1,mkx-1
        do i=istart,iend
           dellat(i,k)=0.
           if(ierr(i).eq.0)then

              xhe   (i,k) =            dellah(i,k)*mbdt + heo(i,k)
              xq    (i,k) =	       dellaq(i,k)*mbdt +  qo(i,k)
              dellat(i,k) =(1./1004.)*(dellah(i,k)-2.5e06*dellaq(i,k))
              xt    (i,k) =	       dellat(i,k)*mbdt +  tn(i,k)
              if(xq(i,k).le.0.) xq(i,k)=1.e-08
              !if(i.eq.ipr.and.j.eq.jpr)then
              !   print *,k,DELLAH(I,K),DELLAQ(I,K),DELLAT(I,K)
              !endif

	   endif
        enddo
     enddo

     do i=istart,iend
        if(ierr(i).eq.0)then
           xhe(i,mkx)  = heo(i,mkx)
           xq (i,mkx)  =  qo(i,mkx)
           xt (i,mkx)  =  tn(i,mkx)
           if(xq(i,mkx).le.0.) xq(i,mkx)  =  1.e-08
        endif
     enddo

     !--- Calculate moist static energy, heights, qes
     call cup_env(j,ipr,jpr,xz,xqes,xhe,xhes,xt,xq,po,z1,mix,mgmxp, &
                  mkx,mgmzp,istart,iend,psur,ierr,tcrit,2)

     !--- Environmental values on cloud levels

     call cup_env_clev(j,ipr,jpr,xt,xqes,xq,xhe,xhes,xz,po,xqes_cup,    &
                       xq_cup,xhe_cup,xhes_cup,xz_cup,po_cup,gamma_cup, &
		       xt_cup,psur,mix,mgmxp,mkx,mgmzp,istart,iend,ierr,z1)

     !**************************** Static Control ************
     !--- Moist static energy inside cloud
     do i=istart,iend
        if(ierr(i).eq.0)then
           xhkb(i)=xhe(i,k22(i))
        endif
     enddo
     call cup_up_he(k22,xhkb,xz_cup,cd,mentr_rate2,xhe_cup,xhc,  &
                    mix,mgmxp,mkx,mgmzp,kbcon,ierr,istart,iend, &
		    xdby,xhe,xhes_cup)

     !--- Normalized mass flux profile
     call cup_up_nms(xzu,xz_cup,mentr_rate2,cd,kbcon,ktop, mix,mgmxp,&
                     mkx,mgmzp,istart,iend,ierr,k22)

     !--- Moisture updraft
   !  call cup_up_moisture(ierr,xz_cup,xqc,xqrc,xpw,xpwav,kbcon,ktop,&
   !                       mix,mgmxp,mkx,mgmzp,istart,iend,cd,xdby,  &
   !                       mentr_rate2,xq,GAMMA_cup,xzu,xqes_cup,k22,xq_cup)

     call cup_up_moisture(ierr,xz_cup,xqc,xqrc,xpw,xpwav,kbcon,ktop,&
                          mix,mgmxp,mkx,mgmzp,istart,iend,cd,xdby,  &
                          mentr_rate2,xq,GAMMA_cup,xzu,xqes_cup,k22, &
			  xq_cup,  wu, dn0_2,ccn,trigg,actual_iens,autoconv )


     !--- Workfunctions for updraft
     call cup_up_aa0(xaa0,xz,xzu,xdby,GAMMA_CUP,xt_cup,kbcon, &
                     ktop,mix,mgmxp,mkx,mgmzp,istart,iend,ierr)

     !srf-feb-2003
     do i=istart,iEND
        if(ierr(i).eq.0) xaa0_ens(i,nens)=xaa0(i)
     enddo

  enddo

  !--------- LARGE SCALE FORCING  -----------------------------------
  call cup_forcing_ens_shal(aa0,aa1,xaa0_ens,mbdt_ens,dtime,xmb, &
       ierr,mix,mgmxp,mjx,mgmyp,mkx,mgmzp,istart,iend,xf_ens,j, &
       'shallow',xland,maxens,iens,maxens2,ipr,jpr,maxens3, &
       ensdim,p_cup,ktop,icoic,iedt,&!sgrell2_2d,
        m2,m3, tscl_KF, &
       dhdt,kbcon,hc,he_cup,z_cup)

  do k=1,mkx
     do i=istart,iend
        if(ierr(i).eq.0)then
           dellat_ens(i,k,iedt) =  dellat(i,k)
           dellaq_ens(i,k,iedt) =  dellaq(i,k)
        else
           dellat_ens(i,k,iedt) = 0.
           dellaq_ens(i,k,iedt) = 0.
        endif
     enddo
  enddo
  !250  continue

  !------------------------------ FEEDBACK ------------------------------------
  call cup_output_ens_shal(xf_ens,ierr,dellat_ens,dellaq_ens, &
       outt,outq,xmb,ktop,mix,mgmxp,mjx,mgmyp,mkx,mgmzp,istart, &
       iend,j,'shallow',maxens2,maxens,ipr,jpr,iens, &
       maxens3,ensdim,xfac1,&!sgrell
       m1,m2,m3,dti)

  !---------------------------done shallow cumulus scheme -------------

  !-------Salva parametros nas analises DO RAMS
  do i=istart,iend
     xierr(i,j)=float(ierr(i))
     if(ierr(i).eq.0)then
          upmf(i,j)= xmb(i)! Shallow updraft mass flux
          xk22(i,j)= float(k22(i))
         xktop(i,j)= float(ktop(i))
        xkbcon(i,j)= float(kbcon(i))
     elseif(ierr(i).ne.0.and.ierr(i).ne.20)then
          upmf(i,j) = 0.
          xk22(i,j) = 0.
         xktop(i,j) = 0.
        xkbcon(i,j) = 0.
     endif
  enddo


  !-------Salva parametros DO Shallow CUP para uso no transporte convectivo:
  if(CCATT == 1) then
    do i=istart,iend
       ierr4d(i) =  ierr(i)
       jmin4d(i) =  1 ! no DOwndrafts
       kdet4d(i) =  1 ! no DOwndrafts
       k224d(i)  =  k22(i)
       kbcon4d(i)=  kbcon(i)
       ktop4d(i) =  ktop(i)
       kstabi4d(i)= kstabi(i)
       kstabm4d(i)= kstabm(i)
       !   kpbl4d(i) =  kpbl(i)
       kpbl4d(i) =  k22(i)
       xmb4d(i)  =  xmb(i)
       edt4d(i)  =  0.  ! no DOwndrafts

   			upmf(i,j) =  0.
       if(ierr(i).eq.0) upmf(i,j) = xmb(i)

       do k=1,mkx
   	  !srf- neste ponto iens=iens_tmp = 1
   	  if(iens.eq.1) then
   	     zcup5d(k,i) = zo_cup(i,k)
   	     pcup5d(k,i) = po_cup(i,k)
   	  endif
   	  enup5d(k,i) = mentr_rate2(i,k)
   	  endn5d(k,i) = 0.	    ! no downdrafts
   	  deup5d(k,i) =  cd(i,k)
   	  dedn5d(k,i) = 0.	    ! no downdrafts
   	  zup5d(k,i) = zuo(i,k)
   	  zdn5d(k,i) = 0.	   ! no downdrafts
   	  !?? it's save for use at wet-deposition scheme (in-cloud)
   	  prup5d(k,i) = 0.	    ! no precip for shallow
   	  clwup5d(k,i) = qrco(i,k)   ! cloud liq  water  - only for upfradt
   	  tup5d(k,i) = t_cup(i,k)
      !>>> em verdade deveria ser a temperatura da parcela
       enddo
    enddo
  endif
end subroutine CUP_enss_shal

!--------------------------------------------------------------------

subroutine cup_dellas_shallow(ierr,z_cup,p_cup,he,mix,mgmxp,mkx,mgmzp, &
     istart,iend,della,itest,j,zu,cd,hc,ktop,k22,kbcon,mentr_rate2, &
     he_cup,ipr,jpr,name)
  use Phys_const, only: g
  implicit none
  character *(*) name
  integer mix,mgmxp,mkx,mgmzp,i,k,istart,iend,itest,j
  real z_cup(mgmxp,mgmzp), p_cup(mgmxp,mgmzp),he(mgmxp,mgmzp)
  real della(mgmxp,mgmzp), hc(mgmxp,mgmzp), cd(mgmxp,mgmzp)
  real zu(mgmxp,mgmzp),he_cup(mgmxp,mgmzp)
  integer kbcon(mgmxp),ktop(mgmxp),k22(mgmxp),ierr(mgmxp),ipr,jpr
  real dp,dz,mentr_rate,subin,entup,detup,subDOwn,entupk,detupk,totmas
  real mentr_rate2(mgmxp,mgmzp)
  !real xsum,xsumt

  !PRINT-----------
  !      IF(j.eq.jpr)THEN
  !      IF(itest.eq.1)THEN
  !        PRINT *,'this one for h '
  !        PRINT *,'in dellas k22(i),kbcon(i),ktop(i)'
  !        PRINT *,k22(i),kbcon(i),ktop(i)
  !      ELSE
  !        PRINT *,'this one for q '
  !        PRINT *,'in dellas k22(i),kbcon(i),ktop(i)'
  !        PRINT *,k22(i),kbcon(i),ktop(i)
  !      ENDIF
  !      ENDIF
  !PRINT-----------

  do K=2,MKX
     do I=ISTART,IEND
        della(i,k)=0.
     enddo
  enddo
  !       xsum=0.
  !       xsumt=0.
  !
  do k=2,mkx-1
     do i=istart,iend
        if (ierr(i).ne.0) cycle
        if (k.gt.ktop(i)) cycle
        !
        !--- specify detrainment of downdraft, has to be consistent
        !--- with zd calculations in soundd.
        !
        dz    = Z_cup(i,k+1)-Z_cup(i,k)
        subin = zu(i,k+1)
        entup = 0.
        detup = 0.
        if(k.ge.kbcon(i).and.k.lt.ktop(i))then
           entup = mentr_rate2(i,k)*dz*zu(i,k)
           detup = cd(i,k+1) *dz*zu(i,k)
        endif
        subdown = zu(i,k)
        entupk  = 0.
        detupk  = 0.

        if(k.eq.k22(i)-1)then
           !srf-fev2003  entupk  = zu(i,kpbl(i))
           entupk  = zu(i,k22(i))
        endif

        if(k.eq.ktop(i)-0)then
           detupk  = zu(i,ktop(i))
           subin   = 0.
        endif

        if(k.lt.kbcon(i))then
           detup   = 0.
        endif
        !
        !--- Changed due to subsidence and entrainment
        !
        totmas=subin-subdown+detup-entup-entupk+detupk
        !         IF(j.eq.jpr.and.i.eq.ipr)PRINT *,'k,totmas,sui,sud = ',k,
        !     1   totmas,subin,subDOwn
        !         IF(j.eq.jpr.and.i.eq.ipr)PRINT *,'updr stuff = ',detup,
        !     1      entup,entupk,detupk
        if(abs(totmas).gt.1.e-6)then
           print *, '**totmass shallow********', i, j, k, totmas, name
           print *,k22(i),kbcon(i),ktop(i)
           print *,'updr stuff = ',subin,subDOwn,detup,entup,entupk,detupk
           stop
        endif

        dp =  100.*( p_cup(i,k)-p_cup(i,k+1) )
        della(i,k)=(subin*he_cup(i,k+1)-subdown*he_cup(i,k)+ &
             detup*.5*( HC(i,K+1)+ HC(i,K))-entup*he(i,k)-&
             entupk*he_cup(i,k22(i))+detupk*hc(i,ktop(i)))*g/dp

        !     IF(j.eq.3.and.i.eq.120)xsumt=xsumt+totmas
        !     IF(j.eq.3.and.i.eq.120)hesum=hesum+he(i,k)*dp
        !     IF(j.eq.3.and.i.eq.120)xsum=xsum+della(i,k)*dp
        !     IF(j.eq.3.and.i.eq.120)PRINT *,'xsum = ',xsum
        !      IF(i.eq.ipr.and.j.eq.jpr)THEN
        !        PRINT *,k,della(i,k),subin*he_cup(i,k+1),subDOwn*he_cup(i,k)
        !        PRINT *,k,detup*.5*(HC(i,K+1)+HC(i,K)),detupk*hc(i,ktop(i)),
        !     1         entup*he(i,k)
        !        PRINT *,k,he_cup(i,k+1),he_cup(i,k),entupk*he_cup(i,k)
        !      ENDIF
     enddo
  enddo
end subroutine cup_dellas_shallow

!------------------------------------------------------------

subroutine cup_forcing_ens_shal(aa0,aa1,xaa0,mbdt,dtime,xmb,ierr, &
       	   mix,mgmxp,mjx,mgmyp,mkx,mgmzp,istart,iend,xf,j,        &
     	   name,xland,maxens,iens,maxens2,ipr,jpr,maxens3,        &
     	   ensdim,p_cup,ktop,icoic,iedt,m2,m3, tscl_KF,&
	   dhdt,kbcon,hc,he_cup,z_cup)

  use Phys_const, only: g
  implicit none
  character *(*) name

  integer k,i,istart,iend,mix,mgmxp,mjx,mgmyp,mkx,mgmzp,j,m2,m3
  integer maxens,maxens3,ensdim,iens,nall,maxens2,ipr,jpr

!24062007 -not using
!------ ensemble 3 dimension = 10
!  integer kclim
!  integer,parameter :: mkxcrt=25
!  real pcrit(mkxcrt),acrit(mkxcrt),acritt(mkxcrt),aclim1
!  real aclim2,aclim3,aclim4
!  data pcrit/ 850., 837.5, 825., 812.5, 800., 787.5, &
!       775., 762.5, 750., 737.5, 725., 712.5, &
!       700., 687.5, 675., 662.5, 650., 637.5, &
!       625., 612.5, 600., 550. , 500., 450., &
!       400./
!  data acrit/ 6.323E-02,  5.795E-02, 5.390E-02,  5.236E-02, &
!       4.450E-02,  4.965E-02, 5.000E-02,  4.983E-02, &
!       5.530E-02,  5.289E-02, 6.080E-02,  5.883E-02, &
!       6.640E-02,  6.766E-02, 7.070E-02,  7.937E-02, &
!       7.500E-02,  9.396E-02, 0.108,      0.111, &
!       0.130,      0.152,     0.221,      0.315, &
!       0.368/
!  data acritt/0.203,  0.299, 0.359,  0.403, 0.515,  0.478, &
!       0.518,  0.530, 0.521,  0.565, 0.543,  0.588, &
!       0.566,  0.602, 0.596,  0.611, 0.625,  0.619, &
!       0.645,  0.627, 0.665,  0.659, 0.688,  0.743, &
!       0.813/

  integer ktop(mgmxp),iedt
  real p_cup(mgmxp,mgmzp),z_cup(mgmxp,mgmzp)
  integer ierr(mgmxp)
  real aa0(mgmxp),aa1(mgmxp),xaa0(mgmxp,maxens),xmb(mgmxp)
  real mbdt(maxens),dtime
  real     xf(mgmxp,mgmyp,ensdim),xland(mgmxp,mgmyp),dhdt(mgmxp,mgmzp)
  real xff_ens3(maxens3),xk(maxens),xff0
  integer nens,ne,n,icoic,  nens3
!
  real  tscl_KF
  real, parameter :: dec_fudge = 0.
!srf- for new closure bound layer quase-equi
  real blqe,dp,dz
  integer kbcon(mgmxp)
  real hc(mgmxp,mgmzp),he_cup(mgmxp,mgmzp)

  nens=0

  !--- LARGE SCALE FORCING
  !
  do I=ISTART,IEND

    xmb(i)=0.

     if(ierr(i).eq.0)then
!

! 24/06/2007
!-srf - not using this forcing ---
!        kclim=0
!        do k=mkxcrt,1,-1
!           if(p_cup(i,ktop(i)).lt.pcrit(k))then
!              kclim=k
!              GO TO 9
!           endif
!        enddo
!        if(p_cup(i,ktop(i)).gt.pcrit(1))kclim=1
!9       continue
!        k= max(kclim-1,1)
!        aclim1= (acrit(kclim) -dec_fudge*acrit(kclim) )*1.e3
!        aclim2= (acrit(k)     -dec_fudge*acrit(k)     )*1.e3
!        aclim3= (acritt(kclim)-dec_fudge*acritt(kclim))*1.e3
!        aclim4= (acritt(k)    -dec_fudge* acritt(k)   )*1.e3
!


        !
        !---- Grell's closure
        !
        xff0	   =  (AA1(I)-AA0(I))/dtime


!----------
!-srf out/2004: only for xff0 output/ shallow
!!srf-2009	sgrell2_2d(i,j) = xff0
!----------

        xff_ens3(1)=    xff0
        xff_ens3(2)= .9*xff_ens3(1)
        xff_ens3(3)=1.1*xff_ens3(1)

! 24/06/2007
!-srf - not using this forcing ---
!        !--- More original Arakawa-Schubert
!        xff_ens3(4)=max(0.,(AA1(I)-aclim1)/dtime)
!        xff_ens3(5)=max(0.,(AA1(I)-aclim2)/dtime)
!        xff_ens3(6)=max(0.,(AA1(I)-aclim3)/dtime)
!        xff_ens3(7)=max(0.,(AA1(I)-aclim4)/dtime)
!
        !- boundary layer quasi-equilibrium closure applied to kbcon level
	blqe = 0.
	!dz = 0. ! to check if the integral is done correctly
	do k=1,kbcon(i)-1 ! kconb is defined at Z_cup levels
         dp =  100.* ( p_cup(i,k) - p_cup(i,k+1) )
	 blqe = blqe + dhdt(i,k)*dp/g
	 !dz =dz+ z_cup(i,k+1)-z_cup(i,k)
	 !print*,k,z_cup(i,k),z_cup(i,k+1)-z_cup(i,k),dp/g,z_cup(i,kbcon(i))
	enddo

	xff_ens3(4)=max(0.,blqe/( hc(i,kbcon(i)) - he_cup(i,kbcon(i)) ))
	!print*,'blqe',xff_ens3(4),blqe,hc(i,kbcon(i)) , he_cup(i,kbcon(i))
	!print*,'dz=', dz

	xff_ens3(5)=0.75 * xff_ens3(4)
	xff_ens3(6)=1.25 * xff_ens3(4)
	xff_ens3(7)=       xff_ens3(4)


        !--- More like Fritsch Chappel or Kain Fritsch (plus triggers)
        xff_ens3(8) =AA1(I)/(60.*20.)
        xff_ens3(9) =AA1(I)/(60.*30.)
        xff_ens3(10)=AA1(I)/(60.*40.)
        !- A. Betts suggestion
        !xff_ens3(8) =AA1(I)/(tscl_KF	)! tscl_KF is already in seconds
        !xff_ens3(9) =AA1(I)/(tscl_KF*0.9)
        !xff_ens3(10)=AA1(I)/(tscl_KF*0.8)
           !

        do nens=1,maxens
           XK(nens)=(XAA0(I,nens)-AA1(I))/MBDT(nens)
           if(xk(nens).le.0.and.xk(nens).gt.-1.e-9) xk(nens)=-1.e-9
           if(xk(nens).gt.0.and.xk(nens).lt.+1.e-9) xk(nens)=+1.e-9
        enddo

        !--- Add up all ensembles
        do ne=1,maxens
           !
           !--- for every xk, we have maxens3 xffs
           !--- iens is from outermost ensemble (most expensive!
           !
           !--- iedt (maxens2 belongs to it)
           !--- is from second, next outermost, not so expensive
           !
           !--- so, for every outermost loop, we have maxens*maxens2*3
           !--- ensembles!!! nall would be 0, IF everything is on first
           !--- loop index, then ne would start counting, THEN iedt, THEN iens....
           !

           nall= (iens-1)*maxens3*maxens*maxens2 &
	        +(iedt-1)*maxens3*maxens         &
                +(ne  -1)*maxens3

           !
           !--- Special treatment for stability closures
           !

           if(xff0.gt.0  .and.  xk(ne).lt.0.)then
              xf(i,j,nall+1) =max(0.,-xff_ens3(1) /xk(ne))
              xf(i,j,nall+2) =max(0.,-xff_ens3(2) /xk(ne))
              xf(i,j,nall+3) =max(0.,-xff_ens3(3) /xk(ne))
           endif
           !
           if(XK(ne).lt.0.                  )then

! 24/06/2007
!-srf - not using this forcing --- (Arakawa&Schubert)
!              xf(i,j,nall+4) =max(0.,-xff_ens3(4) /xk(ne))
!              xf(i,j,nall+5) =max(0.,-xff_ens3(5) /xk(ne))
!              xf(i,j,nall+6) =max(0.,-xff_ens3(6) /xk(ne))
!              xf(i,j,nall+7) =max(0.,-xff_ens3(7) /xk(ne))

              xf(i,j,nall+8) =max(0.,-xff_ens3(8) /xk(ne))
              xf(i,j,nall+9) =max(0.,-xff_ens3(9) /xk(ne))
              xf(i,j,nall+10)=max(0.,-xff_ens3(10)/xk(ne))
           endif
 	   !
           !--- boundary layer quasi-equilibrium closure
	   !
	   xf(i,j,nall+4) =xff_ens3(4)
	   xf(i,j,nall+5) =xff_ens3(5)
	   xf(i,j,nall+6) =xff_ens3(6)
	   xf(i,j,nall+7) =xff_ens3(7)

           !----- 1d closure ensemble -------------
           if(icoic.ge.1)then
               do nens3=1,maxens3
                     xf(i,j,nall+nens3)=xf(i,j,nall+icoic)
               enddo
           endif

        end do

!srf-nov/2004
       ! exit
        cycle
!
     elseif(ierr(i).ne.20.and.ierr(i).ne.0)then
!
        do n=1,ensdim
           xf(i,j,n)=0.
        enddo
!
     endif

  end do

end subroutine cup_forcing_ens_shal

!--------------------------------------------------------------------

subroutine cup_output_ens_shal(xf,ierr,dellat,dellaq,outt,outq,xmb, &
           ktop,mix,mgmxp,mjx,mgmyp,mkx,mgmzp, &
           istart,iend,j,name,maxens2,maxens,ipr, &
           jpr,iens,maxens3,ensdim,xfac1,m1,m2,m3,dti)

  use cup_output_vars, only: &
       xmb_ave,xmb_std,xmb_ske,xmb_cur, &
        pr_ave, pr_std, pr_ske, pr_cur, &
	 x_ave,  x_std,  x_ske,  x_cur, &
	 x_ave_cap,                     &
	 x_ave_cap1,x_ave_cap2,x_ave_cap3,x_ave_cap4,x_ave_cap5, &
        cup_output_vars_alloc,          &
        alloc_cup_output_vars

  implicit none
  character (LEN=*) name

  integer mix,mjx,mkx,istart,iend,mgmxp,mgmyp,mgmzp,ensdim,i,k,j,n,ipr, &
          jpr,m1,m2,m3,maxens,maxens2
  integer ktop(mgmxp),ierr(mgmxp),ncount,iens,maxens3  !,nens3

  real outtes,ddtes,dtt,dtq
  real xmb(mgmxp),xfac1(mgmxp)
  real outt(mgmxp,mgmzp),outq(mgmxp,mgmzp)
  real dellat(mgmxp,mgmzp,maxens2),dellaq(mgmxp,mgmzp,maxens2)
  real xf(mgmxp,mgmyp,ensdim)

  real dti,qmemf,qmem1,qmem2,thresh

  !--- use of statistics properties
  integer, parameter :: i_use_stat_prop = 0
  !data  i_use_stat_prop/0/
  real, parameter :: tunning = 0.
  !data tunning /0.0/

  do K=1,MKX
     do I=ISTART,IEND
        outt(i,k) = 0.
        outq(i,k) = 0.
     enddo
  enddo
  do I=ISTART,IEND
     xmb(i)  =0.
     xfac1(i)=1.
  enddo

  !--- calculate ensemble average mass fluxes
  if(i_use_stat_prop == 1 ) then
     !      IF(i_use_stat_prop == 1 .and. icoic == 0 ) THEN
     if(.not. cup_output_vars_alloc)  &
          call alloc_cup_output_vars(mgmxp,maxens,maxens3)
     !
     !--- calculate ensemble average mass fluxes
     !
!srf-agosto2004- por enquanto ate resolver o problema da conveccao profunda
!inclua os parametros dti e x_ave_cap*
! inclua sgrell1_3d, sgrell2_3d (veja chamada da rotina 'cup_output_ens_catt'
!     call massflx_stats_catt(xf,mix,mgmxp,mjx,mgmyp,ensdim,maxens,maxens2, &
!          maxens3,xmb_ave,xmb_std,xmb_cur,xmb_ske,x_ave,                   &
!          x_std,x_ske,x_cur,x_ave_cap,j,istart,iend,ierr,                  &
!          m1,m2,m3,sgrell,3,&
!	   x_ave_cap1, &
!	   x_ave_cap2, &
!	   x_ave_cap3, &
!	   x_ave_cap4, &
!	   x_ave_cap5, &
!	   dti)

     do I=ISTART,IEND
        if(ierr(i).eq.0)then
           !-- tunning process
           xmb(i)=xmb_ave(i)-tunning*xmb_std(i)
           if(xmb(i).lt.0.) xmb(i)=.1*xmb_ave(i)
           if(xmb_ave(i) .lt. 1.e-10) ierr(i) = 11
        endif
        xfac1(i)=xmb(i)
     enddo

  else
     !----  Simple average
     do I=ISTART,IEND
        ncount=0
        xmb(i)=0.
        if(ierr(i).eq.0)then
           do n=(iens-1)*maxens2*maxens*maxens3+1, &
                    iens*maxens2*maxens*maxens3
              if(xf(i,j,n).gt.0.)then
                 xmb(i) = xmb(i) + xf(i,j,n)
                 ncount = ncount + 1
              endif
           enddo
           if(ncount.gt.0)then
              xmb(i)=xmb(i)/float(ncount)
!srf-20/11/2005
! limita fluxo de massa ao maximo de 0.1 kgm^2/s^2
              xmb(i)=min(xmb(i),0.1)
!
           else
              xmb(i)=0.
              ierr(i)=13
           endif
!srf - print------------------
!	   if(j==7 .and. i==48) then
!	     !if( k==1)
!	     print*,'BRAMS: XMB'
!	     print*,xmb(i)
!	     print*,xf(i,j,:)
!	   endif
!srf - print------------------


        endif
     enddo
  endif

  !--- now do feedback
!  go to 400
!  ddtes=250.
!  if(name.eq.'shallow')ddtes=500.
!
!  do K=1,MKX
!     do I=ISTART,IEND
!	 dtt=0.
!	 dtq=0.
!	 if(ierr(i).eq.0.and.k.le.ktop(i))then
!	    do n=1,maxens2
!	       dtt = dtt + dellat(i,k,n)
!	       dtq = dtq + dellaq(i,k,n)
!	    enddo
!	    outtes=dtt*xmb(i)*86400./float(maxens2)
!	    if(outtes .gt. 2.*ddtes .and. k.gt.2) then
!	       xmb(i)= 2.*ddtes/outtes * xmb(i)
!	       outtes=    ddtes
!	    endif
!	    if(outtes .lt. -ddtes) then
!	       xmb(i)= -ddtes/outtes * xmb(i)
!	       outtes= -ddtes
!	    endif
!	    if(outtes .gt. .5*ddtes.and.k.le.2) then
!	       xmb(i)=   ddtes/outtes * xmb(i)
!	       outtes=.5*ddtes
!	    endif
!	    !if(i.eq.ipr.and.j.eq.jpr)print *,'xmb = ',i,j,k,xmb(i)
!	    outt(i,k)= xmb(i) *dtt / float(maxens2)
!	    outq(i,k)= xmb(i) *dtq / float(maxens2)
!
!	 endif
!     enddo
!  enddo
!
!-srf- reformulated heating/cooling limits according with G3d version
!
!-srf  Thresh=200. ! max heating/cooling rate allowed  K/day
  Thresh=200. ! max heating/cooling rate allowed  K/day
  do I=ISTART,IEND
     if(ierr(i) /= 0)cycle
     do K=1,ktop(i)
	dtt=0.
	dtq=0.
	do n=1,maxens2
	   dtt = dtt + dellat(i,k,n)
	   dtq = dtq + dellaq(i,k,n)
	enddo
	outt(i,k)= xmb(i) *dtt / float(maxens2)
	outq(i,k)= xmb(i) *dtq / float(maxens2)
     enddo

     Qmemf=1.
     Qmem1=0.

     do K=1,ktop(i)
	Qmem1=abs(outt(i,k))*86400.

	 if(Qmem1 > Thresh)then
	   Qmem2 = Thresh/Qmem1
	   Qmemf = min(Qmemf,Qmem2)
	   Qmemf = max(0.,Qmemf)
	 endif
     enddo
      ! - limiting heating/cooling
     do K=1,ktop(i)
       outt(i,k) =  outt(i,k) *Qmemf
       outq(i,k) =  outq(i,k) *Qmemf
      enddo ! end k
      xmb(i)=xmb(i)*qmemf
  enddo
!-srf
!
  do i=istart,iend
     if(ierr(i).eq.0)then
        xfac1(i) = xmb(i)/(xfac1(i)+1.e-8)
     endif
  enddo

end subroutine cup_output_ens_shal
!---------------------------------------------------------------------------------------------
