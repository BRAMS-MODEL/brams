!---------------------------GRELL CUMULUS SCHEME---------------------------
subroutine CUPARTH_CATT(CCATT,  &
     mynum   	   , &
     mgmxp   	   , &
     mgmyp   	   , &
     mgmzp   	   , &
     m1      	   , &
     m2      	   , &
     m3      	   , &
     ia      	   , &
     iz      	   , &
     ja      	   , &
     jz      	   , &
     i0      	   , &
     j0      	   , &
     maxiens       , &
     iens          , &
     ngrid         , &
     ngrids_cp     , &
     dtime         , &
     time          , &
     ua            , &
     va            , &
     wa            , &
     theta         , &
     thetail       , &
     pp            , &
     pi0   	   , &
     dn0   	   , &
     rv    	   , &
     tke   	   , &
     tkmin 	   , &
     rcp           , &
     rrp           , &
     rpp           , &
     rsp           , &
     rap           , &
     rgp           , &
     rhp           , &
     topt          , &
     rtgt          , &
     tht           , &
     rtt           , &
     pt            , &
     outtem        , &
     outrt         , &
     outcl         , &
     precip        , &
!     sgrell1_3d    , &
!     sgrell2_3d    , &
!     sgrell9_3d    , &
!     sgrell1_2d    , &
     ierr4d        , &
     jmin4d        , &
     !
     kdet4d        , &
     k224d         , &
     kbcon4d       , &
     ktop4d        , &
     kpbl4d        , &
     kstabi4d      , &
     kstabm4d      , &
     xmb4d 	   , &
     edt4d 	   , &
     zcup5d 	   , &
     pcup5d 	   , &
     enup5d    	   , &
     endn5d 	   , &
     deup5d 	   , &
     dedn5d 	   , &
     zup5d 	   , &
     zdn5d 	   , &
     prup5d 	   , &
     clwup5d 	   , &
     tup5d 	   , &
     upmf   	   , &
     dnmf   	   , &
     xierr 	   , &
     xktop 	   , &
     xkbcon 	   , &
     xk22 	   , &
     xjmin 	   , &
     xkdt  	   , &
     xiact_p 	   , &
     xiact_c  	   , &
     confrq   	   , &
     frqanl   	   , &
     deltaxy2 	   , &
     patch_area    , &
     npat          , &
     level         , &
     glat          , &
     glon          , &
     sflux_r       , &
     sflux_t       , &
     trigg 	   , &
     autoconv        &
      )


  !  Modules for Grell Parameterization
  use mem_grell_param, only : maxens,  & !INTENT(IN)
       maxens2,                        & !INTENT(IN)
       maxens3,                        & !INTENT(IN)
       ensdim,                         & !INTENT(IN)
       icoic                             !INTENT(IN)

 use mem_scratch2_grell, only : &
       massflx,                  &       !INTENT(OUT)
       iact_old_gr,              &       !INTENT(OUT)
       aa0,                      &       !INTENT(OUT)
       xland,                    &       !INTENT(OUT)
       kdt,                      &       !INTENT(OUT)
       iact_gr,                  &       !INTENT(OUT)
       kdet,                     &       !INTENT(OUT)
       pret,                     &       !INTENT(OUT)
       mconv,                    &       !INTENT(OUT)
       umean,                    &       !INTENT(OUT)
       vmean,                    &       !INTENT(OUT)
       pmean,                    &       !INTENT(OUT)
       TER11,                    &       !INTENT(OUT)
       glatg,                    &       !INTENT(OUT)
       glong,                    &       !INTENT(OUT)
       PSUR,                     &       !INTENT(OUT)
       PO,                       &       !INTENT(OUT)
       US_Grell,                 &       !INTENT(OUT)
       VS_Grell,                 &       !INTENT(OUT)
       OMEG,                     &       !INTENT(OUT)
       T,                        &       !INTENT(OUT)
       Q,                        &       !INTENT(OUT)
       TKEG,                     &       !INTENT(OUT)
       RCPG,                     &       !INTENT(OUT)
       TN,                       &       !INTENT(OUT)
       QO,                       &       !INTENT(OUT)
       P,                        &       !INTENT(OUT)
       OUTT,                     &       !INTENT(OUT)
       OUTQ,                     &       !INTENT(OUT)
       OUTQC,                    &       !INTENT(OUT)
       DIRECTION,                &       !INTENT(OUT)
       massfln                           !INTENT(?)

  use Phys_const, only: cp, p00, tcrit, g, cpor !INTENT(IN)

  !use extras            , only: extra3d,extra2d,na_EXTRA3D
  use mem_varinit, only: nudlat
  use mem_grid, only: nxtnest, initial, jdim
  use node_mod, only: ibcon

  !- incluindo o efeito de aerosois na precipita��o
  use mem_carma, only: carma


  !  this is the lowest level, then comes ensemble 1!
  !       ensdim=maxiens*maxens*maxens2*maxens3 !Ensemble dimension
  !
  !srf- ICOIC is used for choice a specific closure
  ! icoic = 0 -> ensemble (all closures)  [en]
  ! icoic = 1 -> Grell                    [gr]
  ! icoic = 4 -> low level omega          [lo]
  ! icoic = 7 -> moisture convergence     [mc]
  ! icoic =10 -> like Fritsch Chappel or Kain Fritsch [sc]
  ! icoic =13 -> Arakawa-Schubert         [as]

  implicit none
  integer, intent(in) ::CCATT, mgmxp, mgmyp, mgmzp, ngrid, ngrids_cp
  integer, intent(in) :: maxiens ! maxens,maxens2,maxens3,ensdim
  integer, intent(in) :: iens,trigg ,autoconv

  integer, intent(in) :: m1, m2, m3, ia, iz, ja, jz, i0, j0, mynum, npat, level

  integer :: j1, j2 !Local

  real, intent(in) :: time

  real, dimension(m1,m2,m3), intent(in) :: ua, va, wa, theta, pp, pi0, dn0,   &
                               rv, tht, rtt, pt,   &
			       tke,rcp,rrp,rpp,rsp,rap,rgp,rhp,thetail

  real, dimension(m1,m2,m3), intent(inout) :: outtem, outrt,outcl

  !- For use without CATT
  !real, dimension(m1,m2,m3) :: sgrell1_3d,sgrell2_3d,sgrell9_3d
  !real, dimension(m2,m3)    :: sgrell1_2d
  !real, pointer :: sgrell1_3d(:,:,:), sgrell2_3d(:,:,:), sgrell3_3d(:,:,:)
  !real, pointer :: sgrell1_2d(:,:)


  real, dimension(m2,m3,npat), intent(in)  :: patch_area

  real, intent(in) :: tkmin, confrq, frqanl, deltaxy2


  real, dimension(m2,m3), intent(in)    :: topt,glat,glon,sflux_r,sflux_t
  real, dimension(m2,m3), intent(inout) :: PRECIP
  real, dimension(m2,m3), intent(in)    :: rtgt
  !real, dimension(m2,m3)     :: train !Not used

  real, dimension(m2,m3) :: upmf, dnmf, xierr, xktop, xkbcon, xjmin,  &
                            xkdt, xiact_p, xiact_c,xk22

  !-------salva parametros da CUP para uso no transporte convectivo:
  integer, dimension(mgmxp,mgmyp,maxiens,ngrids_cp) :: ierr4d, jmin4d,  &
       kdet4d, k224d, kbcon4d, ktop4d, kpbl4d, kstabi4d, kstabm4d

  real,dimension(mgmxp,mgmyp,maxiens,ngrids_cp) :: xmb4d, edt4d

  real,dimension(mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp) :: enup5d,  &
                                                         endn5d,  &
							 deup5d,  &
							 dedn5d,  &
							 zup5d,   &
							 zdn5d,   &
            						 prup5d,  &
            						 clwup5d, &
            						 tup5d,   &
            						 zcup5d,  &
            						 pcup5d
  ! no futuro, coloque estas var no mem_scratch_grell ....
  real, dimension(mgmxp) :: sflux_rg, sflux_tg

  real, dimension(mgmxp) :: dnsup  ! densidade em superficie, para converter o
                                   ! fluxo de calor sensivel

  real, dimension(mgmxp) :: ccn    ! calculado atraves da AOT que vem do Carma


  real, dimension(mgmxp,mgmzp) :: mentr_rate2, wu, dn0_2,xcoast

  integer :: istart, iend, i, j, k, mix, mjx, mkx, kr, m, ipr, jpr !kk
  real    :: vspeed, dp, dtime, dq, cpdTdt, exner
  real    :: r_liq, r_sol,tempk,fxc,dfxcdt

  real    :: aot500 ! variavel local para a aot
  real :: dti, tscl_KF !Local


  !----------------------------------------------------------------------
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

!----- tmp
!istart=19; iend=19; j1=27;j2=27!----- tmp
!----- tmp

  !---- Coordenadas para escrita ascii
  !      ipr=19 - i0
  !      jpr=19 - j0
  ipr=0 !- i0
  jpr=0 !- j0

!- for ensemble average
  if (frqanl /= 0.) then
      dti = confrq/frqanl
  else
      dti = 0.
  endif

  do j=1,m3  ! loop em todo dominio para passar informacoes da fronteira
     do i=1,m2 ! dos nodes
        massflx(i,j)     = dnmf(i,j)
        iact_old_gr(i,j) = int(xiact_p(i,j))
     enddo
  enddo

  !-not used anymore
  ! A. Betts suggestion for the KF timescale < DELTAX  / 50 m/s
  !  tscl_KF =  sqrt(deltaxy2) / 50.  !units: sec
                                     ! use 25 m/s for shallow
  !- loop externo : j
  do j=j1,j2


     do I = ISTART,IEND
        aa0(i)           =0.
        xland(i,j)       = patch_area(i,j,1) !flag < 1 para land
	                                     !flag  =1 para water
        !if(xland(i,j) < 0.8) xland(i,j) = 0.
	!sgrell1_2d(i,j)=xland(i,j)
        !xcoast(i,j) = 0.
	!if(xland(i,j) >= 0.8) then
	!  if(0.25*(xland(i+1,j)+xland(i-1,j) + xland(i,j+1)+xland(i,j-1)) &
	!                .lt.0.7)  xcoast(i,j) = 1.
	!endif
	!sgrell1_2d(i,j)=xcoast(i,j)

	iact_gr(i,j)     = 0
        kdt(i,j)         = 0
        precip(i,j)      = 0.
     enddo


     !--- Prepare input, erase output

     do I = ISTART,IEND
        kdet(i)  =2
        pret(i)  =0.
        mconv(i) =0.
        umean(i) =0.
        vmean(i) =0.
        pmean(i) =0.
     enddo


    if(autoconv == 2) then
          do I = ISTART,IEND

            !aot500 = .1 + MAX(carma(ngrid)%aot(i,j,11),0.)
            aot500 = .01! + MAX(carma(ngrid)%aot(i,j,11),0.)
            !ccn(i) = max( 300., ( 370.37*aot500 )**1.555 )
            ccn(i) = max( 150., ( 370.37*aot500 )**1.555 )

          enddo
     endif
     !------- Transfere valores do RAMS para o esquema de cumuls
     do K=1,MKX
        kr = K + 1          ! nivel K da grade do Grell corresponde ao
                            ! nivel K + 1 do RAMS
        do I = ISTART,IEND

           TER11(I)= topt(i,j)
           glatg(i) = glat(i,j)
           glong(i) = glon(i,j)

           sflux_rg(i) = sflux_r(i,j)
	   sflux_tg(i) = sflux_t(i,j)
	   dnsup(i)    = (dn0(1,i,j)+dn0(2,i,j))*0.5

	   dn0_2(i,k) = 0.001*dn0(kr,i,j) ! unidades: g/cm�

           ! Pressure in mbar
           PSUR(I) = .5*( ((pp(1,i,j)+pi0(1,i,j))/cp)**cpor*p00 +  &
                          ((pp(2,i,j)+pi0(2,i,j))/cp)**cpor*p00 )*1.e-2
           PO(I,K) = ((pp(kr,i,j)+pi0(kr,i,j))/cp)**cpor*p00*1.e-2

           US_Grell(I,K) = .5*( ua(kr,i,j) + ua(kr,i-1,j) )
           VS_Grell(I,K) = .5*( va(kr,i,j) + va(kr,i,j-1) )
           OMEG(I,K)   = -g*dn0(kr,i,j)*.5*( wa(kr,i,j)+wa(kr-1,i,j) )

           T(I,K)  = theta(kr,i,j)*(pp(kr,i,j)+pi0(kr,i,j))/cp
           Q(I,K)  = rv(kr,i,j)

           !- variables for PBL top height
	   TKEG(I,K) = TKE(kr,i,j)
	   RCPG(I,K) = RCP(kr,i,j)
           !        Calcula tendencia projetada na temperatura em funcao
           !        das tendencias de theta e PI : cp*T=Pi*Theta
           exner= pp(kr,i,j)+pi0(kr,i,j)
           cpdTdt= exner*tht(kr,i,j) + theta(kr,i,j)*pt(kr,i,j)
           !cpdTdt  = exner*tht(kr,i,j) ! assumindo PT(KR,I,J) << exner*THT(KR,I,J)/theta

           !        Temperatura projetada se a conveccao nao ocorrer
           TN(I,K) = T(I,K) + ( cpdTdt/cp )*dtime

           !        Umidade projetada se a conveccao nao ocorrer
           QO(I,K) = Q(I,K) +   rtt(kr,i,j)*dtime

           P(I,K)  = PO(I,K)
!print*,'1gd',i,j,k,US_grell(I,K), VS_grell(I,K) , T(I,K),q(I,K),tn(i,k),qo(i,k),po(i,k)

           if((PSUR(I)-P(I,K)).gt.150.and.P(I,K).gt.300.)then
              DP       = -.5*(P(I,K+1)-P(I,K-1))
              UMEAN(I) = UMEAN(I)+US_Grell(I,K)*DP
              VMEAN(I) = VMEAN(I)+VS_Grell(I,K)*DP
              PMEAN(I) = PMEAN(I)+DP
           endif

           if(TN(I,K).lt.200.)    TN(I,K) = T(I,K)
           if(QO(I,K).lt.1.E-08)  QO(I,K) = 1.E-08

           OUTT(I,K)  = 0. !- Tendencia no campo de temperatura
	                   !  associada aos cumulus
           OUTQ(I,K)  = 0. !- Tendencia na razao de mist. de vapor d'agua
	                   !  assoc. aos cumulus
           OUTQC(I,K) = 0. !- Tendencia na razao de mistura de agua de nuvem e/ou gelo
                           !  associada aos cumulus
        enddo
     enddo

     do I = ISTART,IEND
        UMEAN(I)=UMEAN(I)/PMEAN(I)
        VMEAN(I)=VMEAN(I)/PMEAN(I)
        VSPEED=sqrt(UMEAN(I)*UMEAN(I)+VMEAN(I)*VMEAN(I))
        DIRECTION(I)=(atan2(UMEAN(I),VMEAN(I))+3.1415926)*57.29578
        if(DIRECTION(I).gt.360.)DIRECTION(I)=DIRECTION(I)-360.
        if(VSPEED.lt.5.)DIRECTION(I)=9999.
        !sgrell1_3d(24,i,j)=UMEAN(I)
        !sgrell1_3d(25,i,j)=VMEAN(I)
     enddo

     do K=2,MKX-1
        do I = ISTART,IEND
           dq=.5*(q(i,k+1)-q(i,k-1))
           !- convergencia de umidade da coluna (omega em pa/s)
           mconv(i)=mconv(i)+omeg(i,k)*dq/g
        enddo
     enddo
     do I = ISTART,IEND
        if(mconv(i) .lt. 0.)  mconv(i) = 0.
     enddo

     !---  CUMULUS PARAMETERIZATION
     !srf- aqui se deve colocal o loop no ensemble dependente do tipo de cumulus
     ! iens =1
     call CUP_enss_catt(ccatt,ngrid,mynum,m1,m2,m3,i0,j0,ipr,jpr,             &
          mgmxp,mgmyp,mgmzp,maxiens,maxens,maxens2,maxens3,             &
	  ensdim,icoic,j,iens,istart,iend,mix,mjx,mkx,                  &
	  massfln,massflx,iact_gr,iact_old_gr,xland,ter11,              &
	  aa0, t, q, tn, qo, po, pret, p, outt, outq, outqc,            &
          dtime, psur, us_grell, vs_grell,kdet,                         &
	  tcrit,time,mconv, omeg, direction, tkeg,rcpg,tkmin,           &
            ierr4d(1,j,iens,ngrid),   jmin4d(1,j,iens,ngrid),           &
            kdet4d(1,j,iens,ngrid),    k224d(1,j,iens,ngrid),           &
           kbcon4d(1,j,iens,ngrid),   ktop4d(1,j,iens,ngrid),           &
            kpbl4d(1,j,iens,ngrid),                                     &
          kstabi4d(1,j,iens,ngrid), kstabm4d(1,j,iens,ngrid),           &
             xmb4d(1,j,iens,ngrid),    edt4d(1,j,iens,ngrid),           &
          zcup5d(1,1,j,iens,ngrid), pcup5d(1,1,j,iens,ngrid),           &
          enup5d(1,1,j,iens,ngrid), endn5d(1,1,j,iens,ngrid),           &
          deup5d(1,1,j,iens,ngrid), dedn5d(1,1,j,iens,ngrid),           &
           zup5d(1,1,j,iens,ngrid),  zdn5d(1,1,j,iens,ngrid),           &
          prup5d(1,1,j,iens,ngrid),clwup5d(1,1,j,iens,ngrid),           &
           tup5d(1,1,j,iens,ngrid),                                     &
          upmf,dnmf,xierr,xktop,xkbcon,xk22,xjmin,xkdt,xiact_p,xiact_c, &
!          sgrell1_3d,&
!	  sgrell2_3d,&
!	  sgrell9_3d,&
!	  sgrell1_2d,&
          dti,tscl_KF      ,&
          glatg,glong, sflux_rg,sflux_tg, dnsup, dn0_2,ccn,trigg,autoconv,xcoast)

     !--- Output
     !srf out/2004: coupling with RAMS microphysics
    if(level <= 2) then
      do K=1,MKX-1
        kr = K + 1
         do I = ISTART,IEND
           ! Converte tend da temperatura (OUTT) em tend de theta (OUTTEM)
           ! cp*T=Pi*Theta => cp dT/dt = Theta*dPi/dt + Pi*dTheta/dt,
           ! Exner's function = pp(kr,i,j)+pi0(kr,i,j)
           exner          = pp(kr,i,j) + pi0(kr,i,j)
           ! tendencia do theta devida a conv profunda
	   !-srf 26062012 - corrigindo taxas de convers�o
           if(OUTT(I,K) /= 0.) &
           outtem    (kr,i,j) = cp/exner * OUTT(I,K) - theta(kr,i,j)*pt(kr,i,j)/exner

           ! tendencia do Rtotal  devida a conv profunda
           outrt(kr,i,j)  = OUTQ(I,K)  + OUTQC(I,K)
	   outcl(kr,i,j)  = 0.0
         enddo
       enddo
    elseif(level > 2) then
      do K=1,MKX-1
        kr = K + 1
         do I = ISTART,IEND
           ! converte tend da temperatura (outt) em tend de theta (outtem)
           ! cp*T=Pi*Theta => cp dT/dt = Theta*dPi/dt + Pi*dTheta/dt,
           ! Exner's function = pp(kr,i,j)+pi0(kr,i,j)
	   !-srf 26062012 - corrigindo taxas de convers�o
           if(OUTT(I,K) /= 0.) then!-srf 26062012

	     	exner= pp(kr,i,j) + pi0(kr,i,j)
             	! tendencia do theta  devida a conv profunda
             	outtem (kr,i,j) = cp/exner * outt(i,k) - theta(kr,i,j)*pt(kr,i,j)/exner

             	! tendencia do theta_il devida a conv profunda
             	r_liq= max(0.,rcp(kr,i,j) + rrp(kr,i,j))

             	r_sol= max(0.,rsp(kr,i,j)+rpp(kr,i,j)+      &
             		     rap(kr,i,j)+rgp(kr,i,j)+  &
             		     rhp(kr,i,j))

	     	tempk = theta(kr,i,j)*(exner)/cp ! air temp (Kelvin)

	     	if(tempk.le.253) then
	     	 fxc =   (2.5e6*r_liq+2.83e6*r_sol)/(cp*amax1(tempk,253.))
	     	 dfxcdt = 2.83e6*OUTQC(I,K)/(cp*amax1(tempk,253.))
             	 outtem (kr,i,j) = (1./(1+fxc))*( outtem (kr,i,j) - thetail(kr,i,j)*dfxcdt )

	     	else
	     	 fxc =   (2.5e6*r_liq+2.83e6*r_sol)/(cp*amax1(tempk,253.))
	     	 dfxcdt = 2.5e6*OUTQC(I,K)/(cp*amax1(tempk,253.)) - &
	     		  fxc/(cp*amax1(tempk,253.)) * cp * OUTT(I,K)

             	 outtem (kr,i,j) = (1./(1+fxc))*( outtem (kr,i,j) - thetail(kr,i,j)*dfxcdt )

	     	endif

	   endif!-srf 26062012

           ! tendencia da vapor d'agua devida a conv profunda
           outrt     (kr,i,j) = OUTQ(I,K) !+ OUTQC(I,K)
           ! tendencia da agua condensada devida a conv profunda
	   outcl(kr,i,j)= OUTQC(I,K)
         enddo
       enddo
      endif

    do I = ISTART,IEND
        PRECIP(I,J)=PRET(I)
!	 if(precip(i,j).gt.0.)	 print*,'pretgd =', i,j,pret(i),outtem(:,i,j); call flush(6)

    enddo

     do I = ISTART,IEND
        if(precip(i,j).le.0.)then
           iact_gr(i,j) =0
           precip(i,j)  =0.
           do k=1,mkx
              kr = k + 1
              outtem(kr,i,j)    = 0.! comente p/ testes com conservacao (c0=0.)
              outrt(kr,i,j)     = 0.! comente p/ testes com conservacao (c0=0.)
	      outcl(kr,i,j)= 0.! comente p/ testes com conservacao (c0=0.)
           enddo
           do k=1,ensdim
              massfln(i,j,k)=0.
           enddo
        else
           iact_gr(i,j)=1
        endif
     enddo

     !--- Salva nas analises

     do I = ISTART,IEND
       !massflx(i,j) = dnmf(i,j)
        xiact_c(i,j) = float(IACT_GR(I,J))
       !xiact_p(i,j) = float(iact_old_gr(I,J))
     enddo

  enddo     ! loop externo - j -

end subroutine CUPARTH_CATT

!--------------------------------------------------------------------

subroutine CUP_enss_catt(ccatt,ngrid, mynum, m1, m2, m3, i0, j0, ipr, jpr,            &
     	       mgmxp, mgmyp, mgmzp, maxiens, maxens, maxens2, maxens3, ensdim,  &
     	       icoic, j, iens, ISTART, IEND, mix, mjx, mkx, massfln, massflx,	&
     	       iact_gr, iact_old_gr, xland, Z1,                          &
     	       AAEQ, T, Q, TN, QO, PO, PRE, P, OUTT, OUTQ, OUTQC, DTIME, &
     	       PSUR, US, VS, KDET, TCRIT, time, mconv, omeg, direction,  &
	       tkeg,rcpg,tkmin,  &
     	       ierr4d, jmin4d, kdet4d, k224d, kbcon4d,ktop4d, kpbl4d,    &
     	       kstabi4d, kstabm4d, xmb4d, edt4d,                         &
	       zcup5d,pcup5d, enup5d, endn5d, deup5d, dedn5d,  &
     	       zup5d, zdn5d, prup5d,clwup5d,tup5d,             &
     	       upmf, dnmf, xierr, xktop, xkbcon, xk22, xjmin,xkdt, xiact_p, xiact_c, &
!               sgrell1_3d,&
!	       sgrell2_3d,&
!	       sgrell9_3d,&
!	       sgrell1_2d,&
	       dti,tscl_KF,&
               glatg,glong,sflux_rg,sflux_tg, dnsup, dn0_2,ccn,trigg,autoconv,xcoast)

!- USE Modules for Grell Cumulus Parameterization
  use mem_scratch3_grell
  use Phys_const, only: rgas,cp,rm,p00,g,cpor,pkdcut
  use kbcon_ecmwf
  !tmp
  !use extras            , only: extra3d,extra2d,na_EXTRA3D
  !USE mod_therm_lib, ONLY:  esat               ! Function

  implicit none
  integer ccatt,maxiens,maxens,maxens2,maxens3,ensdim
  integer mix,mjx,mkx,mgmxp, mgmyp, mgmzp
  integer nall,nens,iens,iedt,trigg,autoconv        !,ktau
  integer nens3
  integer icoic
  real time

  !--- Input variables -----------------------------
  !
  ! basic environmental input includes moisture convergence (mconv)
  ! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
  ! convection for this call only and at that particular gridpoint
  !

  real mconv(mgmxp), Z1(mgmxp), direction(mgmxp), AAEQ(mgmxp),  &
  pre(mgmxp), PSUR(mgmxp), glatg(mgmxp),glong(mgmxp), sflux_rg(mgmxp),sflux_tg(mgmxp), dnsup(mgmxp)

  real T(mgmxp,mgmzp), Q(mgmxp,mgmzp), TN(mgmxp,mgmzp), QO(mgmxp,mgmzp),   &
       P(mgmxp,mgmzp), PO(mgmxp,mgmzp), US(mgmxp,mgmzp), VS(mgmxp,mgmzp),  &
       omeg(mgmxp,mgmzp), tkeg(mgmxp,mgmzp),rcpg(mgmxp,mgmzp)

  real mentr_rate2(mgmxp,mgmzp)

  real massfln(mgmxp,mgmyp,ensdim)
  real massflx(mgmxp,mgmyp), xland(mgmxp,mgmyp)
  integer iact_gr(mgmxp,mgmyp), iact_old_gr(mgmxp,mgmyp), kdet(mgmxp)
  integer ngrid, mynum, i0, j0, ipr, jpr, m1, m2, m3, kr

  !-------CU_P parameters for the convective transport:
  integer, dimension(mgmxp) :: ierr4d, jmin4d, kdet4d, k224d, kbcon4d,  &
       ktop4d, kpbl4d, kstabi4d, kstabm4d

  real, dimension(mgmxp) :: xmb4d, edt4d

  real, dimension(mgmzp,mgmxp) ::            &
                    zcup5d, pcup5d,          &
                    enup5d, endn5d, deup5d,  &
                    dedn5d,  zup5d,  zdn5d,  &
                    prup5d,clwup5d,   tup5d

  real, dimension(mgmxp,mgmzp) :: wu

  !------Variables saved in RAMS Analisys
  ! use (m1,m2,m3) para dimensionar os vetores que sao
  ! escritos nas analises do RAMS
  real, dimension(m2,m3) :: upmf, dnmf, xierr, xktop, xkbcon, xjmin,    &
       xkdt, xiact_p, xiact_c, xk22
  !real, dimension(m1,m2,m3) :: sgrell1_3d,sgrell2_3d,sgrell9_3d
  !real, dimension(m2,m3)    :: train !sgrell1_2d
  real tkmin,dti,tscl_KF
  !
  !--- Work variables - Allocatable in this point --
  !
  integer fquasi, fstab, fmconv, iresult
  integer ki, m
  integer I, J, K, ISTART, IEND
  real mbdt

  !--- Output variables ----------------------------
  ! outt   = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip

  real OUTQC(mgmxp,mgmzp), OUTT(mgmxp,mgmzp), OUTQ(mgmxp,mgmzp)

  real day, dz, tcrit, dtime
  real dellaqsum, dellaqcsum, dp

  !--- New entrainment/detrainment related stuff --------------------

  real mentr_rate, mentrd_rate, entr_rate, radius,              &  !entrd_rate,  &
       massfld, zcutdown, edtmax, edtmin, depth_min, zkbmax,    &
       z_detr, zktop, dh, cap_maxs,masscon

  real ::  fdiur, tke_start,tkeminnot,qs1,qsk
  integer :: kstart,kstart_way
  integer,parameter :: i_cup_dir_flag = 0

  real, dimension (m1,m2,m3) :: dn0

  real, dimension(mgmxp,mgmzp) :: dn0_2,xcoast
  real, dimension(mgmxp):: ccn,closure_n

  !---   the following are your basic environmental
  !	 variables. They carry a "_cup" if they are
  !	 on model cloud levels (staggered). They carry
  !	 an "o"-ending (z becomes zo), if they are the forced
  !	 variables. They are preceded by x (z becomes xz)
  !	 to indicate modification by some typ of cloud
  !
  ! z  	        = heights of model levels
  ! q  	        = environmental mixing ratio
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
  ! he_cup      = moist static energy on model cloud levels
  ! hes_cup     = saturation moist static energy on model cloud levels
  ! gamma_cup   = gamma on model cloud levels
  ! hcd         = moist static energy in downdraft
  ! zd          = normalized downdraft mass flux
  ! dby         = buoancy term
  ! entr        = entrainment rate
  ! zd          = downdraft normalized mass flux
  ! entr        = entrainment rate
  ! hcd         = h in model cloud
  ! bu          = buoancy term
  ! zd          = normalized downdraft mass flux
  ! gamma_cup   = gamma on model cloud levels
  ! mentr_rate  = entrainment rate
  ! qcd         = cloud q (including liquid water) after entrainment
  ! qrch        = saturation q in cloud
  ! pwd         = evaporate at that level
  ! pwev        = total normalized integrated evaoprate (I2)
  ! entr        = entrainment rate
  ! z1          = terrain elevation
  ! entr        = downdraft entrainment rate
  ! jmin        = downdraft originating level
  ! kdet        = level above ground where downdraft start detraining
  ! psur        = surface pressure
  ! z1          = terrain elevation
  ! pr_ens      = precipitation ensemble
  ! xf          = mass flux ensembles
  ! massfln     = downdraft mass flux ensembles used in next timestep
  ! omeg        = omega from large scale model
  ! mconv       = moisture convergence from large scale model
  ! zd          = downdraft normalized mass flux
  ! zu          = updraft normalized mass flux
  ! dir         = "storm motion"
  ! mbdt        = arbitrary numerical parameter
  ! dtime       = dt over which forcing is applied
  ! iact_gr_old = flag to tell where convection was active
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! dby         = buoancy term
  ! ktop        = cloud top (output)
  ! xmb         = total base mass flux
  ! hc          = cloud moist static energy
  ! hkb         = moist static energy at originating level
  ! mentr_rate  = entrainment rate
  ! cd          = detrainment function for updraft
  ! cdd         = detrainment function for downdraft
  ! dellat      = change of temperature per unit mass flux of cloud ensemble
  ! dellaq      = change of q per unit mass flux of cloud ensemble
  ! dellaqc     = change of qc per unit mass flux of cloud ensemble
  ! aa0         = cloud work function for downdraft
  ! edt         = epsilon
  ! aa0         = cloud work function without forcing effects
  ! aa1         = cloud work function with forcing effects
  ! xaa0        = cloud work function with cloud effects (ensemble dependent)
  ! edt         = epsilon

  !----------------------------------------------------------------
  !
  !     if(ktau.gt.3.and.ktau.lt.7) ...
  !
  day=86400.
  !
  !--- specify entrainmentrate and detrainmentrate
  !
  !     radius=14000.-float(iens)*2000.
  radius=12000.
  !      radius=5000.
  fquasi=1
  fstab=0
  fmconv=0
  !
  !--- gross entrainment rate (these may be changed later on in the
  !--- program, depending what your detrainment is!!)
  !
  !entr_rate=.2/radius
  !entr_rate=1.0E-4 ! ECMWF
  entr_rate=5.0E-5
  !
  !--- entrainment of mass
  !
  !      mentrd_rate=0.

  mentrd_rate=entr_rate
  mentr_rate =entr_rate

  !--- initial detrainmentrates
  !
  do k=1,mkx
     do i=istart,iend
        cd(i,k)  = 0.1*entr_rate
	!qs1= (0.622 * ESAT (T(i,2))  ) / (P(i,2) - ESAT (T(i,2)))
	!qsk= (0.622 * ESAT (T(i,k))  ) / (P(i,k) - ESAT (T(i,k)))

	mentr_rate2(i,k) = entr_rate!*4.*(qsk/qs1)**2

        cdd(i,k) = 0.
     enddo
  enddo
  !
  !
  !--- max/min allowed value for epsilon (ratio downdraft base mass flux/updraft)
  !    base mass flux
  !
  edtmax=.8
  edtmin=.2
  !
  !--- minimum depth (m), clouds must have
  !
  depth_min=500.
  !
  !--- maximum depth (mb) of capping
  !--- inversion (larger cap = no convection)
  !
  !---- for kbcon_cin test
  !  cap_maxs=10.
  !  IF (iens.EQ.2) cap_maxs=20.
  !  IF (iens.EQ.3) cap_maxs=30.
  !  IF (iens.EQ.4) cap_maxs=40.
  !  cap_maxs = 20.
  !
  ! cap_maxs trigger parameter
  !srf- 01/11/2008
  !cap_maxs =  90.
  cap_maxs = 110.
  cap_max_increment(:)=20.

  !--- initialize cap_max
  tkeminnot = 0.08
  tke_start = 0.1 + tkeminnot

  do i=istart,iend
  ! ciclo diurno off
    !cap_max(i)=cap_maxs
  !
  !--- diurnal cycle on cap_maxs:
     fdiur = 1.0
     if( tkeg(i,1) - tkeminnot <  tke_start .and. tkeg(i,1) - tkeminnot > tkmin) fdiur = 0.5
     cap_max(i)= 70. + fdiur * max(0.,cap_maxs - 70.) &
     	      *( tkeg(i,1) - tkmin - tkeminnot) / max(tkeg(i,1), tkmin + tkeminnot)
     !- for output only
     !sgrell9_3d(2,i,j)=cap_max(i)
  end do

  do I=ISTART,IEND
     aa0(i)=0.
     aa1(i)=0.
     aad(i)=0.
     closure_n(i)=16.! = maxens3
     kstabm(i)=mkx-2
     if (aaeq(i).lt.0.) then
        ierr(i)=20
     else
        ierr(i)=0
        XIERR(i,j)=0.
        cwf(i,j)=0.
        pwf(i,j)=0.
        pwdf(i,j)=0.
        eddt(i,j)=0.
        ! xktop(i,j)=0.
        ! xkbas(i,j)=0.
        xmass(i,j)=0.
        predb(i,j)=0.
     endif
       ierr2(i)=ierr(i)
       ierr3(i)=ierr(i)
  enddo


  !--- first check for upstream convection
  if(i_cup_dir_flag == 1 ) then

   do i=istart,iend
     if (ierr(i).eq.0) then
        iresult=0
        massfld=0.
        call cup_direction2(i,j,direction,iact_old_gr,mix,mjx,  &
             mgmxp,mgmyp,massflx,iresult,ensdim,0,0,maxens3,massfld)

  !--- increase cap_max if is there upstream convection
        if (iresult.eq.1) then
           cap_max(i)=cap_max(i)+20.
        endif
     endif
   enddo

  endif

  !--- max height(m) above ground where updraft air can originate

  zkbmax=4000.

  !--- height(m) above which no downdrafts are allowed to originate

  zcutdown=3000.

  !--- depth(m) over which downdraft detrains all its mass

  z_detr=1250.

  !--- MBDT parameter
  !srf - for capmax ensemble mbdt_ens is constant (=4e-3*dtime)
  !      independent of nens
  do nens=1,maxens
     mbdt_ens(nens)=(float(nens)-3.)*dtime*1.e-3+dtime*5.E-03
  enddo
  do nens=1,maxens2
    !edt_ens(nens)=.7-float(nens)*.1
     edt_ens(nens)=.95-float(nens)*.01
  enddo
  !
  !--- environmental conditions, FIRST HEIGHTS
  !
  do i=istart,iend
     if(ierr(i).ne.20)then
        do k=1,maxens*maxens2*maxens3
           xf_ens(  i,j,(iens-1)*maxens*maxens2*maxens3+k)= 0.
           pr_ens(  i,j,(iens-1)*maxens*maxens2*maxens3+k)= 0.
           outt_ens(i,j,(iens-1)*maxens*maxens2*maxens3+k)= 0.
        enddo
     endif
  enddo

  !--- calculate moist static energy, heights, qes
  !
  call cup_env(j, ipr, jpr, z, qes, he, hes, t, q, p, z1, mix, mgmxp,  &
       mkx, mgmzp, istart, iend, psur, ierr, tcrit, 0)
  call cup_env(j, ipr, jpr, zo, qeso, heo, heso, tn, qo, po, z1, mix,mgmxp, &
       mkx, mgmzp, istart, iend, psur, ierr, tcrit, 0)
  !
  !--- environmental values on cloud levels
  !
  call cup_env_clev(j, ipr, jpr, t, qes, q, he, hes, z, p, qes_cup,  &
       q_cup, he_cup, hes_cup, z_cup, p_cup, gamma_cup, t_cup, psur, &
       mix, mgmxp, mkx, mgmzp, istart, iend, ierr, z1)
  call cup_env_clev(j, ipr, jpr, tn, qeso, qo, heo, heso, zo, po,    &
       qeso_cup, qo_cup, heo_cup, heso_cup, zo_cup, po_cup,          &
       gammao_cup, tn_cup, psur, mix, mgmxp, mkx, mgmzp, istart,     &
       iend, ierr, z1)

  do i=istart,iend
     if (ierr(i).eq.0) then

        do k=1,mkx
           if (zo_cup(i,k).gt.zkbmax+z1(i)) then
              kbmax(i)=k
              exit
           endif
        enddo

        !--- level where detrainment for downdraft starts

        do k=1,mkx
           if (zo_cup(i,k) .gt. z_detr+z1(i)) then
              kdet(i)=k
              !GO TO 26
              exit
           endif
        enddo
        !26      CONTINUE

     endif
  enddo

  !--- DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22

  !   here, ,should just take top of PBL for shallow clouds, also maybe for
  !   deep clouds in tropics: kstart=level(pbltop)
  !srf-fev2003
  kstart_way = 0

!srf- TESTAR NO FUTURO -
  if (kstart_way == 1) then
  !--- Determine PBL top using TKE (TKEG) and liquid water mixing ratio (RCPG)
       call get_zi(mix,mgmxp,mkx,mgmzp,istart,iend,j,ierr,kzi,tkeg  &
                ,rcpg,zo,z1,tkmin)
       do i=istart,iend
        !srf-14-fev-2003
         if(ierr(i).eq.0) k22(i) = kzi(i)
        !Esta produz uma cobertura de shallow mais representativa
        !Em caso de alta resolucao vertical tente a versao original (k22=kzi)
        !if(ierr(i).eq.0) k22(i) = max(2, kzi(i) - 1)
        !	    IF(ierr(i).eq.0) PRINT*,k22(i)
       enddo

  else
  ! - old way to define k22
      kzi(1:mgmxp) = 1
      kstart = 3
      CALL MAXIMI(HEO_CUP,mix,mgmxp,mkx,mgmzp,kstart,KBMAX,K22, &
                 ISTART,IEND,ierr)
  endif


  do I=ISTART,IEND
     if (ierr(I).eq.0) then
        if (K22(I).ge.KBMAX(i)) ierr(i)=2
     endif
  enddo

!srf-print----------------
!	do I=ISTART,IEND
!	   if(j==16 .and. i==95) then
!	      print*,' HES & HES_CUP profiles'
!	    do k=1,mkx
!	     print*,k,heo_cup(i,k),heso_cup(i,k)
!	    enddo
!	   endif
!	 enddo
!srf-print----------------



  !--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON

   if (trigg == 1) then

      call cup_kbcon_catt(cap_max_increment,1,k22,kbcon,heo_cup,heso_cup, &
           hkbo,kzi,mix,mgmxp,mkx,mgmzp,istart,iend,ierr,kbmax,po_cup,cap_max, j)

   elseif(trigg == 2) then

      CALL trigg_ecmwf(mentr_rate2,psur,dnsup, istart,iend,mgmxp,mgmzp,mkx,t_cup, q_cup, qes_cup,z_cup,&
	    p_cup, sflux_rg,sflux_tg,t,q,qes,z,ierr,k22,kbcon,ktop,m1,m2,m3,j,wu)
      !- sgrell9_3d = velocidade vertical no interior das nuvens (equacao 6, Jakob e Siebesma 2003, MWR).

      !if(mynum == 1)print*,'tr1=',ktop(2+14),kbcon(2+14),k22(2+14),ierr(2+14);call flush(6)
      !if(mynum == 2)print*,'tr2=',ktop(2),kbcon(2),k22(2),ierr(2);call flush(6)

      !do k=1,mkx
      !  kr = k + 1  ! nivel k da grade do grell corresponde ao nivel k + 1 do rams
      !  do i = istart,iend
      !    wu(i,k) = sgrell9_3d(kr,i,j)
      !  enddo
      !enddo
   else

    stop 'trigger function must be 1 or 2'

   endif


   do i=istart,iend
    if(ierr(i).eq.0) hkb(i)=hkbo(i)
   enddo

  !srf 30-jan-2002. This routine blow up the model
  !CALL cup_kbcon_cin(1, k22,kbcon, heo_cup, heso_cup, z, tn_cup,  &
  !     qeso_cup, mix, mgmxp, mkx, mgmzp, istart, iend, ierr,      &
  !     kbmax, po_cup, cap_max)

  !--- Increase detrainment in stable layers
  call MINIMI(HEso_cup, mix, mgmxp, mkx, mgmzp, Kbcon, kstabm,  &
       kstabi, ISTART, IEND, ierr)

  do I=ISTART,IEND
     if (ierr(I).eq.0.) then
        if (kstabm(i)-1.gt.kstabi(i)) then
           do k=kstabi(i),kstabm(i)-1
	     ! cd(i,k)=cd(i,k-1)+1.5*entr_rate
              cd(i,k)=cd(i,k-1)+1.5*mentr_rate2(i,k)

	      !if (iens.gt.4) then
              !  ! cd(i,k) = cd(i,k-1)+float(iens-4)*entr_rate/  &
              !  !           float(kstabm(i)-kstabi(i))
              !
	      !	 cd(i,k) = cd(i,k-1)+float(iens-4)*mentr_rate2(i,k)/  &
              !             float(kstabm(i)-kstabi(i))
	      !else
              !   cd(i,k)=cd(i,k)
              !endif
              !if (cd(i,k).gt.10.0*entr_rate) cd(i,k)=10.0*entr_rate
               if (cd(i,k).gt.10.0*mentr_rate2(i,k)) cd(i,k)=10.0*mentr_rate2(i,k)
	   enddo
        endif
     endif
  enddo


  !--- Calculate incloud moist static energy

  call cup_up_he(k22, hkb, z_cup, cd, mentr_rate2, he_cup, hc,   &
       mix, mgmxp, mkx, mgmzp, kbcon, ierr, istart, iend, dby,  &
       he, hes_cup)
  call cup_up_he(k22, hkbo, zo_cup, cd, mentr_rate2, heo_cup, hco,  &
       mix, mgmxp, mkx, mgmzp, kbcon, ierr, istart, iend, dbyo,    &
       heo, heso_cup)

  !--- Determine cloud top - KTOP
  IF (trigg == 1) then
     call cup_ktop(1, dbyo, kbcon, ktop, mix, mgmxp, mkx, mgmzp, istart, &
          iend, ierr)
  ELSE
      DO I=ISTART,IEND
        if (ierr(I).eq.0) then
          do k=ktop(i)+1,mkx
             dby(i,k)=0.
          enddo
        endif
    enddo
  ENDIF

  !srf--print--------------------
  !DO I=ISTART,IEND
  !   if(i.eq.ipr.and.j.eq.jpr) then
  !     print*,'-----------K T O P--2----------'
  !     print*,'MYNUM I J=',mynum,i,j
  !     print*,'k22-kbcon-ierr-Ktop=',k22(i),kbcon(i),ierr(i),ktop(i)
  !     do K=KBCON(I)+1,MKX-2
  !        print*,'k-heo_cup-dbyo=',k,heo_cup(i,k),dbyo(i,k)
  !     enddo
  !     print*,'-------------------------------'
  !   endif
  !enddo
  !srf--print--------------------

  do I=ISTART,IEND
     kzdown(i)=0
     if (ierr(i).eq.0) then
        zktop = (zo_cup(i,ktop(i))-z1(i))*.6
        zktop = min(zktop+z1(i),zcutdown+z1(i))
        do k=1,mkx
           if (zo_cup(i,k).gt.zktop) then
              kzdown(i)=k
              exit
           endif
        enddo
     endif
  enddo


  !--- Downdraft originating level - JMIN
  call MINIMI(HEso_cup, mix, mgmxp, mkx, mgmzp, K22, kzdown,  &
       JMIN, ISTART, IEND, ierr)

  !srf--print-------------------
  !      DO I=ISTART,IEND
  !       if(i.eq.ipr.and.j.eq.jpr) then
  !       print*,'-------------------------------'
  !          write(mynum+2,111) mynum,i0+i,j0+j,Jmin(i),Ktop(i)
  ! 111      format(1x,5i5)
  !       print*,'MYNUM I J=',mynum,i,j
  !       print*,'k22-kbcon-ierr=',k22(i),kbcon(i),ierr(i)
  !       print*,'ktop-he_cup-dbyo=',ktop(i),he_cup(i,12),dbyo(i,12)
  !       print*,'Ktop Jmin=',ktop(i),Jmin(i)
  !       print*,'-------------------------------'
  !       endif
  !      enddo
  !srf--print--------------------

  do I=ISTART,IEND
     if (ierr(I).eq.0.) then

        !--- Check whether it would have buoyancy, if there where
        !--- no entrainment/detrainment

101     continue

	!if(i==2 .and. j==3) print*,'x',jmin(i),kdet(i),ktop(i)

	if (jmin(i)-1 .lt. kdet(i)  ) kdet(i) = jmin(i)-1
        if (jmin(i)   .ge. ktop(i)-1) jmin(i) = ktop(i)-2
        ki=jmin(i)
	!

        if (JMIN(I) == 0) then
           ierr(i)=9
	    GO TO 100
        endif

	!
	!
        hcdo(i,ki) = heso_cup(i,ki)
        DZ         = Zo_cup(i,Ki+1)-Zo_cup(i,Ki)
        !dh         = dz*(HCDo(i,Ki)-heso_cup(i,ki))
        dh         = 0.

        do k=ki-1,1,-1
           hcdo(i,k) = heso_cup(i,jmin(i))
           DZ        = Zo_cup(i,K+1)-Zo_cup(i,K)
           dh        = dh+dz*(HCDo(i,K)-heso_cup(i,k))
           if (dh .gt. 0.) then
              jmin(i)=jmin(i)-1
              if (jmin(i) .gt. 3     ) then
                 GO TO 101
              else if (jmin(i) .le. 3) then
                 ierr(i)=9
                 GO TO 100
              endif
           endif
        enddo


     endif
100  continue
  enddo


  !--- Must have at least depth_min m between cloud convective
  !    base and cloud top.

  do i=istart,iend
     if (ierr(i).eq.0.) then
        if (-zo_cup(i,kbcon(i))+zo_cup(i,ktop(i)) .lt. depth_min) ierr(i)=6
     endif
  enddo

  !--- Normalized updraft mass flux profile

  call cup_up_nms(zu, z_cup, mentr_rate2, cd, kbcon, ktop,  &
       mix, mgmxp, mkx, mgmzp, istart, iend, ierr, k22)
  call cup_up_nms(zuo, zo_cup, mentr_rate2, cd, kbcon, ktop,  &
       mix, mgmxp, mkx, mgmzp, istart, iend, ierr, k22)

  !--- Normalized downdraft mass flux profile,also work on
  !    bottom detrainment
  !--- in this routine

  call cup_dd_nms(zd, z_cup, cdd, mentrd_rate, jmin, ierr,    &
       mix, mgmxp, mkx, mgmzp, istart, iend, 0, kdet, z1)
  call cup_dd_nms(zdo, zo_cup, cdd, mentrd_rate, jmin, ierr,  &
       mix, mgmxp, mkx, mgmzp, istart, iend, 1, kdet, z1)

  !--- Downdraft moist static energy

  call cup_dd_he(hes_cup, zd, hcd, z_cup, cdd, mentrd_rate,  &
       jmin, ierr, mix, mgmxp, mkx, mgmzp, istart, iend,he,  &
       kdet, dbyd, he_cup)
  call cup_dd_he(heso_cup, zdo, hcdo, zo_cup, cdd,           &
       mentrd_rate, jmin, ierr, mix, mgmxp, mkx, mgmzp,      &
       istart, iend,heo, kdet, dbydo, he_cup)

  !--- Calculate moisture properties of downdraft

  call cup_dd_moisture(j, zd, hcd, hes_cup, qcd, qes_cup,  &
       pwd, q_cup, z_cup, cdd, mentrd_rate, jmin, ierr,    &
       gamma_cup, pwev, mix, mgmxp, mkx, mgmzp, istart,    &
       iend, bu, qrcd, q, he, hc, t_cup, 2)
  call cup_dd_moisture(j, zdo, hcdo, heso_cup, qcdo,       &
       qeso_cup, pwdo, qo_cup, zo_cup, cdd, mentrd_rate,   &
       jmin, ierr, gammao_cup, pwevo, mix, mgmxp, mkx,     &
       mgmzp, istart, iend, bu, qrcdo, qo, heo, hco, tn_cup, 1)

  !--- Calculate moisture properties of updraft

  !call cup_up_moisture(ierr, z_cup, qc, qrc, pw, pwav,     &
  !     kbcon, ktop, mix, mgmxp, mkx, mgmzp, istart, iend,  &
  !     cd, dby, mentr_rate2, q, GAMMA_cup, zu, qes_cup,     &
  !     k22, q_cup)
  !call cup_up_moisture(ierr, zo_cup, qco, qrco, pwo, pwavo,  &
  !     kbcon, ktop, mix, mgmxp, mkx, mgmzp, istart, iend,    &
  !     cd, dbyo,                                             &
  !     mentr_rate2, qo, GAMMAo_cup, zuo, qeso_cup, k22, qo_cup)

 call cup_up_moisture(ierr, z_cup, qc, qrc, pw, pwav,     &
       kbcon, ktop, mix, mgmxp, mkx, mgmzp, istart, iend,  &
       cd, dby, mentr_rate2, q, GAMMA_cup, zu, qes_cup,     &
       k22, q_cup, wu, dn0_2,ccn,TRIGG,iens,autoconv)

 call cup_up_moisture(ierr, zo_cup, qco, qrco, pwo, pwavo,  &
       kbcon, ktop, mix, mgmxp, mkx, mgmzp, istart, iend,    &
       cd, dbyo,mentr_rate2, qo, GAMMAo_cup, zuo, qeso_cup, k22,&
       qo_cup,  wu, dn0_2,ccn,TRIGG,iens,autoconv)



  !--- Calculate workfunctions for updrafts

  call cup_up_aa0(aa0, z, zu, dby, GAMMA_CUP, t_cup, kbcon,  &
       ktop, mix, mgmxp, mkx, mgmzp, istart, iend, ierr)
  call cup_up_aa0(aa1, zo, zuo, dbyo, GAMMAo_CUP, tn_cup,    &
       kbcon, ktop, mix, mgmxp, mkx, mgmzp, istart, iend, ierr)
  do i=istart,iend
     if (ierr(i).eq.0) then
        if (aa1(i).eq.0.) ierr(i)=17
     endif
  enddo

  !--- Determine downdraft strength in terms of windshear

  call cup_dd_edt(ierr, us, vs, zo, ktop, kbcon, edt, po,   &
       pwavo, pwevo, mix, mgmxp, mkx, mgmzp, istart, iend,  &
       edtmax, edtmin, maxens2, edtc, vshear, sdp, vws)

!- srf - big loop starts here! ---------------------------------------------

  do iedt=1,maxens2 !orig: DO 250 iedt=1,maxens2

     do i=istart,iend
        if (ierr(i).eq.0) then
           edt(i)  = edtc(i,iedt)
           edto(i) = edtc(i,iedt)
           edtx(i) = edtc(i,iedt)
        endif
     enddo
     do k=1,mkx
        do i=istart,iend
            dellat_ens(i,k,iedt) = 0.
            dellaq_ens(i,k,iedt) = 0.
           dellaqc_ens(i,k,iedt) = 0.
               pwo_ens(i,k,iedt) = 0.
        enddo
     enddo
!      if(j.eq.jpr)then
!      if(j.eq.jpr)then
!         i=ipr
!        write(0,*)'in 250 loop ',iedt,edt(ipr),ierr(ipr)
!       if(ierr(i).eq.0.or.ierr(i).eq.3)then
!         write(0,*)'250',k22(I),kbcon(i),ktop(i),jmin(i)
!         write(0,*)edt(i),aa0(i),aa1(i)
!         do k=1,mkx
!           write(0,*)k,z(i,k),he(i,k),hes(i,k)
!         enddo
!         write(0,*)'end 250 loop ',iedt,edt(ipr),ierr(ipr)
!         do k=1,ktop(i)+1
!           write(0,*)zu(i,k),zd(i,k),pw(i,k),pwd(i,k)
!         enddo
!        endif
!      endif

     !--- downdraft workfunction
!srf- comentado em mar/2004
     !x     do I=ISTART,IEND
     !x      aad(i)=0.
     !x     enddo
     !x     call cup_dd_aa0(edto, ierr, aad, jmin, gammao_cup, tn_cup, &
     !x          hcdo, heso_cup, zo, mix, mgmxp, mkx, mgmzp, istart,   &
     !x          iend, zdo)
!srf-

     !--- Change per unit mass that a model cloud would modify the environment

     !--- 1. in bottom layer

     call cup_dellabot_catt(0, 0, heo_cup, ierr, zo_cup, po_cup,                &
                            hcdo, edto, zdo, cdd, heo, mix, mgmxp, mkx, mgmzp,  &
                            istart, iend, dellah, 1, j, mentrd_rate, zo)
     call cup_dellabot_catt(ipr, jpr, qo_cup, ierr, zo_cup,                 &
                            po_cup, qrcdo, edto, zdo, cdd, qo, mix, mgmxp,  &
                            mkx, mgmzp, istart, iend, dellaq, 2, j,         &
                            mentrd_rate, zo)

     !--- 2. everywhere else

     call cup_dellas_catt(ierr, zo_cup, po_cup, hcdo, edto, zdo,             &
                          cdd, heo, mix, mgmxp, mkx, mgmzp, istart, iend,    &
                          dellah, 1, j, mentrd_rate, zuo, cd, hco, ktop,     &
                          k22, kbcon, mentr_rate2, jmin, heo_cup, kdet, k22,  &
                          0, 0, 'deep')

     !-- Take out cloud liquid water for detrainment

     do k=1,mkx
        do i=istart,iend
           scr1(i,k)   = 0.
           dellaqc(i,k)= 0.

	   if (ierr(i).eq.0) then

	      scr1(i,k)=qco(i,k)-qrco(i,k)

              if (k.eq.ktop(i)-0) dellaqc(i,k) = .01*zuo(i,ktop(i))*  &
                                  qrco(i,ktop(i))*g/(po_cup(i,k)-po_cup(i,k+1))

              if (k.lt.ktop(i).and.k.gt.kbcon(i)) then

		 dz = zo_cup(i,k+1)-zo_cup(i,k)
                 dellaqc(i,k) = .01*g*cd(i,k)*dz*zuo(i,k)*            &
                                .5*(  qrco(i,k) +   qrco(i,k+1) )/    &
                                   (po_cup(i,k) - po_cup(i,k+1) )
              endif
           endif
        enddo
     enddo

     call cup_dellas_catt(ierr, zo_cup, po_cup, qrcdo, edto, zdo,  &
          cdd, qo, mix, mgmxp, mkx, mgmzp, istart, iend,      &
          dellaq, 2, j, mentrd_rate, zuo, cd, scr1, ktop,     &
          k22, kbcon, mentr_rate2, jmin, qo_cup, kdet, k22,    &
          ipr, jpr, 'deep')

     !srf-----print-------
     !do i=istart,iend
     !    if (ierr(i).eq.0) then
     !       if (j.eq.jpr.and.i.eq.ipr) then
     !  	dellaqsum  =0.
     !  	dellaqcsum =0.
     !  	do k=1,mkx
     !  	   dp	      = -100.*(p_cup(i,k+1)-p_cup(i,k))
     !  	   !	     dp=100.*(po(i,k-1)-po(i,k))
     !  	   !	     if(k.eq.1) dp=996.2952
     !  	   dellaqsum  = dellaqsum + dellaq(i,k) * dp/9.81
     !  	   dellaqcsum = dellaqcsum + dellaqc(i,k)* dp/9.81
     !
     !  	   if (k.eq.1) then
     !  	      write(6,'(a1,78a1)') ' ',('-',m=1,78)
     !  	    print *, 'i k Z_cup P_cup dp dellaq dellaq',  &
     !  		   ' dellaqsum dellaqcsum'
     !  	   endif
     !  	   write(6,'(2i4,5F12.4,3f16.4)') i, k, 	  &
     !  		Z_cup(i,k),p_cup(i,k), dp,		  &
     !  		1.e+3*86400.*dellaq(i,k),		  &
     !  		1.e+3*86400.*dellaqc(i,k),		  &
     !  		1.e+3*86400.*dellaqsum, 		  &
     !  		1.e+3*86400.*dellaqcsum,		  &
     !  		1.e+3*86400.*(dellaqsum+dellaqcsum)
     !  	   !& 100.*(p_cup(i,k-1)-p_cup(i,k))/
     !  	   !& (g*(z_cup(i,k)-z_cup(i,k-1)))
     !  	   ! densidade
     !  	enddo
     !      endif
     !   endif
     !enddo
     !srf-----print-------

     !--- Using dellas, calculate changed environmental profiles

        mbdt=mbdt_ens(2)
        do i=istart,iend
           do k=1,maxens
              xaa0_ens(i,k)=0.
          enddo
        enddo

        do k=1,mkx-1
           do i=istart,iend
              dellat(i,k)=0.
              if (ierr(i).eq.0) then
                 XHE     (I,K) = DELLAH(I,K)*MBDT + HEO(I,K)
                 XQ      (I,K) = DELLAQ(I,K)*MBDT +  QO(I,K)
                 DELLAT  (I,K) = (1./1004.)*(DELLAH(I,K)-2.5E06*DELLAQ(I,K))
                 XT_Grell(I,K) = DELLAT(I,K)*MBDT +  TN(I,K) !(XT_Grell == old XT)
                 if (XQ(I,K).le.0.) XQ(I,K)=1.E-08
                 !if (i.eq.ipr.and.j.eq.jpr) then
                 !   print *, k, DELLAH(I,K), DELLAQ(I,K), DELLAT(I,K)
                 !endif
              endif
           enddo
        enddo
        !
        do i=istart,iend
           if (ierr(i).eq.0) then
              XHE     (I,mkx)  = HEO(I,mkx)
              XQ      (I,mkx)  = QO (I,mkx)
              XT_Grell(I,mkx)  = TN (I,mkx)
              if (XQ(I,mkx).le.0.) XQ(I,mkx) = 1.E-08
           endif
        enddo

        !--- Calculate moist static energy, heights, qes

        call cup_env(j, ipr, jpr, xz, xqes, xhe, xhes, xt_grell,  &
             xq, po, z1, mix, mgmxp, mkx, mgmzp, istart,    &
             iend, psur, ierr, tcrit, 2)

        !--- Environmental values on cloud levels

        call cup_env_clev(j, ipr, jpr, xt_grell, xqes, xq, xhe,   &
             xhes, xz, po, xqes_cup, xq_cup, xhe_cup,       &
             xhes_cup, xz_cup, po_cup, gamma_cup, xt_cup,   &
             psur, mix, mgmxp, mkx, mgmzp, istart, iend,    &
             ierr, z1)

        !**************************** Static Control

        !--- Moist static energy inside cloud

        do i=istart,iend
           if (ierr(i).eq.0) then
              xhkb(i)=xhe(i,k22(i))
           endif
        enddo
        call cup_up_he(k22, xhkb, xz_cup, cd, mentr_rate2,  &
             xhe_cup, xhc, mix, mgmxp, mkx, mgmzp, kbcon,  &
             ierr, istart, iend, xdby, xhe, xhes_cup)

        !--- Normalized mass flux profile

        call cup_up_nms(xzu, xz_cup, mentr_rate2, cd, kbcon,  &
             ktop, mix, mgmxp, mkx, mgmzp, istart, iend,     &
             ierr, k22)
        call cup_dd_nms(xzd, xz_cup, cdd, mentrd_rate, jmin, &
             ierr, mix, mgmxp, mkx, mgmzp, istart, iend, 1,  &
             kdet, z1)

        !--- Moisture downdraft

        call cup_dd_he(xhes_cup, xzd, xhcd, xz_cup, cdd,     &
             mentrd_rate, jmin, ierr, mix, mgmxp, mkx,       &
             mgmzp, istart, iend, xhe, kdet, dbyd, xhe_cup)
        call cup_dd_moisture(j, xzd, xhcd, xhes_cup, xqcd,   &
             xqes_cup, xpwd, xq_cup, xz_cup, cdd,            &
             mentrd_rate, jmin, ierr, gamma_cup, xpwev,mix,  &
             mgmxp, mkx, mgmzp, istart, iend, bu, xqrcd, xq, &
             xhe, xhc, xt_cup,3)

        !--- Moisture updraft

        !call cup_up_moisture(ierr, xz_cup, xqc, xqrc, xpw,   &
        !     xpwav, kbcon, ktop, mix, mgmxp, mkx, mgmzp,     &
        !     istart, iend, cd, xdby, mentr_rate2, xq,         &
        !     GAMMA_cup, xzu, xqes_cup, k22, xq_cup)

	call cup_up_moisture(ierr, xz_cup, xqc, xqrc, xpw,   &
             xpwav, kbcon, ktop, mix, mgmxp, mkx, mgmzp,     &
             istart, iend, cd, xdby, mentr_rate2, xq,         &
             GAMMA_cup, xzu, xqes_cup, k22, xq_cup, &
	     wu, dn0_2,ccn,TRIGG,iens,autoconv)


	!
        !--- Workfunctions for updraft
        !
        call cup_up_aa0(xaa0, xz, xzu, xdby, GAMMA_CUP,   &
             xt_cup, kbcon, ktop, mix, mgmxp, mkx, mgmzp, &
             istart, iend,ierr)

        !--- Workfunctions for downdraft
        !call cup_dd_aa0(edtx,ierr,xaa0,jmin,gamma_cup,xt_cup, &
        !     xhcd,xhes_cup,xz,mix,mgmxp,mkx,mgmzp,istart,iend,xzd)

!------------------------- loop at  cap_max ensemble ----------
      do nens=1,maxens


        do i=istart,iend
           if (ierr(i).eq.0) then

              xaa0_ens(i,nens) = xaa0(i)
              nall=  (iens-1)*maxens3*maxens*maxens2  &
     		   + (iedt-1)*maxens3*maxens	    &
      		   + (nens-1)*maxens3

              do k=1,mkx
                 if (k.le.ktop(i)) then
                    do nens3=1,maxens3

!- pr_ens is the total precipitation on the surface associated to each ensemble member

                       if (nens3.eq.7) then
                   !--- b=0
                          pr_ens(i,j,nall+nens3) = pr_ens(i,j,nall+nens3) +  &
                                                   pwo(i,k)+edto(i)*pwdo(i,k)
                   !--- b=beta
                       else if (nens3.eq.8) then
                          pr_ens(i,j,nall+nens3) = pr_ens(i,j,nall+nens3) +  &
                                                   pwo(i,k)
                   !--- b=beta/2
                       else if (nens3.eq.9) then
                          pr_ens(i,j,nall+nens3) = pr_ens(i,j,nall+nens3) +  &
                                                   pwo(i,k)+.5*edto(i)*pwdo(i,k)
                       else
                          pr_ens(i,j,nall+nens3) = pr_ens(i,j,nall+nens3) +  &
                                                   pwo(i,k)+edto(i)*pwdo(i,k)
                       endif

                    enddo
                 endif
              enddo

              if(pr_ens(i,j,nall+7) .lt. 0.) then
                 ierr(i) = 18
                 do nens3=1,maxens3
                    pr_ens(i,j,nall+nens3)= 0.
                 enddo
	      endif

	      do nens3=1,maxens3
	         if(pr_ens(i,j,nall+nens3) .lt. 0.) then
                    pr_ens(i,j,nall+nens3)= 0.
	         endif
	      enddo
              do nens3=1,maxens3
                 outt_ens(i,j,nall+nens3) = dellat(i,1)
              enddo

           endif
        enddo

     enddo
!---------------------end of loop at  cap_max ensemble ----------

     !--- Large scale forcing

     !   here, ,should just take top of PBL for shallow clouds, also maybe for
     !   deep clouds in tropics: kstart=level(pbltop)
     !kstart = 3
     !call MAXIMI(HE_CUP, mix, mgmxp, mkx, mgmzp,kstart,   &
     !     KBMAX, K22x, ISTART, IEND, ierr)

     do I=ISTART,IEND
        if (ierr(i).eq.0) then
           !if (K22x(I).ge.KBMAX(i)) ierr(i)=998
	   k22x(i)=k22(i)
        endif
        ierr2(i)=ierr(i)
        ierr3(i)=ierr(i)
     enddo

     !--- Determine the level of convective cloud base - KBCON

     if (trigg == 1) then

       call cup_kbcon_catt(cap_max_increment,2,k22x,kbconx,heo_cup,heso_cup, &
          hkbo,kzi,mix,mgmxp,mkx,mgmzp,istart,iend,ierr2,kbmax,po_cup,cap_max,&
	  j)

       call cup_kbcon_catt(cap_max_increment,3,k22x,kbconx,heo_cup,heso_cup, &
          hkbo,kzi,mix,mgmxp,mkx,mgmzp,istart,iend,ierr3,kbmax,po_cup,cap_max,&
	  j)

     endif



     if(ensdim.eq.1)then

        call cup_forcing_ens_1_catt(aa0,aa1,xaa0_ens,mbdt_ens,dtime,    &
             xmb,ierr,mix,mgmxp,mjx,mgmyp,mkx,mgmzp,istart,iend,xf_ens, &
             j,fquasi,fstab,'deeps',xland,maxens,iens,iedt,maxens2,ipr, &
             jpr,maxens3,mconv,omeg,zdo,kbcon,zuo,pr_ens,edto,aad,kbcon, &
             massflx,iact_old_gr,direction,ensdim, &
             massfln,xff_ens3, xk,ierr2,ierr3, i_cup_dir_flag)

     elseif(maxens3.eq.16) then

        call cup_forcing_ens_16_catt(aa0,aa1,xaa0_ens,mbdt_ens,dtime,     &
             xmb,ierr,mix,mgmxp,mjx,mgmyp,mkx,mgmzp,istart,iend,xf_ens,j, &
             fquasi,fstab,'deeps',xland,maxens,iens,iedt,maxens2,ipr,jpr, &
             maxens3,mconv,omeg,zdo,kbcon,zuo,pr_ens,edto,aad,kbcon,      &
             massflx,iact_old_gr,direction,ensdim,                        &
             massfln,massfld,iresult,xff_ens3, xk,p_cup,ktop,icoic,ierr2, &
             ierr3,ngrid,                                                 &
!	     sgrell1_2d, &
	     m2,m3, tscl_KF, i_cup_dir_flag,                   &
	     glatg,glong,closure_n,xcoast)

     endif

     do k=1,mkx
        do i=istart,iend
           if (ierr(i).eq.0) then
              dellat_ens (i,k,iedt)  =  dellat(i,k)
              dellaq_ens (i,k,iedt)  =  dellaq(i,k)
              dellaqc_ens(i,k,iedt) =  dellaqc(i,k)
              pwo_ens    (i,k,iedt)  = pwo(i,k) + edt(i)*pwdo(i,k)

          else

              dellat_ens (i,k,iedt) = 0.
              dellaq_ens (i,k,iedt) = 0.
              dellaqc_ens(i,k,iedt) = 0.
              pwo_ens    (i,k,iedt) = 0.
           endif
        enddo
     enddo
     !250  CONTINUE
  enddo

  !--- FEEDBACK

  call cup_output_ens_catt(xf_ens,ierr,dellat_ens,dellaq_ens,dellaqc_ens,   &
       outt,outq,outqc,pre,pwo_ens,xmb,ktop,mix,mgmxp,mjx,mgmyp,mkx,mgmzp,  &
       istart,iend,j,'deep',maxens2,maxens,ipr,jpr,iens,pr_ens,outt_ens,    &
       maxens3,ensdim,icoic,massfln,xfac1,xfac_for_dn, &
!       sgrell1_3d,sgrell2_3d,&
       m1,m2,m3,ngrid,dti,closure_n,xland)


!  do I=ISTART,IEND
!     if(ierr(i).eq.0)then
!      masscon=0.
!      do k=1,ktop(i)
!       masscon=masscon+(outq(i,k)+outqc(i,k))*100.*(p_cup(i,k+1)-p_cup(i,k))/9.81
!      enddo
!
!      print*,'masscon pre=' ,i,masscon,pre(i),100.*(masscon-pre(i))/(1.e-10+masscon)
!      endif
!  enddo


  do i=istart,iend
     pre(i) = max(pre(i),0.)
  enddo


  !
  !---------------------done----------------------------------------------------------------
  if(CCATT == 1) then

  !---Salva parametros para uso no transporte convectivo:
    do i=istart,iend
       ierr4d(i)   = ierr(i)
       jmin4d(i)   = jmin(i)
       kdet4d(i)   = kdet(i)
       k224d(i)    = k22(i)
       kbcon4d(i)  = kbcon(i)
       ktop4d(i)   = ktop(i)
       kstabi4d(i) = kstabi(i)
       kstabm4d(i) = kstabm(i)
       !kpbl4d(i)   = kpbl(i)
       kpbl4d(i)   = k22(i) ! por enquanto kpbl==k22
       xmb4d(i)    = xmb(i)

       ! downdraft mass flux averaged
       dnmf(i,j)  = 0.
       do k=1,ensdim
          dnmf(i,j) = dnmf(i,j) + massfln(i,j,k)
       enddo
       dnmf(i,j) = dnmf(i,j)/float(ensdim)
        !recalcula o parametro edt4d para o transporte convectivo
        !modifique posteriormente para uso direto do fluxo de massa do downdraft
        !dnmf(i,j) = edt4d(i)*xmb(i)    ! downdraft mass flux averaged
       edt4d(i)  = dnmf(i,j) / ( xmb(i) + 1.e-16 )

       do k=1,mkx
 	  if (iens.eq.1) then
 	     zcup5d(k,i) = zo_cup(i,k)
 	     pcup5d(k,i) = po_cup(i,k)
 	  endif
 	  enup5d(k,i) = mentr_rate2(i,k)
 	  endn5d(k,i) = mentrd_rate
 	  deup5d(k,i) =   cd(i,k)
 	  dedn5d(k,i) =  cdd(i,k)
 	   zup5d(k,i)  = zuo(i,k)
 	   zdn5d(k,i)  = zdo(i,k)


   !--tmp
   !if (iens.eq.1) then
   !extra3d(10,ngrid)%d3(k,i,j)= mentr_rate2(i,k)!
   !extra3d(11,ngrid)%d3(k,i,j)= mentrd_rate
   !extra3d(12,ngrid)%d3(k,i,j)=   cd(i,k)
   !extra3d(13,ngrid)%d3(k,i,j)=  cdd(i,k)
   !extra3d(14,ngrid)%d3(k,i,j)= zuo(i,k)
   !extra3d(15,ngrid)%d3(k,i,j)= zdo(i,k)
   !endif
   !tmp-

 	   ! p-lw5d is the ratio between precip
 	   ! rate/liq water content after rainout which
 	   ! unit is:  kg[ar]/(m^2s)
 	   ! It's save for use at wet-deposition scheme
 	   ! (in-cloud)
 	   prup5d (k,i) = xmb(i)*pwo(i,k) !only for upfradt
 	   clwup5d(k,i) =  qrco(i,k)	  !only for upfradt
 	   tup5d  (k,i) = t_cup(i,k)  !  em verdade deveria ser a temperatura da parcela
 				      !>>> de ar no updraft e _NAO_ a temperatura ambiente
 	enddo
     enddo
  endif
  !-------Salva parametros nas analises do RAMS

  do I=ISTART,IEND
     xierr(i,j)=float(ierr(i))
     if (ierr(i).eq.0) then
        upmf(i,j)   = xmb(i)

        ! downdraft mass flux averaged
        dnmf(i,j)  = 0.
        do k=1,ensdim
           dnmf(i,j) = dnmf(i,j) + massfln(i,j,k)
        enddo
        dnmf(i,j) = dnmf(i,j)/float(ensdim)
        !
        xktop(i,j ) = float(ktop(i) )
        xkbcon(i,j) = float(kbcon(i))
        xkdt(i,j  ) = float(kdet(i) )
        xjmin(i,j ) = float(jmin(i) )
	xk22(i,j  ) = float(k22(i)  )
     elseif (ierr(i).ne.0.and.ierr(i).ne.20) then
        upmf(i,j)   = 0.
        dnmf(i,j)   = 0.
        xktop(i,j)  = 0.
        xkbcon(i,j) = 0.
        xkdt(i,j)   = 0.
        xjmin(i,j)  = 0.
        xk22(i,j)   = 0.
	!sgrell9_3d(:,i,j)=0.  !<<<<<<<<<<<<<<<<<< 15/12/2010 Boulder
     endif
  enddo

end subroutine CUP_enss_catt
!--------------------------------------------------------------------------------

subroutine cup_output_ens_catt(xf_ens,ierr,dellat,dellaq,dellaqc,           &
           outt,outq,outqc,pre,pwo_ens,xmb,ktop,mix,mgmxp,mjx,mgmyp,        &
           mkx,mgmzp,istart,iend,j,name,maxens2,maxens,ipr,jpr,iens,pr_ens, &
           outt_ens,maxens3,ensdim,icoic,massfln,xfac1,xfac_for_dn,         &
!           sgrell1_3d,sgrell2_3d,
	   m1,m2,m3,ngrid, dti,closure_n,xland)


 use cup_output_vars, only: &
       xmb_ave,xmb_std,xmb_ske,xmb_cur, &
        pr_ave, pr_std, pr_ske, pr_cur, &
	 x_ave,  x_std,  x_ske,  x_cur, &
	 x_ave_cap,                     &
	 x_ave_cap1,x_ave_cap2,x_ave_cap3,x_ave_cap4,x_ave_cap5, &
        cup_output_vars_alloc,          &
        alloc_cup_output_vars

  implicit none
  ! xf = ensemble mass fluxes
  ! pr_ens = surface precipitation ensembles
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble
  ! outt = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
  ! xmb    = total base mass flux
  ! xfac1  = correction factor
  ! pwo_ens = pw -epsilon*pd (ensemble dependent)
  ! ierr error value, maybe modified in this routine
  !
  character (LEN=*):: name  !CHARACTER *(*) name

  integer mix,mjx,mkx,istart,iend,mgmxp,mgmyp,mgmzp,    &
          ensdim,i,k,j,maxens2,n,maxens,ipr,jpr,        &
          icoic,m1,m2,m3,ngrid,ncount,iens,maxens3

  integer:: ktop(mgmxp),ierr(mgmxp)

  real ::outtes,ddtes,dtt,dtq,dtqc,dtpw

  real :: xf_ens(mgmxp,mgmyp,ensdim), pr_ens(mgmxp,mgmyp,ensdim) &
       ,outt_ens(mgmxp,mgmyp,ensdim),massfln(mgmxp,mgmyp,ensdim)

  real ::outt(mgmxp,mgmzp),outq(mgmxp,mgmzp),outqc(mgmxp,mgmzp),       &
       dellat(mgmxp,mgmzp,maxens2), dellaq(mgmxp,mgmzp,maxens2),     &
       pwo_ens(mgmxp,mgmzp,maxens2),dellaqc(mgmxp,mgmzp,maxens2),     &
       pre(mgmxp),xmb(mgmxp),xfac1(mgmxp),xfac_for_dn(mgmxp)

  !real, dimension(m1,m2,m3) :: sgrell1_3d,sgrell2_3d
  !real, dimension(m2,m3) :: train
  real ::dti
  !srf -dec/2003
  ! Tunning parameter = use it to calibrate the precipitation field.
  real, parameter :: tunning = 0. ! tunning = 0.7

  ! Statistics properties
  integer, parameter :: i_use_stat_prop = 1
  integer, parameter :: calc_xmb_with_statis_prop = 0 ! = 0 use simple aritmetic mean
                                                      ! = 1 use ensemble averages
  real:: closure_n(mgmxp),clos_wei
  real:: xland(mgmxp,mgmyp)
  !real,parameter:: train= 1.0

  do K=1,MKX
     do I=ISTART,IEND
        outt(i,k) = 0.
        outq(i,k) = 0.
       outqc(i,k) = 0.
     enddo
  enddo

  do I=ISTART,IEND
     pre(i)  =0.
     xmb(i)  =0.
     xfac1(i)=1.
     xfac_for_dn(i) = 1.
  enddo
  !
  !--- Calculate ensemble average mass fluxes
  !
  if(i_use_stat_prop == 1                  ) then

     if(.not. cup_output_vars_alloc)  &
          call alloc_cup_output_vars(mgmxp,maxens,maxens3)


     !1st time: Updraft mass fluxes statistics properties
     !          sgrell1_3d and sgrell2_3d arrays will store them for output
     call massflx_stats_catt( &
          xf_ens , &	      !01
          mix    , &	    !02
          mgmxp  , &	    !03
          mjx    , &	    !04
          mgmyp  , &	    !05
          ensdim , &	    !06
          maxens , &	    !07
          maxens2, & 	     !08
          maxens3, &	     !09
          xmb_ave, &	     !10
          !
          xmb_std, &	     !11
          xmb_cur, &	     !12
          xmb_ske, &	     !13
          x_ave  , &	   !14
          x_std  , &	   !15
          x_ske  , &	   !16
          x_cur  , &	   !17
          x_ave_cap, &	     !18
          j,       &	     !19
          istart,  &	     !20
          !
          iend,    &	     !21
          ierr,    &	     !22
          m1,      & 	     !23
          m2,      & 	     !24
          m3,      & 	     !25
!          sgrell1_3d, &	     !26
!	  sgrell2_3d, &
	  1,	      &
	  x_ave_cap1, &
	  x_ave_cap2, &
	  x_ave_cap3, &
	  x_ave_cap4, &
	  x_ave_cap5,  &
	  dti)

     !2nd time: precipitation statistics properties
     !sgrell1_3d array will store them for output
     call massflx_stats_catt( &
          pr_ens, &          !01
          mix, &	     !02
          mgmxp, &	     !03
          mjx, &	     !04
          mgmyp, &	     !05
          ensdim, &	     !06
          maxens, &	     !07
          maxens2, & 	     !08
          maxens3, &	     !09
          pr_ave, &	     !10
          !
          pr_std, &	     !11
          pr_cur, &	     !12
          pr_ske, &	     !13
          x_ave, &	     !14
          x_std, &	     !15
          x_ske, &	     !16
          x_cur, &	     !17
          x_ave_cap, &	     !18
          j, &		     !19
          istart, &	     !20
          !
          iend, &	     !21
          ierr, &	     !22
          m1, & 	     !23
          m2, & 	     !24
          m3, & 	     !25
!          sgrell1_3d, &	     !26
!	  sgrell2_3d, &
	  2,	      &
	  x_ave_cap1, &
	  x_ave_cap2, &
	  x_ave_cap3, &
	  x_ave_cap4, &
	  x_ave_cap5, &
	  dti)

  ENDIF
  IF(calc_xmb_with_statis_prop == 1 .and. i_use_stat_prop == 1) then

     do i=istart,iend

        if(ierr(i).eq.0)then

           ncount = 0
           do n=(iens-1)*maxens2*maxens*maxens3+1, &
                    iens*maxens2*maxens*maxens3
                pr_ens(i,j,n) =   pr_ens(i,j,n)*xf_ens(i,j,n)
              outt_ens(i,j,n) = outt_ens(i,j,n)*xf_ens(i,j,n)

              xmb(i) = xmb(i) + xf_ens(i,j,n)
              ncount = ncount + 1
           enddo
           xfac_for_dn(i) = xmb(i)/float(ncount)

           !tunning process
           xmb(i) = xmb_ave(i) - tunning * xmb_std(i)
           if(xmb(i).lt.0.) xmb(i)= .1 * xmb_ave(i)

          ! --- Now use proper count of how many closures were actually
          ! used in cup_forcing_ens (including screening of some
          ! closures over water) to properly normalize xmb
           clos_wei=16./max(1.,closure_n(i))
           !if (xland(i,j).gt.0.98)xmb(i)=xmb(i)*clos_wei
	   xmb(i)=xmb(i)*clos_wei

           xfac1(i) = xmb(i)
           !srf - inconsistency at downdraft mass flux due tunning process
           !      on xmb - correction factor:
           xfac_for_dn(i) = xmb(i)/(xfac_for_dn(i) + 1.e-16)

           if(xmb_ave(i) .lt. 1.e-10) ierr(i) = 11

        endif
     enddo

  else

     !----  Simple average
     do I=ISTART,IEND
        ncount=0
        xmb(i)=0.

        if(ierr(i).eq.0)then

           do n=(iens-1)*maxens2*maxens*maxens3+1,&
	         iens   *maxens2*maxens*maxens3
                pr_ens(i,j,n) =   pr_ens(i,j,n)*xf_ens(i,j,n)
              outt_ens(i,j,n) = outt_ens(i,j,n)*xf_ens(i,j,n)
              !
              !srf- esta restricao gera incompatibilidade com o calculo de
              !     do fluxo de massa do  downdraft calculado fora
              !IF(xf_ens(i,j,n).GT.0.)THEN !Lufla
              if(xf_ens(i,j,n).ge.0.)then
                 xmb(i) = xmb(i) + xf_ens(i,j,n)
                 ncount = ncount + 1
                 !if(i.eq.ipr.and.j.eq.jpr) print *,'XF =',n,xf_ens(i,j,n)
              endif
           enddo

           if(ncount.gt.0)then
              xmb(i)=xmb(i)/float(ncount)
           else
              xmb(i)=0.
              ierr(i)=13
           endif

          ! --- Now use proper count of how many closures were actually
          ! used in cup_forcing_ens (including screening of some
          ! closures over water) to properly normalize xmb
           clos_wei=16./max(1.,closure_n(i))
           !if (xland(i,j).gt.0.98)xmb(i)=xmb(i)*clos_wei
	   xmb(i)=xmb(i)*clos_wei
	  !rmf xmb(i)=xmb(i)*train(i,j)
!if(i==19 .and. j==27)  write(0,*)'out1',xmb(i),closure_n(i)


        endif
        xfac1(i)=xmb(i)
     enddo
  end if
  !
  !-- Now do feedback
  !
  ddtes=250.
  !     if(name.eq.'shal')ddtes=500.
!-srf from WRF
! ddtes=500.
!     if(name.eq.'shal')ddtes=200.
!
  do k=1,mkx
     do i=istart,iend
        dtt  =0.
        dtq  =0.
        dtqc =0.
        dtpw =0.

        if(ierr(i).eq.0.and.k.le.ktop(i))then
           do n=1,maxens2
              dtt  = dtt  +  dellat(i,k,n)
              dtq  = dtq  +  dellaq(i,k,n)
              dtqc = dtqc + dellaqc(i,k,n)
              dtpw = dtpw + pwo_ens(i,k,n)

           enddo

           outtes = dtt*XMB(I)*86400./float(maxens2)

           if (outtes .gt. 2.*ddtes .and. k.gt.2) then
              XMB(I) = 2.*ddtes/outtes * xmb(i)
              outtes = 1.*ddtes
           endif

           if (outtes .lt. -ddtes)                then
              XMB(I) = -ddtes/outtes * xmb(i)
              outtes = -ddtes
           endif

           if (outtes .gt. .5*ddtes .and. k.le.2) then
              XMB(I) =    ddtes/outtes * xmb(i)
              outtes = .5*ddtes
           endif

            outt(i,k) = outt(i,k) +xmb(i)*dtt /float(maxens2)
            outq(i,k) = outq(i,k) +xmb(i)*dtq /float(maxens2)
           outqc(i,k) =outqc(i,k) +xmb(i)*dtqc/float(maxens2)
           pre(i)      = pre(i)   +xmb(i)*dtpw/float(maxens2) ! unit : kg[liq water]/(m^2 s)
	  ! if(i==19 .and. j==27)  write(0,*)'out2', xmb(i),dtpw,OUTQ(I,K),pre(i)

           !srf-----print-------
           !          if(j.eq.jpr   .and. i.eq.ipr) then
           !	   if(k.eq.1) then
           !            write(6,'(a1,78a1)') ' ',('-',m=1,78)
           !	    print*,'maxens2 j i k outtP  OUTQ OUTQC PREC dtpw xmb pwo_ens'
           !	   endif
           !          write(6,'(4i4,8e12.4)') maxens2,j,i,k
           !     &   ,86400.*outt(I,K),1000.*86400.*OUTQ(I,K)
           !     &   ,1000.*86400.*OUTQC(I,K)
           !     &   ,pre(i),dtpw,xmb(i),pwo_ens(i,k,n)
           !
           !
           !	   if(k.eq.ktop(i)) then
           !            write(6,'(a1,78a1)') ' ',('-',m=1,78)
           !           endif
           !           if(k.eq.mkx) write(6,'(a1,78a1)') ' ',('-',m=1,78)
           !          endif
           !srf-----print-------

        endif
     enddo
  enddo


!srf-----print-------
!     if(i==47 .and. j==21) then
!	    dtpw=0.
!	   dtt=0.
!	    do n=(iens-1)*maxens2*maxens*maxens3+1,&
!		 iens	*maxens2*maxens*maxens3
!		 !if(pr_ens(i,j,n) == 0.)
!		 !print*,'pr_ens n',pr_ens(i,j,n)*3600. ,n
!		  dtpw=dtpw+pr_ens(i,j,n)/float(maxens2*maxens*maxens3)
!
!!		 dtt =dtt+xf_ens(i,j,n)/float(maxens2*maxens*maxens3)
!	    enddo
!	 print*,'==========================================================='
!	 print*,'fluxo de massa e prec final =',xmb(i),pre(i)*3600.
!	 print*,'fluxo de massa e prec final=', dtt,dtpw*3600.
!	    do n=1,maxens2
!	      dtq  =0.
!	      do k=1,ktop(i)
!		 dtq = dtq + pwo_ens(i,k,n)
!	      enddo
!	      print*,'prec maxens2 n=',n,dtq*3600.
!	    enddo
!
!
!	print*,'==========================================================='
!
!	!stop 43434
!     endif
!     endif
!  enddo
!!srf-----print-------
!


  do I=ISTART,IEND
     if(ierr(i).eq.0)then
        xfac1(i)=xmb(i)/(xfac1(i)+1.e-16)
        do k=(iens-1)*maxens2*maxens*maxens3+1,iens*maxens2*maxens*maxens3
            massfln(i,j,k)=  massfln(i,j,k)*xfac1(i)*xfac_for_dn(i)
             pr_ens(i,j,k)=   pr_ens(i,j,k)*xfac1(i)
           outt_ens(i,j,k)= outt_ens(i,j,k)*xfac1(i)
        enddo
     endif
  enddo

end subroutine cup_output_ens_catt
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine massflx_stats_catt( &
     xf_ens, &       !01
     mix, &          !02
     mgmxp, &        !03
     mjx, &          !04
     mgmyp, &        !05
     ensdim, &       !06
     maxens, &       !07
     maxens2, &      !08
     maxens3, &      !09
    !
     xt_ave, &       !10
     xt_std, &       !11
     xt_cur, &	  !12
     xt_ske, &	  !13
     x_ave,  &	  !14
     x_std,  &	  !15
     x_ske,  &	  !16
     x_cur,  &	  !17
     x_ave_cap, &	  !18
     j, &		  !19
     istart, &	  !20
     !
     iend, &         !21
     ierr, &	  !22
     m1, & 	  !23
     m2, & 	  !24
     m3, & 	  !25
!     sgrell1_3d, &	  !26
!     sgrell2_3d, &
     itest,      &  !27
     x_ave_cap1, &
     x_ave_cap2, &
     x_ave_cap3, &
     x_ave_cap4, &
     x_ave_cap5, &
     dti)
!--(DMK-CCATT-INI)-----------------------------------------------------
  use ccatt_start, only:        &
       CCATT           ! intent(in)
!--(DMK-CCATT-OLD)-----------------------------------------------------
!  use catt_start, only:        &
!       CATT           ! intent(in)
!--(DMK-CCATT-FIM)-----------------------------------------------------

  implicit none

  integer :: mix,mgmxp,mjx,mgmyp,ensdim,maxens   !,mkx
  integer :: maxens2,maxens3,i,j,k,istart,iend,num,ielem,ncount
  integer :: kk,m1,m2,m3,num2,num3,iedt,itest
  integer :: ierr(mgmxp)
  real :: xf_ens(mgmxp,mgmyp,ensdim)
  real :: xt_ave(mgmxp),xt_std(mgmxp),xt_ske(mgmxp),xt_cur(mgmxp)
  real :: x_ave(mgmxp,maxens3),x_std(mgmxp,maxens3)
  real :: x_ske(mgmxp,maxens3),x_cur(mgmxp,maxens3)
  real :: x_ave_cap(mgmxp,maxens)

  real :: small_number,rxxx


!08-03-2004
!   extra arrays for outputting precip from cap ensembles in dependence of closures
!
  real :: x_ave_cap1(mgmxp,maxens), x_ave_cap2(mgmxp,maxens) &
         ,x_ave_cap3(mgmxp,maxens), x_ave_cap4(mgmxp,maxens) &
         ,x_ave_cap5(mgmxp,maxens)
  real dti,pw_cap

!- for output at analysis only
!  real, dimension(m1,m2,m3) :: sgrell1_3d,sgrell2_3d

  num =maxens *maxens2 ! old way: ensdim/maxens3
  num2=maxens2*maxens3 ! old way: ensdim/maxens
  num3=maxens2         ! old way: num2/maxens3
  small_number=1.e-6


  do k=1,maxens
     do i=istart,iend
        x_ave_cap (i,k)=0.
        x_ave_cap1(i,k)=0.
        x_ave_cap2(i,k)=0.
        x_ave_cap3(i,k)=0.
        x_ave_cap4(i,k)=0.
        x_ave_cap5(i,k)=0.
     enddo
  enddo

  do k=1,maxens3
     do i=istart,iend
        x_ave(i,k)=0.
        x_std(i,k)=0.
        x_ske(i,k)=0.
        x_cur(i,k)=0.
     enddo
  enddo

  do i=istart,iend
     xt_ave(i)=0.
     xt_std(i)=0.
     xt_ske(i)=0.
     xt_cur(i)=0.
  enddo


!!!!!!!!!!!!!!!!!!
!   if(itest == 1) xf_ens(:,j,:)=0.01
!   if(itest == 2) xf_ens(:,j,:)=0.001
!!!!!!!!!!!!!!!!!!



!srf-                The mean for each cap_max (1,2 and 3) :
!     average in maxen2 and maxens3 ensemble elements for each
!     maxens (cap_max) ensemble element
  do iedt=1,maxens2
     do   k=1,maxens
        do kk=1,maxens3
           do i=istart,iend
              if(ierr(i).eq.0) &
                 x_ave_cap(i,k) = x_ave_cap(i,k) + &
		                  xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+kk)



!===========print
!      if(itest==1) then
!      if(i==47 .and. j==21) then
!       if(iedt==1 .and. kk==1 .and. k==1)print*,'=== itest=',itest
!     	print*,iedt,k,kk,&
!	x_ave_cap(i,k)&
!	,xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+kk)
!!     	 enddo
!!	 print*,'PW MEAN CAP=',pw_cap/15.
!!       print*,'==========================================================='
!      endif
!      endif
!===========print


           enddo
        enddo
     enddo
  enddo

!srf-  The mean for each cap_max mean on the closure group (grell, AS, KF, etc)
!                                     and the precip efficiency
!      average in maxen2 and maxens3 ensemble elements for each
!      maxens (cap_max) ensemble element
! x_ave_cap1(k) = fluxo de massa media sobre os tres elementos
!                 de fechamento grell e sobre os tres elementos
!                 da eficiencia de precipitacao para cada
!                 cap_max (k=1,2,3)
  do iedt=1,maxens2
    do   k=1,maxens
      do i=istart,iend
       if(ierr(i).eq.0) then
          x_ave_cap1(i,k) = x_ave_cap1(i,k) +                        &
     	   (xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+1)+     &
     	    xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+2)+     &
     	    xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+3))*.33333

          x_ave_cap2(i,k) = x_ave_cap2(i,k) +                        &
     	    (xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+4)+    &
     	     xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+5)+    &
     	     xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+6))*.33333

          x_ave_cap3(i,k) = x_ave_cap3(i,k) +                       &
     	    (xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+7)+   &
     	     xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+8)+   &
     	     xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+9))*.33333

          x_ave_cap4(i,k) = x_ave_cap4(i,k) +                       &
     	    (xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+10)+  &
     	     xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+11)+  &
     	     xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+12))*.33333

          x_ave_cap5(i,k) = x_ave_cap5(i,k) +                        &
     	   (xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+13)+    &
!    	     xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+14)+
     	    xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+15)+    &
!    	     xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+16))*.25
     	    xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+16))*.33333


!===========print
!       if(itest==1) then
!       if(i==47 .and. j==21) then
!
!     	!if(iedt==1 .and. k==1)
!	print*,'======= 	FLUXO DE MASSA ',iedt,k
!          print*,x_ave_cap1(i,k),xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+1),     &
!				  xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+2),     &
!				  xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+3)
!	   print*,x_ave_cap2(i,k),xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+4),     &
!				  xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+5),     &
!				  xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+6)
!	   print*,x_ave_cap3(i,k),xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+7),     &
!				  xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+8),     &
!				  xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+9)
!	   print*,x_ave_cap4(i,k),xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+10),     &
!				  xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+11),     &
!				  xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+12)
!	   print*,x_ave_cap5(i,k),xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+13),     &
!				  xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+14),     &
!				  xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+15)
!
!	  print*,'=='
!
!      endif
!      endif
!===========print



       endif
    enddo
   enddo
  enddo

  do k=1,maxens
     do i=istart,iend
        if(ierr(i).eq.0) then
	   x_ave_cap (i,k)=x_ave_cap (i,k)/float(num2)
           x_ave_cap1(i,k)=x_ave_cap1(i,k)/float(num3)
           x_ave_cap2(i,k)=x_ave_cap2(i,k)/float(num3)
           x_ave_cap3(i,k)=x_ave_cap3(i,k)/float(num3)
           x_ave_cap4(i,k)=x_ave_cap4(i,k)/float(num3)
           x_ave_cap5(i,k)=x_ave_cap5(i,k)/float(num3)




!===========print
!	if(itest==1) then
!	if(i==47 .and. j==21) then
!
!	 !if(iedt==1 .and. k==1)
!	 print*,'=======	 FLUXO DE MASSA M ',k,num2,num3
!	   print*,x_ave_cap (i,k)
!	   print*,x_ave_cap1(i,k)
!
!
!	   print*,x_ave_cap2(i,k)
!
!
!	   print*,x_ave_cap3(i,k)
!
!
!	   print*,x_ave_cap4(i,k)
!
!
!	   print*,x_ave_cap5(i,k)
!
!
!
!	  print*,'=='
!
!      endif
!      endif
!===========print




	 endif
     enddo
  enddo


!srf- mass flux average for each closure:
!     average in maxens and maxens2 ensemble elements for each
!     maxens3 (closure) ensemble element
!
  do kk=1,num ! num = maxens * maxens2
     do k=1,maxens3
        do i=istart,iend
           if(ierr(i).eq.0)then
              x_ave(i,k)=x_ave(i,k) + xf_ens(i,j,maxens3*(kk-1)+k)
           endif
        enddo
     enddo
  enddo

  do k=1,maxens3
     do i=istart,iend
        if(ierr(i).eq.0) x_ave(i,k) = x_ave(i,k)/float(num)
     enddo
  enddo

!srf- total average in maxens, maxen2 and maxens3 ensemble elements
!new way
  do k=1,maxens3
     do i=istart,iend
        if(ierr(i).eq.0) xt_ave(i) = xt_ave(i) + x_ave(i,k)/float(maxens3)
     enddo
  enddo

!old way
!  do k=1,maxens3
!     do i=istart,iend
!        if(ierr(i).eq.0) xt_ave(i)=xt_ave(i)+x_ave(i,k)
!     enddo
!  enddo
!  do i=istart,iend
!     if(ierr(i).eq.0) xt_ave(i)=xt_ave(i)/float(maxens3)
!  enddo

!srf-----print-------
!  do I=ISTART,IEND
!     if(ierr(i).eq.0)then
!     if(i==47 .and. j==21) then!
!
!       print*,'==========================================================='
!       print*,'fluxo de massa 1=',xt_ave(i)
!       print*,'==========================================================='
!     endif
!     endif
!  enddo

!  do I=ISTART,IEND
!     if(ierr(i).eq.0)then
!     if(i==47 .and. j==21) then
!       xt_ave(i)=0.
!       ncount =0
!         do kk=1,num ! num = maxens * maxens2
!            do k=1,maxens3
!	      ncount =ncount +1
!               xt_ave(i) = xt_ave(i) + xf_ens(i,j,maxens3*(kk-1)+k)
!	       print*,k,kk,ncount,xf_ens(i,j,maxens3*(kk-1)+k),xt_ave(i)
!	     enddo
!          enddo
!
!	print*,'==========================================================='
!	print*,'fluxo de massa 2 =',xt_ave(i),ncount,xt_ave(i)/float(ncount)
!	print*,'==========================================================='
!	stop 4442
!     endif
!     endif
!  enddo
!     !srf-----print-------
!
!
!




  !--- now do std, skewness,curtosis
  !srf- stats properties in maxens and maxens2 ensemble elements for each
  !     maxens3 (closure) ensemble element (for each closure)
 do kk=1,num
     do k=1,maxens3
        do i=istart,iend
           if(ierr(i).eq.0.and.x_ave(i,k).gt.0.)then
!-----------------------    bug- underflow
              rxxx       = max( small_number, xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))
             !x_std(i,k) = x_std(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**2
              x_std(i,k) = x_std(i,k)+(rxxx)**2

             !x_ske(i,k) = x_ske(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**3
              x_ske(i,k) = x_ske(i,k)+(rxxx)**3

	     !x_cur(i,k) = x_cur(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**4
	      x_cur(i,k) = x_cur(i,k)+(rxxx)**4
!-orig
!             x_std(i,k) = x_std(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**2
!             x_ske(i,k) = x_ske(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**3
!	      x_cur(i,k) = x_cur(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**4
!-----------------------
              !
              xt_std(i)  = xt_std(i) +(xf_ens(i,j,maxens3*(kk-1)+k)-xt_ave(i) )**2

           endif
        enddo
     enddo
  enddo



  do k=1,maxens3
     do i=istart,iend
        if(ierr(i).eq.0.and.xt_ave(i).gt.0.)then
           !        xt_std(i) = xt_std(i)+(x_ave(i,k)-xt_ave(i))**2
           xt_ske(i) = xt_ske(i)+(x_ave(i,k)-xt_ave(i))**3

!bug- underflow
            if(x_ave(i,k)-xt_ave(i) > small_number) &
            xt_cur(i) = xt_cur(i)+(x_ave(i,k)-xt_ave(i))**4
!           xt_cur(i) = xt_cur(i)+(x_ave(i,k)-xt_ave(i))**4
        endif
     enddo
  enddo
  do k=1,maxens3
     do i=istart,iend
!srf-21/02/2005
!       if(ierr(i).eq.0.and.x_std(i,k).gt.0.)then
        if(ierr(i).eq.0.and.x_std(i,k).gt.small_number)then
           x_std(i,k) = x_std(i,k)/float(num)
           x_std(i,k) = sqrt(x_std(i,k))
           if(xt_std(i) .gt. 0. ) then
              x_ske(i,k) = x_ske(i,k)/float(num)/x_std(i,k)**3
              x_cur(i,k) = x_cur(i,k)/float(num)/x_std(i,k)**4

           endif
        endif
        !       print*,'                               '
        !       print*,'Some statistics at gridpoint i,j, ierr',i,j,ierr(i)
        !       print*,'statistics for closure number ',k
        !       print*,'Average= ',x_ave(i,k),'  Std= ',x_std(i,k)
        !       print*,'Skewness= ',x_ske(i,k),' Curtosis= ',x_cur(i,k)
        !       print*,'                               '
     enddo
  enddo

  do i=istart,iend
     if(ierr(i).eq.0.and.xt_ave(i).gt.0.)then
        !           xt_std(i)=xt_std(i)/float(maxens3)
        xt_std(i)=xt_std(i)/float(num*maxens3)
        xt_std(i)=sqrt(xt_std(i))
!srf-21/02/2005
!       if(xt_std(i) .gt. 0. ) then
        if(xt_std(i) .gt. small_number ) then
           xt_ske(i)=xt_ske(i)/float(maxens3)/xt_std(i)**3
           xt_cur(i)=xt_cur(i)/float(maxens3)/xt_std(i)**4
        endif
        !       print*,' 			      '
        !       print*,'Total ensemble independent statistics at i j =',i,j
        !       print*,'Average= ',xt_ave(i),'  Std= ',xt_std(i)
        !       print*,'Skewness= ',xt_ske(i),' Curtosis= ',xt_cur(i)
        !       print*,' 			      '
     endif
     !srf -- store at analysis files (for output)
     !-------------------For Deep Cumulus -----------------------
!   if (CATT==1) then
!     if(itest.eq.1)then
!        if(maxens3.ne.16 .or. m1 .lt. 30) then
!           print*,'**** maxens3 for deep different of 16 *****'
!           print*,'****               OR                 *****'
!           print*,'****       NZP   less than 30         *****'
!           stop
!        endif
!        !this the total mean mass flux (average over all ensemble)
!        sgrell1_3d(2,i,j) = xt_ave(i)
!        sgrell1_3d(3,i,j) = xt_std(i)
!        sgrell1_3d(4,i,j) = xt_ske(i)
!        sgrell1_3d(5,i,j) = xt_cur(i)

!        !this the mass flux for the closure type 1  - Grell
!        sgrell1_3d(6,i,j) = .333*(x_ave(i,1)+x_ave(i,2)+x_ave(i,3))
!        !this the mass flux for the closure type 4  - low level omega
!        sgrell1_3d(7,i,j) = .333*(x_ave(i,4)+x_ave(i,5)+x_ave(i,6))
!        !this the mass flux for the closure type 7  - moisture convergence
!        sgrell1_3d(8,i,j) = .333*(x_ave(i,7)+x_ave(i,8)+x_ave(i,9))
!        !this the mass flux for the closure type 10  - Stability closure (FC-KF)
!        sgrell1_3d(9,i,j) = .333*(x_ave(i,10)+x_ave(i,11)+x_ave(i,12))
!        !this the mass flux for the closure type 13  - Arakawa-Schubert
!        sgrell1_3d(10,i,j) = .25*(x_ave(i,13)+x_ave(i,14)+x_ave(i,15)+x_ave(i,16))
!
!        !this the average mass flux for each cap_max ensemble element
!        do k=1,maxens
!           sgrell1_3d(16+(k-1),i,j) = x_ave_cap(i,k)
!        enddo

!         sgrell2_3d(1,i,j)  =x_ave_cap1(i,1)
!         sgrell2_3d(2,i,j)  =x_ave_cap1(i,2)
!         sgrell2_3d(3,i,j)  =x_ave_cap1(i,3)
!         sgrell2_3d(4,i,j)  =x_ave_cap2(i,1)
!         sgrell2_3d(5,i,j)  =x_ave_cap2(i,2)
!         sgrell2_3d(6,i,j)  =x_ave_cap2(i,3)
!         sgrell2_3d(7,i,j)  =x_ave_cap3(i,1)
!         sgrell2_3d(8,i,j)  =x_ave_cap3(i,2)
!         sgrell2_3d(9,i,j)  =x_ave_cap3(i,3)
!         sgrell2_3d(10,i,j) =x_ave_cap4(i,1)
!         sgrell2_3d(11,i,j) =x_ave_cap4(i,2)
!         sgrell2_3d(12,i,j) =x_ave_cap4(i,3)
!         sgrell2_3d(13,i,j) =x_ave_cap5(i,1)
!         sgrell2_3d(14,i,j) =x_ave_cap5(i,2)
!         sgrell2_3d(15,i,j) =x_ave_cap5(i,3)

!===========print
!       if(i==47 .and. j==21) then
!
!     	print*,'==========================================================='
!     	print*,'======= 	FLUXO DE MASSA  	   ================'
!     	print*,'ens elem - inst  - mean : dti=',dti
!     	  do ielem=1,15
!     	    print*,ielem,sgrell2_3d(ielem ,i,j)
!     	 enddo
!	 !stop 4343
!       print*,'==========================================================='
!      endif
!===========print
!===========print
 !	if(i==47 .and. j==21) then
 !	 print*,'==========================================================='
 !	 print*,'=======  fluxo de massa medio =',xt_ave(i)
 !	 print*,'X_AVE_CAP=',x_ave_cap(i,1),x_ave_cap(i,2),x_ave_cap(i,3)
 !
 !	   do ielem=1,15,3
 !	     print*,'X_AVE_CAP=',ielem,sgrell2_3d(ielem,i,j  ) &
 !				      ,sgrell2_3d(ielem+1,i,j) &
 !				      ,sgrell2_3d(ielem+2,i,j)
 !        enddo
 !      print*,'==========================================================='
 !     endif
 !===========print
        !
!     elseif(itest.eq.2)then
!     sgrell1_3d(11,i,j) = .333*(x_ave(i,1)+x_ave(i,2)+x_ave(i,3))
!     sgrell1_3d(12,i,j) = .333*(x_ave(i,4)+x_ave(i,5)+x_ave(i,6))
!     sgrell1_3d(13,i,j) = .333*(x_ave(i,7)+x_ave(i,8)+x_ave(i,9))
!     sgrell1_3d(14,i,j) = .333*(x_ave(i,10)+x_ave(i,11)+x_ave(i,12))
!     sgrell1_3d(15,i,j) = .25*(x_ave(i,13)+x_ave(i,14)+x_ave(i,15)+x_ave(i,16))
!
!     sgrell1_3d(11,i,j) = .333*(x_ave(i,1)+x_ave(i,2)+x_ave(i,3)) &
!     			      *sgrell1_3d(6,i,j)*3600.
!     sgrell1_3d(12,i,j) = .333*(x_ave(i,4)+x_ave(i,5)+x_ave(i,6)) &
!     			      *sgrell1_3d(7,i,j)*3600.
!     sgrell1_3d(13,i,j) = .333*(x_ave(i,7)+x_ave(i,8)+x_ave(i,9)) &
!     			      *sgrell1_3d(8,i,j)*3600.
!     sgrell1_3d(14,i,j) = .333*(x_ave(i,10)+x_ave(i,11)+x_ave(i,12)) &
!     			      *sgrell1_3d(9,i,j)*3600.
!     sgrell1_3d(15,i,j) = .250*(x_ave(i,13)+x_ave(i,14)+x_ave(i,15)+x_ave(i,16)) &
!     			      *sgrell1_3d(10,i,j)*3600.
!
!!    convert from kg/m2/s to kg/m2/hour
!     sgrell2_3d(1,i,j)  =x_ave_cap1(i,1)*sgrell2_3d(1 ,i,j)*3600.
!     sgrell2_3d(2,i,j)  =x_ave_cap1(i,2)*sgrell2_3d(2 ,i,j)*3600.
!     sgrell2_3d(3,i,j)  =x_ave_cap1(i,3)*sgrell2_3d(3 ,i,j)*3600.
!     sgrell2_3d(4,i,j)  =x_ave_cap2(i,1)*sgrell2_3d(4 ,i,j)*3600.
!     sgrell2_3d(5,i,j)  =x_ave_cap2(i,2)*sgrell2_3d(5 ,i,j)*3600.
!     sgrell2_3d(6,i,j)  =x_ave_cap2(i,3)*sgrell2_3d(6 ,i,j)*3600.
!     sgrell2_3d(7,i,j)  =x_ave_cap3(i,1)*sgrell2_3d(7 ,i,j)*3600.
!     sgrell2_3d(8,i,j)  =x_ave_cap3(i,2)*sgrell2_3d(8 ,i,j)*3600.
!     sgrell2_3d(9,i,j)  =x_ave_cap3(i,3)*sgrell2_3d(9 ,i,j)*3600.
!     sgrell2_3d(10,i,j) =x_ave_cap4(i,1)*sgrell2_3d(10,i,j)*3600.
!     sgrell2_3d(11,i,j) =x_ave_cap4(i,2)*sgrell2_3d(11,i,j)*3600.
!     sgrell2_3d(12,i,j) =x_ave_cap4(i,3)*sgrell2_3d(12,i,j)*3600.
!     sgrell2_3d(13,i,j) =x_ave_cap5(i,1)*sgrell2_3d(13,i,j)*3600.
!     sgrell2_3d(14,i,j) =x_ave_cap5(i,2)*sgrell2_3d(14,i,j)*3600.
!     sgrell2_3d(15,i,j) =x_ave_cap5(i,3)*sgrell2_3d(15,i,j)*3600.
!
!!-time average (hardwire to franl=10800 confrq=600)
!!         dti = = franl/confrq
!!     dti=0.05555 (for franl=10800 confrq=600)
!
!!- precipitation closure 1 (grell, mm/h)
!     sgrell2_3d(16,i,j)  =sgrell2_3d(16,i,j)+ sgrell2_3d(1 ,i,j)*dti ! large cap
!     sgrell2_3d(17,i,j)  =sgrell2_3d(17,i,j)+ sgrell2_3d(2 ,i,j)*dti ! medium cap
!     sgrell2_3d(18,i,j)  =sgrell2_3d(18,i,j)+ sgrell2_3d(3 ,i,j)*dti ! low cap'
!
!!- precipitation closure 2 (omega, mm/h)
!     sgrell2_3d(19,i,j)  =sgrell2_3d(19,i,j)+ sgrell2_3d(4 ,i,j)*dti ! large cap
!     sgrell2_3d(20,i,j)  =sgrell2_3d(20,i,j)+ sgrell2_3d(5 ,i,j)*dti ! medium cap
!     sgrell2_3d(21,i,j)  =sgrell2_3d(21,i,j)+ sgrell2_3d(6 ,i,j)*dti ! low cap'
!
!!- precipitation closure 3 (Kuo, mm/h)
!     sgrell2_3d(22,i,j)  =sgrell2_3d(22,i,j)+ sgrell2_3d(7 ,i,j)*dti ! large cap
!     sgrell2_3d(23,i,j)  =sgrell2_3d(23,i,j)+ sgrell2_3d(8 ,i,j)*dti ! medium cap
!     sgrell2_3d(24,i,j)  =sgrell2_3d(24,i,j)+ sgrell2_3d(9 ,i,j)*dti ! low cap'
!
!!- precipitation closure 4 (K&F, mm/h)
!     sgrell2_3d(25,i,j)  =sgrell2_3d(25,i,j)+ sgrell2_3d(10,i,j)*dti ! large cap
!     sgrell2_3d(26,i,j)  =sgrell2_3d(26,i,j)+ sgrell2_3d(11,i,j)*dti ! medium cap
!     sgrell2_3d(27,i,j)  =sgrell2_3d(27,i,j)+ sgrell2_3d(12,i,j)*dti ! low cap'
!
!!- precipitation closure 5 (A&S, mm/h)
!     sgrell2_3d(28,i,j)  =sgrell2_3d(28,i,j)+ sgrell2_3d(13,i,j)*dti ! large cap
!     sgrell2_3d(29,i,j)  =sgrell2_3d(29,i,j)+ sgrell2_3d(14,i,j)*dti ! medium cap
!     sgrell2_3d(30,i,j)  =sgrell2_3d(30,i,j)+ sgrell2_3d(15,i,j)*dti ! low cap'

!===========print
!       if(i==47 .and. j==21) then
!        pw_cap = 0.
!     	print*,'==========================================================='
!     	print*,'======= 	PRECIPITACAO		   ================'
!     	print*,'ens elem - inst  - mean : dti=',dti
!     	  do ielem=1,15
!	    pw_cap=pw_cap+sgrell2_3d(ielem ,i,j)
!     	    print*,ielem,sgrell2_3d(ielem ,i,j),sgrell2_3d(ielem+15 ,i,j)
!     	 enddo
!	 print*,'PW MEAN CAP=',pw_cap/15.
!       print*,'==========================================================='
!      endif
!===========print
!        !---------------For Shallow Cumulus ---------------------------------
!    elseif(itest.eq.3)then
!        if(maxens3.ne.10) then
!           print*,'**** maxens3 for shallow different of 10 *****'
!           stop
!        endif
        !this the mass flux for the closure type 1  - Grell

!        sgrell1_3d(21,i,j) = .333*(x_ave(i,1)+x_ave(i,2)+x_ave(i,3))
!        !this the mass flux for the closure type 4 -  Arakawa-Schubert
!        sgrell1_3d(22,i,j) = .025*(x_ave(i,4)+x_ave(i,5)+x_ave(i,6)+x_ave(i,7) )
!        !this the mass flux for the closure type 7 - Stability closure (FC-KF)
!        sgrell1_3d(23,i,j) = .333*(x_ave(i,8)+x_ave(i,9)+x_ave(i,10))

        !       if(ierr(i).eq.0) then
        !	print*,'STATISTIC - SHALLOW XMB=',xt_ave(i)
        !	print*,'GR=',sgrell1_3d(21,i,j),'AS=',sgrell1_3d(22,i,j)
        !     &        ,'SC=',sgrell1_3d(23,i,j)
        !	print*,'AS',x_ave(i,4),x_ave(i,5),x_ave(i,6),x_ave(i,7)
        !       endif

        !this the total mean mass flux (average ovell all ensemble)
        !       sgrell1_3d(24,i,j) = xt_ave(i)
        !       sgrell1_3d(25,i,j) = xt_std(i)
        !       sgrell1_3d(26,i,j) = xt_ske(i)
        !       sgrell1_3d(27,i,j) = xt_cur(i)

!     endif
!    endif !CATT
  enddo

end subroutine massflx_stats_catt

!--------------------------------------------------------------------------
subroutine cup_forcing_ens_1_catt(aa0,aa1,xaa0,mbdt,dtime,xmb,ierr, &
     mix,mgmxp,mjx,mgmyp,mkx,mgmzp,istart,iend,xf,j,fquasi, &
     fstab,name,xland,maxens,iens,iedt,maxens2,ipr,jpr,maxens3, &
     mconv,omeg,zd,k22,zu,pr_ens,edt,aad,kbcon,massflx, &
     iact_old_gr,dir,ensdim,massfln,xff_ens3, xk,ierr2,ierr3,i_cup_dir_flag)
  implicit none
  character *(*) name

  integer i,istart,iend,mix,mgmxp,mjx,mgmyp,mkx,mgmzp,j,maxens,maxens3
  integer ensdim,ierr(mgmxp),iens,nall,iedt,maxens2,ipr,jpr,i_cup_dir_flag
  integer ierr2(mgmxp),ierr3(mgmxp)
  integer k22(mgmxp),kbcon(mgmxp),iact_old_gr(mgmxp,mgmyp)
  real aa0(mgmxp),aa1(mgmxp),xaa0(mgmxp,maxens),xmb(mgmxp)
  real mbdt(maxens),dtime,edt(mgmxp),aad(mgmxp),dir(mgmxp)  !,dxxf
  real xf(mgmxp,mgmyp,ensdim),xland(mgmxp,mgmyp)
  real pr_ens(mgmxp,mgmyp,ensdim)
  real mconv(mgmxp),omeg(mgmxp,mgmzp),zd(mgmxp,mgmzp),zu(mgmxp,mgmzp)
  real xff_ens3(maxens3),xk(maxens),xff0      !,xff1,xff2,xff3,xff
  real massflx(mgmxp,mgmyp)
  real massfln(mgmxp,mgmyp,ensdim)
  real massfld
  integer fquasi,fstab,nens,ne,n,nens3,iresult,iresultd,iresulte
  nens=0

  !  if(j.eq.jpr)then
  !     print *,'!!!!!!! IN CUP_FORCING !!!!!!!!!!!'
  !     print *,massflx(ipr,jpr)
  !  endif

  !--- LARGE SCALE FORCING
  !
  !  DO 100 I=ISTART,IEND
  do i=istart,iend
     xmb(i)=0.
     if(name.eq.'deeps'.and.ierr(i).gt.995)then
        aa0(i)=0.
        ierr(i)=0
     endif
     if(ierr(i).eq.0)then
        !--- treatment different for this closure
        if(name.eq.'deeps')then
           xff0= (AA1(I)-AA0(I))/DTIME
           xff_ens3(1)=(AA1(I)-AA0(I))/dtime
           !--- more original Arakawa-Schubert (climatologic value of aa0)
           !
           !     xff_ens3(2)=max(0.,(AA1(I)-1500.)/dtime)
           !     xff_ens3(2)=xff_ens3(1)
           !
           !--- omeg is in bar/s, mconv done with omeg in Pa/s
           !     more like Brown (1979), or Frank-Cohen (199?)
           !     xff_ens3(3)=-omeg(i,k22(i))/9.81
           !     xff_ens3(4)=-omeg(i,kbcon(i))/9.81
           !--- more like Krishnamurti et al.
           !     xff_ens3(5)=mconv(i)
           !     xff_ens3(6)=mconv(i)
           !--- more like Fritsch Chappel or Kain Fritsch (plus triggers)
           !     xff_ens3(7)=AA1(I)/(60.*30.)
           !     xff_ens3(8)=AA1(I)/(60.*40.)
           do nens=1,maxens
              XK(nens)=(XAA0(I,1)-AA1(I))/MBDT(2)
              !print*,nens,XK(nens),XAA0(I,nens)
              if(xk(nens).le.0.and.xk(nens).gt.-1.e-9) xk(nens)=-1.e-9
              if(xk(nens).gt.0.and.xk(nens).lt.1.e-9) xk(nens)=1.e-9
           enddo
           !--- add up all ensembles
           !DO 350 ne=1,maxens
           do ne=1,maxens
              !--- for every xk, we have maxens3 xffs
              !--- iens is from outermost ensemble (most expensive!
              !--- iedt (maxens2 belongs to it)
              !--- is from second, next outermost, not so expensive
              !--- so, for every outermost loop, we have maxens*maxens2*3
              !--- ensembles!!! nall would be 0, if everything is on first
              !--- loop index, then ne would start counting, then iedt, then iens....
              iresultd=0
              iresulte=0
              iresulte=1
              nall= (iens-1)*maxens3*maxens*maxens2+ &
	           (iedt-1)*maxens3*maxens+ &
                   (  ne-1)*maxens3
              !--- check for upwind convectionc
              if(i_cup_dir_flag == 1 ) then
               iresult=0
               massfld=0.
               call cup_direction2(i,j,dir,iact_old_gr,mix,mjx, &
                    mgmxp,mgmyp,massflx,iresult,ensdim,1,nall, &
                    maxens3,massfld)

               if(XK(ne).lt.0.and.xff0.gt.0.) iresultd=1
               iresulte=max(iresult,iresultd)
	      endif

	      if(iresulte.eq.1)then
                 !--- special treatment for stability closures
                 if(xff0.gt.0.)then
                    xf(i,j,nall+1)=max(0.,-xff_ens3(1)/xk(ne))+massfld
                 else
                    xf(i,j,nall+1)=massfld
                 endif
                 !--- store new for next time step
                 do nens3=1,maxens3
                    massfln(i,j,nall+nens3)=edt(i)*xf(i,j,nall+nens3)
                    massfln(i,j,nall+nens3)=max(0.,massfln(i,j,nall+nens3))

                 enddo
              endif
           end do
           !350       continue
           !go to 100
           cycle
        end if
     elseif(ierr(i).ne.20.and.ierr(i).ne.0)then
        do n=1,ensdim
           xf(i,j,n)=0.
           massfln(i,j,n)=0.
        enddo
     end if
  end do
  !100   continue
end subroutine cup_forcing_ens_1_catt

!--------------------------------------------------------------------------

subroutine cup_dellas_catt(ierr, z_cup, p_cup, hcd, edt, zd, cdd, he, mix,   &
     mgmxp, mkx, mgmzp, istart, iend, della, itest, j, mentrd_rate, zu, &
     cd, hc, ktop, k22, kbcon, mentr_rate2, jmin, he_cup, kdet, kpbl,    &
     ipr, jpr, name)
  use Phys_const, only: g
  implicit none
  character (LEN=*) name  !CHARACTER *(*) name
  integer mix, mgmxp, mkx, mgmzp, i, k, istart, iend,  &
       itest, j
  real z_cup(mgmxp,mgmzp), p_cup(mgmxp,mgmzp), hcd(mgmxp,mgmzp),  &
       zd(mgmxp,mgmzp), cdd(mgmxp,mgmzp), he(mgmxp,mgmzp),        &
       della(mgmxp,mgmzp), hc(mgmxp,mgmzp), cd(mgmxp,mgmzp),      &
       zu(mgmxp,mgmzp), he_cup(mgmxp,mgmzp)
  real edt(mgmxp)
  integer kbcon(mgmxp), ktop(mgmxp), k22(mgmxp), jmin(mgmxp)
  integer ierr(mgmxp), kdet(mgmxp), kpbl(mgmxp), ipr, jpr
  real entdo, dp, dz, mentrd_rate,                     &  !detdo1, detdo2,
       mentr_rate, subin, detdo, entup, detup,         &
       subdown, entdoj, entupk, detupk, totmas
  real xsum,xsumt,hesum

  real mentr_rate2(mgmxp,mgmzp)


  do k=2,mkx
     do I=ISTART,IEND
        della(i,k) = 0.
     enddo
  enddo
  !  xsum=0.
  !  hesum=0.
  !  xsumt=0.

  do K=2,MKX-1
     !DO 100 K=2,MKX-1
     do I=ISTART,IEND
        !DO 100 I=ISTART,IEND
        !IF (ierr(i).NE.0) GO TO 100
        if (ierr(i).ne.0) cycle
        !IF (K.GT.KTOP(I)) GO TO 100
        if (K.gt.KTOP(I)) cycle
        !
        !--- Specify detrainment of downdraft, has to be consistent
        !--- with zd calculations in soundd.
        !
        dz    = Z_cup(I,K+1)-Z_cup(I,K)
        detdo = edt(i)*CDD(i,K)   *dz*zd(i,k+1)
        entdo = edt(i)*mentrd_rate*dz*zd(i,k+1)
        subin = zu(i,k+1)-zd(i,k+1)*edt(i)
        entup = 0.
        detup = 0.
        if (k.ge.kbcon(i).and.k.lt.ktop(i)) then
           entup = mentr_rate2(i,k)*dz*zu(i,k)
           detup = CD(i,K+1) *dz*zu(i,k)
        endif
        subdown = ( zu(i,k)-zd(i,k)*edt(i) )
        entdoj  = 0.
        entupk  = 0.
        detupk  = 0.

        if (k.eq.jmin(i)) then
           entdoj  = zd(i,k)*edt(i)
        endif

        if (k.eq.k22(i)-1) then
!       if (k.eq.kpbl(i) ) then
           entupk  = zu(i,kpbl(i))
        endif

        if (k.gt.kdet(i)) then
           detdo   = 0.
        endif

        if (k.eq.ktop(i)-0) then
           detupk  = zu(i,ktop(i))
           subin   = 0.
        endif
        if (k.lt.kbcon(i)) then
           detup   = 0.
        endif
        !
        !--- Changed due to subsidence and entrainment
        !
        totmas =subin-subdown+detup-entup-entdo + detdo-entupk-entdoj+detupk
        if (abs(totmas).gt.1.e-6) then
           print *, '**totmass deep********', i, j, k, totmas, name
           !print *, kpbl(i), k22(i), kbcon(i), ktop(i)
           !          print *,'updr stuff = ',subin,
           !    1      subdown,detup,entup,entupk,detupk
           !          print *,'dddr stuff = ',entdo,
           !    1      detdo,entdoj
           stop
        endif

        !srf         dp =  100.*( p_cup(i,k-1)-p_cup(i,k) )
        dp =  100.*( p_cup(i,k)-p_cup(i,k+1) )
!        della(i,k)=(subin  *he_cup(i,k+1) - subdown*he_cup(i,k  ) +         &
!             	    detup*.5*( HC(i,K+1)+ HC(i,K)) +			    &
!             	    detdo*.5*(HCD(i,K+1)+HCD(i,K)) -			    &
!             	    entup*he(i,k) - entdo*he(i,k) -			    &
!             	    entupk*he_cup(i,k22(i)) -  entdoj*he_cup(i,jmin(i)) +   &
!             	    detupk*hc(i,ktop(i))                                    &
!		    )*g/dp
         della(i,k)=(                                 &
     		    subin  *he_cup(i,k+1)	      &
     		  - subdown*he_cup(i,k  )	      &
     		  + detup*.5*( HC(i,K+1)+ HC(i,K))    &
     		  + detdo*.5*(HCD(i,K+1)+HCD(i,K))    &
     		  - entup*he(i,k)		      &
     		  - entdo*he(i,k)		      &
     		  - entupk*he_cup(i,k22(i))	      &
     		  - entdoj*he_cup(i,jmin(i))	      &
     		  + detupk*hc(i,ktop(i))	      &
     		   )*g/dp
!      if(i.eq.ipr.and.j.eq.jpr)then
!         write(0,*)'d',k,della(i,k),(zu(i,k+1)*he_cup(i,k+1) &
!                    -zu(i,k)*he_cup(i,k))*g/dp,subin,subdown
!        write(0,*)'d',detup,entup,entdo,entupk,entdoj
!        print *,k,della(i,k),subin*he_cup(i,k+1),subdown*he_cup(i,k),
!     1            detdo*.5*(HCD(i,K+1)+HCD(i,K))
!        print *,k,detup*.5*(HC(i,K+1)+HC(i,K)),detupk*hc(i,ktop(i)),
!     1         entup*he(i,k),entdo*he(i,k)
!        print *,k,he_cup(i,k+1),he_cup(i,k),entupk*he_cup(i,k)
!       endif

        !     if(j.eq.3.and.i.eq.120)xsumt=xsumt+totmas
        !     if(j.eq.3.and.i.eq.120)hesum=hesum+he(i,k)*dp
        !     if(j.eq.3.and.i.eq.120)xsum=xsum+della(i,k)*dp
        !     if(j.eq.3.and.i.eq.120)print *,'xsum = ',xsum
        !IF (i.EQ.ipr.AND.j.EQ.jpr) THEN
        !   PRINT *, k, della(i,k), subin*he_cup(i,k+1),                  &
        !        subdown*he_cup(i,k), detdo*.5*(HCD(i,K+1)+HCD(i,K))
        !      PRINT *, k, detup*.5*(HC(i,K+1)+HC(i,K)),                  &
        !           detupk*hc(i,ktop(i)), entup*he(i,k),entdo*he(i,k)
        !      PRINT *, k, he_cup(i,k+1), he_cup(i,k), entupk*he_cup(i,k)
        !   ENDIF

        !100     CONTINUE
     enddo
  enddo
end subroutine cup_dellas_catt

!--------------------------------------------------------------------------

subroutine cup_dellabot_catt(ipr, jpr, he_cup, ierr, z_cup, p_cup, hcd, edt,   &
     zd, cdd, he, mix, mgmxp, mkx, mgmzp, istart, iend, della, itest, j,  &
     mentrd_rate, z)
  use Phys_const, only: g
  implicit none
  integer mix, mgmxp, mkx, mgmzp, i, istart, iend, itest, j
  real z_cup(mgmxp,mgmzp), p_cup(mgmxp,mgmzp), hcd(mgmxp,mgmzp),   &
       zd(mgmxp,mgmzp), cdd(mgmxp,mgmzp), he(mgmxp,mgmzp),         &
       della(mgmxp,mgmzp), he_cup(mgmxp,mgmzp), z(mgmxp,mgmzp), edt(mgmxp)
  integer ierr(mgmxp), ipr, jpr, m

!-local variables in this routine
  real detdo1, detdo2, entdo, dp, dz, mentrd_rate, subin, detdo

  do i=istart,iend
     !DO 100 i=istart,iend
     della(i,1)=0.
     !IF (ierr(i).NE.0) GO TO 100
     if (ierr(i).ne.0) cycle
     dz        =       z_cup(i,2)-z_cup(i,1)
     dp        = 100.*(p_cup(i,1)-p_cup(i,2))
     detdo1    = edt(i)*zd(i,2)*cdd(i,1)*dz
     detdo2    = edt(i)*zd(i,1)
     entdo     = edt(i)*zd(i,2)*mentrd_rate*dz
     subin     =-edt(I)*zd(i,2)
     detdo     = detdo1+detdo2-entdo+subin

     !check for mass conservation
     if(abs(detdo).gt.1.e-6)then
        print*,edt(i),zd(i,1),cdd(i,1),dz
        print*,detdo1,detdo2,-entdo,subin
        print*,'DELLABOT: totmas = ',detdo
        stop
     endif

     DELLA(I,1)= (detdo1*.5*(hcd(i,1)+hcd(i,2)) +  &
          detdo2*hcd(i,1) + subin*he_cup(i,2) -    &
          entdo*he(i,1))*g/dp

!      if(i.eq.ipr.and.j.eq.jpr)then
!       write(0,*)'db1',della(i,1),subin,entdo
!       write(0,*)'db2',detdo1,detdo2,detdo1+detdo2-entdo+subin
!      endif
!100  CONTINUE
  enddo
end subroutine cup_dellabot_catt

!------------------------------------------------------------

subroutine cup_forcing_ens_16_catt(aa0,aa1,xaa0,mbdt,dtime,xmb,ierr,   &
     mix,mgmxp,mjx,mgmyp,mkx,mgmzp,istart,iend,xf,j,fquasi,       &
     fstab,name,xland,maxens,iens,iedt,maxens2,ipr,jpr,maxens3,   &
     mconv,omeg,zd,k22,zu,pr_ens,edt,aad,kbcon,massflx,		  &
     iact_old_gr,dir,ensdim,massfln,massfld,iresult,xff_ens3,xk,  &
     p_cup,ktop,icoic,ierr2,ierr3,ngrid, &
!     sgrell1_2d,
     m2,m3,tscl_KF,i_cup_dir_flag,&
     glatg,glong,closure_n,xcoast)
  !
  ! ierr error value, maybe modified in this routine
  ! pr_ens = precipitation ensemble
  ! xf     = mass flux ensembles
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

  use Phys_const, only: g
  implicit none
  character (LEN=*) name

  integer k,i,istart,iend,mix,mgmxp,mjx,mgmyp,mkx,mgmzp,j,        &
          maxens,maxens3,ngrid

  !------ ensemble 3 dimension = 16
  integer mkxcrt,kclim
  parameter (mkxcrt=15)
  real pcrit(mkxcrt),acrit(mkxcrt),acritt(mkxcrt),aclim1,         &
       aclim2,aclim3,aclim4
  data pcrit/850.,800.,750.,700.,650.,600.,550.,500.,450.,400.,   &
       350.,300.,250.,200.,150./
  data acrit/.0633,.0445,.0553,.0664,.075,.1082,.1521,.2216,      &
       .3151,.3677,.41,.5255,.7663,1.1686,1.6851/
  !  GDAS derived acrit
  data acritt/.203,.515,.521,.566,.625,.665,.659,.688,            &
       .743,.813,.886,.947,1.138,1.377,1.896/
  integer ktop(mgmxp)
  real p_cup(mgmxp,mgmzp)

  integer fquasi,fstab,nens,ne,n,nens3,iresult,iresultd,          &
          iresulte,icoic
  integer ensdim,iens,nall,iedt,maxens2,ipr,jpr,m2,m3
  integer k22(mgmxp),kbcon(mgmxp),ierr(mgmxp),ierr2(mgmxp),ierr3(mgmxp) !Lufla
  integer iact_old_gr(mgmxp,mgmyp)

  real aa0(mgmxp),aa1(mgmxp),xaa0(mgmxp,maxens),xmb(mgmxp)
  real mbdt(maxens),dtime,edt(mgmxp),aad(mgmxp),dir(mgmxp)
  real xf(mgmxp,mgmyp,ensdim),xland(mgmxp,mgmyp)
  real pr_ens(mgmxp,mgmyp,ensdim)
  real mconv(mgmxp),omeg(mgmxp,mgmzp),zd(mgmxp,mgmzp),            &
       zu(mgmxp,mgmzp)
  real xff_ens3(maxens3),xk(maxens),xff0   !,xff1,xff2,xff3,xff
  real massflx(mgmxp,mgmyp)
  real massfln(mgmxp,mgmyp,ensdim)
  real xomg,massfld,a1

  real xcoast(mgmxp,mgmzp)

!  real, dimension(m2,m3)    :: sgrell1_2d
  real tscl_KF
  integer i_cup_dir_flag

  !srf  training
  real train_geo
  real glatg(mgmxp),glong(mgmxp),closure_n(mgmxp)
  integer training
  training =0
  !---

  nens=0

  !--- LARGE SCALE FORCING
  !
  do I=ISTART,IEND
     xmb(i)=0.
     if(name.eq.'deeps'.and.ierr(i).gt.995)then
        aa0(i) =0.
        ierr(i)=0
     endif
     if(ierr(i).eq.0)then
        !Added for ensemble 3 with dimension = 16
        kclim=0
        do k=mkxcrt,1,-1
           if(p_cup(i,ktop(i)).lt.pcrit(k))then
              kclim=k
              GO TO 9
           endif
        enddo
        if(p_cup(i,ktop(i)).gt.pcrit(1))kclim=1
9       continue
         kclim=max(kclim,1)
         k=max(kclim-1,1)
         aclim1=acrit(kclim)*1.e3
         aclim2=acrit(k)*1.e3
         aclim3=acritt(kclim)*1.e3
         aclim4=acritt(k)*1.e3
!
        !
        !--- Treatment different for this closure
        !
        if(name.eq.'deeps')then
        !
        !---- Grell's closure
        !
           xff0       =  (AA1(I)-AA0(I))/dtime
           xff_ens3(1)=    xff0
           xff_ens3(2)= .9*xff_ens3(1)
           xff_ens3(3)=1.1*xff_ens3(1)

!----------
!-srf out/2004: only for xff0 output
!	   sgrell1_2d(i,j) = xff0
!----------

           !
           !     More like Brown (1979), or Frank-Cohen (199?)
           !
           !---  omeg is in bar/s, mconv done with omeg in Pa/s
           xff_ens3(4)=     -omeg(i,k22(i))/g
           xff_ens3(5)=     -omeg(i,kbcon(i))/g
           xff_ens3(6)=     -omeg(i,1)/g

           do k=2,kbcon(i)-1
              xomg = -omeg(i,k)/g
              if(xomg .gt. xff_ens3(6)) xff_ens3(6)=xomg
           enddo
           !
           !--- More like Krishnamurti et al.
           !
           xff_ens3(7)=    mconv(i)
           xff_ens3(8)= .9*mconv(i)
           xff_ens3(9)=1.1*mconv(i)
           !
           !--- More like Fritsch Chappel or Kain Fritsch (plus triggers)
           !srf - changed at dec/2002 - greater timescale instab. removal
           xff_ens3(10)=AA1(I)/(60.*20.)
           xff_ens3(11)=AA1(I)/(60.*30.)
           xff_ens3(12)=AA1(I)/(60.*40.)
           !srf: A. Betts suggestion
           !xff_ens3(10)=AA1(I)/(tscl_KF    )! tscl_KF is already in seconds
           !xff_ens3(11)=AA1(I)/(tscl_KF*0.9)
           !xff_ens3(12)=AA1(I)/(tscl_KF*0.8)
	   !
           !--- More original Arakawa-Schubert (climatologic value of aa0)
           !
           xff_ens3(13)=max(0.,(AA1(I)-aclim1)/dtime)
           !srf- dec/2002 - later xff_ens3(14) will be equal to xff_ens3(13)
           xff_ens3(14)=max(0.,(AA1(I)-aclim2)/dtime)
           xff_ens3(15)=max(0.,(AA1(I)-aclim3)/dtime)
           xff_ens3(16)=max(0.,(AA1(I)-aclim4)/dtime)
!if(i==19 .and. j==27) write(0,*) 'xff_ens3 =',xff_ens3

           ! observe that XK(nens) is the same for any nens.
           do nens=1,maxens
!srf          XK(nens)=(XAA0(I,nens)-AA1(I))/MBDT(nens)
              XK(nens)=(XAA0(I,nens)-AA1(I))/MBDT(2  )
              if(xk(nens).le.0.and.xk(nens).gt.-1.e-9) xk(nens)=-1.e-9
              if(xk(nens).gt.0.and.xk(nens).lt.+1.e-9) xk(nens)=+1.e-9
           enddo
           !
           !--- Add up all ensembles
           !
           do ne=1,maxens
              nall= (iens-1)*maxens3*maxens*maxens2  &
                   +(iedt-1)*maxens3*maxens          &
                   +(ne  -1)*maxens3

              !          observe the mass flux calculation:
              !--------------------------------------------------------------!
              !           ne   |     ierr     | mass flux		     !
              !           1    |     ierr =0  |  mf1 = xff_ens3 / xk (ne)    !
              !           1    |     ierr >0  |  mf1 =  0		     !
              !           2    |     ierr2=0  |  mf2 = mf1		     !
              !           2    |     ierr2>0  |  mf2 =  0		     !
              !           3    |     ierr3=0  |  mf3 = mf1		     !
              !           3    |     ierr3>0  |  mf3 =  0		     !
              !							             !
              ! xk(ne) is the same for any 'ne'.			     !
              !--------------------------------------------------------------!
              ! if ierr2 > 0 (convection was not permited for that cap_max)
              ! then equal to zero the mass flux for the second member of the ensemble (maxens)
              if(ne.eq.2 .and. ierr2(i).gt.0)then
                 do nens3=1,maxens3
                         xf(i,j,nall+nens3)=0.
                    massfln(i,j,nall+nens3)=0.
                 enddo
                 cycle
              endif
              ! if ierr3 > 0 (convection was not permited for that cap_max)
              ! then equal to zero the mass flux for the third member of the ensemble (maxens)
              if(ne.eq.3 .and. ierr3(i).gt.0)then
                 do nens3=1,maxens3
                         xf(i,j,nall+nens3)=0.
                    massfln(i,j,nall+nens3)=0.
                 enddo
                 cycle
              endif
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
              iresultd=0
              iresulte=0

 go to 505

!-------------------------------------------------------------
! over water, enfor!e small cap for some of the closures
!
              if(xland(i,j).gt.0.98)then
!                 if(ierr2(i).gt.0.or.ierr3(i).gt.0)then
!       - ierr2 - 75 mb cap thickness, ierr3 - 125 cap thickness
!
! - for larger cap, set Grell closure to zero
                      xff_ens3(1) =0.
                      massfln(i,j,nall+1)=0.
                      xff_ens3(2) =0.
                      massfln(i,j,nall+2)=0.
                      xff_ens3(3) =0.
                      massfln(i,j,nall+3)=0.
                      if(ne.eq.1) closure_n(i)=closure_n(i)-3.

!                      xff_ens3(7) =0.
!                      massfln(i,j,nall+7)=0.
!                      xff_ens3(8) =0.
!                      massfln(i,j,nall+8)=0.
!                      xff_ens3(9) =0.
!!!                     massfln(i,j,nall+9)=0.
!                      closure_n(i)=closure_n(i)-1.
!                endif
!
!   also take out some closures in general
!
!                      xff_ens3(4) =0.
!                      massfln(i,j,nall+4)=0.
!                      xff_ens3(5) =0.
!                      massfln(i,j,nall+5)=0.
!                      xff_ens3(6) =0.
!                      massfln(i,j,nall+6)=0.
!                      closure_n(i)=closure_n(i)-3.
!
! taking out KF closures:
!
!              endif
		     xff_ens3(10)=0.
		     massfln(i,j,nall+10)=0.
		     xff_ens3(11)=0.
		     massfln(i,j,nall+11)=0.
		     xff_ens3(12)=0.
		     massfln(i,j,nall+12)=0.
		     if(ne.eq.1)closure_n(i)=closure_n(i)-3
             endif
! end water treatment
!- srf
!- over the coast, take out some trash closures
!             if( xcoast(i,j) == 1.)then
!             if( xland(i,j) > 0.2 .and. xland(i,j) <= 0.98 )then
!
!                xff_ens3(1) =0.
!                massfln(i,j,nall+1)=0.
!		xff_ens3(10)=0.
!		massfln(i,j,nall+10)=0.
!                if(ne.eq.1) closure_n(i)=closure_n(i)-2.
!
!             endif
! taking out AS closures:
505 continue
		     xff_ens3(13)=0.
		     massfln(i,j,nall+13)=0.
		     xff_ens3(14)=0.
		     massfln(i,j,nall+14)=0.
		     xff_ens3(15)=0.
		     massfln(i,j,nall+15)=0.
		     xff_ens3(16)=0.
		     massfln(i,j,nall+16)=0.
		     if(ne.eq.1 .and. iedt==1)closure_n(i)=closure_n(i)-4
!
!
!
!505 continue
!-------------------------------------------------------------              !
!Testar no futuro : ver sigficado dos param. iresultd,iresulte,iresult
              !--- check for upwind convection
              if(i_cup_dir_flag == 1 ) then
	      !
                  iresult=0
                  massfld=0.
                  call cup_direction2(i,j,dir,iact_old_gr,mix,mjx,  &
                       mgmxp,mgmyp,massflx,iresult,ensdim,1,nall,   &
                       maxens3,massfld)
                  print*,'cup_dir: i j massfld= ',i,j,massfld
	          !if(i.eq.ipr.and.j.eq.jpr.and.iedt.eq.1.and.ne.eq.1)then
                  !   print *,massfld,ne,iedt,iens
                  !   print *,xk(ne),xff_ens3(1),xff_ens3(2),xff_ens3(3)
                  !endif
                  if(XK(ne).lt.0.and.xff0.gt.0.)iresultd=1
                  iresulte=max(iresult,iresultd)
	          iresulte=1
	      else
	          massfld=0.
	          iresulte=1
	      endif

	      if(iresulte.eq.1)then
                 !
                 !--- Special treatment for stability closures
                 !

                 if(xff0.gt.0.)then
                    xf(i,j,nall+1) =max(0., -xff_ens3( 1) /xk(ne))+massfld
                    xf(i,j,nall+2) =max(0., -xff_ens3( 2) /xk(ne))+massfld
                    xf(i,j,nall+3) =max(0., -xff_ens3( 3) /xk(ne))+massfld
                    xf(i,j,nall+13)=max(0., -xff_ens3(13) /xk(ne))+massfld
                    xf(i,j,nall+14)=max(0., -xff_ens3(14) /xk(ne))+massfld
                    xf(i,j,nall+15)=max(0., -xff_ens3(15) /xk(ne))+massfld
                    xf(i,j,nall+16)=max(0., -xff_ens3(16) /xk(ne))+massfld
                 else
                    xf(i,j,nall+1) =massfld
                    xf(i,j,nall+2) =massfld
                    xf(i,j,nall+3) =massfld
                    xf(i,j,nall+13)=massfld
                    xf(i,j,nall+14)=massfld
                    xf(i,j,nall+15)=massfld
                    xf(i,j,nall+16)=massfld
                 endif
                 !
                 !--- if iresult.eq.1, following independent of xff0
                 !
                 xf(i,j,nall+4)=max(0.,xff_ens3(4)+massfld)
                 xf(i,j,nall+5)=max(0.,xff_ens3(5)+massfld)
                 xf(i,j,nall+6)=max(0.,xff_ens3(6)+massfld)


                !============================================
                !for physical initialization================
                !xff_ens3(7) = pcp rate observed
                !icoic = 7
                !only for the initial time
                !check units for moist convergence
                !at first 20 minutes
                !============================================


                a1 = max(1.e-9,pr_ens(i,j,nall+7))
                xf(i,j,nall+7)=max(0.,xff_ens3(7)/a1)

                a1 = max(1.e-9,pr_ens(i,j,nall+8))
                xf(i,j,nall+8)=max(0.,xff_ens3(8)/a1)

                a1 = max(1.e-9,pr_ens(i,j,nall+9))
                xf(i,j,nall+9)=max(0.,xff_ens3(9)/a1)

!               xf(i,j,nall+7)=max(0.,xff_ens3(7)/pr_ens(i,j,nall+7))
!               xf(i,j,nall+8)=max(0.,xff_ens3(8)/pr_ens(i,j,nall+8))
!               xf(i,j,nall+9)=max(0.,xff_ens3(9)/pr_ens(i,j,nall+9))

!testar a opcao abaixo
!mar2004???      if(XK(ne).lt.0.and.xff0.gt.0.)then
                 if(XK(ne).lt.0.)then
                    xf(i,j,nall+10)=max(0.,-xff_ens3(10)/xk(ne))+massfld
                    xf(i,j,nall+11)=max(0.,-xff_ens3(11)/xk(ne))+massfld
                    xf(i,j,nall+12)=max(0.,-xff_ens3(12)/xk(ne))+massfld
                 else
                    xf(i,j,nall+10)=massfld
                    xf(i,j,nall+11)=massfld
                    xf(i,j,nall+12)=massfld
                 endif
                 !==============
                 !05-12-2002
                 !srf - forcing 14 is too bad, use the same for 13:
                 !A&S (14) = A&S (13)
                 xf(i,j,nall+14)=xf(i,j,nall+13)
                 !
                 !-----------------------------------------
		 !srf 02-jul 2006
		 ! training
                 if(training==1 .and. icoic == 0) then
		    ! 80% GRELL
                    xf(i,j,nall+4)=xf(i,j,nall+1)
                    xf(i,j,nall+5)=xf(i,j,nall+2)
                    xf(i,j,nall+6)=xf(i,j,nall+3)
                    xf(i,j,nall+7)=xf(i,j,nall+1)
                    xf(i,j,nall+8)=xf(i,j,nall+2)
                    xf(i,j,nall+9)=xf(i,j,nall+3)
                    !xf(i,j,nall+13) = xf(i,j,nall+10)
                    !xf(i,j,nall+14) = xf(i,j,nall+10)
                    !xf(i,j,nall+15) = xf(i,j,nall+11)
                    !xf(i,j,nall+16) = xf(i,j,nall+12)
		    xf(i,j,nall+13) = xf(i,j,nall+2)
		    xf(i,j,nall+14) = xf(i,j,nall+3)
		    xf(i,j,nall+15) = xf(i,j,nall+1)
		    xf(i,j,nall+16) = xf(i,j,nall+2)
		 endif
		 !---------------------------------------
                 !----- training geografico-------------
                 !if(training==1) then
		 !   train_geo=1.
		 !   if(glatg(i) .lt. -20.) train_geo=min(2.,1+0.2*(abs(glatg(i))-20.))
                 !   do nens3=1,maxens3
                 !      xf(i,j,nall+nens3)=train_geo*xf(i,j,nall+nens3)
                 !   enddo
                 !endif

                 !----- 1d closure ensemble -------------
                 if(icoic.ge.1)then
		    closure_n(i)=0.
                    do nens3=1,maxens3
                       xf(i,j,nall+nens3)=xf(i,j,nall+icoic)
                    enddo
                 endif
                 !
                 !--- store new for next time step
                 !
                 do nens3=1,maxens3
                    massfln(i,j,nall+nens3)=edt(i)*xf(i,j,nall+nens3)
                    massfln(i,j,nall+nens3)=max(0.,massfln(i,j,nall+nens3))
                 enddo
              endif
           enddo
           cycle !go to 100
        endif
     elseif(ierr(i).ne.20.and.ierr(i).ne.0)then
        do n=1,ensdim
           xf(i,j,n)=0.
           massfln(i,j,n)=0.
        enddo
     endif
     !100  CONTINUE
  enddo
  return
end subroutine cup_forcing_ens_16_catt
!--------------------------------------------------------------------
