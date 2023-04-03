!----------------------------------------------------------------------
subroutine cup_up_he(k22, hkb, z_cup, cd, entr2, he_cup, hc,               &
     mix, mgmxp, mkx, mgmzp, kbcon, ierr, istart, iend, dby, he, hes_cup)
  implicit none
  integer i, j, k, mix, mgmxp, mkx, mgmzp, istart, iend
  real he_cup(mgmxp,mgmzp), hc(mgmxp,mgmzp), hkb(mgmxp),       &
       z_cup(mgmxp,mgmzp), cd(mgmxp,mgmzp), dby(mgmxp,mgmzp),  &
       he(mgmxp,mgmzp), hes_cup(mgmxp,mgmzp)
       
  integer kbcon(mgmxp), ierr(mgmxp), k22(mgmxp)
  real entr, dz
  real entr2(mgmxp,mgmzp)
  
  
  ! hc = cloud moist static energy
  ! hkb = moist static energy at originating level
  ! he = moist static energy on model levels
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function 
  ! z_cup = heights of model cloud levels 
  ! entr = entrainment rate
  !--- Moist static energy inside cloud

  do i=istart,iend
     if (ierr(i).eq.0.) then
        hkb(i)=he_cup(i,k22(i))
        do k=1,k22(i)
           hc(i,k)=he_cup(i,k)
           DBY(I,K)=0.
        enddo
        do k=k22(i),kbcon(i)-1
           hc(i,k)=hkb(i)
           DBY(I,K)=0.
        enddo
        k=kbcon(i)
        hc(i,k)=hkb(i)
        DBY(I,Kbcon(i))=Hkb(I)-HES_cup(I,K)
     endif
  enddo
  do k=2,mkx-1
     do i=istart,iend
        if (k.gt.kbcon(i).and.ierr(i).eq.0.) then
           DZ=Z_cup(i,K)-Z_cup(i,K-1)
           HC(i,K)=(HC(i,K-1)*(1.-.5*CD(i,K)*DZ)+entr2(i,k)*      &
                DZ*HE(i,K-1))/(1.+entr2(i,k)*DZ-.5*cd(i,k)*dz)
           DBY(I,K)=HC(I,K)-HES_cup(I,K)
        endif
     enddo
  enddo
  return
end subroutine cup_up_he


!----------------------------------------------------------------------
subroutine cup_up_moisture(ierr, z_cup, qc, qrc, pw, pwav,            &
     kbcon, ktop, mix, mgmxp, mkx, mgmzp, istart, iend, cd, dby,      &
     mentr_rate2, q, GAMMA_cup, zu, qes_cup, k22, qe_cup, &
      w_up, rho, ccn, TRIGG, iens,autoconv)
     
  implicit none
  INTEGER, INTENT(in) :: trigg, iens,autoconv
  integer istart, iend, mix, mgmxp, mkx, mgmzp, i, k
  real q(mgmxp,mgmzp), zu(mgmxp,mgmzp), GAMMA_cup(mgmxp,mgmzp),       &
       qe_cup(mgmxp,mgmzp), dby(mgmxp,mgmzp), cd(mgmxp,mgmzp),        &
       z_cup(mgmxp,mgmzp), qes_cup(mgmxp,mgmzp)
       
  real qc(mgmxp,mgmzp), qrc(mgmxp,mgmzp), pw(mgmxp,mgmzp),pwav(mgmxp)
  integer kbcon(mgmxp), ktop(mgmxp), ierr(mgmxp), k22(mgmxp), iall
 
  real radius, xl, dz, dh, qrch, c0, mentr_rate
  real mentr_rate2(mgmxp,mgmzp)
  
  real ccn (mgmxp)
  real qtc(mgmxp,mgmzp)
  real w_up(mgmxp,mgmzp)
  real con, q1
  real rho(mgmxp,mgmzp) ! g/cm^3
  
 ! real, parameter :: ANBASE =  50. !*1.e+6 !Berry-number at cloud base #/m^3(maritime)
 ! REAL, PARAMETER :: CCN =1000.  !*1.e+6 !Berry-number at cloud base #/m^3(continental)
  real, parameter :: BDISPM = 0.366       !Berry--size dispersion (maritime)
  REAL, PARAMETER :: BDISPC = 0.146       !Berry--size dispersion (continental)
  !integer, parameter :: berry = 2     ! if berry = 2, Berry parameterization on, else berry = 0 to 
                                      ! default autoconversion param.
				      
  integer, parameter:: deep=1,shallow=2
  
  IF(autoconv == 2 .AND. TRIGG /= 2) STOP "BERRY MUST BE USED WITH TRIGG = 2"

  iall=0
  xl=2.5e6
  c0=.002
  
  ! c0=.000   !para teste de conservacao (prec=0)
  ! cd= detrainment function 
  ! q = environmental q on model levels
  ! qe_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function 
  ! zu = normalized updraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! mentr_rate = entrainment rate
  ! qc = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! qrc = liquid water content in cloud after rainout
  ! pw = condensate that will fall out at that level
  ! pwav = totan normalized integrated condensate (I1)
  ! c0 = conversion rate (cloud to rain)

  !--- No precip for small clouds
  IF(iens.eq.shallow)  c0=0.

!do k=1,mkx
!     do i=istart,iend
!        if(mentr_rate2(i,k).gt.0.)then
!          radius = 0.2/(0.2/1200)
!	 if(radius.lt.900.)
!     !         if(radius.lt.900.)iall=0
!         endif
!     enddo
!enddo

 
  do i=istart,iend
     pwav(i)=0.
  enddo
  do k=1,mkx
     do i=istart,iend
        pw(i,k) =0.
        !_srf     qc(i,k) =qes_cup(i,k)
        qc(i,k) =qe_cup(i,k)
        qrc(i,k)=0.
     enddo
  enddo
  do i=istart,iend
     if(ierr(i).eq.0.)then
        do k=k22(i),kbcon(i)-1
           qc(i,k)=qe_cup(i,k22(i))
        enddo
     endif
  enddo

  do K=2,MKX-1
     do I=ISTART,IEND
      if (ierr(i).ne.0)  cycle
      if (K.lt.KBCON(I)) cycle
      if (K.gt.KTOP(I))  cycle

      DZ=Z_cup(i,K)-Z_cup(i,K-1)

      !--- 1. Steady state plume equation, for what could
      !---    be in cloud without condensation

      QC(i,K)=(QC(i,K-1)*(1.-.5*CD(i,K)*DZ)+mentr_rate2(i,k)*	     &
     	   DZ*Q(i,K-1))/(1.+mentr_rate2(i,k)*DZ-.5*cd(i,k)*dz)

      !--- 3.Condensation
      !   QTC  = total condensed water (var local)
      !--- Saturation  in cloud, this is what is allowed to be in it

      QRCH=QES_cup(I,K)+(1./XL)*(GAMMA_cup(i,k)/(1.+GAMMA_cup(i,k)))*DBY(I,K)


!------------------------- comeco do trecho modificado----------------------------
     IF(autoconv.eq.2 .and. iens.eq.deep) then

      
      ! total condensed water (QTC)
      QTC(i,k) = max(0., QC(I,K)-QRCH) ! kg[h2o]/kg[ar]= g[h2o]/g[ar]
        
      if(w_up(i,k) <= 0. .or. QTC(i,k) == 0. ) then
            CON = 0.
      
      else
      
            !- rho = air density (g[ar]/cm^3)
            q1 = rho(i,k)*QTC(i,k)  ! g[h2o]/cm^3
     
            !- Berry's formulation for autoconversion
	    CON = 1.e+6*q1*q1/(60.0*(5.0*rho(i,k) + 0.0366*CCN(i)/ &
	          ( 1.e6 * qtc(i,k) * BDISPC)  ) ) ! kg[h2o]/ ( kg[ar] s)
            
            !- w_up is the vertical velocity of the updraft (m/s)
 	    CON =      CON/(0.75*min(10.,w_up(i,k))) ! kg[h2o]/ ( kg[ar] m)
      endif
            
      !- rain water production (pw)
      pw(i,k)=CON*dz*zu(i,k)  !unit: kg[liq water]/kg[air]
      
      !- limit pw to the max allowed (the total condensed)
      pw(i,k)=min(pw(i,k),QTC(i,k))
      
      !- condensed water remained in cloud
      QRC(I,K) = QTC(i,k) - pw(i,k)

      !!!!pw(i,k)=pw(i,k)*zu(i,k)  !unit: kg[liq water]/kg[air]
     else
        !--- 3.Condensation
        !--- Liquid water content in cloud after rainout

        QRC(I,K)=(QC(I,K)-QRCH)/(1.+C0*DZ*zu(i,k))
        if (qrc(i,k).lt.0.) then
           qrc(i,k)=0.
        endif
  
        PW(i,k)=c0*dz*QRC(I,K)*zu(i,k)  !unit: kg[liq water]/kg[air]
                                        !unit of c0 is m^-1
        !if (iall.eq.1) then
        !   qrc(i,k)=0.
        !   pw(i,k)=(QC(I,K)-QRCH)*zu(i,k)
        !   if (pw(i,k).lt.0.) pw(i,k)=0.
        !endif
     endif
     !
     !print*,'cons=',i,k,j,pw(i,k)*1000.!,100.*(pw(i,k)/zu(i,k)+qrc(i,k)+qrch-QC(I,K))/(1.e-13+QC(I,K));call flush(6)
     
     
     !--- Set next level for the in cloud total water
     QC(I,K)=QRC(I,K)+qrch
     !
     !--- Integrated normalized ondensate
     !
     PWAV(I)=PWAV(I)+PW(I,K)

    enddo
  enddo
  return
end subroutine cup_up_moisture

!----------------------------------------------------------------------
subroutine cup_up_nms(zu, z_cup, entr2, cd, kbcon, ktop,               &
     mix, mgmxp, mkx, mgmzp, istart, iend, ierr, k22)
  implicit none
  integer i, k, mix, mgmxp, mkx, mgmzp, istart, iend
  real zu(mgmxp,mgmzp), z_cup(mgmxp,mgmzp), cd(mgmxp,mgmzp)
  integer kbcon(mgmxp), ktop(mgmxp), k22(mgmxp), ierr(mgmxp)
  real entr, dz
  real entr2(mgmxp,mgmzp)
  

  
  do k=1,mkx
     do i=istart,iend
        zu(i,k)=0.
     enddo
  enddo
  do i=istart,iend
     if (ierr(I).eq.0) then
        do k=k22(i),kbcon(i)
           zu(i,k)=1.
        enddo
        do K=KBcon(i)+1,KTOP(i)
           DZ=Z_cup(i,K)-Z_cup(i,K-1)
           ZU(i,K)=ZU(i,K-1)*(1.+(entr2(i,k)-cd(i,k))*DZ)
        enddo
     endif
  enddo
  return
end subroutine cup_up_nms


!----------------------------------------------------------------------
subroutine cup_up_aa0(aa0, z, zu, dby, GAMMA_CUP, t_cup,              &
     kbcon, ktop, mix, mgmxp, mkx, mgmzp, istart, iend, ierr)
  implicit none
  integer i, k, mix, mgmxp, mkx, mgmzp, istart, iend
  real aa0(mgmxp), z(mgmxp,mgmzp), zu(mgmxp,mgmzp),                   &
       gamma_cup(mgmxp,mgmzp), t_cup(mgmxp,mgmzp), dby(mgmxp,mgmzp)
  integer kbcon(mgmxp), ktop(mgmxp), ierr(mgmxp)
  real  dz, da
  
  
  do I=ISTART,IEND
     aa0(i)=0.
  enddo
  do K=2,MKX-1
  !DO 100 K=2,MKX-1
     do I=ISTART,IEND
     !DO 100 I=ISTART,IEND
        !IF (ierr(i).NE.0)  GO TO 100
        if (ierr(i).ne.0)  cycle
        !IF (K.LE.KBCON(I)) GO TO 100
        if (K.le.KBCON(I)) cycle
        !IF (K.GT.KTOP(I))  GO TO 100
        if (K.gt.KTOP(I))  cycle
        DZ = Z(I,K)-Z(I,K-1)
        da = zu(i,k)*DZ*(9.81/(1004.*((T_cup(I,K)))))*DBY(I,K-1)/  &
             (1.+GAMMA_CUP(I,K))
        !IF (K.EQ.KTOP(I).AND.da.LE.0.) go to 100
        if (K.eq.KTOP(I).and.da.le.0.) cycle
        AA0(I)=AA0(I)+da
        if (aa0(i).lt.0.) aa0(i)=0.

!100     CONTINUE
     enddo
  enddo

  return
end subroutine cup_up_aa0

