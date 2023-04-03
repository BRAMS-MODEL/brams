!--------------------------------------------------------------------
subroutine cup_env(j, ipr, jpr, z, qes, he, hes, t, q, p, z1,        &
     mix, mgmxp, mkx, mgmzp, istart, iend, psur, ierr, tcrit, itest)

  implicit none

  ! Arguments:
  integer, intent(in) :: j, ipr, jpr, mix, mgmxp, mkx, mgmzp, istart, &
       iend, itest
  integer, intent(in) :: ierr(mgmxp)
  real, intent(in)    :: tcrit
  real, intent(in)    :: t(mgmxp,mgmzp), p(mgmxp,mgmzp), &
       z1(mgmxp), psur(mgmxp)
  real, intent(inout) :: qes(mgmxp,mgmzp), q(mgmxp,mgmzp), z(mgmxp,mgmzp), &
       he(mgmxp,mgmzp), hes(mgmxp,mgmzp)

  ! Local Variables:
  integer :: i, k, iph,  m
  real    :: tv(mgmxp,mgmzp)
  real    :: AE(2), BE(2), HT(2), e, tvbar
  real    :: xl, cp

  XL=2.5E06
  cp=1004.
  HT(1)=XL/CP
  HT(2)=2.834E6/CP
  BE(1)=.622*HT(1)/.286
  AE(1)=BE(1)/273.+ALOG(610.71)
  BE(2)=.622*HT(2)/.286
  AE(2)=BE(2)/273.+ALOG(610.71)
  do K=1,MKX
     do I=ISTART,IEND
        if(ierr(i).eq.0)then
           !sgb - IPH is for phase, dependent on TCRIT (water or ice)
           IPH = 1
           if (T(I,K).le.TCRIT) IPH = 2
           E = exp(AE(IPH)-BE(IPH)/T(I,K))
           QES(I,K) = .622*E/(100.*P(I,K)-E)
           if (QES(I,K) .le. 1.E-08)   QES(I,K)=1.E-08
           if (Q(I,K)   .gt. QES(I,K)) Q(I,K)=QES(I,K)
           TV(I,K) = T(I,K)+.608*Q(I,K)*T(I,K)
        endif
     enddo
  enddo

  !--- z's are calculated with changed h's and q's and t's
  !--- if itest=2

  if (itest.ne.2) then
     do I=ISTART,IEND
        if (ierr(i).eq.0) then
           Z(I,1)=max(0.,Z1(I))-(ALOG(P(I,1))-ALOG(PSUR(I)))*287.*TV(I,1)/9.81

           ! **(JP)** passivel de eliminar print, permitindo vetorizacao
           ! **(JP)** se eliminar, muda resultado (pois vetoriza ALOG acima)
           !srf-----print-------
!!$           if (j.eq.jpr   .and. i.eq.ipr) then
!!$              k=1
!!$              if(k.eq.1) then
!!$                 write(6,'(a1,78a1)') ' ',('-',m=1,78)
!!$                 print*,'i k Z1(I) Z(I,k) P(I,k) PSUR(I) TV(I,k)'
!!$              endif
!!$              write(6,'(2i4,5F12.4)') i,k,Z1(I),Z(I,1),P(I,1),PSUR(I),TV(I,1)
!!$           endif
           !srf-----print-------
           ! **(JP)** fim de modificacao

        endif
     enddo

     !--- Calculate heights
     do K=2,MKX
        do I=ISTART,IEND
           if (ierr(i).eq.0) then
              TVBAR =  .5*TV(I,K)+.5*TV(I,K-1)
              Z(I,K) = Z(I,K-1)-(ALOG(P(I,K))-ALOG(P(I,K-1)))*287.*TVBAR/9.81

              ! **(JP)** passivel de eliminar print, permitindo vetorizacao
              ! **(JP)** se eliminar, muda resultado (pois vetoriza ALOG acima)
              !srf-----print-------
!!$              if (i.eq.ipr.and. j.eq.jpr) then
!!$                 if (k.eq.1) then
!!$                    write(6,'(a1,78a1)') ' ',('-',m=1,78)
!!$                    print*,'i k Z1(I) Z(I,k) P(I,k) PSUR(I) TV(I,k)'
!!$                 endif
!!$                 write(6,'(2i4,5F12.4)') i, k, Z1(I), Z(I,k), P(I,k),   &
!!$                      PSUR(I), TV(I,k)
!!$                 if (k.eq.mkx) write(6,'(a1,78a1)') ' ',('-',m=1,78)
!!$              endif
              !srf-----print-------
              ! **(JP)** fim de modificacao

           endif
        enddo
     enddo
  else
     do k=1,mkx
        do i=istart,iend
           if (ierr(i).eq.0) then
              z(i,k) = (he(i,k)-1004.*t(i,k)-2.5e6*q(i,k))/9.81
              z(i,k) = max(1.e-3,z(i,k))
           endif
        enddo
     enddo
  endif

  !--- calculate moist static energy - HE
  !--- Saturated moist static energy - HES

  do K=1,MKX
     do I=ISTART,IEND
        if (ierr(i).eq.0) then
           if (itest.eq.0) HE(I,K) = 9.81*Z(I,K)+1004.*T(I,K)+2.5E06*Q(I,K)
           HES(I,K) = 9.81*Z(I,K)+1004.*T(I,K)+2.5E06*QES(I,K)

           if (HE(I,K).ge.HES(I,K)) HE(I,K) = HES(I,K)

           ! **(JP)** elimita print, permitindo vetorizacao
           !srf-----print-------
!!$           if (i.eq.ipr.and.k.eq.jpr)then
!!$              if(k.eq.1) then
!!$                 write(6,'(a1,78a1)') ' ',('-',m=1,78)
!!$                 print*,'i,k,T(I,K),Q(I,K),Z(I,K),HE(I,K),HES(I,K)'
!!$              endif
!!$              write(6,'(2i4,5F12.4)') i,k,T(I,K),Q(I,K),Z(I,K),     &
!!$                   HE(I,K),HES(I,K)
!!$              if (k.eq.mkx) write(6,'(a1,78a1)') ' ',('-',m=1,78)
!!$           endif
           !srf-----print-------
           ! **(JP)** fim de modificacao

        endif
     enddo
  enddo

  return
end subroutine cup_env

!--------------------------------------------------------------------
subroutine cup_env_clev(j, ipr, jpr, t, qes, q, he, hes, z, p, qes_cup,  &
     q_cup, he_cup, hes_cup, z_cup, p_cup, gamma_cup, t_cup, psur,       &
     mix, mgmxp, mkx, mgmzp, istart, iend, ierr, z1)
  implicit none
  integer i, j, k, mix, mgmxp, mkx, mgmzp, istart, iend, ipr, jpr, m
  real qes_cup(mgmxp,mgmzp), q_cup(mgmxp,mgmzp), he_cup(mgmxp,mgmzp),  &
       hes_cup(mgmxp,mgmzp), z_cup(mgmxp,mgmzp), p_cup(mgmxp,mgmzp),   &
       gamma_cup(mgmxp,mgmzp), t_cup(mgmxp,mgmzp)
  real qes(mgmxp,mgmzp), q(mgmxp,mgmzp), he(mgmxp,mgmzp),              &
       hes(mgmxp,mgmzp), z(mgmxp,mgmzp), p(mgmxp,mgmzp),               &
       t(mgmxp,mgmzp), psur(mgmxp), z1(mgmxp)
  integer ierr(mgmxp)
  real xl, rv, cp
  xl=2.5e6
  rv=461.9
  cp=1004.
  do k=2,mkx
     do i=istart,iend
        if (ierr(i).eq.0)then
           qes_cup(i,k) = .5*(qes(i,k-1) + qes(i,k))
           q_cup(i,k)   = .5*(  q(i,k-1) +   q(i,k))
           hes_cup(i,k) = .5*(hes(i,k-1) + hes(i,k))
           he_cup(i,k)  = .5*( he(i,k-1) +  he(i,k))
           if (he_cup(i,k).gt.hes_cup(i,k)) he_cup(i,k) = hes_cup(i,k)

           z_cup(i,k) = .5*(z(i,k-1) + z(i,k))
           p_cup(i,k) = .5*(p(i,k-1) + p(i,k))
           t_cup(i,k) = .5*(t(i,k-1) + t(i,k))

           gamma_cup(i,k) =(xl/cp)*(xl/(rv*t_cup(i,k)*t_cup(i,k)))*qes_cup(i,k)
        endif
     enddo
  enddo
  do i=istart,iend
     if (ierr(i).eq.0)then
        qes_cup(i,1) = qes(i,1)
        q_cup(i,1)   = q(i,1)
        hes_cup(i,1) = hes(i,1)
        he_cup(i,1)  = he(i,1)

        !srf
        !        z_cup(i,1) = .5*( z(i,1) +   z1(i))
        !        p_cup(i,1) = .5*( p(i,1) + psur(i))
        !        t_cup(i,1) =      t(i,1)
        z_cup(i,1)   = z1(i)
        p_cup(i,1)   = psur(i)
        t_cup(i,1)   = t(i,1)
        !srf	
        gamma_cup(i,1) = xl/cp*(xl/(rv*t_cup(i,1)*t_cup(i,1)))*qes_cup(i,1)
     endif
  enddo

  ! **(JP)** elimita print, permitindo vetorizacao
  !srf-----print-------
!!$  do i=istart,iend
!!$     if (ierr(i).eq.0) then
!!$        do k=1,mkx
!!$           if (j.eq.jpr   .and. i.eq.ipr) then
!!$              if (k.eq.1) then
!!$                 write(6,'(a1,78a1)') ' ',('-',m=1,78)
!!$                 print*,'i k topo z t p z_cup t_cup p_cup'
!!$                 !	   print*,'i k Qes_cup Q_cup He_cup z_cup t_cup p_cup'
!!$              endif
!!$              !WRITE(6,'(2i4,6F12.4)') i,k,QES_CUP(I,K),Q_CUP(I,K),  &
!!$              !     He_cup(i,k),Z_CUP(I,K),T_CUP(I,K),P_CUP(i,k)
!!$              write(6,'(2i4,7F12.4)') i,k,z1(i),z(I,K)-z1(i),t(I,K),  &
!!$                   p(i,k),Z_CUP(I,K)-z1(i),T_CUP(I,K),P_CUP(i,k)
!!$
!!$              if (k.eq.mkx) write(6,'(a1,78a1)') ' ',('-',m=1,78)
!!$           endif
!!$       enddo
!!$     endif
!!$  enddo
  !srf-----print-------
  ! **(JP)** fim de modificacao

  return
end subroutine cup_env_clev

!--------------------------------------------------------------------
subroutine cup_direction2(i, j, dir, id, mix, mjx, mgmxp, mgmyp,  &
     massflx, iresult, num, imass, nall, maxens3, massfld)

  implicit none

  ! Arguments:
  integer, intent(in)  :: i, j, mix, mjx, mgmxp, mgmyp, num, imass, nall, maxens3
  integer, intent(out) :: iresult
  !srf      integer id(mix,mjx)
  !srf      real dir(mgmxp), massflx(mix,mjx)
  integer, intent(in)  ::  id(mgmxp,mgmyp)
  real, intent(inout)  ::  dir(mgmxp)
  real, intent(inout)  ::  massflx(mgmxp,mgmyp)
  real, intent(out)    ::  massfld

  ! Local Variables:
  integer :: k
  integer :: ia, ib, ja, jb
  real    ::  diff

  if (imass.eq.1) then
     massfld = massflx(i,j)
  endif
  iresult=0
  !      return

  diff = 22.5
  if (dir(i).lt.22.5) dir(i)=360.+dir(i)
  if (id(i,j).eq.1)   iresult=1
  !      ja=max(2,j-1)
  !      ia=max(2,i-1)
  !      jb=min(mjx-1,j+1)
  !      ib=min(mix-1,i+1)
  ja=j-1
  ia=i-1
  jb=j+1
  ib=i+1
  if (dir(i).gt.90.-diff.and.dir(i).le.90.+diff) then

     !--- Steering flow from the east
     if (id(ib,j).eq.1) then
        iresult=1
        if (imass.eq.1) then
           massfld = max(massflx(ib,j),massflx(i,j))
        endif
        return
     endif
  else if (dir(i).gt.135.-diff.and.dir(i).le.135.+diff)then

     !--- Steering flow from the south-east
     if (id(ib,ja).eq.1) then
        iresult=1
        if (imass.eq.1) then
           massfld = max(massflx(ib,ja),massflx(i,j))
        endif
        return
     endif

     !--- Steering flow from the south
  else if (dir(i).gt.180.-diff.and.dir(i).le.180.+diff) then
     if (id(i,ja).eq.1) then
        iresult=1
        if (imass.eq.1) then
           massfld = max(massflx(i,ja),massflx(i,j))
        endif
        return
     endif

     !--- Steering flow from the south west
  else if (dir(i).gt.225.-diff.and.dir(i).le.225.+diff) then
     if (id(ia,ja).eq.1) then
        iresult=1
        if (imass.eq.1) then
           massfld = max(massflx(ia,ja),massflx(i,j))
        endif
        return
     endif

     !--- Steering flow from the west
  else if (dir(i).gt.270.-diff.and.dir(i).le.270.+diff) then
     if (id(ia,j).eq.1) then
        iresult=1
        if (imass.eq.1)then
           massfld = max(massflx(ia,j),massflx(i,j))
        endif
        return
     endif

     !--- Steering flow from the north-west
  else if (dir(i).gt.305.-diff.and.dir(i).le.305.+diff) then
     if (id(ia,jb).eq.1) then
        iresult=1
        if (imass.eq.1) then
           massfld = max(massflx(ia,jb),massflx(i,j))
        endif
        return
     endif

     !--- Steering flow from the north
  else if (dir(i).gt.360.-diff.and.dir(i).le.360.+diff) then
     if (id(i,jb).eq.1) then
        iresult=1
        if (imass.eq.1) then
           massfld = max(massflx(i,jb),massflx(i,j))
        endif
        return
     endif

     !--- Steering flow from the north-east
  else if (dir(i).gt.45.-diff.and.dir(i).le.45.+diff) then
     if (id(ib,jb).eq.1) then
        iresult=1
        if (imass.eq.1) then
           massfld = max(massflx(ib,jb),massflx(i,j))
        endif
        return
     endif
  endif

  !srf
  !      if(massfld.gt.0.) print*,'---- MASS FLD=',massfld        
  !srf

  return
end subroutine cup_direction2

!--------------------------------------------------------------------
subroutine cup_ktop(ilo, dby, kbcon, ktop, mix, mgmxp, mkx, mgmzp,   &
     istart, iend, ierr)
  implicit none
  integer mix, mgmxp, mkx, mgmzp, i, k, istart, iend, ilo
  real dby(mgmxp,mgmzp)
  integer ierr(mgmxp), kbcon(mgmxp), ktop(mgmxp)

  do I=ISTART,IEND
  !DO 42 I=ISTART,IEND
     ktop(i)=1
     if (ierr(I).eq.0) then

        do K=KBCON(I)+1,MKX-2
        !DO 40 K=KBCON(I)+1,MKX-2
           if (DBY(I,K).le.0.) then
              KTOP(I) = K-1
              GO TO 41
           endif
!40         CONTINUE
        enddo
        if (ilo.eq.1) ierr(i)=5
        if (ilo.eq.2) ierr(i)=998
        GO TO 42
41      continue
        do k=ktop(i)+1,mkx
           dby(i,k)=0.
        enddo
     endif
42   continue
  enddo
  return
end subroutine cup_ktop

!--------------------------------------------------------------------
subroutine cup_kbcon(iloop, k22, kbcon, he_cup, hes_cup, mix, mgmxp,  &
     mkx, mgmzp, istart, iend, ierr, kbmax, p_cup, cap_max)
  implicit none
  integer i, mix, mgmxp, mkx, mgmzp, istart, iend, iloop
  integer kbcon(mgmxp), k22(mgmxp), ierr(mgmxp), kbmax(mgmxp)
  real he_cup(mgmxp,mgmzp), hes_cup(mgmxp,mgmzp), p_cup(mgmxp,mgmzp)
  real pbcdif, cap_max(mgmxp)
  !srf
  !      integer k
  !srf

  !--- Determine the level of convective cloud base  - KBCON

  do I=ISTART,IEND
  !DO 27 I=ISTART,IEND
     kbcon(i)=1
     !IF (ierr(I).NE.0) GO TO 27
     if (ierr(I).ne.0) cycle
     KBCON(I)=K22(I)

     !srf
     !       print*,'-----------------------------------------'
     !       print*,'1',i,ierr(I),k22(i),kbcon(i),cap_max(i)
     !       do k=mkx,k22(i),-1
     !       do k=mkx,1,-1
     !       print*,k,HE_cup(I,K),HES_cup(I,K),p_cup(i,k)
     !       enddo
     !srf

     GO TO 32
31   continue
     KBCON(I)=KBCON(I)+1
     if (KBCON(I).gt.KBMAX(i)+2) then
        if(iloop.eq.1)ierr(i)=3
        if(iloop.eq.2)ierr(i)=997
        !srf
        !       print*,'2',i,ierr(I),k22(i),kbcon(i),cap_max(i)
        !srf

        !GO TO 27
        cycle
     endif
32   continue
     if (HE_cup(I,K22(I)).lt.HES_cup(I,KBCON(I))) GO TO 31

     !     Cloud base pressure and max moist static energy pressure
     !     i.e., the depth (in mb) of the layer of negative buoyancy

     if (KBCON(I)-K22(I).eq.1) then
        !srf
        !       print*,'2',i,ierr(I),k22(i),kbcon(i),cap_max(i)
        !srf
        !GO TO 27
        cycle
     endif

     PBCDIF = -P_cup(I,KBCON(I)) + P_cup(I,K22(I))
     if (PBCDIF.gt.cap_max(i)) then
        K22(I)   = K22(I)+1
        KBCON(I) = K22(I)
        GO TO 32
     endif
!27   CONTINUE
  enddo

  return
end subroutine cup_kbcon


!--------------------------------------------------------------------
subroutine cup_kbcon_cin(iloop, k22, kbcon, he_cup, hes_cup, z, tmean, qes,  &
            mix, mgmxp, mkx, mgmzp, istart, iend, ierr, kbmax, p_cup, cap_max)
  implicit none
  !srf      include 'METCON'
  integer i, mix, mgmxp, mkx, mgmzp, istart, iend, iloop
  integer kbcon(mgmxp), k22(mgmxp), ierr(mgmxp), kbmax(mgmxp)
  real he_cup(mgmxp,mgmzp), hes_cup(mgmxp,mgmzp), p_cup(mgmxp,mgmzp)
  real z(mgmxp,mgmzp), tmean(mgmxp,mgmzp), qes(mgmxp,mgmzp)
  real pbcdif, cap_max(mgmxp)
  real cin, cin_max, dh, tprim, gamma

  !srf - Substitue o METCON include file
  real LOVCP_P
  real xl, rv, cp
  xl=2.5e6
  rv=461.525
  cp=1004.
  LOVCP_P = xl/cp
  !srf- -----------------------------------
  !
  !--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
  !
  do I=ISTART,IEND
  !DO 27 I=ISTART,IEND
     cin_max=-cap_max(i)

     kbcon(i)=1
     cin = 0.
     !IF (ierr(I).NE.0) GO TO 27
     if (ierr(I).ne.0) cycle
     KBCON(I)=K22(I)
     GO TO 32
31   continue
     KBCON(I)=KBCON(I)+1
     if (KBCON(I).gt.KBMAX(i)+2) then
        if (iloop.eq.1) ierr(i)=3
        if (iloop.eq.2) ierr(i)=997

        !srf
        !      print*,'KCON_CIN ',i,k22(i),kbcon(i),ierr(i),cin
        !srf

        !GO TO 27
        cycle
     endif
32   continue
     dh = HE_cup(I,K22(I)) - HES_cup(I,KBCON(I))
     if (dh.lt. 0.) then
        GAMMA = LOVCP_P*(xl/(rv*(Tmean(I,K22(i))**2)))*QES(I,K22(i))
        tprim = dh/(cp*(1.+gamma))

        cin   = cin + 9.8066*tprim* &
	        (z(i,k22(i))-z(i,k22(i)-1))/tmean(i,k22(i))

        go to 31
     end if

     !     If negative energy in negatively buoyant layer
     !       exceeds convective inhibition (CIN) threshold,
     !       then set K22 level one level up and see if that
     !       will work.

     if (cin.lt.cin_max) then
        K22(I)=K22(I)+1
        KBCON(I)=K22(I)
        GO TO 32
     endif

     !srf
     !      print*,'KCON_CIN ',i,k22(i),kbcon(i),ierr(i),cin
     !srf

!27   CONTINUE
  enddo

  return
end subroutine cup_kbcon_cin


!--------------------------------------------------------------------
subroutine MINIMI(ARRAY, MXX, mgmxp, MZX, mgmzp, KS, KEND,  &
     KT, ISTART, IEND, ierr)

  implicit none

  integer MXX, MZX, ISTART, IEND, I, K, mgmxp, mgmzp
  real ARRAY(mgmxp,mgmzp), X(mgmxp)
  integer KT(mgmxp), KS(mgmxp), KEND(mgmxp), KSTOP, ierr(mgmxp)

  do I=ISTART,IEND
  !DO 200 I=ISTART,IEND
     KT(I)=KS(I)
     if (ierr(i).eq.0) then
        X(I)=ARRAY(I,KS(I))
        KSTOP=max(KS(I)+1,KEND(I))

        do K=KS(I)+1,KSTOP
        !DO 100 K=KS(I)+1,KSTOP
           if(ARRAY(I,K).lt.X(I)) then
              X(I)=ARRAY(I,K)
              KT(I)=K
           endif
!100        CONTINUE
        enddo
     endif
!200  CONTINUE
  enddo

  return
end subroutine MINIMI

!--------------------------------------------------------------------
subroutine MAXIMI(ARRAY, MXX, mgmxp, MZX, mgmzp, KS, KE, MAXX, ISTART, &
     IEND, ierr)

  implicit none

  integer mgmxp, mgmzp, MXX, MZX, KS, KE(mgmxp), ISTART, IEND, I, K
  real ARRAY(mgmxp,mgmzp), X(mgmxp)
  real XAR
  integer MAXX(mgmxp), ierr(mgmxp)

  do I=ISTART,IEND
  !DO 200 I=ISTART,IEND
     MAXX(I)=KS
     if(ierr(i).eq.0)then
        X(I)=ARRAY(I,KS)

        do K=KS,KE(i)
        !DO 100 K=KS,KE(i)
           XAR=ARRAY(I,K)
           if(XAR.ge.X(I)) then
              X(I)=XAR
              MAXX(I)=K
           endif
!100        CONTINUE
        enddo
     endif
!200  CONTINUE
  enddo

  return
end subroutine MAXIMI
