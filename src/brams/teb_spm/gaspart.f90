!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

SUBROUTINE le_fontes(ng, n1, n2, n3, np, ia, iz, ja, jz, time)

  USE mem_grid, ONLY: &    !INTENT(IN)
       ngrids,        &    !INTENT(IN)
       grid_g,        &    !INTENT(IN)
!!$       time,          &    !INTENT(IN)
!!$       zt,            &    !INTENT(IN)
       dzt                 !INTENT(IN)
  USE mem_leaf, ONLY: &
       leaf_g              !INTENT(IN)
  USE mem_basic, ONLY: &
       basic_g             !INTENT(IN)
  USE mem_gaspart, ONLY: &
       gaspart_g           !INTENT(INOUT)

!!$  USE node_mod, ONLY: mchnum, master_num !INTENT(IN) - DEBUG
  USE node_mod, ONLY: mynum !INTENT(IN) - DEBUG

  IMPLICIT NONE
  ! Arguments:
  INTEGER, INTENT(in) :: ng, n1, n2, n3, np, ia, iz, ja, jz
  REAL, INTENT(in)    :: time
  ! Local Variables:
  INTEGER :: ig, kl, i, j
  INTEGER, PARAMETER :: ngases=6
  INTEGER :: kgas(ngases)
  INTEGER :: len1
  CHARACTER(len=10) :: gas(ngases), tracer

  print *, "DEBUG-ALF:le_fontes:mynum,ng,n1,n2,n3,np,ia,iz,ja,jz=", &
       mynum,ng,n1,n2,n3,np,ia,iz,ja,jz
  call flush(6)

  ! Init arrays: gas, kgas
  gas(:)  = (/'NO  ', 'NO2 ', 'PM25', 'CO  ', 'SO2 ', 'VOCS'/)
  kgas(:) = (/  1,    4,     7,     10,   13,     16/)
!!$  data ( gas(ig),ig=1,ngases) /'NO', 'NO2', 'PM25', 'CO', 'SO2', 'VOCS'/
!!$  data (kgas(ig),ig=1,ngases) /  1,    4,     7,     10,   13,     16/

  !------------ Allocation table of qsc = dum1------------------------
  ! Distribute sources in order to avoid highly localized sources:
  ! (somente definidas em 1 ponto de grade)
  ! kgas(1) = 1,2,3    => NO
  ! kgas(2) = 4,5,6    => NO2
  ! kgas(3) = 7,8,9    => PM25
  ! kgas(4) = 10,11,12 => CO
  ! kgas(5) = 13,14,15 => SO2 
  ! kgas(6) = 16,17,18 => VOCS
  !----------------------------------------------------------------------


  !sources only in the inner grid
  IF (ng/=ngrids) RETURN

  !Allocate memory in DUM1 vector for sources:
  ! first kl 3D vector gaspart_g(ng)%gasr(1,1,1) levels are left for
  ! emission sources calculation
  kl = 3*ngases 

  gaspart_g(ng)%gasr(:,:,:) = 0.

  DO ig=1,ngases

     !Call dum1_zero(n1,n2,n3,ia,iz,ja,jz,kl,scalar_g(ig,ng)%dum1(1,1,1))  

     tracer = gas(ig)

     len1   = LEN_TRIM(tracer) + 1

     CALL read_sources_teb(ng, n1, n2, n3, np, ia, iz, ja, jz, &
          gaspart_g(ng)%gasr(1,1,1), gaspart_g(ng)%fusog(1,1), &
          tracer(1:len1), kgas(ig), leaf_g(ng)%G_URBAN(1,1,1), &
          grid_g(ng)%dxt(1,1), grid_g(ng)%dyt(1,1), time       )

     CALL reorganize_sources_teb(ng, n1, n2, n3, ia, iz, ja, jz, &
          gaspart_g(ng)%gasr(1,1,1), kgas(ig)                )

     CALL convert_to_misture_ratio_teb(ng, n1, n2, n3, ia, iz, ja, jz,     &
          kgas(ig), gaspart_g(ng)%gasr(1,1,1),                             &
          basic_g(ng)%dn0(1,1,1), grid_g(ng)%rtgt(1,1),                    &
          grid_g(ng)%dxt(1,1), grid_g(ng)%dyt(1,1), dzt                    )

  ENDDO !end of gases' looping

!!$  if (mchnum==master_num) then
!!$     print *, "DEBUG-ALF:le_fontes:mchnum,time=", mchnum,time
!!$     call flush(6)
!!$  endif

END SUBROUTINE le_fontes

!--------------------------------------------------

SUBROUTINE read_sources_teb(ng, n1, n2, n3, np, ia, iz, ja, jz, &
     qsc, fuso, gas, kgas, schar, dxt, dyt, time)

  USE mem_emiss, ONLY: &
       EINDNO,         & !INTENT(IN)
       EINDNO2,        & !INTENT(IN)
       EINDPM,         & !INTENT(IN)
       EINDCO,         & !INTENT(IN)
       EINDSO2,        & !INTENT(IN)
       EINDVOC,        & !INTENT(IN)
       EVEINO,         & !INTENT(IN)
       EVEINO2,        & !INTENT(IN)
       EVEIPM,         & !INTENT(IN)
       EVEICO,         & !INTENT(IN)
       EVEISO2,        & !INTENT(IN)
       EVEIVOC,        & !INTENT(IN)
       EFSAT,          & !INTENT(IN)
       EFSUN,          & !INTENT(IN)
       WEEKDAYIN         !INTENT(IN)
       
  USE teb_vars_const, ONLY : &
       RUSHH1,               & !INTENT(IN)
       RUSHH2,               & !INTENT(IN)
       DAYLIGHT                !INTENT(IN)

  IMPLICIT NONE
  ! Arguments:
  INTEGER, INTENT(IN)          :: ng, n1, n2, n3, np, ia, iz, ja, jz
  REAL, INTENT(INOUT)          :: qsc(n1,n2,n3)
  REAL, INTENT(IN)             :: fuso(n2,n3)
  CHARACTER(len=*), INTENT(IN) :: gas
  INTEGER, INTENT(IN)          :: kgas
  REAL, INTENT(IN)             :: schar(n2,n3,np)
  REAL, INTENT(IN)             :: dxt(n2,n3)
  REAL, INTENT(IN)             :: dyt(n2,n3)
  REAL, INTENT(IN)             :: time
  ! Local Variables:
  INTEGER :: i, j, idays
  CHARACTER(len=3) :: cday
  REAL :: pfat, pfat2, emiss, emiind, area, pft, r_q, tign, ax1, ax2, bx1, bx2
  REAL :: timeq1, timeq2

!!$  real :: sumqsc !DEBUG-ALF

  !****************************************************************************
  !defining the emission rate for each gas/particle
  !****************************************************************************

  !The following values are provided by CETESB's report for 2001
  !these emission rates are relative to an area of 8 million m2, wich 
  !represent the total area of RMSP. However, the urbanized part of this
  !area is about 1.5 and it will be considered here.

  !    veicular emissions rate in ug/min/m2 (for an area of 8 million of m2)

  !             CO   = 390.6
  !             NOx  =  89.3
  !             NO   = NOx*0.9
  !             NO2  = NOx*0.1
  !             VOC's=  91.0
  !             SO2  =   5.3

  ! In order to consider the diurnal cycle of veicular emission, the emission
  ! rates will have units of kg/m2/day

  !    veicular emissions rate in kg/day/m2 (for an area of 1.5 thousand of m2)

  !             CO   = 3.0173515E-03
  !             NOx  = 6.8566209E-04
  !             NO   = NOx*0.9 =  6.1709585E-04
  !             NO2  = NOx*0.1 =  6.8566209E-05
  !             VOC's= 6.9954334E-04
  !             SO2  = 4.0730592E-05
  !             PM15 = 6.2648396E-06
  !    industrial emissions rate in kg/s/m2 (for an area of 1.5 thousand of m2)
  !
  !  emiind
  !             CO   = 8.1599860E-10
  !             NOx  = 6.8566209E-04
  !             NO   = 2.6636227E-10
  !             NO2  = 2.9595805E-11
  !             VOC's= 2.5367833E-10
  !             SO2  = 3.6149164E-10
  !             PM25 = 4.3421278E-10

  IF     (gas=='CO')   THEN
     emiss  = EVEICO 
     emiind = EINDCO
!!$     print *, "DEBUG-ALF:read_sources_teb:gas=", gas
!!$     call flush(6)
  ELSEIF (gas=='NO')   THEN
     emiss  = EVEINO
     emiind = EINDNO
!!$     print *, "DEBUG-ALF:read_sources_teb:gas=", gas
!!$     call flush(6)
  ELSEIF (gas=='NO2')  THEN
     emiss  = EVEINO2
     emiind = EINDNO2
!!$     print *, "DEBUG-ALF:read_sources_teb:gas=", gas
!!$     call flush(6)
  ELSEIF (gas=='VOCS') THEN
     emiss  = EVEIVOC
     emiind = EINDVOC
!!$     print *, "DEBUG-ALF:read_sources_teb:gas=", gas
!!$     call flush(6)
  ELSEIF (gas=='SO2')  THEN
     emiss  = EVEISO2
     emiind = EINDSO2
!!$     print *, "DEBUG-ALF:read_sources_teb:gas=", gas
!!$     call flush(6)
  ELSEIF (gas=='PM25') THEN
     emiss  = EVEIPM
     emiind = EINDPM
!!$     print *, "DEBUG-ALF:read_sources_teb:gas=", gas
!!$     call flush(6)
  ENDIF

!!$  emiss=emiss

  pft   = emiss/273234.9
  !the value of 273234.9 was obtained in order to have the integral for
  !one day = emiss
  !it is used to distribute emissions following the diurnal cycle 

  tign  = 0.0*3600.                            !UTC time of ignition 

  !print*,'valor de time=',time

  idays = INT((time/3600.)/24.0 + 0.00001)  !number of days of simulation
  tign  = tign + REAL(idays)*24.0*3600.0

  pfat  = 1.0

  CALL EMFACTOR(WEEKDAYIN,idays,cday)
  IF     (cday=='SAT') THEN
     pfat=EFSAT
  ELSEIF (cday=='SUN') THEN
     pfat=EFSUN
  ENDIF

  !******************************************************************
  !gaussian distribution                                            *
  !                                     timeqi                      *
  !                                      ____                       *
  !                               _     |    |      _               *
  !              1               |      (x-mi)**2    |              *
  !f(x)= ----------------- * EXP | -  -------------- |  i=1,2       *
  !      sigmai* sqrt(2*PI)      |    2*(sigmai**2)  |              *
  !     |__________________|      -  |____________| -               *
  !             axi                    8.5 or 10.6                  *
  !******************************************************************

  ax1 = 4.5
  ax2 = 9.2

  DO i = 1,n2
     DO j= 1,n3

        IF (NINT(SCHAR(I,J,2))/=0) THEN

           !bx1=10.81
           !bx2=19.0
           bx1 = RUSHH1 - fuso(i,j) + DAYLIGHT
           bx2 = RUSHH2 - fuso(i,j) + DAYLIGHT

           timeq1 = (time/3600.0 - tign/3600.0) - bx1
           timeq2 = (time/3600.0 - tign/3600.0) - bx2

           r_q = pft*(ax1*EXP(-(timeq1)**2/8.5) + &
                ax2*EXP(-(timeq2)**2/10.6))*pfat 

           area = 1./(dxt(i,j)*dyt(i,j))

           !    Identificando pontos de grade com a classificacao urbano 1 
           IF     (NINT(SCHAR(I,J,2))==1) THEN
              pfat2 = 1.
              ! Identificando pontos de grade com a classificacao urbano 2 
           ELSEIF (NINT(SCHAR(I,J,2))==2) THEN
              pfat2 = 0.33333333
              ! Identificando pontos de grade com a classificacao urbano 3 
           ELSEIF (NINT(SCHAR(I,J,2))==3) THEN
              pfat2 = 0.2
           ENDIF

           qsc(kgas,i,j) = (r_q*area*pfat*pfat2) + (emiind*area)

        ENDIF

     ENDDO
  ENDDO

!!$  sumqsc = sum(qsc(kgas,:,:))
!!$  print *, "DEBUG-ALF:read_sources_teb:kgas,sumqsc=", kgas,sumqsc

END SUBROUTINE read_sources_teb

!-------------------------------------------------------

SUBROUTINE reorganize_sources_teb(ngrid, n1, n2, n3, ia, iz, ja, jz, qsc, kgas)

  USE node_mod, ONLY: &
       nodeibcon,     & !INTENT(IN)
       mynum            !INTENT(IN)

  IMPLICIT NONE
  ! Arguments:
  INTEGER, INTENT(IN) :: ngrid, n1, n2, n3, ia, iz, ja, jz
  REAL, INTENT(INOUT) :: qsc(n1,n2,n3)
  INTEGER, INTENT(IN) :: kgas
  ! Local Variables:
  INTEGER :: i,j, ii, jj
  REAL    :: f
  INTEGER         :: jInit, jEnd, iInit, iEnd

!!$  real :: sumqsc !DEBUG-ALF
  !fator de distribuicao 20% para cada um do 9 primeiros vizinhos
  !( incluindo o proprio site i,j)

  f = 0.2

  ! Defining bounds
  if (btest(nodeibcon(mynum,ngrid),0)) then !OESTE
     iInit = 3
  else
     iInit = 1
  endif
  if (btest(nodeibcon(mynum,ngrid),1)) then !LESTE
     iEnd = n2-2
  else
     iEnd = n2
  endif
  if (btest(nodeibcon(mynum,ngrid),2)) then !NORTE
     jInit = 3
  else
     jInit = 1
  endif
  if (btest(nodeibcon(mynum,ngrid),3)) then !SUL
     jEnd = n3-2
  else
     jEnd = n3
  endif

  !do j=3,n3-2
     !do i=3,n2-2
  !DO j=1,n3
     !DO i=1,n2
  DO j=jInit,jEnd
     DO i=iInit,iEnd

        ! ponto grade i,j,k(=2)
        qsc(kgas+1,i,j) = qsc(kgas+1,i,j) + (1.-f)*qsc(kgas,i,j)

        !distribuicao nos 9 sites em torno de i,j    !  j+1  .   .   .
        !do jj = j-1,j+1                          !   j   .   .   .  K=2     
           !do ii = i-1,i+1                        !  j-1  .   .   .
              !      i-1  i  i+1
        DO jj = MAX(j-1,1),MIN(j+1,n3)
           DO ii = MAX(i-1,1),MIN(i+1,n2)

              qsc(kgas+1,ii,jj) = qsc(kgas+1,ii,jj) + &
                   (1./9.)*f*qsc(kgas,i,j)  ! ponto grade i,j,k= 2 e 3

           ENDDO
        ENDDO

     ENDDO
  ENDDO

!!$  sumqsc = sum(qsc(kgas,:,:))
!!$  print *, "DEBUG-ALF:reorganize_sources_teb:kgas,sumqsc=", kgas,sumqsc
!!$  sumqsc = sum(qsc(kgas+1,:,:))
!!$  print *, "DEBUG-ALF:reorganize_sources_teb:kgas+1,sumqsc=", kgas+1,sumqsc
!!$  call flush(6)

END SUBROUTINE reorganize_sources_teb

!-------------------------------------------------------------

SUBROUTINE convert_to_misture_ratio_teb(ng, n1, n2, n3, ia, iz, ja, jz, &
     kgas, qsc, dn0, rtgt, dxt, dyt, dzt) !zt

  IMPLICIT NONE
  ! Arguments:
  INTEGER, INTENT(IN) :: ng, n1, n2, n3, ia,iz,ja,jz, kgas
  REAL, INTENT(INOUT) :: qsc(n1,n2,n3)
  REAL, INTENT(IN)    :: dn0(n1,n2,n3)
  REAL, INTENT(IN)    :: rtgt(n2,n3)
  REAL, INTENT(IN)    :: dxt(n2,n3)
  REAL, INTENT(IN)    :: dyt(n2,n3)
  REAL, INTENT(IN)    :: dzt(n1)
  ! Local variables:
!!$  real, dimension(n1) ::  zt
  REAL    :: fcu, vol	
  INTEGER :: i, j, k
  !

  !Fator de conversao de unidades
  fcu = 1.        !=> kg [gas/part] /kg [ar]
  !!fcu =1.e+12   !=> ng [gas/part] /kg [ar]
  !fcu =1.e+6      !=> mg [gas/part] /kg [ar]  

  DO j = 1,n3
     DO i = 1,n2
        DO k = 2,2

           vol             = 1./(dxt(i,j)*dyt(i,j)*dzt(k)*rtgt(i,j))

           qsc(kgas+1,i,j) = fcu*qsc(kgas+1,i,j)/(vol*dn0(k,i,j))

        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE convert_to_misture_ratio_teb

!------------------------------------------------------------

SUBROUTINE sources_teb(n1, n2, n3, ia, iz, ja, jz, ig, ngrids)

  USE mem_gaspart, ONLY : gaspart_g

  IMPLICIT NONE
  ! Arguments:
  INTEGER, INTENT(in) :: n1, n2, n3, ia, iz, ja, jz, ig, ngrids

  IF (ig==ngrids) THEN
     CALL EMISSAO(n1, n2, n3, ia, iz, ja, jz, gaspart_g(ig))
  ENDIF

END SUBROUTINE sources_teb

!------------------------------------------------------------

SUBROUTINE EMISSAO(n1, n2, n3, ia, iz, ja, jz, gaspart)

  USE mem_gaspart, ONLY : gaspart_vars ! Type

  IMPLICIT NONE
  ! Arguments:
  INTEGER, INTENT(in)               :: n1, n2, n3, ia, iz, ja, jz
  TYPE(gaspart_vars), INTENT(inout) :: gaspart

  CALL tendgas(n1,n2,n3,ia,iz,ja,jz,gaspart%pnot  (1),gaspart%gasr(1,1,1),2)
  CALL tendgas(n1,n2,n3,ia,iz,ja,jz,gaspart%pno2t (1),gaspart%gasr(1,1,1),5)
  CALL tendgas(n1,n2,n3,ia,iz,ja,jz,gaspart%ppm25t(1),gaspart%gasr(1,1,1),8)
  CALL tendgas(n1,n2,n3,ia,iz,ja,jz,gaspart%pcot  (1),gaspart%gasr(1,1,1),11)
  CALL tendgas(n1,n2,n3,ia,iz,ja,jz,gaspart%pso2t (1),gaspart%gasr(1,1,1),14)
  CALL tendgas(n1,n2,n3,ia,iz,ja,jz,gaspart%pvoct (1),gaspart%gasr(1,1,1),17)

END SUBROUTINE EMISSAO

!---------------------------------------------------------------

SUBROUTINE TENDGAS(n1, n2, n3, ia, iz, ja, jz, tende, qsc, k)

  IMPLICIT NONE
  ! Arguments:
  INTEGER, INTENT(in) :: n1, n2, n3, ia, iz, ja, jz, k
  REAL, INTENT(inout) :: tende(n1,n2,n3)
  REAL, INTENT(in)    :: qsc(n1,n2,n3)
  ! Local Variables:
  INTEGER :: i, j

!!$  real :: sumtendea, sumtendeb !DEBUG-ALF
!!$
!!$  sumtendea = sum(tende(2,ia:iz,ja:jz))

  DO i=ia,iz
     DO j=ja,jz
        tende(2,i,j) = tende(2,i,j) + qsc(k,i,j)
     ENDDO
  ENDDO

!!$  sumtendeb = sum(tende(2,ia:iz,ja:jz))
!!$  print *, "DEBUG-ALF:TENDGAS:k,sumtendea,b=", k,sumtendea,sumtendeb

END SUBROUTINE TENDGAS

!!$!---------------------------------------------------------------
!!$
!!$subroutine dum1_zero(n1,n2,n3,ia,iz,ja,jz,kl,d)
!!$
!!$  implicit none
!!$  integer :: n1,n2,n3,ia,iz,ja,jz,kl,i,j,k
!!$  real, dimension(n1,n2,n3) :: d
!!$
!!$  do j=ja,jz
!!$     do i=ia,iz
!!$        do k=1,kl 
!!$           d(k,i,j) = 0.
!!$        enddo
!!$     enddo
!!$  enddo
!!$  return
!!$end subroutine dum1_zero

!----------------------------------------------------------------

SUBROUTINE init_conc1(ictrl, ngrid, n1, n2, n3, np, &
     G_URBAN, no,                                   &
     no2, pm25,                                     &
     co, voc,                                       &
     so2, so4, aer, zt                              )

  USE node_mod, ONLY: &
       nodeibcon,     & !INTENT(IN)
       mynum            !INTENT(IN)

  IMPLICIT NONE
  ! Arguments:
  INTEGER, INTENT(in) :: ictrl, ngrid, n1, n2, n3, np
  REAL, INTENT(inout) :: G_URBAN(n2,n3,np)
  REAL, INTENT(inout) :: no(n1,n2,n3)
  REAL, INTENT(inout) :: no2(n1,n2,n3)
  REAL, INTENT(inout) :: pm25(n1,n2,n3)
  REAL, INTENT(inout) :: co(n1,n2,n3)
  REAL, INTENT(inout) :: voc(n1,n2,n3)
  REAL, INTENT(inout) :: so2(n1,n2,n3)
  REAL, INTENT(inout) :: so4(n1,n2,n3)
  REAL, INTENT(inout) :: aer(n1,n2,n3)
  REAL, INTENT(inout) :: zt(n1)
  ! Local Variables:
  REAL, PARAMETER :: PMAR=28.97
  INTEGER         :: i, j, ii, jj, k
  REAL            :: pfat, f, expo
  REAL            :: nox(n1,n2,n3)
  REAL            :: no2x(n1,n2,n3)
  REAL            :: pm25x(n1,n2,n3)
  REAL            :: cox(n1,n2,n3)
  REAL            :: vocx(n1,n2,n3)
  REAL            :: so2x(n1,n2,n3)
  REAL            :: no0, no20, pm250, co0, voc0, so20
  INTEGER         :: jInit, jEnd, iInit, iEnd

  no(:,:,:)   = 0.			       !NO
  no2(:,:,:)  = 0.			       !NO2
  pm25(:,:,:) = 0.			       !PM25
  co(:,:,:)   = 0.			       !CO
  so2(:,:,:)  = 0.			       !SO2
  so4(:,:,:)  = 0.			       !SO4
  voc(:,:,:)  = 0.			       !VOCS
  aer(:,:,:)  = 0.			       !VOCS

!!$  print *, "DEBUG-ALF:init_conc1:ictrl=", ictrl
!!$  print *, "DEBUG-ALF:init_conc1:sum(G_URBAN(:,:,2))=", sum(G_URBAN(:,:,2))
!!$  call flush(6)

  IF (ictrl==1) THEN
     !--------
     !Flag which defines initial condition (IF RUNTYPE = INITIAL)
     DO j=1,n3
        DO i=1,n2

           pfat = 0.0

           !identify surface type by using G_URBAN 
           IF (NINT(G_URBAN(i,j,2))/=0.) THEN

              IF (NINT(G_URBAN(i,j,2))==1) THEN
                 pfat = 1.
              ENDIF

              IF (NINT(G_URBAN(i,j,2))==2) THEN
                 pfat = 0.3     !using only 30% of values for urban regions
              ENDIF
              IF (NINT(G_URBAN(i,j,2))==3) THEN
                 pfat = 0.2     !using only 20% of values for urban regions
              ENDIF
           ELSE
              pfat = 0.1
           ENDIF

           !defining surface backgroud concentrations (typical values for 00 Z in urban regions)
           !can be better difined if using real values (eg. read a concentration file)

           no0   = 248.04 * pfat  !NO  (estacao congonhas dia 31/07/99 ug/m3)
           no20  =  62.30 * pfat  !NO2 (estacao congonhas dia 31/07/99 ug/m3)
           pm250 =  17.28 * pfat  !PM25(estacao santana   dia 31/07/99 ug/m3)
           co0   =   1.20 * pfat  !CO  (estacao congonhas dia 31/07/99 ppm)
           so20  =  17.50 * pfat  !SO2 (estacao congonhas dia 31/07/99 ug/m3)
           voc0  =   0.31 * pfat  !VOCS (Leila)

           DO k=2,n1-1
              expo         = EXP(-zt(k)/zt(n1-1))

              nox  (k,i,j) = (no0*expo)  *(1.e-9)/1.275	     !NO
              no2x (k,i,j) = (no20*expo) *(1.e-9)/1.275	     !NO2
              pm25x(k,i,j) = (pm250*expo)*(1.e-9)/1.275	     !PM25
              cox  (k,i,j) = (co0*expo)  *(28.00/PMAR)*1.e-6 !CO
              so2x (k,i,j) = (so20*expo) *(1.e-9)/1.275	     !SO2
              vocx (k,i,j) = (voc0*expo) *(42.08/PMAR)*1.e-6 !VOCS
           ENDDO  !end of looping in k(vertical levels)

        ENDDO  !end of looping in i (longitude)
     ENDDO   !end of looping in j (latitute)

  ENDIF  !end of if in ictrl

  f = 0.2

  !re-distributing sources
  ! Defining bounds
  if (btest(nodeibcon(mynum,ngrid),0)) then !OESTE
     iInit = 3
  else
     iInit = 1
  endif
  if (btest(nodeibcon(mynum,ngrid),1)) then !LESTE
     iEnd = n2-2
  else
     iEnd = n2
  endif
  if (btest(nodeibcon(mynum,ngrid),2)) then !NORTE
     jInit = 3
  else
     jInit = 1
  endif
  if (btest(nodeibcon(mynum,ngrid),3)) then !SUL
     jEnd = n3-2
  else
     jEnd = n3
  endif

  !do j=3,n3-2
     !do i=3,n2-2
  !DO j=1,n3
     !DO i=1,n2
  DO j=jInit,jEnd
     DO i=iInit,iEnd

        DO k=2,n1-1

           no  (k,i,j) = no  (k,i,j) + (1. - f)*nox  (k,i,j)
           no2 (k,i,j) = no2 (k,i,j) + (1. - f)*no2x (k,i,j)
           pm25(k,i,j) = pm25(k,i,j) + (1. - f)*pm25x(k,i,j)
           co  (k,i,j) = co  (k,i,j) + (1. - f)*cox  (k,i,j)
           so2 (k,i,j) = so2 (k,i,j) + (1. - f)*so2x (k,i,j)
           voc (k,i,j) = voc (k,i,j) + (1. - f)*vocx (k,i,j)

           !distribution into the 9 sites around i,j    !  j+1  .   .   .
           !do jj = j-1,j+1                          !   j   .   .   .       
              !do ii = i-1,i+1                        !  j-1  .   .   .
                 !      i-1  i  i+1
           DO jj = MAX(j-1,1),MIN(j+1,n3)
              DO ii = MAX(i-1,1),MIN(i+1,n2)

                 no  (k,ii,jj) = no  (k,ii,jj) + (1./9.)*f*nox  (k,i,j)
                 no2 (k,ii,jj) = no2 (k,ii,jj) + (1./9.)*f*no2x (k,i,j)
                 pm25(k,ii,jj) = pm25(k,ii,jj) + (1./9.)*f*pm25x(k,i,j)
                 co  (k,ii,jj) = co  (k,ii,jj) + (1./9.)*f*cox  (k,i,j)
                 so2 (k,ii,jj) = so2 (k,ii,jj) + (1./9.)*f*so2x (k,i,j)
                 voc (k,ii,jj) = voc (k,ii,jj) + (1./9.)*f*vocx (k,i,j)
                 so4 (k,ii,jj) = so2 (k,ii,jj)*0.2  
                 aer (k,ii,jj) = so2 (k,ii,jj)*0.0 

              ENDDO
           ENDDO

        ENDDO !vertical
     ENDDO  !longitude
  ENDDO   !latitude

!!$  print *, "DEBUG-ALF:init_conc1:mynum,sum(no,no2,pm25,co,so2,voc,so4,aer)=", &
!!$       mynum, sum(no), sum(no2), sum(pm25), sum(co), sum(so2), sum(voc), &
!!$       sum(so4), sum(aer)
!!$  call flush(6)
!!$
!!$  print *, "DEBUG-ALF:init_conc1:mynum,n2,sum(dados(:,1,:)),sum(dados(:,18,:)),sum(dados(:,n2,:))=", &
!!$       mynum, n2, &
!!$       (sum(no(:,1,:))+sum(no2(:,1,:))+sum(pm25(:,1,:))+sum(co(:,1,:))+sum(so2(:,1,:))+sum(voc(:,1,:))+sum(so4(:,1,:))+sum(aer(:,1,:))), &
!!$       (sum(no(:,18,:))+sum(no2(:,18,:))+sum(pm25(:,18,:))+sum(co(:,18,:))+sum(so2(:,18,:))+sum(voc(:,18,:))+sum(so4(:,18,:))+sum(aer(:,18,:))), &
!!$       (sum(no(:,n2,:))+sum(no2(:,n2,:))+sum(pm25(:,n2,:))+sum(co(:,n2,:))+sum(so2(:,n2,:))+sum(voc(:,n2,:))+sum(so4(:,n2,:))+sum(aer(:,n2,:)))
!!$  call flush(6)

END SUBROUTINE init_conc1

!----------------------------------------------------------------

SUBROUTINE init_conc2(ictrl, ngrid, n1, n2, n3, np,    &
     G_URBAN, s7p, s8p, s9p, s10p, s11p, s12p, s13p, zt)

  USE node_mod, ONLY: &
       nodeibcon,     & !INTENT(IN)
       mynum            !INTENT(IN)

  IMPLICIT NONE

  ! Arguments:
  INTEGER, INTENT(in) :: ictrl, ngrid
  INTEGER, INTENT(in) :: n1, n2, n3, np
  REAL, INTENT(in)    :: G_URBAN(n2,n3,np)
  REAL, INTENT(out)   :: s7p(n1,n2,n3)
  REAL, INTENT(out)   :: s8p(n1,n2,n3)
  REAL, INTENT(out)   :: s9p(n1,n2,n3)
  REAL, INTENT(out)   :: s10p(n1,n2,n3)
  REAL, INTENT(out)   :: s11p(n1,n2,n3)
  REAL, INTENT(out)   :: s12p(n1,n2,n3)
  REAL, INTENT(out)   :: s13p(n1,n2,n3)
  REAL, INTENT(in)    :: zt(n1)
  ! Local variables:
  INTEGER         :: i, j, ii, jj, k
  REAL, PARAMETER :: PMAR=28.97
  REAL            :: pfat, f, expo
  REAL            :: s7p0, s8p0, s9p0, s10p0, s11p0, s12p0, s13p0
  REAL            :: s7p2(n1,n2,n3)
  REAL            :: s8p2(n1,n2,n3)
  REAL            :: s9p2(n1,n2,n3)
  REAL            :: s10p2(n1,n2,n3)
  REAL            :: s11p2(n1,n2,n3)
  REAL            :: s12p2(n1,n2,n3)
  REAL            :: s13p2(n1,n2,n3)
  INTEGER         :: jInit, jEnd, iInit, iEnd


  IF (ictrl==1) THEN
     !--------
     !Flag which defines initial condition (IF RUNTYPE = INITIAL)

     s7p(:,:,:)  = 0.			       !O3
     s8p(:,:,:)  = 0.			       !RHCO
     s9p(:,:,:)  = 0.                          !HO2
     s10p(:,:,:) = 0.                          !O3P
     s11p(:,:,:) = 0.                          !O1D
     s12p(:,:,:) = 0.                          !HO
     s13p(:,:,:) = 0.                          !RO2

     DO j=1,n3
        DO i=1,n2

           !identify surface type by using G_URBAN 
           IF (NINT(G_URBAN(i,j,2))/=0.) THEN

              IF (NINT(G_URBAN(i,j,2))==1) THEN
                 pfat=1.
              ENDIF

              IF (NINT(G_URBAN(i,j,2))==2) THEN
                 pfat=0.3     !using only 30% of values for urban regions
              ENDIF

              IF (NINT(G_URBAN(i,j,2))==3) THEN
                 pfat=0.2     !using only 20% of values for urban regions
              ENDIF
           ELSE
              pfat=0.1     
           ENDIF

           ! defining surface backgroud concentrations
           ! (typical values for 00 Z in urban regions)
           ! can be better difined if using real values
           ! (eg. read a concentration file)

           s7p0  = 8.00  *pfat        !O3 (estacao santana dia 31/07/99 ug/m3)
           s8p0  = 0.0187*pfat        !RHCO (valor aproximado segundo Leila,2003, comunicacao pessoal)
           s9p0  = (4. *1.0e-7) *pfat !HO2 ppm
           s10p0 = (4.8*1.0e-10)*pfat !o3p ppm
           s11p0 = (4.8*1.0e-11)*pfat !o1d ppm
           s12p0 = (8. *1.0e-7) *pfat !HO  ppm
           s13p0 = (4. *1.0e-7) *pfat !RO2 ppm

           DO k=2,n1-1
              expo         = EXP(-zt(k)/zt(n1-1))

              s7p2(k,i,j)  = (s7p0 *expo)*(1.0e-9)/1.275      !O3
              s8p2(k,i,j)  = (s8p0 *expo)*(44.05/PMAR)*1.e-6  !RHCO
              s9p2(k,i,j)  = (s9p0 *expo)*(33.0/PMAR) *1.e-6  !HO2
              s10p2(k,i,j) = (s10p0*expo)*(16.0/PMAR) *1.e-6  !O3P
              s11p2(k,i,j) = (s11p0*expo)*(16.0/PMAR) *1.e-6  !O1D
              s12p2(k,i,j) = (s12p0*expo)*(17.0/PMAR) *1.e-6  !HO
              s13p2(k,i,j) = (s13p0*expo)*(47.0/PMAR) *1.e-6  !RO2

           ENDDO

        ENDDO
     ENDDO

  ENDIF

  f=0.2

  ! Defining bounds
  if (btest(nodeibcon(mynum,ngrid),0)) then !OESTE
     iInit = 3
  else
     iInit = 1
  endif
  if (btest(nodeibcon(mynum,ngrid),1)) then !LESTE
     iEnd = n2-2
  else
     iEnd = n2
  endif
  if (btest(nodeibcon(mynum,ngrid),2)) then !NORTE
     jInit = 3
  else
     jInit = 1
  endif
  if (btest(nodeibcon(mynum,ngrid),3)) then !SUL
     jEnd = n3-2
  else
     jEnd = n3
  endif

  !do j=3,n3-2
     !do i=3,n2-2
  !DO j=1,n3
     !DO i=1,n2
  DO j=jInit,jEnd
     DO i=iInit,iEnd

        DO k=2,n1-1

           s7p(k,i,j)  =  s7p(k,i,j) + (1.-f)*s7p2(k,i,j)
           s8p(k,i,j)  =  s8p(k,i,j) + (1.-f)*s8p2(k,i,j)
           s9p(k,i,j)  =  s9p(k,i,j) + (1.-f)*s9p2(k,i,j)
           s10p(k,i,j) = s10p(k,i,j) + (1.-f)*s10p2(k,i,j)
           s11p(k,i,j) = s11p(k,i,j) + (1.-f)*s11p2(k,i,j)
           s12p(k,i,j) = s12p(k,i,j) + (1.-f)*s12p2(k,i,j)
           s13p(k,i,j) = s13p(k,i,j) + (1.-f)*s13p2(k,i,j)

           !distribution into the 9 sites around i,j    !  j+1  .   .   .
           !do jj = j-1,j+1                          !   j   .   .   .       
              !do ii = i-1,i+1                        !  j-1  .   .   .
                 !      i-1  i  i+1
           DO jj = MAX(j-1,1),MIN(j+1,n3)
              DO ii = MAX(i-1,1),MIN(i+1,n2)

                 s7p(k,ii,jj)  =  s7p(k,ii,jj)+(1./9.)*f* s7p2(k,i,j)
                 s8p(k,ii,jj)  =  s8p(k,ii,jj)+(1./9.)*f* s8p2(k,i,j)
                 s9p(k,ii,jj)  =  s9p(k,ii,jj)+(1./9.)*f* s9p2(k,i,j)
                 s10p(k,ii,jj) = s10p(k,ii,jj)+(1./9.)*f*s10p2(k,i,j)
                 s11p(k,ii,jj) = s11p(k,ii,jj)+(1./9.)*f*s11p2(k,i,j)
                 s12p(k,ii,jj) = s12p(k,ii,jj)+(1./9.)*f*s12p2(k,i,j)
                 s13p(k,ii,jj) = s13p(k,ii,jj)+(1./9.)*f*s13p2(k,i,j)

              ENDDO
           ENDDO

        ENDDO

     ENDDO
  ENDDO

!!$  print *, "DEBUG-ALF:init_conc2:mynum,sum(s7p,s8p,s9p,s10p,s11p,s12p,s13p)=",&
!!$       mynum, sum(s7p), sum(s8p), sum(s9p), sum(s10p), sum(s11p), sum(s12p), &
!!$       sum(s13p)
!!$  call flush(6)
!!$
!!$  print *, "DEBUG-ALF:init_conc2:mynum,n2,sum(dados(:,1,:)),sum(dados(:,18,:)),sum(dados(:,n2,:))=", &
!!$       mynum, n2, &
!!$       (sum(s7p(:,1,:))+sum(s8p(:,1,:))+sum(s9p(:,1,:))+sum(s10p(:,1,:))+sum(s11p(:,1,:))+sum(s12p(:,1,:))+sum(s13p(:,1,:)) ), &
!!$       (sum(s7p(:,18,:))+sum(s8p(:,18,:))+sum(s9p(:,18,:))+sum(s10p(:,18,:))+sum(s11p(:,18,:))+sum(s12p(:,18,:))+sum(s13p(:,18,:)) ), &
!!$       (sum(s7p(:,n2,:))+sum(s8p(:,n2,:))+sum(s9p(:,n2,:))+sum(s10p(:,n2,:))+sum(s11p(:,n2,:))+sum(s12p(:,n2,:))+sum(s13p(:,n2,:)) )
!!$  call flush(6)


END SUBROUTINE init_conc2

!----------------------------------------------------------------

SUBROUTINE EMFACTOR(dayin, idays, dayout)
  IMPLICIT NONE
  ! Arguments:
  CHARACTER(len=3), INTENT(IN)  :: dayin
  INTEGER, intent(IN)           :: idays
  CHARACTER(len=3), INTENT(OUT) :: dayout
  ! Local Variables:
  INTEGER ::i, j, k, id, ndays, irest
  CHARACTER(len=3) :: cday(7)

  ! Init array: cday
  cday(:) = (/'SUN', 'MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT'/)
!!$  DATA (cday(id),id=1,7) /'SUN', 'MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT'/
!!$  !                                               0      1      2      3   
!!$  !                         4      5      6      7      8      9      10
!!$  !                         11     12     13     14     15     16     17
!!$  !                         18     19     20     21     22     23     24

  ndays = (idays/7)*7    ! (21)
  irest = idays - ndays  !  (3)

  DO i=1,7
     IF (dayin==cday(i)) j = i
  ENDDO

  k = j + irest

  IF (k>7) k = k-7

  dayout = cday(k)

END SUBROUTINE EMFACTOR

!###########################################################################

SUBROUTINE init_conc_prev()

  ! This routine initializes gas variables from a previous (day-1) history file

  USE mem_grid, ONLY: &
       maxgrds,       & ! INTENT(IN)
       ngridsh          ! INTENT(OUT)
  USE mem_emiss, ONLY: &
       chemdata_in      ! INTENT(IN)
  USE node_mod, ONLY: &
       mchnum,        & ! INTENT(IN)
       master_num       ! INTENT(IN)
  USE ParLib, ONLY: &
       parf_bcast       ! Subroutine
  IMPLICIT NONE
  INCLUDE "i8.h"
  include "files.h"

  ! Local variables:
  INTEGER :: ngrids1, ioutput1,  &
       nnxp1(maxgrds), nnyp1(maxgrds), nnzp1(maxgrds), nzg1, nzs1, npatch1
  INTEGER :: ie, maxarr, ngr, nc
  CHARACTER (len=f_name_length) :: hnameinh
  INTEGER, EXTERNAL  :: cio_i
  INTEGER, PARAMETER :: iunhd=11


  ! Open the input history header file and read some of the info.

  IF (mchnum==master_num) THEN !io-process only

     WRITE(*,*)'chemdata_in=',chemdata_in
     !pause

     nc       = LEN_TRIM(chemdata_in)
     hnameinh = chemdata_in(1:nc-9)//'.vfm'

     CALL rams_f_open(iunhd, chemdata_in(1:len_trim(chemdata_in)), 'FORMATTED', 'OLD', 'READ', 0)

     ie      = cio_i(iunhd, 1, 'ngrids',  ngrids1,  1)
     ngridsh = ngrids1
     ie      = cio_i(iunhd, 1, 'nnxp',    nnxp1,    ngrids1)
     ie      = cio_i(iunhd, 1, 'nnyp',    nnyp1,    ngrids1)
     ie      = cio_i(iunhd, 1, 'nnzp',    nnzp1,    ngrids1)
     ie      = cio_i(iunhd, 1, 'npatch',  npatch1,  1)
     ie      = cio_i(iunhd, 1, 'nzg',     nzg1,     1)
     ie      = cio_i(iunhd, 1, 'nzs',     nzs1,     1)
     ie      = cio_i(iunhd, 1, 'ioutput', ioutput1, 1)

     ! Find maximum size of any array on history file.
     ! Allocate scratch array of this size.

     maxarr = 0
     DO ngr=1,ngridsh
        maxarr = MAX(                            &
             maxarr,                             &
             nnxp1(ngr)*nnyp1(ngr)*nnzp1(ngr),   &
             nnxp1(ngr)*nnyp1(ngr)*nzg1*npatch1, &
             nnxp1(ngr)*nnyp1(ngr)*nzs1*npatch1  )
     ENDDO

  ENDIF

  ! Broadcasting Data
  CALL parf_bcast(ngrids1, master_num)
  CALL parf_bcast(nnxp1, INT(SIZE(nnxp1,1),i8), master_num)
  CALL parf_bcast(nnyp1, INT(SIZE(nnyp1,1),i8), master_num)
  CALL parf_bcast(nnzp1, INT(SIZE(nnzp1,1),i8), master_num)
  CALL parf_bcast(npatch1, master_num)
  CALL parf_bcast(nzg1, master_num)
  CALL parf_bcast(nzs1, master_num)
  CALL parf_bcast(ioutput1, master_num)
  CALL parf_bcast(maxarr, master_num)

  ! read stuff here

  CALL hist_pol_read(maxarr, hnameinh(1:len_trim(hnameinh)), iunhd)

  IF (mchnum==master_num) THEN !io-process only
     PRINT*, 'back from read'
     CLOSE(iunhd)
  ENDIF


END SUBROUTINE init_conc_prev


!******************************************************************************

SUBROUTINE hist_pol_read(maxarr, hnamein, iunhd)

  USE an_header, ONLY: &
       head_table,     &  ! Type
       nvbtab             ! INTENT(OUT)
  USE var_tables, ONLY: &
       nvgrids,         & ! INTENT(IN)
       num_var,         & ! INTENT(IN)
       vtab_r             ! INTENT(INOUT)
  USE node_mod, ONLY: &
       mchnum,        & ! INTENT(IN)
       master_num       ! INTENT(IN)
  USE ParLib, ONLY: &
       parf_bcast       ! Subroutine

  IMPLICIT NONE
  INCLUDE "i8.h"
  include "files.h"

  ! Arguments:
  INTEGER, INTENT(in)           :: maxarr
  CHARACTER(len=f_name_length), INTENT(in) :: hnamein
  INTEGER, INTENT(in)           :: iunhd
  ! Local variables:
  INTEGER            :: ngr, npts, nptsh, nv, nvh, i
  REAL, ALLOCATABLE  :: scr(:)
  INTEGER, PARAMETER :: inhunt=10
  TYPE(head_table), ALLOCATABLE :: hr_table(:)
  INTEGER            :: iCopyFlg


  !  Read variable header info
  IF (mchnum==master_num) THEN !io-process only
     ALLOCATE (scr(maxarr))
     REWIND(iunhd)
     READ(iunhd,*) nvbtab
     ALLOCATE (hr_table(nvbtab))
     DO nv=1,nvbtab
        READ(iunhd,*)                &
             hr_table(nv)%string,    &
             hr_table(nv)%npointer,  &
             hr_table(nv)%idim_type, &
             hr_table(nv)%ngrid,     &
             hr_table(nv)%nvalues
     ENDDO

  ENDIF

  ! Broadcasting data
  CALL parf_bcast(nvbtab, master_num)

  IF (mchnum==master_num) THEN !io-process only
     ! Open history data file
     CALL rams_f_open(inhunt, hnamein(1:len_trim(hnamein)), 'UNFORMATTED', 'OLD', 'READ', 0)
  ENDIF

  ! Loop through all variables
  DO nvh=1,nvbtab

     IF (mchnum==master_num) THEN !io-process only
        ! Read a variable
        nptsh = hr_table(nvh)%nvalues
        READ(inhunt) (scr(i), i=1,nptsh)
        
        !  See if this variable is active in the current run
        ngr = hr_table(nvh)%ngrid
        IF (ngr>nvgrids) CYCLE
        
        DO nv=1,num_var(ngr)
           npts = vtab_r(nv,ngr)%npts
           IF (hr_table(nvh)%string==vtab_r(nv,ngr)%name) THEN
              IF (nptsh/=npts) THEN
                 PRINT*, 'Grid point number mismatch on history field:',  &
                      vtab_r(nv,ngr)%name,npts,nptsh
                 CALL fatal_error('History read number points error')
              ENDIF

              iCopyFlg = 0
              
              IF(  vtab_r(nv,ngr)%name=='PNO'   .OR. &
                   vtab_r(nv,ngr)%name=='PNO2'  .OR. &
                   vtab_r(nv,ngr)%name=='PPM25' .OR. &
                   vtab_r(nv,ngr)%name=='PCO'   .OR. &
                   vtab_r(nv,ngr)%name=='PVOC'  .OR. &
                   vtab_r(nv,ngr)%name=='PSO2'  .OR. &
                   vtab_r(nv,ngr)%name=='PSO4'  .OR. &
                   vtab_r(nv,ngr)%name=='PAER'  .OR. &
                   vtab_r(nv,ngr)%name=='PVOC'  .OR. &
                   vtab_r(nv,ngr)%name=='PSO2'  .OR. &
                   vtab_r(nv,ngr)%name=='PO3'   .OR. &
                   vtab_r(nv,ngr)%name=='PRHCO' .OR. &
                   vtab_r(nv,ngr)%name=='PHO2'  .OR. &
                   vtab_r(nv,ngr)%name=='PO3P'  .OR. &
                   vtab_r(nv,ngr)%name=='PO1D'  .OR. &
                   vtab_r(nv,ngr)%name=='PHO'   .OR. &
                   vtab_r(nv,ngr)%name=='PROO'       ) THEN
                 
                 WRITE (UNIT=6, FMT='(a25,2i5,3x,a18,i10)') &
                      'Polutants History filling grid: ',   &
                      ngr, nv, vtab_r(nv,ngr)%name, npts
                 !	 pause
                 iCopyFlg = 1
!!$                 call atob(npts, scr(1), vtab_r(nv,ngr)%var_p)

                 EXIT
              ENDIF

           ENDIF
        ENDDO

     ENDIF

     ! Broadcasting last array read
     IF (iCopyFlg==1) THEN
        CALL parf_bcast(nv, master_num)
        CALL parf_bcast(ngr, master_num)
        CALL parf_bcast(npts, master_num)
        CALL parf_bcast(scr, INT(npts, i8), master_num)

        CALL atob(npts, scr(1), vtab_r(nv,ngr)%var_p)
     ENDIF

  ENDDO
  
  IF (mchnum==master_num) THEN !io-process only
     ! Close the input history file
     CLOSE(inhunt)
     DEALLOCATE(scr, hr_table)
  ENDIF

END SUBROUTINE hist_pol_read


!!$!-------------------------------
!!$!DEBUG-ALF
!!$subroutine debug_gas()
!!$
!!$  use mem_gaspart
!!$
!!$  implicit none
!!$
!!$  print *, "DEBUG-ALF:debug_gas:sum(so2)=", sum(gaspart_g(1)%pso2(:,:,:))
!!$
!!$
!!$end subroutine debug_gas
