MODULE mod_chem_orage

  USE mem_scratch1_grell, ONLY: &
        ierr4d                     ! (IN)

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: eclair_driver ! Subroutine


CONTAINS

  SUBROUTINE eclair_driver(mzp,mxp,myp,ia,iz,ja,jz,dtlt,rtgt,dxt,dyt,dzt,zt, &
                           nzpmax,mgmxp,mgmyp,maxiens,nnqparm,wp,rtp,rv, &
                           pp,pi0,theta,sc_t,weight,nspecies,no,cp,cpor)

    ! node_mod
    INTEGER , INTENT(IN) :: mzp
    INTEGER , INTENT(IN) :: mxp  
    INTEGER , INTENT(IN) :: myp
    INTEGER , INTENT(IN) :: ia
    INTEGER , INTENT(IN) :: iz
    INTEGER , INTENT(IN) :: ja
    INTEGER , INTENT(IN) :: jz

    ! mem_grid
    REAL    , INTENT(IN) :: dtlt
    REAL    , INTENT(IN) :: rtgt(mxp,myp)
    REAL    , INTENT(IN) :: dxt(mxp,myp)
    REAL    , INTENT(IN) :: dyt(mxp,myp)
    REAL    , INTENT(IN) :: dzt(nzpmax)
    REAL    , INTENT(IN) :: zt(nzpmax)

    ! grid_dims
    INTEGER , INTENT(IN) :: nzpmax

    ! mem_grell_param
    INTEGER , INTENT(IN) :: mgmxp
    INTEGER , INTENT(IN) :: mgmyp
    INTEGER , INTENT(IN) :: maxiens

!    ! mem_scratch1_grell
!    INTEGER , INTENT(IN) :: ierr4d(mgmxp,mgmyp,maxiens)

    ! mem_cuparm
    INTEGER , INTENT(IN) :: nnqparm

    ! mem_basic
    REAL    , INTENT(IN)    :: wp(mzp,mxp,myp)
    REAL    , INTENT(IN)    :: rtp(mzp,mxp,myp)
    REAL    , INTENT(IN)    :: rv(mzp,mxp,myp)
    REAL    , INTENT(IN)    :: pp(mzp,mxp,myp)
    REAL    , INTENT(IN)    :: pi0(mzp,mxp,myp)
    REAL    , INTENT(IN)    :: theta(mzp,mxp,myp)

    ! mem_chem1
    REAL    , INTENT(INOUT) :: sc_t(mzp*mxp*myp)

    ! chem1_list
    REAL    , INTENT(IN)    :: weight(nspecies)
    INTEGER , INTENT(IN)    :: nspecies
    INTEGER , INTENT(IN)    :: no

    ! rconstants
    REAL, INTENT(IN) :: cp
    REAL, INTENT(IN) :: cpor

    INTEGER ierr(mxp,myp)

    !  if(.NOT. allocated(chem1_g(NO)%sc_t)) return

    ! general control for activation or not of the NO production by 'eclair'
    ! only for Grell scheme and deep convection
    IF( NNQPARM == 2) THEN
       ierr(1:mxp,1:myp) = ierr4d(1:mxp,1:myp,1,1)
    ELSE
       RETURN
       !ierr(1:mxp,1:myp) = 0
    ENDIF
    CALL eclair(mzp,mxp,myp,ia,iz,ja,jz,sc_t(1),wp(1,1,1),rtp(1,1,1),rv(1,1,1),pp(1,1,1), & 
                pi0(1,1,1),theta(1,1,1),rtgt(1,1),dxt(1,1),dyt(1,1),dzt,zt,dtlt,ierr, &
                weight,nspecies,no,cp,cpor)

  END SUBROUTINE eclair_driver
  !------------------------------------------------------------------------------------

  SUBROUTINE eclair(mzp,mxp,myp,ia,iz,ja,jz,tend_NO,wp,rtp,rv,pp,pi0,theta,rtgt,dxt, &
                    dyt,dzt,zt,dtlt,ierr,weight,nspecies,no,cp,cpor) 		      

    INTEGER , INTENT(IN)    :: mzp
    INTEGER , INTENT(IN)    :: mxp
    INTEGER , INTENT(IN)    :: myp
    INTEGER , INTENT(IN)    :: ia
    INTEGER , INTENT(IN)    :: iz
    INTEGER , INTENT(IN)    :: ja
    INTEGER , INTENT(IN)    :: jz
    REAL    , INTENT(IN)    :: rtgt(mxp,myp)
    REAL    , INTENT(IN)    :: dxt(mxp,myp)
    REAL    , INTENT(IN)    :: dyt(mxp,myp)
    REAL    , INTENT(IN)    :: zt(mzp)
    REAL    , INTENT(IN)    :: dzt(mzp)

    REAL    , INTENT(IN)    :: wp(mzp,mxp,myp)
    REAL    , INTENT(IN)    :: rtp(mzp,mxp,myp)
    REAL    , INTENT(IN)    :: rv(mzp,mxp,myp) 
    REAL    , INTENT(IN)    :: pp(mzp,mxp,myp)
    REAL    , INTENT(IN)    :: pi0(mzp,mxp,myp)
    REAL    , INTENT(IN)    :: theta(mzp,mxp,myp)
    REAL    , INTENT(INOUT) :: tend_NO(mzp,mxp,myp)

    INTEGER , INTENT(IN)    :: ierr(mxp,myp)
    REAL    , INTENT(IN)    :: dtlt
    REAL    , INTENT(IN)    :: weight(nspecies)
    INTEGER , INTENT(IN)    :: nspecies
    INTEGER , INTENT(IN)    :: no

    ! rconstants
    REAL, INTENT(IN) :: cp
    REAL, INTENT(IN) :: cpor

    REAL, DIMENSION(mzp) :: w,xtempc,rtot,rvap,NOIC,NOGC,xpress
    REAL, PARAMETER      :: pmar=28.96  
    REAL                 :: v,ycunit,emi,xtemp
    INTEGER              :: i,j,k  
    !
    !
    DO j = ja,jz
       DO i = ia,iz

          IF(ierr(i,j) > 0) CYCLE ! if there is not convection, cycle

          DO k = 2, mzp-1  
             xtempc(k) = theta(k,i,j)*(pi0(k,i,j)+pp(k,i,j))/cp
             xtempc(k) = xtempc(k) - 273.16  
             rtot  (k) = rtp(k,i,j)
             rvap  (k) = rv(k,i,j)
             w     (k) = wp(k,i,j)
             xpress(k) = (((pi0(k,i,j)+pp(k,i,j))/cp)**cpor)*1.e5
          END DO

          CALL orage_NO(mzp,w,rtot,rvap,zt,rtgt(i,j),xtempc,NOIC,NOGC)

          ! conversion en molecules par cm3
          DO k = 2, mzp-1
             !srf    v= deltax*deltay*(zt(k+1)-zt(k)) * 1.e6
             v= (1./(dxt(i,j)*dyt(i,j)*dzt(k))*rtgt(i,j)) * 1.e6

             NOIC(k) = NOIC(k) /V     
             NOGC(k) = NOGC(k) /V     
             ! conversion en ppbm /sec
             xtemp = xtempc(k) + 273.16
             ycunit =(6.02e23*1e-15*xpress(k)/(8.314e0*xtemp)*(pmar/weight(NO)))
             !lpce    ycunit = 2.45e13*(298.15/xtemp)*(xpress(k)/101325.)
             NOIC(k) = NOIC(k)/ycunit 
             NOGC(k) = NOGC(k)/ycunit   
             NOIC(k) = NOIC(k)/60.
             NOGC(k) = NOGC(k)/60.
             emi = (NOIC(k)+NOGC(k)) 
             !ajout emission a chaque pas de temps dtlt
             ! la parametrisation des taux de NOX est calculee pour 1min soit 60s
             ! on calcule a chaque pas de temps et on fait la proportion par rapport au temps.
             ! NOX(k,i,j)     = NOX(k,i,j)    + emi*dtlt
             !srf-   converting to tendency
             !the parameterization for the NOx rate is calculated for 1 minutes so you
             !have to take into account the fact that the timestep you are using
             !is not 60 s so that you need to make a ratio :
             !noxrate(1 time step) = noxrate (1 minute)*timestep/60.
             tend_NO(k,i,j) = tend_NO(k,i,j) + emi*dtlt/60.
          END DO
       END DO
    END DO
  END SUBROUTINE Eclair

  !--------------------------------------------------------
  SUBROUTINE orage_NO(mzp,w,rtot,rvap,zt,rtgt,xtempc,noic,nogc)
  !--------------------------------------------------------

    INTEGER , INTENT(IN)    :: mzp
    REAL    , INTENT(IN)    :: w(mzp)
    REAL    , INTENT(IN)    :: rtot(mzp)
    REAL    , INTENT(IN)    :: rvap(mzp)
    REAL    , INTENT(IN)    :: zt(mzp)
    REAL    , INTENT(IN)    :: xtempc(mzp)
    REAL    , INTENT(IN)    :: rtgt
    REAL    , INTENT(INOUT) :: noic(mzp) ! (DMK) alterado (OUT) para (INOUT)
    REAL    , INTENT(INOUT) :: nogc(mzp) ! (DMK) alterado (OUT) para (INOUT)

    INTEGER :: ick1,ick2,ick3
    REAL    :: wmax,rnug,f,xr,dh,z,p,xn1,xn2,fic,fgc

    INTEGER :: k

    ! le calcul est effectue sur une colonne dans une maille i,j
    ! pour un nuage en developpement (w >0.)
    ! Pickering et al,1998 JGR
    ! F : Nombre total de flashs/mn  
    ! Zone IC : phase glace uniquement 
    ! Zone GC : phase  mixte rtot-rvap > 10-3 g/kg soit > 20 dBZ = 10-6 kg/kg
    ! ICK1 : altitude de iso 0C
    ! ICK2 : sommet de zone GC temp > -15C
    ! ICK3 : sommet de zone IC rtp > 10-5 kg/kg
    ! P : fraction de flashs dans GC

    ! calcul de la vitesse max dans la colonne
    wmax = 0.
    rnug=0.
    DO k = 2, mzp-1
       NOIC(k) = 0.
       NOGC(k)=0.
       IF(w(k).GT.wmax) THEN
          wmax = w(k)
       END IF
       ! calcul de l eau condensee max
       xr = rtot(k) - rvap(k)
       IF(xr.GT.rnug) THEN
          rnug = xr
       END IF
    END DO
    IF(wmax.LE.0.) THEN
       RETURN
    END IF
    IF(rnug.LT.1.e-6) THEN
       RETURN
    END IF
    !  print*, 'je suis dans orage_no'  

    ! Calcul de F
    F = 5.E-6*(wmax**4.54 )
    !
    ! determination des zones IC et GC
    ICK1 = 0
    ICK2 = 0
    ICK3 = 0
    DO K = 1, mzp
       ! altitude de l iso 0
       IF (xtempc(k).GE.0.) THEN
          ICK1 = K
       END IF
       ! altitude de l iso -15
       IF (xtempc(k).GE.-15.) THEN
          ICK2 = K
       END IF
       ! sommet zone IC
       xr = rtot(k) - rvap(k)
       IF (xr.GE.1.e-5) THEN
          ICK3 = K
       END IF
    END DO
    !   print*, 'ICK1=', ICK1, 'ICK2=', ICK2, 'ICK3=',ICK3
    !  pas de possibilite a differencier zone IC et GC
    IF (ICK3.LE.ICK2) THEN
       RETURN
    END IF
    ! calcul de la proportion de flashs dans GC
    ! DH en km
    !srf  DH = (zt(ICK3)-zt(ICK1))/1000.0
    DH = (zt(ICK3)-zt(ICK1))*rtgt/1000.0
    !  print*, 'DH=',DH
    z = 0.02*(DH**4) - 0.648*(DH**3) + 7.493*(DH**2) &
         -36.54*DH + 63.09
    P = 1./(z+1.)
    ! intervalle en metre des regions IC et GC
    !srf  xn1 = ZT(ICK3)-ZT(ICK2)
    xn1 = (ZT(ICK3)-ZT(ICK2))*rtgt
    !srf  xn2 = ZT(ICK2)
    xn2 = ZT(ICK2)*rtgt
    !  print*, 'xn1=',xn1,'xn2=',xn2
    ! nombre de flashs dans IC par metre et par min
    FIC = F*(1.-P)/xn1
    ! nombre de flashs dans GC par metre et par min
    FGC = F*P/xn2
    !  print*, 'fic=',fic,'fgc=',fgc
    ! calcul production de NO
    DO k = 2, mzp-1
       ! production de molecules de NO pour zone IC  /min
       IF(k.GE.ICK2.AND.k.LT.ICK3) THEN
          !srf      NOIC(K) = 6.7e25 * FIC* (zt(k+1) - zt(k))
          NOIC(K) = 6.7e25 * FIC* (zt(k+1) - zt(k))*rtgt
          !  print*, 'NOIC=',NOIC(K),'k=', k
       END IF
       ! production de molecules de NO pour zone GC  /min
       xr = rtot(k) - rvap(k)
       IF(xr.GE.1.e-6) THEN
          IF (k.LT.ICK2) THEN
             !srf        NOGC(K) = 6.7e26 * FGC * (zt(k+1) - zt(k))
             NOGC(K) = 6.7e26 * FGC * (zt(k+1) - zt(k))*rtgt
             !  print*, 'NOGC=',NOGC(K),'k=', k
          END IF
       END IF
    END DO

  END SUBROUTINE orage_NO


END MODULE mod_chem_orage
