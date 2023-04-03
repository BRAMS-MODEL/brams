!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE ISRP1F
!! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FOREWARD PROBLEM OF
!!     AN AMMONIUM-SULFATE AEROSOL SYSTEM.
!!     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY
!!     THE AMBIENT RELATIVE HUMIDITY.
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!  REVISION HISTORY:                                                   *
!!   Modified by Yang Zhang of AER in Sept., 2001 to incorporate the CMU*
!!           hybrid mass transfer module into MADRID1                   *
!!           The code for the CMU hybrid mass transfer module was       *
!!           provided by Bonyoung Koo/Spyros Pandis of CMU on           *
!!           May 18, 2001                                               *
!!   Modified by Betty Pun, AER, February/March 2002 to work with       *
!!           organic anion and water associated with hydrophilic        *
!!           organic SOA (MADRID 2)                                     *
!!***********************************************************************
SUBROUTINE isrp1f (wi, rhi, tempi)
   USE Isorropia_Module, ONLY: sulrat,w,zero,rh,drnh42s4,drnh4hs4,drlc, &
                               ncomp,metstbl,scase,init1,calcnh3
   implicit none
      
   DOUBLE precision, INTENT(IN OUT)                     :: wi(ncomp)
   DOUBLE precision, INTENT(IN OUT)                     :: rhi
   DOUBLE precision, INTENT(IN OUT)                     :: tempi
   double precision :: dc
   
   ! *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
   
   CALL init1 (wi, rhi, tempi)
   
   ! *** CALCULATE SULFATE RATIO *******************************************
   
   sulrat = w(3)/w(2)
   
   ! *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
   
   ! *** SULFATE POOR
   
   IF (2.0 <= sulrat) THEN
     dc   = w(3) - 2.001D0*w(2)  ! For numerical stability
     w(3) = w(3) + MAX(-dc, zero)
     
     IF(metstbl == 1) THEN
       scase = 'A2'
       CALL calca2                 ! Only liquid (metastable)
     ELSE
       
       IF (rh < drnh42s4) THEN
         scase = 'A1'
         CALL calca1              ! NH42SO4              ; case A1
         
       ELSE IF (drnh42s4 <= rh) THEN
         scase = 'A2'
         CALL calca2              ! Only liquid          ; case A2
       END IF
     END IF
     
   ! *** SULFATE RICH (NO ACID)
     
   ELSE IF (1.0 <= sulrat .AND. sulrat < 2.0) THEN
     
     IF(metstbl == 1) THEN
       scase = 'B4'
       CALL calcb4                 ! Only liquid (metastable)
     ELSE
       
       IF (rh < drnh4hs4) THEN
         scase = 'B1'
         CALL calcb1              ! NH4HSO4,LC,NH42SO4   ; case B1
         
       ELSE IF (drnh4hs4 <= rh .AND. rh < drlc) THEN
         scase = 'B2'
         CALL calcb2              ! LC,NH42S4            ; case B2
         
       ELSE IF (drlc <= rh .AND. rh < drnh42s4) THEN
         scase = 'B3'
         CALL calcb3              ! NH42S4               ; case B3
         
       ELSE IF (drnh42s4 <= rh) THEN
         scase = 'B4'
         CALL calcb4              ! Only liquid          ; case B4
       END IF
     END IF
     CALL calcnh3
     
   ! *** SULFATE RICH (FREE ACID)
     
   ELSE IF (sulrat < 1.0) THEN
     
     IF(metstbl == 1) THEN
       scase = 'C2'
       CALL calcc2                 ! Only liquid (metastable)
     ELSE
       
       IF (rh < drnh4hs4) THEN
         scase = 'C1'
         CALL calcc1              ! NH4HSO4              ; case C1
         
       ELSE IF (drnh4hs4 <= rh) THEN
         scase = 'C2'
         CALL calcc2              ! Only liquid          ; case C2
         
       END IF
     END IF
     CALL calcnh3
   END IF
   
END SUBROUTINE  isrp1f

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE ISRP2F
!! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FOREWARD PROBLEM OF
!!     AN AMMONIUM-SULFATE-NITRATE AEROSOL SYSTEM.
!!     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY
!!     THE AMBIENT RELATIVE HUMIDITY.
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE isrp2f (wi, rhi, tempi)
   USE Isorropia_Module, ONLY: sulrat,w,metstbl,scase,rh,drnh4no3,drnh42s4, &
                               drnh4hs4,drlc,ncomp,init2,CalcNa,water
   IMPLICIT NONE

   DOUBLE precision, INTENT(IN OUT)                     :: wi(ncomp)
   DOUBLE precision, INTENT(IN OUT)                     :: rhi
   DOUBLE precision, INTENT(IN OUT)                     :: tempi
   
   
   ! *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
   !print*,"CASE0=",scase,water
   
   CALL init2 (wi, rhi, tempi)
   !print*,"CASE1=",scase,water
   
   ! *** CALCULATE SULFATE RATIO *******************************************
   
   sulrat = w(3)/w(2)
   
   ! *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
   
   ! *** SULFATE POOR
   
   IF (2.0 <= sulrat) THEN
     
     IF(metstbl == 1) THEN
       scase = 'D3'
       CALL calcd3                 ! Only liquid (metastable)
     ELSE
       
       IF (rh < drnh4no3) THEN
         scase = 'D1'
         CALL calcd1              ! NH42SO4,NH4NO3       ; case D1
         
       ELSE IF (drnh4no3 <= rh .AND. rh < drnh42s4) THEN
         scase = 'D2'
         CALL calcd2              ! NH42S4               ; case D2
         
       ELSE IF (drnh42s4 <= rh) THEN
         scase = 'D3'
         CALL calcd3              ! Only liquid          ; case D3
       END IF
     END IF
     
   ! *** SULFATE RICH (NO ACID)
   !     FOR SOLVING THIS CASE, NITRIC ACID IS ASSUMED A MINOR SPECIES,
   !     THAT DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM.
   !  SUBROUTINES CALCB? ARE CALLED, AND THEN THE NITRIC ACID IS DISSOLVED
   !     FROM THE HNO3(G) -> (H+) + (NO3-) EQUILIBRIUM.
     
   ELSE IF (1.0 <= sulrat .AND. sulrat < 2.0) THEN
     
     IF(metstbl == 1) THEN
       scase = 'B4'
       CALL calcb4                 ! Only liquid (metastable)
       scase = 'E4'
     ELSE
       
       IF (rh < drnh4hs4) THEN
         scase = 'B1'
         CALL calcb1              ! NH4HSO4,LC,NH42SO4   ; case E1
         scase = 'E1'
         
       ELSE IF (drnh4hs4 <= rh .AND. rh < drlc) THEN
         scase = 'B2'
         CALL calcb2              ! LC,NH42S4            ; case E2
         scase = 'E2'
         
       ELSE IF (drlc <= rh .AND. rh < drnh42s4) THEN
         scase = 'B3'
         CALL calcb3              ! NH42S4               ; case E3
         scase = 'E3'
         
       ELSE IF (drnh42s4 <= rh) THEN
         scase = 'B4'
         CALL calcb4              ! Only liquid          ; case E4
         scase = 'E4'
       END IF
     END IF
     
     CALL calcna                 ! HNO3(g) DISSOLUTION
     
   ! *** SULFATE RICH (FREE ACID)
   !     FOR SOLVING THIS CASE, NITRIC ACID IS ASSUMED A MINOR SPECIES,
   !     THAT DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM
   !  SUBROUTINE CALCC? IS CALLED, AND THEN THE NITRIC ACID IS DISSOLVED
   !     FROM THE HNO3(G) -> (H+) + (NO3-) EQUILIBRIUM.
     
   ELSE IF (sulrat < 1.0) THEN
     
     IF(metstbl == 1) THEN
       scase = 'C2'
       CALL calcc2                 ! Only liquid (metastable)
       scase = 'F2'
     ELSE
       
       IF (rh < drnh4hs4) THEN
         scase = 'C1'
         CALL calcc1              ! NH4HSO4              ; case F1
         scase = 'F1'
         
       ELSE IF (drnh4hs4 <= rh) THEN
         scase = 'C2'
         CALL calcc2              ! Only liquid          ; case F2
         scase = 'F2'
       END IF
     END IF
     
     CALL calcna                 ! HNO3(g) DISSOLUTION
   END IF
   !print*,"CASE2=",scase,water
END SUBROUTINE  isrp2f

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE ISRP3F
!! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FORWARD PROBLEM OF
!!     AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM AEROSOL SYSTEM.
!!     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE & SODIUM
!!     RATIOS AND BY THE AMBIENT RELATIVE HUMIDITY.
!! *** COPYRIGHT 1996-2000 UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE isrp3f (wi, rhi, tempi)
   USE Isorropia_Module
   IMPLICIT NONE
   
   double precision, INTENT(OUT)                        :: wi(ncomp)
   double precision, INTENT(IN OUT)                     :: rhi
   double precision, INTENT(IN OUT)                     :: tempi
   
   double precision :: rest
   
   ! *** ADJUST FOR TOO LITTLE AMMONIUM AND CHLORIDE ***********************
   
   wi(3) = MAX (wi(3), 1.d-10)  ! NH4+ : 1e-4 umoles/m3
   wi(5) = MAX (wi(5), 1.d-10)  ! Cl-  : 1e-4 umoles/m3
   
   ! *** ADJUST FOR TOO LITTLE SODIUM, SULFATE AND NITRATE COMBINED ********
   
   IF (wi(1)+wi(2)+wi(4) <= 1D-10) THEN
     wi(1) = 1.d-10  ! Na+  : 1e-4 umoles/m3
     wi(2) = 1.d-10  ! SO4- : 1e-4 umoles/m3
   END IF
   
   ! *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
   
   CALL isoinit3 (wi, rhi, tempi)
   
   ! *** CHECK IF TOO MUCH SODIUM ; ADJUST AND ISSUE ERROR MESSAGE *********
   
   rest = 2.d0*w(2) + w(4) + w(5)
   IF (w(1) > rest) THEN            ! NA > 2*SO4+CL+NO3 ?
     w(1) = (one-1D-6)*rest         ! Adjust Na amount
     CALL pusherr (0050, 'ISRP3F')  ! Warning error: Na adjusted
   END IF
   
   ! *** CALCULATE SULFATE & SODIUM RATIOS *********************************
   
   sulrat = (w(1)+w(3))/w(2)
   sodrat = w(1)/w(2)
   
   ! *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
   
   ! *** SULFATE POOR ; SODIUM POOR
   
   IF (2.0 <= sulrat .AND. sodrat < 2.0) THEN
     
     IF(metstbl == 1) THEN
       scase = 'G5'
       CALL calcg5                 ! Only liquid (metastable)
     ELSE
       
       IF (rh < drnh4no3) THEN
         scase = 'G1'
         CALL calcg1              ! NH42SO4,NH4NO3,NH4CL,NA2SO4
         
       ELSE IF (drnh4no3 <= rh .AND. rh < drnh4cl) THEN
         scase = 'G2'
         CALL calcg2              ! NH42SO4,NH4CL,NA2SO4
         
       ELSE IF (drnh4cl <= rh  .AND. rh < drnh42s4) THEN
         scase = 'G3'
         CALL calcg3              ! NH42SO4,NA2SO4
         
       ELSE IF (drnh42s4 <= rh  .AND. rh < drna2so4) THEN
         scase = 'G4'
         CALL calcg4              ! NA2SO4
         
       ELSE IF (drna2so4 <= rh) THEN
         scase = 'G5'
         CALL calcg5              ! Only liquid
       END IF
     END IF
     
   ! *** SULFATE POOR ; SODIUM RICH
     
   ELSE IF (sulrat >= 2.0 .AND. sodrat >= 2.0) THEN
     
     IF(metstbl == 1) THEN
       scase = 'H6'
       CALL calch6                 ! Only liquid (metastable)
     ELSE
       
       IF (rh < drnh4no3) THEN
         scase = 'H1'
         CALL calch1              ! NH4NO3,NH4CL,NA2SO4,NACL,NANO3
         
       ELSE IF (drnh4no3 <= rh .AND. rh < drnano3) THEN
         scase = 'H2'
         CALL calch2              ! NH4CL,NA2SO4,NACL,NANO3
         
       ELSE IF (drnano3 <= rh  .AND. rh < drnacl) THEN
         scase = 'H3'
         CALL calch3              ! NH4CL,NA2SO4,NACL
         
       ELSE IF (drnacl <= rh   .AND. rh < drnh4cl) THEN
         scase = 'H4'
         CALL calch4              ! NH4CL,NA2SO4
         
       ELSE IF (drnh4cl <= rh .AND. rh < drna2so4) THEN
         scase = 'H5'
         CALL calch5              ! NA2SO4
         
       ELSE IF (drna2so4 <= rh) THEN
         scase = 'H6'
         CALL calch6              ! NO SOLID
       END IF
     END IF
     
   ! *** SULFATE RICH (NO ACID)
     
   ELSE IF (1.0 <= sulrat .AND. sulrat < 2.0) THEN
     
     IF(metstbl == 1) THEN
       scase = 'I6'
       CALL calci6                 ! Only liquid (metastable)
     ELSE
       
       IF (rh < drnh4hs4) THEN
         scase = 'I1'
         CALL calci1              ! NA2SO4,(NH4)2SO4,NAHSO4,NH4HSO4,LC
         
       ELSE IF (drnh4hs4 <= rh .AND. rh < drnahso4) THEN
         scase = 'I2'
         CALL calci2              ! NA2SO4,(NH4)2SO4,NAHSO4,LC
         
       ELSE IF (drnahso4 <= rh .AND. rh < drlc) THEN
         scase = 'I3'
         CALL calci3              ! NA2SO4,(NH4)2SO4,LC
         
       ELSE IF (drlc <= rh     .AND. rh < drnh42s4) THEN
         scase = 'I4'
         CALL calci4              ! NA2SO4,(NH4)2SO4
         
       ELSE IF (drnh42s4 <= rh .AND. rh < drna2so4) THEN
         scase = 'I5'
         CALL calci5              ! NA2SO4
         
       ELSE IF (drna2so4 <= rh) THEN
         scase = 'I6'
         CALL calci6              ! NO SOLIDS
       END IF
     END IF
     
     CALL calcnha                ! MINOR SPECIES: HNO3, HCl
     CALL calcnh3                !                NH3
     
   ! *** SULFATE RICH (FREE ACID)
     
   ELSE IF (sulrat < 1.0) THEN
     
     IF(metstbl == 1) THEN
       scase = 'J3'
       CALL calcj3                 ! Only liquid (metastable)
     ELSE
       
       IF (rh < drnh4hs4) THEN
         scase = 'J1'
         CALL calcj1              ! NH4HSO4,NAHSO4
         
       ELSE IF (drnh4hs4 <= rh .AND. rh < drnahso4) THEN
         scase = 'J2'
         CALL calcj2              ! NAHSO4
         
       ELSE IF (drnahso4 <= rh) THEN
         scase = 'J3'
         CALL calcj3
       END IF
     END IF
     
     CALL calcnha                ! MINOR SPECIES: HNO3, HCl
     CALL calcnh3                !                NH3
   END IF
   
END SUBROUTINE  isrp3f

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCA2
!! *** CASE A2
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0)
!!     2. LIQUID AEROSOL PHASE ONLY POSSIBLE
!!     FOR CALCULATIONS, A BISECTION IS PERFORMED TOWARDS X, THE
!!     AMOUNT OF HYDROGEN IONS (H+) FOUND IN THE LIQUID PHASE.
!!     FOR EACH ESTIMATION OF H+, FUNCTION FUNCB2A CALCULATES THE
!!     CONCENTRATION OF IONS FROM THE NH3(GAS) - NH4+(LIQ) EQUILIBRIUM.
!!     ELECTRONEUTRALITY IS USED AS THE OBJECTIVE FUNCTION.
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calca2
   USE Isorropia_Module
   implicit none
   
   double precision :: omelo,omehi,x1,y1,funca2,dx,x2,y2,x3,y3
   INTEGER :: i
   
   ! *** SETUP PARAMETERS ************************************************
   
   calaou    =.true.       ! Outer loop activity calculation flag
   omelo     = tiny        ! Low  limit: SOLUTION IS VERY BASIC
   omehi     = 2.0D0*w(2)  ! High limit: FROM NH4+ -> NH3(g) + H+(aq)
   
   ! *** CALCULATE WATER CONTENT *****************************************
   
   molal(5) = w(2)
   molal(6) = zero
   CALL calcmr
   
   ! *** INITIAL VALUES FOR BISECTION ************************************
   
   x1 = omehi
   y1 = funca2 (x1)
   IF (ABS(y1) <= eps) RETURN
   
   ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
   
   dx = (omehi-omelo)/FLOAT(ndiv)
   DO  i=1,ndiv
     x2 = MAX(x1-dx, omelo)
     y2 = funca2 (x2)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y2) < zero) GO TO 20  ! (Y1*Y2.LT.ZERO)
     x1 = x2
     y1 = y2
   END DO
   IF (ABS(y2) <= eps) THEN
     RETURN
   ELSE
     CALL pusherr (0001, 'CALCA2')    ! WARNING ERROR: NO SOLUTION
     RETURN
   END IF
   
   ! *** PERFORM BISECTION ***********************************************
   
   20    DO  i=1,maxit
     x3 = 0.5*(x1+x2)
     y3 = funca2 (x3)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y3) <= zero) THEN  ! (Y1*Y3 .LE. ZERO)
       y2    = y3
       x2    = x3
     ELSE
       y1    = y3
       x1    = x3
     END IF
     IF (ABS(x2-x1) <= eps*x1) GO TO 40
   END DO
   CALL pusherr (0002, 'CALCA2')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CONVERGED ; RETURN **********************************************
   
   40    x3 = 0.5*(x1+x2)
   y3 = funca2 (x3)

END SUBROUTINE  calca2
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** FUNCTION FUNCA2
!! *** CASE A2
!!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE A2 ;
!!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCA2.
!!=======================================================================
DOUBLE PRECISION FUNCTION funca2 (omegi)
   USE Isorropia_Module
   implicit none
   
   double precision, INTENT(IN)                         :: omegi
   DOUBLE PRECISION :: lamda,psi,a1,a2,a3,zeta,denom
   INTEGER :: i
   
   ! *** SETUP PARAMETERS ************************************************
   
   frst   = .true.
   calain = .true.
   psi    = w(2)         ! INITIAL AMOUNT OF (NH4)2SO4 IN SOLUTION
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     a1    = xk1*water/gama(7)*(gama(8)/gama(7))**2.
     a2    = xk2*r*temp/xkw*(gama(8)/gama(9))**2.
     a3    = xkw*rh*water*water
     
     lamda = psi/(a1/omegi+one)
     zeta  = a3/omegi
     
   ! *** SPECIATION & WATER CONTENT ***************************************
     
     molal (1) = omegi                        ! HI
     molal (3) = w(3)/(one/a2/omegi + one)    ! NH4I
     molal (5) = MAX(psi-lamda,tiny)          ! SO4I
     molal (6) = lamda                        ! HSO4I
     gnh3      = MAX (w(3)-molal(3), tiny)    ! NH3GI
     coh       = zeta                         ! OHI
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   
   ! *** CALCULATE OBJECTIVE FUNCTION ************************************
   
   20    CONTINUE
   denom = (2.0*molal(5)+molal(6)+organion)
   funca2= (molal(3)/denom - one) + molal(1)/denom

END FUNCTION funca2

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCA1
!! *** CASE A1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0)
!!     2. SOLID AEROSOL ONLY
!!     3. SOLIDS POSSIBLE : (NH4)2SO4
!!     A SIMPLE MATERIAL BALANCE IS PERFORMED, AND THE SOLID (NH4)2SO4
!!     IS CALCULATED FROM THE SULFATES. THE EXCESS AMMONIA REMAINS IN
!!     THE GAS PHASE.
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calca1
   USE Isorropia_Module
   implicit none
   
   cnh42s4 = w(2)
   gnh3    = MAX (w(3)-2.0*cnh42s4, zero)

END SUBROUTINE  calca1
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCB4
!! *** CASE B4
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!!     2. LIQUID AEROSOL PHASE ONLY POSSIBLE
!!     FOR CALCULATIONS, A BISECTION IS PERFORMED WITH RESPECT TO H+.
!!     THE OBJECTIVE FUNCTION IS THE DIFFERENCE BETWEEN THE ESTIMATED H+
!!     AND THAT CALCULATED FROM ELECTRONEUTRALITY.
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcb4
   USE Isorropia_Module
   implicit none
   
   double precision :: ak1,bet,gam,bb,cc,dd
   INteger :: i
   
   ! *** SOLVE EQUATIONS **************************************************
   
   frst       = .true.
   calain     = .true.
   calaou     = .true.
   
   ! *** CALCULATE WATER CONTENT ******************************************
   
   CALL calcb1a         ! GET DRY SALT CONTENT, AND USE FOR WATER.
   molalr(13) = clc
   molalr(9)  = cnh4hs4
   molalr(4)  = cnh42s4
   clc        = zero
   cnh4hs4    = zero
   cnh42s4    = zero
 !srf  water      = molalr(13)/(m0(13)+tinydenom)+molalr(9)/(m0(9)+tinydenom)+molalr(4)/(m0(4)+tinydenom)
     water      = molalr(13)/m0(13)+molalr(9)/m0(9)+molalr(4)/m0(4)

!if(m0(13)< tinydenom .or. m0(4)< tinydenom .or.m0(9)< tinydenom ) then
!  print*,"calcb4@water =",water;call flush(6)
!  stop 3333
!endif 
 !  if(m0(13)>0.0D+00 .and. m0(9)>0.0D+00 .and.m0(4)>0.0D+00) THEN
 !    water      = molalr(13)/m0(13)+molalr(9)/m0(9)+molalr(4)/m0(4)
 !  else
 !    water      =tiny
 !  endif
 
   water      = water +watorg
   
 
   molal(3)   = w(3)   ! NH4I
   
   DO  i=1,nsweep
     ak1   = xk1*((gama(8)/gama(7))**2.)*(water/gama(7))
     bet   = w(2)
     gam   = molal(3)
     
     bb    = bet + ak1 - gam + organion
     cc    =-ak1*bet
     dd    = MAX(bb*bb-4.d0*cc, 0.d0)       ! US patch 12/20/01
     
   ! *** SPECIATION & WATER CONTENT ***************************************
     
     molal (5) = MAX(tiny,MIN(0.5*(-bb + SQRT(dd)), w(2))) ! SO4I
     molal (6) = MAX(tiny,MIN(w(2)-molal(5),w(2)))         ! HSO4I
     molal (1) = MAX(tiny,MIN(ak1*molal(6)/molal(5),w(2))) ! HI
     CALL calcmr                                           ! Water content
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (.NOT.calain) GO TO 30
     CALL calcact
   END DO
   
   30    RETURN
   
END SUBROUTINE  calcb4

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCB3
!! *** CASE B3
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!!     2. BOTH LIQUID & SOLID PHASE IS POSSIBLE
!!     3. SOLIDS POSSIBLE: (NH4)2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcb3
   USE Isorropia_Module
   implicit none
   
   double precision :: x,y,tlc,tnh42s4,tnh4hs4
   
   ! *** CALCULATE EQUIVALENT AMOUNT OF HSO4 AND SO4 ***********************
   
   x = MAX(2*w(2)-w(3), zero)   ! Equivalent NH4HSO4
   y = MAX(w(3)  -w(2), zero)   ! Equivalent NH42SO4
   
   ! *** CALCULATE SPECIES ACCORDING TO RELATIVE ABUNDANCE OF HSO4 *********
   
   IF (x < y) THEN             ! LC is the MIN (x,y)
     scase   = 'B3 ; SUBCASE 1'
     tlc     = x
     tnh42s4 = y-x
     CALL calcb3a (tlc,tnh42s4)      ! LC + (NH4)2SO4
   ELSE
     scase   = 'B3 ; SUBCASE 2'
     tlc     = y
     tnh4hs4 = x-y
     CALL calcb3b (tlc,tnh4hs4)      ! LC + NH4HSO4
   END IF
   
END SUBROUTINE  calcb3
   
!!=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCB3A
!! *** CASE B3 ; SUBCASE 1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH (1.0 < SULRAT < 2.0)
!!     2. BOTH LIQUID & SOLID PHASE IS POSSIBLE
!!     3. SOLIDS POSSIBLE: (NH4)2SO4
!!     FOR CALCULATIONS, A BISECTION IS PERFORMED TOWARDS ZETA, THE
!!     AMOUNT OF SOLID (NH4)2SO4 DISSOLVED IN THE LIQUID PHASE.
!!     FOR EACH ESTIMATION OF ZETA, FUNCTION FUNCB3A CALCULATES THE
!!     AMOUNT OF H+ PRODUCED (BASED ON THE SO4 RELEASED INTO THE
!!     SOLUTION). THE SOLUBILITY PRODUCT OF (NH4)2SO4 IS USED AS THE
!!     OBJECTIVE FUNCTION.
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcb3a (tlc, tnh42s4)
   USE Isorropia_Module
   implicit none
   
   double precision :: zlo,zli,tnh42s4,z1,y1,ylo,dz,y2,z2,zk,tlc
   double precision :: zhi,yhi,z3,y3
   double precision, external :: funcb3a
   integer :: i
   
   calaou = .true.         ! Outer loop activity calculation flag
   zlo    = zero           ! MIN DISSOLVED (NH4)2SO4
   zhi    = tnh42s4        ! MAX DISSOLVED (NH4)2SO4
   
   ! *** INITIAL VALUES FOR BISECTION (DISSOLVED (NH4)2SO4 ****************
   
   z1 = zlo
   y1 = funcb3a (z1, tlc, tnh42s4)
   IF (ABS(y1) <= eps) RETURN
   ylo= y1
   
   ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO ***********************
   
   dz = (zhi-zlo)/FLOAT(ndiv)
   DO  i=1,ndiv
     z2 = z1+dz
     y2 = funcb3a (z2, tlc, tnh42s4)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y2) < zero) GO TO 20  ! (Y1*Y2.LT.ZERO)
     z1 = z2
     y1 = y2
   END DO
   
   ! *** NO SUBDIVISION WITH SOLUTION FOUND
   
   yhi= y1                      ! Save Y-value at HI position
   IF (ABS(y2) < eps) THEN   ! x2 IS A SOLUTION
     RETURN
     
   ! *** { YLO, YHI } < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH LC
     
   ELSE IF (ylo < zero .AND. yhi < zero) THEN
     z1 = zhi
     z2 = zhi
     GO TO 40
     
   ! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH LC
     
   ELSE IF (ylo > zero .AND. yhi > zero) THEN
     z1 = zlo
     z2 = zlo
     GO TO 40
   ELSE
     CALL pusherr (0001, 'CALCB3A')    ! WARNING ERROR: NO SOLUTION
     RETURN
   END IF
   
   ! *** PERFORM BISECTION ***********************************************
   
   20    DO  i=1,maxit
     z3 = 0.5*(z1+z2)
     y3 = funcb3a (z3, tlc, tnh42s4)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y3) <= zero) THEN  ! (Y1*Y3 .LE. ZERO)
       y2    = y3
       z2    = z3
     ELSE
       y1    = y3
       z1    = z3
     END IF
     IF (ABS(z2-z1) <= eps*z1) GO TO 40
   END DO
   CALL pusherr (0002, 'CALCB3A')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CONVERGED ; RETURN ************************************************
   
   40    zk = 0.5*(z1+z2)
   y3 = funcb3a (zk, tlc, tnh42s4)
   
END SUBROUTINE  calcb3a
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** FUNCTION FUNCB3A
!! *** CASE B3 ; SUBCASE 1
!!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE B3
!!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCA3.
!!=======================================================================
DOUBLE PRECISION FUNCTION funcb3a (zk, y, x)
   USE Isorropia_Module
   implicit none  
   
   double precision, INTENT(IN)                         :: zk
   double precision, INTENT(IN)                         :: y
   double precision, INTENT(IN)                         :: x
   DOUBLE PRECISION :: kk,dd,grat1
   integer :: i
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   frst   = .true.
   calain = .true.
   DO  i=1,nsweep
     grat1 = xk1*water/gama(7)*(gama(8)/gama(7))**2.
     dd    = SQRT( (zk+grat1+y+organion)**2. + 4.0*(y*grat1  &
         - organion* (y+zk))  )
     kk    = 0.5*(-(zk+grat1+y+organion) + dd )
     
   ! *** SPECIATION & WATER CONTENT ***************************************
     
     molal (1) = kk + organion         ! HI
     molal (5) = kk+zk+y           ! SO4I
     molal (6) = MAX (y-kk, tiny)  ! HSO4I
     molal (3) = 3.0*y+2*zk        ! NH4I
     cnh42s4   = x-zk              ! Solid (NH4)2SO4
     CALL calcmr                   ! Water content
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 30
     END IF
   END DO
   
   ! *** CALCULATE OBJECTIVE FUNCTION ************************************
   
   !CC30    FUNCB3A= ( SO4I*NH4I**2.0 )/( XK7*(WATER/GAMA(4))**3.0 )
   30    funcb3a= molal(5)*molal(3)**2.0
   funcb3a= funcb3a/(xk7*(water/gama(4))**3.0) - one

END FUNCTION funcb3a
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCB3B
!! *** CASE B3 ; SUBCASE 2
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH (1.0 < SULRAT < 2.0)
!!     2. LIQUID PHASE ONLY IS POSSIBLE
!!     SPECIATION CALCULATIONS IS BASED ON THE HSO4 <--> SO4 EQUILIBRIUM.
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcb3b (y, x)
   USE Isorropia_Module
   implicit none   
   
   DOUBLE PRECISION, INTENT(IN)                         :: y
   DOUBLE PRECISION, INTENT(IN)                         :: x
   DOUBLE PRECISION :: kk,dd,grat1
   integer :: i
   
   calaou = .false.        ! Outer loop activity calculation flag
   frst   = .false.
   calain = .true.
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     grat1 = xk1*water/gama(7)*(gama(8)/gama(7))**2.
     dd    = SQRT( (grat1+y+organion)**2. + 4.0* ( (x+y)*grat1 - y * organion))
     kk    = 0.5*(-(grat1+y+organion) + dd )
     
   ! *** SPECIATION & WATER CONTENT ***************************************
     
     molal (1) = kk + organion        ! HI
     molal (5) = y+kk                 ! SO4I
     molal (6) = MAX (x+y-kk, tiny)   ! HSO4I
     molal (3) = 3.0*y+x              ! NH4I
     CALL calcmr                      ! Water content
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (.NOT.calain) GO TO 30
     CALL calcact
   END DO
   
   30    RETURN
   
END SUBROUTINE  calcb3b

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCB2
!! *** CASE B2
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : LC, (NH4)2SO4
!!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON THE SULFATE RATIO:
!!     1. WHEN BOTH LC AND (NH4)2SO4 ARE POSSIBLE (SUBROUTINE CALCB2A)
!!     2. WHEN ONLY LC IS POSSIBLE (SUBROUTINE CALCB2B).
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcb2
   USE Isorropia_Module
   implicit none
   
   double precision :: x,y
   ! *** CALCULATE EQUIVALENT AMOUNT OF HSO4 AND SO4 ***********************
   
   x = MAX(2*w(2)-w(3), tiny)   ! Equivalent NH4HSO4
   y = MAX(w(3)  -w(2), tiny)   ! Equivalent NH42SO4
   
   ! *** CALCULATE SPECIES ACCORDING TO RELATIVE ABUNDANCE OF HSO4 *********
   
   IF (x <= y) THEN             ! LC is the MIN (x,y)
     scase = 'B2 ; SUBCASE 1'
     CALL calcb2a (x,y-x)      ! LC + (NH4)2SO4 POSSIBLE
   ELSE
     scase = 'B2 ; SUBCASE 2'
     CALL calcb2b (y,x-y)      ! LC ONLY POSSIBLE
   END IF
   
END SUBROUTINE  calcb2
   
!!=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCB2
!! *** CASE B2 ; SUBCASE A.
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH (1.0 < SULRAT < 2.0)
!!     2. SOLID PHASE ONLY POSSIBLE
!!     3. SOLIDS POSSIBLE: LC, (NH4)2SO4
!!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
!!     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!!     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE
!!     FOR SOLID CALCULATIONS, A MATERIAL BALANCE BASED ON THE STOICHIMETRIC
!!     PROPORTION OF AMMONIUM AND SULFATE IS DONE TO CALCULATE THE AMOUNT
!!     OF LC AND (NH4)2SO4 IN THE SOLID PHASE.
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcb2a (tlc, tnh42s4)
   USE Isorropia_Module
   implicit none
   
   double precision :: tlc,tnh42s4
   ! *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************
   
   IF (rh < drmlcas) THEN
     scase   = 'B2 ; SUBCASE A1'    ! SOLIDS POSSIBLE ONLY
     clc     = tlc
     cnh42s4 = tnh42s4
     scase   = 'B2 ; SUBCASE A1'
   ELSE
     scase = 'B2 ; SUBCASE A2'
     CALL calcb2a2 (tlc, tnh42s4)   ! LIQUID & SOLID PHASE POSSIBLE
     scase = 'B2 ; SUBCASE A2'
   END IF
   
END SUBROUTINE  calcb2a
   
!!======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCB2A2
!! *** CASE B2 ; SUBCASE A2.
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH (1.0 < SULRAT < 2.0)
!!     2. SOLID PHASE ONLY POSSIBLE
!!     3. SOLIDS POSSIBLE: LC, (NH4)2SO4
!!     THIS IS THE CASE WHERE THE RELATIVE HUMIDITY IS IN THE MUTUAL
!!     DRH REGION. THE SOLUTION IS ASSUMED TO BE THE SUM OF TWO WEIGHTED
!!     SOLUTIONS ; THE SOLID PHASE ONLY (SUBROUTINE CALCB2A1) AND THE
!!     THE SOLID WITH LIQUID PHASE (SUBROUTINE CALCB3).
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcb2a2 (tlc, tnh42s4)
   USE Isorropia_Module
   implicit none
   
   double precision :: wf,onemwf,clco,cnh42so,tnh42s4,tlc
   ! *** FIND WEIGHT FACTOR **********************************************
   
   IF (wftyp == 0) THEN
     wf = zero
   ELSE IF (wftyp == 1) THEN
     wf = 0.5D0
   ELSE
     wf = (drlc-rh)/(drlc-drmlcas)
   END IF
   onemwf  = one - wf
   
   ! *** FIND FIRST SECTION ; DRY ONE ************************************
   
   clco     = tlc                     ! FIRST (DRY) SOLUTION
   cnh42so  = tnh42s4
   
   ! *** FIND SECOND SECTION ; DRY & LIQUID ******************************
   
   clc     = zero
   cnh42s4 = zero
   CALL calcb3                        ! SECOND (LIQUID) SOLUTION
   
   ! *** FIND SOLUTION AT MDRH BY WEIGHTING DRY & LIQUID SOLUTIONS.
   
   molal(1)= onemwf*molal(1)                                   ! H+
   molal(3)= onemwf*(2.d0*(cnh42so-cnh42s4) + 3.d0*(clco-clc)) ! NH4+
   molal(5)= onemwf*(cnh42so-cnh42s4 + clco-clc)               ! SO4--
   molal(6)= onemwf*(clco-clc)                                 ! HSO4-
   
   water   = onemwf*water
   
   clc     = wf*clco    + onemwf*clc
   cnh42s4 = wf*cnh42so + onemwf*cnh42s4
   
END SUBROUTINE  calcb2a2
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCB2
!! *** CASE B2 ; SUBCASE B
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH (1.0 < SULRAT < 2.0)
!!     2. BOTH LIQUID & SOLID PHASE IS POSSIBLE
!!     3. SOLIDS POSSIBLE: LC
!!     FOR CALCULATIONS, A BISECTION IS PERFORMED TOWARDS ZETA, THE
!!     AMOUNT OF SOLID LC DISSOLVED IN THE LIQUID PHASE.
!!     FOR EACH ESTIMATION OF ZETA, FUNCTION FUNCB2A CALCULATES THE
!!     AMOUNT OF H+ PRODUCED (BASED ON THE HSO4, SO4 RELEASED INTO THE
!!     SOLUTION). THE SOLUBILITY PRODUCT OF LC IS USED AS THE OBJECTIVE
!!     FUNCTION.
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcb2b (tlc,tnh4hs4)
   USE Isorropia_Module
   implicit none
   
   double precision ::  zlo,zhi,tlc,x1,y1,yhi,dx,x2,y2,x3,y3
   double precision ::  tnh4hs4,ylo
   double precision, external :: funcb2b
   integer :: i
   
   calaou = .true.       ! Outer loop activity calculation flag
   zlo    = zero
   zhi    = tlc          ! High limit: all of it in liquid phase
   
   ! *** INITIAL VALUES FOR BISECTION **************************************
   
   x1 = zhi
   y1 = funcb2b (x1,tnh4hs4,tlc)
   IF (ABS(y1) <= eps) RETURN
   yhi= y1                        ! Save Y-value at Hi position
   
   ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO ************************
   
   dx = (zhi-zlo)/ndiv
   DO  i=1,ndiv
     x2 = x1-dx
     y2 = funcb2b (x2,tnh4hs4,tlc)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y2) < zero) GO TO 20  ! (Y1*Y2.LT.ZERO)
     x1 = x2
     y1 = y2
   END DO
   
   ! *** NO SUBDIVISION WITH SOLUTION FOUND
   
   ylo= y1                      ! Save Y-value at LO position
   IF (ABS(y2) < eps) THEN   ! x2 IS A SOLUTION
     RETURN
     
   ! *** { YLO, YHI } < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH LC
     
   ELSE IF (ylo < zero .AND. yhi < zero) THEN
     x1 = zhi
     x2 = zhi
     GO TO 40
     
   ! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH LC
     
   ELSE IF (ylo > zero .AND. yhi > zero) THEN
     x1 = zlo
     x2 = zlo
     GO TO 40
   ELSE
     CALL pusherr (0001, 'CALCB2B')    ! WARNING ERROR: NO SOLUTION
     RETURN
   END IF
   
   ! *** PERFORM BISECTION *************************************************
   
   20    DO  i=1,maxit
     x3 = 0.5*(x1+x2)
     y3 = funcb2b (x3,tnh4hs4,tlc)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y3) <= zero) THEN  ! (Y1*Y3 .LE. ZERO)
       y2    = y3
       x2    = x3
     ELSE
       y1    = y3
       x1    = x3
     END IF
     IF (ABS(x2-x1) <= eps*x1) GO TO 40
   END DO
   CALL pusherr (0002, 'CALCB2B')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CONVERGED ; RETURN ************************************************
   
   40    x3 = 0.5*(x1+x2)
   y3 = funcb2b (x3,tnh4hs4,tlc)
   
END SUBROUTINE  calcb2b
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** FUNCTION FUNCB2B
!! *** CASE B2 ;
!!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE B2 ; SUBCASE 2
!!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCB2B.
!!=======================================================================
DOUBLE PRECISION FUNCTION funcb2b (x,tnh4hs4,tlc)
   USE Isorropia_Module
   implicit none
   
   double precision :: grat2,parm,x,delta, tnh4hs4,tlc,omega
   integer :: i
   
   ! *** SOLVE EQUATIONS **************************************************
   
   frst   = .true.
   calain = .true.
   DO  i=1,nsweep
     grat2 = xk1*water*(gama(8)/gama(7))**2./gama(7)
     parm  = x+grat2+organion
     delta = parm*parm + 4.0*((x+tnh4hs4)*grat2-x*organion)
     omega = 0.5*(-parm + SQRT(delta))         ! Thetiki riza (ie:H+>0)
     
   ! *** SPECIATION & WATER CONTENT ***************************************
     
     molal (1) = omega + organion              ! HI
     molal (3) = 3.0*x+tnh4hs4                 ! NH4I
     molal (5) = x+omega                       ! SO4I
     molal (6) = MAX (x+tnh4hs4-omega, tiny)   ! HSO4I
     clc       = MAX(tlc-x,zero)               ! Solid LC
     CALL calcmr                               ! Water content
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP ******************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 30
     END IF
   END DO
   
   ! *** CALCULATE OBJECTIVE FUNCTION **************************************
   
   !CC30    FUNCB2B= ( NH4I**3.*SO4I*HSO4I )/( XK13*(WATER/GAMA(13))**5. )
   30    funcb2b= (molal(3)**3.)*molal(5)*molal(6)
   funcb2b= funcb2b/(xk13*(water/gama(13))**5.) - one

END FUNCTION funcb2b
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCB1
!! *** CASE B1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : LC, (NH4)2SO4, NH4HSO4
!!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
!!     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!!     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCB1A)
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcb1
   USE Isorropia_Module
   implicit none
   
   ! *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************
   
   IF (rh < drmlcab) THEN
     scase = 'B1 ; SUBCASE 1'
     CALL calcb1a              ! SOLID PHASE ONLY POSSIBLE
     scase = 'B1 ; SUBCASE 1'
   ELSE
     scase = 'B1 ; SUBCASE 2'
     CALL calcb1b              ! LIQUID & SOLID PHASE POSSIBLE
     scase = 'B1 ; SUBCASE 2'
   END IF
   
END SUBROUTINE  calcb1
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCB1A
!! *** CASE B1 ; SUBCASE 1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH
!!     2. THERE IS NO LIQUID PHASE
!!     3. SOLIDS POSSIBLE: LC, { (NH4)2SO4  XOR  NH4HSO4 } (ONE OF TWO
!!                         BUT NOT BOTH)
!!     A SIMPLE MATERIAL BALANCE IS PERFORMED, AND THE AMOUNT OF LC
!!     IS CALCULATED FROM THE (NH4)2SO4 AND NH4HSO4 WHICH IS LEAST
!!     ABUNDANT (STOICHIMETRICALLY). THE REMAINING EXCESS OF SALT
!!     IS MIXED WITH THE LC.
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcb1a
   USE Isorropia_Module
   implicit none
   
   DOUBLE precision :: x,y
   ! *** SETUP PARAMETERS ************************************************
   
   x = 2*w(2)-w(3)       ! Equivalent NH4HSO4
   y = w(3)-w(2)         ! Equivalent (NH4)2SO4
   
   ! *** CALCULATE COMPOSITION *******************************************
   
   IF (x <= y) THEN      ! LC is the MIN (x,y)
     clc     = x        ! NH4HSO4 >= (NH4)2S04
     cnh4hs4 = zero
     cnh42s4 = y-x
   ELSE
     clc     = y        ! NH4HSO4 <  (NH4)2S04
     cnh4hs4 = x-y
     cnh42s4 = zero
   END IF

END SUBROUTINE  calcb1a
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCB1B
!! *** CASE B1 ; SUBCASE 2
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE: LC, { (NH4)2SO4  XOR  NH4HSO4 } (ONE OF TWO
!!                         BUT NOT BOTH)
!!     THIS IS THE CASE WHERE THE RELATIVE HUMIDITY IS IN THE MUTUAL
!!     DRH REGION. THE SOLUTION IS ASSUMED TO BE THE SUM OF TWO WEIGHTED
!!     SOLUTIONS ; THE SOLID PHASE ONLY (SUBROUTINE CALCB1A) AND THE
!!     THE SOLID WITH LIQUID PHASE (SUBROUTINE CALCB2).
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcb1b
   USE Isorropia_Module
   implicit none
   
   double precision :: wf,onemwf,clco,cnh42so,cnh4hso
   ! *** FIND WEIGHT FACTOR **********************************************
   
   IF (wftyp == 0) THEN
     wf = zero
   ELSE IF (wftyp == 1) THEN
     wf = 0.5D0
   ELSE
     wf = (drnh4hs4-rh)/(drnh4hs4-drmlcab)
   END IF
   onemwf  = one - wf
   
   ! *** FIND FIRST SECTION ; DRY ONE ************************************
   
   CALL calcb1a
   clco     = clc               ! FIRST (DRY) SOLUTION
   cnh42so  = cnh42s4
   cnh4hso  = cnh4hs4
   
   ! *** FIND SECOND SECTION ; DRY & LIQUID ******************************
   
   clc     = zero
   cnh42s4 = zero
   cnh4hs4 = zero
   CALL calcb2                  ! SECOND (LIQUID) SOLUTION
   
   ! *** FIND SOLUTION AT MDRH BY WEIGHTING DRY & LIQUID SOLUTIONS.
   
   molal(1)= onemwf*molal(1)                                   ! H+
   molal(3)= onemwf*(2.d0*(cnh42so-cnh42s4) + (cnh4hso-cnh4hs4)  &
       + 3.d0*(clco-clc))                          ! NH4+
   molal(5)= onemwf*(cnh42so-cnh42s4 + clco-clc)               ! SO4--
   molal(6)= onemwf*(cnh4hso-cnh4hs4 + clco-clc)               ! HSO4-
   
   water   = onemwf*water
   
   clc     = wf*clco    + onemwf*clc
   cnh42s4 = wf*cnh42so + onemwf*cnh42s4
   cnh4hs4 = wf*cnh4hso + onemwf*cnh4hs4
   
END SUBROUTINE  calcb1b
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCC2
!! *** CASE C2
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
!!     2. THERE IS ONLY A LIQUID PHASE
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcc2
   USE Isorropia_Module
   implicit none
   
   DOUBLE PRECISION :: lamda, kapa,psi,parm,bb,cc
   integer :: i
   
   calaou =.true.         ! Outer loop activity calculation flag
   frst   =.true.
   calain =.true.
   
   ! *** SOLVE EQUATIONS **************************************************
   
   lamda  = w(3)           ! NH4HSO4 INITIALLY IN SOLUTION
   psi    = w(2)-w(3)      ! H2SO4 IN SOLUTION
   DO  i=1,nsweep
     parm  = water*xk1/gama(7)*(gama(8)/gama(7))**2.
     bb    = psi+parm+organion
     cc    =-parm*(lamda+psi)
     kapa  = 0.5*(-bb+SQRT(bb*bb-4.0*cc))
     
   ! *** SPECIATION & WATER CONTENT ***************************************
     
     molal(1) = psi+kapa+organion                      ! HI
     molal(3) = lamda                                  ! NH4I
     molal(5) = kapa                                   ! SO4I
     molal(6) = MAX(lamda+psi-kapa, tiny)              ! HSO4I
     ch2so4   = MAX(molal(5)+molal(6)-molal(3), zero)  ! Free H2SO4
     CALL calcmr                                       ! Water content
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (.NOT.calain) GO TO 30
     CALL calcact
   END DO
   
   30    RETURN
   
END SUBROUTINE  calcc2
   
!>======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCC1
!! *** CASE C1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE: NH4HSO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcc1
   USE Isorropia_Module
   implicit none
   
   DOUBLE PRECISION :: klo, khi,x1,y1,ylo,dx,x2,y2,yhi,x3,y3
   integer :: i
   double precision, external :: funcc1
   
   calaou = .true.    ! Outer loop activity calculation flag
   klo    = tiny
   khi    = w(3)
   
   ! *** INITIAL VALUES FOR BISECTION *************************************
   
   x1 = klo
   y1 = funcc1 (x1)
   IF (ABS(y1) <= eps) GO TO 50
   ylo= y1
   
   ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO ***********************
   
   dx = (khi-klo)/FLOAT(ndiv)
   DO  i=1,ndiv
     x2 = x1+dx
     y2 = funcc1 (x2)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y2) < zero) GO TO 20 ! (Y1*Y2 .LT. ZERO)
     x1 = x2
     y1 = y2
   END DO
   
   ! *** NO SUBDIVISION WITH SOLUTION FOUND
   
   yhi= y2                 ! Save Y-value at HI position
   IF (ABS(y2) < eps) THEN   ! x2 IS A SOLUTION
     GO TO 50
     
   ! *** { YLO, YHI } < 0.0  SOLUTION IS ALWAYS UNDERSATURATED WITH NH4HS04
     
   ELSE IF (ylo < zero .AND. yhi < zero) THEN
     GO TO 50
     
   ! *** { YLO, YHI } > 0.0 SOLUTION IS ALWAYS SUPERSATURATED WITH NH4HS04
     
   ELSE IF (ylo > zero .AND. yhi > zero) THEN
     x1 = klo
     x2 = klo
     GO TO 40
   ELSE
     CALL pusherr (0001, 'CALCC1')    ! WARNING ERROR: NO SOLUTION
     GO TO 50
   END IF
   
   ! *** PERFORM BISECTION OF DISSOLVED NH4HSO4 **************************
   
   20    DO  i=1,maxit
     x3 = 0.5*(x1+x2)
     y3 = funcc1 (x3)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y3) <= zero) THEN  ! (Y1*Y3 .LE. ZERO)
       y2    = y3
       x2    = x3
     ELSE
       y1    = y3
       x1    = x3
     END IF
     IF (ABS(x2-x1) <= eps*x1) GO TO 40
   END DO
   CALL pusherr (0002, 'CALCC1')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CONVERGED ; RETURN ***********************************************
   
   40    x3 = 0.5*(x1+x2)
   y3 = funcc1 (x3)
   
   50    RETURN
   
END SUBROUTINE  calcc1
   
!!=======================================================================
!! *** ISORROPIA CODE
!! *** FUNCTION FUNCC1
!! *** CASE C1 ;
!!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE C1
!!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCC1.
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CANREGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
DOUBLE PRECISION FUNCTION funcc1 (kapa)
   USE Isorropia_Module
   implicit none
   
   DOUBLE PRECISION, INTENT(IN)             :: kapa
   DOUBLE PRECISION :: lamda,psi,par1,par2,bb,cc
   integer :: i
   
   ! *** SOLVE EQUATIONS **************************************************
   
   frst   = .true.
   calain = .true.
   
   psi = w(2)-w(3)
   DO  i=1,nsweep
     par1  = xk1*water/gama(7)*(gama(8)/gama(7))**2.0
     par2  = xk12*(water/gama(9))**2.0
     bb    = psi + par1 + organion
     cc    =-par1*(psi+kapa)
     lamda = 0.5*(-bb+SQRT(bb*bb-4*cc))
     
   ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY *******************************
     
     molal(1) = psi+lamda+organion           ! HI
     molal(3) = kapa                         ! NH4I
     molal(5) = lamda                        ! SO4I
     molal(6) = MAX (zero, psi+kapa-lamda)   ! HSO4I
     cnh4hs4  = MAX(w(3)-molal(3), zero)     ! Solid NH4HSO4
     ch2so4   = MAX(psi, zero)               ! Free H2SO4
     CALL calcmr                             ! Water content
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 30
     END IF
   END DO
   
   ! *** CALCULATE ZERO FUNCTION *******************************************
   
   !CC30    FUNCC1= (NH4I*HSO4I/PAR2) - ONE
   30    funcc1= (molal(3)*molal(6)/par2) - one

END FUNCTION funcc1
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCD3
!! *** CASE D3
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0)
!!     2. THERE IS OLNY A LIQUID PHASE
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcd3
   USE Isorropia_Module
   use solut
   implicit none
   
   double precision :: delta,psi4lo,psi4hi,x1,y1,ylo,dx,x2,y2,yhi,p4,yy
   double precision :: x3,y3
   double precision, external :: funcd3
   integer :: i
   
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** FIND DRY COMPOSITION **********************************************
 !print*,"calcd3@water0", water;call flush(6)
   CALL calcd1a
   
 ! print*,"calcd3@water1", water;call flush(6)
  ! *** SETUP PARAMETERS ************************************************
   
   chi1 = cnh4no3               ! Save from CALCD1 run
   chi2 = cnh42s4
   chi3 = ghno3
   chi4 = gnh3
   
   psi1 = cnh4no3               ! ASSIGN INITIAL PSI's
   psi2 = chi2
   psi3 = zero
   psi4 = zero
   
   molal(5) = zero
   molal(6) = zero
   molal(3) = psi1
   molal(7) = psi1
   CALL calcmr                  ! Initial water
   !print*,"calcd3@water2", water;call flush(6)
  
   calaou = .true.              ! Outer loop activity calculation flag
   psi4lo = tiny                ! Low  limit
   psi4hi = chi4                ! High limit
   
   ! *** INITIAL VALUES FOR BISECTION ************************************
   
   60    x1 = psi4lo
   y1 = funcd3 (x1)
   IF (ABS(y1) <= eps) RETURN
   ylo= y1                 ! Save Y-value at HI position
   
   ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
   
   dx = (psi4hi-psi4lo)/FLOAT(ndiv)
   DO  i=1,ndiv
     x2 = x1+dx
     y2 = funcd3 (x2)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y2) < zero) GO TO 20  ! (Y1*Y2.LT.ZERO)
     x1 = x2
     y1 = y2
   END DO
   
   ! *** NO SUBDIVISION WITH SOLUTION FOUND
   
   yhi= y1                      ! Save Y-value at Hi position
   IF (ABS(y2) < eps) THEN   ! x2 IS A SOLUTION
     RETURN
     
   ! *** { YLO, YHI } < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH3
   ! Physically I dont know when this might happen, but I have put this
   ! branch in for completeness.
     
   ELSE IF (ylo < zero .AND. yhi < zero) THEN
     p4 = chi4
     yy = funcd3(p4)
     GO TO 50
     
   ! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH3
   ! This happens when Sul.Rat. = 2.0, so some NH4+ from sulfate evaporates
   ! and goes to the gas phase ; so I redefine the LO and HI limits of PSI4
   ! and proceed again with root tracking.
     
   ELSE IF (ylo > zero .AND. yhi > zero) THEN
     psi4hi = psi4lo
     psi4lo = psi4lo - 0.1*(psi1+psi2) ! No solution; some NH3 evaporates
     IF (psi4lo < -(psi1+psi2)) THEN
       CALL pusherr (0001, 'CALCD3')  ! WARNING ERROR: NO SOLUTION
       RETURN
     ELSE
       GO TO 60                        ! Redo root tracking
     END IF
   END IF
   
   ! *** PERFORM BISECTION ***********************************************
   
   20    DO  i=1,maxit
     x3 = 0.5*(x1+x2)
     y3 = funcd3 (x3)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y3) <= zero) THEN  ! (Y1*Y3 .LE. ZERO)
       y2    = y3
       x2    = x3
     ELSE
       y1    = y3
       x1    = x3
     END IF
     IF (ABS(x2-x1) <= eps*ABS(x1)) GO TO 40
   END DO
   CALL pusherr (0002, 'CALCD3')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CONVERGED ; RETURN **********************************************
   
   40    x3 = 0.5*(x1+x2)
   y3 = funcd3 (x3)
   
   ! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
   
   50    CONTINUE
   IF (molal(1) > tiny) THEN
     CALL calchs4 (molal(1), molal(5), zero, delta)
     molal(1) = molal(1) - delta                     ! H+   EFFECT
     molal(5) = molal(5) - delta                     ! SO4  EFFECT
     molal(6) = delta                                ! HSO4 EFFECT
   END IF
  ! print*,"calcd3@water4", water;call flush(6)

END SUBROUTINE  calcd3
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** FUNCTION FUNCD3
!! *** CASE D3
!!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE D3 ;
!!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCD3.
!!=======================================================================
DOUBLE PRECISION FUNCTION funcd3 (p4)
   USE Isorropia_Module
   use solut   
   implicit none
   
   double precision, INTENT(IN)                         :: p4
   double precision :: bb,denm,abb,ahi
   integer :: i
   
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   frst   = .true.
   calain = .true.
   psi4   = p4
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     a2   = xk7*(water/gama(4))**3.0
     a3   = xk4*r*temp*(water/gama(10))**2.0
     a4   = (xk2/xkw)*r*temp*(gama(10)/gama(5))**2.0
     a7   = xkw *rh*water*water
     
     psi3 = a3*a4*chi3*(chi4-psi4) - psi1*(2.d0*psi2+psi1+psi4)
     psi3 = psi3/(a3*a4*(chi4-psi4) + 2.d0*psi2+psi1+psi4)
     psi3 = MIN(MAX(psi3, zero), chi3)
     
     bb   = psi4 - psi3 - organion
   !CCOLD         AHI  = 0.5*(-BB + SQRT(BB*BB + 4.d0*A7)) ! This is correct also
   !CC         AHI  =2.0*A7/(BB+SQRT(BB*BB + 4.d0*A7)) ! Avoid overflow when HI->0
     denm = bb+SQRT(bb*bb + 4.d0*a7)
     IF (denm <= tiny) THEN       ! Avoid overflow when HI->0
       abb  = ABS(bb)
       denm = (bb+abb) + 2.0*a7/abb ! Taylor expansion of SQRT
     END IF
     ahi = 2.0*a7/denm
     
   ! *** SPECIATION & WATER CONTENT ***************************************
     
     molal (1) = ahi                             ! HI
     molal (3) = psi1 + psi4 + 2.d0*psi2         ! NH4I
     molal (5) = psi2                            ! SO4I
     molal (6) = zero                            ! HSO4I
     molal (7) = psi3 + psi1                     ! NO3I
     cnh42s4   = chi2 - psi2                     ! Solid (NH4)2SO4
     cnh4no3   = zero                            ! Solid NH4NO3
     ghno3     = chi3 - psi3                     ! Gas HNO3
     gnh3      = chi4 - psi4                     ! Gas NH3
     CALL calcmr                                 ! Water content
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   
   ! *** CALCULATE OBJECTIVE FUNCTION ************************************
   
   20    CONTINUE
   !CC      FUNCD3= NH4I/HI/MAX(GNH3,TINY)/A4 - ONE
   funcd3= molal(3)/molal(1)/MAX(gnh3,tiny)/a4 - one

END FUNCTION funcd3

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCD2
!! *** CASE D2
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcd2
   USE Isorropia_Module
   use solut
   implicit none
   
   double precision :: delta,psi4lo,psi4hi,x1,y1,ylo,dx,x2,y2,yhi,p4,yy
   double precision :: x3,y3
   double precision, external :: funcd2
   integer :: i
   
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** FIND DRY COMPOSITION **********************************************
   
   CALL calcd1a
   
   ! *** SETUP PARAMETERS ************************************************
   
   chi1 = cnh4no3               ! Save from CALCD1 run
   chi2 = cnh42s4
   chi3 = ghno3
   chi4 = gnh3
   
   psi1 = cnh4no3               ! ASSIGN INITIAL PSI's
   psi2 = cnh42s4
   psi3 = zero
   psi4 = zero
   
   molal(5) = zero
   molal(6) = zero
   molal(3) = psi1
   molal(7) = psi1
   CALL calcmr                  ! Initial water
   
   calaou = .true.              ! Outer loop activity calculation flag
   psi4lo = tiny                ! Low  limit
   psi4hi = chi4                ! High limit
   
   ! *** INITIAL VALUES FOR BISECTION ************************************
   
   60    x1 = psi4lo
   y1 = funcd2 (x1)
   IF (ABS(y1) <= eps) RETURN
   ylo= y1                 ! Save Y-value at HI position
   
   ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
   
   dx = (psi4hi-psi4lo)/FLOAT(ndiv)
   DO  i=1,ndiv
     x2 = x1+dx
     y2 = funcd2 (x2)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y2) < zero) GO TO 20  ! (Y1*Y2.LT.ZERO)
     x1 = x2
     y1 = y2
   END DO
   
   ! *** NO SUBDIVISION WITH SOLUTION FOUND
   
   yhi= y1                      ! Save Y-value at Hi position
   IF (ABS(y2) < eps) THEN   ! x2 IS A SOLUTION
     RETURN
     
   ! *** { YLO, YHI } < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH3
   ! Physically I dont know when this might happen, but I have put this
   ! branch in for completeness.
     
   ELSE IF (ylo < zero .AND. yhi < zero) THEN
     p4 = chi4
     yy = funcd2(p4)
     GO TO 50
     
   ! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH3
   ! This happens when Sul.Rat. = 2.0, so some NH4+ from sulfate evaporates
   ! and goes to the gas phase ; so I redefine the LO and HI limits of PSI4
   ! and proceed again with root tracking.
     
   ELSE IF (ylo > zero .AND. yhi > zero) THEN
     psi4hi = psi4lo
     psi4lo = psi4lo - 0.1*(psi1+psi2) ! No solution; some NH3 evaporates
     IF (psi4lo < -(psi1+psi2)) THEN
       CALL pusherr (0001, 'CALCD2')  ! WARNING ERROR: NO SOLUTION
       RETURN
     ELSE
       GO TO 60                        ! Redo root tracking
     END IF
   END IF
   
   ! *** PERFORM BISECTION ***********************************************
   
   20    DO  i=1,maxit
     x3 = 0.5*(x1+x2)
     y3 = funcd2 (x3)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y3) <= zero) THEN  ! (Y1*Y3 .LE. ZERO)
       y2    = y3
       x2    = x3
     ELSE
       y1    = y3
       x1    = x3
     END IF
     IF (ABS(x2-x1) <= eps*ABS(x1)) GO TO 40
   END DO
   CALL pusherr (0002, 'CALCD2')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CONVERGED ; RETURN **********************************************
   
   40    x3 = 0.5*(x1+x2)
   y3 = funcd2 (x3)
   
   ! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
   
   50    CONTINUE
   IF (molal(1) > tiny) THEN
     CALL calchs4 (molal(1), molal(5), zero, delta)
     molal(1) = molal(1) - delta                     ! H+   EFFECT
     molal(5) = molal(5) - delta                     ! SO4  EFFECT
     molal(6) = delta                                ! HSO4 EFFECT
   END IF

END SUBROUTINE  calcd2
   
!!=======================================================================
!! *** ISORROPIA CODE
!! *** FUNCTION FUNCD2
!! *** CASE D2
!!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE D2 ;
!!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCD2.
!!=======================================================================
DOUBLE PRECISION FUNCTION funcd2 (p4)
   USE Isorropia_Module
   use solut
   implicit none
   
   double precision, INTENT(IN)                         :: p4
   double precision :: psi14,bb,denm,abb,ahi
   integer :: i,islv
   
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   frst   = .true.
   calain = .true.
   psi4   = p4
   psi2   = chi2
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     a2  = xk7*(water/gama(4))**3.0
     a3  = xk4*r*temp*(water/gama(10))**2.0
     a4  = (xk2/xkw)*r*temp*(gama(10)/gama(5))**2.0
     a7  = xkw *rh*water*water
     
     IF (chi2 > tiny .AND. water > tiny) THEN
       psi14 = psi1+psi4
       CALL poly3 (psi14,0.25*psi14**2.,-a2/4.d0, psi2, islv)  ! PSI2
       IF (islv == 0) THEN
         psi2 = MIN (psi2, chi2)
       ELSE
         psi2 = zero
       END IF
     END IF
     
     psi3  = a3*a4*chi3*(chi4-psi4) - psi1*(2.d0*psi2+psi1+psi4)
     psi3  = psi3/(a3*a4*(chi4-psi4) + 2.d0*psi2+psi1+psi4)
     psi3  = MIN(MAX(psi3, zero), chi3)
     
     bb   = psi4 - psi3 - organion
   !CCOLD         AHI  = 0.5*(-BB + SQRT(BB*BB + 4.d0*A7)) ! This is correct also
     denm = bb+SQRT(bb*bb + 4.d0*a7)
     IF (denm <= tiny) THEN       ! Avoid overflow when HI->0
       abb  = ABS(bb)
       denm = (bb+abb) + 2.0*a7/abb ! Taylor expansion of SQRT
     END IF
     ahi = 2.0*a7/denm
     
   ! *** SPECIATION & WATER CONTENT ***************************************
     
     molal (1) = ahi                              ! HI
     molal (3) = psi1 + psi4 + 2.d0*psi2          ! NH4
     molal (5) = psi2                             ! SO4
     molal (6) = zero                             ! HSO4
     molal (7) = psi3 + psi1                      ! NO3
     cnh42s4   = chi2 - psi2                      ! Solid (NH4)2SO4
     cnh4no3   = zero                             ! Solid NH4NO3
     ghno3     = chi3 - psi3                      ! Gas HNO3
     gnh3      = chi4 - psi4                      ! Gas NH3
     CALL calcmr                                  ! Water content
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   
   ! *** CALCULATE OBJECTIVE FUNCTION ************************************
   
   20    CONTINUE
   !CC      FUNCD2= NH4I/HI/MAX(GNH3,TINY)/A4 - ONE
   funcd2= molal(3)/molal(1)/MAX(gnh3,tiny)/a4 - one

END FUNCTION funcd2

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCD1
!! *** CASE D1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0)
!!     2. SOLID AEROSOL ONLY
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3
!!     THERE ARE TWO REGIMES DEFINED BY RELATIVE HUMIDITY:
!!     1. RH < MDRH ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCD1A)
!!     2. RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcd1
   USE Isorropia_Module
   implicit none
   
   EXTERNAL calcd1a, calcd2
   
   ! *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************
   
   IF (rh < drmasan) THEN
     scase = 'D1 ; SUBCASE 1'   ! SOLID PHASE ONLY POSSIBLE
     CALL calcd1a
     scase = 'D1 ; SUBCASE 1'
   ELSE
     scase = 'D1 ; SUBCASE 2'   ! LIQUID & SOLID PHASE POSSIBLE
     CALL calcmdrh (rh, drmasan, drnh4no3, calcd1a, calcd2)
     scase = 'D1 ; SUBCASE 2'
   END IF
   
END SUBROUTINE  calcd1
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCD1A
!! *** CASE D1 ; SUBCASE 1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0)
!!     2. SOLID AEROSOL ONLY
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3
!!     THE SOLID (NH4)2SO4 IS CALCULATED FROM THE SULFATES, WHILE NH4NO3
!!     IS CALCULATED FROM NH3-HNO3 EQUILIBRIUM. 'ZE' IS THE AMOUNT OF
!!     NH4NO3 THAT VOLATIZES WHEN ALL POSSILBE NH4NO3 IS INITIALLY IN
!!     THE SOLID PHASE.
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcd1a
   USE Isorropia_Module
   implicit none
   
   double precision :: parm,x,ps,om,omps,diak,ze
   
   ! *** SETUP PARAMETERS ************************************************
   parm    = xk10/(r*temp)/(r*temp)
   
   ! *** CALCULATE NH4NO3 THAT VOLATIZES *********************************
   
   cnh42s4 = w(2)
   x       = MAX(zero, MIN(w(3)-2.0*cnh42s4, w(4)))  ! MAX NH4NO3
   ps      = MAX(w(3) - x - 2.0*cnh42s4, zero)
   om      = MAX(w(4) - x, zero)
   
   omps    = om+ps
   diak    = SQRT(omps*omps + 4.0*parm)              ! DIAKRINOUSA
   ze      = MIN(x, 0.5*(-omps + diak))              ! THETIKI RIZA
   
   ! *** SPECIATION *******************************************************
   
   cnh4no3 = x  - ze    ! Solid NH4NO3
   gnh3    = ps + ze    ! Gas NH3
   ghno3   = om + ze    ! Gas HNO3
   
END SUBROUTINE  calcd1a

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCG5
!! *** CASE G5
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcg5
   USE Isorropia_Module
   USE caseg
   implicit none
   
   DOUBLE PRECISION :: delta,psi6lo,psi6hi,x1,y1,dx,x2,y2,x3,y3
   integer :: i
   double precision, external :: funcg5a
   
   !LFR  COMMON /caseg/ chi1, chi2, chi3, chi4, chi5, chi6, lamda,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7
   
   ! *** SETUP PARAMETERS ************************************************
   
   calaou = .true.
   chi1   = 0.5*w(1)
   chi2   = MAX (w(2)-chi1, zero)
   chi3   = zero
   chi4   = MAX (w(3)-2.d0*chi2, zero)
   chi5   = w(4)
   chi6   = w(5)
   
   psi1   = chi1
   psi2   = chi2
   psi6lo = tiny
   psi6hi = chi6-tiny    ! MIN(CHI6-TINY, CHI4)
   
   water  = chi2/m0(4) + chi1/m0(2)
!srf   water  = chi2/(m0(4)+tinydenom) + chi1/(m0(2)+tinydenom)
   
   ! *** INITIAL VALUES FOR BISECTION ************************************
   
   x1 = psi6lo
   y1 = funcg5a (x1)
   IF (chi6 <= tiny) GO TO 50
   !cc      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50
   !cc      IF (WATER .LE. TINY) RETURN                    ! No water
   
   ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
   
   dx = (psi6hi-psi6lo)/FLOAT(ndiv)
   DO  i=1,ndiv
     x2 = x1+dx
     y2 = funcg5a (x2)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y2) < zero) GO TO 20  ! (Y1*Y2.LT.ZERO)
     x1 = x2
     y1 = y2
   END DO
   
   ! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
   
   IF (ABS(y2) > eps) y2 = funcg5a (psi6hi)
   GO TO 50
   
   ! *** PERFORM BISECTION ***********************************************
   
   20    DO  i=1,maxit
     x3 = 0.5*(x1+x2)
     y3 = funcg5a (x3)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y3) <= zero) THEN  ! (Y1*Y3 .LE. ZERO)
       y2    = y3
       x2    = x3
     ELSE
       y1    = y3
       x1    = x3
     END IF
     IF (ABS(x2-x1) <= eps*x1) GO TO 40
   END DO
   CALL pusherr (0002, 'CALCG5')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CONVERGED ; RETURN **********************************************
   
   40    x3 = 0.5*(x1+x2)
   y3 = funcg5a (x3)
   
   ! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
   
   50    CONTINUE
   IF (molal(1) > tiny .AND. molal(5) > tiny) THEN  ! If quadrat.called
     CALL calchs4 (molal(1), molal(5), zero, delta)
     molal(1) = molal(1) - delta                    ! H+   EFFECT
     molal(5) = molal(5) - delta                    ! SO4  EFFECT
     molal(6) = delta                               ! HSO4 EFFECT
   END IF
   
END SUBROUTINE  calcg5
   
!!=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE FUNCG5A
!! *** CASE G5
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
DOUBLE PRECISION FUNCTION funcg5a (x)
   USE Isorropia_Module
   USE caseg
   implicit none
   
   double precision, INTENT(IN)                         :: x
   double precision :: smin,hi,ohi,akk,bb,cc,dd
   integer :: i
   
   !DOUBLE PRECISION :: lamda
   !LFR  COMMON /caseg/ chi1, chi2, chi3, chi4, chi5, chi6, lamda,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7
   
   ! *** SETUP PARAMETERS ************************************************
   
   psi6   = x
   frst   = .true.
   calain = .true.
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     
     a1  = xk5 *(water/gama(2))**3.0
     a2  = xk7 *(water/gama(4))**3.0
     a4  = (xk2/xkw)*r*temp*(gama(10)/gama(5))**2.0
     a5  = xk4 *r*temp*(water/gama(10))**2.0
     a6  = xk3 *r*temp*(water/gama(11))**2.0
     akk = a4*a6
     
   !  CALCULATE DISSOCIATION QUANTITIES
     
     IF (chi5 >= tiny) THEN
       psi5 = psi6*chi5/(a6/a5*(chi6-psi6) + psi6)
     ELSE
       psi5 = tiny
     END IF
     
   !CC      IF(CHI4.GT.TINY) THEN
     IF(w(2) > tiny) THEN       ! Accounts for NH3 evaporation
       bb   =-(chi4 + psi6 + psi5 + 1.d0/a4)
       cc   = chi4*(psi5+psi6) - 2.d0*psi2/a4
       dd   = MAX(bb*bb-4.d0*cc, 0.d0)       ! US patch 11/26/01
       psi4 =0.5D0*(-bb - SQRT(dd))
     ELSE
       psi4 = tiny
     END IF
     
   ! *** CALCULATE SPECIATION ********************************************
     
     molal (2) = 2.0D0*psi1                          ! NAI
     molal (3) = 2.0*psi2 + psi4                     ! NH4I
     molal (4) = psi6                                ! CLI
     molal (5) = psi2 + psi1                         ! SO4I
     molal (6) = zero
     molal (7) = psi5                                ! NO3I
     
     smin  = 2.d0*molal(5)+molal(7)+molal(4)-molal(2)-molal(3) + organion
     CALL calcph (smin, hi, ohi)
     molal (1) = hi
     
     gnh3      = MAX(chi4 - psi4, tiny)              ! Gas NH3
     ghno3     = MAX(chi5 - psi5, tiny)              ! Gas HNO3
     ghcl      = MAX(chi6 - psi6, tiny)              ! Gas HCl
     
     cnh42s4   = zero                                ! Solid (NH4)2SO4
     cnh4no3   = zero                                ! Solid NH4NO3
     cnh4cl    = zero                                ! Solid NH4Cl
     
     CALL calcmr                                     ! Water content
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   
   ! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
   
   20    funcg5a = molal(1)*molal(4)/ghcl/a6 - one
   !CC         FUNCG5A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
   
END FUNCTION funcg5a
   
!!=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCG4
!! *** CASE G4
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcg4
   USE Isorropia_Module
   USE caseg
   implicit none

   DOUBLE PRECISION :: delta,psi6lo,psi6hi,x1,y1,dx,x2,y2,x3,y3
   double precision, external :: funcg4a
   integer :: i
   
   !LFR  COMMON /caseg/ chi1, chi2, chi3, chi4, chi5, chi6, lamda,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7
   
   ! *** SETUP PARAMETERS ************************************************
   
   calaou = .true.
   chi1   = 0.5*w(1)
   chi2   = MAX (w(2)-chi1, zero)
   chi3   = zero
   chi4   = MAX (w(3)-2.d0*chi2, zero)
   chi5   = w(4)
   chi6   = w(5)
   
   psi2   = chi2
   psi6lo = tiny
   psi6hi = chi6-tiny    ! MIN(CHI6-TINY, CHI4)
   
   water  = chi2/m0(4) + chi1/m0(2)
 !srf  water  = chi2/(m0(4)+tinydenom) + chi1/(m0(2)+tinydenom)
   
   ! *** INITIAL VALUES FOR BISECTION ************************************
   
   x1 = psi6lo
   y1 = funcg4a (x1)
   IF (chi6 <= tiny) GO TO 50
   !CC      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY .OR. WATER .LE. TINY) GOTO 50
   !CC      IF (WATER .LE. TINY) RETURN                    ! No water
   
   ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
   
   dx = (psi6hi-psi6lo)/FLOAT(ndiv)
   DO  i=1,ndiv
     x2  = x1+dx
     y2  = funcg4a (x2)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y2) < zero) GO TO 20  ! (Y1*Y2.LT.ZERO)
     x1  = x2
     y1  = y2
   END DO
   
   ! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
   
   IF (ABS(y2) > eps) y2 = funcg4a (psi6lo)
   GO TO 50
   
   ! *** PERFORM BISECTION ***********************************************
   
   20    DO  i=1,maxit
     x3 = 0.5*(x1+x2)
     y3 = funcg4a (x3)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y3) <= zero) THEN  ! (Y1*Y3 .LE. ZERO)
       y2    = y3
       x2    = x3
     ELSE
       y1    = y3
       x1    = x3
     END IF
     IF (ABS(x2-x1) <= eps*x1) GO TO 40
   END DO
   CALL pusherr (0002, 'CALCG4')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CONVERGED ; RETURN **********************************************
   
   40    x3 = 0.5*(x1+x2)
   y3 = funcg4a (x3)
   
   ! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
   
   50    CONTINUE
   IF (molal(1) > tiny .AND. molal(5) > tiny) THEN
     CALL calchs4 (molal(1), molal(5), zero, delta)
     molal(1) = molal(1) - delta                     ! H+   EFFECT
     molal(5) = molal(5) - delta                     ! SO4  EFFECT
     molal(6) = delta                                ! HSO4 EFFECT
   END IF
   
END SUBROUTINE  calcg4
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE FUNCG4A
!! *** CASE G4
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
DOUBLE PRECISION FUNCTION funcg4a (x)
   USE Isorropia_Module
   USE caseg
   implicit none
   
   double precision, INTENT(IN)                         :: x
   DOUBLE PRECISION :: nai, nh4i, no3i, hi,ohi,bb,cc,dd,cli,so4i
   integer :: i,islv
   
   !LFR  COMMON /caseg/ chi1, chi2, chi3, chi4, chi5, chi6, lamda,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7
   
   ! *** SETUP PARAMETERS ************************************************
   
   psi6   = x
   psi1   = chi1
   frst   = .true.
   calain = .true.
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     
     a1  = xk5 *(water/gama(2))**3.0
     a2  = xk7 *(water/gama(4))**3.0
     a4  = (xk2/xkw)*r*temp*(gama(10)/gama(5))**2.0
     a5  = xk4 *r*temp*(water/gama(10))**2.0
     a6  = xk3 *r*temp*(water/gama(11))**2.0
     
   !  CALCULATE DISSOCIATION QUANTITIES
     
     IF (chi5 >= tiny) THEN
       psi5 = psi6*chi5/(a6/a5*(chi6-psi6) + psi6)
     ELSE
       psi5 = tiny
     END IF
     
   !CC      IF(CHI4.GT.TINY) THEN
     IF(w(2) > tiny) THEN       ! Accounts for NH3 evaporation
       bb   =-(chi4 + psi6 + psi5 + 1.d0/a4)
       cc   = chi4*(psi5+psi6) - 2.d0*psi2/a4
       dd   = MAX(bb*bb-4.d0*cc, 0.d0)       ! US patch 12/20/01
       psi4 =0.5D0*(-bb - SQRT(dd))
     ELSE
       psi4 = tiny
     END IF
     
   !  CALCULATE CONCENTRATIONS
     
     nh4i = 2.0*psi2 + psi4
     cli  = psi6
     so4i = psi2 + psi1
     no3i = psi5
     nai  = 2.0D0*psi1
     
     CALL calcph(2.d0*so4i+no3i+cli-nai-nh4i+organion, hi, ohi)
     
   ! *** Na2SO4 DISSOLUTION
     
     IF (chi1 > tiny .AND. water > tiny) THEN        ! PSI1
       CALL poly3 (psi2, zero, -a1/4.d0, psi1, islv)
       IF (islv == 0) THEN
         psi1 = MIN (psi1, chi1)
       ELSE
         psi1 = zero
       END IF
     ELSE
       psi1 = zero
     END IF
     
   ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
     
     molal (1) = hi
     molal (2) = nai
     molal (3) = nh4i
     molal (4) = cli
     molal (5) = so4i
     molal (6) = zero
     molal (7) = no3i
     
   ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
     
     gnh3      = MAX(chi4 - psi4, tiny)
     ghno3     = MAX(chi5 - psi5, tiny)
     ghcl      = MAX(chi6 - psi6, tiny)
     
     cnh42s4   = zero
     cnh4no3   = zero
     cnh4cl    = zero
     
   ! *** CALCULATE MOLALR ARRAY, WATER AND ACTIVITIES **********************
     
     CALL calcmr
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   
   ! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
   
   20    funcg4a = molal(1)*molal(4)/ghcl/a6 - one
   !CC         FUNCG4A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
   
END FUNCTION funcg4a
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCG3
!! *** CASE G3
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!!     2. LIQUID & SOLID PHASE ARE BOTH POSSIBLE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcg3
   USE Isorropia_Module
   implicit none
   
   EXTERNAL calcg1a, calcg4
   integer :: i
   
   ! *** REGIME DEPENDS ON THE EXISTANCE OF WATER AND OF THE RH ************
   
   IF (w(4) > tiny .AND. w(5) > tiny) THEN ! NO3,CL EXIST, WATER POSSIBLE
     scase = 'G3 ; SUBCASE 1'
     CALL calcg3a
     scase = 'G3 ; SUBCASE 1'
   ELSE                                      ! NO3, CL NON EXISTANT
     scase = 'G1 ; SUBCASE 1'
     CALL calcg1a
     scase = 'G1 ; SUBCASE 1'
   END IF
   
   IF (water <= tiny) THEN
     IF (rh < drmg3) THEN        ! ONLY SOLIDS
       water = tiny
       DO  i=1,nions
         molal(i) = zero
       END DO
       CALL calcg1a
       scase = 'G3 ; SUBCASE 2'
       RETURN
     ELSE
       scase = 'G3 ; SUBCASE 3'  ! MDRH REGION (NA2SO4, NH42S4)
       CALL calcmdrh (rh, drmg3, drnh42s4, calcg1a, calcg4)
       scase = 'G3 ; SUBCASE 3'
     END IF
   END IF
   
END SUBROUTINE  calcg3
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCG3A
!! *** CASE G3 ; SUBCASE 1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcg3a
   USE Isorropia_Module
   USE caseg
   IMPLICIT NONE
   DOUBLE PRECISION :: delta
   DOUBLE PRECISION :: dx, psi6hi, psi6lo, x1, x2, x3, y1, y2, y3
   DOUBLE PRECISION, EXTERNAL :: funcg3a
   INTEGER :: i,islv
   !LFR  COMMON /caseg/ chi1, chi2, chi3, chi4, chi5, chi6, lamda,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7
   
   ! *** SETUP PARAMETERS ************************************************
   
   calaou = .true.
   chi1   = 0.5*w(1)
   chi2   = MAX (w(2)-chi1, zero)
   chi3   = zero
   chi4   = MAX (w(3)-2.d0*chi2, zero)
   chi5   = w(4)
   chi6   = w(5)
   
   psi6lo = tiny
   psi6hi = chi6-tiny    ! MIN(CHI6-TINY, CHI4)
   
   water  = tiny
   
   ! *** INITIAL VALUES FOR BISECTION ************************************
   
   x1 = psi6lo
   y1 = funcg3a(x1)
   IF (chi6 <= tiny) GO TO 50
   !CC      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY .OR. WATER .LE. TINY) GOTO 50
   !CC      IF (WATER .LE. TINY) RETURN                    ! No water
   
   ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
   
   dx = (psi6hi-psi6lo)/FLOAT(ndiv)
   DO  i=1,ndiv
     x2  = x1+dx
     y2  = funcg3a (x2)
     
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y2) < zero) GO TO 20  ! (Y1*Y2.LT.ZERO)
     x1  = x2
     y1  = y2
   END DO
   
   ! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
   
   IF (ABS(y2) > eps) y2 = funcg3a (psi6lo)
   GO TO 50
   
   ! *** PERFORM BISECTION ***********************************************
   
   20    DO  i=1,maxit
     x3 = 0.5*(x1+x2)
     y3 = funcg3a (x3)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y3) <= zero) THEN  ! (Y1*Y3 .LE. ZERO)
       y2    = y3
       x2    = x3
     ELSE
       y1    = y3
       x1    = x3
     END IF
     IF (ABS(x2-x1) <= eps*x1) GO TO 40
   END DO
   CALL pusherr (0002, 'CALCG3A')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CONVERGED ; RETURN **********************************************
   
   40    x3 = 0.5*(x1+x2)
   y3 = funcg3a (x3)
   
   ! *** FINAL CALCULATIONS *************************************************
   
   50    CONTINUE
   
   ! *** Na2SO4 DISSOLUTION
   
   IF (chi1 > tiny .AND. water > tiny) THEN        ! PSI1
     CALL poly3 (psi2, zero, -a1/4.d0, psi1, islv)
     IF (islv == 0) THEN
       psi1 = MIN (psi1, chi1)
     ELSE
       psi1 = zero
     END IF
   ELSE
     psi1 = zero
   END IF
   molal(2) = 2.0D0*psi1               ! Na+  EFFECT
   molal(5) = molal(5) + psi1          ! SO4  EFFECT
   cna2so4  = MAX(chi1 - psi1, zero)   ! NA2SO4(s) depletion
   
   ! *** HSO4 equilibrium
   
   IF (molal(1) > tiny .AND. molal(5) > tiny) THEN
     CALL calchs4 (molal(1), molal(5), zero, delta)
     molal(1) = molal(1) - delta                     ! H+   EFFECT
     molal(5) = molal(5) - delta                     ! SO4  EFFECT
     molal(6) = delta                                ! HSO4 EFFECT
   END IF
   
END SUBROUTINE  calcg3a
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE FUNCG3A
!! *** CASE G3 ; SUBCASE 1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
DOUBLE PRECISION FUNCTION funcg3a (x)
   USE Isorropia_Module
   USE caseg
   implicit none
   
   double precision, INTENT(IN)                         :: x
   DOUBLE PRECISION :: psi20, smin, hi, ohi,bb,cc,dd
   integer :: i,islv
   
   !LFR  COMMON /caseg/ chi1, chi2, chi3, chi4, chi5, chi6, lamda,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7
   
   ! *** SETUP PARAMETERS ************************************************
   
   psi6   = x
   psi2   = chi2
   frst   = .true.
   calain = .true.
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     
     a1  = xk5 *(water/gama(2))**3.0
     a2  = xk7 *(water/gama(4))**3.0
     a4  = (xk2/xkw)*r*temp*(gama(10)/gama(5))**2.0
     a5  = xk4 *r*temp*(water/gama(10))**2.0
     a6  = xk3 *r*temp*(water/gama(11))**2.0
     
   !  CALCULATE DISSOCIATION QUANTITIES
     
     IF (chi5 >= tiny) THEN
       psi5 = psi6*chi5/(a6/a5*(chi6-psi6) + psi6)
     ELSE
       psi5 = tiny
     END IF
     
   !CC      IF(CHI4.GT.TINY) THEN
     IF(w(2) > tiny) THEN       ! Accounts for NH3 evaporation
       bb   =-(chi4 + psi6 + psi5 + 1.d0/a4)
       cc   = chi4*(psi5+psi6) - 2.d0*psi2/a4
       dd   = MAX(bb*bb-4.d0*cc, 0.d0)       ! US patch 12/20/01
       psi4 =0.5D0*(-bb - SQRT(dd))
     ELSE
       psi4 = tiny
     END IF
     
     IF (chi2 > tiny .AND. water > tiny) THEN
       CALL poly3 (psi4, psi4*psi4/4.d0, -a2/4.d0, psi20, islv)
       IF (islv == 0) psi2 = MIN (psi20, chi2)
     END IF
     
   ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
     
     molal (2) = zero                                ! Na
     molal (3) = 2.0*psi2 + psi4                     ! NH4I
     molal (4) = psi6                                ! CLI
     molal (5) = psi2                                ! SO4I
     molal (6) = zero                                ! HSO4
     molal (7) = psi5                                ! NO3I
     
     smin  = 2.d0*molal(5)+molal(7)+molal(4)-molal(2)-molal(3) + organion
     CALL calcph (smin, hi, ohi)
     molal (1) = hi
     
     gnh3      = MAX(chi4 - psi4, tiny)              ! Gas NH3
     ghno3     = MAX(chi5 - psi5, tiny)              ! Gas HNO3
     ghcl      = MAX(chi6 - psi6, tiny)              ! Gas HCl
     
     cnh42s4   = chi2 - psi2                         ! Solid (NH4)2SO4
     cnh4no3   = zero                                ! Solid NH4NO3
     cnh4cl    = zero                                ! Solid NH4Cl
     
     CALL calcmr                                     ! Water content
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   
   ! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
   
   20    funcg3a = molal(1)*molal(4)/ghcl/a6 - one
   !CC         FUNCG3A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
   
END FUNCTION funcg3a
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCG2
!! *** CASE G2
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!!     2. LIQUID & SOLID PHASE ARE BOTH POSSIBLE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcg2
   USE Isorropia_Module
   implicit none
   
   EXTERNAL calcg1a, calcg3a, calcg4
   integer :: i
   
   ! *** REGIME DEPENDS ON THE EXISTANCE OF NITRATES ***********************
   
   IF (w(4) > tiny) THEN        ! NO3 EXISTS, WATER POSSIBLE
     scase = 'G2 ; SUBCASE 1'
     CALL calcg2a
     scase = 'G2 ; SUBCASE 1'
   ELSE                          ! NO3 NON EXISTANT, WATER NOT POSSIBLE
     scase = 'G1 ; SUBCASE 1'
     CALL calcg1a
     scase = 'G1 ; SUBCASE 1'
   END IF
   
   ! *** REGIME DEPENDS ON THE EXISTANCE OF WATER AND OF THE RH ************
   
   IF (water <= tiny) THEN
     IF (rh < drmg2) THEN             ! ONLY SOLIDS
       water = tiny
       DO  i=1,nions
         molal(i) = zero
       END DO
       CALL calcg1a
       scase = 'G2 ; SUBCASE 2'
     ELSE
       IF (w(5) > tiny) THEN
         scase = 'G2 ; SUBCASE 3'    ! MDRH (NH4CL, NA2SO4, NH42S4)
         CALL calcmdrh (rh, drmg2, drnh4cl, calcg1a, calcg3a)
         scase = 'G2 ; SUBCASE 3'
       END IF
       IF (water <= tiny .AND. rh >= drmg3) THEN
         scase = 'G2 ; SUBCASE 4'    ! MDRH (NA2SO4, NH42S4)
         CALL calcmdrh (rh, drmg3, drnh42s4, calcg1a, calcg4)
         scase = 'G2 ; SUBCASE 4'
       ELSE
         water = tiny
         DO  i=1,nions
           molal(i) = zero
         END DO
         CALL calcg1a
         scase = 'G2 ; SUBCASE 2'
       END IF
     END IF
   END IF
   
END SUBROUTINE  calcg2
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCG2A
!! *** CASE G2 ; SUBCASE 1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcg2a
   USE Isorropia_Module
   USE caseg
   implicit none
   
   DOUBLE PRECISION :: delta,x1,x2,x3,y1,y2,y3,psi6lo,psi6hi,dx
   double precision, external :: funcg2a
   integer :: i,islv
   
   !LFR  COMMON /caseg/ chi1, chi2, chi3, chi4, chi5, chi6, lamda,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7
   
   ! *** SETUP PARAMETERS ************************************************
   
   calaou = .true.
   chi1   = 0.5*w(1)
   chi2   = MAX (w(2)-chi1, zero)
   chi3   = zero
   chi4   = MAX (w(3)-2.d0*chi2, zero)
   chi5   = w(4)
   chi6   = w(5)
   
   psi6lo = tiny
   psi6hi = chi6-tiny
   
   water  = tiny
   
   ! *** INITIAL VALUES FOR BISECTION ************************************
   
   x1 = psi6lo
   y1 = funcg2a (x1)
   IF (chi6 <= tiny) GO TO 50
   !CC      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50
   !CC      IF (WATER .LE. TINY) GOTO 50               ! No water
   
   ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
   
   dx = (psi6hi-psi6lo)/FLOAT(ndiv)
   DO  i=1,ndiv
     x2 = x1+dx
     y2 = funcg2a (x2)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y2) < zero) GO TO 20  ! (Y1*Y2.LT.ZERO)
     x1 = x2
     y1 = y2
   END DO
   
   ! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
   
   IF (ABS(y2) > eps) water = tiny
   GO TO 50
   
   ! *** PERFORM BISECTION ***********************************************
   
   20    DO  i=1,maxit
     x3 = 0.5*(x1+x2)
     y3 = funcg2a (x3)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y3) <= zero) THEN  ! (Y1*Y3 .LE. ZERO)
       y2    = y3
       x2    = x3
     ELSE
       y1    = y3
       x1    = x3
     END IF
     IF (ABS(x2-x1) <= eps*x1) GO TO 40
   END DO
   CALL pusherr (0002, 'CALCG2A')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CONVERGED ; RETURN **********************************************
   
   40    x3 = 0.5*(x1+x2)
   IF (x3 <= tiny2) THEN   ! PRACTICALLY NO NITRATES, SO DRY SOLUTION
     water = tiny
   ELSE
     y3 = funcg2a (x3)
   END IF
   
   ! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
   
   50    CONTINUE
   
   ! *** Na2SO4 DISSOLUTION
   
   IF (chi1 > tiny .AND. water > tiny) THEN        ! PSI1
     CALL poly3 (psi2, zero, -a1/4.d0, psi1, islv)
     IF (islv == 0) THEN
       psi1 = MIN (psi1, chi1)
     ELSE
       psi1 = zero
     END IF
   ELSE
     psi1 = zero
   END IF
   molal(2) = 2.0D0*psi1               ! Na+  EFFECT
   molal(5) = molal(5) + psi1          ! SO4  EFFECT
   cna2so4  = MAX(chi1 - psi1, zero)   ! NA2SO4(s) depletion
   
   ! *** HSO4 equilibrium
   
   IF (molal(1) > tiny .AND. molal(5) > tiny) THEN
     CALL calchs4 (molal(1), molal(5), zero, delta)
     molal(1) = molal(1) - delta     ! H+   AFFECT
     molal(5) = molal(5) - delta     ! SO4  AFFECT
     molal(6) = delta                ! HSO4 AFFECT
   END IF
   
END SUBROUTINE  calcg2a
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE FUNCG2A
!! *** CASE G2
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
DOUBLE PRECISION FUNCTION funcg2a (x)
   USE Isorropia_Module
   USE caseg
   implicit none
   
   double precision, INTENT(IN)                         :: x
   DOUBLE PRECISION :: psi20,smin,hi,ohi,deno,delt,bb,cc,dd,psi31,psi32
   integer :: i,islv
   
   !LFR  COMMON /caseg/ chi1, chi2, chi3, chi4, chi5, chi6, lamda,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7
   
   ! *** SETUP PARAMETERS ************************************************
   
   psi6   = x
   psi2   = chi2
   psi3   = zero
   frst   = .true.
   calain = .true.
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     
     a1  = xk5 *(water/gama(2))**3.0
     a2  = xk7 *(water/gama(4))**3.0
     a4  = (xk2/xkw)*r*temp*(gama(10)/gama(5))**2.0
     a5  = xk4 *r*temp*(water/gama(10))**2.0
     a6  = xk3 *r*temp*(water/gama(11))**2.0
     
     deno = MAX(chi6-psi6-psi3, zero)
     psi5 = chi5/((a6/a5)*(deno/psi6) + one)
     
     psi4 = MIN(psi5+psi6,chi4)
     
     IF (chi2 > tiny .AND. water > tiny) THEN
       CALL poly3 (psi4, psi4*psi4/4.d0, -a2/4.d0, psi20, islv)
       IF (islv == 0) psi2 = MIN (psi20, chi2)
     END IF
     
   ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
     
     molal (2) = zero                             ! NA
     molal (3) = 2.0*psi2 + psi4                  ! NH4I
     molal (4) = psi6                             ! CLI
     molal (5) = psi2                             ! SO4I
     molal (6) = zero                             ! HSO4
     molal (7) = psi5                             ! NO3I
     
   !CC      MOLAL (1) = MAX(CHI5 - PSI5, TINY)*A5/PSI5   ! HI
     smin = 2.d0*molal(5)+molal(7)+molal(4)-molal(2)-molal(3) + organion
     CALL calcph (smin, hi, ohi)
     molal (1) = hi
     
   ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
     
     gnh3      = MAX(chi4 - psi4, tiny)
     ghno3     = MAX(chi5 - psi5, tiny)
     ghcl      = MAX(chi6 - psi6, tiny)
     
     cnh42s4   = MAX(chi2 - psi2, zero)
     cnh4no3   = zero
     
   ! *** NH4Cl(s) calculations
     
     a3   = xk6 /(r*temp*r*temp)
     IF (gnh3*ghcl > a3) THEN
       delt = MIN(gnh3, ghcl)
       bb = -(gnh3+ghcl)
       cc = gnh3*ghcl-a3
       dd = MAX(bb*bb-4.d0*cc, 0.d0)       ! US patch 12/20/01
       psi31 = 0.5D0*(-bb + SQRT(dd))
       psi32 = 0.5D0*(-bb - SQRT(dd))
       IF (delt-psi31 > zero .AND. psi31 > zero) THEN
         psi3 = psi31
       ELSE IF (delt-psi32 > zero .AND. psi32 > zero) THEN
         psi3 = psi32
       ELSE
         psi3 = zero
       END IF
     ELSE
       psi3 = zero
     END IF
     
   ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
     
     gnh3    = MAX(gnh3 - psi3, tiny)
     ghcl    = MAX(ghcl - psi3, tiny)
     cnh4cl  = psi3
     
   ! *** CALCULATE MOLALR ARRAY, WATER AND ACTIVITIES **********************
     
     CALL calcmr
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   
   ! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
   
   20    IF (chi4 <= tiny) THEN
     funcg2a = molal(1)*molal(4)/ghcl/a6 - one
   ELSE
     funcg2a = molal(3)*molal(4)/ghcl/gnh3/a6/a4 - one
   END IF
   
END FUNCTION funcg2a
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCG1
!! *** CASE G1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!!     2. SOLID AEROSOL ONLY
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3, NH4CL, NA2SO4
!!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
!!     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!!     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCG1A)
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcg1
   USE Isorropia_Module
   implicit none
   
   EXTERNAL calcg1a, calcg2a
   
   ! *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************
   
   IF (rh < drmg1) THEN
     scase = 'G1 ; SUBCASE 1'
     CALL calcg1a              ! SOLID PHASE ONLY POSSIBLE
     scase = 'G1 ; SUBCASE 1'
   ELSE
     scase = 'G1 ; SUBCASE 2'  ! LIQUID & SOLID PHASE POSSIBLE
     CALL calcmdrh (rh, drmg1, drnh4no3, calcg1a, calcg2a)
     scase = 'G1 ; SUBCASE 2'
   END IF
   
END SUBROUTINE  calcg1
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCG1A
!! *** CASE G1 ; SUBCASE 1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0)
!!     2. SOLID AEROSOL ONLY
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3
!!     SOLID (NH4)2SO4 IS CALCULATED FROM THE SULFATES, WHILE NH4NO3
!!     IS CALCULATED FROM NH3-HNO3 EQUILIBRIUM. 'ZE' IS THE AMOUNT OF
!!     NH4NO3 THAT VOLATIZES WHEN ALL POSSILBE NH4NO3 IS INITIALLY IN
!!     THE SOLID PHASE.
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcg1a
   USE Isorropia_Module
   implicit none
   
   DOUBLE PRECISION :: lamda,lamda1, lamda2, kapa, kapa1, kapa2
   double precision :: alf,bet,gam,a1,a2,theta1,theta2,bb,cc,dd,sqdd
   double precision :: dd1,dd2,sqdd1,sqdd2,rtsq
   
   
   ! *** CALCULATE NON VOLATILE SOLIDS ***********************************
   
   cna2so4 = 0.5*w(1)
   cnh42s4 = w(2) - cna2so4
   
   ! *** CALCULATE VOLATILE SPECIES **************************************
   
   alf     = w(3) - 2.0*cnh42s4
   bet     = w(5)
   gam     = w(4)
   
   rtsq    = r*temp*r*temp
   a1      = xk6/rtsq
   a2      = xk10/rtsq
   
   theta1  = gam - bet*(a2/a1)
   theta2  = a2/a1
   
   ! QUADRATIC EQUATION SOLUTION
   
   bb      = (theta1-alf-bet*(one+theta2))/(one+theta2)
   cc      = (alf*bet-a1-bet*theta1)/(one+theta2)
   dd      = MAX( bb*bb - 4.0D0*cc, 0.d0 )  ! US patch 12/20/01
   IF (dd < zero) GO TO 100   ! Solve each reaction seperately
   
   ! TWO ROOTS FOR KAPA, CHECK AND SEE IF ANY VALID
   
   sqdd    = SQRT(dd)
   kapa1   = 0.5D0*(-bb+sqdd)
   kapa2   = 0.5D0*(-bb-sqdd)
   lamda1  = theta1 + theta2*kapa1
   lamda2  = theta1 + theta2*kapa2
   
   IF (kapa1 >= zero .AND. lamda1 >= zero) THEN
     IF (alf-kapa1-lamda1 >= zero .AND.  &
           bet-kapa1 >= zero .AND. gam-lamda1 >= zero) THEN
       kapa = kapa1
       lamda= lamda1
       GO TO 200
     END IF
   END IF
   
   IF (kapa2 >= zero .AND. lamda2 >= zero) THEN
     IF (alf-kapa2-lamda2 >= zero .AND.  &
           bet-kapa2 >= zero .AND. gam-lamda2 >= zero) THEN
       kapa = kapa2
       lamda= lamda2
       GO TO 200
     END IF
   END IF
   
   ! SEPERATE SOLUTION OF NH4CL & NH4NO3 EQUILIBRIA
   
   100   kapa  = zero
   lamda = zero
   dd1   = (alf+bet)*(alf+bet) - 4.0D0*(alf*bet-a1)
   dd2   = (alf+gam)*(alf+gam) - 4.0D0*(alf*gam-a2)
   
   ! NH4CL EQUILIBRIUM
   
   IF (dd1 >= zero) THEN
     sqdd1 = SQRT(dd1)
     kapa1 = 0.5D0*(alf+bet + sqdd1)
     kapa2 = 0.5D0*(alf+bet - sqdd1)
     
     IF (kapa1 >= zero .AND. kapa1 <= MIN(alf,bet)) THEN
       kapa = kapa1
     ELSE IF (kapa2 >= zero .AND. kapa2 <= MIN(alf,bet)) THEN
       kapa = kapa2
     ELSE
       kapa = zero
     END IF
   END IF
   
   ! NH4NO3 EQUILIBRIUM
   
   IF (dd2 >= zero) THEN
     sqdd2 = SQRT(dd2)
     lamda1= 0.5D0*(alf+gam + sqdd2)
     lamda2= 0.5D0*(alf+gam - sqdd2)
     
     IF (lamda1 >= zero .AND. lamda1 <= MIN(alf,gam)) THEN
       lamda = lamda1
     ELSE IF (lamda2 >= zero .AND. lamda2 <= MIN(alf,gam)) THEN
       lamda = lamda2
     ELSE
       lamda = zero
     END IF
   END IF
   
   ! IF BOTH KAPA, LAMDA ARE > 0, THEN APPLY EXISTANCE CRITERION
   
   IF (kapa > zero .AND. lamda > zero) THEN
     IF (bet < lamda/theta1) THEN
       kapa = zero
     ELSE
       lamda= zero
     END IF
   END IF
   
   ! *** CALCULATE COMPOSITION OF VOLATILE SPECIES ***********************
   
   200   CONTINUE
   cnh4no3 = lamda
   cnh4cl  = kapa
   
   gnh3    = MAX(alf - kapa - lamda, zero)
   ghno3   = MAX(gam - lamda, zero)
   ghcl    = MAX(bet - kapa, zero)
   
END SUBROUTINE  calcg1a

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCH6
!! *** CASE H6
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calch6
   USE Isorropia_Module
   use solut
   implicit none
   
   double precision :: delta,frna,psi6lo,psi6hi,x1,y1,dx,x2,y2,x3,y3
   double precision, external :: funch6a
   integer :: i
   
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   calaou = .true.
   chi1   = w(2)                                ! CNA2SO4
   chi2   = zero                                ! CNH42S4
   chi3   = zero                                ! CNH4CL
   frna   = MAX (w(1)-2.d0*chi1, zero)
   chi8   = MIN (frna, w(4))                    ! CNANO3
   chi4   = w(3)                                ! NH3(g)
   chi5   = MAX (w(4)-chi8, zero)               ! HNO3(g)
   chi7   = MIN (MAX(frna-chi8, zero), w(5))    ! CNACL
   chi6   = MAX (w(5)-chi7, zero)               ! HCL(g)
   
   psi6lo = tiny
   psi6hi = chi6-tiny    ! MIN(CHI6-TINY, CHI4)
   
   ! *** INITIAL VALUES FOR BISECTION ************************************
   
   x1 = psi6lo
   y1 = funch6a (x1)
   IF (ABS(y1) <= eps .OR. chi6 <= tiny) GO TO 50
   
   ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
   
   dx = (psi6hi-psi6lo)/FLOAT(ndiv)
   DO  i=1,ndiv
     x2 = x1+dx
     y2 = funch6a (x2)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y2) < zero) GO TO 20  ! (Y1*Y2.LT.ZERO)
     x1 = x2
     y1 = y2
   END DO
   
   ! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
   
   IF (ABS(y2) > eps) y2 = funch6a (psi6lo)
   GO TO 50
   
   ! *** PERFORM BISECTION ***********************************************
   
   20    DO  i=1,maxit
     x3 = 0.5*(x1+x2)
     y3 = funch6a (x3)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y3) <= zero) THEN  ! (Y1*Y3 .LE. ZERO)
       y2    = y3
       x2    = x3
     ELSE
       y1    = y3
       x1    = x3
     END IF
     IF (ABS(x2-x1) <= eps*x1) GO TO 40
   END DO
   CALL pusherr (0002, 'CALCH6')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CONVERGED ; RETURN **********************************************
   
   40    x3 = 0.5*(x1+x2)
   y3 = funch6a (x3)
   
   ! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
   
   50    CONTINUE
   IF (molal(1) > tiny .AND. molal(5) > tiny) THEN
     CALL calchs4 (molal(1), molal(5), zero, delta)
     molal(1) = molal(1) - delta                     ! H+   EFFECT
     molal(5) = molal(5) - delta                     ! SO4  EFFECT
     molal(6) = delta                                ! HSO4 EFFECT
   END IF
   
END SUBROUTINE  calch6
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE FUNCH6A
!! *** CASE H6
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
DOUBLE PRECISION FUNCTION funch6a (x)
   USE Isorropia_Module
   use solut
   implicit none
   
   double precision, INTENT(IN)                         :: x
   double precision :: smin,hi,ohi,a9,bb,cc,dd
   integer :: i
   
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   psi6   = x
   psi1   = chi1
   psi2   = zero
   psi3   = zero
   psi7   = chi7
   psi8   = chi8
   frst   = .true.
   calain = .true.
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     
     a1  = xk5 *(water/gama(2))**3.0
     a4  = (xk2/xkw)*r*temp*(gama(10)/gama(5))**2.0
     a5  = xk4 *r*temp*(water/gama(10))**2.0
     a6  = xk3 *r*temp*(water/gama(11))**2.0
     a7  = xk8 *(water/gama(1))**2.0
     a8  = xk9 *(water/gama(3))**2.0
     a9  = xk1*water/gama(7)*(gama(8)/gama(7))**2.
     
   !  CALCULATE DISSOCIATION QUANTITIES
     
     psi5 = chi5*(psi6+psi7) - a6/a5*psi8*(chi6-psi6-psi3)
     psi5 = psi5/(a6/a5*(chi6-psi6-psi3) + psi6 + psi7)
     psi5 = MAX(psi5, tiny)
     
     IF (chi1 > tiny .AND. water > tiny) THEN  ! First try 3rd order soln
       bb   =-(chi4 + psi6 + psi5 + 1.d0/a4)
       cc   = chi4*(psi5+psi6)
       dd   = MAX(bb*bb-4.d0*cc, 0.d0)       ! US patch 12/20/01
       psi4 =0.5D0*(-bb - SQRT(dd))
       psi4 = MIN(psi4,chi4)
     ELSE
       psi4 = tiny
     END IF
     
   ! *** CALCULATE SPECIATION ********************************************
     
     molal (2) = psi8 + psi7 + 2.d0*psi1               ! NAI
     molal (3) = psi4                                  ! NH4I
     molal (4) = psi6 + psi7                           ! CLI
     molal (5) = psi2 + psi1                           ! SO4I
     molal (6) = zero                                  ! HSO4I
     molal (7) = psi5 + psi8                           ! NO3I
     
     smin = 2.d0*molal(5)+molal(7)+molal(4)-molal(2)-molal(3) + organion
     CALL calcph (smin, hi, ohi)
     molal (1) = hi
     
     gnh3      = MAX(chi4 - psi4, tiny)
     ghno3     = MAX(chi5 - psi5, tiny)
     ghcl      = MAX(chi6 - psi6, tiny)
     
     cnh42s4   = zero
     cnh4no3   = zero
     cnacl     = MAX(chi7 - psi7, zero)
     cnano3    = MAX(chi8 - psi8, zero)
     cna2so4   = MAX(chi1 - psi1, zero)
     
     CALL calcmr                                    ! Water content
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   
   ! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
   
   20    funch6a = molal(3)*molal(4)/ghcl/gnh3/a6/a4 - one
   
END FUNCTION funch6a
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCH5
!! *** CASE H5
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calch5
   USE Isorropia_Module
   use solut
   implicit none
   
   double precision :: delta,frna,psi6lo,psi6hi,x1,y1,dx,x2,y2,x3,y3
   double precision, external :: funch5a
   integer :: i
   
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** REGIME DEPENDS ON THE EXISTANCE OF NITRATES ***********************
   
   IF (w(4) <= tiny .AND. w(5) <= tiny) THEN
     scase = 'H5'
     CALL calch1a
     scase = 'H5'
     RETURN
   END IF
   
   ! *** SETUP PARAMETERS ************************************************
   
   calaou = .true.
   chi1   = w(2)                                ! CNA2SO4
   chi2   = zero                                ! CNH42S4
   chi3   = zero                                ! CNH4CL
   frna   = MAX (w(1)-2.d0*chi1, zero)
   chi8   = MIN (frna, w(4))                    ! CNANO3
   chi4   = w(3)                                ! NH3(g)
   chi5   = MAX (w(4)-chi8, zero)               ! HNO3(g)
   chi7   = MIN (MAX(frna-chi8, zero), w(5))    ! CNACL
   chi6   = MAX (w(5)-chi7, zero)               ! HCL(g)
   
   psi6lo = tiny
   psi6hi = chi6-tiny    ! MIN(CHI6-TINY, CHI4)
   
   ! *** INITIAL VALUES FOR BISECTION ************************************
   
   x1 = psi6lo
   y1 = funch5a (x1)
   IF (ABS(y1) <= eps .OR. chi6 <= tiny) GO TO 50
   
   ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
   
   dx = (psi6hi-psi6lo)/FLOAT(ndiv)
   DO  i=1,ndiv
     x2 = x1+dx
     y2 = funch5a (x2)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y2) < zero) GO TO 20  ! (Y1*Y2.LT.ZERO)
     x1 = x2
     y1 = y2
   END DO
   
   ! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
   
   IF (ABS(y2) > eps) y2 = funch5a (psi6lo)
   GO TO 50
   
   ! *** PERFORM BISECTION ***********************************************
   
   20    DO  i=1,maxit
     x3 = 0.5*(x1+x2)
     y3 = funch5a (x3)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y3) <= zero) THEN  ! (Y1*Y3 .LE. ZERO)
       y2    = y3
       x2    = x3
     ELSE
       y1    = y3
       x1    = x3
     END IF
     IF (ABS(x2-x1) <= eps*x1) GO TO 40
   END DO
   CALL pusherr (0002, 'CALCH5')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CONVERGED ; RETURN **********************************************
   
   40    x3 = 0.5*(x1+x2)
   y3 = funch5a (x3)
   
   ! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
   
   50    CONTINUE
   IF (molal(1) > tiny .AND. molal(5) > tiny) THEN
     CALL calchs4 (molal(1), molal(5), zero, delta)
     molal(1) = molal(1) - delta                     ! H+   EFECT
     molal(5) = molal(5) - delta                     ! SO4  EFFECT
     molal(6) = delta                                ! HSO4 EFFECT
   END IF
   
END SUBROUTINE  calch5
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE FUNCH5A
!! *** CASE H5
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : NONE
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
DOUBLE PRECISION FUNCTION funch5a (x)
   USE Isorropia_Module
   use solut
   implicit none
   
   double precision, INTENT(IN)                         :: x
   double precision :: smin,hi,ohi,aa,bb,cc,a9,dd
   integer :: i,islv
   
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   psi6   = x
   psi1   = chi1
   psi2   = zero
   psi3   = zero
   psi7   = chi7
   psi8   = chi8
   frst   = .true.
   calain = .true.
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     
     a1  = xk5 *(water/gama(2))**3.0
     a4  = (xk2/xkw)*r*temp*(gama(10)/gama(5))**2.0
     a5  = xk4 *r*temp*(water/gama(10))**2.0
     a6  = xk3 *r*temp*(water/gama(11))**2.0
     a7  = xk8 *(water/gama(1))**2.0
     a8  = xk9 *(water/gama(3))**2.0
     a9  = xk1*water/gama(7)*(gama(8)/gama(7))**2.
     
   !  CALCULATE DISSOCIATION QUANTITIES
     
     psi5 = chi5*(psi6+psi7) - a6/a5*psi8*(chi6-psi6-psi3)
     psi5 = psi5/(a6/a5*(chi6-psi6-psi3) + psi6 + psi7)
     psi5 = MAX(psi5, tiny)
     
     IF (chi1 > tiny .AND. water > tiny) THEN  ! First try 3rd order soln
       bb   =-(chi4 + psi6 + psi5 + 1.d0/a4)
       cc   = chi4*(psi5+psi6)
       dd   = MAX(bb*bb-4.d0*cc, 0.d0)       ! US patch 12/20/01
       psi4 =0.5D0*(-bb - SQRT(dd))
       psi4 = MIN(psi4,chi4)
     ELSE
       psi4 = tiny
     END IF
     
     IF (chi1 > tiny .AND. water > tiny) THEN     ! NA2SO4 DISSOLUTION
       aa = psi7+psi8
       bb = aa*aa
       cc =-a1/4.d0
       CALL poly3 (aa, bb, cc, psi1, islv)
       IF (islv == 0) THEN
         psi1 = MIN (psi1, chi1)
       ELSE
         psi1 = zero
       END IF
     END IF
     
   ! *** CALCULATE SPECIATION ********************************************
     
     molal (2) = psi8 + psi7 + 2.d0*psi1                ! NAI
     molal (3) = psi4                                   ! NH4I
     molal (4) = psi6 + psi7                            ! CLI
     molal (5) = psi2 + psi1                            ! SO4I
     molal (6) = zero
     molal (7) = psi5 + psi8                            ! NO3I
     
     smin = 2.d0*molal(5)+molal(7)+molal(4)-molal(2)-molal(3) + organion
     CALL calcph (smin, hi, ohi)
     molal (1) = hi
     
     gnh3      = MAX(chi4 - psi4, tiny)
     ghno3     = MAX(chi5 - psi5, tiny)
     ghcl      = MAX(chi6 - psi6, tiny)
     
     cnh42s4   = zero
     cnh4no3   = zero
     cnacl     = MAX(chi7 - psi7, zero)
     cnano3    = MAX(chi8 - psi8, zero)
     cna2so4   = MAX(chi1 - psi1, zero)
     
     CALL calcmr                               ! Water content
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   
   ! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
   
   20    funch5a = molal(3)*molal(4)/ghcl/gnh3/a6/a4 - one
   
END FUNCTION funch5a
   
!!=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCH4
!! *** CASE H4
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calch4
   USE Isorropia_Module
   use solut
   implicit none
         
   double precision :: delta,frna,psi6lo,psi6hi,x1,y1,dx,x2,y2,x3,y3
   double precision, external :: funch4a
   integer :: i
   
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** REGIME DEPENDS ON THE EXISTANCE OF NITRATES ***********************
   
   IF (w(4) <= tiny .AND. w(5) <= tiny) THEN
     scase = 'H4'
     CALL calch1a
     scase = 'H4'
     RETURN
   END IF
   
   ! *** SETUP PARAMETERS ************************************************
   
   calaou = .true.
   chi1   = w(2)                                ! CNA2SO4
   chi2   = zero                                ! CNH42S4
   chi3   = zero                                ! CNH4CL
   frna   = MAX (w(1)-2.d0*chi1, zero)
   chi8   = MIN (frna, w(4))                    ! CNANO3
   chi4   = w(3)                                ! NH3(g)
   chi5   = MAX (w(4)-chi8, zero)               ! HNO3(g)
   chi7   = MIN (MAX(frna-chi8, zero), w(5))    ! CNACL
   chi6   = MAX (w(5)-chi7, zero)               ! HCL(g)
   
   psi6lo = tiny
   psi6hi = chi6-tiny    ! MIN(CHI6-TINY, CHI4)
   
   ! *** INITIAL VALUES FOR BISECTION ************************************
   
   x1 = psi6lo
   y1 = funch4a (x1)
   IF (ABS(y1) <= eps .OR. chi6 <= tiny) GO TO 50
   
   ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
   
   dx = (psi6hi-psi6lo)/FLOAT(ndiv)
   DO  i=1,ndiv
     x2 = x1+dx
     y2 = funch4a (x2)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y2) < zero) GO TO 20  ! (Y1*Y2.LT.ZERO)
     x1 = x2
     y1 = y2
   END DO
   
   ! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
   
   IF (ABS(y2) > eps) y2 = funch4a (psi6lo)
   GO TO 50
   
   ! *** PERFORM BISECTION ***********************************************
   
   20    DO  i=1,maxit
     x3 = 0.5*(x1+x2)
     y3 = funch4a (x3)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y3) <= zero) THEN  ! (Y1*Y3 .LE. ZERO)
       y2    = y3
       x2    = x3
     ELSE
       y1    = y3
       x1    = x3
     END IF
     IF (ABS(x2-x1) <= eps*x1) GO TO 40
   END DO
   CALL pusherr (0002, 'CALCH4')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CONVERGED ; RETURN **********************************************
   
   40    x3 = 0.5*(x1+x2)
   y3 = funch4a (x3)
   
   ! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
   
   50    CONTINUE
   IF (molal(1) > tiny .AND. molal(5) > tiny) THEN
     CALL calchs4 (molal(1), molal(5), zero, delta)
     molal(1) = molal(1) - delta                      ! H+   EFFECT
     molal(5) = molal(5) - delta                      ! SO4  EFFECT
     molal(6) = delta                                 ! HSO4 EFFECT
   END IF
   
END SUBROUTINE  calch4
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE FUNCH4A
!! *** CASE H4
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
DOUBLE PRECISION FUNCTION funch4a (x)
   USE Isorropia_Module
   use solut
   implicit none
   
   double precision , INTENT(IN)                         :: x
   double precision :: smin,hi,ohi,aa,bb,cc,a9,dd,delt,psi31,psi32
   integer :: i,islv
   
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   psi6   = x
   psi1   = chi1
   psi2   = zero
   psi3   = zero
   psi7   = chi7
   psi8   = chi8
   frst   = .true.
   calain = .true.
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     
     a1  = xk5 *(water/gama(2))**3.0
     a4  = (xk2/xkw)*r*temp*(gama(10)/gama(5))**2.0
     a5  = xk4 *r*temp*(water/gama(10))**2.0
     a6  = xk3 *r*temp*(water/gama(11))**2.0
     a7  = xk8 *(water/gama(1))**2.0
     a8  = xk9 *(water/gama(3))**2.0
     a9  = xk1*water/gama(7)*(gama(8)/gama(7))**2.
     
   !  CALCULATE DISSOCIATION QUANTITIES
     
     psi5 = chi5*(psi6+psi7) - a6/a5*psi8*(chi6-psi6-psi3)
     psi5 = psi5/(a6/a5*(chi6-psi6-psi3) + psi6 + psi7)
     psi5 = MAX(psi5, tiny)
     
     IF (chi1 > tiny .AND. water > tiny) THEN  ! First try 3rd order soln
       bb   =-(chi4 + psi6 + psi5 + 1.d0/a4)
       cc   = chi4*(psi5+psi6)
       dd   = MAX(bb*bb-4.d0*cc, 0.d0)       ! US patch 12/20/01
       psi4 =0.5D0*(-bb - SQRT(dd))
       psi4 = MIN(psi4,chi4)
     ELSE
       psi4 = tiny
     END IF
     
     IF (chi1 > tiny .AND. water > tiny) THEN     ! NA2SO4 DISSOLUTION
       aa = psi7+psi8
       bb = aa*aa
       cc =-a1/4.d0
       CALL poly3 (aa, bb, cc, psi1, islv)
       IF (islv == 0) THEN
         psi1 = MIN (psi1, chi1)
       ELSE
         psi1 = zero
       END IF
     END IF
     
   ! *** CALCULATE SPECIATION ********************************************
     
     molal (2) = psi8 + psi7 + 2.d0*psi1                ! NAI
     molal (3) = psi4                                   ! NH4I
     molal (4) = psi6 + psi7                            ! CLI
     molal (5) = psi2 + psi1                            ! SO4I
     molal (6) = zero
     molal (7) = psi5 + psi8                            ! NO3I
     
     smin = 2.d0*molal(5)+molal(7)+molal(4)-molal(2)-molal(3) + organion
     CALL calcph (smin, hi, ohi)
     molal (1) = hi
     
   ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
     
     gnh3      = MAX(chi4 - psi4, tiny)
     ghno3     = MAX(chi5 - psi5, tiny)
     ghcl      = MAX(chi6 - psi6, tiny)
     
     cnh42s4   = zero
     cnh4no3   = zero
     cnacl     = MAX(chi7 - psi7, zero)
     cnano3    = MAX(chi8 - psi8, zero)
     cna2so4   = MAX(chi1 - psi1, zero)
     
   ! *** NH4Cl(s) calculations
     
     a3   = xk6 /(r*temp*r*temp)
     delt = MIN(gnh3, ghcl)
     bb = -(gnh3+ghcl)
     cc = gnh3*ghcl-a3
     dd = MAX(bb*bb-4.d0*cc, 0.d0)       ! US patch 12/20/01
     psi31 = 0.5D0*(-bb + SQRT(dd))
     psi32 = 0.5D0*(-bb - SQRT(dd))
     IF (delt-psi31 > zero .AND. psi31 > zero) THEN
       psi3 = psi31
     ELSE IF (delt-psi32 > zero .AND. psi32 > zero) THEN
       psi3 = psi32
     ELSE
       psi3 = zero
     END IF
     
   ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
     
     gnh3    = MAX(gnh3 - psi3, tiny)
     ghcl    = MAX(ghcl - psi3, tiny)
     cnh4cl  = psi3
     
     CALL calcmr                           ! Water content
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   
   ! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
   
   20    funch4a = molal(3)*molal(4)/ghcl/gnh3/a6/a4 - one
   
END FUNCTION funch4a
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCH3
!! *** CASE H3
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : NH4CL, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!======================================================================
SUBROUTINE calch3
   USE Isorropia_Module
   use solut
   implicit none
   
   double precision :: delta,frna,psi6lo,psi6hi,x1,y1,dx,x2,y2,x3,y3
   double precision, external :: funch3a
   integer :: i
   
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** REGIME DEPENDS ON THE EXISTANCE OF NITRATES ***********************
   
   IF (w(4) <= tiny) THEN        ! NO3 NOT EXIST, WATER NOT POSSIBLE
     scase = 'H3'
     CALL calch1a
     scase = 'H3'
     RETURN
   END IF
   
   ! *** SETUP PARAMETERS ************************************************
   
   calaou = .true.
   chi1   = w(2)                                ! CNA2SO4
   chi2   = zero                                ! CNH42S4
   chi3   = zero                                ! CNH4CL
   frna   = MAX (w(1)-2.d0*chi1, zero)
   chi8   = MIN (frna, w(4))                    ! CNANO3
   chi4   = w(3)                                ! NH3(g)
   chi5   = MAX (w(4)-chi8, zero)               ! HNO3(g)
   chi7   = MIN (MAX(frna-chi8, zero), w(5))    ! CNACL
   chi6   = MAX (w(5)-chi7, zero)               ! HCL(g)
   
   psi6lo = tiny
   psi6hi = chi6-tiny    ! MIN(CHI6-TINY, CHI4)
   
   ! *** INITIAL VALUES FOR BISECTION ************************************
   
   x1 = psi6lo
   y1 = funch3a (x1)
   IF (ABS(y1) <= eps .OR. chi6 <= tiny) GO TO 50
   
   ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
   
   dx = (psi6hi-psi6lo)/FLOAT(ndiv)
   DO  i=1,ndiv
     x2 = x1+dx
     y2 = funch3a (x2)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y2) < zero) GO TO 20  ! (Y1*Y2.LT.ZERO)
     x1 = x2
     y1 = y2
   END DO
   
   ! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
   
   IF (ABS(y2) > eps) y2 = funch3a (psi6lo)
   GO TO 50
   
   ! *** PERFORM BISECTION ***********************************************
   
   20    DO  i=1,maxit
     x3 = 0.5*(x1+x2)
     y3 = funch3a (x3)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y3) <= zero) THEN  ! (Y1*Y3 .LE. ZERO)
       y2    = y3
       x2    = x3
     ELSE
       y1    = y3
       x1    = x3
     END IF
     IF (ABS(x2-x1) <= eps*x1) GO TO 40
   END DO
   CALL pusherr (0002, 'CALCH3')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CONVERGED ; RETURN **********************************************
   
   40    x3 = 0.5*(x1+x2)
   y3 = funch3a (x3)
   
   ! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
   
   50    CONTINUE
   IF (molal(1) > tiny .AND. molal(5) > tiny) THEN
     CALL calchs4 (molal(1), molal(5), zero, delta)
     molal(1) = molal(1) - delta                     ! H+   EFFECT
     molal(5) = molal(5) - delta                     ! SO4  EFFECT
     molal(6) = delta                                ! HSO4 EFFECT
   END IF
   
END SUBROUTINE  calch3
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE FUNCH3A
!! *** CASE H3
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
DOUBLE PRECISION FUNCTION funch3a (x)
   USE Isorropia_Module
   use solut
   implicit none
   
   double precision, INTENT(IN)                         :: x
   double precision :: smin,hi,ohi,aa,bb,cc,a9,dd,diak,delt,psi31,psi32
   integer :: i,islv
   
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   psi6   = x
   psi1   = chi1
   psi2   = zero
   psi3   = zero
   psi7   = chi7
   psi8   = chi8
   frst   = .true.
   calain = .true.
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     
     a1  = xk5 *(water/gama(2))**3.0
     a4  = (xk2/xkw)*r*temp*(gama(10)/gama(5))**2.0
     a5  = xk4 *r*temp*(water/gama(10))**2.0
     a6  = xk3 *r*temp*(water/gama(11))**2.0
     a7  = xk8 *(water/gama(1))**2.0
     a8  = xk9 *(water/gama(3))**2.0
     a9  = xk1*water/gama(7)*(gama(8)/gama(7))**2.
     
   !  CALCULATE DISSOCIATION QUANTITIES
     
     psi5 = chi5*(psi6+psi7) - a6/a5*psi8*(chi6-psi6-psi3)
     psi5 = psi5/(a6/a5*(chi6-psi6-psi3) + psi6 + psi7)
     psi5 = MAX(psi5, tiny)
     
     IF (chi1 > tiny .AND. water > tiny) THEN  ! First try 3rd order soln
       bb   =-(chi4 + psi6 + psi5 + 1.d0/a4)
       cc   = chi4*(psi5+psi6)
       dd   = MAX(bb*bb-4.d0*cc, 0.d0)       ! US patch 12/20/01
       psi4 =0.5D0*(-bb - SQRT(dd))
       psi4 = MIN(psi4,chi4)
     ELSE
       psi4 = tiny
     END IF
     
     IF (chi7 > tiny .AND. water > tiny) THEN     ! NACL DISSOLUTION
       diak = (psi8-psi6)**2.d0 + 4.d0*a7
       psi7 = 0.5D0*( -(psi8+psi6) + SQRT(diak) )
       psi7 = MAX(MIN(psi7, chi7), zero)
     END IF
     
     IF (chi1 > tiny .AND. water > tiny) THEN     ! NA2SO4 DISSOLUTION
       aa = psi7+psi8
       bb = aa*aa
       cc =-a1/4.d0
       CALL poly3 (aa, bb, cc, psi1, islv)
       IF (islv == 0) THEN
         psi1 = MIN (psi1, chi1)
       ELSE
         psi1 = zero
       END IF
     END IF
     
   ! *** CALCULATE SPECIATION ********************************************
     
     molal (2) = psi8 + psi7 + 2.d0*psi1             ! NAI
     molal (3) = psi4                                ! NH4I
     molal (4) = psi6 + psi7                         ! CLI
     molal (5) = psi2 + psi1                         ! SO4I
     molal (6) = zero
     molal (7) = psi5 + psi8                         ! NO3I
     
     smin = 2.d0*molal(5)+molal(7)+molal(4)-molal(2)-molal(3) + organion
     CALL calcph (smin, hi, ohi)
     molal (1) = hi
     
   ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
     
     gnh3      = MAX(chi4 - psi4, tiny)
     ghno3     = MAX(chi5 - psi5, tiny)
     ghcl      = MAX(chi6 - psi6, tiny)
     
     cnh42s4   = zero
     cnh4no3   = zero
     cnacl     = MAX(chi7 - psi7, zero)
     cnano3    = MAX(chi8 - psi8, zero)
     cna2so4   = MAX(chi1 - psi1, zero)
     
   ! *** NH4Cl(s) calculations
     
     a3   = xk6 /(r*temp*r*temp)
     delt = MIN(gnh3, ghcl)
     bb = -(gnh3+ghcl)
     cc = gnh3*ghcl-a3
     dd = MAX(bb*bb-4.d0*cc, 0.d0)       ! US patch 12/20/01
     psi31 = 0.5D0*(-bb + SQRT(dd))
     psi32 = 0.5D0*(-bb - SQRT(dd))
     IF (delt-psi31 > zero .AND. psi31 > zero) THEN
       psi3 = psi31
     ELSE IF (delt-psi32 > zero .AND. psi32 > zero) THEN
       psi3 = psi32
     ELSE
       psi3 = zero
     END IF
     
   ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
     
     gnh3    = MAX(gnh3 - psi3, tiny)
     ghcl    = MAX(ghcl - psi3, tiny)
     cnh4cl  = psi3
     
     CALL calcmr                                 ! Water content
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   
   ! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
   
   20    funch3a = molal(3)*molal(4)/ghcl/gnh3/a6/a4 - one
   
END FUNCTION funch3a
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCH2
!! *** CASE H2
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!!     2. SOLID & LIQUID AEROSOL POSSIBLE
!!     3. SOLIDS POSSIBLE : NH4Cl, NA2SO4, NANO3, NACL
!!     THERE ARE THREE REGIMES IN THIS CASE:
!!     1. NH4NO3(s) POSSIBLE. LIQUID & SOLID AEROSOL (SUBROUTINE CALCH2A)
!!     2. NH4NO3(s) NOT POSSIBLE, AND RH < MDRH. SOLID AEROSOL ONLY
!!     3. NH4NO3(s) NOT POSSIBLE, AND RH >= MDRH. (MDRH REGION)
!!     REGIMES 2. AND 3. ARE CONSIDERED TO BE THE SAME AS CASES H1A, H2B
!!     RESPECTIVELY (BECAUSE MDRH POINTS COINCIDE).
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calch2
   USE Isorropia_Module
   implicit none
   
   EXTERNAL calch1a, calch3
   
   ! *** REGIME DEPENDS ON THE EXISTANCE OF NITRATES ***********************
   
   IF (w(4) > tiny) THEN        ! NO3 EXISTS, WATER POSSIBLE
     scase = 'H2 ; SUBCASE 1'
     CALL calch2a
     scase = 'H2 ; SUBCASE 1'
   ELSE                          ! NO3 NON EXISTANT, WATER NOT POSSIBLE
     scase = 'H2 ; SUBCASE 1'
     CALL calch1a
     scase = 'H2 ; SUBCASE 1'
   END IF
   
   IF (water <= tiny .AND. rh < drmh2) THEN      ! DRY AEROSOL
     scase = 'H2 ; SUBCASE 2'
     
   ELSE IF (water <= tiny .AND. rh >= drmh2) THEN  ! MDRH OF H2
     scase = 'H2 ; SUBCASE 3'
     CALL calcmdrh (rh, drmh2, drnano3, calch1a, calch3)
     scase = 'H2 ; SUBCASE 3'
   END IF
   
END SUBROUTINE  calch2
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCH2A
!! *** CASE H2 ; SUBCASE 1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : NH4CL, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calch2a
   USE Isorropia_Module
   use solut
   implicit none
   
   double precision :: delta,frna,psi6lo,psi6hi,x1,y1,dx,x2,y2,x3,y3
   double precision, external :: funch2a
   integer :: i
   
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   calaou = .true.
   chi1   = w(2)                                ! CNA2SO4
   chi2   = zero                                ! CNH42S4
   chi3   = zero                                ! CNH4CL
   frna   = MAX (w(1)-2.d0*chi1, zero)
   chi8   = MIN (frna, w(4))                    ! CNANO3
   chi4   = w(3)                                ! NH3(g)
   chi5   = MAX (w(4)-chi8, zero)               ! HNO3(g)
   chi7   = MIN (MAX(frna-chi8, zero), w(5))    ! CNACL
   chi6   = MAX (w(5)-chi7, zero)               ! HCL(g)
   
   psi6lo = tiny
   psi6hi = chi6-tiny    ! MIN(CHI6-TINY, CHI4)
   
   ! *** INITIAL VALUES FOR BISECTION ************************************
   
   x1 = psi6lo
   y1 = funch2a (x1)
   IF (ABS(y1) <= eps .OR. chi6 <= tiny) GO TO 50
   
   ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
   
   dx = (psi6hi-psi6lo)/FLOAT(ndiv)
   DO  i=1,ndiv
     x2 = x1+dx
     y2 = funch2a (x2)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y2) < zero) GO TO 20  ! (Y1*Y2.LT.ZERO)
     x1 = x2
     y1 = y2
   END DO
   
   ! *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
   
   IF (y2 > eps) y2 = funch2a (psi6lo)
   GO TO 50
   
   ! *** PERFORM BISECTION ***********************************************
   
   20    DO  i=1,maxit
     x3 = 0.5*(x1+x2)
     y3 = funch2a (x3)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y3) <= zero) THEN  ! (Y1*Y3 .LE. ZERO)
       y2    = y3
       x2    = x3
     ELSE
       y1    = y3
       x1    = x3
     END IF
     IF (ABS(x2-x1) <= eps*x1) GO TO 40
   END DO
   CALL pusherr (0002, 'CALCH2A')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CONVERGED ; RETURN **********************************************
   
   40    x3 = 0.5*(x1+x2)
   y3 = funch2a (x3)
   
   ! *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
   
   50    CONTINUE
   IF (molal(1) > tiny .AND. molal(5) > tiny) THEN
     CALL calchs4 (molal(1), molal(5), zero, delta)
     molal(1) = molal(1) - delta                    ! H+   EFFECT
     molal(5) = molal(5) - delta                    ! SO4  EFFECT
     molal(6) = delta                               ! HSO4 EFFECT
   END IF
   
END SUBROUTINE  calch2a
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE FUNCH2A
!! *** CASE H2 ; SUBCASE 1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
DOUBLE PRECISION FUNCTION funch2a (x)
   USE Isorropia_Module
   use solut
   implicit none
   
   double precision, INTENT(IN)                         :: x
   double precision :: smin,hi,ohi,aa,bb,cc,a64,a9,dd,diak,delt,psi31,psi32
   integer :: i,islv
   
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   psi6   = x
   psi1   = chi1
   psi2   = zero
   psi3   = zero
   psi7   = chi7
   psi8   = chi8
   frst   = .true.
   calain = .true.
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     
     a1  = xk5 *(water/gama(2))**3.0
     a4  = (xk2/xkw)*r*temp*(gama(10)/gama(5))**2.0
     a5  = xk4 *r*temp*(water/gama(10))**2.0
     a6  = xk3 *r*temp*(water/gama(11))**2.0
     a7  = xk8 *(water/gama(1))**2.0
     a8  = xk9 *(water/gama(3))**2.0
     a64 = (xk3*xk2/xkw)*(gama(10)/gama(5)/gama(11))**2.0
     a64 = a64*(r*temp*water)**2.0
     a9  = xk1*water/gama(7)*(gama(8)/gama(7))**2.
     
   !  CALCULATE DISSOCIATION QUANTITIES
     
     psi5 = chi5*(psi6+psi7) - a6/a5*psi8*(chi6-psi6-psi3)
     psi5 = psi5/(a6/a5*(chi6-psi6-psi3) + psi6 + psi7)
     psi5 = MAX(psi5, tiny)
     
     IF (chi1 > tiny .AND. water > tiny) THEN  ! First try 3rd order soln
       bb   =-(chi4 + psi6 + psi5 + 1.d0/a4)
       cc   = chi4*(psi5+psi6)
       dd   = MAX(bb*bb-4.d0*cc, 0.d0)       ! US patch 12/20/01
       psi4 =0.5D0*(-bb - SQRT(dd))
       psi4 = MIN(psi4,chi4)
     ELSE
       psi4 = tiny
     END IF
     
     IF (chi7 > tiny .AND. water > tiny) THEN     ! NACL DISSOLUTION
       diak = (psi8-psi6)**2.d0 + 4.d0*a7
       psi7 = 0.5D0*( -(psi8+psi6) + SQRT(diak) )
       psi7 = MAX(MIN(psi7, chi7), zero)
     END IF
     
     IF (chi8 > tiny .AND. water > tiny) THEN     ! NANO3 DISSOLUTION
       diak = (psi7-psi5)**2.d0 + 4.d0*a8
       psi8 = 0.5D0*( -(psi7+psi5) + SQRT(diak) )
       psi8 = MAX(MIN(psi8, chi8), zero)
     END IF
     
     IF (chi1 > tiny .AND. water > tiny) THEN     ! NA2SO4 DISSOLUTION
       aa = psi7+psi8
       bb = aa*aa
       cc =-a1/4.d0
       CALL poly3 (aa, bb, cc, psi1, islv)
       IF (islv == 0) THEN
         psi1 = MIN (psi1, chi1)
       ELSE
         psi1 = zero
       END IF
     END IF
     
   ! *** CALCULATE SPECIATION ********************************************
     
     molal (2) = psi8 + psi7 + 2.d0*psi1                 ! NAI
     molal (3) = psi4                                    ! NH4I
     molal (4) = psi6 + psi7                             ! CLI
     molal (5) = psi2 + psi1                             ! SO4I
     molal (6) = zero                                    ! HSO4I
     molal (7) = psi5 + psi8                             ! NO3I
     
     smin = 2.d0*molal(5)+molal(7)+molal(4)-molal(2)-molal(3) +organion
     CALL calcph (smin, hi, ohi)
     molal (1) = hi
     
     gnh3      = MAX(chi4 - psi4, tiny)
     ghno3     = MAX(chi5 - psi5, tiny)
     ghcl      = MAX(chi6 - psi6, tiny)
     
     cnh42s4   = zero
     cnh4no3   = zero
     cnacl     = MAX(chi7 - psi7, zero)
     cnano3    = MAX(chi8 - psi8, zero)
     cna2so4   = MAX(chi1 - psi1, zero)
     
   ! *** NH4Cl(s) calculations
     
     a3   = xk6 /(r*temp*r*temp)
     delt = MIN(gnh3, ghcl)
     bb = -(gnh3+ghcl)
     cc = gnh3*ghcl-a3
     dd = MAX(bb*bb-4.d0*cc, 0.d0)       ! US patch 12/20/01
     psi31 = 0.5D0*(-bb + SQRT(dd))
     psi32 = 0.5D0*(-bb - SQRT(dd))
     IF (delt-psi31 > zero .AND. psi31 > zero) THEN
       psi3 = psi31
     ELSE IF (delt-psi32 > zero .AND. psi32 > zero) THEN
       psi3 = psi32
     ELSE
       psi3 = zero
     END IF
     
   ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
     
     gnh3    = MAX(gnh3 - psi3, tiny)
     ghcl    = MAX(ghcl - psi3, tiny)
     cnh4cl  = psi3
     
     CALL calcmr                        ! Water content
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   
   ! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
   
   20    funch2a = molal(3)*molal(4)/ghcl/gnh3/a64 - one
   
END FUNCTION funch2a
   
!=======================================================================
! *** ISORROPIA CODE
! *** SUBROUTINE CALCH1
! *** CASE H1
!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!     2. SOLID AEROSOL ONLY
!     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NA2SO4
!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
!     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCH1A)
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! @author Athanasios Nenes
!=======================================================================
SUBROUTINE calch1
   USE Isorropia_Module
   implicit none

   EXTERNAL calch1a, calch2a
   
   ! *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************
   
   IF (rh < drmh1) THEN
     scase = 'H1 ; SUBCASE 1'
     CALL calch1a              ! SOLID PHASE ONLY POSSIBLE
     scase = 'H1 ; SUBCASE 1'
   ELSE
     scase = 'H1 ; SUBCASE 2'  ! LIQUID & SOLID PHASE POSSIBLE
     CALL calcmdrh (rh, drmh1, drnh4no3, calch1a, calch2a)
     scase = 'H1 ; SUBCASE 2'
   END IF
   
END SUBROUTINE  calch1
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCH1A
!! *** CASE H1 ; SUBCASE 1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!!     2. SOLID AEROSOL ONLY
!!     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NANO3, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calch1a
   USE Isorropia_Module
   implicit none

   DOUBLE PRECISION :: lamda, lamda1, lamda2, kapa, kapa1, kapa2, nafr, no3fr
   DOUBLE PRECISION :: clfr,alf,bet,gam,rtsq,a1,a2,theta1,theta2,bb,cc,dd
   double precision :: sqdd,dd1,dd2,sqdd1,sqdd2
   
   ! *** CALCULATE NON VOLATILE SOLIDS ***********************************
   
   cna2so4 = w(2)
   cnh42s4 = zero
   nafr    = MAX (w(1)-2*cna2so4, zero)
   cnano3  = MIN (nafr, w(4))
   no3fr   = MAX (w(4)-cnano3, zero)
   cnacl   = MIN (MAX(nafr-cnano3, zero), w(5))
   clfr    = MAX (w(5)-cnacl, zero)
   
   ! *** CALCULATE VOLATILE SPECIES **************************************
   
   alf     = w(3)                     ! FREE NH3
   bet     = clfr                     ! FREE CL
   gam     = no3fr                    ! FREE NO3
   
   rtsq    = r*temp*r*temp
   a1      = xk6/rtsq
   a2      = xk10/rtsq
   
   theta1  = gam - bet*(a2/a1)
   theta2  = a2/a1
   
   ! QUADRATIC EQUATION SOLUTION
   
   bb      = (theta1-alf-bet*(one+theta2))/(one+theta2)
   cc      = (alf*bet-a1-bet*theta1)/(one+theta2)
   dd      = MAX( bb*bb - 4.0D0*cc, 0.d0 )  ! US patch 12/20/01
   IF (dd < zero) GO TO 100   ! Solve each reaction seperately
   
   ! TWO ROOTS FOR KAPA, CHECK AND SEE IF ANY VALID
   
   sqdd    = SQRT(dd)
   kapa1   = 0.5D0*(-bb+sqdd)
   kapa2   = 0.5D0*(-bb-sqdd)
   lamda1  = theta1 + theta2*kapa1
   lamda2  = theta1 + theta2*kapa2
   
   IF (kapa1 >= zero .AND. lamda1 >= zero) THEN
     IF (alf-kapa1-lamda1 >= zero .AND.  &
           bet-kapa1 >= zero .AND. gam-lamda1 >= zero) THEN
       kapa = kapa1
       lamda= lamda1
       GO TO 200
     END IF
   END IF
   
   IF (kapa2 >= zero .AND. lamda2 >= zero) THEN
     IF (alf-kapa2-lamda2 >= zero .AND.  &
           bet-kapa2 >= zero .AND. gam-lamda2 >= zero) THEN
       kapa = kapa2
       lamda= lamda2
       GO TO 200
     END IF
   END IF
   
   ! SEPERATE SOLUTION OF NH4CL & NH4NO3 EQUILIBRIA
   
   100   kapa  = zero
   lamda = zero
   dd1   = (alf+bet)*(alf+bet) - 4.0D0*(alf*bet-a1)
   dd2   = (alf+gam)*(alf+gam) - 4.0D0*(alf*gam-a2)
   
   ! NH4CL EQUILIBRIUM
   
   IF (dd1 >= zero) THEN
     sqdd1 = SQRT(dd1)
     kapa1 = 0.5D0*(alf+bet + sqdd1)
     kapa2 = 0.5D0*(alf+bet - sqdd1)
     
     IF (kapa1 >= zero .AND. kapa1 <= MIN(alf,bet)) THEN
       kapa = kapa1
     ELSE IF (kapa2 >= zero .AND. kapa2 <= MIN(alf,bet)) THEN
       kapa = kapa2
     ELSE
       kapa = zero
     END IF
   END IF
   
   ! NH4NO3 EQUILIBRIUM
   
   IF (dd2 >= zero) THEN
     sqdd2 = SQRT(dd2)
     lamda1= 0.5D0*(alf+gam + sqdd2)
     lamda2= 0.5D0*(alf+gam - sqdd2)
     
     IF (lamda1 >= zero .AND. lamda1 <= MIN(alf,gam)) THEN
       lamda = lamda1
     ELSE IF (lamda2 >= zero .AND. lamda2 <= MIN(alf,gam)) THEN
       lamda = lamda2
     ELSE
       lamda = zero
     END IF
   END IF
   
   ! IF BOTH KAPA, LAMDA ARE > 0, THEN APPLY EXISTANCE CRITERION
   
   IF (kapa > zero .AND. lamda > zero) THEN
     IF (bet < lamda/theta1) THEN
       kapa = zero
     ELSE
       lamda= zero
     END IF
   END IF
   
   ! *** CALCULATE COMPOSITION OF VOLATILE SPECIES ***********************
   
   200   CONTINUE
   cnh4no3 = lamda
   cnh4cl  = kapa
   
   gnh3    = alf - kapa - lamda
   ghno3   = gam - lamda
   ghcl    = bet - kapa
   
END SUBROUTINE  calch1a

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCI6
!! *** CASE I6
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!!     2. SOLID & LIQUID AEROSOL POSSIBLE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calci6
   USE Isorropia_Module
   use solut
   implicit none
   
   DOUBLE PRECISION :: bb,cc,dd
   integer :: i
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** FIND DRY COMPOSITION **********************************************
   
   CALL calci1a
   
   ! *** SETUP PARAMETERS ************************************************
   
   chi1 = cnh4hs4               ! Save from CALCI1 run
   chi2 = clc
   chi3 = cnahso4
   chi4 = cna2so4
   chi5 = cnh42s4
   
   psi1 = cnh4hs4               ! ASSIGN INITIAL PSI's
   psi2 = clc
   psi3 = cnahso4
   psi4 = cna2so4
   psi5 = cnh42s4
   
   calaou = .true.              ! Outer loop activity calculation flag
   frst   = .true.
   calain = .true.
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     
     a6 = xk1 *water/gama(7)*(gama(8)/gama(7))**2.
     
   !  CALCULATE DISSOCIATION QUANTITIES
     
     bb   = psi2 + psi4 + psi5 + a6 +organion        ! PSI6
     cc   =-a6*(psi2+psi3+psi1) +organion*(psi2+psi4+psi5)
     dd   = MAX(bb*bb-4.d0*cc, 0.d0)       ! US patch 12/20/01
     psi6 = 0.5D0*(-bb + SQRT(dd))
     
   ! *** CALCULATE SPECIATION ********************************************
     
     molal (1) = psi6 + organion                         ! HI
     molal (2) = 2.d0*psi4 + psi3                        ! NAI
     molal (3) = 3.d0*psi2 + 2.d0*psi5 + psi1            ! NH4I
     molal (5) = psi2 + psi4 + psi5 + psi6               ! SO4I
     molal (6) = psi2 + psi3 + psi1 - psi6               ! HSO4I
     clc       = zero
     cnahso4   = zero
     cna2so4   = chi4 - psi4
     cnh42s4   = zero
     cnh4hs4   = zero
     CALL calcmr                                         ! Water content
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   
   20    RETURN

END SUBROUTINE  calci6
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCI5
!! *** CASE I5
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!!     2. SOLID & LIQUID AEROSOL POSSIBLE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!=======================================================================
SUBROUTINE calci5
   USE Isorropia_Module
   use solut
   implicit none
   
   double precision :: psi4lo,psi4hi,y1,x1,yhi,dx,x2,y2,x3,y3,ylo
   double precision, external :: funci5a
   integer :: i
   
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** FIND DRY COMPOSITION **********************************************
   
   CALL calci1a
   
   ! *** SETUP PARAMETERS ************************************************
   
   chi1 = cnh4hs4               ! Save from CALCI1 run
   chi2 = clc
   chi3 = cnahso4
   chi4 = cna2so4
   chi5 = cnh42s4
   
   psi1 = cnh4hs4               ! ASSIGN INITIAL PSI's
   psi2 = clc
   psi3 = cnahso4
   psi4 = zero
   psi5 = cnh42s4
   
   calaou =.true.               ! Outer loop activity calculation flag
   psi4lo = zero                ! Low  limit
   psi4hi = chi4                ! High limit
   
   ! *** IF NA2SO4(S) =0, CALL FUNCI5B FOR Y4=0 ***************************
   
   IF (chi4 <= tiny) THEN
     y1 = funci5a (zero)
     GO TO 50
   END IF
   
   ! *** INITIAL VALUES FOR BISECTION ************************************
   
   x1 = psi4hi
   y1 = funci5a (x1)
   yhi= y1                      ! Save Y-value at HI position
   
   ! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NA2SO4 **
   
   IF (ABS(y1) <= eps .OR. yhi < zero) GO TO 50
   
   ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
   
   dx = (psi4hi-psi4lo)/FLOAT(ndiv)
   DO  i=1,ndiv
     x2 = x1-dx
     y2 = funci5a (x2)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y2) < zero) GO TO 20  ! (Y1*Y2.LT.ZERO)
     x1 = x2
     y1 = y2
   END DO
   
   ! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH4CL
   
   ylo= y1                      ! Save Y-value at Hi position
   IF (ylo > zero .AND. yhi > zero) THEN
     y3 = funci5a (zero)
     GO TO 50
   ELSE IF (ABS(y2) < eps) THEN   ! x2 IS A SOLUTION
     GO TO 50
   ELSE
     CALL pusherr (0001, 'CALCI5')    ! WARNING ERROR: NO SOLUTION
     GO TO 50
   END IF
   
   ! *** PERFORM BISECTION ***********************************************
   
   20    DO  i=1,maxit
     x3 = 0.5*(x1+x2)
     y3 = funci5a (x3)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y3) <= zero) THEN  ! (Y1*Y3 .LE. ZERO)
       y2    = y3
       x2    = x3
     ELSE
       y1    = y3
       x1    = x3
     END IF
     IF (ABS(x2-x1) <= eps*x1) GO TO 40
   END DO
   CALL pusherr (0002, 'CALCI5')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CONVERGED ; RETURN **********************************************
   
   40    x3 = 0.5*(x1+x2)
   y3 = funci5a (x3)
   
   50    RETURN
   
END SUBROUTINE  calci5
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE FUNCI5A
!! *** CASE I5
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!!     2. SOLID & LIQUID AEROSOL POSSIBLE
!!     3. SOLIDS POSSIBLE : NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
DOUBLE PRECISION FUNCTION funci5a (p4)
   USE Isorropia_Module
   use solut
   implicit none
   
   DOUBLE PRECISION, INTENT(IN)                         :: p4
   double precision :: bb,cc,dd
   integer :: i
   
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   psi4   = p4     ! PSI3 already assigned in FUNCI5A
   frst   = .true.
   calain = .true.
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     
     a4 = xk5 *(water/gama(2))**3.0
     a5 = xk7 *(water/gama(4))**3.0
     a6 = xk1 *water/gama(7)*(gama(8)/gama(7))**2.
     
   !  CALCULATE DISSOCIATION QUANTITIES
     
     bb   = psi2 + psi4 + psi5 + a6 +organion       ! PSI6
     cc   =-a6*(psi2 + psi3 + psi1)+organion*(psi2+psi4+psi5)
     dd   = MAX(bb*bb-4.d0*cc, 0.d0)       ! US patch 12/20/01
     psi6 = 0.5D0*(-bb + SQRT(dd))
     
   ! *** CALCULATE SPECIATION ********************************************
     
     molal (1) = psi6+ organion                  ! HI
     molal (2) = 2.d0*psi4 + psi3                ! NAI
     molal (3) = 3.d0*psi2 + 2.d0*psi5 + psi1    ! NH4I
     molal (5) = psi2 + psi4 + psi5 + psi6       ! SO4I
     molal (6) = psi2 + psi3 + psi1 - psi6       ! HSO4I
     clc       = zero
     cnahso4   = zero
     cna2so4   = chi4 - psi4
     cnh42s4   = zero
     cnh4hs4   = zero
     CALL calcmr                                 ! Water content
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   
   ! *** CALCULATE OBJECTIVE FUNCTION ************************************
   
   20    a4     = xk5 *(water/gama(2))**3.0
   funci5a= molal(5)*molal(2)*molal(2)/a4 - one

END FUNCTION funci5a

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCI4
!! *** CASE I4
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!!     2. SOLID & LIQUID AEROSOL POSSIBLE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calci4
   USE Isorropia_Module
   use solut
   implicit none
   
   double precision :: psi4lo,psi4hi,x1,y1,yhi,dx,x2,y2,ylo,x3,y3
   double precision, external :: funci4a
   integer :: i
   
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** FIND DRY COMPOSITION **********************************************
   
   CALL calci1a
   
   ! *** SETUP PARAMETERS ************************************************
   
   chi1 = cnh4hs4               ! Save from CALCI1 run
   chi2 = clc
   chi3 = cnahso4
   chi4 = cna2so4
   chi5 = cnh42s4
   
   psi1 = cnh4hs4               ! ASSIGN INITIAL PSI's
   psi2 = clc
   psi3 = cnahso4
   psi4 = zero
   psi5 = zero
   
   calaou = .true.              ! Outer loop activity calculation flag
   psi4lo = zero                ! Low  limit
   psi4hi = chi4                ! High limit
   
   ! *** IF NA2SO4(S) =0, CALL FUNCI4B FOR Y4=0 ***************************
   
   IF (chi4 <= tiny) THEN
     y1 = funci4a (zero)
     GO TO 50
   END IF
   
   ! *** INITIAL VALUES FOR BISECTION ************************************
   
   x1 = psi4hi
   y1 = funci4a (x1)
   yhi= y1                      ! Save Y-value at HI position
   
   ! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NA2SO4 **
   
   IF (ABS(y1) <= eps .OR. yhi < zero) GO TO 50
   
   ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
   
   dx = (psi4hi-psi4lo)/FLOAT(ndiv)
   DO  i=1,ndiv
     x2 = x1-dx
     y2 = funci4a (x2)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y2) < zero) GO TO 20  ! (Y1*Y2.LT.ZERO)
     x1 = x2
     y1 = y2
   END DO
   
   ! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH4CL
   
   ylo= y1                      ! Save Y-value at Hi position
   IF (ylo > zero .AND. yhi > zero) THEN
     y3 = funci4a (zero)
     GO TO 50
   ELSE IF (ABS(y2) < eps) THEN   ! x2 IS A SOLUTION
     GO TO 50
   ELSE
     CALL pusherr (0001, 'CALCI4')    ! WARNING ERROR: NO SOLUTION
     GO TO 50
   END IF
   
   ! *** PERFORM BISECTION ***********************************************
   
   20    DO  i=1,maxit
     x3 = 0.5*(x1+x2)
     y3 = funci4a (x3)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y3) <= zero) THEN  ! (Y1*Y3 .LE. ZERO)
       y2    = y3
       x2    = x3
     ELSE
       y1    = y3
       x1    = x3
     END IF
     IF (ABS(x2-x1) <= eps*x1) GO TO 40
   END DO
   CALL pusherr (0002, 'CALCI4')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CONVERGED ; RETURN **********************************************
   
   40    x3 = 0.5*(x1+x2)
   y3 = funci4a (x3)
   
   50    RETURN
   
END SUBROUTINE  calci4
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE FUNCI4A
!! *** CASE I4
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!!     2. SOLID & LIQUID AEROSOL POSSIBLE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
DOUBLE PRECISION FUNCTION funci4a (p4)
   USE Isorropia_Module
   use solut
   implicit none
   
   double precision, INTENT(IN)                         :: p4
   double precision :: bb,cc,dd
   integer :: i
   
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   psi4   = p4     ! PSI3 already assigned in FUNCI4A
   frst   = .true.
   calain = .true.
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     
     a4 = xk5 *(water/gama(2))**3.0
     a5 = xk7 *(water/gama(4))**3.0
     a6 = xk1 *water/gama(7)*(gama(8)/gama(7))**2.
     a7 = SQRT(a4/a5)
     
   !  CALCULATE DISSOCIATION QUANTITIES
     
     bb   = psi2 + psi4 + psi5 + a6 + organion         ! PSI6
     cc   =-a6*(psi2 + psi3 + psi1)+organion*(psi2+psi4+psi5)
     dd   = MAX(bb*bb-4.d0*cc, 0.d0)       ! US patch 12/20/01
     psi6 = 0.5D0*(-bb + SQRT(dd))
     
     psi5 = (psi3 + 2.d0*psi4 - a7*(3.d0*psi2 + psi1))/2.d0/a7
     psi5 = MIN (psi5, chi5)
     
   ! *** CALCULATE SPECIATION ********************************************
     
     molal (1) = psi6 + organion                 ! HI
     molal (2) = 2.d0*psi4 + psi3                ! NAI
     molal (3) = 3.d0*psi2 + 2.d0*psi5 + psi1    ! NH4I
     molal (5) = psi2 + psi4 + psi5 + psi6       ! SO4I
     molal (6) = psi2 + psi3 + psi1 - psi6       ! HSO4I
     clc       = zero
     cnahso4   = zero
     cna2so4   = chi4 - psi4
     cnh42s4   = chi5 - psi5
     cnh4hs4   = zero
     CALL calcmr                                 ! Water content
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   
   ! *** CALCULATE OBJECTIVE FUNCTION ************************************
   
   20    a4     = xk5 *(water/gama(2))**3.0
   funci4a= molal(5)*molal(2)*molal(2)/a4 - one

END FUNCTION funci4a

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCI3
!! *** CASE I3
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!!     2. SOLID & LIQUID AEROSOL POSSIBLE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, NAHSO4, LC
!!     THERE ARE THREE REGIMES IN THIS CASE:
!!     1.(NA,NH4)HSO4(s) POSSIBLE. LIQUID & SOLID AEROSOL (SUBROUTINE CALCI3A)
!!     2.(NA,NH4)HSO4(s) NOT POSSIBLE, AND RH < MDRH. SOLID AEROSOL ONLY
!!     3.(NA,NH4)HSO4(s) NOT POSSIBLE, AND RH >= MDRH. SOLID & LIQUID AEROSOL
!!     REGIMES 2. AND 3. ARE CONSIDERED TO BE THE SAME AS CASES I1A, I2B
!!     RESPECTIVELY
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calci3
   USE Isorropia_Module
   implicit none
   
   EXTERNAL calci1a, calci4
   Integer :: i
   
   ! *** FIND DRY COMPOSITION **********************************************
   
   CALL calci1a
   
   ! *** REGIME DEPENDS UPON THE POSSIBLE SOLIDS & RH **********************
   
   IF (cnh4hs4 > tiny .OR. cnahso4 > tiny) THEN
     scase = 'I3 ; SUBCASE 1'
     CALL calci3a                     ! FULL SOLUTION
     scase = 'I3 ; SUBCASE 1'
   END IF
   
   IF (water <= tiny) THEN
     IF (rh < drmi3) THEN         ! SOLID SOLUTION
       water = tiny
       DO  i=1,nions
         molal(i) = zero
       END DO
       CALL calci1a
       scase = 'I3 ; SUBCASE 2'
       
     ELSE IF (rh >= drmi3) THEN     ! MDRH OF I3
       scase = 'I3 ; SUBCASE 3'
       CALL calcmdrh (rh, drmi3, drlc, calci1a, calci4)
       scase = 'I3 ; SUBCASE 3'
     END IF
   END IF
   
END SUBROUTINE  calci3
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCI3A
!! *** CASE I3 ; SUBCASE 1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!!     2. SOLID & LIQUID AEROSOL POSSIBLE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, LC
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calci3a
   USE Isorropia_Module
   use solut
   implicit none
   
   double precision :: psi2lo,psi2hi,x1,y1,yhi,dx,x2,y2,x3,y3
   double precision, external :: funci3a
   integer :: i
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** FIND DRY COMPOSITION **********************************************
   
   CALL calci1a         ! Needed when called from CALCMDRH
   
   ! *** SETUP PARAMETERS ************************************************
   
   chi1 = cnh4hs4               ! Save from CALCI1 run
   chi2 = clc
   chi3 = cnahso4
   chi4 = cna2so4
   chi5 = cnh42s4
   
   psi1 = cnh4hs4               ! ASSIGN INITIAL PSI's
   psi2 = zero
   psi3 = cnahso4
   psi4 = zero
   psi5 = zero
   
   calaou = .true.              ! Outer loop activity calculation flag
   psi2lo = zero                ! Low  limit
   psi2hi = chi2                ! High limit
   
   ! *** INITIAL VALUES FOR BISECTION ************************************
   
   x1 = psi2hi
   y1 = funci3a (x1)
   yhi= y1                      ! Save Y-value at HI position
   
   ! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH LC *********
   
   IF (yhi < eps) GO TO 50
   
   ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
   
   dx = (psi2hi-psi2lo)/FLOAT(ndiv)
   DO  i=1,ndiv
     x2 = MAX(x1-dx, psi2lo)
     y2 = funci3a (x2)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y2) < zero) GO TO 20  ! (Y1*Y2.LT.ZERO)
     x1 = x2
     y1 = y2
   END DO
   
   ! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH LC
   
   IF (y2 > eps) y2 = funci3a (zero)
   GO TO 50
   
   ! *** PERFORM BISECTION ***********************************************
   
   20    DO  i=1,maxit
     x3 = 0.5*(x1+x2)
     y3 = funci3a (x3)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y3) <= zero) THEN  ! (Y1*Y3 .LE. ZERO)
       y2    = y3
       x2    = x3
     ELSE
       y1    = y3
       x1    = x3
     END IF
     IF (ABS(x2-x1) <= eps*x1) GO TO 40
   END DO
   CALL pusherr (0002, 'CALCI3A')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CONVERGED ; RETURN **********************************************
   
   40    x3 = 0.5*(x1+x2)
   y3 = funci3a (x3)
   
   50    RETURN
   
END SUBROUTINE  calci3a
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE FUNCI3A
!! *** CASE I3 ; SUBCASE 1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!!     2. SOLID & LIQUID AEROSOL POSSIBLE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, LC
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
DOUBLE PRECISION FUNCTION funci3a (p2)
   USE Isorropia_Module
   use solut
   IMPLICIT NONE
   
   double precision , INTENT(IN)                         :: p2
   double precision :: psi4lo,psi4hi,x1,y1,yhi,dx,x2,y2,x3,y3
   double precision, external :: funci3b
   integer :: i
   
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   psi2   = p2                  ! Save PSI2 in COMMON BLOCK
   psi4lo = zero                ! Low  limit for PSI4
   psi4hi = chi4                ! High limit for PSI4
   
   ! *** IF NH3 =0, CALL FUNCI3B FOR Y4=0 ********************************
   
   IF (chi4 <= tiny) THEN
     funci3a = funci3b (zero)
     GO TO 50
   END IF
   
   ! *** INITIAL VALUES FOR BISECTION ************************************
   
   x1 = psi4hi
   y1 = funci3b (x1)
   IF (ABS(y1) <= eps) GO TO 50
   yhi= y1                      ! Save Y-value at HI position
   
   ! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NA2SO4 *****
   
   IF (yhi < zero) GO TO 50
   
   ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
   
   dx = (psi4hi-psi4lo)/FLOAT(ndiv)
   DO  i=1,ndiv
     x2 = MAX(x1-dx, psi4lo)
     y2 = funci3b (x2)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y2) < zero) GO TO 20  ! (Y1*Y2.LT.ZERO)
     x1 = x2
     y1 = y2
   END DO
   
   ! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NA2SO4
   
   IF (y2 > eps) y2 = funci3b (psi4lo)
   GO TO 50
   
   ! *** PERFORM BISECTION ***********************************************
   
   20    DO  i=1,maxit
     x3 = 0.5*(x1+x2)
     y3 = funci3b (x3)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y3) <= zero) THEN  ! (Y1*Y3 .LE. ZERO)
       y2    = y3
       x2    = x3
     ELSE
       y1    = y3
       x1    = x3
     END IF
     IF (ABS(x2-x1) <= eps*x1) GO TO 40
   END DO
   CALL pusherr (0004, 'FUNCI3A')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** INNER LOOP CONVERGED **********************************************
   
   40    x3 = 0.5*(x1+x2)
   y3 = funci3b (x3)
   
   ! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
   
   50    a2      = xk13*(water/gama(13))**5.0
   funci3a = molal(5)*molal(6)*molal(3)**3.d0/a2 - one

END FUNCTION funci3a
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** FUNCTION FUNCI3B
!! *** CASE I3 ; SUBCASE 2
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!!     2. SOLID & LIQUID AEROSOL POSSIBLE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, LC
!!     SOLUTION IS SAVED IN COMMON BLOCK /CASE/
!!=======================================================================
DOUBLE PRECISION FUNCTION funci3b (p4)
   USE Isorropia_Module
   use solut
   IMPLICIT NONE
   
   double precision , INTENT(IN)                         :: p4
   double precision :: bb,cc,dd
   integer :: i
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   psi4   = p4
   
   ! *** SETUP PARAMETERS ************************************************
   
   frst   = .true.
   calain = .true.
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     
     a4 = xk5*(water/gama(2))**3.0
     a5 = xk7*(water/gama(4))**3.0
     a6 = xk1*water/gama(7)*(gama(8)/gama(7))**2.
     a7 = SQRT(a4/a5)
     
   !  CALCULATE DISSOCIATION QUANTITIES
     
     bb   = psi2 + psi4 + psi5 + a6 + organion       ! PSI6
     cc   =-a6*(psi2 + psi3 + psi1) +organion*(psi2+psi4+psi5)
     dd   = MAX(bb*bb-4.d0*cc, 0.d0)       ! US patch 12/20/01
     psi6 = 0.5D0*(-bb + SQRT(dd))
     
     psi5 = (psi3 + 2.d0*psi4 - a7*(3.d0*psi2 + psi1))/2.d0/a7
     psi5 = MIN (psi5, chi5)
     
   ! *** CALCULATE SPECIATION ********************************************
     
     molal(1) = psi6 +organion                        ! HI
     molal(2) = 2.d0*psi4 + psi3                      ! NAI
     molal(3) = 3.d0*psi2 + 2.d0*psi5 + psi1          ! NH4I
     molal(5) = psi2 + psi4 + psi5 + psi6             ! SO4I
     molal(6) = MAX(psi2 + psi3 + psi1 - psi6, tiny)  ! HSO4I
     clc      = MAX(chi2 - psi2, zero)
     cnahso4  = zero
     cna2so4  = MAX(chi4 - psi4, zero)
     cnh42s4  = MAX(chi5 - psi5, zero)
     cnh4hs4  = zero
     CALL calcmr                                       ! Water content
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   
   ! *** CALCULATE OBJECTIVE FUNCTION ************************************
   
   20    a4     = xk5 *(water/gama(2))**3.0
   funci3b= molal(5)*molal(2)*molal(2)/a4 - one

END FUNCTION funci3b

!!=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCI2
!! *** CASE I2
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!!     2. SOLID & LIQUID AEROSOL POSSIBLE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, NAHSO4, LC
!!     THERE ARE THREE REGIMES IN THIS CASE:
!!     1. NH4HSO4(s) POSSIBLE. LIQUID & SOLID AEROSOL (SUBROUTINE CALCI2A)
!!     2. NH4HSO4(s) NOT POSSIBLE, AND RH < MDRH. SOLID AEROSOL ONLY
!!     3. NH4HSO4(s) NOT POSSIBLE, AND RH >= MDRH. SOLID & LIQUID AEROSOL
!!     REGIMES 2. AND 3. ARE CONSIDERED TO BE THE SAME AS CASES I1A, I2B
!!     RESPECTIVELY
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calci2
   USE Isorropia_Module
   IMPLICIT NONE
   
   EXTERNAL calci1a, calci3a
   integer :: i
   
   ! *** FIND DRY COMPOSITION **********************************************
   
   CALL calci1a
   
   ! *** REGIME DEPENDS UPON THE POSSIBLE SOLIDS & RH **********************
   
   IF (cnh4hs4 > tiny) THEN
     scase = 'I2 ; SUBCASE 1'
     CALL calci2a
     scase = 'I2 ; SUBCASE 1'
   END IF
   
   IF (water <= tiny) THEN
     IF (rh < drmi2) THEN         ! SOLID SOLUTION ONLY
       water = tiny
       DO  i=1,nions
         molal(i) = zero
       END DO
       CALL calci1a
       scase = 'I2 ; SUBCASE 2'
       
     ELSE IF (rh >= drmi2) THEN     ! MDRH OF I2
       scase = 'I2 ; SUBCASE 3'
       CALL calcmdrh (rh, drmi2, drnahso4, calci1a, calci3a)
       scase = 'I2 ; SUBCASE 3'
     END IF
   END IF
   
END SUBROUTINE  calci2
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCI2A
!! *** CASE I2 ; SUBCASE A
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!!     2. SOLID & LIQUID AEROSOL POSSIBLE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, NAHSO4, LC
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calci2a
   USE Isorropia_Module
   use solut
   IMPLICIT NONE
   
   double precision :: psi2lo,psi2hi,x1,y1,yhi,dx,x2,y2,x3,y3
   double precision, external :: funci2a,funci3a
   integer :: i
   
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** FIND DRY COMPOSITION **********************************************
   
   CALL calci1a    ! Needed when called from CALCMDRH
   
   ! *** SETUP PARAMETERS ************************************************
   
   chi1 = cnh4hs4               ! Save from CALCI1 run
   chi2 = clc
   chi3 = cnahso4
   chi4 = cna2so4
   chi5 = cnh42s4
   
   psi1 = cnh4hs4               ! ASSIGN INITIAL PSI's
   psi2 = zero
   psi3 = zero
   psi4 = zero
   psi5 = zero
   
   calaou = .true.              ! Outer loop activity calculation flag
   psi2lo = zero                ! Low  limit
   psi2hi = chi2                ! High limit
   
   ! *** INITIAL VALUES FOR BISECTION ************************************
   
   x1 = psi2hi
   y1 = funci2a (x1)
   yhi= y1                      ! Save Y-value at HI position
   
   ! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH LC *********
   
   IF (yhi < eps) GO TO 50
   
   ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
   
   dx = (psi2hi-psi2lo)/FLOAT(ndiv)
   DO  i=1,ndiv
     x2 = MAX(x1-dx, psi2lo)
     y2 = funci2a (x2)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y2) < zero) GO TO 20  ! (Y1*Y2.LT.ZERO)
     x1 = x2
     y1 = y2
   END DO
   
   ! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH LC
   
   IF (y2 > eps) y2 = funci3a (zero)
   GO TO 50
   
   ! *** PERFORM BISECTION ***********************************************
   
   20    DO  i=1,maxit
     x3 = 0.5*(x1+x2)
     y3 = funci2a (x3)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y3) <= zero) THEN  ! (Y1*Y3 .LE. ZERO)
       y2    = y3
       x2    = x3
     ELSE
       y1    = y3
       x1    = x3
     END IF
     IF (ABS(x2-x1) <= eps*x1) GO TO 40
   END DO
   CALL pusherr (0002, 'CALCI2A')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CONVERGED ; RETURN **********************************************
   
   40    x3 = 0.5*(x1+x2)
   y3 = funci2a (x3)
   
   50    RETURN
   
   ! *** END OF SUBROUTINE CALCI2A *****************************************
   
END SUBROUTINE  calci2a
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE FUNCI2A
!! *** CASE I2 ; SUBCASE 1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!!     2. SOLID & LIQUID AEROSOL POSSIBLE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4, NAHSO4, LC
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!======================================================================
DOUBLE PRECISION FUNCTION funci2a (p2)
   USE Isorropia_Module
   use solut
   IMPLICIT NONE
   
   double precision , INTENT(IN)                         :: p2
   double precision :: aa,bb,cc,dd
   integer :: i,islv
   
   !LFR  COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
       !LFR  psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
       !LFR  a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   frst   = .true.
   calain = .true.
   psi2   = p2                  ! Save PSI2 in COMMON BLOCK
   psi3   = chi3
   psi4   = chi4
   psi5   = chi5
   psi6   = zero
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     
     a3 = xk11*(water/gama(9))**2.0
     a4 = xk5 *(water/gama(2))**3.0
     a5 = xk7 *(water/gama(4))**3.0
     a6 = xk1 *water/gama(7)*(gama(8)/gama(7))**2.
     a7 = SQRT(a4/a5)
     
   !  CALCULATE DISSOCIATION QUANTITIES
     
     IF (chi5 > tiny .AND. water > tiny) THEN
       psi5 = (psi3 + 2.d0*psi4 - a7*(3.d0*psi2 + psi1))/2.d0/a7
       psi5 = MAX(MIN (psi5, chi5), tiny)
     END IF
     
     IF (chi4 > tiny .AND. water > tiny) THEN
       aa   = psi2+psi5+psi6+psi3
       bb   = psi3*aa
       cc   = 0.25D0*(psi3*psi3*(psi2+psi5+psi6)-a4)
       CALL poly3 (aa, bb, cc, psi4, islv)
       IF (islv == 0) THEN
         psi4 = MIN (psi4, chi4)
       ELSE
         psi4 = zero
       END IF
     END IF
     
     IF (chi3 > tiny .AND. water > tiny) THEN
       aa   = 2.d0*psi4 + psi2 + psi1 - psi6
       bb   = 2.d0*psi4*(psi2 + psi1 - psi6) - a3
       cc   = zero
       CALL poly3 (aa, bb, cc, psi3, islv)
       IF (islv == 0) THEN
         psi3 = MIN (psi3, chi3)
       ELSE
         psi3 = zero
       END IF
     END IF
     
     bb   = psi2 + psi4 + psi5 + a6 +organion              ! PSI6
     cc   =-a6*(psi2 + psi3 + psi1)+ organion*(psi2+psi4+psi5)
     dd   = MAX(bb*bb-4.d0*cc, 0.d0)       ! US patch 12/20/01
     psi6 = 0.5D0*(-bb + SQRT(dd))
     
   ! *** CALCULATE SPECIATION ********************************************
     
     molal (1) = psi6 + organion                ! HI
     molal (2) = 2.d0*psi4 + psi3               ! NAI
     molal (3) = 3.d0*psi2 + 2.d0*psi5 + psi1   ! NH4I
     molal (5) = psi2 + psi4 + psi5 + psi6      ! SO4I
     molal (6) = psi2 + psi3 + psi1 - psi6      ! HSO4I
     clc       = chi2 - psi2
     cnahso4   = chi3 - psi3
     cna2so4   = chi4 - psi4
     cnh42s4   = chi5 - psi5
     cnh4hs4   = zero
     CALL calcmr                                ! Water content
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   
   ! *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
   
   20    a2      = xk13*(water/gama(13))**5.0
   funci2a = molal(5)*molal(6)*molal(3)**3.d0/a2 - one

END FUNCTION funci2a
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCI1
!! *** CASE I1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!!     2. SOLID AEROSOL ONLY
!!     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NA2SO4
!!     THERE ARE TWO POSSIBLE REGIMES HERE, DEPENDING ON RELATIVE HUMIDITY:
!!     1. WHEN RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!!     2. WHEN RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCI1A)
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calci1
   USE Isorropia_Module
   IMPLICIT NONE
   
   EXTERNAL calci1a, calci2a
   
   ! *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************
   
   IF (rh < drmi1) THEN
     scase = 'I1 ; SUBCASE 1'
     CALL calci1a              ! SOLID PHASE ONLY POSSIBLE
     scase = 'I1 ; SUBCASE 1'
   ELSE
     scase = 'I1 ; SUBCASE 2'  ! LIQUID & SOLID PHASE POSSIBLE
     CALL calcmdrh (rh, drmi1, drnh4hs4, calci1a, calci2a)
     scase = 'I1 ; SUBCASE 2'
   END IF
   
   ! *** AMMONIA IN GAS PHASE **********************************************
   
   CALL calcnh3
   
END SUBROUTINE  calci1
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCI1A
!! *** CASE I1 ; SUBCASE 1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
!!     2. SOLID AEROSOL ONLY
!!     3. SOLIDS POSSIBLE : NH4HSO4, NAHSO4, (NH4)2SO4, NA2SO4, LC
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calci1a
   USE Isorropia_Module
   IMPLICIT NONE
   
   double precision :: frso4,frnh4
   ! *** CALCULATE NON VOLATILE SOLIDS ***********************************
   
   cna2so4 = 0.5D0*w(1)
   cnh4hs4 = zero
   cnahso4 = zero
   cnh42s4 = zero
   frso4   = MAX(w(2)-cna2so4, zero)
   
   clc     = MIN(w(3)/3.d0, frso4/2.d0)
   frso4   = MAX(frso4-2.d0*clc, zero)
   frnh4   = MAX(w(3)-3.d0*clc,  zero)
   
   IF (frso4 <= tiny) THEN
     clc     = MAX(clc - frnh4, zero)
     cnh42s4 = 2.d0*frnh4
     
   ELSE IF (frnh4 <= tiny) THEN
     cnh4hs4 = 3.d0*MIN(frso4, clc)
     clc     = MAX(clc-frso4, zero)
     IF (cna2so4 > tiny) THEN
       frso4   = MAX(frso4-cnh4hs4/3.d0, zero)
       cnahso4 = 2.d0*frso4
       cna2so4 = MAX(cna2so4-frso4, zero)
     END IF
   END IF
   
   ! *** CALCULATE GAS SPECIES *********************************************
   
   ghno3 = w(4)
   ghcl  = w(5)
   gnh3  = zero
   
END SUBROUTINE  calci1a

!!=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCJ3
!! *** CASE J3
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
!!     2. THERE IS ONLY A LIQUID PHASE
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcj3
   USE Isorropia_Module
   IMPLICIT NONE
   
   DOUBLE PRECISION :: lamda, kapa,chi1,chi2,psi1,psi2,a3,bb,cc,dd
   integer :: i
   
   ! *** SETUP PARAMETERS ************************************************
   
   calaou = .true.              ! Outer loop activity calculation flag
   frst   = .true.
   calain = .true.
   
   lamda  = MAX(w(2) - w(3) - w(1), tiny)  ! FREE H2SO4
   chi1   = w(1)                           ! NA TOTAL as NaHSO4
   chi2   = w(3)                           ! NH4 TOTAL as NH4HSO4
   psi1   = chi1
   psi2   = chi2                           ! ALL NH4HSO4 DELIQUESCED
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     
     a3 = xk1  *water/gama(7)*(gama(8)/gama(7))**2.0
     
   !  CALCULATE DISSOCIATION QUANTITIES
     
     bb   = a3+lamda+organion                 ! KAPA
     cc   =-a3*(lamda + psi1 + psi2)
     dd   = MAX(bb*bb-4.d0*cc, 0.d0)       ! US PATCH 12/20/01
     kapa = 0.5D0*(-bb+SQRT(dd))
     
   ! *** CALCULATE SPECIATION ********************************************
     
     molal (1) = lamda + kapa + organion      ! HI
     molal (2) = psi1                         ! NAI
     molal (3) = psi2                         ! NH4I
     molal (4) = zero                         ! CLI
     molal (5) = kapa                         ! SO4I
     molal (6) = lamda + psi1 + psi2 - kapa   ! HSO4I
     molal (7) = zero                         ! NO3I
     CALL calcmr                              ! Water content
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 50
     END IF
   END DO
   
   50    RETURN
   
END SUBROUTINE  calcj3

!!=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCJ2
!! *** CASE J2
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : NAHSO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcj2
   USE Isorropia_Module
   USE casej
   IMPLICIT NONE
   
   DOUBLE PRECISION :: psi1lo,psi1hi,x1,y1,yhi,dx,x2,y2,ylo,x3,y3
   double precision, external :: funcj2
   integer :: i
   
   !LFR  COMMON /casej/ chi1, chi2, chi3, lamda, kapa, psi1, psi2, psi3, a1,   a2,   a3
   
   ! *** SETUP PARAMETERS ************************************************
   
   calaou = .true.              ! Outer loop activity calculation flag
   chi1   = w(1)                ! NA TOTAL
   chi2   = w(3)                ! NH4 TOTAL
   psi1lo = tiny                ! Low  limit
   psi1hi = chi1                ! High limit
   
   ! *** INITIAL VALUES FOR BISECTION ************************************
   
   x1 = psi1hi
   y1 = funcj2 (x1)
   yhi= y1                      ! Save Y-value at HI position
   
   ! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH42SO4 ****
   
   IF (ABS(y1) <= eps .OR. yhi < zero) GO TO 50
   
   ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
   
   dx = (psi1hi-psi1lo)/FLOAT(ndiv)
   DO  i=1,ndiv
     x2 = x1-dx
     y2 = funcj2 (x2)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y2) < zero) GO TO 20  ! (Y1*Y2.LT.ZERO)
     x1 = x2
     y1 = y2
   END DO
   
   ! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH42SO4
   
   ylo= y1                      ! Save Y-value at Hi position
   IF (ylo > zero .AND. yhi > zero) THEN
     y3 = funcj2 (zero)
     GO TO 50
   ELSE IF (ABS(y2) < eps) THEN   ! x2 IS A SOLUTION
     GO TO 50
   ELSE
     CALL pusherr (0001, 'CALCJ2')    ! WARNING ERROR: NO SOLUTION
     GO TO 50
   END IF
   
   ! *** PERFORM BISECTION ***********************************************
   
   20    DO  i=1,maxit
     x3 = 0.5*(x1+x2)
     y3 = funcj2 (x3)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y3) <= zero) THEN  ! (Y1*Y3 .LE. ZERO)
       y2    = y3
       x2    = x3
     ELSE
       y1    = y3
       x1    = x3
     END IF
     IF (ABS(x2-x1) <= eps*x1) GO TO 40
   END DO
   CALL pusherr (0002, 'CALCJ2')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CONVERGED ; RETURN **********************************************
   
   40    x3 = 0.5*(x1+x2)
   y3 = funcj2 (x3)
   
   50    RETURN

END SUBROUTINE  calcj2
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE FUNCJ2
!! *** CASE J2
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
DOUBLE PRECISION FUNCTION funcj2 (p1)
   USE Isorropia_Module
   USE casej
   IMPLICIT NONE
   
   double precision, INTENT(IN)                         :: p1
   double precision :: bb,cc,dd
   integer :: i
   !LFR  DOUBLE PRECISION :: lamda, kapa
   !LFR  COMMON /casej/ chi1, chi2, chi3, lamda, kapa, psi1, psi2, psi3, a1,   a2,   a3
   
   ! *** SETUP PARAMETERS ************************************************
   
   frst   = .true.
   calain = .true.
   
   lamda  = MAX(w(2) - w(3) - w(1), tiny)  ! FREE H2SO4
   psi1   = p1
   psi2   = chi2                           ! ALL NH4HSO4 DELIQUESCED
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     
     a1 = xk11 *(water/gama(12))**2.0
     a3 = xk1  *water/gama(7)*(gama(8)/gama(7))**2.0
     
   !  CALCULATE DISSOCIATION QUANTITIES
     
     bb   = a3+lamda+organion                  ! KAPA
     cc   =-a3*(lamda + psi1 + psi2)
     dd   = MAX(bb*bb-4.d0*cc, 0.d0)       ! US PATCH 12/20/01
     kapa = 0.5D0*(-bb+SQRT(dd))
     
   ! *** CALCULATE SPECIATION ********************************************
     
     molal (1) = lamda + kapa + organion       ! HI
     molal (2) = psi1                          ! NAI
     molal (3) = psi2                          ! NH4I
     molal (4) = zero                          ! CLI
     molal (5) = kapa                          ! SO4I
     molal (6) = lamda + psi1 + psi2 - kapa    ! HSO4I
     molal (7) = zero                          ! NO3I
     CALL calcmr                               ! Water content
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   
   ! *** CALCULATE OBJECTIVE FUNCTION ************************************
   
   20    funcj2 = molal(2)*molal(6)/a1 - one
   
END FUNCTION funcj2
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCJ1
!! *** CASE J1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : NH4HSO4, NAHSO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcj1
   USE Isorropia_Module
   USE casej
   IMPLICIT NONE
   
   DOUBLE PRECISION :: psi1lo,psi1hi,x1,y1,yhi,dx,x2,y2,ylo,x3,y3
   double precision, external :: funcj1
   integer :: i
   
   !LFR  COMMON /casej/ chi1, chi2, chi3, lamda, kapa, psi1, psi2, psi3, a1,   a2,   a3
   
   ! *** SETUP PARAMETERS ************************************************
   
   calaou =.true.               ! Outer loop activity calculation flag
   chi1   = w(1)                ! Total NA initially as NaHSO4
   chi2   = w(3)                ! Total NH4 initially as NH4HSO4
   
   psi1lo = tiny                ! Low  limit
   psi1hi = chi1                ! High limit
   
   ! *** INITIAL VALUES FOR BISECTION ************************************
   
   x1 = psi1hi
   y1 = funcj1 (x1)
   yhi= y1                      ! Save Y-value at HI position
   
   ! *** YHI < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH42SO4 ****
   
   IF (ABS(y1) <= eps .OR. yhi < zero) GO TO 50
   
   ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
   
   dx = (psi1hi-psi1lo)/FLOAT(ndiv)
   DO  i=1,ndiv
     x2 = x1-dx
     y2 = funcj1 (x2)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y2) < zero) GO TO 20  ! (Y1*Y2.LT.ZERO)
     x1 = x2
     y1 = y2
   END DO
   
   ! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH42SO4
   
   ylo= y1                      ! Save Y-value at Hi position
   IF (ylo > zero .AND. yhi > zero) THEN
     y3 = funcj1 (zero)
     GO TO 50
   ELSE IF (ABS(y2) < eps) THEN   ! x2 IS A SOLUTION
     GO TO 50
   ELSE
     CALL pusherr (0001, 'CALCJ1')    ! WARNING ERROR: NO SOLUTION
     GO TO 50
   END IF
   
   ! *** PERFORM BISECTION ***********************************************
   
   20    DO  i=1,maxit
     x3 = 0.5*(x1+x2)
     y3 = funcj1 (x3)
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y3) <= zero) THEN  ! (Y1*Y3 .LE. ZERO)
       y2    = y3
       x2    = x3
     ELSE
       y1    = y3
       x1    = x3
     END IF
     IF (ABS(x2-x1) <= eps*x1) GO TO 40
   END DO
   CALL pusherr (0002, 'CALCJ1')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CONVERGED ; RETURN **********************************************
   
   40    x3 = 0.5*(x1+x2)
   y3 = funcj1 (x3)
   
   50    RETURN
   
END SUBROUTINE  calcj1
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE FUNCJ1
!! *** CASE J1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
DOUBLE PRECISION FUNCTION funcj1 (p1)
   USE Isorropia_Module
   USE casej
   IMPLICIT NONE
   
   double precision, INTENT(IN)                         :: p1
   double precision :: bb,cc,dd
   integer :: i
   !LFR  DOUBLE PRECISION :: lamda, kapa
   !LFR  COMMON /casej/ chi1, chi2, chi3, lamda, kapa, psi1, psi2, psi3, a1,   a2,   a3
   
   ! *** SETUP PARAMETERS ************************************************
   
   frst   = .true.
   calain = .true.
   
   lamda  = MAX(w(2) - w(3) - w(1), tiny)  ! FREE H2SO4
   psi1   = p1
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     
     a1 = xk11 *(water/gama(12))**2.0
     a2 = xk12 *(water/gama(09))**2.0
     a3 = xk1  *water/gama(7)*(gama(8)/gama(7))**2.0
     
     psi2 = 0.5*(-(lamda+psi1) + SQRT((lamda+psi1)**2.d0+4.d0*a2))  ! PSI2
     psi2 = MIN (psi2, chi2)
     
     bb   = a3+lamda+organion                   ! KAPA
     cc   =-a3*(lamda + psi2 + psi1)
     dd   = MAX(bb*bb-4.d0*cc, 0.d0)       ! US PATCH 12/20/01
     kapa = 0.5D0*(-bb+SQRT(dd))
     
   ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
     
     molal (1) = lamda + kapa + organion       ! HI
     molal (2) = psi1                          ! NAI
     molal (3) = psi2                          ! NH4I
     molal (4) = zero
     molal (5) = kapa                          ! SO4I
     molal (6) = lamda + psi1 + psi2 - kapa    ! HSO4I
     molal (7) = zero
     CALL calcmr                               ! Water content
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   
   ! *** CALCULATE OBJECTIVE FUNCTION ************************************
   
   20    funcj1 = molal(2)*molal(6)/a1 - one
   
END FUNCTION funcj1
   
   
