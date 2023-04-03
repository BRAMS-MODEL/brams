!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE ISRP1R
!! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE REVERSE PROBLEM OF
!!     AN AMMONIUM-SULFATE AEROSOL SYSTEM.
!!     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY
!!     THE AMBIENT RELATIVE HUMIDITY.
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!!  @author Athanasios Nenes
!!  REVISION HISTORY:                                                   *
!!   Modified by Yang Zhang of AER in Sept., 2001 to incorporate the CMU*
!!           hybrid mass transfer module into MADRID1                   *
!!           The code for the CMU hybrid mass transfer module was       *
!!           provided by Bonyoung Koo/Spyros Pandis of CMU on           *
!!           May 18, 2001                                               *
!!                                                                       *
!!***********************************************************************
SUBROUTINE isrp1r (wi, rhi, tempi)
   USE Isorropia_Module
   IMPLICIT NONE
    
   double precision, INTENT(IN OUT)                     :: wi(ncomp)
   double precision, INTENT(IN OUT)                     :: rhi
   double precision, INTENT(IN OUT)                     :: tempi
      
      
   ! *** INITIALIZE COMMON BLOCK VARIABLES *********************************
   
   CALL init1 (wi, rhi, tempi)
   
   ! *** CALCULATE SULFATE RATIO *******************************************
   
   IF (rh >= drnh42s4) THEN         ! WET AEROSOL, NEED NH4 AT SRATIO=2.0
     sulratw = getasr(waer(2), rhi)     ! AEROSOL SULFATE RATIO
   ELSE
     sulratw = 2.0D0                    ! DRY AEROSOL SULFATE RATIO
   END IF
   sulrat  = waer(3)/waer(2)         ! SULFATE RATIO
   
   ! *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
   
   ! *** SULFATE POOR
   
   IF (sulratw <= sulrat) THEN
     
     IF(metstbl == 1) THEN
       scase = 'K2'
       CALL calck2                 ! Only liquid (metastable)
     ELSE
       
       IF (rh < drnh42s4) THEN
         scase = 'K1'
         CALL calck1              ! NH42SO4              ; case K1
         
       ELSE IF (drnh42s4 <= rh) THEN
         scase = 'K2'
         CALL calck2              ! Only liquid          ; case K2
       END IF
     END IF
     
   ! *** SULFATE RICH (NO ACID)
     
   ELSE IF (1.0 <= sulrat .AND. sulrat < sulratw) THEN
     w(2) = waer(2)
     w(3) = waer(3)
     
     IF(metstbl == 1) THEN
       scase = 'B4'
       CALL calcb4                 ! Only liquid (metastable)
       scase = 'L4'
     ELSE
       
       IF (rh < drnh4hs4) THEN
         scase = 'B1'
         CALL calcb1              ! NH4HSO4,LC,NH42SO4   ; case B1
         scase = 'L1'
         
       ELSE IF (drnh4hs4 <= rh .AND. rh < drlc) THEN
         scase = 'B2'
         CALL calcb2              ! LC,NH42S4            ; case B2
         scase = 'L2'
         
       ELSE IF (drlc <= rh .AND. rh < drnh42s4) THEN
         scase = 'B3'
         CALL calcb3              ! NH42S4               ; case B3
         scase = 'L3'
         
       ELSE IF (drnh42s4 <= rh) THEN
         scase = 'B4'
         CALL calcb4              ! Only liquid          ; case B4
         scase = 'L4'
       END IF
     END IF
     
     CALL calcnh3p
     
   ! *** SULFATE RICH (FREE ACID)
     
   ELSE IF (sulrat < 1.0) THEN
     w(2) = waer(2)
     w(3) = waer(3)
     
     IF(metstbl == 1) THEN
       scase = 'C2'
       CALL calcc2                 ! Only liquid (metastable)
       scase = 'M2'
     ELSE
       
       IF (rh < drnh4hs4) THEN
         scase = 'C1'
         CALL calcc1              ! NH4HSO4              ; case C1
         scase = 'M1'
         
       ELSE IF (drnh4hs4 <= rh) THEN
         scase = 'C2'
         CALL calcc2              ! Only liquid          ; case C2
         scase = 'M2'
       END IF
     END IF
     
     CALL calcnh3p
     
   END IF
END SUBROUTINE isrp1r
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE ISRP2R
!! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE REVERSE PROBLEM OF
!!     AN AMMONIUM-SULFATE-NITRATE AEROSOL SYSTEM.
!!     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY
!!     THE AMBIENT RELATIVE HUMIDITY.
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE isrp2r (wi, rhi, tempi)
   USE Isorropia_Module
   IMPLICIT NONE
   
   double precision, INTENT(IN OUT)                     :: wi(ncomp)
   double precision, INTENT(IN OUT)                     :: rhi
   double precision, INTENT(IN OUT)                     :: tempi
   
   LOGICAL :: tryliq
   
   ! *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
   
   tryliq = .true.             ! Assume liquid phase, sulfate poor limit
   
   10    CALL init2 (wi, rhi, tempi)
   
   ! *** CALCULATE SULFATE RATIO *******************************************
   
   IF (tryliq .AND. rh >= drnh4no3) THEN ! *** WET AEROSOL
     sulratw = getasr(waer(2), rhi)     ! LIMITING SULFATE RATIO
   ELSE
     sulratw = 2.0D0                    ! *** DRY AEROSOL
   END IF
   sulrat = waer(3)/waer(2)
   
   ! *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
   
   ! *** SULFATE POOR
   
   IF (sulratw <= sulrat) THEN
     
     IF(metstbl == 1) THEN
       scase = 'N3'
       CALL calcn3                 ! Only liquid (metastable)
     ELSE
       
       IF (rh < drnh4no3) THEN
         scase = 'N1'
         CALL calcn1              ! NH42SO4,NH4NO3       ; case N1
         
       ELSE IF (drnh4no3 <= rh .AND. rh < drnh42s4) THEN
         scase = 'N2'
         CALL calcn2              ! NH42S4               ; case N2
         
       ELSE IF (drnh42s4 <= rh) THEN
         scase = 'N3'
         CALL calcn3              ! Only liquid          ; case N3
       END IF
     END IF
     
   ! *** SULFATE RICH (NO ACID)
     
   !     FOR SOLVING THIS CASE, NITRIC ACID AND AMMONIA IN THE GAS PHASE ARE
   !     ASSUMED A MINOR SPECIES, THAT DO NOT SIGNIFICANTLY AFFECT THE
   !     AEROSOL EQUILIBRIUM.
     
   ELSE IF (1.0 <= sulrat .AND. sulrat < sulratw) THEN
     w(2) = waer(2)
     w(3) = waer(3)
     w(4) = waer(4)
     
     IF(metstbl == 1) THEN
       scase = 'B4'
       CALL calcb4                 ! Only liquid (metastable)
       scase = 'O4'
     ELSE
       
       IF (rh < drnh4hs4) THEN
         scase = 'B1'
         CALL calcb1              ! NH4HSO4,LC,NH42SO4   ; case O1
         scase = 'O1'
         
       ELSE IF (drnh4hs4 <= rh .AND. rh < drlc) THEN
         scase = 'B2'
         CALL calcb2              ! LC,NH42S4            ; case O2
         scase = 'O2'
         
       ELSE IF (drlc <= rh .AND. rh < drnh42s4) THEN
         scase = 'B3'
         CALL calcb3              ! NH42S4               ; case O3
         scase = 'O3'
         
       ELSE IF (drnh42s4 <= rh) THEN
         scase = 'B4'
         CALL calcb4              ! Only liquid          ; case O4
         scase = 'O4'
       END IF
     END IF
     
     CALL calcnap                ! HNO3, NH3 dissolved
     CALL calcnh3p
     
   ! *** SULFATE RICH (FREE ACID)
     
   !     FOR SOLVING THIS CASE, NITRIC ACID AND AMMONIA IN THE GAS PHASE ARE
   !     ASSUMED A MINOR SPECIES, THAT DO NOT SIGNIFICANTLY AFFECT THE
   !     AEROSOL EQUILIBRIUM.
     
   ELSE IF (sulrat < 1.0) THEN
     w(2) = waer(2)
     w(3) = waer(3)
     w(4) = waer(4)
     
     IF(metstbl == 1) THEN
       scase = 'C2'
       CALL calcc2                 ! Only liquid (metastable)
       scase = 'P2'
     ELSE
       
       IF (rh < drnh4hs4) THEN
         scase = 'C1'
         CALL calcc1              ! NH4HSO4              ; case P1
         scase = 'P1'
         
       ELSE IF (drnh4hs4 <= rh) THEN
         scase = 'C2'
         CALL calcc2              ! Only liquid          ; case P2
         scase = 'P2'
       END IF
     END IF
     
     CALL calcnap                ! HNO3, NH3 dissolved
     CALL calcnh3p
   END IF
   
   ! *** IF SULRATW < SULRAT < 2.0 and WATER = 0 => SULFATE RICH CASE.
   
   IF (sulratw <= sulrat .AND. sulrat < 2.0 .AND. water <= tiny) THEN
     tryliq = .false.
     GO TO 10
   END IF
   
END SUBROUTINE isrp2r

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE ISRP3R
!! *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE REVERSE PROBLEM OF
!!     AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM AEROSOL SYSTEM.
!!     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE & SODIUM
!!     RATIOS AND BY THE AMBIENT RELATIVE HUMIDITY.
!! *** COPYRIGHT 1996-2000 UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE isrp3r (wi, rhi, tempi)
   USE Isorropia_Module
   IMPLICIT NONE
   
   double precision, INTENT(IN OUT)                     :: wi(ncomp)
   double precision, INTENT(IN OUT)                     :: rhi
   double precision, INTENT(IN OUT)                     :: tempi
   
   LOGICAL :: tryliq
   double precision :: frso4,sri
   INTEGER :: i
   
   !cC
   !cC *** ADJUST FOR TOO LITTLE AMMONIUM AND CHLORIDE ***********************
   !cC
   !c      WI(3) = MAX (WI(3), 1.D-10)  ! NH4+ : 1e-4 umoles/m3
   !c      WI(5) = MAX (WI(5), 1.D-10)  ! Cl-  : 1e-4 umoles/m3
   
   ! *** INITIALIZE ALL VARIABLES ******************************************
   
   tryliq = .true.             ! Use liquid phase sulfate poor limit
   
   10    CALL isoinit3 (wi, rhi, tempi) ! COMMON block variables
   !cC
   !cC *** CHECK IF TOO MUCH SODIUM ; ADJUST AND ISSUE ERROR MESSAGE *********
   !cC
   !c      REST = 2.D0*WAER(2) + WAER(4) + WAER(5)
   !c      IF (WAER(1).GT.REST) THEN            ! NA > 2*SO4+CL+NO3 ?
   !c         WAER(1) = (ONE-1D-6)*REST         ! Adjust Na amount
   !c         CALL PUSHERR (0050, 'ISRP3R')     ! Warning error: Na adjusted
   !c      ENDIF
   
   ! *** CALCULATE SULFATE & SODIUM RATIOS *********************************
   
   IF (tryliq .AND. rh >= drnh4no3) THEN  ! ** WET AEROSOL
     frso4   = waer(2) - waer(1)/2.0D0     ! SULFATE UNBOUND BY SODIUM
     frso4   = MAX(frso4, tiny)
     sri     = getasr(frso4, rhi)          ! SULFATE RATIO FOR NH4+
     sulratw = (waer(1)+frso4*sri)/waer(2) ! LIMITING SULFATE RATIO
     sulratw = MIN (sulratw, 2.0D0)
   ELSE
     sulratw = 2.0D0                     ! ** DRY AEROSOL
   END IF
   sulrat = (waer(1)+waer(3))/waer(2)
   sodrat = waer(1)/waer(2)
   
   ! *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
   
   ! *** SULFATE POOR ; SODIUM POOR
   
   IF (sulratw <= sulrat .AND. sodrat < 2.0) THEN
     
     IF(metstbl == 1) THEN
       scase = 'Q5'
       CALL calcq5                 ! Only liquid (metastable)
       scase = 'Q5'
     ELSE
       
       IF (rh < drnh4no3) THEN
         scase = 'Q1'
         CALL calcq1              ! NH42SO4,NH4NO3,NH4CL,NA2SO4
         
       ELSE IF (drnh4no3 <= rh .AND. rh < drnh4cl) THEN
         scase = 'Q2'
         CALL calcq2              ! NH42SO4,NH4CL,NA2SO4
         
       ELSE IF (drnh4cl <= rh  .AND. rh < drnh42s4) THEN
         scase = 'Q3'
         CALL calcq3              ! NH42SO4,NA2SO4
         
       ELSE IF (drnh42s4 <= rh  .AND. rh < drna2so4) THEN
         scase = 'Q4'
         CALL calcq4              ! NA2SO4
         scase = 'Q4'
         
       ELSE IF (drna2so4 <= rh) THEN
         scase = 'Q5'
         CALL calcq5              ! Only liquid
         scase = 'Q5'
       END IF
     END IF
     
   ! *** SULFATE POOR ; SODIUM RICH
     
   ELSE IF (sulrat >= sulratw .AND. sodrat >= 2.0) THEN
     
     IF(metstbl == 1) THEN
       scase = 'R6'
       CALL calcr6                 ! Only liquid (metastable)
       scase = 'R6'
     ELSE
       
       IF (rh < drnh4no3) THEN
         scase = 'R1'
         CALL calcr1              ! NH4NO3,NH4CL,NA2SO4,NACL,NANO3
         
       ELSE IF (drnh4no3 <= rh .AND. rh < drnano3) THEN
         scase = 'R2'
         CALL calcr2              ! NH4CL,NA2SO4,NACL,NANO3
         
       ELSE IF (drnano3 <= rh  .AND. rh < drnacl) THEN
         scase = 'R3'
         CALL calcr3              ! NH4CL,NA2SO4,NACL
         
       ELSE IF (drnacl <= rh   .AND. rh < drnh4cl) THEN
         scase = 'R4'
         CALL calcr4              ! NH4CL,NA2SO4
         
       ELSE IF (drnh4cl <= rh .AND. rh < drna2so4) THEN
         scase = 'R5'
         CALL calcr5              ! NA2SO4
         scase = 'R5'
         
       ELSE IF (drna2so4 <= rh) THEN
         scase = 'R6'
         CALL calcr6              ! NO SOLID
         scase = 'R6'
       END IF
     END IF
     
   ! *** SULFATE RICH (NO ACID)
     
   ELSE IF (1.0 <= sulrat .AND. sulrat < sulratw) THEN
     DO  i=1,ncomp
       w(i) = waer(i)
     END DO
     
     IF(metstbl == 1) THEN
       scase = 'I6'
       CALL calci6                 ! Only liquid (metastable)
       scase = 'S6'
     ELSE
       
       IF (rh < drnh4hs4) THEN
         scase = 'I1'
         CALL calci1              ! NA2SO4,(NH4)2SO4,NAHSO4,NH4HSO4,LC
         scase = 'S1'
         
       ELSE IF (drnh4hs4 <= rh .AND. rh < drnahso4) THEN
         scase = 'I2'
         CALL calci2              ! NA2SO4,(NH4)2SO4,NAHSO4,LC
         scase = 'S2'
         
       ELSE IF (drnahso4 <= rh .AND. rh < drlc) THEN
         scase = 'I3'
         CALL calci3              ! NA2SO4,(NH4)2SO4,LC
         scase = 'S3'
         
       ELSE IF (drlc <= rh     .AND. rh < drnh42s4) THEN
         scase = 'I4'
         CALL calci4              ! NA2SO4,(NH4)2SO4
         scase = 'S4'
         
       ELSE IF (drnh42s4 <= rh .AND. rh < drna2so4) THEN
         scase = 'I5'
         CALL calci5              ! NA2SO4
         scase = 'S5'
         
       ELSE IF (drna2so4 <= rh) THEN
         scase = 'I6'
         CALL calci6              ! NO SOLIDS
         scase = 'S6'
       END IF
     END IF
     
     CALL calcnhp                ! HNO3, NH3, HCL in gas phase
     CALL calcnh3p
     
   ! *** SULFATE RICH (FREE ACID)
     
   ELSE IF (sulrat < 1.0) THEN
     DO  i=1,ncomp
       w(i) = waer(i)
     END DO
     
     IF(metstbl == 1) THEN
       scase = 'J3'
       CALL calcj3                 ! Only liquid (metastable)
       scase = 'T3'
     ELSE
       
       IF (rh < drnh4hs4) THEN
         scase = 'J1'
         CALL calcj1              ! NH4HSO4,NAHSO4
         scase = 'T1'
         
       ELSE IF (drnh4hs4 <= rh .AND. rh < drnahso4) THEN
         scase = 'J2'
         CALL calcj2              ! NAHSO4
         scase = 'T2'
         
       ELSE IF (drnahso4 <= rh) THEN
         scase = 'J3'
         CALL calcj3
         scase = 'T3'
       END IF
     END IF
     
     CALL calcnhp                ! HNO3, NH3, HCL in gas phase
     CALL calcnh3p
     
   END IF
   
   ! *** IF AFTER CALCULATIONS, SULRATW < SULRAT < 2.0
   !                            and WATER = 0          => SULFATE RICH CASE.
   
   IF (sulratw <= sulrat .AND. sulrat < 2.0 .AND. water <= tiny) THEN
     tryliq = .false.
     GO TO 10
   END IF
   
END SUBROUTINE isrp3r

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCK2
!! *** CASE K2
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0)
!!     2. LIQUID AEROSOL PHASE ONLY POSSIBLE
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calck2
   USE Isorropia_Module
   implicit none
   
   DOUBLE PRECISION :: nh4i, nh3gi, nh3aq, hi,ohi,del,so4i,a2,akw,hso4i
   integer :: i
   
   ! *** SETUP PARAMETERS ************************************************
   
   calaou   =.true.     ! Outer loop activity calculation flag
   frst     =.true.
   calain   =.true.
   
   ! *** CALCULATE WATER CONTENT *****************************************
   
   molalr(4)= MIN(waer(2), 0.5D0*waer(3))
   water    = molalr(4)/m0(4)  ! ZSR correlation
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
   !C         A21  = XK21*WATER*R*TEMP
     a2   = xk2 *r*temp/xkw/rh*(gama(8)/gama(9))**2.
     akw  = xkw *rh*water*water
     
     nh4i = waer(3)
     so4i = waer(2)
     hso4i= zero
     
     CALL calcph (2.d0*so4i - nh4i, hi, ohi)    ! Get pH
     
     nh3aq = zero                               ! AMMONIA EQUILIBRIUM
     IF (hi < ohi) THEN
       CALL calcamaq (nh4i, ohi, del)
       nh4i  = MAX (nh4i-del, zero)
       ohi   = MAX (ohi -del, tiny)
       nh3aq = del
       hi    = akw/ohi
     END IF
     
     CALL calchs4 (hi, so4i, zero, del)         ! SULFATE EQUILIBRIUM
     so4i  = so4i - del
     hi    = hi   - del
     hso4i = del
     
     nh3gi = nh4i/hi/a2   !    NH3AQ/A21
     
   ! *** SPECIATION & WATER CONTENT ***************************************
     
     molal(1) = hi
     molal(3) = nh4i
     molal(5) = so4i
     molal(6) = hso4i
     coh      = ohi
     gasaq(1) = nh3aq
     gnh3     = nh3gi
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   
   20    RETURN
   
END SUBROUTINE calck2

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCK1
!! *** CASE K1
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
SUBROUTINE calck1
   USE Isorropia_Module
   implicit none
   
   cnh42s4 = MIN(waer(2),0.5D0*waer(3))  ! For bad input problems
   gnh3    = zero
   
   w(2)    = cnh42s4
   w(3)    = 2.d0*cnh42s4 + gnh3
   
END SUBROUTINE calck1
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCN3
!! *** CASE N3
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0)
!!     2. THERE IS ONLY A LIQUID PHASE
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcn3
   USE Isorropia_Module
   USE solut
   implicit none
   
   DOUBLE PRECISION :: nh4i, no3i, nh3aq, no3aq,hi,ohi,gg,so4i,del
   double precision :: aml5,akw,hso4i
   integer :: i
!   COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
!       psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
!       a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   calaou =.true.              ! Outer loop activity calculation flag
   frst   =.true.
   calain =.true.
   
   ! *** AEROSOL WATER CONTENT
   
   molalr(4) = MIN(waer(2),0.5D0*waer(3))       ! (NH4)2SO4
   aml5      = MAX(waer(3)-2.d0*molalr(4),zero) ! "free" NH4
   molalr(5) = MAX(MIN(aml5,waer(4)), zero)     ! NH4NO3=MIN("free",NO3)
   water     = molalr(4)/m0(4) + molalr(5)/m0(5)
   water     = MAX(water, tiny)
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     a2    = xk2 *r*temp/xkw/rh*(gama(8)/gama(9))**2.
   !C         A21   = XK21*WATER*R*TEMP
     a3    = xk4*r*temp*(water/gama(10))**2.0
     a4    = xk7*(water/gama(4))**3.0
     akw   = xkw *rh*water*water
     
   ! ION CONCENTRATIONS
     
     nh4i  = waer(3)
     no3i  = waer(4)
     so4i  = waer(2)
     hso4i = zero
     
     CALL calcph (2.d0*so4i + no3i - nh4i, hi, ohi)
     
   ! AMMONIA ASSOCIATION EQUILIBRIUM
     
     nh3aq = zero
     no3aq = zero
     gg    = 2.d0*so4i + no3i - nh4i
     IF (hi < ohi) THEN
       CALL calcamaq2 (-gg, nh4i, ohi, nh3aq)
       hi    = akw/ohi
     ELSE
       hi    = zero
       CALL calcniaq2 (gg, no3i, hi, no3aq) ! HNO3
       
   ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
       
       CALL calchs4 (hi, so4i, zero, del)
       so4i  = so4i  - del
       hi    = hi    - del
       hso4i = del
       ohi   = akw/hi
     END IF
     
   ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
     
     molal (1) = hi
     molal (3) = nh4i
     molal (5) = so4i
     molal (6) = hso4i
     molal (7) = no3i
     coh       = ohi
     
     cnh42s4   = zero
     cnh4no3   = zero
     
     gasaq(1)  = nh3aq
     gasaq(3)  = no3aq
     
     ghno3     = hi*no3i/a3
     gnh3      = nh4i/hi/a2   !   NH3AQ/A21
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP ******************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   
   ! *** RETURN ***********************************************************
   
   20    RETURN
   
END SUBROUTINE calcn3

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCN2
!! *** CASE N2
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0)
!!     2. THERE IS BOTH A LIQUID & SOLID PHASE
!!     3. SOLIDS POSSIBLE : (NH4)2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcn2
   USE Isorropia_Module
   USE solut
   
   implicit none
   
   double precision :: psi1lo,psi1hi,x1,y1,yhi,x2,y2,ylo,p4,yy,x3,y3,dx
   double precision, external :: funcn2
   integer :: i
   
!   COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
!       psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
!       a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   chi1   = MIN(waer(2),0.5D0*waer(3))     ! (NH4)2SO4
   chi2   = MAX(waer(3) - 2.d0*chi1, zero) ! "Free" NH4+
   chi3   = MAX(waer(4) - chi2, zero)      ! "Free" NO3
   
   psi2   = chi2
   psi3   = chi3
   
   calaou = .true.              ! Outer loop activity calculation flag
   psi1lo = tiny                ! Low  limit
   psi1hi = chi1                ! High limit
   
   ! *** INITIAL VALUES FOR BISECTION ************************************
   
   x1 = psi1hi
   y1 = funcn2 (real(x1))
   IF (y1 <= eps) RETURN   ! IF (ABS(Y1).LE.EPS .OR. Y1.LE.ZERO) RETURN
   yhi= y1                 ! Save Y-value at HI position
   
   ! *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
   
   dx = (psi1hi-psi1lo)/FLOAT(ndiv)
   DO  i=1,ndiv
     x2 = MAX(x1-dx, zero)
     y2 = funcn2 (real(x2))
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y2) < zero) GO TO 20  ! (Y1*Y2.LT.ZERO)
     x1 = x2
     y1 = y2
   END DO
   
   ! *** NO SUBDIVISION WITH SOLUTION FOUND
   
   ylo= y1                      ! Save Y-value at Hi position
   IF (ABS(y2) < eps) THEN   ! x2 IS A SOLUTION
     RETURN
     
   ! *** { YLO, YHI } < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH3
     
   ELSE IF (ylo < zero .AND. yhi < zero) THEN
     p4 = chi4
     yy = funcn2(real(p4))
     GO TO 50
     
   ! *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH3
     
   ELSE IF (ylo > zero .AND. yhi > zero) THEN
     p4 = tiny
     yy = funcn2(real(p4))
     GO TO 50
   ELSE
     CALL pusherr (0001, 'CALCN2')    ! WARNING ERROR: NO SOLUTION
     RETURN
   END IF
   
   ! *** PERFORM BISECTION ***********************************************
   
   20    DO  i=1,maxit
     x3 = 0.5*(x1+x2)
     y3 = funcn2 (real(x3))
     IF (SIGN(1.d0,y1)*SIGN(1.d0,y3) <= zero) THEN  ! (Y1*Y3 .LE. ZERO)
       y2    = y3
       x2    = x3
     ELSE
       y1    = y3
       x1    = x3
     END IF
     IF (ABS(x2-x1) <= eps*x1) GO TO 40
   END DO
   CALL pusherr (0002, 'CALCN2')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CONVERGED ; RETURN **********************************************
   
   40    x3 = 0.5*(x1+x2)
   y3 = funcn2 (real(x3))
   50    CONTINUE

END SUBROUTINE calcn2
   
!>======================================================================
!! *** ISORROPIA CODE
!! *** FUNCTION FUNCN2
!! *** CASE D2
!!     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE D2 ;
!!     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCN2.
!!=======================================================================
DOUBLE PRECISION FUNCTION funcn2 (p1)
   USE Isorropia_Module
   use solut
   implicit none
   
   REAL, INTENT(IN)                         :: p1
   DOUBLE PRECISION :: nh4i, no3i, nh3aq, no3aq, hi, ohi,gg,so4i,del
   double precision :: akw,hso4i
   integer :: i

!   COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
!       psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
!       a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   frst   = .true.
   calain = .true.
   psi1   = p1
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     a2    = xk2 *r*temp/xkw/rh*(gama(8)/gama(9))**2.
   !C         A21   = XK21*WATER*R*TEMP
     a3    = xk4*r*temp*(water/gama(10))**2.0
     a4    = xk7*(water/gama(4))**3.0
     akw   = xkw *rh*water*water
     
   ! ION CONCENTRATIONS
     
     nh4i  = 2.d0*psi1 + psi2
     no3i  = psi2 + psi3
     so4i  = psi1
     hso4i = zero
     
     CALL calcph (2.d0*so4i + no3i - nh4i, hi, ohi)
     
   ! AMMONIA ASSOCIATION EQUILIBRIUM
     
     nh3aq = zero
     no3aq = zero
     gg    = 2.d0*so4i + no3i - nh4i
     IF (hi < ohi) THEN
       CALL calcamaq2 (-gg, nh4i, ohi, nh3aq)
       hi    = akw/ohi
     ELSE
       hi    = zero
       CALL calcniaq2 (gg, no3i, hi, no3aq) ! HNO3
       
   ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
       
       CALL calchs4 (hi, so4i, zero, del)
       so4i  = so4i  - del
       hi    = hi    - del
       hso4i = del
       ohi   = akw/hi
     END IF
     
   ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
     
     molal (1) = hi
     molal (3) = nh4i
     molal (5) = so4i
     molal (6) = hso4i
     molal (7) = no3i
     coh       = ohi
     
     cnh42s4   = chi1 - psi1
     cnh4no3   = zero
     
     gasaq(1)  = nh3aq
     gasaq(3)  = no3aq
     
     ghno3     = hi*no3i/a3
     gnh3      = nh4i/hi/a2   !   NH3AQ/A21
     
   ! *** CALCULATE MOLALR ARRAY, WATER AND ACTIVITIES **********************
     
     CALL calcmr
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   
   ! *** CALCULATE OBJECTIVE FUNCTION ************************************
   
   20    funcn2= nh4i*nh4i*so4i/a4 - one

END FUNCTION funcn2

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCN1
!! *** CASE N1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0)
!!     2. SOLID AEROSOL ONLY
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3
!!     THERE ARE TWO REGIMES DEFINED BY RELATIVE HUMIDITY:
!!     1. RH < MDRH  ; ONLY SOLID PHASE POSSIBLE (SUBROUTINE CALCN1A)
!!     2. RH >= MDRH ; LIQUID PHASE POSSIBLE (MDRH REGION)
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcn1
   USE Isorropia_Module
   implicit none
   
   EXTERNAL calcn1a, calcn2
   
   ! *** REGIME DEPENDS UPON THE AMBIENT RELATIVE HUMIDITY *****************
   
   IF (rh < drmasan) THEN
     scase = 'N1 ; SUBCASE 1'
     CALL calcn1a              ! SOLID PHASE ONLY POSSIBLE
     scase = 'N1 ; SUBCASE 1'
   ELSE
     scase = 'N1 ; SUBCASE 2'
     CALL calcmdrp (rh, drmasan, drnh4no3, calcn1a, calcn2)
     scase = 'N1 ; SUBCASE 2'
   END IF
   
END SUBROUTINE calcn1
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCN1A
!! *** CASE N1 ; SUBCASE 1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0)
!!     2. SOLID AEROSOL ONLY
!!     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcn1a
   USE Isorropia_Module
   implicit none
   
   double precision :: psi1,psi2
   ! *** SETUP PARAMETERS *************************************************
   
   !CC      A1      = XK10/R/TEMP/R/TEMP
   
   ! *** CALCULATE AEROSOL COMPOSITION ************************************
   
   !CC      CHI1    = 2.D0*WAER(4)        ! Free parameter ; arbitrary value.
   psi1    = waer(4)
   
   ! *** The following statment is here to avoid negative NH4+ values in
   !     CALCN? routines that call CALCN1A
   
   psi2    = MAX(MIN(waer(2),0.5D0*(waer(3)-psi1)),tiny)
   
   cnh4no3 = psi1
   cnh42s4 = psi2
   
   !CC      GNH3    = CHI1 + PSI1 + 2.0*PSI2
   !CC      GHNO3   = A1/(CHI1-PSI1) + PSI1
   gnh3    = zero
   ghno3   = zero
   
   w(2)    = psi2
   w(3)    = gnh3  + psi1 + 2.0*psi2
   w(4)    = ghno3 + psi1
   
END SUBROUTINE calcn1a
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCQ5
!! *** CASE Q5
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0); SODIUM POOR (SODRAT < 2.0)
!!     2. LIQUID AND SOLID PHASES ARE POSSIBLE
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcq5
   USE Isorropia_Module
   use solut
   implicit none
   
   DOUBLE PRECISION :: nh4i, nai, no3i, nh3aq, no3aq, claq,ohi,ggcl,cli,hi
   double precision :: ggno3,so4i,del,gg,slout,akw
   double precision :: bb,cc,dd,hso4i
   integer :: i
   
!   COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
!       psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
!       a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   frst    =.true.
   calain  =.true.
   calaou  =.true.
   
   ! *** CALCULATE INITIAL SOLUTION ***************************************
   
   CALL calcq1a
   
   psi1   = cna2so4      ! SALTS DISSOLVED
   psi4   = cnh4cl
   psi5   = cnh4no3
   psi6   = cnh42s4
   
   CALL calcmr           ! WATER
   
   nh3aq  = zero
   no3aq  = zero
   claq   = zero
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     akw = xkw*rh*water*water               ! H2O       <==> H+
     
   ! ION CONCENTRATIONS
     
     nai    = waer(1)
     so4i   = waer(2)
     nh4i   = waer(3)
     no3i   = waer(4)
     cli    = waer(5)
     
   ! SOLUTION ACIDIC OR BASIC?
     
     gg   = 2.d0*so4i + no3i + cli - nai - nh4i
     IF (gg > tiny) THEN                        ! H+ in excess
       bb =-gg
       cc =-akw
       dd = bb*bb - 4.d0*cc
       hi = 0.5D0*(-bb + SQRT(dd))
       ohi= akw/hi
     ELSE                                        ! OH- in excess
       bb = gg
       cc =-akw
       dd = bb*bb - 4.d0*cc
       ohi= 0.5D0*(-bb + SQRT(dd))
       hi = akw/ohi
     END IF
     
   ! UNDISSOCIATED SPECIES EQUILIBRIA
     
     IF (hi < ohi) THEN
       CALL calcamaq2 (-gg, nh4i, ohi, nh3aq)
       hi    = akw/ohi
       hso4i = zero
     ELSE
       ggno3 = MAX(2.d0*so4i + no3i - nai - nh4i, zero)
       ggcl  = MAX(gg-ggno3, zero)
       IF (ggcl > tiny) CALL calcclaq2 (ggcl, cli, hi, claq) ! HCl
       IF (ggno3 > tiny) THEN
         IF (ggcl <= tiny) hi = zero
         CALL calcniaq2 (ggno3, no3i, hi, no3aq)              ! HNO3
       END IF
       
   ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
       
       CALL calchs4 (hi, so4i, zero, del)
       so4i  = so4i  - del
       hi    = hi    - del
       hso4i = del
       ohi   = akw/hi
     END IF
     
   ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
     
     molal(1) = hi
     molal(2) = nai
     molal(3) = nh4i
     molal(4) = cli
     molal(5) = so4i
     molal(6) = hso4i
     molal(7) = no3i
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   !cc      CALL PUSHERR (0002, 'CALCQ5')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
   
   20    a2      = (xk2/xkw)*r*temp*(gama(10)/gama(5))**2. ! NH3  <==> NH4+
   a3      = xk4 *r*temp*(water/gama(10))**2.        ! HNO3 <==> NO3-
   a4      = xk3 *r*temp*(water/gama(11))**2.        ! HCL  <==> CL-
   
   gnh3    = nh4i/hi/a2
   ghno3   = hi*no3i/a3
   ghcl    = hi*cli /a4
   
   gasaq(1)= nh3aq
   gasaq(2)= claq
   gasaq(3)= no3aq
   
   cnh42s4 = zero
   cnh4no3 = zero
   cnh4cl  = zero
   cnacl   = zero
   cnano3  = zero
   cna2so4 = zero
   
END SUBROUTINE calcq5

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCQ4
!! *** CASE Q4
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0); SODIUM POOR (SODRAT < 2.0)
!!     2. LIQUID AND SOLID PHASES ARE POSSIBLE
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=====================================================================
SUBROUTINE calcq4
   USE Isorropia_Module
   use solut
   implicit none
   
   LOGICAL :: psconv1
   DOUBLE PRECISION :: nh4i, nai, no3i, nh3aq, no3aq, claq,bb,cc,dd
   double precision :: root3,ohi,ggcl,cli,hi,ggno3,so4i,del,gg
   double precision :: psi1o,hso4i,akw
   integer :: i,islv

!   COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
!       psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
!       a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   frst    =.true.
   calain  =.true.
   calaou  =.true.
   
   psconv1 =.true.
   psi1o   =-great
   root3   = zero
   
   ! *** CALCULATE INITIAL SOLUTION ***************************************
   
   CALL calcq1a
   
   chi1   = cna2so4      ! SALTS
   
   psi1   = cna2so4      ! AMOUNT DISSOLVED
   psi4   = cnh4cl
   psi5   = cnh4no3
   psi6   = cnh42s4
   
   CALL calcmr           ! WATER
   
   nai    = waer(1)      ! LIQUID CONCENTRATIONS
   so4i   = waer(2)
   nh4i   = waer(3)
   no3i   = waer(4)
   cli    = waer(5)
   hso4i  = zero
   nh3aq  = zero
   no3aq  = zero
   claq   = zero
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     a5  = xk5 *(water/gama(2))**3.         ! Na2SO4    <==> Na+
     akw = xkw*rh*water*water               ! H2O       <==> H+
     
   ! SODIUM SULFATE
     
     IF (nai*nai*so4i > a5) THEN
       bb =-(waer(2) + waer(1))
       cc = waer(1)*waer(2) + 0.25*waer(1)*waer(1)
       dd =-0.25*(waer(1)*waer(1)*waer(2) - a5)
       CALL poly3(bb, cc, dd, root3, islv)
       IF (islv /= 0) root3 = tiny
       root3 = MIN (root3, waer(1)/2.0, waer(2), chi1)
       root3 = MAX (root3, zero)
       psi1  = chi1-root3
     END IF
     psconv1 = ABS(psi1-psi1o) <= eps*psi1o
     psi1o   = psi1
     
   ! ION CONCENTRATIONS ; CORRECTIONS
     
     nai = waer(1) - 2.d0*root3
     so4i= waer(2) - root3
     nh4i   = waer(3)
     no3i   = waer(4)
     cli    = waer(5)
     
   ! SOLUTION ACIDIC OR BASIC?
     
     gg   = 2.d0*so4i + no3i + cli - nai - nh4i
     IF (gg > tiny) THEN                        ! H+ in excess
       bb =-gg
       cc =-akw
       dd = bb*bb - 4.d0*cc
       hi = 0.5D0*(-bb + SQRT(dd))
       ohi= akw/hi
     ELSE                                        ! OH- in excess
       bb = gg
       cc =-akw
       dd = bb*bb - 4.d0*cc
       ohi= 0.5D0*(-bb + SQRT(dd))
       hi = akw/ohi
     END IF
     
   ! UNDISSOCIATED SPECIES EQUILIBRIA
     
     IF (hi < ohi) THEN
       CALL calcamaq2 (-gg, nh4i, ohi, nh3aq)
       hi    = akw/ohi
       hso4i = zero
     ELSE
       ggno3 = MAX(2.d0*so4i + no3i - nai - nh4i, zero)
       ggcl  = MAX(gg-ggno3, zero)
       IF (ggcl > tiny) CALL calcclaq2 (ggcl, cli, hi, claq) ! HCl
       IF (ggno3 > tiny) THEN
         IF (ggcl <= tiny) hi = zero
         CALL calcniaq2 (ggno3, no3i, hi, no3aq)              ! HNO3
       END IF
       
   ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
       
       CALL calchs4 (hi, so4i, zero, del)
       so4i  = so4i  - del
       hi    = hi    - del
       hso4i = del
       ohi   = akw/hi
     END IF
     
   ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
     
     molal(1) = hi
     molal(2) = nai
     molal(3) = nh4i
     molal(4) = cli
     molal(5) = so4i
     molal(6) = hso4i
     molal(7) = no3i
     
   ! *** CALCULATE WATER **************************************************
     
     CALL calcmr
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       IF (psconv1) GO TO 20
     END IF
   END DO
   !cc      CALL PUSHERR (0002, 'CALCQ4')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
   
   20    a2      = (xk2/xkw)*r*temp*(gama(10)/gama(5))**2. ! NH3  <==> NH4+
   a3      = xk4 *r*temp*(water/gama(10))**2.        ! HNO3 <==> NO3-
   a4      = xk3 *r*temp*(water/gama(11))**2.        ! HCL  <==> CL-
   
   gnh3    = nh4i/hi/a2
   ghno3   = hi*no3i/a3
   ghcl    = hi*cli /a4
   
   gasaq(1)= nh3aq
   gasaq(2)= claq
   gasaq(3)= no3aq
   
   cnh42s4 = zero
   cnh4no3 = zero
   cnh4cl  = zero
   cnacl   = zero
   cnano3  = zero
   cna2so4 = chi1 - psi1
   
END SUBROUTINE calcq4

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCQ3
!! *** CASE Q3
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!!     2. SOLID & LIQUID AEROSOL POSSIBLE
!!     3. SOLIDS POSSIBLE : NH4CL, NA2SO4, NANO3, NACL
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcq3
   USE Isorropia_Module
   implicit none
   
   LOGICAL :: exno, excl
   EXTERNAL calcq1a, calcq4
   
   ! *** REGIME DEPENDS ON AMBIENT RELATIVE HUMIDITY & POSSIBLE SPECIES ***
   
   exno = waer(4) > tiny
   excl = waer(5) > tiny
   
   IF (exno .OR. excl) THEN             ! *** NITRATE OR CHLORIDE EXISTS
     scase = 'Q3 ; SUBCASE 1'
     CALL calcq3a
     scase = 'Q3 ; SUBCASE 1'
     
   ELSE                                 ! *** NO CHLORIDE AND NITRATE
     IF (rh < drmg3) THEN
       scase = 'Q3 ; SUBCASE 2'
       CALL calcq1a             ! SOLID
       scase = 'Q3 ; SUBCASE 2'
     ELSE
       scase = 'Q3 ; SUBCASE 3' ! MDRH (NH4)2SO4, NA2SO4
       CALL calcmdrp (rh, drmg3, drnh42s4, calcq1a, calcq4)
       scase = 'Q3 ; SUBCASE 3'
     END IF
   END IF
   
END SUBROUTINE calcq3
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCQ3A
!! *** CASE Q3 ; SUBCASE A
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0); SODIUM POOR (SODRAT < 2.0)
!!     2. LIQUID AND SOLID PHASES ARE POSSIBLE
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcq3a
   USE Isorropia_Module
   use solut
   implicit none
   
   LOGICAL :: psconv1, psconv6
   DOUBLE PRECISION :: nh4i, nai, no3i, nh3aq, no3aq, claq, gg
   double precision :: bb,cc,dd,root3,root1,ohi,ggcl,cli,hi,ggno3,so4i,del
   double precision :: psi1o,psi6o,hso4i,akw
   integer :: i,islv
   
!   COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
!       psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
!       a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   frst    =.true.
   calain  =.true.
   calaou  =.true.
   
   psconv1 =.true.
   psconv6 =.true.
   
   psi1o   =-great
   psi6o   =-great
   
   root1   = zero
   root3   = zero
   
   ! *** CALCULATE INITIAL SOLUTION ***************************************
   
   CALL calcq1a
   
   chi1   = cna2so4      ! SALTS
   chi4   = cnh4cl
   chi6   = cnh42s4
   
   psi1   = cna2so4      ! AMOUNT DISSOLVED
   psi4   = cnh4cl
   psi5   = cnh4no3
   psi6   = cnh42s4
   
   CALL calcmr           ! WATER
   
   nai    = waer(1)      ! LIQUID CONCENTRATIONS
   so4i   = waer(2)
   nh4i   = waer(3)
   no3i   = waer(4)
   cli    = waer(5)
   hso4i  = zero
   nh3aq  = zero
   no3aq  = zero
   claq   = zero
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     a5  = xk5 *(water/gama(2))**3.         ! Na2SO4    <==> Na+
     a7  = xk7 *(water/gama(4))**3.         ! (NH4)2SO4 <==> Na+
     akw = xkw*rh*water*water               ! H2O       <==> H+
     
   ! SODIUM SULFATE
     
     IF (nai*nai*so4i > a5) THEN
       bb =-(waer(2) + waer(1) - root1)
       cc = waer(1)*(waer(2) - root1) + 0.25*waer(1)*waer(1)
       dd =-0.25*(waer(1)*waer(1)*(waer(2) - root1) - a5)
       CALL poly3(bb, cc, dd, root3, islv)
       IF (islv /= 0) root3 = tiny
       root3 = MIN (root3, waer(1)/2.0, waer(2) - root1, chi1)
       root3 = MAX (root3, zero)
       psi1  = chi1-root3
     END IF
     psconv1 = ABS(psi1-psi1o) <= eps*psi1o
     psi1o   = psi1
     
   ! AMMONIUM SULFATE
     
     IF (nh4i*nh4i*so4i > a4) THEN
       bb =-(waer(2)+waer(3)-root3)
       cc =  waer(3)*(waer(2)-root3+0.5D0*waer(3))
       dd =-((waer(2)-root3)*waer(3)**2.d0 + a4)/4.d0
       CALL poly3(bb, cc, dd, root1, islv)
       IF (islv /= 0) root1 = tiny
       root1 = MIN(root1, waer(3), waer(2)-root3, chi6)
       root1 = MAX(root1, zero)
       psi6  = chi6-root1
     END IF
     psconv6 = ABS(psi6-psi6o) <= eps*psi6o
     psi6o   = psi6
     
   ! ION CONCENTRATIONS
     
     nai = waer(1) - 2.d0*root3
     so4i= waer(2) - root1 - root3
     nh4i= waer(3) - 2.d0*root1
     no3i= waer(4)
     cli = waer(5)
     
   ! SOLUTION ACIDIC OR BASIC?
     
     gg   = 2.d0*so4i + no3i + cli - nai - nh4i
     IF (gg > tiny) THEN                        ! H+ in excess
       bb =-gg
       cc =-akw
       dd = bb*bb - 4.d0*cc
       hi = 0.5D0*(-bb + SQRT(dd))
       ohi= akw/hi
     ELSE                                        ! OH- in excess
       bb = gg
       cc =-akw
       dd = bb*bb - 4.d0*cc
       ohi= 0.5D0*(-bb + SQRT(dd))
       hi = akw/ohi
     END IF
     
   ! UNDISSOCIATED SPECIES EQUILIBRIA
     
     IF (hi < ohi) THEN
       CALL calcamaq2 (-gg, nh4i, ohi, nh3aq)
       hi    = akw/ohi
       hso4i = zero
     ELSE
       ggno3 = MAX(2.d0*so4i + no3i - nai - nh4i, zero)
       ggcl  = MAX(gg-ggno3, zero)
       IF (ggcl > tiny) CALL calcclaq2 (ggcl, cli, hi, claq) ! HCl
       IF (ggno3 > tiny) THEN
         IF (ggcl <= tiny) hi = zero
         CALL calcniaq2 (ggno3, no3i, hi, no3aq)              ! HNO3
       END IF
       
   ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
       
       CALL calchs4 (hi, so4i, zero, del)
       so4i  = so4i  - del
       hi    = hi    - del
       hso4i = del
       ohi   = akw/hi
     END IF
     
   ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
     
     molal(1) = hi
     molal(2) = nai
     molal(3) = nh4i
     molal(4) = cli
     molal(5) = so4i
     molal(6) = hso4i
     molal(7) = no3i
     
   ! *** CALCULATE WATER **************************************************
     
     CALL calcmr
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       IF (psconv1 .AND. psconv6) GO TO 20
     END IF
   END DO
   !cc      CALL PUSHERR (0002, 'CALCQ3A')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
   
   20    a2      = (xk2/xkw)*r*temp*(gama(10)/gama(5))**2. ! NH3  <==> NH4+
   a3      = xk4 *r*temp*(water/gama(10))**2.        ! HNO3 <==> NO3-
   a4      = xk3 *r*temp*(water/gama(11))**2.        ! HCL  <==> CL-
   
   gnh3    = nh4i/hi/a2
   ghno3   = hi*no3i/a3
   ghcl    = hi*cli /a4
   
   gasaq(1)= nh3aq
   gasaq(2)= claq
   gasaq(3)= no3aq
   
   cnh42s4 = chi6 - psi6
   cnh4no3 = zero
   cnh4cl  = zero
   cnacl   = zero
   cnano3  = zero
   cna2so4 = chi1 - psi1
   
END SUBROUTINE calcq3a

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCQ2
!! *** CASE Q2
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!!     2. SOLID & LIQUID AEROSOL POSSIBLE
!!     3. SOLIDS POSSIBLE : NH4CL, NA2SO4, NANO3, NACL
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcq2
   USE Isorropia_Module
   implicit none
   
   LOGICAL :: exno, excl
   EXTERNAL calcq1a, calcq3a, calcq4
   
   ! *** REGIME DEPENDS ON AMBIENT RELATIVE HUMIDITY & POSSIBLE SPECIES ***
   
   exno = waer(4) > tiny
   excl = waer(5) > tiny
   
   IF (exno) THEN                       ! *** NITRATE EXISTS
     scase = 'Q2 ; SUBCASE 1'
     CALL calcq2a
     scase = 'Q2 ; SUBCASE 1'
     
   ELSE IF (.NOT.exno .AND. excl) THEN   ! *** ONLY CHLORIDE EXISTS
     IF (rh < drmg2) THEN
       scase = 'Q2 ; SUBCASE 2'
       CALL calcq1a             ! SOLID
       scase = 'Q2 ; SUBCASE 2'
     ELSE
       scase = 'Q2 ; SUBCASE 3' ! MDRH (NH4)2SO4, NA2SO4, NH4CL
       CALL calcmdrp (rh, drmg2, drnh4cl, calcq1a, calcq3a)
       scase = 'Q2 ; SUBCASE 3'
     END IF
     
   ELSE                                 ! *** NO CHLORIDE AND NITRATE
     IF (rh < drmg3) THEN
       scase = 'Q2 ; SUBCASE 2'
       CALL calcq1a             ! SOLID
       scase = 'Q2 ; SUBCASE 2'
     ELSE
       scase = 'Q2 ; SUBCASE 4' ! MDRH (NH4)2SO4, NA2SO4
       CALL calcmdrp (rh, drmg3, drnh42s4, calcq1a, calcq4)
       scase = 'Q2 ; SUBCASE 4'
     END IF
   END IF

END SUBROUTINE calcq2
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCQ2A
!! *** CASE Q2 ; SUBCASE A
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0); SODIUM POOR (SODRAT < 2.0)
!!     2. LIQUID AND SOLID PHASES ARE POSSIBLE
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcq2a
   USE Isorropia_Module
   use solut
   implicit none
   
   LOGICAL :: psconv1, psconv4, psconv6
   DOUBLE PRECISION :: nh4i, nai, no3i, nh3aq, no3aq, claq
   double precision :: bb,cc,dd,root3,root1,ohi,ggcl,cli,hi,ggno3,so4i,del,gg
   double precision :: psi1o,psi4o,psi6o,root2,hso4i,a14,akw,root2a,root2b
   integer :: i,islv
   
!   COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
!       psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
!       a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   frst    =.true.
   calain  =.true.
   calaou  =.true.
   
   psconv1 =.true.
   psconv4 =.true.
   psconv6 =.true.
   
   psi1o   =-great
   psi4o   =-great
   psi6o   =-great
   
   root1   = zero
   root2   = zero
   root3   = zero
   
   ! *** CALCULATE INITIAL SOLUTION ***************************************
   
   CALL calcq1a
   
   chi1   = cna2so4      ! SALTS
   chi4   = cnh4cl
   chi6   = cnh42s4
   
   psi1   = cna2so4      ! AMOUNT DISSOLVED
   psi4   = cnh4cl
   psi5   = cnh4no3
   psi6   = cnh42s4
   
   CALL calcmr           ! WATER
   
   nai    = waer(1)      ! LIQUID CONCENTRATIONS
   so4i   = waer(2)
   nh4i   = waer(3)
   no3i   = waer(4)
   cli    = waer(5)
   hso4i  = zero
   nh3aq  = zero
   no3aq  = zero
   claq   = zero
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     a5  = xk5 *(water/gama(2))**3.         ! Na2SO4    <==> Na+
     a14 = xk14*(water/gama(6))**2.         ! NH4Cl     <==> NH4+
     a7  = xk7 *(water/gama(4))**3.         ! (NH4)2SO4 <==> Na+
     akw = xkw*rh*water*water               ! H2O       <==> H+
     
   ! AMMONIUM CHLORIDE
     
     IF (nh4i*cli > a14) THEN
       bb    =-(waer(3) + waer(5) - 2.d0*root1)
       cc    = waer(5)*(waer(3) - 2.d0*root1) - a14
       dd    = bb*bb - 4.d0*cc
       IF (dd < zero) THEN
         root2 = zero
       ELSE
         dd    = SQRT(dd)
         root2a= 0.5D0*(-bb+dd)
         root2b= 0.5D0*(-bb-dd)
         IF (zero <= root2a) THEN
           root2 = root2a
         ELSE
           root2 = root2b
         END IF
         root2 = MIN(root2, waer(5), waer(3) - 2.d0*root1, chi4)
         root2 = MAX(root2, zero)
         psi4  = chi4 - root2
       END IF
     END IF
     psconv4 = ABS(psi4-psi4o) <= eps*psi4o
     psi4o   = psi4
     
   ! SODIUM SULFATE
     
     IF (nai*nai*so4i > a5) THEN
       bb =-(waer(2) + waer(1) - root1)
       cc = waer(1)*(waer(2) - root1) + 0.25*waer(1)*waer(1)
       dd =-0.25*(waer(1)*waer(1)*(waer(2) - root1) - a5)
       CALL poly3(bb, cc, dd, root3, islv)
       IF (islv /= 0) root3 = tiny
       root3 = MIN (root3, waer(1)/2.0, waer(2) - root1, chi1)
       root3 = MAX (root3, zero)
       psi1  = chi1-root3
     END IF
     psconv1 = ABS(psi1-psi1o) <= eps*psi1o
     psi1o   = psi1
     
   ! AMMONIUM SULFATE
     
     IF (nh4i*nh4i*so4i > a4) THEN
       bb =-(waer(2)+waer(3)-root2-root3)
       cc = (waer(3)-root2)*(waer(2)-root3+0.5D0*(waer(3)-root2))
       dd =-((waer(2)-root3)*(waer(3)-root2)**2.d0 + a4)/4.d0
       CALL poly3(bb, cc, dd, root1, islv)
       IF (islv /= 0) root1 = tiny
       root1 = MIN(root1, waer(3)-root2, waer(2)-root3, chi6)
       root1 = MAX(root1, zero)
       psi6  = chi6-root1
     END IF
     psconv6 = ABS(psi6-psi6o) <= eps*psi6o
     psi6o   = psi6
     
   ! ION CONCENTRATIONS
     
     nai = waer(1) - 2.d0*root3
     so4i= waer(2) - root1 - root3
     nh4i= waer(3) - root2 - 2.d0*root1
     no3i= waer(4)
     cli = waer(5) - root2
     
   ! SOLUTION ACIDIC OR BASIC?
     
     gg   = 2.d0*so4i + no3i + cli - nai - nh4i
     IF (gg > tiny) THEN                        ! H+ in excess
       bb =-gg
       cc =-akw
       dd = bb*bb - 4.d0*cc
       hi = 0.5D0*(-bb + SQRT(dd))
       ohi= akw/hi
     ELSE                                        ! OH- in excess
       bb = gg
       cc =-akw
       dd = bb*bb - 4.d0*cc
       ohi= 0.5D0*(-bb + SQRT(dd))
       hi = akw/ohi
     END IF
     
   ! UNDISSOCIATED SPECIES EQUILIBRIA
     
     IF (hi < ohi) THEN
       CALL calcamaq2 (-gg, nh4i, ohi, nh3aq)
       hi    = akw/ohi
       hso4i = zero
     ELSE
       ggno3 = MAX(2.d0*so4i + no3i - nai - nh4i, zero)
       ggcl  = MAX(gg-ggno3, zero)
       IF (ggcl > tiny) CALL calcclaq2 (ggcl, cli, hi, claq) ! HCl
       IF (ggno3 > tiny) THEN
         IF (ggcl <= tiny) hi = zero
         CALL calcniaq2 (ggno3, no3i, hi, no3aq)              ! HNO3
       END IF
       
   ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
       
       CALL calchs4 (hi, so4i, zero, del)
       so4i  = so4i  - del
       hi    = hi    - del
       hso4i = del
       ohi   = akw/hi
     END IF
     
   ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
     
     molal(1) = hi
     molal(2) = nai
     molal(3) = nh4i
     molal(4) = cli
     molal(5) = so4i
     molal(6) = hso4i
     molal(7) = no3i
     
   ! *** CALCULATE WATER **************************************************
     
     CALL calcmr
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       IF (psconv1 .AND. psconv4 .AND. psconv6) GO TO 20
     END IF
   END DO
   !cc      CALL PUSHERR (0002, 'CALCQ2A')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
   
   20    a2      = (xk2/xkw)*r*temp*(gama(10)/gama(5))**2. ! NH3  <==> NH4+
   a3      = xk4 *r*temp*(water/gama(10))**2.        ! HNO3 <==> NO3-
   a4      = xk3 *r*temp*(water/gama(11))**2.        ! HCL  <==> CL-
   
   gnh3    = nh4i/hi/a2
   ghno3   = hi*no3i/a3
   ghcl    = hi*cli /a4
   
   gasaq(1)= nh3aq
   gasaq(2)= claq
   gasaq(3)= no3aq
   
   cnh42s4 = chi6 - psi6
   cnh4no3 = zero
   cnh4cl  = chi4 - psi4
   cnacl   = zero
   cnano3  = zero
   cna2so4 = chi1 - psi1
   
END SUBROUTINE calcq2a

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCQ1
!! *** CASE Q1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!!     2. SOLID AEROSOL ONLY
!!     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, (NH4)2SO4, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcq1
   USE Isorropia_Module
   implicit none
   
   LOGICAL :: exno, excl
   EXTERNAL calcq1a, calcq2a, calcq3a, calcq4
   
   ! *** REGIME DEPENDS ON AMBIENT RELATIVE HUMIDITY & POSSIBLE SPECIES ***
   
   exno = waer(4) > tiny
   excl = waer(5) > tiny
   
   IF (exno .AND. excl) THEN           ! *** NITRATE & CHLORIDE EXIST
     IF (rh < drmg1) THEN
       scase = 'Q1 ; SUBCASE 1'
       CALL calcq1a             ! SOLID
       scase = 'Q1 ; SUBCASE 1'
     ELSE
       scase = 'Q1 ; SUBCASE 2' ! MDRH (NH4)2SO4, NA2SO4, NH4CL, NH4NO3
       CALL calcmdrp (rh, drmg1, drnh4no3, calcq1a, calcq2a)
       scase = 'Q1 ; SUBCASE 2'
     END IF
     
   ELSE IF (exno .AND. .NOT.excl) THEN ! *** ONLY NITRATE EXISTS
     IF (rh < drmq1) THEN
       scase = 'Q1 ; SUBCASE 1'
       CALL calcq1a             ! SOLID
       scase = 'Q1 ; SUBCASE 1'
     ELSE
       scase = 'Q1 ; SUBCASE 3' ! MDRH (NH4)2SO4, NA2SO4, NH4NO3
       CALL calcmdrp (rh, drmq1, drnh4no3, calcq1a, calcq2a)
       scase = 'Q1 ; SUBCASE 3'
     END IF
     
   ELSE IF (.NOT.exno .AND. excl) THEN ! *** ONLY CHLORIDE EXISTS
     IF (rh < drmg2) THEN
       scase = 'Q1 ; SUBCASE 1'
       CALL calcq1a             ! SOLID
       scase = 'Q1 ; SUBCASE 1'
     ELSE
       scase = 'Q1 ; SUBCASE 4' ! MDRH (NH4)2SO4, NA2SO4, NH4CL
       CALL calcmdrp (rh, drmg2, drnh4cl, calcq1a, calcq3a)
       scase = 'Q1 ; SUBCASE 4'
     END IF
     
   ELSE                                ! *** NO CHLORIDE AND NITRATE
     IF (rh < drmg3) THEN
       scase = 'Q1 ; SUBCASE 1'
       CALL calcq1a             ! SOLID
       scase = 'Q1 ; SUBCASE 1'
     ELSE
       scase = 'Q1 ; SUBCASE 5' ! MDRH (NH4)2SO4, NA2SO4
       CALL calcmdrp (rh, drmg3, drnh42s4, calcq1a, calcq4)
       scase = 'Q1 ; SUBCASE 5'
     END IF
   END IF

END SUBROUTINE calcq1
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCQ1A
!! *** CASE Q1 ; SUBCASE 1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
!!     2. SOLID AEROSOL ONLY
!!     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, (NH4)2SO4, NA2SO4
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!======================================================================
SUBROUTINE calcq1a
   USE Isorropia_Module
   implicit none
   
   double precision :: frso4,frnh3
   ! *** CALCULATE SOLIDS **************************************************
   
   cna2so4 = 0.5D0*waer(1)
   frso4   = MAX (waer(2)-cna2so4, zero)
   
   cnh42s4 = MAX (MIN(frso4,0.5D0*waer(3)), tiny)
   frnh3   = MAX (waer(3)-2.d0*cnh42s4, zero)
   
   cnh4no3 = MIN (frnh3, waer(4))
   !CC      FRNO3   = MAX (WAER(4)-CNH4NO3, ZERO)
   frnh3   = MAX (frnh3-cnh4no3, zero)
   
   cnh4cl  = MIN (frnh3, waer(5))
   !CC      FRCL    = MAX (WAER(5)-CNH4CL, ZERO)
   frnh3   = MAX (frnh3-cnh4cl, zero)
   
   ! *** OTHER PHASES ******************************************************
   
   water   = zero
   
   gnh3    = zero
   ghno3   = zero
   ghcl    = zero
   
END SUBROUTINE calcq1a

!!=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCR6
!! *** CASE R6
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0); SODIUM RICH (SODRAT >= 2.0)
!!     2. THERE IS ONLY A LIQUID PHASE
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcr6
   USE Isorropia_Module
   use solut
   implicit none
   
   DOUBLE PRECISION :: nh4i, nai, no3i, nh3aq, no3aq, claq
   double precision :: gg,ohi,ggcl,hi,ggno3,so4i,del,cli
   double precision :: hso4i,akw,bb,cc,dd
   integer :: i
   
!   COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
!       psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
!       a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   CALL calcr1a
   
   psi1   = cna2so4
   psi2   = cnano3
   psi3   = cnacl
   psi4   = cnh4cl
   psi5   = cnh4no3
   
   frst   = .true.
   calain = .true.
   calaou = .true.
   
   ! *** CALCULATE WATER **************************************************
   
   CALL calcmr
   
   ! *** SETUP LIQUID CONCENTRATIONS **************************************
   
   hso4i  = zero
   nh3aq  = zero
   no3aq  = zero
   claq   = zero
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     akw = xkw*rh*water*water                        ! H2O    <==> H+
     
     nai    = waer(1)
     so4i   = waer(2)
     nh4i   = waer(3)
     no3i   = waer(4)
     cli    = waer(5)
     
   ! SOLUTION ACIDIC OR BASIC?
     
     gg  = 2.d0*waer(2) + no3i + cli - nai - nh4i
     IF (gg > tiny) THEN                        ! H+ in excess
       bb =-gg
       cc =-akw
       dd = bb*bb - 4.d0*cc
       hi = 0.5D0*(-bb + SQRT(dd))
       ohi= akw/hi
     ELSE                                        ! OH- in excess
       bb = gg
       cc =-akw
       dd = bb*bb - 4.d0*cc
       ohi= 0.5D0*(-bb + SQRT(dd))
       hi = akw/ohi
     END IF
     
   ! UNDISSOCIATED SPECIES EQUILIBRIA
     
     IF (hi < ohi) THEN
       CALL calcamaq2 (-gg, nh4i, ohi, nh3aq)
       hi    = akw/ohi
     ELSE
       ggno3 = MAX(2.d0*so4i + no3i - nai - nh4i, zero)
       ggcl  = MAX(gg-ggno3, zero)
       IF (ggcl > tiny) CALL calcclaq2 (ggcl, cli, hi, claq) ! HCl
       IF (ggno3 > tiny) THEN
         IF (ggcl <= tiny) hi = zero
         CALL calcniaq2 (ggno3, no3i, hi, no3aq)              ! HNO3
       END IF
       
   ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
       
       CALL calchs4 (hi, so4i, zero, del)
       so4i  = so4i  - del
       hi    = hi    - del
       hso4i = del
       ohi   = akw/hi
     END IF
     
   ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
     
     molal(1) = hi
     molal(2) = nai
     molal(3) = nh4i
     molal(4) = cli
     molal(5) = so4i
     molal(6) = hso4i
     molal(7) = no3i
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       GO TO 20
     END IF
   END DO
   !cc      CALL PUSHERR (0002, 'CALCR6')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
   
   20    a2       = (xk2/xkw)*r*temp*(gama(10)/gama(5))**2. ! NH3  <==> NH4+
   a3       = xk4 *r*temp*(water/gama(10))**2.        ! HNO3 <==> NO3-
   a4       = xk3 *r*temp*(water/gama(11))**2.        ! HCL  <==> CL-
   
   gnh3     = nh4i/hi/a2
   ghno3    = hi*no3i/a3
   ghcl     = hi*cli /a4
   
   gasaq(1) = nh3aq
   gasaq(2) = claq
   gasaq(3) = no3aq
   
   cnh42s4  = zero
   cnh4no3  = zero
   cnh4cl   = zero
   cnacl    = zero
   cnano3   = zero
   cna2so4  = zero
   
END SUBROUTINE calcr6

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCR5
!! *** CASE R5
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0); SODIUM RICH (SODRAT >= 2.0)
!!     2. LIQUID AND SOLID PHASES ARE POSSIBLE
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcr5
   USE Isorropia_Module
   use solut
   implicit none
   
   LOGICAL :: psconv
   DOUBLE PRECISION :: nh4i, nai, no3i, nh3aq, no3aq, claq,gg
   double precision :: bb,cc,dd,root,ohi,ggcl,cli,hi,ggno3,so4i,del
   double precision :: psio,hso4i,akw
   integer :: i,islv
   
!   COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
!       psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
!       a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   LOGICAL :: nean, neac, nesn, nesc
   
   ! *** SETUP PARAMETERS ************************************************
   
   CALL calcr1a                             ! DRY SOLUTION
   
   nean = cnh4no3 <= tiny    ! NH4NO3       ! Water exists?
   neac = cnh4cl <= tiny    ! NH4CL
   nesn = cnano3 <= tiny    ! NANO3
   nesc = cnacl  <= tiny    ! NACL
   IF (nean .AND. neac .AND. nesn .AND. nesc) RETURN
   
   chi1   = cna2so4
   
   psi1   = cna2so4
   psi2   = cnano3
   psi3   = cnacl
   psi4   = cnh4cl
   psi5   = cnh4no3
   
   psio   =-great
   
   ! *** CALCULATE WATER **************************************************
   
   CALL calcmr
   
   frst   = .true.
   calain = .true.
   calaou = .true.
   psconv = .false.
   
   ! *** SETUP LIQUID CONCENTRATIONS **************************************
   
   nai    = waer(1)
   so4i   = waer(2)
   nh4i   = waer(3)
   no3i   = waer(4)
   cli    = waer(5)
   hso4i  = zero
   nh3aq  = zero
   no3aq  = zero
   claq   = zero
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     a5  = xk5*(water/gama(2))**3.                   ! Na2SO4 <==> Na+
     akw = xkw*rh*water*water                        ! H2O    <==> H+
     
   ! SODIUM SULFATE
     
     root = zero
     IF (nai*nai*so4i > a5) THEN
       bb =-3.d0*chi1
       cc = 3.d0*chi1**2.0
       dd =-chi1**3.0 + 0.25D0*a5
       CALL poly3(bb, cc, dd, root, islv)
       IF (islv /= 0) root = tiny
       root = MIN (MAX(root,zero), chi1)
       psi1 = chi1-root
     END IF
     psconv = ABS(psi1-psio) <= eps*psio
     psio   = psi1
     
   ! ION CONCENTRATIONS
     
     nai  = waer(1) - 2.d0*root
     so4i = waer(2) - root
     nh4i = waer(3)
     no3i = waer(4)
     cli  = waer(5)
     
   ! SOLUTION ACIDIC OR BASIC?
     
     gg   = 2.d0*so4i + no3i + cli - nai - nh4i
     IF (gg > tiny) THEN                        ! H+ in excess
       bb =-gg
       cc =-akw
       dd = bb*bb - 4.d0*cc
       hi = 0.5D0*(-bb + SQRT(dd))
       ohi= akw/hi
     ELSE                                        ! OH- in excess
       bb = gg
       cc =-akw
       dd = bb*bb - 4.d0*cc
       ohi= 0.5D0*(-bb + SQRT(dd))
       hi = akw/ohi
     END IF
     
   ! UNDISSOCIATED SPECIES EQUILIBRIA
     
     IF (hi < ohi) THEN
       CALL calcamaq2 (-gg, nh4i, ohi, nh3aq)
       hi    = akw/ohi
     ELSE
       ggno3 = MAX(2.d0*so4i + no3i - nai - nh4i, zero)
       ggcl  = MAX(gg-ggno3, zero)
       IF (ggcl > tiny) CALL calcclaq2 (ggcl, cli, hi, claq) ! HCl
       IF (ggno3 > tiny) THEN
         IF (ggcl <= tiny) hi = zero
         CALL calcniaq2 (ggno3, no3i, hi, no3aq)              ! HNO3
       END IF
       
   ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
       
       CALL calchs4 (hi, so4i, zero, del)
       so4i  = so4i  - del
       hi    = hi    - del
       hso4i = del
       ohi   = akw/hi
     END IF
     
   ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
     
     molal(1) = hi
     molal(2) = nai
     molal(3) = nh4i
     molal(4) = cli
     molal(5) = so4i
     molal(6) = hso4i
     molal(7) = no3i
     
   ! *** CALCULATE WATER **************************************************
     
     CALL calcmr
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       IF (psconv) GO TO 20
     END IF
   END DO
   !cc      CALL PUSHERR (0002, 'CALCR5')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
   
   20    a2       = (xk2/xkw)*r*temp*(gama(10)/gama(5))**2. ! NH3  <==> NH4+
   !C      A21      = XK21*WATER*R*TEMP
   a3       = xk4 *r*temp*(water/gama(10))**2.        ! HNO3 <==> NO3-
   a4       = xk3 *r*temp*(water/gama(11))**2.        ! HCL  <==> CL-
   
   gnh3     = nh4i/hi/a2  ! NH4I*OHI/A2/AKW
   ghno3    = hi*no3i/a3
   ghcl     = hi*cli /a4
   
   gasaq(1) = nh3aq
   gasaq(2) = claq
   gasaq(3) = no3aq
   
   cnh42s4  = zero
   cnh4no3  = zero
   cnh4cl   = zero
   cnacl    = zero
   cnano3   = zero
   cna2so4  = chi1 - psi1
   
END SUBROUTINE calcr5

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCR4
!! *** CASE R4
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!!     2. SOLID AEROSOL ONLY
!!     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NA2SO4, NANO3, NACL
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcr4
   USE Isorropia_Module
   implicit none
   
   LOGICAL :: exan, exac, exsn, exsc
   EXTERNAL calcr1a, calcr5
   
   double precision :: bb,cc,dd,root,ohi,ggcl,cli,ggno3,hi,so4i,del,gg
   
   ! *** SOLVE FOR DRY CASE AND SEE WHICH SOLIDS ARE POSSIBLE **************
   
   scase = 'R4 ; SUBCASE 2'
   CALL calcr1a              ! SOLID
   scase = 'R4 ; SUBCASE 2'
   
   exan = cnh4no3 > tiny    ! NH4NO3
   exac = cnh4cl > tiny    ! NH4CL
   exsn = cnano3 > tiny    ! NANO3
   exsc = cnacl  > tiny    ! NACL
   
   ! *** REGIME DEPENDS ON RELATIVE HUMIDITY AND POSSIBLE SPECIES **********
   
   IF (exan .OR. exsn .OR. exsc) THEN   ! *** NH4NO3,NANO3 EXIST
     IF (rh >= drmh1) THEN
       scase = 'R4 ; SUBCASE 1'
       CALL calcr4a
       scase = 'R4 ; SUBCASE 1'
     END IF
     
   ELSE IF (exac) THEN                  ! *** NH4CL EXISTS ONLY
     IF (rh >= drmr5) THEN
       scase = 'R4 ; SUBCASE 3'
       CALL calcmdrp (rh, drmr5, drnh4cl, calcr1a, calcr5)
       scase = 'R4 ; SUBCASE 3'
     END IF
   END IF
   
END SUBROUTINE calcr4
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCR4A
!! *** CASE R4A
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0); SODIUM RICH (SODRAT >= 2.0)
!!     2. LIQUID AND SOLID PHASES ARE POSSIBLE
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcr4a
   USE Isorropia_Module
   use solut
   implicit none
   
   LOGICAL :: psconv1, psconv4
   DOUBLE PRECISION :: nh4i, nai, no3i, nh3aq, no3aq, claq,gg
   double precision :: bb,cc,dd,root,ohi,ggcl,cli,hi,ggno3,so4i,del
   double precision :: psio1,psio4,hso4i,a14,akw
   integer :: i,islv
   
!   COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
!       psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
!       a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   frst    = .true.
   calain  = .true.
   calaou  = .true.
   psconv1 = .false.
   psconv4 = .false.
   psio1   =-great
   psio4   =-great
   
   ! *** CALCULATE INITIAL SOLUTION ***************************************
   
   CALL calcr1a
   
   chi1   = cna2so4      ! SALTS
   chi4   = cnh4cl
   
   psi1   = cna2so4
   psi2   = cnano3
   psi3   = cnacl
   psi4   = cnh4cl
   psi5   = cnh4no3
   
   CALL calcmr           ! WATER
   
   nai    = waer(1)      ! LIQUID CONCENTRATIONS
   so4i   = waer(2)
   nh4i   = waer(3)
   no3i   = waer(4)
   cli    = waer(5)
   hso4i  = zero
   nh3aq  = zero
   no3aq  = zero
   claq   = zero
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     a5  = xk5 *(water/gama(2))**3.                  ! Na2SO4 <==> Na+
     a14 = xk14*(water/gama(6))**2.                  ! NH4Cl  <==> NH4+
     akw = xkw*rh*water*water                        ! H2O    <==> H+
     
   ! SODIUM SULFATE
     
     root = zero
     IF (nai*nai*so4i > a5) THEN
       bb =-3.d0*chi1
       cc = 3.d0*chi1**2.0
       dd =-chi1**3.0 + 0.25D0*a5
       CALL poly3(bb, cc, dd, root, islv)
       IF (islv /= 0) root = tiny
       root = MIN (MAX(root,zero), chi1)
       psi1 = chi1-root
       nai  = waer(1) - 2.d0*root
       so4i = waer(2) - root
     END IF
     psconv1 = ABS(psi1-psio1) <= eps*psio1
     psio1   = psi1
     
   ! AMMONIUM CHLORIDE
     
     root = zero
     IF (nh4i*cli > a14) THEN
       bb   =-(nh4i + cli)
       cc   =-a14 + nh4i*cli
       dd   = bb*bb - 4.d0*cc
       root = 0.5D0*(-bb-SQRT(dd))
       IF (root > tiny) THEN
         root    = MIN(root, chi4)
         psi4    = chi4 - root
         nh4i    = waer(3) - root
         cli     = waer(5) - root
       END IF
     END IF
     psconv4 = ABS(psi4-psio4) <= eps*psio4
     psio4   = psi4
     
     no3i   = waer(4)
     
   ! SOLUTION ACIDIC OR BASIC?
     
     gg   = 2.d0*so4i + no3i + cli - nai - nh4i
     IF (gg > tiny) THEN                        ! H+ in excess
       bb =-gg
       cc =-akw
       dd = bb*bb - 4.d0*cc
       hi = 0.5D0*(-bb + SQRT(dd))
       ohi= akw/hi
     ELSE                                        ! OH- in excess
       bb = gg
       cc =-akw
       dd = bb*bb - 4.d0*cc
       ohi= 0.5D0*(-bb + SQRT(dd))
       hi = akw/ohi
     END IF
     
   ! UNDISSOCIATED SPECIES EQUILIBRIA
     
     IF (hi < ohi) THEN
       CALL calcamaq2 (-gg, nh4i, ohi, nh3aq)
       hi    = akw/ohi
     ELSE
       ggno3 = MAX(2.d0*so4i + no3i - nai - nh4i, zero)
       ggcl  = MAX(gg-ggno3, zero)
       IF (ggcl > tiny) CALL calcclaq2 (ggcl, cli, hi, claq) ! HCl
       IF (ggno3 > tiny) THEN
         IF (ggcl <= tiny) hi = zero
         CALL calcniaq2 (ggno3, no3i, hi, no3aq)              ! HNO3
       END IF
       
   ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
       
       CALL calchs4 (hi, so4i, zero, del)
       so4i  = so4i  - del
       hi    = hi    - del
       hso4i = del
       ohi   = akw/hi
     END IF
     
   ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
     
     molal(1) = hi
     molal(2) = nai
     molal(3) = nh4i
     molal(4) = cli
     molal(5) = so4i
     molal(6) = hso4i
     molal(7) = no3i
     
   ! *** CALCULATE WATER **************************************************
     
     CALL calcmr
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       IF (psconv1 .AND. psconv4) GO TO 20
     END IF
   END DO
   !cc      CALL PUSHERR (0002, 'CALCR4A')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
   
   20    a2      = (xk2/xkw)*r*temp*(gama(10)/gama(5))**2. ! NH3  <==> NH4+
   a3      = xk4 *r*temp*(water/gama(10))**2.        ! HNO3 <==> NO3-
   a4      = xk3 *r*temp*(water/gama(11))**2.        ! HCL  <==> CL-
   
   gnh3    = nh4i/hi/a2
   ghno3   = hi*no3i/a3
   ghcl    = hi*cli /a4
   
   gasaq(1)= nh3aq
   gasaq(2)= claq
   gasaq(3)= no3aq
   
   cnh42s4 = zero
   cnh4no3 = zero
   cnh4cl  = chi4 - psi4
   cnacl   = zero
   cnano3  = zero
   cna2so4 = chi1 - psi1
   
END SUBROUTINE calcr4a

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCR3
!! *** CASE R3
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!!     2. SOLID AEROSOL ONLY
!!     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NA2SO4, NANO3, NACL
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcr3
   USE Isorropia_Module
   implicit none
   
   LOGICAL :: exan, exac, exsn, exsc
   EXTERNAL calcr1a, calcr4a, calcr5
   
   ! *** SOLVE FOR DRY CASE AND SEE WHICH SOLIDS ARE POSSIBLE **************
   
   scase = 'R3 ; SUBCASE 2'
   CALL calcr1a              ! SOLID
   scase = 'R3 ; SUBCASE 2'
   
   exan = cnh4no3 > tiny    ! NH4NO3
   exac = cnh4cl > tiny    ! NH4CL
   exsn = cnano3 > tiny    ! NANO3
   exsc = cnacl  > tiny    ! NACL
   
   ! *** REGIME DEPENDS ON RELATIVE HUMIDITY AND POSSIBLE SPECIES **********
   
   IF (exan .OR. exsn) THEN                   ! *** NH4NO3,NANO3 EXIST
     IF (rh >= drmh1) THEN
       scase = 'R3 ; SUBCASE 1'
       CALL calcr3a
       scase = 'R3 ; SUBCASE 1'
     END IF
     
   ELSE IF (.NOT.exan .AND. .NOT.exsn) THEN   ! *** NH4NO3,NANO3 = 0
     IF      (     exac .AND.      exsc) THEN
       IF (rh >= drmr4) THEN
         scase = 'R3 ; SUBCASE 3'
         CALL calcmdrp (rh, drmr4, drnacl, calcr1a, calcr4a)
         scase = 'R3 ; SUBCASE 3'
       END IF
       
     ELSE IF (.NOT.exac .AND.      exsc) THEN
       IF (rh >= drmr2) THEN
         scase = 'R3 ; SUBCASE 4'
         CALL calcmdrp (rh, drmr2, drnacl, calcr1a, calcr4a)
         scase = 'R3 ; SUBCASE 4'
       END IF
       
     ELSE IF (     exac .AND. .NOT.exsc) THEN
       IF (rh >= drmr5) THEN
         scase = 'R3 ; SUBCASE 5'
         CALL calcmdrp (rh, drmr5, drnacl, calcr1a, calcr5)
         scase = 'R3 ; SUBCASE 5'
       END IF
     END IF
     
   END IF
   
END SUBROUTINE calcr3
   
!!=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCR3A
!! *** CASE R3A
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0); SODIUM RICH (SODRAT >= 2.0)
!!     2. LIQUID AND SOLID PHASES ARE POSSIBLE
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcr3a
   USE Isorropia_Module
   use solut
   implicit none
   
   LOGICAL :: psconv1, psconv3, psconv4
   DOUBLE PRECISION :: nh4i, nai, no3i, nh3aq, no3aq, claq,gg
   double precision :: bb,cc,dd,root1,ohi,ggcl,cli,hi,ggno3,so4i,del
   double precision :: psi1o,psi3o,psi4o,root2,root3,hso4i,a14,akw
   double precision :: root2a,root2b,root3a,root3b
   integer :: i,islv
   
!   COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
!       psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
!       a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   frst    =.true.
   calain  =.true.
   calaou  =.true.
   psconv1 =.true.
   psconv3 =.true.
   psconv4 =.true.
   psi1o   =-great
   psi3o   =-great
   psi4o   =-great
   root1   = zero
   root2   = zero
   root3   = zero
   
   ! *** CALCULATE INITIAL SOLUTION ***************************************
   
   CALL calcr1a
   
   chi1   = cna2so4      ! SALTS
   chi4   = cnh4cl
   chi3   = cnacl
   
   psi1   = cna2so4
   psi2   = cnano3
   psi3   = cnacl
   psi4   = cnh4cl
   psi5   = cnh4no3
   
   CALL calcmr           ! WATER
   
   nai    = waer(1)      ! LIQUID CONCENTRATIONS
   so4i   = waer(2)
   nh4i   = waer(3)
   no3i   = waer(4)
   cli    = waer(5)
   hso4i  = zero
   nh3aq  = zero
   no3aq  = zero
   claq   = zero
   
   molal(1) = zero
   molal(2) = nai
   molal(3) = nh4i
   molal(4) = cli
   molal(5) = so4i
   molal(6) = hso4i
   molal(7) = no3i
   
   CALL calcact          ! CALCULATE ACTIVITY COEFFICIENTS
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     a5  = xk5 *(water/gama(2))**3.                  ! Na2SO4 <==> Na+
     a8  = xk8 *(water/gama(1))**2.                  ! NaCl   <==> Na+
     a14 = xk14*(water/gama(6))**2.                  ! NH4Cl  <==> NH4+
     akw = xkw*rh*water*water                        ! H2O    <==> H+
     
   ! AMMONIUM CHLORIDE
     
     IF (nh4i*cli > a14) THEN
       bb    =-(waer(3) + waer(5) - root3)
       cc    =-a14 + nh4i*(waer(5) - root3)
       dd    = MAX(bb*bb - 4.d0*cc, zero)
       root2a= 0.5D0*(-bb+SQRT(dd))
       root2b= 0.5D0*(-bb-SQRT(dd))
       IF (zero <= root2a) THEN
         root2 = root2a
       ELSE
         root2 = root2b
       END IF
       root2 = MIN(MAX(zero, root2), MAX(waer(5)-root3,zero), chi4, waer(3))
       psi4  = chi4 - root2
     END IF
     psconv4 = ABS(psi4-psi4o) <= eps*psi4o
     psi4o   = psi4
     
   ! SODIUM SULFATE
     
     IF (nai*nai*so4i > a5) THEN
       bb =-(chi1 + waer(1) - root3)
       cc = 0.25D0*(waer(1) - root3)*(4.d0*chi1+waer(1)-root3)
       dd =-0.25D0*(chi1*(waer(1)-root3)**2.d0 - a5)
       CALL poly3(bb, cc, dd, root1, islv)
       IF (islv /= 0) root1 = tiny
       root1 = MIN (MAX(root1,zero), MAX(waer(1)-root3,zero), chi1, waer(2))
       psi1  = chi1-root1
     END IF
     psconv1 = ABS(psi1-psi1o) <= eps*psi1o
     psi1o   = psi1
     
   ! ION CONCENTRATIONS
     
     nai = waer(1) - (2.d0*root1 + root3)
     so4i= waer(2) - root1
     nh4i= waer(3) - root2
     cli = waer(5) - (root3 + root2)
     no3i= waer(4)
     
   ! SODIUM CHLORIDE  ; To obtain new value for ROOT3
     
     IF (nai*cli > a8) THEN
       bb    =-((chi1-2.d0*root1) + (waer(5) - root2))
       cc    = (chi1-2.d0*root1)*(waer(5) - root2) - a8
       dd    = SQRT(MAX(bb*bb - 4.d0*cc, tiny))
       root3a= 0.5D0*(-bb-SQRT(dd))
       root3b= 0.5D0*(-bb+SQRT(dd))
       IF (zero <= root3a) THEN
         root3 = root3a
       ELSE
         root3 = root3b
       END IF
       root3   = MIN(MAX(root3, zero), chi3)
       psi3    = chi3-root3
     END IF
     psconv3 = ABS(psi3-psi3o) <= eps*psi3o
     psi3o   = psi3
     
   ! SOLUTION ACIDIC OR BASIC?
     
     gg   = 2.d0*so4i + no3i + cli - nai - nh4i
     IF (gg > tiny) THEN                        ! H+ in excess
       bb =-gg
       cc =-akw
       dd = bb*bb - 4.d0*cc
       hi = 0.5D0*(-bb + SQRT(dd))
       ohi= akw/hi
     ELSE                                        ! OH- in excess
       bb = gg
       cc =-akw
       dd = bb*bb - 4.d0*cc
       ohi= 0.5D0*(-bb + SQRT(dd))
       hi = akw/ohi
     END IF
     
   ! UNDISSOCIATED SPECIES EQUILIBRIA
     
     IF (hi < ohi) THEN
       CALL calcamaq2 (-gg, nh4i, ohi, nh3aq)
       hi    = akw/ohi
     ELSE
       ggno3 = MAX(2.d0*so4i + no3i - nai - nh4i, zero)
       ggcl  = MAX(gg-ggno3, zero)
       IF (ggcl > tiny) CALL calcclaq2 (ggcl, cli, hi, claq) ! HCl
       IF (ggno3 > tiny) THEN
         IF (ggcl <= tiny) hi = zero
         CALL calcniaq2 (ggno3, no3i, hi, no3aq)              ! HNO3
       END IF
       
   ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
       
       CALL calchs4 (hi, so4i, zero, del)
       so4i  = so4i  - del
       hi    = hi    - del
       hso4i = del
       ohi   = akw/hi
     END IF
     
   ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
     
     molal(1) = hi
     molal(2) = nai
     molal(3) = nh4i
     molal(4) = cli
     molal(5) = so4i
     molal(6) = hso4i
     molal(7) = no3i
     
   ! *** CALCULATE WATER **************************************************
     
     CALL calcmr
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       IF (psconv1.AND.psconv3.AND.psconv4) GO TO 20
     END IF
   END DO
   !cc      CALL PUSHERR (0002, 'CALCR3A')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
   
   20    IF (cli <= tiny .AND. waer(5) > tiny) THEN !No disslv Cl-;solid only
     DO  i=1,nions
       molal(i) = zero
     END DO
     DO  i=1,ngasaq
       gasaq(i) = zero
     END DO
     CALL calcr1a
   ELSE
     a2      = (xk2/xkw)*r*temp*(gama(10)/gama(5))**2. ! NH3  <==> NH4+
     a3      = xk4 *r*temp*(water/gama(10))**2.        ! HNO3 <==> NO3-
     a4      = xk3 *r*temp*(water/gama(11))**2.        ! HCL  <==> CL-
     
     gnh3    = nh4i/hi/a2
     ghno3   = hi*no3i/a3
     ghcl    = hi*cli /a4
     
     gasaq(1)= nh3aq
     gasaq(2)= claq
     gasaq(3)= no3aq
     
     cnh42s4 = zero
     cnh4no3 = zero
     cnh4cl  = chi4 - psi4
     cnacl   = chi3 - psi3
     cnano3  = zero
     cna2so4 = chi1 - psi1
   END IF
   
END SUBROUTINE calcr3a

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCR2
!! *** CASE R2
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!!     2. SOLID AEROSOL ONLY
!!     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NA2SO4, NANO3, NACL
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcr2
   USE Isorropia_Module
   implicit none
   
   LOGICAL :: exan, exac, exsn, exsc
   EXTERNAL calcr1a, calcr3a, calcr4a, calcr5
   
   ! *** SOLVE FOR DRY CASE AND SEE WHICH SOLIDS ARE POSSIBLE **************
   
   scase = 'R2 ; SUBCASE 2'
   CALL calcr1a              ! SOLID
   scase = 'R2 ; SUBCASE 2'
   
   exan = cnh4no3 > tiny    ! NH4NO3
   exac = cnh4cl > tiny    ! NH4CL
   exsn = cnano3 > tiny    ! NANO3
   exsc = cnacl  > tiny    ! NACL
   
   ! *** REGIME DEPENDS ON RELATIVE HUMIDITY AND POSSIBLE SPECIES **********
   
   IF (exan) THEN                             ! *** NH4NO3 EXISTS
     IF (rh >= drmh1) THEN
       scase = 'R2 ; SUBCASE 1'
       CALL calcr2a
       scase = 'R2 ; SUBCASE 1'
     END IF
     
   ELSE IF (.NOT.exan) THEN                   ! *** NH4NO3 = 0
     IF      (     exac .AND.      exsn .AND.      exsc) THEN
       IF (rh >= drmh2) THEN
         scase = 'R2 ; SUBCASE 3'
         CALL calcmdrp (rh, drmh2, drnano3, calcr1a, calcr3a)
         scase = 'R2 ; SUBCASE 3'
       END IF
       
     ELSE IF (.NOT.exac .AND.      exsn .AND.      exsc) THEN
       IF (rh >= drmr1) THEN
         scase = 'R2 ; SUBCASE 4'
         CALL calcmdrp (rh, drmr1, drnano3, calcr1a, calcr3a)
         scase = 'R2 ; SUBCASE 4'
       END IF
       
     ELSE IF (.NOT.exac .AND. .NOT.exsn .AND.      exsc) THEN
       IF (rh >= drmr2) THEN
         scase = 'R2 ; SUBCASE 5'
         CALL calcmdrp (rh, drmr2, drnacl, calcr1a, calcr4a)
         scase = 'R2 ; SUBCASE 5'
       END IF
       
     ELSE IF (.NOT.exac .AND.      exsn .AND. .NOT.exsc) THEN
       IF (rh >= drmr3) THEN
         scase = 'R2 ; SUBCASE 6'
         CALL calcmdrp (rh, drmr3, drnano3, calcr1a, calcr3a)
         scase = 'R2 ; SUBCASE 6'
       END IF
       
     ELSE IF (     exac .AND. .NOT.exsn .AND.      exsc) THEN
       IF (rh >= drmr4) THEN
         scase = 'R2 ; SUBCASE 7'
         CALL calcmdrp (rh, drmr4, drnacl, calcr1a, calcr4a)
         scase = 'R2 ; SUBCASE 7'
       END IF
       
     ELSE IF (     exac .AND. .NOT.exsn .AND. .NOT.exsc) THEN
       IF (rh >= drmr5) THEN
         scase = 'R2 ; SUBCASE 8'
         CALL calcmdrp (rh, drmr5, drnh4cl, calcr1a, calcr5)
         scase = 'R2 ; SUBCASE 8'
       END IF
       
     ELSE IF (     exac .AND.      exsn .AND. .NOT.exsc) THEN
       IF (rh >= drmr6) THEN
         scase = 'R2 ; SUBCASE 9'
         CALL calcmdrp (rh, drmr6, drnano3, calcr1a, calcr3a)
         scase = 'R2 ; SUBCASE 9'
       END IF
     END IF
     
   END IF
   
END SUBROUTINE calcr2
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCR2A
!! *** CASE R2A
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0); SODIUM RICH (SODRAT >= 2.0)
!!     2. LIQUID AND SOLID PHASES ARE POSSIBLE
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcr2a
   USE Isorropia_Module
   use solut
   implicit none
   
   LOGICAL :: psconv1, psconv2, psconv3, psconv4
   DOUBLE PRECISION :: nh4i, nai, no3i, nh3aq, no3aq, claq
   double precision :: bb,cc,dd,gg,root1,ohi,ggcl,cli,hi,ggno3,so4i,del
   double precision :: psi1o,psi2o,psi3o,psi4o,root2,root3,root4,hso4i
   double precision :: a9,a14,akw,root2a,root2b,root3a,root3b,root4a,root4b
   integer :: i,islv
   
!   COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8,  &
!       psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8,  &
!       a1,   a2,   a3,   a4,   a5,   a6,   a7,   a8
   
   ! *** SETUP PARAMETERS ************************************************
   
   frst    =.true.
   calain  =.true.
   calaou  =.true.
   
   psconv1 =.true.
   psconv2 =.true.
   psconv3 =.true.
   psconv4 =.true.
   
   psi1o   =-great
   psi2o   =-great
   psi3o   =-great
   psi4o   =-great
   
   root1   = zero
   root2   = zero
   root3   = zero
   root4   = zero
   
   ! *** CALCULATE INITIAL SOLUTION ***************************************
   
   CALL calcr1a
   
   chi1   = cna2so4      ! SALTS
   chi2   = cnano3
   chi3   = cnacl
   chi4   = cnh4cl
   
   psi1   = cna2so4
   psi2   = cnano3
   psi3   = cnacl
   psi4   = cnh4cl
   psi5   = cnh4no3
   
   CALL calcmr           ! WATER
   
   nai    = waer(1)      ! LIQUID CONCENTRATIONS
   so4i   = waer(2)
   nh4i   = waer(3)
   no3i   = waer(4)
   cli    = waer(5)
   hso4i  = zero
   nh3aq  = zero
   no3aq  = zero
   claq   = zero
   
   molal(1) = zero
   molal(2) = nai
   molal(3) = nh4i
   molal(4) = cli
   molal(5) = so4i
   molal(6) = hso4i
   molal(7) = no3i
   
   CALL calcact          ! CALCULATE ACTIVITY COEFFICIENTS
   
   ! *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
   
   DO  i=1,nsweep
     a5  = xk5 *(water/gama(2))**3.                  ! Na2SO4 <==> Na+
     a8  = xk8 *(water/gama(1))**2.                  ! NaCl   <==> Na+
     a9  = xk9 *(water/gama(3))**2.                  ! NaNO3  <==> Na+
     a14 = xk14*(water/gama(6))**2.                  ! NH4Cl  <==> NH4+
     akw = xkw*rh*water*water                        ! H2O    <==> H+
     
   ! AMMONIUM CHLORIDE
     
     IF (nh4i*cli > a14) THEN
       bb    =-(waer(3) + waer(5) - root3)
       cc    = nh4i*(waer(5) - root3) - a14
       dd    = MAX(bb*bb - 4.d0*cc, zero)
       dd    = SQRT(dd)
       root2a= 0.5D0*(-bb+dd)
       root2b= 0.5D0*(-bb-dd)
       IF (zero <= root2a) THEN
         root2 = root2a
       ELSE
         root2 = root2b
       END IF
       root2 = MIN(MAX(root2, zero), chi4)
       psi4  = chi4 - root2
     END IF
     psconv4 = ABS(psi4-psi4o) <= eps*psi4o
     psi4o   = psi4
     
   ! SODIUM SULFATE
     
     IF (nai*nai*so4i > a5) THEN
       bb =-(waer(2) + waer(1) - root3 - root4)
       cc = waer(1)*(2.d0*root3 + 2.d0*root4 - 4.d0*waer(2) - one)  &
           -(root3 + root4)**2.0 + 4.d0*waer(2)*(root3 + root4)
       cc =-0.25*cc
       dd = waer(1)*waer(2)*(one - 2.d0*root3 - 2.d0*root4) +  &
           waer(2)*(root3 + root4)**2.0 - a5
       dd =-0.25*dd
       CALL poly3(bb, cc, dd, root1, islv)
       IF (islv /= 0) root1 = tiny
       root1 = MIN (MAX(root1,zero), chi1)
       psi1  = chi1-root1
     END IF
     psconv1 = ABS(psi1-psi1o) <= eps*psi1o
     psi1o   = psi1
     
   ! SODIUM NITRATE
     
     IF (nai*no3i > a9) THEN
       bb    =-(waer(4) + waer(1) - 2.d0*root1 - root3)
       cc    = waer(4)*(waer(1) - 2.d0*root1 - root3) - a9
       dd    = SQRT(MAX(bb*bb - 4.d0*cc, tiny))
       root4a= 0.5D0*(-bb-dd)
       root4b= 0.5D0*(-bb+dd)
       IF (zero <= root4a) THEN
         root4 = root4a
       ELSE
         root4 = root4b
       END IF
       root4 = MIN(MAX(root4, zero), chi2)
       psi2  = chi2-root4
     END IF
     psconv2 = ABS(psi2-psi2o) <= eps*psi2o
     psi2o   = psi2
     
   ! ION CONCENTRATIONS
     
     nai = waer(1) - (2.d0*root1 + root3 + root4)
     so4i= waer(2) - root1
     nh4i= waer(3) - root2
     no3i= waer(4) - root4
     cli = waer(5) - (root3 + root2)
     
   ! SODIUM CHLORIDE  ; To obtain new value for ROOT3
     
     IF (nai*cli > a8) THEN
       bb    =-(waer(1) - 2.d0*root1 + waer(5) - root2 - root4)
       cc    = (waer(5) + root2)*(waer(1) - 2.d0*root1 - root4) - a8
       dd    = SQRT(MAX(bb*bb - 4.d0*cc, tiny))
       root3a= 0.5D0*(-bb-dd)
       root3b= 0.5D0*(-bb+dd)
       IF (zero <= root3a) THEN
         root3 = root3a
       ELSE
         root3 = root3b
       END IF
       root3   = MIN(MAX(root3, zero), chi3)
       psi3    = chi3-root3
     END IF
     psconv3 = ABS(psi3-psi3o) <= eps*psi3o
     psi3o   = psi3
     
   ! SOLUTION ACIDIC OR BASIC?
     
     gg   = 2.d0*so4i + no3i + cli - nai - nh4i
     IF (gg > tiny) THEN                        ! H+ in excess
       bb =-gg
       cc =-akw
       dd = bb*bb - 4.d0*cc
       hi = 0.5D0*(-bb + SQRT(dd))
       ohi= akw/hi
     ELSE                                        ! OH- in excess
       bb = gg
       cc =-akw
       dd = bb*bb - 4.d0*cc
       ohi= 0.5D0*(-bb + SQRT(dd))
       hi = akw/ohi
     END IF
     
   ! UNDISSOCIATED SPECIES EQUILIBRIA
     
     IF (hi < ohi) THEN
       CALL calcamaq2 (-gg, nh4i, ohi, nh3aq)
       hi    = akw/ohi
     ELSE
       ggno3 = MAX(2.d0*so4i + no3i - nai - nh4i, zero)
       ggcl  = MAX(gg-ggno3, zero)
       IF (ggcl > tiny) CALL calcclaq2 (ggcl, cli, hi, claq) ! HCl
       IF (ggno3 > tiny) THEN
         IF (ggcl <= tiny) hi = zero
         CALL calcniaq2 (ggno3, no3i, hi, no3aq)              ! HNO3
       END IF
       
   ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
       
       CALL calchs4 (hi, so4i, zero, del)
       so4i  = so4i  - del
       hi    = hi    - del
       hso4i = del
       ohi   = akw/hi
     END IF
     
   ! *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
     
     molal(1) = hi
     molal(2) = nai
     molal(3) = nh4i
     molal(4) = cli
     molal(5) = so4i
     molal(6) = hso4i
     molal(7) = no3i
     
   ! *** CALCULATE WATER **************************************************
     
     CALL calcmr
     
   ! *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
     
     IF (frst.AND.calaou .OR. .NOT.frst.AND.calain) THEN
       CALL calcact
     ELSE
       IF (psconv1.AND.psconv2.AND.psconv3.AND.psconv4) GO TO 20
     END IF
   END DO
   !cc      CALL PUSHERR (0002, 'CALCR2A')    ! WARNING ERROR: NO CONVERGENCE
   
   ! *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
   
   20    IF (cli <= tiny .AND. waer(5) > tiny) THEN !No disslv Cl-;solid only
     DO  i=1,nions
       molal(i) = zero
     END DO
     DO  i=1,ngasaq
       gasaq(i) = zero
     END DO
     CALL calcr1a
   ELSE                                     ! OK, aqueous phase present
     a2      = (xk2/xkw)*r*temp*(gama(10)/gama(5))**2. ! NH3  <==> NH4+
     a3      = xk4 *r*temp*(water/gama(10))**2.        ! HNO3 <==> NO3-
     a4      = xk3 *r*temp*(water/gama(11))**2.        ! HCL  <==> CL-
     
     gnh3    = nh4i/hi/a2
     ghno3   = hi*no3i/a3
     ghcl    = hi*cli /a4
     
     gasaq(1)= nh3aq
     gasaq(2)= claq
     gasaq(3)= no3aq
     
     cnh42s4 = zero
     cnh4no3 = zero
     cnh4cl  = chi4 - psi4
     cnacl   = chi3 - psi3
     cnano3  = chi2 - psi2
     cna2so4 = chi1 - psi1
   END IF
   
END SUBROUTINE calcr2a

!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCR1
!! *** CASE R1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!!     2. SOLID AEROSOL ONLY
!!     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NA2SO4, NANO3, NACL
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcr1
   USE Isorropia_Module
   implicit none
   
   LOGICAL :: exan, exac, exsn, exsc
   EXTERNAL calcr1a, calcr2a, calcr3a, calcr4a, calcr5
   
   ! *** SOLVE FOR DRY CASE AND SEE WHICH SOLIDS ARE POSSIBLE **************
   
   scase = 'R1 ; SUBCASE 1'
   CALL calcr1a              ! SOLID
   scase = 'R1 ; SUBCASE 1'
   
   exan = cnh4no3 > tiny    ! NH4NO3
   exac = cnh4cl > tiny    ! NH4CL
   exsn = cnano3 > tiny    ! NANO3
   exsc = cnacl  > tiny    ! NACL
   
   ! *** REGIME DEPENDS ON RELATIVE HUMIDITY AND POSSIBLE SPECIES **********
   
   IF (exan.AND.exac.AND.exsc.AND.exsn) THEN  ! *** ALL EXIST
     IF (rh >= drmh1) THEN
       scase = 'R1 ; SUBCASE 2'  ! MDRH
       CALL calcmdrp (rh, drmh1, drnh4no3, calcr1a, calcr2a)
       scase = 'R1 ; SUBCASE 2'
     END IF
     
   ELSE IF (.NOT.exan) THEN                   ! *** NH4NO3 = 0
     IF      (     exac .AND.      exsn .AND.      exsc) THEN
       IF (rh >= drmh2) THEN
         scase = 'R1 ; SUBCASE 3'
         CALL calcmdrp (rh, drmh2, drnano3, calcr1a, calcr3a)
         scase = 'R1 ; SUBCASE 3'
       END IF
       
     ELSE IF (.NOT.exac .AND.      exsn .AND.      exsc) THEN
       IF (rh >= drmr1) THEN
         scase = 'R1 ; SUBCASE 4'
         CALL calcmdrp (rh, drmr1, drnano3, calcr1a, calcr3a)
         scase = 'R1 ; SUBCASE 4'
       END IF
       
     ELSE IF (.NOT.exac .AND. .NOT.exsn .AND.      exsc) THEN
       IF (rh >= drmr2) THEN
         scase = 'R1 ; SUBCASE 5'
         CALL calcmdrp (rh, drmr2, drnacl, calcr1a, calcr3a) !, CALCR4A)
         scase = 'R1 ; SUBCASE 5'
       END IF
       
     ELSE IF (.NOT.exac .AND.      exsn .AND. .NOT.exsc) THEN
       IF (rh >= drmr3) THEN
         scase = 'R1 ; SUBCASE 6'
         CALL calcmdrp (rh, drmr3, drnano3, calcr1a, calcr3a)
         scase = 'R1 ; SUBCASE 6'
       END IF
       
     ELSE IF (     exac .AND. .NOT.exsn .AND.      exsc) THEN
       IF (rh >= drmr4) THEN
         scase = 'R1 ; SUBCASE 7'
         CALL calcmdrp (rh, drmr4, drnacl, calcr1a, calcr3a) !, CALCR4A)
         scase = 'R1 ; SUBCASE 7'
       END IF
       
     ELSE IF (     exac .AND. .NOT.exsn .AND. .NOT.exsc) THEN
       IF (rh >= drmr5) THEN
         scase = 'R1 ; SUBCASE 8'
         CALL calcmdrp (rh, drmr5, drnh4cl, calcr1a, calcr3a) !, CALCR5)
         scase = 'R1 ; SUBCASE 8'
       END IF
       
     ELSE IF (     exac .AND.      exsn .AND. .NOT.exsc) THEN
       IF (rh >= drmr6) THEN
         scase = 'R1 ; SUBCASE 9'
         CALL calcmdrp (rh, drmr6, drnano3, calcr1a, calcr3a)
         scase = 'R1 ; SUBCASE 9'
       END IF
     END IF
     
   ELSE IF (.NOT.exac) THEN                   ! *** NH4CL  = 0
     IF      (     exan .AND.      exsn .AND.      exsc) THEN
       IF (rh >= drmr7) THEN
         scase = 'R1 ; SUBCASE 10'
         CALL calcmdrp (rh, drmr7, drnh4no3, calcr1a, calcr2a)
         scase = 'R1 ; SUBCASE 10'
       END IF
       
     ELSE IF (     exan .AND. .NOT.exsn .AND.      exsc) THEN
       IF (rh >= drmr8) THEN
         scase = 'R1 ; SUBCASE 11'
         CALL calcmdrp (rh, drmr8, drnh4no3, calcr1a, calcr2a)
         scase = 'R1 ; SUBCASE 11'
       END IF
       
     ELSE IF (     exan .AND. .NOT.exsn .AND. .NOT.exsc) THEN
       IF (rh >= drmr9) THEN
         scase = 'R1 ; SUBCASE 12'
         CALL calcmdrp (rh, drmr9, drnh4no3, calcr1a, calcr2a)
         scase = 'R1 ; SUBCASE 12'
       END IF
       
     ELSE IF (     exan .AND.      exsn .AND. .NOT.exsc) THEN
       IF (rh >= drmr10) THEN
         scase = 'R1 ; SUBCASE 13'
         CALL calcmdrp (rh, drmr10, drnh4no3, calcr1a, calcr2a)
         scase = 'R1 ; SUBCASE 13'
       END IF
     END IF
     
   ELSE IF (.NOT.exsn) THEN                  ! *** NANO3  = 0
     IF      (     exan .AND.      exac .AND.      exsc) THEN
       IF (rh >= drmr11) THEN
         scase = 'R1 ; SUBCASE 14'
         CALL calcmdrp (rh, drmr11, drnh4no3, calcr1a, calcr2a)
         scase = 'R1 ; SUBCASE 14'
       END IF
       
     ELSE IF (     exan .AND.      exac .AND. .NOT.exsc) THEN
       IF (rh >= drmr12) THEN
         scase = 'R1 ; SUBCASE 15'
         CALL calcmdrp (rh, drmr12, drnh4no3, calcr1a, calcr2a)
         scase = 'R1 ; SUBCASE 15'
       END IF
     END IF
     
   ELSE IF (.NOT.exsc) THEN                  ! *** NACL   = 0
     IF      (     exan .AND.      exac .AND.      exsn) THEN
       IF (rh >= drmr13) THEN
         scase = 'R1 ; SUBCASE 16'
         CALL calcmdrp (rh, drmr13, drnh4no3, calcr1a, calcr2a)
         scase = 'R1 ; SUBCASE 16'
       END IF
     END IF
   END IF
   
END SUBROUTINE calcr1
   
!>=======================================================================
!! *** ISORROPIA CODE
!! *** SUBROUTINE CALCR1A
!! *** CASE R1 ; SUBCASE 1
!!     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
!!     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
!!     2. SOLID AEROSOL ONLY
!!     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NANO3, NA2SO4, NANO3, NACL
!! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
!!!! @author Athanasios Nenes
!!=======================================================================
SUBROUTINE calcr1a
   USE Isorropia_Module
   implicit none
   
   double precision :: frna,frno3,frcl,frnh3
   ! *** CALCULATE SOLIDS **************************************************
   
   cna2so4 = waer(2)
   frna    = MAX (waer(1)-2*cna2so4, zero)
   
   cnh42s4 = zero
   
   cnano3  = MIN (frna, waer(4))
   frno3   = MAX (waer(4)-cnano3, zero)
   frna    = MAX (frna-cnano3, zero)
   
   cnacl   = MIN (frna, waer(5))
   frcl    = MAX (waer(5)-cnacl, zero)
   frna    = MAX (frna-cnacl, zero)
   
   cnh4no3 = MIN (frno3, waer(3))
   frno3   = MAX (frno3-cnh4no3, zero)
   frnh3   = MAX (waer(3)-cnh4no3, zero)
   
   cnh4cl  = MIN (frcl, frnh3)
   frcl    = MAX (frcl-cnh4cl, zero)
   frnh3   = MAX (frnh3-cnh4cl, zero)
   
   ! *** OTHER PHASES ******************************************************
   
   water   = zero
   
   gnh3    = zero
   ghno3   = zero
   ghcl    = zero
   
END SUBROUTINE calcr1a
   
   