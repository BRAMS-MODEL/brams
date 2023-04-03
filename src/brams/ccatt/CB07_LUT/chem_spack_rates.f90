  
 MODULE mod_chem_spack_rates
  
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: rates ! subroutine
 CONTAINS
  
   SUBROUTINE rates(rk,y,w,ngas,ijkbeg,ijkend,maxblock_size,nr)
 
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This routine computes the reaction rates.
!     This routine is automatically generated by SPACK.
!     Mechanism: ../Mechanism/CB07   
!     Species: ../Mechanism/ciCB07 
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     RK: kinetic rates.
!     Y: chemical concentrations.
!
!     -- INPUT/OUTPUT VARIABLES
!
!     -- OUTPUT VARIABLES
!
!     W: reaction rates.
!
!------------------------------------------------------------------------
!
!     -- REMARKS
!
!------------------------------------------------------------------------
!
!     -- MODIFICATIONS
!
!------------------------------------------------------------------------
!
!     -- AUTHOR(S)
!
!     SPACK.
!
!------------------------------------------------------------------------
 
      IMPLICIT NONE
 
 
 
      INTEGER  	, INTENT(IN)  :: ngas	               
      INTEGER  	, INTENT(IN)  :: ijkbeg 	       
      INTEGER  	, INTENT(IN)  :: ijkend 	       
      INTEGER  	, INTENT(IN)  :: maxblock_size	       
      INTEGER  	, INTENT(IN)  :: nr		       
      DOUBLE PRECISION , INTENT(IN)  :: rk(maxblock_size,nr)  
      DOUBLE PRECISION , INTENT(IN)  :: y(maxblock_size,NGAS) 
      DOUBLE PRECISION , INTENT(OUT) :: w(maxblock_size,nr)   
      INTEGER  	              :: ijk		       
 
      DO ijk=ijkbeg,ijkend
      w(ijk,  1) =  rk(ijk,  1) * Y(ijk,  1)
      w(ijk,  2) =  rk(ijk,  2) * Y(ijk,  4)
      w(ijk,  3) =  rk(ijk,  3) * Y(ijk,  4)
      w(ijk,  4) =  rk(ijk,  4) * Y(ijk, 11)
      w(ijk,  5) =  rk(ijk,  5) * Y(ijk, 10)
      w(ijk,  6) =  rk(ijk,  6) * Y(ijk, 12)
      w(ijk,  7) =  rk(ijk,  7) * Y(ijk,  5)
      w(ijk,  8) =  rk(ijk,  8) * Y(ijk,  5)
      w(ijk,  9) =  rk(ijk,  9) * Y(ijk, 13)
      w(ijk, 10) =  rk(ijk, 10) * Y(ijk, 15)
      w(ijk, 11) =  rk(ijk, 11) * Y(ijk, 15)
      w(ijk, 12) =  rk(ijk, 12) * Y(ijk, 16)
      w(ijk, 13) =  rk(ijk, 13) * Y(ijk, 33)
      w(ijk, 14) =  rk(ijk, 14) * Y(ijk, 28)
      w(ijk, 15) =  rk(ijk, 15) * Y(ijk, 30)
      w(ijk, 16) =  rk(ijk, 16) * Y(ijk,  3)
      w(ijk, 17) =  rk(ijk, 17) * Y(ijk,  3) * Y(ijk,  4)
      w(ijk, 18) =  rk(ijk, 18) * Y(ijk,  6)
      w(ijk, 19) =  rk(ijk, 19) * Y(ijk,  6)
      w(ijk, 20) =  rk(ijk, 20) * Y(ijk,  6)
      w(ijk, 21) =  rk(ijk, 21) * Y(ijk,  4) * Y(ijk,  7)
      w(ijk, 22) =  rk(ijk, 22) * Y(ijk,  4) * Y(ijk,  8)
      w(ijk, 23) =  rk(ijk, 23) * Y(ijk,  7) * Y(ijk,  8)
      w(ijk, 24) =  rk(ijk, 24) * Y(ijk, 13) * Y(ijk,  7)
      w(ijk, 25) =  rk(ijk, 25) * Y(ijk,  8) * Y(ijk,  8)
      w(ijk, 26) =  rk(ijk, 26) * Y(ijk,  8) * Y(ijk,  8)
      w(ijk, 27) =  rk(ijk, 27) * Y(ijk,  3) * Y(ijk,  2)
      w(ijk, 28) =  rk(ijk, 28) * Y(ijk,  3) * Y(ijk,  1)
      w(ijk, 29) =  rk(ijk, 29) * Y(ijk,  3) * Y(ijk,  1)
      w(ijk, 30) =  rk(ijk, 30) * Y(ijk,  7) * Y(ijk,  2)
      w(ijk, 31) =  rk(ijk, 31) * Y(ijk,  7) * Y(ijk,  1)
      w(ijk, 32) =  rk(ijk, 32) * Y(ijk,  7) * Y(ijk,  5)
      w(ijk, 33) =  rk(ijk, 33) * Y(ijk,  8) * Y(ijk,  2)
      w(ijk, 34) =  rk(ijk, 34) * Y(ijk,  8) * Y(ijk,  1)
      w(ijk, 35) =  rk(ijk, 35) * Y(ijk, 12)
      w(ijk, 36) =  rk(ijk, 36) * Y(ijk,  8) * Y(ijk,  5)
      w(ijk, 37) =  rk(ijk, 37) * Y(ijk,  7) * Y(ijk, 11)
      w(ijk, 38) =  rk(ijk, 38) * Y(ijk,  7) * Y(ijk, 10)
      w(ijk, 39) =  rk(ijk, 39) * Y(ijk,  7) * Y(ijk, 12)
      w(ijk, 40) =  rk(ijk, 40) * Y(ijk,  4) * Y(ijk,  2)
      w(ijk, 41) =  rk(ijk, 41) * Y(ijk,  4) * Y(ijk,  1)
      w(ijk, 42) =  rk(ijk, 42) * Y(ijk,  2) * Y(ijk,  2)
      w(ijk, 43) =  rk(ijk, 43) * Y(ijk,  2) * Y(ijk,  1)
      w(ijk, 44) =  rk(ijk, 44) * Y(ijk, 11) * Y(ijk, 11)
      w(ijk, 45) =  rk(ijk, 45) * Y(ijk,  5) * Y(ijk,  2)
      w(ijk, 46) =  rk(ijk, 46) * Y(ijk,  5) * Y(ijk,  1)
      w(ijk, 47) =  rk(ijk, 47) * Y(ijk,  5) * Y(ijk,  1)
      w(ijk, 48) =  rk(ijk, 48) * Y(ijk,  9)
      w(ijk, 49) =  rk(ijk, 49) * Y(ijk,  9)
      w(ijk, 50) =  rk(ijk, 50) * Y(ijk,  9)
      w(ijk, 51) =  rk(ijk, 51) * Y(ijk,  5) * Y(ijk,  5)
      w(ijk, 52) =  rk(ijk, 52) * Y(ijk, 14) * Y(ijk,  7)
      w(ijk, 53) =  rk(ijk, 53) * Y(ijk, 14) * Y(ijk,  7)
      w(ijk, 54) =  rk(ijk, 54) * Y(ijk, 15) * Y(ijk,  3)
      w(ijk, 55) =  rk(ijk, 55) * Y(ijk, 15) * Y(ijk,  7)
      w(ijk, 56) =  rk(ijk, 56) * Y(ijk, 15) * Y(ijk,  5)
      w(ijk, 57) =  rk(ijk, 57) * Y(ijk, 16) * Y(ijk,  3)
      w(ijk, 58) =  rk(ijk, 58) * Y(ijk, 16) * Y(ijk,  7)
      w(ijk, 59) =  rk(ijk, 59) * Y(ijk, 16) * Y(ijk,  5)
      w(ijk, 60) =  rk(ijk, 60) * Y(ijk, 17) * Y(ijk,  2)
      w(ijk, 61) =  rk(ijk, 61) * Y(ijk, 17) * Y(ijk,  1)
      w(ijk, 62) =  rk(ijk, 62) * Y(ijk, 19)
      w(ijk, 63) =  rk(ijk, 63) * Y(ijk, 17) * Y(ijk,  8)
      w(ijk, 64) =  rk(ijk, 64) * Y(ijk, 17) * Y(ijk, 17)
      w(ijk, 65) =  rk(ijk, 65) * Y(ijk,  7) * Y(ijk, 20)
      w(ijk, 66) =  rk(ijk, 66) * Y(ijk, 22)
      w(ijk, 67) =  rk(ijk, 67) * Y(ijk, 22)
      w(ijk, 68) =  rk(ijk, 68) * Y(ijk, 22) * Y(ijk,  1)
      w(ijk, 69) =  rk(ijk, 69) * Y(ijk,  3) * Y(ijk, 24)
      w(ijk, 70) =  rk(ijk, 70) * Y(ijk,  7) * Y(ijk, 24)
      w(ijk, 71) =  rk(ijk, 71) * Y(ijk,  4) * Y(ijk, 24)
      w(ijk, 72) =  rk(ijk, 72) * Y(ijk,  3) * Y(ijk, 23)
      w(ijk, 73) =  rk(ijk, 73) * Y(ijk,  7) * Y(ijk, 23)
      w(ijk, 74) =  rk(ijk, 74) * Y(ijk,  4) * Y(ijk, 23)
      w(ijk, 75) =  rk(ijk, 75) * Y(ijk,  5) * Y(ijk, 23)
      w(ijk, 76) =  rk(ijk, 76) * Y(ijk,  3) * Y(ijk, 32)
      w(ijk, 77) =  rk(ijk, 77) * Y(ijk,  7) * Y(ijk, 32)
      w(ijk, 78) =  rk(ijk, 78) * Y(ijk,  4) * Y(ijk, 32)
      w(ijk, 79) =  rk(ijk, 79) * Y(ijk,  5) * Y(ijk, 32)
      w(ijk, 80) =  rk(ijk, 80) * Y(ijk,  1) * Y(ijk, 32)
      w(ijk, 81) =  rk(ijk, 81) * Y(ijk, 33) * Y(ijk,  7)
      w(ijk, 82) =  rk(ijk, 82) * Y(ijk, 33) * Y(ijk,  4)
      w(ijk, 83) =  rk(ijk, 83) * Y(ijk, 33) * Y(ijk,  5)
      w(ijk, 84) =  rk(ijk, 84) * Y(ijk,  7) * Y(ijk, 25)
      w(ijk, 85) =  rk(ijk, 85) * Y(ijk, 27) * Y(ijk,  2)
      w(ijk, 86) =  rk(ijk, 86) * Y(ijk, 27)
      w(ijk, 87) =  rk(ijk, 87) * Y(ijk,  7) * Y(ijk, 26)
      w(ijk, 88) =  rk(ijk, 88) * Y(ijk,  5) * Y(ijk, 26)
      w(ijk, 89) =  rk(ijk, 89) * Y(ijk, 29) * Y(ijk,  1)
      w(ijk, 90) =  rk(ijk, 90) * Y(ijk, 28) * Y(ijk,  7)
      w(ijk, 91) =  rk(ijk, 91) * Y(ijk, 28) * Y(ijk,  4)
      w(ijk, 92) =  rk(ijk, 92) * Y(ijk,  7) * Y(ijk, 31)
      w(ijk, 93) =  rk(ijk, 93) * Y(ijk,  7) * Y(ijk, 30)
      w(ijk, 94) =  rk(ijk, 94) * Y(ijk, 18) * Y(ijk,  2)
      w(ijk, 95) =  rk(ijk, 95) * Y(ijk, 21) * Y(ijk,  2)
      w(ijk, 96) =  rk(ijk, 96) * Y(ijk, 18) * Y(ijk, 18)
      w(ijk, 97) =  rk(ijk, 97) * Y(ijk, 18) * Y(ijk,  8)
      w(ijk, 98) =  rk(ijk, 98) * Y(ijk, 21) * Y(ijk,  8)
      w(ijk, 99) =  rk(ijk, 99) * Y(ijk, 21) * Y(ijk, 21)
      w(ijk,100) =  rk(ijk,100) * Y(ijk, 18) * Y(ijk, 21)
      w(ijk,101) =  rk(ijk,101) * Y(ijk, 34) * Y(ijk,  7)
      w(ijk,102) =  rk(ijk,102) * Y(ijk, 35) * Y(ijk,  7)
      w(ijk,103) =  rk(ijk,103) * Y(ijk, 36) * Y(ijk,  7)
       END DO
 
   END SUBROUTINE rates
 
  END MODULE mod_chem_spack_rates
 
