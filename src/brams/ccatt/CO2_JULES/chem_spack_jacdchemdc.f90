                                                                                                                                                      
 MODULE mod_chem_spack_jacdchemdc                                                                                                                     
                                                                                                                                                      
   USE mod_chem_spack_dratedc, ONLY: dratedc  ! subroutine                                                                                            
   IMPLICIT NONE                                                                                                                                      
   PRIVATE                                                                                                                                            
   PUBLIC :: jacdchemdc ! subroutine                                                                                                                  
 CONTAINS                                                                                                                                             
                                                                                                                                                      
   SUBROUTINE jacdchemdc(y,rk,JacC,ngas,ijkbeg,ijkend,maxblock_size,nr)                                                                               
                                                                                                                                                      
!------------------------------------------------------------------------                                                                             
!                                                                                                                                                     
!     -- DESCRIPTION                                                                                                                                  
!                                                                                                                                                     
!     This routine computes the Jacobian matrix for the gas-phase.                                                                                    
!     This routine is automatically generated by SPACK.                                                                                               
!     Mechanism: ../Mechanism/CB07                                                                                                                    
!     Species: ../Mechanism/ciCB07                                                                                                                    
!                                                                                                                                                     
!------------------------------------------------------------------------                                                                             
!                                                                                                                                                     
!     -- INPUT VARIABLES                                                                                                                              
!                                                                                                                                                     
!     Y: chemical concentrations.                                                                                                                     
!     RK: kinetic rates.                                                                                                                              
!                                                                                                                                                     
!     -- INPUT/OUTPUT VARIABLES                                                                                                                       
!                                                                                                                                                     
!     -- OUTPUT VARIABLES                                                                                                                             
!                                                                                                                                                     
!     JACC: Jacobian matrix.                                                                                                                          
!                                                                                                                                                     
!------------------------------------------------------------------------                                                                             
!                                                                                                                                                     
!     -- REMARKS                                                                                                                                      
!                                                                                                                                                     
!     The matrix JACC could be stored in a low-dimensional vector.                                                                                    
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
                                                                                                                                                      
                                                                                                                                                      
       INTEGER 	 , INTENT(IN)  :: ngas                                                                                                                
       INTEGER 	 , INTENT(IN)  :: ijkbeg			                                                                                                           
       INTEGER 	 , INTENT(IN)  :: ijkend			                                                                                                           
       INTEGER 	 , INTENT(IN)  :: maxblock_size 		                                                                                                    
       INTEGER 	 , INTENT(IN)  :: nr 				                                                                                                             
       DOUBLE PRECISION , INTENT(IN)  :: rk(maxblock_size,nr)		                                                                                       
       DOUBLE PRECISION , INTENT(IN)  :: y(maxblock_size,NGAS) 	                                                                                      
       DOUBLE PRECISION , INTENT(OUT) :: JacC(maxblock_size,NGAS,NGAS)                                                                                
       								                                                                                                                                       
       DOUBLE PRECISION :: dw(maxblock_size,nr,NGAS)			                                                                                               
       INTEGER :: ijk							                                                                                                                          
                                                                                                                                                      
                                                                                                                                                      
      CALL dratedc(rk,y,dw,ngas,ijkbeg,ijkend,maxblock_size,nr)                                                                                       
print*,"asd",maxblock_size,NGAS,nr,ijkbeg,ijkend                                                                                
                                                                                                                                                      
      DO ijk=ijkbeg,ijkend                                                                                                                            
      JacC(ijk,  1,  1) =  - dw(ijk,  1,  1) &
                           - dw(ijk, 28,  1) &
                           - dw(ijk, 29,  1) &
                           - dw(ijk, 31,  1) &
                           - dw(ijk, 34,  1) &
                           - dw(ijk, 41,  1) &
                           - dw(ijk, 43,  1) &
                           - dw(ijk, 47,  1) &
                           - dw(ijk, 61,  1) &
                           - dw(ijk, 68,  1) &
                           - dw(ijk, 80,  1) &
                           - dw(ijk, 89,  1)
      JacC(ijk,  2,  1) =  + dw(ijk,  1,  1) &
                           + dw(ijk, 28,  1) &
                           - dw(ijk, 43,  1) &
                           + dw(ijk, 46,  1) &
                          + 0.2000000000000000D+00*dw(ijk, 80,  1)
      JacC(ijk,  3,  1) =  + dw(ijk,  1,  1) &
                           - dw(ijk, 28,  1) &
                           - dw(ijk, 29,  1)
      JacC(ijk,  4,  4) =  - dw(ijk,  2,  4) &
                           - dw(ijk,  3,  4) &
                           - dw(ijk, 17,  4) &
                           - dw(ijk, 21,  4) &
                           - dw(ijk, 22,  4) &
                           - dw(ijk, 40,  4) &
                           - dw(ijk, 41,  4) &
                           - dw(ijk, 71,  4) &
                           - dw(ijk, 74,  4) &
                           - dw(ijk, 78,  4) &
                           - dw(ijk, 82,  4) &
                           - dw(ijk, 91,  4)
      JacC(ijk,  6,  4) =  + dw(ijk,  2,  4)
      JacC(ijk,  3,  4) =  + dw(ijk,  3,  4) &
                           - dw(ijk, 17,  4)
      JacC(ijk,  2,  5) =  + dw(ijk,  7,  5) &
                           - dw(ijk, 45,  5) &
                           + dw(ijk, 46,  5)
      JacC(ijk,  5,  5) =  - dw(ijk,  7,  5) &
                           - dw(ijk,  8,  5) &
                           - dw(ijk, 32,  5) &
                           - dw(ijk, 36,  5) &
                           - dw(ijk, 45,  5) &
                           - dw(ijk, 46,  5) &
                           - dw(ijk, 47,  5) &
                          - 0.2000000000000000D+01*dw(ijk, 51,  5) &
                          - 0.2000000000000000D+01*dw(ijk, 51,  5) &
                           - dw(ijk, 56,  5) &
                           - dw(ijk, 59,  5) &
                           - dw(ijk, 75,  5) &
                           - dw(ijk, 79,  5) &
                           - dw(ijk, 83,  5) &
                           - dw(ijk, 88,  5)
      JacC(ijk,  1,  5) =  + dw(ijk,  8,  5) &
                           + dw(ijk, 32,  5) &
                          + 0.7000000000000000D+00*dw(ijk, 36,  5) &
                          + 0.2000000000000000D+01*dw(ijk, 45,  5) &
                           - dw(ijk, 47,  5) &
                          + 0.2000000000000000D+01*dw(ijk, 51,  5) &
                          + 0.2000000000000000D+01*dw(ijk, 51,  5) &
                           + dw(ijk, 75,  5) &
                          + 0.2000000000000000D+00*dw(ijk, 79,  5)
      JacC(ijk,  3,  5) =  + dw(ijk,  8,  5)
      JacC(ijk,  3,  3) =  - dw(ijk, 16,  3) &
                           - dw(ijk, 17,  3) &
                           - dw(ijk, 27,  3) &
                           - dw(ijk, 28,  3) &
                           - dw(ijk, 29,  3) &
                           - dw(ijk, 54,  3) &
                           - dw(ijk, 57,  3) &
                           - dw(ijk, 69,  3) &
                           - dw(ijk, 72,  3) &
                           - dw(ijk, 76,  3)
      JacC(ijk,  4,  3) =  + dw(ijk, 16,  3) &
                           - dw(ijk, 17,  3)
      JacC(ijk,  3,  6) =  + dw(ijk, 18,  6) &
                           + dw(ijk, 19,  6)
      JacC(ijk,  6,  6) =  - dw(ijk, 18,  6) &
                           - dw(ijk, 19,  6) &
                           - dw(ijk, 20,  6)
      JacC(ijk,  1,  3) =  + dw(ijk, 27,  3) &
                           - dw(ijk, 28,  3) &
                           - dw(ijk, 29,  3)
      JacC(ijk,  1,  2) =  + dw(ijk, 27,  2) &
                           + dw(ijk, 33,  2) &
                           + dw(ijk, 40,  2) &
                          + 0.2000000000000000D+01*dw(ijk, 42,  2) &
                          + 0.2000000000000000D+01*dw(ijk, 42,  2) &
                           - dw(ijk, 43,  2) &
                          + 0.2000000000000000D+01*dw(ijk, 45,  2) &
                           + dw(ijk, 60,  2) &
                          + 0.9000000000000000D+00*dw(ijk, 85,  2) &
                           + dw(ijk, 94,  2)
      JacC(ijk,  2,  3) =  - dw(ijk, 27,  3) &
                           + dw(ijk, 28,  3)
      JacC(ijk,  2,  2) =  - dw(ijk, 27,  2) &
                           - dw(ijk, 30,  2) &
                           - dw(ijk, 33,  2) &
                           - dw(ijk, 40,  2) &
                          - 0.2000000000000000D+01*dw(ijk, 42,  2) &
                          - 0.2000000000000000D+01*dw(ijk, 42,  2) &
                           - dw(ijk, 43,  2) &
                           - dw(ijk, 45,  2) &
                           - dw(ijk, 60,  2) &
                           - dw(ijk, 85,  2) &
                           - dw(ijk, 94,  2) &
                           - dw(ijk, 95,  2)
      JacC(ijk,  3,  2) =  - dw(ijk, 27,  2)
      JacC(ijk,  5,  3) =  + dw(ijk, 29,  3)
      JacC(ijk,  5,  1) =  + dw(ijk, 29,  1) &
                           + dw(ijk, 41,  1) &
                           - dw(ijk, 46,  1) &
                           - dw(ijk, 47,  1)
      JacC(ijk,  1,  4) =  + dw(ijk, 40,  4) &
                           - dw(ijk, 41,  4)
      JacC(ijk,  2,  4) =  - dw(ijk, 40,  4)
      JacC(ijk,  4,  2) =  - dw(ijk, 40,  2)
      JacC(ijk,  4,  1) =  - dw(ijk, 41,  1)
      JacC(ijk,  5,  4) =  + dw(ijk, 41,  4)
      JacC(ijk,  5,  2) =  - dw(ijk, 45,  2)
       END DO                                                                                                                                         
print*,"asd22",maxblock_size,NGAS,nr,ijkbeg,ijkend                                                                                
                                                                                                                                                      
   END SUBROUTINE jacdchemdc                                                                                                                          
                                                                                                                                                      
  END MODULE mod_chem_spack_jacdchemdc                                                                                                                
                                                                                                                                                      
