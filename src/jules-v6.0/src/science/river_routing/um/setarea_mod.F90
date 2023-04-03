#if defined(UM_JULES)
!Huw Lewis (MO), Jan 2015
!DEPRECATED CODE
!This code was transferred from the UM repository at UM vn9.2 / JULES vn 4.1.
!Future developments will supercede these subroutines, and as such they
!should be considered deprecated. They will be retained in the codebase to
!maintain backward compatibility with functionality prior to
!UM vn10.0 / JULES vn 4.2, until such time as they become redundant.
!
!
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: River Routing

MODULE setarea_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SETAREA_MOD'

CONTAINS

SUBROUTINE setarea(nx, ny, area, jmax, offset_ny)

USE arealat1_mod, ONLY: arealat1

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
!
!     set area [m^2] of each grid box
!
INTEGER, INTENT(IN) :: nx, ny, jmax, offset_ny
REAL(KIND=real_jlslsm), INTENT(OUT) :: area(nx, ny)
REAL(KIND=real_jlslsm) :: atmp
INTEGER :: j_offset
INTEGER :: i, j
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SETAREA'

j_offset = offset_ny


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
DO j = 1, jmax
  IF ((j + j_offset) <= 90) THEN
    atmp = arealat1(INT(ABS(91.0 - (j + j_offset)))) * 1.0e06
  ELSE
    atmp = arealat1(INT(ABS((j + j_offset) - 90.0))) * 1.0e06
  END IF

  DO i = 1, nx
    area(i,j) = atmp
  END DO
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE setarea

END MODULE setarea_mod
#endif
