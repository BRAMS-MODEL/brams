#if defined(UM_JULES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: river routing

MODULE init_riv_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='INIT_RIV_MOD'

CONTAINS
SUBROUTINE init_riv(nstep_since_riv, icode, cmessage)

! Purpose: To initialise variables for river routing when used with MetUM. 
! Method: Sets up variables in the dump

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE um_latlon_mod, ONLY: um_latlon

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
INTEGER :: icode                 ! OUT Internal return code
CHARACTER(LEN=errormessagelength) :: cmessage
                              ! OUT Internal error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_RIV'

INTEGER, INTENT(OUT) :: nstep_since_riv


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,                &
                        zhook_handle)
icode = 0
cmessage=""
nstep_since_riv = 0

! Set global (full domain) latitude, longitude and mask variables
! using UM module variables
CALL um_latlon()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,               &
                        zhook_handle)
RETURN
END SUBROUTINE init_riv
END MODULE init_riv_mod
#endif
