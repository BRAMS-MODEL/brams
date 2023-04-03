! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE calc_direct_albsoil_mod

USE um_types, ONLY: real_jlslsm

USE parkind1,                 ONLY: jprb, jpim
USE yomhook,                  ONLY: lhook, dr_hook

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CALC_DIRECT_ALBSOIL_MOD'

CONTAINS

FUNCTION calc_direct_albsoil(albdif, cosz) RESULT (albdir)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  albdif,                                                                     &
    ! Diffuse albedo of soil
  cosz
    ! Cosine of the solar zenith angle
REAL(KIND=real_jlslsm) ::                                                     &
  albdir
    ! Direct albedo, RESULT of function


REAL(KIND=real_jlslsm) ::                                                     &
  gamma_sim,                                                                  &
    ! Isotropic similarity parameter
  albdif_app
    ! Approximate diffuse albedo

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_DIRECT_ALBSOIL'

!End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


!           Pade approximant to the isotropic similarity parameter
!           using the diffuse albedo.
gamma_sim = (1.0 - albdif) /                                                  &
            (1.0 + (4.0 / 3.0) * albdif)
!           Calculate approximate diffuse albedo to scale direct
!           to a more accurate value
albdif_app = (1.0 - gamma_sim) *                                              &
             (1.0 - LOG(1.0+2.0 * gamma_sim) /                                &
             (2.0 * gamma_sim) ) / gamma_sim
!           Equations 24 and 40 from Hapke (1981) with scaling correction.
albdir = (albdif / albdif_app) *                                              &
                 (1.0 - gamma_sim) /                                          &
                 (1.0+2.0 * gamma_sim * cosz)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION calc_direct_albsoil


END MODULE calc_direct_albsoil_mod


