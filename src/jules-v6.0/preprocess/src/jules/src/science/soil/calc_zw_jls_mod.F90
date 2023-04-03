! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   Subroutine CALC_ZW-------------------------------------------------

!   Purpose: To calculate the mean water table depth from the soil
!            moisture deficit as described in Koster et al., 2000.,
!            using the Newton-Raphson method.

! Documentation: UNIFIED MODEL DOCUMENTATION PAPER NO 25

MODULE calc_zw_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CALC_ZW_MOD'

CONTAINS

SUBROUTINE calc_zw(npnts, nshyd, soil_pts, soil_index,                        &
                   bexp, sathh, smcl, smclzw, smclsat, smclsatzw, v_sat, zw)

USE water_constants_mod,  ONLY: rho_water
USE jules_hydrology_mod,  ONLY: zw_max

USE parkind1,             ONLY: jprb, jpim
USE yomhook,              ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!Subroutine arguments
!-----------------------------------------------------------------------------
! Scalar arguments with intent(IN):
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  npnts,                                                                      &
    ! Number of gridpoints.
  nshyd,                                                                      &
    ! Number of soil moisture levels.
  soil_pts
    ! Number of soil points.

!-----------------------------------------------------------------------------
! Array arguments with intent(IN):
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  soil_index(npnts)
    ! Array of soil points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  bexp(npnts),                                                                &
    ! Clapp-Hornberger exponent.
  sathh(npnts),                                                               &
    ! Saturated soil water pressure (m).
  smcl(npnts,nshyd),                                                          &
    ! Total soil moisture contents of each layer (kg/m2).
  smclsat(npnts,nshyd),                                                       &
    ! Soil moisture contents of each layer at saturation (kg/m2).
  smclzw(npnts),                                                              &
    ! Moisture content in deep layer.
  smclsatzw(npnts),                                                           &
    ! Moisture content in deep layer at saturation.
  v_sat(npnts)
    ! Volumetric soil moisture.

!-----------------------------------------------------------------------------
! Array arguments with intent(INOUT):
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
  zw(npnts)
    ! Water table depth (m).

!-----------------------------------------------------------------------------
! Local parameters.
!-----------------------------------------------------------------------------
INTEGER, PARAMETER :: niter = 3
    ! Number of iterations.

!-----------------------------------------------------------------------------
! Local variables:
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  i, j, n, it
    !Loop counters

REAL(KIND=real_jlslsm) ::                                                     &
  zw_old,                                                                     &
    !zw from last timestep/iteration.
  fn,                                                                         &
    !fn calc in Newton-Raphson iteration.
  dfn,                                                                        &
    !derivative of fn.
  zwest,                                                                      &
  smd,                                                                        &
    !soil moisture deficit
  psisat
    !Saturated soil water pressure (m) (negative).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_ZW'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i,zw_old,smd,psisat,zwest,it,fn,dfn)                          &
!$OMP SHARED(soil_pts,soil_index,zw,smclsat,smcl,smclsatzw,smclzw,sathh,bexp, &
!$OMP        v_sat,nshyd,zw_max)
DO j = 1,soil_pts
  i = soil_index(j)
  zw_old = zw(i)

  !---------------------------------------------------------------------------
  ! Calculate total soil moisture deficit:
  !---------------------------------------------------------------------------
  smd = 0.0
  DO n = 1,nshyd
    smd = smd + ((smclsat(i,n) - smcl(i,n)) / rho_water)
  END DO

  smd = smd + ((smclsatzw(i) - smclzw(i)) / rho_water)

  !---------------------------------------------------------------------------
  ! Estimate water table depth assuming an equilibrium profile
  ! above water table (i.e. pressure head gradient and gravity balance):
  !---------------------------------------------------------------------------
  psisat = -sathh(i)
  zwest  = zw_old

  DO it = 1,niter
    IF (zwest <  0.0) THEN
      zwest = 0.0
    END IF

    IF (zwest >  zw_max) THEN
      zwest = zw_max
    END IF

    zw_old = zwest
    zwest  = zw_old

    ! Newton-Raphson. zw(next)=zw-f(zw)/f'(zw).
    fn = zwest * v_sat(i) - smd                                               &
         - v_sat(i) * bexp(i) / (bexp(i) - 1.0) * psisat                      &
         * ( 1.0 - ( 1.0 - zwest / psisat )**                                 &
           ( 1.0 - 1.0 / bexp(i) ) )       !   f(zw)
    dfn =  v_sat(i) + v_sat(i) * ( 1.0 - zwest / psisat )** (-1.0 / bexp(i))

    IF (ABS(dfn) >  1.0e-10) THEN
      zwest = zwest - fn / dfn
    END IF

    IF (zwest <  0) THEN
      zwest = 0.0
    END IF

    IF (zwest >  zw_max) THEN
      zwest = zw_max
    END IF
  END DO  !  iterations

  zw(i) = zwest

  IF (zw(i) <  0) THEN
    zw(i) = 0.0
  END IF

  IF (zw(i) >  zw_max) THEN
    zw(i) = zw_max
  END IF

END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE calc_zw
END MODULE calc_zw_mod
