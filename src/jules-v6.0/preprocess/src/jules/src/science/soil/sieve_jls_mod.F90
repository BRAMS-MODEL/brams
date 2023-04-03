! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  SUBROUTINE SIEVE-----------------------------------------------------------

!  PURPOSE : TO CALCULATE THE THROUGHFALL OF WATER FALLING
!            THROUGH THE SURFACE CANOPY

!  SUITABLE FOR SINGLE COLUMN MODEL USE

!  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER NO 25
!                  SECTION (3B(II)), EQN(P252.9)

!-----------------------------------------------------------------------------

MODULE sieve_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SIEVE_MOD'

CONTAINS

SUBROUTINE sieve ( npnts, surft_pts, surft_index, area, can_cpy, r, frac,     &
                   timestep, can_wcnt, tot_tfall, tot_tfall_surft )

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(IN):
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  npnts,                                                                      &
    ! Total number of land points.
  surft_pts,                                                                  &
    ! Number of tile points.
  surft_index(npnts)
    ! Index of tile points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  area(npnts),                                                                &
    ! Fractional area of gridbox over which water falls (%).
  can_cpy(npnts),                                                             &
    ! Canopy capacity (kg/m2).
  r(npnts),                                                                   &
    ! Water fall rate (kg/m2/s).
  frac(npnts),                                                                &
    ! Tile fraction.
  timestep
    ! Timestep (s).

!-----------------------------------------------------------------------------
! Arguments with INTENT(INOUT):
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
  can_wcnt(npnts),                                                            &
    ! Canopy water content (kg/m2).
  tot_tfall(npnts),                                                           &
    ! Cummulative canopy throughfall (kg/m2/s).
  tot_tfall_surft(npnts)
    ! Throughfall contributions for this tile (kg/m2/s).

!-----------------------------------------------------------------------------
! Local variables:
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  i, j
    ! Loop counters:
    ! i for land point
    ! j for tile point

REAL(KIND=real_jlslsm) ::                                                     &
  aexp,                                                                       &
    ! Used in calculation of exponential in throughfall formula.
  can_ratio,                                                                  &
    ! CAN_WCNT / CAN_CPY
  tfall(npnts),                                                               &
    ! Local throughfall (kg/m2/s).
  smallp,                                                                     &
    ! Small positive number << 1.
  smallestp
    ! Smallest +ve real which can be represented.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SIEVE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!-----------------------------------------------------------------------------

smallestp = TINY(1.0)
smallp    = EPSILON( r )

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(j, i, aexp, can_ratio)                &
!$OMP SHARED(surft_pts, surft_index, can_cpy, r, smallp, area, timestep,      &
!$OMP        smallestp, can_wcnt, tfall, tot_tfall, frac, tot_tfall_surft)    &
!$OMP SCHEDULE(STATIC)
DO j = 1,surft_pts
  i = surft_index(j)
  IF ( can_cpy(i)  >  0.0 .AND. r(i) > smallp ) THEN
    aexp = area(i) * can_cpy(i) / (r(i) * timestep)
    ! Only calculate if AEXP is small enough to avoid underflow
    IF (aexp < -LOG(smallestp)) THEN
      aexp = EXP(-aexp)
    ELSE
      aexp = 0.0
    END IF

    can_ratio = can_wcnt(i) / can_cpy(i)

    can_ratio = MIN(can_ratio,1.0)
    tfall(i) = r(i) * ((1.0 - can_ratio) * aexp + can_ratio)
  ELSE
    tfall(i) = r(i)
  END IF

  can_wcnt(i)  = can_wcnt(i) + (r(i) - tfall(i)) * timestep
  tot_tfall(i) = tot_tfall(i) + frac(i) * tfall(i)

  ! Increment the tile throughfall for this precip type
  tot_tfall_surft(i) = tot_tfall_surft(i) + tfall(i)

END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE sieve
END MODULE sieve_mod
