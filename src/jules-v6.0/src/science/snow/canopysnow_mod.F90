! *****************************COPYRIGHT*******************************

! (c) [University of Edinburgh]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
!  SUBROUTINE CANOPYSNOW------------------------------------------------------

! Description:
!     Partition snowfall into canopy interception, throughfall and unloading.

! Subroutine Interface:
MODULE canopysnow_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CANOPYSNOW_MOD'

CONTAINS

SUBROUTINE canopysnow ( land_pts, surft_pts, timestep, cansnowtile,           &
                        surft_index, catch_snow, con_snow, ls_snow, ls_graup, &
                        unload_backgrnd_surft,                                &
                        melt_surft, snow_can, snowfall, graupfall )

USE jules_snow_mod, ONLY:                                                     &
  ! imported scalars
  snowinterceptfact,                                                          &
  ! Constant in relationship between mass of intercepted snow and
  ! snowfall rate.
  snowunloadfact
  ! Constant in relationship between canopy snow unloading and canopy
  ! snow melt rate.

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER,INTENT(IN) ::                                                         &
  land_pts,                                                                   &
    !  Number of land points.
  surft_pts
    !  Number of tile points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  timestep                 ! Timestep (s).

LOGICAL, INTENT(IN) ::                                                        &
  cansnowtile              ! Switch for canopy snow model.

!-----------------------------------------------------------------------------
! Array arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  surft_index(land_pts)    ! Index of tile points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  catch_snow(land_pts),                                                       &
    ! Canopy snow capacity (kg/m2).
  con_snow(land_pts),                                                         &
    ! Convective snowfall rate (kg/m2/s).
  ls_snow(land_pts),                                                          &
    ! Large-scale frozen precip fall rate (kg/m2/s).
  ls_graup(land_pts),                                                         &
    ! Large-scale graupel fall rate (kg/m2/s).
  melt_surft(land_pts),                                                       &
    ! Canopy snow melt rate (kg/m2/s).
  unload_backgrnd_surft(land_pts)
    ! Background canopy unloading rate.

!-----------------------------------------------------------------------------
! Array arguments with intent(inout)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
  snow_can(land_pts)
    ! Canopy snow load, if the canopy snow model is selected (kg/m2).

!-----------------------------------------------------------------------------
! Array arguments with intent(out)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  snowfall(land_pts),                                                         &
    ! Frozen precip reaching the ground in timestep (kg/m2).
  graupfall(land_pts)
    ! Graupel reaching the ground in timestep (kg/m2).

!-----------------------------------------------------------------------------
! Local scalars
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  i,                                                                          &
    ! Land point index
  k
    ! Tile point index

REAL(KIND=real_jlslsm) ::                                                     &
  intercept,                                                                  &
    ! Snow intercepted by canopy in timestep (kg/m2).
 unload
    ! Canopy snow unloaded in timestep (kg/m2).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CANOPYSNOW'

!-----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( cansnowtile ) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,intercept,unload)                                           &
!$OMP SHARED(surft_pts,surft_index,snowfall,ls_snow,ls_graup,con_snow,        &
!$OMP        timestep,graupfall,catch_snow,snow_can,melt_surft,               &
!$OMP        unload_backgrnd_surft,snowinterceptfact,snowunloadfact)
  DO k = 1,surft_pts
    i = surft_index(k)
    snowfall(i)  = ( ls_snow(i) - ls_graup(i) + con_snow(i) ) * timestep
    graupfall(i) = ls_graup(i) * timestep
    intercept    = snowinterceptfact * (catch_snow(i) - snow_can(i))          &
                   * (1.0 - EXP(-snowfall(i) / catch_snow(i)))
    unload       = snowunloadfact * melt_surft(i) * timestep                  &
                   + unload_backgrnd_surft(i) * snow_can(i) * timestep
    !-----------------------------------------------------------------------
    ! At this point, the value of unload can be larger than the
    ! amount of snow on the canopy (which has already had melt
    ! and sublimation removed), so we need to limit to the amount
    ! of snow available. However, snow_can can be <0 (with small
    ! absolute value) because of "issues" in the surface flux code,
    ! so we also restrict unload to be >=0.
    !-----------------------------------------------------------------------
    unload      = MAX( MIN( unload, snow_can(i) ), 0.0 )
    snow_can(i) = snow_can(i) + intercept - unload
    snowfall(i) = snowfall(i) + graupfall(i) - intercept + unload
  END DO
!$OMP END PARALLEL DO
ELSE

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i)                                                            &
!$OMP SHARED(surft_pts,surft_index,snowfall,ls_snow,con_snow,timestep,        &
!$OMP        graupfall,ls_graup)
  DO k = 1,surft_pts
    i = surft_index(k)
    snowfall(i)  = ( ls_snow(i) + con_snow(i) ) * timestep
    graupfall(i) = ls_graup(i) * timestep
  END DO
!$OMP END PARALLEL DO

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE canopysnow
END MODULE canopysnow_mod
