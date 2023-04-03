! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE SOIL_HYD_WT---------------------------------------------------

! Description: Updates water table depth, and calculates drainage and surface
!              runoff.

! Documentation : UM Documentation Paper 25


MODULE soil_hyd_wt_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SOIL_HYD_WT_MOD'

CONTAINS

SUBROUTINE soil_hyd_wt (npnts, nshyd, soil_pts, soil_index,                   &
                        bexp, fw, sathh, timestep, v_sat,                     &
                        slow_runoff, smcl, surf_roff, w_flux,                 &
                        stf_slow_runoff,                                      &
                        zw, sthzw, qbase, qbase_l,                            &
                        drain, l_top, smclzw, smclsatzw, smclsat,             &
                        surf_roff_inc)

! Use relevant subroutines
USE calc_zw_mod,  ONLY: calc_zw

USE parkind1,     ONLY: jprb, jpim
USE yomhook,      ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(IN):
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  npnts,                                                                      &
    ! Number of gridpoints.
  nshyd,                                                                      &
    ! Number of soil moisture levels.
  soil_pts
    ! Number of soil points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  timestep
    ! Model timestep (s).

LOGICAL, INTENT(IN) ::                                                        &
  stf_slow_runoff,                                                            &
    ! Stash flag for sub-surface runoff.
  l_top
    ! Flag for TOPMODEL-based hydrology.

INTEGER, INTENT(IN) ::                                                        &
  soil_index(npnts)
    ! Array of soil points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  bexp(npnts,nshyd),                                                          &
    ! Clapp-Hornberger exponent.
  fw(npnts),                                                                  &
    ! Throughfall from canopy plus snowmelt minus surface runoff (kg/m2/s).
  sathh(npnts,nshyd),                                                         &
    ! Saturated soil water pressure (m).
  v_sat(npnts,nshyd),                                                         &
    ! Volumetric soil moisture concentration at saturation (m3 H2O/m3 soil).
  w_flux(npnts,0:nshyd)
    ! The fluxes of water between layers (kg/m2/s).

!-----------------------------------------------------------------------------
! Arguments with INTENT(OUT):
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  slow_runoff(npnts),                                                         &
    ! Drainage from the base of the soil profile (kg/m2/s).
  drain(npnts),                                                               &
    ! Drainage out of nshyd'th level (kg/m2/s).
  qbase(npnts),                                                               &
    ! Base flow (kg/m2/s).
  surf_roff_inc(npnts)
    ! Increment to surface runoff (kg m-2 s-1).

!-----------------------------------------------------------------------------
! Arguments with INTENT(INOUT):
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
  smcl(npnts,nshyd),                                                          &
    ! Total soil moisture contents of each layer (kg/m2).
  surf_roff(npnts),                                                           &
    ! Surface runoff (kg/m2/s).
  sthzw(npnts),                                                               &
    ! Soil moisture fraction in deep layer.
  qbase_l(npnts,nshyd+1),                                                     &
    ! Base flow from each level (kg/m2/s).
  smclzw(npnts),                                                              &
    ! Moisture content in deep layer (kg/m2).
  smclsatzw(npnts),                                                           &
    ! Moisture content in deep layer at saturation (kg/m2).
  smclsat(npnts,nshyd),                                                       &
    ! The saturation moisture content of each layer (kg/m2).
  zw(npnts)
    ! Mean water table depth (m).

!-----------------------------------------------------------------------------
! Local scalars:
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  i, j, n                ! WORK Loop counters.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SOIL_HYD_WT'

!End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!-----------------------------------------------------------------------------

IF (l_top) THEN

  !---------------------------------------------------------------------------
  ! Diagnose the new water table depth.
  ! Assume local equilibrium psi profile.
  !---------------------------------------------------------------------------

  DO j = 1,soil_pts
    i = soil_index(j)
    smclzw(i) = smclzw(i) - (qbase_l(i,nshyd+1) - w_flux(i,nshyd)) * timestep
    ! Update prognostic deep layer soil moisture fraction:
    sthzw(i) = smclzw(i) / smclsatzw(i)
  END DO

  CALL calc_zw(npnts, nshyd, soil_pts, soil_index,                            &
               bexp, sathh, smcl, smclzw, smclsat, smclsatzw, v_sat, zw)

  !---------------------------------------------------------------------------
  ! Dont allow negative base flows:
  !---------------------------------------------------------------------------
  DO j = 1,soil_pts
    i = soil_index(j)
    qbase(i) = 0.0
    DO n = 1,nshyd+1
      qbase_l(i,n) = MAX(qbase_l(i,n),0.0)
      qbase(i)     = qbase(i) + qbase_l(i,n)
    END DO
  END DO

END IF
!-----------------------------------------------------------------------------
! Output slow runoff (drainage) diagnostic.
!-----------------------------------------------------------------------------
IF (stf_slow_runoff) THEN
  DO i = 1,npnts
    slow_runoff(i) = 0.0
  END DO
  DO j = 1,soil_pts
    i = soil_index(j)
    ! Ensure correct field is output as deep runoff: dependant on L_TOP
    IF (l_top) THEN
      slow_runoff(i) = qbase(i)
    ELSE
      slow_runoff(i) = w_flux(i,nshyd)
    END IF
    drain(i) = w_flux(i,nshyd)
  END DO
END IF

!-----------------------------------------------------------------------------
! Update surface runoff diagnostic.
!-----------------------------------------------------------------------------
DO j = 1,soil_pts
  i = soil_index(j)
  surf_roff(i)     = surf_roff(i) + (fw(i) - w_flux(i,0))
  surf_roff_inc(i) = fw(i) - w_flux(i,0)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE soil_hyd_wt
END MODULE soil_hyd_wt_mod
