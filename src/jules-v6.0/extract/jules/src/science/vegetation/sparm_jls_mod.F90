! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! *********************************************************************
! SUBROUTINE sparm
!
! Routine to calculate the gridbox mean land surface parameters from
! the areal fractions of the surface types and the structural
! properties of the plant functional types.
!
! This routine no longer calculates the max infiltration rate. This
! is now done in science/soil/infiltration_rate.F90 and called seperately
! as required.
!
! Future developments make calling infiltration_rate from within sparm
! undesirable (eg smcl-dependence)
!
! *********************************************************************
MODULE sparm_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SPARM_MOD'

CONTAINS

SUBROUTINE sparm (land_pts, nsurft, surft_pts, surft_index,                   &
                  frac_surft, canht_pft, lai_pft, z0m_soil_gb,                &
                  catch_snow_surft, catch_surft, z0_surft, z0h_bare_surft,    &
                  ztm_gb)

!Use in relevant subroutines
USE pft_sparm_mod,            ONLY: pft_sparm

!Use in relevant variables
USE jules_surface_types_mod,  ONLY: lake, npft, ntype, urban_canyon,          &
                                    urban_roof, soil
USE jules_vegetation_mod,     ONLY: can_model
USE blend_h,                  ONLY: lb
USE nvegparm,                 ONLY: catch_nvg, z0_nvg
USE jules_surface_mod,        ONLY: i_aggregate_opt, l_vary_z0m_soil,         &
                                    l_aggregate
USE jules_snow_mod,           ONLY: cansnowtile, snowloadlai
USE c_z0h_z0m,                ONLY: z0h_z0m
USE switches_urban,           ONLY: l_moruses
USE dust_param,               ONLY: z0_soil

USE stochastic_physics_run_mod, ONLY: l_rp2, i_rp_scheme, i_rp2b, z0hm_pft_rp

USE parkind1,                 ONLY: jprb, jpim
USE yomhook,                  ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                        &
  land_pts,                                                                   &
    ! Number of land points to be processed.
  nsurft,                                                                     &
    ! Number of surface tiles.
  surft_pts(ntype),                                                           &
    ! Number of land points which include the nth surface type.
  surft_index(land_pts,ntype)
    ! Indices of land points which include the nth surface type.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  frac_surft(land_pts,ntype),                                                 &
    ! Fractional cover of each surface type.
  canht_pft(land_pts,npft),                                                   &
    ! Vegetation height (m).
  lai_pft(land_pts,npft),                                                     &
    ! Leaf area index.
  z0m_soil_gb(land_pts),                                                      &
    ! z0m of bare soil (m)
  ztm_gb(land_pts)
    ! Roughness length

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  catch_snow_surft(land_pts,nsurft),                                          &
    ! Snow capacity for tile (kg/m2)
  catch_surft(land_pts,nsurft),                                               &
    ! Canopy capacity for each tile (kg/m2).
  z0_surft(land_pts,nsurft),                                                  &
    ! Roughness length for each tile (m).
  z0h_bare_surft(land_pts,nsurft)
    ! Snow-free thermal roughness length for each tile (m).


!-----------------------------------------------------------------------------
! Local variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  j,l,n  ! Loop counters

REAL(KIND=real_jlslsm) ::                                                     &
  catch(land_pts),                                                            &
    ! GBM canopy capacity (kg/m2).
  catch_t(land_pts,ntype),                                                    &
    ! Canopy capacities for types (kg/m2).
  fz0(land_pts),                                                              &
    ! Aggregation function of Z0.
  fz0h(land_pts),                                                             &
    ! Aggregation function of Z0H.
  z0(land_pts),                                                               &
    ! GBM roughness length (m).
  z0h(land_pts),                                                              &
    ! GBM thermal roughness length (m).
  z0_t(land_pts,ntype)
    ! Roughness lengths for types (m).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SPARM'

!-----------------------------------------------------------------------------
!end of header
!-----------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Set parameters for vegetated surface types
!-----------------------------------------------------------------------------
IF ( l_rp2 .AND. i_rp_scheme == i_rp2b) THEN
  DO n = 1,npft
    z0h_z0m(n) = z0hm_pft_rp(n)
  END DO
END IF

DO n = 1,npft
  CALL pft_sparm (land_pts, n, surft_pts(n), surft_index(:,n),                &
                  canht_pft(:,n), lai_pft(:,n), catch_t(:,n), z0_t(:,n))
END DO

! cansnowtile is only used when l_aggregate is .FALSE. - needs to be
! consistent with logic where cansnowtile is set.
IF (can_model  ==  4 .AND. .NOT. l_aggregate) THEN
  DO n = 1,npft
    IF ( cansnowtile(n) ) THEN
      DO j = 1,surft_pts(n)
        l = surft_index(j,n)
        catch_snow_surft(l,n) = snowloadlai * lai_pft(l,n)
      END DO
    END IF
  END DO
END IF

!-----------------------------------------------------------------------------
! Set parameters for non-vegetated surface types
!-----------------------------------------------------------------------------
DO n = npft+1,ntype
  DO j = 1,surft_pts(n)
    l = surft_index(j,n)
    catch_t(l,n) = catch_nvg(n - npft)
  END DO
END DO

IF (l_vary_z0m_soil) THEN
  DO n = npft+1,ntype
    IF (n == soil) THEN
      ! for soil, get the z0m from the input data:
      DO l = 1,land_pts
        z0_t(l,n) = z0m_soil_gb(l)
      END DO
    ELSE
      ! for non soil types, get the z0m from the input namelist:
      DO j = 1,surft_pts(n)
        l = surft_index(j,n)
        z0_t(l,n) = z0_nvg(n - npft)
      END DO
    END IF
  END DO
ELSE
  ! get all the non veg types:
  DO n = npft+1,ntype
    DO j = 1,surft_pts(n)
      l = surft_index(j,n)
      z0_t(l,n) = z0_nvg(n - npft)
    END DO
  END DO
END IF

! MORUSES Set canyon & roof roughness length
IF ( l_moruses ) THEN
  n = urban_canyon
  DO j = 1,surft_pts(n)
    l = surft_index(j,n)
    z0_t(l,n)          = ztm_gb(l)
    z0_t(l,urban_roof) = ztm_gb(l)
  END DO
END IF

IF ( l_aggregate ) THEN
  !---------------------------------------------------------------------------
  ! Form means and copy to tile arrays if required for aggregate tiles
  !---------------------------------------------------------------------------
  DO l = 1,land_pts
    catch(l)  = 0.0
    fz0(l)    = 0.0
    fz0h(l)   = 0.0
    z0(l)     = 0.0
  END DO

  DO n = 1,ntype
    DO j = 1,surft_pts(n)
      l = surft_index(j,n)
      fz0(l) = fz0(l) + frac_surft(l,n) / (LOG(lb / z0_t(l,n)))**2
      ! Explicit aggregation of z0h if required.
      IF (i_aggregate_opt == 1) THEN
        fz0h(l) = fz0h(l) + frac_surft(l,n) /                                 &
                  ( LOG(lb / z0_t(l,n)) * LOG(lb / (z0h_z0m(n) * z0_t(l,n))) )
      END IF
    END DO
  END DO

  DO l = 1,land_pts
    z0(l) = lb * EXP(-SQRT(1.0 / fz0(l)))
    ! Explicit aggregation of z0h if required.
    IF (i_aggregate_opt == 1) THEN
      z0h(l) = lb * EXP(-1.0 / (fz0h(l) * LOG(lb / z0(l))) )
    END IF
  END DO

  DO n = 1,ntype
    DO j = 1,surft_pts(n)
      l = surft_index(j,n)
      catch(l) = catch(l) + frac_surft(l,n) * catch_t(l,n)
    END DO
  END DO

  DO l = 1,land_pts
    ! Canopy capacity is average over non-lake surface types
    catch_surft(l,1) = 0.0
    IF ( lake > 0 ) THEN
      IF ( frac_surft(l,lake) < 1.0 ) THEN
        catch_surft(l,1) = catch(l) / (1.0 - frac_surft(l,lake))
      END IF
    END IF
    z0_surft(l,1) = z0(l)
    IF (i_aggregate_opt == 1) z0h_bare_surft(l,1) = z0h(l)
  END DO

ELSE
  !---------------------------------------------------------------------------
  ! .NOT. l_aggregate
  ! Copy surface-type arrays to tiles if separate tiles used
  !---------------------------------------------------------------------------
  DO n = 1,ntype
    DO j = 1,surft_pts(n)
      l = surft_index(j,n)
      catch_surft(l,n) = catch_t(l,n)
      z0_surft(l,n)    = z0_t(l,n)
    END DO
  END DO

END IF  !  l_aggregate

!-----------------------------------------------------------------------------
! Set bare soil roughness for use in 1 tile dust scheme
!-----------------------------------------------------------------------------
! The following line has been added to readlsta so to fix CRUNs
! and is now a duplication.
! When l_vary_z0m_soil is true, the variable is also passed down to the
! dust emission, so this isn't used.

z0_soil = z0_nvg(soil - npft)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE sparm
END MODULE sparm_mod
