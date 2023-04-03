#if defined(UM_JULES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt and UM file CodeOwners.txt

MODULE check_jules_unavailable_options_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE check_jules_unavailable_options()

USE ereport_mod, ONLY: ereport
USE jules_print_mgr, ONLY:                                                    &
    jules_message, jules_print, PrNorm, newline
USE jules_vegetation_mod, ONLY: photo_acclim_model, photo_collatz,            &
    photo_farquhar, photo_model, stomata_jacobs, stomata_model,               &
    l_gleaf_fix, l_o3_damage, l_use_pft_psi, fsmc_shape, l_prescsow,          &
    l_croprotate
USE jules_irrig_mod, ONLY: l_irrig_dmd, l_irrig_limit
USE switches_urban, ONLY: l_urban_empirical
USE jules_rivers_mod, ONLY: i_river_vn, rivers_trip
USE jules_soil_biogeochem_mod, ONLY: l_ch4_microbe
USE jules_soil_mod, ONLY: l_tile_soil
USE jules_surface_mod, ONLY: l_surface_type_ids
USE jules_surface_types_mod, ONLY: ncpft
USE overbank_inundation_mod, ONLY: l_riv_overbank
USE jules_water_resources_mod, ONLY: l_water_resources
USE jules_model_environment_mod, ONLY: lsm_id, cable

IMPLICIT NONE

!Local variables
INTEGER :: errcode, error_sum
CHARACTER(LEN=*), PARAMETER :: RoutineName='CHECK_JULES_UNAVAILABLE_OPTIONS'


error_sum = 0
! jules_model_environment_mod
IF ( lsm_id == cable ) THEN
  error_sum = error_sum + 1
  WRITE(jules_message,'(I0,A,I0)') error_sum,                                 &
     ": CABLE can only be used in standalone mode. lsm_id = ", lsm_id
  CALL jules_print(RoutineName, jules_message, level = PrNorm)
END IF

! switches_urban
! An ereport has not been added for l_cosz as it is .true. by default. If the
! default is changed to .false. it can be added here.
IF ( l_urban_empirical ) THEN
  error_sum = error_sum + 1
  WRITE(jules_message,'(I0,A)') error_sum,                                    &
     ": l_urban_empirical should be .false."
  CALL jules_print(RoutineName, jules_message, level = PrNorm)
END IF

! jules_rivers
IF ( i_river_vn == rivers_trip ) THEN
  error_sum = error_sum + 1
  WRITE(jules_message,'(I0,A,I0)') error_sum,                                 &
     ": Rivers trip can only be run in standalone mode. i_river_vn = ",       &
     rivers_trip
  CALL jules_print(RoutineName, jules_message, level = PrNorm)
END IF

! jules_soil_biogeochem_mod
IF ( l_ch4_microbe ) THEN
  error_sum = error_sum + 1
  WRITE(jules_message,'(I0,A,L1)') error_sum,                                 &
     ": Microbial methane scheme is only available to standalone JULES. " //  &
     "l_ch4_microbe = ", l_ch4_microbe
  CALL jules_print(RoutineName, jules_message, level = PrNorm)
END IF

! jules_soil_mod
IF ( l_tile_soil ) THEN
  error_sum = error_sum + 1
  WRITE(jules_message,'(I0,A,L1)') error_sum,                                 &
     ": Soil tiling is only available to standalone JULES. " //               &
     "l_tile_soil = ", l_tile_soil
  CALL jules_print(RoutineName, jules_message, level = PrNorm)
END IF

! jules_surface_types
IF ( ncpft > 0 ) THEN
  error_sum = error_sum + 1
  WRITE(jules_message,'(I0,A,I0)') error_sum,                                 &
     ": JULES crop can only be run in standalone. ncpft = ", ncpft
  CALL jules_print(RoutineName, jules_message, level = PrNorm)
END IF

! jules_vegetation
IF ( fsmc_shape == 1 ) THEN
  error_sum = error_sum + 1
  WRITE(jules_message,'(I0,A,I0)') error_sum,                                 &
     ": Piece-wise linear in soil potential (fscm_shape = 1) is only " //     &
     "available to standalone. fscm_shape= ", fsmc_shape
  CALL jules_print(RoutineName, jules_message, level = PrNorm)
END IF

! photo_farquhar is not allowed yet, but if it is allowed further work will
! still be required to allow photo_acclim_model /= 0.
!*******************************************************************************
! Clarification: IF photo_model == photo_collatz
! * photo_acclim_model does nothing regardless of value of photo_acclim_model.
! * photo_acclim_model == imdi when photo_model == photo_collatz, therefore the
!   photo_model == photo_farquhar condition is required here. It appears
!   inconsistent with the metadata, however it is not.
IF ( photo_model == photo_farquhar .AND. photo_acclim_model /= 0 ) THEN
  error_sum = error_sum + 1
  WRITE(jules_message,'(I0,A)') error_sum,                                    &
     ": photo_acclim_model must be set to zero (no acclimation)."
  CALL jules_print(RoutineName, jules_message, level = PrNorm)
END IF

IF ( photo_model /= photo_collatz ) THEN
  error_sum = error_sum + 1
  WRITE(jules_message,'(I0,A)') error_sum,                                    &
     ": photo_model must be set to use the Collatz model."
  CALL jules_print(RoutineName, jules_message, level = PrNorm)
END IF

IF ( stomata_model /= stomata_jacobs ) THEN
  error_sum = error_sum + 1
  WRITE(jules_message,'(I0,A)') error_sum,                                    &
     ": stomata_model must be set to use the original (Jacobs) model."
  CALL jules_print(RoutineName, jules_message, level = PrNorm)
END IF

IF ( l_irrig_dmd ) THEN
  error_sum = error_sum + 1
  WRITE(jules_message,'(I0,A,L1)') error_sum,                                 &
     ": Irrigation demand model is only available to standalone. " //         &
     "l_irrig_dmd = ", l_irrig_dmd
  CALL jules_print(RoutineName, jules_message, level = PrNorm)
END IF

IF ( l_irrig_limit ) THEN
  error_sum = error_sum + 1
  WRITE(jules_message,'(I0,A,L1)') error_sum,                                 &
     ": Limiting irrigation supply is only available to standalone. " //      &
     "l_irrig_dmd = ", l_irrig_limit
  CALL jules_print(RoutineName, jules_message, level = PrNorm)
END IF

IF ( l_gleaf_fix ) THEN
  error_sum = error_sum + 1
  WRITE(jules_message,'(I0,A,L1)') error_sum,                                 &
     ": Bug fix for accumuating g_leaf_phen_acc between calls to " //         &
     "TRIFFID is for standalone only. l_gleaf_fix = ", l_gleaf_fix
  CALL jules_print(RoutineName, jules_message, level = PrNorm)
END IF

IF ( l_o3_damage ) THEN
  error_sum = error_sum + 1
  WRITE(jules_message,'(I0,A,L1)') error_sum,                                 &
     ": Ozone damage for vegetation is only available to standalone. " //     &
     "l_o3_damage = ", l_o3_damage
  CALL jules_print(RoutineName, jules_message, level = PrNorm)
END IF

IF ( l_prescsow ) THEN
  error_sum = error_sum + 1
  WRITE(jules_message,'(I0,A,L1)') error_sum,                                 &
     ": Prescribed sowing dates for crops is only available to " //           &
     "standalone. l_prescsow = ", l_prescsow
  CALL jules_print(RoutineName, jules_message, level = PrNorm)
END IF

IF ( l_use_pft_psi ) THEN
  error_sum = error_sum + 1
  WRITE(jules_message,'(I0,A,L1)') error_sum,                                 &
     ": Use psi_open_io and psi_close_io to calculate the soil moisture " //  &
     "stress function on vegetation is only available to standalone. " //     &
     "l_use_pft_psi = ", l_use_pft_psi
  CALL jules_print(RoutineName, jules_message, level = PrNorm)
END IF

IF ( l_croprotate ) THEN
  error_sum = error_sum + 1
  WRITE(jules_message,'(I0,A,L1)') error_sum,                                 &
     ": Double cropping is only available to standalone. " //                 &
     "l_croprotate = ", l_croprotate
  CALL jules_print(RoutineName, jules_message, level = PrNorm)
END IF

IF ( l_riv_overbank ) THEN
  error_sum = error_sum + 1
  WRITE(jules_message,'(I0,A,L1)') error_sum,                                 &
     ": Overbank inundation is only available to standalone. " //             &
     "l_riv_overbank = ", l_riv_overbank
  CALL jules_print(RoutineName, jules_message, level = PrNorm)
END IF

! jules_water_resources
IF ( l_water_resources ) THEN
  error_sum = error_sum + 1
  WRITE(jules_message,'(I0,A,L1)') error_sum,                                 &
     ": Water resource modelling is only available to standalone. " //        &
     "l_water_resources = ", l_water_resources
  CALL jules_print(RoutineName, jules_message, level = PrNorm)
END IF

IF ( error_sum > 0 ) THEN
  errcode = 30
  WRITE(jules_message,'(A,I0,A)') "One or more JULES options (", error_sum,   &
     ") have been incorrectly set for use in UM-JULES. " //                   &
     "Please see job output for details."
  CALL ereport(RoutineName, errcode, jules_message)
END IF

! Warnings for information
IF ( l_surface_type_ids ) THEN
  errcode = -10
  WRITE(jules_message,'(A,L1)') newline //                                    &
     " l_surface_type_ids has no effect in the UM and should be " //newline// &
     "triggered off, however it may indicate that something else "//newline// &
     "is wrong. Was this intended? l_surface_type_ids = ", l_surface_type_ids
  CALL ereport(RoutineName, errcode, jules_message)
END IF


RETURN
END SUBROUTINE check_jules_unavailable_options
END MODULE check_jules_unavailable_options_mod
#endif
