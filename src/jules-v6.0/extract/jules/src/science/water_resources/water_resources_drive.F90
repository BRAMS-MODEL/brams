#if !defined(UM_JULES)
!******************************COPYRIGHT**************************************
! (c) UK Centre for Ecology & Hydrology.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT**************************************

MODULE water_resources_drive_mod

!-----------------------------------------------------------------------------
! Description:
!   Driver routine for water resource management.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in HYDROLOGY
!
! Code Description:
!   Language: Fortran 90.
!-----------------------------------------------------------------------------

IMPLICIT NONE

PRIVATE  !  private scope by default
PUBLIC water_resources_drive

! Module parameters.
CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                       &
  ModuleName = 'WATER_RESOURCES_DRIVE_MOD'

CONTAINS

!#############################################################################

SUBROUTINE water_resources_drive( ntype, land_index, plant_n_gb,              &
             dvi_cpft, flandg, frac_irr_soilt, frac_soilt,                    &
             irrig_eff, grid_area_ij, smvccl_soilt, smvcst_soilt, sthf_soilt, &
             demand_accum, smcl_soilt, sthu_irr_soilt, sthu_soilt,            &
             irrig_water_gb )

USE ancil_info, ONLY: land_pts, nsoilt

USE atm_fields_bounds_mod, ONLY: pdims_s

USE irrigation_water_mod, ONLY: add_irrigation_to_soil, calc_irrigation_demand

USE jules_soil_mod, ONLY: sm_levels

USE jules_surface_types_mod, ONLY: ncpft

USE jules_water_resources_mod, ONLY: l_water_irrigation, nwater_use,          &
      use_irrigation

USE theta_field_sizes, ONLY: row_length=>t_i_length, rows=>t_j_length

USE um_types, ONLY: real_jlslsm

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook


IMPLICIT NONE

!-----------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  ntype
    ! Number of surface types.

!-----------------------------------------------------------------------------
! Array arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  land_index(land_pts),                                                       &
    ! Index of land points.
  plant_n_gb(land_pts)
    ! Best planting date for non-rice crops.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  dvi_cpft(land_pts,ncpft),                                                   &
    ! Development index for crop tiles.
  flandg(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),        &
    ! Land fraction.
  frac_irr_soilt(land_pts,nsoilt),                                            &
    ! Irrigation fraction.
  frac_soilt(land_pts,nsoilt),                                                &
    !  Fraction of gridbox for each soil tile.
  irrig_eff(land_pts),                                                        &
    ! Irrigation efficiency i.e. the fraction of the water withdrawn for
    ! irrigation that is demanded by the crop scheme.
  grid_area_ij(row_length,rows),                                              &
    ! Area of each gridbox (m2).
  smvccl_soilt(land_pts,nsoilt,sm_levels),                                    &
    ! Critical volumetric SMC (cubic m per cubic m of soil).
  smvcst_soilt(land_pts,nsoilt,sm_levels),                                    &
    ! Volumetric saturation point (m3/m3 of soil).
  sthf_soilt(land_pts,nsoilt,sm_levels)
    ! Frozen soil moisture content of each layer as a fraction of saturation.

!-----------------------------------------------------------------------------
! Array arguments with intent(in out)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                     &
  demand_accum(land_pts,nwater_use),                                          &
    ! Demands for water accumulated over the water resource timestep (kg).
  smcl_soilt(land_pts,nsoilt,sm_levels),                                      &
    ! Soil moisture content of each layer (kg/m2).
  sthu_irr_soilt(land_pts,nsoilt,sm_levels),                                  &
    ! Unfrozen soil moisture content of each layer as a fraction of
    ! saturation in irrigated fraction
  sthu_soilt(land_pts,nsoilt,sm_levels)
    ! Unfrozen soil moisture content of each layer as a fraction of
    ! saturation.

!-----------------------------------------------------------------------------
! Array arguments with intent(out)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  irrig_water_gb(land_pts)
    ! Water added to soil via irrigation (kg m-2 s-1).

!-----------------------------------------------------------------------------
! Local parameters
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'WATER_RESOURCES_DRIVE'

!-----------------------------------------------------------------------------
! Local scalar variables
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  i, j, l
    ! Loop counters and indices.

!-----------------------------------------------------------------------------
! Local array variables
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                     &
  demand_irrig_layer(land_pts,nsoilt,sm_levels),                              &
    ! Demand for irrigation water for each soil layer (kg m-2).
    ! This has sm_levels layers (though we only need nlayer_irrig) for ease of
    ! ensuring compatability with the code used with subroutine
    ! irrigation_demand when l_water_resources=F.
  demand_irrig_soilt(land_pts,nsoilt),                                        &
    ! Demand for water for irrigation for each soil tile (kg m-2).
  land_area(land_pts),                                                        &
    ! Area of land in gridbox (m2).
  supply_irrig(land_pts)
    ! Water supplied for irrigation (kg).

! Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------------
!end of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Calculate the demand for water for irrigation.
!-----------------------------------------------------------------------------
IF ( l_water_irrigation ) THEN

  ! Get the area of land in each land gridbox.
  DO l = 1, land_pts
    j = (land_index(l) - 1) / row_length + 1
    i = land_index(l) - (j-1) * row_length
    land_area(l) = grid_area_ij(i,j) * flandg(i,j)
  END DO

  CALL calc_irrigation_demand( plant_n_gb, dvi_cpft, frac_irr_soilt,          &
         frac_soilt, irrig_eff, land_area, smcl_soilt, smvccl_soilt,          &
         smvcst_soilt,                                                        &
         sthf_soilt, sthu_irr_soilt, sthu_soilt,                              &
         demand_accum(:,use_irrigation), demand_irrig_layer,                  &
         demand_irrig_soilt )
END IF

!----------------------------------------------------------------------------
! Here we will allocate global arrays (with size global_land_pts), then
! gather land_pts fields into fields with size global_land_pts.
! The "science" subroutines will then be run, on a single processor because
! they need knowledge of global fields.
! Finally information will then be scattered across tasks and global arrays
! deallocated.
! This code does not exist yet!
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! To allow testing during development we will assume that all irrigation
! requirements can be met. Note that the extra water demanded to account for
! the limited efficiency of irrigation will be accounted for in the conveyance
! loss - but this is not yet coded.
!-----------------------------------------------------------------------------
DO l = 1,land_pts
  supply_irrig(l) = demand_accum(l,use_irrigation) / ( 2.0 - irrig_eff(l) )
END DO

!-----------------------------------------------------------------------------
! Add irrigation to the soil store.
!-----------------------------------------------------------------------------
IF ( l_water_irrigation ) THEN
  CALL add_irrigation_to_soil( demand_accum(:,use_irrigation),                &
                               demand_irrig_soilt,                            &
                               demand_irrig_layer,                            &
                               frac_irr_soilt, land_area, smvcst_soilt,       &
                               sthf_soilt, supply_irrig,                      &
                               smcl_soilt, sthu_irr_soilt, sthu_soilt,        &
                               irrig_water_gb )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE water_resources_drive

END MODULE water_resources_drive_mod
#endif
