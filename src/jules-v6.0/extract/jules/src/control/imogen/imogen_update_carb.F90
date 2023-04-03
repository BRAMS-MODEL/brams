#if !defined(UM_JULES)
!******************************COPYRIGHT**************************************
! (c) Centre for Ecology and Hydrology. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT**************************************

SUBROUTINE imogen_update_carb(progs)

USE model_time_mod, ONLY: current_time, main_run_end

USE ancil_info, ONLY: dim_cs1, dim_cslayer

USE parallel_mod, ONLY: is_master_task, gather_land_field,                    &
      scatter_land_field

USE model_grid_mod, ONLY: global_land_pts

USE ancil_info, ONLY: land_pts, nsoilt

USE jules_hydrology_mod, ONLY: l_top

USE imogen_run, ONLY:                                                         &
  include_co2,c_emissions,anom,anlg,co2_init_ppmv,land_feed_co2,              &
  nyr_emiss, ocean_feed, initial_co2_ch4_year, land_feed_ch4

USE imogen_clim, ONLY:                                                        &
  dctot_co2, dctot_ch4, lat

USE aero, ONLY: co2_mmr

USE imogen_constants, ONLY:                                                   &
  conv,ocean_area,nfarray

USE imogen_progs, ONLY:                                                       &
  co2_ppmv, co2_start_ppmv, dtemp_o, co2_change_ppmv, fa_ocean,               &
  d_land_atmos_co2, d_land_atmos_ch4, d_ocean_atmos, c_emiss_out

USE imogen_anlg_vals, ONLY:                                                   &
  t_ocean_init

USE imogen_io_vars, ONLY:                                                     &
  yr_emiss,c_emiss

USE logging_mod, ONLY: log_fatal

!TYPE definitions
USE jules_fields_mod, ONLY: toppdm, trifctltype
USE prognostics, ONLY: progs_type

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Updates the ocean and atmosphere carbon stores based on the accumulated
!   Carbon uptake of the land and prescribed emissions/concentrations.
!   Can include the natural methane feedback when land_feed_ch4 is TRUE
!
! Code Owner: Please refer to ModuleLeaders.txt
!             This file belongs in IMOGEN
!
! Code Description:
!   Language: Fortran 90.
!
!-----------------------------------------------------------------------------

! Arguments

!TYPES containing the data
TYPE(progs_type), INTENT(IN OUT) :: progs

INTEGER ::                                                                    &
  emiss_tally ! Checks that datafile of emissions includes
              ! year of interest

INTEGER :: l, n, nn, m ! loop counters

REAL ::                                                                       &
  d_land_atmos_field(land_pts)
              !d_land_atmos scattered on to each processor

REAL, ALLOCATABLE ::                                                          &
  global_data_dctot(:),                                                       &
              !  dctot on full grid
  global_data_lat(:),                                                         &
              !  lat on full grid
  global_data_tmp(:)  ! d_land_atmos on full grid
                      ! same value for each cell
!-----------------------------------------------------------------


!-----------------------------------------------------------------
! Calculate land CO2 atmosphere carbon exchange based in the difference
! between the current land carbon and the land carbon at the start of the
! year (this is rather confusingly stored as dctot_co2 in imogen_update_clim).
!-----------------------------------------------------------------
IF (land_feed_co2) THEN
  DO l = 1,land_pts
    dctot_co2(l) = trifctltype%cv_gb(l) - dctot_co2(l)
  END DO

  DO l = 1,land_pts
    DO m = 1,nsoilt
      DO n = 1,dim_cs1
        DO nn = 1,dim_cslayer
          dctot_co2(l) = dctot_co2(l) + progs%cs_pool_soilt(l,m,nn,n)
        END DO
      END DO
    END DO
  END DO

  !need to gather dctot_co2 and lat to global grid
  IF ( is_master_task() ) ALLOCATE(global_data_dctot(global_land_pts))
  IF ( is_master_task() ) ALLOCATE(global_data_lat(global_land_pts))
  IF ( is_master_task() ) ALLOCATE(global_data_tmp(global_land_pts))
  CALL gather_land_field(dctot_co2, global_data_dctot)
  CALL gather_land_field(lat, global_data_lat)

  IF ( is_master_task() ) THEN
    CALL diffcarb_land_co2(global_land_pts, d_land_atmos_co2,                 &
              global_data_lat, global_data_dctot, conv)
    global_data_tmp(:) = d_land_atmos_co2
  END IF

  !need to scatter d_land_atmos_co2 to each processor
  CALL scatter_land_field(global_data_tmp,d_land_atmos_field)
  d_land_atmos_co2 = d_land_atmos_field(1)

  IF ( ALLOCATED(global_data_dctot) ) DEALLOCATE(global_data_dctot)
  IF ( ALLOCATED(global_data_lat) ) DEALLOCATE(global_data_lat)
  IF ( ALLOCATED(global_data_tmp) ) DEALLOCATE(global_data_tmp)

END IF

!-----------------------------------------------------------------
! Now do the carbon cycling update.
!-----------------------------------------------------------------
IF (include_co2 .AND. c_emissions .AND. anom .AND. anlg) THEN

  ! Include anthropogenic carbon emissions.
  emiss_tally = 0
  DO n = 1,nyr_emiss
    IF (yr_emiss(n) == current_time%year) THEN
      c_emiss_out = c_emiss(n)
      co2_ppmv = co2_ppmv + conv * c_emiss_out
      emiss_tally = emiss_tally + 1
      ! We have found the right year so we can exit the loop
      EXIT
    END IF
  END DO

  IF (emiss_tally /= 1)                                                       &
    CALL log_fatal("imogen_update_carb",                                      &
                   'IMOGEN: Emission dataset does not match run')

  IF (land_feed_co2) THEN
    ! Update with land feedbacks if required
    co2_ppmv = co2_ppmv + d_land_atmos_co2
  END IF

  IF (ocean_feed) THEN
    ! Update with ocean feedbacks if required
    CALL ocean_co2(                                                           &
      current_time%year - initial_co2_ch4_year+1,                             &
      1,co2_ppmv,co2_init_ppmv,dtemp_o(1),                                    &
      fa_ocean,ocean_area,co2_change_ppmv,                                    &
      main_run_end%year - initial_co2_ch4_year,                               &
      t_ocean_init,nfarray,d_ocean_atmos                                      &
    )
    co2_ppmv = co2_ppmv + d_ocean_atmos
  END IF
END IF

IF (include_co2) co2_change_ppmv = co2_ppmv - co2_start_ppmv

!Unit conversion ppm to mmr (g/g): mol mass co2 / mol mass dry air * 1e-6
co2_mmr = co2_ppmv * 44.0 / 28.97 * 1.0e-6

!-----------------------------------------------------------------
! Calculate ch4 feedbacks
!-----------------------------------------------------------------
IF (land_feed_ch4) THEN
  ! This following is an attempt to future proof the format of the ch4
  ! emissions. At some point fire CH4 emissions could also be included however
  ! this will require a better treatment of the fch4_ref  and it may also be
  ! advisable to use a different logical than l_top to identify wetland
  ! emissions
  IF (l_top) THEN
    DO l = 1,land_pts
      ! Sum flux of accumulated ch4 over land
      dctot_ch4(l) = SUM(toppdm%fch4_wetl_acc_soilt(l,:)) !kg C / year
      toppdm%fch4_wetl_acc_soilt(l,:) = 0.0 !reset
    END DO
  ELSE
    CALL log_fatal("imogen_update_carb",                                      &
                   "top model is switched off (l_top) " //                    &
                   " therefore no CH4 feedbacks are included" )
  END IF

  ! Gather dctot_ch4 and lat to global grid
  IF ( is_master_task() ) ALLOCATE(global_data_dctot(global_land_pts))
  IF ( is_master_task() ) ALLOCATE(global_data_lat(global_land_pts))
  IF ( is_master_task() ) ALLOCATE(global_data_tmp(global_land_pts))
  CALL gather_land_field(dctot_ch4, global_data_dctot)
  CALL gather_land_field(lat, global_data_lat)

  ! Calculate the additional global flux of ch4 from land to atmosphere
  ! compared to a reference flux which is applicable to the atmospheric
  ! baseline ch4 concentration .
  IF ( is_master_task() ) THEN
    CALL diffcarb_land_ch4(global_land_pts, d_land_atmos_ch4,                 &
                global_data_lat, global_data_dctot)
    global_data_tmp(:) = d_land_atmos_ch4
  END IF

  ! Need to scatter d_land_atmos_ch4 to each processor
  CALL scatter_land_field(global_data_tmp, d_land_atmos_field)
  d_land_atmos_ch4 = d_land_atmos_field(1)

  ! Update atmospheric methane at end of year. This accounts for the change in
  ! methane global flux from natural source from a reference amount, and the
  ! change in decay rate due to a changed atmospheric concentration
  CALL diff_atmos_ch4(d_land_atmos_ch4, conv)

  IF ( ALLOCATED(global_data_dctot) ) DEALLOCATE(global_data_dctot)
  IF ( ALLOCATED(global_data_lat) )   DEALLOCATE(global_data_lat)
  IF ( ALLOCATED(global_data_tmp) )   DEALLOCATE(global_data_tmp)

END IF

RETURN

END SUBROUTINE imogen_update_carb
#endif
