! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT****************************************
MODULE river_control_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RIVER_CONTROL_MOD'

CONTAINS
SUBROUTINE river_control(                                                     &
   !INTEGER, INTENT(IN)
   land_pts,                                                                  &
   !REAL, INTENT(IN)
   sub_surf_roff, surf_roff,                                                  &
   !INTEGER, INTENT(INOUT)
   a_steps_since_riv,                                                         &
   !REAL, INTENT (INOUT)
   tot_surf_runoff_gb, tot_sub_runoff_gb, acc_lake_evap_gb,                   &
   !REAL, INTENT (OUT)
   rivers_sto_per_m2_on_landpts, rflow, rrun,                                 &
   !Arguments for the UM-----------------------------------------
   !INTEGER, INTENT(IN)
   n_proc, row_length, rows, river_row_length, river_rows,                    &
   land_index, ntype, aocpl_row_length, aocpl_p_rows, g_p_field,              &
   g_r_field, global_row_length, global_rows, global_river_row_length,        &
   global_river_rows, nsurft,                                                 &
   !REAL, INTENT(IN)
   fqw_surft, delta_lambda, delta_phi, xx_cos_theta_latitude,                 &
   xpa, xua, xva, ypa, yua, yva, flandg, trivdir,                             &
   trivseq, r_area, slope, flowobs1, r_inext, r_jnext, r_land,                &
   smvcst_soilt, smvcwt_soilt, frac_surft,                                    &
   !REAL, INTENT(INOUT)
   substore, surfstore, flowin, bflowin, twatstor, smcl_soilt, sthu_soilt,    &
   !REAL, INTENT(OUT)
   inlandout_atm_gb, inlandout_riv, riverout, box_outflow, box_inflow,        &
   riverout_rgrid                                                             &
   )

!Module imports

!Common modules
USE rivers_route_mod,         ONLY: rivers_drive, scatter_land_from_riv_field
USE jules_rivers_mod,         ONLY: i_river_vn, nstep_rivers, rivers_rfm,     &
                                    rivers_trip, rivers_call, np_rivers,      &
                                    rivers_sto_rp, rivers_boxareas_rp,        &
                                    ! UM only
                                    l_inland, rivers_um_trip, rivers_speed,   &
                                    rivers_meander
USE timestep_mod,             ONLY: timestep

USE ereport_mod,              ONLY: ereport

USE atm_fields_bounds_mod,    ONLY: tdims_s, pdims_s, pdims

USE ancil_info,               ONLY: nsoilt

USE jules_soil_mod,           ONLY: sm_levels

USE overbank_update_mod,      ONLY: overbank_update

USE jules_surface_types_mod,  ONLY: lake

USE overbank_inundation_mod,  ONLY: frac_fplain_lp, frac_fplain_rp,           &
                                    l_riv_overbank

USE jules_irrig_mod,     ONLY: l_irrig_limit

!Module imports - Variables required only in UM-mode
#if defined(UM_JULES)
USE riv_intctl_mod_1A,        ONLY: riv_intctl_1a

USE timestep_mod,             ONLY: timestep_number

USE level_heights_mod,        ONLY: r_theta_levels

USE atm_land_sea_mask,        ONLY: global_land_pts => atmos_number_of_landpts
USE um_parallel_mod,          ONLY: is_master_task, gather_land_field
USE um_riv_to_jules_mod,      ONLY: um_riv_to_jules, jules_riv_to_um
#else
!Variables required only in JULES standalone-mode
!Imported module routines
USE model_grid_mod,           ONLY: global_land_pts
USE parallel_mod,             ONLY: is_master_task, gather_land_field,        &
                                    scatter_land_field
USE diag_swchs,               ONLY: srflow, srrun
#endif

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Control routine for rivers
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Subroutine arguments

INTEGER, INTENT(IN) ::                                                        &
  n_proc,                                                                     &
  land_pts,                                                                   &
  ! Number of land points
  row_length,                                                                 &
  rows,                                                                       &
  river_row_length,                                                           &
  river_rows,                                                                 &
  land_index(land_pts),                                                       &
  ntype,                                                                      &
  aocpl_row_length,                                                           &
  aocpl_p_rows,                                                               &
  g_p_field,                                                                  &
  g_r_field,                                                                  &
  global_row_length,                                                          &
  global_rows,                                                                &
  global_river_row_length,                                                    &
  global_river_rows,                                                          &
  nsurft

INTEGER, INTENT(INOUT)  ::                                                    &
  a_steps_since_riv

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  surf_roff(land_pts),                                                        &
  ! Surface runoff (kg m-2 s-1)
  sub_surf_roff(land_pts),                                                    &
  ! Sub-surface runoff (kg m-2 s-1)
  fqw_surft(land_pts,nsurft),                                                 &
  delta_lambda,                                                               &
  delta_phi,                                                                  &
  xx_cos_theta_latitude(tdims_s%i_start:tdims_s%i_end,                        &
                        tdims_s%j_start:tdims_s%j_end),                       &
  xpa(aocpl_row_length+1),                                                    &
  xua(0:aocpl_row_length),                                                    &
  xva(aocpl_row_length+1),                                                    &
  ypa(aocpl_p_rows),                                                          &
  yua(aocpl_p_rows),                                                          &
  yva(0:aocpl_p_rows),                                                        &
  flandg(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),        &
  trivdir(river_row_length, river_rows),                                      &
  trivseq(river_row_length, river_rows),                                      &
  r_area(row_length, rows),                                                   &
  slope(row_length, rows),                                                    &
  flowobs1(row_length, rows),                                                 &
  r_inext(row_length, rows),                                                  &
  r_jnext(row_length, rows),                                                  &
  r_land(row_length, rows),                                                   &
  smvcst_soilt(land_pts,nsoilt),                                              &
  smvcwt_soilt(land_pts,nsoilt),                                              &
  frac_surft(land_pts,ntype)

REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
  tot_surf_runoff_gb(land_pts),                                               &
  tot_sub_runoff_gb(land_pts),                                                &
  acc_lake_evap_gb(row_length,rows),                                          &
  twatstor(river_row_length, river_rows),                                     &
  smcl_soilt(land_pts,nsoilt,sm_levels),                                      &
  sthu_soilt(land_pts,nsoilt,sm_levels),                                      &
  substore(row_length, rows),                                                 &
  surfstore(row_length, rows),                                                &
  flowin(row_length, rows),                                                   &
  bflowin(row_length, rows)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
   rivers_sto_per_m2_on_landpts(land_pts),                                    &
   rflow(land_pts),                                                           &
   ! River flow diagnostic on land points
   rrun(land_pts),                                                            &
   ! Runoff diagnostic on land points
   inlandout_atm_gb(land_pts),                                                &
   inlandout_riv(river_row_length,river_rows),                                &
   riverout(row_length, rows),                                                &
   box_outflow(river_row_length, river_rows),                                 &
   box_inflow(river_row_length, river_rows),                                  &
   riverout_rgrid(river_row_length, river_rows)

! Local array variables.
!-----------------------------------------------------------------------------
! Jules-standalone
REAL(KIND=real_jlslsm) :: riverout_rgrid_1d(np_rivers),                       &
                          ! River outflow into the ocean on river points
                          ! Units = kg s-1. This doesn't look like it's used.
                          rivers_sto_per_m2_rgrid(np_rivers)
                          ! Water storage on river points in kg m-2

REAL(KIND=real_jlslsm), ALLOCATABLE :: global_tot_sub_runoff(:)
REAL(KIND=real_jlslsm), ALLOCATABLE :: global_tot_surf_runoff(:)
REAL(KIND=real_jlslsm), ALLOCATABLE :: global_rrun(:)
REAL(KIND=real_jlslsm), ALLOCATABLE :: global_rflow(:)

#if defined(UM_JULES)
REAL(KIND=real_jlslsm) ::                                                     &
  a_boxareas(row_length,rows),                                                &
  inlandout_atmos(row_length,rows)
#endif

!Local variables
INTEGER ::                                                                    &
  error,                                                                      &
  ! Error status from each call to ALLOCATE.
  error_sum,                                                                  &
  ! Accumulated error status.
  gather_pe_rivers,                                                           &
  l,i,j,ip

#if defined(UM_JULES)
REAL(KIND=real_jlslsm) ::                                                     &
   riv_step
   ! River timestep (secs)

LOGICAL ::                                                                    &
   first_routing,                                                             &
   invert_atmos,                                                              &
   invert_ocean = .FALSE.
#endif

!Error reporting
CHARACTER(LEN=256)       :: message
INTEGER                  :: errorstatus
CHARACTER(LEN=*), PARAMETER  :: RoutineName = 'RIVER_CONTROL'

!Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------------
!end of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

#if defined(UM_JULES)
!Pass inland flow to soil moisture every timestep
IF ( .NOT. l_inland) THEN
  DO j = 1, rows
    DO i = 1, row_length
      inlandout_atmos(i,j) = 0.0
      inlandout_riv(i,j)   = 0.0
    END DO
  END DO
END IF

!Set rrun and rflow to zero
DO l = 1, land_pts
  rrun(l)  = 0.0
  rflow(l) = 0.0
END DO

!Set the river routing to run on the 'last' PE as PE0 is very busy
gather_pe_rivers = n_proc - 1
#endif

!Initialise the accumulated surface and subsurface runoff to zero
!at the beginning of river routing timestep
IF ( a_steps_since_riv == 0 ) THEN
  tot_surf_runoff_gb = 0.0
  tot_sub_runoff_gb  = 0.0
  acc_lake_evap_gb   = 0.0
END IF

! Increment counters.
a_steps_since_riv = a_steps_since_riv + 1

IF (a_steps_since_riv == nstep_rivers) THEN
  rivers_call = .TRUE.
ELSE
  rivers_call = .FALSE.
END IF

!Accumulate the runoff as Kg/m2/s over the River Routing period
DO l = 1, land_pts
  IF (surf_roff(l) <  0.0) THEN
    !      WRITE(umMessage,*)'surf_roff(',l,')= ',surf_roff(l)
    !      CALL umPrint(umMessage,src='river_control')
  ELSE
    tot_surf_runoff_gb(l) = tot_surf_runoff_gb(l) +                           &
                         (surf_roff(l) / REAL(nstep_rivers))
  END IF
  IF (sub_surf_roff(l) <  0.0) THEN
    !      WRITE(umMessage,*)'sub_surf_roff(',l,')= ',sub_surf_roff(L)
    !      CALL umPrint(umMessage,src='river_control')
  ELSE
    tot_sub_runoff_gb(l) = tot_sub_runoff_gb(l) +                             &
                        (sub_surf_roff(l) / REAL(nstep_rivers))
  END IF
END DO

! Could this UM ifdef be replaced with the i_rivers_vn switch?
#if defined(UM_JULES)
DO l = 1, land_pts
  j = (land_index(l) - 1) / row_length +1
  i = land_index(l) - (j-1) * row_length
  acc_lake_evap_gb(i,j) = acc_lake_evap_gb(i,j) +                             &
     frac_surft(l,lake) * fqw_surft(l,lake) * timestep
END DO

!Detect first entry into river routing and initialise diagnostics
first_routing = .FALSE.
IF (timestep_number == nstep_rivers) THEN
  first_routing = .TRUE.
  riverout = 0.0
  box_outflow = 0.0
  box_inflow  = 0.0
END IF
#endif

!-------------------------------------------------------------------------
! Main call to river routing
!-------------------------------------------------------------------------
IF ( rivers_call ) THEN

#if defined(UM_JULES)
  !If ATMOS fields are as Ocean (i.e. inverted NS) set invert_atmos
  invert_atmos = .FALSE.

  IF ( .NOT. invert_ocean) THEN
    invert_atmos = .TRUE.
  END IF

  !Calculate the Atmosphere gridbox areas
  DO j = 1, rows
    DO i = 1, row_length
      a_boxareas(i,j) = r_theta_levels(i,j,0)                                 &
                    * r_theta_levels(i,j,0)                                   &
                    * delta_lambda * delta_phi                                &
                    * xx_cos_theta_latitude(i,j)
    END DO
  END DO
#endif

  SELECT CASE ( i_river_vn )
#if defined(UM_JULES)
  CASE ( rivers_um_trip )
    IF (nsoilt == 1) THEN
      riv_step = REAL(nstep_rivers) * timestep
      !This subroutine is deprecated and has not been adapted to work with
      !soil tiling.
      CALL riv_intctl_1a(                                                     &
        xpa, xua, xva, ypa, yua, yva,                                         &
        g_p_field, g_r_field, n_proc,                                         &
        gather_pe_rivers,land_pts,land_index,                                 &
        invert_atmos, row_length, rows,                                       &
        global_row_length, global_rows,                                       &
        river_row_length, river_rows,                                         &
        global_river_row_length, global_river_rows,                           &
        flandg(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
        riv_step, rivers_speed, rivers_meander,                               &
        trivdir, trivseq, twatstor, riverout_rgrid, a_boxareas,               &
        delta_phi,first_routing,                                              &
        r_area, slope, flowobs1,r_inext,r_jnext,r_land,                       &
        substore,surfstore,flowin,bflowin,                                    &
        !in/out accumulated runoff
        tot_surf_runoff_gb, tot_sub_runoff_gb,                                &
        !out
        box_outflow, box_inflow, riverout,                                    &
        !add inland basin arguments in call to rivctl
        inlandout_atmos,inlandout_riv,                                        &
        !required for soil moisture correction for water conservation
        sm_levels,acc_lake_evap_gb,smvcst_soilt,smvcwt_soilt,                 &
        smcl_soilt(1:,1,sm_levels),sthu_soilt(1:,1,sm_levels)                 &
        )
    ELSE
      errorstatus = 10
      WRITE (message,*) 'riv_intctl_1a cannot be used when nsoilt > 1'
      CALL Ereport ( RoutineName, errorstatus, message)
    END IF
#endif

  CASE ( rivers_rfm, rivers_trip )
#if defined(UM_JULES)
    ! Translate UM variables to jules_rivers_mod variables
    CALL um_riv_to_jules(g_p_field, global_river_row_length,                  &
                         global_river_rows,land_pts, land_index,              &
                         delta_phi, a_boxareas,                               &
                         r_area, flowobs1, r_inext, r_jnext, r_land,          &
                         substore, surfstore, flowin, bflowin)
#endif
    !--------------------------------------------------------------------------
    !   Gather runoff information from all processors
    !--------------------------------------------------------------------------

    IF ( is_master_task() ) THEN
      ALLOCATE(global_tot_sub_runoff(global_land_pts), stat = error)
      error_sum = error
      ALLOCATE(global_tot_surf_runoff(global_land_pts), stat = error)
      error_sum = error_sum + error
      ALLOCATE(global_rrun(global_land_pts), stat = error)
      error_sum = error_sum + error
      ALLOCATE(global_rflow(global_land_pts), stat = error)
      error_sum = error_sum + error
    ELSE
      ALLOCATE(global_tot_sub_runoff(1), stat = error)
      error_sum = error
      ALLOCATE(global_tot_surf_runoff(1), stat = error)
      error_sum = error_sum + error
      ALLOCATE(global_rrun(1), stat = error)
      error_sum = error_sum + error
      ALLOCATE(global_rflow(1), stat = error)
      error_sum = error_sum + error
    END IF

    IF ( error_sum /= 0 ) THEN
      errorstatus = 10
      CALL ereport( RoutineName, errorstatus,                                 &
                     "Error related to allocation of runoff variables." )
    END IF

    CALL gather_land_field(tot_sub_runoff_gb, global_tot_sub_runoff)
    CALL gather_land_field(tot_surf_runoff_gb, global_tot_surf_runoff)

    !-------------------------------------------------------------------------
    ! Call RFM or TRIP routing driver on single processor
    !-------------------------------------------------------------------------
    IF ( is_master_task() ) THEN

      CALL rivers_drive( global_land_pts,                                     &
                         global_tot_sub_runoff,                               &
                         global_tot_surf_runoff,                              &
                         global_rrun, global_rflow, riverout_rgrid_1d )

      !-----------------------------------------------------------------------
      ! Compute overbank inundation
      !-----------------------------------------------------------------------
      IF ( l_riv_overbank ) THEN
        CALL overbank_update()
      END IF
    END IF    ! end is_master

    !-------------------------------------------------------------------------
    ! Update output diagnostics
    !-------------------------------------------------------------------------
#if defined(UM_JULES)

    CALL jules_riv_to_um(global_rflow, riverout,                              &
                         substore, surfstore, flowin, bflowin)

#else
    IF (srflow .OR. srrun) THEN
      CALL scatter_land_field(global_rrun, rrun)
      CALL scatter_land_field(global_rflow, rflow)

      !-------------------------------------------------------------------------
      ! Update riv storage on land points used to calculate water for irrigation
      !-------------------------------------------------------------------------
      IF ( l_irrig_limit ) THEN
        DO ip = 1, np_rivers
          rivers_sto_per_m2_rgrid(ip) = rivers_sto_rp(ip) /                   &
                                        rivers_boxareas_rp(ip)
        END DO
        CALL scatter_land_from_riv_field(rivers_sto_per_m2_rgrid,             &
                                         rivers_sto_per_m2_on_landpts)
      END IF
    END IF

    !-------------------------------------------------------------------------
    ! Update fraction of inundated floodplain from overbank inundation routine
    !-------------------------------------------------------------------------
    IF ( l_riv_overbank ) THEN
      CALL scatter_land_from_riv_field( frac_fplain_rp, frac_fplain_lp )
    END IF
#endif

    DEALLOCATE(global_rflow)
    DEALLOCATE(global_rrun)
    DEALLOCATE(global_tot_surf_runoff)
    DEALLOCATE(global_tot_sub_runoff)

  CASE DEFAULT

    errorstatus = 10
    WRITE (message,'(A,I6,A)') 'River model type option ',                    &
                       i_river_vn,' not recognised.'
    CALL Ereport ( RoutineName, errorstatus, message)

  END SELECT

#if defined(UM_JULES)
  !compress inland basin outputs to land points only
  IF (l_inland) THEN
    DO l = 1,land_pts
      j = (land_index(l) - 1) / row_length +1
      i = land_index(l) - (j-1) * row_length
      inlandout_atm_gb(l) = inlandout_atmos(i,j)
    END DO
  END IF
#endif

  !-------------------------------------------------------------------------
  !   Reset counters after a call to routing.
  !-------------------------------------------------------------------------
  !Mult RIVEROUT by the number of physics timesteps per River routing
  !timestep as DAGHYD stores RIVEROUT every timestep. Non-routing
  !timestep vals are passed in as 0.0
  a_steps_since_riv = 0

END IF ! rivers_call

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE river_control

END MODULE river_control_mod
