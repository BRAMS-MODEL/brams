#if defined(UM_JULES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!     Driver and science routines for calculating river flow routing
!     using a kinematic wave model
!     RFM: see Bell et al. 2007 Hydrol. Earth Sys. Sci. 11. 532-549
!     TRIP: see Oki et al 1999 J.Met.Soc.Japan, 77, 235-255.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

MODULE um_riv_to_jules_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CONTAINS

!###############################################################################
! subroutine um_riv_to_jules
! Conversion routine for UM variables to JULES runoff routing

SUBROUTINE um_riv_to_jules(g_p_field, global_river_row_length,                &
                           global_river_rows, land_points, land_index,        &
                           delta_phi, a_boxareas,                             &
                           r_area, flowobs1, r_inext,r_jnext,r_land,          &
                           substore,surfstore,flowin,bflowin)

USE UM_ParVars,        ONLY: glsize
USE Field_Types,       ONLY: fld_type_p

USE ancil_info,        ONLY: land_pts
USE theta_field_sizes, ONLY: t_i_length, t_j_length

USE atm_land_sea_mask, ONLY: global_land_pts => atmos_number_of_landpts
USE um_parallel_mod,   ONLY: is_master_task, um_gather_field
USE um_latlon_mod,     ONLY: global_land_index
USE jules_rivers_mod,  ONLY: i_river_vn

USE jules_rivers_mod,  ONLY:                                                  &
  !  imported scalar parameters
  np_rivers, rivers_dlat, rivers_first, rivers_dx                             &
  ,i_river_vn, rivers_trip, rivers_rfm, rivers_um_trip                        &
  ,rivers_regrid                                                              &
  !  imported arrays with intent(in)
  ,rfm_flowobs1_rp, il_river_grid                                             &
  !  imported arrays with intent(inout)
  ,rivers_next_rp, rfm_iarea_rp, rfm_land_rp, rivers_sto_rp                   &
  ,rfm_substore_rp, rfm_surfstore_rp, rfm_flowin_rp, rfm_bflowin_rp           &
  ,rfm_rivflow_rp, rfm_baseflow_rp, rivers_boxareas_rp

USE planet_constants_mod,     ONLY: planet_radius

USE ereport_mod,              ONLY:                                           &
  ereport

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                        &
  land_points                                                                 &
                ! number of landpoints
  ,land_index (land_points)                                                   &
                ! index of land to global points
  ,g_p_field                                                                  &
                ! size of global ATMOS field
  ,global_river_row_length                                                    &
                ! size of global RIV field
  ,global_river_rows
                ! size of global RIV field

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  delta_phi
                ! RCM gridsize (radians)

! ancillary variables for river routing model
! should be defined on river_row_length, river_rows
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  r_area(t_i_length,t_j_length)                                               &
                ! accumulated areas file
  ,r_inext(t_i_length,t_j_length)                                             &
                ! x-coordinate of downstream grid pt
  ,r_jnext(t_i_length,t_j_length)                                             &
                ! y-coordinate of downstream grid pt
  ,flowobs1(t_i_length,t_j_length)                                            &
                ! optional initialisation for flows
  ,r_land(t_i_length,t_j_length)                                              &
                ! land/river depends on value of a_thresh
  ,a_boxareas(t_i_length,t_j_length)
                ! ATMOS gridbox areas

! prognostic variables for river routing model
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  substore(t_i_length,t_j_length)                                             &
                ! routing sub_surface store (mm)
  ,surfstore(t_i_length,t_j_length)                                           &
                ! routing surface store (mm)
  ,flowin(t_i_length,t_j_length)                                              &
                ! surface lateral inflow (mm)
  ,bflowin(t_i_length,t_j_length)
                ! sub-surface lateral inflow (mm)

! local variables
INTEGER :: i,j,l,gf,ip
INTEGER :: inext, jnext, global_row_length
INTEGER :: error, error_sum, errcode

INTEGER, ALLOCATABLE :: mapfr(:,:)                   ! map full to river

! Gathering and Scattering variables:
REAL(KIND=real_jlslsm) :: gather_riverout_ATMOS(g_p_field)
                ! river outflow at seapoints on the ATMOS grid
REAL(KIND=real_jlslsm) :: gather_r_area(g_p_field)
                ! global field of accumulated area
REAL(KIND=real_jlslsm) :: gather_r_inext(g_p_field)
                ! global field of x-flow directions
REAL(KIND=real_jlslsm) :: gather_r_jnext(g_p_field)
                ! global field of y-flow directions
REAL(KIND=real_jlslsm) :: gather_slope(g_p_field)
                ! global field of slope
REAL(KIND=real_jlslsm) :: gather_flowobs1(g_p_field)
                ! global field initial flow values
REAL(KIND=real_jlslsm) :: gather_r_land(g_p_field)
                ! global field of land-type
REAL(KIND=real_jlslsm) :: gather_substore(g_p_field)
                ! global field of surface storage
REAL(KIND=real_jlslsm) :: gather_surfstore(g_p_field)
                ! global field sub-surface storage
REAL(KIND=real_jlslsm) :: gather_flowin(g_p_field)
                ! global field of flowin
REAL(KIND=real_jlslsm) :: gather_bflowin(g_p_field)
                ! global field of bflowin
REAL(KIND=real_jlslsm) :: gather_flandg(g_p_field)
                ! field for gather to land/sea mask
REAL(KIND=real_jlslsm) :: gather_boxareas(g_p_field)
                ! field for gather to grid box areas

! 2.  Initialise river routing parameters
!-----------------------------------------------------------------------------
! River routing algorithm
! Lat/lon grid spacing
rivers_dlat = delta_phi

! Regular grid spacing (m)
! rivers_dx is read from namelist jules_rivers_props in standalone mode;
! in UM it will have to be properly calculated for RFM
! The value of 1500 is what is needed for UKV
rivers_dx = planet_radius * delta_phi

! Flag to specify if regridding required from land to river grids
rivers_regrid = .FALSE.

! 3.   Compute number of valid river points
!-----------------------------------------------------------------------------
np_rivers = global_land_pts     ! assume 1:1 if not regrid

! 4.   Allocate river routing variables
!-----------------------------------------------------------------------------
IF ( rivers_first ) THEN
  IF ( is_master_task() ) THEN

    ALLOCATE( mapfr(global_river_row_length,global_river_rows), stat = error )
    error_sum = error
    ALLOCATE( rfm_flowobs1_rp(np_rivers),    stat = error )
    error_sum = error
    ALLOCATE( rivers_next_rp(np_rivers),     stat = error )
    error_sum = error_sum + error
    ALLOCATE( rfm_iarea_rp(np_rivers),       stat = error )
    error_sum = error_sum +error
    ALLOCATE( rfm_land_rp(np_rivers),        stat = error )
    error_sum = error_sum + error
    ALLOCATE( rivers_sto_rp(np_rivers),      stat = error )
    error_sum = error_sum + error
    ALLOCATE( rivers_boxareas_rp(np_rivers), stat = error )
    error_sum = error_sum + error
    ALLOCATE( rfm_rivflow_rp(np_rivers),     stat = error )
    error_sum = error_sum + error
    ALLOCATE( il_river_grid(np_rivers),      stat = error )
    error_sum = error_sum + error
    ALLOCATE( rfm_surfstore_rp(np_rivers),   stat = error )
    error_sum = error_sum + error
    ALLOCATE( rfm_substore_rp(np_rivers),    stat = error )
    error_sum = error_sum + error
    ALLOCATE( rfm_flowin_rp(np_rivers),      stat = error )
    error_sum = error_sum + error
    ALLOCATE( rfm_bflowin_rp(np_rivers),     stat = error )
    error_sum = error_sum + error
    ALLOCATE( rfm_baseflow_rp(np_rivers),    stat = error )
    error_sum = error_sum + error

    IF ( error_sum /= 0 )                                                     &
      errcode = 150
    CALL ereport("um_riv_to_jules", errcode,                                  &
                 "Error allocating JULES model arrays")

    DO ip = 1, np_rivers
      rfm_land_rp(ip)      = 2       ! Initialise all points to land
      rfm_iarea_rp(ip)     = 0

      rfm_flowobs1_rp(ip)  = 0.0
      rfm_surfstore_rp(ip) = 0.0
      rfm_substore_rp(ip)  = 0.0
      rfm_flowin_rp(ip)    = 0.0
      rfm_bflowin_rp(ip)   = 0.0
      rfm_rivflow_rp(ip)   = 0.0
      rfm_baseflow_rp(ip)  = 0.0
      rivers_next_rp(ip)   = 0
      il_river_grid(ip)    = 0
      rivers_sto_rp(ip)    = 0.0
    END DO

    DO j = 1, global_river_rows
      DO i = 1, global_river_row_length
        mapfr(i,j) = 0
      END DO
    END DO

  END IF       !! end master_task
END IF    !! end rivers_first

! 5.   Gather river routing ancillary and prognostic variables
!-----------------------------------------------------------------------------
CALL um_gather_field(r_area,     gather_r_area)
CALL um_gather_field(flowobs1,   gather_flowobs1)
CALL um_gather_field(r_inext,    gather_r_inext)
CALL um_gather_field(r_jnext,    gather_r_jnext)
CALL um_gather_field(r_land,     gather_r_land)
CALL um_gather_field(a_boxareas, gather_boxareas)
CALL um_gather_field(substore,   gather_substore)
CALL um_gather_field(surfstore,  gather_surfstore)
CALL um_gather_field(flowin,     gather_flowin)
CALL um_gather_field(bflowin,    gather_bflowin)

! 6.   Translate gathered variables to river points vector variables
!-----------------------------------------------------------------------------
IF ( rivers_first ) THEN
  IF ( is_master_task() ) THEN

    global_row_length = glsize(1,fld_type_p)
    DO l = 1,global_land_pts
      j = (global_land_index(l) - 1) / global_row_length + 1
      i = global_land_index(l) - (j-1) * global_row_length
      gf = i + (j-1) * global_row_length

      rfm_iarea_rp(l) = NINT(gather_r_area(gf))
      rfm_land_rp(l) = NINT(gather_r_land(gf))
      rfm_substore_rp(l) = gather_substore(gf)
      rivers_boxareas_rp(l) = gather_boxareas(gf)
      rfm_surfstore_rp(l) = gather_surfstore(gf)
      rfm_flowin_rp(l) = gather_flowin(gf)
      rfm_bflowin_rp(l) = gather_bflowin(gf)
      rfm_flowobs1_rp(l) = gather_flowobs1(gf)

      mapfr(i,j) = l      ! full 2D rivgrid to rivpts
    END DO

    DO l = 1,global_land_pts
      j = (global_land_index(l) - 1) / global_row_length + 1
      i = global_land_index(l) - (j-1) * global_row_length
      gf = i + (j-1) * global_row_length

      IF ( gather_r_inext(gf) >= -1 .AND. gather_r_jnext(gf) >= -1 ) THEN
        inext =  i + NINT(gather_r_inext(gf))
        jnext =  j + NINT(gather_r_jnext(gf))

        ! cyclic bcs for rivers_reglatlon    !! check this (cf init_river_props)
        IF ( inext > global_river_row_length ) inext = 1
        IF ( inext < 1 ) inext = i !inext + global_river_row_length

        IF ( jnext > global_river_rows ) jnext = 1
        IF ( jnext < 1 ) jnext = j !jnext + global_river_rows

        IF ( mapfr(inext, jnext) > 0 ) THEN
          rivers_next_rp(l) = mapfr(inext,jnext)
        END IF
      END IF
    END DO

    DEALLOCATE(mapfr)

  END IF    !! master task
END IF    !! rivers_first

END SUBROUTINE um_riv_to_jules

!###############################################################################
! subroutine jules_riv_to_um
! Conversion routine for JULES runoff routing variables to UM variables

SUBROUTINE jules_riv_to_um(global_rflow, riverout_atmos,                      &
                           substore, surfstore, flowin, bflowin)

USE rivers_regrid_mod, ONLY: rivpts_to_landpts
USE atm_land_sea_mask, ONLY: global_land_pts => atmos_number_of_landpts
USE theta_field_sizes, ONLY: t_i_length, t_j_length
USE um_parallel_mod,   ONLY: is_master_task, scatter_land2d_field

USE jules_rivers_mod,  ONLY:                                                  &
  !  imported scalar parameters
  np_rivers                                                                   &
  !  imported arrays with intent(inout)
  ,rivers_sto_rp, rfm_substore_rp, rfm_surfstore_rp, rfm_flowin_rp            &
  ,rfm_bflowin_rp, rfm_rivflow_rp, rfm_baseflow_rp

IMPLICIT NONE

REAL(KIND=real_jlslsm), INTENT(IN) :: global_rflow(:)
                                           ! River flow diagnostic

! prognostic variables for river routing model
REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
  riverout_atmos(t_i_length, t_j_length)                                      &
                ! ROUTING RIVER FLOW on atmosphere grid (kg/m2/s)
  ,substore(t_i_length,t_j_length)                                            &
                ! ROUTING SUB_SURFACE STORE (MM)
  ,surfstore(t_i_length,t_j_length)                                           &
                ! ROUTING SURFACE STORE (MM)
  ,flowin(t_i_length,t_j_length)                                              &
                !SURFACE LATERAL INFLOW (MM)
  ,bflowin(t_i_length,t_j_length)
                ! SUB-SURFACE LATERAL INFLOW (MM)

REAL(KIND=real_jlslsm) ::                                                     &
  substore_lp(global_land_pts)                                                &
                ! ROUTING SUB_SURFACE STORE (MM)
  ,surfstore_lp(global_land_pts)                                              &
                ! ROUTING SURFACE STORE (MM)
  ,flowin_lp(global_land_pts)                                                 &
                !SURFACE LATERAL INFLOW (MM)
  ,bflowin_lp(global_land_pts)
                ! SUB-SURFACE LATERAL INFLOW (MM)

INTEGER :: l    ! Loop counter

! Initialization
DO l = 1, global_land_pts
  flowin_lp = 0.0
  bflowin_lp = 0.0
  surfstore_lp = 0.0
  substore_lp = 0.0
END DO

IF ( is_master_task() ) THEN
  CALL rivpts_to_landpts( np_rivers, rfm_flowin_rp, global_land_pts, flowin_lp )
  CALL rivpts_to_landpts( np_rivers, rfm_bflowin_rp, global_land_pts, bflowin_lp )
  CALL rivpts_to_landpts( np_rivers, rfm_surfstore_rp, global_land_pts, surfstore_lp )
  CALL rivpts_to_landpts( np_rivers, rfm_substore_rp, global_land_pts, substore_lp )
END IF

! Update river flow diagnostics for UM dump
CALL scatter_land2d_field(global_rflow, riverout_atmos)

! Update prognostics for UM dump
CALL scatter_land2d_field(flowin_lp, flowin)
CALL scatter_land2d_field(bflowin_lp, bflowin)
CALL scatter_land2d_field(surfstore_lp, surfstore)
CALL scatter_land2d_field(substore_lp, substore)

END SUBROUTINE jules_riv_to_um

END MODULE um_riv_to_jules_mod
#endif
