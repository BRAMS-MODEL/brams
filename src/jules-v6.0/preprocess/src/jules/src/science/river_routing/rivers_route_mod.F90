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

MODULE rivers_route_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CONTAINS

!###############################################################################
! subroutine rivers_drive
! Driver routine for runoff routing by a kinematic wave model.

SUBROUTINE rivers_drive( global_land_pts, sub_runoffin, surf_runoffin,        &
                         runoff_out, rivflow, riverout_rgrid )
!-------------------------------------------------------------------------------
!
! Description:
!   Perform the routing of surface and sub-surface runoff defined on landpts
!
!   This routine regrids the total surface runoff to the river grid and passes
!   it to the RFM or TRIP routines to be routed.
!
!-------------------------------------------------------------------------------
! Modules used:

USE jules_rivers_mod, ONLY:                                                   &
!  imported scalars with intent(in)
     i_river_vn, rivers_rfm, rivers_trip, np_rivers

USE rivers_regrid_mod, ONLY:                                                  &
!  imported procedures
     landpts_to_rivpts, rivpts_to_landpts

USE rivers_route_rfm_mod, ONLY:                                               &
!  imported procedures
     rivers_route_rfm

USE rivers_route_trip_mod, ONLY:                                              &
!  imported procedures
     rivers_route_trip

USE jules_print_mgr, ONLY:                                                    &
  jules_message,                                                              &
  jules_print

IMPLICIT NONE

!-------------------------------------------------------------------------------
! Scalar arguments with intent(in)
INTEGER, INTENT(IN) :: global_land_pts
                             ! Size of GLOBAL runoff arrays on land points
! Array arguments with intent(in)
REAL(KIND=real_jlslsm), INTENT(IN) :: sub_runoffin(global_land_pts)
                             ! Average rate of sub surface runoff since
                             ! last rivers call on land_pts in kg m-2 s-1
REAL(KIND=real_jlslsm), INTENT(IN) :: surf_runoffin(global_land_pts)
                             ! Average rate of surface runoff since last
                             ! rivers call on land_pts in kg m-2 s-1

! Array arguments with intent(out)
REAL(KIND=real_jlslsm), INTENT(OUT) :: runoff_out(global_land_pts)
                             ! Total runoff diagnostic on land_pts
                             ! in kg m-2 s-1
REAL(KIND=real_jlslsm), INTENT(OUT) :: rivflow(global_land_pts)
                             ! River flow diagnostic on land_pts
                             ! in kg m-2 s-1
REAL(KIND=real_jlslsm), INTENT(OUT) :: riverout_rgrid(np_rivers)
                             ! River outflow into the ocean on river grid
                             ! in kg s-1

! Local scalar variables
INTEGER :: l, ip             !  loop counter

! Local array variables.
REAL(KIND=real_jlslsm) :: outflow(np_rivers)
                             !  rate of channel surface flow leaving gridbox
                             !  (kg m-2 s-1)
REAL(KIND=real_jlslsm) :: baseflow(np_rivers)
                             !  rate of channel base flow leaving
                             !  gridbox (kg m-2 s-1)
REAL(KIND=real_jlslsm) :: surf_runoff_rivers(np_rivers)
                             !  average rate of surface runoff since last
                             !  rivers call (kg m-2 s-1) on rivers grid.
REAL(KIND=real_jlslsm) :: sub_runoff_rivers(np_rivers)
                             !  average rate of sub-surface runoff since last
                             !  rivers call (kg m-2 s-1) on rivers grid.
REAL(KIND=real_jlslsm) :: tot_runoff_rivers(np_rivers)
                             ! runoff diagnostic on rivers grid (kg m-2 s-1)

!-------------------------------------------------------------------------------
! Initialisation
!-------------------------------------------------------------------------------
DO ip = 1, np_rivers
  outflow(ip)    = 0.0
  baseflow(ip)   = 0.0
  riverout_rgrid(ip) = 0.0
  sub_runoff_rivers(ip) = 0.0
  surf_runoff_rivers(ip) = 0.0
END DO

DO l = 1, global_land_pts
  rivflow(l)    = 0.0
  runoff_out(l) = 0.0
END DO

!-------------------------------------------------------------------------------
! Regrid surface and subsurface runoff from land points to rivers points
!-------------------------------------------------------------------------------
CALL landpts_to_rivpts( global_land_pts, sub_runoffin,                        &
                        np_rivers, sub_runoff_rivers )
CALL landpts_to_rivpts( global_land_pts, surf_runoffin,                       &
                        np_rivers, surf_runoff_rivers )

!-------------------------------------------------------------------------------
! Call the routing science routine.
!-------------------------------------------------------------------------------
SELECT CASE ( i_river_vn )

CASE ( rivers_rfm )
  CALL rivers_route_rfm( surf_runoff_rivers, sub_runoff_rivers,               &
                         outflow, baseflow, riverout_rgrid )
CASE ( rivers_trip )
  CALL rivers_route_trip( surf_runoff_rivers, sub_runoff_rivers,              &
                          outflow, baseflow, riverout_rgrid )
CASE DEFAULT
  WRITE(jules_message,*) 'ERROR: rivers_drive: ' //                           &
                         'do not recognise i_river_vn=', i_river_vn
  CALL jules_print('rivers_route_drive',jules_message)
END SELECT

!-------------------------------------------------------------------------------
!   Update diagnostics
!------------------------------------------------------------------------------
! Sum the surface and subsurface runoffs
DO ip = 1, np_rivers
  tot_runoff_rivers(ip) = surf_runoff_rivers(ip) + sub_runoff_rivers(ip)
END DO

! Compute river storage per m2 (used with irrigation scheme)
!rivers_sto_per_m2_rgrid(:) = rivers_sto_rp(:) / rivers_boxareas_rp(:)

!-------------------------------------------------------------------------------
!   Regrid from rivers to land grid
!------------------------------------------------------------------------------
CALL rivpts_to_landpts( np_rivers, outflow, global_land_pts, rivflow )
CALL rivpts_to_landpts( np_rivers, tot_runoff_rivers,                         &
                        global_land_pts, runoff_out )

END SUBROUTINE rivers_drive

!###############################################################################

SUBROUTINE scatter_land_from_riv_field( global_var_riv, var_land )
!-----------------------------------------------------------------------------
! Description:
!   Uses to the river variable on global river points array (e.g. a prognostic)
!   to fill a variable on land points array, scattered across all processors
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

USE jules_rivers_mod, ONLY:                                                   &
!  imported scalars with intent(in)
    np_rivers
USE rivers_regrid_mod, ONLY: rivpts_to_landpts

USE parallel_mod, ONLY: is_master_task, scatter_land_field
USE model_grid_mod, ONLY: global_land_pts

USE ancil_info, ONLY: land_pts
USE ereport_mod, ONLY: ereport

IMPLICIT NONE

! Array arguments with intent(in)
REAL(KIND=real_jlslsm), INTENT(IN) :: global_var_riv(np_rivers)
                             ! Input variable on riverpoints
! Array arguments with intent(in)
REAL(KIND=real_jlslsm), INTENT(OUT) :: var_land(land_pts)
                             ! Input variable on riverpoints
! Local arrays
REAL(KIND=real_jlslsm), ALLOCATABLE :: global_rivers_var_on_landpts(:)
          ! River routing output variable defined on global model land points
! Local variables
INTEGER:: error, errorstatus    ! Error status from ALLOCATE call
INTEGER:: l                     ! Loop counter

!-------------------------------------------------------------------------------
! Initialise
DO l = 1, land_pts
  var_land(l) = 0.0
END DO

!-------------------------------------------------------------------------------
! Allocate global variable on land points
IF ( is_master_task() ) THEN
  ALLOCATE(global_rivers_var_on_landpts(global_land_pts), stat = error)
  IF ( error == 0 ) THEN
    DO l = 1, global_land_pts
      global_rivers_var_on_landpts(l) = 0.0
    END DO
  END IF
ELSE
  ALLOCATE(global_rivers_var_on_landpts(1), stat = error)
  IF ( error == 0 ) THEN
    global_rivers_var_on_landpts = 0.0
  END IF
END IF

IF ( error /= 0 ) THEN
  errorstatus = 10
  CALL ereport( 'scatter_land_from_riv_field', errorstatus,                   &
                "Error related to allocation in scatter_land_from_riv." )
END IF


!-------------------------------------------------------------------------------
! Regrid/remap variable from rivers to land grid
IF ( is_master_task() ) THEN

  CALL rivpts_to_landpts( np_rivers, global_var_riv,                          &
                          global_land_pts, global_rivers_var_on_landpts )

END IF

!-------------------------------------------------------------------------------
! Scatter land variable onto separate processors
CALL scatter_land_field(global_rivers_var_on_landpts, var_land)

DEALLOCATE(global_rivers_var_on_landpts)

END SUBROUTINE scatter_land_from_riv_field

!###############################################################################

SUBROUTINE adjust_routestore
!-----------------------------------------------------------------------------
! Description:
!   Uses to the river adjustment factor on land points (which accounts for
!   water extracted from the river storage for irrigation)
!   to adjust the river storage on river points array (a prognostic).
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

USE jules_rivers_mod, ONLY:                                                   &
!  imported scalars with intent(in)
    np_rivers                                                                 &
!  imported array with intent(in)
    ,rivers_sto_rp, rivers_adj_on_landpts

USE rivers_regrid_mod, ONLY: landpts_to_rivpts

USE parallel_mod, ONLY: is_master_task, gather_land_field
USE model_grid_mod, ONLY: global_land_pts

IMPLICIT NONE

REAL(KIND=real_jlslsm), ALLOCATABLE :: global_rivers_adj_on_landpts(:)
          ! rivers adjustment factor on global model land points

REAL(KIND=real_jlslsm) :: rivers_adj_rgrid(np_rivers)
          ! rivers adjustment factor on routing grid
INTEGER :: ip ! array indices

IF ( is_master_task() ) THEN
  ALLOCATE(global_rivers_adj_on_landpts(global_land_pts))
ELSE
  ALLOCATE(global_rivers_adj_on_landpts(1))
END IF

CALL gather_land_field(rivers_adj_on_landpts, global_rivers_adj_on_landpts)

IF ( is_master_task() ) THEN

  DO ip = 1, np_rivers
    rivers_adj_rgrid(ip) = 1.0
  END DO

  !-------------------------------------------------------------------------------
  !   Regrid from land to rivers grid
  !-------------------------------------------------------------------------------
  CALL landpts_to_rivpts( global_land_pts, global_rivers_adj_on_landpts,      &
                          np_rivers, rivers_adj_rgrid )

  !-------------------------------------------------------------------------------
  ! Apply the correction factor to the route storage, which is a prognostic
  !-------------------------------------------------------------------------------
  DO ip = 1,np_rivers
    rivers_sto_rp(ip) = rivers_adj_rgrid(ip) * rivers_sto_rp(ip)
  END DO

END IF

DEALLOCATE(global_rivers_adj_on_landpts)

END SUBROUTINE adjust_routestore

!###############################################################################

END MODULE rivers_route_mod
