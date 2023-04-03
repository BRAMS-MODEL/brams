!###############################################################################
!###############################################################################

MODULE rivers_regrid_mod

!-----------------------------------------------------------------------------
! Description:
!   Contains river routing regridding/remapping functions for standalone
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CONTAINS

!###############################################################################
! subroutine rivers_get_xy_pos
!
! Internal procedure in module grid_utils
!    Given a point number and the extents (size) of a 2-D (xy) grid, returns
!    the x and y indices (coords) of the point. All coords are relative to (1,1)
!    at bottom left of grid, with numbering running left to right, bottom to
!    top.

SUBROUTINE rivers_get_xy_pos( i, nx, ny, ix, iy )

IMPLICIT NONE

! Scalar arguments with intent (in)
INTEGER, INTENT(IN) ::                                                        &
                       i,                                                     &
                   !  the point number
                       nx,                                                    &
                   !  the extent (size) of the grid in the x direction
                       ny
                   !  the extent (size) of the grid in the y direction
                   !  NB ny is only required for (cursory) error checking.

! Scalars with intent (out)
INTEGER, INTENT(OUT) ::  ix,iy   !  the x and y coordinates of the point

! Assume that i>=1!
iy = ( i - 1 ) / nx + 1
ix = i - (iy-1) * nx

! If locations are out of range (grid is too small for this point number),
! set to -1.
IF ( ix > nx .OR. iy > ny ) THEN
  ix = -1
  iy = -1
END IF

END SUBROUTINE rivers_get_xy_pos
!###############################################################################


!###############################################################################
! subroutine rivers_remap_match
! Computes regridding between land points and main model grid where both
! grids overlap (i.e. rivers_regrid=False)
SUBROUTINE rivers_remap_match
!-----------------------------------------------------------------------------
!
! Description:
!   Computes remapping between land points and river routing grid
!
!   If the "main" model grid is 2-D, this target grid is the 2-D grid.
!   If the "main" grid is a vector (in offline applications of JULES this is
!   possible if points from a larger grid have been compressed - e.g. land
!   points selected from a larger grid.), the target grid is the larger grid,
!   across which the points are to be scattered.
!
!------------------------------------------------------------------------------
! Modules used:

USE jules_rivers_mod, ONLY:                                                   &
!  imported arrays with intent(in)
     rivers_lat_rp, rivers_lon_rp, np_rivers, ir_land_grid, il_river_grid

#if defined(UM_JULES)
USE atm_land_sea_mask, ONLY: global_land_pts => atmos_number_of_landpts
USE um_latlon_mod, ONLY: latitude_of_land_pts => latitude,                    &
                         longitude_of_land_pts => longitude
USE um_parallel_mod, ONLY: is_master_task,                                    &
                           gather_land_field => gather_land2d_field
#else
USE model_grid_mod, ONLY: global_land_pts, latitude_of_land_pts,              &
    longitude_of_land_pts
USE parallel_mod, ONLY: is_master_task, gather_land_field
#endif

IMPLICIT NONE

REAL(KIND=real_jlslsm), ALLOCATABLE :: global_lat_of_land_pts(:)
REAL(KIND=real_jlslsm), ALLOCATABLE :: global_lon_of_land_pts(:)

! Local scalar variables.

INTEGER :: l      !  loop counter (land point)
INTEGER :: ip     !  work

!------------------------------------------------------------------------------
! Compute full land_pts grid lat/lon
!------------------------------------------------------------------------------

IF ( is_master_task() ) THEN
  ALLOCATE(global_lat_of_land_pts(global_land_pts))
  ALLOCATE(global_lon_of_land_pts(global_land_pts))
ELSE
  ALLOCATE(global_lat_of_land_pts(1))
  ALLOCATE(global_lon_of_land_pts(1))
END IF

global_lat_of_land_pts(:) = 0.0
global_lon_of_land_pts(:) = 0.0

CALL gather_land_field(latitude_of_land_pts, global_lat_of_land_pts)
CALL gather_land_field(longitude_of_land_pts, global_lon_of_land_pts)

!------------------------------------------------------------------------------
! Compute remap between land points and full river grid
!         general case, covers all grid types when rivers_regrid=False
!------------------------------------------------------------------------------

IF ( is_master_task() ) THEN
  DO l = 1,global_land_pts
    DO ip = 1,np_rivers
      IF (rivers_lat_rp(ip) == global_lat_of_land_pts(l) .AND.                &
         rivers_lon_rp(ip) == global_lon_of_land_pts(l)) THEN
        il_river_grid(ip) = l
        ir_land_grid(l) = ip
        EXIT
      END IF
    END DO
  END DO
END IF

DEALLOCATE(global_lat_of_land_pts)
DEALLOCATE(global_lon_of_land_pts)

END SUBROUTINE rivers_remap_match
!###############################################################################
!###############################################################################
! subroutine rivers_remap_unmatch
! Setup initial arrays for regridding between land points and main model grid
! where both grids DO NOT overlap
SUBROUTINE rivers_remap_unmatch
!-----------------------------------------------------------------------------
!
! Description:
!
!------------------------------------------------------------------------------
! Modules used:

USE jules_rivers_mod, ONLY:                                                   &
!  imported scalars with intent(in) - defining grid
     nx_grid, reg_dlat, reg_dlon, reg_lon1, reg_lat1,                         &
     rivers_reglatlon, rivers_regrid,                                         &
!  imported arrays with intent(in)
     ir_land_grid

#if defined(UM_JULES)
USE atm_land_sea_mask, ONLY: global_land_pts => atmos_number_of_landpts
USE um_latlon_mod, ONLY: latitude_of_land_pts => latitude,                    &
                         longitude_of_land_pts => longitude
USE um_parallel_mod, ONLY: is_master_task,                                    &
                           gather_land_field => gather_land2d_field
#else
USE model_grid_mod, ONLY: global_land_pts, latitude_of_land_pts,              &
    longitude_of_land_pts
USE parallel_mod, ONLY: is_master_task, gather_land_field
#endif

IMPLICIT NONE

REAL(KIND=real_jlslsm), ALLOCATABLE :: global_lat_of_land_pts(:)
REAL(KIND=real_jlslsm), ALLOCATABLE :: global_lon_of_land_pts(:)

! Local scalar variables.

INTEGER :: l      !  loop counter (land point)
INTEGER :: i,j    !  work

!------------------------------------------------------------------------------
! Compute full land_pts grid lat/lon
!------------------------------------------------------------------------------

IF ( is_master_task() ) THEN
  ALLOCATE(global_lat_of_land_pts(global_land_pts))
  ALLOCATE(global_lon_of_land_pts(global_land_pts))
ELSE
  ALLOCATE(global_lat_of_land_pts(1))
  ALLOCATE(global_lon_of_land_pts(1))
END IF

global_lat_of_land_pts(:) = 0.0
global_lon_of_land_pts(:) = 0.0

CALL gather_land_field(latitude_of_land_pts, global_lat_of_land_pts)
CALL gather_land_field(longitude_of_land_pts, global_lon_of_land_pts)

!------------------------------------------------------------------------------
! Precalculate landpoints mapping to full 2D grid if river grid and land
! grids do not coincide
!------------------------------------------------------------------------------
IF ( is_master_task() ) THEN
  IF ( rivers_regrid .AND. rivers_reglatlon ) THEN

    DO l = 1,global_land_pts
      i = NINT( (global_lon_of_land_pts(l) - reg_lon1) / reg_dlon ) + 1
      j = NINT( (global_lat_of_land_pts(l) - reg_lat1) / reg_dlat ) + 1
      ir_land_grid(l) = (j-1) * nx_grid + i
    END DO

  END IF
END IF

DEALLOCATE(global_lat_of_land_pts)
DEALLOCATE(global_lon_of_land_pts)

END SUBROUTINE rivers_remap_unmatch
!###############################################################################

!###############################################################################
! subroutine landpts_to_rivpts
! Handles regridding/remapping/no change of variables from landpts to rivpts

SUBROUTINE landpts_to_rivpts( land_pts, var_land, riv_pts, var_riv )

!------------------------------------------------------------------------------
! Description:
!   Checks what method of conversion required between land_pts and riv_pts, and
!   invokes regridding/remapping/no change for a given input variable.

USE jules_rivers_mod, ONLY:                                                   &
     rivers_regrid, il_river_grid

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts  !  a vector length
INTEGER, INTENT(IN) :: riv_pts  !  a vector length

! Array arguments with intent(in)
REAL(KIND=real_jlslsm), INTENT(IN) :: var_land(land_pts)
                                          ! input variable on land_pts

! Array arguments with intent(out)
REAL(KIND=real_jlslsm), INTENT(OUT) :: var_riv(riv_pts)
                                          ! output variable on riv_pts

! Local counter
INTEGER :: ip

! If regridding required, call routine
IF ( rivers_regrid ) THEN

  CALL rivers_regrid_from_land( land_pts, var_land, riv_pts, var_riv )

  ! If not regridding, translate between global_land and river point vectors
ELSE IF (land_pts /= riv_pts) THEN

  DO ip = 1,riv_pts
    IF (il_river_grid(ip) > 0) THEN
      var_riv(ip) = var_land(il_river_grid(ip))
    END IF
  END DO

  ! If grids identical, including land/riv points in same order, no regrid
ELSE

  DO ip = 1,riv_pts
    var_riv(ip)  = var_land(ip)
  END DO

END IF

END SUBROUTINE landpts_to_rivpts

!###############################################################################
! subroutine rivpts_to_landpts
! Handles regridding/remapping/no change of variables from landpts to rivpts

SUBROUTINE rivpts_to_landpts( riv_pts, var_riv, land_pts, var_land )

!------------------------------------------------------------------------------
! Description:
!   Checks what method of conversion required between land_pts and riv_pts, and
!   invokes regridding/remapping/no change for a given input variable.

USE jules_rivers_mod, ONLY:                                                   &
     rivers_regrid, il_river_grid

IMPLICIT NONE

INTEGER, INTENT(IN) :: riv_pts  !  a vector length
INTEGER, INTENT(IN) :: land_pts  !  a vector length

! Array arguments with intent(in)
REAL(KIND=real_jlslsm), INTENT(IN) :: var_riv(riv_pts)
                                         ! input variable on riv_pts

! Array arguments with intent(out)
REAL(KIND=real_jlslsm), INTENT(OUT) :: var_land(land_pts)
                                         ! output variable on land_pts

! Local counter
INTEGER :: ip

! If regridding required, call routine
IF ( rivers_regrid ) THEN

  CALL rivers_regrid_to_land( riv_pts, var_riv, land_pts, var_land )

  ! If not regridding, translate between global_land and river point vectors
ELSE IF (land_pts /= riv_pts) THEN

  DO ip = 1,riv_pts
    IF (il_river_grid(ip) > 0) THEN
      var_land(il_river_grid(ip)) = var_riv(ip)
    END IF
  END DO

  ! If grids identical, including land/riv points in same order, no regrid
ELSE

  DO ip = 1,riv_pts
    var_land(ip) = var_riv(ip)
  END DO

END IF

END SUBROUTINE rivpts_to_landpts

!###############################################################################
! subroutine rivers_regrid_from_land
! Handles regridding of runoff from "main" to rivers grids.

SUBROUTINE rivers_regrid_from_land( land_pts, runoff_land, riv_pts,           &
                                    runoff_riv )
!------------------------------------------------------------------------------
! Description:
!   Regrids runoff from a source grid to a target grid (the rivers grid).
!   Both grids must be regular in latitude and longitude.

USE jules_rivers_mod, ONLY:                                                   &
!  imported scalars with intent(in)
     nx_rivers,ny_rivers, nx_grid, ny_grid, rivers_index_rp, ir_land_grid

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts  !  a vector length
INTEGER, INTENT(IN) :: riv_pts  !  a vector length

! Array arguments with intent(in)

REAL(KIND=real_jlslsm), INTENT(IN) :: runoff_land(land_pts)
!  runoff rate on model grid (kg m-2 s-1)

! Array arguments with intent(out)
REAL(KIND=real_jlslsm), INTENT(OUT) :: runoff_riv(riv_pts)
!  runoff rate on rivers grid (kg m-2 s-1)

INTEGER :: maxlin, l, ip, ix, iy

REAL(KIND=real_jlslsm) :: runoff_grid(nx_grid,ny_grid)
REAL(KIND=real_jlslsm) :: runoff_out(nx_rivers,ny_rivers)

! Convert from land_pts to full grid
runoff_grid(:,:) = 0.0
DO l = 1,land_pts
  CALL rivers_get_xy_pos( ir_land_grid(l), nx_grid, ny_grid, ix, iy )
  runoff_grid(ix,iy) = runoff_land(l)
END DO

! Call to lat/lon grid-based regridding routine
maxlin = ( nx_grid + nx_rivers ) * ( ny_grid + ny_rivers )
CALL rivers_route_regrid ( maxlin, runoff_grid, runoff_out )

! Convert from rivers grid to riv_pts
DO ip = 1,riv_pts
  CALL rivers_get_xy_pos(rivers_index_rp(ip),nx_rivers,ny_rivers,ix,iy)
  runoff_riv(ip) = runoff_out(ix,iy)
END DO

END SUBROUTINE rivers_regrid_from_land

!###############################################################################
! subroutine rivers_regrid_to_land
! Handles regridding of runoff from rivers to "main" grids.

SUBROUTINE rivers_regrid_to_land( riv_pts, runoff_riv, land_pts, runoff_land )
!------------------------------------------------------------------------------
! Description:
!   Regrids runoff from a source grid to a target grid (the rivers grid).
!   Both grids must be regular in latitude and longitude.

USE jules_rivers_mod, ONLY:                                                   &
!  imported scalars with intent(in)
     nx_rivers,ny_rivers, nx_grid, ny_grid, rivers_index_rp, ir_land_grid

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts  !  a vector length
INTEGER, INTENT(IN) :: riv_pts  !  a vector length

! Array arguments with intent(in)
REAL(KIND=real_jlslsm), INTENT(IN) :: runoff_riv(riv_pts)
!  runoff rate on rivers grid (kg m-2 s-1)

! Array arguments with intent(out)
REAL(KIND=real_jlslsm), INTENT(OUT) :: runoff_land(land_pts)
!  runoff rate on model grid (kg m-2 s-1)

INTEGER :: maxlin, l, ip, ix, iy

REAL(KIND=real_jlslsm) :: runoff_grid(nx_grid,ny_grid)
REAL(KIND=real_jlslsm) :: runoff_out(nx_rivers,ny_rivers)

! Convert from riv_pts to rivers grid
runoff_out(:,:) = 0.0
DO ip = 1,riv_pts
  CALL rivers_get_xy_pos(rivers_index_rp(ip),nx_rivers,ny_rivers,ix,iy)
  runoff_out(ix,iy) = runoff_riv(ip)
END DO

! Call to lat/lon grid-based regridding routine
maxlin = ( nx_grid + nx_rivers ) * ( ny_grid + ny_rivers )
CALL rivers_route_regrid_invert ( maxlin, runoff_out, runoff_grid )

! Convert from full grid to land_pts
DO l = 1,land_pts
  CALL rivers_get_xy_pos( ir_land_grid(l), nx_grid, ny_grid, ix, iy )
  runoff_land(l) = runoff_grid(ix,iy)
END DO

END SUBROUTINE rivers_regrid_to_land

!###############################################################################
! subroutine rivers_route_regrid
! Handles regridding of runoff from "main" to rivers grids.

SUBROUTINE rivers_route_regrid( maxlin, runoff_grid, runoff_out )
!------------------------------------------------------------------------------
! Description:
!   Regrids runoff from a source grid to a target grid (the rivers grid).
!   Both grids must be regular in latitude and longitude.
!
!------------------------------------------------------------------------------
! Modules used:

USE jules_rivers_mod, ONLY:                                                   &
!  imported scalars with intent(in)
     nx_rivers,ny_rivers,rivers_dlat,rivers_dlon,rivers_lat1,rivers_lon1      &
     , nx_grid, ny_grid, reg_dlat, reg_dlon, reg_lat1, reg_lon1

USE areaver_mod, ONLY: pre_areaver, do_areaver

IMPLICIT NONE

! Scalar arguments with intent(in)

INTEGER, INTENT(IN) :: maxlin  !  a vector length

! Array arguments with intent(in)

REAL(KIND=real_jlslsm), INTENT(IN) :: runoff_grid(nx_grid,ny_grid)
                               !  runoff rate on model grid (kg m-2 s-1)

! Array arguments with intent(out)

REAL(KIND=real_jlslsm), INTENT(OUT) :: runoff_out(nx_rivers,ny_rivers)
                               !  runoff rate on rivers grid (kg m-2 s-1)

! Local scalar variables

INTEGER :: adjust       !  mode (option number) of area averaging
INTEGER :: ix           !  loop counter
INTEGER :: iy           !  loop counter
INTEGER :: maxl         !  number of entries in output from pre_areaver

LOGICAL :: cyclic_srce  !  TRUE source (model) grid is cyclic in x
LOGICAL :: cyclic_targ  !  TRUE target (model) grid is cyclic in x
LOGICAL :: invert_srce  !  TRUE source (model) grid runs N to S
                        !  FALSE source grid runs S to N
LOGICAL :: spherical    !  TRUE  coordinates are lat/lon on a sphere
                        !  FALSE Cartesian axes.
LOGICAL :: want         !  Value of masks at locations for which data are
                        !  required

! Local array variables.

INTEGER :: count_targ(nx_rivers,ny_rivers)
                        !  number of model gridboxes that
                        !  contribute to each rivers gridbox
INTEGER :: base_targ(nx_rivers,ny_rivers)
                        !  base (starting) index for each rivers gridbox
                        !  (i.e. location of first element in lists)
INTEGER :: index_srce(maxlin)
                        !  list of source (model) gridboxes
                        !  that contribute to each rivers gridbox

REAL(KIND=real_jlslsm) :: sourcelat(ny_grid+1)
                        !  latitudes of edges of source (model) gridboxes.
                        !  First value is the S edge of first gridbox,
                        !  all other values are the N edge of each gridbox.
REAL(KIND=real_jlslsm) :: sourcelon(nx_grid+1)
                        !  longitudes of edges of source (model) gridboxes
                        !  First value is the W edge of first gridbox, all
                        !  other values are the E edge of each gridbox.
REAL(KIND=real_jlslsm) :: targetlat(ny_rivers+1)
                        !  latitudes of edges of target (mrivers) gridboxes.
                        !  First value is the S edge of first gridbox, all
                        !  other values are the N edge of each gridbox.
REAL(KIND=real_jlslsm) :: targetlon(nx_rivers+1)
                        !  longitudes of edges of target (rivers) gridboxes
                        !  First value is the W edge of first gridbox, all
                        !  other values are the E edge of each gridbox.

REAL(KIND=real_jlslsm) :: adjust_targ(nx_rivers,ny_rivers)
                        !  adjustment factors (not used)
REAL(KIND=real_jlslsm) :: weight(maxlin)
                        !  lists of weights for each source (model) gridbox

LOGICAL :: mask_srce(nx_grid,ny_grid)
LOGICAL :: mask_targ(nx_rivers,ny_rivers)

!------------------------------------------------------------------------------

! Set values

spherical = .TRUE.    !  Calculations are for lat/lon coordinates.
adjust    = 0         !  "normal" area averaging
maxl      = maxlin

! Decide if grids run S-N or N-S in latitude

IF ( reg_dlat > 0.0) THEN
  invert_srce = .FALSE.                   !  model grid runs S to N
ELSE
  invert_srce = .TRUE.                    !  model grid runs N to S
END IF

! Decide if grids are cyclic in longitude.

IF ( REAL(nx_grid) * reg_dlon > 359.9 ) THEN
  cyclic_srce = .TRUE.
ELSE
  cyclic_srce = .FALSE.
END IF
IF ( REAL(nx_rivers) * rivers_dlon > 359.9 ) THEN
  cyclic_targ = .TRUE.
ELSE
  cyclic_targ = .FALSE.
END IF

!------------------------------------------------------------------------------
! Set coordinates of edges of model gridboxes
!------------------------------------------------------------------------------

DO ix = 1,nx_grid+1
  sourcelon(ix) = reg_lon1 + (REAL(ix-1) - 0.5) * reg_dlon
END DO
DO iy = 1,ny_grid+1
  sourcelat(iy) = reg_lat1 + (REAL(iy-1) - 0.5) * reg_dlat
END DO

!------------------------------------------------------------------------------
! Set coordinates of edges of rivers gridboxes
!------------------------------------------------------------------------------

DO ix = 1,nx_rivers+1
  targetlon(ix) = rivers_lon1 + (REAL(ix-1) - 0.5) * rivers_dlon
END DO
DO iy = 1,ny_rivers+1
  targetlat(iy) = rivers_lat1 + (REAL(iy-1) - 0.5) * rivers_dlat
END DO

!------------------------------------------------------------------------------
! Set masks to indicate that all points in both grids are to be used
!------------------------------------------------------------------------------

want           = .TRUE.
mask_srce(:,:) = want
mask_targ(:,:) = want

!------------------------------------------------------------------------------
! Call setup routing for averaging
!------------------------------------------------------------------------------

CALL pre_areaver( nx_grid, sourcelon, ny_grid, sourcelat                      &
                 ,cyclic_srce, nx_grid, want, mask_srce                       &
                 ,nx_rivers, targetlon, ny_rivers, targetlat                  &
                 ,cyclic_targ, spherical, maxl, count_targ, base_targ         &
                 ,index_srce, weight )

!------------------------------------------------------------------------------
! Call averaging routine
!------------------------------------------------------------------------------

CALL do_areaver( nx_grid, ny_grid, nx_grid                                    &
                ,invert_srce, runoff_grid, nx_rivers, ny_rivers               &
                ,count_targ, base_targ, nx_rivers, want, mask_targ            &
                ,index_srce, weight, adjust, runoff_out, nx_grid, 2           &
                ,adjust_targ )

END SUBROUTINE rivers_route_regrid
!###############################################################################

!###############################################################################
! subroutine rivers_route_regrid_invert
! Handles regridding of runoff from "main" to rivers grids.

SUBROUTINE rivers_route_regrid_invert( maxlin, runoff_out, runoff_grid )
!------------------------------------------------------------------------------
! Description:
!   Regrids runoff from a source grid to a target grid (the rivers grid).
!   Both grids must be regular in latitude and longitude.
!   Inverse of rivers_route_regrid
!
!------------------------------------------------------------------------------
! Modules used:

USE jules_rivers_mod, ONLY:                                                   &
!  imported scalars with intent(in)
     nx_rivers,ny_rivers,rivers_dlat,rivers_dlon,rivers_lat1,rivers_lon1      &
     , nx_grid, ny_grid, reg_dlat, reg_dlon, reg_lat1, reg_lon1

USE areaver_mod, ONLY: pre_areaver, do_areaver

IMPLICIT NONE

! Scalar arguments with intent(in)

INTEGER, INTENT(IN) :: maxlin    !  a vector length

! Array arguments with intent(in)

REAL(KIND=real_jlslsm), INTENT(IN) :: runoff_out(nx_rivers,ny_rivers)
                                 !  runoff rate on rivers grid (kg m-2 s-1)

! Array arguments with intent(out)

REAL(KIND=real_jlslsm), INTENT(OUT) :: runoff_grid(nx_grid,ny_grid)
                                 !  runoff rate on model grid (kg m-2 s-1)

! Local scalar variables

INTEGER :: adjust        !  mode (option number) of area averaging
INTEGER :: ix            !  loop counter
INTEGER :: iy            !  loop counter
INTEGER :: maxl          !  number of entries in output from pre_areaver

LOGICAL :: cyclic_srce   !  TRUE source (model) grid is cyclic in x
LOGICAL :: cyclic_targ   !  TRUE target (model) grid is cyclic in x
LOGICAL :: invert_srce   !  TRUE source (model) grid runs N to S
                         !  FALSE source grid runs S to N
LOGICAL :: spherical     !  TRUE  coordinates are lat/lon on a sphere
                         !  FALSE Cartesian axes.
LOGICAL :: want
                         !  Value of masks at locations for which data
                         !  are required

! Local array variables

INTEGER :: count_targ(nx_grid,ny_grid)
                         !  number of model gridboxes that contribute to
                         !  each rivers gridbox
INTEGER :: base_targ(nx_grid,ny_grid)
                         !  base (starting) index for each rivers gridbox
                         !  (i.e. location of first element in lists)
INTEGER :: index_srce(maxlin)
                         !  list of source (model) gridboxes that contribute
                         !  to each rivers gridbox

REAL(KIND=real_jlslsm) :: sourcelat(ny_rivers+1)
                         !  latitudes of edges of target (model) gridboxes.
                         !  First value is the S edge of first gridbox, all
                         !  other values are the N edge of each gridbox.
REAL(KIND=real_jlslsm) :: sourcelon(nx_rivers+1)
                         !  longitudes of edges of target (model) gridboxes
                         !  First value is the W edge of first gridbox, all
                         !  other values are the E edge of each gridbox.
REAL(KIND=real_jlslsm) :: targetlat(ny_grid+1)
                         !  latitudes of edges of source (rivers) gridboxes.
                         !  First value is the S edge of first gridbox, all
                         !  other values are the N edge of each gridbox.
REAL(KIND=real_jlslsm) :: targetlon(nx_grid+1)
                         !  longitudes of edges of source (rivers) gridboxes
                         !  First value is the W edge of first gridbox, all
                         !  other values are the E edge of each gridbox.

REAL(KIND=real_jlslsm) :: adjust_targ(nx_grid,ny_grid)
                         !  adjustment factors (not used)
REAL(KIND=real_jlslsm) :: weight(maxlin)
                         !  lists of weights for each source (model) gridbox

LOGICAL :: mask_srce(nx_rivers,ny_rivers)
LOGICAL :: mask_targ(nx_grid,ny_grid)

!------------------------------------------------------------------------------

! Set values
spherical   = .TRUE.    !  Calculations are for lat/lon coordinates.
adjust      = 0         !  "normal" area averaging
maxl        = maxlin
invert_srce = .FALSE.   !  rivers grid runs S to N

! Decide if grids are cyclic in longitude

IF ( REAL(nx_rivers) * rivers_dlon > 359.9 ) THEN
  cyclic_srce = .TRUE.
ELSE
  cyclic_srce = .FALSE.
END IF
IF ( REAL(nx_grid) * reg_dlon > 359.9 ) THEN
  cyclic_targ = .TRUE.
ELSE
  cyclic_targ = .FALSE.
END IF

!------------------------------------------------------------------------------
! Set coordinates of edges of rivers gridboxes
!------------------------------------------------------------------------------

DO ix = 1,nx_rivers+1
  sourcelon(ix) = rivers_lon1 + (REAL(ix-1) - 0.5) * rivers_dlon
END DO
DO iy = 1,ny_rivers+1
  sourcelat(iy) = rivers_lat1 + (REAL(iy-1) - 0.5) * rivers_dlat
END DO

!------------------------------------------------------------------------------
! Set coordinates of edges of model gridboxes
!------------------------------------------------------------------------------

DO ix = 1,nx_grid+1
  targetlon(ix) = reg_lon1 + (REAL(ix-1) - 0.5) * reg_dlon
END DO
DO iy = 1,ny_grid+1
  targetlat(iy) = reg_lat1 + (REAL(iy-1) - 0.5) * reg_dlat
END DO

!------------------------------------------------------------------------------
! Set masks to indicate that all points in both grids are to be used
!------------------------------------------------------------------------------

want           = .TRUE.
mask_srce(:,:) = want
mask_targ(:,:) = want

!------------------------------------------------------------------------------
! Call setup rivers for averaging.
!------------------------------------------------------------------------------

CALL pre_areaver( nx_rivers, sourcelon, ny_rivers, sourcelat, cyclic_srce     &
                 ,nx_rivers, want, mask_srce, nx_grid, targetlon, ny_grid     &
                 ,targetlat, cyclic_targ, spherical, maxl, count_targ         &
                 ,base_targ, index_srce, weight )

!------------------------------------------------------------------------------
! Call averaging routine.
!------------------------------------------------------------------------------

CALL do_areaver( nx_rivers, ny_rivers, nx_rivers, invert_srce, runoff_out     &
                ,nx_grid, ny_grid, count_targ, base_targ, nx_grid, want       &
                ,mask_targ, index_srce, weight, adjust, runoff_grid           &
                ,nx_rivers, 2, adjust_targ )

END SUBROUTINE rivers_route_regrid_invert
!###############################################################################

END MODULE rivers_regrid_mod
