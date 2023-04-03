#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE init_grid_mod

USE logging_mod, ONLY: log_info, log_warn, log_fatal

IMPLICIT NONE

PRIVATE
PUBLIC init_grid

CONTAINS

! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.
#include "init_grid.inc"
#include "init_input_grid.inc"
#include "init_latlon.inc"
#include "init_model_grid.inc"
#include "init_land_frac.inc"
#include "init_surf_hgt.inc"
#include "init_z_land.inc"


END MODULE init_grid_mod

#endif
