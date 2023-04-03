#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE initial_conditions_mod

USE logging_mod, ONLY: log_info, log_warn, log_error, log_fatal

IMPLICIT NONE

PRIVATE
PUBLIC init_ic

CONTAINS

! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.
#include "init_ic.inc"
#include "total_snow_init.inc"
#include "topmodel_init.inc"
#include "get_default_ic_values.inc"
#include "calc_fit_fsat.inc"

END MODULE initial_conditions_mod

#endif
