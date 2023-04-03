#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE dump_mod

USE io_constants, ONLY: format_len, format_ascii, format_ncdf,                &
                         max_sdf_name_len

USE logging_mod, ONLY: log_info, log_fatal

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------------

INTEGER, PARAMETER :: max_dim_dump = 14
INTEGER, PARAMETER :: max_var_dump = 84

CHARACTER(LEN=max_sdf_name_len), PARAMETER ::                                 &
  land_dim_name    = "land",                                                  &
  pft_dim_name     = "pft",                                                   &
  cpft_dim_name    = "cpft",                                                  &
  sc_pool_dim_name = "scpool",                                                &
  soil_n_pool_dim_name = "snpool",                                            &
  sc_layer_dim_name = "sclayer",                                              &
  ch4layer_dim_name  = "ch4layer",                                            &
  snow_dim_name    = "snow",                                                  &
  soil_dim_name    = "soil",                                                  &
  tile_dim_name    = "tile",                                                  &
  type_dim_name    = "type",                                                  &
  soilt_dim_name   = "soilt",                                                 &
  scalar_dim_name  = "scalar",                                                &
  nolevs_dim_name  = "olevs",                                                 &
  nfarray_dim_name = "nfarray",                                               &
  seed_dim_name    = "seed",                                                  &
  bedrock_dim_name = "bedrock",                                               &
  p_rivers_dim_name = "p_rivers"

#if defined(NCDF_DUMMY)
CHARACTER(LEN=format_len), PARAMETER :: dump_format = format_ascii
#else
CHARACTER(LEN=format_len), PARAMETER :: dump_format = format_ncdf
#endif

  !Create a defined type to store the flags for which ancils are read
  !from the dump file
TYPE ancil_flags
  LOGICAL :: latlon
  LOGICAL :: frac
  LOGICAL :: vegetation_props
  LOGICAL :: soil_props
  LOGICAL :: top
  LOGICAL :: agric
  LOGICAL :: crop_props
  LOGICAL :: irrig
  LOGICAL :: rivers_props
  LOGICAL :: water_resources_props
  LOGICAL :: co2
  LOGICAL :: flake
END TYPE

TYPE(ancil_flags) :: ancil_dump_read

!-----------------------------------------------------------------------------
! Visibility declarations
!-----------------------------------------------------------------------------
PRIVATE
PUBLIC  max_var_dump, required_vars_for_configuration,                        &
        read_dump, write_dump,                                                &
        ancil_dump_read

CONTAINS


! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.
#include "required_vars_for_configuration.inc"
#include "read_dump.inc"
#include "write_dump.inc"
#include "get_dim_info.inc"

END MODULE dump_mod
#endif
