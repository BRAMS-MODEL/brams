! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module defines absolute maximum values for array dimensions
! that are used in IO

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.


MODULE max_dimensions

IMPLICIT NONE

INTEGER, PARAMETER ::                                                         &
  nnpft_max         = 32,                                                     &
      ! Maximum number of natural pfts
  elev_tile_max = 25,                                                         &
      ! Maximum number of elevated ice tiles
  ncpft_max         = 12,                                                     &
      ! Maximum number of crop pfts
  npft_max          = nnpft_max + ncpft_max,                                  &
  nnvg_max          = 10 + (elev_tile_max * 2),                               &
  ntype_max         = npft_max + nnvg_max,                                    &
  nsurft_max        = ntype_max,                                              &
  sm_levels_max     = 20,                                                     &
  cs_layer_max     = 60,                                                      &
    ! Maximum allowed number of soil carbon layers.
  snow_layers_max   = 10,                                                     &
  ndep_species_max  = 50
    ! Maximum allowed number of trace gas and aerosol components for
    ! deposition.

END MODULE max_dimensions
