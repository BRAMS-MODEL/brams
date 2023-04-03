! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module for elevation values.

MODULE c_elevate

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

! Declarations:

USE max_dimensions, ONLY: nsurft_max
USE missing_data_mod, ONLY: rmdi

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

REAL(KIND=real_jlslsm) :: z_land_io = rmdi
                                   ! Single point gridbox mean elevation
REAL(KIND=real_jlslsm) :: surf_hgt_band(nsurft_max) = rmdi
                                        ! Spatially invariant elevation bands

LOGICAL :: l_elev_absolute_height(nsurft_max)
                                   ! F - tile surf_hgt are
                                   ! anomalies from forcing data
                                   ! altitude (default JULES behaviour)
                                   ! T - tile surf_hgt are absolute
                                   ! elevations above sea level

REAL(KIND=real_jlslsm) ::                                                     &
 surf_hgt_io(nsurft_max)

! In example JULES control files, these are all 0, so initialise at that
DATA surf_hgt_io / nsurft_max * 0.0 /
DATA l_elev_absolute_height / nsurft_max * .FALSE. /

! Namelist used in the UM only
NAMELIST  / jules_elevate/ surf_hgt_io, l_elev_absolute_height

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='C_ELEVATE'

CONTAINS

SUBROUTINE print_nlist_jules_elevate()
USE jules_print_mgr, ONLY: jules_print
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer

CALL jules_print('c_elevate',                                                 &
    'Contents of namelist jules_elevate')

WRITE(lineBuffer,*)' surf_hgt_io = ',surf_hgt_io
CALL jules_print('c_elevate',lineBuffer)

WRITE(lineBuffer,*)' l_elev_absolute_height = ',l_elev_absolute_height
CALL jules_print('c_elevate',lineBuffer)

CALL jules_print('c_elevate',                                                 &
    '- - - - - - end of namelist - - - - - -')

END SUBROUTINE print_nlist_jules_elevate


END MODULE c_elevate
