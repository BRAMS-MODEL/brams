! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This module contains variables used for reading in triffid data
! and initialisations


! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.


MODULE trif_io

USE max_dimensions, ONLY:                                                     &
  npft_max
USE missing_data_mod, ONLY: imdi, rmdi
USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!---------------------------------------------------------------------
! Set up variables to use in IO (a fixed size version of each array
! in trif that we want to initialise).
!---------------------------------------------------------------------
INTEGER ::                                                                    &
  crop_io(npft_max) = imdi

REAL(KIND=real_jlslsm) ::                                                     &
  g_area_io(npft_max) = rmdi,                                                 &
  g_grow_io(npft_max) = rmdi,                                                 &
  g_root_io(npft_max) = rmdi,                                                 &
  g_wood_io(npft_max) = rmdi,                                                 &
  lai_max_io(npft_max) = rmdi,                                                &
  lai_min_io(npft_max) = rmdi,                                                &
  alloc_fast_io(npft_max) = rmdi,                                             &
  alloc_med_io(npft_max) = rmdi,                                              &
  alloc_slow_io(npft_max) = rmdi,                                             &
  dpm_rpm_ratio_io(npft_max) = rmdi,                                          &
  retran_l_io(npft_max) = rmdi,                                               &
  retran_r_io(npft_max) = rmdi

!---------------------------------------------------------------------
! Set up a namelist for reading and writing these arrays
!---------------------------------------------------------------------
NAMELIST  / jules_triffid/ crop_io,g_area_io,g_grow_io,g_root_io,             &
                         g_wood_io,lai_max_io,lai_min_io,                     &
                         alloc_fast_io,alloc_med_io,alloc_slow_io,            &
                         dpm_rpm_ratio_io,retran_l_io,retran_r_io

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TRIF_IO'

CONTAINS
SUBROUTINE print_nlist_jules_triffid()
USE jules_print_mgr, ONLY: jules_print
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer

CALL jules_print('trif_io',                                                   &
    'Contents of namelist jules_triffid')

WRITE(lineBuffer,*)' crop_io = ',crop_io
CALL jules_print('trif_io',lineBuffer)
WRITE(lineBuffer,*)' g_area_io = ',g_area_io
CALL jules_print('trif_io',lineBuffer)
WRITE(lineBuffer,*)' g_grow_io = ',g_grow_io
CALL jules_print('trif_io',lineBuffer)
WRITE(lineBuffer,*)' g_root_io = ',g_root_io
CALL jules_print('trif_io',lineBuffer)
WRITE(lineBuffer,*)' g_wood_io = ',g_wood_io
CALL jules_print('trif_io',lineBuffer)
WRITE(lineBuffer,*)' lai_max_io = ',lai_max_io
CALL jules_print('trif_io',lineBuffer)
WRITE(lineBuffer,*)' lai_min_io = ',lai_min_io
CALL jules_print('trif_io',lineBuffer)
WRITE(lineBuffer,*)' alloc_fast_io = ',alloc_fast_io
CALL jules_print('trif_io',lineBuffer)
WRITE(lineBuffer,*)' alloc_med_io = ',alloc_med_io
CALL jules_print('trif_io',lineBuffer)
WRITE(lineBuffer,*)' alloc_slow_io = ',alloc_slow_io
CALL jules_print('trif_io',lineBuffer)
WRITE(lineBuffer,*)' dpm_rpm_ratio_io = ',dpm_rpm_ratio_io
CALL jules_print('trif_io',lineBuffer)
WRITE(lineBuffer,*)' retran_r_io = ',retran_r_io
CALL jules_print('trif_io',lineBuffer)
WRITE(lineBuffer,*)' retran_l_io = ',retran_l_io
CALL jules_print('trif_io',lineBuffer)

CALL jules_print('trif_io',                                                   &
    '- - - - - - end of namelist - - - - - -')

END SUBROUTINE print_nlist_jules_triffid


END MODULE trif_io
