#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_triffid(nml_dir)

USE missing_data_mod, ONLY:                                                   &
!  imported scalar parameters
     rmdi

USE io_constants, ONLY: namelist_unit

USE string_utils_mod, ONLY: to_string

USE jules_surface_types_mod, ONLY: npft, nnpft

USE jules_vegetation_mod, ONLY: l_triffid, l_phenol

USE trif, ONLY: crop, g_area, g_grow, g_root, g_wood, lai_max, lai_min,       &
                alloc_fast, alloc_med, alloc_slow, dpm_rpm_ratio,             &
                retran_r, retran_l

USE trif_io, ONLY: jules_triffid,                                             &
                         crop_io,g_area_io,g_grow_io,g_root_io,               &
                         g_wood_io,lai_max_io,lai_min_io,                     &
                         alloc_fast_io,alloc_med_io,alloc_slow_io,            &
                         dpm_rpm_ratio_io,retran_l_io,retran_r_io
  
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises the TRIFFID parameters
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Arguments
CHARACTER(LEN=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                         ! namelists

! Work variables
INTEGER :: error  ! Error indicator
CHARACTER(LEN=errormessagelength) :: iomessage


!-----------------------------------------------------------------------------


! Nothing to do if neither triffid or phenology are selected
IF ( .NOT. l_triffid .AND. .NOT. l_phenol ) RETURN

!-----------------------------------------------------------------------------
! Initialise namelist values before reading them
!-----------------------------------------------------------------------------
crop_io(:) = -9

g_area_io(:)  = rmdi
g_grow_io(:)  = rmdi
g_root_io(:)  = rmdi
g_wood_io(:)  = rmdi
lai_max_io(:) = rmdi
lai_min_io(:) = rmdi
alloc_fast_io(:) = rmdi
alloc_med_io(:)  = rmdi
alloc_slow_io(:) = rmdi
dpm_rpm_ratio_io(:) = rmdi
retran_r_io(:)   = rmdi
retran_l_io(:)   = rmdi

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
CALL log_info("init_triffid", "Reading JULES_TRIFFID namelist...")

! Open the pft parameters namelist file
OPEN(namelist_unit, FILE=(TRIM(nml_dir) // '/' // 'triffid_params.nml'),      &
               STATUS='old', POSITION='rewind', ACTION='read', IOSTAT = error,&
               IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_triffid",                                              &
                 "Error opening namelist file triffid_params.nml " //         &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

READ(namelist_unit, NML = jules_triffid, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_triffid",                                              &
                 "Error reading namelist JULES_TRIFFID " //                   &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

! Close the namelist file
CLOSE(namelist_unit, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_triffid",                                              &
                 "Error closing namelist file triffid_params.nml " //         &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")


!-----------------------------------------------------------------------------
! Copy values from fixed length arrays used in namelist into allocated arrays
!-----------------------------------------------------------------------------
crop(:)    = crop_io(1:npft)
g_area(:)  = g_area_io(1:npft)
g_grow(:)  = g_grow_io(1:npft)
g_root(:)  = g_root_io(1:npft)
g_wood(:)  = g_wood_io(1:npft)
lai_max(:) = lai_max_io(1:npft)
lai_min(:) = lai_min_io(1:npft)
alloc_fast(:) = alloc_fast_io(1:npft)
alloc_med(:)  = alloc_med_io(1:npft)
alloc_slow(:) = alloc_slow_io(1:npft)
dpm_rpm_ratio(:) = dpm_rpm_ratio_io(1:npft)
retran_r(:)   = retran_r_io(1:npft)
retran_l(:)   = retran_l_io(1:npft)
 
!-----------------------------------------------------------------------------
! Check that all required variables were present in the namelist.
! The namelist variables were initialised to rmdi.
!-----------------------------------------------------------------------------
error = 0
IF ( ANY( crop(1:nnpft) < 0 ) ) THEN  ! crop was initialised to < 0
  error = 1
  CALL log_error("init_triffid", "No value for crop")
END IF
IF ( ANY( ABS( g_area(1:nnpft) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error("init_triffid", "No value for g_area")
END IF
IF ( ANY( ABS( g_grow(1:nnpft) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error("init_triffid", "No value for g_grow")
END IF
IF ( ANY( ABS( g_root(1:nnpft) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error("init_triffid", "No value for g_root")
END IF
IF ( ANY( ABS( g_wood(1:nnpft) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error("init_triffid", "No value for g_wood")
END IF
IF ( ANY( ABS( lai_max(1:nnpft) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error("init_triffid", "No value for lai_max")
END IF
IF ( ANY( ABS( lai_min(1:nnpft) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error("init_triffid", "No value for lai_min")
END IF
IF ( ANY( ABS( alloc_fast(1:nnpft) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error("init_triffid", "No value for alloc_fast")
END IF
IF ( ANY( ABS( alloc_med(1:nnpft) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error("init_triffid", "No value for alloc_med")
END IF
IF ( ANY( ABS( alloc_slow(1:nnpft) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error("init_triffid", "No value for alloc_slow")
END IF
IF ( ANY( ABS( dpm_rpm_ratio(1:nnpft) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error("init_triffid", "No value for dpm_rpm_ratio")
END IF
IF ( ANY( ABS( retran_r(1:nnpft) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error("init_triffid", "No value for retran_r")
END IF
IF ( ANY( ABS( retran_l(1:nnpft) - rmdi ) < EPSILON(1.0) ) ) THEN
  error = 1
  CALL log_error("init_triffid", "No value for retran_l")
END IF

IF ( error /= 0 )                                                             &
  CALL log_fatal("init_triffid",                                              &
                 "Variable(s) missing from namelist - see earlier " //        &
                 "error message(s)")

RETURN

END SUBROUTINE init_triffid
#endif
