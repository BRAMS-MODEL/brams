! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE jules_irrig_mod

!-----------------------------------------------------------------------------
! Description:
!   Module containing irrigation options and a namelist for setting them
!-----------------------------------------------------------------------------


USE missing_data_mod, ONLY: imdi, rmdi
USE max_dimensions, ONLY: npft_max

IMPLICIT NONE

! Parameters identifying alternative ways of calculating the irrigation
! season. There are used with variable irr_crop and should have unique values.
INTEGER, PARAMETER ::                                                         &
  irr_crop_all_year = 0,                                                      &
    ! Irrigation season lasts all year.
  irr_crop_doell = 1,                                                         &
    ! Irrigation season is determined from driving data according to Doll and
    ! Siebert (2002) method.
    ! Doell, P., and Siebert, S., 2002, Global modeling of irrigation water
    ! requirements, Water Resour. Res., 38(4).
  irr_crop_dvimax = 2
    ! Irrigation season is determined by maximum crop development index
    ! across all tiles. Requires ncpft > 0. 

LOGICAL :: frac_irrig_all_tiles = .FALSE.
    ! T - assign irr fraction to all tiles
    ! F - assign irr fraction to specfied tiles only
INTEGER :: irrigtiles(npft_max) = imdi
    ! Indices of the tiles to be irrigated.
    ! Only used if frac_irrig_all_tiles = F.
    ! This is the fixed length equivalent of irrtiles where size is defined
    ! after namelist have been read in and set by irrig_vars_alloc

LOGICAL :: set_irrfrac_on_irrtiles = .FALSE.
    ! Switch to set irrigation fraction for irrigated tiles,
    ! set to FALSE to reproduce original results using frac_irrig_all_tiles

LOGICAL ::                                                                    &
  l_irrig_dmd = .FALSE.,                                                      &
    !   Switch for using irrigation demand code
  l_irrig_limit = .FALSE.
    !   Switch for limiting irrigation supply

INTEGER :: nirrtile = imdi
   !  Number of tiles that can have irrigated fraction

INTEGER :: irr_crop = irr_crop_all_year
   !  Switch for irrigation cropping model.

! Non-namelist variable but arrays used elsewhere after namelists are
! read in - moved out of crop_vars_mod (Hopefully this is the best place
! to put these).
INTEGER, ALLOCATABLE :: irrtiles(:)
      ! Tiles that can have irrigated fraction
      ! Only used when frac_irrig_all_tiles = .FALSE.
!-----------------------------------------------------------------------------
! Single namelist definition for UM and standalone
!-----------------------------------------------------------------------------
NAMELIST  / jules_irrig/  l_irrig_dmd, l_irrig_limit, irr_crop,               &
                          frac_irrig_all_tiles, nirrtile, irrigtiles,         &
                          set_irrfrac_on_irrtiles

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='JULES_IRRIG_MOD'

CONTAINS

SUBROUTINE irrig_vars_alloc(npft, l_irrig_dmd)

!No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
INTEGER, INTENT(IN) :: npft
LOGICAL, INTENT(IN) :: l_irrig_dmd

!Local variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='IRRIG_VARS_ALLOC'
!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( l_irrig_dmd ) THEN
  ALLOCATE(irrtiles(npft))
  irrtiles(:) = 0
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE irrig_vars_alloc


SUBROUTINE check_jules_irrig()

USE ereport_mod, ONLY: ereport
USE jules_print_mgr, ONLY: jules_message, jules_print, PrNorm
USE jules_hydrology_mod, ONLY: l_top
USE jules_vegetation_mod, ONLY: l_crop
USE jules_surface_types_mod, ONLY: c3_grass, c4_grass
!-----------------------------------------------------------------------------
! Description:
!   Checks JULES_IRRIG namelist for consistency
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IMPLICIT NONE

INTEGER :: errcode, error_sum   ! error code to pass to ereport
CHARACTER(LEN=*), PARAMETER :: RoutineName='CHECK_JULES_IRRIG'


! Check that irr_crop is reasonable.
IF ( l_irrig_dmd ) THEN
  SELECT CASE ( irr_crop )
  CASE ( irr_crop_all_year, irr_crop_doell, irr_crop_dvimax )
    ! These are valid, nothing more to do.
  CASE DEFAULT
    errcode = 101
    CALL ereport(RoutineName, errcode,                                        &
               'Invalid value for irr_crop.')
  END SELECT
END IF

! Verify that if switches are on that other switches are also available
IF ( l_irrig_limit .AND. .NOT. l_irrig_dmd ) THEN
  errcode = 101
  CALL ereport(RoutineName, errcode,                                          &
               'l_irrig_limit=T requires l_irrig_dmd=T ')
END IF

IF ( l_irrig_limit .AND. .NOT. l_top ) THEN
  errcode = 101
  CALL ereport(RoutineName, errcode,                                          &
               'l_irrig_limit=T requires l_top=T ')
END IF

IF ( l_irrig_dmd .AND. irr_crop ==  irr_crop_dvimax .AND. .NOT. l_crop ) THEN
  errcode = 101
  CALL ereport(RoutineName, errcode,                                          &
               'Irrigation triggered by crop DVI (irr_crop = 2) ' //          &
               'requires crop model to be active')
END IF

IF ( set_irrfrac_on_irrtiles .AND. frac_irrig_all_tiles ) THEN
  errcode = 101
  CALL ereport(RoutineName, errcode,                                          &
               "Cannot set both frac_irrig_all_tiles and  " //                &
               "set_irrfrac_on_irrtiles to be true ")
END IF

END SUBROUTINE check_jules_irrig


SUBROUTINE print_nlist_jules_irrig()

USE jules_print_mgr, ONLY: jules_print

IMPLICIT NONE

CHARACTER(LEN=50000) :: lineBuffer

!-----------------------------------------------------------------------------

CALL jules_print('jules_irrig', 'Contents of namelist jules_irrig')

WRITE(lineBuffer,*) ' l_irrig_dmd = ',l_irrig_dmd
CALL jules_print('jules_irrig',lineBuffer)

WRITE(lineBuffer,*) ' l_irrig_limit = ',l_irrig_limit
CALL jules_print('jules_irrig',lineBuffer)

WRITE(lineBuffer,*) ' irr_crop = ',irr_crop
CALL jules_print('jules_irrig',lineBuffer)

WRITE(lineBuffer,*) ' frac_irrig_all_tiles = ',frac_irrig_all_tiles
CALL jules_print('jules_irrig',lineBuffer)

WRITE(lineBuffer,*) ' nirrtile = ',nirrtile
CALL jules_print('jules_irrig',lineBuffer)

WRITE(lineBuffer,*) ' irrigtiles = ',irrigtiles
CALL jules_print('jules_irrig',lineBuffer)

WRITE(lineBuffer,*) ' set_irrfrac_on_irrtiles = ',set_irrfrac_on_irrtiles
CALL jules_print('jules_irrig',lineBuffer)

CALL jules_print('jules_irrig',                                               &
    '- - - - - - end of namelist - - - - - -')


END SUBROUTINE print_nlist_jules_irrig


SUBROUTINE read_nml_jules_irrig(nml_dir)

! Description:
!  Read the JULES_IRRIG namelist (standalone)

USE io_constants, ONLY: max_file_name_len, namelist_unit
USE string_utils_mod, ONLY: to_string
USE logging_mod, ONLY: log_info, log_fatal
USE errormessagelength_mod, ONLY: errormessagelength
USE missing_data_mod, ONLY: imdi, rmdi


IMPLICIT NONE

! Arguments
CHARACTER(LEN=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                         ! namelists
INTEGER :: error  ! Error indicator
CHARACTER(LEN=errormessagelength) :: iomessage

CHARACTER(LEN=max_file_name_len) :: FILE
                      ! The name of the file (or variable name template) to
                      ! use for variables that need to be filled from file

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_JULES_IRRIG'

! Open the irrigation namelist file
OPEN(namelist_unit, FILE=(TRIM(nml_dir) // '/' // 'jules_irrig.nml'),         &
               STATUS='old', POSITION='rewind', ACTION='read', IOSTAT = error,&
               IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal(RoutineName,                                                 &
                 "Error opening namelist file jules_irrig.nml " //            &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")
! Read namelist
CALL log_info(RoutineName, "Reading JULES_IRRIG namelist...")
READ(namelist_unit, NML = jules_irrig, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal(RoutineName,                                                 &
                 "Error reading namelist JULES_IRRIG " //                     &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")


! Close the namelist file
CLOSE(namelist_unit, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal(RoutineName,                                                 &
                 "Error closing namelist file irrig.nml " //                  &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")



END SUBROUTINE read_nml_jules_irrig




END MODULE jules_irrig_mod
