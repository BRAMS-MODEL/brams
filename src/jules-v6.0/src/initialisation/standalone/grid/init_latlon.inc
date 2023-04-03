#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

SUBROUTINE init_latlon()

USE mem_brams_jules, ONLY: glatB,glonB !JWA

USE dump_mod, ONLY: ancil_dump_read

USE errormessagelength_mod, ONLY: errormessagelength

USE input_mod, ONLY: input_grid => grid, fill_variables_from_file

USE io_constants, ONLY: max_sdf_name_len, max_file_name_len, namelist_unit

USE jules_water_resources_mod, ONLY: l_water_irrigation

USE missing_data_mod, ONLY: rmdi

USE model_grid_mod, ONLY: grid_area_ij, grid_lat => latitude,                 &
      grid_lon => longitude

USE model_interface_mod, ONLY: identifier_len, populate_var, get_var_id

USE string_utils_mod, ONLY: to_string

USE templating_mod, ONLY: tpl_has_var_name, tpl_substitute_var

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises the latitude and longitude of each gridbox.
!   Optionally initialises gridbox area.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Local parameters.
!-----------------------------------------------------------------------------

INTEGER, PARAMETER :: max_vars = 3
  ! The maximum possible number of variables that can be given.

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_LATLON'

!-----------------------------------------------------------------------------
! Local scalar variables:
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  error,                                                                      &
    ! Error code.
  error_sum,                                                                  &
    ! Error indicator.
  i,                                                                          &
    ! Loop counter.
  nvars_file,                                                                 &
    ! The number of variables that will be set from the given file.
  nvars_required
    ! The number of variables that are required in this configuration.

CHARACTER(LEN=errormessagelength) :: iomessage

!-----------------------------------------------------------------------------
! Local array variables:
!-----------------------------------------------------------------------------
CHARACTER(LEN=max_sdf_name_len) :: file_tpl_name(max_vars)
     ! The name to substitute in a template for each variable.
CHARACTER(LEN=identifier_len) :: file_var(max_vars)
     ! The variable identifiers of the variables to set from file.
CHARACTER(LEN=max_sdf_name_len) :: file_var_name(max_vars)
     ! The name of each variable in the file
CHARACTER(LEN=identifier_len) :: required_vars(max_vars)
     ! The variable identifiers of the required variables.

!-----------------------------------------------------------------------------
! Definition of the jules_latlon namelist
!-----------------------------------------------------------------------------
INTEGER :: nvars
  ! The number of variables to be set.

LOGICAL :: read_from_dump
  ! Switch to read values from dump file.

CHARACTER(LEN=max_file_name_len) :: FILE
  ! The name of the file (or variable name template) to use for variables
  ! that need to be filled from file.

REAL(KIND=real_jlslsm) :: const_val(max_vars)
  ! The constant value to use for each variable if use_file = F for that
  ! variable.

LOGICAL :: use_file(max_vars)
  ! T - the variable uses the file
  ! F - the variable is set using a constant value

CHARACTER(LEN=identifier_len) :: var(max_vars)
  ! The identifiers of the variables.

CHARACTER(LEN=max_sdf_name_len) :: var_name(max_vars)
  ! The name of each variable in the file.

CHARACTER(LEN=max_sdf_name_len) :: tpl_name(max_vars)
  ! The name to substitute in a template for each variable.

NAMELIST  / jules_latlon/ read_from_dump, FILE, nvars, var,                   &
                          use_file, var_name, tpl_name, const_val

!-----------------------------------------------------------------------------
!end of header

!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
read_from_dump = .FALSE.
nvars_required = 0
nvars_file     = 0
nvars          = 0
use_file(:)    = .TRUE.  ! Default is for every var to be read from file
FILE=''                  ! Empty file name.
var(:)         = ''      ! Empty identifiers.
var_name(:)    = ''      ! Empty variable names.
tpl_name(:)    = ''      ! Empty template name.
const_val(:)   = rmdi    ! Missing data value.

CALL log_info(RoutineName,                                                    &
              "Getting latitude and longitude for the full input grid...")

! Allocate the model lat and lon arrays to the full input grid (for now).
! For ease we always allocate grid_area_ij, though it is not always needed.
error_sum = 0
ALLOCATE(grid_area_ij(input_grid%nx,input_grid%ny), stat = error)
error_sum = error_sum + error
ALLOCATE(grid_lat(input_grid%nx,input_grid%ny), stat = error)
error_sum = error_sum + error
ALLOCATE(grid_lon(input_grid%nx,input_grid%ny), stat = error)
error_sum = error_sum + error
IF ( error_sum > 0 ) THEN
  CALL log_fatal(RoutineName,                                                 &
                 "Error allocating arrays for latitude and longitude")
END IF

!-----------------------------------------------------------------------------
! Read the namelist.
!-----------------------------------------------------------------------------
CALL log_info(RoutineName, "Reading JULES_LATLON namelist...")

READ(namelist_unit, NML = jules_latlon, IOSTAT = error, IOMSG = iomessage)

IF ( error /= 0 ) THEN
 CALL log_fatal(RoutineName,                                                 &
                "Error reading namelist JULES_LATLON " //                    &
                "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                TRIM(iomessage) // ")")
END IF

ancil_dump_read%latlon = read_from_dump

!-----------------------------------------------------------------------------
! Identify which variables will be read here.
!-----------------------------------------------------------------------------
! First, set up the variables that are always required.
! latitude and longitude are always required here, regardless of the value of
! read_from_dump.
nvars_required = 2
required_vars(1:nvars_required) = (/ 'latitude ',                             &
                                     'longitude' /)

! If irrigation is used, add gridbox area (unless it will be read from dump).
IF ( l_water_irrigation ) THEN
  IF ( .NOT. ancil_dump_read%latlon ) THEN
    nvars_required = nvars_required + 1
    required_vars(nvars_required) =    'grid_area'
  END IF
ELSE
  ! Clarify that no variables from this section need to be read from dump.
  ancil_dump_read%latlon = .FALSE.
END IF

!-----------------------------------------------------------------------------
! Set variables.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Check that all the required variables are provided.
!-----------------------------------------------------------------------------
DO i = 1,nvars_required
  IF ( .NOT. ANY(var(1:nvars) == required_vars(i)) ) THEN
    CALL log_fatal(RoutineName,                                               &
                   "No value given for required variable '" //                &
                   TRIM(required_vars(i)) // "'")
  END IF
END DO

!-----------------------------------------------------------------------------
! Check which variables we will be using and partition them into variables
! set to constant values and variables set from file.
!-----------------------------------------------------------------------------
DO i = 1,nvars

  !---------------------------------------------------------------------------
  ! If the variable is required, we will be using it.
  !---------------------------------------------------------------------------
  IF ( ANY( required_vars(1:nvars_required) == TRIM(var(i)) ) ) THEN

    IF ( use_file(i) ) THEN
      ! This variable will be read from file.
      CALL log_info(RoutineName,                                              &
                    "'" // TRIM(var(i)) // "' will be read from file")

      ! If the variable will be filled from file, register it here.
      nvars_file = nvars_file + 1
      file_var(nvars_file) = var(i)
      file_var_name(nvars_file) = var_name(i)
      file_tpl_name(nvars_file) = tpl_name(i)
    ELSE
      ! THis variables will be set as a constant; populate it here.
      ! First check that a value has been provided.
      IF ( ABS( const_val(i) - rmdi ) < EPSILON(1.0) ) THEN
        CALL log_fatal(RoutineName,                                           &
                       "No constant value provided for variable '"            &
                       // TRIM(var(i)) // "'" )
      END IF
      CALL log_info(RoutineName,                                              &
                    "'" // TRIM(var(i)) // "' will be set to a " //           &
                    "constant = " // to_string(const_val(i)))
      CALL populate_var(get_var_id(var(i)), const_val = const_val(i))
    END IF

  ELSE

    ! The variable is not a required variable. Warn about not using it.
    CALL log_warn(RoutineName,                                                &
                  "Provided variable '" // TRIM(var(i)) //                    &
                  "' is not required, so will be ignored")
  END IF

END DO  !  i (loop over variables)

!-----------------------------------------------------------------------------
! Set variables from file
!-----------------------------------------------------------------------------
!DSM IF ( nvars_file > 0 ) THEN

!DSM  ! Check that a file name was provided.
!DSM  IF ( LEN_TRIM(FILE) == 0 ) THEN
!DSM    CALL log_fatal(RoutineName, "No file name provided")
!DSM  END IF

!DSM  IF ( tpl_has_var_name(FILE) ) THEN
!DSM    ! We are using a file name template, so loop through the variables
!DSM    ! setting one from each file.
!DSM    DO i = 1,nvars_file
!DSM      ! Check that a template string was provided for this variable.
!DSM      IF ( LEN_TRIM(file_tpl_name(i)) == 0 ) THEN
!DSM        CALL log_fatal( RoutineName,                                          &
!DSM                        "No variable name template substitution " //          &
!DSM                        "provided for " // TRIM(file_var(i)) )
!DSM      END IF
!DSM      CALL fill_variables_from_file(                                          &
!DSM                               tpl_substitute_var(FILE, file_tpl_name(i)),    &
!DSM                               (/ file_var(i) /), (/ file_var_name(i) /) )
!DSM    END DO
!DSM  ELSE
!DSM    ! We are not using a file name template, so set all variables from the
!DSM    ! same file.
!DSM    CALL fill_variables_from_file( FILE, file_var(1:nvars_file),              &
!DSM                                   file_var_name(1:nvars_file) )
!DSM  END IF


!DSM END IF  !  nvars_file

!DSM !-----------------------------------------------------------------------------
!DSM ! Advise if some values will be read from the dump file.
!DSM !-----------------------------------------------------------------------------
!DSM IF ( ancil_dump_read%latlon ) THEN
!DSM  CALL log_info(RoutineName,                                                  &
!DSM                "grid_area will be read from the dump file. " //              &
!DSM                "Namelist values will be ignored.")

!DSM END IF

!print*, 'Data is on a grid - reading latitude and longitude from BRAMS'
grid_lat=glatB
grid_lon=glonB

!-----------------------------------------------------------------------------
! Check that the values of latitude and longitude seem sensible.
!-----------------------------------------------------------------------------
IF ( ANY(grid_lat < -90.0) .OR. ANY(grid_lat > 90.0) )                        &
  CALL log_fatal(RoutineName,                                                 &
                 "Latitude is out of range - allowed range is -90.0 to " //   &
                 "90.0, given range is " //                                   &
                 TRIM(to_string(MINVAL(grid_lat))) // " to " //               &
                 TRIM(to_string(MAXVAL(grid_lat))))

IF ( ANY(grid_lon < -180.0) .OR. ANY(grid_lon > 360.0) )                      &
  CALL log_fatal(RoutineName,                                                 &
                 "Longitude is out of range - allowed range is -180.0 " //    &
                 "to 360.0, given range is " //                               &
                 TRIM(to_string(MINVAL(grid_lon))) // " to " //               &
                 TRIM(to_string(MAXVAL(grid_lon))))

RETURN

END SUBROUTINE init_latlon
#endif
