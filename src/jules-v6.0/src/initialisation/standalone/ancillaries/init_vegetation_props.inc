#if !defined(UM_JULES)
!******************************COPYRIGHT**************************************
! (c) Centre for Ecology and Hydrology.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237] 
!******************************COPYRIGHT**************************************
SUBROUTINE init_vegetation_props()

USE io_constants, ONLY: max_file_name_len, max_sdf_name_len, namelist_unit

USE string_utils_mod, ONLY: to_string

USE templating_mod, ONLY: tpl_has_var_name, tpl_substitute_var

USE model_interface_mod, ONLY: identifier_len, populate_var, get_var_id

USE jules_vegetation_mod, ONLY: photo_adapt, photo_acclim_model

USE dump_mod, ONLY: ancil_dump_read

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in ancils required for vegetation and checks them for consistency.
!-----------------------------------------------------------------------------

! Local parameters.
INTEGER, PARAMETER :: max_vars = 1
  ! The maximum possible number of variables that can be given.

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'INIT_VEGETATION_PROPS'

! Local scalar variables.
INTEGER ::                                                                    &
  error,                                                                      &
    ! Error indicator.
  i,                                                                          &
    ! Loop counter.
  nvars_file,                                                                 &
    ! The number of variables that will be set from the given file.
  nvars_required
    ! The number of variables that are required in this configuration.

CHARACTER(LEN=errormessagelength) :: iomessage
    ! I/O error message string.

! Local array variables.
CHARACTER(LEN=identifier_len) :: required_vars(max_vars)
  ! The variable identifiers of the required variables.

!-----------------------------------------------------------------------------
! Variables for and definition of the jules_vegetation_props namelist.
!-----------------------------------------------------------------------------
INTEGER :: nvars
  ! The number of variables listed.

LOGICAL :: read_from_dump
  ! T means read variables from the dump file (not an ancillary file).

CHARACTER(LEN=max_file_name_len) :: FILE
  ! The name of the file (or variable name template) to use for variables that
  ! need to be filled from file.

REAL :: const_val(max_vars)
  ! The constant value to use for a variable if use_file = F for that
  ! variable.

LOGICAL :: use_file(max_vars)
  ! T - the variable is in the file
  ! F - the variable is set using a constant value.

CHARACTER(LEN=identifier_len) :: var(max_vars)
  ! The variable identifiers of the variables.

CHARACTER(LEN=max_sdf_name_len) :: var_name(max_vars)
  ! The name of each variable in the file.

CHARACTER(LEN=max_sdf_name_len) :: tpl_name(max_vars)
  ! The name to substitute in a template for each variable.

NAMELIST  / jules_vegetation_props/ read_from_dump, FILE, nvars, var, use_file, &
                           var_name, tpl_name, const_val

!-----------------------------------------------------------------------------
! If thermal adaptation is not selected, we have nothing more to do.
!-----------------------------------------------------------------------------
IF ( photo_acclim_model /= photo_adapt ) RETURN

!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
read_from_dump = .FALSE.
nvars_required = 0
nvars_file     = 0
nvars          = 0
use_file(:)    = .TRUE.  ! Default is for every variable to be read from file
FILE           = ''      ! Empty file name
tpl_name(:)    = ''      ! Empty template string

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
CALL log_info(RoutineName, "Reading JULES_VEGETATION_PROPS namelist...")

READ(namelist_unit, NML = jules_vegetation_props, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 ) THEN
  CALL log_fatal(RoutineName,                                                 &
                 "Error reading namelist JULES_VEGETATION_PROPS " //          &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")
END IF

ancil_dump_read%vegetation_props = read_from_dump

IF ( .NOT. ancil_dump_read%vegetation_props) THEN
  !---------------------------------------------------------------------------
  ! We read from the ancil file.
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  ! Set up properties using namelist values.
  !---------------------------------------------------------------------------
  ! Set up the required variables.
  ! All the variables are required.
  nvars_required   = max_vars
  required_vars(:) = (/ 't_growth_gb' /)

  !---------------------------------------------------------------------------
  ! Check that all the required variables are provided.
  !---------------------------------------------------------------------------
  DO i = 1,nvars_required
    IF ( .NOT. ANY(var(1:nvars) == required_vars(i)) ) THEN
      CALL log_fatal(RoutineName,                                             &
                     "No value given for required variable '" //              &
                     TRIM(required_vars(i)) // "'")
    END IF
  END DO

  !---------------------------------------------------------------------------
  ! Check which variables we will be using and partition them into variables
  ! set to constant values and variables set from file.
  !---------------------------------------------------------------------------
  DO i = 1,nvars

    IF ( ANY(required_vars(1:nvars_required) == var(i)) ) THEN
      ! This is a required variable.
      IF ( use_file(i) ) THEN
        ! Variable will be filled from file; register it here.
        CALL log_info(RoutineName,                                            &
                      "'" // TRIM(var(i)) // "' will be read from file")

        nvars_file = nvars_file + 1
        ! Since nvars_file <= i (so we will not overwrite unprocessed values)
        ! and we do not need the values from these arrays for any non-file
        ! variables from now on, we can just compress them down onto variables
        ! that are in the file.
        var(nvars_file)      = var(i)
        var_name(nvars_file) = var_name(i)
        tpl_name(nvars_file) = tpl_name(i)
      ELSE
        ! Variable will not be read from file. Populate using a constant.
        CALL log_info(RoutineName,                                            &
                      "'" // TRIM(var(i)) // "' will be set to a " //         &
                      "constant = " // to_string(const_val(i)))

        CALL populate_var(get_var_id(var(i)), const_val = const_val(i))
      END IF

    ELSE

      ! This is not a required variable. Warn about not using it.
      CALL log_warn(RoutineName,                                              &
                    "Provided variable '" // TRIM(var(i)) //                  &
                    "' is not required, so will be ignored")
    END IF

  END DO

  !---------------------------------------------------------------------------
  ! Set variables from file
  !---------------------------------------------------------------------------
  IF ( nvars_file > 0 ) THEN
    ! Check that a file name was provided.
    IF ( LEN_TRIM(FILE) == 0 ) THEN
      CALL log_fatal(RoutineName, "No file name provided")
    END IF

    IF ( tpl_has_var_name(FILE) ) THEN
      ! We are using a file name template, so loop through the variables setting
      ! one from each file.
      DO i = 1,nvars_file
        ! If using a variable name template, check that a template string was
        ! provided for the current variable.
        IF ( LEN_TRIM(tpl_name(i)) == 0 ) THEN
          CALL log_fatal(RoutineName,                                         &
                         "No variable name template substitution " //         &
                         "(tpl_name) provided for " // TRIM(var(i)))
        END IF

        CALL fill_variables_from_file( tpl_substitute_var(FILE, tpl_name(i)), &
                                       (/ var(i) /), (/ var_name(i) /) )

      END DO

    ELSE

      ! We are not using a file name template, so set all variables from the
      ! same file.
      CALL fill_variables_from_file( FILE, var(1:nvars_file),                 &
                                     var_name(1:nvars_file) )

    END IF  !  tpl_has_var_name

  END IF  !  nvars_file

ELSE

  !---------------------------------------------------------------------------
  ! We read from the dump file
  !---------------------------------------------------------------------------
  CALL log_info(RoutineName,                                                  &
                "Veg property ancils will be read from the dump file.  " //   &
                "Namelist values ignored")

END IF !  ancil_dump_read%vegetation_props

RETURN

END SUBROUTINE init_vegetation_props
#endif
