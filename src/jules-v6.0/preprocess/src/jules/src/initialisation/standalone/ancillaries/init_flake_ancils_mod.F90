MODULE init_flake_ancils_mod

! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT****************************************
IMPLICIT NONE

PRIVATE

PUBLIC :: init_flake_ancils

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='INIT_FLAKE_ANCILS_MOD'

CONTAINS

!===============================================================================
! Public subroutine
!===============================================================================

SUBROUTINE init_flake_ancils()

!Module imports

USE dump_mod, ONLY: ancil_dump_read

USE errormessagelength_mod, ONLY: errormessagelength

USE input_mod, ONLY: fill_variables_from_file

USE io_constants, ONLY: max_file_name_len, max_sdf_name_len, namelist_unit

USE jules_surface_mod, ONLY: l_flake_model

USE logging_mod, ONLY: log_info, log_warn, log_fatal

USE missing_data_mod, ONLY: rmdi

USE model_interface_mod, ONLY: identifier_len, populate_var, get_var_id

USE string_utils_mod, ONLY: to_string

USE templating_mod, ONLY: tpl_has_var_name, tpl_substitute_var

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Subroutine for importing FLake ancillaries. 
!
! Current Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
INTEGER, PARAMETER :: max_flakemodel_vars = 1
       ! The maximum number of FLake model variables that can be given
       ! (lake depth is the only parameter set from namelist at the moment
       ! other paremeters, e.g. fetch, may be added in the future). 

INTEGER :: nvars_required      ! The number of variables that are
                               ! required in this configuration

CHARACTER(LEN=identifier_len) :: required_vars(max_flakemodel_vars)
                               ! The variable identifiers of the required
                               ! variables

INTEGER :: nvars_file          ! The number of variables that will be set
                               ! from the given file (template?)

INTEGER :: error  ! Error indicator

INTEGER :: i      ! Loop counter

!-----------------------------------------------------------------------------
! Definition of the jules_flake namelist
!-----------------------------------------------------------------------------

LOGICAL :: read_from_dump
CHARACTER(LEN=max_file_name_len) :: FILE
                      ! The name of the file (or variable name template) to
                      ! use for variables that need to be filled from file
INTEGER :: nvars      ! The number of variables in this section
CHARACTER(LEN=identifier_len) :: var(max_flakemodel_vars)
                      ! The variable identifiers of the variables
CHARACTER(LEN=max_sdf_name_len) :: var_name(max_flakemodel_vars)
                      ! The name of the lake depth variable in the file
CHARACTER(LEN=max_sdf_name_len) :: tpl_name(max_flakemodel_vars)
                      ! The name to substitute in a template for each
                      ! variable
LOGICAL :: use_file(max_flakemodel_vars)
                      !   T - the variable uses the file
                      !   F - the variable is set using a constant value
CHARACTER(LEN=errormessagelength) :: iomessage
                      ! I/O error message string
REAL :: const_val(max_flakemodel_vars)
                      ! The constant value to use for lake depth if
                      ! use_file = F

NAMELIST  / jules_flake/ read_from_dump, FILE, nvars, var, var_name, use_file,&
                       tpl_name, const_val

!-----------------------------------------------------------------------------

! If FLake model is not on, we have nothing to do
IF ( .NOT. l_flake_model ) RETURN

!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------

read_from_dump  = .FALSE.
nvars_required  = 0
nvars_file      = 0
nvars           = 0
use_file(:)     = .FALSE. ! Default is for every variable to be read from file
FILE(:)         = ''      ! Empty file names
var_name(:)     = ''      ! Empty variable names. 
tpl_name(:)     = ''      ! Empty template string
const_val(:)    = 5.0     ! Default lake depth is 5.0m everywhere

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
CALL log_info("init_flake_ancils", "Reading JULES_FLAKE namelist...")

READ(namelist_unit, NML = jules_flake, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_flake_ancils",                                         &
                 "Error reading namelist JULES_FLAKE " //                     &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
 
                TRIM(iomessage) // ")")

ancil_dump_read%flake = read_from_dump

IF ( .NOT. ancil_dump_read%flake ) THEN 
  !---------------------------------------------------------------------------
  ! Set up FLake properties using namelist values
  !---------------------------------------------------------------------------
  ! Set up the required variables
  ! All the FLake model variables are always required
  nvars_required = max_flakemodel_vars
  required_vars(:) = (/ 'lake_depth' /)

  !-------------------------------------------------------------------------
  ! Check that variable identifiers are not empty.
  ! Although we might later decide that the identifier is not required, for
  ! clarity we check here whether the claimed amount of information was
  ! provided.
  !-------------------------------------------------------------------------
  DO i = 1,nvars
    IF ( LEN_TRIM(var(i)) == 0 )                                              &
      CALL log_fatal("init_flake_ancils",                                     &
                     "Insufficient values for var. " //                       &
                     "No name provided for var at position #" //              &
                     TRIM(to_string(i)) )
  END DO

  !---------------------------------------------------------------------------
  ! Check that all the required variables are there
  !---------------------------------------------------------------------------
  DO i = 1,nvars_required
    IF ( .NOT. ANY(var(1:nvars) == required_vars(i)) )                        &
      CALL log_fatal("init_flake_ancils",                                     &
                     "No value given for required variable '" //              &
                     TRIM(required_vars(i)) // "'")
  END DO

  !---------------------------------------------------------------------------
  ! Check which variables we will be using and partition them into variables
  ! set to constant values and variables set from file
  !---------------------------------------------------------------------------
  DO i = 1,nvars
    !---------------------------------------------------------------------------
    ! If the variable is one of the required vars, then we will be using it
    !---------------------------------------------------------------------------
    IF ( ANY(required_vars(1:nvars_required) == var(i)) ) THEN
      IF ( use_file(i) ) THEN
        CALL log_info("init_flake_ancils",                                    &
                      "'" // TRIM(var(i)) // "' will be read from file")

        ! If the variable will be filled from file, register it here
        nvars_file = nvars_file + 1
        ! Since nvars_file <= i (so we will not overwrite unprocessed values)
        ! and we do not need the values from these arrays for any non-file
        !variables from now on, we can just compress them down onto variables
        !that are in the file
        var(nvars_file) = var(i)
        var_name(nvars_file) = var_name(i)
        tpl_name(nvars_file) = tpl_name(i)
      ELSE
        ! If the variable is being set as a constant, populate it here.
        ! First check that a value has been provided.
        IF ( ABS( const_val(i) - rmdi ) < EPSILON(1.0) )                      &
          CALL log_fatal("init_flake_ancils",                                 &
                         "No constant value provided for variable '"          &
                         // TRIM(var(i)) // "'" )

        CALL log_info("init_flake_ancils",                                    &
                      "'" // TRIM(var(i)) // "' will be set to a " //         &
                      "constant = " // to_string(const_val(i)))

        CALL populate_var(get_var_id(var(i)), const_val = const_val(i))
      END IF !use_file(i)
    ELSE
      ! If the variable is not a required variable, warn about not using it
      CALL log_warn("init_flake_ancils",                                      &
                    "Provided variable '" // TRIM(var(i)) //                  &
                    "' is not required, so will be ignored")
    END IF
  END DO !loop over nvars

  !---------------------------------------------------------------------------
  ! Set variables from file
  !---------------------------------------------------------------------------
  IF ( nvars_file > 0 ) THEN
    ! Check that a file name was provided.
    IF ( LEN_TRIM(FILE) == 0 )                                                &
      CALL log_fatal("init_flake_ancils", "No file name provided")

    IF ( tpl_has_var_name(FILE) ) THEN
      ! We are using a file name template, so loop through the variables setting
      ! one from each file
      DO i = 1,nvars_file
        ! If using a variable name template, check that a template string was
        !provided for the current variable
        IF ( LEN_TRIM(tpl_name(i)) == 0 )                                     &
          CALL log_fatal("init_flake_ancils",                                 &
                         "No variable name template substitution " //         &
                         "(tpl_name) provided for " // TRIM(var(i)))

        CALL fill_variables_from_file(                                        &
          tpl_substitute_var(FILE, tpl_name(i)),                              &
          (/ var(i) /), (/ var_name(i) /)                                     &
        )
      END DO
    ELSE
      ! We are not using a file name template, so set all variables from the same
      ! file
      CALL fill_variables_from_file(                                          &
        FILE, var(1:nvars_file), var_name(1:nvars_file)                       &
      )
    END IF !tpl_has_var_name
  END IF !nvars_file

ELSE !We read from the dump file
  CALL log_info("init_flake_ancils",                                          &
                "flake ancils will be read from the dump file.  " //          &
                "Namelist values ignored")

END IF !.NOT. ancil_dump_read%flake


RETURN
END SUBROUTINE init_flake_ancils
END MODULE init_flake_ancils_mod
