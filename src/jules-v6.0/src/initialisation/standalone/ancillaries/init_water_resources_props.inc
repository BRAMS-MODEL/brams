#if !defined(UM_JULES)

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Reads in JULES_WATER_RESOURCES_PROPS namelist if necessary

SUBROUTINE init_water_resources_props

USE dump_mod, ONLY: ancil_dump_read

USE io_constants, ONLY: max_file_name_len, max_sdf_name_len, namelist_unit

USE jules_water_resources_mod, ONLY: l_water_resources, l_water_irrigation

USE missing_data_mod, ONLY: imdi, rmdi

USE model_interface_mod, ONLY: identifier_len, populate_var, get_var_id

USE string_utils_mod, ONLY: to_string

USE templating_mod, ONLY: tpl_has_var_name, tpl_substitute_var

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!    Reads ancillary fields related to  water resource modelling.
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
INTEGER, PARAMETER :: max_vars = 4
  ! The maximum possible number of variables that can be given.

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_WATER_RESOURCES_PROPS'

!-----------------------------------------------------------------------------
! Local scalar variables:
!-----------------------------------------------------------------------------

INTEGER :: nvars_required      ! The number of water resource variables that 
                               ! are required in this configuration
INTEGER :: nvars_file       ! The number of variables that will be set
                            ! from the given file 

INTEGER :: error  ! Error indicator

INTEGER :: i  ! Loop counter

!-----------------------------------------------------------------------------
! Local array variables:
!-----------------------------------------------------------------------------
CHARACTER(LEN=identifier_len) :: required_vars(max_vars)
                               ! The variable identifiers of the required
                               ! variables
CHARACTER(LEN=identifier_len) :: file_var(max_vars)
                      ! The variable identifiers of the variables to set
                      ! from file
CHARACTER(LEN=max_sdf_name_len) :: file_var_name(max_vars)
                      ! The name of each variable in the file

CHARACTER(LEN=max_sdf_name_len) :: file_tpl_name(max_vars)
                      ! The name to substitute in a template for each
                      ! variable

!-----------------------------------------------------------------------------
! Variables that are in the namelist.
!-----------------------------------------------------------------------------
INTEGER :: nvars
  ! The number of variables to be set.

LOGICAL :: read_from_dump

REAL(KIND=real_jlslsm) :: const_val(max_vars)
  ! The constant value to use for each variable if use_file = F for that
  ! variable.

LOGICAL :: use_file(max_vars)
  ! T - the variable uses the file
  ! F - the variable is set using a constant value

CHARACTER(LEN=max_file_name_len) :: FILE
  ! The name of the file (or variable name template) to use for variables
  ! that need to be filled from file.

CHARACTER(LEN=identifier_len) :: var(max_vars)
  ! The identifiers of the variables.
CHARACTER(LEN=max_sdf_name_len) :: var_name(max_vars)
  ! The name of each variable in the file.
CHARACTER(LEN=max_sdf_name_len) :: tpl_name(max_vars)
  ! The name to substitute in a template for each variable.

NAMELIST  / jules_water_resources_props/ read_from_dump, FILE, nvars, var,    &
                            use_file, var_name, tpl_name, const_val

!-----------------------------------------------------------------------------
!end of header

! Nothing to do if water resource model is not selected.
IF ( .NOT. l_water_resources ) RETURN

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

!-----------------------------------------------------------------------------
! Read the namelist.
!-----------------------------------------------------------------------------
CALL log_info(RoutineName, "Reading JULES_WATER_RESOURCES_PROPS namelist...")

READ(namelist_unit, NML = jules_water_resources_props, IOSTAT = error)

IF ( error /= 0 ) THEN
  CALL log_fatal(RoutineName,                                                 &
                "Error reading namelist JULES_WATER_RESOURCES_PROPS " //      &
                "(IOSTAT=" // TRIM(to_string(error)) // ")")
END IF

ancil_dump_read%water_resources_props = read_from_dump

IF ( .NOT. ancil_dump_read%water_resources_props) THEN 
  !---------------------------------------------------------------------------
  ! Read from the ancil file.
  !---------------------------------------------------------------------------
  
  ! First, set up the variables that are always required.
  nvars_required = 2
  required_vars(1:nvars_required) = (/ 'conveyance_loss',                     &
                                       'sfc_water_frac ' /)

  ! If irrigation is used, add irrigation efficiency and gridbox area.
  IF ( l_water_irrigation ) THEN
    nvars_required = nvars_required + 1
    required_vars(nvars_required) =    'irrig_eff      '
    !    nvars_required = nvars_required + 1
    !    required_vars(nvars_required) =    'land_area      '
  END IF

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

    !-------------------------------------------------------------------------
    ! If the variable is required, we will be using it.
    !-------------------------------------------------------------------------
    IF ( ANY( required_vars(1:nvars_required) == TRIM(var(i)) ) ) THEN

      ! Decide if variable will come from a file or use a constant.
      IF ( use_file(i) ) THEN
        CALL log_info(RoutineName,                                            &
                      "'" // TRIM(var(i)) // "' will be read from file")

        ! If the variable will be filled from file, register it here.
        nvars_file = nvars_file + 1
        file_var(nvars_file) = var(i)
        file_var_name(nvars_file) = var_name(i)
        file_tpl_name(nvars_file) = tpl_name(i)
      ELSE
        ! If the variable is being set as a constant, populate it here.
        ! First check that a value has been provided.
        IF ( ABS( const_val(i) - rmdi ) < EPSILON(1.0) ) THEN
          CALL log_fatal(RoutineName,                                         &
                         "No constant value provided for variable '"          &
                         // TRIM(var(i)) // "'" )
        END IF

        CALL log_info(RoutineName,                                            &
                      "'" // TRIM(var(i)) // "' will be set to a " //         &
                      "constant = " // to_string(const_val(i)))
        CALL populate_var(get_var_id(var(i)), const_val = const_val(i))
      END IF

    ELSE

      ! The variable is not a required variable. Warn about not using it.
      CALL log_warn(RoutineName,                                              &
                    "Provided variable '" // TRIM(var(i)) //                  &
                    "' is not required, so will be ignored")
    END IF

  END DO  !  i (loop over variables)

  !---------------------------------------------------------------------------
  ! Set variables from file
  !---------------------------------------------------------------------------
  IF ( nvars_file > 0 ) THEN

    ! Check that a file name was provided.
    IF ( LEN_TRIM(FILE) == 0 ) THEN
      CALL log_fatal(RoutineName, "No file name provided")
    END IF

    IF ( tpl_has_var_name(FILE) ) THEN
      ! We are using a file name template, so loop through the variables
      ! setting one from each file.
      DO i = 1,nvars_file
        ! Check that a template string was provided for this variable.
        IF ( LEN_TRIM(file_tpl_name(i)) == 0 ) THEN
          CALL log_fatal( RoutineName,                                        &
                          "No variable name template substitution " //        &
                          "provided for " // TRIM(file_var(i)) )
        END IF
        CALL fill_variables_from_file(                                        &
                                 tpl_substitute_var(FILE, file_tpl_name(i)),  &
                                 (/ file_var(i) /), (/ file_var_name(i) /) )
      END DO
    ELSE
      ! We are not using a file name template, so set all variables from the
      ! same file.
      CALL fill_variables_from_file( FILE, file_var(1:nvars_file),            &
                                     file_var_name(1:nvars_file) )
    END IF
  END IF  !  nvars_file

ELSE

  !---------------------------------------------------------------------------
  ! ancil_dump_read%water_resources_props = .TRUE.
  ! Values will be read from the dump file.
  !---------------------------------------------------------------------------
  CALL log_info(RoutineName,                                                  &
                "water resources properties will be read from the dump " //   &
                "file. Namelist values will be ignored.")

END IF ! ancil_dump_read%water_resources_props

END SUBROUTINE init_water_resources_props
#endif
