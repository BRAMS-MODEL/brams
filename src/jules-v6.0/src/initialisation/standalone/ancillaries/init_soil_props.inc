#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_soil_props()

USE io_constants, ONLY: max_sdf_name_len, max_file_name_len, namelist_unit

USE missing_data_mod, ONLY:                                                   &
!  imported scalar parameters
     rmdi

USE string_utils_mod, ONLY: to_string

USE templating_mod, ONLY: tpl_has_var_name, tpl_substitute_var

USE model_interface_mod, ONLY: identifier_len, populate_var, get_var_id

USE dump_mod, ONLY: ancil_dump_read

USE errormessagelength_mod, ONLY: errormessagelength

USE jules_soil_biogeochem_mod, ONLY:                                          &
 ! imported scalar parameters
   soil_model_ecosse, soil_model_rothc,                                       &
 ! imported scalar variables (IN)
   soil_bgc_model

USE jules_soil_mod, ONLY: soil_props_const_z, l_tile_soil, l_broadcast_ancils

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in the soil properties and checks them for consistency
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER :: RoutineName='init_soil_props'

! Work variables
INTEGER, PARAMETER :: max_soil_vars = 11 ! The maximum possible number of
                                         ! soil variables that can be given

INTEGER :: nvars_required      ! The number of soil variables that are
                               ! required in this configuration
CHARACTER(LEN=identifier_len) :: required_vars(max_soil_vars)
                               ! The variable identifiers of the required
                               ! variables

INTEGER :: nvars_file       ! The number of variables that will be set
                            ! from the given file (template?)
! Variables passed to fill_variables_from_file
CHARACTER(LEN=identifier_len) :: file_var(max_soil_vars)
                      ! The variable identifiers of the variables to set
                      ! from file
CHARACTER(LEN=max_sdf_name_len) :: file_var_name(max_soil_vars)
                      ! The name of each variable in the file
CHARACTER(LEN=max_sdf_name_len) :: file_tpl_name(max_soil_vars)
                      ! The name to substitute in a template for each
                      ! variable

INTEGER :: i    ! Index variables

INTEGER :: error  ! Error indicator

!-----------------------------------------------------------------------------
! Definition of the jules_soil_props namelist
!-----------------------------------------------------------------------------
LOGICAL :: read_from_dump
LOGICAL :: const_z            ! T - the same properties are used for each
                              !     soil layer
                              ! F - properties for each layer are read from
                              !     file
CHARACTER(LEN=max_file_name_len) :: FILE
                      ! The name of the file (or variable name template) to
                      ! use for variables that need to be filled from file

INTEGER :: nvars      ! The number of variables in this section
CHARACTER(LEN=identifier_len) :: var(max_soil_vars)
                      ! The variable identifiers of the variables
LOGICAL :: use_file(max_soil_vars)
                      !   T - the variable uses the file
                      !   F - the variable is set using a constant value
CHARACTER(LEN=max_sdf_name_len) :: var_name(max_soil_vars)
                      ! The name of each variable in the file
CHARACTER(LEN=max_sdf_name_len) :: tpl_name(max_soil_vars)
                      ! The name to substitute in a template for each
                      ! variable
CHARACTER(LEN=errormessagelength) :: iomessage
                      ! I/O error message string
REAL(KIND=real_jlslsm) :: const_val(max_soil_vars)
                      ! The constant value to use for each variable if
                      ! use_file = F for that variable
NAMELIST  / jules_soil_props/ read_from_dump, const_z, FILE, nvars, var,      &
                            use_file, var_name, tpl_name, const_val

!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
read_from_dump = .FALSE.
nvars_required = 0
nvars_file     = 0
const_z        = .FALSE. ! Default is to read a value for each soil level
nvars          = 0
use_file(:)    = .TRUE.  ! Default is for every var to be read from file
FILE=''      ! Empty file name.
var(:)         = ''      ! Empty identifiers.
var_name(:)    = ''      ! Empty variable names.
tpl_name(:)    = ''      ! Empty template name.
const_val(:)   = rmdi    ! Missing data value.

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
CALL log_info(RoutineName, "Reading JULES_SOIL_PROPS namelist...")

! Read the soil properties namelist
READ(namelist_unit, NML = jules_soil_props, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal(RoutineName,                                                 &
                 "Error reading namelist JULES_SOIL_PROPS " //                &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

ancil_dump_read%soil_props = read_from_dump

IF ( .NOT. ancil_dump_read%soil_props) THEN
  !-------------------------------------------------------------------------
  ! Read from the ancil file.
  !-------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  ! Set soil properties using namelist values.
  !---------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! Check that variable identifiers are not empty.
  ! Although we might later decide that the identifier is not required, for
  ! clarity we check here whether the claimed amount of information was
  ! provided.
  !-------------------------------------------------------------------------
  DO i = 1,nvars
    IF ( LEN_TRIM(var(i)) == 0 )                                              &
      CALL log_fatal(RoutineName,                                             &
                     "Insufficient values for var. " //                       &
                     "No name provided for var at position #" //              &
                     TRIM(to_string(i)) )
  END DO

  ! Set up the required variables

  ! First, set up the variables that are always required.
  nvars_required = 9
  required_vars(1:nvars_required) = (/ 'b      ', 'sathh  ', 'satcon ',       &
                                       'sm_sat ', 'sm_crit', 'sm_wilt',       &
                                       'hcap   ', 'hcon   ', 'albsoil' /)

  ! If RothC is selected, clay is required.
  IF ( soil_bgc_model == soil_model_rothc ) THEN
    nvars_required                = nvars_required + 1
    required_vars(nvars_required) = 'clay'
    ! In older versions the soil clay content was set to zero instead of
    ! being read in. If it is not available in ancillaries give a stern
    ! warning and set it to zero.
    IF ( .NOT. ANY(var(1:nvars) == 'clay') ) THEN
      CALL log_warn(RoutineName,                                              &
                    "No value given for soil clay content. "            //    &
                    "Soil clay content is required with RothC model. "  //    &
                    "It will be set to 0.0 as for previous versions. "  //    &
                    "This is WRONG - please try and find values for clay.")
      ! Add clay to the list of variables, so that it will be set to zero.
      nvars                     = nvars + 1
      var(nvars)                = 'clay'
      use_file(nvars_required)  = .FALSE.
      const_val(nvars_required) = 0.0
    END IF
  END IF

  ! Variables used with ECOSSE.
  IF ( soil_bgc_model == soil_model_ecosse ) THEN
    nvars_required                = nvars_required + 1
    required_vars(nvars_required) = 'clay'
    nvars_required                = nvars_required + 1
    required_vars(nvars_required) = 'soil_ph'
  END IF

  !-------------------------------------------------------------------------
  ! Check that all the required variables are provided.
  !-------------------------------------------------------------------------
  DO i = 1,nvars_required
    IF ( .NOT. ANY(var(1:nvars) == required_vars(i)) )                        &
      CALL log_fatal(RoutineName,                                             &
                     "No value given for required variable '" //              &
                     TRIM(required_vars(i)) // "'")
  END DO

  !---------------------------------------------------------------------------
  ! Determine whether to append _soilt as well to tell model_interface_mod
  ! that the variables read in will have a soil tile dimension.
  ! This is required when l_tile_soil = T, and l_broadcast_ancils = F.
  !---------------------------------------------------------------------------
  IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils ) THEN
    DO i = 1,nvars
      var(i) = TRIM(var(i)) // "_soilt"
    END DO

    DO i = 1,nvars_required
      required_vars(i) = TRIM(required_vars(i)) // "_soilt"
    END DO
  END IF

  !-------------------------------------------------------------------------
  ! Constant Z (i.e. spatially varying but constant through vertical levels)
  ! is implemented by having a separate input variable in
  ! model_interface_mod called <var>_const_z that has no vertical levels.
  ! Hence, once the previous check is done, we add _const_z to both
  ! required and provided variable identifiers if asked for.
  !-------------------------------------------------------------------------
  soil_props_const_z = const_z
  IF ( soil_props_const_z ) THEN
    DO i = 1,nvars
      ! Don't change variables that do not have multiple levels.
      IF ( var(i) /= 'albsoil' ) THEN
        var(i) = TRIM(var(i)) // "_const_z"
      END IF
    END DO

    DO i = 1,nvars_required
      ! Don't change variables that do not have multiple levels.
      IF ( required_vars(i) /= 'albsoil' ) THEN
        required_vars(i) = TRIM(required_vars(i)) // "_const_z"
      END IF
    END DO
  END IF  !  soil_ancil_const_z

  !-------------------------------------------------------------------------
  ! Check which variables we will be using and partition them into variables
  ! set to constant values and variables set from file
  !-------------------------------------------------------------------------
  DO i = 1,nvars

    !-------------------------------------------------------------------------
    ! If the variable is one of the required vars, then we will be using it
    !-------------------------------------------------------------------------
    IF ( ANY( required_vars(1:nvars_required) == TRIM(var(i)) ) ) THEN
      IF ( use_file(i) ) THEN
        CALL log_info(RoutineName,                                            &
                      "'" // TRIM(var(i)) // "' will be read from file")

        ! If the variable will be filled from file, register it here
        nvars_file = nvars_file + 1
        file_var(nvars_file) = var(i)
        file_var_name(nvars_file) = var_name(i)
        file_tpl_name(nvars_file) = tpl_name(i)
      ELSE
        ! If the variable is being set as a constant, populate it here.
        ! First check that a value has been provided.
        IF ( ABS( const_val(i) - rmdi ) < EPSILON(1.0) )                      &
          CALL log_fatal(RoutineName,                                         &
                         "No constant value provided for variable '"          &
                         // TRIM(var(i)) // "'" )

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

  !-------------------------------------------------------------------------
  ! Set variables from file
  !-------------------------------------------------------------------------
  IF ( nvars_file > 0 ) THEN
    ! Check that a file name was provided.
    IF ( LEN_TRIM(FILE) == 0 )                                                &
      CALL log_fatal(RoutineName, "No file name provided")

    IF ( tpl_has_var_name(FILE) ) THEN
      ! We are using a file name template, so loop through the variables
      ! setting one from each file.
      DO i = 1,nvars_file
        ! Check that a template string was provided for this variable.
        IF ( LEN_TRIM(file_tpl_name(i)) == 0 )                                &
          CALL log_fatal( RoutineName,                                        &
                          "No variable name template substitution " //        &
                          "provided for " // TRIM(file_var(i)) )
        CALL fill_variables_from_file(                                        &
          tpl_substitute_var(FILE, file_tpl_name(i)),                         &
          (/ file_var(i) /), (/ file_var_name(i) /)                           &
        )
      END DO
    ELSE
      ! We are not using a file name template, so set all variables from the
      ! same file.
      CALL fill_variables_from_file(                                          &
        FILE,file_var(1:nvars_file), file_var_name(1:nvars_file)              &
      )
    END IF
  END IF

ELSE

  !-------------------------------------------------------------------------
  ! ancil_dump_read%soil_props = .TRUE.
  ! Values will be read from the dump file.
  !-------------------------------------------------------------------------
  CALL log_info(RoutineName,                                                  &
                "soil properties will be read from the dump file.  " //       &
                "Namelist values ignored")

END IF  ! ancil_dump_read%soil_props

RETURN

END SUBROUTINE init_soil_props
#endif
