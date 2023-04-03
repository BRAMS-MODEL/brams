#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_top()

USE io_constants, ONLY: max_file_name_len, max_sdf_name_len, namelist_unit

USE missing_data_mod, ONLY:                                                   &
!  imported scalar parameters
     rmdi

USE string_utils_mod, ONLY: to_string

USE templating_mod, ONLY: tpl_has_var_name, tpl_substitute_var

USE model_interface_mod, ONLY: identifier_len, populate_var, get_var_id

USE jules_hydrology_mod, ONLY: l_top

USE jules_soil_mod, ONLY: l_tile_soil, l_broadcast_ancils

USE dump_mod, ONLY: ancil_dump_read

USE errormessagelength_mod, ONLY: errormessagelength

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in the TOPMODEL properties and checks them for consistency
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
INTEGER, PARAMETER :: max_topmodel_vars = 3
       ! The maximum possible number of TOPMODEL variables that can be given

INTEGER :: nvars_required      ! The number of variables that are
                               ! required in this configuration
CHARACTER(LEN=identifier_len) :: required_vars(max_topmodel_vars)
                               ! The variable identifiers of the required
                               ! variables

INTEGER :: nvars_file       ! The number of variables that will be set
                            ! from the given file (template?)

INTEGER :: error  ! Error indicator

INTEGER :: i  ! Loop counter

!-----------------------------------------------------------------------------
! Definition of the jules_top namelist - this combines local variables with
! some from c_topog
!-----------------------------------------------------------------------------
LOGICAL :: read_from_dump
CHARACTER(LEN=max_file_name_len) :: FILE
                      ! The name of the file (or variable name template) to
                      ! use for variables that need to be filled from file

INTEGER :: nvars      ! The number of variables in this section
CHARACTER(LEN=identifier_len) :: var(max_topmodel_vars)
                      ! The variable identifiers of the variables
LOGICAL :: use_file(max_topmodel_vars)
                      !   T - the variable uses the file
                      !   F - the variable is set using a constant value
CHARACTER(LEN=max_sdf_name_len) :: var_name(max_topmodel_vars)
                      ! The name of each variable in the file
CHARACTER(LEN=max_sdf_name_len) :: tpl_name(max_topmodel_vars)
                      ! The name to substitute in a template for each
                      ! variable
CHARACTER(LEN=errormessagelength) :: iomessage
                      ! I/O error message string
REAL(KIND=real_jlslsm) :: const_val(max_topmodel_vars)
                      ! The constant value to use for each variable if
                      ! use_file = F for that variable
NAMELIST  / jules_top/ read_from_dump, FILE, nvars, var, use_file, var_name,  &
                     tpl_name, const_val


!-----------------------------------------------------------------------------

! If TOPMODEL is not on, we have nothing to do
IF ( .NOT. l_top ) RETURN


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
read_from_dump = .FALSE.
nvars_required = 0
nvars_file     = 0
nvars          = 0
use_file(:)    = .TRUE.  ! Default is for every variable to be read from file
FILE=''      ! Empty file name
var(:)         = ''      ! Empty identifiers.
var_name(:)    = ''      ! Empty variable names.
tpl_name(:)    = ''      ! Empty template string
const_val(:)   = rmdi    ! Missing data value.

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
CALL log_info("init_top", "Reading JULES_TOP namelist...")

READ(namelist_unit, NML = jules_top, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_top",                                                  &
                 "Error reading namelist JULES_TOP " //                       &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

ancil_dump_read%top = read_from_dump

IF ( .NOT. ancil_dump_read%top) THEN !we read from the ancil file
  !---------------------------------------------------------------------------
  ! Set up TOPMODEL properties using namelist values
  !---------------------------------------------------------------------------
  ! Set up the required variables
  ! All the TOPMODEL variables are always required
  nvars_required = max_topmodel_vars
  ! ti_skew may be available in later releases
  required_vars(:) = (/ 'fexp   ', 'ti_mean', 'ti_sig ' /)

  !---------------------------------------------------------------------------
  ! Determine whether to append _soilt as well to tell model_interface_mod
  ! whether the variables read in will have a soil tile dimension.
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
  ! Check that variable identifiers are not empty.
  ! Although we might later decide that the identifier is not required, for
  ! clarity we check here whether the claimed amount of information was
  ! provided.
  !-------------------------------------------------------------------------
  DO i = 1,nvars
    IF ( LEN_TRIM(var(i)) == 0 )                                              &
      CALL log_fatal("init_top",                                              &
                     "Insufficient values for var. " //                       &
                     "No name provided for var at position #" //              &
                     TRIM(to_string(i)) )
  END DO

  !---------------------------------------------------------------------------
  ! Check that all the required variables are there
  !---------------------------------------------------------------------------
  DO i = 1,nvars_required
    IF ( .NOT. ANY(var(1:nvars) == required_vars(i)) )                        &
      CALL log_fatal("init_top",                                              &
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
        CALL log_info("init_top",                                             &
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
          CALL log_fatal("init_top",                                          &
                         "No constant value provided for variable '"          &
                         // TRIM(var(i)) // "'" )

        CALL log_info("init_top",                                             &
                      "'" // TRIM(var(i)) // "' will be set to a " //         &
                      "constant = " // to_string(const_val(i)))

        CALL populate_var(get_var_id(var(i)), const_val = const_val(i))
      END IF
    ELSE
      ! If the variable is not a required variable, warn about not using it
      CALL log_warn("init_top",                                               &
                    "Provided variable '" // TRIM(var(i)) //                  &
                    "' is not required, so will be ignored")
    END IF
  END DO

  !---------------------------------------------------------------------------
  ! Set variables from file
  !---------------------------------------------------------------------------
  IF ( nvars_file > 0 ) THEN
    !   Check that a file name was provided.
    IF ( LEN_TRIM(FILE) == 0 )                                                &
      CALL log_fatal("init_top", "No file name provided")

    IF ( tpl_has_var_name(FILE) ) THEN
      ! We are using a file name template, so loop through the variables setting
      ! one from each file
      DO i = 1,nvars_file
        ! If using a variable name template, check that a template string was
        !provided for the current variable
        IF ( LEN_TRIM(tpl_name(i)) == 0 )                                     &
          CALL log_fatal("init_top",                                          &
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
    END IF
  END IF

ELSE !We read from the dump file
  CALL log_info("init_top",                                                   &
                "topmodel ancils will be read from the dump file.  " //       &
                "Namelist values ignored")

END IF !.NOT. ancil_dump_read%top

RETURN

END SUBROUTINE init_top
#endif
