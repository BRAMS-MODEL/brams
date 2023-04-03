#if !defined(UM_JULES)
!******************************COPYRIGHT**************************************
! (c) Centre for Ecology and Hydrology. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT**************************************

SUBROUTINE init_pdm()

USE io_constants, ONLY: max_file_name_len, max_sdf_name_len, namelist_unit

USE string_utils_mod, ONLY: to_string

USE templating_mod, ONLY: tpl_has_var_name, tpl_substitute_var

USE model_interface_mod, ONLY: identifier_len, populate_var, get_var_id

USE jules_hydrology_mod, ONLY: l_spdmvar

USE errormessagelength_mod, ONLY: errormessagelength

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in the PDM properties and checks them for consistency
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
INTEGER, PARAMETER :: max_pdm_vars = 1
       ! The maximum possible number of PDM variables that can be given

INTEGER :: nvars_required      ! The number of variables that are
                               ! required in this configuration
CHARACTER(LEN=identifier_len) :: required_vars(max_pdm_vars)
                               ! The variable identifiers of the required
                               ! variables

INTEGER :: nvars_file       ! The number of variables that will be set
                            ! from the given file (template?)

INTEGER :: error  ! Error indicator

INTEGER :: i  ! Loop counter

!-----------------------------------------------------------------------------
! Definition of the jules_pdm namelist - this combines local variables with
! some from c_topog
!-----------------------------------------------------------------------------
LOGICAL :: read_from_dump
CHARACTER(LEN=max_file_name_len) :: FILE
                      ! The name of the file (or variable name template) to
                      ! use for variables that need to be filled from file

INTEGER :: nvars      ! The number of variables in this section
CHARACTER(LEN=identifier_len) :: var(max_pdm_vars)
                      ! The variable identifiers of the variables
LOGICAL :: use_file(max_pdm_vars)
                      !   T - the variable uses the file
                      !   F - the variable is set using a constant value
CHARACTER(LEN=max_sdf_name_len) :: var_name(max_pdm_vars)
                      ! The name of each variable in the file
CHARACTER(LEN=max_sdf_name_len) :: tpl_name(max_pdm_vars)
                      ! The name to substitute in a template for each
                      ! variable
CHARACTER(LEN=errormessagelength) :: iomessage
                      ! I/O error message string
REAL(KIND=real_jlslsm) :: const_val(max_pdm_vars)
                      ! The constant value to use for each variable if
                      ! use_file = F for that variable
NAMELIST  / jules_pdm/ read_from_dump, FILE, nvars, var, use_file, var_name,  &
                     tpl_name, const_val


!-----------------------------------------------------------------------------

! If s_pdm as a slope dependent parameter is not on, we have nothing to do
IF ( .NOT. l_spdmvar ) RETURN


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
read_from_dump = .FALSE.
nvars_required = 0
nvars_file     = 0
nvars          = 0
use_file(:)    = .TRUE.  ! Default is for every variable to be read from file
FILE=''      ! Empty file name
tpl_name(:)    = ''      ! Empty template string

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
CALL log_info("init_pdm", "Reading JULES_PDM namelist...")

READ(namelist_unit, NML = jules_pdm, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_pdm",                                                  &
                 "Error reading namelist JULES_PDM " //                       &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

IF ( .NOT. read_from_dump) THEN !we read from the ancil file
  !---------------------------------------------------------------------------
  ! Set up PDM properties using namelist values
  !---------------------------------------------------------------------------
  ! Set up the required variables
  ! All the PDM variables are always required
  nvars_required = max_pdm_vars
  required_vars(:) = (/ 'slope' /)


  !---------------------------------------------------------------------------
  ! Check that all the required variables are there
  !---------------------------------------------------------------------------
  DO i = 1,nvars_required
    IF ( .NOT. ANY(var(1:nvars) == required_vars(i)) )                        &
      CALL log_fatal("init_pdm",                                              &
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
        CALL log_info("init_pdm",                                             &
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
        ! If the variable is being set as a constant, just populate it here
        CALL log_info("init_pdm",                                             &
                      "'" // TRIM(var(i)) // "' will be set to a " //         &
                      "constant = " // to_string(const_val(i)))

        CALL populate_var(get_var_id(var(i)), const_val = const_val(i))
      END IF
    ELSE
      ! If the variable is not a required variable, warn about not using it
      CALL log_warn("init_pdm",                                               &
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
      CALL log_fatal("init_pdm", "No file name provided")

    IF ( tpl_has_var_name(FILE) ) THEN
      ! We are using a file name template, so loop through the variables setting
      ! one from each file
      DO i = 1,nvars_file
        ! If using a variable name template, check that a template string was
        !provided for the current variable
        IF ( LEN_TRIM(tpl_name(i)) == 0 )                                     &
          CALL log_fatal("init_pdm",                                          &
                         "No variable name template substitution " //         &
                         "(tpl_name) provided for " // TRIM(var(i)))

        CALL fill_variables_from_file(                                        &
                                      tpl_substitute_var(FILE, tpl_name(i)),  &
                                      (/ var(i) /), (/ var_name(i) /))
      END DO
    ELSE
      ! We are not using a file name template, so set all variables from the same
      ! file
      CALL fill_variables_from_file(FILE, var(1:nvars_file),                  &
                                    var_name(1:nvars_file))
    END IF
  END IF

ELSE !We read from the dump file
  CALL log_info("init_pdm",                                                   &
                "pdm ancils will be read from the dump file.  " //            &
                "Namelist values ignored")

END IF !.NOT. ancil_dump_read%frac

RETURN

END SUBROUTINE init_pdm
#endif
