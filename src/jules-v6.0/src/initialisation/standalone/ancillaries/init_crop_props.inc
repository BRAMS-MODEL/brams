#if !defined(UM_JULES)

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Reads in JULES_CROP_PROPS namelist if necessary

SUBROUTINE init_crop_props()

USE io_constants, ONLY: max_file_name_len, max_sdf_name_len, namelist_unit

USE missing_data_mod, ONLY:                                                   &
!  imported scalar parameters
     rmdi

USE templating_mod, ONLY: tpl_has_var_name, tpl_substitute_var

USE model_interface_mod, ONLY: identifier_len, populate_var, get_var_id

USE string_utils_mod, ONLY: to_string

USE jules_vegetation_mod, ONLY: l_crop, l_prescsow, l_croprotate

USE dump_mod, ONLY: ancil_dump_read

USE errormessagelength_mod, ONLY: errormessagelength

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!    Reads in TTveg, TTrep
!    Reads in the sowing date for each of the crop pfts if l_prescsow=T
!    and the latest possible harvest date for the crop pfts if provided
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
INTEGER, PARAMETER :: max_crop_vars = 4
       ! The maximum possible number of crop variables that can be given
       ! includes the croplatestharvdate which is only a required variable
       ! when l_croprotate=T and optional when l_croprotate=F.

INTEGER :: nvars_required      ! The number of variables that are
                               ! required in this configuration

CHARACTER(LEN=identifier_len), ALLOCATABLE :: required_vars(:)
                               ! The variable identifiers of the required
                               ! variables

INTEGER :: nvars_file       ! The number of variables that will be set
                            ! from the given file (template?)

INTEGER :: error  ! Error indicator

INTEGER :: i  ! Loop counter

!-----------------------------------------------------------------------------
! Definition of the jules_crop_props namelist
!-----------------------------------------------------------------------------
LOGICAL :: read_from_dump
CHARACTER(LEN=max_file_name_len) :: FILE
                      ! The name of the file (or variable name template) to
                      ! use for variables that need to be filled from file
INTEGER :: nvars      ! The number of variables in this section
CHARACTER(LEN=identifier_len) :: var(max_crop_vars)
                      ! The variable identifiers of the variables
CHARACTER(LEN=errormessagelength) :: iomessage
                      ! Error message string for I/O errors
LOGICAL :: use_file(max_crop_vars)
                      !   T - the variable uses the file
                      !   F - the variable is set using a constant value
CHARACTER(LEN=max_sdf_name_len) :: var_name(max_crop_vars)
                      ! The name of each variable in the file
CHARACTER(LEN=max_sdf_name_len) :: tpl_name(max_crop_vars)
                      ! The name to substitute in a template for each
                      ! variable
REAL(KIND=real_jlslsm) :: const_val(max_crop_vars)
                      ! The constant value to use for each variable if
                      ! use_file = F for that variable
NAMELIST  / jules_crop_props/                                                 &
                     read_from_dump, FILE, nvars, var, use_file, var_name,    &
                     tpl_name, const_val

!-----------------------------------------------------------------------------

! Nothing to do if crop model is not on
IF ( .NOT. l_crop ) RETURN

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
CALL log_info("init_crop_props", "Reading JULES_CROP_PROPS namelist...")

READ(namelist_unit, NML = jules_crop_props, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_crop_props",                                           &
                "Error reading namelist JULES_CROP_PROPS " //                 &
                "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //          &
                TRIM(iomessage) // ")")

ancil_dump_read%crop_props = read_from_dump

IF ( .NOT. ancil_dump_read%crop_props) THEN !we read from the ancil file

  !-------------------------------------------------------------------------
  ! Set up crop properties using namelist values
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! Check that variable identifiers are not empty.
  ! Although we might later decide that the identifier is not required, for
  ! clarity we check here whether the claimed amount of information was
  ! provided.
  !-------------------------------------------------------------------------
  DO i = 1,nvars
    IF ( LEN_TRIM(var(i)) == 0 )                                              &
      CALL log_fatal("init_crop_props",                                       &
                     "Insufficient values for var. " //                       &
                     "No name provided for var at position #" //              &
                     TRIM(to_string(i)) )
  END DO

  ! Set up the required variables

  IF (l_prescsow) THEN
    IF ( (l_croprotate) .AND.                                                 &
         .NOT. (ANY(var(1:nvars) == 'croplatestharvdate')) )                  &

      CALL log_fatal("init_crop_props",                                       &
                     "If l_croprotate is T then value for "    //             &
                     "latestharvdate must be provided"         //             &
                     TRIM(to_string(i)) )

    IF (ANY(var(1:nvars) == 'croplatestharvdate')) THEN
      nvars_required = max_crop_vars
      ALLOCATE(required_vars(nvars_required))
      required_vars(:) = (/ 'cropsowdate       ', 'cropttveg         ',       &
                            'cropttrep         ', 'croplatestharvdate'     /)
    ELSE
      nvars_required = max_crop_vars - 1
      ALLOCATE(required_vars(nvars_required))
      required_vars(:) = (/ 'cropsowdate', 'cropttveg  ', 'cropttrep  ' /)
    END IF
  ELSE
    nvars_required = max_crop_vars - 2
    ALLOCATE(required_vars(nvars_required))
    required_vars(:) = (/ 'cropttveg', 'cropttrep' /)
  END IF

  !---------------------------------------------------------------------------
  ! Check that all the required variables are there
  !---------------------------------------------------------------------------
  DO i = 1,nvars_required
    IF ( .NOT. ANY(var(1:nvars) == required_vars(i)) )                        &
      CALL log_fatal("init_crop_props",                                       &
                     "No value given for required variable '" //              &
                     TRIM(required_vars(i)) // "'")
  END DO

  !---------------------------------------------------------------------------
  ! Check which variables we will be using and partition them into variables
  ! set to constant values and variables set from file
  !---------------------------------------------------------------------------
  DO i = 1,nvars
    !-------------------------------------------------------------------------
    ! If the variable is one of the required vars, then we will be using it
    !-------------------------------------------------------------------------
    IF ( ANY(required_vars(1:nvars_required) == var(i)) ) THEN
      IF ( use_file(i) ) THEN
        CALL log_info("init_crop_props",                                      &
                      "'" // TRIM(var(i)) // "' will be read from file")

        ! If the variable will be filled from file, register it here
        nvars_file = nvars_file + 1
        ! Since nvars_file <= i (so we will not overwrite unprocessed values)
        ! and we do not need the values from these arrays for any non-file
        ! variables from now on, we can just compress them down onto variables
        !that are in the file
        var(nvars_file) = var(i)
        var_name(nvars_file) = var_name(i)
        tpl_name(nvars_file) = tpl_name(i)
      ELSE
        ! If the variable is being set as a constant, populate it here.
        ! First check that a value has been provided.
        IF ( ABS( const_val(i) - rmdi ) < EPSILON(1.0) )                      &
          CALL log_fatal("init_crop_props",                                   &
                         "No constant value provided for variable '"          &
                         // TRIM(var(i)) // "'" )

        CALL log_info("init_crop_props",                                      &
                      "'" // TRIM(var(i)) // "' will be set to a " //         &
                      "constant = " // to_string(const_val(i)))

        CALL populate_var(get_var_id(var(i)), const_val = const_val(i))
      END IF
    ELSE
      ! If the variable is not a required variable, warn about not using it
      CALL log_warn("init_crop_props",                                        &
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
      CALL log_fatal("init_crop_props", "No file name provided")

    IF ( tpl_has_var_name(FILE) ) THEN
      ! We are using a file name template, so loop through the variables
      !setting one from each file
      DO i = 1,nvars_file
        ! If using a variable name template, check that a template string was
        !provided for the current variable
        IF ( LEN_TRIM(tpl_name(i)) == 0 )                                     &
          CALL log_fatal("init_crop_props",                                   &
                         "No variable name template substitution " //         &
                         "(tpl_name) provided for " // TRIM(var(i)))

        CALL fill_variables_from_file(                                        &
          tpl_substitute_var(FILE, tpl_name(i)),                              &
          (/ var(i) /), (/ var_name(i) /)                                     &
        )
      END DO
    ELSE
      ! We are not using a file name template, so set all variables from the
      ! same file
      CALL fill_variables_from_file(                                          &
        FILE, var(1:nvars_file), var_name(1:nvars_file)                       &
      )
    END IF
  END IF

  DEALLOCATE(required_vars)

ELSE !We read from the dump file
  CALL log_info("init_crop_props",                                            &
                "crop properties will be read from the dump file.  " //       &
                "Namelist values ignored")

END IF !.NOT. ancil_dump_read%crop_props

END SUBROUTINE init_crop_props
#endif
