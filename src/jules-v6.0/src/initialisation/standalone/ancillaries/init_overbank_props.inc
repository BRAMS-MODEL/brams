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

SUBROUTINE init_overbank_props()

USE mpi, ONLY: mpi_comm_world

USE io_constants, ONLY: max_sdf_name_len, max_file_name_len, namelist_unit

USE missing_data_mod, ONLY: rmdi

USE ancil_info, ONLY: land_pts

USE string_utils_mod, ONLY: to_string

USE templating_mod, ONLY: tpl_has_var_name, tpl_substitute_var

USE model_interface_mod, ONLY: identifier_len, populate_var, get_var_id

USE model_grid_mod, ONLY: global_land_pts

USE jules_rivers_mod, ONLY:                                                   &
  np_rivers, nx_rivers, ny_rivers, l_rivers

USE overbank_inundation_mod, ONLY:                                            &
  l_riv_overbank, l_riv_hypsometry,                                           &
  logn_mean, logn_stdev, logn_mean_rp, logn_stdev_rp,                         &
  frac_fplain_lp, frac_fplain_rp,                                             &
  qbf, dbf, wbf

USE logging_mod, ONLY: log_info

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in river overbank inundation parameters and properties, and
!   allocates various overbank variables.
!
!  Code Owner: Please refer to ModuleLeaders.txt
!  This file belongs in section: Hydrology
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER :: RoutineName='init_overbank_props'

! Work variables

INTEGER ::                                                                    &
  error        = 0,                                                           &
                     ! Variable for trapping the error from each
                     ! individual call to allocate
  error_sum    = 0
                     ! Variable to track the sum of all errors
                     ! resulting from calls to allocate. Hence we
                     ! know that everything was successful if and
                     ! only if this is zero at the end

INTEGER, PARAMETER :: max_inundation_vars = 2 ! The maximum possible number
                                              ! of inundation variables

INTEGER :: nvars_required      ! The number of inundation variables that are
                               ! required in this configuration
INTEGER :: nvars_optional      ! The number of optional inundation variables
                               ! in this configuration
CHARACTER(LEN=identifier_len) :: required_vars(max_inundation_vars)
                               ! The variable identifiers of the required
                               ! variables
CHARACTER(LEN=identifier_len) :: optional_vars(max_inundation_vars)
                               ! The variable identifiers of any optional
                               ! variables

INTEGER :: nvars_file          ! The number of variables that will be set
                               ! from the given file (template?)

! Variables passed to fill_variables_from_file

CHARACTER(LEN=identifier_len) :: file_var(max_inundation_vars)
                      ! The variable identifiers of the variables to set
                      ! from file
CHARACTER(LEN=max_sdf_name_len) :: file_var_name(max_inundation_vars)
                      ! The name of each variable in the file
CHARACTER(LEN=max_sdf_name_len) :: file_tpl_name(max_inundation_vars)
                      ! The name to substitute in a template for each
                      ! variable

CHARACTER(LEN=max_file_name_len) :: FILE
                      ! The name of the file (or variable name template) to
                      ! use for variables that need to be filled from file

INTEGER :: nvars      ! The number of variables in this section
CHARACTER(LEN=identifier_len) :: var(max_inundation_vars)
                      ! The variable identifiers of the variables
LOGICAL :: use_file(max_inundation_vars)
                      !   T - the variable uses the file
                      !   F - the variable is set using a constant value
CHARACTER(LEN=max_sdf_name_len) :: var_name(max_inundation_vars)
                      ! The name of each variable in the file
CHARACTER(LEN=max_sdf_name_len) :: tpl_name(max_inundation_vars)
                      ! The name to substitute in a template for each
                      ! variable
REAL(KIND=real_jlslsm) :: const_val(max_inundation_vars)
                      ! The constant value to use for each variable if
                      ! use_file = F for that variable

INTEGER :: i  ! Index variables

NAMELIST  / jules_overbank_props/                                             &
   FILE, nvars, var, use_file, var_name, tpl_name, const_val

! If rivers or overbank inundation are not on, there is nothing to do
IF (l_rivers .AND. l_riv_overbank) THEN

  !---------------------------------------------------------------------------
  ! Initialise
  !---------------------------------------------------------------------------

  nvars_required = 0
  nvars_optional = 0
  nvars_file     = 0
  nvars          = 0
  use_file(:)    = .TRUE. ! Default is for every variable to be read from file
  FILE=''     ! Empty file name.

  IF ( global_land_pts <= 1 ) THEN
    l_rivers = .FALSE.
    CALL log_warn(RoutineName,                                                &
                  "River routing not appropriate for single point runs. " //  &
                  "Overbank inundation disabled.")

  ELSE

    !-------------------------------------------------------------------------
    ! Read overbank inundation namelist
    !--------------------------------------------------------------------------

    CALL log_info(RoutineName,"Reading JULES_OVERBANK_PROPS namelist...")

    READ(namelist_unit, NML = jules_overbank_props, IOSTAT = error)
    IF ( error /= 0 ) CALL log_fatal(RoutineName,                             &
                  "Error reading namelist JULES_OVERBANK_PROPS " //           &
                  "(IOSTAT=" // TRIM(to_string(error)) // ")")

    !--------------------------------------------------------------------------
    ! Set up river overbank inundation properties
    !--------------------------------------------------------------------------

    IF ( l_riv_hypsometry ) THEN
      nvars_required = 2
      required_vars(1:nvars_required) = (/ 'logn_mean ', 'logn_stdev' /)

      !-----------------------------------------------------------------------
      ! Check that all the required variables are there
      !-----------------------------------------------------------------------
      DO i = 1,nvars_required
        IF ( .NOT. ANY(var(1:nvars) == required_vars(i)) )                    &
        CALL log_fatal(RoutineName,                                           &
                       "No value given for required variable '" //            &
                       TRIM(required_vars(i)) // "'")
      END DO

      !-----------------------------------------------------------------------
      ! Allocate inundation-specific arrays from ancillary information
      !-----------------------------------------------------------------------

      IF ( nx_rivers > 0 .AND. ny_rivers > 0) THEN

        ALLOCATE( logn_mean(nx_rivers, ny_rivers), stat = error )
        error_sum = error
        ALLOCATE( logn_stdev(nx_rivers, ny_rivers), stat = error )
        error_sum = error_sum + error
        IF ( error_sum == 0 ) THEN
          ! Initialise to impossible values.
          logn_mean(:,:)  = rmdi
          logn_stdev(:,:) = rmdi
        ELSE
          CALL log_fatal(RoutineName,                                         &
                         "Error allocating arrays from ancillaries")
        END IF

      ELSE

        CALL log_fatal(RoutineName,                                           &
                  "No overbank inundation with invalid routing dimensions" )
      END IF
    END IF

    !-------------------------------------------------------------------------
    ! Check which variables we will be using and partition them into variables
    ! set to constant values and variables set from file
    !-------------------------------------------------------------------------
    DO i = 1,nvars

      !-----------------------------------------------------------------------
      ! If the variable is one of the required vars, then we will be using it
      !-----------------------------------------------------------------------

      IF ( ANY(required_vars(1:nvars_required) == TRIM(var(i)))               &
           .OR. (nvars_optional > 0 .AND.                                     &
           ANY(optional_vars(1:nvars_optional) == TRIM(var(i))) )) THEN
        IF ( use_file(i) ) THEN
          CALL log_info(RoutineName,                                          &
                        "'" // TRIM(var(i)) // "' will be read from file")

          ! If the variable will be filled from file, register it here
          nvars_file = nvars_file + 1
          file_var(nvars_file) = TRIM(var(i))
          file_var_name(nvars_file) = var_name(i)
          file_tpl_name(nvars_file) = TRIM(tpl_name(i))

        ELSE

          ! If the variable is being set as a constant, just populate it here
          CALL log_info(RoutineName,                                          &
                        "'" // TRIM(var(i)) // "' will be set to a " //       &
                        "constant = " // to_string(const_val(i)))

          CALL populate_var(get_var_id(var(i)), const_val = const_val(i))
        END IF

      ELSE

        ! If the variable is not a required variable, warn about not using it
        CALL log_warn(RoutineName,                                            &
                      "Provided variable '" // TRIM(var(i)) //                &
                      "' is not required, so will be ignored")
      END IF

    END DO ! nvars

    !-------------------------------------------------------------------------
    ! Set variables from file
    !-------------------------------------------------------------------------

    IF ( nvars_file > 0 ) THEN

      ! Check that a file name was provided
      IF ( LEN_TRIM(FILE) == 0 ) THEN
        CALL log_fatal(RoutineName, "No file name provided")
      END IF

      CALL log_info(RoutineName,                                              &
                    "Reading overbank inundation information from file " //   &
                    TRIM(FILE))

      IF ( tpl_has_var_name(FILE) ) THEN
        ! We are using a file name template, so loop through the variables
        ! setting one from each file.
        DO i = 1,nvars_file
          CALL fill_variables_from_file(tpl_substitute_var(FILE,              &
                                        file_tpl_name(i)),                    &
                                        (/ file_var(i) /),                    &
                                        (/ file_var_name(i) /) )
        END DO

      ELSE
        ! We are not using a file name template, so set all variables from
        ! the same file.
        CALL fill_variables_from_file(                                        &
          FILE,file_var(1:nvars_file), file_var_name(1:nvars_file))

      END IF
    END IF

    !-------------------------------------------------------------------------
    ! Allocate inundation variables defined on river routing points
    !-------------------------------------------------------------------------

    ALLOCATE(logn_mean_rp(np_rivers),       stat = error)
    error_sum = error
    ALLOCATE(logn_stdev_rp(np_rivers),      stat = error)
    error_sum = error_sum + error
    ALLOCATE(qbf(np_rivers),                stat = error)
    error_sum = error_sum + error
    ALLOCATE(dbf(np_rivers),                stat = error)
    error_sum = error_sum + error
    ALLOCATE(wbf(np_rivers),                stat = error)
    error_sum = error_sum + error
    ALLOCATE(frac_fplain_rp(np_rivers), stat = error)
    error_sum = error_sum + error
    IF ( error_sum == 0 ) THEN
      ! Initialise array values
      logn_mean_rp(:)       = rmdi
      logn_stdev_rp(:)      = rmdi
      qbf(:)                = rmdi
      dbf(:)                = rmdi
      wbf(:)                = rmdi
      frac_fplain_rp(:)     = rmdi
    ELSE
      CALL log_fatal(RoutineName,                                             &
                     "Error allocating river grid arrays")
    END IF

    !-------------------------------------------------------------------------
    ! Allocate inundation variables defined on land points
    !-------------------------------------------------------------------------

    ALLOCATE(frac_fplain_lp(land_pts), stat = error)
    error_sum = error
    IF ( error_sum == 0 ) THEN
      ! Initialise array values
      frac_fplain_lp(:) = 0.0
    ELSE
      CALL log_fatal(RoutineName,                                             &
                     "Error allocating land grid arrays")
    END IF

  END IF ! (check global_landpoint > 1)

END IF

END SUBROUTINE init_overbank_props
#endif
