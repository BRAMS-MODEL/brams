#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_urban_props(ainfo,urban_param)

  !Use in subroutines
USE tilepts_mod, ONLY: tilepts

!Use in variables
USE io_constants, ONLY: max_sdf_name_len, max_file_name_len, namelist_unit

USE missing_data_mod, ONLY:                                                   &
!  imported scalar parameters
     rmdi

USE string_utils_mod, ONLY: to_string

USE templating_mod, ONLY: tpl_has_var_name, tpl_substitute_var

USE model_interface_mod, ONLY: identifier_len, populate_var, get_var_id

USE input_mod, ONLY: fill_variables_from_file

USE ancil_info, ONLY: land_pts, surft_pts

USE jules_surface_types_mod, ONLY: ice, urban, urban_canyon

USE switches_urban, ONLY: l_urban_empirical, l_moruses_macdonald, l_moruses

USE urban_param_mod, ONLY: cdz, kappa2, a, z0m_mat

USE logging_mod, ONLY: log_info, log_warn, log_fatal

USE errormessagelength_mod, ONLY: errormessagelength

USE um_types, ONLY: real_jlslsm

!TYPE definitions
USE ancil_info,    ONLY: ainfo_type
USE urban_param_mod, ONLY: urban_param_type

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises urban parameters and properties
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

!Arguments

TYPE(ainfo_type), INTENT(IN OUT) :: ainfo
TYPE(urban_param_type), INTENT(IN OUT) :: urban_param

! Work variables
INTEGER, PARAMETER :: max_urban_vars = 9
       ! The maximum possible number of TOPMODEL variables that can be given

INTEGER :: nvars_required      ! The number of variables that are
                               ! required in this configuration
CHARACTER(LEN=identifier_len) :: required_vars(max_urban_vars)
                               ! The variable identifiers of the required
                               ! variables

INTEGER :: nvars_file       ! The number of variables that will be set
                            ! from the given file (template?)

REAL(KIND=real_jlslsm) :: sc_hwr(land_pts), d_h(land_pts)  ! Work variables
REAL(KIND=real_jlslsm) :: lambdaf, lambdap

INTEGER :: i,l  ! Index variables

INTEGER :: error  ! Error indicator
CHARACTER(LEN=errormessagelength) :: iomessage

!-----------------------------------------------------------------------------
! Definition of the urban_properties namelist - this specifies how urban
! properties are set
!-----------------------------------------------------------------------------
CHARACTER(LEN=max_file_name_len) :: FILE
                      ! The name of the file (or variable name template) to
                      ! use for variables that need to be filled from file

INTEGER :: nvars      ! The number of variables in this section
CHARACTER(LEN=identifier_len) :: var(max_urban_vars)
                      ! The variable identifiers of the variables
LOGICAL :: use_file(max_urban_vars)
                      !   T - the variable uses the file
                      !   F - the variable is set using a constant value
CHARACTER(LEN=max_sdf_name_len) :: var_name(max_urban_vars)
                      ! The name of each variable in the file
CHARACTER(LEN=max_sdf_name_len) :: tpl_name(max_urban_vars)
                      ! The name to substitute in a template for each
                      ! variable
REAL(KIND=real_jlslsm) :: const_val(max_urban_vars)
                      ! The constant value to use for each variable if
                      ! use_file = F for that variable
CHARACTER(LEN=*), PARAMETER :: RoutineName='init_urban_props'

NAMELIST  /urban_properties/ FILE, nvars, var, use_file, var_name, tpl_name,  &
                            const_val

!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
nvars_required = 0
nvars_file     = 0
nvars          = 0
use_file(:)    = .TRUE.  ! Default is to read every variable from file
FILE=''      ! Empty file name.
var(:)         = ''      ! Empty identifiers.
var_name(:)    = ''      ! Empty variable names.
tpl_name(:)    = ''      ! Empty template name.
const_val(:)   = rmdi    ! Missing data value.

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
CALL log_info(RoutineName, "Reading URBAN_PROPERTIES namelist...")
READ(namelist_unit, NML = urban_properties, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal(RoutineName,                                                 &
                 "Error reading namelist URBAN_PROPERTIES " //                &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

!-----------------------------------------------------------------------------
! Process the namelist values and set derived variables
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Check that the run actually has urban and that urban schemes are not run
! in error
!-----------------------------------------------------------------------------
CALL tilepts(land_pts, ainfo%frac_surft, surft_pts, ainfo%surft_index,        &
             ainfo%l_lice_point)
IF ( surft_pts(urban_canyon) == 0 )                                           &
  CALL log_warn(RoutineName,                                                  &
                "URBAN-2T or MORUSES is selected but there are no " //        &
                "urban land points - extra calculations may be being " //     &
                "performed that will not impact on results")

!-----------------------------------------------------------------------------
! Read values for urban properties
!-----------------------------------------------------------------------------
! Set up the required variables
! First, variables that are always required and can't be derived
nvars_required = 4
required_vars(1:nvars_required) = (/ 'albwl', 'albrd', 'emisw', 'emisr' /)

IF ( l_moruses ) THEN
  IF ( .NOT. l_urban_empirical ) THEN
    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'wrr'
    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'hwr'
    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'hgt'
  END IF
  IF ( .NOT. l_moruses_macdonald ) THEN
    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'ztm'
    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'disp'
  END IF
ELSE
  ! For urban2t, we only need wrr
  nvars_required = 1
  required_vars(1) = 'wrr'
END IF

!-----------------------------------------------------------------------------
! Check that variable identifiers are not empty.
! Although we might later decide that the identifier is not required, for
! clarity we check here whether the claimed amount of information was
! provided.
!-----------------------------------------------------------------------------
DO i = 1,nvars
  IF ( LEN_TRIM(var(i)) == 0 )                                                &
    CALL log_fatal(RoutineName,                                               &
                   "Insufficient values for var. " //                         &
                   "No name provided for var at position #" //                &
                   TRIM(to_string(i)) )
END DO

!-----------------------------------------------------------------------------
! Check that all the required variables are there
!-----------------------------------------------------------------------------
DO i = 1,nvars_required
  IF ( .NOT. ANY(var(1:nvars) == required_vars(i)) )                          &
    CALL log_fatal(RoutineName,                                               &
                   "No value given for required variable '" //                &
                   TRIM(required_vars(i)) // "'")
END DO

!-----------------------------------------------------------------------------
! Check which variables we will be using and partition them into variables
! set to constant values and variables set from file
!-----------------------------------------------------------------------------
DO i = 1,nvars
  !-----------------------------------------------------------------------------
  ! If the variable is one of the required vars, then we will be using it
  !-----------------------------------------------------------------------------
  IF ( ANY(required_vars(1:nvars_required) == var(i)) ) THEN
    IF ( use_file(i) ) THEN
      CALL log_info(RoutineName,                                              &
                    "'" // TRIM(var(i)) // "' will be read from file")

      ! If the variable will be filled from file, register it here
      nvars_file = nvars_file + 1
      ! Since nvars_file <= i (so we will not overwrite unprocessed values)
      ! and we do not need the values from these arrays for any non-file variables
      ! from now on, we can just compress them down onto variables that are in the
      ! file
      var(nvars_file) = var(i)
      var_name(nvars_file) = var_name(i)
      tpl_name(nvars_file) = tpl_name(i)
    ELSE
      ! If the variable is being set as a constant, populate it here.
      ! First check that a value has been provided.
      IF ( ABS( const_val(i) - rmdi ) < EPSILON(1.0) )                        &
        CALL log_fatal(RoutineName,                                           &
                       "No constant value provided for variable '"            &
                       // TRIM(var(i)) // "'" )

      CALL log_info(RoutineName,                                              &
                    "'" // TRIM(var(i)) // "' will be set to a constant")

      CALL populate_var(get_var_id(var(i)), const_val = const_val(i))
    END IF
  ELSE
    ! If the variable is not a required variable, warn about not using it
    CALL log_warn(RoutineName,                                                &
                  "Provided variable '" // TRIM(var(i)) //                    &
                  "' is not required, so will be ignored")
  END IF
END DO

!-----------------------------------------------------------------------------
! Set variables from file
!-----------------------------------------------------------------------------
IF ( nvars_file > 0 ) THEN

  ! Check that a file name was provided
  IF ( LEN_TRIM(FILE) == 0 )                                                  &
     CALL log_fatal(RoutineName, "No file name provided")

  CALL log_info(RoutineName,                                                  &
              "Reading urban information from file " // TRIM(FILE))

  IF ( tpl_has_var_name(FILE) ) THEN
    ! We are using a file name template, so loop through the variables setting
    ! one from each file
    DO i = 1,nvars_file
      ! Check that a template string was provided for this variable.
      IF ( LEN_TRIM(tpl_name(i)) == 0 )                                       &
        CALL log_fatal( RoutineName,                                          &
                        "No variable name template substitution " //          &
                        "provided for " // TRIM(var(i)) )
      CALL fill_variables_from_file(                                          &
        tpl_substitute_var(FILE, tpl_name(i)),                                &
        (/ var(i) /), (/ var_name(i) /)                                       &
      )
    END DO
  ELSE
    ! We are not using a file name template, so set all variables from the same
    ! file
    CALL fill_variables_from_file(                                            &
      FILE, var(1:nvars_file), var_name(1:nvars_file)                         &
    )
  END IF
END IF


!-----------------------------------------------------------------------------
! Empirical relationships derived from correlating CEH urban fraction and
! LUCID urban geometry data for London. Obtained from collaboration with the
! University of Reading. See:
!     Bohnenstengel, S.I., Evans, S., Clark, P., Belcher, S.E. (2010);
!     Simulations of the London urban heat island, Q.J.R.Meteorol. Soc., to
!     be submitted.
! for more information

! Check for ice has been left in to be consistent with UM, but is not actually
! required here
!-----------------------------------------------------------------------------
IF ( l_urban_empirical ) THEN
  CALL log_info(RoutineName,                                                  &
                "Using empirical relationships for urban geometry: wrr")

  IF ( l_moruses )                                                            &
    CALL log_info(RoutineName,                                                &
                  "Using empirical relationships for urban geometry: hwr")

  DO l = 1, land_pts
    IF ( ainfo%frac_surft(l,urban) > 0.0 .AND.                                &
         ABS(ainfo%frac_surft(l,ice)) < EPSILON(1.0) ) THEN
      lambdap = 22.878 * ainfo%frac_surft(l,urban)**6 -                       &
                59.473 * ainfo%frac_surft(l,urban)**5 +                       &
                57.749 * ainfo%frac_surft(l,urban)**4 -                       &
                25.108 * ainfo%frac_surft(l,urban)**3 +                       &
                4.3337 * ainfo%frac_surft(l,urban)**2 +                       &
                0.1926 * ainfo%frac_surft(l,urban)    +                       &
                0.036
      urban_param%wrr_gb(l) = 1.0 - lambdap

      IF ( l_moruses ) THEN
        lambdaf = 16.412 * ainfo%frac_surft(l,urban)**6 -                     &
                  41.855 * ainfo%frac_surft(l,urban)**5 +                     &
                  40.387 * ainfo%frac_surft(l,urban)**4 -                     &
                  17.759 * ainfo%frac_surft(l,urban)**3 +                     &
                  3.2399 * ainfo%frac_surft(l,urban)**2 +                     &
                  0.0626 * ainfo%frac_surft(l,urban)    +                     &
                  0.0271
        urban_param%hwr_gb(l) = 4.0 * ATAN(1.0) / 2.0 * lambdaf               &
                                / ( 1.0 - lambdap )
      END IF
    END IF
  END DO
END IF

IF ( l_moruses ) THEN
  IF ( l_urban_empirical ) THEN
    CALL log_info(RoutineName,                                                &
                  "Using empirical relationships for urban geometry: hgt")

    DO l = 1,land_pts
      IF (ainfo%frac_surft(l,urban) > 0.0 .AND.                               &
          ABS(ainfo%frac_surft(l,ice)) < EPSILON(1.0) ) THEN
        urban_param%hgt_gb(l) = 167.409 * ainfo%frac_surft(l,urban)**5 -      &
                    337.853 * ainfo%frac_surft(l,urban)**4 +                  &
                    247.813 * ainfo%frac_surft(l,urban)**3 -                  &
                    76.3678 * ainfo%frac_surft(l,urban)**2 +                  &
                    11.4832 * ainfo%frac_surft(l,urban)  +                    &
                    4.48226
      END IF
    END DO
  END IF

  IF ( l_moruses_macdonald ) THEN
    ! Macdonald Formulation
    CALL log_info(RoutineName, "Using MacDonald formulation")

    sc_hwr(:) = 0.5 * ( urban_param%hwr_gb(:) / (2.0 * ATAN(1.0)) )
    d_h(:)    = 1.0 - urban_param%wrr_gb(:) *                                 &
                ( a**(urban_param%wrr_gb(:) - 1.0) )
    DO l = 1,land_pts
      IF ( urban_param%wrr_gb(l) > 0.0 .AND. urban_param%wrr_gb(l) < 1.0 ) THEN
        urban_param%disp_gb(l) = d_h(l) * urban_param%hgt_gb(l)
        urban_param%ztm_gb(l)  = (cdz * (1.0 - d_h(l)) *                      &
                     sc_hwr(l) * urban_param%wrr_gb(l) / kappa2)**(-0.5)
        urban_param%ztm_gb(l)  = (1.0 - d_h(l)) * EXP(-urban_param%ztm_gb(l))
        urban_param%ztm_gb(l)  = urban_param%ztm_gb(l) * urban_param%hgt_gb(l)
        urban_param%ztm_gb(l)  = MAX(urban_param%ztm_gb(l), z0m_mat)
      ELSE
        urban_param%disp_gb(l) = 0.0
        urban_param%ztm_gb(l)  = 0.0
      END IF
    END DO
  END IF
END IF

RETURN

END SUBROUTINE init_urban_props
#endif
