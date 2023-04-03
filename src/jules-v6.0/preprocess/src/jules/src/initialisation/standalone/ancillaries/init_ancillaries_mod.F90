! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE init_ancillaries_mod

USE input_mod, ONLY: fill_variables_from_file

USE logging_mod, ONLY: log_info, log_warn, log_fatal

IMPLICIT NONE

PRIVATE
PUBLIC init_ancillaries

CONTAINS

! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_ancillaries(nml_dir, crop_vars, ainfo, trif_vars, urban_param,&
                            trifctltype)

USE io_constants, ONLY: namelist_unit

USE string_utils_mod, ONLY: to_string

USE aero, ONLY: co2_mmr_mod => co2_mmr

USE dump_mod, ONLY: ancil_dump_read

USE init_flake_ancils_mod, ONLY: init_flake_ancils

USE errormessagelength_mod, ONLY: errormessagelength

USE jules_surface_mod, ONLY: l_urban2t

USE switches_urban, ONLY: l_moruses

USE um_types, ONLY: real_jlslsm

USE jules_model_environment_mod, ONLY: lsm_id, cable

USE allocate_cable_arrays_mod, ONLY: allocate_cable_progs

!TYPE definitions
USE crop_vars_mod, ONLY: crop_vars_type
USE ancil_info,    ONLY: ainfo_type
USE trif_vars_mod, ONLY: trif_vars_type
USE urban_param_mod, ONLY: urban_param_type
USE trifctl,   ONLY: trifctl_type

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in information about the model ancillaries
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Arguments
CHARACTER(LEN=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                         ! namelists

!TYPES containing field data (IN OUT)
TYPE(crop_vars_type), INTENT(IN OUT) :: crop_vars
TYPE(ainfo_type), INTENT(IN OUT) :: ainfo
TYPE(trif_vars_type), INTENT(IN OUT) :: trif_vars
TYPE(urban_param_type), INTENT(IN OUT) :: urban_param
TYPE(trifctl_type), INTENT(IN OUT) :: trifctltype

CHARACTER(LEN=errormessagelength) :: iomessage

! Work variables
INTEGER :: error  ! Error indicator

LOGICAL :: read_from_dump
REAL(KIND=real_jlslsm)    :: co2_mmr

CHARACTER(LEN=*), PARAMETER :: RoutineName='init_ancillaries'

NAMELIST  / jules_co2 / read_from_dump, co2_mmr

!-----------------------------------------------------------------------------
! Open the ancillaries namelist file
OPEN(namelist_unit, FILE=(TRIM(nml_dir) // '/' // 'ancillaries.nml'),         &
               STATUS='old', POSITION='rewind', ACTION='read', IOSTAT = error,&
               IOMSG = iomessage)
IF ( error /= 0 ) THEN
  CALL log_fatal(RoutineName,                                                 &
                 "Error opening namelist file ancillaries.nml " //            &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")
END IF

! Defer to specialist routines to process each namelist
CALL init_frac()
CALL init_vegetation_props()
CALL init_soil_props()
CALL init_top()
CALL init_pdm()
CALL init_agric(trif_vars%frac_past_gb,trifctltype)
CALL init_crop_props()
CALL init_irrig_props()
CALL init_rivers_props()
CALL init_water_resources_props()
CALL init_flake_ancils()
IF ( l_urban2t .OR. l_moruses ) CALL init_urban_props(ainfo,urban_param)
IF (lsm_id == cable) THEN
  CALL allocate_cable_progs()
  CALL init_cable_progs()
END IF

! Read the JULES_CO2 namelist
! This is so simple it doesn't really need its own subroutine
read_from_dump = .FALSE.

!Copy across the default value set up in the module.
co2_mmr = co2_mmr_mod

CALL log_info(RoutineName, "Reading JULES_CO2 namelist...")
READ(namelist_unit, NML = jules_co2, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 ) THEN
  CALL log_fatal(RoutineName,                                                 &
                 "Error reading namelist JULES_CO2 " //                       &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")
END IF

ancil_dump_read%co2 = read_from_dump

IF ( .NOT. ancil_dump_read%co2) THEN
  co2_mmr_mod = co2_mmr
ELSE
  ! We read from the dump file.
  CALL log_info(RoutineName,                                                  &
                "co2_mmr will be read from the dump file.  " //               &
                "Namelist value ignored")
END IF

! Close the namelist file
CLOSE(namelist_unit, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 ) THEN
  CALL log_fatal(RoutineName,                                                 &
                 "Error closing namelist file ancillaries.nml " //            &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")
END IF

RETURN

END SUBROUTINE init_ancillaries
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_frac()

USE data_cube_mod, ONLY: data_cube, cube_create

USE mem_brams_jules, ONLY: fracB,nxB,nyB !JWS

USE jules_surface_types_mod,       ONLY : ntype

USE io_constants, ONLY: max_file_name_len, max_sdf_name_len, namelist_unit

USE string_utils_mod, ONLY: to_string

USE jules_vegetation_mod, ONLY: l_veg_compete
  
USE dump_mod, ONLY: ancil_dump_read

USE errormessagelength_mod, ONLY: errormessagelength

USE model_interface_mod, ONLY: get_var_id, populate_var

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in the tile fractions and checks for consistency
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables

TYPE(data_cube) :: DATA


INTEGER :: error,i,j,k,l  ! Error indicator

!-----------------------------------------------------------------------------
! Definition of the jules_frac namelist
!-----------------------------------------------------------------------------
LOGICAL :: read_from_dump
CHARACTER(LEN=max_file_name_len) :: FILE
CHARACTER(LEN=max_sdf_name_len) :: frac_name
CHARACTER(LEN=errormessagelength) :: iomessage
NAMELIST  / jules_frac/ read_from_dump, FILE, frac_name


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
read_from_dump = .FALSE.
FILE=''
frac_name      = ''

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
CALL log_info("init_frac", "Reading JULES_FRAC namelist...")

!First, we read the namelist
READ(namelist_unit, NML = jules_frac, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_frac",                                                 &
                "Error reading namelist JULES_FRAC " //                      &
                "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                TRIM(iomessage) // ")")

ancil_dump_read%frac = read_from_dump

IF ( .NOT. ancil_dump_read%frac) THEN !we read from the ancil file
  !-------------------------------------------------------------------------
  ! Set frac using namelist values
  !-------------------------------------------------------------------------
  ! Frac is prognostic if competing veg is on, so must be read as an initial
  ! condition
  IF ( l_veg_compete ) THEN
    CALL log_info("init_frac",                                                &
                  "Competing vegetation is enabled - frac will be read " //   &
                  "as an initial condition")
    RETURN
  END IF

  !Check that file name was provided
IF ( LEN_TRIM(FILE) == 0 )                                                  &
  CALL log_fatal("init_frac", "No file name provided")

!JWS CALL log_info("init_frac",                                                  &
!JWS              "Reading tile fractions from file " // TRIM(FILE))


open(unit=66,file='jules.log',position='append',status='old',action='write')
write(unit=66,fmt='(A)') 'Reading tile fractions from BRAMS'
close(unit=66)

!JWS CALL fill_variables_from_file(FILE, (/ 'frac' /), (/ frac_name /))

DATA = cube_create((/nxB, nyB, ntype/))

   DATA%SHAPE(1)=nxB
   DATA%SHAPE(2)=nyB
   DATA%SHAPE(3)=ntype
   l=0
   do k=1,ntype
      do j=1,nyB
         do i=1,nxB
            l=l+1
            DATA%values(l)=fracB(i,j,k)
         enddo
      enddo
   enddo

CALL populate_var(get_var_id('frac'), DATA)

ELSE !We read from the dump file
  CALL log_info("init_frac",                                                  &
                "frac will be read from the dump file.  " //                  &
                "Namelist values ignored")

END IF !.NOT. ancil_dump_read%frac

RETURN

END SUBROUTINE init_frac
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
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_agric(frac_past_grid,trifctltype)

USE io_constants, ONLY: mdi, max_sdf_name_len, max_file_name_len,             &
                         namelist_unit

USE string_utils_mod, ONLY: to_string

USE input_mod, ONLY: input_grid => grid

USE jules_vegetation_mod, ONLY: l_trif_crop

USE dump_mod, ONLY: ancil_dump_read

USE errormessagelength_mod, ONLY: errormessagelength

USE um_types, ONLY: real_jlslsm

!TYPE Definitions
USE trifctl, ONLY: trifctl_type

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises the agricultural fraction
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Arguments
REAL(KIND=real_jlslsm), INTENT(OUT) :: frac_past_grid(:)
                               !Fraction of pasture
! TYPE Definitions
TYPE(trifctl_type), INTENT(IN OUT) :: trifctltype

! Work variables
INTEGER :: error  ! Error indicator

!-----------------------------------------------------------------------------
! Definition of the jules_agric namelist
! For data at a single point, a single fraction is specified in the namelist
! In all other cases, agricultural fraction is read from a file
!-----------------------------------------------------------------------------
LOGICAL :: read_from_dump
LOGICAL :: zero_agric           ! T - set agr. frac. at all points to 0.0
                                ! F - read agr. frac. from input
LOGICAL :: zero_past            ! T - set agr. frac. at all points to 0.0
                                ! F - read agr. frac. from input
REAL(KIND=real_jlslsm) :: frac_agr                ! Single point fraction
REAL(KIND=real_jlslsm) :: frac_past               ! Single point fraction

CHARACTER(LEN=max_file_name_len) :: FILE
                                ! The file to read fraction from
CHARACTER(LEN=max_sdf_name_len) :: agric_name
                                ! The name of the variable in the file
CHARACTER(LEN=max_file_name_len) :: file_past
                                ! The file to read fraction from
CHARACTER(LEN=max_sdf_name_len) :: past_name
                                ! The name of the variable in the file
CHARACTER(LEN=errormessagelength) :: iomessage
                                ! I/O error string
NAMELIST  / jules_agric/ read_from_dump, zero_agric, zero_past, frac_agr,     &
                       frac_past,FILE, agric_name, file_past, past_name


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
read_from_dump = .FALSE.
zero_agric = .TRUE.  ! Default is to set agricultural fraction to 0 everywhere
zero_past  = .TRUE.  ! Default is to set pasture fraction to 0 everywhere
frac_agr   = mdi     ! Initialised to missing data so we can tell if it is set
                     ! using the namelist
frac_past  = mdi     ! Initialised to missing data so we can tell if it is set
                     ! using the namelist
FILE=''
agric_name = ''
file_past  = ''
past_name = ''

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
CALL log_info("init_agric", "Reading JULES_AGRIC namelist...")

! First, we read the namelist
READ(namelist_unit, NML = jules_agric, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_agric",                                                &
                 "Error reading namelist JULES_AGRIC " //                     &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

ancil_dump_read%agric = read_from_dump

IF ( .NOT. l_trif_crop) zero_past = .TRUE.

IF ( ancil_dump_read%agric) THEN !We read from the dump file
  CALL log_info("init_agric",                                                 &
                "agric ancils will be read from the dump file.  " //          &
                "Namelist values ignored")
ELSE !we read from the ancil file
  !---------------------------------------------------------------------------
  ! Set values derived from namelist and verify for consistency
  !---------------------------------------------------------------------------
  ! If zero fraction is selected, then that is all we have to do
  IF ( zero_agric ) THEN
    CALL log_info("init_agric", "Zero agricultural fraction indicated")
    trifctltype%frac_agr_gb(:) = 0.0
    frac_past_grid(:) = 0.0
    RETURN
  END IF

  IF ( input_grid%nx * input_grid%ny <= 1 ) THEN
    !-------------------------------------------------------------------------
    ! If we are running a single point, read fraction from the namelist
    !-------------------------------------------------------------------------
    ! Check if frac_agr was set via the namelist.
    IF ( ABS(frac_agr - mdi) < EPSILON(frac_agr) )                            &
      CALL log_fatal("init_agric",                                            &
                     "For data at a single point, agricultural fraction "//   &
                     "is read from the namelist JULES_AGRIC")

    ! Now we know they have been set, copy their values into the model arrays
    CALL log_info("init_agric",                                               &
                  "Data is at a single point - reading agricultural " //      &
                  "fraction from namelist JULES_AGRIC")
    trifctltype%frac_agr_gb(1) = frac_agr
  ELSE
    !-------------------------------------------------------------------------
    ! If we have a grid, set fraction from the specified file
    !-------------------------------------------------------------------------
    CALL log_info("init_agric",                                               &
                  "Data is on a grid - reading agricultural fraction " //     &
                  "from file " // TRIM(FILE))

    ! Check that file name was provided
    IF ( LEN_TRIM(FILE) == 0 )                                                &
      CALL log_fatal("init_agric", "No file name provided")

    CALL fill_variables_from_file(FILE, (/ 'frac_agr' /), (/ agric_name /))
  END IF

  !---------------------------------------------------------------------------
  ! Check that the values seem sensible
  !---------------------------------------------------------------------------
  IF ( ANY(trifctltype%frac_agr_gb < 0.0) .OR.                                &
       ANY(trifctltype%frac_agr_gb > 1.0) )                                   &
    CALL log_fatal("init_agric",                                              &
                   "Agricultural fraction should be in range 0.0 to 1.0, "//  &
                   "given range is " //                                       &
                   TRIM(to_string(MINVAL(trifctltype%frac_agr_gb))) //        &
                   " to " //                                                  &
                   TRIM(to_string(MAXVAL(trifctltype%frac_agr_gb))))
END IF ! ancil_dump_read%agric

!---------------------------------------------------------------------------
! Read Pasture Ancil
!---------------------------------------------------------------------------

IF ( .NOT. ancil_dump_read%agric) THEN !we read from the ancil file
  !---------------------------------------------------------------------------
  ! Set values derived from namelist and verify for consistency
  !---------------------------------------------------------------------------
  ! If zero fraction is selected, then that is all we have to do
  IF ( zero_past ) THEN
    CALL log_info("init_agric", "Zero pasture fraction indicated")
    frac_past_grid(:) = 0.0
    RETURN
  END IF

  IF ( input_grid%nx * input_grid%ny <= 1 ) THEN
    !-------------------------------------------------------------------------
    ! If we are running a single point, read fraction from the namelist
    !-------------------------------------------------------------------------
    ! Check if frac_past was set via the namelist.
    IF ( ABS(frac_past - mdi) < EPSILON(frac_past) )                          &
      CALL log_fatal("init_agric",                                            &
                     "For data at a single point, pasture fraction "//        &
                     "is read from the namelist JULES_AGRIC")

    ! Now we know they have been set, copy their values into the model arrays
    CALL log_info("init_agric",                                               &
                  "Data is at a single point - reading pasture " //           &
                  "fraction from namelist JULES_AGRIC")
    frac_past_grid(1) = frac_past
  ELSE
    !-------------------------------------------------------------------------
    ! If we have a grid, set fraction from the specified file
    !-------------------------------------------------------------------------
    CALL log_info("init_agric",                                               &
                  "Data is on a grid - reading pasture fraction " //          &
                  "from file " // TRIM(file_past))

    ! Check that file name was provided
    IF ( LEN_TRIM(file_past) == 0 )                                           &
      CALL log_fatal("init_agric", "No file name provided")

    CALL fill_variables_from_file(file_past, (/ 'frac_past' /), (/ past_name /))
  END IF

  !---------------------------------------------------------------------------
  ! Check that the values seem sensible
  !---------------------------------------------------------------------------
  IF ( ANY(frac_past_grid < 0.0) .OR. ANY(frac_past_grid > 1.0) )             &
    CALL log_fatal("init_agric",                                              &
                   "Pasture fraction should be in range 0.0 to 1.0, "//       &
                   "given range is " //                                       &
                   TRIM(to_string(MINVAL(frac_past_grid))) // " to " //       &
                   TRIM(to_string(MAXVAL(frac_past_grid))))

ELSE !We read from the dump file
  CALL log_info("init_agric",                                                 &
                "agric ancils will be read from the dump file.  " //          &
                "Namelist values ignored")

END IF !.NOT. ancil_dump_read%agric

RETURN

END SUBROUTINE init_agric

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
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

SUBROUTINE init_irrig_props()

USE io_constants, ONLY: max_file_name_len, max_sdf_name_len, namelist_unit

USE missing_data_mod, ONLY:                                                   &
!  imported scalar parameters
     imdi, rmdi

USE model_interface_mod, ONLY: populate_var, get_var_id

USE string_utils_mod, ONLY: to_string
                         
USE jules_irrig_mod, ONLY: l_irrig_dmd,  irrtiles, frac_irrig_all_tiles,      &
                           nirrtile, irrigtiles,                              &
                           set_irrfrac_on_irrtiles

USE dump_mod, ONLY: ancil_dump_read

USE jules_surface_types_mod, ONLY: npft

USE um_types, ONLY: real_jlslsm

USE jules_print_mgr, ONLY: jules_print

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!    Initialises irrigation parameters and properties
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
INTEGER :: error  ! Error indicator

INTEGER :: i  ! Loop counter

CHARACTER(LEN=50000) :: lineBuffer

!-----------------------------------------------------------------------------
! Definition of the jules_irrig_props namelist
!-----------------------------------------------------------------------------
LOGICAL :: read_from_dump
CHARACTER(LEN=max_file_name_len) :: irrig_frac_file
                      ! The name of the file to use for irrigation fraction
                      ! if read_file = .false.
LOGICAL :: read_file
                      !   T - read from file
                      !   F - use a constant value for all points
CHARACTER(LEN=max_sdf_name_len) :: var_name
                      ! The name of the variable in the file
REAL(KIND=real_jlslsm) :: const_frac_irr = rmdi
                      ! The constant value to use for each variable if
                      ! read_file = F and frac_irrig_all_tiles = T
                      ! for that variable
REAL :: const_irrfrac_irrtiles = rmdi
                      ! The constant value to use for each variable if
                      ! read_file = F and set_irrig_on_irrtiles = T
                      ! for that variable

NAMELIST  / jules_irrig_props/ read_from_dump, irrig_frac_file, read_file, var_name, const_frac_irr, &
                       const_irrfrac_irrtiles

!-----------------------------------------------------------------------------

! Nothing to do if irrigation demand model is not on
IF ( .NOT. l_irrig_dmd ) RETURN

!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
read_from_dump = .FALSE.
read_file = .TRUE.       ! Default is to read every variable from file
irrig_frac_file =''      ! Empty file name
var_name  = ''           ! Empty variable name.

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
CALL log_info("init_irrig_props", "Reading JULES_IRRIG_PROPS namelist...")

READ(namelist_unit, NML = jules_irrig_props, IOSTAT = error)

IF ( error /= 0 )                                                             &
  CALL log_fatal("init_irrig_props",                                          &
                "Error reading namelist JULES_IRRIG_PROPS " //                &
                "(IOSTAT=" // TRIM(to_string(error)) // ")")

! Print the contents of the namelist_unit
CALL jules_print('init_irrig_props', 'Contents of namelist init_irrig_props')

WRITE(lineBuffer,*) ' read_from_dump = ',read_from_dump
CALL jules_print('init_irrig_props',lineBuffer)

WRITE(lineBuffer,*) ' irrig_frac_file = ',irrig_frac_file
CALL jules_print('init_irrig_props',lineBuffer)

WRITE(lineBuffer,*) ' read_file = ',read_file
CALL jules_print('init_irrig_props',lineBuffer)

WRITE(lineBuffer,*) ' var_name = ',var_name
CALL jules_print('init_irrig_props',lineBuffer)

WRITE(lineBuffer,*) ' const_frac_irr = ',const_frac_irr
CALL jules_print('init_irrig_props',lineBuffer)

WRITE(lineBuffer,*) ' const_irrfrac_irrtiles = ',const_irrfrac_irrtiles
CALL jules_print('init_irrig_props',lineBuffer)

ancil_dump_read%irrig = read_from_dump

IF ( .NOT. ancil_dump_read%irrig) THEN !we read from the ancil file

  !---------------------------------------------------------------------------
  ! Determine if variable to be read from file
  !---------------------------------------------------------------------------
  IF ( read_file ) THEN
    CALL log_info("init_irrig_props", "'frac_irrig' will be read from file")

    ! If the variable will be filled from file, register it here
    IF ( LEN_TRIM(irrig_frac_file) == 0 )                                     &
      CALL log_fatal("init_irrig_props", "No file name provided")

    CALL fill_variables_from_file(irrig_frac_file, (/ 'frac_irrig' /), (/ var_name /))
  ELSE
    ! If the variable is being set as a constant, populate it here.
    IF ( set_irrfrac_on_irrtiles ) THEN
       ! The irrigation fraction is specified for irrigated tiles
       ! using irrfrac_irrtiles
       ! First check that a value has been provided.
      IF ( ABS( const_irrfrac_irrtiles - rmdi ) < EPSILON(1.0) )              &
        CALL log_fatal("init_irrig_props",                                    &
                  "No constant value provided for irrfrac_irrtiles.")
      CALL log_info("init_irrig_props",                                       &
               "'const_irrfrac_irrtiles' will be set to a constant = " //     &
                        to_string(const_irrfrac_irrtiles))
      CALL populate_var(get_var_id('irrfrac_irrtiles'),                       &
                                        const_val = const_irrfrac_irrtiles)
    ELSE
       ! The irrigated fraction is set as const_frac_irr from namelist.
       ! First check that a value has been provided.
      IF ( ABS( const_frac_irr - rmdi ) < EPSILON(1.0) )                      &
        CALL log_fatal("init_irrig_props",                                    &
                     "No constant value provided for frac_irrig.")

      CALL log_info("init_irrig_props",                                       &
                  "'frac_irrig' will be set to a constant = " //              &
                  to_string(const_frac_irr))

      CALL populate_var(get_var_id('frac_irrig'), const_val = const_frac_irr)
    END IF
  END IF

  !---------------------------------------------------------------------------
  ! Process pft names to be assigned irr fraction (if not all)
  !---------------------------------------------------------------------------
  IF ( .NOT. frac_irrig_all_tiles ) THEN
    IF ( ( nirrtile > npft ) .OR. ( nirrtile < 1 ) ) THEN
      CALL log_info("init_irrig_props",                                       &
                    "nirrtile should be greater than 1 and less than or" //   &
                    "equal to the number of pfts.")
    END IF
    ! Copy pft indices from namelist into allocated array.
    ! First check that values were given.
    IF ( ANY( irrigtiles(1:nirrtile) == imdi ) )                              &
      CALL log_fatal("init_irrig_props",                                      &
                     "Insufficient values provided for irrigtiles.")
                     
    ! irrtiles is now allocated with SIZE(npft) instead of SIZE(npft_max).
    ! It should really be allocated with SIZE(nirrtile) though, this will 
    ! be changed through ticket jules:#1065
    irrtiles(1:nirrtile) = irrigtiles(1:nirrtile)
  END IF

ELSE !We read from the dump file
  CALL log_info("init_frac",                                                  &
                "frac will be read from the dump file.  " //                  &
                "Namelist values ignored")

END IF !.NOT. ancil_dump_read%frac

END SUBROUTINE init_irrig_props
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

SUBROUTINE init_rivers_props()

USE mpi, ONLY: mpi_comm_world

USE io_constants, ONLY: max_sdf_name_len, max_file_name_len, namelist_unit

USE grid_utils_mod, ONLY: grid_create, grid_info

USE input_mod, ONLY:                                                          &
    od_grid => grid,                                                          &
    dummy_grid => grid,                                                       &
    fill_variables_from_file, use_subgrid

USE missing_data_mod, ONLY: imdi, rmdi

USE ancil_info, ONLY: land_pts

USE string_utils_mod, ONLY: to_string

USE templating_mod, ONLY: tpl_has_var_name, tpl_substitute_var

USE model_interface_mod, ONLY: identifier_len, populate_var, get_var_id

USE model_grid_mod, ONLY: latitude, longitude, global_land_pts, model_grid

USE conversions_mod, ONLY: pi_over_180

USE planet_constants_mod, ONLY: planet_radius  ! the Earth's radius (m)

USE jules_rivers_mod, ONLY:                                                   &
     i_river_vn, rivers_rfm, rivers_trip,                                     &
     flow_dir_delta, np_rivers, nx_rivers, ny_rivers,                         &
     rivers_dlat, rivers_dlon, rivers_lat1, rivers_lon1, nx_grid, ny_grid,    &
     nseqmax, reg_lon1, reg_lat1, reg_dlon, reg_dlat, rivers_dx,              &
     rivers_reglatlon, rivers_regrid, l_rivers,                               &
     rivers_index_rp, rivers_next_rp, rivers_seq_rp, rivers_dir_rp,           &
     rivers_sto_rp, rivers_dra_rp, rivers_lat_rp, rivers_lon_rp,              &
     rivers_seq, rivers_dir, rivers_dra,                                      &
     rivers_lat2d, rivers_lon2d, rivers_xgrid, rivers_ygrid,                  &
     rfm_flowobs1_rp, rfm_surfstore_rp, rfm_substore_rp, rfm_flowin_rp,       &
     rfm_bflowin_rp, rfm_rivflow_rp, rfm_baseflow_rp,                         &
     rfm_iarea_rp, rfm_land_rp, rivers_boxareas_rp,                           &
     rivers_sto_per_m2_on_landpts, rivers_adj_on_landpts,                     &
     il_river_grid, ir_land_grid, a_thresh

USE rivers_utils, ONLY: rivers_earth_area

USE rivers_regrid_mod, ONLY: rivers_remap_match, rivers_remap_unmatch,        &
                             rivers_get_xy_pos

USE overbank_inundation_mod, ONLY:                                            &
     l_riv_overbank, l_riv_hypsometry, use_rosgen,                            &
     logn_mean, logn_stdev, logn_mean_rp, logn_stdev_rp,                      &
     qbf, wbf, dbf,                                                           &
     riv_a, riv_b, riv_c, riv_f, coef_b, exp_c

USE jules_irrig_mod, ONLY: l_irrig_limit

USE jules_print_mgr, ONLY:                                                    &
  jules_message,                                                              &
  jules_print

USE string_utils_mod, ONLY: to_string
USE logging_mod, ONLY: log_info

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!    Reads in details of the river routing grid and river properties,
!    and allocates various river variables.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Work variables

INTEGER :: error, error_sum  ! Error indicators

INTEGER :: ntasks ! Parallel mode indicator

INTEGER, PARAMETER :: max_rivers_vars = 6 ! The maximum possible number of
                                          ! routing variables

INTEGER :: nvars_required      ! The number of routing variables that are
                               ! required in this configuration
INTEGER :: nvars_optional      ! The number of optional routing variables
                               ! in this configuration
CHARACTER(LEN=identifier_len) :: required_vars(max_rivers_vars)
                               ! The variable identifiers of the required
                               ! variables
CHARACTER(LEN=identifier_len) :: optional_vars(max_rivers_vars)
                               ! The variable identifiers of any optional
                               ! variables

INTEGER :: nvars_file          ! The number of variables that will be set
                               ! from the given file (template?)

! Variables passed to fill_variables_from_file

CHARACTER(LEN=identifier_len) :: file_var(max_rivers_vars)
                      ! The variable identifiers of the variables to set
                      ! from file
CHARACTER(LEN=max_sdf_name_len) :: file_var_name(max_rivers_vars)
                      ! The name of each variable in the file
CHARACTER(LEN=max_sdf_name_len) :: file_tpl_name(max_rivers_vars)
                      ! The name to substitute in a template for each
                      ! variable

CHARACTER(LEN=max_file_name_len) :: FILE
                      ! The name of the file (or variable name template) to
                      ! use for variables that need to be filled from file

INTEGER :: nvars      ! The number of variables in this section
CHARACTER(LEN=identifier_len) :: var(max_rivers_vars)
                      ! The variable identifiers of the variables
LOGICAL :: use_file(max_rivers_vars)
                      !   T - the variable uses the file
                      !   F - the variable is set using a constant value
CHARACTER(LEN=max_sdf_name_len) :: var_name(max_rivers_vars)
                      ! The name of each variable in the file
CHARACTER(LEN=max_sdf_name_len) :: tpl_name(max_rivers_vars)
                      ! The name to substitute in a template for each
                      ! variable
REAL(KIND=real_jlslsm) :: const_val(max_rivers_vars)
                      ! The constant value to use for each variable if
                      ! use_file = F for that variable

INTEGER :: i,ip,ix,iy,j  ! Index variables
INTEGER :: inext, jnext
INTEGER :: i1, i2, ilat, ilon, step
REAL(KIND=real_jlslsm) :: reg_lat2, reg_lon2, min_lon, dlat, dlon, tiny_0
REAL(KIND=real_jlslsm) :: contribarea
                      ! (upstream contributing drainag area (km^2)

!Conversion from m2 to km2
REAL(KIND=real_jlslsm), PARAMETER :: m2tokm2 = 1.0e-6

LOGICAL :: use_sub_local
LOGICAL :: set_dra

! Remapping full grid to vector points
INTEGER, ALLOCATABLE :: mapfr(:,:)                   ! map full to river

! The names and sizes of the x and y dimensions

CHARACTER(LEN=max_sdf_name_len) :: x_dim_name, y_dim_name
INTEGER :: nx, ny

TYPE(grid_info) :: local_grid, rivers_grid

NAMELIST  / jules_rivers_props/                                               &
   x_dim_name, y_dim_name, nx, ny, FILE,                                      &
   nvars, var, use_file, var_name, tpl_name, const_val,                       &
   rivers_regrid, rivers_reglatlon,                                           &
   nx_grid, ny_grid, reg_lat1, reg_lon1, reg_dlon, reg_dlat, rivers_dx

!-----------------------------------------------------------------------------
! Notes on river routing grid assumptions (see also jules_rivers_mod):
!
!        The JULES river routing implementation assumes that river routing
!        variables are defined on a river grid with cell (1,1) in the lower
!        left corner (i.e. W-E, S-N). Ancillary information may be read in
!        on alternative grids, but this routine will translate these points
!        to the anticipated ordering.
!
!        The JULES river routing algorithms also require information on the
!        dimensions of the full model grid, which may not be typically defined
!        for some standalone applications where only land_pts are specified.
!        Information on the full model grid is therefore also read in and
!        checked here.
!
!-----------------------------------------------------------------------------

! If rivers are not on, there is nothing to do
IF ( .NOT. l_rivers ) RETURN

!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------

nvars_required = 0
nvars_optional = 0
nvars_file     = 0
nvars          = 0
use_file(:)    = .TRUE. ! Default is for every variable to be read from file
FILE=''     ! Empty file name.

IF ( global_land_pts <= 1 ) THEN
  l_rivers = .FALSE.
  CALL log_warn("init_rivers_props",                                          &
                "River routing not appropriate for single point runs. " //    &
                "River routing disabled.")

ELSE

  !---------------------------------------------------------------------------
  ! Read river routing properties namelist
  !---------------------------------------------------------------------------
  CALL log_info("init_rivers_props", "Reading JULES_RIVERS_PROPS namelist...")

  READ(namelist_unit, NML = jules_rivers_props, IOSTAT = error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_rivers_props",                                       &
                   "Error reading namelist JULES_RIVERS_PROPS " //            &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  !---------------------------------------------------------------------------
  ! Define regular routing grid lat/lon
  !---------------------------------------------------------------------------

  use_sub_local = use_subgrid
  use_subgrid   = .FALSE.

  ! Temporarily copy saved grid (set for full model grid) to a local variable
  ! before overwriting to define the river routing grid
  local_grid = grid_create(dummy_grid%is_1d,dummy_grid%dim_name,              &
                    dummy_grid%nx,dummy_grid%x_name,dummy_grid%nx,            &
                    dummy_grid%y_name,dummy_grid%ny)

  ! Read latitude dimension from file
  od_grid = grid_create( grid_is_1d = .TRUE.,                                 &
                         dim_name = TRIM(y_dim_name), npoints = ny )
  ALLOCATE(rivers_ygrid(ny), stat = error)
  CALL fill_variables_from_file(FILE,                                         &
                                (/ 'rivers_ygrid' /),(/ TRIM(y_dim_name) /))

  ! Read longitude dimension from file
  od_grid = grid_create( grid_is_1d = .TRUE.,                                 &
                         dim_name = TRIM(x_dim_name), npoints = nx)
  ALLOCATE(rivers_xgrid(nx), stat = error)
  CALL fill_variables_from_file(FILE,                                         &
                                (/ 'rivers_xgrid' /),(/ TRIM(x_dim_name) /))
 
  !---------------------------------------------------------------------------
  ! Build the routing grid object from the namelist values
  !---------------------------------------------------------------------------
  dummy_grid = grid_create( .FALSE., "", 0, x_dim_name, nx, y_dim_name, ny)
  rivers_grid = grid_create( .FALSE., "", 0, x_dim_name, nx, y_dim_name, ny)

  nx_rivers = rivers_grid%nx
  ny_rivers = rivers_grid%ny

  ! Set all river routing arrays to be S to N irrespective of how latitudes
  ! are stored in the ancillary files
  IF ( rivers_ygrid(1) > MINVAL(rivers_ygrid) ) THEN
    i1   = ny_rivers
    i2   = 1
    step = -1
  ELSE
    i1   = 1
    i2   = ny_rivers
    step = 1
  END IF

  !---------------------------------------------------------------------------
  ! Set up routing properties using namelist values
  !---------------------------------------------------------------------------

  ! Set up the required and optional routing variables

  SELECT CASE ( i_river_vn )

  CASE ( rivers_trip )

    nvars_required = 2
    required_vars(1:nvars_required) = (/ 'direction', 'sequence ' /)
    nvars_optional = 2
    optional_vars(1:nvars_optional) = (/ 'latitude_2d ', 'longitude_2d' /)

  CASE ( rivers_rfm )

    IF ( rivers_reglatlon ) THEN
      nvars_required = 1
      required_vars(1:nvars_required) = (/ 'direction' /)
      nvars_optional = 4
      optional_vars(1:nvars_optional) =                                       &
        (/ 'latitude_2d ', 'longitude_2d',                                    &
           'area        ', 'sequence    ' /)
    ELSE
      nvars_required = 3
      required_vars(1:nvars_required) = (/ 'direction   ',                    &
                                  'latitude_2d ', 'longitude_2d' /)
      nvars_optional = 2
      optional_vars(1:nvars_optional) =                                       &
        (/ 'area        ', 'sequence    ' /)
    END IF

  CASE DEFAULT
    CALL log_fatal("init_rivers_props",                                       &
                   "Do not recognise i_river_vn = '" //                       &
                   TRIM(to_string(i_river_vn)) // "'")
  END SELECT

  !---------------------------------------------------------------------------
  ! Check that all the required variables are there
  !---------------------------------------------------------------------------
  DO i = 1,nvars_required
    IF ( .NOT. ANY(var(1:nvars) == required_vars(i)) )                        &
      CALL log_fatal("init_rivers_props",                                     &
                     "No value given for required variable '" //              &
                     TRIM(required_vars(i)) // "'")
  END DO

  !---------------------------------------------------------------------------
  ! Allocate routing-specific arrays from ancillary information
  !---------------------------------------------------------------------------

  IF ( nx_rivers > 0 .AND. ny_rivers > 0 ) THEN

    ALLOCATE( rivers_dir(nx_rivers, ny_rivers),      stat = error )
    error_sum = error
    ALLOCATE( rivers_seq(nx_rivers, ny_rivers),      stat = error )
    error_sum = error_sum + error
    ALLOCATE( rivers_dra(nx_rivers, ny_rivers),      stat = error )
    error_sum = error_sum + error
    ALLOCATE( rivers_lat2d(nx_rivers, ny_rivers),    stat = error )
    error_sum = error_sum + error
    ALLOCATE( rivers_lon2d(nx_rivers, ny_rivers),    stat = error )
    error_sum = error_sum + error
    ALLOCATE( mapfr(nx_rivers, ny_rivers),           stat = error )
    error_sum = error_sum + error
    IF ( error_sum /= 0 ) THEN
      CALL log_fatal("init_rivers_props",                                     &
                     "Error allocating for rivers arrays.")
    END IF

    ! Initialise to impossible values.

    rivers_dir(:,:)      = rmdi
    rivers_seq(:,:)      = rmdi
    rivers_dra(:,:)      = rmdi
    rivers_lat2d(:,:)    = rmdi
    rivers_lon2d(:,:)    = rmdi
    mapfr(:,:)           = 0

  ELSE

    CALL log_fatal("init_rivers_props",                                       &
                   "Invalid routing grid dimensions provided (nx_rivers = "// &
                    TRIM(to_string(nx_rivers)) // ", ny_rivers = " //         &
                    TRIM(to_string(ny_rivers)) // "). Check inputs. ")

  END IF

  !---------------------------------------------------------------------------
  ! Check which variables we will be using and partition them into variables
  ! set to constant values and variables set from file
  !---------------------------------------------------------------------------
  DO i = 1,nvars

    !-------------------------------------------------------------------------
    ! If the variable is one of the required vars, then we will be using it
    !-------------------------------------------------------------------------

    IF ( ANY(required_vars(1:nvars_required) == TRIM(var(i)))                 &
       .OR. (nvars_optional > 0 .AND.                                         &
       ANY(optional_vars(1:nvars_optional) == TRIM(var(i)))) ) THEN
      IF ( use_file(i) ) THEN
        CALL log_info("init_rivers_props",                                    &
                      "'" // TRIM(var(i)) // "' will be read from file")

        ! If the variable will be filled from file, register it here
        nvars_file = nvars_file + 1
        file_var(nvars_file) = TRIM(var(i))
        file_var_name(nvars_file) = var_name(i)
        file_tpl_name(nvars_file) = TRIM(tpl_name(i))

      ELSE

        ! If the variable is being set as a constant, just populate it here
        CALL log_info("init_rivers_props",                                    &
                      "'" // TRIM(var(i)) // "' will be set to a " //         &
                      "constant = " // to_string(const_val(i)))

        CALL populate_var(get_var_id(var(i)), const_val = const_val(i))
      END IF

    ELSE

      ! If the variable is not a required variable, warn about not using it
      CALL log_warn("init_rivers_props",                                      &
                    "Provided variable '" // TRIM(var(i)) //                  &
                    "' is not required, so will be ignored")
    END IF

  END DO ! nvars

  !---------------------------------------------------------------------------
  ! Set variables from file
  !---------------------------------------------------------------------------

  IF ( nvars_file > 0 ) THEN

    ! Check that a file name was provided
    IF ( LEN_TRIM(FILE) == 0 ) THEN
      CALL log_fatal("init_rivers_props", "No file name provided")
    END IF

    CALL log_info("init_rivers_props",                                        &
                  "Reading routing information from file " // TRIM(FILE))

    IF ( tpl_has_var_name(FILE) ) THEN
      ! We are using a file name template, so loop through the variables
      ! setting one from each file.
      DO i = 1,nvars_file
        CALL fill_variables_from_file(                                        &
          tpl_substitute_var(FILE, file_tpl_name(i)),                         &
          (/ file_var(i) /), (/ file_var_name(i) /)                           &
        )
      END DO

    ELSE
      ! We are not using a file name template, so set all variables from the
      ! same file.
      CALL fill_variables_from_file(                                          &
        FILE,file_var(1:nvars_file), file_var_name(1:nvars_file))

    END IF
  END IF

  !---------------------------------------------------------------------------
  ! Setup regular river routing grid (assuming same as model grid if not set
  ! from namelist or model grid information)
  !---------------------------------------------------------------------------

  ! Wrap longitude inputs -180 to +180
  DO ix = 1,nx
    IF ( rivers_xgrid(ix) > 180.0 ) THEN
      rivers_xgrid(ix) = rivers_xgrid(ix) - 360.0
    END IF
  END DO

  ! Set regular Lat/Lon grid variables (S to N)
  rivers_lat1 = MINVAL( rivers_ygrid )
  rivers_lon1 = MINVAL( rivers_xgrid )
  rivers_dlat = ABS( rivers_ygrid(2) - rivers_ygrid(1) )
  rivers_dlon = ABS( rivers_xgrid(2) - rivers_xgrid(1) )

  CALL mpi_comm_size(mpi_comm_world,ntasks,error)
  IF ( ntasks == 1 ) THEN
    IF ( reg_dlat > -900.0 .OR. reg_dlon > -900.0 .OR.                        &
        ny_grid > 0 .OR. nx_grid > 0 .OR.                                     &
        reg_lat1 > -900.0 .OR. reg_lon1  > -900.0 ) THEN

      CALL log_info("init_rivers_props",                                      &
                   "Running in serial mode and with grid_is_1d==.FALSE. " //  &
                   "means that the values; nx_grid, ny_grid, reg_dlat, "  //  &
                   "reg_dlon, reg_lat1 and reg_lon1 will be ignored "     //  &
                   "in the init_rivers_props namelist.")
    END IF
  END IF

  IF ( reg_lon1 <= -900.0 .OR. reg_lat1 <= -900.0 ) THEN

    ! Infer model grid settings for 1-dimensional grid input
    IF ( model_grid%nx == 1 .OR. .NOT. rivers_reglatlon ) THEN

      CALL log_info("init_rivers_props",                                      &
              "No regular model grid set - assuming same as routing grid")
      nx_grid  = nx_rivers
      ny_grid  = ny_rivers
      reg_lat1 = rivers_lat1
      reg_lon1 = rivers_lon1
      reg_dlon = rivers_dlon
      reg_dlat = rivers_dlat

    ELSE

      CALL log_info("init_rivers_props",                                      &
               "Using 2D input model grid to define grid parameters")
      nx_grid  = SIZE( latitude,1 )
      ny_grid  = SIZE( latitude,2 )
      reg_lat1 = MINVAL( latitude )
      reg_lon1 = MINVAL( longitude )
      reg_dlon = ABS( longitude(2,2) - longitude(1,1) )
      reg_dlat = ABS( latitude(2,2) - latitude(1,1) )

    END IF
  ELSE
    IF ( reg_dlat <= -900.0 .OR. reg_dlon <= -900.0 .OR.                      &
         ny_grid <= 0 .OR. nx_grid <= 0 ) THEN
      CALL log_info("init_rivers_props",                                      &
                    "If reg_lon1 and  reg_lon1 are specified in namelist," // &
                    "specify reg_dlat, reg_dlon, ny_grid, nx_grid too." )
    END IF
  END IF

  CALL log_info("init_rivers_props","Setting regular routing grid " //        &
            " Minimum " // TRIM(x_dim_name) // ": " //                        &
              TRIM(to_string(MINVAL(rivers_xgrid)))  //                       &
            " Minimum " // TRIM(y_dim_name) // ": " //                        &
              TRIM(to_string(MINVAL(rivers_ygrid))) //                        &
            " Maximum " // TRIM(x_dim_name) // ": " //                        &
              TRIM(to_string(MAXVAL(rivers_xgrid))) //                        &
            " Maximum " // TRIM(y_dim_name) // ": " //                        &
              TRIM(to_string(MAXVAL(rivers_ygrid))))

  !---------------------------------------------------------------------------
  ! Process the namelist values and set derived variables
  !---------------------------------------------------------------------------

  ! Consider routing points on valid domain only - otherwise unset rivers_dir
  !---------------------------------------------------------------------------

  reg_lat2 = reg_lat1 + ny_grid * reg_dlat
  reg_lon2 = reg_lon1 + nx_grid * reg_dlon

  DO iy = 1,ny_rivers
    DO ix = 1,nx_rivers
      IF ( (rivers_xgrid(ix) < reg_lon1) .OR.                                 &
           (rivers_ygrid(iy) < reg_lat1) .OR.                                 &
           (rivers_xgrid(ix) > reg_lon2) .OR.                                 &
           (rivers_ygrid(iy) > reg_lat2) ) THEN
        rivers_dir(ix,iy) = -1
      END IF
    END DO
  END DO

  ! Detect number of river and discharge points in river grid
  !--------------------------------------------------------------------------- 

  np_rivers = COUNT( rivers_dir >= 1 .AND. rivers_dir <=10 )

  CALL log_info("init_rivers_props",                                          &
                "River routing points = " // TRIM(to_string(np_rivers)))
  CALL log_info("init_rivers_props",                                          &
                "Global land pts = " // TRIM(to_string(global_land_pts)))

  IF ( np_rivers <= 0 ) THEN
    CALL log_fatal("init_rivers_props",                                       &
                   "Invalid number of valid routing grid points " //          &
                   "(np_rivers = " // TRIM(to_string(np_rivers)) //           &
                   "). Check inputs. ")
  END IF

  ! Allocate land point arrays and initialise
  !---------------------------------------------------------------------------

  ALLOCATE( ir_land_grid(global_land_pts), stat = error )
  error_sum = error
  IF ( error_sum /= 0 ) THEN
    CALL log_fatal("init_rivers_props",                                       &
                   "Error allocating overbank variables.")
  END IF


  ! If required, read overbank inundation ancillary information from file
  ! before translating to the river routing grid
  IF ( l_riv_overbank ) CALL init_overbank_props()

  !---------------------------------------------------------------------------
  ! Allocate routing point arrays
  !---------------------------------------------------------------------------

  ALLOCATE(rivers_index_rp(np_rivers),    stat = error)
  error_sum = error
  ALLOCATE(rivers_next_rp(np_rivers),     stat = error)
  error_sum = error_sum + error
  ALLOCATE(rivers_seq_rp(np_rivers),      stat = error)
  error_sum = error_sum + error
  ALLOCATE(rivers_sto_rp(np_rivers),      stat = error)
  error_sum = error_sum + error
  ALLOCATE(rivers_dir_rp(np_rivers),      stat = error)
  error_sum = error_sum + error
  ALLOCATE(rivers_dra_rp(np_rivers),      stat = error)
  error_sum = error_sum + error
  ALLOCATE(rivers_lat_rp(np_rivers),      stat = error)
  error_sum = error_sum + error
  ALLOCATE(rivers_lon_rp(np_rivers),      stat = error)
  error_sum = error_sum + error
  ALLOCATE(rivers_boxareas_rp(np_rivers), stat = error)
  error_sum = error_sum + error
  ALLOCATE(il_river_grid(np_rivers),      stat = error)
  error_sum = error_sum + error
  ALLOCATE(rfm_rivflow_rp(np_rivers),      stat = error)
  error_sum = error_sum + error

  ! Initialise array values

  rivers_index_rp(:) = imdi
  rivers_next_rp(:)  = 0
  rivers_seq_rp(:)   = rmdi
  rivers_dir_rp(:)   = rmdi
  rivers_sto_rp(:)   = 0.0
  il_river_grid(:)   = 0
  rfm_rivflow_rp(:)  = 0.0

  IF ( i_river_vn == rivers_rfm ) THEN
    ALLOCATE( rfm_flowobs1_rp(np_rivers),  stat = error )
    error_sum = error_sum + error
    ALLOCATE( rfm_surfstore_rp(np_rivers), stat = error )
    error_sum = error_sum + error
    ALLOCATE( rfm_substore_rp(np_rivers),  stat = error )
    error_sum = error_sum + error
    ALLOCATE( rfm_flowin_rp(np_rivers),    stat = error )
    error_sum = error_sum + error
    ALLOCATE( rfm_bflowin_rp(np_rivers),   stat = error )
    error_sum = error_sum + error
    ALLOCATE( rfm_iarea_rp(np_rivers),     stat = error )
    error_sum = error_sum + error
    ALLOCATE( rfm_land_rp(np_rivers),      stat = error )
    error_sum = error_sum + error
    ALLOCATE( rfm_baseflow_rp(np_rivers),  stat = error )
    error_sum = error_sum + error
  END IF
  IF ( error_sum /= 0 ) THEN
    CALL log_fatal("init_rivers_props",                                       &
                   "Error allocating for routing point arrays.")
  END IF

  ! Initialise array values
  IF ( i_river_vn == rivers_rfm ) THEN
    rfm_land_rp(:)      = 2       ! Initialise all points to land
    rfm_iarea_rp(:)     = imdi

    rfm_flowobs1_rp(:)  = 0.0
    rfm_surfstore_rp(:) = 0.0
    rfm_substore_rp(:)  = 0.0
    rfm_flowin_rp(:)    = 0.0
    rfm_bflowin_rp(:)   = 0.0
    rfm_baseflow_rp(:)  = 0.0

  END IF

  inext   = imdi
  jnext   = imdi
  nseqmax = imdi
  ip      = 0

  ! Initialise iarea from drainage?
  set_dra = .FALSE.
  IF ( MAXVAL(rivers_dra) > 0.0 ) set_dra = .TRUE.

  min_lon = MINVAL(rivers_xgrid)

  ! Set up routing point arrays, and correct 2d lat/lon grid (S to N)

  DO ix = 1,nx_rivers

    ilat = 0
    IF ( rivers_reglatlon ) THEN
      ilon = NINT( (rivers_xgrid(ix) - min_lon) / rivers_dlon ) + 1
      IF ( ilon > nx_rivers ) THEN
        ilon = ilon - nx_rivers
      END IF
    ELSE
      ilon = ix
    END IF

    DO iy = i1,i2,step

      ilat = ilat + 1

      IF ( rivers_reglatlon ) THEN
        ! Assume 1d x and y grid variables define lat / lon
        rivers_lat2d(ilon,ilat) = rivers_ygrid(iy)
        rivers_lon2d(ilon,ilat) = rivers_xgrid(ix)
      END IF

      ! Set appropriate routing flags etc
      IF ( rivers_dir(ix,iy) >= 1 .AND. rivers_dir(ix,iy) <=10 ) THEN

        ip = ip + 1
        rivers_index_rp(ip) = (ilat-1) * nx_rivers + ilon

        ! Update river vector information from 2D input ancillaries
        rivers_lon_rp(ip) = rivers_lon2d(ilon,ilat)
        rivers_lat_rp(ip) = rivers_lat2d(ilon,ilat)
        rivers_dir_rp(ip) = rivers_dir(ix,iy)
        rivers_dra_rp(ip) = rivers_dra(ix,iy)
        rivers_seq_rp(ip) = rivers_seq(ix,iy)
        IF ( NINT(rivers_seq_rp(ip)) > nseqmax ) THEN
          nseqmax = NINT(rivers_seq_rp(ip))
        END IF

        ! Initialise drainage area from iarea ancillary
        IF ( i_river_vn == rivers_rfm ) THEN
          IF ( set_dra ) THEN
            rfm_iarea_rp(ip) = rivers_dra(ix,iy)
          ELSE
            ! if not present, set all to river points
            rfm_iarea_rp(ip) = a_thresh + 1
          END IF
        END IF

        ! Define mapping between ip river vector and 2D x,y routing grid
        mapfr(ilon,ilat) = ip

        ! Define overbank inundation ancillary variables, if required
        IF ( l_riv_overbank .AND. l_riv_hypsometry ) THEN
          logn_mean_rp(ip) = logn_mean(ix,iy)
          logn_stdev_rp(ip) = logn_stdev(ix,iy)
        END IF

      END IF  !  rivers_dir

    END DO     ! iy river rows
  END DO    ! ix river columns

  !----------------------------------------------------------------------------
  ! Set next variables from flow directions
  !----------------------------------------------------------------------------

  DO ip = 1,np_rivers

    IF ( rivers_dir_rp(ip) >= 1 .AND. rivers_dir_rp(ip) <= 10 ) THEN
      CALL rivers_get_xy_pos(rivers_index_rp(ip),nx_rivers,ny_rivers,ix,iy)
      inext = ix + flow_dir_delta( NINT(rivers_dir_rp(ip)),1 )
      jnext = iy + flow_dir_delta( NINT(rivers_dir_rp(ip)),2 )

      ! cyclic bcs
      IF ( rivers_reglatlon ) THEN
        IF ( inext > nx_rivers ) inext = 1
        IF ( inext < 1 ) inext = inext + nx_rivers
      END IF

      IF ( mapfr(inext, jnext) > 0 ) THEN
        rivers_next_rp(ip) = mapfr(inext,jnext)
      END IF

      ! HL: Note this line is currently included for trip only to preserve
      !     rose stem testing between 'xy' and 'vector' RFM approaches
      IF ( i_river_vn == rivers_trip ) THEN
        ! If the downstream point is the same as the current point, set to
        ! zero.
        IF ( rivers_next_rp(ip) == ip ) THEN
          rivers_next_rp(ip) = 0
        END IF
      END IF

    END IF  !  rivers_dir_rp
  END DO

  DEALLOCATE(mapfr)

  !---------------------------------------------------------------------------
  ! Check routing and model grid settings are appropriate
  !---------------------------------------------------------------------------
  
  ! Check that full routing latitude/longitude grid within range
  IF ( ANY(rivers_lat2d < -90.0) .OR. ANY(rivers_lat2d > 90.0) )              &
    CALL log_fatal("init_rivers_props",                                       &
                   "Latitude is out of range - allowed range is " //          &
                   "-90.0 to 90.0, given range is " //                        &
                   TRIM(to_string(MINVAL(rivers_lat2d))) // " to " //         &
                   TRIM(to_string(MAXVAL(rivers_lat2d))))

  IF ( ANY(rivers_lon2d < -180.0) .OR. ANY(rivers_lon2d > 360.0) )            &
     CALL log_fatal("init_rivers_props",                                      &
                   "Longitude is out of range - allowed range is " //         &
                   "-180.0 to 360.0, given range is " //                      &
                   TRIM(to_string(MINVAL(rivers_lon2d))) // " to " //         &
                   TRIM(to_string(MAXVAL(rivers_lon2d))))

  ! Checks for regular latitude/longitude grid
  IF ( rivers_reglatlon ) THEN

    ! Check cells are 'square' if grid defined as regular
    tiny_0 = TINY(0.0)
    IF ( ABS(rivers_dlat - rivers_dlon) > tiny_0 ) THEN
      CALL log_fatal("init_rivers_props",                                     &
         "Error in routing - latitude and longitude grid steps " //           &
         " must be equal if set as regular grid. rivers_dlat =" //            &
         TRIM(to_string(rivers_dlat)) // " , rivers_dlon =" //                &
         TRIM(to_string(rivers_dlon)))
    END IF

  ELSE

    ! Must be regular for any regridding methods
    IF ( rivers_regrid ) THEN
      WRITE(jules_message,*) "ERROR: init_rivers: routing codes with " //     &
                             "regridding only currently exist for " //        &
                             "regular lat/lon grids."
      CALL jules_print('init_rivers_props',jules_message)
      CALL log_fatal("init_rivers_props",                                     &
                     "Error in routing - non-regular lat/lon river grid " //  &
                     "defined and regridding requested. This is not a " //    &
                     "valid option currently.")

    END IF

    ! Check if rivers_dx set
    IF ( rivers_dx <= 0 ) THEN
      WRITE(jules_message,*) "ERROR: init_rivers: routing codes with " //     &
                           "non-regular lat/lon grid require setting " //     &
                           "rivers_dx in the namelist."
      CALL jules_print('init_rivers_props',jules_message)
      CALL log_fatal("init_rivers_props",                                     &
                   "Error in routing - non-regular lat/lon river grid " //    &
                   "defined but no grid dimension size provided.")
    END IF

  END IF
   
  !---------------------------------------------------------------------------
  ! Check if routing and model grids are same (required for RFM implementation)
  ! Note some lat/lon rotation may still be required
  !---------------------------------------------------------------------------

  IF ( rivers_reglatlon ) THEN
    ! n.b. in parallel mode, this check is not comprehensive (there might be
    ! a different gap around the edges of the regions given to each task)
    DO i = 1,SIZE(latitude, dim = 1)
      DO j = 1,SIZE(longitude, dim = 2)
        IF ( MOD(latitude(i,j) - reg_lat1, ABS(reg_dlat) ) > tiny_0 .OR.      &
             MOD(longitude(i,j) - reg_lon1, ABS(reg_dlon)) > tiny_0 ) THEN
          WRITE(jules_message,*) "Land grid points do not coincide " //       &
                              "with calculated main model grid."
          CALL jules_print('init_rivers_props',jules_message)
          CALL log_fatal("init_rivers_props",                                 &
            "Error in routing - model latitude and longitude grid " //        &
            "steps must be equal if set as regular grid. reg_dlat =" //       &
            TRIM(to_string(reg_dlat)) // " , reg_dlon =" //                   &
            TRIM(to_string(reg_dlon)))
        END IF
      END DO
    END DO
  END IF

  ! Checks for compatible dimension sizes if not regridding

  IF ( .NOT. rivers_regrid ) THEN

    IF ( nx_rivers /= nx_grid .OR. ny_rivers /= ny_grid ) THEN
      CALL jules_print('init_rivers_props',jules_message)
      CALL log_fatal("init_rivers_props",                                     &
         "Error in routing - model grid and river grid dimensions " //        &
         " must be equal size if no regridding attempted.")
    END IF
  END IF

  !---------------------------------------------------------------------------
  ! Initialise river gridbox areas (m2) and set grid spacing (m)
  !---------------------------------------------------------------------------

  IF ( rivers_dlat == 0 ) THEN
    rivers_dlat = rivers_lat2d(2,2) - rivers_lat2d(1,1)
  END IF
  IF ( rivers_dlon == 0 ) THEN
    rivers_dlon = rivers_lon2d(2,2) - rivers_lon2d(1,1)
  END IF

  IF ( rivers_reglatlon ) THEN

    IF ( rivers_dx < 0 ) THEN
      rivers_dx = planet_radius * (ABS(rivers_dlat) * pi_over_180)
    END IF

    dlat = 0.5 * rivers_dlat
    dlon = 0.5 * rivers_dlon

    DO ip = 1,np_rivers
      rivers_boxareas_rp(ip) = ABS( rivers_earth_area(                        &
                                               rivers_lat_rp(ip) - dlat,      &
                                               rivers_lat_rp(ip) + dlat,      &
                                               rivers_lon_rp(ip) - dlon,      &
                                               rivers_lon_rp(ip) + dlon ))
    END DO

  ELSE

    rivers_boxareas_rp(:) = rivers_dx * rivers_dx

  END IF

  !---------------------------------------------------------------------------
  ! Define overbank inundation initial variables for bankfull on river points
  ! if using Rosgen entrenchment option.
  ! Bankfull discharge (qbf) from power-law relationship from
  ! "Flood Modeling, Prediction and Mitigation" by Şen 2018.
  ! Bankfull width and depth (wbf, dbf) from Leopold & Maddock (1953).
  !---------------------------------------------------------------------------
  IF ( l_riv_overbank .AND. use_rosgen ) THEN
    DO ip = 1,np_rivers
      contribarea = (rivers_boxareas_rp(ip) * m2tokm2) *                      &
                    ( 1.0 + MAX(0.0, rfm_iarea_rp(ip)) )
      qbf(ip) = coef_b * ( contribarea**exp_c )
      wbf(ip) = riv_a * ( qbf(ip)**riv_b )
      dbf(ip) = riv_c * ( qbf(ip)**riv_f )
    END DO
  END IF

  !---------------------------------------------------------------------------
  ! Calculate land_pts to river_pts remappings ahead of main run
  !---------------------------------------------------------------------------
  IF ( rivers_regrid ) THEN
    CALL rivers_remap_unmatch()
  ELSE IF ( global_land_pts /= np_rivers ) THEN
    CALL rivers_remap_match()
  END IF

  !---------------------------------------------------------------------------
  ! Reset saved JULES grid from river to full land model
  !---------------------------------------------------------------------------

  dummy_grid = grid_create(local_grid%is_1d,local_grid%dim_name,              &
                           local_grid%nx,local_grid%x_name,local_grid%nx,     &
                           local_grid%y_name,local_grid%ny)
  use_subgrid = use_sub_local

END IF ! (check global_landpoint > 1)

END SUBROUTINE init_rivers_props
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

SUBROUTINE init_cable_progs()

USE io_constants, ONLY: max_sdf_name_len, max_file_name_len, namelist_unit

USE string_utils_mod, ONLY: to_string

USE templating_mod, ONLY: tpl_has_var_name, tpl_substitute_var

USE model_interface_mod, ONLY: identifier_len, populate_var, get_var_id

USE ancil_info, ONLY: land_pts!, soil_pts, soil_index, sm_levels


IMPLICIT NONE

!------------------------------------------------------------------------------
! Description:
!   
!   Reads in information about CABLE prognostic variables for their
!   initialisation
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!------------------------------------------------------------------------------
! Work variables
INTEGER, PARAMETER :: nCABLE_VARS  = 10

INTEGER :: nvars_required      ! The number of prognostic variables that are
                               ! required in this configuration
  
CHARACTER(LEN=identifier_len) :: required_vars(nCABLE_VARS)
                               ! The variable identifiers of the required
                               ! variables

INTEGER :: nvars_file       ! The number of variables that will be set
                            ! from the given file (template?)

! Variables passed to fill_variables_from_file
CHARACTER(LEN=identifier_len) :: file_var(nCABLE_VARS)
                      ! The variable identifiers of the variables to set
                      ! from file
CHARACTER(LEN=max_sdf_name_len) :: file_var_name(nCABLE_VARS)
                      ! The name of each variable in the file

CHARACTER(LEN=max_sdf_name_len) :: file_tpl_name(nCABLE_VARS)
                      ! The name to substitute in a template for each
                      ! variable

INTEGER :: i,l  ! Index variables

INTEGER :: error  ! Error indicator

!-----------------------------------------------------------------------------
! Definition of the cable_progs namelist
!-----------------------------------------------------------------------------
CHARACTER(LEN=max_file_name_len) :: FILE
                      ! The name of the file (or variable name template) to
                      ! use for variables that need to be filled from file

INTEGER :: nvars      ! The number of variables in this section

CHARACTER(LEN=identifier_len) :: var(nCABLE_VARS)
                      ! The variable identifiers of the variables
LOGICAL :: use_file(nCABLE_VARS)
                      !   T - the variable uses the file
                      !   F - the variable is set using a constant value
CHARACTER(LEN=max_sdf_name_len) :: var_name(nCABLE_VARS)
                      ! The name of each variable in the file
CHARACTER(LEN=max_sdf_name_len) :: tpl_name(nCABLE_VARS)
                      ! The name to substitute in a template for each
                      ! variable
REAL :: const_val(nCABLE_VARS)
INTEGER:: iconst_val(nCABLE_VARS)
                      ! The constant value to use for each variable if
                      ! use_file = F for that variable

NAMELIST  / cable_progs/ FILE, nvars, use_file, var, var_name

!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
nvars_required = 0
nvars_file     = 0
nvars          = 0
use_file(:)    = .TRUE.  ! Default is for every variable to be read from file

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
CALL log_info("init_cable_progs", "Reading CABLE_PROGS namelist...")

! First, we read the cable_progs namelist
READ(namelist_unit, NML = cable_progs, IOSTAT = error)

IF ( error /= 0 )                                                             &
   CALL log_fatal("init_cable_progs",                                         &
                  "Error reading namelist CABLE_PROGS" //                     &
                  "(IOSTAT=" // TRIM(to_string(error)) // ")")

!-----------------------------------------------------------------------------
! Set up CABLE prognostics using namelist values
!-----------------------------------------------------------------------------
! Set up the required variables
! All the CABLE variables are always required for CABLE runs
nvars_required = nCABLE_VARS
required_vars(:) = (/                                                         &
                        'ThreeLayerSnowFlag_CABLE',                           &
                        'OneLyrSnowDensity_CABLE ',                           &
                        'SnowAge_CABLE           ',                           &
                        'SnowDensity_CABLE       ',                           &
                        'SnowMass_CABLE          ',                           &
                        'SnowDepth_CABLE         ',                           &
                        'SnowTemp_CABLE          ',                           &
                        'FrozenSoilFrac_CABLE    ',                           &
                        'SoilMoisture_CABLE      ',                           &
                        'SoilTemp_CABLE          '                            &
                        /)

!-----------------------------------------------------------------------------
! Check that all the required variables are there
!-----------------------------------------------------------------------------

DO i = 1,nvars_required
  IF ( .NOT. ANY(var(1:nvars) == TRIM(required_vars(i))) )                    &
  !IF ( trim( var(1) ) /= TRIM( required_vars(1) ) )                         &
    CALL log_fatal("init_cable_progs",                                        &
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
      CALL log_info("init_cable_progs",                                       &
                    "'" // TRIM(var(i)) // "' will be read from file")

      ! If the variable will be filled from file, register it here
      nvars_file = nvars_file + 1
      file_var(nvars_file) = var(i)
      file_var_name(nvars_file) = var_name(i)
      file_tpl_name(nvars_file) = tpl_name(i)
    ELSE
      ! If the variable is being set as a constant, just populate it here
      CALL log_info("init_cable_progs",                                       &
                    "'" // TRIM(var(i)) // "' will be set to a constant")

      CALL populate_var(get_var_id(var(i)), const_val = const_val(i))
    END IF
  ELSE
    ! If the variable is not a required variable, warn about not using it
    CALL log_warn("init_cable_progs",                                         &
                  "Provided variable '" // TRIM(var(i)) //                    &
                  "' is not required, so will be ignored")
  END IF
END DO

!-----------------------------------------------------------------------------
! Set variables from file
!-----------------------------------------------------------------------------
IF ( nvars_file > 0 ) THEN
  IF ( tpl_has_var_name(FILE) ) THEN
    ! We are using a file name template, so loop through the variables setting
    ! one from each file
    DO i = 1,nvars_file
        
      CALL fill_variables_from_file(                                          &
        tpl_substitute_var(FILE, file_tpl_name(i)),                           &
        (/ file_var(i) /), (/ file_var_name(i) /)                             &
      )
    END DO
  ELSE
    ! We are not using a file name template, so set all variables from the same
    ! file
       
    CALL fill_variables_from_file(                                            &
      FILE,file_var(1:nvars_file), file_var_name(1:nvars_file)                &
    )
  END IF
END IF

RETURN

END SUBROUTINE init_cable_progs

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

END MODULE init_ancillaries_mod

