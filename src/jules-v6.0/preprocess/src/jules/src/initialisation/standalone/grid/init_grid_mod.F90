! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE init_grid_mod

USE logging_mod, ONLY: log_info, log_warn, log_fatal

IMPLICIT NONE

PRIVATE
PUBLIC init_grid

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


SUBROUTINE init_grid(nml_dir, crop_vars_data,psparms_data,top_pdm_data,       &
                              fire_vars_data,ainfo_data,trif_vars_data,       &
                              soil_ecosse_vars_data, aero_data,               &
                              urban_param_data,progs_data,trifctl_data,       &
                              coastal_data,jules_vars_data)

USE io_constants, ONLY: namelist_unit

USE string_utils_mod, ONLY: to_string

USE errormessagelength_mod, ONLY: errormessagelength

!TYPE definitions
USE crop_vars_mod, ONLY: crop_vars_data_type
USE p_s_parms,     ONLY: psparms_data_type
USE top_pdm,       ONLY: top_pdm_data_type
USE fire_vars_mod, ONLY: fire_vars_data_type
USE ancil_info,    ONLY: ainfo_data_type
USE trif_vars_mod, ONLY: trif_vars_data_type
USE soil_ecosse_vars_mod, ONLY: soil_ecosse_vars_data_type
USE aero,          ONLY: aero_data_type
USE urban_param_mod, ONLY:urban_param_data_type
USE prognostics, ONLY: progs_data_type
USE trifctl,       ONLY: trifctl_data_type
USE coastal,       ONLY: coastal_data_type
USE jules_vars_mod, ONLY: jules_vars_data_type

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in information about the model grids
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
TYPE(crop_vars_data_type), INTENT(IN OUT) :: crop_vars_data
TYPE(psparms_data_type), INTENT(IN OUT) :: psparms_data
TYPE(top_pdm_data_type), INTENT(IN OUT) :: top_pdm_data
TYPE(fire_vars_data_type), INTENT(IN OUT) :: fire_vars_data
TYPE(ainfo_data_type), INTENT(IN OUT) :: ainfo_data
TYPE(trif_vars_data_type), INTENT(IN OUT) :: trif_vars_data
TYPE(soil_ecosse_vars_data_type), INTENT(IN OUT) :: soil_ecosse_vars_data
TYPE(aero_data_type), INTENT(IN OUT) :: aero_data
TYPE(urban_param_data_type), INTENT(IN OUT) :: urban_param_data
TYPE(progs_data_type), INTENT(IN OUT) :: progs_data
TYPE(coastal_data_type), INTENT(IN OUT) :: coastal_data
TYPE(jules_vars_data_type), INTENT(IN OUT) :: jules_vars_data

! Work variables
INTEGER :: error  ! Error indicator
CHARACTER(LEN=errormessagelength) :: iomessage

TYPE(trifctl_data_type), INTENT(IN OUT) :: trifctl_data

!-----------------------------------------------------------------------------


! Open the grid namelist file
OPEN(namelist_unit, FILE=(TRIM(nml_dir) // '/' // 'model_grid.nml'),          &
               STATUS='old', POSITION='rewind', ACTION='read', IOSTAT = error,&
               IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_grid",                                                 &
                 "Error opening namelist file model_grid.nml " //             &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

! Defer to specialised routines to initialise the different aspect of the grid
CALL init_input_grid()
CALL init_latlon()
CALL init_land_frac()
CALL init_model_grid(crop_vars_data,psparms_data,top_pdm_data,                &
                     fire_vars_data,ainfo_data,trif_vars_data,                &
                     soil_ecosse_vars_data, aero_data, urban_param_data,      &
                     progs_data,trifctl_data, coastal_data,jules_vars_data)
CALL init_surf_hgt(jules_vars_data)
CALL init_z_land(ainfo_data,jules_vars_data)


CLOSE(namelist_unit, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_grid",                                                 &
                 "Error closing namelist file model_grid.nml " //             &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

RETURN

END SUBROUTINE init_grid
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_input_grid()

USE mem_brams_jules, ONLY: nxB,nyB  !DSM

USE io_constants, ONLY: namelist_unit, max_sdf_name_len

USE grid_utils_mod, ONLY: grid_create

USE string_utils_mod, ONLY: to_string

USE input_mod, ONLY: input_grid => grid

USE model_interface_mod, ONLY: pft_dim_name, cpft_dim_name, nvg_dim_name,     &
                                type_dim_name,                                &
                                tile_dim_name, snow_dim_name, soil_dim_name,  &
                                scpool_dim_name, bedrock_dim_name,            &
                                sclayer_dim_name

USE time_varying_input_mod, ONLY: time_dim_name
  
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in information about the grid for input data and checks for
!   consistency
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
INTEGER :: error  ! Error indicator
CHARACTER(LEN=errormessagelength) :: iomessage

!-----------------------------------------------------------------------------
! Definition of the jules_input_grid namelist - this combines variables
! from input_mod with local variables
!-----------------------------------------------------------------------------
LOGICAL :: grid_is_1d  ! T - grid is a vector
                       ! F - grid is 2D (e.g. lat/lon, but not necessarily)

! Used if grid_is_1d = T - the name and size of the single grid dimension
CHARACTER(LEN=max_sdf_name_len) :: grid_dim_name
INTEGER :: npoints

! Used if grid_is_1d = F - the names and sizes of the x and y dimensions
CHARACTER(LEN=max_sdf_name_len) :: x_dim_name, y_dim_name
INTEGER :: nx, ny

NAMELIST  / jules_input_grid/ grid_is_1d,                                     &
! Required options for input on a 1d grid
                              grid_dim_name, npoints,                         &
! Required options for input on a 2d grid
                              x_dim_name, y_dim_name, nx, ny,                 &
! Levels dimension names
                              pft_dim_name, cpft_dim_name, nvg_dim_name,      &
                              type_dim_name,                                  &
                              tile_dim_name, snow_dim_name, soil_dim_name,    &
                              scpool_dim_name, bedrock_dim_name,              &
                              sclayer_dim_name,                               &
! Name of the time dimension for time-varying files
                              time_dim_name


!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
grid_is_1d = .FALSE.

grid_dim_name = "land"
npoints       = 0

x_dim_name = "x"
y_dim_name = "y"
nx         = 0
ny         = 0

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
CALL log_info("init_input_grid", "Reading JULES_INPUT_GRID namelist...")

! First, we read the namelist
!print*,'1sssssnx=',nx
READ(namelist_unit, NML = jules_input_grid, IOSTAT = error, IOMSG = iomessage)
!print*,'2sssssnx=',nx
nx=nxB  !DSM
ny=nyB  !DSM
!print*,'Reading nx and ny from BRAMS',' nx=',nx,' ny=',ny
!stop

IF ( error /= 0 )                                                             &
  CALL log_fatal("init_input_grid",                                           &
                 "Error reading namelist JULES_INPUT_GRID " //                &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

!-----------------------------------------------------------------------------
! Build the input grid object from the namelist values
!-----------------------------------------------------------------------------
input_grid = grid_create(                                                     &
  grid_is_1d, grid_dim_name, npoints, x_dim_name, nx, y_dim_name, ny          &
)

! Check that a grid of non-zero size has been specified
IF ( input_grid%nx * input_grid%ny <= 0 )                                     &
  CALL log_fatal("init_input_grid",                                           &
                 "Invalid dimensions specified for grid - " //                &
                 "nx =  " // TRIM(to_string(input_grid%nx)) // "; " //        &
                 "ny = " // TRIM(to_string(input_grid%ny)))

CALL log_info("init_input_grid",                                              &
              "Size of input grid - " // TRIM(to_string(input_grid%nx)) //    &
              " x " // TRIM(to_string(input_grid%ny)))

RETURN

END SUBROUTINE init_input_grid
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

SUBROUTINE init_latlon()

USE mem_brams_jules, ONLY: glatB,glonB !JWA

USE dump_mod, ONLY: ancil_dump_read

USE errormessagelength_mod, ONLY: errormessagelength

USE input_mod, ONLY: input_grid => grid, fill_variables_from_file

USE io_constants, ONLY: max_sdf_name_len, max_file_name_len, namelist_unit

USE jules_water_resources_mod, ONLY: l_water_irrigation

USE missing_data_mod, ONLY: rmdi

USE model_grid_mod, ONLY: grid_area_ij, grid_lat => latitude,                 &
      grid_lon => longitude

USE model_interface_mod, ONLY: identifier_len, populate_var, get_var_id

USE string_utils_mod, ONLY: to_string

USE templating_mod, ONLY: tpl_has_var_name, tpl_substitute_var

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises the latitude and longitude of each gridbox.
!   Optionally initialises gridbox area.
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

INTEGER, PARAMETER :: max_vars = 3
  ! The maximum possible number of variables that can be given.

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_LATLON'

!-----------------------------------------------------------------------------
! Local scalar variables:
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  error,                                                                      &
    ! Error code.
  error_sum,                                                                  &
    ! Error indicator.
  i,                                                                          &
    ! Loop counter.
  nvars_file,                                                                 &
    ! The number of variables that will be set from the given file.
  nvars_required
    ! The number of variables that are required in this configuration.

CHARACTER(LEN=errormessagelength) :: iomessage

!-----------------------------------------------------------------------------
! Local array variables:
!-----------------------------------------------------------------------------
CHARACTER(LEN=max_sdf_name_len) :: file_tpl_name(max_vars)
     ! The name to substitute in a template for each variable.
CHARACTER(LEN=identifier_len) :: file_var(max_vars)
     ! The variable identifiers of the variables to set from file.
CHARACTER(LEN=max_sdf_name_len) :: file_var_name(max_vars)
     ! The name of each variable in the file
CHARACTER(LEN=identifier_len) :: required_vars(max_vars)
     ! The variable identifiers of the required variables.

!-----------------------------------------------------------------------------
! Definition of the jules_latlon namelist
!-----------------------------------------------------------------------------
INTEGER :: nvars
  ! The number of variables to be set.

LOGICAL :: read_from_dump
  ! Switch to read values from dump file.

CHARACTER(LEN=max_file_name_len) :: FILE
  ! The name of the file (or variable name template) to use for variables
  ! that need to be filled from file.

REAL(KIND=real_jlslsm) :: const_val(max_vars)
  ! The constant value to use for each variable if use_file = F for that
  ! variable.

LOGICAL :: use_file(max_vars)
  ! T - the variable uses the file
  ! F - the variable is set using a constant value

CHARACTER(LEN=identifier_len) :: var(max_vars)
  ! The identifiers of the variables.

CHARACTER(LEN=max_sdf_name_len) :: var_name(max_vars)
  ! The name of each variable in the file.

CHARACTER(LEN=max_sdf_name_len) :: tpl_name(max_vars)
  ! The name to substitute in a template for each variable.

NAMELIST  / jules_latlon/ read_from_dump, FILE, nvars, var,                   &
                          use_file, var_name, tpl_name, const_val

!-----------------------------------------------------------------------------
!end of header

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

CALL log_info(RoutineName,                                                    &
              "Getting latitude and longitude for the full input grid...")

! Allocate the model lat and lon arrays to the full input grid (for now).
! For ease we always allocate grid_area_ij, though it is not always needed.
error_sum = 0
ALLOCATE(grid_area_ij(input_grid%nx,input_grid%ny), stat = error)
error_sum = error_sum + error
ALLOCATE(grid_lat(input_grid%nx,input_grid%ny), stat = error)
error_sum = error_sum + error
ALLOCATE(grid_lon(input_grid%nx,input_grid%ny), stat = error)
error_sum = error_sum + error
IF ( error_sum > 0 ) THEN
  CALL log_fatal(RoutineName,                                                 &
                 "Error allocating arrays for latitude and longitude")
END IF

!-----------------------------------------------------------------------------
! Read the namelist.
!-----------------------------------------------------------------------------
CALL log_info(RoutineName, "Reading JULES_LATLON namelist...")

READ(namelist_unit, NML = jules_latlon, IOSTAT = error, IOMSG = iomessage)

IF ( error /= 0 ) THEN
 CALL log_fatal(RoutineName,                                                 &
                "Error reading namelist JULES_LATLON " //                    &
                "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                TRIM(iomessage) // ")")
END IF

ancil_dump_read%latlon = read_from_dump

!-----------------------------------------------------------------------------
! Identify which variables will be read here.
!-----------------------------------------------------------------------------
! First, set up the variables that are always required.
! latitude and longitude are always required here, regardless of the value of
! read_from_dump.
nvars_required = 2
required_vars(1:nvars_required) = (/ 'latitude ',                             &
                                     'longitude' /)

! If irrigation is used, add gridbox area (unless it will be read from dump).
IF ( l_water_irrigation ) THEN
  IF ( .NOT. ancil_dump_read%latlon ) THEN
    nvars_required = nvars_required + 1
    required_vars(nvars_required) =    'grid_area'
  END IF
ELSE
  ! Clarify that no variables from this section need to be read from dump.
  ancil_dump_read%latlon = .FALSE.
END IF

!-----------------------------------------------------------------------------
! Set variables.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Check that all the required variables are provided.
!-----------------------------------------------------------------------------
DO i = 1,nvars_required
  IF ( .NOT. ANY(var(1:nvars) == required_vars(i)) ) THEN
    CALL log_fatal(RoutineName,                                               &
                   "No value given for required variable '" //                &
                   TRIM(required_vars(i)) // "'")
  END IF
END DO

!-----------------------------------------------------------------------------
! Check which variables we will be using and partition them into variables
! set to constant values and variables set from file.
!-----------------------------------------------------------------------------
DO i = 1,nvars

  !---------------------------------------------------------------------------
  ! If the variable is required, we will be using it.
  !---------------------------------------------------------------------------
  IF ( ANY( required_vars(1:nvars_required) == TRIM(var(i)) ) ) THEN

    IF ( use_file(i) ) THEN
      ! This variable will be read from file.
      CALL log_info(RoutineName,                                              &
                    "'" // TRIM(var(i)) // "' will be read from file")

      ! If the variable will be filled from file, register it here.
      nvars_file = nvars_file + 1
      file_var(nvars_file) = var(i)
      file_var_name(nvars_file) = var_name(i)
      file_tpl_name(nvars_file) = tpl_name(i)
    ELSE
      ! THis variables will be set as a constant; populate it here.
      ! First check that a value has been provided.
      IF ( ABS( const_val(i) - rmdi ) < EPSILON(1.0) ) THEN
        CALL log_fatal(RoutineName,                                           &
                       "No constant value provided for variable '"            &
                       // TRIM(var(i)) // "'" )
      END IF
      CALL log_info(RoutineName,                                              &
                    "'" // TRIM(var(i)) // "' will be set to a " //           &
                    "constant = " // to_string(const_val(i)))
      CALL populate_var(get_var_id(var(i)), const_val = const_val(i))
    END IF

  ELSE

    ! The variable is not a required variable. Warn about not using it.
    CALL log_warn(RoutineName,                                                &
                  "Provided variable '" // TRIM(var(i)) //                    &
                  "' is not required, so will be ignored")
  END IF

END DO  !  i (loop over variables)

!-----------------------------------------------------------------------------
! Set variables from file
!-----------------------------------------------------------------------------
!DSM IF ( nvars_file > 0 ) THEN

!DSM  ! Check that a file name was provided.
!DSM  IF ( LEN_TRIM(FILE) == 0 ) THEN
!DSM    CALL log_fatal(RoutineName, "No file name provided")
!DSM  END IF

!DSM  IF ( tpl_has_var_name(FILE) ) THEN
!DSM    ! We are using a file name template, so loop through the variables
!DSM    ! setting one from each file.
!DSM    DO i = 1,nvars_file
!DSM      ! Check that a template string was provided for this variable.
!DSM      IF ( LEN_TRIM(file_tpl_name(i)) == 0 ) THEN
!DSM        CALL log_fatal( RoutineName,                                          &
!DSM                        "No variable name template substitution " //          &
!DSM                        "provided for " // TRIM(file_var(i)) )
!DSM      END IF
!DSM      CALL fill_variables_from_file(                                          &
!DSM                               tpl_substitute_var(FILE, file_tpl_name(i)),    &
!DSM                               (/ file_var(i) /), (/ file_var_name(i) /) )
!DSM    END DO
!DSM  ELSE
!DSM    ! We are not using a file name template, so set all variables from the
!DSM    ! same file.
!DSM    CALL fill_variables_from_file( FILE, file_var(1:nvars_file),              &
!DSM                                   file_var_name(1:nvars_file) )
!DSM  END IF


!DSM END IF  !  nvars_file

!DSM !-----------------------------------------------------------------------------
!DSM ! Advise if some values will be read from the dump file.
!DSM !-----------------------------------------------------------------------------
!DSM IF ( ancil_dump_read%latlon ) THEN
!DSM  CALL log_info(RoutineName,                                                  &
!DSM                "grid_area will be read from the dump file. " //              &
!DSM                "Namelist values will be ignored.")

!DSM END IF

!print*, 'Data is on a grid - reading latitude and longitude from BRAMS'
grid_lat=glatB
grid_lon=glonB

!-----------------------------------------------------------------------------
! Check that the values of latitude and longitude seem sensible.
!-----------------------------------------------------------------------------
IF ( ANY(grid_lat < -90.0) .OR. ANY(grid_lat > 90.0) )                        &
  CALL log_fatal(RoutineName,                                                 &
                 "Latitude is out of range - allowed range is -90.0 to " //   &
                 "90.0, given range is " //                                   &
                 TRIM(to_string(MINVAL(grid_lat))) // " to " //               &
                 TRIM(to_string(MAXVAL(grid_lat))))

IF ( ANY(grid_lon < -180.0) .OR. ANY(grid_lon > 360.0) )                      &
  CALL log_fatal(RoutineName,                                                 &
                 "Longitude is out of range - allowed range is -180.0 " //    &
                 "to 360.0, given range is " //                               &
                 TRIM(to_string(MINVAL(grid_lon))) // " to " //               &
                 TRIM(to_string(MAXVAL(grid_lon))))

RETURN

END SUBROUTINE init_latlon
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_model_grid(crop_vars_data,psparms_data,top_pdm_data,          &
                           fire_vars_data,ainfo_data,trif_vars_data,          &
                           soil_ecosse_vars_data, aero_data, urban_param_data,&
                           progs_data,trifctl_data, coastal_data,             &
                           jules_vars_data)

USE allocate_jules_arrays_mod, ONLY: allocate_jules_arrays

USE io_constants, ONLY: max_file_name_len, namelist_unit, points_file_unit

USE string_utils_mod, ONLY: to_string

USE grid_utils_mod, ONLY: subgrid_info, grid_create, subgrid_create,          &
                           subgrid_extract, subgrid_restrict

USE input_mod, ONLY: grid_in => grid,                                         &
! Variables to do with the extraction of subgrids if specified
! use_subgrid does not have an _in so that it can appear as use_subgrid in the
! namelist
                        use_subgrid, subgrid_in => subgrid

USE output_mod, ONLY: grid_out => grid,                                       &
! Variables to do with writing subgrids if multiple tasks are used
                         use_subgrid_out => use_subgrid,                      &
                         subgrid_out => subgrid

USE parallel_mod, ONLY: decompose_domain

USE switches, ONLY: l_co2_interactive

USE jules_soil_biogeochem_mod, ONLY: soil_model_1pool, soil_model_ecosse,     &
                                      soil_model_rothc, l_layeredC,           &
                                      soil_bgc_model, dim_ch4layer

USE jules_surface_mod, ONLY: l_point_data, l_aggregate

USE jules_vegetation_mod, ONLY: l_triffid, l_crop

USE model_grid_mod, ONLY: model_grid, global_land_pts, global_land_mask,      &
                          grid_area_ij, latitude, longitude,                  &
                          latitude_of_land_pts, longitude_of_land_pts

USE coastal, ONLY: flandg

USE atm_fields_bounds_mod, ONLY: atm_fields_bounds_init
USE theta_field_sizes, ONLY: t_i_length, t_j_length,                          &
                             u_i_length,u_j_length,                           &
                             v_i_length,v_j_length

USE ancil_info, ONLY: land_pts, co2_dim_len,                                  &
                       co2_dim_row, dim_cs1, dim_cs2,                         &
                       nsurft, dim_cslayer,                                   &
                       row_length, rows, dim_soil_n_pool, n_rows, nsoilt

USE jules_deposition_mod, ONLY: ndry_dep_species

USE jules_surface_types_mod, ONLY: npft, ncpft, nnvg, ntype

USE jules_snow_mod, ONLY: nsmax

USE jules_sea_seaice_mod, ONLY: nice, nice_use

USE jules_soil_mod, ONLY: sm_levels, ns_deep, l_tile_soil

USE model_interface_mod, ONLY: bl_level_dim_size, bedrock_dim_size,           &
                                pft_dim_size, cpft_dim_size, nvg_dim_size,    &
                                type_dim_size, tile_dim_size, snow_dim_size,  &
                                soil_dim_size, soilt_dim_size,                &
                                scpool_dim_size, ch4layer_dim_size,           &
                                sclayer_dim_size, tracer_dim_size

USE nlsizes_namelist_mod, ONLY: bl_levels, jules_nlsizes

USE errormessagelength_mod, ONLY: errormessagelength

USE um_types, ONLY: real_jlslsm

!TYPE definitions
USE crop_vars_mod, ONLY: crop_vars_data_type
USE p_s_parms,     ONLY: psparms_data_type
USE top_pdm,       ONLY: top_pdm_data_type
USE fire_vars_mod, ONLY: fire_vars_data_type
USE ancil_info,    ONLY: ainfo_data_type
USE trif_vars_mod, ONLY: trif_vars_data_type
USE soil_ecosse_vars_mod, ONLY: soil_ecosse_vars_data_type
USE aero,          ONLY: aero_data_type
USE urban_param_mod, ONLY:urban_param_data_type
USE prognostics,   ONLY: progs_data_type
USE trifctl,       ONLY: trifctl_data_type
USE coastal,       ONLY: coastal_data_type
USE jules_vars_mod, ONLY: jules_vars_data_type

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in information about how the model grid is defined, sets it up and
!   checks for consistency
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
!Arguments
!TYPES containing field data (IN OUT)
TYPE(crop_vars_data_type), INTENT(IN OUT) :: crop_vars_data
TYPE(psparms_data_type), INTENT(IN OUT) :: psparms_data
TYPE(top_pdm_data_type), INTENT(IN OUT) :: top_pdm_data
TYPE(fire_vars_data_type), INTENT(IN OUT) :: fire_vars_data
TYPE(ainfo_data_type), INTENT(IN OUT) :: ainfo_data
TYPE(trif_vars_data_type), INTENT(IN OUT) :: trif_vars_data
TYPE(soil_ecosse_vars_data_type), INTENT(IN OUT) :: soil_ecosse_vars_data
TYPE(aero_data_type), INTENT(IN OUT) :: aero_data
TYPE(urban_param_data_type), INTENT(IN OUT) :: urban_param_data
TYPE(progs_data_type), INTENT(IN OUT) :: progs_data
TYPE(trifctl_data_type), INTENT(IN OUT) :: trifctl_data
TYPE(coastal_data_type), INTENT(IN OUT) :: coastal_data
TYPE(jules_vars_data_type), INTENT(IN OUT) :: jules_vars_data

! Work variables
LOGICAL :: input_mask(grid_in%nx, grid_in%ny)
                                      ! Defined on the full input grid
                                      !   T - the point will be modelled
                                      !   F - the point will not be modelled
REAL(KIND=real_jlslsm) :: point_lat, point_lon
                              ! Latitude and longitude of points read from
                              ! points file

TYPE(subgrid_info) :: task_subgrid  ! The subgrid that the current MPI task
                                    ! will model, as a subgrid of the full
                                    ! model grid

REAL(KIND=real_jlslsm)                                                        &
    , ALLOCATABLE :: global_area(:,:), global_lat(:,:), global_lon(:,:), global_land_frac(:,:)
                             ! Latitude, longitude and land fraction on
                             ! the full model grid
                             ! These are not required outside of this
                             ! routine

INTEGER :: error, error_sum  ! Error indicators
CHARACTER(LEN=errormessagelength) :: iomessage

INTEGER :: i, i2, j, j2, l  ! Index variables


!-----------------------------------------------------------------------------
! Definition of the jules_model_grid namelist - this combines variables
! from input_mod with some local variables
!-----------------------------------------------------------------------------
LOGICAL :: land_only      ! T - only model land points
                          ! F - model all specified points


LOGICAL :: force_1d_grid  ! T - model grid is 1D
                          ! F - model grid is 1D unless input grid is 2D
                          !     and model grid is the whole input grid

LOGICAL :: latlon_region  ! T - subgrid is to be selected with latitude and
                          !     longitude bounds
                          ! F - subgrid is to be selected using a list
                          !     of latitudes and longitudes
REAL(KIND=real_jlslsm) :: lat_bounds(2)   ! USED IF latlon_region=T
REAL(KIND=real_jlslsm) :: lon_bounds(2)
                        ! Upper and lower bounds for latitude and longitude

INTEGER :: npoints = 0 ! The number of points to read from file
CHARACTER(LEN=max_file_name_len) :: points_file
                        ! The file to read latitudes and longitudes for
                        ! specified points

NAMELIST  / jules_model_grid/ land_only, use_subgrid, latlon_region,          &
                            lat_bounds, lon_bounds, npoints, points_file,     &
                            force_1d_grid


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
land_only = .TRUE.  ! Default is to model only land points from whatever
                    ! points are selected
force_1d_grid = .FALSE. ! Default is model grid is 1D unless input grid is 2D
                        ! and model grid is the whole input grid
points_file = ''    ! Empty file name.

!-----------------------------------------------------------------------------
! Read the init_model_grid namelist.
!-----------------------------------------------------------------------------
CALL log_info("init_model_grid", "Reading JULES_MODEL_GRID namelist...")

READ(namelist_unit, NML = jules_model_grid, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_model_grid",                                           &
                 "Error reading namelist JULES_MODEL_GRID " //                &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

!-----------------------------------------------------------------------------
! Read the jules_nlsizes namelist.
!-----------------------------------------------------------------------------
CALL log_info("init_model_grid", "Reading JULES_NLSIZES namelist...")

READ(namelist_unit, NML = jules_nlsizes, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 ) THEN
  CALL log_fatal("init_model_grid",                                           &
                 "Error reading namelist JULES_NLSIZES " //                   &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")
END IF

!-----------------------------------------------------------------------------
! Set up the model grid variables
!-----------------------------------------------------------------------------
CALL log_info("init_model_grid", "Setting up model grid variables...")

!*****************************************************************************
!*****************************************************************************
! Work out what points in the input grid will make up the model grid
!
! By the time we get here, we know that we have latitude, longitude and
! land fraction on the full input grid already read in
! Use these in combination with the values in the jules_model_grid namelist
! to work out what points to use
!*****************************************************************************
!*****************************************************************************
! To start with, assume we are using all the points
input_mask(:,:) = .TRUE.

IF ( use_subgrid ) THEN
  IF ( latlon_region ) THEN
    !-----------------------------------------------------------------------------
    ! The subgrid will be selected using latitude and longitude bounds
    !-----------------------------------------------------------------------------
    CALL log_info("init_model_grid",                                          &
                  "Subgrid will be selected using latitude and " //           &
                  "longitude bounds")
    CALL log_info("init_model_grid",                                          &
                  "Latitude range - " // TRIM(to_string(lat_bounds(1))) //    &
                  " to " // TRIM(to_string(lat_bounds(2))))
    CALL log_info("init_model_grid",                                          &
                  "Longitude range - " // TRIM(to_string(lon_bounds(1))) //   &
                  " to " // TRIM(to_string(lon_bounds(2))))

    input_mask = ( lat_bounds(1) <= latitude ) .AND.                          &
                 ( latitude <= lat_bounds(2) ) .AND.                          &
                 ( lon_bounds(1) <= longitude ) .AND.                         &
                 ( longitude <= lon_bounds(2) )
  ELSE
    !-----------------------------------------------------------------------------
    ! The subgrid will be selected using a list of latitudes and longitudes
    ! from the given file
    !-----------------------------------------------------------------------------
    CALL log_info("init_model_grid",                                          &
                  "Subgrid will be selected using a list of points " //       &
                  "from " // TRIM(points_file))

    !     Check that a file name was provided.
    IF ( LEN_TRIM(points_file) == 0 )                                         &
      CALL log_fatal("init_model_grid", "No points file provided" )

    ! No points have been specified yet, so set input_mask accordingly
    input_mask(:,:) = .FALSE.

    !-----------------------------------------------------------------------------
    ! Read the point latitudes and longitudes from file
    !-----------------------------------------------------------------------------
    OPEN(points_file_unit, FILE=points_file,                                  &
               STATUS='old', POSITION='rewind', ACTION='read', IOSTAT = error,&
               IOMSG = iomessage)
    IF ( error /= 0 )                                                         &
    CALL log_fatal("init_model_grid",                                         &
                   "Error opening points file " //                            &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")

    DO i = 1,npoints
      ! Try to read the next latitude/longitude pair
      READ(points_file_unit, *, IOSTAT = error, IOMSG = iomessage) point_lat, &
                                                               point_lon
      IF ( error /= 0 )                                                       &
        CALL log_fatal("init_model_grid",                                     &
                       "Error reading lat/lon pair from points file " //      &
                       "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //   &
                       TRIM(iomessage) // ")")

      ! Set input_mask to .TRUE. at that point
      WHERE ( ABS(latitude - point_lat) < EPSILON(1.0) .AND.                  &
              ABS(longitude - point_lon) < EPSILON(1.0) )
        input_mask = .TRUE.
      END WHERE
    END DO

    CLOSE(points_file_unit, IOSTAT = error, IOMSG = iomessage)
    IF ( error /= 0 )                                                         &
      CALL log_fatal("init_model_grid",                                       &
                     "Error closing points file " //                          &
                     "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //     &
                     TRIM(iomessage) // ")")
  END IF  ! latlon_region
END IF  ! use_subgrid

IF ( land_only ) THEN
  !-----------------------------------------------------------------------------
  ! If requested, select land points only from the points we have left
  !-----------------------------------------------------------------------------
  CALL log_info("init_model_grid",                                            &
                "From the points specified, only land points will be " //     &
                "modelled")

  input_mask = input_mask .AND. ( flandg > 0.0 )
END IF


!*****************************************************************************
!-----------------------------------------------------------------------------
! At this point, the .TRUE. points in input_mask will comprise the FULL
! model grid
! Issue an error if all points have been excluded
!-----------------------------------------------------------------------------
IF ( .NOT. ANY(input_mask) )                                                  &
  CALL log_fatal("init_model_grid",                                           &
                 "All points in input grid have been excluded from " //       &
                 "model grid")
!*****************************************************************************

!-----------------------------------------------------------------------------
! Build a subgrid from the input grid and mask
!-----------------------------------------------------------------------------
subgrid_in = subgrid_create(grid_in, input_mask, force_1d_grid)

! Construct the model grid to be the same size as the input subgrid
! It doesn't matter about dimension names as it will not be used for I/O
! It is always created as a 2D grid
model_grid = grid_create(                                                     &
  .FALSE., x_name = "", y_name = "", nx = subgrid_in%nx, ny = subgrid_in%ny   &
)

!-----------------------------------------------------------------------------
! Store the latitude, longitude and land fraction on the full model grid
!-----------------------------------------------------------------------------
ALLOCATE(global_area(model_grid%nx,model_grid%ny))
global_area(:,:) = subgrid_extract(subgrid_in, grid_area_ij)

ALLOCATE(global_lat(model_grid%nx,model_grid%ny))
global_lat(:,:) = subgrid_extract(subgrid_in, latitude)

ALLOCATE(global_lon(model_grid%nx,model_grid%ny))
global_lon(:,:) = subgrid_extract(subgrid_in, longitude)

ALLOCATE(global_land_frac(model_grid%nx,model_grid%ny))
global_land_frac(:,:) = subgrid_extract(subgrid_in, flandg)

! Work out how many land points are in the full model grid and set the land
! mask
ALLOCATE(global_land_mask(model_grid%nx,model_grid%ny))
global_land_mask(:,:) = ( global_land_frac(:,:) > EPSILON(1.0) )
global_land_pts = COUNT(global_land_mask)

!-----------------------------------------------------------------------------
! Check that we do not have any repeated coordinates.
!-----------------------------------------------------------------------------
DO i = 1,model_grid%nx
  DO j = 1,model_grid%ny

    ! Scan the rest of this row.
    DO i2 = i+1,model_grid%nx
      IF ( global_lat(i,j) == global_lat(i2,j) .AND.                          &
           global_lon(i,j) == global_lon(i2,j) ) THEN
        CALL log_fatal( "init_model_grid",                                    &
                        "Duplicate coordinates have been specified."  //      &
                        " Latitude=" // to_string(global_lat(i,j))    //      &
                        " Longitude=" // to_string(global_lon(i,j)) )
      END IF
    END DO

    ! Scan all higher rows.
    DO i2 = 1,model_grid%nx
      DO j2 = j+1,model_grid%ny
        IF ( global_lat(i,j) == global_lat(i2,j2) .AND.                       &
             global_lon(i,j) == global_lon(i2,j2) ) THEN
          CALL log_fatal( "init_model_grid",                                  &
                          "Duplicate coordinates have been specified."  //    &
                          " Latitude=" // to_string(global_lat(i,j))    //    &
                          " Longitude=" // to_string(global_lon(i,j)) )
        END IF
      END DO
    END DO

  END DO  !  j
END DO  !  i

!-----------------------------------------------------------------------------
! Decompose the model domain
!
! Returns the subgrid of the full model grid that the current task will model
!-----------------------------------------------------------------------------
task_subgrid = decompose_domain(model_grid)

!-----------------------------------------------------------------------------
! Set up the extraction of the input grid
!-----------------------------------------------------------------------------
! First, we combine the task subgrid with the subgrid generated from input_mask
! to get the input subgrid for this task
subgrid_in = subgrid_restrict(subgrid_in, task_subgrid)

! Do we actually have a subgrid, or does the subgrid cover the whole input grid?
use_subgrid = ( subgrid_in%nx /= grid_in%nx ) .OR.                            &
              ( subgrid_in%ny /= grid_in%ny )

IF ( use_subgrid ) THEN
  ! Regrid the latitude and longitude onto the model grid for this task

  DEALLOCATE(grid_area_ij)
  ALLOCATE(grid_area_ij(task_subgrid%nx,task_subgrid%ny))
  grid_area_ij(:,:) = subgrid_extract(task_subgrid, global_area)

  DEALLOCATE(latitude)
  ALLOCATE(latitude(task_subgrid%nx,task_subgrid%ny))
  latitude(:,:) = subgrid_extract(task_subgrid, global_lat)

  DEALLOCATE(longitude)
  ALLOCATE(longitude(task_subgrid%nx,task_subgrid%ny))
  longitude(:,:) = subgrid_extract(task_subgrid, global_lon)

  DEALLOCATE(flandg)
  ALLOCATE(flandg(task_subgrid%nx,task_subgrid%ny))
  flandg(:,:) = subgrid_extract(task_subgrid, global_land_frac)
END IF

! Deallocate some arrays that are no longer required
DEALLOCATE(global_area)
DEALLOCATE(global_lat)
DEALLOCATE(global_lon)
DEALLOCATE(global_land_frac)

!-----------------------------------------------------------------------------
! Set up the output grid
!
! The full output grid is the same as the full model grid
! The subgrid that this task is responsible for outputting is the same as
! the subgrid it is responsible for modelling
!-----------------------------------------------------------------------------
grid_out    = model_grid
subgrid_out = task_subgrid

! Do we actually have a subgrid, or does the subgrid cover the whole output grid?
use_subgrid_out = ( subgrid_out%nx /= grid_out%nx ) .OR.                      &
                  ( subgrid_out%ny /= grid_out%ny )

!-----------------------------------------------------------------------------
! Set up the model grid for this task, i.e. the grid that the science code
! knows about. Only control code has to deal with the concept of a 'task
! model grid' and a 'full model grid'
!-----------------------------------------------------------------------------
t_i_length = task_subgrid%nx
t_j_length = task_subgrid%ny
v_i_length = task_subgrid%nx
v_j_length = task_subgrid%ny
u_i_length = task_subgrid%nx
u_j_length = task_subgrid%ny

! At the moment, JULES uses no halos, and the t, p, u and v grids are the same
CALL atm_fields_bounds_init(0, 0, 0, 0, t_i_length, t_j_length, t_j_length, 1,&
                            bl_levels_opt = bl_levels)
! Copy the values into row_length, rows and n_rows for legacy purposes
row_length = t_i_length
rows   = t_j_length
n_rows = t_j_length

CALL log_info("init_model_grid",                                              &
              "Size of model grid - " // TRIM(to_string(t_i_length)) //       &
              " x " // TRIM(to_string(t_j_length)))

! Calculate the number of land points
land_pts = COUNT(flandg > EPSILON(1.0))
CALL log_info("init_model_grid",                                              &
              "Selected grid contains " // TRIM(to_string(land_pts)) //       &
              " land points")

CALL log_info("init_model_grid",                                              &
              "bl_levels = " // TRIM(to_string(bl_levels)))

! Calculate number of surface types based on whether we are using the aggregate
! surface type
IF ( l_aggregate ) THEN
  nsurft = 1
ELSE
  nsurft = ntype
END IF

! Set up the number of soil tiles
IF ( l_tile_soil ) THEN
  nsoilt = nsurft
ELSE
  nsoilt = 1
  !This is also the default value, but reaffirming it here for clarity
END IF

! Set up the ice categories - there is only one
nice = 1
! nice_use MUST EQUAL NICE FOR STANDALONE JULES
nice_use = nice

! Set dimensions for soil carbon and nitrogen model.
SELECT CASE ( soil_bgc_model )
CASE ( soil_model_1pool )
  dim_cs1         = 1  !  1 pool
  dim_soil_n_pool = 0  !  Not used, set to zero so no space used.
  dim_cs2         = 1  !  1 as this model is not used with TRIFFID
  IF (l_layeredC) THEN
    dim_cslayer = sm_levels
  ELSE
    dim_cslayer = 1
  END IF
CASE ( soil_model_rothc )
  dim_cs1         = 4   !  4 soil pools
  dim_soil_n_pool = 0   !  Not used, set to zero so no space used.
  dim_cs2         = land_pts
  IF (l_layeredC) THEN
    dim_cslayer = sm_levels
  ELSE
    dim_cslayer = 1
  END IF
CASE ( soil_model_ecosse )
  dim_cs1         = 4    !  4 soil organic C pools: DPM, RPM, BIO, HUM
  dim_soil_n_pool = 6    !  4 soil organic N pools, NH4, NO3
  dim_cs2         = land_pts
  ! dim_cslayer was set in jules_soil_biogeochem_mod.
END SELECT

IF ( l_co2_interactive ) THEN
  co2_dim_len = t_i_length
  co2_dim_row = t_j_length
ELSE
  co2_dim_len = 1
  co2_dim_row = 1
END IF

! Warn if using l_point_data for a non-single-point run
IF ( l_point_data .AND. t_i_length * t_j_length > 1 )                         &
  CALL log_warn("init_model_grid",                                            &
                "l_point_data is selected but there is more than one point")

!-----------------------------------------------------------------------------
! Allocate model arrays
!
! Note: CABLE model arrays are only allocated later on in init.F90 due to a
!       dependency on surft_pts, which has not yet been initialised
!-----------------------------------------------------------------------------

ALLOCATE(latitude_of_land_pts(land_pts))
ALLOCATE(longitude_of_land_pts(land_pts))

CALL allocate_jules_arrays(crop_vars_data,psparms_data,top_pdm_data,          &
                            fire_vars_data,ainfo_data,trif_vars_data,         &
                            soil_ecosse_vars_data, aero_data,                 &
                            urban_param_data, progs_data,trifctl_data,        &
                            coastal_data,jules_vars_data)

ainfo_data%land_mask(:,:) = ( flandg(:,:) > EPSILON(1.0) )

l = 0
DO j = 1,t_j_length
  DO i = 1,t_i_length
    IF ( ainfo_data%land_mask(i,j) ) THEN
      l = l + 1
      ainfo_data%land_index(l) = (j-1) * t_i_length + i
      coastal_data%fland(l) = flandg(i,j)
      latitude_of_land_pts(l)  = latitude(i,j)
      longitude_of_land_pts(l) = longitude(i,j)
    END IF
  END DO
END DO

!-----------------------------------------------------------------------------
! Now we know all the required info, set up the z-dimension sizes for IO
!-----------------------------------------------------------------------------
bl_level_dim_size = bl_levels
pft_dim_size      = npft
cpft_dim_size     = ncpft
nvg_dim_size      = nnvg
type_dim_size     = ntype
tile_dim_size     = nsurft
soilt_dim_size    = nsoilt
snow_dim_size     = nsmax
soil_dim_size     = sm_levels
scpool_dim_size   = dim_cs1
sclayer_dim_size  = dim_cslayer
ch4layer_dim_size = dim_ch4layer
tracer_dim_size   = ndry_dep_species
bedrock_dim_size  = ns_deep

RETURN

END SUBROUTINE init_model_grid
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_land_frac()

USE mem_brams_jules, ONLY: land_fracB !JWA

USE io_constants, ONLY: max_sdf_name_len, max_file_name_len, namelist_unit

USE string_utils_mod, ONLY: to_string

USE coastal, ONLY: flandg

USE input_mod, ONLY: input_grid => grid, fill_variables_from_file
  
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in the land fraction on the full grid and checks for consistency
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
INTEGER :: error  ! Error indicator
CHARACTER(LEN=errormessagelength) :: iomessage

!-----------------------------------------------------------------------------
! Definition of the jules_land namelist - this combines variables
! from input_mod with some local variables
!-----------------------------------------------------------------------------
CHARACTER(LEN=max_file_name_len) :: FILE
CHARACTER(LEN=max_sdf_name_len) :: land_frac_name
NAMELIST  / jules_land_frac/ FILE, land_frac_name

!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
!JWS FILE=''  !  Empty file name.
land_frac_name = ''  !  Empty variable name.

!JWS CALL log_info("init_land_frac", "Reading JULES_LAND_FRAC namelist...")

! First, we read the namelist
!JWS READ(namelist_unit, NML = jules_land_frac, IOSTAT = error, IOMSG = iomessage)
!JWS IF ( error /= 0 )                                                             &
!JWS  CALL log_fatal("init_land_frac",                                            &
!JWS                 "Error reading namelist JULES_LAND_FRAC " //                 &
!JWS                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
!JWS                 TRIM(iomessage) // ")")


!-----------------------------------------------------------------------------
! Verify for consistency
!-----------------------------------------------------------------------------
CALL log_info("init_land_frac",                                               &
              "Getting land fraction for the full input grid...")

! For now, we allocate flandg on the full input grid
ALLOCATE(flandg(input_grid%nx,input_grid%ny), stat = error)
IF ( error > 0 )                                                              &
  CALL log_fatal("init_land_frac", "Error allocating flandg")

IF ( input_grid%nx * input_grid%ny <= 1 ) THEN
  ! For single point data, just assume the single point is 100% land
  CALL log_info("init_land_frac",                                             &
                "Data is at a single point - setting land fraction to 1.0")
  flandg(1,1) = 1.0
ELSE
  ! For any other grid, read land fraction from file
!JWS  CALL log_info("init_land_frac",                                             &
!JWS                "Data is on a grid - reading land fraction from file " //     &
!JWS                TRIM(FILE))

Print*, '
open(unit=66,file='jules.log',position='append',status='old',action='write')
write(unit=66,fmt='(A)') 'Data is on a grid - reading land fraction from BRAMS'
close(unit=66)'

  !   Check that a file name was provided.
  
!JWS  IF ( LEN_TRIM(FILE) == 0 )                                                  &
!JWS    CALL log_fatal("init_land_frac", "No file name provided")

!JWS  CALL fill_variables_from_file(FILE, (/ 'land_fraction' /),                  &
!JWS                                        (/  land_frac_name   /))

   flandg(:,:)=land_fracB(:,:)   !JWS 
END IF

! For now, make sure that land fraction is either 1.0 or 0.0 - until JULES
! can deal with coastal tiling
WHERE ( flandg > EPSILON(1.0) )
  flandg = 1.0
ELSEWHERE
  flandg = 0.0
END WHERE

RETURN

END SUBROUTINE init_land_frac
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_surf_hgt(jules_vars_data)

USE io_constants, ONLY: max_sdf_name_len, max_file_name_len, namelist_unit

USE string_utils_mod, ONLY: to_string

USE input_mod, ONLY: fill_variables_from_file

USE jules_surface_mod, ONLY: l_aggregate

USE c_elevate, ONLY: surf_hgt_io, l_elev_absolute_height

USE ancil_info, ONLY: nsurft

USE errormessagelength_mod, ONLY: errormessagelength

!TYPE definitions
USE jules_vars_mod, ONLY: jules_vars_data_type

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises tile heights relative to the gridbox mean.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

!TYPES containing field data (IN OUT)
TYPE(jules_vars_data_type), INTENT(IN OUT) :: jules_vars_data

! Work variables
INTEGER :: error  ! Error indicator

!-----------------------------------------------------------------------------
! Definition of the jules_surf_hgt namelist
!-----------------------------------------------------------------------------
LOGICAL :: zero_height  ! T - set height at all points to 0.0
                        ! F - read height from file
LOGICAL :: use_file     ! T - the variable uses the file
                        ! F - the variable is set using a constant value
CHARACTER(LEN=max_file_name_len) :: FILE  ! If input grid has more than one
                                          ! point, read heights from this
                                          ! file
CHARACTER(LEN=max_sdf_name_len) :: surf_hgt_name  ! Use the variable with
                                                  ! this name
CHARACTER(LEN=errormessagelength) :: iomessage

NAMELIST  / jules_surf_hgt/ zero_height, l_elev_absolute_height, use_file,    &
                          FILE, surf_hgt_name, surf_hgt_io

!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
zero_height               = .TRUE.  ! Default is to use zero height everywhere
l_elev_absolute_height(:) = .FALSE. ! T - tiles have absolute heights above
                                    ! sea-level
                                    ! F - tiles have heights relative to the
                                    ! gridbox mean
use_file                  = .TRUE.  ! Default is to read variable from a file
FILE=''      ! Empty file name.
surf_hgt_name             = ''      ! Empty variable name.

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
CALL log_info("init_surf_hgt", "Reading JULES_SURF_HGT namelist...")

! First, we read the namelist
READ(namelist_unit, NML = jules_surf_hgt, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_surf_hgt",                                             &
                 "Error reading namelist JULES_SURF_HGT " //                  &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

!-----------------------------------------------------------------------------
! Set values derived from namelist and verify for consistency
!-----------------------------------------------------------------------------
! Zero height must be used with the aggregate surface scheme
IF ( l_aggregate ) THEN
  CALL log_warn("init_surf_hgt",                                              &
                "Aggregate surface selected, so using zero height")
  zero_height = .TRUE.
END IF

! If zero height is selected, then that is all we have to do
IF ( zero_height ) THEN
  CALL log_info("init_surf_hgt",                                              &
                "Zero height selected - setting all heights to 0.0")
  jules_vars_data%surf_hgt_surft(:,:) = 0.0
  RETURN
END IF

!-----------------------------------------------------------------------------
! All tiles have heights relative to the gridbox mean
!-----------------------------------------------------------------------------
IF ( .NOT. ANY(l_elev_absolute_height) ) THEN

  IF ( use_file ) THEN
    !-----------------------------------------------------------------------------
    ! If we have a grid, set heights from the specified file
    !-----------------------------------------------------------------------------
    CALL log_info("init_surf_hgt",                                            &
                  "Data is on a grid - reading surf_hgt from file " //        &
                  TRIM(FILE))

    !   Check that a file name was provided
    IF ( LEN_TRIM(FILE) == 0 )                                                &
      CALL log_fatal("init_surf_hgt", "No file name provided")

    CALL fill_variables_from_file(FILE, (/ 'surf_hgt' /), (/ surf_hgt_name /))

  ELSE
    !-----------------------------------------------------------------------------
    ! If we are reading data at a single point, read height from the namelist
    !-----------------------------------------------------------------------------
    ! Copy their values into the model arrays
    CALL log_info("init_surf_hgt",                                            &
                  "Data is at a single point - reading surf_hgt from " //     &
                  "surf_hgt_io in namelist JULES_SURF_HGT")
    jules_vars_data%surf_hgt_surft(1,:) = surf_hgt_io(1:nsurft)

  END IF
END IF

RETURN

END SUBROUTINE init_surf_hgt
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_z_land(ainfo_data,jules_vars_data)

USE io_constants, ONLY: max_sdf_name_len, max_file_name_len, namelist_unit

USE string_utils_mod, ONLY: to_string

USE input_mod, ONLY: fill_variables_from_file

USE missing_data_mod, ONLY: rmdi

USE c_elevate, ONLY: l_elev_absolute_height, z_land_io, surf_hgt_band

USE ancil_info, ONLY: nsurft, land_pts

USE theta_field_sizes, ONLY: t_i_length

USE errormessagelength_mod, ONLY: errormessagelength

!TYPE definitions
USE ancil_info,    ONLY: ainfo_data_type
USE jules_vars_mod, ONLY: jules_vars_data_type

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises the elevation of the forcing data elevation if any tile has
!   an absolute height above sea-level.
!   Tile heights are set to be spatially invarient.
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Arguments
TYPE(ainfo_data_type), INTENT(IN OUT) :: ainfo_data
TYPE(jules_vars_data_type), INTENT(IN OUT) :: jules_vars_data

! Work variables
INTEGER :: error  ! Error indicator

!-----------------------------------------------------------------------------
! Definition of the jules_z_land namelist
!-----------------------------------------------------------------------------
LOGICAL :: use_file
                      !   T - the variable uses the file
                      !   F - the variable is set using a constant value

CHARACTER(LEN=max_file_name_len) :: FILE  ! If input grid has more than one
                                          ! point, read heights from this
                                          ! file
CHARACTER(LEN=max_sdf_name_len) :: z_land_name    ! Use the variable with
                                                  ! this name for forcing
                                                  ! altitude
CHARACTER(LEN=errormessagelength) :: iomessage

NAMELIST  / jules_z_land/ use_file, FILE, z_land_io, z_land_name,             &
                        surf_hgt_band

INTEGER :: i, j, l, n

!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
use_file = .TRUE.      ! Default is for every variable to be read from file
FILE=''       ! Empty file name
z_land_name = 'z_land' ! Default variable name

!-----------------------------------------------------------------------------
! If some tiles have absolute heights then the gridbox mean height must also
! be provided. For a single point, this can come from the namelist, otherwise
! we need to use the file option
!-----------------------------------------------------------------------------
jules_vars_data%z_land_ij(:,:) = 0.0

IF ( ANY(l_elev_absolute_height) ) THEN

  !-----------------------------------------------------------------------------
  !   Read namelist
  !-----------------------------------------------------------------------------
  CALL log_info("init_z_land", "Reading JULES_Z_LAND namelist...")

  !   First, we read the namelist
  READ(namelist_unit, NML = jules_z_land, IOSTAT = error, IOMSG = iomessage)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_z_land",                                             &
                 "Error reading namelist JULES_Z_LAND " //                    &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

  !-----------------------------------------------------------------------------
  !   Set values derived from namelist and verify for consistency
  !-----------------------------------------------------------------------------
  CALL log_info("init_z_land","Some tiles have  " //                          &
         "absolute heights above sea-level, see l_elev_absolute_height " //   &
         "for which. Where l_elev_absolute_height is false, surf_hgt "   //   &
         "offsets can only be applied as global to that tile "           //   &
         "type (usually these are 0, indicating no offset from "         //   &
         "the gridbox mean).")

  !-----------------------------------------------------------------------------
  !   Set the heights (relative or absolute) to be constant across a domain.
  !-----------------------------------------------------------------------------
  DO n = 1,nsurft
    jules_vars_data%surf_hgt_surft(:,n) = surf_hgt_band(n)
  END DO

  !-----------------------------------------------------------------------------
  !   If we have a grid, set gridbox mean heights from the specified file
  !-----------------------------------------------------------------------------
  IF ( use_file ) THEN

    CALL log_info("init_z_land",                                              &
                  "Data is on a grid - reading z_land from file " //          &
                  TRIM(FILE))

    IF ( error > 0 )                                                          &
      CALL log_fatal("init_z_land", "Error allocating z_land")

    !     Check that a file name was provided
    IF ( LEN_TRIM(FILE) == 0 )                                                &
      CALL log_fatal("init_z_land", "No file name provided for gridbox " //   &
                     "mean heights")

    CALL fill_variables_from_file(FILE,                                       &
                                  (/ 'z_land_land' /), (/ z_land_name /))

    DO l = 1,land_pts
      j = ( ainfo_data%land_index(l) - 1 ) / t_i_length + 1
      i = ainfo_data%land_index(l) - (j-1) * t_i_length
      jules_vars_data%z_land_ij(i,j) = jules_vars_data%z_land_land(l)
    END DO

    !-----------------------------------------------------------------------------
    !     If we are reading data at a single point, read height from the namelist
    !-----------------------------------------------------------------------------
  ELSE

    CALL log_info("init_z_land",                                              &
                  "Data is at a single point - reading z_land from " //       &
                  "z_land_io in namelist JULES_Z_LAND")
    jules_vars_data%z_land_ij(1,1) = z_land_io

    !-----------------------------------------------------------------------------
    !     Check that a value for the gridbox height has been set
    !-----------------------------------------------------------------------------
    IF ( z_land_io == rmdi )                                                  &
      CALL log_fatal("init_z_land", "Some tiles have absolute "   //          &
                     "heights above sea-level but no value for "  //          &
                     "z_land has been provided. Set a value for " //          &
                     "z_land in the JULES_Z_LAND namelist")

    !-----------------------------------------------------------------------------
    !     Check that a value for the elevation bands have been set
    !-----------------------------------------------------------------------------
    IF ( ANY(surf_hgt_band(1:nsurft)  == rmdi)  )                             &
      CALL log_fatal("init_z_land", "Some tiles have absolute "   //          &
                     "heights above sea-level but no values for " //          &
                     "elevation bands have been provided. Set a " //          &
                     "value for surf_hgt_band in the "//                      &
                     "JULES_Z_LAND namelist")

  END IF

END IF

RETURN

END SUBROUTINE init_z_land


END MODULE init_grid_mod

