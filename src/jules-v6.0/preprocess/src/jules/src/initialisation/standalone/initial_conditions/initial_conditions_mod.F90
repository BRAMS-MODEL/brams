! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE initial_conditions_mod

USE logging_mod, ONLY: log_info, log_warn, log_error, log_fatal

IMPLICIT NONE

PRIVATE
PUBLIC init_ic

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


SUBROUTINE init_ic(nml_dir,crop_vars,psparms,toppdm,ainfo,trif_vars,          &
                   soilecosse,urban_param,progs,jules_vars)

  !Use in relevant subroutines
USE freeze_soil_mod, ONLY: freeze_soil
USE flake_init_mod,  ONLY: flake_init
USE tilepts_mod,     ONLY: tilepts

USE missing_data_mod, ONLY:                                                   &
!  imported scalar parameters
     rmdi

USE mpi, ONLY: mpi_integer, mpi_comm_world

USE dictionary_mod, ONLY: dict, dict_get, dict_free

USE io_constants, ONLY: max_sdf_name_len, max_file_name_len, namelist_unit

USE conversions_mod, ONLY: zerodegc

USE water_constants_mod, ONLY: rho_water, dpsidt

USE string_utils_mod, ONLY: to_string

USE veg_param, ONLY: litc_norm

USE dump_mod, ONLY: max_var_dump, required_vars_for_configuration,            &
                     read_dump

USE input_mod, ONLY: fill_variables_from_file

USE templating_mod, ONLY: tpl_has_var_name, tpl_substitute_var

USE model_interface_mod, ONLY: identifier_len, populate_var, get_var_id

USE ecosse_init_mod, ONLY: ecosse_init

USE jules_rivers_mod, ONLY: np_rivers, rivers_sto_rp, rivers_boxareas_rp,     &
                            rivers_sto_per_m2_on_landpts

USE jules_hydrology_mod, ONLY: l_top

USE jules_vegetation_mod, ONLY: l_phenol, l_triffid, frac_min, l_crop,        &
                                 l_recon, l_nitrogen, l_croprotate

USE jules_irrig_mod, ONLY: l_irrig_dmd, l_irrig_limit

USE update_mod, ONLY: l_daily_disagg, precip_disagg_method

USE jules_surface_mod, ONLY: l_urban2t

USE jules_soil_biogeochem_mod, ONLY:                                          &
  ! imported scalar parameters
  soil_model_ecosse, soil_model_rothc,                                        &
  ! imported scalar variables
  bio_hum_cn, l_layeredC, soil_bgc_model, tau_lit

USE jules_soil_ecosse_mod, ONLY:                                              &
  ! imported array variables (IN)
  dz_soilc

USE jules_surface_mod, ONLY: l_elev_land_ice, l_flake_model, l_aggregate

USE jules_surface_types_mod, ONLY: npft, ncpft, nnpft,                        &
                                    ice, urban_canyon, urban_roof

USE ancil_info, ONLY: land_pts, lice_pts, soil_pts, surft_pts, nsoilt,        &
                      soilt_pts, dim_cslayer

USE prognostics, ONLY: l_broadcast_soilt_in_mod => l_broadcast_soilt

USE crop_vars_mod, ONLY: startyr, startmon, startday,                         &
                         starttime

USE crop_utils_mod, ONLY: croplai_min, cropcanht_min, croprootc_min

USE pftparm, ONLY: rootd_ft

USE jules_soil_mod, ONLY: dzsoil, sm_levels, l_tile_soil

USE update_mod, ONLY: assign_irrig_fraction

USE time_info_mod, ONLY: current_model_time

USE root_frac_mod, ONLY: root_frac

USE rivers_route_mod, ONLY: scatter_land_from_riv_field

USE soil_biogeochem_control_mod, ONLY:                                        &
  ! imported scalar parameters
  call_initial,                                                               &
  ! imported procedures
  increment_soil_drivers

USE parallel_mod, ONLY: master_task_id, is_master_task

USE errormessagelength_mod, ONLY: errormessagelength

USE um_types, ONLY: real_jlslsm

!TYPE definitions
USE prognostics, ONLY: progs_type
USE crop_vars_mod, ONLY: crop_vars_type
USE p_s_parms,  ONLY: psparms_type
USE top_pdm, ONLY: top_pdm_type
USE ancil_info, ONLY: ainfo_type
USE trif_vars_mod, ONLY: trif_vars_type
USE soil_ecosse_vars_mod, ONLY: soil_ecosse_vars_type
USE urban_param_mod, ONLY: urban_param_type
USE jules_vars_mod, ONLY: jules_vars_type
USE mem_brams_jules, ONLY: runtypeB, mynumB, dir_run_idB, dump_periodB,hfilinB  !DSM


IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Sets up the initial conditions for the run
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
TYPE(psparms_type), INTENT(IN OUT) :: psparms
TYPE(top_pdm_type), INTENT(IN OUT) :: toppdm
TYPE(ainfo_type), INTENT(IN OUT) :: ainfo
TYPE(trif_vars_type), INTENT(IN OUT) :: trif_vars
TYPE(soil_ecosse_vars_type), INTENT(IN OUT) :: soilecosse
TYPE(urban_param_type), INTENT(IN OUT) :: urban_param
TYPE(progs_type), INTENT(IN OUT) :: progs
TYPE(jules_vars_type), INTENT(IN OUT) :: jules_vars

!-----------------------------------------------------------------------------
! Local scalar parameters.
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER ::                                                &
     RoutineName = 'init_ic'   ! Name of this routine.

! Work variables
INTEGER :: nvars_required     ! The number of variables that are
                              ! required in this configuration

CHARACTER(LEN=identifier_len) :: required_vars(max_var_dump)
                               ! The variable identifiers of the required
                               ! variables

INTEGER :: nvars_from_ancil
CHARACTER(LEN=identifier_len) :: vars_from_ancil(max_var_dump)
                               ! The variable identifiers of the ancil
                               ! variables

TYPE(dict) :: default_values  ! Dictionary mapping identifier => default value
                              ! Values should be REAL
REAL(KIND=real_jlslsm) :: default_value
                       ! Variable to contain values read from the dict

INTEGER :: nvars_file   ! The number of variables that will be set
                        ! from the given file (template?)

LOGICAL :: reset_done  ! Indicates if a reset of frac to frac_min was
                       ! performed

LOGICAL :: firstcall = .TRUE.

REAL(KIND=real_jlslsm) :: urban_fraction
                        ! Used in partitioning of urban fraction into
                        ! canyon and roof

REAL(KIND=real_jlslsm) :: f_root_pft(sm_levels)
                              ! Root fraction in each soil layer
                              ! used to initialise plant available inorg N

REAL(KIND=real_jlslsm) :: f_root_pft_dz(sm_levels)
                              ! Normalised roots in each soil layer
                              ! used to initialise plant available inorg N

REAL(KIND=real_jlslsm) :: rivers_sto_per_m2_rgrid(np_rivers)
          ! water storage on routing grid in kg m-2

REAL ::                                                                       &
 qbase_l_soilt_dummy(land_pts,nsoilt,sm_levels+1),                            &
    ! Lateral flux of water from each soil layer (kg m-2 s-1).
    ! Used as a dummy argument.
 w_flux_soilt_dummy(land_pts,nsoilt,0:sm_levels)
    ! Downward water flux at bottom of each soil layer (kg m-2 s-1).
    ! Used as a dummy argument.

INTEGER :: nstep_dummy        ! Used as a dummy argument.
LOGICAL :: l_run_model_dummy  ! Used as a dummy argument.

INTEGER :: i,j,l,n,m,ip  ! Loop counters

INTEGER :: n_counted ! Counter for assigning soil_pts

INTEGER :: year, month, day, time

INTEGER :: error  ! Error indicator
CHARACTER(LEN=errormessagelength) :: iomessage

CHARACTER(LEN=256) :: aux,aux2   !DSM
INTEGER ::  hh, mm

!-----------------------------------------------------------------------------
! Definition of the jules_initial namelist
!-----------------------------------------------------------------------------
LOGICAL :: dump_file  ! T - the given file is a dump file
                      ! F - the given file is not a dump file

LOGICAL :: total_snow
                      ! Switch indicating how the snow model is initialised
                      !   T - only snow_tile is needed
                      !       If nsmax>0, the layer values are determined
                      !       by the model
                      !   F - all snow variables are supplied directly (what
                      !       these are depends on nsmax)

CHARACTER(LEN=max_file_name_len) :: FILE
                      ! The name of the file (or variable name template) to
                      ! use for variables that need to be filled from file

INTEGER :: nvars      ! The number of variables in this section
CHARACTER(LEN=identifier_len) :: var(max_var_dump)
                      ! The variable identifiers of the variables
LOGICAL :: use_file(max_var_dump)
                      !   T - the variable uses the file
                      !   F - the variable is set using a constant value
                      ! Defaults to T for every variable
CHARACTER(LEN=max_sdf_name_len) :: var_name(max_var_dump)
                      ! The name of each variable in the file
CHARACTER(LEN=max_sdf_name_len) :: tpl_name(max_var_dump)
                      ! The name to substitute in a template for each
                      ! variable
REAL(KIND=real_jlslsm) :: const_val(max_var_dump)
                      ! The constant value to use for each variable if
                      ! use_file = F for that variable
LOGICAL :: in_init_nml, loaded_from_ancil
                      !Used for checking whether all req vars have been
                      !accounted for

LOGICAL :: l_broadcast_soilt = .FALSE.
                      !IO local variable.
                      !Value then held in prognostics module.
                      !Switch to broadcast model state around all soil tiles.
                      !Only has an effect if l_tile_soil is true.
                      !Does not do anything with ancils- this is done by
                      !l_broadcast_ancils

LOGICAL, PARAMETER :: l_output_mode = .FALSE.

NAMELIST  / jules_initial/ total_snow, dump_file, FILE,                       &
                         nvars, var, use_file, var_name, tpl_name,            &
                         const_val, l_broadcast_soilt

!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
nvars_required = 0
nvars_file     = 0
dump_file      = .FALSE.
total_snow     = .FALSE.
nvars          = 0
use_file(:)    = .TRUE.  ! Default is to set every variable from file
FILE=''      ! Empty file name.
var(:)         = ''      ! Empty identifiers.
var_name(:)    = ''      ! Empty variable names.
tpl_name(:)    = ''      ! Empty template names.
const_val(:)   = rmdi    ! Missing data value. This might later be replaced
                         ! by a default value from a dictionary.

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
CALL log_info(RoutineName, "Reading JULES_INITIAL namelist...")

OPEN(namelist_unit, FILE=(TRIM(nml_dir) // '/' // 'initial_conditions.nml'),  &
               STATUS='old', POSITION='rewind', ACTION='read', IOSTAT = error,&
               IOMSG = iomessage)
IF ( error /= 0 ) THEN
  CALL log_fatal(RoutineName,                                                 &
                 "Error opening namelist file initial_conditions.nml " //     &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")
END IF

READ(namelist_unit, NML = jules_initial, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 ) THEN
  CALL log_fatal(RoutineName,                                                 &
                 "Error reading namelist JULES_INITIAL " //                   &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")
END IF

CLOSE(namelist_unit, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 ) THEN
  CALL log_fatal(RoutineName,                                                 &
                 "Error closing namelist file initial_conditions.nml " //     &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")
END IF

!DSM{
if (trim(runtypeB)=='INITIAL') then
   dump_file=.false.
   file='./dummy.dat'
elseif (trim(runtypeB)=='HISTORY') then
   dump_file=.true.
   nvars=21
   use_file(:)=.true.
   aux=dir_run_idB(1:index(dir_run_idB,'/',BACK = .TRUE.)-1)
   aux2=trim(dir_run_idB(index(dir_run_idB,'/',BACK = .TRUE.)+1:180))
   write(aux2,'(a,i5.5)') trim(aux2),mynumB
   file=trim(aux)//'/'//trim(aux2)//'.dump.'
   
   aux=trim(hfilinB(index(hfilinB,'-head.txt',BACK = .TRUE.)-17:index(hfilinB,'-head.txt',BACK = .TRUE.)-1))
   
   read(aux(12:13),*) hh
   read(aux(14:15),*) mm
   write(aux2,*) hh*3600+mm*60
   aux2=adjustl(aux2)
   
   file=trim(file)//aux(1:4)//aux(6:7)//aux(9:10)//'.'//trim(aux2)//'.nc'
else
   print*, 'ERROR. The code is not prepared for RUNTYPE='//trim(runtypeB)
   stop
endif
!DSM}

!-----------------------------------------------------------------------------
! Check that variable identifiers are not empty.
! Although we might later decide that the identifier is not required, for
! clarity we check here whether the claimed amount of information was
! provided.
!-----------------------------------------------------------------------------
DO i = 1,nvars
  IF ( LEN_TRIM(var(i)) == 0 ) THEN
    CALL log_fatal(RoutineName,                                               &
                   "Insufficient values for var. " //                         &
                   "No name provided for var at position #" //                &
                   TRIM(to_string(i)) )
  END IF
END DO

!Copy across _io variable to its proper place
l_broadcast_soilt_in_mod = l_broadcast_soilt

!-----------------------------------------------------------------------------
! Check switches against soil tiling
!-----------------------------------------------------------------------------
IF ( l_broadcast_soilt ) THEN
  IF ( l_tile_soil ) THEN
    CALL log_info(RoutineName,                                                &
                  "Non-soil tiled initial conditions will be broadcast " //   &
                  "to all soil tiles. Users should consider appropriate" //   &
                  " spinup period")
  ELSE
    CALL log_warn(RoutineName,                                                &
                  "l_broadcast_soilt will have no effect: l_tile_soilt = F")
  END IF
END IF

!-----------------------------------------------------------------------------
! Set up initial conditions using namelist values
!-----------------------------------------------------------------------------

! Set up the required variables - we get the list by calling a procedure in
! dump_mod that tells us the required prognostic variables for the current
! model configuration
! We indicate whether or not we are using total_snow, as this affects the
! required variables
! We also indicate that we do not want IMOGEN prognostics in the list, even
! if they are required, since they are initialised in init_imogen
! We do request the ancillaries
CALL required_vars_for_configuration( nvars_required, required_vars,          &
                                      nvars_from_ancil, vars_from_ancil,      &
                                      l_output_mode,                          &
                                      total_snow, .FALSE., dump_file)

!Get a dictionary of default values in case they are needed
default_values = get_default_ic_values(total_snow)

! If we are initialising from a dump and no variables were specified, then
! we assume that all variables will be initialised from the dump file
IF ( dump_file .AND. nvars < 1 ) THEN
  CALL log_info(RoutineName,                                                  &
                "No variables given - will attempt to initialise all " //     &
                "required variables from specified dump file")
  nvars = nvars_required
  var(1:nvars) = required_vars(1:nvars)
  ! Every variable will use the file
  use_file(:) = .TRUE.
  ! We don't need to set var_name, tpl_name or const_val since they are never
  ! used in the case of a dump file anyway

ELSE
  !Check that all required variables are accounted for and fill in with
  !default values as needed
  !** NB- Default values are not available when using a dump file **
  DO i = 1, nvars_required

    !Test if the required variable is in the initial conditions namelist
    IF ( ANY(var(1:nvars) == required_vars(i)) ) THEN
      in_init_nml = .TRUE.
    ELSE
      in_init_nml = .FALSE.
    END IF

    !Test if the required variable has been loaded from an ancil
    IF ( ANY(vars_from_ancil(1:nvars_from_ancil) == required_vars(i)) ) THEN
      loaded_from_ancil = .TRUE.
    ELSE
      loaded_from_ancil = .FALSE.
    END IF

    !Log our assessment of whether the variable is accounted for
    IF ( in_init_nml .OR. loaded_from_ancil ) THEN
      CALL log_info(RoutineName,                                              &
           "'" // TRIM(required_vars(i)) // "' accounted for- OK")
    ELSE
      CALL log_info(RoutineName,                                              &
           "'" // TRIM(required_vars(i)) // "' not accounted for- " //        &
           "will attempt to set to default value")

      !Add the req variable to the var array along with with its default value
      nvars            = nvars + 1
      var(nvars)       = required_vars(i)
      use_file(nvars)  = .FALSE.
      CALL dict_get(default_values, required_vars(i), default_value)
      const_val(nvars) = default_value

    END IF !in_init_nml .OR. loaded_from_ancil
  END DO !nvars_required
END IF !dump_file .AND. nvars < 1

!-----------------------------------------------------------------------------
! Check which variables we will be using and partition them into variables
! set to constant values and variables set from file
!-----------------------------------------------------------------------------
DO i = 1,nvars
  !---------------------------------------------------------------------------
  ! If the variable is one of the required vars, then we will be using it
  !---------------------------------------------------------------------------
  IF ( ANY(required_vars(1:nvars_required) == var(i)) ) THEN
    IF ( use_file(i) ) THEN
      CALL log_info(RoutineName,                                              &
                    "'" // TRIM(var(i)) // "' will be read from file")

      ! If the variable will be filled from file, register it here
      nvars_file = nvars_file + 1
      ! Since nvars_file <= i (so we will not overwrite unprocessed values)
      ! and we do not need the values from these arrays for any non-file
      ! variables from now on, we can just compress them down onto variables
      ! that are in the file.
      var(nvars_file) = var(i)
      var_name(nvars_file) = var_name(i)
      tpl_name(nvars_file) = tpl_name(i)
    ELSE
      ! If the variable is being set as a constant, populate it here.
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
    END IF  !  use_file
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
  !   Check that a file name was provided.
  IF ( LEN_TRIM(FILE) == 0 ) THEN
    CALL log_fatal(RoutineName, "No file name provided")
  END IF

  IF ( dump_file ) THEN
    ! If we are using a dump file, use read_dump to fill the variables
    CALL read_dump(FILE, var(1:nvars_file))
  ELSE IF ( tpl_has_var_name(FILE) ) THEN
    ! If we are using a non-dump file with a variable name template, loop
    ! through the variables setting one from each file.
    DO i = 1,nvars_file
      ! Check that a template string was provided for this variable.
      IF ( LEN_TRIM(tpl_name(i)) == 0 ) THEN
        CALL log_fatal( RoutineName,                                          &
                        "No variable name template substitution " //          &
                        "provided for " // TRIM(var(i)) )
      END IF
      CALL fill_variables_from_file(                                          &
        tpl_substitute_var(FILE, tpl_name(i)),                                &
        (/ var(i) /), (/ var_name(i) /)                                       &
      )
    END DO
  ELSE
    ! We are not using a file name template, so set all variables from the
    ! same file.
    CALL fill_variables_from_file(                                            &
      FILE, var(1:nvars_file), var_name(1:nvars_file)                         &
    )
  END IF  !  dump_file
END IF  !  nvars_file

! Free the dictionary of default values as it is no longer required
CALL dict_free(default_values)

!*****************************************************************************
! Further processing depending on options specified
!*****************************************************************************

!-----------------------------------------------------------------------------
! Set up derived soil values
!-----------------------------------------------------------------------------
! Set surface values.
psparms%hcon_soilt(:,:,0)   = psparms%hcon_soilt(:,:,1)
psparms%satcon_soilt(:,:,0) = psparms%satcon_soilt(:,:,1)

! Check that psparms%sathh_soilt>=0 - a common error!
IF ( ANY( psparms%sathh_soilt(:,:,:) < 0.0 ) )                                &
  CALL log_fatal("init_soil",                                                 &
                 "sathh < 0.0 detected - for JULES, sathh is " //             &
                 "abs(saturated head)")

! Check that clay_soilt>=0 and <=1 - a common error!
SELECT CASE ( soil_bgc_model )
CASE ( soil_model_ecosse, soil_model_rothc )
  IF ( ANY( psparms%clay_soilt(:,:,:) < 0.0 ) .OR.                            &
       ANY( psparms%clay_soilt(:,:,:) > 1.0 ) ) THEN
    CALL log_fatal("init_soil", "clay <0 and/or >1.0 detected - for JULES")
  END IF
END SELECT

!-----------------------------------------------------------------------------
! Detect soil points.
! If top layer saturation moisture content > 0, this is a soil point.
! Note that land ice points are no longer assigned here.
!-----------------------------------------------------------------------------
soil_pts        = 0
ainfo%soil_index(:)   = 0
ainfo%l_soil_point(:) = .FALSE.
DO l = 1,land_pts
  IF (nsoilt == 1) THEN
    m = 1
    !Use existing logic
    IF ( psparms%smvcst_soilt(l,m,1) >                                        &
         EPSILON(psparms%smvcst_soilt(l,m,1)) ) THEN
      soil_pts             = soil_pts + 1
      ainfo%soil_index(soil_pts) = l
      ainfo%l_soil_point(l)      = .TRUE.
    END IF
  ELSE
    !Use soil-tile aware logic

    !Test every soil tile for the smvcst condition
    n_counted = 0
    DO m = 1,nsoilt
      IF (psparms%smvcst_soilt(l,m,1) >                                       &
          EPSILON(psparms%smvcst_soilt(l,m,1))) THEN
        n_counted = n_counted + 1
      END IF
    END DO

    !If there are any tiles that pass the above, make sure it is all of them
    IF (n_counted == nsoilt) THEN
      soil_pts             = soil_pts + 1
      ainfo%soil_index(soil_pts) = l
      ainfo%l_soil_point(l)      = .TRUE.
    ELSE IF (n_counted > 0 .AND. n_counted /= nsoilt) THEN
      CALL log_fatal(RoutineName,                                             &
                     "Gridboxes must be entirely soil or land ice")
      !ELSE
        !It's not a soil point. Do nothing
    END IF
  END IF !nsoilt == 1
END DO !land_pts

CALL log_info(RoutineName,                                                    &
              "Number of soil points = " // TRIM(to_string(soil_pts)))

IF ( soil_pts <= 0 ) THEN
  CALL log_warn(RoutineName,                                                  &
                "There are no soil points - any land points are land ice")
END IF

!Set up up soilt tile points, index and frac
IF ( nsoilt == 1 ) THEN
  m = 1
  !Direct mapping of soil tile points to soil points
  soilt_pts(1)     = soil_pts
  ainfo%soilt_index(:,m) = ainfo%soil_index

  !All fractions must equal 1 by definition
  ainfo%frac_soilt(:,m)  = 1.0

ELSE !nsoilt == nsurft

  DO m = 1,nsoilt

    !Set to zero initially to ensure we don't make a mess of gbm calculations
    soilt_pts(m) = 0
    DO l = 1, land_pts
      ainfo%frac_soilt(l,m)  = 0.0
    END DO

    !Now work through each soil point, updating the number of points, index
    !and fraction as appropriate.
    !Assumes all non-zero fractions are non-negligible
    DO l = 1, soil_pts
      i = ainfo%soil_index(l)
      IF (ainfo%frac_surft(i,m) > 0.0) THEN
        soilt_pts(m)     = soilt_pts(m) + 1
        ainfo%soilt_index(i,m) = i
        ainfo%frac_soilt(i,m)  = ainfo%frac_surft(i,m)
      END IF
    END DO

  END DO
END IF


! This accounts for some of the water that remains unfrozen when the soil is
! very cold
DO i = 1,soil_pts
  l = ainfo%soil_index(i)
  psparms%sthu_min_soilt(l,:,:) =                                             &
    ( dpsidt * zerodegc / psparms%sathh_soilt(l,:,:) )**                      &
    ( -1.0 / psparms%bexp_soilt(l,:,:) )
END DO

! Calculate soil moisture content from wetness
DO i = 1,sm_levels
  progs%smcl_soilt(:,:,i) = rho_water * dzsoil(i) *                           &
    jules_vars%sthuf_soilt(:,:,i) * psparms%smvcst_soilt(:,:,i)
END DO

!-----------------------------------------------------------------------------
! Set the normalisation factor for vertical profile of litter input.
!-----------------------------------------------------------------------------
IF ( l_layeredc .AND. soil_bgc_model == soil_model_rothc ) THEN
  ! Calculate dz * exp(-tau*z) for mid-point depth in each layer.
  litc_norm =  dzsoil(1) * EXP( -tau_lit * 0.5 * dzsoil(1) )
  DO j = 2,sm_levels
    litc_norm = litc_norm + dzsoil(j) *                                       &
                EXP( -tau_lit * ( SUM(dzsoil(1:j-1)) + 0.5 * dzsoil(j) ) )
  END DO
ELSE IF ( soil_bgc_model == soil_model_ecosse ) THEN
  ! Calculate dz * exp(-tau*z) for mid-point depth in each layer.
  litc_norm = dz_soilc(1) * EXP(-tau_lit * dz_soilc(1) * 0.5)
  DO j = 2,dim_cslayer
    litc_norm = litc_norm + dz_soilc(j) *                                     &
                EXP( -tau_lit *                                               &
                ( SUM(dz_soilc(1:j-1)) + 0.5 * dz_soilc(j) ) )
  END DO
END IF

!-----------------------------------------------------------------------------
! If using the two-tile urban schemes and only a combined urban fraction is
! given then split the fraction between the canyon and roof. This has to be
! done here instead of in init_urban to be consistent with triffid.
!-----------------------------------------------------------------------------
IF ( l_urban2t ) THEN
  CALL log_info(RoutineName,                                                  &
                "Either URBAN-2T or MORUSES is in use - splitting urban " //  &
                "tile into canyon/roof if needed")

  DO l = 1, land_pts
    IF ( ainfo%frac_surft(l,urban_canyon) > 0.0 .AND.                         &
         ainfo%frac_surft(l,urban_roof) == 0.0 ) THEN
      urban_fraction             = ainfo%frac_surft(l,urban_canyon)
      ainfo%frac_surft(l,urban_canyon) = urban_fraction * urban_param%wrr_gb(l)
      ainfo%frac_surft(l,urban_roof)   = urban_fraction -                     &
                                   ainfo%frac_surft(l,urban_canyon)
    ELSE IF ( ainfo%frac_surft(l,urban_canyon) > 0.0 .AND.                    &
              ainfo%frac_surft(l,urban_roof) > 0.0 .AND.                      &
              firstcall ) THEN
      CALL log_warn(RoutineName,                                              &
         "WARNING: URBAN-2T/M dump being used to initialise? " //             &
         " Splitting not done: Roof fraction already exists")
      firstcall = .FALSE.
    END IF
  END DO
END IF

!-----------------------------------------------------------------------------
!   If using TRIFFID (with or without competing veg), ensure that fractions of
!   PFTs are not below minimum. Only do this over soil points - land ice
!   points should have zero fractions.
!-----------------------------------------------------------------------------
! Set up a flag to see if any points were reset to frac_min
reset_done = .FALSE.

IF ( l_triffid .AND. l_recon ) THEN
  DO j = 1,soil_pts
    i = ainfo%soil_index(j)
    IF ( ANY( ainfo%frac_surft(i,:) < frac_min ) ) THEN
      ! Reset all small values. Renormalisation is done later, but will fail
      ! if frac_min is sufficiently large. We only reset natural PFT tiles.
      WHERE ( ainfo%frac_surft(i,1:nnpft) < frac_min )
        ainfo%frac_surft(i,1:nnpft) = frac_min
      END WHERE
      reset_done = .TRUE.
    END IF
  END DO

  IF ( reset_done ) THEN
    CALL log_warn(RoutineName,                                                &
                  "frac < frac_min at one or more points - reset to " //      &
                  "frac_min at those points")
  END IF
END IF

!-----------------------------------------------------------------------------
!   If using l_croprotate ensure that fractions of
!   crop PFTs are not below minimum. Only do this over soil points- land ice
!   points should have zero fractions.
!   This ensures that each crop fraction is initialised properly so it can be
!   used later. If a fraction is zero at the start it cannot be used later.
!-----------------------------------------------------------------------------

IF (l_croprotate) THEN
  DO j = 1, soil_pts
    i = ainfo%soil_index(j)
    IF ( ANY( ainfo%frac_surft(i,nnpft+1:nnpft + ncpft) < frac_min ) ) THEN
      ! Reset all small values. Renormalisation is done later, but will fail
      ! if frac_min is sufficiently large. We only reset crop PFT tiles.
      DO n = 1,ncpft
        WHERE ( ainfo%frac_surft(:,n + nnpft) < frac_min )
          ainfo%frac_surft(:,n + nnpft) = frac_min
        END WHERE
      END DO
      reset_done = .TRUE.
    END IF
  END DO

  IF ( reset_done ) THEN
    CALL log_warn(RoutineName,                                                &
                "frac < frac_min at one or more points - reset to " //        &
                "frac_min at those points")
  END IF
END IF

!-----------------------------------------------------------------------------
! Initialise soil bio and hum pools from soil carbon.
!-----------------------------------------------------------------------------
IF ( soil_bgc_model == soil_model_rothc ) THEN
  !Triffid is incompatible with soil tiling at present. Hard code soil tile
  !index to 1 using m = 1
  m = 1
  DO n = 1,dim_cslayer
    progs%ns_pool_gb(:,n,3) = progs%cs_pool_soilt(:,m,n,3) / bio_hum_cn
    progs%ns_pool_gb(:,n,4) = progs%cs_pool_soilt(:,m,n,4) / bio_hum_cn
  END DO
END IF

!-----------------------------------------------------------------------------
! Reconfigure ECOSSE.
! For now we always call this (regardless of l_recon) so as to deal with
! zero (or tiny) pools.
!-----------------------------------------------------------------------------
IF ( soil_bgc_model == soil_model_ecosse ) THEN
  CALL ecosse_init(psparms%soil_pH_soilt,ainfo%frac_surft,progs%cs_pool_soilt)
END IF

!-----------------------------------------------------------------------------
! Initialise soil biogeochemistry drivers.
! At present this need only be done for ECOSSE.
!-----------------------------------------------------------------------------
IF ( soil_bgc_model == soil_model_ecosse ) THEN
  ! Most of these arguments are not used in this call (but must be present).
  ! However under some circumstances the values of sthf_soilt, sthu_soilt and
  ! t_soil_soilt passed here are used to set initial values - so give sensible
  ! values!
  CALL increment_soil_drivers( call_initial, land_pts, nstep_dummy,           &
                               l_run_model_dummy, trif_vars%deposition_n_gb,  &
                               qbase_l_soilt_dummy, psparms%sthf_soilt,       &
                               psparms%sthu_soilt, progs%t_soil_soilt,        &
                               w_flux_soilt_dummy,                            &
                               ainfo%soil_index,                              &
                               ! TYPES
                               soilecosse )
END IF

!-----------------------------------------------------------------------------
! If using crop model, fill lai and canht from croplai_cpft and cropcanht if
! lai and
! canht_pft were not prognostics (whether they're prognostics is specified in
! required_vars_for_configuration)
!-----------------------------------------------------------------------------
IF ( l_crop ) THEN
  IF ( l_phenol ) THEN
    DO n = 1,ncpft
      crop_vars%croplai_cpft(:,n) = progs%lai_pft(:,n + nnpft)
    END DO
  ELSE
    DO n = 1,ncpft
      progs%lai_pft(:,n + nnpft) = crop_vars%croplai_cpft(:,n)
    END DO
  END IF

  IF ( l_triffid ) THEN
    DO n = 1,ncpft
      crop_vars%cropcanht_cpft(:,n) = progs%canht_pft(:,n + nnpft)
    END DO
  ELSE
    DO n = 1,ncpft
      progs%canht_pft(:,n + nnpft) = crop_vars%cropcanht_cpft(:,n)
    END DO
  END IF

  ! ensure that crop leaf area index is not below minimum
  reset_done = .FALSE.
  DO i = 1,land_pts
    IF ( ANY( crop_vars%croplai_cpft(i,:) < croplai_min ) ) THEN
      ! Reset all small values. We only reset crop PFT tiles
      WHERE ( crop_vars%croplai_cpft(i,:) < croplai_min )
        crop_vars%croplai_cpft(i,:) = croplai_min
      END WHERE
      reset_done = .TRUE.
    END IF
  END DO

  IF ( reset_done ) THEN
    CALL log_warn(RoutineName,                                                &
                "croplai_cpft < croplai_min at one or more points - " //      &
                "reset to croplai_min at those points")
    DO n = 1,ncpft
      progs%lai_pft(:,n + nnpft) = crop_vars%croplai_cpft(:,n)
    END DO
  END IF

  ! ensure that crop canopy height is not below minimum
  reset_done = .FALSE.
  DO i = 1,land_pts
    IF ( ANY( crop_vars%cropcanht_cpft(i,:) < cropcanht_min ) ) THEN
      ! Reset all small values. We only reset crop PFT tiles
      WHERE ( crop_vars%cropcanht_cpft(i,:) < cropcanht_min )
        crop_vars%cropcanht_cpft(i,:) = cropcanht_min
      END WHERE
      reset_done = .TRUE.
    END IF
  END DO

  IF ( reset_done ) THEN
    CALL log_warn(RoutineName,                                                &
                "cropcanht < cropcanht_min at one or more points - " //       &
                "reset to cropcanht_min at those points")
    DO n = 1,ncpft
      progs%canht_pft(:,n + nnpft) = crop_vars%cropcanht_cpft(:,n)
    END DO
  END IF

  ! ensure that crop root carbon is not below minimum
  reset_done = .FALSE.
  DO i = 1,land_pts
    IF ( ANY( crop_vars%rootc_cpft(i,:) < croprootc_min ) ) THEN
      ! Reset all small values.
      WHERE ( crop_vars%rootc_cpft(i,:) < croprootc_min )
        crop_vars%rootc_cpft(i,:) = croprootc_min
      END WHERE
      reset_done = .TRUE.
    END IF
  END DO

  IF ( reset_done ) THEN
    CALL log_warn(RoutineName,                                                &
                "rootc_cpft < croprootc_min at one or more points - " //      &
                "reset to croprootc_min at those points")
  END IF
END IF

!-----------------------------------------------------------------------------
! Check that frac_surft sums to 1.0 (with a bit of leeway).
!-----------------------------------------------------------------------------
IF ( l_recon .AND. .NOT. l_croprotate) THEN
  DO i = 1,land_pts
    ! If the discrepancy is big enough, bail
    IF ( ABS( SUM(ainfo%frac_surft(i,:)) - 1.0 ) >= 1.0e-2 ) THEN
      CALL log_fatal(RoutineName,                                             &
                     "frac does not sum to 1 at point " //                    &
                     TRIM(to_string(i)) //                                    &
                     " and the discrepancy is too big to be removed")
    ELSE IF ( ABS( SUM(ainfo%frac_surft(i,:)) - 1.0 ) > 1.0e-4 ) THEN
      ! Correct a small discrepancy
      CALL log_warn(RoutineName,                                              &
                    "frac does not sum to 1 at point " //                     &
                     TRIM(to_string(i)) //                                    &
                     " - removing small discrepancy")
      ! Ignore small discrepancies and (re)normalise
      ainfo%frac_surft(i,:) = ainfo%frac_surft(i,:) / SUM(ainfo%frac_surft(i,:))
    END IF

    ! Ignore any discrepancy below the threshold completely
  END DO
END IF

!-----------------------------------------------------------------------------
! Process the ice fraction field.
! Identify land ice points, only if the ice surface type is specified.
!-----------------------------------------------------------------------------
lice_pts        = 0
ainfo%lice_index(:)   = 0
ainfo%l_lice_point(:) = .FALSE.

IF ( ice > 0 .AND. .NOT. l_elev_land_ice ) THEN
  DO l = 1,land_pts
    IF ( ainfo%frac_surft(l,ice) > 0.0 ) THEN
      ! This is a land ice point.
      lice_pts             = lice_pts + 1
      ainfo%lice_index(lice_pts) = l
      ainfo%l_lice_point(l)      = .TRUE.

      ! If reconfiguring, ensure pools are set (mostly to zero) over ice.
      IF ( l_recon ) THEN
        IF ( l_triffid ) THEN
          progs%canht_pft(l,:)            = 0.0
          progs%lai_pft(l,:)              = 0.0
          ainfo%frac_surft(l,1:nnpft)     = 0.0
        END IF
        IF ( soil_bgc_model == soil_model_rothc ) THEN
          progs%n_inorg_gb(l)             = 0.0
          progs%n_inorg_soilt_lyrs(l,:,:) = 0.0
        END IF
        IF ( soil_bgc_model == soil_model_ecosse .OR.                         &
             soil_bgc_model == soil_model_rothc ) THEN
          progs%cs_pool_soilt(l,:,:,:)    = 1.0e-6
          progs%ns_pool_gb(l,:,:)         = 1.0e-6
        END IF
      END IF

      ! At present, land ice and soil points are mutually exclusive.
      ! Check this is not a soil point
      IF ( ANY(ainfo%soil_index == l) ) THEN
        CALL log_fatal(RoutineName,                                           &
                       "Land ice points and soil points are mutually " //     &
                       "exclusive")
      END IF

      ! Check that ice fraction is one (cannot have partial ice coverage).
      IF ( ABS(ainfo%frac_surft(l,ice) - 1.0) > EPSILON(1.0) ) THEN
        CALL log_error(RoutineName,                                           &
                       "Ice fraction must be 1 at an ice point")
      END IF

    END IF  !  land ice points
  END DO

END IF  !  ice > 0 .AND. .NOT.l_elev_land_ice

!-----------------------------------------------------------------------------
! fall back to the old smvcst (NOT soil) criteria when using the tiled land ice
! subsurface - consistency with the UM
!-----------------------------------------------------------------------------
IF ( l_elev_land_ice ) THEN
  !Only need to check the 1st soil tile as any problems will have been
  !picked up by the soil_pts test above. Set m = 1
  m = 1
  DO l = 1,land_pts
    IF ( psparms%smvcst_soilt(l,m,1) <=                                       &
         EPSILON(psparms%smvcst_soilt(l,m,1)) ) THEN
      lice_pts = lice_pts + 1
      ainfo%lice_index(lice_pts) = l
      ainfo%l_lice_point(l)      = .TRUE.
    END IF
    !
    !Reverse of soil criteria, so automatically exclusive.
    !Partial ice fractions (e.g in each elevation class) allowed for land ice
    !points here
    !
  END DO
END IF

CALL log_info(RoutineName,                                                    &
              "Number of land ice points = " // TRIM(to_string(lice_pts)))

!-----------------------------------------------------------------------------
! Check that all land points have been identified as either soil or ice
!-----------------------------------------------------------------------------
IF ( soil_pts + lice_pts /= land_pts ) THEN
  CALL log_fatal(RoutineName,                                                 &
                 "All points should be either soil or land ice points - " //  &
                 "have land_pts = " // TRIM(to_string(land_pts)) //           &
                 " and soil_pts + lice_pts = " //                             &
                 TRIM(to_string(soil_pts + lice_pts)))
END IF

! Additional check using ainfo%l_soil_point and ainfo%l_lice_point
DO l = 1, land_pts
  IF (ainfo%l_soil_point(l) .AND. ainfo%l_lice_point(l)) THEN
    CALL log_fatal(RoutineName,                                               &
                   "Point " // TRIM(to_string(l)) // " has been assigned" //  &
                   " as both a soil and land ice point. Please review the"//  &
                   " consistency for the soil and frac ancils.")
  ELSE IF ( .NOT. ainfo%l_soil_point(l) .AND. .NOT. ainfo%l_lice_point(l)) THEN
    CALL log_fatal(RoutineName,                                               &
                   "Point " // TRIM(to_string(l)) // " has been assigned" //  &
                   " as neither a soil or land ice point. Please review"  //  &
                   " the consistency for the soil and frac ancils.")
  END IF
END DO

!-----------------------------------------------------------------------------
! Set up tile index
!-----------------------------------------------------------------------------
CALL tilepts( land_pts, ainfo%frac_surft, surft_pts, ainfo%surft_index,       &
              ainfo%l_lice_point )

! For URBAN-2T or MORUSES: Check that urban canyons also have roofs
IF ( l_urban2t ) THEN
  IF ( surft_pts(urban_canyon) /= surft_pts(urban_roof) ) THEN
    CALL log_fatal(RoutineName,                                               &
                   "URBAN-2T or MORUSES - # canyons /= # roofs")
  END IF
END IF

!-----------------------------------------------------------------------------
! Deal with "simple" initialisation of snow variables.
!-----------------------------------------------------------------------------
IF ( total_snow ) CALL total_snow_init(ainfo, progs)

!-----------------------------------------------------------------------------
! Calculate frozen and unfrozen fractions of soil moisture.
! freeze_soil assumes bexp_soilt, sathh_soilt, smvcst_soilt are
! constant in soil column, so loop around soil layers to allow
! depth varying soil properties and to maintain compatibility with the UM
!-----------------------------------------------------------------------------
DO i = 1,sm_levels
  DO m = 1,nsoilt
    CALL freeze_soil (land_pts, 1, psparms%bexp_soilt(:,m,i), dzsoil(i:i),    &
                      psparms%sathh_soilt(:,m,i), progs%smcl_soilt(:,m,i),    &
                      progs%t_soil_soilt(:,m,i), psparms%smvcst_soilt(:,m,i), &
                      psparms%sthu_soilt(:,m,i), psparms%sthf_soilt(:,m,i))
  END DO
END DO
!-----------------------------------------------------------------------------
! Finish initialising TOPMODEL
!-----------------------------------------------------------------------------
IF ( l_top ) THEN
  !Calculate fitting parameters
  CALL calc_fit_fsat(toppdm,ainfo)
  CALL topmodel_init(psparms,toppdm,ainfo)
END IF

!-----------------------------------------------------------------------------
! Initialise FLake.
!-----------------------------------------------------------------------------

IF ( l_flake_model                                                            &
   .AND. ( .NOT. l_aggregate)                                                 &
   .AND. (land_pts > 0 ) ) THEN
  CALL flake_init(ainfo, progs)
END IF ! FLake

!---------------------------------------------------------------------------
! Finish initialising irrigation
! Process pft names to be assigned irr fraction (if not all)
!---------------------------------------------------------------------------
CALL current_model_time(year,month,day,time)
startyr   = year
startmon  = month
startday  = day
starttime = time

IF ( l_irrig_dmd ) THEN
  CALL assign_irrig_fraction(crop_vars,ainfo)
  IF ( l_irrig_limit ) THEN
    DO ip = 1, np_rivers
      rivers_sto_per_m2_rgrid(ip) = rivers_sto_rp(ip) / rivers_boxareas_rp(ip)
    END DO
    CALL scatter_land_from_riv_field(rivers_sto_per_m2_rgrid,                 &
                                     rivers_sto_per_m2_on_landpts)
  END IF
END IF
!-----------------------------------------------------------------------------
! If we are using the disaggregator with random rainfall, initialise the
! random seed
!-----------------------------------------------------------------------------
IF ( l_daily_disagg .AND. precip_disagg_method > 1 ) THEN
  IF ( .NOT. dump_file ) THEN
    ! If not using a dump file, we need to initialise a new random seed on the
    ! master task and broadcast it to all other tasks
    IF ( is_master_task() ) THEN
      CALL RANDOM_SEED()
      CALL RANDOM_SEED( get = progs%seed_rain )
    END IF
    CALL mpi_bcast(progs%seed_rain, SIZE(progs%seed_rain), mpi_integer,       &
                   master_task_id, mpi_comm_world, error)
  END IF

  ! Make sure that we put the seed, whether it came from the dump or was
  ! broadcast by the master task
  CALL RANDOM_SEED( put = progs%seed_rain )
END IF

IF (l_nitrogen .AND. l_layeredC) THEN
  ! Initialise n_inorg_avail_pft so it can be leached before the first call
  ! to TRIFFID.
  DO n = 1,npft
    CALL root_frac(n,sm_levels,dzsoil,rootd_ft(n),f_root_pft)
    f_root_pft_dz = f_root_pft / dzsoil / (f_root_pft(1) / dzsoil(1))
    DO i = 1,soil_pts
      l = ainfo%soil_index(i)
      DO j = 1,dim_cslayer
        IF (nsoilt == 1) THEN
          m = 1
          progs%n_inorg_avail_pft(l,n,j) = progs%n_inorg_soilt_lyrs(l,m,j) *  &
                                     f_root_pft_dz(j)
        ELSE
          progs%n_inorg_avail_pft(l,n,j) = progs%n_inorg_soilt_lyrs(l,n,j) *  &
                                     f_root_pft_dz(j)
        END IF
      END DO
    END DO
  END DO
END IF

RETURN

END SUBROUTINE init_ic
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE total_snow_init(ainfo,progs)

USE layersnow_mod,   ONLY: layersnow

USE ancil_info, ONLY: land_pts, nsurft, surft_pts, lice_pts, nsoilt

USE jules_snow_mod, ONLY: nsmax, rho_snow_const, rho_snow_fresh, canSnowTile

USE jules_surface_mod, ONLY: l_elev_land_ice

USE jules_radiation_mod, ONLY: l_snow_albedo, l_embedded_snow

USE um_types, ONLY: real_jlslsm

!JULES TYPEs
USE ancil_info,    ONLY: ainfo_type
USE prognostics,   ONLY: progs_type

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Sets up the initial snow conditions from snow_surft only
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Arguments

!JULES TYPEs
TYPE(ainfo_type), INTENT(IN OUT) :: ainfo
TYPE(progs_type), INTENT(IN OUT) :: progs

! Work variables
INTEGER :: i, j, k, n, m  ! Index variables

REAL(KIND=real_jlslsm) :: snow_on_ground(land_pts,nsurft)
              ! Snow considered to be on the ground (not in canopy) (kg m-2)
              ! This is all the snow

!-------------------------------------------------------------------------------
! Put all snow onto the ground and zero canopy snow.
! Currently all snow is held in snow_surft.
! For can_model=4 tiles, put snow into snow_grnd_surft and zero snow_surft.
!-------------------------------------------------------------------------------

! Save input value.
snow_on_ground(:,:) = progs%snow_surft(:,:)

! Initialise stores to zero.
progs%snow_grnd_surft(:,:) = 0.0
progs%snow_surft(:,:) = 0.0

! Initialise other variables with values that will be retained where there is
! no tile - using "sensible" values for when these are printed.
progs%snowdepth_surft(:,:) = 0.0
IF ( nsmax < 1 ) THEN
  progs%rho_snow_grnd_surft(:,:) = rho_snow_const
ELSE
  progs%rho_snow_grnd_surft(:,:) = rho_snow_fresh
  progs%tsnow_surft(:,:,:) = 273.15
  progs%ds_surft(:,:,:) = 0.0
  IF ( l_snow_albedo .OR. l_embedded_snow ) progs%rgrainl_surft(:,:,:) = 0.0
END IF

DO n = 1,nsurft
  IF ( canSnowTile(n) ) THEN
    DO j = 1,surft_pts(n)
      i = ainfo%surft_index(j,n)
      progs%snow_grnd_surft(i,n) = snow_on_ground(i,n)
    END DO
  ELSE
    DO j = 1,surft_pts(n)
      i = ainfo%surft_index(j,n)
      progs%snow_surft(i,n) = snow_on_ground(i,n)
    END DO
  END IF
END DO

!-------------------------------------------------------------------------------
! Set snow density, calculate snow depth and set temperature of snow to equal
! that of soil.
!-------------------------------------------------------------------------------

!==============================================================================
!**NOTICE REGARDING SOIL TILING**
!
!The following section facilitates the use of soil tiling. As implemented,
!there are two soil tiling options:
!
!nsoilt == 1
!Operate as with a single soil tile, functionally identical to JULES upto
! at least vn4.7 (Oct 2016)
! This means that a soilt variable being passed 'up' to the surface is
! broadcast to the surft variable (with weighting by frac if requred)
!
!nsoilt > 1
!Operate with nsoilt = nsurft, with a direct mapping between them
! This means that a soilt variable being passed 'up' to the surface is simply
! copied into the surft variable
!
! This will need to be refactored for other tiling approaches. This note
! will be replicated elsewhere in the code as required
!
!These comments apply until **END NOTICE REGARDING SOIL TILING**
!==============================================================================

DO n = 1,nsurft

  !Set the current soil tile (see notice above)
  IF (nsoilt == 1) THEN
    !There is only 1 soil tile
    m = 1
  ELSE ! nsoilt == nsurft
    !Soil tiles map directly on to surface tiles
    m = n
  END IF !nsoilt

  DO j = 1,surft_pts(n)
    i = ainfo%surft_index(j,n)
    !     Use the constant (snowpack) density for nsmax=0 and if there is an
    !     existing pack. If nsmax>0 and there is no pack, initialise the density
    !     to the fresh snow value so that this value is used when/if a snowpack
    !     next develops.
    IF ( nsmax == 0 .OR.                                                      &
         ( snow_on_ground(i,n) > EPSILON(snow_on_ground) ) ) THEN
      progs%rho_snow_grnd_surft(i,n) = rho_snow_const
    ELSE
      progs%rho_snow_grnd_surft(i,n) = rho_snow_fresh
    END IF
    progs%snowdepth_surft(i,n) = snow_on_ground(i,n) /                        &
                                progs%rho_snow_grnd_surft(i,n)
    IF ( nsmax > 0 ) THEN

      progs%tsnow_surft(i,n,:) = progs%t_soil_soilt(i,m,1)
      IF ( l_snow_albedo .OR. l_embedded_snow )                               &
        progs%rgrainl_surft(i,n,:) = progs%rgrain_surft(i,n)
    END IF
  END DO
END DO

!==============================================================================
!**END NOTICE REGARDING SOIL TILING**
!==============================================================================

IF (l_elev_land_ice) THEN
  DO j = 1,lice_pts
    i = ainfo%lice_index(j)
    IF ( nsmax > 0 ) THEN
      DO n = 1,nsurft
        progs%tsnow_surft(i,n,:) = progs%tsurf_elev_surft(i,n)
      END DO
    END IF
  END DO
END IF

progs%nsnow_surft(:,:) = 0
IF ( nsmax > 0 ) THEN
  !-------------------------------------------------------------------------------
  ! Calculate snow layer thicknesses.
  !-------------------------------------------------------------------------------
  DO n = 1,nsurft
    CALL layersnow(land_pts, surft_pts(n), ainfo%surft_index(:,n),            &
                   progs%snowdepth_surft(:,n), progs%nsnow_surft(:,n),        &
                   progs%ds_surft(:,n,:))
  END DO

  !-------------------------------------------------------------------------------
  ! Set layer frozen and liquid contents.
  !-------------------------------------------------------------------------------
  progs%sice_surft(:,:,:) = 0.0
  progs%sliq_surft(:,:,:) = 0.0
  DO n = 1,nsurft
    DO j = 1,surft_pts(n)
      i = ainfo%surft_index(j,n)
      DO k = 1,progs%nsnow_surft(i,n)
        progs%sice_surft(i,n,k) = snow_on_ground(i,n) * progs%ds_surft(i,n,k) &
                              / progs%snowdepth_surft(i,n)
      END DO
    END DO
  END DO
END IF

RETURN

END SUBROUTINE total_snow_init
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE topmodel_init(psparms,toppdm,ainfo)

  !Use in relevant subroutines
USE calc_fsat_mod,           ONLY: calc_fsat
USE calc_zw_inund_mod,       ONLY: calc_zw_inund
USE calc_baseflow_jules_mod, ONLY: calc_baseflow_jules

!Use in relevant varaibles
USE ancil_info, ONLY:                                                         &
  land_pts, soil_pts, nsoilt

USE jules_soil_mod,       ONLY:                                               &
  dzsoil, sm_levels

USE jules_hydrology_mod,  ONLY:                                               &
  l_wetland_unfrozen

USE water_constants_mod,  ONLY:                                               &
  rho_water  !  density of pure water (kg/m3)

USE um_types, ONLY: real_jlslsm

!TYPE definitions
USE p_s_parms, ONLY: psparms_type
USE top_pdm, ONLY: top_pdm_type
USE ancil_info, ONLY: ainfo_type

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Finish initialising TOPMODEL by calculating surface saturated fraction
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Arguments

!TYPES containing the data
TYPE(psparms_type), INTENT(IN OUT) :: psparms
TYPE(top_pdm_type), INTENT(IN OUT) :: toppdm
TYPE(ainfo_type), INTENT(IN OUT) :: ainfo

! Work variables
LOGICAL, PARAMETER :: l_gamtot = .FALSE.  !  Switch for calculation of gamtot.
                                          !  .FALSE. so that calc_fsat
                                          !  calculates toppdm%fsat_soilt.

INTEGER :: i, j, n, m                    ! Index variables

REAL(KIND=real_jlslsm) :: qbase_l_soilt(land_pts,nsoilt,sm_levels+1)
  !Base flow from each layer (kg/m2/s).
REAL(KIND=real_jlslsm) :: top_crit_soilt(land_pts,nsoilt)
  !Critical topographic index required to calculate the surface saturation
  !fraction.
REAL(KIND=real_jlslsm) :: zdepth(0:sm_levels)
  !Lower soil layer boundary depth (m).
REAL(KIND=real_jlslsm) :: wutot_soilt(land_pts,nsoilt)
  !UNFROZEN to TOTAL fraction at ZW.
REAL(KIND=real_jlslsm) :: dumwutot_soilt(land_pts,nsoilt)
  !Dummy UNFROZEN to TOTAL fraction at ZW (always set to 1).
REAL(KIND=real_jlslsm) :: ksz_soilt(land_pts,nsoilt,0:sm_levels)
  !Saturated hydraulic conductivity for each layer (kg/m2/s).
REAL(KIND=real_jlslsm) :: dumsthf_soilt(land_pts,nsoilt,sm_levels)
  !Dummy Frozen soil moisture content of each layer as a fraction of
  !saturation (always set to 0).
REAL(KIND=real_jlslsm) :: smcl_soilt(land_pts,nsoilt,sm_levels)
  !Soil moisture content of each layer (kg/m2).
REAL(KIND=real_jlslsm) :: smclsat_soilt(land_pts,nsoilt,sm_levels)
  !Saturated soil moisture content of each layer (kg/m2).
REAL(KIND=real_jlslsm) :: zw_inund_soilt(land_pts,nsoilt)
  ! Adjusted Water table depth (m).

!-------------------------------------------------------------------------------
zdepth(0) = 0.0
DO n = 1,sm_levels
  zdepth(n) = zdepth(n-1) + dzsoil(n)
END DO

! Set values that are retained at non-soil points.
toppdm%fsat_soilt(:,:)     = 0.0
toppdm%fwetl_soilt(:,:)    = 0.0
dumwutot_soilt(:,:) = 1.0

IF ( soil_pts /= 0 ) THEN
  DO j = 1,soil_pts
    i = ainfo%soil_index(j)
    DO n = 0,sm_levels
      ksz_soilt(i,:,n) = psparms%satcon_soilt(i,:,n)
    END DO
    DO n = 1,sm_levels
      smclsat_soilt(i,:,n) = rho_water * dzsoil(n) * psparms%smvcst_soilt(i,:,n)
      smcl_soilt(i,:,n)    =                                                  &
        (psparms%sthu_soilt(i,:,n) + psparms%sthf_soilt(i,:,n))               &
        * smclsat_soilt(i,:,n)
    END DO
  END DO

  IF (L_wetland_unfrozen) THEN
    DO m = 1,nsoilt
      CALL calc_zw_inund(                                                     &
        land_pts,sm_levels,soil_pts,ainfo%soil_index,zdepth,                  &
        psparms%bexp_soilt(:,m,1),psparms%sathh_soilt(:,m,1),                 &
        smclsat_soilt(:,m,:),                                                 &
        smcl_soilt(:,m,:),psparms%sthu_soilt(:,m,:),toppdm%sthzw_soilt(:,m),  &
        toppdm%zw_soilt(:,m),                                                 &
        zw_inund_soilt(:,m),wutot_soilt(:,m))

    END DO
    dumsthf_soilt(:,:,:) = 0.0
  ELSE
    zw_inund_soilt(:,:)  = toppdm%zw_soilt(:,:)
    dumsthf_soilt(:,:,:) = psparms%sthf_soilt(:,:,:)
  END IF

  !   We need top_crit_soilt - get this from calc_baseflow.
  DO m = 1,nsoilt
    CALL calc_baseflow_jules(                                                 &
      soil_pts, ainfo%soil_index, land_pts, sm_levels, zdepth,                &
      ksz_soilt(:,m,:), psparms%bexp_soilt(:,m,:), toppdm%fexp_soilt(:,m),    &
      toppdm%ti_mean_soilt(:,m), zw_inund_soilt(:,m), dumsthf_soilt(:,m,:),   &
      top_crit_soilt(:,m), toppdm%qbase_soilt(:,m), qbase_l_soilt(:,m,:))

    !   Call calc_fsat with 1st argument (l_gamtot)=.FALSE. so as to
    !   get toppdm%fsat_soilt.
    CALL calc_fsat(                                                           &
      l_gamtot, soil_pts, ainfo%soil_index, land_pts,                         &
      toppdm%ti_mean_soilt(:,m), toppdm%ti_sig_soilt(:,m),                    &
      dumwutot_soilt(:,m), top_crit_soilt(:,m),                               &
      toppdm%gamtot_soilt(:,m), toppdm%fsat_soilt(:,m), toppdm%fwetl_soilt(:,m))
  END DO

  IF (L_wetland_unfrozen) THEN
    DO j = 1,soil_pts
      i = ainfo%soil_index(j)
      toppdm%fsat_soilt(i,:)  = wutot_soilt(i,:) * toppdm%fsat_soilt(i,:)
      toppdm%fwetl_soilt(i,:) = wutot_soilt(i,:) * toppdm%fwetl_soilt(i,:)
    END DO
  END IF

END IF  !  soil_pts

RETURN

END SUBROUTINE topmodel_init
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

FUNCTION get_default_ic_values(total_snow) RESULT(defaults_dict)

! USE statments
USE dictionary_mod, ONLY: dict, dict_set, dict_create
USE logging_mod,    ONLY: log_info

USE dump_mod,       ONLY: max_var_dump, required_vars_for_configuration

USE model_interface_mod, ONLY: identifier_len

USE fire_mod,        ONLY: fire_inis

USE lake_mod,        ONLY: lake_fetch_0, lake_h_mxl_0, lake_shape_0,          &
                           lake_T_mean_0, lake_T_mxl_0, lake_T_ice_0,         &
                           lake_H_ice_0

USE metstats_mod,    ONLY: metstats_inis

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Gets default values for initialisation
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Arguments
LOGICAL, INTENT(IN) :: total_snow

CHARACTER(LEN=identifier_len) :: identifiers(max_var_dump)
                             ! The variable identifiers of the required
                             ! variables

CHARACTER(LEN=identifier_len) :: vars_from_ancil(max_var_dump)
                             ! The variable identifiers of the ancil
                             ! variables (not used in this subroutine)

! Work variables
INTEGER :: i !Counter
INTEGER :: nvars_required
INTEGER :: nvars_from_ancil

! Return type
TYPE(dict) :: defaults_dict

LOGICAL, PARAMETER :: l_output_mode = .FALSE.

!End of header
!-----------------------------------------------------------------------------

  !Get the list of variables we need to have default variables for
CALL required_vars_for_configuration(nvars_required, identifiers,             &
                                     nvars_from_ancil, vars_from_ancil,       &
                                     l_output_mode,                           &
                                     total_snow, .FALSE.)

!Create the dictionary containing the default values
defaults_dict = dict_create(nvars_required, 1.0)

!Set the default values one by one
DO i = 1, nvars_required
  SELECT CASE ( identifiers(i) )

    !      CASE ( 'canopy' )
    !      CASE ( 'cs' )
    !      CASE ( 'gs' )
    !      CASE ( 'snow_tile' )
    !      CASE ( 'sthuf' )
    !      CASE ( 't_soil' )
    !      CASE ( 'tstar_tile' )
    !      CASE ( 'lai' )
    !      CASE ( 'canht' )
    !      CASE ( 'frac' )
    !      CASE ( 'sthzw' )
    !      CASE ( 'zw' )
    !      CASE ( 'rgrain' )

          !Metstat module variables
  CASE ( 'temp_max_00h_r' )
    CALL dict_set(defaults_dict, 'temp_max_00h_r', metstats_inis%temp_max_00h)
  CASE ( 'temp_ave_00h_r' )
    CALL dict_set(defaults_dict, 'temp_ave_00h_r', metstats_inis%temp_ave_00h)
  CASE ( 'prec_tot_00h_r' )
    CALL dict_set(defaults_dict, 'prec_tot_00h_r', metstats_inis%prec_tot_00h)
  CASE ( 'prec_tot_12h_r' )
    CALL dict_set(defaults_dict, 'prec_tot_12h_r', metstats_inis%prec_tot_12h)
  CASE ( 'rhum_min_00h_r' )
    CALL dict_set(defaults_dict, 'rhum_min_00h_r', metstats_inis%rhum_min_00h)
  CASE ( 'dewp_ave_00h_r' )
    CALL dict_set(defaults_dict, 'dewp_ave_00h_r', metstats_inis%dewp_ave_00h)
  CASE ( 'wind_ave_00h_r' )
    CALL dict_set(defaults_dict, 'wind_ave_00h_r', metstats_inis%wind_ave_00h)
  CASE ( 'temp_max_00h' )
    CALL dict_set(defaults_dict, 'temp_max_00h',   metstats_inis%temp_max_00h)
  CASE ( 'temp_ave_00h' )
    CALL dict_set(defaults_dict, 'temp_ave_00h',   metstats_inis%temp_ave_00h)
  CASE ( 'temp_ave_nday' )
    CALL dict_set(defaults_dict, 'temp_ave_nday',  metstats_inis%temp_ave_nday)
  CASE ( 'temp_pnt_12h' )
    CALL dict_set(defaults_dict, 'temp_pnt_12h',   metstats_inis%temp_pnt_12h)
  CASE ( 'prec_tot_00h' )
    CALL dict_set(defaults_dict, 'prec_tot_00h',   metstats_inis%prec_tot_00h)
  CASE ( 'prec_tot_12h' )
    CALL dict_set(defaults_dict, 'prec_tot_12h',   metstats_inis%prec_tot_12h)
  CASE ( 'rhum_min_00h' )
    CALL dict_set(defaults_dict, 'rhum_min_00h',   metstats_inis%rhum_min_00h)
  CASE ( 'rhum_pnt_12h' )
    CALL dict_set(defaults_dict, 'rhum_pnt_12h',   metstats_inis%rhum_pnt_12h)
  CASE ( 'dewp_ave_00h' )
    CALL dict_set(defaults_dict, 'dewp_ave_00h',   metstats_inis%dewp_ave_00h)
  CASE ( 'wind_ave_00h' )
    CALL dict_set(defaults_dict, 'wind_ave_00h',   metstats_inis%wind_ave_00h)
  CASE ( 'wind_pnt_12h' )
    CALL dict_set(defaults_dict, 'wind_pnt_12h',   metstats_inis%wind_pnt_12h)

    !Fire module variables
  CASE ( 'fire_canadian_ffmc' )
    CALL dict_set(defaults_dict, 'fire_canadian_ffmc'     , fire_inis%canadian_ffmc )
  CASE ( 'fire_canadian_ffmc_mois' )
    CALL dict_set(defaults_dict, 'fire_canadian_ffmc_mois', fire_inis%canadian_ffmc_mois)
  CASE ( 'fire_canadian_dmc' )
    CALL dict_set(defaults_dict, 'fire_canadian_dmc'      , fire_inis%canadian_dmc )
  CASE ( 'fire_canadian_dc' )
    CALL dict_set(defaults_dict, 'fire_canadian_dc'       , fire_inis%canadian_dc )
  CASE ( 'fire_mcarthur_r_dr' )
    CALL dict_set(defaults_dict, 'fire_mcarthur_r_dr'     , fire_inis%mcarthur_r_dr )
  CASE ( 'fire_mcarthur_n_dr' )
    CALL dict_set(defaults_dict, 'fire_mcarthur_n_dr'     , fire_inis%mcarthur_n_dr)
  CASE ( 'fire_nesterov' )
    CALL dict_set(defaults_dict, 'fire_nesterov'          , fire_inis%nesterov_index )

    !FLake variables
  CASE ( 'lake_fetch_gb' )
    CALL dict_set(defaults_dict, 'lake_fetch_gb'       , lake_fetch_0)
  CASE ( 'lake_h_mxl_gb' )
    CALL dict_set(defaults_dict, 'lake_h_mxl_gb'       , lake_h_mxl_0)
  CASE ( 'lake_shape_factor_gb' )
    CALL dict_set(defaults_dict, 'lake_shape_factor_gb', lake_shape_0)
  CASE ( 'lake_t_mean_gb' )
    CALL dict_set(defaults_dict, 'lake_t_mean_gb'      , lake_T_mean_0)
  CASE ( 'lake_t_mxl_gb' )
    CALL dict_set(defaults_dict, 'lake_t_mxl_gb'       , lake_T_mxl_0)
  CASE ( 'lake_t_ice_gb' )
    CALL dict_set(defaults_dict, 'lake_t_ice_gb'       , lake_T_ice_0)
  CASE ( 'lake_h_ice_gb' )
    CALL dict_set(defaults_dict, 'lake_h_ice_gb'       , lake_H_ice_0)

    !Snow variables
    !      CASE ( 'rho_snow' )
    !      CASE ( 'snow_depth' )
    !      CASE ( 'snow_grnd' )
    !      CASE ( 'nsnow' )
    !      CASE ( 'snow_ds' )
    !      CASE ( 'snow_ice' )
    !      CASE ( 'snow_liq' )
    !      CASE ( 'tsnow' )

  CASE DEFAULT
    CALL log_info("get_default_ic_values",                                    &
                   "No default value defined. Value must be defined in namelist - " &
                   // TRIM(identifiers(i)))
  END SELECT
END DO

RETURN

END FUNCTION get_default_ic_values
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   Subroutine CALC_FIT_FSAT-------------------------------------------
!
!   Purpose: To speed up the large scale hydrology code (LTOP=TRUE)
!            dramatically. This is done by fitting exponential
!            functions to the incomplete gamma function for each grid
!            box and the complete range of possible "water table"
!            (top_crit) cases - see documentation.
!            Estimates the fitted parameters for Fsat=function(ZW)
!            and  Fwet=function(ZW) for each land grid point.
!            (Calculating the incomplete gamma function for each grid
!            box at each time step was very time consuming).
!                                                             !
! Documentation: UNIFIED MODEL DOCUMENTATION PAPER NO 25
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  6.4    10/01/07   New Deck         Nic Gedney
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.

SUBROUTINE calc_fit_fsat(toppdm,ainfo)

  !Use in relevant subroutines
USE calc_fsat_mod, ONLY: calc_fsat

!Use in relevant varaibles
USE ancil_info,           ONLY:                                               &
  land_pts, soil_pts, nsoilt

USE jules_hydrology_mod,  ONLY:                                               &
  l_top, zw_max, nfita

USE jules_soil_mod,       ONLY:                                               &
  dzsoil, sm_levels

USE jules_print_mgr,      ONLY:                                               &
  jules_message, jules_print

USE um_types, ONLY: real_jlslsm

! Type Definitions
USE top_pdm, ONLY: top_pdm_type
USE ancil_info,    ONLY: ainfo_type

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Subroutine arguments
! TYPE Definitions
TYPE(top_pdm_type), INTENT(IN OUT) :: toppdm
TYPE(ainfo_type), INTENT(IN OUT) :: ainfo

!-----------------------------------------------------------------------------

! Work variables
REAL(KIND=real_jlslsm) :: zdepth           ! Standard Soil model DEPTH.

! Local scalars:
INTEGER, PARAMETER :: nzw = 100  ! Number of ZW values used in fit.
                                 ! Maximum value for a significant
                                 ! improvement in the fit.

INTEGER :: i,j,iz,n,ifita,m  ! Loop counters.

REAL(KIND=real_jlslsm) ::                                                     &
  dzw,                                                                        &
                 ! WORK ZW increment ; defined by zw_max and nzw.
  rms,                                                                        &
                 ! WORK rms errors for given fsat fit values.
  rmsw,                                                                       &
                 ! WORK rms errors for given fwet fit values.
  rmsold,                                                                     &
                 ! WORK rms errors for given fsat fit values.
                 !      for best fit so far.
  rmswold,                                                                    &
                 ! WORK rms errors for given fwet fit values
                 !      for best fit so far.
  cfit,                                                                       &
                 ! WORK CFit value for given loop.
  cfitmax
                 ! Maximum possible value for Cfit.

REAL(KIND=real_jlslsm), PARAMETER ::                                          &
  thr_err = 5.0e-3,                                                           &
                 ! Error threshold value

  cfitmin = 0.0
                 ! Minimum possible value for Cfit.

! Local arrays:
REAL(KIND=real_jlslsm) ::                                                     &
  fsat_calc_soilt(land_pts,nsoilt,nzw),                                       &
                            ! WORK Surface saturation fraction.
  fsat_fit(nzw),                                                              &
                            ! WORK Fitted surface saturation fraction.
  fwet_calc_soilt(land_pts,nsoilt,nzw),                                       &
                            ! WORK Wetland fraction.
  fwet_fit(nzw),                                                              &
                            ! WORK Fitted wetland fraction.
  dumzw(nzw),                                                                 &
                            ! WORK Dummy water table depth (m).
  dumfsat_soilt(land_pts,nsoilt),                                             &
                            ! WORK Dummy surface saturation fraction.
  dumfwetl_soilt(land_pts,nsoilt)
                            ! WORK Dummy wetland fraction.
! DBC We could use this local version of gamtot in this version,
!     but using the module version to match other versions of code.
!     &,toppdm%gamtot_soilt(land_pts)  ! WORK Integrated complete Gamma function

REAL(KIND=real_jlslsm) ::                                                     &
  top_crit_soilt(land_pts,nsoilt,nzw),                                        &
                            ! WORK LOG(QBASE_MAX/QBASE) -see document.
  top_crit1z_soilt(land_pts,nsoilt),                                          &
                            ! WORK As above but for an individual zw.
  top_min_soilt(land_pts,nsoilt),                                             &
                            ! WORK value for when zw=zw_max.
  wutot_soilt(land_pts,nsoilt)           ! WORK Dummy (set to 1.0).

INTEGER :: errorstatus
CHARACTER(LEN=80) :: cmessage
CHARACTER(LEN=13) :: routinename


!-----------------------------------------------------------------------------
routinename='CALC_FIT_FSAT'

IF ( l_top ) THEN

  cfitmax = 0.15 * nfita

  ! Define the water table depths to be used in the fitting process:
  dzw = 1.0 / REAL(nzw) * zw_max
  DO iz = 1,nzw
    dumzw(iz) = REAL(iz-1) * dzw
  END DO

  DO i = 1,land_pts        ! initialise to zero
    wutot_soilt(i,:)    = 1.0
    dumfsat_soilt(i,:)  = 0.0
    dumfwetl_soilt(i,:) = 0.0
  END DO
  zdepth = 0.0

  ! Calculate total soil depth
  DO n = 1,sm_levels
    zdepth = zdepth + dzsoil(n)
  END DO

  ! Calculate Gamtot
  DO m = 1,nsoilt
    toppdm%gamtot_soilt(:,m) = 0.0
    CALL calc_fsat( .TRUE., soil_pts, ainfo%soil_index, land_pts,             &
                  toppdm%ti_mean_soilt(:,m),                                  &
                  toppdm%ti_sig_soilt(:,m), wutot_soilt(:,m),                 &
                  top_crit1z_soilt(:,m), toppdm%gamtot_soilt(:,m),            &
                  dumfsat_soilt(:,m), dumfwetl_soilt(:,m))
  END DO

  ! Calculate top_crit for the water table depths:
  DO iz = 1,nzw
    DO m = 1, nsoilt
      DO j = 1,soil_pts
        i = ainfo%soil_index(j)

        fsat_calc_soilt(i,m,iz) = 0.0
        fwet_calc_soilt(i,m,iz) = 0.0
        top_crit_soilt(i,m,iz)  = 0.0

        IF ( toppdm%ti_mean_soilt(i,m) > 0.0 .AND.                            &
             toppdm%ti_sig_soilt(i,m) > 0.0 ) THEN
          top_min_soilt(i,m) = 1.0 / toppdm%fexp_soilt(i,m)                   &
                              * EXP(-toppdm%fexp_soilt(i,m) * (zw_max - zdepth))

          IF ( dumzw(iz) <= zdepth ) THEN
            top_crit1z_soilt(i,m) = LOG(zdepth+1.0                            &
                               / toppdm%fexp_soilt(i,m) - top_min_soilt(i,m)) &
                               - LOG(zdepth - dumzw(iz) + 1.0                 &
                               / toppdm%fexp_soilt(i,m) - top_min_soilt(i,m))
          END IF

          IF ( dumzw(iz) >  zdepth ) THEN
            top_crit1z_soilt(i,m) = LOG(zdepth+1.0                            &
                                    / toppdm%fexp_soilt(i,m)                  &
                                    - top_min_soilt(i,m))                     &
                                    - LOG(1 / toppdm%fexp_soilt(i,m)          &
                                    * EXP(-toppdm%fexp_soilt(i,m)             &
                                    * (dumzw(iz) - zdepth))                   &
                                    - top_min_soilt(i,m))
          END IF
        END IF
      END DO !soil_pts

      ! Calculate FSAT and FWET for one ZW at all soil land_pts:
      CALL calc_fsat( .FALSE., soil_pts, ainfo%soil_index, land_pts,          &
                     toppdm%ti_mean_soilt(:,m),                               &
                     toppdm%ti_sig_soilt(:,m), wutot_soilt(:,m),              &
                     top_crit1z_soilt(:,m), toppdm%gamtot_soilt(:,m),         &
                     dumfsat_soilt(:,m), dumfwetl_soilt(:,m))

      DO j = 1,soil_pts
        i = ainfo%soil_index(j)

        IF ( toppdm%ti_mean_soilt(i,m) > 0.0 .AND.                            &
             toppdm%ti_sig_soilt(i,m) > 0.0 ) THEN
          fsat_calc_soilt(i,m,iz) = dumfsat_soilt(i,m)
          fwet_calc_soilt(i,m,iz) = dumfwetl_soilt(i,m)
          top_crit_soilt(i,m,iz)  = top_crit1z_soilt(i,m)

          IF ( dumzw(iz) <  dzw ) THEN ! Values at zw=0m
            toppdm%a_fsat_soilt(i,m) = fsat_calc_soilt(i,m,iz)
            toppdm%a_fwet_soilt(i,m) = fwet_calc_soilt(i,m,iz)
          END IF
        END IF
      END DO
    END DO !nsoilt
  END DO                  !ZW calc_fsat loop

  ! Now carry out fit for FSAT, where FSAT=function(ZW). (Likewise FWET)
  DO m = 1, nsoilt
    DO j = 1,soil_pts
      i = ainfo%soil_index(j)

      IF ( toppdm%ti_mean_soilt(i,m) > 0.0 .AND.                              &
           toppdm%ti_sig_soilt(i,m) > 0.0 ) THEN
        rmsold  = 1.0e10
        rmswold = 1.0e10

        DO ifita = 1,nfita
          cfit = cfitmax * (ifita) / REAL(nfita)

          rms  = 0.0
          rmsw = 0.0
          ! top_crit=TI_MAX when zw=zw_max
          DO iz = 1,nzw
            fsat_fit(iz) = toppdm%a_fsat_soilt(i,m)                           &
                            * EXP(-cfit * top_crit_soilt(i,m,iz))
            fwet_fit(iz) = toppdm%a_fwet_soilt(i,m)                           &
                            * EXP(-cfit * top_crit_soilt(i,m,iz))
            rms  = rms  + (fsat_calc_soilt(i,m,iz) - fsat_fit(iz))**2
            rmsw = rmsw + (fwet_calc_soilt(i,m,iz) - fwet_fit(iz))**2
          END DO         !ZW

          rms  = SQRT(rms) / REAL(nzw)
          rmsw = SQRT(rmsw) / REAL(nzw)

          IF ( rms < rmsold ) THEN
            rmsold = rms
            toppdm%c_fsat_soilt(i,m) = cfit
          END IF
          IF ( rmsw < rmswold ) THEN
            rmswold = rmsw
            toppdm%c_fwet_soilt(i,m) = cfit
          END IF
        END DO

        DO iz = 1,nzw
          fsat_fit(iz) = toppdm%a_fsat_soilt(i,m) *                           &
                         EXP(-toppdm%c_fsat_soilt(i,m) *                      &
                         top_crit_soilt(i,m,iz))
          fwet_fit(iz) = toppdm%a_fwet_soilt(i,m) *                           &
                         EXP(-toppdm%c_fwet_soilt(i,m) *                      &
                         top_crit_soilt(i,m,iz))
        END DO            !ZW

        IF ( rmsold >= thr_err ) THEN
          IF ( toppdm%c_fsat_soilt(i,m) <= cfitmin .OR.                       &
              toppdm%c_fsat_soilt(i,m) >= cfitmax ) THEN

            WRITE(jules_message,*) 'ERROR cfit FSAT', i,                      &
                toppdm%c_fsat_soilt(i,m), cfitmin, cfitmax
            CALL jules_print('calc_fit_fsat',jules_message)
            WRITE(jules_message,*)                                            &
              'If c_fsat=cfitmax try increasing nfita'
            CALL jules_print('calc_fit_fsat',jules_message)
            WRITE(jules_message,*) 'fsat_calc=',                              &
                fsat_calc_soilt(i,m,1), fsat_calc_soilt(i,m,3),               &
                fsat_calc_soilt(i,m,5)
            CALL jules_print('calc_fit_fsat',jules_message)
            WRITE(jules_message,*) 'fsat_fit=',                               &
                fsat_fit(1), fsat_fit(3), fsat_fit(5)
            CALL jules_print('calc_fit_fsat',jules_message)
            WRITE(jules_message,*) 'rms=', rmsold
            CALL jules_print('calc_fit_fsat',jules_message)
            ErrorStatus = 35
            WRITE(CMessage, *) 'Error in cfit FSAT in LSH model setup'
          END IF
        END IF

        IF ( rmswold >= thr_err ) THEN
          IF ( toppdm%c_fwet_soilt(i,m) <= cfitmin .OR.                       &
              toppdm%c_fwet_soilt(i,m) >= cfitmax ) THEN

            WRITE(jules_message,*) 'ERROR cfit FWET', i,                      &
                toppdm%c_fwet_soilt(i,m), cfitmin, cfitmax
            CALL jules_print('calc_fit_fsat',jules_message)
            WRITE(jules_message,*) 'fwet_calc=',                              &
                fwet_calc_soilt(i,m,1), fwet_calc_soilt(i,m,3),               &
                fwet_calc_soilt(i,m,5)
            CALL jules_print('calc_fit_fsat',jules_message)
            WRITE(jules_message,*) 'fwet_fit=',                               &
                fwet_fit(1), fwet_fit(3), fwet_fit(5)
            CALL jules_print('calc_fit_fsat',jules_message)
            WRITE(jules_message,*) 'rmsw=', rmswold
            CALL jules_print('calc_fit_fsat',jules_message)
            WRITE(jules_message,*) '(fsat_calc=)',                            &
                fsat_calc_soilt(i,m,1), fsat_calc_soilt(i,m,3),               &
                fsat_calc_soilt(i,m,5)
            CALL jules_print('calc_fit_fsat',jules_message)
            WRITE(jules_message,*) '(fsat_fit=)',                             &
                fsat_fit(1), fsat_fit(3), fsat_fit(5)
            CALL jules_print('calc_fit_fsat',jules_message)
            WRITE(jules_message,*) '(rms=)', rmsold
            CALL jules_print('calc_fit_fsat',jules_message)
            ErrorStatus = 35
            ErrorStatus = 40
            WRITE(CMessage,*) 'Error in cfit FWET in LSH model setup'
          END IF
        END IF

        IF ( rmsold > thr_err ) THEN
          WRITE(jules_message,*) 'Warning LSH rms Error in fit:',             &
              rmsold, rmswold
          CALL jules_print('calc_fit_fsat',jules_message)
        END IF
      END IF

    END DO                  ! land_pts
  END DO !nsoilt
ELSE
  DO m = 1,nsoilt
    DO j = 1,soil_pts
      i = ainfo%soil_index(j)

      toppdm%a_fsat_soilt(i,m) = 0.0
      toppdm%c_fsat_soilt(i,m) = 0.0
      toppdm%a_fwet_soilt(i,m) = 0.0
      toppdm%c_fwet_soilt(i,m) = 0.0
      toppdm%gamtot_soilt(i,m) = 0.0
    END DO
  END DO
END IF

RETURN

END SUBROUTINE calc_fit_fsat

END MODULE initial_conditions_mod

