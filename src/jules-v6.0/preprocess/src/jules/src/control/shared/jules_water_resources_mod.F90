!******************************COPYRIGHT**************************************
! (c) UK Centre for Ecology & Hydrology.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT**************************************

MODULE jules_water_resources_mod

USE um_types, ONLY: real_jlslsm

USE missing_data_mod, ONLY:                                                   &
  ! imported scalar parameters
  imdi, rmdi

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Contains water resource management options, parameters and variables,
!   and a namelist for setting some of them.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Public scope by default.

!-----------------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                       &
  ModuleName = 'JULES_WATER_RESOURCES_MOD'

INTEGER, PARAMETER, PRIVATE ::                                                &
  name_len = 3,                                                               &
    ! Length of sector names.
  nwater_use_max = 6
    ! Maximum possible number of water sectors.

! Parameters identifying alternative models for non-renewable groundwater.
! These should have unique, non-zero values. A value of zero is used to
! indicate that non-renewable groundwater is not used.
INTEGER, PARAMETER ::                                                         &
  nr_gwater_model_last = 1,                                                   &
    ! Indicates that non-renewable groundwater is used as a last resort.
  nr_gwater_model_use = 2
    ! Indicates that non-renewable groundwater is used as part of the mix of
    ! water sources.

CHARACTER(LEN=name_len), PARAMETER ::                                         &
  name_domestic = 'dom',                                                      &
    ! Name used to identify domestic use.
  name_environment = 'env',                                                   &
    ! Name used to identify environmental use.
  name_industry = 'ind',                                                      &
    ! Name used to identify industrial use.
  name_irrigation = 'irr',                                                    &
    ! Name used to identify irrigation use.
  name_livestock = 'liv',                                                     &
    ! Name used to identify livestock use.
  name_transfers = 'tra'
    ! Name used to identify water transfers use.

!-----------------------------------------------------------------------------
! Module variables
!-----------------------------------------------------------------------------

! Items set in namelist jules_water_resources.

! Top-level switch for the parameterisation.
LOGICAL ::                                                                    &
  l_water_resources = .FALSE.
    ! Switch to select water resource management modelling.
    ! .TRUE.  = represent water resources
    ! .FALSE. = no water resources

! Switches that control which sectors are considered.
LOGICAL ::                                                                    &
  l_water_domestic = .FALSE.,                                                 &
    ! .TRUE. =  consider demand for water for domestic use
    ! .FALSE. = do not consider domestic demand
  l_water_environment = .FALSE.,                                              &
    ! .TRUE.  = consider demand for water for environmental flow requirements
    ! .FALSE. = do not consider environmental demand
  l_water_industry = .FALSE.,                                                 &
    ! .TRUE.  = consider demand for water for industrial use
    ! .FALSE. = do not consider industrial demand
  l_water_irrigation = .FALSE.,                                               &
    ! .TRUE.  = consider demand for water for irrigation
    ! .FALSE. = do not consider irrigation demand
    ! In future this will be made consistent with existing irrigation switches
    ! such as l_irrig_dmd; at present there is no conflicting functionality.
  l_water_livestock = .FALSE.,                                                &
    ! .TRUE.  = consider demand for water for livestock
    ! .FALSE. = do not consider demand from livestock
  l_water_transfers = .FALSE.
    ! .TRUE.  = consider (explicit) water transfers
    ! .FALSE. = do not consider transfers

LOGICAL ::                                                                    &
  l_prioritise = .FALSE.
    ! Switch controlling prioritisation beween demands.
    ! .TRUE.  = rank demands in priority order
    ! .FALSE. = no prioritisation

! Other parameters of the model.
INTEGER ::                                                                    &
  nr_gwater_model = imdi,                                                     &
    ! Chosen model for non-renewable groundwater.
  nstep_water_res = imdi
    ! Timestep length for water resource model (number of "main model"
    ! timesteps).

REAL(KIND=real_jlslsm) ::                                                     &
  rf_domestic = rmdi,                                                         &
    ! Fraction of water that is returned after abstraction for domestic use
    ! (via sewage systems etc.).
  rf_industry = rmdi,                                                         &
    ! Fraction of water that is returned after abstraction for industrial
    ! purposes.
  rf_livestock = rmdi
    ! Fraction of water that is returned after abstraction for livestock.

CHARACTER(LEN=name_len) ::                                                    &
  priority(nwater_use_max) = 'xxx'
    ! Water sector names, in order of decreasing priority.
    ! This fixed-length variable allows for all possible sectors.

!-----------------------------------------------------------------------------
! Declare the namelist.
!-----------------------------------------------------------------------------
NAMELIST  / jules_water_resources /                                           &
! Shared
    l_prioritise, l_water_domestic, l_water_environment, l_water_industry,    &
    l_water_irrigation, l_water_livestock, l_water_resources,                 &
    l_water_transfers, nr_gwater_model, nstep_water_res, priority,            &
    rf_domestic, rf_industry, rf_livestock

!-----------------------------------------------------------------------------
! Variables below here are not in the namelist.
!-----------------------------------------------------------------------------
! Scalar variables.
INTEGER ::                                                                    &
  nwater_use,                                                                 &
    ! Number of water resource sectors that are considered.
  use_domestic = 0,                                                           &
    ! Index of domestic use in multi-use arrays.
  use_environment = 0,                                                        &
    ! Index of environmental use in multi-use arrays.
  use_industry = 0,                                                           &
    ! Index of industrial use in multi-use arrays.
  use_irrigation = 0,                                                         &
    ! Index of irrigation use in multi-use arrays.
  use_livestock = 0,                                                          &
    ! Index of livestock use in multi-use arrays.
  use_transfers = 0,                                                          &
    ! Index of transfers in multi-use arrays.
  water_res_count = 0
    ! Counter of timesteps done in current water resource timestep.

! Array variables.
INTEGER, ALLOCATABLE ::                                                       &
  priority_order(:,:)
    ! Priorities of water demands at each gridpoint, in order of decreasing
    ! priority. Values are the index in multi-sector arrays.
    ! e.g. priority_order(l,1) = 3 indicates that the first priority use is
    !      in slice 3 of multi-sector arrays.

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  demand_accum(:,:),                                                          &
    ! Demands for water accumulated over the water resource timestep (kg).
    ! Note that in general this should be written to restart files (dumps)
    ! but this is not done yet.
  !---------------------------------------------------------------------------
  ! Demands that can be prescribed.
  !---------------------------------------------------------------------------
  demand_rate_domestic(:),                                                    &
    ! Demand for water for domestic use (kg s-1).
  demand_rate_industry(:),                                                    &
    ! Demand for water for industrial use (kg s-1).
  demand_rate_livestock(:),                                                   &
    ! Demand for water for livestock (kg s-1).
  demand_rate_transfers(:),                                                   &
    ! Demand for water for (explicit) transfers (kg s-1).
  !---------------------------------------------------------------------------
  ! Ancillary fields.
  !---------------------------------------------------------------------------
  conveyance_loss(:),                                                         &
    ! Fraction of water that is lost during conveyance from source to user.
  irrig_eff(:),                                                               &
    ! Irrigation efficiency i.e. the fraction of the water withdrawn for
    ! irrigation that is demanded by the crop scheme.
  sfc_water_frac(:)
    ! Target for the fraction of demand that will be met from surface water
    ! (as opposed to groundwater).

CONTAINS

!#############################################################################
!#############################################################################

SUBROUTINE check_jules_water_resources()

USE ereport_mod, ONLY: ereport

USE jules_rivers_mod, ONLY: l_rivers, nstep_rivers, np_rivers

USE jules_soil_mod, ONLY: sm_levels

USE jules_irrig_mod, ONLY: l_irrig_dmd, l_irrig_limit

!-----------------------------------------------------------------------------
! Description:
!   Checks that values in JULES_WATER_RESOURCES namelist were provided and
!   are acceptable.
!-----------------------------------------------------------------------------

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Local scalar parameters.
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER ::                                                &
   RoutineName = 'CHECK_JULES_WATER_RESOURCES'   ! Name of this procedure.

!-----------------------------------------------------------------------------
! Local scalar variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  error_status,                                                               &
    ! Error status.
  i,                                                                          &
    ! Loop counter.
  nsector
    ! Number of sectors selected.

LOGICAL ::                                                                    &
  in_list
    ! T when a string appears in a list.

!-----------------------------------------------------------------------------
! Set error status to show a fatal error.
error_status = 101

! If water resources are not to be modelled, there's nothing to do here.
IF ( .NOT. l_water_resources ) RETURN

!-----------------------------------------------------------------------------
! Check that a timestep was provided and is reasonable.
!-----------------------------------------------------------------------------
IF ( nstep_water_res == imdi ) THEN
  CALL ereport ( RoutineName, error_status,                                   &
                 "nstep_water_res not found." )
END IF
IF ( nstep_water_res < 1 ) THEN
  CALL ereport ( RoutineName, error_status,                                   &
                 "nstep_water_res must be at least 1." )
END IF

! Water resource and river models must be in synchrony.
! It is likely that future code additions will include the possibility that
! the configuration includes water resources but with no need to access
! rivers. At that time this test would be revised.
IF ( l_rivers .AND. nstep_water_res /= nstep_rivers ) THEN
  CALL ereport ( RoutineName, error_status,                                   &
                 "Water resource and river modules must be called on the " // &
                 "same timestep, i.e. nstep_water_res = nstep_rivers." )
END IF

!-----------------------------------------------------------------------------
! Count the number of sectors to be considered, and check that at least one is
! selected.
!-----------------------------------------------------------------------------
nwater_use = COUNT( (/ l_water_domestic, l_water_environment,                 &
                       l_water_industry, l_water_irrigation,                  &
                       l_water_livestock, l_water_transfers /) )
IF ( nwater_use == 0 ) THEN
  CALL ereport ( RoutineName, error_status,                                   &
                 "At least one resource sector must be selected." )
END IF

IF ( nwater_use > nwater_use_max ) THEN
  CALL ereport ( RoutineName, error_status,                                   &
                 "Increase size of nwater_use_max." )
END IF

!-----------------------------------------------------------------------------
! Check that a return flow value is given for all relevant and active sectors,
! and that it is reasonable.
!-----------------------------------------------------------------------------
IF ( l_water_domestic ) THEN
  IF ( ABS( rf_domestic - rmdi ) < EPSILON(1.0) ) THEN
    CALL ereport ( RoutineName, error_status,                                 &
                   "rf_domestic not found." )
  ELSE IF ( rf_domestic < 0.0 .OR. rf_domestic > 1.0 ) THEN
    CALL ereport ( RoutineName, error_status,                                 &
                   "rf_domestic must lie in the range 0 to 1." )
  END IF
END IF

IF ( l_water_livestock ) THEN
  IF ( ABS( rf_livestock - rmdi ) < EPSILON(1.0) ) THEN
    CALL ereport ( RoutineName, error_status,                                 &
                   "rf_livestock not found." )
  ELSE IF ( rf_livestock < 0.0 .OR. rf_livestock > 1.0 ) THEN
    CALL ereport ( RoutineName, error_status,                                 &
                   "rf_livestock must lie in the range 0 to 1." )
  END IF
END IF

IF ( l_water_industry ) THEN
  IF ( ABS( rf_industry - rmdi ) < EPSILON(1.0) ) THEN
    CALL ereport ( RoutineName, error_status,                                 &
                   "rf_industry not found." )
  ELSE IF ( rf_industry < 0.0 .OR. rf_industry > 1.0 ) THEN
    CALL ereport ( RoutineName, error_status,                                 &
                   "rf_industry must lie in the range 0 to 1." )
  END IF
END IF

!-----------------------------------------------------------------------------
! Check that the names provided for prioritisation are valid.
!-----------------------------------------------------------------------------
IF ( l_prioritise ) THEN
  ! Only consider the first nwater_use values.

  ! Check there are no duplicate names.
  DO i = 1,nwater_use-1
    IF ( ANY( priority(i+1:nwater_use) == priority(i) ) ) THEN
      CALL ereport ( RoutineName, error_status,                               &
                     "Duplicate priority: " // TRIM(priority(i)) )
    END IF
  END DO

  ! Check that names are valid.
  DO i = 1,nwater_use
    SELECT CASE ( priority(i) )
    CASE ( name_domestic, name_environment, name_industry, name_irrigation,   &
           name_livestock, name_transfers )
      ! These are valid names, nothing more to do.
    CASE ( 'xxx' )
      CALL ereport ( RoutineName, error_status,                               &
                     "Insufficient values provided for priority - " //        &
                     "each active sector must be listed." )
    CASE ( '' )
      CALL ereport ( RoutineName, error_status,                               &
                     "Priority name empty." )
    CASE DEFAULT
      CALL ereport ( RoutineName, error_status,                               &
                     "Priority name not valid: " // TRIM(priority(i)) )
    END SELECT
  END DO

  !---------------------------------------------------------------------------
  ! Check that all active sectors are given a priority, and inactive sectors
  ! are not given a priority.
  !---------------------------------------------------------------------------
  ! Domestic
  in_list = ANY( priority(1:nwater_use) == name_domestic )
  IF ( l_water_domestic ) THEN
    IF ( .NOT. in_list ) THEN
      CALL ereport ( RoutineName, error_status,                               &
                     "Domestic use is not given a priority." )
    END IF
  ELSE IF ( in_list ) THEN
    CALL ereport ( RoutineName, error_status,                                 &
                   "Domestic use is not modelled and so should not " //       &
                   "given a priority." )
  END IF

  ! Environmental
  in_list = ANY( priority(1:nwater_use) == name_environment )
  IF ( l_water_environment ) THEN
    IF ( .NOT. in_list ) THEN
      CALL ereport ( RoutineName, error_status,                               &
                     "Environmental use is not given a priority." )
    END IF
  ELSE IF ( in_list ) THEN
    CALL ereport ( RoutineName, error_status,                                 &
                   "Environmental use is not modelled and so should not " //  &
                   "be given a priority." )
  END IF

  ! Industrial
  in_list = ANY( priority(1:nwater_use) == name_industry )
  IF ( l_water_industry ) THEN
    IF ( .NOT. in_list ) THEN
      CALL ereport ( RoutineName, error_status,                               &
                     "Industrial use is not given a priority." )
    END IF
  ELSE IF ( in_list ) THEN
    CALL ereport ( RoutineName, error_status,                                 &
                   "Industrial use is not modelled and so should not " //     &
                   "be given a priority." )
  END IF

  ! Irrigation
  in_list = ANY( priority(1:nwater_use) == name_irrigation )
  IF ( l_water_irrigation ) THEN
    IF ( .NOT. in_list ) THEN
      CALL ereport ( RoutineName, error_status,                               &
                     "Irrigation use is not given a priority." )
    END IF
  ELSE IF ( in_list ) THEN
    CALL ereport ( RoutineName, error_status,                                 &
                   "Irrigation use is not modelled and so should not " //     &
                   "be given a priority." )
  END IF

  ! Livestock
  in_list = ANY( priority(1:nwater_use) == name_livestock )
  IF ( l_water_livestock ) THEN
    IF ( .NOT. in_list ) THEN
      CALL ereport ( RoutineName, error_status,                               &
                     "Livestock use is not given a priority." )
    END IF
  ELSE IF ( in_list ) THEN
    CALL ereport ( RoutineName, error_status,                                 &
                   "Livestock use is not modelled and so should not " //      &
                   "be given a priority." )
  END IF

  ! Transfers
  in_list = ANY( priority(1:nwater_use) == name_transfers )
  IF ( l_water_transfers ) THEN
    IF ( .NOT. in_list ) THEN
      CALL ereport ( RoutineName, error_status,                               &
                     "Water transfers are not given a priority." )
    END IF
  ELSE IF ( in_list ) THEN
    CALL ereport ( RoutineName, error_status,                                 &
                   "Water transfers are not modelled and so should not " //   &
                   "be given a priority." )
  END IF

END IF  !  l_prioritise

!-----------------------------------------------------------------------------
! Check option for non-renewable groundwater.
!-----------------------------------------------------------------------------
SELECT CASE ( nr_gwater_model )
CASE ( 0, nr_gwater_model_last, nr_gwater_model_use )
  ! Valid options, nothing more to do here.
CASE ( imdi )
  CALL ereport ( RoutineName, error_status,                                   &
                 "No value given for nr_gwater_model." )
CASE DEFAULT
  CALL ereport ( RoutineName, error_status,                                   &
                 "Invalid value for nr_gwater_model (non-renewable " //       &
                 "groundwater option)." )
END SELECT

!-----------------------------------------------------------------------------
! Ensure compatability with older irrigation options.
!-----------------------------------------------------------------------------
IF ( l_water_irrigation ) THEN
  ! Irrigation requires l_irrig_dmd to be TRUE (as that ensures that
  ! irrigation is considered in the surface flux and soil code).
  IF ( .NOT. l_irrig_dmd ) THEN
    CALL ereport ( RoutineName, error_status,                                 &
                   "l_water_irrigation also requires l_irrig_dmd." )
  END IF

  IF ( l_irrig_limit ) THEN
    CALL ereport ( RoutineName, error_status,                                 &
                   "l_water_irrigation requires l_irrig_limit=.false." )
  END IF
END IF

END SUBROUTINE check_jules_water_resources

!#############################################################################
!#############################################################################

SUBROUTINE set_jules_water_resources()

!-----------------------------------------------------------------------------
! Description:
!   Sets values related to water resource code. This should be called after
!   the namelist-processing routines.
!-----------------------------------------------------------------------------

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Local scalar parameters.
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER ::                                                &
   RoutineName = 'SET_JULES_WATER_RESOURCES'   ! Name of this procedure.

!-----------------------------------------------------------------------------
! Local scalar variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  n
    ! Counters.
   
!-----------------------------------------------------------------------------

IF ( l_water_resources ) THEN

  !---------------------------------------------------------------------------
  ! Set indices for each use to show position in multi-use arrays.
  !--------------------------------------------------------------------------- 
  n = 0
  IF ( l_water_domestic ) THEN
    n = n + 1
    use_domestic = n
  END IF

  IF ( l_water_environment ) THEN
    n = n + 1
    use_environment = n
  END IF

  IF ( l_water_industry ) THEN
    n = n + 1
    use_industry = n
  END IF

  IF ( l_water_irrigation ) THEN
    n = n + 1
    use_irrigation = n
  END IF

  IF ( l_water_livestock ) THEN
    n = n + 1
    use_livestock = n
  END IF

  IF ( l_water_transfers ) THEN
    n = n + 1
    use_transfers = n
  END IF

ELSE
  ! l_water_resources = .FALSE.

  ! Ensure the sector-specific switches are FALSE (namelists might have set
  ! them to TRUE), to simplify later checking.
  l_water_domestic    = .FALSE.
  l_water_environment = .FALSE.
  l_water_industry    = .FALSE.
  l_water_irrigation  = .FALSE.
  l_water_livestock   = .FALSE.
  l_water_transfers   = .FALSE.

END IF  !  l_water_resources 

END SUBROUTINE set_jules_water_resources

!#############################################################################
!#############################################################################

!#############################################################################
!#############################################################################

SUBROUTINE water_resources_alloc( land_pts )

! Allocate arrays.
! No USE statements other than error reporting and Dr Hook (if required).

IMPLICIT NONE

! Arguments.
INTEGER, INTENT(IN) :: land_pts
  ! Number of land points.

! Local parameters.
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'WATER_RESOURCES_ALLOC'

! Local scalar variables.
INTEGER :: land_pts_dim, nwater_use_dim

!-----------------------------------------------------------------------------
! Arrays are always allocated, but with minimal size if the science is not
! selected. Decide on sizes.
!-----------------------------------------------------------------------------
IF ( l_water_resources ) THEN
  land_pts_dim   = land_pts
  nwater_use_dim = nwater_use
ELSE
  land_pts_dim   = 1
  nwater_use_dim = 1
END IF

!-----------------------------------------------------------------------------
! Individual demands (which can be prescibed).
! We allocate a minimum size if a sector is not being used.
!-----------------------------------------------------------------------------
IF ( l_water_domestic ) THEN
  ALLOCATE( demand_rate_domestic(land_pts_dim) )
ELSE
  ALLOCATE( demand_rate_domestic(1) )
END IF

IF ( l_water_industry ) THEN
  ALLOCATE( demand_rate_industry(land_pts_dim) )
ELSE
  ALLOCATE( demand_rate_industry(1) )
END IF

IF ( l_water_livestock ) THEN
  ALLOCATE( demand_rate_livestock(land_pts_dim) )
ELSE
  ALLOCATE( demand_rate_livestock(1) )
END IF

IF ( l_water_transfers ) THEN
  ALLOCATE( demand_rate_transfers(land_pts_dim) )
ELSE
  ALLOCATE( demand_rate_transfers(1) )
END IF

demand_rate_domestic(:)  = 0.0
demand_rate_industry(:)  = 0.0
demand_rate_livestock(:) = 0.0
demand_rate_transfers(:) = 0.0

!-----------------------------------------------------------------------------
! Accumulated demands.
!-----------------------------------------------------------------------------
ALLOCATE( demand_accum(land_pts_dim,nwater_use_dim) )
demand_accum(:,:) = 0.0

!-----------------------------------------------------------------------------
! Ancillary fields.
!-----------------------------------------------------------------------------
ALLOCATE( conveyance_loss(land_pts_dim) )
IF ( l_water_irrigation ) THEN
  ALLOCATE( irrig_eff(land_pts_dim) )
ELSE
  ALLOCATE( irrig_eff(1) )
END IF
ALLOCATE( sfc_water_frac(land_pts_dim) )
conveyance_loss(:) = 0.0
irrig_eff(:)       = 0.0
sfc_water_frac(:)  = 0.0

!-----------------------------------------------------------------------------
! Other variables.
!-----------------------------------------------------------------------------
ALLOCATE( priority_order(land_pts_dim,nwater_use_dim) )
priority_order(:,:) = 0

RETURN
END SUBROUTINE water_resources_alloc

!#############################################################################

END MODULE jules_water_resources_mod
