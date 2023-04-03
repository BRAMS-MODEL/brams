! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE jules_rivers_mod

!-----------------------------------------------------------------------------
! Description:
!   Contains river routing options and a namelist for setting them
!   This currently holds both "physics" and control variables
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

USE missing_data_mod, ONLY: imdi, rmdi

USE um_types, ONLY: real_jlslsm

USE jules_irrig_mod, ONLY: l_irrig_limit

USE ereport_mod, ONLY: ereport

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Scalar parameters.
INTEGER, PARAMETER ::                                                         &
   rivers_um_trip = 1
                            ! value of i_river_vn indicating that the UM
                            ! TRIP linear model is to be used.
                            ! TRIP really only refers to the river network
                            ! data, but is used here to refer to a linear
                            ! model that has often been used with the data.
INTEGER, PARAMETER ::                                                         &
   rivers_rfm = 2
                            ! value of i_river_vn indicating that the RFM
                            ! Kinematic Wave model is to be used.
INTEGER, PARAMETER ::                                                         &
   rivers_trip = 3
                            ! value of i_river_vn indicating that the TRIP
                            ! linear model is to be used.
                            ! TRIP really only refers to the river network
                            ! data, but is used here to refer to a linear
                            ! model that has often been used with the data.

!-----------------------------------------------------------------------------
! Array parameters.

INTEGER,PARAMETER :: flow_dir_delta(0:10,2) = RESHAPE( (/                     &
!      missing  N    NE    E    SE    S    SW     W   NW   no_dir no_dir
!  x offsets:
         -9,    0,   1,    1,    1,   0,   -1,   -1,  -1,    0,   0,          &
!  y offsets:
         -9,    1,   1,    0,   -1,  -1,   -1,    0,   1,    0,   0           &
         /), (/ 11,2 /) )

!      Components of displacement (number of gridboxes) in x and
!      y directions to the immediately downstream gridbox. These
!      correspond to the directions in flow_dir_set.
!      The elements in the 2nd dimension are x and y respectively.
!             e.g. flow_dir_delta(1,1:2)=(/0,1/)
!             flow_dir_set(1) refers to a displacement to the N
!             flow_dir_delta(1,2)=1 means 1 gridbox in the N (y) direction
!             flow_dir_delta(1,1)=0 means no displacement in the W-E (x)
!                                 direction
!     The zeroth column (flow_dir_delta(0,:)) is not currently used, but is
!     included so as to make size(1) of flow_dir_delta and flow_dir_set equal
!     in an attempt for clarity.
!     The "no defined direction" must be the last column.

! Note on flow direction assumptions:
!       values of flow direction code that represent the
!       displacement given by the corresponding positions in flow_dir_delta.
!       The following are the meanings of each element of the array:
!        0: no data value (e.g. over sea)
!        1-8: N,NE,E,SE,S,SW,W,NW  i.e. clockwise from N
!             Although referred to via these compass directions, they are used
!             as "grid-relative" directions, i.e. it is assumed that columns
!             run S-N, and rows W-E, so "N" means "same column, one row up".
!             If the grid is actually rotated (so that columns do not run S-N
!             on the earth), the point that is "same column, one row up" in
!             fact does not lie immediately N.
!        9: undefined flow direction (i.e. no outflow, or outflow to sea)
!             For some encoding schemes, there is not a single value that
!             represents undefined flow (e.g. ARC combines values of all
!             directions with equal slope - so "undefined" flow direction
!             depends upon slope to neighbouring points). This is not a problem
!             here, since none of these schemes is currently encoded, BUT may
!             have to be accounted for in future - code could be added below
!             where currently we stop with "unexpected value".
!        Note that the values of flow_dir_set that correspond to "real" flow
!        directions (as opposed to missing data or no defined flow direction)
!        must be >0. Further, ALL values of flow_dir_set (incl sea and no flow)
!        must be >=0. In that case, values <0 can be used to indicate points
!        that have a flow direction that points off the edge of the grid (for
!        non-global applications).  These restrictions are necessary to make
!        the current algorithm work. In particular, a run with a smaller grid
!        (e.g. regional) should give the same routing pathways as a run with a
!        larger grid (at the points that are on both grids).

!-----------------------------------------------------------------------------
! Items set in namelist
!-----------------------------------------------------------------------------

LOGICAL ::                                                                    &
   l_rivers = .FALSE.                                                         &
                            ! Switch for runoff routing
   ,l_inland = .FALSE.
                            ! Control rerouting of inland basin water

INTEGER ::                                                                    &
   nstep_rivers = imdi                                                        &
                            ! Timestep for runoff routing
                            ! (number of model timesteps)
   ,i_river_vn  = imdi                                                        &
                            ! integer representation of river routing type
                            !   1 == 'um_trip'
                            !   2 == 'rfm'
                            !   3 == 'trip'
   ,a_thresh = 1
                            ! threshold area (TRIP: pixels; RFM: km)

REAL(KIND=real_jlslsm) ::                                                     &
   rivers_meander = 1.4                                                       &
                            ! meander ratio for rivers - the ratio of
                            ! actual river length to the calculated length
   ,rivers_speed = 0.4                                                        &
                            ! flow speed for rivers (m s-1)
   ,runoff_factor = 1.0                                                       &
                            ! runoff volume factor
   ,cland  = 0.2                                                              &
                            ! land wave speed (m/s)
   ,criver = 0.62                                                             &
                            ! subsurf river wave speed (m/s)
   ,cbland  = 0.1                                                             &
                            ! subsurf land wave speed (m/s)
   ,cbriver = 0.15                                                            &
                            ! subsurf river wave speed (m/s)
   ,retl = 0.0                                                                &
                            ! return flow (land squares) (<1)
   ,retr = 0.005
                            ! return flow (river squares) (<1)

!-----------------------------------------------------------------------------
! Rivers parameters
!-----------------------------------------------------------------------------

! Scalar variables (general)
INTEGER ::                                                                    &
   np_rivers = 0                                                              &
                            ! number of points in the rivers grid at which
                            ! routing is calculated
   ,nx_rivers = 0                                                             &
                            ! row length for rivers grid
   ,ny_rivers = 0
                            ! column length for rivers grid

!-----------------------------------------------------------------------------
! Definition of the river routing grid
!-----------------------------------------------------------------------------

REAL(KIND=real_jlslsm) ::                                                     &
   rivers_dlat  = rmdi                                                        &
                            ! size of gridbox of (regular) rivers grid
                            ! in latitude (degrees)
   ,rivers_dlon = rmdi                                                        &
                            ! size of gridbox of (regular) rivers grid
                            ! in longitude (degrees)
   ,rivers_lat1 = rmdi                                                        &
                            ! latitude of southernmost row of gridpoints
                            ! on a regular rivers grid (degrees)
   ,rivers_lon1 = rmdi                                                        &
                            !  longitude of westernmost (first) column of
                            ! gridpoints on a regular rivers grid (degrees)
   ,rivers_dx   = rmdi
                            ! size of gridbox of rivers grid in m (for
                            ! non-regular lat/lon grids)

INTEGER ::                                                                    &
   nx_grid      = imdi                                                        &
                            ! row length of full land model grid
                            ! (only needed for river routing)
   ,ny_grid     = imdi                                                        &
                            ! column length of model grid
                            ! (only needed for river routing)
   ,nseqmax     = imdi
                            ! maximum value of routing grid sequence

REAL(KIND=real_jlslsm) ::                                                     &
   reg_lon1     = rmdi                                                        &
                            ! longitude of westernnmost row of gridpoints
                            ! on a regular full model grid (degrees)
   ,reg_lat1    = rmdi                                                        &
                            ! latitude of southernnmost row of gridpoints
                            ! on a regular full model grid (degrees)
   ,reg_dlon    = rmdi                                                        &
                            ! size of gridbox of (regular) full model grid
                            ! in longitude (degrees)
   ,reg_dlat    = rmdi
                            ! size of gridbox of (regular) full model grid
                            ! in latitude (degrees)

LOGICAL ::                                                                    &
   rivers_reglatlon = .TRUE.                                                  &
                            ! flag indicating if rivers grid is regular in
                            ! latitude and longitude See above for
                            ! definition of a regular grid.
   ,rivers_regrid   = .TRUE.                                                  &
                            ! flag indicating if model and rivers grids
                            ! are identical
                            !     FALSE grids are identical
                            !     TRUE grids differ and regridding required
   ,rivers_first    = .TRUE.                                                  &
                            ! TRUE indicates first river rivers timestep
   ,rivers_call     = .FALSE.

!-----------------------------------------------------------------------------
! Array variables defined on land points required as input for river routing
!-----------------------------------------------------------------------------

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  tot_surf_runoff_gb(:),                                                      &
                           !  accumulated surface runoff (production)
                           !  between calls to rivers (kg m-2 s-1)
  tot_sub_runoff_gb(:),                                                       &
                           !  accumulatd sub-surface runoff (production)
                           !  rate between calls to rivers (kg m-2 s-1)
  acc_lake_evap_gb(:),                                                        &
                           !  accumulated lake evap over river routing
                           !  timestep (Kg/m2) - on land points
  rivers_sto_per_m2_on_landpts(:),                                            &
                           ! Water storage (kg m-2) on land points
  rivers_adj_on_landpts(:)
                           ! adjustment factor for water storage on landpts

INTEGER, ALLOCATABLE ::                                                       &
  il_river_grid(:),                                                           &
                           ! map of land point index on river routing grid
  ir_land_grid(:)
                           ! map of river point index on land model grid

!-----------------------------------------------------------------------------
! Array variables defined on river routing grid used in routing calculations.
!-----------------------------------------------------------------------------

! Arrays defined on rivers points only

INTEGER, ALLOCATABLE ::                                                       &
   rivers_index_rp(:)                                                         &
                            ! Index of points where routing is calculated
   ,rivers_next_rp(:)
                            ! Index of the next downstream point.

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
   rivers_seq_rp(:),                                                          &
                            ! River routing pathway sequence
   rivers_dir_rp(:),                                                          &
                            ! River routing direction index
   rivers_sto_rp(:),                                                          &
                            ! Water storage (kg)
   rivers_dra_rp(:),                                                          &
                            ! Catchment area draining to a grid cell
                            !    (no. of grid cells)
   rivers_lat_rp(:),                                                          &
                            ! River routing point latitude
   rivers_lon_rp(:)
                            ! River routing point longitude

! Ancillary arrays defined on full 2D rivers grid
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
   rivers_seq(:,:),                                                           &
                            ! River routing pathway sequence
   rivers_dir(:,:),                                                           &
                            ! River routing direction index
   rivers_dra(:,:),                                                           &
                            ! Catchment area draining to a grid cell
                            !    (no. of grid cells)
   rivers_lat2d(:,:),                                                         &
                            ! Full 2D latitude field
                            ! (enables non-regular lat-lon river grids)
   rivers_lon2d(:,:),                                                         &
                            ! Full 2D longitude field
                            ! (enables non-regular lat-lon river grids)
   rivers_xgrid(:),                                                           &
                            ! 1D x-dimension of rivers grid
   rivers_ygrid(:)
                            ! 1D y-dimension of rivers grid

! Arrays defined on 1d river points vectors
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
   rfm_flowobs1_rp(:),                                                        &
                            ! Initial (observed) river flow (kg m-2 s-1)
   rfm_surfstore_rp(:),                                                       &
                            ! Surface storage (m3)
   rfm_substore_rp(:),                                                        &
                            ! Sub-surface storage (m3)
   rfm_flowin_rp(:),                                                          &
                            ! Surface lateral inflow (m3)
   rfm_bflowin_rp(:),                                                         &
                            ! Sub-surface lateral inflow (m3)
   rfm_rivflow_rp(:),                                                         &
                            ! Surface river flow (m3 s-1)
   rfm_baseflow_rp(:)
                            ! Sub-surface flow (m3 s-1)

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
   rivers_boxareas_rp(:),                                                     &
                            ! Gridbox area of each river grid pixel (m2)
   rfm_iarea_rp(:),                                                           &
                            ! Number of pixels draining to a pixel
   rfm_land_rp(:)
                            ! Flag to indicate river grid pixel type
                            !    0 = sea, 1 = river, 2 = land
                            ! (n.b. declared REAL not INTEGER)

!-----------------------------------------------------------------------------
! Single namelist definition for UM and standalone
!-----------------------------------------------------------------------------

NAMELIST  / jules_rivers/                                                     &
  l_rivers, l_inland, i_river_vn, nstep_rivers,                               &
  cland, criver, cbland, cbriver, runoff_factor, retl, retr,                  &
  a_thresh, rivers_meander, rivers_speed

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='JULES_RIVERS_MOD'

!-----------------------------------------------------------------------------

CONTAINS

SUBROUTINE jules_rivers_alloc(land_pts)

!No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
INTEGER, INTENT(IN) :: land_pts

!Local variables
INTEGER :: temp_size, temp_tiles, temp_layers
INTEGER :: error, error_sum  ! Error indicators

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='JULES_RIVERS_ALLOC'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE( tot_surf_runoff_gb(land_pts))
ALLOCATE( tot_sub_runoff_gb(land_pts))
ALLOCATE( acc_lake_evap_gb(land_pts))

tot_surf_runoff_gb(:) = 0.0
tot_sub_runoff_gb(:)  = 0.0
acc_lake_evap_gb(:)   = 0.0

IF ( l_irrig_limit ) THEN
  temp_size = land_pts
ELSE
  temp_size = 1
END IF

ALLOCATE( rivers_sto_per_m2_on_landpts(temp_size), stat = error )
error_sum = error
ALLOCATE( rivers_adj_on_landpts(temp_size), stat = error )
error_sum = error_sum + error
IF ( error_sum /= 0 ) THEN
  CALL ereport( RoutineName, error_sum,                                       &
                'Error allocating for l_irrig_limit.' )
END IF
rivers_sto_per_m2_on_landpts(:) = 0.0
rivers_adj_on_landpts(:)        = 1.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE jules_rivers_alloc

#if !defined(LFRIC)
SUBROUTINE check_jules_rivers()

USE ereport_mod, ONLY: ereport

USE jules_irrig_mod, ONLY: l_irrig_limit

!-----------------------------------------------------------------------------
! Description:
!   Checks JULES_RIVERS namelist for consistency
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IMPLICIT NONE

INTEGER :: errcode

IF ( .NOT. l_rivers ) THEN
  IF (l_irrig_limit) THEN
    errcode = 101
    CALL ereport("check_jules_rivers", errcode,                               &
                 'l_irrig_limit=T requires l_rivers=T ')
  END IF

  ! If rivers are not enabled, there is nothing further to check
  RETURN
END IF

! Check that a valid integer timestep was given
IF ( nstep_rivers <= 0 ) THEN
  errcode = 101
  CALL ereport("check_jules_rivers", errcode, 'nstep_rivers must be > 0')
END IF

! Check that parameter values are appropriate for the selected algorithm
! This also serves as a check that we have a recognised rivers type
SELECT CASE ( i_river_vn )
CASE ( rivers_rfm )
  IF ( cland <= 0.0 .OR. criver <= 0.0 ) THEN
    errcode = 101
    CALL ereport("check_jules_rivers", errcode,                               &
                 "Surface wave speeds must be > 0")
  END IF
  IF ( cbland <= 0.0 .OR. cbriver <= 0.0) THEN
    errcode = 102
    CALL ereport("check_jules_rivers", errcode,                               &
                 "Sub surface wave speeds must be > 0")
  END IF
  IF ( runoff_factor <= 0.0 ) THEN
    errcode = 103
    CALL ereport("check_jules_rivers", errcode,                               &
                 "Runoff factor must be > 0")
  END IF

CASE ( rivers_trip )
  IF ( rivers_speed <= 0.0 .OR. rivers_meander <= 0.0 ) THEN
    errcode = 104
    CALL ereport("check_jules_rivers", errcode,                               &
                 "River speed and meander ratio must be > 0")
  END IF
CASE ( rivers_um_trip )
  IF ( rivers_speed <= 0.0 .OR. rivers_meander <= 0.0 ) THEN
    errcode = 104
    CALL ereport("check_jules_rivers", errcode,                               &
                 "River speed and meander ratio must be > 0")
  END IF
CASE DEFAULT
  errcode = 101
  CALL ereport("check_jules_rivers", errcode,                                 &
               'Unrecognised river routing algorithm (i_river_vn)')
END SELECT

IF ( l_irrig_limit .AND. ( i_river_vn /= rivers_trip ) ) THEN
  errcode = 101
  CALL ereport("check_jules_rivers", errcode,                                 &
               'l_irrig_limit=T requires i_river_vn=3')
END IF

END SUBROUTINE check_jules_rivers
#endif

SUBROUTINE print_nlist_jules_rivers()

USE jules_print_mgr, ONLY: jules_print

IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_JULES_RIVERS'

CALL jules_print('jules_rivers_inputs_mod',                                   &
   'Contents of namelist jules_rivers')

WRITE(lineBuffer,*)' l_rivers = ',l_rivers
CALL jules_print('jules_rivers',lineBuffer)
WRITE(lineBuffer,*)' l_inland = ',l_inland
CALL jules_print('jules_rivers',lineBuffer)
WRITE(lineBuffer,*)' i_river_vn = ',i_river_vn
CALL jules_print('jules_rivers',lineBuffer)
WRITE(lineBuffer,*)' nstep_rivers = ',nstep_rivers
CALL jules_print('jules_rivers',lineBuffer)

SELECT CASE ( i_river_vn )

CASE ( rivers_trip, rivers_um_trip )
  WRITE(lineBuffer,*)' RIVERS_SPEED = ',rivers_speed
  CALL jules_print('jules_rivers',lineBuffer)
  WRITE(lineBuffer,*)' RIVERS_MEANDER = ',rivers_meander
  CALL jules_print('jules_rivers',lineBuffer)

CASE ( rivers_rfm )
  WRITE(lineBuffer,*)' CLAND = ',cland
  CALL jules_print('jules_rivers',lineBuffer)
  WRITE(lineBuffer,*)' CRIVER = ',criver
  CALL jules_print('jules_rivers',lineBuffer)
  WRITE(lineBuffer,*)' CBLAND = ',cbland
  CALL jules_print('jules_rivers',lineBuffer)
  WRITE(lineBuffer,*)' CBRIVER = ',cbriver
  CALL jules_print('jules_rivers',lineBuffer)
  WRITE(lineBuffer,*)' RETL = ',retl
  CALL jules_print('jules_rivers',lineBuffer)
  WRITE(lineBuffer,*)' RETR = ',retr
  CALL jules_print('jules_rivers',lineBuffer)
  WRITE(lineBuffer,*)' A_THRESH = ',a_thresh
  CALL jules_print('jules_rivers',lineBuffer)
  WRITE(lineBuffer,*)' RUNOFF_FACTOR = ',runoff_factor
  CALL jules_print('jules_rivers',lineBuffer)

END SELECT

CALL jules_print('jules_rivers_mod',                                          &
   '- - - - - - end of namelist - - - - - -')

END SUBROUTINE print_nlist_jules_rivers


#if defined(UM_JULES) && !defined(LFRIC)
SUBROUTINE read_nml_jules_rivers(unit_in)

USE setup_namelist,   ONLY: setup_nml_type
USE check_iostat_mod, ONLY: check_iostat
USE UM_parcore,       ONLY: mype
USE parkind1,         ONLY: jprb, jpim
USE yomhook,          ONLY: lhook, dr_hook
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER,INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus, errcode
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_JULES_RIVERS'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 3
INTEGER, PARAMETER :: n_real = 9
INTEGER, PARAMETER :: n_log = 2

TYPE my_namelist
  SEQUENCE
  INTEGER :: nstep_rivers
  INTEGER :: a_thresh
  INTEGER :: i_river_vn
  REAL(KIND=real_jlslsm) :: cland
  REAL(KIND=real_jlslsm) :: criver
  REAL(KIND=real_jlslsm) :: cbland
  REAL(KIND=real_jlslsm) :: cbriver
  REAL(KIND=real_jlslsm) :: runoff_factor
  REAL(KIND=real_jlslsm) :: retl
  REAL(KIND=real_jlslsm) :: retr
  REAL(KIND=real_jlslsm) :: rivers_meander
  REAL(KIND=real_jlslsm) :: rivers_speed
  LOGICAL :: l_rivers
  LOGICAL :: l_inland
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in = n_int,              &
                    n_real_in = n_real, n_log_in = n_log)

IF (mype == 0) THEN
  READ(UNIT = unit_in, NML = jules_rivers, IOSTAT = ErrorStatus,              &
                                           IOMSG = iomessage)
  CALL check_iostat(errorstatus, "namelist JULES_RIVERS", iomessage)

  my_nml % nstep_rivers = nstep_rivers
  my_nml % a_thresh = a_thresh
  my_nml % i_river_vn = i_river_vn
  my_nml % cland = cland
  my_nml % criver = criver
  my_nml % cbland = cbland
  my_nml % cbriver = cbriver
  my_nml % runoff_factor = runoff_factor
  my_nml % retl = retl
  my_nml % retr = retr
  my_nml % rivers_meander = rivers_meander
  my_nml % rivers_speed = rivers_speed
  my_nml % l_rivers = l_rivers
  my_nml % l_inland = l_inland
END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN
  nstep_rivers = my_nml % nstep_rivers
  a_thresh = my_nml % a_thresh
  i_river_vn = my_nml % i_river_vn
  cland = my_nml % cland
  criver = my_nml % criver
  cbland = my_nml % cbland
  cbriver = my_nml % cbriver
  runoff_factor = my_nml % runoff_factor
  retl = my_nml % retl
  retr = my_nml % retr
  a_thresh = my_nml % a_thresh
  rivers_meander = my_nml % rivers_meander
  rivers_speed = my_nml % rivers_speed
  l_rivers = my_nml % l_rivers
  l_inland = my_nml % l_inland
END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_jules_rivers
#endif

END MODULE jules_rivers_mod

!-----------------------------------------------------------------------------
