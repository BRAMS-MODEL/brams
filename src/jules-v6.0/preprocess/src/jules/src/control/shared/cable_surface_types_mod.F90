MODULE cable_surface_types_mod

!-----------------------------------------------------------------------------
! Description:
!   Contains CABLE surface type information and a namelist for setting them
!   The CABLE equivalent of jules_surface_types_mod
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE 
!-----------------------------------------------------------------------------

USE missing_data_mod, ONLY: imdi

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  nnpft_cable,                                                                &
                ! Number of natural pfts
                !   Derived from namelist inputs as npft_cable - ncpft_cable
  ncpft_cable = 0,                                                            &
                ! Number of crop pfts
  npft_cable  = imdi,                                                         &
                ! Number of plant functional types
  nnvg_cable  = imdi,                                                         &
                ! Number of non-vegetation surface types
  ntype_cable,                                                                &
                ! Total number of surface types
  lakes_cable  = imdi,                                                        &
                ! Index of the CABLE lakes surface type
  ice_cable  = imdi,                                                          &
                ! Index of the CABLE ice surface type
  urban_cable  = imdi,                                                        &
                ! Index of the CABLE urban surface type
  barren_cable = imdi
                ! Index of the barren surface type

!-----------------------------------------------------------------------------
! UM only: The following are only required for UM flexible tiles.
!-----------------------------------------------------------------------------
! The UM uses the following to identify the original PFTs in the UM output.

INTEGER ::                                                                    &
   evergreen_needleleaf = imdi,                                               &
                  ! Index of surface type 'evergreen needleleaf'
   evergreen_broadleaf  = imdi,                                               &
                  ! Index of surface type 'evergreen broadleaf'
   deciduous_needleleaf = imdi,                                               &
                  ! Index of surface type 'deciduous needleleaf'
   deciduous_broadleaf  = imdi,                                               &
                  ! Index of surface type 'deciduous broadleaf'
   shrub_cable          = imdi,                                               &
                  ! Index of surface type 'shrub'
   c3_grassland = imdi,                                                       &
                  ! Index of surface type 'C3 grassland'
   c4_grassland = imdi,                                                       &
                  ! Index of surface type 'C4 grassland'
   tundra = imdi,                                                             &
                  ! Index of surface type 'tundra'
   c3_cropland = imdi,                                                        &
                  ! Index of surface type 'C3 cropland'
   c4_cropland = imdi,                                                        &
                  ! Index of surface type 'C4 cropland'
   wetland = imdi,                                                            &
                  ! Index of surface type 'wetland'
   empty1 = imdi,                                                             &
                  ! Index of surface type 'empty 1'
   empty2 = imdi 
                  ! Index of surface type 'empty 2'


!-----------------------------------------------------------------------------
! Single namelist definition for UM and standalone
!-----------------------------------------------------------------------------
NAMELIST  / cable_surface_types/                                              &
    npft_cable, ncpft_cable, nnvg_cable,                                      &
    evergreen_needleleaf, evergreen_broadleaf, deciduous_needleleaf,          &
    deciduous_broadleaf, shrub_cable, c3_grassland, c4_grassland,             &
    tundra, c3_cropland, c4_cropland, wetland, empty1, empty2,                &
    barren_cable, urban_cable, lakes_cable, ice_cable

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CABLE_SURFACE_TYPES_MOD'

CONTAINS

SUBROUTINE check_cable_surface_types()

USE max_dimensions, ONLY: npft_max, ncpft_max, nnvg_max

USE ereport_mod, ONLY: ereport

!-----------------------------------------------------------------------------
! Description:
!   Checks CABLE_SURFACE_TYPES namelist for consistency
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE 
!-----------------------------------------------------------------------------

IMPLICIT NONE

INTEGER :: i ! Loop counter

INTEGER :: errorstatus

!-----------------------------------------------------------------------------
! Check that the given values are less than the fixed values for IO
!-----------------------------------------------------------------------------
IF ( npft_cable > npft_max ) THEN
  errorstatus = 101
  CALL ereport("check_jules_surface_types", errorstatus,                      &
               "npft_cable > npft_max - increase npft_max and recompile")
END IF
IF ( ncpft_cable > ncpft_max ) THEN
  errorstatus = 101
  CALL ereport("check_jules_surface_types", errorstatus,                      &
               "ncpft_cable > ncpft_max - increase ncpft_max and recompile")
END IF
IF ( nnvg_cable > nnvg_max ) THEN
  errorstatus = 101
  CALL ereport("check_jules_surface_types", errorstatus,                      &
               "nnvg_cable > nnvg_max - increase nnvg_max and recompile")
END IF
IF ( ncpft_cable > npft_cable ) THEN
  errorstatus = 101
  CALL ereport("check_jules_surface_types", errorstatus,                      &
               "ncpft_cable > npft_cable - total number of PFTs must be >= "  &
               // "number of crop PFTs")
END IF

! If these are true, then we automatically have:
!   ntype_cable <= ntype_max (since ntype_cable(_max) = npft_cable(_max) +
!   nnvg_cable(_max))
!   nnpft_cable <= nnpft_max (since nnpft_cable(_max) = npft_cable(_max) -
!   ncpft_cable(_max))

!-----------------------------------------------------------------------------
! Check values for the specific surface types are sensible
!-----------------------------------------------------------------------------
! PFT surface types must come before non-veg types, so if urban_cable,
! lakes_cable or ice_cable are given (i.e. > 0) then they must be > npft_cable
! A soil type is required
IF ( urban_cable > 0 .AND.                                                    &
     ( urban_cable <= npft_cable .OR. urban_cable > ntype_cable ) ) THEN
  errorstatus = 101
  CALL ereport("check_cable_surface_types", errorstatus,                      &
               "urban_cable tile is given but is out of range")
END IF

IF ( lakes_cable > 0 .AND.                                                    &
     ( lakes_cable <= npft_cable .OR. lakes_cable > ntype_cable ) ) THEN
  errorstatus = 101
  CALL ereport("check_cable_surface_types", errorstatus,                      &
               "lakes_cable tile is given but is out of range")
END IF

IF ( ice_cable > 0 .AND. ( ice_cable <= npft_cable .OR. ice_cable > ntype_cable ) ) THEN
  errorstatus = 101
  CALL ereport("check_cable_surface_types", errorstatus,                      &
               "ice_cable tile is given but is out of range")
END IF

END SUBROUTINE check_cable_surface_types


SUBROUTINE print_nlist_cable_surface_types()

USE jules_print_mgr, ONLY: jules_print

IMPLICIT NONE

INTEGER :: i ! Loop counter

CHARACTER(LEN=50000) :: lineBuffer


!-----------------------------------------------------------------------------

CALL jules_print('cable_surface_types',                                       &
                 'Contents of namelist cable_surface_types')

WRITE(lineBuffer, *) '  npft_cable = ', npft_cable
CALL jules_print('jules_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  nnvg_cable = ', nnvg_cable
CALL jules_print('jules_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  evergreen_needleleaf = ', evergreen_needleleaf
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  evergreen_broadleaf = ', evergreen_broadleaf
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  deciduous_needleleaf = ', deciduous_broadleaf
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  shrub_cable = ', shrub_cable
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  c3_grassland = ', c3_grassland
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  c4_grassland = ', c4_grassland
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  tundra = ', tundra
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  c3_cropland = ', c3_cropland
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  c4_cropland = ', c4_cropland
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  wetland = ', wetland
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  empty1 = ', empty1
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  empty2 = ', empty2
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  barren_cable = ', barren_cable
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  urban_cable = ', urban_cable
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  lakes_cable = ', lakes_cable
CALL jules_print('cable_surface_types', lineBuffer)

WRITE(lineBuffer, *) '  ice_cable = ', ice_cable
CALL jules_print('cable_surface_types', lineBuffer)

CALL jules_print('cable_surface_types',                                       &
    '- - - - - - end of namelist - - - - - -')

END SUBROUTINE print_nlist_cable_surface_types


SUBROUTINE set_derived_variables_cable_surface_types()

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Derive ntype_cable and nnpft_cable from the namelist values
!-----------------------------------------------------------------------------
ntype_cable = npft_cable + nnvg_cable
nnpft_cable = npft_cable - ncpft_cable

RETURN

END SUBROUTINE set_derived_variables_cable_surface_types

END MODULE cable_surface_types_mod
