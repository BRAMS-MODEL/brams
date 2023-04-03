#if !defined(UM_JULES)
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

SUBROUTINE init_deposition_species(nml_dir)

USE deposition_species_io_mod, ONLY:                                          &
! imported scalars
  ch4_scaling_io, cuticle_o3_io, dd_ice_coeff_io,                             &
  dep_species_name_io, diffusion_coeff_io, diffusion_corr_io,                 &
  r_tundra_io, r_wet_soil_o3_io,                                              &
! imported arrays
  ch4dd_tundra_io, h2dd_c_io, h2dd_m_io, h2dd_q_io, rsurf_std_io,             &
! imported namelist
  jules_deposition_species

USE deposition_species_mod, ONLY:                                             &
! imported procedures
  check_jules_deposition_species,                                             &
! imported scalars
  ch4_scaling, cuticle_o3, r_wet_soil_o3,                                     &
! imported arrays
  ch4dd_tundra, dd_ice_coeff, dep_species_name, diffusion_coeff,              &
  diffusion_corr, h2dd_c, h2dd_m, h2dd_q, r_tundra, rsurf_std

USE jules_deposition_mod, ONLY:                                               &
! imported scalars
  l_deposition, ndry_dep_species

USE jules_surface_types_mod, ONLY:                                            &
! imported scalars
  ntype

USE io_constants, ONLY: namelist_unit
USE missing_data_mod, ONLY:                                                   &
! imported scalar parameters
  rmdi
USE string_utils_mod, ONLY: to_string
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises parameters related to atmospheric deposition.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Scalar arguments with INTENT(IN).
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                         ! namelists

!-----------------------------------------------------------------------------
! Parameters.
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_DEPOSITION_SPECIES'

!-----------------------------------------------------------------------------
! Scalar variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  error,                                                                      &
    ! Error indicator
  i,                                                                          &
    ! Loop counter.
  j
    ! Loop counter.

CHARACTER(LEN=errormessagelength) :: iomessage

!-----------------------------------------------------------------------------
! There's nothing to do if deposition is not selected.
!-----------------------------------------------------------------------------
IF ( .NOT. l_deposition ) RETURN

!-----------------------------------------------------------------------------
! Read deposition species namelists, one per species.
!-----------------------------------------------------------------------------
CALL log_info(RoutineName, "Reading JULES_DEPOSITION namelist...")

! Open the namelist file.
OPEN(namelist_unit, FILE=(TRIM(nml_dir) // '/' // 'jules_deposition.nml'),    &
     STATUS='old', POSITION='rewind', ACTION='read', IOSTAT = error,          &
     IOMSG = iomessage)

IF ( error /= 0 ) THEN
  CALL log_fatal(RoutineName,                                                 &
                 "Error opening namelist file jules_deposition.nml " //       &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")
END IF

!-----------------------------------------------------------------------------
! Loop over the species.
!-----------------------------------------------------------------------------
DO i = 1,ndry_dep_species

  ! Initialise namelist values before reading the next species.
  dep_species_name_io = 'unset'
  ch4_scaling_io      = rmdi
  ch4dd_tundra_io(:)  = rmdi
  cuticle_o3_io       = rmdi
  dd_ice_coeff_io(:)  = rmdi
  diffusion_coeff_io  = rmdi
  diffusion_corr_io   = rmdi
  h2dd_c_io(:)        = rmdi
  h2dd_m_io(:)        = rmdi
  h2dd_q_io(:)        = rmdi
  r_tundra_io         = rmdi
  r_wet_soil_o3_io    = rmdi
  rsurf_std_io(:)     = rmdi

  ! Read namelist for this species.
  READ( namelist_unit, NML = jules_deposition_species, IOSTAT = error,        &
       IOMSG = iomessage)

  IF ( error /= 0 ) THEN
    CALL log_fatal(RoutineName,                                               &
                   "Error reading namelist JULES_DEPOSITION_SPECIES " //      &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")
  END IF

  !---------------------------------------------------------------------------
  ! Load values into the final variables.
  !---------------------------------------------------------------------------
  ! Variables that have a value for every species (even if not necessarily
  ! used for all species).
  dep_species_name(i) = dep_species_name_io 
  dd_ice_coeff(i,:)   = dd_ice_coeff_io(:)
  diffusion_coeff(i)  = diffusion_coeff_io 
  diffusion_corr(i)   = diffusion_corr_io 
  r_tundra(i)         = r_tundra_io 
  rsurf_std(:,i)      = rsurf_std_io(1:ntype)

  ! Deal with variables that only apply to one species.
  ! If values are given for other species they will be ignored.
  SELECT CASE ( dep_species_name(i) )

  CASE ( 'CH4' )

    ch4_scaling     = ch4_scaling_io
    ch4dd_tundra(:) = ch4dd_tundra_io(:)

  CASE ( 'H2' )
    h2dd_c(:) = h2dd_c_io(1:ntype)
    h2dd_m(:) = h2dd_m_io(1:ntype)
    h2dd_q(:) = h2dd_q_io(1:ntype)

  CASE ( 'O3' )
    cuticle_o3    = cuticle_o3_io
    r_wet_soil_o3 = r_wet_soil_o3_io

  END SELECT

END DO  !  species

!-----------------------------------------------------------------------------
! Close the namelist file
!-----------------------------------------------------------------------------
CLOSE(namelist_unit, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 ) THEN
  CALL log_fatal(RoutineName,                                                 &
                 "Error closing namelist file jules_deposition.nml " //       &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")
END IF

!-----------------------------------------------------------------------------
! Check the values.
!-----------------------------------------------------------------------------
CALL check_jules_deposition_species()

!-----------------------------------------------------------------------------
! Print information about the parameters.
!-----------------------------------------------------------------------------
DO i = 1,ndry_dep_species
  CALL log_info( RoutineName, "Species: " // TRIM(dep_species_name(i)) )
  CALL log_info( RoutineName, "  diffusion_coeff = " //                       &
                              to_string(diffusion_coeff(i)) )
  CALL log_info( RoutineName, "  diffusion_corr = " //                        &
                              to_string(diffusion_corr(i)) )
  CALL log_info( RoutineName, "  r_tundra = " //                              &
                              to_string(r_tundra(i)) )

  SELECT CASE ( dep_species_name(i) )
  CASE ( 'CH4' )
    CALL log_info( RoutineName, "  ch4_scaling = " //                         &
                                to_string(ch4_scaling) )
  CASE ( 'O3' )
    CALL log_info( RoutineName, "  r_wet_soil_o3 = " //                       &
                                to_string(r_wet_soil_o3) )
    CALL log_info( RoutineName, "  cuticle_o3 = " //                          &
                                to_string(cuticle_o3) )
  END SELECT
  
END DO

RETURN

END SUBROUTINE init_deposition_species

#endif
