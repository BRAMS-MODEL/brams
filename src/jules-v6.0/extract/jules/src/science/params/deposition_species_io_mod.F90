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

MODULE deposition_species_io_mod

USE deposition_species_mod, ONLY: species_name_len
USE max_dimensions, ONLY: ndep_species_max, ntype_max

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Module containing variables used for reading deposition species
!   parameters, such as fixed-size versions of parameters, and an associated
!   namelist.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in BIOGENIC FLUXES.
!
! Code Description:
!   Language: Fortran 90.
!-----------------------------------------------------------------------------

! Public module variables

! Scalar variables.
REAL(KIND=real_jlslsm) ::                                                     &
  ch4_scaling_io,                                                             &
    ! Scaling applied to CH4 soil uptake (dimensionless).
    ! Originally this was used to match the value from the IPCC TAR.
  cuticle_o3_io,                                                              &
    ! Constant for calculation of cuticular resistance for ozone (s m-1).
  diffusion_coeff_io,                                                         &
    ! Diffusion coefficient (m2 s-1).
  diffusion_corr_io,                                                          &
    ! Diffusion correction for stomatal resistance, accounting for the
    ! different diffusivities of water and other species (dimensionless).
  r_tundra_io,                                                                &
    ! Surface resistance used in tundra region (s m-1).
  r_wet_soil_o3_io
    ! Wet soil surface resistance for ozone (s m-1).

CHARACTER(LEN=species_name_len) ::                                            &
  dep_species_name_io

! Array variables.
REAL(KIND=real_jlslsm) ::                                                     &
  ch4dd_tundra_io(4),                                                         &
    ! Coefficients of cubic polynomial relating CH4 loss for tundra to
    ! temperature.
  dd_ice_coeff_io(3),                                                         &
    ! Coefficients in quadratic function relating dry deposition over ice to
    ! temperature.
  h2dd_c_io(ntype_max),                                                       &
    ! Constant in quadratic function relating hydrogen deposition to soil
    ! moisture, for each surface type (s m-1).
  h2dd_m_io(ntype_max),                                                       &
    ! Coefficient of first order term in quadratic function relating hydrogen
    ! deposition to soil moisture, for each surface type (s m-1).
  h2dd_q_io(ntype_max),                                                       &
    ! Coefficient of second order term in quadratic function relating
    ! hydrogen deposition to soil moisture, for each surface type (s m-1).
  rsurf_std_io(ntype_max)
    ! Surface resistance used in tundra region for each surface type (s m-1).

!-----------------------------------------------------------------------------
! Namelist.
!-----------------------------------------------------------------------------
NAMELIST / jules_deposition_species /                                         &
  ch4_scaling_io, ch4dd_tundra_io, cuticle_o3_io,                             &
  dd_ice_coeff_io, dep_species_name_io, diffusion_coeff_io,                   &
  diffusion_corr_io, h2dd_c_io, h2dd_m_io, h2dd_q_io,                         &
  r_tundra_io, r_wet_soil_o3_io, rsurf_std_io

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DEPOSITION_SPECIES_IO'

CONTAINS

!#############################################################################

SUBROUTINE print_nlist_jules_deposition_species()
USE jules_print_mgr, ONLY: jules_print

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER ::                                                &
  RoutineName = 'PRINT_NLIST_JULES_DEPOSITION_SPECIES'
    ! Name of this procedure.

CHARACTER(LEN=50000) :: lineBuffer

CALL jules_print(RoutineName,                                                 &
                 'Contents of namelist jules_deposition_species')

WRITE(lineBuffer,"(2A)")' dep_species_name_io = ',dep_species_name_io
CALL jules_print(RoutineName,lineBuffer)

WRITE(lineBuffer,"(A,G11.4E2)")' ch4_scaling = ',ch4_scaling_io
CALL jules_print(RoutineName,lineBuffer)

WRITE(lineBuffer,*)' ch4dd_tundra_io = ',ch4dd_tundra_io(:)
CALL jules_print(RoutineName,lineBuffer)

WRITE(lineBuffer,"(A,G11.4E2)")' cuticle_o3_io = ',cuticle_o3_io
CALL jules_print(RoutineName,lineBuffer)

WRITE(lineBuffer,"(A,G11.4E2)")' r_wet_soil_o3_io = ',r_wet_soil_o3_io
CALL jules_print(RoutineName,lineBuffer)

WRITE(lineBuffer,"(A,G11.4E2)")' diffusion_coeff_io = ',diffusion_coeff_io
CALL jules_print(RoutineName,lineBuffer)

WRITE(lineBuffer,"(A,G11.4E2)")' diffusion_corr_io = ',diffusion_corr_io
CALL jules_print(RoutineName,lineBuffer)

WRITE(lineBuffer,*)' dd_ice_coeff_io = ',dd_ice_coeff_io(:)
CALL jules_print(RoutineName,lineBuffer)

WRITE(lineBuffer,*)' h2dd_c_io = ',h2dd_c_io(:)
CALL jules_print(RoutineName,lineBuffer)

WRITE(lineBuffer,*)' h2dd_m_io = ',h2dd_c_io(:)
CALL jules_print(RoutineName,lineBuffer)

WRITE(lineBuffer,*)' h2dd_q_io = ',h2dd_c_io(:)
CALL jules_print(RoutineName,lineBuffer)

WRITE(lineBuffer,*)' rsurf_std_io = ',rsurf_std_io(:)
CALL jules_print(RoutineName,lineBuffer)

CALL jules_print(RoutineName,                                                 &
                 '- - - - - - end of namelist - - - - - -')

END SUBROUTINE print_nlist_jules_deposition_species

!#############################################################################

END MODULE deposition_species_io_mod
#endif
