! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE jules_hydrology_mod

!-----------------------------------------------------------------------------
! Description:
!   Contains hydrology options and a namelist for setting them
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Switches
!-----------------------------------------------------------------------------
LOGICAL ::                                                                    &
  l_hydrology = .TRUE.,                                                       &
      ! Turns off hydrology code for UM_JULES
  l_top  = .FALSE.,                                                           &
      ! Switch for TOPMODEL-based hydrology
  l_pdm = .FALSE.,                                                            &
      ! Switch for PDM hydrology
  l_spdmvar = .FALSE.,                                                        &
      ! Switch for slope dependent s_pdm in PDM hydrology
  l_baseflow_corr = .TRUE.,                                                   &
      ! Switch for using a correction to the calculation of baseflow
      ! Only used if l_top = T
  l_var_rainfrac = .FALSE.,                                                   &
      ! Switch for using convective and large scale rain fractions
      ! as calculated in the UM atmosphere
      ! Has no effect in standalone JULES
  l_wetland_unfrozen = .FALSE.,                                               &
      ! Switch for TOPMODEL-based hydrology with unfrozen wetland inundation
      ! Only used if l_top=.T.
  l_limit_gsoil = .FALSE.
      ! Switch for limiting gsoil above theta_crit

!-----------------------------------------------------------------------------
! PDM parameters
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                     &
  dz_pdm = 1.0,                                                               &
      ! Soil layer thickness for PDM (m)
  b_pdm  = 1.0,                                                               &
      ! Shape factor for PDM
  s_pdm = 0.0,                                                                &
! So/Smax factor for PDM
   slope_pdm_max = 6.0
! Maximum topographic slope (deg) in the slope dependent s_pdm
! linear function
!-----------------------------------------------------------------------------
! TOPMODEL parameters
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  nfita = 20
      ! Number of loops for fitting.
REAL(KIND=real_jlslsm) ::                                                     &
  ti_max = 10.0,                                                              &
      ! Maximum topographic index considered
  ti_wetl = 1.5,                                                              &
      ! Parameter to remove very high water tables from the calculated
      ! wetland fraction
  zw_max = 6.0
      ! Maximum allowed water table depth (m)
!------------------------------------------------------------------------------
! Single namelist definition for UM and standalone
!------------------------------------------------------------------------------
NAMELIST  / jules_hydrology/                                                  &
  l_hydrology, l_top, l_pdm, l_spdmvar, l_baseflow_corr, l_var_rainfrac,      &
  l_wetland_unfrozen, l_limit_gsoil,                                          &
  dz_pdm, b_pdm, s_pdm, slope_pdm_max, ti_max, ti_wetl, zw_max, nfita


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='JULES_HYDROLOGY_MOD'

CONTAINS

SUBROUTINE check_jules_hydrology()

USE ereport_mod, ONLY: ereport

!-----------------------------------------------------------------------------
! Description:
!   Checks JULES_HYDROLOGY namelist for consistency
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IMPLICIT NONE

INTEGER :: errorstatus

! PDM and TOPMODEL are not allowed together
IF ( l_top .AND. l_pdm ) THEN
  errorstatus = 101
  CALL ereport("check_jules_hydrology", errorstatus,                          &
               "PDM and TOPMODEL cannot be used together")
END IF

END SUBROUTINE check_jules_hydrology


SUBROUTINE print_nlist_jules_hydrology()

USE jules_print_mgr, ONLY: jules_print

IMPLICIT NONE

CHARACTER(LEN=50000) :: lineBuffer

CALL jules_print('jules_hydrology',                                           &
                 'Contents of namelist jules_hydrology')

WRITE(lineBuffer, *) '  l_top = ', l_top
CALL jules_print('jules_hydrology', lineBuffer)

WRITE(lineBuffer, *) '  l_pdm = ', l_pdm
CALL jules_print('jules_hydrology', lineBuffer)

WRITE(lineBuffer, *) '  l_spdmvar = ', l_spdmvar
CALL jules_print('jules_hydrology', lineBuffer)

WRITE(lineBuffer, *) '  l_baseflow_corr = ', l_baseflow_corr
CALL jules_print('jules_hydrology', lineBuffer)

WRITE(lineBuffer, *) '  l_var_rainfrac = ', l_var_rainfrac
CALL jules_print('jules_hydrology', lineBuffer)

WRITE(lineBuffer, *) '  l_wetland_unfrozen = ', l_wetland_unfrozen
CALL jules_print('jules_hydrology', lineBuffer)

WRITE(lineBuffer, *) '  l_limit_gsoil = ', l_limit_gsoil
CALL jules_print('jules_hydrology', lineBuffer)

WRITE(lineBuffer, *) '  dz_pdm = ', dz_pdm
CALL jules_print('jules_hydrology', lineBuffer)

WRITE(lineBuffer, *) '  b_pdm = ', b_pdm
CALL jules_print('jules_hydrology', lineBuffer)

WRITE(lineBuffer, *) '  s_pdm = ', s_pdm
CALL jules_print('jules_hydrology', lineBuffer)

WRITE(lineBuffer, *) '  slope_pdm_max = ', s_pdm
CALL jules_print('jules_hydrology', lineBuffer)

WRITE(lineBuffer, *) '  ti_max = ', ti_max
CALL jules_print('jules_hydrology', lineBuffer)

WRITE(lineBuffer, *) '  ti_wetl = ', ti_wetl
CALL jules_print('jules_hydrology', lineBuffer)

WRITE(lineBuffer, *) '  zw_max = ', zw_max
CALL jules_print('jules_hydrology', lineBuffer)

WRITE(lineBuffer, *) '  nfita = ', nfita
CALL jules_print('jules_hydrology', lineBuffer)

CALL jules_print('jules_hydrology',                                           &
    '- - - - - - end of namelist - - - - - -')

END SUBROUTINE print_nlist_jules_hydrology


END MODULE jules_hydrology_mod
