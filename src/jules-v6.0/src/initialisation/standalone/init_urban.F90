#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_urban(nml_dir)

!Use in variables
USE jules_surface_mod, ONLY: l_urban2t

USE switches_urban, ONLY: read_nml_jules_urban_switches,                      &
                          print_nlist_jules_urban_switches,                   &
                          check_jules_urban_switches

USE urban_param_mod, ONLY:  print_nlist_jules_urban2t_param

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
! Arguments
CHARACTER(LEN=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                         ! namelists

IF ( l_urban2t ) CALL read_nml_jules_urban_switches(nml_dir)
CALL print_nlist_jules_urban_switches()
CALL check_jules_urban_switches()
CALL print_nlist_jules_urban2t_param()

RETURN

END SUBROUTINE init_urban
#endif
