! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

MODULE time_info_mod

! Information about the 360 day year comes from different places in the UM
! and standalone
! We unify them and expose them as a variable in this module
USE datetime_mod, ONLY: l_360, l_leap

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Module for unifying access to date / time information in the UM and
!   standalone
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------



CONTAINS


SUBROUTINE current_model_time(year, month, day, time)

! Note: Standalone JULES and the UM define their model time information in
! modules which have the same name (their content is different) - include the
! appropriate module here
USE model_time_mod, ONLY: current_time

IMPLICIT NONE

INTEGER, OPTIONAL, INTENT(OUT) :: year, month, day, time
    ! The time components as integers

!-----------------------------------------------------------------------------
! Description:
!   Gets the current model time as integer parts
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IF ( PRESENT(year) ) year = current_time%year
IF ( PRESENT(month) ) month = current_time%month
IF ( PRESENT(day) ) day = current_time%day
IF ( PRESENT(time) ) time = current_time%time

END SUBROUTINE current_model_time

END MODULE time_info_mod
