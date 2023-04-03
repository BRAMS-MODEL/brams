#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE internal_init_profile_vars(profile)

USE io_constants, ONLY: max_dim_var

USE data_cube_mod, ONLY: cube_create, cube_free

USE datetime_mod, ONLY: datetime_advance

USE model_interface_mod, ONLY: get_var_levs_dims

USE model_time_mod, ONLY: is_spinup, spinup_start, timestep_len

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises fields required to start sampling data for the current section
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(output_profile), INTENT(INOUT) :: profile  ! The profile to initialise
                                                ! fields for


! Work variables
INTEGER :: grid_shape(2)  ! The size of the two grid dimensions (depending
                          ! on whether a subgrid is used)

INTEGER :: var_ndims   ! The number of levels dimensions the variable has
INTEGER :: lev_dim_sizes(max_dim_var)  ! The sizes of the levels dimensions

TYPE(datetime) :: section_start  ! The start time of this 'section' of output

INTEGER :: i  ! Index variable


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Allocate the data cubes for the output fields
!-----------------------------------------------------------------------------
! Work out what shape the grid is
IF ( use_subgrid ) THEN
  grid_shape(:) = (/ subgrid%nx, subgrid%ny /)
ELSE
  grid_shape(:) = (/ grid%nx, grid%ny /)
END IF

DO i = 1,profile%nfields
  ! Get the number and sizes of the levels dimensions for the variable
  CALL get_var_levs_dims(                                                     &
    profile%fields(i)%var_id, ndims = var_ndims, dim_sizes = lev_dim_sizes    &
  )

  ! Deallocate any existing cube first
  CALL cube_free(profile%fields(i)%field_data)

  ! Allocate a cube of the correct shape for the variable
  profile%fields(i)%field_data =                                              &
                   cube_create((/ grid_shape, lev_dim_sizes(1:var_ndims) /))
END DO


!-----------------------------------------------------------------------------
! Initialise the profile so that the sampling and output periods start at
! the correct time
!-----------------------------------------------------------------------------
! Work out the start time for the current section of output
IF ( is_spinup ) THEN
  section_start = spinup_start
ELSE
  section_start = profile%output_start
END IF

! First sample will be at time=sample_period after start of section. We
! reduce by one timestep so as to match the time at the start of the
! timestep.
profile%next_sample_time    = datetime_advance( section_start,                &
                                                profile%sample_period         &
                                                - timestep_len )
profile%current_output_time = section_start
profile%next_output_time    = datetime_advance(section_start, profile%output_period)
profile%samples_in_period   = 0

RETURN

END SUBROUTINE internal_init_profile_vars
#endif
