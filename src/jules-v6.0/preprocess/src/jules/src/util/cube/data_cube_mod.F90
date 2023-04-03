
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

MODULE data_cube_mod

USE io_constants, ONLY: mdi

USE logging_mod, ONLY: log_fatal

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Type definitions
!-----------------------------------------------------------------------------
TYPE data_cube
  !-----------------------------------------------------------------------------
  ! This type encapsulates information about a cube of data
  !-----------------------------------------------------------------------------
  INTEGER, POINTER :: SHAPE(:) => NULL()  ! The size of each dimension of
                                          ! the cube
  REAL, POINTER :: values(:) => NULL()  ! The data contained in the cube

END TYPE data_cube

!-----------------------------------------------------------------------------
! Overloads
!-----------------------------------------------------------------------------
! We want to use a single name to convert arrays of different dimensionalities
! to cubes
INTERFACE cube_from_array
MODULE PROCEDURE cube_from_array_1d,     cube_from_array_2d,                  &
                 cube_from_array_3d,     cube_from_array_4d,                  &
                 cube_from_array_5d,     cube_from_array_6d,                  &
                 cube_from_array_7d
END INTERFACE cube_from_array

! We want to use a single name to read data from the cube in different
! dimensionalities
INTERFACE cube_get_data
MODULE PROCEDURE cube_get_data_scalar, cube_get_data_1d,                      &
                 cube_get_data_2d,     cube_get_data_3d,                      &
                 cube_get_data_4d,     cube_get_data_5d,                      &
                 cube_get_data_6d,     cube_get_data_7d
END INTERFACE cube_get_data

! Operators
INTERFACE operator(+)
MODULE PROCEDURE cube_add_scalar, cube_add_cube
END INTERFACE

INTERFACE operator(-)
MODULE PROCEDURE cube_sub_scalar, cube_sub_cube
END INTERFACE

INTERFACE operator( * )
MODULE PROCEDURE cube_mul_scalar, cube_mul_cube
END INTERFACE

INTERFACE operator(/)
MODULE PROCEDURE cube_div_scalar, cube_div_cube
END INTERFACE


!-----------------------------------------------------------------------------
! Visibility declarations
!-----------------------------------------------------------------------------
PRIVATE
PUBLIC                                                                        &
! Public types
    data_cube,                                                                &
! Public operators
    operator(+), operator(-), operator( * ), operator(/),                     &
! Public procedures
    cube_create, cube_safe_copy, cube_from_array, cube_get_data, cube_free,   &
    cube_min, cube_max


CONTAINS


! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

FUNCTION cube_create(SHAPE) RESULT(cube)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Create a data cube
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
INTEGER, INTENT(IN) :: SHAPE(:)  ! The shape of the cube

! Return type
TYPE(data_cube) :: cube  ! The (empty) cube


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Allocate and populate the cube
!-----------------------------------------------------------------------------
ALLOCATE(cube%SHAPE(SIZE(SHAPE)))
ALLOCATE(cube%values(PRODUCT(SHAPE)))

cube%SHAPE(:)  = SHAPE(:)
cube%values(:) = 0.0  ! Initialise values to 0

RETURN

END FUNCTION cube_create
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

SUBROUTINE cube_safe_copy(c1, c2)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Copies cube c2 to cube c1 in a memory safe way
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(data_cube), INTENT(INOUT) :: c1  ! The cube to copy to
TYPE(data_cube), INTENT(IN) :: c2     ! The cube to copy


!-----------------------------------------------------------------------------


! First, make sure the cube we are copying to is deallocated
CALL cube_free(c1)

! Create a new cube of the correct shape
c1 = cube_create(c2%SHAPE)

! Copy the values across
c1%values(:) = c2%values(:)

RETURN

END SUBROUTINE cube_safe_copy
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

FUNCTION cube_from_array_1d(values) RESULT(cube)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Creates a cube with the same shape and values as the given array
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
REAL, INTENT(IN) :: values(:)  ! The array to create the cube from


! Return type
TYPE(data_cube) :: cube  ! The cube containing the data


!-----------------------------------------------------------------------------


cube = cube_create(SHAPE(values))
cube%values(:) = RESHAPE(values, (/ SIZE(values) /))

RETURN

END FUNCTION cube_from_array_1d


FUNCTION cube_from_array_2d(values) RESULT(cube)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Creates a cube with the same shape and values as the given array
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
REAL, INTENT(IN) :: values(:,:)  ! The array to create the cube from


! Return type
TYPE(data_cube) :: cube  ! The cube containing the data


!-----------------------------------------------------------------------------


cube = cube_create(SHAPE(values))
cube%values(:) = RESHAPE(values, (/ SIZE(values) /))

RETURN

END FUNCTION cube_from_array_2d


FUNCTION cube_from_array_3d(values) RESULT(cube)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Creates a cube with the same shape and values as the given array
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
REAL, INTENT(IN) :: values(:,:,:)  ! The array to create the cube from


! Return type
TYPE(data_cube) :: cube  ! The cube containing the data


!-----------------------------------------------------------------------------


cube = cube_create(SHAPE(values))
cube%values(:) = RESHAPE(values, (/ SIZE(values) /))

RETURN

END FUNCTION cube_from_array_3d


FUNCTION cube_from_array_4d(values) RESULT(cube)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Creates a cube with the same shape and values as the given array
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
REAL, INTENT(IN) :: values(:,:,:,:)  ! The array to create the cube from


! Return type
TYPE(data_cube) :: cube  ! The cube containing the data


!-----------------------------------------------------------------------------


cube = cube_create(SHAPE(values))
cube%values(:) = RESHAPE(values, (/ SIZE(values) /))

RETURN

END FUNCTION cube_from_array_4d


FUNCTION cube_from_array_5d(values) RESULT(cube)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Creates a cube with the same shape and values as the given array
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
REAL, INTENT(IN) :: values(:,:,:,:,:)  ! The array to create the cube from


! Return type
TYPE(data_cube) :: cube  ! The cube containing the data


!-----------------------------------------------------------------------------


cube = cube_create(SHAPE(values))
cube%values(:) = RESHAPE(values, (/ SIZE(values) /))

RETURN

END FUNCTION cube_from_array_5d


FUNCTION cube_from_array_6d(values) RESULT(cube)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Creates a cube with the same shape and values as the given array
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
REAL, INTENT(IN) :: values(:,:,:,:,:,:)  ! The array to create the cube from


! Return type
TYPE(data_cube) :: cube  ! The cube containing the data


!-----------------------------------------------------------------------------


cube = cube_create(SHAPE(values))
cube%values(:) = RESHAPE(values, (/ SIZE(values) /))

RETURN

END FUNCTION cube_from_array_6d


FUNCTION cube_from_array_7d(values) RESULT(cube)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Creates a cube with the same shape and values as the given array
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
REAL, INTENT(IN) :: values(:,:,:,:,:,:,:)  ! The array to create the cube from


! Return type
TYPE(data_cube) :: cube  ! The cube containing the data


!-----------------------------------------------------------------------------


cube = cube_create(SHAPE(values))
cube%values(:) = RESHAPE(values, (/ SIZE(values) /))

RETURN

END FUNCTION cube_from_array_7d
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*****************************************************************************
! Note that these operators are all 'mdi-propagating', i.e. any values set
! to mdi in incoming cubes will propagate into outgoing cubes as mdi
!*****************************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


PURE LOGICAL FUNCTION is_mdi(x)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Checks if a value is (epsilon close to) mdi
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 95.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
REAL, INTENT(IN) :: x

REAL, PARAMETER :: small_val = EPSILON(1.0)


!-----------------------------------------------------------------------------


! Return TRUE if the result is epsilon close to mdi, otherwise return false
is_mdi = ABS(x - mdi) < small_val

RETURN

END FUNCTION is_mdi


!*****************************************************************************


FUNCTION cube_add_scalar(c1, add) RESULT(c2)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes a cube and adds a scalar to it's contents
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 95.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(data_cube), INTENT(IN) :: c1  ! The cube to add to
REAL, INTENT(IN) :: add  ! The scalar to add


! Return type
TYPE(data_cube) :: c2  ! The cube containing the summed data


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Create a new cube of the correct shape
!-----------------------------------------------------------------------------
c2 = cube_create(c1%SHAPE)

!-----------------------------------------------------------------------------
! Set the data for the outgoing cube in a way that preserves mdis from the
! incoming cubes
!-----------------------------------------------------------------------------
c2%values = mdi_safe_add(c1%values, add)

RETURN

END FUNCTION cube_add_scalar


FUNCTION cube_add_cube(c1, c2) RESULT(c3)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes two cubes and adds their contents together element-wise
!   If the cubes do not have the same shape, an error is thrown
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 95.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(data_cube), INTENT(IN) :: c1, c2  ! The cubes to add


! Return type
TYPE(data_cube) :: c3  ! The cube containing the summed data


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that the two cubes to add have compatible shapes
!-----------------------------------------------------------------------------
IF ( SIZE(c1%SHAPE) /= SIZE(c2%SHAPE) .OR. .NOT. ALL(c1%SHAPE == c2%SHAPE) )  &
  CALL log_fatal("cube_add_cube",                                             &
                 "Error adding cubes - shapes do not match")

!-----------------------------------------------------------------------------
! Create a cube of the correct shape
!-----------------------------------------------------------------------------
c3 = cube_create(c1%SHAPE)

!-----------------------------------------------------------------------------
! Set the data for the outgoing cube in a way that preserves mdis from the
! incoming cubes
!-----------------------------------------------------------------------------
c3%values = mdi_safe_add(c1%values, c2%values)

RETURN

END FUNCTION cube_add_cube


ELEMENTAL FUNCTION mdi_safe_add(x, y) RESULT(z)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes two values and adds them together in a mdi-preserving way
!   This function is declared as ELEMENTAL, meaning that it can operate on
!   arrays
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 95.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
REAL, INTENT(IN) :: x, y  ! The values to add together

! Return type
REAL :: z


!-----------------------------------------------------------------------------


! If either of the incoming values is mdi, then the result is mdi
IF ( is_mdi(x) .OR. is_mdi(y) ) THEN
  z = mdi
  RETURN
END IF

! Otherwise, just add the values together to produce the result
z = x + y

RETURN

END FUNCTION mdi_safe_add


!*****************************************************************************


FUNCTION cube_sub_scalar(c1, sub) RESULT(c2)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes a cube and subtracts a scalar from it's contents
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 95.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(data_cube), INTENT(IN) :: c1  ! The cube to subtract from
REAL, INTENT(IN) :: sub  ! The scalar to subtract


! Return type
TYPE(data_cube) :: c2  ! The cube containing the subtracted data


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Create a new cube of the correct shape
!-----------------------------------------------------------------------------
c2 = cube_create(c1%SHAPE)

!-----------------------------------------------------------------------------
! Set the data for the outgoing cube in a way that preserves mdis from the
! incoming cubes
!-----------------------------------------------------------------------------
c2%values = mdi_safe_sub(c1%values, sub)

RETURN

END FUNCTION cube_sub_scalar


FUNCTION cube_sub_cube(c1, c2) RESULT(c3)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes two cubes and subtracts the contents of the second from the first,
!   element-wise
!   If the cubes do not have the same shape, an error is thrown
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 95.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(data_cube), INTENT(IN) :: c1, c2  ! The cubes to add


! Return type
TYPE(data_cube) :: c3  ! The cube containing the summed data


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that the two cubes to add have compatible shapes
!-----------------------------------------------------------------------------
IF ( SIZE(c1%SHAPE) /= SIZE(c2%SHAPE) .OR. .NOT. ALL(c1%SHAPE == c2%SHAPE) )  &
  CALL log_fatal("cube_sub_cube",                                             &
                 "Error adding cubes - shapes do not match")

!-----------------------------------------------------------------------------
! Create a cube of the correct shape
!-----------------------------------------------------------------------------
c3 = cube_create(c1%SHAPE)

!-----------------------------------------------------------------------------
! Set the data for the outgoing cube in a way that preserves mdis from the
! incoming cubes
!-----------------------------------------------------------------------------
c3%values = mdi_safe_sub(c1%values, c2%values)

RETURN

END FUNCTION cube_sub_cube


ELEMENTAL FUNCTION mdi_safe_sub(x, y) RESULT(z)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes two values x and y and subtracts y from x in a mdi-preserving way
!   This function is declared as ELEMENTAL, meaning that it can operate on
!   arrays
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 95.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
REAL, INTENT(IN) :: x, y  ! The values to subtract

! Return type
REAL :: z


!-----------------------------------------------------------------------------


! If either of the incoming values is mdi, then the result is mdi
IF ( is_mdi(x) .OR. is_mdi(y) ) THEN
  z = mdi
  RETURN
END IF

! Otherwise, just subtract y from x to produce the result
z = x - y

RETURN

END FUNCTION mdi_safe_sub


!*****************************************************************************


FUNCTION cube_mul_scalar(c1, mult) RESULT(c2)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes a cube and multiplies its contents by a scalar
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 95.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(data_cube), INTENT(IN) :: c1  ! The cube to multiply
REAL, INTENT(IN) :: mult  ! The scalar to multiply by


! Return type
TYPE(data_cube) :: c2  ! The cube containing the multiplied data


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Create a cube of the correct shape
!-----------------------------------------------------------------------------
c2 = cube_create(c1%SHAPE)

!-----------------------------------------------------------------------------
! Set the data for the outgoing cube in a way that preserves mdis from the
! incoming cubes
!-----------------------------------------------------------------------------
c2%values = mdi_safe_mul(c1%values, mult)

RETURN

END FUNCTION cube_mul_scalar


FUNCTION cube_mul_cube(c1, c2) RESULT(c3)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes two cubes and multiplies their contents together element-wise
!   If the cubes do not have the same shape, an error is thrown
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 95.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(data_cube), INTENT(IN) :: c1, c2  ! The cubes to multiply


! Return type
TYPE(data_cube) :: c3  ! The cube containing the multiplied data


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that the two cubes to add have compatible shapes
!-----------------------------------------------------------------------------
IF ( SIZE(c1%SHAPE) /= SIZE(c2%SHAPE) .OR. .NOT. ALL(c1%SHAPE == c2%SHAPE) )  &
  CALL log_fatal("cube_mul_cube",                                             &
                 "Error multiplying cubes - shapes do not match")

!-----------------------------------------------------------------------------
! Create a cube of the correct shape
!-----------------------------------------------------------------------------
c3 = cube_create(c1%SHAPE)

!-----------------------------------------------------------------------------
! Set the data for the outgoing cube in a way that preserves mdis from the
! incoming cubes
!-----------------------------------------------------------------------------
c3%values = mdi_safe_mul(c1%values, c2%values)

RETURN

END FUNCTION cube_mul_cube


ELEMENTAL FUNCTION mdi_safe_mul(x, y) RESULT(z)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes two values and multiplies them in a mdi-preserving way
!   This function is declared as ELEMENTAL, meaning that it can operate on
!   arrays
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 95.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
REAL, INTENT(IN) :: x, y  ! The values to multiply

! Return type
REAL :: z


!-----------------------------------------------------------------------------


! If either of the incoming values is mdi, then the result is mdi
IF ( is_mdi(x) .OR. is_mdi(y) ) THEN
  z = mdi
  RETURN
END IF

! Otherwise, just subtract y from x to produce the result
z = x * y

RETURN

END FUNCTION mdi_safe_mul


!*****************************************************************************


FUNCTION cube_div_scalar(c1, divisor) RESULT(c2)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes a cube and divides its contents by a scalar
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 95.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(data_cube), INTENT(IN) :: c1  ! The cube to divide
REAL, INTENT(IN) :: divisor  ! The scalar to divide by


! Return type
TYPE(data_cube) :: c2  ! The cube containing the divided data


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Create a cube of the correct shape
!-----------------------------------------------------------------------------
c2 = cube_create(c1%SHAPE)

!-----------------------------------------------------------------------------
! Set the data for the outgoing cube in a way that preserves mdis from the
! incoming cubes
!-----------------------------------------------------------------------------
c2%values = mdi_safe_div(c1%values, divisor)

RETURN

END FUNCTION cube_div_scalar


FUNCTION cube_div_cube(c1, c2) RESULT(c3)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes two cubes and divides the contents of the first by the contents of
!   the second, element-wise
!   Where the second cube has near-zero values, the resulting cube will contain
!   zeroes
!   If the cubes do not have the same shape, an error is thrown
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 95.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(data_cube), INTENT(IN) :: c1, c2  ! The cubes to divide


! Return type
TYPE(data_cube) :: c3  ! The cube containing the divided data


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Create a cube of the correct shape
!-----------------------------------------------------------------------------
IF ( SIZE(c1%SHAPE) /= SIZE(c2%SHAPE) .OR. .NOT. ALL(c1%SHAPE == c2%SHAPE) )  &
  CALL log_fatal("cube_div_cube",                                             &
                 "Error dividing cubes - shapes do not match")

!-----------------------------------------------------------------------------
! Create a cube containing divided data
!-----------------------------------------------------------------------------
c3 = cube_create(c1%SHAPE)

!-----------------------------------------------------------------------------
! Set the data for the outgoing cube in a way that preserves mdis from the
! incoming cubes
!-----------------------------------------------------------------------------
c3%values = mdi_safe_div(c1%values, c2%values)

RETURN

END FUNCTION cube_div_cube


ELEMENTAL FUNCTION mdi_safe_div(x, y) RESULT(z)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes two values x and y and divides x by y in a mdi-preserving way
!   This function is declared as ELEMENTAL, meaning that it can operate on
!   arrays
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 95.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
REAL, INTENT(IN) :: x, y  ! The values to divide

! Return type
REAL :: z


!-----------------------------------------------------------------------------


! If either of the incoming values is mdi, then the result is mdi
IF ( is_mdi(x) .OR. is_mdi(y) ) THEN
  z = mdi
  RETURN
END IF

! If the denominator is 0, return 0
IF ( ABS(y) < EPSILON(1.0) ) THEN
  z = 0
  RETURN
END IF

! Otherwise, just divide x by y to produce the result
z = x / y

RETURN

END FUNCTION mdi_safe_div
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*****************************************************************************
! Note that these operators are all 'mdi-propagating', i.e. any values set
! to mdi in incoming cubes will propagate into outgoing cubes as mdi
!*****************************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

FUNCTION cube_min(c1, c2) RESULT(c3)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes two cubes and finds the minimum values element-wise
!   If the cubes do not have the same shape, an error is thrown
!
! Code Author: Karina Williams
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 95.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(data_cube), INTENT(IN) :: c1, c2  ! The cubes to compare


! Return type
TYPE(data_cube) :: c3  ! The cube containing the minimum values


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that the two input cubes have compatible shapes
!-----------------------------------------------------------------------------
IF ( SIZE(c1%SHAPE) /= SIZE(c2%SHAPE) .OR. .NOT. ALL(c1%SHAPE == c2%SHAPE) )  &
  CALL log_fatal("min_of_cube_and_cube",                                      &
                 "Error finding min of cubes - shapes do not match")

!-----------------------------------------------------------------------------
! Create a cube of the correct shape
!-----------------------------------------------------------------------------
c3 = cube_create(c1%SHAPE)

!-----------------------------------------------------------------------------
! Set the data for the outgoing cube in a way that preserves mdis from the
! incoming cubes
!-----------------------------------------------------------------------------
c3%values = mdi_safe_min(c1%values, c2%values)

RETURN

END FUNCTION cube_min


!*****************************************************************************


ELEMENTAL FUNCTION mdi_safe_min(x, y) RESULT(z)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes two values and finds minimum in a mdi-preserving way
!   This function is declared as ELEMENTAL, meaning that it can operate on
!   arrays
!
! Code Author: Karina Williams
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 95.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
REAL, INTENT(IN) :: x, y  ! The values to find the minimum of

! Return type
REAL :: z


!-----------------------------------------------------------------------------


! If either of the incoming values is mdi, then the result is mdi
IF ( is_mdi(x) .OR. is_mdi(y) ) THEN
  z = mdi
  RETURN
END IF

! Otherwise, just find the minimum
z = MIN(x, y)

RETURN

END FUNCTION mdi_safe_min


!-----------------------------------------------------------------------------


FUNCTION cube_max(c1, c2) RESULT(c3)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes two cubes and finds the maximum values element-wise
!   If the cubes do not have the same shape, an error is thrown
!
! Code Author: Karina Williams
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 95.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(data_cube), INTENT(IN) :: c1, c2  ! The cubes to compare


! Return type
TYPE(data_cube) :: c3  ! The cube containing the maximum values


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that the two input cubes have compatible shapes
!-----------------------------------------------------------------------------
IF ( SIZE(c1%SHAPE) /= SIZE(c2%SHAPE) .OR. .NOT. ALL(c1%SHAPE == c2%SHAPE) )  &
  CALL log_fatal("max_of_cube_and_cube",                                      &
                 "Error finding max of cubes - shapes do not match")

!-----------------------------------------------------------------------------
! Create a cube of the correct shape
!-----------------------------------------------------------------------------
c3 = cube_create(c1%SHAPE)

!-----------------------------------------------------------------------------
! Set the data for the outgoing cube in a way that preserves mdis from the
! incoming cubes
!-----------------------------------------------------------------------------
c3%values = mdi_safe_max(c1%values, c2%values)

RETURN

END FUNCTION cube_max


!*****************************************************************************


ELEMENTAL FUNCTION mdi_safe_max(x, y) RESULT(z)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes two values and finds maximum in a mdi-preserving way
!   This function is declared as ELEMENTAL, meaning that it can operate on
!   arrays
!
! Code Author: Karina Williams
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 95.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
REAL, INTENT(IN) :: x, y  ! The values to find the maximum of

! Return type
REAL :: z


!-----------------------------------------------------------------------------


! If either of the incoming values is mdi, then the result is mdi
IF ( is_mdi(x) .OR. is_mdi(y) ) THEN
  z = mdi
  RETURN
END IF

! Otherwise, just find the maximum
z = MAX(x, y)

RETURN

END FUNCTION mdi_safe_max

! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

SUBROUTINE cube_get_data_scalar(cube, VALUE)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Get the data from a cube where the data is a scalar
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(data_cube), INTENT(IN) :: cube  ! The cube containing the data

REAL, INTENT(OUT) :: VALUE  ! The scalar value to retrieve


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that the data in the cube is of the correct size to extract a scalar
!-----------------------------------------------------------------------------
IF ( SIZE(cube%values) /= 1 )                                                 &
  CALL log_fatal("cube_get_data_scalar",                                      &
                 "Cube must have size 1 to extract a scalar")

!-----------------------------------------------------------------------------
! Extract the scalar and return it
!-----------------------------------------------------------------------------
VALUE = cube%values(1)

RETURN

END SUBROUTINE cube_get_data_scalar


SUBROUTINE cube_get_data_1d(cube, values)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Get the data from a cube where the data is a 1d array
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(data_cube), INTENT(IN) :: cube  ! The cube containing the data

REAL, INTENT(OUT) :: values(:)  ! The array to put values in


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that values is the right shape for the data in the cube
!-----------------------------------------------------------------------------
! First make sure that cube%shape has the correct number of dimensions
IF ( SIZE(cube%SHAPE) /= 1 )                                                  &
  CALL log_fatal("cube_get_data_1d",                                          &
                 "values has a different rank to the cube")

! Check that values is the correct size
IF ( SIZE(values) /= cube%SHAPE(1) )                                          &
  CALL log_fatal("cube_get_data_1d",                                          &
                 "values is not the same size as the cube")

!-----------------------------------------------------------------------------
! The cube's values are already stored as a 1d array
!-----------------------------------------------------------------------------
values(:) = cube%values(:)

RETURN

END SUBROUTINE cube_get_data_1d


SUBROUTINE cube_get_data_2d(cube, values)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Get the data from a cube where the data is a 2d array
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(data_cube), INTENT(IN) :: cube  ! The cube containing the data

REAL, INTENT(OUT) :: values(:,:)  ! The array to put values in

INTEGER :: shape2D(2) ! Constant size array required by RESHAPE intrinsic

!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that values is the right shape for the data in the cube
!-----------------------------------------------------------------------------
! First make sure that cube%shape has the correct number of dimensions
IF ( SIZE(cube%SHAPE) /= 2 )                                                  &
  CALL log_fatal("cube_get_data_2d",                                          &
                 "values has a different rank to the cube")

! Check that all the dimensions of values have the correct size
IF ( .NOT. ALL(SHAPE(values) == cube%SHAPE(1:2)) )                            &
  CALL log_fatal("cube_get_data_2d",                                          &
                 "values is not the same shape as the cube")

!-----------------------------------------------------------------------------
! Reshape the cube's data to the expected shape
!-----------------------------------------------------------------------------
shape2D(1:2) = cube%SHAPE(1:2)
values = RESHAPE(cube%values, shape2D)

RETURN

END SUBROUTINE cube_get_data_2d


SUBROUTINE cube_get_data_3d(cube, values)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Get the data from a cube where the data is a 3d array
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(data_cube), INTENT(IN) :: cube  ! The cube containing the data

REAL, INTENT(OUT) :: values(:,:,:)  ! The array to put values in

INTEGER :: shape3D(3) ! Constant size array required by RESHAPE intrinsic

!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that values is the right shape for the data in the cube
!-----------------------------------------------------------------------------
! First make sure that cube%shape has the correct number of dimensions
IF ( SIZE(cube%SHAPE) /= 3 )                                                  &
  CALL log_fatal("cube_get_data_3d",                                          &
                 "values has a different rank to the cube")

! Check that all the dimensions of values have the correct size
IF ( .NOT. ALL(SHAPE(values) == cube%SHAPE(1:3)) )                            &
  CALL log_fatal("cube_get_data_3d",                                          &
                 "values is not the same shape as the cube")

!-----------------------------------------------------------------------------
! Reshape the cube's data to the expected shape
!-----------------------------------------------------------------------------
shape3D(1:3) = cube%SHAPE(1:3)
values = RESHAPE(cube%values, shape3D)

RETURN

END SUBROUTINE cube_get_data_3d


SUBROUTINE cube_get_data_4d(cube, values)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Get the data from a cube where the data is a 4d array
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(data_cube), INTENT(IN) :: cube  ! The cube containing the data

REAL, INTENT(OUT) :: values(:,:,:,:)  ! The array to put values in

INTEGER :: shape4D(4) ! Constant size array required by RESHAPE intrinsic

!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that values is the right shape for the data in the cube
!-----------------------------------------------------------------------------
! First make sure that cube%shape has the correct number of dimensions
IF ( SIZE(cube%SHAPE) /= 4 )                                                  &
  CALL log_fatal("cube_get_data_4d",                                          &
                 "values has a different rank to the cube")

! Check that all the dimensions of values have the correct size
IF ( .NOT. ALL(SHAPE(values) == cube%SHAPE(1:4)) )                            &
  CALL log_fatal("cube_get_data_4d",                                          &
                 "values is not the same shape as the cube")

!-----------------------------------------------------------------------------
! Reshape the cube's data to the expected shape
!-----------------------------------------------------------------------------
shape4D(1:4) = cube%SHAPE(1:4)
values = RESHAPE(cube%values, shape4D)

RETURN

END SUBROUTINE cube_get_data_4d


SUBROUTINE cube_get_data_5d(cube, values)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Get the data from a cube where the data is a 5d array
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(data_cube), INTENT(IN) :: cube  ! The cube containing the data

REAL, INTENT(OUT) :: values(:,:,:,:,:)  ! The array to put values in

INTEGER :: shape5D(5) ! Constant size array required by RESHAPE intrinsic

!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that values is the right shape for the data in the cube
!-----------------------------------------------------------------------------
! First make sure that cube%shape has the correct number of dimensions
IF ( SIZE(cube%SHAPE) /= 5 )                                                  &
  CALL log_fatal("cube_get_data_5d",                                          &
                 "values has a different rank to the cube")

! Check that all the dimensions of values have the correct size
IF ( .NOT. ALL(SHAPE(values) == cube%SHAPE(1:5)) )                            &
  CALL log_fatal("cube_get_data_5d",                                          &
                 "values is not the same shape as the cube")

!-----------------------------------------------------------------------------
! Reshape the cube's data to the expected shape
!-----------------------------------------------------------------------------
shape5D(1:5) = cube%SHAPE(1:5)
values = RESHAPE(cube%values, shape5D)

RETURN

END SUBROUTINE cube_get_data_5d


SUBROUTINE cube_get_data_6d(cube, values)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Get the data from a cube where the data is a 6d array
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(data_cube), INTENT(IN) :: cube  ! The cube containing the data

REAL, INTENT(OUT) :: values(:,:,:,:,:,:)  ! The array to put values in

INTEGER :: shape6D(6) ! Constant size array required by RESHAPE intrinsic

!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that values is the right shape for the data in the cube
!-----------------------------------------------------------------------------
! First make sure that cube%shape has the correct number of dimensions
IF ( SIZE(cube%SHAPE) /= 6 )                                                  &
  CALL log_fatal("cube_get_data_6d",                                          &
                 "values has a different rank to the cube")

! Check that all the dimensions of values have the correct size
IF ( .NOT. ALL(SHAPE(values) == cube%SHAPE(1:6)) )                            &
  CALL log_fatal("cube_get_data_6d",                                          &
                 "values is not the same shape as the cube")

!-----------------------------------------------------------------------------
! Reshape the cube's data to the expected shape
!-----------------------------------------------------------------------------
shape6D(1:6) = cube%SHAPE(1:6)
values = RESHAPE(cube%values, shape6D)

RETURN

END SUBROUTINE cube_get_data_6d


SUBROUTINE cube_get_data_7d(cube, values)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Get the data from a cube where the data is a 7d array
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(data_cube), INTENT(IN) :: cube  ! The cube containing the data

REAL, INTENT(OUT) :: values(:,:,:,:,:,:,:)  ! The array to put values in

INTEGER :: shape7D(7) ! Constant size array required by RESHAPE intrinsic

!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that values is the right shape for the data in the cube
!-----------------------------------------------------------------------------
! First make sure that cube%shape has the correct number of dimensions
IF ( SIZE(cube%SHAPE) /= 7 )                                                  &
  CALL log_fatal("cube_get_data_7d",                                          &
                 "values has a different rank to the cube")

! Check that all the dimensions of values have the correct size
IF ( .NOT. ALL(SHAPE(values) == cube%SHAPE(1:7)) )                            &
  CALL log_fatal("cube_get_data_7d",                                          &
                 "values is not the same shape as the cube")

!-----------------------------------------------------------------------------
! Reshape the cube's data to the expected shape
!-----------------------------------------------------------------------------
shape7D(1:7) = cube%SHAPE(1:7)
values = RESHAPE(cube%values, shape7D)

RETURN

END SUBROUTINE cube_get_data_7d
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

SUBROUTINE cube_free(cube)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Frees the memory used by the given cube. Data will no longer be able to
!   be read
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
TYPE(data_cube), INTENT(INOUT) :: cube  ! The cube to free


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Just deallocate the arrays
!-----------------------------------------------------------------------------
IF ( ASSOCIATED(cube%SHAPE) ) THEN
  DEALLOCATE(cube%SHAPE)
  NULLIFY(cube%SHAPE)
END IF

IF ( ASSOCIATED(cube%values) ) THEN
  DEALLOCATE(cube%values)
  NULLIFY(cube%values)
END IF

RETURN

END SUBROUTINE cube_free


END MODULE data_cube_mod
