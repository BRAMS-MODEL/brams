!******************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT****************************************

MODULE file_ascii_generic_sync_mod

!-----------------------------------------------------------------------------
! Description:
!   Provides an abstracted subroutine (file_sync()) which forces a
!   synchronisation of file contents with the disk (if possible). This requires
!   compiler specific extensions, so only actually attempts a sync on compilers
!   for which a method has been specifically written. Currently, these are: -
!
!     * GNU GFORTRAN Compiler
!     * Intel IFORT Compiler
!     * Cray Fortran Compiler
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90. (+ compiler specific extensions)
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Define list of compilers we have specifically supported synchronisation for: -


USE, INTRINSIC :: ISO_C_BINDING

IMPLICIT NONE

PRIVATE

PUBLIC :: file_sync


! Use iso_c_binding to CALL the standard libc function fsync() to force
! synchronisation with the disk.

INTERFACE
FUNCTION bind_c_fsync(fd) BIND(c,NAME="fsync")
IMPORT :: C_INT
INTEGER(KIND=C_INT), VALUE :: fd
INTEGER(KIND=C_INT) :: bind_c_fsync
END FUNCTION bind_c_fsync
END INTERFACE


!------------------------------------------------------------------------------!
CONTAINS
!------------------------------------------------------------------------------!

SUBROUTINE file_sync(unitno)

!-----------------------------------------------------------------------------
! Description:
!   Sync the file with disk (if we can).
!
! Method:
!   Call bind_c_fsync(get_fd_from_unit(unitno)) if that is defined.
!   Warn on failure, rather than fail - this allows the subroutine to do nothing
!   when built with a compiler which does not have an appropriate extension to
!   implement the unit number to file descriptor conversion in get_fd_from_unit()
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

USE logging_mod, ONLY: log_info, log_warn

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: unitno

! Local variables
INTEGER :: ret = -1
CHARACTER(LEN=40) :: message

!-----------------------------------------------------------------------------

! Only sync on compilers we have specificially implemented:

ret = bind_c_fsync(get_fd_from_unit(unitno))

IF (ret == 0) THEN
  WRITE(message,'(A,I0,A)') 'Succeeded to sync unit ',unitno,' with disk'
  CALL log_info("file_sync", message)
ELSE
  WRITE(message,'(A,I0,A)') 'Failed to sync unit ',unitno,' with disk'
  CALL log_warn("file_sync", message)
END IF

END SUBROUTINE file_sync

!------------------------------------------------------------------------------!

FUNCTION get_fd_from_unit(unitno)

!-----------------------------------------------------------------------------
! Description:
!   Convert Fortran unit number to C file descriptor
!
! Method:
!   If there is not a compiler-specific implementation of this function for
!   the current compiler, it will return -1 (which is the error code for
!   implemented compilers)
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------


IMPLICIT NONE

! Function arguments
INTEGER, INTENT(IN) :: unitno

! Return type
INTEGER :: get_fd_from_unit

! Local variables

INTEGER :: tmp_unitno=-1


!-----------------------------------------------------------------------------


! use GNU extension FNUM()
tmp_unitno = fnum(unitno)


get_fd_from_unit = tmp_unitno

END FUNCTION get_fd_from_unit

!------------------------------------------------------------------------------!

END MODULE file_ascii_generic_sync_mod
